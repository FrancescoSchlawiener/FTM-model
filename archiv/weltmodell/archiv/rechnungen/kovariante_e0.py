# -*- coding: utf-8 -*-
"""
FTM - kovariante Fassung der Einzeltrichter-Hierarchie E0 (Wirkung FTM-AS-E0).
Rechnungen zu forschung/37-kovariante-e0.md
Aufruf:  python kovariante_e0.py > kovariante_e0.txt      (numpy, scipy)

Abschnitte
  0  Grundgroessen, Sonnenumgebung (a0 = 1,2e-10 gemessen)
  1  Gezeitentensor und Glaettungsradius: Invariante, Normierung, effektive und baryonische Gezeit
  2  Zonenmittel der Gezeit: Identitaet, Sonne selbst, Planeten, Nachbarsterne, Gravitationswellen
  3  Komar-Masse, Fixpunkt fuer l
  4  Newton-Grenzfall: diskrete Wirkung auf 3D-Gitter, Variationsprobe, Energie = E0, Kraftgesetz
  5  Sonne im Milchstrassenfeld: gamma(r), Zusatz an den Planeten, Q2 (Cassini), Periheldrehungen,
     weite Doppelsterne, Oortsche Wolke
  6  Groesse der kovarianten Korrekturen
  7  Gezeit aus der vollen Riemann-Projektion: Kopplung an Gravitationswellen (Abschaetzung)
  8  Schwellenuebergang: Energiesprung beim Einsetzen der Glaettung, l-Kraft
Eintraege mit [G] sind Literatur-/Beobachtungswerte aus dem Gedaechtnis (nicht verifiziert).
"""
import io, math, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import k1

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(all="ignore")

# ------------------------------------------------------------------ Konstanten
G    = 6.67430e-11
c    = 299792458.0
AU   = 1.495978707e11
PC   = 3.0856775814913673e16
KPC  = 1e3 * PC
MPC  = 1e6 * PC
YR   = 3.15576e7
MAS  = 180 / math.pi * 3600e3
MSUN = 1.98847e30
GMs  = 1.32712440018e20
a0   = 1.2e-10                      # gemessener Wert (durchgehend)
a0H  = c * (70e3 / MPC) / (2 * math.pi)   # Modellwert cH0/2π bei H0 = 70 (nur Vergleich)
KK   = 1e3 / KPC                    # km/s/kpc -> 1/s

def H(t): print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

def nu(x):
    x = np.maximum(np.asarray(x, dtype=float), 1e-300)
    return 1.0 + 1.0 / np.expm1(np.sqrt(x))

def xdnudx(x):
    """dnu/dln g = -nu(nu-1) s/2, s = sqrt(g/a0)"""
    s = math.sqrt(x); n = float(nu(x)); return -n * (n - 1) * s / 2

_XL, _WL = np.polynomial.legendre.leggauss(48)
def Qs_of_y(y):
    """Q_s(y) = ∫_0^y (nu(√y')-1) dy' = 4∫_0^s t³/(e^t-1) dt,  s = y^(1/4); glatt (Gauss-Legendre)"""
    y = np.asarray(y, dtype=float); s = np.maximum(y, 0.0) ** 0.25
    t = 0.5 * s[..., None] * (1 + _XL)
    f = np.where(t > 0, 4 * t ** 3 / np.expm1(np.maximum(t, 1e-300)), 0.0)
    return 0.5 * s * np.sum(_WL * f, axis=-1)

def inv_rar(gob, a0_):
    lo, hi = 1e-8, 1e4
    for _ in range(200):
        m = math.sqrt(lo * hi)
        if m * float(nu(m)) < gob / a0_: lo = m
        else: hi = m
    return math.sqrt(lo * hi) * a0_

def plummer_frac(x):
    """Anteil der Plummer-Masse innerhalb r = x·l"""
    return x ** 3 / (1 + x ** 2) ** 1.5

def tf(T): return T - np.trace(T) / 3 * np.eye(3)
def nrm(T): return math.sqrt(float(np.sum(T ** 2)))
SQ32 = math.sqrt(1.5)

# ============================================================================ 0
H("0   Grundgroessen und Sonnenumgebung (a0 = 1,2e-10 m/s², gemessen)")
v0, R0 = 233e3, 8.178 * KPC                               # [G]
gobs = v0 ** 2 / R0
gNe = inv_rar(gobs, a0); be = gNe / a0; nue = float(nu(be)); Ke = xdnudx(be) / nue
print(f"   v0 = 233 km/s, R0 = 8,178 kpc [G]  ->  g_obs = {gobs:.4e} m/s² = {gobs/a0:.3f} a0")
print(f"   Newtonsches Elternfeld aus RAR-Umkehr: g_Ne = {gNe:.4e} = {be:.4f} a0,  nu_e = {nue:.4f},  K = dln nu/dln g = {Ke:.4f}")
rs = math.sqrt(GMs / a0); re = math.sqrt(GMs / gNe)
print(f"   Sonne: r* = {rs/AU:.1f} AE,  r_e = √(GM/g_Ne) = {re/AU:.1f} AE")
ce = nue * (1 + Ke / 3) - 1
print(f"   lineare Antwort des kollektiven Feldes auf ein schwaches eigenes Feld (Kugelmittel):"
      f"  c_e = nu_e(1+K/3) - 1 = {ce:.4f}")

# ============================================================================ 1
H("1   Gezeitentensor und Glaettungsradius  l³ = G M / (√(3/2)·|Ē^TF|)")
A_, B_ = 15.3 * KK, -11.9 * KK                            # [G] Oort-Konstanten
rho0 = 0.10 * MSUN / PC ** 3                              # [G] dynamische Dichte in der Scheibenmitte
rhob = 0.084 * MSUN / PC ** 3                             # [G] baryonische Dichte in der Scheibenmitte
Om2 = (A_ - B_) ** 2
T_eff = np.diag([(A_ - B_) * (-3 * A_ - B_), Om2, 4 * math.pi * G * rho0 + 2 * (A_ ** 2 - B_ ** 2)])
dlnv = -(A_ + B_) / (A_ - B_)
dlngobs = 2 * dlnv - 1
dlngN = dlngobs / (1 + Ke)
PhiRR_N = gNe / R0 * dlngN; Phipp_N = gNe / R0
T_bar = np.diag([PhiRR_N, Phipp_N, 4 * math.pi * G * rhob - PhiRR_N - Phipp_N])
tideJ = 4 * A_ * (A_ - B_)
print(f"   Normierung: Punktmasse als Eltern, Kreisbahn: Ē^TF = (GM/R³)·diag(-2,1,1) -> √(3/2)|Ē^TF| = √(3/2)·√6 = {SQ32*math.sqrt(6):.4f}·GM/R³")
print(f"               = Jacobi-Wert 4A(A-B) = 3Ω² = 3GM/R³ der Kepler-Bahn (A = 3Ω/4, B = -Ω/4).")
print(f"   effektive Gezeit (Oort [G], ρ0 = 0,10 M☉/pc³): diag = ({T_eff[0,0]:.3e}, {T_eff[1,1]:.3e}, {T_eff[2,2]:.3e}) s⁻²")
print(f"      spurfrei: diag = ({tf(T_eff)[0,0]:.3e}, {tf(T_eff)[1,1]:.3e}, {tf(T_eff)[2,2]:.3e}),  |TF| = {nrm(tf(T_eff)):.3e},  √(3/2)|TF| = {SQ32*nrm(tf(T_eff)):.3e}")
print(f"   baryonische Newton-Gezeit (aus g_Ne, dln v/dln R = {dlnv:.4f}, dln g_N/dln R = {dlngN:.4f}, ρ_b = 0,084 M☉/pc³ [G]):")
print(f"      diag = ({T_bar[0,0]:.3e}, {T_bar[1,1]:.3e}, {T_bar[2,2]:.3e}) s⁻²;  spurfrei |TF| = {nrm(tf(T_bar)):.3e},  √(3/2)|TF| = {SQ32*nrm(tf(T_bar)):.3e}")
print(f"   Jacobi-Wert 4A(A-B) (Bahnrahmen, enthaelt Ω²) = {tideJ:.3e} s⁻²")
def ell_of(GM, Tn): return (GM / Tn) ** (1 / 3)
EB = SQ32 * nrm(tf(T_bar)); EE = SQ32 * nrm(tf(T_eff))
ELLS = [("kovariant: √(3/2)|∂∂Φ̄_N^TF| (baryonisch, gewaehlt)", EB),
        ("√(3/2)|Φ_ij^TF| effektiv", EE),
        ("volle Frobenius-Norm effektiv (mit Spur)", nrm(T_eff)),
        ("Jacobi 4A(A-B) (Bahnrahmen)", tideJ)]
print(f"\n   {'Gezeitenmass':52s} {'|E| [s⁻²]':>11s} {'l [pc]':>8s} {'l [AE]':>10s} {'l/r*':>7s} {'l/r_e':>7s}")
for nm, Tn in ELLS:
    l_ = ell_of(GMs, Tn)
    print(f"   {nm:52s} {Tn:11.3e} {l_/PC:8.3f} {l_/AU:10.4g} {l_/rs:7.2f} {l_/re:7.2f}")
ellB = ell_of(GMs, EB); ellJ = ell_of(GMs, tideJ); ellE = ell_of(GMs, EE)

# ============================================================================ 2
H("2   Zonenmittel der Gezeit  ⟨Ē^TF⟩_l = ∫ K_l(x) ∂∂Φ̄_N^TF(x) d³x   (Plummer-Gewicht um den Schwerpunkt)")
# Identitaet: Quelle im Abstand d -> Schalen r>d tragen 0, Schalen r<d den Zentralwert -> ⟨T⟩ = T(0)·F(d/l)
_xg64, _wg64 = np.polynomial.legendre.leggauss(64)
_bp = [-1.0, 0.0] + [1 - 10.0 ** (-k) for k in range(1, 15)] + [1.0]
_MU4 = np.concatenate([0.5 * (b - a) * _xg64 + 0.5 * (b + a) for a, b in zip(_bp[:-1], _bp[1:])])
_W4 = np.concatenate([0.5 * (b - a) * _wg64 for a, b in zip(_bp[:-1], _bp[1:])])
def shell_avg_Tzz(r, d):
    """Kugelmittel von ∂z²(1/|x - d ẑ|) ueber die Schale |x| = r (Gauss-Legendre in μ, bei μ -> 1 verdichtet)"""
    mu, w = _MU4, _W4
    dz = r * mu - d; R2 = r ** 2 * (1 - mu ** 2) + dz ** 2
    Tzz = 3 * dz ** 2 / R2 ** 2.5 - 1 / R2 ** 1.5
    return 0.5 * np.sum(w * Tzz)
_XL4, _WL4 = np.polynomial.legendre.leggauss(200)
def plummer_avg_Tzz(d, ell):
    tot = 0.0; big = max(d, ell)
    segs = [(0, 0.5 * d), (0.5 * d, 0.9 * d), (0.9 * d, d), (d, 1.1 * d), (1.1 * d, 2 * d), (2 * d, 50 * big), (50 * big, 5000 * big)]
    for a, b in segs:
        r = 0.5 * (b - a) * _XL4 + 0.5 * (b + a); wr = 0.5 * (b - a) * _WL4
        dens = 4 * math.pi * r ** 2 * 3 / (4 * math.pi * ell ** 3) * (1 + r ** 2 / ell ** 2) ** -2.5
        tot += np.sum(wr * dens * np.array([shell_avg_Tzz(ri, d) for ri in r]))
    return tot
print("   Probe der Identitaet ⟨∂z²(1/|x-d|)⟩_l = (2/d³)·F(d/l),  F(x) = x³/(1+x²)^(3/2):")
for d, ell in [(0.3, 1.0), (1.0, 1.0), (3.0, 1.0)]:
    num = plummer_avg_Tzz(d, ell); ana = 2 / d ** 3 * plummer_frac(d / ell)
    print(f"      d/l = {d/ell:4.1f}:  Quadratur {num:.6f},  Formel {ana:.6f}")
print("   Eigenes Feld (d = 0, kugelsymmetrisch): Zonenmittel des spurfreien Anteils = 0 exakt.")

print(f"\n   Beitraege zur Gezeit der Sonne (l = {ellB/PC:.3f} pc; Mass √(3/2)|TF| = 3GM/d³ fuer eine Punktquelle):")
print(f"   {'Quelle':22s} {'d':>10s} {'Punktwert [s⁻²]':>16s} {'Zonenmittel [s⁻²]':>18s} {'rel. zu Galaxie':>16s}")
PL = [("Merkur", 0.387, 2.2032e13), ("Venus", 0.723, 3.2486e14), ("Erde+Mond", 1.0, 4.0350e14), ("Mars", 1.524, 4.2828e13),
      ("Jupiter", 5.2026, 1.26686534e17), ("Saturn", 9.5549, 3.7931187e16), ("Uranus", 19.218, 5.793939e15), ("Neptun", 30.11, 6.836529e15)]
pt_tot = 0.0; av_tot = 0.0
for nm, dA, GMp in PL:
    d = dA * AU; pv = 3 * GMp / d ** 3; av = pv * plummer_frac(d / ellB)
    pt_tot += pv; av_tot += av
    print(f"   {nm:22s} {dA:8.3f} AE {pv:16.3e} {av:18.3e} {av/EB:16.2e}")
print(f"   {'Summe Planeten':22s} {'':>10s} {pt_tot:16.3e} {av_tot:18.3e} {av_tot/EB:16.2e}")
print(f"   {'Galaxie (baryonisch)':22s} {'':>10s} {EB:16.3e} {EB:18.3e} {1.0:16.2e}")
lpt = ell_of(GMs, pt_tot + EB)
print(f"   Mit dem Punktwert statt des Zonenmittels waere l = {lpt/AU:.1f} AE = {lpt/rs:.4f} r*  (Cassini verlangt l ≥ 1,30 r*).")
print(f"   Planetenbeitrag im Zonenmittel ≈ 3·ΣGM_p/l³ (unabhaengig von d): {3*sum(p[2] for p in PL)/ellB**3:.3e} s⁻²"
      f"  -> δl/l ≤ {av_tot/EB/3:.1e}")

# Nachbarsterne (Koernigkeit der geglaetteten Sterndichte)
print("\n   Nachbarsterne: Poisson-Feld n = 0,1 Sterne/pc³ [G], Massen log-gleichverteilt 0,1-1 M☉, je ueber l_j geglaettet;")
print("   Gezeit eines Plummer-Sterns im Abstand d, gemittelt ueber die Sonnenzone ≈ Hesse-Matrix des Potentials mit ε² = l_j² + l_s².")
rng = np.random.default_rng(1)
def star_tide(pos, GMj, eps):
    d2 = np.sum(pos ** 2, axis=1) + eps ** 2
    w1 = GMj / d2 ** 1.5; w2 = 3 * GMj / d2 ** 2.5
    return np.sum(w1) * np.eye(3) - np.einsum("n,ni,nj->ij", w2, pos, pos)
Rbox = 15 * PC; ratios = []
for it in range(600):
    N = rng.poisson(0.1 * 4 / 3 * math.pi * (Rbox / PC) ** 3)
    u = rng.random(N) ** (1 / 3) * Rbox; v = rng.normal(size=(N, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    pos = u[:, None] * v
    m = 10 ** rng.uniform(-1, 0, N); GMj = m * GMs
    lj = ell_of(GMj, EB); ls = ellB
    for _ in range(3):
        Ttot = T_bar + star_tide(pos, GMj, np.sqrt(lj ** 2 + ls ** 2))
        ls = ell_of(GMs, SQ32 * nrm(tf(Ttot)))
    ratios.append(ls / ellB)
ratios = np.array(ratios)
print(f"   l_Sonne / l_glatt: Median {np.median(ratios):.3f},  5 % {np.percentile(ratios,5):.3f},  95 % {np.percentile(ratios,95):.3f},"
      f"  Minimum {ratios.min():.3f}  ->  l_min = {ratios.min()*ellB/rs:.1f} r*")

# Gravitationswellen: Zonenmittel einer ebenen Welle = Fourier-Transformierte der Plummer-Dichte
def Pi(kl): return kl * k1(kl)
def Pi_num(kl):
    r = np.linspace(1e-6, 400, 400001)
    dens = 3 * r ** 2 * (1 + r ** 2) ** -2.5
    return np.trapezoid(dens * np.sinc(kl * r / math.pi), r)
print("\n   Zonenmittel einer ebenen Welle e^{ik·x}: Π(kl) = kl·K1(kl);  Probe (Quadratur):  " +
      ",  ".join(f"kl={x}: {Pi(x):.6f} / {Pi_num(x):.6f}" for x in [0.3, 1.0, 3.0, 10.0]))
print(f"   Gezeit einer Welle ½ω²h_c gegen die galaktische Gezeit (h_c = 2,4e-15 (f·yr)^(-2/3) [G, PTA-Spektrum, unter 2 nHz extrapoliert]):")
print(f"   {'f [Hz]':>10s} {'λ [pc]':>9s} {'kl':>9s} {'½ω²h_c':>11s} {'Π(kl)':>10s} {'Zonenmittel/Gal.':>17s}")
for f in [5e-10, 1e-9, 2e-9, 3e-9, 1e-8, 3.17e-8, 1e-7]:
    hc = 2.4e-15 * (f * YR) ** (-2 / 3); om = 2 * math.pi * f; kl = om / c * ellB
    Eg = 0.5 * om ** 2 * hc
    print(f"   {f:10.2e} {c/f/PC:9.3f} {kl:9.3f} {Eg:11.3e} {Pi(kl):10.3e} {Eg*Pi(kl)/EB:17.3e}")
kl_ligo = 2 * math.pi * 100 / c * ellB
print(f"   LIGO-Band (100 Hz, h = 1e-21): ½ω²h = {0.5*(2*math.pi*100)**2*1e-21:.2e} s⁻² (Punktwert {0.5*(2*math.pi*100)**2*1e-21/EB:.1e}-fach der Galaxie),"
      f" kl = {kl_ligo:.2e}, Π ~ exp(-{kl_ligo:.1e}) = 0")

# ============================================================================ 3
H("3   Masse als Komar-Fluss und Fixpunkt fuer l")
print("   Statisch, σ = ln N (q -> 0 im Innern kompakter Koerper), Schwarzschild: M_σ = (c²/4πG)∮ n·Dσ dA = M/N(r),")
print("   M_K = (c²/4πG)∮ e^σ n·Dσ dA = M (unabhaengig von der Flaeche).")
for nm, GM, R in [("Sonne, Oberflaeche", GMs, 6.957e8), ("Weisser Zwerg 0,6 M☉, 8700 km", 0.6 * GMs, 8.7e6), ("Neutronenstern 1,4 M☉, 12 km", 1.4 * GMs, 1.2e4)]:
    N_ = math.sqrt(1 - 2 * GM / (R * c ** 2))
    print(f"   {nm:32s}: N = {N_:.8f},  M_σ/M = {1/N_:.8f},  M_K/M = 1,  Unterschied in l: {(1/N_)**(1/3)-1:.2e}")
ell = 0.1 * ellB
for k in range(40):
    Et = EB + sum(3 * GMp / (dA * AU) ** 3 * plummer_frac(dA * AU / ell) for _, dA, GMp in PL)
    ell_new = ell_of(GMs, Et)
    if abs(ell_new / ell - 1) < 1e-14: break
    ell = ell_new
print(f"   Fixpunkt-Iteration l -> (GM/√(3/2)|⟨E⟩_l|)^(1/3) mit Planeten (Start 0,1 l): {k} Schritte, l = {ell/PC:.6f} pc,"
      f"  l/l_ohne Planeten - 1 = {ell/ellB-1:.2e}")

# ============================================================================ 4
H("4   Newton-Grenzfall der Wirkung: diskrete Probe auf einem 3D-Gitter (Einheiten G = a0 = 1, Gitterweite 1)")
print("   Diskrete Fassung von  L_NR = (1/8πG)[(∇Ψ)² - 2∇Φ·∇Ψ + (∇Φ-∇Φ_N)² + a0²Q_s(|∇Φ̄_N|²/a0²) + 2φ(∇²Φ̄_N - 𝒢[∇²Φ_N])] - ρΦ")
print("   𝒢[s] = (1-χ)s + K_l·∫χ s  (Koerpergebiet χ, Plummer-Kern K_l um den Schwerpunkt, Masse aus dem Fluss)")
n = 14; Nn = n ** 3
I1 = sp.identity(n, format="csr")
D1 = sp.diags([np.ones(n), -np.ones(n)], [0, -1], shape=(n + 1, n), format="csr")       # Vorwaertsdifferenz, Dirichlet 0
C1 = sp.diags([0.5 * np.ones(n - 1), -0.5 * np.ones(n - 1)], [1, -1], shape=(n, n), format="csr")
Dop = sp.vstack([sp.kron(sp.kron(D1, I1), I1), sp.kron(sp.kron(I1, D1), I1), sp.kron(sp.kron(I1, I1), D1)]).tocsr()
Mop = sp.vstack([sp.kron(sp.kron(C1, I1), I1), sp.kron(sp.kron(I1, C1), I1), sp.kron(sp.kron(I1, I1), C1)]).tocsr()
Lap = (Dop.T @ Dop).tocsc(); LU = spla.splu(Lap)
ix = np.indices((n, n, n)).reshape(3, -1).T.astype(float)
cb = np.array([5.0, 6.0, 6.0])
chi = (np.linalg.norm(ix - cb, axis=1) <= 1.01).astype(float)
ell_g = 2.5
Kp = (1 + np.sum((ix - cb) ** 2, axis=1) / ell_g ** 2) ** -2.5; Kp /= Kp.sum()
def rho_of(eps=0.0, eta=None):
    rb = 20.0 * chi / chi.sum()
    rg = np.exp(-np.sum((ix - np.array([9.0, 8.0, 7.0])) ** 2, axis=1) / (2 * 1.5 ** 2)); rg *= 8.0 / rg.sum(); rg *= (1 - chi)
    r_ = rb + rg
    return r_ if eta is None else r_ + eps * eta
def Gop(s): return (1 - chi) * s + Kp * np.dot(chi, s)
def GTop(f): return (1 - chi) * f + chi * np.dot(Kp, f)
def solve_fields(rho):
    PhiN = -4 * math.pi * LU.solve(rho)
    src = -(Lap @ PhiN)
    Phib = -LU.solve(Gop(src))
    g = (Mop @ Phib).reshape(3, -1); y = np.sum(g ** 2, axis=0)
    Qp = nu(np.sqrt(y)) - 1
    phi = LU.solve(Mop.T @ (np.tile(Qp, 3) * (Mop @ Phib)))
    Phi = PhiN + GTop(phi)
    return dict(Psi=Phi.copy(), Phi=Phi, PhiN=PhiN, Phib=Phib, phi=phi)
def action(F, rho):
    DPs, DPh, DPN = Dop @ F["Psi"], Dop @ F["Phi"], Dop @ F["PhiN"]
    g = (Mop @ F["Phib"]).reshape(3, -1); y = np.sum(g ** 2, axis=0)
    con = -(Lap @ F["Phib"]) - Gop(-(Lap @ F["PhiN"]))
    return (np.dot(DPs, DPs) - 2 * np.dot(DPh, DPs) + np.dot(DPh - DPN, DPh - DPN) + np.sum(Qs_of_y(y))
            + 2 * np.dot(F["phi"], con)) / (8 * math.pi) - np.dot(rho, F["Phi"])
def V_E0(rho):
    F = solve_fields(rho); DPN = Dop @ F["PhiN"]
    g = (Mop @ F["Phib"]).reshape(3, -1); y = np.sum(g ** 2, axis=0)
    return -(np.dot(DPN, DPN) + np.sum(Qs_of_y(y))) / (8 * math.pi), F
rho = rho_of(); F0 = solve_fields(rho); S0 = action(F0, rho)
gb = (Mop @ F0["Phib"]).reshape(3, -1); yb = np.sqrt(np.sum(gb ** 2, axis=0))
print(f"   Gitter {n}³, Koerper 20 (7 Knoten), Gas 8, l = {ell_g};  |∇Φ̄_N|/a0 im Gitter: {yb.min():.3f} ... {yb.max():.3f}")
rng2 = np.random.default_rng(7)
print(f"   Richtungsableitungen der Wirkung in der Loesung (zufaellige Richtungen, zentrale Differenzen):")
for key in ["Psi", "Phi", "PhiN", "Phib", "phi"]:
    e = rng2.normal(size=Nn); e /= np.linalg.norm(e); h_ = 1e-4 * max(1.0, np.abs(F0[key]).max())
    Fp = {k: v.copy() for k, v in F0.items()}; Fm = {k: v.copy() for k, v in F0.items()}
    Fp[key] += h_ * e; Fm[key] -= h_ * e
    dS = (action(Fp, rho) - action(Fm, rho)) / (2 * h_)
    ref = abs(np.dot(Dop @ F0["PhiN"], Dop @ e)) / (4 * math.pi)
    print(f"      δS/δ{key:5s}·e = {dS: .3e}   (Vergleichsgroesse (1/4π)|∇Φ_N·∇e| = {ref:.3e})")
Vs, Fv = V_E0(rho)
print(f"   Wirkung in der Loesung S = {S0:.10f},  -V_E0 = {-Vs:.10f},  Differenz {S0+Vs:.2e}")
GT = GTop(F0["phi"])
print(f"   Strukturpotential im Koerpergebiet (𝒢†φ): min {GT[chi>0].min():.10f}, max {GT[chi>0].max():.10f},  K·φ = {np.dot(Kp, F0['phi']):.10f}")
eta = np.exp(-np.sum((ix - np.array([6.0, 6.5, 6.0])) ** 2, axis=1) / (2 * 1.2 ** 2)); eta /= eta.sum()
e_ = 1e-4
dV = (V_E0(rho_of(e_, eta))[0] - V_E0(rho_of(-e_, eta))[0]) / (2 * e_)
print(f"   Kraftgesetz: dV_E0/dε (neu geloest) = {dV:.10f}")
print(f"                Σ Φ·η mit Φ = Φ_N + 𝒢†φ  = {np.dot(F0['Phi'], eta):.10f}")
print(f"                Σ (Φ_N + φ)·η (ohne Zonenmittel im Koerper) = {np.dot(F0['PhiN'] + F0['phi'], eta):.10f}")

# ============================================================================ 5
H("5   Sonne im Milchstrassenfeld mit dem kovarianten Glaettungsradius")
MU, WMU = np.polynomial.legendre.leggauss(2000)
def g_avg(r, GM, gNe_, a0_, sm):
    """kugelgemittelte Einwaertsbeschleunigung (Flussidentitaet): Newton + Fluss der kollektiven Struktur,
       Koerper im kollektiven Feld ueber sm geglaettet (E0)"""
    r = np.atleast_1d(np.asarray(r, dtype=float)); gN = GM / r ** 2
    gs = (GM * r / (r ** 2 + sm ** 2) ** 1.5)[:, None]; m = MU[None, :]
    gm = np.sqrt(np.maximum(gs ** 2 + gNe_ ** 2 + 2 * gs * gNe_ * m, 0.0))
    return gN + 0.5 * np.sum(WMU[None, :] * (nu(gm / a0_) - 1) * (gs + gNe_ * m), axis=1)
def J_quad(beta, lam, Nn=240):
    xg, wg = np.polynomial.legendre.leggauss(Nn)
    bps = {-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 3.5, 6.0, 9.0, 12.0}
    if lam > 0:
        L = math.log(lam); bps |= {L - 1.5, L - 0.7, L - 0.3, L, L + 0.3, L + 0.7, L + 1.5}
    bps = sorted(b for b in bps if -10.0 <= b <= 12.0)
    tseg = [(0.0, 0.5 * math.pi), (0.5 * math.pi, 0.9 * math.pi), (0.9 * math.pi, 0.99 * math.pi), (0.99 * math.pi, math.pi)]
    tot = 0.0
    for u0, u1 in zip(bps[:-1], bps[1:]):
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0); wu = 0.5 * (u1 - u0) * wg; rho_ = np.exp(u)
        for t0, t1 in tseg:
            th = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0); wt = 0.5 * (t1 - t0) * wg
            R_, TH = np.meshgrid(rho_, th, indexing="ij"); cc = np.cos(TH)
            h = R_ / (R_ ** 2 + lam ** 2) ** 1.5
            gm = np.sqrt(np.maximum(h ** 2 + 2 * cc * h + 1, 0.0))
            kern = -3 * (3 * cc ** 2 - 1) * h / R_ ** 2 + (9 * cc - 15 * cc ** 3) / R_ ** 2
            F = 0.5 * (nu(beta * gm) - 1.0) * kern * R_ * np.sin(TH)
            tot += np.einsum("i,j,ij->", wu, wt, F)
    return tot + (4 / 15) * xdnudx(beta) / lam ** 3
Q2tide = -1.5 * tf(T_eff)[0, 0]
print(f"   Quadrupol Q2 = -(3/2)(g_e/r_e)·J(β, λ = l/r_e);  Cassini: Q2 = (3 ± 3)·10⁻²⁷ s⁻² (Hees et al. 2014); galaktische Gezeit: {Q2tide:.2e} s⁻²")
print(f"   {'Glaettungsradius':44s} {'l [pc]':>7s} {'λ':>7s} {'J':>11s} {'Q2 [s⁻²]':>11s} {'Q2/σ_Cassini':>13s} {'Q2 (a0=cH0/2π)':>15s}")
gNeH = inv_rar(gobs, a0H); reH = math.sqrt(GMs / gNeH)
Q2c = {}
for nm, l_ in [("kovariant, baryonisch (gewaehlt)", ellB), ("kovariant, effektiv", ellE), ("Jacobi 4A(A-B)", ellJ),
               ("Minimum Nachbarsterne (0,90 l)", 0.90 * ellB)]:
    J = J_quad(be, l_ / re); Q2 = -1.5 * (gNe / re) * J
    lH = l_ * (a0H / a0) ** 0; QH = -1.5 * (gNeH / reH) * J_quad(gNeH / a0H, lH / reH)
    Q2c[nm] = Q2
    print(f"   {nm:44s} {l_/PC:7.3f} {l_/re:7.2f} {J:11.4e} {Q2:11.3e} {Q2/3e-27:13.1e} {QH:15.3e}")
print(f"   Oberflaechenterm allein (4/15)(dν/dln g)/λ³·(-3/2)g_e/r_e bei l_bar: {-1.5*(gNe/re)*(4/15)*xdnudx(be)/(ellB/re)**3:.3e} s⁻²")

print(f"\n   γ(r) = ⟨g⟩/g_N der Sonne (l = {ellB/PC:.3f} pc) und Zusatz ohne den gleichfoermigen Anteil (der Sonne und Planeten gleich beschleunigt):")
print(f"   {'r [AE]':>9s} {'g_N [m/s²]':>11s} {'γ - 1 (Quadratur)':>18s} {'γ - 1 (linear: c_e r³/(r²+l²)^1,5)':>36s} {'δa [m/s²]':>11s}")
for rA in [1, 5.2, 30, 100, 1e3, 5e3, 1e4, 3e4, 5e4, 1e5, 2e5]:
    r_ = rA * AU; gN = GMs / r_ ** 2
    gq = float(g_avg(r_, GMs, gNe, a0, ellB)[0]) / gN - 1
    gl = ce * r_ ** 3 / (r_ ** 2 + ellB ** 2) ** 1.5
    print(f"   {rA:9.4g} {gN:11.3e} {gq:18.3e} {gl:36.3e} {gl*gN:11.3e}")
kap = ce * GMs / ellB ** 3
print(f"   Innen: isotroper Hooke-Zusatz δa = -κ r, κ = c_e GM/l³ = {kap:.3e} s⁻²  (zum Vergleich Λc²/3 bei Λ_obs = 1,1e-52 m⁻²: {1.1e-52*c**2/3:.2e} s⁻²)")

# Periheldrehungen (Gauss-Gleichungen, Bahnmittel)
eps_ = math.radians(23.4393); ra, de = math.radians(266.405), math.radians(-28.936)
xq = np.array([math.cos(de) * math.cos(ra), math.cos(de) * math.sin(ra), math.sin(de)])
egc = np.array([xq[0], xq[1] * math.cos(eps_) + xq[2] * math.sin(eps_), -xq[1] * math.sin(eps_) + xq[2] * math.cos(eps_)])
def dvarpi(aAU, e, inc, Om, om, acc, GM=GMs, N=40000):
    a = aAU * AU; n_ = math.sqrt(GM / a ** 3); p = a * (1 - e * e)
    M = np.linspace(0, 2 * math.pi, N, endpoint=False); E = M.copy()
    for _ in range(50): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    f = 2 * np.arctan2(math.sqrt(1 + e) * np.sin(E / 2), math.sqrt(1 - e) * np.cos(E / 2)); r = a * (1 - e * np.cos(E))
    i, O, w = map(math.radians, (inc, Om, om)); u = w + f
    rh = np.stack([math.cos(O) * np.cos(u) - math.sin(O) * np.sin(u) * math.cos(i),
                   math.sin(O) * np.cos(u) + math.cos(O) * np.sin(u) * math.cos(i), np.sin(u) * math.sin(i)], 1)
    wh = np.array([math.sin(O) * math.sin(i), -math.cos(O) * math.sin(i), math.cos(i)])
    sh = np.cross(wh[None, :], rh)
    av = acc(rh * r[:, None])
    Rc, Sc, Wc = np.sum(av * rh, 1), np.sum(av * sh, 1), av @ wh
    d = math.sqrt(1 - e * e) / (n_ * a * e) * (-Rc * np.cos(f) + Sc * (1 + r / p) * np.sin(f)) \
        + math.tan(i / 2) * r * np.sin(u) * Wc / (n_ * a * a * math.sqrt(1 - e * e))
    return float(np.mean(d))
planets = [("Merkur", 0.387098, 0.205630, 7.0050, 48.331, 29.124, 0.030), ("Venus", 0.723332, 0.006773, 3.3946, 76.680, 54.884, 0.016),
           ("Erde", 1.000000, 0.016709, 0.0, 0.0, 102.937, 0.0019), ("Mars", 1.523679, 0.093405, 1.8497, 49.558, 286.502, 0.00037),
           ("Jupiter", 5.2026, 0.048498, 1.3033, 100.464, 273.867, 0.283), ("Saturn", 9.5549, 0.055546, 2.4845, 113.665, 339.392, 0.0047),
           ("Uranus", 19.2184, 0.046381, 0.7732, 74.006, 96.999, 3.90), ("Neptun", 30.1104, 0.009456, 1.7700, 131.784, 276.336, 5.40)]
Q2B = Q2c["kovariant, baryonisch (gewaehlt)"]
accT = lambda X: -kap * X + (Q2B + Q2tide) * (np.outer(X @ egc, egc) - X / 3)
print(f"\n   Periheldrehungen aus Hooke-Zusatz (κ) + Quadrupol (Glaettung {Q2B:.2e} + Gezeit {Q2tide:.2e} s⁻²); Fehler EPM2011 [G]:")
print(f"   {'Planet':8s} {'Δϖ̇ [mas/Jahr]':>15s} {'σ_EPM [mas/Jahr]':>17s} {'Verhaeltnis':>12s} {'δa_Hooke [m/s²]':>16s} {'δa/g_N':>10s}")
for nm, a, e, i, O, w, err in planets:
    dv = dvarpi(a, e, i, O, w, accT) * YR * MAS
    print(f"   {nm:8s} {dv:15.3e} {err:17.2g} {abs(dv)/err:12.1e} {kap*a*AU:16.3e} {kap*(a*AU)**3/GMs:10.2e}")

# weite Doppelsterne und Oortsche Wolke
Mb = 1.5; GMb = Mb * GMs; ellb = ell_of(GMb, EB)
sv = np.array([2, 5, 7, 10, 15, 20, 30]) * 1e3 * AU
gb_ = g_avg(sv, GMb, gNe, a0, ellb) / (GMb / sv ** 2)
print(f"\n   Weite Doppelsterne (1,5 M☉, l = {ellb/PC:.3f} pc):  " + ",  ".join(f"{s/AU/1e3:.0f} kAE: γ = {g:.5f}" for s, g in zip(sv, gb_)))
ro = np.array([1e4, 2e4, 5e4, 1e5, 2e5]) * AU
go_ = g_avg(ro, GMs, gNe, a0, ellB) / (GMs / ro ** 2); goJ = g_avg(ro, GMs, gNe, a0, ellJ) / (GMs / ro ** 2)
print("   Oortsche Wolke (Kreisbahnen), γ mit l_bar / mit l_Jacobi:  " + ",  ".join(f"{r/AU:.0e} AE: {g:.4f} / {h:.4f}" for r, g, h in zip(ro, go_, goJ)))

# ============================================================================ 6
H("6   Groesse der kovarianten Korrekturen (Sonne, l = l_bar)")
vcmb = 369.8e3                                             # [G] Geschwindigkeit gegen das CMB-Ruhesystem
b2 = (vcmb / c) ** 2
print(f"   Vorzugssystem: Kern im Ruhesystem des Koerpers, Gradient im Ruhesystem des Mediums; v = 369,8 km/s [G], (v/c)² = {b2:.2e}")
print(f"      anisotroper Anteil des Hooke-Zusatzes ~ (v/c)²·κ = {b2*kap:.2e} s⁻²,  Korrektur an Q2 ~ (v/c)²·Q2 = {b2*Q2B:.2e} s⁻²")
print(f"      Projektion der Gezeit mit U statt u: relative Aenderung ~ (v/c)² = {b2:.1e}  ->  δl/l ~ {b2/3:.1e}")
print(f"   Metrische Laenge im Kern: Potentialtiefe (v0/c)² = {(v0/c)**2:.2e};  eigenes Feld am Zonenrand GM/(c² l) = {GMs/(c**2*ellB):.2e}")
print(f"   Masse: σ-Fluss statt Komar-Fluss an der Sonnenoberflaeche: 2,1e-6 (-> δl/l = 7e-7); Komar-Fluss: exakt")
print(f"   Instantanes σ̄ auf dem Ruheschnitt: Lichtlaufzeit ueber die Zone l/c = {ellB/c/YR:.2f} Jahre; Retardierungsfehler ~ (v/c)² = {b2:.1e} relativ")
rho_ph = 0.01 * MSUN / PC ** 3                             # [G]
Lsc = R0
lnf = math.log(2 * Lsc / ellB) - 4 / 3
dg = 2 * math.pi * G * ellB ** 2 * lnf * rho_ph / Lsc
print(f"   η - 1 (Stern spuert ⟨Φ_s⟩_l, Licht und Gas Φ_s): Δg ≈ 2πG l² [ln(2L/l) - 4/3] ∇ρ_ph = {dg:.2e} m/s²,  Δg/g_obs = {dg/gobs:.1e}")
print(f"      (ρ_ph = 0,01 M☉/pc³ [G], Skala L = R0; Plummer-Kern: ⟨r²⟩ logarithmisch, abgeschnitten bei L)")
print(f"   Planeten in der Zone: δl/l = {ell/ellB-1:.1e} (Fixpunkt), Schwankung mit den Umlaeufen von Jupiter/Saturn in derselben Groesse")

# ============================================================================ 7
H("7   Zum Vergleich: Gezeit aus der vollen Riemann-Projektion -> Kopplung der Zonen an Gravitationswellen")
MU4, WMU4 = np.polynomial.legendre.leggauss(400)
XQ, WQ = np.polynomial.legendre.leggauss(12)
def V_coll(GM, prof, ell, gNe_, a0_, Nr=1500):
    """-(a0²/8πG)∫[Q_s(y)-Q_s(y_e)] d³x fuer das kollektive Feld g_e + eigenes geglaettetes Feld prof(r)"""
    r = np.logspace(math.log10(ell) - 5, math.log10(ell) + 5, Nr)
    gs = prof(r)[:, None]; m = MU4[None, :]
    ye = (gNe_ / a0_) ** 2; nue_ = float(nu(gNe_ / a0_))
    y = (gs ** 2 + gNe_ ** 2 + 2 * gs * gNe_ * m) / a0_ ** 2
    yp = ye + (y - ye)[..., None] * (1 + XQ) / 2
    dQ = (y - ye) / 2 * np.sum(WQ * (nu(np.sqrt(np.maximum(yp, 0))) - nue_), axis=-1)
    dens = (0.5 * np.sum(WMU4[None, :] * dQ, axis=1) + (nue_ - 1) * gs[:, 0] ** 2 / a0_ ** 2) * 4 * math.pi * r ** 2
    return -(a0_ ** 2 / (8 * math.pi * G)) * np.trapezoid(dens, r)
Ms = GMs / G
def Vpl(l_): return V_coll(GMs, lambda r: GMs * r / (r ** 2 + l_ ** 2) ** 1.5, ellB, gNe, a0)
dl = 0.02
V0, Vp, Vm = Vpl(ellB), Vpl(ellB * (1 + dl)), Vpl(ellB * (1 - dl))
Vl = (Vp - Vm) / (2 * dl * ellB); Vll = (Vp - 2 * V0 + Vm) / (dl * ellB) ** 2
VE = Vl * (-ellB / (3 * EB)); VEE = Vll * (ellB / (3 * EB)) ** 2 + Vl * 4 * ellB / (9 * EB ** 2)
print(f"   Strukturenergie der geglaetteten Sonne: V(l) = {V0:.4e} J,  dV/dl = {Vl:.4e} J/m,  d²V/dl² = {Vll:.4e} J/m²")
print(f"   l = (GM/|E|)^(1/3):  dV/d|E| = {VE:.3e} J s²,  d²V/d|E|² = {VEE:.3e} J s⁴")
print(f"   Welle mit Gezeit ½ω²h: Zusatz in der Wellen-Lagrangedichte -⅛ n V_EE ω² ḣ² Π(kl)²  ->  δ = 4πG n V_EE ω² Π²/c²,  n = 0,1 pc⁻³ [G]")
nst = 0.1 / PC ** 3
print(f"   {'f [Hz]':>10s} {'kl':>8s} {'Π(kl)':>10s} {'|δc_T/c| ≈ |δ|/2':>18s}")
for f in [1e-10, 1e-9, 3e-9, 1e-8, 3.17e-8, 1e-7, 1e-3, 100.0]:
    om = 2 * math.pi * f; kl = om / c * ellB; P = Pi(kl) if kl < 700 else 0.0
    dlt = 4 * math.pi * G * nst * VEE * om ** 2 * P ** 2 / c ** 2
    print(f"   {f:10.2e} {kl:8.3g} {P:10.3e} {abs(dlt)/2:18.3e}")
print("   Mit der Hesse-Form des Faltungspotentials σ̄ (gewaehlt) ist die Kopplung null: σ koppelt linear nicht an TT-Moden.")

# ============================================================================ 8
H("8   Schwellenuebergang (Koerper = Zusammenhangskomponente mit ρ_σ > ρ_J = 3√(3/2)|Ē^TF|/(4πG)) und l-Kraft")
rhoJ = 3 * EB / (4 * math.pi * G)
print(f"   ρ_J an der Sonne = {rhoJ:.3e} kg/m³ = {rhoJ/(MSUN/PC**3):.3f} M☉/pc³ = {rhoJ/1.6726e-27/1e6*0.73:.1f} H-Atome/cm³ (X = 0,73)")
print(f"   An der Schwelle gilt mittlere Dichte = ρ_J, also Radius R = (3M/4πρ_J)^(1/3) = l: das Objekt ist so gross wie seine Zone.")
Vpt = V_coll(GMs, lambda r: GMs / r ** 2, ellB, gNe, a0)
Vun = V_coll(GMs, lambda r: np.where(r < ellB, GMs * r / ellB ** 3, GMs / r ** 2), ellB, gNe, a0)
for nm, Vx in [("Punkt -> Plummer l", Vpt), ("homogene Kugel R = l -> Plummer l", Vun)]:
    dVx = V0 - Vx
    print(f"   Sprung {nm:36s}: ΔV = {dVx:.3e} J,  ΔV/M = {dVx/Ms:.3e} m²/s²,  √(2|ΔV|/M) = {math.sqrt(2*abs(dVx)/Ms):.2f} m/s")
Fl = -Vl / Ms * (2 / 3) * ellB / R0
print(f"   l-Kraft (l ∝ R^(2/3) bei flacher Rotationskurve): F/M = {Fl:.3e} m/s² = {Fl/gobs:.1e} g_obs")
