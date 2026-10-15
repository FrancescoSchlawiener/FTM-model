# -*- coding: utf-8 -*-
"""
FTM - tensorielle (kovariante) Feldgleichung.  Rechnungen zu forschung/06-tensorielle-feldgleichung.md
Aufruf:  python tensor.py > tensor.txt        (benoetigt nur numpy)

Abschnitte
  1  Kovariante Form: Kinetikfunktion Q~, NR-Grenzfall (numerische Pruefung), Metrik-Einheit von
     Faltung (g_rr) und Takt (g_tt), kovariante Integrodifferentialgleichung, Vakuumterm
  2  Klassische Tests in FTM-Groessen, q und Zusatzterm je Test, PPN, Quadrupol durch das
     aeussere Feld der Milchstrasse (Cassini)
  3  Linsenmasse gegen dynamische Masse: eta fuer Option A und B
"""
import math, sys
from fractions import Fraction
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
np.seterr(over="ignore", under="ignore", divide="ignore", invalid="ignore")

# ------------------------------------------------------------------ Konstanten
G    = 6.67430e-11
c    = 299792458.0
AU   = 1.495978707e11
PC   = 3.0856775814913673e16
KPC  = 1e3 * PC
MPC  = 1e6 * PC
LY   = 9.4607304725808e15
YR   = 3.15576e7                    # julianisches Jahr
CY   = 100 * YR
AS   = 180 / math.pi * 3600         # rad -> Bogensekunden
MAS  = AS * 1e3
H0   = 70e3 / MPC
RH   = c / H0
a0H  = c * H0 / (2 * math.pi)       # a0 = c H0 / 2pi  (FTM)
a0E  = 1.2e-10                      # empirischer RAR-Wert
A0S  = [("a0=cH0/2π", a0H), ("a0=1,2e-10", a0E)]
Msun = 1.98847e30
GMs  = 1.32712440018e20             # Sonne
GMe  = 3.986004418e14               # Erde
Rs   = 6.957e8
Re   = 6.3781e6                     # Erde, Aequator
Me   = 5.9722e24

def H(t): print("\n" + "=" * 90 + "\n" + t + "\n" + "=" * 90)
def ds(GM): return GM / (3 * c * c)                  # Faltungstiefe
def rstar(GM, a0): return math.sqrt(GM / a0)         # Schwellenradius

# ------------------------------------------------------------------ nu, q, Q~
def nu_m1_s(s):
    """nu-1 = sum_{n>=1} q^n = q/(1-q) = 1/(e^s-1),  s = sqrt(a_N/a0) = -ln q"""
    return 1.0 / np.expm1(s)
def nu_x(x):
    return 1.0 + nu_m1_s(np.sqrt(x))

# Bernoulli-Zahlen (B1 = -1/2), exakt
def bernoulli(N):
    B = [Fraction(1)]
    for m in range(1, N + 1):
        B.append(-sum(Fraction(math.comb(m + 1, k)) * B[k] for k in range(m)) / (m + 1))
    return [float(b) for b in B]
BN = bernoulli(40)
Q0 = -4 * math.pi ** 4 / 15                         # Q~(y=0) = -4 zeta(4) * 3! = -4 pi^4/15

def Qt_s(s):
    """Q~(y), y = (a_N/a0)^2 = s^4:  Q~ = -int_y^inf (nu-1) dy' = -int_s^inf 4 t^3/(e^t-1) dt"""
    s = np.asarray(s, dtype=float)
    out = np.empty_like(s)
    lo = s < 2.0
    sl = s[lo]
    acc = np.zeros_like(sl)
    for k in range(0, 40):
        acc += BN[k] * sl ** (k + 3) / (math.factorial(k) * (k + 3))
    out[lo] = Q0 + 4 * acc
    sh = s[~lo]
    acc = np.zeros_like(sh)
    for n in range(1, 80):
        acc += np.exp(-n * sh) * (sh ** 3 / n + 3 * sh ** 2 / n ** 2 + 6 * sh / n ** 3 + 6 / n ** 4)
    out[~lo] = -4 * acc
    return out

def Qt_polylog(s, N=200000):
    """dieselbe Funktion ueber Polylogarithmen: -4[ s^3 Li1(q) + 3 s^2 Li2(q) + 6 s Li3(q) + 6 Li4(q) ]"""
    q = math.exp(-s); n = np.arange(1, N + 1, dtype=float); qn = q ** n
    Li = lambda k: float(np.sum(qn / n ** k))
    return -4 * (s ** 3 * (-math.log1p(-q)) + 3 * s ** 2 * Li(2) + 6 * s * Li(3) + 6 * Li(4))

def Qt_int(s0):
    t = np.linspace(s0, 80.0, 800001)
    f = np.where(t > 0, 4 * t ** 3 / np.expm1(t), 0.0)
    return -np.trapezoid(f, t)

# ============================================================================ 1
H("1a  Kinetikfunktion des Strukturterms:  Q~(y) mit dQ~/dy = nu-1 = Σ q^n,  q = e^{-s},  s = y^{1/4} = sqrt(a_N/a0)")
print("   Q~(y) = -∫_y^∞ (nu-1) dy'  =  -4 [ s³Li₁(q) + 3s²Li₂(q) + 6sLi₃(q) + 6Li₄(q) ]   (Normierung Q~ -> 0 fuer q -> 0)")
print(f"   Q~(0)  = -4π⁴/15 = {Q0:.10f}")
print(f"   {'s':>6s} {'q':>11s} {'Q~ (Reihe)':>16s} {'Q~ (Polylog)':>16s} {'Q~ (Integral)':>16s} {'dQ~/dy num.':>13s} {'nu-1':>13s}")
for s in [0.05, 0.3, 1.0, 2.0, 3.0, 6.0, 12.0]:
    y = s ** 4; h = 1e-5 * y
    d = (Qt_s(np.array([(y + h) ** 0.25]))[0] - Qt_s(np.array([(y - h) ** 0.25]))[0]) / (2 * h)
    print(f"   {s:6.2f} {math.exp(-s):11.4e} {Qt_s(np.array([s]))[0]:16.10f} {Qt_polylog(s):16.10f} {Qt_int(s):16.10f} {d:13.6e} {float(nu_m1_s(s)):13.6e}")
print("   Grenzverhalten q -> 0:  Q~ ≈ -4q(s³+3s²+6s+6)  ->  0  (schneller als jede Potenz von 1/s);  nu-1 ≈ q")
print("   Tiefes Regime s -> 0:  nu ≈ 1/s = sqrt(a0/a_N)  ->  a = sqrt(a_N a0)  (Tully-Fisher)")

# ---------------------------------------------------------------------------- 1b
H("1b  Numerische Pruefung: statischer schwacher Grenzfall der Einstein-Hilbert-Wirkung")
print("   Metrik ds² = -(1+2εf)c²dt² + (1-2εh)dx⃗²,  f,h nur von x abhaengig.  Behauptung (Text):")
print("   √-g R = 4h'' + ε²[2h'² - 4f'h'] + Divergenzen + O(ε³)  =>  L_EH,NR = (1/8πG)[(∇Ψ)² - 2∇Φ·∇Ψ]")

def ricci_scalar_1d(g, x):
    N = len(x)
    gi = np.linalg.inv(g.transpose(2, 0, 1)).transpose(1, 2, 0)
    dg = np.zeros((4, 4, 4, N)); dg[1] = np.gradient(g, x, axis=2, edge_order=2)
    # Gam[l,m,n] = 1/2 g^{ls}(d_m g_sn + d_n g_sm - d_s g_mn)
    Gam = np.zeros((4, 4, 4, N))
    for l in range(4):
        for m in range(4):
            for n in range(4):
                acc = 0
                for s_ in range(4):
                    acc = acc + gi[l, s_] * (dg[m, s_, n] + dg[n, s_, m] - dg[s_, m, n])
                Gam[l, m, n] = 0.5 * acc
    dGam = np.zeros((4, 4, 4, 4, N)); dGam[1] = np.gradient(Gam, x, axis=3, edge_order=2)
    R = np.zeros(N)
    for m in range(4):
        for n in range(4):
            Rmn = 0
            for l in range(4):
                Rmn = Rmn + dGam[l, l, m, n] - dGam[n, l, l, m]
                for s_ in range(4):
                    Rmn = Rmn + Gam[l, l, s_] * Gam[s_, m, n] - Gam[l, n, s_] * Gam[s_, l, m]
            R = R + gi[m, n] * Rmn
    sq = np.sqrt(-np.linalg.det(g.transpose(2, 0, 1)))
    return R, sq, Gam

x = np.linspace(-14, 14, 56001)
f = np.exp(-(x - 0.7) ** 2 / 2.0)
hh = 0.6 * np.exp(-(x + 0.4) ** 2 / 3.0) + 0.3 * np.exp(-(x - 1.5) ** 2 / 0.8)
def metric(eps):
    g = np.zeros((4, 4, len(x)))
    g[0, 0] = -(1 + 2 * eps * f); g[1, 1] = g[2, 2] = g[3, 3] = 1 - 2 * eps * hh
    return g
Svals = {}
for eps in [1e-3, -1e-3]:
    R, sq, Gam = ricci_scalar_1d(metric(eps), x)
    Svals[eps] = np.trapezoid(sq * R, x)
    if eps > 0:
        a1 = Gam[0, 0, 1]                     # a_x = Γ^0_{01} fuer u = N^{-1} ∂_t
        lnN = 0.5 * np.log(1 + 2 * eps * f)
        err_a = np.max(np.abs(a1 - np.gradient(lnN, x))) / np.max(np.abs(np.gradient(lnN, x)))
S2 = (Svals[1e-3] + Svals[-1e-3]) / (2 * 1e-6)
fp, hp = np.gradient(f, x), np.gradient(hh, x)
ref = np.trapezoid(2 * hp ** 2 - 4 * fp * hp, x)
print(f"   Koeffizient ε² aus numerischem √-g R:  {S2:.6f}")
print(f"   ∫(2h'² - 4f'h') dx:                    {ref:.6f}     relative Abweichung {abs(S2/ref-1):.2e}")
print(f"   Aether-Beschleunigung a_x = Γ⁰₀₁ gegen ∂_x ln N (u = N⁻¹∂_t): max. rel. Abweichung {err_a:.2e}")
print("   => mit Φ=c²εf, Ψ=c²εh: L_EH,NR = (1/8πG)[(∇Ψ)² - 2∇Φ·∇Ψ];  c²|a| = |∇Φ| (Neigung des Mediums)")

# ---------------------------------------------------------------------------- 1c
H("1c  Variationspruefung der NR-Wirkung (G = a0 = M = 1, also r* = 1), Plummer-Kugel b = 0,2")
print("   L = (1/8πG)[(∇Ψ)² - 2∇Φ·∇Ψ + (∇Φ-∇Φ_N)² + a0² Q~(|∇Φ_N|²/a0²)] - ρΦ")
b = 0.2
r = np.linspace(1e-4, 60.0, 600001); dr = r[1] - r[0]
m = r ** 3 / (r ** 2 + b ** 2) ** 1.5
rho = 3 * b ** 2 / (4 * math.pi * (r ** 2 + b ** 2) ** 2.5)
gN = m / r ** 2                                      # dΦ_N/dr
nuv = 1 + nu_m1_s(np.sqrt(gN))
def integ(d): return np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * dr)])
PhiN = integ(gN)
confs = {"A: Φ'=Ψ'=nu·Φ_N'": (integ(nuv * gN), integ(nuv * gN), PhiN),
         "B: Φ'=nu·Φ_N', Ψ=Φ_N": (integ(nuv * gN), PhiN, PhiN),
         "GR: Φ=Ψ=Φ_N": (PhiN, PhiN, PhiN)}
def action(Phi, Psi, PN):
    rm = 0.5 * (r[1:] + r[:-1])
    d = lambda F: (F[1:] - F[:-1]) / dr
    P, S_, N_ = d(Phi), d(Psi), d(PN)
    L = (S_ ** 2 - 2 * P * S_ + (P - N_) ** 2 + Qt_s(np.sqrt(np.abs(N_)))) / (8 * math.pi)
    kin = np.sum(4 * math.pi * rm ** 2 * L) * dr
    src = np.trapezoid(4 * math.pi * r ** 2 * rho * Phi, r)
    return kin - src
print(f"   {'Konfiguration':24s} {'δS/δΦ':>11s} {'δS/δΨ':>11s} {'δS/δΦ_N':>11s}   (Bump bei r = 0,5 / 1 / 3; Werte skaliert)")
for name, (Phi, Psi, PN) in confs.items():
    row = []
    for which in range(3):
        vals = []
        for rc in [0.5, 1.0, 3.0]:
            bump = np.exp(-(r - rc) ** 2 / (0.15 * rc) ** 2)
            e = 1e-4
            F = [Phi, Psi, PN]
            Fp = [F[i] + (e * bump if i == which else 0) for i in range(3)]
            Fm = [F[i] - (e * bump if i == which else 0) for i in range(3)]
            dS = (action(*Fp) - action(*Fm)) / (2 * e)
            norm = np.trapezoid(4 * math.pi * r ** 2 * rho * bump, r) + np.trapezoid(r ** 2 * gN * np.abs(np.gradient(bump, r)), r)
            vals.append(abs(dS) / norm)
        row.append(max(vals))
    print(f"   {name:24s} {row[0]:11.2e} {row[1]:11.2e} {row[2]:11.2e}")
print("   => nur Konfiguration A macht die Wirkung stationaer: ∇²Φ_N = 4πGρ, Ψ = Φ, ∇²Φ = ∇·[nu ∇Φ_N]")

# ---------------------------------------------------------------------------- 1d
H("1d  FTM-Groessen im Strukturterm")
for nm, a0 in A0S:
    l0 = c * c / a0
    print(f"   {nm:11s}: a0 = {a0:.4e} m/s²   ℓ0 = c²/a0 = {l0:.4e} m = {l0/MPC:8.1f} Mpc   ℓ0/(2π R_H) = {l0/(2*math.pi*RH):.6f}")
GM = GMs; r_test = 5000 * AU
for nm, a0 in A0S:
    l0 = c * c / a0; aN = GM / r_test ** 2; Y = (aN / c ** 2) ** 2
    s_cov = (l0 ** 2 * Y) ** 0.25
    rs = math.sqrt(3 * c * c * ds(GM) / a0)
    print(f"   Sonne, r = 5000 AE, {nm}: (ℓ0²𝒴)^{{1/4}} = {s_cov:.6f}   r*/r = {rs/r_test:.6f}   τ = i·{s_cov/(2*math.pi):.6f}   q = e^{{2πiτ}} = {math.exp(-s_cov):.6f}")

# ---------------------------------------------------------------------------- 1e
H("1e  Faltung und Takt aus EINER Metrik: Δs aus g_rr, dτ/dt aus g_tt")
print("   (i) Excess-Radius einer homogenen Kugel aus g_rr = 1/(1-2Gm(r)/(c²r)), m(r) = M r³/R³:")
print("       δR = ∫0^R (√g_rr - 1) dr = R[arcsin(√x)/√x - 1],  x = 2GM/(c²R);   1. Ordnung: GM/(3c²) = Δs")
bodies = [("Erde", GMe, 6.371e6), ("Sonne", GMs, Rs), ("Sirius B", 1.018 * G * Msun, 5.8e6), ("Neutronenstern", 1.4 * G * Msun, 1.2e4)]
print(f"   {'Koerper':15s} {'Δs = GM/3c²':>14s} {'δR exakt':>14s} {'δR numerisch':>14s} {'δR/Δs':>10s}")
for nm, GM, R in bodies:
    xx = 2 * GM / (c * c * R)
    ex = R * (math.asin(math.sqrt(xx)) / math.sqrt(xx) - 1)
    rr = np.linspace(0, R, 200001)
    num = np.trapezoid(1 / np.sqrt(1 - xx * rr ** 2 / R ** 2) - 1, rr)
    print(f"   {nm:15s} {ds(GM):14.6e} {ex:14.6e} {num:14.6e} {ex/ds(GM):10.6f}")
print("   (ii) Vakuum-Einsteingleichungen, numerisch (Einheiten GM/c² = 1, Punkt r = 10, θ = 1):")

def metric_sph(kind):
    def gfun(X):
        t, rr_, th, ph = X
        g = np.zeros((4, 4))
        A = 1 - 2 / rr_
        if kind == "Schwarzschild": gt, gr = -A, 1 / A
        elif kind == "Takt ohne Faltung": gt, gr = -A, 1.0
        elif kind == "Takt quadriert": gt, gr = -A * A, 1 / A
        g[0, 0], g[1, 1] = gt, gr
        g[2, 2] = rr_ ** 2; g[3, 3] = (rr_ * math.sin(th)) ** 2
        return g
    return gfun
def christ_pt(gfun, X, h=1e-4):
    g = gfun(X); gi = np.linalg.inv(g); dg = np.zeros((4, 4, 4))
    for a in range(4):
        e = np.zeros(4); e[a] = h
        dg[a] = (gfun(X + e) - gfun(X - e)) / (2 * h)
    Gm = np.zeros((4, 4, 4))
    for l in range(4):
        for m_ in range(4):
            for n in range(4):
                Gm[l, m_, n] = 0.5 * sum(gi[l, s_] * (dg[m_, s_, n] + dg[n, s_, m_] - dg[s_, m_, n]) for s_ in range(4))
    return Gm
def einstein_pt(gfun, X, h=1e-3):
    Gm = christ_pt(gfun, X); dG = np.zeros((4, 4, 4, 4))
    for a in range(4):
        e = np.zeros(4); e[a] = h
        dG[a] = (christ_pt(gfun, X + e) - christ_pt(gfun, X - e)) / (2 * h)
    Ric = np.zeros((4, 4))
    for m_ in range(4):
        for n in range(4):
            Ric[m_, n] = sum(dG[l, l, m_, n] - dG[n, l, l, m_] for l in range(4)) + \
                sum(Gm[l, l, s_] * Gm[s_, m_, n] - Gm[l, n, s_] * Gm[s_, l, m_] for l in range(4) for s_ in range(4))
    g = gfun(X); gi = np.linalg.inv(g); Rs_ = np.sum(gi * Ric)
    Gt = Ric - 0.5 * g * Rs_
    return gi @ Gt                     # G^mu_nu
X0 = np.array([0.0, 10.0, 1.0, 0.3])
for kind in ["Schwarzschild", "Takt ohne Faltung", "Takt quadriert"]:
    Gmix = einstein_pt(metric_sph(kind), X0)
    gf = metric_sph(kind)(X0)
    print(f"       {kind:18s}: g_tt·g_rr = {gf[0,0]*gf[1,1]:+.6f}   G^t_t = {Gmix[0,0]:+.3e}   G^r_r = {Gmix[1,1]:+.3e}   G^θ_θ = {Gmix[2,2]:+.3e}")
print("       -> nur g_tt·g_rr = -1 loest die Vakuumgleichung; G^t_t - G^r_r = -(2/r) e^{-2β}(α'+β') erzwingt es.")
print("   (iii) Aussen: g_rr = 1/(1-6Δs/r) (aus G^t_t)  und  g_tt = -1/g_rr (aus G^t_t = G^r_r)  =>  dτ/dt = √(1-6Δs/r)")
for nm, GM, R in bodies[:2]:
    rr_ = 2 * R
    print(f"       {nm:6s} r = 2R: √(-g_tt) aus g_rr: {math.sqrt(1/(1/(1-6*ds(GM)/rr_))):.12f}   √(1-6Δs/r): {math.sqrt(1-6*ds(GM)/rr_):.12f}")

# ---------------------------------------------------------------------------- 1f
H("1f  Kovariante Integrodifferentialgleichung: Kurzgedaechtnis-Grenzfall auf der Kreisbahn")
print("   Kreisbahn (07, 1.3): a = a_N + κ0·M_k·a/(λ²+Ω²),  Ω² = a/r.  Mit M_k = q, κ0 = λ²:  λ -> ∞  =>  a = a_N/(1-q) = nu·a_N")
GMmw = 6e10 * G * Msun; a0 = a0E
rs_mw = rstar(GMmw, a0)
print(f"   Milchstrasse 6e10 M☉, a0 = 1,2e-10:  r* = {rs_mw/KPC:.3f} kpc")
print(f"   {'r/r*':>6s} {'nu':>9s}   a/a_N fuer λ/Ω_K = 1, 3, 10, 100, 1000")
for fr in [0.5, 1.0, 2.0, 5.0]:
    rr_ = fr * rs_mw; aN = GMmw / rr_ ** 2; q = math.exp(-math.sqrt(aN / a0)); OmK = math.sqrt(aN / rr_)
    out = []
    for k in [1, 3, 10, 100, 1000]:
        lam = k * OmK; lo, hi = aN, 1e3 * aN
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            F = mid - aN - lam ** 2 * q * mid / (lam ** 2 + mid / rr_)
            lo, hi = (mid, hi) if F < 0 else (lo, mid)
        out.append(mid / aN)
    print(f"   {fr:6.2f} {1/(1-q):9.5f}   " + "  ".join(f"{v:9.5f}" for v in out))

# ---------------------------------------------------------------------------- 1g
H("1g  Nebenbefund: Vakuumwert der Kinetikfunktion wirkt wie eine kosmologische Konstante")
OmL = 0.685
for nm, a0 in A0S:
    Lam_s = -Q0 * (a0 / c ** 2) ** 2
    Lam_o = 3 * OmL * H0 ** 2 / c ** 2
    print(f"   {nm:11s}: Λ_s = (4π⁴/15)(a0/c²)² = {Lam_s:.4e} m⁻²   Λ_obs(Ω_Λ={OmL}, H0=70) = {Lam_o:.4e} m⁻²   Verhaeltnis {Lam_s/Lam_o:.3f}")
print(f"   mit a0 = cH0/2π exakt: Λ_s = (π²/15)·H0²/c²,  Λ_s/Λ_obs = π²/(45 Ω_Λ) = {math.pi**2/(45*OmL):.3f}")

# ============================================================================ 2
H("2a  Klassische Tests in FTM-Groessen (γ = β = 1)")
dS = ds(GMs); dE = ds(GMe)
print(f"   Δs☉ = GM☉/3c² = {dS:.4f} m      Δs⊕ = {dE*1e3:.4f} mm")
alpha = 12 * dS / Rs
print(f"   Lichtablenkung am Sonnenrand: α = 12Δs/R☉ = {alpha:.6e} rad = {alpha*AS:.4f}″   (4GM/c²b = {4*GMs/(c*c*Rs)*AS:.4f}″)")
aM, eM, PM = 5.7909050e10, 0.205630, 87.9691 * 86400
dphi = 18 * math.pi * dS / (aM * (1 - eM ** 2))
print(f"   Merkur: Δφ = 18πΔs/(a(1-e²)) = {dphi:.6e} rad/Umlauf = {dphi*AS*CY/PM:.3f}″/Jh   (6πGM/(c²a(1-e²)): {6*math.pi*GMs/(c*c*aM*(1-eM**2))*AS*CY/PM:.3f}″/Jh)")
b_c, rE_, rC_ = 1.6 * Rs, AU, 7.4 * AU
dt_sh = 12 * dS / c * math.log(4 * rE_ * rC_ / b_c ** 2)
print(f"   Shapiro (Cassini-Geometrie b = 1,6 R☉, r⊕ = 1 AE, r_C = 7,4 AE): Δt_hin+rück = (12Δs/c)·((1+γ)/2)·ln(4r⊕r_C/b²) = {dt_sh*1e6:.1f} μs")
print("   Messung Cassini (Bertotti, Iess, Tortora 2003): γ - 1 = (2,1 ± 2,3)·10⁻⁵")
# Rotation / Drehfaltung
C_e, w_e = 8.034e37, 7.2921159e-5
J_e = C_e * w_e
dSJ_e = G * J_e / (3 * c ** 3)
aK_e = J_e / (Me * c)
J_s = 1.92e41
dSJ_s = G * J_s / (3 * c ** 3)
print("\n   Drehfaltung Δs_J = GJ/(3c³) = Δs·a_K,  a_K = J/(Mc) (Kerr-Laenge); gravitomagnetisch: g_0i = -6Δs_J (n×x)_i/r³")
print(f"   Erde : J = {J_e:.4e} kg m²/s   a_K = {aK_e:.4f} m   Δs_J = {dSJ_e:.4e} m²   (Δs⊕·a_K = {dE*aK_e:.4e} m²)")
print(f"   Sonne: J = {J_s:.3e} kg m²/s   a_K = {J_s/(Msun*c):.2f} m   Δs_J = {dSJ_s:.4e} m²")
aGPB, decIM = 7027.4e3, math.radians(16.84)
Om_orb = math.sqrt(GMe / aGPB ** 3)
Om_geo = 4.5 * dE / aGPB * Om_orb
Om_LT = 1.5 * c * dSJ_e / aGPB ** 3
print(f"   GP-B (a = 7027,4 km): geodaetisch Ω = (9/2)(Δs⊕/a)Ω_Bahn = {Om_geo*YR*MAS:.1f} mas/Jahr    Messung 6601,8 ± 18,3")
print(f"                         Frame-Dragging ⟨Ω⟩ = (3/2)cΔs_J/a³ = {Om_LT*YR*MAS:.2f} mas/Jahr; projiziert (IM Peg, δ = 16,84°): {Om_LT*math.cos(decIM)*YR*MAS:.2f} mas/Jahr    Messung 37,2 ± 7,2")
for nm, a_, e_ in [("LAGEOS", 12270e3, 0.0045), ("LAGEOS 2", 12163e3, 0.014)]:
    Om_node = 6 * c * dSJ_e / (a_ ** 3 * (1 - e_ ** 2) ** 1.5)
    print(f"   {nm:9s} Knoten: dΩ/dt = 6cΔs_J/(a³(1-e²)^1,5) = {Om_node*YR*MAS:.2f} mas/Jahr")

H("2b  Wert von q und Groesse des Zusatzterms je Test (isolierte Quelle)")
tests = [("Lichtablenkung, b = R☉", GMs, Rs, 1e-4),
         ("Merkur-Perihel, r = a", GMs, aM, 1e-4),
         ("Cassini, r bis 7,4 AE", GMs, 7.4 * AU, 2.3e-5),
         ("Erdbahn 1 AE", GMs, AU, None),
         ("Neptun 30,07 AE", GMs, 30.07 * AU, None),
         ("GP-B, r = 7027 km", GMe, aGPB, 0.19),
         ("LAGEOS, r = 12270 km", GMe, 12270e3, 0.05),
         ("Sonne bei r*", GMs, None, None)]
print(f"   {'Test':26s} {'a_N [m/s²]':>11s} {'s=√(a_N/a0)':>13s} {'log10 q (cH0/2π)':>17s} {'log10 q (1,2e-10)':>18s} {'rel. Zusatz ≈ q':>16s} {'Messgenauigkeit':>15s}")
for nm, GM, rr_, prec in tests:
    if rr_ is None: rr_ = rstar(GM, a0E)
    aN = GM / rr_ ** 2
    sH, sE = math.sqrt(aN / a0H), math.sqrt(aN / a0E)
    lqH, lqE = -sH / math.log(10), -sE / math.log(10)
    rel = f"10^{lqE:.0f}" if lqE < -3 else f"{math.exp(-sE)/(1-math.exp(-sE)):.3f}"
    pr = f"{prec:.1e}" if prec else "-"
    print(f"   {nm:26s} {aN:11.3e} {sE:13.4e} {lqH:17.1f} {lqE:18.1f} {rel:>16s} {pr:>15s}")
print("   (rel. Zusatz = (nu-1) = Zusatzbeschleunigung/Newton; konservativ mit dem groesseren a0 = 1,2e-10)")
print("   PPN: Strukturterm ∝ exp(-c·√(U/(a0 r))) hat keine Entwicklung in 1/c²  =>  γ = β = 1 in jeder PN-Ordnung")

H("2c  Quadrupol durch das aeussere Feld der Milchstrasse (QUMOND-Grenzfall der kovarianten Fassung)")
print("   Φ_Q = -(Q2/2) x^i x^j (e_i e_j - δ_ij/3),  Q2 = -(3/2) ∂z²Φ_s(0),  ∂z²Φ_s(0) = (1/4π)∫ (nu-1)∇Φ_N·∇[∂z²(1/r)] d³x")
print("   Einheiten r_e = √(GM☉/g_e), g_e: I(β) = ½∫∫ (nu(β|ĝ|)-nu(β)) [-3(3c²-1)/ρ⁴ + (9c-15c³)/ρ²] dρ dc,  β = g_e/a0")
def I_quad(beta, Nn=400, nuf=None):
    if nuf is None: nuf = lambda xx: 1 + nu_m1_s(np.sqrt(xx))
    xg, wg = np.polynomial.legendre.leggauss(Nn)
    useg = [(-7.0, -2.0), (-2.0, -0.3), (-0.3, 0.0), (0.0, 0.3), (0.3, 2.0), (2.0, 9.5)]
    tseg = [(0.0, 0.5 * math.pi), (0.5 * math.pi, 0.9 * math.pi), (0.9 * math.pi, math.pi)]
    nue = float(nuf(np.array(beta)))
    tot = 0.0
    for (u0, u1) in useg:
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0); wu = 0.5 * (u1 - u0) * wg
        rho_ = np.exp(u)
        for (t0, t1) in tseg:
            th = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0); wt = 0.5 * (t1 - t0) * wg
            R_, TH = np.meshgrid(rho_, th, indexing="ij")
            cc = np.cos(TH)
            gm = np.sqrt(np.maximum(1 / R_ ** 4 + 2 * cc / R_ ** 2 + 1, 0.0))
            dnu = nuf(beta * gm) - nue
            kern = -3 * (3 * cc ** 2 - 1) / R_ ** 4 + (9 * cc - 15 * cc ** 3) / R_ ** 2
            F = 0.5 * dnu * kern * R_ * np.sin(TH)
            tot += np.einsum("i,j,ij->", wu, wt, F)
    return tot
def I_direct(beta, Nu=2000, Nt=2000):
    """unabhaengige Methode: Phantomdichte S = ∇·[(nu-1)∇Φ_N] per Differenzen, dann -1/(4π)∫ S ∂z²(1/r) d³x"""
    u = np.linspace(-6.0, 9.0, Nu); th = np.linspace(1e-6, math.pi - 1e-6, Nt)
    R_, TH = np.meshgrid(np.exp(u), th, indexing="ij"); cc, sn = np.cos(TH), np.sin(TH)
    gr = 1 / R_ ** 2 + cc; gt = -sn                        # ĝ = r̂/ρ² + ẑ
    gm = np.sqrt(np.maximum(gr ** 2 + gt ** 2, 1e-300))
    nue = 1 + float(nu_m1_s(math.sqrt(beta)))
    dn = nu_m1_s(np.sqrt(beta * gm)) + 1 - nue                # nu - nu_e (Konstante traegt nicht bei)
    Vr, Vt = dn * gr, dn * gt
    dVr = np.gradient(R_ ** 2 * Vr, u, axis=0) / R_           # ∂ρ(ρ²Vr) = (1/ρ)∂u(...)
    S = dVr / R_ ** 2 + np.gradient(sn * Vt, th, axis=1) / (R_ * sn)
    # nu_e ist konstant: ∇·(nu_e ĝ) = 0 ausser im Ursprung; dort ist nu-1 = 0 => Korrekturterm (nu_e-1)·∇·ĝ entfaellt
    K = (3 * cc ** 2 - 1) / R_ ** 3                          # ∂z²(1/r)
    integrand = S * K * R_ ** 3 * sn                         # d³x = 2π ρ³ sinθ du dθ
    return -(1 / (4 * math.pi)) * 2 * math.pi * np.trapezoid(np.trapezoid(integrand, th, axis=1), u)
print("   Konvergenz (β = 1):  " + "  ".join(f"N={n}: {I_quad(1.0, n):.6f}" for n in [100, 200, 400, 800]))
print(f"   unabhaengige Methode (Phantomdichte, Gitter): I(1) = {I_direct(1.0):.6f}")
print(f"   {'g_e [m/s²]':>11s} {'a0':>11s} {'β=g_e/a0':>9s} {'I(β)':>10s} {'r_e [AE]':>9s} {'Q2 [s⁻²]':>11s} {'Q2/σ_Cassini':>13s}")
Q2list = []
for ge in [0.8e-10, 1.0e-10, 1.2e-10, 1.6e-10, 2.0e-10]:
    for nm, a0 in A0S:
        beta = ge / a0; I = I_quad(beta, 400); re = math.sqrt(GMs / ge)
        Q2 = -1.5 * (ge / re) * I
        Q2list.append(Q2)
        print(f"   {ge:11.2e} {nm:>11s} {beta:9.3f} {I:10.6f} {re/AU:9.0f} {Q2:11.3e} {(Q2-3e-27)/3e-27:13.1f}")
print("   Cassini-Ranging (Hees et al. 2014): Q2 = (3 ± 3)·10⁻²⁷ s⁻²")

# ============================================================================ 3
H("3   Linsenmasse gegen dynamische Masse: η = M_dyn/M_lens")
print("   M_dyn = r²|∇Φ|/G,  M_lens = r²|∇(Φ+Ψ)/2|/G.  Allgemein mit Strukturschlupf γ_s = Ψ_s/Φ_s:")
print("   η = nu / [1 + (1+γ_s)(nu-1)/2];   A: γ_s = 1 -> η = 1;   B: γ_s = 0 -> η = 2nu/(1+nu);   konform: γ_s = -1 -> η = nu")
print(f"   {'a_N/a0':>8s} {'nu':>8s} {'η_A':>5s} {'η_B':>7s} {'log10(g_lens/g_RAR) B':>22s} {'η_konform':>10s}")
for xv in [100, 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 1e-3, 1e-4]:
    n_ = float(nu_x(xv)); eB = 2 * n_ / (1 + n_)
    print(f"   {xv:8.0e} {n_:8.3f} {1:5d} {eB:7.3f} {math.log10((1+n_)/(2*n_)):22.3f} {n_:10.3f}")
print("\n   Isolierte Galaxien (Punktmasse, a0 = 1,2e-10):")
print(f"   {'M_b [M☉]':>9s} {'r [kpc]':>8s} {'g_bar [m/s²]':>13s} {'nu':>8s} {'η_B':>6s}")
for Mb in [1e10, 6e10, 1e11]:
    for rk in [30, 100, 300, 1000]:
        aN = G * Mb * Msun / (rk * KPC) ** 2; n_ = float(nu_x(aN / a0E))
        print(f"   {Mb:9.0e} {rk:8d} {aN:13.3e} {n_:8.2f} {2*n_/(1+n_):6.3f}")
print("\n   Linsen-Tully-Fisher (tief): v_lens⁴ = G M a0 (A)  bzw.  G M a0/4 (B)  ->  bei festem v: M_B/M_A = 4, Δlog M = "
      f"{math.log10(4):.3f} dex")
AK, nK, AE_, nE = 1006.0, 6.17, 14000.0, 30.0
print(f"   Statistik-Skalierung Euclid/KiDS-1000 (Flaeche·Quellendichte): √({AE_:.0f}·{nE:.0f} / ({AK:.0f}·{nK})) = {math.sqrt(AE_*nE/(AK*nK)):.1f}-fach kleinere Fehler")

H("2d  Welche Schaerfe des Uebergangs verlangt Cassini?  Familie nu_n(x) = (1 - e^{-x^n})^{-1/(2n)}  (n = 1/2: FTM-Rueckkopplung)")
print("   g_e = 1,2e-10 m/s², a0 = 1,2e-10 (β = 1)")
print(f"   {'n':>5s} {'nu(0,1)':>8s} {'nu(1)':>7s} {'nu(10)':>8s} {'I(1)':>10s} {'Q2 [s⁻²]':>11s} {'Q2/σ_Cassini':>13s}")
for n_ in [0.5, 1.0, 2.0, 3.0]:
    nuf = (lambda nn: (lambda xx: (-np.expm1(-np.power(xx, nn))) ** (-1 / (2 * nn))))(n_)
    I = I_quad(1.0, 400, nuf); re = math.sqrt(GMs / 1.2e-10); Q2 = -1.5 * (1.2e-10 / re) * I
    print(f"   {n_:5.1f} {float(nuf(np.array(0.1))):8.3f} {float(nuf(np.array(1.0))):7.3f} {float(nuf(np.array(10.0))):8.4f} {I:10.6f} {Q2:11.3e} {(Q2-3e-27)/3e-27:13.1f}")
R0 = 8.2 * KPC
print(f"   Zum Vergleich: reiner Gezeitenanteil des galaktischen Feldes g_e/R0 (R0 = 8,2 kpc, g_e = 1,2e-10) = {1.2e-10/R0:.2e} s⁻²")
print("   (so klein waere Q2, wenn nur Gezeiten und nicht das homogene aeussere Feld in nu eingingen)")
