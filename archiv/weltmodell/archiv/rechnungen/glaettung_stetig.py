# -*- coding: utf-8 -*-
"""
Stetiger Einsatz der Glättung in der Hierarchie der Trichter.
Regel (stetig): geglättet wird nur der Überschuss der Dichte über die Gezeitendichte ρ_J = 3T/(4πG),
    ρ̄ = min(ρ, ρ_J) + Σ_i M_i^c·K_{ℓ_i^c},   M_i^c = ∫_{V_i}(ρ − ρ_J) dV,   (ℓ_i^c)³ = G·M_i^c / T_i,
T = √(3/2)|⟨Ē^TF⟩| (Elterngezeit), K = Plummer-Kern. Punktmasse: M^c = M, also die bisherige Hierarchie.
Abschnitte
  1  Einsatz: Zentraldichte der Plummer-Wolke = ρ_J; Überschuss eines Plummer-Körpers; Energie beim Einsetzen
  2  Sonne und weite Doppelsterne (SI)
  3  Sternsysteme: Kugelsternhaufen, Crater II, Draco, Fornax, DF2, DF4 (Jeans, Feldform im äußeren Feld)
  4  Konsequente Stufendichte: Milchstraße als Körper in der Lokalen Gruppe
  5  Kovariante Nebenbedingung: Sprung beim Verschmelzen zweier Komponenten
[A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.  a₀ = 1,2·10⁻¹⁰ m/s².
"""
import io, os, sys, math
import numpy as np
from scipy.special import k1
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(all="ignore")
out = []
def P(s=""): out.append(s)
def H(t): P(); P("=" * 100); P(t); P("=" * 100)

# ------------------------------------------------------------------ gemeinsame Funktionen
def nu(x):
    """ν(g/a₀) = 1/(1 − e^{−s}), s = √(g/a₀)"""
    x = np.maximum(np.asarray(x, float), 1e-300); return 1.0 / (-np.expm1(-np.sqrt(x)))
def dlnnu(x):
    """d ln ν / d ln g"""
    return -(nu(x) - 1) * np.sqrt(x) / 2
MU, WMU = np.polynomial.legendre.leggauss(400)
def phantom(g, ge):
    """Monopol der Strukturbeschleunigung bei innerem Newton-Feld g (radial) und homogenem äußeren Feld ge
       (Flussmittel über die Kugel = Monopol der Feldform, Satz von Gauß). Einheit wie g, a₀ = 1."""
    g = np.atleast_1d(np.asarray(g, float)); gs = g[:, None]; m = MU[None, :]
    gm = np.sqrt(np.maximum(gs**2 + ge**2 + 2 * gs * ge * m, 0))
    return 0.5 * np.sum(WMU * (nu(gm) - 1) * (gs + ge * m), axis=1)

# Fourier-Transformierte des Plummer-Kerns: K̂(κ) = κ·K1(κ), κ = kℓ
KAP = np.linspace(1e-6, 80, 400001); KH = KAP * k1(KAP)
def FT_enclosed(u, n):
    """M(<u·ℓ)/M für die n-fache Plummer-Faltung (Skala ℓ)"""
    return (2 / np.pi) * np.trapezoid(KH**n * (np.sin(KAP * u) - KAP * u * np.cos(KAP * u)) / KAP, KAP)
def FT_density(u, n):
    """ρ(u·ℓ)·ℓ³/M für die n-fache Plummer-Faltung"""
    return (1 / (2 * np.pi**2)) * np.trapezoid(KH**n * KAP**2 * np.sinc(KAP * u / np.pi), KAP)

# ============================================================================ 1
H("1   Einsatz der Glättung")
P("   a) Identität: Plummer-Wolke der Masse M mit ℓ³ = GM/T hat die Zentraldichte 3M/(4πℓ³) = 3T/(4πG) = ρ_J.")
P("      Die Hierarchie senkt jeden Körper auf die Gezeitendichte ab; ein Körper, dessen größte Dichte schon ≤ ρ_J ist,")
P("      hätte nichts abzugeben. Stetige Regel: nur der Überschuss (ρ − ρ_J)_+ wird geglättet.")
P("   Probe Fourier-Kern: Plummer M(<ℓ)/M = " + f"{FT_enclosed(1.0, 1):.6f} (exakt 2^(-3/2) = {2**-1.5:.6f}),"
  + f"  ρ(0)ℓ³/M = {FT_density(0.0, 1):.6f} (exakt 3/4π = {3/(4*np.pi):.6f})")
def excess_plummer(bl):
    """Überschussanteil M^c/M eines Plummer-Körpers mit Skala b bei b/ℓ = bl (ℓ³ = GM/T mit der ganzen Masse)"""
    if bl >= 1: return 0.0
    x = math.sqrt(bl**(-1.2) - 1)
    return x**3 / (1 + x * x)**1.5 - (bl * x)**3
def jacobi_frac(bl):
    """Anteil innerhalb des Radius, in dem die mittlere umschlossene Dichte = ρ_J ist (umschlossene Form)"""
    return max(1 - bl * bl, 0.0)**1.5
# numerische Probe des Überschusses
def excess_num(bl):
    r = np.linspace(0, 20, 400001); rho = 3 / (4 * np.pi * bl**3) * (1 + (r / bl)**2)**-2.5; rJ = 3 / (4 * np.pi)
    return np.trapezoid(4 * np.pi * r**2 * np.maximum(rho - rJ, 0), r)
P("   b) Überschussanteil M^c/M eines Plummer-Körpers (Skala b), Einsatz bei b = ℓ (Zentraldichte = ρ_J):")
P("      b/ℓ        " + "  ".join(f"{x:6.2f}" for x in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.2)))
P("      Überschuss " + "  ".join(f"{(1.0 if x == 0 else excess_plummer(x)):6.3f}" for x in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.2)))
P("      numerisch  " + "  ".join(f"{(1.0 if x == 0 else excess_num(x)):6.3f}" for x in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.2)))
P("      umschlossen" + "  ".join(f"{jacobi_frac(x):6.3f}" for x in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.2)))
P("      Nahe am Einsatz: Überschuss ≈ " + f"{excess_plummer(0.99)/(0.01)**2.5:.2f}·(1 − b/ℓ)^(5/2) (bei b/ℓ = 0,99),"
  + f" umschlossen ≈ {jacobi_frac(0.99)/(0.01)**1.5:.2f}·(1 − b/ℓ)^(3/2)")
P("      Scharfe Ersetzung Körper → Plummer(ℓ): für ein Plummer-Profil bei b = ℓ identisch (kein Sprung), für andere Profile ein Sprung.")
P("   Poisson-Form (Tiefe w, e^{−|k|w}, P_a·P_b = P_{a+b}): Punktmasse in Tiefe w hat Zentraldichte M/(π²w³);")
P(f"      Tiefe bis zur Gezeitendichte w_J = (4/3π)^(1/3)·ℓ = {(4/(3*np.pi))**(1/3):.4f}·ℓ; ein Körper der Tiefe w_b erhält (w_J − w_b)_+,")
P("      Einsatz ebenfalls bei Zentraldichte = ρ_J. Die Kernform (Plummer in 𝒢, Poisson auf der Halbgeraden) legt das Modell nicht fest.")

# ---- Sonnenumgebung (SI), für Energie und Abschnitt 2
G = 6.67430e-11; c = 299792458.0; AU = 1.495978707e11; PC = 3.0856775814913673e16; KPC = 1e3 * PC
MSUN = 1.98847e30; GMs = 1.32712440018e20; a0 = float(__import__("os").environ.get("FTM_A0", "1.2e-10")); KK = 1e3 / KPC
def inv_rar(gob):
    lo, hi = 1e-8, 1e4
    for _ in range(200):
        m = math.sqrt(lo * hi)
        if m * float(nu(m)) < gob / a0: lo = m
        else: hi = m
    return math.sqrt(lo * hi) * a0
def tf(T): return T - np.trace(T) / 3 * np.eye(3)
def nrm(T): return math.sqrt(float(np.sum(T**2)))
SQ32 = math.sqrt(1.5)
v0, R0 = 233e3, 8.178 * KPC                                   # [G]
gobs = v0**2 / R0; gNe = inv_rar(gobs); be = gNe / a0; nue = float(nu(be)); Ke = float(dlnnu(be))
A_, B_ = 15.3 * KK, -11.9 * KK; rhob = 0.084 * MSUN / PC**3   # [G]
dlnv = -(A_ + B_) / (A_ - B_); dlngN = (2 * dlnv - 1) / (1 + Ke)
T_bar = np.diag([gNe / R0 * dlngN, gNe / R0, 4 * math.pi * G * rhob - gNe / R0 * dlngN - gNe / R0])
EB = SQ32 * nrm(tf(T_bar)); ellS = (GMs / EB)**(1 / 3); rhoJ = 3 * EB / (4 * math.pi * G)
ce = nue * (1 + Ke / 3) - 1
P(f"   c) Energie beim Einsetzen, Wolke 1 M☉ am Sonnenort: g_Ne = {be:.3f} a₀ [G], ν_e = {nue:.4f}, ℓ = {ellS/PC:.3f} pc,"
  f" ρ_J = {rhoJ/(MSUN/PC**3):.3f} M☉/pc³")
P(f"      Strukturenergie linear im äußeren Feld: V_s = −c_e·W, W = (1/8πG)∫|∇φ̄|², c_e = ν_e(1 + K/3) − 1 = {ce:.4f}")
# Probe der linearen Näherung gegen die volle Strukturenergie (Q_s) für Kugel R = ℓ → Plummer ℓ
XQ, WQ = np.polynomial.legendre.leggauss(12)
def Vs_full(prof, ell):
    r = np.logspace(math.log10(ell) - 5, math.log10(ell) + 5, 3000)
    gs = prof(r)[:, None]; m = MU[None, :]; ye = be**2
    y = (gs**2 + gNe**2 + 2 * gs * gNe * m) / a0**2
    yp = ye + (y - ye)[..., None] * (1 + XQ) / 2
    dQ = (y - ye) / 2 * np.sum(WQ * (nu(np.sqrt(np.maximum(yp, 0))) - nue), axis=-1)
    dens = (0.5 * np.sum(WMU[None, :] * dQ, axis=1) + (nue - 1) * gs[:, 0]**2 / a0**2) * 4 * math.pi * r**2
    return -(a0**2 / (8 * math.pi * G)) * np.trapezoid(dens, r)
def W_of(prof, ell):
    r = np.logspace(math.log10(ell) - 6, math.log10(ell) + 6, 20000)
    return np.trapezoid(prof(r)**2 * r**2, r) / (2 * G)
Ms_kg = GMs / G
hom = lambda R, GM: (lambda r: np.where(r < R, GM * r / R**3, GM / r**2))
plu = lambda l, GM: (lambda r: GM * r / (r**2 + l**2)**1.5)
dV_full = Vs_full(plu(ellS, GMs), ellS) - Vs_full(hom(ellS, GMs), ellS)
dV_lin = -ce * (W_of(plu(ellS, GMs), ellS) - W_of(hom(ellS, GMs), ellS))
P(f"      Probe Kugel R = ℓ → Plummer ℓ (scharfer Einsatz): ΔV/M voll = {dV_full/Ms_kg:.1f} m²/s², linear = {dV_lin/Ms_kg:.1f} m²/s²")
P("      Homogene Wolke (1 M☉), Dichte ρ₀ = (ℓ/R)³·ρ_J; ΔV/M gegen die ungeglättete Wolke [m²/s²]:")
P("      ρ₀/ρ_J        scharf (Körper → Plummer ℓ)    stetig (Überschuss → Plummer ℓ^c)")
for q in (0.99, 1.0 + 1e-9, 1.001, 1.01, 1.1, 1.5, 3.0, 10.0, 100.0):
    R = ellS / q**(1 / 3); Wh = W_of(hom(R, GMs), R)
    sharp = 0.0 if q <= 1 else -ce * (W_of(plu(ellS, GMs), ellS) - Wh)
    if q <= 1: cont = 0.0
    else:
        GMr = GMs / q; GMc = GMs - GMr; lc = (GMc / EB)**(1 / 3)
        cont = -ce * (W_of(lambda r: hom(R, GMr)(r) + plu(lc, GMc)(r), R) - Wh)
    P(f"      {q:10.4f}   {sharp/Ms_kg:14.1f}                  {cont/Ms_kg:14.3f}")

# ============================================================================ 2
H("2   Sonne und weite Doppelsterne (a₀ = 1,2·10⁻¹⁰, Sonnenumgebung [G])")
rho_sun = MSUN / (4 / 3 * math.pi * 6.957e8**3)
P(f"   Sonne: ρ_J/ρ̄☉ = {rhoJ/rho_sun:.2e}  ->  M^c/M − 1 = −{rhoJ/rho_sun:.1e}; ℓ = {ellS/PC:.4f} pc = {ellS/math.sqrt(GMs/a0):.2f} r* (wie bisher)")
re = math.sqrt(GMs / gNe)
def J_quad(beta, lam, Nn=240):
    xg, wg = np.polynomial.legendre.leggauss(Nn)
    bps = {-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 3.5, 6.0, 9.0, 12.0}
    L = math.log(lam); bps |= {L - 1.5, L - 0.7, L - 0.3, L, L + 0.3, L + 0.7, L + 1.5}
    bps = sorted(b for b in bps if -10.0 <= b <= 12.0)
    tseg = [(0.0, 0.5 * math.pi), (0.5 * math.pi, 0.9 * math.pi), (0.9 * math.pi, 0.99 * math.pi), (0.99 * math.pi, math.pi)]
    tot = 0.0
    for u0, u1 in zip(bps[:-1], bps[1:]):
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0); wu = 0.5 * (u1 - u0) * wg; rho_ = np.exp(u)
        for t0, t1 in tseg:
            th = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0); wt = 0.5 * (t1 - t0) * wg
            R_, TH = np.meshgrid(rho_, th, indexing="ij"); cc = np.cos(TH)
            h = R_ / (R_**2 + lam**2)**1.5
            gm = np.sqrt(np.maximum(h**2 + 2 * cc * h + 1, 0.0))
            kern = -3 * (3 * cc**2 - 1) * h / R_**2 + (9 * cc - 15 * cc**3) / R_**2
            F = 0.5 * (nu(beta * gm) - 1.0) * kern * R_ * np.sin(TH)
            tot += np.einsum("i,j,ij->", wu, wt, F)
    return tot + (4 / 15) * float(dlnnu(beta)) * float(nu(beta)) / lam**3
Q2 = -1.5 * (gNe / re) * J_quad(be, ellS / re)
P(f"   Quadrupol aus der Glättung Q₂ = {Q2:.2e} s⁻²  (Cassini-Grenze ~3·10⁻²⁷ s⁻² [G])")
MU2, WMU2 = np.polynomial.legendre.leggauss(2000)
def g_avg(r, GM, sm):
    r = np.atleast_1d(np.asarray(r, float)); gN = GM / r**2
    gs = (GM * r / (r**2 + sm**2)**1.5)[:, None]; m = MU2[None, :]
    gm = np.sqrt(np.maximum(gs**2 + gNe**2 + 2 * gs * gNe * m, 0.0))
    return gN + 0.5 * np.sum(WMU2[None, :] * (nu(gm / a0) - 1) * (gs + gNe * m), axis=1)
rr = np.array([1e3, 1e4, 3e4]) * AU
P("   Sonne γ − 1 = ⟨g⟩/g_N − 1:  " + ",  ".join(f"{r/AU:.0e} AE: {g-1:.1e}" for r, g in zip(rr, g_avg(rr, GMs, ellS) / (GMs / rr**2))))
# Doppelstern 2 × 0,75 M☉ [A]: jeder Stern ist Körper (M^c = M); Fixpunkt ℓ mit der geglätteten Gezeit des Begleiters
m1 = 0.75 * GMs
def comp_tide(d, l):
    """√(3/2)|TF| der Gezeit des über ℓ geglätteten Begleiters, gemittelt über die eigene Zone ℓ (zweifache Plummer-Faltung)"""
    u = d / l
    return abs(4 * math.pi * m1 * FT_density(u, 2) / l**3 - 3 * m1 * FT_enclosed(u, 2) / d**3)
P("   Weiter Doppelstern 2 × 0,75 M☉ [A]: Fixpunkte von ℓ³·√(3/2)|T_gal + T_Begleiter(ℓ)| = G·m (Achse längs R / φ / z)")
for dk in (5, 10, 30):
    d = dk * 1e3 * AU; res = []
    for ax in range(3):
        nn = np.zeros(3); nn[ax] = 1.0
        def h(l):
            Tc = comp_tide(d, l) * (np.outer(nn, nn) - np.eye(3) / 3)
            return l**3 * SQ32 * nrm(tf(T_bar) + Tc) - m1
        ls = np.logspace(math.log10(d) - 1.5, math.log10(ellS) + 0.5, 300); hv = [h(l) for l in ls]
        roots = [brentq(h, ls[i], ls[i + 1]) for i in range(len(ls) - 1) if hv[i] * hv[i + 1] < 0]
        res.append(roots)
    P(f"      d = {dk:2d} kAE:  " + ";  ".join("ℓ = " + ", ".join(f"{x/PC:.4f}" for x in rs) + " pc" for rs in res)
      + f"   (ohne Begleiter {(m1/EB)**(1/3)/PC:.4f} pc)")
Mb = 1.5 * GMs; lb = (Mb / EB)**(1 / 3); sv = np.array([5, 10, 20, 30]) * 1e3 * AU
P(f"   Paar als Ganzes (1,5 M☉, ℓ = {lb/PC:.3f} pc): γ − 1 = " + ",  ".join(f"{s/AU/1e3:.0f} kAE: {g-1:.1e}" for s, g in zip(sv, g_avg(sv, Mb, lb) / (Mb / sv**2))))

# ============================================================================ 3
H("3   Sternsysteme: Plummer (b = R_e [A]), isotrope Jeans-Gleichung, σ_los im Fenster R < R_e; Wirt als Punktmasse")
A0 = float(__import__("os").environ.get("FTM_A0", "1.2e-10")) * 3.0856775814913673e19 / 1e6                 # a₀ in (km/s)²/kpc
GK = 4.30091e-6
GMu = lambda M: GK * M / A0                                 # G·M in kpc²·a₀
kms = lambda s2: np.sqrt(np.maximum(s2, 0) * A0)
def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return kms(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r))
def tides(Mh, D):
    """Elterngezeit √(3/2)|E^TF| am Ort: baryonisch (Hesse von Φ̄_N, wie im Modell) und wirksam (Hesse des verstärkten Felds)"""
    GMh = GMu(Mh); Tb = 3 * GMh / D**3
    gD = lambda x: nu(GMh / x**2) * GMh / x**2; hh = 1e-4 * D
    gp = (gD(D + hh) - gD(D - hh)) / (2 * hh); Te = abs(gp - gD(D) / D)
    return Tb, float(Te)
def model(M, Re, Mh, D, T=None):
    """σ: Newton, Feldform (Sterne einzeln geglättet, System kollektiv), stetige Regel auf der Stufe des Systems (Elterngezeit T)"""
    b = Re; GM = GMu(M); r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, 5000)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5; Grho = 3 * GM / (4 * np.pi * b**3) * (1 + (r / b)**2)**-2.5
    ge = GMu(Mh) / D**2; gN = GMenc / r**2
    sN = jeans(r, Grho, gN, Re)
    sF = jeans(r, Grho, gN + phantom(gN, ge), Re)
    if T is None: return sN, sF
    GrJ = 3 * T / (4 * np.pi); Grc = np.maximum(Grho - GrJ, 0.0)
    GMc_enc = np.concatenate([[0.0], np.cumsum(0.5 * (4 * np.pi * r[1:]**2 * Grc[1:] + 4 * np.pi * r[:-1]**2 * Grc[:-1]) * np.diff(r))])
    GMc = GMc_enc[-1]
    if GMc <= 0: return sN, sF, sF, 0.0, 0.0
    lc = (GMc / T)**(1 / 3)
    gbar = (GMenc - GMc_enc) / r**2 + GMc * r / (r**2 + lc**2)**1.5
    a_s = phantom(gbar, ge)                                  # Monopol der Strukturbeschleunigung (einwärts) aus dem geglätteten Feld
    # Kraftgesetz aus der Energie: 𝒢†f = (1 − θ)f + θ·[⟨f⟩_K + (ℓ/3)∂_ℓ⟨f⟩_K],  θ = Stufe {ρ > ρ_J}
    Phis = -np.concatenate([np.cumsum((0.5 * (a_s[1:] + a_s[:-1]) * np.diff(r))[::-1])[::-1], [0.0]]) - a_s[-1] * r[-1]
    def zmean(l):
        K = 3 / (4 * np.pi * l**3) * (1 + (r / l)**2)**-2.5
        return np.trapezoid(4 * np.pi * r**2 * K * Phis, r)
    Phit = zmean(lc) + (lc / 3) * (zmean(lc * 1.001) - zmean(lc * 0.999)) / (0.002 * lc)
    theta = Grc > 0; ic = np.nonzero(theta)[0][-1]; rc = r[ic]
    Delta = np.interp(rc, r, Phis) - Phit                    # Potentialsprung am Rand von V (Zonenmittel innen, Φ_s außen)
    gv = gN + np.where(theta, 0.0, a_s)
    f = Grho * gv
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]]) + np.where(theta, GrJ * Delta, 0.0)
    w_ = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    sV = kms(np.trapezoid(r**2 * p * w_, r) / np.trapezoid(r**2 * Grho * w_, r))
    p0 = p - np.where(theta, GrJ * Delta, 0.0)
    sV0 = kms(np.trapezoid(r**2 * p0 * w_, r) / np.trapezoid(r**2 * Grho * w_, r))
    w = Grc / Grho                                           # Anteilsform: nur der Überschuss spürt das Zonenmittel
    sS = jeans(r, Grho, gN + (1 - w) * a_s, Re)
    STORE[(M, Re, T)] = dict(rc=rc, Delta=Delta, sV0=sV0, sS=sS)
    return sN, sF, sV, GMc / GM, lc
STORE = {}
SYS = [  # Name, M [M☉], R_e [kpc], Wirt [M☉], Abstand [kpc], σ gemessen, Fehler (+, −)
    ("kompakter KSH [A]", 2e5, 0.003, 6e10, 10.0, None, None),
    ("Pal 14 [G]", 2e4, 0.0247, 6e10, 71.0, 0.38, (0.12, 0.12)),
    ("NGC 2419 [G]", 9e5, 0.019, 6e10, 90.0, 4.1, (0.5, 0.5)),
    ("Crater II [G]", 3.2e5, 1.066, 6e10, 117.5, 2.7, (0.3, 0.3)),
    ("Draco [G]", 5.8e5, 0.221, 6e10, 76.0, 9.1, (1.2, 1.2)),
    ("Fornax [G]", 4.0e7, 0.71, 6e10, 147.0, 11.7, (0.9, 0.9)),
    ("DF2 [G]", 2.0e8, 2.2, 1e11, 80.0, 8.5, (2.3, 3.1)),
    ("DF4 [G]", 1.5e8, 1.6, 1e11, 165.0, 4.2, (4.4, 2.2))]
P("   Eingaben: Masse mit Υ = 2 [A] (Zwerge: L_V [G]), R_e [G], Wirt Milchstraße 6·10¹⁰ M☉ bzw. NGC 1052 10¹¹ M☉ [A], Abstand [G]")
P("   ρ̄/ρ_J = mittlere Dichte in r_½ gegen Gezeitendichte (baryonisch / wirksam); f = Überschussanteil M^c/M (baryonisch / wirksam)")
P(f"   {'System':18s} {'g_i/a₀':>8s} {'g_e/a₀':>8s} {'r½/r*':>6s} {'ρ̄/ρ_J':>13s} {'f':>11s} {'ℓ_b[kpc]':>8s} "
  f"{'Newton':>7s} {'Feldform':>8s} {'stetig b':>8s} {'stetig w':>8s}  {'gemessen':>14s}  {'σ²/σ_N² verlangt/Feldform':>26s}")
ROWS = {}
for nm, M, Re, Mh, D, so, er in SYS:
    GM = GMu(M); rh = Re / math.sqrt(2**(2 / 3) - 1); gi = GM / 2 / rh**2; ge = GMu(Mh) / D**2
    Tb, Te = tides(Mh, D)
    sN, sF, sSb, fb, lcb = model(M, Re, Mh, D, Tb)
    _, _, sSw, fw, lcw = model(M, Re, Mh, D, Te)
    rJb, rJw = GM / 2 / rh**3 / Tb, GM / 2 / rh**3 / Te
    ms = f"{so} +{er[0]}/−{er[1]}" if so else "–"
    ver = f"{(so/sN)**2:5.2f} / {(sF/sN)**2:5.2f}" if so else ""
    ROWS[nm] = dict(sN=sN, sF=sF, sSb=sSb, sSw=sSw, fb=fb, fw=fw, gi=gi, ge=ge, rJ=(rJb, rJw), rs=rh / math.sqrt(GM), lb=(GM / Tb)**(1 / 3))
    P(f"   {nm:18s} {gi:8.4f} {ge:8.4f} {rh/math.sqrt(GM):6.2f} {rJb:6.1f} / {rJw:5.1f} {fb:5.2f} / {fw:4.2f} {(GM/Tb)**(1/3):8.2f} "
      f"{sN:7.2f} {sF:8.2f} {sSb:8.2f} {sSw:8.2f}  {ms:>14s}  {ver:>26s}")
P("   Kraftgesetz der stetigen Regel (baryonische Gezeit): Rand r_c von V = {ρ > ρ_J}, Potentialsprung Δ = Φ_s(r_c) − Φ̃ dort;")
P("   σ mit Sprung (Tabelle) / ohne Sprung / Anteilsform (nur der Überschuss spürt das Zonenmittel):")
for nm, M, Re, Mh, D, so, er in SYS:
    Tb, _ = tides(Mh, D); s = STORE.get((M, Re, Tb))
    if s: P(f"      {nm:18s} r_c = {s['rc']:8.4f} kpc (r_c/R_e = {s['rc']/Re:5.2f}),  Δ = {s['Delta']*A0:9.2f} (km/s)²,  σ = {ROWS[nm]['sSb']:.2f} / {s['sV0']:.2f} / {s['sS']:.2f} km/s")
def sigma_c1(M, Re, Mh, D, T):
    """glatte Form ohne freie Zahl: geglättet wird E(ρ) = (ρ − ρ_J)²_+/ρ, Gewicht im Kraftgesetz E'(ρ) = 1 − (ρ_J/ρ)² (stetig)"""
    b = Re; GM = GMu(M); r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, 5000)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5; Grho = 3 * GM / (4 * np.pi * b**3) * (1 + (r / b)**2)**-2.5
    ge = GMu(Mh) / D**2; gN = GMenc / r**2; GrJ = 3 * T / (4 * np.pi)
    E = np.where(Grho > GrJ, (Grho - GrJ)**2 / Grho, 0.0); Ep = np.where(Grho > GrJ, 1 - (GrJ / Grho)**2, 0.0)
    GMe = np.concatenate([[0.0], np.cumsum(0.5 * (4 * np.pi * r[1:]**2 * E[1:] + 4 * np.pi * r[:-1]**2 * E[:-1]) * np.diff(r))])
    GMc = GMe[-1]; lc = (GMc / T)**(1 / 3)
    gbar = (GMenc - GMe) / r**2 + GMc * r / (r**2 + lc**2)**1.5; a_s = phantom(gbar, ge)
    Phis = -np.concatenate([np.cumsum((0.5 * (a_s[1:] + a_s[:-1]) * np.diff(r))[::-1])[::-1], [0.0]]) - a_s[-1] * r[-1]
    zm = lambda l: np.trapezoid(4 * np.pi * r**2 * 3 / (4 * np.pi * l**3) * (1 + (r / l)**2)**-2.5 * Phis, r)
    Phit = zm(lc) + (lc / 3) * (zm(lc * 1.001) - zm(lc * 0.999)) / (0.002 * lc)
    g = gN + (1 - Ep) * a_s + np.gradient(Ep, r) * (Phit - Phis)
    return jeans(r, Grho, g, Re), GMc / GM
P("   Glatte Form (ρ − ρ_J)²_+/ρ, Kraftgewicht 1 − (ρ_J/ρ)² (baryonische Gezeit): der Randterm wird zur verteilten Kraft ∂_r(Gewicht)·(Φ̃ − Φ_s):")
P("      " + ",  ".join(f"{nm.split(' [')[0]} σ = {sigma_c1(M, Re, Mh, D, tides(Mh, D)[0])[0]:.2f} km/s (f = {sigma_c1(M, Re, Mh, D, tides(Mh, D)[0])[1]:.2f})"
                         for nm, M, Re, Mh, D, so, er in SYS if nm.split(' [')[0] in ("Crater II", "Fornax", "DF2", "DF4")))
P("   Feldform = bisherige Regel (Sterne sind Körper, das System ist kollektiv) = stetige Regel mit der Dichte der Quelle (Lesart i).")
P("   stetig b / w = stetige Regel mit der Dichte der Stufe (Lesart ii, Sterne vorher geglättet), Elterngezeit baryonisch / wirksam.")
d4, fx = ROWS["DF4 [G]"], ROWS["Fornax [G]"]
P(f"   DF4 gegen Fornax: g_i/a₀ {d4['gi']:.4f} / {fx['gi']:.4f},  g_e/a₀ {d4['ge']:.4f} / {fx['ge']:.4f},  r½/r* {d4['rs']:.2f} / {fx['rs']:.2f},"
  f"  ρ̄/ρ_J {d4['rJ'][0]:.0f} / {fx['rJ'][0]:.0f},  ℓ/r½ {d4['lb']/(1.6/math.sqrt(2**(2/3)-1)):.1f} / {fx['lb']/(0.71/math.sqrt(2**(2/3)-1)):.1f}")
P("   Verlangte Richtung je System: (σ_gem/σ_Feldform)², aufsteigend nach der Kenngröße geordnet (< 1: weniger Verstärkung als die Feldform):")
MEAS = {nm: (so / ROWS[nm]["sF"])**2 for nm, M, Re, Mh, D, so, er in SYS if so}
for lab, fn in (("g_i/a₀", lambda d: d["gi"]), ("g_e/a₀", lambda d: d["ge"]), ("g_i/g_e", lambda d: d["gi"] / d["ge"]),
                ("ρ̄/ρ_J", lambda d: d["rJ"][0]), ("r½/r*", lambda d: d["rs"])):
    srt = sorted(MEAS, key=lambda n: fn(ROWS[n]))
    P(f"      {lab:8s}: " + "  ".join(f"{n.split(' [')[0]} ({fn(ROWS[n]):.3g}: {MEAS[n]:.2f})" for n in srt))

# ============================================================================ 4
H("4   Konsequente Stufendichte: auch Wirte sind Körper (Lesart ii ohne Schnitt)")
Mmw, bmw, M31, D31 = 6e10, 3.0, 1e11, 785.0                   # [A] Plummer b = 3 kpc; M31 baryonisch [A], Abstand [G]
Tlg = 3 * GMu(M31) / D31**3; lmw = (GMu(Mmw) / Tlg)**(1 / 3)
r = np.logspace(-2, 3.5, 4000); GM = GMu(Mmw)
GMenc = GM * r**3 / (r**2 + bmw**2)**1.5; Grho = 3 * GM / (4 * np.pi * bmw**3) * (1 + (r / bmw)**2)**-2.5; gN = GMenc / r**2
ge = GMu(M31) / D31**2
GrJ = 3 * Tlg / (4 * np.pi); Grc = np.maximum(Grho - GrJ, 0)
GMc_enc = np.concatenate([[0.0], np.cumsum(0.5 * (4 * np.pi * r[1:]**2 * Grc[1:] + 4 * np.pi * r[:-1]**2 * Grc[:-1]) * np.diff(r))])
GMc = GMc_enc[-1]; lc = (GMc / Tlg)**(1 / 3)
gbar = (GMenc - GMc_enc) / r**2 + GMc * r / (r**2 + lc**2)**1.5
gF = gN + phantom(gN, ge); gS = gN + np.where(Grc > 0, 0.0, phantom(gbar, ge))   # Kraftgesetz: in V nur das Zonenmittel (v_c ohne Randsprung)
rc_edge = r[np.argmax(Grc <= 0)]
P(f"   Milchstraße 6·10¹⁰ M☉ (Plummer b = 3 kpc [A]) im Feld von M31 (10¹¹ M☉ [A], 785 kpc [G]): ℓ = {lmw:.0f} kpc,"
  f" Überschussanteil {GMc/GM:.4f}, Überschuss bis r = {rc_edge:.0f} kpc")
P("   v_c [km/s]      " + "  ".join(f"{x:>6.0f} kpc" for x in (8, 20, 50, 100)))
for lab, g in (("Newton", gN), ("Feldform", gF), ("stetig (ii)", gS)):
    P(f"   {lab:14s}  " + "  ".join(f"{math.sqrt(np.interp(x, r, g)*x*A0):10.1f}" for x in (8, 20, 50, 100)))
P("   Alle Satelliten liegen in dieser Zone: σ fällt für Crater II, Draco, Fornax auf die Newton-Spalte.")

# ============================================================================ 5
H("5   Kovariante Nebenbedingung D²σ̄ = 𝒢[σ] mit Überschuss: Verschmelzen zweier Komponenten (Sonnenumgebung)")
P("   𝒢[σ] = D²σ − (4πG/c²)[(ρ_σ − ρ_J)_+ − Σ_i γ_i M_i^c 𝒦_{ℓ_i^c}(ρ_i)],  M_i^c = ∫_{V_i}(ρ_σ − ρ_J) dV,  V_i Komponente von {ρ_σ > ρ_J}")
I0 = (2 / np.pi) * np.trapezoid(KH**2, KAP)                  # ⟨1/|x−y|⟩·ℓ für zwei Plummer-Wolken gleicher Skala im selben Zentrum
P(f"   Probe: ∫(κK1)²dκ = {np.trapezoid(KH**2, KAP):.6f} (exakt 3π²/32 = {3*np.pi**2/32:.6f})")
m = GMs                                                       # je 1 M☉ Überschuss, kompakte Kerne (b ≪ ℓ) [A]
lm = (m / EB)**(1 / 3); l2 = (2 * m / EB)**(1 / 3)
W_before = 2 * (3 * np.pi / 32) * m**2 / lm / G + m**2 * I0 / lm / G
W_after = (3 * np.pi / 32) * (2 * m)**2 / l2 / G
dVm = -ce * (W_after - W_before) / (2 * Ms_kg)
P(f"   Zwei kompakte Kerne je 1 M☉ verschmelzen (Sattel erreicht ρ_J bei Abstand ≪ ℓ): vorher 2 × Plummer(ℓ = {lm/PC:.3f} pc),"
  f" nachher Plummer(ℓ = {l2/PC:.3f} pc)")
P(f"   Sprung ΔV/M = {dVm:.1f} m²/s² (zum Vergleich: scharfer Einsatz einer homogenen Wolke {dV_lin/Ms_kg:.1f} m²/s²)")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "glaettung_stetig.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
