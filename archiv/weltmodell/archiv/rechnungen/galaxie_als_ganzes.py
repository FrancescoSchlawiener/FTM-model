# -*- coding: utf-8 -*-
"""
Beschleunigung eines Systems als Ganzes im homogenen äußeren Newton-Feld g_e.
Jede Masse folgt der örtlichen Regel im Gesamtfeld: a(x) = ν(|g_i + g_e|)(g_i + g_e).
Das System als Ganzes: massengewichtetes Mittel A = (1/M)∫ρ·a dV (kein zusätzliches Gesetz).
Im Bezug des Systems: a_int(x) = a(x) − A, ∫ρ·a_int dV = 0 nach Konstruktion.
Abschnitte
  1  Plummer-Kugel: A/g_e gegen X = GM/(b²a₀) und g_e/a₀; Grenzfälle g_e und ν(g_e)g_e; Symmetrie A ∥ g_e (3D-Quadratur, Monte Carlo)
  2  Milchstraßen-ähnliche Masse (6·10¹⁰ M☉), g_e = 0,03·a₀ in der Ebene: A und Einseitigkeit m = 1 für Punktmasse, Hernquist, Exponentialscheibe
  3  Druckgestützte Systeme (Eingaben aus glaettung_stetig.py): Jeans mit a_int im Bezug A; Virial; Rand, an dem a_int nach außen zeigt
  4  Ordnung der Abweichungen nach A
  5  Unabhängige Prüfung: erste Ordnung g_i ≫ g_e (Beta-Funktion), zweite Ordnung g_i ≪ g_e, Monopol eines konstanten Vektors
Stellen, an denen Newton-Felder von Sternen vorkämen, sind geglättet: g_i ist das Feld der glatten Verteilung.
[A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.  a₀ = A0 aus dem Kern.
"""
import io, os, sys, math
import numpy as np
from scipy.special import i0e, i1e, k0e, k1e, gamma
from ftm_kern import G, M_sun, kpc, A0, nu

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(all="ignore")
out = []
def P(s=""): out.append(s)
def H(t): P(); P("=" * 100); P(t); P("=" * 100)

def nuu(x):
    """ν in Einheiten a₀ = 1 (x = |g|/a₀), sicher bei x = 0"""
    return nu(np.maximum(np.asarray(x, float), 1e-300), 1.0)

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def gl_split(edges, n):
    xs, ws = zip(*[gl(a, b, n) for a, b in zip(edges[:-1], edges[1:])])
    return np.concatenate(xs), np.concatenate(ws)

# ------------------------------------------------------------------ Schalenmittel (a₀ = 1)
# Kugelschale mit innerem Feld g (einwärts), g_e = ge·ẑ: S(g) = ½∫dμ ν(|g_tot|)(g_tot·ẑ).
# Mit |g_tot| = u²:  μ = (g² + ge² − u⁴)/(2g·ge),  g_tot·ẑ = (u⁴ + ge² − g²)/(2ge),  ν(u²)u² = u²/(1 − e^{−u}) glatt in u:
#   S(g) = ∫_{√|g−ge|}^{√(g+ge)} u³/(1 − e^{−u})·(u⁴ + ge² − g²) du / (2g·ge²)
#   für g < ge/2 direkt in μ; für g > 2ge stabiler mit |g_tot| = g + ge·τ:  S(g) = (1/4g)∫₋₁¹ ν(g_tot)g_tot(2gτ + ge τ² + ge) dτ
UX, UW = np.polynomial.legendre.leggauss(64)
def shell(g, ge):
    g = np.atleast_1d(np.asarray(g, float)); S = np.empty_like(g)
    k = g > 2 * ge
    if np.any(k):
        gk = g[k][:, None]; gt = gk + ge * UX[None, :]
        S[k] = (nuu(gt) * gt * (2 * gk * UX[None, :] + ge * UX[None, :]**2 + ge)) @ UW / (4 * g[k])
    j = g < 0.5 * ge                                          # g ≪ ge: μ-Form ohne Nullstelle von |g_tot|
    if np.any(j):
        gk = g[j][:, None]; gm = np.sqrt(gk**2 + ge**2 - 2 * gk * ge * UX[None, :])
        S[j] = 0.5 * (nuu(gm) * (ge - gk * UX[None, :])) @ UW
    k = ~k & ~j
    if np.any(k):
        gk = g[k]; u1, u2 = np.sqrt(np.abs(gk - ge)), np.sqrt(gk + ge)
        u = 0.5 * (u2 - u1)[:, None] * UX[None, :] + 0.5 * (u2 + u1)[:, None]
        f = u**3 / (-np.expm1(-u)) * (u**4 + ge**2 - gk[:, None]**2)
        S[k] = 0.5 * (u2 - u1) * (f @ UW) / (2 * gk * ge**2)
    return S
def cuts_of(gfun, ge):
    """Stellen y ∈ (0, 1) mit gfun(y) = ge (Knick des Integranden)"""
    from scipy.optimize import brentq
    y = np.unique(np.concatenate([np.logspace(-14, -0.3, 3000), 1 - np.logspace(-14, -0.3, 3000), np.linspace(0.01, 0.99, 3000)]))
    d = gfun(y) - ge; i = np.nonzero(np.sign(d[1:]) != np.sign(d[:-1]))[0]
    return [brentq(lambda x: gfun(x) - ge, y[j], y[j + 1], xtol=1e-16) for j in i]
def nodes(gfun, dmdy, ge, n=32):
    """Gauß-Knoten in y ∈ (0, 1), geteilt an den Knicken und geometrisch an beiden Enden"""
    e = np.concatenate([[0.0, 1.0], np.logspace(-10, -0.4, 25), 1 - np.logspace(-10, -0.4, 25), cuts_of(gfun, ge)])
    y, wt = gl_split(np.unique(e), n)
    return y, wt * dmdy(y)
# Plummer (b = 1): m = w³, g_i = X·w(1 − w²); Variable v mit w = 1 − v² (glatt am Rand): dm = 3w²·2v dv
def plummer_A(X, ge, n=32):
    """A/g_e für Plummer-Kugel, X = GM/(b²a₀)"""
    gf = lambda v: X * (1 - v * v) * (1 - (1 - v * v)**2)
    v, wt = nodes(gf, lambda v: 6 * v * (1 - v * v)**2, ge, n)
    return float(np.sum(wt * shell(gf(v), ge))) / ge

def plummer_A3d(X, ge, n_hat, nw=240, nth=240, nph=256):
    """volle 3D-Quadratur (w, θ, φ), g_e in beliebiger Richtung n_hat; gibt Vektor A/g_e"""
    v, vw = gl_split(np.unique(np.concatenate([[0.0, 1.0], np.logspace(-6, -0.3, 11)])), nw // 12)
    w = 1 - v * v; wm = 6 * v * w * w * vw
    th, wth = gl(0, math.pi, nth); ph = np.arange(nph) * 2 * math.pi / nph; wph = 2 * math.pi / nph
    g = X * w * (1 - w * w)
    er = np.stack([np.sin(th)[:, None] * np.cos(ph)[None, :], np.sin(th)[:, None] * np.sin(ph)[None, :],
                   np.cos(th)[:, None] * np.ones_like(ph)[None, :]])            # (3, nth, nph)
    wang = (wth * np.sin(th))[:, None] * wph / (4 * math.pi)
    gev = ge * np.asarray(n_hat, float)
    A = np.zeros(3)
    for gi, wi in zip(g, wm):
        gt = -gi * er + gev[:, None, None]
        A += wi * np.sum(nuu(np.linalg.norm(gt, axis=0)) * gt * wang, axis=(1, 2))
    return A / ge

def plummer_A_mc(X, ge, n_hat, N=2_000_000, seed=1):
    """Monte Carlo (anderer Weg): Stichprobe der Plummer-Masse (Paare ±x), g_e in Richtung n_hat; Mittel und Fehler je Komponente"""
    rng = np.random.default_rng(seed)
    m = rng.random(N); t = m**(2 / 3); g = X * m / (t / (1 - t))
    v = rng.normal(size=(N, 3)); er = v / np.linalg.norm(v, axis=1)[:, None]
    gev = ge * np.asarray(n_hat, float)[None, :]
    a = 0.0
    for sg in (1, -1):
        gt = -sg * g[:, None] * er + gev
        a = a + 0.5 * nuu(np.linalg.norm(gt, axis=1))[:, None] * gt
    return a.mean(0) / ge, a.std(0) / math.sqrt(N) / ge

# ============================================================================ 1
H("1   Plummer-Kugel (Masse M, Skala b) im homogenen Feld g_e: A/g_e,  X = GM/(b²a₀),  größtes inneres Feld 0,385·X·a₀")
GES = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
XS = [0.026, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 260.0]
P("   Zeilen X (g_i,max/a₀), Spalten g_e/a₀.  Eintrag A/g_e;  Grenzfälle: Newtonsches Zentrum A/g_e = 1, Feldform A/g_e = ν(g_e)")
P("   " + f"{'X':>7s} {'g_i,max':>8s} {'g_i(r½)':>8s} | " + " ".join(f"{x:>8.3f}" for x in GES))
P("   " + f"{'ν(g_e)':>25s} | " + " ".join(f"{float(nuu(x)):8.3f}" for x in GES))
rh_b = 1 / math.sqrt(2**(2 / 3) - 1)
TAB1 = {}
for X in XS:
    row = [plummer_A(X, ge) for ge in GES]; TAB1[X] = row
    P("   " + f"{X:7.3f} {0.3849*X:8.3f} {X*0.5/rh_b**2:8.4f} | " + " ".join(f"{v:8.3f}" for v in row))
P("   Anteil zwischen den Grenzfällen f = (A/g_e − 1)/(ν(g_e) − 1)  (0 = Newtonsches Zentrum, 1 = Feldform):")
for X in XS:
    P("   " + f"{X:7.3f} {'':8s} {'':8s} | " + " ".join(f"{(v-1)/(float(nuu(ge))-1):8.3f}" for v, ge in zip(TAB1[X], GES)))
nh = np.array([1.0, 2.0, 3.0]) / math.sqrt(14)
P("   Symmetrie: g_e in Richtung (1,2,3)/√14, volle 3D-Quadratur (r, θ, φ) und Monte Carlo (2·10⁶ Massenpunkte):")
for X, ge in ((0.1, 0.03), (1.0, 0.1), (10.0, 0.01)):
    A3 = plummer_A3d(X, ge, nh); par = A3 @ nh; perp = np.linalg.norm(A3 - par * nh)
    Amc, emc = plummer_A_mc(X, ge, nh); pmc = Amc @ nh; qmc = np.linalg.norm(Amc - pmc * nh); emc_par = math.sqrt(np.sum((emc * nh)**2))
    P(f"      X = {X:5.2f}, g_e = {ge:5.2f}: A∥/g_e = {par:.5f} (3D), {plummer_A(X, ge):.5f} (axial), {pmc:.4f} ± {emc_par:.4f} (MC);"
      f"  |A⊥|/|A| = {perp/abs(par):.1e} (3D), {qmc/abs(pmc):.1e} (MC, Fehler {np.linalg.norm(emc)/abs(pmc):.1e})")

# ============================================================================ 2
H("2   Milchstraßen-ähnliche Masse 6·10¹⁰ M☉, g_e = 0,03·a₀ in der Ebene: A und Einseitigkeit m = 1 von v_c(φ)")
Mmw = 6e10 * M_sun; GE2 = 0.03
XK = G * Mmw / (A0 * kpc**2)                                  # g/a₀ = XK/r² (r in kpc)
P(f"   Punktmasse als Grenzfall b → 0 der Plummer-Kugel (A/g_e):  " +
  ",  ".join(f"b = {b:g} kpc: {plummer_A(XK/b**2, GE2):.4f}" for b in (3.0, 1.0, 0.3, 0.1, 0.01)))
P(f"   (ν(g_e) = {float(nuu(GE2)):.3f}; bei b → 0 liegt die Masse im Newtonschen Zentrum, A → g_e)")
# Hernquist [A] a = 3 kpc: M(<r)/M = r²/(r+a)²  ->  r = a√m/(1 − √m)
def hernquist_A(Xk, a, ge, n=32):
    """Hernquist (a): m = p², g_i = (GM/a²)(1 − p)², dm = 2p dp"""
    X = Xk / a**2; gf = lambda p: X * (1 - p)**2
    p, wt = nodes(gf, lambda p: 2 * p, ge, n)
    return float(np.sum(wt * shell(gf(p), ge))) / ge
def hern_g(Xk, a, R): return Xk / (R + a)**2
# dünne Exponentialscheibe R_d = 3 kpc: g_R = GM R/(2R_d³)·[I0K0 − I1K1](R/2R_d)
def expd_g(Xk, Rd, R):
    y = np.asarray(R, float) / (2 * Rd)
    return Xk * R / (2 * Rd**3) * (i0e(y) * k0e(y) - i1e(y) * k1e(y))
RD = 3.0
def expd_A(Xk, Rd, ge, n=400, nph=720):
    # Masse in der Ebene: dM/M = (R/R_d²)e^{−R/R_d} dR dφ/2π;  u = R/R_d
    u, wu = gl_split(np.array([0, 0.25, 0.5, 1, 2, 3, 5, 8, 12, 20, 35, 60]), n // 10)
    ph = np.arange(nph) * 2 * math.pi / nph
    R = u * Rd; g = expd_g(Xk, Rd, R)
    gx = ge - g[:, None] * np.cos(ph)[None, :]; gy = -g[:, None] * np.sin(ph)[None, :]
    a = nuu(np.hypot(gx, gy))
    w = wu * u * np.exp(-u)
    Ax = np.einsum("i,ij->", w, a * gx) / nph; Ay = np.einsum("i,ij->", w, a * gy) / nph
    return Ax / ge, Ay / ge, np.sum(w)
Ah = hernquist_A(XK, 3.0, GE2); Axe, Aye, norm_e = expd_A(XK, RD, GE2)
P(f"   Hernquist a = 3 kpc [A]: A/g_e = {Ah:.4f};  Exponentialscheibe R_d = 3 kpc: A/g_e = {Axe:.4f} (A_y/g_e = {Aye:.1e}, Masse-Norm {norm_e:.6f})")
P(f"   Probe Scheibenfeld bei 50 kpc: g_R/(GM/R²) = {expd_g(XK, RD, 50.0)/(XK/2500):.4f}; zentrale Flächendichte M/(2πR_d²) = "
  f"{6e10/(2*math.pi*9e6):.0f} M☉/pc²")
phi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
def m1(gfun, R, ref):
    """m = 1-Amplitude [%] von v_c(φ) bei R (kpc), a_int = ν(|g|)g − ref, g_e = GE2 in x-Richtung"""
    gi = gfun(R); eR = np.stack([np.cos(phi), np.sin(phi)])
    gt = -gi * eR + np.array([GE2, 0.0])[:, None]
    a = nuu(np.linalg.norm(gt, axis=0)) * gt - np.array([ref, 0.0])[:, None]
    vc = np.sqrt(np.clip(-R * np.sum(a * eR, axis=0), 0, None))
    F = np.fft.rfft(vc) / len(vc); return 2 * abs(F[1]) / F[0].real * 100, F[0].real
REFS = lambda Aval: (("g_e", GE2), ("ν(g_e)g_e", float(nuu(GE2)) * GE2), ("A", Aval * GE2))
P("   m = 1 [%] je Bezug:  g_e (Newtonsches Zentrum) | ν(g_e)g_e (Feldform) | A (Massenmittel der eigenen Massen)")
for lab, gfun, Aval in (("Punktmasse (A = g_e)", lambda R: XK / R**2, 1.0),
                        ("Hernquist a = 3 kpc", lambda R: hern_g(XK, 3.0, R), Ah),
                        ("Exponentialscheibe 3 kpc", lambda R: float(expd_g(XK, RD, R)), Axe)):
    P(f"   {lab}  (A/g_e = {Aval:.3f})")
    for Rk in (10, 20, 30, 50):
        vals = [m1(gfun, Rk, ref)[0] for _, ref in REFS(Aval)]
        vbar = m1(gfun, Rk, Aval * GE2)[1] * math.sqrt(A0 * kpc) / 1e3
        P(f"      R = {Rk:2d} kpc, g_i/a₀ = {gfun(Rk):6.3f}, v̄_c = {vbar:6.1f} km/s:   " + " | ".join(f"{v:6.2f}" for v in vals))

# ============================================================================ 3
H("3   Druckgestützte Systeme: Jeans (isotrop, Plummer b = R_e) mit a_int im Bezug A; Eingaben SYS aus glaettung_stetig.py")
KMS2 = A0 * kpc / 1e6                                          # (a₀·kpc) -> (km/s)²
GMu = lambda M: G * M * M_sun / (A0 * kpc**2)
SYS = [  # Name, M [M☉], R_e [kpc], Wirt [M☉], Abstand [kpc], σ gemessen, Fehler  (wie glaettung_stetig.py)
    ("Pal 14 [G]", 2e4, 0.0247, 6e10, 71.0, 0.38, 0.12),
    ("NGC 2419 [G]", 9e5, 0.019, 6e10, 90.0, 4.1, 0.5),
    ("Crater II [G]", 3.2e5, 1.066, 6e10, 117.5, 2.7, 0.3),
    ("Draco [G]", 5.8e5, 0.221, 6e10, 76.0, 9.1, 1.2),
    ("Fornax [G]", 4.0e7, 0.71, 6e10, 147.0, 11.7, 0.9),
    ("DF2 [G]", 2.0e8, 2.2, 1e11, 80.0, 8.5, 2.7),
    ("DF4 [G]", 1.5e8, 1.6, 1e11, 165.0, 4.2, 3.3)]
MU, WMU = np.polynomial.legendre.leggauss(400)
def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return math.sqrt(max(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r), 0) * KMS2)
def inward(g, ge, Aref):
    """einwärts gerichtete Komponente −e_r·a_int je Richtung μ = e_r·ĝ_e; g (N,), Rückgabe (N, nμ)"""
    G_, M_ = g[:, None], MU[None, :]
    gm = np.sqrt(np.maximum(G_**2 + ge**2 - 2 * G_ * ge * M_, 0))
    # a = ν(gm)(−g e_r + g_e);  −e_r·a = ν(gm)(g − ge μ);  −e_r·(−A) = +A μ
    return nuu(gm) * (G_ - ge * M_) + Aref * M_
P("   σ_los im Fenster R < R_e [km/s]. Feldform: Monopol von ν(|g|)g (wie glaettung_stetig.py). Bezug A: Monopol von ν(|g|)g − A.")
P("   Zusammenhalt: Massenanteil, an dem a_int radial nach außen zeigt; in Klammern a_int(0)/g_i(r½) (Schub des Kerns gegen das Ganze, längs g_e).")
P(f"   {'System':14s} {'g_i/a₀':>7s} {'g_e/a₀':>7s} {'A/g_e':>6s} {'ν(g_e)':>6s} {'Newton':>7s} {'Feldform':>8s} {'Bezug A':>8s} "
  f"{'Virial A':>8s} {'gemessen':>11s} | außen (Kern): Bezug g_e | ν(g_e)g_e | A")
ROWS = {}
for nm, M, Re, Mh, D, so, er in SYS:
    b = Re; GM = GMu(M); ge = GMu(Mh) / D**2; X = GM / b**2
    rh = b * rh_b; gi = GM / 2 / rh**2
    r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, 5000)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5; Grho = 3 * GM / (4 * np.pi * b**3) * (1 + (r / b)**2)**-2.5; gN = GMenc / r**2
    Aval = plummer_A(X, ge) * ge
    mono = lambda Aref: 0.5 * np.sum(WMU * inward(gN, ge, Aref), axis=1)
    sN = jeans(r, Grho, gN, Re)
    sF = jeans(r, Grho, mono(0.0), Re)
    sA = jeans(r, Grho, mono(Aval), Re)
    # Virial (3D, Bezug A): 3σ²_3D-Mittel = ∫ρ r·(−e_r·a_int) dV / M;  σ_1D² = das Drittel davon (ganze Kugel, kein Fenster)
    vir = lambda Aref: math.sqrt(np.trapezoid(4 * np.pi * r**2 * Grho * r * mono(Aref), r) / GM / 3 * KMS2)
    sVA, sV0, sVN = vir(Aval), vir(0.0), math.sqrt(np.trapezoid(4 * np.pi * r**2 * Grho * r * gN, r) / GM / 3 * KMS2)
    rout = []
    wq, wwq = gl_split(np.linspace(0, 1, 41), 20); gq = X * wq * (1 - wq * wq); dmq = 3 * wq * wq * wwq
    for Aref in (ge, float(nuu(ge)) * ge, Aval):
        frac = np.sum(dmq[:, None] * 0.5 * WMU[None, :] * (inward(gq, ge, Aref) < 0))
        rout.append(f"{frac:5.3f} ({(float(nuu(ge)) * ge - Aref) / gi:+6.3f})")
    ROWS[nm] = dict(gi=gi, ge=ge, A=Aval / ge, nue=float(nuu(ge)), sN=sN, sF=sF, sA=sA, so=so, er=er)
    P(f"   {nm:14s} {gi:7.4f} {ge:7.4f} {Aval/ge:6.3f} {float(nuu(ge)):6.3f} {sN:7.2f} {sF:8.2f} {sA:8.2f} {sVA:8.2f} "
      f"{so:5.2f} ± {er:4.2f} | " + " | ".join(rout))
    ROWS[nm].update(sVA=sVA, sV0=sV0, sVN=sVN)
P("   Virial ohne Bezug (A = 0) gleich Virial mit A: " + ",  ".join(f"{n.split(' [')[0]} {d['sV0']:.3f}/{d['sVA']:.3f}" for n, d in ROWS.items()))
P("   Größte Differenz Feldform − Bezug A in σ: " + f"{max(abs(d['sF'] - d['sA']) for d in ROWS.values()):.1e} km/s")

# ============================================================================ 4
H("4   Ordnung: verlangte Verstärkung (σ_gem/σ_Bezug A)² gegen A/g_e, A/(ν(g_e)g_e), g_i/g_e")
MEAS = {n: (d["so"] / d["sA"])**2 for n, d in ROWS.items()}
for lab, fn in (("A/g_e", lambda d: d["A"]), ("A/ν(g_e)g_e", lambda d: d["A"] / d["nue"]), ("g_i/g_e", lambda d: d["gi"] / d["ge"])):
    srt = sorted(MEAS, key=lambda n: fn(ROWS[n]))
    P(f"   {lab:12s}: " + "  ".join(f"{n.split(' [')[0]} ({fn(ROWS[n]):.3g}: {MEAS[n]:.2f})" for n in srt))

# ============================================================================ 5
H("5   Unabhängige Prüfung")
P("   a) g_i ≫ g_e, erste Ordnung: A/g_e = ⟨ν(g_i)(1 + K/3)⟩_ρ, K = d ln ν/d ln g. Tiefbereich ν ≈ 1/s + 1/2 + s/12:")
P("      ν(1 + K/3) = ν + (s/6)ν' ≈ (5/6)/s + 1/2 + 7s/72.  Plummer: ⟨(g_i b²/GM)^{−1/2}⟩ = ∫₀¹ m^{−1/6}(1 − m^{2/3})^{−1/2} dm")
Bv = 1.5 * gamma(1.25) * gamma(0.5) / gamma(1.75); Bs = 1.5 * gamma(1.75) * gamma(1.5) / gamma(3.25)
P(f"      = (3/2)·B(5/4, 1/2) = {Bv:.5f}   (Handwert: Γ(5/4) = 0,90640, Γ(1/2) = 1,77245, Γ(7/4) = 0,91906 → 2,6220);"
  f"  ⟨(g_i b²/GM)^{{1/2}}⟩ = (3/2)B(7/4, 3/2) = {Bs:.5f}")
P("      Rest bei endlichem g_e: dort, wo g_i < g_e (Kern, Rand), fällt ν(g_i) weg; Ordnung √(g_e/X)·(5/3)/√X")
for X, ge in ((0.03, 1e-5), (0.03, 1e-7), (0.03, 1e-9), (0.01, 1e-9)):
    pred = (5 / 6) * Bv / math.sqrt(X) + 0.5 + (7 / 72) * math.sqrt(X) * Bs
    P(f"      X = {X}, g_e = {ge:g}: Näherung (5/6)·B/√X + 1/2 + (7/72)√X·B' = {pred:.4f};  Quadratur {plummer_A(X, ge):.4f}")
P("   b) g_i ≪ g_e, zweite Ordnung: A = ν(g_e)g_e·[1 + ⟨g_i²⟩(g_e h'' + 4h')/(6 ν g_e)], h = ν(g); Plummer ⟨g_i²⟩ = X²·(3/2)B(5/2, 3)")
g2 = 1.5 * 2 / (2.5 * 3.5 * 4.5)
P(f"      (3/2)B(5/2, 3) = (3/2)·2/(5/2·7/2·9/2) = {g2:.5f}")
for X, ge in ((0.01, 0.3), (0.003, 0.1), (0.03, 1.0)):
    hh = 1e-4 * ge; h = lambda x: float(nuu(x))
    h1 = (h(ge + hh) - h(ge - hh)) / (2 * hh); h2 = (h(ge + hh) - 2 * h(ge) + h(ge - hh)) / hh**2
    pred = h(ge) * (1 + g2 * X**2 * (ge * h2 + 4 * h1) / (6 * h(ge) * ge))
    P(f"      X = {X}, g_e = {ge}: ν(g_e) = {h(ge):.5f}, Näherung {pred:.5f}, Quadratur {plummer_A(X, ge):.5f}")
P("   c) Monopol eines konstanten Vektors: ½∫₋₁¹ μ dμ = 0, daher gibt jeder Bezug A dieselbe isotrope Jeans-Lösung wie die Feldform;")
P("      Virial: ∫ρ x·A dV = A·∫ρ x dV = 0 (Schwerpunkt im Ursprung). Zahlen in 3 stimmen damit überein.")
P("   d) Schalenformel (Substitution |g_tot| = u²) gegen direkte μ-Quadratur (4000 Gauß-Knoten):")
MU4, WMU4 = np.polynomial.legendre.leggauss(4000)
for g, ge in ((0.05, 0.03), (0.0301, 0.03), (2.0, 0.01), (0.001, 0.3)):
    gm = np.sqrt(np.maximum(g * g + ge * ge - 2 * g * ge * MU4, 0))
    P(f"      g = {g}, g_e = {ge}: {float(shell(g, ge)[0]):.8f}  /  {0.5 * np.sum(WMU4 * nuu(gm) * (ge - g * MU4)):.8f}")
P("   e) Punktmasse, Hand: A = g_e, weil ν(|g_i + g_e|) → 1 für die ganze Masse (g_i → ∞ im Zentrum); Plummer b → 0 in 2 zeigt es.")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "galaxie_als_ganzes.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
