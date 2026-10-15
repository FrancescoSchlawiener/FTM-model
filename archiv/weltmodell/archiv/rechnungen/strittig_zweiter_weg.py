# -*- coding: utf-8 -*-
"""
Strittige Punkte, zweiter Rechenweg (unabhängig von hopf_knoten, info_geometrie, modulraum_zahlen, kontakt_atmen).
Rechnung, noch nicht im Modell.

1 Hopf-Zahl: axiale Eichung A_z = 0 (statt Coulomb/FFT), zentrale Differenzen 2. Ordnung; a_J = 0,5/1/2/5,
  Vorzeichen von a_J, Windung p (Ladung); dazu Quadratur auf S³ in Torus-Koordinaten (a_J tritt dort nicht auf).
2 Verlinkte Tori, irrationale Wicklung: geschlossene Kettenbruch-Näherungen p_n/q_n (Gauß-Integral numerisch),
  Helizität aus dem Fluss (Arnold/Moffatt), Milnor-Invarianten.
3 χ als Geodäte: zweiparametrige Familien (p ∝ e^{−θ₁n−θ₂n²}; negativ-binomial), Christoffel-Symbole numerisch,
  Geodätengleichung integriert.
4 Clausius-Bilanz dQ = T dS auf der Halbgeraden: Laplace-Integral der Rückkehrdichte exakt, Faktor der Periode.
5 Q aus Phasenraum: Wirkungsraum-Volumen, Kapazität (Fläche), Fenster aus der Zelle, fensterfreie Kandidaten.
6 Stabilisatoren in SL(2,ℤ) und PSL(2,ℤ) durch Aufzählen.
7 Kontakt: Stokes-Beweis im Fourierraum; Kontaktformen mit u transversal bzw. S¹ = Zeit: freie Funktionen.
8 Resonator: Rückkehrdichte aus N Maxwell-/Lorentz-Moden gegen Lévy-Dichte T^{−3/2}; lokaler Exponent.
Konstanten aus ftm_kern.py. Ausgabe: strittig_zweiter_weg.txt
"""
import math
import itertools
import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation
from scipy.interpolate import CubicSpline
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, nnls
from ftm_kern import G, c, hbar, kB, A0, M_sun, kpc, Q_KOSM, q_KOSM, zelle_aus_a0, schliesslaenge, tau_D

out = []


def P(s=""):
    out.append(s)
    print(s, flush=True)


def kopf(t):
    P("=" * 100)
    P(t)
    P("=" * 100)


# ================================================================== 1
kopf("1  Hopf-Zahl, zweiter Weg: axiale Eichung A_z = 0, Differenzen 2. Ordnung; Quadratur auf S³")


def nfeld(p, q, a, L, N):
    x = np.linspace(-L, L, N)
    h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = X**2 + Y**2 + Z**2
    Z0 = np.conj(X**2 + Y**2 + (Z - 1j * a)**2) / (r2 + a * a)     # inverse Stereographie, Radius a
    Z1 = 2 * a * (X + 1j * Y) / (r2 + a * a)
    del X, Y, Z, r2
    num = Z1**p if p >= 0 else np.conj(Z1)**(-p)
    den = Z0**q
    s = np.abs(num)**2 + np.abs(den)**2
    cc = num * np.conj(den)
    return np.array([2 * cc.real / s, 2 * cc.imag / s, (np.abs(den)**2 - np.abs(num)**2) / s]), h


def hopf_axial(p, q, a, L, N):
    n, h = nfeld(p, q, a, L, N)
    d = [np.gradient(n, h, axis=i + 1) for i in range(3)]
    F = lambda i, j: np.einsum('a...,a...->...', n, np.cross(d[i], d[j], axis=0))
    B = np.array([F(1, 2), F(2, 0), F(0, 1)])
    del d, n
    # axiale Eichung: A_x = ∫_{-L}^z B_y dz', A_y = −∫_{-L}^z B_x dz', A_z = 0  (Trapez kumulativ)
    def cum(f):
        g = np.cumsum(f, axis=2) * h
        return g - 0.5 * h * (f + f[:, :, :1])
    Ax, Ay = cum(B[1]), -cum(B[0])
    H = -np.sum(Ax * B[0] + Ay * B[1]) * h**3 / (16 * np.pi**2)   # Vorzeichen wie hopf_knoten: w = Z1/Z0, a_J > 0 ⇒ +1
    # Probe: rot A = B (z-Komponente)
    Bz_rek = np.gradient(Ay, h, axis=0) - np.gradient(Ax, h, axis=1)
    fehler = np.max(np.abs(Bz_rek - B[2])[4:-4, 4:-4, 4:-4]) / np.max(np.abs(B[2]))
    return H, fehler


N1 = 150
P(f"1a feste Box ±8, N = {N1} (Gitter unabhängig von a_J), (p,q) = (1,1):")
H_ref = None
for a in (0.5, 1.0, 2.0, 5.0):
    H, f = hopf_axial(1, 1, a, 8.0, N1)
    P(f"   a_J = {a:4.1f}: H = {H:+.4f}   (h/a = {16/(N1-1)/a:.3f}, Box/a = {8/a:5.1f}; Probe rot A = B: {f:.1e})")
P(f"1b Box mitskaliert ±8a_J, N = {N1}: a_J = 0,5: H = {hopf_axial(1, 1, 0.5, 4.0, N1)[0]:+.4f}, a_J = 5: H = {hopf_axial(1, 1, 5.0, 40.0, N1)[0]:+.4f}"
  "  ⇒ Abweichungen in 1a nur Gitter (a_J = 0,5) bzw. Box (a_J = 5)")
P("1c Vorzeichen und Windung (Box ±8, a_J = ±1):")
for p, q, a in ((1, 1, 1.0), (1, 1, -1.0), (-1, 1, 1.0), (2, 1, 1.0), (2, 1, -1.0), (1, 2, 1.0)):
    H, _ = hopf_axial(p, q, a, 8.0, N1)
    P(f"   (p,q) = ({p:+d},{q}), a_J = {a:+.0f}: H = {H:+.4f}   p·q·sign(a_J) = {p*q*int(np.sign(a)):+d}")

# Quadratur auf S³: Z0 = cos η e^{iξ1}, Z1 = sin η e^{iξ2}; α = [(f−1)p dξ2 − f q dξ1]/2π, f = 1/(1+|w|²)
P("1d S³-Quadratur (Torus-Koordinaten η, ξ1, ξ2; Eichung α = [(f−1)p dξ2 − f q dξ1]/2π, glatt an beiden Kernkreisen):")
for p, q in ((1, 1), (2, 1), (1, 2), (3, 2)):
    f = lambda e: 1 / (1 + (np.sin(e)**p / np.cos(e)**q)**2)
    fp = lambda e: (f(e + 1e-6) - f(e - 1e-6)) / 2e-6
    # α∧dα = −(pq f'/4π²) dξ1 dη dξ2  →  H = −pq ∫ f' dη
    Hq = -p * q * quad(fp, 1e-9, np.pi / 2 - 1e-9)[0]
    P(f"   (p,q) = ({p},{q}): H = {Hq:+.6f}   (a_J kommt in dieser Rechnung nicht vor; er ist nur der Radius der Stereographie)")
P("   ⇒ H = p·q·sign(a_J): Betrag von a_J geht nicht ein; a_J = 0 ⇒ Abbildung konstant, H = 0; Umkehr des Drehsinns (a_J → −a_J)")
P("     kehrt H um. Mit Q_el nur über die Windung p (Ladungswindung), linear.")

# ================================================================== 2
kopf("2  Verlinkte Tori bei irrationaler Wicklung: Kettenbruch-Näherungen, Helizität, Milnor")


def torus_kurve(rt, m, k, npts, R0=1.0, ph=0.0):
    t = np.linspace(0, 2 * np.pi, npts, endpoint=False)
    phi, th = m * t, k * t + ph                  # m Umläufe um die z-Achse, k Umläufe um den Ring
    x = np.stack([(R0 + rt * np.cos(th)) * np.cos(phi), (R0 + rt * np.cos(th)) * np.sin(phi), rt * np.sin(th)], 1)
    dx = np.roll(x, -1, 0) - x
    return x + 0.5 * dx, dx


def gauss_lk(c1, c2):
    x1, d1 = c1
    x2, d2 = c2
    s = 0.0
    for i in range(0, len(x1), 400):
        r = x1[i:i + 400, None, :] - x2[None, :, :]
        cr = np.cross(d1[i:i + 400, None, :], d2[None, :, :])
        s += np.sum(np.einsum('ijk,ijk->ij', r, cr) / np.linalg.norm(r, axis=2)**3)
    return s / (4 * np.pi)


W = (math.sqrt(5) - 1) / 2
fib = [(1, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 13), (13, 21)]
P(f"W = (√5−1)/2 = {W:.6f} (Ring-Umläufe je Achsenumlauf) auf zwei Tori (Rohrradius 0,3 und 0,6 um den Ring R₀ = 1).")
P("Geschlossene Näherung p_n/q_n: q_n Achsenumläufe, p_n Ringumläufe; Gauß-Integral numerisch:")
P(f"   {'p/q':>7s} {'Lk (Gauß)':>12s} {'p·q':>6s} {'Lk/q²':>9s} {'|p/q − W|':>10s}")
for pn, qn in fib:
    npts = 260 * qn
    Lk = gauss_lk(torus_kurve(0.3, qn, pn, npts), torus_kurve(0.6, qn, pn, npts, ph=0.7))
    P(f"   {pn:>3d}/{qn:<3d} {Lk:12.4f} {pn*qn:6d} {abs(Lk)/qn**2:9.5f} {abs(pn/qn-W):10.2e}")
P("   Die ganzen Zahlen p_n·q_n hängen an der gewählten Näherung (Abbruch des Kettenbruchs), nicht an W allein.")
# Helizität aus Flüssen (Arnold/Moffatt): Feld tangential an verschachtelten Tori, Wicklung W(ψ)
# H = 2∫ W(Φ) Φ dΦ (Φ = Ringfluss innerhalb der Fläche), asymptotische Verlinkung je Flusspaar = W
for Wf, name in ((lambda x: W + 0 * x, "W konstant"), (lambda x: W * (1 + 0.3 * x), "W(Φ) = W(1+0,3Φ)")):
    Hh = 2 * quad(lambda x: Wf(x) * x, 0, 1)[0]
    P(f"   Helizität je Ringfluss², {name}: H/Φ² = 2∫W(Φ)Φ dΦ = {Hh:.6f}")
P("   Milnor: μ̄(12) = Lk (nur für geschlossene Kurven); μ̄(123) ist nur modulo ggT(Lk12, Lk13, Lk23) bestimmt,")
P("   bei Kurven auf verschachtelten Tori sind alle Lk_ij ≠ 0 ⇒ keine neue ganze Zahl. Irrationales W: keine Kurve")
P("   schließt, μ̄ undefiniert; es bleibt die asymptotische Verlinkung (Arnold) lim Lk(T₁,T₂)/(T₁T₂) = W, reell.")
P("   Entkopplung: Lk ist unter dem Fluss erhalten (Isotopie); Lk(T) = W·T² wächst, fällt nie auf 0 ⇒ keine")
P("   Umlaufzahl bis zur Entkopplung. Die einzige bestimmte Zahl ist W selbst (und der Kettenbruch [0;1,1,1,…]).")

# ================================================================== 3
kopf("3  χ als Geodäte: zweiparametrige Familien, Geodätengleichung explizit")
nn = np.arange(0, 4000, dtype=float)


T2 = np.log1p(nn)


def fam_quad(th):            # p(n) ∝ exp(−θ1 n − θ2 ln(1+n)): Maximum Entropie bei festem <n>, <ln(1+n)>
    e = -th[0] * nn - th[1] * T2
    e -= e.max()
    p = np.exp(e)
    return p / p.sum()


def fam_nb(th):              # negativ-binomial: s = θ1 (q = e^{−s}), r = θ2 (r = 1: geometrisch)
    s, r = th
    from scipy.special import gammaln
    lp = gammaln(nn + r) - gammaln(r) - gammaln(nn + 1) + r * math.log(-math.expm1(-s)) - s * nn
    p = np.exp(lp)
    return p / p.sum()


def fisher(fam, th, h=1e-4):
    p = fam(th)
    if fam is fam_quad:      # Exponentialfamilie: g = Cov(T) exakt
        Tm = np.array([nn, T2])
        m = Tm @ p
        D = Tm - m[:, None]
        return (D * p) @ D.T
    lp = np.log(np.maximum(p, 1e-300))
    sc = []
    for i in range(2):
        d = np.zeros(2); d[i] = h
        sc.append((np.log(np.maximum(fam(th + d), 1e-300)) - np.log(np.maximum(fam(th - d), 1e-300))) / (2 * h))
    return np.array([[np.sum(p * sc[i] * sc[j]) for j in range(2)] for i in range(2)])


def christoffel(fam, th, h=1e-3):
    dg = []
    for k in range(2):
        d = np.zeros(2); d[k] = h
        dg.append((fisher(fam, th + d) - fisher(fam, th - d)) / (2 * h))
    g = fisher(fam, th)
    gi = np.linalg.inv(g)
    Gam1 = np.zeros((2, 2, 2))   # Γ_{ij,k}
    for i, j, k in itertools.product(range(2), repeat=3):
        Gam1[i, j, k] = 0.5 * (dg[i][j, k] + dg[j][i, k] - dg[k][i, j])
    return np.einsum('kl,ijl->kij', gi, Gam1), g


P("Familie A: p(n) ∝ e^{−θ₁n−θ₂ln(1+n)} (Maximum Entropie bei festem <n>, <ln(1+n)>; θ₂ beider Vorzeichen erlaubt).")
P("   Rückkehrverteilung: θ₂ = 0, θ₁ = s. (Mit n² statt ln(1+n) wäre θ₂ = 0 der Rand der Familie, θ₂ < 0 nicht normierbar.)")
P("Kurve θ₂ = 0, Parameter s: Levi-Civita-Beschleunigung Γ^k_11; Geodäte (bis auf Umparametrisierung) ⇔ Γ²_11 = 0.")
P("Für e-Verbindung (α = 1) Γ = 0 in θ exakt (Exponentialfamilie): Gerade θ₂ = 0 ist e-Geodäte mit affinem Parameter s.")
P(f"   {'s':>5s} {'Γ¹_11':>11s} {'Γ²_11 (LC)':>12s} {'Normalkrümmung k_g':>19s}  m-Verb. (α=−1): 2·Γ²_11")
for s in (0.3, 1.0, 2.0):
    Gm, g = christoffel(fam_quad, np.array([s, 0.0]))
    # geodätische Krümmung: |Normalanteil der Beschleunigung| / |Tangente|²
    acc = Gm[:, 0, 0]
    t = np.array([1.0, 0.0])
    tn = t @ g @ t
    an = acc - (acc @ g @ t) / tn * t
    kg = math.sqrt(an @ g @ an) / tn
    P(f"   {s:5.1f} {Gm[0,0,0]:11.4e} {Gm[1,0,0]:12.4e} {kg:19.4e}  {2*Gm[1,0,0]:.4e}")
P("Geodätengleichung (Levi-Civita) integriert ab (θ₁,θ₂) = (1, 0), Richtung ∂_θ₁, Einheitsgeschwindigkeit:")


def geod_rhs(fam):
    def f(t, y):
        th, v = y[:2], y[2:]
        Gm, _ = christoffel(fam, th)
        return np.concatenate([v, -np.einsum('kij,i,j->k', Gm, v, v)])
    return f


g0 = fisher(fam_quad, np.array([1.0, 0.0]))
v0 = np.array([1 / math.sqrt(g0[0, 0]), 0.0])
sol = solve_ivp(geod_rhs(fam_quad), (0, 0.8), np.concatenate([[1.0, 0.0], v0]), t_eval=[0.2, 0.4, 0.8], rtol=1e-7)
for tt, a1, a2 in zip(sol.t, sol.y[0], sol.y[1]):
    P(f"   Bogenlänge {tt:.1f}: θ₁ = {a1:.5f}, θ₂ = {a2:+.3e}   (θ₂ ≠ 0: die LC-Geodäte verlässt die Rückkehrfamilie)")
P("Familie B: negativ-binomial (s, r), r = 1 geometrisch (r unabhängige Rückkehrkanäle):")
for s in (0.5, 1.0, 2.0):
    Gm, g = christoffel(fam_nb, np.array([s, 1.0]))
    P(f"   s = {s:.1f}: Γ²_11 (LC) = {Gm[1,0,0]:+.4e}   (≠ 0 ⇒ r = 1 keine LC-Geodäte)")
P("Parameter entlang der e-Geodäte: affin ist s (und jede Funktion a·s + b). χ(s) = s/(e^s−1) ist nicht affin in s")
for s in (0.5, 1.0, 2.0):
    P(f"   s = {s:.1f}: χ = {s/math.expm1(s):.5f}, dχ/ds = {(math.expm1(s)-s*math.exp(s))/math.expm1(s)**2:+.5f}, d²χ/ds² ≠ 0")
P("   ⇒ Formulierbar ist: die Rückkehrfamilie ist e-Geodäte (e-flach) in jeder Exponentialfamilie, die sie enthält, und steht")
P("     senkrecht (Fisher) auf den Blättern festen Erwartungswerts <n> (Pythagoras). Nicht wahr: χ als Geodäte oder affiner")
P("     Parameter; Levi-Civita-Geodäte ist sie in keiner der beiden Familien.")

# ================================================================== 4
kopf("4  dQ = T dS auf der Halbgeraden: Rückkehrpfade, Verlust, Periode")
P("Rückkehrdichte der Halbgeraden (Diffusion, erste Rückkehr bei Zeit t): f(t) = √(τ_D/4π) t^{−3/2} e^{−τ_D/4t}.")
P("Gewicht eines Pfades: Überleben e^{−Γt} (Wärme = Verlust); Zahl der Pfade e^{S(t)}, S = ln f.")
P("Clausius am Sattel: Γ = dS/dt ⇒ Γ = τ_D/4t² − 3/2t; ohne Vorfaktor: t* = ½√(τ_D/Γ), Exponent Γt*+τ_D/4t* = √(Γτ_D) = s.")
P(f"   {'Γτ_D':>8s} {'∫f e^{−Γt}dt':>14s} {'e^{−√(Γτ_D)}':>14s} {'Γt*':>8s} {'τ_D/4t*':>8s}")
tauD = 1.0
for gt in (0.01, 0.1, 1.0, 10.0):
    Gam = gt / tauD
    f = lambda t: math.sqrt(tauD / (4 * math.pi)) * t**-1.5 * math.exp(-tauD / (4 * t) - Gam * t)
    tm = 0.5 * math.sqrt(tauD / Gam)
    I = quad(f, 0, tm, limit=400)[0] + quad(f, tm, np.inf, limit=400)[0]
    S = lambda t: -tauD / (4 * t) - 1.5 * math.log(t) + 0.5 * math.log(tauD / (4 * math.pi))
    ts = brentq(lambda t: tauD / (4 * t * t) - 1.5 / t - Gam, 1e-6, 1e6)
    S2 = -tauD / (2 * ts**3) + 1.5 / ts**2
    gauss = math.exp(S(ts) - Gam * ts) * math.sqrt(2 * math.pi / -S2)
    t0 = 0.5 * math.sqrt(tauD / Gam)
    P(f"   {gt:8.2f} {I:14.8f} {math.exp(-math.sqrt(gt)):14.8f} {Gam*t0:8.4f} {tauD/(4*t0):8.4f}   (Clausius: Wärme Γt* = Entropie τ_D/4t*)")
P("   ⇒ Das Laplace-Integral gibt e^{−s} exakt (der Vorfaktor t^{−3/2} hebt sich genau); keine freie Konstante im Exponenten.")
P("Faktor in v⁴: s² = Γ·τ_D mit Γ = a_N/c. Welche Periode τ_D, welcher Faktor in v⁴ = k·GMa₀ (tiefer Bereich a → √(a_N c/τ_D)):")
Lam = schliesslaenge(A0)
for name, tD in (("Umlauf 2Λ/c (Kac, Λ = c²/2a₀)", 2 * Lam / c), ("einfach Λ/c", Lam / c),
                 ("ħ/kT mit T = ħa₀/(2πck) (2π, Unruh-Form)", 2 * math.pi * c / A0), ("Diffusion Λ²/D, D = cΛ", Lam / c)):
    k = (c / tD) / A0
    P(f"   {name:42s}: τ_D = {tD:.4e} s, k = c/(τ_D a₀) = {k:.4f}")
M = 1e11 * M_sun
for rk in (100, 1000):
    r = rk * kpc
    aN = G * M / r**2
    s = math.sqrt(aN * tau_D(A0) / c)
    a = aN / -math.expm1(-s)
    P(f"   Punktmasse 10¹¹ M☉, r = {rk} kpc, τ_D = 2Λ/c: v⁴/(GMa₀) = {(a*r)**2/(G*M*A0):.4f}  (→ 1 für r → ∞)")
P("   ⇒ Bilanz dQ = T dS ist dieselbe Rechnung wie s = min_t[Γt + τ_D/4t] (07): Faktor 1 genau dann, wenn die Periode der")
P("     Umlauf 2Λ/c über die Schließlänge ist (T_I = ħa₀/ck, ohne 2π). Die Clausius-Relation setzt den Faktor nicht; ihn setzt Λ.")

# ================================================================== 5
kopf("5  Q aus dem Phasenraum mit anderen Normierungen (GM = 1, r* = 1)")


def aR(r):
    return (1 / r**2) / -math.expm1(-1 / r)


_lr = np.linspace(math.log(1e-4), math.log(2e3), 3000)
_dl = lambda r: (1 / r**2) * math.exp(-1 / r) / -math.expm1(-1 / r)          # a_R − a_N
_cum = [0.0]
for i in range(1, len(_lr)):
    _cum.append(_cum[-1] + quad(lambda u: _dl(math.exp(u)) * math.exp(u), _lr[i - 1], _lr[i])[0])
_sp = CubicSpline(_lr, np.array(_cum) - np.interp(0.0, _lr, _cum))


def PhiR(r):
    return 1.0 - 1.0 / r + float(_sp(max(math.log(r), _lr[0])))


def PhiN(r):
    return 1.0 - 1.0 / r


def Jr(PhiF, E, L, R=None):
    g = lambda r: 2 * (E - PhiF(r)) - L * L / (r * r)
    if L == 0:
        r1, r2 = 1e-12, R
    else:
        rp = brentq(lambda r: -L * L / r**3 + (aR(r) if PhiF is PhiR else 1 / r**2), 1e-12, R)
        if g(rp) <= 0:
            return 0.0
        r1 = brentq(g, 1e-13, rp)
        r2 = brentq(g, rp, R)
    return quad(lambda r: math.sqrt(max(g(r), 0)), r1, r2, limit=200)[0] / math.pi


def vol_wirkung(PhiF, R, nL=60):
    """Volumen {E < Φ(R)} im Wirkungsraum: (2π)³ ∫ 2L·J_r(E,L) dL (L_z ∈ [−L, L])."""
    E = PhiF(R)
    # Lmax: Kreisbahn mit Energie E
    def Ecirc(r):
        v2 = r * (aR(r) if PhiF is PhiR else 1 / r**2)
        return PhiF(r) + v2 / 2
    rc = brentq(lambda r: Ecirc(r) - E, 1e-4, R)
    Lmax = rc * math.sqrt(rc * (aR(rc) if PhiF is PhiR else 1 / rc**2))
    xs, ws = np.polynomial.legendre.leggauss(nL)
    Ls = 0.5 * Lmax * (xs + 1)
    tot = sum(w * 2 * L * Jr(PhiF, E, L, R=R) for L, w in zip(Ls, ws)) * 0.5 * Lmax
    return (2 * math.pi)**3 * tot, 2 * math.pi * Jr(PhiF, E, 0.0, R=R), 2 * math.pi * Lmax


P("(a) Wirkungs-Winkel-Volumen (2π)³∫dJ_r dL dL_z, (b) Kapazität als Fläche: kleinste Wirkung ∮p dq geschlossener")
P("    Charakteristiken des Randes (radial 2πJ_r(E,0) oder kreisförmig 2πL_c), (c) Volumen^{1/3} (Fläche je Freiheitsgrad).")
P(f"   {'R/r*':>6s} {'Q_Vol(Wirkung)':>15s} {'Q_Kap':>9s} {'Q_Vol^{1/3}':>12s}")
res = {}
for R in (1.0, 3.0, 5.0, 10.0, 30.0):
    vR, jrR, lcR = vol_wirkung(PhiR, R)
    vN, jrN, lcN = vol_wirkung(PhiN, R)
    capR, capN = min(jrR, lcR), min(jrN, lcN)
    res[R] = (vR / vN - 1, capR / capN - 1, (vR / vN)**(1 / 3) - 1)
    P(f"   {R:6.1f} {vR/vN-1:15.4f} {capR/capN-1:9.4f} {(vR/vN)**(1/3)-1:12.4f}")
P("   Probe (a) gegen Ortsraum-Volumen (modulraum_zahlen, R = 3: 2,4223): gleiche Zahl ⇒ Liouville, Normierung ändert nichts.")
P("(d) Fenster aus der Zelle R = L/2 bzw. R = Λ = c²/2a₀ (R/r* hängt an M):")
L_z = zelle_aus_a0(A0)
for Ms in (1e9, 1e11, 1e13):
    rs = math.sqrt(G * Ms * M_sun / A0)
    P(f"   M = {Ms:.0e} M☉: r* = {rs/kpc:7.2f} kpc, (L/2)/r* = {L_z/2/rs:.3e}, Λ/r* = {Lam/rs:.3e}  ⇒ Q ≈ (R/r*)^{{3/2}}·(ln)^{{3/2}} ≫ 5,36")
P("(e) Fensterfreie Kandidaten: R = r* (a_N = a₀): Q_Vol = {:.4f}, Q_Kap = {:.4f}; R → ∞: Q → ∞ (jede Bahn gebunden).".format(*res[1.0][:2]))
for idx, name in ((0, "Q_Vol"), (1, "Q_Kap"), (2, "Q_Vol^{1/3}")):
    try:
        f = lambda R: (lambda t: (t[0] / t[1] - 1) if idx == 0 else ((min(t[2], t[3]) / min(t[4], t[5]) - 1) if idx == 1 else (t[0] / t[1])**(1 / 3) - 1))(
            (lambda a, b: (a[0], b[0], a[1], a[2], b[1], b[2]))(vol_wirkung(PhiR, R, 30), vol_wirkung(PhiN, R, 30))) - Q_KOSM
        Rq = brentq(f, 2.0, 400.0, xtol=1e-3)
        P(f"   {name} = 5,364 bei R/r* = {Rq:.3f} (je Normierung ein anderes Fenster)")
    except ValueError:
        P(f"   {name}: 5,364 im Bereich R/r* ∈ [2, 400] nicht erreicht")
P("   ⇒ M und a₀ kürzen sich; es bleibt in jeder Normierung eine Funktion von R/r*, jede mit anderem R für 5,36.")

# ================================================================== 6
kopf("6  Stabilisatoren von i und ρ in SL(2,ℤ) und PSL(2,ℤ) (Aufzählen, Einträge |·| ≤ 4)")
rho = complex(-0.5, math.sqrt(3) / 2)
for z0, name in ((1j, "i"), (rho, "ρ")):
    st = []
    for a_, b_, c_, d_ in itertools.product(range(-4, 5), repeat=4):
        if a_ * d_ - b_ * c_ == 1 and abs((a_ * z0 + b_) / (c_ * z0 + d_) - z0) < 1e-12:
            st.append(np.array([[a_, b_], [c_, d_]]))
    ords = []
    for Mx in st:
        k, X = 1, Mx.copy()
        while not np.array_equal(X, np.eye(2, dtype=int)):
            X = X @ Mx; k += 1
        ords.append(k)
    psl = {tuple((Mx if (Mx[0, 0], Mx[0, 1]) > (0, 0) or (Mx[0, 0] == 0 and Mx[0, 1] > 0) else -Mx).ravel()) for Mx in st}
    P(f"   Stab({name}) in SL(2,ℤ): {len(st)} Elemente, Ordnungen {sorted(ords)} ⇒ zyklisch ℤ{len(st)} (enthält −I);"
      f" in PSL = SL/{{±I}}: {len(psl)} Elemente ⇒ ℤ{len(psl)}")
P("   Von Hand: S = [[0,−1],[1,0]], S² = −I, S⁴ = I (ℤ₄); ST = [[0,−1],[1,1]], (ST)³ = −I, (ST)⁶ = I (ℤ₆).")
P("   ⇒ Beide richtig: ℤ₄, ℤ₆ in SL(2,ℤ); ℤ₂, ℤ₃ in PSL(2,ℤ), das auf der oberen Halbebene treu wirkt (−I wirkt dort als 1).")
P("     Welche Gruppe gilt, entscheidet, ob −I im Modell wirkt: auf τ nicht; auf dem Gitter (ω₁, ω₂) als (−ω₁, −ω₂), z → −z.")

# ================================================================== 7
kopf("7  Kontaktgeometrie: Stokes unabhängig (Fourier), andere Wahlen")
P("Stokes im Fourierraum: β = Σ_k b_k e^{ik·x} auf T⁴; ∫_{T⁴}(dβ)² hat nur Paare (k, −k):")
P("   ∫(dβ)² ∝ Σ_k ε^{ijlm} k_i k_l b_{k,j} b_{−k,m} = 0, da k_i k_l symmetrisch, ε antisymmetrisch in (i,l).")
rng = np.random.default_rng(7)
Nf = 10
xg = np.arange(Nf) * 2 * np.pi / Nf
Xg = np.meshgrid(xg, xg, xg, xg, indexing='ij')
for trial in range(3):
    beta = [np.zeros((Nf,) * 4) for _ in range(4)]
    for _ in range(6):
        kv = rng.integers(-2, 3, 4)
        ph = sum(kv[i] * Xg[i] for i in range(4))
        for j in range(4):
            beta[j] += rng.normal() * np.cos(ph + rng.uniform(0, 6.3))
    kk = np.fft.fftfreq(Nf, 1 / Nf)
    Kg = np.meshgrid(kk, kk, kk, kk, indexing='ij')
    der = lambda f, i: np.real(np.fft.ifftn(1j * Kg[i] * np.fft.fftn(f)))
    F = [[der(beta[j], i) - der(beta[i], j) for j in range(4)] for i in range(4)]
    pf = F[0][1] * F[2][3] - F[0][2] * F[1][3] + F[0][3] * F[1][2]
    P(f"   Zufalls-β {trial+1} (spektrale Ableitung): Mittel Pf(dβ) = {pf.mean():+.2e}, Bereich [{pf.min():.2f}, {pf.max():.2f}] ⇒ Nullstelle")
P("   Beweis bestätigt: Reeb = u = ∂_t verlangt α = dt + β, β global auf T³×S¹, (dβ)² exakt ⇒ ∫ = 0 ⇒ entartet.")

x, y, z, w, t, vp = sp.symbols('x y z w t varpi', real=True)


def kontakt(alpha, X):
    n = len(X)
    A = [alpha[i] for i in range(n)]
    Fm = [[sp.diff(A[j], X[i]) - sp.diff(A[i], X[j]) for j in range(n)] for i in range(n)]
    tot = 0
    for perm in itertools.permutations(range(n)):
        sg = Permutation(list(perm)).signature()
        a_, b_, c_, d_, e_ = perm
        tot += sg * A[a_] * Fm[b_][c_] * Fm[d_][e_]
    return sp.simplify(tot / 8)


def reeb(alpha, X):
    n = len(X)
    Fm = sp.Matrix(n, n, lambda i, j: sp.diff(alpha[j], X[i]) - sp.diff(alpha[i], X[j]))
    R = sp.Matrix(sp.symbols('R0:%d' % n))
    eqs = list(Fm.T * R) + [sum(alpha[i] * R[i] for i in range(n)) - 1]
    sol = sp.solve(eqs, list(R), dict=True)
    return [sp.simplify(sol[0][Ri]) for Ri in R] if sol else None


Af, Bf, Cf = sp.Function('A')(t), sp.Function('B')(t), sp.Function('C')(t)
X5 = [x, y, z, w, t]
al = [Af * sp.cos(z), Af * sp.sin(z), 0, Bf, 1]          # α = dt + A(t)(cos z dx + sin z dy) + B(t) dw ; w ∈ S¹ innere Richtung
P("(i) T³×S¹_w×ℝ_t, u = ∂_t transversal (α(∂_t) = 1), α = dt + A(t)(cos z dx + sin z dy) + B(t)dw:")
P(f"   α∧(dα)² ∝ {sp.factor(kontakt(al, X5))}   ⇒ Kontakt für jedes A, B mit A·(AB' − BA') ≠ 0 (zwei freie Funktionen)")
Rb = reeb(al, X5)
P(f"   Reeb-Feld R = {Rb}  (R_t = 0: Reeb liegt immer in T³×S¹, nie längs u)")
alS = [sp.exp(vp) * sp.cos(z), sp.exp(vp) * sp.sin(z), 0, 0, 1]   # S¹ = Zeit t, ϖ ∈ ℝ offen
Ff = sp.Function('f')(vp)
alS = [Ff * sp.cos(z), Ff * sp.sin(z), 0, 0, 1]
X5b = [x, y, z, vp, t]
P("(ii) S¹ = Zeitrichtung t, innere Richtung ϖ offen: α = dt + f(ϖ)(cos z dx + sin z dy):")
P(f"   α∧(dα)² ∝ {sp.factor(kontakt(alS, X5b))}; Reeb = {reeb(alS, X5b)}  ⇒ Reeb = ∂_t für JEDES f mit f·f' ≠ 0")
P("   Darboux/Gray: jede Kontaktform ist lokal dt + Σ p_i dq_i, Kontaktstrukturen haben keine lokalen Invarianten ⇒")
P("   keine Funktion w(z) folgt aus der Kontaktbedingung allein; A(t), B(t) bzw. f(ϖ) bleiben frei (eine Funktion ≠ Parameter).")
for name, Aexpr in (("A = t^{2/3} (Materie)", t**sp.Rational(2, 3)), ("A = e^{t} (de Sitter)", sp.exp(t))):
    Bex = t if Aexpr != sp.exp(t) else sp.exp(2 * t)
    val = sp.factor(kontakt([Aexpr * sp.cos(z), Aexpr * sp.sin(z), 0, Bex, 1], X5))
    P(f"   Beispiel {name}, B = {Bex}: α∧(dα)² ∝ {val}  (beide Kontakt, verschiedenes Atmen)")

# ================================================================== 8
kopf("8  Analogon: Resonatoren (Maxwell-/Lorentz-Kern) gegen Lévy-Rückkehrdichte")
P("Rückkehranteil q(Γ) = f̂(Γ) = ∫f(t)e^{−Γt}dt; Regel verlangt 1 − f̂ ∝ √(Γτ) für kleines Γ (lokaler Exponent ½).")
tau = 1.0
lev = lambda G_: math.exp(-math.sqrt(G_ * tau))
P("   Ein Maxwell-Kern f = γe^{−γt}: f̂ = γ/(γ+Γ), 1 − f̂ = Γ/(γ+Γ) ⇒ a = a_N/(1−q) = a_N + γc (Versatz, keine Wurzel).")
P("   Ein Lorentz-Kern f ∝ e^{−γt}sin²(ω₀t/2): f̂ rational, 1 − f̂ = c₁Γ + … (Exponent 1).")
g0_, w0_ = 1.0, 3.0
fl = lambda t: math.exp(-g0_ * t) * math.sin(w0_ * t / 2)**2
nl = quad(fl, 0, 60, limit=400)[0]
Gs = np.logspace(-4, 2, 13)
lo = [1 - quad(lambda t: fl(t) * math.exp(-G_ * t), 0, 60, limit=400)[0] / nl for G_ in Gs]
# N Maxwell-Moden, Raten log-verteilt, Gewichte nicht negativ (Foster-Netz), angepasst an 1 − e^{−√Γτ} auf Γτ ∈ [1e-4, 1e2]
Gfit = np.logspace(-4, 2, 200)
P(f"   {'N Moden':>8s} {'max rel. Fehler 1−f̂':>20s} {'Exponent bei Γτ=1e-3':>21s} {'Exponent bei Γτ=1e-6':>21s}")
for Nm in (1, 2, 4, 8, 16):
    lam = np.logspace(-4, 2, Nm) if Nm > 1 else np.array([1.0])
    Amat = np.array([[Gv / (l + Gv) for l in lam] for Gv in Gfit])
    tgt = np.array([1 - lev(Gv) for Gv in Gfit])
    wts, _ = nnls(Amat / tgt[:, None], np.ones_like(tgt))
    fN = lambda Gv: sum(wk * Gv / (l + Gv) for wk, l in zip(wts, lam))
    err = np.max(np.abs(np.array([fN(Gv) for Gv in Gfit]) / tgt - 1))
    ex = lambda Gv: (math.log(fN(Gv * 1.01)) - math.log(fN(Gv / 1.01))) / (2 * math.log(1.01))
    P(f"   {Nm:8d} {err:20.3e} {ex(1e-3):21.3f} {ex(1e-6):21.3f}")
exl = lambda Gv: (math.log(1 - lev(Gv * 1.01)) - math.log(1 - lev(Gv / 1.01))) / (2 * math.log(1.01))
exo = (math.log(lo[1]) - math.log(lo[0])) / (math.log(Gs[1]) - math.log(Gs[0]))
P(f"   Lévy selbst: Exponent bei Γτ = 1e-6: {exl(1e-6):.3f};  ein Lorentz-Resonator (γ = 1, ω₀ = 3) bei Γ = 1e-4: {exo:.3f}")
P("   ⇒ Endlich viele Resonatoren (Maxwell oder Lorentz, parallel oder als endliche Kette): f̂ rational, Exponent 1 unterhalb der")
P("     kleinsten Rate; Wurzel nur zwischen kleinster und größter Rate. Exakt nur mit unendlich vielen Moden: Dichte der Raten")
P("     ∝ λ^{−1/2} (Foster) bzw. unendliche Mori-Kettenbruchkette (Cauer) = Leitung. Zwanzig-Mori-Kettenbruch IST eine Kette.")

kopf("Von Hand")
P("  1: (p,q) = (1,1): H = −pq∫₀^{π/2} f' dη = −(f(π/2) − f(0)) = −(0 − 1) = 1; a_J nicht darin.")
P("  2: 8/13: Lk = 13·8 = 104, 104/169 = 0,6154; W = 0,6180.")
P("  3: χ(1) = 1/(e−1) = 0,58198; affin wäre χ = a·s + b, aber χ(0) = 1, χ(1) = 0,582, χ(2) = 0,313: Differenzen 0,418 ≠ 0,269.")
P("  4: Γτ_D = 1: t* = ½, Γt* + τ_D/4t* = ½ + ½ = 1 = s; e^{−1} = 0,36788. 2Λ/c = c/a₀ = 2,998e8/1,1e-10 = 2,725e18 s.")
P("  6: S² = [[0,−1],[1,0]]² = [[−1,0],[0,−1]] = −I ⇒ Ordnung 4 in SL, 2 in PSL.")
P("  7: α = dt + A(cos z dx + sin z dy) + B dw: α∧(dα)² ∝ A(AB' − BA'); A = t^{2/3}, B = t: t^{2/3}(t^{2/3} − ⅔t^{2/3}) = t^{4/3}/3.")
P("  8: Maxwell γ = Γ: f̂ = ½; Lévy Γτ = 1: e^{−1} = 0,368.")

with open(__file__.replace('.py', '.txt'), 'w') as fh:
    fh.write("\n".join(out) + "\n")
