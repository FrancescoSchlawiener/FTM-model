# -*- coding: utf-8 -*-
"""
Gewicht der Schließungsregel als Rückkehr (Kac/Lévy) statt Maxwell: alle Folgen.
  w_M = x²/(1 + x²)          (jetzige Form, Kern gewicht_w)
  w_R = 1 − e^{−√x}          (Rückkehr mit der Ausdehnungsrate als Verlust), x = θ̃·τ_D, τ_D = c/a₀
1 Herleitung (Lévy-Dichte, Halbgerade bei Λ geschlossen: cosh, Reihe, Differenzen)
2 1 − w(z) im Hintergrund, θ̃ = 3H(z), a₀(z) aus dem Kern
3 lineares Wachstum: Übergabe (1 − w)ω/N* als nicht klumpender Anteil (Lesart 03), Kontrolle 5·10⁻⁴
4 Übergangsbreiten w 0,9 → 0,1 (Zykloide und Energieintegral), Einsatz bei hohem z
5 Umkehrkante (Energieintegral der Schale)
6 Folgen: Ω_s/(1 − w), frisch kollabierte Galaxien, N*, Linse um isolierte Galaxien, Filament, w = ½
7 Änderungsliste, 8 Handprüfungen
Nur Kern-Formeln (ftm_kern); Top-Hat und Galaxie hier auf eigenem Weg (Energieintegral statt ODE).
"""
import math, os, sys, io
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K
from ftm_kern import G, c, M_sun, kpc, Gyr, A0, Q_KOSM, OMEGA_M, OMEGA_R, E_hintergrund, a0_frueh_relativ

OUT = []
def P(s=""):
    OUT.append(s); print(s)
def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)

H0 = K.H0_si()
tauD = c / A0
thc = (A0 / c) / H0                      # θ_c in H0
X9R, X5R, X1R = math.log(10)**2, math.log(2)**2, math.log(1 / 0.9)**2
X9M, X5M, X1M = 3.0, 1.0, 1 / 3

def wM(x):
    x = max(x, 0.0); return x * x / (1 + x * x)
def wR(x):
    return -math.expm1(-math.sqrt(max(x, 0.0)))
def l10_1mw(x, form):                    # log10(1 − w), ohne Unterlauf
    return -math.sqrt(x) / math.log(10) if form == "R" else -math.log10(1 + x * x)

# ======================================================================================== 0
hdr("0  Grundgrößen")
P("   a₀ = %.2e, τ_D = c/a₀ = %.3f Gyr, θ_c = a₀/c = %.4e s⁻¹ = %.5f H₀, 3H₀τ_D = %.4f" % (A0, tauD / Gyr, A0 / c, thc, 3 / thc))
P("   Schwellen: w_R = 0,9 / 0,5 / 0,1 bei x = (ln 10)² = %.5f / (ln 2)² = %.5f / (ln(10/9))² = %.6f;  w_M bei 3 / 1 / 1/3" % (X9R, X5R, X1R))

# ======================================================================================== 1
hdr("1  Herleitung: Laplace-Transformierte der Dichte erster Rückkehrzeiten")
P("   Halbgerade [0, ∞), Start in der Tiefe ϖ, Diffusion D, Rückkehr bei 0: f(t) = ϖ/√(4πDt³)·e^{−ϖ²/(4Dt)}.")
P("   Von Hand: E[e^{−θT}] = ∫f e^{−θt}dt; mit t = ϖ²/(4Du²) wird daraus (2/√π)∫e^{−u²−k²ϖ²/(4u²)}du = e^{−kϖ}, k = √(θ/D)")
P("   (Integral ∫₀^∞ e^{−u²−b²/u²}du = (√π/2)e^{−2b}); Sattelpunkt: min_t[θt + ϖ²/(4Dt)] = ϖ√(θ/D) — der Exponent ist exakt.")
P("   Mit ϖ² = D·τ_D: kϖ = √(θτ_D) = √x  ⇒  1 − w = E[e^{−θ̃T}] = e^{−√x},  w_R = 1 − e^{−√x}.")
P("   Numerisch (quad über ln t, τ_D = 1, eigenes Integral):")
for x in (0.01, 0.4805, 1.0, 17.86, 100.0, 1e4, 1e6):
    # Integrand in u = ln t, Exponent zentriert am Sattel t* = 1/(2√x)
    ts = 0.5 / math.sqrt(x)
    ex0 = -math.sqrt(x)
    g = lambda u: math.exp(-x * math.exp(u) - 1 / (4 * math.exp(u)) - ex0) / math.sqrt(4 * math.pi * math.exp(u))
    lo, hi = math.log(ts) - 12, math.log(ts) + 12
    val = quad(g, lo, hi, points=[math.log(ts)], limit=400, epsabs=0, epsrel=1e-12)[0]
    P("     x = %-9g ∫ = e^{−√x}·%.10f   (Rückkehr: log10(1 − w) = %.4f, Formel %.4f)" % (x, val, (math.log(val) + ex0) / math.log(10), -math.sqrt(x) / math.log(10)))

P()
P("   Halbgerade bei Λ geschlossen (spiegelnd, gleich Kreis 2Λ): 1 − w = cosh(k(Λ−ϖ))/cosh(kΛ);  offen bei Λ (absorbierend): sinh(k(Λ−ϖ))/sinh(kΛ)")
P("   Verhältnis zur offenen Halbgeraden: cosh-Form/e^{−kϖ} = (1 + e^{−2k(Λ−ϖ)})/(1 + e^{−2kΛ}), kΛ = √x·Λ/ϖ.")
Lam = K.schliesslaenge(A0)
P("   ϖ = r*(M) (ϖ² = Dτ_D, D = GM/c), Λ = %.4f Gpc:" % (Lam / (1e9 * K.pc)))
P("   M [M☉]      ϖ/Λ          rel. Änderung von 1 − w durch den Abschluss bei x = 1 | 17,86 | 10³ | 10⁶")
def rel_cosh(x, r):
    sx = math.sqrt(x)
    return math.expm1(math.log1p(math.exp(-2 * sx * (1 / r - 1))) - math.log1p(math.exp(-2 * sx / r)))
for M in (1e8, 1e12, 1e15, 1e18, 1e21, 1e22):
    r = K.r_stern(M * M_sun) / Lam
    P("   %.0e     %.3e    %s" % (M, r, "  ".join("%+.2e" % rel_cosh(x, r) for x in (1.0, 17.86, 1e3, 1e6))))
Sb = Lam / (1 + Q_KOSM)                  # ΣΔs_b der Zelle (Λ = ν·ΣΔs_b)
rz = math.sqrt(6 * Sb * Lam) / Lam       # r* einer Masse mit Δs = ΣΔs_b: r*² = r₀Λ = 6ΣΔs_b·Λ
P("   ganze Zelle als eine Falte (Δs = ΣΔs_b = Λ/ν_kosm): ϖ/Λ = √(6/ν) = %.5f:  %s" % (rz, "  ".join("x=%g: %+.3e" % (x, rel_cosh(x, rz)) for x in (1.0, 17.86, 1e3))))
P("   ⇒ Für jeden Bereich unterhalb der Zelle (ϖ/Λ ≤ 10⁻⁴ bis 10¹⁸ M☉) ändert der Abschluss bei Λ 1 − w um weniger als e^{−2√x·Λ/ϖ}: null.")

# Zweiter Weg: Eigenfunktionsreihe auf [0, Λ] (absorbierend bei 0, spiegelnd bei Λ) und Differenzenlösung D G'' = θ G
P()
P("   Zweiter Weg (Reihe der Überlebenswahrscheinlichkeit): S(t) = Σ 4/((2n+1)π)·sin(k_nϖ)e^{−Dk_n²t}, k_n = (2n+1)π/(2Λ);")
P("   E[e^{−θT}] = 1 − θ·Σ 4/((2n+1)π)·sin(k_nϖ)/(θ + Dk_n²).  Dritter Weg: Differenzen D G'' = θG, G(0) = 1, G'(Λ) = 0.")
def reihe(xx, r, N=400000):
    # Einheiten: Λ = 1, ϖ = r, D·τ_D = ϖ² ⇒ D = r² (τ_D = 1), θ = x
    n = np.arange(N); kn = (2 * n + 1) * math.pi / 2
    return 1 - xx * np.sum(4 / ((2 * n + 1) * math.pi) * np.sin(kn * r) / (xx + r * r * kn * kn))
def differenzen(xx, r, M=20000):
    h = 1.0 / M; D = r * r
    # tridiagonal: D(G_{i+1} − 2G_i + G_{i−1})/h² − θG_i = 0, G_0 = 1, Spiegel bei M
    from scipy.linalg import solve_banded
    n = M
    ab = np.zeros((3, n)); rhs = np.zeros(n)
    main = -2 * D / h**2 - xx; off = D / h**2
    ab[1, :] = main; ab[0, 1:] = off; ab[2, :-1] = off
    ab[2, n - 2] = 2 * off                         # Spiegel G_{M+1} = G_{M−1}
    rhs[0] = -off
    Gs = solve_banded((1, 1), ab, rhs)
    i = int(round(r * M)) - 1
    return Gs[i]
P("   ϖ/Λ     x       cosh-Form      Reihe          Differenzen    offen e^{−√x}")
for r in (0.2, 0.5, 0.97):
    for xx in (0.1, 1.0, 17.86):
        ch = math.cosh(math.sqrt(xx) / r * (1 - r)) / math.cosh(math.sqrt(xx) / r)
        P("   %.2f   %6.2f   %.8f     %.8f     %.8f     %.8f" % (r, xx, ch, reihe(xx, r), differenzen(xx, r), math.exp(-math.sqrt(xx))))

# ======================================================================================== 2
hdr("2  1 − w im Hintergrund: θ̃ = 3H(z) (fällt monoton, Preisach = jetziger Wert), a₀(z) = a₀·a0_frueh_relativ(z)")
P("   z        a₀(z)/a₀   x = 3H·c/a₀(z)   1 − w_M          1 − w_R            (1−w_R)/(1−w_M)   w_R")
ZL = (0, 0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 300, 1090, 1100)
bgtab = {}
for z in ZL:
    az = a0_frueh_relativ(z)
    x = 3 * H0 * E_hintergrund(z) * c / (A0 * az)
    lM, lR = l10_1mw(x, "M"), l10_1mw(x, "R")
    bgtab[z] = (x, lM, lR)
    P("   %-7g  %.4f     %.5e      %.4e       10^%9.3f       %.3e          %s" % (z, az, x, 10**lM, lR, 10**(lR - lM), ("%.5f" % (1 - 10**lR))))
P("   mit a₀ = cH₀/2π (x = 6π): 1 − w_M = %.4e, 1 − w_R = %.4e" % (1 / (1 + 36 * math.pi**2), math.exp(-math.sqrt(6 * math.pi))))
P("   z = 1090: √x = %.3f, 1 − w_R = e^{−%.2f} = 10^%.2f (Kontrolle: 10⁻³²¹); 1 − w_M = %.2e" % (math.sqrt(bgtab[1090][0]), math.sqrt(bgtab[1090][0]), bgtab[1090][2], 10**bgtab[1090][1]))

# ======================================================================================== 3
hdr("3  Lineares Wachstum mit Übergabe λ = (1 − w)ω/N*")
P("   Lesart (03, 05): w ist der Speicheranteil — mitatmend, die Mode wirkt wie druckfreier Staub; 1 − w ist der Anteil, der nicht")
P("   mitatmet und mit der Rate λ(w) = (1 − w)·ω_orb/N* in die Operatorform übergeht (03, Zeile der Gleichung λ(w)).")
P("   Im homogenen Hintergrund zählt die Operatorform nicht (kein Gefälle, q = 1): konservativ klumpt der übergegangene Anteil nicht.")
P("   dε/dln a = (1 − w)·√(Ω_m(a)/2)·(1 − ε)/N*  (ω_orb = H√(Ω_m(a)/2)),   D'' + (2 + dlnE/dlna)D' = 1,5Ω_m(a)(1 − ε)D")
P("   Obere Schranke (sofort): Anteil (1 − w) der Mode klumpt nicht: Quelle 1,5Ω_m(a)(1 − (1 − w)·Q/(1+Q))D.")

def wachstum(Om, h, a0, form, Ns=1.0, lesart="uebergabe", a0z=True):
    OL = 1 - Om
    H0s = h * 100e3 / (1e6 * K.pc)
    def one_m_w(a):
        E = math.sqrt(Om * a**-3 + OL)
        x = 3 * H0s * E * c / (a0 * (a0_frueh_relativ(1 / a - 1) if a0z else 1.0))
        return 10**l10_1mw(x, form)
    def rhs(n, y):
        a = math.exp(n)
        E2 = Om * a**-3 + OL; Oma = Om * a**-3 / E2
        D, Dp, eps = y
        omw = one_m_w(a)
        if lesart == "uebergabe":
            deps = omw * math.sqrt(Oma / 2) * (1 - eps) / Ns if Ns else 0.0
            src = 1.5 * Oma * (1 - eps) * D
        elif lesart == "sofort":
            deps = 0.0
            src = 1.5 * Oma * (1 - omw * Q_KOSM / (1 + Q_KOSM)) * D
        else:
            deps = 0.0; src = 1.5 * Oma * D
        return [Dp, -(2 - 1.5 * Oma) * Dp + src, deps]
    # ε ab a = 10⁻⁴ (vor a = 10⁻³ vernachlässigbar), D ab 10⁻³ mit D = a
    s0 = solve_ivp(lambda n, y: [0, 0, rhs(n, [1, 1, y[2]])[2]], (math.log(1e-4), math.log(1e-3)), [0, 0, 0.0], rtol=1e-10, atol=1e-14)
    e3 = s0.y[2][-1]
    s = solve_ivp(rhs, (math.log(1e-3), 0.0), [1e-3, 1e-3, e3], rtol=1e-11, atol=1e-16, dense_output=True)
    return s

def dev(Om, h, a0, form, Ns=1.0, lesart="uebergabe", zs=(0.0, 0.5, 1.0), a0z=True):
    ref = wachstum(Om, h, a0, form, None, "keine")
    mod = wachstum(Om, h, a0, form, Ns, lesart, a0z)
    out = []
    for z in zs:
        n = math.log(1 / (1 + z))
        D0, Dp0, _ = ref.sol(n); D1, Dp1, e1 = mod.sol(n)
        out.append((z, D1 / D0 - 1, Dp1 / Dp0 - 1, e1))   # Δσ₈/σ₈ = ΔD/D; Δ(fσ₈)/fσ₈ = Δ(D')/D'
    return out

P()
P("   Kontrolle (Eingaben wie die Rechnung zur Schließungsregel: Ω_m = 0,3138, h = 0,6736, a₀ = 1,2·10⁻¹⁰, a₀ fest, N* = 1, w_M):")
Omc = (0.02237 + 0.1200) / 0.6736**2
k = dev(Omc, 0.6736, 1.2e-10, "M", 1.0, zs=(0.0, 0.5, 1.0, 2.0), a0z=False)
P("     " + "  ".join("z=%.1f: Δ(fσ₈)/fσ₈ = %+.3e" % (z, dfs) for z, dD, dfs, e in k) + "   (Ziel −4,65e-04 / −2,19e-04 / −1,09e-04 / −3,60e-05)")
KON = k[0][2]

P()
P("   Kern (Ω_m = %.4f, h = 0,674, a₀ = 1,1·10⁻¹⁰, a₀(z) aus dem Kern):" % OMEGA_M)
P("   Form  Lesart       N*     ε(z=0)       Δσ₈/σ₈ z=0 / 0,5 / 1                 Δ(fσ₈)/fσ₈ z=0 / 0,5 / 1")
WT = {}
for form in ("M", "R"):
    for les, Ns in (("uebergabe", 1.0), ("uebergabe", 0.3), ("uebergabe", 2.4), ("sofort", None)):
        r = dev(OMEGA_M, 0.674, A0, form, Ns if Ns else 1.0, les)
        WT[(form, les, Ns)] = r
        P("   w_%s   %-11s  %-5s  %.3e    %+.2e %+.2e %+.2e          %+.2e %+.2e %+.2e"
          % (form, les, "–" if Ns is None else "%.1f" % Ns, r[0][3], r[0][1], r[1][1], r[2][1], r[0][2], r[1][2], r[2][2]))
P("   (Übergabe ∝ 1/N*: die Abweichung skaliert fast genau mit 1/N*.)")
P("   Hintergrundstrahlung: 1 − w_R(1090) = 10^%.1f; Schallgeschwindigkeit c_s² = 1/(3(1 + R + (1−w)QR)) unverändert; Peaks wie bisher." % bgtab[1090][2])

# ======================================================================================== 4
hdr("4  Übergangsbreite w 0,9 → 0,1 vor der Umkehr (Preisach: θ̃ = θ bis zur Umkehr)")
Om = OMEGA_M; OL = 1 - Om
def t_von_a(a):                           # ohne Strahlung, in 1/H₀
    return 2 / (3 * math.sqrt(OL)) * math.asinh(math.sqrt(OL / Om) * a**1.5)
def alter_kern(z):                        # mit Strahlung, in s
    return quad(lambda zz: 1 / ((1 + zz) * H0 * E_hintergrund(zz)), z, np.inf, limit=400)[0]

# (a) Zykloide (EdS): R = A(1 − cos η), t = B(η − sin η), t_c = 2πB = Alter bei z_c (Kern, mit Strahlung)
def zykloide(zc, xs):
    B = alter_kern(zc) / (2 * math.pi)
    th = lambda eta: 3 * math.sin(eta) / (B * (1 - math.cos(eta))**2)
    tt = []
    for x in xs:
        e = brentq(lambda eta: th(eta) - x * A0 / c, 1e-6, math.pi)
        tt.append(B * (e - math.sin(e)))
    return tt, B
# (b) Energieintegral im Kern-Hintergrund ohne Strahlung: Ṙ² = 2μ(1/R − 1/R_ta) − Ω_L(R_ta² − R²), H₀ = 1, R in R_L
MU = Om / 2
def Rdot(R, Rta):
    return math.sqrt(max(2 * MU * (1 / R - 1 / Rta) - OL * (Rta**2 - R**2), 0.0))
def t_bis(s1, s2, Rta):                   # Zeit von R = s1·R_ta bis s2·R_ta (s1 < s2 ≤ 1), Substitution s = 1 − u²
    f = lambda u: 2 * u * Rta / max(Rdot(Rta * (1 - u * u), Rta), 1e-300)
    u1, u2 = math.sqrt(1 - s2), math.sqrt(1 - s1)
    return quad(f, u1, u2, limit=400, epsrel=1e-11)[0]
def t_ta(Rta):
    return t_bis(0.0, 1.0, Rta)
def Rta_fuer(zc):
    tc = t_von_a(1 / (1 + zc))
    return brentq(lambda R: 2 * t_ta(R) - t_bis(0.0, 0.01, R) - tc, 1e-3, 1.5, xtol=1e-13)
def theta_s(s, Rta):                      # θ/H₀ bei R = s·R_ta (Expansion)
    return 3 * Rdot(s * Rta, Rta) / (s * Rta)
def s_fuer_x(x, Rta):
    return brentq(lambda s: theta_s(s, Rta) / thc - x, 1e-6, 1 - 1e-15, xtol=1e-15)
def breite_energie(zc, x9, x1):
    Rta = Rta_fuer(zc)
    s9, s1 = s_fuer_x(x9, Rta), s_fuer_x(x1, Rta)
    return t_bis(s9, s1, Rta) / H0 / Gyr * 1e3, t_bis(s9, 1.0, Rta) / H0 / Gyr * 1e3, Rta

P("   (a) Zykloide (Einstein-de Sitter, t_c = Alter bei z_c im Kern);  (b) Energieintegral im Kern-Hintergrund (Ω_Λ, ohne Strahlung)")
P("   z_c    Maxwell (a)   (b)   [gewicht_preisach 7,24/22,48/72,00/156,40]    Rückkehr (a)   (b)   [kac 14,36/44,58/142,64/309,26]   R/M   w_R=0,9 vor Umkehr (b)")
BR = {}
for zc in (3.4, 2.0, 1.0, 0.5):
    (t9, t1), B = zykloide(zc, (X9M, X1M)); bMa = (t1 - t9) / Gyr * 1e3
    (t9r, t1r), _ = zykloide(zc, (X9R, X1R)); bRa = (t1r - t9r) / Gyr * 1e3
    bMb, _, Rta = breite_energie(zc, X9M, X1M)
    bRb, vorR, _ = breite_energie(zc, X9R, X1R)
    BR[zc] = (bMa, bMb, bRa, bRb, vorR)
    P("   %4.1f   %8.2f   %8.2f                                             %8.2f   %8.2f                                  %.3f   %8.2f" % (zc, bMa, bMb, bRa, bRb, bRb / bMb, vorR))
P("   Einsatz in frühen Dichtespitzen (Zykloide):  z_c   Maxwell [Mio. J.]   Rückkehr [Mio. J.]   (stoerungen_rekombination: 0,02/0,06/0,46/2,85)")
for zc in (30.0, 20.0, 10.0, 5.0):
    (a9, a1), _ = zykloide(zc, (X9M, X1M)); (r9, r1), _ = zykloide(zc, (X9R, X1R))
    P("                                                %5.1f    %8.3f            %8.3f" % (zc, (a1 - a9) / Gyr * 1e3, (r1 - r9) / Gyr * 1e3))
P("   Hand (Zykloide, linear um die Umkehr): θ̇_ta = −3/(8B²), Δt = Δx·θ_c·8B²/3; z_c = 1: B = %.4f Gyr" % (zykloide(1.0, (1.0,))[1] / Gyr))
B1 = zykloide(1.0, (1.0,))[1]
P("     Maxwell Δx = 8/3: Δt = %.2f Mio. J.;  Rückkehr Δx = %.4f: Δt = %.2f Mio. J. (linear; exakt oben, Rückkehr reicht weiter vor die Umkehr)"
  % (8 / 3 * (A0 / c) * 8 * B1**2 / 3 / Gyr * 1e3, X9R - X1R, (X9R - X1R) * (A0 / c) * 8 * B1**2 / 3 / Gyr * 1e3))

# ======================================================================================== 5
hdr("5  Umkehrkante heute (Energieintegral der Schale im Kern-Hintergrund ohne Strahlung)")
T0 = t_von_a(1.0)
def alter_schale(mu, Kc):                 # Zeit bis R = 1 (expandierend), Ṙ² = 2μ/R + Ω_L R² − K
    Rm = (mu / OL)**(1 / 3) if OL > 0 else 2.0
    if Rm < 1 and 2 * mu / Rm + OL * Rm * Rm - Kc <= 0:
        return 1e3                         # Schale kehrt vor R = 1 um
    f = lambda R: 1 / math.sqrt(max(2 * mu / R + OL * R * R - Kc, 1e-300))
    return quad(f, 0, 1, limit=400, epsrel=1e-12)[0]
def mu_fuer_y(y):                         # y = θ/3H₀ heute
    return brentq(lambda mu: alter_schale(mu, 2 * mu + OL - y * y) - T0, 0.05, 20, xtol=1e-14)
mu0 = mu_fuer_y(1.0)
P("   Probe: y = θ/3H₀ = 1 ⇒ 1 + Δ = 2μ/Ω_m = %.8f (soll 1)" % (2 * mu0 / Om))
# Umkehr heute: y = 0 (K = 2μ + Ω_L), R = R_ta = 1: Zeit bis zur Umkehr = t₀
def t_ta_unit(mu):                        # R_ta = 1, Ṙ² = 2μ(1/R − 1) − Ω_L(1 − R²)
    f = lambda u: 2 * u / math.sqrt(max(2 * mu * (1 / (1 - u * u) - 1) - OL * (1 - (1 - u * u)**2), 1e-300))
    return quad(f, 0, 1, limit=400, epsrel=1e-12)[0]
mu_ta = brentq(lambda mu: t_ta_unit(mu) - T0, 0.5, 50, xtol=1e-14)
D_ta = 2 * mu_ta / Om
P("   Umkehr heute: 1 + Δ_ta = %.4f (trichterkette_gewicht mit Strahlung: 11,816)" % D_ta)
P("   w       Maxwell: θ̃/3H₀   1 + Δ      Rückkehr: θ̃/3H₀   1 + Δ      (trichterkette Maxwell: 0,99→4,599, 0,9→9,345, 0,5→10,965, 0,1→11,529; kac: 7,619/11,404/11,806)")
KANTE = {}
for wv, xm, xr in ((0.99, math.sqrt(0.99 / 0.01), math.log(100)**2), (0.9, X9M, X9R), (0.5, X5M, X5R), (0.1, X1M, X1R)):
    ym, yr = xm * thc / 3, xr * thc / 3
    Dm = 2 * mu_fuer_y(ym) / Om if ym < 1 else float("nan")
    Dr = 2 * mu_fuer_y(yr) / Om if yr < 1 else float("nan")
    KANTE[wv] = (Dm, Dr)
    P("   %.2f    %.5f        %8.4f       %.5f          %s" % (wv, ym, Dm, yr, ("%8.4f" % Dr) if yr < 1 else "– (θ̃ > 3H₀: nur in Leeren)"))
P("   w_R = 0,99 verlangt θ̃ = %.3f·3H₀: nur Bereiche, die schneller als der Hintergrund atmen (Leeren)." % (math.log(100)**2 * thc / 3))
P("   w im Hintergrund: w_M = %.5f, w_R = %.5f" % (wM(3 / thc), wR(3 / thc)))
P("   1 + Δ     θ/3H₀      w_M        w_R        (Punktmasse: R/R_ta = ((1+Δ_ta)/(1+Δ))^{1/3})")
def y_fuer_D(D):
    mu = D * Om / 2
    Kc = brentq(lambda Kc: alter_schale(mu, Kc) - T0, -50, 2 * mu + OL - 1e-12, xtol=1e-14)
    return math.sqrt(2 * mu + OL - Kc)
for D in (0.3, 0.5, 1.0, 1.5, 2.0, 3.0, D_ta / 1.5**3, 5.0, 8.0, 10.0):
    y = y_fuer_D(D)
    P("   %6.3f    %.5f    %.5f    %.5f     R/R_ta = %.3f" % (D, y, wM(3 * y / thc), wR(3 * y / thc), (D_ta / D)**(1 / 3)))
P("   Schale w 0,9 → 0,1 in R/R_ta (Punktmasse): Maxwell %.3f → %.3f, Rückkehr %.3f → %.3f"
  % ((D_ta / KANTE[0.9][0])**(1 / 3), (D_ta / KANTE[0.1][0])**(1 / 3), (D_ta / KANTE[0.9][1])**(1 / 3), (D_ta / KANTE[0.1][1])**(1 / 3)))
for name, M in (("Lokale Gruppe 2·10¹²", 2e12), ("Virgo 1,2·10¹⁵", 1.2e15)):
    rta = (3 * M * M_sun / (4 * math.pi * K.rho_m_mitbewegt() * D_ta))**(1 / 3) / (1e6 * K.pc)
    P("   %-22s R_ta = %.2f Mpc; w = ½ bei %.3f (M) / %.3f Mpc (R); w = 0,9 bei %.3f (M) / %.3f Mpc (R)"
      % (name, rta, rta * (D_ta / KANTE[0.5][0])**(1 / 3), rta * (D_ta / KANTE[0.5][1])**(1 / 3), rta * (D_ta / KANTE[0.9][0])**(1 / 3), rta * (D_ta / KANTE[0.9][1])**(1 / 3)))

# ======================================================================================== 6
hdr("6  Folgen")
P("6a Ω_s/(1 − w):  x = 3Hτ_D = 6π/Θ (Θ = 2πa₀/(cH));  Ω_s = (π²/45)Θ²")
P("   Maxwell: 1 − w = 1/(1 + (6π/Θ)²) ≈ (Θ/6π)² ⇒ Ω_s/(1 − w) = (π²/45)(Θ² + 36π²) ≈ 36π⁴/45 = %.2f" % (36 * math.pi**4 / 45))
P("   Rückkehr: 1 − w = e^{−√(6π/Θ)} ⇒ Ω_s/(1 − w) = (π²/45)·Θ²·e^{√(6π/Θ)}")
P("   (Ω_s hier mit festem a₀ zur Epoche gerechnet wie im Modell: Θ(z) = 2πa₀(z)/(cH(z)))")
P("   z       Θ          Ω_s/(1−w_M)    Ω_s/(1−w_R)")
for z in (0, 0.5, 1, 2, 10, 1090):
    Th = 2 * math.pi * A0 * a0_frueh_relativ(z) / (c * H0 * E_hintergrund(z))
    x = 6 * math.pi / Th
    rM = K.K_S * Th**2 * (1 + x * x)
    lrR = math.log10(K.K_S * Th**2) + math.sqrt(x) / math.log(10)
    P("   %-6g  %.4e   %8.3f       %s" % (z, Th, rM, ("%.3f" % 10**lrR) if lrR < 6 else "10^%.1f" % lrR))
P("   heute direkt: Ω_s = %.4f, 1 − w_R = %.4e ⇒ %.3f (kac: 16,72)" % (K.omega_vakuum(A0), 10**bgtab[0][2], K.omega_vakuum(A0) / 10**bgtab[0][2]))

# ---- 6b frisch kollabierte Galaxien
P()
P("6b Frisch kollabierte Galaxien (M_b = 6·10¹⁰, Hernquist a = 1,5 kpc·6^0,3, Top-Hat (1+Q)M_b; Energieintegral statt ODE)")
Mb = 6e10 * M_sun; ab = 1.5 * kpc * 6**0.3
RL = (3 * (1 + Q_KOSM) * Mb / (4 * math.pi * K.rho_m_mitbewegt()))**(1 / 3)
def I_pre(zc, form):
    """I bis t_c: ∫(1 − w)ω dt (Expansion) + ∫ω(max(R, R_ta/2))dt (Einfall, w = 0); ω = √(Ω_m/(2R³)) in H₀."""
    Rta = Rta_fuer(zc)
    om = lambda R: math.sqrt(Om / (2 * R**3))
    def f_up(u):                                   # s = 1 − u², dt = 2u R_ta du/Ṙ
        s = 1 - u * u; R = s * Rta
        rd = Rdot(R, Rta)
        if rd <= 0: return 0.0
        x = 3 * rd / R / thc
        omw = 10**l10_1mw(x, form)
        return omw * om(R) * 2 * u * Rta / rd
    up = quad(f_up, 0, math.sqrt(1 - 1e-5), limit=800, epsrel=1e-10, points=[math.sqrt(1 - s_fuer_x(1.0, Rta))])[0]
    def f_dn(u):
        s = 1 - u * u; R = s * Rta
        rd = Rdot(R, Rta)
        return om(max(R, Rta / 2)) * 2 * u * Rta / max(rd, 1e-300)
    dn = quad(f_dn, 0, math.sqrt(1 - 0.01), limit=400, epsrel=1e-11)[0]
    tc = 2 * t_ta(Rta) - t_bis(0.0, 0.01, Rta)
    return up, dn, tc, Rta
def ratio_rar(rk, W, rv):
    rr = rk * kpc
    mb = K.hernquist_masse(Mb, ab, rr); gN = G * mb / rr**2; n_ = float(K.nu(gN))
    g = gN + (1 - W) * (n_ - 1) * gN + W * G * Q_KOSM * Mb * min(rr / rv, 1.0) / rr**2
    return g / (n_ * gN)
RQ = brentq(lambda R: (float(K.nu(G * K.hernquist_masse(Mb, ab, R) / R**2)) - 1) * K.hernquist_masse(Mb, ab, R) - Q_KOSM * Mb, 1 * kpc, 500 * kpc)
P("   R_Q = %.2f kpc (gewicht_preisach 53,09)" % (RQ / kpc))
FG = {}
for form in ("M", "R"):
    res = {}
    for zc in (3.4, 2.0, 1.0, 0.5):
        up, dn, tc, Rta = I_pre(zc, form)
        res[zc] = dict(up=up, dn=dn, tc=tc, ov=math.sqrt(Om / (2 * (Rta / 2)**3)), rv=Rta / 2 * RL)
    Iat = lambda zc, t: res[zc]["up"] + res[zc]["dn"] + res[zc]["ov"] * max(t - res[zc]["tc"], 0.0)
    cf05 = max(abs(ratio_rar(rk, 1.0, res[0.5]["rv"]) - 1) for rk in (10, 20, 30, 50))
    Ns = Iat(0.5, T0) / math.log(cf05 / 0.05)
    P("   w_%s: I vor der Umkehr (Expansion) z_c = 3,4/2/1/0,5: %s;  Einfall %s" % (form, " / ".join("%.4f" % res[z]["up"] for z in (3.4, 2.0, 1.0, 0.5)), " / ".join("%.3f" % res[z]["dn"] for z in (3.4, 2.0, 1.0, 0.5))))
    P("        N* aus der Enge der Relation (z_c = 0,5, ≤ 5 %% bei 10–50 kpc heute) = %.3f   (gewicht_preisach, Maxwell: 2,416)" % Ns)
    P("        z_c → z_obs   W = e^{−I/N*}   g/g_RAR − 1 bei 10 / 20 / 50 kpc      Δlog M_b bei festem v (10 / 50 kpc) [dex]")
    for zc, zo in ((3.4, 2.0), (2.0, 1.0), (1.0, 0.5), (0.5, 0.0)):
        W = math.exp(-Iat(zc, t_von_a(1 / (1 + zo))) / Ns)
        d = [ratio_rar(rk, W, res[zc]["rv"]) - 1 for rk in (10, 20, 50)]
        dm = []
        for rk, dd in ((10, d[0]), (50, d[2])):
            gN = G * K.hernquist_masse(Mb, ab, rk * kpc) / (rk * kpc)**2
            dm.append(math.log10(1 + dd) / (1 - float(K.chi(gN)) / 2))
        P("        %4.1f → %3.1f      %.4e      %+.4f  %+.4f  %+.4f              %+.4f / %+.4f" % (zc, zo, W, d[0], d[1], d[2], dm[0], dm[1]))
        FG[(form, zc)] = (W, d, dm)
    FG[(form, "N")] = Ns
P("   nötiger Massennullpunkt (3σ, z_c = 3,4 → 2, 10 kpc) = Δ/3: Maxwell %.4f dex, Rückkehr %.4f dex"
  % (abs(FG[("M", 3.4)][2][0]) / 3, abs(FG[("R", 3.4)][2][0]) / 3))
NsR = FG[("R", "N")]; NsM = FG[("M", "N")]
P("   1/λ_m heute bei 200–500 ρ_krit mit N* = 0,3–N*(Relation): Maxwell %.2f–%.2f Gyr, Rückkehr %.2f–%.2f Gyr"
  % (K.uebergabezeit(500, 0.0, (0.3, NsM))[1] / Gyr, K.uebergabezeit(200, 0.0, (0.3, NsM))[2] / Gyr,
     K.uebergabezeit(500, 0.0, (0.3, NsR))[1] / Gyr, K.uebergabezeit(200, 0.0, (0.3, NsR))[2] / Gyr))

# ---- 6c Übergabe im Hintergrund, λ = (1 − w)ω_orb/N*
P()
P("6c Übergabe im Hintergrund λ/H = (1 − w)√(Ω_m(a)/2)/N* und aufsummiert ε (N* = 1, Kern):")
for z in (0, 0.5, 1, 2, 10):
    a = 1 / (1 + z); Oma = Om * a**-3 / (Om * a**-3 + OL)
    x = bgtab[z][0] if z in bgtab else None
    P("   z = %-4g λ/H: w_M %.3e, w_R %.3e" % (z, 10**bgtab[z][1] * math.sqrt(Oma / 2), 10**bgtab[z][2] * math.sqrt(Oma / 2)))
P("   ε(z = 0) mit N* = 1: w_M %.3e, w_R %.3e;  mit N* aus der Relation (%.2f / %.2f): w_M %.3e, w_R %.3e"
  % (WT[("M", "uebergabe", 1.0)][0][3], WT[("R", "uebergabe", 1.0)][0][3], NsM, NsR,
     WT[("M", "uebergabe", 1.0)][0][3] / NsM, WT[("R", "uebergabe", 1.0)][0][3] / NsR))
rM = dev(OMEGA_M, 0.674, A0, "M", NsM); rR = dev(OMEGA_M, 0.674, A0, "R", NsR)
P("   Δ(fσ₈)/fσ₈ mit N* aus der Relation: w_M (N* = %.2f) %s;  w_R (N* = %.2f) %s"
  % (NsM, " ".join("z=%.1f: %+.2e" % (z, x) for z, _, x, _ in rM), NsR, " ".join("z=%.1f: %+.2e" % (z, x) for z, _, x, _ in rR)))
P("   Δσ₈/σ₈ mit N* aus der Relation: w_M %s;  w_R %s" % (" ".join("%+.2e" % x for _, x, _, _ in rM), " ".join("%+.2e" % x for _, x, _, _ in rR)))

# ---- 6d Filament, w = ½, Einsatz
P()
P("6d w = ½ und Filament")
P("   w_M = ½ bei θ̃ = θ_c;  w_R = ½ bei θ̃ = (ln 2)²θ_c = %.4f θ_c = %.4f H₀;  w_R(θ_c) = 1 − 1/e = %.4f" % (X5R, X5R * thc, wR(1.0)))
P("   Filament (eine Achse mit H₀, zwei Querachsen r): θ = H₀ + 2r; w = ½ ⇒ r = (θ_½ − H₀)/2: Maxwell %.3f H₀, Rückkehr %.3f H₀"
  % ((thc - 1) / 2, (X5R * thc - 1) / 2))
P("   Wand (zwei Achsen mit H₀): r = θ_½ − 2H₀: Maxwell %.3f H₀, Rückkehr %.3f H₀" % (thc - 2, X5R * thc - 2))

# ---- 6e Linse um isolierte Galaxien (linsen_isoliert, w ersetzt)
P()
P("6e Linse um isolierte Galaxien (Milchstraße 6·10¹⁰, z = 0, Ansatz lokal; Funktionen aus linsen_isoliert, w_R statt gewicht_w):")
import linsen_isoliert as LI
w_alt = K.gewicht_w
def run_li(tag):
    LI._cache.clear()
    Mmw = 6e10 * M_sun
    for prof in ("flach", "NFW"):
        row = []
        for Ns in (0.3, 2.4, NsR):
            for rk in (200, 300):
                r = np.array([rk * kpc])
                row.append(LI.g_3D(Mmw, r, Ns=Ns, prof=prof, z=0.0, form="lokal")[0] / LI.g_3D(Mmw, r)[0] - 1)
        P("   %s %-5s: N* = 0,3: %+.3f %+.3f | N* = 2,4: %+.3f %+.3f | N* = %.2f: %+.3f %+.3f   (200/300 kpc)" % ((tag, prof, row[0], row[1], row[2], row[3], NsR, row[4], row[5])))
run_li("w_M")
K.gewicht_w = lambda th, a0=A0: -math.expm1(-math.sqrt(max(th, 0.0) * c / a0))
try:
    run_li("w_R")
finally:
    K.gewicht_w = w_alt
P("   (linsen_isoliert.txt, Maxwell: flach −0,031/−0,192, NFW −0,072/−0,394 bei N* = 2,4)")

# ======================================================================================== 7
hdr("7  Änderungsliste (alt → neu mit w_R)")
L = [
 ("03:13", "Filament mitatmend, bis Querschnitt schneller als 0,42·H₀ schrumpft (w < ½)", "0,42·H₀", "%.2f·H₀" % ((1 - X5R * thc) / 2)),
 ("03:37", "w wechselt in dünner Schale vor der Umkehr, heute bei 10- bis 12-facher Dichte", "1+Δ 9,3–11,5 (w 0,9→0,1)", "1+Δ %.1f–%.1f" % (KANTE[0.9][1], KANTE[0.1][1])),
 ("03:37", "außerhalb w ≥ 0,99", "w ≥ 0,99", "w_R ≥ %.3f nur im Hintergrund, bei 1+Δ = 3,5 (1,5 R_ta) w_R = %.3f" % (wR(3 / thc), wR(3 * y_fuer_D(D_ta / 1.5**3) / thc))),
 ("03:51", "w = θ̃²/(θ̃² + θ_c²)", "Maxwell", "w = 1 − e^{−√(θ̃/θ_c)}"),
 ("03:56", "1 − w heute", "3,1·10⁻³", "%.2e" % 10**bgtab[0][2]),
 ("03:56", "1 − w mit cH₀/2π", "2,8·10⁻³", "%.2e" % math.exp(-math.sqrt(6 * math.pi))),
 ("03:56", "1 − w Rekombination", "3·10⁻¹²", "10^%.0f" % bgtab[1090][2]),
 ("03:56", "fσ₈ auf 5·10⁻⁴ (N* = 1, a₀ = 1,2)", "%.1e (Kontrolle)" % abs(KON), "%.1e (N* = 1, a₀ = 1,1); %.1e (N* = 0,3)" % (abs(WT[("R", "uebergabe", 1.0)][0][2]), abs(WT[("R", "uebergabe", 0.3)][0][2]))),
 ("03:57", "Übergangsbreiten 7/22/72/156 Mio. J.", "7/22/72/156", "/".join("%.0f" % BR[z][3] for z in (3.4, 2.0, 1.0, 0.5))),
 ("03:57", "N* ≈ 0,3–2,4 (obere Grenze aus der Enge der Relation)", "2,4", "%.2f" % NsR),
 ("03:57", "1/λ_m heute 0,28–3,5 Mrd. J.", "0,28–3,5", "%.2f–%.2f" % (K.uebergabezeit(500, 0.0, (0.3, NsR))[1] / Gyr, K.uebergabezeit(200, 0.0, (0.3, NsR))[2] / Gyr)),
 ("03:74", "1 − w = 3·10⁻¹² bei der Rekombination", "3·10⁻¹²", "10^%.0f" % bgtab[1090][2]),
 ("03:105", "1 − w ≈ (Θ/6π)², Ω_s = 78(1 − w) in jeder Epoche", "78, fest", "(π²/45)Θ²e^{√(6π/Θ)}, heute %.1f" % (K.omega_vakuum(A0) / 10**bgtab[0][2])),
 ("03:131", "frisch kollabiert: +1,9 % (3,4→2, 10 kpc)", "+1,9 %", "%+.1f %%" % (100 * FG[("R", 3.4)][1][0])),
 ("03:131", "frisch kollabiert: −2 bis −6 % (1→0,5)", "−2 bis −6 %", "%+.1f bis %+.1f %%" % (100 * FG[("R", 1.0)][1][0], 100 * FG[("R", 1.0)][1][2])),
 ("03:131", "0,013 dex (10 kpc) bis 0,05 dex (50 kpc); Nullpunkt 0,004 dex", "0,013/0,05/0,004", "%.3f/%.3f/%.4f" % (abs(FG[("R", 3.4)][2][0]), abs(FG[("R", 1.0)][2][1]), abs(FG[("R", 3.4)][2][0]) / 3)),
 ("05:89", "w = x²/(1 + x²) als Speicheranteil eines Maxwell-Elements", "Maxwell", "Rückkehranteil 1 − e^{−√x} (Lévy-Dichte, Kac)"),
 ("07:35", "w = ½ bei θ = θ_c (gesetzt, θ_c = 1/τ_D)", "½ bei θ_c", "½ bei (ln 2)²θ_c = 0,48θ_c; θ_c = 1/τ_D folgt"),
]
for a, b, o, n in L:
    P("   %-7s %-72s  alt: %-26s neu: %s" % (a, b, o, n))
P("   03:137 Linse um isolierte Galaxien (N* = 2,4: 3–7 %% bei 200, 19–43 %% bei 300 kpc) → siehe 6e mit N* = %.2f" % NsR)
P()
P("   Prüfstand (ftm_pruefstand.py) und Gesetze (ftm_gesetze.py):")
PL = [
 ("pruefstand:69", "1 − w heute (Standardwert a₀)", "3,1e-3", "%.3e" % 10**bgtab[0][2]),
 ("pruefstand:70", "1 − w heute (a₀ = cH₀/2π)", "2,8e-3", "%.3e" % math.exp(-math.sqrt(6 * math.pi))),
 ("pruefstand:312-315", "w 0,9 → 0,1 bei z_c = 3,4/2/1/0,5 (liest gewicht_preisach.txt)", "7/22/72/156", "/".join("%.1f" % BR[z][3] for z in (3.4, 2.0, 1.0, 0.5))),
 ("pruefstand:316", "Preisach = jetzige Form im Hintergrund", "1,0", "1,0 (formunabhängig)"),
 ("pruefstand:353", "1 − w bei der Rekombination", "3,4e-12", "10^%.2f (Gleitkomma subnormal: log10 prüfen)" % bgtab[1090][2]),
 ("pruefstand:491", "w = ½ bei θ = a₀/c", "0,5", "%.4f (½ erst bei (ln 2)²a₀/c)" % wR(1.0)),
 ("pruefstand:494", "Ω_s/(1 − w) heute", "78", "%.2f" % (K.omega_vakuum(A0) / 10**bgtab[0][2])),
 ("pruefstand:562", "Ω_s/(1 − w) = 36π⁴/45", "78", "entfällt (nicht fest)"),
 ("pruefstand:569", "Preisach Δlog M_b (3,4→2, 10 kpc) [dex] (liest fruehe_galaxien_a0.txt)", "0,013", "%.4f" % abs(FG[("R", 3.4)][2][0])),
 ("pruefstand:570", "Preisach nötiger Massennullpunkt [dex]", "0,004", "%.4f" % (abs(FG[("R", 3.4)][2][0]) / 3)),
 ("pruefstand:580", "Umkehr heute 1 + Δ_ta", "11,8", "11,8 (unverändert)"),
 ("gesetze:666", "w = Speicheranteil eines Maxwell-Elements (sim M)", "trägt", "entfällt für w (Maxwell bleibt für den Kern in 05)"),
]
for a, b, o, n in PL:
    P("   %-20s %-68s  alt: %-12s neu: %s" % (a, b, o, n))

# ======================================================================================== 8
hdr("8  Handprüfungen")
P("   1) x heute = 3H₀c/a₀ = 3·2,1843e-18·2,9979e8/1,1e-10 = %.4f; √x = %.4f; e^{−√x} = %.5f; Maxwell 1/(1+x²) = %.5f"
  % (3 * H0 * c / A0, math.sqrt(3 * H0 * c / A0), math.exp(-math.sqrt(3 * H0 * c / A0)), 1 / (1 + (3 * H0 * c / A0)**2)))
P("      Von Hand: 17,859² = 318,9; 1/319,9 = 3,126e-3; √17,859 = 4,2260; e^{−4,226} = e^{−4}·e^{−0,226} = 0,018316·0,7977 = 0,014610")
P("   2) Ω_s/(1 − w_R) heute = 0,2443·e^{4,226} = 0,2443·68,44 = 16,72")
P("   3) Rekombination: x = 5,4526e5, √x = 738,42, log10 e^{−738,42} = −738,42·0,43429 = −320,69")
P("   4) Filament: θ_½ = 0,4805·0,16798 H₀ = 0,0807 H₀; r = (0,0807 − 1)/2 = −0,4596 H₀")
P("   5) Umkehrkante EdS-Probe: Energieintegral mit Ω_Λ = 0 muss 9π²/16 = 5,5517 geben:")
OL_s = OL; OL = 0.0; Om_s = Om
T0_eds = 2 / 3.0
def t_ta_eds(mu):
    f = lambda u: 2 * u / math.sqrt(max(2 * mu * (1 / (1 - u * u) - 1), 1e-300))
    return quad(f, 0, 1, limit=400, epsrel=1e-12)[0]
mu_e = brentq(lambda mu: t_ta_eds(mu) - T0_eds, 0.5, 50, xtol=1e-14)
OL = OL_s
P("      1 + Δ_ta (EdS, Ω_m = 1) = 2μ = %.5f (9π²/16 = %.5f)" % (2 * mu_e, 9 * math.pi**2 / 16))

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
