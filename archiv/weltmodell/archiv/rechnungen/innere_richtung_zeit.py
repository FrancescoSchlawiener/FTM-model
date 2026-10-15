"""Innere Richtung als imaginäre Zeit?

Frage: Die innere Richtung ist ein Kreis vom Umfang 2Λ ohne Umlauf (spiegelsymmetrisch die
bei Λ geschlossene Halbgerade), Laufzeit τ_D = 2Λ/c. Ist sie eine euklidische (imaginäre) Zeit
der Periode β = 2Λ/c, und was folgt daraus?

Blickrichtung: T_GH je Beobachter und Epoche (G, hängt an H); T_I innen (I, hängt an Λ);
T_Unruh je beschleunigtem Ort (L).
Rechnet nur mit ftm_kern. Ausgabe: innere_richtung_zeit.txt
"""
import math
import numpy as np
from scipy.optimize import brentq

import ftm_kern as K

out = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


a0 = K.A0
Lam = K.schliesslaenge(a0)          # Λ = c²/(2a₀)
tauD = K.tau_D(a0)                  # c/a₀ = 2Λ/c
H0 = K.H0_si()


def T_I(a0_=a0):
    """Temperatur der inneren Richtung: Periode 2Λ/c in imaginärer Zeit, kT = ħc/(2Λ) = ħa₀/c."""
    return K.hbar * c_ / (2 * K.schliesslaenge(a0_) * K.kB)


def T_GH(H):
    return K.hbar * H / (2 * math.pi * K.kB)


def T_U(a):
    return K.hbar * a / (2 * math.pi * K.c * K.kB)


c_ = K.c

p("=" * 78)
p("0  Grundgrößen")
p(f"a₀ = {a0:.3e} m/s²,  Λ = c²/(2a₀) = {Lam:.4e} m = {Lam/K.Gpc if hasattr(K,'Gpc') else Lam/(1e9*K.pc):.3f} Gpc")
p(f"τ_D = c/a₀ = 2Λ/c = {tauD:.4e} s = {tauD/K.Gyr:.2f} Mrd. J.;   2π/H₀ = {2*math.pi/H0/K.Gyr:.2f} Mrd. J.")

# ------------------------------------------------------------------ 1
p("=" * 78)
p("1  T_I = ħc/(2Λk) = ħa₀/(ck) gegen T_GH = ħH/(2πk);  T_I/T_GH = 2πa₀/(cH) = Θ")
TI0 = T_I()
p(f"T_I (a₀ fest) = {TI0:.4e} K;   T_GH(H₀) = {T_GH(H0):.4e} K;   Θ₀ = {K.theta_(a0, H0):.4f}")
p(f"Probe: T_I/T_GH − Θ = {TI0/T_GH(H0) - K.theta_(a0, H0):.1e}")


def E(z):
    return K.E_hintergrund(z)


p("")
p(f"{'Epoche':<28}{'H/H₀':>10}{'a₀(z)/a₀':>10}{'T_I [K]':>12}{'T_GH [K]':>12}{'Θ=T_I/T_GH':>12}")
for z in [0.0, 1.0, 10.0, 1090.0]:
    Hz = H0 * E(z)
    a0z = a0 * K.a0_frueh_relativ(z)
    TI = K.hbar * a0z / (K.c * K.kB)
    p(f"{'z = %g' % z:<28}{E(z):>10.4g}{a0z/a0:>10.4f}{TI:>12.4e}{T_GH(Hz):>12.4e}{TI/T_GH(Hz):>12.4g}")
# Zukunft
Ls = K.lambda_vakuum(Lam)
Hs = math.sqrt(Ls * K.c**2 / 3)
Hinf = H0 * math.sqrt(K.OMEGA_DE)
for name, Hz in [("Zukunft H → H_s (Vakuumteil)", Hs), ("Zukunft H → H₀√Ω_DE", Hinf)]:
    p(f"{name:<28}{Hz/H0:>10.4g}{1.0:>10.4f}{TI0:>12.4e}{T_GH(Hz):>12.4e}{TI0/T_GH(Hz):>12.4g}")
p(f"Θ_s = √45/π = {math.sqrt(45)/math.pi:.4f} = 1/√k_s (exakt, unabhängig von a₀)")
# Kreuzung Θ = 1
zx = brentq(lambda z: K.theta_(a0 * K.a0_frueh_relativ(z), H0 * E(z)) - 1, 0, 2)
p(f"Θ = 1 bei z = {zx:.4f};  |Θ − 1| ≤ 10 % für z ≤ "
  f"{brentq(lambda z: K.theta_(a0*K.a0_frueh_relativ(z), H0*E(z)) - 0.9, 0, 3):.3f} (Vergangenheit)")
tz = brentq(lambda y: K.theta_(a0, H0 * math.sqrt(K.OMEGA_M * y**-3 + K.OMEGA_R * y**-4 + K.OMEGA_DE)) - 1.1, 1, 5)
p(f"Θ = 1,1 in der Zukunft bei Skalenfaktor a = {tz:.3f}")
p("Blickrichtung: T_GH gehört zum Horizont des Beobachters (G, fällt mit H); T_I zur inneren")
p("Richtung (I, fest wie Λ). Gleich sind sie nur, solange Θ ≈ 1.")

# ------------------------------------------------------------------ 2
p("=" * 78)
p("2  Folgt ν = 1 + n_B(s) aus T_I?  Nötig: ħω/(kT_I) = s  ⇔  ω* = s·a₀/c = √(a_N a₀)/c = √(Γ/τ_D)")
M = 6e10 * K.M_sun
rs = K.r_stern(M, a0)
vinf = K.v_flach(M, a0)
p(f"Milchstraße 6·10¹⁰ M☉: r* = {rs/(1e3*K.pc):.3f} kpc, v∞ = {vinf/1e3:.1f} km/s")
p(f"{'Ort':<14}{'s':>7}{'ω* [1/s]':>11}{'Γ/ω*':>9}{'Ω/ω*':>9}{'(c/r)/ω*':>11}{'(a/c)/ω*':>10}{'(a₀/c)/ω*':>11}{'H₀/ω*':>9}")
orte = [("MW 0,3 r*", M, 0.3 * rs), ("MW r*", M, rs), ("MW 3 r*", M, 3 * rs), ("MW 30 r*", M, 30 * rs),
        ("Sonne 1 AE", K.M_sun, K.AU)]
for name, m, r in orte:
    aN = K.G * m / r**2
    s = float(K.s_(aN, a0))
    w = s * a0 / K.c
    Gam = K.verlust_gamma(aN)
    a = float(K.a_regel(aN, a0))
    Om = math.sqrt(a / r)
    p(f"{name:<14}{s:>7.3g}{w:>11.3e}{Gam/w:>9.3g}{Om/w:>9.3g}{(K.c/r)/w:>11.3e}{(a/K.c)/w:>10.4g}{(a0/K.c)/w:>11.3g}{H0/w:>9.3g}")
p("Γ/ω* = s, (a₀/c)/ω* = 1/s: keine der Raten ist konstant ω*, außer a/c im Tiefbereich (a → √(a_N a₀)).")
p("a/c = ω* gilt nur für s → 0; exakt ist (a/c)/ω* = s·ν = s/(1 − e^{−s}) (Newton-Grenze → s).")
p("Sattelpunkt der Rückkehr t* = ½√(τ_D/Γ) gibt ω* = 1/(2t*): hängt selbst an s, keine Eigenfrequenz des Orts.")
p("")
p("Räumliche Lesung (Welle in ϖ mit c, Mode k = 1/r der Kreisbahn): s = k·ϖ_R = r*/r (Punktmasse),")
p("also Periode r* in ϖ, kT_r* = ħc/r*, nicht 2Λ.  r*² = r₀·Λ (r₀ = 6Δs):")
r0 = K.r0(M)
p(f"  r* = {rs:.5e} m,  √(r₀Λ) = {math.sqrt(r0*Lam):.5e} m;  T_r* = {K.hbar*K.c/(rs*K.kB):.3e} K,  T_r*/T_I = 2Λ/r* = {2*Lam/rs:.3e}")
p(f"  Sonne: T_r*/T_I = {2*Lam/K.r_stern(K.M_sun, a0):.3e}  (massenabhängig ∝ M^(−1/2))")

# ------------------------------------------------------------------ 3
p("=" * 78)
p("3  Unruh: T_U(a) = ħa/(2πck) = T_I  ⇔  a_U = 2π·a₀")
aU = brentq(lambda a: T_U(a) - TI0, 1e-12, 1e-6)
p(f"a_U = {aU:.4e} m/s² = {aU/a0:.4f}·a₀ (2π = {2*math.pi:.4f})")
p(f"T_U(a₀) = {T_U(a0):.4e} K = T_I/{TI0/T_U(a0):.4f}")
p(f"Hubble-Kugel-Neigung cH₀/2 = {K.c*H0/2:.4e} = {K.c*H0/2/a0:.4f}·a₀ = (π/Θ₀)·a₀;  a_U/(cH₀/2) = 2Θ₀ = {aU/(K.c*H0/2):.4f}")
p(f"T_U(cH) = T_GH exakt; a_U = Θ·cH₀ = {aU/(K.c*H0):.4f}·cH₀")
aN = aU
p(f"Die Regel bei a_N = 2πa₀: s = {float(K.s_(aN,a0)):.4f}, ν = {float(K.nu(aN,a0)):.4f}, χ = {float(K.chi(aN,a0)):.4f} (kein ausgezeichneter Punkt)")
p("a₀ selbst ist kT_I·c/ħ = c/τ_D (Rate der Rückkehr), nicht die Unruh-Beschleunigung von T_I (Faktor 2π).")
p(f"Unruh-Lesung mit ω = Γ = a_N/c: ħΓ/(kT_I) = a_N/a₀ = s², n_B(s²): bei a_N = a₀ ν = {1+float(K.n_B(1.0)):.4f} (gleich),")
for x in [0.1, 0.01]:
    p(f"   a_N = {x}a₀: ν(s²) = {1+float(K.n_B(x)):.3f} gegen ν(s) = {float(K.nu(x*a0,a0)):.3f};  a → a_N + a₀ statt √(a_N a₀)")

# ------------------------------------------------------------------ 4
p("=" * 78)
p("4  Kac und Periode: Rückkehr-Laplace-Transformierte  E[e^{−ΓT}], k = √(Γ/D)")


def L_intervall(k, w, Lm):
    """[0, Λ], Ort absorbierend, Λ reflektierend."""
    return np.cosh(k * (Lm - w)) / np.cosh(k * Lm)


def L_kreis(k, w, Lm):
    """Kreis 2Λ ohne Drift, Ort absorbierend von beiden Seiten."""
    return (np.sinh(k * (2 * Lm - w)) + np.sinh(k * w)) / np.sinh(2 * k * Lm)


def L_thermisch(k, w, Lm):
    """normierter euklidischer Propagator G_β(τ)/G_β(0), β = 2Λ (in ϖ), ω = k:
    [(1 + n)e^{−kϖ} + n e^{kϖ}]/(1 + 2n), n = n_B(2kΛ)."""
    n = float(K.n_B(2 * k * Lm))
    return ((1 + n) * np.exp(-k * w) + n * np.exp(k * w)) / (1 + 2 * n)


def L_kreis_drift(k, w, Lm, Pe, D=1.0):
    """Kreis 2Λ mit Drift v (Umlauf), Pe = v·ϖ_R/D (ϖ_R = 1): D u'' + v u' − Γ u = 0, u(0) = u(2Λ) = 1."""
    v = Pe * D
    Gam = k * k * D
    disc = math.sqrt(v * v + 4 * D * Gam)
    rp, rm = (-v + disc) / (2 * D), (-v - disc) / (2 * D)
    Lt = 2 * Lm
    # u = A e^{rp ϖ} + B e^{rm ϖ}, numerisch stabil über Verschiebung
    # Randwerte: A + B = 1, A e^{rp Lt} + B e^{rm Lt} = 1
    ep, em = math.exp(rp * Lt) if rp * Lt < 700 else float('inf'), math.exp(rm * Lt)
    if math.isinf(ep):
        # A klein: A = (1 − B e^{rm Lt}) e^{−rp Lt}, B ≈ 1
        B = 1.0
        return math.exp(rm * w) + math.exp(rp * (w - Lt)) * (1 - em)
    A = (1 - em) / (ep - em)
    B = 1 - A
    return A * math.exp(rp * w) + B * math.exp(rm * w)


dmax = 0.0
for kL in [0.1, 1.0, 3.0, 10.0]:
    for x in [0.1, 0.5, 0.9]:
        a_, b_, t_ = L_intervall(kL, x, 1.0), L_kreis(kL, x, 1.0), L_thermisch(kL, x, 1.0)
        dmax = max(dmax, abs(a_ - b_), abs(a_ - t_))
p(f"Intervall [0,Λ] reflektierend = Kreis 2Λ ohne Drift = euklidischer Propagator mit β = 2Λ: max. Abweichung {dmax:.1e}")
p(f"Beispiel kΛ = 1, ϖ = Λ/2: {L_intervall(1.0,0.5,1.0):.6f} / {L_kreis(1.0,0.5,1.0):.6f} / {L_thermisch(1.0,0.5,1.0):.6f}")
p("Also: streng periodisch (Kreis) ändert gegenüber reflektierend nichts, solange kein Umlauf (Drift).")
p("")
p("Rückkehranteil mit ϖ_R = r*, k = 1/r (Punktmasse): q_Λ = cosh(s(Λ/r* − 1))/cosh(sΛ/r*)")
p("   = e^{−s} + thermische Beimischung mit n_B(2Λ/r) = n_B(T_r/T_I)-artig, T_r = ħc/(kr)")
for Lr in [10.0, 100.0]:
    p(f" Λ/r* = {Lr:g}:")
    for s in [3.0, 1.0, 0.3, 0.1, 0.03, 0.01]:
        qL = float(L_intervall(s, 1.0, Lr))
        n = float(K.n_B(2 * s * Lr))
        p(f"   s = {s:<5g} 1−q_Λ = {1-qL:.4e}  (1−q_∞ = {-math.expm1(-s):.4e}),  ν_Λ/ν_∞ = {(-math.expm1(-s))/(1-qL):.4f},  n_B(2sΛ/r*) = {n:.3e}")
Lr = Lam / rs
x100 = 2 * Lam / (100e3 * K.pc)
p(f" Milchstraße heute: Λ/r* = {Lr:.3e}; bei 100 kpc 2Λ/r = {x100:.3e}, n_B ≈ e^(−{x100:.3e}) (unmessbar)")
p(f" 10 % Abweichung (Kac-Budget) erst bei r ≈ 0,657·Λ = {0.657*Lam/K.Gpc:.2f} Gpc, dort 2Λ/r = {2/0.657:.2f}, n_B = {float(K.n_B(2/0.657)):.3f}")
p("")
p("Mit Umlauf (Drift v, Pe = vϖ_R/D) auf dem Kreis 2Λ, Λ/ϖ_R = 1000: lokale Steigung d ln(1−q)/d ln s")
for Pe in [0.0, 0.01, 0.1, 1.0, 10.0]:
    row = []
    for s in [1.0, 0.1, 0.01, 0.001]:
        f = lambda ss: 1 - L_kreis_drift(ss, 1.0, 1000.0, Pe)
        h = 1.02
        slope = (math.log(f(s*h)) - math.log(f(s/h))) / (2*math.log(h))
        row.append(f"{slope:6.3f}")
    p(f"  Pe = {Pe:<5g}  s = 1 / 0,1 / 0,01 / 0,001:  " + "  ".join(row))
p("Pe = 0: Wurzel (Steigung 1) bis s ≈ ϖ_R/Λ, danach 2 (Kac-Budget). Drift bricht die Wurzel ab s ≲ Pe.")
p("Kac: mittlere Rückkehr 2Λ/c ist bei Diffusion ein Mittel über eine breite Verteilung (∝ T^{−3/2}),")
p("bei Umlauf wäre sie eine feste Periode – dann aber keine Wurzel. Die Periode 2Λ/c ist Geometrie (Umfang), kein Takt.")

# ------------------------------------------------------------------ 5
p("=" * 78)
p("5  Wick-Drehung: □δσ = (−∂_t²/c² + ∇²)δσ, t = −iϖ/c  ⇒  (∇² + ∂_ϖ²)δσ = 0 (erweiterte Faltungsgleichung)")
# Probe numerisch: f = e^{ikx} e^{−kϖ} erfüllt Laplace; Welle e^{ikx}e^{−iωt} mit ω = ck geht bei t = −iϖ/c in e^{ikx}e^{−kϖ} über
k = 1.3; x = 0.4; w = 0.7
lap = (-(k**2) + k**2) * np.exp(1j*k*x) * np.exp(-k*w)
t = -1j * w / K.c
wave = np.exp(1j*k*x) * np.exp(-1j*(K.c*k)*t)
p(f"Probe: Laplace-Rest {abs(lap):.1e};  Welle bei t = −iϖ/c minus Glättungsmode: {abs(wave - np.exp(1j*k*x)*np.exp(-k*w)):.1e}")
p("Die Glättung e^{−|k|ϖ} ist die Welle ω = c|k| in imaginärer Zeit ϖ/c. Reelle Zeit (hyperbolisch) im Raum,")
p("imaginäre (elliptisch) in der Rückkehr; V1 der Halbgeraden (ϖ raumartig) ist genau diese Signatur.")
p("Bahntorus: q = e^{2πiτ_M} = e^{2πi·Re τ_M}·e^{−s}: Phase aus dem Umlauf (reelle Zeit), Betrag aus der Rückkehr.")

with open(__file__.replace(".py", ".txt"), "w") as f:
    f.write("\n".join(out) + "\n")
