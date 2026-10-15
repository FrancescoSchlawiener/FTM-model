# -*- coding: utf-8 -*-
"""
Kern der Kopfgleichung (07) aus derselben Rückkehr wie Regel und Gewicht.

Frage: 07 schreibt den Kern noch als Maxwell-Relaxation e^{−λ_K(t−t′)}; Regel (q = e^{−s}) und Gewicht
(w = 1 − e^{−√(θ̃τ_D)}) kommen aus der Dichte erster Rückkehrzeiten. Welcher Zeitkern folgt aus ihr, ist
λ_K = Γ = a_N/c, und was wird aus der Schranke λ_K ≲ 3·10⁻⁴ Ω (kern_gedaechtnis.py)?

1  Dichte erster Rückkehrzeiten (Lévy, einseitig stabil ½):  f(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/(4T)},
   E[e^{−pT}] = e^{−√(pτ_D)}  (numerisch gegen geschlossen).
2  Mit Verlust Γ:  K(T) = f(T)·e^{−ΓT},  K̂(p) = e^{−√((p+Γ)τ_D)},  K̂(0) = e^{−s} = q.
   Summe über Rückläufe ν̂(p) = 1/(1 − K̂(p)); tiefer Bereich: ν̂ ≈ 1/√((p+Γ)τ_D), gedämpfte Halbintegration.
3  Sonnenkreis (Punktmasse 6·10¹⁰ M☉, 8 kpc, wie kern_gedaechtnis.py): Γ/Ω zwei Wege (Kern; Hand v_N²/(c·v)).
4  Gewicht des Kerns bei Bahnfrequenz: |K̂(iΩ)|/K̂(0) gegen Maxwell λ/√(λ²+Ω²) (zwei Wege: komplexe Wurzel; geschlossene Form).
5  Anteil des Kerns innerhalb t (Weltalter, Bahnperiode): numerisch und geschlossen (erfc-Formel).
"""
import sys, io, math, cmath
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import *

out = []
def P(s=""):
    out.append(s); print(s)

TD = tau_D()

def f_levy(T, td=1.0):
    return math.sqrt(td / (4 * math.pi)) * T**-1.5 * math.exp(-td / (4 * T)) if T > 0 else 0.0

# ---------------------------------------------------------------- 1
P("1) Dichte erster Rückkehrzeiten, τ_D = 1:  E[e^{−pT}] numerisch gegen e^{−√p}")
for p in (0.01, 0.3, 1.0, 3.0):
    num = quad(lambda T: f_levy(T) * math.exp(-p * T), 0, np.inf, limit=400)[0]
    P(f"   p = {p:5.2f}:  numerisch {num:.8f}   geschlossen {math.exp(-math.sqrt(p)):.8f}")
norm = quad(f_levy, 0, np.inf, limit=400)[0]
P(f"   Norm ∫f = {norm:.6f} (Rückkehr sicher);  Modus T = τ_D/6;  Mittel unendlich (Schwanz T^{{−3/2}})")
P(f"   τ_D = c/a₀ = {TD:.4e} s = {TD/Gyr:.1f} Mrd. Jahre;  Modus τ_D/6 = {TD/6/Gyr:.1f} Mrd. Jahre")

# ---------------------------------------------------------------- 2
P()
P("2) Kern mit Verlust Γ:  K(T) = f(T)e^{−ΓT},  K̂(p) = e^{−√((p+Γ)τ_D)}")
for s in (0.3, 1.0, 3.0):
    G_ = s * s     # Γτ_D = s²
    k0 = quad(lambda T: f_levy(T) * math.exp(-G_ * T), 0, np.inf, limit=400)[0]
    P(f"   s = {s}:  K̂(0) numerisch {k0:.6f},  e^{{−s}} = {math.exp(-s):.6f}")
P("   Summe über Rückläufe ν̂(p) = 1/(1 − K̂(p)); für s → 0: ν̂·√((p+Γ)τ_D) →")
for x in (0.3, 0.1, 0.01):
    P(f"      √((p+Γ)τ_D) = {x}:  {1/(-math.expm1(-x))*x:.5f}  (1 + x/2 + …)")
P("   Tiefbereich: ν̂ ≈ ((p+Γ)τ_D)^{−1/2} ⇒ Zeitkern e^{−Γt}·t^{−1/2}/√(πτ_D): gedämpfte Halbintegration")
P("   (d/dt + Γ)^{−1/2}; exakt ist K̂ = exp[−√τ_D·(d/dt + Γ)^{1/2}], keine Mittag-Leffler-Form.")
P("   Der Exponentialfaktor des Kerns ist genau der Verlust: λ_K = Γ, aber mal der Lévy-Dichte, nicht allein.")

# ---------------------------------------------------------------- 3
P()
P("3) Sonnenkreis: Γ/Ω")
M = 6e10 * M_sun
r = 8 * kpc
aN = G * M / r**2
Gam = verlust_gamma(aN)
v_rule = v_kreis_punkt(M, r)
Om_rule = v_rule / r
Om_220 = 220e3 / r
s_sun = float(s_(aN))
P(f"   a_N = {aN:.4e} m/s² ({aN/A0:.3f} a₀), s = {s_sun:.4f}, Γ = a_N/c = {Gam:.4e} s⁻¹, Γτ_D = s² = {Gam*TD:.4f}")
P(f"   Weg A (Kern): Ω Regel = {Om_rule:.4e} s⁻¹ (v = {v_rule/1e3:.1f} km/s):  Γ/Ω = {Gam/Om_rule:.3e};  mit Ω = 220 km/s/8 kpc: {Gam/Om_220:.3e}")
vN2 = GM_sun * 6e10 / r
P(f"   Weg B (Hand): Γ/Ω = v_N²/(c·v),  v_N² = GM/r = {vN2:.4e} m²/s²,  c·v = {c*v_rule:.4e}  →  {vN2/(c*v_rule):.3e}")
P(f"   Tiefbereich allgemein: Γ/Ω = (v_N/v)²·v/c;  v/c = {v_rule/c:.3e}")
P(f"   a_N, die 5,2·10⁻⁴ gäbe: {5.2e-4*Om_rule*c:.3e} m/s² (≈ {5.2e-4*Om_rule*c/A0:.2f} a₀, z. B. anderes Baryonenmodell)")
P("   Schranke 3,07·10⁻⁴ Ω: eigene Rechnung kern_gedaechtnis.py für den reinen Exponentialkern,")
P("   Hub der Radialbewegung nach 10 Mrd. Jahren ≥ 90 %; die Annahme 'alte Scheiben behalten ihre Radialstreuung'")
P("   und Ω = 220 km/s/8 kpc sind Messgrößen der Literatur (nicht aus eigenen Daten).")

# ---------------------------------------------------------------- 4
P()
P("4) Gewicht bei Bahnfrequenz, R = |K̂(iΩ)|/K̂(0) (Maxwell: λ/√(λ²+Ω²))")
for name, Om in (("Ω Regel", Om_rule), ("Ω 220 km/s", Om_220)):
    zA = cmath.sqrt((Gam + 1j * Om) * TD)
    RA = math.exp(-(zA.real - s_sun))
    reB = math.sqrt(TD) * math.sqrt((math.hypot(Gam, Om) + Gam) / 2)     # Re√(a+ib) = √((|z|+a)/2)
    RB = math.exp(-(reB - s_sun))
    P(f"   {name}: Ωτ_D = {Om*TD:.1f};  Weg A (komplexe Wurzel) Re = {zA.real:.4f}, R = {RA:.3e};  "
      f"Weg B (√((|z|+a)/2)) Re = {reB:.4f}, R = {RB:.3e}")
    P(f"      Maxwell mit λ_K = Γ: R = {Gam/math.hypot(Gam, Om):.3e};  Schranke 3,07·10⁻⁴")
P("   Hand: Re ≈ √(Ωτ_D/2) = " + f"{math.sqrt(Om_rule*TD/2):.3f}" + ", minus s → e^{−33,8} ≈ 2·10⁻¹⁵.")

# ---------------------------------------------------------------- 5
P()
P("5) Anteil des Kerns bis zur Zeit t:  F(t)/q,  F(t) = ½[e^{−s}erfc(√(τ_D/4t) − √(Γt)) + e^{s}erfc(√(τ_D/4t) + √(Γt))]")
def F_closed(t, s, td):
    Gm = s * s / td
    a, b = math.sqrt(td / (4 * t)), math.sqrt(Gm * t)
    return 0.5 * (math.exp(-s) * erfc(a - b) + math.exp(s) * erfc(a + b))
for s in (0.3, 1.0, s_sun):
    for tname, t in (("13,8 Mrd. J.", 13.8 * Gyr), ("10 Mrd. J.", 10 * Gyr), ("Umlauf 2π/Ω", 2 * math.pi / Om_rule)):
        Gm = s * s / TD
        num = quad(lambda T: f_levy(T, TD) * math.exp(-Gm * T), 0, t, limit=400, points=[t / 10])[0]
        cl = F_closed(t, s, TD)
        P(f"   s = {s:.3f}, {tname:12s}: numerisch {num/math.exp(-s):.4e}, geschlossen {cl/math.exp(-s):.4e}  (Anteil an q)")

P()
P("Ergebnis: K(T) = √(τ_D/4π)T^{−3/2}e^{−τ_D/4T − ΓT}, K̂(p) = e^{−√((p+Γ)τ_D)}. λ_K = Γ ist der Dämpfungsfaktor,")
P("die Lévy-Dichte legt das Gewicht aber auf T ~ τ_D (Modus τ_D/6 ≈ 14 Mrd. J.); bei Bahnfrequenz ist der Kern um ~10⁻¹⁵")
P("unterdrückt. In Echtzeit gelesen hätte der Kern im Weltalter nur einen Bruchteil von q aufgebaut: T ist die Rückkehrzeit in")
P("der inneren Richtung (ϖ = c·imaginäre Zeit), der Kern wirkt stationär auf q, nicht auf ṙ; die Exponential-Schranke gilt nicht.")

with open("kopfgleichung_kern.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
