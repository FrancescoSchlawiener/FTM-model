# -*- coding: utf-8 -*-
# Zählt der Vakuumteil einmal oder je geschlossener Richtung des T³ (Ω_DE = 3·Ω_s?)?
# A Q̃(s) aus n_B des Kerns (Q̃ fehlt im Kern); wo steckt die 3 der Raumdimension in π⁴/15?
# B Friedmann: die 3 in Ω_s = Λ_s c²/(3H²) ist die Spur über die drei Richtungen (G₀₀ = 3H²/c²)
# C Topologie: welcher Klasse gehört ein Vakuumterm an (b₀, b₁, b₃)? Rückstellung = Volumenmode
# D Alternative „je Schleife eine eigene Rückkehr“: Λ_i ∝ L_i ∝ a, oder Λ_i fest
# E Zahlen, falls Faktor 3: H₀ aus Flachheit, a₀ aus H₀, Q-Fenster, w
import sys, io, math
import numpy as np
from scipy.integrate import quad
from scipy.special import zeta, gamma
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ftm_kern import *

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []; P = lambda s="": out.append(s)
pi = math.pi
H100 = H0_si(100.0)
Gpc = 1e9 * pc

# ---------------------------------------------------------------- A
P("A) Q̃(s) = −4∫_s^∞ x³ n_B(x) dx (aus n_B des Kerns; Q̃ selbst steht nicht im Kern)")
def Qt(s):
    return -4 * quad(lambda x: x**3 * n_B(x), s, np.inf, limit=200)[0] if s > 0 else \
           -4 * quad(lambda x: x**3 * n_B(x), 1e-12, np.inf, limit=200)[0]
try:
    import mpmath as mp
    def Qt_li(s):
        z = mp.e**(-s)
        return float(-4 * (s**3 * mp.polylog(1, z) + 3 * s**2 * mp.polylog(2, z) + 6 * s * mp.polylog(3, z) + 6 * mp.polylog(4, z)))
    for s in (0.3, 1.0, 3.0):
        P(f"   s = {s}: Integral {Qt(s):.8f}, Polylog-Form (07) {Qt_li(s):.8f}")
except ImportError:
    P("   mpmath fehlt: Polylog-Form nicht verglichen")
P(f"   Q̃(0) = {Qt(0):.8f};  −4π⁴/15 = {-4*pi**4/15:.8f};  Γ(4)ζ(4) = {gamma(4)*zeta(4):.8f}")
P("   Bose-Leiter ∫₀^∞ s^d n_B ds = Γ(d+1)ζ(d+1) (Kern-n_B):")
for d in (1, 2, 3):
    I = quad(lambda x: x**d * n_B(x), 1e-12, np.inf, limit=200)[0]
    P(f"     d = {d}: {I:.6f}  Γ(d+1)ζ(d+1) = {gamma(d+1)*zeta(d+1):.6f}")
P("   Planck-Lesart von d = 3: ∫d³k |k| n_B(|k|) = 4π·∫k³n_B dk; die drei Richtungen sind in x² (Zustandsdichte) und 4π schon integriert.")
P("   Allgemein ∫d^dk |k| n_B = S_{d−1}·Γ(d+1)ζ(d+1):")
for d in (1, 2, 3):
    S = 2 * pi**(d / 2) / gamma(d / 2)
    P(f"     d = {d}: S_(d−1) = {S:.4f},  Γ(d+1)ζ(d+1) = {gamma(d+1)*zeta(d+1):.4f}")
P(f"   Je Richtung (d = 1) gezählt: 3·Γ(2)ζ(2) = {3*zeta(2):.4f} = π²/2, nicht 3·π⁴/15 = {3*pi**4/15:.4f}.")
P("   ⇒ 3·(π⁴/15) wäre die dreidimensionale Größe dreimal, also die Richtungen doppelt gezählt.")
P("   Herkunft der Potenz: s⁴ = ℓ₀²Y (Y quadratisch im Gradienten, s² = a_N/a₀) ⇒ dQ̃/ds = 4s³n_B; Gewicht 1/n⁴ = 1/n^(d+1), d = 3.")
P()

# ---------------------------------------------------------------- B
P("B) Friedmann: H² = (8πG/3)ρ + Λc²/3; G₀₀ = 3(ȧ/a)²/c² ist die Summe der drei Raumrichtungen.")
Lam = schliesslaenge(A0); Ls = lambda_vakuum(Lam); H = H0_si()
Om_s = Ls * c**2 / (3 * H**2)
P(f"   a₀ = {A0:.2e}: Λ = {Lam/Gpc:.3f} Gpc, Λ_s = {Ls:.4e} m⁻², Ω_s = Λ_s c²/(3H₀²) = {Om_s:.4f}  (Kern omega_vakuum: {omega_vakuum():.4f})")
P(f"   3Ω_s = Λ_s c²/H₀² = {Ls*c**2/H**2:.4f}: der Faktor 3 hebt genau die 3 aus G₀₀ auf (die Spur über die drei Richtungen).")
P("   Vakuumterm T_μν = −ρ_s g_μν: je Richtung Druck p = −ρ_s; die Richtungen stehen schon in ρ + 3p = −2ρ_s (ä/a = +(8πG/3)ρ_s).")
P("   Die Energie ρ_s ist die 00-Komponente, einmal; 3 Druckkomponenten als 3 Energien zu zählen ist die Spur statt der Dichte.")
P()

# ---------------------------------------------------------------- C
P("C) Topologie (Betti aus topologie_ganzes: T³ [1, 3, 3, 1], T³ × [0, Λ] [1, 3, 3, 1, 0]).")
P("   ℒ_𝒮 enthält (2/ℓ₀²)·Q̃(0) als Faktor der Volumenform √−g d⁴x: eine Zahl je Punkt, integriert über V = L₁L₂L₃;")
P("   die drei Schleifen gehen nur über das Volumen ein, der Term ist eine 0-Form (Skalar), keine 1-Form.")
P("   b₁ = 3 zählt harmonische 1-Formen (gleichförmige Vektorfelder, Verschiebungen) = Hintergrund der Triade, nicht die Rückstellung.")
P("   Rückstellung = gleichförmige Dehnung = Volumenmode (1 von 6 metrischen Nullmoden; 5 Scherung < 10⁻⁶ der Atemrate).")
P("   Die innere Richtung bringt keine Klasse (zusammenziehbar, H¹ = 0): Λ ist eine Länge, einmal; tiefste Mode des Ganzen k = 0 im Raum.")
P()

# ---------------------------------------------------------------- D
P("D) Alternative: je Schleife L_i eine eigene Rückkehr mit eigenem Schluss Λ_i.")
Lmin = KREISPAAR_GRENZE * d_LSS()
P(f"   (a) Λ_i ∝ L_i ∝ a (Schluss an der Raumschleife): Vakuumteil wirkt wie w_s = −⅓; heute 3·Ω_s = {3*Om_s:.3f} als Krümmung (gemessen Ω_K = {OMEGA_K_MESS[0]} ± {OMEGA_K_MESS[1]}) – ausgeschlossen.")
P(f"   (b) Λ_i = Λ fest, drei gleiche Kopien: verlangt drei innere Richtungen. Die Rückkehr ist aber weder Querrichtung noch Wicklung")
P("       des Torus (07: irrationale Drehung gibt höchstens drei Rückkehrzeiten, keine Wurzel; SPARC d_s = 1,03 ± 0,07 ⇒ eine eindimensionale Rückkehr).")
P("       Die innere Richtung ist orthogonal zum ganzen Raum, nicht zu einer Schleife: eine Richtung, ein Λ = ν_kosm·ΣΔs_b (nur das Volumen).")
P("   (c) je Falte eine Halbgerade: alle bis dasselbe Λ; der Vakuumterm hängt an keinem einzelnen Trichter (sonst ∝ Zahl der Falten).")
P(f"   Nebenbefund: c²/L_min = {c**2/Lmin:.3e} m/s² (L_min = {Lmin/Gpc:.2f} Gpc) – Λ = L/2 wäre zahlengleich, aber L atmet, Λ nicht.")
P()

# ---------------------------------------------------------------- E
P("E) Zahlen, falls Ω_DE = 3Ω_s (nur zur Prüfung).")
wm = OMEGA_B_H2 * (1 + Q_KOSM); wr = OMEGA_R * h_**2
P(f"   Ω_DE (Kern, Planck-Flachheit) = {OMEGA_DE:.4f};  ω_m = (1+Q)ω_b = {wm:.5f}")
for a0 in (1.00e-10, 1.06e-10, 1.10e-10, 1.13e-10, 1.20e-10):
    P(f"   a₀ = {a0:.2e}: Θ = {theta_(a0, H):.4f}, Ω_s = {omega_vakuum(a0):.4f}, 3Ω_s = {3*omega_vakuum(a0):.4f}  (H₀ = {H0_KMS})")
a0_req = math.sqrt(OMEGA_DE / (3 * K_S)) * c * H / (2 * pi)
P(f"   verlangt: 3Ω_s = Ω_DE ⇒ a₀ = {a0_req:.4e}, Θ = {theta_(a0_req, H):.4f}")
P("   Flachheit mit ω_m fest: H₀² = H₁₀₀²(ω_m + ω_r) + (4π⁴/15)(a₀/c)²  [1× Vakuumteil: 4π⁴/45]")
for a0 in (1.00e-10, 1.06e-10, 1.10e-10, 1.13e-10):
    h3 = math.sqrt(H100**2 * (wm + wr) + (4 * pi**4 / 15) * (a0 / c)**2) / H0_si(1.0)
    h1 = math.sqrt(H100**2 * (wm + wr) + (4 * pi**4 / 45) * (a0 / c)**2) / H0_si(1.0)
    P(f"   a₀ = {a0:.2e}: H₀(3Ω_s) = {h3:.2f} km/s/Mpc, Ω_m = {wm/(h3/100)**2:.4f};   H₀(nur Ω_s) = {h1:.2f}")
a0lo, a0hi = 1.06e-10 - math.hypot(0.09e-10, 0.04e-10), 1.06e-10 + math.hypot(0.09e-10, 0.04e-10)
hs = [math.sqrt(H100**2 * (wm + wr) + (4 * pi**4 / 15) * (a / c)**2) / H0_si(1.0) for a in (a0lo, 1.06e-10, a0hi)]
P(f"   SPARC a₀ = 1,06 ± 0,10: H₀ = {hs[1]:.1f} (+{hs[2]-hs[1]:.1f} / −{hs[1]-hs[0]:.1f}) km/s/Mpc – eine Gleichung, zwei Unbekannte: gibt H₀(a₀) oder a₀(H₀), nicht beide.")
P("   Die CMB-Winkel legen bei konstantem Ω_DE H₀ = 67,4 fest (ΛCDM); dann a₀ = %.4e. Unabhängige Vorhersage entsteht nur mit a₀ aus Galaxien (±9 %%)." % a0_req)
qlo, qhi = q_fenster()
P(f"   Q-Fenster (a₀ = 1,1): Q ∈ [{qlo:.2f}, {qhi:.2f}]; mit a₀ = {a0_req:.3e}: Q ≤ {q_fenster(a0_req)[1]:.2f}; Q = {Q_KOSM:.3f} liegt drin.")
P(f"   (1+Q)a₀ = {(1+Q_KOSM)*a0_req:.3e} ≤ 7,17·10⁻¹⁰: erfüllt. Faktor 3 berührt Λ, Q, Zelle nicht (Λ = ν ΣΔs_b ist unabhängig von der Zählung).")
P(f"   a₀ fest seit z ≈ 10 ⇒ w = −1 exakt, Hintergrund = ΛCDM mit Ω_DE; Ω_DE(z = 1090)/Ω_m-Verhältnis wie ΛCDM. Gegen w₀ > −1 [G] kein Raum für die Bindungsrate.")
P()

# ---------------------------------------------------------------- Hand
P("Handprobe: a₀ = 1,06e-10, H₀ = 67,4 ⇒ H = 2,1843e-18 s⁻¹, cH = 6,5484e-10, 2πa₀ = 6,6602e-10, Θ = 1,0171, Θ² = 1,0345,")
P("           Ω_s = 0,21932·1,0345 = 0,2269, 3Ω_s = 0,6807.  Code: %.4f" % (3 * omega_vakuum(1.06e-10)))
P("Handprobe H₀: ω_m H₁₀₀² = 0,14237·(3,2408e-18)² = 1,4953e-36; (4π⁴/15)(1,1e-10/c)² = 25,976·1,3463e-37 = 3,4971e-36;")
P("           Summe 4,992e-36 ⇒ H₀ = 2,2344e-18 s⁻¹ = 68,95 km/s/Mpc.")

txt = "\n".join(out)
print(txt)
open(__file__.replace('.py', '.txt'), 'w', encoding='utf-8').write(txt + "\n")
