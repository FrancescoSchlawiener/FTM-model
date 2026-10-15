# -*- coding: utf-8 -*-
"""
Echte Kontraktion oder nur fallende Atemrate? Minimalmodell des Ganzen (G).
Blick: L(t) bzw. a(t) ist das Atmen des Ganzen (G), H = ȧ/a die Rate je Epoche; Λ ist die Schließlänge in der
inneren Richtung (I), fest oder mit neuer Faltung an der Engstelle wachsend; der Horizont c/H und seine Entropie
gehören je einem Beobachter. Nichts davon ist eine Dichte am Ort.

Ausdehnungsgesetz auf dem flachen T³ (nur die Volumen-Nullmode atmet, keine Krümmung):
  E² = H²/H0² = Ω_m a⁻³ (Bremsglied der Materie, σ = ΣΔs − ΔL = 0) + Ω_r a⁻⁴ + Ω_s(Λ) + Ω_x(a)
  Ω_s = (π²/45)Θ² = Ω_s0·(Λ0/Λ)²  (Vakuumteil, Rückstellung des Ganzen),  Ω_x = übrige zwei Drittel (offen)
  Ḣ = −(3/2)H² Σ(1 + w_i)Ω_i/E²,   q = −1 − Ḣ/H²
Umkehr des Atmens (echte Kontraktion) ⇔ E² = 0 bei endlichem a.
1 Umkehrbedingung.  2 Λ fest oder wachsend: Endrate, Kac-Maß.  3 Vorzeichen von Ḣ mit Λ(a) aus dem Kern.
4 Kreislaufbedingung aus 04 bei fester Dunkler Energie.  5 von Hand.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, c, G, l_P, Mpc, Gpc, yr, Gyr, OMEGA_M, OMEGA_R, OMEGA_DE, H0_KMS, H0_si,
                      schliesslaenge, omega_vakuum, lambda_vakuum, horizontentropie, kreislauf_mindestrate)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def H_(t):
    P(); P("=" * 100); P(t); P("=" * 100)

H0 = H0_si()
Lam0 = schliesslaenge(A0)
Om, Or = OMEGA_M, OMEGA_R
Os = omega_vakuum(A0)
Ox = 1 - Om - Or - Os
SDOT = 7.1e86                       # Ṡ_i je Horizont heute, k/s (Wert des Modells, 03)
KMS = Mpc / 1e3

H_("1 Umkehrbedingung auf dem flachen T³")
P(f"Ω_m = {Om:.4f}, Ω_r = {Or:.2e}, Ω_s = (π²/45)Θ² = {Os:.4f} (Λ = {Lam0/Gpc:.2f} Gpc), übrige Ω_x = {Ox:.4f}")
P("Jeder Summand ist ≥ 0: Ω_s ∝ Γ(4)Σ1/n⁴ > 0 (Zählung der Rückläufe), Materie und Strahlung ≥ 0, T³ ohne Krümmung.")
P("⇒ E² > 0 für alle a: H wechselt das Vorzeichen nie. Numerisch, min E² über a ∈ [1, 10⁶]:")
aa = np.geomspace(1, 1e6, 4001)
for name, Ex in [("Ω_x fest (w = −1)", lambda a: Ox + 0 * a), ("Ω_x ∝ a⁻³·⁽¹⁺ʷ⁾, w = −0,9", lambda a: Ox * a**-0.3),
                 ("Ω_x ∝ a⁻², Bindung eingefroren", lambda a: Ox * a**-2.0), ("Ω_x = 0 ab heute", lambda a: 0 * a)]:
    E2 = Om * aa**-3 + Or * aa**-4 + Os + Ex(aa)
    P(f"  {name:32s} min E² = {E2.min():.4f}  → H∞ = {H0_KMS*math.sqrt(E2[-1]):.2f} km/s/Mpc")
Hs = H0 * math.sqrt(Os)
P(f"Endrate bei festem Λ, Rest verklungen: H_s = H0·√Ω_s = {Hs*KMS:.2f} km/s/Mpc = π²c/(√45·Λ) = "
  f"{math.pi**2*c/(math.sqrt(45)*Lam0)*KMS:.2f}")
P()
P("Was eine Umkehr bräuchte: einen Summanden Ω_n < 0 mit |Ω_n(a)| = Ω_m a⁻³ + Ω_s bei endlichem a.")
P("Bei festem Λ bleibt Ω_s; ein negativer Summand ∝ a^{−3(1+w_n)} holt ihn nur ein, wenn er wächst (w_n < −1).")
for wn in [-0.5, -1.0, -1.1]:
    On = 0.01
    E2 = Om * aa**-3 + Os - On * aa**(-3 * (1 + wn))
    i = np.argmax(E2 <= 0) if np.any(E2 <= 0) else None
    P(f"  Ω_n = −{On}, w_n = {wn:+.1f}: " + (f"Umkehr bei a = {aa[i]:.1f}" if i else f"keine Umkehr (min E² = {E2.min():.3f})"))
P("Mit wachsendem Λ fällt Ω_s ∝ Λ⁻² gegen null, bleibt aber > 0: auch dann keine Umkehr, H → 0 von oben.")
P("Das Modell hat keinen negativen Summanden: die Rückstellung σ = 0 ist das Bremsglied (w = 0, ≥ 0), der Vakuumteil > 0.")

H_("2 Λ fest oder wachsend: Endzustand, Kac-Maß, Rate")
S_L = 45 / math.pi**3 * (Lam0 / l_P)**2
P(f"Horizontentropie zum Vakuumteil S_Λ = (45/π³)(Λ/ℓ_P)² = {S_L:.3e} k  (Probe πc²/(ℓ_P²H_s²) = {horizontentropie(Hs):.3e})")
P(f"Λ fest: Maß des Zustandsraums je Horizont endlich, e^{{S_Λ}}; Kac: Rückkehr des Zustands nach ~ e^{{S_Λ}} "
  f"(log10 ≈ {S_L/math.log(10):.2e}). Das ist ein Endzustand (H = H_s fest) mit Rückkehr als Schwankung.")
eps0, _, _ = kreislauf_mindestrate(Lam0, SDOT)
P(f"Λ wächst mindestens mit ε = Λ̇/Λ = Ṡ_i/(2S_Λ) = {eps0:.3e} s⁻¹ (1/ε = {1/eps0/yr:.2e} a) ⇒ dS_Λ/dt ≥ Ṡ_i:")
P("  das Maß wächst mindestens so schnell wie Entropie entsteht; Kac greift nicht, der Zustand kehrt nicht zurück (Spirale),")
P("  die Form kehrt zurück (flache Matrix → Wirbel → Faltung).")
K = math.pi**3 * l_P**2 * SDOT / 45           # dΛ²/dt bei festem Ṡ_i je Horizont
t2 = 3 * Lam0**2 / K
P(f"Mit Ṡ_i je Horizont fest: Λ² = Λ0² + (π³ℓ_P²Ṡ_i/45)·t ⇒ Λ ∝ t^½ spät, H_s ∝ 1/Λ; H_s halbiert sich nach {t2/yr:.2e} a")
P(f"  Rate: Ḣ_s/H_s = −ε,  q = −1 + ε/H_s = −1 + {eps0/Hs:.2e};  w_s = −1 + (2/3)ε/H = −1 + {2/3*eps0/H0:.2e} (heute)")
P("  ⇒ Λ fest: H → H_s = const (Endzustand, verletzt 04). Λ wächst: H → 0 asymptotisch, ohne Umkehr.")

H_("3 Vorzeichen von Ḣ mit Λ(a) = Λ_m(1 + Ω_r/(Ω_m a)) (Kern) und neuer Faltung ε")
P("Ḣ/H² = −(3/2)[Ω_m a⁻³ + (4/3)Ω_r a⁻⁴ + (1 + w_s)Ω_s]/E²,  1 + w_s = (2/3)·d ln Λ/d ln a,")
P("d ln Λ/d ln a = −(Ω_r/(Ω_m a))/(1 + Ω_r/(Ω_m a)) + ε/H   (Strahlung verliert Faltung; neue Faltung an der Engstelle)")
def Ea(a, Oxf):
    r = Or / (Om * a); Lr = (1 + r) / (1 + Or / Om)
    return math.sqrt(Om * a**-3 + Or * a**-4 + Os / Lr**2 + Oxf(a))
def Hdot_H2(a, Oxf, eps):
    r = Or / (Om * a); Lr = (1 + r) / (1 + Or / Om)
    E = Ea(a, Oxf); H = H0 * E
    dlnL = -r / (1 + r) + eps / H
    ws1 = 2 / 3 * dlnL
    return -1.5 * (Om * a**-3 + 4 / 3 * Or * a**-4 + ws1 * Os / Lr**2) / E**2
Ox0 = lambda a: 0.0
P("  (Rest Ω_x = 0, nur Ω_m, Ω_r, Ω_s(Λ(a)); bei a = 1 daher nicht der heutige Wert)")
for a in [1, 10, 50, 81, 100, 1e3, 1e6, 1e12, 1e14, 1e15, 1e16]:
    P(f"  a = {a:9.3g}:  Ḣ/H² = {Hdot_H2(a, Ox0, eps0):+.3e}   q = {-1 - Hdot_H2(a, Ox0, eps0):+.6f}")
a_x1 = brentq(lambda a: Hdot_H2(a, Ox0, eps0), 10, 1e4)
a_x2 = brentq(lambda a: Hdot_H2(a, Ox0, eps0), 1e10, 1e18)
tdS = lambda a1, a2: quad(lambda la: 1 / (H0 * Ea(math.exp(la), Ox0)), math.log(a1), math.log(a2), limit=400)[0]
P(f"  Ḣ = 0 bei a₁ = {a_x1:.1f} (Materie gegen Strahlungsverlust der Faltung, von Hand a² ≈ 1,5Ω_m²/(Ω_rΩ_s) = "
  f"{math.sqrt(1.5*Om**2/(Or*Os)):.1f}), t ≈ {tdS(1, a_x1)/Gyr:.0f} Gyr nach heute")
P(f"  Ḣ = 0 bei a₂ = {a_x2:.2e} (neue Faltung überwiegt, von Hand a ≈ Ω_rH_s/(Ω_mε) = {Or*Hs/(Om*eps0):.2e}), "
  f"t ≈ {tdS(1, a_x2)/yr:.2e} a  [bei Ṡ_i fest]")
P(f"  dazwischen Ḣ > 0, Betrag: H steigt insgesamt um den Faktor 1 + Ω_r/(Ω_m a₁) = 1 + {Or/(Om*a_x1):.1e}")
P("  Mit Rest Ω_x (w_x ≥ −1) kommt nur ein negativer Beitrag zu Ḣ hinzu. Heute (Rest w = −1):")
P(f"  Ḣ0/H0² = {Hdot_H2(1, lambda a: Ox, eps0):+.4f},  q0 = {-1 - Hdot_H2(1, lambda a: Ox, eps0):+.4f}; "
  f"minimal (Rest weg, nur Ω_m + Ω_s): Beschleunigung ab Ω_m a⁻³ = 2Ω_s, z = {(2*Os/Om)**(1/3)-1:.3f}")

H_("4 Kreislaufbedingung (04) −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k) bei fester Dunkler Energie (w = −1, Ω_DE = %.3f)" % OMEGA_DE)
def cond(a):
    E = math.sqrt(Om * a**-3 + Or * a**-4 + OMEGA_DE); H = H0 * E
    mHd = 1.5 * H0**2 * (Om * a**-3 + 4 / 3 * Or * a**-4)
    return mHd - l_P**2 * H**3 * SDOT / (2 * math.pi * c**2)
a_v = brentq(cond, 1, 1e9)
tv = quad(lambda la: 1 / (H0 * math.sqrt(Om * math.exp(-3*la) + OMEGA_DE)), -12, math.log(a_v), limit=400)[0]
P(f"  verletzt ab a = {a_v:.1f}, Weltalter t = {tv/yr:.2e} a  (04: ≈ 2,5·10¹¹ a)")
P("  Mit wachsendem Λ (ε ≥ Ṡ_i/(2S_Λ)) gilt −Ḣ = εH am Ende: die Bedingung ist genau die Mindestrate aus 03, erfüllt ohne Umkehr.")
P("  Mit Λ(a) aus dem Kern (Strahlung verliert Faltung) verlangt dieselbe Bedingung im Fenster a₁ < a < a₂ mehr:")
P("  ε ≥ Ṡ_i/(2S_Λ) + (Ω_r/(Ω_m a))·H  (sonst Ḣ > 0, die Horizontfläche schrumpft):")
for a in [a_x1, 1e3, 1e6, 1e10, 1e14, a_x2]:
    r = Or / (Om * a); need = eps0 + r / (1 + r) * Hs
    P(f"    a = {a:9.3g}: ε_nötig = {need:.2e} s⁻¹ = {need/eps0:.2e}·ε_min")

H_("5 Von Hand")
P("  Λ = c²/(2a₀) = 8,988·10¹⁶/2,2·10⁻¹⁰ = 4,085·10²⁶ m;  H_s = π²c/(√45Λ) = 9,870·2,998·10⁸/(6,708·4,085·10²⁶)")
P("      = 1,080·10⁻¹⁸ s⁻¹ = 33,3 km/s/Mpc;  H_s/H0 = 1,080/2,184 = 0,494 = √0,244 = √Ω_s.")
P(f"  Kern: {Hs*KMS:.2f} km/s/Mpc, √Ω_s = {math.sqrt(Os):.4f}")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kontraktion_atemrate.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
