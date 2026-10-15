"""Schluss oben: Das Ganze (Zelle) als oberste Ebene, deren Rückkehr sich schließt.

Fixpunktbedingungen am oberen Ende der Leiter, VOR dem Rechnen festgelegt (aus dem Modell begründet):
  Bezeichnungen: Δs_b = ΣΔs_b(Zelle) = GM_b/(3c²) (M_b mitbewegt fest, G), Λ (I, seit z ≈ 10 fest), a₀ = c²/2Λ,
  r*² = GM_b/a₀ = 6Δs_b·Λ = r₀Λ (Kac-Beziehung, jede Masse), r₀ = 6Δs_b, ℓ₀ = 2Λ.
  Jede Bedingung der Form Λ = m·Δs_b gibt mit der Lesart Λ = ν_kosm·Δs_b: ν = m, Q = m − 1.

  F1  Leiter in ϖ (07/08: r₀ : r* : Λ = (2δ*)² : 2δ* : 1, Halbgerade, Schluss reflektierend bei Λ):
      oberstes Ende m = r*/r₀ = 1, d. h. r₀ = r* = Λ.  Wegen r*² = r₀Λ zieht jede der drei Gleichheiten die anderen nach.
      Diffusionsform: ϖ_R² = D⟨T⟩ = Λ² (die Rückkehr erreicht in einer Kac-Zeit genau den Schluss) – dieselbe Bedingung.
  F1s dieselbe als Schranke: ϖ_R(Zelle) ≤ Λ (Rückkehrtiefe nicht jenseits der Schließung; sonst Form ungültig).
  F2  Triade in Raumlängen (01/08, M_U): 3Δs = r* = ℓ₀ = 2Λ (auch hier zieht eine Gleichheit die andere nach).
  F3  Kreis 2Λ statt Halbgerade (ℓ₀ = Umfang):  (a) r₀ = ℓ₀,  (b) r* = ℓ₀ (= F2),  (c) r₀ = r* (= F1).
  F4  Kac am Ort der Zelle: 1/ν = μ(Ort)/μ(X), X = [0, Λ] bzw. Kreis 2Λ, Ort = r₀, 3Δs_b oder Δs_b.
  F5  Rückkehr ohne Gefälle über die tiefste Mode des Ganzen: räumlich k = 0 (Γ = 0, q = 1); in I k₁ = π/2Λ
      (Ort Dirichlet, Λ reflektierend) bzw. π/Λ (Kreis 2Λ); Rückkehrtiefe ϖ_R = r*(M_b); q = e^{−kϖ_R} bzw. geschlossen
      cosh(k(Λ − ϖ_R))/cosh(kΛ); Fixpunkt ν = 1/(1 − q) mit r*/Λ = √(6/ν).
  F2', F1' mit gezählter Faltung Δs_eff = ν·Δs_b = Λ statt Δs_b.
a₀: Jede Bedingung verknüpft nur M_b (G, mitbewegt) mit Λ (I) und ist daher zeitfest; L tritt nur über M_b = ρ_b L³ auf.
Eine Gleichung, zwei Unbekannte (Λ, L): a₀ folgt erst mit L; Bedingungen mit L_phys oder H wären nicht zeitfest.
"""
import math
import numpy as np
from scipy.optimize import brentq
import ftm_kern as K

c, G, Gpc = K.c, K.G, K.Mpc * 1e3
Q_P, dQ_P = K.Q_KOSM, 0.0646            # Planck 2018 TT,TE,EE+lowE+lensing [G], wie herkunft_q
Q_B, dQ_B = 5.3225, 0.0525              # Planck 2018 + BAO [G], wie herkunft_q
A0_BAND = (1.04e-10, 1.10e-10, 1.19e-10)

H0 = K.H0_si()
rho_b = K.OMEGA_B * 3 * H0**2 / (8 * math.pi * G)
dLSS = K.d_LSS()
Lmin = K.KREISPAAR_GRENZE * dLSS


def ds_b(L):
    return G * rho_b * L**3 / (3 * c**2)


def a0_bei(m, L):
    return c**2 / (2 * m * ds_b(L))


def L_fuer(m, a0):
    return (3 * c**4 / (2 * m * G * rho_b * a0)) ** (1 / 3)


def sig(Q):
    return (Q - Q_P) / dQ_P, (Q - Q_B) / dQ_B


def linie(ch="="):
    print(ch * 112)


linie()
print("0  Grundgrößen")
linie()
print(f"   ρ_b = {rho_b:.4e} kg/m³ (Ω_b = {K.OMEGA_B:.5f});  d_LSS = {dLSS/Gpc:.3f} Gpc;  L_min = 0,985·d_LSS = {Lmin/Gpc:.3f} Gpc")
print(f"   Δs_b(d_LSS) = {ds_b(dLSS)/Gpc:.4f} Gpc;  Δs_b(L_min) = {ds_b(Lmin)/Gpc:.4f} Gpc")
print(f"   gemessen Q = {Q_P:.4f} ± {dQ_P} (Planck) bzw. {Q_B} ± {dQ_B} (+BAO) [G]; a₀-Band 1,04–1,19·10⁻¹⁰")
nu_m = 1 + Q_P
print(f"   heute mit der Lesart (ν = {nu_m:.4f}): r*/Λ = √(6/ν) = {math.sqrt(6/nu_m):.4f},  r₀/Λ = 6/ν = {6/nu_m:.4f},  3Δs_b/ℓ₀ = 1,5/ν = {1.5/nu_m:.4f}")

# ------------------------------------------------------------------ 1
linie()
print("1  Bedingungen Λ = m·Δs_b  ⇒  ν = m, Q = m − 1; a₀ bei würfelförmiger Zelle")
linie()
kand = [
    ("F1  ϖ-Leiter: r₀ = r* = Λ (Halbgerade)", 6.0),
    ("F2  Raum-Triade: 3Δs = r* = ℓ₀ = 2Λ", 1.5),
    ("F3a Kreis 2Λ: r₀ = ℓ₀", 3.0),
    ("F3b Kreis 2Λ: r* = ℓ₀ (= F2)", 1.5),
    ("F3c Kreis 2Λ: r₀ = r* (= F1)", 6.0),
]
print(f"   {'Bedingung':42s} {'ν':>6s} {'Q':>6s} {'σ Planck':>9s} {'σ +BAO':>8s}  a₀(L_min)  a₀(d_LSS)   L für a₀ = 1,04 / 1,10 / 1,19 [Gpc]")
for name, m in kand:
    sP, sB = sig(m - 1)
    Ls = [L_fuer(m, a) / Gpc for a in A0_BAND]
    print(f"   {name:42s} {m:6.2f} {m-1:6.2f} {sP:+9.1f} {sB:+8.1f}  {a0_bei(m, Lmin)*1e10:8.4f}   {a0_bei(m, dLSS)*1e10:8.4f}   "
          + " / ".join(f"{x:6.2f}" for x in Ls) + f"   (L_min = {Lmin/Gpc:.2f})")
print("   Kreis 2Λ: r₀ = ℓ₀ gibt r*² = r₀Λ = ℓ₀²/2, r* = ℓ₀/√2 ≠ ℓ₀: auf dem Kreis fallen die drei Längen nie zusammen.")
print("   Halbgerade (F1) und Raum-Triade (F2) haben je genau einen Tripelpunkt (r*² = r₀Λ bzw. r*² = 3Δs·ℓ₀).")

# ------------------------------------------------------------------ 2
linie()
print("2  Gezählte Faltung statt Baryonen (Δs_eff = ν·Δs_b = Λ nach der Lesart)")
linie()
print("   F1': 6Δs_eff = Λ  ⇒ 6Λ = Λ: Widerspruch.   F2': 3Δs_eff = 2Λ ⇒ 3Λ = 2Λ: Widerspruch.")
print(f"   r*_eff/Λ = √(6ν·Δs_b/Λ) = √6 = {math.sqrt(6):.4f} > 1 für jedes ν: Die gezählte Masse liegt stets jenseits des Schlusses;")
print("   die Lesart setzt die Zelle mit gezählter Faltung auf ΣΔs_eff = Λ, d. h. M = 1,5·M_U (keine freie Bedingung).")
print("   Nur die Baryonenfaltung (wie in der Regel r* = √(GM_b/a₀)) gibt eine Bedingung.")

# ------------------------------------------------------------------ 3
linie()
print("3  F4  Kac am Ort der Zelle: 1/ν = μ(Ort)/μ(X) mit Λ = ν·Δs_b")
linie()
for X, fX in (("Halbgerade Λ", 1.0), ("Kreis 2Λ", 2.0)):
    for Ort, fO in (("r₀ = 6Δs_b", 6.0), ("3Δs_b", 3.0), ("Δs_b", 1.0)):
        # 1/ν = fO·Δs_b/(fX·ν·Δs_b) ⇔ 1 = fO/fX
        print(f"   X = {X:12s} Ort = {Ort:11s}: 1/ν = {fO/fX:.1f}/ν  ⇒ " + ("gilt für jedes ν (Umschreibung)" if abs(fO/fX-1) < 1e-12 else f"Widerspruch ({fO/fX:.1f} = 1)"))
print("   Kac legt ν nie fest: entweder Umschreibung der Lesart oder Widerspruch.")

# ------------------------------------------------------------------ 4
linie()
print("4  F5  Fixpunkt der eigenen Rückkehr der Zelle mit der tiefsten Mode der inneren Richtung")
linie()
print("   räumlich: tiefste Mode k = 0 ⇒ s = 0, q = 1, ν = ∞ (kein Fixpunkt).")


def q_offen(nu, kfak):
    s = kfak * math.sqrt(6 / nu)          # s = k·r*, k = kfak/Λ
    return math.exp(-s)


def q_zu(nu, kfak):
    x = math.sqrt(6 / nu)                 # r*/Λ
    if x > 1:
        return float("nan")
    return math.cosh(kfak * (1 - x)) / math.cosh(kfak)


for name, kfak in (("Halbgerade k₁ = π/2Λ", math.pi / 2), ("Kreis 2Λ k₁ = π/Λ", math.pi), ("k = 1/Λ (Λ als Radius: L/I verwechselt)", 1.0)):
    for form, qf in (("offen e^{−kϖ_R}", q_offen), ("geschlossen cosh", q_zu)):
        f = lambda nu: nu - 1 / (1 - qf(nu, kfak))
        grid = np.linspace(1.0001, 400, 400000)
        vals = np.array([f(v) for v in grid])
        ok = np.isfinite(vals)
        wz = np.where(np.sign(vals[ok][:-1]) != np.sign(vals[ok][1:]))[0]
        loes = [brentq(f, grid[ok][i], grid[ok][i+1]) for i in wz]
        txt = ", ".join(f"ν = {v:.4f} (Q = {v-1:.3f}, {sig(v-1)[0]:+.1f}σ; r*/Λ = {math.sqrt(6/v):.3f}{' > 1: jenseits des Schlusses, ungültig' if v < 6 else ''})" for v in loes) or "keine Lösung mit r* ≤ Λ"
        print(f"   {name:40s} {form:18s}: {txt}")

# ------------------------------------------------------------------ 5
linie()
print("5  Schranke statt Gleichung: ϖ_R(Zelle) = r*(M_b) ≤ Λ  ⇔  ν ≥ 6, Q ≥ 5 (untere Schranke, zeitfest)")
linie()
print("   Obere Schranke aus L ≥ L_min (keine Kreispaare): ν ≤ 3c⁴/(2Gρ_bL_min³a₀) – hängt an a₀ und Zellform (nicht zwingend)")
g = {"würfelförmig": 1.0, "sechseckig": 1.0491, "raumzentriert": 1.0911, "flächenzentriert": 1.1225}
for form, gg in g.items():
    V = (Lmin / gg) ** 3
    zeile = []
    for a in A0_BAND:
        nmax = 3 * c**4 / (2 * G * rho_b * V * a)
        zeile.append(f"{nmax:6.3f}")
    print(f"   {form:17s} ν_max bei a₀ = 1,04 / 1,10 / 1,19: " + " / ".join(zeile))
print(f"   Würfel, a₀ = 1,1: Fenster 6 ≤ ν ≤ {3*c**4/(2*G*rho_b*Lmin**3*1.1e-10):.3f}, also 5 ≤ Q ≤ {3*c**4/(2*G*rho_b*Lmin**3*1.1e-10)-1:.3f}; gemessen 5,364")
print(f"   Würfel: Fenster leer (ν_max < 6) für a₀ > {3*c**4/(2*G*rho_b*Lmin**3*6)*1e10:.4f}·10⁻¹⁰ (= a₀ bei F1 und L = L_min)")
print(f"   Abstand zur Gleichheit heute: r*/Λ = {math.sqrt(6/nu_m):.4f}, 1 − r*/Λ = {1-math.sqrt(6/nu_m):.4f}; Q − 5 = {Q_P-5:.3f} = {(Q_P-5)/dQ_P:.1f}σ")

# ------------------------------------------------------------------ 6
linie()
print("6  Zeitfestigkeit (Beobachterfehler): r*/Λ mit Λ(z)/Λ(0) = 1/a0_frueh_relativ(z), M_b mitbewegt fest")
linie()
for z in (0, 1, 10, 100, 1090, 1e4):
    fl = K.a0_frueh_relativ(0) / K.a0_frueh_relativ(z)
    print(f"   z = {z:7g}: Λ(z)/Λ(0) = {fl:.4f},  r*/Λ = {math.sqrt(6/nu_m/fl):.4f}  (F1 als Gleichung gälte nur bei einem Λ; die Schranke gilt in jeder Epoche)")
print("   Nicht zeitfest und daher ausgeschlossen: Bedingungen mit L_phys = L̄·a (G atmet), H, Θ oder d_LSS des heutigen Beobachters.")

# ------------------------------------------------------------------ 7
linie()
print("7  Handwert")
linie()
m = 6.0
Lam = m * ds_b(dLSS)
print(f"   F1, L = d_LSS: Δs_b = {ds_b(dLSS)/Gpc:.4f} Gpc, Λ = 6Δs_b = {Lam/Gpc:.3f} Gpc, a₀ = c²/2Λ = {c**2/(2*Lam):.4e} m/s²")
