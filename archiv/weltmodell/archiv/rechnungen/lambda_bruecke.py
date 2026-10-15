# -*- coding: utf-8 -*-
"""
Brücke Λ = ν_kosm·ΣΔs_b(Zelle): Gibt es einen Weg, Λ ohne ν festzulegen?

Werkzeug: Das Gewicht ist dieselbe Rückkehr wie die Regel. Für die ganze Zelle als Bereich im Atmen
(θ̃ = 3H, Hintergrund) gilt w_Zelle = 1 − e^{−√x}, x = 3H·τ_D = 6HΛ/c.

Kandidaten, VORAB festgelegt (keine nachgeschobenen):
  K1  Gewicht = Kac-Maß („der mitatmende Anteil trägt die nicht in Zentren gefaltete Länge“):
      Der Anteil der inneren Länge, der in den Zentren gefaltet liegt (Kac-Maß des offenen Bereichs, ΣΔs_b/Λ),
      ist der Anteil der Rückkehr der Zelle, der im Atmen ankommt:   ΣΔs_b/Λ = 1 − w_Zelle = e^{−√(6HΛ/c)}.
  K2  Leiterende an der Zelle („Zentrum, Wirbel und Schwelle auf einer Länge“, 01/08):
      Die Zelle als eine Falte (baryonisch, Schritt r₀ = 6ΣΔs_b) erreicht mit ihrer Rückkehrtiefe die Schließung:
      (a) ϖ_R = √(r₀Λ) = Λ (Halbgerade [0, Λ])  ⇒  Λ = 6·ΣΔs_b
      (b) r* = ℓ₀ = 2Λ (obere Schnittmasse M_U, Kreis 2Λ)  ⇒  Λ = 1,5·ΣΔs_b
  K3  Die Rückstellung schließt die Rückkehr:
      relative Dehnung der Raumschleife durch die Rückstellung (G, ΔL_Raum/L mit |ΔL_Raum| = Λ)
      = Anteil der Rückkehr der Zelle, der im Atmen schließt (I):   Λ/L = 1 − w_Zelle = e^{−√(6HΛ/c)}.
Jeder Kandidat ist ν-frei. K1 und K2 brauchen ΣΔs_b(L), K3 braucht L: Die Zelle L legt das Modell nicht fest
(nur L ≥ 0,985·d_LSS ohne Kreispaare). Verglichen wird ohne Anpassung mit a₀ = 1,1·10⁻¹⁰ (Kern) und mit
ν = 1 + Q_KOSM (eingesetzte Zahl, keine Vorhersage), bei L = L_min (Kreispaare), L = d_LSS und L aus a₀.
Nebenform: 1 − w der Zelle mit Abschluss bei Λ (cosh-Form), weil die Zelle als eine Falte ϖ ≈ Λ hat.

Blickrichtung: Λ, ϖ_R, τ_D [I]; ΣΔs_b [L, über die Zelle summiert]; ΔL_Raum, L [G, Raumschleife]; H [G, Epoche].
Rechnet nur mit ftm_kern. Ausgabe: lambda_bruecke.txt
"""
import io, os, sys, math
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

out = io.StringIO()


def p(*a):
    print(*a, file=out)


def kopf(t):
    p("=" * 110)
    p(t)
    p("=" * 110)


H0 = K.H0_si()
Gpc = K.Gpc
c = K.c
A0 = K.A0
NU = 1 + K.Q_KOSM                      # eingesetzt
dLSS = K.d_LSS()
Lmin = K.KREISPAAR_GRENZE * dLSS
L11 = K.zelle_aus_a0(A0)               # Zelle der Lesart bei a₀ = 1,1 (benutzt ν)
Lam11 = K.schliesslaenge(A0)
rho_b0 = K.OMEGA_B * 3 * H0**2 / (8 * math.pi * K.G)   # ν-frei (Ω_b direkt)
FORMEN = [("würfelförmig", 1.0), ("sechseckig", math.sqrt(3) / 2), ("raumzentriert", 4 / (3 * math.sqrt(3))),
          ("flächenzentriert", 1 / math.sqrt(2))]


def sds_b(L, f=1.0):
    """ΣΔs_b der Zelle (mitbewegt, heute): G·ρ_b·f·L³/(3c²)."""
    return K.G * rho_b0 * f * L**3 / (3 * c**2)


def eins_minus_w(H, Lam):
    """1 − w_Zelle, offene Form (Falte ≪ Λ): e^{−√(3H·τ_D)}, τ_D = 2Λ/c."""
    return math.exp(-math.sqrt(6 * H * Lam / c))


def eins_minus_w_cosh(H, Lam, sds):
    """1 − w der Zelle als eine Falte mit Abschluss bei Λ: cosh(k(Λ−ϖ))/cosh(kΛ), ϖ = √(6ΣΔs·Λ)."""
    x = 6 * H * Lam / c
    varpi = math.sqrt(6 * sds * Lam)
    kv = math.sqrt(x)
    kL = kv * Lam / varpi
    wert = (math.exp(-kv) + math.exp(kv - 2 * kL)) / (1 + math.exp(-2 * kL))   # = cosh(kL − kv)/cosh(kL), stabil
    return wert, varpi / Lam


kopf("0  Grundgrößen (heute)")
p(f"   a₀ = {A0:.2e} m/s² ⇒ Λ = c²/2a₀ = {Lam11/Gpc:.4f} Gpc [I];  H₀ = {H0:.5e} s⁻¹ [G];  x = 6H₀Λ/c = {6*H0*Lam11/c:.4f}")
p(f"   1 − w_Zelle (offen) = {eins_minus_w(H0, Lam11):.5f};  ν = 1 + Q_KOSM = {NU:.4f} (eingesetzt), 1/ν = {1/NU:.5f}")
p(f"   d_LSS = {dLSS/Gpc:.3f} Gpc;  L_min = 0,985·d_LSS = {Lmin/Gpc:.3f} Gpc;  L(a₀ = 1,1, Lesart) = {L11/Gpc:.3f} Gpc [G]")
p(f"   ΣΔs_b(L) = Gρ_bL³/(3c²) mit Ω_b = {K.OMEGA_B:.5f} (ν-frei):  L_min {sds_b(Lmin)/Gpc:.4f}, d_LSS {sds_b(dLSS)/Gpc:.4f}, "
  f"L(1,1) {sds_b(L11)/Gpc:.4f} Gpc")
p(f"   Probe Lesart: ν·ΣΔs_b(L(1,1)) = {NU*sds_b(L11)/Gpc:.4f} Gpc = Λ {Lam11/Gpc:.4f} Gpc")
p(f"   Probe Kern: a0_aus_zelle(L_min) = {K.a0_aus_zelle(Lmin):.4e} (Lesart, würfelförmig)")

# ------------------------------------------------------------------ K1
kopf("1  K1  ΣΔs_b/Λ = 1 − w_Zelle = e^{−√(6HΛ/c)}   (Gewicht = Kac-Maß)")
p("   Mit u = Λ/ΣΔs_b und b = 6H·ΣΔs_b/c:  ln u = √(b·u).  ln u/√u hat das Maximum 2/e bei u = e²;")
p("   Lösung nur für b ≤ 4/e² = 0,5413 (Berührung bei u = e² = 7,389).")
b_krit = 4 / math.e**2
p(f"   b_krit = {b_krit:.5f}")
p("   L               ΣΔs_b [Gpc]   b          Lösungen u = Λ/ΣΔs_b          Λ [Gpc]        a₀ [10⁻¹⁰]")
for name, L in [("L_min", Lmin), ("d_LSS", dLSS), ("L(1,1)", L11)]:
    S = sds_b(L)
    b = 6 * H0 * S / c
    if b > b_krit:
        p(f"   {name:7s} {L/Gpc:6.2f}   {S/Gpc:9.4f}   {b:8.4f}   keine (b > b_krit)")
    else:
        g = lambda u: math.log(u) - math.sqrt(b * u)
        u1 = brentq(g, 1.0 + 1e-12, math.e**2); u2 = brentq(g, math.e**2, 1e8)
        p(f"   {name:7s} {L/Gpc:6.2f}   {S/Gpc:9.4f}   {b:8.4f}   {u1:.4f}, {u2:.4f}   {u1*S/Gpc:.3f}, {u2*S/Gpc:.3f}")
S_krit = b_krit * c / (6 * H0)
L_krit = (3 * c**2 * S_krit / (K.G * rho_b0)) ** (1 / 3)
p(f"   Grenze: ΣΔs_b ≤ {S_krit/Gpc:.4f} Gpc ⇔ L ≤ {L_krit/Gpc:.3f} Gpc (würfelförmig) = {L_krit/Lmin:.3f}·L_min: durch Kreispaare ausgeschlossen")
Lam_b = math.e**2 * S_krit
p(f"   an der Grenze: Λ = e²·ΣΔs_b = {Lam_b/Gpc:.4f} Gpc, a₀ = {K.a0_aus_schliesslaenge(Lam_b):.4e}")
# nötiges ΣΔs_b für a₀ = 1,1
S_noetig = Lam11 * eins_minus_w(H0, Lam11)
L_noetig = (3 * c**2 * S_noetig / (K.G * rho_b0)) ** (1 / 3)
p(f"   verlangt für a₀ = 1,1: ΣΔs_b = Λ·e^(−√17,86) = {S_noetig/Gpc:.4f} Gpc ⇒ L = {L_noetig/Gpc:.3f} Gpc < L_min; ν_impl = 1/(1 − w) = {1/eins_minus_w(H0, Lam11):.2f} (gegen {NU:.3f})")
# mit der Lesart (ν eingesetzt): 1/ν = e^{−√x}
a0_K1L = 3 * H0 * c / math.log(NU)**2
p(f"   mit der Lesart Λ = ν·ΣΔs_b (ν eingesetzt): 1/ν = e^(−√(3Hτ_D)) ⇒ a₀ = 3H₀c/(ln ν)² = {a0_K1L:.4e} "
  f"({a0_K1L/A0:.2f}·a₀), L = {K.zelle_aus_a0(a0_K1L)/Gpc:.2f} Gpc < L_min")
# cosh-Form
p("   Nebenform mit Abschluss bei Λ (Zelle als eine Falte): gleiche Gleichung mit cosh(k(Λ−ϖ))/cosh(kΛ)")
for name, L in [("L_min", Lmin), ("d_LSS", dLSS)]:
    S = sds_b(L)
    us = np.geomspace(1.0001, 1e4, 20001)
    diff = np.array([S / (u * S) - eins_minus_w_cosh(H0, u * S, S)[0] for u in us])
    wechsel = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    p(f"   {name}: min(ΣΔs_b/Λ − (1 − w_cosh)) über u = 1…10⁴ = {diff.min():.4f};  Vorzeichenwechsel: {len(wechsel)}")
    for i in wechsel:
        u0 = brentq(lambda u: 1 / u - eins_minus_w_cosh(H0, u * S, S)[0], us[i], us[i + 1])
        vr = eins_minus_w_cosh(H0, u0 * S, S)[1]
        p(f"      Lösung u = {u0:.4f}: Λ = {u0*S/Gpc:.4f} Gpc, a₀ = {K.a0_aus_schliesslaenge(u0*S):.4e}, ϖ/Λ = √(6/u) = {vr:.4f}"
          f" {'> 1: Tiefe jenseits der Schließung, Form ungültig' if vr > 1 else ''}")
# Epoche
p("   Epoche (ΣΔs_b mitbewegt fest, H = H(z)): b(z) = b(0)·E(z) wächst ⇒ früher erst recht keine Lösung:")
for z in (0, 1, 10, 1090):
    p(f"     z = {z:5d}: b = {6*H0*K.E_hintergrund(z)*sds_b(dLSS)/c:.4e} (b_krit {b_krit:.4f}), d_LSS-Zelle")

# ------------------------------------------------------------------ K2
kopf("2  K2  Leiterende an der Zelle: (a) ϖ_R(Zelle) = Λ ⇒ Λ = 6ΣΔs_b;  (b) r* = 2Λ ⇒ Λ = 1,5ΣΔs_b")
p("   a₀ = c²/(2Λ) = c⁴/(2·m·G·ρ_b·V) mit m = 6 bzw. 1,5;  Lesart: m = ν.  Verhältnis zur Lesart bei gleichem L: ν/m")
for tag, m in [("(a)", 6.0), ("(b)", 1.5)]:
    p(f"   {tag} m = {m}:  ν_impl = {m} (Q = {m-1:.2f} gegen Q_KOSM = {K.Q_KOSM:.4f});  a₀/a₀(Lesart) = ν/m = {NU/m:.4f}")
    p("       Form                f         a₀(L_min) [10⁻¹⁰]   L(a₀ = 1,1) [Gpc]   L(1,1)/L_min")
    for name, f in FORMEN:
        a_min = c**2 / (2 * m * sds_b(Lmin, f))
        L_req = (3 * c**2 * (c**2 / (2 * A0 * m)) / (K.G * rho_b0 * f)) ** (1 / 3)
        p(f"       {name:18s} {f:.4f}    {a_min*1e10:8.4f}             {L_req/Gpc:8.3f}           {L_req/Lmin:.4f}")
    p(f"       Kac-Lemma μ(A)·⟨T_A⟩ ≤ μ(X): ν·ΣΔs_b ≤ Λ = {m}·ΣΔs_b verlangt ν ≤ {m}; mit ν = {NU:.4f}: "
      f"{'verletzt' if NU > m else 'erfüllt'} (Überschuss {NU/m:.4f})")
p(f"   Lesart heute: ϖ_R(Zelle)/Λ = √(6/ν) = {math.sqrt(6/NU):.5f}; r*(Zelle)/2Λ = √(6/ν)/2 = {math.sqrt(6/NU)/2:.5f}")
p("   Epoche: ΣΔs_b mitbewegt fest ⇒ Λ fest wie in der Lesart (Faktor m statt ν); a₀(z) und Ω_s(z) wie dort.")

# ------------------------------------------------------------------ K3
kopf("3  K3  Λ/L = 1 − w_Zelle = e^{−√(6HΛ/c)}   (Rückstellung schließt die Rückkehr)")
p("   v = Λ/L, B = 6HL/c:  −ln v = √(B·v)  hat für jedes B > 0 genau eine Lösung (links +∞, bei v = 1 negativ, monoton).")


def K3_loes(L, H):
    B = 6 * H * L / c
    v = brentq(lambda v: -math.log(v) - math.sqrt(B * v), 1e-300, 1.0)
    return v * L, v, B


p("   L                B           v = Λ/L      Λ [Gpc]     a₀ [10⁻¹⁰]   ν_impl = Λ/ΣΔs_b")
for name, L in [("L_min", Lmin), ("d_LSS", dLSS), ("L(1,1)", L11)]:
    Lam, v, B = K3_loes(L, H0)
    p(f"   {name:7s} {L/Gpc:6.2f}   {B:9.3f}   {v:.5f}    {Lam/Gpc:8.4f}    {K.a0_aus_schliesslaenge(Lam)*1e10:8.4f}     {Lam/sds_b(L):.4f}")
L_noetig3 = Lam11 / eins_minus_w(H0, Lam11)
p(f"   verlangt für a₀ = 1,1: L = Λ·e^(√17,86) = {L_noetig3/Gpc:.1f} Gpc ({L_noetig3/dLSS:.1f}·d_LSS); dann ν_impl = Λ/ΣΔs_b = "
  f"{Lam11/sds_b(L_noetig3):.3e} (gegen {NU:.3f})")
# mit der Lesart: Λ = ν·ΣΔs_b(L) = Gρ_mL³/(3c²) und Λ/L = e^{−√(6HΛ/c)}
g3 = lambda L: (K.G * K.rho_m_mitbewegt() * L**3 / (3 * c**2)) / L - eins_minus_w(H0, K.G * K.rho_m_mitbewegt() * L**3 / (3 * c**2))
L_fp = brentq(g3, 1 * Gpc, 100 * Gpc)
Lam_fp = K.G * K.rho_m_mitbewegt() * L_fp**3 / (3 * c**2)
p(f"   mit der Lesart (ν eingesetzt): Fixpunkt L = {L_fp/Gpc:.3f} Gpc (= {L_fp/Lmin:.3f}·L_min), Λ = {Lam_fp/Gpc:.4f} Gpc, "
  f"a₀ = {K.a0_aus_schliesslaenge(Lam_fp):.4e}")
# cosh-Form
p("   Nebenform mit Abschluss bei Λ (ϖ = √(6ΣΔs_b·Λ)):")
for name, L in [("L_min", Lmin), ("d_LSS", dLSS)]:
    S = sds_b(L)
    h = lambda Lam: Lam / L - eins_minus_w_cosh(H0, Lam, S)[0]
    Lam_c = brentq(h, 1e-6 * Gpc, L)
    vr = eins_minus_w_cosh(H0, Lam_c, S)[1]
    p(f"   {name}: Λ = {Lam_c/Gpc:.4f} Gpc, ϖ/Λ = {vr:.4f}, a₀ = {K.a0_aus_schliesslaenge(Lam_c):.4e}"
      f"{'  (ϖ > Λ: Tiefe jenseits der Schließung, Form ungültig)' if vr > 1 else ''}")
# Epoche: mitbewegte Zelle, L_phys = L̄/(1+z), H(z)
p("   Epoche (Zelle mitbewegt, L_phys = L̄/(1+z), H = H(z)), L̄ = d_LSS; Θ = 2πa₀/(cH), Ω_s = (π²/45)Θ²:")
p("   z        Λ(z)/Λ(0)     a₀(z)/a₀(0)    Θ(z)          Ω_s(z)        Lesart: a₀(z)/a₀(0), Ω_s(z)")
Lam0 = K3_loes(dLSS, H0)[0]
for z in (0, 0.5, 1, 3, 10, 100, 1090):
    Hz = H0 * K.E_hintergrund(z)
    Lz = K3_loes(dLSS / (1 + z), Hz)[0]
    a0z = K.a0_aus_schliesslaenge(Lz)
    Th = K.theta_(a0z, Hz)
    a0L = A0 * K.a0_frueh_relativ(z)
    ThL = K.theta_(a0L, Hz)
    p(f"   {z:<7}  {Lz/Lam0:.4e}    {a0z/K.a0_aus_schliesslaenge(Lam0):.4e}    {Th:.4e}    {K.K_S*Th**2:.4e}    {a0L/A0:.4f}, {K.K_S*ThL**2:.3e}")
Hz = H0 * K.E_hintergrund(1090); Lz = K3_loes(dLSS / 1091, Hz)[0]
f1090 = K.K_S * K.theta_(K.a0_aus_schliesslaenge(Lz), Hz)**2
p(f"   Vakuumteil bei z = 1090 als fester Anteil gelesen: ΔN_eff ≈ {K.delta_N_eff(f1090):.2f} (gemessen N_eff = 2,99 ± 0,17 [Kern])")

# ------------------------------------------------------------------ Blickrichtung
kopf("4  Beobachterfehler und Brücken")
p("   K1: ΣΔs_b [L summiert, Zelle, epochenfest] / Λ [I]  =  1 − w_Zelle [I-Rückkehr mit G-Verlust 3H, Epoche].")
p("       Brücke: Kac-Maß (zeitfrei, Maß) = Laplace-Gewicht (zustandsabhängig über H). Nicht hergeleitet;")
p("       setzt eine feste Größe einer Epochengröße gleich (Zustandsfehler). Widerspricht zudem 03: Die Mode zählt")
p("       im mitatmenden Hintergrund (w ≈ 1) mit 1 + Q wie Staub; das Gewicht trennt nicht gefaltet/ungefaltet.")
p("   K2: ϖ_R(Zelle) [I, EINER Falte] = Λ [I]. Richtung stimmt, aber die Summe ΣΔs_b [L, viele Falten] wird als")
p("       eine Falte mit einer Halbgeraden gelesen (G als L); nach 07 hat jede Falte ihre eigene Halbgerade bis Λ.")
p("       Brücke: Leiterende (01/08) ist dort eine Aussage über M_U, nicht über die Zelle. Nicht hergeleitet;")
p("       dazu verletzt (a) mit ν = 1 + Q_KOSM das Kac-Lemma (ν·ΣΔs_b > Λ).")
p("   K3: Λ [I] / L [G, Raumschleife] = 1 − w_Zelle. Liest Λ als Bruchteil einer Raumlänge (I als G), das")
p("       verbietet 03 (Λ keine Raumlänge); |ΔL_Raum| = Λ gilt nur dem Betrag nach aus der Bilanz. Mit H(z) und L ∝ a")
p("       atmet Λ mit (früh a₀/H nicht → 0). Nicht hergeleitet.")
p("   Lesart Λ = ν·ΣΔs_b: Kac-Lemma + kein Außen (I, G), ν = Rückkehrzeit des Bereichs ΣΔs_b (Zusatzannahme X_d).")

# ------------------------------------------------------------------ Ergebnis
kopf("5  Ergebnis")
p("   K1  keine Lösung für L ≥ 15,9 Gpc (Kreispaare verlangen L ≥ 27,4); mit ν eingesetzt a₀ = 5,7·10⁻¹⁰, L = 15,9 Gpc;")
p("       cosh-Form: eine Lösung (a₀ = 2,2·10⁻¹⁰), aber mit ϖ/Λ = 1,40 > 1 (ungültig).")
p("   K2  (a) a₀ = 1,195·10⁻¹⁰ bei L_min (würfelförmig), L(1,1) = 28,1 Gpc; ν_impl = 6 (Q = 5) verletzt mit Q_KOSM das")
p("       Kac-Lemma; (b) ν_impl = 1,5. Legt Λ nicht ohne L fest; ersetzt nur ν durch eine Zahl.")
p("   K3  a₀ = 4,4·10⁻¹⁰ (d_LSS; Ω_s heute 3,9 > 1); für 1,1 verlangt L = 906 Gpc mit ν_impl ≈ 2·10⁻⁴; mit ν eingesetzt")
p("       L = 15,9 Gpc < L_min, a₀ = 5,8·10⁻¹⁰; cosh-Form 2,0–2,1·10⁻¹⁰ mit ϖ > Λ (ungültig); atmet mit der Epoche")
p("       (a₀(1)/a₀(0) = 1,9, Ω_s(1090) = 0,32, ΔN_eff ≈ 3,6).")
p("   Kein Kandidat trägt. Die Lesart Λ = ν_kosm·ΣΔs_b bleibt Lesart; ν und L bleiben offen.")

# ------------------------------------------------------------------ Prüfungen
kopf("6  Unabhängige Prüfungen (zweiter Weg)")
# K1: dichtes Gitter statt Nullstellensuche
S = sds_b(dLSS)
us = np.geomspace(1.0 + 1e-9, 1e6, 200001)
m1 = np.max(np.log(us) / np.sqrt(us))
p(f"   K1 Gitter: max ln u/√u = {m1:.6f} (2/e = {2/math.e:.6f}); b(d_LSS) = {6*H0*S/c:.4f} > b_krit ⇒ keine Lösung")
# K3: Fixpunktiteration Λ ← L·e^{−√(6HΛ/c)} gedämpft
L = dLSS; Lam = 0.1 * L
for _ in range(2000):
    Lam = 0.5 * Lam + 0.5 * L * math.exp(-math.sqrt(6 * H0 * Lam / c))
p(f"   K3 gedämpfte Iteration (d_LSS): Λ = {Lam/Gpc:.4f} Gpc (brentq {K3_loes(dLSS, H0)[0]/Gpc:.4f})")
# K2: aus a0_aus_zelle mal ν/6
p(f"   K2(a) über den Kern: a0_aus_zelle(L_min)·ν/6 = {K.a0_aus_zelle(Lmin)*NU/6:.4e}")

txt = out.getvalue()
print(txt)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "lambda_bruecke.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt)
