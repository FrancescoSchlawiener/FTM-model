# -*- coding: utf-8 -*-
"""
Rate der neuen Faltung an der Engstelle: Schwarze Löcher (lokale Engstelle r₀ = 6Δs) und Λ.
Blick: Λ = ν_kosm·ΣΔs_b(Zelle) ist I/G (innere Richtung, Zelle); die Falte r ≤ r₀ ist L (am Ort).
Entropie S = 4πk·α_G(M) je Falte, S_Λ = (45/π³)k(Λ/ℓ_P)² je Horizont. Nichts davon ist Dichte am Ort.
1 Λ unter Akkretion und Verschmelzung (Δs additiv, ν, Abstrahlung).  2 Raten und w(z) im Vergleich.
3 Rückgabe (Hawking): Zeitskalen und Bilanz.  4 Beobachterfehler.  5 von Hand.
Literaturwerte [G] = aus dem Gedächtnis, nicht nachgeprüft, nur Größenordnung.
"""
import math, os, sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, c, G, hbar, kB, l_P, Mpc, Gpc, yr, Gyr, M_sun, T_CMB, OMEGA_M, OMEGA_R, OMEGA_B,
                      OMEGA_DE, Q_KOSM, H0_si, schliesslaenge, omega_vakuum, horizontentropie,
                      kreislauf_mindestrate, isco_kerr, alpha_G, ds, r0, r_stern, hawking_lebensdauer,
                      lebensdauer_falte, grenzmasse_hawking, E_hintergrund)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def H_(t):
    P(); P("=" * 110); P(t); P("=" * 110)

H0 = H0_si()
Lam0 = schliesslaenge(A0)
Os = omega_vakuum(A0)
Ox = OMEGA_DE - Os
nu_k = 1 + Q_KOSM
rho_c = 3 * H0**2 / (8 * math.pi * G)
rho_b = OMEGA_B * rho_c
eps_min, S_Lam, H_s = kreislauf_mindestrate(Lam0)
SDOT = 7.1e86                                   # Ṡ_i je Horizont heute (Wert des Modells, 03)
S_OBS = 3.1e104                                 # Entropie der Falten im beobachtbaren Kosmos [G] (Egan & Lineweaver 2010)

# [G] Masse der großen Falten je Volumen heute und Form der Wachstumsgeschichte (nur Größenordnung)
RHO_BH_MSUN_MPC3 = (3e5, 4.5e5, 6e5)            # [G] Galaxienkern-Falten, lokale Dichte
def psi_form(z):                                # [G] Form der Stern-/Akkretionsgeschichte (Madau & Dickinson 2014), nur Form
    return (1 + z)**2.7 / (1 + ((1 + z) / 2.9)**5.6)

def Hz(z):
    return H0 * E_hintergrund(z)

# Normierung: ∫ Ṁ dt von z = 20 bis 0 = ρ_BH(0)
_I = quad(lambda z: psi_form(z) / ((1 + z) * Hz(z)), 0, 20, limit=200)[0]
def rhodot_bh(z, rho0):                         # kg/m³/s, mitbewegt
    return rho0 * psi_form(z) / _I
def rho_bh(z, rho0):
    return rho0 * quad(lambda x: psi_form(x) / ((1 + x) * Hz(x)), z, 20, limit=200)[0] / _I
def kgm3(x):
    return x * M_sun / Mpc**3

H_("1  Λ = ν_kosm·ΣΔs_b(Zelle): was ändert sich, wenn Baryonen in Falten (r ≤ r₀) fallen?")
P(f"   Λ = c²/2a₀ = {Lam0/Gpc:.2f} Gpc,  ν_kosm = 1 + Q = {nu_k:.3f},  ρ_b = {rho_b:.3e} kg/m³")
M = 1e8 * M_sun
P(f"   (a) Δs additiv: ΣΔs = GΣM/3c². Umlagerung von Masse in die Falte ändert ΣM der Zelle nicht: dΣΔs = 0.")
P(f"   (b) Rückkehrtiefe: r*(M_b) der Galaxie hängt nur an ihrer Gesamtmasse (Gauß, kein Außen); ob ein Teil in r₀ liegt,")
P(f"       sieht die Ebene darüber nicht. Beispiel 10⁸ M☉: r₀ = {r0(M)/1.496e11:.2f} AE, r* = {r_stern(M)/(1e3*3.0857e16):.2f} kpc, r*/r₀ = {r_stern(M)/r0(M):.2e}")
P(f"       ν_kosm ist Kac-Zählung in der inneren Richtung der Zelle (I/G), kein Mittel örtlicher ν: dν_kosm = 0.")
aN_r0 = G * M / r0(M)**2
P(f"       [L-Lesart, falsch] ν am Ort r₀: a_N = c⁴/(36GM) = {aN_r0:.2e} m/s² = {aN_r0/A0:.1e}·a₀ ⇒ ν = 1;")
P(f"       setzte man dort ν = 1 statt ν_kosm, fiele Λ um (Q/(1+Q))·dM_BH/M_b = {Q_KOSM/nu_k:.3f}·dM_BH/M_b.")
P(f"   (c) Abstrahlung: Beim Einfall wird der Anteil η der Ruhemasse Strahlung; Strahlung verliert ihre Faltung durch")
P(f"       Rotverschiebung (03) ⇒ dΛ/Λ = −η·dM_BH/M_b (Obergrenze, der Verlust läuft mit 1/a nach).")
for j in (0.0, 0.7, 0.998):
    r_isco = isco_kerr(j)
    # Bindung auf der ISCO (Kerr, mitläufig): E = (r^1.5 − 2r^0.5 + j)/(r^0.75·√(r^1.5 − 3r^0.5 + 2j))
    E = (r_isco**1.5 - 2 * r_isco**0.5 + j) / (r_isco**0.75 * math.sqrt(r_isco**1.5 - 3 * r_isco**0.5 + 2 * j))
    P(f"       j = {j:5.3f}: r_ISCO = {r_isco:.3f}·3Δs = {r_isco/2:.3f}·r₀,  η = 1 − E_ISCO = {1-E:.4f}")
P(f"       Verschmelzung gleicher Falten: S ∝ M² und S_f ≥ S₁ + S₂ ⇒ M_f ≥ √2·M₁, abgestrahlt ≤ 1 − 1/√2 = {1-1/math.sqrt(2):.3f}")
P(f"   ⇒ Einfall und Verschmelzung vergrößern Λ nie; sie senken es um den abgestrahlten Anteil.")

H_("2  Raten im Vergleich und w(z)")
P(f"   Mindestrate aus der Entropiebilanz (04): ε = Ṡ_i/(2S_Λ) = {eps_min:.3e} s⁻¹ (S_Λ = {S_Lam:.3e} k)")
for rb in RHO_BH_MSUN_MPC3:
    r0b = kgm3(rb)
    P(f"   ρ_BH(0) = {rb:.1e} M☉/Mpc³ [G]: ρ_BH/ρ_b = {r0b/rho_b:.2e},  Ṁ/M_b heute = {rhodot_bh(0, r0b)/rho_b:.2e} s⁻¹,"
      f"  Mittel ρ_BH/(ρ_b t₀) = {r0b/rho_b/(13.8*Gyr):.2e} s⁻¹")
r0b = kgm3(RHO_BH_MSUN_MPC3[1])
fdot0 = rhodot_bh(0, r0b) / rho_b
eta = 0.1
P(f"   Abstrahlung (η = {eta}): dΛ/Λdt = −η·Ṁ/M_b = {-eta*fdot0:.2e} s⁻¹  (Betrag = {eta*fdot0/eps_min:.1e}·ε, Vorzeichen falsch)")
P(f"   Obergrenze, als ob jede eingefallene Masse Λ neu um ν·Δs verlängerte (Doppelzählung): +{fdot0:.2e} s⁻¹ = {fdot0/eps_min:.1e}·ε")
# Entropie-Weg: Ṡ der Falten, S ∝ M²
sdot_bh = 2 * S_OBS * rhodot_bh(0, r0b) / r0b
P(f"   Entropie-Weg: Ṡ_Falten ≈ 2S_obs·Ṁ_BH/M_BH = {sdot_bh:.2e} k/s [G] (Modellwert Ṡ_i = {SDOT:.1e}: von Falten getragen)")
P(f"     ⇒ ε_BH = Ṡ_Falten/(2S_Λ) = {sdot_bh/(2*S_Lam):.2e} s⁻¹ = Mindestrate; keine Faltung, nur Maß des Zustandsraums.")
# nötige Rate für die gemessene Form, wenn der Vakuumteil allein variiert (Ω_x fest)
P(f"   Nötig, wenn der Vakuumteil (Ω_s = {Os:.3f} von Ω_DE = {OMEGA_DE:.3f}) allein die Form w₀, w_a [G] trüge (Ω_x fest):")
P(f"     d ln Λ/d ln a = −½·d ln ρ_s/d ln a = (3/2)(Ω_DE(a)/Ω_s)(1 + w(a)),  w(a) = w₀ + w_a(1 − a)")
for w0, wa in ((-0.75, -0.9), (-0.8, -0.75), (-0.85, -0.6)):
    need0 = 1.5 * (OMEGA_DE / Os) * (1 + w0)
    P(f"     w₀ = {w0}, w_a = {wa}: heute d ln Λ/d ln a = {need0:.3f} ⇒ Λ̇/Λ = {need0*H0:.2e} s⁻¹ = {need0*H0/eps_min:.1e}·ε;"
      f"  ρ_DE-Maximum verlangt Vorzeichenwechsel bei 1 + w = 0, z = {1/(1 - (-1 - w0)/wa) - 1:.2f}")
P(f"   Die übrigen zwei Drittel (Ω_x = {Ox:.3f}) wachsen mit Λ nicht: Ω_s ∝ Λ⁻² kann nur fallen, wenn Λ wächst.")
P()
P("   w(z) mit Λ an die Akkretion gebunden (Form ψ(z) [G], ρ_BH(0) = 4,5·10⁵ M☉/Mpc³ [G]);")
P("   w_s = −1 + (2/3)d ln Λ/d ln a,  w_DE = −1 + (Ω_s/Ω_DE)(1 + w_s) mit Ω_x fest:")
P("        z    Ṁ/(M_b H)     w_s − (−1): Abstrahlung η=0,1   Obergrenze (Doppelzählung)    w_DE − (−1) Obergrenze")
for z in (0.0, 0.4, 1.0, 2.0, 3.0, 5.0):
    x = rhodot_bh(z, r0b) / rho_b / Hz(z)            # d(M_BH/M_b)/d ln a
    a = 1 / (1 + z)
    os_a = Os / (Os + Ox)                            # Ω_s-Anteil an der Dunklen Energie (beide fest)
    P(f"   {z:6.1f}   {x:.3e}      {-(2/3)*eta*x:+.3e}                 {(2/3)*x:+.3e}                  {os_a*(2/3)*x:+.3e}")
P("   Die Form verlangt vor z ≈ 0,36 schrumpfendes, danach wachsendes Λ; die Akkretion hat ein Vorzeichen.")
P("   Maximum der Akkretion bei z ≈ 2, nicht bei 0,4; Betrag 10⁻⁵ statt 10⁻¹; Vorzeichen im Modell negativ.")

H_("3  Rückgabe: Faltung → Engstelle → Strahlung. Zeitskalen und Bilanz")
T_s = hbar * H_s / (2 * math.pi * kB)
M_gr_s = grenzmasse_hawking(T_s)
P(f"   Horizonttemperatur im Endzustand bei festem Λ: T_s = ħH_s/(2πk) = {T_s:.2e} K; Falten mit M < {M_gr_s/M_sun:.1e} M☉ geben zurück")
P(f"   heute (T_CMB) erst unter {grenzmasse_hawking(T_CMB):.2e} kg; Beginn der Rückgabe 4·10¹¹–10¹² a (04)")
for m in (10, 1e6, 1e9, 1e11):
    Mm = m * M_sun
    t1, t2 = hawking_lebensdauer(Mm), lebensdauer_falte(Mm)
    P(f"   M = {m:8.0e} M☉: S/k = 4πα_G = {4*math.pi*alpha_G(Mm):.2e},  Lebensdauer {t1/yr:.2e} a (3840·(S/k)·Δs/c: {t2/yr:.2e} a)")
P(f"   1/ε = {1/eps_min/yr:.2e} a  ≪ Lebensdauer jeder Falte ab {(1/eps_min*hbar*c**4/(5120*math.pi*G**2))**(1/3)/M_sun:.1e} M☉")
P("   Entropie der Rückgabe (Schwarzkörper, dE = c²dM, S_rad = (4/3)∫dE/T_H):")
Mm = 1e6 * M_sun
Srad = quad(lambda m: (4/3) * c**2 * 8 * math.pi * G * kB * m / (hbar * c**3) / kB, 0, Mm)[0]
P(f"     S_rad/S_Falte = {Srad/(4*math.pi*alpha_G(Mm)):.4f} (= 4/3): Rückgabe erzeugt Entropie, sie kehrt die Einfaltung nicht um.")
P("   Bilanz der Faltung: Einfall ΔΣΔs = 0; Rückgabe macht die Masse zu Strahlung, deren Faltung mit 1/a vergeht;")
P(f"   Λ fällt am Ende höchstens um den Anteil der Falten an den Baryonen (heute {r0b/rho_b:.1e}, obere Grenze 1).")
P(f"   Die Entropiebilanz verlangt in derselben Zeit Wachstum Λ² = Λ0²(1 + 2εt): nach 10⁶⁷ a Λ/Λ0 = {math.sqrt(1+2*eps_min*1e67*yr):.1e}.")
P("   ⇒ Der lokale Kreislauf schließt in der Form (flach → Wirbel → Falte → Strahlung → flach), nicht in Λ.")

H_("4  Beobachterfehler")
P("   Λ, ν_kosm: I/G (innere Richtung, Zelle), zeitfest außer über ΣΔs_b. Falte, r₀, T_H, Lebensdauer: L.")
P("   S_Falte = 4πkα_G(M) gilt je Falte; S_Λ und S_H je Horizont eines Beobachters; verglichen wird nur ΣS_Falten ≤ S_Λ im selben Horizont.")
P(f"   ΣS_Falten ≤ Ω_b²·S_H: heute S_obs/S_Λ = {S_OBS/S_Lam:.1e}; die Falten können S_Λ nie füllen, Λ muss nicht ihretwegen wachsen,")
P("   sondern damit das Maß des Zustandsraums schneller wächst als die Entropie (Mindestrate ε).")
P("   ν = 1 am Ort r₀ ist örtlich (L); als Beitrag zu ν_kosm gelesen wäre es die Verwechslung L statt I.")

H_("5  Von Hand")
P("   ρ_BH = 4,5·10⁵·1,989·10³⁰ kg/(3,0857·10²² m)³ = 8,95·10³⁵/2,938·10⁶⁷ = 3,05·10⁻³² kg/m³")
P(f"   ρ_b = 4,20·10⁻²⁸ ⇒ ρ_BH/ρ_b = 7,3·10⁻⁵; Mittel über t₀ = 4,35·10¹⁷ s: 1,7·10⁻²² s⁻¹; mal η = 0,1: 1,7·10⁻²³ s⁻¹")
P(f"   Mittel/ε = 1,7·10⁻²³/3,83·10⁻³⁷ = 4,4·10¹³ (heute nach ψ: 1,1·10¹³).   Kern: ρ_BH/ρ_b = {r0b/rho_b:.3e}, Mittel·η = {0.1*r0b/rho_b/(13.8*Gyr):.3e} s⁻¹")
P("   η (j = 0): E_ISCO = √(8/9) = 0,9428 ⇒ η = 0,0572")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
