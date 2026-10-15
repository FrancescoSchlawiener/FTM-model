# -*- coding: utf-8 -*-
"""
Bilanz aller Glieder der Ausdehnungsgleichung und Beobachterfehler in der Messung der Dunklen Energie.
Frage: Gibt es im Modell EIN Glied mit positiver, heute fester oder wachsender Dichte neben dem Vakuumteil?
Trägt das Mittel der Strukturwirkung ⟨ℒ_𝒮⟩ über gebundene Bereiche? Verzerrt der Ort des Beobachters
(Takt δ, gebundener Bereich w = 0, Licht sieht den wirbelfreien Teil) die gelesenen Abstände so, dass w(z) entsteht?
Blick: Ausdehnung = G (Zelle, Nullmode); Q_s am Ort = L, gemittelt über die Zelle = G; Λ, Kac-Zählung = I.
1 Inventur aller Glieder (Vorzeichen, Skalierung, heutiger Wert).
2 ⟨ℒ_𝒮⟩: Mittel von Q_s über gebundene Bereiche (Punktmassen bis R_ta), lineares Feld dazwischen.
3 Messung: Takt am Ort des Beobachters, gebundener Bereich, mittlere Konvergenz; gelesenes Δw₀, Δw_a.
4 Von Hand.
Formeln/Konstanten aus ftm_kern. Nicht im Kern (hier nur als Lesart der Messung bzw. Integral über Kernformeln):
  Q_s(s) = ∫₀^s 4t³n_B(t)dt (Integral über n_B des Kerns), CPL-Lesung E(z; w₀, w_a) der Messung, Linsengewicht.
"""
import math, os, sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, kpc, Mpc, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B, Q_KOSM, H0_KMS, H0_si,
                      omega_vakuum, a0_modell, n_B, nu, r_stern, v_flach, gewicht_w, tau_D, umkehrradius_heute,
                      rho_m_mitbewegt, schliesslaenge, K_S)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 110); P(t); P("=" * 110)

H0 = H0_si()
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
QS_INF = 4 * math.pi**4 / 15
rho_c = 3 * H0**2 / (8 * math.pi * G)

def Q_s(s):
    """Q_s(s) = ∫₀^s 4t³ n_B(t) dt; Q̃ = Q_s − 4π⁴/15."""
    if s <= 0:
        return 0.0
    if s > 60:
        return QS_INF
    return quad(lambda t: 4 * t**3 * float(n_B(t)) if t > 0 else 0.0, 0, s, limit=200)[0]

# ============================================================== 1 Inventur
kopf("1  Inventur der Ausdehnungsgleichung  H²/H₀² = Σ Glieder (flacher Torus, ε = −ℒ für statische Glieder)")
th_bg = 3 * H0                                  # mitbewegt θ̃ = 3H
w_bg = gewicht_w(th_bg)
s_bg = math.sqrt(th_bg * tau_D())
nu_bg = 1 / (-math.expm1(-s_bg))                  # ν(s) = 1 + n_B(s) mit s² = θ̃τ_D
rest_w = (1 - w_bg) * (nu_bg - 1)
P("   Hintergrund: θ̃τ_D = %.2f, 1 − w = %.4f, (1 − w)(ν − 1) = %.2e (s² = θ̃τ_D, a_N = 0)" % (th_bg * tau_D(), 1 - w_bg, rest_w))
P("   Ω_s = (π²/45)Θ² = %.4f (a₀ = %.2e), gemessene Lücke Ω_DE = %.4f, fehlend Ω_x = %.4f (Anteil %.1f %%)"
  % (OM_S, A0, OMEGA_DE, OM_X, 100 * OM_X / OMEGA_DE))
P("   Nur mit Ω_s allein: Θ = %.3f nötig, a₀ = %.3e m/s² (cH₀/2π = %.3e); SPARC 1,04–1,19, Zelle ≤ 1,13"
  % (math.sqrt(OMEGA_DE / K_S), math.sqrt(OMEGA_DE / K_S) * a0_modell(), a0_modell()))

# Werte für Teil 2, hier vorweg für die Tabelle
v_sun = 200e3
glieder = []   # (Name, Vorzeichen, Skalierung, Wert heute, Blick)
P()
P("   %-44s %-6s %-26s %-14s %s" % ("Glied", "Vorz.", "Skalierung mit a", "Ω heute", "Blick / Zustand"))

# ============================================================== 2 ⟨ℒ_𝒮⟩
def J_punkt(sR):
    """∫_{r<R} Q_s(r*/r) 4πr²dr = 4π r*³ J, J = ∫_{s_R}^∞ Q_s(s)/s⁴ ds."""
    f = lambda s: Q_s(s) / s**4
    return quad(f, sR, 1.0, limit=200)[0] + quad(f, 1.0, 60.0, limit=200)[0] + QS_INF / (3 * 60.0**3)

rows_Q = []
for lgM in (8, 10, 11, 12, 13, 14, 15):
    Mb = 10**lgM * M_sun
    rs = r_stern(Mb)
    Rta = umkehrradius_heute((1 + Q_KOSM) * Mb)
    J = J_punkt(rs / Rta)
    vf = v_flach(Mb)
    # Energie je Objekt a₀²/(8πG)·4πr*³J = M_b v_f² J/2; alle gebundenen Baryonen in dieser Masse (F = 1):
    om_Q = OMEGA_B * (vf / c)**2 * J / 2
    rows_Q.append((lgM, rs / kpc, Rta / kpc, J, vf / 1e3, om_Q, om_Q / OM_S))

# lineares Feld zwischen Strukturen: a_N = (1–4)·10⁻⁴ a₀ (Modell 07, Forschungsstand)
QS_lin = [(x, Q_s(math.sqrt(x))) for x in (1e-4, 4e-4)]

om_Q_max = max(r[5] for r in rows_Q)
glieder = [
    ("Staub (1 + Q)ρ_b (Kac-Zählung)", "+", "a⁻³", "%.4f" % OMEGA_M, "G; Maß in I, fällt"),
    ("Strahlung", "+", "a⁻⁴", "%.2e" % OMEGA_R, "G; fällt"),
    ("Vakuumteil π⁴/(15Λ²)", "+", "a⁰ (Λ fest seit z≈10)", "%.4f" % OM_S, "I→G; fest, trägt 1/3"),
    ("Rest der Schließungsregel (1 − w)", "0", "in (1+Q)ρ_b, a⁻³", "0 (1−w=%.3f)" % (1 - w_bg), "Mode = Staub (c_s²=0)"),
    ("Gegenquelle/Nullmode σ = ΣΔs − ΔL = 0", "0", "ist das Staubglied", "0", "G; ⟨ρ − ρ̄⟩ = 0"),
    ("Scherung/Wirbel der Zelle", "+", "a⁻⁶", "≤ 1e-12", "G; (σ/H)² < 10⁻¹²"),
    ("⟨ℒ_𝒮⟩: −a₀²⟨Q_s⟩/8πG gebunden", "−", "∝ F(a)·a⁻³", "≥ −%.1e" % om_Q_max, "L gemittelt → G; wächst mit F"),
    ("⟨ℒ_𝒮⟩: lineares Feld a_N ≤ 4e-4 a₀", "−", "∝ (D/a²)^{3/2}, fällt", "−%.1e" % (OM_S * QS_lin[1][1] / QS_INF), "L gemittelt → G"),
    ("2h^{μν}B_μB_ν (Strom J)", "±", "∝ gebunden", "≤ %.0e" % (1e-9 * OMEGA_M), "Spannung 10⁻⁹ρc², Bilanz"),
    ("Arbeitsaustausch der Regel", "±", "mittelt weg", "≤ %.0e" % (3e-4 * (v_sun / c)**2 * OMEGA_M), "L; 3·10⁻⁴·σ² je 10 Mrd. J."),
    ("Gewicht (1 − w)(ν − 1)", "0", "nur Störungen", "0 (Faktor %.1e)" % rest_w, "wirkt auf a_N ≠ 0, Hintergrund a_N = 0"),
]
for g in glieder:
    P("   %-44s %-6s %-26s %-14s %s" % g)
P("   Summe der Glieder außer Staub/Strahlung/Vakuumteil: |ΣΩ| ≤ %.1e gegen nötig +%.3f." % (om_Q_max + 1e-9, OM_X))

kopf("2  ⟨ℒ_𝒮⟩ über gebundene Bereiche: ε_𝒮 = −a₀²Q̃/(8πG) = ε_s·(1 − ⟨Q_s⟩/(4π⁴/15))")
P("   Normierung: ℒ_𝒮 ⊃ (c⁴/16πG)(2/ℓ₀²)Q̃, ℓ₀ = c²/a₀; Q̃(0) = −4π⁴/15 gibt genau Λ_s = 4π⁴/(15ℓ₀²) (Vakuumteil).")
P("   Wo Feld ist (s > 0), steigt Q̃ gegen 0: Die Strukturwirkung SENKT den Vakuumteil am Ort (im Newton-Bereich ganz).")
P("   Punktmasse: s = r*/r, ∫Q_s dV = 4πr*³·J, J = ∫_{s_R}^∞ Q_s/s⁴ ds ≈ (4/3)ln(R/r*) + C; Energie je Objekt M_b·v_f²·J/2.")
P("   Obergrenze: alle Baryonen gebunden (F = 1) in Objekten der Masse M_b, Bereich bis R_ta.")
P("   %-8s %-9s %-10s %-7s %-9s %-12s %s" % ("M_b/M☉", "r* [kpc]", "R_ta [kpc]", "J", "v_f km/s", "|Ω_𝒮|", "|Ω_𝒮|/Ω_s"))
for r in rows_Q:
    P("   1e%-6d %-9.1f %-10.0f %-7.2f %-9.0f %-12.2e %.2e" % r)
C_fit = rows_Q[3][3] - 4 / 3 * math.log(rows_Q[3][2] / rows_Q[3][1])
P("   J − (4/3)ln(R_ta/r*) = %.3f (Konstante des inneren Teils, M_b = 10¹² M☉)" % C_fit)
for x, qs in QS_lin:
    P("   Lineares Feld a_N = %.0e·a₀: Q_s = %.2e, Anteil am Vakuumteil %.1e" % (x, qs, qs / QS_INF))
P("   Vorzeichen: negativ (senkt ε). Verlauf: ∝ F(a)·a⁻³ je Objektmasse, Betrag wächst mit der Bindung → Dichte der")
P("   Dunklen Energie fiele mit der Zeit (w > −1, richtiger Quadrant), aber um ≤ %.0e·Ω_s. Trägt nicht (nötig +%.2f·Ω_s)."
  % (om_Q_max / OM_S, OM_X / OM_S))
P("   Auch mit umgekehrtem Vorzeichen gelesen bliebe der Betrag ≤ %.0e: keine Lesart der Strukturwirkung trägt." % om_Q_max)

# ============================================================== 3 Messung
kopf("3  Beobachterfehler in der Messung: Takt am Ort, gebundener Bereich (w = 0), mittlere Konvergenz")

def E_cpl(z, w0=-1.0, wa=0.0, Om=OMEGA_M, Ode=OMEGA_DE):
    a = 1 / (1 + z)
    return math.sqrt(Om / a**3 + OMEGA_R / a**4 + Ode * a**(-3 * (1 + w0 + wa)) * math.exp(-3 * wa * (1 - a)))

ZG = np.linspace(0.02, 1.5, 149)
def dL(z, **k):
    return (1 + z) * quad(lambda x: 1 / E_cpl(x, **k), 0, z)[0]
mu0 = np.array([5 * math.log10(dL(z)) for z in ZG])
eps = 1e-3
dmu_w0 = (np.array([5 * math.log10(dL(z, w0=-1 + eps)) for z in ZG]) - mu0) / eps
dmu_wa = (np.array([5 * math.log10(dL(z, wa=eps)) for z in ZG]) - mu0) / eps
A = np.vstack([np.ones_like(ZG), dmu_w0, dmu_wa]).T
def lese_w(dmu):
    x, *_ = np.linalg.lstsq(A, dmu, rcond=None)
    return x[1], x[2]

# (a) Takt am Ort des Beobachters: Tiefe im gebundenen Bereich (Regel bis R_ta, außen (1 + Q)·Newton)
def tiefe(Mb, r):
    Rta = umkehrradius_heute((1 + Q_KOSM) * Mb)
    innen = quad(lambda x: float(nu(G * Mb / x**2)) * G * Mb / x**2, r, Rta, limit=400)[0]
    aussen = (1 + Q_KOSM) * G * Mb / Rta
    return (innen + aussen) / c**2, Rta
d_MW, Rta_MW = tiefe(6e10 * M_sun, 8.2 * kpc)
R31 = umkehrradius_heute((1 + Q_KOSM) * 1e11 * M_sun)
d_M31 = tiefe(1e11 * M_sun, 780 * kpc)[0] if 780 * kpc < R31 else (1 + Q_KOSM) * G * 1e11 * M_sun / (780 * kpc * c**2)
delta_obs = d_MW + d_M31
P("   (a) Tiefe Φ/c² am Beobachter: Milchstraße (6·10¹⁰ M☉, 8,2 kpc, Regel bis R_ta = %.0f kpc) %.2e, Andromeda (10¹¹ M☉, 780 kpc) %.2e"
  % (Rta_MW / kpc, d_MW, d_M31))
P("       Summe δ_obs = %.2e. Wirkung: 1 + z_gel = (1 + z)(1 − δ_obs + δ_Quelle); H₀ in eigener Zeit um δ verschoben." % delta_obs)
dz = -delta_obs * (1 + ZG)
dmu_a = np.array([5 * math.log10(dL(z + d)) for z, d in zip(ZG, dz)]) - mu0
w0a, waa = lese_w(dmu_a)
P("       Δμ(z = 0,02 / 0,5 / 1,5) = %.1e / %.1e / %.1e mag; gelesen Δw₀ = %+.1e, Δw_a = %+.1e"
  % (dmu_a[0], dmu_a[np.argmin(abs(ZG - 0.5))], dmu_a[-1], w0a, waa))
P("       Quellen sitzen selbst in Mulden (δ_Quelle ~ δ_obs, zufällig): Streuung, kein Gang mit z über 10⁻⁶ hinaus.")

# (b) gebundener Bereich w = 0 um den Beobachter: Rand läuft mit H·R_ta, Kompensation bis 2. Ordnung
x2 = (H0 * Rta_MW / c)**2
mulde = 0.03 * (1 + Q_KOSM) * G * 6e10 * M_sun / (Rta_MW * c**2)
P("   (b) Bereich mit w = 0 (R_ta = %.2f Mpc): Rotverschiebung des nicht mitatmenden Wegs wird am Rand durch dessen Lauf H·R_ta"
  % (Rta_MW / Mpc))
P("       ausgeglichen; es bleiben (H₀R_ta/c)² = %.1e und die Mulde der Linse (3 %% von (1+Q)M_b an R_ta): %.1e." % (x2, mulde))
P("       Als Abstandsfehler bei z = 0,02: %.1e relativ; bei z = 0,5 %.1e. Kein Gang wie w(z)." % ((x2 + mulde) / 0.02, (x2 + mulde) / 0.5))

# (c) mittlere Konvergenz: Licht sieht den wirbelfreien Teil; Gauß auf dem Torus
def kappa_mittel(zs, epsL, F_z=lambda z: 1.0):
    chi = lambda z: quad(lambda x: 1 / E_cpl(x), 0, z)[0]
    xs = chi(zs)
    f = lambda z: (lambda x: x * (xs - x) / xs)(chi(z)) * (1 + z) / E_cpl(z) * epsL * F_z(z)
    return 1.5 * OMEGA_M * quad(f, 0, zs, limit=100)[0]
epsL = -0.03 * 0.85                       # Obergrenze: 3 % Mulde × gebundener Anteil 0,85, ohne Schließung
Fz = lambda z: 1 / (1 + z)                # Bindung wächst (Form nur zur Größenordnung)
dmu_c = np.array([-5 / math.log(10) * kappa_mittel(z, epsL, Fz) for z in ZG])
w0c, wac = lese_w(dmu_c)
P("   (c) Mittlere Konvergenz, falls das vom Licht gesehene Mittel um ε_L = %.3f unter (1 + Q)ρ̄_b läge (Obergrenze):" % epsL)
P("       Δμ(z = 0,5 / 1,5) = %.1e / %.1e mag; gelesen Δw₀ = %+.1e, Δw_a = %+.1e."
  % (dmu_c[np.argmin(abs(ZG - 0.5))], dmu_c[-1], w0c, wac))
P("       Auf dem geschlossenen Torus ist aber nur ρ − ρ̄ lösbar (Gauß, Σ Quelle = 0): ⟨κ⟩ = 0 genau; ein Fehlbetrag")
P("       fiele in die Nullmode, also ins Staubglied der Ausdehnung (Staubverlust, falsches Vorzeichen). Scheiben (η > 1)")
P("       berühren nur die Streuung. Ergebnis (c): 0 im Mittel.")

# (d) Größe, die eine Messverzerrung liefern müsste
h_mod = math.sqrt(OMEGA_M + OMEGA_R + OM_S)
def dL_mod(z):   # physikalische Dichten fest (ω_m aus der Strahlung, ω_s aus a₀): H₀,mod = h_mod·H₀
    return (1 + z) * quad(lambda x: 1 / math.sqrt(OMEGA_M * (1 + x)**3 + OMEGA_R * (1 + x)**4 + OM_S) , 0, z)[0]
dmu_d = np.array([5 * math.log10(dL_mod(z)) for z in ZG]) - mu0
w0d, wad = lese_w(dmu_d)
P("   (d) Nötig: Modell nur mit Ω_s (ω_m, ω_s fest) gibt H₀,mod = %.3f·H₀ = %.1f km/s/Mpc; Δμ(0,5) − Δμ(0,02) = %.3f mag,"
  % (h_mod, h_mod * H0_KMS, dmu_d[np.argmin(abs(ZG - 0.5))] - dmu_d[0]))
P("       linear projiziert Δw₀ = %+.2f, Δw_a = %+.2f (nur Richtung). Messeffekte: (a) |Δμ| ≤ 2·10⁻⁴ mag, |Δw| ≤ 10⁻³; (b) 10⁻⁸; (c) 0 → 3 Größenordnungen und mehr zu klein." % (w0d, wad))
P("       Einen Faktor %.2f in H₀ gäbe nur δ_obs ≈ %.2f (gemessen %.0e)." % (1 / h_mod, 1 - h_mod, delta_obs))

# ============================================================== 3b Blick je Glied
kopf("3b Blick (L/I/G) und Zustand je Glied")
for t in [
    "Staub (1+Q)ρ_b: G (Zelle, Nullmode); 1+Q ist Kac-Zählung in I, kein Stoff; Epoche: seit Rekombination w≈1 → Staub.",
    "Vakuumteil: Λ liegt in I, als Rückstellung des Ganzen in G gelesen; zustandsfest, solange Λ fest (seit z≈10).",
    "Schließungsregel (1−w): Z/G, zustandsabhängig (θ̃, Preisach); ändert Dynamik der Störungen, nicht das Maß.",
    "Q_s gebunden: am Ort (L) eine Regelgröße, keine Dichte; als Mittel über die Zelle (G) ins Atmen gelesen.",
    "B², Arbeitsaustausch: L, Bilanz über 𝒫L auf der Schließung; keine Energie am Ort.",
    "Scherung: G; Schleifen atmen gleich (01).",
    "Messung: δ_obs ist L am Beobachter (gebunden, w = 0), H und d_L sind G; Verwechslung nur auf 10⁻⁶-Stufe.",
]:
    P("   " + t)

# ============================================================== 4 Von Hand
kopf("4  Von Hand")
Mb = 1e15 * M_sun
vf = (G * Mb * A0)**0.25
P("   Haufen 10¹⁵ M☉: v_f = (GMa₀)^{1/4} = (1,327·10³⁵·1,1·10⁻¹⁰)^{1/4} = %.3e m/s; (v_f/c)² = %.2e" % (vf, (vf / c)**2))
Rta = umkehrradius_heute((1 + Q_KOSM) * Mb); rs = r_stern(Mb)
Jh = 4 / 3 * math.log(Rta / rs) + C_fit
P("   J ≈ (4/3)ln(R_ta/r*) + C = (4/3)·ln(%.0f/%.0f) + %.2f = %.2f (Code %.2f)" % (Rta / kpc, rs / kpc, C_fit, Jh, rows_Q[-1][3]))
P("   |Ω_𝒮| = Ω_b·(v_f/c)²·J/2 = 0,0492·%.2e·%.2f/2 = %.2e (Code %.2e)" % ((vf / c)**2, Jh, 0.0492 * (vf / c)**2 * Jh / 2, rows_Q[-1][5]))
qs_poly = None
try:
    import mpmath as mp
    s = 2.0; x = mp.e**(-s)
    qt = -4 * (s**3 * mp.polylog(1, x) + 3 * s**2 * mp.polylog(2, x) + 6 * s * mp.polylog(3, x) + 6 * mp.polylog(4, x))
    P("   Q̃(2) aus Li: %.6f; Q_s(2) − 4π⁴/15 per Integral: %.6f" % (float(qt), Q_s(2.0) - QS_INF))
except Exception:
    pass
P("   Tiefe Milchstraße: v_f²/c²·[1 + ln(R_ta/r) + (1+Q)v_N²/v_f²...] ≈ (%.0f km/s/c)²·(1 + ln(%.0f/8,2)) = %.1e (Code %.1e)"
  % (v_flach(6e10 * M_sun) / 1e3, Rta_MW / kpc, (v_flach(6e10 * M_sun) / c)**2 * (1 + math.log(Rta_MW / (8.2 * kpc))), d_MW))
P("   a₀ für Ω_s = Ω_DE: Θ = √(0,6865·45/π²) = √3,130 = 1,769; a₀ = 1,769·1,042·10⁻¹⁰ = 1,84·10⁻¹⁰ m/s²")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_bilanz_glieder.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
