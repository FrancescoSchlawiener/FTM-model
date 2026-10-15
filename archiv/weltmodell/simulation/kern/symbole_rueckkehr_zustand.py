# -*- coding: utf-8 -*-
"""
Symbole der Rückkehr (I), der Faltung am Ort (L) und des Zustandsraums (Z), die die Formeln in rueckkehr.py,
faltung.py und zustandsraum.py zusätzlich brauchen. Felder wie in symbole.py; eingebunden dort mit einer Importzeile.
Zeichen nach Register: Epizykelfrequenz κ_E und Vertikalfrequenz ω_z (nicht κ, ν_z), Entropie der Rücklaufzahl S_n
(nicht H), der bei Λ geschlossene Rückkehranteil als q_Λ (derselbe Nome q auf der Ebene I, geschlossen).
"""
from .symbole import SYMBOLE

REGISTERMODUL = True
_s = SYMBOLE.neu

# ------------------------------------------------------------------ Ebene I: Kern, Glättung, Schluss bei Λ (07)
_s("T_R", "T_rueck", "Rückkehrzeit in der inneren Richtung (Variable des Kerns, keine Zeit längs der Bahn)", "s", "I",
   "je_system", "Baustein 07")
_s("K(T)", "K_zeit", "Kern im Zeitbereich: Dichte erster Rückkehrzeiten mit dem Verlust Γ", "s⁻¹", "I", "je_system",
   "Baustein 07")
_s("Ω", "Omega_bahn", "Winkelgeschwindigkeit der Kreisbahn, Ω² = a/r", "s⁻¹", "L", "je_system", "Baustein 07")
_s("|K̂(iΩ)|/K̂(0)", "K_bahn", "Anteil des Kerns bei Bahnfrequenz gegen den stationären Kern", "1", "I", "je_system",
   "Baustein 07")
_s("k_ϖ", "k_welle", "Wellenzahl einer Mode des Feldes am Ort", "m⁻¹", "L", "je_system", "Baustein 07")
_s("P_ϖ", "P_glatt", "Glättung einer Mode über die Tiefe ϖ, e^{−|k|ϖ} (Halbgruppe P_a·P_b = P_{a+b})", "1", "I",
   "je_system", "Baustein 07")
_s("λ_Λ", "lambda_schluss", "kleinster Eigenwert der bei Λ geschlossenen Rückkehr, Dπ²/(4Λ²)", "s⁻¹", "I",
   "je_system", "Baustein 07")
_s("Λ/ϖ_R", "Lam_zu_varpiR", "Schließlänge gegen Rückkehrtiefe", "1", "I", "je_system", "Baustein 07")
_s("q_Λ", "q_schluss", "Rückkehranteil der bei Λ geschlossenen Rückkehr cosh(s(Λ/ϖ_R − 1))/cosh(sΛ/ϖ_R)", "1", "I",
   "je_system", "Baustein 07", stamm="q")
_s("kΛ", "k_Lam", "Schließlänge in Einheiten der Abschirmlänge, Λ·√(Γ/D)", "1", "I", "je_system", "Baustein 07")
_s("ε_Λ", "eps_schluss", "Abweichung der Regel von Tully-Fisher durch den Schluss bei Λ, 1/tanh(kΛ) − 1", "1", "I",
   "je_system", "Baustein 07")
_s("r_Λ", "r_schluss", "Radius, an dem Γ auf den kleinsten Eigenwert der geschlossenen Rückkehr fällt", "m", "LI",
   "je_system", "Baustein 07")
_s("T_I", "T_innen", "Periode 2Λ/c der inneren Richtung als Temperatur geschrieben, kT_I = ħa₀/c (keine Wärme)", "K",
   "I", "seit_z10", "Baustein 07")
_s("x_V", "x_feld", "Argument des Strukturteils im Energiefunktional, |∇Φ̄_N|²/a₀² = s⁴", "1", "L", "je_system",
   "Baustein 07")

# ------------------------------------------------------------------ Ebene L: Glättung und Hierarchie, Flussgesetz (07, 01)
_s("Ē^TF", "E_TF", "spurfreie Hesse-Form des geglätteten Elternfelds (3×3), gemittelt über die eigene Zone", "s⁻²",
   "L", "je_system", "Baustein 07")
_s("M×", "M_kreuz", "Masse, an der die Glättungsgrenze ℓ = R die Schwelle r* = R schneidet, a₀³/(G𝒯²)", "kg", "L",
   "je_system", "Baustein 07")
_s("R×", "R_kreuz", "Radius dieses Schnitts, a₀/𝒯", "m", "L", "je_system", "Baustein 07")
_s("ℓ/r*", "ell_zu_rstern", "Glättungslänge gegen Schwellenradius, (M×/M)^{1/6}", "1", "L", "je_system", "Baustein 07")
_s("δ*", "delta_stern", "Tiefe an der Schwelle (v_flach/c)², Schritt der Leiter 3Δs : r* : ℓ₀", "1", "LI", "je_system",
   "Baustein 01")
_s("d_eff", "d_eff", "Flussdimension einer Form der Dimension d unter der Regel 1 + (d − 1)(1 − χ/2); Kugel 1 − d ln a/d ln r = 3 − χ", "1", "L", "je_system", "Baustein 07")
_s("κ_E", "kappa_E", "Epizykelfrequenz der Kreisbahn", "s⁻¹", "L", "je_system", "Baustein 07")
_s("ω_z", "omega_z", "Vertikalfrequenz der Kreisbahn (im Modelltext ν_z)", "s⁻¹", "L", "je_system", "Baustein 07")
_s("(κ_E²+ω_z²)/Ω²", "epizyklen", "Summe der Epizyklen im Maß der Umlauffrequenz, 2 + χ", "1", "L", "je_system",
   "Baustein 07")

# ------------------------------------------------------------------ Ebene Z: Bose-Leiter, Rücklaufzahl, Schleifen, Gedächtnis (05, 06, 07)
_s("ℬ_d", "bose_stufe", "Stufe d der Bose-Leiter s^d·n_B(s)", "1", "Z", "je_system", "Baustein 07")
_s("ℐ_d", "bose_integral", "Integral der Stufe d über s, Γ(d+1)ζ(d+1)", "1", "Z", "fest", "Baustein 07")
_s("g_F", "g_fisher", "Fisher-Metrik der Rücklaufzahl −dν/ds = ν(ν − 1)", "1", "Z", "je_system", "Baustein 07")
_s("S_n", "S_n", "Entropie der Rücklaufzahl ln ν + χ (im Modelltext H)", "1", "Z", "je_system", "Baustein 07")
_s("a_F", "a_falte", "Parameter der Falte in der Normalform ẋ = −(x³ − a_F·x − c)", "1", "Z", "je_system", "Baustein 06")
_s("δ_F", "delta_falte", "Abstand des Arbeitspunkts zur Falte", "1", "Z", "je_system", "Baustein 06")
_s("k_K", "k_kopplung", "Stärke der Kopplung zweier Phasen", "s⁻¹", "Z", "je_system", "Baustein 06")
_s("φ_K", "phi_kopplung", "Verzögerungswinkel der Kopplung", "rad", "Z", "je_system", "Baustein 06")
_s("Δω_R", "d_omega_rast", "größter Frequenzabstand, bei dem zwei Zyklen einrasten, 2k|cos φ|", "s⁻¹", "Z",
   "je_system", "Baustein 06")
_s("u_P", "preisach_folge", "Folge der Eingangswerte eines Preisach-Gedächtnisses", "1", "Z", "je_system",
   "Baustein 06")
_s("α_P", "preisach_auf", "Aufwärtsschwelle je Hysteron", "1", "Z", "je_system", "Baustein 06")
_s("β_P", "preisach_ab", "Abwärtsschwelle je Hysteron", "1", "Z", "je_system", "Baustein 06")
_s("μ_P", "preisach_mass", "Gewicht je Hysteron (Preisach-Maß)", "1", "Z", "je_system", "Baustein 06")
_s("y_P", "preisach_ausgang", "Ausgang Σ μ_P·Zustand des Preisach-Gedächtnisses", "1", "Z", "je_system", "Baustein 06")
_s("n_R", "n_abgriff", "Glied der RC-Kette, dessen Spannung auf den Eingang zurückgeführt wird", "1", "Z", "je_system",
   "Baustein 06")
_s("RG", "RG_kette", "Längswiderstand mal Querleitwert je Glied der RC-Kette", "1", "Z", "je_system", "Baustein 06")
_s("k_w", "k_modulgewicht", "Gewicht einer Modulform (gerade, ≥ 4)", "1", "Z", "fest", "Baustein 05")
_s("G_k", "G_eisenstein", "Eisenstein-Reihe Σ′(m + nτ)^{−k}, 𝒫 auf dem Torus angewandt auf z^{−k}", "1", "Z",
   "je_system", "Baustein 05")

# ------------------------------------------------------------------ Labor im Zustandsraum: Kupferdraht, Hüllkurve (06)
_s("Δf_d", "df_anregung", "Verstimmung der Anregung gegen die Eigenfrequenz", "Hz", "L", "je_system", "Baustein 06")
_s("V_D", "V_draht", "Amplitude der Antriebsspannung", "V", "L", "je_system", "Baustein 06")
_s("ΔT_D", "dT_draht", "Erwärmung des heißen Asts im stationären Punkt", "K", "L", "je_system", "Baustein 06")
_s("Δf_F", "df_falte", "Verstimmung, an der der heiße Ast endet (Falte)", "Hz", "L", "je_system", "Baustein 06")
_s("Δf_H", "df_hopf", "Abstand der Hopf-Grenze vor der Falte des heißen Asts", "Hz", "L", "je_system", "Baustein 06")
_s("f_ring", "f_ring", "Ringfrequenz an der Hopf-Grenze", "Hz", "L", "je_system", "Baustein 06")
_s("τ_th/τ_m", "tau_th_zu_m", "Wärmezeit C_th/G_th gegen Ringzeit 2Q_l/ω₀ der Amplitude", "1", "L", "je_system",
   "Baustein 06")

# ------------------------------------------------------------------ Spielbeispiele und Rechenpunkte (04, 05, 06)
# Eigene Zeichen, weil γ und β im Modelltext die PPN-Koeffizienten sind, v die Kreisgeschwindigkeit, D die
# Rückkehrdiffusion.
_s("β_Page", "beta_page", "Parameter der Page-Kurve der Abstrahlung einer Falte (Anteil der Emission)", "1", "L",
   "fest", "Messwert [G] (Page 1976)")
_s("γ_K", "gamma_spielkern", "Kopplung des Maxwell-Spielkerns (Feder γ_K, Dämpfer γ_K·t_K)", "1", "Z", "je_system",
   "Baustein 05")
_s("v_K", "v_spielkern", "Anfangsgeschwindigkeit im Spielbeispiel (Einheiten GM = 1, x₀ = 1)", "1", "Z", "je_system",
   "Baustein 05")
_s("Δh", "delta_h", "Höhenunterschied im Vergleich der beiden Kanäle", "m", "L", "je_system", "Baustein 05")
_s("D_Z", "dissipation_Z", "Dissipation D(Z) der Torus-Dynamik T² × S (transversale Rate −D_Z)", "1", "Z", "je_system",
   "Baustein 06")

# ------------------------------------------------------------------ Labor: Kupferdraht, Kriterium der Hopf-Grenze und Zerlegung der Schleife (06, 09)
_s("ℓ₁", "lyapunov_1", "erster Lyapunov-Koeffizient der Hüllkurve an der Hopf-Grenze (Koordinaten Re X/|X₀|, "
   "Im X/|X₀|, ΔT/K); > 0 unterkritisch, < 0 überkritisch", "s⁻¹", "L", "je_system", "Baustein 06")
_s("q_σ", "q_sigma", "Anteil der Schleifenverstärkung am Draht über σ_el(T) allein", "1", "L", "je_system",
   "Baustein 09", stamm="q")
_s("q_E", "q_E", "Anteil der Schleifenverstärkung am Draht über E(T) allein", "1", "L", "je_system", "Baustein 09",
   stamm="q")
_s("q_N", "q_N", "Anteil der Schleifenverstärkung am Draht über die Wärmespannung N(T) allein", "1", "L", "je_system",
   "Baustein 09", stamm="q")

# ------------------------------------------------------------------ Zustandsraum: Resolvente, Spielkern, Rastpunkt, Diskriminante (05, 06)
_s("α_D", "alpha_drehung", "Drehzahl einer irrationalen Drehung des Kreises (golden: (√5 − 1)/2)", "1", "Z", "fest",
   "Baustein 05")
_s("K_M", "k_moden", "Zahl der Fluchtmoden je Umlaufsinn (Wellenzahlen 1 … K_M)", "1", "Z", "je_system", "Baustein 05")
_s("ν_F", "nu_flucht", "größte Verstärkung der Fluchtmoden 1/|1 − q·e^{2πikα_D}| der Resolvente (I − qT)⁻¹", "1", "Z",
   "je_system", "Baustein 05")
_s("x_U", "x_umkehr", "Umkehrpunkt des Spielbeispiels mit Maxwell-Spielkern (Einheiten GM = 1, x₀ = 1)", "1", "Z",
   "je_system", "Baustein 05")
_s("c_Δ", "c_delta", "Koeffizienten der q-Reihe der Diskriminante Δ = Σ c_Δ(n)qⁿ (Ramanujan)", "1", "Z", "fest",
   "Baustein 05")
_s("Δ_τ", "delta_modul", "Diskriminante Δ(τ) = q∏(1 − qⁿ)²⁴ des Torus ℂ/(ℤ + τℤ)", "1", "Z", "je_system",
   "Baustein 05")
_s("λ_ψ", "lambda_rast", "Rate am Rastpunkt der Phasendifferenz, −√((2k cos φ)² − Δω²)", "s⁻¹", "Z", "je_system",
   "Baustein 06")
_s("N_K", "n_zyklen", "Zahl gekoppelter Zyklen an einem gemeinsamen Medium (Kuramoto)", "1", "Z", "je_system",
   "Baustein 06")
_s("(N−1)k₁/W", "einrasten_erst", "Kopplung beim ersten Einrasten im Maß der Bandbreite W, mal N − 1", "1", "Z",
   "je_system", "Baustein 06")
_s("λ_⊥", "lambda_quer", "transversale Rate des Torus T² × S, −D_Z", "1", "Z", "je_system", "Baustein 06")

# ------------------------------------------------------------------ Rückkehr: Gitter-Rückkehr, Drehung, Windungskorrelation (04, 07)
_s("G_d(1)", "watson_G", "größte Verstärkung der Gitter-Rückkehr in d Dimensionen, G_d(1) = 1/(1 − p_d)", "1", "I",
   "fest", "Baustein 07")
_s("T_μ", "t_drehung", "Rückkehrzeiten (Schritte) einer irrationalen Drehung in einen Bogen vom Maß μ(A)", "1", "I",
   "je_system", "Baustein 07")
_s("N_H", "n_harmonische", "Zahl der Harmonischen des Dreiecks in der Windungskorrelation", "1", "I", "je_system",
   "Baustein 04")
_s("S_CHSH", "s_chsh", "größte CHSH-Größe einer drehinvarianten Korrelation E(Δ)", "1", "I", "je_system",
   "Baustein 04")

# ------------------------------------------------------------------ Faltung: Rückgabe, zwei Kanäle (04, 05)
_s("Ṁ_F", "m_punkt_falte", "Rückgabe der Falte dM/dt", "kg/s", "L", "je_system", "Baustein 04")
_s("R_K", "r_kanal", "Verhältnis differentieller Kanal (Vorspannung → Resonanz) gegen Takt-Kanal bei radialer "
   "Verschiebung", "1", "L", "je_system", "Baustein 05")

# ------------------------------------------------------------------ Labor: VO₂ auf einer Mikrobrücke im Magnetfeld (06)
_s("𝔄_VO₂", "aufbau_vo2", "Satz der Annahmen des VO₂-Aufbaus (Übergang, Widerstände, Wärme, Kreis, Brücke, Feld; "
   "Einheiten µs, K, V, kΩ, mA, nF, nJ/K, µm)", "gemischt", "L", "je_system", "Baustein 06")
_s("V_H", "v_hopf", "Quellspannung an der Hopf-Grenze des VO₂-Elements", "V", "L", "je_system", "Baustein 06")
_s("f_H", "f_hopf", "Frequenz des geborenen Zyklus an der Hopf-Grenze", "Hz", "L", "je_system", "Baustein 06")

# ------------------------------------------------------------------ Selbsterregter Umlauf mit thermischer Rückwirkung (05)
_s("g_E", "g_entdaempfung", "Entdämpfung (negativer Widerstand) eines selbsterregten Umlaufs", "Ω", "L", "je_system",
   "Baustein 05")
_s("b_I", "b_saettigung", "kubische Sättigung des Umlaufs (Glied b_I·I³)", "Ω/A²", "L", "je_system", "Baustein 05")
_s("dR/dT", "dR_dT", "Temperaturgang des Widerstands", "Ω/K", "L", "je_system", "Baustein 05")
