# -*- coding: utf-8 -*-
"""
Symbole der Gruppe Teilchen, Sterne und Licht (02, 08, 09; Ebene L mit den Brücken der Längenleiter). Ergänzt das
Register in symbole.py (dort eine Importzeile); Regeln und Felder dort.
"""
from .symbole import SYMBOLE

REGISTERMODUL = True
_s = SYMBOLE.neu

# ------------------------------------------------------------------ Einheiten und Zählgrößen im Stern (02)
_s("n_P", "n_poly", "Index der Polytrope (1,5 konvektiv, 3 strahlend)", "1", "L", "je_system", "Baustein 02")
_s("x₁_LE", "x1_LE", "erste Nullstelle der Lane-Emden-Lösung θ(x₁) = 0", "1", "Z", "fest", "Baustein 02")
_s("m₁_LE", "m1_LE", "Massenzahl der Polytrope −x₁²θ′(x₁)", "1", "Z", "fest", "Baustein 02")
_s("ρ_c/ρ̄", "rho_c_mittel", "Zentraldichte gegen mittlere Dichte der Polytrope x₁³/(3m₁)", "1", "Z", "fest", "Baustein 02")
_s("f_T", "f_tiefe", "Anteil der Tiefe in kT_c = f·δ·m̄c², x₁/((n + 1)m₁)", "1", "Z", "fest", "Baustein 02")
_s("μ_m", "mu_mol", "mittlere Teilchenmasse in Nukleonmassen (m̄ = μ_m·m_u; im Gas des Halos μ_m·m_p)", "1", "L", "je_system", "Baustein 02")
_s("μ_e", "mu_e", "Masse je Elektron in m_u", "1", "L", "je_system", "Baustein 02")
_s("T_c", "T_kern", "Kerntemperatur eines nicht entarteten Trichters, am Ort gemessen", "K", "L", "je_system",
   "Baustein 02")
_s("M_Grenze", "M_grenze", "Grenzmasse der Zählregel der Elektronen (Tiefe gegen hochrelativistische Zählregel)", "kg",
   "L", "fest", "Baustein 02")
_s("T_max", "T_max", "höchste Kerntemperatur eines Trichters, den die Zählregel der Elektronen hält", "K", "L",
   "je_system", "Baustein 02")
_s("g_T", "g_decke", "Decke kT_max/((μ_m/μ_e)m_ec²) an der Grenzmasse, max_x[x/4 − Π(x)]", "1", "Z", "fest",
   "Baustein 02")
_s("𝒞", "C_carter", "Carter-Faktor α_G/(α¹²(m_e/m_p)⁴)", "1", "Z", "fest", "Baustein 02")

# ------------------------------------------------------------------ Leiter der Längen und Schnittmassen (01, 08)
_s("r*/R", "r_stern_zu_R", "Lage der Schwelle zur Hülle √(g(R)/a₀); trennt kompakte und kollektive Klasse", "1", "L",
   "je_system", "Baustein 08")
_s("M̃", "M_spiegel", "Spiegelmasse m_Pl²/M (vertauscht λ̄ und 3Δs)", "kg", "L", "je_system", "Baustein 08")
_s("r_sp", "r_spiegel", "Spiegelort r*²/r (bildet r₀ auf Λ ab, den Raum auf die innere Richtung)", "m", "LI",
   "je_system", "Autor (FTM-PDF)")

# ------------------------------------------------------------------ Raten und Kopplungen (09)
_s("α/α_G", "alpha_zu_alphaG", "elektrische gegen gravitative Wechselwirkung zweier Teilchen", "1", "L", "je_system",
   "Baustein 09")
_s("J", "J_dreh", "Drehimpuls eines Körpers", "J·s", "L", "je_system", "Baustein 08")
_s("μ_mag", "mu_mag", "magnetisches Moment Q_el·c·a_J", "J/T", "L", "je_system", "Baustein 09")
_s("Γ_G", "Gamma_G", "Faltungsrate α_G(M)·ω_C", "s⁻¹", "L", "je_system", "Baustein 09")
_s("ω_F", "omega_F", "Faltungstakt c/(3Δs) (Faltung im Zentrum, keine Raumlänge)", "s⁻¹", "L", "je_system",
   "Baustein 09")
_s("ω_H", "omega_H", "Umlauf des Elektrons im Wasserstoff αc/a_B", "s⁻¹", "L", "fest", "Baustein 09")
_s("Δω", "d_omega", "Frequenzunterschied zweier Zyklen", "s⁻¹", "Z", "je_system", "Baustein 09")
_s("K_E", "K_einrast", "Kopplung zweier Zyklen beim Einrasten", "s⁻¹", "Z", "je_system", "Baustein 09")
_s("Δθ", "d_theta", "Phasenlage zweier eingerasteter Zyklen", "rad", "Z", "je_system", "Baustein 09")
_s("t_φ", "t_dephas", "thermische Dephasierungszeit ħ/(k_BT)", "s", "L", "je_system", "Baustein 09")

# ------------------------------------------------------------------ Licht und Starkfeld (08)
_s("r_ph", "r_photon", "Photonensphäre 9Δs", "m", "L", "je_system", "Baustein 08")
_s("r̃", "r_tilde", "Bahnradius in Einheiten 3Δs", "1", "L", "je_system", "Baustein 08")
_s("ω_φ", "omega_phi", "azimutale Grundfrequenz einer Kerr-Bahn", "s⁻¹", "L", "je_system", "Baustein 08")
_s("ω_r", "omega_r", "radiale Grundfrequenz einer Kerr-Bahn", "s⁻¹", "L", "je_system", "Baustein 08")
_s("ω_θ", "omega_theta", "vertikale Grundfrequenz einer Kerr-Bahn", "s⁻¹", "L", "je_system", "Baustein 08")
_s("ζ", "zeta_knoten", "Knotendrehung der Bahnebene (ω_φ − ω_θ)/ω_φ", "1", "L", "je_system", "Baustein 08")
_s("ϑ", "theta_pol", "Polwinkel gegen die Drehachse", "rad", "L", "je_system", "Baustein 08")
_s("r_E", "r_ergo", "Nullfläche des Takts ruhender Uhren 3Δs + √(9Δs² − a_J²cos²ϑ)", "m", "L", "je_system",
   "Baustein 08")
_s("r_±", "r_horizont", "Horizonte 3Δs ± √(9Δs² − a_J²)", "m", "L", "je_system", "Baustein 08")
_s("ΔT_pr", "dT_prograd", "Umlaufzeit mitläufig gegen gegenläufig 4π·a_J/c", "s", "L", "je_system", "Baustein 08")
_s("cos ψ", "cos_psi", "Kosinus des Winkels zwischen Drehachse und Ort", "1", "L", "je_system", "Baustein 08")
_s("Ω_LT", "Omega_mitf", "Mitführung eines Kreisels (Frame-Dragging), Betrag", "s⁻¹", "L", "je_system", "Baustein 08")
_s("δ_L", "delta_leitstern", "Winkel des Leitsterns zur Bahnebene (Deklination bei polarer Bahn)", "rad", "L", "je_system", "Baustein 08")
_s("Ω_LS", "Omega_mitf_leitstern", "Mitführung vom Kreisel zum Leitstern gesehen, Ω_LT·cos δ_L", "s⁻¹", "L", "je_system", "Baustein 08")
_s("Ω_geo", "Omega_geod", "geodätische Drehung eines Kreisels auf Kreisbahn", "s⁻¹", "L", "je_system", "Baustein 08")
_s("P_GW", "P_gw", "Quadrupol-Leistung einer Kreisbahn", "W", "L", "je_system", "Baustein 08")
_s("f_gw", "f_gw", "Frequenz einer Gravitationswelle im Detektor", "Hz", "L", "je_system", "Baustein 08")
_s("ḟ_gw", "fdot_gw", "Anstieg der Frequenz (Chirp)", "Hz/s", "L", "je_system", "Baustein 08")
_s("Δs_c", "ds_chirp", "Faltung der Chirp-Masse G·M_c/(3c²)", "m", "L", "je_system", "Baustein 08", stamm="Δs")
_s("F_j", "F_ringdown", "Eigenfrequenz des Trichters in Einheiten c/(2π·3Δs), hängt nur am Endspin", "1", "L",
   "je_system", "Baustein 08")
_s("f_RD", "f_ringdown", "Ringdown-Frequenz c/(2π·3Δs)·F(j)", "Hz", "L", "je_system", "Baustein 08")
_s("E_rad", "E_rad", "abgestrahlte Energie einer Quelle", "J", "L", "je_system", "Baustein 08")
_s("Δs_rad", "ds_rad", "Faltung der abgestrahlten Energie G·E_rad/(3c⁴)", "m", "L", "je_system", "Baustein 08",
   stamm="Δs")
_s("D_L", "D_L", "Leuchtkraftabstand einer Quelle", "m", "G", "epoche", "Baustein 08")
_s("κ_Θ", "kappa_richtung", "Richtungsfaktor des Gedächtnis-Effekts", "1", "L", "je_system", "Baustein 08")
_s("h_mem", "h_mem", "bleibende Verschiebung nach einer Gravitationswelle", "1", "L", "je_system", "Baustein 08")
_s("Φ_X", "Phi_X", "Zusatzpotential des kollektiven Trichters v_flach²·ln ν(s)", "m²/s²", "LI", "je_system",
   "Baustein 08")
_s("x_F", "x_fermi", "Fermi-Impuls der Elektronen in m_ec, λ̄_e·(3π²n_e)^{1/3}", "1", "L", "je_system", "Baustein 02")
_s("Π_e", "Pi_e", "Druck der Zählregel je Elektron in m_ec², P_e/(n_e m_ec²)", "1", "L", "je_system", "Baustein 02")
_s("∂P/∂T", "dP_dT", "Änderung der eingespeisten Heizleistung mit der Temperatur (Schleife Feld → Strom → Wärme → Widerstand)",
   "W/K", "L", "je_system", "Baustein 09")

# ------------------------------------------------------------------ Starke Kopplung und Ladung (02)
_s("α_s", "alpha_s", "starke Kopplung bei der Skala μ_R (läuft mit der Skala)", "1", "L", "je_system", "Baustein 02")
_s("μ_R", "mu_skala", "Energieskala der laufenden Kopplung", "GeV", "L", "je_system", "Baustein 02")
_s("M_Z", "M_Z", "Masse des Z-Bosons (Bezugsskala der starken Kopplung)", "GeV", "L", "fest", "Baustein 02")
_s("m_q", "m_quark", "Quarkmasse, Schwelle des Laufs (m_c, m_b, m_t)", "GeV", "L", "fest", "Baustein 02")
_s("n_f", "n_flavour", "Zahl der leichten Quarksorten unter der Skala", "1", "L", "je_system", "Baustein 02")
_s("q_QCD", "q_qcd", "Nome der starken Kopplung Λ_QCD/M = e^{−2π/(β₀α_s)} (im Modelltext |q̃|)", "1", "L", "je_system",
   "Baustein 02", stamm="q")
_s("g_D/e", "dirac_verhaeltnis", "magnetische gegen elektrische Elementarladung 1/(2α) (Dirac-Bedingung)", "1", "L",
   "fest", "Baustein 02")
_s("r_Qel", "r_ladungskreis", "Radius des Ladungskreises 2ℓ_P/√α", "m", "L", "fest", "Baustein 02")
_s("E_Qel", "E_ladungskreis", "Energie eines Umlaufs um den Ladungskreis ħc/r_Qel", "J", "L", "fest", "Baustein 02")

# ------------------------------------------------------------------ Zustand des Stoffs und Sterne (02)
_s("m_A", "m_atom", "Masse eines Atoms", "kg", "L", "je_system", "Baustein 02")
_s("δ_Z", "delta_zustand", "Zustandsschwelle 2k_BT/(m_Ac²) = 2(v_th/c)²", "1", "L", "je_system", "Baustein 02")
_s("f_s", "f_saha", "Oberflächentemperatur in Einheiten α²m_ec²: k_BT_eff/(α²m_ec²)", "1", "L", "je_system",
   "Baustein 02")
_s("β_E", "beta_gas", "Anteil des Gasdrucks am Druck (Eddington-Standardmodell, Polytrope n = 3)", "1", "L", "je_system",
   "Baustein 02")
_s("M_β", "M_eddington", "Masse eines Sterns mit Gasdruckanteil β_E im Eddington-Standardmodell", "kg", "L", "je_system",
   "Baustein 02")

# ------------------------------------------------------------------ Syntropie und Falte (04)
_s("S_F", "S_falte", "Entropie der Falte 36π(Δs/ℓ_P)² in Einheiten von k_B", "1", "L", "je_system", "Baustein 04")
_s("ΣS_F", "S_falten_summe", "Summe der Entropien aller Falten im Horizont, in Einheiten von k_B", "1", "G", "epoche",
   "Baustein 04")
_s("t_P/t_F", "page_anteil", "Page-Zeit als Anteil der Lebensdauer der Falte", "1", "Z", "fest", "Baustein 04")
_s("z_f", "z_faktor", "Rotverschiebungsfaktor a₂/a₁ der Strahlung bis zur Rückverdichtung", "1", "G", "je_system",
   "Baustein 04")
_s("N_RV", "n_rueckverdichtung", "Mindestzahl der Falten, deren Strahlung zu einer neuen Falte zurückverdichtet", "1", "Z",
   "je_system", "Baustein 04")
_s("P_abs", "P_absorbiert", "absorbierte Strahlungsleistung eines Körpers", "W", "L", "je_system", "Baustein 04")
_s("A_alb", "albedo", "Albedo (zurückgeworfener Anteil der einfallenden Strahlung)", "1", "L", "je_system", "Baustein 04")
_s("T_A", "T_abstrahlung", "Abstrahlungstemperatur eines Körpers im Strahlungsgleichgewicht", "K", "L", "je_system",
   "Baustein 04")
_s("J_S", "J_entropie", "Entropieexport eines Körpers im Strahlungsstrom", "W/K", "L", "je_system", "Baustein 04")

# ------------------------------------------------------------------ Kerne (06)
_s("a_Vol", "a_volumen", "Volumenkoeffizient des Tröpfchenmodells", "J", "L", "fest", "Baustein 06")
_s("a_Obf", "a_oberflaeche", "Oberflächenkoeffizient des Tröpfchenmodells (E_S = a_Obf·A^{2/3})", "J", "L", "fest",
   "Baustein 06")
_s("a_Cou", "a_coulomb", "Coulomb-Koeffizient des Tröpfchenmodells (E_C = a_Cou·Z²/A^{1/3})", "J", "L", "fest",
   "Baustein 06")
_s("a_Asy", "a_asymmetrie", "Asymmetriekoeffizient des Tröpfchenmodells", "J", "L", "fest", "Baustein 06")
_s("a_Paar", "a_paarung", "Paarungskoeffizient des Tröpfchenmodells", "J", "L", "fest", "Baustein 06")
_s("x_Sp", "x_spalt", "Spaltbarkeit E_C/(2E_S) = a_Cou·Z²/(2a_Obf·A); Kippstelle bei 1", "1", "L", "je_system",
   "Baustein 06")

# ------------------------------------------------------------------ Wellen (08)
_s("Q_RD", "Q_ringdown", "Güte des Ringdowns f_RD/(2·Dämpfungsrate), hängt nur am Endspin", "1", "L", "je_system",
   "Baustein 08")
_s("τ_RD", "tau_ringdown", "Abklingzeit des Ringdowns Q_RD/(πf_RD)", "s", "L", "je_system", "Baustein 08")
_s("ι", "iota_quelle", "Neigung der Bahnachse einer Wellenquelle gegen die Sichtlinie (Kante: π/2)", "rad", "L",
   "je_system", "Baustein 08")

# ------------------------------------------------------------------ Kohärenz (09)
_s("ε_r", "eps_rel", "relative Permittivität eines Mediums", "1", "L", "je_system", "Baustein 09")
_s("I_ion", "I_ion", "Ionenstärke einer Lösung ½Σc_iz_i²", "mol/m³", "L", "je_system", "Baustein 09")
_s("λ_D", "lambda_debye", "Debye-Länge √(ε₀ε_rk_BT/(2N_Ae²I_ion))", "m", "L", "je_system", "Baustein 09")
