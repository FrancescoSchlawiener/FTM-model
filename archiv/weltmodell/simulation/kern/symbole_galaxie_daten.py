# -*- coding: utf-8 -*-
"""
Symbole der Gruppe Galaxie und Daten (Regeln in symbole.py; eingebunden dort mit einer Importzeile).
"""
from .symbole import SYMBOLE

REGISTERMODUL = True
_s = SYMBOLE.neu
SP = "SPARC (Lelli, McGaugh, Schombert 2016)"
AME = "AME2020 (Huang, Wang, Kondev, Audi, Naimi 2021)"
PUB = "PubChem, Periodensystem der Elemente"

# ------------------------------------------------------------------ Daten: SPARC
_s("V_obs", "V_obs", "gemessene Rotationsgeschwindigkeit", "m/s", "L", "je_system", SP)
_s("e_V", "e_V", "Messfehler von V_obs", "m/s", "L", "je_system", SP)
_s("V_gas", "V_gas", "Beitrag des Gases (mit Helium) zur Rotation", "m/s", "L", "je_system", SP)
_s("V_disk", "V_disk", "Beitrag der Sternscheibe zur Rotation bei Υ = 1", "m/s", "L", "je_system", SP)
_s("V_bul", "V_bul", "Beitrag des Bulges zur Rotation bei Υ = 1", "m/s", "L", "je_system", SP)
_s("V_f", "V_f", "flache Rotationsgeschwindigkeit, gemessen", "m/s", "L", "je_system", SP)
_s("Υ_d", "Upsilon_d", "Masse-zu-Licht der Scheibe bei 3,6 μm", "M☉/L☉[3,6]", "L", "je_system", SP)
_s("Υ_b", "Upsilon_b", "Masse-zu-Licht des Bulges bei 3,6 μm", "M☉/L☉[3,6]", "L", "je_system", SP)
_s("R_d", "R_d", "Skalenlänge der Scheibe", "m", "L", "je_system", SP)
_s("SB₀", "SB0", "zentrale Flächenhelligkeit der Scheibe bei 3,6 μm", "L☉[3,6]/pc²", "L", "je_system", SP)
_s("L_3,6", "L_36", "Leuchtkraft bei 3,6 μm", "L☉[3,6]", "L", "je_system", SP)
_s("M_HI", "M_HI", "Masse des neutralen Wasserstoffs (im Datensatz in M☉)", "kg", "L", "je_system", SP)
_s("i_inc", "i_inc", "Neigung (Scheibe gegen die Himmelsebene, Ring gegen die Scheibe)", "rad", "L", "je_system", SP)
_s("Q_fl", "Q_fl", "Güte einer SPARC-Kurve (1 hoch, 2 mittel, 3 niedrig)", "1", "L", "je_system", SP)

# ------------------------------------------------------------------ Daten: AME2020, PubChem
_s("Z_A", "Z_kern", "Kernladungszahl", "1", "L", "fest", AME)
_s("A_A", "A_kern", "Massenzahl", "1", "L", "fest", AME)
_s("E_B/A", "EB_je_A", "Bindungsenergie je Nukleon", "J", "L", "fest", AME)
_s("Δ_M", "mass_excess", "Massenexzess", "J", "L", "fest", AME)
_s("T_schm", "T_schmelz", "Schmelztemperatur bei Normaldruck", "K", "L", "fest", PUB)
_s("T_sied", "T_siede", "Siedetemperatur bei Normaldruck", "K", "L", "fest", PUB)
_s("E_ion", "E_ion", "erste Ionisierungsenergie", "J", "L", "fest", PUB)

# ------------------------------------------------------------------ Galaxie: Lesung an Rotationskurven (03, 07)
_s("V_bar", "V_bar", "baryonische Rotation, V_bar² = V_gas|V_gas| + Υ_d·V_disk² + Υ_b·V_bul²", "m/s", "L", "je_system",
   "Baustein 03")
_s("g_obs", "g_obs", "gemessene Zentripetalbeschleunigung (V_obs² + σ_v²·R/R_d)/R, mit dem Druck des Gases", "m/s²",
   "L", "je_system", "Baustein 03")
_s("ε_R", "eps_R", "Rest je Galaxie: Median von log₁₀(g_obs/(ν(g_N)·g_N)) in dex", "1", "L", "je_system", "Baustein 03")
_s("σ_ε", "sigma_eps", "Streuung einer Lesung aus den Messfehlern von Abstand, Neigung und Υ (dex)", "1", "L",
   "je_system", "Baustein 03")
_s("a₀_gal", "a0_gal", "Schwelle je Galaxie gelesen (Rest null)", "m/s²", "LI", "je_system", "Baustein 03", stamm="a₀")
_s("a₀_prof", "a0_prof", "Schwelle aus allen SPARC-Kurven zugleich, Abstand, Neigung und Υ je Galaxie als Störgrößen",
   "m/s²", "LI", "je_system", "Baustein 03", stamm="a₀")
_s("a₀_pkt", "a0_pkt", "Schwelle an den SPARC-Punkten, Bose-Form an Bin-Medianen", "m/s²", "LI", "je_system",
   "Baustein 07", stamm="a₀")
_s("Ȳ_bin", "y_bin", "Median einer Größe je Bin in log g_N (Rechenhilfe)", "1", "L", "je_system", "Baustein 07")
_s("F_q", "F_q", "Fehlerquadratsumme einer Anpassung an Bin-Mediane", "1", "L", "je_system", "Baustein 07")
_s("d_s", "d_s", "Spektraldimension der Rückkehr, 2·d ln a/d ln a_N", "1", "LI", "je_system", "Baustein 07")
_s("Σ", "Sigma", "Flächendichte einer Scheibe am Ort", "kg/m²", "L", "je_system", "Baustein 07")
_s("Σ₀", "Sigma0", "zentrale Flächendichte einer Scheibe (SPARC: Υ_d·SB₀)", "kg/m²", "L", "je_system", "Baustein 03")
_s("Σ_dyn", "Sigma_dyn", "zentrale dynamische Flächendichte (Toomre-Integral) (1/2πG)∫V²/R² dR", "kg/m²", "L",
   "je_system", "Baustein 03")

# ------------------------------------------------------------------ Galaxie: Scheibe, Feldform, Frequenzen (07)
_s("g_R", "g_R", "Newton-Feld in der Scheibenebene, radial (Betrag)", "m/s²", "L", "je_system", "Baustein 07")
_s("g_z", "g_z", "Newton-Feld senkrecht zur Ebene (Betrag; an der Scheibe 2πGΣ)", "m/s²", "L", "je_system",
   "Baustein 07")
_s("a_R", "a_R", "radiale Beschleunigung in der Ebene nach der Regel", "m/s²", "L", "je_system", "Baustein 07")
_s("f_Σ", "f_Sigma", "Gewicht des senkrechten Scheibenfelds in der Regel (Feldform 1, Bahnregel in der Ebene 0)", "1",
   "L", "je_system", "Baustein 07")
_s("b_K", "b_K", "Skala der Kuzmin-Scheibe", "m", "L", "je_system", "Baustein 07")
_s("a_MN", "a_MN", "radiale Skala der Miyamoto-Nagai-Scheibe", "m", "L", "je_system", "Baustein 07")
_s("b_MN", "b_MN", "vertikale Skala der Miyamoto-Nagai-Scheibe", "m", "L", "je_system", "Baustein 07")
_s("ω_p", "omega_p", "Knotenpräzession eines geneigten Rings", "s⁻¹", "L", "je_system", "Baustein 07")
_s("Ω_N", "Omega_N", "Kreisfrequenz nach Newton √(a_N/r), Punktmasse √(GM/r³), von der fernen Uhr gezählt", "s⁻¹", "L", "je_system", "Baustein 07", stamm="Ω")
_s("ω_z,N", "omega_zN", "Vertikalfrequenz nach Newton", "s⁻¹", "L", "je_system", "Baustein 07")
_s("ρ", "rho_ort", "Dichte der Baryonen am Ort", "kg/m³", "L", "je_system", "Baustein 07")
_s("ρ_X", "rho_X", "gelesene Dichte eines Beobachters, der Newton erwartet (Ausgabe seiner Rechnung, kein Stoff)",
   "kg/m³", "L", "je_system", "Baustein 07")

# ------------------------------------------------------------------ System als Ganzes im äußeren Feld (03)
_s("A", "A_mittel", "Massenmittel der Regel über das System ⟨ν(|g_i + g_e|)(g_i + g_e)⟩_M (Betrag, parallel zu g_e)",
   "m/s²", "L", "je_system", "Baustein 03")
_s("A_m1", "A_m1", "Amplitude der Einseitigkeit m = 1 der Rotation (Anteil am Mittel)", "1", "L", "je_system",
   "Baustein 03")
_s("R_e", "R_e", "Halblichtradius, Skala der Plummer-Kugel", "m", "L", "je_system", "Baustein 03")

# ------------------------------------------------------------------ Gezeit und Hierarchie der Trichter (02, 07)
_s("A_O", "A_Oort", "Oort-Konstante A = −½R·dΩ/dR", "s⁻¹", "L", "je_system", "Baustein 07")
_s("B_O", "B_Oort", "Oort-Konstante B = A − Ω", "s⁻¹", "L", "je_system", "Baustein 07")
_s("𝒯_H", "T_hill", "Gezeit des Hill-Radius, 4Ω² − κ² = 4A_O(A_O − B_O) (ungeglättet, andere Größe als 𝒯)", "s⁻²",
   "L", "je_system", "Baustein 02")

# ------------------------------------------------------------------ Filament, Faltungsmesser, Temperatur (03)
_s("μ_fil", "mu_fil", "Linienmasse eines Filaments", "kg/m", "L", "je_system", "Baustein 03")
_s("r*_fil", "r_stern_fil", "Schwelle eines Filaments 2Gμ/a₀", "m", "L", "je_system", "Baustein 03", stamm="r*")
_s("Δs_eff", "ds_eff", "gezählte Faltung ν·Δs (Faltung mal Rückkehrzahl; die Regel faltet nicht)", "m", "L",
   "je_system", "Baustein 03", stamm="Δs")
_s("M_b", "M_b", "baryonische Masse eines Systems, aus der Gastemperatur gelesen", "kg", "L", "je_system",
   "Baustein 03", stamm="M")
_s("S_T", "S_T", "logarithmisches Gefälle von Dichte und Temperatur, f = 1/S_T", "1", "L", "je_system", "Baustein 03")
_s("f_He", "f_helium", "Gasmasse gegen M_HI (Helium), M_gas = f_He·M_HI", "1", "L", "fest", "Setzung SPARC")
_s("σ_Υ", "sigma_ups", "Streuung von Υ_d je Galaxie im Fehlermodell (dex)", "1", "L", "fest", "Setzung SPARC")

# ------------------------------------------------------------------ Vorhersagen und Offenes in 03 (Einsetzen als Kernformel)
_s("Ṗ/P", "P_dot_rel", "relative Änderung der Pulsperiode durch die Beschleunigung längs der Sichtlinie", "s⁻¹", "L",
   "je_system", "Baustein 03")
_s("t_W", "t_warp", "Zeit, bis die Knotenlinie des Warps bei gemessener Rate um 90° dreht", "s", "L", "je_system",
   "Baustein 03")
_s("I_F", "I_F", "Übergabeintegral ∫ω dt bis zur ersten Faltung des Phasenraumblatts", "rad", "L", "je_system",
   "Baustein 03")
_s("t_T", "t_trenn", "Zeit, bis sich Linsenmaximum und Gas beim Haufenstoß trennen", "s", "L", "je_system", "Baustein 03")
_s("σ_v,N", "sigma_v_newton", "Geschwindigkeitsstreuung mit Newton allein (ohne Regel)", "m/s", "L", "je_system",
   "Baustein 03")
_s("V_σ", "verstaerkung_sigma", "Verstärkung der Streuung (σ_v/σ_v,N)² druckgestützter Systeme", "1", "L", "je_system",
   "Baustein 03")

# ------------------------------------------------------------------ Feldform in Scheiben (07)
_s("x_gl", "x_gleichstand", "Lage des Gleichstands von Feldform und radialer Regel in der Scheibenebene, in R_d", "1",
   "L", "je_system", "Baustein 07")
_s("ε_F", "eps_feldform", "gelöste Feldform gegen die örtliche Näherung ν(√(g_R² + (2πGΣ)²))·g_R in der Ebene", "1",
   "L", "je_system", "Baustein 07")
_s("a₀_gem", "a0_gem", "Schwelle gemeinsam je Drittel der Flächendichte angepasst (Υ_d fest)", "m/s²", "LI", "je_system",
   "Baustein 03", stamm="a₀")
_s("a₀_e", "a0_e", "Schwelle, die man ohne äußeres Feld liest, wo die Regel mit äußerem Feld wirkt", "m/s²", "LI",
   "je_system", "Baustein 03", stamm="a₀")
_s("k_a", "k_a0", "Anstieg von log₁₀ a₀ je Dekade der Flächendichte Σ₀", "1", "L", "je_system", "Baustein 03")
_s("ε_Σ", "eps_Sigma", "Rest der zentralen Flächendichte log₁₀(Σ_dyn,gem/Σ_dyn,Regel) (dex)", "1", "L", "je_system",
   "Baustein 03")
_s("k_Σ", "k_Sigma", "Anstieg des Rests ε_Σ je Dekade der Flächendichte Σ₀", "1", "L", "je_system", "Baustein 03")
_s("k_gas", "k_gas", "Änderung des Rests der Regel mit dem örtlichen Gasanteil innerhalb der Galaxien (dex)", "1", "L",
   "je_system", "Baustein 03")
_s("ε_Q", "eps_Q", "Rest der Regel an den Punkten jenseits R_Q (a_N < (ln q_kosm)²a₀), dex", "1", "L", "je_system",
   "Baustein 03")
_s("k_aus", "k_aussen", "Steigung des Rests gegen log a_N an den Außenpunkten (a_N < 0,1·a₀)", "1", "L", "je_system",
   "Baustein 03")
_s("σ_ds", "sigma_ds", "Fehler der Spektraldimension aus Ziehungen der SPARC-Kurven", "1", "LI", "je_system",
   "Baustein 07")
_s("B_a", "B_boden", "Boden der Regel a = ν(a_N)·a_N + B, in a₀", "1", "LI", "je_system", "Baustein 07")
_s("Pe", "Peclet", "Péclet-Zahl einer Drift der Rückkehr (Drift gegen Diffusion)", "1", "I", "je_system", "Baustein 07")
_s("b_q", "b_quer", "Bandbreite einer kontrahierenden Querrichtung der Rückkehr", "m", "I", "je_system", "Baustein 07")
_s("t_Wick", "t_aufwickel", "Zeit, bis die Knotenlinien zweier Ringe um Δφ auseinanderlaufen (Aufwicklung des Warps)", "s",
   "L", "je_system", "Baustein 03")
_s("𝒯_Q", "T_quadrupol", "Quadrupol des Restfelds der Regel um einen geglätteten Körper im Feld der Elternebene", "s⁻²",
   "L", "je_system", "Baustein 07")
_s("Σ_z", "Sigma_z", "gelesene Flächendichte der Säule |z| < z₁ am Ort, a_z(z₁)/(2πG) (Ausgabe einer Lesung, kein Stoff)",
   "kg/m²", "L", "je_system", "Baustein 07")
