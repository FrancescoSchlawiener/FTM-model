# -*- coding: utf-8 -*-
"""
Symbole der Module kosmos und naht_haufen, die das Grundregister nicht führt (Regeln und Klasse in symbole.py;
eingebunden am Ende von symbole.py). Dieselbe Sache auf einer anderen Ebene trägt das Grundzeichen mit Index
(Δs_H zu Δs, H_eff zu H, z_c zu z).
"""
from .symbole import SYMBOLE

REGISTERMODUL = True
_s = SYMBOLE.neu

# ------------------------------------------------------------------ Ebene G: Horizonte, Hintergrund, Vakuumteil (03, 04)
_s("a_GH", "a_GH", "Temperatur des Hubble-Horizonts als Beschleunigung cH/2π (heute nahe a₀, Θ = a₀/a_GH)", "m/s²", "G",
   "epoche", "Baustein 03")
_s("d_H", "d_H", "Teilchenhorizont eines Beobachters, mitbewegt", "m", "G", "epoche", "Baustein 03")
_s("ΔL/R", "dL_R", "Dehnung der Raumschleifen je Radius = ΣΔs_b/R im eigenen Teilchenhorizont", "1", "G", "epoche",
   "Baustein 03")
_s("Δs_H", "ds_H", "Faltung des eigenen Horizonts r_H/6 = c/(6H) (globale Engstelle)", "m", "G", "epoche", "Baustein 04",
   stamm="Δs")
_s("S_Λ", "S_Lam", "Horizontentropie des Vakuumteils allein (45/π³)(Λ/ℓ_P)² in Einheiten von k_B", "1", "IG", "seit_z10",
   "Baustein 03")
_s("d ln a₀/d ln a", "dlna0_dlna", "Gang der Schwelle mit dem Skalenfaktor (heute null)", "1", "IG", "epoche",
   "Baustein 03")
_s("L·H/c", "L_H_c", "Zelle in Hubble-Längen der Epoche", "1", "G", "epoche", "Baustein 03")
_s("f³", "f3_zelle", "Formfaktor der Zelle ohne Kreispaare (Würfel 1, sechseckig 1,15, raumzentriert 1,30, "
   "flächenzentriert 1,41)", "1", "Z", "fest", "Baustein 03")
_s("(1+Q)a₀", "a0_nu_grenze", "Grenze von (1 + Q)·a₀ aus der Zelle ohne Kreispaare, f³·3c⁴/(2Gρ_bL_min³)", "m/s²", "G",
   "epoche", "Baustein 03")
_s("H_b", "H_b", "Atemrate aus Bremsglied der Materie und Vakuumteil π⁴/(15Λ²), ohne Fluchtmaß", "s⁻¹", "G", "epoche",
   "Baustein 04")
_s("Ẋ_an", "X_an_dot", "Rate des über die Rückkehr angekommenen Fluchtmaßes X_an", "s⁻¹", "ZG", "epoche", "Baustein 03")
_s("H_eff", "H_eff", "Atemrate mit angekommenem Fluchtmaß H_b + ⅓Ẋ_an", "s⁻¹", "G", "epoche", "Baustein 03", stamm="H")
_s("ln a", "ln_a", "Atemmaß des Ganzen ln(L(t)/L(t₀)) = −ln(1 + z): die Epoche, auch in der Zukunft (z → −1)", "1", "G",
   "epoche", "Baustein 03")
_s("t_A", "t_welt", "Weltzeit längs der Atemrate bis zur Epoche, ∫d ln a/H", "s", "G", "epoche", "Baustein 04")
_s("Ṡ_H", "S_H_dot", "Wachstum der Horizontentropie −2S_H·Ḣ/H in k_B/s", "s⁻¹", "G", "epoche", "Baustein 03")
_s("Ω_X", "Omega_X", "fester Rest der Atemrate im Vergleichshintergrund (Vakuumteil Ω_s, oder im Standardmodell "
   "gelesen 1 − Ω_m − Ω_r)", "1", "G", "epoche", "Baustein 03")
_s("X", "X_flucht", "inaktives Fluchtmaß je Baryon: ln-Maß der ganzen Tori im größten gebundenen Bereich, Regel gegen "
   "Newton (Zustandsraum, kein Stoff)", "1", "ZG", "epoche", "Baustein 03")
_s("X_an", "X_an", "über die Rückkehr angekommener Teil des Fluchtmaßes (Atmen der Raumschleifen δ ln L)", "1", "ZG",
   "epoche", "Baustein 03")
_s("t_F/t_ta", "tF_zu_tta", "Übergabe bei der ersten Faltung in Einheiten der Umkehrzeit", "1", "LG", "fest",
   "Baustein 03")
_s("σ_M", "sigma_M", "lineare Streuung der Dichte heute auf der Masse M (Kugelfilter)", "1", "G", "epoche",
   "Baustein 03")
_s("D_+", "D_plus", "lineares Wachstum der Störungen D(a)", "1", "G", "epoche", "Baustein 03")
_s("δ_ta", "delta_ta", "lineare Schwelle der Umkehr: auf die Umkehrepoche fortgesetzter linearer Anfangskontrast",
   "1", "LG", "epoche", "Baustein 03")
_s("A_s", "A_s", "Amplitude des Anfangsspektrums", "1", "G", "fest", "Messwert Planck 2018 [G]")
_s("n_s", "n_s", "Neigung des Anfangsspektrums", "1", "G", "fest", "Messwert Planck 2018 [G]")
_s("τ_re", "tau_re", "optische Tiefe der Reionisation", "1", "G", "fest", "Messwert Planck 2018 [G]")
_s("Σm_ν", "m_nu", "Summe der Neutrinomassen im Spektrum", "eV", "G", "fest", "Messwert [G]")
_s("Δ(fσ₈)/fσ₈", "d_fsigma8", "relative Änderung des linearen Wachstums fσ₈ durch die Schließungsregel gegen dasselbe "
   "Wachstum ohne Übergabe", "1", "G", "epoche", "Baustein 03")
_s("Φ/c²", "phi_c2", "Potentialtiefe einer linearen Störung in Einheiten von c²", "1", "G", "epoche", "Baustein 03")

# ------------------------------------------------------------------ Ebene Z und G: Form der Zelle, Himmel des Beobachters (03)
_s("Gitter", "gitter_zelle", "Gitterform der Zelle (Würfel, sechseckig, raumzentriert, flächenzentriert)", "1", "Z",
   "fest", "Baustein 03")
_s("φ_D", "phi_dreh", "Drehwinkel einer Drehsymmetrie der Zelle (flacher Torus)", "°", "Z", "fest", "Baustein 03")
_s("ℓ_W", "ell_winkel", "Ordnung einer Winkelmode am Himmel (2 Quadrupol, 3 Oktupol)", "1", "G", "fest", "Baustein 03")
_s("Anteil_ℓ", "anteil_ell", "Anteil der Leistung einer Winkelordnung in einer Auswahl ihrer Formen", "1", "G", "fest",
   "Baustein 03")
_s("C_ℓ/C_ℓ∞", "C_ell_rel", "Leistung einer Winkelordnung in der geschlossenen Zelle gegen den unendlichen Raum, "
   "für den Beobachter heute", "1", "G", "epoche", "Baustein 03")
_s("σ_ℓ", "sigma_ell", "kosmische Streuung: relative Streuung der Leistung einer Winkelordnung an einem Himmel",
   "1", "G", "fest", "Baustein 03")

# ------------------------------------------------------------------ Brücke L–G und Z: Naht (03)
_s("z_c", "z_c", "Rotverschiebung, bei der ein Bereich kollabiert", "1", "LG", "je_system", "Baustein 03", stamm="z")
_s("t_w", "t_w", "Dauer, in der w vor der Umkehr einer Schale von 0,9 auf 0,1 fällt", "s", "LG", "je_system", "Baustein 03")
_s("X_b", "X_b", "Fluchtmaß je Bereich: ln des Maßes der ganzen Tori in der Naht, Regel gegen Newton am selben Rand",
   "1", "LZ", "je_system", "Baustein 03")
_s("R/R_ta", "R_zu_Rta", "Radius einer Schale in Einheiten des Umkehrradius", "1", "LG", "je_system", "Baustein 03")
_s("I_vor", "I_vor", "Übergabeintegral ∫(1 − w)ω dt einer heute umkehrenden Schale bis zur Umkehr", "rad", "LG", "je_system",
   "Baustein 03")
_s("e_M", "e_profil", "Anfangsprofil eines Bereichs δM/M ∝ M^{−e} (Keim e = 1)", "1", "L", "je_system", "Baustein 03")
_s("ε_j", "eps_j", "Drehimpuls der Schale ab der Umkehr, j² = ε_j·GM·r_ta", "1", "L", "je_system", "Baustein 03")
_s("ΔΣ", "delta_Sigma", "Linsensignal Σ̄(<R) − Σ(R) im Abstand R (gelesene Faltung, projiziert)", "kg/m²", "L",
   "je_system", "Baustein 03")
_s("M_b(<r)", "M_b_innen", "eingeschlossene Baryonenmasse eines Aufbaus (Gas, Zentralgalaxie, Galaxien)", "kg", "L",
   "je_system", "Baustein 03")
_s("M_L(<r)", "M_linse_innen", "gelesene Linsenmasse innerhalb r nach einem Referenzprofil (gelesen, kein Stoff)", "kg",
   "L", "je_system", "Baustein 03")
_s("F_H", "F_haufen", "Rest im Haufen: gelesene Linsenmasse durch die Masse nach der örtlichen Regel, M_L/(ν·M_b)", "1",
   "L", "je_system", "Baustein 03")
_s("c_NFW", "c_nfw", "Konzentration des Referenzprofils der Linse (r₂₀₀/r_s)", "1", "L", "je_system", "Baustein 03")
