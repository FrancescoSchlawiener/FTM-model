# -*- coding: utf-8 -*-
"""
Eingaben der Gruppe Kosmos (Zelle, Weltzeit, Linse jenseits der Naht, Rest im Haufen; Klasse und Regeln in
eingabe.py). Gewählte Rechenpunkte und Aufbauten, kein [G].
"""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D = "03-von-der-galaxie-bis-zum-kosmos.md"

# ------------------------------------------------------------------ Störungen bei der Rekombination
_e("PHI_STOERUNG", "Φ/c²", 1e-5, "1", "Potentialtiefe linearer Störungen, untere Lesung der Spanne a_N ≈ 2,6–12·a₀(z) "
   "(Zahl nicht im Text)")

# ------------------------------------------------------------------ Linse der Milchstraße jenseits der Naht
_e("R_LINSE_1MPC", "R", 1.0, "Mpc", "Abstand des Linsensignals der Milchstraße jenseits 1,05·R_ta",
   (_D, "Milchstraße bei 1 Mpc"))

# ------------------------------------------------------------------ Haufen: Aufbau der Baryonen (Notiz 28) und Rechenpunkt
_e("HAUFEN_GAS_M", "M", 1e14, "M☉", "Gas des Haufens innerhalb HAUFEN_GAS_RNORM (β-Profil, Aufbau Notiz 28)")
_e("HAUFEN_GAS_RC", "r", 150.0, "kpc", "Kernradius des β-Gases (Aufbau Notiz 28)")
_e("HAUFEN_GAS_RNORM", "r", 1.0, "Mpc", "Radius, innerhalb dessen das Gas HAUFEN_GAS_M trägt (Aufbau Notiz 28)")
_e("HAUFEN_GAS_RMAX", "r", 3.0, "Mpc", "äußerer Rand des Gases (Aufbau Notiz 28)")
_e("HAUFEN_BCG_M", "M", 1.5e12, "M☉", "Zentralgalaxie des Haufens (Hernquist, Aufbau Notiz 28)")
_e("HAUFEN_BCG_B", "b_H", 20.0, "kpc", "Hernquist-Skala der Zentralgalaxie (Aufbau Notiz 28)")
_e("HAUFEN_GAL_M", "M", 1e13, "M☉", "Galaxien des Haufens (Hernquist, Aufbau Notiz 28)")
_e("HAUFEN_GAL_B", "b_H", 300.0, "kpc", "Hernquist-Skala der Galaxien des Haufens (Aufbau Notiz 28)")
_e("NSTERN_EINS", "N*", 1.0, "rad", "Rechenpunkt des linearen Wachstums mit Übergabe innerhalb der Spanne N* ≈ 0,3–2,4 "
   "(Rechnung gewicht_rueckkehr; Zahl nicht im Satz)")

# ------------------------------------------------------------------ Schale außerhalb der Naht
_e("R_SCHALE_AUSSEN", "R/R_ta", 1.33, "1", "Rechenpunkt außerhalb der Naht, an dem w gelesen wird", (_D, "0,95 ab 1,33·R_ta"))

# ------------------------------------------------------------------ Übergabe bei der ersten Faltung (selbstähnlicher Einfall)
_e("IF_E_KEIM", "e_M", 1.0, "1", "Keim: δM/M ∝ M^{−1} (δ̄ ∝ r^{−3}), Anfangsprofil der Spanne 1,1–1,25 (Rechnung "
   "uebergabe_faltung)")
_e("IF_E_FLACH", "e_M", 0.2, "1", "flachstes gerechnetes Anfangsprofil δ̄ ∝ r^{−0,6} (Rechnung uebergabe_faltung)")
_e("IF_EJ_KLEIN", "ε_j", 0.01, "1", "kleinster gerechneter Drehimpuls ab der Umkehr (Rechnung uebergabe_faltung)")
_e("IF_EJ_GROSS", "ε_j", 0.2, "1", "größter gerechneter Drehimpuls ab der Umkehr (Rechnung uebergabe_faltung)")

# ------------------------------------------------------------------ Lauf des eigenen Atmens (forschung/kosmos.py)
_e("T_F_ZU_T_TA", "t_F/t_ta", 1.8, "1", "Übergabe bei der ersten Faltung nach der Umkehr, t_F/t_ta = 1,71–1,89 "
   "(Rechnung uebergabe_faltung); Mitte als Rechenpunkt")
_e("FLUCHT_M_MIN", "M", 1e6, "M☉", "kleinste gezählte Masse M_ges der gebundenen Bereiche (Zählung ohne untere Schranke "
   "des Modells; Gitter der Massen ab 10⁶ M☉)")
_e("R_HAUFEN_REST_1MPC", "r", 1.0, "Mpc", "äußerer Rechenpunkt des Rests im Haufen", (_D, "1,6 bei 1 Mpc"))
