# -*- coding: utf-8 -*-
"""Eingaben zu 03, Teil a (Galaxie, Dunkle Materie, Schließungsregel; Klasse und Regeln in eingabe.py)."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D = "03-von-der-galaxie-bis-zum-kosmos.md"

# ------------------------------------------------------------------ Tabelle am Kopf von 03
_e("R_MILCHSTRASSE_RAND", "R", 15.3, "kpc", "Rand der Milchstraße, Bezug von Δs/R in der Tabelle (Zahl nicht im Text)")
_e("R_HAUFEN", "R", 2.0, "Mpc", "Radius des Haufens, Bezug von Δs/R in der Tabelle (Zahl nicht im Text)")

# ------------------------------------------------------------------ Filament
_e("MU_FILAMENT", "μ_fil", 1e12, "M☉/Mpc", "Linienmasse des Beispielfilaments", (_D, "10¹² M☉/Mpc"))

# ------------------------------------------------------------------ Milchstraße als Punktmasse und als Scheibe
_e("R_5KPC", "r", 5.0, "kpc", "Punkt der Rotationskurve", (_D, "bei 5 kpc"))
_e("R_15KPC", "r", 15.0, "kpc", "Punkt der Rotationskurve", (_D, "bei 15 kpc"))
_e("R_50KPC", "r", 50.0, "kpc", "Punkt der Rotationskurve, des Außenfelds und der Einseitigkeit", (_D, "bei 50 kpc"))
_e("R_10KPC", "r", 10.0, "kpc", "innerster Punkt der Einseitigkeit der Scheibe (Zahl nicht im Text)")
_e("R_20KPC", "r", 20.0, "kpc", "Punkt der Einseitigkeit", (_D, "20–50 kpc"))
_e("R_30KPC", "r", 30.0, "kpc", "Punkt der Einseitigkeit", (_D, "bis 30 kpc"))
_e("RD_MILCHSTRASSE", "R_d", 3.0, "kpc", "Skalenlänge der Milchstraße als Exponentialscheibe", (_D, "R_d = 3 kpc"))

# ------------------------------------------------------------------ druckgestützte Kugeln, Einzeltrichter oder kollektiv
_e("M_KUGEL_KLEIN", "M", 1e7, "M☉", "druckgestützte Kugel im Tiefbereich", (_D, "bei 10⁷ M☉"))
_e("M_KUGEL_GROSS", "M", 2e8, "M☉", "druckgestützte Kugel im Tiefbereich", (_D, "bei 2·10⁸ M☉"))
_e("G_KUGELSTERNHAUFEN", "g", 25.0, "a₀", "Feld eines kompakten Kugelsternhaufens", (_D, "g/a₀ ≈ 25"))
_e("G_ZWERG", "g", 0.13, "a₀", "Feld sehr dünner Kugelsternhaufen und Zwerggalaxien", (_D, "g/a₀ ≈ 0,13"))

# ------------------------------------------------------------------ äußeres Feld
_e("GE_002", "a_e", 0.02, "a₀", "gemeinsames äußeres Feld senkrecht zur Scheibe, untere Grenze", (_D, "g_e = 0,02–0,03·a₀"))
_e("GE_003", "a_e", 0.03, "a₀", "äußeres Feld (SPARC, Milchstraße bei 50 kpc, Scheibe, Außenpunkte)",
   ((_D, "im Feld 0,03a₀"), (_D, "ein Feld von 0,03·a₀")))
_e("GE_01", "a_e", 0.1, "a₀", "äußeres Feld der Plummer-Kugel", (_D, "bei g_e = 0,1·a₀"))
_e("R_E_PLUMMER", "R_e", 1.0, "kpc", "Skala der Plummer-Kugel; A/g_e hängt nur an GM/b² (Zahl nicht im Text)")
_e("SIG_GAS_8", "σ_v", 8.0, "km/s", "Geschwindigkeitsstreuung des Gases (Druck) in g_obs der SPARC-Punkte, untere Grenze; "
   "dieselbe Eingabe in 07 (Spektraldimension)", (_D, "8–10 km/s"))
_e("R_100KPC", "r", 100.0, "kpc", "Radius im Haufenkern (Rest der Regel gegen die Linse; Rechnung fehlt im Kern)",
   (_D, "bei 100 kpc"))
_e("R_80KPC", "r", 80.0, "kpc", "Radius der Außenbahn einer Galaxie im Haufen", (_D, "bei 80 kpc"))
_e("HERNQUIST_A", "b_H", 3.0, "kpc", "Hernquist-Skala der Galaxie im Haufen (Zahl nicht im Text)")
_e("M_HAUFEN_UMGEBUNG", "M", 1e14, "M☉", "Masse des Haufens innerhalb 1 Mpc, Aufbau der Umgebung (Zahl nicht im Text)")
_e("HAUFEN_R_ISO", "r", 1.0, "Mpc", "isothermer Haufen: M(<r) ∝ r bis 1 Mpc, außen fest (Zahl nicht im Text)")
_e("D_HAUFEN_NAH", "r", 0.5, "Mpc", "Abstand der Galaxie vom Haufenzentrum, nah (Zahl nicht im Text)")
_e("D_HAUFEN_FERN", "r", 2.0, "Mpc", "Abstand der Galaxie vom Haufenzentrum, fern (Zahl nicht im Text)")

# ------------------------------------------------------------------ Schwelle
_e("H0_70", "H₀", 70.0, "km/s/Mpc", "obere Grenze der Atemrate heute für cH₀/2π", (_D, "H₀ = 67–70 km/s/Mpc"))
_e("THETA_UNTEN", "Θ", 0.9, "1", "T_I/T_GH innerhalb von 10 %: untere Grenze 1 − 0,1", (_D, "innerhalb von 10 %"))
_e("UPS_D_04", "Υ_d", 0.4, "M☉/L☉[3,6]", "Masse-zu-Licht der Scheibe, untere Lesung", (_D, "Υ_d = 0,4–0,6"))
_e("UPS_D_06", "Υ_d", 0.6, "M☉/L☉[3,6]", "Masse-zu-Licht der Scheibe, obere Lesung", (_D, "Υ_d = 0,4–0,6"))

# ------------------------------------------------------------------ Schließungsregel
_e("W_OBEN", "w", 0.9, "1", "Beginn des Übergangs von w vor der Umkehr und äußerer Rand der w-Schale (w = 0,9)", (_D, "von 0,9 auf 0,1"))
_e("W_UNTEN", "w", 0.1, "1", "Ende des Übergangs von w vor der Umkehr", (_D, "von 0,9 auf 0,1"))
_e("Z_C_34", "z_c", 3.4, "1", "Kollaps bei z = 3,4", (_D, "bei Kollaps bei z = 3,4"))
_e("Z_C_2", "z_c", 2.0, "1", "Kollaps bei z = 2", (_D, "bei z = 2 / 1 / 0,5"))
_e("Z_C_1", "z_c", 1.0, "1", "Kollaps bei z = 1", (_D, "bei z = 2 / 1 / 0,5"))
_e("Z_C_05", "z_c", 0.5, "1", "Kollaps bei z = 0,5", (_D, "bei z = 2 / 1 / 0,5"))
_e("DELTA_200", "Δ_B", 200.0, "1", "Bindungsgrad des Bereichs, untere Grenze", (_D, "Δ = 200–500"))
_e("DELTA_500", "Δ_B", 500.0, "1", "Bindungsgrad des Bereichs, obere Grenze", (_D, "Δ = 200–500"))
_e("M_B_KLEIN", "M", 1e10, "M☉", "baryonische Masse für R_ta", (_D, "10¹⁰ / 6·10¹⁰ / 10¹² M☉"))
_e("M_B_GROSS", "M", 1e12, "M☉", "baryonische Masse für R_ta", (_D, "10¹⁰ / 6·10¹⁰ / 10¹² M☉"))
