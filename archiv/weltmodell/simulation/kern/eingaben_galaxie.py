# -*- coding: utf-8 -*-
"""Eingaben der Gruppe Galaxie (Auftrag F): gewählte Rechenpunkte der Prüfzahlen in pruefung/zahlen_galaxie.py
(Klasse in eingabe.py)."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D03 = "03-von-der-galaxie-bis-zum-kosmos.md"
_D07 = "07-feldgleichung.md"

# ------------------------------------------------------------------ Feldform in der Exponentialscheibe (07)
_e("RD_RSTERN_HOCH", "R_d", 0.3, "r*", "Exponentialscheibe hoher Flächendichte (Skala in r* der Gesamtmasse)",
   (_D07, "R_d = 0,3 r*"))
_e("RD_RSTERN_MITTEL", "R_d", 1.0, "r*", "Exponentialscheibe mittlerer Flächendichte, Rechenpunkt des Gleichstands")
_e("RD_RSTERN_TIEF", "R_d", 3.0, "r*", "Exponentialscheibe niedriger Flächendichte", (_D07, "R_d = 3 r*"))
_e("X_RD_VIERTEL", "r", 0.25, "R_d", "Punkt R_d/4 in der Scheibenebene")
_e("X_RD_EINS", "r", 1.0, "R_d", "Punkt R_d in der Scheibenebene")
_e("X_NAEHERUNG", "r", (0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0), "R_d",
   "Punkte in der Ebene, an denen die örtliche Näherung mit der Feldform verglichen wird (R_d/4 bis 8R_d)")

# ------------------------------------------------------------------ Scheibenmodell der Milchstraße (Sonnenort, Ringe)
_e("RD_MILCHSTRASSE_26", "R_d", 2.6, "kpc", "Skalenlänge der Milchstraße als dünne Exponentialscheibe (6·10¹⁰ M☉) "
   "für Sonnenort und Ringpräzession")
_e("H_SCHICHT_03", "z_h", 0.3, "kpc", "Höhe der Schicht, Mittelebenendichte Σ/(2h), für die Vertikalfrequenz")
_e("R_RINGE_AUSSEN", "r", (20.0, 25.0), "kpc", "Ringe der Außenscheibe, an denen Kern- und Feldform verglichen werden")
_e("R_SONNE_SCHEIBE", "r", 8.2, "kpc", "Sonnenort im Scheibenmodell (Säule und Feldform)")
_e("Z_SAEULE", "z_h", 1.1, "kpc", "Höhe der Säule der Dynamikmasse am Sonnenort", (_D07, "(±1,1 kpc)"))
_e("MN_SUMME_SONNE", "M", ((1.2e10, 0.0, 0.6), (5.3e10, 3.0, 0.45), (1.2e10, 3.3, 1.0), (1.3e10, 7.0, 0.15)),
   "M☉, kpc, kpc", "Scheibenmodell aus vier Miyamoto-Nagai-Scheiben (M, a, b), angepasst an g_N, ρ_b und 𝒯 der "
   "Sonnenumgebung (Annahme des Aufbaus, wep_glaettung)")

# ------------------------------------------------------------------ Feldform der Milchstraße (McMillan, Gitter)
_e("GITTER_MILCHSTRASSE", "r", (300.0, 2000.0, 0.005), "1, kpc, kpc",
   "Gitter des Feldform-Lösers für das Massenmodell: Zellen je Richtung, Rand, erste Zelle")

# ------------------------------------------------------------------ SPARC-Lesungen
_e("SIGMA_INT", "σ_ε", (0.10, 0.20), "1", "innere Streuung von ln g_obs in der profilierten Lesung (zwei Annahmen)")
_e("R_3RSTERN", "r", 3.0, "r*", "Punktmasse bei 3r*: gegenläufige Bahnen, wenn die Rückkehr am Drehsinn hinge")
_e("S_MAX_BODEN", "s", 0.3, "1", "tiefe Punkte der Bodenanpassung (s < 0,3)")
