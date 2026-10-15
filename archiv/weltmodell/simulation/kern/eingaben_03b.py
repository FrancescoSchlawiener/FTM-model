# -*- coding: utf-8 -*-
"""
Eingaben zu Baustein 03, zweiter Teil (Expansion, Kaleidoskop, Dunkle Energie, Vorhersagen, Offen). Klasse und Regeln in
eingabe.py. Gewählte Rechenpunkte und Annahmen des Aufbaus, kein [G].
"""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D = "03-von-der-galaxie-bis-zum-kosmos.md"

# ------------------------------------------------------------------ Epochen und Bereiche
_e("Z_SCHWELLE_MESSUNG", "z", (0.5, 1.0), "1", "Bereich, in dem die Beschleunigungsrelation die Wanderung von a₀ misst",
   (_D, "bei z = 0,5–1"))
_e("KAE_DOPPELSTERNE", "r", (5.0, 30.0), "kAE", "Abstandsbereich weiter Doppelsterne der Vorhersage (Gaia)",
   (_D, "bei 5–30 kAE"))

# ------------------------------------------------------------------ Punktmasse Milchstraße: Sonnenort, Rückkopplung
_e("R_SONNE_PUNKTMASSE", "r", 8.2, "kpc", "Abstand der Sonne vom Zentrum als Rechenpunkt der Punktmasse (Binärpulsare)")
_e("R_RUECKKOPPLUNG", "r", 50.0, "kpc", "Rechenpunkt der Rückkopplungszeit ν·r/c", (_D, "ν·r/c ≈ 1 Mio. Jahre bei 50 kpc"))

# ------------------------------------------------------------------ Übergabe nach Haufenstößen
_e("I_F_HAUFENSTOSS", "N*", (1.0, 1.25), "rad",
   "Übergabe bei der ersten Faltung (Rechnung uebergabe_faltung); keine Kernformel, als Wert eingesetzt")

# ------------------------------------------------------------------ Vorhersagen: gewählte Punkte
_e("R_FRUEH_50", "r", 50.0, "kpc", "äußerer Rechenpunkt früher Galaxien", (_D, "(50 kpc)"))
_e("Z_FRUEH_34", "z_c", (3.4, 2.0), "1", "Kollaps- und Beobachtungsepoche einer früh kollabierten Galaxie",
   (_D, "Kollaps bei z = 3,4, gesehen bei z = 2"))
_e("Z_FRUEH_05", "z_c", (1.0, 0.5), "1", "Kollaps- und Beobachtungsepoche einer spät kollabierten Galaxie",
   (_D, "gesehen bei z = 0,5"))
_e("R_LINSE_100", "r", 100.0, "kpc", "Radius, bis zu dem die Mode um isolierte Galaxien ganz übergegangen ist",
   (_D, "bis 100 kpc ganz übergegangen"))
_e("R_LINSE_200", "r", 200.0, "kpc", "Radius der Linse um isolierte Galaxien", (_D, "(200 kpc)"))
_e("R_LINSE_300", "r", 300.0, "kpc", "Radius der Linse um isolierte Galaxien", (_D, "(300 kpc)"))

# ------------------------------------------------------------------ SPARC-Außenpunkte, äußeres Feld
_e("SPARC_RAND_AUSWAHL", "a_N", 0.05, "a₀", "Auswahl: äußerster Punkt einer Kurve mit a_N ≤ 0,05·a₀")
_e("GE_RAND_001", "a_e", 0.01, "a₀", "geprüftes äußeres Feld an den Außenpunkten", (_D, "0,01·a₀ kaum"))
_e("GE_RAND_0005", "a_e", 0.005, "a₀", "Grenze der Empfindlichkeit der Außenpunkte", (_D, "unter 0,005·a₀"))
_e("GE_RAND_0001", "a_e", 0.001, "a₀", "geprüftes kleines äußeres Feld (Linse)", (_D, "0,001·a₀ machte"))

# ------------------------------------------------------------------ Warp: Aufbau der Rechnung
_e("WARP_RADIUS", "r", 14.0, "kpc", "Radius der Warp-Rechnung (Drehmomente)", (_D, "bei 14 kpc)"))
_e("WARP_NEIGUNG", "ϑ", 0.1, "rad", "Neigung des Rings in der Warp-Rechnung (Drehmomente)", (_D, "Neigung 0,1,"))
_e("WARP_RINGE", "r", (12.0, 16.0), "kpc", "Ringe der Aufwicklung", (_D, "bei 12 und 16 kpc"))
_e("WARP_AUSEINANDER", "Δθ", 180.0, "°", "Maß der Aufwicklung: Knotenlinien um 180° auseinander", (_D, "um 180° auseinander"))

# ------------------------------------------------------------------ Druckgestützte Systeme: Wirt und Punkte
_e("M_NGC1052", "M", 1e11, "M☉", "baryonische Masse von NGC 1052 als Punktmasse (Wirt von DF2, DF4), Annahme des Aufbaus")

# ------------------------------------------------------------------ Exponentialscheibe gegen Lelli-Verlauf
_e("SIGMA0_EXPSCHEIBE", "Σ₀", (10.0, 1000.0), "M☉/pc²", "Rechenpunkte niedriger und hoher Flächendichte")
