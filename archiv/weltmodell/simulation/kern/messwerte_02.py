# -*- coding: utf-8 -*-
"""Messwerte zu Baustein 02 (Vom Teilchen bis zum Stern); Klasse und Regeln in messwert.py."""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
_D = "02-vom-teilchen-bis-zum-stern.md"

# ------------------------------------------------------------------ Körper
_m("BETEIGEUZE_M", "M", 18.0, None, "M☉", "Beteigeuze, Masse", "G", modell_texte=((_D, "18 M☉"),))
_m("BETEIGEUZE_R", "R", 764.0, None, "R☉", "Beteigeuze, Radius", "G", modell_texte=((_D, "764 R☉"),))
_m("SATURN_R", "R", 58232.0, None, "km", "Saturn, mittlerer Radius", "G")
# Massen der Planeten (Erde über GM⊕); Faltung des Sonnensystems
_m("MERKUR_M", "M", 3.3011e23, None, "kg", "Merkur, Masse", "G")
_m("VENUS_M", "M", 4.8675e24, None, "kg", "Venus, Masse", "G")
_m("MARS_M", "M", 6.4171e23, None, "kg", "Mars, Masse", "G")
_m("JUPITER_M", "M", 1.89813e27, None, "kg", "Jupiter, Masse", "G")
_m("SATURN_M", "M", 5.6832e26, None, "kg", "Saturn, Masse", "G")
_m("URANUS_M", "M", 8.6811e25, None, "kg", "Uranus, Masse", "G")
_m("NEPTUN_M", "M", 1.02409e26, None, "kg", "Neptun, Masse", "G")

# ------------------------------------------------------------------ Sterne
_m("CARTER_GRENZE_BEOB", "M", (1.2, 1.4), None, "M☉",
   "Grenze strahlende/konvektive Hülle auf der Hauptreihe (beobachtet)", "G", modell_texte=((_D, "1,2–1,4 M☉"),))
_m("MU_SONNE_HOMOGEN", "μ_m", 0.62, None, "1", "mittlere Teilchenmasse der Sonne, homogen, solare Mischung", "G")
_m("MU_SONNE_KERN", "μ_m", 0.85, None, "1", "mittlere Teilchenmasse im Sonnenkern heute", "G")
_m("T_ZUEND_PP", "T", 4.0e6, None, "K", "Zündtemperatur des Wasserstoffbrennens (pp), Mindestwert am unteren Massenende",
   "G")
_m("T_ZUEND_3A", "T", 1.0e8, None, "K", "Zündtemperatur des Heliumbrennens (3α)", "G")
_m("T_ZUEND_SI", "T", 3.0e9, None, "K", "Zündtemperatur des Siliziumbrennens (Photozerfall)", "G")

# ------------------------------------------------------------------ Kopplungen
_m("ALPHA_MZ", "α", 1 / 127.95, None, "1", "PDG, 1/α(M_Z) = 127,95 (laufende Kopplung an der Z-Masse)", "G",
   modell_texte=((_D, "1/α(M_Z) = 127,95"),))
