# -*- coding: utf-8 -*-
"""Eingaben der Bausteine 08 und 09 (Klasse und Regeln in eingabe.py): gewählte Rechenpunkte, je Aufbau ein Abschnitt."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
D8 = "08-licht-und-trichter.md"
D9 = "09-elektromagnetismus-und-frequenzleiter.md"

# ------------------------------------------------------------------ Trichter und Schwelle (08)
_e("ISS_HOEHE", "r", 400.0, "km", "Bahnhöhe der ISS über dem mittleren Erdradius; im Modelltext nicht genannt", wahl=True)
_e("R_MILCHSTRASSE_HUELLE", "R", 15.0, "kpc", "Ausdehnung der Milchstraße für r*/R; im Modelltext nicht genannt", wahl=True)
_e("R_HAUFEN_HUELLE", "R", 1.0, "Mpc", "Ausdehnung des Beispielhaufens für r*/R; im Modelltext nicht genannt", wahl=True)

# ------------------------------------------------------------------ Rotation und Verwindung (08)
_e("J_ISCO", "j", 0.9, "1", "Spinzahl des Beispiels der innersten stabilen Bahn", (D8, "liegt bei j = 0,9"))
_e("J_MILCHSTRASSE", "j", 1329.0, "1", "Spinzahl der Milchstraße für die Knotenzahl an r*; im Modelltext nicht genannt", wahl=True)
_e("GPB_HOEHE", "r", 642.0, "km", "Bahnhöhe von Gravity Probe B, wie im Modelltext", (D8, "642 km Höhe"))
_e("GPB_LEITSTERN_DEKL", "δ_L", 16.84, "°", "Deklination des Leitsterns IM Peg (+16°50′), Winkel zur Erdachse bei polarer Bahn", (D8, "Leitstern IM Peg"))
_e("GPB_A", "r", 7027.4, "km", "Bahnhalbachse von Gravity Probe B für die Rechnung; 642 km Höhe über dem mittleren Erdradius "
   "gäbe 7013 km", wahl=True)

# ------------------------------------------------------------------ Frequenzleiter und Kohärenz (09)
_e("T_KOERPER", "T", 310.0, "K", "Körpertemperatur für die thermische Dephasierungszeit; im Modelltext ohne Zahl", wahl=True)
_e("PROBE_MASSE", "m", 1.0, "kg", "Probekörper des Coulomb-Beispiels", (D9, "1 kg mit 1 pC bei 110 V/m"))
_e("PROBE_LADUNG", "Q_el", 1e-12, "C", "Ladung des Probekörpers im Coulomb-Beispiel", (D9, "1 kg mit 1 pC bei 110 V/m"))
