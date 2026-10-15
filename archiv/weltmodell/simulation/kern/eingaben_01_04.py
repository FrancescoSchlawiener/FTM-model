# -*- coding: utf-8 -*-
"""Eingaben der Bausteine 01 und 04 (Klasse und Regeln in eingabe.py): gewählte Rechenpunkte, je Aufbau ein Abschnitt."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D04 = "04-zeit-verschraenkung-syntropie-gedaechtnis.md"

# ------------------------------------------------------------------ Gangrate an Punkten in Einheiten von r₀ (04)
_e("TAKT_R_100", "r", 100.0, "r₀", "Punkt der Gangrate weit außen, Beispiel im Modelltext", (_D04, "0,995 bei 100 r₀"))
_e("TAKT_R_2", "r", 2.0, "r₀", "Punkt der Gangrate nahe der Falte, Beispiel im Modelltext", (_D04, "0,707 bei 2 r₀"))

# ------------------------------------------------------------------ Kalter Kohlenstoff-Sauerstoff-Zwerg (04)
_e("WZ_M", "M", 0.6, "M☉", "Beispiel des Gedächtnisses im Stern", (_D04, "Kohlenstoff-Sauerstoff-Zwerg von 0,6 M☉"))
_e("WZ_R", "R", 8.696e6, "m", "Radius des kalten Zwergs von 0,6 M☉, Wert des alten Prüfstands (Masse-Radius-Beziehung "
   "der Zählregel fehlt im Kern)")

# ------------------------------------------------------------------ Einsatz der Rückgabe (04)
_e("RUECKGABE_M_KLEIN", "M", 10.0, "M☉", "untere Masse der Spanne, Beispiel im Modelltext", (_D04, "für Falten von 10 bis 10⁹ M☉"))
_e("RUECKGABE_M_GROSS", "M", 1e9, "M☉", "obere Masse der Spanne, Beispiel im Modelltext", (_D04, "für Falten von 10 bis 10⁹ M☉"))
