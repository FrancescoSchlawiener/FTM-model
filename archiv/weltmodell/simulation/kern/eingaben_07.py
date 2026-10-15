# -*- coding: utf-8 -*-
"""Eingaben Baustein 07 (Feldgleichung): gewählte Rechenpunkte der Prüfzahlen in zahlen_07 (Klasse in eingabe.py)."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D = "07-feldgleichung.md"

# ------------------------------------------------------------------ Punkte der Rückkehr
_e("S_EINS", "s", 1.0, "1", "Punkt a_N = a₀ (Schwelle, Punktmasse bei r*)", (_D, "(s = 1)"))
_e("S_NULL_KOMMA_DREI", "s", 0.3, "1", "Punkt im Tiefbereich für den Anteil der Rückkehr als Vorgang", (_D, "(s = 0,3)"))
_e("S_TIEF", "s", 1e-9, "1", "Grenzfall s → 0 (χ = 1, volle Rückkehr) als Rechenpunkt")
_e("LAM_ZU_VARPI_GROSS", "Λ/ϖ_R", 1e4, "1", "Schluss weit außen, Λ ≫ ϖ_R (Grenzfall der bei Λ geschlossenen Rückkehr)",
   (_D, "Ist sie bei Λ ≫ ϖ_R geschlossen"))
_e("ABWEICHUNG_ZEHN", "ε_Λ", 0.10, "1", "Schwelle der Abweichung von Tully-Fisher, die der Satz nennt",
   (_D, "um mehr als 10 %"))

# ------------------------------------------------------------------ Milchstraße als Punktmasse, Sonnenkreis
_e("R_ACHT_KPC", "r", 8.0, "kpc", "Sonnenkreis der Milchstraße als Punktmasse (Beispiel im Modelltext)",
   (_D, "Milchstraße bei 8 kpc"))
_e("V_TOY_PAPER", "v", 220.0, "km/s", "Datenpaar des Toy Papers (a = v₀²/r)", (_D, "220 km/s gehören zu"))

# ------------------------------------------------------------------ Geneigte Ringe in der Miyamoto-Nagai-Scheibe
_e("MN_MASSE", "M", 6e10, "M☉", "Masse der Miyamoto-Nagai-Scheibe der Ringrechnung", (_D, "Miyamoto-Nagai-Scheibe 6·10¹⁰ M☉"))
_e("MN_A", "a_MN", 3.0, "kpc", "radiale Skala der Miyamoto-Nagai-Scheibe", (_D, "Skalen 3 und 0,3 kpc"))
_e("MN_B", "b_MN", 0.3, "kpc", "vertikale Skala der Miyamoto-Nagai-Scheibe", (_D, "Skalen 3 und 0,3 kpc"))
_e("RING_R12", "r", 12.0, "kpc", "Ringradius", (_D, "bei 12 kpc"))
_e("RING_R20", "r", 20.0, "kpc", "Ringradius", (_D, "bei 20 kpc"))
_e("RING_R30", "r", 30.0, "kpc", "Ringradius", (_D, "bei 30 kpc"))
_e("RING_I2", "i_inc", 2.0, "°", "schwach geneigter Ring", (_D, "bei 2°"))
_e("RING_I89", "i_inc", 89.0, "°", "fast polarer Ring", (_D, "bei 89°"))

# ------------------------------------------------------------------ SPARC-Punkte der Spektraldimension (Gasdruck: SIG_GAS_8)
_e("S_TRENN", "s", 0.1, "1", "Grenze, unter der die Punkte gezählt werden (die Grenze selbst ist eine Behauptung ohne "
   "Kernformel, siehe ausnahmen_07)")
