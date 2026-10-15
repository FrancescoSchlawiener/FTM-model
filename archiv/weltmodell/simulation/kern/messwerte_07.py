# -*- coding: utf-8 -*-
"""Messwerte Baustein 07 (Feldgleichung): Mond, Weltalter; der Sonnenort steht in messwerte.py (Klasse in messwert.py)."""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
_D = "07-feldgleichung.md"

# ------------------------------------------------------------------ Körper und Kosmos
_m("M_MOND", "M", 7.342e22, None, "kg", "Masse des Mondes", "G")
_m("WELTALTER", "t", 13.8, None, "Mrd. Jahre", "Alter des Kosmos (Planck 2018)", "G", modell_texte=((_D, "13,8 Mrd. Jahren"),))
