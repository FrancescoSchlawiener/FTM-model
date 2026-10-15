# -*- coding: utf-8 -*-
"""Messwerte zu 03, Teil a (Klasse und Regeln in messwert.py)."""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
_BAO = "Planck Collaboration 2020, A&A 641, A6, Tab. 2 (TT,TE,EE+lowE+lensing+BAO), 1σ"

# ------------------------------------------------------------------ Kosmos: Q aus Hintergrundstrahlung mit BAO
_m("OMEGA_B_H2_BAO", "ω_b", 0.02242, 0.00014, "1", _BAO, "literatur")
_m("OMEGA_C_H2_BAO", "ω_c", 0.11933, 0.00091, "1", _BAO, "literatur")
