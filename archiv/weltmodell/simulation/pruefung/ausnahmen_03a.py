# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung zu 03, Teil a (Klasse in ausnahme.py): Beiwerte einer Formel, die der Kern trägt.
Offene Zahlen stehen als offene Prüfzahlen in zahlen_offen.py, Messwerte mit ihrem Ort in messwerte.py.
"""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_x = AUSNAHMEN.neu
D = datei("03")

# ------------------------------------------------------------------ Beiwerte
_x(D, "81", "Beiwert der Formel σ_v⁴ = (4/81)·GMa₀ (Kern sigma_v_tief)")
