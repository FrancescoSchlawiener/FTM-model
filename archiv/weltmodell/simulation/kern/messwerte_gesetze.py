# -*- coding: utf-8 -*-
"""
Messwerte der Gesetze (Klasse und Regeln in messwert.py; alle aus dem Gedächtnis [G]): Steigung der baryonischen
Tully-Fisher-Relation und Fehler der Doppelpotenz von Lelli. Die Werte von Lelli 2016 (Σ_M, Σ_C, α, β) stehen in
messwerte.py (LELLI_SIGMA_DYN); hier nur ihre Fehler.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu

_m("LELLI_BTF_STEIGUNG", "v_flach", 3.85, 0.09, "1",
   "Lelli, McGaugh, Schombert, Desmond, Katz 2019 (MNRAS 484, 3267): Steigung der baryonischen Tully-Fisher-Relation "
   "log M_b gegen log V_f, 1σ", "G")
_m("LELLI_BETA_FEHLER", "Σ_dyn", 0.04, None, "1",
   "Lelli, McGaugh, Schombert, Pawlowski 2016 (ApJL 827, L19): Fehler von β = 0,62 der Doppelpotenz", "G")
_m("LELLI_STREUUNG", "Σ_dyn", 0.21, None, "dex",
   "Lelli, McGaugh, Schombert, Pawlowski 2016 (ApJL 827, L19): Streuung um die Doppelpotenz", "G")
