# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung für die Bausteine 01 und 04 (Klasse in ausnahme.py).
"Koeffizient:" fester Faktor einer Formel, die im Kern steht; geprüft wird die Formel (Gegenprobe, Gesetz), nicht
die Zahl. Offene Zahlen (fehlt Formel, nur Forschung) stehen als offene Prüfzahlen in zahlen_offen.py.
"""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_D04 = datei("04")
_a = AUSNAHMEN.neu

# ------------------------------------------------------------------ 04 Koeffizienten
_a(_D04, "24", "Koeffizient: T_H = ħc/(24π·k·Δs) in grenzmasse_hawking")
_a(_D04, "3840", "Koeffizient: t = 3840·(S/k)·Δs/c in lebensdauer_falte")
_a(_D04, "15", "Koeffizient: Vakuumteil π⁴/(15Λ²) in lambda_vakuum")
_a(_D04, "45", "Koeffizient: H_s = π²c/(√45·Λ) in H_vakuum")
