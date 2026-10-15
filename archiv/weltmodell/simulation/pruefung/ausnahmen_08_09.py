# -*- coding: utf-8 -*-
"""Ausnahmen der Rückrichtung in 08 und 09 (Klasse in ausnahme.py): Zahlenfaktoren in Formeln."""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_a = AUSNAHMEN.neu
D8, D9 = datei("08"), datei("09")

# ------------------------------------------------------------------ 08: Zahlenfaktoren in Formeln
_a(D8, "27", "Zahlenfaktor der Reihe 27/2·(Δs/r)³ und von b = 3√27·Δs")
_a(D8, "48", "Zahlenfaktor im Quellterm □h̄_μν = −48π·σ_μν")
_a(D8, "32", "Zahlenfaktor 32/5 der Quadrupol-Leistung")
_a(D8, "192", "Zahlenfaktor 192π/5 der Abnahme der Umlaufzeit")
_a(D8, "96", "Zahlenfaktor 96/5 des Faltungsmessers")
