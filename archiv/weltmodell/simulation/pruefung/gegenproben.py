# -*- coding: utf-8 -*-
"""
Gegenproben: je Formel des Modells mindestens ein unabhängiger Weg (Felder in gegenprobe.py). So füllen:
GEGENPROBEN.neu(kennung, Formel-Kennung, "zweiter_weg" | "grenzfall" | "handwert", [Argument-Tupel], unabhaengig).
Der unabhängige Weg ruft die geprüfte Formel nicht auf (Wächter). Kennungen: "<Formel>_<Weg>". Die Einträge stehen
je Gruppe in gegenproben_*.py; hier stehen zweite Wege, die mehrere Gruppen nutzen.
"""
import math

import kern as K

from .gegenprobe import GEGENPROBEN  # noqa: F401


def hawking_lebensdauer(M):
    """Lebensdauer der Falte über die Hawking-Form t = 5120πG²M³/(ħc⁴) (Page 1976): zweiter Weg zu lebensdauer_falte
    (Gegenprobe lebensdauer_falte_page)."""
    return 5120 * math.pi * K.G**2 * M**3 / (K.hbar * K.c**4)
