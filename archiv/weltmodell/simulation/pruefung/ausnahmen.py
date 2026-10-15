# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung (Klasse in ausnahme.py): Zahlen im Modelltext, die nichts zu prüfen tragen.
So füllen: AUSNAHMEN.neu(Datei, Zahl wie im Text, Grund[, stichwort=…]). Eine gerechnete oder behauptete Zahl ist nie
eine Ausnahme. Mit stichwort gilt die Ausnahme nur an dieser Fundstelle (dieselbe Zahl meint anderswo etwas anderes).
"""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

AUSNAHMEN.neu(datei("03"), "2574", "Name der Galaxie IC 2574")

# ------------------------------------------------------------------ Fundstellen: Namen und Zahlenfaktoren
AUSNAHMEN.neu(datei("03"), "10⁴", "Zahlenfaktor der Formel M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴)", stichwort="Oberhalb M_b,Q = (1 + Q)²")
AUSNAHMEN.neu(datei("03"), "14", "Name des Kugelsternhaufens Pal 14", stichwort="Pal 14")
AUSNAHMEN.neu(datei("05"), "24", "Koeffizient der Reihe Δ = q − 24q² + … (Ramanujan τ(2))", stichwort="q − 24q²")
AUSNAHMEN.neu(datei("05"), "24", "Ordnung 24 der Grundverschiebung q^{−1/24} (Dedekind-Eta), dieselbe 24 wie in E₂", stichwort="q^{−1/24}")
AUSNAHMEN.neu(datei("05"), "24", "Ordnung 24 der Grundverschiebung (Multiplikatorordnung), dieselbe 24 wie in E₂", stichwort="Multiplikatorordnung")
AUSNAHMEN.neu(datei("07"), "16", "Zahlenfaktor der Wirkung c⁴R/16πG", stichwort="c⁴R/16πG")
AUSNAHMEN.neu(datei("08"), "15", "Zahlenfaktor π⁴/15 der Formel Λ_s·ℓ_P²", stichwort="4·(π⁴/15)")
AUSNAHMEN.neu(datei("08"), "18", "Zahlenfaktor 18π der Periheldrehung", stichwort="18π·Δs")

# ------------------------------------------------------------------ Teile
from . import ausnahmen_01_04  # noqa: F401,E402  Ausnahmen der Bausteine 01 und 04
from . import ausnahmen_02  # noqa: E402,F401
from . import ausnahmen_03a  # noqa: E402,F401  (Ausnahmen 03, Teil a)
from . import ausnahmen_03b  # noqa: E402,F401  (Baustein 03, zweiter Teil)
from . import ausnahmen_05_06  # noqa: E402,F401
from . import ausnahmen_07  # noqa: E402,F401
from . import ausnahmen_08_09  # noqa: E402,F401
