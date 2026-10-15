# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung Baustein 07 (Klasse in ausnahme.py): Name, Formelteil, Vertrauensniveau – nichts zu
prüfen. Offene Zahlen (fehlt Formel, nur Forschung) stehen als offene Prüfzahlen in zahlen_offen.py.
"""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_D = datei("07")
_a = AUSNAHMEN.neu

# ------------------------------------------------------------------ ausnahme
_a(_D, "24", "ausnahme: Koeffizient der Reihe Δ = q − 24q² + 252q³ − … (Ramanujan τ(2), Definition der Spitzenform)")
_a(_D, "252", "ausnahme: Koeffizient der Reihe Δ (Ramanujan τ(3), Definition der Spitzenform)")
_a(_D, "95", "ausnahme: Vertrauensniveau 95 % der SPARC-Anpassungen")
_a(_D, "2574", "ausnahme: Name der Galaxie IC 2574")
_a(_D, "720", "ausnahme: Bernoulli-Koeffizient der Reihe r/r* + ½ + r*/(12r) − (r*/r)³/720 (Formelteil)")
