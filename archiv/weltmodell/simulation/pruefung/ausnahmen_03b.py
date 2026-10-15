# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung zu Baustein 03, zweiter Teil (Klasse in ausnahme.py): Koeffizienten in Formeln, Namen,
Vertrauensniveaus. Keine gerechnete oder behauptete Zahl.
"""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_D = datei("03")
AUSNAHMEN.neu(_D, "15", "Koeffizient der Formel Λ_s = π⁴/(15·Λ²)")
AUSNAHMEN.neu(_D, "45", "Koeffizient der Formeln Ω_s = (π²/45)·Θ² und S_Λ = (45/π³)·k·(Λ/ℓ_P)²")
AUSNAHMEN.neu(_D, "2419", "Name des Kugelsternhaufens NGC 2419")
AUSNAHMEN.neu(_D, "95", "Vertrauensniveau 95 % der SPARC-Anpassung")
