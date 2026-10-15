# -*- coding: utf-8 -*-
"""Eingaben (Klasse und Regeln in eingabe.py). Neue Eingaben kommen hierher, je Aufbau ein Abschnitt."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu

# ------------------------------------------------------------------ Systeme
_D09 = "09-elektromagnetismus-und-frequenzleiter.md"
_e("NS_M", "M", 1.4, "M☉", "typischer Neutronenstern: Tabelle 02, Uhr auf dem Neutronenstern 04, Frequenzleiter 09",
   (("02-vom-teilchen-bis-zum-stern.md", "Neutronenstern (1,4 M☉)"), (_D09, "Neutronenstern (1,4 M☉, 12 km)")))
_e("NS_R", "R", 12.0, "km", "typischer Neutronensternradius zu NS_M (02, 04, 09)", ((_D09, "Neutronenstern (1,4 M☉, 12 km)"),))
_e("M_HAUFEN", "M", 1.7e14, "M☉", "baryonische Masse eines Galaxienhaufens, Beispiel der Tabelle 03 und des Haufens in 08",
   (("03-von-der-galaxie-bis-zum-kosmos.md", "1,7·10¹⁴ M☉"),))
_e("Z_REKOMBINATION", "z", 1090.0, "1", "Epoche der letzten Streufläche (Rekombination), Rechenpunkt des frühen "
   "Vakuumteils und der Rückkehr bei der Rekombination (03)", (("03-von-der-galaxie-bis-zum-kosmos.md", "Ω_s(1090)"),))
_e("M_MILCHSTRASSE", "M", 6e10, "M☉", "baryonische Masse der Milchstraße als Punktmasse, Beispiel der Bausteine 03, 08, 09",
   (("03-von-der-galaxie-bis-zum-kosmos.md", "6·10¹⁰ M☉"), (_D09, "Milchstraße (6·10¹⁰ M☉)")))

# ------------------------------------------------------------------ Labor: eingespannter Kupferdraht (06), Annahmen des Aufbaus
_e("DRAHT_A", "a_D", 50e-6, "m", "Drahtradius, Annahme des Aufbaus (torus_wechselkraefte, Teil 4a)")
_e("DRAHT_ELL", "ℓ_D", 0.02, "m", "Drahtlänge, Annahme des Aufbaus")
_e("DRAHT_B", "B", 0.02, "T", "Magnetfeld, Annahme des Aufbaus")
_e("DRAHT_QM", "Q_m", 1000.0, "1", "mechanische Güte, Annahme des Aufbaus")
_e("DRAHT_SIG0", "σ_0", 20e6, "Pa", "Vorspannung, Annahme des Aufbaus")
from . import eingaben_01_04  # noqa: F401,E402  Eingaben der Bausteine 01 und 04
from . import eingaben_02  # noqa: E402,F401
from . import eingaben_03a  # noqa: E402,F401  (Eingaben 03, Teil a)
from . import eingaben_03b  # noqa: E402,F401  (Baustein 03, zweiter Teil)

from . import eingaben_05_06  # noqa: E402,F401
from . import eingaben_07  # noqa: E402,F401
from . import eingaben_08_09  # noqa: E402,F401
from . import eingaben_teilchen_licht  # noqa: E402,F401  (Gruppe Teilchen, Sterne und Licht)
from . import eingaben_zustand  # noqa: E402,F401
from . import eingaben_galaxie  # noqa: E402,F401  (Gruppe Galaxie, Auftrag F)
from . import eingaben_kosmos  # noqa: E402,F401  (Gruppe Kosmos)
