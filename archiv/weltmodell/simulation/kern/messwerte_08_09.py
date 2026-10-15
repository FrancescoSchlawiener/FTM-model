# -*- coding: utf-8 -*-
"""
Messwerte der Bausteine 08 und 09 (Klasse und Regeln in messwert.py): Messungen, gegen die der Modelltext die Rechnung
stellt, und Eingänge der Wellenquellen. Messwerte aus messwerte.py tragen ihre Orte in 08 dort am Eintrag.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
D8 = "08-licht-und-trichter.md"

# ------------------------------------------------------------------ Messungen im Modelltext (08)
_m("GPB_GEO", "Ω_geo", 6602.0, 18.0, "mas/Jahr", "Gravity Probe B (Everitt et al. 2011), geodätische Drehung, 1σ", "literatur",
   modell_texte=((D8, "6602 ± 18"),))
_m("GPB_LT", "Ω_LT", 37.2, 7.2, "mas/Jahr", "Gravity Probe B (Everitt et al. 2011), Mitführung, 1σ", "literatur",
   modell_texte=((D8, "37,2 ± 7,2"),))
_m("HT_VERHAELTNIS", "Ṗ_b", 0.998, 0.002, "1",
   "PSR B1913+16: gemessene zu erwarteter Abnahme der Umlaufzeit, galaktisch korrigiert (Weisberg, Huang 2016)", "G",
   modell_texte=((D8, "0,998 ± 0,002"),))
_m("EHT_RING", "α_S", 42.0, 3.0, "μas", "EHT Collaboration 2019, ApJL 875, L1 (arXiv:1906.11238): Ringdurchmesser 42 ± 3 μas", "literatur", modell_texte=((D8, "42 ± 3 μas"),))
_m("GW170817_DT", "t", 1.74, 0.05, "s", "GW170817 / GRB 170817A: Gammablitz nach dem Verschmelzen (LIGO, Virgo, Fermi 2017)",
   "G", modell_texte=((D8, "um 1,74 s versetzt"),))

# ------------------------------------------------------------------ Eingänge der Wellenquellen (08)
_m("GW150914_Z", "z", 0.09, 0.03, "1", "GW150914, Rotverschiebung der Quelle (LIGO 2016)", "literatur")
_m("GW150914_MC_DET", "M", 30.2, 2.5, "M☉", "GW150914, Chirp-Masse im Detektor (LIGO 2016, PRL 116, 241102)", "G")
_m("GW170817_MC_DET", "M", 1.1977, 0.0008, "M☉", "GW170817, Chirp-Masse im Detektor (LIGO, Virgo 2017)", "literatur")
