# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 03 (Von der Galaxie bis zum Kosmos). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "03_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
"""
import kern as K
from kern.eingabe import EINGABEN
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "03"
D = datei(B)

# Beispiel (ein echter Eintrag; weitere darunter im selben Aufbau)
AUSSAGEN.neu("03_milchstrasse_punktmasse", (D, "Milchstraße als Punktmasse: v = 250 km/s bei 5 kpc"), "Zahl")
PRUEFZAHLEN.neu("03_v_flach_milchstrasse", B, "03_milchstrasse_punktmasse", "172,0",
                lambda p: K.v_flach(EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun, p.a0) / 1e3, "km/s",
                "Milchstraße als Punktmasse", haengt_an=("M_MILCHSTRASSE", "GM_sun"))
