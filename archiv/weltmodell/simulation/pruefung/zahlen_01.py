# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 01 (Grundlagen). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "01_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
"""
import kern as K
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "01"
D = datei(B)

# Beispiel (ein echter Eintrag; weitere darunter im selben Aufbau)


AUSSAGEN.neu("01_A2_je_erdmasse", (D, "= x · 1,4783 mm"), "Zahl")
PRUEFZAHLEN.neu("01_ds_je_erdmasse", B, "01_A2_je_erdmasse", "1,4783", lambda p: K.ds(K.M_earth) * 1e3, "mm", "(A2)",
                haengt_an=("GM_earth",))

AUSSAGEN.neu("01_grenzradius", (D, "Am Grenzradius r₀ steht der Takt still (Erde 8,87 mm, Sonne 2,95 km)"), "Zahl")
PRUEFZAHLEN.neu("01_r0_erde", B, "01_grenzradius", "8,87", lambda p: K.r0(K.M_earth) * 1e3, "mm", "Grenzradius r₀",
                haengt_an=("GM_earth",))
PRUEFZAHLEN.neu("01_r0_sonne", B, "01_grenzradius", "2,95", lambda p: K.r0(K.M_sun) / 1e3, "km", "Grenzradius r₀",
                haengt_an=("GM_sun",))
