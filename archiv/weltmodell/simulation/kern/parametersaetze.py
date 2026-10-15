# -*- coding: utf-8 -*-
"""Benannte Parametersätze (Klasse in parametersatz.py). "standard" gilt, wenn eine Rechnung keinen anderen nennt.
Die Werte der Sätze kommen aus Messwerten und Eingaben (je eine Quelle); jede Prüfzahl nennt ihren Satz."""
from .eingaben import EINGABEN
from .messwerte import MESSWERTE
from .parametersatz import PARAMETERSAETZE

REGISTERMODUL = True
_p = PARAMETERSAETZE.neu

_p("standard", "standard", MESSWERTE["H0_KMS"].wert, 0.5, 0.0, "gebunden", "hier, heute", H0_quelle="H0_KMS")
_H0 = MESSWERTE["H0_KMS"].wert
# a₀ als Übereinstimmung von heute cH₀/2π (Lesart des Hubble-Horizonts; gilt nur heute, für den eigenen Horizont)
_p("horizont", "horizont", _H0, 0.5, 0.0, "gebunden", "hier, heute (eigener Horizont)", H0_quelle="H0_KMS")
# Epoche der letzten Streufläche
_p("z1090", "standard", _H0, 0.5, EINGABEN["Z_REKOMBINATION"].wert, "gebunden", "hier, bei der Rekombination",
   H0_quelle="H0_KMS")
# Masse-zu-Licht der Scheibe: untere und obere Lesung an SPARC
_p("Upsilon_0.4", "standard", _H0, EINGABEN["UPS_D_04"].wert, 0.0, "gebunden", "hier, heute", H0_quelle="H0_KMS")
_p("Upsilon_0.6", "standard", _H0, EINGABEN["UPS_D_06"].wert, 0.0, "gebunden", "hier, heute", H0_quelle="H0_KMS")
# Υ_b mit Υ_d mitskaliert (Lesart des Archivs): der Modelltext setzt Υ_b = 0,7 fest, dies ist eine Wahl
_p("Upsilon_skaliert", "standard", _H0, 0.5, 0.0, "gebunden", "hier, heute", H0_quelle="H0_KMS", ups_b_skaliert=True,
   wahl="Υ_b mit Υ_d mitskaliert (Lesart des Archivs, im Modelltext nicht gesetzt)")
# Hintergrund zwischen den Strukturen: atmet mit (θ̃ = 3H heute)
_p("mitatmend_heute", "standard", _H0, 0.5, 0.0, "mitatmend", "zwischen den Strukturen, heute", H0_quelle="H0_KMS")
# Sonnenkreis: Milchstraße als Punktmasse (Rechenpunkt 8 kpc) oder gemessener Sonnenort (v₀, R₀, Oort A, B, ρ_b)
_p("Sonnenkreis_Punktmasse", "standard", _H0, 0.5, 0.0, "gebunden", "Sonnenkreis, Milchstraße als Punktmasse",
   H0_quelle="H0_KMS", sonnenkreis="punktmasse")
_p("Sonnenkreis_gemessen", "standard", _H0, 0.5, 0.0, "gebunden", "Sonnenort, gemessen", H0_quelle="H0_KMS",
   sonnenkreis="gemessen")
