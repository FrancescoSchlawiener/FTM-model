# -*- coding: utf-8 -*-
"""Wächter (Beobachterfehler): Formel mit Ausgabe im Zustand „epoche“ hat die Epoche als Eingang; zählend: Sätze, die
eine feste Größe einer Uhrgröße gleichsetzen, ohne „heute“."""
import waechter as W
from kern.formeln import Formel


def test_bestand():
    assert W.epochen_eingabe() == []
    assert isinstance(W.uhr_muster(), list)


def test_negativ_epoche():
    f = Formel("probe", "Θ = 2πa₀/(cH)", ("Θ",), ("a₀",), "G", "hilfe", "Probe")
    assert len(W.epochen_eingabe({"probe": f})) == 1
    g = Formel("probe", "Θ = 2πa₀/(cH)", ("Θ",), ("a₀", "H"), "G", "hilfe", "Probe")
    assert W.epochen_eingabe({"probe": g}) == []


def test_negativ_uhr_muster():
    funde = W.uhr_muster({"03-probe.md": "Die Schwelle ist a₀ ≈ cH₀/2π. Heute gilt Λ ≈ πc/H₀."})
    assert len(funde) == 1 and "a₀" in funde[0]
