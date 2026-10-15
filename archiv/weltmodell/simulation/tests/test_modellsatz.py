# -*- coding: utf-8 -*-
"""Wächter: jede Prüfzahl steht mit ihrem Stichwort im selben Satz der Modelldatei; Orte von Aussagen, Setzungen,
Eingaben und Messwerten stehen wörtlich im Modelltext; Toleranzregel aus der Stellenzahl."""
import pytest

import waechter as W
from pruefung.aussage import AussageRegister
from pruefung.modelltext import zahl_aus_text
from pruefung.pruefzahl import Pruefzahl, PRUEFZAHLEN


def test_bestand():
    assert W.modellsatz() == []


def test_toleranzregel():
    assert zahl_aus_text("492,2") == (492.2, 0.05)
    w, h = zahl_aus_text("5,906·10⁻³⁹")
    assert abs(w / 5.906e-39 - 1) < 1e-12 and abs(h / 5e-43 - 1) < 1e-12
    assert zahl_aus_text("−0,13")[0] == -0.13 and zahl_aus_text("10¹²")[1] is None


def test_toleranzregel_nullen_am_ende():
    """Nullen am Ende einer ganzen Zahl gelten nicht als Stelle; hinter dem Komma schon (Negativprobe: 2780,0)."""
    assert zahl_aus_text("2780")[1] == 5 and zahl_aus_text("7 340")[1] == 5 and zahl_aus_text("1000")[1] == 500
    assert zahl_aus_text("7342")[1] == 0.5 and zahl_aus_text("2780,0")[1] == 0.05
    assert zahl_aus_text("20·10³")[1] == 5e3


def test_negativ_zahl_in_anderem_satz():
    a = next(iter(PRUEFZAHLEN.values()))
    falsch = Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, a.rechnung, a.einheit, "Grundannahme")
    leer = AussageRegister()
    assert len(W.modellsatz({"probe": falsch}, leer)) == 1


def test_negativ_eigene_toleranz_ohne_grund():
    a = next(iter(PRUEFZAHLEN.values()))
    with pytest.raises(ValueError):
        Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, a.rechnung, a.einheit, a.modell_satz, toleranz=0.05)


def test_rund_etwa_bloße_mantisse_halbe_stelle():
    """„rund 4·10¹¹“: halbe letzte Stelle der Mantisse, keine Sondertoleranz (04_rueckgabe_klein)."""
    from pruefung.modelltext import zahl_aus_text
    from pruefung.pruefzahl import PRUEFZAHLEN
    assert zahl_aus_text("4·10¹¹") == (4e11, 0.5e11)
    p = PRUEFZAHLEN["04_rueckgabe_klein"]
    assert p.toleranz is None and abs(p.toleranz_absolut() - 0.5e11) < 1
