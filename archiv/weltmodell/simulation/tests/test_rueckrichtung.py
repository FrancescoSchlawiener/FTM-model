# -*- coding: utf-8 -*-
"""Wächter (zählend): jede Zahl in modell/01–09 ist an ihrer Fundstelle Prüfzahl (auch offen), Messwert, Eingabe,
Setzung oder Ausnahme. Offene Prüfzahlen erfassen, decken aber nicht. Hier Lauf und Negativproben."""
import pytest

import waechter as W
from pruefung.ausnahme import Ausnahme
from pruefung.modelltext import zahlen

_D = "01-probe.md"


def test_bestand_laeuft():
    assert isinstance(W.rueckrichtung(), dict)
    b = W.rueckrichtung_bilanz()
    assert all(v["erfasst"] == v["gedeckt"] + v["offen"] <= v["zahlen"] for v in b.values())


def test_negativ():
    texte = {_D: "Die Faltung ist 492,2 m, die Neigung 274 m/s², die Masse 3,5·10²² kg."}
    luecken = W.rueckrichtung(texte, [W.Deckung(_D, "492,2", "Faltung")])[_D]
    assert luecken == ["274", "3,5·10²²"]


def test_negativ_je_fundstelle():
    """Dieselbe Zahl in zwei Sätzen: gedeckt nur im Satz mit dem Stichwort, die andere Stelle bleibt offen."""
    texte = {_D: "Die Schale endet bei 0,53. Die Scheibe gibt bei R_d/4 den Wert 0,53."}
    assert W.rueckrichtung(texte, [W.Deckung(_D, "0,53", "Schale")])[_D] == ["0,53"]
    assert W.rueckrichtung(texte, [W.Deckung(_D, "0,53", "Schale"), W.Deckung(_D, "0,53", "Scheibe")])[_D] == []
    # Ausnahme ohne Stichwort gilt an jeder Stelle (Namen, Faktoren)
    assert W.rueckrichtung(texte, [W.Deckung(_D, "0,53", ueberall=True)])[_D] == []


def test_negativ_offen_erfasst_nicht_gedeckt():
    texte = {_D: "Die Erde exportiert 6,1·10¹⁴ W/K. Die Sonne 492,2 m."}
    e = [W.Deckung(_D, "6,1·10¹⁴", "exportiert", art="offen"), W.Deckung(_D, "492,2", "Sonne")]
    assert W.rueckrichtung(texte, e)[_D] == []
    assert W.rueckrichtung_bilanz(texte, e)[_D] == dict(zahlen=2, erfasst=2, gedeckt=1, offen=1)


def test_negativ_blosse_potenz():
    """Bloße Zehnerpotenzen mitten im Satz sind Zahlen (vorher las die erste Form „10“ und verwarf es)."""
    assert zahlen("heute nur zu 10⁻¹⁸, rund 10⁻²⁹ je Jahr, 10¹² Jahre, z ≈ 20 (10⁸ M☉)") == \
        ["10⁻¹⁸", "10⁻²⁹", "10¹²", "20", "10⁸"]
    assert zahlen("von 10 bis 10⁹ M☉, 3,5·10²² kg") == ["10⁹", "3,5·10²²"]


def test_negativ_ausnahme_offen():
    """Eine offene Zahl ist keine Ausnahme."""
    with pytest.raises(ValueError):
        Ausnahme(_D, "0,54", "fehlt Formel: Page-Zeit")
    with pytest.raises(ValueError):
        Ausnahme(_D, "0,54", "offen, nur Forschung: …")
