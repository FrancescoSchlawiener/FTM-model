# -*- coding: utf-8 -*-
"""Offene Prüfzahlen, Wiederholungen und Parametersätze: je eine Negativprobe."""
import pytest

import laden  # noqa: F401
import waechter as W
from kern.parametersatz import PARAMETERSAETZE
from pruefung.pruefzahl import Pruefzahl, PRUEFZAHLEN
from pruefung.wiederholung import Wiederholung, WiederholungListe

_D = "01-probe.md"


def _eine():
    laden.lade_alles()
    return next(p for p in PRUEFZAHLEN.values() if not p.offen)


def test_offen_ergebnis():
    a = _eine()
    o = Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, None, None, a.modell_satz, offen="fehlt_formel",
                  beschreibung="Formel X")
    r = o.pruefe()
    assert r["ergebnis"] == "OFFEN" and r["ok"] is None and r["wert"] is None


def test_negativ_offen_mit_rechnung_oder_ohne_beschreibung():
    a = _eine()
    with pytest.raises(ValueError):
        Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, a.rechnung, a.einheit, a.modell_satz,
                  offen="fehlt_formel", beschreibung="X")
    with pytest.raises(ValueError):
        Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, None, None, a.modell_satz, offen="fehlt_formel")
    with pytest.raises(ValueError):
        Pruefzahl("probe", a.baustein, a.aussage, a.modell_text, None, None, a.modell_satz, offen="vielleicht",
                  beschreibung="X")


def test_negativ_vorzeichen_fundstelle():
    """Die Stelle von „−0,13“ ist „0,13“: ein Eintrag ohne Vorzeichen deckt sie."""
    texte = {_D: "Die Kurven fallen mit −0,13 bis −0,27."}
    assert W.rueckrichtung(texte, [W.Deckung(_D, "0,13", "Kurven")])[_D] == ["0,27"]


def test_negativ_wiederholung():
    laden.lade_alles()
    liste = WiederholungListe()
    liste.neu(_D, "9,81", "Tabelle", "G_ERDE")
    with pytest.raises(ValueError):
        liste.neu(_D, "9,81", "Tabelle", "G_ERDE")
    assert Wiederholung(_D, "9,81", "Tabelle", "GIBT_ES_NICHT").eintrag() == (None, None)


def test_negativ_parametersatz_sonnenort():
    """Die Gezeit am Sonnenort rechnet nur im Satz mit gemessenem Sonnenort."""
    from pruefung.sonnenort import gezeit
    with pytest.raises(ValueError):
        gezeit(PARAMETERSAETZE["standard"])
    assert gezeit(PARAMETERSAETZE["Sonnenkreis_gemessen"]) > 0


def test_parametersaetze_vorhanden():
    for k in ("horizont", "z1090", "Upsilon_0.4", "Upsilon_0.6", "mitatmend_heute", "Sonnenkreis_Punktmasse",
              "Sonnenkreis_gemessen"):
        assert k in PARAMETERSAETZE
    h = PARAMETERSAETZE["horizont"]
    assert abs(h.a0 / 1.04e-10 - 1) < 0.01 and PARAMETERSAETZE["z1090"].z == 1090
