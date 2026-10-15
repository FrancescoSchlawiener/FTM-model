# -*- coding: utf-8 -*-
"""Wächter: Katalogeinträge und Versuche nennen Vorhandenes (streng); zählend: Archivskripte ohne Katalogeintrag,
Forschungszeilen mit „Rechnung“ ohne Versuch."""
import waechter as W
from forschung.thema import Invariante, Thema, Versuch
from katalog_eintrag import Katalog


def test_bestand():
    assert W.katalog() == [] and W.themen() == []
    assert isinstance(W.katalog_luecken(), list) and isinstance(W.themen_luecken(), list)


def test_negativ_katalog():
    k = Katalog()
    k.neu("gibt_es_nicht", "nutzlos", None, "Probe")
    k.neu("hopf_knoten", "übernommen", "keine_kennung", "Probe")
    assert len(W.katalog(k)) == 2
    assert W.katalog_luecken(k, ["hopf_knoten", "a0_lesungen"]) == ["a0_lesungen"]


def test_negativ_themen():
    t = Thema("probe", "Frage?", Invariante("Form", "Z", "keine"),
              [Versuch("v1", "2026-09-30", "Weg", "offen", forschungszeile="Gibt es nicht")], {}, "Messung", "offen")
    assert len(W.themen({"probe": t}, titel={"Hopf-Knoten"})) == 1
    zeilen = [dict(nr=1, kurztitel="Hopf-Knoten"), dict(nr=2, kurztitel="Ohne Versuch")]
    assert W.themen_luecken({"probe": t}, zeilen) == ["Z. 1: Hopf-Knoten", "Z. 2: Ohne Versuch"]


def test_urteil_auf_gedaechtnis_wird_offen():
    v = Versuch("v", "2026-09-30", "Weg", "trägt nicht", haengt_an=("DF2_SIGMA",))
    assert v.urteil.startswith("offen, hängt an DF2_SIGMA")
    assert Versuch("w", "2026-09-30", "Weg", "trägt nicht", haengt_an=("a₀",)).urteil == "trägt nicht"
