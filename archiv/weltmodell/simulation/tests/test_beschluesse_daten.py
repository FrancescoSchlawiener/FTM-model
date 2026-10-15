# -*- coding: utf-8 -*-
"""Wächter: Wiederöffnen eines Beschlusses nur mit neuer Grundlage; Hashes der Datendateien unverändert."""
import waechter as W
from forschung.beschluss import BeschlussRegister


def test_bestand():
    assert W.beschluesse() == [] and W.daten() == []


def test_negativ_beschluss():
    r = BeschlussRegister()
    r.neu("alt", "Aussage", ("03-probe.md", "x"), "Rechnung A", "2026-09-29", "hinein")
    r.neu("neu", "Aussage", ("03-probe.md", "x"), "Rechnung A", "2026-09-30", "heraus", ersetzt="alt")
    assert len(W.beschluesse(r)) == 1


def test_negativ_daten():
    soll = {"sparc/a.dat": "11", "elemente/b.csv": "22"}
    assert W.daten({"sparc/a.dat": "11", "elemente/b.csv": "23"}, soll) == ["Daten geändert: elemente/b.csv"]


def test_doku_abschnitt_9_zeigt_offene_beschluesse():
    """Abschnitt 9 der Doku führt jeden offenen Beschluss und jeden Vorschlag mit Kennung und Satz."""
    import doku_erzeugen
    from forschung.beschluesse import OFFEN, VORSCHLAEGE
    text = "\n".join(doku_erzeugen.offen({"themen": {}, "pruefzahlen": {}, "gesetze": {}, "beschluesse": {}}, None))
    for k in list(OFFEN) + list(VORSCHLAEGE):
        assert f"| {k} |" in text
    assert all(v.satz.replace("|", "\\|") in text for v in VORSCHLAEGE.values())
