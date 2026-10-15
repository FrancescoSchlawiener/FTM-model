# -*- coding: utf-8 -*-
"""Wächter Katalog (streng): jede Archivdatei, Skript oder verwaiste Datei, hat genau einen gültigen Eintrag; jede im
Grund genannte Kennung gibt es. Mit Negativproben."""
import os
import re

import waechter as W
from katalog_eintrag import KATALOG, Katalog, archivdateien, verwaiste

_KENNUNG = re.compile(r"(?<![\w/.])([a-z0-9]+(?:_[a-z0-9]+)+)(?![\w/.])")


def _fremde_kennungen(katalog_):
    """Kleingeschriebene Namen mit Unterstrich (ab 6 Zeichen) im Grund, die keine Kennung und keine Archivdatei sind."""
    from kern.formeln import FORMELN
    from pruefung.gegenprobe import GEGENPROBEN
    from pruefung.gesetz import GESETZE
    from pruefung.pruefzahl import PRUEFZAHLEN
    bekannt = set(FORMELN) | set(PRUEFZAHLEN) | set(GESETZE) | set(GEGENPROBEN) | set(archivdateien())
    return [f"{e.skript}: {t}" for e in katalog_.values() for t in _KENNUNG.findall(e.grund)
            if len(t) >= 6 and t not in bekannt]


def test_vollstaendig_und_gueltig():
    assert "katalog_luecken" in W.STRENG and "katalog_luecken" not in W.ZAEHLEND
    assert W.katalog_luecken() == [] and W.katalog() == []
    assert set(KATALOG) == set(archivdateien())
    assert len([d for d in archivdateien() if "." not in d]) == 187


def test_grund_nennt_nur_vorhandenes():
    assert _fremde_kennungen(KATALOG) == []
    k = Katalog()
    k.neu("atmen_w", "nutzlos", None, "ersetzt durch lambda_vakuum und gibt_es_nicht")
    assert _fremde_kennungen(k) == ["atmen_w: gibt_es_nicht"]


def test_negativ_luecke():
    k = Katalog()
    for e in list(KATALOG.values())[1:]:
        k.neu(e.skript, e.klasse, e.ziel, e.grund)
    erste = next(iter(KATALOG))
    assert W.katalog_luecken(k) == [erste]
    k.neu("gibt_es_nicht.txt", "nutzlos", None, "Probe")
    assert W.katalog(k) == ["Katalog gibt_es_nicht.txt: kein Archivskript"]


def test_verwaiste(tmp_path):
    for d in ("skript.py", "skript.txt", "skript_a011.txt", "ohne.npy"):
        (tmp_path / d).write_text("", encoding="utf-8")
    os.mkdir(tmp_path / "__pycache__")
    for d in ("skript.cpython-313.pyc", "weg.cpython-313.pyc"):
        (tmp_path / "__pycache__" / d).write_bytes(b"")
    assert verwaiste(str(tmp_path)) == ["__pycache__/weg.cpython-313.pyc", "ohne.npy", "skript_a011.txt"]
    assert archivdateien(str(tmp_path)) == ["__pycache__/weg.cpython-313.pyc", "ohne.npy", "skript",
                                            "skript_a011.txt"]
