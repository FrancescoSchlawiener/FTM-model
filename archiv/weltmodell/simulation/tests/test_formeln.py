# -*- coding: utf-8 -*-
"""Wächter: jede Funktion in kern/ und forschung/ trägt genau ihren Formel-Eintrag, Status passt zum Ort, kein Name
doppelt; Symbole von aus/ein im Register; Alternativen gültig; Stichwort des Modellorts steht in der Modelldatei."""
import os
import types

import waechter as W
from kern.formeln import FORMELN, Formel, formel, pruefe_alternativen
from kern.symbole import SYMBOLE
from pruefung.modelltext import lies


def test_bestand():
    assert W.registrierung() == []


def test_symbole_alternativen_modellort():
    assert pruefe_alternativen() == []
    for f in FORMELN.values():
        assert all(z in SYMBOLE for z in f.aus + f.ein), f.kennung
        if f.status == "modell":
            assert f.modell_ort[1] in lies(f.modell_ort[0]), (f.kennung, f.modell_ort)
        else:
            assert f.modell_ort is None


def test_negativ_unregistriert():
    mod = types.ModuleType("kern.probe")

    def ohne_eintrag(x):
        return x
    ohne_eintrag.__module__ = "kern.probe"
    mod.ohne_eintrag = ohne_eintrag
    assert any("ohne Formel-Eintrag" in f for f in W.registrierung([mod], register={}))


def test_negativ_status_am_falschen_ort():
    register = {}
    f = Formel("probe_forschung", "x = y", ("s",), ("a_N",), "LI", "forschung", "Probe")
    mod = types.ModuleType("kern.probe")

    def probe_forschung(x):
        return x
    probe_forschung.__module__ = "kern.probe"
    probe_forschung.formel, f.funktion = f, probe_forschung
    register["probe_forschung"] = f
    mod.probe_forschung = probe_forschung
    assert any("gehört nicht nach kern/" in x for x in W.registrierung([mod], register=register))


def test_negativ_dekorator_zweimal():
    import pytest
    with pytest.raises(ValueError):
        @formel("ds", "doppelt", ("Δs",), ("M",), "L", "hilfe", "Probe")
        def ds(M):
            return M
    assert os.path.basename(FORMELN["ds"].funktion.__code__.co_filename) == "faltung.py"
