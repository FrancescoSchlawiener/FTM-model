# -*- coding: utf-8 -*-
"""Wächter: Symbole eindeutig (Zeichen, Code-Name, Grundzeichen mit Index), Ebene und Zustand aus der festen Liste."""
import pytest

import waechter as W
from kern.symbole import EBENEN, SYMBOLE, ZUSTAENDE, Symbol, SymbolRegister


def test_bestand():
    assert W.symbole_eindeutig() == []
    for s in SYMBOLE.values():
        assert set(s.ebene) <= EBENEN and s.ebene and s.zustand in ZUSTAENDE
        assert s.bedeutung and s.einheit_si and s.herkunft


def test_q_ein_objekt():
    assert "q_τ" not in SYMBOLE and "q_S" not in SYMBOLE
    assert SYMBOLE["q_Δs"].stamm == "q" and SYMBOLE["q_kosm"].stamm == "q"
    assert SYMBOLE["τ_M"].stamm == "τ" and SYMBOLE["τ_Z"].stamm == "τ"


def test_negativ_register_wirft():
    r = SymbolRegister()
    r.neu("q", "q", "Nome", "1", "LIZG", "je_system", "Autor")
    with pytest.raises(ValueError):
        r.neu("q", "q2", "anders", "1", "Z", "fest", "Baustein 05")
    with pytest.raises(ValueError):
        r.neu("q_x", "q", "anders", "1", "Z", "fest", "Baustein 05")
    with pytest.raises(ValueError):
        r.neu("p_Δs", "p_ds", "Index an fremdem Grundzeichen", "1", "Z", "fest", "Autor", stamm="q")


def test_negativ_waechter_findet():
    r = SymbolRegister()
    r.neu("q", "q", "Nome", "1", "LIZG", "je_system", "Autor")
    dict.__setitem__(r, "q2", Symbol("q2", "q", "Code doppelt", "1", "Z", "fest", "Autor"))
    assert any("Code-Name doppelt" in f for f in W.symbole_eindeutig(r))
