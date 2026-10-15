# -*- coding: utf-8 -*-
"""Wächter (zählend): Zeichen im Modelltext, deren Grundbuchstabe im Register eine andere Bedeutung hat."""
import waechter as W


def test_bestand_laeuft():
    assert isinstance(W.symbolkollision(), list)


def test_negativ():
    funde = W.symbolkollision({"07-probe.md": "Die Güte χ² der Anpassung und der Schalter χ; dazu q_Δs."})
    assert len(funde) == 0 or all("χ" not in f for f in funde)
    funde = W.symbolkollision({"07-probe.md": "Die Güte χ_fit der Anpassung und κ als Epizykelfrequenz."})
    assert any("χ_fit kollidiert" in f for f in funde) and any("κ nicht im Register" in f for f in funde)
