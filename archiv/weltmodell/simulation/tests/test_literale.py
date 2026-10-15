# -*- coding: utf-8 -*-
"""Wächter: kein Zahlenliteral einer Naturkonstante oder Einheit außerhalb von konstanten.py und messwerte.py."""
import waechter as W


def test_bestand():
    assert W.literale() == []


def test_negativ(tmp_path):
    p = tmp_path / "probe.py"
    p.write_text("M = 1.989e30\nG = 6.674e-11\nx = 2.5\n", encoding="utf-8")
    funde = W.literale([str(p)])
    assert len(funde) == 2 and "M☉" in funde[0] and "G" in funde[1]
