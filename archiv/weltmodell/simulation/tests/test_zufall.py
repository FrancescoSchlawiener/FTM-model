# -*- coding: utf-8 -*-
"""Wächter: Zufall nur mit eigenem Generator mit Seed, je Aufruf (nicht auf Modulebene geteilt)."""
import waechter as W


def test_bestand():
    assert W.zufall() == []


def test_negativ(tmp_path):
    p = tmp_path / "probe.py"
    p.write_text("import random\nimport numpy as np\nRNG = np.random.default_rng(1)\n"
                 "def f():\n    a = np.random.normal()\n    b = np.random.default_rng()\n"
                 "    c = np.random.default_rng(7)\n    return a, b, c\n", encoding="utf-8")
    funde = W.zufall([str(p)])
    assert len(funde) == 4, funde
