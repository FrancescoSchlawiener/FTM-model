# -*- coding: utf-8 -*-
"""Forschung: Das Tunnelbild der Rückkehr (nicht im Modell) rechnet weiter; ψ(r*) = q auf 10⁻⁴. Kein Gesetz des Modells."""
from forschung.tunnelbild import tunnelbild_gegenprobe
from pruefung.gesetz import GESETZE
from pruefung.modelltext import MODELL
import glob


def test_tunnelbild_rechnung():
    psi, q = tunnelbild_gegenprobe()
    assert abs(psi / q - 1) < 1e-4


def test_tunnelbild_nicht_im_modell():
    assert "M_rueckkehr_tunneln" not in GESETZE
    text = " ".join(open(f, encoding="utf-8").read() for f in glob.glob(MODELL + "/*.md"))
    assert "Tunnel" not in text and "tunnel" not in text
