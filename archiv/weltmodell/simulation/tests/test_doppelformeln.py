# -*- coding: utf-8 -*-
"""Wächter: keine inhaltlichen Doppelformeln (gleiche Einheit, gleiche Werte an Zufallseingaben, ohne Alternative)."""
import math

import pytest
import waechter as W
from kern.formeln import FORMELN, Formel


@pytest.mark.schwer
def test_bestand():
    assert W.doppelformeln() == []


def _formel(kennung, fn, alternative=None, ebene="L", gleiche_form=()):
    f = Formel(kennung, "Probe", ("r*",), ("M", "a₀"), ebene, "hilfe", "Probe", alternative=alternative,
               gleiche_form=gleiche_form)
    f.funktion = fn
    return f


def test_negativ():
    a = _formel("probe_a", lambda M: math.sqrt(M) * 3)
    b = _formel("probe_b", lambda M: 3 * M ** 0.5)
    assert len(W.doppelformeln({"probe_a": a, "probe_b": b})) == 1
    b_alt = _formel("probe_b", lambda M: 3 * M ** 0.5, alternative="probe_a")
    assert W.doppelformeln({"probe_a": a, "probe_b": b_alt}) == []
    assert "ds" in FORMELN


def test_gleiche_form_andere_ebene():
    """Dieselbe Form auf anderer Ebene ist erlaubt, wenn sie genannt ist; auf derselben Ebene bleibt sie Doppelformel."""
    from kern.formeln import pruefe_alternativen
    a = _formel("probe_a", lambda M: math.sqrt(M) * 3)
    z = _formel("probe_z", lambda M: 3 * M ** 0.5, ebene="Z", gleiche_form="probe_a")
    assert W.doppelformeln({"probe_a": a, "probe_z": z}) == []
    assert pruefe_alternativen({"probe_a": a, "probe_z": z}) == []
    ohne = _formel("probe_z", lambda M: 3 * M ** 0.5, ebene="Z")
    assert len(W.doppelformeln({"probe_a": a, "probe_z": ohne})) == 1
    gleich = _formel("probe_l", lambda M: 3 * M ** 0.5, ebene="L", gleiche_form="probe_a")
    assert len(W.doppelformeln({"probe_a": a, "probe_l": gleich})) == 1
    assert any("derselben Ebene" in f for f in pruefe_alternativen({"probe_a": a, "probe_l": gleich}))
    assert "kac_rueckkehrzeit" in FORMELN["schleife_thermisch"].gleiche_form
