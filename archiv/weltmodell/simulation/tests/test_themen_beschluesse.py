# -*- coding: utf-8 -*-
"""Themenregister und Beschlüsse: Bestand vollständig, Orte stehen im Modelltext, Urteile rechenbar, Auftragskopf."""
import pytest

import laden
import waechter as W
from forschung.beschluesse import ARTEN_OFFEN, OFFEN, VORSCHLAEGE, Vorschlag
from forschung.beschluss import BESCHLUESSE
from forschung.thema import THEMEN, auftrag
from pruefung.modelltext import lies

PFLICHT = ("dunkle_energie", "herkunft_q", "a0_schwelle", "zwerge", "warp", "naht_uebergabe", "linse_haufen",
           "mikroskala", "zellform", "zustandsraum_labor", "zeit_gedaechtnis", "feldform_scheiben")


def test_themen_vollstaendig():
    laden.lade_alles()
    assert all(k in THEMEN for k in PFLICHT)
    assert W.themen() == [] and W.themen_luecken() == []
    for t in THEMEN.values():
        for v in t.versuche:
            assert v.urteil          # jede Abhängigkeit ist bekannt (sonst KeyError)


def test_urteil_frei_gewaehlt():
    t = THEMEN["dunkle_energie"]
    v = next(v for v in t.versuche if v.kennung == "de_ratenschwelle")
    assert v.urteil == "offen, hängt an frei:Glättungsmasse"
    assert next(v for v in t.versuche if v.kennung == "de_phasenseite").urteil == "trägt nicht"


def test_auftrag_kopf():
    kopf = auftrag("dunkle_energie")
    t = THEMEN["dunkle_energie"]
    assert t.invariante.form in kopf and all(v.kennung in kopf for v in t.versuche)
    assert all(w in kopf for w in t.gesperrt)


def test_beschluesse_orte_im_modell():
    laden.lade_alles()
    assert W.beschluesse() == []
    for b in BESCHLUESSE.values():
        if b.richtung == "hinein":
            assert b.modell_ort[1] in lies(b.modell_ort[0]), b.kennung
    for k, (art, grund) in OFFEN.items():
        assert k in BESCHLUESSE and art in ARTEN_OFFEN and grund, k
    for v in VORSCHLAEGE.values():
        if v.modell_ort is not None:
            assert v.modell_ort[1] in lies(v.modell_ort[0]), v.kennung


def test_negativ_vorschlag():
    with pytest.raises(ValueError):
        Vorschlag("x", ("03-probe.md", "a"), "1", "2", "Laune", "Satz")
    with pytest.raises(ValueError):
        Vorschlag("x", ("03-probe.md", "a"), "1", "2", "Rundung", "")
