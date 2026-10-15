# -*- coding: utf-8 -*-
"""[W] (Wahl): Eingaben mit wahl=True, Parametersätze und Gesetze mit wahl vererben [W] wie [G]; Negativproben zeigen,
dass Nicht-Wahlen und Konstanten keines tragen."""
import laden  # noqa: F401
from kern.eingabe import EINGABEN, Eingabe, leseprotokoll
from kern.herkunft_g import messwerte_g
from kern.herkunft_w import abhaengigkeit_w, eingaben_w, rechnung_w, sind_wahl
from kern.messwert import MESSWERTE
from kern.parametersatz import PARAMETERSAETZE
from pruefung.gesetz import GESETZE
from pruefung.pruefzahl import PRUEFZAHLEN, Pruefzahl

WAHLEN = {"DRAHT_Q_MERKLICH", "TORUS_K", "TORUS_DOMEGA", "ISS_HOEHE", "R_MILCHSTRASSE_HUELLE", "R_HAUFEN_HUELLE",
          "J_MILCHSTRASSE", "GPB_A", "T_KOERPER"}


def test_wahlen_des_bauplans():
    laden.lade_alles()
    assert set(eingaben_w()) == WAHLEN
    assert sind_wahl() == ["Upsilon_skaliert"]
    assert not EINGABEN["R_HAUFEN"].wahl and not EINGABEN["NS_M"].wahl     # Beispiele des Modelltexts sind keine Wahl


def test_pruefzahlen_erben_w():
    laden.lade_alles()
    assert PRUEFZAHLEN["06_rate_rastpunkt"].pruefe()["w"] == ["TORUS_DOMEGA", "TORUS_K"]
    assert PRUEFZAHLEN["09_draht_kanaele_min"].pruefe()["w"] == ["DRAHT_Q_MERKLICH"]
    assert "GPB_A" in PRUEFZAHLEN[next(k for k, p in PRUEFZAHLEN.items() if "GPB_A" in p.haengt_an)].pruefe()["w"]


def test_negativ_keine_wahl():
    laden.lade_alles()
    assert PRUEFZAHLEN["02_ds_erde"].pruefe()["w"] == []
    assert abhaengigkeit_w("R_HAUFEN") == set() and abhaengigkeit_w("GM_sun") == set()
    e = Eingabe("PROBE", "M", 1.0, "M☉", "Probe")                            # ohne wahl
    assert e.wahl is False
    assert Eingabe("PROBE2", "M", 1.0, "M☉", "Probe", wahl=True).wahl is True


def test_leseprotokoll_findet_gelesene_eingaben():
    laden.lade_alles()
    with leseprotokoll() as g:
        _ = EINGABEN["TORUS_K"].wert + EINGABEN["R_HAUFEN"].wert
    assert g == {"TORUS_K", "R_HAUFEN"}
    assert rechnung_w(g) == {"TORUS_K"}
    with leseprotokoll() as g2:
        pass
    assert g2 == set()                                                        # nach dem Block wird nichts mehr notiert
    _ = EINGABEN["TORUS_K"]
    assert g2 == set()


def test_parametersatz_und_gesetz_erben_w():
    laden.lade_alles()
    kopf = PRUEFZAHLEN["02_ds_erde"]
    probe = Pruefzahl("probe", kopf.baustein, kopf.aussage, "1,5", lambda p: 1.5, "1", kopf.modell_satz,
                      satz="Upsilon_skaliert")
    assert probe.pruefe()["w"] == ["Upsilon_skaliert"]
    assert Pruefzahl("probe", kopf.baustein, kopf.aussage, "1,5", lambda p: 1.5, "1", kopf.modell_satz).pruefe()["w"] == []
    assert PRUEFZAHLEN["03_sparc_a0_varianten_unten"].satz == "standard"
    assert PARAMETERSAETZE["Upsilon_skaliert"].ups_b_skaliert and not PARAMETERSAETZE["standard"].ups_b_skaliert
    tafel5 = [k for k, g in GESETZE.items() if g.wahl]
    assert tafel5 and all(GESETZE[k].quelle.startswith("Autor, Tafel 5") for k in tafel5)
    # [W] nur für Gedächtnis, Maxwell-Element und Operator; Poisson-Inversion und Additivität sind bestätigt
    ohne = {k for k, g in GESETZE.items() if g.quelle.startswith("Autor, Tafel 5") and not g.wahl}
    assert ohne == {"Z_poisson_inversion", "Z_faltung_additiv"}
    assert all(k.startswith(("Z_gedaechtnis", "Z_maxwell", "Z_bindungsoperator", "Z_spiel")) or "gedaechtnis" in k for k in tafel5), tafel5
    r = GESETZE["Z_bindungsoperator"].pruefe()
    assert r["w"] and "Tafel 5" in r["w"][0]
    assert GESETZE["Z_poynting_joule"].pruefe()["w"] == []


def test_g_nur_aus_messwerten():
    """[G] zählt nur Messwerte ohne Konstanten (CODATA, IAU); horizontwachstum und kreislauf_versiegen erben nichts
    über ℓ_P und G."""
    laden.lade_alles()
    gz = set(messwerte_g())
    assert not gz & {"G", "m_p", "m_e", "GM_sun", "R_sun", "alpha", "GM_earth"}
    assert gz and all(MESSWERTE[k].ist_g for k in gz)
    for k in ("03_horizontwachstum_heute", "04_kreislauf_versiegen"):
        g = set(PRUEFZAHLEN[k].pruefe()["g"])
        assert g <= gz and not g & {"G", "m_p", "m_e", "GM_sun", "R_sun"}, (k, g)
