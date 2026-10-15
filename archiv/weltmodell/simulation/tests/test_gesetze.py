# -*- coding: utf-8 -*-
"""Gesetze: alle laufen zusammen im Budget von 40 s (Überschreitung nur Warnung); die Ergebnisse außer OK sind genau die bekannten; die
Gegenprobe-Regel hält (die Vergleichsseite ruft nicht nur die Formeln der Modellseite); jedes Gesetz trägt Quelle,
Ebene und Triadenrolle; kein Messwert steht als Literal in den Gesetzen."""
import io
import os
import time
import warnings
import tokenize

import pytest

import kern as K
import laden
import waechter
from kern.messwert import MESSWERTE
from pruefung.gesetz import GESETZE, Gesetz, ident, regel_gegenprobe

GRENZE_S = 40.0
# Gesetze mit einem Ergebnis außer OK (Stand der Übertragung aus ftm_gesetze.py: 175/2/1/1)
ERWARTET = {"B_flaechendichte_beta": "SPANNUNG", "B_flaechendichte_lelli_fit": "SPANNUNG",
            "B_crater2_schaetzer": "HINWEIS", "E_crater2": "WIDERSPRUCH"}
# Zahlen in den Gesetzen, die zufällig einem Messwert gleichen, aber Rechenpunkte sind (mit Grund)
RECHENPUNKTE = {0.06: "Kopplung bzw. Schrittweite der Testfälle (gleich dem Messwert M_NU nur der Zahl nach)",
                1.8: "Antriebsamplitude der Kapitza-Probe (gleich dem Formparameter MCMILLAN_BULGE nur der Zahl nach)",
                0.9: "Koeffizient β₂ der Koexistenz zweier Zyklen bzw. Spin j = 0,9 (gleich MCMILLAN_DICK nur der Zahl nach)",
                0.3: "Testpunkte x = 0,3 (Amplitude, Nome, Skala) der Identitäten, gleich dem Messwert ALBEDO_ERDE nur der Zahl nach",
                0.4: "Hernquist-Skala b_H = 0,4 r* (Gesetz gilt für jede Skala)",
                1.4: "Doppelstern 1,4 + 1,3 M☉ der Quadrupolleistung (Identität gilt für jede Masse)",
                1.2: "Koeffizient β₁ der Koexistenz zweier Zyklen",
                0.04: "Toleranz bzw. Kopplung k der Einrastzunge",
                0.62: "Ort des zweiten Körpers auf dem Gittertorus 1 × 1",
                0.12: "Radius des geschlossenen Wegs auf dem Gittertorus",
                3e9: "Probemasse der Additivität (M☉)",
                42: "Exponent der Massenstichprobe 10⁴² kg",
                80: "Grenze der Nullstellensuche in log₁₀ M",
                20000: "Schrittzahl einer Simulation"}


@pytest.fixture(scope="module")
def ergebnisse():
    laden.lade_alles()
    t = time.time()
    erg = {k: g.pruefe() for k, g in GESETZE.items()}
    return erg, time.time() - t


@pytest.mark.schwer
def test_laufzeit(ergebnisse):
    """Laufzeit ist Warnung mit Budget, kein Testabbruch (hängt an der Last der Maschine)."""
    w = waechter.laufzeit_warnung(ergebnisse[1], GRENZE_S, "Gesetze")
    if w:
        warnings.warn(w)


def test_laufzeit_warnung_kein_abbruch():
    assert waechter.laufzeit_warnung(10.0, 40.0) is None
    assert "Budget von 40 s" in waechter.laufzeit_warnung(50.0, 40.0, "Gesetze")


@pytest.mark.schwer
def test_ergebnisse(ergebnisse):
    erg, _ = ergebnisse
    anders = {k: r["ergebnis"] for k, r in erg.items() if r["ergebnis"] != "OK"}
    assert anders == ERWARTET


@pytest.mark.schwer
def test_klassen_und_zahl(ergebnisse):
    erg, _ = ergebnisse
    zahl = {kl: sum(1 for k in erg if GESETZE[k].klasse == kl) for kl in "ABCZME"}
    assert zahl == {"A": 24, "B": 12, "C": 16, "Z": 49, "M": 124, "E": 9}


@pytest.mark.schwer
def test_gegenprobe_regel(ergebnisse):
    assert regel_gegenprobe(ergebnisse[0]) == []


def test_gegenprobe_regel_negativ():
    """Beide Seiten rufen r0: derselbe Code. Eine Zahl ohne Seite: keine Verfolgung."""
    laden.lade_alles()
    g = Gesetz("A_probe", "A", "Probe", lambda p: ident(lambda: K.r0(K.M_sun), lambda: 2 * K.r0(K.M_sun) / 2),
               "Probe", "L", "Zentrum")
    assert any("derselbe Code" in f for f in regel_gegenprobe({"A_probe": g.pruefe()}))
    g = Gesetz("A_probe2", "A", "Probe", lambda p: ident(K.r0(K.M_sun), 2 * K.G * K.M_sun / K.c**2), "Probe", "L",
               "Zentrum")
    assert any("keine Verfolgung" in f for f in regel_gegenprobe({"A_probe2": g.pruefe()}))
    g = Gesetz("A_probe3", "A", "Probe", lambda p: ident(lambda: K.r0(K.M_sun), lambda: 2 * K.G * K.M_sun / K.c**2),
               "Probe", "L", "Zentrum")
    assert regel_gegenprobe({"A_probe3": g.pruefe()}) == []


def test_felder_pflicht():
    laden.lade_alles()
    for fehlt in (dict(triade=None), dict(ebene="X"), dict(quelle="")):
        kw = dict(kennung="A_x", klasse="A", aussage="x", pruefung=lambda p: None, quelle="q", ebene="L",
                  triade="Zentrum")
        kw.update(fehlt)
        with pytest.raises((ValueError, KeyError)):
            Gesetz(**kw)
    with pytest.raises(ValueError):
        Gesetz("B_x", "A", "x", lambda p: None, "q", "L", "Zentrum")      # Kennung passt nicht zur Klasse
    for g in GESETZE.values():
        assert g.quelle and g.ebene and g.triade


def _messwert_zahlen():
    z = set()
    for m in MESSWERTE.values():
        for w in (m.wert if isinstance(m.wert, tuple) else (m.wert,)):
            z.add(float(w))
    return z


def test_keine_messwert_literale():
    laden.lade_alles()
    werte = _messwert_zahlen()
    ordner = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pruefung")
    funde = []
    for d in sorted(os.listdir(ordner)):
        if not (d.startswith("gesetze_") and d.endswith(".py")):
            continue
        with open(os.path.join(ordner, d), encoding="utf-8") as f:
            for tok in tokenize.generate_tokens(io.StringIO(f.read()).readline):
                if tok.type != tokenize.NUMBER:
                    continue
                try:
                    x = float(tok.string.replace("_", ""))
                except ValueError:
                    continue                                     # imaginäre Zahl (1j)
                if x in werte and x not in RECHENPUNKTE and not (2 * x == int(2 * x) and abs(x) <= 20):
                    funde.append(f"{d}:{tok.start[0]} {tok.string}")
    assert funde == []


def test_messwert_literal_negativ():
    """Crater II 2,7 als Literal fiele auf."""
    assert 2.7 in _messwert_zahlen()
