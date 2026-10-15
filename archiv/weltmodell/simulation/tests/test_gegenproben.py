# -*- coding: utf-8 -*-
"""Wächter: Gegenproben bestehen, und der unabhängige Weg ruft die geprüfte Formel nicht auf (derselbe Code zählt
nicht). Zählend: Formeln des Modells ohne Gegenprobe."""
import kern as K
import pytest

import waechter as W
from pruefung.gegenprobe import Gegenprobe


@pytest.mark.schwer
def test_bestand():
    assert W.gegenproben() == []
    assert isinstance(W.gegenproben_abdeckung(), list)


def test_negativ_derselbe_code():
    probe = Gegenprobe("probe", "r0", "zweiter_weg", [(K.M_sun,)], lambda M: 2 * K.r0(M) / 2)
    assert any("derselbe Code" in f for f in W.gegenproben({"probe": probe}))


def test_negativ_verschachtelt():
    def umweg(M):
        return K.r_stern(M) ** 2 / K.ell0() / 3 * 6   # r*²/ℓ₀ = 3Δs, also 6Δs = r₀, ohne r0 zu rufen
    probe = Gegenprobe("probe", "ds", "zweiter_weg", [(K.M_sun,)], lambda M: K.r0(M) / 6)
    assert any("derselbe Code" in f for f in W.gegenproben({"probe": probe}))
    ok = Gegenprobe("probe2", "r0", "zweiter_weg", [(K.M_sun,)], umweg, toleranz=1e-9)
    assert W.gegenproben({"probe2": ok}) == []


def test_negativ_abstand():
    probe = Gegenprobe("probe", "r0", "handwert", [(K.M_sun,)], lambda M: 2953.0, toleranz=1e-6, grund="Probe")
    assert any("Abstand" in f for f in W.gegenproben({"probe": probe}))


def test_mehrere_werte():
    """Formeln mit mehreren Werten: jeder Wert wird verglichen, Zahl und Form müssen gleich sein."""
    import math

    def lambert(d):
        from scipy.special import lambertw
        s = d + float(lambertw(-d * math.exp(-d)).real)
        return s, s * s, s**d / math.expm1(s)
    ok = Gegenprobe("probe", "bose_stufe_max", "zweiter_weg", [(3,)], lambert, toleranz=1e-10, grund="Probe")
    assert W.gegenproben({"probe": ok}) == []
    falsch = Gegenprobe("probe", "bose_stufe_max", "zweiter_weg", [(3,)],
                        lambda d: (lambert(d)[0], lambert(d)[1], 1.01 * lambert(d)[2]), toleranz=1e-6, grund="Probe")
    assert any("Abstand" in f for f in W.gegenproben({"probe": falsch}))
    kurz = Gegenprobe("probe", "bose_stufe_max", "zweiter_weg", [(3,)], lambda d: lambert(d)[:2], toleranz=1e-6,
                      grund="Probe")
    assert any("Abstand" in f for f in W.gegenproben({"probe": kurz}))


def _formel_paar():
    """Formeln, die r_stern_leiter transitiv aufruft."""
    from kern.formeln import FORMELN
    from kern.verfolgung import verfolge
    _, innen = verfolge(FORMELN["r_stern_leiter"].funktion, K.M_sun, K.A0)
    return innen


def test_negativ_teilt_aufgerufene_formel():
    """Ruft die geprüfte Formel transitiv die Modellformel F und der unabhängige Weg ebenfalls, ist die Probe nicht unabhängig."""
    innen = _formel_paar()
    assert "ds" in innen and "ell0" in innen
    teilend = Gegenprobe("probe", "r_stern_leiter", "zweiter_weg", [(K.M_sun, K.A0)],
                         lambda M, a0: (3 * K.ds(M) * K.ell0(a0)) ** 0.5)
    r = teilend.pruefe()
    assert r["teilt"] and not W._unabhaengig(teilend, None)
    frei = Gegenprobe("probe2", "r_stern_leiter", "zweiter_weg", [(K.M_sun, K.A0)],
                      lambda M, a0: (K.G * M / a0) ** 0.5)
    assert frei.pruefe()["teilt"] == [] and W._unabhaengig(frei, None)


def test_umkehr_zaehlt_nicht_fuer_abdeckung():
    from kern.formeln import FORMELN
    fk = "r0"
    unabh = Gegenprobe("p1", fk, "umkehr", [(K.M_sun,)], lambda M: 2 * K.G * M / K.c**2 * 3)
    zwei = Gegenprobe("p2", fk, "zweiter_weg", [(K.M_sun,)], lambda M: 2 * K.G * M / K.c**2 * 3)
    assert fk in W.gegenproben_abdeckung({"p1": unabh})
    assert fk not in W.gegenproben_abdeckung({"p2": zwei})
    z = W.gegenproben_abdeckung_zahlen({"p1": unabh})
    assert z["nur_umkehr"] == 1 and z["umkehrproben"] == 1 and z["ohne"] >= 1


def test_art_umkehr_erlaubt():
    Gegenprobe("p", "r0", "umkehr", [(K.M_sun,)], lambda M: 6 * K.G * M / K.c**2)
