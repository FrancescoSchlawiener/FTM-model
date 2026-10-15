# -*- coding: utf-8 -*-
"""[G] nur aus Messwert-Einträgen: Naturkonstanten (CODATA, IAU) und die daraus abgeleiteten ℓ_P, m_Pl tragen kein [G];
ein Messwert der Art "G" (Planck, Sonnenort, GP-B …) vererbt [G] an jede Rechnung, die ihn nutzt."""
import laden  # noqa: F401
from kern import konstanten as K
from kern.herkunft_g import abhaengigkeit_g, code_wurzeln, formel_g, messwert_g, rechnung_g, symbol_wurzeln
from kern.parametersatz import PARAMETERSAETZE


def test_konstanten_ohne_g():
    for z in ("G", "m_p", "m_e", "α", "M☉", "GM☉", "ℓ_P", "m_Pl"):
        assert symbol_wurzeln(z) == set(), z
    for k in ("G", "m_p", "GM_sun", "R_sun", "alpha"):
        assert messwert_g(k) == set(), k


def test_keine_vererbung_ueber_planck():
    laden.lade_alles()
    for k in ("leiterzahl", "spiegel_masse", "r0", "ds"):
        assert formel_g(k) == set(), k
    assert rechnung_g(rechnung=lambda p: K.l_P**2 * K.m_Pl * K.G) == set()


def test_messwerte_vererben():
    laden.lade_alles()
    assert abhaengigkeit_g("DF2_SIGMA") == {"DF2_SIGMA"}
    assert abhaengigkeit_g("DF4_SIGMA") == {"DF4_SIGMA"}
    assert abhaengigkeit_g("CRATER_II_SIGMA") == set()  # Art literatur
    assert code_wurzeln(lambda p: K.H0_KMS * K.T_CMB) == set()  # beide Art literatur
    assert symbol_wurzeln("H₀") == set() and symbol_wurzeln("M⊕") == set()
    p = PARAMETERSAETZE["standard"]
    p.zuruecksetzen()
    _ = p.H0
    assert rechnung_g(satz=p) == set()
    p.zuruecksetzen()


def test_negativ_konstante_in_haengt_an():
    """Eine genannte Konstante (GM☉) gibt kein [G], ein genannter Messwert daneben schon."""
    assert rechnung_g(haengt_an=("GM_sun", "DF2_SIGMA")) == {"DF2_SIGMA"}


def test_art_literatur_ohne_g():
    from kern.messwert import Messwert
    m = Messwert("X_TEST", "H₀", 67.4, 0.5, "km/s/Mpc", "Planck 2018, Abstract", "literatur")
    assert not m.ist_g and m.datei is None
