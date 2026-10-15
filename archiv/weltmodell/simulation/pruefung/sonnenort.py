# -*- coding: utf-8 -*-
"""
Sonnenort in Prüfzahlen: die Messwerte des Sonnenorts und die Gezeit 𝒯 des geglätteten Felds der Milchstraße daraus
(Einsetzen in die Kernformel gezeit_sonnenort). Nur im Parametersatz mit sonnenkreis = "gemessen".
"""
import kern as K
from kern.messwert import MESSWERTE

SONNENORT = ("V0_SONNE", "R0_SONNE", "OORT_A", "OORT_B", "RHO_B_SONNE")


def _mw(k):
    return MESSWERTE[k].wert


def gezeit(p):
    """𝒯 [s⁻²] am Sonnenort aus v₀, R₀, Oort A und B, ρ_b im Satz p."""
    if p.sonnenkreis != "gemessen":
        raise ValueError(f"Gezeit am Sonnenort nur im Satz mit gemessenem Sonnenort, nicht in {p.kennung}")
    kk = 1e3 / K.kpc
    return K.gezeit_sonnenort(_mw("V0_SONNE") * 1e3, _mw("R0_SONNE") * K.kpc, _mw("OORT_A") * kk, _mw("OORT_B") * kk,
                              _mw("RHO_B_SONNE") * K.M_sun / K.pc**3, p.a0)


def gezeit_hill(p):
    """𝒯_H [s⁻²] der Hill-Gezeit aus Oort A und B im Satz p."""
    if p.sonnenkreis != "gemessen":
        raise ValueError(f"Hill-Gezeit nur im Satz mit gemessenem Sonnenort, nicht in {p.kennung}")
    return K.gezeit_hill(_mw("OORT_A") * 1e3 / K.kpc, _mw("OORT_B") * 1e3 / K.kpc)
