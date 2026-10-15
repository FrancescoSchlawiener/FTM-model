# -*- coding: utf-8 -*-
"""
Gesetze Klasse E: Messwerte mit Quelle gegen das Modell (Abstand in Standardabweichungen). Aufbau wie gesetze_A.py;
Kennungen beginnen mit "E_". Der Modellwert ist immer gerechnet (Kern, Eingaben, Parametersatz), der Messwert kommt
nur aus dem Register; Messwerte mit Fehler bei 95 % gehen mit dem halben Fehler als 1σ ein.
"""
import math

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE

from .gesetz import GESETZE, mess, mess_wert, schranke

_g = GESETZE.neu
MAS_JAHR = K.yr * 180 / math.pi * 3600e3        # rad/s → mas/Jahr


def _mw(k):
    return MESSWERTE[k].wert


_g("E_a0_gegen_relation", "E", "Schwelle a₀ gegen die gemessene Relation g† (Fehler stat und sys)",
   lambda p: mess(lambda: p.a0, "A0_LIT"), "McGaugh, Lelli, Schombert 2016", "I", "Hintergrund")
_g("E_zelle_ohne_kreispaare", "E", "Zelle ohne Kreispaare: a₀ ≤ a₀(L = 0,985·d_LSS) (kubischer Torus)",
   lambda p: schranke(lambda: p.a0, lambda: K.a0_aus_zelle(K.KREISPAAR_GRENZE * K.d_LSS()), ("KREISPAAR_GRENZE",)),
   "COMPACT, Petersen et al. 2023 (JCAP 01, 030)", "IG", "Hintergrund",
   haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "OMEGA_R"))


def _n_eff(p):
    """N_eff = Standardwert + ΔN_eff des Vakuumteils bei der Rekombination: Ω_s(z) mit a₀(z) und H(z)."""
    a0z = p.a0 * K.a0_frueh_relativ(p.z)
    omega = K.omega_vakuum(a0z, p.H0 * K.E_hintergrund(p.z))
    return K.N_EFF_STD + K.delta_N_eff(omega)


_g("E_n_eff", "E", "N_eff mit dem Vakuumteil des Modells bei der Rekombination (a₀(z), H(z)): früh fast null",
   lambda p: mess(lambda: _n_eff(p), "N_EFF_MESS"), "Planck 2018 + BAO", "G", "Hintergrund",
   haengt_an=("N_EFF_STD", "OMEGA_R", "OMEGA_B_H2", "OMEGA_C_H2"), satz="z1090")


def _omega_k():
    """Krümmung der flachen Zelle: 1 − Ω_m − Ω_r − Ω_DE (Ausdehnungsgesetz als Standardmodell gelesen)."""
    return 1 - K.OMEGA_M - K.OMEGA_R - K.omega_de()


_g("E_omega_k", "E", "Krümmung Ω_K der flachen Zelle (Zelle atmet, Λ fest) gegen die Messung",
   lambda p: mess_wert(_omega_k, _mw("OMEGA_K_MESS"), MESSWERTE["OMEGA_K_MESS"].fehler / 2, ("OMEGA_K_MESS",)),
   "Planck 2018 + BAO (Fehler 95 %)", "G", "Hintergrund", haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "OMEGA_R"))


def _crater2(p):
    """σ_v von Crater II: Plummer-Kugel mit Jeans und der örtlichen Regel im Feld der Milchstraße (Punktmasse)."""
    M = _mw("CRATER2_M") * K.M_sun
    re = _mw("CRATER2_RE") * K.kpc
    ge = K.neigung(EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun, _mw("CRATER2_ABSTAND") * K.kpc)
    return K.sigma_plummer_jeans(M, re, ge, p.a0) / 1e3


_g("E_crater2", "E", "Crater II σ_v, örtliche Regel (Jeans, Außenfeld der Milchstraße)",
   lambda p: mess(lambda: _crater2(p), "CRATER_II_SIGMA"), "Caldwell et al. 2017", "L", "Hintergrund",
   haengt_an=("CRATER2_M", "CRATER2_RE", "CRATER2_ABSTAND", "M_MILCHSTRASSE"))
_GPB = ("GPB_A", "GM_earth")
_g("E_gpb_geodaetisch", "E", "Gravity Probe B geodätische Drehung [mas/Jahr]",
   lambda p: mess(lambda: K.geodaetisch(K.M_earth, EINGABEN["GPB_A"].wert * 1e3) * MAS_JAHR, "GPB_GEO"),
   "Everitt et al. 2011", "L", "Wirbel", haengt_an=_GPB)
_g("E_gpb_mitfuehrung", "E", "Gravity Probe B Mitführung [mas/Jahr] (Bahnmittel längs der Achse)",
   lambda p: mess(lambda: K.mitfuehrung_polarbahn(K.M_earth, K.verwindung(_mw("ERDE_J"), K.M_earth),
                                                  EINGABEN["GPB_A"].wert * 1e3) * MAS_JAHR, "GPB_LT"),
   "Everitt et al. 2011", "L", "Wirbel", haengt_an=_GPB + ("ERDE_J",))
_PSR = ("PSR1913_M1", "PSR1913_M2", "PSR1913_PB", "PSR1913_E")


def _ht_verhaeltnis():
    """Ṗ_b der Faltungsform geteilt durch die Formel der ART (Peters-Mathews von Hand)."""
    m1, m2, Pb, e = (_mw(k) for k in _PSR)
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e**2)**3.5
    art = (-192 * math.pi * K.G**(5 / 3) / (5 * K.c**5) * (Pb / (2 * math.pi))**(-5 / 3) * fe
           * m1 * m2 * K.M_sun**2 / ((m1 + m2) * K.M_sun)**(1 / 3))
    return K.hulse_taylor(m1, m2, Pb, e) / art


_g("E_hulse_taylor", "E", "Hulse-Taylor: gemessene zu erwarteter Abnahme der Umlaufzeit (Modell gleich ART)",
   lambda p: mess(_ht_verhaeltnis, "HT_VERHAELTNIS"), "Weisberg, Huang 2016", "L", "Wirbel", haengt_an=_PSR)
_g("E_schatten_m87", "E", "Schatten M87* (Masse und Abstand EHT) [μas]",
   lambda p: mess(lambda: K.schatten_winkel(_mw("M87_M") * K.M_sun, _mw("M87_D") * K.Mpc) * 180 / math.pi * 3600e6,
                  "EHT_RING"), "EHT 2019", "L", "Zentrum", haengt_an=("M87_M", "M87_D"))
