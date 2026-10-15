# -*- coding: utf-8 -*-
"""
Gesetze Klasse A: Grenzgesetze der Physik, Lehrbuchform (Vergleichsseite, von Hand) gegen die Faltungsform des Kerns
(Modellseite). Aufbau: GESETZE.neu(kennung, "A", Aussage, pruefung, Quelle, Ebene, Triade, …) mit den Feldern aus
gesetz.py; beide Seiten als lambda (Gegenprobe-Regel). Kennungen beginnen mit "A_".
Rechenpunkte ohne Eingabe (Sonnenmasse, 1 AE, Doppelstern 1,4 + 1,3 M☉ bei 1000 km) sind Stellen, an denen die
Identität gilt; sie gilt für jeden Wert. Körper und Bahnen aus Messungen kommen aus dem Messwert-Register.
"""
import math

import numpy as np
from scipy.optimize import minimize_scalar

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE

from .gesetz import GESETZE, ident

_g = GESETZE.neu
M = K.M_sun
GM = K.GM_sun
SIGMA_T = 8 * math.pi / 3 * K.r_e**2          # Thomson-Querschnitt (Lehrbuch)
DOPPEL_M1, DOPPEL_M2, DOPPEL_R = 1.4, 1.3, 1e6  # Rechenpunkt der Quadrupolleistung (M☉, M☉, m)


def _mw(k):
    return MESSWERTE[k].wert


_g("A_schwarzschild_radius", "A", "Schwarzschild-Radius 2GM/c² = r₀ = 6Δs",
   lambda p: ident(lambda: K.r0(M), lambda: 2 * GM / K.c**2), "Schwarzschild 1916", "L", "Zentrum")
_g("A_zeitdilatation", "A", "Zeitdilatation √(1 − 2GM/rc²) = √(1 − 6Δs/r)",
   lambda p: ident(lambda: K.takt(M, K.AU), lambda: math.sqrt(1 - 2 * GM / (K.AU * K.c**2))), "Einstein 1911/1916",
   "L", "Zentrum")
_g("A_lichtablenkung", "A", "Lichtablenkung 4GM/(c²b) = 12Δs/b (Sonnenrand)",
   lambda p: ident(lambda: K.ablenkung(M, K.R_sun), lambda: 4 * GM / (K.c**2 * K.R_sun)), "Einstein 1915", "L",
   "Wirbel")
_g("A_periheldrehung", "A", "Periheldrehung 6πGM/(c²a(1−e²)) = 18πΔs/(a(1−e²)) (Merkur)",
   lambda p: ident(lambda: K.perihel_pro_umlauf(M, _mw("MERKUR_A"), _mw("MERKUR_E")),
                   lambda: 6 * math.pi * GM / (K.c**2 * _mw("MERKUR_A") * (1 - _mw("MERKUR_E")**2))),
   "Einstein 1915", "L", "Wirbel", haengt_an=("MERKUR_A", "MERKUR_E"))
_g("A_photonensphaere", "A", "Photonensphäre 3GM/c² = 9Δs",
   lambda p: ident(lambda: K.photonensphaere(M), lambda: 3 * GM / K.c**2), "Schwarzschild", "L", "Zentrum")
_g("A_schattenradius", "A", "Schattenradius 3√3·GM/c² = 3√27·Δs",
   lambda p: ident(lambda: K.schatten_winkel(M, 1.0) / 2, lambda: 3 * math.sqrt(3) * GM / K.c**2), "Synge 1966", "L",
   "Zentrum")


def _isco_energie(j, mit=True):
    """ISCO als Minimum der Bahnenergie längs der Kreisbahnen (zweiter Weg, nicht die Bardeen-Formel), Einheit GM/c²."""
    s = 1 if mit else -1

    def E(x):
        return (x**1.5 - 2 * x**0.5 + s * j) / (x**0.75 * math.sqrt(x**1.5 - 3 * x**0.5 + s * 2 * j))
    xs = np.linspace(1.3 if mit else 3.3, 12, 4000)
    with np.errstate(invalid="ignore"):
        ev = np.array([E(x) if x**1.5 - 3 * x**0.5 + s * 2 * j > 0 else np.nan for x in xs])
    k = int(np.nanargmin(ev))
    return minimize_scalar(E, bounds=(xs[max(k - 5, 0)], xs[min(k + 5, len(xs) - 1)]), method="bounded",
                           options={"xatol": 1e-12}).x


for _j, _mit, _k in ((0.0, True, "A_isco_schwarzschild"), (None, True, "A_isco_kerr_mitlaeufig"),
                     (None, False, "A_isco_kerr_gegenlaeufig")):
    _g(_k, "A", f"ISCO Kerr ({'j = 0' if _j == 0.0 else 'j = 0,9'}, {'mitläufig' if _mit else 'gegenläufig'}) als "
       "Minimum der Bahnenergie",
       lambda p, j=_j, mit=_mit: ident(lambda: K.isco_kerr(EINGABEN["J_ISCO"].wert if j is None else j, mit),
                                       lambda: _isco_energie(EINGABEN["J_ISCO"].wert if j is None else j, mit), 1e-6),
       "Bardeen, Press, Teukolsky 1972", "L", "Wirbel")

_g("A_hawking_temperatur", "A", "Hawking-Temperatur ħc³/(8πGMk) = ħc/(24πkΔs)",
   lambda p: ident(lambda: K.hbar * K.c / (24 * math.pi * K.kB * K.ds(M)), lambda: K.hbar * K.c**3 / (8 * math.pi * GM * K.kB)),
   "Hawking 1974", "L", "Zentrum")
_g("A_bekenstein_hawking", "A", "Bekenstein-Hawking S = 4πkGM²/(ħc) = 36πk(Δs/ℓ_P)²",
   lambda p: ident(lambda: 36 * math.pi * (K.ds(M) / K.l_P)**2, lambda: 4 * math.pi * K.G * M**2 / (K.hbar * K.c)),
   "Bekenstein 1973", "L", "Zentrum")
_g("A_hawking_lebensdauer", "A", "Hawking-Lebensdauer 5120πG²M³/(ħc⁴) = 3840(S/k)Δs/c",
   lambda p: ident(lambda: K.lebensdauer_falte(M), lambda: 5120 * math.pi * K.G**2 * M**3 / (K.hbar * K.c**4)),
   "Page 1976", "L", "Zentrum")
_g("A_eddington", "A", "Eddington-Grenze 4πGMm_pc/σ_T = (9/2)m_pc³Δs/r_e²",
   lambda p: ident(lambda: K.eddington(M), lambda: 4 * math.pi * GM * K.m_p * K.c / SIGMA_T), "Eddington 1926", "L",
   "Zentrum")
_g("A_gibbons_hawking", "A", "Gibbons-Hawking ħH/(2πk) = ħc/(12πkΔs_H), Δs_H = c/6H (heute)",
   lambda p: ident(lambda: K.hbar * K.c / (12 * math.pi * K.kB * K.faltung_horizont(K.H0_si(p.H0))),
                   lambda: K.hbar * (p.H0 * 1e3 / K.Mpc) / (2 * math.pi * K.kB)),
   "Gibbons, Hawking 1977", "G", "Hintergrund")
_g("A_kepler", "A", "Kepler T² = 4π²a³/(GM) = 4π²a³/(3c²Δs) (1 AE)",
   lambda p: ident(lambda: 2 * math.pi * math.sqrt(K.AU**3 / (3 * K.c**2 * K.ds(M))),
                   lambda: 2 * math.pi * math.sqrt(K.AU**3 / GM)), "Kepler 1619, Newton", "L", "Wirbel")


def _quadrupol_lehrbuch():
    m1, m2 = DOPPEL_M1 * M, DOPPEL_M2 * M
    return 32 / 5 * K.G**4 * (m1 * m2)**2 * (m1 + m2) / (K.c**5 * DOPPEL_R**5)


_g("A_quadrupolleistung", "A", "Quadrupolleistung (32/5)G⁴(m₁m₂)²M/(c⁵r⁵) = (32/5)(c⁵/G)η_m²(3Δs/r)⁵",
   lambda p: ident(lambda: K.quadrupol_leistung((DOPPEL_M1 + DOPPEL_M2) * M,
                                                DOPPEL_M1 * DOPPEL_M2 / (DOPPEL_M1 + DOPPEL_M2)**2, DOPPEL_R),
                   _quadrupol_lehrbuch), "Peters, Mathews 1963", "L", "Wirbel")


def _peters_pdot():
    m1, m2, Pb, e = (_mw(k) for k in ("PSR1913_M1", "PSR1913_M2", "PSR1913_PB", "PSR1913_E"))
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e**2)**3.5
    return (-192 * math.pi * K.G**(5 / 3) / (5 * K.c**5) * (Pb / (2 * math.pi))**(-5 / 3) * fe
            * m1 * m2 * M**2 / ((m1 + m2) * M)**(1 / 3))


_PSR = ("PSR1913_M1", "PSR1913_M2", "PSR1913_PB", "PSR1913_E")
_g("A_hulse_taylor", "A", "Abnahme der Umlaufzeit Ṗ_b (Peters-Mathews) = Faltungsform (PSR B1913+16)",
   lambda p: ident(lambda: K.hulse_taylor(*(_mw(k) for k in _PSR)), _peters_pdot), "Peters 1964", "L", "Wirbel",
   haengt_an=_PSR)
_GPB = ("GPB_A", "GM_earth")
_g("A_geodaetisch", "A", "geodätische Präzession (3/2)(GM/c²r)ω = (9/2)(Δs/r)ω (Bahn von Gravity Probe B)",
   lambda p: ident(lambda: K.geodaetisch(K.M_earth, EINGABEN["GPB_A"].wert * 1e3),
                   lambda: 1.5 * K.GM_earth / (K.c**2 * EINGABEN["GPB_A"].wert * 1e3)
                   * math.sqrt(K.GM_earth / (EINGABEN["GPB_A"].wert * 1e3)**3)),
   "de Sitter 1916", "L", "Wirbel", haengt_an=_GPB)
_g("A_mitfuehrung", "A", "Mitführung im Polarorbit GJ/(2c²a³) = 3Δs·a_J·c/(2a³)",
   lambda p: ident(lambda: K.mitfuehrung_polarbahn(K.M_earth, K.verwindung(_mw("ERDE_J"), K.M_earth),
                                                   EINGABEN["GPB_A"].wert * 1e3),
                   lambda: K.G * _mw("ERDE_J") / (2 * K.c**2 * (EINGABEN["GPB_A"].wert * 1e3)**3)),
   "Lense, Thirring 1918; Schiff 1960", "L", "Wirbel", haengt_an=_GPB + ("ERDE_J",))
_g("A_sternmasse", "A", "Sternmasse m_Pl³/m_p² (Landau-Chandrasekhar-Skala)",
   lambda p: ident(lambda: K.sternmasse(), lambda: (K.hbar * K.c / K.G)**1.5 / K.m_p**2), "Landau 1932", "L",
   "Zentrum")
_g("A_alpha_G", "A", "α_G(M) = 3Δs/λ̄ = GM²/(ħc) (M = 10⁻²⁰ kg)",
   lambda p: ident(lambda: K.alpha_G(1e-20), lambda: K.G * 1e-40 / (K.hbar * K.c)), "Dirac, Eddington", "L", "Zentrum")


# ------------------------------------------------------------------ PPN und Newton-Grenze (07)
def _gamma_ppn():
    """γ aus der Lichtablenkung α = 2(1 + γ)GM/(c²b)."""
    return K.ablenkung(M, K.R_sun) * K.c**2 * K.R_sun / (2 * GM) - 1


def _beta_ppn():
    """β aus der Periheldrehung Δφ = (2 + 2γ − β)/3 · 6πGM/(c²a(1 − e²)), mit γ aus der Ablenkung."""
    a, e = _mw("MERKUR_A"), _mw("MERKUR_E")
    einheit = 6 * math.pi * GM / (K.c**2 * a * (1 - e * e))
    return 2 + 2 * _gamma_ppn() - 3 * K.perihel_pro_umlauf(M, a, e) / einheit


_g("A_ppn_gamma", "A", "γ = 1: Raumkrümmung je Masse wie in der ART (aus der Lichtablenkung der Faltungsform)",
   lambda p: ident(_gamma_ppn, lambda: 1.0, 1e-12), "Eddington, Robertson; Will 1993 (PPN)", "L", "Zentrum")
_g("A_ppn_beta", "A", "β = 1: Nichtlinearität der Überlagerung wie in der ART (aus der Periheldrehung, γ aus der Ablenkung)",
   lambda p: ident(_beta_ppn, lambda: 1.0, 1e-12), "Eddington, Robertson; Will 1993 (PPN)", "L", "Wirbel",
   haengt_an=("MERKUR_A", "MERKUR_E"))
_g("A_newton_erdbahn", "A", "Newton-Grenze: an der Erdbahn a = a_N (Zusatz e^{−s}, s = r*/r ≈ 7342, unter der Rechengenauigkeit)",
   lambda p: ident(lambda: K.a_regel(GM / K.AU**2, p.a0) / (GM / K.AU**2), lambda: 1.0, 1e-15), "Newton 1687", "L",
   "Wirbel")
