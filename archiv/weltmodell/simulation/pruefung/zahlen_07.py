# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 07 (Feldgleichung). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "07_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Alle Prüfzahlen laufen im Satz "standard" (a₀ = 1,1·10⁻¹⁰, H₀ heute, Υ_d = 0,5, gebunden, hier und heute).
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN
from .sonnenort import SONNENORT, gezeit

B = "07"
D = datei(B)
_A, _P = AUSSAGEN.neu, PRUEFZAHLEN.neu
MW_E = ("M_MILCHSTRASSE", "GM_sun")


def _e(k):
    return EINGABEN[k].wert


def _m(k):
    return MESSWERTE[k].wert


def _mw():
    """Milchstraße als Punktmasse [kg]."""
    return _e("M_MILCHSTRASSE") * K.M_sun


def _r8(p):
    """Sonnenkreis der Milchstraße als Punktmasse (Rechenpunkt 8 kpc) [m]; nur im Satz Sonnenkreis_Punktmasse."""
    if p.sonnenkreis != "punktmasse":
        raise ValueError(f"Sonnenkreis als Punktmasse nur im Satz mit sonnenkreis = punktmasse, nicht in {p.kennung}")
    return _e("R_ACHT_KPC") * K.kpc


# ------------------------------------------------------------------ Bose-Leiter (Beispiel, erster Eintrag)
_A("07_maxima_bose_leiter", (D, "Die Maxima liegen bei s → 0 (d = 1)"), "Zahl")
_P("07_maximum_d3", B, "07_maxima_bose_leiter", "7,96", lambda p: K.bose_stufe_max(3)[1], "a₀", "(d = 3)")
_P("07_maximum_d2_ort", B, "07_maxima_bose_leiter", "2,54", lambda p: K.bose_stufe_max(2)[1], "a₀",
   "größte Zusatzbeschleunigung")
_P("07_maximum_d2_wert", B, "07_maxima_bose_leiter", "0,648", lambda p: K.bose_stufe_max(2)[2], "a₀",
   "größte Zusatzbeschleunigung")
_P("07_maximum_d2_punktmasse", B, "07_maxima_bose_leiter", "0,628", lambda p: 1 / K.bose_stufe_max(2)[0], "r*",
   "Punktmasse bei")
_P("07_maximum_d2_milchstrasse", B, "07_maxima_bose_leiter", "5,5",
   lambda p: K.r_stern(_mw(), p.a0) / K.bose_stufe_max(2)[0] / K.kpc, "kpc", "Milchstraße", haengt_an=MW_E)

_A("07_leiter_bei_rstern", (D, "Bei r* sind alle Stufen gleich"), "Zahl")
_P("07_stufen_bei_rstern", B, "07_leiter_bei_rstern", "0,582", lambda p: K.bose_leiter(_e("S_EINS"), 2), "1",
   "Bei r* sind alle Stufen gleich", haengt_an=("S_EINS",))
_A("07_leiter_stufe3", (D, "Vakuumteil −Q̃(0) = 4π⁴/15"), "Zahl")
_P("07_leiter_nenner_15", B, "07_leiter_stufe3", "15", lambda p: math.pi**4 / K.bose_leiter_integral(3), "1",
   "Vakuumteil −Q̃(0) = 4π⁴/15")

# ------------------------------------------------------------------ Spitzenform klein (Δs/r)
_A("07_spitzenform_klein", (D, "In allen realen Systemen ist Δs/r klein"), "Zahl")
_P("07_ds_r_erdbahn", B, "07_spitzenform_klein", "3,3·10⁻⁹", lambda p: K.ds(K.M_sun) / K.AU, "1",
   "Erdbahn im Sonnentrichter", haengt_an=("GM_sun",))
_P("07_ds_r_milchstrasse", B, "07_spitzenform_klein", "1,2·10⁻⁷", lambda p: K.ds(_mw()) / _r8(p), "1",
   "Milchstraße bei 8 kpc", haengt_an=MW_E + ("R_ACHT_KPC",), satz="Sonnenkreis_Punktmasse")

# ------------------------------------------------------------------ Schwelle und Rückkehr
_A("07_schalter_schwelle", (D, "a₀ = 1,1·10⁻¹⁰ m/s²"), "Zahl")
_P("07_a0_standard", B, "07_schalter_schwelle", "1,1·10⁻¹⁰", lambda p: p.a0, "m/s²", "a₀ = 1,1·10⁻¹⁰")

_A("07_rate_theta_c", (D, "Dieselbe Rate a₀/c"), "Zahl")
_P("07_rate_a0_c", B, "07_rate_theta_c", "3,7·10⁻¹⁹", lambda p: K.verlust_gamma(p.a0), "s⁻¹", "Dieselbe Rate a₀/c")

_A("07_tau_D", (D, "τ_D hängt nicht von M ab"), "Zahl")
_P("07_tau_D_mrd", B, "07_tau_D", "86", lambda p: K.tau_D(p.a0) / K.Gyr, "Mrd. Jahre", "τ_D hängt nicht von M ab")

_A("07_periode_im_tau", (D, "Eine Periode P in Im τ"), "Zahl")
_P("07_periode_sonnenoberflaeche", B, "07_periode_im_tau", "2,5·10⁵",
   lambda p: K.s_(K.neigung(K.M_sun, K.R_sun), p.a0) / (2 * math.pi), "1", "Sonnenoberfläche: P >",
   haengt_an=("GM_sun", "R_sun"))

# ------------------------------------------------------------------ SPARC: Spektraldimension
def _sparc_punkte(p):
    """(g_N, g_obs) aller Punkte der SPARC-Auswahl mit Υ_d des Satzes und dem Gasdruck σ = 8 km/s."""
    S = K.SPARC_SATZ
    gN, g = [], []
    for d in K.sparc_tafel(S["q_max"], S["inc_min"]):
        vb2 = K.v_bar_quadrat(d["Vg"], d["Vd"], d["Vb"], p.Upsilon_d, S["ups_b"])
        m = (vb2 > 0) & (d["Vo"] > 0)
        gN.append(vb2[m] / d["R"][m])
        g.append(K.g_gemessen(d["R"][m], d["Vo"][m], d["Rd"], _e("SIG_GAS_8") * 1e3))
    return np.concatenate(gN), np.concatenate(g)


_A("07_spektraldimension_sparc", (D, "An den SPARC-Kurven gibt die Form der Regel d_s"), "Zahl")
_P("07_ds_sparc", B, "07_spektraldimension_sparc", "1,03", lambda p: K.ds_anpassung(*_sparc_punkte(p), p.a0), "1",
   "An den SPARC-Kurven gibt die Form der Regel d_s", haengt_an=("SIG_GAS_8",))
_P("07_punkte_unter_s_trenn", B, "07_spektraldimension_sparc", "16",
   lambda p: float(np.sum(K.s_(_sparc_punkte(p)[0], p.a0) < _e("S_TRENN"))), "Punkte", "dort liegen 16 Punkte",
   haengt_an=("SIG_GAS_8", "S_TRENN"))

# ------------------------------------------------------------------ Umlaufsinn
_A("07_umlaufsinn", (D, "Zählen beide Sinne voll, entsteht bei gleichem Tully-Fisher-Nullpunkt coth s"), "Zahl")
# coth s = 1 + 2n_B(2s) (Einsetzen in n_B)
_P("07_coth_schwelle", B, "07_umlaufsinn", "1,313", lambda p: 1 + 2 * K.n_B(2 * _e("S_EINS")), "1",
   "bei a_N = a₀ 1,313", haengt_an=("S_EINS",))


def _beide_sinne():
    """(1 + 2n_B)/(1 + n_B) bei s = 1: beide Sinne voll gegen einen Sinn."""
    nb = K.n_B(_e("S_EINS"))
    return (1 + 2 * nb) / (1 + nb)


_A("07_beide_sinne_tief", (D, "Beide Sinne gäben im tiefen Bereich χ → 2"), "Zahl")
_P("07_beide_sinne_faktor", B, "07_beide_sinne_tief", "1,37", lambda p: _beide_sinne(), "1", "das 1,37-Fache",
   haengt_an=("S_EINS",))
_P("07_milchstrasse_rstern_beide", B, "07_beide_sinne_tief", "253",
   lambda p: K.v_kreis_punkt(_mw(), K.r_stern(_mw(), p.a0), p.a0) * math.sqrt(_beide_sinne()) / 1e3, "km/s",
   "Milchstraße bei r* 253 statt 216 km/s", haengt_an=MW_E + ("S_EINS",))
_P("07_milchstrasse_rstern", B, "07_beide_sinne_tief", "216",
   lambda p: K.v_kreis_punkt(_mw(), K.r_stern(_mw(), p.a0), p.a0) / 1e3, "km/s",
   "Milchstraße bei r* 253 statt 216 km/s", haengt_an=MW_E)

# ------------------------------------------------------------------ Stationär: Rückkehr als Vorgang
def _anteil_vorgang(s, p):
    """Anteil des stationären q, der nach dem Weltalter stünde: ∫₀^t K(T)dT / e^{−s} (K mit Γ = s²a₀/c)."""
    Gam = K.verlust_gamma(s * s * p.a0)
    t = _m("WELTALTER") * K.Gyr
    return quad(lambda T: K.kern_zeit(T, Gam, p.a0), 0, t, epsabs=0, limit=400)[0] / math.exp(-s)


_A("07_rueckkehr_als_vorgang", (D, "liefe die Rückkehr als Vorgang in der Zeit ab"), "Zahl")
_P("07_vorgang_s1", B, "07_rueckkehr_als_vorgang", "19", lambda p: 100 * _anteil_vorgang(_e("S_EINS"), p), "%",
   "liefe die Rückkehr als Vorgang in der Zeit ab", haengt_an=("WELTALTER", "S_EINS"))
_P("07_vorgang_s03", B, "07_rueckkehr_als_vorgang", "10",
   lambda p: 100 * _anteil_vorgang(_e("S_NULL_KOMMA_DREI"), p), "%",
   "liefe die Rückkehr als Vorgang in der Zeit ab", haengt_an=("WELTALTER", "S_NULL_KOMMA_DREI"))

# ------------------------------------------------------------------ Eine Gleichung: Verluste unter einer Wurzel
def _theta_heute(p):
    """θ̃ = 3H₀ (Ausdehnungsrate des Hintergrunds heute); nur im mitatmenden Satz (zwischen den Strukturen)."""
    if p.zustand != "mitatmend":
        raise ValueError(f"θ̃ des Hintergrunds nur im mitatmenden Satz, nicht in {p.kennung}")
    return 3 * K.H0_si(p.H0)


_A("07_verluste_eine_wurzel", (D, "stehen unter einer Wurzel, nicht als Produkt"), "Zahl")
_P("07_theta_tau_D", B, "07_verluste_eine_wurzel", "17,9", lambda p: _theta_heute(p) * K.tau_D(p.a0), "1",
   "θτ_D = 17,9", satz="mitatmend_heute")
_P("07_wurzel_summe", B, "07_verluste_eine_wurzel", "1,30·10⁻²",
   lambda p: K.rueckkehr(0.0, K.verlust_gamma(p.a0), _theta_heute(p), p.a0), "1", "stehen unter einer Wurzel", satz="mitatmend_heute")
_P("07_wurzel_produkt", B, "07_verluste_eine_wurzel", "5,4·10⁻³",
   lambda p: K.rueckkehr(0.0, K.verlust_gamma(p.a0), 0.0, p.a0) * K.rueckkehr(0.0, 0.0, _theta_heute(p), p.a0), "1",
   "nicht als Produkt", satz="mitatmend_heute")

_A("07_schrittlaenge", (D, "die Schrittlänge fällt heraus"), "Zahl")
_P("07_q_milchstrasse_8kpc", B, "07_schrittlaenge", "0,33624", lambda p: K.q_(K.neigung(_mw(), _r8(p)), p.a0), "1",
   "Milchstraße bei 8 kpc: q = 0,33624", haengt_an=MW_E + ("R_ACHT_KPC",), satz="Sonnenkreis_Punktmasse")

_A("07_gewicht_zwischen_strukturen", (D, "Zwischen Strukturen sieht die Regel die Summe der Verluste"), "Zahl")
_P("07_gewicht_rueckkehr_heute", B, "07_gewicht_zwischen_strukturen", "2·10⁻⁴",
   lambda p: (1 - K.gewicht_w(_theta_heute(p), p.a0)) * K.n_B(math.sqrt(_theta_heute(p) * K.tau_D(p.a0))), "1",
   "(1 − w)(ν − 1) ≈ (1 − w)²/w", satz="mitatmend_heute")

# ------------------------------------------------------------------ Schluss bei Λ
def _rand_zehn_prozent():
    """r/(r*·Λ/ϖ_R), ab dem die bei Λ ≫ ϖ_R geschlossene Rückkehr (ϖ_R = r*) um den gewählten Anteil über der offenen
    liegt: ν_Λ/ν − 1 = ε."""
    L, eps = _e("LAM_ZU_VARPI_GROSS"), _e("ABWEICHUNG_ZEHN")
    s = brentq(lambda x: (1 - math.exp(-x)) / (1 - K.rueckkehr_geschlossen(x, L)) - 1 - eps, 1e-8, 1.0, xtol=1e-14)
    return 1 / (s * L)


_HS = ("LAM_ZU_VARPI_GROSS", "ABWEICHUNG_ZEHN")
_A("07_schluss_bei_lambda", (D, "Ist sie bei Λ ≫ ϖ_R geschlossen"), "Zahl")
_P("07_rand_065", B, "07_schluss_bei_lambda", "0,66", lambda p: _rand_zehn_prozent(), "r*·Λ/ϖ_R",
   "r ≈ 0,66·r*·Λ/ϖ_R", haengt_an=_HS)
_P("07_flach_bis_N", B, "07_schluss_bei_lambda", "1,5", lambda p: 1 / _rand_zehn_prozent(), "N",
   "Λ/ϖ_R ≳ 1,5·N", haengt_an=_HS)
_P("07_abweichung_eigenwert", B, "07_schluss_bei_lambda", "9,0", lambda p: 100 * K.abweichung_schluss(math.pi / 2), "%",
   "1/tanh(π/2) − 1 = 9,0 %")

_A("07_thermischer_propagator", (D, "sie trägt nur die Abweichung ab 0,657·Λ"), "Zahl")
_P("07_rand_0657", B, "07_thermischer_propagator", "0,657", lambda p: _rand_zehn_prozent(), "Λ",
   "sie trägt nur die Abweichung ab 0,657·Λ", haengt_an=_HS)

_A("07_netz_aller_falten", (D, "Ein gemeinsames Netz aller Falten der Länge Λ trüge die Rückkehrtiefen nicht"), "Zahl")

# ------------------------------------------------------------------ Offene Wahl des Imaginärteils
_A("07_wahl_imaginaerteil", (D, "Offene Wahl des Imaginärteils"), "Zahl")
_P("07_im_tau_erdbahn", B, "07_wahl_imaginaerteil", "3,1", lambda p: K.im_tau_aus_nome(K.ds(K.M_sun) / K.AU), "1",
   "an der Erdbahn bei 3,1", haengt_an=("GM_sun",))
_P("07_im_tau_milchstrasse", B, "07_wahl_imaginaerteil", "2,5", lambda p: K.im_tau_aus_nome(K.ds(_mw()) / _r8(p)), "1",
   "in der Milchstraße bei 2,5", haengt_an=MW_E + ("R_ACHT_KPC",), satz="Sonnenkreis_Punktmasse")

# ------------------------------------------------------------------ Auswertung: Toy Paper
_A("07_toy_paper", (D, "Die Form a = v₀²/r des Toy Papers"), "Zahl")
_P("07_masse_220", B, "07_toy_paper", "1,60·10¹¹", lambda p: K.faltungsmesser(_e("V_TOY_PAPER") * 1e3, p.a0) / K.ds(K.M_sun),
   "M☉", "220 km/s gehören zu", haengt_an=("V_TOY_PAPER", "GM_sun"))

# ------------------------------------------------------------------ Realteil: Umlaufphase und Resonanzen
_A("07_rosette", (D, "χ = 1 gibt Re τ_M = √2"), "Zahl")
_P("07_rosette_winkel", B, "07_rosette", "254,6", lambda p: 360 / math.sqrt(K.realteil_quadrat_punkt(_e("S_TIEF"))),
   "°", "Rosetten mit 254,6°", haengt_an=("S_TIEF",))

_A("07_resonanzen_punktmasse", (D, "örtliche Resonanzen bei rationalem κ/Ω"), "Zahl")
for _pk, _txt in ((9 / 8, "0,45"), (6 / 5, "0,68"), (5 / 4, "0,94"), (4 / 3, "2,1"), (7 / 5, "12")):
    _P(f"07_resonanz_rstern_{_txt.replace(',', '_')}", B, "07_resonanzen_punktmasse", _txt,
       lambda p, pk=_pk: K.radius_resonanz(pk, 1.0), "r*", f"bei {_txt} r*")

_A("07_resonanzen_milchstrasse", (D, "In der Milchstraße liegen diese Resonanzen bei"), "Zahl")
for _pk, _txt in ((6 / 5, "5,9"), (5 / 4, "8,2"), (9 / 7, "10,9"), (4 / 3, "18,0"), (7 / 5, "108")):
    _P(f"07_resonanz_mw_{_txt.replace(',', '_')}", B, "07_resonanzen_milchstrasse", _txt,
       lambda p, pk=_pk: K.radius_resonanz(pk, K.r_stern(_mw(), p.a0)) / K.kpc, "kpc", f"{_txt} kpc",
       haengt_an=MW_E)

# ------------------------------------------------------------------ Kern bei Bahnfrequenz
_A("07_kern_bahnfrequenz", (D, "Bei Bahnfrequenz trägt er"), "Zahl")
_P("07_kern_sonnenkreis", B, "07_kern_bahnfrequenz", "2·10⁻¹⁵",
   lambda p: K.kern_bahnfrequenz(K.v_kreis_punkt(_mw(), _r8(p), p.a0) / _r8(p), K.verlust_gamma(K.neigung(_mw(), _r8(p))),
                                 p.a0), "1", "am Sonnenkreis 2·10⁻¹⁵", haengt_an=MW_E + ("R_ACHT_KPC",), satz="Sonnenkreis_Punktmasse")

# ------------------------------------------------------------------ Klassische Tests
_A("07_zusatzterm_sonnensystem", (D, "Im Sonnensystem ist der Zusatzterm verschwindend"), "Zahl")
_P("07_s_erdbahn", B, "07_zusatzterm_sonnensystem", "7342", lambda p: K.r_stern(K.M_sun, p.a0) / K.AU, "1",
   "an der Erdbahn e^{−7342}", haengt_an=("GM_sun",))

# ------------------------------------------------------------------ Geneigte Ringe (örtliche Regel gegen Newton)
def _ring(rk, ik, p):
    M = _e("MN_MASSE") * K.M_sun
    a, b = _e("MN_A") * K.kpc, _e("MN_B") * K.kpc
    R, i = _e(rk) * K.kpc, math.radians(_e(ik))
    return K.praezession_ring(M, a, b, R, i, True, p.a0) / K.praezession_ring(M, a, b, R, i, False, p.a0)


_A("07_ringe_geneigt", (D, "Stark geneigte Ringe präzedieren absolut viel langsamer"), "Zahl")
for _rk, _ik, _txt, _satz in (("RING_R12", "RING_I2", "1,42", "bei 12 kpc 1,42"), ("RING_R12", "RING_I89", "1,66", "1,66 bei 89°"),
                              ("RING_R20", "RING_I2", "1,69", "bei 20 kpc 1,69"), ("RING_R20", "RING_I89", "1,91", "1,69 und 1,91"),
                              ("RING_R30", "RING_I2", "2,00", "bei 30 kpc 2,00"), ("RING_R30", "RING_I89", "2,19", "2,00 und 2,19")):
    _P(f"07_ring_{_rk[-2:]}_{_ik[6:]}", B, "07_ringe_geneigt", _txt, lambda p, rk=_rk, ik=_ik: _ring(rk, ik, p), "1",
       _satz, haengt_an=("MN_MASSE", "MN_A", "MN_B", _rk, _ik))

# ------------------------------------------------------------------ Hierarchie der Trichter: Gezeit am Sonnenort
_GEZEIT = SONNENORT


_A("07_glaettungswolke", (D, "Die Glättungswolke mit ℓ³ = GM/𝒯"), "Zahl")
_P("07_ell_erde", B, "07_glaettungswolke", "2780", lambda p: K.glaettungslaenge(K.M_earth, gezeit(p)) / K.AU, "AE",
   "Erde 2780 AE", haengt_an=_GEZEIT + ("GM_earth",), satz="Sonnenkreis_gemessen")
_P("07_ell_mond", B, "07_glaettungswolke", "641", lambda p: K.glaettungslaenge(_m("M_MOND"), gezeit(p)) / K.AU, "AE",
   "Mond 641 AE", haengt_an=_GEZEIT + ("M_MOND",), satz="Sonnenkreis_gemessen")
_P("07_umlauf_an_ell", B, "07_glaettungswolke", "84,5", lambda p: 2 * math.pi / math.sqrt(gezeit(p)) / K.yr / 1e6,
   "Mio. Jahre", "Umlauf an ℓ 84,5 Mio. Jahre", haengt_an=_GEZEIT, satz="Sonnenkreis_gemessen")

_A("07_schnitt_glaettung_schwelle", (D, "M× = a₀³/(G𝒯²)"), "Zahl")
_P("07_schnittmasse", B, "07_schnitt_glaettung_schwelle", "3,3·10⁸",
   lambda p: K.schnittmasse_glaettung(gezeit(p), p.a0) / K.M_sun, "M☉", "am Sonnenort 3,3·10⁸ M☉", haengt_an=_GEZEIT, satz="Sonnenkreis_gemessen")
_P("07_schnittradius", B, "07_schnitt_glaettung_schwelle", "642",
   lambda p: K.schnittradius_glaettung(gezeit(p), p.a0) / K.pc, "pc", "642 pc", haengt_an=_GEZEIT, satz="Sonnenkreis_gemessen")
_P("07_ell_zu_rstern_sonne", B, "07_schnitt_glaettung_schwelle", "26,2",
   lambda p: K.glaettung_zu_schwelle(K.M_sun, gezeit(p), p.a0), "1", "(Sonne 26,2)", haengt_an=_GEZEIT + ("GM_sun",), satz="Sonnenkreis_gemessen")

# Form der Rückkehr: Gitter-Rückkehr in drei Dimensionen, goldene Drehung
_A("07_watson", (D, "in drei ist ν ≤ 1,516 (Watsons Integral)"), "Zahl")
_P("07_watson_drei", B, "07_watson", "1,516", lambda p: K.watson_rueckkehr(_e("WATSON_D")), "1", "Watsons Integral",
   haengt_an=("WATSON_D",))
_GOLDEN = (math.sqrt(5) - 1) / 2   # Drehzahl der goldenen Drehung (im Satz genannt, kein Rechenpunkt)
_A("07_drehung_golden", (D, "golden, Zelle 1 %: 55, 89, 144"), "Zahl")
for _i, _z in enumerate(("55", "89", "144")):
    _P(f"07_drehung_golden_{_z}", B, "07_drehung_golden", _z,
       lambda p, i=_i: K.rueckkehrzeiten_drehung(_GOLDEN, _e("DREHUNG_ZELLE"))[i], "1", "Zelle 1 %",
       haengt_an=("DREHUNG_ZELLE",))

# Punkt 1 bei r*: s = 1, Im τ_M = s/2π; Δ gegen die Reihe Σqⁿ = q/(1 − q)
_A("07_delta_bei_rstern", (D, "bei r* nur noch 4,7·10⁻⁸ der Reihe"), "Zahl")
_P("07_delta_gegen_reihe", B, "07_delta_bei_rstern", "4,7·10⁻⁸",
   lambda p: K.delta_modul(1j * _e("S_EINS") / (2 * math.pi)) / (K.schleifensumme(math.exp(-_e("S_EINS"))) - 1), "1",
   "nur noch 4,7·10⁻⁸ der Reihe", haengt_an=("S_EINS",))
