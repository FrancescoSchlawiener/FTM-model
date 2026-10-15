# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 02 (Vom Teilchen bis zum Stern). Je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung,
(D, Stichwort), "Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py.
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Parametersatz überall "standard" (a₀ = 1,1·10⁻¹⁰, gebunden, hier und heute); a₀ geht nur über p.a0 ein.
"""
import math

from scipy.optimize import brentq

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN
from .sonnenort import SONNENORT, gezeit, gezeit_hill

B = "02"
D = datei(B)
_A = AUSSAGEN.neu
_P = PRUEFZAHLEN.neu


def _mw(k):
    return MESSWERTE[k].wert


def _ein(k):
    return EINGABEN[k].wert


_SONNENORT = SONNENORT + ("GM_sun",)


def _m_planeten():
    return sum(_mw(k) for k in ("MERKUR_M", "VENUS_M", "MARS_M", "JUPITER_M", "SATURN_M", "URANUS_M",
                                "NEPTUN_M")) + K.M_earth


_PLANETEN = ("MERKUR_M", "VENUS_M", "MARS_M", "JUPITER_M", "SATURN_M", "URANUS_M", "NEPTUN_M", "GM_earth", "GM_sun")


def _steilheit_h():
    """Bindungssteilheit der Wasserstoff-Hülle: Kreisbahn Ω² = a/r mit Ω = ω_H, r = a_B."""
    return K.umlauf_wasserstoff()**2 * K.a_B


def _m_ns():
    return _ein("NS_M") * K.M_sun


def _r_ns():
    return _ein("NS_R") * 1e3


def _m_rr():
    return _mw("BETEIGEUZE_M") * K.M_sun


def _r_rr():
    return _mw("BETEIGEUZE_R") * K.R_sun


def _delta_zuend(T, mu):
    """Tiefe δ, bei der die Kerntemperatur (n = 3) die Zündtemperatur T erreicht: T_c ist linear in δ = 3Δs/R."""
    n = _ein("N_STRAHLEND")
    return T / K.temperatur_kern(K.M_sun, K.R_sun, mu, n) * 3 * K.ds(K.M_sun) / K.R_sun


def _m_he_zuend(n):
    """Kleinste Kernmasse, deren höchste Temperatur (Tiefe gegen Zählregel der Elektronen) Helium zündet."""
    T, mu, mue = _mw("T_ZUEND_3A"), _ein("MU_HE"), _ein("MU_E")
    Mg = K.grenzmasse_zaehlregel(mue, n)
    return brentq(lambda m: K.temperatur_max(m * K.M_sun, mu, mue, n) - T, 0.01, 0.999 * Mg / K.M_sun, xtol=1e-6)


# ------------------------------------------------------------------ Tabelle
_A("02_tabelle_proton", (D, "| Proton | 4,14·10⁻⁵⁵ m | 4,9·10⁻⁴⁰ | 1,6·10⁻⁷ |"), "Zahl")
_P("02_ds_proton", B, "02_tabelle_proton", "4,14·10⁻⁵⁵", lambda p: K.ds(K.m_p), "m", "| Proton |",
   haengt_an=("m_p",))
_P("02_ds_R_proton", B, "02_tabelle_proton", "4,9·10⁻⁴⁰", lambda p: K.ds(K.m_p) / K.R_p, "1", "| Proton |",
   haengt_an=("m_p", "R_p"))
_P("02_g_proton", B, "02_tabelle_proton", "1,6·10⁻⁷", lambda p: K.neigung(K.m_p, K.R_p), "m/s²", "| Proton |",
   haengt_an=("m_p", "R_p"))

_A("02_tabelle_wasserstoff", (D, "| Wasserstoff-Atom | 4,14·10⁻⁵⁵ m | 7,8·10⁻⁴⁵ | Bindung 9·10²² |"), "Zahl")
_P("02_ds_R_wasserstoff", B, "02_tabelle_wasserstoff", "7,8·10⁻⁴⁵", lambda p: K.ds(K.m_p) / K.a_B, "1",
   "| Wasserstoff-Atom |", haengt_an=("m_p", "a_B"))
_P("02_bindung_wasserstoff", B, "02_tabelle_wasserstoff", "9·10²²", lambda p: _steilheit_h(), "m/s²",
   "| Wasserstoff-Atom |", haengt_an=("a_B",))

_A("02_tabelle_erde", (D, "| Erde | 1,478 mm | 2,3·10⁻¹⁰ | 9,82 |"), "Zahl")
_P("02_ds_erde", B, "02_tabelle_erde", "1,478", lambda p: K.ds(K.M_earth) * 1e3, "mm", "| Erde |",
   haengt_an=("GM_earth",))
_P("02_ds_R_erde", B, "02_tabelle_erde", "2,3·10⁻¹⁰", lambda p: K.ds(K.M_earth) / K.R_earth, "1", "| Erde |",
   haengt_an=("GM_earth", "R_earth"))
_P("02_g_erde", B, "02_tabelle_erde", "9,82", lambda p: K.neigung(K.M_earth, K.R_earth), "m/s²", "| Erde |",
   haengt_an=("GM_earth", "R_earth"))

_A("02_tabelle_sonne", (D, "| Sonne | 492,2 m | 7,1·10⁻⁷ | 274 |"), "Zahl")
_P("02_ds_sonne", B, "02_tabelle_sonne", "492,2", lambda p: K.ds(K.M_sun), "m", "| Sonne |", haengt_an=("GM_sun",))
_P("02_ds_R_sonne", B, "02_tabelle_sonne", "7,1·10⁻⁷", lambda p: K.ds(K.M_sun) / K.R_sun, "1", "| Sonne |",
   haengt_an=("GM_sun", "R_sun"))
_P("02_g_sonne", B, "02_tabelle_sonne", "274", lambda p: K.neigung(K.M_sun, K.R_sun), "m/s²", "| Sonne |",
   haengt_an=("GM_sun", "R_sun"))

_A("02_tabelle_neutronenstern", (D, "| Neutronenstern (1,4 M☉) | 689 m | 0,057 | 1,3·10¹² |"), "Zahl")
_P("02_ds_neutronenstern", B, "02_tabelle_neutronenstern", "689", lambda p: K.ds(_m_ns()), "m",
   "| Neutronenstern (1,4 M☉) |", haengt_an=("NS_M", "GM_sun"))
_P("02_ds_R_neutronenstern", B, "02_tabelle_neutronenstern", "0,057", lambda p: K.ds(_m_ns()) / _r_ns(), "1",
   "| Neutronenstern (1,4 M☉) |", haengt_an=("NS_M", "NS_R", "GM_sun"))
_P("02_g_neutronenstern", B, "02_tabelle_neutronenstern", "1,3·10¹²", lambda p: K.neigung(_m_ns(), _r_ns()), "m/s²",
   "| Neutronenstern (1,4 M☉) |", haengt_an=("NS_M", "NS_R", "GM_sun"))

_A("02_tabelle_roter_riese", (D, "| Roter Riese (18 M☉, 764 R☉) | 8,86 km | 1,7·10⁻⁸ | 0,0085 |"), "Zahl")
_P("02_ds_roter_riese", B, "02_tabelle_roter_riese", "8,86", lambda p: K.ds(_m_rr()) / 1e3, "km", "| Roter Riese",
   haengt_an=("BETEIGEUZE_M", "GM_sun"))
_P("02_ds_R_roter_riese", B, "02_tabelle_roter_riese", "1,7·10⁻⁸", lambda p: K.ds(_m_rr()) / _r_rr(), "1",
   "| Roter Riese", haengt_an=("BETEIGEUZE_M", "BETEIGEUZE_R", "GM_sun", "R_sun"))
_P("02_g_roter_riese", B, "02_tabelle_roter_riese", "0,0085", lambda p: K.neigung(_m_rr(), _r_rr()), "m/s²",
   "| Roter Riese", haengt_an=("BETEIGEUZE_M", "BETEIGEUZE_R", "GM_sun", "R_sun"))

# ------------------------------------------------------------------ Teilchen
_A("02_vakuum_flach", (D, "δ ≈ 10⁻³⁹ bei 1 fm neben einem Proton"), "Zahl")
# δ = 1 − dτ/dt verliert unter 10⁻¹⁶ alle Stellen (Auslöschung); hier die führende Ordnung δ = 3Δs/r
_P("02_delta_1fm", B, "02_vakuum_flach", "10⁻³⁹", lambda p: K.delta(K.m_p, _ein("ABSTAND_1FM")), "1",
   "bei 1 fm neben einem Proton", haengt_an=("m_p", "ABSTAND_1FM"), toleranz=0.5,
   toleranz_grund="Größenordnung (≈ 10ⁿ): ±50 %")

_A("02_alpha_em", (D, "α = r_e/λ̄_e = 1/137,036"), "Zahl")
_P("02_alpha_em_laengen", B, "02_alpha_em", "137,036", lambda p: 1 / (K.r_e / K.compton_red(K.m_e)), "1",
   "α = r_e/λ̄_e", haengt_an=("r_e", "m_e"))

# ------------------------------------------------------------------ Atom
_A("02_wasserstoff_huelle", (D, "Im Wasserstoff läuft die Hülle mit 6,58·10¹⁵ Hz"), "Zahl")
_P("02_umlauf_wasserstoff", B, "02_wasserstoff_huelle", "6,58·10¹⁵", lambda p: K.umlauf_wasserstoff() / (2 * math.pi),
   "Hz", "Im Wasserstoff läuft die Hülle", haengt_an=("alpha", "a_B"))
_P("02_bindung_wasserstoff_g", B, "02_wasserstoff_huelle", "9·10²¹", lambda p: _steilheit_h() / K.g_earth, "g⊕",
   "ihre Bindungssteilheit", haengt_an=("alpha", "a_B", "G_ERDE"))

_A("02_steilheit_faltung", (D, "würde eine Faltung von 3,8·10¹² kg verlangen"), "Zahl")
_P("02_steilheit_masse", B, "02_steilheit_faltung", "3,8·10¹²", lambda p: _steilheit_h() / K.neigung(1.0, K.a_B), "kg",
   "würde eine Faltung von", haengt_an=("alpha", "a_B"))

_A("02_schwelle_proton", (D, "Der Schwellenradius des Protons liegt bei r* = 32 fm"), "Zahl")
_P("02_r_stern_proton", B, "02_schwelle_proton", "32", lambda p: K.r_stern(K.m_p, p.a0) * 1e15, "fm",
   "Der Schwellenradius des Protons", haengt_an=("m_p",))
_P("02_bohr_zu_r_stern", B, "02_schwelle_proton", "1660", lambda p: K.a_B / K.r_stern(K.m_p, p.a0), "1",
   "a_B = 1660·r*", haengt_an=("m_p", "a_B"))

_A("02_glaettung_atom", (D, "in der Gezeit der Milchstraße bis zu seiner Glättungslänge von 2,7 mm"), "Zahl")
_P("02_glaettung_wasserstoff", B, "02_glaettung_atom", "2,7", lambda p: K.glaettungslaenge(K.m_p, gezeit(p)) * 1e3,
   "mm", "Glättungslänge von", haengt_an=("m_p",) + _SONNENORT, satz="Sonnenkreis_gemessen")

# ------------------------------------------------------------------ Planeten und Sonnensystem
_A("02_saturn", (D, "Saturn trägt die 95-fache Faltung der Erde"), "Zahl")
_P("02_saturn_faltung", B, "02_saturn", "95", lambda p: K.ds(_mw("SATURN_M")) / K.ds(K.M_earth), "1",
   "Saturn trägt die", haengt_an=("SATURN_M", "GM_earth"))
_P("02_saturn_kante", B, "02_saturn", "14",
   lambda p: 100 * (K.neigung(_mw("SATURN_M"), _mw("SATURN_R") * 1e3) / K.neigung(K.M_earth, K.R_earth) - 1), "%",
   "seine Kante ist nur", haengt_an=("SATURN_M", "SATURN_R", "GM_earth", "R_earth"))

_A("02_erde_bahn", (D, "(Erde: 29,78 km/s, ein Jahr)"), "Zahl")
_P("02_v_erde", B, "02_erde_bahn", "29,78", lambda p: K.v_kreis_punkt(K.M_sun, K.AU, p.a0) / 1e3, "km/s",
   "Erde: 29,78 km/s", haengt_an=("GM_sun",))

_A("02_faltung_system", (D, "Die Faltungen des Systems addieren sich zu 492,9 m, davon 99,87 % aus der Sonne"), "Zahl")
_P("02_ds_sonnensystem", B, "02_faltung_system", "492,9", lambda p: K.ds(K.M_sun + _m_planeten()), "m",
   "Die Faltungen des Systems", haengt_an=_PLANETEN)
_P("02_anteil_sonne", B, "02_faltung_system", "99,87",
   lambda p: 100 * K.ds(K.M_sun) / K.ds(K.M_sun + _m_planeten()), "%", "aus der Sonne", haengt_an=_PLANETEN)

_A("02_schwelle_sonne", (D, "Die Neigung fällt bei r* = √(GM☉/a₀) auf die Schwelle a₀ = 1,1·10⁻¹⁰ m/s²"), "Zahl")
_P("02_a0_schwelle", B, "02_schwelle_sonne", "1,1·10⁻¹⁰", lambda p: p.a0, "m/s²", "auf die Schwelle a₀")

_A("02_glaettung_sonne", (D, "Glättungslänge ℓ ≈ 0,93 pc (Hill-Radius 1,37 pc)"), "Zahl")
_P("02_glaettung_sonne_l", B, "02_glaettung_sonne", "0,93", lambda p: K.glaettungslaenge(K.M_sun, gezeit(p)) / K.pc,
   "pc", "Glättungslänge ℓ", haengt_an=_SONNENORT, satz="Sonnenkreis_gemessen")
_P("02_hill_sonne", B, "02_glaettung_sonne", "1,37",
   lambda p: K.glaettungslaenge(K.M_sun, gezeit_hill(p))
   / K.pc, "pc", "Hill-Radius", haengt_an=("OORT_A", "OORT_B", "GM_sun"), satz="Sonnenkreis_gemessen")

# ------------------------------------------------------------------ Sterne
_A("02_sternmasse", (D, "Δs_* = 912 m,     M_* = m_Pl³/m_p² = 1,853 M☉"), "Zahl")
_P("02_ds_stern", B, "02_sternmasse", "912", lambda p: K.ds(K.sternmasse()), "m", "Δs_* =", haengt_an=("m_p",))
_P("02_m_stern", B, "02_sternmasse", "1,853", lambda p: K.sternmasse() / K.M_sun, "M☉", "M_* = m_Pl³/m_p²",
   haengt_an=("m_p", "GM_sun"))

_A("02_carter_grenze", (D, "(beobachtet bei 1,2–1,4 M☉ ≈ 0,7 M_*)"), "Zahl")
_P("02_carter_grenze_mstern", B, "02_carter_grenze", "0,7",
   lambda p: sum(_mw("CARTER_GRENZE_BEOB")) / 2 * K.M_sun / K.sternmasse(), "M_*", "≈ 0,7 M_*",
   haengt_an=("CARTER_GRENZE_BEOB", "m_p", "GM_sun"))

_A("02_carter_faktor", (D, "Der Faktor 𝒞 = 2,944"), "Zahl")
_P("02_carter", B, "02_carter_faktor", "2,944", lambda p: K.carter_faktor(), "1", "Der Faktor 𝒞 =",
   haengt_an=("m_p", "m_e", "alpha"))
_P("02_carter_eins", B, "02_carter_faktor", "1402", lambda p: K.m_p / K.m_e * K.carter_faktor() ** -0.25, "1",
   "mit 𝒞 = 1 käme m_p/m_e", haengt_an=("m_p", "m_e", "alpha"))

_A("02_carter_sonne", (D, "L_Edd/L_* = 3,28·10⁴"), "Zahl")
_P("02_eddington_zu_L", B, "02_carter_sonne", "3,28·10⁴", lambda p: K.eddington(K.M_sun) / K.L_sun, "1",
   "L_Edd/L_* =", haengt_an=("GM_sun", "L_sun", "m_p", "r_e"))

_A("02_oberflaeche_alpha", (D, "δ_* = 3Δs_*/R_* ∝ (r_e/λ̄_e)² (Sonne 0,04·α²)"), "Zahl")
_P("02_delta_sonne_alpha", B, "02_oberflaeche_alpha", "0,04", lambda p: K.delta(K.M_sun, K.R_sun) / K.alpha**2, "α²",
   "Sonne 0,04·α²", haengt_an=("GM_sun", "R_sun", "alpha"))

_A("02_neutronenstern_huelle", (D, "rückt beim Neutronenstern die Hülle bis auf 2,9 r₀ an das Zentrum heran"), "Zahl")
_P("02_ns_huelle_r0", B, "02_neutronenstern_huelle", "2,9", lambda p: _r_ns() / K.r0(_m_ns()), "r₀",
   "Neutronenstern die Hülle", haengt_an=("NS_M", "NS_R", "GM_sun"))
_P("02_ns_takt", B, "02_neutronenstern_huelle", "19", lambda p: 100 * K.delta(_m_ns(), _r_ns()), "%",
   "der Takt ist dort um", haengt_an=("NS_M", "NS_R", "GM_sun"))

_A("02_tiefe_faktor", (D, "f = 0,54 (konvektiv, n = 1,5),  0,85 (strahlend, n = 3)"), "Zahl")
_P("02_f_konvektiv", B, "02_tiefe_faktor", "0,54", lambda p: K.tiefe_faktor(_ein("N_KONVEKTIV")), "1", "konvektiv",
   haengt_an=("N_KONVEKTIV",))
_P("02_f_strahlend", B, "02_tiefe_faktor", "0,85", lambda p: K.tiefe_faktor(_ein("N_STRAHLEND")), "1", "strahlend",
   haengt_an=("N_STRAHLEND",))

_A("02_kerntemperatur_sonne", (D, "(Sonne 1,2–1,7·10⁷ K)"), "Zahl")
_P("02_T_c_sonne_homogen", B, "02_kerntemperatur_sonne", "1,2",
   lambda p: K.temperatur_kern(K.M_sun, K.R_sun, _mw("MU_SONNE_HOMOGEN"), _ein("N_STRAHLEND")) / 1e7, "10⁷ K",
   "Sonne 1,2–1,7·10⁷ K", haengt_an=("GM_sun", "R_sun", "MU_SONNE_HOMOGEN", "N_STRAHLEND"))
_P("02_T_c_sonne_kern", B, "02_kerntemperatur_sonne", "1,7·10⁷",
   lambda p: K.temperatur_kern(K.M_sun, K.R_sun, _mw("MU_SONNE_KERN"), _ein("N_STRAHLEND")), "K",
   "Sonne 1,2–1,7·10⁷ K", haengt_an=("GM_sun", "R_sun", "MU_SONNE_KERN", "N_STRAHLEND"))

_A("02_brennstufen", (D, "von δ = 7·10⁻⁷ (pp) über 8·10⁻⁶ (3α) bis 1,7·10⁻⁴ (Si)"), "Zahl")
_P("02_delta_pp", B, "02_brennstufen", "7·10⁻⁷", lambda p: _delta_zuend(_mw("T_ZUEND_PP"), _mw("MU_SONNE_HOMOGEN")),
   "1", "(pp)", haengt_an=("T_ZUEND_PP", "MU_SONNE_HOMOGEN", "N_STRAHLEND"))
_P("02_delta_3a", B, "02_brennstufen", "8·10⁻⁶", lambda p: _delta_zuend(_mw("T_ZUEND_3A"), _ein("MU_HE")), "1",
   "(3α)", haengt_an=("T_ZUEND_3A", "MU_HE", "N_STRAHLEND"))
_P("02_delta_si", B, "02_brennstufen", "1,7·10⁻⁴", lambda p: _delta_zuend(_mw("T_ZUEND_SI"), _ein("MU_SI")), "1",
   "(Si)", haengt_an=("T_ZUEND_SI", "MU_SI", "N_STRAHLEND"))

_A("02_decke", (D, "kT ≈ 0,1·(μ/μ_e)·m_ec² ≈ 5·10⁸ K"), "Zahl")
_P("02_decke_faktor", B, "02_decke", "0,1", lambda p: K.decke_faktor(), "(μ/μ_e)·m_ec²", "kT ≈ 0,1·(μ/μ_e)·m_ec²")
_P("02_decke_kohlenstoff", B, "02_decke", "5·10⁸",
   lambda p: K.decke_faktor() * _ein("MU_CO") / _ein("MU_E") * K.m_e * K.c**2 / K.kB, "K", "Zündschwelle des Kohlenstoffs",
   haengt_an=("MU_CO", "MU_E", "m_e"))

_A("02_grenzmasse", (D, "ohne Wärme gibt dieselbe Bilanz die Grenzmasse 1,456 M☉"), "Zahl")
_P("02_grenzmasse_zaehlregel", B, "02_grenzmasse", "1,456",
   lambda p: K.grenzmasse_zaehlregel(_ein("MU_E"), _ein("N_STRAHLEND")) / K.M_sun, "M☉", "die Grenzmasse",
   haengt_an=("MU_E", "N_STRAHLEND", "GM_sun"))

_A("02_helium_zuendung", (D, "Helium zündet im entarteten Kern ab ≈ 0,41·M_Grenze (0,40–0,60 M☉ je nach Profil)"),
   "Zahl")
_P("02_helium_zu_grenze", B, "02_helium_zuendung", "0,41",
   lambda p: _m_he_zuend(_ein("N_STRAHLEND")) * K.M_sun / K.grenzmasse_zaehlregel(_ein("MU_E"), _ein("N_STRAHLEND")),
   "M_Grenze", "Helium zündet", haengt_an=("T_ZUEND_3A", "MU_HE", "MU_E", "N_STRAHLEND", "GM_sun"))
_P("02_helium_konvektiv", B, "02_helium_zuendung", "0,40", lambda p: _m_he_zuend(_ein("N_KONVEKTIV")), "M☉",
   "je nach Profil", haengt_an=("T_ZUEND_3A", "MU_HE", "MU_E", "N_KONVEKTIV", "GM_sun"))
_P("02_helium_strahlend", B, "02_helium_zuendung", "0,60", lambda p: _m_he_zuend(_ein("N_STRAHLEND")), "M☉",
   "je nach Profil", haengt_an=("T_ZUEND_3A", "MU_HE", "MU_E", "N_STRAHLEND", "GM_sun"))

# ------------------------------------------------------------------ Offen: Massenverhältnis und Nome
AUSSAGEN.neu("02_alpha_G", (D, "α_G (exponentiell, 5,906·10⁻³⁹)"), "Zahl")
PRUEFZAHLEN.neu("02_alpha_G_proton", B, "02_alpha_G", "5,906·10⁻³⁹", lambda p: K.alpha_G(K.m_p), "1",
                "α_G (exponentiell", haengt_an=("G", "m_p"))

_A("02_massenverhaeltnis", (D, "treffen 1836,15 und 1/137,036 nicht öfter"), "Zahl")
_P("02_mp_me", B, "02_massenverhaeltnis", "1836,15", lambda p: K.m_p / K.m_e, "1", "treffen 1836,15",
   haengt_an=("m_p", "m_e"))

_A("02_nome_tau", (D, "Die Nome-Form (Im τ_G = 14,0096, Im τ_α = 0,7831)"), "Zahl")
_P("02_im_tau_G", B, "02_nome_tau", "14,0096", lambda p: K.im_tau_aus_nome(K.alpha_G(K.m_p)), "1", "Im τ_G =",
   haengt_an=("G", "m_p"))
_P("02_im_tau_alpha", B, "02_nome_tau", "0,7831", lambda p: K.im_tau_aus_nome(K.alpha), "1", "Im τ_α =",
   haengt_an=("alpha",))

# ------------------------------------------------------------------ Starke Kopplung: Lauf und Nome (Skalen in GeV)
_QCD = ("ALPHA_S_MZ", "M_Z_GEV", "M_CHARM", "M_BOTTOM", "M_TOP")


def _schwellen():
    return tuple(_mw(k) for k in ("M_CHARM", "M_BOTTOM", "M_TOP"))


def _m_pl_gev():
    return K.m_Pl * K.c**2 / K.e_el / 1e9


def _alpha_s_planck(schleifen):
    return K.alpha_s_lauf(_mw("ALPHA_S_MZ"), _mw("M_Z_GEV"), _m_pl_gev(), _schwellen(), schleifen)


def _nome_planck(p):
    """|q̃| = Λ₃/M_Pl: Λ₃ = m_c·q_QCD(α_s(m_c), n_f = 3); die Nome-Form ist die Lösung einer Schleife, der Lauf deshalb
    ebenfalls mit einer Schleife."""
    m_c = _mw("M_CHARM")
    a_c = K.alpha_s_lauf(_mw("ALPHA_S_MZ"), _mw("M_Z_GEV"), m_c, _schwellen(), 1)
    return m_c * K.nome_qcd(a_c, 3) / _m_pl_gev()


def _empfindlichkeit(vz):
    """Änderung von m_p in %, wenn α_s(M_Pl) um vz·1 % wächst: m_p ∝ Skala, an der die abwärts laufende Kopplung
    (zwei Schleifen) denselben starken Wert erreicht (Dynamik fest)."""
    a, m_pl, stark, sw = _alpha_s_planck(2), _m_pl_gev(), _ein("ALPHA_S_STARK"), _schwellen()
    mu = K.skala_bei_kopplung(a, m_pl, stark, sw, 2)
    return 100 * (K.skala_bei_kopplung(a * (1 + vz * _ein("ALPHA_S_AENDERUNG")), m_pl, stark, sw, 2) / mu - 1)


_A("02_nome_starke_kopplung", (D, "1/α_s ≈ 53 an der Planck-Skala, |q̃| = 1,2·10⁻²⁰"), "Zahl")
_P("02_alpha_s_planck", B, "02_nome_starke_kopplung", "53", lambda p: 1 / _alpha_s_planck(2), "1", "1/α_s ≈ 53",
   haengt_an=_QCD)
_P("02_nome_qcd_planck", B, "02_nome_starke_kopplung", "1,2·10⁻²⁰", _nome_planck, "1", "|q̃| = 1,2·10⁻²⁰",
   haengt_an=_QCD)
_P("02_alpha_G_groessenordnung", B, "02_nome_starke_kopplung", "5,9·10⁻³⁹", lambda p: K.alpha_G(K.m_p), "1",
   "und so klein wie 5,9·10⁻³⁹", haengt_an=("G", "m_p"), toleranz=0.5, toleranz_grund="Größenordnung (≈ 10ⁿ): ±50 %")
_A("02_empfindlichkeit_mp", (D, "ändert m_p um +41 % bzw. −30 %"), "Zahl")
_P("02_empfindlichkeit_plus", B, "02_empfindlichkeit_mp", "41", lambda p: _empfindlichkeit(+1), "%",
   "ändert m_p um", haengt_an=_QCD + ("ALPHA_S_AENDERUNG", "ALPHA_S_STARK"))
_P("02_empfindlichkeit_minus", B, "02_empfindlichkeit_mp", "−30", lambda p: _empfindlichkeit(-1), "%",
   "ändert m_p um", haengt_an=_QCD + ("ALPHA_S_AENDERUNG", "ALPHA_S_STARK"))

# ------------------------------------------------------------------ Ladung als Windung
_A("02_dirac", (D, "Verhältnis magnetischer zu elektrischer Elementarladung 1/(2α) = 68,5"), "Zahl")
_P("02_dirac_verhaeltnis", B, "02_dirac", "68,5", lambda p: K.dirac_verhaeltnis(), "1", "1/(2α) = 68,5",
   haengt_an=("alpha",))
_A("02_ladungskreis", (D, "Sein Radius ist 2ℓ_P/√α = 23,4·ℓ_P, ein Umlauf trüge 5·10¹⁷ GeV"), "Zahl")
_P("02_ladungskreis_radius", B, "02_ladungskreis", "23,4", lambda p: K.radius_ladungskreis() / K.l_P, "ℓ_P",
   "Sein Radius ist", haengt_an=("alpha",))
_P("02_ladungskreis_umlauf", B, "02_ladungskreis", "5·10¹⁷", lambda p: K.energie_ladungskreis() / K.e_el / 1e9, "GeV",
   "ein Umlauf trüge", haengt_an=("alpha",))

# ------------------------------------------------------------------ Zustand des Stoffs (PubChem)
def _eisen_zustand(p):
    """δ der Erde (3Δs/R) gegen die Zustandsschwelle von Eisen bei seiner Schmelztemperatur."""
    fe = K.pubchem(_ein("Z_EISEN"))
    return K.delta(K.M_earth, K.R_earth) / K.zustandsschwelle(fe["T_schm"], fe["masse_u"] * K.m_u)


_A("02_zustand_eisen", (D, "Die Erde liegt für Eisen bei dessen Schmelztemperatur 116-fach darüber"), "Zahl")
_P("02_zustand_eisen_erde", B, "02_zustand_eisen", "116", _eisen_zustand, "1", "Die Erde liegt für Eisen",
   haengt_an=("Z_EISEN", "GM_earth", "R_earth"))

# ------------------------------------------------------------------ Oberfläche und Druck der Sterne
_A("02_f_s_sonne", (D, "Sonne: f_s = 0,0183"), "Zahl")
_P("02_f_s_sonne", B, "02_f_s_sonne", "0,0183", lambda p: K.f_saha(K.TEFF_sun), "1", "Sonne: f_s = 0,0183",
   haengt_an=("TEFF_sun",))
_A("02_gleichstand_druck", (D, "Strahlungsdruck und Gasdruck sind bei 72,6 M_* = 134,6 M☉ gleich"), "Zahl")
_P("02_gleichstand_druck_mstern", B, "02_gleichstand_druck", "72,6",
   lambda p: K.eddington_quartik(_ein("BETA_GLEICHSTAND"), _mw("MU_SONNE_HOMOGEN")) / K.sternmasse(), "M_*",
   "Strahlungsdruck und Gasdruck", haengt_an=("BETA_GLEICHSTAND", "MU_SONNE_HOMOGEN"))
_P("02_gleichstand_druck_msonne", B, "02_gleichstand_druck", "134,6",
   lambda p: K.eddington_quartik(_ein("BETA_GLEICHSTAND"), _mw("MU_SONNE_HOMOGEN")) / K.M_sun, "M☉",
   "Strahlungsdruck und Gasdruck", haengt_an=("BETA_GLEICHSTAND", "MU_SONNE_HOMOGEN", "GM_sun"))
