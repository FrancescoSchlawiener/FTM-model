# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 09 (Elektromagnetismus und Frequenzleiter): je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung,
(D, Stichwort), "Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen
mit "09_". Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Spin ½: J = ħ/2, also a_J = verwindung(ħ/2, M).
"""
import math

import kern as K
from kern.eingabe import EINGABEN as E
from kern.messwert import MESSWERTE
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "09"
D = datei(B)
_A = AUSSAGEN.neu
_P = PRUEFZAHLEN.neu
M_MW = E["M_MILCHSTRASSE"].wert * K.M_sun
M_NS, R_NS = E["NS_M"].wert * K.M_sun, E["NS_R"].wert * 1e3


def _aJ_halb(m):
    return K.verwindung(K.hbar / 2, m)


# ------------------------------------------------------------------ Rolle im Zustandsraum
_A("09_em_gegen_grav", (D, "α/α_G = 1,23555·10³⁶"), "Zahl")
_P("09_em_gegen_grav", B, "09_em_gegen_grav", "1,23555·10³⁶", lambda p: K.verhaeltnis_em_grav(K.m_p), "1", "α/α_G",
   haengt_an=("m_p", "alpha"))

# ------------------------------------------------------------------ Verwindung und magnetisches Moment
_A("09_verwindung_elektron", (D, "a_J = λ̄_e/2 = r_e/(2α) = 193,08 fm"), "Zahl")
_P("09_a_J_elektron", B, "09_verwindung_elektron", "193,08", lambda p: _aJ_halb(K.m_e) * 1e15, "fm", "a_J = λ̄_e/2",
   haengt_an=("m_e",))
_A("09_bohr_magneton", (D, "μ_B = eħ/(2m_e) = 9,274·10⁻²⁴ J/T"), "Zahl")
_P("09_bohr_magneton", B, "09_bohr_magneton", "9,274·10⁻²⁴", lambda p: K.moment_magnetisch(K.e_el, _aJ_halb(K.m_e)),
   "J/T", "μ_B = eħ/(2m_e)", haengt_an=("m_e",))
_A("09_spinzahl_teilchen", (D, "(Elektron 2,9·10⁴⁴, Proton 8,5·10³⁷)"), "Zahl")
_P("09_spinzahl_elektron", B, "09_spinzahl_teilchen", "2,9·10⁴⁴", lambda p: K.spinzahl(_aJ_halb(K.m_e), K.m_e), "1",
   "Elektron", haengt_an=("m_e",))
_P("09_spinzahl_proton", B, "09_spinzahl_teilchen", "8,5·10³⁷", lambda p: K.spinzahl(_aJ_halb(K.m_p), K.m_p), "1",
   "Proton", haengt_an=("m_p",))

# ------------------------------------------------------------------ Frequenzleiter (Tafel)
_A("09_frequenzleiter", (D, "| Ebene (Blickrichtung) | Frequenz | Wert |"), "Zahl")
_P("09_compton_takt_proton", B, "09_frequenzleiter", "1,425·10²⁴", lambda p: K.compton_takt(K.m_p), "rad/s",
   "ω₁ = c/λ̄_p", haengt_an=("m_p",))
_P("09_umlauf_wasserstoff", B, "09_frequenzleiter", "4,134·10¹⁶", lambda p: K.umlauf_wasserstoff(), "rad/s",
   "ω₂ = αc/a_B", haengt_an=("alpha", "a_B"))
_P("09_umlauf_wasserstoff_welle", B, "09_frequenzleiter", "45,6", lambda p: 2 * math.pi * K.c / K.umlauf_wasserstoff() * 1e9,
   "nm", "ω₂ = αc/a_B", haengt_an=("alpha", "a_B"))
_P("09_faltungsrate_proton", B, "09_frequenzleiter", "8,42·10⁻¹⁵", lambda p: K.faltungsrate(K.m_p), "s⁻¹",
   "Γ_G = α_G·ω₁", haengt_an=("m_p",))
_P("09_faltungstakt_proton", B, "09_frequenzleiter", "2,41·10⁶²", lambda p: K.faltungstakt(K.m_p), "s⁻¹",
   "c/(3Δs_p)", haengt_an=("m_p",))
_P("09_neutronenstern_takt", B, "09_frequenzleiter", "0,81", lambda p: K.takt(M_NS, R_NS), "1", "Neutronenstern",
   haengt_an=("NS_M", "NS_R", "GM_sun"))
_P("09_neutronenstern_umlauf", B, "09_frequenzleiter", "1,04·10⁴", lambda p: K.umlauf_kepler(M_NS, R_NS), "rad/s",
   "Neutronenstern", haengt_an=("NS_M", "NS_R", "GM_sun"))
_P("09_erdbahn", B, "09_frequenzleiter", "1,99·10⁻⁷", lambda p: K.umlauf_kepler(K.M_sun, K.AU), "rad/s", "| Erdbahn |",
   haengt_an=("GM_sun",))
_P("09_milchstrasse_umlauf", B, "09_frequenzleiter", "8,0·10⁻¹⁶",
   lambda p: K.v_kreis_punkt(M_MW, K.r_stern(M_MW, p.a0), p.a0) / K.r_stern(M_MW, p.a0), "rad/s",
   "Milchstraße (6·10¹⁰ M☉)", haengt_an=("M_MILCHSTRASSE", "GM_sun"))
_P("09_milchstrasse_umlaufzeit", B, "09_frequenzleiter", "248",
   lambda p: 2 * math.pi * K.r_stern(M_MW, p.a0) / K.v_kreis_punkt(M_MW, K.r_stern(M_MW, p.a0), p.a0) / K.yr / 1e6,
   "Mio. Jahre", "Milchstraße (6·10¹⁰ M☉)", haengt_an=("M_MILCHSTRASSE", "GM_sun"))

# ------------------------------------------------------------------ Kohärenz und Phasenseite
_A("09_dephasierung", (D, "ħ/kT = 2,46·10⁻¹⁴ s"), "Zahl")
_P("09_dephasierung_koerper", B, "09_dephasierung", "2,46·10⁻¹⁴", lambda p: K.dephasierzeit(E["T_KOERPER"].wert), "s",
   "ħ/kT", haengt_an=("T_KOERPER",))
_A("09_debye", (D, "die Debye-Länge in physiologischer Salzlösung 0,78 nm"), "Zahl")
_P("09_debye_physiologisch", B, "09_debye", "0,78",
   lambda p: K.debye_laenge(E["T_KOERPER"].wert, MESSWERTE["EPS_R_WASSER_37"].wert, E["I_PHYSIOLOGISCH"].wert * 1e3)
   * 1e9, "nm", "Debye-Länge", haengt_an=("T_KOERPER", "EPS_R_WASSER_37", "I_PHYSIOLOGISCH"))
_A("09_coulomb_unter_a0", (D, "1 kg mit 1 pC bei 110 V/m"), "Zahl")
_P("09_coulomb_feld_a0", B, "09_coulomb_unter_a0", "110", lambda p: p.a0 * E["PROBE_MASSE"].wert / E["PROBE_LADUNG"].wert,
   "V/m", "1 kg mit 1 pC", haengt_an=("PROBE_MASSE", "PROBE_LADUNG"))


def _draht_kanaele(p):
    """|q_N/q_σ| an den Arbeitspunkten der Tabelle (kälteste stationäre Lösung), an denen die Schleife merklich trägt."""
    Qm, grenze = E["DRAHT_QM"].wert, E["DRAHT_Q_MERKLICH"].wert
    aus = []
    for V, dfd in E["DRAHT_ARBEITSPUNKTE"].wert:
        q, q_s, _, q_n = K.draht_schleife(K.kupfer_huellkurve(dfd, V, Qm)['stationaer'][0], dfd, V, Qm)
        if abs(q) >= grenze:
            aus.append(abs(q_n / q_s))
    return aus


# Stoffwerte Kupfer [G] gehen über die Voreinstellung von draht() ein, nicht über einen Formelaufruf
_KD = ("DRAHT_ARBEITSPUNKTE", "DRAHT_Q_MERKLICH", "DRAHT_QM", "CU_RHO_EL", "CU_AL_R", "CU_E", "CU_AL_E", "CU_AL_TH",
       "CU_RHO", "CU_CP", "CU_K_TH", "DRAHT_A", "DRAHT_ELL", "DRAHT_B", "DRAHT_SIG0")
_A("09_draht_kanaele", (D, "über die Wärmespannung in der Struktur 200- bis 1900-mal stärker als über σ_el(T)"), "Zahl")
_P("09_draht_kanaele_min", B, "09_draht_kanaele", "200", lambda p: min(_draht_kanaele(p)), "1",
   "200- bis 1900-mal stärker", haengt_an=_KD)
_P("09_draht_kanaele_max", B, "09_draht_kanaele", "1900", lambda p: max(_draht_kanaele(p)), "1",
   "200- bis 1900-mal stärker", haengt_an=_KD)
