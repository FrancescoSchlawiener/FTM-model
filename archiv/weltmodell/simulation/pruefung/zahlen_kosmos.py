# -*- coding: utf-8 -*-
"""
Prüfzahlen der Gruppe Kosmos (Bausteine 03 und 04), gerechnet aus den Kernformeln, die für offene Zahlen hinzukamen:
Form und Drehungen der Zelle, kosmische Streuung, Gang der Schwelle aus der Strahlung, Störungen bei der Rekombination,
Horizontwachstum und Kreislauf, Weltzeit (Einsatz der Rückgabe), Linse jenseits der Naht, Rest im Haufen.
Aufbau wie zahlen_01.py; Messwerte und Eingaben nur über ihre Register, genannt in haengt_an. Satz je Prüfzahl:
Größen des Ganzen im Satz "mitatmend_heute", die Zelle (Formgröße) und gebundene Bereiche im Satz "standard", die
Störung bei der Rekombination im Satz "z1090".
"""
import cmath
import math

import kern as K
from kern.eingabe import EINGABEN as E
from kern.messwert import MESSWERTE as M
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

_A = AUSSAGEN.neu
_P = PRUEFZAHLEN.neu
D03, D04 = datei("03"), datei("04")
KOSMOS = ("H0_KMS", "OMEGA_B_H2", "OMEGA_C_H2")
QUADRAT, SECHSECK = 1j, cmath.exp(2j * math.pi / 3)     # τ = i und τ = e^{2πi/3} stehen im Satz selbst
FORMEN = ("sechseckig", "raumzentriert", "flächenzentriert")


def _e(k, einheit_si=1.0):
    return E[k].wert * einheit_si


# ------------------------------------------------------------------ 03 Kaleidoskop: Drehungen der Zelle
_A("03_drehungen_kompakt", (D03, "In einem kompakten flachen Raum sind nur Drehungen um"), "Zahl")
for _i, _w in enumerate(("60", "90", "120", "180")):
    _P(f"03_drehung_moeglich_{_w}", "03", "03_drehungen_kompakt", _w,
       lambda p, i=_i: sorted(set(K.drehwinkel_zelle(QUADRAT)) | set(K.drehwinkel_zelle(SECHSECK)))[i], "°",
       "sind nur Drehungen um")
_A("03_drehungen_tori", (D03, "Das sind genau die Drehsymmetrien flacher Tori"), "Zahl")
_P("03_drehung_jeder", "03", "03_drehungen_tori", "180",
   lambda p: min(set(K.drehwinkel_zelle(QUADRAT)) & set(K.drehwinkel_zelle(SECHSECK))), "°", "hat jeder")
_P("03_drehung_quadrat", "03", "03_drehungen_tori", "90", lambda p: min(K.drehwinkel_zelle(QUADRAT)), "°",
   "nur der quadratische")
_P("03_drehung_sechseck_60", "03", "03_drehungen_tori", "60", lambda p: min(K.drehwinkel_zelle(SECHSECK)), "°",
   "nur der sechseckige")
_P("03_drehung_sechseck_120", "03", "03_drehungen_tori", "120", lambda p: sorted(K.drehwinkel_zelle(SECHSECK))[1], "°",
   "nur der sechseckige")

# ------------------------------------------------------------------ 03 Kaleidoskop: Quadrupol an einem Himmel
_A("03_quadrupol_streuung", (D03, "weniger als seine kosmische Streuung (63 %)"), "Zahl")
_P("03_quadrupol_kosmische_streuung", "03", "03_quadrupol_streuung", "63", lambda p: 100 * K.streuung_kosmisch(2), "%",
   "kosmische Streuung")

_A("03_ausrichtung_isotrop", (D03, "(unendlicher Raum 40 %)"), "Zahl")
_P("03_quadrupol_isotrop_achsen", "03", "03_ausrichtung_isotrop", "40", lambda p: 100 * K.anteil_isotrop(2, 2), "%",
   "unendlicher Raum")          # Achsenformen des Würfels: 2 der 5 Formen (E_g)
_P("03_quadrupol_isotrop_m0", "03", "03_ausrichtung_isotrop", "20", lambda p: 100 * K.anteil_isotrop(1, 2), "%",
   "um ihre Achse")             # Form um die Achse (m = 0): 1 der 5 Formen

# ------------------------------------------------------------------ 03 Zelle: Grenze von a₀ und Formfaktor je Form
_A("03_zellgrenze_formen", (D03, "sechseckig bis 1,30, raumzentriert bis 1,46, flächenzentriert als dichteste"), "Zahl")
for _f, _z in zip(FORMEN, ("1,30", "1,46", "1,59")):
    _P(f"03_a0_max_{_f}", "03", "03_zellgrenze_formen", _z,
       lambda p, f=_f: K.zellgrenze(K.formfaktor_zelle(f)) / (1 + K.Q_KOSM) / 1e-10, "10⁻¹⁰ m/s²",
       "raumzentriert bis", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)
_A("03_formfaktor_formen", (D03, "(Würfel f³ = 1, sechseckig 1,15, raumzentriert 1,30, flächenzentriert 1,41)"), "Zahl")
for _f, _z in zip(FORMEN, ("1,15", "1,30", "1,41")):
    _P(f"03_formfaktor_{_f}", "03", "03_formfaktor_formen", _z, lambda p, f=_f: K.formfaktor_zelle(f), "1",
       "Würfel f³ = 1")

# ------------------------------------------------------------------ 03 Gebunden: Übergabe bei der ersten Faltung
_A("03_uebergabe_faltung", (D03, "I_F = ∫ω dt ≈ 1,1–1,25 Radiant nach der Umkehr"), "Zahl")
_P("03_I_F_keim_klein", "03", "03_uebergabe_faltung", "1,1",
   lambda p: K.uebergabe_faltung(_e("IF_E_KEIM"), _e("IF_EJ_KLEIN"), n_t=1000, n_l=400, grenze=3e-3)[0], "rad",
   "I_F = ∫ω dt",       # gröberes Gitter: die Zahl hat nur zwei Stellen, der Fixpunkt läuft hier am langsamsten
   haengt_an=("IF_E_KEIM", "IF_EJ_KLEIN"))
_P("03_I_F_keim_gross", "03", "03_uebergabe_faltung", "1,25",
   lambda p: K.uebergabe_faltung(_e("IF_E_KEIM"), _e("IF_EJ_GROSS"))[0], "rad", "I_F = ∫ω dt",
   haengt_an=("IF_E_KEIM", "IF_EJ_GROSS"))
_P("03_I_F_flach", "03", "03_uebergabe_faltung", "1,75",
   lambda p: K.uebergabe_faltung(_e("IF_E_FLACH"), _e("IF_EJ_GROSS"))[0], "rad", "flachere Anfangsprofile bis",
   haengt_an=("IF_E_FLACH", "IF_EJ_GROSS"))

# ------------------------------------------------------------------ 03 Naht: Schale außerhalb R_ta, Übergabe vor der Umkehr
_A("03_naht_schale", (D03, "w wechselt in einer Schale bis 1,16·R_ta"), "Zahl")
_P("03_schale_rand", "03", "03_naht_schale", "1,16", lambda p: K.rand_schale_aussen(_e("W_OBEN"), p.a0), "R_ta",
   "w wechselt in einer Schale bis", haengt_an=("W_OBEN",) + KOSMOS)
_A("03_uebergabe_vor_umkehr", (D03, "(∫(1 − w)ω dt ≈ 0,1)"), "Zahl")
_P("03_I_vor_umkehr", "03", "03_uebergabe_vor_umkehr", "0,1", lambda p: K.uebergabe_vor_umkehr(p.a0), "rad",
   "Die Übergabe füllt die Schale bis heute nicht", haengt_an=KOSMOS)
_A("03_schale_aussen", (D03, "Außerhalb steigt w rasch (0,95 ab 1,33·R_ta"), "Zahl")
_P("03_w_aussen_133", "03", "03_schale_aussen", "0,95", lambda p: K.gewicht_schale_aussen(_e("R_SCHALE_AUSSEN"), p.a0),
   "1", "Außerhalb steigt w rasch", haengt_an=("R_SCHALE_AUSSEN",) + KOSMOS)

# ------------------------------------------------------------------ 03 Linear: Wachstum mit der Schließungsregel
_A("03_wachstum_fsigma8", (D03, "(fσ₈ auf 1,2·10⁻³, höchstens 5·10⁻³)"), "Zahl")
_P("03_fsigma8_uebergabe", "03", "03_wachstum_fsigma8", "1,2·10⁻³",
   lambda p: -K.wachstum_uebergabe(p.z, _e("NSTERN_EINS"), p.a0), "1", "fσ₈ auf",
   haengt_an=("NSTERN_EINS", "OMEGA_R") + KOSMOS, satz="mitatmend_heute")
_P("03_fsigma8_sofort", "03", "03_wachstum_fsigma8", "5·10⁻³",
   lambda p: -K.wachstum_uebergabe(p.z, a0=p.a0, sofort=True), "1", "höchstens",
   haengt_an=("OMEGA_R",) + KOSMOS, satz="mitatmend_heute")

# ------------------------------------------------------------------ 03 Störungen bei der Rekombination
_A("03_stoerung_rekombination", (D03, "(147 Mpc: a_N ≈ 2,6–12·a₀(z))"), "Zahl")
_P("03_stoerung_aN_rekombination", "03", "03_stoerung_rekombination", "2,6",
   lambda p: K.neigung_stoerung(_e("PHI_STOERUNG"), M["BAO_SKALA"].wert * K.Mpc, p.z)
   / (p.a0 * K.a0_frueh_relativ(p.z)), "a₀(z)", "lineare Störungen über der Schwelle",
   haengt_an=("BAO_SKALA", "PHI_STOERUNG", "OMEGA_R") + KOSMOS, satz="z1090")

# ------------------------------------------------------------------ 03 Kreislauf: Horizontwachstum, Zerstrahlung
_A("03_horizont_verduennung", (D03, "heute 6,6·10¹⁷-mal schneller als die Entropieproduktion"), "Zahl")
_P("03_horizontwachstum_heute", "03", "03_horizont_verduennung", "6,6·10¹⁷",
   lambda p: K.horizontwachstum(-math.log1p(p.z), H0=p.H0) / M["ENTROPIEPRODUKTION"].wert, "1",
   "Solange die Verdünnung der Materie", haengt_an=("ENTROPIEPRODUKTION", "OMEGA_R") + KOSMOS, satz="mitatmend_heute")
_A("03_zerstrahlung", (D03, "Zerstrahlung senkt Λ, weil Strahlung ihre Faltung durch Rotverschiebung verliert"), "Zahl")
_P("03_ws_zerstrahlung", "03", "03_zerstrahlung", "−1,0002",
   lambda p: -1 + K.zustand_vakuum(K.gang_schwelle_strahlung(p.z)), "1", "heute w_s =",
   haengt_an=("OMEGA_R",) + KOSMOS, satz="mitatmend_heute")

# ------------------------------------------------------------------ 03 Linse jenseits der Naht, Rest im Haufen
_A("03_linse_jenseits_naht", (D03, "Milchstraße bei 1 Mpc ΔΣ = 0,12 M☉/pc²"), "Zahl")
_P("03_linse_milchstrasse_1mpc", "03", "03_linse_jenseits_naht", "0,12",
   lambda p: K.linse_punktmasse((1 + K.Q_KOSM) * _e("M_MILCHSTRASSE", K.M_sun), _e("R_LINSE_1MPC", K.Mpc))
   / (K.M_sun / K.pc**2), "M☉/pc²", "Milchstraße bei 1 Mpc",
   haengt_an=("M_MILCHSTRASSE", "R_LINSE_1MPC", "OMEGA_B_H2", "OMEGA_C_H2"))

_A("03_haufen_rest", (D03, "Mit der Regel allein bleibt in Haufen ein Rest"), "Zahl")
_HAUFEN = ("HAUFEN_GAS_M", "HAUFEN_GAS_RC", "HAUFEN_GAS_RNORM", "HAUFEN_GAS_RMAX", "HAUFEN_BCG_M", "HAUFEN_BCG_B",
           "HAUFEN_GAL_M", "HAUFEN_GAL_B", "NFW_M200", "NFW_C", "NFW_Z", "H0_KMS")


def _haufen_rest(p, r):
    kpc, Ms = K.kpc, K.M_sun
    mb = K.haufen_baryonen(r, _e("HAUFEN_GAS_M", Ms), _e("HAUFEN_GAS_RC", kpc), _e("HAUFEN_GAS_RNORM", K.Mpc),
                           _e("HAUFEN_GAS_RMAX", K.Mpc), _e("HAUFEN_BCG_M", Ms), _e("HAUFEN_BCG_B", kpc),
                           _e("HAUFEN_GAL_M", Ms), _e("HAUFEN_GAL_B", kpc))
    ml = K.masse_linse_nfw(r, M["NFW_M200"].wert * Ms, M["NFW_C"].wert, M["NFW_Z"].wert)
    return K.haufen_rest(r, mb, ml, p.a0)


_P("03_haufen_rest_100kpc", "03", "03_haufen_rest", "3,1", lambda p: _haufen_rest(p, _e("R_100KPC", K.kpc)), "1",
   "Faktor 3,1 bei 100 kpc", haengt_an=_HAUFEN + ("R_100KPC",))
_P("03_haufen_rest_1mpc", "03", "03_haufen_rest", "1,6", lambda p: _haufen_rest(p, _e("R_HAUFEN_REST_1MPC", K.Mpc)), "1",
   "1,6 bei 1 Mpc", haengt_an=_HAUFEN + ("R_HAUFEN_REST_1MPC",))

# ------------------------------------------------------------------ 04 Einsatz der Rückgabe (Weltzeit), Kreislauf
_A("04_rueckgabe_einsatz", (D04, "für Falten von 10 bis 10⁹ M☉ nach rund 4·10¹¹ bis 10¹² Jahren"), "Zahl")


def _rueckgabe(p, m_kennung, Ox):
    """Weltzeit, zu der T_CMB/a die Temperatur der Falte erreicht: a = M/M_T(T_CMB), weil T_H ∝ 1/M."""
    return K.weltzeit(math.log(_e(m_kennung, K.M_sun) / K.grenzmasse_hawking()), Ox=Ox, H0=p.H0) / K.yr


_P("04_rueckgabe_klein", "04", "04_rueckgabe_einsatz", "4·10¹¹",
   lambda p: _rueckgabe(p, "RUECKGABE_M_KLEIN", K.omega_de()), "Jahr", "nach rund",
   haengt_an=("RUECKGABE_M_KLEIN", "T_CMB", "OMEGA_R") + KOSMOS, satz="mitatmend_heute")
_P("04_rueckgabe_gross", "04", "04_rueckgabe_einsatz", "10¹²",
   lambda p: _rueckgabe(p, "RUECKGABE_M_GROSS", K.omega_vakuum(p.a0, p.H0)), "Jahr", "nach rund",
   haengt_an=("RUECKGABE_M_GROSS", "T_CMB", "OMEGA_R") + KOSMOS, satz="mitatmend_heute", toleranz=0.5,
   toleranz_grund="bloße Zehnerpotenz mit „rund“: Größenordnung (0,5 … 1,5·10¹²)")

_A("04_kreislauf_verletzt", (D04, "Mit konstanter Dunkler Energie ist diese Bedingung ab etwa 2,6·10¹¹ Jahren verletzt"),
   "Zahl")
_P("04_kreislauf_versiegen", "04", "04_kreislauf_verletzt", "2,6·10¹¹",
   lambda p: K.weltzeit(K.kreislauf_versiegen(M["ENTROPIEPRODUKTION"].wert, Ox=K.omega_de(), H0=p.H0),
                        Ox=K.omega_de(), H0=p.H0) / K.yr, "Jahr", "Mit konstanter Dunkler Energie",
   haengt_an=("ENTROPIEPRODUKTION", "OMEGA_R") + KOSMOS, satz="mitatmend_heute")
