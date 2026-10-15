# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 03, zweiter Teil: Expansion und Rotverschiebung, optisches Kaleidoskop, Dunkle Energie, prüfbare
Vorhersagen, Offen (Aufbau wie zahlen_01.py). Rechnung frisch aus dem Kern; Messwerte und Eingaben nur über ihre
Register. Satz überall "standard" (a₀ = 1,1·10⁻¹⁰, H₀ = 67,4, Υ_d = 0,5, z = 0); einen mitatmenden Satz für Zelle,
Vakuumteil und Kreislauf gibt es noch nicht.
"""
import math

import kern as K
from kern.eingabe import EINGABEN as E
from kern.messwert import MESSWERTE as M
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN
from .sonnenort import SONNENORT, gezeit

B = "03"
D = datei(B)
_A = AUSSAGEN.neu
_P = PRUEFZAHLEN.neu
KOSMOS = ("H0_KMS", "OMEGA_B_H2", "OMEGA_C_H2")
Q_MESS = ("OMEGA_B_H2", "OMEGA_C_H2")
M_MW = E["M_MILCHSTRASSE"].wert * K.M_sun
MIO_J = 1e6 * K.yr

# ------------------------------------------------------------------ Kaleidoskop: Zelle ohne Kreispaare
_A("03_zelle_ohne_kreispaare", (D, "L ≥ 27,4 Gpc"), "Zahl")
_P("03_zelle_mindestens", B, "03_zelle_ohne_kreispaare", "27,4", lambda p: K.KREISPAAR_GRENZE * K.d_LSS() / K.Gpc,
   "Gpc", "kürzeste geschlossene Schleife", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)
_A("03_zelle_fenster", (D, "Prüfbar ist sie nur im Fenster 27–42 Gpc"), "Zahl")
_P("03_zelle_fenster_unten", B, "03_zelle_fenster", "27", lambda p: K.KREISPAAR_GRENZE * K.d_LSS() / K.Gpc, "Gpc",
   "Prüfbar ist sie nur im Fenster", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)

# ------------------------------------------------------------------ Rückstellung: Bilanz der Faltung
_A("03_rueckstellbilanz", (D, "ΔL_Raum/R = ΣΔs/R = 0,083"), "Zahl")
_P("03_rueckstell_baryonen", B, "03_rueckstellbilanz", "0,083", lambda p: K.rueckstellbilanz(p.z), "1",
   "ΔL_Raum/R = ΣΔs/R", haengt_an=KOSMOS)
_P("03_rueckstell_ruecklaeufe", B, "03_rueckstellbilanz", "0,53",
   lambda p: K.rueckstellbilanz(p.z) * (1 + K.Q_KOSM), "1", "ν_kosm·ΣΔs/R", haengt_an=KOSMOS)

# ------------------------------------------------------------------ Vakuumteil und Θ
_A("03_vakuumteil_heute", (D, "Heute (Θ ≈ 1) ist Ω_s = π²/45 = 0,22"), "Zahl")
_P("03_omega_s_theta1", B, "03_vakuumteil_heute", "0,22", lambda p: K.omega_vakuum(p.a0, p.H0), "1",
   "Heute (Θ ≈ 1)", satz="horizont")
_P("03_omega_s_standard", B, "03_vakuumteil_heute", "0,24", lambda p: K.omega_vakuum(p.a0, p.H0), "1",
   "Heute (Θ ≈ 1)")
_A("03_theta_steuert", (D, "Θ steuert alles Kosmische der Schwelle"), "Zahl")
_P("03_omega_s_zu_rest", B, "03_theta_steuert", "16,7",
   lambda p: K.omega_vakuum(p.a0, p.H0) / (1 - K.gewicht_hintergrund(K.theta_(p.a0, K.H0_si(p.H0)))), "1",
   "heute Ω_s = 16,7·(1 − w)")
_P("03_theta_heute", B, "03_theta_steuert", "1,06", lambda p: K.theta_(p.a0, K.H0_si(p.H0)), "1",
   "Heute ist Θ = 1,06")

# ------------------------------------------------------------------ Λ heute fest, früh
_A("03_lambda_fest", (D, "Λ ist heute fest."), "Zahl")
_P("03_a0_fest_seit_z10", B, "03_lambda_fest", "0,3", lambda p: 100 * (1 - K.a0_frueh_relativ(10.0)), "%",
   "Seit z ≈ 10 ist a₀ bis auf", haengt_an=("OMEGA_R",) + KOSMOS)
_P("03_a0_rekombination", B, "03_lambda_fest", "0,76",
   lambda p: K.a0_frueh_relativ(p.z), "1", "bei der Rekombination 0,76·a₀",
   haengt_an=("Z_REKOMBINATION", "OMEGA_R") + KOSMOS, satz="z1090")
_P("03_omega_s_rekombination", B, "03_lambda_fest", "2,6·10⁻¹⁰",
   lambda p: K.omega_vakuum(p.a0 * K.a0_frueh_relativ(p.z),
                            p.H0 * K.E_hintergrund(p.z)), "1", "Ω_s(1090)",
   haengt_an=("Z_REKOMBINATION", "OMEGA_R") + KOSMOS, satz="z1090")

# ------------------------------------------------------------------ Schwelle und Zelle
_A("03_schwelle_zelle", (D, "Der Wert bindet die Schwelle an die Zelle."), "Zahl")
_P("03_zelle_a0_unten", B, "03_schwelle_zelle", "28,1", lambda p: K.zelle_aus_a0(M["A0_SPARC"].wert[0]) / K.Gpc,
   "Gpc", "verlangt L = 28,1–26,9 Gpc", haengt_an=("A0_SPARC",) + KOSMOS)
_P("03_zelle_a0_oben", B, "03_schwelle_zelle", "26,9", lambda p: K.zelle_aus_a0(M["A0_SPARC"].wert[1]) / K.Gpc,
   "Gpc", "verlangt L = 28,1–26,9 Gpc", haengt_an=("A0_SPARC",) + KOSMOS)
_P("03_zelle_a0_standard", B, "03_schwelle_zelle", "27,6", lambda p: K.zelle_aus_a0(p.a0) / K.Gpc, "Gpc",
   "1,1·10⁻¹⁰: 27,6 Gpc", haengt_an=KOSMOS)
_P("03_d_lss_heute", B, "03_schwelle_zelle", "27,8", lambda p: K.d_LSS() / K.Gpc, "Gpc",
   "Durchmesser der letzten Streufläche (27,8 Gpc", haengt_an=KOSMOS)
_P("03_zelle_beobachter_z1", B, "03_schwelle_zelle", "1,32", lambda p: K.d_LSS() / K.d_LSS_beobachter(1.0), "1",
   "L/d_LSS = 1,32", haengt_an=KOSMOS)
_P("03_a0_max_wuerfel", B, "03_schwelle_zelle", "1,13·10⁻¹⁰",
   lambda p: K.zellgrenze(1.0) / (1 + K.Q_KOSM), "m/s²", "damit ist a₀ ≤ 1,13·10⁻¹⁰",
   haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)

# ------------------------------------------------------------------ Kehrseite der Bindung
_A("03_fluchtmass_milchstrasse", (D, "Milchstraße ln 284 = 5,65"), "Zahl")
_P("03_fluchtmass_verhaeltnis", B, "03_fluchtmass_milchstrasse", "284",
   lambda p: math.exp(K.fluchtmass_bereich(M_MW, p.a0)), "1", "Milchstraße ln 284", haengt_an=("M_MILCHSTRASSE",
                                                                                                "GM_sun") + KOSMOS)
_P("03_fluchtmass_ln", B, "03_fluchtmass_milchstrasse", "5,65", lambda p: K.fluchtmass_bereich(M_MW, p.a0), "1",
   "Milchstraße ln 284 = 5,65", haengt_an=("M_MILCHSTRASSE", "GM_sun") + KOSMOS)

# ------------------------------------------------------------------ Kreislauf
_A("03_kreislauf_entropie", (D, "heute S_H/S_Λ = 0,24"), "Zahl")
_P("03_entropie_verhaeltnis", B, "03_kreislauf_entropie", "0,24",
   lambda p: K.horizontentropie(K.H0_si(p.H0)) / K.entropie_vakuum(K.schliesslaenge(p.a0)), "1",
   "heute S_H/S_Λ = 0,24")
_A("03_kreislauf_mindestrate", (D, "Mit der heutigen Entropieproduktion"), "Zahl")
_P("03_schliesslaenge_standard", B, "03_kreislauf_mindestrate", "13,2", lambda p: K.schliesslaenge(p.a0) / K.Gpc,
   "Gpc", "Λ = c²/(2a₀) = 13,2 Gpc")


def _rate(p):
    return K.kreislauf_mindestrate(K.schliesslaenge(p.a0), M["ENTROPIEPRODUKTION"].wert)[0]


_P("03_mindestrate", B, "03_kreislauf_mindestrate", "3,8·10⁻³⁷", _rate, "s⁻¹", "Mit der heutigen Entropieproduktion",
   haengt_an=("ENTROPIEPRODUKTION",))
_P("03_mindestrate_jahr", B, "03_kreislauf_mindestrate", "10⁻²⁹", lambda p: _rate(p) * K.yr, "1/Jahr",
   "rund 10⁻²⁹ je Jahr", haengt_an=("ENTROPIEPRODUKTION",), toleranz=0.5,
   toleranz_grund="im Text als Größenordnung (rund)")
_P("03_lambda_zuwachs", B, "03_kreislauf_mindestrate", "1,6·10⁻¹⁰", lambda p: _rate(p) * K.schliesslaenge(p.a0),
   "m/s", "Λ wächst um", haengt_an=("ENTROPIEPRODUKTION",))
# 1 + w_s = −⅔·d ln a₀/d ln a mit d ln a₀/d ln a = −(Λ̇/Λ)/H (zustand_vakuum gibt 1 + w_s direkt).
_P("03_ws_heute", B, "03_kreislauf_mindestrate", "1·10⁻¹⁹",
   lambda p: K.zustand_vakuum(-_rate(p) / K.H0_si(p.H0)), "1", "w_s = −1 + 1·10⁻¹⁹",
   haengt_an=("ENTROPIEPRODUKTION",))

# ------------------------------------------------------------------ Vorhersagen
_A("03_vorhersage_binaerpulsare", (D, "Beschleunigung am Sonnenort aus Binärpulsaren"), "Vorhersage")
# Milchstraße als Punktmasse
_P("03_pulsar_regel", B, "03_vorhersage_binaerpulsare", "6,3·10⁻¹⁹",
   lambda p: K.pulsar_rate(K.a_regel(K.neigung(M_MW, E["R_SONNE_PUNKTMASSE"].wert * K.kpc), p.a0)), "s⁻¹", "radial Ṗ/P",
   haengt_an=("M_MILCHSTRASSE", "R_SONNE_PUNKTMASSE", "GM_sun"))
_P("03_pulsar_newton", B, "03_vorhersage_binaerpulsare", "4,1·10⁻¹⁹",
   lambda p: K.pulsar_rate(K.neigung(M_MW, E["R_SONNE_PUNKTMASSE"].wert * K.kpc)), "s⁻¹", "radial Ṗ/P",
   haengt_an=("M_MILCHSTRASSE", "R_SONNE_PUNKTMASSE", "GM_sun"))
_A("03_vorhersage_massenaenderung", (D, "Nach einer Massenänderung folgt die Außenbeschleunigung"), "Vorhersage")
_P("03_rueckkopplung_zeit", B, "03_vorhersage_massenaenderung", "1",
   lambda p: K.nu(K.neigung(M_MW, E["R_RUECKKOPPLUNG"].wert * K.kpc), p.a0) * E["R_RUECKKOPPLUNG"].wert * K.kpc
   / K.c / MIO_J, "Mio. Jahre", "ν·r/c ≈ 1 Mio. Jahre", haengt_an=("M_MILCHSTRASSE", "R_RUECKKOPPLUNG", "GM_sun"))
_A("03_vorhersage_haufenstoss", (D, "Linsenmaximum und Gas trennen sich bei Haufenstößen"), "Vorhersage")
# unten kleinstes I_F bei größtem Δ, oben größtes I_F bei kleinstem Δ; uebergabezeit gibt 1/ω_orb
_P("03_haufenstoss_unten", B, "03_vorhersage_haufenstoss", "0,9",
   lambda p: K.haufenstoss_trennzeit(E["I_F_HAUFENSTOSS"].wert[0],
                                     1 / K.uebergabezeit(E["DELTA_500"].wert, p.z)[0]) / K.Gyr,
   "Mrd. Jahre", "Mrd. Jahre bei Δ = 200–500", haengt_an=("I_F_HAUFENSTOSS", "DELTA_200", "DELTA_500") + KOSMOS)
_P("03_haufenstoss_oben", B, "03_vorhersage_haufenstoss", "1,8",
   lambda p: K.haufenstoss_trennzeit(E["I_F_HAUFENSTOSS"].wert[1],
                                     1 / K.uebergabezeit(E["DELTA_200"].wert, p.z)[0]) / K.Gyr,
   "Mrd. Jahre", "Mrd. Jahre bei Δ = 200–500", haengt_an=("I_F_HAUFENSTOSS", "DELTA_200", "DELTA_500") + KOSMOS)

_A("03_vorhersage_flaechendichte", (D, "Zentrale Flächendichte aus der Rotationskurve (Toomre-Integral)"),
   "Vorhersage")
_SIG = K.M_sun / K.pc**2


def _lelli(S):
    """Gemessener Verlauf Σ_dyn(Σ₀) (Messwert LELLI_SIGMA_DYN), Σ₀ in M☉/pc²."""
    sm, sc, al, be = M["LELLI_SIGMA_DYN"].wert
    return sm * (1 + S / sc) ** al * (S / sc) ** be


def _exp_gegen_lelli(S, a0):
    return math.log10(K.sigma_dyn_expscheibe(S * _SIG, a0) / _SIG / _lelli(S))


_P("03_expscheibe_hoch", B, "03_vorhersage_flaechendichte", "0,08",
   lambda p: _exp_gegen_lelli(E["SIGMA0_EXPSCHEIBE"].wert[1], p.a0), "dex", "eine ideale Exponentialscheibe läge",
   haengt_an=("LELLI_SIGMA_DYN", "SIGMA0_EXPSCHEIBE"))
_P("03_expscheibe_tief", B, "03_vorhersage_flaechendichte", "0,22",
   lambda p: _exp_gegen_lelli(E["SIGMA0_EXPSCHEIBE"].wert[0], p.a0), "dex", "eine ideale Exponentialscheibe läge",
   haengt_an=("LELLI_SIGMA_DYN", "SIGMA0_EXPSCHEIBE"))
_P("03_sparc_auswahl", B, "03_vorhersage_flaechendichte", "153",
   lambda p: len(K.sparc_tafel(K.SPARC_SATZ["q_max"], K.SPARC_SATZ["inc_min"])), "Galaxien", "SPARC, 153 Galaxien",
   haengt_an=("SPARC_Q_MAX", "SPARC_INC_MIN"))

# ------------------------------------------------------------------ Offen: Herkunft von Q und Zelle
_A("03_offen_herkunft_q", (D, "Herkunft von Q: Q ist dynamisch, nicht topologisch"), "Zahl")
_P("03_q_abstand_fuenf", B, "03_offen_herkunft_q", "5,6", lambda p: (K.Q_KOSM - 5) / K.Q_fehler(), "σ",
   "liegt 5,6σ daneben", haengt_an=Q_MESS)
_P("03_q_max_wuerfel", B, "03_offen_herkunft_q", "5,5", lambda p: K.q_fenster(p.a0)[1], "1",
   "Q ≤ 5,5 (Würfel)", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)
_P("03_rstern_zu_lambda", B, "03_offen_herkunft_q", "0,97",
   lambda p: K.r_stern(K.rho_b_mitbewegt() * K.zelle_aus_a0(p.a0) ** 3, p.a0) / K.schliesslaenge(p.a0), "1",
   "r*/Λ = √(6/ν_kosm) = 0,97", haengt_an=KOSMOS)
_P("03_zellgrenze_wuerfel", B, "03_offen_herkunft_q", "7,17·10⁻¹⁰", lambda p: K.zellgrenze(1.0), "m/s²",
   "(1 + Q)·a₀ ≤ 7,17·10⁻¹⁰", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)
_P("03_q_fehler_planck", B, "03_offen_herkunft_q", "1,2", lambda p: 100 * K.Q_fehler() / K.Q_KOSM, "%",
   "(heute 1,0–1,2 %)", haengt_an=Q_MESS)

# ------------------------------------------------------------------ Offen: Warp
_A("03_offen_warp", (D, "Der Warp der Milchstraße präzediert gemessen rechtläufig"), "Zahl")
_P("03_warp_viertel", B, "03_offen_warp", "0,14",
   lambda p: K.warp_viertelzeit(M["WARP_PRAEZESSION"].wert * 1e3 / K.kpc) / K.Gyr, "Mrd. Jahre",
   "nach 0,14 Mrd. Jahren", haengt_an=("WARP_PRAEZESSION",))

# ------------------------------------------------------------------ Offen: SPARC-Außenpunkte
_A("03_offen_aussenpunkte", (D, "prüfen die äußersten SPARC-Punkte"), "Zahl")


def _aussen_aN(p):
    """a_N/a₀ am äußersten Punkt jeder Kurve (Tafel mit SPARC_SATZ, Υ_d aus dem Parametersatz)."""
    t = K.sparc_tafel(K.SPARC_SATZ["q_max"], K.SPARC_SATZ["inc_min"])
    return [float(K.v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], p.Upsilon_d)[-1] / g["R"][-1] / p.a0) for g in t]


_P("03_aussen_zahl", B, "03_offen_aussenpunkte", "74",
   lambda p: sum(x <= E["SPARC_RAND_AUSWAHL"].wert for x in _aussen_aN(p)), "Galaxien", "74 Galaxien bis a_N",
   haengt_an=("SPARC_RAND_AUSWAHL", "SPARC_Q_MAX", "SPARC_INC_MIN", "SPARC_UPS_B"))
_P("03_aussen_kleinstes", B, "03_offen_aussenpunkte", "0,008", lambda p: min(_aussen_aN(p)), "a₀",
   "74 Galaxien bis a_N = 0,008·a₀", haengt_an=("SPARC_Q_MAX", "SPARC_INC_MIN", "SPARC_UPS_B"))

# ------------------------------------------------------------------ Offen: druckgestützte Systeme im äußeren Feld
_A("03_offen_druckgestuetzt", (D, "Druckgestützte Systeme im äußeren Feld"), "Zahl")
_WIRT = {"PAL14": "M_MILCHSTRASSE", "NGC2419": "M_MILCHSTRASSE", "CRATER2": "M_MILCHSTRASSE",
         "DRACO": "M_MILCHSTRASSE", "FORNAX": "M_MILCHSTRASSE", "DF2": "M_NGC1052", "DF4": "M_NGC1052"}


def _sys(k):
    """(M, R_e, g_e) eines Systems in SI; g_e = Newton-Feld des Wirts am heutigen Abstand."""
    g_e = K.neigung(E[_WIRT[k]].wert * K.M_sun, M[f"{k}_ABSTAND"].wert * K.kpc)
    return M[f"{k}_M"].wert * K.M_sun, M[f"{k}_RE"].wert * K.kpc, g_e


def _sigma(k, p, regel=True, aussen=True):
    m, re, g_e = _sys(k)
    return K.sigma_plummer_jeans(m, re, g_e if aussen else 0.0, p.a0, regel=regel) / 1e3


def _dep(k):
    return (f"{k}_M", f"{k}_RE", f"{k}_ABSTAND", _WIRT[k], "GM_sun")


for _k, _t in (("PAL14", "1,8"), ("NGC2419", "6,5"), ("FORNAX", "14,0"), ("DF2", "17,7"), ("DF4", "19,2")):
    _P(f"03_druck_regel_{_k.lower()}", B, "03_offen_druckgestuetzt", _t, lambda p, k=_k: _sigma(k, p), "km/s",
       "für Pal 14 1,8", haengt_an=_dep(_k))
for _k, _t in (("CRATER2", "1,4"), ("DRACO", "3,3")):
    _P(f"03_druck_regel_{_k.lower()}", B, "03_offen_druckgestuetzt", _t, lambda p, k=_k: _sigma(k, p), "km/s",
       "Für Crater II und Draco gibt sie", haengt_an=_dep(_k))
for _k, _t in (("DF2", "6,3"), ("DF4", "7,2")):
    _P(f"03_druck_verstaerkung_{_k.lower()}", B, "03_offen_druckgestuetzt", _t,
       lambda p, k=_k: K.verstaerkung_druck(_sigma(k, p), _sigma(k, p, regel=False)), "1", "statt 6,3 bzw. 7,2",
       haengt_an=_dep(_k))
_P("03_druck_fornax_verlangt", B, "03_offen_druckgestuetzt", "4,5",
   lambda p: K.verstaerkung_druck(M["FORNAX_SIGMA"].wert, _sigma("FORNAX", p, regel=False)), "1",
   "die dreifache Verstärkung (4,5)",
   haengt_an=_dep("FORNAX") + ("FORNAX_SIGMA",))
for _k, _t in (("FORNAX", "5,5"), ("DRACO", "1,2")):
    _P(f"03_druck_newton_geglaettet_{_k.lower()}", B, "03_offen_druckgestuetzt", _t,
       lambda p, k=_k: _sigma(k, p, regel=False), "km/s", "aber auch Fornax (5,5) und Draco (1,2 km/s)",
       haengt_an=_dep(_k))
for _k, _t in (("PAL14", "0,665"), ("NGC2419", "5,09"), ("DF4", "7,15")):
    _P(f"03_druck_newton_{_k.lower()}", B, "03_offen_druckgestuetzt", _t, lambda p, k=_k: _sigma(k, p, regel=False),
       "km/s", "unter dem Newton-Wert (0,665; 5,09; 7,15 km/s)", haengt_an=_dep(_k))
_P("03_druck_draco_isoliert", B, "03_offen_druckgestuetzt", "4,74", lambda p: _sigma("DRACO", p, aussen=False), "km/s",
   "Draco über dem Wert ohne äußeres Feld (4,74 km/s)", haengt_an=_dep("DRACO"))

# ------------------------------------------------------------------ Vorhersage: Quadrupol im Sonnensystem
_A("03_quadrupol_sonnensystem", (D, "Kein anomaler Quadrupol im Sonnensystem über"), "Vorhersage")
_P("03_gezeit_quadrupol", B, "03_quadrupol_sonnensystem", "5,6·10⁻³⁰", gezeit, "s⁻²",
   "Kein anomaler Quadrupol im Sonnensystem über", haengt_an=SONNENORT, satz="Sonnenkreis_gemessen")
