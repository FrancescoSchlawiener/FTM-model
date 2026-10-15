# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 04 (Zeit, Verschränkung, Syntropie, Gedächtnis). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "04_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
"""
import math

import kern as K
from kern.eingabe import EINGABEN as E
from kern.messwert import MESSWERTE as M
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "04"
D = datei(B)

# Beispiel (ein echter Eintrag; weitere darunter im selben Aufbau)
AUSSAGEN.neu("04_grenze_rueckgabe", (D, "heute wächst jede astrophysikalische Falte netto (Grenze 4,5·10²² kg)"), "Zahl")
PRUEFZAHLEN.neu("04_grenzmasse_T_CMB", B, "04_grenze_rueckgabe", "4,5·10²²", lambda p: K.grenzmasse_hawking(K.T_CMB), "kg",
                "astrophysikalische Falte netto", haengt_an=("T_CMB", "G"))
_TAG = 86400.0   # Sekunden je Tag (Einheit)

AUSSAGEN.neu("04_gangrate", (D, "Die Gangrate hängt nur von r/Δs ab"), "Zahl")
PRUEFZAHLEN.neu("04_takt_100r0", B, "04_gangrate", "0,995",
                lambda p: K.takt(K.M_sun, E["TAKT_R_100"].wert * K.r0(K.M_sun)), "1", "Gangrate",
                haengt_an=("TAKT_R_100", "GM_sun"))
PRUEFZAHLEN.neu("04_takt_2r0", B, "04_gangrate", "0,707",
                lambda p: K.takt(K.M_sun, E["TAKT_R_2"].wert * K.r0(K.M_sun)), "1", "Gangrate",
                haengt_an=("TAKT_R_2", "GM_sun"))

AUSSAGEN.neu("04_uhren", (D, "Eine Uhr an der Erdoberfläche geht pro Jahr 0,022 s nach"), "Zahl")
PRUEFZAHLEN.neu("04_uhr_erde", B, "04_uhren", "0,022", lambda p: K.delta(K.M_earth, K.R_earth) * K.yr, "s",
                "Erdoberfläche", haengt_an=("GM_earth", "R_earth"))
PRUEFZAHLEN.neu("04_uhr_sonne", B, "04_uhren", "67", lambda p: K.delta(K.M_sun, K.R_sun) * K.yr, "s",
                "Sonnenoberfläche", haengt_an=("GM_sun", "R_sun"))
PRUEFZAHLEN.neu("04_uhr_neutronenstern", B, "04_uhren", "69,5",
                lambda p: K.delta(E["NS_M"].wert * K.M_sun, E["NS_R"].wert * 1e3) * K.yr / _TAG, "d", "Neutronenstern",
                haengt_an=("NS_M", "NS_R", "GM_sun"))

AUSSAGEN.neu("04_endrate", (D, "H_s = π²c/(√45·Λ) = H₀√Ω_s = 33 km/s/Mpc"), "Zahl")
PRUEFZAHLEN.neu("04_endrate_H_s", B, "04_endrate", "33", lambda p: K.H_vakuum(K.schliesslaenge(p.a0)) * K.Mpc / 1e3,
                "km/s/Mpc", "Endzustand")

AUSSAGEN.neu("04_faltung_horizont", (D, "Δs_H = c/6H (heute 0,74 Gpc = 0,056·Λ)"), "Zahl")
PRUEFZAHLEN.neu("04_ds_H_heute", B, "04_faltung_horizont", "0,74", lambda p: K.faltung_horizont(K.H0_si(p.H0)) / K.Gpc,
                "Gpc", "Δs_H = c/6H")
PRUEFZAHLEN.neu("04_ds_H_zu_Lambda", B, "04_faltung_horizont", "0,056",
                lambda p: K.faltung_horizont(K.H0_si(p.H0)) / K.schliesslaenge(p.a0), "1", "Δs_H = c/6H")

AUSSAGEN.neu("04_weisser_zwerg", (D, "Kohlenstoff-Sauerstoff-Zwerg von 0,6 M☉ ist mit δ = 1,0·10⁻⁴ tiefer"), "Zahl")
PRUEFZAHLEN.neu("04_delta_weisser_zwerg", B, "04_weisser_zwerg", "1,0·10⁻⁴",
                lambda p: K.delta(E["WZ_M"].wert * K.M_sun, E["WZ_R"].wert), "1", "Kohlenstoff-Sauerstoff-Zwerg",
                haengt_an=("WZ_M", "WZ_R", "GM_sun"))

# ------------------------------------------------------------------ Syntropie: Export der Erde


def _export_erde(p):
    """Erde im Strahlungsstrom der Sonne: aufgenommen mit T_eff☉, abgegeben mit der Gleichgewichtstemperatur."""
    P = K.leistung_absorbiert(K.L_sun, K.AU, K.R_earth, M["ALBEDO_ERDE"].wert)
    return K.entropieexport(P, K.temperatur_abstrahlung(P, K.R_earth), K.TEFF_sun)


AUSSAGEN.neu("04_export_erde", (D, "Die Erde exportiert 6,1·10¹⁴ W/K"), "Zahl")
PRUEFZAHLEN.neu("04_entropieexport_erde", B, "04_export_erde", "6,1·10¹⁴", _export_erde, "W/K", "Die Erde exportiert",
                haengt_an=("ALBEDO_ERDE", "L_sun", "R_earth", "TEFF_sun"))

# ------------------------------------------------------------------ Entropie, Page-Zeit, Rückverdichtung der Falte
AUSSAGEN.neu("04_entropie_falte", (D, "S = 36π·k·(Δs/ℓ_P)² = 4π·k·α_G(M)"), "Zahl")
PRUEFZAHLEN.neu("04_entropie_falte_36", B, "04_entropie_falte", "36",
                lambda p: K.entropie_falte(K.M_sun) / (math.pi * (K.ds(K.M_sun) / K.l_P)**2), "1",
                "S = 36π·k·(Δs/ℓ_P)²", haengt_an=("GM_sun",))
AUSSAGEN.neu("04_page_zeit", (D, "Die Information kehrt mit der Strahlung zurück (Page-Zeit 0,54 der Lebensdauer)"), "Zahl")
PRUEFZAHLEN.neu("04_page_zeit", B, "04_page_zeit", "0,54", lambda p: K.page_zeit(M["BETA_PAGE"].wert), "1",
                "Page-Zeit 0,54", haengt_an=("BETA_PAGE",))
AUSSAGEN.neu("04_fuellung_horizont", (D, "Die Falten füllen den Entropiehaushalt des Horizonts heute nur zu 10⁻¹⁸"), "Zahl")
PRUEFZAHLEN.neu("04_fuellung_horizont", B, "04_fuellung_horizont", "10⁻¹⁸",
                lambda p: M["S_FALTEN_HEUTE"].wert / K.horizontentropie(K.H0_si(p.H0)), "1",
                "Entropiehaushalt des Horizonts heute", haengt_an=("S_FALTEN_HEUTE",), toleranz=0.5,
                toleranz_grund="Größenordnung (≈ 10ⁿ): ±50 %")
AUSSAGEN.neu("04_rueckverdichtung", (D, "Rückverdichtung gelingt nur hierarchisch aufwärts (N ≥ 1,5·z_f² Teile)"), "Zahl")
PRUEFZAHLEN.neu("04_rueckverdichtung_faktor", B, "04_rueckverdichtung", "1,5",
                lambda p: K.rueckverdichtung_teile(E["Z_F_EINS"].wert, M["BETA_PAGE"].wert), "1", "N ≥ 1,5·z_f²",
                haengt_an=("BETA_PAGE", "Z_F_EINS"))

AUSSAGEN.neu("04_chsh_oberschwingung", (D, "jede Oberschwingung senkt S (zwei Harmonische des Dreiecks: 2,36)"), "Zahl")
PRUEFZAHLEN.neu("04_chsh_zwei_harmonische", B, "04_chsh_oberschwingung", "2,36",
                lambda p: K.chsh_oberschwingungen(E["CHSH_HARMONISCHE"].wert), "1", "zwei Harmonische des Dreiecks",
                haengt_an=("CHSH_HARMONISCHE",))


def _rueckgabe_nenner(p):
    """Nenner N in dM/dt = −ħc⁴/(Nπ·G²·M²) aus der Lebensdauer: dM/dt = −M/(3t_F) (t_F ∝ M³); unabhängig von M."""
    M = K.M_sun
    return K.hbar * K.c**4 / (math.pi * K.G**2 * M**2) / (M / (3 * K.lebensdauer_falte(M)))


AUSSAGEN.neu("04_rueckgabe", (D, "Rückgabe dM/dt = −ħc⁴/(15360π·G²·M²)"), "Zahl")
PRUEFZAHLEN.neu("04_rueckgabe_nenner", B, "04_rueckgabe", "15360", _rueckgabe_nenner, "1", "Rückgabe dM/dt")
