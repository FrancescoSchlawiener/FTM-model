# -*- coding: utf-8 -*-
"""
Prüfzahlen der Gruppe Galaxie (Auftrag F), Bausteine 03 und 07: Zahlen, deren Kernformel bisher fehlte (Feldform in
Scheiben, Sonnenort, Ringpräzession, Warp, SPARC-Lesungen, druckgestützte Systeme). Aufbau wie zahlen_03a.py;
Rechnung frisch aus dem Kern, Messwerte und Eingaben nur über ihre Register (in haengt_an genannt). Große Lesungen
(profilierte Tafeln, Löser) hält der Kern je Lauf im Speicher; die Prüfzahl ruft die Kernformel trotzdem jedes Mal,
damit die Aufrufverfolgung sie sieht. Zweite Fundstellen derselben Größe stehen als Wiederholung am Ende.
"""
import math

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE
from kern.parametersatz import PARAMETERSAETZE

from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN
from .wiederholung import WIEDERHOLUNGEN

D03, D07 = datei("03"), datei("07")
_a = AUSSAGEN.neu
_p = PRUEFZAHLEN.neu
_w = WIEDERHOLUNGEN.neu
_S = K.SPARC_SATZ
_SP = ("SPARC_Q_MAX", "SPARC_INC_MIN", "SPARC_UPS_B", "SPARC_UPS_STREUUNG")
MSPC = K.M_sun / K.pc**2


def E(k):
    return EINGABEN[k].wert


def M(k):
    return MESSWERTE[k].wert


def _tafel():
    return K.sparc_tafel(_S["q_max"], _S["inc_min"])


def _satz():
    return K.sparc_profil_satz(_tafel())


# ================================================================== 07: Feldform gegen radiale Regel in der Scheibe
_a("07_feldform_scheibe", (D07, "Exponentialscheibe, Feldform gegen radiale Regel"), "Zahl")
_FF = "Exponentialscheibe, Feldform"
_p("07_feldform_viertel_hoch", "07", "07_feldform_scheibe", "0,82",
   lambda p: 1 / K.eta_scheibe("exp", E("RD_RSTERN_HOCH"), E("X_RD_VIERTEL")), "1", _FF,
   haengt_an=("RD_RSTERN_HOCH", "X_RD_VIERTEL"))
_p("07_feldform_viertel_tief", "07", "07_feldform_scheibe", "0,53",
   lambda p: 1 / K.eta_scheibe("exp", E("RD_RSTERN_TIEF"), E("X_RD_VIERTEL")), "1", _FF,
   haengt_an=("RD_RSTERN_TIEF", "X_RD_VIERTEL"))
_p("07_feldform_rd_tief", "07", "07_feldform_scheibe", "0,79",
   lambda p: 1 / K.eta_scheibe("exp", E("RD_RSTERN_TIEF"), E("X_RD_EINS")), "1", _FF,
   haengt_an=("RD_RSTERN_TIEF", "X_RD_EINS"))
_p("07_feldform_rd_hoch", "07", "07_feldform_scheibe", "0,92",
   lambda p: 1 / K.eta_scheibe("exp", E("RD_RSTERN_HOCH"), E("X_RD_EINS")), "1", _FF,
   haengt_an=("RD_RSTERN_HOCH", "X_RD_EINS"))
_p("07_feldform_gleichstand", "07", "07_feldform_scheibe", "3,4",
   lambda p: K.gleichstand_scheibe(E("RD_RSTERN_MITTEL")), "R_d", _FF, haengt_an=("RD_RSTERN_MITTEL",))

_a("07_feldform_naeherung", (D07, "die Feldform auf 0,02 dex"), "Zahl")
_p("07_feldform_naeherung_dex", "07", "07_feldform_naeherung", "0,02",
   lambda p: max(abs(math.log10(K.feldform_naeherung("exp", h, x))) for h in
                 (E("RD_RSTERN_HOCH"), E("RD_RSTERN_MITTEL"), E("RD_RSTERN_TIEF")) for x in E("X_NAEHERUNG")),
   "dex", "Für die Kuzmin-Scheibe gilt", haengt_an=("RD_RSTERN_HOCH", "RD_RSTERN_MITTEL", "RD_RSTERN_TIEF", "X_NAEHERUNG"))

# ================================================================== 07: SPARC, Gewicht des senkrechten Scheibenfelds
_a("07_gewicht_f", (D07, "geht das senkrechte Scheibenfeld mit dem Gewicht"), "Zahl")
_FS = "An den SPARC-Rotationskurven, mit"
_p("07_gewicht_f_unten", "07", "07_gewicht_f", "0,03", lambda p: K.gewicht_f(_tafel(), p.Upsilon_d)[0], "1", _FS,
   haengt_an=_SP, satz="Upsilon_0.4")
_p("07_gewicht_f_oben", "07", "07_gewicht_f", "0,13", lambda p: K.gewicht_f(_tafel(), p.Upsilon_d)[0], "1", _FS,
   haengt_an=_SP, satz="Upsilon_0.6")
_a("07_gewicht_f_ohne", (D07, "Das positive f trägt eine Galaxie"), "Zahl")
_p("07_gewicht_f_ohne_ic2574", "07", "07_gewicht_f_ohne", "0,00",
   lambda p: K.gewicht_f(_tafel(), p.Upsilon_d, ("IC2574",))[0], "1", "Das positive f", haengt_an=_SP)
_p("07_gewicht_f_ohne_fehler", "07", "07_gewicht_f_ohne", "0,03",
   lambda p: K.gewicht_f(_tafel(), p.Upsilon_d, ("IC2574",))[1], "1", "Das positive f", haengt_an=_SP)

# ================================================================== 07: Dynamikmasse am Sonnenort (Scheibenmodell)


def _saeule(p, form):
    """Σ_z [M☉/pc²] der Säule |z| < 1,1 kpc am Sonnenort im Scheibenmodell (Milchstraße als Exponentialscheibe)."""
    M_, Rd, R0 = E("M_MILCHSTRASSE") * K.M_sun, E("RD_MILCHSTRASSE_26") * K.kpc, E("R_SONNE_SCHEIBE") * K.kpc
    gR = float(K.g_expscheibe(M_, Rd, R0))
    return K.dynamikmasse_ort(M("SIGMA_BAR_SONNE") * MSPC, gR, E("Z_SAEULE") * K.kpc, R0, form, p.a0) / MSPC


_H_SAEULE = ("SIGMA_BAR_SONNE", "M_MILCHSTRASSE", "RD_MILCHSTRASSE_26", "R_SONNE_SCHEIBE", "Z_SAEULE", "GM_sun")
_a("07_dynamikmasse", (D07, "Die Dynamikmasse am Sonnenort"), "Zahl")
_DM = "Die Dynamikmasse am Sonnenort"
_p("07_dynamikmasse_kern", "07", "07_dynamikmasse", "57", lambda p: _saeule(p, "kern"), "M☉/pc²", _DM,
   haengt_an=_H_SAEULE)
_p("07_dynamikmasse_ort_schicht", "07", "07_dynamikmasse", "69", lambda p: _saeule(p, "ort_schicht"), "M☉/pc²", _DM,
   haengt_an=_H_SAEULE)
_p("07_dynamikmasse_ort_ebene", "07", "07_dynamikmasse", "70", lambda p: _saeule(p, "ort_ebene"), "M☉/pc²", _DM,
   haengt_an=_H_SAEULE)


def _milchstrasse(p):
    """Feldform-Lösung des baryonischen Massenmodells (McMillan 2017) in Einheiten G = a₀ = 1, Länge kpc."""
    teile = dict(bulge=M("MCMILLAN_BULGE"), duenn=M("MCMILLAN_DUENN"), dick=M("MCMILLAN_DICK"), HI=M("MCMILLAN_HI"),
                 H2=M("MCMILLAN_H2"))
    q = K.milchstrasse_quelle(teile)
    n, L, h0 = E("GITTER_MILCHSTRASSE")
    return q, K.feldform_loeser(q, int(n), L, h0)


_H_MCM = ("MCMILLAN_BULGE", "MCMILLAN_DUENN", "MCMILLAN_DICK", "MCMILLAN_HI", "MCMILLAN_H2", "GITTER_MILCHSTRASSE",
          "GM_sun")
_p("07_dynamikmasse_feldform", "07", "07_dynamikmasse", "92",
   lambda p: K.saeule_loesung(_milchstrasse(p)[1], E("R_SONNE_SCHEIBE"), E("Z_SAEULE")) * (p.a0 / K.G) / MSPC,
   "M☉/pc²", _DM, haengt_an=_H_MCM + ("R_SONNE_SCHEIBE", "Z_SAEULE"))

# ================================================================== 07: Präzession geneigter Ringe, Kern- gegen Feldform


def _faktor_ring(p, R):
    M_, Rd, h = E("M_MILCHSTRASSE") * K.M_sun, E("RD_MILCHSTRASSE_26") * K.kpc, E("H_SCHICHT_03") * K.kpc
    return (K.knotenpraezession_scheibe(M_, Rd, h, R * K.kpc, "ort", p.a0)
            / K.knotenpraezession_scheibe(M_, Rd, h, R * K.kpc, "kern", p.a0))


_a("07_ringe_kern_feld", (D07, "Präzessionsraten geneigter Ringe unterscheiden Kern- und Feldform"), "Zahl")
_H_RING = ("M_MILCHSTRASSE", "RD_MILCHSTRASSE_26", "H_SCHICHT_03", "R_RINGE_AUSSEN", "GM_sun")
_p("07_ringe_faktor_innen", "07", "07_ringe_kern_feld", "2,3", lambda p: _faktor_ring(p, E("R_RINGE_AUSSEN")[0]), "1",
   "Präzessionsraten geneigter", haengt_an=_H_RING)
_p("07_ringe_faktor_aussen", "07", "07_ringe_kern_feld", "3,0", lambda p: _faktor_ring(p, E("R_RINGE_AUSSEN")[1]), "1",
   "Präzessionsraten geneigter", haengt_an=_H_RING)

# ================================================================== 07: Hill-Radius statt Gezeit (Scheibenmodell, Sonnenort)


def _g_sonne_mn(p):
    """Newton-Feld des Scheibenmodells aus vier Miyamoto-Nagai-Scheiben am Sonnenort R₀ (gemessen)."""
    if p.sonnenkreis != "gemessen":
        raise ValueError("Sonnenort nur im Satz mit gemessenem Sonnenort")
    R0 = M("R0_SONNE") * K.kpc
    return R0, sum(K.g_miyamoto_nagai(m * K.M_sun, a * K.kpc, b * K.kpc, R0, 0.0, "R") for m, a, b in E("MN_SUMME_SONNE"))


_a("07_hill_sonne", (D07, "die Sonne liefe dann im selben Scheibenmodell"), "Zahl")
_p("07_hill_sonne_newton", "07", "07_hill_sonne", "188",
   lambda p: math.sqrt(_g_sonne_mn(p)[1] * _g_sonne_mn(p)[0]) / 1e3, "km/s", "Mit dem Hill-Radius",
   haengt_an=("MN_SUMME_SONNE", "R0_SONNE", "GM_sun"), satz="Sonnenkreis_gemessen")
_p("07_hill_sonne_regel", "07", "07_hill_sonne", "229",
   lambda p: math.sqrt(K.a_regel(_g_sonne_mn(p)[1], p.a0) * _g_sonne_mn(p)[0]) / 1e3, "km/s", "Mit dem Hill-Radius",
   haengt_an=("MN_SUMME_SONNE", "R0_SONNE", "GM_sun"), satz="Sonnenkreis_gemessen")

# ================================================================== 07: Restfeld im Planetensystem (Quadrupol der Glättung)


def _quadrupol(p):
    """Quadrupol des Zusatzfelds um die über ℓ☉ geglättete Sonne im geglätteten Feld der Galaxie am Sonnenort."""
    from .sonnenort import gezeit
    g_e = K.a_N_aus_kreisbahn(M("V0_SONNE") * 1e3, M("R0_SONNE") * K.kpc, p.a0)
    return K.quadrupol_glaettung(K.GM_sun, g_e, K.glaettungslaenge(K.M_sun, gezeit(p)), p.a0)


_a("07_quadrupol_planetensystem", (D07, "Im Planetensystem bleibt ein Restfeld"), "Zahl")
_p("07_quadrupol_glaettung", "07", "07_quadrupol_planetensystem", "8,6·10⁻³¹", _quadrupol, "s⁻²",
   "Im Planetensystem bleibt", haengt_an=("V0_SONNE", "R0_SONNE", "OORT_A", "OORT_B", "RHO_B_SONNE", "GM_sun"),
   satz="Sonnenkreis_gemessen")

# ================================================================== 07: Spektraldimension und Boden an SPARC
_a("07_ds_fehler", (D07, "An den SPARC-Kurven gibt die Form der Regel d_s"), "Zahl")
_p("07_ds_fehler_ziehung", "07", "07_ds_fehler", "0,07",
   lambda p: K.ds_streuung(_tafel(), p.Upsilon_d, E("SIG_GAS_8") * 1e3, a0=p.a0), "1", "Messbar ist d_s",
   haengt_an=_SP + ("SIG_GAS_8",))


def _boden(p):
    return K.boden_grenze(_tafel(), p.Upsilon_d, E("SIG_GAS_8") * 1e3, E("S_MAX_BODEN"), a0=p.a0)


_a("07_boden", (D07, "an SPARC liegt ein Boden unter"), "Zahl")
_H_BODEN = _SP + ("SIG_GAS_8", "S_MAX_BODEN")
_p("07_boden_grenze", "07", "07_boden", "0,03", _boden, "a₀", "Eine Drift zum Torus hin", haengt_an=_H_BODEN)
_p("07_boden_peclet", "07", "07_boden", "0,016", lambda p: _boden(p) * p.a0 / K.boden_drift(1.0, p.a0), "1",
   "Eine Drift zum Torus hin", haengt_an=_H_BODEN)
_p("07_boden_bandbreite", "07", "07_boden", "25", lambda p: K.boden_quer(1.0, 1.0, p.a0) / (_boden(p) * p.a0),
   "ϖ_R", "Eine Drift zum Torus hin", haengt_an=_H_BODEN)

# ================================================================== 03: Schwelle an SPARC, profiliert


def _ups_b(p, ups_d=None):
    """Υ_b: fest 0,7 (Setzung SPARC_UPS_B). Im Satz Upsilon_skaliert (Wahl) wird Υ_b mit Υ_d mitskaliert: Υ_b/Υ_d wie
    im Standardsatz (Lesart des Archivs)."""
    if not p.ups_b_skaliert:
        return _S["ups_b"]
    return _S["ups_b"] * (p.Upsilon_d if ups_d is None else ups_d) / PARAMETERSAETZE["standard"].Upsilon_d


def _a0(p, sig_int, verlust="linear", ups_d=None):
    ud = p.Upsilon_d if ups_d is None else ups_d
    return K.a0_profiliert(_satz(), ud, sig_int, verlust, ups_b=_ups_b(p, ud))


def _a0_varianten(p):
    """Lesungen der Varianten: Υ_d 0,4 / 0,5 / 0,6, σ_int 0,10 / 0,20 (linear), dazu Υ_d 0,5 robust (soft_l1); Υ_b
    fest oder mitskaliert nach dem Satz p."""
    s1, s2 = E("SIGMA_INT")
    lin = {(k, s): _a0(p, s, ups_d=PARAMETERSAETZE[k].Upsilon_d)[0]
           for k in ("Upsilon_0.4", "standard", "Upsilon_0.6") for s in (s1, s2)}
    return lin, _a0(p, s1, "soft_l1", ups_d=PARAMETERSAETZE["standard"].Upsilon_d)[0]


def _a0_syst(p):
    lin, robust = _a0_varianten(p)
    s1, s2 = E("SIGMA_INT")
    b0 = lin[("standard", s1)]
    sY = max(abs(lin[(k, s1)] - b0) for k in ("Upsilon_0.4", "Upsilon_0.6"))
    return math.sqrt(sY**2 + (lin[("standard", s2)] - b0)**2 + (robust - b0)**2)


_a("03_sparc_a0_profil", (D03, "mit Abstand, Neigung und Υ_d als Störgrößen"), "Zahl")
_AS = "An den SPARC-Rotationskurven (149"
_H_A0 = _SP + ("SIGMA_INT",)
_p("03_sparc_a0_profil_wert", "03", "03_sparc_a0_profil", "1,06", lambda p: _a0(p, E("SIGMA_INT")[0])[0] / 1e-10,
   "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0)
_p("03_sparc_a0_profil_stat", "03", "03_sparc_a0_profil", "0,085", lambda p: _a0(p, E("SIGMA_INT")[0])[1] / 1e-10,
   "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0)
_p("03_sparc_a0_profil_syst", "03", "03_sparc_a0_profil", "0,11", lambda p: _a0_syst(p) / 1e-10, "10⁻¹⁰ m/s²", _AS,
   haengt_an=_H_A0)
_p("03_sparc_a0_varianten_unten", "03", "03_sparc_a0_profil", "0,96", lambda p: min(_a0_varianten(p)[0].values()) / 1e-10,
   "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0)
_p("03_sparc_a0_varianten_oben", "03", "03_sparc_a0_profil", "1,19", lambda p: max(_a0_varianten(p)[0].values()) / 1e-10,
   "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0)
# dieselbe Spanne mit Υ_b mitskaliert (Lesart des Archivs, Wahl): das Archiv gab 1,00 bis 1,13; mit Υ_b = 0,7 fest
# (Setzung) liegt sie höher (Prüfzahlen oben)

_a("03_sparc_a0_drittel", (D03, "Gemeinsam je Drittel der Flächendichte angepasst"), "Zahl")
_DR = "Gemeinsam je Drittel der"
_p("03_sparc_a0_drittel_unten", "03", "03_sparc_a0_drittel", "0,75", lambda p: K.a0_je_drittel(_satz(), p.Upsilon_d)[0] / 1e-10,
   "10⁻¹⁰ m/s²", _DR, haengt_an=_SP)
_p("03_sparc_a0_drittel_mitte", "03", "03_sparc_a0_drittel", "1,51", lambda p: K.a0_je_drittel(_satz(), p.Upsilon_d)[1] / 1e-10,
   "10⁻¹⁰ m/s²", _DR, haengt_an=_SP)
_p("03_sparc_a0_drittel_oben", "03", "03_sparc_a0_drittel", "1,06·10⁻¹⁰", lambda p: K.a0_je_drittel(_satz(), p.Upsilon_d)[2],
   "m/s²", _DR, haengt_an=_SP)
_p("03_sparc_a0_drittel_anstieg", "03", "03_sparc_a0_drittel", "0,09", lambda p: K.a0_je_drittel(_satz(), p.Upsilon_d)[3],
   "1", _DR, haengt_an=_SP)
_p("03_sparc_a0_drittel_fehler", "03", "03_sparc_a0_drittel", "0,08", lambda p: K.a0_je_drittel(_satz(), p.Upsilon_d)[4],
   "1", _DR, haengt_an=_SP)


def _senkung(p, lage, gewichtet):
    return 100 * (1 - K.a0_senkung_aussenfeld(_satz(), p.Upsilon_d, E("GE_003") * p.a0, lage, gewichtet, p.a0) / p.a0)


_a("03_sparc_aussenfeld_senkung", (D03, "Das äußere Feld senkt den gelesenen Wert"), "Zahl")
_AF = "Das äußere Feld senkt den"
_H_AF = _SP + ("GE_003",)
_p("03_aussenfeld_senkrecht_gew", "03", "03_sparc_aussenfeld_senkung", "7", lambda p: _senkung(p, "senkrecht", True),
   "%", _AF, haengt_an=_H_AF)
_p("03_aussenfeld_senkrecht_ungew", "03", "03_sparc_aussenfeld_senkung", "14", lambda p: _senkung(p, "senkrecht", False),
   "%", _AF, haengt_an=_H_AF)
_p("03_aussenfeld_gleich_gew", "03", "03_sparc_aussenfeld_senkung", "22", lambda p: _senkung(p, "gleich", True),
   "%", _AF, haengt_an=_H_AF)
_p("03_aussenfeld_gleich_ungew", "03", "03_sparc_aussenfeld_senkung", "31", lambda p: _senkung(p, "gleich", False),
   "%", _AF, haengt_an=_H_AF)

# ================================================================== 03: Vorhersagen und Offenes (SPARC, Sonnenort, Warp)
_a("03_sigma_dyn_rest", (D03, "Zentrale Flächendichte aus der Rotationskurve"), "Zahl")
_ZF = "| Zentrale Flächendichte"
_p("03_sigma_dyn_rest_mittel", "03", "03_sigma_dyn_rest", "0,05", lambda p: -K.sigma_dyn_rest(_tafel(), p.Upsilon_d, p.a0)[0],
   "dex", _ZF, haengt_an=_SP)
_p("03_sigma_dyn_rest_ups04", "03", "03_sigma_dyn_rest", "0,01", lambda p: -K.sigma_dyn_rest(_tafel(), p.Upsilon_d, p.a0)[0],
   "dex", _ZF, haengt_an=_SP, satz="Upsilon_0.4")
_p("03_sigma_dyn_rest_ups06", "03", "03_sigma_dyn_rest", "0,09", lambda p: -K.sigma_dyn_rest(_tafel(), p.Upsilon_d, p.a0)[0],
   "dex", _ZF, haengt_an=_SP, satz="Upsilon_0.6")
_p("03_sigma_dyn_rest_anstieg", "03", "03_sigma_dyn_rest", "0,05", lambda p: K.sigma_dyn_rest(_tafel(), p.Upsilon_d, p.a0)[2],
   "1", _ZF, haengt_an=_SP)
_p("03_sigma_dyn_rest_anstieg_fehler", "03", "03_sigma_dyn_rest", "0,025",
   lambda p: K.sigma_dyn_rest(_tafel(), p.Upsilon_d, p.a0)[3], "1", _ZF, haengt_an=_SP)

_a("03_rest_gasanteil", (D03, "Die Regel hängt nicht am Zustand der Materie"), "Zahl")
_GA = "| Die Regel hängt nicht"
_p("03_rest_gasanteil_unten", "03", "03_rest_gasanteil", "−0,015", lambda p: K.rest_gasanteil(_tafel(), p.Upsilon_d, 0.0, p.a0)[0],
   "dex", _GA, haengt_an=_SP, satz="Upsilon_0.4")
_p("03_rest_gasanteil_oben", "03", "03_rest_gasanteil", "0,018",
   lambda p: K.rest_gasanteil(_tafel(), p.Upsilon_d, E("SIG_GAS_8") * 1e3, p.a0)[0], "dex", _GA,
   haengt_an=_SP + ("SIG_GAS_8",), satz="Upsilon_0.6")
_p("03_rest_gasanteil_fehler", "03", "03_rest_gasanteil", "0,013", lambda p: K.rest_gasanteil(_tafel(), p.Upsilon_d, 0.0, p.a0)[1],
   "dex", _GA, haengt_an=_SP, satz="Upsilon_0.4")


def _jenseits(p):
    return K.rest_jenseits(_tafel(), p.Upsilon_d, E("SIG_GAS_8") * 1e3, math.log(K.q_KOSM)**2 * p.a0, p.a0)


_a("03_rest_jenseits_rq", (D03, "jenseits folgen die SPARC-Kurven der Regel"), "Zahl")
_H_RQ = _SP + ("SIG_GAS_8", "OMEGA_B_H2", "OMEGA_C_H2")
_p("03_rest_jenseits_rq_wert", "03", "03_rest_jenseits_rq", "−0,02", lambda p: _jenseits(p)[0], "dex", "- **Naht:** Bis",
   haengt_an=_H_RQ)
_p("03_rest_jenseits_rq_fehler", "03", "03_rest_jenseits_rq", "0,05", lambda p: _jenseits(p)[1], "dex", "- **Naht:** Bis",
   haengt_an=_H_RQ)


def _reststeigung(p):
    return K.reststeigung_aussen(_tafel(), p.Upsilon_d, E("SIG_GAS_8") * 1e3, E("GE_003") * p.a0, "senkrecht",
                                 E("SPARC_RAND_AUSWAHL") * p.a0, p.a0)


_a("03_reststeigung_aussen", (D03, "Reststeigung"), "Zahl")
_H_RS = _SP + ("SIG_GAS_8", "GE_003", "SPARC_RAND_AUSWAHL")
_p("03_reststeigung_wert", "03", "03_reststeigung_aussen", "−0,30", lambda p: _reststeigung(p)[0], "1", "Ob das mitatmende",
   haengt_an=_H_RS)
_p("03_reststeigung_fehler", "03", "03_reststeigung_aussen", "0,07", lambda p: _reststeigung(p)[1], "1", "Ob das mitatmende",
   haengt_an=_H_RS)

_a("03_sattel_sonne_erde", (D03, "Am Sattelpunkt Sonne–Erde"), "Zahl")
_SA = "| Am Sattelpunkt Sonne–Erde"
_p("03_sattel_abstand", "03", "03_sattel_sonne_erde", "2,59·10⁸", lambda p: K.sattelpunkt(K.GM_sun, K.GM_earth, K.AU), "m",
   _SA, haengt_an=("GM_sun", "GM_earth"))
_p("03_sattel_feld_galaxie", "03", "03_sattel_sonne_erde", "1,34",
   lambda p: K.a_N_aus_kreisbahn(M("V0_SONNE") * 1e3, M("R0_SONNE") * K.kpc, p.a0) / p.a0, "a₀", _SA,
   haengt_an=("V0_SONNE", "R0_SONNE"), satz="Sonnenkreis_gemessen")
_p("03_sattel_ungeglaettet", "03", "03_sattel_sonne_erde", "7·10⁻¹¹", lambda p: K.bose_stufe_max(2)[2] * p.a0, "m/s²", _SA)

_a("03_warp_aufwicklung", (D03, "Die freie Präzession wickelt den Warp zudem auf"), "Zahl")


def _aufwicklung(p, form):
    """Zeit [Mrd. Jahre], bis die Knotenlinien der Ringe bei 12 und 16 kpc um 180° auseinanderlaufen (Massenmodell)."""
    q, f = _milchstrasse(p)
    s = math.sqrt(p.a0 / K.kpc)                                             # Frequenzeinheit des Lösers
    w = []
    for R in E("WARP_RINGE"):
        rho0 = float(q["rho_punkt"](R, 0.0))
        Om, oz, g = K.frequenzen_loesung(f, rho0, R)
        nu_ = 1.0 if form == "newton" else float(K.nu(g, 1.0))
        w.append(K.knotenpraezession_frequenzen(Om, oz, nu_, "ort") * s)
    return K.warp_aufwicklung(w[0], w[1], math.radians(E("WARP_AUSEINANDER"))) / (1e9 * K.yr)


_H_WARP = _H_MCM + ("WARP_RINGE", "WARP_AUSEINANDER")
_p("03_warp_aufwicklung_regel", "03", "03_warp_aufwicklung", "0,26", lambda p: _aufwicklung(p, "ort"), "Mrd. Jahre",
   "Die freie Präzession wickelt", haengt_an=_H_WARP)
_p("03_warp_aufwicklung_newton", "03", "03_warp_aufwicklung", "0,28", lambda p: _aufwicklung(p, "newton"), "Mrd. Jahre",
   "Die freie Präzession wickelt", haengt_an=_H_WARP)


def _df2(p, lage):
    ge = 0.0 if lage == "ohne" else K.neigung(E("M_NGC1052") * K.M_sun, M("DF2_ABSTAND") * K.kpc)
    return K.sigma_plummer_jeans(M("DF2_M") * K.M_sun, M("DF2_RE") * K.kpc, ge, p.a0,
                                 lage="kugel" if lage == "ohne" else lage) / 1e3


_a("03_druck_lage_aussenfeld", (D03, "DF2 und DF4 je nach Lage des Außenfelds"), "Zahl")
_H_DF2 = ("DF2_M", "DF2_RE", "DF2_ABSTAND", "M_NGC1052", "GM_sun")
_p("03_druck_df2_gleichgerichtet", "03", "03_druck_lage_aussenfeld", "17", lambda p: _df2(p, "gleich"), "km/s",
   "- Druckgestützte Systeme", haengt_an=_H_DF2)
_p("03_druck_df2_ohne_feld", "03", "03_druck_lage_aussenfeld", "21", lambda p: _df2(p, "ohne"), "km/s",
   "- Druckgestützte Systeme", haengt_an=_H_DF2)

_a("03_gegenlaeufig", (D03, "Gegenläufige Sternscheiben derselben Galaxie"), "Zahl")


def _gegenlaeufig(p):
    """v_gegen/v_mit bei 3r* einer Punktmasse (a_N = a₀·(r*/r)²), wenn nur Rückläufe im Drehsinn zählten: Newton
    gegen die Regel, √(a_N/a)."""
    aN = p.a0 / E("R_3RSTERN")**2
    return math.sqrt(aN / K.a_regel(aN, p.a0))


_p("03_gegenlaeufig_drehsinn", "03", "03_gegenlaeufig", "0,53", _gegenlaeufig, "1", "| Gegenläufige Sternscheiben",
   haengt_an=("R_3RSTERN",))

# ================================================================== Zweite Fundstellen derselben Größe
_w(D07, "0,53", "In Scheiben innen", "07_feldform_viertel_tief")
_w(D07, "0,82", "In Scheiben innen", "07_feldform_viertel_hoch")
_w(D07, "0,79", "In Scheiben innen", "07_feldform_rd_tief")
_w(D07, "0,92", "In Scheiben innen", "07_feldform_rd_hoch")
_w(D07, "3,4", "In Scheiben innen", "07_feldform_gleichstand")
_w(D03, "3,4", "| Scheiben innen", "07_feldform_gleichstand")
_w(D03, "0,03", "Gemessen (SPARC, mit Fehlern", "07_gewicht_f_unten")
_w(D03, "0,13", "Gemessen (SPARC, mit Fehlern", "07_gewicht_f_oben")
_w(D03, "0,00", "Gemessen (SPARC, mit Fehlern", "07_gewicht_f_ohne_ic2574")
_w(D07, "0,00", "Das positive f", "07_gewicht_f_ohne_ic2574")
_w(D07, "0,4", "An den SPARC-Rotationskurven, mit", "UPS_D_04")
_w(D07, "0,6", "An den SPARC-Rotationskurven, mit", "UPS_D_06")
_w(D07, "0,4", "Die Grenzwerte 0 und 1 von", "UPS_D_04")
_w(D07, "0,6", "Die Grenzwerte 0 und 1 von", "UPS_D_06")
_w(D07, "69", "In der Mittelebene ist das", "07_dynamikmasse_ort_schicht")
_w(D07, "70", "In der Mittelebene ist das", "07_dynamikmasse_ort_ebene")
_w(D03, "69", "| Beschleunigung am Sonnenort", "07_dynamikmasse_ort_schicht")
_w(D03, "70", "| Beschleunigung am Sonnenort", "07_dynamikmasse_ort_ebene")
_w(D03, "0,09", "- Der Anstieg der Schwelle", "03_sparc_a0_drittel_anstieg")
_w(D03, "0,08", "- Der Anstieg der Schwelle", "03_sparc_a0_drittel_fehler")
_w(D03, "0,03", "Mit a₀ ≤ 1,13·10⁻¹⁰ aus", "GE_003")
_w(D03, "5,5", "das Fenster 5 ≤ Q ≤ 5,5", "03_q_max_wuerfel")
_w(D03, "0,7", "Υ_b = 0,7 fest", "SPARC_UPS_B")
