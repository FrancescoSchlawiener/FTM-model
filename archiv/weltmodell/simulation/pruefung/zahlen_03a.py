# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 03, Teil a: Kopftabelle, „Galaxie, Haufen, Filamente“, „Dunkle Materie: geänderte Bewegungsregel“,
„Haufen und Hintergrundstrahlung: die Schließungsregel“. Aufbau wie zahlen_01.py; Kennungen beginnen mit "03a_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
SPARC: Auswahl, Υ_b und Helium aus SPARC_SATZ; Υ_d aus dem Parametersatz. Zufall nur mit festem Seed je Aufruf.
Große SPARC-Lesungen (Monte Carlo über die Messfehler) rechnet _einmal() je Lauf; wer sie liest, nennt ihre Formeln in
haengt_an, weil die Aufrufverfolgung einen gespeicherten Wert nicht sieht.
"""
import math

import numpy as np

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE

from . import zahlen_03  # noqa: F401  (Aussage 03_milchstrasse_punktmasse)
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "03"
D = datei(B)
_a = AUSSAGEN.neu
_p = PRUEFZAHLEN.neu
_S = K.SPARC_SATZ


def E(k):
    """Wert einer Eingabe."""
    return EINGABEN[k].wert


def M(k):
    """Masse einer Eingabe in M☉ als kg."""
    return EINGABEN[k].wert * K.M_sun


MW = "M_MILCHSTRASSE"
_GM = ("GM_sun",)
_Q = ("OMEGA_B_H2", "OMEGA_C_H2")


# ================================================================== Kopftabelle
_a("03a_tabelle", (D, "| System | baryonische Masse | Δs | Δs/R |"), "Zahl")
_p("03a_ds_milchstrasse", B, "03a_tabelle", "197,4", lambda p: K.ds(M(MW)) / K.AU, "AE", "| Milchstraße |",
   haengt_an=(MW,) + _GM)
_p("03a_dsR_milchstrasse", B, "03a_tabelle", "6,3·10⁻⁸", lambda p: K.ds(M(MW)) / (E("R_MILCHSTRASSE_RAND") * K.kpc),
   "1", "| Milchstraße |", haengt_an=(MW, "R_MILCHSTRASSE_RAND") + _GM)
_p("03a_ds_haufen", B, "03a_tabelle", "8,8", lambda p: K.ds(M("M_HAUFEN")) / K.ly, "Lj", "| Galaxienhaufen |",
   haengt_an=("M_HAUFEN",) + _GM)
_p("03a_dsR_haufen", B, "03a_tabelle", "1,4·10⁻⁶", lambda p: K.ds(M("M_HAUFEN")) / (E("R_HAUFEN") * K.Mpc), "1",
   "| Galaxienhaufen |", haengt_an=("M_HAUFEN", "R_HAUFEN") + _GM)


def _masse_horizont():
    """Baryonen in der Kugel des eigenen Teilchenhorizonts heute: ρ_b·(4π/3)·d_H³."""
    return K.rho_b_mitbewegt() * 4 / 3 * math.pi * K.teilchenhorizont(0.0) ** 3


_p("03a_masse_universum", B, "03a_tabelle", "7,4·10²²", lambda p: _masse_horizont() / K.M_sun, "M☉",
   "| beobachtbares Universum", haengt_an=_GM)
_p("03a_ds_universum", B, "03a_tabelle", "3,9", lambda p: K.ds(_masse_horizont()) / K.ly / 1e9, "Mrd Lj",
   "| beobachtbares Universum")
_p("03a_dsR_universum", B, "03a_tabelle", "0,083", lambda p: K.rueckstellbilanz(0.0), "1", "| beobachtbares Universum")

# ================================================================== Galaxie, Haufen, Filamente
_a("03a_straffung", (D, "in der Milchstraße fällt sie vom Zentrum"), "Zahl")
_p("03a_straffung_rand", B, "03a_straffung", "6·10⁻⁸", lambda p: K.ds(M(MW)) / (E("R_MILCHSTRASSE_RAND") * K.kpc),
   "1", "bis zum Rand", haengt_an=(MW, "R_MILCHSTRASSE_RAND") + _GM)

_a("03a_flussdimension", (D, "d_eff = 1 + (d − 1)(1 − χ/2)"), "Zahl")
_p("03a_deff_kugel", B, "03a_flussdimension", "2", lambda p: K.flussdimension(1e-9 * p.a0, p.a0, 3), "1",
   "Die Kugel geht von 3 auf 2")
_p("03a_deff_filament", B, "03a_flussdimension", "1,5", lambda p: K.flussdimension(1e-9 * p.a0, p.a0, 2), "1",
   "ein Filament von 2 auf 1,5")

_a("03a_filament", (D, "r*_fil = 2Gμ/a₀"), "Zahl")
_p("03a_r_stern_filament", B, "03a_filament", "2,5",
   lambda p: K.r_stern_filament(M("MU_FILAMENT") / K.Mpc, p.a0) / K.kpc, "kpc", "r*_fil = 2Gμ/a₀",
   haengt_an=("MU_FILAMENT",) + _GM)


def _filament_querschnitt(p):
    """Schrumpfrate des Querschnitts (in H₀), bei der w = ½: θ̃ = H₀(1 − 2x) aus Dehnung längs H₀ und quer −xH₀."""
    from scipy.optimize import brentq
    H = K.H0_si(p.H0)
    th = brentq(lambda t: K.gewicht_w(t * H, p.a0) - 0.5, 1e-9, 1.0, xtol=1e-15)
    return (1 - th) / 2


_p("03a_filament_querschnitt", B, "03a_filament", "0,46", _filament_querschnitt, "H₀", "schneller als 0,46·H₀")

# ================================================================== Dunkle Materie: Punktmasse, Streuung
_a("03a_punktmasse", (D, "Milchstraße als Punktmasse: v = 250 km/s bei 5 kpc"), "Zahl")
for _k, _t, _r in (("5", "250", "R_5KPC"), ("15", "198", "R_15KPC"), ("50", "180", "R_50KPC")):
    _p(f"03a_v_milchstrasse_{_k}kpc", B, "03a_punktmasse", _t,
       lambda p, r=_r: K.v_kreis_punkt(M(MW), E(r) * K.kpc, p.a0) / 1e3, "km/s", "Milchstraße als Punktmasse",
       haengt_an=(MW, _r) + _GM)
_p("03a_v_milchstrasse_rstern", B, "03a_punktmasse", "216",
   lambda p: K.v_kreis_punkt(M(MW), K.r_stern(M(MW), p.a0), p.a0) / 1e3, "km/s", "Milchstraße als Punktmasse",
   haengt_an=(MW,) + _GM)

_a("03a_streuung_tief", (D, "σ_v⁴ = (4/81)·GMa₀"), "Zahl")
_p("03a_sigma_1e7", B, "03a_streuung_tief", "9,2", lambda p: K.sigma_v_tief(M("M_KUGEL_KLEIN"), p.a0) / 1e3, "km/s",
   "σ_v⁴ = (4/81)·GMa₀", haengt_an=("M_KUGEL_KLEIN",) + _GM)
_p("03a_sigma_2e8", B, "03a_streuung_tief", "19,5", lambda p: K.sigma_v_tief(M("M_KUGEL_GROSS"), p.a0) / 1e3, "km/s",
   "σ_v⁴ = (4/81)·GMa₀", haengt_an=("M_KUGEL_GROSS",) + _GM)


# ================================================================== Universalität (SPARC)
def _tafel():
    return K.sparc_tafel(_S["q_max"], _S["inc_min"])


def _reste(p, sig=0.0, ge=0.0, teil="alle", ups_d=None):
    """Rest je Galaxie gegen die Regel mit den eigenen Baryonen; sig in m/s, ge in a₀ (senkrecht)."""
    u = p.Upsilon_d if ups_d is None else ups_d
    return np.array([K.rest_galaxie(g["R"], g["Vo"], g["Vg"], g["Vd"], g["Vb"], g["Rd"], u, sig, ge * p.a0, p.a0,
                                    teil=teil) for g in _tafel()])


def _gasanteil(p):
    """Gasanteil f_He·M_HI/(Υ_d·L + f_He·M_HI) je Galaxie der Tafel."""
    return np.array([_S["helium"] * g["MHI"] / (p.Upsilon_d * g["L"] + _S["helium"] * g["MHI"]) for g in _tafel()])


GASREICH = 0.5    # Festlegung „gasreich“: Gasanteil über ½ (ein Symbol für den Gasanteil fehlt im Register)


def _rang_gasreich(p, ge):
    """Spearman ρ und p von Rest und log SB₀ in gasreichen Galaxien, Gasdruck σ = 8 km/s, äußeres Feld ge (a₀)."""
    from scipy.stats import spearmanr
    r = _reste(p, E("SIG_GAS_8") * 1e3, ge)
    lsb = np.log10([g["SB0"] for g in _tafel()])
    m = (_gasanteil(p) > GASREICH) & np.isfinite(r)
    return spearmanr(lsb[m], r[m])


def _innen_aussen(p):
    d = _reste(p, teil="innen") - _reste(p, teil="aussen")
    d = d[np.isfinite(d)]
    return d.mean(), d.std() / math.sqrt(d.size)


_a("03a_universal", (D, "hängt der Rest je Galaxie an keiner der geprüften Größen"), "Zahl")
_p("03a_sparc_galaxien", B, "03a_universal", "149", lambda p: float(np.isfinite(_reste(p)).sum()), "1",
   "SPARC, 149 Galaxien")
_p("03a_upsilon_d", B, "03a_universal", "0,5", lambda p: p.Upsilon_d, "M☉/L☉[3,6]", "SPARC, 149 Galaxien")
_p("03a_innen_aussen", B, "03a_universal", "−0,02", lambda p: _innen_aussen(p)[0], "dex", "innen und außen")
_p("03a_innen_aussen_fehler", B, "03a_universal", "0,013", lambda p: _innen_aussen(p)[1], "dex", "innen und außen")
_p("03a_rang_sb_gasreich", B, "03a_universal", "0,36", lambda p: _rang_gasreich(p, 0.0)[0], "1",
   "Rangkorrelation 0,36", haengt_an=("SIG_GAS_8",))
_p("03a_rang_sb_ge002", B, "03a_universal", "0,28", lambda p: _rang_gasreich(p, E("GE_002"))[0], "1",
   "senkt das auf 0,28–0,23", haengt_an=("SIG_GAS_8", "GE_002"))
_p("03a_rang_sb_ge003", B, "03a_universal", "0,23", lambda p: _rang_gasreich(p, E("GE_003"))[0], "1",
   "senkt das auf 0,28–0,23", haengt_an=("SIG_GAS_8", "GE_003"))
_p("03a_rang_sb_ge003_p", B, "03a_universal", "0,07", lambda p: _rang_gasreich(p, E("GE_003"))[1], "1",
   "nicht mehr gesichert, p = 0,07", haengt_an=("SIG_GAS_8", "GE_003"))

_SPEICHER = {}
_F_REST = ("sparc_tafel", "sparc_kurve", "sparc_tabelle", "rest_galaxie", "v_bar_quadrat", "g_gemessen", "nu",
           "streuung_messfehler")
_F_SCHWELLE = _F_REST + ("a0_je_galaxie",)
N_MC_REST, N_MC_SCHWELLE, N_BOOT, SEED = 200, 100, 500, 20260930


def _einmal(p, name):
    """Monte Carlo über die Messfehler je Galaxie (eigener Seed je Galaxie), einmal je Lauf gerechnet:
    rest(σ_Υ) = (Reste σ = 8 km/s, Streuung je Galaxie aus den Fehlern); schwelle(σ_Υ) = (log a₀ je Galaxie, Streuung)."""
    schluessel = (p.kennung,) + name
    if schluessel in _SPEICHER:
        return _SPEICHER[schluessel]
    art, sups = name
    sig = E("SIG_GAS_8") * 1e3
    u = p.Upsilon_d
    if art == "rest":
        def lesung(R, Vo, Vg, Vd, Vb, Rd):
            return K.rest_galaxie(R, Vo, Vg, Vd, Vb, Rd, u, sig, 0.0, p.a0)
        n = N_MC_REST
    else:
        def lesung(R, Vo, Vg, Vd, Vb, Rd):
            return math.log10(K.a0_je_galaxie(R, Vo, Vg, Vd, Vb, Rd, u, sig))
        n = N_MC_SCHWELLE
    gem, sim = [], []
    for i, g in enumerate(_tafel()):
        x = lesung(g["R"], g["Vo"], g["Vg"], g["Vd"], g["Vb"], g["Rd"])
        if not np.isfinite(x):
            continue
        gem.append(x)
        sim.append(K.streuung_messfehler(g, lesung, n, SEED + i, sups))
    _SPEICHER[schluessel] = (np.array(gem), np.array(sim))
    return _SPEICHER[schluessel]


def _streuung_regel(p):
    """Streuung der Regel selbst: √(Var(Rest) − ⟨σ_Fehler²⟩) bei 0,1 dex Streuung von Υ."""
    g, s = _einmal(p, ("rest", _S["ups_streuung_dex"]))
    return math.sqrt(max(g.var() - np.mean(s**2), 0.0))


_p("03a_streuung_rest", B, "03a_universal", "0,14", lambda p: _einmal(p, ("rest", 0.0))[0].std(), "dex",
   "Die Streuung des Rests (0,14 dex)", haengt_an=("SIG_GAS_8",) + _F_REST)
_p("03a_streuung_regel", B, "03a_universal", "0,05", _streuung_regel, "dex", "Die Streuung des Rests (0,14 dex)",
   haengt_an=("SIG_GAS_8",) + _F_REST)

# ================================================================== Einzeltrichter oder kollektiv
_a("03a_kollektiv", (D, "Kompakte Kugelsternhaufen (g/a₀ ≈ 25)"), "Zahl")
_p("03a_nu_kugelsternhaufen", B, "03a_kollektiv", "1,007", lambda p: K.nu(E("G_KUGELSTERNHAUFEN") * p.a0, p.a0), "1",
   "Kompakte Kugelsternhaufen", haengt_an=("G_KUGELSTERNHAUFEN",))
_p("03a_nu_zwerg", B, "03a_kollektiv", "3,3", lambda p: K.nu(E("G_ZWERG") * p.a0, p.a0), "1",
   "Kompakte Kugelsternhaufen", haengt_an=("G_ZWERG",))
_p("03a_grenzflaechendichte", B, "03a_kollektiv", "251", lambda p: K.sigma_kollektiv(p.a0), "M☉/pc²",
   "Σ < a₀/(πG)", haengt_an=_GM)

# ================================================================== Umgebung
_a("03a_umgebung", (D, "Die Regel wirkt nicht additiv pro Teilchen"), "Zahl")


def _anteil_im_haufen(p, d):
    """v im Haufen / v frei bei 80 kpc: Galaxie Hernquist, Haufen isotherm (M ∝ r bis HAUFEN_R_ISO), Feld senkrecht."""
    r = E("R_80KPC") * K.kpc
    gi = K.G * K.hernquist_masse(M(MW), E("HERNQUIST_A") * K.kpc, r) / r**2
    D_ = E(d) * K.Mpc
    ge = K.G * M("M_HAUFEN_UMGEBUNG") * min(1.0, D_ / (E("HAUFEN_R_ISO") * K.Mpc)) / D_**2
    return math.sqrt(K.nu_aussenfeld(gi, ge, "senkrecht", p.a0) / K.nu(gi, p.a0)) * 100


_H = (MW, "R_80KPC", "HERNQUIST_A", "M_HAUFEN_UMGEBUNG", "HAUFEN_R_ISO") + _GM
_p("03a_haufen_nah", B, "03a_umgebung", "50", lambda p: _anteil_im_haufen(p, "D_HAUFEN_NAH"), "%",
   "des freien Werts", haengt_an=_H + ("D_HAUFEN_NAH",))
_p("03a_haufen_fern", B, "03a_umgebung", "77", lambda p: _anteil_im_haufen(p, "D_HAUFEN_FERN"), "%",
   "des freien Werts", haengt_an=_H + ("D_HAUFEN_FERN",))


def _a_i50():
    return K.G * M(MW) / (E("R_50KPC") * K.kpc) ** 2


_HM = (MW, "R_50KPC", "GE_003") + _GM
_p("03a_nu_frei_50", B, "03a_umgebung", "6,25", lambda p: K.nu(_a_i50(), p.a0), "1", "ν = 6,25 frei", haengt_an=_HM)
_p("03a_nu_senkrecht_50", B, "03a_umgebung", "5,36",
   lambda p: K.nu_aussenfeld(_a_i50(), E("GE_003") * p.a0, "senkrecht", p.a0), "1", "5,36 senkrecht", haengt_an=_HM)
_p("03a_nu_gleich_50", B, "03a_umgebung", "4,59",
   lambda p: K.nu_aussenfeld(_a_i50(), E("GE_003") * p.a0, "gleich", p.a0), "1", "4,59 gleichgerichtet",
   haengt_an=_HM)


def _plummer(p):
    """A/g_e der Plummer-Kugel mit GM/b² = a₀ bei g_e = 0,1·a₀ (A hängt nur an GM/b² und g_e)."""
    b = E("R_E_PLUMMER") * K.kpc
    ge = E("GE_01") * p.a0
    return K.massenmittel_A_plummer(p.a0 * b * b / K.G, b, ge, p.a0) / ge


_p("03a_plummer_A", B, "03a_umgebung", "2,37", _plummer, "g_e", "das 2,37-Fache von g_e",
   haengt_an=("R_E_PLUMMER", "GE_01"))
_p("03a_nu_ge01", B, "03a_umgebung", "3,69", lambda p: K.nu(E("GE_01") * p.a0, p.a0), "1", "ν(g_e) = 3,69",
   haengt_an=("GE_01",))


def _A_scheibe(p):
    return K.massenmittel_A_scheibe(M(MW), E("RD_MILCHSTRASSE") * K.kpc, E("GE_003") * p.a0, p.a0)


_HS = (MW, "RD_MILCHSTRASSE", "GE_003") + _GM
_p("03a_scheibe_A", B, "03a_umgebung", "1,27", lambda p: _A_scheibe(p) / (E("GE_003") * p.a0), "g_e",
   "A = 1,27·g_e", haengt_an=_HS)


def _m1_scheibe(p, r):
    R = E(r) * K.kpc
    gi = float(K.g_expscheibe(M(MW), E("RD_MILCHSTRASSE") * K.kpc, R))
    return K.einseitigkeit_m1(gi, R, E("GE_003") * p.a0, _A_scheibe(p), p.a0) * 100


def _m1_punkt(p, r, bezug):
    R, ge = E(r) * K.kpc, E("GE_003") * p.a0
    gi = K.G * M(MW) / R**2
    ref = ge if bezug == "g_e" else K.nu(ge, p.a0) * ge
    return K.einseitigkeit_m1(gi, R, ge, ref, p.a0) * 100


_p("03a_m1_scheibe_30", B, "03a_umgebung", "4,4",
   lambda p: max(_m1_scheibe(p, r) for r in ("R_10KPC", "R_20KPC", "R_30KPC")), "%", "um höchstens 4,4 %",
   haengt_an=_HS + ("R_10KPC", "R_20KPC", "R_30KPC"))
_p("03a_m1_scheibe_50", B, "03a_umgebung", "31", lambda p: _m1_scheibe(p, "R_50KPC"), "%", "bei 50 kpc um 31 %",
   haengt_an=_HS + ("R_50KPC",))
_p("03a_m1_punkt_20_ge", B, "03a_umgebung", "2", lambda p: _m1_punkt(p, "R_20KPC", "g_e"), "%", "2–36 % mit g_e",
   haengt_an=_HM + ("R_20KPC",))
_p("03a_m1_punkt_50_ge", B, "03a_umgebung", "36", lambda p: _m1_punkt(p, "R_50KPC", "g_e"), "%", "2–36 % mit g_e",
   haengt_an=_HM)
_p("03a_m1_punkt_20_nu", B, "03a_umgebung", "13", lambda p: _m1_punkt(p, "R_20KPC", "nu"), "%",
   "13–20 % mit ν(g_e)·g_e", haengt_an=_HM + ("R_20KPC",))
_p("03a_m1_punkt_50_nu", B, "03a_umgebung", "20", lambda p: _m1_punkt(p, "R_50KPC", "nu"), "%",
   "13–20 % mit ν(g_e)·g_e", haengt_an=_HM)

# ================================================================== Die Schwelle
_a("03a_schwelle", (D, "Heute ist a₀ die Temperatur des Hubble-Horizonts"), "Zahl")
_p("03a_theta_heute", B, "03a_schwelle", "1,06", lambda p: K.theta_(p.a0, K.hubble_rate(p)), "1", "(Θ = 1,06)")


def _z_theta_unten(p):
    """z, bis zu dem Θ(z) = 2πa₀/(cH(z)) nicht unter 1 − 10 % fällt."""
    from scipy.optimize import brentq
    return brentq(lambda z: K.theta_(p.a0, K.H0_si(p.H0) * K.E_hintergrund(z)) - E("THETA_UNTEN"), 0.0, 3.0,
                  xtol=1e-12)


_p("03a_z_theta", B, "03a_schwelle", "0,3", _z_theta_unten, "1", "innerhalb von 10 % nur für z ≲ 0,3",
   haengt_an=("THETA_UNTEN",))
_p("03a_lambda_horizont", B, "03a_schwelle", "14", lambda p: K.schliesslaenge(p.a0) / K.Gpc, "Gpc",
   "heute Λ ≈ πc/H₀ (12–14 Gpc)", satz="horizont")
_p("03a_lambda_literatur", B, "03a_schwelle", "12", lambda p: K.schliesslaenge(K.A0_LIT) / K.Gpc, "Gpc",
   "heute Λ ≈ πc/H₀ (12–14 Gpc)", haengt_an=("A0_LIT",))


def _vakuum_theta1(p):
    """Vakuumteil bei Θ = 1 (a₀ = cH/2π der Epoche; Satz "horizont")."""
    return K.omega_vakuum(p.a0, p.H0)


_p("03a_vakuum_atmend", B, "03a_schwelle", "22", lambda p: _vakuum_theta1(p) * 100, "%", "zu jeder Zeit bei 22 %", satz="horizont")
_p("03a_delta_neff", B, "03a_schwelle", "2,09", lambda p: K.delta_N_eff(_vakuum_theta1(p)), "1", "ΔN_eff = 2,09", satz="horizont")
_p("03a_vakuum_kruemmung", B, "03a_schwelle", "0,22", _vakuum_theta1, "1", "wie eine Krümmung von 0,22", satz="horizont")
_p("03a_a0_horizont_674", B, "03a_schwelle", "1,04", lambda p: p.a0 * 1e10, "10⁻¹⁰ m/s²",
   "Der heutige Wert cH₀/2π ist 1,04–1,08·10⁻¹⁰ m/s²", satz="horizont")
_p("03a_a0_horizont_70", B, "03a_schwelle", "1,08·10⁻¹⁰", lambda p: K.a0_modell(E("H0_70")), "m/s²",
   "Der heutige Wert cH₀/2π ist 1,04–1,08·10⁻¹⁰ m/s²", haengt_an=("H0_70",))
_p("03a_unter_literatur_70", B, "03a_schwelle", "10", lambda p: (1 - K.a0_modell(E("H0_70")) / K.A0_LIT) * 100, "%",
   "10–13 % unter dem Literaturwert", haengt_an=("H0_70", "A0_LIT"))
_p("03a_unter_literatur_674", B, "03a_schwelle", "13", lambda p: (1 - p.a0 / K.A0_LIT) * 100, "%",
   "10–13 % unter dem Literaturwert", haengt_an=("A0_LIT",), satz="horizont")


def _h0_literatur(p):
    from scipy.optimize import brentq
    return brentq(lambda h: K.a0_modell(h) - K.A0_LIT, 10.0, 200.0, xtol=1e-12)


_p("03a_h0_literatur", B, "03a_schwelle", "77,6", _h0_literatur, "km/s/Mpc", "H₀ = 77,6 km/s/Mpc",
   haengt_an=("A0_LIT",))


def _tully_fisher(p):
    """Äußerster Punkt mit V_bar² > 0 je Galaxie mit V_flach > 0: (a₀ aus V⁴/(G M_b), s dort, a₀ mit der Regel)."""
    from scipy.optimize import brentq
    Vf, Mb, R, aN = [], [], [], []
    for g in _tafel():
        vb2 = K.v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], p.Upsilon_d)
        if g["Vf"] <= 0 or not (vb2[-3:] > 0).any():
            continue
        k = int(np.nonzero(vb2 > 0)[0][-1])
        Vf.append(g["Vf"])
        Mb.append((p.Upsilon_d * g["L"] + _S["helium"] * g["MHI"]) * K.M_sun)
        R.append(g["R"][k])
        aN.append(vb2[k] / g["R"][k])
    Vf, Mb, R, aN = map(np.array, (Vf, Mb, R, aN))
    a0_tf = p.a0 * 10 ** np.median(4 * np.log10(Vf / np.array([K.v_flach(m, p.a0) for m in Mb])))
    s = np.median(K.s_(aN, p.a0))
    la = brentq(lambda x: np.median(np.log10(Vf / np.sqrt(K.a_regel(aN, 10**x) * R))), -10.6, -9.4, xtol=1e-12)
    return a0_tf, s, 10**la


_p("03a_tf_vflach", B, "03a_schwelle", "1,5·10⁻¹⁰", lambda p: _tully_fisher(p)[0], "m/s²",
   "Die Tully-Fisher-Relation mit V_flach liest 1,5·10⁻¹⁰", haengt_an=_GM)
_p("03a_tf_s_median", B, "03a_schwelle", "0,23", lambda p: _tully_fisher(p)[1], "1", "s im Median 0,23",
   haengt_an=_GM)
_p("03a_tf_regel", B, "03a_schwelle", "0,97", lambda p: _tully_fisher(p)[2] * 1e10, "10⁻¹⁰ m/s²",
   "mit der Regel an denselben Punkten gelesen gibt sie 0,97", haengt_an=_GM)

_p("03a_a0_je_galaxie", B, "03a_schwelle", "1,09·10⁻¹⁰", lambda p: 10 ** np.median(_einmal(p, ("schwelle", 0.0))[0]),
   "m/s²", "Je Galaxie gelesen", haengt_an=("SIG_GAS_8",) + _F_SCHWELLE)
_p("03a_a0_streuung_gemessen", B, "03a_schwelle", "0,28",
   lambda p: 1.4826 * np.median(np.abs(_einmal(p, ("schwelle", 0.0))[0] - np.median(_einmal(p, ("schwelle", 0.0))[0]))),
   "dex", "wie die Messfehler es erwarten lassen (0,27–0,30 dex)", haengt_an=("SIG_GAS_8",) + _F_SCHWELLE)
_p("03a_a0_streuung_unten", B, "03a_schwelle", "0,27", lambda p: np.median(_einmal(p, ("schwelle", 0.0))[1]), "dex",
   "wie die Messfehler es erwarten lassen (0,27–0,30 dex)", haengt_an=("SIG_GAS_8",) + _F_SCHWELLE)
_p("03a_a0_streuung_oben", B, "03a_schwelle", "0,30",
   lambda p: np.median(_einmal(p, ("schwelle", _S["ups_streuung_dex"]))[1]), "dex",
   "wie die Messfehler es erwarten lassen (0,27–0,30 dex)", haengt_an=("SIG_GAS_8",) + _F_SCHWELLE)


def _anstieg(p):
    """Theil-Sen-Steigung von log a₀ je Galaxie gegen log Σ₀ (Σ₀ = Υ_d·SB₀), Fehler per Bootstrap über Galaxien."""
    from scipy.stats import theilslopes
    y = []
    x = []
    for g in _tafel():
        a = K.a0_je_galaxie(g["R"], g["Vo"], g["Vg"], g["Vd"], g["Vb"], g["Rd"], p.Upsilon_d, E("SIG_GAS_8") * 1e3)
        if np.isfinite(a):
            x.append(math.log10(p.Upsilon_d * g["SB0"]))
            y.append(math.log10(a))
    x, y = np.array(x), np.array(y)
    rng = np.random.default_rng(SEED)
    bs = [theilslopes(y[i], x[i])[0] for i in (rng.integers(0, x.size, x.size) for _ in range(N_BOOT))]
    return theilslopes(y, x)[0], float(np.std(bs))


_p("03a_anstieg_schwelle", B, "03a_schwelle", "0,06", lambda p: _anstieg(p)[0], "dex/Dekade",
   "steigt nicht messbar mit der Flächendichte", haengt_an=("SIG_GAS_8",))
_p("03a_anstieg_schwelle_fehler", B, "03a_schwelle", "0,04", lambda p: _anstieg(p)[1], "dex/Dekade",
   "steigt nicht messbar mit der Flächendichte", haengt_an=("SIG_GAS_8",))
_p("03a_a0_zelle", B, "03a_schwelle", "1,13·10⁻¹⁰", lambda p: K.a0_aus_zelle(K.KREISPAAR_GRENZE * K.d_LSS()), "m/s²",
   "Mit a₀ ≤ 1,13·10⁻¹⁰ aus der Zelle", haengt_an=_Q + ("KREISPAAR_GRENZE",))
_p("03a_a0_standard", B, "03a_schwelle", "1,1·10⁻¹⁰", lambda p: p.a0, "m/s²", "Das Modell rechnet mit")

# ================================================================== Schließungsregel: mitatmend
_a("03a_mitatmend", (D, "**Mitatmend (w ≈ 1):**"), "Zahl")
_p("03a_rest_w_heute", B, "03a_mitatmend", "1,5·10⁻²", lambda p: 1 - K.gewicht_w(3 * K.hubble_rate(p), p.a0), "1",
   "1 − w = 1,5·10⁻² heute", satz="mitatmend_heute")
_p("03a_rest_w_horizont", B, "03a_mitatmend", "1,3·10⁻²",
   lambda p: 1 - K.gewicht_w(3 * K.hubble_rate(p), p.a0), "1", "1,3·10⁻² mit cH₀/2π", satz="horizont")


def _exponent_rekombination(p):
    """√(θ̃/θ_c) bei der Rekombination: θ̃ = 3H(z), a₀(z) mit dem Strahlungsanteil der inneren Faltung."""
    z = p.z
    a0z = p.a0 * K.a0_frueh_relativ(z)
    return float(K.rueckkehr_exponent(0.0, 0.0, 3 * K.H0_si(p.H0) * K.E_hintergrund(z), a0z))


_p("03a_rest_w_rekombination", B, "03a_mitatmend", "738", _exponent_rekombination, "1", "(e^{−738})",
   haengt_an=("Z_REKOMBINATION",), satz="z1090")
_p("03a_Q_kosmisch", B, "03a_mitatmend", "5,36", lambda p: K.Q_KOSM, "1", "Ω_c/Ω_b = Q = q/(1 − q) = 5,36",
   haengt_an=_Q)
_p("03a_q_kosmisch", B, "03a_mitatmend", "0,843", lambda p: K.q_KOSM, "1", "kosmisch q = 0,843", haengt_an=_Q)
_p("03a_Q_radiant", B, "03a_mitatmend", "5,283", lambda p: K.kac_rueckkehrzeit(2 * math.pi, 1.0) - 1, "1",
   "Q = 2π − 1 = 5,283")
_p("03a_Q_bao", B, "03a_mitatmend", "5,32",
   lambda p: MESSWERTE["OMEGA_C_H2_BAO"].wert / MESSWERTE["OMEGA_B_H2_BAO"].wert, "1",
   "Gemessen (Hintergrundstrahlung mit BAO bzw. allein)", haengt_an=("OMEGA_B_H2_BAO", "OMEGA_C_H2_BAO"))
_p("03a_Q_planck_fehler", B, "03a_mitatmend", "0,065", lambda p: K.Q_fehler(), "1",
   "Gemessen (Hintergrundstrahlung mit BAO bzw. allein)", haengt_an=_Q)
_p("03a_Q_planck_abstand", B, "03a_mitatmend", "1,3",
   lambda p: (K.Q_KOSM - (K.kac_rueckkehrzeit(2 * math.pi, 1.0) - 1)) / K.Q_fehler(), "σ",
   "Gemessen (Hintergrundstrahlung mit BAO bzw. allein)", haengt_an=_Q)

# ================================================================== Schließungsregel: gebunden
_a("03a_gebunden", (D, "**Gebunden (w → 0):**"), "Zahl")
_ZC = (("Z_C_34", "14"), ("Z_C_2", "45"), ("Z_C_1", "143"), ("Z_C_05", "309"))
for _z, _t in _ZC:
    _p(f"03a_schale_{_z.lower()}", B, "03a_gebunden", _t,
       lambda p, z=_z: K.schalen_dauer(E(z), p.a0, E("W_OBEN"), E("W_UNTEN")) / K.yr / 1e6, "Mio. Jahre",
       "fällt w vor dem Umkehrpunkt von 0,9 auf 0,1", haengt_an=(_z, "W_OBEN", "W_UNTEN"))


def _nachlauf(p, delta, z, grenze):
    """1/λ_m = N*/ω_orb in Mrd. Jahren; grenze "min" (N* unten) oder "max" (N* oben)."""
    _, lo, hi = K.uebergabezeit(E(delta), z)
    return (lo if grenze == "min" else hi) / K.Gyr


_N = ("NSTERN_MIN", "NSTERN_MAX")
_p("03a_nachlauf_heute_min", B, "03a_gebunden", "0,28", lambda p: _nachlauf(p, "DELTA_500", 0.0, "min"),
   "Mrd. Jahre", "ist 1/λ_m heute 0,28–3,5 Mrd. Jahre", haengt_an=_N + ("DELTA_500",))
_p("03a_nachlauf_heute_max", B, "03a_gebunden", "3,5", lambda p: _nachlauf(p, "DELTA_200", 0.0, "max"),
   "Mrd. Jahre", "ist 1/λ_m heute 0,28–3,5 Mrd. Jahre", haengt_an=_N + ("DELTA_200",))
_p("03a_nachlauf_z2_min", B, "03a_gebunden", "0,1", lambda p: _nachlauf(p, "DELTA_500", E("Z_C_2"), "min"),
   "Mrd. Jahre", "für eine bei z = 2 kollabierte Galaxie 0,1–1,2 Mrd. Jahre", haengt_an=_N + ("DELTA_500", "Z_C_2"))
_p("03a_nachlauf_z2_max", B, "03a_gebunden", "1,2", lambda p: _nachlauf(p, "DELTA_200", E("Z_C_2"), "max"),
   "Mrd. Jahre", "für eine bei z = 2 kollabierte Galaxie 0,1–1,2 Mrd. Jahre", haengt_an=_N + ("DELTA_200", "Z_C_2"))

# ================================================================== Schließungsregel: Galaxien, Naht, Haufen
_a("03a_galaxien_RQ", (D, "Die Bilanz Q·M_b = (ν − 1)·M_b gilt genau am Radius R_Q"), "Zahl")


def _aN_RQ(p):
    return K.G * M(MW) / K.r_Q(M(MW), p.a0) ** 2


_p("03a_aN_RQ", B, "03a_galaxien_RQ", "0,029", lambda p: _aN_RQ(p) / p.a0, "a₀", "(a_N = 0,029 a₀",
   haengt_an=(MW,) + _GM + _Q)
_p("03a_RQ_milchstrasse", B, "03a_galaxien_RQ", "51", lambda p: K.r_Q(M(MW), p.a0) / K.kpc, "kpc", "Milchstraße 51 kpc",
   haengt_an=(MW,) + _GM + _Q)
_p("03a_chi_RQ", B, "03a_galaxien_RQ", "0,917", lambda p: K.chi(_aN_RQ(p), p.a0), "1", "= 0,917",
   haengt_an=(MW,) + _GM + _Q)

_a("03a_naht", (D, "Bis zum Umkehrradius (θ̃ = 0) ist w = 0"), "Zahl")


def _R_ta(k):
    return K.umkehrradius_heute((1 + K.Q_KOSM) * M(k))


_p("03a_Rta_milchstrasse", B, "03a_naht", "0,58", lambda p: _R_ta(MW) / K.Mpc, "Mpc", "Milchstraße an R_ta = 0,58 Mpc",
   haengt_an=(MW,) + _GM + _Q)
_p("03a_nu_Rta", B, "03a_naht", "67", lambda p: K.nu(K.G * M(MW) / _R_ta(MW) ** 2, p.a0), "1", "ν ≈ 67",
   haengt_an=(MW,) + _GM + _Q)
for _k, _t in (("M_B_KLEIN", "319"), (MW, "580"), ("M_B_GROSS", "1482")):
    _p(f"03a_Rta_{_k.lower()}", B, "03a_naht", _t, lambda p, k=_k: _R_ta(k) / K.kpc, "kpc", "an R_ta ∝ M_b^{1/3}",
       haengt_an=(_k,) + _GM + _Q)
_p("03a_w_mitbewegt", B, "03a_naht", "0,985", lambda p: K.gewicht_w(3 * K.hubble_rate(p), p.a0), "1",
   "mitbewegt mit θ̃ = 3H 0,985")

_a("03a_haufen", (D, "Gebunden ist die Mode nach wenigen Radiant der Bahn übergeben"), "Zahl")
_p("03a_masse_bQ", B, "03a_haufen", "4,45·10¹³", lambda p: K.masse_bQ(p.a0, K.hubble_rate(p)) / K.M_sun, "M☉",
   "M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴)", haengt_an=_GM + _Q)
