# -*- coding: utf-8 -*-
"""
Gesetze Klasse B: Gesetze der Galaxien aus der Literatur, mit der örtlichen Regel des Kerns nachgerechnet. Aufbau wie
gesetze_A.py; Kennungen beginnen mit "B_". Messwerte nur über das Register, Eingaben über ihr Register; Rechenpunkte
reiner Grenzfälle (tiefe Plummer-Kugel, Skalierung um 3) stehen hier als Namen.
SPARC: Auswahl, Υ_b aus SPARC_SATZ, Υ_d aus dem Parametersatz und den Eingaben der unteren und oberen Lesung.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE

from .gesetz import GESETZE, hinweis, ident, mess, mess_wert, sim, wahr, Ergebnis

_g = GESETZE.neu
MW = "M_MILCHSTRASSE"
MSPC = K.M_sun / K.pc**2


def _E(k):
    return EINGABEN[k].wert


def _m_mw():
    return _E(MW) * K.M_sun


# ------------------------------------------------------------------ Form, Tiefbereich, Virialsatz, Tully-Fisher
GN_PUNKTE = np.logspace(-14, -8, 13)          # m/s², über die ganze Relation
_g("B_beschleunigungsrelation", "B", "Beschleunigungsrelation: Form g/(1 − e^{−√(g/g†)}) = örtliche Regel ν(g_N)·g_N",
   lambda p: ident(lambda: K.a_regel(GN_PUNKTE, p.a0), lambda: GN_PUNKTE / (1 - np.exp(-np.sqrt(GN_PUNKTE / p.a0))),
                   1e-12),
   "McGaugh, Lelli, Schombert 2016 (PRL 117, 201101)", "LI", "Wirbel")

TIEF_M, TIEF_R, SKALA = 1e7, 100.0, 3.0     # M☉, kpc: a_N = 10⁻⁴·a₀; Streckung um 3


def _skaleninvarianz(p):
    M = TIEF_M * K.M_sun
    r = TIEF_R * K.kpc
    return K.a_regel(K.G * M / (SKALA * r)**2, p.a0) * SKALA / K.a_regel(K.G * M / r**2, p.a0)


_g("B_tiefbereich_skaleninvariant", "B", "Tiefbereich skaleninvariant: a(λr)·λ = a(r) (a_N = 10⁻⁴·a₀, λ = 3)",
   lambda p: sim(lambda: 1.0, lambda: _skaleninvarianz(p), 0.02), "Milgrom 2009 (ApJ 698, 1630)", "LI", "Wirbel")

VIRIAL_M, VIRIAL_B = 1e5, 20.0               # M☉, kpc: Plummer-Kugel ganz im Tiefbereich


def _virial_sigma2(p):
    """Massengewichtetes σ² (eindimensional) der Plummer-Kugel aus dem Virialsatz mit der örtlichen Regel."""
    Mv, b = VIRIAL_M * K.M_sun, VIRIAL_B * K.kpc

    def f(u):
        x = math.exp(u)
        rho = 3 * Mv / (4 * math.pi * b**3) * (1 + (x / b)**2)**-2.5
        menc = Mv * x**3 / (x * x + b * b)**1.5
        return 4 * math.pi * x**3 * rho * float(K.a_regel(K.G * menc / x**2, p.a0)) * x
    return quad(f, math.log(b) - 20, math.log(b) + 20, limit=400)[0] / (3 * Mv)


_g("B_virialsatz_tief", "B", "Virialsatz im Tiefbereich σ⁴ = (4/81)GMa₀ (Plummer, örtliche Regel, jedes Profil)",
   lambda p: sim(lambda: K.sigma_v_tief(VIRIAL_M * K.M_sun, p.a0)**2, lambda: _virial_sigma2(p), 0.02),
   "Milgrom 1995; McGaugh, Milgrom 2013", "LI", "Wirbel")

WEIT_AUSSEN = 1e3                            # r = 10³·r*: dort trägt der Rest GM/(2r) nur 2,5·10⁻⁴ in v


def _flache_kurve(p):
    M = _m_mw()
    r = WEIT_AUSSEN * K.r_stern(M, p.a0)
    return math.sqrt(K.a_regel(K.G * M / r**2, p.a0) * r)


_g("B_tully_fisher_grenzfall", "B", "Tully-Fisher v⁴ = GMa₀ als Grenzfall der örtlichen Regel weit außen",
   lambda p: sim(lambda: K.v_flach(_m_mw(), p.a0), lambda: _flache_kurve(p), 1e-3), "Milgrom 1983; McGaugh 2016",
   "LI", "Wirbel", haengt_an=(MW,))

BTF_MASSEN = (1e9, 1e11)                     # M☉: zwei Massen weit außen (10³·r*)


def _btf_steigung(p):
    """d log M/d log v der äußeren Geschwindigkeit, gerechnet mit der örtlichen Regel (nicht gesetzt)."""
    lv = []
    for m in BTF_MASSEN:
        M = m * K.M_sun
        r = WEIT_AUSSEN * K.r_stern(M, p.a0)
        lv.append(math.log(math.sqrt(K.a_regel(K.G * M / r**2, p.a0) * r)))
    return math.log(BTF_MASSEN[1] / BTF_MASSEN[0]) / (lv[1] - lv[0])


_g("B_btf_steigung", "B", "Baryonische Tully-Fisher: Steigung log M gegen log v aus der Regel (Grenzfall 4) gegen die Messung",
   lambda p: mess(lambda: _btf_steigung(p), "LELLI_BTF_STEIGUNG"), "Lelli et al. 2019 (MNRAS 484, 3267)", "LI",
   "Wirbel")


# ------------------------------------------------------------------ Zentrale Flächendichte (Toomre)
EXP_M, EXP_RD, EXP_PUNKTE = 1.0, 1.0, 20001   # Newton-Kurve der Exponentialscheibe, 10⁻⁵ … 10³ R_d logarithmisch


def _toomre_newton():
    R = np.geomspace(1e-5, 1e3, EXP_PUNKTE) * EXP_RD
    V = np.sqrt(R * K.g_expscheibe(EXP_M, EXP_RD, R))
    return K.toomre_sigma_dyn(R, V)


_g("B_toomre_exponentialscheibe", "B",
   "Toomre: (1/2πG)∫V²/R² dR = Σ(0) für die Exponentialscheibe (Newton, abgetastete Kurve)",
   lambda p: ident(_toomre_newton, lambda: EXP_M / (2 * math.pi * EXP_RD**2), 1e-4), "Toomre 1963, Freeman 1970",
   "L", "Zentrum")

SIGMA0_FIT = np.logspace(0.5, 4.5, 25)       # M☉/pc², über den Bereich der Messung


def _lelli(S):
    SM, SC, al, be = MESSWERTE["LELLI_SIGMA_DYN"].wert
    return SM * (1 + S / SC)**al * (S / SC)**be


def _beta_regel(p):
    """β der Doppelpotenz Σ_dyn = Σ_M(1 + Σ₀/Σ_C)^{1−β}(Σ₀/Σ_C)^β, angepasst an die ideale Scheibe mit der Regel."""
    sd = np.array([K.sigma_dyn_expscheibe(s * MSPC, p.a0) / MSPC for s in SIGMA0_FIT])

    def dpl(x, s0, sc, bet):
        return np.log10(s0 * (1 + 10**x / sc)**(1 - bet) * (10**x / sc)**bet)
    return curve_fit(dpl, np.log10(SIGMA0_FIT), np.log10(sd), p0=(1500, 1500, 0.6), maxfev=20000)[0][2]


_g("B_flaechendichte_beta", "B", "Zentrale Flächendichte, ideale Exponentialscheibe: Steigung β der Doppelpotenz (α = 1 − β)",
   lambda p: mess_wert(lambda: _beta_regel(p), MESSWERTE["LELLI_SIGMA_DYN"].wert[3],
                       MESSWERTE["LELLI_BETA_FEHLER"].wert, ("LELLI_SIGMA_DYN", "LELLI_BETA_FEHLER")),
   "Lelli et al. 2016 (ApJL 827, L19)", "L", "Zentrum")

SIGMA0_PUNKTE = (EINGABEN["SIGMA0_EXPSCHEIBE"].wert[0], 100.0, EINGABEN["SIGMA0_EXPSCHEIBE"].wert[1])   # M☉/pc²


def _lelli_abstand(p):
    streu = MESSWERTE["LELLI_STREUUNG"].wert
    dd = [math.log10(K.sigma_dyn_expscheibe(s * MSPC, p.a0) / MSPC / _lelli(s)) for s in SIGMA0_PUNKTE]
    x = max(dd, key=abs)
    return Ergebnis("OK" if abs(x) <= streu else "SPANNUNG", x, 0.0, " / ".join(f"{d:+.2f}" for d in dd),
                    ("LELLI_SIGMA_DYN", "LELLI_STREUUNG"), None, None, "mess")


_g("B_flaechendichte_lelli_fit", "B",
   "Zentrale Flächendichte, ideale Exponentialscheibe gegen den Lelli-Fit bei Σ₀ = 10 / 100 / 1000 M☉/pc² [dex], "
   "innerhalb der Streuung", _lelli_abstand, "Lelli et al. 2016 (ApJL 827, L19)", "L", "Zentrum",
   haengt_an=("SIGMA0_EXPSCHEIBE",))


def _sparc_toomre(ups_d, a0):
    """Je Galaxie (Σ₀ = Υ_d·SB₀, Σ_dyn gemessen, Σ_dyn der Regel mit den eigenen Baryonen), in M☉/pc²; radiale Form."""
    S = K.SPARC_SATZ
    aus = []
    for g in K.sparc_tafel(S["q_max"], S["inc_min"]):
        vb2 = K.v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], ups_d, S["ups_b"])
        gN = np.clip(vb2, 1e-6, None) / g["R"]
        vm = np.sqrt(K.nu(gN, a0) * gN * g["R"])
        aus.append((ups_d * g["SB0"], K.toomre_sigma_dyn(g["R"], g["Vo"]) / MSPC, K.toomre_sigma_dyn(g["R"], vm) / MSPC))
    a = np.array(aus)
    return a[(a > 0).all(axis=1)]


def _rest(ups_d, a0):
    a = _sparc_toomre(ups_d, a0)
    y = np.log10(a[:, 1] / a[:, 2])
    return a, y


def _sparc_rest(p):
    """Mittel des Rests log Σ_gem/Σ_Regel; Fehler aus Streuung und der Spanne Υ_d = 0,4–0,6 (halbe Differenz)."""
    _, y = _rest(p.Upsilon_d, p.a0)
    spanne = abs(_rest(_E("UPS_D_06"), p.a0)[1].mean() - _rest(_E("UPS_D_04"), p.a0)[1].mean()) / 2
    f = math.hypot(y.std(ddof=1) / math.sqrt(len(y)), spanne)
    n = y.mean() / f
    return Ergebnis("OK" if abs(n) <= 2 else ("SPANNUNG" if abs(n) <= 3 else "WIDERSPRUCH"), float(y.mean()), 0.0,
                    f"{n:+.1f}σ", (), None, None, "mess")


_SP = ("SPARC_Q_MAX", "SPARC_INC_MIN", "SPARC_UPS_B", "UPS_D_04", "UPS_D_06")
_g("B_flaechendichte_sparc_rest", "B",
   "Zentrale Flächendichte Galaxie für Galaxie (eigene Baryonen): Rest log Σ_gem/Σ_Regel [dex] null",
   _sparc_rest, "SPARC, Lelli et al. 2016 (Tafel und Kurven)", "L", "Zentrum", haengt_an=_SP)


def _sparc_trend(p):
    a, y = _rest(p.Upsilon_d, p.a0)
    x = np.log10(a[:, 0])
    A = np.vstack([x - x.mean(), np.ones_like(x)]).T
    cf, rs, *_ = np.linalg.lstsq(A, y, rcond=None)
    ek = math.sqrt(rs[0] / (len(x) - 2) / np.sum((x - x.mean())**2))
    n = cf[0] / ek
    return Ergebnis("OK" if abs(n) <= 2 else ("SPANNUNG" if abs(n) <= 3 else "WIDERSPRUCH"), float(cf[0]), 0.0,
                    f"{n:+.1f}σ", (), None, None, "mess")


_g("B_flaechendichte_sparc_trend", "B", "Zentrale Flächendichte Galaxie für Galaxie: Trend des Rests je Dekade Σ₀ null",
   _sparc_trend, "SPARC, Lelli et al. 2016 (Tafel und Kurven)", "L", "Zentrum", haengt_an=_SP)


# ------------------------------------------------------------------ Druckgestützte Systeme und äußeres Feld
def _crater2(p):
    """Crater II im Feld der Milchstraße: Plummer-Virialsumme (σ² = (π/32)·νGM/R_e) gegen den Schätzer von McGaugh
    (σ² = νGM/(4R_e), ganze statt halber Masse) mit derselben Verstärkung ν(g_e)."""
    M = MESSWERTE["CRATER2_M"].wert * K.M_sun
    re = MESSWERTE["CRATER2_RE"].wert * K.kpc
    nue = float(K.nu(K.G * _m_mw() / (MESSWERTE["CRATER2_ABSTAND"].wert * K.kpc)**2, p.a0))
    s_mcg = math.sqrt(nue * K.G * M / (3 * 4 / 3 * re)) / 1e3
    s_vir = math.sqrt(math.pi / 32 * nue * K.G * M / re) / 1e3
    return hinweis(lambda: s_vir, lambda: s_mcg, f"Faktor {s_mcg**2 / s_vir**2:.2f} in σ²",
                   ("CRATER2_M", "CRATER2_RE", "CRATER2_ABSTAND"))


_g("B_crater2_schaetzer", "B",
   "Crater II im Außenfeld: Plummer-Virial gegen den Schätzer McGaugh (σ² ∝ M statt M/2) [km/s]", _crater2,
   "McGaugh 2016 (ApJL 832, L8): 2,1 (+0,9/−0,6)", "L", "Hintergrund", haengt_an=(MW,))


def _aussenfeld(p):
    gi = K.G * _m_mw() / (_E("R_50KPC") * K.kpc)**2
    ge = _E("GE_003") * p.a0
    frei, senk, gleich = float(K.nu(gi, p.a0)), K.nu_aussenfeld(gi, ge, "senkrecht", p.a0), K.nu_aussenfeld(gi, ge, "gleich", p.a0)
    return wahr(frei > senk > gleich, f"ν {frei:.2f} > {senk:.2f} > {gleich:.2f}")


_g("B_aussenfeld_richtung", "B",
   "Äußeres Feld senkt ν, die Richtung zählt: frei > senkrecht > gleichgerichtet (Milchstraße, 50 kpc, 0,03·a₀)",
   _aussenfeld, "Chae et al. 2020 (ApJ 904, 51)", "L", "Hintergrund", haengt_an=(MW, "R_50KPC", "GE_003"))
