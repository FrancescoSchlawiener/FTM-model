# -*- coding: utf-8 -*-
"""
Gesetze Klasse M: unsere Gesetze aus den Bausteinen 01–09 und den Notizen 71/72, zweiter Rechenweg (Nullstelle,
Metrik, Differenzenrechnung, Reihe, Monte Carlo). Aufbau wie gesetze_A.py; Kennungen beginnen mit "M_".
Dazu die Brücken der Ebenen (Daten in kern/zustandsraum.BRUECKEN), die Budgets je Ebene, die Grenzfälle des
Sonnensystems und die Beziehungen, die als alte Prüfungen des Prüfstands kamen (Beziehung → Gesetz).
Rechenpunkte (Masse, Radius, s) ohne Eingabe sind Stellen, an denen die Beziehung gilt; sie gilt für jeden Wert.
"""
import functools
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import solve_banded
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gamma, gammaln, zeta

import kern as K
from forschung.rueckkehr import gewicht_w_maxwell
from kern.eingabe import EINGABEN
from kern.zustandsraum import (BRUECKEN, GROESSEN_ROLLEN, GRENZLAENGEN, GRENZLAENGEN_ORT, ROLLEN_ZUSTANDSRAUM,
                               TRIADE_GLIEDER)
from kern.symbole import SYMBOLE

from .aussage import AUSSAGEN
from .gesetz import GESETZE, ident, sim, wahr
from .modelltext import datei
from .sonnenort import SONNENORT, gezeit

_g = GESETZE.neu
SEED = 29092026
MW = "M_MILCHSTRASSE"
M = K.M_sun
GM = K.GM_sun


def _E(k):
    return EINGABEN[k].wert


def _m_mw():
    return _E(MW) * K.M_sun


def _lam(m):
    return K.hbar / (m * K.c)


# ================================================================== M1 Formeln A1–A3
_g("M_neigung_A1", "M", "(A1) g = x/y²·g⊕ gegen GM/R² (Sonne; g⊕ = 9,81 gesetzt)",
   lambda p: ident(lambda: K.neigung_xy(K.M_sun / K.M_earth, K.R_sun / K.R_earth), lambda: GM / K.R_sun**2, 2e-3),
   "Baustein 01", "L", "Wirbel", haengt_an=("G_ERDE", "R_sun", "R_earth"))
_g("M_zeitdifferential_bahn", "M", "Zeitdifferential einer Bahn δ = 3Δs/r = v²/c² (Erdbahn, v aus dem Jahr)",
   lambda p: ident(lambda: K.delta(M, K.AU), lambda: (2 * math.pi * K.AU / K.yr)**2 / K.c**2, 2e-3), "Baustein 01",
   "L", "Zentrum")
_g("M_brechzahl_strahl", "M", "Brechzahl n = 1 + 6Δs/r als Strahlintegral: α = 12Δs/b (Takt und Faltung je die Hälfte)",
   lambda p: ident(lambda: K.ablenkung(M, K.R_sun),
                   lambda: 6 * K.ds(M) / K.R_sun * quad(lambda u: (1 + u * u)**-1.5, -np.inf, np.inf)[0], 1e-8),
   "Baustein 08", "L", "Wirbel")


# ================================================================== M2 Produkt fest, Übergang bei Gleichheit (Notiz 72)
def _wurzel_log(f, lo, hi):
    return math.exp(brentq(lambda u: math.log(f(math.exp(u))), lo, hi))


def _m_pl_wurzel():
    return _wurzel_log(lambda m: 3 * K.ds(m) / _lam(m), -40, 0)


def _m_c_wurzel(a0):
    return _wurzel_log(lambda m: K.r_stern(m, a0) / _lam(m), -200, 0)


def _m_u_wurzel(a0):
    return _wurzel_log(lambda m: K.r_stern(m, a0) / (3 * K.ds(m)), 100, 200)


def _m_min_wurzel(a0):
    """λ̄ = ℓ₀: Nullstelle in log₁₀ M."""
    return 10**brentq(lambda lm: math.log(K.hbar / (10**lm * K.c) / (K.c**2 / a0)), -80, -60)


_g("M_schnitt_planck", "M", "Schnitt 3Δs = λ̄ bei m_Pl (Nullstelle)", lambda p: ident(_m_pl_wurzel, lambda: K.m_Pl, 1e-9),
   "Bausteine 02, 08; Notiz 72", "L", "Zentrum")
_g("M_schnitt_compton", "M", "Schnitt r* = λ̄ bei M_c = (ħ²a₀/(Gc²))^{1/3} (Nullstelle)",
   lambda p: ident(lambda: K.masse_r_stern_gleich_compton(p.a0), lambda: _m_c_wurzel(p.a0), 1e-9),
   "Baustein 08; Notiz 72", "LI", "Zentrum")
_g("M_schnitt_oben", "M", "Schnitt r* = 3Δs bei M_U = c⁴/(Ga₀) = c³τ_D/G = m_Pl/ε_L (Nullstelle)",
   lambda p: ident(lambda: K.schnittmasse_oben(p.a0), lambda: _m_u_wurzel(p.a0), 1e-9), "Baustein 08; Notiz 72", "LI",
   "Zentrum")
_g("M_schnitt_unten", "M", "Schnitt λ̄ = ℓ₀ bei M_min = ħ/(c²τ_D) = ε_L·m_Pl (Nullstelle)",
   lambda p: ident(lambda: K.schnittmasse_unten(p.a0), lambda: _m_min_wurzel(p.a0), 1e-9), "Baustein 08", "LI",
   "Zentrum")
_g("M_laengenleiter", "M", "Längenleiter M_c³·M_U = m_Pl⁴ (aus den Nullstellen)",
   lambda p: ident(lambda: _m_c_wurzel(p.a0)**3 * _m_u_wurzel(p.a0), lambda: _m_pl_wurzel()**4, 1e-8), "Baustein 08",
   "LI", "Zentrum")
PROBE_MASSEN = np.logspace(-20, 42, 9)            # kg, von unter der Planck-Masse bis über die Zelle
_g("M_produkt_fest", "M", "r*²·λ̄ = ℓ_P²·ℓ₀ für jede Masse (Stichprobe 10⁻²⁰ … 10⁴² kg)",
   lambda p: ident(lambda: np.array([K.r_stern(m, p.a0)**2 * K.compton_red(m) for m in PROBE_MASSEN]),
                   lambda: K.l_P**2 * K.c**2 / p.a0 * np.ones(PROBE_MASSEN.size), 1e-12), "Baustein 08; Notiz 72",
   "LI", "Zentrum")
_g("M_produkt_faltung_compton", "M", "Produkt fest: 3Δs·λ̄ = ℓ_P² für jede Masse (Δs·λ̄/ℓ_P² = ⅓)",
   lambda p: ident(lambda: np.array([K.ds(m) * K.compton_red(m) for m in PROBE_MASSEN]),
                   lambda: K.hbar * K.G / K.c**3 / 3 * np.ones(PROBE_MASSEN.size), 1e-12), "Baustein 02", "L", "Zentrum")
_g("M_compton_mittel", "M", "Geometrisches Mittel λ̄_e² = r_e·a_B",
   lambda p: ident(lambda: K.compton_red(K.m_e)**2, lambda: K.r_e * K.a_B, 1e-9), "Baustein 08", "L", "Zentrum")


def _sternskala():
    D = brentq(lambda x: x * x * K.ds(K.m_p) - (_lam(K.m_p) / 3)**3, 1, 1e4)
    return 3 * K.c**2 * D / K.G


_g("M_sternskala", "M", "Sternskala: Δs_*²·Δs_p = (λ̄_p/3)³ ⇒ M_* = m_Pl³/m_p²",
   lambda p: ident(lambda: K.sternmasse(), _sternskala, 1e-9), "Baustein 02", "L", "Zentrum")

# ================================================================== M3 Verwindung (Kerr aus der Metrik)
KERR_A, KERR_TH = 0.6, 0.7                       # a_J in GM/c², Polwinkel


def _nullflaeche_metrik():
    gtt = lambda r: -(1 - 2 * r / (r * r + KERR_A**2 * math.cos(KERR_TH)**2))           # noqa: E731
    return brentq(gtt, 1.0001, 5) * K.G * M / K.c**2


_g("M_nullflaeche_takt", "M", "Nullfläche des Takts r_E(θ) = 3Δs + √(9Δs² − a_J²cos²θ) (Nullstelle von g_tt, θ = 0,7)",
   lambda p: ident(lambda: K.nullflaeche_takt(M, KERR_A * K.G * M / K.c**2, KERR_TH), _nullflaeche_metrik, 1e-10),
   "Baustein 08", "L", "Zentrum")


def _umlauf_metrik(r):
    """T_mit − T_gegen aus den Kreisbahnen der Metrik (Äquator, Boyer-Lindquist), in s (GM/c³ als Zeiteinheit)."""
    tt, tp, pp = -2 / r**2, 2 * KERR_A / r**2, 2 * r - 2 * KERR_A**2 / r**2
    om = [(-tp + s * math.sqrt(tp * tp - tt * pp)) / pp for s in (1, -1)]
    return (2 * math.pi / om[0] - 2 * math.pi / abs(om[1])) * K.G * M / K.c**3


for _r in (10.0, 1000.0):
    _g(f"M_umlaufdifferenz_{int(_r)}", "M", f"T_mitläufig − T_gegenläufig = 4πa_J/c unabhängig von r (Metrik, r = {_r:g} GM/c²)",
       lambda p, r=_r: ident(lambda: K.umlauf_differenz(KERR_A * K.G * M / K.c**2), lambda: _umlauf_metrik(r), 1e-6),
       "Baustein 08", "L", "Wirbel")
_g("M_spin_halb", "M", "Spin ½: a_J·Δs = ℓ_P²/6 (Elektron); j = ½(m_Pl/M)²",
   lambda p: ident(lambda: K.hbar / 2 / (K.m_e * K.c) * K.ds(K.m_e), lambda: K.l_P**2 / 6, 1e-9), "Baustein 09", "L",
   "Wirbel")
_g("M_magnetisches_moment", "M", "Magnetisches Moment μ = Q_el·c·a_J (g = 2) = μ_B = eħ/(2m_e) beim Elektron",
   lambda p: ident(lambda: K.moment_magnetisch(K.e_el, K.verwindung(K.hbar / 2, K.m_e)),
                   lambda: K.e_el * K.hbar / (2 * K.m_e), 1e-9), "Baustein 09", "L", "Wirbel")


# ================================================================== M4 Rückkehr offen und bei Λ geschlossen (07, 71)
def _bvp(L, versiegelt=True, h=1e-3):
    n = int(round(L / h))
    ab = np.zeros((3, n))
    ab[0, 1:] = 1.0
    ab[1] = -(2 + h * h)
    ab[2, :-1] = 1.0
    if versiegelt:
        ab[2, -2] = 2.0
    rhs = np.zeros(n)
    rhs[0] = -1.0
    return solve_banded((1, 1), ab, rhs)      # Eintrag i ist G(x = (i + 1)h)


_g("M_rueckkehr_geschlossen", "M", "Rückkehr bei Λ geschlossen: G(x) = cosh((L − x)/λ)/cosh(L/λ) (L = 3λ, x = 1λ)",
   lambda p: ident(lambda: K.rueckkehr_geschlossen(1.0, 3.0), lambda: _bvp(3.0)[999], 1e-5), "Baustein 07; Notiz 71",
   "I", "Wirbel")
_g("M_rueckkehr_offen", "M", "Rückkehr offen: q = e^{−s} (Differenzenrechnung, x = 1λ)",
   lambda p: ident(lambda: K.q_(p.a0, p.a0), lambda: _bvp(40.0, False)[999], 1e-5), "Baustein 07", "I", "Wirbel")


# ================================================================== M5 Schalter χ
def _elastizitaet(aN, a0, nach):
    h = 1e-5
    if nach == "a0":
        return (math.log(K.a_regel(aN, a0 * (1 + h))) - math.log(K.a_regel(aN, a0 * (1 - h)))) / (math.log1p(h) - math.log1p(-h))
    return (math.log(K.a_regel(aN * (1 + h), a0)) - math.log(K.a_regel(aN * (1 - h), a0))) / (math.log1p(h) - math.log1p(-h))


for _s in (0.3, 1.0, 3.0):
    _n = str(_s).replace(".", "_")
    _g(f"M_chi_elastizitaet_a0_{_n}", "M", f"d ln a/d ln a₀ = χ/2 (s = {_s:g})",
       lambda p, s=_s: ident(lambda: float(K.chi(s * s * p.a0, p.a0)) / 2, lambda: _elastizitaet(s * s * p.a0, p.a0, "a0"),
                             1e-7), "Bausteine 03, 07", "LI", "Wirbel")
    _g(f"M_chi_flussdimension_{_n}", "M", f"d ln a/d ln a_N = 1 − χ/2, Flussdimension 3 − χ (s = {_s:g})",
       lambda p, s=_s: ident(lambda: 1 - float(K.chi(s * s * p.a0, p.a0)) / 2,
                             lambda: _elastizitaet(s * s * p.a0, p.a0, "aN"), 1e-7), "Baustein 03; Notiz 72", "LI",
       "Wirbel")
_g("M_spektraldimension_tief", "M", "Spektraldimension d_s = 2·d ln a/d ln a_N → 1 im Tiefbereich (s = 10⁻³)",
   lambda p: ident(lambda: K.spektraldimension_regel(1e-6 * p.a0, p.a0), lambda: 2 * _elastizitaet(1e-6 * p.a0, p.a0, "aN"),
                   1e-3), "Notiz 71 (Pólya)", "LI", "Wirbel")
S_SONNE = 0.7


def _zusatz_sonne(a0):
    r = K.r_stern(M, a0) / S_SONNE
    return float(K.a_regel(GM / r**2, a0)) - GM / r**2


_g("M_zusatz_kern_nome", "M", "Zusatz der Rückkehr a − a_N = χ·(Δs/r)(3c²/r*) (Sonne, s = 0,7)",
   lambda p: ident(lambda: _zusatz_sonne(p.a0),
                   lambda: float(K.chi(GM / (K.r_stern(M, p.a0) / S_SONNE)**2, p.a0)) * K.ds(M)
                   / (K.r_stern(M, p.a0) / S_SONNE) * 3 * K.c**2 / K.r_stern(M, p.a0), 1e-12), "Baustein 07", "LI",
   "Wirbel")


# ================================================================== M6 Zusatzpotential, M7 Epizyklen, M8 Entropie
def _zusatz_integral(M_, r, a0, verschiebung=0.0):
    f = lambda x: float(K.a_regel(K.G * M_ / (x + verschiebung)**2, a0)) - K.G * M_ / (x + verschiebung)**2   # noqa: E731
    return quad(f, 1e-9 * r if verschiebung == 0 else 0.0, r, limit=400)[0]


for _s in (0.5, 2.0):
    _g(f"M_zusatzpotential_{str(_s).replace('.', '_')}", "M",
       f"Φ_X = v_f²·ln ν(s) als Integral von a − a_N (Punktmasse, s = {_s:g})",
       lambda p, s=_s: ident(lambda: K.zusatzpotential(_m_mw(), K.r_stern(_m_mw(), p.a0) / s, p.a0),
                             lambda: _zusatz_integral(_m_mw(), K.r_stern(_m_mw(), p.a0) / s, p.a0), 1e-7),
       "Baustein 08", "LI", "Wirbel", haengt_an=(MW,))
HERNQ_B, HERNQ_R = 0.4, 1.5                      # b_H und r in r*


def _hernquist_kern(a0):
    rs = K.r_stern(_m_mw(), a0)
    return K.zusatzpotential(_m_mw(), (HERNQ_R + HERNQ_B) * rs, a0) - K.zusatzpotential(_m_mw(), HERNQ_B * rs, a0)


_g("M_zusatzpotential_hernquist", "M", "Φ_X Hernquist = v_f²·ln[ν(s)/ν(s₀)] (b_H = 0,4 r*, r = 1,5 r*)",
   lambda p: ident(lambda: _hernquist_kern(p.a0),
                   lambda: _zusatz_integral(_m_mw(), HERNQ_R * K.r_stern(_m_mw(), p.a0), p.a0,
                                            HERNQ_B * K.r_stern(_m_mw(), p.a0)), 1e-7),
   "Baustein 08", "LI", "Wirbel", haengt_an=(MW,))


def _kappa2_om2(aNf, r):
    """κ²/Ω² aus der Kreisbahn der Regel, numerisch (Einheiten GM = a₀ = r* = 1)."""
    O2 = lambda x: float(K.a_regel(aNf(x), 1.0)) / x                        # noqa: E731
    return (r * (O2(r * 1.0001) - O2(r / 1.0001)) / (r * 1.0001 - r / 1.0001) + 4 * O2(r)) / O2(r)


for _r in (0.5, 2.0):
    _g(f"M_re_tau_punkt_{str(_r).replace('.', '_')}", "M", f"(Re τ_M)² = κ²/Ω² = 1 + χ, Punktmasse (r = {_r:g} r*)",
       lambda p, r=_r: ident(lambda: K.realteil_quadrat_punkt(1 / r), lambda: _kappa2_om2(lambda x: 1 / x**2, r), 1e-6),
       "Baustein 07; Notiz 72", "Z", "tangential")
_g("M_re_tau_hernquist", "M", "(Re τ_M)² = 3 + γ_N(1 − χ/2), Hernquist (r = 1 r*)",
   lambda p: ident(lambda: 3 + (-2 / (1 + HERNQ_B)) * (1 - float(K.chi(1 / (1 + HERNQ_B)**2, 1.0)) / 2),
                   lambda: _kappa2_om2(lambda x: 1 / (x + HERNQ_B)**2, 1.0), 1e-6), "Baustein 07", "Z", "tangential")


def _entropie_summe(s):
    q = math.exp(-s)
    n = np.arange(int(600 / s))
    P = (1 - q) * q**n
    return -np.sum(P * np.log(P))


for _s in (0.3, 1.0):
    _g(f"M_entropie_ruecklauf_{str(_s).replace('.', '_')}", "M", f"Entropie der Rücklaufzahl S_n = ln ν + χ (s = {_s:g})",
       lambda p, s=_s: ident(lambda: K.entropie_ruecklauf(s), lambda: _entropie_summe(s), 1e-9), "Notiz 71 (Jaynes)",
       "Z", "transversal")
_g("M_fisher_metrik", "M", "Fisher-Metrik der Rücklaufzahl g_F = −dν/ds = ν(ν − 1) = 1/(4 sinh²(s/2)) (s = 1)",
   lambda p: ident(lambda: K.fisher_metrik(1.0), lambda: 1 / (4 * math.sinh(0.5)**2), 1e-9), "Baustein 07", "Z",
   "transversal")
BERNOULLI_X = 3.0


def _bernoulli_reihe(x):
    return x + 0.5 + 1 / (12 * x) - (1 / x)**3 / 720


_g("M_bernoulli_reihe", "M", "M_schein/M = 1/(1 − e^{−r*/r}) = r/r* + ½ + r*/(12r) − (r*/r)³/720 + … (r = 3r*)",
   lambda p: ident(lambda: float(K.nu(p.a0 / BERNOULLI_X**2, p.a0)), lambda: _bernoulli_reihe(BERNOULLI_X), 1e-5),
   "Baustein 07", "LI", "Wirbel")


# ================================================================== M9 Kritische Verzweigung (71)
def _lPB(n):
    return -n + (n - 1) * np.log(n) - gammaln(n + 1)          # Borel: Gesamtzahl S der kritischen Poisson-Verzweigung


def _tfix(Gm):
    return brentq(lambda T: T - math.exp(-Gm) * math.exp(T - 1), 0.5, 1 - 1e-13, xtol=1e-15)


def _borel_summe():
    n = np.arange(1, 400001, dtype=float)
    return 1 - np.sum(np.exp(_lPB(n) - 1e-4 * n))


_g("M_verzweigung_borel", "M", "Kritische Verzweigung (Borel): Summe über P(S = n) gegen T = z·e^{T−1} (Γ = 10⁻⁴)",
   lambda p: ident(_borel_summe, lambda: 1 - _tfix(1e-4), 1e-6), "Notiz 71", "I", "Wirbel")
_g("M_verzweigung_wurzel", "M", "Kritische Verzweigung: 1 − E[e^{−ΓS}] = √(2Γ) (Γ = 10⁻⁸), der Exponent ½ der Rückkehr",
   lambda p: sim(lambda: 1 - float(K.rueckkehr(0.0, 1e-8, 0.0, K.c / 2)), lambda: 1 - _tfix(1e-8), 1e-3),
   "Notiz 71; Baustein 07", "I", "Wirbel")
_g("M_verzweigung_potenz", "M", "Kritische Verzweigung: P(S = n)·n^{3/2}·√(2π) → 1 (n = 10⁷)",
   lambda p: ident(lambda: math.exp(_lPB(1e7)) * 1e7**1.5 * math.sqrt(2 * math.pi), lambda: 1.0, 2e-6), "Notiz 71", "I",
   "Wirbel")


# ================================================================== M10 Dimension als Faktor
def _sehnen():
    rng = np.random.default_rng(SEED)
    U = rng.uniform(0.0, 1.0, 400000)
    return np.mean(2 * np.sqrt(1 - (2 * U - 1)**2)), np.mean(2 * np.sqrt(U))


_g("M_sehne_kreis", "M", "Mittlere Sehne im Kreis πA/P (η₂ = π; Monte Carlo)",
   lambda p: sim(lambda: math.pi / 2, lambda: _sehnen()[0], 5e-3), "Notizen 71, 72 (Cauchy)", "L", "Zentrum")
_g("M_sehne_kugel", "M", "Mittlere Sehne in der Kugel 4V/S (η₃ = 4; Monte Carlo)",
   lambda p: sim(lambda: 4 / 3, lambda: _sehnen()[1], 5e-3), "Notizen 71, 72 (Cauchy)", "L", "Zentrum")
_g("M_dimension_drei", "M", "V/(S·R) = 1/d für die d-Kugel ⇒ Δs = GM/(3c²) für d = 3",
   lambda p: ident(lambda: K.ds(M) * K.c**2 / GM, lambda: quad(lambda r: r * r, 0, 1)[0], 1e-12), "Notiz 72", "L",
   "Zentrum")


# ================================================================== M11 Modulformen (05)
def _sig(k, n):
    return sum(d**k for d in range(1, n + 1) if n % d == 0)


NS = 120
_S1 = [_sig(1, n) for n in range(NS)]


def _modulformen():
    S3 = [_sig(3, n) for n in range(NS)]
    S5 = [_sig(5, n) for n in range(NS)]
    Pf = lambda q: 1 - 24 * sum(_S1[n] * q**n for n in range(1, NS))       # noqa: E731
    Qf = lambda q: 1 + 240 * sum(S3[n] * q**n for n in range(1, NS))       # noqa: E731
    Rf = lambda q: 1 - 504 * sum(S5[n] * q**n for n in range(1, NS))       # noqa: E731
    qd = lambda f, q: q * (f(q * (1 + 1e-5)) - f(q * (1 - 1e-5))) / (2e-5 * q)   # noqa: E731
    q = 0.05
    P, Q, R = Pf(q), Qf(q), Rf(q)
    D = lambda x: x * np.prod([(1 - x**n)**24 for n in range(1, NS)])     # noqa: E731
    ok = (abs(qd(Pf, q) - (P * P - Q) / 12) < 1e-8 and abs(qd(Qf, q) - (P * Q - R) / 3) < 1e-7
          and abs(qd(Rf, q) - (P * R - Q * Q) / 2) < 1e-6 and abs((Q**3 - R**2) / 1728 / D(q) - 1) < 1e-10
          and abs(qd(D, q) / (P * D(q)) - 1) < 1e-8)
    return wahr(ok, "q = 0,05")


_g("M_modulformen", "M", "Modulformen: q dE₂/dq = (E₂² − E₄)/12 usw., Δ = (E₄³ − E₆²)/1728 = q∏(1 − qⁿ)²⁴, D₁₂Δ = 0",
   lambda p: _modulformen(), "Baustein 05; Ramanujan 1916", "Z", "tangential")


def _herkunft_zahlen():
    ok = (all(abs((2 * math.pi)**k / (gamma(k) * zeta(k)) - v) < 1e-9 for k, v in ((2, 24), (4, 240), (6, 504)))
          and 3 * 240 + 2 * 504 == 1728 and abs(2 / (1 - 1 / 2 - 1 / 3) - 12) < 1e-12
          and abs(4 * math.pi**4 / 15 - 4 * (2 * math.pi)**4 / 240) < 1e-12)
    return wahr(ok, "24, 240, 504, 1728, 12, 4π⁴/15")


_g("M_herkunft_zahlen", "M",
   "Herkunft der Zahlen: 24, 240, 504 = (2π)^k/(Γ(k)ζ(k)), 1728 = 3·240 + 2·504, 12 = 2/(1 − ½ − ⅓), 4π⁴/15 = 4(2π)⁴/240",
   lambda p: _herkunft_zahlen(), "Baustein 05", "Z", "tangential")
LAMBERT_S = 0.1
_g("M_lambert_reihe", "M", "Lambert-Reihe Σσ₁(n)e^{−ns} = π²/(6s²) − 1/(2s) + 1/24 (Modultransformation, s = 0,1)",
   lambda p: ident(lambda: sum(d * math.exp(-d * LAMBERT_S) / (1 - math.exp(-d * LAMBERT_S)) for d in range(1, 2000)),
                   lambda: math.pi**2 / (6 * LAMBERT_S**2) - 1 / (2 * LAMBERT_S) + 1 / 24, 1e-12), "Baustein 05", "Z",
   "tangential")


def _hecke():
    lq = lambda q: sum(_S1[n] * q**n for n in range(1, NS))                 # noqa: E731
    return wahr(all(abs(sum(_sig(1, p * m) * 0.3**m for m in range(1, 100)) - ((1 + p) * lq(0.3) - p * lq(0.3**p)))
                    < 1e-9 for p in (2, 3, 5)), "p = 2, 3, 5")


_g("M_hecke", "M", "Hecke: Σσ₁(pm)q^m = (1 + p)ℓ(q) − p·ℓ(q^p) (Übertrag zwischen Skalenebenen)", lambda p: _hecke(),
   "Baustein 05", "Z", "parametrisch")
RR_Q = 0.5
_g("M_rogers_ramanujan", "M", "Rogers-Ramanujan Σq^{n²}/(q;q)_n = ∏1/((1 − q^{5n+1})(1 − q^{5n+4})) (q = ½)",
   lambda p: ident(lambda: sum(RR_Q**(n * n) / np.prod([1 - RR_Q**k for k in range(1, n + 1)]) for n in range(0, 30)),
                   lambda: np.prod([1 / ((1 - RR_Q**(5 * n + 1)) * (1 - RR_Q**(5 * n + 4))) for n in range(0, 200)]),
                   1e-12), "Baustein 05", "Z", "tangential")


def _fibonacci():
    f = [1, 1]
    for _ in range(40):
        f.append(f[-1] + f[-2])
    return f


_g("M_ausschlussregel", "M", "Ausschlussregel Abstand ≥ 2: Zustände je Platz → φ (Fibonacci)",
   lambda p: ident(lambda: _fibonacci()[-1] / _fibonacci()[-2], lambda: (1 + math.sqrt(5)) / 2, 1e-12), "Baustein 05",
   "Z", "tangential")


# ================================================================== M12 Ramanujan-Summen, Resolvente 𝒫_q (05)
def _ramanujan_summen():
    N = 12
    T = np.roll(np.eye(N), 1, 0)
    cq = lambda q, n: sum(math.cos(2 * math.pi * a * n / q) for a in range(1, q + 1) if math.gcd(a, q) == 1)   # noqa: E731
    E = {q: sum(cq(q, n) * np.linalg.matrix_power(T, n) for n in range(N)) / N for q in (1, 2, 3, 4, 6, 12)}
    ok = (np.allclose(sum(E.values()), np.eye(N)) and np.allclose(E[1], np.ones((N, N)) / N)
          and all(np.allclose(E[a] @ E[b], E[a] if a == b else 0) for a in E for b in E))
    return wahr(ok, "ℤ₁₂")


_g("M_ramanujan_projektoren", "M", "Ramanujan-Summen: E_q Projektoren, ΣE_q = I, E_qE_q' = 0, E₁ = 𝒫 (ℤ₁₂)",
   lambda p: _ramanujan_summen(), "Baustein 05", "Z", "transversal")
PQ = 0.37


def _resolvente_eigen():
    ev = np.linalg.eigvals((1 - PQ) * np.linalg.inv(np.eye(64) - PQ * np.roll(np.eye(64), 1, 0)))
    return max(abs(ev - 1 / (1 + PQ)))


_g("M_resolvente_kreis", "M", "𝒫_q = (1 − q)(I − qT)⁻¹: Eigenwerte auf dem Kreis um 1/(1 + q) mit Radius q/(1 + q)",
   lambda p: ident(_resolvente_eigen, lambda: PQ / (1 + PQ), 1e-10), "Baustein 05", "Z", "transversal")

# ================================================================== M13 Wicklung
GOLD = (math.sqrt(5) - 1) / 2
_g("M_goldene_wicklung", "M", "Goldene Wicklung: n·‖nφ‖ → 1/√5 (n = F₂₀)",
   lambda p: ident(lambda: _fibonacci()[19] * abs(_fibonacci()[19] * GOLD - round(_fibonacci()[19] * GOLD)),
                   lambda: 1 / math.sqrt(5), 1e-6), "Baustein 01", "Z", "tangential")


def _drei_abstaende():
    ok, mw = True, []
    for w in (GOLD, math.sqrt(2) - 1, math.pi - 3):
        orb = (np.arange(200000) * w) % 1.0
        hit = np.where(orb < 0.137)[0]
        d = np.unique(np.diff(hit))
        ok &= len(d) <= 3 and (len(d) < 3 or d[2] == d[0] + d[1])
        mw.append(np.mean(np.diff(hit)) * 0.137)
    return wahr(ok and max(abs(np.array(mw) - 1)) < 1e-3, "Kac " + "/".join(f"{x:.4f}" for x in mw))


_g("M_drei_abstaende", "M", "Drei-Abstände-Satz: Rückkehrzeiten nehmen ≤ 3 Werte an, der größte = Summe; Mittel = Kac",
   lambda p: _drei_abstaende(), "Notiz 71; Baustein 01; Slater 1967", "Z", "tangential")


# ================================================================== M14 Kapitza, M15 Äquivalenzschicht
def _kapitza(Aw):
    om = 60.0
    A = Aw / om
    s = solve_ivp(lambda t, z: [z[1], (1.0 - A * om**2 * math.cos(om * t)) * math.sin(z[0])], (0, 30), [0.2, 0.0],
                  max_step=0.004, rtol=1e-8)
    return max(abs(s.y[0])) < 1.0


_g("M_kapitza", "M", "Kapitza-Pendel: zeitgemitteltes Einschließen oben genau für (Aω)² > 2gℓ ((Aω)²/gℓ = 3,24 / 1,21)",
   lambda p: wahr(_kapitza(1.8) and not _kapitza(1.1), "3,24 / 1,21"), "Baustein 05 (Paul-Falle); Kapitza 1951", "Z",
   "parametrisch")


def _aequivalenzschicht():
    R, rq = 1.0, 2.5
    schale = quad(lambda t: 0.5 * math.sin(t) * (rq - R * math.cos(t)) / (rq**2 + R**2 - 2 * rq * R * math.cos(t))**1.5,
                  0, math.pi)[0]
    kugel = quad(lambda r: 3 * r**2 * quad(lambda t: 0.5 * math.sin(t) * (rq - r * math.cos(t))
                                          / (rq**2 + r**2 - 2 * rq * r * math.cos(t))**1.5, 0, math.pi)[0], 0, 1)[0]
    return wahr(abs(schale * rq**2 - 1) < 1e-9 and abs(kugel * rq**2 - 1) < 1e-9,
                f"{schale * rq**2:.10f} / {kugel * rq**2:.10f}")


_g("M_aequivalenzschicht", "M", "Äquivalenzschicht: Schale und Vollkugel geben außen dasselbe Feld GM/r² (Dichte nicht eindeutig)",
   lambda p: _aequivalenzschicht(), "Baustein 05; Notiz 71", "L", "Zentrum")


# ================================================================== M16 Engstelle (04)
def _s_bh(m):
    return 4 * math.pi * (2 * K.G * m / K.c**2)**2 / (4 * K.l_P**2)          # Horizontfläche/4ℓ_P², in k


M1_, M2_ = 3 * M, 7 * M
_g("M_kreuzterm", "M", "Kreuzterm beim Verschmelzen 8πG·M₁M₂k/(ħc) (aus der Horizontfläche)",
   lambda p: ident(lambda: 4 * math.pi * K.alpha_G(M1_ + M2_) - 4 * math.pi * K.alpha_G(M1_) - 4 * math.pi * K.alpha_G(M2_),
                   lambda: _s_bh(M1_ + M2_) - _s_bh(M1_) - _s_bh(M2_), 1e-9), "Baustein 04", "L", "Zentrum")
_g("M_kreuzterm_halb", "M", "Gleiche Massen: Kreuzterm = ½ der Endentropie",
   lambda p: ident(lambda: (K.alpha_G(2 * M) - 2 * K.alpha_G(M)) / K.alpha_G(2 * M),
                   lambda: (_s_bh(2 * M) - 2 * _s_bh(M)) / _s_bh(2 * M), 1e-12), "Baustein 04", "L", "Zentrum")
_g("M_entropie_alpha_G", "M", "Entropie der Falte S = 36πk(Δs/ℓ_P)² = 4πk·α_G(M) (Sonne)",
   lambda p: ident(lambda: 4 * math.pi * K.alpha_G(M), lambda: 36 * math.pi * (K.ds(M) / K.l_P)**2, 1e-12),
   "Baustein 04", "L", "Zentrum")
RUECKGABE_M = 1e12                               # kg


def _rueckgabe_ableitung():
    tau = K.lebensdauer_falte(RUECKGABE_M)
    Mt = lambda t: RUECKGABE_M * (1 - t / tau)**(1 / 3)                    # noqa: E731
    return (Mt(1e-6 * tau) - Mt(0)) / (1e-6 * tau)


_g("M_rueckgabe_rate", "M", "Rückgabe dM/dt = −ħc⁴/(15360πG²M²) (Ableitung von M(t) aus der Lebensdauer)",
   lambda p: ident(_rueckgabe_ableitung, lambda: -K.hbar * K.c**4 / (15360 * math.pi * K.G**2 * RUECKGABE_M**2), 1e-5),
   "Baustein 04", "L", "Zentrum")
_g("M_horizont_doppelt", "M", "Horizont bei gleichem Radius doppelt so heiß: (ħH/2πk)/(ħc/(4πk·r₀)), r₀ = c/H",
   lambda p: ident(lambda: (K.hbar * K.c / (12 * math.pi * K.kB * K.faltung_horizont(K.H0_si(p.H0))))
                   / (K.hbar * K.c**3 / (8 * math.pi * K.kB * K.G * (K.c / K.H0_si(p.H0)) * K.c**2 / (2 * K.G))),
                   lambda: 2.0, 1e-12), "Baustein 04", "G", "Hintergrund")


def _nariai():
    f = lambda r, m: 1 - 2 * m / r - r * r                                # noqa: E731
    mN = brentq(lambda m: minimize_scalar(lambda r: -f(r, m), bounds=(0.01, 1), method="bounded").fun, 0.1, 0.3,
                xtol=1e-14)
    return minimize_scalar(lambda r: -f(r, mN), bounds=(0.01, 1), method="bounded", options={"xatol": 1e-12}).x**2


_g("M_nariai_drittel", "M", "Größte mögliche Falte (Nariai-Grenze, numerisch) hat ⅓ der Horizontentropie",
   lambda p: sim(lambda: 1 / 3, _nariai, 1e-6), "Baustein 04", "G", "Zentrum")


def _falten_schranke(p):
    """Eine Falte aus allen Baryonen der Hubble-Kugel (ΣΔs_b): S/S_H = Ω_b² (Faltung nur der Baryonen)."""
    H = K.H0_si(p.H0)
    M_b = K.rho_b_mitbewegt() * 4 / 3 * math.pi * (K.c / H)**3
    return 4 * math.pi * K.alpha_G(M_b) / K.horizontentropie(H)


_g("M_falten_schranke", "M", "Alle Falten zusammen höchstens Ω_b²·S_H (eine Falte aus ΣΔs_b der Hubble-Kugel)",
   lambda p: ident(lambda: _falten_schranke(p), lambda: K.OMEGA_B**2, 1e-9), "Baustein 04", "G", "Zentrum",
   haengt_an=("OMEGA_B_H2",))
S_I_PROBE = 1e90                                 # k/s, Rechenpunkt: die Ungleichung gilt für jede Rate


def _endzustand(p):
    H = K.H0_si(p.H0)
    dSdH = -2 * K.horizontentropie(H) / H          # dS_H/dH
    return S_I_PROBE / (-dSdH)


_g("M_kein_endzustand", "M", "Kein Endzustand: −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k) aus Ṡ_H ≥ Ṡ_i",
   lambda p: ident(lambda: _endzustand(p),
                   lambda: K.l_P**2 * (p.H0 * 1e3 / K.Mpc)**3 * S_I_PROBE / (2 * math.pi * K.c**2), 1e-9),
   "Baustein 04", "G", "Hintergrund")


# ================================================================== M17 Vakuumteil und Strukturwirkung (03, 07)
def _omega_aus_lambda(p):
    return K.lambda_vakuum(K.schliesslaenge(p.a0)) * K.c**2 / (3 * K.H0_si(p.H0)**2)


_g("M_vakuumteil_theta", "M", "Λ_s = π⁴/(15Λ²) ⇒ Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH)",
   lambda p: ident(lambda: _omega_aus_lambda(p),
                   lambda: math.pi**2 / 45 * (2 * math.pi * p.a0 / (K.c * p.H0 * 1e3 / K.Mpc))**2, 1e-12),
   "Baustein 03", "IG", "Hintergrund")
_g("M_entropie_vakuum", "M", "S_Λ = S_H(H_s) = (45/π³)(Λ/ℓ_P)² bei H_s² = Λ_s c²/3",
   lambda p: ident(lambda: K.entropie_vakuum(K.schliesslaenge(p.a0)),
                   lambda: 45 / math.pi**3 * (K.c**2 / (2 * p.a0) / K.l_P)**2, 1e-12), "Baustein 03", "IG", "Hintergrund")
_g("M_entropie_vakuum_heute", "M", "S_Λ = S_H/Ω_s: die Entropie des Endzustands aus der heutigen Horizontentropie",
   lambda p: ident(lambda: K.entropie_vakuum(K.schliesslaenge(p.a0)),
                   lambda: math.pi * K.c**2 / (K.l_P**2 * (p.H0 * 1e3 / K.Mpc)**2)
                   / (math.pi**2 / 45 * (2 * math.pi * p.a0 / (K.c * p.H0 * 1e3 / K.Mpc))**2), 1e-12),
   "Baustein 03", "IG", "Hintergrund")
_g("M_endzustand_H", "M", "Endzustand H_s/H₀ = √Ω_s (nur der Vakuumteil)",
   lambda p: ident(lambda: K.H_vakuum(K.schliesslaenge(p.a0)) / K.H0_si(p.H0),
                   lambda: math.sqrt(math.pi**2 / 45) * 2 * math.pi * p.a0 / (K.c * p.H0 * 1e3 / K.Mpc), 1e-12),
   "Baustein 03", "G", "Hintergrund")
_g("M_theta_kuenftig", "M", "Θ im Endzustand (H → H_s) → √45/π",
   lambda p: ident(lambda: K.theta_(p.a0, K.H_vakuum(K.schliesslaenge(p.a0))), lambda: math.sqrt(45) / math.pi, 1e-12),
   "Baustein 03", "G", "Hintergrund")


def _li(n, x):
    return sum(x**k / k**n for k in range(1, 400))


for _s in (0.5, 1.0, 3.0):
    _g(f"M_strukturwirkung_{str(_s).replace('.', '_')}", "M",
       f"Strukturwirkung: Q̃(s) + 4π⁴/15 = Q_s(s⁴), der Vakuumteil ist ihre Konstante (s = {_s:g})",
       lambda p, s=_s: ident(lambda: K.Q_tilde(s) + 4 * math.pi**4 / 15,
                             lambda: quad(lambda x: 1 / (1 - math.exp(-math.sqrt(math.sqrt(x)))) - 1, 0, s**4, limit=400)[0],
                             1e-8), "Bausteine 07, 03", "LG", "Hintergrund")


# ================================================================== M18 Kern (06)
def _flaeche(e):
    a, b = 1 + e, (1 + e)**-0.5
    return quad(lambda x: 2 * math.pi * b * math.sqrt(1 - x * x / a**2 + (b * x / a**2)**2), -a, a, limit=200,
                epsabs=1e-14, epsrel=1e-13)[0]


def _coulomb(e):
    a, b = 1 + e, (1 + e)**-0.5
    return 0.3 * quad(lambda l: 1 / ((b * b + l) * math.sqrt(a * a + l)), 0, np.inf, limit=200, epsabs=1e-14,
                      epsrel=1e-13)[0]


def _spaltbarkeit():
    k2 = lambda f, e: (f(e) + f(-e) - 2 * f(0)) / (f(0) * e * e)         # noqa: E731
    cS = (4 * k2(_flaeche, 0.01) - k2(_flaeche, 0.02)) / 3
    cC = (4 * k2(_coulomb, 0.01) - k2(_coulomb, 0.02)) / 3
    return -cS / (2 * cC)


_g("M_spaltbarkeit", "M", "Spaltbarkeit: Rückstellung ∝ (1 − x), Kippstelle bei x = E_C/(2E_S) = 1 (Ellipsoid numerisch)",
   lambda p: sim(lambda: 1.0, _spaltbarkeit, 1e-4), "Baustein 06; Bohr, Wheeler 1939", "L", "Zentrum")
Z_A_PROBE = 0.45                                 # Z/A; die Lage des Maximums hängt an keinem Koeffizienten


def _bindung_viertel():
    aS, aC = 1.0, 1.0                            # Oberfläche und Coulomb in beliebigen Einheiten
    A = minimize_scalar(lambda A: aS * A**(-1 / 3) + aC * Z_A_PROBE**2 * A**(2 / 3), bounds=(1e-3, 1e3),
                        method="bounded", options={"xatol": 1e-12}).x
    return aC * Z_A_PROBE**2 * A**(5 / 3) / (2 * aS * A**(2 / 3))


_g("M_bindung_viertel", "M", "Größte Bindung je Nukleon (Oberfläche + Coulomb) bei x = ¼, E_C = E_S/2 (für alle Koeffizienten)",
   lambda p: ident(_bindung_viertel, lambda: 0.25, 1e-6), "Baustein 06", "L", "Zentrum")


# ================================================================== M19 Aufbauzeit der Rückkehr (71)
def _aufbauzeit():
    dx, dt = 0.05, 0.001
    x = np.arange(0, 12, dx)
    cm = np.zeros_like(x)
    cm[0] = 1.0
    acc = np.zeros_like(x)
    ab = np.zeros((3, x.size - 2))
    ab[0, 1:] = 1
    ab[1] = -(2 + dx**2)
    ab[2, :-1] = 1
    rr = np.zeros(x.size - 2)
    rr[0] = -1
    cs = np.concatenate([[1.0], solve_banded((1, 1), ab, rr), [0.0]])
    cs[-1] = 1.0
    for _ in range(int(25 / dt)):
        lap = np.zeros_like(cm)
        lap[1:-1] = (cm[2:] - 2 * cm[1:-1] + cm[:-2]) / dx**2
        cm = cm + dt * (lap - cm)
        cm[0], cm[-1] = 1.0, 0.0
        acc += (1 - cm / cs) * dt
    return acc[int(3 / dx)]


_g("M_aufbauzeit", "M", "Aufbauzeit der Rückkehr τ(x) = x/(2√(Dk)) = τ_D/(2s) (x = 3λ)",
   lambda p: sim(lambda: 1.5, _aufbauzeit, 5e-3), "Notiz 71 (Morphogen)", "I", "Wirbel")


# ================================================================== M21 Eine Zahl trägt die Leiter
_g("M_leiterzahl", "M", "ε_L = ℓ_P/ℓ₀ = a₀/a_P (a_P = c²/ℓ_P)",
   lambda p: ident(lambda: K.leiterzahl(p.a0), lambda: p.a0 / (K.c**2 / K.l_P), 1e-12), "Baustein 08", "Z",
   "tangential")
_g("M_spiegel_leiter", "M", "Spiegel der Leiter: M_min·M_U = m_Pl² (aus den Nullstellen)",
   lambda p: ident(lambda: _m_min_wurzel(p.a0) * _m_u_wurzel(p.a0), lambda: K.m_Pl**2, 1e-9), "Baustein 08", "LI",
   "Zentrum")
_g("M_min_planck_laenge", "M", "Bei M_min zugleich r* = ℓ_P (Nullstelle von λ̄ = ℓ₀ eingesetzt)",
   lambda p: ident(lambda: K.r_stern(_m_min_wurzel(p.a0), p.a0), lambda: K.l_P, 1e-9), "Baustein 08", "LI", "Zentrum")
_g("M_zelle_masse", "M", "M_U = ⅔·ρ_m·L³ (Zelle aus a₀ = 3c⁴/(2Gρ_mL³))",
   lambda p: ident(lambda: 2 / 3 * K.rho_m_mitbewegt() * K.zelle_aus_a0(p.a0)**3, lambda: _m_u_wurzel(p.a0), 1e-9),
   "Bausteine 08, 03", "IG", "Zentrum", haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "H0_KMS"))
_g("M_entropie_leiter", "M", "S_max(Λ)/k = (45/π³)(Λ/ℓ_P)² = (45/4π³)·M_U/M_min",
   lambda p: ident(lambda: K.entropie_vakuum(K.schliesslaenge(p.a0)),
                   lambda: 45 / (4 * math.pi**3) * _m_u_wurzel(p.a0) / _m_min_wurzel(p.a0), 1e-9), "Bausteine 08, 03",
   "IG", "Zentrum")
_g("M_faltungstakt_MU", "M", "Faltungstakt der Zellmasse: c/(3Δs(M_U)) = a₀/c",
   lambda p: ident(lambda: K.faltungstakt(K.schnittmasse_oben(p.a0)), lambda: p.a0 / K.c, 1e-12), "Baustein 08", "IG",
   "Zentrum")
_g("M_ds_MU_lambda", "M", "Faltung der Zellmasse: Δs(M_U)/Λ = ⅔",
   lambda p: ident(lambda: K.ds(_m_u_wurzel(p.a0)) / K.schliesslaenge(p.a0), lambda: 2 / 3, 1e-9), "Baustein 08", "IG",
   "Zentrum")


# ================================================================== M22 Spiegel der Diffusion
def _diffusion_spiegel(a0):
    return max(abs((K.r_stern(m, a0)**2 / K.tau_D(a0)) * (K.hbar / m) / (K.l_P * K.c)**2 - 1)
               for m in np.logspace(-30, 42, 9))


_g("M_diffusion_spiegel", "M", "D·ħ/M = (ℓ_P c)² für jede Masse (D = r*²/τ_D aus der Rückkehr; Rechnung, nicht im Modell)",
   lambda p: ident(lambda: 1.0 + _diffusion_spiegel(p.a0), lambda: 1.0, 1e-12), "Rechnung, nicht im Modell; Nelson 1966",
   "LI", "Wirbel")
_g("M_quantendiffusion_spiegel", "M", "Quantendiffusion der Spiegelmasse ħ/(2M′) = D/2, M′ = m_Pl²/M (Sonne)",
   lambda p: ident(lambda: K.hbar / (2 * K.spiegel_masse(M)), lambda: K.r_stern(M, p.a0)**2 * p.a0 / K.c / 2, 1e-12),
   "Baustein 08", "LI", "Wirbel")
_g("M_faltungstakt_spiegel", "M", "Faltungstakt c/(3Δs) = Compton-Takt der Spiegelmasse (Proton)",
   lambda p: ident(lambda: K.faltungstakt(K.m_p), lambda: K.compton_takt(K.spiegel_masse(K.m_p)), 1e-12), "Baustein 09",
   "L", "Zentrum")
_g("M_compton_faltung", "M", "ω_C² = Γ_G·c/(3Δs) (Proton): Compton-Takt als Mittel von Faltungsrate und Faltungstakt",
   lambda p: ident(lambda: K.faltungsrate(K.m_p) * K.faltungstakt(K.m_p), lambda: (K.m_p * K.c**2 / K.hbar)**2, 1e-12),
   "Baustein 09", "L", "Zentrum")
_g("M_carter_frequenzen", "M", "Carter als Frequenzleiter: Γ_G·ω₁³/(αω₂)⁴ = 𝒞",
   lambda p: ident(lambda: K.carter_faktor(), lambda: K.faltungsrate(K.m_p) * K.compton_takt(K.m_p)**3
                   / (K.alpha * K.umlauf_wasserstoff())**4, 1e-9), "Bausteine 02, 09", "Z", "tangential")


# ================================================================== M23 Barriere (das Tunnelbild ist Forschung: forschung/tunnelbild.py)
_g("M_barriere_min", "M", "Barriere bei r* = Ruheenergie der untersten Schnittmasse: kT_I = ħa₀/c = M_min c²",
   lambda p: ident(lambda: K.temperatur_innen(p.a0), lambda: _m_min_wurzel(p.a0) * K.c**2 / K.kB, 1e-9), "Bausteine 07, 08",
   "I", "Zentrum")

# ================================================================== M24 Schluss bei Λ, M25 Vakuum in Planck-Einheiten
KAC_L = 2000.0


def _abweichung_numerisch():
    x = 2 * KAC_L / math.pi
    s = 1 / x
    q = math.exp(-s) * (1 + math.exp(-2 * s * (KAC_L - 1))) / (1 + math.exp(-2 * s * KAC_L))
    return (1 - math.exp(-s)) / (1 - q) - 1


_g("M_schluss_abweichung", "M", "Geschlossene Rückkehr: Abweichung bei r = 2Λ/π gleich 1/tanh(π/2) − 1 (Γ = kleinster Eigenwert)",
   lambda p: sim(lambda: K.abweichung_schluss(math.pi / 2), _abweichung_numerisch, 2e-3), "Baustein 07", "I", "Wirbel")
_g("M_vakuum_planck", "M", "Λ_s·ℓ_P² = (4π⁴/15)·ε_L²",
   lambda p: ident(lambda: K.lambda_vakuum(K.schliesslaenge(p.a0)) * K.l_P**2,
                   lambda: 4 * math.pi**4 / 15 * (K.l_P * p.a0 / K.c**2)**2, 1e-12), "Baustein 08", "IG", "Hintergrund")


def _einseitigkeit_linear(p):
    R = _E("R_20KPC") * K.kpc
    gi = K.G * _m_mw() / R**2
    ge = _E("GE_003") * p.a0
    s = math.sqrt(gi / p.a0)
    n = float(K.nu(gi, p.a0))
    dn = -n * (n - 1) * s / 2                    # dν/d ln g
    return 0.5 * (-n - dn + float(K.nu(ge, p.a0))) * ge / (n * gi)


def _einseitigkeit_voll(p):
    R = _E("R_20KPC") * K.kpc
    ge = _E("GE_003") * p.a0
    return K.einseitigkeit_m1(K.G * _m_mw() / R**2, R, ge, float(K.nu(ge, p.a0)) * ge, p.a0)


_g("M_einseitigkeit", "M", "Einseitigkeit m = 1 bei 20 kpc, g_e = 0,03·a₀ in der Ebene: erste Ordnung gegen volle Rechnung",
   lambda p: sim(lambda: _einseitigkeit_voll(p), lambda: _einseitigkeit_linear(p), 0.03), "Baustein 03", "L",
   "Hintergrund", haengt_an=(MW, "R_20KPC", "GE_003"))


# ================================================================== M26 Γ = a/c, M20 Gewicht und Gedächtnis
def _a_selbst(aN, a0):
    return a0 * brentq(lambda y: y * (1 - math.exp(-math.sqrt(y))) - aN / a0, aN / a0, 1e3 * aN / a0 + 10, xtol=1e-15,
                       rtol=1e-14)


def _steigung_selbst(p):
    M_ = _m_mw()
    r = 1e5 * K.r_stern(M_, p.a0)
    v = lambda x: math.sqrt(_a_selbst(K.G * M_ / x**2, p.a0) * x)          # noqa: E731
    return math.log(v(1.1 * r) / v(r)) / math.log(1.1)


_g("M_verlust_selbstkonsistent", "M", "Γ = a/c selbstkonsistent: d ln v/d ln r → −1/6 (fallende Kurven; die Regel nutzt Γ = a_N/c)",
   lambda p: sim(lambda: -1 / 6, lambda: _steigung_selbst(p), 5e-3), "Baustein 07", "LI", "Wirbel", haengt_an=(MW,))


@functools.lru_cache(maxsize=None)
def _maxwell_sim(x, n=2):
    """Maxwell-Element im Zeitbereich (t in t_K): σ̇ = ε̇ − σ, ε = sin(xt); (Speicher, Verlust) aus einer Periode."""
    T = 2 * math.pi / x
    s = solve_ivp(lambda t, y: [x * math.cos(x * t) - y[0]], (0, (n + 1) * T + 20), [0.0], rtol=1e-10, atol=1e-13,
                  dense_output=True)
    tt = np.linspace(20 + n * T - T, 20 + n * T, 4001)
    sg = s.sol(tt)[0]
    return 2 * np.trapezoid(sg * np.sin(x * tt), tt) / T, 2 * np.trapezoid(sg * np.cos(x * tt), tt) / T


MAXWELL_X = (0.1, 0.5, 1.0, 3.0, 17.86)


def _maxwell_abweichung(p, welche):
    d = 0.0
    for x in MAXWELL_X:
        sp, vl = _maxwell_sim(x)
        if welche == "speicher":
            d = max(d, abs(sp - gewicht_w_maxwell(x / K.tau_D(p.a0), p.a0)))
        else:
            d = max(d, abs(vl - x / (1 + x * x)))
    return 1.0 + d


_g("M_gewicht_maxwell", "M",
   "Frühere Form des Gewichts (Rechnung): Speicheranteil eines Maxwell-Elements mit t_K = τ_D (x = θτ_D; Zeitbereich)",
   lambda p: sim(lambda: 1.0, lambda: _maxwell_abweichung(p, "speicher"), 1e-6), "Rechnung, nicht im Modell", "G",
   "Hintergrund")
_g("M_gewicht_laplace", "M",
   "Gewicht der Schließungsregel aus der Rückkehr: 1 − w = Laplace-Transformierte der Dichte erster Rückkehrzeiten = e^{−√x} (x = 3)",
   lambda p: sim(lambda: 1 - K.gewicht_w(3 * p.a0 / K.c, p.a0),
                 lambda: quad(lambda t: (1 / math.sqrt(4 * math.pi * t**3)) * math.exp(-1 / (4 * t)) * math.exp(-3.0 * t),
                              0, np.inf, limit=400)[0], 1e-8), "Bausteine 03, 05", "IG", "Hintergrund")
_g("M_maxwell_verlust_zeit", "M", "Verlustanteil des Maxwell-Elements x/(1 + x²) (Zeitbereich)",
   lambda p: sim(lambda: 1.0, lambda: _maxwell_abweichung(p, "verlust"), 1e-6), "Baustein 05", "Z", "transversal")
KREIS_OM, KREIS_LAM, KREIS_R, RADIAL_K = 1.0, 0.3, 2.0, 1.3
_g("M_gedaechtnis_kreisbahn", "M", "Gedächtniskern auf der Kreisbahn: rω²/(λ_m² + ω²) = r·(Speicheranteil, x = ω/λ_m); numerisch",
   lambda p: sim(lambda: KREIS_R * KREIS_OM**2 / (KREIS_LAM**2 + KREIS_OM**2),
                 lambda: quad(lambda u: math.exp(-KREIS_LAM * u) * KREIS_R * KREIS_OM * math.sin(KREIS_OM * u), 0, np.inf,
                              limit=400)[0], 1e-8), "Baustein 07", "Z", "tangential")
_g("M_gedaechtnis_radial", "M", "Radialschwingung im Kern: Speicher κ²/(λ²+κ²), Verlust κλ/(λ²+κ²) (rundet Bahnen); numerisch",
   lambda p: sim(lambda: (RADIAL_K**2 + RADIAL_K * KREIS_LAM) / (KREIS_LAM**2 + RADIAL_K**2),
                 lambda: quad(lambda u: math.exp(-KREIS_LAM * u) * RADIAL_K * math.sin(RADIAL_K * u), 0, np.inf, limit=400)[0]
                 + quad(lambda u: math.exp(-KREIS_LAM * u) * RADIAL_K * math.cos(RADIAL_K * u), 0, np.inf, limit=400)[0],
                 1e-8), "Baustein 07", "Z", "transversal")
def _umlauf_takt_hand(a0):
    """Bei r* ist a_N = a₀, also a = a₀/(1 − e^{−1}) und v = √(a·r*): (a/a₀)·(c/v)."""
    rs = K.r_stern(_m_mw(), a0)
    a = a0 / (1 - math.exp(-1))
    return a / a0 * K.c / math.sqrt(a * rs)


_g("M_umlauf_takt", "M", "Ω·τ_D = (a/a₀)·(c/v) (Milchstraße bei r*)",
   lambda p: ident(lambda: K.v_kreis_punkt(_m_mw(), K.r_stern(_m_mw(), p.a0), p.a0) / K.r_stern(_m_mw(), p.a0)
                   * K.tau_D(p.a0), lambda: _umlauf_takt_hand(p.a0), 1e-12),
   "Baustein 07", "LI", "Wirbel", haengt_an=(MW,))


# ================================================================== Beziehungen aus dem alten Prüfstand (Beziehung → Gesetz)
_g("M_hubble_kugel_nu", "M", "An der Hubble-Kugel ist g = cH/2 = πa₀ (heute, a₀ = cH₀/2π): ν(πa₀) = 1/(1 − e^{−√π})",
   lambda p: ident(lambda: float(K.nu(K.c * K.H0_si(p.H0) / 2, K.a0_modell(p.H0))),
                   lambda: 1 / (1 - math.exp(-math.sqrt(math.pi))), 1e-12), "Baustein 03", "G", "Hintergrund",
   satz="horizont")
_g("M_preisach_hintergrund", "M", "Preisach: 1 − w im mitatmenden Hintergrund gleich der jetzigen Form (θ fällt stetig)",
   lambda p: ident(lambda: 1 - K.gewicht_w_preisach([6 * K.H0_si(p.H0), 3 * K.H0_si(p.H0)], p.a0),
                   lambda: math.exp(-math.sqrt(3 * (p.H0 * 1e3 / K.Mpc) * K.c / p.a0)), 1e-12), "Baustein 03", "G",
   "Hintergrund", satz="mitatmend_heute")
_g("M_nu_ohne_gefaelle", "M", "Ohne Gefälle ist die Rückkehr das Gewicht: ν = 1/w (θ̃ = 3H₀)",
   lambda p: ident(lambda: 1 / (1 - float(K.rueckkehr(0.0, 0.0, 3 * K.H0_si(p.H0), p.a0))),
                   lambda: 1 / K.gewicht_w(3 * K.H0_si(p.H0), p.a0), 1e-12), "Baustein 03", "IG", "Hintergrund",
   satz="mitatmend_heute")
_g("M_gewicht_eins", "M", "w bei θ = a₀/c: 1 − 1/e",
   lambda p: ident(lambda: K.gewicht_w(p.a0 / K.c, p.a0), lambda: 1 - math.exp(-1), 1e-12), "Baustein 07", "IG",
   "Hintergrund")
_g("M_gewicht_halb", "M", "w = ½ bei θ = (ln 2)²·a₀/c (Nullstelle)",
   lambda p: ident(lambda: (math.log(2))**2 * p.a0 / K.c,
                   lambda: brentq(lambda t: K.gewicht_w(t, p.a0) - 0.5, 1e-3 * p.a0 / K.c, p.a0 / K.c, xtol=1e-40,
                                  rtol=1e-14), 1e-10), "Baustein 07", "IG", "Hintergrund")
RING_R = "RING_R20"


def _ring_verhaeltnis(p):
    M_ = _E("MN_MASSE") * K.M_sun
    a, b, R = _E("MN_A") * K.kpc, _E("MN_B") * K.kpc, _E(RING_R) * K.kpc
    return K.knotenpraezession_klein(M_, a, b, R, True, p.a0) / K.knotenpraezession_klein(M_, a, b, R, False, p.a0)


def _ring_wurzel_nu(p):
    M_ = _E("MN_MASSE") * K.M_sun
    a, b, R = _E("MN_A") * K.kpc, _E("MN_B") * K.kpc, _E(RING_R) * K.kpc
    return math.sqrt(float(K.nu(K.g_miyamoto_nagai(M_, a, b, R, 0.0, "R"), p.a0)))


_g("M_ring_wurzel_nu", "M", "Ring bei Neigung → 0: ω_p örtlich/Newton = √ν (Miyamoto-Nagai, 20 kpc)",
   lambda p: ident(lambda: _ring_verhaeltnis(p), lambda: _ring_wurzel_nu(p), 1e-4), "Baustein 07", "L", "Wirbel",
   haengt_an=("MN_MASSE", "MN_A", "MN_B", RING_R))
_g("M_leiter_vflach", "M", "Leiterschritt δ* = 3Δs/r* = (v_flach/c)² (Milchstraße)",
   lambda p: ident(lambda: K.leiter_schritt(_m_mw(), p.a0), lambda: math.sqrt(K.G * _m_mw() * p.a0) / K.c**2, 1e-12),
   "Bausteine 01, 08", "LI", "Zentrum", haengt_an=(MW,))
_g("M_rstern_r0", "M", "r*/r₀ = 1/(2δ*) (Milchstraße)",
   lambda p: ident(lambda: K.r_stern(_m_mw(), p.a0) / K.r0(_m_mw()), lambda: 1 / (2 * K.leiter_schritt(_m_mw(), p.a0)),
                   1e-12), "Baustein 08", "LI", "Zentrum", haengt_an=(MW,))
_g("M_gw_schwelle", "M", "Schwellenfrequenz der Wellen f*_gw·π = Ω_N(r*) (Milchstraße)",
   lambda p: ident(lambda: K.schwellenfrequenz_gw(_m_mw(), p.a0) * math.pi,
                   lambda: K.umlauf_kepler(_m_mw(), K.r_stern(_m_mw(), p.a0)), 1e-12), "Baustein 08", "LI", "Wirbel",
   haengt_an=(MW,))
_g("M_diffusion_r0", "M", "Diffusion der Rückkehr D = r*²/τ_D = r₀c/2 = GM/c (Sonne)",
   lambda p: ident(lambda: K.r_stern(M, p.a0)**2 / K.tau_D(p.a0), lambda: K.r0(M) * K.c / 2, 1e-12), "Baustein 07", "LI",
   "Wirbel")


_g("M_zelle_hubble", "M", "Zelle in Hubble-Längen L·H₀/c = (8π²/(ν_kosm Ω_b Θ))^{1/3} (L = d_LSS, a₀ aus der Zelle)",
   lambda p: ident(lambda: K.zelle_hubble(K.theta_(K.a0_aus_zelle(K.d_LSS()), K.H0_si(p.H0))),
                   lambda: K.d_LSS() * (p.H0 * 1e3 / K.Mpc) / K.c, 1e-9), "Baustein 03", "G", "Hintergrund",
   haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "OMEGA_R"))
Z_ZUKUNFT = -1 + 1e-9                            # Beobachter im Endzustand (a → ∞)


def _lss_zukunft_hand(p):
    """2(c/H₀)∫ da/(a²E(a)) von der Rekombination bis a → ∞ (E von Hand)."""
    Om, Or = K.OMEGA_M, K.OMEGA_R
    E = lambda a: math.sqrt(Om / a**3 + Or / a**4 + 1 - Om - Or)             # noqa: E731
    a1 = 1 / (1 + EINGABEN["Z_REKOMBINATION"].wert)
    teile = quad(lambda u: 1 / (math.exp(u) * E(math.exp(u))), math.log(a1), 0, limit=400, epsrel=1e-12)[0]
    teile += quad(lambda a: 1 / (a * a * E(a)), 1, np.inf, limit=400, epsrel=1e-12)[0]
    return 2 * K.c / K.H0_si(p.H0) * teile


_g("M_zelle_kuenftig", "M", "Durchmesser der letzten Streufläche im Endzustand (Beobachter bei a → ∞): L/d_LSS fällt auf seinen Grenzwert",
   lambda p: ident(lambda: K.zelle_aus_a0(p.a0) / K.d_LSS_beobachter(Z_ZUKUNFT),
                   lambda: K.zelle_aus_a0(p.a0) / _lss_zukunft_hand(p), 1e-6), "Baustein 03", "G", "Hintergrund",
   haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "OMEGA_R", "Z_REKOMBINATION"))

# ================================================================== Leiter und Brücken der Ebenen
AUSSAGEN.neu("M_leiter", (datei("08"), "r*² = 3Δs·ℓ₀"), "Beziehung")


def _leiter_massen():
    return [K.m_p, K.M_sun, _m_mw(), _E("M_HAUFEN") * K.M_sun]


_g("M_leiter_ebenen", "M", "Leiter r*² = 3Δs·ℓ₀ auf jeder Ebene (Proton, Sonne, Milchstraße, Haufen): r* ist das "
   "geometrische Mittel von Faltung und Wirkungslänge",
   lambda p: ident(lambda: np.array([K.r_stern(m, p.a0)**2 for m in _leiter_massen()]),
                   lambda: np.array([3 * K.ds(m) * K.ell0(p.a0) for m in _leiter_massen()]), 1e-12),
   "Bausteine 01, 08", "LI", "Zentrum", haengt_an=(MW, "M_HAUFEN"), bezug="M_leiter")
for _i, (_von, _nach, _gl, _f, _ort) in enumerate(BRUECKEN):
    AUSSAGEN.neu(f"M_bruecke_{_von}{_nach}_{_i}", _ort, "Beziehung")
_g("M_bruecke_LI_r0_lambda", "M", "Brücke L–I: r*² = r₀·Λ (Raum und innere Richtung)",
   lambda p: ident(lambda: K.r_stern(_m_mw(), p.a0)**2, lambda: K.r0(_m_mw()) * K.schliesslaenge(p.a0), 1e-12),
   "Baustein 07", "LI", "Zentrum", haengt_an=(MW,), bezug="M_bruecke_LI_1")
_g("M_bruecke_LI_gamma", "M", "Brücke L–I: Γ = a_N/c = D·k² (Glättungsrate mit D = GM/c, k = 1/r; Sonne bei 1 AE)",
   lambda p: ident(lambda: K.verlust_gamma(GM / K.AU**2), lambda: (GM / K.c) / K.AU**2, 1e-12), "Baustein 07", "LI",
   "Wirbel", bezug="M_bruecke_LI_2")


# ================================================================== Rollen der Größen und Grenzlängen (Register des Kerns)
for _i, (_z, _tri, _o1, _rol, _o2) in enumerate(GROESSEN_ROLLEN):
    for _j, _o in enumerate((_o1, _o2)):
        if _o:
            AUSSAGEN.neu(f"M_rolle_{_i}_{_j}", _o, "Beziehung")
AUSSAGEN.neu("M_grenzlaengen", GRENZLAENGEN_ORT, "Beziehung")
_g("M_groessen_rollen", "M", "Jede Größe des Modells trägt höchstens ein Glied der Triade und höchstens eine Rolle im Zustandsraum; "
   "Zeichen stehen im Symbolregister",
   lambda p: ident(lambda: np.array([float(_z in SYMBOLE or _z in ("|q|", "arg q", "Q_el", "Re τ_M", "Δω_R"))
                                     * float(_tri is None or _tri in TRIADE_GLIEDER)
                                     * float(_rol is None or _rol in ROLLEN_ZUSTANDSRAUM)
                                     for _z, _tri, _o1, _rol, _o2 in GROESSEN_ROLLEN]),
                   lambda: np.ones(len(GROESSEN_ROLLEN)), 0.0),
   "Bausteine 01, 05", "LIZG", "Wirbel", bezug="M_rolle_0_0")
_g("M_grenzlaengen_ebenen", "M", "Grenzlängen jeder Ebene 3Δs, r*, ℓ₀: Faltung und Schwellenradius im Raum (L), Wirkungslänge in der "
   "inneren Richtung (I); ihre Schritte sind gleich (δ* = 3Δs/r* = r*/ℓ₀)",
   lambda p: ident(lambda: np.array([K.leiter_schritt(m, p.a0) for m in _leiter_massen()]),
                   lambda: np.array([K.r_stern(m, p.a0) / K.ell0(p.a0) for m in _leiter_massen()]), 1e-12),
   "Baustein 08", "LI", "Zentrum", haengt_an=(MW, "M_HAUFEN"), bezug="M_grenzlaengen")


def _zelle_faltung(p):
    """ν_kosm·ΣΔs_b der Zelle L = zelle_aus_a0(a₀): Faltung aller Baryonen der Zelle mal ν_kosm = 1 + Q."""
    L = K.zelle_aus_a0(p.a0)
    return (1 + K.Q_KOSM) * K.ds(K.rho_b_mitbewegt() * L**3)


_g("M_bruecke_IG_lambda", "M", "Brücke I–G: Λ = ν_kosm·ΣΔs_b der Zelle, ν_kosm = 1 + Q",
   lambda p: ident(lambda: K.schliesslaenge(p.a0), lambda: _zelle_faltung(p), 1e-12), "Baustein 03", "IG", "Zentrum",
   haengt_an=("OMEGA_B_H2", "OMEGA_C_H2", "H0_KMS"), bezug="M_bruecke_IG_3")
_g("M_bruecke_LG_RQ", "M", "Brücke L–G: bei R_Q = r*/(−ln q_kosm) erreicht das örtliche q das kosmische",
   lambda p: ident(lambda: K.q_(K.G * _m_mw() / K.r_Q(_m_mw(), p.a0)**2, p.a0), lambda: K.q_KOSM, 1e-12), "Baustein 03",
   "LG", "Hintergrund", haengt_an=(MW, "OMEGA_B_H2", "OMEGA_C_H2"), bezug="M_bruecke_LG_4")


def _bruecke_zl():
    """Brücke Z–L: 𝒫_q = (1 − q)(I − qT)⁻¹ (Regel, Mittel über die Translation T) vertauscht mit dem Erzeuger L = T − Tᵀ,
    also L_eff = 𝒫L = L𝒫; im Grenzfall q → 1 mittelt 𝒫 die Nullmode heraus (𝒫L → 0, die Regel löscht dort nur den
    Erzeuger, keine Freiheitsgrade). Zweiter Weg: Eigenwerte im Fourier-Raum, 𝒫 → (1 − q)/(1 − q e^{−ik})."""
    n, q = 32, 0.6
    T = np.roll(np.eye(n), 1, 0)
    L = T - T.T
    P = (1 - q) * np.linalg.inv(np.eye(n) - q * T)
    k = 2 * np.pi * np.arange(n) / n
    ew = (1 - q) / (1 - q * np.exp(-1j * k)) * (np.exp(-1j * k) - np.exp(1j * k))
    ew_leff = np.linalg.eigvals(P @ L)
    P1 = (1 - 0.999999) * np.linalg.inv(np.eye(n) - 0.999999 * T)
    return wahr(np.allclose(P @ L, L @ P) and all(np.min(np.abs(ew_leff - e)) < 1e-9 for e in ew) and np.linalg.norm(P1 @ L) < 1e-4,
                f"‖[𝒫, L]‖ = {np.linalg.norm(P @ L - L @ P):.1e}")


_g("M_bruecke_ZL_regel", "M", "Brücke Z–L: Der Zustandsraum wirkt auf den Ort über die Regel, L_eff = 𝒫L "
   "(𝒫 vertauscht mit dem Erzeuger, löscht keine Freiheitsgrade)", lambda p: _bruecke_zl(), "Baustein 05", "ZL",
   "transversal", bezug="M_bruecke_ZL_5")
_g("M_bruecke_ZG_kac", "M", "Brücke Z–G: Der Zustandsraum bindet sich ans Ganze über das Kac-Budget, "
   "Σ_i 1/⟨T_{A_i}⟩ = 1 für jede Zerlegung des geschlossenen Raums",
   lambda p: ident(lambda: K.kac_budget(7.0, [1.0, 2.5, 3.5]), lambda: 1.0, 1e-12),
   "Baustein 06", "ZG", "transversal", bezug="M_bruecke_ZG_6")


# ================================================================== Budgets je Ebene (01: kein Außen)
def _gauss_torus():
    """Poisson auf dem Torus (periodisch, FFT): lösbar genau für Σρ = 0 (Windungen heben sich auf); mit Σρ ≠ 0 bleibt
    der Nullmode als Rest."""
    n = 64
    rng = np.random.default_rng(SEED + 1)
    k = np.fft.fftfreq(n) * 2 * math.pi
    kx, ky = np.meshgrid(k, k, indexing="ij")
    k2 = kx**2 + ky**2
    rest = []
    for mittel in (0.0, 0.1):
        rho = rng.normal(size=(n, n))
        rho += mittel - rho.mean()
        rh = np.fft.fft2(rho)
        phi = np.where(k2 == 0, 0.0, -rh / np.where(k2 == 0, 1.0, k2))
        lap = np.real(np.fft.ifft2(-k2 * phi))
        rest.append(np.abs(lap - rho).max())
    return wahr(rest[0] < 1e-12 and abs(rest[1] - 0.1) < 1e-12, f"Rest {rest[0]:.0e} / {rest[1]:.2f}")


_g("M_budget_ladung", "M", "Kein Außen: auf dem geschlossenen Raum verlangt Gauß Σ Quelle = 0 (Σ Q_el = 0)",
   lambda p: _gauss_torus(), "Baustein 01", "G", "Hintergrund")


def _faltung_horizont_kugel():
    R = K.teilchenhorizont(0.0)
    return K.ds(K.rho_b_mitbewegt() * 4 / 3 * math.pi * R**3) / R


_g("M_budget_rueckstellung", "M", "Rückstellung σ = ΣΔs − ΔL = 0: Dehnung der Raumschleifen je Radius = Faltung der Baryonen je Radius",
   lambda p: ident(lambda: K.rueckstellbilanz(0.0), _faltung_horizont_kugel, 1e-12), "Bausteine 01, 03", "G", "Zentrum",
   haengt_an=("OMEGA_B_H2", "H0_KMS", "OMEGA_C_H2", "OMEGA_R"))


# ================================================================== Grenzfälle Sonnensystem (03, 07)
def _sattel(p):
    """Sattelpunkt Sonne–Erde liegt in beiden Glättungszonen: dort ist das Feld Newtonsch (kein Zusatz)."""
    T = gezeit(p)
    w = math.sqrt(K.M_earth / K.M_sun)
    d = K.AU * w / (1 + w)                                 # Abstand von der Erde
    lE, lS = K.glaettungslaenge(K.M_earth, T), K.glaettungslaenge(K.M_sun, T)
    return wahr(d < lE and K.AU - d < lS, f"Sattel {d:.3g} m, ℓ_Erde {lE / K.AU:.0f} AE, ℓ_Sonne {lS / K.AU:.0f} AE")


_g("M_sattelpunkt_newton", "M", "Sattelpunkt Sonne–Erde ohne Zusatzbeschleunigung: er liegt in den Glättungszonen beider Körper",
   _sattel, "Bausteine 03, 07", "L", "Hintergrund", haengt_an=SONNENORT, satz="Sonnenkreis_gemessen")
_g("M_glaettung_schwelle", "M", "Glättung: ℓ³ = GM/𝒯 und ℓ/r* = (M×/M)^{1/6} mit M× = a₀³/(G𝒯²) (Sonne im Feld der Galaxie)",
   lambda p: ident(lambda: K.glaettungslaenge(K.M_sun, gezeit(p)) / K.r_stern(K.M_sun, p.a0),
                   lambda: K.glaettung_zu_schwelle(K.M_sun, gezeit(p), p.a0), 1e-12), "Baustein 07", "L", "Hintergrund",
   haengt_an=SONNENORT, satz="Sonnenkreis_gemessen")


# ================================================================== Erhaltung in der Feldform (07): V, ∂_i𝒮_ij = 0
FF_N, FF_BREITE, FF_ORT1, FF_ORT2, FF_MASSE2 = 128, 0.06, (0.35, 0.45), (0.62, 0.55), 3.0   # Gitter auf dem Torus 1 × 1


class _Feldform:
    """Zwei glatte Körper auf dem periodischen Gitter (2D, FFT): Newton-Feld ∇Φ_N, Feldform ∇²Φ = ∇·[ν(|∇Φ_N|)∇Φ_N]
    und die Kräfte auf beide Körper; die Schwelle liegt mitten im Feld (Tief- und Newton-Bereich zugleich)."""

    def __init__(self):
        n = FF_N
        x = (np.arange(n) + 0.5) / n
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")
        k = np.fft.fftfreq(n, 1 / n) * 2 * math.pi
        self.kx, self.ky = np.meshgrid(k, k, indexing="ij")
        self.k2 = np.where(self.kx**2 + self.ky**2 == 0, 1.0, self.kx**2 + self.ky**2)
        self.r1 = self.koerper(*FF_ORT1, 1.0)
        gx, gy = self.grad(self.poisson(self.r1 + self.koerper(*FF_ORT2, FF_MASSE2)))
        self.a0 = float(np.median(np.hypot(gx, gy)))

    def koerper(self, x0, y0, m):
        return m * np.exp(-((self.X - x0)**2 + (self.Y - y0)**2) / (2 * FF_BREITE**2))

    def poisson(self, f):
        h = np.fft.fft2(f)
        h[0, 0] = 0.0                              # Torus: nur die Abweichung vom Mittel hat ein Feld
        return -h / self.k2

    def grad(self, ph):
        return np.real(np.fft.ifft2(1j * self.kx * ph)), np.real(np.fft.ifft2(1j * self.ky * ph))

    def kraefte(self, x2, y2):
        """(Feldform auf 1, auf 2, örtliche Regel auf 1, auf 2) mit Körper 2 bei (x2, y2)."""
        r2 = self.koerper(x2, y2, FF_MASSE2)
        gx, gy = self.grad(self.poisson(self.r1 + r2))
        nu = K.nu(np.maximum(np.hypot(gx, gy), 1e-300), self.a0)
        div = np.real(np.fft.ifft2(1j * self.kx * np.fft.fft2(nu * gx) + 1j * self.ky * np.fft.fft2(nu * gy)))
        Fx, Fy = self.grad(self.poisson(div))
        f = [np.array([np.sum(r * Fx), np.sum(r * Fy)]) for r in (self.r1, r2)]
        o = [np.array([np.sum(r * nu * gx), np.sum(r * nu * gy)]) for r in (self.r1, r2)]
        return f[0], f[1], o[0], o[1]


def _feldform_impuls(p):
    """Die Kräfte der Feldform heben sich auf (Impuls erhalten, bis auf die Auflösung des Gitters: 128² gibt 10⁻⁴,
    256² 10⁻⁵); die örtliche Regel je Körper tut es nicht."""
    f1, f2, o1, o2 = _Feldform().kraefte(*FF_ORT2)
    feld = np.linalg.norm(f1 + f2) / np.linalg.norm(f1)
    oertlich = np.linalg.norm(o1 + o2) / np.linalg.norm(o1)
    return wahr(feld < 1e-3 and oertlich > 0.05, f"Feldform {feld:.0e}, örtliche Regel {oertlich:.0e}")


_g("M_feldform_impuls", "M",
   "Feldform V: Impuls erhalten (Gesamtkraft null, ∂_i𝒮_ij = 0 im Newton-Grenzfall); die örtliche Regel allein nicht",
   _feldform_impuls, "Baustein 07", "L", "Wirbel")
FF_KREIS, FF_SCHRITTE = 0.12, 32


def _feldform_energie(p):
    """Arbeit der Feldform-Kraft auf Körper 2 längs eines geschlossenen Kreises um seinen Ort: null (Kraft = −∇V)."""
    ff = _Feldform()
    arbeit, betrag = 0.0, 0.0
    for i in range(FF_SCHRITTE):
        ph = 2 * math.pi * i / FF_SCHRITTE
        f2 = ff.kraefte(FF_ORT2[0] + FF_KREIS * math.cos(ph), FF_ORT2[1] + FF_KREIS * math.sin(ph))[1]
        arbeit += f2 @ (np.array([-math.sin(ph), math.cos(ph)]) * FF_KREIS * 2 * math.pi / FF_SCHRITTE)
        betrag += np.linalg.norm(f2) / FF_SCHRITTE
    rel = abs(arbeit) / (betrag * 2 * math.pi * FF_KREIS)
    return wahr(rel < 1e-3, f"Umlaufarbeit / (|F|·Umfang) {rel:.0e}")


_g("M_feldform_energie", "M", "Feldform V: Energie erhalten, die Kraft ist ein Gradient (Arbeit auf geschlossenem Weg null)",
   _feldform_energie, "Baustein 07", "L", "Wirbel")
