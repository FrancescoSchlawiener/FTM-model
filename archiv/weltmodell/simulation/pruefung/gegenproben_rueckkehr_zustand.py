# -*- coding: utf-8 -*-
"""
Gegenproben der Rückkehr (I), der Faltung am Ort (L) und des Zustandsraums (Z): je Formel aus kern/rueckkehr.py,
kern/faltung.py und kern/zustandsraum.py mindestens ein unabhängiger Weg (zweiter Rechenweg, Grenzfall, Handwert).
Der unabhängige Weg ruft die geprüfte Formel nicht auf (Wächter). Zufall nur mit eigenem Generator mit Seed je Aufruf.
Handwerte aus dem Modelltext tragen die Toleranz ihrer letzten Stelle.
"""
import math

import numpy as np

import kern as K

from .gegenprobe import GEGENPROBEN

_neu = GEGENPROBEN.neu
A0 = K.A0
M_MW = 6e10 * K.M_sun
GOLDEN = (math.sqrt(5) - 1) / 2


def _nu(s):
    """ν = 1/(1 − e^{−s}) ohne Kern (für die Differenzenwege)."""
    return 1 / -math.expm1(-s)


def _a_punkt(r, GM=1.0, a0=1.0):
    """Punktmasse mit der örtlichen Regel, ohne Kern: a = a_N/(1 − e^{−√(a_N/a₀)})."""
    aN = GM / r**2
    return aN * _nu(math.sqrt(aN / a0))


def _laplace_dichte(lam):
    """∫₀^∞ √(1/4π)·x^{−3/2}·e^{−1/(4x)}·e^{−λx} dx (Dichte erster Rückkehrzeiten in Einheiten τ_D), über x = e^u."""
    from scipy.integrate import quad
    f = lambda u: math.exp(-u / 2 - math.exp(-u) / 4 - lam * math.exp(u)) / math.sqrt(4 * math.pi)
    return quad(f, -40, 40, limit=400, points=[-math.log(max(lam, 1e-30)) / 2], epsabs=0, epsrel=1e-12)[0]


# ================================================================== Rückkehr (I)
def _taktgefaelle(aN):
    """Γ = c·|∇ln(dτ/dt)| an der Sonne im Abstand mit a_N, dτ/dt = √(1 − 2GM/(rc²)), Differenzenquotient."""
    GM = K.G * K.M_sun
    r = math.sqrt(GM / aN)
    h = 1e-4 * r
    ln_takt = lambda x: 0.5 * math.log1p(-2 * GM / (x * K.c**2))
    return K.c * (ln_takt(r + h) - ln_takt(r - h)) / (2 * h)


_neu("verlust_gamma_taktgefaelle", "verlust_gamma", "zweiter_weg", [(1e-12,), (1.1e-10,), (1e-8,)], _taktgefaelle,
     toleranz=1e-6, grund="Γ = c·|∇ln(dτ/dt)| als Differenzenquotient (h = 10⁻⁴r), schwaches Feld 6Δs/r ≪ 10⁻⁶")

_neu("tau_D_rueckkehrtiefe", "tau_D", "zweiter_weg", [(A0,), (1.2e-10,)],
     lambda a0: K.r_stern(K.M_sun, a0)**2 / (K.G * K.M_sun / K.c))


def _exponent_minimum(p, Gamma, theta, a0):
    """s = min_t[(p + Γ + θ)t + τ_D/(4t)]: Verlust gegen Rückkehrzeit, Minimum über ln t."""
    from scipy.optimize import minimize_scalar
    rate, tD = p + Gamma + theta, K.c / a0
    u0 = 0.5 * math.log(tD / (4 * rate))
    r = minimize_scalar(lambda u: rate * math.exp(u) + tD / (4 * math.exp(u)), bounds=(u0 - 3, u0 + 3),
                        method="bounded", options=dict(xatol=1e-10))
    return r.fun


_neu("rueckkehr_exponent_minimum", "rueckkehr_exponent", "zweiter_weg",
     [(0.0, 1e-18, 0.0, A0), (1e-19, 3e-19, 2e-19, A0), (0.0, 4e-17, 0.0, 1.2e-10)], _exponent_minimum,
     toleranz=1e-9, grund="numerisches Minimum (xatol 10⁻¹⁰ in ln t)")

_neu("rueckkehr_laplace", "rueckkehr", "zweiter_weg",
     [(0.0, 1e-18, 0.0, A0), (1e-18, 1e-19, 3e-19, A0), (0.0, 4 * A0 / K.c, 0.0, A0), (0.0, 0.01 * A0 / K.c, 0.0, A0)],
     lambda p, Gamma, theta, a0: _laplace_dichte((p + Gamma + theta) * K.c / a0),
     toleranz=1e-9, grund="Laplace-Integral der Dichte erster Rückkehrzeiten mit quad")

_neu("s_punktmasse", "s_", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     lambda aN, a0: K.r_stern(K.M_sun, a0) / math.sqrt(K.G * K.M_sun / aN))

_neu("q_bose", "q_", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     lambda aN, a0: (lambda nb: nb / (1 + nb))(float(K.n_B(math.sqrt(aN / a0)))))

_neu("nu_bose", "nu", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10), (1e-5, A0)],
     lambda aN, a0: 1 + float(K.n_B(math.sqrt(aN / a0))))
_neu("nu_grenzfall_tief", "nu", "grenzfall", [(1e-12 * A0, A0), (1e-14 * A0, A0)],
     lambda aN, a0: 1 / math.sqrt(aN / a0), toleranz=1e-6, grund="s → 0: ν → 1/s, Rest s/2 (s ≤ 10⁻⁶)")
_neu("nu_grenzfall_newton", "nu", "grenzfall", [(1600 * A0, A0), (2500 * A0, A0)], lambda aN, a0: 1.0)

_neu("chi_nome", "chi", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     lambda aN, a0: (lambda s, q: s * q / (1 - q))(math.sqrt(aN / a0), float(K.q_(aN, a0))))
_neu("chi_elastizitaet", "chi", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     lambda aN, a0: (lambda s, h: -(math.log(_nu(s * (1 + h))) - math.log(_nu(s * (1 - h))))
                     / (math.log(1 + h) - math.log(1 - h)))(math.sqrt(aN / a0), 1e-5),
     toleranz=1e-8, grund="χ = −d ln ν/d ln s als Differenzenquotient (h = 10⁻⁵)")

_neu("n_B_rueckkehrsumme", "n_B", "zweiter_weg", [(0.3,), (1.0,), (3.0,)],
     lambda s: sum(math.exp(-n * s) for n in range(1, int(45 / s) + 2)), toleranz=1e-12,
     grund="Reihe Σ_{n≥1} e^{−ns} abgeschnitten bei e^{−45}")

_neu("Q_tilde_vakuum", "Q_tilde", "handwert", [(0.0,)], lambda s: -4 * math.pi**4 / 15)
def _q_tilde_reihe(s):
    """Q̃ = −4Σ_{n≥1} Γ(4, ns)/n⁴ (obere unvollständige Gammafunktion)."""
    from scipy.special import gammaincc
    n = np.arange(1, 4000)
    return float(-4 * np.sum(6 * gammaincc(4, n * s) / n**4.0))


_neu("Q_tilde_gammareihe", "Q_tilde", "zweiter_weg", [(0.5,), (1.0,), (3.0,)], _q_tilde_reihe,
     toleranz=1e-9, grund="Reihe −4Σ Γ(4, ns)/n⁴ bis n = 4000 gegen quad")

_neu("schliesslaenge_kac", "schliesslaenge", "zweiter_weg", [(A0,), (1.2e-10,)], lambda a0: K.c * K.tau_D(a0) / 2)
_neu("ell0_leiter", "ell0", "umkehr", [(A0,), (1.2e-10,)],
     lambda a0: K.r_stern(K.M_sun, a0)**2 / (3 * K.ds(K.M_sun)))
_neu("a0_aus_schliesslaenge_kac", "a0_aus_schliesslaenge", "zweiter_weg", [(4.1e26,), (3.7e26,)],
     lambda Lam: K.c / K.kac_zeit_schluss(Lam))

_neu("gewicht_w_laplace", "gewicht_w", "zweiter_weg", [(3.67e-19, A0), (1e-19, A0), (6.6e-18, A0), (1e-18, 1.2e-10)],
     lambda theta, a0: 1 - _laplace_dichte(theta * K.c / a0),
     toleranz=1e-9, grund="1 − w = E[e^{−θT}] als Laplace-Integral mit quad")
_neu("gewicht_w_halb", "gewicht_w", "handwert", [(math.log(2)**2 * A0 / K.c, A0)], lambda theta, a0: 0.5)
_neu("gewicht_w_preisach_minimum", "gewicht_w_preisach", "zweiter_weg",
     [([1e-17, 3e-19, 5e-18], A0), ([-1e-18, 1e-17], A0), ([4e-18, 2e-18, 9e-18, 2.5e-18], 1.2e-10)],
     lambda folge, a0: -math.expm1(-math.sqrt(max(min(folge), 0.0) * K.c / a0)))

_neu("kern_zeit_ueberleben", "kern_zeit", "zweiter_weg",
     [(1e17, 1e-19, A0), (2.7e18, 3e-19, A0), (5e19, 1e-20, 1.2e-10)],
     lambda T, Gamma, a0: (lambda tD, h: -(math.erf(math.sqrt(tD / (4 * (T + h)))) - math.erf(math.sqrt(tD / (4 * (T - h)))))
                           / (2 * h) * math.exp(-Gamma * T))(K.c / a0, 1e-5 * T),
     toleranz=1e-8, grund="Dichte = −d/dT erf(ϖ_R/√(4DT)) als Differenzenquotient (h = 10⁻⁵T)")

_neu("kern_bahnfrequenz_wurzel", "kern_bahnfrequenz", "zweiter_weg",
     [(7.8e-16, 2e-18, A0), (1e-17, 1e-18, A0), (1e-18, 4e-19, 1.2e-10)],
     lambda Om, Gamma, a0: (lambda x, y: math.exp(-math.sqrt((math.hypot(x, y) + x) / 2) + math.sqrt(x)))(
         Gamma * K.c / a0, Om * K.c / a0),
     toleranz=1e-9, grund="Re√(x + iy) = √((|z| + x)/2), Exponent bis 50")


def _glaettung_randwert(k, varpi, n=4000):
    """Beschränkte Lösung von f'' = k²f, f(0) = 1, auf [0, 40/k] (f(40/k) = 0), Differenzengleichung."""
    from scipy.linalg import solve_banded
    L = 40 / abs(k)
    h = L / n
    ab = np.zeros((3, n - 1))
    ab[0, 1:] = 1.0
    ab[1, :] = -(2 + (k * h)**2)
    ab[2, :-1] = 1.0
    b = np.zeros(n - 1)
    b[0] = -1.0
    f = np.concatenate([[1.0], solve_banded((1, 1), ab, b), [0.0]])
    return float(np.interp(varpi, np.linspace(0, L, n + 1), f))


_neu("glaettung_innen_randwert", "glaettung_innen", "zweiter_weg", [(1.0, 0.5), (2.0, 1.0), (1e-20, 3e19)],
     _glaettung_randwert, toleranz=2e-4, grund="Differenzengleichung h = 40/(4000k), Fehler (kh)²/12 je Länge")


def _erneuerung_randwert(gN, varpi, M, a0):
    """∂_ϖ²u = (Γ/D)u auf [0, W], u(W) = 0, Erneuerung u(0) − u(ϖ_R) = g_N; ϖ_R² = D·τ_D mit D = GM/c."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    k = math.sqrt(abs(gN) * K.c / (K.G * M) / K.c)
    vR = math.sqrt(K.G * M / K.c * K.c / a0)
    m = 2000
    h = vR / m
    W = max(40 / k, 2 * vR, 1.2 * varpi)
    n = int(W / h) + 1
    A = diags([np.ones(n), -(2 + (k * h)**2) * np.ones(n + 1), np.ones(n)], [-1, 0, 1], format="lil")
    A[0, :] = 0
    A[0, 0], A[0, m] = 1.0, -1.0
    A[n, :] = 0
    A[n, n] = 1.0
    b = np.zeros(n + 1)
    b[0] = gN
    u = spsolve(A.tocsr(), b)
    return float(np.interp(varpi, np.arange(n + 1) * h, u))


_neu("erneuerung_randwert", "erneuerung", "zweiter_weg",
     [(-1e-10, 0.0, M_MW, A0), (-1e-10, 1e19, M_MW, A0), (-5e-11, 0.0, 1e10 * K.M_sun, 1.2e-10)],
     _erneuerung_randwert, toleranz=1e-5, grund="Differenzengleichung mit 2000 Schritten bis ϖ_R")


def _eigenwert_matrix(D, Lam, n=3000):
    """Kleinster Eigenwert von −D∂² auf [0, Λ], Ort fest (Dirichlet), bei Λ reflektierend (Neumann in der Zellmitte,
    Knoten x_i = i·h, Λ = (n + ½)h)."""
    from scipy.linalg import eigvalsh_tridiagonal
    h = Lam / (n + 0.5)
    d = np.full(n, 2.0)
    d[-1] = 1.0
    return D / h**2 * float(eigvalsh_tridiagonal(d, -np.ones(n - 1), select="i", select_range=(0, 0))[0])


_neu("eigenwert_schluss_matrix", "eigenwert_schluss", "zweiter_weg", [(1.0, 1.0), (4.4e22, 4e26)], _eigenwert_matrix,
     toleranz=1e-6, grund="Differenzenmatrix mit 3000 Punkten, Fehler (πh/2Λ)²/12")


def _matsubara(s, L, n=4000):
    """Thermischer Propagator der Periode β = 2Λ (Einheit ϖ_R), G(τ) = (1/β)Σ_n cos(ω_nτ)/(ω_n² + k²), k = s;
    Glied ohne k als Bernoulli-Polynom, Rest ∝ n⁻⁴. Gibt G(ϖ_R)/G(0)."""
    beta, k = 2 * L, s
    w = 2 * math.pi * np.arange(1, n) / beta

    def G(tau):
        x = tau / beta
        grund = beta**2 / 2 * (x * x - x + 1 / 6)
        rest = 2 * np.sum(np.cos(w * tau) / (w**2 * (w**2 + k**2)))
        return 1 / k**2 + grund - k**2 * rest
    return G(1.0) / G(0.0)


_neu("rueckkehr_geschlossen_matsubara", "rueckkehr_geschlossen", "zweiter_weg",
     [(1.0, 3.0), (0.2, 10.0), (0.05, 31.4), (3.162e-3, 100.0)], _matsubara,
     toleranz=1e-8, grund="Matsubara-Summe bis n = 4000, Rest ∝ n⁻³")
_neu("rueckkehr_geschlossen_irrweg", "rueckkehr_geschlossen", "handwert", [(math.sqrt(1e-7 / 0.5) * 10, 100.0)],
     lambda s, L: 0.998133, toleranz=6e-7,
     grund="Irrweg ϖ = 10, Λ = 1000, Γ = 10⁻⁷, D = ½ (Rechnung wicklung_rueckkehr, 6 Stellen)")

_neu("abweichung_schluss_grenzfall", "abweichung_schluss", "grenzfall", [(math.pi / 2,), (1.0,), (3.0,)],
     lambda kL: (1 - math.exp(-1e-6)) / (1 - K.rueckkehr_geschlossen(1e-6, kL / 1e-6)) - 1,
     toleranz=1e-4, grund="s = 10⁻⁶ statt s → 0 bei festem kΛ")
_neu("abweichung_schluss_modelltext", "abweichung_schluss", "handwert", [(math.pi / 2,)], lambda kL: 0.090,
     toleranz=0.006, grund="9,0 % (halbe letzte Stelle 0,05 %)")


def _rc_offen(kL):
    """RC-Kette mit offenem Ende bei n_Λ = kΛ·λ Gliedern (λ = 1/√(RG) = 2000), Abgriff n_R = 2 auf den Eingang
    zurückgeführt: Abweichung ν_offen/ν_unendlich − 1, beide Ketten als Netz gelöst."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    RG, nR = 1 / 2000**2, 2

    def nu_netz(n):
        A = diags([np.ones(n), -(2 + RG) * np.ones(n + 1), np.ones(n)], [-1, 0, 1], format="lil")
        A[0, :] = 0
        A[0, 0], A[0, nR] = 1.0, -1.0
        A[n, n - 1], A[n, n] = 1.0, -(1 + RG)
        b = np.zeros(n + 1)
        b[0] = 1.0
        return spsolve(A.tocsr(), b)[0]
    return nu_netz(int(round(kL * 2000))) / nu_netz(40 * 2000) - 1


_neu("abweichung_schluss_rc_kette", "abweichung_schluss", "zweiter_weg", [(math.pi / 2,)], _rc_offen,
     toleranz=2e-2, grund="Kette mit Abgriff n_R = 2 bei λ = 2000 (s = 10⁻³), Ende auf ganze Glieder gerundet")


def _radius_schluss(Lam, rstern, vR):
    from scipy.optimize import brentq
    GM, tD = rstern**2 * A0, K.c / A0
    lam = vR**2 / tD * math.pi**2 / (4 * Lam**2)
    x = brentq(lambda u: math.log(GM / (math.exp(2 * u) * K.c)) - math.log(lam), -200, 200, xtol=1e-14)
    return math.exp(x)


_neu("radius_schluss_eigenwert", "radius_schluss", "zweiter_weg",
     [(4.1e26, 2.7e20, 2.7e20), (4.1e26, 1e15, 3e15)], _radius_schluss,
     toleranz=1e-9, grund="Nullstelle Γ(r) = Dπ²/(4Λ²) mit brentq")

_neu("temperatur_innen_periode", "temperatur_innen", "zweiter_weg", [(A0,), (1.2e-10,)],
     lambda a0: K.hbar / (K.kB * 2 * K.schliesslaenge(a0) / K.c))
_neu("temperatur_innen_handwert", "temperatur_innen", "handwert", [(A0,)], lambda a0: 2.803e-30,
     toleranz=2e-4, grund="2,803·10⁻³⁰ K (Rechnung triade_ebenen, 4 Stellen)")


def _kac_halbgerade(Lam, n=200):
    """Irrweg auf 0 … N, an beiden Enden reflektierend, Schrittdauer ℓ/c mit ℓ = Λ/N: stationäres Maß als Eigenvektor,
    ⟨T₀⟩ = 1/π₀ Schritte."""
    P = np.zeros((n + 1, n + 1))
    P[0, 1] = P[n, n - 1] = 1.0
    for i in range(1, n):
        P[i, i - 1] = P[i, i + 1] = 0.5
    w, v = np.linalg.eig(P.T)
    pi = np.real(v[:, np.argmin(abs(w - 1))])
    pi /= pi.sum()
    return 1 / pi[0] * (Lam / n) / K.c


_neu("kac_zeit_schluss_markow", "kac_zeit_schluss", "zweiter_weg", [(4.1e26,), (1.0,)], _kac_halbgerade,
     toleranz=1e-9, grund="Eigenvektor der Übergangsmatrix (201 Zustände)")

_neu("strukturteil_Q_tilde", "strukturteil", "zweiter_weg", [(0.01,), (1.0,), (50.0,)],
     lambda x: K.Q_tilde(x**0.25) + 4 * math.pi**4 / 15, toleranz=1e-8, grund="zwei quad-Integrale")

# ================================================================== Faltung am Ort (L)
_neu("neigung_xy_newton", "neigung_xy", "zweiter_weg", [(1.0, 1.0), (317.8, 11.2)],
     lambda x, y: K.G * x * K.M_earth / (y * K.R_earth)**2, toleranz=2e-3,
     grund="g⊕ = 9,81 m/s² gesetzt, GM⊕/R⊕² = 9,798 m/s²")
_neu("ds_leiter", "ds", "umkehr", [(K.m_p,), (K.M_sun,), (M_MW,)],
     lambda M: K.r_stern(M)**2 / K.ell0() / 3)
_neu("neigung_faltung", "neigung", "zweiter_weg", [(K.M_sun, K.R_sun), (K.M_earth, K.R_earth)],
     lambda M, R: 3 * K.c**2 * K.ds(M) / R**2)
_neu("takt_metrik", "takt", "zweiter_weg", [(K.M_sun, K.R_sun), (K.M_sun, 1e4)],
     lambda M, r: math.sqrt(1 - 2 * K.G * M / (r * K.c**2)))
_neu("delta_metrik", "delta", "zweiter_weg", [(K.M_sun, K.R_sun), (K.M_sun, 1e4), (K.m_p, 1e-15)],
     lambda M, r: -math.expm1(0.5 * math.log1p(-2 * K.G * M / (r * K.c**2))))
_neu("delta_1fm_hand", "delta", "handwert", [(K.m_p, 1e-15)], lambda M, r: 1.2421e-39, toleranz=1e-4,
     grund="bei 1 fm neben dem Proton von Hand δ ≈ GM_p/(rc²) = 6,6743·10⁻¹¹·1,67262·10⁻²⁷/(10⁻¹⁵·8,98755·10¹⁶) "
           "= 1,2421·10⁻³⁹ (führende Ordnung, der Rest ist 10⁻³⁹-mal kleiner); vier Stellen")
_neu("r0_schwarzschild", "r0", "zweiter_weg", [(K.M_sun,), (K.M_earth,)], lambda M: 2 * K.G * M / K.c**2)
_neu("r_stern_leiterweg", "r_stern", "umkehr", [(K.M_sun, A0), (M_MW, 1.2e-10)], K.r_stern_leiter)
_neu("r_stern_leiter_direkt", "r_stern_leiter", "zweiter_weg", [(K.M_sun, A0), (M_MW, 1.2e-10)],
     lambda M, a0: math.sqrt(K.G * M / a0))
_neu("leiter_schritt_innen", "leiter_schritt", "zweiter_weg", [(K.m_p, A0), (K.M_sun, A0), (M_MW, 1.2e-10)],
     lambda M, a0: K.r_stern(M, a0) / K.ell0(a0))
_neu("leiter_schritt_flach", "leiter_schritt", "zweiter_weg", [(K.M_sun, A0), (M_MW, A0)],
     lambda M, a0: (K.v_flach(M, a0) / K.c)**2)


def _gezeit_eigen(E):
    """√(3/2)·|Ē^TF| über die Eigenwerte des spurfreien Teils."""
    lam = np.linalg.eigvalsh(np.asarray(E, float))
    lam = lam - lam.mean()
    return math.sqrt(1.5 * float(np.sum(lam**2)))


def _drehung(E, w):
    c_, s_ = math.cos(w), math.sin(w)
    R = np.array([[c_, -s_, 0], [s_, c_, 0], [0, 0, 1]])
    return R @ np.asarray(E) @ R.T


_neu("gezeit_betrag_eigenwerte", "gezeit_betrag", "zweiter_weg",
     [(np.diag([2e-30, -1e-30, -1e-30]),), (_drehung(np.diag([3e-30, -0.5e-30, 1e-30]), 0.7),)], _gezeit_eigen,
     toleranz=1e-12, grund="Eigenwerte statt Frobenius-Norm")
_neu("glaettungslaenge_dichte", "glaettungslaenge", "umkehr", [(K.M_sun, 5.553e-30), (K.M_earth, 5.553e-30)],
     lambda M, T: (3 * M / (4 * math.pi * K.gezeitendichte(T)))**(1 / 3))
_neu("gezeitendichte_wolke", "gezeitendichte", "umkehr", [(5.553e-30,), (1e-28,)],
     lambda T: K.M_sun / (4 * math.pi / 3 * K.glaettungslaenge(K.M_sun, T)**3))


def _schnittmasse(T, a0):
    from scipy.optimize import brentq
    u = brentq(lambda u: math.log(K.glaettungslaenge(math.exp(u), T) / K.r_stern(math.exp(u), a0)), 0, 250, xtol=1e-13)
    return math.exp(u)


_neu("schnittmasse_glaettung_gleichheit", "schnittmasse_glaettung", "zweiter_weg", [(5.553e-30, A0), (1e-28, 1.2e-10)],
     _schnittmasse, toleranz=1e-9, grund="Nullstelle ℓ(M) = r*(M) mit brentq")
_neu("schnittradius_glaettung_gleichheit", "schnittradius_glaettung", "zweiter_weg", [(5.553e-30, A0), (1e-28, 1.2e-10)],
     lambda T, a0: K.r_stern(_schnittmasse(T, a0), a0), toleranz=1e-9, grund="Nullstelle ℓ(M) = r*(M) mit brentq")
_neu("glaettung_zu_schwelle_direkt", "glaettung_zu_schwelle", "zweiter_weg",
     [(K.M_sun, 5.553e-30, A0), (1e4 * K.M_sun, 1e-28, 1.2e-10)],
     lambda M, T, a0: K.glaettungslaenge(M, T) / K.r_stern(M, a0))
_neu("glaettung_zu_schwelle_modelltext", "glaettung_zu_schwelle", "handwert",
     [(K.M_sun, A0 / (642 * K.pc), A0)], lambda M, T, a0: 26.2,
     toleranz=2e-3, grund="Sonne 26,2 mit R× = a₀/𝒯 = 642 pc (Stellen des Modelltexts)")


def _flussdimension(aN, a0):
    """1 − d ln a/d ln r an der Punktmasse (GM = 1, a₀ = 1 skaliert), Differenzenquotient."""
    r, h = math.sqrt(a0 / aN), 1e-5
    return 1 - (math.log(_a_punkt(r * (1 + h))) - math.log(_a_punkt(r * (1 - h)))) / (math.log(1 + h) - math.log(1 - h))


_neu("flussdimension_ableitung", "flussdimension", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     _flussdimension, toleranz=1e-8, grund="Differenzenquotient h = 10⁻⁵")


def _epizyklen(aN, a0):
    """Punktmasse mit Regel (GM = a₀ = 1): Ω² = a/r, κ² = r dΩ²/dr + 4Ω², ω_z² = ∂(a_z)/∂z bei z = 0."""
    r, h = math.sqrt(a0 / aN), 1e-5
    Om2 = lambda x: _a_punkt(x) / x
    kap2 = r * (Om2(r * (1 + h)) - Om2(r * (1 - h))) / (2 * h * r) + 4 * Om2(r)
    dz = h * r
    az = lambda z: _a_punkt(math.hypot(r, z)) * z / math.hypot(r, z)
    wz2 = (az(dz) - az(-dz)) / (2 * dz)
    return (kap2 + wz2) / Om2(r)


_neu("epizyklen_summe_bahn", "epizyklen_summe", "zweiter_weg", [(1e-12, A0), (1.1e-10, A0), (3e-9, 1.2e-10)],
     _epizyklen, toleranz=1e-8, grund="Differenzenquotienten h = 10⁻⁵r")

# ================================================================== Zustandsraum (Z)
_neu("realteil_quadrat_punkt_epizykel", "realteil_quadrat_punkt", "zweiter_weg", [(0.1,), (1.0,), (4.0,)],
     lambda s: (lambda r, h: (r * ((_a_punkt(r * (1 + h)) / (r * (1 + h)) - _a_punkt(r * (1 - h)) / (r * (1 - h)))
                                   / (2 * h * r)) + 4 * _a_punkt(r) / r) / (_a_punkt(r) / r))(1 / s, 1e-5),
     toleranz=1e-8, grund="(κ/Ω)² als Differenzenquotient, r* = 1")
_neu("realteil_quadrat_punkt_triade", "realteil_quadrat_punkt", "handwert", [(1.0,)], lambda s: 1.25777**2,
     toleranz=1e-5, grund="Re τ_M(r*) = 1,25777 (Rechnung triade_ebenen, 6 Stellen)")


def _resonanz(p_zu_k, rstern):
    from scipy.optimize import brentq
    s = brentq(lambda x: 1 + x * float(K.n_B(x)) - p_zu_k**2, 1e-6, 50, xtol=1e-15)
    return rstern / s


_neu("radius_resonanz_bose", "radius_resonanz", "zweiter_weg", [(1.2, 1.0), (1.25, 2.7e20), (4 / 3, 1.0)], _resonanz,
     toleranz=1e-9, grund="Nullstelle 1 + s·n_B(s) = (p/k)² mit brentq")
_neu("radius_resonanz_modelltext", "radius_resonanz", "handwert", [(1.2, 1.0), (1.25, 1.0), (9 / 8, 1.0)],
     lambda pk, r: {1.2: 0.68, 1.25: 0.95, 9 / 8: 0.45}[pk], toleranz=1.2e-2,
     grund="0,68 / 0,95 / 0,45 r* (halbe letzte Stelle)")
_neu("radius_resonanz_modelltext_aussen", "radius_resonanz", "handwert", [(4 / 3, 1.0), (7 / 5, 1.0)],
     lambda pk, r: {4 / 3: 2.1, 7 / 5: 12.0}[pk], toleranz=4.2e-2, grund="2,1 / 12 r* (halbe letzte Stelle)")


def _im_tau(q):
    from scipy.optimize import brentq
    return brentq(lambda y: math.log(q) + 2 * math.pi * y, 1e-12, 1e3, xtol=1e-15)


_neu("im_tau_aus_nome_umkehr", "im_tau_aus_nome", "umkehr", [(0.843,), (0.5,), (1e-3,)], _im_tau,
     toleranz=1e-12, grund="Nullstelle von |e^{2πiτ}| = q mit brentq")


def _kac_drehung(mu_X, mu_A, n_start=20000):
    """Goldene Drehung x → x + α (mod 1), Zelle [0, μ(A)/μ(X)): mittlere erste Rückkehrzeit über gleichmäßig verteilte
    Startpunkte der Zelle."""
    mu = mu_A / mu_X
    x0 = (np.arange(n_start) + 0.5) / n_start * mu
    x, T, offen = x0.copy(), np.zeros(n_start), np.ones(n_start, bool)
    for n in range(1, int(10 / mu) + 10):
        x = (x + GOLDEN) % 1.0
        zurueck = offen & (x < mu)
        T[zurueck] = n
        offen &= ~zurueck
        if not offen.any():
            break
    return float(T.mean())


_neu("kac_rueckkehrzeit_drehung", "kac_rueckkehrzeit", "zweiter_weg", [(1.0, 0.01), (100.0, 5.0)], _kac_drehung,
     toleranz=1e-4, grund="Mittel über 2·10⁴ Startpunkte (Drei-Abstände-Satz: Rückkehrzeiten 55, 89, 144)")


def _kac_markow(mu_X, mu_A_liste, n_schritte=400000):
    """Markow-Kette: träger Irrweg auf dem Ring mit μ(X) Zuständen (doppelt stochastisch), Zustände zufällig auf die
    Bereiche verteilt (Größen μ(A_i)); Rückkehrzeiten je Bereich aus einer Bahn gemessen (eigener Generator mit Seed)."""
    rng = np.random.default_rng(29092026)
    n = int(mu_X)
    bereich = rng.permutation(np.repeat(np.arange(len(mu_A_liste)), np.asarray(mu_A_liste, int)))
    schritt = rng.choice([-1, 0, 1], size=n_schritte, p=[0.3, 0.4, 0.3])
    zustand = np.cumsum(schritt) % n
    summe = 0.0
    for i in range(len(mu_A_liste)):
        t = np.nonzero(bereich[zustand] == i)[0]
        summe += 1 / float(np.mean(np.diff(t)))
    return summe


_neu("kac_budget_markow", "kac_budget", "zweiter_weg", [(12, [3, 4, 5]), (20, [1, 2, 7, 10])], _kac_markow,
     toleranz=1e-2, grund="Markow-Bahn mit 4·10⁵ Schritten, Seed fest")


def _kac_drehung_rueckwaerts(mu_X, mu_A_liste, n=200000):
    """Wicklung im umgekehrten Umlaufsinn (x → x − α): Zerlegung des Kreises in Zellen, Σ 1/⟨T_i⟩ aus einer Bahn."""
    grenzen = np.cumsum([0.0] + [m / mu_X for m in mu_A_liste])
    x = (0.1234 - GOLDEN * np.arange(n)) % 1.0
    zelle = np.searchsorted(grenzen, x, side="right") - 1
    return sum(1 / float(np.mean(np.diff(np.nonzero(zelle == i)[0]))) for i in range(len(mu_A_liste)))


_neu("kac_budget_wicklung_rueckwaerts", "kac_budget", "zweiter_weg", [(7, [1, 1, 1, 1, 1, 1, 1]), (10, [1, 2, 3, 4])],
     _kac_drehung_rueckwaerts, toleranz=1e-3, grund="Bahn mit 2·10⁵ Schritten")

_neu("bose_leiter_summe", "bose_leiter", "zweiter_weg", [(0.5, 1), (1.0, 2), (1.6, 2), (3.0, 3)],
     lambda s, d: s**d * sum(math.exp(-n * s) for n in range(1, int(45 / s) + 2)), toleranz=1e-12,
     grund="Reihe abgeschnitten bei e^{−45}")
_neu("bose_leiter_integral_quad", "bose_leiter_integral", "zweiter_weg", [(1,), (2,), (3,)],
     lambda d: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda s: s**d / math.expm1(s), 0, 200,
                                                                    epsabs=0, epsrel=1e-12)[0],
     toleranz=1e-10, grund="quad")


def _leiter_punktmasse(d):
    """∫(a − a_N)(r*/r)^{d−1} d ln r / a₀ an der Punktmasse (GM = a₀ = r* = 1), über u = ln r."""
    from scipy.integrate import quad
    f = lambda u: math.exp(-2 * u) / math.expm1(min(math.exp(-u), 700.0)) * math.exp(-(d - 1) * u)   # a − a_N = a_N/(e^s − 1)
    return quad(f, -12, 30, limit=400, epsabs=0, epsrel=1e-11)[0]


_neu("bose_leiter_integral_flussmass", "bose_leiter_integral", "zweiter_weg", [(1,), (2,), (3,)], _leiter_punktmasse,
     toleranz=1e-8, grund="Ortsintegral der Punktmasse mit quad, abgeschnitten bei r = e⁻¹² und e³⁰")

_neu("fisher_metrik_ableitung", "fisher_metrik", "zweiter_weg", [(0.2,), (1.0,), (4.0,)],
     lambda s: -(_nu(s * (1 + 1e-5)) - _nu(s * (1 - 1e-5))) / (2e-5 * s), toleranz=1e-8,
     grund="−dν/ds als Differenzenquotient")
_neu("fisher_metrik_sinh", "fisher_metrik", "zweiter_weg", [(0.2,), (1.0,), (4.0,)],
     lambda s: 1 / (4 * math.sinh(s / 2)**2))


def _entropie_shannon(s):
    q = math.exp(-s)
    n = np.arange(0, int(60 / s) + 2)
    p = (1 - q) * q**n
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


_neu("entropie_ruecklauf_shannon", "entropie_ruecklauf", "zweiter_weg", [(0.1,), (1.0,), (5.0,)], _entropie_shannon,
     toleranz=1e-10, grund="Summe −Σp_n ln p_n abgeschnitten bei qⁿ = e^{−60}")

_neu("schleifensumme_rueckläufe", "schleifensumme", "zweiter_weg", [(0.3,), (0.843,), (-0.05,), (0.99,)],
     lambda q: sum(q**n for n in range(0, int(45 / -math.log(abs(q))) + 2)), toleranz=1e-12,
     grund="Reihe Σ qⁿ abgeschnitten bei |q|ⁿ = e^{−45}")


def _falte(a, delta):
    """Normalform ẋ = −(x³ − a·x − c), c = 2(a/3)^{3/2} − δ: stabiler Ast vor der Falte, 1 − q = f'(x*)."""
    r = np.real(np.roots([1.0, 0.0, -a, -(2 * (a / 3)**1.5 - delta)]))
    x = min(r)
    return 1 - (3 * x * x - a)


_neu("kipp_schleife_normalform", "kipp_schleife", "grenzfall", [(1.0, 1e-6), (2.0, 1e-6), (0.5, 1e-7)], _falte,
     toleranz=1e-5, grund="Asymptote vor der Falte (δ → 0), Rest ∝ δ")


def _adler(k, phi):
    """ψ̇ = Δω − 2k cos φ sin ψ bei Δω = 2k + 1 (sicher außerhalb): Schwebung Ω_b aus der Zeit eines Umlaufs
    (solve_ivp), dann Rastbreite √(Δω² − Ω_b²)."""
    from scipy.integrate import solve_ivp
    dw = 2 * k + 1.0
    ereignis = lambda t, y: y[0] - 2 * math.pi
    ereignis.terminal = True
    sol = solve_ivp(lambda t, y: [dw - 2 * k * math.cos(phi) * math.sin(y[0])], (0, 1e4), [0.0], events=ereignis,
                    rtol=1e-12, atol=1e-12)
    Ob = 2 * math.pi / sol.t_events[0][0]
    return math.sqrt(dw * dw - Ob * Ob)


_neu("rastbreite_schwebung", "rastbreite", "zweiter_weg", [(0.5, 0.0), (0.5, math.pi / 3), (0.2, 2.0)], _adler,
     toleranz=1e-7, grund="Umlaufzeit aus solve_ivp (rtol 10⁻¹²)")


def _preisach_eingaben():
    """Hysteronen und Folgen mit eigenem Generator mit Seed."""
    rng = np.random.default_rng(1935)
    a, b = rng.uniform(-1, 1, 3000), rng.uniform(-1, 1, 3000)
    al, be, mass = np.maximum(a, b), np.minimum(a, b), rng.uniform(0, 1, 3000)
    folgen = [np.concatenate([[-1.0], np.cumsum(rng.normal(0, 0.2, 300)) * 0.3]), [-1, 0.8, -0.3, 0.6, 0.2],
              [-1, 0.8, -0.1, 0.4, -0.3, 0.6, 0.2]]
    return [(list(f), al, be, mass) for f in folgen]


def _preisach_voll(folge, alpha, beta, mass):
    z = -np.ones(alpha.shape)
    for u in folge:
        z = np.where(u >= alpha, 1.0, np.where(u <= beta, -1.0, z))
    return float(mass @ z)


_neu("preisach_ausgang_volle_folge", "preisach_ausgang", "zweiter_weg", _preisach_eingaben(), _preisach_voll)

def _rc_netz(nR, RG):
    """Leiter aus Längswiderstand R und Querleitwert G je Glied, Eingang V(0) = V_ein + V(n_R), Ende nach e^{−60}:
    Knotengleichungen als Netz gelöst, ν = V(0)/V_ein."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    n = int(60 / math.acosh(1 + RG / 2)) + nR + 10
    A = diags([np.ones(n), -(2 + RG) * np.ones(n + 1), np.ones(n)], [-1, 0, 1], format="lil")
    A[0, :] = 0
    A[0, 0], A[0, nR] = 1.0, -1.0
    A[n, :] = 0
    A[n, n] = 1.0
    b = np.zeros(n + 1)
    b[0] = 1.0
    return float(spsolve(A.tocsr(), b)[0])


_neu("rc_kette_nu_netz", "rc_kette_nu", "zweiter_weg", [(20, 0.01), (5, 1e-4), (3, 0.5)], _rc_netz,
     toleranz=1e-10, grund="Netz mit Ende nach e^{−60}")
_neu("rc_kette_nu_modelltext", "rc_kette_nu", "handwert", [(20, 0.01)], lambda nR, RG: 1.1567, toleranz=5e-5,
     grund="1,1567 (halbe letzte Stelle)")


def _gittersumme(k, tau, N=400):
    """Σ′(m + nτ)^{−k} über |m|, |n| ≤ N, dazu der Rest außerhalb als Integral 2π/((k − 2)N^{k−2}·Im τ) für k = 4
    nicht nötig: Toleranz nach N."""
    m, n = np.meshgrid(np.arange(-N, N + 1), np.arange(-N, N + 1))
    z = m + n * complex(tau)
    z = z[(m != 0) | (n != 0)]
    return float(np.sum(z**(-int(k))).real)


_neu("eisenstein_G_gittersumme", "eisenstein_G", "zweiter_weg", [(4, 1j), (4, 1.5j), (6, 1.5j), (8, 1j)], _gittersumme,
     toleranz=2e-5, grund="Gitter |m|, |n| ≤ 400, Rest ∝ N^{2−k}")
_neu("eisenstein_G_lemniskate", "eisenstein_G", "handwert", [(4, 1j)],
     lambda k, tau: math.gamma(0.25)**8 / (960 * math.pi**2))

_neu("nu_aus_modul_nome", "nu_aus_modul", "umkehr", [(0.5j,), (2j,), (0.3 + 0.2j,)],
     lambda tau: float(K.nu((2 * math.pi * complex(tau).imag)**2 * A0, A0)))
_neu("nu_aus_modul_modelltext", "nu_aus_modul", "handwert", [(0.5j,)], lambda tau: 1.045, toleranz=5e-4,
     grund="1,045 (halbe letzte Stelle)")
_neu("nu_aus_modul_modelltext_spitze", "nu_aus_modul", "handwert", [(2j,)], lambda tau: 1.0000035, toleranz=5e-8,
     grund="1,0000035 (halbe letzte Stelle)")


# ------------------------------------------------------------------ Kupferdraht: Impedanz und Routh-Hurwitz
def _kupfer_routh(Qm, V):
    """Heißer Ast über die Impedanz Z_m = Γ²/(c_m + i(mω − k(ΔT)/ω)); Falte = Ende der drei Lösungen; Hopf-Grenze über
    Routh-Hurwitz a₁a₂ − a₃ = 0 des charakteristischen Polynoms der Hüllkurve (Jacobi-Matrix mit Schritt 10⁻⁵).
    Gibt (Abstand vor der Falte [Hz], Ringfrequenz √a₂/2π [Hz])."""
    from scipy.optimize import brentq
    p = K.draht(K.CU, dict(K.DRAHT, Qm=Qm))
    R0, alR, Gam, m, cm, kk, Gth, Cth, w0 = (p[x] for x in ('R0', 'alR', 'Gam', 'm', 'c_m', 'k', 'G_th', 'C_th', 'w0'))

    def f_imp(dT, wd):
        Rt = R0 * (1 + alR * dT)
        Zm = Gam**2 / (cm + 1j * (m * wd - kk(dT) / wd))
        return 0.5 * np.abs(V / (Rt + Zm))**2 * (Rt + Zm.real) - Gth * dT

    def loesungen(d):
        wd = w0 + 2 * math.pi * d
        xs = np.linspace(0.0, 3.0, 30001)
        fv = f_imp(xs, wd)
        idx = np.nonzero(np.sign(fv[:-1]) != np.sign(fv[1:]))[0]
        return [brentq(f_imp, xs[i], xs[i + 1], args=(wd,), xtol=1e-14) for i in idx]

    def hurwitz(d):
        wd = w0 + 2 * math.pi * d
        dT = loesungen(d)[-1]
        Rt = R0 * (1 + alR * dT)
        X = Gam * V / Rt / (kk(dT) - m * wd**2 + 1j * wd * cm + 1j * wd * Gam**2 / Rt)

        def F(y):
            Xc, T = y[0] + 1j * y[1], y[2]
            R = R0 * (1 + alR * T)
            I = (V - 1j * wd * Gam * Xc) / R
            dX = (Gam * I - (kk(T) - m * wd**2 + 1j * wd * cm) * Xc) / (2j * wd * m)
            return np.array([dX.real, dX.imag, (0.5 * R * abs(I)**2 + 0.5 * cm * wd**2 * abs(Xc)**2 - Gth * T) / Cth])
        y = np.array([X.real, X.imag, dT])
        J = np.column_stack([(F(y + e) - F(y - e)) / (2 * e.max()) for e in
                             np.diag([1e-5 * abs(X), 1e-5 * abs(X), 1e-5])])
        a = np.poly(J)
        return a[1] * a[2] - a[3], a[2]

    d = np.arange(0.0, -400.0, -1.0)
    ja = [x for x in d if len(loesungen(x)) >= 3]
    lo, hi = ja[-1], ja[-1] - 1.0
    for _ in range(34):
        mi = 0.5 * (lo + hi)
        lo, hi = (mi, hi) if len(loesungen(mi)) >= 3 else (lo, mi)
    falte = lo
    abst = np.geomspace(1e-3, 5.0, 40)
    h = [hurwitz(falte + x)[0] for x in abst]
    i = max(j for j in range(len(abst)) if h[j] < 0)
    a, b = abst[i], abst[i + 1]
    for _ in range(30):
        mi = 0.5 * (a + b)
        a, b = (mi, b) if hurwitz(falte + mi)[0] < 0 else (a, mi)
    return a, math.sqrt(hurwitz(falte + a)[1]) / (2 * math.pi)


def _v_kupfer(Qm, Vf):
    """Antrieb mit gleicher Schleifenverstärkung wie Q_m = 1000 bei 10 mV (belastete Güte 851)."""
    p = K.draht(K.CU, dict(K.DRAHT, Qm=Qm))
    Ql = p['m'] * p['w0'] / p['c_m'] / (1 + p['Gam']**2 / p['c_m'] / p['R0'])
    return 10e-3 * Vf * math.sqrt(851.0 / Ql)


_KUPFER_FAELLE = [(1e4, _v_kupfer(1e4, 2.0)), (3e3, _v_kupfer(3e3, 2.0))]
_neu("kupfer_hopf_abstand_routh", "kupfer_hopf_abstand", "zweiter_weg", _KUPFER_FAELLE,
     lambda Qm, V: _kupfer_routh(Qm, V)[0], toleranz=2e-3,
     grund="Falte auf dem Gitter der stationären Punkte (10⁻⁴ K), Routh-Hurwitz statt Eigenwerte")
_neu("kupfer_ringfrequenz_routh", "kupfer_ringfrequenz", "zweiter_weg", _KUPFER_FAELLE,
     lambda Qm, V: _kupfer_routh(Qm, V)[1], toleranz=2e-3,
     grund="Ringfrequenz √a₂/2π an der Routh-Hurwitz-Grenze")
_neu("kupfer_ringfrequenz_alt", "kupfer_ringfrequenz", "handwert", [_KUPFER_FAELLE[0]], lambda Qm, V: 1.323,
     toleranz=1e-3, grund="1,323 Hz bei Q_m = 10⁴, V = 9,676 mV (Rechnung kupfer_selbstpulsieren, 4 Stellen)")
_neu("kupfer_zeitverhaeltnis_mode", "kupfer_zeitverhaeltnis", "zweiter_weg", [(1e3,), (3e3,), (1e4,)],
     lambda Qm: (lambda p: (p['C_th'] / p['G_th']) / (2 * p['m'] / (p['c_m'] + p['Gam']**2 / p['R0'])))(
         K.draht(K.CU, dict(K.DRAHT, Qm=Qm))))
_neu("kupfer_zeitverhaeltnis_modelltext", "kupfer_zeitverhaeltnis", "handwert", [(3e3,)], lambda Qm: 0.74,
     toleranz=7e-3, grund="0,74 bei Güte 3000 (halbe letzte Stelle)")


def _bose_stufe_lambert(d):
    """Maximum von s^d/(e^s − 1) geschlossen: d(1 − e^{−s}) = s ⇒ s = d + W₀(−d·e^{−d}); gibt (s, s², Wert)."""
    from scipy.special import lambertw
    s = d + float(lambertw(-d * math.exp(-d)).real)
    return s, s * s, s**d / math.expm1(s)


_neu("bose_stufe_max_lambert", "bose_stufe_max", "zweiter_weg", [(2,), (3,), (4,)], _bose_stufe_lambert,
     toleranz=1e-10, grund="Nullstelle mit brentq (xtol 2·10⁻¹²) gegen Lambert-W")


# ------------------------------------------------------------------ Kupferdraht: Lyapunov-Koeffizient und Zerlegung der Schleife
def _draht_symbolisch(Qm, V, dfd):
    """Fluss der Hüllkurve (Re Ẋ, Im Ẋ, ΔṪ) und Wärmeleistung P(ΔT_R, ΔT_E, ΔT_N) je Kanal als sympy-Ausdrücke in
    den Stoffwerten des Drahts (Hüllkurve und Impedanz der Mode, ohne Kernfunktion des Flusses)."""
    import sympy as sp
    d = dict(K.DRAHT, Qm=Qm)
    p = K.draht(K.CU, d)
    cu = K.CU
    A, I, ell = math.pi * d['a']**2, math.pi * d['a']**4 / 4, d['ell']
    x, y, T, tR, tE, tN = sp.symbols('x y T tR tE tN', real=True)
    R0, alR, Gam, m, cm, Gth, Cth = (sp.Float(p[k], 30) for k in ('R0', 'alR', 'Gam', 'm', 'c_m', 'G_th', 'C_th'))
    kf = lambda tE_, tN_: (8 * sp.pi**4 * cu['E'] * (1 + cu['alE'] * tE_) * I / ell**3
                           + 2 * sp.pi**2 * (d['sig0'] * A - cu['E'] * A * cu['alth'] * tN_) / ell)
    wd = p['w0'] + 2 * sp.pi * dfd
    Rt = R0 * (1 + alR * T)
    Ic = sp.expand_complex((V - sp.I * wd * Gam * (x + sp.I * y)) / Rt)
    dX = sp.expand_complex((Gam * Ic - (kf(T, T) - m * wd**2 + sp.I * wd * cm) * (x + sp.I * y)) / (2 * sp.I * wd * m))
    fluss = [sp.re(dX), sp.im(dX), (Rt * (sp.re(Ic)**2 + sp.im(Ic)**2) / 2 + cm * wd**2 * (x**2 + y**2) / 2
                                     - Gth * T) / Cth]
    R = R0 * (1 + alR * tR)
    Zm = sp.expand_complex(Gam**2 / (cm + sp.I * (m * wd - kf(tE, tN) / wd)))
    Zr, Zi = sp.re(Zm), sp.im(Zm)
    P = V**2 / 2 * (R + Zr) / ((R + Zr)**2 + Zi**2)
    return (x, y, T, tR, tE, tN), fluss, P, Gth


def _kupfer_lyapunov_symbolisch(Qm, V):
    """ℓ₁ mit exakten Ableitungen (sympy) bis zur dritten Ordnung; Gleichgewicht aus der Impedanzform, Hopf-Grenze aus
    kupfer_hopf."""
    import sympy as sp
    abst, falte, _ = K.kupfer_hopf(Qm, V)
    dfd = falte + abst
    (x, y, T, *_), fl, _, _ = _draht_symbolisch(Qm, V, dfd)
    p = K.draht(K.CU, dict(K.DRAHT, Qm=Qm))
    wd = p['w0'] + 2 * math.pi * dfd
    T0 = K.kupfer_huellkurve(dfd, V, Qm)['heiss']
    Rt = p['R0'] * (1 + p['alR'] * T0)
    X0 = p['Gam'] * V / Rt / (p['k'](T0) - p['m'] * wd**2 + 1j * wd * p['c_m'] + 1j * wd * p['Gam']**2 / Rt)
    sk = [abs(X0), abs(X0), 1.0]
    u = sp.symbols('u0:3', real=True)
    F = [f.subs({x: u[0] * sk[0], y: u[1] * sk[1], T: u[2] * sk[2]}) / sk[i] for i, f in enumerate(fl)]
    pt = {u[0]: X0.real / sk[0], u[1]: X0.imag / sk[1], u[2]: T0}
    J = np.array([[float(sp.diff(F[i], u[j]).subs(pt)) for j in range(3)] for i in range(3)])
    B2 = np.array([[[float(sp.diff(F[i], u[j], u[k]).subs(pt)) for k in range(3)] for j in range(3)] for i in range(3)])
    C3 = np.array([[[[float(sp.diff(F[i], u[j], u[k], u[l]).subs(pt)) for l in range(3)] for k in range(3)]
                    for j in range(3)] for i in range(3)])
    ev, vr = np.linalg.eig(J)
    i = int(np.argmax(ev.imag))
    om, q = ev[i].imag, vr[:, i]
    evl, vl = np.linalg.eig(J.T)
    pv = vl[:, int(np.argmin(np.abs(evl - np.conj(ev[i]))))]
    pv = pv / np.conj(np.vdot(pv, q))
    B = lambda a, b: np.einsum('ijk,j,k->i', B2, a, b)
    C = lambda a, b, c: np.einsum('ijkl,j,k,l->i', C3, a, b, c)
    h11 = np.linalg.solve(J, B(q, np.conj(q)))
    h20 = np.linalg.solve(2j * om * np.eye(3) - J, B(q, q))
    return float((np.vdot(pv, C(q, q, np.conj(q))) - 2 * np.vdot(pv, B(q, h11))
                  + np.vdot(pv, B(np.conj(q), h20))).real / (2 * om))


_neu("kupfer_lyapunov_symbolisch", "kupfer_lyapunov", "zweiter_weg", [(3e3, 6.577e-3), (1e4, _v_kupfer(1e4, 2.0))],
     _kupfer_lyapunov_symbolisch, toleranz=1e-4,
     grund="exakte Ableitungen (sympy) gegen zentrale Differenzen dritter Ordnung (Schritt 2·10⁻³)")


def _draht_schleife_symbolisch(dT, dfd, V, Qm):
    """q und seine drei Kanäle als exakte Ableitungen der Wärmeleistung der Impedanzform (sympy)."""
    import sympy as sp
    (_, _, _, tR, tE, tN), _, P, Gth = _draht_symbolisch(Qm, V, dfd)
    wert = {tR: dT, tE: dT, tN: dT}
    teile = [float(sp.diff(P, t).subs(wert)) / float(Gth) for t in (tR, tE, tN)]
    return (sum(teile), *teile)


_neu("draht_schleife_symbolisch", "draht_schleife", "zweiter_weg",
     [(0.1450, -10.0, 5e-3, 1e3), (0.0222, -2.0, 2e-3, 1e3), (0.6060, -40.0, 10e-3, 1e3)], _draht_schleife_symbolisch,
     toleranz=1e-6, grund="exakte Ableitung (sympy) gegen zentrale Differenz (Schritt 10⁻⁶ K); q als Summe der Kanäle")


# ------------------------------------------------------------------ Diskriminante, Resolvente, Spielkern, Rastpunkt (05, 06)
def _sigma(n, k):
    return sum(d**k for d in range(1, n + 1) if n % d == 0)


def _delta_eisenstein(n_max):
    """c_Δ(n) aus (E₄³ − E₆²)/1728 mit E₄ = 1 + 240Σσ₃(n)qⁿ, E₆ = 1 − 504Σσ₅(n)qⁿ (ganzzahlige Reihen)."""
    n_max = int(n_max)
    E4 = [1] + [240 * _sigma(n, 3) for n in range(1, n_max + 1)]
    E6 = [1] + [-504 * _sigma(n, 5) for n in range(1, n_max + 1)]
    mal = lambda a, b: [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(n_max + 1)]
    d = [x - y for x, y in zip(mal(mal(E4, E4), E4), mal(E6, E6))]
    return np.array([x // 1728 for x in d[1:]], dtype=float)


_neu("delta_reihe_eisenstein", "delta_reihe", "zweiter_weg", [(8,), (20,)], _delta_eisenstein)


def _delta_gespiegelt(tau):
    """Δ(it) = t^{−12}·Δ(i/t) (Gewicht 12), Δ(i/t) als q-Reihe aus E₄ und E₆ im Nome e^{−2π/t}."""
    t = complex(tau).imag
    qq = math.exp(-2 * math.pi / t)
    c = _delta_eisenstein(20)
    return t**-12 * sum(cn * qq**(n + 1) for n, cn in enumerate(c))


_neu("delta_modul_gespiegelt", "delta_modul", "zweiter_weg", [(1j / (2 * math.pi),), (0.5j,), (0.8j,)],
     _delta_gespiegelt, toleranz=1e-12, grund="Produkt bis 60/(2π Im τ) Glieder gegen 20 Glieder der gespiegelten Reihe")


def _resolvente_fibonacci(q, K_M):
    """Nächste Rückkehr der goldenen Drehung nur an Fibonacci-Zahlen k (beste Näherungen);
    |1 − q·e^{iθ}|⁻¹ = (1 − 2q cos θ + q²)^{−1/2}."""
    f, fib = [1, 2], []
    while f[-1] <= K_M:
        f.append(f[-1] + f[-2])
    fib = [k for k in f if k <= K_M]
    th = [2 * math.pi * ((k * GOLDEN) % 1.0) for k in fib]
    return max((1 - 2 * q * math.cos(x) + q * q)**-0.5 for x in th)


_neu("resolvente_flucht_fibonacci", "resolvente_flucht", "zweiter_weg", [(0.9999, 2000), (0.99, 2000), (0.9, 300)],
     _resolvente_fibonacci, toleranz=1e-9)


def _spiel_bahn(v, gamma):
    """Bahn ẍ = −1/x² − γ_K·y, ẏ = ẋ (unendliches Gedächtnis), Umkehr bei ẋ = 0."""
    from scipy.integrate import solve_ivp
    ereignis = lambda t, z: z[1]
    ereignis.terminal, ereignis.direction = True, -1
    s = solve_ivp(lambda t, z: [z[1], -1 / z[0]**2 - gamma * z[2], z[1]], (0, 1e4), [1.0, v, 0.0], rtol=1e-12,
                  atol=1e-13, events=ereignis)
    return s.y_events[0][0][0]


_neu("spiel_umkehrpunkt_bahn", "spiel_umkehrpunkt", "zweiter_weg", [(1.6, 0.01), (1.2, 0.05), (2.0, 0.001)],
     _spiel_bahn, toleranz=1e-8, grund="Bahnintegration (rtol 10⁻¹²) gegen Energieerhaltung")
_neu("spiel_umkehrpunkt_modelltext", "spiel_umkehrpunkt", "handwert", [(1.6, 0.01)], lambda v, g: 9.75,
     toleranz=5.2e-4, grund="9,75 im Modelltext (halbe letzte Stelle)")


def _rastpunkt_ableitung(k, phi, domega):
    """d ψ̇/dψ am stabilen Nullpunkt ψ* = arcsin(Δω/(2k cos φ)) von ψ̇ = Δω − 2k cos φ sin ψ, Differenzenquotient."""
    b = 2 * k * math.cos(phi)
    psi = math.asin(domega / b)
    f = lambda x: domega - b * math.sin(x)
    h = 1e-6
    return (f(psi + h) - f(psi - h)) / (2 * h)


_neu("rate_rastpunkt_ableitung", "rate_rastpunkt", "zweiter_weg", [(0.05, 0.0, 0.02), (0.3, 0.5, 0.1)],
     _rastpunkt_ableitung, toleranz=1e-8, grund="Differenzenquotient (Schritt 10⁻⁶)")


def _vo2_hurwitz(aufbau):
    """Hopf-Grenze des VO₂-Elements über Hurwitz: Fixpunkt als Nullstelle des vollen Flusses (fsolve), Jacobi-Matrix
    exakt (sympy), charakteristisches Polynom λ⁴ + a₁λ³ + a₂λ² + a₃λ + a₄; das elektrothermische Paar kreuzt die
    Achse, wo a₁a₂a₃ − a₃² − a₁²a₄ das Vorzeichen wechselt; dort ω² = a₃/a₁."""
    import sympy as sp
    from scipy.optimize import brentq, fsolve
    a = dict(aufbau)
    V, T, x, v, Vs = sp.symbols('V T x v Vs', real=True)
    s = lambda sh: 1 / (1 + sp.exp((T - a['T_C'] - sh) / a['breite']))
    R = sp.exp(sp.log(a['R_M']) + sp.log(a['R_I'] / a['R_M']) * s(a['kappa'] * x))
    I = (V - a['Gamma'] * a['B'] * v) / R
    om = 2 * sp.pi * a['f_m']
    y = x - (a['D_tr'] * (1 - s(0)) + a['lambda_B'] * a['B'] * I)
    F = sp.Matrix([((Vs - V) / a['R_L'] - I) / a['C'], (I * I * R - a['G_th'] * (T - a['T_0'])) / a['C_th'], v,
                   -om / a['Q'] * v - om**2 * y * (1 + (y / a['x_D'])**2)])
    Fn = sp.lambdify(((V, T, x, v), Vs), list(F), 'numpy')
    Jn = sp.lambdify(((V, T, x, v), Vs), F.jacobian([V, T, x, v]), 'numpy')
    start = [np.array([0.85, 340.2, 0.124, 0.0])]

    def poly(u):
        z = fsolve(lambda w: Fn(w, u), start[0], xtol=1e-12)
        start[0] = z
        c = np.poly(np.array(Jn(z, u), dtype=float))
        return c

    def hurwitz(u):
        c = poly(u)
        return c[1] * c[2] * c[3] - c[3]**2 - c[1]**2 * c[4]
    VH = brentq(hurwitz, *a['V_suche'], xtol=1e-12)
    c = poly(VH)
    return VH, math.sqrt(c[3] / c[1]) / (2 * math.pi) * 1e6


_neu("vo2_hopf_hurwitz", "vo2_hopf", "zweiter_weg", [(K.eingabe.EINGABEN["VO2_AUFBAU"].wert,)], _vo2_hurwitz,
     toleranz=1e-6, grund="Fixpunkt auf dem Temperaturgitter (70 mK, linear eingeschnürt) und Differenzen gegen fsolve "
                          "und exakte Jacobi-Matrix; Bisektion gegen brentq")


def _selbsterregt_bahn(g_E, b_I, dR_dT, G_th):
    """Grenzzyklus q̈ + (R − g_E)q̇ + b_I·q̇³ + q = 0 (L = C = 1) gerechnet, P = ⟨R·I² + b_I·I⁴⟩ über 200 Zeiteinheiten
    nach dem Einschwingen, bei R = 0,1 ± 10⁻⁴; Differenzenquotient nach T."""
    from scipy.integrate import solve_ivp
    R0, h = 0.1, 1e-4

    def P(R):
        s = solve_ivp(lambda t, z: [z[1], -z[0] - (R - g_E) * z[1] - b_I * z[1]**3], (0, 400), [0.5, 0.0],
                      rtol=1e-10, atol=1e-12, dense_output=True)
        tt = np.linspace(200, 400, 100001)
        I = s.sol(tt)[1]
        return np.trapezoid(R * I**2 + b_I * I**4, tt) / 200
    return (P(R0 + h) - P(R0 - h)) / (2 * h) * dR_dT / G_th


_neu("schleife_selbsterregt_bahn", "schleife_selbsterregt", "zweiter_weg", [(0.30, 0.2, 0.05, 1.0)],
     _selbsterregt_bahn, toleranz=0.03, grund="Mittelung über den Umlauf (erste Ordnung in g_E − R = 0,2) gegen den "
                                                "gerechneten Grenzzyklus")


def _kuramoto_klein(N_liste):
    """N = 2: ψ̇ = 1 − 2k sin ψ rastet ab k = ½; N = 3 (ω = −½, 0, ½, symmetrisch): ½ = k·max(sin φ + sin 2φ).
    Auf das Gitter k·N = 0,30 … 1,20 (Schritt 0,01) aufgerundet, wie es der Lauf liest."""
    from scipy.optimize import minimize_scalar
    aus = []
    for n in N_liste:
        if n == 2:
            kN = 2 * 0.5
        else:
            m = -minimize_scalar(lambda x: -(math.sin(x) + math.sin(2 * x)), bounds=(0, math.pi / 2),
                                 method='bounded', options=dict(xatol=1e-12)).fun
            kN = 3 * 0.5 / m
        kN = math.ceil(round(kN / 0.01, 9)) * 0.01
        aus.append((n - 1) * kN / n)
    return np.array(aus)


_neu("kuramoto_einrasten_klein", "kuramoto_einrasten", "zweiter_weg", [((2, 3),)], _kuramoto_klein, toleranz=1e-9)


def _rate_quer_fluss(D):
    """Kontraktion ρ̇ = −D·ρ der Querrichtung, integriert: ln(ρ(T)/ρ(0))/T."""
    from scipy.integrate import solve_ivp
    s = solve_ivp(lambda t, r: -D * r, (0, 10.0), [1.0], rtol=1e-12, atol=1e-15)
    return math.log(s.y[0, -1]) / 10.0


_neu("rate_quer_fluss", "rate_quer", "zweiter_weg", [(0.5,), (0.1,)], _rate_quer_fluss, toleranz=1e-8,
     grund="Integration rtol 10⁻¹²")


# ------------------------------------------------------------------ Gitter-Rückkehr, Drehung, Windungskorrelation (04, 07)
def _watson_gamma(d):
    """Watson (1939): G₃(1) = √6/(32π³)·Γ(1/24)Γ(5/24)Γ(7/24)Γ(11/24)."""
    from scipy.special import gamma
    return math.sqrt(6) / (32 * math.pi**3) * gamma(1 / 24) * gamma(5 / 24) * gamma(7 / 24) * gamma(11 / 24)


_neu("watson_rueckkehr_gamma", "watson_rueckkehr", "zweiter_weg", [(3,)], _watson_gamma, toleranz=1e-9)


def _drehung_je_start(alpha, mu):
    """Erste Rückkehrzeit je Startpunkt x ∈ [0, μ) (Gitter aus 20000 Punkten): kleinstes n ≥ 1 mit
    (x + nα) mod 1 < μ; gibt die verschiedenen Werte."""
    x = (np.arange(20000) + 0.5) / 20000 * mu
    t = np.zeros_like(x)
    offen = np.ones(x.size, bool)
    n = 0
    while offen.any():
        n += 1
        treffer = offen & (np.mod(x + n * alpha, 1.0) < mu)
        t[treffer] = n
        offen &= ~treffer
    return np.unique(t)


_neu("rueckkehrzeiten_drehung_je_start", "rueckkehrzeiten_drehung", "zweiter_weg",
     [(GOLDEN, 0.01), (GOLDEN, 0.001), (math.sqrt(2) - 1, 0.01)], _drehung_je_start)


def _chsh_frei(N_H):
    """CHSH mit freien Winkeln a = 0, b, a′, b′ (Nelder-Mead aus mehreren Starts), E aus den Harmonischen."""
    from scipy.optimize import minimize
    k = 2 * np.arange(int(N_H)) + 1.0
    w = (1 / k**2) / np.sum(1 / k**2)
    E = lambda d: -np.dot(w, np.cos(k * d))
    S = lambda x: -abs(E(-x[0]) - E(-x[2]) + E(x[1] - x[0]) + E(x[1] - x[2]))
    beste = 0.0
    for start in ((0.7, 1.6, -0.7), (0.5, 1.2, -0.5), (0.3, 1.0, -0.3), (1.0, 2.0, -1.0)):
        r = minimize(S, start, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000))
        beste = max(beste, -r.fun)
    return beste


_neu("chsh_oberschwingungen_frei", "chsh_oberschwingungen", "zweiter_weg", [(1,), (2,), (3,)], _chsh_frei,
     toleranz=1e-7, grund="freie Winkel (Nelder-Mead, xatol 10⁻¹²) gegen den drehinvarianten Ansatz")
_neu("chsh_oberschwingungen_grenzen", "chsh_oberschwingungen", "grenzfall", [(1,), (10000,)],
     lambda N: 2 * math.sqrt(2) if N == 1 else 2.0, toleranz=1e-4,
     grund="N_H = 1: 2√2 (Tsirelson); N_H → ∞: Dreieck, 2 (Rest ∝ 1/N_H)")


# ------------------------------------------------------------------ Zwei Kanäle, Rückgabe (04, 05)
def _kanal_ableitungen(M, r, f_g):
    """Differenzenquotienten: ½f_g·|d ln(GM/r²)/dr| gegen d ln√(1 − 2GM/(rc²))/dr."""
    GM = K.G * M
    h = 1e-5 * r
    dg = (math.log(GM / (r + h)**2) - math.log(GM / (r - h)**2)) / (2 * h)
    takt = lambda x: 0.5 * math.log1p(-2 * GM / (x * K.c**2))
    dt = (takt(r + h) - takt(r - h)) / (2 * h)
    return 0.5 * f_g * abs(dg) / dt


_neu("kanal_verhaeltnis_ableitungen", "kanal_verhaeltnis", "zweiter_weg",
     [(K.M_earth, K.R_earth, 1.0), (K.M_sun, K.R_sun, 0.5), (1.4 * K.M_sun, 1.2e4, 1.0)], _kanal_ableitungen,
     toleranz=1e-8, grund="Differenzenquotient mit Schritt 10⁻⁵r")
_neu("rueckgabe_rate_lebensdauer", "rueckgabe_rate", "zweiter_weg", [(K.M_sun,), (1e12,), (1e40,)],
     lambda M: -M / (3 * K.lebensdauer_falte(M)))
