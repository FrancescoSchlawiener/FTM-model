# -*- coding: utf-8 -*-
"""
Gegenproben der Gruppe Galaxie und Daten (Felder in gegenprobe.py): je Formel des Modells mindestens ein Weg, der die
geprüfte Funktion nicht aufruft. Wo die Eingaben einer Formel schon abgeleitete Größen sind (Frequenzen, Bin-Mediane),
rechnet der zweite Weg beim Laden aus dem Feld selbst und vergleicht mit dem Wert, den die Formel aus den Größen gibt.
"""
import math

import numpy as np

import kern as K

from .gegenprobe import GEGENPROBEN

_n = GEGENPROBEN.neu
A0, kpc, G = K.A0, K.kpc, K.G
MW = 6e10 * K.M_sun


def _nu(x, a0=A0):
    """ν = 1/(1 − e^{−√(x/a₀)}) von Hand, ohne den Kern."""
    return 1.0 / (1.0 - math.exp(-math.sqrt(x / a0)))


def _chi(x, a0=A0):
    s = math.sqrt(x / a0)
    return s / math.expm1(s)


# ------------------------------------------------------------------ Lesung an Rotationskurven
_n("v_bar_quadrat_handwert", "v_bar_quadrat", "handwert", [(-10e3, 100e3, 50e3, 0.5, 0.7)],
   lambda *a: -1e8 + 5e9 + 1.75e9)
_n("g_gemessen_ohne_druck", "g_gemessen", "grenzfall", [(3 * kpc, 150e3, 2 * kpc, 0.0)], lambda R, V, Rd, s: V * V / R)
_n("g_gemessen_handwert", "g_gemessen", "handwert", [(2 * kpc, 100e3, 2 * kpc, 10e3)], lambda R, V, Rd, s: 1.01e10 / R)

_R = np.linspace(1, 20, 15) * kpc
_Vd = np.full(_R.size, 100e3)
_gN = 0.5 * _Vd**2 / _R
_Vo = np.sqrt(np.array([_nu(x) for x in _gN]) * _gN * _R) * 10**0.05            # Rest genau 0,1 dex
_ge = 0.03 * A0
_Vo_e = np.sqrt(np.array([_nu(math.hypot(x, _ge)) for x in _gN]) * _gN * _R)    # mit äußerem Feld senkrecht, Rest 0
_Vo_p = np.sqrt(np.array([_nu(x) for x in _gN]) * _gN * _R * 10**0.1 - (8e3)**2 * _R / (3 * kpc))  # Gasdruck σ = 8 km/s
_Z = np.zeros(_R.size)
_n("rest_galaxie_synthetisch", "rest_galaxie", "grenzfall",
   [(_R, _Vo, _Z, _Vd, _Z, 3 * kpc, 0.5), (_R, _Vo_e * 10**0.05, _Z, _Vd, _Z, 3 * kpc, 0.5, 0.0, _ge),
    (_R, _Vo_p, _Z, _Vd, _Z, 3 * kpc, 0.5, 8e3)],
   lambda *a: 0.1, toleranz=1e-9)
_Vo13 = np.sqrt(np.array([_nu(x, 1.3e-10) for x in _gN]) * _gN * _R)
_n("a0_je_galaxie_synthetisch", "a0_je_galaxie", "grenzfall", [(_R, _Vo13, _Z, _Vd, _Z, 3 * kpc, 0.5)],
   lambda *a: 1.3e-10, toleranz=1e-8, grund="Nullstelle mit brentq, xtol 2·10⁻¹² in log a₀")

_GAL = dict(D=10 * K.Mpc, eD=0.1 * K.Mpc, inc=60.0, eInc=0.0, R=_R, Vo=_Vo, Vg=_Z, Vd=_Vd, Vb=_Z, Rd=3 * kpc)
_n("streuung_messfehler_abstand", "streuung_messfehler", "handwert",
   [(_GAL, lambda R, Vo, Vg, Vd, Vb, Rd: math.log10(Rd / (3 * kpc)), 4000, 7, 0.0)],
   lambda *a: 0.01 / math.log(10), toleranz=0.05,
   grund="Streuung von log₁₀(1 + δ), δ ~ N(0; 0,01): 0,01/ln 10; Stichprobe 4000 Ziehungen (≈ 1 %)")

_Rf = np.linspace(1, 40, 40001) * kpc
_n("toomre_flache_kurve", "toomre_sigma_dyn", "handwert", [(_Rf, np.full(_Rf.size, 200e3))],
   lambda R, V: (200e3)**2 / (math.pi * G * kpc), toleranz=1e-6, grund="Trapez auf 40001 Punkten")


# ------------------------------------------------------------------ Scheiben
_n("g_expscheibe_fernfeld", "g_expscheibe", "grenzfall", [(MW, 3 * kpc, 600 * kpc)], lambda M, Rd, R: G * M / R**2,
   toleranz=2e-4, grund="R = 200 R_d, Korrektur ∝ (R_d/R)², hier 1,1·10⁻⁴")
_n("sigma_dyn_expscheibe_newton", "sigma_dyn_expscheibe", "grenzfall", [(100 * K.M_sun / K.pc**2, 1e-40)],
   lambda S0, a0: S0, toleranz=1e-6, grund="Toomre: Newton gibt Σ_dyn(0) = Σ₀; Quadratur")


def _dphi(phi, R, z, k, h):
    return -(phi(R + h, z) - phi(R - h, z)) / (2 * h) if k == "R" else -(phi(R, z + h) - phi(R, z - h)) / (2 * h)


_n("kuzmin_feld_potential", "kuzmin_feld", "zweiter_weg",
   [(MW, 3 * kpc, 8 * kpc, "R"), (MW, 3 * kpc, 8 * kpc, "z")],
   lambda M, b, R, k: _dphi(lambda r, z: G * M / math.hypot(r, b + z), R, 1e-7 * kpc, k, 1e-8 * kpc),
   toleranz=1e-6, grund="Differenzen des Potentials −GM/√(R² + (b + z)²) knapp über der Ebene")
_n("g_miyamoto_nagai_potential", "g_miyamoto_nagai", "zweiter_weg",
   [(MW, 3 * kpc, 0.3 * kpc, 12 * kpc, 0.0, "R"), (MW, 3 * kpc, 0.3 * kpc, 12 * kpc, 0.5 * kpc, "z")],
   lambda M, a, b, R, z, k: _dphi(lambda r, zz: G * M / math.hypot(r, a + math.hypot(zz, b)), R, z, k, 1e-4 * kpc),
   toleranz=1e-6, grund="zentrale Differenzen des Potentials")
_n("a_scheibe_ebene_radial", "a_scheibe_ebene", "grenzfall", [(3e-11, 50 * K.M_sun / K.pc**2, 0.0)],
   lambda g, S, f: _nu(g) * g, toleranz=1e-12)
_gk, _Sk = 1e-16, 1e-16 / (2 * math.pi * G) * 0.75            # tief, g_z/g_R = 0,75
_n("a_scheibe_ebene_kuzmin_tief", "a_scheibe_ebene", "grenzfall", [(_gk, _Sk, 1.0)],
   lambda g, S, f: math.sqrt(A0 * g) * math.sqrt(g / math.hypot(g, 2 * math.pi * G * S)), toleranz=1e-3,
   grund="Tiefbereich: a_Feld/a_radial → √(cos θ), Korrektur ∝ s")

_h = 1e-5 * kpc
_Rv = 12 * kpc
_gRv = K.g_miyamoto_nagai(MW, 3 * kpc, 0.3 * kpc, _Rv, 0.0, "R")
_gzv = K.g_miyamoto_nagai(MW, 3 * kpc, 0.3 * kpc, _Rv, _h, "z")
_wzN, _OmN, _nuv = math.sqrt(_gzv / _h), math.sqrt(_gRv / _Rv), _nu(_gRv)
_az_ort = _nu(math.hypot(_gRv, _gzv)) * _gzv
_az_kern = _gzv + (_nu(math.hypot(_gRv, _gzv)) - 1) * math.hypot(_gRv, _gzv) * _h / math.hypot(_Rv, _h)
_WERT_V = {"ort": math.sqrt(_az_ort / _h), "kern": math.sqrt(_az_kern / _h)}
_n("vertikalfrequenz_feld", "vertikalfrequenz", "zweiter_weg",
   [(_wzN, _OmN, _nuv, "ort"), (_wzN, _OmN, _nuv, "kern")], lambda w, O, n_, form: _WERT_V[form],
   toleranz=1e-6, grund="ω_z² = a_z(h)/h aus dem Feld der Miyamoto-Nagai-Scheibe, h = 10⁻⁵ kpc")

_aM, _bM = 3 * kpc, 0.3 * kpc


def _praez_analytisch(M, a, b, R, regel=True):
    gR = G * M * R / (R * R + (a + b)**2)**1.5
    wzN2 = G * M * (a + b) / (b * (R * R + (a + b)**2)**1.5)           # ∂g_z/∂z in der Ebene
    n_ = _nu(gR) if regel else 1.0
    return math.sqrt(n_ * gR / R) - math.sqrt(n_ * wzN2)


_n("knotenpraezession_analytisch", "knotenpraezession_klein", "zweiter_weg",
   [(MW, _aM, _bM, 20 * kpc, True), (MW, _aM, _bM, 20 * kpc, False)], _praez_analytisch, toleranz=1e-6,
   grund="Differenzen mit h = 10⁻⁴ kpc gegen die Ableitung von Hand")
_n("praezession_ring_klein", "praezession_ring", "grenzfall", [(MW, _aM, _bM, 12 * kpc, math.radians(0.1), True)],
   lambda M, a, b, R, i, regel: K.knotenpraezession_klein(M, a, b, R, regel), toleranz=1e-3,
   grund="Neigung 0,1°: Bahnintegration gegen die Frequenzen am Ort, Rest ∝ i² (0,5°: 0,3 %)")


# ------------------------------------------------------------------ System als Ganzes im äußeren Feld
def _m1_quad(g_i, R, g_e, bezug, a0=A0):
    from scipy.integrate import quad

    def v(p):
        gx, gy = -g_i * math.cos(p) + g_e, -g_i * math.sin(p)
        n_ = _nu(math.hypot(gx, gy), a0)
        return math.sqrt(max(-R * ((n_ * gx - bezug) * math.cos(p) + n_ * gy * math.sin(p)), 0.0))
    F0 = quad(v, 0, 2 * math.pi, limit=200)[0] / (2 * math.pi)
    F1 = quad(lambda p: v(p) * math.cos(p), 0, 2 * math.pi, limit=200)[0] / (2 * math.pi)
    return 2 * abs(F1) / F0


_gi20 = G * MW / (20 * kpc)**2
_n("einseitigkeit_m1_quadratur", "einseitigkeit_m1", "zweiter_weg",
   [(_gi20, 20 * kpc, 0.03 * A0, 0.03 * A0), (_gi20 / 6.25, 50 * kpc, 0.03 * A0, _nu(0.03 * A0) * 0.03 * A0)],
   _m1_quad, toleranz=1e-6, grund="Fourier-Koeffizient mit quad statt FFT über 720 Punkte")
_n("massenmittel_A_plummer_schwach", "massenmittel_A_plummer", "grenzfall", [(1e5 * K.M_sun, 10 * kpc, 0.1 * A0)],
   lambda M, b, ge: _nu(ge) * ge, toleranz=1e-4, grund="g_i ≪ g_e: A → ν(g_e)·g_e, Korrektur ∝ g_i/g_e")
_n("massenmittel_A_scheibe_schwach", "massenmittel_A_scheibe", "grenzfall", [(1e5 * K.M_sun, 10 * kpc, 0.1 * A0)],
   lambda M, Rd, ge: _nu(ge) * ge, toleranz=1e-4, grund="g_i ≪ g_e: A → ν(g_e)·g_e, Korrektur ∝ g_i/g_e")


def _sigma_newton(M, Re, ge, a0=A0, regel=False):
    """Plummer isotrop: σ_r² = GM/(6√(r² + b²)); Fenster R < R_e über den Anteil jeder Schale im Zylinder."""
    from scipy.integrate import quad
    w = lambda r: 1.0 if r <= Re else 1 - math.sqrt(1 - (Re / r)**2)                    # noqa: E731
    rho = lambda r: (1 + (r / Re)**2) ** -2.5                                              # noqa: E731
    z = quad(lambda r: r * r * rho(r) * w(r) * G * M / (6 * math.hypot(r, Re)), 0, 1e4 * Re, points=[Re], limit=400)[0]
    n = quad(lambda r: r * r * rho(r) * w(r), 0, 1e4 * Re, points=[Re], limit=400)[0]
    return math.sqrt(z / n)


_n("sigma_plummer_jeans_newton", "sigma_plummer_jeans", "zweiter_weg",
   [(1e5 * K.M_sun, 0.02 * kpc, 0.0, A0, False), (1.5e8 * K.M_sun, 1.6 * kpc, 0.0, A0, False)], _sigma_newton,
   toleranz=1e-3, grund="Trapez auf 1500 logarithmischen Radien gegen σ_r² = GM/(6√(r² + b²)) mit quad")
_n("sigma_v_tief_handwert", "sigma_v_tief", "handwert", [(1e7 * K.M_sun,), (2e8 * K.M_sun,)],
   lambda M: 9.2e3 if M < 1e8 * K.M_sun else 19.5e3, toleranz=5e-3, grund="Modelltext 9,2 und 19,5 km/s (drei Stellen)")


# ------------------------------------------------------------------ Gezeit und Hierarchie der Trichter
_n("gezeit_betrag_handwert", "gezeit_betrag", "handwert",
   [(np.diag([1e-30, 2e-30, 3e-30]),), (np.array([[0, 2e-30, 0], [2e-30, 0, 0], [0, 0, 0]]),)],
   lambda E: math.sqrt(3) * (1e-30 if E[0, 0] else 2e-30))


def _gezeit_newton(v0, R0, A, B, rho, a0):
    g = v0 * v0 / R0
    dlnv = -(A + B) / (A - B)
    lam = [g / R0 * (2 * dlnv - 1), g / R0]
    lam.append(4 * math.pi * G * rho - lam[0] - lam[1])
    m = sum(lam) / 3
    return math.sqrt(1.5 * sum((x - m)**2 for x in lam))


_KK = 1e3 / kpc
_n("gezeit_sonnenort_newton", "gezeit_sonnenort", "grenzfall",
   [(233e3, 8.2 * kpc, 15.3 * _KK, -11.9 * _KK, 0.08 * K.M_sun / K.pc**3, 1e-15)], _gezeit_newton, toleranz=1e-9)
_n("gezeit_hill_frequenzen", "gezeit_hill", "zweiter_weg", [(15.3 * _KK, -11.9 * _KK)],
   lambda A, B: 4 * (A - B)**2 - (-4 * B * (A - B)))


def _brent(f, a, b):
    from scipy.optimize import brentq
    return brentq(f, a, b, xtol=1e-300, rtol=1e-14)


_T = 5.5e-30
_n("glaettungslaenge_frequenz", "glaettungslaenge", "zweiter_weg", [(K.M_sun, _T), (K.M_earth, _T)],
   lambda M, T: _brent(lambda r: G * M / r**3 - T, 1e10, 1e20), toleranz=1e-9)
_n("gezeitendichte_mittel", "gezeitendichte", "umkehr", [(_T,)],
   lambda T: K.M_sun / (4 / 3 * math.pi * K.glaettungslaenge(K.M_sun, T)**3), toleranz=1e-9)
_n("schnittmasse_glaettung_schnitt", "schnittmasse_glaettung", "zweiter_weg", [(_T, A0)],
   lambda T, a0: _brent(lambda M: math.log(K.glaettungslaenge(M, T) / K.r_stern(M, a0)), 1e30, 1e45), toleranz=1e-9)
_n("schnittradius_glaettung_schnitt", "schnittradius_glaettung", "zweiter_weg", [(_T, A0)],
   lambda T, a0: K.r_stern(_brent(lambda M: math.log(K.glaettungslaenge(M, T) / K.r_stern(M, a0)), 1e30, 1e45), a0), toleranz=1e-9)


# ------------------------------------------------------------------ Bose-Form und Spektraldimension
_xb = np.array([-1.5, -1.0, -0.5, 0.0])
_sb = np.array([0.01, 0.02, 0.01, 0.03])
_kb = np.array([1.0, -2.0, 0.5, 0.0])
_Bq = np.stack([_xb, [math.log10(_nu(10**x * A0, 1.2e-10)) for x in _xb] + _kb * _sb, np.full(4, 50), _sb], axis=1)
_n("fq_form_handwert", "fq_form", "handwert", [(_Bq, 1.2e-10)], lambda B, a0: 1.0 + 4.0 + 0.25, toleranz=1e-9)
_Ba = _Bq.copy()
_Ba[:, 1] = [math.log10(_nu(10**x * A0, 1.2e-10)) for x in _xb]
_n("a0_an_punkten_synthetisch", "a0_an_punkten", "grenzfall", [(_Ba,)], lambda B: 1.2e-10, toleranz=1e-4,
   grund="begrenzte Minimierung, xatol 10⁻⁵ in log a₀")
_n("nu_ds_bose", "nu_ds", "grenzfall", [(3e-11, 1.0), (5e-10, 1.0)], lambda x, d: _nu(x), toleranz=1e-12)
_gs = np.logspace(-13, -9, 60)
_gd = np.array([x / (-math.expm1(-(x / 1.65e-10) ** (1 - 1.3 / 2))) for x in _gs])
_n("ds_anpassung_synthetisch", "ds_anpassung", "grenzfall", [(_gs, _gd)], lambda gN, g: 1.3, toleranz=1e-5,
   grund="Anpassung (least_squares) an genaue Werte mit d_s = 1,3, a₀ = 1,65·10⁻¹⁰")


def _steigung(aN, a0=A0, h=1e-5):
    f = lambda x: math.log(_nu(x, a0) * x)                    # noqa: E731
    return 2 * (f(aN * math.exp(h)) - f(aN * math.exp(-h))) / (2 * h)


_n("spektraldimension_steigung", "spektraldimension_regel", "zweiter_weg", [(1e-12,), (1.1e-10,), (5e-9,)], _steigung,
   toleranz=1e-8, grund="zentrale Differenz in ln a_N")


# ------------------------------------------------------------------ Flussdimension, gelesene Dichte, Faltungsmesser
_n("flussdimension_handwert", "flussdimension", "handwert", [(1e-40, A0, 3), (1e-40, A0, 2), (1e-12, A0, 1)],
   lambda aN, a0, d: {3: 2.0, 2: 1.5, 1: 1.0}[d])


def _dichte_numerisch(aN, r, rho, a0=A0, h=1e-5):
    M = aN * r * r / G
    zusatz = lambda rr: (_nu(G * M / rr**2, a0) - 1) * G * M / rr**2 * rr * rr      # noqa: E731  r²(a − a_N)
    div = (zusatz(r * math.exp(h)) - zusatz(r * math.exp(-h))) / (2 * h) / r / (r * r)
    f = lambda x: _nu(x, a0) * x                                                      # noqa: E731
    dadaN = (f(aN * math.exp(h)) - f(aN * math.exp(-h))) / (aN * 2 * math.sinh(h))
    return div / (4 * math.pi * G) + rho * (dadaN - 1)


_n("dichte_gelesen_divergenz", "dichte_gelesen", "zweiter_weg",
   [(1e-11, 30 * kpc, 0.0), (5e-10, 3 * kpc, 1e-21), (2.8e-10, 5 * kpc, 1e-20)], _dichte_numerisch, toleranz=1e-6,
   grund="Divergenz r⁻²d(r²(a − a_N))/dr (Punktmasse außen) und d(νa_N)/da_N mit zentralen Differenzen")
_n("r_stern_filament_schnitt", "r_stern_filament", "zweiter_weg", [(1e12 * K.M_sun / K.Mpc,)],
   lambda mu: _brent(lambda r: 2 * G * mu / r - A0, 1e15, 1e25), toleranz=1e-9)
_n("r_stern_filament_handwert", "r_stern_filament", "handwert", [(1e12 * K.M_sun / K.Mpc,)], lambda mu: 2.5 * kpc,
   toleranz=0.02, grund="Modelltext 2,5 kpc (zwei Stellen)")
_n("faltungsmesser_masse", "faltungsmesser", "zweiter_weg", [(172e3,), (50e3,)], lambda v: K.ds(v**4 / (G * A0)),
   toleranz=1e-12)
_n("ds_eff_beschleunigung", "ds_eff", "zweiter_weg", [(MW, 30 * kpc), (K.M_sun, K.AU)],
   lambda M, r: K.a_regel(G * M / r**2) * r * r / (3 * K.c**2), toleranz=1e-12)
_n("masse_aus_temperatur_kreis", "masse_aus_temperatur", "umkehr", [(1.16e7, 2.0, 0.6)],
   lambda T, S, mu: _brent(lambda M: K.v_flach(M)**2 / (S * K.kB * T / (mu * K.m_p)) - 1, 1e35, 1e50), toleranz=1e-9)


# ------------------------------------------------------------------ Einsetzen als Kernformel (03)
_n("pulsar_rate_hand", "pulsar_rate", "handwert", [(1e-10,)], lambda a: 3.33564e-19, toleranz=1e-5,
   grund="10⁻¹⁰ m/s² / 2,99792458·10⁸ m/s = 3,33564·10⁻¹⁹ s⁻¹ von Hand, sechs Stellen")
_n("warp_viertelzeit_hand", "warp_viertelzeit", "handwert", [(1.0,), (2.0,)], lambda w: {1.0: 1.5707963, 2.0: 0.78539816}[w],
   toleranz=1e-7, grund="π/2 = 1,5707963 und π/4 = 0,78539816 von Hand")
_n("verstaerkung_druck_hand", "verstaerkung_druck", "handwert", [(2.0, 1.0), (3.0, 2.0)],
   lambda s, sN: {(2.0, 1.0): 4.0, (3.0, 2.0): 2.25}[(s, sN)])


# ------------------------------------------------------------------ Regel, Kreisbahn, Grenzflächendichte (von Hand)
_n("a_regel_hand", "a_regel", "handwert", [(0.873 * A0,)], lambda x: 1.43785 * A0, toleranz=1e-5,
   grund="a_N = 0,873a₀: s = 0,93434, e^{−s} = 0,39284, ν = 1,6470, a = 1,43785a₀ von Hand (sechs Stellen)")
_n("v_kreis_punkt_hand", "v_kreis_punkt", "zweiter_weg", [(MW, 8 * kpc), (MW, 50 * kpc)],
   lambda M, r: math.sqrt(_nu(G * M / r**2) * G * M / r))
_n("v_flach_schnitt", "v_flach", "zweiter_weg", [(MW,)], lambda M: _brent(lambda v: v**4 / (G * M * A0) - 1, 1e3, 1e7),
   toleranz=1e-9)
_n("nu_aussenfeld_hand", "nu_aussenfeld", "zweiter_weg", [(0.1 * A0, 0.03 * A0, "senkrecht"), (0.1 * A0, 0.03 * A0, "gleich")],
   lambda gi, ge, lage: _nu((gi**4 + ge**4 + 2 * gi**2 * ge**2) ** 0.25 if lage == "senkrecht" else gi + ge))
_n("sigma_kollektiv_hand", "sigma_kollektiv", "handwert", [(1.1e-10,)], lambda a0: 251.20, toleranz=5e-4,
   grund="1,1·10⁻¹⁰/(π·6,6743·10⁻¹¹) kg/m² = 0,52461 kg/m², M☉/pc² = 2,0884·10⁻³ kg/m²: 251,20 von Hand")


# ------------------------------------------------------------------ Feldform in Scheiben (Löser, Einheiten G = a₀ = 1)
def _kuzmin_eta(profil, h, x):
    """Kuzmin: Feldform exakt ν(|g_N|)·g_R, |g_N| = 1/(R² + h²) in der Ebene; η = ν(g_R)/ν(|g_N|)."""
    R = x * h
    return _nu(R / (R * R + h * h) ** 1.5, 1.0) / _nu(1 / (R * R + h * h), 1.0)


_n("eta_scheibe_kuzmin", "eta_scheibe", "grenzfall", [("kuzmin", 3.0, 1.0), ("kuzmin", 0.3, 2.0)], _kuzmin_eta,
   toleranz=3e-3, grund="Kuzmin-Scheibe: Feldform exakt ν(|g_N|)·g_N; Gitter 200 × 200 (Kuzmin h = 3 von Hand 1,168)")
_n("feldform_naeherung_kuzmin", "feldform_naeherung", "grenzfall", [("kuzmin", 3.0, 1.0), ("kuzmin", 0.3, 2.0)],
   lambda profil, h, x: 1.0, toleranz=3e-3,
   grund="Kuzmin-Scheibe: in der Ebene ist |g_N| = √(g_R² + (2πGΣ)²) exakt, die Näherung ist die Feldform; Gitter")
_n("gleichstand_scheibe_nullstelle", "gleichstand_scheibe", "umkehr", [(1.0,)],
   lambda h: _brent(lambda x: K.eta_scheibe("exp", h, x) - 1, 2.0, 6.0), toleranz=2e-3,
   grund="Nullstelle von η − 1 mit brentq über die interpolierte η statt Interpolation zwischen zwei Zellmitten")
_LK = K.feldform_loeser(K.scheibe_quelle("kuzmin", 1.0))


def _kuzmin_saeule(los, R, z1, feld):
    """Kuzmin (M = 1, h = 1) oberhalb der Ebene: Feld einer Punktmasse bei z = −h; Feldform a = ν(|g|)·g exakt."""
    d2 = R * R + (z1 + 1.0)**2
    gz = (z1 + 1.0) / d2**1.5
    return (gz if feld == "newton" else _nu(1 / d2, 1.0) * gz) / (2 * math.pi)


_n("saeule_loesung_kuzmin", "saeule_loesung", "grenzfall", [(_LK, 3.0, 1.0, "newton"), (_LK, 3.0, 1.0, "feldform")],
   _kuzmin_saeule, toleranz=3e-3, grund="Kuzmin-Scheibe oberhalb der Ebene exakt (Punktmasse gespiegelt); Gitter")


# ------------------------------------------------------------------ Sonnenort, Sattel, Ringe, Warp
_aN0 = 1.3e-10
_vN0 = math.sqrt(_nu(_aN0) * _aN0 * 8 * kpc)
_n("a_N_aus_kreisbahn_vorwaerts", "a_N_aus_kreisbahn", "zweiter_weg", [(_vN0, 8 * kpc)], lambda v, r: _aN0,
   toleranz=1e-9, grund="v aus a_N = 1,3·10⁻¹⁰ vorwärts mit ν von Hand")
_n("sattelpunkt_nullstelle", "sattelpunkt", "zweiter_weg", [(K.GM_sun, K.GM_earth, K.AU)],
   lambda G1, G2, D: _brent(lambda x: G2 / x**2 - G1 / (D - x)**2, 1e6, 1e10), toleranz=1e-9)
_Sb = 47.1 * K.M_sun / K.pc**2


def _saeule_hand(S, gR, z, R, form):
    gz = 2 * math.pi * G * S
    if form == "kern":
        return S + (_nu(gR) - 1) * gR * z / R / (2 * math.pi * G)
    return (_nu(math.hypot(gR, gz)) if form == "ort_schicht" else _nu(gR)) * S


_n("dynamikmasse_ort_hand", "dynamikmasse_ort", "zweiter_weg",
   [(_Sb, 1.374e-10, 1.1 * kpc, 8.2 * kpc, f) for f in ("kern", "ort_schicht", "ort_ebene")], _saeule_hand)
_n("knotenpraezession_frequenzen_hand", "knotenpraezession_frequenzen", "handwert",
   [(1e-15, 3e-15, 1.8, "ort"), (1e-15, 3e-15, 1.8, "kern")],
   lambda O, w, n_, form: math.sqrt(n_) * (O - w) if form == "ort" else math.sqrt(n_) * O - math.sqrt(w * w + (n_ - 1) * O * O))


def _praez_freeman(M, Rd, h, R, form):
    """Freeman-Scheibe mit den Bessel-Funktionen von scipy, Ableitung von v² über die Bessel-Ableitungen."""
    from scipy.special import i0, i1, k0, k1
    y = R / (2 * Rd)
    S0 = M / (2 * math.pi * Rd * Rd)
    B = i0(y) * k0(y) - i1(y) * k1(y)
    dB = i1(y) * k0(y) - i0(y) * k1(y) - (i0(y) - i1(y) / y) * k1(y) - i1(y) * (-k0(y) - k1(y) / y)
    v2 = 4 * math.pi * G * S0 * Rd * y * y * B
    dv2 = 4 * math.pi * G * S0 * Rd * (2 * y * B + y * y * dB) / (2 * Rd)
    OmN = math.sqrt(v2) / R
    wzN = math.sqrt(4 * math.pi * G * S0 * math.exp(-R / Rd) / (2 * h) - dv2 / R)
    n_ = _nu(v2 / R)
    wz = math.sqrt(n_) * wzN if form == "ort" else math.sqrt(wzN**2 + (n_ - 1) * OmN**2)
    return math.sqrt(n_) * OmN - wz


_n("knotenpraezession_scheibe_bessel", "knotenpraezession_scheibe", "zweiter_weg",
   [(MW, 2.6 * kpc, 0.3 * kpc, 20 * kpc, "ort"), (MW, 2.6 * kpc, 0.3 * kpc, 20 * kpc, "kern")], _praez_freeman,
   toleranz=1e-6, grund="dv²/dR mit zentraler Differenz (10⁻⁴R) gegen die Bessel-Ableitung von Hand")
_n("warp_aufwicklung_hand", "warp_aufwicklung", "handwert", [(-3e-16, -1e-16, math.pi)], lambda a, b, d: math.pi / 2e-16)


# ------------------------------------------------------------------ SPARC: synthetische Galaxien aus der Regel
def _syn(nm, Md, rd, a0=A0, f=1.0, B=0.0, gas=0.0, SB0=100.0, D=10.0, inc=60.0, n=20, rmax=12.0, sig_rel=0.02):
    """Galaxie der Tafel aus der Regel: Scheibe (Kuzmin in der Ebene, Υ = 0,5), Gas als Anteil gas·R/(R + r_d) der Scheibe;
    V_obs² = (f·ν(g_N; a₀)·g_N + B·a₀)·R. Fehler e_D, e_i null, e_V = sig_rel·V."""
    R = np.linspace(0.3, rmax, n) * rd * kpc
    Vd2 = G * Md * K.M_sun * R**2 / (R**2 + (rd * kpc)**2)**1.5 / 0.5
    Vg = np.sqrt(gas * 0.5 * Vd2 * R / (R + rd * kpc))
    gN = (Vg * Vg + 0.5 * Vd2) / R
    Vo = np.sqrt((f * np.array([_nu(x, a0) for x in gN]) * gN + B * a0) * R)
    return dict(nm=nm, D=D * K.Mpc, eD=0.0, fD=2, inc=inc, eInc=0.0, R=R, Vo=Vo, eV=sig_rel * Vo, Vg=Vg, Vd=np.sqrt(Vd2),
                Vb=np.zeros(n), Rd=rd * kpc, SB0=SB0, SBd=SB0 * np.exp(-R / (rd * kpc)), SBb=np.zeros(n))


_SYN = [_syn(f"s{i}", m, rd) for i, (m, rd) in enumerate(((1e9, 1.0), (5e9, 2.0), (2e10, 3.0), (6e10, 3.5), (2.5e8, 0.8)))]
_SYN = [_syn(g["nm"], m, rd, a0=1.25e-10) for g, (m, rd) in zip(_SYN, ((1e9, 1.0), (5e9, 2.0), (2e10, 3.0), (6e10, 3.5),
                                                                      (2.5e8, 0.8)))]
_n("a0_profiliert_synthetisch", "a0_profiliert", "grenzfall", [(K.sparc_profil_satz(_SYN), 0.5, 0.1, "linear", None, 0)],
   lambda *a: 1.25e-10, toleranz=2e-3, grund="Galaxien genau nach der Regel mit a₀ = 1,25·10⁻¹⁰; Parabel im Gitter "
   "(Schritt 0,01 dex)")
_SAT4 = K.sparc_profil_satz(K.sparc_tafel(K.SPARC_SATZ["q_max"], K.SPARC_SATZ["inc_min"])[:4])


def _profil_ls(satz, a0, ups_d, sig_int):
    """Profilierte Fehlerquadratsumme je Galaxie mit scipy least_squares (Trust Region), eigene Residuen."""
    from scipy.optimize import least_squares
    su = math.log(10) * K.SPARC_SATZ["ups_streuung_dex"]
    aus = []
    for k in range(len(satz["name"])):
        m = satz["maske"][k]
        lgo, sv = satz["lgo"][k, m], satz["sv"][k, m]
        gas, disk, bul = satz["gas"][k, m], satz["disk"][k, m], satz["bul"][k, m]
        i0, eD, ei = satz["inc"][k], satz["eD"][k], satz["ei"][k]

        def res(p):
            uD, di, uY = p
            inc = min(max(i0 + di, 0.2), math.pi / 2)
            gb = np.maximum(gas + math.exp(uY) * (ups_d * disk + K.SPARC_SATZ["ups_b"] * bul), 1e-14)
            a = np.array([_nu(x, a0) for x in gb]) * gb
            r = (lgo - uD + 2 * math.log(math.sin(i0) / math.sin(inc)) - np.log(a)) / np.sqrt(sv**2 + sig_int**2)
            return np.r_[r, uD / eD, di / ei, uY / su]
        aus.append(2 * least_squares(res, [0.0, 0.0, 0.0], xtol=1e-14, ftol=1e-14, gtol=1e-14).cost)
    return np.array(aus)


_n("profil_chi2_least_squares", "profil_chi2", "zweiter_weg", [(_SAT4, 1.1e-10, 0.5, 0.1)], _profil_ls, toleranz=1e-6,
   grund="Levenberg-Marquardt je Galaxie gegen scipy least_squares (Trust Region), gleiche Minima auf 10⁻⁶")


def _syn3():
    """Drei Gruppen zu je fünf Galaxien, je Gruppe eine Flächendichte und eine Schwelle."""
    aus = []
    for j, (a0, sb) in enumerate(((0.8e-10, 20.0), (1.0e-10, 200.0), (1.25e-10, 2000.0))):
        for i, (m, rd) in enumerate(((1e9, 1.0), (5e9, 2.0), (2e10, 3.0), (6e10, 3.5), (2.5e8, 0.8))):
            aus.append(_syn(f"d{j}{i}", m, rd, a0=a0, SB0=sb))
    return aus


def _drittel_hand(*a):
    x = np.log10(0.5 * np.array([20.0, 200.0, 2000.0]))
    y = np.log10([0.8e-10, 1.0e-10, 1.25e-10])
    return 0.8e-10, 1.0e-10, 1.25e-10, np.polyfit(x, y, 1)[0]


_n("a0_je_drittel_synthetisch", "a0_je_drittel", "grenzfall", [(K.sparc_profil_satz(_syn3()), 0.5, 0)], _drittel_hand,
   toleranz=3e-3, grund="je Gruppe Galaxien genau nach der Regel mit eigener Schwelle; Parabel im Gitter (49 Punkte)")
_SATZ_ALLE = K.sparc_profil_satz(K.sparc_tafel(K.SPARC_SATZ["q_max"], K.SPARC_SATZ["inc_min"]))
_n("a0_senkung_ohne_feld", "a0_senkung_aussenfeld", "grenzfall", [(_SATZ_ALLE, 0.5, 0.0, "senkrecht", True)],
   lambda *a: A0, toleranz=1e-8, grund="ohne äußeres Feld liest man a₀ selbst; Minimierung mit xatol 10⁻¹⁰")


def _syn_toomre():
    """Galaxien mit V_obs = V_Regel·10^{ε/2}: der Rest der zentralen Flächendichte ist genau ε je Galaxie."""
    aus, eps = [], []
    for i, (m, rd, sb) in enumerate(((1e9, 1.0, 20.0), (5e9, 2.0, 60.0), (2e10, 3.0, 200.0), (6e10, 3.5, 700.0),
                                     (2.5e8, 0.8, 10.0), (1e10, 2.5, 300.0))):
        e = 0.1 + 0.05 * (math.log10(0.5 * sb) - 2) + (0.01 if i % 2 else -0.01)
        g = _syn(f"t{i}", m, rd, SB0=sb)
        g["Vo"] = g["Vo"] * 10 ** (e / 2)
        aus.append(g)
        eps.append((math.log10(0.5 * sb), e))
    return aus, np.array(eps)


_TT, _TE = _syn_toomre()


def _toomre_hand(tafel, ups_d, a0=A0):
    x, y = _TE[:, 0], _TE[:, 1]
    k, c = np.polyfit(x, y, 1)
    s2 = np.sum((y - (k * x + c))**2) / (y.size - 2)
    return y.mean(), y.std(ddof=1) / math.sqrt(y.size), k, math.sqrt(s2 / np.sum((x - x.mean())**2))


_n("sigma_dyn_rest_synthetisch", "sigma_dyn_rest", "grenzfall", [(_TT, 0.5)], _toomre_hand, toleranz=1e-9,
   grund="Toomre-Integral linear in V²: der Rest ist das eingesetzte ε")


def _syn_rest(k_gas, c_aussen=0.0, k_log=0.0):
    """Galaxien mit Rest c_g + k_gas·f_gas + k_log·log₁₀(g_N/a₀) (σ = 0), Gas als halbe Scheibe, weit hinaus."""
    aus = []
    for i, (m, rd) in enumerate(((1e8, 1.0), (2.5e8, 1.5), (1e9, 2.0), (4e8, 1.2))):
        g = _syn(f"r{i}", m, rd, gas=0.3 + 0.2 * i, rmax=40.0, n=30)
        Vb2 = g["Vg"]**2 + 0.5 * g["Vd"]**2
        gN = Vb2 / g["R"]
        f = g["Vg"]**2 / Vb2
        r = 0.02 * i + k_gas * f + k_log * np.log10(gN / A0) + c_aussen
        g["Vo"] = np.sqrt(10**r * np.array([_nu(x) for x in gN]) * gN * g["R"])
        aus.append(g)
    return aus


def _gas_hand(tafel, *a):
    """Rest = c_g + 0,03·f: Nullpunkt je Galaxie (Median) abgezogen, Steigung gegen f − ⟨f⟩ über alle Punkte."""
    F = [g["Vg"]**2 / (g["Vg"]**2 + 0.5 * g["Vd"]**2) for g in tafel]
    x = np.concatenate(F)
    y = np.concatenate([0.03 * (f - np.median(f)) for f in F])
    return np.sum((x - x.mean()) * y) / np.sum((x - x.mean())**2)


_n("rest_gasanteil_synthetisch", "rest_gasanteil", "zweiter_weg", [(_syn_rest(0.03), 0.5, 0.0, A0, 0)], _gas_hand,
   toleranz=1e-9, grund="Rest genau linear im Gasanteil; Schätzer von Hand")
_aQ = math.log(K.q_KOSM)**2 * A0
_n("rest_jenseits_synthetisch", "rest_jenseits", "grenzfall", [(_syn_rest(0.0), 0.5, 0.0, _aQ, A0, 0, 0)],
   lambda *a: 0.03, toleranz=1e-9, grund="Rest je Galaxie 0,02·i konstant: Median der vier Galaxien 0,03")
_n("reststeigung_synthetisch", "reststeigung_aussen", "grenzfall",
   [(_syn_rest(0.0, 0.0, -0.3), 0.5, 0.0, 0.0, "senkrecht", 0.05 * A0, A0, 0)], lambda *a: -0.3, toleranz=1e-9,
   grund="Rest genau linear in log₁₀ a_N mit Steigung −0,3")
_GAS = [_syn(f"g{i}", 5e8, 1.5, gas=1.0) for i in range(4)]
for _g in _GAS:
    _g["Vd"] = np.zeros_like(_g["Vd"])
    _g["Vo"] = np.sqrt((np.array([_nu(x) for x in _g["Vg"]**2 / _g["R"]]) * _g["Vg"]**2 / _g["R"] + 0.02 * A0) * _g["R"])
_n("ds_streuung_eine_ziehung", "ds_streuung", "grenzfall", [(_GAS, 0.5, 0.0, 1, 1)], lambda *a: 0.0,
   grund="eine Ziehung gleicher Galaxien nur aus Gas, ohne Fehler von Abstand und Neigung: Streuung 0")
_n("boden_grenze_ohne_fehler", "boden_grenze", "grenzfall", [(_GAS, 0.5, 0.0, 0.3, 5, 1)], lambda *a: 0.02,
   toleranz=1e-6, grund="gleiche Galaxien nur aus Gas mit Boden 0,02·a₀, ohne Fehler: Grenze = Boden")
_gb = np.logspace(-13, -10, 50)
_n("boden_anpassung_synthetisch", "boden_anpassung", "grenzfall",
   [(_gb, np.array([_nu(x) * x for x in _gb]) + 0.02 * A0, 0.3)], lambda *a: 0.02, toleranz=1e-6,
   grund="Punkte genau nach a = ν(a_N)·a_N + 0,02·a₀")
_n("boden_drift_hand", "boden_drift", "handwert", [(0.016,)], lambda Pe: 0.032 * A0)
_n("boden_quer_hand", "boden_quer", "handwert", [(1.0, 25.0)], lambda v, b: 0.7978846 / 25 * A0, toleranz=1e-7,
   grund="√(2/π) = 0,7978846 von Hand")


def _syn_f(f, a0):
    """Galaxien ohne Gas nach a = ν(√(g_R² + f(2πGΣ)²))·g_R mit Σ = Υ_d·SB_disk (Υ_d = 0,5)."""
    aus = []
    for i, (m, rd, sb) in enumerate(((1e9, 1.0, 50.0), (5e9, 2.0, 200.0), (2e10, 3.0, 500.0), (2.5e8, 0.8, 30.0))):
        g = _syn(f"f{i}", m, rd, SB0=sb)
        gR = 0.5 * g["Vd"]**2 / g["R"]
        gz = 2 * math.pi * G * 0.5 * g["SBd"] * K.M_sun / K.pc**2
        g["Vo"] = np.sqrt(np.array([_nu(math.sqrt(x * x + f * z * z), a0) for x, z in zip(gR, gz)]) * gR * g["R"])
        aus.append(g)
    return aus


_AGf = np.exp(np.linspace(math.log(0.7e-10), math.log(1.9e-10), 61))
_n("gewicht_f_synthetisch", "gewicht_f", "grenzfall", [(_syn_f(0.3, _AGf[30]), 0.5, (), False)], lambda *a: 0.3,
   toleranz=0.02, grund="Galaxien genau nach der Regel mit f = 0,3 und a₀ auf dem Gitter; Parabel im Gitter von f (0,02)")



def _quadrupol_simpson(GM, g_e, ell, a0=A0):
    """Dasselbe Integral mit Simpson auf gleichmäßigen Gittern (2001 Punkte in ln ρ, 801 in θ), eigener Code."""
    from scipy.integrate import simpson
    re = math.sqrt(GM / g_e)
    beta, lam = g_e / a0, ell / re
    u, th = np.linspace(-10, 12, 2001), np.linspace(0, math.pi, 801)
    R_, TH = np.meshgrid(np.exp(u), th, indexing="ij")
    cc = np.cos(TH)
    h = R_ / (R_**2 + lam**2)**1.5
    gm = np.sqrt(np.maximum(h**2 + 2 * cc * h + 1, 0.0))
    k = -3 * (3 * cc**2 - 1) * h / R_**2 + (9 * cc - 15 * cc**3) / R_**2
    nu_ = 1 / (-np.expm1(-np.sqrt(np.maximum(beta * gm, 1e-300))))
    J = simpson(simpson(0.5 * (nu_ - 1) * k * R_ * R_ * np.sin(TH), x=th, axis=1), x=u)
    s = math.sqrt(beta)
    J += (4 / 15) * (-s / math.expm1(s) / 2) * _nu(beta, 1.0) / lam**3
    return -1.5 * g_e / re * J


_n("quadrupol_glaettung_simpson", "quadrupol_glaettung", "zweiter_weg", [(K.GM_sun, 1.34 * A0, 0.93 * K.pc)],
   _quadrupol_simpson, toleranz=1e-5, grund="Simpson auf gleichmäßigen Gittern gegen Gauß-Legendre mit Stützstellen")
