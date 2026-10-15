# -*- coding: utf-8 -*-
"""
Gegenproben der Gruppe Teilchen, Sterne und Licht (02, 04, 08, 09): je Formel ein unabhängiger Weg, der die geprüfte
Formel nicht aufruft (zweiter Rechenweg, Grenzfall der ART bzw. Newtons, Handwert aus dem Modelltext).
Kennungen: "<Formel>_<Weg>". Toleranz über 1e-9 nur mit Grund (Handwert: halbe letzte Stelle des Modelltexts).
"""
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

import kern as K
from kern.messwert import MESSWERTE

from .gegenprobe import GEGENPROBEN
from .gegenproben import hawking_lebensdauer

_g = GEGENPROBEN.neu
G, c, hbar, kB = K.G, K.c, K.hbar, K.kB
ARCSEK = 180 / math.pi * 3600

# Handwerte ([G]) aus dem Messwert-Register
_W = {k: MESSWERTE[k].wert for k in ("MERKUR_A", "MERKUR_E", "MERKUR_P", "ERDE_J", "M87_M", "M87_D", "GW150914_MF",
                                     "GW150914_F", "PSR1913_M1", "PSR1913_M2", "PSR1913_PB", "PSR1913_E")}
M_SUN, M_EARTH, R_SUN = K.M_sun, K.M_earth, K.R_sun
M_MW = 6e10 * M_SUN


def _nu(aN, a0):
    """Regel ν = 1/(1 − e^{−√(a_N/a₀)}), hier unabhängig vom Kern geschrieben."""
    return 1 / (-math.expm1(-math.sqrt(aN / a0)))


def _polytrope(n):
    """Eigene Lane-Emden-Lösung (Taylor-Start, dichte Ausgabe): (x₁, m₁, Lösung)."""
    x0 = 1e-5
    y0 = [1 - x0**2 / 6 + n * x0**4 / 120, -x0 / 3 + n * x0**3 / 30]

    def rand(x, y):
        return y[0]
    rand.terminal, rand.direction = True, -1
    s = solve_ivp(lambda x, y: [y[1], -max(y[0], 0.0)**n - 2 * y[1] / x], [x0, 30], y0, method="DOP853",
                  events=rand, rtol=1e-12, atol=1e-14, dense_output=True)
    x1 = s.t_events[0][0]
    return x1, -x1**2 * s.y_events[0][0][1], s


# ------------------------------------------------------------------ Faltung (01): Handwerte der Sonne
_g("ds_sonne_handwert", "ds", "handwert", [(M_SUN,)], lambda M: 492.2, toleranz=1.1e-4,
   grund="Modelltext 492,2 m, halbe letzte Stelle")
_g("r0_sonne_handwert", "r0", "handwert", [(M_SUN,)], lambda M: 2.95e3, toleranz=1.7e-3,
   grund="Modelltext 2,95 km, halbe letzte Stelle")

# ------------------------------------------------------------------ Teilchen (02, 04)
_g("compton_red_mittel", "compton_red", "zweiter_weg", [(K.m_e,)], lambda m: math.sqrt(K.r_e * K.a_B),
   toleranz=1e-10, grund="λ̄_e² = r_e·a_B mit drei CODATA-Werten (Rundung 10⁻¹¹)")
_g("alpha_G_faltung", "alpha_G", "zweiter_weg", [(K.m_p,), (K.m_e,), (M_SUN,)],
   lambda M: 3 * (G * M / (3 * c**2)) / (hbar / (M * c)))
_g("sternmasse_faltung", "sternmasse", "zweiter_weg", [()],
   lambda: 3 * c**2 / G * math.sqrt((hbar / (K.m_p * c) / 3)**3 / (G * K.m_p / (3 * c**2))))
_g("grenzmasse_hawking_temperatur", "grenzmasse_hawking", "zweiter_weg", [(K.T_CMB,), (1e-8,), (300.0,)],
   lambda T: 3 * c**2 / G * hbar * c / (24 * math.pi * kB * T))
_g("masse_r_stern_gleich_compton_nullstelle", "masse_r_stern_gleich_compton", "zweiter_weg", [(K.A0,), (1.2e-10,)],
   lambda a0: 10**brentq(lambda lm: 0.5 * math.log(G * 10**lm / a0) - math.log(hbar / (10**lm * c)), -60, 10,
                         xtol=1e-15))
_g("eddington_thomson", "eddington", "zweiter_weg", [(M_SUN,), (10 * M_SUN,)],
   lambda M: 4 * math.pi * G * M * K.m_p * c / (8 * math.pi / 3 * K.r_e**2))
_g("eddington_sonne_handwert", "eddington", "handwert", [(M_SUN,)], lambda M: 1.257e31, toleranz=4e-4,
   grund="Wert 1,257·10³¹ W, halbe letzte Stelle")
_g("lebensdauer_falte_page", "lebensdauer_falte", "zweiter_weg", [(K.m_p,), (1e12,), (M_SUN,)], hawking_lebensdauer)


def _carter_sonne():
    """√𝒞 = (2π³/45)·f_s⁴·(Δs/Δs_*)·(α²/δ)²·L_Edd/L_* mit L_* = 4πR²σT_eff⁴ (Identität der Sternzahlen)."""
    sigma_sb = math.pi**2 * kB**4 / (60 * hbar**3 * c**2)
    L = 4 * math.pi * R_SUN**2 * sigma_sb * K.TEFF_sun**4
    fs = kB * K.TEFF_sun / (K.alpha**2 * K.m_e * c**2)
    M_st = (hbar * c / G)**1.5 / K.m_p**2
    delta = G * M_SUN / (R_SUN * c**2)
    L_edd = 1.5 * G * M_SUN * K.m_p * c / K.r_e**2
    return ((2 * math.pi**3 / 45) * fs**4 * (M_SUN / M_st) * (K.alpha**2 / delta)**2 * L_edd / L)**2


_g("carter_faktor_sonne", "carter_faktor", "zweiter_weg", [()], _carter_sonne)
_g("carter_faktor_frequenzleiter", "carter_faktor", "zweiter_weg", [()],
   lambda: (G * K.m_p**3 * c / hbar**2) * (K.m_p * c**2 / hbar)**3 / (K.alpha * K.alpha * c / K.a_B)**4)

# ------------------------------------------------------------------ Tiefe und Temperatur (02)


def _f_potential(n):
    """f = Δδ_innen/((n + 1)δ): Potentialunterschied Hülle – Zentrum aus dem Dichteprofil, in GM/R."""
    x1, m1, s = _polytrope(n)
    dphi = quad(lambda x: -s.sol(x)[1] if x > 1e-5 else x / 3, 0, x1, limit=200, epsabs=0, epsrel=1e-12)[0]
    return dphi / (m1 / x1) / (n + 1)


_g("tiefe_faktor_potential", "tiefe_faktor", "zweiter_weg", [(1.5,), (3.0,)], _f_potential, toleranz=1e-8,
   grund="eigene Integration und Quadratur (rtol 10⁻¹²)")
_g("tiefe_faktor_analytisch", "tiefe_faktor", "grenzfall", [(0.0,), (1.0,)], lambda n: 0.5, toleranz=1e-8,
   grund="n = 0: x₁ = √6, m₁ = 2√6; n = 1: x₁ = m₁ = π; beide f = ½ (Integration rtol 10⁻¹¹)")
_g("temperatur_kern_sonne_hand", "temperatur_kern", "handwert", [(M_SUN, R_SUN, 0.62, 3.0)],
   lambda M, R, mu, n: 0.8543 * G * M / R * mu * K.m_u / kB, toleranz=1e-4,
   grund="θ₃ = 0,8543 auf vier Stellen (Tafel, von Hand)")


def _chandrasekhar(mu_e, n):
    """M_Ch = (√(3π)/2)·m₁(3)·(ħc/G)^{3/2}/(μ_e m_u)² (nur n = 3)."""
    return math.sqrt(3 * math.pi) / 2 * _polytrope(3.0)[1] * (hbar * c / G)**1.5 / (mu_e * K.m_u)**2


_g("grenzmasse_zaehlregel_chandrasekhar", "grenzmasse_zaehlregel", "zweiter_weg", [(2.0, 3.0), (2.15, 3.0)],
   _chandrasekhar, toleranz=1e-9)
_g("grenzmasse_zaehlregel_handwert", "grenzmasse_zaehlregel", "handwert", [(2.0, 3.0)], lambda mu_e, n: 1.456 * M_SUN,
   toleranz=3.5e-4, grund="Modelltext 1,456 M☉, halbe letzte Stelle")


def _tmax_langsam(M, mu, mu_e, n):
    """Kleine Masse: Π = x²/5, kT/(μ m_u c²) = a/R − B/R², Maximum a²/(4B) geschlossen."""
    m_u = K.m_u
    x1, m1, _ = _polytrope(n)
    D = x1**3 / (3 * m1)
    f = x1 / ((n + 1) * m1)
    lam_e = hbar / (K.m_e * c)
    a = f * G * M / c**2
    B = K.m_e / (mu_e * m_u) / 5 * lam_e**2 * (3 * math.pi**2 * D * 3 * M / (4 * math.pi * mu_e * m_u))**(2 / 3)
    return mu * m_u * c**2 * a * a / (4 * B) / kB


_g("temperatur_max_langsam", "temperatur_max", "grenzfall", [(1e-3 * M_SUN, 0.62, 1.143, 1.5), (3e-4 * M_SUN, 1.0, 2.0, 3.0)],
   _tmax_langsam, toleranz=3e-4, grund="nichtrelativistische Näherung Π = x²/5 (x ≈ 0,02)")


def _decke():
    """max_x[x/4 − Π(x)] über die Nullstelle der Ableitung, Π eigenständig geschrieben."""
    def Pi(x):
        return (x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x)) / (8 * x**3)

    def ab(x, h=1e-6):
        return 0.25 - (Pi(x + h) - Pi(x - h)) / (2 * h)
    x = brentq(ab, 0.3, 5.0, xtol=1e-13)
    return x / 4 - Pi(x)


_g("decke_faktor_ableitung", "decke_faktor", "zweiter_weg", [()], _decke, toleranz=1e-9)

# ------------------------------------------------------------------ Leiter der Längen und Schnittmassen (01, 08)
_LP = math.sqrt(hbar * G / c**3)
_MPL = math.sqrt(hbar * c / G)
_g("leiterzahl_beschleunigung", "leiterzahl", "zweiter_weg", [(K.A0,), (1.2e-10,)], lambda a0: a0 / (c**2 / _LP))
_g("leiterzahl_handwert", "leiterzahl", "handwert", [(K.A0,)], lambda a0: 2.0e-62, toleranz=2.5e-2,
   grund="Modelltext 2,0·10⁻⁶², halbe letzte Stelle")
_g("schnittmasse_unten_planck", "schnittmasse_unten", "zweiter_weg", [(K.A0,), (1.2e-10,)],
   lambda a0: a0 * _LP / c**2 * _MPL)
_g("schnittmasse_oben_produkt", "schnittmasse_oben", "zweiter_weg", [(K.A0,), (1.2e-10,)],
   lambda a0: _MPL**2 / (hbar * a0 / c**3))
_g("schnittmasse_oben_handwert", "schnittmasse_oben", "handwert", [(K.A0,)], lambda a0: 1.1e54, toleranz=4.6e-2,
   grund="Modelltext 1,1·10⁵⁴ kg, halbe letzte Stelle")
_g("schwelle_zu_huelle_flaeche", "schwelle_zu_huelle", "zweiter_weg",
   [(M_SUN, R_SUN, K.A0), (M_MW, 15 * K.kpc, K.A0)],
   lambda M, R, a0: math.sqrt(math.pi * G * (M / (math.pi * R**2)) / a0))
_g("schwelle_zu_huelle_sonne_handwert", "schwelle_zu_huelle", "handwert", [(M_SUN, R_SUN, K.A0)],
   lambda M, R, a0: 1.6e6, toleranz=3.2e-2, grund="Modelltext 1,6·10⁶, halbe letzte Stelle")
_g("spiegel_masse_faltung", "spiegel_masse", "zweiter_weg", [(K.m_p,), (M_SUN,)],
   lambda M: c**2 / G * hbar / (M * c))   # 3Δs(M̃) = λ̄(M)
_g("spiegel_radius_schliesslaenge", "spiegel_radius", "grenzfall",
   [(2 * G * M / c**2, M, a0) for M in (M_SUN, M_MW) for a0 in (K.A0, 1.2e-10)],
   lambda r, M, a0: c**2 / (2 * a0), grund="r = r₀ = 6Δs wird auf Λ = c²/(2a₀) abgebildet", toleranz=1e-12)

# ------------------------------------------------------------------ Raten und Kopplungen (09)
_g("verhaeltnis_em_grav_laengen", "verhaeltnis_em_grav", "zweiter_weg", [(K.m_p,)],
   lambda m: (K.r_e / (hbar / (K.m_e * c))) / (G * m / c**2 / (hbar / (m * c))), toleranz=1e-10,
   grund="α = r_e/λ̄_e aus CODATA-Werten (Rundung 10⁻¹¹)")
_g("compton_takt_energie", "compton_takt", "zweiter_weg", [(K.m_p,), (K.m_e,)], lambda M: M * c**2 / hbar)
_g("faltungsrate_mittel", "faltungsrate", "zweiter_weg", [(K.m_p,), (M_SUN,)],
   lambda M: (M * c**2 / hbar)**2 / (c**3 / (G * M)))   # ω_C² = Γ_G·c/(3Δs)
_g("faltungsrate_proton_handwert", "faltungsrate", "handwert", [(K.m_p,)], lambda M: 8.42e-15, toleranz=6e-4,
   grund="Modelltext 8,42·10⁻¹⁵ s⁻¹, halbe letzte Stelle")
_g("faltungstakt_proton_handwert", "faltungstakt", "handwert", [(K.m_p,)], lambda M: 2.41e62, toleranz=2.1e-3,
   grund="Modelltext 2,41·10⁶² s⁻¹, halbe letzte Stelle")
_g("faltungstakt_faltung", "faltungstakt", "zweiter_weg", [(K.m_p,), (M_SUN,)], lambda M: c**3 / (G * M))
_g("umlauf_wasserstoff_hartree", "umlauf_wasserstoff", "zweiter_weg", [()],
   lambda: K.alpha**2 * K.m_e * c**2 / hbar, toleranz=1e-10, grund="α, m_e, a_B aus CODATA (Rundung 10⁻¹¹)")
_g("umlauf_kepler_neutronenstern_handwert", "umlauf_kepler", "handwert", [(1.4 * M_SUN, 12e3)],
   lambda M, r: 1.04e4, toleranz=4.9e-3, grund="Modelltext 1,04·10⁴ rad/s, halbe letzte Stelle")
_g("umlauf_kepler_tempo", "umlauf_kepler", "zweiter_weg", [(M_SUN, K.AU), (1.4 * M_SUN, 12e3)],
   lambda M, r: c * math.sqrt(3 * (G * M / (3 * c**2)) / r) / r)   # v = c·√(3Δs/r)
_g("spinzahl_spin_halb", "spinzahl", "zweiter_weg", [(hbar / (2 * m * c), m) for m in (K.m_e, K.m_p, 1e3)],
   lambda aJ, M: 0.5 * (hbar * c / G) / M**2)   # j = ½(m_Pl/M)²
_g("spinzahl_elektron_handwert", "spinzahl", "handwert", [(hbar / (2 * K.m_e * c), K.m_e)], lambda aJ, M: 2.9e44,
   toleranz=1.8e-2, grund="Modelltext 2,9·10⁴⁴, halbe letzte Stelle")
_g("moment_magnetisch_bohr", "moment_magnetisch", "handwert", [(K.e_el, hbar / (2 * K.m_e * c))],
   lambda Q, aJ: 9.274e-24, toleranz=6e-5, grund="Bohr-Magneton 9,274·10⁻²⁴ J/T, halbe letzte Stelle")


def _einrasten(dw, K_):
    s = solve_ivp(lambda t, y: [dw - K_ * math.sin(y[0])], [0, 400 / K_], [0.0], rtol=1e-12, atol=1e-14)
    return s.y[0][-1]


_g("einrast_phase_bewegung", "einrast_phase", "zweiter_weg", [(0.3, 1.0), (-0.5, 0.8)], _einrasten, toleranz=1e-8,
   grund="Integration der Phasengleichung bis zum Stillstand")
_g("dephasierzeit_koerper_handwert", "dephasierzeit", "handwert", [(310.0,)], lambda T: 2.46e-14, toleranz=2.1e-3,
   grund="Modelltext 2,46·10⁻¹⁴ s bei 310 K, halbe letzte Stelle")

# ------------------------------------------------------------------ Licht (08)


def _ablenkung_exakt(M, b):
    """Ablenkung auf der exakten Schwarzschild-Bahn: α + π = 2∫₀^{π/2} dφ/√(1 − r_s(u₀² + u₀u + u²)/(u₀ + u)),
    u = u₀ sin φ, u₀² − r_s u₀³ = 1/b²."""
    rs = 2 * G * M / c**2
    u0 = brentq(lambda u: u * u - rs * u**3 - 1 / b**2, 0.5 / b, 1.5 / b, xtol=1e-30, rtol=1e-15)

    def f(p):   # 1/√(1 − ε) − 1, damit der flache Anteil π nicht abgezogen werden muss
        u = u0 * math.sin(p)
        return math.expm1(-0.5 * math.log1p(-rs * (u0**2 + u0 * u + u * u) / (u0 + u)))
    return 2 * quad(f, 0, math.pi / 2, epsabs=0, epsrel=1e-12)[0]


_g("ablenkung_art", "ablenkung", "grenzfall", [(M_SUN, R_SUN), (M_SUN, 1e3 * R_SUN)], _ablenkung_exakt,
   toleranz=1e-5, grund="exakte Bahn; zweite Ordnung (15π/16)·GM/(c²b) = 6·10⁻⁶ am Sonnenrand")
_g("ablenkung_sonnenrand_handwert", "ablenkung", "handwert", [(M_SUN, R_SUN)], lambda M, b: 1.75 / ARCSEK,
   toleranz=2.9e-3, grund="1,75″ am Sonnenrand, halbe letzte Stelle")
_g("ablenkung_kollektiv_rstern", "ablenkung_kollektiv", "zweiter_weg", [(M_MW, K.A0), (1.7e14 * M_SUN, K.A0)],
   lambda M, a0: 6 * math.pi * (G * M / (3 * c**2)) / math.sqrt(G * M / a0))


def _linse_sichtlinie(M, b, a0):
    """α = (2/c²)∫ a(r)·b/r dz längs der Sichtlinie, a = ν(GM/r²)·GM/r²: Quadratur in ln z, Rest jenseits Z über
    ν ≃ r/r* + ½ geschlossen."""
    def fz(z):
        r2 = b * b + z * z
        return _nu(G * M / r2, a0) * G * M * b / r2**1.5
    rs = math.sqrt(G * M / a0)
    z0, Z = 1e-6 * b, 1e9 * max(b, rs)
    innen = fz(0.0) * z0
    mitte = quad(lambda u: fz(math.exp(u)) * math.exp(u), math.log(z0), math.log(Z), limit=800, epsabs=0,
                 epsrel=1e-12)[0]
    rest = G * M * b * (1 / (rs * Z) + 1 / (4 * Z * Z))
    return 2 / c**2 * 2 * (innen + mitte + rest)


_g("ablenkung_linse_sichtlinie", "ablenkung_linse", "zweiter_weg",
   [(M_MW, x * math.sqrt(G * M_MW / K.A0), K.A0) for x in (0.1, 1.0, 10.0)], _linse_sichtlinie, toleranz=1e-7,
   grund="Quadratur längs der Sichtlinie")
_g("ablenkung_linse_I1_handwert", "ablenkung_linse", "handwert", [(M_MW, math.sqrt(G * M_MW / K.A0), K.A0)],
   lambda M, b, a0: 12 * G * M / (3 * c**2) / b * 2.135, toleranz=2.4e-4, grund="Modelltext I(1) = 2,135")
_g("zusatzpotential_quadratur", "zusatzpotential", "zweiter_weg",
   [(M_MW, x * math.sqrt(G * M_MW / K.A0), K.A0) for x in (0.3, 1.0, 10.0)],
   lambda M, r, a0: quad(lambda u: (_nu(G * M / u**2, a0) - 1) * G * M / u**2, 0, r, limit=400, epsabs=0,
                         epsrel=1e-12)[0], toleranz=1e-9)


def _b_min(M):
    """Kleinster Stoßparameter des Lichts: Minimum von b(r) = r/√(1 − 2GM/(rc²)) (Nullstelle der Ableitung)."""
    rs = 2 * G * M / c**2
    r = brentq(lambda r: 1 - rs / r - rs / (2 * r), 1.01 * rs, 10 * rs, xtol=1e-12 * rs)
    return r, r / math.sqrt(1 - rs / r)


_g("photonensphaere_bahn", "photonensphaere", "zweiter_weg", [(M_SUN,), (1e9 * M_SUN,)], lambda M: _b_min(M)[0],
   toleranz=1e-11, grund="Nullstelle mit xtol 10⁻¹²")
_g("schatten_winkel_bahn", "schatten_winkel", "zweiter_weg", [(4.3e6 * M_SUN, 8.2 * K.kpc)],
   lambda M, D: 2 * _b_min(M)[1] / D, toleranz=1e-11, grund="Nullstelle mit xtol 10⁻¹²")
_g("schatten_winkel_m87_handwert", "schatten_winkel", "handwert", [(_W["M87_M"] * M_SUN, _W["M87_D"] * K.Mpc)],
   lambda M, D: 40e-6 / ARCSEK, toleranz=1.3e-2, grund="Modelltext Ø ≈ 40 μas, halbe letzte Stelle")
_g("perihel_pro_umlauf_art", "perihel_pro_umlauf", "zweiter_weg", [(M_SUN, _W["MERKUR_A"], _W["MERKUR_E"])],
   lambda M, a, e: 6 * math.pi * G * M / (c**2 * a * (1 - e * e)))
_g("perihel_merkur_handwert", "perihel_pro_umlauf", "handwert", [(M_SUN, _W["MERKUR_A"], _W["MERKUR_E"])],
   lambda M, a, e: 43 / ARCSEK * _W["MERKUR_P"] / 36525, toleranz=1.2e-2,
   grund="Modelltext 43″ je Jahrhundert, halbe letzte Stelle")
_g("schwellenfrequenz_gw_umlauf", "schwellenfrequenz_gw", "zweiter_weg", [(M_SUN, K.A0), (M_MW, K.A0)],
   lambda M, a0: math.sqrt(G * M / math.sqrt(G * M / a0)**3) / math.pi)


def _isco_nullstelle(j, mitlaeufig=True):
    """ω_r = 0: 1 − 6/r̃ + 8j r̃^{−3/2} − 3j²/r̃² = 0 (äußerste Nullstelle)."""
    jj = j if mitlaeufig else -j
    return brentq(lambda r: 1 - 6 / r + 8 * jj * r**-1.5 - 3 * jj * jj / r**2, 1 + 1e-9 if jj > 0.99 else 1.5, 9.5,
                  xtol=1e-14)


_g("isco_kerr_nullstelle", "isco_kerr", "zweiter_weg", [(0.0, True), (0.5, True), (0.9, True), (0.9, False)],
   _isco_nullstelle, toleranz=1e-10, grund="Nullstelle mit xtol 10⁻¹⁴")
_g("hulse_taylor_peters", "hulse_taylor", "zweiter_weg",
   [(_W["PSR1913_M1"], _W["PSR1913_M2"], _W["PSR1913_PB"], _W["PSR1913_E"])],
   lambda m1, m2, Pb, e: -192 * math.pi * G**(5 / 3) / (5 * c**5) * (Pb / (2 * math.pi))**(-5 / 3)
   * m1 * m2 * M_SUN**2 * ((m1 + m2) * M_SUN)**(-1 / 3)
   * (1 + 73 / 24 * e * e + 37 / 96 * e**4) / (1 - e * e)**3.5)
_g("hulse_taylor_handwert", "hulse_taylor", "handwert",
   [(_W["PSR1913_M1"], _W["PSR1913_M2"], _W["PSR1913_PB"], _W["PSR1913_E"])], lambda m1, m2, Pb, e: -2.40e-12,
   toleranz=2.1e-3, grund="Modelltext −2,40·10⁻¹², halbe letzte Stelle")

# ------------------------------------------------------------------ Verwindung (08)
_g("verwindung_flaeche", "verwindung", "zweiter_weg", [(_W["ERDE_J"], M_EARTH)],
   lambda J, M: (G * J / c**3) / (G * M / c**2))   # GJ/c³ = 3Δs·a_J
_g("verwindung_erde_handwert", "verwindung", "handwert", [(_W["ERDE_J"], M_EARTH)], lambda J, M: 3.27,
   toleranz=1.6e-3, grund="Modelltext 3,27 m, halbe letzte Stelle")
_g("umlauf_differenz_frequenzen", "umlauf_differenz", "zweiter_weg", [(3.27,), (1e3,)],
   lambda aJ: (lambda M, rt, j: 2 * math.pi * G * M / c**3 * ((rt**1.5 + j) - (rt**1.5 - j)))(
       M_EARTH, 2e6, aJ / (G * M_EARTH / c**2)))   # 2π/ω_φ(−) − 2π/ω_φ(+) mit ω_φ = c³/(GM(r̃^{3/2} ± j)), vertauscht
_g("umlauf_differenz_erde_handwert", "umlauf_differenz", "handwert", [(3.27,)], lambda aJ: 137e-9, toleranz=3.7e-3,
   grund="Modelltext 137 ns, halbe letzte Stelle")
_g("kerr_omega_phi_kepler", "kerr_omega_phi", "grenzfall", [(M_SUN, rt, 0.0) for rt in (6.0, 1e4)],
   lambda M, rt, j: math.sqrt(G * M / (rt * G * M / c**2)**3), toleranz=1e-12, grund="j = 0: Kepler, Rundung")
_g("kerr_omega_r_schwarzschild", "kerr_omega_r", "grenzfall", [(M_SUN, rt, 0.0) for rt in (7.0, 50.0)],
   lambda M, rt, j: math.sqrt(G * M / (rt * G * M / c**2)**3) * math.sqrt(1 - 6 / rt), toleranz=1e-12,
   grund="j = 0: Schwarzschild-Epizykel, Rundung")
_g("kerr_omega_theta_schwarzschild", "kerr_omega_theta", "grenzfall", [(M_SUN, rt, 0.0) for rt in (7.0, 50.0)],
   lambda M, rt, j: math.sqrt(G * M / (rt * G * M / c**2)**3), toleranz=1e-12, grund="j = 0: ω_θ = ω_φ, Rundung")
_g("knotendrehung_frequenzen", "knotendrehung", "grenzfall",
   [(M_SUN, 1e8 * G * M_SUN / c**2, j) for j in (0.3, 0.9)],
   lambda M, r, j: -math.expm1(0.5 * math.log1p(-4 * j * (r * c**2 / (G * M))**-1.5 + 3 * j * j / (r * c**2 / (G * M))**2)),
   toleranz=1e-4, grund="führende Ordnung; nächste Ordnung ¾·j·r̃^{−1/2} ≤ 7·10⁻⁵ bei r̃ = 10⁸")
_g("nullflaeche_takt_aequator", "nullflaeche_takt", "grenzfall", [(M_SUN, 0.7 * G * M_SUN / c**2, math.pi / 2)],
   lambda M, aJ, th: 2 * G * M / c**2, toleranz=1e-12, grund="Äquator 6Δs, cos(π/2) numerisch 6·10⁻¹⁷")
_g("nullflaeche_takt_pol", "nullflaeche_takt", "grenzfall", [(M_SUN, 0.7 * G * M_SUN / c**2, 0.0)],
   lambda M, aJ, th: max(np.roots([1.0, -2 * G * M / c**2, aJ * aJ]).real), toleranz=1e-12,
   grund="Pol = äußerer Horizont, Nullstellen des Polynoms")
_g("horizont_kerr_wurzeln", "horizont_kerr", "zweiter_weg",
   [(M_SUN, a * G * M_SUN / c**2, au) for a in (0.3, 0.9) for au in (True, False)],
   lambda M, aJ, au: sorted(np.roots([1.0, -2 * G * M / c**2, aJ * aJ]).real)[1 if au else 0], toleranz=1e-12,
   grund="Nullstellen von r² − 6Δs·r + a_J² (numpy.roots)")


def _mitf_dipol(M, aJ, r, cos_psi):
    J = aJ * M * c
    return G * J / (c**2 * r**3) * math.sqrt(1 + 3 * cos_psi**2)


_g("mitfuehrung_dipol", "mitfuehrung", "zweiter_weg", [(M_EARTH, 3.27, 7.0e6, x) for x in (0.0, 0.5, 1.0)],
   _mitf_dipol)


def _mitf_bahnmittel(M, aJ, r):
    """Mittel des Vektors [3(ĵ·r̂)r̂ − ĵ] über eine polare Kreisbahn (Gauß-Legendre), Anteil längs ĵ."""
    x, w = np.polynomial.legendre.leggauss(64)
    phi = math.pi * (x + 1)
    rz = np.sin(phi)                          # ĵ·r̂ auf der polaren Bahn
    feld_z = 3 * rz * rz - 1
    return G * M / c**2 * aJ * c / r**3 * float(np.sum(w * feld_z) / 2)


_g("mitfuehrung_polarbahn_mittel", "mitfuehrung_polarbahn", "zweiter_weg", [(M_EARTH, 3.27, 7.0274e6)],
   _mitf_bahnmittel, toleranz=1e-12, grund="Gauß-Legendre, Rundung")
_g("geodaetisch_art", "geodaetisch", "zweiter_weg", [(M_EARTH, 7.0274e6), (M_SUN, K.AU)],
   lambda M, r: 1.5 * (G * M)**1.5 / (c**2 * r**2.5))

# ------------------------------------------------------------------ Wellen (08)
_g("quadrupol_leistung_peters", "quadrupol_leistung", "zweiter_weg", [(2.83 * M_SUN, 0.2495, 1e9)],
   lambda M, eta, r: 32 / 5 * G**4 * (eta * M)**2 * M**3 / (c**5 * r**5))


def _chirp(f, fdot):
    Mc = c**3 / G * (5 / 96 * math.pi**(-8 / 3) * f**(-11 / 3) * fdot)**0.6
    return G * Mc / (3 * c**2)


_g("faltungsmesser_chirp_masse", "faltungsmesser_chirp", "zweiter_weg", [(35.0, 300.0), (100.0, 20.0)], _chirp)
_g("ringdown_frequenz_masse", "ringdown_frequenz", "zweiter_weg", [(_W["GW150914_MF"] * M_SUN, _W["GW150914_F"])],
   lambda M, F: c**3 * F / (2 * math.pi * G * M))
_g("ringdown_quelle_handwert", "ringdown_frequenz", "handwert", [(_W["GW150914_MF"] * M_SUN, _W["GW150914_F"])],
   lambda M, F: 272.4, toleranz=2e-4, grund="Modelltext Quellsystem 272,4 Hz, halbe letzte Stelle")
_g("gedaechtnis_welle_energie", "gedaechtnis_welle", "zweiter_weg", [(3 * M_SUN * c**2, 410 * K.Mpc, 17 / 96)],
   lambda E, D, k: 4 * k * G * E / (c**4 * D))


def _schleife_sprung(dP_dT, G_th, dP=1.0):
    """Linearisierte Wärmebilanz: offener Sprung ΔT₀ = ΔP/G, geschlossener ΔT = ΔP/(G − ∂P/∂T); q = 1 − ΔT₀/ΔT."""
    return 1 - (dP / G_th) / (dP / (G_th - dP_dT))


_g("schleife_thermisch_sprung", "schleife_thermisch", "zweiter_weg", [(0.2, 1.0), (-3e-3, 2e-2), (0.9, 1.0)],
   _schleife_sprung, toleranz=1e-12, grund="Rundung in 1 − ΔT₀/ΔT")

# ------------------------------------------------------------------ Wellen (08): Abklingzeit und Richtungsfaktor
_g("ringdown_abklingzeit_daempfung", "ringdown_abklingzeit", "zweiter_weg",
   [(_W["GW150914_MF"] * M_SUN, _W["GW150914_F"], 3.1669), (10 * M_SUN, 0.3737, 2.1)],
   lambda M, F, Q: (G * M / c**3) / (F / (2 * Q)))   # τ = (GM/c³)/(M·ω_I), M·ω_I = F/(2Q)
_g("richtungsfaktor_gedaechtnis_kosinus", "richtungsfaktor_gedaechtnis", "zweiter_weg",
   [(x,) for x in (0.3, 1.0, math.pi / 2, 2.6)],
   lambda i: (lambda u: (1 - u * u) * (17 + u * u) / 96)(math.cos(i)), toleranz=1e-12, grund="Rundung")
_g("richtungsfaktor_gedaechtnis_kante", "richtungsfaktor_gedaechtnis", "grenzfall", [(math.pi / 2,)],
   lambda i: 17 / 96, toleranz=1e-12, grund="Kante: sin² = 1, cos² = 0 (Rundung von cos(π/2))")

# ------------------------------------------------------------------ Falte und Syntropie (04)
_g("entropie_falte_alpha_G", "entropie_falte", "zweiter_weg", [(K.m_p,), (M_SUN,), (1e12,)],
   lambda M: 4 * math.pi * G * M * M / (hbar * c))


def _entropie_smarr(M):
    """Smarr: Mc² = 2·T_H·S mit T_H = ħc³/(8πGMk_B), S in k_B."""
    T_H = hbar * c**3 / (8 * math.pi * G * M * kB)
    return M * c * c / (2 * kB * T_H)


_g("entropie_falte_smarr", "entropie_falte", "zweiter_weg", [(M_SUN,), (1e12,)], _entropie_smarr, toleranz=1e-12,
   grund="Rundung")


def _page_gleichstand(beta):
    """Gleichstand S_rad = S_F mit S_F/S₀ = (1 − x)^{2/3}, S_rad = β(S₀ − S_F), numerisch."""
    return brentq(lambda x: (1 - x)**(2 / 3) - beta * (1 - (1 - x)**(2 / 3)), 1e-12, 1 - 1e-12, xtol=1e-15)


_g("page_zeit_gleichstand", "page_zeit", "zweiter_weg", [(1.0,), (1.4847,), (3.0,)], _page_gleichstand,
   toleranz=1e-10, grund="Nullstelle (brentq, xtol 10⁻¹⁵)")
_g("rueckverdichtung_teile_bilanz", "rueckverdichtung_teile", "zweiter_weg", [(1.0, 1.4847), (10.0, 1.4847), (0.3, 2.0)],
   lambda z, b: brentq(lambda N: (N / z)**2 - N * b, 0.5 * b * z * z, 2 * b * z * z, xtol=1e-14),
   toleranz=1e-10, grund="Nullstelle der Entropiebilanz (brentq)")


def _stefan_boltzmann_planck():
    """σ_SB aus dem Planck-Spektrum: σ = (2πk⁴/(h³c²))·∫x³/(eˣ − 1)dx."""
    return 2 * math.pi * kB**4 / ((2 * math.pi * hbar)**3 * c * c) * quad(lambda x: x**3 / math.expm1(x), 0, 200)[0]


_g("temperatur_abstrahlung_planck", "temperatur_abstrahlung", "zweiter_weg", [(1.215e17, K.R_earth), (1.0, 1.0)],
   lambda P, R: (P / (4 * math.pi * R * R * _stefan_boltzmann_planck()))**0.25, toleranz=1e-10,
   grund="Quadratur des Planck-Spektrums")


def _entropie_je_energie(T):
    """S/E der Hohlraumstrahlung aus der Bose-Besetzung n: s ∝ ∫x²[(1+n)ln(1+n) − n ln n]dx, u ∝ kT·∫x³n dx."""
    def s(x):
        n = 1 / math.expm1(x)
        return x * x * ((1 + n) * math.log1p(n) - n * math.log(n))
    return quad(s, 1e-9, 200, limit=200)[0] / (kB * T * quad(lambda x: x**3 / math.expm1(x), 0, 200)[0]) * kB


_g("entropieexport_bose", "entropieexport", "zweiter_weg", [(1.215e17, 255.0, 5772.0), (1.0, 3.0, 300.0)],
   lambda P, Ta, Te: P * (_entropie_je_energie(Ta) - _entropie_je_energie(Te)), toleranz=1e-8,
   grund="Quadratur der Entropie der Bose-Besetzung")

# ------------------------------------------------------------------ Starke Kopplung und Ladung (02)
_SCHWELLEN = (1.27, 4.18, 172.7)


def _lauf_eine_schleife(a0, mu0, mu1, schwellen=(), schleifen=1):
    """Eine Schleife geschlossen: 1/α(μ₁) = 1/α(μ₀) + Σ β₀(n_f)/(2π)·Δln μ über die Abschnitte."""
    L0, L1 = math.log(mu0), math.log(mu1)
    kanten = sorted({L0, L1} | {math.log(m) for m in schwellen if min(mu0, mu1) < m < max(mu0, mu1)})
    inv = 1 / a0
    for a, b in zip(kanten[:-1], kanten[1:]):
        nf = 3 + sum(1 for m in schwellen if math.log(m) < (a + b) / 2)
        inv += (11 - 2 * nf / 3) / (2 * math.pi) * (b - a) * (1 if L1 > L0 else -1)
    return 1 / inv


_g("alpha_s_lauf_geschlossen", "alpha_s_lauf", "zweiter_weg",
   [(0.118, 91.1876, mu, _SCHWELLEN, 1) for mu in (1.5, 10.0, 1e3, 1.22e19)], _lauf_eine_schleife, toleranz=1e-9,
   grund="Integration (rtol 10⁻¹¹) gegen die geschlossene Form einer Schleife")


def _skala_eine_schleife(a0, mu0, a_ziel, schwellen=(), schleifen=1):
    """Eine Schleife geschlossen abwärts: je Abschnitt Δln μ = −(2π/β₀)(1/α_Ziel − 1/α) bis zum Ziel."""
    L, inv = math.log(mu0), 1 / a0
    for m in sorted([m for m in schwellen if m < mu0], reverse=True) + [0.0]:
        nf = 3 + sum(1 for x in schwellen if x < math.exp(L) * (1 - 1e-12))
        b = (11 - 2 * nf / 3) / (2 * math.pi)
        L_ziel = L - (inv - 1 / a_ziel) / b
        if m == 0.0 or L_ziel > math.log(m):
            return math.exp(L_ziel)
        inv -= b * (L - math.log(m))
        L = math.log(m)


_g("skala_bei_kopplung_geschlossen", "skala_bei_kopplung", "zweiter_weg",
   [(0.02, 1.22e19, 0.5, _SCHWELLEN, 1), (0.118, 91.1876, 0.3, _SCHWELLEN, 1)], _skala_eine_schleife, toleranz=1e-8,
   grund="Ereignis der Integration (rtol 10⁻¹¹) gegen die geschlossene Form einer Schleife")


def _nome_pol(alpha_s, nf):
    """Λ/M ist die Skala, an der die eine Schleife 1/α(μ) = 1/α_s + (β₀/2π)ln(μ/M) null wird (Nullstelle)."""
    b0 = (11 * 3 - 2 * nf) / 3
    return math.exp(brentq(lambda L: 1 / alpha_s + b0 / (2 * math.pi) * L, -1e4, 0.0, xtol=1e-14))


_g("nome_qcd_pol", "nome_qcd", "zweiter_weg", [(0.1, 3), (0.3, 5), (0.0189, 6)], _nome_pol, toleranz=1e-9,
   grund="Nullstelle (brentq, xtol 10⁻¹⁴ in ln μ)")
_g("dirac_verhaeltnis_si", "dirac_verhaeltnis", "zweiter_weg", [()],
   lambda: 2 * math.pi * hbar * K.eps0 * c / K.e_el**2, toleranz=1e-9,
   grund="h/(μ₀ce²) mit μ₀ = 1/(ε₀c²): ε₀, e und α aus CODATA stimmen auf 10⁻¹⁰")
_g("radius_ladungskreis_kaluza_klein", "radius_ladungskreis", "zweiter_weg", [()],
   lambda: math.sqrt(16 * math.pi * G * K.eps0) * hbar / (c * K.e_el), toleranz=1e-9,
   grund="Kaluza-Klein q₁ = √(16πGε₀)ħ/(cR) = e; ε₀, e und α aus CODATA stimmen auf 10⁻¹⁰")
_g("energie_ladungskreis_planck", "energie_ladungskreis", "zweiter_weg", [()],
   lambda: math.sqrt(K.alpha) / 2 * math.sqrt(hbar * c / G) * c * c)

# ------------------------------------------------------------------ Zustand und Sterne (02)
def _zustand_maxwell(T, m):
    """2⟨v_x²⟩/c² aus der Maxwell-Verteilung e^{−mv²/2kT} (Quadratur, v in Einheiten √(kT/m) nur zur Skalierung)."""
    v0 = math.sqrt(kB * T / m)
    z = quad(lambda v: v * v * math.exp(-m * v * v / (2 * kB * T)), -40 * v0, 40 * v0, points=[0.0])[0]
    n = quad(lambda v: math.exp(-m * v * v / (2 * kB * T)), -40 * v0, 40 * v0, points=[0.0])[0]
    return 2 * z / n / c**2


_g("zustandsschwelle_maxwell", "zustandsschwelle", "zweiter_weg", [(1811.0, 55.845 * K.m_u), (300.0, K.m_p)],
   _zustand_maxwell, toleranz=1e-9, grund="Quadratur der Maxwell-Verteilung")
_g("f_saha_hartree", "f_saha", "zweiter_weg", [(5772.0,), (3e4,)],
   lambda T: kB * T * K.m_e * K.a_B**2 / hbar**2, toleranz=1e-9,
   grund="α²m_ec² = ħ²/(m_ea_B²); a_B und α aus CODATA stimmen auf 10⁻¹⁰")


def _eddington_druck(beta, mu):
    """Polytrope n = 3 mit festem β: P = K ρ^{4/3}, K = [(k/(μm_u))⁴·3(1 − β)/(aβ⁴)]^{1/3}, a = π²k⁴/(15ħ³c³);
    M = 4π·m₁·(K/(πG))^{3/2} (eigene Lane-Emden-Lösung)."""
    a_rad = math.pi**2 * kB**4 / (15 * hbar**3 * c**3)
    Kp = ((kB / (mu * K.m_u))**4 * 3 * (1 - beta) / (a_rad * beta**4))**(1 / 3)
    return 4 * math.pi * _polytrope(3.0)[1] * (Kp / (math.pi * G))**1.5


_g("eddington_quartik_druck", "eddington_quartik", "zweiter_weg", [(0.5, 0.62), (0.9, 0.6), (0.99, 1.0)],
   _eddington_druck, toleranz=1e-8, grund="eigene Integration der Polytrope (rtol 10⁻¹²)")
_g("eddington_quartik_chandrasekhar", "eddington_quartik", "handwert", [(0.5, 0.62), (0.9, 0.6)],
   lambda beta, mu: M_SUN * math.sqrt((1 - beta) / (0.002979 * mu**4 * beta**4)), toleranz=2e-3,
   grund="Tafelwert 1 − β = 0,002979·(M/M☉)²μ⁴β⁴ (Chandrasekhar), vier Stellen; hier 0,002990")

# ------------------------------------------------------------------ Kerne (06)
def _troepfchen_normalgleichung(A_min):
    """Dieselbe Anpassung über die Normalgleichungen (XᵀX)a = XᵀB, Tafel eigens gelesen."""
    X, B = [], []
    for Z in range(0, 121):
        for A in range(max(int(math.ceil(A_min)), 1), 301):
            try:
                k = K.ame2020(Z, A)
            except KeyError:
                continue
            if not k["gemessen"] or not math.isfinite(k["BA"]):
                continue
            N = A - Z
            d = 1.0 if Z % 2 == 0 and N % 2 == 0 else (-1.0 if Z % 2 and N % 2 else 0.0)
            X.append([A, -A**(2 / 3), -Z * (Z - 1) / A**(1 / 3), -(N - Z)**2 / A, d / math.sqrt(A)])
            B.append(k["BA"] * A)
    X, B = np.array(X), np.array(B)
    return tuple(np.linalg.solve(X.T @ X, X.T @ B))


_g("troepfchen_anpassung_normalgleichung", "troepfchen_anpassung", "zweiter_weg", [(16,)], _troepfchen_normalgleichung,
   toleranz=1e-8, grund="Normalgleichungen gegen Kleinste-Quadrate-Löser (Kondition)")
_g("spaltbarkeit_energien", "spaltbarkeit", "zweiter_weg", [(92, 238, 1.1e-13, 2.7e-12), (26, 56, 1.1e-13, 2.7e-12)],
   lambda Z, A, aC, aS: (aC * Z * Z / A**(1 / 3)) / (2 * aS * A**(2 / 3)), toleranz=1e-12, grund="Rundung")

# ------------------------------------------------------------------ Kohärenz (09)
def _debye_bjerrum(T, eps_r, I):
    """λ_D = (8π·ℓ_B·N_A·I)^{−1/2} mit der Bjerrum-Länge ℓ_B = e²/(4πε₀ε_rk_BT)."""
    l_B = K.e_el**2 / (4 * math.pi * K.eps0 * eps_r * kB * T)
    return (8 * math.pi * l_B * K.N_A * I)**-0.5


_g("debye_laenge_bjerrum", "debye_laenge", "zweiter_weg", [(310.0, 74.2, 150.0), (298.15, 78.5, 1.0)], _debye_bjerrum,
   toleranz=1e-12, grund="Rundung")
_g("debye_laenge_handwert", "debye_laenge", "handwert", [(298.15, 78.5, 150.0)], lambda T, e, I: 0.304e-9 / math.sqrt(0.15),
   toleranz=2e-3, grund="Faustformel 0,304 nm/√(I/(mol/l)) bei 25 °C, drei Stellen")
