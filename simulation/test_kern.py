# -*- coding: utf-8 -*-
"""Schnelle Prüfungen des Kerns gegen Handrechnungen (Aufruf: python -m pytest -q simulation)."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kern as K  # noqa: E402


def test_faltung_sonne_von_hand():
    assert abs(K.faltung(K.M_sun) - 1.327124400e20 / (3 * 299792458.0**2)) < 1e-6


def test_takt_am_grenzradius_null():
    assert K.takt(K.M_sun, K.grenzradius(K.M_sun)) == 0.0


def test_regel_grenzfaelle():
    assert abs(K.verstaerkung(1e4 * K.a0) - 1) < 1e-40       # Newton
    aN = 1e-6 * K.a0                                          # tief: a → √(a_N a₀)
    assert abs(K.verstaerkung(aN) * aN / math.sqrt(aN * K.a0) - 1) < 1e-3
    assert abs(K.verstaerkung(K.a0) - 1 / (1 - math.exp(-1))) < 1e-12


def test_schwelle_und_leiter():
    M = K.M_sun
    assert abs(K.schwellenradius(M)**2 / (3 * K.faltung(M) * 2 * K.schliesslaenge()) - 1) < 1e-12
    assert abs(K.rueckkehrzeit() / (K.c / K.a0) - 1) < 1e-12


def test_tully_fisher_von_hand():
    M = 6e10 * K.M_sun
    assert abs(K.v_flach(M) / (6e10 * 1.327124400e20 * 1.06e-10)**0.25 - 1) < 1e-9


def test_gehaltener_trichter_spreizung():
    """Hülle R = λ_H·N^{1/3} auf der Massenskala: δ = GM_*/(Rc²) = λ̄_p/λ_H (Weißer Zwerg m_e/m_p)."""
    M = K.sternmasse()
    R = K.hbar / (K.m_e * K.c) * (M / K.m_p)**(1 / 3)
    assert abs(K.G * M / (R * K.c**2) / (9.1093837015e-31 / 1.67262192369e-27) - 1) < 1e-9

def test_spiegelmasse_quantenlaenge_ist_grenzradius():
    for M in (K.m_e, K.M_earth, K.M_sun, 1e12 * K.M_sun):
        assert abs(K.compton(K.spiegelmasse(M)) / K.grenzradius(M) - 1) < 1e-12
        assert abs(K.spiegelmasse(M) * M / (K.m_Pl**2 / 2) - 1) < 1e-9


def test_modulgruppe_spiegel_und_dreier():
    """PSL(2,ℤ) = ℤ₂ ∗ ℤ₃: S² = −I, (ST)³ = −I; Fixpunkte i (S) und e^{2πi/3} (ST)."""
    import numpy as np
    S, T = np.array([[0, -1], [1, 0]]), np.array([[1, 1], [0, 1]])
    assert (S @ S == -np.eye(2)).all()
    ST = S @ T
    assert (ST @ ST @ ST == -np.eye(2)).all()
    wirk = lambda m, z: (m[0, 0] * z + m[0, 1]) / (m[1, 0] * z + m[1, 1])
    rho = complex(-0.5, math.sqrt(3) / 2)
    assert abs(wirk(S, 1j) - 1j) < 1e-12 and abs(wirk(ST, rho) - rho) < 1e-12

def test_regel_gleichgewicht_nullpunkt_und_schleifenspiegel():
    """ν = 1 + n_B, q = n_B/ν, ν(s) + ν(−s) = 1, ν = ½ + ½coth(s/2); ganze Schleifensumme Z = ∏1/(1 − e^{−ns}):
    η(i/t) = √t·η(it) mit s = 2πt, also s ↔ 4π²/s, und ln Z → π²/(6s) + ½ln(s/2π) − s/24 für kleine s."""
    import mpmath as mp
    for s in (0.3, 1.0, 2.5, 7.0):
        nB, nu = 1 / math.expm1(s), 1 / (1 - math.exp(-s))
        assert abs(nB / nu - math.exp(-s)) < 1e-12 and abs(nu + 1 / (1 - math.exp(s)) - 1) < 1e-12
        assert abs(nu - (0.5 + 0.5 / math.tanh(s / 2))) < 1e-12
        eta = lambda t: mp.e**(-mp.pi * t / 12) * mp.nprod(lambda n: 1 - mp.e**(-2 * mp.pi * n * t), [1, mp.inf])
        t = mp.mpf(s) / (2 * mp.pi)
        assert abs(eta(1 / t) / (mp.sqrt(t) * eta(t)) - 1) < 1e-12
    s = mp.mpf("0.05")
    lnZ = -mp.nsum(lambda n: mp.log(1 - mp.e**(-n * s)), [1, mp.inf])
    assert abs(lnZ - (mp.pi**2 / (6 * s) + mp.log(s / (2 * mp.pi)) / 2 - s / 24)) < 1e-12

def test_zeta_spiegel_ganz_und_oertlich():
    """n_B(s) = 1/s + Σ ζ(−n)(−s)ⁿ/n! (örtlich), ∫ x^d n_B = Γ(d+1)ζ(d+1) (ganz); ζ(1−k) = 2(2π)^{−k}cos(πk/2)Γ(k)ζ(k)."""
    import mpmath as mp
    s = mp.mpf("0.3")
    reihe = 1 / s + mp.nsum(lambda n: mp.zeta(-n) * (-s)**n / mp.factorial(n), [0, 30])
    assert abs(reihe - 1 / mp.expm1(s)) < 1e-12
    for d in (1, 2, 3):
        k = d + 1
        assert abs(mp.quad(lambda x: x**d / mp.expm1(x), [0, mp.inf]) - mp.gamma(k) * mp.zeta(k)) < 1e-12
        assert abs(mp.zeta(1 - k) - 2 * (2 * mp.pi)**(-k) * mp.cos(mp.pi * k / 2) * mp.gamma(k) * mp.zeta(k)) < 1e-12
    assert abs(mp.zeta(-3) - mp.mpf(1) / 120) < 1e-15 and abs(mp.gamma(4) * mp.zeta(4) - mp.pi**4 / 15) < 1e-15


def test_drei_arten_der_torussymmetrie():
    """Spur < 2 Fixpunkt innen (Spiegel i, Dreier ρ), = 2 Rand (voller Umlauf), > 2 durchlaufend (golden: φ)."""
    import numpy as np
    phi = (1 + 5**0.5) / 2
    M = np.array([[2, 1], [1, 1]])
    assert abs((M[0, 0] * phi + M[0, 1]) / (M[1, 0] * phi + M[1, 1]) - phi) < 1e-12
    spuren = [abs(np.trace(m)) for m in (np.array([[0, -1], [1, 0]]), np.array([[0, -1], [1, 1]]), np.array([[1, 1], [0, 1]]), M)]
    assert spuren == [0, 1, 2, 3]


def test_netz_auf_dem_geschlossenen_raum():
    """Morse auf T³ (Euler-Charakteristik 0): Trichter − Fadensättel + Wandsättel − Leerräume = 0."""
    import itertools
    import numpy as np
    from scipy.optimize import fsolve
    rng = np.random.default_rng(5)
    K = np.array([k for k in itertools.product(range(-1, 2), repeat=3) if k > (0, 0, 0)])
    A, P = rng.normal(size=len(K)), rng.uniform(0, 2 * np.pi, len(K))
    grad = lambda x: -(A * np.sin(K @ x + P)) @ K
    hess = lambda x: -(K.T * (A * np.cos(K @ x + P))) @ K
    pkt = []
    for x0 in rng.uniform(0, 2 * np.pi, (6000, 3)):
        x, _, ier, _ = fsolve(grad, x0, full_output=True)
        if ier != 1 or np.linalg.norm(grad(x)) > 1e-9:
            continue
        x = np.mod(x, 2 * np.pi)
        if all(np.linalg.norm((x - q + np.pi) % (2 * np.pi) - np.pi) > 1e-5 for q in pkt):
            pkt.append(x)
    n = [sum(int((np.linalg.eigvalsh(hess(x)) < 0).sum()) == i for x in pkt) for i in range(4)]
    assert min(n) > 0 and n[0] - n[1] + n[2] - n[3] == 0

def test_kepler_grosskreis_und_kraftspiegel():
    """Fock: die Impulse einer Keplerbahn (Kreis, Hamilton) liegen auf der Einheitskugel auf einem Großkreis.
    Bohlin: w = z² bildet die Schwingerbahn (r¹) auf eine Keplerellipse (r⁻²) mit Brennpunkt im Ursprung ab.
    Newton–Kasner: r^a ↔ r^A mit (a + 3)(A + 3) = 4, Spiegelpunkt a = −1 (tiefer Bereich der Regel)."""
    import numpy as np
    e, h = 0.6, 1.0
    th = np.linspace(0, 2 * np.pi, 50)
    p = np.stack([-np.sin(th), e + np.cos(th)], 1) / h
    p0 = np.sqrt((1 - e**2) / h**2)                       # p0² = −2E
    u = np.column_stack([2 * p0 * p, (p**2).sum(1) - p0**2]) / ((p**2).sum(1) + p0**2)[:, None]
    assert np.allclose((u**2).sum(1), 1) and np.linalg.svd(u, compute_uv=False)[2] < 1e-10
    a, b = 2.0, 0.7
    t = np.linspace(0, 2 * np.pi, 50)
    w = (a * np.cos(t) + 1j * b * np.sin(t))**2
    assert np.allclose(abs(w) + abs(w - (a**2 - b**2)), a**2 + b**2)
    dual = lambda x: 4 / (x + 3) - 3
    assert dual(-2) == 1 and dual(1) == -2 and dual(-1) == -1

def test_vakuumteil_als_spiegel():
    """ρ_s c² = (π³/30)·√(c⁷/(ħG²) · ħc/ℓ₀⁴), ℓ₀ = 2Λ: geometrisches Mittel aus Planck-Dichte und Quantendichte der Schließlänge."""
    l0 = 2 * K.schliesslaenge()
    mittel = math.sqrt(K.c**7 / (K.hbar * K.G**2) * K.hbar * K.c / l0**4)
    assert abs(K.dichte_vakuumteil() * K.c**2 / mittel / (math.pi**3 / 30) - 1) < 1e-12

def test_kleinstes_netz_betti():
    """Φ = cos x + cos y + cos z auf T³: kritische Punkte bei x, y, z ∈ {0, π}, Index = Zahl negativer Krümmungen:
    1 Trichter, 3 Fäden, 3 Wände, 1 Leerraum = Betti-Zahlen (1, 3, 3, 1); Φ → −Φ tauscht Index k und 3 − k."""
    import itertools
    n = [0, 0, 0, 0]
    for p in itertools.product((0.0, math.pi), repeat=3):
        n[sum(-math.cos(x) < 0 for x in p)] += 1
    assert n == [1, 3, 3, 1] and n == n[::-1] and n == [math.comb(3, k) for k in range(4)]

def test_bericht_wird_nie_ueberschrieben(tmp_path):
    import lauf
    (tmp_path / "bericht_0001.md").write_text("alt")
    (tmp_path / "bericht_0007.md").write_text("alt")
    assert lauf.naechste_nummer(str(tmp_path)) == 8


def test_falte_spiegel_und_temperatur():
    """Falte: r₀ = λ̄(m̃), r₀(m̃) = λ̄(M); kT = m̃c²/4π ist die Hawking-Temperatur (6,17·10⁻⁸ K je Sonnenmasse)."""
    M = 7.3e35
    mt = K.spiegelmasse(M)
    assert abs(K.hbar / (mt * K.c) / (2 * K.G * M / K.c**2) - 1) < 1e-12
    assert abs(2 * K.G * mt / K.c**2 / (K.hbar / (M * K.c)) - 1) < 1e-12
    assert abs(K.falte_temperatur(K.M_sun) / 6.17e-8 - 1) < 1e-3
    assert abs(K.falte_temperatur(M) * 24 * math.pi * K.k_B * K.faltung(M) / (K.hbar * K.c) - 1) < 1e-12


def test_schallhorizont_zweiter_weg():
    """r_d unabhängig von der Anpassung: Rekombination (Peebles, Faktor 1,14), Drag-Tiefe τ = 1, Integral c_s/H."""
    import numpy as np
    from scipy.integrate import solve_ivp, quad
    h, me, mH, sT, eV, Y = 6.62607015e-34, K.m_e, 1.6735575e-27, 6.6524587e-29, K.e, 0.245
    rc = K.rho_krit100
    rr, _, rm = K.dichte_heute()
    rg = math.pi**2 / 15 * (K.k_B * K.T_CMB)**4 / (K.hbar * K.c)**3 / K.c**2
    for Q in (None, 2 * math.pi - 1):
        wc = K.omega_c if Q is None else Q * K.omega_b
        rho_m = (K.omega_b + wc) * rc
        H = lambda z: math.sqrt(8 * math.pi * K.G / 3 * (rr * (1 + z)**4 + rho_m * (1 + z)**3))
        nH0 = (1 - Y) * K.omega_b * rc / mH
        aB = lambda T: 1.14e-19 * 4.309 * (T / 1e4)**-0.6166 / (1 + 0.6703 * (T / 1e4)**0.53)
        bB = lambda T: aB(T) * (2 * math.pi * me * K.k_B * T / h**2)**1.5 * math.exp(-3.4 * eV / (K.k_B * T))

        def dx(z, x):
            T, n, Hz = K.T_CMB * (1 + z), nH0 * (1 + z)**3, H(z)
            k = 121.567e-9**3 / (8 * math.pi * Hz)
            C = (1 + k * 8.22458 * n * (1 - x[0])) / (1 + k * (8.22458 + bB(T)) * n * (1 - x[0]))
            return [C * (aB(T) * n * x[0]**2 - bB(T) * (1 - x[0]) * math.exp(-10.2 * eV / (K.k_B * T))) / (Hz * (1 + z))]
        z0 = 1700
        T0 = K.T_CMB * (1 + z0)
        S = (2 * math.pi * me * K.k_B * T0 / h**2)**1.5 * math.exp(-13.6057 * eV / (K.k_B * T0)) / (nH0 * (1 + z0)**3)
        sol = solve_ivp(dx, [z0, 200], [(-S + math.sqrt(S * S + 4 * S)) / 2], method="LSODA", dense_output=True, rtol=1e-8, atol=1e-12)
        R = lambda z: 3 * K.omega_b * rc / (4 * rg) / (1 + z)
        zs = np.linspace(200, z0, 15001)
        f = np.array([nH0 * (1 + z)**3 * sol.sol(z)[0] * sT * K.c / ((1 + z) * H(z)) / R(z) for z in zs])
        tau = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(zs))])
        zd = np.interp(1.0, tau, zs)
        rd = quad(lambda z: K.c / math.sqrt(3 * (1 + R(z))) / H(z), zd, np.inf, limit=200)[0] / K.Mpc
        assert abs(rd / K.schallhorizont(Q) - 1) < 1e-3


def test_elemente_zweiter_weg():
    """Füllfolge von Hand; Kippstelle im Tal kontinuierlich und ganzzahlig (zwei Verfahren derselben Simulation) stimmen überein."""
    import kerne as KE
    folge = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p", "8s", "5g", "6f", "7d", "8p"]
    platz = {"s": 2, "p": 6, "d": 10, "f": 14, "g": 18}
    z, ende = 0, {}
    for o in folge:
        ende[o] = (z + 1, z + platz[o[1]])
        z += platz[o[1]]
    assert ende["7p"][1] == 118 and ende["8s"] == (119, 120) and ende["5g"] == (121, 138) and ende["8p"][1] == 168
    assert tuple(e[1] for o, e in ende.items() if o[1] == "p" or o == "1s") == K.schliesszahlen_huelle(8)
    for f in (KE.TROPFEN, KE.SCHALE):
        zk, ak = KE.kippstelle(f)
        zg, a1, a2 = KE.kippstelle_ganz(f)
        assert abs(zk - zg) <= 1.5 and a1 - 10 <= ak <= a2 + 10


def test_systemnamen():
    """Systematische IUPAC-Namen der neuen Elemente."""
    import elemente as E
    soll = {119: ("Ununennium", "Uue"), 120: ("Unbinilium", "Ubn"), 121: ("Unbiunium", "Ubu"), 122: ("Unbibium", "Ubb"),
            123: ("Unbitrium", "Ubt"), 124: ("Unbiquadium", "Ubq"), 125: ("Unbipentium", "Ubp"), 126: ("Unbihexium", "Ubh"),
            127: ("Unbiseptium", "Ubs"), 128: ("Unbioctium", "Ubo"), 129: ("Unbiennium", "Ube"), 130: ("Untrinilium", "Utn")}
    assert all(E.systemname(z) == v for z, v in soll.items())


def test_dirac_1s():
    """Dirac 1s: Punktkern gibt die geschlossene Formel; Grenze 6c⁴·Δs_e/G mit dem Kernradius des Modells auf zweitem Weg (Schuss bei E = −mc²)."""
    from scipy.integrate import solve_ivp
    assert abs(K.z_kritisch(lambda Z: 1e-6, E=0.5, lo=100, hi=130) - math.sqrt(0.75) / K.alpha) < 0.01
    lam_e = K.hbar / (K.m_e * K.c) * 1e15
    import kerne as KE
    lam_H = KE.halterlaenge(KE.a_C(KE.TROPFEN))

    def knoten(Z):
        R = lam_H * KE.a_tal(Z, KE.TROPFEN)**(1 / 3) / lam_e
        z, E = Z * K.alpha, -1.0 + 1e-6
        V = lambda r: -z * (3 - r * r / (R * R)) / (2 * R) if r < R else -z / r
        s = solve_ivp(lambda r, y: [y[0] / r + (E - V(r) + 1) * y[1], -y[1] / r - (E - V(r) - 1) * y[0]],
                      [R * 1e-4, 300], [R * 1e-4, 0.0], rtol=1e-9, atol=1e-30, method="DOP853")
        G = s.y[0]
        return int((G[1:] * G[:-1] < 0).sum())
    zk = KE.huellen_grenze(KE.TROPFEN)
    assert knoten(zk - 0.5) == 0 and knoten(zk + 0.5) >= 1


def test_alpha_lebensdauer_u238():
    """Tunnelform des Alpha-Zerfalls trifft U-238 (gemessen 4,463 Gy) auf besser als 0,2 Dekaden."""
    import kerne as KE
    assert abs(KE.alpha_lebensdauer()[3]) < 0.2 and KE.alpha_lebensdauer()[1] < 0.6


def test_takt_aus_dem_tempo_gibt_dirac():
    """1s und 2s im Trichter −Z/r mit dem Takt aus dem eigenen Tempo gegen die Dirac-Energie von Hand, Z = 80."""
    import numpy as np
    import huelle_2_entstehung as E
    from huelle_1_loeser import _gitter, c
    r, h = _gitter()
    Z = 80
    for n in (1, 2):
        e, _ = E._bahn(Z, -Z / r, r, h, (n, 0))
        g = math.sqrt(1 - (Z / c)**2)
        dirac = c * c * (1 / math.sqrt(1 + (Z / c)**2 / (n - 1 + g)**2) - 1)
        assert abs(e / dirac - 1) < 1e-4


def test_magnetismus_windung_mal_verwindung():
    """μ = e·c·λ̄/2 ist eħ/(2m_e) von Hand; j = 1 bei m_Pl/√2; Langevin steckt in der Regel; Spin ½ ist die Fermi-Form;
    Kippstelle: m ≈ √(3(1 − t)), h_c ≈ (2/3)(1 − t)^{3/2}."""
    assert abs(K.magneton(K.m_e) / (K.e * K.hbar / (2 * K.m_e)) - 1) < 1e-12
    assert abs(K.spinzahl_teilchen(K.m_Pl / math.sqrt(2)) - 1) < 1e-12
    for x in (0.1, 1.0, 5.0):
        assert abs(K.verstaerkung(K.a0 * x * x) - (0.5 + 1 / x + 0.5 * K.langevin(x / 2))) < 1e-10
        assert abs(K.ausrichtung_windung(x) - (1 - 2 / (math.exp(2 * x) + 1))) < 1e-12
    assert abs(K.einrasten_magnet(0.9999) / math.sqrt(3 * 1e-4) - 1) < 1e-3
    assert abs(K.hysteresefenster(0.9999) / (2 / 3 * 1e-4**1.5) - 1) < 1e-3
    assert K.einrasten_magnet(1.5) < 1e-6                     # über der Kippstelle keine Ausrichtung ohne Feld



def test_schallwinkel_liest_das_atmen():
    """Mit dem gemessenen ω_c liest der Schallwinkel dieselbe Schwelle wie das Atmen nach Planck; z_* und r_* liegen
    bei den Werten der Standardrechnung (z_* ≈ 1090, r_* ≈ 144 Mpc)."""
    assert abs(K.a0_aus_schallwinkel() / K.a0_aus_atmen() - 1) < 1e-6
    assert abs(K.z_entkopplung() / K.z_stern - 1) < 3e-3
    assert 143 < K.schallhorizont_stern(lambda a: K.dichte_de_planck()) / K.Mpc < 146


def test_windung_gegen_faltung():
    """An der Grenze M = m_p·(α/α_G)^{3/2} ist die Bindung je Nukleon aus der Faltung gleich der aus der Windung, wenn Atome
    den Körper halten (R = a_B·N^{1/3}); seine Tiefe ist dort α²·m_e/m_p."""
    M, d = K.masse_windung_faltung()
    N = M / K.m_p
    R = K.bohr_radius() * N**(1 / 3)
    assert abs(K.G * M * K.m_p / R / (K.alpha * K.hbar * K.c / K.bohr_radius()) - 1) < 1e-9
    assert abs(K.trichter(M, R)["tiefe"] / d - 1) < 1e-9


def test_offene_windung_faellt():
    """Die offene Windung des Wasserstoffs, Umlauf für Umlauf mit dem Verlust (8π/3)(r_e/r)(v/c) je Umlauf geführt, fällt
    in derselben Zeit wie nach der geschlossenen Formel a_B³/(4r_e²c)."""
    import math
    r_e = K.coulomb_konstante() * K.e**2 / (K.m_e * K.c**2)
    r, t = K.bohr_radius(), 0.0
    while r > 1e-3 * K.bohr_radius():
        v = math.sqrt(K.coulomb_konstante() * K.e**2 / (K.m_e * r))
        f = K.verlust_windung(r, v)
        n = max(1, int(1e-3 / f))
        t += n * 2 * math.pi * r / v
        r *= (1 - f)**n
    assert abs(t / K.fall_offene_windung() - 1) < 5e-3


def test_halo_des_standardmodells():
    """NFW mit c = 10: größte Geschwindigkeit bei c·x = 2,16258 mit 1,2051·V₂₀₀ (von Hand); dort ist der Rest null."""
    import numpy as np
    import daten as D
    r200 = 100 / (10 * 70 / 1e3)                              # V₂₀₀ = 100 km/s, H₀ = 70: R₂₀₀ [kpc]
    null = np.zeros(1)
    k = dict(R=np.array([0.216258 * r200]), V=np.array([120.51]), eV=np.ones(1), Vgas=null, Vdisk=null, Vbul=null)
    d, _ = D._reste([0, 0, 0, 2.0, 1.0], None, None, dict(Inc=90.0, D=10.0), k, 0.1, 1.0, K.kpc, 0.5, halo=70.0, roh=True)
    assert abs(d[0]) < 1e-4


def test_schalen_des_gehaltenen_trichters():
    """Abstand der Schalen bei λ_H = 1,2 fm von Hand: (5/6)(3/2)^{4/3}·(197,327 MeV fm)²/(938,27 MeV·1,44 fm²) = 41,2 MeV;
    reiner Oszillator: geschlossene Schalen 8, 20, 40, 70 binden, die Mitten 14, 30, 55 nicht."""
    import kerne as KE
    assert abs(KE.schalenabstand(1.2) - 41.24) < 0.05
    d = KE.schalen_niveaus([8, 20, 40, 70, 14, 30, 55], 0.0, 0.0)
    assert all(d[:4] < 0) and all(d[4:] > 0)


def test_alpha_tunnel_geschlossen():
    """Tunnel im Coulomb-Wall geschlossen gerechnet: ∫_R^b √(1/r − 1/b) dr = √b·(arccos √x − √(x(1 − x))), x = R/b."""
    import numpy as np
    import kerne as KE
    lam = KE.halterlaenge(KE.a_C(KE.TROPFEN))
    hc, mu = K.hbar * K.c / KE.MEV * 1e15, K.m_u * K.c**2 / KE.MEV
    d = []
    for Z, A, Q, lgT in KE._alpha_strahler():
        m, R = mu * 4 * (A - 4) / A, lam * ((A - 4)**(1 / 3) + 4**(1 / 3))
        k = 2 * (Z - 2) * K.alpha * hc
        b = k / Q
        x = R / b
        I = math.sqrt(2 * m * k) * math.sqrt(b) * (math.acos(math.sqrt(x)) - math.sqrt(x * (1 - x))) / hc
        d.append(math.log10(math.log(2) * 2 * R / (K.c * 1e15 * math.sqrt(2 * Q / m))) + 2 * I / math.log(10) - lgT)
    assert abs(np.mean(d) - KE.alpha_tunnel()[0]) < 1e-3


def test_masse_aus_bahn_und_turm():
    """Erdbahn (1 AE, Sternjahr 365,256363 d) gibt die Sonnenmasse; im Turm ist 3Δs·h/R² = g·h/c² (von Hand)."""
    assert abs(K.masse_aus_bahn(K.AU, 365.256363 * 86400) / K.M_sun - 1) < 1e-4
    g = K.G * K.M_earth / K.R_ERDE**2
    assert abs(3 * K.faltung(K.M_earth) * 22.5 / K.R_ERDE**2 / (g * 22.5 / K.c**2) - 1) < 1e-12


def test_linie_des_wasserstoffs_von_hand():
    """Fermi-Form von Hand: ν = (16/3)·α²·cR∞·(m_e/m_p)·(g_p/2)·(1 + m_e/m_p)⁻³, cR∞ = α²m_ec²/(2h); dazu (1 + α/2π)(1 + 3α²/2)."""
    cR = K.alpha**2 * K.m_e * K.c**2 / (2 * K.h)
    r = K.m_e / K.m_p
    nu = 16 / 3 * K.alpha**2 * cR * r * K.g_p_CODATA / 2 * (1 + r)**-3 * (1 + K.alpha / (2 * math.pi)) * (1 + 1.5 * K.alpha**2)
    assert abs(K.hyperfein_wasserstoff() / nu - 1) < 1e-9


def test_zwei_umlaufsinne():
    """Reihe der Rückläufe von Hand: Σqⁿ und Σ(−1)ⁿ⁺¹qⁿ gegen Bose- und Fermi-Form; über die Stufen 7/8 = 1 − 2⁻³ (η(4)/ζ(4))."""
    for s in (0.3, 1.0, 2.5):
        q = math.exp(-s)
        assert abs(sum(q**n for n in range(1, 400)) - K.rueckkehr_summe(s, 1)) < 1e-12
        assert abs(sum((-1)**(n + 1) * q**n for n in range(1, 400)) - K.rueckkehr_summe(s, 2)) < 1e-12
    assert abs(K.strahlung_je_sinn(1) - math.pi**4 / 15) < 1e-8
    assert abs(K.strahlung_je_sinn(2) / K.strahlung_je_sinn(1) - 7 / 8) < 1e-9


def test_halbe_drehung_von_hand():
    """Stufen der Leiter mit halber Drehung (η(d+1)/ζ(d+1) = 1 − 2^{−d}), ∫n_F = ln 2, Kugelmittel 2/3, (11/4)^{1/3}."""
    for d in (1, 2, 3):
        assert abs(K.strahlung_je_sinn(2, d) / K.strahlung_je_sinn(1, d) - (1 - 2**-d)) < 1e-9
    assert abs(K.strahlung_je_sinn(2, 0) - math.log(2)) < 1e-9
    assert abs(K.kugelmittel_quer() - 2 / 3) < 1e-12
    assert abs(K.licht_gegen_neutrinos() - (11 / 4)**(1 / 3)) < 1e-9


def test_euler_sperre():
    """Sperrende Windungen gegen Rückläufe in ungeraden Moden (Euler), zwei getrennte Produkte."""
    for q in (0.1, 0.5, 0.9):
        a, b = K.sperre_und_ungerade(q)
        assert abs(a / b - 1) < 1e-12


def test_spiegel_halbe_drehung():
    """Produkte direkt gegen die Entwicklung im Tiefbereich; der Fehler der Windung folgt −e^{−2π²/s} (Spiegel mit halber
    Drehung), der der Faltung e^{−4π²/s}."""
    for s in (0.05, 0.2, 0.5):
        for z in (1, 2):
            a, b = K.alle_vielfachen(s, z)
            assert abs(a - b) < 1e-9
    a, b = K.alle_vielfachen(2.0, 2)
    assert abs((a - b) / -math.exp(-2 * math.pi**2 / 2.0) - 1) < 0.01
    for s in (0.5, 1.0, 3.0):                                  # Spiegelung der Sperre auf die ungeraden Moden
        q, qs = math.exp(-s), math.exp(-2 * math.pi**2 / s)
        links, rechts = 1.0, (q * qs)**(-1 / 24) / math.sqrt(2)
        for n in range(1, 3000):
            links *= 1 + q**n
            rechts *= 1 - qs**(2 * n - 1)
        assert abs(links / rechts - 1) < 1e-12


def test_atmen_mit_eigenem_omega_b():
    """H₀ aus a₀ mit vorgegebenem ω_b gleich dem Weg über das gemessene ω_b; von Hand h² = ω_b(1 + Q) + ω_ν + ω_r + ω_DE."""
    Q = K.q_je_radiant()
    assert abs(K.atmen_heute(1.08e-10, Q=Q, wb=K.omega_b)[0] / K.atmen_heute(1.08e-10, Q=Q)[0] - 1) < 1e-12
    wb, wnu = 0.0225, 0.06 / 93.14
    h2 = (wb * (1 + Q) + wnu + K.dichte_heute()[0] / K.rho_krit100
          + 3 * math.pi**3 * 1.08e-20 * 1.08 / (30 * K.G * K.c**2) / K.rho_krit100)
    assert abs(K.atmen_heute(1.08e-10, Q=Q, wb=wb, wnu=wnu)[0] / (100 * math.sqrt(h2)) - 1) < 1e-12


def test_maximum_dunkle_energie_von_hand():
    """Dichte der Dunklen Energie im Gang w₀ + w_a(1 − a): ρ ∝ a^{−3(1 + w₀ + w_a)}·e^{−3w_a(1 − a)}; ihr Maximum numerisch."""
    from scipy.optimize import minimize_scalar
    for w0, wa in ((-0.84, -0.56), (-0.8, -0.7)):
        lnrho = lambda z: -3 * (1 + w0 + wa) * math.log(1 / (1 + z)) - 3 * wa * (1 - 1 / (1 + z))
        z = minimize_scalar(lambda z: -lnrho(z), bounds=(0, 3), method="bounded", options=dict(xatol=1e-9)).x
        assert abs(z - K.maximum_dunkle_energie(w0, wa)) < 1e-6
        assert abs(math.log(K.gang_dunkle_energie(w0, wa, 0.7)) - lnrho(0.7)) < 1e-12


def test_streuung_im_feld_grenzen():
    """Plummer-Kugel nach der Regel: ganz im Tiefbereich σ⁴ = (4/81)·GMa₀ (jedes Profil), ganz Newtonsch σ² = πGM/(32R_e)
    (Virialsatz der Plummer-Kugel von Hand); ein äußeres Feld senkt σ."""
    M, R = 1e3 * K.M_sun, 20 * K.kpc
    assert abs(K.streuung_im_feld(M, R) / (4 / 81 * K.G * M * K.a0)**0.25 - 1) < 1e-3
    M, R = 1e13 * K.M_sun, 0.05 * K.kpc
    assert abs(K.streuung_im_feld(M, R) / math.sqrt(math.pi * K.G * M / (32 * R)) - 1) < 1e-4
    M, R = 1e7 * K.M_sun, 0.3 * K.kpc
    assert K.streuung_im_feld(M, R, 0.1 * K.a0) < K.streuung_im_feld(M, R)


def test_bahn_wendepunkte_kepler():
    """Weit über der Schwelle ist die Regel Newton: Bahn um 1 M☉ bei 1 AE mit 0,8·v_Kreis rein tangential ist am
    Apozentrum, das Perizentrum nach Kepler 2a − r mit a = 1/(2/r − v²/GM)."""
    M, r = K.M_sun, K.AU
    v = 0.8 * math.sqrt(K.G * M / r)
    peri, apo = K.bahn_wendepunkte(M, r, 0.0, v)
    a = 1 / (2 / r - v**2 / (K.G * M))
    assert abs(apo / r - 1) < 1e-9 and abs(peri / (2 * a - r) - 1) < 1e-6
