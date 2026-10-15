# -*- coding: utf-8 -*-
"""Kern des FTM-Modells: alle Konstanten und Formeln genau einmal.

Konstanten: CODATA 2018, GM☉ nach IAU 2015. Kosmische Messwerte (H₀, Ω) nach Planck 2018.
Jede Formel trägt den Baustein des Modells, in dem sie steht.
"""
import math

# ---------------------------------------------------------------- Konstanten
c = 299_792_458.0                 # m/s
G = 6.674_30e-11                  # m³/(kg s²)
h = 6.626_070_15e-34              # J s
hbar = h / (2 * math.pi)
k_B = 1.380_649e-23               # J/K
e = 1.602_176_634e-19             # C
eps0 = 8.854_187_8128e-12         # F/m
m_e = 9.109_383_7015e-31          # kg
m_p = 1.672_621_923_69e-27        # kg
alpha = e**2 / (4 * math.pi * eps0 * hbar * c)

GM_sun = 1.327_124_400_18e20      # m³/s² (IAU)
M_sun = GM_sun / G
GM_earth = 3.986_004_418e14       # m³/s²
M_earth = GM_earth / G

AU = 1.495_978_707e11             # m
R_ERDE, R_SONNE = 6.371e6, 6.957e8  # m, Radien von Erde und Sonne (Eingaben)
pc = 3.085_677_581_491_367e16     # m
kpc, Mpc, Gpc = 1e3 * pc, 1e6 * pc, 1e9 * pc
yr = 365.25 * 86400.0             # s (julianisches Jahr)
Gyr = 1e9 * yr

H0_kms = 67.4                     # km/s/Mpc (Planck 2018)
H0 = H0_kms * 1e3 / Mpc           # 1/s
Omega_b = 0.02237 / (H0_kms / 100)**2
Omega_c = 0.1200 / (H0_kms / 100)**2
Q_kosm = Omega_c / Omega_b        # aus der Messung eingesetzt (03)

a0 = 1.06e-10                     # m/s², Standardwert des Modells (03): SPARC und Atmen lesen dieselbe Schwelle
a0_literatur = 1.2e-10            # nur zum Vergleich


# ---------------------------------------------------------------- 01 Grundlagen
def faltung(M):
    """Δs = GM/(3c²): Faltungstiefe."""
    return G * M / (3 * c**2)


def energie(M):
    """E = 3c⁴·Δs/G: Energie als gebundene Faltung."""
    return 3 * c**4 * faltung(M) / G


def neigung(M, r):
    """g = 3c²Δs/r² = GM/r²: Neigung an der Kante des Trichters."""
    return 3 * c**2 * faltung(M) / r**2


def takt(M, r):
    """dτ/dt = √(1 − 6Δs/r): Abwicklungstakt am Ort r."""
    return math.sqrt(1 - 6 * faltung(M) / r)


def grenzradius(M):
    """r₀ = 6Δs: dort steht der Takt still."""
    return 6 * faltung(M)


def trichter(M, R, traeger=None, a0_=None):
    """Der Trichter der Masse M mit Hülle R in seinen drei Richtungen.
    parametrisch (Einbettung): Tiefe δ = 3Δs/R und Halter λ_H = R/N^{1/3}, N = M/Träger (Vorgabe Proton);
    transversal (innere Richtung): Schwelle r*, Lage R/r*, Rückkehranteil q = e^{−r*/R} an der Hülle;
    tangential (um den Knoten): Tempo der Hülle v/c = √δ."""
    a0_ = a0 if a0_ is None else a0_
    traeger = m_p if traeger is None else traeger
    d = 3 * faltung(M) / R
    rs = math.sqrt(G * M / a0_)
    return dict(tiefe=d, halter=R / (M / traeger)**(1 / 3), schwelle=rs, lage=R / rs, rueckkehr=math.exp(-rs / R),
                tempo=math.sqrt(d))


def masse_windung_faltung():
    """Windungen heben sich auf, Faltungen addieren sich: Die Faltung von N Trägern erreicht die Windung bei
    N^{2/3}·α_G = α (Bindung je Nukleon, unabhängig vom Halter), M = m_p·(α/α_G)^{3/2}; die Tiefe eines von Atomen
    gehaltenen Körpers (λ_H = a_B) ist dort δ = α²·m_e/m_p. Gibt (M, δ)."""
    return m_p * (alpha / (3 * faltung(m_p) / compton(m_p)))**1.5, alpha**2 * m_e / m_p


# ---------------------------------------------------------------- 03/07 Rückkehr und Schwelle
def schliesslaenge(a0_=a0):
    """Λ = c²/(2a₀): Schließlänge der Rückkehr, Länge der inneren Richtung."""
    return c**2 / (2 * a0_)


def rueckkehrzeit(a0_=a0):
    """τ_D = 2Λ/c = c/a₀: Kac-Rückkehrzeit, für jede Masse gleich."""
    return 2 * schliesslaenge(a0_) / c


def schwellenradius(M, a0_=a0):
    """r* = √(GM/a₀) = √(3Δs·ℓ₀), ℓ₀ = 2Λ."""
    return math.sqrt(G * M / a0_)


def s_(aN, a0_=a0):
    """s = √(a_N/a₀)."""
    return math.sqrt(aN / a0_)


def rueckkehranteil(aN, a0_=a0):
    """q = e^{−s}: Anteil, der je Umlauf zurückkehrt."""
    return math.exp(-s_(aN, a0_))


def verstaerkung(aN, a0_=a0):
    """ν = Σqⁿ = 1/(1 − q): Summe der Rückläufe."""
    if aN <= 0:
        return float("inf")
    return 1 / (-math.expm1(-s_(aN, a0_)))


def schalter(aN, a0_=a0):
    """χ = s/(e^s − 1)."""
    s = s_(aN, a0_)
    return s / math.expm1(s) if s > 0 else 1.0


def bahntempo_regel(M, r, a0_=a0):
    """Kreisbahn nach der Regel um die Punktmasse M: v = √(ν(a_N)·a_N·r)."""
    aN = G * M / r**2
    return math.sqrt(verstaerkung(aN, a0_) * aN * r)


def v_flach(M, a0_=a0):
    """v⁴ = GMa₀: flache Rotationskurve im tiefen Bereich."""
    return (G * M * a0_)**0.25


def theta_(a0_=a0, H=H0):
    """Θ = 2πa₀/(cH)."""
    return 2 * math.pi * a0_ / (c * H)


def vakuumteil(a0_=a0, H=H0):
    """Ω_s = (π²/45)·Θ²: Rückstellung des Ganzen, als Energie gelesen."""
    return math.pi**2 / 45 * theta_(a0_, H)**2


# ---------------------------------------------------------------- 02 Teilchen bis Stern
def compton(M):
    """λ̄ = ħ/(Mc)."""
    return hbar / (M * c)


def spiegelmasse(M):
    """Die Rückkehrgleichung mit D = GM/c ist die Schrödinger-Gleichung in imaginärer Zeit mit D = ħ/(2m̃): m̃ = ħc/(2GM) = m_Pl²/(2M)."""
    return hbar * c / (2 * G * M)


def falte_temperatur(M):
    """Temperatur der Falte: Ruheenergie der Spiegelmasse durch 4π, kT = m̃c²/(4π) = ħc/(24πΔs)."""
    return spiegelmasse(M) * c**2 / (4 * math.pi * k_B)


def kopplung_grav(M):
    """α_G(M) = 3Δs/λ̄ = GM²/(ħc)."""
    return 3 * faltung(M) / compton(M)


def sternmasse():
    """M_* = m_Pl³/m_p², m_Pl = √(ħc/G)."""
    m_pl = math.sqrt(hbar * c / G)
    return m_pl**3 / m_p**2


def bahntempo(M, r):
    """v = c·√(3Δs/r) auf der Kreisbahn."""
    return c * math.sqrt(3 * faltung(M) / r)


# ---------------------------------------------------------------- 08 Licht
def ablenkung(M, b):
    """α = 6Δs/b (Takt) + 6Δs/b (Faltung) = 12Δs/b."""
    return 12 * faltung(M) / b


def perihel_je_umlauf(M, a, ecc):
    """Δφ = 18πΔs/(a(1 − e²)) je Umlauf: (12π + 12π − 6π)·Δs/p."""
    return 18 * math.pi * faltung(M) / (a * (1 - ecc**2))


# ---------------------------------------------------------------- 09 Elektromagnetismus
def coulomb_konstante():
    return 1 / (4 * math.pi * eps0)


def bohr_radius():
    """a_B = ħ/(α m_e c)."""
    return hbar / (alpha * m_e * c)


def verlust_windung(r, v, M=m_e):
    """Verlust einer offenen Windung (Ladung e um eine Gegenladung e) je Umlauf gegen ihre Bindung: Zwei Umlaufsinne,
    der Schwerpunkt der Ladung schwingt, Abgabe als Dipol: ΔE/E = (8π/3)·(r_e/r)·(v/c), r_e = e²/(4πε₀Mc²)."""
    return 8 * math.pi / 3 * coulomb_konstante() * e**2 / (M * c**2) / r * v / c


def verlust_faltung(m, M, r):
    """Verlust einer Bahn der Masse m um M je Umlauf gegen ihre Bindung: ein Umlaufsinn, der Schwerpunkt der Faltung
    steht, Abgabe erst als Quadrupol: ΔE/E = (128π/5)·(m/M)·(v/c)⁵, v = √(GM/r)."""
    return 128 * math.pi / 5 * m / M * (math.sqrt(G * M / r) / c)**5


def fall_offene_windung():
    """Zeit, in der die Windung des Wasserstoffs fiele, wenn sie nicht schlösse: dr/dt = −(4/3)·r_e²c/r², t = a_B³/(4r_e²c)."""
    r_e = coulomb_konstante() * e**2 / (m_e * c**2)
    return bohr_radius()**3 / (4 * r_e**2 * c)


# ---------------------------------------------------------------- 09 Magnetismus: Windung mal Verwindung
g_e_CODATA = 2 * 1.001_159_652_181_28   # g-Faktor des Elektrons, CODATA 2018 (Messung; nicht in archiv/weltmodell/daten)
g_p_CODATA = 5.585_694_6893              # g-Faktor des Protons, CODATA 2018 (Inhalt, nicht Form)
NU_21CM = 1_420_405_751.768             # Hz, Hyperfeinlinie des Wasserstoffs (Hellwig u. a. 1970)
QED_BEIWERTE = (0.5, -0.328_478_965, 1.181_241_456, -1.9097, 9.16)  # g−2 des Elektrons in (α/π)ⁿ, n = 1 … 5 (Aoyama u. a. 2012)


def rueckkehr_windung_wechselnd(n):
    """Beiwert der n-ten Rückkehr der Windung in (α/π)ⁿ, wenn jede Rückkehr mit halber Drehung kommt (zwei Umlaufsinne,
    q → −q): (−1)^{n+1}·(α/2π)ⁿ/(α/π)ⁿ = (−1)^{n+1}/2ⁿ."""
    return (-1)**(n + 1) / 2**n


def hyperfein_wasserstoff(rueckkehr=True):
    """Linie des Wasserstoffs aus Windung mal Verwindung [Hz]: Moment der Windung g_e·μ_B (g_e = 2, mit der ersten Rückkehr
    2(1 + α/2π)) gegen das des Protons g_p·μ_N am Ort des Kerns, ΔE = (2μ₀/3)·g_e μ_B·g_p μ_N·|ψ(0)|², |ψ(0)|² = 1/(π a³) mit
    der Bahn um den gemeinsamen Schwerpunkt, a = a_B(1 + m_e/m_p), und dem Takt aus dem Tempo (Dirac) 1 + (3/2)α²."""
    ge = 2 * erste_rueckkehr_windung() if rueckkehr else 2.0
    a = bohr_radius() * (1 + m_e / m_p)
    return 2 * mu0 / 3 * ge * magneton(m_e) * g_p_CODATA * magneton(m_p) / (math.pi * a**3) * (1 + 1.5 * alpha**2) / h


def magneton(M):
    """μ = Q_el·c·a_J mit der Verwindung a_J = λ̄/2 des Spins ½: e·c·λ̄/2 = eħ/(2M); Elektron μ_B, Proton μ_N."""
    return e * c * compton(M) / 2


def spinzahl_teilchen(M):
    """j = a_J/(3Δs) eines Teilchens mit Spin ½: a_J = ħ/(2Mc) gegen 3Δs = GM/c², also m_Pl²/(2M²); j = 1 bei m_Pl/√2."""
    return hbar / (2 * M * c) / (3 * faltung(M))


def erste_rueckkehr_windung():
    """Die Windung kehrt je Radiant mit ihrem Gewicht α auf sich zurück: μ/μ_B = 1 + α/(2π) (erste Rückkehr)."""
    return 1 + alpha / (2 * math.pi)


def langevin(x):
    """L(x) = coth x − 1/x: Ausrichtung einer Faltung (alle Richtungen, ein Umlaufsinn) im Feld x = μB/kT."""
    return 1 / math.tanh(x) - 1 / x


def ausrichtung_windung(x):
    """Ausrichtung einer Windung mit zwei Umlaufsinnen (Spin ½) im Feld x = μB/kT: tanh x = 1 − 2f(2x), Fermi-Form."""
    return math.tanh(x)


def einrasten_magnet(t, h=0.0):
    """Eingerastete Windungen (Weiss): m = tanh((m + h)/t), t = T/T_C, h = μB/(kT_C); Kippstelle t = 1 (q = T_C/T = 1).
    Ohne Feld über der Kippstelle null, sonst der Ast mit dem Vorzeichen des Feldes."""
    from scipy.optimize import brentq
    if h == 0.0 and t >= 1.0:
        return 0.0
    g = lambda m: m - math.tanh((m + h) / t)
    lo, hi = (1e-12, 1.0) if h >= 0 else (-1.0, -1e-12)
    return brentq(g, lo, hi, xtol=1e-15)


def hysteresefenster(t):
    """Kippfeld der eingerasteten Windungen, t < 1: h_c = t·artanh(√(1 − t)) − √(1 − t) [kT_C/μ], ≈ (2/3)(1 − t)^{3/2}."""
    w = math.sqrt(1 - t)
    return abs(t * math.atanh(w) - w)


def rueckkehr_windung_reihe(n=2):
    """Wiederholte Rückkehr der Windung als geometrische Reihe: Σ_{k≤n} (α/2π)^k; zum Vergleich mit der Messung."""
    return sum((alpha / (2 * math.pi))**k for k in range(n + 1))


mu0 = 4 * math.pi * 1e-7 * 1.000_000_000_55   # Vs/(Am), CODATA 2018
J_ERDE = 5.86e33                              # kg m²/s, Drehimpuls der Erde (Eingabe, Trägheitsmoment 0,3307·M·R²)
R_GPB = 7.027e6                               # m, Bahnradius von Gravity Probe B (Eingabe)


def verwindung(J, M):
    """a_J = J/(Mc): die Verwindung eines Trichters als Länge."""
    return J / (M * c)


def mitfuehrung(J, r):
    """Mitführung der Bahnebene durch die Verwindung (Lense-Thirring), polare Bahn im Mittel: Ω = GJ/(2c²r³) = 3Δs·c·a_J/(2r³)."""
    return G * J / (2 * c**2 * r**3)



# ---------------------------------------------------------------- 08 An den Messungen (Messwerte mit Quelle; nicht in archiv/weltmodell/daten)
GPS_BAHN = 26_561.75e3                # m, große Halbachse der GPS-Bahnen (IS-GPS-200)
GPS_WERK = 4.4647e-10                 # Werkseinstellung: Die Bahnuhren gehen um diesen Anteil langsamer (10,22999999543 statt 10,23 MHz; IS-GPS-200)
CASSINI_GAMMA = (2.1e-5, 2.3e-5)      # γ − 1 aus der Laufzeit zur Sonde Cassini (Bertotti, Iess, Tortora 2003)
MERKUR_GEMESSEN = (43.11, 0.45)       # ″ im Jahrhundert, Drehung des Perihels über die Störungen hinaus (Clemence 1947)
GPB_KREISEL = (6601.8, 18.3)          # mas/a, Präzession der Kreisel in der Bahnebene (Gravity Probe B, Everitt u. a. 2011)
GPB_MITFUEHRUNG = (37.2, 7.2)         # mas/a, Mitführung (Gravity Probe B, Everitt u. a. 2011)
PSR_1913 = dict(m1=1.438, m2=1.390, Pb=0.322997448918 * 86400, e=0.6171340, verhaeltnis=(0.9983, 0.0016))  # Weisberg, Huang 2016
M87 = dict(M=6.5e9, D=16.8, ring=(42.0, 3.0))   # M☉, Mpc, µas (Event Horizon Telescope 2019)
GW170817_DC = (-3e-15, 7e-16)         # Grenzen von (v_GW − c)/c (Abbott u. a. 2017)
ALPHA_G = (0.75, 0.13, 0.16)          # Fall von Antiwasserstoff in g, ± (Statistik und System) ± (Simulation) (ALPHA-g, Nature 621, 716, 2023)
S2 = dict(a=0.12540, e=0.88466, P=16.0518, R0=8249.0, f=(1.10, 0.19))  # ″, –, a, pc; f_SP der Drehung je Umlauf (GRAVITY 2020)
SGRA_RING = (51.8, 2.3)               # µas, Ring um Sgr A* (Event Horizon Telescope 2022, Beobachtung 2017)
J0737 = dict(m1=1.3381, m2=1.2489, Pb=0.10225156 * 86400, e=0.087779, genau=1.3e-4)  # M☉, s; Ṗ_b bestätigt auf 1,3·10⁻⁴ (Kramer u. a. 2021)
POUND = dict(h=22.5, genau=0.01)      # m Turmhöhe, Genauigkeit der Rotverschiebung (Pound, Rebka, Snider 1960–1965)


def uhr_bahn(M, R, r):
    """Gang einer Uhr auf der Kreisbahn r gegen eine am Rand R (A3, erste Ordnung): Am Rand fehlt dem Takt 3Δs/R, auf der Bahn
    3Δs/r aus der Tiefe und ½(v/c)² = (3/2)Δs/r aus dem Tempo; die Bahnuhr geht um 3Δs/R − (9/2)Δs/r vor."""
    d = faltung(M)
    return 3 * d / R - 4.5 * d / r


def kreisel_bahn(M, r):
    """Präzession eines Kreisels auf der Kreisbahn r: (3/2)·(3Δs/r)·Ω, Ω = √(GM/r³) [rad/s]."""
    return 1.5 * 3 * faltung(M) / r * math.sqrt(G * M / r**3)


def doppelpulsar(m1, m2, Pb, e):
    """Abnahme der Umlaufzeit eines Paares durch den Faltungsquadrupol (transversal, ein Umlaufsinn):
    Ṗ_b = −(192π/5)·(G/c³)^{5/3}·(2π/P_b)^{5/3}·m₁m₂(m₁ + m₂)^{−1/3}·f(e), f(e) = (1 + 73e²/24 + 37e⁴/96)/(1 − e²)^{7/2}."""
    f = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e**2)**3.5
    return -192 * math.pi / 5 * (G / c**3)**(5 / 3) * (2 * math.pi / Pb)**(5 / 3) * m1 * m2 * (m1 + m2)**(-1 / 3) * f


def masse_aus_bahn(a, P):
    """Masse aus Halbachse und Umlaufzeit: T² = 4π²a³/(3c²Δs) gibt Δs und M = 3c²Δs/G."""
    return 3 * c**2 * (4 * math.pi**2 * a**3 / (3 * c**2 * P**2)) / G


def schatten(M, D):
    """Durchmesser des Schattens einer Falte: Die Photonensphäre 9Δs fängt jedes Licht mit Stoßparameter b ≤ 9√3·Δs;
    am Himmel 18√3·Δs/D [rad]."""
    return 18 * math.sqrt(3) * faltung(M) / D

def feld_treffpunkt(M):
    """Feld, bei dem die Feldrate eB/M den Compton-Takt c/λ̄ erreicht: B = M²c²/(eħ) (Schwinger-Feld), E = B·c."""
    return M**2 * c**2 / (e * hbar)


def diamagnetismus_h():
    """Molare Suszeptibilität des Wasserstoffatoms (Langevin, ⟨r²⟩ = 3a_B²): χ = −μ₀·N_A·e²·a_B²/(2m_e) [m³/mol]."""
    N_A = 6.022_140_76e23
    return -mu0 * N_A * e**2 * bohr_radius()**2 / (2 * m_e)


def flussquant():
    """Φ₀ = h/(2e): Schließung der Phase eines eingerasteten Paares je Umlauf."""
    return h / (2 * e)


def josephson_einrasten(I, I_c, R):
    """Zwei eingerastete Windungen mit Spannung: ψ̇ = Δω − K·sin ψ mit Δω = 2eRI/ħ, K = 2eRI_c/ħ (dieselbe Form wie das
    Einrasten zweier Zyklen); eingerastet (V = 0) für |I| ≤ I_c. Gibt (Δω, K, eingerastet)."""
    dw, K_ = 2 * e * R * I / hbar, 2 * e * R * I_c / hbar
    return dw, K_, abs(dw) <= K_


def london_tiefe(n_e):
    """Schließlänge des Feldes im eingerasteten Zustand: λ_L = √(m_e/(μ₀·n·e²))."""
    return math.sqrt(m_e / (mu0 * n_e * e**2))


# ---------------------------------------------------------------- 03 Umgebung, Naht, Zelle
def verstaerkung_aussenfeld(aN, g_e, lage, a0_=a0):
    """ν im äußeren Feld: Betrag des Gesamtfelds, senkrecht √(a_N² + g_e²), gleichgerichtet a_N + g_e."""
    betrag = math.hypot(aN, g_e) if lage == "senkrecht" else aN + g_e
    return verstaerkung(betrag, a0_)


def naht_radius(M, Q=Q_kosm, a0_=a0):
    """R_Q = r*/(−ln q_kosm): dort erreicht das örtliche q das kosmische."""
    return schwellenradius(M, a0_) / (-math.log(Q / (1 + Q)))


def form_zelle(Q=Q_kosm):
    """r*(M_b der Zelle)/Λ = √(6/ν_kosm), ν_kosm = 1 + Q; der Schluss verlangt ≤ 1, also Q ≥ 5."""
    return math.sqrt(6 / (1 + Q))


# ---------------------------------------------------------------- 02 Grenzmasse
m_u = 1.660_539_066_60e-27        # kg
m_n = 1.008_664_915_95 * m_u       # kg, Neutron (AME2020)
M_MYON = 105.658_3755 * 1e6 * e / c**2   # kg, Myon (CODATA 2018)
TAU_MYON = 2.196_9811e-6                 # s, mittlere Lebensdauer des Myons (PDG)
MCF_FUSIONEN = 150                       # Verschmelzungen je Myon in dichtem D-T (Los Alamos, Jones u. a. 1986)
MCF_HAFTEN = 0.006                       # wirksame Haftwahrscheinlichkeit des Myons am α (Übersicht 2026)
NU_U233 = 2.487                          # Neutronen je Spaltung von ²³³U, thermisch (IAEA Standards 2017)
U233_SPALT, U233_FANG = 531.16, 45.25      # b, Spaltung und Einfang von ²³³U bei 2200 m/s (Tabelle nach ENDF, nuclear-power.com)
QUERSCHNITT = dict(Li6=940.0, C_streu=5.551, C_abs=0.0035, Th232=7.37)   # b, 2200 m/s, gebunden (Sears 1992, periodictable)
SPALLATION = (0.59, 9.6)                   # GeV, Neutronen je Proton in dickem Pb/Bi (PSI, Nucl. Instr. Meth. B 2003)
WIRKUNGSGRAD_THERMOAKUSTIK = 0.30          # Wärme zu Schall, ohne bewegte Teile (Backhaus und Swift, Nature 1999)
m_H = 1.007_825_032_23 * m_u       # kg, Wasserstoffatom (AME2020)
_LANE_EMDEN = {}


def lane_emden(n=3.0, h=1e-3):
    """Form des Trichters mit Druck ∝ Dichte^{1+1/n}: θ'' + (2/ξ)θ' + θⁿ = 0, θ(0) = 1, θ'(0) = 0.
    Gibt den Rand ξ₁ (θ = 0) und den Beiwert ω_n = −ξ₁²·θ'(ξ₁); n = 3: ξ₁ = 6,897, ω₃ = 2,018."""
    if (n, h) in _LANE_EMDEN:
        return _LANE_EMDEN[(n, h)]

    def f(xi, y):
        th, dth = y
        return (dth, -2 * dth / xi - max(th, 0.0)**n)
    xi = h
    y = (1 - xi**2 / 6 + n * xi**4 / 120, -xi / 3 + n * xi**3 / 30)
    while y[0] > 0:
        k1 = f(xi, y)
        k2 = f(xi + h / 2, (y[0] + h / 2 * k1[0], y[1] + h / 2 * k1[1]))
        k3 = f(xi + h / 2, (y[0] + h / 2 * k2[0], y[1] + h / 2 * k2[1]))
        k4 = f(xi + h, (y[0] + h * k3[0], y[1] + h * k3[1]))
        alt = (xi, y)
        y = (y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
             y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))
        xi += h
    t = alt[1][0] / (alt[1][0] - y[0])
    xi1 = alt[0] + t * h
    dth = alt[1][1] + t * (y[1] - alt[1][1])
    _LANE_EMDEN[(n, h)] = (xi1, -xi1**2 * dth)
    return _LANE_EMDEN[(n, h)]


def grenzmasse(mu_e=2.0):
    """M = (ω₃·√(3π)/2)·(ħc/G)^{3/2}/(μ_e·m_u)², ω₃ = 2,018 der Beiwert der Form n = 3 (lane_emden)."""
    return lane_emden(3.0)[1] * math.sqrt(3 * math.pi) / 2 * (hbar * c / G)**1.5 / (mu_e * m_u)**2


def druck_elektronen(xf):
    """Zählregel der Elektronen: Druck durch n_e·m_ec², Π(x) = (3/8x³)·[x(2x²/3 − 1)√(1 + x²) + arsinh x], x = λ̄_e·k_F
    (x²/5 langsam, x/4 schnell)."""
    return 3 / (8 * xf**3) * (xf * (2 * xf**2 / 3 - 1) * math.sqrt(1 + xf**2) + math.asinh(xf))


def decke_entartet(mu_je_mu_e):
    """Höchste Temperatur eines Trichters, dessen Tiefe die Zählregel der Elektronen hält (an der Grenzmasse):
    kT_max = (μ/μ_e)·m_ec²·max_x[x/4 − Π(x)]; gibt (kT_max [J], x am Maximum)."""
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda xf: -(xf / 4 - druck_elektronen(xf)), bounds=(0.01, 10), method="bounded", options=dict(xatol=1e-10))
    return mu_je_mu_e * m_e * c**2 * -r.fun, r.x


def polytrope(n):
    """Lane-Emden θ'' + 2θ'/ξ = −θⁿ: Rand ξ₁, −θ'(ξ₁) und der Faktor der Kerntemperatur f = 1/((n + 1)·ξ₁·|θ'(ξ₁)|)
    in kT_c = f·δ·m̄c², δ = GM/(Rc²)."""
    from scipy.integrate import solve_ivp
    ende = lambda x, y: y[0]
    ende.terminal = True
    s = solve_ivp(lambda x, y: [y[1], -max(y[0], 0)**n - 2 * y[1] / x], [1e-6, 50], [1 - 1e-12 / 6, -1e-6 / 3],
                  events=ende, rtol=1e-12, atol=1e-14)
    x1, dt1 = s.t_events[0][0], s.y_events[0][0][1]
    return x1, -dt1, 1 / ((n + 1) * x1 * -dt1)


# ---------------------------------------------------------------- 06 Labor
def rc_kette_verstaerkung(n_R, RG):
    """RC-Kette mit Rückführung des Abgriffs n_R auf den Eingang: ν = 1/(1 − e^{−s}), s = n_R·arcosh(1 + RG/2)."""
    return 1 / (-math.expm1(-n_R * math.acosh(1 + RG / 2)))


# ---------------------------------------------------------------- 03 Systeme im Tiefbereich
def glaettungslaenge(M, gezeit):
    """Glättungslänge ℓ [m] eines kompakten Trichters der Masse M [kg] in der Gezeit 𝒯 [1/s²] des Elternfelds: ℓ³ = GM/𝒯;
    bis ℓ bleibt er Newtonsch (Einzeltrichter)."""
    return (G * M / gezeit)**(1 / 3)


def gezeit_flach(v, R):
    """Gezeit 𝒯 [1/s²] einer flachen Rotationskurve v [m/s] im Abstand R [m]: −dg/dR = v²/R² = Ω²."""
    return (v / R)**2


BREMSSTRAHLUNG = 1.4e-40 * 1.2                  # W·m³·K^(−1/2), Λ = Wert·√T für vollständig ionisiertes Gas, Gaunt-Faktor 1,2 (Rybicki und Lightman 1979)


WEITE_DOPPELSTERNE = dict(abstand_pc=130, kAU=(1, 30), mond_zuwachs=0.20)   # Gaia DR3: kein Zuwachs, Newton (Cookson u. a. 2026)


def streuung_kugel(M, a0_=a0):
    """Druckgestützte Kugel ganz im Tiefbereich, isotrop, jedes Profil: σ_v⁴ = (4/81)·GMa₀."""
    return (4 / 81 * G * M * a0_)**0.25


SONNE_GALAKTISCH = dict(R0=8.249,       # kpc, Abstand zum Zentrum (GRAVITY, wie bei S2)
                        z=0.0208,       # kpc über der Scheibe (Bennett und Bovy 2019)
                        v=(12.9, 245.6, 7.78))  # km/s gegen das Zentrum (Drimmel und Poggio 2018, nach Reid und Brunthaler 2004)


def bahn_wendepunkte(M, r, v_r, v_t, a0_=a0):
    """Peri- und Apozentrum [m] einer Bahn im Feld der Punktmasse M nach der Regel, g = ν(GM/r²)·GM/r², aus Ort r [m],
    radialer und tangentialer Geschwindigkeit [m/s]: Nullstellen von ½v_r² + L²/2r² − L²/2x² − (Φ(x) − Φ(r))."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    g = lambda x: verstaerkung(G * M / x**2, a0_) * G * M / x**2
    L = r * v_t
    f = lambda x: 0.5 * v_r**2 + 0.5 * v_t**2 - L**2 / (2 * x**2) - (quad(g, r, x, limit=200)[0] if x > r else
                                                                     -quad(g, x, r, limit=200)[0])
    innen, aussen = r * (1 - 1e-9), r * (1 + 1e-9)            # steht der Körper am Wendepunkt, ist r selbst einer
    if f(innen) <= 0:
        peri = r
    else:
        peri = brentq(f, r * 1e-4, innen) if f(r * 1e-4) < 0 else r * 1e-4
    if f(aussen) <= 0:
        return peri, r
    hi = r * 1.5
    while f(hi) > 0 and hi < r * 1e4:
        hi *= 1.5
    return peri, (brentq(f, aussen, hi) if f(hi) < 0 else float("inf"))


def regel_monopol(g_i, g_e=0.0, a0_=a0):
    """Gemittelte Radialkomponente der Regel auf das Gesamtfeld, ½∫ν(|g_i·e_r + g_e·e_z|)·(g_i + g_e·μ)dμ: was ein
    gestapeltes Linsensignal um einen Bereich im gleichförmigen äußeren Feld g_e liest; g_e = 0 gibt ν(g_i)·g_i."""
    import numpy as np
    mu, w = np.polynomial.legendre.leggauss(64)
    g = np.sqrt(g_i**2 + g_e**2 + 2 * g_i * g_e * mu)
    return float(np.sum(w / 2 / -np.expm1(-np.sqrt(g / a0_)) * (g_i + g_e * mu)))


def linse_umgerechnet(r, g, R):
    """Was die Messung um einen kugeligen Bereich als Beschleunigung liest: aus dem Feld g(r) [m/s²] auf dem Gitter r [m]
    die eingeschlossene Masse r²g/G, ihre Dichte, projiziert Σ(R) und M_2D(<R), ΔΣ = M_2D/(πR²) − Σ, umgerechnet wie die
    Messung g = 4G·ΔΣ. Exakt für g ∝ 1/r (dort gleich g), für eine Punktmasse (4/π)·GM/R²."""
    import numpy as np
    M = r**2 * g / G
    dM, lr = np.gradient(M, r), np.log(r)                        # dM/dr = 4πr²ρ
    aus = []
    for Ri in np.atleast_1d(R):
        s = np.concatenate([[0.0], np.geomspace(1e-5 * Ri, math.sqrt(r[-1]**2 - Ri**2), 4000)])
        rr = np.sqrt(Ri**2 + s**2)
        sigma = np.trapezoid(2 * np.interp(np.log(rr), lr, dM) / (4 * math.pi * rr**2), s)
        rs = np.geomspace(Ri, r[-1], 4000)
        m2 = np.interp(math.log(Ri), lr, M) + np.trapezoid(np.interp(np.log(rs), lr, dM) * (1 - np.sqrt(1 - (Ri / rs)**2)), rs)
        aus.append(4 * G * (m2 / (math.pi * Ri**2) - sigma))
    return np.array(aus)


M_B_MILCHSTRASSE = 6.1e10            # M☉, Baryonen bis 50 kpc (Flynn u. a. 2006, nach Bland-Hawthorn und Gerhard 2016)
M_B_M31 = 1.4e11                      # M☉, Sterne (Tamm u. a. 2012) und kaltes Gas 7,7·10⁹ (Chemin u. a. 2009)
M_V_SONNE = 4.83                      # absolute V-Helligkeit der Sonne


def streuung_im_feld(M, R_e, g_e=0.0, a0_=a0, n=200, M_B=None):
    """Sichtlinienstreuung σ [m/s] einer druckgestützten Plummer-Kugel (Masse M, projizierter Halblichtradius R_e) im
    äußeren Newtonfeld g_e der Umgebung, isotrop. Die Regel wirkt auf das Gesamtfeld, g = ν(|g_N,i + g_e|)·(g_N,i + g_e);
    die Bewegung als Ganzes fällt im Virialsatz heraus (∫ρ·r = 0): σ² = −∫ρ r·g dV/(3M), integriert über den Massenanteil x
    und den Winkel μ zum äußeren Feld. Ohne äußeres Feld ganz im Tiefbereich σ⁴ = (4/81)·GMa₀, ganz Newtonsch
    σ² = πGM/(32R_e). M_B: die Baryonen, die der Bereich einmal gebunden hat (Gedächtnis der Bindung, gleiches Profil); die
    Neigung am Ort liest M, die Rückkehr (innere Richtung) M_B: g = g_N + (ν(|g_N,B|) − 1)·g_N,B."""
    import numpy as np
    x, wx = np.polynomial.legendre.leggauss(n)
    mu, wmu = np.polynomial.legendre.leggauss(64)
    x, wx = (x + 1) / 2, wx / 2                                  # Massenanteil 0–1
    r = R_e * x[:, None]**(1 / 3) / np.sqrt(1 - x[:, None]**(2 / 3))
    g_i = G * M * x[:, None] / r**2                              # Newtonsch, nach innen
    g_b = g_i if M_B is None else G * M_B * x[:, None] / r**2   # Faltung, die die Rückkehr liest
    gn = np.sqrt(g_b**2 + g_e**2 - 2 * g_b * g_e * mu[None, :])
    nu = 1 / -np.expm1(-np.sqrt(gn / a0_))
    g = (g_e * mu[None, :] - g_i) + (nu - 1) * (g_e * mu[None, :] - g_b)
    W = np.sum(wx[:, None] * wmu[None, :] / 2 * r * g)
    return math.sqrt(-W / 3)


def flaechendichte_schwelle(a0_=a0):
    """Ein System liegt als Ganzes unter der Schwelle, wenn Σ < a₀/(πG)."""
    return a0_ / (math.pi * G)


def filament_schwelle(mu, a0_=a0):
    """Schwelle eines Filaments der Linienmasse μ: r*_fil = 2Gμ/a₀."""
    return 2 * G * mu / a0_


# ---------------------------------------------------------------- 07 Umlaufphase, Maximum der Zusatzbeschleunigung
def umlaufphase(aN, a0_=a0):
    """Re τ_M = ω_rad/ω_azim, Punktmasse: (Re τ_M)² = 1 + χ."""
    return math.sqrt(1 + schalter(aN, a0_))


def maximum_zusatz():
    """Stufe d = 2 der Bose-Leiter, s²·n_B(s), hat ihr Maximum bei 2(1 − e^{−s}) = s (Wien in d = 2).
    Gibt (a_N/a₀, (a − a_N)/a₀) an dieser Stelle."""
    lo, hi = 1.0, 2.0
    for _ in range(200):
        m = (lo + hi) / 2
        if 2 * (-math.expm1(-m)) - m > 0:
            lo = m
        else:
            hi = m
    s = (lo + hi) / 2
    return s**2, s**2 / math.expm1(s)


# ---------------------------------------------------------------- 08 Leiter der Längen, Licht am kollektiven Trichter
l_P = math.sqrt(hbar * G / c**3)
m_Pl = math.sqrt(hbar * c / G)


def leiterzahl(a0_=a0):
    """ε_L = ℓ_P/ℓ₀, ℓ₀ = c²/a₀ = 2Λ."""
    return l_P / (c**2 / a0_)


def schnittmassen(a0_=a0):
    """M_min = ε_L·m_Pl, M_c = ε_L^{1/3}·m_Pl, m_Pl, M_U = m_Pl/ε_L."""
    e = leiterzahl(a0_)
    return e * m_Pl, e**(1 / 3) * m_Pl, m_Pl, m_Pl / e


def ablenkung_grenz(M, a0_=a0):
    """Ablenkung am kollektiven Trichter für b ≫ r*: α_∞ = 6πΔs/r* = 2π·(v_flach/c)²."""
    return 2 * math.pi * (v_flach(M, a0_) / c)**2


# ---------------------------------------------------------------- 04 Entropie der Engstelle
def entropie_falte(M):
    """S/k = 36π·(Δs/ℓ_P)²."""
    return 36 * math.pi * (faltung(M) / l_P)**2


# ---------------------------------------------------------------- 03 Kosmos: Atmen, Zelle, Fluchtmaß des Ganzen
T_CMB = 2.7255                    # K (Fixsen 2009)
N_eff = 3.046                     # Neutrinos (Standardwert)
omega_b, omega_c = 0.02237, 0.1200  # Ω·h² (Planck 2018)
omega_b_fehler, omega_c_fehler = 0.00015, 0.0012  # 1σ (Planck 2018)
theta_stern_rel = 0.00031 / 1.04110  # relativer 1σ des Schallwinkels 100θ_* = 1,04110 ± 0,00031 (Planck 2018)
z_stern = 1089.92                   # letzte Streufläche (Planck 2018)
H100 = 100e3 / Mpc
rho_krit100 = 3 * H100**2 / (8 * math.pi * G)


def dichte_photonen():
    """Dichte der Photonen heute aus T_CMB."""
    return math.pi**2 / 15 * (k_B * T_CMB)**4 / (hbar * c)**3 / c**2


def dichte_heute(Q=None):
    """Dichten heute als Masse je Volumen: (ρ_r, ρ_b, ρ_m); Strahlung aus T_CMB, Neutrinos mit N_eff.
    Q = None: ω_c gemessen; sonst ω_c = Q·ω_b (Mode ξ mit dem Gewicht Q)."""
    rho_r = dichte_photonen() * (1 + 7 / 8 * (4 / 11)**(4 / 3) * N_eff)
    wc = omega_c if Q is None else Q * omega_b
    return rho_r, omega_b * rho_krit100, (omega_b + wc) * rho_krit100


def schallhorizont(Q=None):
    """r_d [Mpc] am Ende der Schallwelle des frühen Plasmas; die Mode ξ wirkt dort linear wie Staub, also gilt
    die Standardrechnung, hier als Anpassung an sie (Aubourg u. a. 2015, 0,02 %), ω_ν = 0,0006.
    Q = None: ω_c gemessen; sonst ω_c = Q·ω_b."""
    wc = omega_c if Q is None else Q * omega_b
    return 55.154 * math.exp(-72.3 * 0.0006**2) / ((omega_b + wc)**0.25351 * omega_b**0.12807)


def dichte_vakuumteil(a0_=a0):
    """ρ_s = π³a₀²/(30Gc²), zeitfest (Λ fest)."""
    return math.pi**3 * a0_**2 / (30 * G * c**2)


def zelle_laenge(a=1.0, a0_=a0, rho_m=None):
    """Räumliche Zelle L aus a₀ = 3c⁴/(2Gρ_m L³), mitbewegt: L ∝ a; rho_m: Dichte von Baryonen und Mode heute (sonst
    gemessen)."""
    return (3 * c**4 / (2 * G * (dichte_heute()[2] if rho_m is None else rho_m) * a0_))**(1 / 3) * a


def schwelle_zelle(a=1.0, a0_=a0):
    """s der ganzen Zelle am Raum: r*(Baryonen der Zelle)/R, R Kugel gleichen Volumens; wächst die Hülle, fällt s ∝ 1/a."""
    L = zelle_laenge(1.0, a0_)
    r_stern = math.sqrt(G * dichte_heute()[1] * L**3 / a0_)
    return r_stern / (L * a * (3 / (4 * math.pi))**(1 / 3))


def fluchtmass_ganzes(s):
    """Rückkehrzahl des Ganzen nach Kac, ν_G = μ_R/μ_N: Maß aller Tori, die ganz in der homogenen Kugel bleiben,
    unter der Regel gegen Newton am selben Rand; s = r*/R. Je Ort ist das Geschwindigkeitsmaß ein Ellipsoid,
    v_r² + v_t²(1 − r²/R²) ≤ 2(Φ(R) − Φ(r)), Volumen (4π/3)·(2ΔΦ)^{3/2}/(1 − r²/R²)."""
    from scipy.integrate import quad
    a0_ = 1 / s**2                                       # Einheiten R = 1, GM = 1
    beschl = lambda r: verstaerkung(r, a0_) * r if r > 0 else 0.0
    dphi = lambda r: quad(beschl, r, 1, limit=200)[0]
    mu_R = quad(lambda r: r * r * (2 * dphi(r))**1.5 / (1 - r * r), 1e-7, 1 - 1e-10, limit=400)[0]
    mu_N = quad(lambda r: r * r * math.sqrt(1 - r * r), 0, 1)[0]   # Newton: 2ΔΦ = 1 − r²
    return mu_R / mu_N


def abstand_mitbewegt(dichte_de, z=None, Q=None):
    """Mitbewegter Abstand D_M bis z (Vorgabe: letzte Streufläche), flach; dichte_de(a) die Dichte der Dunklen Energie.
    Das ist, was der Beobachter sieht: Jede Schale geht mit ihrer eigenen, früheren Epoche ein."""
    from scipy.integrate import quad
    z = z_stern if z is None else z
    rho_r, _, rho_m = dichte_heute(Q)
    H = lambda a: math.sqrt(8 * math.pi * G / 3 * (rho_r / a**4 + rho_m / a**3 + dichte_de(a)))
    return quad(lambda x: c / (math.exp(x) * H(math.exp(x))), -math.log(1 + z), 0.0, limit=500, epsabs=0, epsrel=1e-10)[0]


_TAFEL = {}


def fluchtmass_ganzes_glatt(s):
    """ν_G(s) aus einer Tafel (120 Stützstellen, kubisch in ln s), für Läufe über die Zeit."""
    import numpy as np
    from scipy.interpolate import CubicSpline
    if "sp" not in _TAFEL:
        sg = np.exp(np.linspace(math.log(0.2), math.log(60.0), 120))
        _TAFEL["sp"] = CubicSpline(np.log(sg), np.log([fluchtmass_ganzes(x) for x in sg]))
    return math.exp(float(_TAFEL["sp"](math.log(s))))


def rueckkehr_summe(s, sinne=1):
    """Summe der Rückläufe Σ_{n≥1} (±1)^{n+1}·e^{−ns}: ein Umlaufsinn gibt die Bose-Form 1/(eˢ − 1), zwei Umlaufsinne (jede
    Rückkehr mit halber Drehung, q → −q) die Fermi-Form 1/(eˢ + 1)."""
    return 1 / (math.exp(s) - (1 if sinne == 1 else -1))


def strahlung_je_sinn(sinne=1, d=3):
    """Stufe d der Leiter über alle s: ∫s^d·(Summe der Rückläufe) ds; bei einem Umlaufsinn Γ(d + 1)ζ(d + 1) (d = 3: π⁴/15),
    bei zweien (1 − 2^{−d})-mal so viel, d = 0: ln 2 (bei einem Sinn unbeschränkt)."""
    from scipy.integrate import quad
    return quad(lambda x: x**d * rueckkehr_summe(x, sinne), 1e-12, 80, limit=400)[0]


def sperre_und_ungerade(q, n=4000):
    """Windungen, die einander sperren (jede Mode höchstens einmal), ∏(1 + qⁿ), und Rückläufe nur in ungeraden Moden,
    ∏1/(1 − q^{2n−1}); nach Euler gleich."""
    a = b = 1.0
    for k in range(1, n + 1):
        a *= 1 + q**k
        b /= 1 - q**(2 * k - 1)
    return a, b


def alle_vielfachen(s, sinne=1, n=200000):
    """Alle Vielfachen einer Schleife, ln Z: Faltung (ein Sinn) −Σ ln(1 − e^{−ns}), Windung mit Sperre (zwei Sinne)
    Σ ln(1 + e^{−ns}); dazu die Entwicklung im Tiefbereich, π²/(6s) + ½·ln(s/2π) − s/24 bzw. π²/(12s) − ½·ln 2 + s/24.
    Gibt (ln Z, Entwicklung)."""
    if sinne == 1:
        z = -sum(math.log1p(-math.exp(-k * s)) for k in range(1, n))
        return z, math.pi**2 / (6 * s) + 0.5 * math.log(s / (2 * math.pi)) - s / 24
    z = sum(math.log1p(math.exp(-k * s)) for k in range(1, n))
    return z, math.pi**2 / (12 * s) - 0.5 * math.log(2) + s / 24


def kugelmittel_quer():
    """Anteil der zwei Richtungen quer zu einer Achse im Mittel über die Kugel, ∫sin²θ dΩ/(4π) = 2/3."""
    from scipy.integrate import quad
    return quad(lambda th: math.sin(th)**2 * math.sin(th) * 2 * math.pi, 0, math.pi)[0] / (4 * math.pi)


def licht_gegen_neutrinos():
    """Licht gegen Neutrinos nach dem Zerstrahlen von Elektron und Positron: Diese zählen mit zwei Sinnen und zwei Ladungen
    4·(7/8) gegen die 2 des Lichts, T_γ/T_ν = (1 + 2·(7/8))^{1/3} = (11/4)^{1/3}."""
    return (1 + 2 * strahlung_je_sinn(2) / strahlung_je_sinn(1))**(1 / 3)


def zustandsteil_je_gas(Q=None):
    """Fluchtmaß der Bindung je gebundenem Baryon: Gezählt werden die Tori innerhalb der Naht, R_Q = r*/(−ln q) mit
    q = Q/(1 + Q); für jeden Bereich ist also s = −ln q, unabhängig von seiner Masse, und das inaktive Maß ν_G(−ln q) − 1.
    Mit ρ = ν·ρ_s ist der Zustandsteil (ν_G(−ln q) − 1)·F·ρ_s, F der als Gas gebundene Anteil der Baryonen."""
    Q = Q_kosm if Q is None else Q
    return fluchtmass_ganzes(-math.log(Q / (1 + Q))) - 1


# Kaltes Gas im heutigen Kosmos, nur Wasserstoff: atomar Ω_HI·h₇₀ (ALFALFA α.100, Jones u. a. 2018, Fehler statistisch und
# systematisch), molekular Ω_H₂·h (xCOLD GASS, Fletcher u. a. 2021)
KALTES_GAS = dict(hi=(3.9e-4, 0.1e-4, 0.6e-4), h2=(5.34e-5, 0.47e-5))


def anteil_kaltes_gas(omega_b_=None, H0=None, Y=None):
    """Anteil der Baryonen im kalten Gas heute, F_kalt = (Ω_HI + Ω_H₂)/((1 − Y)·Ω_b), mit Helium im Anteil Y (sonst Y_p) und
    der Zählung von KALTES_GAS auf H₀ umgerechnet (Ω_HI ∝ 1/h₇₀, Ω_H₂ ∝ 1/h). Gibt (F_kalt, Fehler)."""
    omega_b_ = omega_b if omega_b_ is None else omega_b_
    H0 = H0_kms if H0 is None else H0
    Y = YP_GEMESSEN[0] if Y is None else Y
    h = H0 / 100
    hi, s_stat, s_sys = KALTES_GAS["hi"]
    h2, s_h2 = KALTES_GAS["h2"]
    Ob = omega_b_ / h**2
    n = (1 - Y) * Ob
    F = (hi / (h / 0.7) + h2 / h) / n
    sF = math.hypot(math.hypot(s_stat, s_sys) / (h / 0.7), s_h2 / h) / n
    return F, sF


def anteil_gesehen(dichte_de):
    """Fester Anteil x·ρ_s, der denselben Abstand zur letzten Streufläche gibt wie dichte_de(a): so liest die Messung."""
    from scipy.optimize import brentq
    D = abstand_mitbewegt(dichte_de)
    rho_s = dichte_vakuumteil()
    return brentq(lambda x: abstand_mitbewegt(lambda a: x * rho_s) - D, 0.5, 5.0, xtol=1e-10)


def weltalter(dichte_de):
    """Weltalter t₀ = ∫ d ln a / H, flach, Dichten heute aus ω_b, ω_c, T_CMB."""
    from scipy.integrate import quad
    rho_r, _, rho_m = dichte_heute()
    H = lambda a: math.sqrt(8 * math.pi * G / 3 * (rho_r / a**4 + rho_m / a**3 + dichte_de(a)))
    return quad(lambda x: 1 / H(math.exp(x)), math.log(1e-12), 0.0, limit=400, epsabs=0, epsrel=1e-10)[0]


def rueckkehr_dichte(T, a0_=a0):
    """Dichte erster Rückkehrzeiten ohne Verlust: f(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/4T}; Gipfel bei τ_D/6."""
    tD = rueckkehrzeit(a0_)
    return math.sqrt(tD / (4 * math.pi)) * T**-1.5 * math.exp(-tD / (4 * T))


def rueckkehr_angekommen(T, a0_=a0):
    """Bis zur Zeit T angekommener Anteil: erfc(√(τ_D/4T))."""
    return math.erfc(math.sqrt(rueckkehrzeit(a0_) / (4 * T)))


def bindungsmasse(sigma, a0_=a0):
    """Kleinste Masse, die Gas der Streuung σ im Tiefbereich bindet: M = (81/4)·σ⁴/(G·a₀) (Umkehr von σ⁴ = (4/81)GMa₀)."""
    return 81 / 4 * sigma**4 / (G * a0_)


def z_gleichheit():
    """Gleichheit von Strahlung und Faltung (mit 1 + Q gelesen): 1 + z_eq = ω_m/ω_r."""
    rho_r, _, rho_m = dichte_heute()
    return rho_m / rho_r - 1


def spektrum_neigung(M, form="BBKS"):
    """Neigung n des Spektrums auf der Masse M (gesamt, mit 1 + Q gelesen): σ²(M) ∝ M^{−(n+3)/3}.
    Nur die Form zählt, keine Normierung; Übergangsfunktion nach BBKS oder nach Eisenstein-Hu ohne Baryonen."""
    from scipy.integrate import quad
    h = H0_kms / 100
    om = omega_b + omega_c
    rho_m = om * rho_krit100
    Gam = om / h                                    # Ω_m·h
    th = T_CMB / 2.7

    def T(k):                                       # k in h/Mpc
        if form == "BBKS":
            q = k / Gam
            return math.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q)**2 + (5.46 * q)**3 + (6.71 * q)**4)**-0.25
        q = k * th**2 / Gam
        L0 = math.log(2 * math.e + 1.8 * q)
        return L0 / (L0 + (14.2 + 731 / (1 + 62.5 * q)) * q * q)

    def sig2(R):                                    # R in Mpc/h
        W = lambda x: 3 * (math.sin(x) - x * math.cos(x)) / x**3
        return quad(lambda lk: math.exp(4 * lk) * T(math.exp(lk))**2 * W(math.exp(lk) * R)**2,
                    math.log(1e-4), math.log(1e4), limit=800)[0]
    R = lambda m: (m * M_sun / (4 * math.pi / 3 * rho_m))**(1 / 3) / Mpc * h
    e = 1.05
    dl = (math.log(sig2(R(M * e))) - math.log(sig2(R(M / e)))) / (2 * math.log(e))
    return -3 - 3 * dl


def grenze_verschmelzen(form="BBKS"):
    """Masse, bei der die Neigung −2 ist: darunter senkt Verschmelzen das inaktive Maß je Baryon, darüber hebt es."""
    from scipy.optimize import brentq
    return 10**brentq(lambda lm: spektrum_neigung(10**lm, form) + 2, 9, 15, xtol=1e-4)


def haltegrenze(v_flach):
    """Größte Streuung, die ein Bereich im Tiefbereich als Gas hält: σ = (4/81)^{1/4}·v_flach = 0,47·v_flach."""
    return (4 / 81)**0.25 * v_flach


# ---------------------------------------------------------------- 03 Dunkle Energie aus der Triade des Ganzen
B1_ZELLE = 3                      # harmonische Richtungen (Schleifen) des Torus T³


def dichte_dunkle_energie(a0_=a0, b1=B1_ZELLE):
    """ρ_DE = b₁·ρ_s = b₁·π³a₀²/(30Gc²): jede Schleife des Ganzen trägt die volle Rücklaufsumme der Faltungsmode."""
    return b1 * dichte_vakuumteil(a0_)


def atmen_heute(a0_=a0, b1=B1_ZELLE, Q=None, wb=None, wnu=0.0):
    """H₀ und Ω_DE aus ω_b, ω_c (oder Q·ω_b), T_CMB und a₀ allein (flach, geschlossen): H₀² = (8πG/3)(ρ_r + ρ_m + ρ_DE).
    wb: ω_b statt des gemessenen (dann ω_c = Q·wb); wnu: ω der Neutrinomasse."""
    rho_r, _, rho_m = dichte_heute(Q)
    if wb is not None:
        rho_m = (wb * (1 + Q) + wnu) * rho_krit100
    rho_de = dichte_dunkle_energie(a0_, b1)
    H = math.sqrt(8 * math.pi * G / 3 * (rho_r + rho_m + rho_de))
    return H * Mpc / 1e3, rho_de / (rho_r + rho_m + rho_de)


def a0_aus_atmen(b1=B1_ZELLE):
    """a₀, für das b₁·ρ_s die Dunkle Energie nach Planck gibt (ω_DE aus H₀, Ω_b, Ω_c, ω_r)."""
    h2 = (H0_kms / 100)**2
    omega_de = (1 - Omega_b - Omega_c) * h2 - dichte_heute()[0] / rho_krit100
    return a0 * math.sqrt(omega_de * rho_krit100 / dichte_dunkle_energie(a0, b1))


def z_entkopplung(Q=None):
    """Rotverschiebung der letzten Streufläche aus ω_b und ω_m (Anpassung von Hu und Sugiyama 1996); die Mode ξ wirkt
    dort linear wie Staub, also gilt die Standardrechnung."""
    wb = omega_b
    wm = omega_b + (omega_c if Q is None else Q * omega_b)
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763)
    g2 = 0.560 / (1 + 21.1 * wb**1.81)
    return 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)


def schallhorizont_stern(dichte_de, Q=None):
    """Mitbewegter Schallhorizont an der letzten Streufläche, r_* = ∫ c_s/(aH) d ln a, c_s = c/√(3(1 + R)),
    R = 3ρ_b·a/(4ρ_γ) [m]."""
    from scipy.integrate import quad
    rho_r, rho_b, rho_m = dichte_heute(Q)
    rho_g = dichte_photonen()
    H = lambda a: math.sqrt(8 * math.pi * G / 3 * (rho_r / a**4 + rho_m / a**3 + dichte_de(a)))
    cs = lambda a: c / math.sqrt(3 * (1 + 3 * rho_b * a / (4 * rho_g)))
    a_stern = 1 / (1 + z_entkopplung(Q))
    return quad(lambda x: cs(math.exp(x)) / (math.exp(x) * H(math.exp(x))), math.log(1e-9), math.log(a_stern),
                limit=500, epsabs=0, epsrel=1e-11)[0]


def schallwinkel(dichte_de, Q=None):
    """Winkel des Schallhorizonts an der letzten Streufläche, θ_* = r_*/D_M(z_*)."""
    return schallhorizont_stern(dichte_de, Q) / abstand_mitbewegt(dichte_de, z_entkopplung(Q), Q)


TAU_PLANCK = (0.054, 0.007)           # optische Tiefe der Reionisation (Planck 2018 VI, aus dem Kurzbericht)
A_PLANCK = (1.0, 0.0025)              # Eichung der Planck-Spektren, plik_lite
DH_GEMESSEN = (2.527e-5, 0.030e-5)    # Deuterium D/H (Cooke, Pettini, Steidel 2018)
YP_GEMESSEN = (0.2449, 0.0040)        # Helium Y_p (Aver, Olive, Skillman 2015)


def cmb_spektren(ombh2, omch2, H0, ln10As, ns, tau, lmax=2508, alles=False, wachstum=(), w0=-1.0, wa=0.0):
    """Spektren der Hintergrundstrahlung D_ℓ = ℓ(ℓ + 1)C_ℓ/2π [µK²] von ℓ = 2 bis lmax (TT, TE, EE, gelinst). Die Mode ξ wirkt
    linear wie Staub mit ω_c, die Dunkle Energie als Formteil (w = −1); gerechnet mit CAMB, N_eff und eine Neutrinomasse
    von 0,06 eV wie bei Planck. Die Anfangsamplitude A_s und die Neigung n_s sind Inhalt, keine Form des Modells.
    alles: dazu das Linsenpotential [L(L + 1)]²C_L^φφ/2π von L = 2 bis 2500 und der Schallhorizont r_d [Mpc];
    wachstum: Rotverschiebungen (fallend), für die dazu f·σ₈ und σ₈ heute kommen (lineares Wachstum);
    w₀, w_a: Gang der Dunklen Energie, w(a) = w₀ + w_a(1 − a) (PPF), w = −1 der Formteil allein."""
    import camb
    extra = dict(redshifts=list(wachstum) + [0.0], kmax=2.0) if wachstum else {}
    if (w0, wa) != (-1.0, 0.0):
        extra.update(w=w0, wa=wa, dark_energy_model="ppf")
    par = camb.set_params(ombh2=ombh2, omch2=omch2, H0=H0, As=math.exp(ln10As) * 1e-10, ns=ns, tau=tau, mnu=0.06,
                          lmax=lmax + 300, lens_potential_accuracy=1, **extra)
    erg = camb.get_results(par)
    cl = erg.get_cmb_power_spectra(par, CMB_unit="muK", raw_cl=False)["total"]
    aus = (cl[2:lmax + 1, 0], cl[2:lmax + 1, 3], cl[2:lmax + 1, 1])
    if alles:
        aus += (erg.get_lens_potential_cls(2500)[2:, 0], erg.get_derived_params()["rdrag"])
    if wachstum:
        aus += (erg.get_fsigma8()[:-1], erg.get_sigma8_0())
    return aus


def scherung_cl(ombh2, omch2, H0, ln10As, ns, nz, w0=-1.0, wa=0.0, ell=None):
    """Scherung der Bilder durch die Faltung längs des Lichtwegs, Limber: C_ℓ^{ij} für die Paare der Rotverschiebungsbins
    (nz = (z, [n_1(z), …]) als Histogramm). Die Mode wirkt wie druckfreier Staub (nichtlinear nach CAMB/HMcode 2020 ohne
    Rückwirkung der Baryonen), der Gang der Dunklen Energie wie bei cmb_spektren. Getrennt: GG (Faltung), GI und II (Ausrichtung
    der Bilder am Feld nach dem linearen Modell, je Einheit der Amplitude A, C1·ρ_krit = 0,0134). Gibt (ℓ, {(i, j): (GG, GI,
    II)}, S8)."""
    import camb
    import numpy as np
    z, ns_ = nz
    extra = dict(w=w0, wa=wa, dark_energy_model="ppf") if (w0, wa) != (-1.0, 0.0) else {}
    par = camb.set_params(ombh2=ombh2, omch2=omch2, H0=H0, As=math.exp(ln10As) * 1e-10, ns=ns, mnu=0.06,
                          halofit_version="mead2020", **extra)
    par.set_matter_power(redshifts=[0.0], kmax=2.0)
    erg = camb.get_results(par)
    S8 = erg.get_sigma8_0() * math.sqrt((ombh2 + omch2 + 0.06 / 93.14) / (H0 / 100)**2 / 0.3)
    P = camb.get_matter_power_interpolator(par, nonlinear=True, hubble_units=False, k_hunit=False, kmax=60.0,
                                           zmax=float(z.max()) + 0.1).P
    Pl = camb.get_matter_power_interpolator(par, nonlinear=False, hubble_units=False, k_hunit=False, kmax=60.0,
                                            zmax=float(z.max()) + 0.1).P
    zz = np.linspace(1e-3, float(z.max()), 600)
    chi = erg.comoving_radial_distance(zz)                         # Mpc
    dchi_dz = 299792.458 / erg.hubble_parameter(zz)
    om = (ombh2 + omch2 + 0.06 / 93.14) / (H0 / 100)**2
    D = np.sqrt(Pl(zz, np.full_like(zz, 0.01), grid=False) / Pl(0.0, 0.01))
    n = [np.interp(zz, z, ni, left=0, right=0) for ni in ns_]
    n = [ni / np.trapezoid(ni, zz) for ni in n]                      # je dz normiert
    vor = 1.5 * om * (H0 / 299792.458)**2
    q = [vor * chi * (1 + zz) * np.array([np.trapezoid(ni[k:] * (chi[k:] - chi[k]) / chi[k:], zz[k:]) for k in range(len(zz))])
         for ni in n]
    nchi = [ni / dchi_dz for ni in n]
    F = -0.0134 * om / D                                           # Ausrichtung je Einheit von A
    ell = np.geomspace(1, 1e5, 300) if ell is None else ell
    kk = (ell[:, None] + 0.5) / chi[None, :]
    Pk = np.where(kk < 60.0, P(np.broadcast_to(zz, kk.shape).ravel(), np.minimum(kk, 60.0).ravel(), grid=False)
                  .reshape(kk.shape), 0.0)                         # jenseits k = 60/Mpc trägt kein Bin mehr
    aus = {}
    for i in range(len(n)):
        for j in range(i, len(n)):
            w = dchi_dz / chi**2
            gg = np.trapezoid(Pk * q[i] * q[j] * w, zz, axis=1)
            gi = np.trapezoid(Pk * F * (q[i] * nchi[j] + nchi[i] * q[j]) * w, zz, axis=1)
            ii = np.trapezoid(Pk * F**2 * nchi[i] * nchi[j] * w, zz, axis=1)
            aus[(i, j)] = (gg, gi, ii)
    return ell, aus, S8


def scherung_xi(ell, cl, theta):
    """ξ₊ und ξ₋ (θ in Bogenminuten) aus C_ℓ: ξ± = ∫ ℓ dℓ/(2π)·C_ℓ·J₀/₄(ℓθ)."""
    import numpy as np
    from scipy.special import jv
    L = np.geomspace(ell[0], ell[-1], 40000)
    c = np.exp(np.interp(np.log(L), np.log(ell), np.log(np.maximum(np.abs(cl), 1e-300)))) * np.sign(np.interp(np.log(L), np.log(ell), cl))
    t = np.atleast_1d(theta) * math.pi / 180 / 60
    return (np.array([np.trapezoid(L * L * c * jv(0, L * x), np.log(L)) for x in t]) / (2 * math.pi),
            np.array([np.trapezoid(L * L * c * jv(4, L * x), np.log(L)) for x in t]) / (2 * math.pi))


def daempfung_zelle(ombh2, omch2, H0, ns, tau, L, w0=-1.0, wa=0.0, lmax=29, kc_chi=200.0, sachs_wolfe=False):
    """Verhältnis C_ℓ(Zelle)/C_ℓ(offen) der Temperatur für ℓ = 2 bis lmax in der geschlossenen Würfelzelle der Kante L
    [Mpc], mitbewegt heute. Die Zelle lässt nur die Moden k = 2π·n/L (n ganzzahlig, n ≠ 0) zu, keine Dichte, eine Form:
    C_ℓ = (2π)³/L³·Σ_n P(k)Δ_ℓ(k)²/k³ statt 4π∫P(k)Δ_ℓ(k)² dk/k; über k_c = kc_chi/χ* das Integral. Δ_ℓ(k) aus CAMB im
    Punkt (Sachs-Wolfe, früher und später Takt, Doppler); sachs_wolfe: nur Δ_ℓ = j_ℓ(kχ*) (zweiter Weg, nur Geometrie).
    Gemittelt über die Lage des Beobachters und die Richtung der Zelle gilt das für ⟨C_ℓ⟩ ohne Bezug auf die Achsen.
    Gibt (ℓ, Verhältnis, χ*)."""
    import camb
    import numpy as np
    from scipy.interpolate import CubicSpline
    from scipy.signal import fftconvolve
    from scipy.special import spherical_jn
    ls = np.arange(2, lmax + 1)
    if sachs_wolfe:
        chi = float(sachs_wolfe)
        kk = np.geomspace(1e-7, 2.0, 400000)
        tr = np.array([spherical_jn(l, kk * chi) for l in ls])
    else:
        extra = dict(w=w0, wa=wa, dark_energy_model="ppf") if (w0, wa) != (-1.0, 0.0) else {}
        par = camb.set_params(ombh2=ombh2, omch2=omch2, H0=H0, As=2e-9, ns=ns, tau=tau, mnu=0.06, lmax=lmax + 100, **extra)
        par.Accuracy.lSampleBoost, par.Accuracy.IntkAccuracyBoost = 50, 4
        erg = camb.get_transfer_functions(par)
        t = erg.get_cmb_transfer_data("scalar")
        chi = erg.conformal_time(0) - erg.tau_maxvis
        kk, tr = t.q, np.array([t.delta_p_l_k[0, list(t.L).index(l)] for l in ls])
    P = (kk / 0.05)**(ns - 1)
    kc = kc_chi / chi
    offen = np.trapezoid(P * tr**2 / kk, kk, axis=1)
    o = kk >= kc
    rest = np.trapezoid(P[o] * tr[:, o]**2 / kk[o], kk[o], axis=1)
    N = int(kc * L / (2 * math.pi)) + 1
    eins = np.bincount(np.arange(-N, N + 1)**2, minlength=N * N + 1).astype(float)
    zahl = np.rint(fftconvolve(fftconvolve(eins, eins), eins)[1:N * N + 1])   # Gitterpunkte je |n|² = 1 … N²
    k = 2 * math.pi * np.sqrt(np.arange(1, len(zahl) + 1)) / L
    g = (zahl > 0) & (k < kc)
    k, zahl = k[g], zahl[g]
    Pk = (k / 0.05)**(ns - 1)
    trk = CubicSpline(kk, tr, axis=1)(k)
    zelle = 2 * math.pi**2 / L**3 * np.sum(zahl * Pk * trk**2 / k**3, axis=1)
    return ls, (zelle + rest) / offen, chi


def gang_dunkle_energie(w0, wa, z):
    """ρ_DE(z)/ρ_DE(0) im Gang w = w₀ + w_a·z/(1 + z)."""
    a = 1 / (1 + z)
    return math.exp(-3 * (1 + w0 + wa) * math.log(a) - 3 * wa * (1 - a))


def maximum_dunkle_energie(w0, wa):
    """Rotverschiebung, bei der die Dunkle Energie im Gang w = w₀ + w_a·z/(1 + z) ihr Maximum hat (w = −1)."""
    return (-1 - w0) / (wa + 1 + w0)


def dichte_de_planck():
    """Dunkle Energie in Plancks bestem Modell (H₀, Ω_b, Ω_c, ω_r, flach), als Masse je Volumen."""
    h2 = (H0_kms / 100)**2
    return ((1 - Omega_b - Omega_c) * h2 - dichte_heute()[0] / rho_krit100) * rho_krit100


def a0_aus_schallwinkel(Q=None, b1=B1_ZELLE):
    """a₀, bei dem b₁·ρ_s(a₀) den Schallwinkel aus Plancks bestem Modell trifft; Q = None: ω_c gemessen, sonst ω_c = Q·ω_b.
    Mit dem gemessenen ω_c ist das das Atmen nach Planck."""
    from scipy.optimize import brentq
    th = schallwinkel(lambda a: dichte_de_planck())
    return brentq(lambda x: schallwinkel(lambda a: dichte_dunkle_energie(x * 1e-10, b1), Q) - th, 0.8, 1.4,
                  xtol=1e-7) * 1e-10


def q_je_radiant():
    """Zählung je Radiant: Das Ganze ist ein voller Umlauf, ν_kosm = 2π, der offene Teil ein Radiant; Q = 2π − 1."""
    return 2 * math.pi - 1


def q_messung():
    """Q = ω_c/ω_b mit 1σ aus den Fehlern von ω_b und ω_c."""
    Qm = omega_c / omega_b
    return Qm, Qm * math.hypot(omega_c_fehler / omega_c, omega_b_fehler / omega_b)


def fuellung(n_max=8):
    """Füllfolge der Hülle nach n + l (bei gleichem n + l zuerst kleineres n): (n, l, erstes Z, letztes Z)."""
    unter = sorted(((n + l, n, l) for n in range(1, n_max + 1) for l in range(n)), key=lambda x: (x[0], x[1]))
    z, aus = 0, []
    for _, n, l in unter:
        aus.append((n, l, z + 1, z + 2 * (2 * l + 1)))
        z += 2 * (2 * l + 1)
    return aus


def schliesszahlen_huelle(n_max=7):
    """Schließzahlen der Hülle, wenn die Schalen nach n + l füllen (bei gleichem n + l zuerst kleineres n): Ende jeder p-Schale."""
    unter = sorted(((n + l, n, l) for n in range(1, n_max + 1) for l in range(n)), key=lambda x: (x[0], x[1]))
    summe, enden = 0, []
    for _, n, l in unter:
        summe += 2 * (2 * l + 1)
        if l == 1 or (n, l) == (1, 0):
            enden.append(summe)
    return tuple(enden)


def dirac_1s_abgleich(Z, R, E=-1.0):
    """1s-Zustand (κ = −1) im Feld einer gleichmäßig geladenen Kugel vom Radius R [in λ̄_e], Energie E [mc²]:
    normierte Wronski-Determinante von regulärer Innen- und abfallender Außenlösung (null = Eigenwert)."""
    from scipy.integrate import solve_ivp
    z = Z * alpha
    V = lambda r: -z * (3 - r * r / (R * R)) / (2 * R) if r < R else -z / r
    f = lambda r, y: [y[0] / r + (E - V(r) + 1) * y[1], -y[1] / r - (E - V(r) - 1) * y[0]]
    rm = max(R, 0.3 / z)                                          # Abgleichsort außerhalb des Kerns
    si = solve_ivp(f, [R * 1e-6, rm], [R * 1e-6, 0.0], rtol=1e-12, atol=1e-30, method="DOP853")
    rb = 4000.0 if E <= -1.0 else 40.0 / math.sqrt(1 - E * E)
    if E <= -1.0:
        q = (-math.sqrt(2 * z / rb) - 1 / rb) / (z / rb)        # G ∝ exp(−√(8zr)), F aus G' = −G/r·κ + (z/r)F
    else:
        q = -math.sqrt((1 - E) / (1 + E))                     # freier Abfall
    so = solve_ivp(f, [rb, rm], [1e-30, q * 1e-30], rtol=1e-12, atol=1e-300, method="DOP853")
    gi, fi, go, fo = si.y[0, -1], si.y[1, -1], so.y[0, -1], so.y[1, -1]
    return (gi * fo - fi * go) / (math.hypot(gi, fi) * math.hypot(go, fo))


def z_kritisch(R_von_Z, E=-1.0, lo=150.0, hi=190.0):
    """Kernladung, bei der der 1s-Zustand die Energie E [mc²] erreicht (E = −1: Bindung 2mc²); R_von_Z(Z) in λ̄_e."""
    from scipy.optimize import brentq
    return brentq(lambda Z: dirac_1s_abgleich(Z, R_von_Z(Z), E), lo, hi, xtol=1e-3)
