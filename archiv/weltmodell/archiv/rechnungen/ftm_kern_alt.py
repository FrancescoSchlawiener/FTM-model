# -*- coding: utf-8 -*-
"""
FTM-Kern: alle Formeln des Modells an einer Stelle, mit einheitlichen Konstanten.
Jede Rechnung importiert diesen Kern (from ftm_kern import *), statt Formeln neu zu schreiben.
Bezeichnungen und Standardwerte: Notation in weltmodell/forschung/00-forschungsstand.md (die eine Deklaration).
Nach jeder Änderung am Modell: python ftm_pruefstand.py (rechnet die Zahlen aus modell/01–09 nach).

Abschnitte
  0 Konstanten und Messwerte (Quelle angegeben; [G] = aus dem Gedächtnis, nicht nachgeprüft)
  1 Grundformeln (01): Neigung, Faltung, Takt
  2 Schwelle und Regel (03, 07): r*, s, q, ν, χ, örtliche Regel, äußeres Feld, Bose-Leiter
  3 Rückkehr (07): Verlust, Rückkehrzeit, Schließlänge Λ
  4 Kosmos (03, 04): Hintergrund H(z), Abstände, Vakuumteil, Λ als ganze innere Faltung, Schließungsregel, Gedächtnis, Kreislauf
  5 Licht, Rotation, Wellen (08)
  6 Teilchen und Sterne (02, 09)
  7 Bahnen: Integrator für die örtliche Regel in einem achsensymmetrischen Newton-Feld
"""
import math
import numpy as np

# ------------------------------------------------------------------ 0 Konstanten (SI)
G = 6.67430e-11            # CODATA 2018
c = 299792458.0
hbar = 6.62607015e-34 / (2 * math.pi)   # h exakt (SI 2019)
kB = 1.380649e-23
e_el = 1.602176634e-19
m_p = 1.67262192369e-27
m_e = 9.1093837015e-31
alpha = 1 / 137.035999084
r_e = 2.8179403262e-15
a_B = 5.29177210903e-11
R_p = 0.8414e-15           # Ladungsradius des Protons (CODATA 2018)

GM_sun = 1.32712440018e20  # IAU
M_sun = GM_sun / G
R_sun = 6.957e8            # IAU nominal
L_sun = 3.828e26           # IAU 2015 B3 nominal
TEFF_sun = 5772.0          # IAU 2015 B3 nominal
GM_earth = 3.986004418e14
M_earth = GM_earth / G
R_earth = 6.371e6
g_earth = 9.81             # Wert der Formel (A1)

AU = 1.495978707e11
pc = 3.0856775814913673e16
kpc, Mpc, Gpc = 1e3 * pc, 1e6 * pc, 1e9 * pc
ly = 9.4607304725808e15
yr = 3.15576e7             # julianisches Jahr
Gyr = 1e9 * yr

l_P = math.sqrt(hbar * G / c**3)
m_Pl = math.sqrt(hbar * c / G)

# Messwerte
A0 = 1.1e-10               # a₀ Standardwert des Modells: SPARC 1,04–1,19 (realistische Fehler), Zelle ohne Kreispaare ≤ 1,13
A0_LIT = 1.2e-10           # a₀ Literaturwert der Beschleunigungsrelation (McGaugh 2016: 1,20 ± 0,02 ± 0,24), nur zum Vergleich
A0_SPARC = (1.04e-10, 1.19e-10)   # SPARC mit realistischen Fehlern, Υ_d = 0,5 (eigene Rechnung, Notiz 74 + Gegenprobe)
H0_KMS = 67.4              # Planck 2018 (Abstract)
OMEGA_B_H2 = 0.02237       # Planck 2018 [G]
OMEGA_C_H2 = 0.1200        # Planck 2018 [G]
OMEGA_BC_H2_FEHLER = (0.00015, 0.0012)   # 1σ von ω_b, ω_c, Planck 2018 [G]
T_CMB = 2.7255
N_EFF_STD = 3.046
N_EFF_MESS = (2.99, 0.17)  # Planck 2018 + BAO (Abstract)
OMEGA_K_MESS = (0.0007, 0.0037)   # Planck 2018 + BAO, 95 %
KREISPAAR_GRENZE = 0.985    # COMPACT 2023 (JCAP 01/030): kürzeste Schleife durch unseren Ort > 0,985·d_LSS (Planck 2015, kubischer Torus: 0,97)


def H0_si(H0=H0_KMS):
    return H0 * 1e3 / Mpc


h_ = H0_KMS / 100
Q_KOSM = OMEGA_C_H2 / OMEGA_B_H2                 # 5,364
q_KOSM = Q_KOSM / (1 + Q_KOSM)                   # 0,843
OMEGA_B = OMEGA_B_H2 / h_**2
OMEGA_M = OMEGA_B * (1 + Q_KOSM)                 # 0,3134
OMEGA_R = 9.14e-5                                 # Photonen + Neutrinos
OMEGA_DE = 1 - OMEGA_M - OMEGA_R
K_S = math.pi**2 / 45                             # Vakuumteil bei Θ = 1


def a0_modell(H0=H0_KMS):
    """a₀ = cH/2π (heutiger Wert der Schwelle nach dem Modell)."""
    return c * H0_si(H0) / (2 * math.pi)


# ------------------------------------------------------------------ 1 Grundformeln (01)
def neigung_xy(x, y):
    """(A1) g(x, y) = x/y²·g⊕, x = M/M⊕, y = R/R⊕."""
    return x / y**2 * g_earth


def ds(M):
    """(A2) Faltungstiefe Δs = GM/(3c²)."""
    return G * M / (3 * c**2)


def neigung(M, R):
    return G * M / R**2


def takt(M, r):
    """(A3) dτ/dt = √(1 − 6Δs/r)."""
    return math.sqrt(1 - 6 * ds(M) / r)


def delta(M, r):
    return 1 - takt(M, r)


def r0(M):
    return 6 * ds(M)


# ------------------------------------------------------------------ 2 Schwelle und Regel
def r_stern(M, a0=A0):
    return math.sqrt(G * M / a0)


def s_(aN, a0=A0):
    """s = √(Γ·τ_D) mit Γ = a_N/c: Exponent der einen Rückkehr ohne Ausdehnung (rueckkehr_exponent)."""
    return rueckkehr_exponent(0.0, np.asarray(aN) / c, 0.0, a0)


def q_(aN, a0=A0):
    return np.exp(-s_(aN, a0))


def nu(aN, a0=A0):
    """Verstärkung ν = 1/(1 − e^{−s}) = 1 + n_B(s)."""
    return 1 / (-np.expm1(-s_(aN, a0)))


def chi(aN, a0=A0):
    s = s_(aN, a0)
    return s / np.expm1(s)


def a_regel(aN, a0=A0):
    """Örtliche Regel, Betrag: a = ν(a_N)·a_N."""
    return nu(aN, a0) * aN


def v_kreis_punkt(M, r, a0=A0):
    aN = G * M / r**2
    return math.sqrt(a_regel(aN, a0) * r)


def v_flach(M, a0=A0):
    return (G * M * a0) ** 0.25


def nu_aussenfeld(a_i, a_e, lage, a0=A0):
    """ν bei äußerem Feld: senkrecht |g| = √(a_i² + a_e²) (s⁴ = s_i⁴ + s_e⁴), gleichgerichtet |g| = a_i + a_e."""
    g = math.hypot(a_i, a_e) if lage == "senkrecht" else a_i + a_e
    return float(nu(g, a0))


def Q_tilde(s):
    """Strukturwirkung Q̃(s) = −4∫_s^∞ x³ n_B(x) dx = −4[s³Li₁ + 3s²Li₂ + 6s·Li₃ + 6Li₄](e^{−s}) (07); Q̃(0) = −4π⁴/15."""
    from scipy.integrate import quad
    return -4 * quad(lambda x: x**3 * math.exp(-x) / -math.expm1(-x) if x > 0 else 0.0, s, np.inf)[0]


def n_B(s):
    return 1 / np.expm1(s)


def bose_stufe_max(d):
    """Maximum von s^d·n_B(s): d(e^s − 1) = s·e^s; gibt (s, a_N/a₀, Wert)."""
    from scipy.optimize import brentq
    s = brentq(lambda x: d * np.expm1(x) - x * math.exp(x), 0.1, 20)
    return s, s * s, s**d * n_B(s)


def sigma_kollektiv(a0=A0):
    """Grenzflächendichte a₀/(πG) in M☉/pc²."""
    return a0 / (math.pi * G) * pc**2 / M_sun


# ------------------------------------------------------------------ 3 Rückkehr (07)
def verlust_gamma(aN):
    return aN / c


def tau_D(a0=A0):
    return c / a0


def rueckkehr_exponent(p, Gamma, theta=0.0, a0=A0):
    """Die eine Rückkehrgleichung ∂_t u + (Γ + θ)u = D∂_ϖ²u (07), am Ort projiziert:
    û(ϖ_R) = e^{−√((p + Γ + θ)·τ_D)}. Gibt den Exponenten √((p + Γ + θ)·τ_D) (p komplex erlaubt).
    p = θ = 0: Regel (q = e^{−s}); Γ = 0: Rest 1 − w des Gewichts; p ≠ 0: Kern K̂(p)."""
    x = (p + Gamma + theta) * tau_D(a0)
    return np.sqrt(x)


def rueckkehr(p, Gamma, theta=0.0, a0=A0):
    """û = e^{−√((p + Γ + θ)·τ_D)}: Rückkehranteil q, Rest 1 − w und Kern K̂ in einem Ausdruck."""
    return np.exp(-rueckkehr_exponent(p, Gamma, theta, a0))


def schliesslaenge(a0=A0):
    """Λ = c²/(2a₀)."""
    return c**2 / (2 * a0)


def ell0(a0=A0):
    """ℓ₀ = 2Λ = c²/a₀; Leiter r*² = 3Δs·ℓ₀."""
    return 2 * schliesslaenge(a0)


def a0_aus_schliesslaenge(Lam):
    return c**2 / (2 * Lam)


# ------------------------------------------------------------------ 4 Kosmos
def E_hintergrund(z, Om=OMEGA_M, Or=OMEGA_R):
    """Hintergrund mit festem Vakuumteil und fester Rückstellung (heute w = −1)."""
    a = 1 / (1 + z)
    return math.sqrt(Om * a**-3 + Or * a**-4 + (1 - Om - Or))


def abstand_mitbewegt(z, H0=H0_KMS, n=20000):
    zs = np.linspace(0, z, n)
    Ez = np.array([E_hintergrund(x) for x in zs]) if n < 2000 else np.sqrt(OMEGA_M * (1 + zs)**3 + OMEGA_R * (1 + zs)**4 + (1 - OMEGA_M - OMEGA_R))
    return c / H0_si(H0) * np.trapezoid(1 / Ez, zs)


def d_LSS(z=1090.0):
    """Durchmesser der letzten Streufläche (mitbewegt)."""
    zs = np.concatenate([np.linspace(0, 10, 20001), np.geomspace(10, z, 20001)[1:]])
    Ez = np.sqrt(OMEGA_M * (1 + zs)**3 + OMEGA_R * (1 + zs)**4 + (1 - OMEGA_M - OMEGA_R))
    return 2 * c / H0_si() * np.trapezoid(1 / Ez, zs)


def theta_(a0, H):
    return 2 * math.pi * a0 / (c * H)


def omega_vakuum(a0=A0, H0=H0_KMS):
    """Vakuumteil Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH)."""
    return K_S * theta_(a0, H0_si(H0))**2


def lambda_vakuum(Lam):
    """Λ_s = π⁴/(15Λ²) = 4π⁴/(15ℓ₀²), ℓ₀ = 2Λ."""
    return math.pi**4 / (15 * Lam**2)


def H_vakuum(Lam):
    """Atemrate des Vakuumteils allein: H_s² = Λ_s·c²/3."""
    return math.sqrt(lambda_vakuum(Lam) * c**2 / 3)


def rho_m_mitbewegt():
    return OMEGA_M * 3 * H0_si()**2 / (8 * math.pi * G)


def a0_aus_zelle(L):
    """Λ als ganze innere Faltung der Zelle: Λ = Gρ_mL³/(3c²), a₀ = 3c⁴/(2Gρ_mL³)."""
    return 3 * c**4 / (2 * G * rho_m_mitbewegt() * L**3)


def zelle_aus_a0(a0):
    return (3 * c**4 / (2 * G * rho_m_mitbewegt() * a0)) ** (1 / 3)


def a0_frueh_relativ(z):
    """Strahlung erhöht die innere Faltung: Λ(a) = Λ_m(1 + Ω_r/(Ω_m a))."""
    a = 1 / (1 + z)
    return 1 / (1 + OMEGA_R / (OMEGA_M * a))


def delta_N_eff(f):
    """Fester Anteil f am Gesamtbudget im Strahlungszeitalter als ΔN_eff."""
    return f / (1 - f) * (1 + 0.2271 * N_EFF_STD) / 0.2271


def gewicht_w_preisach(theta_geschichte, a0=A0):
    """Schließungsregel mit Gedächtnis: w(θ̃), θ̃ = max(min θ der Geschichte, 0)."""
    return gewicht_w(min(theta_geschichte), a0)


def gewicht_w(theta, a0=A0):
    """Schließungsregel aus derselben Rückkehr wie die Regel: 1 − w = E[e^{−θ₊T}] über die Dichte
    erster Rückkehrzeiten, also w = 1 − e^{−√(θ₊/θ_c)}, θ_c = a₀/c = 1/τ_D."""
    th = max(theta, 0.0)
    return -math.expm1(-float(rueckkehr_exponent(0.0, 0.0, th, a0)))


def gewicht_w_maxwell(theta, a0=A0):
    """Frühere Form (Rechnung, nicht im Modell): Maxwell-Speicheranteil w = θ₊²/(θ₊² + θ_c²)."""
    th = max(theta, 0.0)
    thc = a0 / c
    return th**2 / (th**2 + thc**2)


def uebergabezeit(dichte_faktor, z=0.0, N_stern=(0.3, 2.4)):
    """1/λ_m = N*/ω_orb, ω_orb² = (4πG/3)·Δ·ρ_krit(z); gibt (1/ω_orb, 1/λ_m min, max) in s."""
    Hz = H0_si() * E_hintergrund(z)
    w = math.sqrt(dichte_faktor / 2) * Hz          # (4πG/3)·Δ·3H²/(8πG) = ΔH²/2
    return 1 / w, N_stern[0] / w, N_stern[1] / w


def horizontentropie(H):
    """S_max = πkc²/(ℓ_P²H²) in Einheiten von k."""
    return math.pi * c**2 / (l_P**2 * H**2)


def kreislauf_mindestrate(Lam, sigma_k=7.1e86):
    """Λ̇/Λ ≥ Ṡ_i/(2S_max) im Endzustand, H² = π⁴c²/(45Λ²); Entropieproduktion Ṡ_i (sigma_k) in k/s je Horizont."""
    H_end = c * math.pi**2 / (math.sqrt(45) * Lam)
    S = horizontentropie(H_end)
    return sigma_k / (2 * S), S, H_end


# ------------------------------------------------------------------ 5 Licht, Rotation, Wellen (08)
def ablenkung(M, b):
    return 12 * ds(M) / b


def ablenkung_kollektiv(M, a0=A0):
    """α_∞ = 6πΔs/r* = 2π(v_flach/c)²."""
    return 2 * math.pi * (v_flach(M, a0) / c)**2


def perihel_pro_umlauf(M, a, e):
    return 18 * math.pi * ds(M) / (a * (1 - e * e))


def schatten_winkel(M, D):
    return 2 * 3 * math.sqrt(3) * G * M / c**2 / D


def hulse_taylor(m1, m2, Pb, e):
    M = (m1 + m2) * M_sun
    eta = m1 * m2 / (m1 + m2)**2
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e * e)**3.5
    return -(192 * math.pi / 5) * eta * (2 * math.pi * 3 * ds(M) / (c * Pb))**(5 / 3) * fe


def masse_r_stern_gleich_compton(a0=A0):
    """r* = λ̄ bei M = (ħ²a₀/(Gc²))^{1/3}."""
    return (hbar**2 * a0 / (G * c**2)) ** (1 / 3)


# ------------------------------------------------------------------ 6 Teilchen und Sterne (02, 09)
def compton_red(m):
    return hbar / (m * c)


def alpha_G(M=m_p):
    return G * M**2 / (hbar * c)


def sternmasse():
    return m_Pl**3 / m_p**2


def grenzmasse_hawking(T=T_CMB):
    """Falte mit T_H = T: M = ħc³/(8πGkT)."""
    return hbar * c**3 / (8 * math.pi * G * kB * T)


# ------------------------------------------------------------------ 7 Bahnen
def bahnen_meridian(gR_fn, gz_fn, R0, z0, vR0, vz0, Lz, T, dt, a0=A0, regel=True):
    """Leapfrog in (R, z) mit L_z fest; gR_fn, gz_fn: Newton-Feld (m/s²) als Funktionen von (R, z) in m.
       regel=True: örtliche Regel a = ν(|g_N|)·g_N; False: Newton."""
    R, z, vR, vz = map(np.array, (R0, z0, vR0, vz0))

    def acc(R, z):
        gR, gz = gR_fn(R, z), gz_fn(R, z)
        f = nu(np.hypot(gR, gz), a0) if regel else 1.0
        return f * gR + Lz**2 / R**3, f * gz
    aR, az = acc(R, z)
    for _ in range(int(T / dt)):
        vR = vR + 0.5 * dt * aR; vz = vz + 0.5 * dt * az
        R = R + dt * vR; z = z + dt * vz
        aR, az = acc(R, z)
        vR = vR + 0.5 * dt * aR; vz = vz + 0.5 * dt * az
    return R, z, vR, vz


# ------------------------------------------------------------------ 8 Ergänzungen für den Prüfstand
def isco_kerr(j, mitlaeufig=True):
    """Innerste stabile Kreisbahn (Bardeen) in Einheiten GM/c² = 3Δs."""
    z1 = 1 + (1 - j * j) ** (1 / 3) * ((1 + j) ** (1 / 3) + (1 - j) ** (1 / 3))
    z2 = math.sqrt(3 * j * j + z1 * z1)
    sg = -1 if mitlaeufig else 1
    return 3 + z2 + sg * math.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def realteil_quadrat_punkt(s):
    """(Re τ_M)² = 1 + χ(s) für die Punktmasse, s = r*/r."""
    return 1 + s / math.expm1(s)


def radius_resonanz(p_zu_q, rstern):
    """Radius, an dem Re τ_M = p/q (Punktmasse): 1 + χ(s) = (p/q)²."""
    from scipy.optimize import brentq
    target = p_zu_q**2 - 1
    s = brentq(lambda x: x / math.expm1(x) - target, 1e-6, 50)
    return rstern / s


def hernquist_masse(M, a, r):
    return M * r * r / (r + a) ** 2


def schwellenfrequenz_gw(M, a0=A0):
    """f*_gw = a₀^{3/4}/(π(3c²Δs)^{1/4})."""
    return a0 ** 0.75 / (math.pi * (3 * c**2 * ds(M)) ** 0.25)


def eddington(M):
    """L_Edd = (9/2)·m_p·c³·Δs/r_e²."""
    return 4.5 * m_p * c**3 * ds(M) / r_e**2


def lebensdauer_falte(M):
    """t = 3840·(S/k)·Δs/c mit S/k = 4π·α_G(M)."""
    return 3840 * 4 * math.pi * alpha_G(M) * ds(M) / c


def hawking_lebensdauer(M):
    return 5120 * math.pi * G**2 * M**3 / (hbar * c**4)


# ------------------------------------------------------------------ 9 Labor: eingespannter Kupferdraht im Magnetfeld (06)
# Stoffwerte Kupfer [G]; Draht, Vorspannung, B und Q_m sind Annahmen (wie torus_wechselkraefte.py, Teil 4a)
CU = dict(rho_el=1.68e-8, alR=3.93e-3, E=128e9, alE=-3.8e-4, alth=16.5e-6, rho=8960.0, cp=385.0, kth=401.0)
DRAHT = dict(a=50e-6, ell=0.02, B=0.02, Qm=1000.0, sig0=20e6)


def draht(cu=CU, d=DRAHT):
    """Größen des Drahts (Galerkin-Mode 1 − cos 2πx/ℓ): R0, G_th, C_th, m_eff, Γ = Bℓ, c_m, k(ΔT)."""
    A = math.pi * d['a']**2
    I = math.pi * d['a']**4 / 4
    ell = d['ell']
    p = dict(R0=cu['rho_el'] * ell / A, G_th=12 * cu['kth'] * A / ell, C_th=cu['rho'] * cu['cp'] * A * ell,
             m=1.5 * cu['rho'] * A * ell, Gam=d['B'] * ell, alR=cu['alR'])

    def k(dT):
        Et = cu['E'] * (1 + cu['alE'] * dT)
        N = d['sig0'] * A - cu['E'] * A * cu['alth'] * dT
        return 8 * math.pi**4 * Et * I / ell**3 + 2 * math.pi**2 * N / ell
    p['k'] = k
    p['w0'] = math.sqrt(k(0.0) / p['m'])
    p['c_m'] = p['m'] * p['w0'] / d['Qm']
    return p


# ------------------------------------------------------------------ 10 Naht, Haufen, Q-Fenster (03)
def r_Q(M, a0=A0, qk=q_KOSM):
    """R_Q = r*/(−ln q_kosm): dort erreicht das örtliche q das kosmische, a_N/a₀ = (ln q)²."""
    return r_stern(M, a0) / (-math.log(qk))


def masse_bQ(a0=A0, H=None, Qk=Q_KOSM):
    """M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴): Baryonenmasse, bei der R_Q die Grenze R₂₀₀ erreicht."""
    H = H0_si() if H is None else H
    lq = math.log(Qk / (1 + Qk))
    return (1 + Qk)**2 * lq**6 * a0**3 / (1e4 * G * H**4)


_UMKEHR = {}


def dichte_umkehr_heute(Om=OMEGA_M):
    """1 + Δ_ta: mittlere Dichte (gegen ρ_m) einer Kugelschale, die heute umkehrt (Top-Hat im flachen
    Hintergrund mit Vakuumteil 1 − Ω_m, Start a = 10⁻³ mit Hubble-Fluss (1 − δ/3)). Unabhängig von der Masse."""
    if Om in _UMKEHR:
        return _UMKEHR[Om]
    from scipy.integrate import solve_ivp, quad
    from scipy.optimize import brentq
    H0, OL = H0_si(), 1 - Om
    Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
    t_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
    t0, ai = t_a(1.0), 1e-3
    ti = t_a(ai)

    def lauf(di):   # Einheit: Radius ri = 1, GM = (H_i²/2)·Ω_m(a_i)·(1 + δ)
        Hi = Hb(ai); GM = 0.5 * Hi**2 * (Om / ai**3) / (Om / ai**3 + OL) * (1 + di)
        f = lambda t, y: [y[1], -GM / y[0]**2 + OL * H0**2 * y[0]]
        ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
        sol = solve_ivp(f, (ti, 2 * t0), [1.0, Hi * (1 - di / 3)], events=ev, rtol=1e-9, atol=1e-12)
        return (sol.t_events[0][0], sol.y_events[0][0][0]) if len(sol.t_events[0]) else (2 * t0, None)
    ld = brentq(lambda x: lauf(10**x)[0] - t0, -4, 0, xtol=1e-10)
    di = 10**ld; R = lauf(di)[1]
    _UMKEHR[Om] = (1 + di) / (R**3 * ai**3)     # Masse/Volumen gegen ρ_m heute: (1 + δ)·ρ_m(a_i)·r_i³/R³
    return _UMKEHR[Om]


def umkehrradius_heute(M_ges):
    """R_ta = (3M/(4πρ_m(1 + Δ_ta)))^{1/3}, M_ges = (1 + Q)·M_b (Maß des Ganzen)."""
    return (3 * M_ges / (4 * math.pi * rho_m_mitbewegt() * dichte_umkehr_heute()))**(1 / 3)


def q_fenster(a0=A0, L_min=None):
    """Fenster der kosmischen Rückkehrzahl: unten r*(M_b,Zelle)/Λ = √(6/ν) ≤ 1 ⇒ Q ≥ 5; oben a₀ = 3c⁴/(2Gνρ_bL³)
    mit L ≥ 0,985·d_LSS (keine Kreispaare) ⇒ ν ≤ 3c⁴/(2Gρ_b a₀ L_min³)."""
    L_min = KREISPAAR_GRENZE * d_LSS() if L_min is None else L_min
    rho_b = OMEGA_B * 3 * H0_si()**2 / (8 * math.pi * G)
    return 5.0, 3 * c**4 / (2 * G * rho_b * a0 * L_min**3) - 1



# ------------------------------------------------------------------ 11 Form des Ganzen (Formverhältnisse, zeitfest)
def Q_fehler():
    """σ(Q) aus σ(ω_b), σ(ω_c), unkorreliert: Q·√((σ_c/ω_c)² + (σ_b/ω_b)²)."""
    sb, sc = OMEGA_BC_H2_FEHLER
    return Q_KOSM * math.hypot(sc / OMEGA_C_H2, sb / OMEGA_B_H2)


def im_tau_aus_nome(q):
    """Flacher Torus ℂ/(ℤ+τℤ), Nome q = e^{2πiτ}: Im τ = −ln|q|/2π."""
    return -math.log(q) / (2 * math.pi)


def torus_streifen_nu(R_zu_r, breite=1.0, lage="aussen"):
    """Eingebetteter Torus (R, r), Flächenelement r(R + r cos θ)dθdφ. Offener Teil = Streifen der Breite
    'breite' (Radiant) im Rohrwinkel θ um θ = 0 (außen) bzw. θ = π (innen): Anteil (b·R ± 2r·sin(b/2))/(2πR);
    ν = 1/Anteil = 2πR/(b·R ± 2r·sin(b/2))."""
    sgn = 1.0 if lage == "aussen" else -1.0
    return 2 * math.pi * R_zu_r / (breite * R_zu_r + sgn * 2 * math.sin(breite / 2))
