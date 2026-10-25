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


def takt(M, r, Q=0.0):
    """dτ/dt = √(1 − 6Δs/r): Abwicklungstakt am Ort r. Mit Windung Q liegt ihr gefaltetes Feld außerhalb des Knotens,
    innerhalb von r bleibt M − Q²/(8πε₀c²r), also (dτ/dt)² = 1 − 6Δs/r + r_Q·3Δs/r² mit r_Q = Q²/(4πε₀Mc²) und
    r_Q·3Δs = GQ²/(4πε₀c⁴) = (ℓ_P·Q/q_P)², ohne Masse."""
    return math.sqrt(1 - 6 * faltung(M) / r + G * Q**2 / (4 * math.pi * eps0 * c**4 * r**2))


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


def neigung_verlangt(v_sicht, v_regel):
    """Neigung [rad], unter der eine drehende Scheibe mit der Geschwindigkeit längs der Sicht v_sicht = v·sin i die Regel
    trifft: sin i = v_sicht/v_Regel. Die Neigung liest nur die Geschwindigkeit, nicht die Baryonen und nicht a_N, der Abstand
    von der Regel ist lg(v_sicht/v_Regel) − lg sin i. Über 1 gibt es keine solche Neigung (π/2)."""
    return math.asin(min(1.0, v_sicht / v_regel))


def v_flach(M, a0_=a0):
    """v⁴ = GMa₀: flache Rotationskurve im tiefen Bereich."""
    return (G * M * a0_)**0.25


def theta_(a0_=a0, H=H0):
    """Θ = 2πa₀/(cH)."""
    return 2 * math.pi * a0_ / (c * H)


def scheibe_achsen(v, th, Q=None, lam=None, a0_=a0, numerisch=False):
    """Größe einer Scheibe, die den linearen Drehimpuls ihres Gases trägt, R_d = λ'·V₂₀₀·R₂₀₀/(√2·v_flach) mit
    V₂₀₀·R₂₀₀ = (G²M²/(10H))^{1/3}, M = (1 + Q)·M_b, gegen ihre Schwelle in den Achsen der Leiter:
    R_d/r* = (λ'/√2)·(1 + Q)^{2/3}·(Θ·j/(20π))^{1/3} mit der Spinzahl j = c/v_flach an der Schwelle und Θ = 2πa₀/(cH) zur
    Zeit der Bindung. v [m/s] = v_flach. lam None: der Anteil der Drehung λ', bei dem die Scheibe an der Schwelle der
    Flächendichte liegt, R_d = r*/√2, λ' = (1 + Q)^{−2/3}·(20π/(Θ·j))^{1/3}, sonst R_d/r*. numerisch: zweiter Weg über
    M_b = v⁴/(Ga₀), H aus Θ, R₂₀₀ bei 200-facher kritischer Dichte und r* = √(GM_b/a₀)."""
    Q = q_je_radiant() if Q is None else Q
    if numerisch:
        Mb = v**4 / (G * a0_)
        H = 2 * math.pi * a0_ / (c * th)
        M = (1 + Q) * Mb
        R2 = (G * M / (100 * H * H))**(1 / 3)
        VR = math.sqrt(G * M / R2) * R2
        rs = math.sqrt(G * Mb / a0_)
        return rs * v / VR if lam is None else lam * VR / (math.sqrt(2) * v) / rs
    f = (1 + Q)**(2 / 3) * (th * (c / v) / (20 * math.pi))**(1 / 3) / math.sqrt(2)
    return 1 / (math.sqrt(2) * f) if lam is None else lam * f


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


def leistung_quadrupol(m1, m2, a):
    """Abstrahlung zweier Massen auf der Kreisbahn mit Abstand a über den Faltungsquadrupol:
    P = (32/5)·G⁴(m₁m₂)²(m₁ + m₂)/(c⁵a⁵)."""
    return 32 / 5 * G**4 * (m1 * m2)**2 * (m1 + m2) / (c**5 * a**5)


def lebensdauer_falte(M):
    """Die Falte strahlt als schwarzer Körper ihrer Temperatur über den Horizont 4πr₀², P = σT⁴·4πr₀² ∝ M⁻²,
    also t = Mc²/(3P)."""
    r0 = 6 * faltung(M)
    kT = k_B * falte_temperatur(M)
    P = math.pi**2 / 60 * kT**4 / (hbar**3 * c**2) * 4 * math.pi * r0**2
    return M * c**2 / (3 * P)


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


def zeit_abgabe(l, E, P, drehung=False, k=None, spiegel=False, d=None):
    """Zeit bis zum Ende oder seit dem Anfang einer Abgabe der kleinsten Ordnung l (Windung l = 1 als Dipol, Faltung l = 2
    als Quadrupol) bei der Energie E und der Leistung P. Wächst die Abgabe wie P ∝ E^k, ist t = E/(|k − 1|·P). Bahn: Bindung
    E ∝ 1/r, P ∝ r^{−(l+3)}, k = l + 3, Fall nach t = E/((l + 2)·P). Drehung: E = ½IΩ², P ∝ Ω^{2l+2}, k = l + 1, Bremsindex
    n = 2l + 1, gedreht seit höchstens t = Ω/((n − 1)·|Ω̇|) = E/(l·P). k: der Exponent direkt, für jede Energie, die mit
    P ∝ E^k ab- oder zunimmt, Falte k = −2 (P ∝ M⁻², t = E/(3P) wie lebensdauer_falte), Wachsen an der Grenze k = 1 (E/P je
    Faktor e, ohne Ende und Anfang in endlicher Zeit). spiegel: dieselbe Rate in der Spiegelmasse, Ẽ = m̃c² =
    (m_Pl·c²)²/(2E), Rate (Ẽ/E)·P ∝ Ẽ^{2 − k}. Der Spiegel kehrt k − 1 um und hält t fest, Fixpunkt k = 1, die Falte
    gespiegelt k = 4 (Ṫ ∝ T⁴). d: Bahn im Feld der Flussdimension d, g ∝ r^{1−d}, E die kinetische Energie, Bindung je
    e-Faktor des Radius (4 − d)·E, P ∝ r^{2l}·Ω^{2l+2}, Fall nach t = (4 − d)·E/((l·(d − 2) + 2)·P), d = 3 wie oben,
    d = 2 E/P für jedes l. Gibt t [s], mit spiegel (2 − k, Ẽ [J], Rate von Ẽ [W], t [s])."""
    if d is not None and not drehung:
        return (4 - d) * E / ((l * (d - 2) + 2) * P)
    if k is None:
        k = (l + 1) if drehung else (l + 3)
    t = E / (abs(k - 1) * P) if k != 1 else E / P
    if not spiegel:
        return t
    Et = (m_Pl * c**2)**2 / (2 * E)
    Pt = Et / E * P
    return 2 - k, Et, Pt, (Et / (abs(1 - k) * Pt) if k != 1 else Et / Pt)


def verschmelzen_gleich(m, delta):
    """Zwei gleiche Massen m auf der Kreisbahn mit δ = 3Δs(m)/a: Bindung ½mc²·δ, Leistung (64/5)·(c⁵/G)·δ⁵, bis zur
    Verschmelzung t = E/(4P) = (5/512)·(3Δs/c)·δ⁻⁴ in N = δ^{−5/2}/(64√2·π) Umläufen. Gibt (t [s], N)."""
    return 5 / 512 * 3 * faltung(m) / c * delta**-4, delta**-2.5 / (64 * math.sqrt(2) * math.pi)


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


def paar_einrasten(lam, numerisch=False):
    """Paar gegensinniger Windungen, das über einen Takt des Gitters (ħω_D = 1) mit der Kopplung λ einrastet (BCS als Einrasten
    nach Weiss für ein Paar): 1 = λ·∫₀¹ dξ·(1 − 2f(E/kT))/E, E = √(ξ² + Δ²), f die Fermi-Form zweier Umlaufsinne
    (rueckkehr_summe mit sinne=2). Gibt (Δ₀, kT_c, Δ₀/kT_c) in Einheiten des Takts. Geschlossen, schwach gekoppelt:
    Δ₀ = 1/sinh(1/λ) ≈ 2·e^{−1/λ}, kT_c = (2e^γ/π)·e^{−1/λ}, Δ₀/kT_c = π·e^{−γ}. numerisch: zweiter Weg, beide Integrale
    über die Quadratur mit der Fermi-Form und Nullstellen in Δ₀ und T_c."""
    import numpy as np
    if not numerisch:
        tc = 2 * math.exp(np.euler_gamma) / math.pi * math.exp(-1 / lam)
        return 1 / math.sinh(1 / lam), tc, math.pi * math.exp(-np.euler_gamma)
    from scipy.integrate import quad
    from scipy.optimize import brentq
    null = lambda d: lam * quad(lambda x: 1 / math.hypot(x, d), 0, 1, limit=200)[0] - 1
    d0 = brentq(null, 1e-300 if lam < 0.05 else 1e-12, 1.0, xtol=1e-16)
    tanh_f = lambda x, t: 1 - 2 * rueckkehr_summe(x / t, sinne=2) if x / t < 700 else 1.0
    warm = lambda t: lam * quad(lambda x: tanh_f(x, t) / x, 0, 1, limit=400, points=[t])[0] - 1
    tc = brentq(warm, d0 / 3, d0, xtol=1e-16)
    return d0, tc, d0 / tc


def zweite_rueckkehr_windung(numerisch=False):
    """Zweite Rückkehr der Windung, in (α/π)² gezählt, wie die Quantenelektrodynamik sie rechnet (Petermann, Sommerfield):
    C₂ = 197/144 + ζ(2)/2 − 3ζ(2)·ln 2 + (3/4)·ζ(3), mit ζ(n) = Σ 1/kⁿ der Summe über volle Umläufe und ln 2 = Σ (−1)^{k+1}/k
    der Summe über halbe Drehungen. Gibt (C₂, rationaler Teil, Teil der Umläufe ζ(2)/2 + (3/4)ζ(3), Teil der halben Drehung
    −3ζ(2)·ln 2). numerisch: zweiter Weg, die drei Summen direkt über k (Euler-Maclaurin am Ende, die wechselnde Summe
    gemittelt über zwei Teilsummen)."""
    import numpy as np
    if numerisch:
        k = np.arange(1, 200001, dtype=float)
        N = k[-1]
        z2 = np.sum(1 / k**2) + 1 / N - 1 / (2 * N**2)
        z3 = np.sum(1 / k**3) + 1 / (2 * N**2) - 1 / (2 * N**3)
        s = np.cumsum((-1)**(k + 1) / k)
        l2 = (s[-1] + s[-2]) / 2
    else:
        from scipy.special import zeta
        z2, z3, l2 = math.pi**2 / 6, float(zeta(3)), math.log(2)
    rat, uml, halb = 197 / 144, z2 / 2 + 0.75 * z3, -3 * z2 * l2
    return rat + uml + halb, rat, uml, halb


def hysteresefenster(t):
    """Kippfeld der eingerasteten Windungen, t < 1: h_c = t·artanh(√(1 − t)) − √(1 − t) [kT_C/μ], ≈ (2/3)(1 − t)^{3/2}."""
    w = math.sqrt(1 - t)
    return abs(t * math.atanh(w) - w)


def besetzung(h, J=0.0):
    """Besetzung einer Bindestelle, die höchstens einen Partner trägt (Sperre): Fermi-Form θ = 1/(1 + e^{−(h + J·θ)}) mit
    h = ln([L]/K_d) und der Kopplung J [kT] benachbarter Stellen im mittleren Feld (J < 4: eine Lösung)."""
    from scipy.optimize import brentq
    return brentq(lambda t: t - 1 / (1 + math.exp(-(h + J * t))), 0.0, 1.0, xtol=1e-15)


def hill_zahl(J):
    """Steilheit der Besetzung bei halber Füllung (h = −J/2), n_H = d ln(θ/(1 − θ))/dh = 1/(1 − q) mit q = J/4: die Summe
    der Rückkopplung, Kippstelle bei J = 4."""
    return 1 / (1 - J / 4)


def rueckkehr_membran(x, numerisch=False, punkte=200000, weite=40.0):
    """Rückbindung an ein Feld von Rezeptoren auf einer Membran: Der Wirkstoff kehrt senkrecht zur Membran zurück, die eine
    Rückkehrgleichung D·G″ = Γ·G auf der Halbgeraden, und das Feld fängt bei jeder Rückkehr mit κ = ρ·k_on ein (Rand
    D·G′(0) = −κ·(1 − G(0))). Mit x = κ/√(DΓ) ist der Anteil, der wieder bindet, bevor er verloren ist, q = x/(1 + x).
    numerisch: zweiter Weg, finite Differenzen in Einheiten √(D/Γ) bis weite."""
    if not numerisch:
        return x / (1 + x)
    import numpy as np
    from scipy.linalg import solve_banded
    h = weite / punkte
    n = punkte                                                   # Unbekannte u_0 … u_{n−1}, u_n = 0
    ab = np.zeros((3, n))
    ab[0, 1:], ab[1, :], ab[2, :-1] = 1.0, -(2 + h * h), 1.0
    ab[1, 0], ab[0, 1] = -(1 + x * h + h * h / 2), 1.0            # Robin-Rand zweiter Ordnung
    rhs = np.zeros(n)
    rhs[0] = -x * h
    return float(solve_banded((1, 1), ab, rhs)[0])


def rueckkehr_dimension(d, r_durch_a, ka, numerisch=False, punkte=400000, weite=60.0):
    """Anteil, der aus dem Abstand r zu einer Stelle vom Radius a zurückfindet, bevor er verloren ist, in d Dimensionen:
    u″ + (d − 1)/r·u′ = k²·u, u(a) = 1, u(∞) = 0, k = √(Γ/D). Geschlossen d = 1: e^{−k(r − a)}, d = 2: K₀(kr)/K₀(ka),
    d = 3: (a/r)·e^{−k(r − a)}. numerisch: zweiter Weg, finite Differenzen von a bis weite/k."""
    import numpy as np
    k = ka                                                       # Längen in Einheiten von a
    r = r_durch_a
    if not numerisch:
        if d == 1:
            return math.exp(-k * (r - 1))
        if d == 2:
            from scipy.special import k0
            return float(k0(k * r) / k0(k))
        return math.exp(-k * (r - 1)) / r
    from scipy.linalg import solve_banded
    ende = 1 + weite / k
    x = np.linspace(1, ende, punkte + 1)
    h = x[1] - x[0]
    xi = x[1:-1]
    n = len(xi)
    unter, mitte, ober = 1 / h**2 - (d - 1) / (2 * h * xi), -2 / h**2 - k * k, 1 / h**2 + (d - 1) / (2 * h * xi)
    ab = np.zeros((3, n))
    ab[0, 1:], ab[1, :], ab[2, :-1] = ober[:-1], mitte, unter[1:]
    rhs = np.zeros(n)
    rhs[0] = -unter[0]                                           # u(a) = 1
    u = solve_banded((1, 1), ab, rhs)
    return float(np.interp(r, xi, u))


def rueckkehr_lesart(q, sinne=1, numerisch=False, umkehr=False):
    """Eine Rückkehr, zwei Lesarten: Mit einem Umlaufsinn (Stapeln) summieren sich die Rückläufe zur Bose-Form Σqⁿ =
    q/(1 − q), Kippstelle bei q = 1. Mit zwei Umlaufsinnen (Sperre) kommt jede Rückkehr mit einer halben Drehung zurück,
    Σ(−1)ⁿ⁺¹qⁿ = q/(1 + q), die Besetzung einer Stelle mit q = [L]/K_d, gespiegelt θ(q) + θ(1/q) = 1. Die beiden Lesarten
    kehren einander um. numerisch: zweiter Weg, die Reihe gliedweise summiert, bei zwei Sinnen und q > 1 über den Spiegel.
    umkehr: bei einem Sinn die Lesart mit zwei Sinnen numerisch umgekehrt, x mit x/(1 + x) = q."""
    if umkehr and sinne == 1:
        from scipy.optimize import brentq
        return brentq(lambda x: x / (1 + x) - q, 0.0, 1e12, xtol=1e-14, rtol=1e-15)
    if not numerisch:
        return q / (1 - q) if sinne == 1 else q / (1 + q)
    if sinne == 2 and q > 1:
        return 1 - rueckkehr_lesart(1 / q, 2, True)
    s, t, n = 0.0, q, 1
    while abs(t) > 1e-18:
        s += t if sinne == 1 or n % 2 else -t
        t *= q
        n += 1
    return s


def besetzung_hysterese(J, numerisch=False):
    """Kooperative Besetzung über der Kippstelle, J > 4: Breite des Hysteresefensters in h zwischen den beiden Falten,
    D = 4u/(1 − u²) − 4·artanh u mit u = √(1 − 4/J), nahe der Kippstelle (J − 4)^{3/2}/3. numerisch: zweiter Weg, die
    Lösungen von θ = 1/(1 + e^{−(h + Jθ)}) auf einem feinen Gitter gezählt und die Stellen, an denen aus drei eine wird,
    durch Halbieren gesucht."""
    import numpy as np
    if not numerisch:
        u = math.sqrt(1 - 4 / J)
        return 4 * u / (1 - u * u) - 4 * math.atanh(u)
    t = np.linspace(1e-9, 1 - 1e-9, 200001)
    zahl = lambda h: int(np.sum(np.diff(np.sign(t - 1 / (1 + np.exp(-(h + J * t))))) != 0))

    def falte(innen, aussen):
        for _ in range(60):
            m = (innen + aussen) / 2
            innen, aussen = (m, aussen) if zahl(m) == 3 else (innen, m)
        return (innen + aussen) / 2
    return falte(-J / 2, 0.0) - falte(-J / 2, -J)


def gewebe_gedaechtnis(folge, numerisch=False, n=300, saat=3):
    """Gewebe aus n kooperativen Gruppen über ihrer Kippstelle, je Gruppe J gleichverteilt zwischen 4,5 und 8 und die
    Mitte um c (normal, Streuung 1) versetzt, θ = 1/(1 + e^{−(h + c + J·θ)}). Die Gabe h läuft die Extremwerte folge
    der Reihe nach ab, alle Gruppen anfangs leer. Gibt je Gruppe 1 (besetzt) oder 0. Geschlossen über die Schwellen der
    Falten, besetzt über −J/2 + D/2 − c, leer unter −J/2 − D/2 − c (D aus besetzung_hysterese). numerisch: zweiter Weg,
    die Besetzung in Schritten von 0,01 nachgeführt, je Schritt 400-mal eingesetzt."""
    import numpy as np
    zufall = np.random.default_rng(saat)
    J, c = zufall.uniform(4.5, 8.0, n), zufall.normal(0.0, 1.0, n)
    weg, halt = [folge[0]], {0}
    for a, b in zip(folge, folge[1:]):
        schritte = max(2, int(round(abs(b - a) / 0.01)))
        weg += list(np.linspace(a, b, schritte + 1)[1:])
        halt.add(len(weg) - 1)
    if numerisch:
        th = np.full(n, 1e-6)
        for i, h in enumerate(weg):
            for _ in range(200000 if i in halt else 400):             # an jedem Extremwert einschwingen lassen
                th = 1 / (1 + np.exp(-(h + c + J * th)))
        return (th > 0.5).astype(int)
    D = np.array([besetzung_hysterese(j) for j in J])
    hoch, tief = -J / 2 + D / 2 - c, -J / 2 - D / 2 - c
    zustand = np.zeros(n, dtype=int)
    for h in weg:
        zustand = np.where(h > hoch, 1, np.where(h < tief, 0, zustand))
    return zustand


def besetzung_rueckkehr(h, J, numerisch=False):
    """Rate, mit der die Besetzung nach einer kleinen Störung zurückkehrt, dθ/dt = −θ + 1/(1 + e^{−(h + Jθ)}) (Grundrate
    eins), auf dem unteren Ast: 1 − q mit der Verstärkung der Schleife q = J·θ(1 − θ). Bei halber Besetzung ist das 1/n_H,
    vor einer Falte fällt sie mit der Wurzel des Abstands. numerisch: zweiter Weg, das Abklingen integriert."""
    import numpy as np
    from scipy.integrate import solve_ivp
    t = np.linspace(1e-9, 1 - 1e-9, 200001)
    g = t - 1 / (1 + np.exp(-(h + J * t)))
    i = int(np.argmax(np.diff(np.sign(g)) != 0))
    from scipy.optimize import brentq
    th = brentq(lambda x: x - 1 / (1 + math.exp(-(h + J * x))), t[i], t[i + 1], xtol=1e-16)
    if not numerisch:
        return 1 - J * th * (1 - th)
    r = solve_ivp(lambda _, y: [-y[0] + 1 / (1 + math.exp(-(h + J * y[0])))], (0, 2.0), [th + 1e-7], rtol=1e-12,
                  atol=1e-18, dense_output=True)
    return float(-math.log((r.sol(2.0)[0] - th) / 1e-7) / 2.0)


def rueckkehr_diffusion(s, punkte=4000, weite=40.0):
    """Zweiter Weg zum Rückkehranteil: G″ = s²·G auf [0, weite] (in Einheiten des Rückkehrwegs) mit G(0) = 1 und G = 0 am
    fernen Rand, finite Differenzen. Gibt G(1), geschlossen e^{−s}."""
    import numpy as np
    from scipy.linalg import solve_banded
    h = weite / punkte
    n = punkte - 1
    ab = np.zeros((3, n))
    ab[0, 1:], ab[1, :], ab[2, :-1] = 1.0, -(2 + (s * h)**2), 1.0
    rhs = np.zeros(n)
    rhs[0] = -1.0
    G = solve_banded((1, 1), ab, rhs)
    return float(np.interp(1.0, np.arange(1, punkte) * h, G))


def einrasten_zeit(dw, K):
    """Zwei Zyklen mit der Verstimmung Δω und der dissipativen Kopplung K = 2k·|cos φ|, ψ̇ = Δω − K·sin ψ: eingerastet
    (|Δω| < K) die Zeit, in der sich die Phase auf die Passung zieht, 1/√(K² − Δω²), sonst die Dauer eines
    Phasenschlupfs, 2π/√(Δω² − K²). Am Rand des Fensters wachsen beide ohne Grenze (Kippstelle)."""
    d = K * K - dw * dw
    return 1 / math.sqrt(d) if d > 0 else 2 * math.pi / math.sqrt(-d)


def einrasten_schwelle(dw, K):
    """Höhe der Schwelle, die die eingerastete Phase hält, im gekippten Potential V(ψ) = −Δω·ψ − K·cos ψ:
    ΔV = 2√(K² − Δω²) − 2Δω·arccos(Δω/K), bei Δω = 0 genau 2K, am Rand des Fensters null."""
    return 2 * math.sqrt(K * K - dw * dw) - 2 * dw * math.acos(dw / K)


def einrasten_numerisch(dw, K, art="zeit"):
    """Zweiter Weg zu einrasten_zeit und einrasten_schwelle: ψ̇ = Δω − K·sin ψ numerisch integriert (Dauer eines
    Schlupfs, Abklingen einer kleinen Auslenkung) und V(ψ) auf einem feinen Gitter abgetastet (Höhe der Schwelle)."""
    import numpy as np
    from scipy.integrate import solve_ivp
    if art == "schwelle":
        x = np.linspace(-math.pi, 2 * math.pi, 600001)
        V = -dw * x - K * np.cos(x)
        i = int(np.argmin(np.where(x < math.pi / 2, V, np.inf)))
        j = i + int(np.argmax(V[i:i + 300000]))
        return float(V[j] - V[i])
    f = lambda t, y: [dw - K * math.sin(y[0])]
    if dw * dw > K * K:                                          # Dauer eines Schlupfs: ψ wächst um 2π
        ev = lambda t, y: y[0] - 2 * math.pi
        ev.terminal = True
        r = solve_ivp(f, (0, 1e6), [0.0], events=ev, rtol=1e-11, atol=1e-12)
        return float(r.t_events[0][0])
    p = math.asin(dw / K)                                       # Abklingen um die Passung: e^{−t/τ}
    r = solve_ivp(f, (0, 5.0), [p + 1e-6], rtol=1e-12, atol=1e-16, dense_output=True)
    return float(-5.0 / math.log((r.sol(5.0)[0] - p) / 1e-6))


def spiegel_deckung(punkte, gitter=0):
    """Wie weit sich benannte Punkte (Kontakte) durch eine Drehung auf ihr Spiegelbild legen lassen: kleinster mittlerer
    Abstand nach der besten eigentlichen Drehung (Kabsch, det = +1) zwischen den Punkten und ihrem Spiegelbild x → −x,
    beide um den Schwerpunkt gelegt. Null heißt: Die Spiegelbilder sind nicht zu unterscheiden. gitter > 0: als zweiter
    Weg alle Drehungen auf einem Gitter der Euler-Winkel abgetastet (gitter Stufen je Winkel)."""
    import numpy as np
    P = np.asarray(punkte, dtype=float)
    P = P - P.mean(axis=0)
    S = P * np.array([-1.0, 1.0, 1.0])
    if gitter:
        def dreh(a, b, c):
            z = lambda w: np.array([[math.cos(w), -math.sin(w), 0], [math.sin(w), math.cos(w), 0], [0, 0, 1]])
            y = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
            return z(a) @ y @ z(c)
        g, h = np.linspace(0, 2 * math.pi, 2 * gitter + 1), np.linspace(0, math.pi, gitter + 1)
        return float(min(np.sqrt(np.mean(np.sum((S @ dreh(a, b, c) - P)**2, axis=1))) for a in g for b in h for c in g))
    U, _, Vt = np.linalg.svd(S.T @ P)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return float(np.sqrt(np.mean(np.sum((S @ R - P)**2, axis=1))))


def allosterie(dw, verschiebung, K=1.0, D=0.25, numerisch=False):
    """Allosterischer Wirkstoff als reaktiver Teil der Kopplung: Er rastet nicht ein, er verschiebt die Frequenz der
    Tasche um verschiebung und damit die Verstimmung jedes anderen Wirkstoffs, Δω → Δω + verschiebung. Gibt den Faktor
    der Verweilzeit des anderen, aus der Rate des Lösens (√(K² − Δω²)/2π)·e^{−ΔV/D} mit ΔV aus einrasten_schwelle(|Δω|).
    numerisch: zweiter Weg, die Schwelle aus dem abgetasteten Potential (einrasten_numerisch)."""
    def rate(d):
        d = abs(d)
        dv = einrasten_numerisch(d, K, "schwelle") if numerisch else einrasten_schwelle(d, K)
        return math.sqrt(K * K - d * d) / (2 * math.pi) * math.exp(-dv / D)
    return rate(dw) / rate(dw + verschiebung)


def wechsel_rate(k_an, k_ab, numerisch=False, schritte=400000):
    """Wechsel einer Bindestelle je Zeit: Lösungen (gleich viele Bindungen) bei der Rate des Bindens k_an = k_on·[L] und
    des Lösens k_ab, geschlossen k_ab·θ mit θ = k_an/(k_an + k_ab). numerisch: zweiter Weg, eine gewürfelte Folge von
    Binden und Lösen mit festem Startwert, die Lösungen gezählt und durch die Dauer geteilt."""
    if not numerisch:
        return k_ab * k_an / (k_an + k_ab)
    import numpy as np
    zufall = np.random.default_rng(7)
    leer = zufall.exponential(1 / k_an, schritte)                # Wartezeiten leer, bis gebunden wird
    voll = zufall.exponential(1 / k_ab, schritte)                # Wartezeiten gebunden, bis gelöst wird
    return float(schritte / (leer.sum() + voll.sum()))


def erster_durchgang(x, numerisch=False):
    """Rückkehranteil der Faltung als erster Durchgang: Die Zeit T, in der eine Irrfahrt in einer Dimension (D = 1) die
    Strecke 1 zurücklegt, ist verteilt wie T^{−3/2}·e^{−1/(4T)}/√(4π), und ihr Mittel ⟨e^{−x·T}⟩ ist geschlossen
    e^{−√x}, mit x = a_N/a₀ die Regel. numerisch: zweiter Weg, das Mittel über die Verteilung integriert."""
    if not numerisch:
        return math.exp(-math.sqrt(x))
    from scipy.integrate import quad
    f = lambda u: math.exp(-x * math.exp(u) - 1 / (4 * math.exp(u))) * math.exp(-0.5 * u) / math.sqrt(4 * math.pi)
    return quad(f, -30, 30, limit=400, epsabs=1e-14, epsrel=1e-12)[0]


def irrfahrt_rueckkehr(dim, numerisch=True):
    """Anteil der Irrfahrten auf dem hyperkubischen Gitter der Dimension dim, die je zu ihrem Anfang zurückkehren (Pólya):
    1 − 1/u, u die mittlere Zahl der Besuche am Anfang, u = ∫₀^∞ (e^{−t/dim}·I₀(t/dim))^dim dt, in einer und zwei
    Dimensionen eins. numerisch=False, nur dim = 3: zweiter Weg, Watsons geschlossene Form
    u = √6/(32π³)·Γ(1/24)·Γ(5/24)·Γ(7/24)·Γ(11/24)."""
    if dim <= 2:
        return 1.0
    if not numerisch:
        g = math.gamma
        return 1 - 1 / (math.sqrt(6) / (32 * math.pi**3) * g(1 / 24) * g(5 / 24) * g(7 / 24) * g(11 / 24))
    from scipy.integrate import quad
    from scipy.special import ive
    u = quad(lambda t: ive(0, t / dim)**dim, 0, math.inf, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    return 1 - 1 / u


def kac_wiederkehr(p_ab=0.1, p_ba=0.3, numerisch=False, schritte=1000000, saat=4):
    """Geschlossenes, stationäres System aus zwei Zuständen, je Schritt von A nach B mit p_ab, von B nach A mit p_ba.
    Gibt (Σ 1/⟨T⟩ der mittleren Wiederkehrzeiten, gezählt von jedem Schritt im Zustand aus, Σ 1/⟨τ⟩ der mittleren
    Verweilzeiten). Geschlossen nach Kac 1, die Verweilzeiten 1/p_ab und 1/p_ba geben p_ab + p_ba. numerisch: zweiter
    Weg, die Folge der Zustände gewürfelt und gezählt."""
    if not numerisch:
        return 1.0, p_ab + p_ba
    import numpy as np
    zufall = np.random.default_rng(saat)
    u = zufall.random(schritte)
    z = np.zeros(schritte, dtype=int)
    for i in range(1, schritte):
        z[i] = (1 if u[i] < p_ab else 0) if z[i - 1] == 0 else (0 if u[i] < p_ba else 1)
    summe_T = sum(1 / np.mean(np.diff(np.where(z == k)[0])) for k in (0, 1))
    wechsel = np.flatnonzero(np.diff(z)) + 1
    laeufe = np.diff(np.concatenate([[0], wechsel, [schritte]]))
    art = z[np.concatenate([[0], wechsel])]
    summe_tau = sum(1 / np.mean(laeufe[1:-1][art[1:-1] == k]) for k in (0, 1))
    return float(summe_T), float(summe_tau)


def regelkreis(q, tau=1.0, numerisch=False):
    """Regelkreis aus Teil und Ganzem: Die Wirkung W = Gabe + R, und der Zustand der Zelle R folgt mit der Rate 1/τ der
    Verstärkung der Schleife, τ·dR/dt = q·W − R (q < 0 Gegenkopplung, 0 < q < 1 Mitkopplung). Gibt (Wirkung im
    Gleichgewicht unter der Gabe 1, Wirkung gleich nach dem Absetzen, Rate des Einschwingens und Abklingens), geschlossen
    1/(1 − q), q/(1 − q), (1 − q)/τ. numerisch: zweiter Weg, die Gleichung integriert, Gabe an, dann ab."""
    if not numerisch:
        return 1 / (1 - q), q / (1 - q), (1 - q) / tau
    import numpy as np
    from scipy.integrate import solve_ivp
    dauer = 60 * tau / (1 - q)
    an = solve_ivp(lambda t, R: [(q * (1 + R[0]) - R[0]) / tau], (0, dauer), [0.0], rtol=1e-12, atol=1e-14)
    R0 = an.y[0, -1]
    ab = solve_ivp(lambda t, R: [(q * R[0] - R[0]) / tau], (0, 1.0), [R0], rtol=1e-12, atol=1e-14, dense_output=True)
    rate = -math.log(ab.sol(1.0)[0] / R0)
    return float(1 + R0), float(R0), float(rate)


def wechsel_grundtakt(theta, S0, k_ab, numerisch=False, dauer=200000.0):
    """Wechsel, die die Zelle liest, wenn der freie Rezeptor selbst mit dem Grundtakt S0 wechselt und gebunden mit dem
    Wirkstoff (k_ab): S = (1 − θ)·S0 + θ·k_ab. numerisch: zweiter Weg, eine gewürfelte Folge, frei mit der Rate des
    Bindens k_ab·θ/(1 − θ), gebunden mit k_ab, in den freien Zeiten Wechsel des Rezeptors mit S0, alle Wechsel gezählt."""
    if not numerisch:
        return (1 - theta) * S0 + theta * k_ab
    import numpy as np
    zufall = np.random.default_rng(11)
    k_an = k_ab * theta / (1 - theta)
    t, zahl = 0.0, 0
    while t < dauer:
        frei = zufall.exponential(1 / k_an)
        zahl += zufall.poisson(S0 * frei) + 1                     # eigene Wechsel in der freien Zeit, dann ein Lösen
        t += frei + zufall.exponential(1 / k_ab)
    return zahl / t


def resonanz_fenster(p, q, K, numerisch=False):
    """Fenster der Resonanz p:q beim Einrasten je Umlauf, θ → θ + Ω − (K/2π)·sin 2πθ (eine Gabe je Umlauf): Breite des
    Bereichs von Ω, in dem die Umlaufzahl genau p/q ist. Gesucht über die Falten der q-fachen Abbildung, an denen
    F^q(θ) − θ − p sein Minimum oder Maximum bei null hat, ausgehend von dem Ω, an dem sein Mittel null ist. numerisch:
    zweiter Weg, je Ω 10⁶ Schritte, im Fenster bleibt der aufgerollte Winkel nahe n·p/q, außerhalb läuft er davon, die
    Ränder durch Halbieren."""
    import numpy as np
    from scipy.optimize import brentq

    def F(t, om, n):
        for _ in range(n):
            t = t + om - K / (2 * math.pi) * np.sin(2 * math.pi * t)
        return t
    weit = 0.5 / q**2
    t = np.linspace(0, 1, 4001)
    mitte = brentq(lambda om: float(np.mean(F(t, om, q) - t - p)), p / q - weit, p / q + weit, xtol=1e-15)
    if numerisch:
        def drin(om):
            x, n = 0.0, 1000000
            for k in range(n):
                x = x + om - K / (2 * math.pi) * math.sin(2 * math.pi * x)
                if abs(x - (k + 1) * p / q) > 2:
                    return False
            return True
        rand = []
        for aussen in (mitte - weit, mitte + weit):
            innen = mitte
            for _ in range(40):
                m = (innen + aussen) / 2
                innen, aussen = (m, aussen) if drin(m) else (innen, m)
            rand.append((innen + aussen) / 2)
        return rand[1] - rand[0]
    tief = lambda om: float(np.min(F(t, om, q) - t - p))
    hoch = lambda om: float(np.max(F(t, om, q) - t - p))
    return brentq(tief, mitte, mitte + weit, xtol=1e-15) - brentq(hoch, mitte - weit, mitte, xtol=1e-15)


def bindung_feld(q, k_ab=1.0, k_bind=1.0, K_nah=1.0, numerisch=False):
    """Rückkehr im Feld der Rezeptoren mit drei Zuständen: gebunden (B), nah im Feld (N), fern (F). B → N mit k_ab,
    N → B mit k_bind, N → F mit k_weg, F → N mit K_nah·k_weg (das Gleichgewicht von nah und fern hängt nur an der Form
    des Feldes), der Rückkehranteil q = k_bind/(k_bind + k_weg). Geschlossen: das Verhältnis P_B/P_F = K_nah·k_bind/k_ab
    hängt nicht an q, die mittleren Zeiten B → F und F → B sind (1/k_ab + 1/(k_bind + k_weg))/(1 − q) und
    (1/(K_nah·k_weg) + 1/(k_bind + k_weg))/q. numerisch: zweiter Weg, stationäre Verteilung und mittlere Übergangszeiten
    aus der Übergangsmatrix. Gibt (P_B/P_F, Zeit B → F, Zeit F → B)."""
    import numpy as np
    k_weg = k_bind * (1 - q) / q
    k_an = K_nah * k_weg
    if not numerisch:
        n = 1 / (k_bind + k_weg)
        return K_nah * k_bind / k_ab, (1 / k_ab + n) / (1 - q), (1 / k_an + n) / q
    Q = np.zeros((3, 3))                                        # 0 gebunden, 1 nah, 2 fern
    Q[0, 1], Q[1, 0], Q[1, 2], Q[2, 1] = k_ab, k_bind, k_weg, k_an
    Q -= np.diag(Q.sum(axis=1))
    pi = np.linalg.lstsq(np.vstack([Q.T, np.ones(3)]), np.array([0, 0, 0, 1.0]), rcond=None)[0]

    def zeit(start, ziel):                                     # mittlere Zeit bis zum ersten Erreichen des Ziels
        rest = [i for i in range(3) if i != ziel]
        T = np.linalg.solve(Q[np.ix_(rest, rest)], -np.ones(len(rest)))
        return float(T[rest.index(start)])
    return float(pi[0] / pi[2]), zeit(0, 2), zeit(2, 0)


def belegung_haushalt(an, ab):
    """Zeitanteile eines Rezeptors, um den mehrere Wirkstoffe werben: Zustand 0 leer, Zustand i mit Wirkstoff i, Übergänge
    0 → i mit der Rate an[i] (k_on·[L_i]) und i → 0 mit ab[i] (k_off). Stationäre Verteilung der Übergangsmatrix, die
    Anteile summieren sich zu eins (Haushalt der Wiederkehr). Gibt [π₀, π₁, …]."""
    import numpy as np
    n = len(an) + 1
    Q = np.zeros((n, n))
    for i, (r1, r2) in enumerate(zip(an, ab), start=1):
        Q[0, i], Q[i, 0] = r1, r2
    Q -= np.diag(Q.sum(axis=1))
    A = np.vstack([Q.T, np.ones(n)])
    return np.linalg.lstsq(A, np.concatenate([np.zeros(n), [1.0]]), rcond=None)[0]


def hemmer_haelfte(L_durch_Kd):
    """Konzentration eines Hemmers [I]/K_I, bei der die Belegung durch den Wirkstoff auf die Hälfte fällt, aus der
    stationären Verteilung (belegung_haushalt). Geschlossen 1 + [L]/K_d (Cheng-Prusoff)."""
    from scipy.optimize import brentq
    th = lambda x: belegung_haushalt([L_durch_Kd, x], [1.0, 1.0])[1]
    t0 = th(0.0)
    return brentq(lambda x: th(x) - t0 / 2, 1e-9, 1e6, xtol=1e-14)


def zyklen_im_band(band_durch_k):
    """Höchstzahl freier Zyklen, die ein Medium mit Paarkopplung k im Frequenzband W trägt: N ≈ 1 + 2W/(πk)."""
    return 1 + 2 * band_durch_k / math.pi


def bindung_huelle(N, a, c, P=math.inf, b=2 / 3):
    """Wirkstoff als gehaltener Trichter: N schwere Atome auf Abstand des Halters bilden die Hülle R = λ_H·N^{1/3}, gebunden
    wird über die Hülle, ∝ N^{2/3}. Obere Grenze der gelesenen Bindung pK = a·N^b − c (c der Eintritt: die Richtungen des
    ganzen Körpers werden inaktiv), am Treffpunkt der Hülle mit der Tasche die Höhe P."""
    import numpy as np
    return np.minimum(a * np.asarray(N, dtype=float)**b - c, P)


def bindung_je_atom_groesse(a, c, b=2 / 3, numerisch=False):
    """Größe mit der stärksten Bindung je Atom, d/dN [(a·N^b − c)/N] = 0: N_e = (c/(a·(1 − b)))^{1/b}, bei b = 2/3
    (3c/a)^{3/2}: Dort ist die Bindung über die Hülle dreimal der Eintritt. numerisch: zweiter Weg, das Maximum gesucht."""
    if numerisch:
        from scipy.optimize import minimize_scalar
        return float(minimize_scalar(lambda N: -(a * N**b - c) / N, bounds=(1e-3, 1e6), method="bounded",
                                     options={"xatol": 1e-10}).x)
    return (c / (a * (1 - b)))**(1 / b)


def rueckkehr_windung_reihe(n=2):
    """Wiederholte Rückkehr der Windung als geometrische Reihe: Σ_{k≤n} (α/2π)^k; zum Vergleich mit der Messung."""
    return sum((alpha / (2 * math.pi))**k for k in range(n + 1))


mu0 = 4 * math.pi * 1e-7 * 1.000_000_000_55   # Vs/(Am), CODATA 2018
J_ERDE = 5.86e33                              # kg m²/s, Drehimpuls der Erde (Eingabe, Trägheitsmoment 0,3307·M·R²)
R_GPB = 7.027e6                               # m, Bahnradius von Gravity Probe B (Eingabe)


def verwindung(J, M):
    """a_J = J/(Mc): die Verwindung eines Trichters als Länge."""
    return J / (M * c)


def horizonte(J, M, Q=0.0):
    """Die beiden Horizonte des verwundenen Trichters als Wurzeln von r² − 6Δs·r + a_J² + r_Q·3Δs = 0, r_Q = Q²/(4πε₀Mc²):
    Summe die Tiefe, Produkt Verwindung und Windung, die Höhe im Thaleskreis die Hypotenuse aus a_J und √(r_Q·3Δs),
    Treffpunkt bei j² + r_Q/(3Δs) = 1. Ohne Windung für j = a_J/(3Δs) > 1 gespiegelt komplex, r± = 3Δs ± i·√(a_J² − 9Δs²),
    auf dem Kreis |r±| = a_J im Winkel cos φ = 1/j. Gibt (r₊, r₋, j), j = (J/ħ)/α_G."""
    import cmath
    d, a = 3 * faltung(M), verwindung(J, M)
    w = cmath.sqrt(d * d - a * a - Q**2 / (4 * math.pi * eps0 * M * c**2) * d)
    return d + w, d - w, a / d


def falte_gedreht(M, J, Q=0.0):
    """Die gedrehte Falte für j² + r_Q/(3Δs) ≤ 1. Ohne Windung ein rechtwinkliges Dreieck im Thaleskreis über r₀ = 6Δs: Höhe
    a_J, der Fußpunkt teilt r₀ in r₊ und r₋ (Höhensatz, a_J² = r₊r₋), Katheten r₀(M_irr)² = r₀·r₊ und r₀(M_rot)² = r₀·r₋
    (Kathetensatz), also M² = M_irr² + M_rot² und M_irr·M_rot = Jc/(2G). Mit Windung M² = (M_irr + M_el)² + M_rot²,
    M_el = Q²/(16πε₀G·M_irr), ohne Drehung M_irr die Faltung innerhalb r₊. M_irr aus der Fläche 4π·(r₊² + a_J²) = 4π·r₀(M_irr)²,
    Entropie S/k = 4π·(M_irr/m_Pl)², herausholbar (M − M_irr)/M, Temperatur ħc·(r₊ − r₋)/(4πk·r₀(M_irr)²). Gibt dict."""
    rp, rm, j = horizonte(J, M, Q)
    rQ = Q**2 / (4 * math.pi * eps0 * M * c**2)
    if j * j + rQ / (3 * faltung(M)) > 1 + 2e-9:
        raise ValueError("j² + r_Q/(3Δs) > 1: kein Dreieck, die Horizonte sind komplex")
    a = verwindung(J, M)
    rp, rm = rp.real, rm.real
    A = rp * rp + a * a
    M_irr = c**2 * math.sqrt(A) / (2 * G)
    return dict(j=j, r_Q=rQ, M_irr=M_irr, M_rot=J * c / (2 * G * M_irr), M_el=Q**2 / (16 * math.pi * eps0 * G * M_irr),
                flaeche=4 * math.pi * A, entropie=4 * math.pi * (M_irr / m_Pl)**2, herausholbar=1 - M_irr / M,
                temperatur=hbar * c * (rp - rm) / (4 * math.pi * k_B * A))


def falte_endzustand(M, a0_=a0):
    """Die Falte im Endzustand: Ohne flachen Außenraum geht die ruhende Uhr mit dτ/dt = √f, f = 1 − 6Δs/r − (H_∞·r/c)²,
    H_∞ = c·√Λ_s. Am schnellsten bei R_max³ = GM/H_∞², wo ein Körper frei ruht, mit dem Takt der Bahnuhr
    √(1 − 9Δs/R_max) = √(1 − (M/M_N)^{2/3}). Die Horizonte sind die Wurzeln von u³ − u + u₀ = 0, u = r·H_∞/c,
    u₀ = 6Δs·H_∞/c: u_k = (2/√3)·cos(θ − 2πk/3), cos 3θ = −M/M_N, M_N = c³/(3√3·G·H_∞). Die Wurzel der Falte folgt aus dem
    Produkt der Wurzeln, −u₀, ohne Auslöschung. Wärme für die frei ruhende Uhr bei R_max, kT = ħc·|∂_r f|/(4π·√f(R_max)),
    bei M_N für beide ħc/(2π·r_N). M ist die Masse, deren Feld die Quelle trägt (beim Haufen (1 + Q)·M_b). Über M_N bleiben
    Takt, Horizonte und Wärme None. Gibt dict(R_max, takt_max, r_falte, r_ganzes, r_negativ, T_falte, T_ganzes, M_N, r_N)."""
    H = groesste_umkehr(M, a0_)[1]
    L = c / H
    M_N = c**3 / (3 * math.sqrt(3) * G * H)
    r_N = L / math.sqrt(3)
    R = (G * M / H**2)**(1 / 3)
    r0 = 6 * faltung(M)
    erg = dict(R_max=R, takt_max=None, r_falte=None, r_ganzes=None, r_negativ=None, T_falte=None, T_ganzes=None,
               M_N=M_N, r_N=r_N)
    if M > M_N:
        return erg
    th = math.acos(max(-1.0, -M / M_N)) / 3
    u_g = 2 / math.sqrt(3) * math.cos(th)
    u_n = 2 / math.sqrt(3) * math.cos(th - 4 * math.pi / 3)
    u_f = r0 / L / (u_g * -u_n)
    f_max = -math.expm1(2 / 3 * math.log(M / M_N)) if M > 0 else 1.0
    if f_max > 1e-12:
        T = lambda r: hbar * c * abs(r0 / r**2 - 2 * r / L**2) / (4 * math.pi * k_B * math.sqrt(f_max))
        T_f, T_g = (T(u_f * L) if u_f > 0 else None), T(u_g * L)
    else:
        T_f = T_g = hbar * c / (2 * math.pi * k_B * r_N)
    erg.update(takt_max=math.sqrt(f_max), r_falte=u_f * L, r_ganzes=u_g * L, r_negativ=u_n * L, T_falte=T_f, T_ganzes=T_g)
    return erg

def mitfuehrung(J, r, bahnebene=False):
    """Mitführung durch die Verwindung: die Achse eines Kreisels auf polarer Bahn im Mittel (Schiff), Ω = GJ/(2c²r³) = 3Δs·c·a_J/(2r³),
    mit bahnebene=True die Ebene der Bahn selbst (Lense-Thirring), Ω = 2GJ/(c²r³) = 6Δs·c·a_J/r³."""
    return G * J * (2 if bahnebene else 0.5) / (c**2 * r**3)


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


def takte_kreisbahn(delta):
    """Kreisbahn um einen Körper ohne Verwindung bei δ = 3Δs/r, jeder Kreis die Nullstelle eines Takts: ruhende Uhr √(1 − 2δ)
    (6Δs), umlaufende √(1 − 2δ − v²/c²) = √(1 − 3δ) mit v² = δ·c² von außen (9Δs, dort v_ort = c), radiale Schwingung gegen
    den Umlauf (Re τ_M)² = 1 − 6δ (18Δs, v_ort = c/2). Dazu v_ort/c = √(δ/(1 − 2δ)), Energie (1 − 2δ)/√(1 − 3δ) (eins bei
    12Δs) und Kreisel je Umlauf 2π·(1 − √(1 − 3δ)). Gibt dict."""
    um = math.sqrt(1 - 3 * delta) if delta <= 1 / 3 else float("nan")
    return dict(ruhend=math.sqrt(max(1 - 2 * delta, 0.0)), umlaufend=um, radial=1 - 6 * delta,
                tempo_ort=math.sqrt(delta / (1 - 2 * delta)) if delta < 0.5 else float("nan"),
                energie=(1 - 2 * delta) / um if delta < 1 / 3 else float("nan"),
                kreisel=2 * math.pi * (1 - um) if delta <= 1 / 3 else float("nan"))


def wachsen_grenze(eps, th=None, a0_=a0):
    """Wachsen einer Falte an der Eddington-Grenze: Gibt, was zufällt, bis zur innersten stabilen Bahn den Anteil ε seiner
    Ruheenergie ab und erreicht diese Abgabe L_Edd, wächst die Falte mit dM/dt = ((1 − ε)/ε)·L_Edd/c², um den Faktor e in
    t_e = (ε/(1 − ε))·t_Edd. Die Zeit der Grenze ist die Tiefe des Protons an r_e gespiegelt, c·t_Edd = σ_T/(4π·3Δs_p) =
    ⅔·r_e²/(3Δs_p) mit σ_T = 4π·⅔·r_e². An r*_p gespiegelt gibt dieselbe Tiefe ℓ₀, also t_Edd = ⅔·(r_e/r*_p)²·τ_D für jede
    Masse, a₀ fällt heraus. ε/(1 − ε) = 1/E − 1 mit der Energie E der innersten Bahn je Ruheenergie (ungedreht
    3/(2√2) − 1, bei j = 1 √3 − 1). Mit th die Zahl der e-Stufen bis zur Uhr Θ, ln(M/M_K) = t(Θ)/t_e =
    ((1 − ε)/ε)·(r*_p/r_e)²·Θ_∞·artanh(Θ/Θ_∞)/2π. Gibt dict(t_edd [s], t_e [s], stufen)."""
    r_e = alpha * compton(m_e)
    t_edd = 2 / 3 * (r_e / schwellenradius(m_p, a0_))**2 * rueckkehrzeit(a0_)
    t_e = eps / (1 - eps) * t_edd
    return dict(t_edd=t_edd, t_e=t_e, stufen=None if th is None else uhr_theta(th, a0_)[0] / t_e)


def gezeiten_falte(m, R, M=None, kreisbahn=False):
    """Stern der Masse m [kg] mit der Hülle R [m] an einer ungedrehten Falte der Masse M [kg]. Er bleibt ganz, solange R
    unter seiner Glättungslänge ℓ³ = Gm/𝒯 liegt. Die Gezeit zählt im Rahmen des Sterns: Quer zum Radius bewegt ist ihr
    streckender Eigenwert 𝒯 = (GM/r³)·(2 + 3L²/(r²c²)) = (3γ² − 1)·GM/r³, γ der Gang der ruhenden Uhr gegen seine. An der
    Umkehr einer Bahn mit Energie eins ist γ² = 1/(1 − 2δ), δ = 3Δs/r, die kleinste Umkehr 12Δs (δ = ¼). Auf der Kreisbahn
    (kreisbahn) kommt der Umlauf hinzu, der mitlaufende Kreisel dreht in der eigenen Zeit mit √(GM/r³) gegen den Radius,
    mit γ² = (1 − 2δ)/(1 − 3δ) zusammen 3γ²·GM/r³, Grenze die innerste stabile Bahn 18Δs (δ = ⅙). r_T [m] für M (Vorgabe
    None): Umkehr oder Kreis, an der ℓ = R wird, weit außen R·(2M/m)^{1/3} oder R·(3M/m)^{1/3}. M_V [kg]: Masse, bei der
    das am Kreis geschieht, M_V/m = √f·(δ_K/δ_m)^{3/2} mit f = 𝒯r³/GM dort (5 oder 4) und δ_m = 3Δs_m/R. Gibt
    dict(r_T, M_V, faktor f, spreizung δ_m)."""
    from scipy.optimize import brentq
    f = (lambda u: 3 * (1 - 2 * u) / (1 - 3 * u)) if kreisbahn else (lambda u: 3 / (1 - 2 * u) - 1)
    dK, u_max = (1 / 6, 1 / 3) if kreisbahn else (1 / 4, 1 / 2)
    d = 3 * faltung(m) / R
    r_T = None
    if M is not None:
        ln_x = math.log(3 * faltung(M) / R)
        u = brentq(lambda u: 3 * (ln_x - math.log(u)) - math.log(f(u) * M / m), 1e-300, u_max * (1 - 1e-12))
        r_T = 3 * faltung(M) / u
    return dict(r_T=r_T, M_V=m * math.sqrt(f(dK)) * (dK / d)**1.5, faktor=f(dK), spreizung=d)


def umlauf_abbildung(k2, lauf=2 * math.pi):
    """Abbildung einer kleinen Abweichung nach einem Lauf, x″ + k²·x = 0 im eigenen Lauf s (08): (x, x′) linear abgebildet,
    Determinante 1, Spur 2·cos(k·s), die Eigenwerte ein Spiegel mit dem festen Produkt eins, λ² − Spur·λ + 1 = 0. Nach der
    Spur von den drei Arten der Symmetrie (01): k² > 0 dreht (Eigenwerte e^{±ik·s} auf dem Einheitskreis, keine feste
    Richtung, außer bei ganzem k·s/π), k² = 0 schert, ((1, s), (0, 1)) mit x′ = 0 als einziger fester Richtung, k² < 0 dehnt
    (e^{±√(−k²)·s}, zwei feste Richtungen). Bahn: k² = (Re τ_M)² = κ²/Ω² im Bahnwinkel, s = 2π, mit der Faltung 1 − 6δ
    (Scherung bei 18Δs, −1 bei 24Δs, Spur 2·cosh 2π an der Photonensphäre), mit der Regel 1 + χ, mit der Rückstellung
    ((4 − d_eff)·g/R − 4H_∞²)/(g/R − H_∞²). Gibt dict(matrix, spur, eigen, art)."""
    if k2 > 0:
        k = math.sqrt(k2)
        co, si = math.cos(k * lauf), math.sin(k * lauf)
        m = ((co, si / k), (-k * si, co))
    elif k2 < 0:
        k = math.sqrt(-k2)
        co, si = math.cosh(k * lauf), math.sinh(k * lauf)
        m = ((co, si / k), (k * si, co))
    else:
        m = ((1.0, lauf), (0.0, 1.0))
    spur = m[0][0] + m[1][1]
    wurzel = complex(spur * spur / 4 - 1) ** 0.5
    art = "Drehung" if k2 > 0 else ("Scherung" if k2 == 0 else "Dehnung")
    return dict(matrix=m, spur=spur, eigen=(spur / 2 + wurzel, spur / 2 - wurzel), art=art)


def flussdimension_bahn(teile):
    """Kreisbahn im Bahnwinkel (08): d²u/dφ² + u = Σ F_i(u), L der Drehimpuls je Masse in der eigenen Zeit. Jedes Glied
    F_i ∝ u^{d_i − 3} trägt auf dem Kreis den Anteil w_i = F_i/u, Σ w_i = 1, eine Abweichung schwingt mit
    (Re τ_M)² = 4 − d_eff, d_eff = Σ w_i·d_i. Für eine Neigung g_i ist F_i = g_i/(L²u²) und d_i = 1 − d ln g_i/d ln r
    (Newton 3, Regel 3 − χ, Rückstellung 0 mit dem negativen Anteil −H_∞²/Ω²). Die Faltung gibt das Glied 3GMu²/c² mit
    d = 5, keine Neigung am Ort: Anteil 3δ, Newton 1 − 3δ (Quadrat der umlaufenden Uhr), also d_eff = 3 + 6δ und
    (Re τ_M)² = 1 − 6δ. Unter 4 dreht die Abbildung, bei 4 schert sie (18Δs, je die Hälfte), darüber dehnt sie
    (Photonensphäre 5). Spiegel der Kräfte: Re τ_M·Re τ_M′ = 2. teile: Folge von (d_i, w_i) oder (d_i, F_i) auf dem Kreis,
    die Anteile werden normiert. Gibt (d_eff, (Re τ_M)²)."""
    s = sum(x for _, x in teile)
    d = sum(di * x for di, x in teile) / s
    return d, 4 - d


def bahn_windung(a_phi, e=0.3, numerisch=False):
    """Bahn um eine Windung (08): Coulomb K_el = |Q₁Q₂|/(4πε₀) mit der Energie des Trägers bei jedem Tempo,
    W = √(p²c² + m²c⁴) − K_el/r. Im Bahnwinkel ist d²u/dφ² + u = K_el·(W + K_el·u)/(p_φ·c)² = K_el·W/(p_φ·c)² + α_φ²·u
    mit α_φ = K_el/(p_φ·c), auf dem Kreis genau sein Tempo gegen c. Coulomb trägt das Quadrat der umlaufenden Uhr,
    1 − α_φ², das zweite Glied den Rest α_φ² mit d = 4 (flussdimension_bahn), keine Neigung, sondern die nach innen
    wachsende Energie des Trägers. Es ist linear in u: Jede Bahn dreht als Rosette um denselben Winkel,
    (Re τ_M)² = 1 − α_φ², schwach je Umlauf πα_φ², bei gleichem Tempo (α_φ² = δ) ein Sechstel der Faltung (6πδ), denn
    die umlaufende Uhr trägt hier keinen Takt. Scherung bei α_φ = 1, K_el = p_φ·c, mit p_φ = ħ bei |Q₁Q₂| = q_P²
    (Zα = 1), darüber Dehnung. numerisch: die Bewegung in der Zeit (m = c = p_φ = 1, K_el = α_φ < 1) vom Perizentrum
    r_k/(1 + e), r_k = √(1/α_φ² − 1) der Kreis, gemessen der Winkel Δφ zwischen zwei Perizentren, (Re τ_M)² = (2π/Δφ)²,
    ohne die Gleichung im Bahnwinkel. Gibt (d_eff, (Re τ_M)², Drehung je Umlauf [rad])."""
    if not numerisch:
        d, k2 = flussdimension_bahn([(3, 1 - a_phi**2), (4, a_phi**2)])
    else:
        from scipy.integrate import solve_ivp
        r_k = math.sqrt(1 / a_phi**2 - 1)

        def f(t, y):
            r, ph, pr = y
            E = math.sqrt(pr * pr + 1 / r**2 + 1)
            return [pr / E, 1 / (E * r * r), 1 / (E * r**3) - a_phi / r**2]
        peri = lambda t, y: y[2]
        peri.direction = 1
        dauer, ph = 8 * math.pi * r_k / (a_phi * (1 - e)**1.5), []
        while len(ph) < 2:
            s = solve_ivp(f, [0, dauer], [r_k / (1 + e), 0.0, 0.0], events=peri, rtol=1e-12, atol=1e-14 * r_k)
            ph, dauer = [y[1] for y in s.y_events[0]], 2 * dauer
        k2 = (2 * math.pi / (ph[1] - ph[0]))**2
        d = 4 - k2
    return d, k2, 2 * math.pi * (1 / math.sqrt(k2) - 1) if k2 > 0 else float("inf")


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


def inaktives_mass(M, H_=H0, a0_=a0, kappa=1.0):
    """Inaktives Maß der Bindung je Baryon im Tiefbereich, (R_ta/r*)^{3/2}, mit dem Rand der Umkehr R_ta³ = GM/(κ·H²).
    In den Größen des Modells ist es a₀/(√κ·H·v_flach) = (v_flach/r*)/(√κ·H): die Umlaufrate an der Schwelle gegen die Rate des
    Atmens. Gibt (R_ta/r*)^{3/2}."""
    R_ta = (G * M / (kappa * H_**2))**(1 / 3)
    return (R_ta / schwellenradius(M, a0_))**1.5


def naht_radius(M, Q=Q_kosm, a0_=a0):
    """R_Q = r*/(−ln q_kosm): dort erreicht das örtliche q das kosmische."""
    return schwellenradius(M, a0_) / (-math.log(Q / (1 + Q)))


def reibung_haufen(M_b, r_J, Q=Q_kosm, a0_=a0):
    """Reibung einer Galaxie, deren Mode übergeben ist, in einem Haufen, dessen Mode es nicht ist. Die Regel des Haufens
    ist die Summe der Rückläufe am Nome des Ganzen, q = Q/(1 + Q), unabhängig von a_N und darum in jeder Richtung gleich:
    Entstehen und Ziehen des Kielwassers zählen Σ(n + 1)·qⁿ = (1 + Q)² mal Chandrasekhar mit den Baryonen. Jenseits ihres
    Jacobi-Radius r_J stört die Galaxie mit ihren Baryonen M_b. Der Beobachter mit Halo rechnet mit Staub (1 + Q)·ρ_b und
    der Masse ν(r_J)·M_b, die er innerhalb r_J aus ihrer Bewegung liest. Gibt ((1 + Q)², Modell durch Beobachter
    (1 + Q)/ν(r_J))."""
    return (1 + Q)**2, (1 + Q) / verstaerkung(G * M_b / r_J**2, a0_)


def welle_haufen(f_gas, alpha_M, gamma_P, gamma=5 / 3):
    """Welle im nie abgekühlten Gas eines Haufens: Die Verdichtung zieht mit (1 + Q)·G nach, ω² = c_s²k²·(1 − q),
    q = 4πG(1 + Q)ρ_Gas/(c_s²k²). Im Gleichgewicht ist 4πG(1 + Q)ρ_b·r² = v_c²·α_M und v_c² = γ_P·kT/(μm_p), darum fällt
    die Regel heraus: Für k = 2π/r ist q = f_Gas·α_M·γ_P/(4π²·γ). Gibt (q, Verstärkung 1/(1 − q), γ_P, bei dem eine
    Welle aus reinem Gas mit α_M = 3 aufschaukelt)."""
    q = f_gas * alpha_M * gamma_P / (4 * math.pi**2 * gamma)
    return q, (1 / (1 - q) if q < 1 else math.inf), 4 * math.pi**2 * gamma / 3


def uebergabe_linse(M200, H0_=H0, Q=Q_kosm, a0_=a0):
    """Bereich, der als Staub M₂₀₀ = (1 + Q)·M_b trägt, R₂₀₀ bei 200·ρ_krit. Ist seine Mode übergeben und liegen seine
    Baryonen innen, liest die Linse in R₂₀₀ die Regel, ν(a_N(R₂₀₀))·M_b statt M₂₀₀. Gibt (R₂₀₀, R_Q, ν/(1 + Q))."""
    Mb = M200 / (1 + Q)
    R = (3 * M200 / (4 * math.pi * 200 * 3 * H0_**2 / (8 * math.pi * G)))**(1 / 3)
    return R, naht_radius(Mb, Q, a0_), verstaerkung(G * Mb / R**2, a0_) / (1 + Q)


def kollaps_kugel(delta_i=1e-5):
    """Kugeliger Kollaps ohne Rückstellung, numerisch: Schale mit wachsender Störung δ_i im flachen Staub (GM = 1,
    r̄³ = 9t²/2), gerechnet bis zur Kreuzung in der Mitte; Rückgabe (linearer Wert δ_i·(t_c/t_i)^{2/3}, Verdichtung
    auf r_U/2 gegen das Mittel bei t_c, Verdichtung an der Umkehr gegen das Mittel bei t_U)."""
    from scipy.integrate import solve_ivp
    rbar = lambda t: (4.5 * t * t)**(1 / 3)
    r0 = rbar(1.0) * (1 - delta_i / 3)
    v0 = 2 / 3 * r0 * (1 - delta_i / 3)
    umkehr = lambda t, y: y[1]
    ende = lambda t, y: y[0] - 1e-3 * r0
    ende.terminal, ende.direction = True, -1
    s = solve_ivp(lambda t, y: [y[1], -1 / y[0]**2], (1.0, 1e14), [r0, v0], events=(umkehr, ende), rtol=1e-12, atol=1e-16,
                  method="DOP853")
    tc = s.t_events[1][0]
    tU, rU = s.t_events[0][0], s.y_events[0][0][0]
    dichte = lambda r, t: 1 / (4 * math.pi / 3 * r**3) * 6 * math.pi * t**2
    return delta_i * tc**(2 / 3), dichte(rU / 2, tc), dichte(rU, tU)


def leerraum_rand(th=None, eps=1e-4):
    """Rand eines Leerraums mit scharfem Rand (außen das Mittel). Er kehrt nie um und atmet mit, die Mode wird wie Staub mit
    1 + Q gelesen: der rollende Kreis mit imaginärem Winkel, r = A(cosh ϑ − 1), t = B(sinh ϑ − ϑ),
    1 + δ = (9/2)(sinh ϑ − ϑ)²/(cosh ϑ − 1)³, linear δ_L = −(3/20)(6(sinh ϑ − ϑ))^{2/3}. Außerhalb des Randes fällt die
    mittlere Störung wie 1/M, die Schalen kreuzen bei d ln(1 + δ)/d ln|δ_L| = −1, also H_v = (1 + f/3)·H mit
    f = d ln D/d ln a. th None: ohne Rückstellung (f = 1) genau, 8(cosh ϑ − 1)² = 9 sinh ϑ (sinh ϑ − ϑ), H_v = (4/3)·H.
    th: Θ, bei dem der Rand kreuzt, im Atmen des Formteils (H_∞ = 1, a³ = sinh²(3t/2), Θ = Θ_∞·tanh(3t/2)), zwei Schalen
    am Rand mit wachsender Mode vom Anfang, die Masse der äußeren um eps größer. Gibt (δ_v linear, 1 + δ, R/R_L, H_v/H, f)."""
    from scipy.optimize import brentq
    if th is None:
        S = lambda x: math.sinh(x) - x
        C = lambda x: math.cosh(x) - 1
        x = brentq(lambda x: 8 * C(x)**2 - 9 * math.sinh(x) * S(x), 1.0, 10.0)
        d = 4.5 * S(x)**2 / C(x)**3
        return -0.15 * (6 * S(x))**(2 / 3), d, d**(-1 / 3), 1.5 * math.sinh(x) * S(x) / C(x)**2, 1.0
    from scipy.integrate import quad, solve_ivp
    t_k = 2 / 3 * math.atanh(th * math.pi / math.sqrt(15))
    a = lambda t: math.sinh(1.5 * t)**(2 / 3)
    E = lambda x: math.sqrt(1 + x**-3)
    D = lambda x: 2.5 * E(x) * quad(lambda u: (u * E(u))**-3, 0, x)[0]
    ti = 1e-5

    def kreuzung(Di):
        start = lambda M, d: (lambda r: [r, r / math.tanh(1.5 * ti) * (1 - d / 3)])(a(ti) * (M / (1 + d))**(1 / 3))
        ev = lambda t, y: y[2] - y[0]
        ev.terminal, ev.direction = True, -1
        f = lambda t, y: [y[1], -0.5 / y[0]**2 + y[0], y[3], -0.5 * (1 + eps) / y[2]**2 + y[2]]
        s = solve_ivp(f, (ti, 4 * t_k), start(1.0, Di) + start(1 + eps, Di / (1 + eps)), events=ev, rtol=1e-11,
                      atol=1e-16, method="DOP853")
        return (s.t_events[0][0], s.y_events[0][0]) if len(s.t_events[0]) else (4 * t_k, None)
    Di = brentq(lambda d: kreuzung(d)[0] - t_k, -2e-2, -1e-4, xtol=1e-12)
    tc, y = kreuzung(Di)
    d = a(tc)**3 / y[0]**3
    f = (math.log(D(a(tc) * 1.0001)) - math.log(D(a(tc) * 0.9999))) / 2e-4
    return Di * D(a(tc)) / D(a(ti)), d, d**(-1 / 3), y[1] / y[0] * math.tanh(1.5 * tc), f


def leerraum_noetig(punkt, H_innen, z_rand):
    """Leerraum um den Beobachter, der innen mit H_innen atmet, wo das Ganze im Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀,
    w_a) mit H₀ atmet, bis zum Rand bei z_rand. Er atmet mit, seine Mode wird wie Staub gelesen: linear ist
    3·(H − H_v) = dδ/dt = f·H·δ_L, also H_v = (1 − f·δ_L/3)·H mit f = d ln D/d ln a heute. Nichtlinear dehnt er sich bei
    gleicher Tiefe weniger, δ_L = −3·(H_innen/H₀ − 1)/f ist die kleinste Tiefe. Sein Radius heute ist der Weg des Lichts
    bis z_rand in seiner eigenen Dehnung, R = c·∫dz/H geteilt durch H_innen/H₀. Seine Masse nähme im Mittel den Radius
    R_L = R·(1 + δ)^{1/3} ein, 1 + δ aus dem rollenden Kreis mit imaginärem Winkel bei δ_L (mit der Rückstellung dieselben
    δ und δ_L, leerraum_rand). Gibt (f, δ_L, R [Mpc], R_L [Mpc], Schwankung σ des linearen Spektrums in einer Kugel R_L,
    |δ_L|/σ)."""
    import numpy as np
    from scipy.integrate import quad
    from scipy.optimize import brentq
    f = float(wachstum_linear(punkt)[1](1.0))
    e = H_innen / punkt[2]
    d = -3 * (e - 1) / f
    H = atmen_punkt(punkt)[4]
    R = c * quad(lambda z: 1 / H(1 / (1 + z)), 0, z_rand)[0] / Mpc / e
    S, C = (lambda x: math.sinh(x) - x), (lambda x: math.cosh(x) - 1)
    x = brentq(lambda x: -0.15 * (6 * S(x))**(2 / 3) - d, 1e-6, 50.0)
    RL = R * (4.5 * S(x)**2 / C(x)**3)**(1 / 3)
    k, P = spektrum_linear(punkt, kmax=10.0)
    y = k * RL
    W = 3 * (np.sin(y) - y * np.cos(y)) / y**3
    s = math.sqrt(np.trapezoid(k**3 * P * W**2 / (2 * math.pi**2), np.log(k)))
    return f, d, R, RL, s, abs(d) / s


def netz_achsen(zaehlung="masse", numerisch=False):
    """Netz aus der Verformung der mitatmenden Mode, die linear wie Staub mit 1 + Q wirkt, x = q − D·∇ψ (Seldowitsch).
    Die Eigenwerte λ₁ ≥ λ₂ ≥ λ₃ von T_ij = ∂_i∂_jψ sind drei Achsen: λ > 0 fällt bei D·λ = 1, λ < 0 dehnt. Gaußsch und
    isotrop ist ⟨T_ij·T_kl⟩ = (σ²/15)·(δ_ij·δ_kl + δ_ik·δ_jl + δ_il·δ_jk), Spur und Scherung sind unabhängig, und von
    ⟨Σλ²⟩ = σ² trägt die Spur 1/3 wie in jeder ebenen Welle T_ij = δ·k̂_i·k̂_j. Nach der Masse fallen alle drei Achsen beim
    Anteil P₃ = 1/2 − arctan(√5)/(2π) − √5/(3π), ebenso viel fällt keine (Spiegel Φ → −Φ, λ_k → −λ_{4−k}), je 1/2 − P₃
    fallen zwei oder eine, ⟨λ₁⟩ = −⟨λ₃⟩ = 3σ/√(10π), ⟨λ₂⟩ = 0. zaehlung "stellen": die Stellen ohne Neigung eines
    isotropen Gaußschen Feldes, deren Krümmung unabhängig von der Neigung dieselbe Verteilung hat, jede gezählt mit
    |λ₁λ₂λ₃| (Kac und Rice): Trichter und Leerraum je (29 − 6√6)/116, Faden und Wandsattel je (29 + 6√6)/116, je Trichter
    (29 + 6√6)/(29 − 6√6) Fadensättel, jeder mit zwei Enden in Trichtern. numerisch: dieselben Anteile aus der Verteilung
    nach Doroschkewitsch (σ = 1), p = (15³/(8π√5))·exp(−3I₁² + (15/2)·I₂)·(λ₁ − λ₂)(λ₂ − λ₃)(λ₁ − λ₃),
    I₁ = Σλ, I₂ = λ₁λ₂ + λ₂λ₃ + λ₁λ₃, jeder Bereich der Vorzeichen (3, 2, 1, 0 positive) einzeln integriert, der Spiegel
    nicht vorausgesetzt. Gibt dict (knoten, faden, wand, leer, bei "masse" dazu l1 [σ] und spur)."""
    if not numerisch:
        if zaehlung == "masse":
            p3 = 0.5 - math.atan(math.sqrt(5)) / (2 * math.pi) - math.sqrt(5) / (3 * math.pi)
            return dict(knoten=p3, faden=0.5 - p3, wand=0.5 - p3, leer=p3, l1=3 / math.sqrt(10 * math.pi), spur=1 / 3)
        w = 6 * math.sqrt(6)
        return dict(knoten=(29 - w) / 116, faden=(29 + w) / 116, wand=(29 + w) / 116, leer=(29 - w) / 116)
    from scipy.integrate import tplquad
    C, L = 15**3 / (8 * math.pi * math.sqrt(5)), 4.0

    def p(l3, l2, l1):
        I1, I2 = l1 + l2 + l3, l1 * l2 + l2 * l3 + l1 * l3
        g = abs(l1 * l2 * l3) if zaehlung == "stellen" else 1.0
        return g * C * math.exp(-3 * I1 * I1 + 7.5 * I2) * (l1 - l2) * (l2 - l3) * (l1 - l3)

    def teil(f, n):                                 # λ₁ ≥ λ₂ ≥ λ₃ mit n positiven Eigenwerten
        a1, b1 = (0, L) if n > 0 else (-L, 0)
        a2, b2 = (lambda x: 0, lambda x: x) if n > 1 else (lambda x: -L, lambda x: min(x, 0))
        a3, b3 = (lambda x, y: 0, lambda x, y: y) if n > 2 else (lambda x, y: -L, lambda x, y: min(y, 0))
        return tplquad(f, a1, b1, a2, b2, a3, b3, epsabs=1e-13, epsrel=1e-10)[0]
    k = [teil(p, n) for n in (3, 2, 1, 0)]
    alle = sum(k)
    aus = dict(knoten=k[0] / alle, faden=k[1] / alle, wand=k[2] / alle, leer=k[3] / alle)
    if zaehlung == "masse":
        aus["l1"] = sum(teil(lambda l3, l2, l1: l1 * p(l3, l2, l1), n) for n in range(4)) / alle
        aus["spur"] = (sum(teil(lambda l3, l2, l1: (l1 + l2 + l3)**2 / 3 * p(l3, l2, l1), n) for n in range(4))
                       / sum(teil(lambda l3, l2, l1: (l1 * l1 + l2 * l2 + l3 * l3) * p(l3, l2, l1), n) for n in range(4)))
    return aus


def netz_drehung(numerisch=False, n=400000, seed=1):
    """Verwindung am Netz. Linear ist J_i ∝ ε_ijk·(T·I)_jk, der Kommutator der Scherung T mit der Form I (09). In den
    Achsen von T, T = diag(λ₁, λ₂, λ₃), ist J_k = (λ_i − λ_j)·I_ij: Die Spur vertauscht mit jeder Form, und jede Achse
    dreht mit dem Unterschied der beiden anderen. Ist die Form gegen die Achsen zufällig verdreht, trägt jede Achse von
    ⟨J²⟩ den Anteil ⟨(λ_i − λ_j)²⟩/⟨Σ(λ_i − λ_j)²⟩. Die Scherung ohne Spur, ⟨Σμ²⟩ = 2σ²/3, ist μ_k = √(2/3)·ρ·cos(φ − 2πk/3)
    mit dem Gewicht |sin 3φ| unabhängig von ρ: ⟨(λ₁ − λ₃)²⟩ = 2⟨ρ²⟩·⟨sin²(φ + π/3)⟩ = (19/15)·σ², die äußeren Paare je
    (11/30)·σ², zusammen 2σ². Die mittlere Achse trägt 19/30, jede äußere 11/60, für jedes Spektrum und jede Skala.
    numerisch: T Gaußsch mit ⟨T_ij·T_kl⟩ = (σ²/15)·(δ_ij·δ_kl + δ_ik·δ_jl + δ_il·δ_jk) und I unabhängig ebenso gezogen,
    J = ε·(T·I) in den Achsen von T gemittelt, ohne die Zerlegung in Eigenwerte vorauszusetzen. Gibt (λ₁, λ₂, λ₃) als
    Anteile von ⟨J²⟩, dazu ⟨(λ₁ − λ₃)²⟩/σ²."""
    if not numerisch:
        from scipy.integrate import quad
        w = quad(lambda f: math.sin(3 * f), 0, math.pi / 3)[0]
        s2 = quad(lambda f: math.sin(3 * f) * math.sin(f + math.pi / 3)**2, 0, math.pi / 3)[0] / w
        d13 = 2 * (2 / 3) * s2
        aussen = (2 - d13) / 2 / 2
        return aussen, d13 / 2, aussen, d13
    import numpy as np
    rng = np.random.default_rng(seed)

    def feld():
        A = rng.standard_normal((n, 3, 3))
        S = (A + np.transpose(A, (0, 2, 1))) / 2
        return math.sqrt(2 / 15) * S + rng.standard_normal(n)[:, None, None] / math.sqrt(15) * np.eye(3)
    T, I = feld(), feld()
    TI = T @ I
    J = np.stack([TI[:, 1, 2] - TI[:, 2, 1], TI[:, 2, 0] - TI[:, 0, 2], TI[:, 0, 1] - TI[:, 1, 0]], axis=1)
    lam, v = np.linalg.eigh(T)
    Jk2 = np.mean(np.einsum("nik,ni->nk", v, J)**2, axis=0)
    teil = Jk2 / Jk2.sum()
    return float(teil[2]), float(teil[1]), float(teil[0]), float(np.mean((lam[:, 2] - lam[:, 0])**2))


def zwei_schwellen(numerisch=False, n=200000, dS=5e-4, seed=2):
    """Bereich und Leerraum als zwei Ränder derselben Irrfahrt: Über die Skalen läuft die lineare Störung eines Baryons in
    S = σ² eine Irrfahrt (scharfe Grenze in k) und kreuzt zuerst die Schwelle des Bereichs δ_c (rollender Kreis ohne
    Rückstellung) oder die des Leerraums δ_v (Kreuzung seines Randes, leerraum_rand). Ein Leerraum im Bereich wird
    zerdrückt, ein Bereich im Leerraum bleibt. Über alle Skalen kreuzt jede Irrfahrt einen der Ränder, den des Leerraums
    nach dem Ruin mit δ_c/(δ_c + |δ_v|), für jedes Spektrum, denn beide Schwellen wachsen mit demselben D. numerisch:
    zweiter Weg, Irrfahrten mit Gaußschen Schritten dS gewürfelt bis zum ersten Rand. Gibt (Leerraum, Bereich, δ_c, δ_v)."""
    dc, dv = 0.15 * (12 * math.pi)**(2 / 3), leerraum_rand()[0]
    if not numerisch:
        return dc / (dc - dv), -dv / (dc - dv), dc, dv
    import numpy as np
    rng = np.random.default_rng(seed)
    x, leer, s = np.zeros(n), 0, math.sqrt(dS)
    while x.size:
        x = x + s * rng.standard_normal(x.size)
        leer += int(np.count_nonzero(x <= dv))
        x = x[(x > dv) & (x < dc)]
    return leer / n, 1 - leer / n, dc, dv


def schale_im_spiegel(E, GM, lauf):
    """Schale ohne Rückstellung, ½ṙ² − GM/r = E je Masse, M die Masse des Feldes innerhalb. Im Spiegel r = z² mit dem
    Lauf dλ = dt/r ist sie der gehaltene Körper z″ = −K·z mit K = −E/2 und (dz/dλ)² + K·z² = GM/2, ab r = 0 bei λ = 0.
    Gebunden (K > 0) z = √(GM/2K)·sin(√K·λ), der rollende Kreis mit θ = 2√K·λ, die mittlere Schale (K = 0)
    z = √(GM/2)·λ mit r³ = (9/2)·GM·t², der Leerraum (K < 0) mit sinh. Gibt (r, t)."""
    K = -E / 2
    if K > 0:
        k = math.sqrt(K)
        z, t = math.sqrt(GM / (2 * K)) * math.sin(k * lauf), GM / (4 * K) * (lauf - math.sin(2 * k * lauf) / (2 * k))
    elif K < 0:
        k = math.sqrt(-K)
        z, t = math.sqrt(GM / (-2 * K)) * math.sinh(k * lauf), GM / (-4 * K) * (math.sinh(2 * k * lauf) / (2 * k) - lauf)
    else:
        z, t = math.sqrt(GM / 2) * lauf, GM * lauf**3 / 6
    return z * z, t


def umkehr_heute(th=None, a0_=a0):
    """Umkehr einer Schale um einen Bereich, dessen Gas nie abkühlte: Bis zur Umkehr atmet die Mode mit und wirkt wie
    Staub, im Atmen des Formteils mit der Rückstellung, r̈ = −(1 + Q)·GM_b/r² + H_∞²·r mit den Baryonen M_b
    innerhalb der Schale, wachsende Mode vom Anfang. In Einheiten H_∞ = 1, (1 + Q)·GM_b = ½ ist R_max³ = ½, das Atmen
    a³ = sinh²(3t/2) und die Laufzeit bis zur Umkehr bei R ∫dr/√(1/r − 1/R + r² − R²), mit r = R(1 − y²) ohne Wurzel am
    Umkehrpunkt. Die Schale, die bei Θ umkehrt (Vorgabe heute, H₀ aus dem Atmen), t = (2/3)·artanh(Θ/Θ_∞), umschließt
    für jede Masse dieselbe Verdichtung. Gibt (R_ta/R_max, Verdichtung gegen das Mittel der Faltung, gegen das Mittel
    in R_max)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    th = theta_(a0_, atmen_heute(a0_)[0] * 1e3 / Mpc) if th is None else th
    t = 2 / 3 * math.atanh(th * math.pi / math.sqrt(15))

    def zeit(R):
        def g(y):
            r = R * (1 - y * y)
            return 0.0 if r <= 0 else 2 * R / math.sqrt(1 / r - R * (R + r))
        return quad(g, 0, 1, limit=200)[0]
    R = brentq(lambda x: zeit(x) - t, 1e-6, 0.5**(1 / 3) * (1 - 1e-9))
    return 2**(1 / 3) * R, math.sinh(1.5 * t)**2 / R**3, 1 / (2 * R**3)


def groesste_umkehr(M, a0_=a0):
    """Größte Umkehr eines gebundenen Bereichs: An der ruhenden Schale (θ̃ = 0, w = 0) hält die Neigung des Tiefbereichs,
    v_flach²/R, der Rückstellung des Endzustands aus drei Vakuumteilen, H_∞²·R mit H_∞ = c·√Λ_s, Λ_s = π⁴/(15Λ²), die Waage:
    R_max = v_flach/H_∞ = (√15/2π²)·√(r*·ℓ₀). Gibt (R_max, H_∞)."""
    H_inf = c * math.pi**2 / (math.sqrt(15) * schliesslaenge(a0_))
    return v_flach(M, a0_) / H_inf, H_inf


def gezeit_ganzes(th=None, a0_=a0):
    """Die Rückstellung als Gezeit des Ganzen: Das Feld H_∞²·r dehnt in jeder Richtung gleich, 𝒯 = H_∞², mit der
    Gezeitendichte ρ_J = 3H_∞²/(4πG) = 2b₁·ρ_s. R_max³ = GM_dyn/H_∞² hat die Form der Glättungslänge ℓ³ = GM/𝒯. Eine Kugel
    vom Mittel der Faltung läuft wie jede Schale, r̈ = −GM/r² + H_∞²·r, mit r̈/r = ä/a. Im Atmen des Formteils ohne Strahlung
    ist ρ_J/ρ̄ = 2x²/(1 − x²) mit x = Θ/Θ_∞, eins bei Θ_∞/√3, wo die Dehnung zu beschleunigen beginnt, und die Kugel liegt
    bei r̄ = (ρ_J/ρ̄)^{1/3}·R_max ihrer Masse. th: Θ (Vorgabe heute aus dem Atmen). Gibt (ρ_J [kg/m³], ρ_J/ρ̄, r̄/R_max)."""
    H_inf = groesste_umkehr(M_sun, a0_)[1]
    th = theta_(a0_, atmen_heute(a0_)[0] * 1e3 / Mpc) if th is None else th
    x2 = (th * math.pi / math.sqrt(15))**2
    v = 2 * x2 / (1 - x2)
    return gezeitendichte(H_inf**2), v, v**(1 / 3)


def schwelle_umkehr(a0_=a0):
    """Treffpunkt von Schwelle und größter Umkehr ohne die Kreise der Falte. Die Rückstellung H_∞²·r erreicht a₀ bei
    r_× = a₀/H_∞² = (Θ_∞/2π)²·ℓ₀, dem Spiegelbild der Länge der Wirkung am Horizont des Endzustands ohne Masse,
    r_×·ℓ₀ = (c/H_∞)². Die größte Umkehr liegt auf dem Weg von r* zu r_× beim Bruchteil 1/d_eff,
    R_max = r*^{1 − 1/d_eff}·r_×^{1/d_eff}: nach der Regel (d_eff = 2) √(r*·r_×), nach Newton (d_eff = 3) (r*²·r_×)^{1/3},
    also R_max/r* = (M_×/M)^{1/(2d_eff)}. Bei M_× = a₀³/(G·H_∞⁴) = (Θ_∞/2π)⁴·M_U ist r* = r_× = R_max in beiden Formen,
    v_flach = a₀/H_∞. Gilt für die beiden Grenzformen durch (r*, a₀), das volle Feld der Regel ist an der Schwelle
    (e/(e − 1))·a₀. Gibt (M_× [kg], r_× [m], v_flach [m/s])."""
    H = groesste_umkehr(M_sun, a0_)[1]
    return a0_**3 / (G * H**4), a0_ / H**2, a0_ / H


def flucht_regel(M, r, a0_=a0, rand=None):
    """Fluchtgeschwindigkeit [m/s] aus dem Abstand r um die Punktmasse M nach der Regel (Galaxie mit übergebener Mode). Das
    Potential der Regel ist Φ = v_flach²·ln n_B(s), n_B = 1/(e^s − 1), s = r*/r (dΦ/dr = ν·a_N), dazu die Rückstellung des
    Endzustands −½H_∞²r². Der Rand ist die ruhende Schale ν(a_N)·a_N = H_∞²·R_e (tief R_max = v_flach/H_∞):
    v² = 2v_flach²·ln(n_B(s_e)/n_B(s)) − H_∞²·(R_e² − r²), ln(e^x − 1) als x + ln(1 − e^{−x}) gegen Überlauf. rand [m]: statt
    der ruhenden Schale eine feste äußere Grenze R_e, bis zu der die Flucht gezählt wird."""
    from scipy.optimize import brentq
    vf, rs = v_flach(M, a0_), schwellenradius(M, a0_)
    R_max, H = groesste_umkehr(M, a0_)
    Re = rand if rand is not None else brentq(lambda x: verstaerkung(G * M / x**2, a0_) * G * M / x**2 - H**2 * x, r,
                                               10 * R_max)
    ln_innen = rs / r + math.log1p(-math.exp(-rs / r))
    return math.sqrt(2 * vf**2 * (ln_innen - math.log(math.expm1(rs / Re))) - H**2 * (Re**2 - r**2))


def verdampfen_anteil(X):
    """Anteil einer Maxwell-Verteilung jenseits der Flucht, X = v_esc/(√2·σ): erfc X + (2X/√π)·e^{−X²}. So viel hebt jede
    Relaxation über den Rand. Newton, mittlere Flucht 4⟨v²⟩ eines Bereichs im Gleichgewicht: X² = 6. Um einen Bereich mit
    übergebener Mode setzt die Gezeit den Rand R_e, tief d_eff·Ω²·r_J² = H_∞²·R_max² = v_flach², mit σ² = v_flach²/γ ist
    X² = (γ/2)·[2·ln(n_B(r*/R_e)/n_B(s)) − 1 + (r/R_e)²], e^{−X²} das Verhältnis der Rückkehrzahlen in der Potenz γ."""
    return math.erfc(X) + 2 * X / math.sqrt(math.pi) * math.exp(-X * X)


def gezeitenrand(m, M, D, m_B=None, a0_=a0):
    """Jacobi-Radius r_J [m] eines Bereichs der Masse m auf der Kreisbahn im Abstand D um die Punktmasse M nach der Regel, in
    erster Ordnung in r/D. Gezeit des Wirts 𝒯 = −dg/dD = (d_eff − 1)·Ω², mit dem Umlauf d_eff·Ω², d_eff = 3 − χ. Das eigene
    Feld: Newton aus den heutigen Baryonen m, die Rückkehr aus den einmal gebundenen m_B (Vorgabe m),
    g = g_N + (ν(g_N,B) − 1)·g_N,B. r_J aus g(r_J) = d_eff·Ω²·r_J, tief r_J = v_flach/(√d_eff·Ω)."""
    from scipy.optimize import brentq
    m_B = m if m_B is None else m_B
    aN = G * M / D**2
    T = (3 - schalter(aN, a0_)) * verstaerkung(aN, a0_) * aN / D
    g = lambda r: G * m / r**2 + (verstaerkung(G * m_B / r**2, a0_) - 1) * G * m_B / r**2
    return brentq(lambda r: g(r) - T * r, 1e-6 * D, D)


def zeitargument(M1, M2, r, v_r, H=None, ein_bereich=False, v_t=0.0, gegenzug=False, numerisch=False):
    """Zeitargument zweier Galaxien mit übergebener Mode aus gemeinsamem Anfang, Abstand r, Annäherung v_r, quer zur
    Verbindung v_t (Vorgabe 0: radial). Jede läuft als Ganzes im verstärkten Feld der anderen, a = ν(g₁)·g₁ + ν(g₂)·g₂ mit
    g_i = GM_i/r² (ein_bereich: das Paar als ein Bereich, ν(g)·g mit g = G(M₁ + M₂)/r², gegenzug: der leichtere Bereich
    fällt als Ganzes im verstärkten Feld des schwereren, und der erfährt den Gegenzug, a = ν(g_s)·g_s·(1 + M_l/M_s) mit
    M_s ≥ M_l, der Impuls bleibt), dazu die Rückstellung −H²·r (Vorgabe H_∞, H = 0 ohne). Die Beschleunigung ist zentral,
    der Drehimpuls L = r·v_t bleibt, die Wendepunkte sind die Nullstellen von ½(v_r² + v_t²) − (U(x) − U(r)) − L²/2x². Tief
    ist a = (v₁² + v₂²)/r, und radial vergeht von der Umkehr R bis zur Mitte √(π/2)·R/v. Außerhalb der Baryonen entsteht
    kein Kielwasser, die Bahn läuft ohne Verlust. Gibt für die zweite Annäherung (Umkehr [m], Zeit seit der größten Nähe
    [s], seit dem Anfang [s], größte Nähe [m], radial 0). Die nächste größte Nähe folgt nach (seit dem Anfang) − 2·(seit der
    größten Nähe). numerisch: zweiter Weg, die Bahn rückwärts in der Zeit integriert (in der Ebene, radial mit Abprall an
    der Mitte) bis zu den Wendepunkten."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    H = groesste_umkehr(M_sun)[1] if H is None else H
    if gegenzug:
        ml, ms = sorted((M1, M2))
        acc = lambda x: verstaerkung(G * ms / x**2) * G * ms / x**2 * (1 + ml / ms) - H**2 * x
    else:
        massen = (M1 + M2,) if ein_bereich else (M1, M2)
        acc = lambda x: sum(verstaerkung(G * m / x**2) * G * m / x**2 for m in massen) - H**2 * x
    if numerisch:
        from scipy.integrate import solve_ivp
        if v_t > 0:
            f2 = lambda t, y: (lambda d: [y[2], y[3], -acc(d) * y[0] / d, -acc(d) * y[1] / d])(math.hypot(y[0], y[1]))
            wende = lambda t, y: y[0] * y[2] + y[1] * y[3]
            s = solve_ivp(f2, (0, 1e19), [r, 0.0, v_r, v_t], events=wende, rtol=1e-11, atol=1e-6, method="DOP853",
                          max_step=1e15, dense_output=False)
            tw, yw = s.t_events[0], s.y_events[0]
            d = [math.hypot(y[0], y[1]) for y in yw]
            return d[0], tw[1], tw[3], d[1]
        eps = 1e-4 * r
        mitte = lambda t, y: y[0] - eps
        mitte.terminal, mitte.direction = True, -1
        umkehr = lambda t, y: y[1]
        umkehr.direction = -1
        f1 = lambda t, y: [y[1], -acc(y[0])]
        a = solve_ivp(f1, (0, 1e19), [r, v_r], events=(mitte, umkehr), rtol=1e-11, atol=1e-6, method="DOP853")
        b = solve_ivp(f1, (0, 1e19), [eps, -a.y_events[0][0][1]], events=(mitte, umkehr), rtol=1e-11, atol=1e-6,
                      method="DOP853")
        return a.y_events[1][0][0], a.t_events[0][0], a.t_events[0][0] + b.t_events[0][0], 0.0
    U = lambda x: quad(acc, r, x, limit=200)[0]
    f = lambda x: 0.5 * (v_r**2 + v_t**2) - U(x) - (r * v_t)**2 / (2 * x * x)
    ra = brentq(f, r, 3 * r)
    rp = brentq(f, 1e-9 * r, r * (1 - 1e-9)) if v_t > 0 else 0.0
    v = lambda x: math.sqrt(max(2 * f(x), 1e-300))
    T = quad(lambda x: 1 / v(x), max(rp, 1e-6 * r), ra, limit=500, points=[r])[0]
    zurueck = quad(lambda x: 1 / v(x), r, ra, limit=500)[0]
    return ra, T + zurueck, 3 * T + zurueck, rp


def fluss_um_bereich(massen, r, t, form="mitatmend", H=None, Q=None):
    """Radiale Geschwindigkeit [m/s] einer Schale im Abstand r [m] um einen Bereich aus Galaxien mit den Baryonen massen
    [kg] zur Zeit t [s] nach dem gemeinsamen Anfang, mit der Rückstellung −H²·r (Vorgabe H_∞). form "mitatmend": bis zur
    ersten Kreuzung trägt die Schale die Mode mit, a = (1 + Q)·G·ΣM/r²; "regel": jede Galaxie mit ihrer Regel,
    a = Σ ν(g_i)·g_i. Aus der Laufzeit ∫dx/√(2(E − U)) vom Anfang: noch vor der Umkehr nach außen (> 0), danach im Fall
    (< 0). Die Wurzel am Umkehrpunkt hebt x = R(1 − y²) heraus."""
    from scipy.integrate import quad, cumulative_trapezoid
    from scipy.optimize import brentq
    import numpy as np
    H = groesste_umkehr(M_sun)[1] if H is None else H
    Q = Q_kosm if Q is None else Q
    if form == "mitatmend":
        U = lambda x: -(1 + Q) * G * sum(massen) / x - 0.5 * H**2 * x**2
    else:
        xs = np.geomspace(1e-6 * r, 50 * r, 6000)
        acc = np.array([sum(verstaerkung(G * m / x**2) * G * m / x**2 for m in massen) for x in xs])
        Ug = np.concatenate([[0.0], cumulative_trapezoid(acc, xs)]) - acc[0] * xs[0] * math.log(xs[0] / r)
        U = lambda x: float(np.interp(math.log(x), np.log(xs), Ug)) - 0.5 * H**2 * x**2

    def bis(R, E, a):
        def g(y):
            d = E - U(R * (1 - y * y))
            return 2 * R * y / math.sqrt(2 * d) if d > 0 else 0.0
        return quad(g, 0, math.sqrt(1 - a / R), limit=400)[0]
    T_um = lambda R: bis(R, U(R), 1e-6 * R)
    if H == 0:
        R_m = math.inf
    elif form == "mitatmend":
        R_m = ((1 + Q) * G * sum(massen) / H**2)**(1 / 3)
    else:
        R_m = brentq(lambda x: sum(verstaerkung(G * m / x**2) * G * m / x**2 for m in massen) - H**2 * x, 1e-3 * r, 1e3 * r)
    if r >= R_m or t < T_um(r):
        E0 = U(min(r, R_m))
        T = lambda ls: bis(r, E0 + math.exp(ls), 1e-6 * r) - t
        return math.sqrt(2 * (E0 + math.exp(brentq(T, math.log(1e2), math.log(1e16))) - U(r)))
    Rt = brentq(lambda R: T_um(R) + bis(R, U(R), r) - t, r * (1 + 1e-6), min(R_m * (1 - 1e-3), 1e3 * r))
    return -math.sqrt(2 * (U(Rt) - U(r)))


def zeitargument_newton(M, r, t, H=0.0):
    """Zeitargument nach Newton mit der Masse M und der Rückstellung −H²·r: Annäherung [m/s] im Abstand r nach der Zeit t
    aus gemeinsamem Anfang, auf dem ersten Fall nach der Umkehr, aus der Laufzeit ∫dx/√(2(E − U)) mit U = −GM/x − ½H²x²."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    U = lambda x: -G * M / x - 0.5 * H**2 * x**2
    R_m = (G * M / H**2)**(1 / 3) if H > 0 else math.inf
    v_max = math.sqrt(2 * ((U(R_m) if H > 0 else 0.0) - U(r)))
    def laufzeit(v):
        E = 0.5 * v**2 + U(r)
        ra = brentq(lambda x: U(x) - E, r, R_m if H > 0 else 1e6 * r)
        # x = ra·(1 − y²): E − U(x) = y²·[GM/(ra(1 − y²)) − ½H²ra²(2 − y²)], die Wurzel am Umkehrpunkt hebt sich heraus
        g = lambda y: 2 * ra / math.sqrt(2 * (G * M / (ra * (1 - y * y)) - 0.5 * H**2 * ra**2 * (2 - y * y)))
        bis = lambda x0: quad(g, 0, math.sqrt(1 - x0 / ra), limit=200)[0]
        return bis(1e-12 * r) + bis(r)
    return brentq(lambda v: laufzeit(v) - t, 1.0, v_max * (1 - 1e-6))


def naht_haufen(Delta=200, H0_=H0, Q=Q_kosm, a0_=a0):
    """Naht eines Haufens in geschlossener Form: Mit g = (1 + Q)·g_bar ist am Rand R_Δ (Δ·ρ_krit) a_N = Δ·H²·R/(2(1 + Q)),
    die Naht (a_N = (ln q)²·a₀, q = Q/(1 + Q)) liegt für jede Masse bei R = 2(1 + Q)·(ln q)²·a₀/(Δ·H²), und R_Δ trifft sie bei
    M = 4(1 + Q)³·(ln q)⁶·a₀³/(Δ²·G·H⁴). Gibt (R, M)."""
    s2 = math.log(Q / (1 + Q))**2
    R = 2 * (1 + Q) * s2 * a0_ / (Delta * H0_**2)
    return R, Delta / 2 * H0_**2 * R**3 / G


def abkuehlung_licht(z=0.0, Y=None):
    """Abkühlung heißen, ionisierten Gases am Licht der Hintergrundstrahlung: Je Stoß gewinnt das Photon im Mittel
    4kT/(m_e c²) seiner Energie, das Elektron gibt je Zeit 4σ_T·c·u_γ·kT/(m_e c²) ab. Gegen die Wärme je Elektron,
    (3/2)·(μ_e/μ)·kT, ist t = 3μ_e·m_e·c/(8μ·σ_T·u_γ), unabhängig von Dichte und Temperatur, u_γ ∝ (1 + z)⁴,
    σ_T = (8π/3)·(α·λ̄_e)². Gas der Urmischung (Y = Y_p), voll ionisiert: μ = 1/(2X + 3Y/4), μ_e = 2/(1 + X), X = 1 − Y.
    Gibt t [s]."""
    Y = YP_GEMESSEN[0] if Y is None else Y
    X = 1 - Y
    mu, mu_e = 1 / (2 * X + 0.75 * Y), 2 / (1 + X)
    sigma_T = 8 * math.pi / 3 * (alpha * compton(m_e))**2
    u = dichte_photonen() * c**2 * (1 + z)**4
    return 3 * mu_e * m_e * c / (8 * mu * sigma_T * u)


def kuehlung_licht(Q=None, x_e=None):
    """Kühlung am Licht gegen das Atmen. Ionisiertes Gas gibt seine Wärme in t = 3μ_e·m_e·c/(8μ·σ_T·u_γ) an das Licht
    (abkuehlung_licht), t ∝ T⁻⁴ in der Temperatur T des Lichts: Windung, σ_T = (8π/3)·(α·λ̄_e)², und die Stufe d = 3 des
    Lichts, u_γ ∝ T⁴. Die Zeit des Atmens 1/H ist, solange die Faltung das Atmen trägt, ∝ T^{−3/2}, denn
    H² = (8πG/3)·(1 + Q)·ε_b·n_γ/c² mit der zeitfesten Faltung je Photon ε_b. Gleich sind beide bei einer Temperatur des
    Ganzen, (kT_ab/m_e c²)^{5/2} = (135/(64π³α²))·(n/n_e)·√((16ζ(3)/3π)·(1 + Q)·α_G(m_e)·ε_b/(m_e c²)), α_G(m_e) = (m_e/m_Pl)²,
    n/n_e = μ_e/μ Teilchen je freiem Elektron der voll ionisierten Urmischung. z_ab gerechnet aus t(z)·H(z) = 1 mit Strahlung
    und Formteil b₁·ρ_s, dazu die geschlossene Form. x_e: freie Elektronen je Wasserstoff, Helium neutral,
    n/n_e = (1 + f_He + x_e)/x_e mit f_He = Y_p/(4X). Das Gas trägt die Wärme des Lichts bis t·H = 1, geschlossen
    (1 + z_t)^{5/2} = (1 + z_ab)^{5/2}·(n/n_e)/(μ_e/μ). Bei festem x_e und H ∝ (1 + z)^{3/2} ist y = T_Gas/T_Licht die Lösung
    von dy/du = −y + e^{−5u/2}·(1 − y), u = ln(a/a_t): bei z_t ist y = e^{2/5}·(2/5)^{2/5}·Γ(3/5, 2/5), spät
    y → (2/5)^{2/5}·Γ(3/5)·e^{−u}. Q = None: gemessen (ω_c/ω_b). Gibt dict (z_ab, kT_ab [eV], T_ab [K], z_ab_form,
    kT_ab_form [eV], n_je_e, y_t, mit x_e dazu z_t und z_t_form)."""
    from scipy.optimize import brentq
    from scipy.special import gamma, gammaincc
    Qw = Q_kosm if Q is None else Q
    Y = YP_GEMESSEN[0]
    X = 1 - Y
    mu, mu_e = 1 / (2 * X + 0.75 * Y), 2 / (1 + X)
    rho_r, _, rho_m = dichte_heute(Q)
    rho_de = dichte_dunkle_energie()
    H = lambda z: math.sqrt(8 * math.pi * G / 3 * (rho_r * (1 + z)**4 + rho_m * (1 + z)**3 + rho_de))
    z_ab = brentq(lambda z: abkuehlung_licht(z) * H(z) - 1, 0.0, 1e3, xtol=1e-12)
    eps_b = faltung_je_photon()["eps_b"] * e
    zeta3 = strahlung_je_sinn(1, 2) / 2
    mec2 = m_e * c**2
    k52 = (135 / (64 * math.pi**3 * alpha**2) * mu_e / mu
           * math.sqrt(16 * zeta3 / (3 * math.pi) * (1 + Qw) * (m_e / m_Pl)**2 * eps_b / mec2))
    kT_form = k52**0.4 * mec2
    aus = dict(z_ab=z_ab, kT_ab=k_B * T_CMB * (1 + z_ab) / e, T_ab=T_CMB * (1 + z_ab),
               z_ab_form=kT_form / (k_B * T_CMB) - 1, kT_ab_form=kT_form / e, n_je_e=mu_e / mu,
               y_t=float(math.exp(0.4) * 0.4**0.4 * gammaincc(0.6, 0.4) * gamma(0.6)))
    if x_e is not None:
        n_je_e = (1 + Y / (4 * X) + x_e) / x_e
        aus["z_t"] = brentq(lambda z: n_je_e * mu / mu_e * abkuehlung_licht(z) * H(z) - 1, 0.0, 1e4, xtol=1e-10)
        aus["z_t_form"] = (1 + aus["z_ab_form"]) * (n_je_e * mu / mu_e)**0.4 - 1
    return aus


def signal_21cm(z=17.0, Q=None):
    """Gas und Licht nach der Rekombination und die Tiefe des 21-cm-Signals. Den Rest der Ionisation x_e(z) gibt die
    Rekombination in CAMB (ω_b, ω_c = Q·ω_b, H₀ aus a₀). Das Gas trägt die Wärme des Lichts, bis seine Kühlung am Licht mit
    n/n_e = (1 + f_He + x_e)/x_e die Zeit des Atmens erreicht (kuehlung_licht), bei z_t, danach fällt T_Gas wie (1 + z)²:
    scharf T_Licht·(1 + z)/(1 + z_t), mit dem Übergang y = (2/5)^{2/5}·Γ(3/5) mal so viel, zweiter Weg die Temperatur der
    Baryonen von CAMB. Die Linie geht mit A₁₀ = μ₀ω³μ_B²/(3πħc³) bei ν₂₁, die Abgabe der Windung des Elektrons. Folgt der
    Spin dem Gas, T_S = T_Gas, ist die Aufnahme τ = 3c³ħA₁₀n_HI/(16kν₂₁²T_S·H) und das Signal
    δT = (T_S − T_Licht)·(1 − e^{−τ})/(1 + z), das tiefste, das ohne fremde Wärme möglich ist. Gibt dict (z_t, x_e dort,
    T_Gas scharf, mit y, CAMB [K], T_Licht [K], A10 [1/s], tau, dT scharf, mit y, CAMB [mK], Frequenz [MHz])."""
    import camb
    from scipy.optimize import brentq
    Qw = Q_kosm if Q is None else Q
    par = camb.set_params(ombh2=omega_b, omch2=Qw * omega_b, H0=atmen_heute(Q=Q)[0], As=2.1e-9, ns=0.965, tau=0.054)
    r = camb.get_background(par)
    ev = lambda zz, n: float(r.get_background_redshift_evolution([zz], [n], format="array")[0, 0])
    Y = YP_GEMESSEN[0]
    X = 1 - Y
    mu, mu_e = 1 / (2 * X + 0.75 * Y), 2 / (1 + X)
    rho_r, rho_b, rho_m = dichte_heute(Q)
    H = lambda zz: math.sqrt(8 * math.pi * G / 3 * (rho_r * (1 + zz)**4 + rho_m * (1 + zz)**3 + dichte_dunkle_energie()))
    zt = brentq(lambda zz: (1 + Y / (4 * X) + ev(zz, "x_e")) / ev(zz, "x_e") * mu / mu_e * abkuehlung_licht(zz) * H(zz) - 1,
                20, 1000, xtol=1e-8)
    TL = T_CMB * (1 + z)
    T_s = TL * (1 + z) / (1 + zt)
    T_y = 0.4**0.4 * math.gamma(0.6) * T_s
    T_c = ev(z, "T_b")
    w = 2 * math.pi * NU_21CM
    A10 = mu0 * w**3 * (e * hbar / (2 * m_e))**2 / (3 * math.pi * hbar * c**3)
    nHI = X * rho_b / m_H * (1 + z)**3 * (1 - ev(z, "x_e"))
    tau = lambda TS: 3 * c**3 * hbar * A10 * nHI / (16 * k_B * NU_21CM**2 * TS * H(z))
    dT = lambda TS: (TS - TL) * (1 - math.exp(-tau(TS))) / (1 + z) * 1e3
    return dict(z_t=zt, x_e=ev(zt, "x_e"), T_scharf=T_s, T_y=T_y, T_camb=T_c, T_licht=TL, A10=A10, tau=tau(T_s),
                dT_scharf=dT(T_s), dT_y=dT(T_y), dT_camb=dT(T_c), MHz=NU_21CM / (1 + z) / 1e6)


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


def grenzmasse(mu_e=2.0, beta=None, mu=None):
    """Halter bei c, P = b·ħc·(ρ/m)^{4/3}: Form n = 3, M = (4ω₃/√π)·b^{3/2}·m_Pl³/m², ω₃ = 2,018 der Beiwert der Form
    (lane_emden). Zählregel der Elektronen nach Fermi, b = (3π²)^{1/3}/4, m = μ_e·m_u: M = (ω₃·√(3π)/2)·m_Pl³/(μ_e·m_u)².
    Mit beta: Licht nach Bose mit der Wärme des Gases, β der Anteil des Gases am Druck im ganzen Trichter,
    b = (45·(1 − β)/(π²β⁴))^{1/3}, m = μ·m_H."""
    if beta is None:
        b, m = (3 * math.pi**2)**(1 / 3) / 4, mu_e * m_u
    else:
        b, m = (45 * (1 - beta) / (math.pi**2 * beta**4))**(1 / 3), mu * m_H
    return 4 * lane_emden(3.0)[1] / math.sqrt(math.pi) * b**1.5 * (hbar * c / G)**1.5 / m**2


def atmen_gehalten(teile, n=2.0, beta=None):
    """Atmen eines gehaltenen Körpers, alle Abstände gleich, r → (1 + ε)·r (01). Die Energie der Träger, ∝ pⁿ, fällt dabei
    wie im Atmen des Ganzen mit (1 + ε)^{−n}, 3∫P dV = n·E_kin. Jede Neigung, die wie (1 + ε)^{1 − d_i} geht (Flussdimension
    d_i), trägt ⟨r·g_i⟩ ∝ (1 + ε)^{2 − d_i}, im Gleichgewicht ist n·E_kin = M·Σ⟨r·g_i⟩. Dann ist ε̈ + K·ε = 0 mit
    K = Σ(n + 2 − d_i)·⟨r·g_i⟩/⟨r²⟩. Über ln(1 + ε) sind Halter (Steigung −n) und Neigung (2 − d_i) Geraden, das
    Gleichgewicht ihr Schnittpunkt, K der Unterschied ihrer Steigungen: dreht für K > 0, schert parallel bei d = n + 2
    (Ruhe nur bei einer Masse, dann bei jedem Radius), dehnt für K < 0. Kreisbahn: Halter der Umlauf L²/r³, n = 2,
    K = κ² = Σ(4 − d_i)·Ω_i², die Rückstellung (d = 0, ⟨r·g⟩/⟨r²⟩ = −H_∞²) zählt vierfach. Newton d = 3:
    K = (n − 1)·⟨r·g⟩/⟨r²⟩, Bindung je Masse −((n − 1)/n)·⟨r·g⟩, gleichmäßig ⟨r·g⟩/⟨r²⟩ = GM/R³. Halter bei c, n = 1,
    parallel zu Newton. Regel d = 3 − χ̄ (χ̄ gewichtet mit r·g), tief d = 2. beta: Licht hält mit der Wärme des Gases
    (Anteil β des Gases am Druck), n − 1 = β·(4 − 3β)/(8 − 7β). teile: Folge von (d_i, ⟨r·g_i⟩/⟨r²⟩ [1/s²]).
    Gibt K [1/s²]."""
    if beta is not None:
        n = 1 + beta * (4 - 3 * beta) / (8 - 7 * beta)
    return sum((n + 2 - d) * x for d, x in teile)


def atmen_ganzes(n=2.0, numerisch=False):
    """Das Ganze als gehaltener Körper (atmen_gehalten): Es hat keine Hülle und keinen Halter außer seiner Rückstellung.
    Ohne Druck ruht es nur, wo die Faltung die Rückstellung hält, GM/R³ = H_∞², und dort ist für jeden Träger
    K = (n − 1)·H_∞² − (n + 2)·H_∞² = −3·H_∞²: Es dehnt mit der Rate √3·H_∞ und hat keine Periode. numerisch: zweiter
    Weg, r̈ = −GM/r² + H_∞²·r um die Ruhe abgeleitet, K = −dr̈/dr. Gibt K/H_∞²."""
    if not numerisch:
        return atmen_gehalten(((3, 1.0), (0, -1.0)), n=n)
    GM, H2 = 1.0, 1.0
    rE = (GM / H2)**(1 / 3)
    f = lambda r: -GM / r**2 + H2 * r
    h = 1e-6 * rE
    return -(f(rE + h) - f(rE - h)) / (2 * h) / H2


def tov_sperre(m_traeger, x_c=None, m_masse=None, faltend=True, traegheit=False):
    """Trichter, den die Zählregel eines Trägers der Masse m_traeger hält (freies Fermi-Gas, x = λ̄·k_F, zwei Sinne).
    Der Druck ist Energie und faltet mit: dP/dr = −(3c²Δs_r/r²)·(ε + P)/c²·(1 + 4πr³P/(M_r·c²))/(1 − 6Δs_r/r),
    dM_r/dr = 4πr²ε/c², Δs_r = GM_r/(3c²), mit u = m_traeger·c²/(8π²λ̄³), P = u·[x(2x²/3 − 1)√(1 + x²) + arsinh x],
    ε = u·[x√(1 + x²)(1 + 2x²) − arsinh x] + n·(m_masse − m_traeger)·c². m_masse ist die Masse je Teilchen des Halters
    (Vorgabe der Träger selbst, beim Weißen Zwerg 2m_u). faltend=False: Newton, ρ = n·m_masse. x_c=None (nur faltend):
    die größte Masse, das erste Maximum über x_c. Gibt (M [kg], R [m], δ = 3Δs/R, x_c), mit traegheit dazu
    k = I/(MR²), I = (8π/3)·∫ρ·r⁴dr mit ρ = ε/c² (faltend) oder n·m_masse (Newton)."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import minimize_scalar
    m_masse = m_traeger if m_masse is None else m_masse
    if x_c is None:
        x, alt, neu = 0.05, 0.0, tov_sperre(m_traeger, 0.05, m_masse)[0]
        while neu > alt:
            x, alt, neu = 1.25 * x, neu, tov_sperre(m_traeger, 1.25 * x, m_masse)[0]
        r = minimize_scalar(lambda lx: -tov_sperre(m_traeger, math.exp(lx), m_masse)[0],
                            bounds=(math.log(x / 1.5625), math.log(x)), method="bounded", options=dict(xatol=1e-9))
        return tov_sperre(m_traeger, math.exp(r.x), m_masse, traegheit=traegheit)
    lam = compton(m_traeger)
    u = m_traeger * c**2 / (8 * math.pi**2 * lam**3)

    def druck(x):
        if x < 1e-2:
            return u * 8 / 3 * (x**5 / 5 - x**7 / 14 + x**9 / 24)
        return u * (x * (2 * x**2 / 3 - 1) * math.sqrt(1 + x**2) + math.asinh(x))

    def energie_dichte(x):
        if x < 1e-2:
            eF = u * 8 * (x**3 / 3 + x**5 / 10 - x**7 / 56)
        else:
            eF = u * (x * math.sqrt(1 + x**2) * (1 + 2 * x**2) - math.asinh(x))
        return eF + x**3 / (3 * math.pi**2 * lam**3) * (m_masse - m_traeger) * c**2

    def f(r, y):
        x, Mr = max(y[0], 1e-12), y[1]
        p, eps = druck(x), energie_dichte(x)
        if faltend:
            dp = -3 * c**2 * faltung(Mr) / r**2 * (eps + p) / c**2 * (1 + 4 * math.pi * r**3 * p / (Mr * c**2)) \
                 / (1 - 6 * faltung(Mr) / r)
            dm = 4 * math.pi * r**2 * eps / c**2
        else:
            rho = x**3 / (3 * math.pi**2 * lam**3) * m_masse
            dp, dm = -G * Mr * rho / r**2, 4 * math.pi * r**2 * rho
        return [dp / (u * 8 / 3 * x**4 / math.sqrt(1 + x**2)), dm, 2 / 3 * r * r * dm]

    ende = lambda r, y: y[0] - 1e-6 * x_c
    ende.terminal = True
    skala = c / math.sqrt(G * energie_dichte(x_c) / c**2)
    r0 = 1e-6 * skala
    m0 = 4 / 3 * math.pi * r0**3 * energie_dichte(x_c) / c**2
    s = solve_ivp(f, [r0, 1e6 * skala], [x_c, m0, 0.4 * m0 * r0**2], method="LSODA",
                  events=ende, rtol=1e-10, atol=[1e-14 * x_c, 1e-12 * M_sun, 1e-12 * M_sun * r0**2], first_step=r0)
    M, R = s.y_events[0][0][1], s.t_events[0][0]
    if traegheit:
        return M, R, 3 * faltung(M) / R, x_c, s.y_events[0][0][2] / (M * R**2)
    return M, R, 3 * faltung(M) / R, x_c

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


def mitte_gehalten(M, R, numerisch=False):
    """Mitte eines ruhenden, gleichmäßig gehaltenen Körpers: Von außen gelesen ist (ρc² + P)·dτ/dt überall gleich, darum
    P/(ρc²) = dτ/dt|_R / dτ/dt|_0 − 1, und die Uhr der Mitte geht mit dτ/dt|_0 = (3/2)·dτ/dt|_R − ½. Für jede nach außen nicht
    steigende Dichte ist das die schnellste Uhr der Mitte. Sie steht bei dτ/dt|_R = ⅓, R = (9/8)·r₀ = 27Δs/4 (Buchdahl).
    numerisch: aus dem Gleichgewicht nach Tolman, Oppenheimer und Volkoff bei fester Dichte (G = c = ρ = 1), ohne die
    geschlossene Form. Gibt (dτ/dt|_0, P/(ρc²))."""
    if not numerisch:
        T = takt(M, R)
        T0 = 1.5 * T - 0.5
        return T0, T / T0 - 1 if T0 > 0 else float("inf")
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    ziel = 3 * faltung(M) / R

    def kugel(pc):
        def f(r, y):
            m = 4 * math.pi / 3 * r**3
            h = (m + 4 * math.pi * r**3 * y[0]) / (r**2 * (1 - 2 * m / r))
            return [-(1 + y[0]) * h, 2 * h]
        rand = lambda r, y: y[0]
        rand.terminal, rand.direction = True, -1
        s = solve_ivp(f, (1e-9, 10), [pc, 0.0], events=rand, rtol=1e-11, atol=1e-14, max_step=1e-3)
        Rk, nu = s.t_events[0][0], s.y_events[0][0][1]
        d = 4 * math.pi / 3 * Rk**2
        return d, math.exp((math.log(1 - 2 * d) - nu) / 2)
    pc = brentq(lambda p: kugel(p)[0] - ziel, 1e-8, 1e7, xtol=1e-14)
    return kugel(pc)[1], pc


def kappe(M, R, r=0.0, H=0.0):
    """Ruhender, gleichmäßig gehaltener Körper als Kappe einer Kugel S³ vom Radius A: c²/A² = (8πG/3)·ρ + H², ohne
    Rückstellung (H = 0) A² = R³/r₀. A misst die Krümmung des Inneren, keine Strecke im Raum. Hülle im Winkel χ_R vom Pol,
    R = A·sin χ_R (H = 0: sin²χ_R = 2δ), der Ort r ≤ R bei r = A·sin χ. Mit dem Weg nach außen wächst der Umfang wie
    dr/dl = cos χ, an der Hülle der Takt, mit H der des Endzustands √(1 − 6Δs/R − (H·R/c)²). Uhr am Ort r (nur H = 0):
    (3/2)·cos χ_R − ½·cos χ, in der Mitte steht sie bei cos χ_R = ⅓. Bogen vom Pol zur Hülle A·χ_R, schwach um R³/(6A²) = Δs
    länger als R. M = 0 mit H: halbe Kugel vom Radius c/H. Dieselbe Kappe hat ein gebundener Bereich im Feld GM/r² im
    Augenblick seiner Umkehr. Gibt dict(A [m], chi_R, chi [rad], takt, uhr, ueberschuss [m]), uhr None für H ≠ 0."""
    rho = M / (4 * math.pi / 3 * R**3)
    A = c / math.sqrt(8 * math.pi * G / 3 * rho + H**2)
    chi_R, chi = math.asin(min(R / A, 1.0)), math.asin(min(r / A, 1.0))
    uhr = 1.5 * math.cos(chi_R) - 0.5 * math.cos(chi) if H == 0 else None
    d = chi_R - math.sin(chi_R) if chi_R > 1e-2 else chi_R**3 / 6 * (1 - chi_R**2 / 20)
    return dict(A=A, chi_R=chi_R, chi=chi, takt=math.cos(chi_R), uhr=uhr, ueberschuss=A * d)

def polytrope(n, traegheit=False):
    """Lane-Emden θ'' + 2θ'/ξ = −θⁿ: Rand ξ₁, −θ'(ξ₁) und der Faktor der Kerntemperatur f = 1/((n + 1)·ξ₁·|θ'(ξ₁)|)
    in kT_c = f·δ·m̄c², δ = GM/(Rc²). traegheit: dazu k = I/(MR²) = (2/3)·∫θⁿξ⁴dξ/(ξ₁²·∫θⁿξ²dξ)."""
    from scipy.integrate import solve_ivp
    ende = lambda x, y: y[0]
    ende.terminal = True
    s = solve_ivp(lambda x, y: [y[1], -max(y[0], 0)**n - 2 * y[1] / x, max(y[0], 0)**n * x * x,
                                max(y[0], 0)**n * x**4], [1e-6, 50], [1 - 1e-12 / 6, -1e-6 / 3, 0.0, 0.0],
                  events=ende, rtol=1e-12, atol=1e-14)
    x1, dt1 = s.t_events[0][0], s.y_events[0][0][1]
    if traegheit:
        m2, m4 = s.y_events[0][0][2], s.y_events[0][0][3]
        return x1, -dt1, 1 / ((n + 1) * x1 * -dt1), 2 * m4 / (3 * x1**2 * m2)
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


def gezeitendichte(gezeit):
    """Gezeitendichte ρ_J = 3𝒯/(4πG) [kg/m³]: Was im Mittel dichter ist, bleibt Einzeltrichter, nur Dünneres wird
    geglättet. Gleichwertig: Ein Körper der Masse M ist kompakt, wenn sein Radius unter ℓ liegt, (4π/3)·ρ_J·ℓ³ = M."""
    return 3 * gezeit / (4 * math.pi * G)


def zerreiss_umlauf(rho, d_eff=3):
    """Umlaufzeit [s], unter der ein Körper der mittleren Dichte ρ [kg/m³] auf der Kreisbahn zerreißt: Er bleibt ganz,
    solange er dichter ist als die Gezeitendichte 3𝒯/(4πG) mit der Gezeit samt Umlauf 𝒯 = d_eff·Ω². Gleich sind beide bei
    Ω² = 4πGρ/(3·d_eff), also P = √(3π·d_eff/(Gρ)), nach Newton 3·√(π/(Gρ)): Die Grenze liest nur die eigene Dichte."""
    return math.sqrt(3 * math.pi * d_eff / (G * rho))


def gezeit_flach(v, R):
    """Gezeit 𝒯 [1/s²] einer flachen Rotationskurve v [m/s] im Abstand R [m]: −dg/dR = v²/R² = Ω²."""
    return (v / R)**2


def gezeit_quelle(M, D, regel=True, a0_=a0):
    """Gezeit eines Bereichs der Masse M [kg] im Abstand D [m] nach dem Zustand seiner Quelle. Ist sein Gas abgekühlt und
    seine Mode übergeben, zieht seine Regel, g = ν(g_N)·g_N, und der Gezeitentensor ist −∂g/∂x = (g/D)·(d_eff·ŝŝ − 1) mit
    der Flussdimension d_eff = 3 − χ (d_eff − 1 = −d ln g/d ln D). Der isotrope Teil übt kein Moment aus, kippen kann nur
    der gerichtete Teil d_eff·g/D, nach Newton 3GM/D³, im Verhältnis ν·(3 − χ)/3. regel=False: Newton. Gibt
    (d_eff·g/D [1/s²], Verhältnis zu Newton)."""
    gN = G * M / D**2
    if not regel:
        return 3 * gN / D, 1.0
    nu, chi = verstaerkung(gN, a0_), schalter(gN, a0_)
    return (3 - chi) * nu * gN / D, nu * (3 - chi) / 3


def kippung_satellit(M, x, v, R, M_wirt, regel=True, a0_=a0):
    """Kippmoment eines umlaufenden Satelliten der Masse M [kg] (Ort x [m], Geschwindigkeit v [m/s], z längs des
    Drehimpulses der Scheibe) auf den Ring R [m] einer Scheibe, deren Baryonen M_wirt [kg] er umläuft. Der Ring mit Umlauf
    Ω (Regel um die Punktmasse) und Achse n dreht sich unter dem Gezeitentensor A der Quelle mit ṅ = (A·n) × n/(2Ω),
    Moment (mR²/2)·(A·n) × n auf den Drehimpuls mR²Ω. Die ungeneigte Scheibe kippt er mit der Rate k·s_z·s_∥ (s_z und s_∥
    die Anteile seiner Richtung senkrecht zur Scheibe und in ihr), k = d_eff·g/(2ΩD). Eine schon geneigte Scheibe
    (Kipprichtung ε̂) dreht er zusätzlich um k·(s_z² − (ŝ_∥·ε̂)²), also höchstens um k. Sein Azimut dreht mit
    φ̇ = (x·v_y − y·v_x)/(x² + y²), rechtläufig für φ̇ > 0. Er läuft im verstärkten Feld des Wirts und der Wirt in seinem.
    Gibt dict(k [1/s], phi_punkt [1/s], neigung [°] der Bahn gegen die Scheibe, peri, apo [m])."""
    D = math.sqrt(sum(a * a for a in x))
    Om = bahntempo_regel(M_wirt, R, a0_) / R
    k = gezeit_quelle(M, D, regel, a0_)[0] / (2 * Om)
    L = (x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0])
    v_r = sum(a * b for a, b in zip(x, v)) / D
    v_t = math.sqrt(max(sum(b * b for b in v) - v_r**2, 0.0))
    peri, apo = bahn_wendepunkte(M_wirt, D, v_r, v_t, a0_, M2=M)
    neigung = math.degrees(math.acos(L[2] / math.sqrt(sum(a * a for a in L))))
    return dict(k=k, phi_punkt=L[2] / (x[0]**2 + x[1]**2), neigung=neigung, peri=peri, apo=apo)


def feld_quelle(M_B, M_h, r, a0_=a0):
    """Feld [m/s²] einer Quelle nach dem Zustand ihrer Quelle im Abstand r [m]: Newton aus den heutigen Baryonen M_h [kg],
    die Rückkehr aus den einmal gebundenen M_B [kg], g = g_N + (ν(g_N,B) − 1)·g_N,B (auch für Felder von r)."""
    import numpy as np
    gB = G * M_B / r**2
    return G * M_h / r**2 + gB / -np.expm1(-np.sqrt(gB / a0_)) - gB


def bahnen_gemeinsam(M, koerper, dauer, schritte=3001, H=None, a0_=a0, wirt=False, eltern=None):
    """Wirt mit den Baryonen M [kg], heute ruhend im Ursprung, und Körper (M_B, M_heute, x, v) [kg, kg, m, m/s] gemeinsam
    rückwärts in der Zeit, dauer [s]. Jeder Körper zieht nach dem Zustand seiner Quelle (feld_quelle), der Wirt mit
    ν(g_N)·g_N. Ein Bereich läuft als Ganzes im verstärkten Feld seines Wirts und aller, die nicht in ihm liegen; von
    seinen eigenen Satelliten erfährt er den Gegenzug zu ihrem Fall in seinem verstärkten Feld, m_Sat/m_Wirt mal so groß,
    der Impuls des Bereichs bleibt.
    eltern: je Körper der Index seines Wirts (0 der Wirt im Ursprung, k der k-te Körper, −1 neben dem Wirt, Vorgabe 0);
    wer neben dem Wirt läuft, spürt dessen Bereich als einen, mit allen Baryonen darin am Schwerpunkt. Körper mit M_B = 0
    sind Probekörper. Dazu die Rückstellung H²·(x − x_Wirt) (Vorgabe H_∞). Gibt (Zeiten [s], Orte und Geschwindigkeiten
    gegen den Wirt [m, m/s] je (n, 3, Zeit), ort(t): Orte gegen den Wirt [m] zu jeder Zeit t [s] als (n, 3)), mit wirt
    dazu Ort und Geschwindigkeit des Wirts im ruhenden Rahmen [m, m/s] je (3, Zeit)."""
    import numpy as np
    from scipy.integrate import solve_ivp
    H = groesste_umkehr(M, a0_)[1] if H is None else H
    q = [(M, M)] + [(k[0], k[1]) for k in koerper]
    n = len(q)
    vorfahren = [set()]
    for e in (eltern or [0] * len(koerper)):
        vorfahren.append(set() if e < 0 else {e} | vorfahren[e])
    innen = np.array([i == 0 or 0 in vorfahren[i] for i in range(n)])
    neben = ~innen
    mB_innen, mh_innen = sum(q[i][0] for i in range(n) if innen[i]), np.array([q[i][1] if innen[i] else 0.0 for i in range(n)])

    def rate(t, y):
        x = y[:3 * n].reshape(n, 3)
        a = np.zeros((n, 3))
        for j, (mB, mh) in enumerate(q):
            if mB <= 0:
                continue
            d = x[j] - x
            r = np.linalg.norm(d, axis=1)
            r[j] = 1.0
            g = feld_quelle(mB, mh, r, a0_)
            g[j] = 0.0
            for i in vorfahren[j]:                               # der Gegenzug im Bereich
                g[i] = mh / q[i][1] * feld_quelle(q[i][0], q[i][1], r[i:i + 1], a0_)[0]
            if innen[j]:
                g[neben] = 0.0                                   # von außen ein Bereich, unten
            a += (g / r)[:, None] * d
        if neben.any():
            s_ = mh_innen @ x / mh_innen.sum()
            d = s_ - x[neben]
            r = np.linalg.norm(d, axis=1)
            a[neben] += (feld_quelle(mB_innen, mh_innen.sum(), r, a0_) / r)[:, None] * d
        a[1:] += H**2 * (x[1:] - x[0])
        return np.concatenate([y[3 * n:], a.ravel()])
    y0 = np.concatenate([np.zeros(3)] + [np.asarray(k[2], float) for k in koerper]
                        + [np.zeros(3)] + [np.asarray(k[3], float) for k in koerper])
    t = np.linspace(0, -dauer, schritte)
    lsg = solve_ivp(rate, (0, -dauer), y0, t_eval=t, rtol=1e-9, atol=1e-3, method="DOP853", dense_output=True)
    x = lsg.y[:3 * n].reshape(n, 3, -1)
    v = lsg.y[3 * n:].reshape(n, 3, -1)
    ort = lambda tt: (lambda z: z[3:3 * n].reshape(n - 1, 3) - z[:3])(lsg.sol(tt))
    if wirt:
        return t, x[1:] - x[0], v[1:] - v[0], ort, x[0], v[0]
    return t, x[1:] - x[0], v[1:] - v[0], ort


def warp_ringe(quellen, M, R_d=None, h=0.5 * kpc, start=None, gekoppelt=True, regel=True, daempfung=0.0, R_aus=14 * kpc,
               R_innen=8 * kpc, N=70, R_rand=30 * kpc, eps=kpc, n_phi=180, a0_=a0, gemeinsam=False, eltern=None,
               kopplung_regel=True):
    """Biegung einer Scheibe aus gekoppelten Ringen unter umlaufenden Satelliten. quellen: je Satellit (M_B, x, v) oder
    (M_B, x, v, M_heute), M_B [kg] die Baryonen, die er einmal gebunden hat, M_heute die heutigen (Vorgabe M_B), Ort x [m]
    und Geschwindigkeit v [m/s] heute, z längs des Drehimpulses der Scheibe. Die Scheibe: Baryonen M [kg], exponentiell mit
    R_d (Vorgabe an der Schwelle der Flächendichte, r*/√2), N Ringe bis R_rand, Dicke h. Die Regel antwortet quer mit ν wie
    der Umlauf: Ring i kippt gegen Ring j mit B_ij = ν_i·G·m_j·R_j/(4π·R_i·Ω_i)·∫cos φ dφ/D³ (Laplace und Lagrange mit ν),
    ẇ_i = i·Σ_j B_ij·(w_j − w_i), w der Kippvektor als komplexe Zahl. Das Feld jedes Satelliten, Newton aus den heutigen
    Baryonen, liest der Ring mit seiner eigenen Regel quer zur Neigung, ν_i·g_N, direkt summiert über n_phi Punkte jedes
    flachen Rings (Glättung eps): Moment τ_i, Antrieb (τ_x + iτ_y)/(m_i·R_i²·Ω_i). Seine Bahn läuft rückwärts im
    verstärkten Feld der Scheibe als Punktmasse, die Scheibe erfährt den Gegenzug, M_heute/M mal so groß, dazu die
    Rückstellung H_∞²·r. Start mit flacher Scheibe zur Zeit start [s] vor heute (Vorgabe: letztes Apozentrum des ersten
    Satelliten). gekoppelt=False: Ringe ohne Kopplung. regel=False: die Satelliten ziehen nach Newton, die Bahnen bleiben.
    daempfung [1/s]: Abklingen der freien Biegungen (Vorgabe 0). kopplung_regel=False: die Ringe koppeln nur mit Newton, der
    Anteil der Rückkehr wirkt rund und gibt keine eigene Präzession, der Umlauf bleibt mit der Regel. gemeinsam: die Bahnen aller Satelliten und der Scheibe
    gemeinsam (bahnen_gemeinsam mit eltern), jeder auch im Feld der anderen, 5 Gyr zurück. Gibt (Neigung bei R_aus gegen die innere Scheibe bis
    R_innen, nach Drehimpuls gemittelt [rad], Drehung ihrer Knotenlinie heute [1/s], rechtläufig positiv)."""
    import numpy as np
    from scipy.integrate import solve_ivp
    R_d = schwellenradius(M, a0_) / math.sqrt(2) if R_d is None else R_d
    H = groesste_umkehr(M, a0_)[1]

    def zug(m_B, m, r, regel=True):                              # Newton aus m, Rückkehr aus m_B
        return feld_quelle(m_B, m, r, a0_) if regel else G * m / r**2
    quellen = [(q[0], q[1], q[2], q[3] if len(q) > 3 else q[0]) for q in quellen]
    bahnen = []
    if gemeinsam:
        ort = bahnen_gemeinsam(M, [(q[0], q[3], q[1], q[2]) for q in quellen], 5 * Gyr, a0_=a0_, eltern=eltern)[3]
        bahnen = [type("Bahn", (), {"sol": staticmethod(lambda t, i=i: np.concatenate([ort(t)[i], np.zeros(3)]))})()
                  for i in range(len(quellen))]
    for M_s, x, v, M_h in (() if gemeinsam else quellen):
        def bahn(t, y, M_s=M_s, M_h=M_h):
            r = math.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
            f = -zug(M, M, r) * (1 + M_h / M) / r + H**2
            return [y[3], y[4], y[5], f * y[0], f * y[1], f * y[2]]
        bahnen.append(solve_ivp(bahn, (0, -20 * Gyr), list(x) + list(v), rtol=1e-11, atol=1e-3, max_step=2e6 * yr,
                                dense_output=True))
    if start is None:                                            # letztes Apozentrum, rückwärts das erste Maximum von r
        tt = np.linspace(0, -(5 if gemeinsam else 20) * Gyr, 20001)
        r = np.array([np.linalg.norm(bahnen[0].sol(z)[:3]) for z in tt]) if gemeinsam else \
            np.linalg.norm(bahnen[0].sol(tt)[:3], axis=0)
        apo = np.nonzero((np.diff(r)[:-1] > 0) & (np.diff(r)[1:] <= 0))[0]
        start = -tt[apo[0] + 1] if len(apo) else 8 * Gyr
    Re = np.linspace(0, R_rand, N + 1)
    R = 0.5 * (Re[1:] + Re[:-1])
    F = lambda a: (1 + a / R_d) * np.exp(-a / R_d)
    m = M * (F(Re[:-1]) - F(Re[1:]))
    phi = np.linspace(0, 2 * math.pi, 721)[:-1]
    Dm = np.sqrt(R[:, None, None]**2 + R[None, :, None]**2 - 2 * R[:, None, None] * R[None, :, None] * np.cos(phi) + h**2)
    w_phi = 2 * math.pi / len(phi)
    gR = np.sum(G * m / (2 * math.pi) * np.sum((R[:, None, None] - R[None, :, None] * np.cos(phi)) / Dm**3, axis=2) * w_phi,
                axis=1)
    nu = np.array([verstaerkung(g, a0_) for g in gR])
    Om = np.sqrt(nu * gR / R)
    B = (nu if kopplung_regel else np.ones(N))[:, None] * G * m[None, :] * R[None, :] / (4 * math.pi * R[:, None] * Om[:, None]) * np.sum(np.cos(phi) / Dm**3,
                                                                                                         axis=2) * w_phi
    np.fill_diagonal(B, 0.0)
    if not gekoppelt:
        B[:] = 0.0
    Bs = B.sum(axis=1)
    ps = (np.arange(n_phi) + 0.5) * 2 * math.pi / n_phi
    P = np.stack([R[:, None] * np.cos(ps), R[:, None] * np.sin(ps), np.zeros((N, n_phi))], 2)

    def antrieb(t):
        f = np.zeros(N, complex)
        for (M_s, _, _, M_h), b in zip(quellen, bahnen):
            d = b.sol(t)[:3] - P
            dd = np.sqrt(np.sum(d * d, axis=2) + eps**2)
            tau = np.mean(np.cross(P, ((nu[:, None] if regel else 1.0) * G * M_h / dd**3)[:, :, None] * d), axis=1)
            f += (tau[:, 0] + 1j * tau[:, 1]) / (R**2 * Om)
        return f

    def ableitung(t, y):
        w = y[:N] + 1j * y[N:]
        dw = 1j * (B @ w - Bs * w) - daempfung * w + antrieb(t)
        return np.concatenate([dw.real, dw.imag])
    y = solve_ivp(ableitung, (-start, 0), np.zeros(2 * N), rtol=1e-9, atol=1e-14, max_step=2e6 * yr).y[:, -1]
    dy = ableitung(0.0, y)
    w, dw = y[:N] + 1j * y[N:], dy[:N] + 1j * dy[N:]
    innen, L = R < R_innen, m * R**2 * Om
    i = int(np.argmin(abs(R - R_aus)))
    rel = w[i] - np.sum(L[innen] * w[innen]) / np.sum(L[innen])
    drel = dw[i] - np.sum(L[innen] * dw[innen]) / np.sum(L[innen])
    return abs(rel), float((drel / rel).imag)


def reibung(m, v, rho_b, sigma, b_min, b_max, aN=None, ell=None, a0_=a0):
    """Bremsung |dv/dt| [m/s²] eines Körpers der Masse m [kg], der mit v [m/s] auf einer Kreisbahn durch Sterne der Dichte
    ρ_b [kg/m³] mit isotroper Streuung σ [m/s] läuft (Chandrasekhar), b_min und b_max [m] die Grenzen des Logarithmus.
    aN [m/s²]: Feld, das die Rückkehr eines Bereichs mit übergebener Mode liest. Das Kielwasser einer Welle k entsteht mit
    G·ν[1 − (χ/2)(k̂·ĝ)²] und zieht quer zur Neigung mit G·ν. Gemittelt über die Wellen mit dem Gewicht der Bremsung
    μ²·e^{−X²μ²}, μ = k̂·v̂, gibt das 1 − (χ/4)(1 − μ̄²), μ̄² von 3/5 (langsam) bis 0 (schnell). Bis zur Glättungslänge ell
    bleibt der Körper Einzeltrichter (Newton):
    |dv/dt| = 4πG²m·ρ_b·F(X)·[ln(ℓ/b_min) + ν²(1 − (χ/4)(1 − μ̄²))·ln(b_max/ℓ)]/v², F(X) = erf X − (2X/√π)·e^{−X²},
    X = v/(√2σ). aN = None: Newton, ln(b_max/b_min), so rechnet auch der Beobachter mit Halo, mit seiner Dichte statt ρ_b."""
    X = v / (math.sqrt(2) * sigma)
    F = math.erf(X) - 2 * X * math.exp(-X * X) / math.sqrt(math.pi)
    lam = math.log(b_max / b_min)
    if aN is not None:
        J2 = math.sqrt(math.pi) * math.erf(X) / (4 * X**3) - math.exp(-X * X) / (2 * X * X)   # ∫₀¹ μ²·e^{−X²μ²} dμ
        mu2 = 3 / (2 * X * X) - math.exp(-X * X) / (2 * X * X * J2)
        ell = b_min if ell is None else min(max(ell, b_min), b_max)
        lam = math.log(ell / b_min) + verstaerkung(aN, a0_)**2 * (1 - schalter(aN, a0_) / 4 * (1 - mu2)) * math.log(b_max / ell)
    return 4 * math.pi * G**2 * m * rho_b * F * lam / v**2


def sinkzeit(m, r0, M, a, sigma, r_koerper, regel=True, r_ende=None, a0_=a0, schritte=4000):
    """Zeit [s], in der ein Körper der Masse m [kg] auf einer Kreisbahn durch Reibung (reibung) von r0 bis r_ende [m]
    (Vorgabe a/10) sinkt, in einer Plummer-Kugel der Baryonen M [kg] mit der Länge a [m] und der Streuung sigma [m/s]:
    ρ_b = 3M/(4πa³)·(1 + r²/a²)^{−5/2}, a_N = GMr/(r² + a²)^{3/2}, mit der Regel a = ν(a_N)·a_N (regel), sonst Newton
    mit den Baryonen allein. b_max = r, b_min der größere von r_koerper [m] und Gm/v², Glättungslänge ℓ³ = Gm/𝒯 mit der
    Gezeit 𝒯 = −da/dr der Kugel. Der Drehimpuls L = r·v fällt mit dL/dt = −r·|dv/dt|, also dr/dt = −r·|dv/dt|/(dL/dr)."""
    import numpy as np
    r_ende = a / 10 if r_ende is None else r_ende
    aN = lambda r: G * M * r / (r * r + a * a)**1.5
    acc = (lambda r: verstaerkung(aN(r), a0_) * aN(r)) if regel else aN
    rho = lambda r: 3 * M / (4 * math.pi * a**3) * (1 + r * r / (a * a))**-2.5
    rs = np.geomspace(r_ende, r0, schritte)
    t = 0.0
    for r1, r2 in zip(rs[1:][::-1], rs[:-1][::-1]):
        r = math.sqrt(r1 * r2)
        h = 1e-4 * r
        v = math.sqrt(acc(r) * r)
        dL = (math.sqrt(acc(r + h) * (r + h)) * (r + h) - math.sqrt(acc(r - h) * (r - h)) * (r - h)) / (2 * h)
        gez = max(-(acc(r + h) - acc(r - h)) / (2 * h), 1e-60)
        ell = (G * m / gez)**(1 / 3)
        b_min = max(r_koerper, G * m / v**2)
        f = reibung(m, v, rho(r), sigma, b_min, max(r, 1.0001 * b_min), aN=aN(r) if regel else None, ell=ell, a0_=a0_)
        t += (r1 - r2) * dL / (r * f)
    return t


def relaxation_regel(aN, a0_=a0):
    """Relaxation durch Begegnungen in einem kollektiven Bereich mit übergebener Mode, ν und χ im Feld a_N [m/s²], das die
    Rückkehr liest (die einmal gebundenen Baryonen). Jede Begegnung im Abstand b gibt den Stoß δv = 2Gm/(bv). Bis zur
    Glättungslänge ℓ, ℓ³ = Gm/𝒯, ist der stoßende Stern Einzeltrichter (Newton), jenseits antwortet die Regel auf sein
    Feld, δa = ν[δg − (χ/2)ĝ(ĝ·δg)], im Kugelmittel (eine Richtung längs, zwei quer) |δa|² = ν²·f·|δg|² mit
    f = (2 + (1 − χ/2)²)/3, tief ¾. In der Zeit R/v ist Δv²/v² = (8/N)·(GM/(v²R))²·[ln(ℓ/b_min) + ν²f·ln(R/ℓ)],
    b_min = Gm/v². Tief ist 𝒯 = Ω² = v²/R², also ℓ³ = b_min·R²: ℓ teilt ln(R/b_min) im Drittel. Ohne Gedächtnis
    (v² = ν·GM/R) ist die Klammer je ln(R/b_min) dann 2/(3ν²) + f/3, tief ¼ (Newton 1). Der Beobachter mit Halo rechnet
    bei demselben Tempo nur mit den Stößen der Sterne nach Newton, (8/N)·(GM/(v²R))²·ln(R/b_min): Der Bereich relaxiert
    (2 + ν²f)/3-mal so schnell, und ohne Gedächtnis bestünde sein Halo als Stoff gelesen aus Körpern von
    (ν²f − 1)/(3(ν − 1)) Sternmassen. Gibt (f, Klammer je ln(R/b_min) ohne Gedächtnis, Verhältnis zum Beobachter mit
    Halo, Masse der Körper seines Halos in Sternmassen)."""
    nu, chi = verstaerkung(aN, a0_), schalter(aN, a0_)
    f = (2 + (1 - chi / 2)**2) / 3
    koerper = (nu**2 * f - 1) / (3 * (nu - 1)) if nu > 1 else 0.0
    return f, 2 / (3 * nu**2) + f / 3, (2 + nu**2 * f) / 3, koerper


BREMSSTRAHLUNG = 1.4e-40 * 1.2                  # W·m³·K^(−1/2), Λ = Wert·√T für vollständig ionisiertes Gas, Gaunt-Faktor 1,2 (Rybicki und Lightman 1979)


WEITE_DOPPELSTERNE = dict(abstand_pc=130, kAU=(1, 30), mond_zuwachs=0.20)   # Gaia DR3: kein Zuwachs, Newton (Cookson u. a. 2026)


def streuung_kugel(M, a0_=a0):
    """Druckgestützte Kugel ganz im Tiefbereich, isotrop, jedes Profil: σ_v⁴ = (4/81)·GMa₀."""
    return (4 / 81 * G * M * a0_)**0.25


def isotherme_kugel(rho0=None, r_ende=1e8):
    """Selbsttragende isotherme Kugel, isotrop, in Einheiten G = a₀ = σ = 1: d ln ρ/dr = −g/σ², dM/dr = 4πr²ρ. rho0 = None:
    ganz im Tiefbereich, g = √M/r. Sonst mit der vollen Regel g = ν(a_N)·a_N ab der Dichte rho0 in der Mitte. Gibt
    (γ = −d ln ρ/d ln r bei r_ende, M·4/81): tief γ = 9/2 und M = (81/4)·σ⁴/(G·a₀), mit Newtonschem Kern γ kleiner."""
    from scipy.integrate import solve_ivp
    if rho0 is None:
        g = lambda r, M: math.sqrt(max(M, 0.0)) / r
        y0, r0 = [0.0, 4 * math.pi / 3 * 1e-18], 1e-6
    else:
        g = lambda r, M: verstaerkung(M / r**2, 1.0) * M / r**2 if M > 0 else 0.0
        r0 = 1e-6 / math.sqrt(rho0)
        y0 = [math.log(rho0), 4 * math.pi / 3 * rho0 * r0**3]
    f = lambda r, y: [-g(r, y[1]), 4 * math.pi * r * r * math.exp(y[0])]
    s = solve_ivp(f, (r0, 10 * r_ende), y0, rtol=1e-11, atol=1e-30, dense_output=True, method="DOP853")
    return -(s.sol(1.01 * r_ende)[0] - s.sol(r_ende)[0]) / math.log(1.01), s.sol(r_ende)[1] * 4 / 81


SONNE_GALAKTISCH = dict(R0=8.249,       # kpc, Abstand zum Zentrum (GRAVITY, wie bei S2)
                        z=0.0208,       # kpc über der Scheibe (Bennett und Bovy 2019)
                        v=(12.9, 245.6, 7.78))  # km/s gegen das Zentrum (Drimmel und Poggio 2018, nach Reid und Brunthaler 2004)


def bahn_wendepunkte(M, r, v_r, v_t, a0_=a0, M2=0.0):
    """Peri- und Apozentrum [m] einer Bahn im Feld der Punktmasse M nach der Regel, g = ν(GM/r²)·GM/r², aus Ort r [m],
    radialer und tangentialer Geschwindigkeit [m/s]: Nullstellen von ½v_r² + L²/2r² − L²/2x² − (Φ(x) − Φ(r)). M2: die
    Baryonen des Körpers selbst, wenn sein Gas abkühlte. Dann läuft er im verstärkten Feld von M und M in seinem, die
    Felder beider addieren sich (Vorgabe 0)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    g = lambda x: sum(verstaerkung(G * m / x**2, a0_) * G * m / x**2 for m in (M, M2) if m > 0)
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


def bahn_phase(M, r, v_r, v_t, a0_=a0):
    """Phase einer Bahn im Feld der Punktmasse M nach der Regel: die Laufzeit vom Perizentrum bis r durch die vom Peri-
    bis zum Apozentrum, ∫dx/|v_r(x)|, 0 im Perizentrum, 1 im Apozentrum (inf, wenn ungebunden). Eine Schar im
    Gleichgewicht verteilt sich gleich über die Phasen, Mittel ½ und Streuung 1/√12."""
    from scipy.integrate import quad
    peri, apo = bahn_wendepunkte(M, r, v_r, v_t, a0_)
    if not math.isfinite(apo):
        return math.inf
    if apo - peri < 1e-9 * r:
        return 0.5
    g = lambda x: verstaerkung(G * M / x**2, a0_) * G * M / x**2
    L, E = r * v_t, 0.5 * (v_r**2 + v_t**2)
    phi = lambda x: quad(g, r, x, limit=200)[0]
    # x = peri + (apo − peri)·sin²u hebt die Wurzeln an beiden Wendepunkten heraus
    def zeit(u1):
        def h(u):
            x = peri + (apo - peri) * math.sin(u)**2
            w = 2 * (E - phi(x)) - L**2 / x**2
            return 2 * (apo - peri) * math.sin(u) * math.cos(u) / math.sqrt(w) if w > 0 else 0.0
        return quad(h, 0, u1, limit=200)[0]
    u_r = math.asin(math.sqrt(min(max((r - peri) / (apo - peri), 0.0), 1.0)))
    return zeit(u_r) / zeit(math.pi / 2)


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


def streuung_im_feld(M, R_e, g_e=0.0, a0_=a0, n=200, M_B=None, einfall=None):
    """Sichtlinienstreuung σ [m/s] einer druckgestützten Plummer-Kugel (Masse M, projizierter Halblichtradius R_e) im
    äußeren Newtonfeld g_e der Umgebung, isotrop. Die Regel wirkt auf das Gesamtfeld, g = ν(|g_N,i + g_e|)·(g_N,i + g_e);
    die Bewegung als Ganzes fällt im Virialsatz heraus (∫ρ·r = 0): σ² = −∫ρ r·g dV/(3M), integriert über den Massenanteil x
    und den Winkel μ zum äußeren Feld. Ohne äußeres Feld ganz im Tiefbereich σ⁴ = (4/81)·GMa₀, ganz Newtonsch
    σ² = πGM/(32R_e). M_B: die Baryonen, die der Bereich einmal gebunden hat (Gedächtnis der Bindung, gleiches Profil); die
    Neigung am Ort liest M, die Rückkehr (innere Richtung) M_B: g = g_N + (ν(|g_N,B|) − 1)·g_N,B. einfall = (c, Q, H₀ [1/s]):
    statt der Regel die nicht übergebene Mode Q·M in der Form nach dem Einfall (NFW mit Konzentration c, M₂₀₀ = (1 + Q)·M bei
    200·ρ_krit), ohne äußeres Feld."""
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
    if einfall is not None:
        c_, Q_, H_ = einfall
        R200 = (3 * (1 + Q_) * M / (4 * math.pi * 200 * 3 * H_**2 / (8 * math.pi * G)))**(1 / 3)
        m = lambda y: np.log(1 + y) - y / (1 + y)
        g = -g_i - G * Q_ * M * m(r * c_ / R200) / m(c_) / r**2 + 0 * mu[None, :]
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


def atom_und_schliessung(a0_=a0):
    """Am Atom die Windung α = r_e/λ̄_e und die Faltung von Elektron und Proton α_G,ep = 3Δs_p/λ̄_e, ihr Verhältnis
    r_e/(3Δs_p), gegen den Abstand des Elektrons von der Schließung ℓ₀/λ̄_e. Gleich hieße r*_p² = λ̄_e·r_e: die Schwelle des
    Protons im Mittel der beiden Längen der Windung des Elektrons, dann a₀ = G·m_p/(α³·a_B²) und
    N_b = (3/(2(1 + Q)))·α·(m_Pl²/(m_e·m_p))². Gibt (ℓ₀/λ̄_e, α/α_G,ep, r*_p/√(λ̄_e·r_e), a₀ aus dem Atom, N_b aus α)."""
    r_e = alpha * compton(m_e)
    x = c**2 / a0_ / compton(m_e)
    y = r_e / (3 * faltung(m_p))
    Q = q_je_radiant()
    return (x, y, schwellenradius(m_p, a0_) / math.sqrt(compton(m_e) * r_e), G * m_p / (alpha**3 * bohr_radius()**2),
            3 / (2 * (1 + Q)) * alpha * (m_Pl**2 / (m_e * m_p))**2)


def exponent_inhalt(f, groesse, eps=1e-3):
    """Wie eine Größe am Inhalt hängt: d ln f/d ln x für x = α ("alpha"), m_e ("m_e") oder m_p ("m_p"), zentral mit
    x·(1 ± eps) gerechnet, danach wieder mit den gemessenen Werten. α wird bei festem e, ħ und c mit ε₀ ∝ 1/α und μ₀ ∝ α
    zusammen geändert, so folgt jede Formel mit, die α über die Ladung liest (Coulomb, μ₀). f: Funktion ohne Argument mit
    positivem Wert. Null heißt, die Rechnung liest x nicht, jede andere Zahl, wie x in ihr steht (a_B: −1 in α und m_e,
    M_*: −2 in m_p). Nicht erfasst: Vorgaben von Argumenten, die beim Laden gebunden sind (verlust_windung, M = m_e), die
    gemessenen m_H und m_u und gespeicherte Zwischenwerte (lru_cache): f darf aus dem Speicher nur Daten lesen, keine
    Rechnung, die x enthält."""
    g = globals()
    namen = {"alpha": (("alpha", 1), ("eps0", -1), ("mu0", 1)), "m_e": (("m_e", 1),), "m_p": (("m_p", 1),)}[groesse]
    alt = {n: g[n] for n, _ in namen}
    werte = []
    try:
        for s in (1 + eps, 1 - eps):
            for n, k in namen:
                g[n] = alt[n] * s**k
            werte.append(f())
    finally:
        g.update(alt)
    return (math.log(werte[0]) - math.log(werte[1])) / (2 * eps)


def schwelle_spiegel(M, a0_=a0):
    """r̃ = √(λ̄·ℓ₀): Mittel der tangentialen und der transversalen Länge, die Schwelle der Spiegelmasse, r̃ = √2·r*(m̃)."""
    return math.sqrt(compton(M) * c**2 / a0_)


def kreuzung_gehalten(halter, a0_=a0, traeger=None):
    """Der gehaltene Körper R = λ_H·N^{1/3}, N = M/Träger (Vorgabe Proton), als Linie in der Ebene aus Masse und Hülle:
    Masse, bei der er r̃ kreuzt, M = m·(λ̄_m·ℓ₀/λ_H²)^{3/5}, und Masse, bei der er r* kreuzt, M = m·(λ_H/r*_m)⁶."""
    m = m_p if traeger is None else traeger
    return m * (compton(m) * c**2 / a0_ / halter**2)**0.6, m * (halter / schwellenradius(m, a0_))**6


def leiter_laenge(k, r, M, a0_=a0):
    """Eine Formel für jede Länge der Leiter: λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1, ein Punkt im Dreieck der drei Richtungen,
    gleichwertig ℓ = ℓ_P·(M/m_Pl)^k·ε_L^{−r} mit k = q − p (parametrisch gegen tangential) und r (transversal).
    λ̄ (−1, 0), 3Δs (1, 0), ℓ₀ (0, 1), ℓ_P (0, 0), r* (½, ½), r̃ (−½, ½), Schwerpunkt (0, ⅓), λ̄/α_G (−3, 0)."""
    return l_P * (M / m_Pl)**k * leiterzahl(a0_)**(-r)


def leiter_ort(M, a0_=a0):
    """Lage einer Masse auf der Leiter, x = ln(M/m_Pl)/ln(1/ε_L): −1 bei M_min, 0 bei m_Pl, 1 bei M_U."""
    return math.log(M / m_Pl) / math.log(1 / leiterzahl(a0_))


def leiter_treffen(punkt1, punkt2, a0_=a0, Z=1.0):
    """Masse, bei der zwei Längen der Leiter, (k₁, r₁) und (k₂, r₂), im Verhältnis ℓ₁/ℓ₂ = Z stehen (Vorgabe Z = 1: gleich).
    In x = ln(M/m_Pl)/ln(1/ε_L) ist jede Länge die Gerade y = r + k·x, sie treffen sich bei x = (r₂ − r₁)/(k₁ − k₂). Die
    Steigung k = 1/n hält GM/R^n fest (Spreizung, Neigung, Gezeit, im Tiefbereich der Umlauf), k → −k tauscht es gegen M·R^n.
    Ein festes Verhältnis Z liegt bei der Masse des Treffpunkts mal Z^{1/(k₁ − k₂)}: 2 zwischen Spreizung und Neigung,
    3/2 zwischen Spreizung und Gezeit, 6 zwischen Neigung und Gezeit, 4 zwischen Neigung und Umlauf im Tiefbereich.
    Gibt (M, x)."""
    (k1, r1), (k2, r2) = punkt1, punkt2
    x = (r2 - r1 + math.log(Z) / math.log(1 / leiterzahl(a0_))) / (k1 - k2)
    return m_Pl * leiterzahl(a0_)**(-x), x


def kopplungen(M, a0_=a0):
    """Die drei Kopplungen eines Trichters, je eine Kante zwischen zwei Richtungen, Verhältnisse seiner Grundlängen:
    α_G = 3Δs/λ̄ (parametrisch zu tangential), α_I = λ̄/ℓ₀ = kT_I/(Mc²) (tangential zu transversal) und 3Δs/ℓ₀
    (parametrisch zu transversal)."""
    l, d, l0 = compton(M), 3 * faltung(M), c**2 / a0_
    return d / l, l / l0, d / l0


def schnittmassen(a0_=a0):
    """M_min = ε_L·m_Pl, M_c = ε_L^{1/3}·m_Pl, m_Pl, M_U = m_Pl/ε_L."""
    e = leiterzahl(a0_)
    return e * m_Pl, e**(1 / 3) * m_Pl, m_Pl, m_Pl / e


def ablenkung_grenz(M, a0_=a0):
    """Ablenkung am kollektiven Trichter für b ≫ r*: α_∞ = 6πΔs/r* = 2π·(v_flach/c)²."""
    return 2 * math.pi * (v_flach(M, a0_) / c)**2



def ablenkung_regel(M, b, a0_=a0):
    """Ablenkung an einer Punktmasse mit der Regel längs des ganzen Strahls, α = (2/c²)·∫ν(a_N)·a_N·b/r dl, beide Hälften wie
    12Δs/b: innen 12Δs/b·(1 + (π²/6)·(b/r*)²), denn ν − 1 = Σqⁿ und die n-te Rückkehr gibt 1/n², außen α_∞ = 6πΔs/r*."""
    from scipy.integrate import quad
    rs = schwellenradius(M, a0_)
    f = lambda th: math.cos(th) * verstaerkung(G * M * math.cos(th)**2 / b**2, a0_) if th < math.pi / 2 else 0.0
    return 4 * G * M / (c**2 * b) * quad(f, 0, math.pi / 2, points=[math.acos(min(1, b / rs))] if b < rs else None, limit=200)[0]

# ---------------------------------------------------------------- 04 Entropie der Engstelle
def entropie_falte(M):
    """S/k = 36π·(Δs/ℓ_P)²."""
    return 36 * math.pi * (faltung(M) / l_P)**2


def waerme_im_takt(M, r, T_aussen):
    """Wärme im Takt (Tolman): Ein ruhender Körper im Gleichgewicht hat eine Wärmelänge, eine Schleife des ganzen Körpers
    in t, und jede ruhende Uhr misst sie mit ihrem Takt, T·dτ/dt = T_∞. Ebenso spürt sie die Neigung, in t die halbe
    Steigung g = ½·d(c²·(dτ/dt)²)/dr, als g/(dτ/dt) = c²·d(dτ/dt)/dr. Im gleichmäßig gehaltenen Körper trägt der Druck
    denselben Takt, T/T_R = 1 + P/(ρc²) (mitte_gehalten). Gibt (T am Ort [K], gespürte Neigung [m/s²])."""
    N = takt(M, r)
    return T_aussen / N, neigung(M, r) / N


def flaeche_strahl(D_E=None, a0_=a0):
    """Die Regel längs des Strahls um eine kompakte Quelle, b ≪ r*: Ihr Zuwachs der Ablenkung ist (π²/3)·b/Λ, linear in b,
    die Ablenkung einer gleichförmigen Fläche. In Flächendichte gelesen ist sie (1/πG)·∫(a − a_N)·d ln r, an der Punktmasse
    a₀/(πG)·∫χ ds = (π²/6)·a₀/(πG), die Flächendichte der Schwelle mal dem Integral des Schalters, gleich für jede Masse.
    Gibt Σ [kg/m²], mit D_E dazu (Σ, κ) mit κ = 4πG·Σ·D_E/c² = (π²/6)·2D_E/Λ."""
    S = math.pi**2 / 6 * flaechendichte_schwelle(a0_)
    return S if D_E is None else (S, 4 * math.pi * G * S * D_E / c**2)


def laufzeit_regel(M, D_E, y, a0_=a0):
    """Punktlinse der Masse M mit der Regel längs des ganzen Strahls (ablenkung_regel), Abstand D_E = D_l·D_ls/D_s, Quelle bei
    y·R_E, R_E = √(2r₀·D_E). Bilder aus b − D_E·α(b) = y·R_E, Laufzeit aus dem Fermat-Potential
    [(b − y·R_E)²/(2D_E) − ∫α db]/c, am Ort der Linse gezählt. Gibt (Δt [s], r₀ der Bilder b₊·|b₋|/(2D_E), y der Bilder
    √(b₊/|b₋|) − √(|b₋|/b₊)). Mit der Fläche der Regel ist Δt = (r₀/c)·F(y der Bilder) mit der Faltung der Linse selbst,
    die Bilder lesen r₀/(1 − κ)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    RE = math.sqrt(2 * grenzradius(M) * D_E)
    al = lambda b: ablenkung_regel(M, b, a0_)
    psi = lambda b: quad(al, RE, b, limit=200)[0] if b > RE else -quad(al, b, RE, limit=200)[0]
    bild = lambda b: b - D_E * math.copysign(al(abs(b)), b) - y * RE
    bp, bm = brentq(bild, 0.5 * RE, 4 * RE), brentq(bild, -4 * RE, -0.05 * RE)
    t = lambda b: ((b - y * RE)**2 / (2 * D_E) - psi(abs(b))) / c
    return t(bm) - t(bp), bp * abs(bm) / (2 * D_E), math.sqrt(bp / abs(bm)) - math.sqrt(abs(bm) / bp)


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


def verdichtung_baryonen(M_b, R):
    """Verdichtung der Baryonen M_b innerhalb des Radius R gegen ihr Mittel heute, M_b/((4π/3)·R³·ρ̄_b), ρ̄_b aus ω_b."""
    return M_b / (4 * math.pi / 3 * R**3 * dichte_heute()[1])


def schallhorizont(Q=None, wb=None):
    """r_d [Mpc] am Ende der Schallwelle des frühen Plasmas; die Mode ξ wirkt dort linear wie Staub, also gilt
    die Standardrechnung, hier als Anpassung an sie (Aubourg u. a. 2015, 0,02 %), ω_ν = 0,0006.
    Q = None: ω_c gemessen; sonst ω_c = Q·ω_b. wb: ω_b statt des gemessenen."""
    wb = omega_b if wb is None else wb
    wc = omega_c if Q is None else Q * wb
    return 55.154 * math.exp(-72.3 * 0.0006**2) / ((wb + wc)**0.25351 * wb**0.12807)


def dichte_vakuumteil(a0_=a0):
    """ρ_s = π³a₀²/(30Gc²), zeitfest (Λ fest)."""
    return math.pi**3 * a0_**2 / (30 * G * c**2)


def zelle_laenge(a=1.0, a0_=a0, rho_m=None):
    """Räumliche Zelle L aus a₀ = 3c⁴/(2Gρ_m L³), mitbewegt: L ∝ a; rho_m: Dichte von Baryonen und Mode heute (sonst
    gemessen)."""
    return (3 * c**4 / (2 * G * (dichte_heute()[2] if rho_m is None else rho_m) * a0_))**(1 / 3) * a


def lichtweg_zelle(a0_=a0, Q=None):
    """Teilchen- und Ereignishorizont gegen die Zelle im Formteil (drei Vakuumteile, Strahlung aus T_CMB, flach): der ganze
    Weg des Lichts c·η_∞ = ∫c·dt/a vom Anfang bis in den Endzustand, ohne Strahlung c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴)
    = 0,70 für jedes a₀ und ω_m. Der Radius der letzten Streufläche wächst bis χ*_∞ = c·(η_∞ − η*), die Kreispaare bis
    cos α_∞ = L/(2χ*_∞), Abbilder über die Diagonalen der Flächen (√2·L) erst bei χ* > L/√2. Q = None: ω_c gemessen, sonst
    ω_c = Q·ω_b. Gibt (Teilchenhorizont/L, Ereignishorizont/L, c·η_∞/L, χ*_∞/L, α_∞ in Grad)."""
    de = lambda a: dichte_dunkle_energie(a0_)
    L = zelle_laenge(a0_=a0_, rho_m=dichte_heute(Q)[2])
    teil = abstand_mitbewegt(de, 1e12, Q)
    ereignis = -abstand_mitbewegt(de, -1 + 1e-15, Q)
    chi = abstand_mitbewegt(de, None, Q) + ereignis
    return teil / L, ereignis / L, (teil + ereignis) / L, chi / L, math.degrees(math.acos(min(1.0, L / (2 * chi))))


def theta_zelle(a0_=a0, strahlung=True, z_stern=1090.0):
    """Θ, bei dem die Würfelzelle die letzte Streufläche erreicht (L = 2χ*), im Formteil (drei Vakuumteile, Q = 2π − 1).
    Aus L³ρ_m = 3c⁴/(2Ga₀) folgt L·H/c = (8π²/(Θ·Ω_m))^{1/3}, und mit drei Vakuumteilen ist Ω_m = 1 − (π²/15)·Θ²: Ohne
    Strahlung hängt L/(2χ*) nur von Θ ab (strahlung=False, ohne a₀). strahlung=True: mit der Strahlung aus T_CMB in der
    Zeit gerechnet, bei a₀. Gibt Θ."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    if not strahlung:
        def F(th):
            om = 1 - math.pi**2 / 15 * th**2
            x = quad(lambda u: 1 / (u * u * math.sqrt(om / u**3 + 1 - om)), 1 / (1 + z_stern), 1.0, limit=200)[0]
            return (8 * math.pi**2 / (th * om))**(1 / 3) / (2 * x) - 1
        return brentq(F, 0.9, 1.15, xtol=1e-10)
    Q = q_je_radiant()
    rho_m = ((1 + Q) * omega_b + 0.06 / 93.14) * rho_krit100
    rho_r, rho_de = dichte_heute()[0], dichte_dunkle_energie(a0_)
    H = lambda a: math.sqrt(8 * math.pi * G / 3 * (rho_m / a**3 + rho_r / a**4 + rho_de))
    L = zelle_laenge(a0_=a0_, rho_m=rho_m)
    chi = lambda a: c * quad(lambda x: 1 / (x * x * H(x)), 1 / (1 + z_stern), a, limit=200)[0]
    a_L = brentq(lambda a: L / (2 * chi(a)) - 1, 0.1, 10.0, xtol=1e-10)
    return 2 * math.pi * a0_ / (c * H(a_L))


def schleifen_ganzes():
    """Die drei Schleifen des Ganzen nach Richtungen: die Zelle L (tangential, Umlauf des Lichts), die Wärmeschleife des
    Horizonts 2πc/H (parametrisch, Atmen) und die innere Schleife ℓ₀ = c²/a₀ (transversal). Im Formteil ohne Strahlung,
    Ω_m = 1 − (π²/15)·Θ², ist Θ = (2πc/H)/ℓ₀, L/ℓ₀ = Θ^{2/3}/(π·Ω_m)^{1/3} und L/(2πc/H) = (π·Θ·Ω_m)^{−1/3}. Gibt Θ, bei dem
    je zwei gleich werden: (Horizont = innen, Zelle = innen, Zelle = Horizont)."""
    from scipy.optimize import brentq
    om = lambda th: 1 - math.pi**2 / 15 * th**2
    return (1.0, brentq(lambda th: th**(2 / 3) / (math.pi * om(th))**(1 / 3) - 1, 0.9, 1.15, xtol=1e-12),
            brentq(lambda th: math.pi * th * om(th) - 1, 0.9, 1.15, xtol=1e-12))


def uhr_theta(th, a0_=a0):
    """Die Uhr Θ im Formteil ohne Strahlung: Ḣ = −(3/2)(H² − H_∞²), also H = H_∞·coth(3H_∞t/2) und Θ = Θ_∞·tanh(3H_∞t/2)
    mit Θ_∞ = √15/π und H_∞·τ_D = 2π/Θ_∞. Gibt (Alter t = (√15/3π²)·τ_D·artanh(Θ/Θ_∞) [s], H·t, Verzögerung
    −ä·a/ȧ² = ½ − (3/2)(Θ/Θ_∞)²)."""
    u = th / (math.sqrt(15) / math.pi)
    return (math.sqrt(15) / (3 * math.pi**2) * rueckkehrzeit(a0_) * math.atanh(u), 2 / 3 * math.atanh(u) / u, 0.5 - 1.5 * u * u)


def zeit_zwischen_theta(th1, th2, a0_=a0):
    """Zeit vom Atmen bei Θ₁ bis Θ₂ im Formteil ohne Strahlung, aus der geschlossenen Form der Uhr Θ [s]."""
    return uhr_theta(th2, a0_)[0] - uhr_theta(th1, a0_)[0]


def zelle_kreis(wb, a0_, wa, Q=None):
    """Kante L der Würfelzelle gegen den Durchmesser der letzten Streufläche, L/(2χ*), im Punkt mit Sättigung (ω_c = Q·ω_b,
    H₀ aus a₀, w₀ = −1 − w_a·(1 − a_s)); ist L < 2χ*, schneidet die Streufläche ihre Abbilder in drei Paaren
    gegenüberliegender Kreise vom Radius α, cos α = L/(2χ*). Gibt (L/(2χ*), α in Grad, 0 ohne Kreise)."""
    import camb
    Q = q_je_radiant() if Q is None else Q
    wm = (1 + Q) * wb + 0.06 / 93.14
    a_s = (wm * rho_krit100 / dichte_dunkle_energie(a0_))**(1 / 3)
    H0 = atmen_heute(a0_, Q=Q, wb=wb, wnu=0.06 / 93.14)[0]
    par = camb.set_params(ombh2=wb, omch2=Q * wb, H0=H0, mnu=0.06, w=-1 - wa * (1 - a_s), wa=wa, dark_energy_model="ppf")
    r = camb.get_background(par)
    x = zelle_laenge(a0_=a0_, rho_m=wm * rho_krit100) / Mpc / (2 * (r.conformal_time(0) - r.tau_maxvis))
    return x, (math.degrees(math.acos(x)) if x < 1 else 0.0)


def kreis_erwartung(punkt, alpha_grad, fwhm_grad, lmax=1000):
    """Übereinstimmung, die ein Kreispaar der Zelle erwartet, S = ⟨T₁T₂⟩/⟨T²⟩. Beide Bilder zeigen denselben Punkt der
    letzten Streufläche, gesehen unter 180° − 2α: Der Monopol (Sachs-Wolfe mit der eigenen Wärme) ist gleich, der
    Dopplerteil v·n̂ kehrt sich um, ⟨D₁D₂⟩ = −cos 2α·⟨D²⟩, der späte Takt läuft auf getrennten Wegen. Die Übergänge von CAMB
    im Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, τ, w₀, w_a) werden je Wellenzahl zerlegt, Δ_ℓ(k) ≈ A(k)·j_ℓ(kχ*) + B(k)·j_ℓ'(kχ*)
    (kleinste Quadrate über ℓ), C^AA, C^BB, der Rest. Mit dem Strahl fwhm_grad gewichtet ist S zwischen
    Σ(2ℓ+1)B_ℓ²(C^AA − cos 2α·C^BB)/ΣB_ℓ²(2ℓ+1)C_ℓ (Rest unkorreliert) und demselben mit dem Rest (korreliert). Gibt
    (untere, obere Schranke, Anteil des Dopplers an C_ℓ bei ℓ = 220 und 400)."""
    import camb
    import numpy as np
    from scipy.special import spherical_jn
    wb, wc, H0, lnAs, ns, tau, w0, wa = punkt
    extra = dict(w=w0, wa=wa, dark_energy_model="ppf") if (w0, wa) != (-1.0, 0.0) else {}
    par = camb.set_params(ombh2=wb, omch2=wc, H0=H0, As=math.exp(lnAs) * 1e-10, ns=ns, tau=tau, mnu=0.06,
                          lmax=lmax + 100, **extra)
    par.Accuracy.lSampleBoost = 50
    r = camb.get_transfer_functions(par)
    t = r.get_cmb_transfer_data("scalar")
    chi = r.conformal_time(0) - r.tau_maxvis
    L, q, D = np.array(t.L), np.array(t.q), t.delta_p_l_k[0]
    w_ = L <= lmax + 100
    L, D = L[w_], D[w_]
    A, B = np.zeros(len(q)), np.zeros(len(q))
    for i, k in enumerate(q):
        if np.any(D[:, i]):
            M = np.stack([spherical_jn(L, k * chi), spherical_jn(L, k * chi, derivative=True)], 1)
            A[i], B[i] = np.linalg.lstsq(M, D[:, i], rcond=None)[0]
    gew = 4 * math.pi * math.exp(lnAs) * 1e-10 * (q / 0.05)**(ns - 1) * np.gradient(np.log(q))
    ls = np.arange(2, lmax + 1)
    stelle = {l: j for j, l in enumerate(L)}
    CT, CA, CB, CAB = (np.zeros(len(ls)) for _ in range(4))
    for n, l in enumerate(ls):
        j0, j1 = spherical_jn(l, q * chi), spherical_jn(l, q * chi, derivative=True)
        CT[n], CA[n] = np.sum(gew * D[stelle[l]]**2), np.sum(gew * (A * j0)**2)
        CB[n], CAB[n] = np.sum(gew * (B * j1)**2), np.sum(gew * A * j0 * B * j1)
    s = math.radians(fwhm_grad) / math.sqrt(8 * math.log(2))
    g = (2 * ls + 1) * np.exp(-ls * (ls + 1) * s * s)
    c2 = math.cos(math.radians(2 * alpha_grad))
    rest = CT - CA - CB - 2 * CAB
    unten = float(np.sum(g * (CA - c2 * CB)) / np.sum(g * CT))
    oben = float(np.sum(g * (CA - c2 * CB + rest)) / np.sum(g * CT))
    return unten, oben, float(CB[218] / CT[218]), float(CB[398] / CT[398])


def himmel_zelle(punkt, L, lmax=150, saat=1, k_faktor=4.0):
    """Himmel der Hintergrundstrahlung in der Würfelzelle der Kante L [Mpc] (mitbewegt), exakt aus den Moden des Gitters
    k = 2π·n/L: je Mode ζ gaußisch mit ⟨|ζ|²⟩ = 2π²·𝒫(k)/(k³L³), a_ℓm = 4π·Σ ζ·i^ℓ·Δ_ℓ(k)·Y*_ℓm(k̂) über |k| ≤ k_faktor·lmax/χ*,
    Δ_ℓ(k) von CAMB im Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, τ, w₀, w_a), Sachs-Wolfe, Doppler und später Takt zusammen. Die
    Achsen der Zelle sind die Achsen der Karte. saat: Zahl oder Folge von Zahlen (je eine Ziehung bei einer Rechnung der
    Übergänge). Gibt (a_ℓm in der Ordnung von healpy, µK, oder ihre Liste, und χ* [Mpc])."""
    import camb
    import numba as nb
    import numpy as np
    wb, wc, H0, lnAs, ns, tau, w0, wa = punkt
    extra = dict(w=w0, wa=wa, dark_energy_model="ppf") if (w0, wa) != (-1.0, 0.0) else {}
    par = camb.set_params(ombh2=wb, omch2=wc, H0=H0, As=math.exp(lnAs) * 1e-10, ns=ns, tau=tau, mnu=0.06,
                          lmax=lmax + 50, **extra)
    par.Accuracy.lSampleBoost = 50
    r = camb.get_transfer_functions(par)
    t = r.get_cmb_transfer_data("scalar")
    chi = r.conformal_time(0) - r.tau_maxvis
    ls = list(t.L)
    tr = np.array([t.delta_p_l_k[0, ls.index(l)] for l in range(lmax + 1) if l >= 2])
    N = int(k_faktor * lmax / chi * L / (2 * math.pi)) + 1
    g = np.arange(-N, N + 1, dtype=np.int32)
    nx, ny, nz = np.meshgrid(g, g, g, indexing="ij")
    n2 = nx * nx + ny * ny + nz * nz
    halb = ((nz > 0) | ((nz == 0) & (ny > 0)) | ((nz == 0) & (ny == 0) & (nx > 0))) & (n2 <= N * N)
    nx, ny, nz, n2 = nx[halb], ny[halb], nz[halb], n2[halb]
    schalen, idx = np.unique(n2, return_inverse=True)
    ks = 2 * math.pi * np.sqrt(schalen) / L
    D = np.zeros((len(ks), lmax + 1))
    for j, l in enumerate(range(2, lmax + 1)):
        D[:, l] = np.interp(ks, t.q, tr[j], left=0.0, right=0.0)
    P = math.exp(lnAs) * 1e-10 * (ks / 0.05)**(ns - 1)
    var = 2 * math.pi**2 * P / (ks**3 * L**3)
    kn = np.sqrt(n2.astype(float))
    x, phi = nz / kn, np.arctan2(ny, nx)

    @nb.njit(parallel=True, cache=False)
    def summe(x, phi, re, im, idx, D, lmax, teile):
        alm = np.zeros((teile, lmax + 1, lmax + 1), np.complex128)
        n = len(x)
        for c in nb.prange(teile):
            P = np.zeros((lmax + 1, lmax + 1))
            for p in range(c * n // teile, (c + 1) * n // teile):
                xx = x[p]
                sx = math.sqrt(max(0.0, 1 - xx * xx))
                P[0, 0] = 1 / math.sqrt(4 * math.pi)
                for m in range(1, lmax + 1):
                    P[m, m] = -math.sqrt((2 * m + 1) / (2 * m)) * sx * P[m - 1, m - 1]
                for m in range(0, lmax):
                    P[m + 1, m] = math.sqrt(2 * m + 3) * xx * P[m, m]
                for m in range(0, lmax + 1):
                    for l in range(m + 2, lmax + 1):
                        a = math.sqrt((4 * l * l - 1) / (l * l - m * m))
                        b = math.sqrt(((l - 1)**2 - m * m) / (4 * (l - 1)**2 - 1))
                        P[l, m] = a * (xx * P[l - 1, m] - b * P[l - 2, m])
                for l in range(2, lmax + 1):
                    d = D[idx[p], l]
                    if d == 0.0:
                        continue
                    f = 2 * re[p] * d if l % 2 == 0 else 2 * im[p] * d          # ζ + (−1)^ℓ ζ*: 2Re oder 2i·Im
                    for m in range(0, l + 1):
                        e = complex(math.cos(m * phi[p]), -math.sin(m * phi[p]))
                        alm[c, l, m] += f * P[l, m] * e
        return alm.sum(axis=0)
    import healpy as hp
    alle = []
    for sa in (saat if isinstance(saat, (list, tuple)) else (saat,)):
        rng = np.random.default_rng(sa)
        z = np.sqrt(var[idx] / 2) * (rng.standard_normal(len(idx)) + 1j * rng.standard_normal(len(idx)))
        A = summe(x, phi, z.real.copy(), z.imag.copy(), idx, D, lmax, nb.get_num_threads())
        out = np.zeros(hp.Alm.getsize(lmax), np.complex128)
        for l in range(2, lmax + 1):
            il = 4 * math.pi * (1j)**l * (1 if l % 2 == 0 else 1j)           # ungerade: der Faktor i aus 2i·Im
            for m in range(l + 1):
                out[hp.Alm.getidx(lmax, l, m)] = il * A[l, m] * T_CMB * 1e6
        alle.append(out)
    return (alle if isinstance(saat, (list, tuple)) else alle[0]), chi


def achse_multipol(alm, l, lmax):
    """Achse eines Multipols: die Richtung n, um die der Drehimpuls des Multipols ℓ am größten ist. Σ m²|a_ℓm(n)|² =
    ⟨a|(n·L)²|a⟩ ist eine quadratische Form in n, die Achse also der Eigenvektor zum größten Eigenwert der Matrix
    Re⟨a|L_i L_j|a⟩ (die Karte ist in der Ebene senkrecht zu n am flachsten). alm in der Ordnung von healpy bis lmax.
    Gibt (Achse als Einheitsvektor, Anteil ⟨m²⟩/ℓ² um sie)."""
    import healpy as hp
    import numpy as np
    m = np.arange(-l, l + 1)
    a = np.array([alm[hp.Alm.getidx(lmax, l, abs(k))] if k >= 0 else (-1)**(-k) * np.conj(alm[hp.Alm.getidx(lmax, l, -k)])
                  for k in m])
    Lz = np.diag(m.astype(complex))
    Lp = np.diag(np.sqrt(l * (l + 1) - m[:-1] * (m[:-1] + 1)).astype(complex), -1)      # L₊|m⟩ = √(…)|m + 1⟩
    Lx, Ly = (Lp + Lp.conj().T) / 2, (Lp - Lp.conj().T) / 2j
    Ls = (Lx, Ly, Lz)
    M = np.array([[np.real(a.conj() @ (Ls[i] @ Ls[j] + Ls[j] @ Ls[i]) @ a) / 2 for j in range(3)] for i in range(3)])
    w, v = np.linalg.eigh(M)
    return v[:, -1], float(w[-1] / np.real(a.conj() @ a)) / l**2


def s_halb(cl):
    """Korrelation über große Winkel: S₁/₂ = ∫ C(θ)² d cos θ von cos θ = −1 bis ½, C(θ) = Σ (2ℓ + 1)/(4π)·C_ℓ·P_ℓ(cos θ),
    ℓ ab 2. cl: C_ℓ ab ℓ = 0 [µK²]. Gibt S₁/₂ [µK⁴]."""
    import numpy as np
    from numpy.polynomial import legendre as Lg
    x = np.linspace(-1.0, 0.5, 3001)
    c = np.array(cl, dtype=float).copy()
    c[:2] = 0.0
    C = Lg.legval(x, (2 * np.arange(len(c)) + 1) / (4 * math.pi) * c)
    return float(np.trapezoid(C**2, x))


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


def einsteinring(z_l, z_s=None, dichte_de=None):
    """Ring einer runden Linse: innen die mittlere Flächendichte c²/(4πG·D_E), D_E = D_l·D_ls/D_s (Winkelabstände, flach),
    also am Ring für die eingeschlossene Masse im Bild G·M/R_E² = πG·Σ_kr = c²/(4D_E), für jede Masse und runde Verteilung; Vorgabe Quelle die letzte Streufläche,
    Dunkle Energie der Formteil 3ρ_s. Gibt (D, g)."""
    de = (lambda a: 3 * dichte_vakuumteil()) if dichte_de is None else dichte_de
    Dl = abstand_mitbewegt(de, z_l)
    Ds = abstand_mitbewegt(de, z_s)
    D = Dl / (1 + z_l) * (1 - Dl / Ds)
    return D, c**2 / (4 * D)


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


def spitzenfaktor(q, n=60):
    """Faktor ∏(1 − qⁿ)²⁴ der Spitzenform Δ(q) = q·∏(1 − qⁿ)²⁴ zum Kern-Nome q = Δs/r (mit Vorzeichen, die halbe Drehung
    ist q → −q). Mit dem Kern-Nome ist die Regel a = a_N + χ·Δ(q)·3c²/r* = a_N·(1 + n_B(s)·∏(1 − qⁿ)²⁴), denn r* = s·r und
    χ/s = n_B: Sie weicht von a_N/(1 − e^{−s}) nur um diesen Faktor ab, ≈ 1 − 8δ, auf die Regel 1 − 8(v/c)²·e^{−s}(1 − e^{−s}),
    höchstens 2(v/c)². Für kleines q ist er 1 − 24q, nach der halben Drehung 1 + 24q."""
    return math.exp(24 * sum(math.log1p(-q**k) for k in range(1, n)))


def kugelmittel_quer():
    """Anteil der zwei Richtungen quer zu einer Achse im Mittel über die Kugel, ∫sin²θ dΩ/(4π) = 2/3."""
    from scipy.integrate import quad
    return quad(lambda th: math.sin(th)**2 * math.sin(th) * 2 * math.pi, 0, math.pi)[0] / (4 * math.pi)


def licht_gegen_neutrinos():
    """Licht gegen Neutrinos nach dem Zerstrahlen von Elektron und Positron: Diese zählen mit zwei Sinnen und zwei Ladungen
    4·(7/8) gegen die 2 des Lichts, T_γ/T_ν = (1 + 2·(7/8))^{1/3} = (11/4)^{1/3}."""
    return (1 + 2 * strahlung_je_sinn(2) / strahlung_je_sinn(1))**(1 / 3)


def zustand_gas(a, Fs, k, a_s):
    """Anteil der Baryonen im kalten Gas mit übergebener Mode über die Zeit, die Form des Kreises: früh null, solange nichts
    gebunden ist, wachsend mit der Bindung bis zur Sättigung a_s (Maximum F_s), danach fallend, wenn Sternbildung und Wärme
    das Gas aus den Bereichen nehmen: F(a) = F_s·u^k·e^{k(1 − u)} mit u = a/a_s, nie negativ."""
    import numpy as np
    u = np.asarray(a, dtype=float) / a_s
    return Fs * u**k * np.exp(k * (1 - u))


def gas_kreis(a, F0, k, a_s):
    """Das gebundene Gas im Kreis, durch seinen heutigen Wert gelegt: F(a) = F₀·a^k·e^{k(1 − a)/a_s}, dieselbe Form wie
    zustand_gas (früh null, Maximum an der Sättigung a_s, danach Rückfluss), F(1) = F₀. Gibt F(a) und den Gang der Dunklen
    Energie ρ_DE(a)/ρ_DE(1) = (1 + j·F(a))/(1 + j·F₀), j = (ν_G(−ln q) − 1)/b₁ bei Q = 2π − 1 (einmal gerechnet): Formteil
    b₁·ρ_s mal dem Zustandsteil des Gases."""
    import numpy as np
    if "je_gas" not in _TAFEL:
        _TAFEL["je_gas"] = zustandsteil_je_gas(q_je_radiant()) / B1_ZELLE
    j = _TAFEL["je_gas"]
    a = np.asarray(a, dtype=float)
    F = F0 * a**k * np.exp(k * (1 - a) / a_s)
    return F, (1 + j * F) / (1 + j * F0)


def gas_gekoppelt(F0, k, a_s, numerisch=False):
    """Der Zustandsteil aus der Mode gespeist: Was er gewinnt, verliert die Mode, die wie Staub mit a⁻³ verdünnt. Gegen die
    erhaltene Mode gelesen wirkt er auf die Dehnung mit seiner über das Volumen gemittelten Geschichte,
    F_eff(a) = 3a⁻³·∫₀^a F(a₁)·a₁²·da₁ mit F aus gas_kreis, geschlossen 3F₀·e^{k/a_s}·γ(k + 3, ka/a_s)/(a³·(k/a_s)^{k+3}).
    Gibt (Gang der Dunklen Energie ρ_DE(z)/ρ_DE(0) = (1 + j·F_eff)/(1 + j·F_eff(1)) als Funktion von z, F_eff(1)).
    numerisch: zweiter Weg, das Integral über die Quadratur."""
    import numpy as np
    from scipy.special import gammainc, gammaln
    if "je_gas" not in _TAFEL:
        _TAFEL["je_gas"] = zustandsteil_je_gas(q_je_radiant()) / B1_ZELLE
    j = _TAFEL["je_gas"]
    if numerisch:
        from scipy.integrate import quad
        def F_eff(a):
            return 3 / a**3 * quad(lambda x: gas_kreis(x, F0, k, a_s)[0] * x**2, 0, a, limit=200)[0]
        F_eff = np.vectorize(F_eff)
    elif k == 0:
        F_eff = lambda a: F0 + 0 * np.asarray(a, dtype=float)
    else:
        b = k / a_s
        F_eff = lambda a: 3 * F0 * np.exp(b + gammaln(k + 3) - (k + 3) * math.log(b)) * gammainc(k + 3, b * np.asarray(a)) \
            / np.asarray(a, dtype=float)**3
    heute = float(F_eff(1.0))
    return (lambda z: (1 + j * F_eff(1 / (1 + np.asarray(z, dtype=float)))) / (1 + j * heute)), heute


def zustandsteil_je_gas(Q=None):
    """Fluchtmaß der Bindung je gebundenem Baryon: Gezählt werden die Tori innerhalb der Naht, R_Q = r*/(−ln q) mit
    q = Q/(1 + Q); für jeden Bereich ist also s = −ln q, unabhängig von seiner Masse, und das inaktive Maß ν_G(−ln q) − 1.
    Mit ρ = ν·ρ_s ist der Zustandsteil (ν_G(−ln q) − 1)·F·ρ_s, F der als Gas gebundene Anteil der Baryonen."""
    Q = Q_kosm if Q is None else Q
    return fluchtmass_ganzes(-math.log(Q / (1 + Q))) - 1


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


def toomre_schwelle(x, M=None, a0_=a0, numerisch=False, dichte=1.0):
    """Scheibe an der Schwelle der Flächendichte, Σ₀ = a₀/(πG), am Ort R = x·R_d: a_N/a₀ = x·[I₀K₀ − I₁K₁](x/2) (Freeman),
    längs antwortet die Regel mit ν(1 − χ/2), darum Q_T = κσ/(πG·ν(1 − χ/2)·Σ), κ² = R·dΩ²/dR + 4Ω². Gibt dict mit der
    Streuung der Gasscheibe an Q_T = 1 in Einheiten des Bahntempos v am Ort und von v_flach, dazu dieselbe für den
    Beobachter mit starrem Halo und mit M [kg] in m/s. Dazu die Arme: An der Spitze der Scherverstärkung stehen die Kämme
    radial, k quer zur Neigung, die Welle entsteht mit G·ν, m_T = 2πR/λ_T = κ²R/(2πG·ν·Σ) = (4 − d_eff)·g_N/(2πGΣ) mit
    d_eff = 4 − κ²/Ω², ν kürzt sich (arme). Newton mit den Baryonen allein (arme_newton), der Beobachter mit starrem
    Halo rechnet ν·m_T (arme_halo). Die Mitte einer Gasschicht nach Spitzer mit G·ν, ρ₀ = πG·ν·Σ²/(2σ²), gegen die
    Gezeitendichte ρ_J = 3(d_eff − 1)Ω²/(4πG): ρ₀/ρ_J = (2/3)·(4 − d_eff)/((d_eff − 1)·Q_T²·ν·(1 − χ/2)²), bei Q_T = 1
    (gezeit, unendlich für d_eff ≤ 1). numerisch: a_N als Hankel-Integral, κ aus Differenzen. dichte: Σ₀ in Einheiten der
    Schwelle a₀/(πG), bei gleicher Masse R_d = r*/√(2·dichte), v_flach bleibt."""
    from scipy import special as S
    if numerisch:
        from scipy.integrate import quad

        def u_(z):
            f = lambda k: k * S.j1(k * z) / (1 + k * k)**1.5
            return 2 * sum(quad(f, n * math.pi / z, (n + 1) * math.pi / z)[0] for n in range(4000))
        a_ = lambda z: verstaerkung(dichte * u_(z) * a0_, a0_) * dichte * u_(z)
        h = 1e-3 * x
        om2 = lambda z: a_(z) / z
        ko2 = (x * (om2(x + h) - om2(x - h)) / (2 * h) + 4 * om2(x)) / om2(x)
        u = dichte * u_(x)
    else:
        y = x / 2
        B = S.i0(y) * S.k0(y) - S.i1(y) * S.k1(y)
        dB = 2 * (S.i1(y) * S.k0(y) - S.i0(y) * S.k1(y)) + 2 * S.i1(y) * S.k1(y) / y
        u = dichte * x * B
    nu, L = verstaerkung(u * a0_, a0_), 1 - schalter(u * a0_, a0_) / 2
    if not numerisch:
        ko2 = 3 + L * (1 + y * dB / B)
    ko = math.sqrt(ko2)
    sv, sh = float(L * dichte * math.exp(-x) / (u * ko)), float(dichte * math.exp(-x) / (u * ko * nu))
    v_vf = math.sqrt(nu * u * x / math.sqrt(2 * dichte))
    m_T = float(ko2 * u / (2 * dichte * math.exp(-x)))
    d_ = float(4 - ko2)
    aus = dict(nu=nu, laengs=nu * L, sigma_v=sv, sigma_vf=sv * v_vf, halo_v=sh, d_eff=d_, arme=m_T,
               arme_newton=float((3 + (ko2 - 3) / L) * u / (2 * dichte * math.exp(-x))), arme_halo=float(nu * m_T),
               gezeit=(2 / 3) * (4 - d_) / ((d_ - 1) * nu * L**2) if d_ > 1 else math.inf)
    if M is not None:
        v = math.sqrt(nu * u * a0_ * x * schwellenradius(M, a0_) / math.sqrt(2 * dichte))
        aus.update(sigma=sv * v, halo=sh * v)
    return aus


def hi_linie(s, h_gas=2.0, x_max=None, a0_=a0, punkte=4000):
    """Linie des neutralen Gases einer Scheibe an der Schwelle der Flächendichte (toomre_schwelle), von der Seite unter der
    Neigung i gesehen, in Einheiten von v_flach·sin i: Das Gas liegt exponentiell mit der Länge h_gas·R_d (bis x_max·R_d,
    Vorgabe 8·h_gas), läuft auf der Kurve der Regel, v/v_flach = √(ν(u)·u·x/√2) mit u = a_N/a₀ = x·[I₀K₀ − I₁K₁](x/2), und
    streut in jeder Richtung gleich mit s = σ/(v_flach·sin i). Breiten wie bei ALFALFA: am Anteil 0,5 und 0,2 der Spitze
    jeder Seite, außen gelesen. Gibt (W50, W20) in Einheiten von v_flach·sin i."""
    import numpy as np
    from scipy import special as S
    x_max = 8 * h_gas if x_max is None else x_max
    x = np.linspace(x_max / punkte, x_max, punkte)
    u = x * (S.i0(x / 2) * S.k0(x / 2) - S.i1(x / 2) * S.k1(x / 2))
    v = np.sqrt(np.array([verstaerkung(w * a0_, a0_) for w in u]) * u * x / math.sqrt(2))
    phi = (np.arange(2000) + 0.5) * math.pi / 2000
    d = 1e-3
    kanten = np.arange(-1.3 - 6 * s, 1.3 + 6 * s + d, d)
    werte = (v[:, None] * np.cos(phi)[None, :]).ravel()
    gew = np.repeat(x * np.exp(-x / h_gas), phi.size)
    P = np.histogram(werte, kanten, weights=gew)[0] + np.histogram(-werte, kanten, weights=gew)[0]
    if s > d:
        k = np.arange(-int(6 * s / d), int(6 * s / d) + 1) * d
        from scipy.signal import fftconvolve
        P = fftconvolve(P, np.exp(-k**2 / (2 * s * s)), mode="same")
    m = 0.5 * (kanten[1:] + kanten[:-1])
    re, li = m >= 0, m < 0

    def breite(f):
        r = m[re][np.nonzero(P[re] >= f * P[re].max())[0][-1]]
        l = m[li][np.nonzero(P[li] >= f * P[li].max())[0][0]]
        return r - l
    return float(breite(0.5)), float(breite(0.2))


_LINIE = {}


def hi_linie_tafel(h_gas=2.0):
    """Tafel der Linie einer Scheibe an der Schwelle (hi_linie): s = σ/(v_flach·sin i) von 0 bis 10 → (s, W50, W20) in
    Einheiten von v_flach·sin i."""
    import numpy as np
    if h_gas not in _LINIE:
        s = np.concatenate([[0.0], np.geomspace(0.005, 10.0, 81)])
        w = np.array([hi_linie(x, h_gas=h_gas) for x in s])
        _LINIE[h_gas] = (s, w[:, 0], w[:, 1])
    return _LINIE[h_gas]


def gewicht_schicht(Sigma, g_R=0.0, sigma=None, a0_=a0, numerisch=False):
    """Gewicht einer ebenen Schicht der Flächendichte Σ [kg/m²] in einem Bereich, dessen Mode übergeben ist, im Feld g_R
    [m/s²] längs der Schicht. Senkrecht wächst das Feld nach Newton mit der Masse zwischen −z und z, u = 4πG∫₀^z ρ dz′,
    bis 2πGΣ, die Regel liest den Betrag, a_z = ν(√(g_R² + u²))·u. Der Druck in der Mitte ist für jedes Profil
    P₀ = ∫a_z du/(4πG) = (a₀²/2πG)·∫_{s_R}^{s_o} x³·ν(x) dx mit s_R² = g_R/a₀, s_o⁴ = s_R⁴ + (2πGΣ/a₀)²: Die Säule
    durchläuft die Leiter von s_R bis s_o, der spontane Teil x³ gibt Spitzer, πGΣ²/2, die Rückkehr x³·n_B die Stufe d = 3,
    höchstens π³a₀²/(30G) = ρ_s·c². Mit sigma [m/s] (isotherm) dazu die Höhe z_½, innerhalb der die Hälfte der Masse liegt,
    aus du/dz = 4πG·P(u)/σ². numerisch: die Schicht als Anfangswertaufgabe dρ/dz = −ρ·a_z/σ², du/dz = 4πGρ, ρ₀ so gesucht,
    dass u außen 2πGΣ erreicht. Gibt dict(P0 [Pa], spitzer = P₀/(πGΣ²/2), z_halb [m])."""
    from scipy.integrate import quad, solve_ivp
    from scipy.optimize import brentq
    u_o = 2 * math.pi * G * Sigma
    s_ = lambda u: (g_R**2 + u**2)**0.25 / math.sqrt(a0_)
    stufe = lambda s: s**4 / 4 + quad(lambda x: x**3 / math.expm1(x), 0, min(s, 60.0))[0] if s > 0 else 0.0
    oben = stufe(s_(u_o))
    P = lambda u: a0_**2 / (2 * math.pi * G) * (oben - stufe(s_(u)))
    zh = None
    if numerisch:
        sig = 1e4 if sigma is None else sigma

        def lauf(r0, ende=False):
            def f(z, y):
                rho, u = max(y[0], 0.0), y[1]
                g = math.hypot(g_R, u)
                return [-rho * (u / -math.expm1(-math.sqrt(g / a0_)) if g > 0 else 0.0) / sig**2, 4 * math.pi * G * rho]
            s = solve_ivp(f, (0, 200 * sig**2 / math.sqrt(a0_ * u_o)), [r0, 0.0], rtol=1e-12, atol=1e-40,
                          method="DOP853", dense_output=True)
            return (s.y[1, -1], s) if ende else s.y[1, -1] - u_o
        r_N = u_o**2 / (8 * math.pi * G * sig**2)
        r0 = brentq(lauf, r_N, r_N * 1e6, xtol=1e-14 * r_N)
        P0 = r0 * sig**2
        if sigma is not None:
            s = lauf(r0, True)[1]
            zh = brentq(lambda z: s.sol(z)[1] - u_o / 2, 0, s.t[-1])
    else:
        P0 = P(0.0)
        if sigma is not None:
            zh = quad(lambda u: sigma**2 / (4 * math.pi * G * P(u)), 0, u_o / 2)[0]
    return dict(P0=P0, spitzer=P0 / (math.pi * G * Sigma**2 / 2), z_halb=zh)

def letzte_stabile_bahn(M, regel=True, a0_=a0):
    """Letzte stabile Kreisbahn um einen gebundenen Bereich mit der Rückstellung des Endzustands: Ω² = g/R − H_∞²,
    κ² = (4 − d_eff)·g/R − 4H_∞², stabil, solange R²Ω nach außen wächst, bis (1 − d_eff/4)^{1/d_eff}·R_max mit
    Ω² = d_eff/(4 − d_eff)·H_∞²: um eine Galaxie mit übergebener Mode (d_eff = 2) R_max/√2 mit Ω = H_∞, um einen Haufen
    ((1 + Q)·Newton, d_eff = 3, M = (1 + Q)·M_b) R_max/∛4 mit Ω = √3·H_∞. Gibt (R, Ω)."""
    R, H = groesste_umkehr(M, a0_)
    d = 2 if regel else 3
    if not regel:
        R = (G * M / H**2)**(1 / 3)
    return ((4 - d) / 4)**(1 / d) * R, math.sqrt(d / (4 - d)) * H

def bindungsmasse(sigma, a0_=a0, traegt_selbst=True):
    """Kleinste Masse, die Gas der Streuung σ im Tiefbereich bindet. Trägt das Gas sein Feld selbst: M = (81/4)·σ⁴/(G·a₀)
    (Umkehr von σ⁴ = (4/81)GMa₀). Ist es Probekörper im Feld der Sterne: endliche Masse nur mit σ < v_flach/√3,
    also M = 9·σ⁴/(G·a₀)."""
    return (81 / 4 if traegt_selbst else 9.0) * sigma**4 / (G * a0_)


def z_gleichheit():
    """Gleichheit von Strahlung und Faltung (mit 1 + Q gelesen): 1 + z_eq = ω_m/ω_r."""
    rho_r, _, rho_m = dichte_heute()
    return rho_m / rho_r - 1


def schwelle_in_der_zeit(z, photon=False):
    """Die Schwelle in der Zeit. Λ = (1 + Q)·ΣΔs_b ist die ganze Faltung der Zelle mit ihren Rückläufen, und jede Energie
    faltet: Mit der Strahlung, die ohne Rückläufe zählt wie in der Faltung je Photon, faltet die Zelle je mitbewegtem
    Volumen ρ_m + ρ_r/a. Früh war Λ länger und a₀ = c²/(2Λ) kleiner, a₀(z)/a₀(heute) = (1 + a_eq)/(1 + (1 + z)·a_eq) mit
    a_eq = 1/(1 + z_eq). photon: z_eq über die Energie je Photon (faltung_je_photon, kT_eq/kT_CMB) statt ω_m/ω_r.
    Gibt a₀(z)/a₀(heute), z = −1 der Endzustand."""
    zeq = faltung_je_photon()["kT_eq"] * e / (k_B * T_CMB) - 1 if photon else z_gleichheit()
    aeq = 1 / (1 + zeq)
    return (1 + aeq) / (1 + (1 + z) * aeq)


def faltung_je_photon(Q=None, mu=None):
    """Das frühe Ganze an einer zeitfesten Zahl. Licht und Baryonen verdünnen sich mit der Zelle gleich, fest ist die Faltung
    der Baryonen je Photon ε_b = ρ_b·c²/n_γ, n_γ = (2ζ(3)/π²)·(kT/ħc)³ die Stufe d = 2 der Bose-Form, die Energie je Photon
    (π⁴/15)/(2ζ(3))·kT das Verhältnis der Stufen d = 3 und d = 2, mit den Neutrinos ρ_r/ρ_γ. Gleichheit von Strahlung und
    Faltung, wo diese Energie die Faltung je Photon mit ihren Rückläufen erreicht, (1 + Q)·ε_b, 1 + z_eq = kT_eq/kT_CMB.
    Nach der Rekombination, solange das Gas die Wärme des Lichts trägt, ist (ℓ/λ_H)³ = (π²/2ζ(3))·n_γ/n fest
    (ℓ = n^{−1/3}, n = ρ_b/(μ·m_H), λ_H = ħc/kT). Mitatmend liest die Neigung (1 + Q)·Newton, ℓ_J² = kT/(4πG·(1 + Q)·ρ_b·m),
    die Jeans-Masse M_J = (1 + Q)·ρ_b·ℓ_J³ = (ℓ/λ_H)^{3/2}·m_Pl³/((4π)^{3/2}·√(1 + Q)·m²) ist zeitfest (mit 1 + Q gelesen).
    Vorher hält das Licht, c_s² = c²/(3(1 + R)), R = 3ρ_b/(4ρ_γ), an der letzten Streufläche M_J = (1 + Q)·ρ_b·ℓ_J³ mit
    ℓ_J² = c_s²/(4πG·(1 + Q)·ρ_b), dort ℓ_J·H/c, H aus Faltung und Strahlung. Sprung (m·c²/(3(1 + R)·kT_*))^{3/2}.
    Q = None: gemessen (ω_c/ω_b). mu = None: neutrales Gas der Urmischung, μ = 1/(X + Y_p/4).
    Gibt dict (eps_b [eV], kT_eq [eV], mu, l_lam = ℓ/λ_H, M_nach [M☉], R, M_vor [M☉], lJ_H = ℓ_J·H/c, sprung)."""
    Q = Q_kosm if Q is None else Q
    Y = YP_GEMESSEN[0]
    mu = 1 / (1 - Y + Y / 4) if mu is None else mu
    m = mu * m_H
    kT = k_B * T_CMB
    n_g = strahlung_je_sinn(1, 2) / math.pi**2 * (kT / (hbar * c))**3
    e_je = strahlung_je_sinn(1, 3) / strahlung_je_sinn(1, 2)
    rho_r, rho_b, _ = dichte_heute()
    f_nu = rho_r / dichte_photonen()
    eps_b = rho_b * c**2 / n_g
    kT_eq = (1 + Q) * eps_b / (e_je * f_nu)
    l3 = math.pi**2 / strahlung_je_sinn(1, 2) * n_g / (rho_b / m)
    M_nach = math.sqrt(l3) * m_Pl**3 / ((4 * math.pi)**1.5 * math.sqrt(1 + Q) * m**2)
    a_s = 1 / (1 + z_stern)
    R = 3 * rho_b * a_s / (4 * dichte_photonen())
    rho_g = (1 + Q) * rho_b / a_s**3
    cs2 = c**2 / (3 * (1 + R))
    lJ = math.sqrt(cs2 / (4 * math.pi * G * rho_g))
    M_vor = rho_g * lJ**3
    H_s = math.sqrt(8 * math.pi * G / 3 * (rho_g + rho_r / a_s**4))
    return dict(eps_b=eps_b / e, kT_eq=kT_eq / e, mu=mu, l_lam=l3**(1 / 3), M_nach=M_nach / M_sun, R=R,
                M_vor=M_vor / M_sun, lJ_H=lJ * H_s / c, sprung=(m * c**2 / (3 * (1 + R) * kT / a_s))**1.5)


def sunjajew_seldowitsch(d=4):
    """Licht der Hintergrundstrahlung im heißen Plasma, erste Ordnung in δ_e = kT_e/(m_e c²) (Kompaneets, T_e ≫ T_CMB):
    Δn = −y·x⁻²·∂_x(x⁴·⟨Δn²⟩), x = hν/kT_CMB, mit der Schwankung ⟨Δn²⟩ = n_B·(1 + n_B) = −∂_x n_B = eˣ/(eˣ − 1)² und
    y = ∫δ_e·n_e·σ_T·dl. Die Stufe d > 2 der Schwankung, s^d·⟨Δn²⟩, hat ihren Gipfel bei s·coth(s/2) = s·(2ν − 1) = d
    (die Stufe d > 1 der Rückläufe s^d·n_B bei s·ν = d), d = 4 ist die Nullstelle. Integriert ∫s^d·⟨Δn²⟩ds = d·Γ(d)ζ(d):
    Die Energie wächst um 4y, die Zahl der Photonen (Stufe d = 2) bleibt. Rayleigh-Jeans ΔT/T = x·Δn → −2y.
    Gibt (Gipfel x₀ der Stufe d, ν₀ = x₀·kT_CMB/h [Hz], Energie je y, Zahl je y, ΔT/T je y bei x → 0)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    x0 = brentq(lambda s: s / math.tanh(s / 2) - d, 1e-6, 2.0 * d, xtol=1e-14)
    schw = lambda x: math.exp(-x) / math.expm1(-x)**2
    dn = lambda x: x * schw(x) * (x / math.tanh(x / 2) - 4)          # Δn/y
    energie = quad(lambda x: x**3 * dn(x), 1e-9, 80, limit=400)[0] / strahlung_je_sinn(1, 3)
    zahl = quad(lambda x: x**2 * dn(x), 1e-9, 80, limit=400)[0] / strahlung_je_sinn(1, 2)
    x = 1e-4
    return x0, x0 * k_B * T_CMB / h, energie, zahl, x * dn(x)


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


def cmb_spektren(ombh2, omch2, H0, ln10As, ns, tau, lmax=2508, alles=False, wachstum=(), w0=-1.0, wa=0.0, w_tafel=None):
    """Spektren der Hintergrundstrahlung D_ℓ = ℓ(ℓ + 1)C_ℓ/2π [µK²] von ℓ = 2 bis lmax (TT, TE, EE, gelinst). Die Mode ξ wirkt
    linear wie Staub mit ω_c, die Dunkle Energie als Formteil (w = −1); gerechnet mit CAMB, N_eff und eine Neutrinomasse
    von 0,06 eV wie bei Planck. Die Anfangsamplitude A_s und die Neigung n_s sind Inhalt, keine Form des Modells.
    alles: dazu das Linsenpotential [L(L + 1)]²C_L^φφ/2π von L = 2 bis 2500 und der Schallhorizont r_d [Mpc];
    wachstum: Rotverschiebungen (fallend), für die dazu f·σ₈ und σ₈ heute kommen (lineares Wachstum);
    w₀, w_a: Gang der Dunklen Energie, w(a) = w₀ + w_a(1 − a) (PPF), w = −1 der Formteil allein; w_tafel: (a, w) als
    Tafel statt w₀, w_a (PPF)."""
    import camb
    import numpy as np
    extra = dict(redshifts=list(wachstum) + [0.0], kmax=2.0) if wachstum else {}
    if (w0, wa) != (-1.0, 0.0) and w_tafel is None:
        extra.update(w=w0, wa=wa, dark_energy_model="ppf")
    par = camb.set_params(ombh2=ombh2, omch2=omch2, H0=H0, As=math.exp(ln10As) * 1e-10, ns=ns, tau=tau, mnu=0.06,
                          lmax=lmax + 300, lens_potential_accuracy=1, **extra)
    if w_tafel is not None:
        par.set_dark_energy_w_a(np.asarray(w_tafel[0]), np.asarray(w_tafel[1]), dark_energy_model="ppf")
    erg = camb.get_results(par)
    cl = erg.get_cmb_power_spectra(par, CMB_unit="muK", raw_cl=False)["total"]
    aus = (cl[2:lmax + 1, 0], cl[2:lmax + 1, 3], cl[2:lmax + 1, 1])
    if alles:
        aus += (erg.get_lens_potential_cls(2500)[2:, 0], erg.get_derived_params()["rdrag"])
    if wachstum:
        aus += (erg.get_fsigma8()[:-1], erg.get_sigma8_0())
    return aus


def profil_fourier(r, M, k):
    """Fourier-Bild einer runden Verteilung aus ihrer eingeschlossenen Masse M(r) bis r[-1]:
    u(k) = ∫ (dM/dr)·sin(kr)/(kr) dr / M(r[-1]), u → 1 für k → 0 (k in Kehrwerten der Länge von r)."""
    import numpy as np
    kr = np.outer(np.atleast_1d(k), r)
    return np.trapezoid(np.gradient(M, r) * np.sinc(kr / math.pi), r, axis=1) / M[-1]


def scherung_cl(ombh2, omch2, H0, ln10As, ns, nz, w0=-1.0, wa=0.0, ell=None, teile=False, daempfung=None):
    """Scherung der Bilder durch die Faltung längs des Lichtwegs, Limber: C_ℓ^{ij} für die Paare der Rotverschiebungsbins
    (nz = (z, [n_1(z), …]) als Histogramm). Die Mode wirkt wie druckfreier Staub (nichtlinear nach CAMB/HMcode 2020 ohne
    Rückwirkung der Baryonen), der Gang der Dunklen Energie wie bei cmb_spektren. Getrennt: GG (Faltung), GI und II (Ausrichtung
    der Bilder am Feld nach dem linearen Modell, je Einheit der Amplitude A, C1·ρ_krit = 0,0134). teile: je Paar getrennt der
    lineare Teil P_lin und der nichtlineare Rest P_nl − P_lin, den vor allem die gebundenen Bereiche tragen; daempfung: (k [1/Mpc], S),
    der Rest mal S(k) (in ln k gelesen, außerhalb fest). Gibt (ℓ, {(i, j): (GG, GI, II)}, S8), mit teile je Paar
    ((GG, GI, II) linear, (GG, GI, II) Rest)."""
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
    lies = lambda f: np.where(kk < 60.0, f(np.broadcast_to(zz, kk.shape).ravel(), np.minimum(kk, 60.0).ravel(), grid=False)
                              .reshape(kk.shape), 0.0)             # jenseits k = 60/Mpc trägt kein Bin mehr
    lin = lies(Pl)
    rest = lies(P) - lin
    if daempfung is not None:
        rest = rest * np.exp(np.interp(np.log(kk), np.log(daempfung[0]), np.log(daempfung[1])))
    stuecke = (lin, rest) if teile else (lin + rest,)
    aus = {}
    w = dchi_dz / chi**2
    for i in range(len(n)):
        for j in range(i, len(n)):
            je = tuple((np.trapezoid(Pk * q[i] * q[j] * w, zz, axis=1),
                        np.trapezoid(Pk * F * (q[i] * nchi[j] + nchi[i] * q[j]) * w, zz, axis=1),
                        np.trapezoid(Pk * F**2 * nchi[i] * nchi[j] * w, zz, axis=1)) for Pk in stuecke)
            aus[(i, j)] = je if teile else je[0]
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


def atmen_punkt(punkt):
    """Hintergrund im Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀, w_a): Faltung (mit 0,06 eV Neutrinomasse), Strahlung glatt,
    Dunkle Energie mit dem Gang w = w₀ + w_a(1 − a), flach. Gibt (ρ_m, ρ_r, ρ_DE(a), w(a), H(a) [1/s])."""
    wb, wc, H0, _, _, w0, wa = punkt
    rho_m = (wb + wc + 0.06 / 93.14) * rho_krit100
    rho_r = dichte_heute()[0]
    rho_de0 = 3 * (H0 * 1e3 / Mpc)**2 / (8 * math.pi * G) - rho_m - rho_r
    rho_de = lambda a: rho_de0 * gang_dunkle_energie(w0, wa, 1 / a - 1)
    w = lambda a: w0 + wa * (1 - a)
    H = lambda a: math.sqrt(8 * math.pi * G / 3 * (rho_m / a**3 + rho_r / a**4 + rho_de(a)))
    return rho_m, rho_r, rho_de, w, H


def spektrum_linear(punkt, kmax=60.0):
    """Lineares Spektrum der Faltung heute im Punkt, die Mode wirkt linear wie Staub (CAMB). Gibt (k [1/Mpc], P [Mpc³])."""
    import camb
    wb, wc, H0, lnAs, ns, w0, wa = punkt
    extra = dict(w=w0, wa=wa, dark_energy_model="ppf") if (w0, wa) != (-1.0, 0.0) else {}
    par = camb.set_params(ombh2=wb, omch2=wc, H0=H0, As=math.exp(lnAs) * 1e-10, ns=ns, mnu=0.06, **extra)
    par.set_matter_power(redshifts=[0.0], kmax=kmax)
    k, _, P = camb.get_results(par).get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    return k, P[0]


def wachstum_ende(th):
    """Ende des linearen Wachstums im Atmen des Formteils ohne Strahlung (die Mode atmet mit und wirkt wie Staub):
    D(a) = (5/2)·a_×·√(1 + (a_×/a)³)·∫₀^{a/a_×} u^{3/2}(1 + u³)^{−3/2} du, früh D = a, a_× die Dehnung bei Formteil gleich
    Faltung, D_∞ = (5/6)·B(5/6, 2/3)·a_× ohne a₀. Der Stand hängt allein an Θ: (a_×/a)³ = (Θ_∞/Θ)² − 1, Θ_∞ = √15/π.
    Gibt (D_∞/a_×, D(a_×)/D_∞, D(Θ)/D_∞)."""
    from scipy.integrate import quad
    I = lambda U: quad(lambda u: u**1.5 / (1 + u**3)**1.5, 0, U)[0]
    I_inf = quad(lambda u: u**1.5 / (1 + u**3)**1.5, 0, math.inf)[0]
    x = (math.sqrt(15) / math.pi / th)**2 - 1
    return 2.5 * I_inf, math.sqrt(2) * I(1.0) / I_inf, math.sqrt(1 + x) * I(x**(-1 / 3)) / I_inf


def wachstum_linear(punkt):
    """Lineares Wachstum D(a) (D(1) = 1) und f = d ln D/d ln a im Atmen des Punkts, Beginn im Strahlungszeitalter mit
    der wachsenden Lösung D ∝ a + ⅔·a_eq. Gibt Funktionen D(a), f(a)."""
    import numpy as np
    from scipy.integrate import solve_ivp
    rho_m, rho_r, rho_de, w, H = atmen_punkt(punkt)
    a_eq = rho_r / rho_m

    def rhs(x, y):
        a = math.exp(x)
        tot = rho_m / a**3 + rho_r / a**4 + rho_de(a)
        dlnH = -0.5 * (3 * rho_m / a**3 + 4 * rho_r / a**4 + 3 * (1 + w(a)) * rho_de(a)) / tot
        return [y[1], -(2 + dlnH) * y[1] + 1.5 * rho_m / a**3 / tot * y[0]]
    x0 = math.log(1e-6)
    s = solve_ivp(rhs, (x0, 0.0), [1e-6 + 2 / 3 * a_eq, 1e-6], rtol=1e-10, atol=1e-14, dense_output=True)
    D1 = s.sol(0.0)[0]
    return (lambda a: s.sol(np.log(a))[0] / D1), (lambda a: s.sol(np.log(a))[1] / s.sol(np.log(a))[0])


def falte_phase():
    """Phase der Bahn von der Umkehr bis zur ersten Kreuzung: Eine kalte Schale fällt frei bis zur Mitte, dort trifft sie die
    inneren, ω·t = π/(2√2) mit ω² = GM/r_U³ (Zykloide r = r_U(1 + cos η)/2, ω·t = (η + sin η)/(2√2) bis η = π)."""
    return math.pi / (2 * math.sqrt(2))


def bereiche_massenfunktion(punkt, lgM):
    """Zahl der Bereiche je Dekade der Masse heute [1/Mpc³] im Punkt, M [M☉] mit der Mode gelesen (ω_b + ω_c): erster
    Durchgang der linearen Faltung durch die Schwelle des rollenden Kreises δ_c = (3/20)·(12π)^{2/3} ohne Rückstellung,
    dn/dlg M = (ρ_m/M)·|d ln σ/d lg M|·√(2/π)·ν·e^{−ν²/2}, ν = δ_c/σ, σ über die Kugel um M (spektrum_linear)."""
    import numpy as np
    k, P = spektrum_linear(punkt)
    rho_m = (punkt[0] + punkt[1]) * rho_krit100 / M_sun * Mpc**3
    lgM = np.asarray(lgM, dtype=float)

    def sig(m):
        x = k * (3 * 10**m / (4 * math.pi * rho_m))**(1 / 3)
        return math.sqrt(np.trapezoid(P * (3 * (np.sin(x) - x * np.cos(x)) / x**3)**2 * k**3, np.log(k)) / (2 * math.pi**2))
    s = np.array([sig(m) for m in lgM])
    nu = 0.15 * (12 * math.pi)**(2 / 3) / s
    return rho_m / 10**lgM * np.abs(np.gradient(np.log(s), lgM)) * math.sqrt(2 / math.pi) * nu * np.exp(-nu**2 / 2)


def linse_umgebung(punkt, M, R, z=0.2, weg="hankel"):
    """Linse der linearen Faltung um einen Bereich der Masse M [M☉] (mit der Mode gelesen) im Abstand R [Mpc, physikalisch]
    bei z, die Mode wie Staub: g = 4G·b·ΔΣ mit der Verstärkung b = 1 + (ν² − 1)/δ_c aus der Massenfunktion (ν = δ_c/σ(M)),
    ΔΣ = ρ̄·∫k·dk/2π·P(k)·J₂(kR) (weg "hankel") oder aus der Projektion von ξ (weg "xi"), mitbewegt gerechnet. Gibt g
    [m/s²]."""
    import numpy as np
    from scipy.special import jv
    k, P = spektrum_linear(punkt)
    P = P * wachstum_linear(punkt)[0](1 / (1 + z))**2
    rho = (punkt[0] + punkt[1]) * rho_krit100 / M_sun * Mpc**3
    x = k * (3 * M / (4 * math.pi * rho))**(1 / 3)
    sig = math.sqrt(np.trapezoid(P * (3 * (np.sin(x) - x * np.cos(x)) / x**3)**2 * k**3, np.log(k)) / (2 * math.pi**2))
    dc = 0.15 * (12 * math.pi)**(2 / 3)
    b = 1 + ((dc / sig)**2 - 1) / dc
    Rc = np.atleast_1d(np.asarray(R, dtype=float)) * (1 + z)
    if weg == "hankel":
        dS = np.array([rho * np.trapezoid(k**2 * P * jv(2, k * x_) / (2 * math.pi), np.log(k)) for x_ in Rc])
    else:
        rr = np.geomspace(1e-2, 200, 600)
        xi = np.array([np.trapezoid(P * k**3 * np.sinc(k * r / math.pi), np.log(k)) / (2 * math.pi**2) for r in rr])
        zz = np.geomspace(1e-4, 150, 800)
        S = lambda X: 2 * rho * np.trapezoid(np.interp(np.sqrt(X**2 + zz**2), rr, xi), zz)
        from scipy.integrate import cumulative_trapezoid
        Rg = np.geomspace(1e-4, 1.01 * Rc.max(), 600)
        Sg = np.array([S(X) for X in Rg])
        innen = np.concatenate([[0.0], cumulative_trapezoid(Sg * Rg, Rg)]) + Sg[0] * Rg[0]**2 / 2
        dS = np.array([2 * np.interp(X, Rg, innen) / X**2 - S(X) for X in Rc])
    return 4 * G * b * dS * M_sun / Mpc**2 * (1 + z)**2


def breitenfunktion(punkt, kanten, anteil=1.0, sigma=None, v_min=None):
    """Galaxien je Dekade der Breite W50 [1/Mpc³] in den Stücken zwischen den Kanten (lg W50 [km/s]), wenn jeder Bereich
    (bereiche_massenfunktion) den Anteil anteil seiner Baryonen M_b = M·ω_b/(ω_b + ω_c) als Galaxie hält: v_flach⁴ =
    G·anteil·M_b·a₀, W50 = 2·v_flach·sin i, cos i gleichverteilt (P(sin i < x) = 1 − √(1 − x²)). sigma [km/s]: statt
    dessen die Linie einer Scheibe an der Schwelle (hi_linie_tafel) mit dieser Streuung, W50 = V·w50(σ/V) mit
    V = v_flach·sin i, über σ/V = 10 hinaus die Glocke 2√(2 ln 2)·σ. v_min [km/s]: nur Bereiche, deren Tiefe vor der
    Übergabe, V₂₀₀ = (10·G·H₀·M)^{1/3} mit der Mode, darüber liegt."""
    import numpy as np
    lgM = np.linspace(6.0, 15.5, 381)
    dn = bereiche_massenfunktion(punkt, lgM) * (lgM[1] - lgM[0])
    if v_min is not None:
        dn = dn * ((10 * G * punkt[2] * 1e3 / Mpc * 10**lgM * M_sun)**(1 / 3) / 1e3 > v_min)
    v = (G * anteil * 10**lgM * punkt[0] / (punkt[0] + punkt[1]) * M_sun * a0)**0.25 / 1e3
    if sigma is not None:
        s_, w5, _ = hi_linie_tafel()
        n = 2000
        V = v[:, None] * np.sqrt(1 - ((np.arange(n) + 0.5) / n)**2)[None, :]
        s = sigma / V
        W50 = np.where(s <= s_[-1], V * np.interp(np.minimum(s, s_[-1]), s_, w5), 2 * math.sqrt(2 * math.log(2)) * sigma)
        lw = np.log10(W50)
        return np.array([(dn[:, None] / n * ((lw >= lo) & (lw < hi))).sum() for lo, hi in zip(kanten[:-1], kanten[1:])]) \
            / np.diff(kanten)
    F = lambda x: 1 - np.sqrt(1 - np.minimum(x, 1.0)**2)
    W = 10**np.asarray(kanten, dtype=float)
    return (dn[:, None] * (F(W[None, 1:] / (2 * v[:, None])) - F(W[None, :-1] / (2 * v[:, None])))).sum(0) / np.diff(kanten)


def faltung_spitze(punkt, R, R_L, spektrum=None):
    """Mittlere lineare Faltung innerhalb des mitbewegten Radius R [Mpc, auch Feld] um eine Spitze des Gaußfelds im Punkt,
    als Teil der Faltung innerhalb R_L [Mpc]: Δ(<R)/Δ(<R_L) = σ²(R_L, R)/σ²(R_L, R_L), σ²(R₁, R₂) = ∫P·W(kR₁)W(kR₂)k²dk/2π²
    (das Mittel des Felds unter der Bedingung). spektrum: (k, P) aus spektrum_linear, sonst gerechnet."""
    import numpy as np
    k, P = spektrum_linear(punkt) if spektrum is None else spektrum
    W = lambda x: 3 * (np.sin(x) - x * np.cos(x)) / x**3
    sig2 = lambda R1, R2: np.trapezoid(P * W(k * R1) * W(k * R2) * k**3, np.log(k)) / (2 * math.pi**2)
    aus = np.array([sig2(R_L, x) for x in np.atleast_1d(R)]) / sig2(R_L, R_L)
    return aus if np.ndim(R) else float(aus[0])


def einfall_haufen(punkt, R_L, n=1000, innen=0.02, aussen=3.0, a_start=0.02, weich=10.0, schritt=0.03, mittel=1.0,
                   mischen=0, quer=None):
    """Einfall der Mode als druckfreier Staub, gerechnet als Entstehung: Um einen Ort, dessen lineare Faltung innerhalb des
    mitbewegten Radius R_L [Mpc] die Schwelle hat, liegt im Mittel die Faltung Δ(<R) = Δ_L·σ²(R_L, R)/σ²(R_L, R_L),
    σ²(R₁, R₂) = ∫P·W(kR₁)W(kR₂)k²dk/2π² (eine Bedingung an das Gaußfeld des Spektrums). Δ_L ist so hoch, dass die Schale
    bei R_L ohne Kreuzung heute die Mitte erreicht. n Schalen von innen·R_L bis aussen·R_L (Abstand gleich in ln R) fallen
    radial im Atmen des Punkts, jede spürt die Masse innerhalb ihres Orts (Gauß), durch die Mitte gehen sie hindurch
    (innerhalb weich [kpc], mitbewegt, wächst die Kraft linear). mischen = k: Einen Radiant der Bahn nach der Umkehr (Zeit 1/ω, ω² = GM/r_U³), wo das Blatt
    am Ort der Schale faltet, mischt die Schale ihre Bahnen, sie teilt sich in k Teile gleicher Masse, deren Richtung gleich
    verteilt ist (μ = cos des Winkels zur Radialen gleichmäßig in k Stufen, Betrag und Energie der Geschwindigkeit bleiben,
    der quere Teil wird Drehimpuls). quer = τ: An derselben Stelle gibt jede Schale den Anteil τ ihrer Bewegungsenergie in
    die Querrichtung, als Drehimpuls, der Rest bleibt radial. quer = "scherung": τ ist der Anteil der Scherung an der
    Verformung der Spitze, 2/(ν² + 2) mit ν = Δ_L/σ(R_L) ihrer Höhe über der Schwankung des Felds (Spur und Scherung sind
    unabhängig, die Spur trägt ν²σ²/3, die Scherung 2σ²/3).
    Gibt dict: r [m] und M(<r) [kg] heute, über die letzten mittel [Gyr] gemittelt, dazu je Schale die Umkehr (t, r) und
    die erste Kreuzung (t, r)."""
    import numpy as np
    from scipy.integrate import quad, solve_ivp
    rho_m, rho_r, rho_de, w, H = atmen_punkt(punkt)
    R = R_L * np.geomspace(innen, aussen, n)
    D, f = wachstum_linear(punkt)
    t_von_a = lambda a: quad(lambda x: 1 / H(math.exp(x)), math.log(1e-8), math.log(a), epsrel=1e-10, limit=200)[0]
    t0, ti = t_von_a(1.0), t_von_a(a_start)
    ah = solve_ivp(lambda t, y: [y[0] * H(y[0])], (ti, 1.05 * t0), [a_start], rtol=1e-11, dense_output=True)
    a_t = lambda t: float(ah.sol(t)[0])
    M = 4 * math.pi / 3 * rho_m * (R * Mpc)**3                          # Masse innerhalb jeder Schale (Lagrange)
    m = np.diff(np.concatenate([[0.0], M]))
    eps = weich * kpc

    Mm = M - 0.5 * m                                                    # Masse innerhalb der Mitte jeder Schale
    Rm = R * (Mm / M)**(1 / 3)
    prof_m = faltung_spitze(punkt, Rm, R_L)

    def zu_beginn(DL):
        d = DL * prof_m * D(a_start)
        r = a_start * Rm * Mpc * (1 + d)**(-1 / 3)
        return r, H(a_start) * r * (1 - f(a_start) * d / (3 * (1 + d)))

    def ankunft(DL):                                                    # Schale R_L allein: Zeit bis zur Mitte
        i = np.argmin(abs(R - R_L))
        r0, v0 = zu_beginn(DL)
        ext = lambda t, r: -4 * math.pi * G / 3 * (2 * rho_r / a_t(t)**4 + rho_de(a_t(t)) * (1 + 3 * w(a_t(t)))) * r
        ev = lambda t, y: y[0] - 1e-4 * r0[i] / a_start
        ev.terminal = True
        s = solve_ivp(lambda t, y: [y[1], -G * Mm[i] / y[0]**2 + ext(t, y[0])], (ti, 1.04 * t0), [r0[i], v0[i]],
                      events=ev, rtol=1e-10)
        return s.t_events[0][0] if len(s.t_events[0]) else 2 * t0
    from scipy.optimize import brentq
    DL = brentq(lambda x: ankunft(x) - t0, 1.0, 3.0, xtol=1e-6)
    if isinstance(quer, str):                                           # Anteil der Scherung an der Spitze
        k_s, P_s = spektrum_linear(punkt)
        Ws = lambda y: 3 * (np.sin(y) - y * np.cos(y)) / y**3
        hoehe = DL / math.sqrt(np.trapezoid(P_s * Ws(k_s * R_L)**2 * k_s**3, np.log(k_s)) / (2 * math.pi**2))
        quer = 2 / (hoehe**2 + 2)
    k_ = max(mischen, 1)
    r, v = (np.repeat(x, k_) for x in zu_beginn(DL))
    mk = np.repeat(m / k_, k_)
    jj = np.zeros(n * k_)                                               # Drehimpuls je Masse, nach dem Mischen
    mu = -1 + (2 * np.arange(k_) + 1) / k_
    t = ti
    t_um, r_um = np.full(n, np.nan), np.full(n, np.nan)
    t_kr, r_kr, t_mi = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    erste = np.arange(n) * k_
    gitter = np.geomspace(1 * kpc, 20 * Mpc, 400)
    summe, zahl = np.zeros_like(gitter), 0

    def beschl(r, t):
        o = np.argsort(r, kind="stable")
        Me = np.empty(n * k_)
        Me[o] = np.cumsum(mk[o]) - 0.5 * mk[o]
        a = a_t(t)
        e_ = eps * a
        return (-G * Me * r / np.maximum(r, e_)**3 + jj**2 / np.maximum(r, e_)**3
                - 4 * math.pi * G / 3 * (2 * rho_r / a**4 + rho_de(a) * (1 + 3 * w(a))) * r), o
    g, o = beschl(r, t)
    while t < t0:
        s_ = np.maximum(r, eps * a_t(t))
        dt = min(1 / H(a_t(t)), float(np.min(np.sqrt(s_ / np.abs(g)))), float(np.min(s_ / np.hypot(v, jj / s_))))
        dt = min(schritt * dt, t0 - t)
        v_alt = v[erste].copy()
        v += 0.5 * dt * g
        r += dt * v
        innen_ = r < 0
        r[innen_], v[innen_] = -r[innen_], -v[innen_]
        t += dt
        g, o = beschl(r, t)
        v += 0.5 * dt * g
        um = np.isnan(t_um) & (v_alt > 0) & (v[erste] <= 0)
        t_um[um], r_um[um] = t, r[erste][um]
        rang = np.empty(n * k_, int)
        rang[o] = np.arange(n * k_)
        kr = np.isnan(t_kr) & (rang[erste] != erste)
        t_kr[kr], r_kr[kr] = t, r[erste][kr]
        falt = np.isnan(t_mi) & (t >= t_um + np.sqrt(r_um**3 / (G * Mm)))
        t_mi[falt] = t
        if (mischen or quer is not None) and falt.any():
            for i in np.flatnonzero(falt):
                s_ = slice(i * k_, (i + 1) * k_)
                b = math.hypot(v[i * k_], jj[i * k_] / r[i * k_])
                if quer is None:
                    v[s_], jj[s_] = b * mu, r[s_] * b * np.sqrt(1 - mu**2)
                else:
                    tq = quer[i] if np.ndim(quer) else quer
                    v[s_], jj[s_] = np.sign(v[s_]) * b * math.sqrt(1 - tq), r[s_] * b * math.sqrt(tq)
            g, o = beschl(r, t)
        if t > t0 - mittel * Gyr:
            i = np.searchsorted(r[o], gitter)
            summe += dt * np.where(i > 0, np.cumsum(mk[o])[np.maximum(i - 1, 0)], 0.0)
            zahl += dt
    return dict(r=gitter, M=summe / zahl, t0=t0, t_um=t_um, r_um=r_um, t_kr=t_kr, r_kr=r_kr, M_schale=Mm, DL=DL,
                quer=quer)


def konzentration_profil(r, M, H0, von=0.05):
    """Konzentration c = R₂₀₀/r_s, wie die Haufen gelesen werden: NFW mit M(R₂₀₀) fest, in lg M von von·R₅₀₀ bis R₂₀₀
    angepasst (mittlere Dichte 500- und 200-mal die kritische heute). Gibt (c, R₂₀₀ [m], M₂₀₀ [kg], mittlerer Betrag des
    Rests in dex)."""
    import numpy as np
    from scipy.optimize import minimize_scalar
    rc = 3 * (H0 * 1e3 / Mpc)**2 / (8 * math.pi * G)
    ok = M > 0
    r, M = r[ok], M[ok]
    rho = 3 * M / (4 * math.pi * r**3)
    lr = np.log(r)

    def rand(D):                                                        # äußerster Schnitt mit D·ρ_krit
        f = np.log(rho / (D * rc))
        j = np.flatnonzero(f >= 0)[-1]
        return math.exp(lr[j] + (lr[j + 1] - lr[j]) * f[j] / (f[j] - f[j + 1]))
    R = {D: rand(D) for D in (500, 200)}
    M2 = math.exp(np.interp(math.log(R[200]), lr, np.log(M)))
    nfw = lambda x: np.log(1 + x) - x / (1 + x)
    w = (r >= von * R[500]) & (r <= R[200])
    rest = lambda c: np.log10(M[w] / (M2 * nfw(c * r[w] / R[200]) / nfw(c)))
    c = minimize_scalar(lambda c: float(np.sum(rest(c)**2)), bounds=(0.5, 40), method="bounded").x
    return float(c), R[200], M2, float(np.mean(np.abs(rest(c))))


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


def w0_saettigung(wa, wm, a0_=a0):
    """w₀ des Gangs w = w₀ + w_a(1 − a), dessen Maximum dort liegt, wo der Kreis sättigt, wo drei Teile die Faltung
    einholen: a_s³ = ω_m·ρ_krit100/ρ_DE(a₀), w₀ = −1 − w_a·(1 − a_s). wm: ω_m, a0_: a₀ [m/s²] der Dunklen Energie heute."""
    a_s = (wm * rho_krit100 / dichte_dunkle_energie(a0_))**(1 / 3)
    return -1 - wa * (1 - a_s)


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
