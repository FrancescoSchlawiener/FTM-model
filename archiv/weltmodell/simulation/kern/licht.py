# -*- coding: utf-8 -*-
"""
Licht, Rotation, Wellen (08; Ebene L, Bahnen von Licht und Körpern in der Faltung): Ablenkung (einzeln, kollektiv,
am kollektiven Trichter), Zusatzpotential, Photonensphäre, Schatten, Perihel, Verwindung (Kerr: Frequenzen,
Knotendrehung, Horizonte, Nullfläche des Takts, Mitführung, geodätische Drehung), Wellen (Leistung, Faltungsmesser,
Ringdown, Gedächtnis, Schwellenfrequenz, Hulse-Taylor). Gilt im Starkfeld wie die ART, in Δs geschrieben.
"""
from scipy.integrate import quad

from .konstanten import *
from .formeln import formel
from .faltung import ds, r_stern
from .galaxie import v_flach
from .rueckkehr import nu

O8 = "08-licht-und-trichter.md"


@formel('ablenkung', 'α_abl = 12Δs/b',
        ('α_abl',), ('M', 'b'), 'L', 'modell', 'Baustein 08', (O8, '= 12Δs/b'))
def ablenkung(M, b):
    return 12 * ds(M) / b


@formel('ablenkung_kollektiv', 'α_∞ = 6πΔs/r* = 2π(v_flach/c)²',
        ('α_∞',), ('M', 'a₀'), 'L', 'modell', 'Baustein 08', (O8, 'α_∞ = 6πΔs/r* = 2π·(v_flach/c)²'))
def ablenkung_kollektiv(M, a0=A0):
    """α_∞ = 6πΔs/r* = 2π(v_flach/c)²."""
    return 2 * math.pi * (v_flach(M, a0) / c)**2


@formel('perihel_pro_umlauf', 'Δφ_P = 18πΔs/(a_b(1 − e_B²))',
        ('Δφ_P',), ('M', 'a_b', 'e_B'), 'L', 'modell', 'Baustein 08', (O8, 'Δφ_Perihel = 18π·Δs/(a(1 − e²))'))
def perihel_pro_umlauf(M, a, e):
    return 18 * math.pi * ds(M) / (a * (1 - e * e))


@formel('schatten_winkel', 'α_S = 2b/D_A, b = 3√27·Δs = 3√3·GM/c²',
        ('α_S',), ('M', 'D_A'), 'L', 'modell', 'Baustein 08', (O8, 'b = 3√27·Δs ≈ 15,6·Δs'))
def schatten_winkel(M, D):
    return 2 * 3 * math.sqrt(3) * G * M / c**2 / D


@formel('isco_kerr', 'r̃_ISCO = 3 + z₂ ∓ √((3 − z₁)(3 + z₁ + 2z₂)) (Kerr, Bardeen), Einheit 3Δs',
        ('r̃_ISCO',), ('j',), 'L', 'modell', 'Baustein 08', (O8, 'Die innerste stabile Bahn liegt bei j = 0,9'))
def isco_kerr(j, mitlaeufig=True):
    """Innerste stabile Kreisbahn (Bardeen) in Einheiten GM/c² = 3Δs."""
    z1 = 1 + (1 - j * j) ** (1 / 3) * ((1 + j) ** (1 / 3) + (1 - j) ** (1 / 3))
    z2 = math.sqrt(3 * j * j + z1 * z1)
    sg = -1 if mitlaeufig else 1
    return 3 + z2 + sg * math.sqrt((3 - z1) * (3 + z1 + 2 * z2))


@formel('schwellenfrequenz_gw', 'f*_gw = a₀^{3/4}/(π(3c²Δs)^{1/4})',
        ('f*_gw',), ('M', 'a₀'), 'L', 'modell', 'Baustein 08', (O8, 'f*_gw = a₀^{3/4}/(π·(3c²Δs)^{1/4})'))
def schwellenfrequenz_gw(M, a0=A0):
    """f*_gw = a₀^{3/4}/(π(3c²Δs)^{1/4})."""
    return a0 ** 0.75 / (math.pi * (3 * c**2 * ds(M)) ** 0.25)


@formel('hulse_taylor', 'Ṗ_b = −(192π/5)·η_m·(2π·3Δs/(cP_b))^{5/3}·f(e_B)',
        ('Ṗ_b',), ('m₁', 'm₂', 'P_b', 'e_B'), 'L', 'modell', 'Baustein 08', (O8, 'dP_b/dt = −(192π/5)·η_m·(2π·3Δs/(c·P_b))^{5/3}·f(e)'))
def hulse_taylor(m1, m2, Pb, e):
    M = (m1 + m2) * M_sun
    eta = m1 * m2 / (m1 + m2)**2
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e * e)**3.5
    return -(192 * math.pi / 5) * eta * (2 * math.pi * 3 * ds(M) / (c * Pb))**(5 / 3) * fe


@formel('ablenkung_linse', 'α(b) = (12Δs/b)·I(b/r*), I(x) = ∫₀^{π/2} cos θ·ν(cos²θ/x²) dθ (ν über a_N/a₀)',
        ('α_abl',), ('M', 'b', 'a₀'), 'LI', 'modell', 'Baustein 08', (O8, 'α(b) = (12Δs/b)·I(b/r*)'), alternative='ablenkung')
def ablenkung_linse(M, b, a0=A0):
    """Licht am kollektiven Trichter: für b ≪ r* 12Δs/b, für b ≫ r* 6πΔs/r* + 6Δs/b."""
    x = b / r_stern(M, a0)
    I = quad(lambda t: math.cos(t) * float(nu(math.cos(t)**2 / x**2 * a0, a0)), 0, math.pi / 2, limit=400)[0]
    return 12 * ds(M) / b * I


@formel('zusatzpotential', 'Φ_X = v_flach²·ln ν(s) = −v_flach²·ln(1 − e^{−s}), s = r*/r (vom Zentrum gezählt)',
        ('Φ_X',), ('M', 'r', 'a₀'), 'LI', 'modell', 'Baustein 08', (O8, 'Φ_X = v_flach²·ln ν(s) = −v_flach²·ln(1 − e^{−s})'))
def zusatzpotential(M, r, a0=A0):
    """Logarithmus der mittleren Zahl der Durchläufe; jenseits von r* → v_flach²·ln(r/r*)."""
    return -v_flach(M, a0)**2 * math.log(-math.expm1(-r_stern(M, a0) / r))


@formel('photonensphaere', 'r_ph = 9Δs = 3GM/c²',
        ('r_ph',), ('M',), 'L', 'modell', 'Baustein 08', (O8, 'Photonensphäre r = 9Δs'))
def photonensphaere(M):
    return 9 * ds(M)


@formel('verwindung', 'a_J = J/(Mc)',
        ('a_J',), ('J', 'M'), 'L', 'modell', 'Baustein 08', (O8, 'a_J = J/(Mc)'))
def verwindung(J, M):
    """Zweite Faltungsgröße eines rotierenden Körpers; GJ/c³ = 3Δs·a_J ist eine Fläche."""
    return J / (M * c)


@formel('umlauf_differenz', 'T_prograd − T_retrograd = 4π·a_J/c (unabhängig von r)',
        ('ΔT_pr',), ('a_J',), 'L', 'modell', 'Baustein 08', (O8, 'T_prograd − T_retrograd = 4π·a_J/c'))
def umlauf_differenz(aJ):
    return 4 * math.pi * aJ / c


@formel('kerr_omega_phi', 'ω_φ = c³/(GM(r̃^{3/2} ± j)), r̃ = r/3Δs',
        ('ω_φ',), ('M', 'r̃', 'j'), 'L', 'modell', 'Baustein 08', (O8, 'ω_φ = c³/(GM(r̃^{3/2} ± j))'))
def kerr_omega_phi(M, rt, j, mitlaeufig=True):
    """Grundfrequenzen der Bahn im Abstand r̃ = r/3Δs (Einheit der Faltung); j < 0 oder mitlaeufig=False: gegenläufig."""
    return c**3 / (G * M * (rt**1.5 + (j if mitlaeufig else -j)))


@formel('kerr_omega_r', 'ω_r² = ω_φ²(1 − 6/r̃ + 8j·r̃^{−3/2} − 3j²/r̃²)',
        ('ω_r',), ('M', 'r̃', 'j'), 'L', 'modell', 'Baustein 08', (O8, 'ω_r² = ω_φ²(1 − 6/r̃ + 8j·r̃^{−3/2} − 3j²/r̃²)'))
def kerr_omega_r(M, rt, j, mitlaeufig=True):
    """Radiale Grundfrequenz; null an der innersten stabilen Bahn, innerhalb davon keine stabile Bahn (Fehler)."""
    jj = j if mitlaeufig else -j
    return kerr_omega_phi(M, rt, j, mitlaeufig) * math.sqrt(1 - 6 / rt + 8 * jj * rt**-1.5 - 3 * jj**2 / rt**2)


@formel('kerr_omega_theta', 'ω_θ² = ω_φ²(1 − 4j·r̃^{−3/2} + 3j²/r̃²)',
        ('ω_θ',), ('M', 'r̃', 'j'), 'L', 'modell', 'Baustein 08', (O8, 'ω_θ² = ω_φ²(1 − 4j·r̃^{−3/2} + 3j²/r̃²)'))
def kerr_omega_theta(M, rt, j, mitlaeufig=True):
    """Vertikale Grundfrequenz; ohne Rotation exakt ω_θ = ω_φ, die Verwindung trennt den dritten Zyklus ab."""
    jj = j if mitlaeufig else -j
    return kerr_omega_phi(M, rt, j, mitlaeufig) * math.sqrt(1 - 4 * jj * rt**-1.5 + 3 * jj**2 / rt**2)


@formel('knotendrehung', 'ζ = (ω_φ − ω_θ)/ω_φ ≈ 2j·(3Δs/r)^{3/2} (führende Ordnung)',
        ('ζ',), ('M', 'r', 'j'), 'L', 'modell', 'Baustein 08', (O8, 'ζ = (ω_φ − ω_θ)/ω_φ = 2j·(3Δs/r)^{3/2}'))
def knotendrehung(M, r, j):
    return 2 * j * (3 * ds(M) / r)**1.5


@formel('nullflaeche_takt', 'r_E(ϑ) = 3Δs + √(9Δs² − a_J²cos²ϑ)',
        ('r_E',), ('M', 'a_J', 'ϑ'), 'L', 'modell', 'Baustein 08', (O8, 'r_E(θ) = 3Δs + √(9Δs² − a_J²cos²θ)'))
def nullflaeche_takt(M, aJ, theta):
    """Nullfläche des Takts ruhender Uhren: am Äquator 6Δs, am Pol der Horizont."""
    d = ds(M)
    return 3 * d + math.sqrt(9 * d * d - (aJ * math.cos(theta))**2)


@formel('horizont_kerr', 'r_± = 3Δs ± √(9Δs² − a_J²); a_J² = r₊r₋, r₊ + r₋ = r₀ = 6Δs',
        ('r_±',), ('M', 'a_J'), 'L', 'modell', 'Baustein 08', (O8, 'a_J² = r₊·r₋,   r₊ + r₋ = r₀ = 6Δs'))
def horizont_kerr(M, aJ, aussen=True):
    d = ds(M)
    w = math.sqrt(9 * d * d - aJ * aJ)
    return 3 * d + (w if aussen else -w)


@formel('mitfuehrung', '|Ω_Mitführung| = (3Δs·a_J·c/r³)·|3(ĵ·r̂)r̂ − ĵ| = (3Δs·a_J·c/r³)·√(1 + 3cos²ψ)',
        ('Ω_LT',), ('M', 'a_J', 'r', 'cos ψ'), 'L', 'modell', 'Baustein 08', (O8, 'Ω_Mitführung = (3Δs·a_J·c/r³)·[3(ĵ·r̂)r̂ − ĵ]'))
def mitfuehrung(M, aJ, r, cos_psi):
    """Betrag der Mitführung eines Kreisels am Ort (Dipolanteil, universell), ψ gegen die Drehachse."""
    return 3 * ds(M) * aJ * c / r**3 * math.sqrt(1 + 3 * cos_psi**2)


@formel('mitfuehrung_polarbahn', 'Ω̄_Mitführung = ½·3Δs·a_J·c/r³ längs ĵ (Mittel über eine polare Kreisbahn)',
        ('Ω_LT',), ('M', 'a_J', 'r'), 'L', 'modell', 'Baustein 08', (O8, 'Gravity Probe B (polarer Erdorbit, 642 km Höhe)'),
        alternative='mitfuehrung')
def mitfuehrung_polarbahn(M, aJ, r):
    """Bahnmittel des ganzen Vektors (zeigt längs der Drehachse). Ein Kreisel, der auf einen Leitstern in der
    Bahnebene zeigt, misst davon nur den Anteil senkrecht zu seiner Achse: |Ω̄ × Ŝ| = Ω̄·cos δ_Leitstern."""
    return 1.5 * ds(M) * aJ * c / r**3


@formel('leitstern_anteil', 'Ω_LS = Ω_LT·cos δ_L', ('Ω_LS',), ('Ω_LT', 'δ_L'), 'L', 'modell', 'Baustein 08',
        (O8, 'vom Kreisel zum Leitstern IM Peg gesehen'))
def leitstern_anteil(omega, delta_l):
    """Anteil des Bahnmittels senkrecht zur Kreiselachse (Leitstern in der Bahnebene unter dem Winkel δ_L zur Achse)."""
    return omega * math.cos(delta_l)


@formel('geodaetisch', 'Ω_geodätisch = (9/2)·(Δs/r)·ω_Bahn, ω_Bahn = √(GM/r³)',
        ('Ω_geo',), ('M', 'r'), 'L', 'modell', 'Baustein 08', (O8, 'Ω_geodätisch = (9/2)·(Δs/r)·ω_Bahn'))
def geodaetisch(M, r):
    return 4.5 * ds(M) / r * math.sqrt(G * M / r**3)


@formel('quadrupol_leistung', 'P = (32/5)·(c⁵/G)·η_m²·(3Δs/r)⁵',
        ('P_GW',), ('M', 'η_m', 'r'), 'L', 'modell', 'Baustein 08', (O8, 'P = (32/5)·(c⁵/G)·η_m²·(3Δs/r)⁵'))
def quadrupol_leistung(M, eta, r):
    """Leistung einer Kreisbahn: fünfte Potenz des Taktunterschieds 3Δs/r."""
    return 32 / 5 * c**5 / G * eta**2 * (3 * ds(M) / r)**5


@formel('faltungsmesser_chirp', 'ḟ = (96/5)·π^{8/3}·(3Δs_c/c)^{5/3}·f^{11/3} ⇒ Δs_c = (c/3)·(ḟ/((96/5)π^{8/3}f^{11/3}))^{3/5}',
        ('Δs_c',), ('f_gw', 'ḟ_gw'), 'L', 'modell', 'Baustein 08', (O8, 'ḟ = (96/5)·π^{8/3}·(3Δs_c/c)^{5/3}·f^{11/3}'))
def faltungsmesser_chirp(f, fdot):
    """Faltung der Chirp-Masse aus Frequenz und Anstieg; im Detektor gemessen (rotverschoben)."""
    return c / 3 * (fdot / (96 / 5 * math.pi**(8 / 3) * f**(11 / 3)))**0.6


@formel('ringdown_frequenz', 'f = c/(2π·3Δs)·F(j)',
        ('f_RD',), ('M', 'F_j'), 'L', 'modell', 'Baustein 08', (O8, 'f = c/(2π·3Δs)·F(j)'))
def ringdown_frequenz(M, F):
    """Eigenschwingung des Trichters der Endmasse zwischen r₀ und dem flachen Außenraum; F hängt nur am Endspin."""
    return c / (2 * math.pi * 3 * ds(M)) * F


@formel('ringdown_abklingzeit', 'τ_RD = Q_RD/(πf_RD) = 6Δs·Q_RD/(c·F(j))',
        ('τ_RD',), ('M', 'F_j', 'Q_RD'), 'L', 'modell', 'Baustein 08', (O8, 'Ringdown GW150914: 249,9 Hz, τ = 4,03 ms'))
def ringdown_abklingzeit(M, F, Q):
    """Abklingzeit der Eigenschwingung des Trichters; Güte und Frequenz hängen nur am Endspin."""
    return Q / (math.pi * ringdown_frequenz(M, F))


@formel('richtungsfaktor_gedaechtnis', 'κ(ι) = sin²ι·(17 + cos²ι)/96 (Kante ι = π/2: 17/96)',
        ('κ_Θ',), ('ι',), 'L', 'modell', 'Baustein 08', (O8, 'κ(Θ) ist der Richtungsfaktor'))
def richtungsfaktor_gedaechtnis(iota):
    """Richtungsfaktor der bleibenden Verschiebung (führende Ordnung), ι Neigung der Bahnachse gegen die Sichtlinie."""
    s2 = math.sin(iota)**2
    return s2 * (17 + (1 - s2)) / 96


@formel('gedaechtnis_welle','h_mem = 12·κ(Θ)·Δs_rad/D_L, Δs_rad = G·E_rad/(3c⁴)',
        ('h_mem',), ('E_rad', 'D_L', 'κ_Θ'), 'L', 'modell', 'Baustein 08', (O8, 'h_mem = 12·κ(Θ)·Δs_rad/D_L'))
def gedaechtnis_welle(E_rad, D_L, kappa):
    return 12 * kappa * G * E_rad / (3 * c**4) / D_L
