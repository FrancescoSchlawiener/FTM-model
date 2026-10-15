# -*- coding: utf-8 -*-
"""
Teilchen und Sterne (02, 09; Ebene L, Faltungen am Ort): Compton-Länge, α_G, Sternmasse, Eddington, Carter,
Tiefe und Temperatur im Stern (Polytrope, Grenzmasse und Decke der Zählregel), Grenzmasse und Lebensdauer der Falte,
Leiter der Längen mit ihren Schnittmassen (Brücke L–I über a₀, Spiegel M ↔ m_Pl²/M), Raten und Kopplungen aus 09.
"""
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from .konstanten import *
from .formeln import formel
from .faltung import ds, r0, r_stern
from .rueckkehr import tau_D, ell0
from .daten import ame2020
O2 = "02-vom-teilchen-bis-zum-stern.md"
O4 = "04-zeit-verschraenkung-syntropie-gedaechtnis.md"
O6 = "06-torus-dynamik.md"
O8 = "08-licht-und-trichter.md"
O9 = "09-elektromagnetismus-und-frequenzleiter.md"


@formel('compton_red', 'λ̄ = ħ/(mc)',
        ('λ̄',), ('m',), 'L', 'modell', 'Baustein 02', (O2, '(λ̄ = ħ/mc)'))
def compton_red(m):
    return hbar / (m * c)


@formel('alpha_G', 'α_G = GM²/(ħc) = 3Δs/λ̄',
        ('α_G',), ('M',), 'L', 'modell', 'Baustein 02', (O2, 'α_G(M) = 3Δs/λ̄ = (M/m_Pl)²'))
def alpha_G(M=m_p):
    return G * M**2 / (hbar * c)


@formel('sternmasse', 'M_* = m_Pl³/m_p²',
        ('M_*',), ('m_Pl', 'm_p'), 'L', 'modell', 'Baustein 02', (O2, 'M_* = m_Pl³/m_p²'))
def sternmasse():
    return m_Pl**3 / m_p**2


@formel('grenzmasse_hawking', 'T_H = T ⇒ M_T = ħc³/(8πGk_BT), T_H = ħc/(24πk_BΔs)',
        ('M_T',), ('T',), 'L', 'modell', 'Baustein 04', ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'T = ħc/(24π·k·Δs)'))
def grenzmasse_hawking(T=T_CMB):
    """Falte mit T_H = T: M = ħc³/(8πGkT)."""
    return hbar * c**3 / (8 * math.pi * G * kB * T)


@formel('masse_r_stern_gleich_compton', 'r* = λ̄ bei M_c = (ħ²a₀/(Gc²))^{1/3} = ε_L^{1/3}·m_Pl',
        ('M_c',), ('a₀',), 'Z', 'modell', 'Baustein 08', (O8, 'M_c = ε_L^{1/3}·m_Pl'))
def masse_r_stern_gleich_compton(a0=A0):
    """r* = λ̄ bei M = (ħ²a₀/(Gc²))^{1/3}."""
    return (hbar**2 * a0 / (G * c**2)) ** (1 / 3)


@formel('eddington', 'P_Edd = (9/2)·m_p·c³·Δs/r_e²',
        ('P_Edd',), ('M',), 'L', 'modell', 'Baustein 02', (O2, 'L_Edd = (9/2)·m_p·c³·Δs/r_e²'))
def eddington(M):
    """L_Edd = (9/2)·m_p·c³·Δs/r_e²."""
    return 4.5 * m_p * c**3 * ds(M) / r_e**2


@formel('entropie_falte', 'S_F = A/(4ℓ_P²) = π·r₀²/ℓ_P² = 36π·(Δs/ℓ_P)² = 4π·α_G(M) (in k_B)',
        ('S_F',), ('M',), 'L', 'modell', 'Baustein 04', (O4, 'S = 36π·k·(Δs/ℓ_P)²'))
def entropie_falte(M):
    """Entropie der Falte: Fläche der Kugel r₀ = 6Δs in Einheiten 4ℓ_P²."""
    return math.pi * r0(M)**2 / l_P**2


@formel('lebensdauer_falte', 't_F = 3840·S_F·Δs/c',
        ('t_F',), ('M',), 'L', 'modell', 'Baustein 04', (O4, 'Lebensdauer t = 3840·(S/k)·Δs/c'))
def lebensdauer_falte(M):
    """t = 3840·(S/k)·Δs/c mit der Entropie der Falte."""
    return 3840 * entropie_falte(M) * ds(M) / c


@formel('rueckgabe_rate', 'dM/dt = −ħc⁴/(15360π·G²·M²)',
        ('Ṁ_F',), ('M',), 'L', 'modell', 'Baustein 04',
        ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'Rückgabe dM/dt = −ħc⁴/(15360π·G²·M²)'))
def rueckgabe_rate(M):
    """Rückgabe der Falte an den Raum (ohne Atmen); integriert über die Masse gibt sie die Lebensdauer."""
    return -hbar * c**4 / (15360 * math.pi * G**2 * M**2)


@formel('page_zeit', 't_P/t_F = 1 − (β/(1 + β))^{3/2}',
        ('t_P/t_F',), ('β_Page',), 'Z', 'modell', 'Baustein 04', (O4, 'Page-Zeit 0,54 der Lebensdauer'))
def page_zeit(beta):
    """Die Strahlung trägt β-mal die Entropie, die die Falte verliert (S_rad = β(S₀ − S_F)); mit S_F/S₀ = (1 − t/t_F)^{2/3}
    ist die Page-Zeit der Gleichstand S_rad = S_F."""
    return 1 - (beta / (1 + beta))**1.5


@formel('rueckverdichtung_teile', 'N_RV = β·z_f²: (N/z_f)²S₀ ≥ N·β·S₀',
        ('N_RV',), ('z_f', 'β_Page'), 'Z', 'modell', 'Baustein 04', (O4, 'N ≥ 1,5·z_f²'))
def rueckverdichtung_teile(z_f, beta):
    """N Falten strahlen N·β·S₀ ab; die um z_f verdünnte Energie gibt eine Falte mit (N/z_f)²S₀. Zweiter Hauptsatz:
    nur hierarchisch aufwärts."""
    return beta * z_f**2


@formel('leistung_absorbiert', 'P_abs = P_L/(4πd²)·πR²·(1 − A_alb)',
        ('P_abs',), ('P_L☉', 'r', 'R', 'A_alb'), 'L', 'hilfe', 'Baustein 04')
def leistung_absorbiert(L, d, R, albedo):
    """Leistung, die eine Kugel vom Radius R im Abstand d einer Quelle der Leuchtkraft L aufnimmt."""
    return L / (4 * math.pi * d * d) * math.pi * R * R * (1 - albedo)


@formel('temperatur_abstrahlung', 'T_A = (P_abs/(4πR²σ_SB))^{1/4}, σ_SB = π²k_B⁴/(60ħ³c²)',
        ('T_A',), ('P_abs', 'R'), 'L', 'hilfe', 'Baustein 04')
def temperatur_abstrahlung(P, R):
    """Temperatur, bei der eine Kugel die aufgenommene Leistung über ihre ganze Fläche wieder abstrahlt."""
    sigma_sb = math.pi**2 * kB**4 / (60 * hbar**3 * c**2)
    return (P / (4 * math.pi * R * R * sigma_sb))**0.25


@formel('entropieexport', 'J_S = (4/3)·P_abs·(1/T_A − 1/T)',
        ('J_S',), ('P_abs', 'T_A', 'T'), 'L', 'modell', 'Baustein 04', (O4, 'Die Erde exportiert'))
def entropieexport(P, T_aus, T_ein):
    """Entropiestrom der Strahlung (Faktor 4/3) bei Aufnahme mit T_ein und Abgabe mit T_aus."""
    return 4 / 3 * P * (1 / T_aus - 1 / T_ein)


# ------------------------------------------------------------------ Starke Kopplung und Ladung (02)
@formel('alpha_s_lauf', 'dα_s/d ln μ_R = −(β₀/2π)α_s² − (β₁/8π²)α_s³, β₀ = 11 − 2n_f/3, β₁ = 102 − 38n_f/3; '
        'n_f = 3 + Zahl der Schwellen m_q unter μ_R, α_s stetig',
        ('α_s',), ('α_s', 'μ_R', 'm_q'), 'L', 'hilfe', 'Baustein 02')
def alpha_s_lauf(alpha0, mu0, mu1, schwellen=(), schleifen=2):
    """α_s(μ₁) aus α_s(μ₀) (Skalen in GeV, auf- oder abwärts), ein oder zwei Schleifen, je Abschnitt zwischen den
    Schwellen mit festem n_f; Fehler, wenn die Kopplung über 5 wächst."""
    if not (0 < alpha0 < 5 and mu0 > 0 and mu1 > 0):
        raise ValueError("α_s oder Skala außerhalb des Bereichs")
    L0, L1 = math.log(mu0), math.log(mu1)
    kanten = sorted({L0, L1} | {math.log(m) for m in schwellen if min(mu0, mu1) < m < max(mu0, mu1)}, reverse=L1 < L0)
    a = alpha0

    def gross(L, y):
        return y[0] - 5.0
    gross.terminal = True
    for A, B in zip(kanten[:-1], kanten[1:]):
        nf = 3 + sum(1 for m in schwellen if math.log(m) < (A + B) / 2)
        b0, b1 = 11 - 2 * nf / 3, 102 - 38 * nf / 3
        s = solve_ivp(lambda L, y: [-(b0 / (2 * math.pi)) * y[0]**2 - (schleifen == 2) * (b1 / (8 * math.pi**2)) * y[0]**3],
                      (A, B), [a], events=gross, rtol=1e-11, atol=1e-14)
        if s.status != 0:
            raise ValueError("Kopplung divergiert")
        a = s.y[0][-1]
    return a


@formel('skala_bei_kopplung', 'μ_R mit α_s(μ_R) = α_Ziel, abwärts gelaufen von (α_s, μ₀)',
        ('μ_R',), ('α_s', 'm_q'), 'L', 'hilfe', 'Baustein 02')
def skala_bei_kopplung(alpha0, mu0, alpha_ziel, schwellen=(), schleifen=2):
    """Skala (GeV), bei der die abwärts laufende Kopplung den Wert α₀ < α_Ziel < 5 erreicht (Nullstelle in ln μ)."""
    if not (0 < alpha0 < alpha_ziel < 5 and mu0 > 0):
        raise ValueError("α_Ziel muss über α₀ liegen")

    def f(L):
        try:
            return alpha_s_lauf(alpha0, mu0, math.exp(L), schwellen, schleifen) - alpha_ziel
        except ValueError:
            return 1.0     # jenseits des Pols: Kopplung über jedem Ziel
    oben, unten = math.log(mu0), math.log(mu0) - 1.0
    while f(unten) < 0:
        oben, unten = unten, unten - 2.0 * (oben - unten)
        if unten < math.log(mu0) - 1e3:
            raise ValueError("α_Ziel nicht erreicht")
    return math.exp(brentq(f, unten, oben, xtol=1e-13, rtol=1e-14))


@formel('nome_qcd', 'q_QCD = Λ_QCD/M = e^{−2π/(β₀α_s(M))}, β₀ = 11 − 2n_f/3',
        ('q_QCD',), ('α_s', 'n_f'), 'L', 'modell', 'Baustein 02', (O2, 'Λ_QCD/M = e^{−2π/(β₀α_s)}'))
def nome_qcd(alpha_s, nf):
    """Nome der starken Kopplung, Im τ̃ = 1/(β₀α_s): Die Protonmasse folgt aus Λ_QCD = M·|q̃| (Lösung einer Schleife)."""
    return math.exp(-2 * math.pi / ((11 - 2 * nf / 3) * alpha_s))


@formel('dirac_verhaeltnis', 'g_D/e = 1/(2α) (Dirac-Bedingung e·g = ħc/2, Gauß)',
        ('g_D/e',), ('α',), 'L', 'modell', 'Baustein 02', (O2, 'Verhältnis magnetischer zu elektrischer Elementarladung 1/(2α)'))
def dirac_verhaeltnis():
    """Kleinste magnetische Ladung gegen die elektrische: Die Windung des Flusses ist ganzzahlig (c₁ ∈ ℤ)."""
    return 1 / (2 * alpha)


@formel('radius_ladungskreis', 'r_Qel = 2ℓ_P/√α (Ladung e als Windung um einen Kreis)',
        ('r_Qel',), ('α',), 'L', 'modell', 'Baustein 02', (O2, 'Sein Radius ist 2ℓ_P/√α'))
def radius_ladungskreis():
    return 2 * l_P / math.sqrt(alpha)


@formel('energie_ladungskreis', 'E_Qel = ħc/r_Qel = (√α/2)·m_Plc²',
        ('E_Qel',), ('r_Qel',), 'L', 'modell', 'Baustein 02', (O2, 'ein Umlauf trüge'))
def energie_ladungskreis():
    """Energie der ersten Umlaufmode des Ladungskreises; weit über m_ec², die Elektronmasse stammt nicht daher."""
    return hbar * c / radius_ladungskreis()


# ------------------------------------------------------------------ Zustand des Stoffs und Oberfläche der Sterne (02)
@formel('zustandsschwelle', 'δ_Z = 2k_BT/(m_Ac²) = 2(v_th/c)²',
        ('δ_Z',), ('T', 'm_A'), 'L', 'modell', 'Baustein 02', (O2, '2kT/(m_A c²) = 2(v_th/c)²'))
def zustandsschwelle(T, m_A):
    """Tiefe, ab der ein Trichter den Zustand seines Stoffs bestimmt (derselbe Quotient wie im Jeans-Kriterium)."""
    return 2 * kB * T / (m_A * c**2)


@formel('f_saha', 'f_s = k_BT_eff/(α²m_ec²)',
        ('f_s',), ('T', 'α', 'm_e'), 'L', 'modell', 'Baustein 02', (O2, 'f_s = kT_eff/(α²m_ec²)'))
def f_saha(T_eff):
    """Oberflächentemperatur in Einheiten der atomaren Bindung α²m_ec²."""
    return kB * T_eff / (alpha**2 * m_e * c**2)


@formel('eddington_quartik', 'M_β = C_E·(ħc/G)^{3/2}/(μ_m m_u)²·√(1 − β_E)/β_E², C_E = 4m₁(3)·√45/π^{3/2}',
        ('M_β',), ('β_E', 'μ_m'), 'L', 'modell', 'Baustein 02', (O2, 'Strahlungsdruck und Gasdruck sind'))
def eddington_quartik(beta, mu):
    """Eddington-Standardmodell (Polytrope n = 3, β_E fest): Masse mit Gasdruckanteil β_E; β_E = ½ gibt den
    Gleichstand von Strahlungs- und Gasdruck."""
    C_E = 4 * lane_emden(3.0)[1] * math.sqrt(45) / math.pi**1.5
    return C_E * (hbar * c / G)**1.5 / (mu * m_u)**2 * math.sqrt(1 - beta) / beta**2


# ------------------------------------------------------------------ Kerne (06)
_TR = {}


@formel('troepfchen_anpassung', 'B = a_Vol·A − a_Obf·A^{2/3} − a_Cou·Z(Z − 1)/A^{1/3} − a_Asy·(N − Z)²/A + a_Paar·δ_P/√A, '
        'kleinste Quadrate an die gemessenen Kerne aus AME2020 mit A ≥ A_min',
        ('a_Vol', 'a_Obf', 'a_Cou', 'a_Asy', 'a_Paar'), ('A_A',), 'L', 'hilfe', 'AME2020')
def troepfchen_anpassung(A_min=16):
    """Koeffizienten (J) des Tröpfchenmodells; δ_P = +1 gerade-gerade, −1 ungerade-ungerade, sonst 0."""
    if A_min not in _TR:
        zeilen, B = [], []
        for Z in range(0, 121):
            for A in range(max(int(math.ceil(A_min)), 1), 301):
                try:
                    k = ame2020(Z, A)
                except KeyError:
                    continue
                if not k["gemessen"] or not np.isfinite(k["BA"]):
                    continue
                N = A - Z
                dP = 1.0 if Z % 2 == 0 and N % 2 == 0 else (-1.0 if Z % 2 and N % 2 else 0.0)
                zeilen.append([A, -A**(2 / 3), -Z * (Z - 1) / A**(1 / 3), -(N - Z)**2 / A, dP / math.sqrt(A)])
                B.append(k["BA"] * A)
        if len(B) < 5:
            raise ValueError("zu wenige Kerne")
        _TR[A_min] = tuple(np.linalg.lstsq(np.array(zeilen), np.array(B), rcond=None)[0])
    return _TR[A_min]


@formel('spaltbarkeit', 'x_Sp = E_C/(2E_S) = a_Cou·Z²/(2a_Obf·A)',
        ('x_Sp',), ('Z_A', 'A_A', 'a_Cou', 'a_Obf'), 'L', 'modell', 'Baustein 06', (O6, 'Spaltbarkeit x = E_C/(2E_S)'))
def spaltbarkeit(Z, A, a_C, a_S):
    """Verstärkung der Kernverformung: Rückstellung ∝ 1 − x, Kippstelle bei x = 1."""
    return a_C * Z * Z / (2 * a_S * A)


# ------------------------------------------------------------------ Kohärenz (09)
@formel('debye_laenge', 'λ_D = √(ε₀ε_rk_BT/(2N_Ae²I_ion))',
        ('λ_D',), ('T', 'ε_r', 'I_ion'), 'L', 'modell', 'Baustein 09', (O9, 'die Debye-Länge in physiologischer Salzlösung'))
def debye_laenge(T, eps_r, I):
    """Abschirmlänge elektrischer Felder in einer Lösung der Ionenstärke I (mol/m³)."""
    return math.sqrt(eps0 * eps_r * kB * T / (2 * N_A * e_el**2 * I))


# ------------------------------------------------------------------ Carter (02, 09)
@formel('carter_faktor', '𝒞 = α_G(m_p)/(α¹²(m_e/m_p)⁴)',
        ('𝒞',), ('α_G', 'α', 'm_e', 'm_p'), 'Z', 'modell', 'Baustein 02', (O2, 'Der Faktor 𝒞 = 2,944'))
def carter_faktor():
    """Carters Beziehung α_G ∝ α¹²(m_e/m_p)⁴; der Faktor ist eine Identität der Sternzahlen, keine Gleichung."""
    return alpha_G(m_p) / (alpha**12 * (m_e / m_p)**4)


# ------------------------------------------------------------------ Tiefe und Temperatur im Stern (02)
_LE = {}


@formel('lane_emden', 'θ″ + 2θ′/x + θⁿ = 0, θ(0) = 1: gibt (x₁, −x₁²θ′(x₁), ρ_c/ρ̄ = x₁³/(3m₁))',
        ('x₁_LE', 'm₁_LE', 'ρ_c/ρ̄'), ('n_P',), 'Z', 'hilfe', 'Baustein 02', None)
def lane_emden(n):
    """Polytrope vom Index n (0 ≤ n < 5), Integration bis zur ersten Nullstelle; zwischengespeichert je n."""
    if n not in _LE:
        def rhs(x, y):
            return [y[1], -max(y[0], 0.0) ** n - 2 * y[1] / x]

        def rand(x, y):
            return y[0]
        rand.terminal, rand.direction = True, -1
        x0 = 1e-6
        s = solve_ivp(rhs, [x0, 40], [1 - x0**2 / 6, -x0 / 3], events=rand, rtol=1e-11, atol=1e-13)
        if not len(s.t_events[0]):
            raise ValueError(f"Polytrope n = {n}: keine Nullstelle (n ≥ 5)")
        x1, d1 = s.t_events[0][0], s.y_events[0][0][1]
        m1 = -x1**2 * d1
        _LE[n] = (x1, m1, x1**3 / (3 * m1))
    return _LE[n]


@formel('tiefe_faktor', 'f = θ_n = x₁/((n + 1)·m₁), kT_c = f·δ·m̄c²',
        ('f_T',), ('n_P',), 'Z', 'modell', 'Baustein 02', (O2, 'f = 0,54 (konvektiv, n = 1,5),  0,85 (strahlend, n = 3)'))
def tiefe_faktor(n=3.0):
    """Anteil der Tiefe: Die Wärme je Teilchen im Zentrum ist der innere Taktunterschied m̄c²·Δδ_innen, geteilt durch
    n + 1 (Enthalpie des idealen Gases); Δδ_innen = (n + 1)·f·δ."""
    x1, m1, _ = lane_emden(n)
    return x1 / ((n + 1) * m1)


@formel('temperatur_kern', 'T_c = f·δ·μ_m·m_u·c²/k_B, δ = 3Δs/R',
        ('T_c',), ('M', 'R', 'μ_m', 'n_P'), 'L', 'modell', 'Baustein 02', (O2, 'kT_c = m̄c²·Δδ_innen/(n + 1) = f·δ·m̄c²'))
def temperatur_kern(M, R, mu, n=3.0):
    """Kerntemperatur eines nicht entarteten Trichters aus seiner Tiefe δ = 3Δs/R, am Ort gemessen (lokale Uhr)."""
    return tiefe_faktor(n) * 3 * ds(M) / R * mu * m_u * c**2 / kB


@formel('druck_elektronen', 'Π(x) = P_e/(n_e m_ec²) = [x(2x² − 3)√(1 + x²) + 3 arsinh x]/(8x³), x = λ̄_e·k_F',
        ('Π_e',), ('x_F',), 'L', 'hilfe', 'Baustein 02', None)
def druck_elektronen(x):
    """Zählregel der Elektronen (kaltes Gas, exakt): Π = x²/5 langsam, x/4 hochrelativistisch."""
    if x < 1e-3:
        return x * x / 5 - x**4 / 14
    fx = x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x)
    return fx / (8 * x**3)


@formel('grenzmasse_zaehlregel', 'M_Grenze = (K_ur·(3ρ_c/(4πρ̄))^{1/3}/(f·G))^{3/2}, K_ur = (3π²)^{1/3}ħc/(4(μ_e m_u)^{4/3})',
        ('M_Grenze',), ('μ_e', 'n_P'), 'L', 'modell', 'Baustein 02', (O2, 'ohne Wärme gibt dieselbe Bilanz die Grenzmasse 1,456 M☉'))
def grenzmasse_zaehlregel(mu_e=2.0, n=3.0):
    """Masse, bei der Tiefe und hochrelativistische Zählregel der Elektronen gleich mit 1/R laufen: R fällt heraus."""
    K_ur = (3 * math.pi**2) ** (1 / 3) * hbar * c / 4 / (mu_e * m_u) ** (4 / 3)
    D = lane_emden(n)[2]
    return (K_ur * (3 * D / (4 * math.pi)) ** (1 / 3) / (tiefe_faktor(n) * G)) ** 1.5


@formel('temperatur_max', 'kT_max = max_R μ_m m_u c²·[f·δ − (m_e/(μ_e m_u))·Π(x)], x = λ̄_e(3π²n_e)^{1/3}',
        ('T_max',), ('M', 'μ_m', 'μ_e', 'n_P'), 'L', 'modell', 'Baustein 02', (O2, 'kT_max = (μ/μ_e)·m_ec²·g(M/M_Grenze)'))
def temperatur_max(M, mu, mu_e=2.0, n=3.0):
    """Höchste Kerntemperatur beim Zusammenziehen: Die Tiefe wächst wie 1/R, die Zählregel wie 1/R² (Π = x²/5)."""
    f, D = tiefe_faktor(n), lane_emden(n)[2]
    lam_e = hbar / (m_e * c)

    def kT(lgR):
        R = 10.0**lgR
        n_e = D * 3 * M / (4 * math.pi * R**3) / (mu_e * m_u)
        x = lam_e * (3 * math.pi**2 * n_e) ** (1 / 3)
        return mu * m_u * c**2 * (f * G * M / (R * c**2) - m_e / (mu_e * m_u) * druck_elektronen(x))
    r = minimize_scalar(lambda lgR: -kT(lgR), bounds=(1.0, 13.0), method="bounded", options=dict(xatol=1e-9))
    return kT(r.x) / kB


@formel('decke_faktor', 'kT_max → (μ_m/μ_e)·m_ec²·max_x[x/4 − Π(x)] an der Grenzmasse',
        ('g_T',), (), 'Z', 'modell', 'Baustein 02', (O2, 'kT ≈ 0,1·(μ/μ_e)·m_ec²'))
def decke_faktor():
    """Decke der Temperatur eines Trichters, den die Zählregel der Elektronen hält (Gegenlänge λ̄_e)."""
    r = minimize_scalar(lambda x: -(x / 4 - druck_elektronen(x)), bounds=(0.1, 10.0), method="bounded", options=dict(xatol=1e-10))
    return r.x / 4 - druck_elektronen(r.x)


# ------------------------------------------------------------------ Leiter der Längen und Schnittmassen (01, 08)
@formel('leiterzahl', 'ε_L = ℓ_P/ℓ₀ = a₀/a_P, a_P = c²/ℓ_P',
        ('ε_L',), ('a₀',), 'Z', 'modell', 'Baustein 08', (O8, 'ε_L = ℓ_P/ℓ₀ = a₀/a_P'))
def leiterzahl(a0=A0):
    return l_P / ell0(a0)


@formel('schnittmasse_unten', 'M_min = ε_L·m_Pl = ħ/(c²τ_D): λ̄ = ℓ₀ und r* = ℓ_P',
        ('M_min',), ('a₀',), 'Z', 'modell', 'Baustein 08', (O8, 'M_min = ε_L·m_Pl = ħ/(c²τ_D)'))
def schnittmasse_unten(a0=A0):
    return hbar / (c**2 * tau_D(a0))


@formel('schnittmasse_oben', 'M_U = m_Pl/ε_L = c³τ_D/G: r* = 3Δs = ℓ₀',
        ('M_U',), ('a₀',), 'Z', 'modell', 'Baustein 08', (O8, 'M_U = m_Pl/ε_L = c³τ_D/G'))
def schnittmasse_oben(a0=A0):
    return c**3 * tau_D(a0) / G


@formel('schwelle_zu_huelle', 'r*/R = √(GM/(a₀R²)) = √(g(R)/a₀); r*/R = 1 ⇔ Σ̄ = a₀/(πG)',
        ('r*/R',), ('M', 'R', 'a₀'), 'L', 'modell', 'Baustein 08', (O8, 'Nach der Lage von r* zur Hülle gibt es zwei Klassen'))
def schwelle_zu_huelle(M, R, a0=A0):
    """Klasse: r*/R ≫ 1 kompakt (Schwelle weit außen), r*/R ≲ 1 kollektiv (als Ganzes unter der Schwelle)."""
    return r_stern(M, a0) / R


@formel('spiegel_masse', 'M̃ = m_Pl²/M: vertauscht λ̄ und 3Δs, bildet M_min auf M_U ab',
        ('M̃',), ('M',), 'L', 'modell', 'Baustein 08', (O8, 'Die Spiegelung M ↔ m_Pl²/M vertauscht λ̄ und 3Δs'))
def spiegel_masse(M):
    return m_Pl**2 / M


@formel('spiegel_radius', 'r ↦ r*²/r: bildet r₀ = 6Δs auf Λ ab (Raum auf innere Richtung)',
        ('r_sp',), ('r', 'M', 'a₀'), 'LI', 'modell', 'Autor (FTM-PDF)', ('08-licht-und-trichter.md', 'r → r*²/r bildet den innersten Ort r₀ auf die Schließlänge Λ ab'))
def spiegel_radius(r, M, a0=A0):
    """Symmetrie der Leiter um r*; die Regel selbst hat sie nicht."""
    return r_stern(M, a0)**2 / r


# ------------------------------------------------------------------ Raten und Kopplungen (09)
@formel('verhaeltnis_em_grav', 'α/α_G(m)',
        ('α/α_G',), ('α', 'm'), 'L', 'modell', 'Baustein 09', (O9, 'α/α_G = 1,23555·10³⁶'))
def verhaeltnis_em_grav(m=m_p):
    return alpha / alpha_G(m)


@formel('compton_takt', 'ω_C = c/λ̄ = Mc²/ħ',
        ('ω_C',), ('M',), 'L', 'modell', 'Baustein 09', (O9, 'Compton-Takt ω_C = c/λ̄'))
def compton_takt(M):
    return c / compton_red(M)


@formel('faltungsrate', 'Γ_G = α_G(M)·ω_C = GM³c/ħ²',
        ('Γ_G',), ('M',), 'L', 'modell', 'Baustein 09', (O9, 'Faltungsrate Γ_G = α_G(M)·ω_C'))
def faltungsrate(M):
    return alpha_G(M) * compton_takt(M)


@formel('faltungstakt', 'ω_F = c/(3Δs)',
        ('ω_F',), ('M',), 'L', 'modell', 'Baustein 09', (O9, 'Faltungstakt c/(3Δs)'))
def faltungstakt(M):
    """Takt der Faltung im Zentrum (keine Raumlänge)."""
    return c / (3 * ds(M))


@formel('umlauf_wasserstoff', 'ω₂ = αc/a_B',
        ('ω_H',), ('α', 'a_B'), 'L', 'modell', 'Baustein 09', (O9, 'ω₂ = αc/a_B'))
def umlauf_wasserstoff():
    return alpha * c / a_B


@formel('umlauf_kepler', 'ω = √(GM/r³) (ferne Uhr)',
        ('Ω_N',), ('M', 'r'), 'L', 'modell', 'Baustein 09', (O9, '√(GM/R³)'))
def umlauf_kepler(M, r):
    return math.sqrt(G * M / r**3)


@formel('spinzahl', 'j = a_J/(3Δs) = (J/ħ)/α_G(M)',
        ('j',), ('a_J', 'M'), 'L', 'modell', 'Baustein 09', (O9, 'Die Spinzahl a_J/(3Δs) = (J/ħ)/α_G(M)'))
def spinzahl(aJ, M):
    return aJ / (3 * ds(M))


@formel('moment_magnetisch', 'μ = Q_el·c·a_J (g = 2 für jede Masse, Ladung, Drehung)',
        ('μ_mag',), ('Q_el', 'a_J'), 'L', 'modell', 'Baustein 09', (O9, 'μ = Q_el·c·a_J'))
def moment_magnetisch(Q, aJ):
    """Faltung und Windung am selben Knoten, beide Multipolreihen mit dem Versatz i·a_J."""
    return Q * c * aJ


@formel('einrast_phase', 'dΔθ/dt = Δω − K sin Δθ = 0 ⇒ Δθ = arcsin(Δω/K) für K ≥ |Δω|',
        ('Δθ',), ('Δω', 'K_E'), 'Z', 'modell', 'Baustein 09', (O9, 'Δθ = arcsin(Δω/K)'))
def einrast_phase(dw, K):
    """Phasenlage zweier eingerasteter Zyklen; für K < |Δω| laufen sie quasi-periodisch weiter (Fehler)."""
    if K < abs(dw):
        raise ValueError("kein Einrasten: K < |Δω|")
    return math.asin(dw / K)


@formel('dephasierzeit', 't_φ = ħ/(k_BT)',
        ('t_φ',), ('T',), 'L', 'modell', 'Baustein 09', (O9, 'ħ/kT = 2,46·10⁻¹⁴ s'))
def dephasierzeit(T):
    return hbar / (kB * T)


@formel('schleife_thermisch', 'q = (∂P/∂T)/G_th',
        ('q',), ('∂P/∂T', 'G_th'), 'L', 'modell', 'Baustein 09', (O9, 'q = (∂P/∂T)/G'), alternative='q_',
        gleiche_form='kac_rueckkehrzeit')
def schleife_thermisch(dP_dT, G_th):
    """Schleifenverstärkung der thermischen Rückwirkung: derselbe Nome q als Regelkreis am Ort; Vorzeichen aus Antrieb
    und Richtung von R(T), bei q ≥ 1 thermisches Durchgehen."""
    return dP_dT / G_th
