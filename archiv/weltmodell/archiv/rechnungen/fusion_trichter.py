# -*- coding: utf-8 -*-
"""
Fusion im Trichter: Temperatur, Zündung, Gedächtnis, Elemente als Spur der Tiefe.
Frage: Ist die Temperatur ein Gedächtnis, gesteuert von der Tiefe δ = 3Δs/R des Trichters?
Rechnet nur mit ftm_kern (Δs, δ, Konstanten). Daten: AME2020 (daten/elemente/mass_1.mas20.txt).
[G] = Messwert aus dem Gedächtnis, nicht nachgeprüft;  [A] = Annahme.
Beobachter: T ist lokal gemessen (am Ort, eigene Uhr); δ ist der Taktunterschied gegen den flachen Außenraum (ferner Beobachter).
Teile: 1 Tiefe → Temperatur   2 Zündschwellen als Tiefen, Mindestmasse   3 Gedächtnis (Zeiten)   4 Elemente   6 Gegenproben
"""
import os, sys, io, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

HERE = os.path.dirname(os.path.abspath(__file__))
DAT = os.path.join(HERE, "..", "daten", "elemente")
out = []
def P(s=""):
    out.append(s)

m_u = 1.66053906660e-27          # atomare Masseneinheit (CODATA 2018), fehlt im Kern
eV = e_el
keV, MeV = 1e3 * eV, 1e6 * eV
muc2 = m_u * c**2                # 931,494 MeV
lam_e = compton_red(m_e)
R_J, M_J = 7.1492e7, 1.89813e27  # Jupiter [G]

P("FUSION IM TRICHTER: Temperatur, Zündung, Gedächtnis, Elemente")
P("Tiefe δ = 3Δs/R (Taktunterschied Hülle gegen flachen Außenraum), kT lokal gemessen, m̄ = μ·m_u")
P()

# ======================================================================= 1 Tiefe → Temperatur
P("1) TIEFE → TEMPERATUR")
P("   Herleitung im Trichter: Längs eines ruhenden Trichters gilt am Ort dh = T ds − dΦ (Hydrostatik, h = Enthalpie je Teilchen).")
P("   Ohne Entropiegefälle (Polytrope, Konvektion) ist h + m̄Φ fest: Die Wärme je Teilchen im Zentrum ist der Taktunterschied")
P("   zwischen Zentrum und Hülle in Ruheenergie,  h_c = m̄c²·Δδ_innen,  Δδ_innen = (Φ(R) − Φ(0))/c².  Beim idealen Gas h = (n+1)kT:")
P("        kT_c = m̄c²·Δδ_innen/(n+1) = f·δ·m̄c²,    f = θ_n = Δδ_innen/((n+1)δ)")
P()

def lane_emden(n):
    def rhs(x, y):
        th, dth = y
        return [dth, -max(th, 0)**n - 2 * dth / x]
    x0 = 1e-6
    ev = lambda x, y: y[0]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(rhs, [x0, 20], [1 - x0**2 / 6, -x0 / 3], events=ev, rtol=1e-11, atol=1e-13, dense_output=True)
    x1 = s.t_events[0][0]; d1 = s.y_events[0][0][1]
    return x1, -x1**2 * d1, s

LE = {}
for n in (1.5, 3.0):
    x1, mw, sol = lane_emden(n)
    theta = x1 / ((n + 1) * mw)               # kT_c = θ_n·GMm̄/R
    D = x1**3 / (3 * mw)                      # ρ_c/ρ̄
    # zweiter Weg: Potential aus dem Dichteprofil integrieren, Φ(R) − Φ(0) = ∫ G m(r)/r² dr
    xs = np.linspace(1e-6, x1, 20001)
    th = np.clip(sol.sol(xs)[0], 0, None)
    mx = -xs**2 * sol.sol(xs)[1]              # m(ξ) in Einheiten 4πα³ρ_c
    dPhi = np.trapezoid(mx / xs**2, xs)       # in Einheiten 4πGα²ρ_c
    # GM/R in denselben Einheiten: M/R = 4πα³ρ_c·mw/(α x1) → mw/x1
    theta2 = dPhi / (mw / x1) / (n + 1)
    LE[n] = dict(x1=x1, mw=mw, theta=theta, D=D)
    P(f"   Polytrope n = {n}: ξ₁ = {x1:.4f}, −ξ₁²θ'(ξ₁) = {mw:.4f}, ρ_c/ρ̄ = {D:.2f};  f = θ_n = {theta:.4f}"
      f"  (Potentialintegral: {theta2:.4f});  Δδ_innen = {(n+1)*theta:.3f}·δ")
P("   [G] Tafelwerte: ξ₁ = 3,6538 / 6,8969; −ξ₁²θ' = 2,7141 / 2,0182")
P()

def beta_edd(M, mu):
    """Eddington-Quartik: 1 − β = 0,003·(M/M☉)²·μ⁴·β⁴ (Standardmodell n = 3)."""
    return brentq(lambda b: 1 - b - 0.002979 * (M / M_sun)**2 * mu**4 * b**4, 1e-6, 1)

# Körper: (Name, M, R, μ, n, T_c[G], Bemerkung)
koerper = [
    ("Sonne (μ homogen 0,62)", M_sun, R_sun, 0.62, 3.0, 1.57e7, "Hauptreihe"),
    ("Sonne (μ_c heute 0,85 [G])", M_sun, R_sun, 0.85, 3.0, 1.57e7, "Hauptreihe"),
    ("Roter Zwerg 0,3 M☉, 0,29 R☉ [G]", 0.3 * M_sun, 0.29 * R_sun, 0.62, 1.5, 7.5e6, "voll konvektiv, T_c [G ±30 %]"),
    ("Roter Zwerg 0,1 M☉, 0,12 R☉ [G]", 0.1 * M_sun, 0.12 * R_sun, 0.62, 1.5, 4.5e6, "teilentartet, T_c [G ±30 %]"),
    ("massiv 15 M☉, 4,7 R☉ ZAMS [G]", 15 * M_sun, 4.7 * R_sun, 0.62, 3.0, 3.4e7, "Strahlungsdruck β"),
    ("massiv 60 M☉, 10 R☉ ZAMS [G]", 60 * M_sun, 10 * R_sun, 0.62, 3.0, 4.0e7, "Strahlungsdruck β"),
    ("Brauner Zwerg 0,05 M☉, 0,1 R☉ (1 Gyr) [G]", 0.05 * M_sun, 0.1 * R_sun, 0.62, 1.5, 1.0e6, "entartet"),
    ("Jupiter", M_J, R_J, 0.62, 1.0, 2.0e4, "entartet/flüssig, n≈1"),
    ("Weißer Zwerg 0,6 M☉, 0,0125 R☉ [G]", 0.6 * M_sun, 0.0125 * R_sun, 14 / 7.5, 1.5, 1.0e7, "entartet, C/O"),
    ("Sirius B 1,02 M☉, 0,0084 R☉ [G]", 1.02 * M_sun, 0.0084 * R_sun, 14 / 7.5, 3.0, 2.0e7, "entartet"),
]
LE[1.0] = dict(x1=math.pi, mw=math.pi, theta=math.pi / (2 * math.pi), D=math.pi**3 / (3 * math.pi))
P(f"   {'Körper':44s} {'δ=3Δs/R':>10s} {'f':>6s} {'β':>5s} {'T_c Tiefe':>10s} {'T_c [G]':>9s} {'T[G]/T_Tiefe':>12s}")
tab1 = []
for name, M, R, mu, n, TG, bem in koerper:
    d = 3 * ds(M) / R
    f = LE[n]["theta"]
    b = beta_edd(M, mu) if "Strahl" in bem else 1.0
    T = f * b * d * mu * muc2 / kB
    tab1.append((name, d, T, TG))
    P(f"   {name:44s} {d:10.3e} {f:6.3f} {b:5.3f} {T:10.3e} {TG:9.2e} {TG/T:12.3f}   {bem}")
P("   Hauptreihe: T_c = f·δ·m̄c²/k trifft auf Faktor 0,8–1,3 (f aus der Polytrope, keine Anpassung).")
P("   Entartete Trichter (Brauner Zwerg, Jupiter, Weißer Zwerg) liegen um Faktor 10–10³ darunter: Dort hält die Zählregel (Pauli) die Tiefe, nicht die Temperatur.")
P()

# Tiefe gegen Zählregel: kT_c/(μ m_u c²) = f·δ − (m_e/(μ_e m_u))·Π(x),  x = p_F/(m_e c) = λ̄_e·(3π²n_e)^{1/3}
def Pe0(rho, mue):
    """Druck des kalten Elektronengases (Chandrasekhar, exakt)."""
    ne = rho / (mue * m_u)
    x = lam_e * (3 * math.pi**2 * ne) ** (1 / 3)
    fx = x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x)
    return m_e * c**2 / (8 * math.pi**2 * lam_e**3) * fx / 3 * 1.0, x

# Normierung prüfen: nichtrelativistisch P = (1/5)·n_e·p_F²/m_e
_rho = 1e5; _P, _x = Pe0(_rho, 2.0)
_Pnr = (1 / 5) * (_rho / (2 * m_u)) * (_x * m_e * c) ** 2 / m_e
P(f"   Kontrolle Elektronendruck bei ρ = 10⁵ kg/m³: exakt {_P:.4e} Pa, nichtrelativistisch {_Pnr:.4e} Pa")

def kT_c(M, R, mu, mue, n):
    t = LE[n]
    rho_c = t["D"] * 3 * M / (4 * math.pi * R**3)
    Pd, x = Pe0(rho_c, mue)
    return mu * m_u * (t["theta"] * G * M / R - Pd / rho_c), x

def Tmax(M, mu, mue, n):
    r = minimize_scalar(lambda lr: -kT_c(M, 10**lr, mu, mue, n)[0], bounds=(3, 11), method="bounded", options=dict(xatol=1e-7))
    R = 10**r.x
    kt, x = kT_c(M, R, mu, mue, n)
    return kt / kB, R, x

P("   Tiefe gegen Zählregel:  kT_c/(μ m_u c²) = f·δ − (m_e/(μ_e m_u))·Π(x),  Π = P_e/(n_e m_e c²),  x = λ̄_e·(3π²n_e)^{1/3}")
P("   (Π = x²/5 nichtrelativistisch, x/4 hochrelativistisch). Beim Zusammenziehen wächst δ ∝ 1/R, Π ∝ 1/R²: T hat ein Maximum.")
MCH_exact = 5.836 / 2.0**2
def M_grenze(n, mue=2.0):
    """Hochrelativistisch: K_ur·ρ_c^{1/3} = θ_n·GM/R mit ρ_c = D_n·3M/(4πR³) → R fällt heraus, eine Masse bleibt."""
    Kur = (3 * math.pi**2) ** (1 / 3) * hbar * c / 4 / (mue * m_u) ** (4 / 3)
    return (Kur * (3 * LE[n]["D"] / (4 * math.pi)) ** (1 / 3) / (LE[n]["theta"] * G)) ** 1.5 / M_sun
P(f"   Kontrolle: Grenzmasse, bei der Tiefe und hochrelativistische Zählregel gleich mit 1/R laufen: n = 3 → {M_grenze(3.0):.3f} M☉"
  f" (Chandrasekhar 5,836/μ_e² = {MCH_exact:.3f} M☉); n = 1,5 → {M_grenze(1.5):.3f} M☉ (Zentralrelation n = 1,5 nur für kleine Massen)")
P()

# ======================================================================= 2 Zündschwellen als Tiefen
P("2) ZÜNDSCHWELLEN ALS TIEFEN")
L = open(os.path.join(DAT, "mass_1.mas20.txt"), encoding="latin-1").read().split("\n")
start = next(i for i, l in enumerate(L) if "micro-u" in l) + 1
MEH, MEn = 7288.971064, 8071.31806
nuc = {}
for l in L[start:]:
    if len(l) < 110:
        continue
    try:
        N_ = int(l[4:9]); Z_ = int(l[9:14]); A_ = int(l[14:19])
    except ValueError:
        continue
    mes = l[28:42].strip(); bas = l[54:67].strip()
    est = ("#" in mes) or ("#" in bas)
    ME = float(mes.replace("#", "."))
    nuc[(Z_, N_)] = dict(A=A_, ME=ME, B=Z_ * MEH + N_ * MEn - ME, meas=not est, el=l[20:23].strip())

def BA(Z, A):
    return nuc[(Z, A - Z)]["B"] / A / 1e3        # MeV
def ME_(Z, A):
    return nuc[(Z, A - Z)]["ME"] / 1e3             # MeV (Atom)

Q_H = 4 * ME_(1, 1) - ME_(2, 4)
Q_3a = 3 * ME_(2, 4) - ME_(6, 12)
Q_C = 2 * ME_(6, 12) - ME_(10, 20) - ME_(2, 4)
Q_Ne = 2 * ME_(10, 20) - ME_(8, 16) - ME_(12, 24)
Q_O = 2 * ME_(8, 16) - ME_(14, 28) - ME_(2, 4)
Q_Si = 2 * ME_(14, 28) - ME_(28, 56)
S_a_Ne = ME_(8, 16) + ME_(2, 4) - ME_(10, 20)
S_a_Si = ME_(12, 24) + ME_(2, 4) - ME_(14, 28)
Q_Fe_dis = 13 * ME_(2, 4) + 4 * ME_(0, 1) - ME_(26, 56) if (0, 1) in nuc else float('nan')
P(f"   AME2020: Q(4¹H→⁴He) = {Q_H:.3f} MeV, Q(3α) = {Q_3a:.3f}, Q(C+C→Ne+α) = {Q_C:.3f}, Q(2Ne→O+Mg) = {Q_Ne:.3f},"
  f" Q(O+O→Si+α) = {Q_O:.3f}, Q(2Si→⁵⁶Ni) = {Q_Si:.3f} MeV")
P(f"           S_α(²⁰Ne) = {S_a_Ne:.3f} MeV, S_α(²⁸Si) = {S_a_Si:.3f} MeV, Q(⁵⁶Fe → 13α + 4n) = {Q_Fe_dis:.2f} MeV")

# Stufen: Name, T_zünd [G], μ der Asche/Brennstoff (voll ionisiert), Reaktion (Z1,A1,Z2,A2) oder photo-S, Q je Nukleon
def mu_ion(Z, A):
    return A / (Z + 1)
stufen = [
    ("pp",  4.0e6, 0.62, (1, 1, 1, 1), None, Q_H / 4, "Mindestwert am unteren Massenende [G]"),
    ("CNO", 1.7e7, 0.62, (7, 14, 1, 1), None, Q_H / 4, "Gleichstand mit pp [G]"),
    ("3α",  1.0e8, mu_ion(2, 4), (2, 4, 4, 8), None, Q_3a / 12, ""),
    ("C",   6.0e8, mu_ion(6, 12), (6, 12, 6, 12), None, Q_C / 24, ""),
    ("Ne",  1.2e9, mu_ion(10, 20), None, S_a_Ne, Q_Ne / 40, "Photozerfall ²⁰Ne(γ,α)"),
    ("O",   1.8e9, mu_ion(8, 16), (8, 16, 8, 16), None, Q_O / 32, ""),
    ("Si",  3.0e9, mu_ion(14, 28), None, S_a_Si, Q_Si / 56, "Photozerfall, Quasi-Gleichgewicht"),
]
f3 = LE[3.0]["theta"]; f15 = LE[1.5]["theta"]
P(f"   Zündtemperaturen [G] (Lehrbuchwerte, ±20 %). δ_zünd = kT/(f·μ m_u c²), f = θ₃ = {f3:.3f} (Kern ~ n = 3).")
P(f"   {'Stufe':5s} {'T_zünd':>8s} {'kT[keV]':>8s} {'μ':>5s} {'δ_zünd':>9s} {'Schritt':>7s} {'E_G[MeV]':>10s} {'E_G/kT':>9s} {'3E₀/kT':>7s} {'E_C[MeV]':>8s} {'kT/E_C':>8s} {'Q/A/m_uc²':>10s} {'Q/A/(m_uc²δ)':>12s}")
dz_prev = None
rows2 = []
for nm, T, mu, rx, S, QA, bem in stufen:
    kT = kB * T
    dz = kT / (f3 * mu * muc2)
    if rx:
        Z1, A1, Z2, A2 = rx
        mr = A1 * A2 / (A1 + A2) * muc2
        EG = 2 * mr * (math.pi * alpha * Z1 * Z2) ** 2
        tau = (27 * EG / (4 * kT)) ** (1 / 3)
        EC = Z1 * Z2 * alpha * hbar * c / (1.2e-15 * (A1 ** (1 / 3) + A2 ** (1 / 3)))   # [A] r₀ = 1,2 fm
        eg_s = f"{EG/MeV:10.3g} {EG/kT:9.3g} {tau:7.1f} {EC/MeV:8.3f} {kT/EC:8.2e}"
    else:
        EG = None; EC = None
        eg_s = f"{'—':>10s} {'S_α/kT=':>9s}{S/(kT/MeV):6.1f}  {'—':>8s} {'—':>8s}"
    step = f"{dz/dz_prev:7.2f}" if dz_prev else f"{'':7s}"
    dz_prev = dz
    rows2.append((nm, T, dz, EG, EC, QA))
    P(f"   {nm:5s} {T:8.1e} {kT/keV:8.2f} {mu:5.2f} {dz:9.2e} {step} {eg_s} {QA/931.494:10.2e} {QA/931.494/dz:12.3g}   {bem}")
# Gesetz T_zünd gegen E_G für die Coulomb-Stufen
cz = [(r[1], r[3]) for r in rows2 if r[3]]
lt = np.log([a for a, b in cz]); le = np.log([b for a, b in cz])
sl, ic, rr, pp_, se = stats.linregress(le, lt)
P(f"   Coulomb-Stufen: ln T_zünd = {sl:.3f}·ln E_G + c (r = {rr:.3f});  T ∝ E_G hieße Steigung 1, feste Gamow-Zahl 3E₀/kT hieße 1.")
P(f"   Gemessen steigt 3E₀/kT von {min((27*r[3]/(4*kB*r[1]))**(1/3) for r in rows2 if r[3]):.0f} auf {max((27*r[3]/(4*kB*r[1]))**(1/3) for r in rows2 if r[3]):.0f}: Die Stufen sind nicht E_G-proportional,")
P("   weil der Querschnittsfaktor S (pp schwach, sonst stark) und die Dichte im Logarithmus stehen. kT/E_C = 10⁻³ … 10⁻² (Tunneln).")
P("   Ne und Si zünden nicht über eine Barriere, sondern über den Photozerfall bei S_α/kT ≈ 30–45.")
P("   Ende der Leiter: Q/A/(m_u c²) fällt von 7·10⁻³ (H) auf 2·10⁻⁴ (Si), δ_zünd steigt von 10⁻⁶ auf 10⁻⁴: Ihr Quotient (Brenn- zu Kontraktionsvorrat)")
P("   fällt von 10³ auf ≈ 1. Bei Si trägt die Kernbindung den Trichter nicht mehr länger, als seine eigene Tiefe ihn trägt; ⁵⁶Ni/⁵⁶Fe ist Q = 0.")
P()

# Längenform der Zündung
P("   Längenform: E_C/(m̄c²) = Z₁Z₂·r_m/R_Kern mit r_m = αħ/(m̄c) (klassischer Radius der Masse m̄), δ = 3Δs/R:")
P("   δ_zünd/(E_C/m_u c²) ist kein fester Bruch (kT/E_C wandert um Faktor 10), also keine reine Längengleichung Kernradius : Faltung.")
r_p_cl = alpha * hbar / (m_p * c)
P(f"   r_p,klass = {r_p_cl:.3e} m;  pp: E_C/(m_p c²) = {r_p_cl/(1.2e-15*2):.2e} (R_Kern = 2,4 fm [A]),  δ_zünd(pp) = {rows2[0][2]:.2e}")
P()

# Mindestmasse
P("   Mindestmasse des Wasserstoffbrennens: T_max(M) aus Tiefe gegen Zählregel (n = 1,5, μ = 0,62, μ_e = 1,143 [A] X = 0,75)")
mu_bd, mue_bd = 0.62, 2 / (1 + 0.75)
for Tign in (2e6, 3e6, 4e6):
    Mm = brentq(lambda m: Tmax(m * M_sun, mu_bd, mue_bd, 1.5)[0] - Tign, 0.005, 0.5)
    Tm, Rm, xm = Tmax(Mm * M_sun, mu_bd, mue_bd, 1.5)
    P(f"   T_zünd = {Tign:.0e} K → M_min = {Mm:.4f} M☉ ({Mm/(sternmasse()/M_sun):.4f} M_*),  R bei T_max = {Rm/R_sun:.3f} R☉,"
      f" δ = {3*ds(Mm*M_sun)/Rm:.2e}, x = λ̄_e k_F = {xm:.3f}")
# geschlossene Form (nichtrelativistisch): kT_max = (θ G M m̄)²/(4B), B = R²·(P/ρ)_deg·m̄
th = LE[1.5]["theta"]; D15 = LE[1.5]["D"]
Knr = (3 * math.pi**2) ** (2 / 3) * hbar**2 / (5 * m_e)
def Mmin_formel(T):
    # (P/ρ)_deg = Knr·(ρ/(μ_e m_u))^{2/3}/(μ_e m_u), ρ_c = D·3M/(4πR³) → Bterm/R², Bterm = Knr·(3DM/(4π μ_e m_u))^{2/3}/(μ_e m_u)
    # kT = μ m_u (θGM/R − Bterm/R²) → max = μ m_u (θGM)²/(4 Bterm)
    g = lambda M: mu_bd * m_u * (th * G * M) ** 2 / (4 * Knr * (3 * D15 * M / (4 * math.pi * mue_bd * m_u)) ** (2 / 3) / (mue_bd * m_u)) - kB * T
    return brentq(g, 1e27, 1e31) / M_sun
P(f"   geschlossene Form (nur nichtrelativistisch, anderer Weg): M_min(3·10⁶ K) = {Mmin_formel(3e6):.4f} M☉;  M_min ∝ M_*·(kT_zünd/m_ec²)^(3/4)")
P(f"   [G] Modelle: 0,075 M☉ bei T_c ≈ 3·10⁶ K. Die Ein-Zonen-Rechnung liegt um Faktor 1,4 darüber (Größenordnung getroffen, nicht der Wert).")
ratio_star = (kB * 3e6 / (m_e * c**2)) ** 0.75
P(f"   Schnittform: M_min/M_* = C·(kT_zünd/m_ec²)^(3/4), (kT/m_ec²)^(3/4) = {ratio_star:.2e}, C = {Mmin_formel(3e6)*M_sun/sternmasse()/ratio_star:.1f}:"
  " eine Potenz der Temperaturskala, kein geometrisches Mittel zweier Leiterlängen.")
P()

# Zündung weiterer Stufen durch Kerne gegebener Masse (T_max, n = 3)
P("   Welche Stufe erreicht ein Kern der Masse M_K (T_max aus Tiefe gegen Zählregel)?")
for nm, mu, mue, T in (("He-Kern → 3α", 4 / 3, 2.0, 1.0e8), ("C/O-Kern → C", 1.75, 2.0, 6.0e8)):
    for n_ in (1.5, 3.0):
        Mg = 0.999 * M_grenze(n_)
        try:
            Mk = brentq(lambda m: Tmax(m * M_sun, mu, mue, n_)[0] - T, 0.01, Mg)
            Tk, Rk, xk = Tmax(Mk * M_sun, mu, mue, n_)
            P(f"   {nm:14s} n = {n_}: T = {T:.1e} K braucht M_K ≥ {Mk:.3f} M☉  (δ = {3*ds(Mk*M_sun)/Rk:.2e}, x = {xk:.2f})")
        except ValueError:
            P(f"   {nm:14s} n = {n_}: nicht erreicht; T_max bei 0,999·M_Grenze = {Tmax(Mg*M_sun, mu, mue, n_)[0]:.2e} K")
P("   [G] Vergleich: kleinster He-brennender Stern ≈ 0,3 M☉, He-Blitz bei Kern 0,47 M☉, C-Zündung ab Kern ≈ 1,05–1,1 M☉.")
P("   T_max hängt nur an M/M_Grenze (beide Terme laufen gleich mit der Strukturzahl): n = 1,5 und n = 3 geben dieselbe Kurve.")
# Decke der Temperatur entarteter Trichter: bei M = M_Grenze ist θGM/R = (m_ec²/μ_e m_u)·x/4 genau, also kT = (μ/μ_e)·m_ec²·[x/4 − Π(x)]
def Pi_(x):
    return (x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x)) / (8 * x**3)
gx = minimize_scalar(lambda x: -(x / 4 - Pi_(x)), bounds=(0.05, 20), method="bounded")
gmax = -gx.fun
P(f"   Decke: Bei M → M_Grenze gilt kT_max = (μ/μ_e)·m_ec²·max_x[x/4 − Π(x)] = (μ/μ_e)·{gmax:.4f}·m_ec² (bei x = {gx.x:.2f});")
P(f"   C/O (μ/μ_e = 0,875): kT_max = {0.875*gmax*m_e*c**2/keV:.1f} keV = {0.875*gmax*m_e*c**2/kB:.2e} K. Kohlenstoff zündet bei {kB*6e8/keV:.1f} keV [G],")
P("   also an der Decke: Ein Trichter, dessen Tiefe die Elektronen-Zählregel hält, erreicht höchstens ≈ 0,1·m_ec²; Ne, O, Si (100–260 keV) nur")
P("   in Kernen über der Grenzmasse, die nie entarten (massive Sterne). Die Elektronen-Compton-Länge setzt die Grenze der Leiter.")
P()

# ======================================================================= 3 Gedächtnis
P("3) TEMPERATUR ALS GEDÄCHTNIS?")
P("   Kelvin-Helmholtz: t_KH = GM²/(RL) = δ·Mc²/L  (Vorrat der Tiefe durch Abstrahlung)")
P("   Kernbrennen:      t_nuc = (Q/A/m_uc²)·X_K·Mc²/L,   t_nuc/t_KH = (Q/A/m_uc²)·X_K/δ")
P("   Entartete Kühlung (Mestel) [G]: t_cool ≈ 8,8·10⁶ a·(A/12)⁻¹·(M/M☉)^(5/7)·(μ/2)^(−2/7)·(L/L☉)^(−5/7)")
Lsun = L_sun
zeiten = [("Sonne", M_sun, R_sun, 1.0 * Lsun, 0.10, Q_H / 4, 4.6),
          ("15 M☉", 15 * M_sun, 4.7 * R_sun, 2.0e4 * Lsun, 0.10, Q_H / 4, 0.01),
          ("0,1 M☉", 0.1 * M_sun, 0.12 * R_sun, 1.0e-3 * Lsun, 1.0, Q_H / 4, 4.6)]
for nm, M, R, Lx, XK, QA, age in zeiten:
    d = 3 * ds(M) / R
    tkh = d * M * c**2 / Lx / Gyr
    tn = QA / 931.494 * XK * M * c**2 / Lx / Gyr
    P(f"   {nm:7s} δ = {d:.2e}: t_KH = {tkh:.3g} Gyr, t_nuc = {tn:.3g} Gyr [L, X_K: G/A], t_nuc/t_KH = {tn/tkh:.0f}")
# Weißer Zwerg: Wärmevorrat gegen Tiefe
Mwd, Rwd, Twd, Awd = 0.6 * M_sun, 0.0125 * R_sun, 1.0e7, 14.0
dwd = 3 * ds(Mwd) / Rwd
Eth_rel = 1.5 * kB * Twd / (Awd * muc2)
for Lx in (1e-2, 1e-3, 1e-4):
    tc = 8.8e6 * (Awd / 12) ** -1 * 0.6 ** (5 / 7) * (1.0) * Lx ** (-5 / 7) / 1e9
    P(f"   Weißer Zwerg 0,6 M☉ bei L = {Lx:.0e} L☉: t_cool ≈ {tc:.3g} Gyr [G-Formel]")
P(f"   Weißer Zwerg: δ = {dwd:.2e}, Wärme je Masse (3/2)kT/(A m_u c²) = {Eth_rel:.2e} bei 10⁷ K: Die Wärme ist {dwd/Eth_rel:.0f}-fach kleiner als die Tiefe.")
P("   Kühlt er, bleibt δ (Zählregel hält R), T fällt. Die Tiefe ist bindend (Struktur), die Temperatur altert (Zustand).")
P("   Im nicht entarteten Trichter umgekehrt: Abstrahlung vertieft (R sinkt), T steigt (negative Wärmekapazität, dE = −dE_th).")
# Grenze nicht entartet / entartet: ψ = kT/E_F entlang der Hauptreihe
P("   Schwelle der Entartung: Hauptreihe mit ρ_c = D_n·ρ̄, kT = fδm̄c², E_F = m_ec²(√(1+x²) − 1): ψ = kT/E_F")
for M in (1.0, 0.3, 0.1, 0.075, 0.05):
    R = {1.0: 1.0, 0.3: 0.29, 0.1: 0.12, 0.075: 0.10, 0.05: 0.10}[M] * R_sun
    kt, x = kT_c(M * M_sun, R, 0.62, mue_bd, 1.5 if M < 0.5 else 3.0)
    EF = m_e * c**2 * (math.sqrt(1 + x * x) - 1)
    d = 3 * ds(M * M_sun) / R
    P(f"   M = {M:5.3f} M☉ (R [G]) δ = {d:.2e}: kT_c = {kt/keV:.3f} keV, E_F = {EF/keV:.3f} keV, ψ = {kt/EF:.2f}")
P("   Die Entartung beginnt dort, wo die Zählregel so viel trägt wie die Tiefe (ψ ≈ 1); ab da hält der Trichter seine Tiefe ohne Wärme.")
P("   Was nicht vergessen wird, ist nicht T, sondern die Asche: Die Zusammensetzung zeichnet die größte je erreichte Tiefe δ_max auf")
P("   (Fusion unter dem Gleichgewicht ist nicht umkehrbar). He-Zwerg, C/O-Zwerg, O/Ne-Zwerg = erreichte Stufe; das Gegenstück zu θ̃ (Minimum der Geschichte).")
P()

# ======================================================================= 4 Elemente
P("4) ELEMENTE ALS SPUR DER TIEFE")
# Häufigkeiten Sonnensystem log ε, H = 12  [G] (Asplund 2009 Photosphäre bzw. Meteorite)
logeps = {1: 12.00, 2: 10.93, 3: 3.26, 4: 1.30, 5: 2.79, 6: 8.43, 7: 7.83, 8: 8.69, 9: 4.56, 10: 7.93, 11: 6.24, 12: 7.60,
          13: 6.45, 14: 7.51, 15: 5.41, 16: 7.12, 17: 5.50, 18: 6.40, 19: 5.03, 20: 6.34, 21: 3.15, 22: 4.95, 23: 3.93,
          24: 5.64, 25: 5.43, 26: 7.50, 27: 4.99, 28: 6.22, 29: 4.19, 30: 4.56, 31: 3.04, 32: 3.65, 33: 2.30, 34: 3.34,
          35: 2.54, 36: 3.25, 37: 2.52, 38: 2.87, 39: 2.21, 40: 2.58, 41: 1.46, 42: 1.88, 44: 1.75, 45: 0.91, 46: 1.57,
          47: 0.94, 48: 1.71, 49: 0.80, 50: 2.04, 51: 1.01, 52: 2.18, 53: 1.55, 54: 2.24, 55: 1.08, 56: 2.18, 57: 1.10,
          58: 1.58, 59: 0.72, 60: 1.42, 62: 0.96, 63: 0.52, 64: 1.07, 65: 0.30, 66: 1.10, 67: 0.48, 68: 0.92, 69: 0.10,
          70: 0.84, 71: 0.10, 72: 0.85, 73: -0.12, 74: 0.85, 75: 0.26, 76: 1.40, 77: 1.38, 78: 1.62, 79: 0.92, 80: 1.17,
          81: 0.90, 82: 1.75, 83: 0.65, 90: 0.02, 92: -0.54}
# B/A je Element: größtes B/A der gemessenen Kerne dieses Z (≈ stabiles Tal)
baZ = {}
for (Z, N), v in nuc.items():
    if v["meas"] and v["A"] > 1:
        b = v["B"] / v["A"] / 1e3
        if Z not in baZ or b > baZ[Z][0]:
            baZ[Z] = (b, v["A"])
baZ[1] = (0.0, 1)
Zs = sorted(logeps)
P("   Häufigkeiten log ε [G] (Asplund 2009), B/A = größte gemessene Bindung je Nukleon des Elements (AME2020)")
def spear(Zsel, diff=False):
    x = np.array([baZ[z][0] for z in Zsel]); y = np.array([logeps[z] for z in Zsel])
    if diff:
        x, y = np.diff(x), np.diff(y)
    r, p = stats.spearmanr(x, y)
    return r, p, len(x)
def partial(Zsel):
    """Spearman log ε gegen B/A bei festgehaltenem Z (Ränge, Residuen)."""
    rz = stats.rankdata(Zsel); rx = stats.rankdata([baZ[z][0] for z in Zsel]); ry = stats.rankdata([logeps[z] for z in Zsel])
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return stats.pearsonr(ex, ey)
for lab, sel in (("Z 6–28 (Brennstufen)", [z for z in Zs if 6 <= z <= 28]), ("Z 29–92 (jenseits Fe)", [z for z in Zs if z >= 29])):
    r, p, nn = spear(sel)
    rp, pp2 = partial(sel)
    P(f"   {lab:24s} Spearman log ε gegen B/A: r_S = {r:+.2f} (p = {p:.2g}, n = {nn});  bei festem Z (Gang abgezogen): {rp:+.2f} (p = {pp2:.2g})")
P(f"   Li, Be, B: log ε = {logeps[3]:.2f}, {logeps[4]:.2f}, {logeps[5]:.2f} gegen He {logeps[2]:.2f}, C {logeps[6]:.2f}; B/A (⁷Li, ⁹Be, ¹¹B) = {BA(3,7):.2f}, {BA(4,9):.2f}, {BA(5,11):.2f} MeV"
  f" unter ⁴He {BA(2,4):.2f} und ¹²C {BA(6,12):.2f}: Die schwach gebundenen Kerne werden in jeder Trichtertiefe über 10⁶ K verbrannt, nicht aufgebaut.")
sel = [z for z in Zs if 6 <= z <= 28 and z + 1 in logeps]
cons = [z for z in range(6, 29)]
r, p, nn = spear(cons, diff=True)
P(f"   Z 6–28 Nachbardifferenzen (gerade/ungerade): r_S = {r:+.2f} (p = {p:.2g}, n = {nn})")
ev = np.mean([logeps[z] for z in range(6, 29) if z % 2 == 0]); od = np.mean([logeps[z] for z in range(6, 29) if z % 2 == 1])
bev = np.mean([baZ[z][0] for z in range(6, 29) if z % 2 == 0]); bod = np.mean([baZ[z][0] for z in range(6, 29) if z % 2 == 1])
P(f"   gerade Z 6–28: ⟨log ε⟩ = {ev:.2f}, ⟨B/A⟩ = {bev:.3f} MeV; ungerade: {od:.2f}, {bod:.3f} MeV")
# Eisengipfel
nb = [z for z in range(15, 41) if z in logeps]
top = sorted(nb, key=lambda z: -logeps[z])[:3]
P(f"   Größte Häufigkeit Z 15–40: {', '.join(f'Z={z} ({logeps[z]:.2f})' for z in top)};  größte B/A: Fe/Ni {BA(26,56):.4f} / {BA(28,62):.4f} MeV")
P(f"   Fe liegt {logeps[26]-logeps[25]:.2f} bzw. {logeps[26]-logeps[27]:.2f} dex über Mn und Co: Gipfel am Maximum der Bindung.")
# Spaltbarkeit x am Fe
ks = [k for k, v in nuc.items() if v["meas"] and v["A"] >= 16]
Za = np.array([k[0] for k in ks], float); Na = np.array([k[1] for k in ks], float); Aa = Za + Na
Ba_ = np.array([nuc[k]["B"] for k in ks]) / 1e3
dlt = np.where((Za % 2 == 0) & (Na % 2 == 0), 1, np.where((Za % 2 == 1) & (Na % 2 == 1), -1, 0)) / np.sqrt(Aa)
X = np.column_stack([Aa, -Aa ** (2 / 3), -Za * (Za - 1) / Aa ** (1 / 3), -(Na - Za) ** 2 / Aa, dlt])
coef = np.linalg.lstsq(X, Ba_, rcond=None)[0]
aV, aS, aC, aA, aP = coef
xFe = aC * 26 * 25 / 56 ** (1 / 3) / (2 * aS * 56 ** (2 / 3))
P(f"   Tröpfchen (AME2020, {len(ks)} Kerne): a_S = {aS:.2f}, a_C = {aC:.3f} MeV; x = E_C/2E_S bei ⁵⁶Fe = {xFe:.3f} (¼ am Maximum)")
P()
P("   Tiefen jenseits der Brennleiter (Klein: bei T = 0 ist das chemische Potential je Baryon × Takt fest, E_F je Baryon = m_uc²·Δδ):")
Qb_Mn56 = ME_(25, 56) - ME_(26, 56)
Qn = ME_(0, 1) - ME_(1, 1)
for lab, E, mue in (("Neutronisierung freier Protonen (m_n − m_H)c²", Qn, 1.0),
                    ("Neutronisierung ⁵⁶Fe → ⁵⁶Mn (AME)", Qb_Mn56, 56 / 26),
                    ("größte Kernbindung B/A (⁶²Ni)", BA(28, 62), 1.0)):
    P(f"   {lab:46s} E = {E:7.3f} MeV → Δδ = E/(μ_e m_uc²) = {E/(mue*931.494):.2e}")
P(f"   Photozerfall des Eisens: Q = {Q_Fe_dis:.1f} MeV = {-Q_Fe_dis/56:.2f} MeV je Nukleon; bei T ≈ 10¹⁰ K [G] → δ_zünd = {kB*1e10/(f3*56/27*muc2):.1e}")
dNS = delta(1.4 * M_sun, 12e3)
P(f"   Neutronenstern 1,4 M☉, 12 km [G]: δ = 1 − dτ/dt = {dNS:.3f}, 3Δs/R = {3*ds(1.4*M_sun)/12e3:.3f};  δ/(B/A_max/m_uc²) = {dNS/(BA(28,62)/931.494):.0f}")
P("   s-Prozess: Neutronen aus ¹³C(α,n), ²²Ne(α,n) in der He-Schale (T ≈ 1–3·10⁸ K [G]) → Stufe 3α der Leiter, δ ≈ 10⁻⁵.")
P("   r-Prozess: Neutronendichte > 10²⁰ cm⁻³ [G] (Verschmelzung/Kollaps von Neutronensternen) → δ ≈ 0,1–0,2.")
# Gipfel jenseits Fe in den Häufigkeiten
ev_z = [z for z in Zs if z > 30 and z % 2 == 0]
peaks = [z for i, z in enumerate(ev_z[1:-1], 1) if logeps[z] > logeps[ev_z[i - 1]] and logeps[z] > logeps[ev_z[i + 1]]]
P(f"   Gipfel jenseits Fe (nur gerade Z, damit der Gerade-ungerade-Gang nicht zählt): Z = {', '.join(str(z) for z in peaks)}")
P("   [G] Zuordnung: s-Gipfel Sr/Ba/Pb (N = 50, 82, 126), r-Gipfel Se–Kr, Te–Xe, Os–Pt (N = 50, 82, 126 bei kleinerem Z). Die Gipfel folgen den")
P("   Schalenzahlen der Neutronen (Zählregel), nicht der Tiefe; die Tiefe entscheidet nur die Neutronendichte (s langsam, r schnell).")
P()

# ======================================================================= 6 Gegenproben
P("6) GEGENPROBEN (anderer Weg)")
# (a) Sonne T_c über den Mittelwert statt Polytrope: Virial ⟨kT⟩ = (1/3)|E_grav|/N, E_grav = −(3/(5−n))GM²/R
n = 3.0; mu = 0.62
Eg = 3 / (5 - n) * G * M_sun**2 / R_sun
kTmean = Eg / 3 / (M_sun / (mu * m_u))
P(f"   (a) Virial n = 3: ⟨kT⟩ = |E_g|/(3N) = {kTmean/eV:.1f} eV = {kTmean/kB:.3e} K (Mittel); T_c/⟨T⟩ bei n = 3 [G] ≈ 2 → {2*kTmean/kB:.2e} K")
P(f"       Polytrope: {tab1[0][2]:.3e} K. Von Hand: δ☉ = 3·492,2 m/6,957·10⁸ m = 2,1225·10⁻⁶; ·0,62·931,494 MeV = 1225,8 eV; ·0,8543 = 1047 eV = 1,215·10⁷ K")
# (b) Zündtiefe C von Hand
P(f"   (b) δ_zünd(C) von Hand: kT = 8,617·10⁻⁵ eV/K·6·10⁸ K = 51,70 keV; /(0,8543·(12/7)·931,494 MeV) = 51,70·10³/1,3642·10⁹ = 3,79·10⁻⁵;  Skript {rows2[3][2]:.3e}")
# (c) Mindestmasse: numerisch (exakter Druck) gegen geschlossene Form
Mnum = brentq(lambda m: Tmax(m * M_sun, mu_bd, mue_bd, 1.5)[0] - 3e6, 0.005, 0.5)
P(f"   (c) M_min(3·10⁶ K): numerisch mit exaktem Druck {Mnum:.4f} M☉, geschlossene Form {Mmin_formel(3e6):.4f} M☉")
P()

txt = "\n".join(out)
print(txt)
with open(os.path.join(HERE, "fusion_trichter.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
