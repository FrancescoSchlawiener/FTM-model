# -*- coding: utf-8 -*-
"""
Radien eines Trichters auf allen Ebenen und ihre Invarianten.
Aufruf: python trichter_invarianten.py > trichter_invarianten.txt
Rechnet nur mit dem Kern (ftm_kern, a₀ = A0 = 1,1·10⁻¹⁰).
[G] = Wert aus dem Gedächtnis, nicht nachgeprüft;  [A] = Annahme.

1  Radien je Ebene: r₀ = 6Δs, 3Δs, R, r*, ℓ, λ̄, ℓ₀; Verhältnisse
2  Identitäten und Skalengesetze; Leiter in m_Pl und ε_L
3  Klassengrenze r*/R = 1 ⇔ Σ̄ = a₀/(πG); Systemklassen
4  Sonne im Feld der Milchstraße: Regel im Gesamtfeld, 100 AE – 5 pc; welche Definition gibt ℓ = 0,93 pc
5  Zusatzpotential Φ_X und Ablenkung α(b): andere Wege (Quadratur in r bzw. längs der Sichtlinie)
6  Frequenzen der Radien; Spiegel M ↔ m_Pl²/M
"""
import math, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from ftm_kern import *

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
np.seterr(all="ignore")


def H(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)


a0 = A0
l0 = c**2 / a0                       # ℓ₀
eps_L = l_P / l0                     # ε_L = a₀/a_P
t_P = l_P / c
M_U = c**4 / (G * a0)
M_c = masse_r_stern_gleich_compton(a0)
M_min = eps_L * m_Pl
MSUN_PC3 = M_sun / pc**3

# ---------------------------------------------------------------- Elterngezeit am Sonnenort (Newtonsche Hesse-Form)
v0, R0 = 233e3, 8.178 * kpc                          # [G]
KK = 1e3 / kpc                                       # km/s/kpc → s⁻¹
A_O, B_O = 15.3 * KK, -11.9 * KK                     # Oort-Konstanten [G]
rho_b = 0.084 * MSUN_PC3                             # baryonische Dichte in der Scheibe am Sonnenort [G]
g_obs = v0**2 / R0
g_Ne = brentq(lambda x: a_regel(x, a0) - g_obs, 1e-13, 1e-8)   # Newtonsches Feld aus der Regel
K_e = float((math.log(nu(g_Ne * 1.0001, a0)) - math.log(nu(g_Ne / 1.0001, a0))) / (2 * math.log(1.0001)))
dlnv = -(A_O + B_O) / (A_O - B_O)
dlngN = (2 * dlnv - 1) / (1 + K_e)
Phi_RR = g_Ne / R0 * dlngN
Phi_pp = g_Ne / R0
Phi_zz = 4 * math.pi * G * rho_b - Phi_RR - Phi_pp
Hess = np.diag([Phi_RR, Phi_pp, Phi_zz])
ETF = Hess - np.trace(Hess) / 3 * np.eye(3)
normF = math.sqrt(float(np.sum(ETF**2)))
T_gal = math.sqrt(1.5) * normF                       # T = √(3/2)·|E^TF|


def ell(M, T):
    return (G * M / T) ** (1 / 3)


# ================================================================ 1
H("1   Radien eines Trichters (a₀ = 1,1·10⁻¹⁰; Längen in m)")
print(f"   ℓ₀ = c²/a₀ = {l0:.4e} m = {l0/Gpc:.2f} Gpc;  ε_L = ℓ_P/ℓ₀ = {eps_L:.4e}")
print(f"   Elterngezeit am Sonnenort (Newton, [G]: v₀ = 233 km/s, R₀ = 8,178 kpc, A = 15,3, B = −11,9 km/s/kpc, ρ_b = 0,084 M☉/pc³):")
print(f"      g_Ne = {g_Ne/a0:.4f} a₀, K = dlnν/dln g = {K_e:.4f}, Hesse diag = ({Phi_RR:.3e}, {Phi_pp:.3e}, {Phi_zz:.3e}) s⁻²")
print(f"      |E^TF| = {normF:.4e} s⁻²,  T = √(3/2)|E^TF| = {T_gal:.4e} s⁻²,  ρ_J = 3T/(4πG) = {3*T_gal/(4*math.pi*G)/MSUN_PC3:.3f} M☉/pc³")
T_sun_at_earth = 3 * GM_sun / AU**3                  # Punktmasse: √(3/2)·√6·GM/D³ = 3GM/D³
M31, D31 = 1e11 * M_sun, 785 * kpc                   # [A] baryonisch, Abstand [G]
T_LG = 3 * G * M31 / D31**3
print(f"   Wirt Erde: Sonne als Punktmasse bei 1 AE, T = 3GM☉/AE³ = {T_sun_at_earth:.3e} s⁻²")
print(f"   Wirt Milchstraße: M31 als Punktmasse (10¹¹ M☉ [A], 785 kpc [G]), T = {T_LG:.3e} s⁻²")
print(f"   Wirt Sgr A*: Zentrum eines kugelnahen Wirts ⇒ spurfreie Gezeit E^TF = 0, ℓ formal unbegrenzt; Haufen: kein Wirt angesetzt")

SYS = [  # Name, M [kg], R [m], Bezeichnung R, T
    ("Proton", m_p, R_p, "R_p", T_gal),
    ("Elektron", m_e, r_e, "r_e [A]", T_gal),
    ("Erde", M_earth, R_earth, "R⊕", T_sun_at_earth),
    ("Sonne", M_sun, R_sun, "R☉", T_gal),
    ("Neutronenstern", 1.4 * M_sun, 12e3, "12 km", T_gal),
    ("Sgr A*", 4.3e6 * M_sun, None, "r₀ (Horizont)", None),
    ("Milchstraße", 6e10 * M_sun, 15 * kpc, "15 kpc", T_LG),
    ("Haufen 1e14", 1e14 * M_sun, 1 * Mpc, "1 Mpc", None),
    ("Haufen 1,7e14", 1.7e14 * M_sun, 1 * Mpc, "1 Mpc", None),
]
print(f"\n   {'System':15s} {'r₀=6Δs':>10s} {'3Δs':>10s} {'R':>10s} {'r*':>10s} {'ℓ':>10s} {'λ̄':>10s}   R")
ROWS = {}
for nm, M, R, Rn, T in SYS:
    d3 = 3 * ds(M)
    if R is None:
        R = 2 * d3
    rs = r_stern(M, a0)
    lam = compton_red(M)
    lg = ell(M, T) if T else float("nan")
    ROWS[nm] = dict(M=M, R=R, d3=d3, rs=rs, lam=lam, ell=lg, T=T)
    print(f"   {nm:15s} {2*d3:10.3e} {d3:10.3e} {R:10.3e} {rs:10.3e} {lg:10.3e} {lam:10.3e}   {Rn}")

print(f"\n   {'System':15s} {'r*/R':>10s} {'r*/r₀':>10s} {'ℓ/r*':>10s} {'R/3Δs':>10s} {'ℓ/R':>10s} {'ρ̄/ρ_J':>10s} {'δ*':>10s} {'r*/λ̄':>10s}")
for nm, d in ROWS.items():
    rsR = d['rs'] / d['R']
    rr0 = d['rs'] / (2 * d['d3'])
    lr = d['ell'] / d['rs']
    Rd = d['R'] / d['d3']
    lR = d['ell'] / d['R']
    rhoJ = (lR**3) if d['T'] else float('nan')
    dstar = d['d3'] / d['rs']
    print(f"   {nm:15s} {rsR:10.3e} {rr0:10.3e} {lr:10.3e} {Rd:10.3e} {lR:10.3e} {rhoJ:10.3e} {dstar:10.3e} {d['rs']/d['lam']:10.3e}")
S = ROWS["Sonne"]
print(f"\n   Sonne: r* = {S['rs']/AU:.0f} AE, ℓ = {S['ell']/pc:.4f} pc = {S['ell']/S['rs']:.2f} r*;  Proton: r* = {ROWS['Proton']['rs']*1e15:.2f} fm,"
      f" ℓ = {ROWS['Proton']['ell']*1e3:.3f} mm;  Erde: ℓ = {ROWS['Erde']['ell']/AU:.4f} AE (= Hill-Radius), r* = {ROWS['Erde']['rs']/AU:.2f} AE")
print(f"   Sgr A*: r* = {ROWS['Sgr A*']['rs']/pc:.1f} pc = {ROWS['Sgr A*']['rs']/ROWS['Sgr A*']['R']:.3e} r₀;  Milchstraße: r* = {ROWS['Milchstraße']['rs']/kpc:.3f} kpc,"
      f" ℓ (M31) = {ROWS['Milchstraße']['ell']/kpc:.0f} kpc;  Haufen 1e14: r* = {ROWS['Haufen 1e14']['rs']/Mpc:.3f} Mpc")

# unabhängige Handprobe (ohne Kern): r* Sonne, Σ-Grenze, ℓ Sonne über Jacobi-Form
rs_hand = math.sqrt(1.32712440018e20 / 1.1e-10)
print(f"   Handprobe: √(GM☉/a₀) = √(1,327e20/1,1e-10) = {rs_hand:.4e} m = {rs_hand/1.495978707e11:.0f} AE")

# ================================================================ 2
H("2   Identitäten und Skalengesetze")
print(f"   {'System':15s} {'r*²/(3Δs·ℓ₀)':>14s} {'3Δs·λ̄/ℓ_P²':>13s} {'r*²λ̄/(ℓ_P²ℓ₀)':>15s} {'(r*/R)²/(g(R)/a₀)':>18s} {'(ℓ/R)³/(ρ̄/ρ_J)':>16s} {'(ℓ/r*)³/(Ω*²/T)':>16s} {'r*/R/(δ_R/δ*)':>14s}")
for nm, d in ROWS.items():
    M, R = d['M'], d['R']
    i1 = d['rs']**2 / (d['d3'] * l0)
    i2 = d['d3'] * d['lam'] / l_P**2
    i3 = d['rs']**2 * d['lam'] / (l_P**2 * l0)
    i4 = (d['rs'] / R)**2 / (G * M / R**2 / a0)
    if d['T']:
        rhobar = M / (4 / 3 * math.pi * R**3); rhoJ = 3 * d['T'] / (4 * math.pi * G)
        i5 = (d['ell'] / R)**3 / (rhobar / rhoJ)
        Om2 = G * M / d['rs']**3
        i6 = (d['ell'] / d['rs'])**3 / (Om2 / d['T'])
    else:
        i5 = i6 = float('nan')
    i7 = (d['rs'] / R) / ((d['d3'] / R) / (d['d3'] / d['rs']))
    print(f"   {nm:15s} {i1:14.10f} {i2:13.10f} {i3:15.10f} {i4:18.10f} {i5:16.10f} {i6:16.10f} {i7:14.10f}")
print("   (alle Spalten = 1: r*² = 3Δs·ℓ₀; ℓ_P² = 3Δs·λ̄; r*²λ̄ = ℓ_P²ℓ₀; (r*/R)² = g(R)/a₀;")
print("    (ℓ/R)³ = ρ̄/ρ_J (Einzeltrichter ⇔ ℓ > R); (ℓ/r*)³ = Ω(r*)²/T; r*/R = δ(R)/δ*, δ* = 3Δs/r* = √(M/M_U))")

print("\n   Leiter: L = ℓ_P·μ^a·ε_L^b, μ = M/m_Pl  (Exponenten exakt; Probe am Proton)")
mu_p = m_p / m_Pl
Lad = [("λ̄", -1, 0, compton_red(m_p)), ("ℓ_P", 0, 0, l_P), ("3Δs", 1, 0, 3 * ds(m_p)),
       ("r*", 0.5, -0.5, r_stern(m_p, a0)), ("ℓ₀", 0, -1, l0), ("r̃ = √(λ̄ℓ₀)", -0.5, -0.5, math.sqrt(compton_red(m_p) * l0))]
for nm, a, b, val in Lad:
    print(f"      {nm:14s} a = {a:+5.2f}, b = {b:+5.2f}:  Formel/Wert = {l_P * mu_p**a * eps_L**b / val:.12f}")
print(f"      ℓ = ℓ_P·μ^(1/3)·(T t_P²)^(−1/3): Probe Sonne {l_P * (M_sun/m_Pl)**(1/3) * (T_gal*t_P**2)**(-1/3) / S['ell']:.12f};"
      f"  T_gal·t_P² = {T_gal*t_P**2:.3e};  T_gal = {T_gal*R0/a0:.2f}·a₀/R₀")
print(f"      R bei fester Dichte: R ∝ μ^(1/3) ⇒ ℓ/R hängt nur an ρ̄/ρ_J, r*/R ∝ M^(1/6)·ρ̄^(1/3)")
print(f"   Schnittmassen: M_min = {M_min:.3e} kg, M_c = {M_c:.3e} kg (ε_L^(1/3)m_Pl = {eps_L**(1/3)*m_Pl:.3e}), m_Pl = {m_Pl:.3e} kg, M_U = {M_U:.3e} kg (m_Pl/ε_L = {m_Pl/eps_L:.3e})")
print("\n   Abhängigkeit der Verhältnisse:")
print("      r*/3Δs = (M_U/M)^(1/2) = 1/δ*;  r*/λ̄ = (M/M_c)^(3/2);  λ̄/3Δs = (m_Pl/M)²;  λ̄/ℓ₀ = M_min/M   (nur M)")
print("      r*/R = √(πGΣ̄/a₀), Σ̄ = M/(πR²);  R/3Δs = 1/δ_R;  ℓ/R = (ρ̄/ρ_J)^(1/3);  ℓ/r* = (Ω*²/T)^(1/3) ∝ M^(−1/6) bei festem T")
for nm in ("Proton", "Elektron", "Sonne", "Neutronenstern"):
    d = ROWS[nm]
    print(f"      Probe ℓ/r* ∝ M^(−1/6): {nm:15s} ℓ/r* = {d['ell']/d['rs']:.4e},  Sonne·(M☉/M)^(1/6) = {S['ell']/S['rs']*(M_sun/d['M'])**(1/6):.4e}")

# ================================================================ 3
H("3   Klassengrenze r*/R = 1")
Sig = sigma_kollektiv(a0)
Sig_hand = 1.1e-10 / (math.pi * 6.6743e-11) * (3.0856775814913673e16)**2 / 1.98840987e30
print(f"   (r*/R)² = GM/(a₀R²) = πGΣ̄/a₀ ⇒ r*/R = 1 ⇔ Σ̄ = a₀/(πG) = {Sig:.2f} M☉/pc²  (Hand: {Sig_hand:.2f})")
print(f"   Bei fester Dichte ρ̄: r*/R = 1 ⇔ M₁ = a₀³/(G³(4πρ̄/3)²):")
for rho, nm in [(1000.0, "1 g/cm³"), (1.41e3, "Sonne mittel"), (0.29 * MSUN_PC3, "ρ_J am Sonnenort"), (0.1 * MSUN_PC3, "0,1 M☉/pc³")]:
    M1 = a0**3 / (G**3 * (4 * math.pi * rho / 3)**2)
    print(f"      ρ̄ = {rho:.3e} kg/m³ ({nm}): M₁ = {M1:.3e} kg = {M1/M_sun:.3e} M☉")
Mx_T = a0**3 / (G * T_gal**2)
print(f"   Schnitt der Linien g(R) = a₀ und ρ̄ = ρ_J (Elterngezeit T): M× = a₀³/(G·T²) = {Mx_T/M_sun:.3e} M☉, R× = r* = ℓ = a₀/T = {a0/T_gal/pc:.1f} pc")
print(f"      ⇒ ℓ/r* = (M×/M)^(1/6): Sonne (M×/M☉)^(1/6) = {(Mx_T/M_sun)**(1/6):.4f} (direkt {S['ell']/S['rs']:.4f})")
print(f"   Systeme (Σ̄ = M/(πR²) innerhalb R; r*/R = √(Σ̄/{Sig:.0f}); ρ̄/ρ_J mit ρ_J am Sonnenort = {3*T_gal/(4*math.pi*G)/MSUN_PC3:.3f} M☉/pc³ nur als Maßstab)")
CL = [("Sonne", 1.0, R_sun / pc, ""),
      ("Kugelsternhaufen kompakt [G]", 1e5, 3.0, "halbe Masse in r_h"),
      ("Kugelsternhaufen dünn (Pal 14) [G]", 1e4, 25.0, "halbe Masse in r_h"),
      ("offener Haufen [G]", 1e3, 2.0, ""),
      ("Molekülwolke (Larson) [G]", 1e5, 20.0, ""),
      ("Scheibe, Freeman Σ₀ [G]", None, None, "Σ₀ = 140 M☉/pc²"),
      ("Milchstraße (Modell)", 6e10, 15e3, ""),
      ("Elliptische 1e11 [G]", 5e10, 5e3, "halbe Masse in R_e"),
      ("Zwerggalaxie Draco [G]", 2.9e5, 221 * 1.3, "halbe Masse in r_h"),
      ("ultradiffus DF2 [G]", 1e8, 2.2e3 * 1.3, "halbe Masse in r_h"),
      ("Haufen 1e14 (1 Mpc)", 1e14, 1e6, ""),
      ("Haufen 1,7e14 (1 Mpc)", 1.7e14, 1e6, "")]
for nm, M, R, note in CL:
    if M is None:
        Sb = 140.0
        print(f"      {nm:36s} Σ̄ = {Sb:9.3g} M☉/pc²  r*/R = {math.sqrt(Sb/Sig):7.3f}  {note}")
        continue
    Sb = M / (math.pi * R**2)
    rb = M / (4 / 3 * math.pi * R**3)
    print(f"      {nm:36s} Σ̄ = {Sb:9.3g} M☉/pc²  r*/R = {math.sqrt(Sb/Sig):7.3f}  ρ̄/ρ_J = {rb/(3*T_gal/(4*math.pi*G)/MSUN_PC3):9.3g}  {note}")

# ================================================================ 4
H("4   Sonne im Feld der Milchstraße: Regel im Gesamtfeld (Vektor), 100 AE – 5 pc")
print(f"   g_e = g_Ne = {g_Ne:.4e} m/s² = {g_Ne/a0:.4f} a₀ (Newton, aus g_obs = v₀²/R₀ über die Regel), ν(g_e) = {float(nu(g_Ne,a0)):.4f}")
print("   a_int = ν(|g☉ + g_e|)(g☉ + g_e) − ν(g_e)g_e, radiale Komponente längs r̂; Verhältnis zu g☉ = GM☉/r²")
print(f"   {'r':>12s} {'g☉/a₀':>10s} {'T·r/g☉':>10s} {'r̂ ∥ g_e':>10s} {'r̂ ∥ −g_e':>10s} {'r̂ ⊥ g_e':>10s} {'Mittel':>10s} {'ν(g☉) frei':>11s}")
e_hat = np.array([1.0, 0.0, 0.0]); ge_v = g_Ne * e_hat
XG, WG = np.polynomial.legendre.leggauss(48)


def a_int_rad(r, n):
    gs = -GM_sun / r**2 * n                            # zur Sonne hin
    gt = gs + ge_v
    at = float(nu(np.linalg.norm(gt), a0)) * gt - float(nu(g_Ne, a0)) * ge_v
    return -float(at @ n) / (GM_sun / r**2)


for r in [100 * AU, 1e3 * AU, 3e3 * AU, 7.34e3 * AU, 2e4 * AU, 5e4 * AU, 0.5 * pc, 0.93 * pc, 1.37 * pc, 2 * pc, 5 * pc]:
    par = a_int_rad(r, np.array([1.0, 0, 0])); anti = a_int_rad(r, np.array([-1.0, 0, 0])); perp = a_int_rad(r, np.array([0, 1.0, 0]))
    mean = 0.5 * sum(w * a_int_rad(r, np.array([x, math.sqrt(1 - x * x), 0])) for x, w in zip(XG, WG))
    gs = GM_sun / r**2
    lab = f"{r/AU:.0f} AE" if r < 0.3 * pc else f"{r/pc:.2f} pc"
    print(f"   {lab:>12s} {gs/a0:10.3e} {T_gal*r/gs:10.3e} {par:10.4f} {anti:10.4f} {perp:10.4f} {mean:10.4f} {float(nu(gs,a0)):11.4f}")
print("   Grenzfall g☉ ≪ g_e: linear a_int ≈ ν_e(1 + K·cos²ψ)g☉ radial: ∥ " + f"{float(nu(g_Ne,a0))*(1+K_e):.4f}, ⊥ {float(nu(g_Ne,a0)):.4f}")

print("\n   Kandidaten für die Grenze:")
r_eq = math.sqrt(GM_sun / g_Ne)
print(f"      |g☉| = g_e (Newton):             r = {r_eq/AU:.0f} AE = {r_eq/pc:.4f} pc = {r_eq/S['rs']:.3f} r*")
g_eobs = g_obs
print(f"      |g☉| = g_obs (verstärkt):        r = {math.sqrt(GM_sun/g_eobs)/pc:.4f} pc")
print(f"      |g☉| = a₀ (r*):                  r = {S['rs']/pc:.4f} pc")
print(f"      Tide T·r = a₀:                   r = a₀/T = {a0/T_gal/pc:.1f} pc")
T_J = 4 * A_O * (A_O - B_O)
cand = [("T = 4A(A−B) (Jacobi, Scheibenebene)", T_J),
        ("T = |E^TF| (Frobenius)", normF),
        ("T = √(3/2)|E^TF| (Modell), Newton-Hesse", T_gal),
        ("T = größter Eigenwert von E^TF", float(np.max(np.abs(np.linalg.eigvalsh(ETF))))),
        ("T = Φ_zz (vertikal)", Phi_zz)]
for nm, T in cand:
    print(f"      GM☉/r² = T·r mit {nm:42s}: r = {ell(M_sun,T)/pc:.4f} pc = {ell(M_sun,T)/S['rs']:.2f} r*  (T = {T:.3e} s⁻²)")
# unabhängige Probe von ℓ: Nullstelle der Gleichgewichtsbedingung statt Kubikwurzel
lr = brentq(lambda r: GM_sun / r**2 - T_gal * r, 0.01 * pc, 10 * pc)
print(f"   Probe (Nullstelle GM/r² − T·r statt Kubikwurzel): ℓ = {lr/pc:.5f} pc;  Umlauf an ℓ 2π/√T = {2*math.pi/math.sqrt(T_gal)/yr/1e6:.1f} Mio. Jahre")

# ================================================================ 5
H("5   Zusatzpotential Φ_X und Ablenkung α(b): andere Wege")
M = 6e10 * M_sun
rs = r_stern(M, a0); vf2 = v_flach(M, a0)**2
print(f"   Φ_X(r) = ∫₀^r (ν−1)·a_N dr' (Quadratur in r) gegen −v_flach²·ln(1 − e^(−r*/r)), Milchstraße als Punktmasse, r* = {rs/kpc:.3f} kpc")
for x in [0.1, 0.3, 1.0, 3.0, 10.0, 100.0]:
    r = x * rs
    f = lambda u: (float(nu(G * M / u**2, a0)) - 1) * G * M / u**2
    num = quad(f, 1e-6 * rs, r, limit=400, points=[min(r, rs)] if r > rs else None)[0]
    ana = -vf2 * math.log(-math.expm1(-rs / r))
    far = vf2 * math.log(r / rs)
    print(f"      r/r* = {x:6.1f}: Quadratur {num/vf2:12.8f} v²,  geschlossen {ana/vf2:12.8f} v²,  v²·ln(r/r*) = {far/vf2:9.5f} v²")
# Ablenkung längs der Sichtlinie: α = (2/c²)∫ a(r)·b/r dz
print("\n   α(b) = (2/c²)∫_{−∞}^{∞} ν(GM/r²)·GM·b/r³ dz  (Quadratur in z), I = α/(12Δs/b)")


def I_z(x, M=M):
    b = x * r_stern(M, a0)
    zs = np.concatenate([[0.0], np.geomspace(1e-6 * b, 1e8 * b, 200001)])
    fz = nu(G * M / (b * b + zs**2), a0) * G * M * b / (b * b + zs**2)**1.5
    val = 2 * np.trapezoid(fz, zs)
    alpha = 2 / c**2 * val
    return alpha / (12 * ds(M) / b)


def I_theta(x):
    return quad(lambda t: math.cos(t) * float(nu(math.cos(t)**2 / x**2 * a0, a0)), 0, math.pi / 2, limit=400)[0]


for x in [0.1, 0.5, 1.0, 2.0, 10.0, 100.0]:
    Iz = I_z(x); It = I_theta(x)
    asym = (math.pi / 2) * x + 0.5            # α ≈ 6πΔs/r* + 6Δs/b ⇒ I ≈ (π/2)·b/r* + 1/2
    print(f"      b/r* = {x:6.1f}: I (z-Quadratur) = {Iz:.6f},  I (θ-Form, quad) = {It:.6f},  Näherung (π/2)x + ½ = {asym:.4f}")
print(f"   α_∞ = 2π(v_flach/c)²: Milchstraße {ablenkung_kollektiv(6e10*M_sun,a0)*206264.806:.3f}″, Haufen 1e14 {ablenkung_kollektiv(1e14*M_sun,a0)*206264.806:.2f}″,"
      f" 1,7e14 {ablenkung_kollektiv(1.7e14*M_sun,a0)*206264.806:.2f}″")

# ================================================================ 6
H("6   Frequenzen der Radien und Spiegel M ↔ m_Pl²/M")
print("   Ω(r*) = √(GM/r*³):  Ω(r*)⁴ = (c/3Δs)·(a₀/c)³  ⇔  Ω(r*) = (c/3Δs)·δ*^(3/2) = (a₀/c)·δ*^(−1/2)")
for nm in ("Proton", "Sonne", "Milchstraße"):
    d = ROWS[nm]; Om = math.sqrt(G * d['M'] / d['rs']**3); dst = d['d3'] / d['rs']
    print(f"      {nm:12s} Ω(r*) = {Om:.4e} s⁻¹;  (c/3Δs)^(1/4)(a₀/c)^(3/4) = {(c/d['d3'])**0.25*(a0/c)**0.75:.4e};"
          f"  (c/3Δs)δ*^(3/2) = {c/d['d3']*dst**1.5:.4e};  π·f*_gw = {math.pi*schwellenfrequenz_gw(d['M'],a0):.4e}")
print("   Drei Raten je Masse (09): ω_C² = Γ_G·c/(3Δs);  hier dazu Ω(r*)⁴ = (c/3Δs)(a₀/c)³ und Ω(ℓ)² = T (Elterngezeit), Ω(R)²/T = (ℓ/R)³")
print("\n   Spiegel M ↔ m_Pl²/M: 3Δs ↔ λ̄, ℓ₀ fest, r* ↦ r̃ = √(λ̄ℓ₀) = √(ħτ_D/M) (τ_D = c/a₀)")
for nm, m in [("Elektron", m_e), ("Proton", m_p), ("Sonne", M_sun)]:
    rt = math.sqrt(compton_red(m) * l0); rt2 = r_stern(m_Pl**2 / m, a0); rt3 = math.sqrt(hbar * tau_D(a0) / m)
    print(f"      {nm:9s} r̃ = {rt:.4e} m;  r*(m_Pl²/M) = {rt2:.4e};  √(ħτ_D/M) = {rt3:.4e};  r*·r̃/(ℓ_Pℓ₀) = {r_stern(m,a0)*rt/(l_P*l0):.10f}")
print(f"   r* = r̃ bei M = m_Pl; r̃ = λ̄ bei M = M_min = {M_min:.3e} kg")
Mx = brentq(lambda lm: math.log(math.sqrt(compton_red(10**lm) * l0) / (3 * ds(10**lm))), -40, 60)
print(f"      r̃ = 3Δs bei M = {10**Mx:.3e} kg = m_Pl²/M_c = {m_Pl**2/M_c:.3e} kg (Spiegel von M_c)")
