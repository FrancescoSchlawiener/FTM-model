# -*- coding: utf-8 -*-
"""
Ebenen als Projektion: das Feld einer Hierarchie, in der inneren Richtung ϖ harmonisch fortgesetzt
((∇² + ∂_ϖ²)u = 0, Symbol e^{−|k|ϖ}, Poisson-Kern P_ϖ(x) = ϖ/(π²(|x|² + ϖ²)²)).
Nur ftm_kern und weltmodell/daten (SPARC, AME2020). Annahmen [A], Gedächtniswerte [G].
 1 Was die Projektion exakt erhält und was sie abbaut (Zentrum / Hintergrund / Wirbel)
 2 Hierarchie Mond – Erde – Sonne – Milchstraße (– Haufen [A]): Tiefe, bei der eine Ebene in die Elternebene übergeht
 3 Rückkehr als Glättung der Kreisbahnmode, Leiter in ϖ, Enden der Leiter
 4 Kompakt und kollektiv: was die Projektion trennt und was nicht (SPARC, AME2020)
 5 Werte für die Prüfung von Hand
"""
import math, os, sys, io, glob
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from ftm_kern import (G, c, hbar, m_p, m_e, M_sun, M_earth, AU, pc, kpc, Mpc, l_P, m_Pl, A0, ds, r0, r_stern,
                      nu, a_regel, schliesslaenge, compton_red, masse_r_stern_gleich_compton)

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ell0 = c**2 / A0
LAM = schliesslaenge(A0)
M_U = c**4 / (G * A0)
M_c = masse_r_stern_gleich_compton(A0)
K_P = (4 / (3 * math.pi)) ** (1 / 3)          # ϖ_J/ℓ für den Poisson-Kern


def H(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)


# ------------------------------------------------------------------ Kern der Projektion
def m_frac(x):
    """Anteil der Masse innerhalb r, Punktmasse in Tiefe ϖ, x = r/ϖ: (2/π)[atan x − x/(1+x²)]."""
    return (2 / math.pi) * (math.atan(x) - x / (1 + x * x))


def dm_frac(x):
    return 4 * x * x / (math.pi * (1 + x * x) ** 2)


def phi1_phi2_poisson(GM, r, w):
    """Φ′ = GM m(r/ϖ)/r², Φ″ = GM[m′/(ϖr²) − 2m/r³] (ϖ = 0: Punkt)."""
    if w <= 0:
        return GM / r**2, -2 * GM / r**3
    x = r / w
    if x < 1e-3:                                   # Reihe m ≈ (4/3π)x³(1 − 9x²/5·…) → Kern isotrop
        k = 4 * GM / (3 * math.pi * w**3)
        return k * r, k
    m = m_frac(x)
    return GM * m / r**2, GM * (dm_frac(x) / (w * r * r) - 2 * m / r**3)


def phi1_phi2_plummer(GM, r, l):
    return GM * r / (r * r + l * l) ** 1.5, GM * (l * l - 2 * r * r) / (r * r + l * l) ** 2.5


def hesse(GM, dvec, w, kern="poisson"):
    """Hesse-Form von Φ einer geglätteten Punktmasse am Ort dvec (Abstand von der Masse)."""
    r = float(np.linalg.norm(dvec)); n = np.asarray(dvec) / r
    p1, p2 = (phi1_phi2_poisson if kern == "poisson" else phi1_phi2_plummer)(GM, r, w)
    return (p1 / r) * (np.eye(3) - np.outer(n, n)) + p2 * np.outer(n, n)


def T_of(Hm):
    E = Hm - np.trace(Hm) / 3 * np.eye(3)
    return math.sqrt(1.5) * math.sqrt(float(np.sum(E * E)))


def eigen_iso(GM, w, kern="poisson"):
    """isotrope Hesse der eigenen, geglätteten Masse im eigenen Zentrum."""
    return 4 * GM / (3 * math.pi * w**3) if kern == "poisson" else GM / w**3


# ------------------------------------------------------------------ Elterngezeit am Sonnenort (wie trichter_invarianten.py)
v0, R0 = 233e3, 8.178 * kpc                          # [G]
KK = 1e3 / kpc
A_O, B_O = 15.3 * KK, -11.9 * KK                     # [G]
rho_b = 0.084 * M_sun / pc**3                        # [G]
g_obs = v0**2 / R0
g_Ne = brentq(lambda x: float(a_regel(x)) - g_obs, 1e-13, 1e-8)
K_e = float((math.log(nu(g_Ne * 1.0001)) - math.log(nu(g_Ne / 1.0001))) / (2 * math.log(1.0001)))
dlnv = -(A_O + B_O) / (A_O - B_O)
Phi_RR = g_Ne / R0 * (2 * dlnv - 1) / (1 + K_e)
Phi_pp = g_Ne / R0
Phi_zz = 4 * math.pi * G * rho_b - Phi_RR - Phi_pp
H_GAL = np.diag([Phi_RR, Phi_pp, Phi_zz])            # harmonischer Teil über ≪ 100 pc: in ϖ unverändert (1a)
T_GAL = T_of(H_GAL)

GM_S, GM_E = G * M_sun, G * M_earth
M_moon = M_earth / 81.3006                           # [G]
GM_M = G * M_moon
D_EM = 3.844e8                                       # [G]

# ================================================================== 1
H("1   Was die Projektion erhält und was sie abbaut")
print("   Poisson-Kern P_ϖ(x) = ϖ/(π²(r² + ϖ²)²), ∫P_ϖ = 1, P_a·P_b = P_{a+b} (Symbol e^{−|k|ϖ})")
# (a) Fluss = Masse
for x in (1, 10, 100, 1e4):
    print(f"   m(r/ϖ = {x:>7g}) = {m_frac(x):.6f}   (Anteil der Masse in r; 1 − m ≈ (4/π)ϖ/r = {4/math.pi/x:.2e})")
print("   ⇒ Zentrum: der Fluss GM = 3c²Δs durch jede Kugel r ≫ ϖ bleibt; die Faltung Δs ist in ϖ invariant.")
# (b) harmonischer Elternanteil: Schalenmittel eines harmonischen Polynoms = Wert (Mittelwertsatz) → P_ϖ*h = h
D6 = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)   # sphärisches 3-Design
def h_tide(p):
    return p[..., 0]**2 - p[..., 1]**2 + 0.3 * p[..., 0] * p[..., 2] + 0.2 * p[..., 1]   # harmonisch: A·x + E^TF
x0 = np.array([0.7, -0.4, 0.25])
for w in (0.1, 1.0, 10.0):
    t = np.linspace(-math.pi / 2, math.pi / 2, 200001)[1:-1]               # r = ϖ·tan t: Gewicht 4πr²P_ϖ dr = (4/π)sin²t dt
    rr = w * np.tan(t); wt = (4 / math.pi) * np.sin(t)**2 * (t[1] - t[0]) * (t > 0)
    sh = np.array([h_tide(x0 + r * D6).mean() for r in rr[::2000]])
    print(f"   ϖ = {w:>5}: Schalenmittel des Gezeitenpolynoms − Wert: max {np.max(np.abs(sh - h_tide(x0))):.1e} (Mittelwertsatz, Kerngewicht "
          f"{wt.sum():.5f});  isotroper Anteil r²: ∫r²P_ϖ bis r = 10⁴ϖ = {np.sum(wt*rr**2*(rr<1e4*w))/w**2:.0f}·ϖ² (divergiert)")
print("   ⇒ Hintergrund: homogenes Feld A und spurfreie Gezeit E^TF der Elternebene (harmonisch) sind in ϖ exakt unverändert;")
print("     nicht invariant ist nur der Spuranteil (Quelldichte am Ort).")
# (c) eigene Struktur baut ab
alpha_max = max(m_frac(x) / x**2 for x in np.linspace(0.3, 3, 20001))
x_max = np.linspace(0.3, 3, 20001)[np.argmax([m_frac(x) / x**2 for x in np.linspace(0.3, 3, 20001)])]
print(f"   Wirbel: eigene Hesse im Zentrum 4GM/(3πϖ³) ∝ ϖ⁻³; größtes eigenes Feld α·GM/ϖ², α = {alpha_max:.4f} bei r = {x_max:.3f}ϖ;")
print(f"     es fällt bei ϖ_a = √α·r* = {math.sqrt(alpha_max):.3f}·r* auf a₀.")
print("   Triade in der Projektion: invariant bleiben Fluss (Zentrum, Δs) und harmonischer Elternteil (Hintergrund);")
print("     abgebaut wird nur die eigene aufgelöste Struktur (Wirbel), Hesse ∝ ϖ⁻³ gegen die feste Elterngezeit.")

# ================================================================== 2
H("2   Hierarchie Mond – Erde – Sonne – Milchstraße: Tiefe des Übergangs in die Elternebene")
print(f"   Elterngezeit am Sonnenort (Newton-Hesse, [G]-Eingaben wie trichter_invarianten): T = {T_GAL:.4e} s⁻²")
print(f"   Kriterium (07, 13): eigene isotrope Hesse = T des Rests. Poisson-Kern: ϖ_J³ = 4GM/(3πT) ⇒ ϖ_J = (4/3π)^(1/3)·ℓ = {K_P:.4f}·ℓ")
ell = lambda GM, T: (GM / T) ** (1 / 3)
lS, lE, lM = ell(GM_S, T_GAL), ell(GM_E, T_GAL), ell(GM_M, T_GAL)
print(f"   ℓ = (GM/T)^(1/3) in der galaktischen Gezeit: Sonne {lS/pc:.4f} pc, Erde {lE/AU:.0f} AE, Mond {lM/AU:.1f} AE (07: 0,93 pc, 2780, 641)")
print(f"   Poisson-Tiefe ϖ_J = {K_P:.4f}·ℓ: Sonne {K_P*lS/pc:.4f} pc, Erde {K_P*lE/AU:.0f} AE, Mond {K_P*lM/AU:.1f} AE")

# Lagen: Sonne im Ursprung; Richtungen zufällig (Ekliptik ist gegen die Galaxie geneigt) → Spanne
def lagen(seed):
    r = np.random.default_rng(seed)
    e1 = r.normal(size=3); e1 /= np.linalg.norm(e1)
    e2 = r.normal(size=3); e2 -= e2.dot(e1) * e1; e2 /= np.linalg.norm(e2)
    xE = AU * e1
    xM = xE + D_EM * (0.6 * e1 + 0.8 * e2)
    return {"Sonne": (GM_S, np.zeros(3)), "Erde": (GM_E, xE), "Mond": (GM_M, xM)}

def T_rest(name, L, tiefen, kern):
    GMi, xi = L[name]
    Hm = H_GAL.copy()
    for k, (GMk, xk) in L.items():
        if k == name:
            continue
        Hm += hesse(GMk, xi - xk, tiefen[k], kern)
    return T_of(Hm)

print("\n   (A) eine gemeinsame Tiefe ϖ für alle Körper (die Gleichung wörtlich): Verhältnis eigene Hesse / T des Rests")
ws = np.geomspace(1e5, 3e17, 1400)
for seed in (0, 1, 2):
    L = lagen(seed)
    txt = []
    for name in ("Mond", "Erde", "Sonne"):
        rat = np.array([eigen_iso(L[name][0], w) / T_rest(name, L, {k: w for k in L}, "poisson") for w in ws])
        sg = np.sign(np.log(rat)); idx = np.where(np.diff(sg) != 0)[0]
        cr = []
        for i in idx:
            f = lambda lw: math.log(eigen_iso(L[name][0], math.exp(lw)) / T_rest(name, L, {k: math.exp(lw) for k in L}, "poisson"))
            cr.append(math.exp(brentq(f, math.log(ws[i]), math.log(ws[i + 1]))))
        txt.append(f"{name}: " + ", ".join(f"{x/AU:.4g} AE" + (" ↓" if sg[i] > 0 else " ↑") for x, i in zip(cr, idx)))
    print(f"      Lage {seed}: " + " | ".join(txt))
print("      (↓ eigene Ebene geht unter, ↑ taucht wieder auf.) Die letzte Grenze der Erde liegt nahe ϖ_J, der Mond geht in der")
print("      Punktgezeit der Erde unter (Hill-artig) und taucht nicht wieder auf: 641 AE entsteht mit einer gemeinsamen Tiefe nicht.")

print("\n   (B) jeder Körper mit eigener Tiefe (Fixpunkt ϖ_i³ = 4GM_i/(3πT_i), T_i aus den übrigen, je mit ihrer Tiefe geglättet)")
for kern in ("poisson", "plummer"):
    for seed in (0, 1, 2):
        L = lagen(seed)
        tf = {k: (K_P if kern == "poisson" else 1.0) * ell(L[k][0], T_GAL) for k in L}
        for _ in range(60):
            Tn = {k: T_rest(k, L, tf, kern) for k in L}
            tf = {k: ((4 / (3 * math.pi)) if kern == "poisson" else 1.0) ** (1 / 3) * (L[k][0] / Tn[k]) ** (1 / 3) for k in L}
        s = ", ".join(f"{k} {tf[k]/AU:.1f} AE (T/T_gal = {Tn[k]/T_GAL:.4f})" for k in ("Mond", "Erde", "Sonne"))
        print(f"      {kern:8s} Lage {seed}: {s}")
print("   ⇒ 2780 AE / 641 AE (07) folgen, wenn jede Ebene ihre eigene Tiefe hat (Plummer: ℓ; Poisson: 0,752·ℓ).")
print("     Die Reihenfolge in ϖ ist ∝ M^(1/3) in derselben Elterngezeit: Mond < Erde < Sonne; alle drei gehen in die")
print("     kollektive Ebene der Galaxie über, nicht ineinander (der geglättete Kern der Sonne ist isotrop).")

# Hierarchie entlang ϖ am Ort des Mondes: wer trägt die spurfreie Hesse (Hintergrund) und wer die isotrope (Wirbel/Quelle)?
print("\n   Entlang ϖ am Ort des Mondes (Lage 0, gemeinsame Tiefe): Anteile an der spurfreien Hesse")
L = lagen(0)
print("        ϖ          Erde      Sonne     Galaxie   | eigene Hesse Mond / T_gesamt")
for w in (1e6, 1e7, 1e8, 1e9, 1e10, 1e11, AU, 10 * AU, 100 * AU, 482 * AU, 2000 * AU, 1e4 * AU, 0.5 * pc):
    xM = L["Mond"][1]
    TE = T_of(hesse(GM_E, xM - L["Erde"][1], w)); TS = T_of(hesse(GM_S, xM, w))
    tot = T_of(hesse(GM_E, xM - L["Erde"][1], w) + hesse(GM_S, xM, w) + H_GAL)
    print(f"   {w:10.3e} m  {TE/tot:8.2e}  {TS/tot:8.2e}  {T_GAL/tot:8.2e}  | {eigen_iso(GM_M, w)/tot:9.2e}")

# Galaxie im Haufen [A]
M_MW = 6e10 * M_sun; M_H = 1e14 * M_sun; d_H = 1 * Mpc          # [A]
T_H = 3 * G * M_H / d_H**3
print(f"\n   Galaxie 6·10¹⁰ M☉ bei 1 Mpc von 10¹⁴ M☉ [A]: T = 3GM_H/d³ = {T_H:.3e} s⁻²; ℓ = {ell(G*M_MW, T_H)/kpc:.1f} kpc, "
      f"ϖ_J = {K_P*ell(G*M_MW, T_H)/kpc:.1f} kpc, r* = {r_stern(M_MW)/kpc:.2f} kpc")

# ================================================================== 3
H("3   Rückkehr als Glättung der Kreisbahnmode; Leiter in ϖ; Enden")
print("   q = e^{−|k|ϖ_R} mit k = 1/r, ϖ_R = r*: der Anteil der Struktur auf der Bahnskala, der die Tiefe r* übersteht (07)")
cases = [("Erdbahn im Sonnentrichter", M_sun, AU), ("Mond um die Erde", M_earth, D_EM),
         ("Sonne in der Milchstraße (Punkt)", M_MW, 8.2 * kpc), ("Galaxie im Haufen [A]", M_H, 1 * Mpc),
         ("Doppelstern 1,5 M☉, 20 kAE", 1.5 * M_sun, 2e4 * AU)]
for (n_, M, r) in cases:
    s = r_stern(M) / r
    print(f"      {n_:34s} r*/r = {s:10.4g},  q = e^(−s) = {math.exp(-s):.3e},  ν = {1/(-math.expm1(-s)):.4f}")
print(f"   Der Doppelstern liegt mit r* = {r_stern(1.5*M_sun)/AU:.0f} AE < 20 kAE < ϖ_J = {K_P*ell(G*1.5*M_sun, T_GAL)/pc:.2f} pc in der Projektion im Tiefbereich "
      f"(ν = {1/(-math.expm1(-r_stern(1.5*M_sun)/(2e4*AU))):.2f} frei);")
print("   Newton bis ℓ (γ = 1,000, 13) verlangt, dass die Regel das über ℓ geglättete Feld liest: das ist gesetzt (E0), nicht projiziert.")

# WKB-Gegenprobe: Rückkehr durch das in ϖ geglättete Monopolfeld statt mit festem Γ
S_inf = quad(lambda u_: math.sqrt(m_frac(1 / u_)) if u_ > 0 else 1.0, 0, 1, limit=200)[0] + \
        quad(lambda u_: math.sqrt(m_frac(1 / u_)), 1, np.inf, limit=400)[0]
print(f"   Gegenprobe: liest die Rückkehr den in ϖ geglätteten Betrag, Γ(ϖ) = |ḡ(r;ϖ)|/c, so ist s_eff = ∫√(Γ/D)dϖ ≤ {S_inf:.4f}")
print(f"     für jedes r ≪ r* (Punktmasse), ν ≥ {1/(-math.expm1(-S_inf)):.4f}: keine Newton-Grenze. Die Koeffizienten der Rückkehr")
print("     dürfen nicht von ϖ abhängen (07): Die Rückkehr trägt Γ des Orts; geglättet wird nur die Bahnmode, nicht der Verlust.")

print("\n   Leiter in ϖ: Schritt r₀ = 6Δs, Rückkehrtiefe r*, Schluss Λ = ℓ₀/2:  r₀ : r* : Λ = (2δ*)² : 2δ* : 1, Stufen m = r*/r₀ = 1/(2δ*)")
for (n_, M) in (("Proton", m_p), ("Erde", M_earth), ("Sonne", M_sun), ("Milchstraße", M_MW), ("Haufen", 1.7e14 * M_sun)):
    rs = r_stern(M); dst = 3 * ds(M) / rs
    print(f"      {n_:12s} r₀/r* = {r0(M)/rs:.4e}, r*/Λ = {rs/LAM:.4e}, 2δ* = {2*dst:.4e}, m = {rs/r0(M):.4e}")
print(f"   Oberes Ende in ϖ: m = 1 (r₀ = r* = Λ) bei M = M_U/4 = {M_U/4:.3e} kg; Zelle 1,5·M_U hat m = {1/(2*math.sqrt(1.5)):.3f} < 1:")
print(f"     keine Halbgerade, die Rückkehr ist geschlossen. In Raumlängen (01, 08): r* = 3Δs = ℓ₀ bei M_U = {M_U:.3e} kg.")
print(f"   Unteres Ende: r* = λ̄ bei M_c = {M_c:.3e} kg verlangt ħ; die Gleichung (∇² + ∂_ϖ²)u = 0 enthält ħ nicht.")
print(f"     Innerhalb der Gleichung ist die kleinste Länge der Schritt r₀; r₀ = ℓ_P bei M = {l_P*c**2/(2*G):.3e} kg (= m_Pl/2), auch nur mit ħ.")
print("   Glättungslänge: nicht auf der Leiter. ℓ/r* = (M×/M)^(1/6), M× = a₀³/(GT²) hängt an der Elterngezeit (Zustand der Umgebung).")

# ================================================================== 4
H("4   Kompakt und kollektiv: was die Projektion trennt")
for (n_, M, T) in (("Sonne (galakt. Gezeit)", M_sun, T_GAL), ("Milchstraße im Haufen [A]", M_MW, T_H),
                   ("Milchstraße bei M31 [A: 10¹¹ M☉, 785 kpc]", M_MW, 3 * G * 1e11 * M_sun / (785 * kpc)**3)):
    print(f"   {n_:42s} ϖ_J/r* = {K_P*ell(G*M, T)/r_stern(M):8.2f}")
print("   ⇒ ϖ_J > r* für Stern und Galaxie: Die Tiefe des Übergangs trennt die Klassen nicht. Die Dichteregel (ρ > ρ_J) würde")
print("     die Galaxie wie die Sonne glätten (70: Milchstraße Newtonsch). Die Trennung trägt r*/R, die Lage der Hülle zur Rückkehrtiefe.")
print("   Überlappen der Glättungswolken: n·ϖ_J³ = (ρ̄/ρ_J)/π² (Poisson); am Sonnenort ρ̄_* ≈ 0,05–0,1 M☉/pc³ [G] gegen ρ_J = "
      f"{3*T_GAL/(4*math.pi*G)*pc**3/M_sun:.3f} M☉/pc³")

# SPARC: größtes baryonisches Feld je Galaxie (Rotmod, Υ_d = 0,5, Υ_b = 0,7 [A])
gmax, rsR = [], []
for f in sorted(glob.glob(os.path.join(HERE, "..", "daten", "sparc", "Rotmod_LTG", "*_rotmod.dat"))):
    a = np.loadtxt(f, comments="#")
    if a.ndim < 2 or len(a) < 3:
        continue
    R = a[:, 0] * kpc; vg, vd, vb = a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3
    v2 = vg * np.abs(vg) + 0.5 * vd**2 + 0.7 * vb**2
    gb = v2 / R
    Mb = max(v2[-1] * R[-1] / G, 1e30)
    gmax.append(np.max(gb) / A0)
    rsR.append(r_stern(Mb) / R[-1])
gmax = np.array(gmax); rsR = np.array(rsR)
q16, q50, q84 = np.percentile(np.log10(gmax), [16, 50, 84])
print(f"   SPARC ({len(gmax)} Kurven): größtes baryonisches Feld g_b,max/a₀: Median {10**q50:.2f}, 16–84 % {10**q16:.2f}–{10**q84:.2f}, "
      f"Spanne {gmax.min():.3f}–{gmax.max():.1f}; unter a₀ am ganzen Ort: {np.mean(gmax<1)*100:.0f} %")
p16, p50, p84 = np.percentile(rsR, [16, 50, 84])
print(f"   SPARC r*/R_letzter Punkt (M_b aus der letzten Kurvenmasse): Median {p50:.3f}, 16–84 % {p16:.3f}–{p84:.3f}")

# AME2020: alle Kerne als kompakte Ebene
ms = []
for l in open(os.path.join(HERE, "..", "daten", "elemente", "mass_1.mas20.txt"), encoding="latin-1"):
    if len(l) < 120 or not l[0] in " 0":
        continue
    try:
        A_ = int(l[16:19]); mu = float(l[106:109].strip() + l[110:123].replace("#", "").strip()) * 1e-6
    except ValueError:
        continue
    ms.append((A_, mu))
ms = np.array(ms)
u_kg = 1.66053906660e-27
Mk = ms[:, 1] * u_kg; Rk = 1.2e-15 * ms[:, 0] ** (1 / 3)        # R = 1,2·A^(1/3) fm [A]
rsk = np.sqrt(G * Mk / A0)
print(f"   AME2020 ({len(ms)} Kerne): r*/R_Kern = {np.min(rsk/Rk):.1f}–{np.max(rsk/Rk):.1f}; r*/λ̄ = (M/M_c)^(3/2) = "
      f"{np.min((Mk/M_c)**1.5):.1f}–{np.max((Mk/M_c)**1.5):.3g}; ϖ_J (galakt.) = {K_P*ell(G*Mk.min(),T_GAL)*1e3:.2f}–{K_P*ell(G*Mk.max(),T_GAL)*1e3:.2f} mm")
print(f"   Elektron: r*/λ̄ = {(m_e/M_c)**1.5:.3e} (unter dem unteren Ende, Triade nur in der inneren Richtung)")

# ================================================================== 5
H("5   Werte für die Prüfung von Hand")
print(f"   ϖ_J(Erde) = (4GM⊕/(3πT))^(1/3) = {K_P*lE:.5e} m = {K_P*lE/AU:.1f} AE;  ℓ(Erde) = {lE:.5e} m = {lE/AU:.1f} AE;  T = {T_GAL:.5e} s⁻²")
import mpmath as mp
S_mp = mp.quad(lambda u_: mp.sqrt((2/mp.pi)*(mp.atan(1/u_) - (1/u_)/(1+1/u_**2))), [0, 1, 10, mp.inf])
print(f"   S_∞ unabhängig (mpmath, Teilung 0|1|10|∞): {float(S_mp):.5f}")
print(f"   α (Poisson, größtes Feld) = {alpha_max:.5f};  (4/3π)^(1/3) = {K_P:.6f};  S_∞ = {S_inf:.5f}")
print(f"   Sonne: 2δ* = 6Δs/r* = {r0(M_sun)/r_stern(M_sun):.5e}, r*/Λ = {r_stern(M_sun)/LAM:.5e}")
