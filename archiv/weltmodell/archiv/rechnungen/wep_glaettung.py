# -*- coding: utf-8 -*-
"""
WEP und Glättungsschwelle: örtliche Regel a = ν(|ḡ_N|)·ḡ_N im Gesamtfeld, geglättet über (ρ − ρ_J)₊.

Quelle der Regel ist allein das geglättete Newton-Feld ḡ_N der Baryonen; keine gelesene Dichte als Quelle.
Kern: Plummer mit Skala ℓ, ℓ³ = GM^c/T, T = √(3/2)|Ē^TF| (Elterngezeit), ρ_J = 3T/(4πG).

  0  Milchstraßenmodell [A] und Gezeitendichte ρ_J am Sonnenort und bei 15 kpc
  1  Glättungslängen: Gezeit aus σ̄ (geglättetes Feld, 07 kovariant) gegen Hill-Radius im wirklichen Elternfeld
  2  Lesart R1: Regel am Schwerpunkt (Weltlinie) – was ein Körper von seinen Nachbarn in ḡ_N sieht
  3  Lesart R2: Regel über die geglättete Zone gemittelt – Gezeitenterm D(ℓ) und Eigenterm E(M)
  4  WEP-Differenzen gegen LLR, MICROSCOPE, Galaxienbeschleunigung
  5  Glättungsschwelle als Phase: Wolke gegen warmes Medium bei R₀ und 15 kpc
[G] = Messwert aus dem Gedächtnis, nicht nachgeprüft; [A] = Annahme.
"""
import math
import numpy as np
from ftm_kern import G, pc, AU, M_sun, GM_sun, GM_earth, R_earth, m_p, A0, nu, a_regel

out = []
def P(s=""): out.append(s); print(s)
def H(t): P(); P("=" * 100); P(t); P("=" * 100)

kpc = 1e3 * pc
MSPC = M_sun / pc**3
R0 = 8.178 * kpc                                      # [G]
V0 = 233e3                                            # [G]
T_SONNE = 5.549e-30                                   # kovariante_e0: √(3/2)|∂∂Φ̄_N^TF| baryonisch am Sonnenort
GM_moon = GM_earth * 0.0123000371                     # [G] Massenverhältnis Mond/Erde
D_MOON = 3.844e8

# ============================================================================ 0
H("0   Milchstraßenmodell [A] und Gezeitendichte")
# Miyamoto-Nagai-Summe (M [M☉], a [kpc], b [kpc]); angepasst an g_N(R₀), ρ_b(R₀, 0) und T(R₀) der Sonnenumgebung
MN = [(1.2e10, 0.0, 0.6), (5.3e10, 3.0, 0.45), (1.2e10, 3.3, 1.0), (1.3e10, 7.0, 0.15)]

def gN(x, y, z):
    gx = gy = gz = 0.0
    for M, a, b in MN:
        GMk = G * M * M_sun; a *= kpc; b *= kpc
        sb = np.sqrt(z * z + b * b); D = (x * x + y * y + (a + sb) ** 2) ** 1.5
        gx = gx - GMk * x / D; gy = gy - GMk * y / D; gz = gz - GMk * z * (a + sb) / (sb * D)
    return np.array([gx, gy, gz]) if np.ndim(gx) == 0 else (gx, gy, gz)

def f_rule(g):                                        # ν − 1
    return nu(g) - 1.0

def g_s(x, y, z):
    """Zusatzfeld der örtlichen Regel (ν(|ḡ_N|) − 1)·ḡ_N"""
    gx, gy, gz = gN(x, y, z)
    gm = np.sqrt(gx * gx + gy * gy + gz * gz)
    f = f_rule(gm)
    return f * gx, f * gy, f * gz

def hesse(R):
    h = 1 * pc; Hm = np.zeros((3, 3))
    for j, e in enumerate(np.eye(3)):
        Hm[:, j] = -(np.array(gN(*(np.array([R, 0, 0]) + h * e))) - np.array(gN(*(np.array([R, 0, 0]) - h * e)))) / (2 * h)
    return Hm

def gezeit(R):
    Hm = hesse(R); TF = Hm - np.trace(Hm) / 3 * np.eye(3)
    return math.sqrt(1.5) * math.sqrt(float(np.sum(TF ** 2))), np.trace(Hm) / (4 * math.pi * G)

ORTE = {}
for nm, R in [("R₀ = 8,18 kpc", R0), ("15 kpc", 15 * kpc)]:
    g = float(-gN(R, 0, 0)[0]); T, rho = gezeit(R)
    ORTE[nm] = dict(R=R, g=g, T=T, rho=rho)
    P(f"   {nm:14s}: g_N = {g:.3e} m/s² = {g/A0:.3f} a₀,  v_Regel = {math.sqrt(a_regel(g)*R)/1e3:.1f} km/s,  "
      f"ρ_b(z=0) = {rho/MSPC:.3f} M☉/pc³,  T = {T:.3e} s⁻²,  ρ_J = {3*T/(4*math.pi*G)/MSPC:.3f} M☉/pc³")
P(f"   Soll am Sonnenort [G]: v₀ = 233 km/s → g_obs = {V0**2/R0:.3e}; g_N (Regel umgekehrt, a₀ = 1,1) = 1,474e-10; ρ_b = 0,084; T = 5,549e-30 (kovariante_e0)")
rhoJ0 = 3 * T_SONNE / (4 * math.pi * G)
P(f"   ρ_J am Sonnenort aus T = 5,549e-30: {rhoJ0:.4e} kg/m³ = {rhoJ0/MSPC:.4f} M☉/pc³ = {rhoJ0/m_p/1e6:.2f} m_p/cm³ "
  f"= {rhoJ0/(1.36*m_p)/1e6:.2f} H/cm³ (mit Helium, 1,36 m_p je H)")
skal15 = ORTE["15 kpc"]["T"] / ORTE["R₀ = 8,18 kpc"]["T"]
T15 = T_SONNE * skal15; rhoJ15 = 3 * T15 / (4 * math.pi * G)
P(f"   bei 15 kpc (Modellverhältnis T₁₅/T₀ = {skal15:.3f}): T = {T15:.3e}, ρ_J = {rhoJ15/MSPC:.4f} M☉/pc³ = {rhoJ15/m_p/1e6:.2f} m_p/cm³ = {rhoJ15/(1.36*m_p)/1e6:.2f} H/cm³")

# ============================================================================ 1
H("1   Glättungslängen ℓ = (GM/T)^(1/3)")
P("   (σ̄) Gezeit aus dem geglätteten Feld (07 kovariant): Die Plummer-Wolke der Sonne ist im Kern isotrop, ihr spurfreier")
P("       Hesse-Anteil verschwindet für r ≪ ℓ☉; am Ort jedes Körpers im Sonnensystem bleibt die galaktische Gezeit T.")
P("   (Hill) Gezeit des wirklichen Elternkörpers 3GM_E/d³ (Erde in der Sonne, Mond in der Erde, Laborkörper an der Erdoberfläche).")
KOERPER = [("Sonne", GM_sun, 3 * T_SONNE / 3),
           ("Erde", GM_earth, 3 * GM_sun / AU**3),
           ("Mond", GM_moon, 3 * GM_earth / D_MOON**3),
           ("Jupiter", 1.26686534e17, 3 * GM_sun / (5.2 * AU)**3),
           ("Laborkörper 1 kg", G * 1.0, 3 * GM_earth / R_earth**3),
           ("Stern 0,1 M☉", 0.1 * GM_sun, T_SONNE),
           ("Pt 0,4 kg", G * 0.4, 3 * GM_earth / (R_earth + 7.1e5)**3),     # [G] MICROSCOPE-Massen
           ("Ti 0,3 kg", G * 0.3, 3 * GM_earth / (R_earth + 7.1e5)**3)]
ELL = {}
P(f"   {'Körper':18s} {'ℓ (σ̄) [m]':>12s} {'= AE':>10s} {'GM/ℓ² (σ̄)':>12s} {'ℓ (Hill) [m]':>13s} {'GM/ℓ² (Hill)':>13s}")
for nm, GMk, Thill in KOERPER:
    ls = (GMk / T_SONNE) ** (1 / 3)
    lh = (GMk / Thill) ** (1 / 3) if nm not in ("Sonne", "Stern 0,1 M☉") else ls
    ELL[nm] = (GMk, ls, lh)
    P(f"   {nm:18s} {ls:12.3e} {ls/AU:10.3g} {GMk/ls**2:12.3e} {lh:13.3e} {GMk/lh**2:13.3e}")
P(f"   Identität: Im Plummer-Kern (d ≪ ℓ) ist das geglättete Feld eines Körpers GM·d/ℓ³ = T·d, für jede Masse gleich.")

# ============================================================================ 2
H("2   Lesart R1: örtliche Regel am Schwerpunkt, ḡ_N(x_i) = Σ_{j≠i} geglättete Felder (eigenes Feld am Schwerpunkt = 0)")
gR0 = ORTE["R₀ = 8,18 kpc"]["g"]
P(f"   Galaktisches Feld am Sonnenort g_N = {gR0:.3e}; Zusatz (ν−1)g_N = {f_rule(gR0)*gR0:.3e} m/s²")
def plummer_feld(GMk, l, d): return GMk * d / (d * d + l * l) ** 1.5
FAELLE = [("Mond am Ort der Erde", "Mond", D_MOON), ("Erde am Ort des Mondes", "Erde", D_MOON),
          ("Jupiter am Ort der Sonne", "Jupiter", 5.2 * AU), ("Sonne am Ort der Erde", "Sonne", AU),
          ("Erde am Ort eines Laborkörpers", "Erde", R_earth)]
P(f"   {'Nachbarfeld in ḡ_N':34s} {'σ̄ [m/s²]':>11s} {'ν−1 (σ̄)':>10s} {'Hill [m/s²]':>12s} {'ν−1 (Hill)':>11s}")
for nm, k, d in FAELLE:
    GMk, ls, lh = ELL[k]
    for_s = plummer_feld(GMk, ls, d); for_h = plummer_feld(GMk, lh, d)
    gs_ = math.hypot(gR0, for_s); gh_ = math.hypot(gR0, for_h)
    P(f"   {nm:34s} {for_s:11.3e} {float(f_rule(gs_)):10.4f} {for_h:12.3e} {float(f_rule(gh_)):11.3e}")
P("   σ̄: Jeder Körper im Sonnensystem sieht ν des galaktischen Feldes (Nachbarn ≤ T·30 AE = "
  f"{T_SONNE*30*AU:.1e} m/s²). Der Empfänger-ℓ kommt in R1 nicht vor: Δa = 0 exakt für Körper am selben Ort.")
P("   Hill: Nachbarfelder ≫ a₀ an jedem Schwerpunkt, ν − 1 → 0; das Sonnensystem erhielte keinen galaktischen Zusatz,")
P(f"         die Sonne liefe mit v = {math.sqrt(gR0*R0)/1e3:.0f} statt {math.sqrt(a_regel(gR0)*R0)/1e3:.0f} km/s (Galaxienbeschleunigung gemessen 2,3e-10 [G]). Hill-ℓ trägt nicht.")
# Gegenseitiger Zusatz Erde–Mond aus T·d (kein WEP-Term, gleich je Masse)
fe = float(f_rule(gR0)); dfe = float((f_rule(gR0 * 1.001) - f_rule(gR0 * 0.999)) / (0.002 * gR0))
P(f"   σ̄, gegenseitiger Zusatz Erde–Mond aus T·d = {T_SONNE*D_MOON:.2e} m/s²: ≤ (ν−1+g·dν/dg)·T·d = {(fe+abs(gR0*dfe))*T_SONNE*D_MOON:.1e} m/s² je Körper (je Masse gleich)")

# ============================================================================ 3
H("3   Lesart R2: Regel über die Plummer-Zone gemittelt,  a_i = ∫K_ℓ(y)·[(ν(|ḡ_N|)−1)ḡ_N](x_i + y) d³y")
P("   (a) Gezeitenterm D(ℓ) = ∫w_ℓ(r)[S(r) − g_s(x₀)]dr, S(r) Kugelmittel des Zusatzfelds, w_ℓ = 3ℓ²r²/(r²+ℓ²)^(5/2)")
NMU, NPH = 64, 64
mu, wmu = np.polynomial.legendre.leggauss(NMU)
ph = (np.arange(NPH) + 0.5) * 2 * np.pi / NPH
MU, PH = np.meshgrid(mu, ph, indexing="ij"); W2 = np.repeat(wmu[:, None], NPH, axis=1) / (2 * NPH)
SN = np.sqrt(1 - MU**2)
def kugelmittel(x0, r):
    dx, dy, dz = r * SN * np.cos(PH), r * SN * np.sin(PH), r * MU
    a = g_s(x0 + dx, dy, dz); b = g_s(x0 - dx, -dy, -dz)
    return np.array([np.sum(W2 * 0.5 * (a[k] + b[k])) for k in range(3)])

def D_tabelle(x0):
    g0 = np.array(g_s(x0, 0.0, 0.0))
    rr = np.logspace(math.log10(3 * pc), math.log10(3e3 * kpc), 260)
    dS = np.array([kugelmittel(x0, r)[0] - g0[0] for r in rr])
    C = dS[0] / rr[0]**2
    return g0[0], rr, dS, C

def D_von(ell, tab):
    g0, rr, dS, C = tab
    lo = np.logspace(math.log10(ell) - 6, math.log10(rr[0]), 400)[:-1]
    r = np.concatenate([lo, rr]); d = np.concatenate([C * lo**2, dS])
    w = 3 * ell**2 * r**2 / (r**2 + ell**2) ** 2.5
    return np.trapezoid(w * d, r)

TAB0 = D_tabelle(R0)
g0s, rr, dS, C0 = TAB0
# Probe: kleines r gegen Laplace-Operator (finite Differenzen, h = 5 pc)
h = 5 * pc
lap = sum((g_s(*(np.array([R0, 0, 0]) + h * e))[0] + g_s(*(np.array([R0, 0, 0]) - h * e))[0] - 2 * g_s(R0, 0, 0)[0]) / h**2 for e in np.eye(3))
P(f"   Zusatzfeld am Sonnenort g_s = {g0s:.4e} m/s² (radial);  S(r)−g_s ≈ C r²: C = {C0:.4e} m⁻¹s⁻²,  ∇²g_s/6 = {lap/6:.4e}")
for rp in [10, 100, 300, 1000, 3000]:
    i = np.argmin(abs(rr - rp * pc)); P(f"      r = {rr[i]/pc:7.0f} pc: (S − g_s)/(C r²) = {dS[i]/(C0*rr[i]**2):.3f}")
# Empfindlichkeit auf die Dicke der dünnen Scheibe
MN_std = list(MN)
for bt in (0.30, 0.60):
    MN[1] = (5.3e10, 3.0, bt)
    Tb = D_tabelle(R0)
    P(f"      Probe dünne Scheibe b = {bt:.2f} kpc: g_N = {-gN(R0,0,0)[0]:.3e}, C = {Tb[3]:.3e}, D(ℓ☉) = {D_von(ELL['Sonne'][1], Tb):.3e} m/s²")
MN[:] = MN_std
P(f"      Standard b = 0,45 kpc: D(ℓ☉) = {D_von(ELL['Sonne'][1], TAB0):.3e} m/s²;  Näherung 3Cℓ²(ln(2L/ℓ) − 4/3) mit L = 400 pc: "
  f"{3*C0*ELL['Sonne'][1]**2*(math.log(2*400*pc/ELL['Sonne'][1]) - 4/3):.3e}")
P("   (b) Eigenterm E(M): das eigene geglättete Feld u(y) = −GM·y/(y²+ℓ²)^(3/2) im Argument von ν, gleichförmiges Außenfeld g")
def E_eigen(GMk, ell, g):
    x = np.logspace(-4, 3, 700) ; mu_, wm = np.polynomial.legendre.leggauss(96)
    X, M_ = np.meshgrid(x, mu_, indexing="ij")
    u = GMk / ell**2 * X / (1 + X**2) ** 1.5
    tot = 0.0
    for sgn in (+1, -1):
        uz = -sgn * u * M_; up = u * np.sqrt(1 - M_**2)
        gz = g + uz; gm = np.sqrt(gz**2 + up**2)
        tot = tot + 0.5 * (f_rule(gm) * gz - f_rule(g) * g)
    wx = 3 * X**2 / (1 + X**2) ** 2.5                  # Plummer-Schalengewicht in x = r/ℓ
    inner = np.sum(tot * wm[None, :], axis=1) / 2
    return np.trapezoid(inner * wx[:, 0], x)
def E_analytisch(GMk, ell, g):
    e = 1e-3 * g
    f1 = (f_rule(g + e) - f_rule(g - e)) / (2 * e); f2 = (f_rule(g + e) - 2 * f_rule(g) + f_rule(g - e)) / e**2
    u2 = 0.0762 * (GMk / ell**2) ** 2                 # ⟨u²⟩_K = 3·½B(5/2,3)·(GM/ℓ²)²
    return u2 * (2 * f1 / 3 + g * f2 / 6), f1, f2
Es = E_eigen(GM_sun, ELL["Sonne"][1], gR0); Ea, f1, f2 = E_analytisch(GM_sun, ELL["Sonne"][1], gR0)
P(f"   Sonne: E = {Es:.3e} m/s² (Quadratur), 2. Ordnung ⟨u²⟩(2f'/3 + g f''/6) = {float(Ea):.3e};  f' = {float(f1):.3e}, f'' = {float(f2):.3e}")
P("   Feldform (Energiefunktional): ∫Q_s(|g + u(x − x_i)|²)d³x hängt bei gleichförmigem g nicht von x_i ab → Eigenkraft 0;")
P("   E ist ein Term der gemittelten örtlichen Regel, nicht der Feldform (dort bleibt nur die ℓ-Kraft der Notiz 13, 9e-19).")

# ============================================================================ 4
H("4   WEP-Differenzen zur Galaxie hin (Lesart R2; Lesart R1 gibt überall 0)")
g_obs = float(a_regel(gR0))
LLR = 1e-13 * GM_sun / AU**2                            # [G] Nordtvedt: |Δa/a| ≲ 1e-13 relativ zur Sonnenbeschleunigung
MIC = 1e-15 * GM_earth / (R_earth + 7.1e5)**2           # [G] MICROSCOPE: η ≲ 1e-15, Bahn 710 km
P(f"   Vergleich: Galaxienbeschleunigung Regel {g_obs:.3e} (gemessen 2,3e-10 [G]); LLR {LLR:.1e} m/s²; MICROSCOPE {MIC:.1e} m/s²")
def zeile(nm1, nm2, var):
    k = 1 if var == "σ̄" else 2
    GM1, l1 = ELL[nm1][0], ELL[nm1][k]; GM2, l2 = ELL[nm2][0], ELL[nm2][k]
    Dt = D_von(l1, TAB0) - D_von(l2, TAB0)
    Ee = E_eigen(GM1, l1, gR0) - E_eigen(GM2, l2, gR0)
    return Dt, Ee
P(f"   {'Paar':28s} {'ℓ-Wahl':6s} {'ΔD Gezeit':>11s} {'ΔE eigen':>11s} {'|Δa|/g_gal':>11s} {'/LLR':>9s} {'/MICRO':>9s}")
PAARE = [("Pt 0,4 kg", "Ti 0,3 kg"), ("Erde", "Mond"), ("Erde", "Sonne"), ("Laborkörper 1 kg", "Erde"), ("Sonne", "Jupiter"), ("Sonne", "Stern 0,1 M☉")]
for a_, b_ in PAARE:
    for var in ("σ̄", "Hill"):
        if var == "Hill" and a_ == "Sonne" and b_ == "Stern 0,1 M☉": continue
        Dt, Ee = zeile(a_, b_, var)
        tot = abs(Dt) + abs(Ee)
        P(f"   {a_+' – '+b_:28s} {var:6s} {Dt:11.2e} {Ee:11.2e} {tot/g_obs:11.2e} {tot/LLR:9.1e} {tot/MIC:9.1e}")
P(f"   Hill-Erde in R2: eigenes Feld GM/ℓ² = {GM_earth/ELL['Erde'][2]**2:.1e} ≫ a₀ in der ganzen Zone → Erde ohne Zusatz, Sonne mit Zusatz:")
P(f"      Δa ≈ (ν−1)g_N = {f_rule(gR0)*gR0:.2e} m/s² Erde gegen Sonne; erzwungene Verschiebung der Erdbahn Δa/Ω² = "
  f"{f_rule(gR0)*gR0/(GM_sun/AU**3):.0f} m (feste galaktische Richtung) – ausgeschlossen [Planetenentfernungen auf m-Niveau, G].")
P(f"   Alte Zahl (Notiz 13, gelesene Dichte 0,01 M☉/pc³): 1,6–2,0e-18 m/s² = 7–9e-9 g.  Neu, Sonne gegen Stern 0,1 M☉ (σ̄, R2): siehe Zeile oben.")

# ============================================================================ 5
H("5   Glättungsschwelle als Phase: Wolke gegen warmes Medium")
P("   Warmes Medium (n < ρ_J): kein Körper, Selbstfeld ungeglättet im Argument von ν → Eigenbindung linear verstärkt,")
P("   Kugelmittel c_e = ν(1 + K/3) − 1, K = dlnν/dlng. Wolke (n > ρ_J): Überschuss über ℓ_c geglättet; in ḡ_N bleibt im")
P("   Wolkeninnern (R < ℓ_c) Boden ρ_J plus Plummer-Kern ρ_J → Eigenbindung 1 + c_e·2ρ_J/ρ̄ (r ≪ ℓ_c, linear).")
for nm, (T_, rJ) in [("R₀ = 8,18 kpc", (T_SONNE, rhoJ0)), ("15 kpc", (T15, rhoJ15))]:
    o = ORTE[nm]; g = o["g"]; R = o["R"]
    TAB = TAB0 if nm.startswith("R₀") else D_tabelle(R)
    e = 1e-3 * g; K = float((np.log(nu(g + e)) - np.log(nu(g - e))) / (np.log(g + e) - np.log(g - e)))
    ce = float(nu(g)) * (1 + K / 3) - 1
    gob = float(a_regel(g)); v = math.sqrt(gob * R)
    P(f"   {nm}: g_N = {g:.3e}, ν = {float(nu(g)):.3f}, K = {K:.3f}, c_e = {ce:.3f}, v = {v/1e3:.1f} km/s, ρ_J = {rJ/m_p/1e6:.2f} m_p/cm³")
    P(f"      {'Wolke':22s} {'n_H':>7s} {'R [pc]':>7s} {'M^c/M':>6s} {'ℓ_c [pc]':>8s} {'D(ℓ_c)':>10s} {'E':>10s} {'Δa/a':>9s} {'Δv [m/s]':>9s} {'Eigenbindung':>12s}")
    for M, nH in [(1e3, 1e3), (1e5, 1e2), (1e6, 1e2), (1e5, 1e4)]:
        rho = 1.36 * nH * 1e6 * m_p
        Rc = (3 * M * M_sun / (4 * math.pi * rho)) ** (1 / 3)
        fc = max(0.0, 1 - rJ / rho); GMc = G * M * M_sun * fc
        lc = (GMc / T_) ** (1 / 3)
        Dt = D_von(lc, TAB); Ee = E_eigen(GMc, lc, g)
        da = abs(Dt) + abs(Ee)
        P(f"      {M:9.0e} M☉ {'':10s} {nH:7.0e} {Rc/pc:7.1f} {fc:6.3f} {lc/pc:8.1f} {Dt:10.2e} {Ee:10.2e} {da/gob:9.1e} {da*R/(2*v):9.2e} {1+ce*2*rJ/rho:12.4f}")
    P(f"      ℓ_c/R = (ρ̄/ρ_J − 1)^(1/3) (homogene Wolke): n_H = 1e2 → {(1.36*1e8*m_p/rJ-1)**(1/3):.2f}, n_H = 1e4 → {(1.36*1e10*m_p/rJ-1)**(1/3):.2f}")
    P(f"      warmes Medium (n_H = 0,1–1): Δa = 0 (kein Körper), Eigenbindung 1 + c_e = {1+ce:.3f}")
P("   Rotationskurven HI gegen CO unterscheiden sich um einige km/s durch Driftgeschwindigkeit und Nichtkreisbewegung [G];")
P("   der Unterschied der Bahnbeschleunigung aus der Regel ist ≤ 1e-4 davon: nicht messbar. Die Phase wirkt auf die Eigenbindung.")

with open(__file__.replace(".py", ".txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
