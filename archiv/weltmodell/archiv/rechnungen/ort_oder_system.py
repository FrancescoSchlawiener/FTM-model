# -*- coding: utf-8 -*-
"""
Ort oder System: Wo zählt ein kollektives System seine Rückläufe?

Eine Gleichung (07): ∂_t u + (Γ + θ)u = D∂_ϖ²u, am Rand ϖ = 0 gezählt (G(0) = 1), projiziert q = e^{−√((Γ + θ)τ_D)}.
Hodge (01): Auf dem Teilgebiet des Kindes ist das Feld der Eltern harmonisch (Hintergrund g_e), das eigene Feld g_i ist Gradient.

Lesung L (örtlich): Rand = Ort jedes Sterns; Γ(x) = |g_i(x) + g_e|/c; θ = θ̃ des Zustands.
Lesung G (System): Rand = Rand des Systems; Γ aus dem eigenen Gradienten allein; der Hintergrund nur über θ̃.
   G-a: am Ort jedes Sterns, Γ(x) = |g_i(x)|/c (Hintergrund nicht im Gefälle)
   G-b: einmal am Rand des Systems, Γ = |g_i(r_Rand)|/c für alle Sterne (ν fest); r_Rand = r½ oder ℓ (Gezeitentiefe ϖ_J)
   G mit θ̃ ≠ 0: s² = s_i² + θ̃τ_D, gleiche Zusatzrate in jeder Richtung = gleichgerichtete Form (aussenfeld_kosmisch).

1  Druckgestützte Systeme: σ_los (Plummer b = R_e, isotrope Jeans-Gleichung, Fenster R < R_e), Eingaben [G] wie glaettung_stetig
2  SPARC (weltmodell/daten/sparc): L mit g_e, G-a, G-b (ν fest am äußersten Punkt)
3  Warp: Knotenpräzession ω_p = Ω − ν_z eines Rings (Miyamoto-Nagai 6·10¹⁰ M☉ [A]) in L (g_e ⟂ und in der Ebene) und G
4  Von Hand
"""
import io, sys, os, math
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import A0, nu, kpc
import sparc_flaechendichte as sf, sparc_universal as su

np.seterr(all="ignore")
GK = 4.30091e-6                                  # G in kpc (km/s)²/M☉
A0K = A0 * kpc / 1e6                             # a₀ in (km/s)²/kpc
MU, WMU = np.polynomial.legendre.leggauss(400)


def nux(x):
    """ν bei g/a₀ = x"""
    return nu(np.maximum(np.asarray(x, float), 1e-300) * A0, A0)


def mono(gi, ge):
    """Monopol von ν(|g_i + g_e|)(g_i + g_e) − A über die Kugel (A fällt heraus: kein Monopol); g in a₀"""
    g = np.atleast_1d(gi)[:, None]; m = MU[None, :]
    gm = np.sqrt(np.maximum(g**2 + ge**2 + 2 * g * ge * m, 0))
    return 0.5 * np.sum(WMU * nux(gm) * (g + ge * m), axis=1)


def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return math.sqrt(max(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r), 0) * A0K)


SYS = [  # Name, M [M☉] (Υ = 2), R_e [kpc], Wirt [M☉], Abstand [kpc], σ gemessen, (+, −)   alles [G]/[A] aus glaettung_stetig
    ("Pal 14", 2e4, 0.0247, 6e10, 71.0, 0.38, (0.12, 0.12)),
    ("Crater II", 3.2e5, 1.066, 6e10, 117.5, 2.7, (0.3, 0.3)),
    ("Draco", 5.8e5, 0.221, 6e10, 76.0, 9.1, (1.2, 1.2)),
    ("Fornax", 4.0e7, 0.71, 6e10, 147.0, 11.7, (0.9, 0.9)),
    ("DF2", 2.0e8, 2.2, 1e11, 80.0, 8.5, (2.3, 3.1)),
    ("DF4", 1.5e8, 1.6, 1e11, 165.0, 4.2, (4.4, 2.2))]

out = []
P = out.append
P("0  Beide Lesungen in der einen Gleichung ∂_t u + (Γ + θ)u = D∂_ϖ²u, stationär, Rand u(0) = 1, u beschränkt:")
P("   L: u_x(ϖ) = e^{−ϖ√((|g_i(x) + g_e|/c + θ̃)/D)},  q(x) = u_x(ϖ_R),  a(x) = ν(x)·(g_i + g_e) − A")
P("   G: u_x(ϖ) = e^{−ϖ√((|g_i(x)|/c + θ̃)/D)};  g_e wirkt nur als gemeinsamer Fall A (Monopol 0) und über θ̃ des Bereichs")
P("   Gebundener Satellit im Bereich des Wirts: θ̃ = 0 (w = 0) ⇒ G = isoliert.  θ̃ > 0 gibt s² = s_i² + θ̃τ_D = gleichgerichtete Form.")
P("")
P(f"1  Druckgestützte Systeme, a₀ = {A0:.2e} m/s²; σ_los(R < R_e) [km/s]; Eingaben und Messwerte [G] aus glaettung_stetig/zwerge_aussenfeld")
P(f"   {'System':10s} {'g_i(r½)':>8s} {'g_e':>7s} {'ℓ[kpc]':>7s} {'Newton':>7s} {'L':>7s} {'G-a':>7s} {'G-b r½':>7s} {'G-b ℓ':>7s}   gemessen [G]      Abstand in σ: L | G-a")
RES = {}
for nm, M, Re, Mh, D, so, er in SYS:
    b = Re; GM = GK * M / A0K                              # kpc²·a₀
    r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, 6000)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5
    rho = (1 + (r / b)**2)**-2.5
    gN = GMenc / r**2
    ge = GK * Mh / D**2 / A0K
    rh = b / math.sqrt(2**(2 / 3) - 1)                     # 3D-Halbmassenradius
    ell = D * (M / (3 * Mh))**(1 / 3)                      # Gezeitentiefe ℓ = ϖ_J (Wirt als Punkt)
    gh = GM / 2 / rh**2
    gl = GM / ell**2
    sN = jeans(r, rho, gN, Re)
    sL = jeans(r, rho, mono(gN, ge), Re)
    sGa = jeans(r, rho, nux(gN) * gN, Re)
    sGr = jeans(r, rho, float(nux(gh)) * gN, Re)
    sGl = jeans(r, rho, float(nux(gl)) * gN, Re)
    e = lambda s: (s - so) / (er[1] if s > so else er[0])
    RES[nm] = dict(N=sN, L=sL, Ga=sGa, Gr=sGr, Gl=sGl, so=so, gh=gh, ge=ge)
    P(f"   {nm:10s} {gh:8.4f} {ge:7.4f} {ell:7.2f} {sN:7.2f} {sL:7.2f} {sGa:7.2f} {sGr:7.2f} {sGl:7.2f}   {so:5.2f} +{er[0]}/−{er[1]}   {e(sL):+6.1f} | {e(sGa):+6.1f}")
P("   (g in a₀; L = Feldform/örtliche Regel mit Monopol; G-a = eigenes Gefälle am Ort; G-b = ν fest vom Rand r½ bzw. ℓ)")
P("   Richtung, die die Messung verlangt, gegen L:  " + ",  ".join(
    f"{k} {'hoch' if v['so'] > v['L'] else 'runter'} (G-a {'hoch' if v['Ga'] > v['L'] + 1e-9 else 'gleich/runter'})" for k, v in RES.items()))
nL = sum(abs(v['L'] - v['so']) <= 2 * max(SYS[i][6]) for i, (k, v) in enumerate(RES.items()))
nG = sum(abs(v['Ga'] - v['so']) <= 2 * max(SYS[i][6]) for i, (k, v) in enumerate(RES.items()))
P(f"   innerhalb 2σ: L {nL} von {len(RES)},  G-a {nG} von {len(RES)}")

# ------------------------------------------------------------------ 2 SPARC
tab = su.tab_voll()
GAL = []
for n, d in tab.items():
    if d['Q'] > 2 or d['inc'] < 30:
        continue
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', n + '_rotmod.dat'), comments='#')
    R, Vo, Vg, Vd, Vb = a[:, 0], a[:, 1], a[:, 3], a[:, 4], a[:, 5]
    Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        continue
    R, Vo, Vb2 = R[m], Vo[m], Vb2[m]
    Vc2 = Vo**2 + 8.0**2 * R / d['Rd']                      # Gasdruck σ = 8 km/s wie sparc_universal
    GAL.append((n, R, Vc2 / R / A0K, Vb2 / R / A0K))


def bewerte(regel):
    res, sl_m, sl_r = [], [], []
    for n, R, go, gn in GAL:
        gp = regel(gn)
        res.append(np.log10(go / gp))
        if len(R) >= 3:
            x = np.log(R[-3:])
            sl_m.append(np.polyfit(x, 0.5 * np.log(go[-3:] * R[-3:]), 1)[0])
            sl_r.append(np.polyfit(x, 0.5 * np.log(gp[-3:] * R[-3:]), 1)[0])
    rr = np.concatenate(res)
    return np.median(rr), np.std(rr), np.median(sl_m), np.median(sl_r), np.median(np.array(sl_m) - np.array(sl_r))


P("")
P(f"2  SPARC (Q ≤ 2, i ≥ 30°, Υ_d = 0,5, Υ_b = 0,7, Gasdruck 8 km/s): {len(GAL)} Galaxien, {sum(len(g[1]) for g in GAL)} Punkte")
P("   Regel                         Rest Median  Streuung [dex] | d ln v/d ln r der äußersten 3 Punkte: gemessen, Regel, Median Δ")
for lab, f in [("L, g_e = 0 (= G-a)", lambda g: nux(g) * g),
               ("L, g_e = 0,01 senkrecht", lambda g: nux(np.hypot(g, 0.01)) * g),
               ("L, g_e = 0,03 senkrecht", lambda g: nux(np.hypot(g, 0.03)) * g),
               ("L, g_e = 0,03 gleichgerichtet", lambda g: nux(g + 0.03) * g),
               ("G-b, ν fest am äußersten Punkt", lambda g: float(nux(g[-1])) * g),
               ("G-b, ν fest am Median-Punkt", lambda g: float(nux(np.median(g))) * g)]:
    m0, sd, sm, sr, dd = bewerte(f)
    P(f"   {lab:30s} {m0:+10.3f} {sd:10.3f}     | {sm:+.3f}  {sr:+.3f}  {dd:+.3f}")
P("   G-b hat dieselbe Form wie Newton (v² ∝ V_bar²): die Außenkurven fallen wie die baryonischen, die Streuung wächst.")

# ------------------------------------------------------------------ 3 Warp
Mw, aw, bw = 6e10, 3.0, 0.3
GMw = GK * Mw


def gMN(R, z):
    """Newton-Feld Miyamoto-Nagai [(km/s)²/kpc], Vektor (g_R, g_z), einwärts negativ"""
    B = math.sqrt(z * z + bw * bw); S = (R * R + (aw + B)**2)**1.5
    return -GMw * R / S, -GMw * z * (aw + B) / (B * S)


def omega_p(R, ge, lage, h=1e-4):
    """Ω − ν_z am Ring R; äußeres Feld ge (in a₀), senkrecht (längs z) oder in der Ebene (Mittel über das Azimut); None = G"""
    gE = ge * A0K
    def a(Rr, z, cphi=1.0):
        gR, gz = gMN(Rr, z)
        if lage == "senkrecht":
            v = np.array([gR, gz + gE]); A = np.array([0.0, float(nux(gE / A0K)) * gE])   # A ≈ Fall als Ganzes (Monopol 0, Ω/ν_z hängen nicht an A)
        elif lage == "Ebene":
            v = np.array([gR + gE * cphi, gz]); A = np.array([0.0, 0.0])
        else:
            v = np.array([gR, gz]); A = np.array([0.0, 0.0])
        return float(nux(np.hypot(*v) / A0K)) * v - A
    if lage == "Ebene":
        ph = np.linspace(0, 2 * np.pi, 73)[:-1]; W, Nz = [], []
        for p in ph:
            aR = a(R, 0, math.cos(p))[0] - gE * math.cos(p) * float(nux(gE / A0K))  # Fall als Ganzes abgezogen (radiale Projektion)
            W.append(math.sqrt(max(-aR / R, 0)))
            Nz.append(math.sqrt(max(-(a(R, h, math.cos(p))[1] - a(R, -h, math.cos(p))[1]) / (2 * h), 0)))
        return float(np.mean(W) - np.mean(Nz)), float(np.mean(W))
    aR = a(R, 0)[0]
    Om = math.sqrt(-aR / R)
    nz = math.sqrt(-(a(R, h)[1] - a(R, -h)[1]) / (2 * h))
    return Om - nz, Om


P("")
P("3  Warp: freie Knotenpräzession eines wenig geneigten Rings, ω_p = Ω − ν_z [km/s/kpc] (> 0 rechtläufig); Messung +10,9 [G]")
P("   Miyamoto-Nagai 6·10¹⁰ M☉, a = 3, b = 0,3 kpc [A] (wie ringe_oertlich); Lokale Gruppe g_e [A]")
P(f"   {'R[kpc]':>6s} {'Newton':>8s} {'G (=L, g_e=0)':>14s} {'L ⟂ 0,01':>9s} {'L ⟂ 0,03':>9s} {'L Ebene 0,03':>13s} {'G-b ν(ℓ_MW)':>12s}")
ellMW = 785.0 * (6e10 / (3 * 1e11))**(1 / 3)              # Milchstraße im Feld von M31 [G/A], wie glaettung_stetig
for R in (12.0, 14.0, 16.0, 20.0):
    gR0 = gMN(R, 0)[0]
    Nw = math.sqrt(-gR0 / R) - math.sqrt(-(gMN(R, 1e-4)[1] - gMN(R, -1e-4)[1]) / 2e-4)
    wG, _ = omega_p(R, 0.0, None)
    w1, _ = omega_p(R, 0.01, "senkrecht")
    w3, _ = omega_p(R, 0.03, "senkrecht")
    wE, _ = omega_p(R, 0.03, "Ebene")
    nb = float(nux(GMw / ellMW**2 / A0K))
    P(f"   {R:6.1f} {Nw:8.2f} {wG:14.2f} {w1:9.2f} {w3:9.2f} {wE:13.2f} {math.sqrt(nb) * Nw:12.2f}")
P(f"   (G-b: ν fest vom Rand ℓ_MW = {ellMW:.0f} kpc, ν = {float(nux(GMw / ellMW**2 / A0K)):.1f}; ω_p = √ν·ω_N)")
P("   Ein homogenes g_e hat über den Ring kein Moment (G: exakt 0); in L bleibt nur l = 2 ∝ g_e² (warp_amplitude: ≤ 0,02 km/s/kpc).")

# ------------------------------------------------------------------ 4 Von Hand
P("")
P("4  Von Hand")
GMsi = 6.6743e-11 * 2e8 * 1.98892e30
s4 = 4 / 81 * GMsi * 1.1e-10
P(f"   DF2 isoliert, tiefer Grenzfall σ⁴ = (4/81)GMa₀: GM = {GMsi:.4e} m³/s², σ⁴ = {s4:.4e} m⁴/s⁴, σ = {s4**0.25 / 1e3:.2f} km/s (G-a Jeans {RES['DF2']['Ga']:.2f}, Fenster R < R_e)")
GMsi = 6.6743e-11 * 3.2e5 * 1.98892e30
s4 = 4 / 81 * GMsi * 1.1e-10
P(f"   Crater II isoliert: σ = {s4**0.25 / 1e3:.2f} km/s (G-a Jeans {RES['Crater II']['Ga']:.2f}); gemessen 2,7 ± 0,3 [G]")
x = 0.0181
P(f"   DF2, L im Grenzfall g_i ≪ g_e: ν(g_e = {x}) = 1/(1 − e^(−√{x})) = 1/(1 − e^(−{math.sqrt(x):.4f})) = {1 / (1 - math.exp(-math.sqrt(x))):.3f}"
  f"; σ ≈ σ_N·√(ν(1 + K/3)), K = −(ν−1)s/2 = {-(1 / (1 - math.exp(-math.sqrt(x))) - 1) * math.sqrt(x) / 2:.3f}"
  f" → {RES['DF2']['N'] * math.sqrt((1 / (1 - math.exp(-math.sqrt(x)))) * (1 - (1 / (1 - math.exp(-math.sqrt(x))) - 1) * math.sqrt(x) / 6)):.2f} km/s (L {RES['DF2']['L']:.2f}; g_i nicht ≪ g_e)")
print("\n".join(out))
