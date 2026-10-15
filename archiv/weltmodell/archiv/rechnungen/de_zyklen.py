# -*- coding: utf-8 -*-
"""
Dunkle Energie: eine Rückstellung je geschlossenem Zyklus der Zelle? (T³: b₁ = 3 ⇒ k_DE = 3·k_s)
Blick: k_DE = ρ_DE·G·c²/a₀² ist ein Koeffizient der Regel des Ganzen (G), keine Energie am Ort.
1 Herkunft der 3 in π⁴/15: Flussdimension (Feld einer Masse, L) gegen Zahl geschlossener Zyklen (Topologie, G);
  was die Normierung von ℒ_𝒮 (∂ℒ/∂Y = 2n_B, s⁴ = ℓ₀²Y) je Regelinstanz gibt; was additive Zyklen am Ort bedeuteten.
2 Gegenlesarten mit ihrer Zahl n = k_DE/k_s und dem a₀, das θ* [G] dafür verlangt.
3 a₀ aus SPARC (eigene Daten) nach der Methode der Gegenprobe zu Notiz 74 (f = 0, Υ fest mit 0,1 dex Fehler je Punkt,
  Abstand, Neigung, errV; Jackknife über Galaxien), Regel ν aus ftm_kern.
4 Zeitverlauf einer Rückstellung je Zyklus.
5 Handprobe.
Aus ftm_kern: c, G, T_CMB, OMEGA_B_H2, Q_KOSM, N_EFF_STD, nu, Q_tilde, n_B, a0_frueh_relativ, OMEGA_M, OMEGA_R.
"""
import math, os, glob
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from ftm_kern import c, G, T_CMB, OMEGA_B_H2, Q_KOSM, N_EFF_STD, nu, Q_tilde, a0_frueh_relativ

OUT = []
def P(s=""):
    OUT.append(s); print(s)

k_s = math.pi**3 / 30
# ---------------------------------------------------------------- θ* → ω_DE (wie de_invariante)
H100 = 100e3 / 3.0856775814913673e22
Z_STERN, THETA_MESS = 1089.9, 1.0411e-2            # [G] Planck
wb = OMEGA_B_H2; wm = wb * (1 + Q_KOSM)
wg = 2.469e-5 * (T_CMB / 2.7255)**4; wr = wg * (1 + 0.2271 * N_EFF_STD)
def E2(z, wde): return wm * (1 + z)**3 + wr * (1 + z)**4 + wde
def theta_stern(wde):
    rs = quad(lambda z: c / math.sqrt(3 * (1 + 0.75 * wb / wg / (1 + z))) / (H100 * math.sqrt(E2(z, wde))), Z_STERN, np.inf, limit=400)[0]
    dm = quad(lambda z: c / (H100 * math.sqrt(E2(z, wde))), 0, Z_STERN, limit=400)[0]
    return rs / dm
wde = brentq(lambda w: theta_stern(w) - THETA_MESS, 0.1, 0.6)
dw = brentq(lambda w: theta_stern(w) - (THETA_MESS + 3e-6), 0.1, 0.6) - wde
rho_de = wde * 3 * H100**2 / (8 * math.pi * G)
def a0_fuer(n, w=wde): return math.sqrt(w * 3 * H100**2 / (8 * math.pi * G) * G * c**2 / (n * k_s))

P("=" * 100); P("1 Die 3 in π⁴/15: Flussdimension, nicht Zahl der Zyklen"); P("=" * 100)
P(f"  Q̃(0) = {Q_tilde(0.0):.6f} = −4Γ(4)ζ(4) = −4π⁴/15 = {-4*math.pi**4/15:.6f}")
P("  Herkunft der Potenz: s⁴ = ℓ₀²Y, Y = |∇σ|² (ein Feld, ein Betrag) ⇒ dY = 4s³ds/ℓ₀², ∂ℒ/∂Y = 2n_B ⇒ dQ̃/ds = 4s³n_B.")
P("  In d Raumdimensionen verlangt die flache Rotationskurve (Tiefbereich ν ≈ 1/s, g_N ∝ r^{1−d}) s ∝ g_N^{(d−2)/(d−1)},")
P("  Y ∝ s^{2(d−1)/(d−2)}, Gewicht x^p n_B mit p = d/(d−2):")
for d in (3, 4, 5, 6):
    p = d / (d - 2)
    I = quad(lambda x: x**p * math.exp(-x) / -math.expm1(-x), 0, np.inf)[0]
    P(f"     d = {d}: p = {p:.3f}, ∫x^p n_B = Γ(p+1)ζ(p+1) = {I:.4f}")
P("  p = 3 nur bei d = 3: die Potenz ist das Flussgesetz des Feldes einer Masse (Gauß über S², L), keine Summe über Richtungen.")
P(f"  Je Richtung (d = 1 der Planck-Leiter) gäbe Γ(2)ζ(2) = {math.pi**2/6:.4f}; 3·π²/6 = {math.pi**2/2:.4f} ≠ π⁴/15 = {math.pi**4/15:.4f}:")
P("  π⁴/15 ist nicht 3 × (Beitrag je Richtung). Die Einwendung „d = 3 zählt die Richtungen schon“ trägt nicht als Zählargument.")
P("  Aber: ℒ_𝒮 enthält (2/ℓ₀²)Q̃(0) als Skalar je Punkt (0-Form mal √−g), ein Y, ein s. Integriert über die Zelle gibt das")
P("  V·Q̃(0) für jede Topologie und jede Form: b₁ kommt in ℒ_𝒮 nicht vor. Je Regelinstanz (eigenes Y_i, gleiches ℓ₀) gibt die")
P("  Normierung genau einmal Q̃(0), also k_s. b₁·k_s verlangt b₁ Regelinstanzen oder einen Term der Zelle (G).")
P("  (a) b₁ volle Instanzen am Ort: ν − 1 → b₁·n_B, im Tiefbereich a = b₁·√(a₀g_N): gelesenes a₀ = b₁²·a₀ = 9a₀ – ausgeschlossen.")
P("  (b) je Zyklus eine Komponente, ℒ = Σ_i Q̃(s_i), s_i⁴ = ℓ₀²B_i²: Tiefbereich |a|² = a₀g_N·Σ|n_i| ⇒ a₀_eff = a₀Σ|n_i|:")
rng = np.random.default_rng(1)
n = rng.normal(size=(200000, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
sa = np.abs(n).sum(1); la = 0.5 * np.log10(sa)            # Verschiebung von log a im Tiefbereich
P(f"      Achse 1, Flächendiagonale √2 = {math.sqrt(2):.3f}, Raumdiagonale √3 = {math.sqrt(3):.3f}; zufällige Lage: log a streut")
P(f"      um {la.std():.3f} dex (an SPARC kaum trennbar), Mittel +{la.mean():.3f} dex: gelesen würde a₀_eff = {10**(2*la.mean()):.2f}·a₀.")
P(f"      n = 3 verlangte dann an SPARC {1.0733*10**(2*la.mean()):.2f}·10⁻¹⁰ statt 1,07. In einer Scheibe schwankt v⁴ mit dem Azimut bis √2 (m = 4).")
P("      Vor allem widerspricht es der Regel selbst: Y ist ein Skalar, der Verlust der Betrag eines Vektors (s⁴ = s_i⁴ + s_e⁴),")
P("      die Regel am Ort isotrop; Σ_i Q̃(s_i) legte Vorzugsachsen der Zelle in jede Galaxie (Beobachterfehler: G als L gelesen).")
P("  (c) Term der Zelle: Σ_{i≤b₁}|ω_i|² mit harmonischen 1-Formen ω_i (⟨|ω_i|²⟩ = 1). Auf dem flachen T³ ist ω_i = dx_i und")
P("      Σ|ω_i|² = h^{ij}h_ij = 3 an jedem Punkt: am Ort nicht von 3 Raumrichtungen zu unterscheiden, für jede Form des T³ gleich.")
P("      Das ist die einzige widerspruchsfreie Fassung der Lesart; sie folgt nicht aus ℒ_𝒮, sie ist eine Setzung (G).")

P(); P("=" * 100); P("2 Lesarten: n = k_DE/k_s und verlangtes a₀ (θ* fest, σ aus ±0,0003 in 100θ*)"); P("=" * 100)
P(f"  ω_DE = {wde:.4f} ± {abs(dw):.4f};  a₀(n) = cH₁₀₀·√(30ω_DE/(8π⁴n))")
LES = [("Volumenmode b₃ = 1 (nur Vakuumteil)", 1), ("nicht-orientierbar flach, b₁ = 2", 2),
       ("b₁ = 3 (T³, jede Form: kubisch, gestreckt, schief)", 3), ("b₂ = 3 (T³; Poincaré b₂ = b₁)", 3),
       ("Vakuumteil + je Zyklus, 1 + b₁ = 4", 4), ("T⁴ mit innerer Richtung, b₁ = 4", 4),
       ("T⁴, b₂ = 6", 6), ("Σ Betti T³ = 8", 8)]
for lab, nn in LES:
    a = a0_fuer(nn); da = abs(a0_fuer(nn, wde + dw) - a)
    P(f"  {lab:52s} n = {nn}:  a₀ = {a*1e10:.4f} ± {da*1e10:.4f}")
P("  Flache Zellen mit Drehung (b₁ = 1: Halb-, Viertel-, Drittel-, Sechsteldrehung) n = 1 wie Volumenmode; Hantzsche-Wendt b₁ = 0: keine.")

# ---------------------------------------------------------------- 3 SPARC
P(); P("=" * 100); P("3 a₀ aus SPARC (Gegenprobe Notiz 74: f = 0, Υ fest, Fehler je Punkt, Jackknife über Galaxien)"); P("=" * 100)
D0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "daten", "sparc")
KPC = 3.0857e19; GK = 1e6 / KPC
tab = {}
for line in open(os.path.join(D0, "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
    p = line.split()
    if len(p) > 17 and p[1].isdigit() and len(p[0]) > 1:
        try: tab[p[0]] = dict(D=float(p[2]), eD=float(p[3]), fD=int(p[4]), inc=float(p[5]), einc=float(p[6]), Q=int(p[17]))
        except ValueError: pass
gal = []
for fn in sorted(glob.glob(os.path.join(D0, "Rotmod_LTG", "*_rotmod.dat"))):
    name = os.path.basename(fn).replace("_rotmod.dat", ""); t = tab.get(name)
    if t is None or t["Q"] > 2 or t["inc"] < 30: continue
    d = np.loadtxt(fn); m = (d[:, 1] > 0) & (d[:, 0] > 0) & (d[:, 2] > 0)
    gal.append(dict(name=name, prim=t["fD"] >= 2, R=d[m, 0], Vo=d[m, 1], eV=d[m, 2], Vg=d[m, 3], Vd=d[m, 4], Vb=d[m, 5],
                    eD=t["eD"] / t["D"], ei=math.radians(t["einc"]), inc=math.radians(t["inc"])))
P(f"  Galaxien {len(gal)} (Primär-Abstände {sum(g['prim'] for g in gal)}), Punkte {sum(len(g['R']) for g in gal)}; Υ_b = 0,7")
AG = np.exp(np.linspace(np.log(0.6e-10), np.log(2.0e-10), 281))
def lna(g, Yd, a0):
    gN = np.maximum((g["Vg"] * np.abs(g["Vg"]) + Yd * g["Vd"] * np.abs(g["Vd"]) + 0.7 * g["Vb"]**2) / g["R"] * GK, 1e-14)
    return np.log(nu(gN, a0) * gN)
def tafel(Yd):
    T = np.zeros((len(gal), len(AG)))
    for i, g in enumerate(gal):
        dat = np.log(g["Vo"]**2 / g["R"] * GK)
        dY = (lna(g, Yd * np.exp(0.23), 1.2e-10) - lna(g, Yd * np.exp(-0.23), 1.2e-10)) / 2
        s2 = (2 * g["eV"] / g["Vo"])**2 + g["eD"]**2 + (2 * g["ei"] / np.tan(g["inc"]))**2 + dY**2
        for k, a0 in enumerate(AG):
            T[i, k] = np.sum((dat - lna(g, Yd, a0))**2 / s2)
    return T
def a0_best(C):
    k = int(np.argmin(C))
    if 0 < k < len(AG) - 1:
        y0, y1, y2 = C[k-1], C[k], C[k+1]; h = np.log(AG[1] / AG[0])
        return math.exp(math.log(AG[k]) + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2)))
    return AG[k]
prim = np.array([g["prim"] for g in gal]); alle = np.ones(len(gal), bool)
res = {}
for Yd in (0.4, 0.5, 0.6):
    T = tafel(Yd)
    for lab, sel in (("alle", alle), ("Primär", prim), ("Hubble", ~prim)):
        idx = np.where(sel)[0]; a = a0_best(T[sel].sum(0))
        jk = np.array([a0_best(T[idx[idx != i]].sum(0)) for i in idx]); nn = len(idx)
        se = math.sqrt((nn - 1) / nn * np.sum((jk - jk.mean())**2))
        res[(Yd, lab)] = (a, se)
        P(f"  Υ_d = {Yd}  {lab:7s}: a₀ = {a*1e10:.3f} ± {se*1e10:.3f} (Jackknife)  ⇒ n = k_DE/k_s = {rho_de*G*c**2/(a**2*k_s):.2f}")
a5, s5 = res[(0.5, "alle")]
syY = (res[(0.4, "alle")][0] - res[(0.6, "alle")][0]) / 2
syD = abs(res[(0.5, "Primär")][0] - res[(0.5, "Hubble")][0]) / 2
tot = math.sqrt(s5**2 + syY**2 + syD**2)
P(f"  Lesung Υ_d = 0,5: a₀ = {a5*1e10:.3f} ± {s5*1e10:.3f} (stat.) ± {syY*1e10:.3f} (Υ_d ± 0,1) ± {syD*1e10:.3f} (Abstände)"
  f" ⇒ gesamt ± {tot*1e10:.3f} ({tot/a5*100:.0f} %)")
nS = rho_de * G * c**2 / (a5**2 * k_s)
P(f"  n = k_DE/k_s = {nS:.2f} ± {2*tot/a5*nS:.2f} (σ_n/n = 2σ_a/a).  Abstände zur Lesung in σ:")
for nn in (1, 2, 3, 4, 6, 8):
    P(f"     n = {nn}: a₀ = {a0_fuer(nn)*1e10:.3f}, Abstand {(a0_fuer(nn)-a5)/tot:+.1f}σ")
P("  Υ_d bestimmt a₀ fast allein (d ln a₀/d ln Υ_d ≈ −1 im Tiefbereich); ohne eigenes Maß für Υ_d bleibt a₀ auf ≈ 15–20 %.")
P("  Zum Trennen von n = 3 gegen n = 4 (Abstand 14 %) auf 2σ braucht es σ(a₀)/a₀ ≈ 7 %, gegen jede Zahl in 2,4–3,2 etwa 1 %.")

# ---------------------------------------------------------------- 4 Zeit
P(); P("=" * 100); P("4 Zeitverlauf (a₀(z) aus a0_frueh_relativ; Rückstellung je Zyklus trägt ℓ₀ = c²/a₀ wie der Vakuumteil)"); P("=" * 100)
P("  b₁ ist topologisch: das Atmen skaliert L₁, L₂, L₃ gleich, die Zahl der Zyklen bleibt. k_DE(z) = 3k_s·(a₀(z)/a₀)².")
for z in (0, 10, 100, 1090, 1e4, 4e8):
    f = a0_frueh_relativ(z)**2
    ratio_r = wde * f / (wr * (1 + z)**4); ratio_m = wde * f / (wm * (1 + z)**3)
    dN = wde * f / (wg * 0.2271 * (1 + z)**4)
    P(f"  z = {z:9.0f}: (a₀(z)/a₀)² = {f:.4f},  ρ_DE/ρ_m = {ratio_m:.2e},  ρ_DE/ρ_r = {ratio_r:.2e},  als ΔN_eff = {dN:.1e} (nur früh sinnvoll)")
P("  Wie der Vakuumteil: zeitfest seit z ≈ 10 (w = −1), früh verschwindend, N_eff unberührt. Mit Λ_i ∝ L_i (Schluss an der")
P("  Raumschleife) liefe jeder Zyklusterm wie a⁻²: w = −⅓, heute als Krümmung Ω ≈ 0,69 (gemessen 0,0007 ± 0,0037): nur Λ fest trägt.")

P(); P("=" * 100); P("5 Handprobe"); P("=" * 100)
P("  a₀² = ω_DE·3H₁₀₀²/(8πG)·Gc²/(3·π³/30) = ω_DE·30·(cH₁₀₀)²/(8π⁴)")
P("  cH₁₀₀ = 2,99792e8·3,24078e-18 = 9,7156e-10;  30·0,3170/(8·97,409) = 9,510/779,27 = 0,012204;  √ = 0,11047")
P(f"  a₀ = 9,7156e-10·0,11047 = 1,0733e-10;  Code {a0_fuer(3)*1e10:.4f}e-10")
P(f"  Deep-Regel Handwert: g_N = 0,01a₀ ⇒ s = 0,1, ν = 1/(1 − e^{{−0,1}}) = 1/0,095163 = 10,508;  Kern {float(nu(0.011e-10, 1.1e-10)):.3f}")

with open(__file__.replace(".py", ".txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
