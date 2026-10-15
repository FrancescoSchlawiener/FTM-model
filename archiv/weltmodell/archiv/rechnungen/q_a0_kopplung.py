# -*- coding: utf-8 -*-
"""Kopplung von Q und a₀ über die Zelle.

Gleichung der Zelle (03): ν·ρ_b·V = 3c⁴/(2G·a₀), ν = 1 + Q.  ρ_b·V = M_b (G, mitbewegt fest), a₀ = c²/2Λ (I, seit z ≈ 10 fest).
Untere Grenze (Schluss):  r*(M_b) ≤ Λ  ⇔  ν ≥ 6, Q ≥ 5.
Obere Grenze (keine Kreispaare, heutiger Beobachter):  L_min ≥ 0,985·d_LSS, L_min = f·V^⅓ (f je Zellform, bestes Glied der Familie)
    ⇒  ν ≤ K·f³/a₀,   K = 3c⁴/(2Gρ_b·(0,985·d_LSS)³).
Zellformen f³: Würfel 1, sechseckig (Höhe = Seite) 2/√3, raumzentriert 3√3/4, flächenzentriert √2 (dichteste Gitterpackung).
SPARC: a₀ global mit der Regel a = a_N/(1 − e^{−√(a_N/a₀)}); je Galaxie Abstand, Neigung und Υ_d (0,1 dex) als Störgrößen
profiliert; Υ_d-Mittel 0,4 / 0,5 / 0,6 (Υ_b = 1,4·Υ_d); Fehler aus Bootstrap über Galaxien.
"""
import io, sys, os, glob, math
import numpy as np
from scipy.optimize import least_squares
import ftm_kern as K

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
c, G, Gpc = K.c, K.G, K.Mpc * 1e3
H0 = K.H0_si()
rho_b = K.OMEGA_B * 3 * H0**2 / (8 * math.pi * G)
dLSS = K.d_LSS()
Lmin = K.KREISPAAR_GRENZE * dLSS
KZ = 3 * c**4 / (2 * G * rho_b * Lmin**3)          # ν·a₀ ≤ KZ·f³
FORMEN = (("Würfel", 1.0), ("sechseckig", 2 / math.sqrt(3)), ("raumzentriert", 3 * math.sqrt(3) / 4), ("flächenzentriert", math.sqrt(2)))
Q_M, dQ_M = 5.36, 0.07


def linie(ch="="):
    print(ch * 112)


def q_max(a0, f3):
    return KZ * f3 / a0 - 1


def a0_max(Q, f3):
    return KZ * f3 / (1 + Q)


linie()
print("0  Grundgrößen")
linie()
print(f"   ρ_b = {rho_b:.4e} kg/m³, d_LSS = {dLSS/Gpc:.3f} Gpc, L_min = {Lmin/Gpc:.3f} Gpc,  K = 3c⁴/(2Gρ_b L_min³) = {KZ:.4e} m/s²")
print("   Zellformen: " + ", ".join(f"{n} f³ = {f3:.4f}" for n, f3 in FORMEN))

# ------------------------------------------------------------------ 1
linie()
print("1  Erlaubtes Gebiet in (a₀, Q): 5 ≤ Q ≤ K·f³/a₀ − 1   ('leer' = obere Grenze unter 5)")
linie()
print("   a₀/10⁻¹⁰ | " + " | ".join(f"{n:>16s}" for n, _ in FORMEN))
for a in np.arange(1.00, 1.2501, 0.01):
    zeile = []
    for _, f3 in FORMEN:
        qm = q_max(a * 1e-10, f3)
        zeile.append(f"   5 ≤ Q ≤ {qm:5.2f}" if qm >= 5 else f"{'leer':>16s}")
    print(f"     {a:5.2f}   | " + " | ".join(zeile))
print("   Fenster schließt (Q_max = 5) bei a₀ = K·f³/6: " + ", ".join(f"{n} {KZ*f3/6*1e10:.4f}" for n, f3 in FORMEN) + " ·10⁻¹⁰")

# ------------------------------------------------------------------ 2a
linie()
print(f"2a Q = {Q_M} ± {dQ_M} [G] eingesetzt: a₀ ≤ K·f³/(1 + Q)")
linie()
for n, f3 in FORMEN:
    print(f"   {n:17s} a₀ ≤ {a0_max(Q_M, f3)*1e10:.4f}  (Q − 1σ: {a0_max(Q_M-dQ_M, f3)*1e10:.4f}, Q + 1σ: {a0_max(Q_M+dQ_M, f3)*1e10:.4f}, "
          f"Q + 2σ: {a0_max(Q_M+2*dQ_M, f3)*1e10:.4f}) ·10⁻¹⁰")

# ------------------------------------------------------------------ 2b SPARC
linie()
print("2b a₀ aus SPARC (eigene Daten), Regel a = a_N/(1 − e^{−√(a_N/a₀)})")
linie()
D0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "daten", "sparc")
GK = 1e6 / K.kpc
tab = {}
for line in open(os.path.join(D0, "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
    p = line.split()
    if len(p) > 17 and p[1].isdigit() and len(p[0]) > 1:
        try:
            tab[p[0]] = dict(D=float(p[2]), eD=float(p[3]), inc=float(p[5]), einc=float(p[6]), Q=int(p[17]))
        except ValueError:
            pass
gal = []
for fn in sorted(glob.glob(os.path.join(D0, "Rotmod_LTG", "*_rotmod.dat"))):
    name = os.path.basename(fn).replace("_rotmod.dat", "")
    t = tab.get(name)
    if t is None or t["Q"] > 2 or t["inc"] < 30:
        continue
    d = np.loadtxt(fn)
    R, Vo, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    m = (Vo > 0) & (R > 0) & (eV > 0)
    if m.sum() < 5:
        continue
    gal.append(dict(R=R[m], Vo=Vo[m], eV=eV[m], Vg=Vg[m], Vd=Vd[m], Vb=Vb[m], eD=max(t["eD"] / t["D"], 0.01),
                    inc=math.radians(t["inc"]), ei=math.radians(max(t["einc"], 1.0))))
npkt = sum(len(g["R"]) for g in gal)
print(f"   {len(gal)} Galaxien (Q ≤ 2, Neigung ≥ 30°), {npkt} Punkte")
sY = math.log(10) * 0.1
AG = np.linspace(0.85e-10, 1.45e-10, 61)


SIG_INT = 0.0                                                            # innere Streuung in ln g (Variante)


def resid(par, g, a0, Y0):
    uD, di, uY = par
    Yd = Y0 * math.exp(uY)
    gbar = (g["Vg"] * np.abs(g["Vg"]) + Yd * g["Vd"] * np.abs(g["Vd"]) + 1.4 * Yd * g["Vb"]**2) / g["R"] * GK
    gbar = np.maximum(gbar, 1e-14)                                      # a_N hängt nicht an D, g_obs ∝ 1/D
    inc = min(max(g["inc"] + di, 0.2), math.pi / 2)
    gobs = np.log(g["Vo"]**2 / g["R"] * GK) - uD + 2 * math.log(math.sin(g["inc"]) / math.sin(inc))
    r = (gobs - np.log(K.a_regel(gbar, a0))) / np.sqrt((2 * g["eV"] / g["Vo"])**2 + SIG_INT**2)
    return np.r_[r, uD / g["eD"], di / g["ei"], uY / sY]


def chi2_tafel(Y0):
    T = np.zeros((len(gal), len(AG)))
    for i, g in enumerate(gal):
        x = np.zeros(3)
        for j, a0 in enumerate(AG):
            s = least_squares(resid, x, args=(g, a0, Y0))
            x = s.x
            T[i, j] = 2 * s.cost
    return T


def bestwert(c2):
    k = int(np.argmin(c2))
    k = min(max(k, 1), len(AG) - 2)
    y0, y1, y2 = c2[k - 1:k + 2]
    h = AG[1] - AG[0]
    return AG[k] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2))


rng = np.random.default_rng(7)
ergeb = {}
for SIG_INT in (0.0, 0.10, 0.20):
  print(f"   innere Streuung σ_int = {SIG_INT:.2f} in ln g ({SIG_INT/math.log(10):.3f} dex):")
  for Y0 in (0.4, 0.5, 0.6):
    T = chi2_tafel(Y0)
    c2 = T.sum(0)
    ab = bestwert(c2)
    boot = np.array([bestwert(T[rng.integers(0, len(gal), len(gal))].sum(0)) for _ in range(2000)])
    sb = boot.std()
    ergeb[(SIG_INT, Y0)] = (ab, sb)
    print(f"     Υ_d = {Y0}: a₀ = {ab*1e10:.4f} ± {sb*1e10:.4f} ·10⁻¹⁰ (Bootstrap), χ²/Punkt = {c2.min()/npkt:.2f}")
SIG_W = 0.10                                                              # Wahl: χ²/Punkt ≈ 1
ergeb = {Y0: v for (si, Y0), v in ergeb.items() if si == SIG_W}
print(f"   weiter mit σ_int = {SIG_W} (χ²/Punkt nahe 1)")
lo = min(a - s for a, s in ergeb.values())
hi = max(a + s for a, s in ergeb.values())
lo2 = min(a - 2 * s for a, s in ergeb.values())
hi2 = max(a + 2 * s for a, s in ergeb.values())
print(f"   Band über Υ_d 0,4–0,6: 1σ {lo*1e10:.3f}–{hi*1e10:.3f}, 2σ {lo2*1e10:.3f}–{hi2*1e10:.3f} ·10⁻¹⁰")

linie()
print("2c Q-Fenster je Zellform aus dem SPARC-Band (obere Grenze am unteren/oberen Rand des Bandes)")
linie()
for n, f3 in FORMEN:
    for tag, (a1, a2) in (("1σ", (lo, hi)), ("2σ", (lo2, hi2))):
        qa, qb = q_max(a2, f3), q_max(a1, f3)
        print(f"   {n:17s} {tag}: Q ≤ {qa:5.2f} (a₀ = {a2*1e10:.3f}) bis Q ≤ {qb:5.2f} (a₀ = {a1*1e10:.3f});  "
              f"mit Q ≥ 5: " + ("leer am oberen Rand" if qa < 5 else "offen") + f";  Q = {Q_M}: " +
              ("verträglich im ganzen Band" if qa >= Q_M else ("nur für a₀ ≤ " + f"{a0_max(Q_M, f3)*1e10:.3f}" if qb >= Q_M else "ausgeschlossen")))
for Y0, (ab, sb) in ergeb.items():
    s_w = (ab - a0_max(Q_M, 1.0)) / math.hypot(sb, a0_max(Q_M, 1.0) * dQ_M / (1 + Q_M))
    print(f"   Würfel, Υ_d = {Y0}: a₀_SPARC − a₀,max(Q) = {(ab - a0_max(Q_M, 1.0))*1e10:+.4f}·10⁻¹⁰ = {s_w:+.2f}σ (beide Fehler)")

# ------------------------------------------------------------------ 4
linie()
print("4  Zeitfestigkeit: ν·a₀ = 3c⁴/(2G·M_b/V) mit M_b, V mitbewegt fest, a₀(z)/a₀(0) aus Λ(z)")
linie()
for z in (0, 1, 10, 1090):
    print(f"   z = {z:5d}: a₀(z)/a₀(0) = {K.a0_frueh_relativ(z)/K.a0_frueh_relativ(0):.4f};  d_LSS des Beobachters bei z wäre kleiner, die Grenze L ≥ 0,985·d_LSS(heute) bleibt mitbewegt")
print("   L_min ist eine Beobachtung von heute an einer mitbewegt festen Zelle: die Grenze gilt für jede Epoche; sie verschärft sich,")
print("   wenn ein späterer Beobachter mit größerem d_LSS keine Kreise sieht.")

linie()
print("5  Handwert")
linie()
print(f"   L_min = {Lmin:.4e} m, L³ = {Lmin**3:.4e}, Gρ_bL³ = {G*rho_b*Lmin**3:.4e}, 3c⁴/2 = {1.5*c**4:.4e}, K = {KZ:.4e}, K/6 = {KZ/6:.4e}")
