# -*- coding: utf-8 -*-
"""
Schwelle Galaxie für Galaxie an SPARC mit der örtlichen Regel (radiale Form, f = 0) aus dem Kern.
Frage: Steigt das je Galaxie gelesene a₀ mit der zentralen Flächendichte Σ₀ (Notiz 68: +0,117 ± 0,042 je Dekade,
nur mit Geschwindigkeitsfehlern)? Hier je Galaxie Abstand, Neigung und Υ_d als Störgrößen mit ihren Fehlern
(Abstand und Neigung aus der SPARC-Tabelle, Υ_d 0,1 dex) profiliert; Steigung mit innerer Streuung und Bootstrap.
"""
import io, sys, glob, os, math
import numpy as np
from scipy.optimize import least_squares
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import a_regel, G, M_sun, pc, kpc

D0 = r"C:\DEV\FTM\weltmodell\daten\sparc"
GK = 1e6 / kpc                                     # (km/s)²/kpc → m/s²
tab = {}
for line in open(os.path.join(D0, "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
    p = line.split()
    if len(p) > 17 and p[1].isdigit() and len(p[0]) > 1:
        try:
            tab[p[0]] = dict(D=float(p[2]), eD=float(p[3]), fD=int(p[4]), inc=float(p[5]), einc=float(p[6]),
                             SB0=float(p[12]), Q=int(p[17]))
        except ValueError:
            pass

gal = []
for fn in sorted(glob.glob(os.path.join(D0, "Rotmod_LTG", "*_rotmod.dat"))):
    name = os.path.basename(fn).replace("_rotmod.dat", "")
    t = tab.get(name)
    if t is None or t["Q"] > 2 or t["inc"] < 30 or t["SB0"] <= 0: continue
    d = np.loadtxt(fn)
    R, Vo, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    m = (Vo > 0) & (R > 0) & (eV > 0)
    if m.sum() < 5: continue
    gal.append(dict(name=name, prim=t["fD"] >= 2, R=R[m], Vo=Vo[m], eV=eV[m], Vg=Vg[m], Vd=Vd[m], Vb=Vb[m],
                    eD=t["eD"] / t["D"], inc=math.radians(t["inc"]), ei=math.radians(max(t["einc"], 1.0)), SB0=t["SB0"]))

Y0, sY = 0.5, math.log(10) * 0.1
AG = np.exp(np.linspace(math.log(0.3e-10), math.log(4e-10), 81))

def resid(par, g, a0, stoer=True):
    uD, di, uY = par if stoer else (0.0, 0.0, 0.0)
    Yd = Y0 * math.exp(uY)
    gbar = np.maximum((g["Vg"] * np.abs(g["Vg"]) + Yd * g["Vd"] * np.abs(g["Vd"]) + 0.7 * g["Vb"]**2) / g["R"] * GK, 1e-15)
    inc = min(max(g["inc"] + di, 0.2), math.pi / 2)
    gobs = np.log(g["Vo"]**2 / g["R"] * GK) - uD - 2 * math.log(math.sin(g["inc"]) / math.sin(inc))
    r = (gobs - np.log(a_regel(gbar, a0))) / (2 * g["eV"] / g["Vo"])
    return np.r_[r, uD / g["eD"], di / g["ei"], uY / sY] if stoer else r

def profil(g, stoer=True):
    c2 = []
    for a0 in AG:
        if stoer:
            s = least_squares(resid, [0.0, 0.0, 0.0], args=(g, a0))
            c2.append(2 * s.cost)
        else:
            c2.append(np.sum(resid(None, g, a0, False)**2))
    c2 = np.array(c2); k = int(np.argmin(c2))
    if k in (0, len(AG) - 1): return None
    x = np.log(AG); y0, y1, y2 = c2[k - 1], c2[k], c2[k + 1]; h = x[1] - x[0]
    curv = (y0 - 2 * y1 + y2) / h**2
    xb = x[k] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2))
    sig = math.sqrt(2 / curv) if curv > 0 else 1.0                 # Δχ² = 1
    red = max(c2[k] / max(len(g["R"]) - 1, 1), 1.0)                 # Streckung bei χ²/FG > 1
    return xb, sig * math.sqrt(red)

def steigung(x, y, s):
    """ln a₀ = α + β·x mit Messfehler s und innerer Streuung τ (Maximum-Likelihood)."""
    best = None
    for tau in np.linspace(0, 1.0, 201):
        w = 1 / (s**2 + tau**2)
        A = np.vstack([np.ones_like(x), x]).T
        C = np.linalg.inv(A.T @ (A * w[:, None])); ab = C @ (A.T @ (w * y))
        L = -0.5 * np.sum(w * (y - A @ ab)**2) + 0.5 * np.sum(np.log(w))
        if best is None or L > best[0]: best = (L, tau, ab, C)
    _, tau, ab, C = best
    return ab[1], math.sqrt(C[1, 1]), tau, ab[0]

rng = np.random.default_rng(1)
VARIANTEN = (("nur Geschwindigkeitsfehler (wie Notiz 68)", False, sY),
             ("Abstand und Neigung als Störgrößen, Υ_d fest", True, 1e-4),
             ("Abstand, Neigung, Υ_d (0,1 dex) als Störgrößen", True, sY))
for lab, stoer, sYv in VARIANTEN:
    sY = sYv
    rows = []
    for g in gal:
        r = profil(g, stoer)
        if r: rows.append((g["name"], g["prim"], math.log(Y0 * g["SB0"]), r[0], r[1]))
    x = np.array([r[2] for r in rows]) / math.log(10); y = np.array([r[3] for r in rows]); s = np.array([r[4] for r in rows])
    pr = np.array([r[1] for r in rows])
    b, sb, tau, a = steigung(x, y, s)
    bs = []
    for _ in range(1000):
        i = rng.integers(0, len(x), len(x)); bs.append(steigung(x[i], y[i], s[i])[0])
    bs = np.array(bs)
    w = 1 / (s**2 + tau**2); a0m = math.exp(np.sum(w * y) / np.sum(w))
    print(f"{lab}: {len(rows)} Galaxien")
    print(f"  d ln a₀ / d log Σ₀ = {b:+.3f} ± {sb:.3f} (Fehler)  ± {bs.std():.3f} (Bootstrap)  → d log a₀/d log Σ₀ = {b/math.log(10):+.3f} ± {bs.std()/math.log(10):.3f}")
    print(f"  innere Streuung von ln a₀: {tau:.2f} ({tau/math.log(10):.2f} dex), gewichtetes Mittel a₀ = {a0m*1e10:.2f}·10⁻¹⁰")
    for sub, sel in (("Primärabstände", pr), ("Hubble-Fluss", ~pr)):
        bb, ss, tt, _ = steigung(x[sel], y[sel], s[sel])
        print(f"  {sub:15s}: d log a₀/d log Σ₀ = {bb/math.log(10):+.3f} ± {ss/math.log(10):.3f} ({sel.sum()} Galaxien)")
