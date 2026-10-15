# -*- coding: utf-8 -*-
"""Äußeres Feld und angepasste Schwelle mit den Punkten der SPARC-Kurven: Daten mit ν(|g_N + g_e|)·g_N (senkrecht/gleichgerichtet),
angepasst ohne äußeres Feld. Gibt a₀′/a₀ (mit Geschwindigkeitsfehlern gewichtet und ungewichtet)."""
import io, sys, math, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from ftm_kern import nu, A0
from scipy.optimize import minimize_scalar
src = open(r"C:\DEV\FTM\weltmodell\rechnungen\gewicht_scheibenfeld_gegenprobe.py", encoding="utf-8").read()
src = src[:src.index("def nu(y)")]
g_ = {}; exec(src, g_); gal, GK = g_["gal"], g_["GK"]
a0 = A0
pts = []
for g in gal:
    gb = np.maximum((g["Vg"]*np.abs(g["Vg"]) + 0.5*g["Vd"]*np.abs(g["Vd"]) + 0.7*g["Vb"]**2)/g["R"]*GK, 1e-15)
    w = (2*g["eV"]/g["Vo"])**-2
    pts.append((gb, w))
gb = np.concatenate([p[0] for p in pts]); w = np.concatenate([p[1] for p in pts])
for ge in (0.01, 0.03, 0.05):
    for lage in ("senkrecht", "gleich"):
        gt = np.hypot(gb, ge*a0) if lage == "senkrecht" else gb + ge*a0
        obs = np.log(nu(gt, a0)*gb)
        f = lambda la: np.sum(w*(obs - np.log(nu(gb, a0*math.exp(la))*gb))**2)
        r = minimize_scalar(f, bounds=(-3, 1), method="bounded")
        # ungewichtet als zweiter Weg
        f2 = lambda la: np.sum((obs - np.log(nu(gb, a0*math.exp(la))*gb))**2)
        r2 = minimize_scalar(f2, bounds=(-3, 1), method="bounded")
        print(f"g_e = {ge:.2f} a₀ {lage:9s}: a₀'/a₀ = {math.exp(r.x):.2f} (gewichtet), {math.exp(r2.x):.2f} (ungewichtet)")
