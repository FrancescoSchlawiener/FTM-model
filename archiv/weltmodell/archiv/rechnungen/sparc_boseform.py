# -*- coding: utf-8 -*-
"""
Die Bose-Form der Regel direkt an den SPARC-Punkten (Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Υ_b = 0,7).
1) Linearisierung: ν = 1 + n_B(s) ⇔ y := ln(g/(g − g_N)) = s = √(g_N/a₀). Gerade durch den Ursprung, Steigung 1/√a₀.
   Gegenprobe: dieselbe Größe für die „einfache“ und die „Standard“-Form der Literatur (sie krümmen die Gerade).
2) Wien-Stelle der Stufe 2: Maximum der Zusatzbeschleunigung (g − g_N)/a₀ über g_N; Regel: 0,648 bei g_N = 2,54·a₀.
Gebinnt in log g_N, Median je Bin (robust gegen einzelne Punkte).
"""
import math, os
import numpy as np
from scipy.optimize import curve_fit, brentq
from ftm_kern import A0, nu, pc
import sparc_flaechendichte as sf
import sparc_universal as su

kpc = 1e3 * pc


def punkte(ups_d=0.5, sig=0.0):
    t = su.tab_voll(); G_, N_ = [], []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        Vb2 = Vg * np.abs(Vg) + ups_d * Vd**2 + 0.7 * Vb**2
        m = (Vb2 > 0) & (Vo > 0)
        G_.append((Vo[m]**2 + (sig * 1e3)**2 * (R[m] / kpc) / d['Rd']) / R[m]); N_.append(Vb2[m] / R[m])
    return np.concatenate(G_), np.concatenate(N_)


def bins(gN, y, lo=-2.0, hi=1.4, n=17, stat=np.median):
    e = np.linspace(lo, hi, n + 1); x = np.log10(gN / A0); out = []
    for a, b in zip(e[:-1], e[1:]):
        m = (x >= a) & (x < b) & np.isfinite(y)
        if m.sum() >= 20:
            out.append((0.5 * (a + b), stat(y[m]), m.sum(), 1.2533 * y[m].std() / math.sqrt(m.sum())))
    return np.array(out)


def nu_einfach(gN, a0=A0):
    y = gN / a0
    return 0.5 + np.sqrt(0.25 + 1 / y)


def nu_standard(gN, a0=A0):
    y = gN / a0
    return np.sqrt(0.5 + np.sqrt(0.25 + 1 / y**2))


if __name__ == "__main__":
    for ud, sg in ((0.5, 0.0), (0.4, 0.0), (0.6, 0.0), (0.5, 8.0)):
        g, gN = punkte(ud, sg)
        print("=== Υ_d = %.1f, Gasdruck σ = %.0f km/s: %d Punkte" % (ud, sg, g.size))
        # 1) Linearisierung, nur wo g > g_N deutlich (ν > 1,05), gebinnt
        y = np.where(g > 1.05 * gN, np.log(g / np.clip(g - gN, 1e-40, None)), np.nan)
        B = bins(gN, y)
        s = np.sqrt(10**B[:, 0])
        k, cov = np.polyfit(s, B[:, 1], 1, w=1 / B[:, 3], cov=True)
        k0 = np.sum(s * B[:, 1] / B[:, 3]**2) / np.sum(s**2 / B[:, 3]**2)
        q2 = np.polyfit(s, B[:, 1], 2, w=1 / B[:, 3])
        print("  Bose-Linearisierung y = ln(g/(g − g_N)) gegen s = √(g_N/a₀): Gerade y = %.3f·s %+.3f (Achsenabschnitt ± %.3f); durch 0: Steigung %.3f ⇒ a₀ = %.3f·10⁻¹⁰; Krümmung (quadratisch) %.3f"
              % (k[0], k[1], math.sqrt(cov[1, 1]), k0, A0 / k0**2 * 1e10, q2[0]))
        for lab, f in (("Bose (Regel)", lambda x: nu(x)), ("einfach", nu_einfach), ("Standard", nu_standard)):
            yy = np.log(f(10**B[:, 0] * A0) / (f(10**B[:, 0] * A0) - 1))
            chi2 = np.sum(((B[:, 1] - yy) / B[:, 3])**2)
            print("     Form %-13s: χ² = %8.1f für %d Bins" % (lab, chi2, len(B)))
        # 2) Wien-Stelle der Stufe 2
        B2 = bins(gN, (g - gN) / A0, -1.2, 1.4, 13)
        j = np.argmax(B2[:, 1])
        p = np.polyfit(B2[max(j - 2, 0):j + 3, 0], B2[max(j - 2, 0):j + 3, 1], 2)
        xm = -p[1] / (2 * p[0])
        print("  Zusatzbeschleunigung (g − g_N)/a₀, Median je Bin: " + "  ".join("%.2f:%.3f" % (10**a, b) for a, b, *_ in B2))
        print("  Maximum (Parabel um den größten Bin): bei g_N = %.2f·a₀, Höhe %.3f·a₀   (Regel: 2,54·a₀, 0,648·a₀)" % (10**xm, np.polyval(p, xm)))
