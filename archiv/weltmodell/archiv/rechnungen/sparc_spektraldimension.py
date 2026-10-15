# -*- coding: utf-8 -*-
"""
Spektraldimension der Rückkehr an den tiefsten SPARC-Punkten (Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Υ_b = 0,7, Gasdruck σ = 8 km/s).
Die Rückkehr ist im Modell eindimensional und driftfrei: 1 − q ∝ Γ^{1−d_s/2}, a ∝ a_N^{d_s/2}, d_s = 2·d ln a/d ln a_N.
Die volle Bose-Form hat die Steigung 1 − χ/2 (χ = s/(e^s − 1)), also d_s,eff = 2 − χ; d_s = 1 erst für s → 0.

1) Anzahl tiefer Punkte/Galaxien (s < 0,3 und s < 0,15).
2) d_s aus (a) Steigungen zwischen Bin-Medianen, (b) Steigungen je Galaxie (frei von Abstand und Neigung), (c) Potenzgesetz
   über s < 0,3 und s < 0,15. Jeweils derselbe Schätzer auf die Regel ν(a_N)·a_N an denselben Punkten angewandt (Modellwert).
   Fehler: Bootstrap über Galaxien, in jeder Ziehung Abstand, Neigung, Υ (0,1 dex) neu gezogen wie in sparc_streuung.
3) Boden: a = ν(a_N)·a_N + B, B frei; obere Grenze 95 %; Drift Pe = B/(2a₀), Querrichtung b/ϖ_R = √(2/π)·a₀/B.
4) Andere Dimensionen: Gitter-Rückkehr ν_d = G_d(1 − ε), ε = κ·a_N/a₀ (G_d(z) = ∫₀^∞ e^{−t} I₀(zt/d)^d dt),
   d = 1, 2 (Logarithmus), 3 (ν ≤ G_3(1) = 1,516, Watson); χ² der tiefen Bins mit Bootstrap-Kovarianz.
Rechnet nur mit ftm_kern (A0, nu, chi).
"""
import math, os
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.integrate import quad
from scipy.special import i0e
from ftm_kern import A0, nu, chi, pc
import sparc_flaechendichte as sf, sparc_universal as su

kpc = 1e3 * pc
rng = np.random.default_rng(29092026)
SIG = 8.0
S_EDGES = np.array([0.1, 0.15, 0.2, 0.3, 0.45, 0.7, 1.0, 1.5, 2.5, 5.0])   # Bins in s = √(a_N/a₀)


def galaxien():
    t = su.tab_voll(); out = []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat'), comments='#')
        out.append(dict(nm=nm, R=a[:, 0] * kpc, Vo=a[:, 1] * 1e3, eV=a[:, 2] * 1e3, Vg=a[:, 3] * 1e3, Vd=a[:, 4] * 1e3, Vb=a[:, 5] * 1e3, **d))
    return out


def punkte(g, fD=1.0, fi=1.0, fU=1.0, ud=0.5, sig=SIG):
    R = g['R'] * fD; Vo = g['Vo'] * fi
    Vb2 = fD * (g['Vg'] * np.abs(g['Vg']) + fU * (ud * g['Vd']**2 + 0.7 * g['Vb']**2))
    m = (Vb2 > 0) & (Vo > 0)
    gN = Vb2[m] / R[m]
    gc = (Vo[m]**2 + (sig * 1e3)**2 * (R[m] / kpc) / (g['Rd'] * fD)) / R[m]
    return gN, gc


def ziehung(gal, boot=True, mc=True, sig=SIG):
    idx = rng.integers(0, len(gal), len(gal)) if boot else np.arange(len(gal))
    N, C, I = [], [], []
    for k, j in enumerate(idx):
        g = gal[j]
        if mc:
            fD = max(1 + rng.normal() * g['eD'] / g['D'], 0.2)
            inc = min(max(g['inc'] + rng.normal() * g['eInc'], 20.0), 90.0)
            fi = math.sin(math.radians(g['inc'])) / math.sin(math.radians(inc)); fU = 10**(rng.normal() * 0.1)
        else:
            fD = fi = fU = 1.0
        gN, gc = punkte(g, fD, fi, fU, sig=sig)
        N.append(gN); C.append(gc); I.append(np.full(gN.size, k))
    return np.concatenate(N), np.concatenate(C), np.concatenate(I)


def modell(gN, a0=A0):
    return nu(gN, a0) * gN


# ---------------------------------------------------------------- Schätzer
def bin_steigungen(gN, gc):
    """Mediane je s-Bin; Steigung zwischen benachbarten Bins für Messung und Regel (gleiche Punkte)."""
    s = np.sqrt(gN / A0); x = np.log(gN); y = np.log(gc); ym = np.log(modell(gN))
    X, Y, YM, S = [], [], [], []
    for a, b in zip(S_EDGES[:-1], S_EDGES[1:]):
        m = (s >= a) & (s < b)
        X.append(np.median(x[m]) if m.sum() >= 5 else np.nan); Y.append(np.median(y[m]) if m.sum() >= 5 else np.nan)
        YM.append(np.median(ym[m]) if m.sum() >= 5 else np.nan); S.append(np.sqrt(np.exp(X[-1]) / A0))
    X, Y, YM, S = map(np.array, (X, Y, YM, S))
    k = np.diff(Y) / np.diff(X); km = np.diff(YM) / np.diff(X); sm = np.sqrt(S[:-1] * S[1:])
    return sm, k, km, Y - YM


def potenz(gN, gc, smax):
    s = np.sqrt(gN / A0); m = s < smax
    x = np.log(gN[m])
    r = least_squares(lambda p: (p[0] + p[1] * x - np.log(gc[m])), [0, 0.5], loss='soft_l1', f_scale=0.1)
    rm = np.polyfit(x, np.log(modell(gN[m])), 1)
    return r.x[1], rm[0], m.sum()


def je_galaxie(gN, gc, I, smax=0.3, nmin=4, dex=0.2):
    s = np.sqrt(gN / A0); k, km, w = [], [], []
    for j in np.unique(I):
        m = (I == j) & (s < smax)
        if m.sum() < nmin:
            continue
        x = np.log(gN[m])
        if np.ptp(x) < dex * math.log(10):
            continue
        k.append(np.polyfit(x, np.log(gc[m]), 1)[0]); km.append(np.polyfit(x, np.log(modell(gN[m])), 1)[0])
    k, km = np.array(k), np.array(km)
    return np.median(k), np.median(km), np.median(k - km), len(k)


def boden(gN, gc, smax=0.3, a0_frei=False, smax_a0=None):
    """a = ν(a_N; a₀)·a_N + B, robust (soft_l1) in ln a; B/a₀ frei (auch < 0 erlaubt, Grenze aus der Verteilung)."""
    s = np.sqrt(gN / A0); m = s < (smax_a0 if a0_frei else smax)
    x, y = gN[m], np.log(gc[m])

    def f(p):
        a0 = A0 * 10**p[1] if a0_frei else A0
        z = modell(x, a0) + p[0] * A0
        return np.log(np.clip(z, 1e-40, None)) - y
    r = least_squares(f, [0.0, 0.0] if a0_frei else [0.0], loss='soft_l1', f_scale=0.1)
    return r.x[0], (10**r.x[1] if a0_frei else 1.0)


# ---------------------------------------------------------------- Gitter-Rückkehr in d Dimensionen
def G_d(eps, d):
    z = 1 - eps
    f = lambda t: np.exp(-t * eps) * i0e(z * t / d)**d   # = e^{−t} I₀(zt/d)^d
    val, _ = quad(lambda u: f(math.exp(u)) * math.exp(u), -30, math.log(60.0 / eps), limit=1000)   # t = e^u
    return val


_EPS = np.logspace(-7, -0.001, 160)
_TAB = {d: np.array([G_d(e, d) for e in _EPS]) for d in (1, 2, 3)}


def nu_d(gN, kappa, d):
    e = np.clip(kappa * gN / A0, _EPS[0], _EPS[-1])
    return np.interp(np.log(e), np.log(_EPS), _TAB[d])


def nu_familie(gN, ds, a0=A0):
    """1 − q = 1 − e^{−s^{2−d_s}}: Bose-Form für d_s = 1, im Tiefen a ∝ a_N^{d_s/2}."""
    return 1 / (-np.expm1(-(gN / a0)**(1 - ds / 2)))


def familie(gN, gc, smax=99.0):
    m = np.sqrt(gN / A0) < smax; x, y = gN[m], np.log(gc[m])
    r = least_squares(lambda p: np.log(nu_familie(x, p[0], A0 * 10**p[1]) * x) - y, [1.0, 0.0], bounds=([0.2, -1], [1.9, 1]), loss='soft_l1', f_scale=0.1)
    return r.x[0], 10**r.x[1]


def tiefe_bins(gN, gc, edges=(0.1, 0.13, 0.16, 0.2, 0.25, 0.3)):
    s = np.sqrt(gN / A0); out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (s >= a) & (s < b)
        out.append(np.median(np.log(gc[m])) if m.sum() >= 5 else np.nan)
    return np.array(out)


def formen_chi2(gN0, gc0, reps):
    edges = (0.1, 0.13, 0.16, 0.2, 0.25, 0.3)
    Y0 = tiefe_bins(gN0, gc0, edges)
    Yr = np.array([tiefe_bins(gN, gc, edges) for gN, gc, _ in reps])
    Yr = Yr[np.all(np.isfinite(Yr), 1)]
    Cinv = np.linalg.inv(np.cov(Yr.T))
    s = np.sqrt(gN0 / A0); ms = [(s >= a) & (s < b) for a, b in zip(edges[:-1], edges[1:])]

    def chi2(fn):
        mod = np.array([np.median(np.log(fn(gN0[m]))) for m in ms]); r = Y0 - mod
        return r @ Cinv @ r
    res = {}
    res['Bose, a₀ frei'] = minimize_scalar(lambda la: chi2(lambda x: modell(x, A0 * 10**la)), bounds=(-1.5, 1.5), method='bounded')
    res['Bose, a₀ = 1,1'] = type('R', (), dict(fun=chi2(modell), x=0.0))()
    for d in (1, 2, 3):
        res['Gitter d = %d, κ frei' % d] = minimize_scalar(lambda lk: chi2(lambda x: nu_d(x, 10**lk, d) * x), bounds=(-6.5, 4), method='bounded')
    best = min(((chi2(lambda x: A * nu_d(x, 10**lk, 2) * x), A, lk) for A in np.linspace(0.3, 6, 115) for lk in np.linspace(-6.5, 4, 106)))
    res['Gitter d = 2, κ, Amplitude frei'] = type('R', (), dict(fun=best[0], x=(best[1], best[2])))()
    # Potenzgesetz frei (2 Parameter)
    best = min(((chi2(lambda x: np.exp(c) * (x / A0)**p * A0), p, c) for p in np.linspace(0.3, 1.2, 91) for c in np.linspace(-1, 2, 151)))
    res['Potenz frei'] = type('R', (), dict(fun=best[0], x=(best[1], best[2])))()
    # Bose mit Boden (a₀ = 1,1), B frei
    res['Bose + Boden B'] = minimize_scalar(lambda B: chi2(lambda x: modell(x) + B * A0), bounds=(-0.05, 0.3), method='bounded')
    return res, len(Y0)


if __name__ == "__main__":
    gal = galaxien()
    gN0, gc0, I0 = ziehung(gal, boot=False, mc=False)
    s0 = np.sqrt(gN0 / A0)
    print("=== 1) Tiefe Punkte (Q ≤ 2, Neigung ≥ 30°, Υ_d 0,5, Gasdruck σ = %.0f km/s): %d Punkte, %d Galaxien gesamt" % (SIG, gN0.size, len(np.unique(I0))))
    for c in (0.3, 0.2, 0.15, 0.1):
        m = s0 < c
        print("   s < %.2f (a_N < %.4f·a₀): %4d Punkte, %3d Galaxien" % (c, c * c, m.sum(), len(np.unique(I0[m]))))
    print("   kleinstes s = %.3f; 1 %%-Quantil %.3f" % (s0.min(), np.percentile(s0, 1)))

    NB = 400
    reps = [ziehung(gal) for _ in range(NB)]
    reps_boot = [ziehung(gal, mc=False) for _ in range(150)]

    print("\n=== 2a) Steigung zwischen Bin-Medianen (Bootstrap über Galaxien + Monte-Carlo Abstand/Neigung/Υ, %d Ziehungen)" % NB)
    sm0, k0, km0, off0 = bin_steigungen(gN0, gc0)
    K = np.array([bin_steigungen(gN, gc)[1] for gN, gc, _ in reps]); KM = np.array([bin_steigungen(gN, gc)[2] for gN, gc, _ in reps])
    D = K - KM
    print("   s(Mitte)  Steigung gem.   σ      Regel(gleiche Punkte)  1−χ/2   gem.−Regel (σ)     d_s gem.     d_s Regel   (Regel−½)/σ")
    for i in range(len(k0)):
        sd = np.nanstd(K[:, i]); dd = np.nanstd(D[:, i])
        print("   %6.3f    %6.3f      %5.3f        %6.3f            %6.3f   %+6.3f (%.3f)   %5.2f±%.2f   %5.2f        %4.1f"
              % (sm0[i], k0[i], sd, km0[i], 1 - float(chi(sm0[i]**2 * A0)) / 2, k0[i] - km0[i], dd, 2 * k0[i], 2 * sd, 2 * km0[i], (km0[i] - 0.5) / sd))
    print("   Versatz ln(g_gem/g_Regel) je Bin (Median): " + "  ".join("%.2f:%+.3f" % (a, b) for a, b in zip(np.sqrt(S_EDGES[:-1] * S_EDGES[1:]), off0)))

    print("\n=== 2b) Steigung je Galaxie über ihre Punkte mit s < s_max (≥ 4 Punkte, ≥ 0,2 dex in a_N); frei von Abstand und Neigung")
    for smax in (0.3, 0.45):
        k, km, dk, n = je_galaxie(gN0, gc0, I0, smax)
        B = np.array([je_galaxie(gN, gc, I, smax)[2] for gN, gc, I in reps])
        print("   s < %.2f: n = %d Galaxien; Median Steigung gem. %.3f (d_s %.2f), Regel %.3f (d_s %.2f); Median(gem.−Regel) %+.3f ± %.3f"
              % (smax, n, k, 2 * k, km, 2 * km, dk, np.nanstd(B)))

    print("\n=== 2c) Potenzgesetz a ∝ a_N^{d_s/2} über alle Punkte unter s_max (robust)")
    for smax in (0.3, 0.2, 0.15):
        p, pm, n = potenz(gN0, gc0, smax)
        P = np.array([potenz(gN, gc, smax)[:2] for gN, gc, _ in reps])
        Pb = np.array([potenz(gN, gc, smax)[:2] for gN, gc, _ in reps_boot])
        print("   s < %.2f (%4d Punkte): d_s gem. = %.3f ± %.3f (nur Bootstrap ± %.3f); Regel an denselben Punkten %.3f; gem.−Regel %+.3f ± %.3f; (Regel − 1)/σ = %.1f"
              % (smax, n, 2 * p, 2 * P[:, 0].std(), 2 * Pb[:, 0].std(), 2 * pm, 2 * (p - pm), 2 * (P[:, 0] - P[:, 1]).std(), (2 * pm - 1) / (2 * P[:, 0].std())))

    print("\n=== 2d) Familie ν = 1/(1 − e^{−(a_N/a₀)^{1−d_s/2}}) (d_s = 1: Bose-Form), d_s und a₀ frei, robust")
    for smax in (99.0, 1.0):
        ds, f = familie(gN0, gc0, smax)
        F = np.array([familie(gN, gc, smax) for gN, gc, _ in reps])
        print("   %-12s d_s = %.3f ± %.3f (a₀-Faktor %.3f ± %.3f);  (d_s − 1)/σ = %+.1f" % ("alle Punkte" if smax > 50 else "s < 1", ds, F[:, 0].std(), f, F[:, 1].std(), (ds - 1) / F[:, 0].std()))

    print("\n=== 3) Boden a = ν(a_N)·a_N + B")
    for lab, smax, fr, sa in (("a₀ = 1,1 fest, s < 0,3", 0.3, False, None), ("a₀ = 1,1 fest, s < 0,45", 0.45, False, None),
                              ("a₀ frei, s < 1", 0.3, True, 1.0), ("a₀ frei, alle Punkte", 0.3, True, 99.0)):
        b0, a0f = boden(gN0, gc0, smax, fr, sa)
        Bs = np.array([boden(gN, gc, smax, fr, sa)[0] for gN, gc, _ in reps])
        up = max(np.percentile(Bs, 95), 0.0)
        print("   %-24s B/a₀ = %+.4f ± %.4f (a₀-Faktor %.3f); 95 %%: B < %.4f·a₀ ⇒ Pe < %.4f, b > %.0f·ϖ_R"
              % (lab, b0, Bs.std(), a0f, up, up / 2, math.sqrt(2 / math.pi) / up if up > 0 else np.inf))
    Bs = np.array([boden(gN, gc, 0.3)[0] for gN, gc, _ in (ziehung(gal, sig=0.0) for _ in range(200))])
    print("   ohne Gasdruck (σ = 0), a₀ fest, s < 0,3: B/a₀ = %+.4f ± %.4f; 95 %%: B < %.4f·a₀" % (boden(*ziehung(gal, False, False, sig=0.0)[:2], 0.3)[0], Bs.std(), max(np.percentile(Bs, 95), 0)))
    print("   Gasdruck als Stellgröße (a₀ = 1,1 fest, s < 0,3, nur Nennwerte): " + "  ".join(
        "σ = %2.0f km/s: B/a₀ = %+.4f" % (sg, boden(*ziehung(gal, False, False, sig=sg)[:2], 0.3)[0]) for sg in (0.0, 5.0, 8.0, 10.0)))
    rd = np.array([g['Rd'] for g in gal])
    tief = np.unique(I0[s0 < 0.3])
    print("   Gasdruckterm σ²/R_d der Galaxien mit tiefen Punkten (σ = 8 km/s): Median %.4f·a₀ (Spanne 16–84 %%: %.4f–%.4f)"
          % tuple(np.percentile((8e3)**2 / (rd[tief] * kpc) / A0, [50, 16, 84])))

    print("\n=== 4) Rückkehr in d Dimensionen gegen die tiefen Bins (s = 0,1–0,3; Bootstrap-Kovarianz)")
    print("   G_1(1 − ε)·√(2ε) → 1: %.4f (ε = 10⁻⁶);  G_3(1) = %.4f (Watson 1,5164);  G_2(1 − 10⁻⁴) = %.3f, (1/π)·ln(8/ε) = %.3f"
          % (G_d(1e-6, 1) * math.sqrt(2e-6), G_d(1e-9, 3), G_d(1e-4, 2), math.log(8e4) / math.pi))
    res, nb = formen_chi2(gN0, gc0, reps)
    ref = res['Bose, a₀ frei'].fun
    for k, r in res.items():
        print("   %-22s χ² = %9.1f (%d Bins)  Δχ² gegen Bose(a₀ frei) = %+9.1f   Parameter %s" % (k, r.fun, nb, r.fun - ref, np.round(r.x, 3)))
