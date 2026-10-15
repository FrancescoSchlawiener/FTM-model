# -*- coding: utf-8 -*-
"""
Zentrale Flächendichte (Toomre-Integral) Galaxie für Galaxie: gemessene Rotationskurve gegen die Kurve der örtlichen Regel
mit den eigenen Baryonen derselben Galaxie (SPARC, Lelli, McGaugh, Schombert 2016).
Σ_dyn(0) = (1/2πG)∫₀^∞ V²/R² dR, innen V ∝ R bis zum ersten Punkt, außen V = V_letzter; gleiche Behandlung für beide Kurven.
Regel: g_N = V_bar²/R mit V_bar² = V_gas|V_gas| + Υ_d V_disk² + Υ_b V_bul² (Υ_d = 0,5, Υ_b = 0,7), a = ν(g_N)·g_N (radiale Form, f = 0).
Auswahl: Q = 1, 2; Neigung ≥ 30°.
"""
import math, os
import numpy as np
from ftm_kern import G, M_sun, pc, A0, nu, nu_aussenfeld

D = os.path.join(os.path.dirname(__file__), '..', 'daten', 'sparc')
kpc = 1e3 * pc
MSPC = M_sun / pc**2
UPS_D, UPS_B = 0.5, 0.7


def tabelle():
    t = {}
    for l in open(os.path.join(D, 'SPARC_Lelli2016c.mrt'), encoding='latin-1'):
        f = l.split()
        if len(f) < 18:
            continue
        try:
            t[f[0]] = dict(inc=float(f[5]), Rd=float(f[11]), SB0=float(f[12]), Q=int(f[17]))
        except ValueError:
            pass
    return t


def kurve(name):
    a = np.loadtxt(os.path.join(D, 'Rotmod_LTG', name + '_rotmod.dat'), comments='#')
    return a[:, 0] * kpc, a[:, 1] * 1e3, a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3


def toomre(R, V):
    f = (V / R) ** 2
    innen = V[0]**2 / R[0]
    mitte = np.trapezoid(f, R)
    aussen = V[-1]**2 / R[-1]
    return (innen + mitte + aussen) / (2 * math.pi * G) / MSPC


def main(ups_d=UPS_D, a0=A0, ge=0.0):
    t = tabelle()
    rows = []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        try:
            R, Vo, Vg, Vd, Vb = kurve(nm)
        except OSError:
            continue
        Vb2 = Vg * np.abs(Vg) + ups_d * Vd**2 + UPS_B * Vb**2
        gN = np.clip(Vb2, 1e-6, None) / R
        nn = nu(gN, a0) if ge == 0 else np.array([nu_aussenfeld(x, ge * a0, 'senkrecht', a0) for x in gN])
        Vm = np.sqrt(nn * gN * R)
        Sb = toomre(R, np.sqrt(np.clip(Vb2, 0, None)))
        rows.append((nm, ups_d * d['SB0'], toomre(R, Vo), toomre(R, Vm), Sb))
    return rows


if __name__ == "__main__":
    rows = main()
    S0 = np.array([r[1] for r in rows]); So = np.array([r[2] for r in rows]); Sm = np.array([r[3] for r in rows]); Sb = np.array([r[4] for r in rows])
    ok = (S0 > 0) & (So > 0) & (Sm > 0)
    S0, So, Sm, Sb = S0[ok], So[ok], Sm[ok], Sb[ok]
    res = np.log10(So / Sm)
    print("SPARC, Q ≤ 2, Neigung ≥ 30°: %d Galaxien, a₀ = %.2e, Υ_d = %.1f, Υ_b = %.1f" % (len(res), A0, UPS_D, UPS_B))
    print("log Σ_dyn(gemessen)/Σ_dyn(Regel): Mittel %+.3f, Median %+.3f, Streuung %.3f dex, Fehler des Mittels %.3f" %
          (res.mean(), np.median(res), res.std(ddof=1), res.std(ddof=1) / math.sqrt(len(res))))
    for lo, hi in ((0, 30), (30, 300), (300, 1e5)):
        m = (S0 >= lo) & (S0 < hi)
        print("   Σ₀ = Υ_d·SB₀ in [%g, %g) M☉/pc²: n = %3d, Rest %+.3f ± %.3f dex" % (lo, hi, m.sum(), res[m].mean(), res[m].std(ddof=1) / math.sqrt(m.sum())))
    for lab, y in (("gemessen", So), ("Regel", Sm), ("Newton (Baryonen)", Sb)):
        k = np.polyfit(np.log10(S0), np.log10(y), 1)
        lowm = S0 < 300
        kl = np.polyfit(np.log10(S0[lowm]), np.log10(y[lowm]), 1)
        print("   Steigung log Σ_dyn gegen log Σ₀, %-18s: alle %.3f, Σ₀ < 300: %.3f" % (lab, k[0], kl[0]))


def trend(rows):
    S0 = np.array([r[1] for r in rows]); So = np.array([r[2] for r in rows]); Sm = np.array([r[3] for r in rows])
    ok = (S0 > 0) & (So > 0) & (Sm > 0)
    x, y = np.log10(S0[ok]), np.log10(So[ok] / Sm[ok])
    A = np.vstack([x - x.mean(), np.ones_like(x)]).T
    cf, rs, *_ = np.linalg.lstsq(A, y, rcond=None)
    s2 = rs[0] / (len(x) - 2)
    return y.mean(), y.std(ddof=1) / math.sqrt(len(y)), cf[0], math.sqrt(s2 / np.sum((x - x.mean())**2)), y[x < math.log10(30)].mean()


if __name__ == "__main__":
    print("\nStörgrößen: Mittel des Rests, Trend d(Rest)/d log Σ₀ (Steigung gemessen − Regel), Rest bei Σ₀ < 30")
    print("  Υ_d | a₀ [10⁻¹⁰] | g_e/a₀ (senkrecht) | Rest [dex]        | Trend           | Rest Σ₀ < 30")
    for ud, a0x, ge in ((0.5, 1.1, 0), (0.4, 1.1, 0), (0.6, 1.1, 0), (0.5, 1.04, 0), (0.5, 1.19, 0), (0.5, 1.1, 0.01), (0.5, 1.1, 0.03)):
        m, em, k, ek, lo = trend(main(ud, a0x * 1e-10, ge))
        print("  %.1f | %10.2f | %18.2f | %+.3f ± %.3f    | %+.3f ± %.3f  | %+.3f" % (ud, a0x, ge, m, em, k, ek, lo))
