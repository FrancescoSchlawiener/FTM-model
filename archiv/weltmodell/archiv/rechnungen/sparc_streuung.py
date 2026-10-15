# -*- coding: utf-8 -*-
"""
Woher kommt die Streuung des Rests je Galaxie? Monte-Carlo über die Fehler des SPARC-Katalogs:
Abstand (e_D), Neigung (e_Inc), Masse-zu-Licht (log-normal, 0,1 dex je Galaxie, Scheibe und Bulge gemeinsam).
Für jede Galaxie werden die gemessenen Größen mit den gezogenen Fehlern neu skaliert und der Rest neu gerechnet;
die Streuung, die allein die Fehler erzeugen, wird mit der gemessenen verglichen (Gasdruck σ = 8 km/s eingerechnet).
Skalierung: R ∝ D, V_obs ∝ 1/sin i, V_disk², V_bul² ∝ Υ·D (Leuchtkraft ∝ D², Radius ∝ D), V_gas² ∝ D.
"""
import math
import numpy as np
from ftm_kern import A0, nu, pc
import sparc_flaechendichte as sf, sparc_universal as su

kpc = 1e3 * pc
rng = np.random.default_rng(29092026)


def rest(R, Vo, Vg, Vd, Vb, Rd, ud=0.5, sig=8.0):
    Vb2 = Vg * np.abs(Vg) + ud * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        return np.nan
    gN = Vb2[m] / R[m]
    go = (Vo[m]**2 + (sig * 1e3)**2 * (R[m] / kpc) / Rd) / R[m]
    return np.median(np.log10(go / (nu(gN, A0) * gN)))


def main(n_mc=200, dups=0.1):
    t = su.tab_voll(); gem, sim = [], []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        r0 = rest(R, Vo, Vg, Vd, Vb, d['Rd'])
        if not np.isfinite(r0):
            continue
        gem.append(r0)
        rs = []
        for _ in range(n_mc):
            fD = max(1 + rng.normal() * d['eD'] / d['D'], 0.2)
            inc = min(max(d['inc'] + rng.normal() * d['eInc'], 20.0), 90.0)
            fi = math.sin(math.radians(d['inc'])) / math.sin(math.radians(inc))
            fU = 10**(rng.normal() * dups)
            rs.append(rest(R * fD, Vo * fi, Vg * math.sqrt(fD), Vd * math.sqrt(fU * fD), Vb * math.sqrt(fU * fD), d['Rd'] * fD))
        sim.append(np.nanstd(rs))
    gem = np.array(gem); sim = np.array(sim)
    return gem, sim


if __name__ == "__main__":
    for dups in (0.0, 0.1, 0.15):
        g, s = main(dups=dups)
        erw = math.sqrt(np.mean(s**2))
        print("Υ-Streuung %.2f dex: gemessene Streuung des Rests %.3f dex, aus den Fehlern erwartet %.3f dex, Rest innen %.3f dex (n = %d)"
              % (dups, g.std(), erw, math.sqrt(max(g.var() - erw**2, 0)), len(g)))


def a0_je_galaxie(R, Vo, Vg, Vd, Vb, Rd, ud=0.5, sig=8.0):
    from scipy.optimize import brentq
    Vb2 = Vg * np.abs(Vg) + ud * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        return np.nan
    gN = Vb2[m] / R[m]; go = (Vo[m]**2 + (sig * 1e3)**2 * (R[m] / kpc) / Rd) / R[m]
    f = lambda la: np.median(np.log10(go / (nu(gN, 10**la) * gN)))
    try:
        return brentq(f, -13, -7)
    except ValueError:
        return np.nan


def schwelle(n_mc=100, dups=0.1):
    t = su.tab_voll(); gem, sim = [], []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        l0 = a0_je_galaxie(R, Vo, Vg, Vd, Vb, d['Rd'])
        if not np.isfinite(l0):
            continue
        gem.append(l0); ls = []
        for _ in range(n_mc):
            fD = max(1 + rng.normal() * d['eD'] / d['D'], 0.2)
            inc = min(max(d['inc'] + rng.normal() * d['eInc'], 20.0), 90.0)
            fi = math.sin(math.radians(d['inc'])) / math.sin(math.radians(inc)); fU = 10**(rng.normal() * dups)
            ls.append(a0_je_galaxie(R * fD, Vo * fi, Vg * math.sqrt(fD), Vd * math.sqrt(fU * fD), Vb * math.sqrt(fU * fD), d['Rd'] * fD))
        sim.append(np.nanstd(ls))
    return np.array(gem), np.array(sim)


if __name__ == "__main__":
    for dups in (0.0, 0.1):
        g, s = schwelle(dups=dups)
        mad = 1.4826 * np.median(np.abs(g - np.median(g)))
        print("Schwelle je Galaxie (Υ-Streuung %.2f dex): Median a₀ = %.2f·10⁻¹⁰, Streuung %.3f dex (robust %.3f), aus den Fehlern erwartet %.3f dex (Median der Einzelfehler %.3f), n = %d"
              % (dups, 10**np.median(g) * 1e10, g.std(), mad, math.sqrt(np.mean(s**2)), np.median(s), len(g)))


def anstieg(n_boot=2000, dups=0.1):
    """Steigung d log a₀ / d log Σ₀ (Σ₀ = 0,5·SB₀) der direkt gelesenen Schwelle; Fehler per Bootstrap über Galaxien, robust (Theil-Sen)."""
    from scipy.stats import theilslopes
    t = su.tab_voll(); x, y = [], []
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        l0 = a0_je_galaxie(R, Vo, Vg, Vd, Vb, d['Rd'])
        if np.isfinite(l0):
            x.append(math.log10(0.5 * d['SB0'])); y.append(l0)
    x, y = np.array(x), np.array(y)
    k = theilslopes(y, x)[0]
    bs = [theilslopes(y[i], x[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(n_boot))]
    return k, np.std(bs), len(x)


if __name__ == "__main__":
    k, ek, n = anstieg()
    print("Anstieg der direkt gelesenen Schwelle mit der Flächendichte (Theil-Sen, Gasdruck 8 km/s): %+.3f ± %.3f je Dekade (n = %d)" % (k, ek, n))
