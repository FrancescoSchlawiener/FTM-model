# -*- coding: utf-8 -*-
"""
Universalität der örtlichen Regel an SPARC: Rest je Galaxie (Median von log g_gem/g_Regel) gegen Eigenschaften der Galaxie.
Ist die Regel universell, hängt der Rest nur an Messgrößen (Abstand, Neigung, Υ), nicht an physikalischen Eigenschaften.
Spearman-Korrelation mit Irrtumswahrscheinlichkeit; dazu Rest innen (R < 2R_d) gegen außen (R > 2R_d).
"""
import math, os
import numpy as np
from scipy.stats import spearmanr
from ftm_kern import A0, nu, pc
import sparc_flaechendichte as sf

kpc = 1e3 * pc


def tab_voll():
    t = {}
    for l in open(os.path.join(sf.D, 'SPARC_Lelli2016c.mrt'), encoding='latin-1'):
        f = l.split()
        if len(f) < 18:
            continue
        try:
            t[f[0]] = dict(T=int(f[1]), D=float(f[2]), eD=float(f[3]), fD=int(f[4]), inc=float(f[5]), eInc=float(f[6]), L=float(f[7]), Reff=float(f[9]),
                           Rd=float(f[11]), SB0=float(f[12]), MHI=float(f[13]), Vf=float(f[15]), Q=int(f[17]))
        except ValueError:
            pass
    return t


def reste(ups_d=0.5, a0=A0, sig=0.0, ge=0.0):
    t = tab_voll(); out = {}
    for nm, d in t.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        Vb2 = Vg * np.abs(Vg) + ups_d * Vd**2 + 0.7 * Vb**2
        m = (Vb2 > 0) & (Vo > 0)
        if m.sum() < 5:
            continue
        gN = Vb2[m] / R[m]; go = (Vo[m]**2 + (sig * 1e3)**2 * (R[m] / kpc) / d['Rd']) / R[m]   # Gasdruck: V_c² = V² + σ²R/R_d
        nn_ = nu(gN, a0) if ge == 0 else np.array([float(nu(math.hypot(x, ge * a0), a0)) for x in gN])   # äußeres Feld senkrecht
        r = np.log10(go / (nn_ * gN))
        Rk = R[m] / kpc
        inn, aus = r[Rk < 2 * d['Rd']], r[Rk > 2 * d['Rd']]
        Mb = ups_d * d['L'] + 1.33 * d['MHI']
        out[nm] = dict(rest=np.median(r), innen=np.median(inn) if inn.size else np.nan, aussen=np.median(aus) if aus.size else np.nan,
                       lVf=math.log10(d['Vf']) if d['Vf'] > 0 else np.nan, lRd=math.log10(d['Rd']), lSB=math.log10(d['SB0']),
                       fgas=1.33 * d['MHI'] / Mb, inc=d['inc'], lD=math.log10(d['D']), T=d['T'], eDrel=d['eD'] / d['D'],
                       lMb=math.log10(Mb), lgmin=math.log10(gN.min() / a0))
    return out


if __name__ == "__main__":
    o = reste()
    names = list(o)
    rest = np.array([o[n]['rest'] for n in names])
    print("SPARC Q ≤ 2, Neigung ≥ 30°: %d Galaxien; Rest je Galaxie Median %.3f dex, Streuung %.3f dex" % (len(rest), np.median(rest), rest.std()))
    print("Spearman ρ des Rests mit (p = Irrtumswahrscheinlichkeit; |ρ| mit p < 0,003 auffällig):")
    for k, lab in (("lMb", "log Baryonenmasse"), ("lVf", "log V_flach"), ("lRd", "log R_d"), ("lSB", "log SB₀ (Flächenhelligkeit)"),
                   ("fgas", "Gasanteil"), ("lgmin", "log kleinstes g_N/a₀"), ("T", "Hubble-Typ"),
                   ("inc", "Neigung [Messgröße]"), ("lD", "log Abstand [Messgröße]"), ("eDrel", "rel. Abstandsfehler [Messgröße]")):
        x = np.array([o[n][k] for n in names]); m = np.isfinite(x)
        rho, p = spearmanr(x[m], rest[m])
        print("   %-34s ρ = %+.3f  p = %.1e  %s" % (lab, rho, p, "AUFFÄLLIG" if p < 0.003 else ""))
    inn = np.array([o[n]['innen'] for n in names]); aus = np.array([o[n]['aussen'] for n in names]); m = np.isfinite(inn) & np.isfinite(aus)
    print("Rest innen (R < 2R_d) − außen (R > 2R_d): Median %+.3f dex, Mittel %+.3f ± %.3f (n = %d)" % (np.median(inn[m] - aus[m]), (inn[m] - aus[m]).mean(), (inn[m] - aus[m]).std() / math.sqrt(m.sum()), m.sum()))
    for ud in (0.4, 0.6):
        oo = reste(ud); rr = np.array([oo[n]['rest'] for n in oo]); x = np.array([oo[n]['lSB'] for n in oo])
        print("   Υ_d = %.1f: Median Rest %+.3f, ρ(Rest, log SB₀) = %+.3f (p = %.1e)" % (ud, np.median(rr), *spearmanr(x, rr)))
    for lab, sel in (("gasdominiert (f_gas > 0,5)", lambda d: d['fgas'] > 0.5), ("sterndominiert (f_gas < 0,2)", lambda d: d['fgas'] < 0.2)):
        nn = [n for n in names if sel(o[n])]
        x = np.array([o[n]['lSB'] for n in nn]); r = np.array([o[n]['rest'] for n in nn])
        print("   %-30s n = %3d: ρ(Rest, log SB₀) = %+.3f (p = %.1e), Median Rest %+.3f" % (lab, len(nn), *spearmanr(x, r), np.median(r)))
    print("Gasdruck (V_c² = V_obs² + σ²·R/R_d): Streuung des Rests und Zusammenhang mit V_flach (alle / gasdominiert)")
    for sg in (0.0, 8.0, 10.0):
        oo = reste(sig=sg); nn = [n for n in oo if oo[n]['lVf'] == oo[n]['lVf']]
        r = np.array([oo[n]['rest'] for n in nn]); v = np.array([oo[n]['lVf'] for n in nn]); gd = np.array([oo[n]['fgas'] > 0.5 for n in nn])
        lowv = gd & (v < math.log10(60))
        print("   σ = %4.1f km/s: Median %+.3f, Streuung %.3f dex; ρ(Rest, log V_f) alle %+.3f (p %.1e), gasdominiert %+.3f (p %.1e); gasdominiert V_f < 60: %+.3f dex (n = %d)"
              % (sg, np.median(r), r.std(), *spearmanr(v, r), *spearmanr(v[gd], r[gd]), np.median(r[lowv]), lowv.sum()))
    for sg in (0.0, 8.0):
        oo = reste(sig=sg); nn = list(oo)
        r = np.array([oo[n]['rest'] for n in nn]); x = np.array([oo[n]['lSB'] for n in nn]); gd = np.array([oo[n]['fgas'] > 0.5 for n in nn])
        print("   σ = %4.1f km/s: ρ(Rest, log SB₀) alle %+.3f (p %.1e), gasdominiert %+.3f (p %.1e), sterndominiert %+.3f; Streuung alle %.3f dex"
              % (sg, *spearmanr(x, r), *spearmanr(x[gd], r[gd]), spearmanr(x[~gd], r[~gd])[0], r.std()))
    print("Gasdruck σ = 8 km/s mit äußerem Feld g_e (senkrecht): Zusammenhang mit log SB₀ und Streuung")
    for ge in (0.0, 0.01, 0.02, 0.03, 0.05, 0.08):
        oo = reste(sig=8.0, ge=ge); nn = list(oo)
        r = np.array([oo[n]['rest'] for n in nn]); x = np.array([oo[n]['lSB'] for n in nn]); gd = np.array([oo[n]['fgas'] > 0.5 for n in nn])
        print("   g_e/a₀ = %.2f: Median %+.3f, Streuung %.4f dex; ρ(Rest, log SB₀) alle %+.3f (p %.1e), gasdominiert %+.3f (p %.1e)"
              % (ge, np.median(r), r.std(), *spearmanr(x, r), *spearmanr(x[gd], r[gd])))
