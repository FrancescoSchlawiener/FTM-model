# -*- coding: utf-8 -*-
"""
Sehen die äußersten SPARC-Punkte ein äußeres Feld, und zeigen sie die Übergabe an R_Q?
Lesart (a): das großräumige, mitatmende Feld (w ≈ 1) geht in |g_N| der Galaxie ein → Rest fällt für a_N < g_e, Außenkurven fallen.
Lesart (b): nur das Feld innerhalb des eigenen umgekehrten Bereichs geht ein → in SPARC wirkt g_e = 0, Außenkurven flach.

Daten nur weltmodell/daten/sparc (SPARC_Lelli2016c.mrt, Rotmod_LTG); Auswahl, Υ und Gasdruck wie sparc_universal:
Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Υ_b = 0,7, V_c² = V_obs² + σ²·R/R_d (σ = 8 km/s; zum Vergleich σ = 0).
Regel: g = ν_eff(g_N; g_e)·g_N mit ν aus ftm_kern; Formen des äußeren Felds wie aussenfeld_kosmisch:
senkrecht |g| = √(g_N² + g_e²), gleichgerichtet |g| = g_N + g_e, Kugelmittel ⟨ν(|g_N + g_e|)⟩ über die Richtung.

1  Auswahl: Galaxien, deren äußerster Punkt a_N ≤ 0,05·a₀ hat; Punkte je Bereich a_N/a₀ < 0,01 / 0,01–0,03 / 0,03–0,1
2  Rest log(g_obs/g_Regel) je Bereich (Median, Bootstrap über Galaxien), roh und relativ zum eigenen Rest der Galaxie bei 0,1–1 a₀;
   Steigung des Rests gegen log a_N innerhalb der Galaxien (a_N < 0,1 a₀)
3  Steigung d ln v/d ln r der äußersten 3 Punkte, gemessen minus Regel (je g_e), gegen a_N/a₀ am letzten Punkt;
   Fehler: Bootstrap über Galaxien und Monte-Carlo über Abstand, Neigung, Υ (0,1 dex) und Geschwindigkeitsfehler
4  R_Q (a_N = s_kosm²·a₀ = 0,029 a₀): wie viele Kurven reichen hin, Sprung oder Knick des Rests dort?
6  Gegenprobe von Hand (ohne ftm_kern) an einer Galaxie
"""
import io, sys, os, math
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import A0, nu, kpc, q_KOSM
import sparc_flaechendichte as sf, sparc_universal as su

rng = np.random.default_rng(29092026)
GL_X, GL_W = np.polynomial.legendre.leggauss(32)
S_K2 = math.log(q_KOSM)**2                     # a_N/a₀ an R_Q
BER = ((0.0, 0.01), (0.01, 0.03), (0.03, 0.1))
GES = (0.003, 0.01, 0.03)
FORMEN = ("senkrecht", "gleich", "kugel")
LAB = dict(senkrecht="senkrecht", gleich="gleichgerichtet", kugel="Kugelmittel")


def nu_eff(gi, ge, form, a0=A0):
    gi = np.asarray(gi, float)
    if ge == 0:
        return nu(gi, a0)
    if form == "gleich":
        return nu(gi + ge, a0)
    if form == "senkrecht":
        return nu(np.hypot(gi, ge), a0)
    g = np.sqrt(np.maximum(gi[..., None]**2 + ge**2 + 2 * gi[..., None] * ge * GL_X, 1e-40))
    return 0.5 * np.sum(GL_W * nu(g, a0), axis=-1)


MODELLE = [("ohne", 0.0, "senkrecht")] + [(f"{LAB[f]} {ge:g}", ge, f) for f in FORMEN for ge in GES]

# ------------------------------------------------------------------ Daten
tab = su.tab_voll()
ALLE = []
for nm, d in tab.items():
    if d['Q'] > 2 or d['inc'] < 30:
        continue
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat'), comments='#')
    R, Vo, eV, Vg, Vd, Vb = a[:, 0] * kpc, a[:, 1] * 1e3, a[:, 2] * 1e3, a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3
    Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        continue
    ALLE.append(dict(name=nm, R=R[m], Vo=Vo[m], eV=np.maximum(eV[m], 1e3), Vg=Vg[m], Vd=Vd[m], Vb=Vb[m], Rd=d['Rd'],
                     eD=d['eD'] / d['D'], inc=d['inc'], eInc=d['eInc']))


def rechne(g, sig=8.0, fD=1.0, fi=1.0, fU=1.0, dV=None):
    """g_N, g_obs, R für eine (gestörte) Galaxie. Skalierung wie sparc_streuung."""
    R = g['R'] * fD
    Vo = g['Vo'] * fi if dV is None else (g['Vo'] + dV) * fi
    Vb2 = g['Vg'] * np.abs(g['Vg']) * fD + fU * fD * (0.5 * g['Vd']**2 + 0.7 * g['Vb']**2)
    gN = np.clip(Vb2, 1e-30, None) / R
    go = (Vo**2 + (sig * 1e3)**2 * (R / kpc) / (g['Rd'] * fD)) / R
    return R, gN, go


def steigung(x, y):
    x = np.log(x); y = np.log(y)
    return np.polyfit(x, y, 1)[0]


# ------------------------------------------------------------------ 1
def auswahl(sig=8.0):
    s = []
    for g in ALLE:
        R, gN, go = rechne(g, sig)
        if gN[-1] / A0 <= 0.05:
            s.append(g)
    return s


SEL = auswahl()
print(f"SPARC Q ≤ 2, Neigung ≥ 30°: {len(ALLE)} Galaxien; davon äußerster Punkt a_N ≤ 0,05·a₀: {len(SEL)}")
print(f"   (a_N am letzten Punkt; Υ_d = 0,5, Υ_b = 0,7; unabhängig von σ)")
print("1  Punkte je Bereich a_N/a₀ (in den ausgewählten Galaxien; in Klammern Zahl der Galaxien mit Punkten dort):")
for lo, hi in BER:
    n = 0; ng = 0
    for g in SEL:
        R, gN, go = rechne(g); k = ((gN / A0 >= lo) & (gN / A0 < hi)).sum(); n += k; ng += k > 0
    print(f"   {lo:.2f}–{hi:.2f}: {n:4d} Punkte ({ng} Galaxien)")
lastN = np.array([rechne(g)[1][-1] / A0 for g in SEL])
print(f"   a_N/a₀ am letzten Punkt: Median {np.median(lastN):.4f}, kleinster {lastN.min():.4f}; < 0,01: {np.sum(lastN<0.01)}, "
      f"0,01–0,03: {np.sum((lastN>=0.01)&(lastN<0.03))}, 0,03–0,05: {np.sum(lastN>=0.03)}")
Rl = np.array([g['R'][-1] / kpc for g in SEL]); RdS = np.array([g['Rd'] for g in SEL])
print(f"   letzter Radius: Median {np.median(Rl):.1f} kpc ({np.median(Rl/RdS):.1f} R_d), größter {Rl.max():.1f} kpc")
for ge in GES:
    print(f"   Zahl der Punkte mit a_N < g_e = {ge:g}·a₀: {sum(int((rechne(g)[1]/A0 < ge).sum()) for g in SEL)} "
          f"({sum(int((rechne(g)[1]/A0 < ge).any()) for g in SEL)} Galaxien)")


# ------------------------------------------------------------------ 2
def reste(gals, ge, form, sig=8.0, fU=1.0):
    """je Galaxie: log a_N/a₀, Rest, Referenzrest (Median bei 0,1–1 a₀)"""
    out = []
    for g in gals:
        R, gN, go = rechne(g, sig, fU=fU)
        r = np.log10(go / (nu_eff(gN, ge * A0, form) * gN))
        x = gN / A0
        ref = r[(x >= 0.1) & (x < 1.0)]
        out.append((np.log10(x), r, np.median(ref) if ref.size >= 2 else np.nan))
    return out


def binmed(out, lo, hi, rel):
    v = []
    for lx, r, ref in out:
        k = (lx >= math.log10(lo) if lo > 0 else np.ones_like(lx, bool)) & (lx < math.log10(hi))
        if rel:
            if not np.isfinite(ref):
                continue
            v.append(r[k] - ref)
        else:
            v.append(r[k])
    v = [x for x in v if x.size]
    return np.median(np.concatenate(v)) if v else np.nan, len(v)


def innere_steigung(out):
    """Steigung d Rest/d log a_N innerhalb der Galaxien, Punkte a_N < 0,1 a₀ (Galaxien mit ≥ 3 solchen Punkten), Mittel über Punkte"""
    X, Y = [], []
    for lx, r, ref in out:
        k = lx < -1
        if k.sum() >= 3:
            X.append(lx[k] - lx[k].mean()); Y.append(r[k] - r[k].mean())
    if not X:
        return np.nan
    X = np.concatenate(X); Y = np.concatenate(Y)
    return np.sum(X * Y) / np.sum(X * X)


def boot(fun, n, NB=1000):
    idx = [rng.integers(0, n, n) for _ in range(NB)]
    return np.nanstd([fun(i) for i in idx])


print("\n2  Rest log(g_obs/g_Regel) [dex], σ = 8 km/s; Median über Punkte je Bereich ± Bootstrap über Galaxien (1000)")
print("   roh | relativ zum eigenen Rest bei 0,1–1 a₀ | Steigung d Rest/d log a_N innerhalb der Galaxien (a_N < 0,1 a₀)")
print("   (a) verlangt: Modell mit dem richtigen g_e flach (relativ ≈ 0, Steigung ≈ 0); ohne g_e fällt der Rest unter g_e (Steigung > 0)")
ERG2 = {}
for sig in (8.0, 0.0):
    if sig == 0.0:
        print("   — zum Vergleich ohne Gasdruck (σ = 0):")
    for lab, ge, form in MODELLE:
        if sig == 0.0 and not (ge == 0 or form == "senkrecht"):
            continue
        out = reste(SEL, ge, form, sig)
        z = []
        for lo, hi in BER:
            m_raw, _ = binmed(out, lo, hi, False)
            e_raw = boot(lambda i: binmed([out[j] for j in i], lo, hi, False)[0], len(out), 300)
            m_rel, nrel = binmed(out, lo, hi, True)
            e_rel = boot(lambda i: binmed([out[j] for j in i], lo, hi, True)[0], len(out), 300)
            z.append(f"{m_raw:+.3f}±{e_raw:.3f} | {m_rel:+.3f}±{e_rel:.3f}")
        k = innere_steigung(out); ek = boot(lambda i: innere_steigung([out[j] for j in i]), len(out), 500)
        ERG2[(sig, lab)] = (k, ek)
        print(f"   {lab:22s}: <0,01: {z[0]};  0,01–0,03: {z[1]};  0,03–0,1: {z[2]};  Steigung {k:+.3f} ± {ek:.3f}")

print("   Steigung d Rest/d log a_N innerhalb der Galaxien mit Υ_d = 0,4 / 0,6 (Υ_b mitskaliert), σ = 8 km/s:")
for lab, ge, form in (MODELLE[0], MODELLE[2], MODELLE[3]):
    z = []
    for fU in (0.8, 1.2):
        out = reste(SEL, ge, form, 8.0, fU); k = innere_steigung(out)
        z.append(f"{k:+.3f} ± {boot(lambda i: innere_steigung([out[j] for j in i]), len(out), 300):.3f}")
    print(f"   {lab:22s}: " + " / ".join(z))
print("   Anzahl Galaxien mit Referenzrest (0,1–1 a₀, ≥ 2 Punkte): "
      f"{sum(np.isfinite(o[2]) for o in reste(SEL, 0, 'senkrecht'))} von {len(SEL)}")


# ------------------------------------------------------------------ 3
def aussen(g, ge, form, sig=8.0, **kw):
    R, gN, go = rechne(g, sig, **kw)
    R3, gN3, go3 = R[-3:], gN[-3:], go[-3:]
    so = steigung(R3, np.sqrt(go3 * R3))
    sm = steigung(R3, np.sqrt(nu_eff(gN3, ge * A0, form) * gN3 * R3))
    return so, sm, gN[-1] / A0


BER3 = ((0.0, 0.01), (0.01, 0.03), (0.03, 0.051))
print("\n3  Steigung d ln v/d ln r der äußersten 3 Punkte, gemessen (V_c mit Gasdruck) und Regel mit eigenen Baryonen, je Bereich a_N/a₀ am letzten Punkt")
print("   Median über Galaxien; Fehler: Bootstrap über Galaxien ⊕ Monte-Carlo (Abstand, Neigung, Υ 0,1 dex, Geschwindigkeitsfehler; 200 Ziehungen)")
nb = [sum(lo <= x < hi for x in lastN) for lo, hi in BER3]
print(f"   Galaxien je Bereich: <0,01: {nb[0]}, 0,01–0,03: {nb[1]}, 0,03–0,05: {nb[2]}")


def mc_zieh(g):
    fD = max(1 + rng.normal() * g['eD'], 0.2)
    inc = min(max(g['inc'] + rng.normal() * g['eInc'], 20.0), 90.0)
    fi = math.sin(math.radians(g['inc'])) / math.sin(math.radians(inc))
    fU = 10**(rng.normal() * 0.1)
    dV = rng.normal(size=g['Vo'].size) * g['eV']
    return dict(fD=fD, fi=fi, fU=fU, dV=dV)


NMC = 200
ZIEH = [[mc_zieh(g) for g in SEL] for _ in range(NMC)]
ERG3 = {}
for sig in (8.0, 0.0):
    if sig == 0.0:
        print("   — ohne Gasdruck (σ = 0):")
    for lab, ge, form in MODELLE:
        if sig == 0.0 and not (ge == 0 or (form == "senkrecht" and ge == 0.03)):
            continue
        A = np.array([aussen(g, ge, form, sig) for g in SEL])      # so, sm, x
        d = A[:, 0] - A[:, 1]
        nz = NMC if (ge == 0 or lab.startswith("senkrecht 0.03")) else 60
        MCD = np.array([[a_[0] - a_[1] for a_ in (aussen(g, ge, form, sig, **zz) for g, zz in zip(SEL, Z))] for Z in ZIEH[:nz]])
        z = []
        for lo, hi in BER3:
            k = (A[:, 2] >= lo) & (A[:, 2] < hi)
            med_o, med_m, med_d = np.median(A[k, 0]), np.median(A[k, 1]), np.median(d[k])
            eb = np.std([np.median(d[k][rng.integers(0, k.sum(), k.sum())]) for _ in range(1000)])
            kk = (lastN >= lo) & (lastN < hi)                      # Bereich nach ungestörtem a_N
            em = np.std(np.median(MCD[:, kk], axis=1))
            e = math.hypot(eb, em)
            z.append((med_o, med_m, med_d, eb, em, e))
        ERG3[(sig, lab)] = z
        print(f"   {lab:22s}: " + ";  ".join(f"gem {o:+.3f} Regel {m:+.3f} Δ {dd:+.3f}±{e:.3f} (B {eb:.3f}, MC {em:.3f})"
                                           for o, m, dd, eb, em, e in z))
# Trend Δ gegen log a_N (letzter Punkt), Theil-Sen, Bootstrap
from scipy.stats import theilslopes, spearmanr
print("   Trend von Δ = gemessen − Regel gegen log(a_N/a₀) am letzten Punkt (Theil-Sen, Bootstrap; Spearman):")
for lab, ge, form in MODELLE:
    A = np.array([aussen(g, ge, form) for g in SEL]); d = A[:, 0] - A[:, 1]; lx = np.log10(A[:, 2])
    k = theilslopes(d, lx)[0]
    bs = [theilslopes(d[i], lx[i])[0] for i in (rng.integers(0, len(d), len(d)) for _ in range(500))]
    rho, p = spearmanr(lx, d)
    print(f"   {lab:22s}: dΔ/d log a_N = {k:+.3f} ± {np.std(bs):.3f}; ρ = {rho:+.3f} (p = {p:.2f}); Median Δ alle {np.median(d):+.3f}")
# Güte: Summe der quadrierten Bereichsmediane Δ/Fehler
print("   χ² der drei Bereichsmediane (Δ/Fehler)², σ = 8 km/s:")
for lab, ge, form in MODELLE:
    z = ERG3[(8.0, lab)]
    print(f"   {lab:22s}: χ² = {sum((dd/e)**2 for o, m, dd, eb, em, e in z):6.2f}")


# ------------------------------------------------------------------ 4
print(f"\n4  R_Q: a_N/a₀ = s_kosm² = (ln q_kosm)² = {S_K2:.4f}; χ dort = s/(e^s − 1) = {math.sqrt(S_K2)/math.expm1(math.sqrt(S_K2)):.4f}")
reicht = [g for g in ALLE if rechne(g)[1][-1] / A0 < S_K2]
print(f"   Kurven, deren letzter Punkt jenseits R_Q liegt (a_N < {S_K2:.3f} a₀): {len(reicht)} von {len(ALLE)}")
spann = [g for g in ALLE if (rechne(g)[1] / A0 < S_K2).sum() >= 2 and ((rechne(g)[1] / A0 >= S_K2) & (rechne(g)[1] / A0 < 0.1)).sum() >= 2]
print(f"   Kurven mit ≥ 2 Punkten beiderseits (a_N < s_kosm² und s_kosm²–0,1): {len(spann)}")
for sig, lab, ge, form in ((8.0,) + MODELLE[0], (8.0, "senkrecht 0.03", 0.03, "senkrecht"), (0.0, "ohne, σ = 0", 0.0, "senkrecht")):
    dd, kn = [], []
    for g in spann:
        R, gN, go = rechne(g, sig); x = gN / A0
        r = np.log10(go / (nu_eff(gN, ge * A0, form) * gN))
        a_ = x < S_K2; b_ = (x >= S_K2) & (x < 0.1)
        dd.append(np.median(r[a_]) - np.median(r[b_]))
        # Knick: Steigung des Rests gegen log x beiderseits
        if a_.sum() >= 3 and b_.sum() >= 3:
            kn.append(np.polyfit(np.log10(x[a_]), r[a_], 1)[0] - np.polyfit(np.log10(x[b_]), r[b_], 1)[0])
    dd = np.array(dd); kn = np.array(kn)
    e1 = np.std([np.median(dd[rng.integers(0, len(dd), len(dd))]) for _ in range(1000)])
    e2 = np.std([np.median(kn[rng.integers(0, len(kn), len(kn))]) for _ in range(1000)]) if kn.size else np.nan
    print(f"   {lab:16s}: Sprung Rest(jenseits R_Q) − Rest(diesseits bis 0,1 a₀) Median {np.median(dd):+.3f} ± {e1:.3f} dex (n = {len(dd)}); "
          f"Knick der Steigung {np.median(kn):+.3f} ± {e2:.3f} (n = {len(kn)})")
# Regel selbst: d ln a/d ln a_N = 1 − χ/2 bei s_kosm² und bei 0,01
for x in (0.1, S_K2, 0.01):
    s = math.sqrt(x); ch = s / math.expm1(s)
    print(f"   Regel bei a_N = {x:.4f} a₀: χ = {ch:.4f}, d ln a/d ln a_N = {1-ch/2:.4f}, Punktmasse d ln v/d ln r = {-(1-ch/2)/1+0.5:+.4f}")


# ------------------------------------------------------------------ 6
print("\n6  Gegenprobe von Hand (eigener Code ohne ftm_kern, eine Galaxie):")
g = min(SEL, key=lambda g: g['R'][-1] * 0 + rechne(g)[1][-1])   # tiefster letzter Punkt
nm = g['name']
rows = [l.split() for l in open(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat')) if not l.startswith('#')]
rows = [[float(v) for v in r] for r in rows if r]
Rd = g['Rd']; KPC = 3.0856775814913673e19; a0 = 1.1e-10
print(f"   {nm}: R_d = {Rd} kpc; letzte 3 Zeilen (R, V_obs, V_gas, V_disk, V_bul):")
hand = []
for r in rows[-3:]:
    Rk, V, Vg_, Vd_, Vb_ = r[0], r[1], r[3], r[4], r[5]
    vb2 = Vg_ * abs(Vg_) + 0.5 * Vd_**2 + 0.7 * Vb_**2            # (km/s)²
    gn = vb2 * 1e6 / (Rk * KPC)
    x = gn / a0
    nu_h = 1.0 / (1.0 - math.exp(-math.sqrt(x)))
    vc2 = V**2 + 64.0 * Rk / Rd
    gobs = vc2 * 1e6 / (Rk * KPC)
    hand.append((Rk, math.sqrt(vc2), math.sqrt(nu_h * vb2)))
    print(f"   R = {Rk:6.2f}: V_obs = {V:6.2f}, V_bar² = {vb2:8.2f} (km/s)², a_N/a₀ = {x:.5f}, ν = {nu_h:.4f}, "
          f"V_c = {math.sqrt(vc2):6.2f}, V_Regel = {math.sqrt(nu_h*vb2):6.2f}, Rest = {math.log10(gobs/(nu_h*gn)):+.4f} dex")
lr = [math.log(h[0]) for h in hand]
for j, lab in ((1, "gemessen"), (2, "Regel")):
    lv = [math.log(h[j]) for h in hand]; mx = sum(lr) / 3; my = sum(lv) / 3
    k = sum((a - mx) * (b - my) for a, b in zip(lr, lv)) / sum((a - mx)**2 for a in lr)
    print(f"   Steigung {lab}: {k:+.4f}")
R, gN, go = rechne(g)
so, sm, x = aussen(g, 0.0, "senkrecht")
print(f"   dasselbe mit dem Skript (ftm_kern): a_N/a₀ letzter Punkt {x:.5f}, Steigung gemessen {so:+.4f}, Regel {sm:+.4f}, "
      f"Rest letzter Punkt {math.log10(go[-1]/(nu(gN[-1])*gN[-1])):+.4f}")
print(f"   Punkt 2 von Hand: s_kosm² = (ln(5,364/6,364))² = {math.log(0.12/0.02237/(1+0.12/0.02237))**2:.5f}")
