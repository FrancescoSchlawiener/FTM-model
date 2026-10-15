# -*- coding: utf-8 -*-
"""
Äußeres Feld als Einbettung in die kosmische Rückkehr?
Vermutung: Das mittlere äußere Feld jeder Galaxie ist das Feld, an dem die örtliche Rückkehr gleich der kosmischen wird,
g_e = s_kosm²·a₀ = (ln q_kosm)²·a₀,  1/(1 − q_kosm) = ν_kosm = 1 + Q.

1  q_kosm, s_kosm, g_e/a₀ für Q aus dem Kern und für Q = 2π − 1
2  SPARC (Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Gasdruck σ = 8 km/s wie sparc_universal): gemeinsames g_e frei (a₀ = 1,1 fest),
   Formen gleichgerichtet (s² = s_i² + s_e²), senkrecht (s⁴ = s_i⁴ + s_e⁴), Kugelmittel, Produkt der Rückkehr (s = s_i + s_e);
   Fehler durch Bootstrap über Galaxien
3  g_e fest (0; 0,02; s_kosm²; 0,04): gemeinsame Anpassung je Drittel von Σ₀ (Abstand, Neigung als Störgrößen, Υ_d fest)
   und Rangkorrelation Rest gegen SB₀ (gasreich)
4  Zusammensetzung aus q_lokal und q_kosm (Verlustraten, Produkt, Fluchtanteile)
5  Milchstraße: R_Q und Rotation jenseits davon
6  Gegenproben von Hand (eigener Code ohne den Kern für ν)
"""
import io, sys, os, math, glob
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import spearmanr
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import A0, nu, G, M_sun, kpc, Q_KOSM, q_KOSM, r_stern
import sparc_flaechendichte as sf, sparc_universal as su

rng = np.random.default_rng(29092026)
GL_X, GL_W = np.polynomial.legendre.leggauss(32)          # Kugelmittel über μ = cos θ


# ------------------------------------------------------------------ 1
def kosm(Q):
    q = Q / (1 + Q); s = -math.log(q)
    return q, s, s * s


print("1  Kosmische Rückkehr: q_kosm = Q/(1 + Q) aus 1/(1 − q) = 1 + Q")
S_K = {}
for lab, Q in (("Q aus dem Kern (Ω_c h²/Ω_b h²)", Q_KOSM), ("Q = 2π − 1 (Zählung je Radiant)", 2 * math.pi - 1)):
    q, s, ge = kosm(Q)
    S_K[lab] = s
    print(f"   {lab:34s}: Q = {Q:.4f}, q_kosm = {q:.5f}, 1/(1−q) = {1/(1-q):.4f}, s_kosm = −ln q = {s:.5f}, g_e/a₀ = s² = {ge:.5f}")
GE_K = kosm(Q_KOSM)[2]
assert abs(q_KOSM - Q_KOSM / (1 + Q_KOSM)) < 1e-12


# ------------------------------------------------------------------ Formen
def nu_eff(gi, ge, form, a0=A0):
    gi = np.asarray(gi, float)
    if ge == 0:
        return nu(gi, a0)
    if form == "gleich":
        return nu(gi + ge, a0)
    if form == "senkrecht":
        return nu(np.hypot(gi, ge), a0)
    if form == "kugel":
        g = np.sqrt(np.maximum(gi[..., None]**2 + ge**2 + 2 * gi[..., None] * ge * GL_X, 1e-40))
        return 0.5 * np.sum(GL_W * nu(g, a0), axis=-1)
    if form == "produkt":
        s = np.sqrt(gi / a0) + math.sqrt(ge / a0)
        return 1 / (-np.expm1(-s))
    raise ValueError(form)


FORMEN = ("gleich", "senkrecht", "kugel", "produkt")
LAB = dict(gleich="gleichgerichtet s²=s_i²+s_e²", senkrecht="senkrecht s⁴=s_i⁴+s_e⁴", kugel="Kugelmittel über μ",
           produkt="Produkt q_i·q_e, s=s_i+s_e")

# ------------------------------------------------------------------ Daten (wie sparc_universal / schwelle_je_galaxie)
tab = su.tab_voll()
gal = []
for nm, d in tab.items():
    if d['Q'] > 2 or d['inc'] < 30:
        continue
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat'), comments='#')
    R, Vo, eV, Vg, Vd, Vb = a[:, 0] * kpc, a[:, 1] * 1e3, a[:, 2] * 1e3, a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3
    Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        continue
    Mb = 0.5 * d['L'] + 1.33 * d['MHI']
    gal.append(dict(name=nm, R=R[m], Vo=Vo[m], eV=np.maximum(eV[m], 1.0), Vg=Vg[m], Vd=Vd[m], Vb=Vb[m], Rd=d['Rd'],
                    gN=Vb2[m] / R[m], lSB=math.log10(d['SB0']) if d['SB0'] > 0 else np.nan, fgas=1.33 * d['MHI'] / Mb,
                    eD=d['eD'] / d['D'], inc=math.radians(d['inc']), ei=math.radians(max(d['eInc'], 1.0)), SB0=d['SB0']))
SIG = 8.0e3
for g in gal:
    g['go'] = (g['Vo']**2 + SIG**2 * (g['R'] / kpc) / g['Rd']) / g['R']
    g['w'] = (np.log(10) / (2 * g['eV'] / g['Vo']))**2              # Gewicht für log10 aus Geschwindigkeitsfehler
print(f"\nSPARC Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Υ_b = 0,7, Gasdruck σ = 8 km/s: {len(gal)} Galaxien, {sum(len(g['R']) for g in gal)} Punkte")

# ------------------------------------------------------------------ 2
GEG = np.round(np.arange(0.0, 0.1501, 0.0025), 4)


def tabellen(form):
    """je Galaxie und g_e: Median des Rests (dex) und gewichtete Fehlerquadratsumme der Punkte"""
    med = np.zeros((len(gal), len(GEG))); c2 = np.zeros_like(med)
    for i, g in enumerate(gal):
        for j, ge in enumerate(GEG):
            r = np.log10(g['go'] / (nu_eff(g['gN'], ge * A0, form) * g['gN']))
            med[i, j] = np.median(r); c2[i, j] = np.sum(g['w'] * r**2)
    return med, c2


def minimum(y):
    k = int(np.argmin(y))
    if 0 < k < len(y) - 1:
        y0, y1, y2 = y[k - 1], y[k], y[k + 1]; h = GEG[1] - GEG[0]
        den = y0 - 2 * y1 + y2
        return GEG[k] + (h * (y0 - y2) / (2 * den) if den > 0 else 0.0)
    return GEG[k]


print("\n2  Gemeinsames äußeres Feld frei (a₀ = 1,1·10⁻¹⁰ fest); Fehler: Bootstrap über Galaxien (1000)")
print("   A: jede Galaxie gleich (Σ Median²)   B: Punkte mit Geschwindigkeitsfehlern (Σ w r²)")
TAB = {}
for form in FORMEN:
    med, c2 = tabellen(form); TAB[form] = (med, c2)
    yA = np.sum(med**2, axis=0); yB = np.sum(c2, axis=0)
    bA, bB = [], []
    for _ in range(1000):
        i = rng.integers(0, len(gal), len(gal))
        bA.append(minimum(np.sum(med[i]**2, axis=0))); bB.append(minimum(np.sum(c2[i], axis=0)))
    gA, gB = minimum(yA), minimum(yB)
    j0 = 0; jk = int(np.argmin(np.abs(GEG - 0.03)))
    print(f"   {LAB[form]:32s}: A g_e/a₀ = {gA:.4f} ± {np.std(bA):.4f} (16–84 %: {np.percentile(bA,16):.4f}–{np.percentile(bA,84):.4f}); "
          f"B g_e/a₀ = {gB:.4f} ± {np.std(bB):.4f};  ΔΣMedian² (0 → 0,03) = {yA[jk]-yA[j0]:+.3f}, Δχ²_B (0 → 0,03) = {yB[jk]-yB[j0]:+.0f}")
print(f"   Vergleich: s_kosm² = {GE_K:.4f} (Q Kern), {kosm(2*math.pi-1)[2]:.4f} (Q = 2π − 1)")
print("   Mittel des Rests (Median über Galaxien, dex) bei g_e = 0 / 0,02 / s_kosm² / 0,04 und Güte A:")
for form in FORMEN:
    med = TAB[form][0]
    z = []
    for ge in (0.0, 0.02, GE_K, 0.04):
        j = int(np.argmin(np.abs(GEG - ge)))
        z.append(f"{np.median(med[:, j]):+.3f} ({np.sum(med[:, j]**2):.3f})")
    print(f"   {LAB[form]:32s}: " + " / ".join(z))

# Gegenweg zu 2: g_e und a₀ zugleich frei (Form gleichgerichtet, Maß A)
AGR = np.exp(np.linspace(math.log(0.7e-10), math.log(1.7e-10), 41))
best = None
for form in ("gleich", "senkrecht"):
    Z = np.zeros((len(AGR), len(GEG)))
    for ia, a0 in enumerate(AGR):
        for j, ge in enumerate(GEG[::2]):
            Z[ia, 2 * j] = sum(np.median(np.log10(g['go'] / (nu_eff(g['gN'], ge * A0, form, a0) * g['gN'])))**2 for g in gal)
    Zs = Z[:, ::2]; ia, j = np.unravel_index(np.argmin(Zs), Zs.shape)
    print(f"   a₀ und g_e zugleich frei ({LAB[form]}, Maß A): a₀ = {AGR[ia]*1e10:.2f}·10⁻¹⁰, g_e/a₀ = {GEG[::2][j]:.3f} "
          f"(Güte {Zs[ia, j]:.3f}; bei a₀ = 1,1, g_e = 0: {sum(np.median(np.log10(g['go']/(nu(g['gN'])*g['gN'])))**2 for g in gal):.3f})")
for form in ("gleich", "senkrecht", "kugel"):
    zz = [sum(np.median(np.log10(g['go'] / (nu_eff(g['gN'], GE_K * A0, form, a0) * g['gN'])))**2 for g in gal) for a0 in AGR]
    k = int(np.argmin(zz))
    print(f"   g_e = s_kosm²·a₀ fest, a₀ frei ({LAB[form]}): a₀ = {AGR[k]*1e10:.2f}·10⁻¹⁰, Güte {zz[k]:.3f}")


# ------------------------------------------------------------------ 3
print("\n3  g_e fest: Anstieg der Schwelle mit Σ₀ (gemeinsam je Drittel) und Rangkorrelation Rest–SB₀ (gasreich)")
GEFEST = (0.0, 0.02, GE_K, 0.04)
for form in ("senkrecht", "gleich", "kugel"):
    med = TAB[form][0]
    z = []
    for ge in GEFEST:
        # Rest bei festem g_e direkt (nicht über Gitter)
        r = np.array([np.median(np.log10(g['go'] / (nu_eff(g['gN'], ge * A0, form) * g['gN']))) for g in gal])
        x = np.array([g['lSB'] for g in gal]); gd = np.array([g['fgas'] > 0.5 for g in gal]); ok = np.isfinite(x)
        ra = spearmanr(x[ok], r[ok])[0]; rg, pg = spearmanr(x[ok & gd], r[ok & gd])
        z.append(f"{ge:.4f}: alle {ra:+.3f}, gasreich {rg:+.3f} (p {pg:.1e})")
    print(f"   {LAB[form]:32s}: " + "; ".join(z))

# gemeinsame Anpassung je Drittel (Abstand, Neigung als Störgrößen, Υ_d fest) wie schwelle_je_galaxie, ohne Gasdruck
sel = [g for g in gal if g['SB0'] > 0]
AG = np.exp(np.linspace(math.log(0.3e-10), math.log(4e-10), 49))


def resid(par, g, a0, ge, form):
    uD, di = par
    inc = min(max(g['inc'] + di, 0.2), math.pi / 2)
    gb = np.maximum(g['gN'], 1e-15)
    gobs = np.log(g['Vo']**2 / g['R']) - uD - 2 * math.log(math.sin(g['inc']) / math.sin(inc))
    r = (gobs - np.log(nu_eff(gb, ge * A0, form, a0) * gb)) / (2 * g['eV'] / g['Vo'])
    return np.r_[r, uD / g['eD'], di / g['ei']]


def chi2_tafel(ge, form):
    T = np.zeros((len(sel), len(AG)))
    for i, g in enumerate(sel):
        p0 = [0.0, 0.0]
        for k, a0 in enumerate(AG):
            s = least_squares(resid, p0, args=(g, a0, ge, form)); p0 = s.x
            T[i, k] = 2 * s.cost
    return T


def amin(y):
    x = np.log(AG); k = int(np.argmin(y))
    if 0 < k < len(y) - 1:
        y0, y1, y2 = y[k - 1], y[k], y[k + 1]
        return x[k] + (x[1] - x[0]) * (y0 - y2) / (2 * (y0 - 2 * y1 + y2))
    return x[k]


xs = np.log10(np.array([0.5 * g['SB0'] for g in sel]))
t1, t2 = np.percentile(xs, [100 / 3, 200 / 3]); grp = np.digitize(xs, [t1, t2])


def drittel(T, idx=None):
    idx = np.arange(len(sel)) if idx is None else idx
    a, xm = [], []
    for k in range(3):
        ii = idx[grp[idx] == k]
        a.append(amin(T[ii].sum(axis=0)) / math.log(10)); xm.append(xs[ii].mean())
    return np.array(a), np.array(xm)


for form in ("senkrecht", "gleich"):
    for ge in GEFEST:
        if ge == 0 and form == "gleich":
            continue
        T = chi2_tafel(ge, form)
        a, xm = drittel(T); k = np.polyfit(xm, a, 1)[0]
        bs = []
        for _ in range(300):
            i = rng.integers(0, len(sel), len(sel)); aa, xx = drittel(T, i); bs.append(np.polyfit(xx, aa, 1)[0])
        print(f"   {('g_e = 0' if ge == 0 else LAB[form]):32s} g_e/a₀ = {ge:.4f}: a₀ je Drittel = " + ", ".join(f"{10**v*1e10:.2f}" for v in a)
              + f" ·10⁻¹⁰ (log Σ₀ = " + ", ".join(f"{v:.2f}" for v in xm) + f"); Anstieg {k:+.3f} ± {np.std(bs):.3f} je Dekade")


# ------------------------------------------------------------------ 4
print("\n4  Zusammensetzung der Rückkehr aus örtlichem und kosmischem Anteil (s_k = s_kosm)")
sk = S_K["Q aus dem Kern (Ω_c h²/Ω_b h²)"]
print("   Verlustraten addieren (Γ = Γ_i + Γ_k, s = √(Γτ_D)):   s² = s_i² + s_k²  = gleichgerichtete Form mit g_e = s_k²·a₀")
print("   Rückkehrwahrscheinlichkeiten multiplizieren (q = q_i·q_k): s = s_i + s_k")
print("   Fluchtanteile addieren (1 − q = (1 − q_i) + (1 − q_k)):   1/ν = 1/ν_i + 1/ν_k")
for si in (3.0, 1.0, 0.3, 0.1, 0.0):
    qi, qk = math.exp(-si), math.exp(-sk)
    n_rate = 1 / (1 - math.exp(-math.sqrt(si**2 + sk**2)))
    n_prod = 1 / (1 - qi * qk)
    n_flucht = 1 / ((1 - qi) + (1 - qk))
    n_frei = 1 / (1 - qi) if si > 0 else float('inf')
    print(f"   s_i = {si:4.1f}: ν frei {n_frei:8.4f} | Raten {n_rate:.4f} | Produkt {n_prod:.4f} | Fluchtanteile {n_flucht:.4f}")
print(f"   Grenze s_i → 0: alle drei → ν_kosm = 1 + Q = {1+Q_KOSM:.4f} (Fluchtanteile nur dort); Newton (s_i → ∞): Raten, Produkt → 1, Fluchtanteile → {1/(1+1-math.exp(-sk)):.3f} < 1")

# ------------------------------------------------------------------ 5
print("\n5  Milchstraße als Punktmasse, 6·10¹⁰ M☉")
M = 6e10 * M_sun; rs = r_stern(M)
for lab, s in S_K.items():
    print(f"   {lab:34s}: r* = {rs/kpc:.3f} kpc, R_Q = r*/s_kosm = {rs/s/kpc:.2f} kpc")
print("   v [km/s] bei r [kpc], g_e = s_kosm²·a₀:  frei | gleichgerichtet | senkrecht | Kugelmittel | Produkt")
for r in (20, 51, 100, 200, 300, 500, 1000):
    gN = G * M / (r * kpc)**2
    v = [math.sqrt(float(nu(gN)) * gN * r * kpc) / 1e3] + [math.sqrt(float(nu_eff(np.array([gN]), GE_K * A0, f)[0]) * gN * r * kpc) / 1e3 for f in FORMEN]
    print(f"   r = {r:5d}: " + " | ".join(f"{x:6.1f}" for x in v) + f"   (ν_kosm-Kepler √((1+Q)GM/r) = {math.sqrt((1+Q_KOSM)*gN*r*kpc)/1e3:.1f})")
for f in ("gleich", "kugel"):
    rr = np.array([100, 200]) * kpc; gN = G * M / rr**2
    v = np.sqrt(nu_eff(gN, GE_K * A0, f) * gN * rr)
    print(f"   {LAB[f]}: d ln v/d ln r zwischen 100 und 200 kpc = {math.log(v[1]/v[0])/math.log(2):+.3f} (flach 0, Kepler −0,5)")

# ------------------------------------------------------------------ 6
print("\n6  Gegenproben (eigener Code, ν von Hand: 1/(1 − exp(−√x)))")
Qh = 0.1200 / 0.02237; qh = Qh / (1 + Qh); print(f"   von Hand: Q = {Qh:.4f}, q = {qh:.5f}, ln q = {math.log(qh):.5f}, (ln q)² = {math.log(qh)**2:.5f}; "
                                              f"2π − 1: q = 1 − 1/2π = {1-1/(2*math.pi):.5f}, (ln q)² = {math.log(1-1/(2*math.pi))**2:.5f}")
GM = 6.6743e-11 * 6e10 * 1.98841e30; rsh = math.sqrt(GM / 1.1e-10) / 3.0857e19
print(f"   von Hand: r* = √(GM/a₀) = {rsh:.3f} kpc, R_Q = {rsh/abs(math.log(qh)):.2f} kpc; a_N(R_Q)/a₀ = {GM/(rsh/abs(math.log(qh))*3.0857e19)**2/1.1e-10:.5f}")


def nu_h(x):
    return 1.0 / (1.0 - math.exp(-math.sqrt(x)))


for nm in ("DDO154", "NGC3198", "UGC00128"):
    g = next((x for x in gal if x['name'] == nm), None)
    if g is None:
        continue
    best = (1e9, None)
    for k in range(0, 201):
        ge = 0.001 * k
        rr = sorted(math.log10(go / (nu_h((gi + ge * 1.1e-10) / 1.1e-10) * gi)) for gi, go in zip(g['gN'], g['go']))
        n = len(rr); md = rr[n // 2] if n % 2 else 0.5 * (rr[n // 2 - 1] + rr[n // 2])
        if abs(md) < best[0]:
            best = (abs(md), ge)
    j = int(np.argmin(np.abs(TAB['gleich'][0][[x['name'] for x in gal].index(nm)])))
    print(f"   {nm:9s} Hand-Gitter (gleichgerichtet, Median-Rest = 0): g_e/a₀ = {best[1]:.3f}; Gitter des Skripts: {GEG[j]:.4f}")
