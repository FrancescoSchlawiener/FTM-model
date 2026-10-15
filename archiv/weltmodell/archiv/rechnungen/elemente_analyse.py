# Elemente und Übergänge: Phasen, Elektronen, Kern an allen Elementen mit Daten; Anschluss an das Modell mit Zufallsschutz
# Daten: daten/elemente/pubchem_periodensystem.csv (PubChem), daten/elemente/mass_1.mas20.txt (AME2020, nur Werte ohne '#' gelten als gemessen)
# 1 Phasen  2 Elektronen  3 Kern  4 Trichtergrößen  5 Regeln des Modells an den Daten  6 Musterscan mit Zufallsschutz
import numpy as np, csv, sys, io, os, math
from itertools import combinations
from scipy import stats
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
DAT = os.path.join(HERE, "..", "daten", "elemente")
rng = np.random.default_rng(20260929)
out = []; P = lambda s="": out.append(s)
NB, NP = 2000, 5000            # Bootstrap- und Permutationszahl

# Konstanten (CODATA 2018)
G = 6.67430e-11; c = 299792458.0; hbar = 1.054571817e-34; kB = 1.380649e-23; NA = 6.02214076e23
me = 9.1093837015e-31; mu = 1.66053906660e-27; u_keV = 931494.10242; re = 2.8179403262e-15
aB = 5.29177210903e-11; lam_e = 3.8615926796e-13; eV = 1.602176634e-19; a0 = 1.2e-10
mPl = math.sqrt(hbar * c / G); lP = math.sqrt(hbar * G / c**3)
MEH, MEn = 7288.971064, 8071.31806   # Massenüberschuss ¹H und n [keV], AME2020
E2 = 1.439964547                      # e²/(4πε₀) [eV·nm]

# ---------------- Hilfen ----------------
def fnum(x):
    try: return float(x)
    except: return np.nan

def spear(x, y, nb=NB, npm=NP):
    """Spearman ρ, 95-%-Bootstrap-Intervall, Permutations-p (zweiseitig)."""
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]; n = len(x)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    r = np.corrcoef(rx, ry)[0, 1]
    bs = []
    for _ in range(nb):
        i = rng.integers(0, n, n)
        a, b = stats.rankdata(x[i]), stats.rankdata(y[i])
        if a.std() > 0 and b.std() > 0: bs.append(np.corrcoef(a, b)[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    ryp = np.array([rng.permutation(ry) for _ in range(npm)])
    rxs = (rx - rx.mean()) / rx.std(); ryps = (ryp - ryp.mean(1, keepdims=True)) / ryp.std(1, keepdims=True)
    rp = (ryps * rxs).mean(1)
    p = (np.sum(np.abs(rp) >= abs(r) - 1e-12) + 1) / (npm + 1)
    return r, lo, hi, p, n

def fmt_p(p):
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}"

def period_group(z):
    ends = [2, 10, 18, 36, 54, 86, 118]; start = [0] + ends
    p = next(i + 1 for i, e in enumerate(ends) if z <= e); k = z - start[p - 1]
    if p == 1: return p, (1 if z == 1 else 18)
    if p in (2, 3): return p, (k if k <= 2 else k + 10)
    if p in (4, 5): return p, k
    if k <= 2: return p, k
    if k == 3: return p, 3            # La, Ac in Gruppe 3 (Konvention [A])
    if k <= 17: return p, 0           # Ce–Lu, Th–Lr: f-Block, keine Gruppe
    return p, k - 14

# ---------------- Daten laden ----------------
rows = list(csv.DictReader(open(os.path.join(DAT, "pubchem_periodensystem.csv"), encoding="utf-8")))
Z = np.array([int(r["AtomicNumber"]) for r in rows]); sym = [r["Symbol"] for r in rows]
S = {s: i for i, s in enumerate(sym)}
M = np.array([fnum(r["AtomicMass"]) for r in rows]); EN = np.array([fnum(r["Electronegativity"]) for r in rows])
RV = np.array([fnum(r["AtomicRadius"]) for r in rows]); IE = np.array([fnum(r["IonizationEnergy"]) for r in rows])
EA = np.array([fnum(r["ElectronAffinity"]) for r in rows]); Tm = np.array([fnum(r["MeltingPoint"]) for r in rows])
Tb = np.array([fnum(r["BoilingPoint"]) for r in rows]); RHO = np.array([fnum(r["Density"]) for r in rows])
STATE = [r["StandardState"] for r in rows]; BLOCK = [r["GroupBlock"] for r in rows]; CONF = [r["ElectronConfiguration"] for r in rows]
PER = np.array([period_group(z)[0] for z in Z]); GRP = np.array([period_group(z)[1] for z in Z])
COND = np.array([s in ("Solid", "Liquid") for s in STATE]) & np.isfinite(RHO)
VA = np.where(COND, M / (RHO * NA) * 1e24, np.nan)          # Volumen je Atom im kondensierten Zustand [Å³]
METAL = np.array([b in ("Alkali metal", "Alkaline earth metal", "Transition metal", "Lanthanide", "Actinide", "Post-transition metal") for b in BLOCK])
NONMET = np.array([b in ("Nonmetal", "Halogen", "Noble gas") for b in BLOCK]); METALLOID = np.array([b == "Metalloid" for b in BLOCK])
QUIRK = {"He": "Schmelzpunkt nur unter Druck", "C": "sublimiert (Schmelz- und Siedewert fast gleich)", "As": "sublimiert (Schmelzpunkt über Siedepunkt)"}
okR = np.isfinite(Tm) & np.isfinite(Tb) & np.array([s not in QUIRK for s in sym])

P("ELEMENTE UND ÜBERGÄNGE: Phasen, Elektronen, Kern; Anschluss an das Modell")
P(f"PubChem: 118 Elemente; mit Schmelzpunkt {np.isfinite(Tm).sum()}, Siedepunkt {np.isfinite(Tb).sum()}, beide {(np.isfinite(Tm)&np.isfinite(Tb)).sum()}, "
  f"Ionisierungsenergie {np.isfinite(IE).sum()}, Elektronenaffinität {np.isfinite(EA).sum()}, Elektronegativität {np.isfinite(EN).sum()}, "
  f"Radius {np.isfinite(RV).sum()}, Dichte kondensiert {COND.sum()}")
P("Radius in PubChem = Van-der-Waals-Radius [pm]. Dichte der Gase ist Gasdichte, darum Atomvolumen nur für feste und flüssige Elemente.")
P("Sonderfälle (aus dem Verhältnis T_b/T_m genommen): " + "; ".join(f"{k}: {v}" for k, v in QUIRK.items()))
P()

# =====================================================================================
P("1) PHASEN: Schmelz- und Siedepunkt")
P("-" * 90)
ratio = np.where(okR, Tb / Tm, np.nan)
P(f"1a) Verhältnis T_b/T_m über {okR.sum()} Elemente: Median {np.nanmedian(ratio):.2f}, 16–84 % {np.nanpercentile(ratio,16):.2f}–{np.nanpercentile(ratio,84):.2f}, "
  f"Min {np.nanmin(ratio):.3f} ({sym[np.nanargmin(ratio)]}), Max {np.nanmax(ratio):.2f} ({sym[np.nanargmax(ratio)]}); Variationskoeffizient {np.nanstd(ratio)/np.nanmean(ratio):.2f}")
groups = {}
for i in np.where(okR)[0]: groups.setdefault(BLOCK[i], []).append(i)
P("    je Block:  n   Median   Min            Max")
kw = []
for b, idx in sorted(groups.items(), key=lambda t: np.median(ratio[t[1]])):
    rr = ratio[idx]; imin, imax = idx[int(np.argmin(rr))], idx[int(np.argmax(rr))]
    P(f"    {b:24s} {len(idx):3d}  {np.median(rr):6.2f}   {rr.min():.2f} ({sym[imin]:2s})     {rr.max():.2f} ({sym[imax]})")
    if len(idx) >= 3: kw.append(np.log(rr))
H, pkw = stats.kruskal(*kw)
P(f"    Kruskal-Wallis über Blöcke mit n ≥ 3: H = {H:.1f}, p = {fmt_p(pkw)}. Das Fenster T_b/T_m hängt an der Bindungsart, es ist keine Konstante.")
P(f"    [G] Richards (ΔS_m ≈ k je Atom) und Trouton (ΔS_b ≈ 10,5 k): ΔH_b/ΔH_m ≈ 10,5·T_b/T_m = {10.5*np.nanmedian(ratio):.0f} im Median;")
P(f"        Schmelzen nimmt nur etwa {100/(10.5*np.nanmedian(ratio)):.0f} % der Bindung, es löst die Ordnung, nicht den Zusammenhalt.")
P()

P("1b) Was geht mit ln T_m und ln T_b einher? Spearman ρ [95 % Bootstrap], Permutations-p; Bonferroni über alle 16 Tests")
feat = [("Ionisierungsenergie", IE), ("Elektronegativität", EN), ("Elektronenaffinität", EA), ("Van-der-Waals-Radius", RV),
        ("Atomvolumen (kondens.)", VA), ("Dichte (kondens.)", np.where(COND, RHO, np.nan)), ("Masse", M), ("T_m bzw. T_b", None)]
for lab, T in (("ln T_m", Tm), ("ln T_b", Tb)):
    for fn, F in feat:
        if F is None: F = Tb if lab == "ln T_m" else Tm
        r, lo, hi, p, n = spear(F, np.log(T))
        sig = "  *" if p * 16 < 0.05 else ""
        P(f"    {lab} ~ {fn:24s} n={n:3d}  ρ = {r:+.2f} [{lo:+.2f}, {hi:+.2f}]  p = {fmt_p(p)}{sig}")
P("    * = nach Bonferroni signifikant.")
# Wie viel legen Werte des einzelnen Atoms fest, wie viel das Atomvolumen im Festkörper? Leave-one-out-Vorhersage von ln T_m
def loo_r2(Xm, y):
    n = len(y); pred = np.empty(n)
    for i in range(n):
        k = np.arange(n) != i; A_ = np.column_stack([np.ones(k.sum()), Xm[k]])
        cf = np.linalg.lstsq(A_, y[k], rcond=None)[0]; pred[i] = cf[0] + Xm[i] @ cf[1:]
    return 1 - np.sum((y - pred)**2) / np.sum((y - y.mean())**2)
isd = ((GRP >= 3) & (GRP <= 12)).astype(float)
dpar = isd * (GRP - 6.0)**2        # Füllungsparabel, Scheitel fest bei halber Füllung (6 von 12), nicht angepasst
sets = [("Atom allein: IE, EN, r", [IE, EN, RV]), ("+ ln V_A (Festkörper)", [IE, EN, RV, np.log(VA)]),
        ("+ d-Füllung", [IE, EN, RV, np.log(VA), isd, dpar])]
mall = np.isfinite(Tm) & np.isfinite(IE) & np.isfinite(EN) & np.isfinite(RV) & np.isfinite(VA)
y = np.log(Tm[mall]); line = []
for lab, F in sets:
    Xm = np.column_stack([f[mall] for f in F]); r2 = loo_r2(Xm, y)
    null = [loo_r2(Xm, rng.permutation(y)) for _ in range(300)]
    line.append(f"{lab}: {r2:.2f} (Zufall ≤ {np.percentile(null, 95):.2f})")
P(f"    Vorhersage von ln T_m (Leave-one-out-R², n = {mall.sum()}): " + "; ".join(line))
P()

P("1c) d-Block: T_m gegen Gruppe 3–12 (Zahl der s+d-Elektronen), je Reihe quadratisch angepasst")
Rsq = []; peaks = []
for per in (4, 5, 6):
    idx = [i for i in range(len(Z)) if PER[i] == per and 3 <= GRP[i] <= 12 and np.isfinite(Tm[i])]
    g = GRP[idx].astype(float); t = Tm[idx]
    cf = np.polyfit(g, t, 2); fit = np.polyval(cf, g); r2 = 1 - np.sum((t - fit)**2) / np.sum((t - t.mean())**2)
    gpk = -cf[1] / (2 * cf[0]); Rsq.append(r2); peaks.append(gpk)
    imax = idx[int(np.argmax(t))]
    P(f"    Reihe {per}: " + " ".join(f"{sym[i]} {Tm[i]:.0f}" for i in idx))
    P(f"            R² = {r2:.2f}, Scheitel bei Gruppe {gpk:.2f}; höchster Wert {sym[imax]} (Gruppe {GRP[imax]})")
# Permutation: Werte innerhalb jeder Reihe mischen, mittleres R² der Parabel
obs = np.mean(Rsq); cnt = 0
rowsd = []
for per in (4, 5, 6):
    idx = [i for i in range(len(Z)) if PER[i] == per and 3 <= GRP[i] <= 12 and np.isfinite(Tm[i])]
    rowsd.append((GRP[idx].astype(float), Tm[idx]))
for _ in range(NP):
    rs = []
    for g, t in rowsd:
        tp = rng.permutation(t); cf = np.polyfit(g, tp, 2); fit = np.polyval(cf, g)
        rs.append(1 - np.sum((tp - fit)**2) / np.sum((tp - tp.mean())**2))
    cnt += np.mean(rs) >= obs
P(f"    mittleres R² = {obs:.2f}, Permutations-p = {fmt_p((cnt+1)/(NP+1))}; Scheitel im Mittel bei Gruppe {np.mean(peaks):.2f} = halbe Füllung von 12 s+d-Plätzen (Friedel, bekannt)")
P()

P("1d) Größte Sprünge zum Nachbarn Z → Z+1")
for lab, T in (("T_m", Tm), ("T_b", Tb)):
    jumps = []
    for i in range(len(Z) - 1):
        if np.isfinite(T[i]) and np.isfinite(T[i + 1]) and sym[i] not in ("He",) and sym[i + 1] not in ("He",):
            jumps.append((abs(np.log(T[i + 1] / T[i])), i))
    jumps.sort(reverse=True)
    P(f"    {lab}: " + ", ".join(f"{sym[i]}→{sym[i+1]} {T[i]:.0f}→{T[i+1]:.0f} K" for _, i in jumps[:8]))
P("    Die größten Sprünge liegen dort, wo ein Netz fester Bindungen (Metall, kovalentes Gitter) in abgeschlossene Einheiten (Edelgasatom, Molekül N₂, O₂, F₂) übergeht,")
P("    oder umgekehrt. Der Übergang zählt die Bindungen zwischen den Einheiten, nicht die in ihnen.")
P()

P("1e) Wer ist wem ähnlicher: waagrechter Nachbar (Z, Z+1, gleiche Periode) oder senkrechter (gleiche Gruppe, nächste Periode)?  mittleres |Δ|, Mann-Whitney")
def neigh_pairs(sel):
    hor = [(i, i + 1) for i in range(len(Z) - 1) if PER[i] == PER[i + 1] and sel(GRP[i]) and sel(GRP[i + 1]) and GRP[i] > 0 and GRP[i + 1] > 0]
    ver = []
    for i in range(len(Z)):
        for j in range(i + 1, len(Z)):
            if GRP[i] == GRP[j] and GRP[i] > 0 and PER[j] == PER[i] + 1 and sel(GRP[i]): ver.append((i, j))
    return hor, ver
for blab, sel in (("Hauptgruppen (1, 2, 13–18)", lambda g: g in (1, 2) or g >= 13), ("d-Block (3–12)", lambda g: 3 <= g <= 12)):
    hor, ver = neigh_pairs(sel)
    line = []
    for xl, X in (("ln T_m", np.log(Tm)), ("ln T_b", np.log(Tb)), ("IE", IE), ("EN", EN), ("r", RV)):
        dh = np.array([abs(X[j] - X[i]) for i, j in hor if np.isfinite(X[i]) and np.isfinite(X[j])])
        dv = np.array([abs(X[j] - X[i]) for i, j in ver if np.isfinite(X[i]) and np.isfinite(X[j])])
        p = stats.mannwhitneyu(dh, dv).pvalue
        line.append(f"{xl}: waagr. {dh.mean():.2f} / senkr. {dv.mean():.2f} (p = {fmt_p(p)})")
    P(f"    {blab}: " + "; ".join(line))
P()

# =====================================================================================
P("2) ELEKTRONEN")
P("-" * 90)
dIE = []
for i in range(len(Z) - 1):
    if np.isfinite(IE[i]) and np.isfinite(IE[i + 1]): dIE.append((np.log(IE[i + 1] / IE[i]), i))
dIE.sort()
nob = {2, 10, 18, 36, 54, 86}
top = [i for _, i in dIE[:6]]
P(f"2a) Größte Abfälle der Ionisierungsenergie Z → Z+1 (von {len(dIE)} Paaren): " + ", ".join(f"{sym[i]}→{sym[i+1]} ({IE[i]:.2f}→{IE[i+1]:.2f} eV)" for i in top))
hit = sum(Z[i] in nob for i in top)
P(f"    {hit} von 6 liegen am Edelgas → Alkalimetall; Zufallswahrscheinlichkeit, dass die 6 größten genau diese 6 Stellen treffen: 1/C({len(dIE)},6) = {1/math.comb(len(dIE),6):.1e}")
dips = []
for per in (2, 3, 4, 5, 6):
    def ie(g):
        k = [i for i in range(len(Z)) if PER[i] == per and GRP[i] == g]; return IE[k[0]] if k else np.nan
    dips.append((per, ie(2) > ie(13), ie(15) > ie(16)))
P("    Knicke innerhalb der Periode (Gruppe 2 > 13: neue Unterschale; 15 > 16: halbvolle p-Schale): " +
  ", ".join(f"P{p}: {'ja' if a else 'nein'}/{'ja' if b else 'nein'}" for p, a, b in dips))
P()

# IE als Längenverhältnis
rIE = E2 / 2 / IE * 1000   # pm
m = np.isfinite(IE) & np.isfinite(RV)
sl = np.polyfit(np.log(RV[m]), np.log(IE[m]), 1)[0]
bs = []
idx_m = np.where(m)[0]
for _ in range(NB):
    i = rng.choice(idx_m, len(idx_m)); bs.append(np.polyfit(np.log(RV[i]), np.log(IE[i]), 1)[0])
r, lo, hi, p, n = spear(RV, IE)
P("2b) Ionisierungsenergie als Längenverhältnis: IE/(m_e c²) = r_e/(2 r_IE), r_IE = e²/(8πε₀·IE); Wasserstoff: r_IE = a_B")
P(f"    r_IE: H {rIE[0]:.1f} pm (a_B = {aB*1e12:.1f}), He {rIE[1]:.1f}, Cs {rIE[S['Cs']]:.1f} pm; gemessener Radius / r_IE: Median {np.nanmedian(RV[m]/rIE[m]):.2f}, 16–84 % {np.nanpercentile(RV[m]/rIE[m],16):.2f}–{np.nanpercentile(RV[m]/rIE[m],84):.2f}")
prod = IE[m] * RV[m] / 1000 / E2
P(f"    IE·r_vdW / (e²/4πε₀): Median {np.median(prod):.2f}, 16–84 % {np.percentile(prod,16):.2f}–{np.percentile(prod,84):.2f}: Bindungsbruchteil IE/(m_e c²) ≈ r_e/r_vdW der Größenordnung nach")
P(f"    ln IE gegen ln r: Steigung {sl:.2f} [{np.percentile(bs,2.5):.2f}, {np.percentile(bs,97.5):.2f}] (Coulomb mit fester Ladung gäbe −1); Spearman ρ = {r:+.2f} [{lo:+.2f}, {hi:+.2f}], p = {fmt_p(p)}, n = {n}")
P(f"    IE/(m_e c²) liegt zwischen {np.nanmin(IE)/511e3:.1e} ({sym[np.nanargmin(IE)]}) und {np.nanmax(IE)/511e3:.1e} ({sym[np.nanargmax(IE)]}); Wasserstoff α²/2 = {(1/137.036)**2/2:.2e}")
P()

# Elektronenaffinität fehlt genau bei geschlossenen oder halbvollen Unterschalen?
import re as _re
def closed_or_half(conf):
    toks = _re.findall(r"(\d)([spdf])(\d+)", conf.split("(")[0])
    if not toks: return False
    full = {"s": 2, "p": 6, "d": 10, "f": 14}
    for _, l, o in toks:
        o = int(o)
        if not (o == full[l] or (l != "s" and o == full[l] // 2)): return False
    return True
lim = 56
CH = np.array([closed_or_half(CONF[i]) for i in range(lim)]); MISS = ~np.isfinite(EA[:lim])
k_in = int(np.sum(MISS & CH)); nm = int(MISS.sum()); nch = int(CH.sum())
p_h = stats.hypergeom.sf(k_in - 1, lim, nch, nm)
P(f"2c) Elektronenaffinität fehlt (Anion nicht oder kaum gebunden) für Z ≤ {lim}: " + ", ".join(sym[i] for i in range(lim) if MISS[i]))
P(f"    Valenzschalen alle voll oder p, d halbvoll (Definition [A]): " + ", ".join(sym[i] for i in range(lim) if CH[i]))
P(f"    {k_in} der {nm} fehlenden liegen in diesen {nch}; hypergeometrisch p = {p_h:.1e}. Ausnahmen mit positiver Affinität: " +
  ", ".join(f"{sym[i]} {EA[i]:.2f}" for i in range(lim) if CH[i] and not MISS[i]))
P()

# Herzfeld (Clausius-Mossotti): Metall, wenn q = (4π/3)·α'/V_A ≥ 1
def best_threshold(q, y):
    o = np.argsort(q); qs, ys = q[o], y[o]; best = (-1, None)
    for k in range(len(qs) + 1):
        acc = (np.sum(~ys[:k]) + np.sum(ys[k:])) / len(ys)
        thr = qs[0] - 1 if k == 0 else (qs[-1] + 1 if k == len(qs) else 0.5 * (qs[k - 1] + qs[k]))
        if acc > best[0]: best = (acc, thr)
    return best
aIE = E2 / 2 / IE * 10            # Å
q_IE = 6 * np.pi * aIE**3 / VA     # α' = (9/2)a³ (Wasserstoff-Form) [A]
q_r = (4 * np.pi / 3) * (RV / 100)**3 / VA   # α' = r_vdW³ [A]
P("2d) Kippstelle der Polarisation (Clausius-Mossotti, Herzfeld): ε = (1 + 2q)/(1 − q), q = (4π/3)·α'/V_A; bei q = 1 divergiert die Rückkopplung der lokalen Felder (Metall)")
P("    α' nicht in den Daten: geschätzt (A) α' = (9/2)·r_IE³ (Wasserstoff-Form) bzw. (B) α' = r_vdW³ [beides A, Faktor 2 unsicher]; V_A aus Dichte und Masse")
for lab, q in (("(A) aus IE", q_IE), ("(B) aus r_vdW", q_r)):
    sel = COND & np.isfinite(q) & (METAL | NONMET)
    y = METAL[sel]; qq = q[sel]
    acc, thr = best_threshold(np.log(qq), y); thr = math.exp(thr)
    acc1 = (np.sum(y & (qq >= 1)) + np.sum(~y & (qq < 1))) / len(y)
    bt = []
    ii = np.where(sel)[0]
    for _ in range(NB):
        j = rng.choice(ii, len(ii))
        if METAL[j].all() or (~METAL[j]).all(): continue
        bt.append(math.exp(best_threshold(np.log(q[j]), METAL[j])[1]))
    nullacc = [best_threshold(np.log(qq), rng.permutation(y))[0] for _ in range(1000)]
    pn = (np.sum(np.array(nullacc) >= acc) + 1) / 1001
    P(f"    {lab}: {y.sum()} Metalle, {(~y).sum()} Nichtmetalle; beste Schwelle q_c = {thr:.2f} [{np.percentile(bt,2.5):.2f}, {np.percentile(bt,97.5):.2f}], Trefferquote {acc:.2f} (Zufall p = {fmt_p(pn)}); bei q = 1: {acc1:.2f}")
    mis = [sym[i] for i in np.where(sel)[0] if (q[i] >= thr) != METAL[i]]
    P(f"       falsch eingeordnet: {', '.join(mis)};  Halbmetalle: " + ", ".join(f"{sym[i]} {q[i]:.2f}" for i in np.where(COND & METALLOID & np.isfinite(q))[0]))
for lab, X, sgn in (("Elektronegativität", EN, -1), ("Ionisierungsenergie", IE, -1)):
    sel = COND & np.isfinite(X) & (METAL | NONMET)
    acc, thr = best_threshold(sgn * X[sel], METAL[sel])
    P(f"    zum Vergleich {lab} allein: Trefferquote {acc:.2f} (Schwelle {abs(thr):.2f})")
P()

# =====================================================================================
P("3) KERN (AME2020)")
P("-" * 90)
L = open(os.path.join(DAT, "mass_1.mas20.txt"), encoding="latin-1").read().split("\n")
start = next(i for i, l in enumerate(L) if "micro-u" in l) + 1
nuc = {}
for l in L[start:]:
    if len(l) < 110: continue
    try:
        N_ = int(l[4:9]); Z_ = int(l[9:14]); A_ = int(l[14:19])
    except ValueError: continue
    mes = l[28:42].strip(); bas = l[54:67].strip()
    est = ("#" in mes) or ("#" in bas)
    ME = float(mes.replace("#", "."))
    B = Z_ * MEH + N_ * MEn - ME          # keV
    nuc[(Z_, N_)] = dict(A=A_, ME=ME, B=B, meas=not est, el=l[20:23].strip())
meas = {k: v for k, v in nuc.items() if v["meas"]}
P(f"3a) Kerne in AME2020: {len(nuc)}, davon gemessen {len(meas)}; schwerstes gemessenes Z = {max(k[0] for k in meas)}, größtes gemessenes A = {max(v['A'] for v in meas.values())}")
ba = sorted(((v["B"] / v["A"], k) for k, v in meas.items() if v["A"] > 1), reverse=True)
P("    größte Bindung je Nukleon B/A: " + ", ".join(f"{nuc[k]['A']}{nuc[k]['el']} {b:.2f} keV" for b, k in ba[:5]))
mpn = sorted(((v["A"] * u_keV + v["ME"] - k[0] * me * c**2 / eV / 1e3) / v["A"] - u_keV, k) for k, v in meas.items() if v["A"] > 4)
P("    kleinste Kernmasse je Nukleon (Atommasse − Z·m_e, relativ zu u): " + ", ".join(f"{nuc[k]['A']}{nuc[k]['el']} {x:.2f} keV" for x, k in mpn[:3]))
P("    Die kleinste Faltung je Nukleon, Δs/A = G·m/(3c²A), hat ⁵⁶Fe, nicht ⁶²Ni: ⁶²Ni hat mehr Neutronen, die schwerer sind als Protonen.")
P()

# Tröpfchenmodell
ks = [k for k, v in meas.items() if v["A"] >= 16]
Zs = np.array([k[0] for k in ks], float); Ns = np.array([k[1] for k in ks], float); As = Zs + Ns
Bs = np.array([meas[k]["B"] for k in ks]) / 1000.0
dP = np.where((Zs % 2 == 0) & (Ns % 2 == 0), 1.0, np.where((Zs % 2 == 1) & (Ns % 2 == 1), -1.0, 0.0))
X = np.column_stack([As, -As**(2/3), -Zs * (Zs - 1) / As**(1/3), -(Ns - Zs)**2 / As, dP / np.sqrt(As)])
coef, *_ = np.linalg.lstsq(X, Bs, rcond=None)
aV, aS, aC, aA, aP = coef; res = Bs - X @ coef
P(f"3b) Tröpfchenmodell (bekannt) an {len(ks)} gemessene Kerne, A ≥ 16: a_V = {aV:.2f}, a_S = {aS:.2f}, a_C = {aC:.3f}, a_A = {aA:.2f}, a_P = {aP:.2f} MeV; rms {np.sqrt(np.mean(res**2)):.2f} MeV")
def B_ldm(z, n):
    a = z + n; d = 1 if (z % 2 == 0 and n % 2 == 0) else (-1 if (z % 2 == 1 and n % 2 == 1) else 0)
    return aV * a - aS * a**(2/3) - aC * z * (z - 1) / a**(1/3) - aA * (n - z)**2 / a + aP * d / math.sqrt(a)
best = []
for a in range(20, 300):
    zz = min(range(1, a), key=lambda z: z * MEH / 1000 + (a - z) * MEn / 1000 - B_ldm(z, a - z) )
    best.append((B_ldm(zz, a - zz) / a, a, zz))
bmax = max(best)
xfun = lambda z, a: aC * z**2 / (2 * aS * a)
P(f"    Maximum von B/A im Tröpfchenmodell längs des Stabilitätstals: A = {bmax[1]}, Z = {bmax[2]}; gemessen ⁶²Ni (Z = 28)")
P(f"    Spaltbarkeit x = E_C/(2E_S) = a_C Z²/(2a_S A): am Maximum {xfun(bmax[2], bmax[1]):.3f}, ⁵⁶Fe {xfun(26,56):.3f}, ⁶²Ni {xfun(28,62):.3f}")
P("    hergeleitet: bei festem Z/A ist d(B/A)/dA = 0 genau bei E_C = E_S/2, also x = ¼; Oberfläche (je Nukleon fallend) gegen Coulomb (steigend).")
P("    Kippstelle: Die Verformungsenergie einer Quadrupolform ist (2/5)E_S·α²·(1 − x); bei x = 1 verschwindet die Rückstellung (Bohr-Wheeler), Antwort ∝ 1/(1 − x).")
heavy = sorted(meas.keys(), key=lambda k: (k[0], k[0] + k[1]))[-1]
xh = max(xfun(k[0], k[0] + k[1]) for k in meas)
kx = max(meas, key=lambda k: xfun(k[0], k[0] + k[1]))
P(f"    größtes x unter den gemessenen Kernen: {xh:.3f} ({kx[0]+kx[1]}{meas[kx]['el']}), 1/(1 − x) = {1/(1-xh):.1f}; schwerster gemessener {heavy[0]+heavy[1]}{meas[heavy]['el']}")
P()

# Magische Zahlen blind aus den Zwei-Nukleonen-Separationsenergien
def gaps(kind, sub):
    D = {}
    for (z, n), v in meas.items():
        if kind == "n":
            if n % 2: continue
            a, b, cc = (z, n - 2), (z, n), (z, n + 2)
        else:
            if z % 2: continue
            a, b, cc = (z - 2, n), (z, n), (z + 2, n)
        if a in meas and cc in meas and min(a) >= 1:
            s_here = meas[b]["B"] - meas[a]["B"]; s_next = meas[cc]["B"] - meas[b]["B"]
            d = (s_here - s_next) / 1000
            if sub:   # glatten Gang des Tröpfchenmodells abziehen
                d -= (B_ldm(*b) - B_ldm(*a)) - (B_ldm(*cc) - B_ldm(*b))
            D.setdefault(n if kind == "n" else z, []).append(d)
    return {k: (np.median(v), len(v)) for k, v in D.items() if len(v) >= 3}
magic = [2, 8, 20, 28, 50, 82, 126]
P("3c) Magische Zahlen blind: Schalenlücke D = S_2n(N) − S_2n(N+2) (gerades N), Median über alle gemessenen Z mit ≥ 3 Werten; ebenso aus S_2p")
P("    Hauptwert: Lücke minus Lücke des Tröpfchenmodells (entfernt den glatten Gang, der leichte Kerne groß macht)")
for sub in (True, False):
    for kind, lab in (("n", "N"), ("p", "Z")):
        D = gaps(kind, sub); keys = sorted(D)
        order = sorted(keys, key=lambda k: -D[k][0])
        mk = [k for k in magic if k in D]; ntop = len(mk)
        hits = len(set(order[:ntop]) & set(mk))
        p_hyp = stats.hypergeom.sf(hits - 1, len(keys), len(mk), ntop)
        tag = "ohne Tröpfchengang" if sub else "roh"
        P(f"    {lab} {tag}: größte " + ", ".join(f"{k} ({D[k][0]:.2f})" for k in order[:9]) + f"  [MeV; {len(keys)} Werte]")
        P(f"       unter den {ntop} größten {hits} magische (mit Daten: {', '.join(str(k) for k in mk)}); Zufall p = {p_hyp:.1e}; Rang der magischen: " +
          ", ".join(f"{k}:{order.index(k)+1}" for k in mk))
P()

# Tropflinien und Separationsenergien
def sep(kind):
    neg = []; tot = 0
    for (z, n), v in meas.items():
        pr = (z, n - 1) if kind == "n" else (z - 1, n)
        if pr in meas and pr[0] >= 0 and pr[1] >= 0:
            s = v["B"] - meas[pr]["B"]; tot += 1
            if s < 0: neg.append((z, n, s))
    return neg, tot
negn, totn = sep("n"); negp, totp = sep("p")
P(f"3d) Schwellen der Teilchenemission: gemessen S_n < 0 bei {len(negn)} von {totn}, S_p < 0 bei {len(negp)} von {totp} Kernen")
P("    S_n < 0: " + ", ".join(f"{z+n}{nuc[(z,n)]['el']}" for z, n, s in sorted(negn)[:12]) + (" …" if len(negn) > 12 else ""))
P("    S_p < 0 (größte Z): " + ", ".join(f"{z+n}{nuc[(z,n)]['el']} ({s:.0f} keV)" for z, n, s in sorted(negp, key=lambda t: -t[0])[:8]))
P("    Der Kern kann ein Nukleon abgeben, sobald die Separationsenergie das Vorzeichen wechselt; die Tropflinie ist diese Schwelle.")
P()

# Mehrfachstabilität bei geradem A (Paarung)
chains = {}
for (z, n), v in nuc.items():
    chains.setdefault(v["A"], {})[z] = v
mult = {"gerade": [], "ungerade": []}
for a, ch in chains.items():
    if a < 4: continue
    zs = sorted(ch); mins = []
    for z in zs:
        me_ = ch[z]["ME"]
        if (z - 1 not in ch or ch[z - 1]["ME"] > me_) and (z + 1 not in ch or ch[z + 1]["ME"] > me_) and ch[z]["meas"] and z - 1 in ch and z + 1 in ch:
            mins.append(z)
    if mins: mult["gerade" if a % 2 == 0 else "ungerade"].append(len(mins))
for k, v in mult.items():
    v = np.array(v); P(f"3e) A {k}: {len(v)} Isobarenketten; Zahl lokaler Massenminima (stabil gegen einfachen β-Zerfall): 1 → {np.sum(v==1)}, 2 → {np.sum(v==2)}, 3 → {np.sum(v==3)}")
P("    Bei geradem A liegen gerade-gerade und ungerade-ungerade Kerne auf zwei um 2·a_P/√A versetzten Parabeln: Mehrfachstabilität,")
P("    Übergang zwischen den Minima nur über doppelten β-Zerfall (Barriere = der Zwischenkern liegt höher).")
P()

# =====================================================================================
P("4) TRICHTERGRÖSSEN JE ELEMENT UND KERN (Modell: Δs = GM/(3c²), Takt 3Δs/r, α_G = 3Δs/λ̄ = (M/m_Pl)²)")
P("-" * 90)
P("    Element   Δs [m]      Δs/r_vdW    3Δs/R_Kern   α_G = (m/m_Pl)²   r* = √(Gm/a₀)   B/A/(m_u c²)  IE/(m_e c²)  kT_m/(m c²)")
Rk = lambda A: 1.2e-15 * A**(1/3)     # Kernradius [A]
BAmaxZ = {}
for (z, n), v in meas.items():
    if v["A"] > 1: BAmaxZ[z] = max(BAmaxZ.get(z, 0), v["B"] / v["A"])
for s in ("H", "C", "Fe", "Ni", "Pb", "U", "Og"):
    i = S[s]; m_ = M[i] * mu; ds = G * m_ / (3 * c**2); A = round(M[i])
    kt = kB * Tm[i] / (m_ * c**2) if np.isfinite(Tm[i]) else np.nan
    ba_ = BAmaxZ.get(Z[i], np.nan) / u_keV
    ie_ = IE[i] / 510998.95 if np.isfinite(IE[i]) else np.nan
    f1 = lambda v, f: ("–" if not np.isfinite(v) else format(v, f))
    P(f"    {s:3s}      {ds:.2e}   {f1(ds/(RV[i]*1e-12) if np.isfinite(RV[i]) else np.nan, '.1e'):7s}     {3*ds/Rk(A):.1e}      {(m_/mPl)**2:.1e}          {math.sqrt(G*m_/a0)*1e15:7.1f} fm      {f1(ba_, '.2e'):8s}     {f1(ie_, '.1e'):7s}      {f1(kt, '.1e')}")
fin = np.isfinite(Tm) & np.isfinite(RV)
ratio_g = (G * (M * mu)**2 / (2 * RV * 1e-12)) / (kB * Tm)
P(f"    Gravitation zweier Nachbaratome G m²/(2r_vdW) gegen k T_m: {np.nanmin(ratio_g[fin]):.1e} bis {np.nanmax(ratio_g[fin]):.1e}; gegen IE: "
  f"{np.nanmin((G*(M*mu)**2/(2*RV*1e-12))/(IE*eV)):.1e} bis {np.nanmax((G*(M*mu)**2/(2*RV*1e-12))/(IE*eV)):.1e}")
P(f"    Takt des Kerns an seiner Oberfläche gegen seine Bindung (3Δs/R)/(B/A/m_u c²): ⁵⁶Fe {3*G*56*mu/(3*c**2)/Rk(56)/(8790.36/u_keV):.1e}")
P(f"    Rückkehranteil im Labor (Erdfeld): q = e^(−√(g⊕/a₀)) = e^(−{math.sqrt(9.81/a0):.3g}); zwischen zwei Atomen allein wäre a_N ≈ {G*56*mu/(2e-10)**2:.1e} m/s² < a₀,")
P("    aber die Regel sieht das Gesamtfeld, und jedes Labor liegt tief im Newton-Bereich. Die Galaxienregel spielt für Elemente keine Rolle.")
P()

P("4b) Wo ein Trichter den Zustand seines Stoffs bestimmt: Druck im Zentrum mal Atomvolumen gegen kT")
P("    homogene Kugel: P_c·V_A = ½·m_A·GM/R = ½·m_A c²·δ,  δ = 3Δs/R = GM/(Rc²)  (Takt an der Oberfläche)")
P("    ⇒ δ_krit = 2kT/(m_A c²) = 2(v_th/c)²;  derselbe Quotient δ·m_A c²/(kT) entscheidet, ob ein Körper sein Gas hält (Jeans, bekannt)")
dk_m = 2 * kB * Tm / (M * mu * c**2); dk_b = 2 * kB * Tb / (M * mu * c**2)
def Mkrit(d, rho):
    return (d * c**2 / G)**1.5 * math.sqrt(3 / (4 * math.pi * rho))
okm = np.isfinite(dk_m) & COND
Ms = np.array([Mkrit(dk_m[i], RHO[i] * 1000) for i in range(len(Z)) if okm[i]]); Rs = (3 * Ms / (4 * np.pi * RHO[okm] * 1000))**(1/3)
P(f"    δ_krit(T_m) über {okm.sum()} feste/flüssige Elemente: {np.nanmin(dk_m[okm]):.1e} bis {np.nanmax(dk_m[okm]):.1e} (Median {np.nanmedian(dk_m[okm]):.1e})")
P(f"    Körpermasse dazu (Kugel aus dem Element): Median {np.median(Ms):.1e} kg, Radius Median {np.median(Rs)/1e3:.0f} km (16–84 %: {np.percentile(Rs,16)/1e3:.0f}–{np.percentile(Rs,84)/1e3:.0f} km)")
i = S["Fe"]; P(f"    Eisen: δ_krit = {dk_m[i]:.1e}, M = {Mkrit(dk_m[i], RHO[i]*1000):.1e} kg, R = {(3*Mkrit(dk_m[i], RHO[i]*1000)/(4*np.pi*RHO[i]*1000))**(1/3)/1e3:.0f} km")
dE = 3 * 1.4782e-3 / 6.371e6
P(f"    Erde δ = {dE:.2e} (Modell 01) = {dE/dk_m[i]:.0f}·δ_krit(Fe): Der Erdtrichter verschiebt die Phasengrenzen seines Eisens; ein Atom allein nie.")
P("    [G] Mond δ ≈ 3·10⁻¹¹ (M = 7,3·10²² kg, R = 1737 km); Kleinkörper werden erst ab einigen hundert km Radius rund (bekannte Grenze).")
P()

# =====================================================================================
P("5) REGELN DES MODELLS AN DEN DATEN")
P("-" * 90)
P("5a) Rückkopplungssumme 1/(1 − q), Kippstelle q = 1: drei Stellen an den Elementen mit derselben Form")
P(f"    Kern: q = x (Spaltbarkeit), Rückstellung ∝ 1 − x; gemessen x = {xfun(28,62):.2f} (⁶²Ni) bis {xh:.2f}; keine Kippstelle erreicht, schwerste Kerne halten nur über Schalen.")
P("    Hülle: q = (4π/3)α'/V_A (lokales Feld), ε = (1+2q)/(1−q); Metall ab q = 1 (Herzfeld) — geprüft in 2d.")
P("    [G] Magnetismus: χ = χ₀/(1 − I·N(E_F)) (Stoner); ferromagnetisch bei I·N(E_F) ≥ 1 (Fe, Co, Ni); keine Daten dazu.")
P("    Schmelzen und Sieden sind Übergänge erster Art: Der Gleichgewichtspunkt liegt zwischen zwei Stellen q = 1 (Unterkühlungs- und Überhitzungsgrenze);")
P("    die Tabellenwerte T_m, T_b sind nicht die Kippstellen, sondern der Gleichstand der freien Energien innerhalb des Hysteresefensters.")
P()
P(f"5b) Kopplungsfenster: fest = eingerastet (Gitterphasen gebunden), flüssig = gekoppelt, nicht eingerastet, Gas = entkoppelt.")
P(f"    Wären beide Ränder von einer Kopplung gesetzt, wäre T_b/T_m fest; gemessen {np.nanmin(ratio):.2f}–{np.nanmax(ratio):.2f}, nach Block verschieden (p = {fmt_p(pkw)}).")
nob_idx = [S[s] for s in ("Ne", "Ar", "Kr", "Xe", "Rn")]
P(f"    Edelgase (eine Paarkopplung): T_b/T_m = " + ", ".join(f"{sym[i]} {Tb[i]/Tm[i]:.3f}" for i in nob_idx) + "; Metalle breit: die Bindung überlebt den Verlust der Ordnung.")
P()
# geometrisches Mittel
P(f"5c) Geometrisches Mittel als Übergangslänge in der Hülle: λ̄_e² = r_e·a_B: {lam_e**2:.6e} = {re*aB:.6e} (exakt, α-Kette r_e = αλ̄_e = α²a_B)")
P(f"    Gegenstück im Modell ℓ_P² = 3Δs·λ̄; für ein Atom: 3Δs·λ̄ = {3*G*56*mu/(3*c**2)*hbar/(56*mu*c):.4e} = ℓ_P² = {lP**2:.4e} (Fe, für jede Masse gleich)")
Hs = math.sqrt(G * mu * 1.008 / a0); Og = math.sqrt(G * mu * 295 / a0)
P(f"    r* eines Atoms liegt zwischen Kern und Hülle ({Hs*1e15:.0f} fm bei H bis {Og*1e12:.2f} pm bei Og); ohne Wirkung, da jedes Labor im Newton-Bereich liegt.")
P()
orders = [n for n in range(1, 13) if abs(2 * math.cos(2 * math.pi / n) - round(2 * math.cos(2 * math.pi / n))) < 1e-9]
P(f"5d) Drehordnungen: Ein Kristall ist ein Zustand auf dem Torus ℝ³/Λ; eine Drehung des Gitters hat ganzzahlige Spur, 2cos(2π/n) ∈ ℤ ⇒ n ∈ {orders}")
P("    Das ist derselbe Satz wie für die flachen Tori (60°, 90°, 120°, 180°). Kristallstrukturen stehen nicht in den Daten; nicht elementweise prüfbar.")
P()
P("5e) Bose-Leiter, Planck-Form: in den Daten keine Wärmekapazität, kein Spektrum, keine Debye-Temperatur; nicht prüfbar.")
P("    Struktur [G]: Lindemann ⟨u²⟩ = (ħ/2mω)·coth(ħω/2kT) = (ħ/2mω)(1 + 2n_B); Schmelzen liegt im Rayleigh-Jeans-Bereich (kT ≫ ħω), die Form zählt beide Umlaufsinne.")
P()

# =====================================================================================
P("6) MUSTERSCAN ÜBER ALLE ELEMENTE MIT ZUFALLSSCHUTZ")
P("-" * 90)
BAZ = np.array([BAmaxZ.get(z, np.nan) for z in Z])
props = [("IE", IE), ("EN", EN), ("EA", EA), ("r", RV), ("V_A", VA), ("lnT_m", np.log(Tm)), ("lnT_b", np.log(Tb)), ("T_b/T_m", ratio), ("B/A_max", BAZ)]
pairs = list(combinations(range(len(props)), 2)); ntest = 2 * len(pairs)
res_raw = []; res_dif = []
for a, b in pairs:
    (na, xa), (nb_, xb) = props[a], props[b]
    r, lo, hi, p, n = spear(xa, xb, nb=300, npm=3000); res_raw.append((p, r, lo, hi, n, na, nb_))
    # Nachbardifferenzen in der Periode: entfernt den glatten Gang mit Z
    da = np.array([xa[i + 1] - xa[i] if PER[i] == PER[i + 1] else np.nan for i in range(len(Z) - 1)])
    db = np.array([xb[i + 1] - xb[i] if PER[i] == PER[i + 1] else np.nan for i in range(len(Z) - 1)])
    r, lo, hi, p, n = spear(da, db, nb=300, npm=3000); res_dif.append((p, r, lo, hi, n, na, nb_))
for lab, R_ in (("Werte", res_raw), ("Nachbardifferenzen", res_dif)):
    sig = [t for t in R_ if t[0] * ntest < 0.05]
    P(f"    {lab}: {len(sig)} von {len(R_)} Paaren nach Bonferroni (×{ntest}) signifikant:")
    for p, r, lo, hi, n, na, nb_ in sorted(sig, key=lambda t: -abs(t[1])):
        P(f"       {na:8s} ~ {nb_:8s} ρ = {r:+.2f} [{lo:+.2f}, {hi:+.2f}], n = {n}")
P("    Kreuzebene (B/A_max gegen Hülle und Phasen) in den Nachbardifferenzen: " +
  ", ".join(f"{na}~{nb_} ρ = {r:+.2f} (p = {fmt_p(p)})" for p, r, lo, hi, n, na, nb_ in res_dif if "B/A_max" in (na, nb_)))
P()
# Kernmagische Z und gerade/ungerade Z in Hülle und Phasen?
def local_res(X):
    return np.array([X[i] - 0.5 * (X[i - 1] + X[i + 1]) if 0 < i < len(Z) - 1 else np.nan for i in range(len(Z))])
P("    Kernschalen in Hülle und Phasen? Lokalresiduum X(Z) − ½[X(Z−1) + X(Z+1)] an Z = 8, 20, 28, 50, 82 gegen 5 zufällige Z (Permutation)")
for xl, X in (("ln T_m", np.log(Tm)), ("ln T_b", np.log(Tb)), ("IE", IE), ("B/A_max", BAZ)):
    lr = local_res(X); good = np.where(np.isfinite(lr))[0]
    mz = [S[s] for s in ("O", "Ca", "Ni", "Sn", "Pb")]
    obs = np.nanmean(lr[mz]); null = np.array([np.mean(lr[rng.choice(good, 5, replace=False)]) for _ in range(NP)])
    p = (np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (NP + 1)
    ev = [i for i in good if Z[i] % 2 == 0]; od = [i for i in good if Z[i] % 2 == 1]
    pe = stats.mannwhitneyu(lr[ev], lr[od]).pvalue
    P(f"       {xl:8s}: magisch {obs:+.3f}, Zufall {null.mean():+.3f} ± {null.std():.3f}, p = {fmt_p(p)};  gerades gegen ungerades Z: {np.mean(lr[ev]):+.3f} / {np.mean(lr[od]):+.3f}, p = {fmt_p(pe)}")
P()
# Gerade-ungerade in IE: nur Schalenschluss oder auch sonst?
P("    Gerade/ungerade Z in der Ionisierungsenergie ohne Edelgase und ihre Nachbarn (Z_Edelgas − 1 … + 1):")
nobZ = {2, 10, 18, 36, 54, 86}; near = {z + d for z in nobZ for d in (-1, 0, 1)}
lr = local_res(IE); lrB = local_res(BAZ)
good = [i for i in range(len(Z)) if np.isfinite(lr[i]) and Z[i] not in near]
ev = [i for i in good if Z[i] % 2 == 0]; od = [i for i in good if Z[i] % 2 == 1]
P(f"       IE: gerade {np.mean(lr[ev]):+.3f} / ungerade {np.mean(lr[od]):+.3f} eV, p = {fmt_p(stats.mannwhitneyu(lr[ev], lr[od]).pvalue)} (n = {len(ev)}/{len(od)})")
for blk, sel in (("s- und p-Block", lambda i: GRP[i] in (1, 2) or GRP[i] >= 13), ("d-Block", lambda i: 3 <= GRP[i] <= 12), ("f-Block", lambda i: GRP[i] == 0)):
    e2 = [i for i in ev if sel(i)]; o2 = [i for i in od if sel(i)]
    if len(e2) >= 3 and len(o2) >= 3:
        P(f"          {blk}: gerade {np.mean(lr[e2]):+.3f} / ungerade {np.mean(lr[o2]):+.3f}, p = {fmt_p(stats.mannwhitneyu(lr[e2], lr[o2]).pvalue)} (n = {len(e2)}/{len(o2)})")
m2 = np.isfinite(lr) & np.isfinite(lrB) & np.array([z not in near for z in Z])
r, lo, hi, p, n = spear(lr[m2], lrB[m2])
P(f"       Lokalresiduen IE gegen B/A_max ohne Edelgasnähe: ρ = {r:+.2f} [{lo:+.2f}, {hi:+.2f}], p = {fmt_p(p)}, n = {n}")
P()

txt = "\n".join(out)
open(os.path.join(HERE, "elemente_analyse.txt"), "w", encoding="utf-8").write(txt); print(txt)
