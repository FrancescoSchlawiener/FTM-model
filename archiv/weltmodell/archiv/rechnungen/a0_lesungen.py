# -*- coding: utf-8 -*-
"""
a₀-Lesungen an SPARC: warum dieselben Daten 0,95, 1,1 und 1,53·10⁻¹⁰ geben, und eine konsistente Lesung.
Alle Formeln aus dem Kern (ftm_kern: a_regel, nu, s_, chi); hier nur Datenlesen und Anpassung.
Blickrichtung: a₀ = c²/2Λ ist eine Größe der inneren Richtung (I), seit z ≈ 10 fest. Die Lesung an Galaxien ist L:
gelesen wird die am Ort wirksame Regel, verschoben durch Υ, Abstand, Neigung, Ort auf der Kurve und äußeres Feld.
  1a Tully-Fisher mit V_flach: Zerlegung log(V_f⁴/(G M_b a₀)) = 4 log(V_f/V_Regel) + 2 log(ν s) + log(M_enc/M_b)
  1b Log-kleinste-Quadrate: Schritte von 0,95 zur profilierten Anpassung
  2  Konsistente Lesung: profiliert (Abstand, Neigung, Υ_d), Bootstrap, Systematik
  3  Handwert
Aufruf: python a0_lesungen.py > a0_lesungen.txt
"""
import io, os, sys, math
import numpy as np
from scipy.optimize import least_squares, brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
import ftm_kern as K
from ftm_kern import G, M_sun, kpc, A0
import sparc_universal as su          # nur Datenlesen (Tabelle)

D0 = os.path.join(HIER, "..", "daten", "sparc", "Rotmod_LTG")
tv = su.tab_voll()
GAL = {}
for nm, t in tv.items():
    fn = os.path.join(D0, nm + "_rotmod.dat")
    if not os.path.exists(fn):
        continue
    a = np.loadtxt(fn, comments="#")
    GAL[nm] = dict(t=t, R=a[:, 0] * kpc, Vo=a[:, 1] * 1e3, eV=a[:, 2] * 1e3, Vg=a[:, 3] * 1e3, Vd=a[:, 4] * 1e3,
                   Vb=a[:, 5] * 1e3, SBb=a[:, 7], Rk=a[:, 0])


def vbar2(g, Yd=0.5, Yb=0.7):
    return g["Vg"] * np.abs(g["Vg"]) + Yd * g["Vd"]**2 + Yb * g["Vb"]**2


def auswahl(t):
    return t["Q"] <= 2 and t["inc"] >= 30


def linie(t):
    print(); print("=" * 110); print(t); print("=" * 110)


print(f"a₀-Lesungen an SPARC, Regel a = ν(a_N)·a_N aus ftm_kern, Referenz a₀ = {A0:.2e} m/s²; Υ_d = 0,5, Υ_b = 0,7 wo nicht anders gesagt")

# ============================================================================ 1a Tully-Fisher
linie("1a Tully-Fisher mit V_flach (Q ≤ 2, Neigung ≥ 30°, V_flach > 0): wie tief sind die äußersten Punkte, was fehlt in M_b?")
rows = []
for nm, g in GAL.items():
    t = g["t"]
    if not (auswahl(t) and t["Vf"] > 0):
        continue
    Vb2 = vbar2(g)
    if not (Vb2[-3:] > 0).any():
        continue
    k = len(g["R"]) - 1
    while Vb2[k] <= 0:
        k -= 1
    R, aN = g["R"][k], Vb2[k] / g["R"][k]
    Mb = (0.5 * t["L"] + 1.33 * t["MHI"]) * 1e9 * M_sun
    # Bulge-Leuchtkraft aus dem Profil SBbul [L☉/pc²] (Trapez in R): in L steckt sie mit Υ = 0,5 statt 0,7
    Lb = np.trapezoid(2 * math.pi * g["Rk"] * 1e3 * g["SBb"], g["Rk"] * 1e3) if g["SBb"].any() else 0.0
    rows.append(dict(nm=nm, Vf=t["Vf"] * 1e3, Mb=Mb, R=R, aN=aN, Menc=Vb2[k] * R / G, Lb=Lb * M_sun, L=t["L"] * 1e9 * M_sun))
n = len(rows)
Vf = np.array([r["Vf"] for r in rows]); Mb = np.array([r["Mb"] for r in rows]); R = np.array([r["R"] for r in rows])
aN = np.array([r["aN"] for r in rows]); Menc = np.array([r["Menc"] for r in rows])
Lb = np.array([r["Lb"] for r in rows]); L = np.array([r["L"] for r in rows])
s = K.s_(aN, A0)
Vr = np.sqrt(K.a_regel(aN, A0) * R)
print(f"   {n} Galaxien; äußerster Punkt mit V_bar² > 0")
q = np.percentile(s, [10, 25, 50, 75, 90])
print(f"   s = √(a_N/a₀) am äußersten Punkt: 10/25/50/75/90 % = " + " / ".join(f"{x:.3f}" for x in q) +
      f";  a_N/a₀ Median {np.median(aN / A0):.3f};  Anteil s > 0,3: {np.mean(s > 0.3):.0%}, s > 0,5: {np.mean(s > 0.5):.0%}")
nus = K.nu(aN, A0) * s
lg_tf = np.log10(Vf**4 / (G * Mb * A0))
tA = 4 * np.log10(Vf / Vr); tB = 2 * np.log10(nus); tC = np.log10(Menc / Mb)
print(f"   Zerlegung (Median, dex) log(V_f⁴/(G M_b a₀)) = 4 log(V_f/V_Regel) + 2 log(ν s) + log(M_enc/M_b) [exakt je Galaxie]:")
print(f"      gesamt {np.median(lg_tf):+.3f}  (a₀ aus V⁴/GM_b: {A0 * 10**np.median(lg_tf) * 1e10:.3f}·10⁻¹⁰)")
print(f"      Daten gegen Regel am Punkt 4 log(V_f/V_Regel) {np.median(tA):+.3f} (Mittel {tA.mean():+.3f})")
print(f"      endliche Tiefe 2 log(ν s) = 2 log(s/(1 − e^{{−s}})) {np.median(tB):+.3f} (Mittel {tB.mean():+.3f}; ν s → 1 nur für s → 0)")
print(f"      Masse log(V_bar²R/(G M_b)): Kugel-Masse am Punkt gegen Tabellenmasse {np.median(tC):+.3f} (Mittel {tC.mean():+.3f})")
print(f"      Summe der Mittel {tA.mean() + tB.mean() + tC.mean():+.3f} = Mittel gesamt {lg_tf.mean():+.3f}")
fb = 0.2 * Lb / Mb
print(f"   Bulge: {np.sum(Lb > 0)} von {n} mit Bulge; M_b = 0,5·L + 1,33·M_HI zählt ihn mit 0,5 statt 0,7: M_b zu klein um "
      f"Median {np.median(fb[Lb > 0]):.1%} (größter {fb.max():.1%}); a₀_TF sinkt dadurch im Median um "
      f"{(A0 * 10**np.median(lg_tf) - A0 * 10**np.median(lg_tf - np.log10(1 + fb))) * 1e10:.3f}·10⁻¹⁰")
# konsistente Tully-Fisher-Lesung: a₀, bei dem die Regel am äußersten Punkt V_flach trifft (Median)
f_tf = lambda la: np.median(np.log10(Vf / np.sqrt(K.a_regel(aN, 10**la) * R)))
la_tf = brentq(f_tf, -10.6, -9.4)
f_tf_m = lambda la: np.mean(np.log10(Vf / np.sqrt(K.a_regel(aN, 10**la) * R)))
la_tfm = brentq(f_tf_m, -10.6, -9.4)
print(f"   Lesung mit der Regel selbst (V_flach = V_Regel am äußersten Punkt): Median {10**la_tf * 1e10:.3f}, Mittel {10**la_tfm * 1e10:.3f} ·10⁻¹⁰")
s_tf = K.s_(aN, 10**la_tf)
print(f"   Tiefgrenze allein verlangt s → 0: nur Galaxien mit s < 0,2 ({np.sum(s < 0.2)}): V⁴/(G M_b) Median "
      f"{np.median(Vf[s < 0.2]**4 / (G * Mb[s < 0.2])) * 1e10:.3f}·10⁻¹⁰, dort 2 log(ν s) = {np.median(tB[s < 0.2]):+.3f}, "
      f"log(M_enc/M_b) = {np.median(tC[s < 0.2]):+.3f}")

# ============================================================================ 1b Log-kleinste-Quadrate
linie("1b Log-kleinste-Quadrate: von der ungewichteten Anpassung zur profilierten (Abstand, Neigung, Υ_d)")
AG = A0 * 10**np.linspace(-0.25, 0.20, 46)


def bestwert(c2, ag=AG):
    k = int(np.argmin(c2)); k = min(max(k, 1), len(ag) - 2)
    x = np.log(ag); y0, y1, y2 = c2[k - 1:k + 2]; h = x[1] - x[0]
    return math.exp(x[k] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2)))


def punkte(sel, Yd=0.5, nur_eV=False):
    out = []
    for nm, g in GAL.items():
        t = g["t"]
        if sel and not auswahl(t):
            continue
        Vb2 = vbar2(g, Yd, 1.4 * Yd)
        m = (Vb2 > 0) & (g["Vo"] > 0) & (g["R"] > 0)
        if nur_eV:
            m &= g["eV"] > 0
        if sel and m.sum() < 5:
            continue
        if m.sum() == 0:
            continue
        out.append(dict(nm=nm, gN=Vb2[m] / g["R"][m], go=g["Vo"][m]**2 / g["R"][m], w=(2 * g["eV"][m] / g["Vo"][m])))
    return out


def ls_frei(pp, gew="gleich", sig=0.1):
    """Σ w·(log g_obs − log a_Regel)² über a₀ ohne Störgrößen."""
    c2 = np.zeros(len(AG))
    for p in pp:
        if gew == "gleich":
            w = np.ones_like(p["gN"])
        elif gew == "Fehler":
            w = 1 / ((p["w"] / math.log(10))**2 + (sig / math.log(10))**2)
        elif gew == "je Galaxie":
            w = np.ones_like(p["gN"]) / p["gN"].size
        for j, a0 in enumerate(AG):
            c2[j] += np.sum(w * np.log10(p["go"] / K.a_regel(p["gN"], a0))**2)
    return bestwert(c2)


def median_lesung(pp, nur_tief=True):
    gN = np.concatenate([p["gN"] for p in pp]); go = np.concatenate([p["go"] for p in pp])
    f = lambda la: np.median(np.log10(go / K.a_regel(gN, 10**la))[(gN < 10**la) if nur_tief else slice(None)])
    return 10**brentq(f, -10.6, -9.4)


P_all, P_sel, P_selE = punkte(False), punkte(True), punkte(True, nur_eV=True)
gN = np.concatenate([p["gN"] for p in P_sel]); go = np.concatenate([p["go"] for p in P_sel])
r = np.log10(go / K.a_regel(gN, A0)); ch = K.chi(gN, A0)
print(f"   Hebel: d ln a/d ln a₀ = χ/2 (Kern: 1 − Steigung); die Anpassung setzt Σ χ·Rest = 0, also das χ-gewichtete Mittel der Reste")
print(f"   Reste bei a₀ = 1,1 (Auswahl, {gN.size} Punkte): Median {np.median(r):+.4f}, Mittel {r.mean():+.4f}, χ-gewichtetes Mittel "
      f"{np.sum(ch * r) / ch.sum():+.4f} dex; Schiefe: Anteil r < −0,3 dex {np.mean(r < -0.3):.1%}, r > +0,3 dex {np.mean(r > 0.3):.1%}")
lo = r < -0.3
print(f"      Punkte mit r < −0,3 dex: {lo.sum()}, davon mit a_N > a₀ {np.mean(gN[lo] > A0):.0%}")
schritte = [
    ("alle 175 Kurven, gleich gewichtet (Gesamtlauf)", ls_frei(P_all)),
    ("Auswahl Q ≤ 2, Neigung ≥ 30°, gleich gewichtet (Gesamtlauf)", ls_frei(P_sel)),
    ("  + nur Punkte mit errV > 0", ls_frei(P_selE)),
    ("  + jede Galaxie gleich (1/N_Punkte)", ls_frei(P_selE, "je Galaxie")),
    ("  + Gewicht 1/(σ_V² + σ_int²), σ_int = 0,10 ln", ls_frei(P_selE, "Fehler", 0.10)),
    ("  + Gewicht 1/(σ_V² + σ_int²), σ_int = 0,20 ln", ls_frei(P_selE, "Fehler", 0.20)),
    ("Median der Reste = 0 (Punkte mit a_N < a₀), Auswahl", median_lesung(P_sel)),
    ("Median der Reste = 0 (alle Punkte), Auswahl", median_lesung(P_sel, False)),
]


# profiliert
def gal_liste(Yd0):
    out = []
    for nm, g in GAL.items():
        t = g["t"]
        if not auswahl(t):
            continue
        m = (g["Vo"] > 0) & (g["R"] > 0) & (g["eV"] > 0)
        if m.sum() < 5:
            continue
        out.append(dict(nm=nm, fD=t["fD"], R=g["R"][m], Vo=g["Vo"][m], eV=g["eV"][m], Vg=g["Vg"][m], Vd=g["Vd"][m], Vb=g["Vb"][m],
                        eD=max(t["eD"] / t["D"], 0.01), inc=math.radians(t["inc"]), ei=math.radians(max(t["eInc"], 1.0))))
    return out


GL = gal_liste(0.5)
sY = math.log(10) * 0.1
H_VERSCH = 0.0      # ln(73/H₀): Hubble-Fluss-Abstände der Tabelle mit H₀ = 73 auf das H₀ des Kerns umgerechnet (Variante); g_obs ∝ 1/D, a_N frei von D


def resid(par, g, a0, Y0, sig, profil):
    uD, di, uY = par if profil else (0.0, 0.0, 0.0)
    Yd = Y0 * math.exp(uY)
    gb = np.maximum((g["Vg"] * np.abs(g["Vg"]) + Yd * g["Vd"]**2 + 1.4 * Yd * g["Vb"]**2) / g["R"], 1e-14)
    inc = min(max(g["inc"] + di, 0.2), math.pi / 2)
    lgo = np.log(g["Vo"]**2 / g["R"]) - uD - (H_VERSCH if g["fD"] == 1 else 0.0) + 2 * math.log(math.sin(g["inc"]) / math.sin(inc))
    rr = (lgo - np.log(K.a_regel(gb, a0))) / np.sqrt((2 * g["eV"] / g["Vo"])**2 + sig**2)
    return np.r_[rr, uD / g["eD"], di / g["ei"], uY / sY] if profil else rr


def tafel(Y0, sig, profil=True, loss="linear"):
    T = np.zeros((len(GL), len(AG)))
    for i, g in enumerate(GL):
        x = np.zeros(3)
        for j, a0 in enumerate(AG):
            if profil:
                sol = least_squares(resid, x, args=(g, a0, Y0, sig, True), loss=loss, f_scale=1.0)
                x = sol.x; T[i, j] = 2 * sol.cost
            else:
                rr = resid(None, g, a0, Y0, sig, False)
                T[i, j] = np.sum(rr**2) if loss == "linear" else np.sum(2 * (np.sqrt(1 + rr**2) - 1))
    return T


rng = np.random.default_rng(7)


def boot(T, idx=None):
    idx = np.arange(T.shape[0]) if idx is None else idx
    b = bestwert(T[idx].sum(0))
    bs = np.array([bestwert(T[idx[rng.integers(0, idx.size, idx.size)]].sum(0)) for _ in range(1000)])
    return b, bs.std()


npk = sum(g["R"].size for g in GL)
T_np = tafel(0.5, 0.10, profil=False)
b_np, e_np = boot(T_np)
schritte.append(("  gleiche Punkte, Fehlergewicht σ_int = 0,10, ohne Profil (Υ 0,5)", b_np))
TAF = {}
for Y0 in (0.4, 0.5, 0.6):
    for sig in (0.10, 0.20):
        TAF[(Y0, sig, "linear")] = tafel(Y0, sig)
TAF[(0.5, 0.10, "soft_l1")] = tafel(0.5, 0.10, loss="soft_l1")
H_VERSCH = math.log(73.0 / K.H0_KMS)
TAF[(0.5, 0.10, "H₀ 67,4")] = tafel(0.5, 0.10)
H_VERSCH = 0.0
schritte.append(("  + Abstand, Neigung, Υ_d profiliert (Υ 0,5, σ_int 0,10)", bestwert(TAF[(0.5, 0.10, 'linear')].sum(0))))
schritte.append(("  + robust (soft_l1) statt Quadrat", bestwert(TAF[(0.5, 0.10, 'soft_l1')].sum(0))))
for lab, v in schritte:
    print(f"   {lab:70s} a₀ = {v * 1e10:.3f}·10⁻¹⁰")
print(f"   ({len(GL)} Galaxien, {npk} Punkte in den gewichteten und profilierten Zeilen)")

# ============================================================================ 2 Konsistente Lesung
linie("2  Konsistente Lesung: profiliert, Bootstrap über Galaxien (1000), Systematik aus Υ_d, σ_int, Verlust, Abstandsmethode")
print("   Υ_d | σ_int | Verlust  | a₀ [10⁻¹⁰] | Bootstrap | χ²/Punkt | Primär-Abstände | Hubble-Fluss")
prim = np.array([i for i, g in enumerate(GL) if g["fD"] != 1]); hub = np.array([i for i, g in enumerate(GL) if g["fD"] == 1])
res = {}
for key, T in TAF.items():
    b, e = boot(T); bp, ep = boot(T, prim); bh, eh = boot(T, hub)
    res[key] = (b, e)
    print(f"   {key[0]:.1f} | {key[1]:.2f}  | {key[2]:8s} | {b * 1e10:10.3f} | {e * 1e10:9.3f} | {T.sum(0).min() / npk:8.2f} | "
          f"{bp * 1e10:.3f} ± {ep * 1e10:.3f} ({prim.size}) | {bh * 1e10:.3f} ± {eh * 1e10:.3f} ({hub.size})")
b0, e0 = res[(0.5, 0.10, "linear")]
lin = [v[0] for k, v in res.items() if k[2] == "linear"]
sysH = abs(res[(0.5, 0.10, "H₀ 67,4")][0] - b0)
sysY = max(abs(res[(Y, 0.10, 'linear')][0] - b0) for Y in (0.4, 0.6))
sysS = abs(res[(0.5, 0.20, 'linear')][0] - b0)
sysL = abs(res[(0.5, 0.10, 'soft_l1')][0] - b0)
sysT = math.sqrt(sysY**2 + sysS**2 + sysL**2)
print(f"   Hubble-Fluss-Abstände mit H₀ = {K.H0_KMS} statt 73 (D × {73 / K.H0_KMS:.3f}): a₀ = {res[(0.5, 0.10, 'H₀ 67,4')][0] * 1e10:.3f}, "
      f"Verschiebung {(res[(0.5, 0.10, 'H₀ 67,4')][0] - b0) * 1e10:+.3f} (in der Systematik nicht enthalten, Wahl des H₀ der Abstände)")
print(f"   Lesung: a₀ = {b0 * 1e10:.2f} ± {e0 * 1e10:.2f} (stat.) ± {sysT * 1e10:.2f} (syst.: Υ_d ±0,1 {sysY * 1e10:.2f}, σ_int {sysS * 1e10:.2f}, "
      f"Verlust {sysL * 1e10:.2f}) ·10⁻¹⁰;  Spanne aller Varianten {min(lin) * 1e10:.2f}–{max(lin) * 1e10:.2f}")
print(f"   Tully-Fisher mit der Regel (1a): {10**la_tf * 1e10:.2f}·10⁻¹⁰ (Median), dieselbe Größe von außen gelesen")
print(f"   Zum Vergleich: cH₀/2π = {K.a0_modell() * 1e10:.3f} (H₀ = 67,4), {K.a0_modell(73.0) * 1e10:.3f} (H₀ = 73); Zelle Würfel ohne Kreispaare ≤ 1,128")
print("   Blickrichtung: gelesen ist die am Ort wirksame Regel (L). Ein äußeres Feld senkt sie unter die Regel ohne Feld;")
print("   die Lesung ist daher eine untere Schranke der inneren Größe a₀ = c²/2Λ (I), nicht a₀ selbst.")

# ============================================================================ 3 Handwert
linie("3  Handwert")
sh = 0.3
print(f"   Punktmasse, s = 0,3: e^{{−0,3}} = {math.exp(-0.3):.5f}, 1 − e^{{−s}} = {1 - math.exp(-0.3):.5f}, ν s = {0.3 / (1 - math.exp(-0.3)):.5f}, "
      f"(ν s)² = {(0.3 / (1 - math.exp(-0.3)))**2:.4f}: V⁴/(GMa₀) = (ν s)² → Tully-Fisher liest a₀ um {((0.3 / (1 - math.exp(-0.3)))**2 - 1):.0%} zu hoch")
print(f"   Kern: ν(0,09 a₀)·0,3 = {float(K.nu(0.09 * A0)) * 0.3:.5f}")
