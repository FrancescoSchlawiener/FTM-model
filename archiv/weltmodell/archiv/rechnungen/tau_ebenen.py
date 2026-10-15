# -*- coding: utf-8 -*-
"""
τ je Ebene und Spiegelung an der Schwelle.
1  Ist r* (s = 1) der selbstduale Punkt einer echten Symmetrie der Regel (s → 1/s, a_N → a₀²/a_N)?
   Test an ν(s), χ(s), a(a_N) und an der gemessenen Relation (SPARC).
2  τ = Re + i·Im je Ebene (Teilchen, Atom, Planet, Stern, Planetenbahn, Galaxie, Haufen, Zelle),
   beide Schreibweisen des Imaginärteils; Vergleich unter SL(2,ℤ) (Reduktion auf den Fundamentalbereich).
3  α und α_G als Betrag bzw. Phase eines τ: legt die Struktur eine Beziehung fest?
Konstanten und Formeln aus ftm_kern.
"""
import io, sys, os, glob, math, cmath
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import (G, c, hbar, m_p, m_e, alpha, r_e, R_p, a_B, M_sun, R_sun, M_earth, R_earth, AU, pc,
                      A0, ds, r_stern, alpha_G, compton_red, schliesslaenge, n_B)

kpc = 1e3 * pc
def nu_s(s):  return 1.0 / (-np.expm1(-np.asarray(s, float)))
def chi_s(s): s = np.asarray(s, float); return s / np.expm1(s)
def x_of_y(y): return np.asarray(y) * nu_s(np.sqrt(y))          # a/a₀ als Funktion von a_N/a₀
def re_tau(s): return math.sqrt(1 + float(chi_s(s)))             # Punktmasse: (Re τ)² = 1 + χ

L = lambda t: print(t)
def kopf(t): L("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

# ================================================================================================
kopf("1  Spiegelung an der Schwelle: s → 1/s  (Punktmasse: s = r*/r, also r → r*²/r;  a_N → a₀²/a_N)")
# ------------------------------------------------------------------------------------------------
Ld = 3 * ds(M_sun); rs = r_stern(M_sun); l0 = c**2 / A0
L(f"1a Leiter Sonne: 3Δs = {Ld:.4e} m, r* = {rs:.4e} m, ℓ₀ = c²/a₀ = {l0:.4e} m;  r*²/(3Δs·ℓ₀) = {rs**2/(Ld*l0):.12f}")
L(f"   a_N(3Δs)·a_N(ℓ₀)/a₀² = {(G*M_sun/Ld**2)*(G*M_sun/l0**2)/A0**2:.12f}   (Spiegelung r → r*²/r ist a_N → a₀²/a_N)")
L(f"   Bild des innersten Ortes r₀ = 6Δs: r*²/r₀ = {rs**2/(6*ds(M_sun)):.4e} m = Λ = c²/2a₀ = {schliesslaenge():.4e} m")

L("\n1b Gesuchte Dualität: s^k·ν(s) = s^{−k}·ν(1/s) für ein festes k (Gewicht).  Führende Ordnung beider Ränder:")
L("   s → 0: ν ≈ 1/s;  s → ∞: ν → 1.  Beide Ränder zugleich nur mit k = 1/2, d. h. ν(1/s) = s·ν(s),")
L("   in Beschleunigungen a(a₀²/a_N)·(a_N/a₀)^{3/2} = a(a_N)  (tiefer Ast √(a_N a₀) ↔ Newton-Ast a_N).")
L(f"   {'s':>7s} {'ν(s)':>10s} {'ν(1/s)':>10s} {'s·ν(s)':>10s} {'ν(1/s)/(s ν(s))':>16s} {'[dex]':>8s}  {'χ(s)':>9s} {'χ(1/s)':>11s}")
for s in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
    a, b = float(nu_s(s)), float(nu_s(1/s))
    L(f"   {s:7.3f} {a:10.5f} {b:10.5f} {s*a:10.5f} {b/(s*a):16.6f} {math.log10(b/(s*a)):+8.4f}  {float(chi_s(s)):9.5f} {float(chi_s(1/s)):11.3e}")
ss = np.exp(np.linspace(math.log(1e-3), 0, 4001))
best = None
for k in np.linspace(0.3, 0.7, 401):
    D = np.abs(k*np.log(ss) + np.log(nu_s(ss)) - (-k*np.log(ss) + np.log(nu_s(1/ss))))
    if best is None or D.max() < best[1]: best = (k, D.max(), ss[np.argmax(D)])
L(f"   bestes festes k (kleinste größte Abweichung über s = 10⁻³…1): k = {best[0]:.3f}, max |Δ ln| = {best[1]:.4f} "
  f"({best[1]/math.log(10):.4f} dex) bei s = {best[2]:.3f}")
L("   Reihen: tiefer Rand ν(s) = 1/s + 1/2 + s/12 − s³/720 … (Potenzen, Bernoulli);")
L("           Newton-Rand ν(1/s) = 1 + e^{−1/s} + e^{−2/s} … (exponentiell).  Eine Seite Potenzreihe, die andere")
L("           exponentiell klein: keine Funktionalgleichung f(1/s) = s^{2k}·f(s) kann beide verbinden.")
L("   χ(s) → 1 (s → 0), χ(1/s) → e^{−1/s}/s: kein Gewicht k gleicht konstant gegen exponentiell aus.")

L("\n1c Exakte Symmetrien der Regelfunktionen (Kontrolle):")
for s in [0.3, 1.0, 2.5]:
    L(f"   s = {s}: ν(s) + ν(−s) = {float(nu_s(s) + nu_s(-s)):.12f};   χ(−s) − χ(s) − s = {float(chi_s(-s) - chi_s(s) - s):+.2e}")
L("   Exakt ist nur die Spiegelung s → −s (|q| → 1/|q|, Spiegelung am Ort |q| = 1), nicht s → 1/s.")

L("\n1d Modulare Spiegelung τ → −1/τ in der Schreibweise Im τ = s/2π (Punkt 1): s → 4π²/s, Fixpunkt s = 2π (r = r*/2π).")
def ln_eta_s(s, N=4000):   # τ = i s/2π rein imaginär, q = e^{−s}
    n = np.arange(1, N + 1); return -s / 24 + np.sum(np.log1p(-np.exp(-n * s)))
for s in [0.7, 1.0, 2.0, 2*math.pi]:
    sd = 4*math.pi**2/s
    lhs = ln_eta_s(sd); rhs = ln_eta_s(s) + 0.5*math.log(s/(2*math.pi))
    L(f"   η: s = {s:7.4f} ↔ 4π²/s = {sd:8.4f}:  ln η(−1/τ) − [ln η(τ) + ½ ln(−iτ)] = {lhs - rhs:+.2e}")
L("   Die Produktform über alle Oberschwingungen ∏(1 − qⁿ) ist exakt modular (Fixpunkt s = 2π); die Regel ist ein")
L("   einzelner Faktor 1/(1 − q) und zählt nur die Grundschwingung. Schwellenpunkt s = 1 und modularer Fixpunkt s = 2π")
L(f"   fallen nicht zusammen (Faktor 2π); bei s = 2π ist Re τ = {re_tau(2*math.pi):.5f}, τ ≡ i + {re_tau(2*math.pi)-1:.4f} (mod 1).")
L("   Test mit Realteil: τ(s) = √(1+χ) + i s/2π gegen −1/(N τ(1/s)), N = 4π² (Fricke-Form der Leiterspiegelung):")
for s in [0.3, 0.5, 1.0, 2.0]:
    t1 = complex(re_tau(s), s/(2*math.pi)); t2 = complex(re_tau(1/s), 1/(2*math.pi*s))
    f = -1/(4*math.pi**2*t2)
    L(f"   s = {s}: τ(s) = {t1.real:.4f}+{t1.imag:.4f}i,  −1/(4π² τ(1/s)) = {f.real:+.4f}+{f.imag:.4f}i  (mod 1: {f.real%1:.4f})")
L("   In Schreibweise Punkt 3 (|q| = Δs/r) ist die Spiegelung |q|·|q′| = (Δs/r*)², Im τ + Im τ′ = 2 Im τ(r*):")
L("   eine Spiegelung an einer waagrechten Geraden, keine Isometrie der oberen Halbebene, also kein Element von SL(2,ℝ).")

# ------------------------------------------------------------------------------------------------
L("\n1e Gemessene Relation (SPARC, Q ≤ 2, i ≥ 30°, Υ_d = 0,5, Υ_b = 0,7): Dualitätsdefekt")
L("    D(y) = lg x(y) − lg x(1/y) − 1,5·lg y  (y = g_bar/a₀, x = g_obs/a₀; Dualität ⇔ D = 0)")
HERE = os.path.dirname(os.path.abspath(__file__))
D0 = os.path.join(HERE, "..", "daten", "sparc")
tab = {}
for line in open(os.path.join(D0, "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
    p = line.split()
    if len(p) > 17 and p[1].isdigit() and len(p[0]) > 1:
        try: tab[p[0]] = dict(inc=float(p[5]), Q=int(p[17]))
        except ValueError: pass
GK = 1e6 / kpc
def lade(Yd, Yb=0.7):
    gb, go = [], []
    for fn in sorted(glob.glob(os.path.join(D0, "Rotmod_LTG", "*_rotmod.dat"))):
        t = tab.get(os.path.basename(fn).replace("_rotmod.dat", ""))
        if t is None or t["Q"] > 2 or t["inc"] < 30: continue
        d = np.loadtxt(fn); R, Vo, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
        g_b = (Vg*np.abs(Vg) + Yd*Vd*np.abs(Vd) + Yb*Vb**2) / R * GK
        m = (R > 0) & (Vo > 0) & (g_b > 0) & (eV < 0.1*Vo)
        gb.append(g_b[m]); go.append((Vo**2/R*GK)[m])
    return np.concatenate(gb), np.concatenate(go)

def bins(gb, go, a0, w=0.1):
    ly, lx = np.log10(gb/a0), np.log10(go/a0)
    cen = np.round(np.arange(-2.0, 2.0 + 1e-9, w), 3); out = {}
    for cc in cen:
        m = np.abs(ly - cc) < w/2
        if m.sum() >= 10: out[cc] = (lx[m].mean(), lx[m].std(ddof=1)/math.sqrt(m.sum()), m.sum(), ly[m].mean())
    return out

for Yd in (0.4, 0.5, 0.6):
    gb, go = lade(Yd); B = bins(gb, go, A0)
    if Yd == 0.5:
        L(f"   N = {len(gb)} Punkte; lg y von {np.log10(gb/A0).min():.2f} bis {np.log10(gb/A0).max():.2f}")
        L(f"   {'lg y':>6s} {'n(y)':>5s} {'n(1/y)':>6s} {'D Messung':>10s} {'±':>6s} {'D Regel':>8s}")
    rows = []
    for cc in sorted(B):
        if cc <= 0 or -cc not in B: continue
        (mx, ex, n1, my), (mx2, ex2, n2, my2) = B[cc], B[-cc]
        Dm = mx - mx2 - 0.75*(my - my2)
        yA, yB = 10**my, 10**my2
        Dr = math.log10(x_of_y(yA)) - math.log10(x_of_y(yB)) - 0.75*(my - my2)
        rows.append((cc, n1, n2, Dm, math.hypot(ex, ex2), Dr))
    if Yd == 0.5:
        for r_ in rows: L(f"   {r_[0]:6.2f} {r_[1]:5d} {r_[2]:6d} {r_[3]:+10.4f} {r_[4]:6.4f} {r_[5]:+8.4f}")
    Dm = np.array([r_[3] for r_ in rows]); eD = np.array([r_[4] for r_ in rows]); Dr = np.array([r_[5] for r_ in rows])
    L(f"   Υ_d = {Yd}: Paare {len(rows)}, χ²(D = 0) = {np.sum((Dm/eD)**2):8.1f},  χ²(D = Regel) = {np.sum(((Dm-Dr)/eD)**2):8.1f},"
      f"  mittleres D = {Dm.mean():+.3f} dex (Regel {Dr.mean():+.3f})")
L("   (Die Gewichtung 1,5·lg y ist mit den gemessenen Bin-Mitten gebildet: lg y − lg(1/y) = my − my2, also 0,75·(my − my2).)")

L("\n   Freier Spiegelpunkt a_m und freies Gewicht k: lg x(y) − lg x(1/y) = k′·(lg y − lg(1/y)), y = g_bar/a_m; Dualität ⇔ Gerade durch 0 (Υ_d = 0,5)")
gb, go = lade(0.5); res = []
for lam in np.linspace(-10.4, -9.4, 51):
    am = 10**lam; B = bins(gb, go, am); P = [(B[cc], B[-cc]) for cc in sorted(B) if cc > 0 and -cc in B]
    if len(P) < 6: continue
    X = np.array([(p[0][3] - p[1][3]) for p in P]); Y = np.array([p[0][0] - p[1][0] for p in P])
    E = np.array([math.hypot(p[0][1], p[1][1]) for p in P])
    k = np.sum(X*Y/E**2)/np.sum(X*X/E**2); c2 = np.sum(((Y - k*X)/E)**2)
    Xr = X; Yr = np.array([math.log10(x_of_y(10**p[0][3]*am/A0)) - math.log10(x_of_y(10**p[1][3]*am/A0)) for p in P])
    kr = np.sum(Xr*Yr/E**2)/np.sum(Xr*Xr/E**2); c2r = np.sum(((Yr - kr*Xr)/E)**2)
    res.append((lam, len(P), k, c2, kr, c2r))
b = min(res, key=lambda t: t[3] / max(t[1] - 2, 1))
L(f"   bester Spiegelpunkt Messung: a_m = {10**b[0]:.2e} m/s², Paare {b[1]}, k′ = {b[2]:.3f} (Gewicht 3/2 ⇔ 0,75), χ² = {b[3]:.1f} (Freiheitsgrade {b[1]-1})")
b2 = min(res, key=lambda t: t[5] / max(t[1] - 2, 1))
L("   Probe der Trennschärfe: dieselbe freie Anpassung an der Regel, die beweisbar nicht dual ist:")
L(f"   dieselbe Prüfung an der Regel (a₀ = 1,1e-10, fehlerfrei, gleiche Fehlerbalken): a_m = {10**b2[0]:.2e}, k′ = {b2[4]:.3f}, χ² = {b2[5]:.1f}")

L("   Freies k′ ≠ 0,75 widerspricht den beiden Ästen (Newton x = y, tief x = √y verlangen Gewicht 3/2); im Messbereich")
L("   (|lg y| ≤ 1,3) ist die freie Form zu schwach, um Dualität von Nicht-Dualität zu trennen.")

# ================================================================================================
kopf("2  τ = Re + i·Im je Ebene  (Re = Umlaufphase ω_rad/ω_azim bzw. Windungszahl; Im: P1 = s/2π, P3 = ln(r/Δs)/2π)")
# ------------------------------------------------------------------------------------------------
def reduziere(t):
    """SL(2,ℤ)-Reduktion auf den Fundamentalbereich |Re| ≤ ½, |τ| ≥ 1."""
    for _ in range(10000):
        t = complex(t.real - math.floor(t.real + 0.5), t.imag)
        if abs(t) < 1 - 1e-15: t = -1/t
        else: return t
    return t

MW = 6e10 * M_sun
ebenen = [  # Name, Masse, Radius (Hülle bzw. Bahn), Art
    ("Proton, Hülle R_p", m_p, R_p),
    ("Proton, Bohr-Radius (nur Gravitation)", m_p, a_B),
    ("Erde, Hülle", M_earth, R_earth),
    ("Sonne, Hülle", M_sun, R_sun),
    ("Planetenbahn 1 AE", M_sun, AU),
    ("Milchstraße 8 kpc", MW, 8*kpc),
    ("Haufen 1e14 M☉, R = r*/0,55", 1e14*M_sun, r_stern(1e14*M_sun)/0.55),
]
L(f"   {'Ebene':30s} {'s = r*/r':>10s} {'Re τ':>7s} {'Im τ P1':>10s} {'Im τ P3':>8s} {'δ* = (v_f/c)²':>14s} {'Im τ P3(r*)':>11s}  {'P1 reduziert':>18s}")
for name, M, r in ebenen:
    s = r_stern(M)/r; re = re_tau(min(s, 700)); i1 = s/(2*math.pi); i3 = math.log(r/ds(M))/(2*math.pi)
    dstar = math.sqrt(G*M*A0)/c**2; i3s = math.log(3/dstar)/(2*math.pi)
    tr = reduziere(complex(re, i1))
    L(f"   {name:30s} {s:10.3e} {re:7.4f} {i1:10.3e} {i3:8.3f} {dstar:14.3e} {i3s:11.3f}  {tr.real:+.4f}{tr.imag:+.4f}i")
t_star = complex(re_tau(1.0), 1/(2*math.pi)); tr = reduziere(t_star)
L(f"\n   An r* jeder Punktmasse (s = 1): τ* = {t_star.real:.5f} + {t_star.imag:.5f}i (P1), reduziert {tr.real:+.4f}{tr.imag:+.4f}i;"
  f"  j(τ*) nicht ausgezeichnet (weder i noch ρ).")
L("   P3 an r*: |q| = Δs/r* = δ*/3, Im τ_P3(r*) = ln(3/δ*)/2π; δ* ∝ M^{1/2}: Im τ_P3(r*) = const − ln M/(4π), stetig in M.")
L("   Leiter 3Δs : r* : ℓ₀ in P3: |q| = 1/3, δ*/3, δ*²/3 — jede Stufe q → δ*·q (Im τ + ln(1/δ*)/2π), eine Skalierung, kein τ → τ + 1.")
L(f"   Zelle (Obergrenze M_U, r* = 3Δs = ℓ₀, δ* = 1): Im τ_P3 = ln 3/2π = {math.log(3)/(2*math.pi):.4f}, Im τ_P1(r*) = 1/2π = {1/(2*math.pi):.4f};"
  "  Re τ = Wicklung W der Zelle, vom Modell nicht festgelegt (nur irrational).")
L("   Realteil: Kepler (Planeten, Hülle kompakter Körper) Re τ = 1, Tiefbereich √2, homogener Kern 2 (Newton) bzw. √3,5 (mit Regel,")
L("   a ∝ r^{1/2}); geschlossene Bahnen genau bei Re τ ∈ ℤ: 1 (Kepler) und 2 (harmonisch) — q gleich, τ um 1 verschoben.")
L("   Teilchen: Windungszahl n ∈ ℤ (Ladung) als Realteil, e^{2πin} = 1: die Ladung ändert q nicht, nur τ.")
L("\n   Paare von Ebenen unter SL(2,ℤ) gleich?  (gleicher reduzierter Punkt)")
pts = []
for name, M, r in ebenen:
    s = r_stern(M)/r; pts.append((name, reduziere(complex(re_tau(min(s, 700)), s/(2*math.pi)))))
dmin = min(abs(pts[i][1]-pts[j][1]) for i in range(len(pts)) for j in range(i+1, len(pts)))
L(f"   kleinster Abstand zweier reduzierter Punkte (P1): {dmin:.3e}  → keine zwei Ebenen SL(2,ℤ)-äquivalent außer über gleiches s.")

# ================================================================================================
kopf("3  α und α_G: eine Zahl Im τ_X?")
# ------------------------------------------------------------------------------------------------
aG = alpha_G(); yG = -math.log(aG)/(2*math.pi); ya = -math.log(alpha)/(2*math.pi)
L(f"   |q_G| = 3Δs_p/λ̄_p = α_G = {aG:.4e}, Im τ_G = {yG:.4f};  |q_α| = r_e/λ̄_e = α = {alpha:.6e}, Im τ_α = {ya:.4f}")
L(f"   Probe: 3Δs_p/λ̄_p = {3*ds(m_p)/compton_red(m_p):.4e};  r_e/λ̄_e = {r_e/compton_red(m_e):.6e}")
L("   Beide stehen im Modell als Beträge (Längenverhältnis Wirkungslänge/Quantenlänge), keiner als Phase.")
L("   Die Phase der Ladung ist die Windungszahl n ∈ ℤ: Re τ ganzzahlig, q unverändert — sie trägt keine reelle Zahl wie α.")
L(f"   Verhältnis Im τ_G/Im τ_α = {yG/ya:.4f};  Produkt = {yG*ya:.4f};  S-Bild von i·Im τ_α: i·{1/ya:.4f};  Fricke N = 1/(y_G y_α) = {1/(yG*ya):.4f}")
L(f"   α_G/α = {aG/alpha:.4e} = (3Δs_p/r_e)·(λ̄_e/λ̄_p) = (3Δs_p/r_e)·(m_p/m_e) = {3*ds(m_p)/r_e*(m_p/m_e):.4e}: enthält m_p/m_e und G.")
L("   Eine gemeinsame Zahl verlangte τ_G und τ_α als Bilder unter einer Abbildung, die die Struktur festlegt; keine gefunden.")

# ================================================================================================
kopf("5  Von Hand nachgerechnet (Werte zum Vergleich)")
L(f"   ν(0,5) = 1/(1 − e^{{−0,5}}) = {float(nu_s(0.5)):.5f};  ν(2) = {float(nu_s(2)):.5f};  ν(2)/(0,5·ν(0,5)) = {float(nu_s(2)/(0.5*nu_s(0.5))):.5f}")
L(f"   Re τ(r*) = √(1 + 1/(e − 1)) = {re_tau(1.0):.5f};  Milchstraße 6e10 M☉: v_f = {(G*MW*A0)**0.25/1e3:.2f} km/s, "
  f"Im τ_P3(r*) = {math.log(3/(math.sqrt(G*MW*A0)/c**2))/(2*math.pi):.4f}")
