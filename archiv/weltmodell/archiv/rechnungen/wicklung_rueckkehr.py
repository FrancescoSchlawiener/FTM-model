# -*- coding: utf-8 -*-
# Wicklung und Rückkehr: Kann die eindimensionale Rückkehr der Regel (q = e^{−√(Γτ_D)}) eine irrationale
# Wicklung des Torus sein, oder ist sie eine eigene innere Richtung?
# 1 Rückkehrzeiten einer Kreisdrehung (Kac, drei Werte)   2 G(Γ) = ⟨e^{−ΓT}⟩, Exponent von 1 − G, mit Rauschen
# 3 Zelle am Rand, von der Drehung getriebener Lauf, Irrweg auf der bei Λ geschlossenen Halbgeraden, Drift + Diffusion
# 4 Kac-Budget, Umlaufsinn, Wicklung   5 unabhängige Prüfungen (analytisch)
import sys, io, math
import numpy as np
from scipy.integrate import quad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ftm_kern import A0, c, q_, s_

out = []; P = lambda s="": out.append(s)
PHI = (math.sqrt(5) - 1) / 2
ALPHAS = {"golden (√5−1)/2": PHI, "√2 − 1": math.sqrt(2) - 1, "e − 2": math.e - 2}
Gam = np.logspace(-8, 2, 401)                      # Verlustrate in 1/Schritt


def steigung(G1):
    """lokaler Exponent d ln(1 − G)/d ln Γ auf dem Gitter Gam"""
    return np.gradient(np.log(G1), np.log(Gam))


def halb_fenster(sl, tol=0.05):
    """Dekaden, in denen der lokale Exponent in [½ − tol, ½ + tol] liegt"""
    return np.sum(np.abs(sl - 0.5) < tol) * (np.log10(Gam[1]) - np.log10(Gam[0]))


def rueckkehr_drehung(alpha, mu, n=200_000):
    """erste Rückkehrzeiten nach [0, μ) für x gleichverteilt in [0, μ)"""
    x0 = (np.arange(n) + 0.5) / n * mu
    T = np.zeros(n, int); x = x0.copy(); offen = np.ones(n, bool); k = 0
    while offen.any():
        k += 1
        x[offen] = (x[offen] + alpha) % 1.0
        treffer = offen & (x < mu)
        T[treffer] = k; offen &= ~treffer
    return T


# ---------------- 1 ----------------
P("1) Kreisdrehung x → x + α mod 1, Zelle [0, μ): erste Rückkehrzeiten (2·10⁵ Startpunkte gleichverteilt in der Zelle)")
rot = {}
for name, al in ALPHAS.items():
    for mu in (0.1, 0.01, 0.001):
        T = rueckkehr_drehung(al, mu)
        w, h = np.unique(T, return_counts=True); p = h / h.sum()
        rot[(name, mu)] = (w, p)
        drei = f"T₃ − T₁ − T₂ = {w[2]-w[0]-w[1]}" if len(w) == 3 else ""
        P(f"   α = {name:15s} μ = {mu:<6}: ⟨T⟩·μ = {np.sum(w*p)*mu:.6f}   Werte {[int(t) for t in w]} Anteile {np.round(p,4).tolist()} {drei}")
P("   Kac ⟨T⟩ = 1/μ erfüllt; höchstens drei Werte, der größte ist die Summe der beiden anderen (Drei-Abstände-Satz).")
P("   Golden: die Werte sind aufeinanderfolgende Fibonacci-Zahlen (μ = 0,01: 55, 89, 144; μ = 0,001: 610, 987, 1597).")
P("   Normierte längste Rückkehr μ·T_max über μ ∈ [10⁻³, 10⁻¹] (40 Werte, 2·10⁴ Startpunkte):")
mus = np.logspace(-3, -1, 40)
for name, al in ALPHAS.items():
    r = [m * rueckkehr_drehung(al, m, 20_000).max() for m in mus]
    P(f"     α = {name:15s}: max μ·T_max = {max(r):.3f}, Mittel {np.mean(r):.3f}")
P()

# ---------------- 2 ----------------
P("2) G(Γ) = ⟨e^{−ΓT}⟩ (T in Schritten, Schrittdauer 1). Lokaler Exponent d ln(1 − G)/d ln Γ.")
P("   Regel: 1 − e^{−√(Γτ_D)} hat für Γτ_D → 0 den Exponenten ½ über beliebig viele Dekaden.")
tauD = 100.0
Gr = -np.expm1(-np.sqrt(Gam * tauD)); slr = steigung(Gr)
P(f"   Regel (τ_D = 100): Exponent bei Γτ_D = 10⁻⁴: {np.interp(1e-6, Gam, slr):.4f}; bei Γτ_D = 1: {np.interp(1e-2, Gam, slr):.3f}; "
  f"Dekaden mit |Exp − ½| < 0,05 im Gitter (Γτ_D ≥ 10⁻⁶): {halb_fenster(slr):.1f}, ohne untere Grenze")
for name in ALPHAS:
    w, p = rot[(name, 0.01)]
    G1 = 1 - np.exp(-np.outer(Gam, w)) @ p
    sl = steigung(G1)
    P(f"   Drehung α = {name:15s} μ = 0,01: Exponent bei Γ = 10⁻⁶: {np.interp(1e-6, Gam, sl):.4f}, 10⁻²: {np.interp(1e-2, Gam, sl):.3f}, "
      f"10⁰: {np.interp(1.0, Gam, sl):.3f}; Dekaden mit |Exp − ½| < 0,05: {halb_fenster(sl):.2f}")
P("   Drehung ohne Rauschen: 1 − G ≈ Γ⟨T⟩ (Exponent 1), bei Γ ≳ 1/T₁ Übergang auf 0. Der Wert ½ wird nur im Durchgang berührt.")

rng = np.random.default_rng(1)


def rueckkehr_rauschen(alpha, sig, mu, n=200_000, tmax=10**6):
    x = rng.uniform(0, mu, n); T = np.zeros(n, np.int64); offen = np.arange(n); k = 0
    while offen.size and k < tmax:
        k += 1
        x[offen] = (x[offen] + alpha + sig * rng.standard_normal(offen.size)) % 1.0
        tr = x[offen] < mu
        T[offen[tr]] = k; offen = offen[~tr]
    return T


P("   Drehung mit Phasenrauschen x → x + α + σξ (ξ normal), μ = 0,01, 2·10⁵ Läufe:")
for al, an, sigs in ((PHI, "golden", (0.003, 0.01, 0.03, 0.1)), (PHI * 1e-2, "golden·10⁻²", (0.003, 0.01, 0.03, 0.1)),
                     (PHI * 1e-5, "golden·10⁻⁵", (0.001, 0.003))):
    for sig in sigs:
        T = rueckkehr_rauschen(al, sig, 0.01)
        w, h = np.unique(T, return_counts=True); p = h / h.sum()
        G1 = 1 - np.exp(-np.outer(Gam, w)) @ p
        sl = steigung(G1)
        Pe = al * 1.0 / sig**2
        P(f"     α = {an:12s} σ = {sig:<6} (Pe = αL/σ² = {Pe:9.3g}): ⟨T⟩μ = {T.mean()*0.01:.3f}, P(T=1) = {np.mean(T==1):.3f}, "
          f"Exp bei Γ=10⁻⁶: {np.interp(1e-6, Gam, sl):.3f}, 10⁻³: {np.interp(1e-3, Gam, sl):.3f}, 10⁻¹: {np.interp(1e-1, Gam, sl):.3f}; "
          f"Dekaden |Exp−½|<0,05: {halb_fenster(sl):.2f}")
P("   Rauschen verschmiert nur die drei Werte; solange die Drehung die Bewegung trägt (Pe ≫ 1), bleibt kein Bereich mit ½.")
P("   Ein ½-Bereich entsteht erst, wenn die Drift fast null ist (Pe ≪ 1) und die Schritte kleiner als die Zelle: dann ist es ein Irrweg.")
P()

# ---------------- 3 ----------------
P("3a) Zelle am Rand: gefaltete Drehung (Spiegelung an den Enden von [0, 1] = Drehung auf dem Kreis der Länge 2), Zelle [0, μ) am Rand.")
for name, al in ALPHAS.items():
    mu = 0.01; n = 400_000
    y = (np.arange(n) + 0.5) / n * mu          # Lift auf dem Kreis der Länge 2; Zelle = (−μ, μ) mit Maß μ je Seite
    y = np.concatenate([y, 2 - y]); T = np.zeros(2 * n, int); off = np.ones(2 * n, bool); k = 0
    while off.any():
        k += 1; y[off] = (y[off] + 2 * al) % 2.0
        tr = off & ((y < mu) | (y > 2 - mu)); T[tr] = k; off &= ~tr
    w, h = np.unique(T, return_counts=True); p = h / h.sum()
    G1 = 1 - np.exp(-np.outer(Gam, w)) @ p; sl = steigung(G1)
    P(f"   α = {name:15s}: ⟨T⟩μ = {np.sum(w*p)*mu:.4f}, Werte {[int(t) for t in w]}, Exp bei Γ = 10⁻⁶: {np.interp(1e-6, Gam, sl):.4f}, "
      f"Dekaden |Exp−½|<0,05: {halb_fenster(sl):.2f}")
P("   Der Rand ändert nichts: gefaltete Drehung ist wieder eine Drehung, drei Werte, Exponent 1 → 0.")
P()

P("3b) Von der Drehung getriebener Lauf (deterministisch): y_{n+1} = y_n + f(x_n), f = +1 auf [0, ½), −1 auf [½, 1), x → x + α;")
P("    Halbgerade |y|, Rückkehr nach y = 0. 2·10⁴ Startphasen x gleichverteilt, bis 2·10⁵ Schritte.")
NMAX = 200_000
for name, al in ALPHAS.items():
    n = 20_000; x = rng.uniform(0, 1, n); yv = np.zeros(n, np.int64); T = np.full(n, -1, np.int64)
    idx = np.arange(n); k = 0; ymax = 0
    while idx.size and k < NMAX:
        k += 1
        yv[idx] += np.where(x[idx] < 0.5, 1, -1)
        x[idx] = (x[idx] + al) % 1.0
        ymax = max(ymax, int(np.abs(yv[idx]).max()))
        tr = yv[idx] == 0
        T[idx[tr]] = k; idx = idx[~tr]
    ret = T[T > 0]
    w, h = np.unique(ret, return_counts=True)
    ts = np.array([10, 100, 1000, 10000])
    surv = [np.mean((T > t) | (T < 0)) for t in ts]
    G1 = 1 - np.exp(-np.outer(Gam, w)) @ (h / h.sum()); sl = steigung(G1)
    P(f"   α = {name:15s}: Exp von 1 − G bei Γ = 10⁻⁶, 10⁻⁴, 10⁻², 10⁰: {[round(float(np.interp(g, Gam, sl)),3) for g in (1e-6,1e-4,1e-2,1)]}, "
      f"Dekaden |Exp−½|<0,05: {halb_fenster(sl):.2f}; Verhältnis aufeinanderfolgender Werte {w[4]/w[3]:.3f}, {w[5]/w[4]:.3f}")
    P(f"   α = {name:15s}: Rückkehr bis 2·10⁵: {ret.size/n:.4f}, ⟨T⟩ = {ret.mean():.2f}, verschiedene Werte {w.size} ({[int(t) for t in w[:6]]}…), "
      f"max|y| = {ymax}, P(T > 10, 10², 10³, 10⁴) = {[round(float(s),4) for s in surv]}")
P("    Irrweg zum Vergleich: P(T > t) ≈ √(2/(πt)) = " + str([round(math.sqrt(2/(math.pi*t)), 4) for t in (10, 100, 1000, 10000)]))
P("    Der von der Wicklung getriebene Lauf streut nur logarithmisch (Ostrowski); die Rückkehrzeiten bilden eine Leiter")
P("    (golden: Verhältnis → φ³ = 4,236; √2 − 1: (1 + √2)² = 5,828), das Mittel ist endlich, kein t^{−3/2}-Schwanz, keine Wurzel.")
P()

P("3c) Irrweg auf der bei Λ geschlossenen Halbgeraden: Orte 0 … N (Schrittlänge ℓ, Schrittdauer ℓ/c = 1), Spiegel bei N = Λ/ℓ.")


def F_erste(N, gam, m, einsminus=False):
    """E_m[z^{τ₀}], z = e^{−Γ}, einfacher Irrweg auf {0 … N}, spiegelnd bei N; mit einsminus: 1 − E_m[…] (genau)."""
    from scipy.linalg import solve_banded
    z = math.exp(-gam); e = -math.expm1(-gam)
    # v_i = 1 − u_i:  v_i − z(½v_{i−1} + ½v_{i+1}) = 1 − z,  v_0 = 0,  v_N − z·v_{N−1} = 1 − z
    ab = np.zeros((3, N)); ab[1] = 1.0; ab[0, 1:] = -z / 2; ab[2, :-1] = -z / 2; ab[2, N - 2] = -z
    v = solve_banded((1, 1), ab, np.full(N, e))
    return v[m - 1] if einsminus else 1 - v[m - 1]


def kac_mittel(N):
    """mittlere Rückkehrzeit zum Ort, exakt: m_i = 1 + ½m_{i−1} + ½m_{i+1}, m_0 = 0, m_N = 1 + m_{N−1}; ⟨T⟩ = 1 + m_1"""
    from scipy.linalg import solve_banded
    ab = np.zeros((3, N)); ab[1] = 1.0; ab[0, 1:] = -0.5; ab[2, :-1] = -0.5; ab[2, N - 2] = -1.0
    return 1 + solve_banded((1, 1), ab, np.ones(N))[0]


GamL = np.logspace(-11, 1, 481)
for N in (100, 1000, 10000):
    E1 = np.array([-math.expm1(-g) + math.exp(-g) * F_erste(N, g, 1, True) for g in GamL])   # 1 − G, erst Schritt nach 1
    sl = np.gradient(np.log(E1), np.log(GamL))
    dek = np.sum(np.abs(sl - 0.5) < 0.05) * (np.log10(GamL[1]) - np.log10(GamL[0]))
    P(f"   N = {N:6d}: Kac ⟨T⟩ = {kac_mittel(N):.4f} (2N = {2*N}); Exp bei Γ = 0,1/N²: {np.interp(0.1/N**2, GamL, sl):.3f}, "
      f"Γ = 10/N²: {np.interp(10/N**2, GamL, sl):.3f}, Γ = 10⁻¹: {np.interp(0.1, GamL, sl):.3f}; Dekaden |Exp−½|<0,05: {dek:.2f}")
P("   Exponent 1 unterhalb Γ ≈ D/Λ² (Kac, endliches Mittel), ½ zwischen D/Λ² und c/ℓ (Pólya), das Fenster wächst wie 2·log₁₀N.")

P("   Rückkehr aus der Tiefe ϖ (Ort m) mit Spiegel bei Λ, kontinuierlich: G = cosh((Λ − ϖ)k)/cosh(Λk), k = √(Γ/D):")
for m, N in ((10, 1000), (10, 10000)):
    D = 0.5
    for gam in (1e-7, 1e-4, 1e-2):
        k = math.sqrt(gam / D); ex = F_erste(N, gam, m)
        ko = math.exp(-m * k) * (1 + math.exp(-2 * (N - m) * k)) / (1 + math.exp(-2 * N * k))
        P(f"     ϖ = {m}, Λ = {N}, Γ = {gam:g}: diskret {ex:.6f}, kontinuierlich {ko:.6f}, offen e^{{−ϖk}} = {math.exp(-m*k):.6f}")
P()

P("3d) Drift und Diffusion (½ liegt zwischen max(D/L², v²/4D) und D/ϖ²) auf dem Kreis (Umfang L, Zelle am Ort, Start in der Tiefe ϖ hinter der Zelle): D u'' + v u' = Γ u, u = 1 an beiden Zellrändern.")


def G_drift(gam, v, D, L, w):
    """exakt: Laplace-Transformierte der ersten Ankunft am Ort von x = ϖ auf (0, L), Drift v, Diffusion D."""
    disc = np.sqrt(v * v + 4 * D * gam)
    rp = (-v + disc) / (2 * D); rm = (-v - disc) / (2 * D)
    # u = A e^{rp x} + B e^{rm x}, u(0) = u(L) = 1; stabil geschrieben
    x = w
    E1 = np.exp(rp * (x - L)); E2 = np.exp(-rp * L); F1 = np.exp(rm * x); F2 = np.exp(rm * L)
    # A e^{rp L} = a', u(x) = a' e^{rp(x−L)} + B e^{rm x}; a' + B e^{rm L} = 1; a' e^{−rp L} + B = 1
    det = 1 - E2 * F2
    a_ = (1 - F2) / det; B = (1 - E2) / det
    return a_ * E1 + B * F1


L = 1.0; wdep = 1e-3
for Pe in (0, 1, 10, 100, 1000):
    D = 1.0; v = Pe * D / L
    G1 = 1 - np.array([G_drift(g, v, D, L, wdep) for g in Gam * 1e4]); sl = steigung(G1)
    P(f"   Pe = vL/D = {Pe:5d}, ϖ/L = 10⁻³: Exp bei ΓL²/D = 10⁻², 10⁰, 10², 10⁴: "
      f"{[round(float(np.interp(g, Gam*1e4, sl)),3) for g in (1e-2, 1, 1e2, 1e4)]}; Dekaden |Exp−½|<0,05: {halb_fenster(sl):.2f}")
P("   Ohne Drift ist der Kreis vom Umfang L = 2Λ genau die bei Λ geschlossene Halbgerade (Spiegelsymmetrie, vgl. 3c).")
P("   Eine Drift (die Wicklung als Umlauf) verkürzt das ½-Fenster; bei Pe ≳ 10³ bleibt keins.")
P()

# ---------------- 4 ----------------
P("4) Kac-Budget, Umlaufsinn, Wicklung")
cuts = np.array([0, 0.05, 0.12, 0.3, 0.41, 0.6, 0.77, 1.0])
for name, al in ALPHAS.items():
    xs = (np.arange(3_000_000) * al) % 1.0
    lab = np.searchsorted(cuts, xs, side="right") - 1
    summe = 0.0
    for i in range(len(cuts) - 1):
        pos = np.nonzero(lab == i)[0]; summe += 1 / np.mean(np.diff(pos))
    P(f"   Zerlegung des Kreises in 7 Zellen, α = {name:15s}: Σ 1/⟨T_i⟩ = {summe:.6f}")
for name, al in ALPHAS.items():
    w1, p1 = rot[(name, 0.01)]
    T2 = rueckkehr_drehung(1 - al, 0.01); w2, h2 = np.unique(T2, return_counts=True)
    P(f"   Umlaufsinn α ↔ −α ({name}, μ = 0,01): Werte {[int(t) for t in w1]} / {[int(t) for t in w2]}, Anteile gleich: {np.allclose(p1, h2/h2.sum(), atol=1e-5)}")
N = 1000
P(f"   Halbgerade N = {N}: Kac-Anteil des Ortes π₀ = 1/(2N) = {1/(2*N):.6f}; Σ π_i = (1 + 2(N−1) + 1)/(2N) = 1")
P()

# ---------------- 5 ----------------
P("5) Unabhängige Prüfungen")
for t0, gam in ((1.0, 0.3), (25.0, 1e-3), (0.25, 4.0)):
    f = lambda t: math.sqrt(t0 / math.pi) * t**-1.5 * math.exp(-t0 / t - gam * t)
    I = quad(f, 0, np.inf, limit=400)[0]
    P(f"   ∫ √(t₀/π) t^{{−3/2}} e^{{−t₀/t}} e^{{−Γt}} dt, t₀ = {t0}, Γ = {gam}: {I:.10f}; e^{{−√(4Γt₀)}} = {math.exp(-math.sqrt(4*gam*t0)):.10f}")
tD = c / A0
for aN in (A0, 0.01 * A0, 4 * A0):
    gam = aN / c; t0 = tD / 4          # in Einheiten von τ_D rechnen: u = t/τ_D
    f = lambda u: math.sqrt(0.25 / math.pi) * u**-1.5 * math.exp(-0.25 / u - gam * tD * u)
    I = quad(f, 0, 1, limit=400)[0] + quad(f, 1, np.inf, limit=400)[0]
    P(f"   Regel: a_N/a₀ = {aN/A0:g}, Γ = a_N/c, t₀ = τ_D/4, τ_D = c/a₀ = {tD:.3e} s: Integral {I:.8f}, ftm_kern q = {float(q_(aN)):.8f}")
P("   Kac ohne Rechner: Halbgerade 0 … N spiegelnd, Stationärmaß ∝ Grad (1, 2, …, 2, 1), Summe 2N ⇒ ⟨T₀⟩ = 2N (3c exakt gelöst).")
P("   Drei Werte von Hand, golden, μ = 0,01: 55α = 33,9919… → Rest 0,99187; 89α = 55,00502 → Rest 0,00502 < μ; 144α = 88,99689 → Rest 0,99689.")
P("   Von x ∈ [0, μ): Rückkehr nach 89 für x < μ − 0,00502, nach 55 für x ≥ 1 − 0,99187 = 0,00813, dazwischen 144 ⇒ Anteile")
g55 = 1 - (55 * PHI) % 1; g89 = (89 * PHI) % 1
p89 = (0.01 - g89) / 0.01; p55 = (0.01 - g55) / 0.01; p144 = 1 - p89 - p55
P(f"   p₅₅ = {p55:.4f}, p₈₉ = {p89:.4f}, p₁₄₄ = {p144:.4f}; Kac 55p₅₅ + 89p₈₉ + 144p₁₄₄ = {55*p55+89*p89+144*p144:.4f} (1/μ = 100)")

txt = "\n".join(out)
print(txt)
open(__file__.replace(".py", ".txt"), "w", encoding="utf-8").write(txt + "\n")
