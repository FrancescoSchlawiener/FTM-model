# -*- coding: utf-8 -*-
"""
Prüfung an den Axiomen: Muss das in gebundenen Bereichen inaktive Fluchtmaß in die parametrische Rolle (Atmen) übergehen?
Axiome als Aussagen über den Zustandsraum des Ganzen: (i) geschlossen, kein Außen; (ii) Maß des Ganzen erhalten
(Hamiltonscher Teil divergenzfrei; GENERIC-Teil 𝓜∇S erzeugt S im Ganzen); (iii) 𝒫 löscht keine Freiheitsgrade.
Minimales geschlossenes Modell: zwei Tori θ₁, θ₂ mit Kopplung −k·cos(θ₁ − θ₂) (Relativphase ψ, Impuls P; Libration = Verriegelung
T² → T¹ = Bindung) und eine parametrische Variable L (Schleifenlänge, Impuls π_L), N Bereiche teilen dasselbe L:
    H = Σ P_i²/(2L²) − k Σ cos ψ_i + π_L²/(2N·M₀) + ½N·K₀(L − 1)²   (L parametrisch: träge, M₀ = 400, K₀ = 0,02)
1 Liouville des Ganzen (Jacobi-Determinante), 2 gebundenes Maß bei festem und bei wanderndem Parameter (Separatrixfläche 16L√k),
3 Rückwirkung auf L: hängt sie am Maß oder an der Energie?  4 Verriegelung: κ/Ω der Kern-Regel, rationale Radien.
Aus ftm_kern: nu (Regel), A0. Sonst reine Modellrechnung ohne Messwerte.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import nu, A0

OUT = []
def P_(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P_(); P_("=" * 100); P_(t); P_("=" * 100)

M0, K0 = 400.0, 0.02            # je Bereich; parametrisch: träge (Eigenzeit 2π√(M/K) ≈ 888 ≫ Umlauf), M = N·M0, K = N·K0

def rhs(t, y, N, k, Lfix=None):
    psi, P, L, pL = y[:N], y[N:2 * N], y[2 * N], y[2 * N + 1]
    if Lfix is not None:
        L = Lfix(t)
    dpsi = P / L**2
    dP = -k(t) * np.sin(psi)
    dL = pL / (N * M0)
    dpL = np.sum(P**2) / L**3 - N * K0 * (L - 1.0)
    return np.concatenate([dpsi, dP, [dL, dpL]])

def lauf(psi0, P0, k, T, Lfix=None, L0=1.0, n_t=400):
    N = len(psi0)
    y0 = np.concatenate([psi0, P0, [L0, 0.0]])
    ts = np.linspace(0, T, n_t)
    s = solve_ivp(rhs, (0, T), y0, t_eval=ts, args=(N, k, Lfix), method="DOP853", rtol=1e-10, atol=1e-12)
    return s.t, s.y

def gebunden(psi, P, L, k):
    return (P**2 / (2 * L**2) - k * np.cos(psi)) < k

# ============================================================ 1 Liouville des Ganzen
kopf("1  Liouville des Ganzen: Jacobi-Determinante des Flusses (N = 2 Bereiche + L, 6 Dimensionen)")
kc = lambda t: 0.05
y0 = np.array([0.3, 2.0, 0.1, -0.05, 1.0, 0.0])
def fluss(y, T=200.0):
    s = solve_ivp(rhs, (0, T), y, args=(2, kc), method="DOP853", rtol=1e-12, atol=1e-13)
    return s.y[:, -1]
h = 1e-6
J = np.zeros((6, 6))
for j in range(6):
    e = np.zeros(6); e[j] = h
    J[:, j] = (fluss(y0 + e) - fluss(y0 - e)) / (2 * h)
P_("   det J(t = 200) = %.8f   (Ganzes Hamiltonsch: Maß erhalten)" % np.linalg.det(J))
P_("   Projektion ż = (I − Π⊥)F auf einem Teil hat div = −κ < 0;")
P_("   im Hamiltonschen Ganzen gilt punktweise div_Teil + div_Rest = 0: wer im Teil zieht, muss im Rest dehnen.")

# ============================================================ 2 gebundenes Maß
kopf("2  Gebundenes Maß (Libration) bei festem und bei wanderndem Parameter; Ensemble gleichverteilt in ψ, |P| < 1,2")
rng = np.random.default_rng(1)
N = 4000
psi0 = rng.uniform(-np.pi, np.pi, N); P0 = rng.uniform(-1.2, 1.2, N)
band = 2 * np.pi * 2.4
for name, k, Lfix in (
        ("k, L fest          ", lambda t: 0.05, lambda t: 1.0),
        ("k wächst, L fest   ", lambda t: 0.05 + 0.10 * t / 3000, lambda t: 1.0),
        ("k fest, L wächst   ", lambda t: 0.05, lambda t: 1.0 + 1.0 * t / 3000)):
    t, y = lauf(psi0, P0, k, 3000.0, Lfix=Lfix, n_t=7)
    s = "   %s" % name
    for i in (0, 3, 6):
        L = Lfix(t[i]); kk = k(t[i])
        f = gebunden(y[:N, i], y[N:2 * N, i], L, kk).mean()
        s += "  t=%4.0f: gebunden %.4f (16L√k/Band %.4f)" % (t[i], f, 16 * L * math.sqrt(kk) / band)
    P_(s)
P_("   Stichprobe ±0,007. Fester Parameter: gebundenes Maß zeitlich fest (Kac: was eintritt, tritt aus). Netto-Bindung nur, solange ein Parameter")
P_("   wandert, und genau um den Zuwachs der Separatrixfläche. Richtung: Parameter (auch L = Atmen) → Bindung, nicht umgekehrt.")

# ============================================================ 3 Rückwirkung auf L
kopf("3  Rückwirkung der Bindung auf L (L dynamisch, geschlossen): Maß oder Energie?")
P_("   Kraft auf L: Σ P²/L³ (nur der kinetische Teil hängt an L). Vergleich jeweils gegen denselben Start mit k = 0.")
N = 2000
fall = {
    "A Start oben  (ψ≈π, P klein)": (np.pi + rng.normal(0, 0.3, N), rng.normal(0, 0.05, N)),
    "B Start unten (ψ≈0, P mittel)": (rng.normal(0, 0.05, N), rng.normal(0, 0.18, N)),
    "C gleichverteilt, |P| < 0,1": (rng.uniform(-np.pi, np.pi, N), rng.uniform(-0.1, 0.1, N)),
}
Tl = 1800.0
for name, (ps, pp) in fall.items():
    L0 = 1.0
    res = {}
    for kk in (0.0, 0.02):
        t, y = lauf(ps, pp, lambda t, kk=kk: kk, Tl, L0=L0, n_t=900)
        L = y[2 * N]
        f = gebunden(y[:N, -1], y[N:2 * N, -1], L[-1], kk).mean() if kk > 0 else 0.0
        Ekin = np.mean(np.sum(y[N:2 * N]**2, axis=0) / (2 * L**2)) / N
        res[kk] = (L.mean(), f, Ekin)
    dL = res[0.02][0] - res[0.0][0]
    P_("   %s: gebunden %.3f  ⟨L⟩ k=0: %.5f  k=0,02: %.5f  Δ⟨L⟩ = %+.5f   ⟨E_kin⟩/Bereich %.5f → %.5f"
       % (name, res[0.02][1], res[0.0][0], res[0.02][0], dL, res[0.0][2], res[0.02][2]))
P_("   L wird verschoben, sobald die Bindung die Energie an L ändert; das Vorzeichen folgt ΔE_kin (A +, B −), nicht dem")
P_("   gebundenen Maß (A und B fast ganz gebunden, Verschiebung entgegengesetzt). Liouville legt die Verschiebung nicht fest.")

# ============================================================ 4 Verriegelung
kopf("4  Verriegelung: κ/Ω der Kern-Regel (Punktmasse, nur x = r/r* zählt) und Radien mit rationalem Verhältnis")
def k_om(x):
    Om2 = lambda u: float(nu(A0 / u**2) * A0 / u**2) / u     # Einheiten r* = 1, a₀ = 1 fällt heraus
    hh = 1e-5
    d = (np.log(Om2(x * (1 + hh))) - np.log(Om2(x * (1 - hh)))) / (2 * hh)
    return math.sqrt(d + 4)
P_("   " + "  ".join("x=%g: %.4f" % (x, k_om(x)) for x in (0.01, 0.1, 0.5, 1, 2, 5, 20, 100, 1000)))
for p, qq in ((9, 8), (6, 5), (5, 4), (4, 3), (7, 5)):
    r = p / qq
    try:
        x = brentq(lambda u: k_om(u) - r, 0.01, 1e5)
        P_("   κ/Ω = %d/%d = %.4f bei r/r* = %.3f (geschlossene Rosette: %d Radialschwingungen je %d Umläufe)" % (p, qq, r, x, p, qq))
    except ValueError:
        P_("   κ/Ω = %d/%d nicht erreicht" % (p, qq))
P_("   Von Hand: Kepler Ω² ∝ r⁻³ → κ² = −3Ω² + 4Ω² → κ/Ω = 1 (alle Bahnen geschlossen, Entartung = Verriegelung);")
P_("   flach Ω² ∝ r⁻² → κ² = 2Ω² → √2. Dazwischen durchläuft κ/Ω stetig jedes Verhältnis in (1, √2).")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
