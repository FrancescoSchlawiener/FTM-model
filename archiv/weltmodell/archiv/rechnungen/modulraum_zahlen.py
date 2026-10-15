# -*- coding: utf-8 -*-
"""
Modulraum und Zahlen: Prüfung dreier Vorschläge gegen das Modell (streng gegen Zahlensuche).
  1 m_p/m_e an CM-Punkten: was das Modell erzwingt; Suchraum einfacher Ausdrücke aus j, Δ, E₄, E₆, E₂* (|D| ≤ 20)
  2 Q über Phasenraumvolumina (Punktmasse unter der Regel gegen Newton); Abhängigkeit von Wahlen
  3 Isotropie der Fixpunkte von PSL(2,ℤ)
  5 Von Hand nachgerechnet
Konstanten aus ftm_kern.py. Ausgabe: modulraum_zahlen.txt
"""
import math
import random
import itertools
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from ftm_kern import m_p, m_e, hbar, c, G, ds, Q_KOSM, q_KOSM, R_p

mp.mp.dps = 40
MU = m_p / m_e
out = []


def P(s=""):
    out.append(s)
    print(s)


def kopf(t):
    P("=" * 100)
    P(t)
    P("=" * 100)


# ------------------------------------------------------------------ Modulformen (q-Reihen)
def sigma(n, k):
    return sum(d**k for d in range(1, n + 1) if n % d == 0)


def formen(tau, N=80):
    q = mp.exp(2j * mp.pi * tau)
    E2 = 1 - 24 * mp.fsum(sigma(n, 1) * q**n for n in range(1, N))
    E4 = 1 + 240 * mp.fsum(sigma(n, 3) * q**n for n in range(1, N))
    E6 = 1 - 504 * mp.fsum(sigma(n, 5) * q**n for n in range(1, N))
    D = q * mp.fprod((1 - q**n)**24 for n in range(1, N))
    j = E4**3 / D
    E2s = E2 - 3 / (mp.pi * mp.im(tau))           # nichtholomorph vervollständigt
    return dict(j=j, j1728=j - 1728, Delta=D, E4=E4, E6=E6, E2s=E2s)


def reduzierte_formen(D):
    """primitive reduzierte Formen (a,b,c), b² − 4ac = D < 0."""
    res = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            cc = (b * b - D) // (4 * a)
            if cc < a:
                continue
            if b < 0 and (abs(b) == a or a == cc):
                continue
            if math.gcd(math.gcd(a, abs(b)), cc) != 1:
                continue
            res.append((a, b, cc))
        a += 1
    return res


# ================================================================== 1
kopf("1  m_p/m_e an CM-Punkten")
P("1a Was das Modell erzwingt (vor jeder Rechnung festgelegt)")
P("   Das Modell ordnet einem Körper genau einen Modulparameter zu: den des Bahntorus bzw. der Hülle,")
P("   τ = Re (Umlaufphase ω_rad/ω_azim oder Windungszahl n) + i·Im, Im τ = ln(r/Δs)/2π (P3) bzw. s/2π (P1).")
P("   Beide Teile sind Zustandsgrößen (Radius, Bahn, Feld); ein fester CM-Punkt (i, ρ, …) ist keinem Teilchen")
P("   zugeordnet. 02: Proton aus dem Nome der starken Kopplung, Elektron aus seiner Yukawa-Kopplung — zwei")
P("   Mechanismen, kein gemeinsamer Torus. Erzwungen sind daher nur die τ, die aus Δs und einer Länge folgen:")
Dsp, Dse = ds(m_p), ds(m_e)
lp, le = hbar / (m_p * c), hbar / (m_e * c)
tauP_R = 1 + 1j * math.log(R_p / Dsp) / (2 * math.pi)
tauP_l = 1j * math.log(lp / Dsp) / (2 * math.pi)
tauE_l = 1j * math.log(le / Dse) / (2 * math.pi)
P(f"   Proton, Hülle R_p (P3):       τ = {tauP_R.real:.0f} + {tauP_R.imag:.4f}i   (|q| = Δs_p/R_p = {Dsp/R_p:.3e})")
P(f"   Proton, Quantenlänge λ̄_p:     τ = {tauP_l.imag:.4f}i   (|q| = Δs_p/λ̄_p = α_G/3 = {Dsp/lp:.4e})")
P(f"   Elektron, Quantenlänge λ̄_e:   τ = {tauE_l.imag:.4f}i   (|q| = Δs_e/λ̄_e = {Dse/le:.4e})")
P("   Elektron, Hülle: im Modell keine (Punktteilchen; 2ℓ_P/√α ist Ladungsradius, trägt die Masse nicht).")
P("   Diese τ liegen tief in der Spitze (Im τ ≈ 14–17), weit von jedem CM-Punkt kleiner Diskriminante.")
qp, qe = Dsp / lp, Dse / le
P(f"   An der Spitze Δ ≈ q, j ≈ 1/q:  Δ(τ_p)/Δ(τ_e) = (Δs_p/λ̄_p)/(Δs_e/λ̄_e) = (m_p/m_e)² = {qp/qe:.6e}  ((m_p/m_e)² = {MU**2:.6e})")
P("   → das erzwungene Verhältnis enthält m_p/m_e als Eingabe (Δs·λ̄ = ℓ_P²/3, Δs/λ̄ ∝ m²): Kreis, keine Bestimmung.")
P("   Ein CM-Punkt für Proton oder Elektron ist nicht erzwungen; jeder Treffer dort wäre gewählt.")
P("")

P("1b Suchraum (look-elsewhere): alle primitiven CM-Punkte mit |D| ≤ 20")
punkte = []
for D in range(-3, -21, -1):
    if D % 4 not in (0, 1):
        continue
    for (a, b, cc) in reduzierte_formen(D):
        tau = (-b + mp.sqrt(D)) / (2 * a)
        punkte.append((D, (a, b, cc), tau))
basis = []   # (Name, Wert > 0)
for D, f, tau in punkte:
    F = formen(tau)
    P(f"   D = {D:4d}  Form {f}  τ = {complex(tau).real:+.4f}{complex(tau).imag:+.4f}i  "
      f"j = {complex(F['j']).real:.6g}{complex(F['j']).imag:+.1e}i   |Δ| = {float(abs(F['Delta'])):.6e}")
    for k in ("j", "j1728", "Delta", "E4", "E6", "E2s"):
        v = float(abs(F[k]))
        if v > 1e-20:
            basis.append((f"|{k}|(D{D},{f[0]})", v))
P(f"   Punkte: {len(punkte)}, von null verschiedene Grundwerte (|j|, |j−1728|, |Δ|, |E₄|, |E₆|, |E₂*|): {len(basis)}")
P("   Probe (h = 1, j muss ganz sein): " + ", ".join(
    f"D{D}: {float(mp.re(formen(tau)['j'])):.4f}" for D, f, tau in punkte if D in (-7, -8, -11, -12, -16, -19)))

EXP = sorted({sg * p / q for p in range(1, 5) for q in (1, 2, 3, 4, 6, 12) for sg in (1, -1)})
P(f"   Exponenten e = ±p/q, p ≤ 4, q ∈ {{1,2,3,4,6,12}}: {len(EXP)} Werte")
lnB = np.array([math.log(v) for _, v in basis])
namen = [n for n, _ in basis]
E = np.array(EXP)

# Klasse K1: x^e;  K2: x^e1·y^e2 (x ≠ y, auch über zwei Punkte = „Proton hier, Elektron dort“);
# K3: n·x^e mit n = 2…12 (kleiner ganzzahliger Vorfaktor)
K1 = (lnB[:, None] * E[None, :]).ravel()
K1n = [(namen[i], EXP[k]) for i in range(len(namen)) for k in range(len(EXP))]
iu = [(i, j) for i in range(len(lnB)) for j in range(i + 1, len(lnB))]
K2 = (np.array([lnB[i] for i, j in iu])[:, None, None] * E[None, :, None]
      + np.array([lnB[j] for i, j in iu])[:, None, None] * E[None, None, :]).ravel()
K3 = (K1[:, None] + np.log(np.arange(2, 13))[None, :]).ravel()
klassen = [("K1  x^e", K1), ("K2  x^e1·y^e2", K2), ("K3  n·x^e (n = 2…12)", K3)]


def treffer(arr, ziel, tol):
    return int(np.sum(np.abs(arr - math.log(ziel)) < math.log1p(tol)))


random.seed(1)
ziele = [MU * math.exp(random.uniform(math.log(0.67), math.log(1.5))) for _ in range(400)]
P(f"   Ziel m_p/m_e = {MU:.5f}, Toleranz 10⁻³ (relativ). Kontrolle: 400 Zufallsziele im Faktor 0,67–1,5 um das Ziel.")
P(f"   {'Klasse':24s} {'Anzahl':>12s} {'Treffer':>8s} {'erwartet (Zufall)':>18s} {'Anteil Zufallsziele ≥ 1 Treffer':>32s}")
for name, arr in klassen:
    arr_s = np.sort(arr)
    def n_hit(z):
        lo, hi = math.log(z) - math.log1p(1e-3), math.log(z) + math.log1p(1e-3)
        return int(np.searchsorted(arr_s, hi) - np.searchsorted(arr_s, lo))
    h = n_hit(MU)
    zr = [n_hit(z) for z in ziele]
    P(f"   {name:24s} {arr.size:12d} {h:8d} {np.mean(zr):18.2f} {np.mean([x > 0 for x in zr]):32.2f}")
# beste Treffer in K1/K2 benennen
bestK1 = int(np.argmin(np.abs(K1 - math.log(MU))))
P(f"   bester K1: {K1n[bestK1][0]}^{K1n[bestK1][1]:+.4f} = {math.exp(K1[bestK1]):.4f}  (δ = {abs(math.exp(K1[bestK1])/MU-1):.1e})")
iK2 = int(np.argmin(np.abs(K2 - math.log(MU))))
pi_, rest = divmod(iK2, len(EXP) * len(EXP))
e1, e2 = divmod(rest, len(EXP))
P(f"   bester K2: {namen[iu[pi_][0]]}^{EXP[e1]:+.4f} · {namen[iu[pi_][1]]}^{EXP[e2]:+.4f} = {math.exp(K2[iK2]):.4f}"
  f"  (δ = {abs(math.exp(K2[iK2])/MU-1):.1e})")
P("   Suchraum zusammen ≈ %.2e Ausdrücke; jeder Treffer auf 10⁻³ ist hier Zufall ohne erzwungene Wahl." %
  (K1.size + K2.size + K3.size))
P("")

# ================================================================== 2
kopf("2  Q als Verhältnis von Phasenraumvolumina (Punktmasse)")
P("   Einheiten: GM = 1, r* = 1 (a₀ = GM/r*² = 1). Regel: a(r) = (1/r²)/(1 − e^{−1/r}), Potential Φ mit Φ' = a.")
P("   Tiefer Ast Φ ≈ v_f² ln r → ∞: unter der Regel ist JEDE Bahn gebunden, einen offenen Kanal gibt es nicht.")
P("   Endlich wird ein Verhältnis nur mit einem Fenster r < R (Wahl):")
P("   V_R(R) = Vol{r < R, v²/2 + Φ(r) < Φ(R)} (Regel, bleibt in R),  V_N(R) = dasselbe mit Newton,")
P("   Q_ph(R) = V_R/V_N − 1 (zusätzlich geschlossen gegen Newton-geschlossen).")


def a_regel_d(r, GM=1.0, rs=1.0):
    return (GM / r**2) / (-math.expm1(-rs / r))


def Phi(r, GM=1.0, rs=1.0):
    # Φ(r) − Φ(1·rs), Integral der Beschleunigung
    return quad(lambda x: a_regel_d(x, GM, rs), rs, r, limit=200)[0]


def V_regel(R, GM=1.0, rs=1.0):
    PR = Phi(R, GM, rs)
    f = lambda r: 4 * math.pi * r**2 * (4 * math.pi / 3) * max(2 * (PR - Phi(r, GM, rs)), 0) ** 1.5
    return quad(f, 1e-9 * rs, R, limit=200, points=[rs] if rs < R else None)[0]


def V_newton(R, GM=1.0):
    f = lambda r: 4 * math.pi * r**2 * (4 * math.pi / 3) * (2 * GM * (1 / r - 1 / R)) ** 1.5
    return quad(f, 0, R, limit=200)[0]


P(f"   {'R/r*':>8s} {'V_R/V_N':>12s} {'Q_ph':>10s} {'Q lokal q/(1−q), q = e^(−r*/R)':>34s}")
tab = []
for R in (0.1, 0.3, 1, 2, 3, 5, 10, 30, 100):
    vr, vn = V_regel(R), V_newton(R)
    ql = math.exp(-1 / R)
    tab.append((R, vr / vn - 1))
    P(f"   {R:8.2f} {vr/vn:12.5f} {vr/vn-1:10.4f} {ql/(1-ql):34.4f}")
# Wo Q_ph = 5,364?
from scipy.optimize import brentq
Rq = brentq(lambda R: V_regel(R) / V_newton(R) - 1 - Q_KOSM, 1, 100, xtol=1e-6)
P(f"   Q_ph = {Q_KOSM:.3f} bei R = {Rq:.4f} r*;  lokales Q = {Q_KOSM:.3f} bei R = r*/(−ln q_kosm) = {1/(-math.log(q_KOSM)):.4f} r*")
P("   Skalenfreiheit (M und a₀ anders, R/r* fest): ", )
for GM, a0 in ((7.0, 3.0), (0.2, 11.0)):
    rs = math.sqrt(GM / a0)
    R = 3 * rs
    vr = V_regel(R, GM, rs)
    vn = V_newton(R, GM)
    P(f"   GM = {GM}, a₀ = {a0}: r* = {rs:.4f}, R = 3r*:  V_R/V_N − 1 = {vr/vn-1:.6f}")
P("   → Q_ph hängt nicht von M, a₀, nur von R/r*, einer Wahl; es wächst ohne Grenze (R → ∞: V_N ∝ R^{3/2},")
P("     V_R ∝ R³(ln R)^{3/2}). Keine ausgezeichnete Zahl. Ebenso andere Fenster (Energieschnitt, Drehimpuls).")
P("   Symplektische Kapazität: Dimension einer Wirkung (Fläche), nicht eines 6-dim. Volumens; Gromovs Satz ist eine")
P("   Ungleichung (Kugel nicht in schmaleren Zylinder), keine Verhältniszahl. Ein Kapazitätsverhältnis zweier")
P("   skalenfreier Gebiete ist wieder eine Funktion von R/r*.")
P("")

# ================================================================== 3
kopf("3  Isotropie der Fixpunkte von PSL(2,ℤ)")
S = np.array([[0, -1], [1, 0]])
T = np.array([[1, 1], [0, 1]])
ST = S @ T


def ordnung_psl(M):
    X = np.eye(2, dtype=int)
    for k in range(1, 13):
        X = X @ M
        if np.array_equal(X, np.eye(2, dtype=int)) or np.array_equal(X, -np.eye(2, dtype=int)):
            return k
    return None


def moeb(M, z):
    return (M[0, 0] * z + M[0, 1]) / (M[1, 0] * z + M[1, 1])


rho = complex(mp.exp(2j * mp.pi / 3))
P(f"   S·i = {moeb(S, 1j):.6f} (Fixpunkt), Ordnung von S in PSL = {ordnung_psl(S)}")
P(f"   ST·ρ = {moeb(ST, rho):.6f} = ρ, Ordnung von ST in PSL = {ordnung_psl(ST)}")
P("   Stab(i) = ℤ₂, Stab(ρ) = ℤ₃: endliche zyklische Gruppen, 0-dimensional. SU(2), SU(3) sind zusammenhängende")
P("   Lie-Gruppen der Dimension 3 und 8. Gemeinsam ist nur: ℤ₂ = Zentrum von SU(2), ℤ₃ = Zentrum von SU(3);")
P("   ℤ₂ und ℤ₃ sind aber auch Untergruppen von U(1), SO(3), jeder Gruppe mit Elementen dieser Ordnung.")
# Abelisierung von PSL(2,Z) = Z6: Probe über Charaktere
P("   PSL(2,ℤ) = ℤ₂ * ℤ₃ (freies Produkt), Abelisierung ℤ₆; das Standardmodell hat (SU(3)×SU(2)×U(1))/ℤ₆.")
P("   Diese ℤ₆ folgt dort aus den Ladungen der Fermionen, nicht aus einer Wirkung; weder Lie-Algebra, Rang,")
P("   Dimension noch Darstellungen (3, 2, Hyperladung) folgen aus ℤ₂, ℤ₃. Ergebnis: Namensgleichheit der Zentren.")
P("   Im Modell: −I ∈ SL(2,ℤ) wirkt als z → −z, kehrt beide Windungen um (Vorzeichen der Ladung, Chiralität);")
P("   Q_el + t/3 ∈ ℤ (02) ist das ℤ₃-Zentrum der Farbe, gegeben, nicht aus Stab(ρ) abgeleitet.")
P("")

# ================================================================== 5
kopf("5  Von Hand nachgerechnet")
P("   Δ(i) = Γ(1/4)^24/(2^24·π^18):  ln = 24·1,288023 − 16,63553 − 20,60514 = −6,32812 → 1,7853·10⁻³")
P("   Reihe: q = e^{−2π} = 0,0018674; Δ ≈ q(1 − 24q + 252q²) = 0,0018674·0,956061 = 1,7853·10⁻³")
Fi = formen(mp.mpc(0, 1))
P(f"   Code: Δ(i) = {float(mp.re(Fi['Delta'])):.6e},  Γ-Formel = {float(mp.gamma(0.25)**24/(2**24*mp.pi**18)):.6e},"
  f"  j(i) = {float(mp.re(Fi['j'])):.6f}, E₆(i) = {float(abs(Fi['E6'])):.1e}")
P("   V_N(R) = (16π²/3)(2GM)^{3/2}·R^{3/2}·B(3/2, 5/2), B = Γ(1,5)Γ(2,5)/Γ(4) = 0,886227·1,329340/6 = 0,196350 = π/16")
P(f"   → V_N = (π³/3)(2GMR)^{{3/2}}; R = 3: {math.pi**3/3*(6)**1.5:.6f}, Code (Quadratur) {V_newton(3.0):.6f}")

with open("modulraum_zahlen.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
