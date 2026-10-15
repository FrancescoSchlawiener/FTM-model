# -*- coding: utf-8 -*-
"""Form der Zelle: welche Gitterformen tragen (a₀, Q) mit dem Abstand zur letzten Streufläche aus dem eigenen Atmen,
und woran trennen sich die Formen am Himmel.

Zelle (03): (1 + Q)·a₀ = 3c⁴/(2Gρ_b·V), V = L_min³/f³, ohne Kreispaare L_min ≥ 0,985·d_LSS ⇒ (1 + Q)·a₀ ≤ K·f³.
d_LSS aus dem eigenen Atmen: E_eff der Variante S2 (de_rest_anstieg, selbstkonsistent), D_M = c/H₀·∫dz/E_eff, ρ_b fest.
Formen (reine Verschiebungsgitter, bestes Glied je Familie): Würfel, sechseckig (Höhe = Seite), raumzentriert, flächenzentriert.
Eigenmoden −∇²ψ = k²ψ auf T³: k ∈ 2π·Γ* (duales Gitter). Winkelmoden der Hintergrundstrahlung nur im Sachs-Wolfe-Grenzfall
mit skaleninvariantem Spektrum (nur Geometrie): C_ℓ ∝ (1/V)Σ_k 2π² j_ℓ²(kχ)/k³, unendlicher Raum ∫dk/k j_ℓ² = 1/(2ℓ(ℓ+1)).
Aus ftm_kern: c, G, Gpc, OMEGA_B, H0_si, d_LSS, KREISPAAR_GRENZE, Q_KOSM. Aus de_rest_anstieg: E_eff (S2), D_M.
"""
import io, os, sys, math, contextlib, warnings
from itertools import product
import numpy as np
from scipy.special import spherical_jn, sph_harm_y
from scipy.stats import norm
from scipy.integrate import quad
warnings.filterwarnings("ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
import ftm_kern as K

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

c, G, Gpc = K.c, K.G, K.Mpc * 1e3
rho_b = K.OMEGA_B * 3 * K.H0_si()**2 / (8 * math.pi * G)
A_SP, sA_stat, sA_sys = 1.06e-10, 0.09e-10, 0.04e-10          # SPARC, eigene Daten (03)
sA = math.hypot(sA_stat, sA_sys)
Q_M, sQ = 5.36, 0.07                                          # [G]
s3 = math.sqrt(3)
GITTER = {"Würfel": np.eye(3),
          "sechseckig": np.array([[1, 0, 0], [0.5, s3 / 2, 0], [0, 0, 1.0]]),
          "raumzentriert": 0.5 * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1.0]]),
          "flächenzentriert": 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0.0]])}
GITTER = {n: A / abs(np.linalg.det(A))**(1 / 3) for n, A in GITTER.items()}          # Zeilen = Gittervektoren, Volumen 1

def vektoren(A, N):
    n = np.array([x for x in product(range(-N, N + 1), repeat=3) if any(x)])
    return n @ A
def f_form(A):
    return np.min(np.linalg.norm(vektoren(A, 2), axis=1))
F3 = {n: f_form(A)**3 for n, A in GITTER.items()}

# ------------------------------------------------------------------ 0 d_LSS aus dem eigenen Atmen
_src = open(os.path.join(HIER, "de_rest_anstieg.py"), encoding="utf-8").read().split("\n")
_stop = next(i for i, l in enumerate(_src) if l.startswith('kopf("2  Rückwirkung'))
DR = {"__file__": os.path.join(HIER, "de_rest_anstieg.py"), "__name__": "dra"}
with contextlib.redirect_stdout(io.StringIO()):
    exec("\n".join(_src[:_stop]), DR)
dmL = DR['dm_L']; E_S2 = DR['ERG']['S2']['E']; E_S1 = DR['ERG']['S1']['E']; E_V0 = DR['ERG']['V0']['E']
r_S2 = DR['D_M'](E_S2) / dmL
OMS = K.omega_vakuum(); OX0 = DR['ERG']['S2']['OX0']; luecke = K.OMEGA_DE - OMS - OX0
E_voll = lambda a: math.sqrt(E_S2(a)**2 + luecke)             # Rest 18 % als feste Rückstellung (nur Spannweite)
r_voll = DR['D_M'](E_voll) / dmL
dL0 = K.d_LSS(); dS2 = dL0 * r_S2

kopf("0  Abstand zur letzten Streufläche (mitbewegt, heutiger Beobachter, ρ_b fest)")
P("   Standardhintergrund: d_LSS = %.3f Gpc" % (dL0 / Gpc))
for nm, E in (("V0", E_V0), ("S1", E_S1), ("S2", E_S2)):
    P("   %s: D_M-Faktor %.4f  ⇒ d_LSS = %.3f Gpc" % (nm, DR['D_M'](E) / dmL, dL0 * DR['D_M'](E) / dmL / Gpc))
P("   S2 + Rest %.3f als feste Rückstellung (Spannweite, nicht Modell): Faktor %.4f" % (luecke, r_voll))
P("   weiter mit S2 (selbstkonsistent): d_LSS = %.3f Gpc, L_min = 0,985·d_LSS = %.3f Gpc" % (dS2 / Gpc, K.KREISPAAR_GRENZE * dS2 / Gpc))

def KZ(d):
    return 3 * c**4 / (2 * G * rho_b * (K.KREISPAAR_GRENZE * d)**3)
K0, K2, Kv = KZ(dL0), KZ(dS2), KZ(dL0 * r_voll)

# ------------------------------------------------------------------ 1 (a₀, Q)-Gebiet
kopf("1  (a₀, Q)-Gebiet je Form:  5 ≤ Q ≤ K·f³/a₀ − 1;  Prüfung mit a₀ = 1,06 ± 0,09 ± 0,04 (SPARC), Q = 5,36 ± 0,07 [G]")
P("   K = %.4e (Standard) → %.4e m/s² (S2, Faktor %.4f = 1/%.4f³)" % (K0, K2, K2 / K0, r_S2))
P("   Form               f³     Q_max(a₀=1,00) (1,06) (1,10) (1,13)  a₀,max(Q=5,36)  Fenster zu (Q=5)  Abstand/σ  P(erfüllt)  [Rest als Rückstellung]")
PV = {}
for n, f3 in F3.items():
    qm = [K2 * f3 / (x * 1e-10) - 1 for x in (1.00, 1.06, 1.10, 1.13)]
    am = K2 * f3 / (1 + Q_M)
    sig = math.hypot((1 + Q_M) * sA, A_SP * sQ)
    t = (K2 * f3 - (1 + Q_M) * A_SP) / sig
    tv = (Kv * f3 - (1 + Q_M) * A_SP) / sig
    PV[n] = norm.cdf(t)
    P("   %-17s %5.3f   %6.2f        %6.2f %6.2f %6.2f     %6.3f           %6.3f         %+5.2f      %5.3f       [%+5.2f]"
      % (n, f3, *qm, am * 1e10, K2 * f3 / 6 * 1e10, t, norm.cdf(t), tv))
P("   Abstand/σ = (K·f³ − (1 + Q)·a₀)/σ, σ² = (1 + Q)²σ_a² + a₀²σ_Q², σ_a = √(0,09² + 0,04²) = %.4f·10⁻¹⁰;" % (sA * 1e10))
P("   P(erfüllt) = Anteil der gemessenen (a₀, Q)-Verteilung innerhalb der Grenze.")
rng = np.random.default_rng(3)
aa = rng.normal(A_SP, sA, 400000); qq = rng.normal(Q_M, sQ, 400000)
P("   Probe Monte Carlo: " + ", ".join("%s %.3f" % (n, np.mean((1 + qq) * aa <= K2 * f3)) for n, f3 in F3.items()))
P("   SPARC-Spanne 1,00–1,13 (Υ_d, Gewichtung): Würfel Q_max %.2f–%.2f, sechseckig %.2f–%.2f."
  % (K2 / 1.13e-10 - 1, K2 / 1.00e-10 - 1, K2 * F3['sechseckig'] / 1.13e-10 - 1, K2 * F3['sechseckig'] / 1.00e-10 - 1))
P("   Modellwert a₀ = 1,1: Würfel Q ≤ %.2f (Q = 5,36 um %.1fσ_Q darüber), sechseckig Q ≤ %.2f."
  % (K2 / 1.1e-10 - 1, (Q_M - (K2 / 1.1e-10 - 1)) / sQ, K2 * F3['sechseckig'] / 1.1e-10 - 1))

# ------------------------------------------------------------------ 2 Kreispaare, Eigenmoden, Winkelmoden
kopf("2  Messbare Folgen der Form bei gleichem Volumen V = 3c⁴/(2Gρ_b(1 + Q)a₀), Beobachter heute, χ = d_LSS/2 (S2)")
chi = dS2 / 2
def volumen(a0, Q=Q_M):
    return 3 * c**4 / (2 * G * rho_b * (1 + Q) * a0)

def Yreal(l, th, ph):
    cols = []
    for m in range(-l, l + 1):
        y = sph_harm_y(l, abs(m), th, ph)
        cols.append(math.sqrt(2) * (-1)**m * (y.imag if m < 0 else y.real) if m != 0 else y.real)
    return np.array(cols)

def moden(A, V, Xmax=150.0):
    B = 2 * math.pi * np.linalg.inv(A).T / V**(1 / 3)
    kmax = Xmax / chi
    N = int(kmax / np.min(np.linalg.norm(B, axis=1)) * 1.8) + 2
    r = np.arange(-N, N + 1)
    n = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    kv = n @ B; k = np.linalg.norm(kv, axis=1)
    m = (k > 0) & (k < kmax)
    return kv[m], k[m], kmax

def schalen(k, tol=1e-6):
    ks = np.sort(k); out = []
    for x in ks:
        if out and abs(x - out[-1][0]) < tol * x: out[-1][1] += 1
        else: out.append([x, 1])
    return out

def tail(l, kmax):
    return quad(lambda x: spherical_jn(l, x)**2 / x, kmax * chi, np.inf, limit=400)[0]

def winkel(A, V, lmax=20, mit_achse=False):
    kv, k, kmax = moden(A, V)
    x = k * chi; w = 2 * math.pi**2 / (V * k**3)
    Cr = {}
    for l in range(2, lmax + 1):
        Cr[l] = 2 * l * (l + 1) * (np.sum(w * spherical_jn(l, x)**2) + tail(l, kmax))
    ach = {}
    th = np.arccos(np.clip(kv[:, 2] / k, -1, 1)); ph = np.arctan2(kv[:, 1], kv[:, 0])
    for l in (2, 3):
        Y = Yreal(l, th, ph)
        g = w * spherical_jn(l, x)**2
        M = (Y * g) @ Y.T * 4 * math.pi / (2 * l + 1) + np.eye(2 * l + 1) * tail(l, kmax) / (2 * l + 1)
        ev = np.linalg.eigvalsh(M); ach[l] = (ev / np.mean(ev), M / np.trace(M) * (2 * l + 1))
        # Anteil der ersten Schale an C_ℓ
    s1 = schalen(k)[0][0]; m1 = abs(k - s1) < 1e-6 * s1
    ant = {l: np.sum((w * spherical_jn(l, x)**2)[m1]) / (Cr[l] / (2 * l * (l + 1))) for l in (2, 3)}
    return Cr, ach, ant, schalen(k)[:4]

def S_halb(Cr, lmax):
    """S_{1/2} = ∫_{−1}^{1/2} C(θ)² dcosθ, C(θ) = Σ(2ℓ+1)/(4π)·C_ℓ·P_ℓ, C_ℓ ∝ R_ℓ/(ℓ(ℓ+1)) (Sachs-Wolfe)."""
    from numpy.polynomial import legendre as Lg
    mu = np.linspace(-1, 0.5, 3001)
    def C(Rf):
        co = np.zeros(lmax + 1)
        for l in range(2, lmax + 1): co[l] = (2 * l + 1) / (4 * math.pi) * Rf(l) / (l * (l + 1))
        return Lg.legval(mu, co)
    return np.trapezoid(C(lambda l: Cr[l])**2, mu) / np.trapezoid(C(lambda l: 1.0)**2, mu)

# Probe gegen die frühere Modensumme: Würfel mit L = d_LSS, C₂-Rest ≈ 0,84 (Forschungsstand 08)
Cw, _, _, _ = winkel(GITTER["Würfel"], dS2**3, lmax=4)
P("   Probe: Würfel mit L = d_LSS: C₂/C₂(∞) = %.3f, C₃ %.3f, C₄ %.3f (frühere Modensumme: ≈ 0,84)" % (Cw[2], Cw[3], Cw[4]))

FAELLE = [("a₀ = 1,06, Q = 5,36 (SPARC)", volumen(A_SP)), ("a₀ = 1,10, Q = 5,36 (Modellwert)", volumen(1.1e-10))]
ERGW = {}
for titel, V in FAELLE:
    P(); P("   %s:  V^⅓ = %.2f Gpc, d_LSS = %.2f Gpc" % (titel, V**(1 / 3) / Gpc, dS2 / Gpc))
    P("   Form               L_min[Gpc] L_min/d  Kreispaare (Radius)          k₁χ    g₁  k₂/k₁ g₂  C₂/C₂∞  C₃/C₃∞  C₄/C₄∞  S½/S½∞  1. Schale an C₂, C₃")
    for n, A in GITTER.items():
        Lm = f_form(A) * V**(1 / 3)
        g = vektoren(A, 2) * V**(1 / 3); gl = np.linalg.norm(g, axis=1)
        inn = gl[gl < dS2]
        kr = "keine" if len(inn) == 0 else "%d (%.1f°)" % (len(inn) // 2, math.degrees(math.acos(min(inn) / dS2)))
        Cr, ach, ant, sch = winkel(A, V)
        ERGW[(titel, n)] = (Cr, ach)
        P("   %-17s  %6.2f    %5.3f   %-26s  %5.2f  %3d  %5.3f %3d   %6.3f  %6.3f  %6.3f  %6.3f   %4.2f, %4.2f"
          % (n, Lm / Gpc, Lm / dS2, kr, sch[0][0] * chi, sch[0][1], sch[1][0] / sch[0][0], sch[1][1],
             Cr[2], Cr[3], Cr[4], S_halb(Cr, 20), ant[2], ant[3]))

P()
P("   Kreispaare erscheinen (kürzeste Schleife = d_LSS) bei a₀ ≥ 3c⁴/(2Gρ_b(1 + Q)(d_LSS/f)³), Q = 5,36, und mit dem Muster")
P("   der kürzesten Gittervektoren (Zentren der Kreise, Winkel zwischen den Achsen):")
for n, A in GITTER.items():
    g = vektoren(A, 2); gl = np.linalg.norm(g, axis=1); m = abs(gl - gl.min()) < 1e-9; gv = g[m] / gl[m][:, None]
    wk = sorted(set(round(math.degrees(math.acos(np.clip(abs(u @ v), 0, 1))), 1) for i, u in enumerate(gv) for v in gv[i + 1:] if abs(abs(u @ v) - 1) > 1e-9))
    a0c = 3 * c**4 / (2 * G * rho_b * (1 + Q_M) * (dS2 / f_form(A))**3)
    P("   %-17s a₀ ≥ %.3f·10⁻¹⁰:  %d Paare zugleich, Achsenwinkel %s°" % (n, a0c * 1e10, m.sum() // 2, ", ".join("%g" % x for x in wk)))

P()
P("   Ausrichtung (a₀ = 1,06): Eigenwerte der Kovarianz ⟨a_ℓm a_ℓm'⟩ im Achsensystem der Zelle, Mittel = 1 (unendlicher Raum: alle 1)")
T = "a₀ = 1,06, Q = 5,36 (SPARC)"
for n in GITTER:
    Cr, ach = ERGW[(T, n)]
    P("   %-17s ℓ = 2: %s   ℓ = 3: %s" % (n, " ".join("%.2f" % x for x in ach[2][0]), " ".join("%.2f" % x for x in ach[3][0])))
Mw = ERGW[(T, "Würfel")][1][2][1]; Ms = ERGW[(T, "sechseckig")][1][2][1]; Ms3 = ERGW[(T, "sechseckig")][1][3][1]
Mw3 = ERGW[(T, "Würfel")][1][3][1]
P("   Würfel ℓ = 2: Diagonale (m = −2 … 2, reell; m = 0 ∝ 3z²−r², m = 2 ∝ x²−y²): %s  ⇒ E_g (m = 0, 2) %.2f gegen T₂g %.2f"
  % (" ".join("%.2f" % Mw[i, i] for i in range(5)), (Mw[2, 2] + Mw[4, 4]) / 2, (Mw[0, 0] + Mw[1, 1] + Mw[3, 3]) / 3))
P("   Würfel ℓ = 3: Diagonale m = −3 … 3: %s" % " ".join("%.2f" % Mw3[i, i] for i in range(7)))
P("   sechseckig ℓ = 2 (Achse = z): Diagonale m = −2 … 2: %s  ⇒ m = 0 %.2f, |m| = 1 %.2f, |m| = 2 %.2f"
  % (" ".join("%.2f" % Ms[i, i] for i in range(5)), Ms[2, 2], (Ms[1, 1] + Ms[3, 3]) / 2, (Ms[0, 0] + Ms[4, 4]) / 2))
P("   sechseckig ℓ = 3: Diagonale m = −3 … 3: %s  ⇒ m = 0 %.2f, |m| = 3 %.2f"
  % (" ".join("%.2f" % Ms3[i, i] for i in range(7)), Ms3[3, 3], (Ms3[0, 0] + Ms3[6, 6]) / 2))
P("   sechseckig ℓ = 3: sin 3φ (m = −3) %.2f gegen cos 3φ (m = +3) %.2f: die Ecken des Sechsecks legen die Phase fest" % (Ms3[0, 0], Ms3[6, 6]))
P("   Trennschärfe eines Himmels (Achsen bekannt): KL = ½Σ(λ − 1 − ln λ) gegen den unendlichen Raum, ℓ = 2 + 3 [nats]: "
  + ", ".join("%s %.2f" % (n, sum(0.5 * np.sum(ERGW[(T, n)][1][l][0] - 1 - np.log(ERGW[(T, n)][1][l][0])) for l in (2, 3))) for n in GITTER))
P("   Anteil des Quadrupols in E_g (Würfel) bzw. m = 0 um die Achse (sechseckig), unendlicher Raum 0,40 bzw. 0,20:")
P("   Würfel %.3f,  sechseckig %.3f;  Oktupol sechseckig m = 0: %.3f (0,143), |m| = 3: %.3f (0,286)"
  % ((Mw[2, 2] + Mw[4, 4]) / 5, Ms[2, 2] / 5, Ms3[3, 3] / 7, (Ms3[0, 0] + Ms3[6, 6]) / 7))

# ------------------------------------------------------------------ 3 Hängt die Form an etwas?
kopf("3  Hängt die Form im Modell an etwas?")
P("   Volumen: fest durch a₀ und Q (Λ = ν·ΣΔs_b sieht nur V). Form: 5 Scherungen des flachen T³, ohne Rückstellung (Nullmoden).")
P("   Drehsymmetrien der Zelle (Kaleidoskop, 03): Ordnung 1, 2, 3, 4, 6; eine irrationale (goldene) Wicklung ist keine Drehung der Zelle.")
P("   Vergleich der Extreme (Rechnung, keine Herleitung):")
phi = (1 + math.sqrt(5)) / 2
nq = np.arange(1, 200001); dist = np.abs(nq * phi - np.round(nq * phi)) * nq
P("   Wicklung: min n·‖nW‖ über n ≤ 2·10⁵, golden: %.5f → 1/√5 = %.5f (Hurwitz: größter Rückkehrabstand, indefinite Form x² + xy − y²)"
  % (dist[1000:].min(), 1 / math.sqrt(5)))
P("   Gitter: kürzeste Schleife bei festem Volumen, f² = γ (Hermite): 2D Quadrat 1,000, Sechseck %.4f = (4/3)^½ (x² + xy + y²);" % (2 / s3)**0.5)
P("   3D: " + ", ".join("%s %.4f" % (n, f_form(A)) for n, A in GITTER.items()) + "; Maximum flächenzentriert = 2^{1/6} = %.4f" % 2**(1 / 6))
P("   Gälte für die Raumzelle dasselbe Prinzip wie für die Wicklung (größte kürzeste Rückkehr), wäre sie flächenzentriert;")
P("   das Modell setzt dieses Prinzip für die Zelle nicht. Die Zellform folgt nicht.")

# ------------------------------------------------------------------ 4 Beobachter
kopf("4  Beobachterfehler")
P("   Zelle, Form, V, L_min: G (Ganzes), mitbewegt fest; L atmet ∝ a, Λ und a₀ fest (I), Q seit z ≈ 10 fest.")
P("   d_LSS: je Beobachter und Epoche (Radius seines Lichtkegels bis z = 1090); die Grenze ist die von heute an einer mitbewegt")
dm1 = lambda E: np.trapezoid(1 / np.array([E(1 / (1 + x)) for x in np.linspace(0, 1, 2001)]), np.linspace(0, 1, 2001))
d_z1 = (DR['D_M'](E_S2) - dm1(E_S2)) / dmL * dL0
P("   festen Zelle: Beobachter bei z = 1 sähe d_LSS = %.2f Gpc (mitbewegt), L_min/d_LSS (sechseckig, a₀ = 1,06) = %.2f."
  % (d_z1 / Gpc, f_form(GITTER['sechseckig']) * volumen(A_SP)**(1 / 3) / d_z1))
P("   Verschiebungsgitter sind homogen: Muster der Kreise und Achsen gleich für jeden Ort, nur die Phasen δ_k hängen am Ort.")
P("   Die Vorzugsachsen sind Achsen des Ganzen (G), keine Richtung am Ort; die Senkung der großen Winkel ist Folge der fehlenden")
P("   Moden k < k₁ (Schließung), nicht einer Dichte. Sachs-Wolfe allein: das späte Atmen (integrierter Anteil) ist nicht gerechnet.")

# ------------------------------------------------------------------ 5 Hand
kopf("5  Von Hand")
P("   K(S2) = 7,1729·10⁻¹⁰/%.4f³ = %.4f·10⁻¹⁰;  Würfel, a₀ = 1,06: Q_max = %.4f/1,06 − 1 = %.3f" % (r_S2, 7.1729 / r_S2**3, 7.1729 / r_S2**3, 7.1729 / r_S2**3 / 1.06 - 1))
Vh = volumen(A_SP); a_h = (Vh / (s3 / 2))**(1 / 3)
P("   sechseckig, h = a: V = (√3/2)a³ ⇒ a = %.3f Gpc; k₁ = 2π/h = %.4f/Gpc, k₁χ = %.3f (g₁ = 2, längs der Achse);"
  % (a_h / Gpc, 2 * math.pi / a_h * Gpc, 2 * math.pi / a_h * chi))
P("   Würfel: L = V^⅓ = %.3f Gpc, k₁χ = 2π·%.3f/%.3f = %.3f (g₁ = 6)" % (Vh**(1 / 3) / Gpc, chi / Gpc, Vh**(1 / 3) / Gpc, 2 * math.pi * chi / Vh**(1 / 3)))

with open(os.path.join(HIER, "zelle_form.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
