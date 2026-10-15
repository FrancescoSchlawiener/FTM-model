# -*- coding: utf-8 -*-
"""
Torus-Dynamik-Schicht des FTM-Weltmodells - Rechenbeispiel
==========================================================

Grundlage (Autor, "Multiphysik Torus Kopplungen"):

    Zustandsraum   M = T^n x S
    Phasen         dθ_i/dt = ω_i(s) + K_i(θ, s)
    Struktur       ds/dt   = G(θ, s) - D(s) s

    Kanal A (Phase-Phase):  K_i = Σ_j k_ij(s) sin(θ_j - θ_i - φ_ij)
    Kanal B (Medium->Phase): ω_i = ω_i(s)
    Kanal C (Phase->Medium): G(θ, s) mit Produkten schneller Moden
    Kanal D (Dissipation):   -D(s) s, transversal zum Torus

    Stabilitätssignatur:  λ_par ≈ 0 (tangential),  λ_perp < 0 (transversal),
                          kein robustes Mode-Locking  (T^n bleibt T^n).
    "Inkommensurable Frequenzen -> echter Torus;
     rationale Einrastung -> Mode-Locking -> Dimensionsverlust T^n -> S^1."

Teil 1  Lyapunov-Spektrum, Wicklungszahl, Umlaufschnitt, Schließung
Teil 2  Arnold-Zunge (Mode-Locking im Parameterraum)
Teil 3  Kopplungsfenster k_eff(ρ, T) = R(ρ) C(T)  (Dichte taktet, Temperatur
        verschmiert oder schärft)
Teil 4  FGMS-Krümmungszyklus  Kompression -> Krümmung -> Flussänderung

Die Ausgabe wird zusätzlich nach torus_dynamik.txt (gleicher Ordner) geschrieben.
"""
import io
import os
from fractions import Fraction

import numpy as np

_BUF = io.StringIO()


def out(*args):
    text = " ".join(str(a) for a in args)
    print(text)
    _BUF.write(text + "\n")


TWO_PI = 2.0 * np.pi
GOLD = (1.0 + 5.0 ** 0.5) / 2.0

# ---------------------------------------------------------------------------
# Modellparameter (gemeinsam für alle Fälle)
# ---------------------------------------------------------------------------
C_MOD = 0.1     # c   : Empfindlichkeit der Eigenfrequenzen, ω_i(s) = ω_i0 (1 + c (s - s0))
G0 = 1.0        # g0  : konstante Einspeisung in s
G1 = 0.5        # g1  : Einspeisung aus Produkt schneller Moden, cos θ1 cos θ2 (Kanal C)
DAMP = 0.5      # D   : Dissipation/Sättigung der Strukturvariable (Kanal D)
S0 = G0 / DAMP  # s0  : Ruhelage der Strukturvariable
PHI = 0.0       # φ_ij: Verzögerungs-/Transformationswinkel


def rhs(x, p):
    """Vektorfeld auf T^2 x S, gebündelt über B Fälle. x: (B,3)."""
    th1, th2, s = x[:, 0], x[:, 1], x[:, 2]
    mod = 1.0 + C_MOD * (s - S0)
    k1 = p["k11"] * np.sin(th2 - th1 - PHI) + p["k32"] * np.sin(2 * th2 - 3 * th1 - PHI)
    k2 = p["k11"] * np.sin(th1 - th2 - PHI) + p["k32"] * np.sin(3 * th1 - 2 * th2 - PHI)
    f = np.empty_like(x)
    f[:, 0] = p["w1"] * mod + k1
    f[:, 1] = p["w2"] * mod + k2
    f[:, 2] = G0 + G1 * np.cos(th1) * np.cos(th2) - DAMP * s
    return f


def jac(x, p):
    """Jacobi-Matrix des Vektorfelds, (B,3,3)."""
    th1, th2 = x[:, 0], x[:, 1]
    k11, k32 = p["k11"], p["k32"]
    c11a = np.cos(th2 - th1 - PHI)
    c32a = np.cos(2 * th2 - 3 * th1 - PHI)
    c11b = np.cos(th1 - th2 - PHI)
    c32b = np.cos(3 * th1 - 2 * th2 - PHI)
    J = np.zeros((x.shape[0], 3, 3))
    J[:, 0, 0] = -k11 * c11a - 3 * k32 * c32a
    J[:, 0, 1] = k11 * c11a + 2 * k32 * c32a
    J[:, 0, 2] = p["w1"] * C_MOD
    J[:, 1, 0] = k11 * c11b + 3 * k32 * c32b
    J[:, 1, 1] = -k11 * c11b - 2 * k32 * c32b
    J[:, 1, 2] = p["w2"] * C_MOD
    J[:, 2, 0] = -G1 * np.sin(th1) * np.cos(th2)
    J[:, 2, 1] = -G1 * np.cos(th1) * np.sin(th2)
    J[:, 2, 2] = -DAMP
    return J


def rk4(x, p, dt):
    a = rhs(x, p)
    b = rhs(x + 0.5 * dt * a, p)
    c = rhs(x + 0.5 * dt * b, p)
    d = rhs(x + dt * c, p)
    return x + dt / 6.0 * (a + 2 * b + 2 * c + d)


def rk4_var(x, Q, p, dt):
    """RK4 für Zustand und Tangentialbasis (Variationsgleichung dQ/dt = J Q)."""
    a = rhs(x, p)
    A = jac(x, p) @ Q
    x2, Q2 = x + 0.5 * dt * a, Q + 0.5 * dt * A
    b = rhs(x2, p)
    Bm = jac(x2, p) @ Q2
    x3, Q3 = x + 0.5 * dt * b, Q + 0.5 * dt * Bm
    c = rhs(x3, p)
    Cm = jac(x3, p) @ Q3
    x4, Q4 = x + dt * c, Q + dt * Cm
    d = rhs(x4, p)
    Dm = jac(x4, p) @ Q4
    return (x + dt / 6.0 * (a + 2 * b + 2 * c + d),
            Q + dt / 6.0 * (A + 2 * Bm + 2 * Cm + Dm))


def best_rational(w, qmax=12):
    fr = Fraction(w).limit_denominator(qmax)
    return fr, abs(w - float(fr))


def circ_dist(a, b):
    d = np.mod(a - b + np.pi, TWO_PI) - np.pi
    return np.abs(d)


# ---------------------------------------------------------------------------
# Teil 1
# ---------------------------------------------------------------------------
def teil1():
    out("=" * 78)
    out("TEIL 1  Stabilitätssignatur auf M = T^2 x S")
    out("=" * 78)
    out("Parameter: c=%.2f  g0=%.2f  g1=%.2f  D=%.2f  s0=g0/D=%.2f  φ=%.2f"
        % (C_MOD, G0, G1, DAMP, S0, PHI))
    out("K_1 = k11 sin(θ2-θ1-φ) + k32 sin(2θ2-3θ1-φ)")
    out("K_2 = k11 sin(θ1-θ2-φ) + k32 sin(3θ1-2θ2-φ)")
    out("(k11: Kanal-A-Form des Autors; k32: Fourier-Anteil von K(θ,s) an der Resonanz 3:2)")
    out("")

    cases = [
        ("A", "inkommensurabel ω2/ω1 = golden", GOLD, 0.05, 0.02),
        ("B", "nahe rational 1:1, ω2 = 1.02", 1.02, 0.05, 0.00),
        ("C", "exakt rational 3:2, ungekoppelt", 1.50, 0.00, 0.00),
        ("D", "nahe rational 3:2, ω2 = 1.505", 1.505, 0.05, 0.02),
        ("E", "inkommensurabel, Kopplung zu stark", GOLD, 0.40, 0.00),
    ]
    B = len(cases)
    p = {
        "w1": np.ones(B),
        "w2": np.array([c[2] for c in cases]),
        "k11": np.array([c[3] for c in cases]),
        "k32": np.array([c[4] for c in cases]),
    }
    dt, t_trans, t_meas = 0.02, 200.0, 1500.0
    n_trans, n_meas = int(t_trans / dt), int(t_meas / dt)

    x = np.zeros((B, 3))
    x[:, 0], x[:, 1], x[:, 2] = 0.3, 1.1, S0
    for _ in range(n_trans):
        x = rk4(x, p, dt)

    Q = np.tile(np.eye(3), (B, 1, 1))
    lsum = np.zeros((B, 3))
    trace_sum = np.zeros(B)
    th_start = x[:, :2].copy()
    sec = [[] for _ in range(B)]
    s_min, s_max = x[:, 2].copy(), x[:, 2].copy()
    for _ in range(n_meas):
        x_old = x
        x, Q = rk4_var(x, Q, p, dt)
        Q, R = np.linalg.qr(Q)
        lsum += np.log(np.abs(np.diagonal(R, axis1=1, axis2=2)))
        trace_sum += np.trace(jac(x, p), axis1=1, axis2=2)
        s_min, s_max = np.minimum(s_min, x[:, 2]), np.maximum(s_max, x[:, 2])
        n0 = np.floor(x_old[:, 0] / TWO_PI)
        n1 = np.floor(x[:, 0] / TWO_PI)
        for b in np.nonzero(n1 > n0)[0]:
            target = n1[b] * TWO_PI
            fr = (target - x_old[b, 0]) / (x[b, 0] - x_old[b, 0])
            th2c = x_old[b, 1] + fr * (x[b, 1] - x_old[b, 1])
            sec[b].append(np.mod(th2c, TWO_PI))

    lyap = lsum / t_meas
    tr_mean = trace_sum / n_meas
    omega = (x[:, :2] - th_start) / t_meas

    out("Fall | Beschreibung                          | λ1        λ2        λ3       | Σλ       <tr J>")
    for b, c in enumerate(cases):
        out("  %s  | %-37s | %+.5f  %+.5f  %+.5f | %+.4f  %+.4f"
            % (c[0], c[1], lyap[b, 0], lyap[b, 1], lyap[b, 2], lyap[b].sum(), tr_mean[b]))
    out("")

    tol = 5e-3
    out("Fall | Ω1       Ω2       W=Ω2/Ω1    nächstes p/q (q<=12)  |W-p/q|     "
        "neutral  kontrah.  Dim. Invariante")
    summary = {}
    for b, c in enumerate(cases):
        W = omega[b, 1] / omega[b, 0]
        fr, err = best_rational(W)
        n_neu = int(np.sum(np.abs(lyap[b]) < tol))
        n_con = int(np.sum(lyap[b] < -tol))
        dim = {2: "T^2", 1: "S^1", 0: "Punkt"}.get(n_neu, "?")
        summary[c[0]] = (W, n_neu, n_con, dim)
        out("  %s  | %.5f  %.5f  %.8f  %-8s              %.2e    %d        %d         %s"
            % (c[0], omega[b, 0], omega[b, 1], W, str(fr), err, n_neu, n_con, dim))
    out("")

    out("Umlaufschnitt θ1 ≡ 0 (mod 2π): Werte von θ2 bei jedem Umlauf von θ1")
    out("Fall | Durchgänge  belegte Bins (von 60)  Schließungsabstand min_n |θ2(n)-θ2(0)|  "
        "Abstand letzter Umläufe")
    for b, c in enumerate(cases):
        pts = np.array(sec[b])
        bins = np.unique(np.floor(pts / TWO_PI * 60).astype(int)).size
        dmin = circ_dist(pts[1:], pts[0]).min()
        tail = pts[len(pts) // 2:]
        # Anzahl distinkter Werte in der zweiten Hälfte (Toleranz 1e-4)
        srt = np.sort(tail)
        gaps = np.diff(np.concatenate([srt, [srt[0] + TWO_PI]]))
        n_distinct = int(np.sum(gaps > 1e-4))
        out("  %s  | %5d        %3d                    %.3e                              "
            "%d distinkte Werte" % (c[0], pts.size, bins, dmin, n_distinct))
    out("")
    out("Strukturvariable s (transversal): Schwankungsbreite [min, max] nach Einschwingen")
    for b, c in enumerate(cases):
        out("  %s  | s in [%.4f, %.4f]" % (c[0], s_min[b], s_max[b]))
    out("")
    out("Analytische Kontrolle (Fall B, 1:1): ψ = θ2-θ1, dψ/dt = Δω - 2 k11 sin ψ,")
    dw = 0.02
    lam_b = -2 * 0.05 * np.sqrt(1 - (dw / (2 * 0.05)) ** 2)
    out("  Rastbedingung |Δω| <= 2 k11 = 0.100, erfüllt (Δω = 0.020);"
        "  λ_ψ = -2 k11 cos ψ* = %.5f" % lam_b)
    out("Analytische Kontrolle (Fall D, 3:2): ψ = 2θ2-3θ1, dψ/dt = Δ - 5 k32 sin ψ, Δ = 2ω2-3ω1")
    dd = 2 * 1.505 - 3
    lam_d = -5 * 0.02 * np.sqrt(1 - (dd / (5 * 0.02)) ** 2)
    out("  Rastbedingung |Δ| <= 5 k32 = 0.100, erfüllt (Δ = %.3f);  λ_ψ ≈ -5 k32 cos ψ* = %.5f"
        % (dd, lam_d))
    out("Analytische Kontrolle (Fall E): Δω = %.3f <= 2 k11 = 0.800 -> Einrastung trotz "
        "inkommensurabler Eigenfrequenzen" % (GOLD - 1))
    out("Transversale Rate: λ_perp ≈ -D = %.3f" % (-DAMP))
    out("")
    return summary


# ---------------------------------------------------------------------------
# Teil 2: Arnold-Zunge um 3:2
# ---------------------------------------------------------------------------
def winding_scan(w2, k11, k32, dt=0.02, t_trans=600.0, t_meas=1000.0):
    n = w2.size
    p = {"w1": np.ones(n), "w2": w2, "k11": np.full(n, k11), "k32": np.full(n, k32)}
    x = np.zeros((n, 3))
    x[:, 0], x[:, 1], x[:, 2] = 0.3, 1.1, S0
    for _ in range(int(t_trans / dt)):
        x = rk4(x, p, dt)
    th0 = x[:, :2].copy()
    for _ in range(int(t_meas / dt)):
        x = rk4(x, p, dt)
    om = (x[:, :2] - th0) / t_meas
    return om[:, 0], om[:, 1]


def teil2():
    out("=" * 78)
    out("TEIL 2  Mode-Locking im Parameterraum (Arnold-Zunge an 3:2, k11 = 0.05)")
    out("=" * 78)
    w2 = np.linspace(1.30, 1.75, 901)
    dw2 = w2[1] - w2[0]
    ks = [0.0, 0.01, 0.02, 0.04]
    out("Rastbreite analytisch (ohne s-Modulation und ohne k11-Anteil): |2ω2 - 3| <= 5 k32"
        "  ->  Breite in ω2 = 5 k32")
    out("k32    | gerastete ω2-Breite (3:2)  analytisch 5*k32  | Anteil gerastet im Fenster "
        "[1.30,1.75] | ω2=golden gerastet?")
    stairs = None
    for k32 in ks:
        o1, o2 = winding_scan(w2, 0.05, k32)
        locked = np.abs(2 * o2 - 3 * o1) < 1e-3
        width = locked.sum() * dw2
        og1, og2 = winding_scan(np.array([GOLD]), 0.05, k32)
        Wg = og2[0] / og1[0]
        frg, errg = best_rational(Wg)
        g_locked = "ja" if errg < 1e-4 else "nein (W=%.6f)" % Wg
        out("%.3f  | %.4f                     %.4f             | %.3f                              "
            "| %s" % (k32, width, 5 * k32, locked.mean(), g_locked))
        if k32 == 0.02:
            stairs = (w2, o2 / o1)
    out("")
    out("Wicklungszahl W(ω2) für k32 = 0.02 (Treppenstufe = Rastplateau):")
    w2s, Ws = stairs
    for target in [1.40, 1.44, 1.46, 1.48, 1.50, 1.52, 1.54, 1.56, 1.60, GOLD, 1.70]:
        i = int(np.argmin(np.abs(w2s - target)))
        out("   ω2 = %.4f   W = %.6f" % (w2s[i], Ws[i]))
    out("")


# ---------------------------------------------------------------------------
# Teil 3: Kopplungsfenster im (ρ, T)-Raum
# ---------------------------------------------------------------------------
def teil3():
    out("=" * 78)
    out("TEIL 3  Kopplungsfenster k_eff(ρ,T) = k0 R(ρ) C(T)")
    out("=" * 78)
    k0, rho_s, n_R, T_a, T_d = 0.3, 1.0, 2, 0.3, 1.5
    dw, k_min = 0.10, 0.01
    out("Rechenannahmen (Parametrisierung der Autoraussagen):")
    out("  R(ρ) = ρ^%d/(ρ^%d + ρ*^%d), ρ* = %.1f          (monoton, sigmoidal, gesättigt)"
        % (n_R, n_R, n_R, rho_s))
    out("  C(T) = (1 - exp(-T/T_a)) exp(-T/T_d), T_a = %.1f, T_d = %.1f   (Maximum bei T*)"
        % (T_a, T_d))
    Tstar = T_a * np.log(1 + T_d / T_a)
    out("  T* = T_a ln(1 + T_d/T_a) = %.3f" % Tstar)
    out("  1:1-Paar mit Δω = %.2f  ->  Verriegelung für k_eff >= Δω/2 = %.3f" % (dw, dw / 2))
    out("  Querbegrenzung bricht weg für k_eff < k_min = %.3f" % k_min)
    rhos = np.logspace(-1, 1, 13)
    Ts = np.logspace(np.log10(0.03), np.log10(6.0), 13)
    RR = rhos ** n_R / (rhos ** n_R + rho_s ** n_R)
    CC = (1 - np.exp(-Ts / T_a)) * np.exp(-Ts / T_d)
    keff = k0 * np.outer(CC, RR)  # Zeilen: T, Spalten: ρ
    out("")
    out("Karte (Zeilen T, Spalten ρ):  V = verriegelt (T^2->S^1),  o = moduliert (Fenster),  "
        ". = zerfasert")
    out("   T \\ ρ  " + " ".join("%5.2f" % r for r in rhos))
    for i, T in enumerate(Ts):
        row = []
        for j in range(rhos.size):
            k = keff[i, j]
            row.append("  V  " if k >= dw / 2 else ("  o  " if k >= k_min else "  .  "))
        out("  %6.3f  " % T + " ".join(row))
    # numerische Kontrolle der Verriegelungsgrenze
    kk = keff.ravel()
    n = kk.size
    p = {"w1": np.ones(n), "w2": np.full(n, 1 + dw), "k11": kk, "k32": np.zeros(n)}
    x = np.zeros((n, 3))
    x[:, 0], x[:, 1], x[:, 2] = 0.3, 1.1, S0
    dt = 0.02
    for _ in range(int(300 / dt)):
        x = rk4(x, p, dt)
    th0 = x[:, :2].copy()
    for _ in range(int(600 / dt)):
        x = rk4(x, p, dt)
    om = (x[:, :2] - th0) / 600.0
    num_locked = np.abs(om[:, 1] - om[:, 0]) < 1e-3
    ana_locked = kk >= dw / 2
    agree = np.mean(num_locked == ana_locked)
    out("")
    out("Numerische Kontrolle der Grenze k_eff = Δω/2 (Simulation aller %d Gitterpunkte): "
        "Übereinstimmung %.1f %%" % (n, 100 * agree))
    mism = np.nonzero(num_locked != ana_locked)[0]
    if mism.size:
        out("  Abweichende Punkte (k_eff): " + ", ".join("%.4f" % kk[m] for m in mism[:8]))
    out("")


# ---------------------------------------------------------------------------
# Teil 4: FGMS-Krümmungszyklus
# ---------------------------------------------------------------------------
def row_softmax(L):
    L = L - L.max(axis=1, keepdims=True)
    E = np.exp(L)
    return E / E.sum(axis=1, keepdims=True)


def fgms_run(zeta, n=12, steps=400, eta=0.05, beta=1.0, alpha_R=2.0, eta_b=0.05,
             K_min=0.05, K_max=20.0, seed=7, log_at=(0, 25, 50, 100, 200, 400)):
    rng = np.random.default_rng(seed)
    eps = 1e-6
    K = np.ones((n, n))
    m = np.full(n, 1.0 / n) * (1 + 0.05 * rng.standard_normal(n))
    m /= m.sum()
    pi = np.full((n, n), 1.0 / n)
    F = m[:, None] * pi
    Fbar = F.copy()
    mu_C = np.mean(np.log(m + eps))
    rows = []

    def metrics(t):
        Kn = K / K.sum(axis=1, keepdims=True)
        H = -(Kn * np.log(Kn)).sum(axis=1).mean() / np.log(n)
        lam = np.linalg.eigvalsh(np.log(K + eps)).max()
        chi = F.sum(axis=1).max()
        capped = np.mean((K >= K_max - 1e-9) | (K <= K_min + 1e-9))
        rows.append((t, chi, H, lam, K.max(), capped))

    metrics(0)
    for t in range(1, steps + 1):
        # Kompression C_t(i) aus Knotendurchsatz m_t(i)
        logm = np.log(m + eps)
        mu_C = 0.9 * mu_C + 0.1 * logm.mean()
        C = np.clip(logm - mu_C, -3.0, 3.0)
        # Diffusion/Entkrümmung D_t = L_K C
        D = C - (K @ C) / K.sum(axis=1)
        # multiplikatives Krümmungsupdate, harte Kappen
        K = K * np.exp(eta * (C[:, None] + C[None, :]) - zeta * (D[:, None] + D[None, :]))
        K = np.clip(K, K_min, K_max)
        g = -np.log(K + eps)
        # Relationaler Impuls R_t = α_R (F - F̄), Flussprognose RowSoftmax(R - β g)
        R = alpha_R * n * n * (F - Fbar)
        pi = row_softmax(R - beta * g)
        m = m @ pi
        m = m * (1 + 0.002 * rng.standard_normal(n))   # kleine externe Störung ΔF_ext
        m = np.abs(m) / np.abs(m).sum()
        F = m[:, None] * pi
        Fbar = (1 - eta_b) * Fbar + eta_b * F
        if t in log_at:
            metrics(t)
    return rows


def teil4():
    out("=" * 78)
    out("TEIL 4  FGMS-Krümmungszyklus: Kompression -> Krümmung -> Flussänderung")
    out("=" * 78)
    out("n = 12 Regionen, η = 0.05, β = 1, K in [0.05, 20], g = -log K, R = 0,")
    out("π = RowSoftmax(-βg), m_{t+1} = m_t π, C = log m - μ_C, D = L_K C,")
    out("K <- clip(K exp(η(C_i+C_j) - ζ(D_i+D_j)))")
    out("χ = max_i Σ_j F_ij (Zentrierung), H_K = mittlere normierte Zeilenentropie von K,")
    out("λ_max = größter Eigenwert von log K, 'gekappt' = Anteil der K-Einträge an K_min/K_max")
    eta = 0.05
    for zeta in [0.0, 0.03, 0.045, 0.05, 0.07]:
        out("")
        out("ζ = %.3f   (ζ/η = %.2f)" % (zeta, zeta / eta))
        out("    t    χ        H_K      λ_max     K_max     gekappt")
        rows = fgms_run(zeta, eta=eta, alpha_R=0.0, steps=1000,
                        log_at=(100, 200, 400, 600, 1000))
        for (t, chi, H, lam, kmax, cap) in rows:
            out("  %4d   %.4f   %.4f   %7.3f   %7.3f   %.3f" % (t, chi, H, lam, kmax, cap))
    out("")
    out("Gleichverteilung: χ = 1/12 = %.4f" % (1 / 12))
    out("Linearisierung um K = 1, m = 1/n: D_i = C_i - <C> ≈ C_i, also")
    out("  δ log K_ij ≈ (η - ζ)(C_i + C_j) pro Schritt  ->  kritische Linie ζ = η")
    out("")


def main():
    out("Torus-Dynamik-Schicht (FTM) - Rechenbeispiel")
    out("numpy %s" % np.__version__)
    out("")
    teil1()
    teil2()
    teil3()
    teil4()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "torus_dynamik.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_BUF.getvalue())


if __name__ == "__main__":
    main()
