"""Eine innere Richtung je Falte oder eine gemeinsame?

(A) Bündel: Raum T³ und je Masse eine eigene Halbgerade [0, Λ], bei Λ geschlossen (reflektierend).
(B1) Stern: N Fasern der Länge Λ, alle bei ϖ = Λ an einem gemeinsamen Knoten verbunden.
(B3) T³ × [0, Λ] mit seitlicher Bewegung (harmonische Fortsetzung (∇² + ∂_ϖ²), isotrop): das Ende bei Λ ist
     für alle Orte gemeinsam, weil dort die seitlichen Moden abgeklungen sind.
Auf jeder Faser e: D_e·u'' = Γ_e·u (Diffusion ohne Drift), D_e = GM_e/c, Γ_e = a_N,e/c, Start in der Tiefe ϖ_R,e = r*_e.
Knoten: u stetig, Fluss Σ D_e·u_e' = 0 (Kette mit c-Schritten h_e = r₀,e, D_e = h_e·c/2).
q_ij = E[e^{−∫Γ dt}; Ankunft am Ort j], Start auf Faser i.

Blickrichtung: Γ lokal (L, Gesamtfeld am Ort, heute); D, ϖ_R, Λ, Knoten, ⟨T⟩ innen (I); ΣΔs_b, q_kosm, d_LSS
global (G, heute, für unseren Horizont); q_ii, ν am Ort gelesen (L, Beobachter an der Bahn).
Rechnet nur mit ftm_kern. Ausgabe: gemeinsame_innere_richtung.txt
"""
import math
import numpy as np
import mpmath as mp
from scipy.integrate import quad
from scipy.optimize import brentq

import ftm_kern as K

mp.mp.dps = 40
out = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def lg(x):
    x = mp.mpf(x)
    return float(mp.log10(x)) if x > 0 else float("-inf")


# ------------------------------------------------------------------ Stern-Graph, analytisch
def _leit(D, k, l):
    """D·k·coth(k·l) bzw. D/l für k = 0."""
    return D / l if k == 0 else D * k * mp.coth(k * l)


def _quer(D, k, l):
    """D·k/sinh(k·l) bzw. D/l für k = 0."""
    return D / l if k == 0 else D * k / mp.sinh(k * l)


def _prof(k, l, x, a, b):
    """Lösung auf [0, l] mit u(0) = a, u(l) = b (k = 0: linear)."""
    if k == 0:
        return a * (1 - x / l) + b * x / l
    return (a * mp.sinh(k * (l - x)) + b * mp.sinh(k * x)) / mp.sinh(k * l)


def stern(D, k, l, i, L):
    """q_ij für alle j, Start auf Faser i in der Tiefe L. D, k, l: Listen je Faser (mp)."""
    S = sum(_leit(D[e], k[e], l[e]) for e in range(len(D)))
    q = []
    for j in range(len(D)):
        U = _quer(D[j], k[j], l[j]) / S
        q.append(_prof(k[i], l[i], L, mp.mpf(1 if i == j else 0), U))
    return q


def halbgerade_A(k, Lam, L):
    """(A) bei Λ reflektierend: cosh(k(Λ − L))/cosh(kΛ)."""
    return mp.cosh(k * (Lam - L)) / mp.cosh(k * Lam)


# ------------------------------------------------------------------ Größen
M_MW = 6e10 * K.M_sun
M_M31 = 1e11 * K.M_sun               # Beispielmasse eines Nachbarn (gesetzt)
Lam = K.schliesslaenge(K.A0)         # I: Λ = c²/(2a₀)
rs_MW = K.r_stern(M_MW)
p("Λ = c²/(2a₀) = %.2f Gpc (I, heute fest);  r*_MW = %.2f kpc;  r*/Λ = %.2e" % (Lam / K.Gpc, rs_MW / K.kpc, rs_MW / Lam))
p("Punktmasse: k = √(Γ/D) = √(a_N/GM) = 1/r (Wellenzahl der Kreisbahnmode), k·ϖ_R = r*/r = s")
p()

# ======================================================================================
p("=" * 100)
p("1  Rückkehrstatistik je Topologie")
p("=" * 100)
p("(A)  q_i = cosh(k(Λ − ϖ_R))/cosh(kΛ); die Fasern sind unabhängig, jeder Rücklauf kommt am eigenen Ort an.")
p("(B1) Stern: u_e(ϖ) = [δ_ej·sinh(k_e(Λ−ϖ)) + U_j·sinh(k_eϖ)]/sinh(k_eΛ),")
p("     U_j = D_jk_j/sinh(k_jΛ) / Σ_e D_ek_e·coth(k_eΛ)   (Wert am Knoten für Ankunft bei j)")
p("     q_ii = sinh(k_i(Λ−ϖ_R))/sinh(k_iΛ) + U_i·sinh(k_iϖ_R)/sinh(k_iΛ),   q_ij = U_j·sinh(k_iϖ_R)/sinh(k_iΛ)")
p("     Summe Σ_j q_ij: Rückkehr irgendwohin; q_ii: Rückkehr an den eigenen Ort.")
p()

# 1a Grenzfall N = 1 und N = 2 analytisch
p("1a) Prüfung der Formel (unabhängige Wege)")
kk, LL, ll = mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("4.0")
q1 = stern([mp.mpf(1)], [kk], [ll], 0, LL)[0]
p("   N = 1:  Stern %.15f   Halbgerade cosh(k(Λ−L))/cosh(kΛ) %.15f" % (q1, halbgerade_A(kk, ll, LL)))
q2 = stern([mp.mpf(1)] * 2, [kk] * 2, [ll] * 2, 0, LL)
p("   N = 2 gleich = Strecke [0, 2Λ], beide Enden Orte:  q_11 = sinh(k(2Λ−L))/sinh(2kΛ), q_12 = sinh(kL)/sinh(2kΛ)")
p("        Stern q_11 %.15f  Strecke %.15f" % (q2[0], mp.sinh(kk * (2 * ll - LL)) / mp.sinh(2 * kk * ll)))
p("        Stern q_12 %.15f  Strecke %.15f" % (q2[1], mp.sinh(kk * LL) / mp.sinh(2 * kk * ll)))


# Differenzenverfahren auf dem Stern (unabhängig vom analytischen Weg)
def stern_fd(D, k, l, i, L, j, n=4000):
    N = len(D)
    h = [l[e] / n for e in range(N)]
    idx = lambda e, m: e * (n - 1) + (m - 1)    # innere Punkte m = 1 … n−1, Knoten zuletzt
    nn = N * (n - 1) + 1
    from scipy.sparse import lil_matrix
    A = lil_matrix((nn, nn)); b = np.zeros(nn); kn = nn - 1
    for e in range(N):
        for m in range(1, n):
            r = idx(e, m)
            A[r, r] = -2 * D[e] / h[e]**2 - D[e] * k[e]**2
            if m > 1:
                A[r, idx(e, m - 1)] = D[e] / h[e]**2
            else:
                b[r] -= D[e] / h[e]**2 * (1.0 if e == j else 0.0)
            if m < n - 1:
                A[r, idx(e, m + 1)] = D[e] / h[e]**2
            else:
                A[r, kn] = D[e] / h[e]**2
    for e in range(N):                          # Fluss am Knoten, 2. Ordnung
        A[kn, kn] += D[e] * 3 / (2 * h[e]); A[kn, idx(e, n - 1)] -= D[e] * 4 / (2 * h[e]); A[kn, idx(e, n - 2)] += D[e] / (2 * h[e])
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import spsolve
    u = spsolve(A.tocsr(), b)
    m = int(round(float(L / l[i]) * n))
    return u[idx(i, m)]


Dl = [1.0, 2.5, 0.4]; kl = [0.7, 0.2, 1.5]; ll3 = [4.0, 4.0, 4.0]; L3 = 1.3
an = stern([mp.mpf(x) for x in Dl], [mp.mpf(x) for x in kl], [mp.mpf(x) for x in ll3], 0, mp.mpf(L3))
p("   N = 3 verschieden (D = 1; 2,5; 0,4, k = 0,7; 0,2; 1,5, Λ = 4, L = 1,3), analytisch gegen Differenzen:")
for j in range(3):
    p("        q_1%d  %.8f   %.8f" % (j + 1, an[j], stern_fd(Dl, kl, ll3, 0, L3, j)))


# Monte-Carlo der Kette mit c-Schritten (dritter Weg), klein
def stern_mc(m_e, surv, start, n_lauf=200000, seed=1):
    """Kette: Faser e hat m_e Schritte, Überleben je Schritt surv_e, Knoten wählt Faser gleich wahrscheinlich
    (bei Schrittlängen h_e ∝ D_e ist das der Fluss Σ D_e u' = 0 nur bei gleichen h; hier gleiche h)."""
    rng = np.random.default_rng(seed)
    N = len(m_e); treff = np.zeros(N); e0, pos0 = start
    for _ in range(n_lauf):
        e, pos, w = e0, pos0, 1.0
        while True:
            if pos == m_e[e]:                    # Knoten
                e = rng.integers(N); pos = m_e[e] - 1; w *= surv[e]; continue
            pos += 1 if rng.random() < 0.5 else -1
            w *= surv[e]
            if pos == 0:
                treff[e] += w; break
    return treff / n_lauf


def stern_kette_exakt(m_e, surv, start):
    """Dieselbe Kette als lineares Gleichungssystem (exakt)."""
    N = len(m_e); nn = sum(m - 1 for m in m_e) + 1; kn = nn - 1
    off = np.cumsum([0] + [m - 1 for m in m_e[:-1]])
    idx = lambda e, x: off[e] + x - 1
    res = []
    for j in range(N):
        A = np.eye(nn); b = np.zeros(nn)
        for e in range(N):
            for x in range(1, m_e[e]):
                r = idx(e, x)
                for y in (x - 1, x + 1):
                    if y == 0:
                        b[r] += 0.5 * surv[e] * (1.0 if e == j else 0.0)
                    elif y == m_e[e]:
                        A[r, kn] -= 0.5 * surv[e]
                    else:
                        A[r, idx(e, y)] -= 0.5 * surv[e]
        for e in range(N):
            A[kn, idx(e, m_e[e] - 1)] -= surv[e] / N
        u = np.linalg.solve(A, b)
        res.append(u[idx(start[0], start[1])])
    return np.array(res)


m_e = [12, 12, 12]; surv = [0.995, 0.999, 0.97]
p("   Kette (3 Fasern à 12 Schritte, Überleben je Schritt 0,995 / 0,999 / 0,97, Start 3 Schritte tief auf Faser 1):")
p("        exakt  " + "  ".join("%.5f" % x for x in stern_kette_exakt(m_e, surv, (0, 3))))
p("        Monte-Carlo (2·10⁵)  " + "  ".join("%.5f" % x for x in stern_mc(m_e, surv, (0, 3))))
p()

# 1b Galaxie in (B1): örtliche Regel?
p("1b) Milchstraße (6·10¹⁰ M☉) in (A) und (B1); Nachbarfasern: M31-artig (10¹¹ M☉ bei ihrem r*) und N_leer Fasern mit Γ = 0")
p("    ln(q/e^{−s}) als log10 des Betrags; e^{−2kΛ} = e^{−2Λ/r} setzt die Größe")
p("    r [kpc]    s      log10|ln(q_A/e^{−s})|   log10|ln(q_ii,B1/e^{−s})| (N_leer = 10⁶)   log10 q_ij (MW→M31)")
rs31 = K.r_stern(M_M31)
for r in (rs_MW, 3 * rs_MW, 50 * K.kpc):
    s = rs_MW / r
    k_i = mp.mpf(1) / r
    D_i = mp.mpf(K.G * M_MW / K.c)
    D_j = mp.mpf(K.G * M_M31 / K.c); k_j = mp.mpf(1) / rs31
    n_leer = 10**6
    # N_leer gleiche Fasern mit Γ = 0: in der Knotensumme als n_leer·D/Λ
    S = _leit(D_i, k_i, Lam) + _leit(D_j, k_j, Lam) + n_leer * D_i / Lam
    U_i = _quer(D_i, k_i, Lam) / S; U_j = _quer(D_j, k_j, Lam) / S
    L = mp.mpf(rs_MW)
    # Abweichungen sind e^{−10⁵}-klein: als Logarithmen der Einzelterme rechnen, nicht als Differenz zu 1
    # (A):  ln(q_A) + s = ln[(1 + e^{−2k(Λ−L)})/(1 + e^{−2kΛ})] ≈ e^{−2k(Λ−L)}
    lA = float(-2 * k_i * (Lam - L) / mp.log(10))
    # (B1): q_ii/e^{−s} − 1 = −e^{−2k(Λ−L)}(…) + U_i·sinh(kL)/(sinh(kΛ)·e^{−s})
    t2 = U_i * mp.sinh(k_i * L) / mp.sinh(k_i * Lam) / mp.exp(-s)
    lB = max(lA, lg(t2))
    qij = _prof(k_i, Lam, L, 0, U_j)
    p("    %6.1f   %5.3f        %12.1f                %12.1f                        %12.1f" % (
        r / K.kpc, s, lA, lB, lg(qij)))
p("    → (B1) gibt am Ort q = e^{−s}; Abweichung und Übertritt ∝ e^{−2Λ/r}, Bedingung Λ ≫ r (nicht nur Λ ≫ ϖ_R).")
p()

# 1c Wo weicht es ab? r₁₀ in A und im Stern mit vielen leeren Fasern
p("1c) 10-%-Radius: |ν/ν_örtlich − 1| = 0,1")


def nu_loc(s):
    return 1 / (1 - math.exp(-s))


def r10(M, art):
    rs = K.r_stern(M)
    def f(x):                                  # x = r/Λ
        r = x * Lam; k = mp.mpf(1) / r; L = mp.mpf(rs)
        if art == "A":
            q = halbgerade_A(k, Lam, L)
        else:                                   # N_leer → ∞: U_i → 0
            q = mp.sinh(k * (Lam - L)) / mp.sinh(k * Lam)
        return abs(float((1 / (1 - q)) / nu_loc(rs / r)) - 1) - 0.1
    return brentq(f, 0.05, 5)


for M, name in ((M_MW, "Milchstraße"), (1e8 * K.M_sun, "Zwerg 10⁸ M☉")):
    p("    %-14s  (A) r₁₀ = %.4f·Λ (ν größer)   (B1, N_leer → ∞) r₁₀ = %.4f·Λ (ν kleiner)" % (name, r10(M, "A"), r10(M, "B")))
p("    Handwert tief: ν/ν_örtlich = coth(Λ/r) bzw. tanh(Λ/r) ⇒ r₁₀ = Λ/artanh(1/1,1) = %.4f·Λ bzw. Λ/artanh(0,9) = %.4f·Λ"
  % (1 / math.atanh(1 / 1.1), 1 / math.atanh(0.9)))
p("    Weit außen (r ≫ Λ): (A) 1 − q → s²(2Λ/ϖ_R − 1)/2 … a → a₀-Boden; (B1, N_leer → ∞) 1 − q → ϖ_R/Λ, ν → Λ/r* = %.2e (Newton mal Konstante)"
  % (Lam / rs_MW))
p()

# ======================================================================================
p("=" * 100)
p("2  Homogene Verteilung, Netzwerk, kosmisches q")
p("=" * 100)
# 2a gleichartige Fasern: Summe = (A)
for N in (1, 2, 10, 1000):
    q = stern([mp.mpf(1)] * N, [kk] * N, [ll] * N, 0, LL)
    p("   N = %4d gleich:  Σ_j q_ij = %.12f   (A) %.12f   q_ii = %.6f" % (N, sum(q), halbgerade_A(kk, ll, LL), q[0]))
p("   → Bei gleichen Fasern ist der Knoten für die Summe ein Spiegel: Σ_j q_ij = q_A für jedes N (analytisch: U_tot = 1/cosh kΛ).")
p("   Gesamtfeld, nicht additiv: In (A) und (B1) kommt es nur über Γ = |g_ges|/c am Ort herein; die Topologie erzeugt es nicht.")
p("   In (B3) ist der seitliche Ankunftsort die Poisson-Glättung des Modells (Abschnitt 4): Rückläufe aus der Tiefe ϖ mitteln")
p("   über ~ϖ seitlich, das ist der kollektive Trichter, der in die Nachbarn übergeht (Netzwerk).")
p()
# 2b Γ = 0: Rückkehr nach Hause
p("2b) Γ = 0 (homogener Hintergrund): (A) q = 1, ν = ∞.  (B1) q_ii = 1 − (L/Λ)(1 − 1/N): Ankunft woanders ohne Verlust.")
for N in (2, 10, 10**6):
    x = mp.mpf("0.157")
    q = stern([mp.mpf(1)] * 2, [mp.mpf(0)] * 2, [mp.mpf(1)] * 2, 0, x) if N == 2 else None
    qf = 1 - x * (1 - mp.mpf(1) / N)
    if q is not None:
        p("   N = 2, L/Λ = 0,157: Stern q_11 = %.6f, Formel %.6f" % (q[0], qf))
    else:
        p("   N = %d, L/Λ = 0,157: q_ii = %.6f" % (N, qf))
p("   Gambler's ruin = Kac: 1 − q_ii (N → ∞) = L/Λ = Maß des Ortsstücks [0, L] / Maß der Faser.")
p()
# 2c mit Verlust aus dem linearen Feld
p("2c) Mit Verlust (s = k·L) und N → ∞: q_ii = sinh(k(Λ−L))/sinh(kΛ). Verlangtes L/Λ für q_ii = q_kosm = %.4f:" % K.q_KOSM)
for s in (0.0, 0.05, 0.1, 0.3):
    if s == 0:
        x = 1 - K.q_KOSM
    elif math.exp(-s) <= K.q_KOSM:
        p("     s = %.2f:  keine Lösung (q_ii ≤ e^{−s} = %.3f < q_kosm; Grenze s < −ln q_kosm = %.3f)" % (s, math.exp(-s), -math.log(K.q_KOSM)))
        continue
    else:
        x = brentq(lambda x: float(mp.sinh(s / x * (1 - x)) / mp.sinh(s / x)) - K.q_KOSM, 1e-3, 0.999)
    p("     s = %.2f:  L/Λ = %.4f" % (s, x))
p("   (s² = a_N/a₀ ≈ 0,01 im linearen Feld heute, Notiz 69: s ≈ 0,1)")
p()
# 2d kosmische Zahlen
rho_c = 3 * K.H0_si()**2 / (8 * math.pi * K.G)
for name, Lz in (("Zelle = d_LSS (Annahme, gerechnet)", K.d_LSS()), ("Zelle aus a₀ (zelle_aus_a0, benutzt ν_kosm)", K.zelle_aus_a0(K.A0))):
    Mb = K.OMEGA_B * rho_c * Lz**3
    Sds = K.ds(Mb)                         # G: ΣΔs_b der Zelle
    x = Sds / Lam
    p("   %s: L_Zelle = %.2f Gpc, ΣΔs_b = %.3f Gpc, ΣΔs_b/Λ = %.4f ⇒ q = %.4f, ν = Λ/ΣΔs_b = %.3f  (gemessen ν_kosm = %.3f)"
      % (name, Lz / K.Gpc, Sds / K.Gpc, x, 1 - x, 1 / x, 1 + K.Q_KOSM))
x_lss = K.ds(K.OMEGA_B * rho_c * K.d_LSS()**3) / Lam
p("   Andere Startiefen des Kosmos: L = r*_Zelle = √(6ΣΔs_b·Λ) ⇒ 1 − q = √(6·ΣΔs_b/Λ) = %.3f;  L = r₀ = 6ΣΔs_b ⇒ 1 − q = %.3f"
  % (math.sqrt(6 * x_lss), 6 * x_lss))
p("   → Nur L = ΣΔs_b (die Faltung selbst) gibt q_kosm; das ist Λ = ν_kosm·ΣΔs_b umgeschrieben, keine Herleitung von ν.")
p()

# ======================================================================================
p("=" * 100)
p("3  Kac-Budget")
p("=" * 100)


def kac_stern(m_e, h_e, i):
    """Mittlere Rückkehrzeit zum Blatt i der Kette (Schrittzeit h_e/c auf Faser e, c = 1), lineares System."""
    N = len(m_e); nn = sum(m - 1 for m in m_e) + 1 + N; kn = nn - 1
    off = np.cumsum([0] + [m for m in m_e[:-1]])
    idx = lambda e, x: off[e] + x              # x = 0 … m_e−1, Blatt x = 0
    A = np.eye(nn); b = np.zeros(nn)
    for e in range(N):
        for x in range(0, m_e[e]):
            r = idx(e, x)
            nb = [x + 1] if x == 0 else [x - 1, x + 1]
            b[r] = h_e[e]
            for y in nb:
                w = 1.0 / len(nb)
                if e == i and y == 0 and x != 0:
                    continue
                if y == m_e[e]:
                    A[r, kn] -= w
                else:
                    A[r, idx(e, y)] -= w
    b[kn] = np.mean(h_e)
    for e in range(N):
        A[kn, idx(e, m_e[e] - 1)] -= 1.0 / N
    T = np.linalg.solve(A, b)
    return T[idx(i, 0)]


m_e = [10, 20, 40]; h_e = [3.0, 1.5, 0.75]           # alle Fasern Länge 30
p("   Kette mit c-Schritten, 3 Fasern (m = 10/20/40 Schritte, h = 3/1,5/0,75, Länge je 30):")
for i in range(3):
    p("     ⟨T_%d⟩ = %.4f   2·Σℓ_e/c = %.4f" % (i + 1, kac_stern(m_e, h_e, i), 2 * 90.0))
p("   Die Schrittzeit h/c macht das Maß zur Länge: ⟨T_i⟩ = 2·(Gesamtlänge des Netzes)/c, für jeden Ort gleich, unabhängig von der Masse.")
p("   (A): ⟨T⟩ = 2Λ/c = τ_D = c/a₀ (Modell).  (B1): ⟨T⟩ = 2NΛ/c = N·τ_D.")
p("   ⇒ Soll in (B) τ_D = ⟨T⟩ bleiben, ist Λ = c²/(2a₀) die Gesamtlänge des gemeinsamen Netzes, nicht die Tiefe jeder Faser.")
p()
p("   Budget: Σ_i μ(Ort_i)/μ(X) ≤ 1.  Mit μ(Ort_i) = Δs_i, μ(X) = Λ: ΣΔs_b ≤ Λ (nur Ungleichung).")
p("     Zahlen (Zelle = d_LSS): ΣΔs_b/Λ = %.4f (Ortsmaß Δs), %.4f (Ortsmaß r₀/2 = 3Δs, Kette mit c-Schritten)" % (x_lss, 3 * x_lss))
p("     Gleichheit mit dem Faktor ν, Λ = ν·ΣΔs_b, verlangt zusätzlich 1 − q_kosm = μ(Orte)/μ(X) (Kac-Lesung des Modells); das Budget selbst liefert sie nicht.")
p()
p("   Eigene Tiefe: Die örtliche Regel verlangt je Ort eine Faser, die bis ϖ_R = r* nicht mit anderen geteilt ist.")
p("   r*_i = √(6Δs_i·Λ) (r*² = r₀Λ) ⇒ Σ r*_i = √(6Λ)·Σ√Δs_i ≥ √(6Λ·ΣΔs_b) = Λ·√(6ΣΔs_b/Λ)")
Mb = K.OMEGA_B * rho_c * K.d_LSS()**3
N_MW = Mb / M_MW
p("     eine Masse mit allen Baryonen der Zelle: Σ r* = %.3f·Λ;  in Stücken der Milchstraße (N = %.2e): Σ r* = %.2e·Λ"
  % (math.sqrt(6 * x_lss), N_MW, math.sqrt(N_MW) * math.sqrt(6 * x_lss)))
p("     Handwert: N·r*_MW/Λ = %.2e" % (N_MW * rs_MW / Lam))
p("   ⇒ Ein gemeinsames Netz der Gesamtlänge Λ kann die eigenen Tiefen nicht tragen (Faktor ~10⁶ zu kurz).")
p("     (B1) trägt die Regel, aber dann ist Λ Tiefe je Faser und ⟨T⟩ = N·τ_D; das Budget je Ort ist dasselbe wie in (A).")
p()

# ======================================================================================
p("=" * 100)
p("4  Knoten und Korrelation zweier Galaxien")
p("=" * 100)
d_MW_M31 = 0.78 * K.Mpc          # [G] Abstand Milchstraße–M31
p("   (B1) Übertritt MW → M31 (1b): q_ij ~ 10^(−2Λ/r·log10 e) ⇒ bei 50 kpc log10 q_ij ≈ %.0f." % (-(Lam / (50 * K.kpc) + Lam / rs31) / math.log(10)))
p("   Obere Grenze ohne Verlust (Γ = 0 auf der Faser): Anteil woanders ≤ ϖ_R/Λ = r*/Λ = %.1e; δν/ν ≤ ν·r*/Λ = %.1e bei 50 kpc (ν = %.2f)"
  % (rs_MW / Lam, K.nu(K.G * M_MW / (50 * K.kpc)**2) * rs_MW / Lam, K.nu(K.G * M_MW / (50 * K.kpc)**2)))
p("   Messgrenze Rotationskurven einige Prozent ⇒ (B1) um ≥ 10⁴ darunter, mit Verlust um e^{−10⁵}.")
p("   Der gemeinsame Knoten ist für alle Paare derselbe: Er unterscheidet kein verschränktes Paar. Die Verschränkung von 04")
p("   (d_innen = 0) ist ein Knoten in der Tiefe 0 für das Paar, nicht der Knoten bei Λ.")
p()
p("   (B3) T³ × [0, Λ] isotrop: ϖ-Bewegung ist die 1D-Halbgerade (Produkt), also q bei festem Γ und ⟨T⟩ = 2Λ/c wie (A).")
p("   Neu ist nur der seitliche Ankunftsort. Ohne Verlust ist er der Poisson-Kern des Modells (Glättung e^{−|k|ϖ}):")
p("   P(ρ) ∝ ϖ/(ρ² + ϖ²)², Anteil innerhalb u = ρ/ϖ: F(u) = (2/π)[arctan u − u/(1+u²)]")
Fnum = lambda u: quad(lambda x: x**2 / (1 + x**2)**2, 0, u)[0] / quad(lambda x: x**2 / (1 + x**2)**2, 0, np.inf)[0]
Fan = lambda u: 2 / math.pi * (math.atan(u) - u / (1 + u * u))
p("     F(1) = %.4f (Integral %.4f);  Median u = %.3f;  Schwanz 1 − F ≈ 4/(πu)" % (Fan(1), Fnum(1), brentq(lambda u: Fan(u) - 0.5, 0.1, 10)))
u = d_MW_M31 / rs_MW
p("     Rücklauf aus ϖ_R = r*_MW, Ankunft jenseits M31 (d = 0,78 Mpc [G], u = %.0f): ohne Verlust %.2e (Formel %.2e)"
  % (u, 1 - Fan(u), 4 / (math.pi * u)))


# mit Verlust: Γ fest (am Ort), seitliches D = D (Annahme); Dichte erster Ankunft in ϖ, seitlich Gauß 3D
def ankunft_mit_verlust(s, u):
    """Anteil (an den überlebenden Rückläufen) mit seitlichem Abstand > u·ϖ_R. Einheiten ϖ_R = 1, τ_D = 1 (D = 1)."""
    Gam = s * s
    f = lambda t: 1 / math.sqrt(4 * math.pi * t**3) * math.exp(-1 / (4 * t))
    from scipy.special import erfc
    def tail3d(t):                              # P(|X| > u), X ~ N(0, 2t·I₃)
        z = u / math.sqrt(4 * t)
        return erfc(z) + 2 * z / math.sqrt(math.pi) * math.exp(-z * z)
    num = quad(lambda t: f(t) * math.exp(-Gam * t) * tail3d(t), 0, np.inf, limit=400)[0]
    den = quad(lambda t: f(t) * math.exp(-Gam * t), 0, np.inf, limit=400)[0]
    return num, den


num, den = ankunft_mit_verlust(1.0, 1.0)
p("     mit Verlust, s = 1 (Γ fest, seitlich gleiches D, Annahme): Überleben %.4f (e^{−1} = %.4f), Anteil jenseits 1·ϖ_R %.3f" % (den, math.exp(-1), num / den))
num0, den0 = ankunft_mit_verlust(1e-6, 1.0)
p("     Gegenprobe Γ → 0: Anteil jenseits 1·ϖ_R %.4f (Poisson-Kern 1 − F(1) = %.4f)" % (num0 / den0, 1 - Fan(1)))
s50 = rs_MW / (50 * K.kpc)
for uu in (10.0, u):
    num, den = ankunft_mit_verlust(s50, uu)
    p("     Beispiel Γ fest = Wert der MW bei 50 kpc (s = %.3f) auf dem ganzen Weg, jenseits %.0f·ϖ_R: %.2e" % (s50, uu, num / den))
p("   → In (B3) hängt der Übertritt zwischen Galaxien am seitlichen D und am Verlust auf dem Weg, die das Modell nicht festlegt;")
p("     zwischen 0 und 4r*/(πd). Die Regel des Modells (1D, Γ am Ort) ist (B3) mit seitlich gemitteltem Feld.")
p()

# ======================================================================================
p("=" * 100)
p("5  Blickrichtung der Größen")
p("=" * 100)
for z in (
    "Γ_i = a_N,ges/c          L  am Ort, Betrag des Gesamtfelds, Zustand gebunden/stationär, heute",
    "D_i = GM_i/c, ϖ_R = r*   I  innere Richtung, je Masse",
    "Λ = c²/(2a₀)             I  innere Richtung; G für die Zelle; heute fest (seit z ≈ 10)",
    "Knoten bei Λ             I  (B1, B3) gemeinsam für alle Orte der Zelle",
    "q_ii, ν                  L  am Ort gelesen (Beobachter an der Bahn)",
    "q_ij                     I → L  Rücklauf über den Knoten, gelesen an einem anderen Ort",
    "⟨T_i⟩ (Kac)              I  Zeit in der inneren Richtung, gezählt mit c-Schritten",
    "ΣΔs_b(Zelle), L_Zelle    G  Zelle, heute; d_LSS für unseren Horizont (Beobachter), Annahme L = d_LSS",
    "q_kosm = 0,843           G  gelesen an der Hintergrundstrahlung (Planck [G]), zeitlich fest",
    "Poisson-Ankunft          I → L  seitlicher Ankunftsort (B3)",
):
    p("   " + z)

with open(__file__.replace(".py", ".txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
