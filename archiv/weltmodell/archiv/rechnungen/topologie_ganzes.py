# -*- coding: utf-8 -*-
# Das Ganze als topologisches System: welche Größen hängen nur an der Form des geschlossenen Mediums, welche an der Dynamik?
# A Homologie von T³ × innere Richtung (Kettenkomplex, Ränge numerisch)
# B Spektrum des Laplace-Operators: Raumzelle (Würfel, sechseckig, bcc, fcc bei gleichem Volumen) und innere Richtung
# C Gauß auf dem geschlossenen Raum: Nullmode, Gegenquelle, Rückstellung; Ladung
# D Vakuumteil als Summe über die Rückläufe n
# E Legt die Topologie Q fest?
# F Windung: Ladung, Spin, Umlaufsinn
# G Triade als Hodge-Zerlegung, Kindebene
import sys, io, math
import numpy as np
from itertools import product
from scipy.special import zeta
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ftm_kern import *

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []; P = lambda s="": out.append(s)
Gpc = 1e3 * 1e6 * pc
pi = math.pi

# ============================ A Homologie ============================
def kette_kreis(n):      # 0-Zellen n, 1-Zellen n, periodisch
    d = np.zeros((n, n))
    for i in range(n): d[i, i] -= 1; d[(i + 1) % n, i] += 1
    return [n, n], [d]
def kette_strecke(m):    # m Punkte, m-1 Kanten
    d = np.zeros((m, m - 1))
    for i in range(m - 1): d[i, i] -= 1; d[i + 1, i] += 1
    return [m, m - 1], [d]
def tensor(A, B):
    (na, da), (nb, db) = A, B
    ka, kb = len(na) - 1, len(nb) - 1
    K = ka + kb
    blocks = {k: [(i, k - i) for i in range(ka + 1) if 0 <= k - i <= kb] for k in range(K + 1)}
    dims = [sum(na[i] * nb[j] for i, j in blocks[k]) for k in range(K + 1)]
    ds = []
    for k in range(1, K + 1):
        M = np.zeros((dims[k - 1], dims[k]))
        off_lo = {}; o = 0
        for (i, j) in blocks[k - 1]: off_lo[(i, j)] = o; o += na[i] * nb[j]
        o = 0
        for (i, j) in blocks[k]:
            sz = na[i] * nb[j]
            if i >= 1:
                r = off_lo[(i - 1, j)]; blk = np.kron(da[i - 1], np.eye(nb[j]))
                M[r:r + blk.shape[0], o:o + sz] += blk
            if j >= 1:
                r = off_lo[(i, j - 1)]; blk = (-1) ** i * np.kron(np.eye(na[i]), db[j - 1])
                M[r:r + blk.shape[0], o:o + sz] += blk
            o += sz
        ds.append(M)
    return dims, ds
def betti(C):
    dims, ds = C
    rk = [np.linalg.matrix_rank(d) for d in ds]
    return [int(dims[k] - (rk[k - 1] if k >= 1 else 0) - (rk[k] if k < len(ds) else 0)) for k in range(len(dims))]

P("A) Homologie (Kettenkomplex, Ränge numerisch; Betti b₀, b₁, …)")
S = kette_kreis(3); I = kette_strecke(3)
T3 = tensor(tensor(S, S), S)
P(f"   T³                         {betti(T3)}   Euler {sum((-1)**k*b for k,b in enumerate(betti(T3)))}")
P(f"   T³ × Halbgerade [0, Λ]     {betti(tensor(T3, I))}   (innere Richtung zusammenziehbar: keine neue Klasse)")
P(f"   T³ × Kreis 2Λ mit Umlauf   {betti(tensor(T3, kette_kreis(3)))}   (4. Windung, 4. Erhaltungsgröße)")
P("   Kreis 2Λ ohne Umlauf = Kreis / Spiegelung = Strecke [0, Λ]: Betti (1, 0), wie die Halbgerade.")
P("   Metrische Nullmoden des flachen T³: 6 = dim GL(3)/O(3): 1 Volumen (Atmen) + 5 Scherung; Killing-Felder 3 (Verschiebungen, b₁ = 3).")
P()

# ============================ B Spektrum ============================
def gitter(name):
    s3 = math.sqrt(3)
    A = {"Würfel": np.eye(3),
         "sechseckig (h = Seite)": np.array([[1, 0, 0], [0.5, s3 / 2, 0], [0, 0, 1.0]]),
         "raumzentriert (bcc)": 0.5 * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1.0]]),
         "flächenzentriert (fcc)": 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0.0]])}[name]
    return A / abs(np.linalg.det(A)) ** (1 / 3)          # Volumen 1
def kurz(A, N=3):
    v = [np.array(n) @ A for n in product(range(-N, N + 1), repeat=3) if any(n)]
    return np.array(sorted(np.linalg.norm(x) for x in v))
V = zelle_aus_a0(A0) ** 3                                # Λ legt das Volumen fest
Lam = schliesslaenge(A0)
P(f"B) Spektrum. Λ = ν·ΣΔs_b = Gρ_mV/(3c²) sieht nur das Volumen: a₀ = {A0:.2e} ⇒ V^(1/3) = {V**(1/3)/Gpc:.2f} Gpc, Λ = {Lam/Gpc:.2f} Gpc")
P("   Form bei gleichem Volumen (L_min kürzeste Schleife, λ₁ kleinster Eigenwert ≠ 0, g₁ Entartung, a₀,max aus L_min ≥ 0,985·d_LSS, Kac 1/Packungsdichte)")
dL = d_LSS()
a0max_w = None
for nm in ["Würfel", "sechseckig (h = Seite)", "raumzentriert (bcc)", "flächenzentriert (fcc)"]:
    A = gitter(nm); l = kurz(A); B = 2 * pi * np.linalg.inv(A).T; k = kurz(B)
    lam1 = k[0] ** 2; g1 = int(np.sum(np.abs(k ** 2 - lam1) < 1e-9))
    dens = (4 / 3) * pi * (l[0] / 2) ** 3                   # Volumen 1
    Vmin = (KREISPAAR_GRENZE * dL / l[0]) ** 3
    a0max = 3 * c**4 / (2 * G * rho_m_mitbewegt() * Vmin)
    P(f"   {nm:24s} L_min/V^⅓ = {l[0]:.4f}  λ₁V^⅔/(2π)² = {lam1/(2*pi)**2:.4f}  g₁ = {g1:2d}  a₀,max = {a0max*1e10:.3f}·10⁻¹⁰  Kac = {1/dens:.3f}")
Lw = V ** (1 / 3)
P(f"   Würfel: Entartungen der Schalen |n|² = 1…8 (Thetareihe Σq^|n|²): {[sum(1 for n in product(range(-3,4),repeat=3) if sum(x*x for x in n)==m) for m in range(1,9)]}")
P("   Innere Richtung [0, Λ], Ort als Rand (G(0) = 1: homogen Dirichlet), bei Λ reflektierend: k_m = (2m+1)π/(2Λ), keine Nullmode.")
k_raum = 2 * pi / Lw; k_inn = pi / (2 * Lam)
P(f"   k₁ Raum = 2π/L = {k_raum*Gpc:.4f}/Gpc,  k₁ innen = π/2Λ = {k_inn*Gpc:.4f}/Gpc  ⇒ tiefste Mode des Ganzen: räumlich konstant (k = 0), nur in der inneren Richtung.")
P(f"   Probe 2Λ/L = {2*Lam/Lw:.4f} gegen Ω_m(LH₀/c)²/4π = {OMEGA_M*(Lw*H0_si()/c)**2/(4*pi):.4f}")
P(f"   Verhältnis 4Λ/L = {4*Lam/Lw:.3f} heute; L atmet (∝ a), Λ fest ⇒ epochenabhängig, nicht topologisch (2Λ/L = Ω_m(LH/c)²/4π, Würfel).")
P()

# ============================ C Gauß ============================
P("C) Gauß auf dem geschlossenen T³: ∫∇²Φ dV = 0 (Grundklasse, b₃ = 1) ⇒ Σ Quelle = 0.")
n = 32; rng = np.random.default_rng(1)
rho = np.zeros((n, n, n)); pos = rng.integers(0, n, (5, 3)); ms = rng.uniform(1, 3, 5)
for p_, m_ in zip(pos, ms): rho[tuple(p_)] += m_
R = np.fft.fftn(rho); kk = 2 * pi * np.fft.fftfreq(n); KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij")
k2 = KX**2 + KY**2 + KZ**2
P(f"   5 Punktmassen, FFT-Poisson: Quellanteil bei k = 0: {R[0,0,0].real:.3f} = Σm = {ms.sum():.3f}; −k²Φ̂ = R̂ ist bei k = 0 nur lösbar mit R̂(0) = 0,")
P(f"   also mit der gleichförmigen Gegenquelle −ρ̄ = −Σm/V = {-ms.sum()/n**3:.3e} (je Zelle). Sie ist die Nullmode: ohne Ort, ohne Außen.")
Phi = np.real(np.fft.ifftn(np.where(k2 > 0, -(R - np.where(k2 == 0, R, 0)) / np.where(k2 > 0, k2, 1), 0)))
lap = (sum(np.roll(Phi, s, a) for a in range(3) for s in (1, -1)) - 6 * Phi)
P(f"   Probe: diskretes ∇²Φ summiert über die Zelle = {lap.sum():.1e}")
rho_m = rho_m_mitbewegt(); H = H0_si()
P(f"   Die Gegenquelle wirkt als Feld a = +(4πG/3)ρ̄·r (gleichförmige Dehnung); (4πG/3)ρ̄_m/H₀² = {(4*pi*G/3)*rho_m/H**2:.4f} = Ω_m/2 = {OMEGA_M/2:.4f}.")
P("   Das ist genau das Bremsglied des Atmens, ä/a = −(4πG/3)ρ̄: Die Gegenquelle ist die Volumen-Nullmode (Atmen), skaliert ∝ a⁻³ (w = 0).")
P(f"   In Faltungsgrößen an der Kugel R: Φ_gegen/c² = (2πG/3c²)ρ̄R² = (3/2)·ΣΔs/R exakt (Probe: {(2*pi*G/3)/((4*pi*G/9)):.4f}).")
P(f"   Vakuumteil zum Vergleich: Ω_s = {omega_vakuum(A0):.3f}, w = −1 (Λ fest). Die Gegenquelle trägt σ = ΣΔs − ΔL = 0, nicht Ω_s.")
P("   Ladung: dieselbe Bedingung gilt für ∇·E = ρ_el/ε₀, und dort nimmt kein Atmen eine Gegenquelle auf ⇒ Σ Q_el = 0 im geschlossenen Raum.")
P("   Unterschied: Windungen ±n heben sich auf, Faltungen Δs ≥ 0 nicht ⇒ nur die Faltung erzwingt eine Rückstellung.")
P()

# ============================ D Vakuumteil ============================
P("D) Vakuumteil Ω_s = (π²/45)Θ²: π⁴/15 = Γ(4)·Σ_{n≥1} 1/n⁴, n = Zahl der Rückläufe (Windungen der Rückkehr), Gewicht 1/n⁴ aus d = 3.")
tot = 6 * zeta(4)
for N in (1, 2, 3, 10):
    P(f"   bis n = {N:2d}: Γ(4)Σ = {6*sum(1/k**4 for k in range(1,N+1)):.5f}  Anteil {6*sum(1/k**4 for k in range(1,N+1))/tot:.4f}")
P(f"   π⁴/15 = {pi**4/15:.5f}. Topologisch: 3 (Dimension), ganze n (Schließung). Dynamisch: Θ = 2πa₀/(cH) = {theta_(A0, H):.3f} heute.")
P()

# ============================ E Q ============================
Qm, dQ = Q_KOSM, 0.065
P(f"E) Legt die Topologie Q fest? gemessen Q = {Qm:.3f} ± {dQ} (Planck 2018 [G]).")
kand = [("Holonomieordnung 1,2,3,4,6 (orientierbare flache 3-Mannigfaltigkeiten; Hantzsche-Wendt 4)", [1, 2, 3, 4, 6]),
        ("Kusszahl (Bildpunkte nächster Nachbarn) Würfel 6, bcc 8, fcc 12, sechseckig 8", [6, 8, 12]),
        ("Σ Betti T³ = 8", [8])]
for nm, vs in kand:
    P(f"   ν = {vs}: Q = ν − 1, bester Abstand {min(abs(v-1-Qm) for v in vs)/dQ:.1f}σ  – {nm}")
P(f"   Kac in die einbeschriebene Kugel: 1,35–1,91 (Tabelle B) ⇒ Q = 0,35–0,91.")
P(f"   Zählung je Radiant ν = 2π: Q = {2*pi-1:.3f}, {abs(2*pi-1-Qm)/dQ:.1f}σ; 2π ist Länge des Einheitskreises (Metrik der Zählung), keine Homologiezahl.")
P("   Topologische Zahlen sind ganz oder rational; Q = 5,364 ± 0,065 ist 5,6σ von 5 und 9,8σ von 6 entfernt.")
P("   Was die Gleichung erzwingt: ν·ρ_b·V = 3c⁴/(2Ga₀). Mit L_min ≥ 0,985·d_LSS folgt nur eine obere Schranke für ν bei festem a₀:")
for nm in ["Würfel", "sechseckig (h = Seite)", "raumzentriert (bcc)", "flächenzentriert (fcc)"]:
    l = kurz(gitter(nm))[0]; Vmin = (KREISPAAR_GRENZE * dL / l) ** 3
    rho_b = rho_m / (1 + Qm)
    numax = 3 * c**4 / (2 * G * A0 * rho_b * Vmin)
    P(f"   {nm:24s} ν ≤ {numax:.2f}  (Q ≤ {numax-1:.2f})")
P("   Die Zellform ändert nur die Schranke; ein Wert für Q folgt nicht.")
P()

# ============================ F Windung ============================
P("F) Windung.")
th = np.linspace(0, 2 * pi, 2001)[:-1]
for nw in (1, -1, 2):
    z = np.exp(1j * nw * th) * (1 + 0.3 * np.cos(3 * th)) + 0.1 * (rng.standard_normal(th.size) + 1j * rng.standard_normal(th.size))
    w = np.sum(np.angle(np.roll(z, -1) / z)) / (2 * pi)
    P(f"   Phase mit Windung {nw:+d} und Störung (Amplitude ≠ 0): gemessen {w:+.4f}")
z = np.exp(1j * th) - 1.2
P(f"   Amplitude durch null geschoben (Kreis um den Knoten verfehlt ihn): Windung {np.sum(np.angle(np.roll(z,-1)/z))/(2*pi):+.4f} ⇒ n ändert sich nur über |ψ| = 0.")
sz = np.array([[1, 0], [0, -1]], complex)
U = np.cos(pi) * np.eye(2) - 1j * np.sin(pi) * sz    # exp(−i·2π·σ_z/2)
P(f"   Spin: Drehung um 2π in SU(2) = {U[0,0].real:+.0f}·𝟙 ⇒ π₁(SO(3)) = ℤ₂: topologisch ist nur J/ħ ∈ ½ℤ, nicht j.")
je = 0.5 * (m_Pl / m_e) ** 2; jp = 0.5 * (m_Pl / m_p) ** 2
P(f"   j = ½(m_Pl/M)²: Elektron {je:.2e}, Proton {jp:.2e} – stetig in M, dynamisch.")
P(f"   Ladung: π₁(U(1)) = ℤ bzw. c₁ ∈ H²(T²; ℤ) = ℤ; Dirac 1/(2α) = {1/(2*alpha):.2f}.")
P("   Innere Richtung: H¹ = 0 (Tabelle A) ⇒ dort keine Windungszahl. Ihr einer Sinn kommt vom Rand: Die Halbgerade hat genau einen Randpunkt (den Ort),")
P("   die Richtung weg vom Rand ist ausgezeichnet (1 Bit, Zeitsinn), der Kern bleibt gerade im Drehsinn.")
P()

# ============================ G Triade ============================
P("G) Triade als Hodge-Zerlegung auf T² (FFT): Feld = harmonisch (Nullmode, b₁ = 2 Zahlen) + Gradient (Quellen) + Rotation (Wirbel).")
n = 256; x = np.arange(n) / n; X, Y = np.meshgrid(x, x, indexing="ij")
def gauss(x0, y0, w=0.02):
    dx = (X - x0 + 0.5) % 1 - 0.5; dy = (Y - y0 + 0.5) % 1 - 0.5
    return np.exp(-(dx**2 + dy**2) / (2 * w * w))
src = 300 * (gauss(0.3, 0.5) - gauss(0.7, 0.5))                    # Quellpaar, Σ = 0
vort = 300 * gauss(0.5, 0.2)                                    # Wirbel
k = 2 * pi * np.fft.fftfreq(n, 1 / n); KX, KY = np.meshgrid(k, k, indexing="ij"); K2 = np.where(KX**2 + KY**2 > 0, KX**2 + KY**2, 1)
ph = -np.fft.fftn(src) / K2; ph[0, 0] = 0; ps = -np.fft.fftn(vort - vort.mean()) / K2; ps[0, 0] = 0
U0 = np.array([0.3, -0.1])
Fx = U0[0] + np.real(np.fft.ifftn(1j * KX * ph)) + np.real(np.fft.ifftn(1j * KY * ps))
Fy = U0[1] + np.real(np.fft.ifftn(1j * KY * ph)) - np.real(np.fft.ifftn(1j * KX * ps))
Fxh, Fyh = np.fft.fftn(Fx), np.fft.fftn(Fy)
harm = np.array([Fxh[0, 0].real, Fyh[0, 0].real]) / n**2
div = (KX * Fxh + KY * Fyh) / K2; div[0, 0] = 0
gx, gy = np.real(np.fft.ifftn(KX * div)), np.real(np.fft.ifftn(KY * div))
rx, ry = Fx - harm[0] - gx, Fy - harm[1] - gy
E = lambda a, b: np.mean(a * a + b * b)
P(f"   zurückgewonnen: harmonisch {harm.round(6)} (gesetzt {U0}); Anteile an ⟨F²⟩: harmonisch {E(harm[0]+0*Fx,harm[1]+0*Fy)/E(Fx,Fy):.3f}, Gradient {E(gx,gy)/E(Fx,Fy):.3f}, Rotation {E(rx,ry)/E(Fx,Fy):.3f}")
P(f"   Quellpaar exakt im Gradienten: max|Gradient − ∇φ| = {np.max(np.abs(gx-np.real(np.fft.ifftn(1j*KX*ph)))):.1e}; Wirbel exakt in der Rotation.")
# Kindebene: Quadrat fern der Quellen
m = (X > 0.45) & (X < 0.55) & (Y > 0.7) & (Y < 0.8)
gxm, gym = gx[m], gy[m]          # nur das Quellfeld der Eltern
ux, uy = gxm.mean(), gym.mean()
P(f"   Kindebene (Quadrat 0,1 fern der Quellen): Dort ist das Feld der Eltern quellfrei und wirbelfrei, also ganz Hintergrund der Kindebene;")
P(f"   gleichförmiger Teil {(ux**2+uy**2)/np.mean(gxm**2+gym**2):.3f}, Rest Gezeit {1-(ux**2+uy**2)/np.mean(gxm**2+gym**2):.3f}. Mit Rand ist der harmonische Teil nicht mehr b₁-dimensional,")
P("   sondern durch die Randwerte gesetzt: der Hintergrund einer Ebene ist das Feld der Elternebene.")

open(__file__.replace(".py", ".txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out))
