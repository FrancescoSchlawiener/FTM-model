# -*- coding: utf-8 -*-
"""
Elektromagnetismus als Phasenseite derselben Rückkehr (Nome q = e^{2πiτ}).
Rechnung, noch nicht im Modell.

A  Geometrie des Nomes: Re τ ist ein Kreis (Windung, beide Sinne), Im τ eine Halbgerade (Verlust, ein Sinn).
   Summe der Rückläufe: ein Sinn mit Verlust → ν (Schwelle); ohne Verlust nur beide Sinne reell → Kamm (ganze Ladung, keine Skala).
B  dF = 0 als Aussage "die Phase ist eine Windung": Gitter-Zusammenhang, Fluss durch jede geschlossene Fläche = 2π·ganz.
   Quellengleichung auf dem geschlossenen T³: lösbar nur mit ΣQ = 0.
C  g-Faktor: Faltung und Windung am selben komplexen Ort (Versatz i·a_J). Multipole M(ia)^l, Q(ia)^l → μ = Q·c·a_J, g = 2.
   Unabhängige Proben: exaktes Kerr-Newman-Potential (Einstein-Maxwell), klassisch starr drehende Ladung (g = 1).
D  Keine Schwelle auf der Phasenseite: Laborfelder mit qE/m < a₀.
E  Beobachterfehler (Tabelle im Text der Meldung).
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from ftm_kern import *

eps0 = 8.8541878128e-12
kC = 1 / (4 * math.pi * eps0)

# ================================================================= A
print("=" * 96)
print("A  Nome q = e^{2πiτ}: Re τ Kreis (Phase), Im τ Halbgerade (Betrag)")
print("=" * 96)
t_im = 0.3
x = np.linspace(0, 1, 2001)
qq = np.exp(2j * np.pi * (x + 1j * t_im))
wind = np.sum(np.diff(np.unwrap(np.angle(qq)))) / (2 * np.pi)
print(f"A1 Re τ: 0 → 1 bei Im τ = {t_im}: Windung von arg q = {wind:+.4f}, |q| konstant: max−min = {np.ptp(np.abs(qq)):.1e}")
print(f"   Re τ: 0 → −1: Windung {-wind:+.4f}  (beide Sinne; Verlust |q| unberührt)")
print("   Im τ ∈ (0, ∞): Halbgerade, ein Rand (|q| = 1), kein geschlossener Weg ⇒ keine Windung, ein Sinn.")
print("   Streifen / (τ ~ τ+1) = Kreis × Halbgerade: H¹ = ℤ nur vom Kreis. Phase = Windung, Betrag = Verlust.")

print("A2 Summe der Rückläufe  Σ q^n:")
for s in [3.0, 1.0, 0.3, 0.1]:
    print(f"   ein Sinn, Verlust s = {s:4.1f}, Phase 0: Σ_n≥0 e^(−ns) = {1/(-math.expm1(-s)):8.4f} = ν(a_N = s²a₀) = {float(nu(s*s*A0)):8.4f}")
# ohne Verlust (s → 0): ein Sinn gibt 1/(1 − e^{iφ}) = ½ + (i/2)cot(φ/2): nicht reell
phi = np.array([0.5, 1.0, 2.0])
for eps in [1e-2, 1e-4]:
    one = 1 / (1 - np.exp(-eps + 1j * phi))
    print(f"   ein Sinn, s = {eps:.0e}, φ = {phi}: Σ = " + ", ".join(f"{z.real:.3f}{z.imag:+.3f}i" for z in one)
          + "  (Imaginärteil ½cot(φ/2): kein reeller Zählwert)")
# beide Sinne: Poisson-Kern sinh s/(cosh s − cos φ) → 2π Σ δ(φ − 2πm)
f = lambda p: np.exp(-((p - 0.0) ** 2) / (2 * 0.4**2)) + 0.5 * np.exp(-((p - 2 * np.pi) ** 2) / (2 * 0.4**2))
pp = np.linspace(-np.pi, 3 * np.pi, 400001)
dp = pp[1] - pp[0]
for s in [1.0, 0.1, 0.01, 0.001]:
    P = np.sinh(s) / (np.cosh(s) - np.cos(pp))
    val = np.sum(P * f(pp)) * dp / (2 * np.pi)
    print(f"   beide Sinne, s = {s:6.3f}: (1/2π)∫P_s(φ)f(φ)dφ = {val:.5f}  → Σ_m f(2πm) = {f(0.0)+f(2*np.pi):.5f} (Kamm: nur ganze Windungen)")
print("   ⇒ Mit Verlust und einem Sinn: ν(s), Skala a₀ über s = √(a_N/a₀). Ohne Verlust ist nur die Zählung in beiden Sinnen reell;")
print("     sie gibt den Kamm Σδ(φ − 2πm): ganze Windungen, keine Kurve, keine Skala.")
print()

# ================================================================= B
print("=" * 96)
print("B  dF = 0 und Quellengleichung")
print("=" * 96)
rng = np.random.default_rng(1)
N = 8


def plaq(Ax, Ay, Az):
    """Plaketten-Winkel F_ij = A_i(x) + A_j(x+i) − A_i(x+j) − A_j(x) auf periodischem N³-Gitter."""
    A = [Ax, Ay, Az]
    F = {}
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        F[(i, j)] = A[i] + np.roll(A[j], -1, axis=i) - np.roll(A[i], -1, axis=j) - A[j]
    return F


def cube_sum(F):
    # orientierte Summe über die 6 Flächen eines Würfels (Rand des Würfels)
    return ((np.roll(F[(0, 1)], -1, axis=2) - F[(0, 1)])
            + (np.roll(F[(1, 2)], -1, axis=0) - F[(1, 2)])
            + (np.roll(F[(2, 0)], -1, axis=1) - F[(2, 0)]))


for amp, name in [(np.pi, "beliebige Phase je Glied"), (0.3, "glatte Phase")]:
    A = [rng.uniform(-amp, amp, (N, N, N)) for _ in range(3)]
    F = plaq(*A)
    exact = cube_sum(F)
    Fbar = {k: v - 2 * np.pi * np.round(v / (2 * np.pi)) for k, v in F.items()}  # Winkel mod 2π
    m = cube_sum(Fbar) / (2 * np.pi)
    print(f"B1 {name}: Σ_Rand F (ohne mod) max = {np.abs(exact).max():.1e} (Bianchi dF = 0 exakt);")
    print(f"   Winkel mod 2π: Fluss je Würfel/2π ganzzahlig: max|m − round m| = {np.abs(m - np.round(m)).max():.1e}, "
          f"Würfel mit m ≠ 0: {int(np.sum(np.round(m) != 0))}, Σm über den geschlossenen Raum = {int(np.round(m).sum())}")
print("   ⇒ Ist A der Zusammenhang einer Phase (nur mod 2π bestimmt), ist der Fluss jeder geschlossenen Fläche 2π·ganz;")
print("     bei glatter Phase null: dF = 0 (∇·B = 0, ∇×E + ∂_tB = 0). m ≠ 0 nur, wo die Phase springt (Amplitude null).")

# Quellengleichung auf T³ (FFT): ∇·E = ρ lösbar nur mit Σρ = 0
N = 32
L = 1.0
k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2
for charges, name in [([(+1, (8, 8, 8)), (-1, (20, 16, 12))], "Windungen +1, −1"),
                      ([(+1, (8, 8, 8)), (+1, (20, 16, 12))], "Windungen +1, +1")]:
    rho = np.zeros((N, N, N))
    for qv, (i, j, l) in charges:
        rho[i, j, l] += qv / (L / N) ** 3
    rk = np.fft.fftn(rho)
    Ktmp = K2.copy(); Ktmp[0, 0, 0] = 1
    phik = rk / Ktmp; phik[0, 0, 0] = 0
    lap = np.real(np.fft.ifftn(-K2 * phik))
    res = np.abs(-lap - rho).max() / np.abs(rho).max()
    print(f"B2 {name}: Σρ·V = {rho.sum()*(L/N)**3:+.1f}, Rest |−∇²φ − ρ|/max ρ = {res:.2e}"
          + ("  (lösbar)" if res < 1e-8 else "  (nicht lösbar: k = 0 bleibt, keine Nullmode nimmt sie auf)"))
print("   ⇒ Aus d⋆F = J auf dem geschlossenen Raum: ΣQ_el = 0; dJ = 0 (Ladungserhaltung) aus d² = 0.")
print("     d⋆F = J selbst braucht eine Wirkung für F (Metrik, Hodge-Stern): die niedrigste lokale, eichinvariante ist ∫F∧⋆F.")
print()

# ================================================================= C
print("=" * 96)
print("C  g-Faktor: Faltung und Windung am selben komplexen Ort z = i·a_J")
print("=" * 96)
# Komplexes Potential eines bei z = i a versetzten Zentrums (Einheiten Q = 1, a = 1): Φ = 1/√(x²+y²+(z−ia)²)
a = 1.0
lmax = 6
xg, wg = leggauss(64)
for r in [20.0, 100.0]:
    th = np.arccos(xg)
    X = r * np.sin(th); Z = r * np.cos(th)
    Phi = 1 / np.sqrt(X**2 + (Z - 1j * a) ** 2 + 0j)
    # Legendre-Projektion: Φ = Σ c_l P_l(cosθ)/r^{l+1}
    from numpy.polynomial import legendre as Lg
    cl = []
    for l in range(lmax):
        Pl = Lg.legval(xg, [0] * l + [1])
        cl.append((2 * l + 1) / 2 * np.sum(wg * Phi * Pl) * r ** (l + 1))
    print(f"C1 r = {r:5.0f}: Multipole c_l = " + ", ".join(f"{z.real:+.4f}{z.imag:+.4f}i" for z in cl)
          + f"   Probe (ia)^l: {', '.join(str((1j*a)**l) for l in range(4))} …")
print("   Realteil: Monopol Q, Quadrupol −Q a², …; Imaginärteil: Dipol Q·a ⇒ magnetisches Moment μ = Q·a (Gauß), SI μ = Q·c·a.")
print("   Dieselbe Verschiebung am Betrag (Faltung): Imaginärteil des Dipols = M·a = J/c, also a = J/(Mc) = a_J.")
print("   g ≡ 2Mμ/(Q J) = 2M·(Q c a_J)/(Q·M c a_J) = 2, für jedes M, Q, a_J.")

# Probe 1: exaktes Kerr-Newman (Boyer-Lindquist, geometrisierte Gauß-Einheiten): A_φ = −Q a r sin²θ/Σ
M_, Q_, a_ = 1.0, 0.5, 0.7
for r in [10.0, 100.0, 1000.0]:
    th = 1.1
    Sig = r**2 + a_**2 * np.cos(th) ** 2
    Aphi = -Q_ * a_ * r * np.sin(th) ** 2 / Sig
    mu_ = -Aphi * r / np.sin(th) ** 2          # Dipol: A_φ = −μ sin²θ/r (Vorzeichen der Drehung)
    gtp = -2 * M_ * a_ * r * np.sin(th) ** 2 / Sig
    J_ = -gtp * r / (2 * np.sin(th) ** 2)      # g_tφ = −2J sin²θ/r
    print(f"C2 Kerr-Newman r = {r:6.0f}: μ = {mu_:.5f}, J = {J_:.5f}, g = 2Mμ/(QJ) = {2*M_*mu_/(Q_*J_):.5f}")
# Probe 2: klassisch, starr drehend, Ladung folgt der Masse
R, w = 1.0, 1.0
mu_shell = Q_ * w * R**2 / 3; J_shell = 2 / 3 * M_ * R**2 * w
mu_ring = Q_ * w * R**2 / 2; J_ring = M_ * R**2 * w
print(f"C3 klassisch starr: Kugelschale g = {2*M_*mu_shell/(Q_*J_shell):.3f}, Ring g = {2*M_*mu_ring/(Q_*J_ring):.3f}")
print("   (Ladung und Masse drehen im Raum; beim Knoten sitzen beide am selben imaginären Versatz ⇒ g = 2.)")
# Probe 3: Faltungsgrößen aus 08: GJ/c³ = 3Δs·a_J
aJ_e = hbar / (2 * m_e * c)
print(f"C4 Elektron: a_J = ħ/(2m_ec) = {aJ_e:.5e} m; μ = e·c·a_J = {e_el*c*aJ_e:.5e} J/T; μ_B = eħ/2m_e = {e_el*hbar/(2*m_e):.5e}")
print(f"   Probe GJ/c³ = 3Δs·a_J: {G*(hbar/2)/c**3:.5e} = {3*ds(m_e)*aJ_e:.5e} m²;  j = a_J/(3Δs) = {aJ_e/(3*ds(m_e)):.3e}")
rQ_e = math.sqrt(G * kC * e_el**2 / c**4)
print(f"   Ladungslänge r_Q = √(G e²/(4πε₀c⁴)) = {rQ_e:.3e} m; r_Q/3Δs = {rQ_e/(3*ds(m_e)):.3e}, a_J/r_Q = {aJ_e/rQ_e:.3e}: überextremal,")
print("   ohne Horizont; der Multipolsatz M(ia)^l, Q(ia)^l hängt nicht an M, Q, a ⇒ g = 2 gilt auch dort.")
print("   Nicht aus dieser Rechnung: die Abweichung g − 2 (gemessen ≈ 2·0,00116, Literatur) und zusammengesetzte Knoten (Proton).")
print()

# ================================================================= D
print("=" * 96)
print("D  Keine Schwelle auf der Phasenseite: qE/m < a₀")
print("=" * 96)
for name, qm in [("Elektron", e_el / m_e), ("Proton", e_el / m_p), ("Körper 1 g, 1 pC", 1e-12 / 1e-3), ("Körper 1 kg, 1 pC", 1e-12)]:
    print(f"D1 {name:18s}: q/m = {qm:.3e} C/kg, qE/m = a₀ bei E* = a₀m/q = {A0/qm:.3e} V/m")
r_ee = math.sqrt(kC * e_el**2 / (m_e * A0))
print(f"   zwei Elektronen: Coulomb-Beschleunigung = a₀ bei r = √(ke²/(m_e a₀)) = {r_ee:.3e} m")
print("D2 Das Modell: Das Argument der Regel ist das Taktgefälle der Faltung, Γ = c|∇ln(dτ/dt)| = a_N/c, nicht die Beschleunigung.")
print("   Eine Phase ändert |q| und damit den Takt nicht ⇒ qE/m geht nicht in s ein. Die Regel wirkt nur auf g_N:")
for gN_rel in [9.81 / A0, 1.0, 0.01]:
    gN = gN_rel * A0
    aE = 0.1 * A0
    tot_modell = aE + float(nu(gN)) * gN
    print(f"   g_N = {gN_rel:9.3g} a₀, qE/m = 0,1 a₀: EM-Anteil der Bahn = qE/m exakt (Faktor 1,000000), Regel ν(g_N) − 1 = {float(nu(gN))-1:.3e}")
print("   Zum Vergleich, wäre ν auf qE/m angewandt (nicht das Modell):")
for aE_rel in [1.0, 0.1, 0.01]:
    print(f"   qE/m = {aE_rel:5.2f} a₀: Faktor ν = {float(nu(aE_rel*A0)):.3f}")
# Betragsseite eines Laborfeldes: Takt aus der Feldenergie (Reissner-Nordström-Term, gravitativ)
for Qlab, rlab in [(1e-6, 1.0), (e_el, 1e-10)]:
    rQ2 = G * kC * Qlab**2 / c**4
    dlt = rQ2 / rlab**2
    grad = c**2 * rQ2 / rlab**3
    print(f"D3 Feldenergie einer Ladung {Qlab:.2e} C bei r = {rlab:g} m: Taktanteil r_Q²/r² = {dlt:.2e}, Gefälle c²r_Q²/r³ = {grad:.2e} m/s²"
          f" (gravitativ, addiert sich zu g_N)")
print("   ⇒ Keine elektromagnetische Abweichung bei beliebig kleinen Feldern (Coulomb 1/r², Lorentzkraft unverändert);")
print("     die Schwelle kommt allein über den Verlust der Faltung. Vorhersage: kein elektromagnetischer MOND-Effekt,")
print("     auch nicht bei qE/m ≪ a₀ und nicht im tiefen Bereich (g_N ≪ a₀): dort ändert sich nur der gravitative Teil.")
print()

# ================================================================= Hand
print("=" * 96)
print("Von Hand: μ = e·c·a_J = 1,602177e-19 · 2,997925e8 · 1,930796e-13 = 4,803205e-11 · 1,930796e-13 = 9,2740e-24 J/T")
print("          E* = a₀·m_e/e = 1,1e-10 · 9,10938e-31 / 1,602177e-19 = 1,00203e-40 / 1,602177e-19 = 6,254e-22 V/m")
