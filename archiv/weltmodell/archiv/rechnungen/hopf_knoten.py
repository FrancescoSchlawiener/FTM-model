# -*- coding: utf-8 -*-
"""
Faltungsknoten als Hopf-Abbildung? Helizität, Verlinkung, Spinorzeichen, Energieskala.
Rechnung, noch nicht im Modell.

A  Knoten des Modells als Feld: Windung (x+iy), komplexer Versatz i·a_J: R_c² = x²+y²+(z−i a_J)².
   Identität: Z0 = conj(R_c²)/(r²+a_J²), Z1 = 2a_J(x+iy)/(r²+a_J²) liegen exakt auf S³ (|Z0|²+|Z1|² = 1).
B  Hopf-Zahl H = (1/16π²)∫A·B (FFT, Coulomb-Eichung) für w = Z1^p/Z0^q; Abhängigkeit von p (Windung), q, a_J.
C  Unabhängig: Gauß-Verschlingungsintegral zweier Urbildkurven (analytisch parametrisiert).
D  Appell-Knoten 1/R_c (Multipole M(ia)^l, Q(ia)^l): Windung um den Ring, Helizität des eigenen EM-Feldes.
E  Spinorhebung: Drehung um 2π → Vorzeichen; Äquivarianz Raumdrehung α ↔ Zieldrehung pα.
F  Verlinkung von Kurven auf verschachtelten Tori: rationale Wicklung ganz, irrationale Wicklung stetig.
G  Eichung: lokale Dichte A·B nicht eindeutig, Integral fest (Beobachterfehler: H ist eine Größe des Ganzen).
H  Energieskala: E2 ∝ a, E4 ∝ 1/a (Derrick); E* = 2√(E2·E4) gegen |H|.
"""
import math
import numpy as np
from ftm_kern import *

rng = np.random.default_rng(1)


def Zkoord(X, Y, Z, a):
    r2 = X**2 + Y**2 + Z**2
    Rc2 = X**2 + Y**2 + (Z - 1j * a)**2
    return np.conj(Rc2) / (r2 + a * a), 2 * a * (X + 1j * Y) / (r2 + a * a), Rc2


def nfeld(p, q, a, L, N, halb=False):
    x = np.linspace(-L, L, N, endpoint=False) + L / N
    h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    Z0, Z1, _ = Zkoord(X, Y, Z, a)
    num = Z1**p if p >= 0 else np.conj(Z1)**(-p)
    den = np.sqrt(Z0) if halb else Z0**q
    s = np.abs(num)**2 + np.abs(den)**2
    cc = num * np.conj(den)
    return np.array([2 * cc.real / s, 2 * cc.imag / s, (np.abs(den)**2 - np.abs(num)**2) / s]), h


def D4(f, h, ax):
    d = (8 * (np.roll(f, -1, ax) - np.roll(f, 1, ax)) - (np.roll(f, -2, ax) - np.roll(f, 2, ax))) / (12 * h)
    g = np.gradient(f, h, axis=ax, edge_order=2)
    sl = [slice(None)] * f.ndim
    for i in [0, 1, -2, -1]:
        sl[ax] = i
        d[tuple(sl)] = g[tuple(sl)]
    return d


def helizitaet(n, h, chi=None):
    d = [D4(n, h, i + 1) for i in range(3)]
    F = lambda i, j: np.einsum('a...,a...->...', n, np.cross(d[i], d[j], axis=0))
    B = np.array([F(1, 2), F(2, 0), F(0, 1)])
    N = n.shape[1]
    k = 2 * np.pi * np.fft.fftfreq(N, h)
    K = np.array(np.meshgrid(k, k, k, indexing='ij'))
    k2 = (K**2).sum(0)
    k2[0, 0, 0] = 1
    A = np.real(np.fft.ifftn(1j * np.cross(K, np.fft.fftn(B, axes=(1, 2, 3)), axis=0) / k2, axes=(1, 2, 3)))
    H = -np.sum(A * B) * h**3 / (16 * np.pi**2)          # Vorzeichen so, dass w = Z1/Z0 H = +1 hat
    E2 = sum(np.sum(di**2) for di in d) * h**3
    E4 = np.sum(B**2) * h**3
    out = dict(H=H, E2=E2, E4=E4)
    if chi is not None:                                   # Eichtransformation A → A + ∇χ
        Ag = A + np.array([D4(chi, h, i) for i in range(3)])
        dens, densg = (A * B).sum(0), (Ag * B).sum(0)
        out.update(Hg=-np.sum(Ag * B) * h**3 / (16 * np.pi**2),
                   dloc=np.max(np.abs(densg - dens)) / np.max(np.abs(dens)))
    return out


def to_R3(X1, X2, X3, X4, a):
    return a * np.array([X1, X2, X3]) / (1 - X4)


def kurve(theta, phi, xi, a):
    Z1 = np.sin(theta) * np.exp(1j * phi)
    Z0 = np.cos(theta) * np.exp(1j * xi)
    return to_R3(Z1.real, Z1.imag, Z0.imag, Z0.real, a)


def gauss_lk(c1, c2):
    """Gauß-Verschlingungsintegral zweier geschlossener Polygonzüge (3×n)."""
    d1 = np.roll(c1, -1, 1) - c1
    d2 = np.roll(c2, -1, 1) - c2
    m1 = c1 + d1 / 2
    m2 = c2 + d2 / 2
    R = m1[:, :, None] - m2[:, None, :]
    cr = np.cross(d1[:, :, None], d2[:, None, :], axis=0)
    return np.sum((R * cr).sum(0) / np.linalg.norm(R, axis=0)**3) / (4 * np.pi)


print("=" * 100)
print("A  Knoten: Windung (x+iy), Verwindung als komplexer Versatz i·a_J  →  Koordinaten auf S³")
print("=" * 100)
P = rng.normal(size=(3, 200000)) * 3
for a in [0.3, 1.0, 4.0]:
    Z0, Z1, Rc2 = Zkoord(*P, a)
    print(f"A1 a_J = {a:4.1f}: max| |Z0|²+|Z1|² − 1 | = {np.max(np.abs(np.abs(Z0)**2 + np.abs(Z1)**2 - 1)):.1e}"
          f"   (|r²−a²+2iaz|² + 4a²ρ² = (r²+a²)²)")
print("   Nullstellen: Z1 = 0 auf der Drehachse (Windung), Z0 = 0 ⇔ R_c² = 0 auf dem Ring ρ = a_J, z = 0 (Verwindung).")
print("   ⇒ Windung und komplexer Versatz sind genau die inverse stereographische Abbildung ℝ³∪∞ → S³ mit Radius a_J;")
print("     w = Z1^p/Z0^q ist die (p,q)-Hopf-Abbildung S³ → S², n = (2Re w, 2Im w, 1−|w|²)/(1+|w|²).")
print("   a_J = 0: Z1 ≡ 0, Z0 ≡ 1 ⇒ konstante Abbildung, kein Knoten.")

print()
print("=" * 100)
print("B  Hopf-Zahl H = (1/16π²)∫A·B,  B_i = ½ε_ijk n·(∂_j n × ∂_k n),  ∇×A = B  (Gitter, FFT)")
print("=" * 100)
tabE = {}
for (p, q) in [(1, 1), (-1, 1), (2, 1), (1, 2), (3, 1), (2, 2), (1, 3)]:
    zeile = []
    for N in [128, 176]:
        n, h = nfeld(p, q, 1.0, 6.0, N)
        r = helizitaet(n, h)
        zeile.append(r['H'])
    tabE[(p, q)] = r
    print(f"B1 p = {p:+d} (Windung), q = {q} : H(N=128) = {zeile[0]:+.4f}, H(N=176) = {zeile[1]:+.4f}   p·q = {p*q:+d}")
x = np.linspace(-6, 6, 128, endpoint=False) + 6 / 128
h = x[1] - x[0]
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
Z0, Z1, _ = Zkoord(X, Y, Z, 1.0)
w = Z1 * Z0
r0 = helizitaet(np.array([2 * w.real, 2 * w.imag, 1 - np.abs(w)**2]) / (1 + np.abs(w)**2), h)
print(f"B2 w = Z0·Z1 (|w| ≤ ½, trifft nicht die ganze Kugel): H = {r0['H']:+.4f}  (Windung allein, ohne Pol am Ring: 0)")
print("   Gegenladung: Windung −1 heißt conj(Z1): H = −q (Tabelle p = −1).")
for a in [0.25, 0.5, 1.0, 2.0]:
    n, h = nfeld(1, 1, a, 6.0, 160)
    print(f"B3 feste Box ±6, N = 160, a_J = {a:4.2f}: H = {helizitaet(n, h)['H']:+.4f}   (Abweichung nur Auflösung/Box: h/a = {h/a:.2f}, L/a = {6/a:.0f})")
print("   ⇒ H = p·q, ganz; hängt linear an der Windung p (Ladung) und am Exponenten q der Verwindungsnullstelle,")
print("     nicht am Wert von a_J (Skala; H ist skalenfrei). a_J = 0 ⇒ H = 0.")

print()
print("=" * 100)
print("C  Unabhängig: Verschlingung zweier Urbildkurven w = w1, w = w2 (Gauß-Integral)")
print("=" * 100)
M = 1500
t = np.linspace(0, 2 * np.pi, M, endpoint=False)
for (p, q) in [(1, 1), (2, 1), (1, 2), (3, 2), (2, 2)]:
    dd = math.gcd(p, q)
    comps = []
    for th, beta in [(0.6, 0.3), (1.1, 2.0)]:
        cs = []
        for k in range(dd):
            phi = (q // dd) * t + (beta + 2 * np.pi * k) / p
            xi = (p // dd) * t
            cs.append(kurve(th, phi, xi, 1.0))
        comps.append(cs)
    lk = sum(gauss_lk(c1, c2) for c1 in comps[0] for c2 in comps[1])
    print(f"C1 (p,q) = ({p},{q}): Σ Verschlingung der Komponenten = {lk:+.4f}   p·q = {p*q}  (Vorzeichen: Umlaufsinn der Kurven)   (Komponenten je Urbild: {dd})")
print("   Probe Faser (1,1): zwei Hopf-Kreise, Verschlingung 1 (von Hand: zwei Fasern der Hopf-Faserung sind einmal verschlungen).")

print()
print("=" * 100)
print("D  Appell-Knoten 1/R_c (Multipole M(ia)^l, Q(ia)^l): Windung um den Ring und Helizität seines Feldes")
print("=" * 100)
a = 1.0
for eps in [0.3, 0.1]:
    u = np.linspace(0, 4 * np.pi, 4001)
    rho = a + eps * np.cos(u)
    z = eps * np.sin(u)
    Rc = np.sqrt(rho**2 + (z - 1j * a)**2 + 0j)               # Hauptzweig: Schnitt über der Scheibe ρ < a, z = 0
    Rc_stetig = np.empty_like(Rc)
    Rc_stetig[0] = Rc[0]
    for i in range(1, len(u)):                                 # stetig fortgesetzt
        Rc_stetig[i] = Rc[i] if abs(Rc[i] - Rc_stetig[i - 1]) < abs(Rc[i] + Rc_stetig[i - 1]) else -Rc[i]
    ph = np.unwrap(np.angle(1 / Rc_stetig))
    i1 = 2000
    print(f"D1 Schleife um den Ring, Abstand {eps}: Phase von 1/R_c nach 1 Umlauf {(ph[i1]-ph[0])/(2*np.pi):+.4f}·2π, "
          f"nach 2 Umläufen {(ph[-1]-ph[0])/(2*np.pi):+.4f}·2π;  1/R_c nach 1 Umlauf = {(1/Rc_stetig[i1])/(1/Rc_stetig[0]):+.4f}·Start")
print("   ⇒ 1/R_c hat um den Ring die Windung ½ (q = ½): eindeutig erst auf zwei Blättern, Rückkehr nach zwei Umläufen.")
for N in [128, 176]:
    n, h = nfeld(1, 1, 1.0, 6.0, N, halb=True)
    print(f"D2 w = Z1/√Z0 auf einem Blatt (Sprung an der Scheibe), N = {N}: H = {helizitaet(n, h)['H']:+.4f}  (nicht ganz, gitterabhängig)")
# eigenes EM-Feld des Appell-Knotens: E + iB = −∇(Q/R_c)
Pp = rng.normal(size=(3, 50000)) * 2
Pp = Pp[:, np.abs(Pp[2]) > 0.05]
hh = 1e-5
def Fc(P):
    return 1 / np.sqrt(P[0]**2 + P[1]**2 + (P[2] - 1j * a)**2 + 0j)
G = np.array([(Fc(Pp + hh * e[:, None]) - Fc(Pp - hh * e[:, None])) / (2 * hh) for e in np.eye(3)])
E, B = -G.real, -G.imag
ph_hat = np.array([-Pp[1], Pp[0], 0 * Pp[0]]) / np.hypot(Pp[0], Pp[1])
Bphi = np.abs((B * ph_hat).sum(0)).max() / np.abs(B).max()
EB = (E * B).sum(0)
Pm = Pp * np.array([[1], [1], [-1]])
Gm = np.array([(Fc(Pm + hh * e[:, None]) - Fc(Pm - hh * e[:, None])) / (2 * hh) for e in np.eye(3)])
EBm = (-Gm.real * -Gm.imag).sum(0)
print(f"D3 Kerr-Newman-Feld E + iB = −∇(Q/R_c): max|B_φ|/max|B| = {Bphi:.1e} (B poloidal) ⇒ A = A_φ φ̂, A·B = 0 punktweise: ∫A·B = 0;")
print(f"   E·B ungerade in z: max|E·B(z) + E·B(−z)|/max|E·B| = {np.abs(EB + EBm).max()/np.abs(EB).max():.1e} ⇒ ∫E·B = 0.")
print("   ⇒ Das eigene Feld des Knotens trägt keine Helizität; H ≠ 0 nur in der Abbildung (Z0,Z1) → S², nicht im Feld 1/R_c.")

print()
print("=" * 100)
print("E  Spinorhebung und Drehung")
print("=" * 100)
sig = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
psi = rng.normal(size=2) + 1j * rng.normal(size=2)
psi /= np.linalg.norm(psi)
nv = lambda s: np.real([np.conj(s) @ S @ s for S in sig])
for al in [2 * np.pi, 4 * np.pi]:
    U = np.diag([np.exp(-1j * al / 2), np.exp(1j * al / 2)])
    print(f"E1 Drehung von n um α = {al/np.pi:.0f}π: n zurück: {np.allclose(nv(U @ psi), nv(psi))};  Spinor Uψ = {np.real_if_close((U @ psi)[0] / psi[0]):+.4f}·ψ")
print("   ⇒ Das Zeichen −1 nach 2π gehört zu jeder Abbildung n = ψ†σψ (Faser der Hopf-Faserung), unabhängig von H.")
x0 = rng.normal(size=(3, 1000))
for (p, q) in [(1, 1), (2, 1), (1, 2), (3, 1)]:
    al = 0.7
    Rz = np.array([[np.cos(al), -np.sin(al), 0], [np.sin(al), np.cos(al), 0], [0, 0, 1]])
    Z0a, Z1a, _ = Zkoord(*x0, 1.0)
    Z0b, Z1b, _ = Zkoord(*(Rz @ x0), 1.0)
    wa, wb = Z1a**p / Z0a**q, Z1b**p / Z0b**q
    dev = np.abs(wb / wa - np.exp(1j * p * al)).max()
    print(f"E2 (p,q) = ({p},{q}), H = {p*q}: w(R_α x) = e^(ipα)·w(x), Abw. {dev:.1e};  SU(2)-Hebung der Zieldrehung 2πp: (−1)^p = {(-1)**p:+d}")
print("   ⇒ In diesem Ansatz folgt das Zeichen nach einer Raumdrehung um 2π der Windung p, nicht H: (1,2) hat H = 2 und −1, (2,1) H = 2 und +1.")
print("   Ob der Knoten als Spin ½ zu quantisieren ist (Klasse der Schleife im Raum der Abbildungen, π₄(S²) = ℤ₂), rechnet diese Datei nicht.")

print()
print("=" * 100)
print("F  Verlinkung von Kurven auf verschachtelten Tori (Tori θ = const um den Ring ρ = a, z = 0)")
print("=" * 100)
M = 1200
t = np.linspace(0, 2 * np.pi, M, endpoint=False)
for (m1, k1), (m2, k2) in [((1, 0), (0, 1)), ((0, 1), (1, 0)), ((1, 1), (1, 1)), ((2, 1), (1, 3)), ((1, 3), (2, 1)), ((3, 2), (2, 5))]:
    cA = kurve(1.2, m1 * t, k1 * t + 0.1, 1.0)   # näher am Ring (θ → π/2)
    cB = kurve(0.5, m2 * t, k2 * t + 0.7, 1.0)   # näher an der Achse
    lk = gauss_lk(cA, cB)
    print(f"F1 innen (φ-Umläufe {m1}, Ring-Umläufe {k1}), außen ({m2}, {k2}): Verschlingung = {lk:+.4f}   m_innen·k_außen = {m1*k2}")
print("   ⇒ ganz genau bei geschlossenen (rationalen) Kurven: Lk = m_innen·k_außen.")
for W in [0.5, (math.sqrt(5) - 1) / 2, 1 / math.sqrt(2), 0.75]:
    res = []
    for T in [20, 40]:
        tt = np.linspace(0, 2 * np.pi * T, 60 * T, endpoint=False)
        cA = kurve(1.2, tt, W * tt + 0.1, 1.0)
        cB = kurve(0.5, tt, W * tt + 0.7, 1.0)
        res.append(gauss_lk(cA, cB) / T**2)                       # Schluss durch die Sehne (Polygon)
    print(f"F2 Wicklung W = {W:.5f} auf beiden Tori, offene Bahn T Umläufe, Lk/T²: T=20 {res[0]:.4f}, T=40 {res[1]:.4f}   (W = {W:.4f})")
print("   ⇒ Bei irrationaler Wicklung schließt keine Kurve; die mittlere Verschlingung je Umlauf² ist W, eine stetige Zahl, keine ganze.")
print("   Die Tori von 01 liegen im Zustandsraum; die Tori dieses Knotens im Ortsraum um den Ring ρ = a_J.")

print()
print("=" * 100)
print("G  Eichung: A → A + ∇χ")
print("=" * 100)
n, h = nfeld(1, 1, 1.0, 6.0, 128)
x = np.linspace(-6, 6, 128, endpoint=False) + 6 / 128
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
chi = 3 * np.exp(-((X - 0.5)**2 + Y**2 + (Z + 0.3)**2) / 0.8) * X
r = helizitaet(n, h, chi=chi)
print(f"G1 H = {r['H']:+.5f}, nach Eichung {r['Hg']:+.5f};  lokale Dichte A·B ändert sich um bis zu {r['dloc']*100:.0f} % ihres Maximums.")
print("   ⇒ A·B ist kein Wert am Ort; nur das Integral über den ganzen Knoten ist eine Größe (G).")

print()
print("=" * 100)
print("H  Energieskala: E2 = ∫|∂n|² ∝ a, E4 = ∫|B|² ∝ 1/a  (Box ±6a, N = 128; Tabelle B: N = 176)")
print("=" * 100)
for a in [0.7, 1.4]:
    n, h = nfeld(1, 1, a, 6.0 * a, 128)
    r = helizitaet(n, h)
    print(f"H1 a = {a}: E2/a = {r['E2']/a:.2f}, E4·a = {r['E4']*a:.2f}   (unendlicher Raum: beide 32π² = {32*np.pi**2:.2f}; E2 fehlt der Rest außerhalb der Box, ≲ 128π/6 = {128*np.pi/6:.0f})")
print("   ⇒ nur E2 (keine Viererterm): E → 0 für a → 0, kein Gleichgewicht (Derrick). Größe setzt im Modell a_J = J/(Mc), keine Energie.")
Hs, Es = [], []
for (p, q), r in tabE.items():
    if p < 0:
        continue
    Es.append(2 * math.sqrt(r['E2'] * r['E4']))
    Hs.append(p * q)
    print(f"H2 (p,q) = ({p},{q}), |H| = {p*q}: E* = 2√(E2·E4) = {Es[-1]:8.1f},  E*/|H|^(3/4) = {Es[-1]/(p*q)**0.75:6.1f}")
k = np.polyfit(np.log(Hs), np.log(Es), 1)[0]
print(f"   Anpassung E* ∝ |H|^k im Ansatz (nicht minimiert): k = {k:.2f}. Eine Schranke E ≥ c|H|^(3/4) braucht den Viererterm E4 in der Energie;")
print("   das Modell hat keine Energie für die Abbildung n, also folgt daraus keine Stabilitätsgrenze.")

print()
print("=" * 100)
print("Von Hand: Punkt (1,1,1), a = 1: |conj R_c²|² = |r²−a²+2iaz|² = |2+2i|² = 8;  |2a(x+iy)|² = 4·2 = 8;  (r²+a²)² = 4² = 16  ⇒ 8+8 = 16 ✓")
print("          E4(1,1) = ∫16·Ω⁴dV, Ω = 2a/(r²+a²): 256·4π·∫r²/(1+r²)⁴dr = 1024π·π/32 = 32π² = 315,83 (Gitter H1: 315,36)")
