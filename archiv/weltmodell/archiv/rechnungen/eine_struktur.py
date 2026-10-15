# -*- coding: utf-8 -*-
"""
Eine Struktur hinter allen Bausteinen: gezählte Rückkehr mit komplexem Nome q = e^{2πiτ}
(|q| = Verlust, arg q = Windung), abgelesen am Ort, Schluss ohne Außen.
Rechnung, noch nicht im Modell.

A  Perioden der Rückkehr je Ebene: Temperatur = ħ/(k·Periode), Schließlänge Λ = π·c²/κ
B  Verschränkung als gemeinsamer Knoten: CHSH für (L) mitgeführte Eigenschaft und
   (I) Knoten, an dem beide Messrichtungen ansetzen; Rolle der Grundschwingung der Windung.
   Unabhängige Proben: Aufzählung der 16 deterministischen Strategien, Eigenwert des
   CHSH-Operators zweier Spins ½, Monte-Carlo der mitgeführten Richtung.
"""
import itertools
import math
import numpy as np
from ftm_kern import *

print("=" * 96)
print("A  Perioden der Rückkehr (imaginäre Zeit) je Ebene; kT = ħ/Periode; Λ = π·c²/κ")
print("=" * 96)
H0 = H0_si()
Lam = schliesslaenge(A0)
tD = tau_D(A0)
print(f"innere Richtung: τ_D = c/a₀ = 2Λ/c = {tD:.4e} s = {tD/Gyr:.1f} Mrd a,  Λ = {Lam/(1e9*pc):.2f} Gpc")
print(f"   kT_I = ħ/τ_D = ħa₀/c = {hbar/tD/kB:.4e} K")
kap_H = c * H0
per_H = 2 * math.pi * c / kap_H
print(f"Horizont (global): κ = cH₀, Periode 2πc/κ = 2π/H₀ = {per_H:.4e} s;  T_GH = {hbar/per_H/kB:.4e} K")
print(f"   Λ_H = π·c²/κ = πc/H₀ = {math.pi*c/H0/(1e9*pc):.2f} Gpc;  Λ/Λ_H = 2πa₀/(cH₀) = Θ = {Lam/(math.pi*c/H0):.4f}")
M = M_sun
r0s = r0(M)
kap_l = c**2 / (2 * r0s)
per_l = 2 * math.pi * c / kap_l
print(f"Falte (lokal, Sonnenmasse): κ = c²/(2r₀), Periode = 24πΔs/c = {per_l:.4e} s  (Probe 24πΔs/c = {24*math.pi*ds(M)/c:.4e})")
print(f"   T = ħc/(24πkΔs) = {hbar/per_l/kB:.4e} K;  Λ_loc = π·c²/κ = 2π·r₀ = {math.pi*c**2/kap_l:.1f} m (Probe 2πr₀ = {2*math.pi*r0s:.1f})")
print("   ⇒ auf jeder Ebene: Periode = 2Λ/c (Kac-Zeit der bei Λ geschlossenen Halbgeraden), Λ = π·c²/κ, a = κ/2π")
print()

print("=" * 96)
print("B  Verschränkung: gemeinsamer Knoten, CHSH  S = E(a,b) − E(a,b') + E(a',b) + E(a',b')")
print("=" * 96)

# B1 unabhängige Probe: 16 deterministische Strategien (mitgeführte Eigenschaft, L-Lesung)
best = 0
for A, A2, B, B2 in itertools.product([-1, 1], repeat=4):
    S = A * B - A * B2 + A2 * B + A2 * B2
    best = max(best, abs(S))
print(f"B1 mitgeführte Eigenschaft (lokal, jede Mischung deterministischer Strategien): max |S| = {best}")

# B2 L-Lesung konkret: Richtung λ im Knoten, von der Quelle mitgeführt, A = sign cos(λ−a), B = −sign cos(λ−b)
rng = np.random.default_rng(1)
lam = rng.uniform(0, 2 * np.pi, 2_000_000)


def E_mc(a, b):
    return np.mean(np.sign(np.cos(lam - a)) * -np.sign(np.cos(lam - b)))


def E_dreieck(d):
    d = abs((d + np.pi) % (2 * np.pi) - np.pi)
    return -(1 - 2 * d / np.pi)


def S_winkel(E, th):
    # a = 0, b = th, a' = 2th, b' = 3th  →  S = 3E(th) − E(3th) für drehinvariantes E
    return abs(E(th) - E(3 * th) + E(th) + E(th))


ths = np.linspace(0, np.pi / 2, 20001)
S_tri = max(S_winkel(E_dreieck, t) for t in ths[::50])
a, b, a2, b2 = 0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4
S_mc = abs(E_mc(a, b) - E_mc(a, b2) + E_mc(a2, b) + E_mc(a2, b2))
print(f"B2 mitgeführte Richtung (Monte-Carlo, Winkel 0/45/90/135°): S = {S_mc:.4f};  analytisch (Dreieck E = −(1 − 2Δ/π)) max S = {S_tri:.4f}")

# B3 I-Lesung: beide Messrichtungen setzen am Knoten an, Korrelation = Grundschwingung E = −cos(kΔ)
S_cos = max(S_winkel(lambda d: -np.cos(d), t) for t in ths)
print(f"B3 Knoten, Grundschwingung E = −cos Δ: max S = {S_cos:.6f}  (2√2 = {2*math.sqrt(2):.6f}), Winkel π/4")
th2 = ths[np.argmax([S_winkel(lambda d: -np.cos(2 * d), t) for t in ths])]
print(f"   Photon (Polarisation kehrt nach π zurück, E = −cos 2Δ): max S = 2√2 bei θ = {np.degrees(th2):.2f}° ≙ 22,5° mod 45°")

# B4 unabhängige Probe: größter Eigenwert des CHSH-Operators zweier Spins ½
sx = np.array([[0, 1], [1, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)


def n_sig(t):
    return np.cos(t) * sz + np.sin(t) * sx


Aop, A2op, Bop, B2op = n_sig(0), n_sig(np.pi / 2), n_sig(np.pi / 4), n_sig(-np.pi / 4)
C = np.kron(Aop, Bop) + np.kron(Aop, B2op) + np.kron(A2op, Bop) - np.kron(A2op, B2op)
print(f"B4 Eigenwert CHSH-Operator (zwei Spins ½): {max(abs(np.linalg.eigvalsh(C))):.6f}")

# B5 kein Signal: Randverteilung am Knoten mit E = −cos Δ, P(A=±1) unabhängig von b
psi = np.array([0, 1, -1, 0], complex) / math.sqrt(2)   # Singulett
for bb in (0.0, 0.7, 2.0):
    PA = []
    for sgn in (+1, -1):
        Pa = (np.eye(2) + sgn * n_sig(0.3)) / 2
        PA.append(np.real(psi.conj() @ np.kron(Pa, np.eye(2)) @ psi))
    Eab = np.real(psi.conj() @ np.kron(n_sig(0.3), n_sig(bb)) @ psi)
    print(f"B5 b = {bb:.1f}: P(A=+1) = {PA[0]:.4f}, P(A=−1) = {PA[1]:.4f}, E(a,b) = {Eab:+.4f}, −cos(a−b) = {-math.cos(0.3-bb):+.4f}")

# B6 Rolle der Grundschwingung: Dreieck = (8/π²)Σ cos((2m+1)Δ)/(2m+1)²; die ersten N Glieder, auf E(0) = −1 normiert
print("B6 Oberschwingungen der Windung (Dreieck abgeschnitten nach N Gliedern, E(0) = −1):")
for N in (1, 2, 3, 5, 10, 100, 1000):
    w = np.array([1 / (2 * m + 1) ** 2 for m in range(N)])
    w /= w.sum()
    ks = np.array([2 * m + 1 for m in range(N)])
    E = lambda d, w=w, ks=ks: -np.sum(w * np.cos(ks * d))
    S = max(S_winkel(E, t) for t in ths[::20])
    print(f"   N = {N:5d}: Anteil Grundschwingung {w[0]:.4f}, max S = {S:.4f}")
print("   Mischung Grundschwingung (1−ε) + Dreieck ε:")
for eps in (0.0, 0.1, 0.25, 0.5, 1.0):
    E = lambda d, eps=eps: (1 - eps) * -np.cos(d) + eps * E_dreieck(d)
    S = max(S_winkel(E, t) for t in ths[::20])
    print(f"   ε = {eps:.2f}: max S = {S:.4f}")
print("   Grenze ohne Einschränkung (nur kein Signal, Popescu-Rohrlich): S = 4")
