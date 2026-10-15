# -*- coding: utf-8 -*-
"""
Dunkle Energie als eine Invariante der Grundwerte (Abstraktion statt Aufbau).
Blick (G, Zelle): Ein Anteil Ω = ρ/ρ_krit hängt am Zustand H; eine Zählung wie Q (Verhältnis zweier Maße, beide ∝ a⁻³)
kann Ω_DE/Ω_m nicht geben, weil ρ_DE/ρ_m ∝ a³ läuft. Zeitfest ist nur ρ_DE gegen die Grundwerte selbst:
   k_DE = ρ_DE·G·c²/a₀²      (Vakuumteil: k_s = π³/30 aus Λ_s = π⁴/(15Λ²), Λ = c²/2a₀)
1 k_DE aus θ* [G] mit den Dichten des Modells (ω_b, ω_c = Q·ω_b, ω_r) und fester Form, je a₀.
2 Umgekehrt: a₀, das jede Lesart von k verlangt (θ* fest), gegen SPARC (eigene Daten).
3 Zählung je Radiant für das Ganze (τ_D·H₀ = 2π; t₀ = τ_D/2π): Zustand, nicht Gesetz.
4 Probe von Hand.
Aus ftm_kern: c, G, T_CMB, OMEGA_B_H2, Q_KOSM, N_EFF_STD, K_S, A0_SPARC, tau_D, lambda_vakuum, schliesslaenge.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from ftm_kern import c, G, T_CMB, OMEGA_B_H2, Q_KOSM, N_EFF_STD, K_S, A0_SPARC, tau_D, lambda_vakuum, schliesslaenge

OUT = []
def P(s=""):
    OUT.append(s); print(s)

H100 = 100e3 / 3.0856775814913673e22          # 100 km/s/Mpc in 1/s
Z_STERN = 1089.9                                # [G] Planck
THETA_MESS = 1.0411e-2                          # 100θ* = 1,0411 ± 0,0003 [G]
wb = OMEGA_B_H2
wm = wb * (1 + Q_KOSM)
wg = 2.469e-5 * (T_CMB / 2.7255)**4
wr = wg * (1 + 0.2271 * N_EFF_STD)

def E2(z, wde):                                  # H²/H100², feste Form
    return wm * (1 + z)**3 + wr * (1 + z)**4 + wde

def theta_stern(wde):
    rs = quad(lambda z: c / math.sqrt(3 * (1 + 0.75 * wb / wg / (1 + z))) / (H100 * math.sqrt(E2(z, wde))), Z_STERN, np.inf, limit=400)[0]
    dm = quad(lambda z: c / (H100 * math.sqrt(E2(z, wde))), 0, Z_STERN, limit=400)[0]
    return rs / dm

def rho_aus_omega(w):
    return w * 3 * H100**2 / (8 * math.pi * G)

wde = brentq(lambda w: theta_stern(w) - THETA_MESS, 0.1, 0.6)
h = math.sqrt(wm + wr + wde)
rho_de = rho_aus_omega(wde)
k_s = math.pi**3 / 30

P("=" * 100); P("1 Zeitfeste Größe der Dunklen Energie aus θ* [G], Dichten des Modells, feste Form"); P("=" * 100)
P(f"ω_DE = {wde:.4f}, h = {h:.4f}, Ω_DE = {wde/h**2:.4f}   (ω_b = {wb}, Q = {Q_KOSM:.3f}, ω_r = {wr:.3e})")
P(f"Vakuumteil k_s = π³/30 = {k_s:.4f};  Probe Kern: ρ_s = Λ_s c²/8πG → k = {lambda_vakuum(schliesslaenge(1.1e-10))*c**2/(8*math.pi*G)*G*c**2/(1.1e-10)**2:.4f}")
for a0 in (1.00e-10, 1.04e-10, 1.06e-10, 1.10e-10, 1.13e-10, 1.19e-10):
    k = rho_de * G * c**2 / a0**2
    P(f"  a₀ = {a0*1e10:.2f}:  k_DE = {k:.3f},  k_DE/k_s = {k/k_s:.3f},  fehlend (k_DE − k_s)/k_s = {k/k_s-1:.3f}")

P(); P("=" * 100); P("2 a₀, das eine Lesart k verlangt (θ* fest; σ aus ±0,0003 in 100θ*)"); P("=" * 100)
dw = brentq(lambda w: theta_stern(w) - (THETA_MESS + 3e-6), 0.1, 0.6) - wde
def a0_fuer(k, w=wde):
    return math.sqrt(rho_aus_omega(w) * G * c**2 / k)
for name, fk in (("nur Vakuumteil", 1.0), ("2·Vakuumteil", 2.0), ("3·Vakuumteil", 3.0)):
    a = a0_fuer(fk * k_s); da = abs(a0_fuer(fk * k_s, wde + dw) - a)
    P(f"  {name:16s}: a₀ = {a*1e10:.4f} ± {da*1e10:.4f}·10⁻¹⁰   (SPARC {A0_SPARC[0]*1e10:.2f}–{A0_SPARC[1]*1e10:.2f})")
lo, hi = [rho_de * G * c**2 / a**2 / k_s for a in (A0_SPARC[1], A0_SPARC[0])]
P(f"  SPARC-Band erlaubt k_DE/k_s = {lo:.2f}–{hi:.2f} (Breite ±{(hi-lo)/(hi+lo)*100:.0f} %): jede Zahl darin passt, keine Lesart ist so entscheidbar.")
P(f"  Entscheidbar erst mit σ(a₀)/a₀ ≈ 1 % (θ* gibt a₀ je Lesart auf {abs(a0_fuer(k_s, wde+dw)-a0_fuer(k_s))/a0_fuer(k_s)*100:.2f} %).")

P(); P("=" * 100); P("3 Zählung je Radiant für das Ganze"); P("=" * 100)
H0 = h * H100
a_theta1 = c * H0 / (2 * math.pi)
GYR = 3.15576e16
t0 = quad(lambda z: 1 / ((1 + z) * H100 * math.sqrt(E2(z, wde))), 0, np.inf, limit=400)[0]
a_t0 = c / (2 * math.pi * t0)
P(f"  τ_D·H₀ = 2π (Θ₀ = 1):  a₀ = cH₀/2π = {a_theta1*1e10:.4f}·10⁻¹⁰")
P(f"  t₀ = τ_D/2π:            a₀ = c/(2πt₀) = {a_t0*1e10:.4f}·10⁻¹⁰   (t₀ = {t0/GYR:.3f} Gyr, H₀t₀ = {H0*t0:.4f})")
P(f"  Beide zählen je Radiant wie Q, gelten aber nur in einer Epoche (a₀ fest, H und t laufen): Zustand des Beobachters,")
P(f"  kein Gesetz; sie unterscheiden sich um H₀t₀ = {H0*t0:.3f}, beide liegen im SPARC-Band. Zeitfest ist nur k_DE.")
for a0 in (1.06e-10, 1.10e-10):
    P(f"  a₀ = {a0*1e10:.2f}: τ_D = {tau_D(a0)/GYR:.2f} Gyr, τ_D/2π = {tau_D(a0)/GYR/2/math.pi:.3f} Gyr, Θ₀ = {2*math.pi*a0/(c*H0):.4f}")

P(); P("=" * 100); P("4 Von Hand"); P("=" * 100)
P(f"  ρ_s = π³a₀²/(30Gc²): a₀ = 1,1 → {math.pi**3/30*(1.1e-10)**2/(G*c**2):.4e} kg/m³;  Λ_s c²/8πG = {lambda_vakuum(schliesslaenge(1.1e-10))*c**2/(8*math.pi*G):.4e}")
P(f"  ω_DE ≈ 0,685·0,674² = {0.685*0.674**2:.4f} (Standardmodell [G]) gegen {wde:.4f}")
P(f"  k_DE(1,1) ≈ ω·3H100²/(8πG)·Gc²/a₀² = {0.3112*3*H100**2/(8*math.pi)*c**2/(1.1e-10)**2:.3f}")

with open(__file__.replace(".py", ".txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
