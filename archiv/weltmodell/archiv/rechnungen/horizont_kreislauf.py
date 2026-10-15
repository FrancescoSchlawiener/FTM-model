# -*- coding: utf-8 -*-
"""Horizontentropie: H-Form πkc²/(ℓ_P²H²) gegen Λ-Form (45/π³)k(Λ/ℓ_P)², Mindestrate der neuen Faltung,
Bedingung −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k). Rechnet nur mit ftm_kern."""
import math, sys
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ftm_kern import *

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

SI = 7.1e86                       # Ṡ_i heute, k/s je Horizont (S_obs/t₀)
H0 = H0_si()
Lam = schliesslaenge(A0)
Os0 = omega_vakuum(A0)
OL = 1 - OMEGA_M - OMEGA_R
H_s = math.sqrt(lambda_vakuum(Lam) * c**2 / 3)          # Atemrate, die der Vakuumteil allein trägt
S_L = 45 / math.pi**3 * (Lam / l_P)**2                    # Λ-Form
M_je_m = 3 * c**2 / G                                      # Δs = GM/(3c²): Masse je Meter Faltung

def E2_std(a):  return OMEGA_M * a**-3 + OMEGA_R * a**-4 + OL          # Standardhintergrund des Kerns
def E2_vak(a):  return OMEGA_M * a**-3 + OMEGA_R * a**-4 + Os0         # nur der Vakuumteil trägt die Dunkle Energie
def alter(a, E2):
    return quad(lambda x: 1 / (x * math.sqrt(E2(x))), 1e-12, a, limit=400)[0] / H0 / yr

P("=== 1 Beide Formen ===")
P(f"Λ = c²/2a₀ = {Lam/Gpc:.2f} Gpc, Ω_s(0) = {Os0:.4f}, H_s/H₀ = {H_s/H0:.4f} (= √Ω_s: {math.sqrt(Os0):.4f})")
P(f"S_Λ = {S_L:.3e} k, horizontentropie(H_s) = {horizontentropie(H_s):.3e} k")
S_H0 = horizontentropie(H0)
P(f"z = 0:        S_H = {S_H0:.3e} k, S_H/S_Λ = {S_H0/S_L:.4f}")
for a in (2, 5, 10, 100):
    P(f"a = {a:>4} (Std): S_H/S_Λ = {Os0/E2_std(a):.4f};   (nur Vakuumteil): {Os0/E2_vak(a):.4f}")
P(f"a → ∞  (Std): S_H/S_Λ = Ω_s/Ω_Λ = {Os0/OL:.4f}  → nie auf 10 % gleich")
P(f"a → ∞  (nur Vakuumteil): 1")
a10 = brentq(lambda a: Os0 / E2_vak(a) - 0.9, 1, 1e3)
P(f"10 %-Gleichheit im Grenzfall: a = {a10:.3f} (z = {1/a10-1:.3f}), Weltalter {alter(a10, E2_vak)/1e9:.1f} Mrd J"
  f" (dieser Hintergrund, heute {alter(1, E2_vak)/1e9:.1f} Mrd J; Std heute {alter(1, E2_std)/1e9:.2f})")
P(f"  Probe analytisch: a³ = 9Ω_m/Ω_s (ohne Strahlung) → a = {(9*OMEGA_M/Os0)**(1/3):.3f}")

P("\n=== 2 Mindestrate ===")
# S_H = πc²/(ℓ_P²H²),  Ṡ_H/S_H = −2Ḣ/H.  H² = H₀²(Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_rest) + H_s²(Λ),  H_s² ∝ Λ⁻²
# ⇒ Ṡ_H = S_H·H·(3Ω_m(a) + 4Ω_r(a)) + 2·Ω_s(a)·S_H·Λ̇/Λ,   Ω_s(a)·S_H = S_Λ  (für jeden Hintergrund)
# ⇒ Λ̇/Λ ≥ [Ṡ_i − S_H·H·(3Ω_m(a)+4Ω_r(a))] / (2 S_Λ)
r_L = SI / (2 * S_L); r_H = SI / (2 * S_H0)
P(f"Λ-Form Ṡ_i/(2S_Λ) = {r_L:.3e} s⁻¹ = {r_L*yr:.2e}/J → {r_L*yr*Lam*M_je_m/M_sun:.2e} M☉/J;  kreislauf_mindestrate: {kreislauf_mindestrate(Lam)[0]:.3e}")
P(f"H-Form Ṡ_i/(2S_H) heute = {r_H:.3e} s⁻¹ = {r_H*yr:.2e}/J → {r_H*yr*Lam*M_je_m/M_sun:.2e} M☉/J")
mat0 = S_H0 * H0 * (3 * OMEGA_M + 4 * OMEGA_R)
P(f"Heute wächst der Horizont durch die Verdünnung der Materie: Ṡ_H = {mat0:.2e} k/s = {mat0/SI:.1e}·Ṡ_i → keine Mindestrate für Λ")

P("\n=== 3 Bedingung −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k) ===")
def rest(a, E2, Si=SI):
    H = H0 * math.sqrt(E2(a)); mHd = H0**2 * (1.5 * OMEGA_M * a**-3 + 2 * OMEGA_R * a**-4)
    return mHd - l_P**2 * H**3 * Si / (2 * math.pi * c**2)
P(f"heute: −Ḣ = {H0**2*(1.5*OMEGA_M+2*OMEGA_R):.2e} s⁻², verlangt {l_P**2*H0**3*SI/(2*math.pi*c**2):.2e} s⁻²")
for name, E2 in (("Std (Ω_Λ fest)", E2_std), ("nur Vakuumteil", E2_vak)):
    ax = brentq(lambda la: rest(10**la, E2), 0, 12); ax = 10**ax
    P(f"{name}: Ṡ_i = 7,1·10⁸⁶ k/s je Horizont fest → verletzt ab a = {ax:.3e}, Weltalter {alter(ax, E2):.3e} J")
    if name.startswith("Std"):
        tA = 2 / (3 * H0 * math.sqrt(OL)) * math.asinh(math.sqrt(OL / OMEGA_M) * ax**1.5) / yr
        P(f"   Probe analytisch (ohne Strahlung): t = 2/(3H₀√Ω_Λ)·arsinh(√(Ω_Λ/Ω_m)a^1,5) = {tA:.3e} J")
# Ṡ_i je Horizont bei fester Produktion je mitbewegtem Volumen: ∝ (a·r_H)⁻³... hier: Horizontvolumen physikalisch fest → Ṡ_i ∝ Dichte der Quellen ∝ a⁻³
ax = 10**brentq(lambda la: rest(10**la, E2_std, SI * 10**(-3*la)), -1, 12) if rest(1e12, E2_std, SI*1e-36) < 0 else None
P(f"Std, Ṡ_i ∝ a⁻³ (Quellen verdünnen sich im festen Horizontvolumen): verletzt ab a = {ax}")
P(f"Im Grenzfall H → H_s ist Ḣ = 0 nur ohne neue Faltung; mit Λ̇/Λ ≥ Ṡ_i/(2S_Λ) = {r_L:.2e} s⁻¹ gilt −Ḣ = H_s·Λ̇/Λ ≥ ℓ_P²H_s³Ṡ_i/(2πc²k), "
  f"w_s + 1 = (2/3)(Λ̇/Λ)/H_s = {2/3*r_L/H_s:.1e}")

P("\n=== 4 Probe von Hand ===")
lp = 1.616255e-35; Lh = 8.98755e16 / 2.2e-10; H0h = 67.4 / 3.085678e19
SLh = 45 / 31.00628 * (Lh / lp)**2; SHh = 3.141593 * (2.99792e8 / H0h / lp)**2
P(f"Λ = {Lh:.4e} m, S_Λ = {SLh:.3e}, S_H = {SHh:.3e}, Verhältnis {SHh/SLh:.4f}, Raten {7.1e86/2/SLh:.3e} / {7.1e86/2/SHh:.3e} s⁻¹")
P(f"Masse je Jahr: 3c²·Λ·rate·yr/G = {3*8.98755e16*Lh*7.1e86/2/SLh*3.15576e7/6.6743e-11/1.98892e30:.2e} M☉")

open(__file__.replace('.py', '.txt'), 'w').write("\n".join(out) + "\n")
