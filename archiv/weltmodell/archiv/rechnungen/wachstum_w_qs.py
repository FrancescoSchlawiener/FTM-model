# -*- coding: utf-8 -*-
"""
Lineares Wachstum mit dem Gewicht w vor Q_s: ν_eff = 1 + (1 − w)(ν − 1). Rechnung, noch nicht im Modell.

Beobachterfrage: Welches Feld sieht die Regel zwischen Strukturen? Nur die Störung selbst (kein Gefälle im Hintergrund):
   a_N,b(R) = Newton-Feld der Baryonen am Rand einer Kugel R mit dem linearen Kontrast δ_R(a):
   a_N,b = ½·Ω_b0·H0²·R_mit·δ_R(a)/a²  (R mitbewegt, δ_R = σ(R)·D(a)/D_Staub(1), selbstkonsistent mit der eigenen Lösung).
Zwei Argumente der Regel:
   (i)  s² = a_N,b/a₀(z)                 (nur das Feld der Störung)
   (ii) s² = a_N,b/a₀(z) + θ̃·τ_D(z)       (Summe der Verluste auf derselben Rückkehr, 03: θ̃ = 3H, τ_D = c/a₀(z))
Zwei Lesarten der Kopplung (Quelle relativ zu Staub mit 1 + Q):
   A  Ersatz:   a = a_N,b·[1 + w·Q + (1 − w)(ν − 1)]  ⇒ μ = 1 + (1 − w)(ν − 1 − Q)/(1 + Q)   (Anteil 1 − w Operator statt Staub)
   B  Zusatz:   a = a_N,b·[1 + Q + (1 − w)(ν − 1)]    ⇒ μ = 1 + (1 − w)(ν − 1)/(1 + Q)       (Regel auf die Bewegung zwischen Strukturen)
Wachstum: D'' + (2 + dlnE/dlna)·D' = 1,5·Ω_m(a)·μ·D, Start z = 100 mit D = a (w dort 1 − e^{−60}).
Spektrum: CAMB (A_s = 2,1·10⁻⁹, n_s = 0,9649 wie stoerungen_rekombination) für Staub mit Ω_c = Q·Ω_b, σ(R) Kugelfilter.
Kern: A0, c, Mpc, H0_si, OMEGA_M, OMEGA_B, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM, E_hintergrund, a0_frueh_relativ, nu, gewicht_w.
"""
import math, os, sys, io
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, c, Mpc, H0_si, H0_KMS, OMEGA_M, OMEGA_B, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM,
                      E_hintergrund, a0_frueh_relativ, nu, gewicht_w)

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)

H0 = H0_si()
RS = [1.0, 3.0, 10.0, 30.0, 100.0]          # Mpc (mitbewegt)

# ------------------------------------------------------------------ Spektrum σ(R)
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=0.0, omk=0.0)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=200.0)
r = camb.get_results(p)
kh, _, pk = r.get_matter_power_spectrum(minkh=1e-4, maxkh=150.0, npoints=4000)
k = kh * H0_KMS / 100; Pk = pk[0] / (H0_KMS / 100)**3             # Mpc⁻¹, Mpc³
def sigR(R):
    x = k * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(np.trapezoid(k**2 * Pk * W**2, k) / (2 * math.pi**2))
sig8 = sigR(8 / (H0_KMS / 100))

# ------------------------------------------------------------------ Hintergrund, Gewicht
def Ea(a): return E_hintergrund(1 / a - 1)
def dlnE(a, h=1e-5):
    return (math.log(Ea(a * math.exp(h))) - math.log(Ea(a * math.exp(-h)))) / (2 * h)
def a0z(a): return A0 * a0_frueh_relativ(1 / a - 1)
def omw(a):
    th = 3 * H0 * Ea(a)
    return 1 - gewicht_w(th, a0z(a)), th * c / a0z(a)
def aNb(a, R, dR):                          # Newton-Feld der Baryonen der Störung, m/s²
    return 0.5 * OMEGA_B * H0**2 * (R * Mpc) * abs(dR) / a**2

def mu(a, R, dR, arg, les):
    om, x = omw(a)
    if om == 0.0: return 1.0
    s2 = aNb(a, R, dR) / a0z(a) + (x if arg == "ii" else 0.0)
    nm1 = 1 / math.expm1(math.sqrt(s2))     # ν − 1 = n_B(s)
    if les == "A": return 1 + om * (nm1 - Q_KOSM) / (1 + Q_KOSM)
    return 1 + om * nm1 / (1 + Q_KOSM)

ai = 1 / 101
def wachstum(R=None, sR=None, arg=None, les=None, a_end=(1.0, 1 / 1.5, 0.5)):
    """y = (D, dD/dlna); D normiert auf D(ai) = ai; δ_R(a) = sR·D(a)/Dst1."""
    def rhs(lna, y):
        a = math.exp(lna); D, Dp = y
        Om_a = OMEGA_M * a**-3 / Ea(a)**2
        m = 1.0 if R is None else mu(a, R, sR * D / Dst1, arg, les)
        return [Dp, -(2 + dlnE(a)) * Dp + 1.5 * Om_a * m * D]
    sol = solve_ivp(rhs, (math.log(ai), 0.0), [ai, ai], rtol=1e-10, atol=1e-14, dense_output=True, method="DOP853")
    return {ae: sol.sol(math.log(ae)) for ae in a_end}
st = wachstum()
Dst1 = st[1.0][0]

hdr("0  Eingaben")
P(f"   a₀ = {A0:.2e}, Q = {Q_KOSM:.4f}, Ω_m = {OMEGA_M:.4f}, Ω_b = {OMEGA_B:.5f}, H₀ = {H0_KMS}; σ₈ (CAMB, A_s = 2,1·10⁻⁹) = {sig8:.4f}")
P(f"   Staub: D(1)/D(z=100)·(1/101) = {Dst1:.5f}, f(0) = {st[1.0][1] / st[1.0][0]:.5f}  (Probe Ω_m(0)^0,55 = {OMEGA_M**0.55:.5f})")
P("   1 − w im Hintergrund (θ̃ = 3H, a₀(z) aus dem Kern):  " + ",  ".join(f"z = {z}: {omw(1 / (1 + z))[0]:.3e}" for z in (0, 0.5, 1, 2, 3, 5)))

hdr("1  Feld der Störung und ν_eff − 1 = (1 − w)(ν − 1) heute (z = 0)")
P("   R [Mpc]   σ(R)      a_N,b/a₀     Staubfeld (1+Q)a_N,b/a₀    ν(i) − 1    (1−w)(ν−1) (i)    ν(ii) − 1    (1−w)(ν−1) (ii)")
om0, x0 = omw(1.0)
for R in RS:
    s = sigR(R); g = aNb(1.0, R, s) / A0
    n1 = 1 / math.expm1(math.sqrt(g)); n2 = 1 / math.expm1(math.sqrt(g + x0))
    P(f"   {R:6.0f}   {s:7.4f}   {g:10.3e}   {g * (1 + Q_KOSM):12.3e}              {n1:8.3f}    {om0 * n1:12.3e}      {n2:10.3e}   {om0 * n2:12.3e}")
P(f"   θ̃τ_D heute = {x0:.3f}: (ii) ν − 1 ≈ 1/(e^√x − 1) = {1 / math.expm1(math.sqrt(x0)):.4e} ≈ (1 − w)/w, also (1 − w)(ν − 1) ≈ (1 − w)²/w")

hdr("2  Wachstum je Skala: relative Änderung gegen Staub mit 1 + Q (z = 0 / 0,5 / 1)")
P("   Lesart / Argument     R [Mpc]   ΔD/D z=0     Δ(fσ₈)/fσ₈ z=0   z=0,5         z=1")
res = {}
for les in ("A", "B"):
    for arg in ("i", "ii"):
        for R in RS:
            sR = sigR(R)
            m = wachstum(R, sR, arg, les)
            vals = []
            for ae in (1.0, 1 / 1.5, 0.5):
                D, Dp = m[ae]; D0, Dp0 = st[ae]
                vals.append(((D - D0) / D0, (Dp - Dp0) / Dp0))   # fσ₈ ∝ f·D = dD/dlna
            res[(les, arg, R)] = vals
            P(f"   {les} / ({arg}){'':<12s} {R:6.0f}   {vals[0][0]:+.3e}   {vals[0][1]:+.3e}       {vals[1][1]:+.3e}    {vals[2][1]:+.3e}")
        P()

# reiner Grenzfall ohne Regel: Anteil (1 − w) der Mode klumpt nicht (μ = 1 − (1 − w)Q/(1+Q)), zum Vergleich mit gewicht_rueckkehr "sofort"
def wachstum_ohne():
    def rhs(lna, y):
        a = math.exp(lna); D, Dp = y
        Om_a = OMEGA_M * a**-3 / Ea(a)**2
        return [Dp, -(2 + dlnE(a)) * Dp + 1.5 * Om_a * (1 - omw(a)[0] * Q_KOSM / (1 + Q_KOSM)) * D]
    sol = solve_ivp(rhs, (math.log(ai), 0.0), [ai, ai], rtol=1e-10, atol=1e-14, dense_output=True, method="DOP853")
    return sol.sol(0.0)
Dn, Dpn = wachstum_ohne()
P(f"   Grenzfall ν − 1 = 0 (Anteil 1 − w klumpt nicht): ΔD/D = {(Dn - st[1.0][0]) / st[1.0][0]:+.3e}, Δ(fσ₈)/fσ₈ = {(Dpn - st[1.0][1]) / st[1.0][1]:+.3e}"
      "  (gewicht_rueckkehr, w_R sofort: −4,92·10⁻³)")

# Vergleich: Argument (ii) ohne Faktor 1 − w (Verlust θ̃ nur einmal gezählt, im Argument der Regel)
def wachstum_ii_ohne_faktor(les):
    def rhs(lna, y):
        a = math.exp(lna); D, Dp = y
        Om_a = OMEGA_M * a**-3 / Ea(a)**2
        x = 3 * H0 * Ea(a) * c / a0z(a)
        nm1 = 1 / math.expm1(math.sqrt(x)) if x < 1e6 else 0.0
        m = 1 + nm1 / (1 + Q_KOSM)
        return [Dp, -(2 + dlnE(a)) * Dp + 1.5 * Om_a * m * D]
    sol = solve_ivp(rhs, (math.log(ai), 0.0), [ai, ai], rtol=1e-10, atol=1e-14, dense_output=True, method="DOP853")
    return sol.sol(0.0)
Dq, Dpq = wachstum_ii_ohne_faktor("B")
P(f"   Zusatz ν(ii) − 1 ohne Faktor 1 − w (θ̃ nur im Argument): ΔD/D = {(Dq - st[1.0][0]) / st[1.0][0]:+.3e}, Δ(fσ₈)/fσ₈ = {(Dpq - st[1.0][1]) / st[1.0][1]:+.3e}")

hdr("3  Handprüfung (von Hand, nicht mit diesem Code)")
P("   R = 10 Mpc, z = 0: H₀ = 67,4/3,0857e19 = 2,1843e-18 s⁻¹, H₀² = 4,7711e-36; Ω_b = 0,02237/0,674² = 0,04924;")
P("   R = 3,0857e23 m; ½·0,04924·4,7711e-36·3,0857e23 = 3,6246e-14 m/s² je Einheit δ.")
s10 = sigR(10.0)
P(f"   mit σ(10 Mpc) = {s10:.4f} (CAMB): a_N,b = {3.6246e-14 * s10:.4e} m/s² = {3.6246e-14 * s10 / 1.1e-10:.4e}·a₀ (Code: {aNb(1.0, 10.0, s10) / A0:.4e})")
g = 3.6246e-14 * s10 / 1.1e-10
P(f"   (i)  s = √{g:.4e} = {math.sqrt(g):.4f}; ν − 1 = 1/(e^s − 1) ≈ 1/s − ½ + s/12 = {1 / math.sqrt(g) - 0.5 + math.sqrt(g) / 12:.3f};"
  f"  mal 1 − w = 0,0146: {0.0146 * (1 / math.sqrt(g) - 0.5 + math.sqrt(g) / 12):.3f}")
P(f"   (ii) s² = 17,86 + {g:.1e}, s = 4,226, e^s = 68,45, ν − 1 = 1/67,45 = 0,01483; mal 0,0146: 2,17e-4")
P("   Lesart A (i), μ − 1 = 0,0146·(ν − 1 − 5,364)/6,364; Lesart A (ii): 0,0146·(0,0148 − 5,364)/6,364 = −1,23e-2 (Grenzfall ohne Regel −1,23e-2)")

hdr("4  Ergebnis")
P("   (i) Feld der Störung allein (a_N,b = 1–4·10⁻⁴·a₀, ν − 1 = 50–100): (1 − w)(ν − 1) = 0,7–1,5 heute; fσ₈ +3,3–7,2 % (A) bzw. +3,7–7,7 % (B):")
P("       unverträglich mit 1,2·10⁻³ (höchstens 5·10⁻³). Das Argument liest nur den Gradienten (L) und übersieht den Zustand: der Bereich atmet mit.")
P("   (ii) Summe der Verluste s² = a_N/a₀ + θ̃τ_D (θ̃τ_D = 17,9 ≫ a_N/a₀): ν − 1 = 0,0148 ≈ (1 − w)/w, (1 − w)(ν − 1) = 2,2·10⁻⁴, skalenfrei.")
P("       B (Mode bleibt voller Staub, (1 − w)·Q_s zusätzlich): fσ₈ +8,8·10⁻⁶;  A (Anteil 1 − w Operator statt Staub) = Grenzfall ohne Regel: −4,9·10⁻³.")
P("       Beides verträglich; die Übergabe λ(w) bleibt die führende Änderung (1,2·10⁻³ bei N* = 1).")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
