# Ist Q̃ der nichtholomorphe Teil (Schatten-Seite) einer harmonischen Maaß-Form? Und was fällt an Q̃ und ν auf?
# Q̃(s) = −4Σ Γ(4,ns)/n⁴, s = 2π·Im τ (Wahl Im τ_M = s/2π)
import numpy as np, sys, io, os
from mpmath import mp, mpf, gammainc, nsum, inf, exp, diff, findroot
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
mp.dps = 30
out=[]; P=lambda s="": out.append(s)
Qt=lambda s: -4*nsum(lambda n: gammainc(4,n*s)/n**4,[1,inf])
nB=lambda s: 1/(exp(s)-1)

P("1) Bauart: nichtholomorpher Teil einer harmonischen Maaß-Form vom Gewicht k")
P("   F⁻ = Σ_{n≥1} c(n)·Γ(1−k, 4πny)·q^{−n};  auf Re τ = 0:  Σ c(n)·Γ(1−k, 4πny)·e^{2πny}")
P("   Q̃ auf Re τ = 0 (s = 2πy):  −4Σ Γ(4, 2πny)/n⁴.  Grad des Polynoms in Γ(a,x) = e^{−x}·Poly_{a−1}(x) erzwingt 1 − k = 4, k = −3.")
P("   Je Exponent e^{−2πny} müsste gelten: c(n)·[(4πny)³ + 3(4πny)² + 6(4πny) + 6] = d(n)·[(2πny)³ + 3(2πny)² + 6(2πny) + 6]")
P("   y³: c = d/8,  y²: c = d/4,  y¹: c = d/2,  y⁰: c = d  → unvereinbar für d ≠ 0.")
# numerische Probe: beste Anpassung c für n = 1 über y
y=np.linspace(0.05,2,200); n=1
lhs=np.exp(-2*np.pi*n*y)*((4*np.pi*n*y)**3+3*(4*np.pi*n*y)**2+6*(4*np.pi*n*y)+6)
rhs=np.exp(-2*np.pi*n*y)*((2*np.pi*n*y)**3+3*(2*np.pi*n*y)**2+6*(2*np.pi*n*y)+6)
c=np.dot(lhs,rhs)/np.dot(lhs,lhs); res=np.max(np.abs(c*lhs-rhs))/np.max(rhs)
P(f"   beste Konstante c für n = 1: {c:.4f}, verbleibender Fehler {res:.1%} des Maximums → keine Gleichheit")
P("   Zweitens: Q̃ hängt nicht von Re τ ab; harmonisch vom Gewicht k ist eine solche Funktion nur als a + b·y^{1−k}.")
P("   Drittens: Der Schatten eines Gewichts −3 hätte Gewicht 2 − k = 5; auf SL(2,ℤ) gibt es keine Formen ungeraden Gewichts")
P("   (−I wirkt mit (−1)^k: f = −f). Ergebnis: Q̃ ist nicht der nichtholomorphe Teil einer Mock-Form auf SL(2,ℤ).")
P()

P("2) Was gilt: der ξ-Operator verbindet Q̃ und die Rückkehrsumme")
P("   ξ_k F = 2i·y^k·conj(∂F/∂τ̄); für F(y) reell: ξ_k F = y^k·F'(y).  Mit k = −3 und s = 2πy:")
for s in (mpf('0.3'),mpf(1),mpf(3)):
    yv=s/(2*mp.pi)
    lhs=yv**(-3)*diff(lambda yy: Qt(2*mp.pi*yy), yv)
    rhs=4*(2*mp.pi)**4*nB(s)
    P(f"   s = {float(s)}:  y⁻³·dQ̃/dy = {float(lhs):.10f}   4(2π)⁴·n_B(s) = {float(rhs):.10f}")
P("   ξ₋₃ Q̃ = 4(2π)⁴·(ν − 1) auf der imaginären Achse: Strukturwirkung und Rückkehrsumme stehen wie nichtholomorpher Teil")
P("   und Schatten zueinander; der „Schatten“ Σqⁿ ist aber keine Modulform.")
P()

P("3) Planck-Muster: Zusatzbeschleunigung und Strukturwirkung")
P("   a − a_N = a_N·n_B(s) = a₀·s²·n_B(s)       Form des Planck-Spektrums in 2 Dimensionen (Dichte s, Energie s)")
P("   dQ̃/ds = 4s³·n_B(s) = 4s·(a − a_N)/a₀     Form des Planck-Spektrums in 3 Dimensionen")
for d in (2,3):
    sw=findroot(lambda s: d*(1-exp(-s))-s, 1.5 if d==2 else 2.8)
    P(f"   Maximum von s^{d}·n_B (Wien-Stelle, d = {d}): s = {float(sw):.4f},  a_N/a₀ = s² = {float(sw**2):.3f},  Wert s^{d}n_B = {float(sw**d*nB(sw)):.4f}")
sw=findroot(lambda s: 2*(1-exp(-s))-s, 1.5)
P(f"   → Die Zusatzbeschleunigung a − a_N ist am größten bei a_N = {float(sw**2):.3f}·a₀ und beträgt dort {float(sw**2*nB(sw)):.4f}·a₀;")
P(f"     Punktmasse bei r = r*/{float(sw):.4f} = {float(1/sw):.3f}·r* (Milchstraße 8,35 kpc: {float(8.35/sw):.2f} kpc).")
P("   Grenzen: a − a_N → a₀·s (s → 0, also √(a_N a₀)) und → a_N·e^{−s} (s → ∞).")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"mock_schatten.txt"),"w",encoding="utf-8").write(txt); print(txt)
