# Wiederkehrende Muster quer durch das Modell und die Rolle der Dimension
# A Bose-Leiter s^d·n_B(s), d = 0…3   B Integrale Γ(d+1)ζ(d+1) als Modellgrößen   C Wien-Stellen je Dimension
# D die 3 in Δs als Dimension (Überschussradius der d-Kugel)   E Übergangslängen als geometrische Mittel
import numpy as np, sys, io, os
from mpmath import mp, mpf, quad, inf, zeta, gamma, exp, findroot, nsum, gammainc, sqrt, pi
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
mp.dps = 25
out=[]; P=lambda s="": out.append(s)
nB=lambda s: 1/(exp(s)-1)
nu=lambda s: 1/(1-exp(-s))
chi=lambda s: s/(exp(s)-1)
Qt=lambda s: -4*quad(lambda x: x**3*nB(x),[s,inf])

P("A) Bose-Leiter: die Hauptgrößen der Regel sind s^d·n_B(s) für d = 0, 1, 2, 3  (s = √(a_N/a₀))")
for s in (mpf('0.3'),mpf(1),mpf(3)):
    a_N=s**2; a=a_N*nu(s)   # in Einheiten a₀
    dQ=4*s**3*nB(s)
    P(f"   s = {float(s)}:  d=0 n_B = {float(nB(s)):.6f} = ν−1 = {float(nu(s)-1):.6f} = M_schein/M − 1 | d=1 s·n_B = {float(s*nB(s)):.6f} = χ = {float(chi(s)):.6f}")
    P(f"            d=2 s²n_B = {float(s**2*nB(s)):.6f} = (a−a_N)/a₀ = {float(a-a_N):.6f} | d=3 s³n_B = {float(s**3*nB(s)):.6f} = ¼·dQ̃/ds = {float(dQ/4):.6f}")
P("   Jede Stufe ist die vorige mal s. d = 1 ist zugleich der Verlust an Flussdimension: d_eff = 3 − s·n_B.")
P()

P("B) Integrale ∫₀^∞ s^d n_B ds = Γ(d+1)ζ(d+1) und wo sie im Modell stehen")
for d in (1,2,3):
    I=quad(lambda s: s**d*nB(s),[0,inf]); P(f"   d = {d}:  {float(I):.10f}  = Γ({d+1})ζ({d+1}) = {float(gamma(d+1)*zeta(d+1)):.10f}")
# d=1: ∫(a − a_N) d ln r für die Punktmasse, s = r*/r
I1=quad(lambda lr: (lambda s: s**2*nB(s))(exp(-lr)),[-40,40])
P(f"   d = 1:  ∫(a − a_N) d ln r = {float(I1):.10f}·a₀ = (π²/6)·a₀ für jede Punktmasse; derselbe Wert ist der Boden (π²/6)a₀ aller Harmonischen")
s_=mpf('0.02'); L=nsum(lambda n: sum(dd for dd in range(1,int(n)+1) if int(n)%dd==0)*exp(-n*s_),[1,4000])
P(f"           Probe Boden: s²(1 + Σσ₁(n)e^{{−ns}}) bei s = 0,02: {float(s_**2*(1+L)):.5f} → π²/6 = {float(pi**2/6):.5f}")
P(f"   d = 3:  Q̃(0) = {float(Qt(mpf(0))):.10f} = −4π⁴/15 = {float(-4*pi**4/15):.10f}: der Vakuumteil der Dunklen Energie")
I2=quad(lambda s: s**2*nB(s),[0,inf])
P(f"   d = 2:  ∫(a − a_N) d(r*/r) = {float(I2):.6f}·a₀ = 2ζ(3)·a₀: im Modell noch ohne Größe (gezielt zu suchen)")
P()

P("C) Wien-Stellen: Maximum von s^d·n_B, d(1 − e^{−s}) = s")
for d in (1,2,3):
    if d==1:
        P("   d = 1: Maximum am Rand s → 0 (χ → 1, voller Verlust einer Flussdimension)")
        continue
    sw=findroot(lambda s: d*(1-exp(-s))-s, 1.5 if d==2 else 2.8)
    P(f"   d = {d}: s = {float(sw):.4f}, a_N/a₀ = {float(sw**2):.3f}, r = {float(1/sw):.3f}·r* (Punktmasse), Wert {float(sw**d*nB(sw)):.4f}")
P()

P("D) Die 3 in Δs = GM/(3c²) ist die Raumdimension")
for d in (1,2,3,4,5):
    I=quad(lambda x: x**(d-1),[0,1])      # ∫ m(r)/r dr mit m ∝ r^d, R = 1
    P(f"   homogene {d}-Kugel: Überschussradius ∫Gm(r)/(rc²)dr = GM/c² · {float(I):.4f} = GM/({d}c²);  V/(S·R) = 1/{d}")
P("   Δs = (GM/c²)·V/(S·R): Faltungstiefe = Schwarzschild-Länge mal Volumen durch (Oberfläche × Radius).")
P("   Damit werden die Zahlen des Modells zu Vielfachen der Dimension d = 3:")
for name,model,gr in (("Takt δ = 3Δs/r","3","GM/(rc²): 1"),("Grenzradius r₀ = 6Δs","6 = 2d","2GM/c²: 2"),("Ablenkung 12Δs/b","12 = 4d","4GM/(bc²): 4"),("Photonensphäre 9Δs","9 = 3d","3GM/c²: 3")):
    P(f"   {name:24s} Faktor {model:8s}  ART {gr}")
P()

P("E) Übergangslängen als geometrisches Mittel aus Faltung und Gegenlänge")
G=6.674e-11; c=2.99792458e8; hbar=1.054571817e-34; a0=1.2e-10; mp_=1.67262192e-27; Ms=1.989e30
lP=np.sqrt(hbar*G/c**3); l0=c**2/a0
for name,M in (("Proton",mp_),("Sonne",Ms),("Milchstraße",6e10*Ms)):
    ds=G*M/(3*c**2); lam=hbar/(M*c)
    P(f"   {name:12s}: √(3Δs·λ̄) = {np.sqrt(3*ds*lam):.4e} m = ℓ_P = {lP:.4e} m (fest)   √(3Δs·ℓ₀) = {np.sqrt(3*ds*l0):.4e} m = r* (∝ √M)")
P("   Klein: Faltung × Quantenlänge = ℓ_P² (Übergang Faltung = Quantenlänge).  Groß: Faltung × Schwellenlänge = r*² (Übergang Newton = Schwelle).")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"wiederkehrende_muster.txt"),"w",encoding="utf-8").write(txt); print(txt)
