# Die Größe der Stufe d = 2 (Wert 2ζ(3)) und das Zusatzpotential in geschlossener Form
# Punktmasse, s = r*/r, a_N = a₀s², a = a_N/(1 − e^{−s}), v_f² = a₀r* = √(GMa₀), Ω² = a/r, Ω_N(r*)² = a₀/r*
import sys, io, os
from mpmath import mp, mpf, quad, inf, zeta, exp, log, pi, gamma
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
mp.dps = 25
out=[]; P=lambda s="": out.append(s)
nu=lambda s: 1/(1-exp(-s))
# Einheiten: a₀ = 1, r* = 1  → v_f² = 1, Ω_N(r*)² = 1
a =lambda r: (1/r**2)*nu(1/r)
aN=lambda r: 1/r**2
lnr=lambda f: quad(lambda t: f(exp(t)), [-12,-3,0,3,12,40])   # ∫ f(r) d ln r

P("1) Jede Stufe, integriert mit dem Maß ihres Flussgesetzes (Gewicht (r*/r)^{d−1} auf d ln r):")
P("   ∫(a − a_N)·(r*/r)^{d−1} d ln r = a₀·Γ(d+1)ζ(d+1)")
for d in (1,2,3):
    I=lnr(lambda r: (a(r)-aN(r))*(1/r)**(d-1))
    P(f"   d = {d}:  {float(I):.10f}   Γ({d+1})ζ({d+1}) = {float(gamma(d+1)*zeta(d+1)):.10f}")
P("   d = 1: Maß d ln r (Skala);  d = 2: Maß d(1/r) (Umfang, Krümmung der Bahnkreise);  d = 3: Maß d(1/r²) (Fläche, Newton-Feld)")
P()

P("2) Die Größe der Stufe 2 in Bahnfrequenzen: (a − a_N)·r*/r = r*·(Ω² − Ω_N²)")
I2=lnr(lambda r: (a(r)-aN(r))/r)
P(f"   ∫(Ω² − Ω_N²) d ln r = {float(I2):.10f} · a₀/r*  = 2ζ(3)·Ω_N(r*)²,  2ζ(3) = {float(2*zeta(3)):.10f}")
P("   Die zusätzliche Umlauffrequenz im Quadrat, über den Logarithmus des Radius summiert, ist 2ζ(3) mal das Quadrat")
P("   der Kepler-Frequenz am Schwellenradius; für jede Masse dieselbe Zahl.")
G=6.674e-11; a0=1.2e-10; Ms=1.989e30; Myr=3.156e13
for name,M in (("Sonne",Ms),("Milchstraße",6e10*Ms)):
    rs=(G*M/a0)**0.5; Om=(a0/rs)**0.5
    P(f"   {name:12s}: Ω_N(r*) = {Om:.3e} s⁻¹ (Umlauf {2*3.14159265/Om/Myr:.3g} Mio. Jahre),  ∫(Ω² − Ω_N²) d ln r = {float(2*zeta(3))*Om**2:.3e} s⁻²")
P()

P("3) Das Zusatzpotential in geschlossener Form (vom Zentrum gezählt):")
P("   Φ_X(r) = ∫₀^r (a − a_N) dr = a₀r*·∫_s^∞ n_B ds = v_f²·ln ν(s) = −v_f²·ln(1 − e^{−s})")
for r in (mpf('0.3'),mpf(1),mpf(3),mpf(30)):
    num=quad(lambda x: a(x)-aN(x),[0,r]); P(f"   r = {float(r):5.1f} r*:  numerisch {float(num):.10f}   v_f²·ln ν = {float(log(nu(1/r))):.10f}")
P("   ln ν = −ln(1 − e^{−s}) = ln Z ist der Logarithmus der Zustandssumme einer Bose-Mode; außen Φ_X → v_f²·ln(r/r*).")
P("   Takt: δ_X = (v_f/c)²·ln ν(s).")
P("   Dieselbe Zahl 2ζ(3) als Potential: ∫₀^∞ Φ_X d(a_N/a₀) = 2·v_f²∫ s·ln ν ds = 2ζ(3)·v_f²:")
I3=quad(lambda s: 2*s*log(nu(s)),[0,inf]); P(f"   {float(I3):.10f}")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"stufe_zwei.txt"),"w",encoding="utf-8").write(txt); print(txt)
