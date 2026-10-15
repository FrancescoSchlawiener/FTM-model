# Feldgleichung des Autors: a = a_N + χ·(Δs/r)·(3c²/r*), Kern-Nome q_Δs = Δs/r, χ glatt
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19
out=[]; P=lambda s="": out.append(s)
chi=lambda s: s*np.exp(-s)/(-np.expm1(-s))
P("1) Punktmasse 6e10 M☉: a = a_N + χ(s)·(Δs/r)(3c²/r*),  χ(s) = s/(e^s−1),  s = √(a_N/a0) = r*/r")
GM=6e10*Msun*G; rs=np.sqrt(GM/a0); ds=GM/(3*c**2)
P("     r/r*     χ          a/a_N        RAR 1/(1−e^{−s})   Im τ_M = ln(r/Δs)/2π")
for x in [0.1,0.3,1,3,10,100]:
    r=x*rs; aN=GM/r**2; s=np.sqrt(aN/a0)
    a=aN+chi(s)*(ds/r)*(3*c**2/rs)
    P(f"   {x:7.2f}   {chi(s):.6f}   {a/aN:.10f}   {1/(1-np.exp(-s)):.10f}      {np.log(r/ds)/(2*np.pi):.4f}")
P()
P("2) Beliebige Verteilung (Hernquist, a = 3 kpc): a = a_N + χ(s)·√(a_N a0)")
ah=3*kpc
for R in [1,3,8,15,30,60]:
    r=R*kpc; aN=GM/(r+ah)**2; s=np.sqrt(aN/a0)
    a=aN+chi(s)*np.sqrt(aN*a0)
    P(f"   R = {R:3d} kpc   a/a_N = {a/aN:.10f}   RAR {1/(1-np.exp(-s)):.10f}")
P()
P("3) Grenzfälle von χ: a_N ≫ a0 → χ ≈ s·e^{−s} → 0 ;  a_N ≪ a0 → χ ≈ 1 − s/2 → 1")
for s in [7031,10,3,1,0.3,0.1,0.01]:
    P(f"   s = {s:8}:  χ = {chi(s):.6e}   (Re τ_M)² = 1+χ (Punktmasse) = {1+chi(s):.6f}")
P()
P("4) Sonnensystem (Erde, s = 7031): Zusatz χ·√(a_N a0) = 7031·e^(−7031)·√(a_N a0) ≈ 10^%d m/s²" % int(np.log10(7031*np.sqrt(5.93e-3*a0))-7031/np.log(10)))
open("kern_chi.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))
