# Kern auf der exakten Kreisbahn als Schleife mit Verstärkung G = κ0 M/(λ²+ω²); Realteil von τ_M
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from scipy.integrate import quad
from scipy.optimize import brentq
G_=6.674e-11; a0=1.2e-10; Msun=1.989e30
GM=6e10*Msun*G_; rs=np.sqrt(GM/a0)
out=[]; P=lambda s="": out.append(s)

P("1) Gedächtnisintegral numerisch, Kreisbahn, Gedächtnis in Umlaufphase: λ_m = k·ω, κ0 = (1+k²)·ω², M = q = e^{−r*/r}")
for k in [0.3,1.0,3.0]:
    P(f"   k = {k}")
    for x in [0.1,0.5,1,3,10,100]:
        r=x*rs; gN=GM/r**2; q=np.exp(-1/x)
        def resid(v):
            w=v/r; lam=k*w; kap=(1+k*k)*w*w
            I=quad(lambda u: np.exp(-k*u)*np.sin(u),0,80/k,limit=2000)[0]/w   # u = ω s; e_r(t)·e_φ(t') = sin(ω s)
            aK=kap*q*r*w*I
            return v*v/r-gN-aK
        vrar=np.sqrt(r*gN/(1-q))
        v=brentq(resid,0.5*np.sqrt(gN*r),3*vrar)
        P(f"     r/r* = {x:7.2f}   v = {v/1e3:9.4f} km/s   v_RAR = {vrar/1e3:9.4f}   Abw. {v/vrar-1:+.1e}")

P(); P("2) Realteil Re τ_M = κ/Ω aus der Kurve gegen die geschlossene Form √(1 + s/(e^s−1)), s = r*/r")
P("     r/r*     κ/Ω numerisch   √(1+s/(e^s−1))   Perizentrenabstand   Präzession je Radialperiode")
vf=lambda r: np.sqrt(GM/(r*(1-np.exp(-rs/r))))
for x in [0.1,0.3,1,3,10,100]:
    r=x*rs; h=1e-5
    b=(np.log(vf(r*(1+h)))-np.log(vf(r*(1-h))))/(np.log(1+h)-np.log(1-h))
    kn=np.sqrt(2*(1+b)); s=1/x; kc=np.sqrt(1+s/np.expm1(s))
    P(f"   {x:7.2f}     {kn:.8f}      {kc:.8f}        {360/kc:8.3f}°            {360/kc-360:+8.3f}°")
P(); P("3) Reihe im tiefen Bereich: Re τ_M² = 1 + Σ B_n s^n/n! = 2 − s/2 + s²/12 − s⁴/720 + …,  s = 2π Im τ_M")
for s in [0.1,0.5,1.0]:
    P(f"   s = {s}: exakt {1+s/np.expm1(s):.8f}, Reihe bis s⁴ {2-s/2+s*s/12-s**4/720:.8f}")
P(); P("4) Sonnensystem: s = r*/r für die Erde = 7031 → Re τ_M − 1 = ½·7031·e^{−7031} ≈ 10^%d" % int(np.log10(0.5*7031)-7031/np.log(10)))
open("kern_schleife.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))

# 5) Allgemeine Massenverteilung: (Re τ_M)² = 3 + γ_N·(1 − s(ν−1)/2), γ_N = d ln g_N/d ln r, s = √(g_N/a0)
out=[]
P("5) Allgemeine Verteilung (Hernquist-Kugel, M = 6e10 M☉, a = 3 kpc): (Re τ_M)² = 3 + γ_N·(1 − s(ν−1)/2)")
kpc=3.0857e19; ah=3*kpc
gN=lambda r: GM/(r+ah)**2
def v(r):
    s=np.sqrt(gN(r)/a0); return np.sqrt(r*gN(r)/(1-np.exp(-s)))
P("     R [kpc]   κ/Ω numerisch   Formel")
for R in [1,3,8,15,30,60]:
    r=R*kpc; h=1e-5
    b=(np.log(v(r*(1+h)))-np.log(v(r*(1-h))))/(np.log(1+h)-np.log(1-h))
    g=(np.log(gN(r*(1+h)))-np.log(gN(r*(1-h))))/(np.log(1+h)-np.log(1-h))
    s=np.sqrt(gN(r)/a0); nu=1/(1-np.exp(-s))
    P(f"     {R:5.0f}     {np.sqrt(2*(1+b)):.8f}      {np.sqrt(3+g*(1-s*(nu-1)/2)):.8f}")
open("kern_schleife.txt","a",encoding="utf-8").write("\n\n"+"\n".join(out)); print("\n".join(out))
