# Bose-Leiter, Zusatzpotential und Phantomdichte für ausgedehnte Massen
# Einheiten GM = a₀ = 1 → r* = 1, v_f = 1.  Regel a = a_N·ν(s), s = √(a_N/a₀) (kugelsymmetrisch exakt; Scheibe: Mittelebene, algebraisch)
import numpy as np, sys, io, os
from scipy.integrate import quad
from scipy.special import i0, i1, k0, k1, zeta, gamma
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
nB=lambda s: (0.0 if s>700 else 1/np.expm1(s)) if s>1e-12 else 1/max(s,1e-300)
GZ={d: gamma(d+1)*zeta(d+1) for d in (1,2,3)}

def prof(name,L):
    if name=="Punktmasse": M=lambda r: 1.0; rho=lambda r: 0.0
    if name=="Plummer":    M=lambda r: r**3/(r*r+L*L)**1.5; rho=lambda r: 3/(4*np.pi*L**3)*(1+r*r/L**2)**-2.5
    if name=="Hernquist":  M=lambda r: r*r/(r+L)**2; rho=lambda r: L/(2*np.pi*r*(r+L)**3)
    if name=="Kugel":      M=lambda r: min(r/L,1.0)**3; rho=lambda r: 3/(4*np.pi*L**3) if r<L else 0.0
    if name=="Scheibe":
        def aN(r):
            y=r/(2*L)
            if y<1e-8: return 0.0
            if y>300: return 1/r**2
            return (y*y*(i0(y)*k0(y)-i1(y)*k1(y)))/(0.5*L*L)/r*L   # v² = GM/(R_d)·2y²[I₀K₀ − I₁K₁]·… auf M = 1 normiert
        return aN, None
    return (lambda r: M(r)/r**2), rho

# Normierung der Scheibe prüfen: v² = 4πGΣ₀R_d·y²[I₀K₀−I₁K₁], M = 2πΣ₀R_d²  → v² = (2GM/R_d)·y²[…],  a_N = v²/R
def aN_disk(r,L):
    y=r/(2*L)
    if y<1e-8: return 0.0
    if y>300: return 1/r**2
    return (2/L)*y*y*(i0(y)*k0(y)-i1(y)*k1(y))/r
def get(name,L):
    if name=="Scheibe": return (lambda r: aN_disk(r,L)), None
    return prof(name,L)

def extra(aN):  # a − a_N
    return lambda r: (lambda x: x*nB(np.sqrt(x)) if x>0 else 0.0)(aN(r))
def lnr_int(f):
    t=np.linspace(-10,12,6001); v=np.array([f(np.exp(x)) for x in t]); return np.trapezoid(v,t)

P("1) Leiter im Ortsmaß: R_d = ∫(a − a_N)(r*/r)^{d−1} d ln r / (a₀Γ(d+1)ζ(d+1)),  r* aus der Gesamtmasse")
def lnr_int2(f,t0):
    t=np.linspace(t0,12,8001); v=np.array([f(np.exp(x)) for x in t]); return np.trapezoid(v,t)
P("   Profil        L/r*    d=1      d=2 (innere Grenze e⁻¹⁰ / e⁻²⁰)")
for name in ("Punktmasse","Plummer","Hernquist","Kugel","Scheibe"):
    for L in ((1.0,) if name=="Punktmasse" else (0.03,0.3,1.0,3.0)):
        aN,_=get(name,L); ex=extra(aN)
        R1=lnr_int2(lambda r: ex(r),-20)/GZ[1]
        R2a=lnr_int2(lambda r: ex(r)/r,-10)/GZ[2]; R2b=lnr_int2(lambda r: ex(r)/r,-20)/GZ[2]
        P(f"   {name:12s} {L if name!='Punktmasse' else 0:5.2f}  {R1:7.4f}   {R2a:12.4f} / {R2b:12.4f}")
P("   d = 1 hängt vom Profil ab (Kern: Faktor 0,5 bis 3,2) und divergiert logarithmisch, wo a − a_N im Zentrum endlich bleibt (Hernquist).")
P("   d = 2 und 3 divergieren für jedes ausgedehnte Profil: Das Gewicht 1/r bzw. 1/r² wächst im Zentrum ohne Grenze, a − a_N verschwindet")
P("   dort nicht exponentiell. Endlich und exakt ist die Leiter im Ortsmaß nur für die Punktmasse.")
P()

P("2) Leiter im Maß der Bose-Variablen s: ∫ s^d n_B |ds| längs des Radius")
P("   Jeder monotone Ast von s liefert ∫₀^{s_max} s^d n_B ds, unabhängig vom Profil. Ein Profil mit Kern (a_N(0) = 0) durchläuft s zweimal.")
for name,L in (("Plummer",0.3),("Hernquist",0.3),("Kugel",1.0),("Scheibe",0.3)):
    aN,_=get(name,L); r=np.exp(np.linspace(-40,14,80001)); s=np.sqrt(np.array([aN(x) for x in r]))
    smax=s.max(); ds=np.abs(np.diff(s)); sm=0.5*(s[1:]+s[:-1])
    for d in (1,2):
        num=np.sum(sm**d/np.expm1(np.maximum(sm,1e-300))*ds)
        inc=quad(lambda x: x**d/np.expm1(x),0,smax)[0]
        branches=2 if s[0]<0.5*smax else 1
        P(f"   {name:10s} L = {L}: s_max = {smax:6.3f}, Äste {branches},  d = {d}: längs r {num:.5f}   Äste × ∫₀^s_max = {branches*inc:.5f}   (Γζ = {GZ[d]:.5f})")
P()

P("3) Zusatzpotential Φ_X(r) = ∫₀^r (a − a_N)dr gegen die Punktmassenform v_f²·ln ν(s(r)):")
P("   Profil        L/r*   Φ_X − v_f² ln ν bei r = 3L, 10 r*, 100 r*  (Einheiten v_f²)")
for name in ("Plummer","Hernquist","Kugel","Scheibe"):
    for L in (0.03,0.3,1.0,3.0):
        aN,_=get(name,L); ex=extra(aN)
        vals=[]
        for R_ in (3*L,10.0,100.0):
            t=np.linspace(-12,np.log(R_),6001); x=np.exp(t); v=np.array([ex(xx)*xx for xx in x]); Phi=np.trapezoid(v,t)
            s=np.sqrt(aN(R_)); vals.append(Phi+np.log(-np.expm1(-s)))
        P(f"   {name:12s} {L:5.2f}   {vals[0]:+8.4f}  {vals[1]:+8.4f}  {vals[2]:+8.4f}")
P("   Außen bleibt eine feste Verschiebung C: Φ_X = v_f²·ln ν(s) + C; C → 0 für kompakte Profile.")
P()

P("4) Phantomdichte in der Masse (kugelsymmetrisch, exakt):")
P("   4πGρ_X = χΩ² + 4πGρ_b·∂(a − a_N)/∂a_N,   ∂(a − a_N)/∂a_N = n_B(1 − sν/2): Vorzeichenwechsel bei sν = 2 ⇔ 2(1 − e^{−s}) = s")
s2=brentq(lambda s: 2*(1-np.exp(-s))-s,1,2); P(f"   Vorzeichenwechsel bei s = {s2:.4f}, a_N = {s2**2:.3f}·a₀: genau die Wien-Stelle der Stufe 2 (größte Zusatzbeschleunigung)")
L=0.3; aN,rho=get("Plummer",L)
for r in (0.05,0.2,0.5,1.0,3.0):
    h=1e-5*r
    F=lambda x: x*x*(aN(x)*nB(np.sqrt(aN(x))))
    num=(F(r+h)-F(r-h))/(2*h)/r**2
    s=np.sqrt(aN(r)); nu=1+nB(s); chi=s*nB(s); Om2=aN(r)*nu/r
    form=chi*Om2+4*np.pi*rho(r)*nB(s)*(1-s*nu/2)
    P(f"   Plummer L = 0,3, r = {r}: numerisch {num:.6f}   Formel {form:.6f}   Materieglied {4*np.pi*rho(r)*nB(s)*(1-s*nu/2):+.5f}  (a_N/a₀ = {aN(r):.3f})")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"ausgedehnte_massen.txt"),"w",encoding="utf-8").write(txt); print(txt)
