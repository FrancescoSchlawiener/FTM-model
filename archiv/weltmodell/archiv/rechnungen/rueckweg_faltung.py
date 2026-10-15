# Rückweg durch die Faltung: q = e^{−√(a_N/a0)} als Green-Funktion eines abgeschirmten Operators D∂²−Γ
# in einer inneren Richtung w. Prüft: (1) was Strukturwirkung/Faltungsgleichung in w liefern, (2) warum Γ ∝ a_N,
# (3) warum offen, (4) Folgen: Newton, Tully-Fisher, Unabhängigkeit von M, prüfbare Aussage.
import numpy as np, sys, io, os
import sympy as sp, mpmath as mp
from scipy import linalg, optimize, special
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; AU=1.496e11; Gyr=3.156e16
H0=70e3/3.0857e22; l0=c**2/a0
out=[]; P=lambda s="": out.append(s)
nu=lambda s: 1/(-np.expm1(-s)); chi=lambda s: s/np.expm1(s)
rng=np.random.default_rng(20260929)

# ---------------- 0) Grundform ----------------
P("0) Faser w ≥ 0, D·G'' = Γ·G, G(0) = 1, Ende bei w = W reflektierend: G(L) = cosh(m(W−L))/cosh(mW), m = √(Γ/D), s = mL")
P("   W → ∞ (offen): G(L) = e^{−s};  W = L (Spiegel am Start = geschlossener Kreis 2L): 1/cosh s. Finite Differenzen (h = 2,5·10⁻⁵ L bzw. 5·10⁻⁴ L):")
def fd(s,W,N=40000):
    w=np.linspace(0,W,N+1); h=w[1]; m2=s*s; n=N
    ab=np.zeros((3,n)); ab[0,1:]=1; ab[2,:-1]=1; ab[1,:]=-2-m2*h*h; ab[2,-2]=2  # Neumann am Ende (Spiegelpunkt)
    b=np.zeros(n); b[0]=-1
    g=linalg.solve_banded((1,1),ab,b); return np.interp(1.0,w[1:],g)
P("      s      W/L = 1 (FD / 1/cosh s)     W/L = 200 (FD / e^{−s})")
for s in [0.1,1,3]:
    P(f"   {s:5.1f}    {fd(s,1):.6f} / {1/np.cosh(s):.6f}        {fd(s,200,400000):.6f} / {np.exp(-s):.6f}")
P("   Tiefgrenzfall s·ν bei s = 0,01 (Tully-Fisher verlangt 1): Ende bei W = ")
row=""
for W in [1,10,100,1e3,1e4,np.inf]:
    q=np.exp(-0.01) if np.isinf(W) else np.cosh(0.01*(W-1))/np.cosh(0.01*W)
    row+=f"  {W:g}: {0.01/(1-q):.3f}"
P("   "+row+"    → offen heißt: kein Ende vor W ≈ L/s.")
P()

# ---------------- 1) Wirkung und Faltungsgleichung in einer inneren Richtung ----------------
P("1) Strukturwirkung ℒ = 2h^{μν}B_μB_ν + (2/ℓ0²)Q̃(q), q = exp(−(ℓ0²Y)^{1/4}), Y = h^{μν}∂σ∂σ: σ tritt nur über ∂σ auf (σ → σ + konst).")
Y,l,sy=sp.symbols('Y ell0 s',positive=True)
s_of_Y=(l**2*Y)**sp.Rational(1,4); nuY=1/(1-sp.exp(-s_of_Y))
ratio=sp.simplify((2*Y*sp.diff(nuY,Y)/nuY).subs(Y,sy**4/l**2))
P("   Linearisierung von ∇·[ν∇σ] (Faltungsgleichung) um ein statisches Gefälle n = ∇σ/|∇σ|:")
P("   μ_ij = ν·δ_ij + 2Y·ν'(Y)·n_i n_j,   2Yν'/ν = "+str(sp.simplify(ratio))+"  =  −χ/2:  "+
  str(sp.simplify(ratio+sy/(2*(sp.exp(sy)-1)))==0))
def jac(sv,eps=1e-7):
    g0=sv**2  # |∇σ| in Einheiten 1/ℓ0, s² = ℓ0|∇σ|
    F=lambda gx,gy: nu(np.sqrt(np.hypot(gx,gy)))*np.array([gx,gy])
    Jxx=(F(g0*(1+eps),0)[0]-F(g0*(1-eps),0)[0])/(2*g0*eps); Jyy=(F(g0,g0*eps)[1]-F(g0,-g0*eps)[1])/(2*g0*eps)
    return Jxx/Jyy
P("   numerisch (Jacobi-Matrix des Flusses): μ_∥/μ_⊥ = "+", ".join(f"s={sv}: {jac(sv):.6f} (1−χ/2 = {1-chi(sv)/2:.6f})" for sv in [0.5,1,3]))
P("   μ_∥/μ_⊥ = 1 − χ/2 = d ln a/d ln a_N (Steigung der Beschleunigungsrelation).")
P("   Mit innerer Richtung w (h^{ww} = 1, a_w = 0), schwaches Feld: δσ = δΦ_N/c² folgt aus der Gleichung ∇²Φ_N = 4πGρ, im Leeren:")
P("   (∇² + ∂_w²)δσ = 0,  δσ = f(w)·e^{ik·x}:  f'' = k²f,  beschränkt für w → ∞: f = e^{−kw} (Poisson-Kern, zweiter Ordnung, offen).")
P("   Die Faltungsgleichung überträgt δσ auf die Neigung: ∂_w²δΦ − k²δΦ = ν f'' − (μ_⊥k_⊥² + μ_∥k_∥²) f; kein Glied ohne Ableitung.")
P("   Faser an einem Punkt (k → 0): f'' = 0, masselos, q = 1 – kein Verlust.")
P("   Kreisbahn als reelle Schleife (Umfang 2πr), niedrigste Mode e^{iφ} (senkrecht zum Gefälle, μ_⊥ = ν): f = e^{−w/r}, δΦ ebenso.")
P("   → Die Wirkung liefert eine abgeschirmte Gleichung in w mit √(Γ/D) = 1/r, aber keinen Verlust ∝ a_N an einem Punkt.")
P("   q = e^{−s} verlangt dann die innere Länge L = r·s = r·√(a_N/a0); Punktmasse: L = r* = √(GM/a0) = √(3Δs·ℓ0), längs r fest.")
P()

# ---------------- 2) L, D, Γ in Modellgrößen ----------------
P("2) Das Argument von q in Q̃: s² = √(ℓ0²Y) = ℓ0·|∇σ| = (ℓ0/c)·(c|∇σ|) = τ_D·Γ.  Γ = c|∇σ| = a_N/c,  τ_D = L²/D = ℓ0/c.")
P("   Starkfeld σ = ln(dτ/dt): Γ = c|∇ln√(1−6Δs/r)| = (a_N/c)/(1 − 6Δs/r) – das Taktgefälle der Faltung, unendlich am Knoten r0 = 6Δs.")
P("   Satz mit L = r*: D = L²c/ℓ0 = GM/c = 3cΔs,  Γ/D = a_N/(GM) = 1/r² (Punktmasse) = Abschirmung der Kreisbahnmode aus 1).")
for name,M in [("Sonne",Msun),("Milchstraße 6·10¹⁰ M☉",6e10*Msun)]:
    GM=G*M; ds=GM/(3*c**2); rs=np.sqrt(GM/a0); D=GM/c
    P(f"   {name:22s} Δs = {ds:.4e} m  r* = {rs:.4e} m = √(3Δs·ℓ0) {np.sqrt(3*ds*l0):.4e}  D = {D:.3e} m²/s  "
      f"Γ(r*) = {a0/c:.3e} s⁻¹  τ_D = {rs**2/D/Gyr:.2f} Mrd. J.  freie Weglänge 2D/c = {2*D/c:.4e} m = 6Δs {6*ds:.4e}")
P(f"   τ_D = c/a0 = {c/a0/Gyr:.2f} Mrd. Jahre; mit a0 = cH0/2π (H0 = 70) ist τ_D = 2π/H0 = {2*np.pi/H0/Gyr:.2f} Mrd. Jahre.")
P()

# ---------------- 3) Warum Γ ∝ a_N ----------------
P("3) Andere Raten bei diffusivem Rückweg (Punktmasse, s so geeicht, dass s = 1 bei r* der Milchstraße):")
Mref=6e10*Msun; rsref=np.sqrt(G*Mref/a0)
def acc(kind,M,r):
    aN=G*M/r**2
    if kind=="a_N/c": s=np.sqrt(aN/a0)
    elif kind=="√(a_N a0)/c": s=(aN/a0)**0.25
    elif kind=="c|a| (a physikalisch)":
        xN=aN/a0; f=lambda x: x-xN*nu(np.sqrt(x)); return a0*optimize.brentq(f,xN,xN+2*np.sqrt(xN)+10,xtol=1e-15,rtol=1e-13)
    elif kind=="ω = √(a_N/r)": s=np.sqrt(np.sqrt(aN/r)*np.sqrt(rsref/a0))
    return aN*nu(s)
P("      Rate Γ                     d ln v/d ln r (30–100 r*)   d ln v/d ln M bei r = 200 r*(M_MW)   v⁴/(GMa0) bei 100 r*")
for kind in ["a_N/c","√(a_N a0)/c","c|a| (a physikalisch)","ω = √(a_N/r)"]:
    v=lambda M,r: np.sqrt(acc(kind,M,r)*r)
    sl=np.log(v(Mref,100*rsref)/v(Mref,30*rsref))/np.log(100/30)
    sm=np.log(v(1e12*Msun,200*rsref)/v(1e10*Msun,200*rsref))/np.log(100)
    P(f"   {kind:26s}  {sl:+.4f}                      {sm:.4f}                              {v(Mref,100*rsref)**4/(G*Mref*a0):.4f}")
P("   Nur Γ ∝ a_N gibt flache Kurven mit v⁴ = GMa0 (Steigung 0, Massenexponent 1/4).")
P()

# ---------------- 4) Warum offen: Umlaufsinn auf der Überlagerung ----------------
P("4) Irrfahrt auf einer geschlossenen inneren Schleife (Umfang L), Verlust mit Rate Γ; Lage auf der Überlagerung (Gerade) mitgeführt.")
P("   A: erste Windung +1 (Punkt +L) vor dem Verlust;  B: erste Windung ±1 (Austritt aus (−L, L)).  D = 1/2, L = 1, Γ = s²/2, 50 000 Läufer:")
def mc(s,N=50000,dt=2e-3):
    Gm=s*s/2; Tk=rng.exponential(1/Gm,N); x=np.zeros(N); A=np.zeros(N,bool); B=np.zeros(N,bool)
    idx=np.arange(N); t=0.0; sq=np.sqrt(dt)
    while idx.size:
        xo=x[idx]; xn=xo+sq*rng.standard_normal(idx.size)
        pA=np.exp(-2*np.clip(1-xo,0,None)*np.clip(1-xn,0,None)/dt)
        pB=np.exp(-2*np.clip(xo+1,0,None)*np.clip(xn+1,0,None)/dt)
        hitA=(xn>=1)|(rng.random(idx.size)<pA); hitB=hitA|(xn<=-1)|(rng.random(idx.size)<pB)
        B[idx]|=hitB; A[idx]|=hitA; x[idx]=xn; t+=dt
        keep=(~hitA)&(Tk[idx]>t); idx=idx[keep]
    return A.mean(),B.mean()
P("      s      A (MC)    e^{−s}     B (MC)    1/cosh s")
for s in [0.5,1,2]:
    a,b=mc(s); P(f"   {s:5.1f}    {a:.4f}    {np.exp(-s):.4f}     {b:.4f}    {1/np.cosh(s):.4f}")
P("   Dieselbe geschlossene Schleife gibt e^{−s}, wenn nur Windungen in einem Sinn Rückläufe sind, und 1/cosh s, wenn beide zählen.")
s=0.01
P(f"   Summe der Windungen beider Sinne (Green-Funktion der Schleife): Σ_n e^{{−|n|s}} = coth(s/2) = 2ν − 1; s·coth(s/2) bei s = 0,01: {s/np.tanh(s/2):.4f} → v⁴ = 4GMa0.")
P("   Drift Pe = vL/2D in der Faser: q = exp(−(√(Pe² + s²) − Pe)); s·ν bei s = 0,01 für Pe = "+
  ", ".join(f"{pe}: {s/(1-np.exp(-(np.hypot(pe,s)-pe))):.3f}" for pe in [0,0.001,0.01,0.1])+"  → Tiefgrenzfall a → 2Pe·a0 ohne Tully-Fisher.")
P()

# ---------------- 5) Folgen und Prüfungen ----------------
P("5a) Nahfeld: s = r*/r, q = e^{−s} (Sonne, r* = %.0f AE):" % (np.sqrt(G*Msun/a0)/AU))
for nm,r in [("Merkur",0.387*AU),("Erde",AU),("Neptun",30.07*AU),("r*",np.sqrt(G*Msun/a0))]:
    s=np.sqrt(G*Msun/r**2/a0); P(f"      {nm:7s} s = {s:10.2f}   log10 q = {-s/np.log(10):12.2f}")
P("5b) Unabhängigkeit von M: L = r*, D = GM/c ⇒ τ_D = L²/D = c/a0 für jede Masse; q hängt nur von a_N ab.")
for Mx in [1e8,1e10,1e12,1e14]:
    GM=G*Mx*Msun; P(f"      M = 10^{int(np.log10(Mx)):2d} M☉:  τ_D = {GM/a0/(GM/c)/Gyr:.4f} Mrd. J.,  q(a_N = a0) = {np.exp(-np.sqrt(a0/a0)):.6f},  s·ν bei a_N = 10⁻⁴a0: {0.01*nu(0.01):.6f}")
P("5c) Endliche freie Weglänge λ_f bei Lichtgeschwindigkeit (Telegraphengleichung, Wendrate μ = c/λ_f, D = cλ_f/2):")
P("    κL = s·√(1 + ε s²), ε = λ_f/(2ℓ0); q = ½(1 + μ/(μ + Γ + cκ))·e^{−κL}.")
mp.mp.dps=40
def qtel(s,lf,L):
    aN=mp.mpf(s)**2*a0; Gm=aN/c; mu=c/mp.mpf(lf); k=mp.sqrt(Gm*(Gm+2*mu))/c
    return (1+mu/(mu+Gm+c*k))/2*mp.exp(-k*L)
for nm,M in [("10⁸ M☉",1e8*Msun),("6·10¹⁰ M☉",6e10*Msun),("10¹⁵ M☉",1e15*Msun)]:
    GM=G*M; ds=GM/(3*c**2); rs=np.sqrt(GM/a0); vf=(GM*a0)**0.25
    d=[float(qtel(s,6*ds,rs)/mp.exp(-s)-1) for s in [0.3,1,3]]
    P(f"    λ_f = 6Δs, L = r*, {nm:10s}: q/e^{{−s}} − 1 = {d[0]:+.2e} / {d[1]:+.2e} / {d[2]:+.2e} (s = 0,3/1/3);  s·(v_f/c)² bei s = 1: {(vf/c)**2:.2e}, (v_f/c)⁴ = {(vf/c)**4:.2e}")
qb=lambda s: float(qtel(s,2*l0,l0))
P("    L = ℓ0 (ein Schritt, λ_f = 2ℓ0, ballistisch): q = "+" / ".join(f"{qb(s):.4g}" for s in [0.3,1,3])+
  " gegen e^{−s} = "+" / ".join(f"{np.exp(-s):.4f}" for s in [0.3,1,3])+f";  s·ν bei s = 0,01: {0.01/(1-qb(0.01)):.4f} → v⁴ = GMa0/4.")
P("5d) Bahn-Lesart mit L = Schwellenlänge der eingeschlossenen Faltung, L² = G·M(<R)/a0, und Abschirmung 1/R (Kreisbahnmode):")
P("    s_II = √(G·M(<R)/a0)/R gegen lokal s = √(a_N/a0). Dünne Exponentialscheibe, Beispiel r*(M) = 3R_d:")
P("      R/R_d    a_N·R²/(G·M(<R))    s lokal    s_II      a_II/a_lokal − 1")
Rd=1.0; GMd=9.0*a0*Rd**2  # r* = 3 R_d
for x in [0.5,1,2.2,3,5,10,30]:
    y=x/2; vN2=2*GMd/Rd*y*y*(special.i0(y)*special.k0(y)-special.i1(y)*special.k1(y)); R=x*Rd; aN=vN2/R
    Menc=GMd*(1-(1+x)*np.exp(-x)); sI=np.sqrt(aN/a0); sII=np.sqrt(Menc/a0)/R
    P(f"   {x:6.1f}    {aN*R*R/Menc:10.4f}          {sI:8.4f}   {sII:8.4f}    {nu(sII)/nu(sI)-1:+.4f}")
P("    Für Kugelsymmetrie ist a_N R² = G M(<R): beide Lesarten gleich. In Scheiben weicht die Bahn-Lesart im Übergang ab.")
P()
P("Ergebnis: Die Form (Wurzel, Exponential, Σqⁿ) folgt aus einem abgeschirmten Rückweg zweiter Ordnung, der in einem Umlaufsinn")
P("gezählt wird. Die Faltungsgleichung liefert in w die Abschirmung 1/r der Kreisbahnmode; Verlust ∝ a_N an einem Punkt liefert sie nicht.")
P("Offen bleibt die innere Länge L = r* = √(3Δs·ℓ0) (gleichwertig: Γ = c|∇σ| mit τ_D = ℓ0/c); sie steckt im Argument von Q̃.")
txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"rueckweg_faltung.txt"),"w",encoding="utf-8").write(txt); print(txt)
