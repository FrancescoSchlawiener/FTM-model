# Herkunft des kosmischen q (Q = q/(1−q) = 5,36) aus derselben Rückkehr, die a₀ = c²/U festlegt
# Halbgerade der Rückkehr [0, Λ], Λ = U/2 = πc/H, Schritte r₀ = 6Δs mit c; Kac: ⟨T_Ort⟩ = 2Λ/c
# 1 Lesarten: L1 Rückkehr auf der geschlossenen Halbgeraden (Hintergrund, lineares Feld)
#             L2 Kac-Maß in der atmenden Zelle (1 − q = μ(Ort)/μ(X))
#             L3 Kriterium R_Q (lokales q = kosmisches q)       L4 Rand der Hubble-Kugel (a_N = πa₀)
# 2 Zeitabhängigkeit von Q      3 Zufallsprobe (Klassen einfacher Ausdrücke, eigene Lesarten)
# Messwerte mit [G] aus dem Gedächtnis, nicht nachgeprüft.
import numpy as np, sys, io, os, itertools
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.special import zeta
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
np.seterr(all="ignore")

c=2.99792458e8; G=6.674e-11; Msun=1.989e30; pc=3.0857e16; kpc=1e3*pc; Mpc=1e6*pc; Gpc=1e9*pc; yr=3.156e7
h=0.674; H0=67.4e3/Mpc                                              # [G]
MESS={"Planck 2018 TT,TE,EE+lowE+lensing":(0.02237,0.00015,0.1200,0.0012),   # ω_b, σ, ω_c, σ  [G]
      "Planck 2018 + BAO":(0.02242,0.00014,0.11933,0.00091)}
sig8=0.811; ns=0.965                                                # [G]
ob,dob,oc,doc=MESS["Planck 2018 TT,TE,EE+lowE+lensing"]
Qm=oc/ob; dQ=Qm*np.hypot(doc/oc,dob/ob); qm=Qm/(1+Qm); nm=1+Qm
Ob=ob/h**2; Om=(ob+oc)/h**2; Or=9.14e-5; ks=np.pi**2/45; RH=c/H0
P("Messwerte [G] (ohne Korrelation von ω_b und ω_c):")
for k,(b,db,cc,dc) in MESS.items():
    Q_=cc/b; d_=Q_*np.hypot(dc/cc,db/b)
    P(f"   {k:36s} ω_b = {b} ± {db}, ω_c = {cc} ± {dc}:  Q = {Q_:.4f} ± {d_:.4f} ({d_/Q_:.2%}),  q = {Q_/(1+Q_):.5f} ± {d_/(1+Q_)**2:.5f},  ν = 1+Q")
P(f"   Ziel: Q = {Qm:.4f} ± {dQ:.4f} ({dQ/Qm:.2%}), ν = {nm:.4f} ({dQ/nm:.2%}), q = {qm:.5f} ({dQ/(1+Qm)**2/qm:.2%});  s = −ln q = {-np.log(qm):.4f}, s² = {np.log(qm)**2:.4f}")
P()

# Hintergrund des Modells mit atmender Zelle (a₀ = cH/2π, Vakuumteil k_s = π²/45 zu jeder Zeit)
Olr=1-Om-Or-ks
E=lambda z: np.sqrt((Om*(1+z)**3+Or*(1+z)**4+Olr)/(1-ks))
a0z=lambda z: c*H0*E(z)/(2*np.pi)
nuB=lambda s: 1/(1-np.exp(-s))

# ======================= 1 Lesarten =======================
P("1) Lesarten")
P("L1 Rückkehr auf der bei Λ geschlossenen Halbgeraden: D·G'' = Γ·G, G(0) = 1, G'(Λ) = 0  ⇒  G(w) = cosh((Λ−w)κ)/cosh(Λκ),  κ = √(Γ/D) = 1/r")
P("   q = G(L) = cosh((Λ−L)s/L)/cosh(Λs/L),  s = L·κ = √(a_N/a₀);  kleines s: 1 − q = (2Λ/L − 1)·s²/2")
# Prüfung an der diskreten Kette mit Verlust (Schritt h mit c, Überleben je Schritt e^{−Γh/c}, Deckel reflektierend)
def chain_q(ratio,s,N=4000):
    Lam=1.0; L=Lam/ratio; hh=Lam/N; kap=s/L; p=np.exp(-kap**2*hh**2/2)       # Γh/c = κ²·D·h/c = κ²h²/2
    nL=int(round(L/hh)); ab=np.zeros((3,N)); b=np.zeros(N)                  # Unbekannte m_1 … m_N
    ab[1,:]=1.0; ab[0,1:]=-p/2; ab[2,:-1]=-p/2; b[0]=p/2                     # m_i − p/2(m_{i−1}+m_{i+1}) = 0, m_0 = 1
    ab[2,N-2]=-p; ab[1,N-1]=1.0                                               # Deckel: m_N = p·m_{N−1}
    m=solve_banded((1,1),ab,b); return m[nL-1]
qcl=lambda ratio,s: np.cosh((ratio-1)*s)/np.cosh(ratio*s)
for ratio,s in ((2.0,0.5),(np.sqrt(np.pi),np.sqrt(np.pi)),(10.0,0.17),(1e3,0.17)):
    P(f"   Kette N = 4000, Λ/L = {ratio:8.3f}, s = {s:.3f}:  Kette {chain_q(ratio,s):.6f},  Formel {qcl(ratio,s):.6f},  offen e^(−s) = {np.exp(-s):.6f}")
P("   a) Homogener Hintergrund: a_N = 0 ⇒ Γ = 0 ⇒ G ≡ 1, q = 1 für jedes Λ; ν = ∞. Kein endliches Q.")
# b) lineares Feld: g_k = (3/2)Ω_m H₀²(1+z)² D(z) δ₀,k/k  (k mitbewegt), Effektivwert über alle k
Gam=Om*h*np.exp(-Ob-np.sqrt(2*h)*Ob/Om); kk=np.logspace(-5,3,8000)                      # h/Mpc
qT=kk/Gam; Tk=np.log(1+2.34*qT)/(2.34*qT)*(1+3.89*qT+(16.1*qT)**2+(5.46*qT)**3+(6.71*qT)**4)**-0.25
Pk=kk**ns*Tk**2; W=lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
Pk*=(sig8/np.sqrt(np.trapezoid(kk**2*Pk*W(kk*8.0)**2,kk)/(2*np.pi**2)))**2
D2=kk**3*Pk/(2*np.pi**2)
sm1=np.sqrt(np.trapezoid(D2/kk**2,np.log(kk)))*Mpc/h                                    # m
def growth():                                                                           # Start im Strahlungszeitalter (Mészáros: D = 1 + 3y/2)
    f=lambda lna,y: [y[1],-(2+dlnE(lna))*y[1]+1.5*Om*np.exp(-3*lna)/E(np.exp(-lna)-1)**2*y[0]]
    y0=1e-5*Om/Or
    return solve_ivp(f,[np.log(1e-5),0],[1+1.5*y0,1.5*y0],dense_output=True,rtol=1e-9,atol=1e-12)
dlnE=lambda lna: (np.log(E(np.exp(-lna-1e-5)-1))-np.log(E(np.exp(-lna+1e-5)-1)))/2e-5
sol=growth(); Dz=lambda z: sol.sol(-np.log(1+z))[0]/sol.sol(0.0)[0]
P("   b) Lineares Feld (Schätzung: Effektivwert über alle k, BBKS, σ₈, n_s [G]): g = (3/2)Ω_m H₀²(1+z)²D(z)·σ₋₁,  a₀(z) = cH(z)/2π")
L1={}
for z in (0.0,1.0,3.0,1090.0):
    g=1.5*Om*H0**2*(1+z)**2*Dz(z)*sm1; y=g/a0z(z); s=np.sqrt(y); L1[z]=nuB(s)
    P(f"      z = {z:6.0f}:  g = {g:.2e} m/s²,  a_N/a₀ = {y:.2e},  s = {s:.4f},  q = {np.exp(-s):.4f},  ν = {nuB(s):6.2f},  Q = {nuB(s)-1:6.2f}")
P(f"      Deckel: Λk = πck/H ≥ {np.pi*RH/Mpc*0.01*h:.0f} für k ≥ 0,01 h/Mpc, der Faktor aus Λ ändert q um < e^(−180).")
P()

P("L2 Kac in der atmenden Zelle: 1 − q = μ(A)/μ(X), ν = ⟨T_A⟩ in Schritten des Ortes (Kac-Budget: Ω_b/Ω_m + Ω_c/Ω_m = 1)")
P("   Halbgerade [0, Λ], Λ = πc/H; Ort des Kosmos = globale Engstelle r₀ = 6Δs_H = c/H (Δs_H = c/6H):  μ(Ort)/μ(X) = r₀/2Λ = 1/2π")
def kac_return(hs):
    N=len(hs)-1; tau=hs/c; A=np.zeros((N,N)); b=np.zeros(N)
    for i in range(1,N+1):
        k=i-1; A[k,k]=1.0; b[k]=tau[i]
        if i<N:
            if i-1>=1: A[k,i-2]-=0.5
            A[k,i]-=0.5
        else:
            if i-1>=1: A[k,i-2]-=1.0
    return tau[0]+np.linalg.solve(A,b)[0]
rng=np.random.default_rng(3)
for N in (3,10,300):
    hs=np.concatenate([[1.0],rng.uniform(0.1,3,N)]); inner=hs[1:-1].sum()+hs[-1]/2
    hs[1:]*= (np.pi-0.5)/inner                                         # Λ = h₀/2 + Σh_innen + h_N/2 = π (Einheit c/H)
    Lam=0.5*np.sum(hs[:-1]+hs[1:]); T=kac_return(hs*RH)
    P(f"   Kette N = {N:3d}, Schritte zufällig, h₀ = c/H, Λ = {Lam:.6f}·c/H:  ⟨T_Ort⟩·H = {T*H0:.10f}  (2π = {2*np.pi:.10f})")
def cycle_return(hs):                                                   # Schleife: Plätze 0 … N−1, Sprung ±1 je ½, Haltezeit h_j/c
    N=len(hs); tau=hs/c; A=np.eye(N-1); b=tau[1:].copy()
    for i in range(1,N):
        for j in ((i-1)%N,(i+1)%N):
            if j!=0: A[i-1,j-1]-=0.5
    m=np.linalg.solve(A,b); return tau[0]+0.5*(m[0]+m[-1])
for N in (4,50):
    hs=np.concatenate([[1.0],rng.uniform(0.1,3,N-1)]); hs[1:]*=(2*np.pi-1)/hs[1:].sum()   # Umfang Σh = U = 2π·c/H
    P(f"   Schleife N = {N:3d}, Schritte zufällig, h₀ = c/H, Umfang {hs.sum():.6f}·c/H:  ⟨T_Ort⟩·H = {cycle_return(hs*RH)*H0:.10f}")
P(f"   Probe: r*_H² = GM_H/a₀ = (c³/2H)/(cH/2π) = π(c/H)² = r₀·Λ  (Kac-Beziehung r*² = r₀Λ gilt auch für die Hubble-Kugel)")
L2=[("Ort = Engstelle r₀ = c/H (Halbgerade oder Schleife U)",2*np.pi),
    ("Ort = Durchmesser 2c/H",np.pi),
    ("Ort = eigener Schwellenradius r*_H = √π·c/H",2*np.sqrt(np.pi)),
    ("dreidimensional: Hubble-Kugel im Würfel U³",6*np.pi**2),
    ("Zeiteinheit Umlauf U/c statt Schritt c/H",1.0)]
for nmL,nu in L2:
    P(f"   {nmL:54s} ν = {nu:8.4f}   q = {1-1/nu:.5f}   Q = {nu-1:8.4f}   Abstand Q: {(nu-1-Qm)/dQ:+7.1f} σ")
U=2*np.pi*RH
P("   Dieselbe Zuordnung für eine einzelne Masse, ν = 2Λ/r₀ = U/6Δs: " + ", ".join(
  f"{nmM} {U/(2*G*M*Msun/c**2):.1e}" for nmM,M in (("Sonne",1.0),("Milchstraße",6e10),("Haufen",1.7e14))) + " (die Regel dort: ν = 1 … 30)")
P()

P("L3 Kriterium R_Q: lokales q = kosmisches q ⇒ a_N(R_Q) = (ln q)²·a₀ = %.4f a₀, R_Q = r*/(−ln q) = %.2f r*  (Umschreibung von q)" % (np.log(qm)**2,-1/np.log(qm)))
P("   Aus der Zone des Kollapses (Virialdichte 18π²·ρ̄_b(z_c), a₀ = cH₀/2π): s = r*/r_vir ∝ M_b^{1/6}(1+z_c), Q = n_B(s) = 1/(e^s − 1)")
rho_c=3*(100e3/Mpc)**2/(8*np.pi*G)
P("   M_b [M☉]     z_c = 1      z_c = 3      z_c = 6")
for Mb in (1e9,6e10,1e12):
    row=[]
    for zc in (1.0,3.0,6.0):
        rv=(3*Mb*Msun/(4*np.pi*18*np.pi**2*ob*rho_c*(1+zc)**3))**(1/3); rs=np.sqrt(G*Mb*Msun/a0z(0)); s=rs/rv
        row.append(f"Q = {1/np.expm1(s):6.2f}")
    P(f"   {Mb:8.0e}    "+"   ".join(row))
P()

P("L4 Rand der Hubble-Kugel: GM_H = c³/2H, a_N(c/H) = cH/2 = π·a₀, s = √π; geschlossen: L = r*_H, Λ/L = √π, q = cosh(π − s)/cosh π")
L4=[]
for nmL,f in (("volle Dichte",1.0),("nur Ω_m [G]",Om),("nur Ω_b [G]",Ob)):
    s=np.sqrt(np.pi*f); qo=np.exp(-s); qc=np.cosh(np.pi-s)/np.cosh(np.pi)
    L4.append((nmL,1/(1-qc)))
    P(f"   {nmL:14s} s = {s:.4f}:  offen q = {qo:.4f} (Q = {qo/(1-qo):.3f}),  geschlossen q = {qc:.4f} (Q = {qc/(1-qc):.3f})")
P(f"   Verlangt wäre am Rand a_N = {np.log(qm)**2:.4f}·a₀, das {np.pi/np.log(qm)**2:.0f}-mal weniger als cH/2.")
P()

# ======================= 2 Zeitabhängigkeit =======================
P("2) Zeitabhängigkeit (Hintergrund mit atmender Zelle)")
zs=(0.0,1.0,1090.0)
P("   Lesart                                   " + "   ".join(f"Q(z = {z:.0f})" for z in zs))
P("   L2 Engstelle, atmende Zelle U = 2πc/H     " + "   ".join(f"{2*np.pi-1:10.4f}" for z in zs))
U0=27*Gpc
P("   L2 feste Zelle U₀ = 27 Gpc, ν = U₀H/c(1+z) " + "   ".join(f"{U0*H0*E(z)/(c*(1+z))-1:10.3f}" for z in zs))
P("   L1 lineares Feld (Effektivwert)          " + "   ".join(f"{L1[z]-1:10.3f}" for z in zs))
for nmL,f0 in (("L4 nur Ω_m(z)",Om),("L4 nur Ω_b(z)",Ob)):
    vals=[]
    for z in zs:
        fz=f0*(1+z)**3/E(z)**2; s=np.sqrt(np.pi*fz); qc=np.cosh(np.pi-s)/np.cosh(np.pi); vals.append(qc/(1-qc))
    P(f"   {nmL:40s} " + "   ".join(f"{v:10.3f}" for v in vals))
tH=lambda z: __import__("scipy").integrate.quad(lambda zz: 1/((1+zz)*H0*E(zz)),z,np.inf)[0]
P("   Kac-Zeit durch Weltalter, ν = 2π/(H·t)    " + "   ".join(f"{2*np.pi/(H0*E(z)*tH(z))-1:10.3f}" for z in zs))
P("   Nur L2 mit atmender Zelle ist zeitlich fest: Λ = πc/H und r₀ = c/H wachsen gleich, ν = 2π zu jeder Zeit.")
Qt=2*np.pi-1
for k,(b,db,cc,dc) in MESS.items():
    Q_=cc/b; d_=Q_*np.hypot(dc/cc,db/b)
    P(f"   {k:36s} ω_c aus Q = 2π − 1: {Qt*b:.5f} gegen {cc} ± {dc} ({(Qt*b-cc)/dc:+.2f} σ in ω_c); Q: {(Qt-Q_)/d_:+.2f} σ")
P("   Nicht feste Q(z) widersprechen dem exakt staubartigen Hintergrund der Schließungsregel (a³(ρ_b + scheinbare Dichte) fest);")
P(f"   die Strahlung liest Q bei z ≈ 1090 ({Qm:.2f} ± {dQ:.2f} [G]); L1 gäbe dort {L1[1090.0]-1:.0f}, die feste Zelle {U0*H0*E(1090)/(c*1091)-1:.0f}.")
P()

# ======================= 3 Zufallsprobe =======================
P("3) Zufallsprobe: einfache Ausdrücke aus Zahlen des Modells im Band der Messung (δ = |x/Ziel − 1|)")
L2_,L3_,L5_,Lp_=np.log(2),np.log(3),np.log(5),np.log(np.pi)
def classA():
    e=np.arange(-12,13); A,B,C,D=np.meshgrid(e,e,e,e,indexing="ij"); base=(A*L2_+B*L3_+C*L5_+D*Lp_).ravel()
    return np.unique(np.round(np.concatenate([base/n for n in (1,2,3,4)]),12))
def classC():
    p=np.arange(1,65); Pq,Qq,K=np.meshgrid(p,p,np.arange(-6,7),indexing="ij")
    return np.unique(np.round((np.log(Pq)-np.log(Qq)+K*Lp_).ravel(),12))
S={"1":1,"2":2,"3":3,"4":4,"6":6,"9":9,"12":12,"24":24,"π":np.pi,"2π":2*np.pi,"√π":np.sqrt(np.pi),"π²":np.pi**2,
   "e":np.e,"φ":(1+5**.5)/2,"√2":2**.5,"ζ(3)":zeta(3),"π²/6":np.pi**2/6,"π⁴/15":np.pi**4/15,"π²/45":np.pi**2/45,"1/(1−e⁻¹)":1/(1-np.exp(-1))}
def classD():
    v=list(S.values()); vals=set(v)
    for x,y in itertools.product(v,repeat=2):
        for r in (x+y,x-y,x*y,x/y):
            if r>0: vals.add(round(r,12))
    return np.unique(np.round(np.log(np.array(sorted(vals))),12))
CL={"A: 2^a·3^b·5^c·π^d, Wurzeln bis 4":classA(),"C: (p/q)·π^k, p,q ≤ 64":classC(),f"D: x∘y, x,y aus {len(S)} Modellzahlen, ∘ ∈ +−×÷":classD()}
def best(lv,T):
    lt=np.log(T); j=np.searchsorted(lv,lt); cand=lv[max(j-1,0):j+1]; k=np.argmin(abs(cand-lt)); return abs(np.expm1(cand[k]-lt))
cnt=lambda lv,T,d: int(np.sum(abs(lv-np.log(T))<=np.log1p(d)+1e-12))
cntl=lambda lv,T,dl: int(np.sum(abs(lv-np.log(T))<=dl+1e-12))
dens=lambda lv,T,w=0.2: np.sum(abs(lv-np.log(T))<w)/(2*w)
rng=np.random.default_rng(11)
for T,nmT,dr in ((Qm,"Q",dQ/Qm),(nm,"ν = 1+Q",dQ/nm),(qm,"q",dQ/(1+Qm)**2/qm)):
    P(f"   Ziel {nmT} = {T:.5f}, 1σ = {dr:.2%}")
    P("   Klasse                                      Werte     Dichte/ln   bester δ   N(1σ) gef./erw.   N(2σ) gef./erw.")
    for cn,lv in CL.items():
        rho=dens(lv,T)
        P(f"   {cn:42s} {len(lv):9d}   {rho:9.0f}   {best(lv,T):.1e}   {cnt(lv,T,dr):5d} / {rho*2*dr:7.1f}   {cnt(lv,T,2*dr):5d} / {rho*4*dr:7.1f}")
dL2=abs(np.log((2*np.pi-1)/Qm))
P(f"   2π − 1 = {2*np.pi-1:.4f} liegt {dL2:.2%} neben Q; gleich gut oder besser (2π − 1 eingeschlossen): Klasse D {cntl(CL[list(CL)[2]],Qm,dL2)}, C {cntl(CL[list(CL)[1]],Qm,dL2)}, A {cntl(CL[list(CL)[0]],Qm,dL2)}.")
# eigene Lesarten: alle zeitlich festen, endlichen ν
own=[nu for _,nu in L2]+[nu for _,nu in L4]
own=np.array(own); dnu=abs(np.log(2*np.pi/nm)); span=np.log(own.max()/own.min())
p_own=1-(1-min(1,2*dnu/span))**len(own)
tg=np.exp(rng.uniform(np.log(own.min()),np.log(own.max()),200000))
p_mc=np.mean(np.min(abs(np.log(own[None,:]/tg[:,None])),axis=1)<=dnu)
P(f"   Eigene Lesarten mit festem, endlichem ν: {len(own)} Werte zwischen {own.min():.2f} und {own.max():.1f}; bester Abstand zu ν = {nm:.3f}: {dnu:.2%} (2π).")
P(f"   Wahrscheinlichkeit, dass eine von {len(own)} Lesarten zufällig so nahe liegt: {p_own:.3f} (Formel), {p_mc:.3f} (Zufallsziele log-gleich im Bereich).")
P(f"   Eine einzige, zwingende Lesart hätte p = {2*dnu/span:.4f}; für p ≤ 10⁻³ müsste sie auf {1e-3*span/2:.2%} treffen und so genau gemessen sein.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"herkunft_q.txt"),"w",encoding="utf-8").write(txt); print(txt)
