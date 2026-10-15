# Maximum der Schwelle a₀(z) bei z ≈ 0,4: welche Größe des Modells wächst im späten Materiezeitalter und kehrt dort um?
# 0 Verlangte Form (DESI-CPL [G], Vakuumteil trägt die Zeitabhängigkeit)
# 1 Kandidaten X(z), jeweils im selbstkonsistenten H(z) mit a₀(z) = a₀(0)·X(z)/X(0), a₀(0) = cH₀/2π (Θ₀ = 1):
#   Gewicht der Schließungsregel, Anteil jenseits Umkehr und Rate, lineares Wachstum, Gleichstand, Bremsen/Beschleunigen,
#   Rückstellbilanz, Bindungsteil des Strukturterms (lineares Feld und Trichter)
# 2 Tabelle: Lage des Maximums, Form, heutige Steigung, CPL w₀/w_a der gesamten Dunklen Energie
import numpy as np, sys, io, os
from scipy.integrate import solve_ivp, quad, cumulative_trapezoid
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)

# ---------------- Eingaben (wie im Modell) ----------------
c=2.99792458e8; G=6.674e-11; Msun=1.989e30; Mpc=3.0857e22
h=0.674; H0=67.4e3/Mpc; Obh2=0.02237; Q=5.364
Ob=Obh2/h**2; Om=Ob*(1+Q); Or=9.14e-5; ks=np.pi**2/45; ODE=1-Om-Or
a0_0=c*H0/(2*np.pi)
sig8=0.811; ns=0.965; lq=-np.log(Q/(1+Q))            # −ln q_kos
P(f"Eingaben: H₀ = 67,4, Ω_b = {Ob:.4f}, Q = {Q}, Ω_m = {Om:.4f}, Ω_r = {Or:.2e}, k_s = π²/45 = {ks:.4f}, a₀(0) = cH₀/2π = {a0_0:.3e} m/s²")
P()

# ---------------- Gitter ----------------
lna=np.linspace(np.log(1e-3),np.log(3.0),3001); a=np.exp(lna); z=1/a-1
i0=np.argmin(abs(lna)); R=Om*a**-3+Or*a**-4

def dlog(y): return np.gradient(np.log(np.abs(y)),lna)

# ---------------- 0 verlangte Form ----------------
P("0) Verlangte Form: ρ_s(z) = ρ_DE,CPL(z) − (Ω_DE − k_s), a₀(z)/a₀(0) = √(ρ_s/k_s)   (DESI DR2 + CMB + SN, CPL-Werte [G])")
P("   Kombination     z_max   max/heute   a₀(1)/a₀(0)   a₀(2)/a₀(0)   heute/max   z=1/max   s₀ = dln a₀/dln a")
REQ={}
for name,w0,wa in (("Pantheon+",-0.838,-0.62),("DESY5",-0.752,-0.86),("Union3",-0.667,-1.09)):
    rde=ODE*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a)); rs=rde-(ODE-ks)
    m=(z>=0)&(z<=3)&(rs>0); r=np.sqrt(np.where(rs>0,rs,np.nan)/ks)
    j=np.nanargmax(np.where(m,r,np.nan)); f1=interp1d(z[::-1],r[::-1])
    s0=np.gradient(np.log(r),lna)[i0]
    REQ[name]=(z[j],r[j],f1(1.0),f1(2.0),s0)
    P(f"   {name:12s}   {z[j]:5.2f}   {r[j]:8.3f}   {f1(1.0):10.3f}   {f1(2.0):11.3f}   {1/r[j]:8.3f}   {f1(1.0)/r[j]:7.3f}   {s0:+.2f}")
P("   Robust: Maximum bei z = 0,35–0,44; heute/max = 0,80–0,91; a₀(1)/max = 0,87; heutige Steigung −0,76 bis −1,56.")
P()

# ---------------- Wachstum und Strukturgrößen ----------------
def growth(E):
    lnE=np.log(E); dlnE=np.gradient(lnE,lna)
    fE=interp1d(lna,E); fd=interp1d(lna,dlnE)
    def rhs(x,y): return [y[1], -(2+fd(x))*y[1]+1.5*Om*np.exp(-3*x)/fE(x)**2*y[0]]
    x0=np.log(1e-2); sol=solve_ivp(rhs,(x0,lna[-1]),[1e-2,1e-2],t_eval=lna[lna>=x0],rtol=1e-9,atol=1e-12)
    D=np.full_like(a,np.nan); D[lna>=x0]=sol.y[0]; fg=np.full_like(a,np.nan); fg[lna>=x0]=sol.y[1]/sol.y[0]
    D[lna<x0]=a[lna<x0]*D[lna>=x0][0]/1e-2; fg[lna<x0]=1.0
    return D/D[i0], fg

Gam=Om*h*np.exp(-Ob-np.sqrt(2*h)*Ob/Om)
kk=np.logspace(-5,3,4000)                                    # h/Mpc
qT=kk/Gam; T=np.log(1+2.34*qT)/(2.34*qT)*(1+3.89*qT+(16.1*qT)**2+(5.46*qT)**3+(6.71*qT)**4)**-0.25
Pk0=kk**ns*T**2
W=lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
sigR=lambda Rr: np.sqrt(np.trapezoid(kk**2*Pk0*W(kk*Rr)**2,kk)/(2*np.pi**2))
Pk0*=(sig8/sigR(8.0))**2
rhom0=Om*2.775e11                                            # (M☉/h)/(Mpc/h)³
Mg=np.logspace(4,17,300); Rg=(3*Mg/(4*np.pi*rhom0))**(1/3)
sigM=np.array([sigR(r) for r in Rg]); lsig=interp1d(np.log(Mg),np.log(sigM),fill_value="extrapolate")
dlsig=interp1d(np.log(Mg),np.gradient(np.log(sigM),np.log(Mg)),fill_value="extrapolate")
Ig=np.trapezoid(Pk0,kk)/(2*np.pi**2)                        # ∫P dk/2π² in (Mpc/h)²
dta=1.062; dc=1.686

def Qs(s):                                                  # Q_s = 4∫₀^s t³/(eᵗ−1) dt
    s=np.atleast_1d(s).astype(float); outv=np.empty_like(s)
    for i,x in enumerate(s):
        outv[i]=4*quad(lambda t: t**3/np.expm1(t),0,x,limit=200)[0] if x>0 else 0.0
    return outv
sQ=np.logspace(-4,2.5,700); QQ=Qs(sQ); fQs=interp1d(np.log(sQ),QQ,bounds_error=False,fill_value=(0.0,4*np.pi**4/15))
Qs_f=lambda s: np.where(s<1e-4,(4/3)*s**3,fQs(np.log(np.maximum(s,1e-4))))
J=lambda u: quad(lambda x: 4*np.pi*x**2*Qs_f(1/x),1e-6,u,limit=400)[0]
Ju=J(1/lq); vac=4*np.pi**4/15

def MQ(E):                                                   # R_Q = R₂₀₀, a₀ = cH/2π; M in M☉/h
    Hs=H0*E; GM=(1+Q)**3*c**3*lq**6/(8*np.pi**3*1e4*Hs)
    return GM/G/Msun*h

# ---------------- Selbstkonsistenz ----------------
jf=int(np.argmin(abs(z+0.5)))                                  # Zukunft jenseits z = −0,5 eingefroren
def ratio_of(X, p=1.0):
    r=np.sign(X/X[i0])*np.abs(X/X[i0])**p; r=r.copy()
    r[lna>0]=np.maximum(r[lna>0],0.0); r[lna>lna[jf]]=r[jf]; return r
def solve(Xfun, p=1.0, it=50):
    """a₀(z) = a₀(0)·(X/X₀)^p, Θ = (X/X₀)^p/E, ρ_s/ρ_c0 = k_s·(X/X₀)^{2p};  E² = R + Ω_Λr + k_s(X/X₀)^{2p}."""
    Olr=1-ks-Om-Or; E=np.sqrt(R+ODE); res=1.0
    for _ in range(it):
        X=Xfun(E); ratio=np.maximum(ratio_of(X,p),0.0)
        if not np.all(np.isfinite(ratio)): return None, "nicht endlich"
        Enew=np.sqrt(R+Olr+ks*ratio**2)
        res=np.max(abs(Enew/E-1))
        if res<1e-9: E=Enew; break
        E=0.5*E+0.5*Enew
    r=ratio_of(Xfun(E),p)
    bad=(z>=0)&(z<=5)&(r<-1e-6)
    if np.any(bad): return None, f"a₀ < 0 bei z = {z[bad].min():.2f} bis {z[bad].max():.2f}"
    if res>1e-6: return None, f"Fixpunkt nicht erreicht (Rest {res:.0e}, Rückkopplung über Ḣ)"
    return E, ("" if res<1e-9 else f"(Rest {res:.0e})")

def cpl(E, Xfun, p=1.0):
    X=Xfun(E); ratio=(X/X[i0])**p; Olr=1-ks-Om-Or
    rde=Olr+ks*ratio**2; wde=-1-np.gradient(np.log(rde),lna)/3
    m=(z>=0)&(z<=2); A=np.vstack([np.ones(m.sum()),1-a[m]]).T
    w0,wa=np.linalg.lstsq(A,wde[m],rcond=None)[0]
    return w0,wa,wde[i0],ks*ratio**2/E**2

def shape(X):
    m=(z>=-0.5)&(z<=20); Xr=X/X[i0]
    j=np.argmax(np.where(m,Xr,-np.inf)); jmin=np.argmin(np.where(m,Xr,np.inf))
    f=interp1d(z[::-1],Xr[::-1])
    return z[j],Xr[j],f(1.0),f(2.0),dlog(X)[i0],z[jmin],Xr[jmin]

# ---------------- Kandidaten ----------------
def c_w_bg_fix(E): return 1/(1+36*np.pi**2*E**2)              # 1 − w, θ = 3H, θ_c = a₀/c mit a₀ fest
def c_w_rate(E): return E*np.gradient(c_w_bg_fix(E),lna)     # d(1−w)/dt (a₀ fest)
def c_Fta(E):
    D,_=growth(E); M=MQ(E)
    return erfc(dta/(np.sqrt(2)*np.exp(lsig(np.log(M)))*D))
def c_Fta_rate_t(E): return E*np.gradient(c_Fta(E),lna)
def c_Fta_rate_lna(E): return np.gradient(c_Fta(E),lna)
def c_fD(E): D,fg=growth(E); return fg*D                     # dD/d ln a
def c_fDH(E): D,fg=growth(E); return fg*D*E                  # dD/dt
def c_OmODE(E):
    Om_a=Om*a**-3/E**2; return Om_a*(1-Om_a-Or*a**-4/E**2)
def c_OmOs(E): return Om*a**-3/E**2*ks/E**2                  # Vakuumteil fest (ρ_s = konst.)
def c_comH(E): return 1/(a*E)                                 # mitbewegter Hubble-Radius
def c_epsH(E): eps=-np.gradient(np.log(E),lna); return eps/E  # ∝ d(H⁻²)/dt
def c_SdsH(E): return Ob*a**-3/E**2/6                          # ΣΔs_b(<c/H)/(c/H) = Ω_b(z)/6
def c_SdsP(E):
    chi=cumulative_trapezoid(1/(a*E),lna,initial=0)+1/(a[0]*E[0])   # mitbewegter Teilchenhorizont (H₀/c), Strahlungsbeginn
    RP=a*chi; return Ob*a**-3/E**2/6*(E*RP)**2
def c_SdsH_abs(E): return (a*E)**-3                           # ΣΔs_b(<c/H) ∝ Ω_b a⁻³ (c/H)³
def c_Dpot(E): D,fg=growth(E); return (1-fg)*D/a             # −d(D/a)/d ln a, Zerfall des Potentials


P("1) Kandidaten X(z) im selbstkonsistenten H(z) (a₀ ∝ X, Θ₀ = 1). Formwerte X/X(heute); s₀ = dln X/dln a heute.")
P("   Kandidat                                        z_max    max/heute  X(1)/X₀  X(2)/X₀   s₀      w₀ / w_a (ges. DE)  Ω_s(1090)")
cands=[("Gewicht 1 − w, a₀ = cH/2π",None),
       ("Gewicht 1 − w, a₀ fest",c_w_bg_fix),
       ("Rate d(1 − w)/dt, a₀ fest",c_w_rate),
       ("Anteil jenseits Umkehr F(>M_Q)",c_Fta),
       ("  Rate dF/d ln a",c_Fta_rate_lna),
       ("  Rate dF/dt",c_Fta_rate_t),
       ("lineares Wachstum dD/d ln a = fD",c_fD),
       ("lineares Wachstum dD/dt",c_fDH),
       ("Gleichstand Ω_m·Ω_DE (ges.)",c_OmODE),
       ("Gleichstand Ω_m·Ω_s (Vakuumteil fest)",c_OmOs),
       ("mitbewegter Hubble-Radius c/(aH)",c_comH),
       ("  d(H⁻²)/dt ∝ ε/H",c_epsH),
       ("Zerfall des Potentials (1 − f)D/a",c_Dpot),
       ("Rückstellung ΣΔs(<c/H)/(c/H) = Ω_b/6",c_SdsH),
       ("Rückstellung ΣΔs(<R_P)/R_P",c_SdsP),
       ("Rückstellung ΣΔs(<c/H) absolut",c_SdsH_abs)]
TAB={}
for name,fn in cands:
    if fn is None:
        P(f"   {name:46s}  konstant: 1 − w = 1/(1 + 36π²) = {1/(1+36*np.pi**2):.2e} zu jeder Zeit (a₀ ∝ H), kein Extremum"); continue
    E,msg=solve(fn)
    if E is None:
        E0=np.sqrt(R+ODE); X=fn(E0); zm,xm,x1,x2,s0,zmin,xmin=shape(X); w0,wa,_,Os=cpl(E0,fn)
        zm_s=f"{zm:5.2f}" if (zm>-0.49 and zm<19.9) else ("Rand" if zm>=19.9 else "Zukunft")
        P(f"   {name:46s}  {zm_s:>7s}  {xm:8.3f}  {x1:8.3f}  {x2:8.3f}  {s0:+6.2f}   {w0:+.3f} / {wa:+.3f}   (1 Schritt im ΛCDM-Hintergrund; Iteration: {msg})")
        TAB[name]=(zm,xm,x1,x2,s0,w0,wa); continue
    X=fn(E); zm,xm,x1,x2,s0,zmin,xmin=shape(X); w0,wa,wz0,Os=cpl(E,fn)
    TAB[name]=(zm,xm,x1,x2,s0,w0,wa)
    zm_s=f"{zm:5.2f}" if (zm>-0.49 and zm<19.9) else ("Rand" if zm>=19.9 else "Zukunft")
    P(f"   {name:46s}  {zm_s:>7s}  {xm:8.3f}  {x1:8.3f}  {x2:8.3f}  {s0:+6.2f}   {w0:+.3f} / {wa:+.3f}   {Os[0]:.1e} {msg}")
P("   (z_max „Rand“: größter Wert am frühen Ende z = 20; „Zukunft“: Maximum bei z < −0,5.)")
P()

# ---------------- Epochen ----------------
E=np.sqrt(R+ODE)
P("2) Epochen im Hintergrund (ΛCDM-Form mit Q; Θ₀ = 1):")
Om_a=Om*a**-3/E**2; qdec=0.5*Om_a+Or*a**-4/E**2-(1-Om_a-Or*a**-4/E**2)
zacc=brentq(interp1d(z[::-1],qdec[::-1]),0.1,2.0)
zeq=brentq(lambda zz: Om*(1+zz)**3-ODE,-0.5,2)
zeqs=brentq(lambda zz: Om*(1+zz)**3-ks,-0.9,2)
P(f"   Bremsen/Beschleunigen ä = 0 bei z = {zacc:.3f};  ρ_m = ρ_DE bei z = {zeq:.3f};  ρ_m = ρ_s (Vakuumteil π²/45, fest) bei z = {zeqs:.3f}")
P(f"   M_Q(z) (R_Q = R₂₀₀, a₀ = cH/2π): heute {MQ(E)[i0]/h:.2e} M☉, z = 0,4: {interp1d(z[::-1],MQ(E)[::-1]/h)(0.4):.2e}, z = 1: {interp1d(z[::-1],MQ(E)[::-1]/h)(1.0):.2e}")
P()

P("   Abhängigkeit der Lage des Maximums von dF/d ln a von Konventionen (ΛCDM-Hintergrund): Schwellwert δ, Überdichte Δ von R_Δ, Schwelle in M_Q")
D0,_=growth(E)
for dd,dn in ((1.062,"Umkehr"),(1.686,"Kollaps")):
    row=[]
    for Dl in (100,200,500):
        for mode in ("atmend","fest"):
            a0z=c*H0*E/(2*np.pi) if mode=="atmend" else c*H0/(2*np.pi)*np.ones_like(E)
            GM=(1+Q)**3*a0z**3*lq**6*(2/(Dl*(H0*E)**2))**2; M=GM/G/Msun*h
            F=erfc(dd/(np.sqrt(2)*np.exp(lsig(np.log(M)))*D0)); X=np.gradient(F,lna)
            m=(z>=-0.5)&(z<=20); row.append(f"Δ={Dl} {mode}: {z[np.argmax(np.where(m,X,-np.inf))]:5.2f}")
    P(f"   δ = {dd} ({dn}):  " + ";  ".join(row))
P()

# ---------------- Bindungsteil ----------------
P("3) Bindungsteil des Strukturterms, Anteil am Vakuumteil ⟨Q_s⟩/(4π⁴/15); a₀ = cH/2π, Hintergrund ΛCDM-Form; wirkt negativ auf ρ_s")
Dfull,_=growth(E)
def Qlin_at(i):
    a0z=c*H0*E[i]/(2*np.pi); sg=1.5*Om*H0**2*a[i]**-2*Dfull[i]*np.sqrt(Ig)*Mpc/h/(1+Q)   # baryonische Newton-Neigung, 3D
    gq=np.linspace(0,8,801); pdf=np.sqrt(2/np.pi)*gq**2*np.exp(-gq**2/2)
    return np.trapezoid(pdf*Qs_f(np.sqrt(gq*sg/np.sqrt(3)/a0z)),gq)/vac, sg/a0z
def Qhalo_at(i):
    lM=np.log(np.logspace(6,16.5,220)); a0z=c*H0*E[i]/(2*np.pi); sg=np.exp(lsig(lM))*Dfull[i]; nu_=dc/sg
    dndlnM=np.sqrt(2/np.pi)*rhom0/np.exp(lM)*nu_*np.abs(dlsig(lM))*np.exp(-nu_**2/2)      # (Mpc/h)⁻³ mitbewegt
    nphys=dndlnM*h**3/a[i]**3/Mpc**3; rs=np.sqrt(G*np.exp(lM)/h/(1+Q)*Msun/a0z)
    return np.trapezoid(nphys*rs**3*Ju,lM)/vac
zlist=[0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,3.0]
idx=[int(np.argmin(abs(z-zz))) for zz in zlist]
QL=[Qlin_at(i) for i in idx]; QH=[Qhalo_at(i) for i in idx]
P(f"   J(R_Q/r*) = J({1/lq:.2f}) = {Ju:.2f};  σ_g/a₀ heute (baryonisch) = {QL[0][1]:.4f}")
P("   z               " + "  ".join(f"{zz:7.1f}" for zz in zlist))
P("   lineares Feld   " + "  ".join(f"{v[0]:.1e}" for v in QL))
P("   Trichter bis R_Q" + "  ".join(f"{v:.1e}" for v in QH))
zf=np.linspace(0,3,121); qlf=np.array([Qlin_at(int(np.argmin(abs(z-zz))))[0] for zz in zf])
qhf=np.array([Qhalo_at(int(np.argmin(abs(z-zz)))) for zz in zf])
P(f"   Maximum: lineares Feld bei z = {zf[np.argmax(qlf)]:.2f} (×{qlf.max()/qlf[0]:.2f} gegen heute, z = 1: ×{qlf[40]/qlf[0]:.2f}, z = 2: ×{qlf[80]/qlf[0]:.2f});"
  f"  Trichter bei z = {zf[np.argmax(qhf)]:.2f} (×{qhf.max()/qhf[0]:.2f}, z = 1: ×{qhf[40]/qhf[0]:.2f}, z = 2: ×{qhf[80]/qhf[0]:.2f})")
P("   → Form mit innerem Maximum, Betrag 10⁻⁷–10⁻⁴ des Vakuumteils, Vorzeichen negativ: an seinem Maximum ist ρ_s am kleinsten.")
P()

# ---------------- Abstand zur verlangten Form ----------------
P("4) Vergleich mit der verlangten Form (DESY5 [G]: z_max = 0,41, max/heute = 1,18, a₀(1)/max = 0,87, s₀ = −1,16; w₀ > −1, w_a < 0)")
for name,(zm,xm,x1,x2,s0,w0,wa) in TAB.items():
    inner=(zm>0.05 and zm<19.9)
    ok_loc=inner and 0.3<=zm<=0.5
    quad_=(w0>-1 and wa<0) if np.isfinite(w0) else False
    P(f"   {name:46s}  inneres Maximum: {'ja' if inner else 'nein'}{'' if not inner else f' bei {zm:.2f}'};  Lage 0,3–0,5: {'ja' if ok_loc else 'nein'};"
      f"  DESI-Quadrant: {'ja' if quad_ else 'nein'}")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"maximum_schwelle.txt"),"w",encoding="utf-8").write(txt); print(txt)
