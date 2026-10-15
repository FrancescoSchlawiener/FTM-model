# Warp der Milchstraße: welches Drehmoment liefert eine rechtläufige Knotendrehung von +10,9 km/s/kpc [G]?
# 1 Baryonisches Massenmodell (McMillan 2017 [G], Annahme), kugelgemitteltes Feld, Feldform ν = 1/(1 − e^{−√y})
# 2 Lineare Antwort der Feldgleichung auf einen Punktstörer (alle geraden l), Drehmoment auf einen gekippten Ring
# 3 Sagittarius-Zwerg, Große Magellansche Wolke (Newton, ×ν am Ort der Scheibe, volle Antwort), M31
# 4 Äußeres Feld der Lokalen Gruppe (Anteil l = 2 der Feldgleichung)
# 5 Schräg einfallendes Gas (Drehimpulszufuhr)   6 Verwindung a_J (Lense-Thirring, verstärkt oder nicht)   7 Bilanz
# Einheiten: kpc, km/s, M☉;  1 km/s/kpc = 1,0227 Gyr⁻¹.  Vorzeichen: > 0 rechtläufig (im Drehsinn der Scheibe, L ∥ +z).
import numpy as np, sys, io, os
np.seterr(over="ignore", divide="ignore", invalid="ignore")
from scipy.linalg import solve_banded
from scipy.special import eval_legendre
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
G=4.30091e-6; C_KMS=299792.458; UACC=1e6/3.0857e19; a0=float(__import__("os").environ.get("FTM_A0", "1.2e-10"))/UACC; GYR=1.0227
nu_of=lambda y: 1/(-np.expm1(-np.sqrt(np.maximum(y,1e-30))))
def dnu_dy(y):
    y=np.maximum(y,1e-30); s=np.sqrt(y); q=np.exp(-s); return -q/(2*s*(1-q)**2)
OMEGA_GAIA=10.9      # km/s/kpc [G]
NCHK=[0,0]
def chk(name,x,y,tol):
    rel=abs(x-y)/max(abs(y),1e-300); NCHK[0]+=1; NCHK[1]+=(rel>=tol)
    P(f"   [{'OK' if rel<tol else 'ABW'}] {name}: {x:.6g} gegen {y:.6g} (rel. {rel:.1e})")

# ---------------- 1 Massenmodell ----------------
P("1) Baryonische Milchstraße (Massenmodell McMillan 2017 [G], Annahme), Feld kugelgemittelt, a₀ = 1,2·10⁻¹⁰ m/s²")
Rg=np.concatenate([[0],np.geomspace(1e-3,400,2500)]); zg=np.concatenate([[0],np.geomspace(1e-4,400,2500)])
Rc=0.5*(Rg[1:]+Rg[:-1]); zc=0.5*(zg[1:]+zg[:-1]); dR=np.diff(Rg); dz=np.diff(zg)
RR,ZZ=np.meshgrid(Rc,zc,indexing="ij"); VOL=2*np.pi*RR*dR[:,None]*dz[None,:]*2      # beide Hälften
def rho_comp(R,z):
    rp=np.sqrt(R**2+(z/0.5)**2); out_={}
    out_["bulge"]=98.4e9/(1+rp/0.075)**1.8*np.exp(-(rp/2.1)**2)
    for k,(S0,Rd,zd) in (("thin",(896e6,2.5,0.3)),("thick",(183e6,3.02,0.9))):
        out_[k]=S0/(2*zd)*np.exp(-R/Rd-np.abs(z)/zd)
    for k,(S0,Rd,Rm,zd) in (("HI",(53.1e6,7.0,4.0,0.085)),("H2",(2180e6,1.5,12.0,0.045))):
        Rs=np.maximum(R,1e-6); e=np.exp(-np.abs(z)/zd)
        out_[k]=S0/(4*zd)*np.exp(-Rm/Rs-Rs/Rd)*4*e/(1+e)**2
    return out_
comp=rho_comp(RR,ZZ); Mc={k:(v*VOL).sum() for k,v in comp.items()}
comp["bulge"]*=9.23e9/Mc["bulge"]; Mc["bulge"]=9.23e9
rho=sum(comp.values()); disc=comp["thin"]+comp["thick"]+comp["HI"]+comp["H2"]
Mb=sum(Mc.values())
P("   " + ", ".join(f"{k} {v:.3e}" for k,v in Mc.items()) + f";  M_b = {Mb:.3e} M☉")
chk("dünne Scheibe gegen 2πΣ₀R_d²",Mc["thin"],2*np.pi*896e6*2.5**2,5e-3)
rr_=np.sqrt(RR**2+ZZ**2).ravel(); o=np.argsort(rr_); rs_=rr_[o]; Mcum=np.cumsum((rho*VOL).ravel()[o])
M_enc=lambda r: np.interp(r,rs_,Mcum)
gN_mw=lambda r: G*M_enc(r)/r**2
def vcirc(R, form="FTM"):
    g=gN_mw(R); return np.sqrt(R*g*(nu_of(g/a0) if form=="FTM" else 1.0))
R_w=14.0; Om_F=vcirc(R_w)/R_w; Om_N=vcirc(R_w,"N")/R_w
P(f"   bei R = {R_w} kpc: g_N/a₀ = {gN_mw(R_w)/a0:.3f}, ν = {nu_of(gN_mw(R_w)/a0):.3f}, v = {vcirc(R_w):.1f} (Newton {vcirc(R_w,'N'):.1f}) km/s, "
  f"Ω = {Om_F:.2f} (Newton {Om_N:.2f}) km/s/kpc;  r* = {np.sqrt(G*Mb/a0):.2f} kpc")
P()

# ---------------- 2 lineare Antwort ----------------
rgr=np.logspace(-1.5,5,8001); lr=np.log(rgr); rf=np.exp(0.5*(lr[1:]+lr[:-1]))
rfe=np.concatenate([[rgr[0]**2/rf[0]],rf,[rgr[-1]**2/rf[-1]]]); drc=rfe[1:]-rfe[:-1]
def yN_from_g(gobs):                                   # g_obs/a₀ → g_N/a₀
    lo,hi=1e-14,1e3
    for _ in range(200):
        mid=np.sqrt(lo*hi)
        if nu_of(mid)*mid<gobs: lo=mid
        else: hi=mid
    return mid
def response(Mp, dp, l, ge_obs=0.02, alpha_const=None):
    """(r²u')' − l(l+1)u = (r²αf')' − l(l+1)(ν−1)f,  α = ν − 1 + yν'(y) radial, ν − 1 tangential; Hintergrund √(g_N² + g_Ne²)."""
    yNe=yN_from_g(ge_obs) if ge_obs>0 else 0.0
    f=lambda x: np.where(x<dp,-G*Mp*x**l/dp**(l+1),-G*Mp*dp**l/x**(l+1))
    fp=lambda x: np.where(x<dp,-l*G*Mp*x**(l-1)/dp**(l+1),(l+1)*G*Mp*dp**l/x**(l+2))
    def bg(x):
        yy=np.sqrt(gN_mw(x)**2+(yNe*a0)**2)/a0; return nu_of(yy)-1, nu_of(yy)-1+yy*dnu_dy(yy)
    n1c,_=bg(rgr); n1f,alf=bg(rf)
    if alpha_const is not None: n1c=np.full_like(rgr,alpha_const); alf=np.full_like(rf,alpha_const)
    wf=rf**2/(rgr[1:]-rgr[:-1]); n=len(rgr)
    main=-l*(l+1)*drc.copy(); up=np.zeros(n); lo=np.zeros(n)
    main[:-1]-=wf; main[1:]-=wf; up[1:]=wf; lo[:-1]=wf
    flux=rf**2*alf*fp(rf); rhs=-l*(l+1)*n1c*f(rgr)*drc; rhs[:-1]+=flux; rhs[1:]-=flux
    main[0]=1; up[1]=0; rhs[0]=0; main[-1]=1; lo[-2]=0; rhs[-1]=0
    ab=np.zeros((3,n)); ab[0]=up; ab[1]=main; ab[2]=lo
    u=solve_banded((1,1),ab,rhs); return u, f(rgr)
u_t,f_t=response(3.2e9,50.0,2,alpha_const=1.5); k14=np.argmin(abs(rgr-14))
P("2) Lineare Antwort der Feldform auf einen Punktstörer, Multipole l (Kugelhintergrund der Milchstraße + äußeres Feld)")
chk("konstantes ν − 1 = 1,5: u/f bei 14 kpc",u_t[k14]/f_t[k14],1.5,3e-3)

def multipoles(Mp, dp, lmax, mode, ge_obs=0.02):
    """Koeffizienten u_l(R_w), u_l'(R_w) des Störpotentials Σ u_l(r)P_l(r̂·ê); mode: 'N', 'nu' (×ν am Ort der Scheibe), 'F' (volle Antwort)."""
    res=[]
    for l in range(2,lmax+1,2):
        f=np.where(rgr<dp,-G*Mp*rgr**l/dp**(l+1),-G*Mp*dp**l/rgr**(l+1))
        if mode=="F": u,_=response(Mp,dp,l,ge_obs)
        elif mode=="nu": u=(nu_of(gN_mw(R_w)/a0)-1)*f
        else: u=0*f
        tot=f+u; val=np.interp(R_w,rgr,tot); der=np.interp(R_w,rgr,np.gradient(tot,rgr))
        res.append((l,val,der))
    return res

def node_rate(mps, beta_deg, theta, vc, nphi=73, nring=720):
    """Knotenrate des um θ gekippten Rings (Knoten längs x) für alle Azimute φ_p des Störers; Rückgabe max rechtläufig, max rückläufig."""
    ph=np.linspace(0,2*np.pi,nring,endpoint=False)
    pts=R_w*np.stack([np.cos(ph),np.sin(ph)*np.cos(theta),np.sin(ph)*np.sin(theta)],1); rh=pts/R_w
    rates=[]
    for php in np.linspace(0,2*np.pi,nphi,endpoint=False):
        b=np.radians(beta_deg); e=np.array([np.cos(b)*np.cos(php),np.cos(b)*np.sin(php),np.sin(b)])
        cg=rh@e; grad=np.zeros_like(pts)
        for l,val,der in mps:
            Pl=eval_legendre(l,cg); dPl=l*(eval_legendre(l-1,cg)-cg*Pl)/np.maximum(1-cg**2,1e-12)
            grad+=der*Pl[:,None]*rh+val*dPl[:,None]*(e[None,:]-cg[:,None]*rh)/R_w
        tau=np.cross(pts,-grad).mean(0)
        rates.append(tau[0]/(R_w*vc*np.sin(theta)))
    rates=np.array(rates); return rates.max(), rates.min()

def node_rate_newton_exact(Mp, dp, beta_deg, theta, vc, nphi=73, nring=720):
    ph=np.linspace(0,2*np.pi,nring,endpoint=False)
    pts=R_w*np.stack([np.cos(ph),np.sin(ph)*np.cos(theta),np.sin(ph)*np.sin(theta)],1); rates=[]
    for php in np.linspace(0,2*np.pi,nphi,endpoint=False):
        b=np.radians(beta_deg); e=np.array([np.cos(b)*np.cos(php),np.cos(b)*np.sin(php),np.sin(b)])*dp
        dvec=pts-e[None,:]; grad=G*Mp*dvec/np.linalg.norm(dvec,axis=1)[:,None]**3
        tau=np.cross(pts,-grad).mean(0); rates.append(tau[0]/(R_w*vc*np.sin(theta)))
    rates=np.array(rates); return rates.max(), rates.min()
P()

# ---------------- 3 Störer ----------------
P("3) Punktstörer: größte Knotenrate eines Rings bei R = 14 kpc, Neigung θ_w, über alle Azimute des Störers (km/s/kpc)")
P("   Annahmen [G]: LMC 3,2·10⁹ M☉ baryonisch, 50 kpc, Breite −33°;  Sgr 1·10⁸–5·10⁸ M☉ baryonisch, 18 kpc, Breite −20°;")
P("                 M31 1·10¹¹ M☉ baryonisch, 780 kpc, Breite −22°;  äußeres Feld g_e = 0,02 a₀;  θ_w = 0,1 (Warp-Höhe ≈ 1,4 kpc bei 14 kpc)")
theta=0.1
PERT=[("LMC",3.2e9,50.0,-33.0,24),("Sgr (5·10⁸)",5e8,18.0,-20.0,60),("Sgr (1·10⁸)",1e8,18.0,-20.0,60),("M31",1e11,780.0,-22.0,2)]
RES={}
P("   Störer          l_max  Newton exakt       Newton (Multipole)   ×ν am Ort der Scheibe   volle Antwort (Feldform)   Verstärkung")
for name,Mp,dp,beta,lmax in PERT:
    ne=node_rate_newton_exact(Mp,dp,beta,theta,vcirc(R_w,"N"))
    mN=node_rate(multipoles(Mp,dp,lmax,"N"),beta,theta,vcirc(R_w,"N"))
    mnu=node_rate(multipoles(Mp,dp,lmax,"nu"),beta,theta,vcirc(R_w))
    mF=node_rate(multipoles(Mp,dp,lmax,"F"),beta,theta,vcirc(R_w))
    RES[name]=(ne,mN,mnu,mF)
    amp=max(abs(mF[0]),abs(mF[1]))/max(abs(mN[0]),abs(mN[1]))*vcirc(R_w)/vcirc(R_w,"N")
    P(f"   {name:14s} {lmax:4d}   +{ne[0]:.2e}/{ne[1]:.2e}   +{mN[0]:.2e}/{mN[1]:.2e}   +{mnu[0]:.2e}/{mnu[1]:.2e}     +{mF[0]:.2e}/{mF[1]:.2e}      ×{amp:.1f} (Feld)")
P("   (+max/−min: rechtläufig/rückläufig je nach Azimut des Störers gegen die Knotenlinie; das Vorzeichen setzt die Lage.)")
ell=lambda M,d: d*(M/(3*Mb))**(1/3)
P(f"   Glättungslängen im Feld der Milchstraße ℓ = d(M/3M_b)^{{1/3}}: LMC {ell(3.2e9,50):.1f} kpc, Sgr {ell(5e8,18):.1f} kpc, M31 {ell(1e11,780):.0f} kpc;"
  f" der Ring liegt außerhalb der Zonen von LMC und Sgr (kollektives Feld)")
lmcF=RES["LMC"][3][0]; Mreq=3.2e9*OMEGA_GAIA/lmcF
P(f"   Für +10,9 km/s/kpc bräuchte die LMC (volle Antwort, θ_w = 0,1) baryonisch {Mreq:.2e} M☉ = {Mreq/Mb:.1f}·M_b der Milchstraße;"
  f" Sgr {5e8*OMEGA_GAIA/RES['Sgr (5·10⁸)'][3][0]:.1e} M☉")
P()

# ---------------- 4 äußeres Feld ----------------
P("4) Äußeres Feld der Lokalen Gruppe: homogenes g_e übt in Newton kein Drehmoment aus; in der Feldform Anteil l = 2 ∝ g_e²")
rg=np.logspace(-1,5,6001); mu_n,mu_w=np.polynomial.legendre.leggauss(400); P2=0.5*(3*mu_n**2-1)
Mcs=M_enc(rg); rho_sph=np.gradient(Mcs,rg)/(4*np.pi*rg**2)
def efe_C(ge):
    gNe=yN_from_g(ge)*a0; r=rg[:,None]; mu=mu_n[None,:]; gr=(G*Mcs/rg**2)[:,None]; dgr=np.gradient(G*Mcs/rg**2,rg)[:,None]
    s2=gr**2+gNe**2-2*gr*gNe*mu; ab=np.sqrt(s2)
    ds2dr=2*gr*dgr-2*dgr*gNe*mu; ds2dmu=-2*gr*gNe
    gdot=((-gr+gNe*mu)*ds2dr+(gNe*(1-mu**2)/r)*ds2dmu)/(2*ab); y=ab/a0
    div=(nu_of(y)-1)*(-4*np.pi*G*rho_sph[:,None])+dnu_dy(y)/a0*gdot
    rho2=2.5*((-div/(4*np.pi*G))*P2[None,:]*mu_w[None,:]).sum(1)
    m=rg<=R_w
    return -(4*np.pi*G/5)*(np.trapezoid(rho2[m]*rg[m]**4,rg[m])/R_w**5+np.trapezoid(rho2[~m]/rg[~m],rg[~m]))
for ge in (0.01,0.02,0.03):
    Cq=efe_C(ge); wmax=1.5*abs(Cq)*0.5/(Om_F*theta)
    P(f"   g_e = {ge} a₀: C(14) = {Cq:+.2e} (km/s/kpc)²,  größte Knotenrate (θ_w = 0,1, cos β sin β = ½) {wmax:.1e} km/s/kpc")
P()

# ---------------- 5 Gas ----------------
P("5) Schräg einfallendes Gas: dL/dt = Ṁ·j_⊥. Knotenrate ω = (Ṁ_o/M_o)·(j_⊥/j_c)/sin θ_w (j_c = R·v_c, Gas mit Kreisbahn-Drehimpuls, Winkel ψ)")
Mo=(disc*VOL)[Rc>12].sum()
P(f"   Außenscheibe R > 12 kpc: M_o = {Mo:.2e} M☉ (Annahme: alles einfallende Gas landet dort)")
P("   Ṁ [M☉/Jahr, Annahme]   θ_w     ω (ψ = 90°) [km/s/kpc]   ω (ψ = 30°)")
for Md in (1.0,3.0):
    for th in (0.05,0.1):
        w=Md*1e9/Mo/np.sin(th)/GYR
        P(f"   {Md:6.1f}                  {th:4.2f}    {w:8.2f}                 {w*0.5:8.2f}")
for th in (0.05,0.1):
    P(f"   Für +10,9 km/s/kpc bei θ_w = {th}: Ṁ·sin ψ = {OMEGA_GAIA*GYR*Mo*np.sin(th)/1e9:.1f} M☉/Jahr")
P(f"   Im Modell ist j_c um v_FTM/v_N = {vcirc(R_w)/vcirc(R_w,'N'):.2f} größer; stammt j_⊥ aus Kreisbahnen, kürzt sich das; ist j_⊥ fest, sinkt ω um diesen Faktor.")
P("   Das Vorzeichen setzt die Richtung von j_⊥ gegen die Knotenlinie; die Drehung endet nach höchstens 90°, wenn die Scheibe zu j_⊥ gekippt ist.")
P()

# ---------------- 6 Verwindung ----------------
P("6) Verwindung der eigenen Rotation (gravitomagnetisch, rechtläufig)")
vc_grid=vcirc(np.maximum(RR,1e-3))
J=(disc*RR*vc_grid*VOL).sum()+0.0
aJ=J/(Mb*C_KMS); ds=G*Mb/(3*C_KMS**2)
wLT=2*G*J/(C_KMS**2*R_w**3); nuR=nu_of(gN_mw(R_w)/a0)
P(f"   J = {J:.3e} M☉ kpc km/s,  a_J = J/(M_b c) = {aJ*1e3:.2f} pc,  Δs = {ds*3.0857e19:.3e} m,  Spinzahl j = a_J/(3Δs) = {aJ/(3*ds):.0f}")
P(f"   ω_LT = 2GJ/(c²R³) bei 14 kpc: {wLT:.2e};  ×ν(R) = {wLT*nuR:.2e};  ×ν² (Strom und Feld verstärkt) = {wLT*nuR**2:.2e} km/s/kpc")
P(f"   Abstand zum Messwert: Faktor {OMEGA_GAIA/(wLT*nuR**2):.1e}")
P()

# ---------------- 7 Bilanz ----------------
P("7) Bilanz bei R = 14 kpc (größte Beträge, θ_w = 0,1): Messwert +10,9 km/s/kpc [G]")
P(f"   Sgr (5·10⁸ M☉):  {max(abs(np.array(RES['Sgr (5·10⁸)'][3]))):.2e}   Vorzeichen: Lage")
P(f"   LMC:             {max(abs(np.array(RES['LMC'][3]))):.2e}   Vorzeichen: Lage")
P(f"   M31:             {max(abs(np.array(RES['M31'][3]))):.2e}   Vorzeichen: Lage (lineare Antwort, obere Schranke)")
P(f"   äußeres Feld:    {1.5*abs(efe_C(0.03))*0.5/(Om_F*theta):.2e}   (g_e = 0,03 a₀)")
P(f"   Gas, 1–3 M☉/Jahr, ψ = 90°: {1e9/Mo/np.sin(theta)/GYR:.2f}–{3e9/Mo/np.sin(theta)/GYR:.2f}   Vorzeichen: Richtung des Einfalls")
P(f"   Verwindung:      {wLT:.1e}–{wLT*nuR**2:.1e}   rechtläufig")
P(f"Prüfungen: {NCHK[0]}, Abweichungen über Toleranz: {NCHK[1]}")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"warp_drehmoment.txt"),"w",encoding="utf-8").write(txt); print(txt)
