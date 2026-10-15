# Dünne Exponentialscheibe mit der Feldform ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N] gegen die algebraische Regel a = ν(a_N/a₀)·a_N
# Achsensymmetrisches Finite-Volumen-Verfahren in (R, z ≥ 0), Einheiten G = M = a₀ = 1 (r* = 1, v_f = 1)
# Prüfung des Lösers an einer Plummer-Kugel, wo Feldform und algebraische Regel exakt gleich sind.
import numpy as np, sys, io, os
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.special import i0, i1, k0, k1
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
nu=lambda y2: 1.0/(-np.expm1(-np.maximum(y2,1e-300)**0.25))   # Argument |g|², s = √(|g|/a₀)

def grid(N,Lmax,h0):
    a=1.0
    for _ in range(200):   # erste Zellbreite ≈ h0
        a=np.arcsinh(Lmax/h0*np.sinh(a)/N*0+Lmax*a/(N*h0)) if False else a
    lo,hi=0.1,20
    for _ in range(100):
        mid=0.5*(lo+hi); w=Lmax*np.sinh(mid/N)/np.sinh(mid)
        if w>h0: lo=mid
        else: hi=mid
    x=Lmax*np.sinh(mid*np.arange(N+1)/N)/np.sinh(mid); return x

def build(Rf,zf):
    NR,NZ=len(Rf)-1,len(zf)-1
    Rc=0.5*(Rf[1:]+Rf[:-1]); zc=0.5*(zf[1:]+zf[:-1])
    return NR,NZ,Rc,zc

def solve(Rf,zf,fdir,bottom_flux,src_vol,flux_weight=None,phiN=None,fN=None,nu_bnd=None):
    """Löst Σ_Flächen c·(Φ_nb − Φ_k)·w = V·S − (bekannte Flüsse). Ohne flux_weight: Laplace-Operator für Φ.
       Mit flux_weight (Feldform): rechte Seite = Σ_Flächen ν_f·c·(Φ_N,nb − Φ_N,k) (Boden hebt sich heraus)."""
    NR,NZ,Rc,zc=build(Rf,zf)
    idx=lambda i,j: i+NR*j
    rows=[];cols=[];vals=[]; rhs=np.zeros(NR*NZ)
    dz=np.diff(zf); dR2=np.diff(Rf**2)
    def add(k,m,c):
        rows.extend([k,k]); cols.extend([k,m]); vals.extend([-c,c])
    for j in range(NZ):
        for i in range(NR):
            k=idx(i,j)
            # radiale Flächen
            if i<NR-1:
                c=2*np.pi*Rf[i+1]*dz[j]/(Rc[i+1]-Rc[i]); add(k,idx(i+1,j),c)
                if flux_weight is not None: rhs[k]+=flux_weight['R'][i,j]*c*(phiN[idx(i+1,j)]-phiN[k])
            else:
                c=2*np.pi*Rf[i+1]*dz[j]/(Rf[i+1]-Rc[i]); rows.append(k); cols.append(k); vals.append(-c)
                rhs[k]-=c*fdir(Rf[i+1],zc[j])
                if flux_weight is not None: rhs[k]+=nu_bnd(Rf[i+1],zc[j])*c*(fN(Rf[i+1],zc[j])-phiN[k])
            if i>0:
                c=2*np.pi*Rf[i]*dz[j]/(Rc[i]-Rc[i-1]); add(k,idx(i-1,j),c)
                if flux_weight is not None: rhs[k]+=flux_weight['R'][i-1,j]*c*(phiN[idx(i-1,j)]-phiN[k])
            # axiale Flächen
            A=np.pi*dR2[i]
            if j<NZ-1:
                c=A/(zc[j+1]-zc[j]); add(k,idx(i,j+1),c)
                if flux_weight is not None: rhs[k]+=flux_weight['z'][i,j]*c*(phiN[idx(i,j+1)]-phiN[k])
            else:
                c=A/(zf[j+1]-zc[j]); rows.append(k); cols.append(k); vals.append(-c)
                rhs[k]-=c*fdir(Rc[i],zf[j+1])
                if flux_weight is not None: rhs[k]+=nu_bnd(Rc[i],zf[j+1])*c*(fN(Rc[i],zf[j+1])-phiN[k])
            if j>0:
                c=A/(zc[j]-zc[j-1]); add(k,idx(i,j-1),c)
                if flux_weight is not None: rhs[k]+=flux_weight['z'][i,j-1]*c*(phiN[idx(i,j-1)]-phiN[k])
            elif flux_weight is None:
                rhs[k]-=bottom_flux[i]
            if flux_weight is None and src_vol is not None:
                rhs[k]+=np.pi*dR2[i]*dz[j]*4*np.pi*src_vol(Rc[i],zc[j])
    L=sp.csc_matrix((vals,(rows,cols)),shape=(NR*NZ,NR*NZ))
    return splu(L).solve(rhs)

def grads(phi,Rf,zf,dphidz0):
    NR,NZ,Rc,zc=build(Rf,zf); F=phi.reshape(NZ,NR).T
    gR=np.zeros_like(F); gz=np.zeros_like(F)
    gR[1:-1,:]=-(F[2:,:]-F[:-2,:])/(Rc[2:]-Rc[:-2])[:,None]; gR[0,:]=-(F[1,:]-F[0,:])/(Rc[1]-Rc[0]); gR[-1,:]=gR[-2,:]
    gz[:,1:-1]=-(F[:,2:]-F[:,:-2])/(zc[2:]-zc[:-2])[None,:]; gz[:,-1]=gz[:,-2]
    gz[:,0]=-0.5*(dphidz0+(F[:,1]-F[:,0])/(zc[1]-zc[0]))
    return gR,gz

def face_nu(phi,Rf,zf,gR,gz):
    NR,NZ,Rc,zc=build(Rf,zf); F=phi.reshape(NZ,NR).T
    nR=np.zeros((NR-1,NZ)); nz=np.zeros((NR,NZ-1))
    gn=-(F[1:,:]-F[:-1,:])/(Rc[1:]-Rc[:-1])[:,None]; gt=0.5*(gz[1:,:]+gz[:-1,:]); nR[:]=nu(gn**2+gt**2)
    gn=-(F[:,1:]-F[:,:-1])/(zc[1:]-zc[:-1])[None,:]; gt=0.5*(gR[:,1:]+gR[:,:-1]); nz[:]=nu(gn**2+gt**2)
    return {'R':nR,'z':nz}

fN=lambda R,z: -1/np.hypot(R,z)
ftot=lambda R,z: -1/np.hypot(R,z)-np.log(-np.expm1(-1/np.hypot(R,z)))
nub=lambda R,z: nu(1/np.hypot(R,z)**4)

N=300; Lmax=80.0
# ---------- Prüfung: Plummer-Kugel b = 0,3 ----------
b=0.3; Rf=grid(N,Lmax,0.004); zf=grid(N,Lmax,0.004); NR,NZ,Rc,zc=build(Rf,zf)
rho=lambda R,z: 3/(4*np.pi*b**3)*(1+(R*R+z*z)/b**2)**-2.5
phiN=solve(Rf,zf,fN,np.zeros(NR),lambda R,z: rho(R,z))
gR,gz=grads(phiN,Rf,zf,np.zeros(NR)); W=face_nu(phiN,Rf,zf,gR,gz)
phi=solve(Rf,zf,ftot,np.zeros(NR),None,flux_weight=W,phiN=phiN,fN=fN,nu_bnd=nub)
GR,GZ=grads(phi,Rf,zf,np.zeros(NR))
P("1) Prüfung des Lösers an der Plummer-Kugel (b = 0,3 r*): Feldform gegen algebraische Regel, Mittelebene")
P("      R/r*    a_N exakt   a_N Gitter   a Feldform   a algebraisch   Verhältnis")
for Rt in (0.1,0.3,1.0,3.0,10.0):
    i=np.argmin(abs(Rc-Rt)); R=Rc[i]; aNex=R/(R*R+b*b)**1.5
    alg=nu(gR[i,0]**2)*gR[i,0]
    P(f"   {R:8.3f}   {aNex:9.5f}   {gR[i,0]:9.5f}   {GR[i,0]:10.5f}   {alg:12.5f}   {GR[i,0]/alg:9.5f}")
P()

# ---------- Prüfung: Kuzmin-Scheibe (Feldform exakt: oberhalb ist ν·g_N wirbelfrei, g = ν(|g_N|)·g_N) ----------
P("2) Prüfung an der Kuzmin-Scheibe Σ = Mh/(2π(R²+h²)^{3/2}): exakt gilt in der Ebene a = ν(|g_N|)·g_N,R mit |g_N| = GM/(R²+h²)")
P("   (das Feld oberhalb ist das einer Punktmasse unter der Ebene, ν·g_N ist dort wirbelfrei)")
for h in (0.3,3.0):
    Rf=grid(N,Lmax,min(0.004,h/60)); zf=grid(N,Lmax,min(0.004,h/60)); NR,NZ,Rc,zc=build(Rf,zf)
    Mann=(h/np.sqrt(Rf[:-1]**2+h*h)-h/np.sqrt(Rf[1:]**2+h*h))
    Sig=h/(2*np.pi*(Rc**2+h*h)**1.5)
    phiN=solve(Rf,zf,fN,-2*np.pi*Mann,None)
    gR,gz=grads(phiN,Rf,zf,2*np.pi*Sig); W=face_nu(phiN,Rf,zf,gR,gz)
    phi=solve(Rf,zf,ftot,None,None,flux_weight=W,phiN=phiN,fN=fN,nu_bnd=nub)
    GR,GZ=grads(phi,Rf,zf,0)
    P(f"   h = {h} r*:     R/h    a exakt (Feldform)   a Gitter    Verhältnis   algebraisch ν(a_N,R)·a_N,R   exakt/alg.")
    for x in (0.25,0.5,1,2,5,10):
        R=x*h
        if R>25: continue
        i=np.argmin(abs(Rc-R)); R=Rc[i]
        gabs=1/(R*R+h*h); gRex=R/(R*R+h*h)**1.5
        ex=nu(gabs**2)*gRex; alg=nu(gRex**2)*gRex
        P(f"              {R/h:6.2f}   {ex:12.5f}      {abs(GR[i,0]):10.5f}   {abs(GR[i,0])/ex:9.5f}   {alg:12.5f}        {ex/alg:8.4f}")
P()

# ---------- dünne Exponentialscheibe ----------
P("3) Dünne Exponentialscheibe: Feldform (volle Lösung) gegen die algebraische Regel in der Mittelebene, a_alg = ν(a_N,R)·a_N,R")
def freeman(R,Rd):
    y=R/(2*Rd); return (2/Rd)*y*y*(i0(y)*k0(y)-i1(y)*k1(y))/R
res={}
for Rd in (0.3,1.0,3.0):
    Rf=grid(N,Lmax,min(0.004,Rd/60)); zf=grid(N,Lmax,min(0.004,Rd/60)); NR,NZ,Rc,zc=build(Rf,zf)
    Mann=((1+Rf[:-1]/Rd)*np.exp(-Rf[:-1]/Rd)-(1+Rf[1:]/Rd)*np.exp(-Rf[1:]/Rd))
    Sig=np.exp(-Rc/Rd)/(2*np.pi*Rd*Rd)
    phiN=solve(Rf,zf,fN,-2*np.pi*Mann,None)          # Bodenfluss ∫∂Φ/∂n dA = −2πG·M_Ring
    gR,gz=grads(phiN,Rf,zf,2*np.pi*Sig); W=face_nu(phiN,Rf,zf,gR,gz)
    phi=solve(Rf,zf,ftot,None,None,flux_weight=W,phiN=phiN,fN=fN,nu_bnd=nub)
    GR,GZ=grads(phi,Rf,zf,2*np.pi*Sig*nu(gR[:,0]**2+(2*np.pi*Sig)**2))
    P(f"   R_d = {Rd} r*:")
    P("      R/R_d    a_N Gitter  a_N Freeman   a Feldform   a algebraisch   Feldform/alg.")
    rows=[]
    for x in (0.25,0.5,1,2,3,5,8,15,40):
        R=x*Rd
        if R>30: continue
        i=np.argmin(abs(Rc-R)); R=Rc[i]
        aF=freeman(R,Rd); alg=nu(gR[i,0]**2)*gR[i,0]
        rows.append((R/Rd,GR[i,0]/alg))
        P(f"   {R/Rd:8.2f}   {gR[i,0]:10.5f}  {aF:10.5f}   {GR[i,0]:10.5f}   {alg:12.5f}   {GR[i,0]/alg:10.4f}")
    res[Rd]=rows
P()
txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"scheibe_feldform.txt"),"w",encoding="utf-8").write(txt); print(txt)
