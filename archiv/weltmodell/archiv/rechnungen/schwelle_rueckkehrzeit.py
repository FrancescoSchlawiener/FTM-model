# Schwelle als Rückkehrzeit im geschlossenen Raum (Kac): τ_D = ⟨T_A⟩ = t·μ(X)/μ(A),  a₀ = c/τ_D = c²/ℓ₀
# 1 Kac auf der Halbgeraden der Hierarchie (Schritte r₀ = 6Δs mit c, Deckel bei Λ): ⟨T_Ort⟩ = 2Λ/c für jede Schrittlänge
# 2 Lesart K1 (A = Ort): Deckel Λ aus der Hubble-Kugel, aus der atmenden Zelle, feste Zelle; Schranke aus dem Vakuumteil
# 3 Lesart K2 (A = gebundene Orte × Ort): a₀ = (c²/2Λ)·F, F = Anteil jenseits der Umkehr; Schar in der Schwellenmasse
# 4 Lesart K3 (X = Volumen, isotrop mit c): Kac für Flüsse ⟨T⟩ = 4V/(c·S_A)
# Jede Lesart selbstkonsistent mit E² = (R + Ω_Λr)/(1 − k_sΘ²), Θ = 2πa₀/(cH), w_s = −1 − (2/3)·dln a₀/dln a
import numpy as np, sys, io, os
from scipy.integrate import solve_ivp
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
np.seterr(all="ignore")

c=2.99792458e8; G=6.674e-11; Msun=1.989e30; pc=3.0857e16; kpc=1e3*pc; Mpc=1e6*pc; Gpc=1e9*pc; yr=3.156e7
h=0.674; H0=67.4e3/Mpc; Obh2=0.02237; Q=5.364
Ob=Obh2/h**2; Om=Ob*(1+Q); Or=9.14e-5; ks=np.pi**2/45; ODE=1-Om-Or
a0m=1.2e-10; RH=c/H0
sig8=0.811; ns=0.965                                     # [G]
P(f"Eingaben: H₀ = 67,4 km/s/Mpc, Ω_m = Ω_b(1+Q) = {Om:.4f}, Ω_r = {Or:.2e}, k_s = π²/45 = {ks:.4f}, c/H₀ = {RH/Gpc:.3f} Gpc; σ₈ = {sig8}, n_s = {ns} [G]")
P()

# ======================= 1 Kac auf der Halbgeraden =======================
P("1) Kette der Ebenen 0 … N (0 = Ort), Schritt h_i mit c (Haltezeit h_i/c), auf und ab gleich häufig, Deckel reflektierend")
P("   Kac: ⟨T_0⟩ = Σ_j deg(j)·(h_j/c) / deg(0) = 2Λ/c,  Λ = Σ_Kanten (h_i + h_{i+1})/2  (exakt aus den mittleren Trefferzeiten gelöst):")
def kac_return(hs):
    N=len(hs)-1; tau=hs/c
    # mittlere Trefferzeit von Ebene i nach 0: m_0 = 0; m_i = τ_i + ½(m_{i−1} + m_{i+1}); m_N = τ_N + m_{N−1}
    A=np.zeros((N,N)); b=np.zeros(N)
    for i in range(1,N+1):
        k=i-1; A[k,k]=1.0; b[k]=tau[i]
        if i<N:
            if i-1>=1: A[k,i-2]-=0.5
            A[k,i]-=0.5
        else:
            if i-1>=1: A[k,i-2]-=1.0
    m=np.linalg.solve(A,b)
    return tau[0]+m[0]                                    # Haltezeit am Ort + Rückweg von Ebene 1
rng=np.random.default_rng(1)
for N in (10,100,1000):
    for kind in ("gleich","zufällig (D springt)"):
        hs=np.ones(N+1) if kind=="gleich" else rng.uniform(0.1,10,N+1)
        T=kac_return(hs); Lam=0.5*np.sum(hs[:-1]+hs[1:])  # Λ = Σ Kanten, Kante = Mittel der Haltelängen

        P(f"   N = {N:5d}, Schritte {kind:22s}: ⟨T_0⟩·c/(2Λ) = {T*c/(2*Lam):.12f}")
P("   ⇒ ⟨T_0⟩ = 2Λ/c hängt nicht von der Schrittlänge ab, also nicht von M (Schritt r₀ = 6Δs ∝ M).")
P("   Die Tiefe, die in dieser Zeit erreicht wird: √(D·⟨T_0⟩) = √((r₀c/2)·(2Λ/c)) = √(r₀·Λ).")
P("   τ_D = ℓ₀/c = ⟨T_0⟩  ⇔  ℓ₀ = 2Λ  ⇔  r*² = 3Δs·ℓ₀ = r₀·Λ  (Schwellenradius = geometrisches Mittel aus Grenzradius und Deckel)")
for nm,M in (("Sonne",1.0),("Milchstraße",6e10)):
    Ds=G*M*Msun/(3*c**2); r0=6*Ds
    for Lam_nm,Lam in (("Λ = c²/2a₀ (a₀ gemessen)",c**2/(2*a0m)),("Λ = πc/H₀",np.pi*RH)):
        P(f"   {nm:12s} r₀ = {r0/(1 if M<10 else 1.496e11):8.4g} {'m ' if M<10 else 'AE'}, {Lam_nm:24s} = {Lam/Gpc:6.2f} Gpc: √(r₀Λ) = {np.sqrt(r0*Lam)/(1.496e11 if M<10 else kpc):8.4g} {'AE' if M<10 else 'kpc'}")
P()

# ======================= Gitter und Hintergrund =======================
lna=np.linspace(np.log(1e-4),np.log(3.0),4001); a=np.exp(lna); z=1/a-1
i0=int(np.argmin(abs(lna))); R=Om*a**-3+Or*a**-4
fz=lambda arr,zz: float(interp1d(z[::-1],arr[::-1])(zz))
dlog=lambda y: np.gradient(np.log(np.abs(y)),lna)

def growth(E):
    """D'' + (2 + dln E/dln a)D' = (3/2)Ω_m a⁻³/E²·D, RK4 auf dem Gitter (Mittelwerte linear), Start D = a bei a = 10⁻³."""
    dlnE=np.gradient(np.log(E),lna); src=1.5*Om*a**-3/E**2
    k0=int(np.searchsorted(lna,np.log(1e-3))); D=np.empty_like(a); D[:k0+1]=a[:k0+1]
    y=np.array([a[k0],a[k0]])
    def f(i,t):
        dd=dlnE[i]+(dlnE[i+1]-dlnE[i])*t; ss=src[i]+(src[i+1]-src[i])*t
        return lambda y: np.array([y[1],-(2+dd)*y[1]+ss*y[0]])
    for i in range(k0,len(a)-1):
        hh=lna[i+1]-lna[i]; f0=f(i,0.0); fm=f(i,0.5); f1=f(i,1.0)
        k1=f0(y); k2=fm(y+0.5*hh*k1); k3=fm(y+0.5*hh*k2); k4=f1(y+hh*k3)
        y=y+hh*(k1+2*k2+2*k3+k4)/6; D[i+1]=y[0]
    return D/D[i0]

# Leistungsspektrum (BBKS, Ω_m = Ω_b(1+Q)), σ(M) in M☉/h
Gam=Om*h*np.exp(-Ob-np.sqrt(2*h)*Ob/Om); kk=np.logspace(-5,4,6000)
qT=kk/Gam; Tk=np.log(1+2.34*qT)/(2.34*qT)*(1+3.89*qT+(16.1*qT)**2+(5.46*qT)**3+(6.71*qT)**4)**-0.25
Pk=kk**ns*Tk**2; W=lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
sigR=lambda Rr: np.sqrt(np.trapezoid(kk**2*Pk*W(kk*Rr)**2,kk)/(2*np.pi**2))
Pk*=(sig8/sigR(8.0))**2; rhom0=Om*2.775e11
sigM=lambda M: sigR((3*M/(4*np.pi*rhom0))**(1/3))

def report(name,E,a0rel,a00):
    """a0rel = a₀(z)/a₀(0) auf dem Gitter, a00 = a₀(0) absolut."""
    Th=2*np.pi*a00*a0rel/(c*H0*E); Os=ks*Th**2
    m=(z>=-0.3)&(z<=10); j=int(np.argmax(np.where(m,a0rel,-np.inf)))
    zraw=z[j]; zmax=zraw if (zraw>-0.29 and zraw<9.9) else np.nan
    zlab=("> 10" if zraw>=9.9 else ("Zuk." if zraw<=-0.29 else f"{zraw:.2f}"))
    if np.all(np.diff(a0rel[m])<0): zlab="–"                                  # monoton fallend mit der Zeit
    sl=dlog(a0rel); ws=-1-(2/3)*sl
    a1=fz(a0rel,1.0); Os1090=ks*(fz(Th,1090.0))**2
    eps0=-dlog(E)[i0]
    return dict(a00=a00,zmax=zmax,zraw=zraw,zlab=zlab,mx=a0rel[j],a1=a1,a1mx=a1/a0rel[j],s0=sl[i0],ws0=ws[i0],ws1=fz(ws,1.0),Os0=Os[i0],Os1090=Os1090,Th=Th,Os=Os,eps0=eps0)

def line(nm,d):
    zs=d["zlab"]
    return (f"   {nm:34s} {d['a00']/1e-10:6.3f}   {d['a1']:6.3f}   {zs:>5s}   {d['s0']:+6.2f}   {d['ws0']:+6.2f} / {d['ws1']:+6.2f}   "
            f"{d['Os0']:6.3f}   {d['Os1090']:.1e}")
HDR="   Lesart                              a₀(0)    a₀(1)/a₀(0)  z_max   dln a₀/dln a   w_s(0) / w_s(1)   Ω_s(0)   Ω_s(1090)"

# ======================= 2 Lesart K1: A = Ort =======================
P("2) Lesart K1: X = Halbgerade der Hierarchie [0, Λ], A = Ort.  a₀ = c²/(2Λ),  Θ = 2πa₀/(cH) = π·(c/H)/Λ")
P(f"   Vakuumteil Ω_s = k_sΘ² < 1 verlangt Λ > π√k_s·c/H = {np.pi*np.sqrt(ks):.3f}·c/H; Ω_s < Ω_DE heute verlangt Λ > {np.pi*np.sqrt(ks/ODE):.3f}·c/H₀.")
P("   K1a Deckel am eigenen Schwellenradius der Hubble-Kugel, Λ = √(GM_H/a₀), GM_H = c³/2H:")
P("       Λ² = GM_H·2Λ/c²  ⇒  Λ = 2GM_H/c² = c/H,  a₀ = cH/2,  Θ = π,  Ω_s = π⁴/45 = %.3f > 1: kein H(z)." % (np.pi**4/45))
P(f"       (Λ = √π·c/H festgehalten, mit a₀ = cH/2π darin: a₀ = cH/(2√π) = {c*H0/(2*np.sqrt(np.pi))/1e-10:.2f}·10⁻¹⁰, Θ = √π, Ω_s = πk_s = {np.pi*ks:.3f} > Ω_DE = {ODE:.3f}.)")
P("       Ebenso jeder Horizont als Deckel, der in de Sitter gegen c/H geht (Hubble, Ereignis): Θ → π.")
res={}
# K1b: atmende Zelle, U = 2πc/H, Deckel beim Gegenpunkt U/2 (gleich: Schleife vom Umfang U, ⟨T⟩ = U/c)
Olr=1-Om-Or-ks; E=np.sqrt((R+Olr)/(1-ks))
res["K1b"]=report("K1b",E,E/E[i0],c*H0/(2*np.pi))
# K1c: feste mitbewegte Zelle U₀: a₀ = c²/(U₀a), ρ_s ∝ a⁻²
for U0nm,U0 in (("27 Gpc",27*Gpc),("c²/a₀,gem",c**2/a0m)):
    th0=2*np.pi*RH/U0; Olr=1-Om-Or-ks*th0**2; E=np.sqrt(R+Olr+ks*th0**2*a**-2)
    res["K1c "+U0nm]=report("K1c",E,1/a,c**2/U0)
P(HDR)
P(line("K1b atmende Zelle, Λ = U/2 = πc/H",res["K1b"]))
P(line("K1c feste Zelle U₀ = 27 Gpc",res["K1c 27 Gpc"]))
P(line("K1c feste Zelle U₀ = c²/a₀ = %.1f Gpc"%(c**2/a0m/Gpc),res["K1c c²/a₀,gem"]))
fN=lambda f: f/(1-f)*(1+3.046*0.2271)/0.2271
P(f"   K1b früh: Ω_s = k_s zu jeder Zeit, im Strahlungszeitalter wie ΔN_eff = {fN(ks):.2f}.  K1c früh: Ω_s ∝ a²·(aH)⁻² → 0.")
P("   Jede Länge, bei der ein geschlossener Raum schließt (Horizont, Zelle, Gegenpunkt), wächst mit der Zeit; a₀ = c²/2Λ fällt dann monoton.")
P()

# ======================= 3 Lesart K2: A = gebundene Orte × Ort =======================
P("3) Lesart K2: X = Orte × Ebenen (Maß: Faltung × Länge in w), A = Ort der gebundenen Bereiche.")
P("   μ(A)/μ(X) = F·(r₀/2Λ)  ⇒  ⟨T_A⟩ = (2Λ/c)/F,  mit Λ = U/2 = πc/H:  a₀ = (cH/2π)·F,  Θ = F ≤ 1.")
P("   F = erfc(δ/(√2·σ(M_min)·D)): Anteil der Faltung in Bereichen jenseits der Umkehr (δ = 1,062) bzw. des Kollapses (1,686),")
P("   lineares Feld wie im Standardmodell mit Ω_m = Ω_b(1+Q); M_min ist die kleinste gebundene Masse (vom Modell nicht festgelegt).")
CACHE={}
def solve_K2(Mmin,dlt,it=80):
    key=(round(np.log10(Mmin),6),dlt)
    if key in CACHE: return CACHE[key]
    E=np.sqrt(R+ODE); sg=sigM(Mmin)
    for _ in range(it):
        D=growth(E); F=erfc(dlt/(np.sqrt(2)*sg*D)); Olr=1-Om-Or-ks*F[i0]**2
        En=np.sqrt((R+Olr)/(1-ks*F**2)); r=np.max(abs(En/E-1)); E=0.6*En+0.4*E
        if r<1e-9: break
    D=growth(E); F=erfc(dlt/(np.sqrt(2)*sg*D))
    CACHE[key]=(E,F); return E,F
P("   δ = 1,062 (Umkehr)")
P("   log M_min [M☉/h]   σ(M_min)   F(0)     a₀(0)    a₀(1)/a₀(0)  z_max   max/heute  a₀(1)/max   dln a₀/dln a   w_s(0)   Ω_s(1090)")
K2={}
for dlt in (1.062,1.686):
    if dlt==1.686:
        P("   δ = 1,686 (Kollaps)")
    for lM in (4,6,8,10,11,12,13,14):
        E,F=solve_K2(10.0**lM,dlt); d=report("K2",E,E*F/(E[i0]*F[i0]),c*H0*F[i0]/(2*np.pi))
        K2[(dlt,lM)]=(d,F[i0])
        zs=d["zlab"]
        P(f"   {lM:8d}           {sigM(10.0**lM):7.3f}   {F[i0]:6.3f}   {d['a00']/1e-10:6.3f}   {d['a1']:9.3f}   {zs:>6s}   {d['mx']:7.3f}   {d['a1mx']:8.3f}   {d['s0']:+9.2f}   {d['ws0']:+7.2f}   {'≈ 0' if d['Os1090']<1e-300 else format(d['Os1090'],'.1e')}")
# M_min mit z_max = 0,40
def zmax_of(lM,dlt):
    E,F=solve_K2(10.0**lM,dlt); d=report("K2",E,E*F/(E[i0]*F[i0]),c*H0*F[i0]/(2*np.pi))
    return d["zraw"], d, F[i0], E, F
P()
P("   Schwellenmasse, bei der das Maximum bei z = 0,35 / 0,40 / 0,44 liegt (verlangte Lage [G]):")
for dlt in (1.062,1.686):
    for zt in (0.35,0.40,0.44):
        lM=brentq(lambda x: zmax_of(x,dlt)[0]-zt,6,14.5,xtol=1e-3)
        zm,d,F0,E,F=zmax_of(lM,dlt)
        Olr=1-Om-Or-ks*F0**2; rde=Olr+ks*F**2*E**2; wde=-1-np.gradient(np.log(rde),lna)/3
        mm=(z>=0)&(z<=2); A_=np.vstack([np.ones(mm.sum()),1-a[mm]]).T; w0,wa=np.linalg.lstsq(A_,wde[mm],rcond=None)[0]
        if dlt==1.062 and zt==0.40: K2EPS=d["eps0"]
        P(f"   δ = {dlt}: z_max = {zm:.2f} bei log M_min = {lM:.2f}:  F(0) = {F0:.3f}, a₀(0) = {d['a00']/1e-10:.3f}·10⁻¹⁰ ({d['a00']/a0m-1:+.0%}),"
          f" max/heute = {d['mx']:.3f}, a₀(1)/max = {d['a1mx']:.3f}, a₀(1)/a₀(0) = {d['a1']:.3f}, a₀(2)/a₀(0) = {fz(E*F/(E[i0]*F[i0]),2.0):.3f}, Steigung {d['s0']:+.2f},"
          f" w_s(0)/w_s(1) = {d['ws0']:+.2f}/{d['ws1']:+.2f}, gesamte DE w₀/w_a = {w0:+.2f}/{wa:+.2f}")
P("   Verlangt [G]: max/heute 1,11–1,26, a₀(1)/max 0,87, Steigung −0,76 bis −1,56, w₀ > −1, w_a < 0.")
P(f"   Wert: a₀(0) = (cH₀/2π)·F(0) ≤ cH₀/2π = {c*H0/(2*np.pi)/1e-10:.3f}·10⁻¹⁰; der Messwert verlangte F(0) = {a0m/(c*H0/(2*np.pi)):.3f} > 1.")
P("   („–“: kein Maximum für z ≤ 10; „> 10“: Maximum früher als z = 10; „Zuk.“: a₀ wächst noch.)")
P()
P("   K2c: gebundene Orte in einer festen mitbewegten Zelle U₀ = 27 Gpc (Λ = U₀a/2): a₀ = c²F/(U₀a), E² = R + Ω_Λr + k_sθ₀²F²(1+z)², θ₀ = 2πc/(U₀H₀)")
U0=27*Gpc; th0=2*np.pi*RH/U0
def solve_K2c(lM,dlt=1.062,it=80):
    E=np.sqrt(R+ODE); sg=sigM(10.0**lM)
    for _ in range(it):
        D=growth(E); F=erfc(dlt/(np.sqrt(2)*sg*D)); Olr=1-Om-Or-ks*th0**2*F[i0]**2
        En=np.sqrt(R+Olr+ks*th0**2*F**2*a**-2); r=np.max(abs(En/E-1)); E=0.6*En+0.4*E
        if r<1e-9: break
    D=growth(E); F=erfc(dlt/(np.sqrt(2)*sg*D))
    return E,F,report("K2c",E,F/(a*F[i0]),c**2*F[i0]/U0)
lMc=brentq(lambda x: solve_K2c(x)[2]["zraw"]-0.40,6,14.5,xtol=1e-3)
E,F,d=solve_K2c(lMc)
P(f"   z_max = 0,40 bei log M_min = {lMc:.2f} (δ = 1,062): F(0) = {F[i0]:.3f}, a₀(0) = {d['a00']/1e-10:.3f}·10⁻¹⁰, max/heute = {d['mx']:.3f}, a₀(1)/max = {d['a1mx']:.3f},"
  f" a₀(1)/a₀(0) = {d['a1']:.3f}, Steigung {d['s0']:+.2f}, w_s(0)/w_s(1) = {d['ws0']:+.2f}/{d['ws1']:+.2f}, a₀(2)/a₀(0) = {fz(F/(a*F[i0]),2.0):.3f}, Ω_s(1090) {'≈ 0' if d['Os1090']<1e-300 else format(d['Os1090'],'.1e')}")
P(f"   Verlangt [G]: a₀(2)/a₀(0) = 0,39–0,51. Der Messwert a₀(0) verlangte U₀ = c²F(0)/a₀ = {c**2*F[i0]/a0m/Gpc:.1f} Gpc < 27 Gpc; heute Steigung −1 + dln F/dln a = {d['s0']:+.2f}.")
P()

P("   Schranke für jede Lesart a₀ = (c²/2Λ)·F mit Λ ∝ c/H und F nicht fallend: dln a₀/dln a = −ε + dln F/dln a ≥ −ε,  ε = −dln H/dln a")
P(f"   ε heute: K1b {res['K1b']['eps0']:.3f}, K2 (z_max = 0,40) {K2EPS:.3f}")
for name,w0,wa,sreq in (("Pantheon+",-0.838,-0.62,-0.76),("DESY5",-0.752,-0.86,-1.16),("Union3",-0.667,-1.09,-1.56)):
    eps=1.5*Om+2*Or+1.5*(1+w0)*ODE
    P(f"   {name:10s} [G]: ε₀ = {eps:.3f}, verlangte Steigung {sreq:+.2f}  ⇒  dln(μ(A)/μ(X))/dln a heute = {sreq+eps:+.2f} (Menge A müsste schrumpfen)")
P()

# ======================= 4 Lesart K3: X = Volumen =======================
P("4) Lesart K3: X = Volumen V, Bewegung isotrop mit c; Kac für Flüsse (mittlere Sehne): ⟨T_A⟩ = 4V/(c·S_A)")
V_H=4*np.pi/3*RH**3
P("   a) A = Kugel vom Radius r* um das System, V = Hubble-Volumen: ⟨T⟩ = V/(π·r*²·c) ∝ 1/M")
for nm,M in (("Sonne",1.0),("Milchstraße",6e10),("Haufen",1.7e14)):
    rs=np.sqrt(G*M*Msun/a0m); T=V_H/(np.pi*rs**2*c)
    P(f"      {nm:12s} r* = {rs/kpc:9.3g} kpc   ⟨T⟩ = {T/yr:9.2e} Jahre   (τ_D = ℓ₀/c = {c/a0m/yr:.2e} Jahre)")
P("      a₀ = c/⟨T⟩ ∝ r*² ∝ M/a₀  ⇒  a₀ ∝ √M,  v⁴ = GMa₀ ∝ M^{3/2}: Tully-Fisher mit falscher Steigung.")
P("   b) A = alle gebundenen Bereiche (Radius bei der Umkehr, Dichte 5,55·ρ̄), ⟨T⟩ = mittlere freie Weglänge/c, ℓ₀ = ℓ_frei:")
lnMg=np.linspace(np.log(1e6),np.log(1e16),400); Mg=np.exp(lnMg)
sg=np.array([sigM(M) for M in Mg]); dls=np.gradient(np.log(sg),lnMg)
E0=np.sqrt(R+ODE); D0=growth(E0)
for zz in (0.0,0.4,1.0):
    Dz=fz(D0,zz); nu_=1.062/(sg*Dz)
    dn=np.sqrt(2/np.pi)*rhom0/Mg*nu_*np.abs(dls)*np.exp(-nu_**2/2)       # (h/Mpc)³ mitbewegt je ln M
    Rta=(3*Mg/(4*np.pi*5.55*rhom0))**(1/3)                                # mitbewegt, Mpc/h
    row=[]
    for lMm in (8,11,14):
        mm=Mg>=10.0**lMm
        sig_com=np.trapezoid((dn*np.pi*Rta**2)[mm],lnMg[mm])               # h/Mpc (mitbewegt)
        l_phys=1/sig_com/h*Mpc/(1+zz)                                      # physikalisch: n ∝ (1+z)³, σ ∝ (1+z)⁻²
        row.append(f"M_min = 10^{lMm}: ℓ = {l_phys/Gpc:7.3f} Gpc, a₀ = c²/ℓ = {c**2/l_phys:.1e}")
    P(f"      z = {zz:3.1f}:  "+";  ".join(row))
P(f"      verlangt ℓ₀ = c²/a₀ = {c**2/a0m/Gpc:.1f} Gpc: die freie Weglänge zwischen gebundenen Bereichen ist 10–10⁴-mal kürzer.")
P()

P("5) Zusammenfassung")
P("   K1 (Ort): τ_D = 2Λ/c für jede Masse, ℓ₀ = 2Λ. Mit Λ = U/2 folgt a₀ = cH/2π (Wert des Modells), Zeitgesetz a₀ ∝ H: monoton, früh ΔN_eff = 2,09.")
P("             Deckel an einem Horizont (→ c/H) gibt Θ = π, Ω_s > 1. Feste Zelle: a₀ ∝ 1 + z; der Messwert verlangte U₀ = 24,3 Gpc < 27 Gpc.")
P("   K2 (gebunden): Maximum im Fenster nur mit gewählter Schwellenmasse; dann F(0) ≈ 0,5, a₀(0) ≈ 0,52·10⁻¹⁰, Maximum flach (1 %).")
P("   K3 (Volumen): a₀ hängt an M (Kugel r*) oder ℓ₀ ist um Größenordnungen zu kurz (gebundene Bereiche).")
P("   K2c (feste Zelle): Breite und z = 2 nahe, Höhe 1,09, Steigung −0,47; Wert verlangte U₀ = 10,8 Gpc < 27 Gpc; U₀ und M_min frei, die Zelle atmet nicht.")
P("   Allgemein: heute verlangt ist ein a₀, das schneller fällt als H; in a₀ ∝ H·μ(A)/μ(X) (atmende Zelle) hieße das eine schrumpfende Menge A.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"schwelle_rueckkehrzeit.txt"),"w",encoding="utf-8").write(txt); print(txt)
