# Halbgerade der Hierarchie: q = G(L) mit D·G'' = Γ·G, G(0) = 1, G beschränkt. Ist die innere Richtung w die Richtung der Hierarchie der Faltungen?
# 1 Irrfahrt auf den Ebenen (Kette, Boden, Richtung, Drift, Glättungslänge, L, Zahl der Ebenen)
# 2 Schließung der Hierarchie (Λ aus Hubble-Radius und Zelle, 10-%-Radius Milchstraße und Zwerg, Spitze der Hierarchie)
# 3 Umgebung (äußeres Feld im Verlust Γ gegen Ebenengrenze mit additivem D)   4 prüfbare Folge
import numpy as np, sys, io, os
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import solve_banded
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30
pc=3.0857e16; kpc=1e3*pc; Mpc=1e6*pc; Gpc=1e9*pc; yr=3.156e7
out=[]; P=lambda s="": out.append(s)
nu=lambda s: 1/(-np.expm1(-s))
rstar=lambda M: np.sqrt(G*M*Msun/a0)          # M in M☉
tauD=c/a0

# ======================= 1 Irrfahrt auf den Ebenen =======================
P("1a) Kette der Ebenen n = 0 … N (0 = Ort, unterste Ebene), Sprungrate κ je Richtung, Verlust Γ, G_n = E_n[e^{−ΓT₀}]")
P("    κ(G_{n+1} − 2G_n + G_{n−1}) = Γ·G_n  ⇒  G_n = λⁿ,  λ + 1/λ − 2 = Γ/κ;  L = m Ebenen, τ_D = m²/κ, s² = Γ·τ_D")
P("    q_m = G_m = exp(−2m·arsinh(s/2m)) → e^{−s} für m → ∞;  Kontinuum D = κh², L = m·h")

def solve_chain(m, s, N, bottom):
    kap=1.0; Gam=kap*(s/m)**2; n=N+1
    up=np.full(n,kap); lo=np.full(n,kap); dg=np.full(n,-(2*kap+Gam)); rhs=np.zeros(n)
    dg[-1]=-(kap+Gam)                                   # oberste Ebene reflektierend
    if bottom=="abs": dg[0]=1.0; up[0]=0.0; rhs[0]=1.0  # Ankunft am Ort = Rückkehr, G_0 = 1
    else: dg[0]=-(kap+Gam); rhs[0]=-1.0                 # Boden reflektierend, Quelle am Ort
    ab=np.zeros((3,n)); ab[0,1:]=up[:-1]; ab[1]=dg; ab[2,:-1]=lo[1:]
    Gv=solve_banded((1,1),ab,rhs); return Gv/Gv[0]
qm=lambda s,m: np.exp(-2*m*np.arcsinh(s/(2*m)))
dev=0
for s in (0.1,1.0,3.0):
    ga=solve_chain(50,s,20000,"abs"); gr=solve_chain(50,s,20000,"ref")
    dev=max(dev,abs(ga[50]/qm(s,50)-1),abs(gr[50]/ga[50]-1))
P(f"    Kette mit 20 000 Ebenen (Deckel reflektierend), m = 50, s = 0,1/1/3: G_m gegen λᵐ und Boden absorbierend gegen reflektierend mit Quelle: max. Abweichung {dev:.1e}")
P("    ⇒ Beide Böden geben dasselbe normierte Profil; der Boden legt nur die Amplitude fest.")
ss=np.logspace(-3,1.5,6000)
P("    m     max|ν_m/ν − 1|   bei s = 1   q_m/e^{−s} bei s = 10")
for m in (1,2,3,5,10,100):
    d=np.abs(nu(-np.log(qm(ss,m)))/nu(ss)-1)
    P(f"   {m:4d}   {d.max()*100:8.3f} %   {(nu(-np.log(qm(1.0,m)))/nu(1.0)-1)*100:7.3f} %   {qm(10.0,m)/np.exp(-10):10.3g}")
P("    ⇒ Schon wenige Ebenen treffen die Regel um a₀ auf Prozente; der Newton-Schwanz ist (m/s)^{2m} statt e^{−s}, e^{−s} bis s verlangt m ≫ s.")
P()

P("1b) „Nur aufwärts“ als Bewegung: zurück nur über eine geschlossene Kette von m Ebenen, q = (κ/(κ+Γ))^m, m/κ = τ_D")
P("       s    a/a₀ Kette (m = 10³)   a/a₀ Regel = s²·ν(s)")
for s in (1,0.3,0.1,0.01):
    m=1000; qd=(1+s*s/m)**(-m)
    P(f"   {s:6g}   {s*s/(1-qd):12.4f}   {s*s*nu(s):12.4f}")
P("    ⇒ 1 − q → s², a → a₀ (Boden, v² = a₀·r): keine Wurzel, kein Tully-Fisher.")
Pe=0.05
s=1e-4; up_=1-np.exp(-(Pe+np.sqrt(Pe**2+s*s))); dn_=-np.expm1(-(np.sqrt(Pe**2+s*s)-Pe))
P(f"    Drift v, Pe = vL/2D = {Pe}: aufwärts q = e^{{−(Pe + √(Pe² + s²))}}, ν ≤ 1/(1 − e^{{−2Pe}}) = {1/(1-np.exp(-2*Pe)):.2f} (außen Kepler);")
P(f"    abwärts q = e^{{−(√(Pe² + s²) − Pe)}}: a/a₀ bei s = 10⁻⁴: {s*s/dn_:.4f} = 2Pe (Boden). Nur die symmetrische Kette gibt 1 − q ≈ s.")
P()

P("1c) w als Glättungslänge: harmonische Fortsetzung (∇² + ∂_w²)u = 0 ist Poisson-Glättung, Symbol e^{−|k|w}")
P("    P_a·P_b = P_{a+b} (Glättungslängen addieren sich), P_0 = 1, |e^{−|k|w}| ≤ 1; w < 0 hieße e^{+|k|w} (unbeschränkt).")
P("    Schleife verlangte P_Λ = 1 für ein Λ > 0, also e^{−|k|Λ} = 1: nur k = 0 (der Mittelwert) kehrt voll zurück.")
u=lambda r,w: 2/(np.pi*r)*np.arctan(r/w)            # Fortsetzung von 1/r
h=1e-3; res=0
for (r,w) in ((1.0,0.5),(0.3,2.0),(3.0,1.0)):
    urr=(u(r+h,w)-2*u(r,w)+u(r-h,w))/h**2; ur=(u(r+h,w)-u(r-h,w))/(2*h); uww=(u(r,w+h)-2*u(r,w)+u(r,w-h))/h**2
    res=max(res,abs(urr+2*ur/r+uww)/abs(uww))
P(f"    Punktmasse: u(r,w) = (2/πr)·arctan(r/w), Laplace-Rest {res:.1e}, u(r, 10⁻⁹)·r = {u(1.0,1e-9):.9f}")
gP=lambda x: (2/np.pi)*(np.arctan(x)/x**2-1/(x*(1+x*x)))
opt=minimize_scalar(lambda x:-gP(x),bounds=(0.05,10),method="bounded"); alP=-opt.fun; alPl=2/(3*np.sqrt(3))
P(f"    Spitze des geglätteten Feldes α·GM/w²: gleichmäßige Kugel α = 1, Plummer α = {alPl:.4f}, Poisson α = {alP:.4f} (bei r = {opt.x:.3f}w)")
P(f"    Schwelle α·GM/L² = a₀ ⇒ L = √α·r*: 1 / {np.sqrt(alPl):.3f} / {np.sqrt(alP):.3f} r*;  Tully-Fisher v⁴ = GMa₀·(r*/L)²: ×1 / ×{1/alPl:.2f} / ×{1/alP:.2f}")
M=6e10; rs=rstar(M); r=10*kpc; Gam=G*M*Msun/r**2/c; D=G*M*Msun/c
P(f"    Punktmasse: D = L²/τ_D = GM/c = 3Δs·c, Γ = a_N/c ⇒ √(Γ/D)·r = {np.sqrt(Gam/D)*r:.6f}: die Abschirmzahl ist die Wellenzahl 1/r der Kreisbahnmode.")
P()

P("1d) Zahl der Ebenen: Schritt mit c ⇒ D = h·c/2 ⇒ h = 2GM/c² = 6Δs = r₀;  m = r*/r₀;  Weg in τ_D: c·τ_D = ℓ₀")
for name,M in (("Milchstraße",6e10),("Zwerg",1e8)):
    rs=rstar(M); r0=2*G*M*Msun/c**2; m=rs/r0
    P(f"    {name:12s} r* = {rs/kpc:7.3f} kpc, r₀ = {r0/1.496e11:9.4g} AE, m = {m:.3g}, Diskretheit s³/24m² bei s = 1: {1/(24*m*m):.1e}")
P()

# ======================= 2 Schließung =======================
P("2) Schließung oben: Deckel reflektierend bei Λ (geschlossener Raum), Ebenen durchlässig (D gleich)")
P("   q = cosh(k(Λ − L))/cosh(kΛ), k = s/L;  tief: 1 − q ≈ s·tanh(Λ/r) ⇒ 10 % bei Λ/r = artanh(1/1,1)")
H=2*np.pi*a0/c; RH=c/H; l0=c**2/a0
P(f"   a₀ = cH/2π: H = {H*Mpc/1e3:.1f} km/s/Mpc, c/H = {RH/Gpc:.2f} Gpc, U = 2πc/H = ℓ₀ = {l0/Gpc:.2f} Gpc; artanh(1/1,1) = {np.arctanh(1/1.1):.4f} ⇒ r₁₀ = {1/np.arctanh(1/1.1):.4f}·Λ")
def lnq_seg(s,L,Lam):
    k=s/L; return -s+np.log1p(np.exp(-2*k*(Lam-L)))-np.log1p(np.exp(-2*k*Lam))
def r10_seg(M,Lam):
    rs=rstar(M)
    f=lambda r: (-np.expm1(-rs/r))/(-np.expm1(lnq_seg(rs/r,rs,Lam)))-1.1
    return brentq(f,0.01*Lam,10*Lam,xtol=1e-6*Lam)
Lams=(("c/H",RH),("πc/H = U/2",np.pi*RH),("13,5 Gpc (U ≥ 27 Gpc)",13.5*Gpc))
P("   Λ                          Λ/L Milchstraße   Λ/L Zwerg   r₁₀ Milchstraße   r₁₀ Zwerg   a_N/a₀ bei r₁₀ (MS)")
for nm,Lam in Lams:
    rm=r10_seg(6e10,Lam); rz=r10_seg(1e8,Lam)
    P(f"   {nm:24s} {Lam/rstar(6e10):12.3g} {Lam/rstar(1e8):12.3g} {rm/Gpc:11.3f} Gpc {rz/Gpc:9.3f} Gpc {(rstar(6e10)/rm)**2:12.2e}")
T=tauD*(2*RH/rstar(6e10)-1)/2
P(f"   Mittlere Rückkehrzeit ⟨T⟩ = τ_D(2Λ/L − 1)/2, Milchstraße, Λ = c/H: {T/yr:.2e} Jahre (endlich; τ_D = {tauD/yr/1e9:.1f} Mrd. Jahre)")
GMH=c**3/(2*H); rsH=np.sqrt(GMH/a0); sH=rsH/RH
k=sH/rsH; Lam=np.pi*RH; qH=np.cosh(k*(Lam-rsH))/np.cosh(k*Lam)
P(f"   Spitze: Hubble-Kugel GM_H = c³/2H (Randneigung cH/2): L_H = r*_H = {rsH/RH:.4f}·c/H = √π·c/H; L_H/Λ = {rsH/RH:.2f} (Λ = c/H) bzw. {rsH/Lam:.3f} (Λ = πc/H)")
P(f"   ν_H offen = {nu(sH):.4f}, mit Deckel bei πc/H = {1/(1-qH):.4f}: für das Ganze ist L ≈ Λ, für jedes Teilsystem L/Λ ≤ 10⁻⁴ (Haufen 1,7·10¹⁴ M☉: {rstar(1.7e14)/RH:.1e}).")
P()

# ======================= 3 Umgebung =======================
P("3) Umgebung. (a) Einbettung im Verlust: Γ = a_N,ges/c (Gesamtfeld), L = r* (eigene Masse) und Λ unverändert.")
P("   (b) Ebenengrenze: Faltungen additiv ⇒ D_Ebene = 3c·ΣΔs = GM_Ebene/c; oberhalb der eigenen Glättungslänge ℓ ist D um μ = M_p/M größer.")
P("       Schritt mit c ⇒ Fluss stetig (D·G′ stetig); Γ je Zeit gleich ⇒ Z = √(DΓ), R = (1 − √μ)/(1 + √μ)")
P("       q = e^{−kL}·(1 + R·e^{−2k(ℓ−L)})/(1 + R·e^{−2kℓ}), k = s/L;  tief: ν_b/ν = (1 + R·e^{−2ℓ/r})/(1 − R·e^{−2ℓ/r}), R → −1: tanh(ℓ/r)")
def q_bd(s,L,ell,mu):
    k=s/L; R=(1-np.sqrt(mu))/(1+np.sqrt(mu))
    return np.exp(-s)*(1+R*np.exp(-2*k*(ell-L)))/(1+R*np.exp(-2*k*ell))
# Prüfung mit Fluss-Kette
s_,ell_,mu_=0.7,5.0,100.0; hh=0.005; W=ell_+400; N=int(W/hh); wn=np.arange(N+1)*hh
Db=np.where(np.arange(N)+0.5<ell_/hh,1.0,mu_)     # D auf den Kanten i+½
Gam=s_**2; n=N+1
dg=np.zeros(n); up=np.zeros(n); lo=np.zeros(n); rhs=np.zeros(n)
dg[0]=1; rhs[0]=1; dg[-1]=1
i=np.arange(1,N); dg[i]=-(Db[i]+Db[i-1])/hh**2-Gam; up[i]=Db[i]/hh**2; lo[i]=Db[i-1]/hh**2
ab=np.zeros((3,n)); ab[0,1:]=up[:-1]; ab[1]=dg; ab[2,:-1]=lo[1:]
Gv=solve_banded((1,1),ab,rhs)
P(f"   Prüfung Fluss-Kette (s = 0,7, ℓ = 5L, μ = 100): G(L) = {Gv[int(1/hh)]:.5f}, Formel {q_bd(s_,1.0,ell_,mu_):.5f}")
P("   Glättungslänge ℓ³ = GM/(√(3/2)·|E^TF|); Elternmasse als Punkt: √(3/2)·|E^TF| = 3GM_p/d³ ⇒ ℓ = d·(M/3M_p)^{1/3}")
Hs=np.diag([-2.0,1.0,1.0]); P(f"   (Prüfung: √(3/2)·|diag(−2,1,1)| = {np.sqrt(1.5)*np.linalg.norm(Hs):.6f})")
def efe(gi,ge,mode):
    if mode=="s": return nu(np.sqrt(np.hypot(gi,ge)/a0))*gi
    return nu(np.sqrt((gi+ge)/a0))*(gi+ge)-nu(np.sqrt(ge/a0))*ge
def setup(M,ge=None,d=None,Mp=None):
    if Mp is None: Mp=ge*d**2/(G*Msun)
    if ge is None: ge=G*Mp*Msun/d**2
    return Mp, ge, d*(M/(3*Mp))**(1/3)
cases=(("Milchstraße, g_e = 0,03a₀, d = 1 Mpc",6e10,dict(ge=0.03*a0,d=1*Mpc)),
       ("Zwerg 10⁸ M☉ bei 100 kpc von 6·10¹⁰ M☉",1e8,dict(Mp=6e10,d=100*kpc)))
def ratios(M,ge,ell,mu,r):
    rs=rstar(M); gi=G*M*Msun/r**2; free=nu(np.sqrt(gi/a0))*gi
    b=gi/(-np.expm1(np.log(q_bd(rs/r,rs,ell,mu))))
    st=np.sqrt(np.hypot(gi,ge)/a0); both=gi/(-np.expm1(np.log(q_bd(st,rs,ell,mu))))
    return efe(gi,ge,"s")/free, efe(gi,ge,"p")/free, b/free, both/free
for nm,M,kw in cases:
    Mp,ge,ell=setup(M,**kw); mu=Mp/M; rs=rstar(M); R=(1-np.sqrt(mu))/(1+np.sqrt(mu))
    f=lambda r: (-np.expm1(-rs/r))/(-np.expm1(np.log(q_bd(rs/r,rs,ell,mu))))-0.9
    r10=brentq(f,0.05*ell,5*ell)
    P(f"   {nm}: M_p = {Mp:.2e} M☉, g_e = {ge/a0:.4f}·a₀, μ = {mu:.0f}, R = {R:.3f}, ℓ = {ell/kpc:.1f} kpc = {ell/rs:.1f}·r*, r₁₀ (Grenze) = {r10/kpc:.1f} kpc = {r10/ell:.3f}·ℓ")
    P("        r/ℓ     r [kpc]   a/a_frei: äußeres Feld ⊥ / ∥    Ebenengrenze    beide (⊥)")
    for x in (0.3,0.5,0.7,1.0):
        e1,e2,b,bb=ratios(M,ge,ell,mu,x*ell)
        P(f"       {x:4.1f}   {x*ell/kpc:8.2f}        {e1:6.3f} / {e2:6.3f}         {b:6.3f}         {bb:6.3f}")
P("   Grenzwert R → −1: r₁₀ = ℓ/artanh(0,9) = %.3f·ℓ;  weit außen (r ≫ ℓ√μ) 1 − q → s·√μ: zweite, tiefere Ebene v_f·μ^{−1/4}" % (1/np.arctanh(0.9)))
def q_R(s,L,ell,R):
    k=s/L; return np.exp(-s)*(1+R*np.exp(-2*k*(ell-L)))/(1+R*np.exp(-2*k*ell))
for nm,M,kw in cases:
    Mp,ge,ell=setup(M,**kw); mu=Mp/M; rs=rstar(M); R=(1-mu)/(1+mu)
    f=lambda r: (-np.expm1(-rs/r))/(-np.expm1(np.log(q_R(rs/r,rs,ell,R))))-0.9
    P(f"   Variante Γ = D/r² je Ebene (gleiches k, Z = Dk): R = (1 − μ)/(1 + μ) = {R:.4f}, r₁₀ = {brentq(f,0.05*ell,5*ell)/ell:.3f}·ℓ ({nm.split(',')[0]})")
# Stern als unterste Ebene unter der Galaxie
M=6e10; rs=rstar(M); r=8*kpc; Gam=G*M*Msun/r**2/c; Ds=G*Msun/c; Dg=G*M*Msun/c; ls=0.93*pc
k1=np.sqrt(Gam/Ds); k2=np.sqrt(Gam/Dg); Z1=Ds*k1; Z2=Dg*k2
lnC=np.log(Z1)-(k1*ls+np.log((Z1+Z2)/2)+np.log1p(np.exp(-2*k1*ls)*(Z1-Z2)/(Z1+Z2)))
lnq=lnC-k2*(rs-ls)
P(f"   Liegt die Ebene des Sterns (1 M☉, ℓ = 0,93 pc) unter der Galaxie, gilt bei 8 kpc: k₁ℓ = {k1*ls:.1f}, q = 10^{lnq/np.log(10):.1f} statt {np.exp(-rs/r):.3f}:")
P("   die Rückkehr wäre gesperrt. Der Ort muss die Ebene des kollektiven (über die Zonen geglätteten) Feldes sein.")
P()

# ======================= 4 prüfbare Folge =======================
P("4) Gleiches äußeres Feld g_e = 0,03a₀, verschiedener Abstand zur Elternmasse (Milchstraße, r = 50 kpc)")
P("     d [Mpc]    M_p [M☉]     μ      ℓ [kpc]   Ebenengrenze a/a_frei   äußeres Feld ⊥ / ∥")
for d in (0.5,1.0,2.0):
    Mp,ge,ell=setup(6e10,ge=0.03*a0,d=d*Mpc); mu=Mp/6e10
    e1,e2,b,bb=ratios(6e10,ge,ell,mu,50*kpc)
    P(f"     {d:5.1f}   {Mp:10.2e}  {mu:7.0f}   {ell/kpc:7.1f}        {b:6.3f}              {e1:6.3f} / {e2:6.3f}")
P("   Das äußere Feld hängt nur an g_e; die Ebenengrenze an ℓ = (GMd/3g_e)^{1/3}, also am Abstand.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"halbgerade_hierarchie.txt"),"w",encoding="utf-8").write(txt); print(txt)
