# Geometrische Muster im Modell.
# 1 Fläche gegen Umfang: Flussdimension, Schnitt r*, scheinbare Masse, Umlaufphase, Linse, Kosmos
# 2 Herkunft der reinen Zahlen (Dimension, Schließung 2π, Ableitungsordnung, Bose-Summen, Modulgruppe)
# 3 Kreis als Unendlich-Eck: Polygon, Archimedes, Resonanzen p/q der Bahn, Kristallographie, symmetrische Tori
# 4 Wachsende Dimension: Kugel und Würfel, Rekursion 2π/n, Bose-Integral nach Dimension, Kac-Rückkehr
import numpy as np, sys, io, os, math
from math import gcd
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gamma, zeta
import mpmath as mp, warnings
np.seterr(over='ignore'); warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; AU=1.495978707e11
hbar=1.054571817e-34; kB=1.380649e-23
out=[]; P=lambda s="": out.append(s)
chi=lambda s: s/np.expm1(s)
nu=lambda s: 1/(-np.expm1(-s))
pi=np.pi

# =============== 1) Fläche gegen Umfang ===============
P("1a) Flussdimension d_eff = 1 − d ln a/d ln r für a = a_N/(1 − e^{−s}); Punktmasse, Einheiten GM = a₀ = 1 (r* = 1, s = r*/r)")
lna=lambda r: -2*np.log(r)+np.log(nu(1/r))
h=1e-5; rr=np.logspace(-2,2,401); dl=np.log1p(h)-np.log1p(-h)
deff=1-(lna(rr*(1+h))-lna(rr*(1-h)))/dl
P(f"   max|d_eff − (3 − χ(s))| = {max(abs(deff-(3-chi(1/rr)))):.1e}   (analytisch: d ln ν/d ln r = s/(e^s − 1) = χ)")
Om2=lambda r: np.exp(lna(r))/r
k2=(Om2(rr*(1+h))-Om2(rr*(1-h)))/(2*h)+4*Om2(rr)
P(f"   κ²/Ω² numerisch gegen 1 + χ: max Abw. {max(abs(k2/Om2(rr)-(1+chi(1/rr)))):.1e};  gegen 4 − d_eff: {max(abs(k2/Om2(rr)-(4-deff))):.1e}")
yy=np.logspace(-4,4,401); la=lambda y: np.log(y)+np.log(nu(np.sqrt(y)))
sl=(la(yy*(1+h))-la(yy*(1-h)))/dl
P(f"   d ln a/d ln a_N gegen 1 − χ/2: max Abw. {max(abs(sl-(1-chi(np.sqrt(yy))/2))):.1e}; Bereich [{sl.min():.4f}, {sl.max():.4f}] ⊂ (½, 1)")
P("   r/r*   d_eff   Re τ = √(1+χ)   ν = M_schein/M   r/r*   r/r* + ½   r/r* + ½ + r*/(12r)")
for x in (0.3,1,3,10,30):
    s=1/x
    P(f"   {x:4g}   {3-chi(s):.3f}   {np.sqrt(1+chi(s)):.4f}          {nu(s):.4f}         {x:.4f}   {x+0.5:.4f}     {x+0.5+1/(12*x):.4f}")
s=1.0
P(f"   Bernoulli: ν − (1/s + ½ + s/12 − s³/720 + s⁵/30240) bei s = 1: {nu(s)-(1/s+0.5+s/12-s**3/720+s**5/30240):.1e}")
M=6e10*Msun; GM=G*M
rs_num=brentq(lambda r: GM/r**2-np.sqrt(GM*a0)/r, 1e18, 1e22)
P(f"   Schnitt Flächengesetz GM/r² = Umfanggesetz √(GMa₀)/r (Milchstraße 6·10¹⁰ M☉): r = {rs_num/kpc:.3f} kpc = √(GM/a₀) = {np.sqrt(GM/a0)/kpc:.3f} kpc, a/a₀ dort = {GM/rs_num**2/a0:.6f}")
P("   Fluss: a·4πr² = 4πGM, a·2πr = 2π√(GMa₀) ⇒ a_Umfang/a_Fläche = (4πr²/2πr)/(4πr*²/2πr*) = r/r* (π kürzt sich)")
P()

P("1b) Potenzgesetz a ∝ r^{1−d}: κ²/Ω² = 4 − d, Winkel zwischen Perizentren 360°/√(4 − d)")
for d in (2,2.5,3,3.5):
    f=lambda r: r**(1-d)/r
    r=1.0; k2p=(f(r*(1+h))-f(r*(1-h)))/(2*h)+4*f(r)
    P(f"   d = {d:3g}: κ²/Ω² = {k2p/f(r):.6f}, Winkel {360/np.sqrt(k2p/f(r)):.2f}°")
P("   d = 4: κ² = 0 (keine stabile Kreisbahn), d > 4 instabil. d_eff ∈ [2, 3] ⇔ Re τ ∈ [1, √2].")
P()

P("1c) Linse: α = (2/c²)∫a_⊥dz. Für a ∝ r^{1−d}: α = (2/c²)·a(b)·b·F_d, F_d = ∫(1+u²)^{−d/2}du = S_{d−1}/S_{d−2}")
Sk=lambda k: 2*pi**((k+1)/2)/gamma((k+1)/2)   # Fläche der Einheitssphäre S^k ⊂ R^{k+1}
for d in (2,3,4,5):
    F=quad(lambda u: (1+u*u)**(-d/2),-np.inf,np.inf)[0]
    P(f"   d = {d}: F_d = {F:.6f},  S_{d-1}/S_{d-2} = {Sk(d-1)/Sk(d-2):.6f}")
P("   d = 3: α = 4GM/(c²b) = 12Δs/b;  d = 2: α = 2π(v/c)² = 6πΔs/r* (überstrichener Winkel π, mal 2 aus Φ + Ψ)")
def I(x):
    f=lambda z: nu(np.sqrt(x*x+z*z)**-1)/(x*x+z*z)**1.5
    return x*x*(quad(f,0,50*x,limit=400)[0]+quad(f,50*x,np.inf,limit=200)[0])
lnI=lambda x: np.log(I(x))
P("   b/r*   I = M_Zyl/M   πx/2 + ½   d_Linse = 3 − d ln I/d ln x   d_eff(r = b)   I/ν(r = b)")
for x in (0.3,1,3,10,30):
    dlens=3-(lnI(x*1.001)-lnI(x/1.001))/(2*np.log(1.001))
    P(f"   {x:4g}   {I(x):9.4f}   {pi*x/2+0.5:9.4f}   {dlens:.3f}                        {3-chi(1/x):.3f}          {I(x)/nu(1/x):.4f}")
P(f"   I/ν → π/2 = {pi/2:.4f}: Projektion von ρ ∝ 1/r² auf den Zylinder, kein eigenes Gesetz")
P()

P("1d) Kosmos: Hubble-Kugel mit kritischer Dichte, a₀ = cH/2π")
H0=67.4e3/3.0857e22; RH=c/H0; rho=3*H0**2/(8*pi*G); MH=4*pi/3*rho*RH**3
aNH=G*MH/RH**2; a0m=c*H0/(2*pi)
P(f"   a_N(R_H) = GM_H/R_H² = c²/(2R_H): {aNH/(c**2/(2*RH)):.6f};  a₀·2πR_H/c² = {a0m*2*pi*RH/c**2:.6f} (Umfanggesetz: c² durch den Hubble-Umfang)")
P(f"   a_N/a₀ = {aNH/a0m:.6f} = π;  r*_H/R_H = {np.sqrt(G*MH/a0m)/RH:.6f} = √π;  s = √π: ν = {nu(np.sqrt(pi)):.4f}, d_eff = {3-chi(np.sqrt(pi)):.3f}")
P(f"   lineare Störung a_N/a₀ = 10⁻³: Umfanggesetz gäbe ν = {nu(np.sqrt(1e-3)):.1f}; die Schließungsregel wirkt dort wie Staub mit Faktor 1 + Q = 6,36 (d = 3)")
P()

# =============== 2) Herkunft der Zahlen ===============
P("2a) 3 in Δs = GM/3c²: Überschussradius δR = ∫₀^R(1/√(1 − 2Gm/rc²) − 1)dr, in Einheiten GM/c² (GM/Rc² = 10⁻⁶)")
def dR(mu,eps=1e-6):
    return quad(lambda x: (1/np.sqrt(1-2*eps*mu(x)/x)-1) if x>0 else 0.0,0,1,limit=200)[0]/eps
for d in (2,3,4):
    P(f"   homogen in d = {d} Dimensionen (m ∝ r^d): δR = {dR(lambda x: x**d):.4f} GM/c²  (1/d = {1/d:.4f})")
P(f"   ρ ∝ 1/r (m ∝ r²): δR = {dR(lambda x: x**2):.4f} GM/c²")
# Polytrope n = 3
from scipy.integrate import solve_ivp
def le(xi,y,n=3): th,dth=y; return [dth,-max(th,0)**n-2*dth/xi]
x0=1e-6; ev=lambda xi,y: y[0]; ev.terminal=True
sol=solve_ivp(le,[x0,20],[1-x0**2/6,-x0/3],events=ev,rtol=1e-11,atol=1e-13,dense_output=True)
xi1=sol.t_events[0][0]; m1=xi1**2*abs(sol.sol(xi1)[1])
mu3=lambda x: (x*xi1)**2*abs(sol.sol(max(x*xi1,x0))[1])/m1
P(f"   Polytrope n = 3 (ξ₁ = {xi1:.5f}): δR = {quad(lambda x: mu3(x)/x,1e-6,1,limit=200)[0]:.4f} GM/c²  → die 3 ist die Dimension bei homogener Dichte")
P("2b) 6 in r₀ = 6Δs: g_tt = 1 − 2GM/(rc²) (der Takt ist die Wurzel der quadratischen Metrik): 6 = 2·3")
F3=quad(lambda z: 1/(1+z*z)**1.5,-np.inf,np.inf)[0]
P(f"2c) 12 in 12Δs/b: Linsenintegral {F3:.6f} (= S₂/S₁) × (1 + γ) = 2 × 3 (Δs) = {F3*2*3:.6f}")
P("2d) Exponent ½: Rückweg D·G'' = Γ·G, G(0) = 1, G beschränkt; Steigung d ln(−ln q)/d ln Γ")
def q2(Gam,L=1.0,D=1.0,X=None):
    k=np.sqrt(Gam/D); X=X or 40/k
    xs=np.linspace(0,X,4001); y0=np.vstack([np.exp(-k*xs),-k*np.exp(-k*xs)])
    s_=solve_bvp(lambda x,y: np.vstack([y[1],Gam/D*y[0]]), lambda ya,yb: np.array([ya[0]-1,yb[0]]), xs, y0, tol=1e-10, max_nodes=200000)
    return s_.sol(L)[0]
Gs=np.array([0.1,0.3,1,3,10]); qs=np.array([q2(g) for g in Gs])
P(f"   zweite Ordnung: q(Γ=1) = {q2(1.0):.8f} (e⁻¹ = {np.exp(-1):.8f}), Steigung {np.polyfit(np.log(Gs),np.log(-np.log(qs)),1)[0]:.5f}")
P("   erste Ordnung (Drift v·G' = −Γ·G): q = e^{−ΓL/v}, Steigung 1; Ordnung m: charakteristische Wurzel (Γ/D)^{1/m}")
f1=lambda t,L=1.0,D=1.0: L/np.sqrt(4*pi*D*t**3)*np.exp(-L*L/(4*D*t))
for Gam in (0.3,1.0,3.0):
    lap=quad(lambda t: f1(t)*np.exp(-Gam*t),0,np.inf,limit=200)[0]
    P(f"   Dichte erster Rückkehr L/√(4πDt³)·e^{{−L²/4Dt}} (∝ t^{{−3/2}}), Laplace bei Γ = {Gam}: {lap:.8f} = e^{{−√Γ}} = {np.exp(-np.sqrt(Gam)):.8f}")
P("2e) 4 in v⁴ = GMa₀: a = a_N^p·a₀^{1−p} ⇒ d ln v/d ln r = (1 − 2p)/2")
for p in (0.4,0.5,0.6):
    v=lambda r: np.sqrt(r*(GM/r**2)**p*a0**(1-p))
    r=10*kpc; P(f"   p = {p}: d ln v/d ln r = {(np.log(v(r*1.001))-np.log(v(r/1.001)))/(2*np.log(1.001)):+.4f}")
P(f"   p = ½: v² = √(GMa₀) ⇒ v⁴ = GMa₀; 4 = 2 (v² = ar) × 2 (Wurzel)")
P("2f) Bose-Integrale und Eisenstein-Zahlen: ∫x^{k−1}/(e^x − 1)dx = Γ(k)ζ(k);  (2π)^k/(Γ(k)ζ(k)) = 2/|ζ(1−k)|")
for k in (2,4,6):
    Ik=quad(lambda x: x**(k-1)/np.expm1(x),0,np.inf)[0]
    P(f"   k = {k}: Γ(k)ζ(k) = {Ik:.6f} ({gamma(k)*zeta(k):.6f});  (2π)^k/Γ(k)ζ(k) = {(2*pi)**k/Ik:.6f};  2/|ζ(1−k)| = {float(2/abs(mp.zeta(1-k))):.6f}")
P("   → 24 (E₂, Δ = η²⁴), 240 (E₄), 504 (E₆): Schließung (2π)^k geteilt durch Bose-Integral mit Zustandsdichte x^{k−1}")
sm=mp.mpf('1e-12'); qq=mp.e**(-sm)
Qt0=-4*(sm**3*mp.polylog(1,qq)+3*sm**2*mp.polylog(2,qq)+6*sm*mp.polylog(3,qq)+6*mp.polylog(4,qq))
P(f"   Q̃(s → 0) = {float(Qt0):.6f};  −4π⁴/15 = {-4*pi**4/15:.6f} = −4Γ(4)ζ(4) = −4(2π)⁴/240 = {-4*(2*pi)**4/240:.6f}")
sq=1.0; Qt1=-4*(sq**3*mp.polylog(1,mp.e**-sq)+3*sq**2*mp.polylog(2,mp.e**-sq)+6*sq*mp.polylog(3,mp.e**-sq)+6*mp.polylog(4,mp.e**-sq))
P(f"   Kontrolle s = 1: Li-Form {float(Qt1):.6f}, −4∫_s^∞x³n_B = {-4*quad(lambda x: x**3/np.expm1(x),1,np.inf)[0]:.6f}")
P(f"   Ω_s = (4π⁴/15)/(3·(2π)²) = {4*pi**4/15/(3*(2*pi)**2):.6f} = π²/45 = {pi**2/45:.6f}  (Bose-Summe π⁴ über Schließung (2π)² aus a₀ = cH/2π)")
P(f"   Nullpunkt: ζ(−1) = {float(mp.zeta(-1)):.6f} = −ζ(2)/(2π²); E₀ = (2πħc/L)·ζ(−1) = −πħc/(6L); 24 = (2π)²/ζ(2) = {(2*pi)**2/zeta(2):.6f}")
def partitions(n):
    p=[1]+[0]*n
    for m in range(1,n+1):
        s_=0; k=1
        while True:
            g1=k*(3*k-1)//2; g2=k*(3*k+1)//2
            if g1>m: break
            sg=1 if k%2 else -1
            s_+=sg*p[m-g1]
            if g2<=m: s_+=sg*p[m-g2]
            k+=1
        p[m]=s_
    return p[n]
for n in (500,5000):
    lp=math.log(partitions(n))
    P(f"   Zustandszahl p({n}): ln p = {lp:.4f};  2π√(n/6) − ln(4√3·n) = 2√(ζ(2)n) − … = {2*pi*np.sqrt(n/6)-np.log(4*np.sqrt(3)*n):.4f}  (6 = π²/ζ(2))")
P("2g) Spitzenform und Modulgruppe")
N=10
def sig(n,k): return sum(d**k for d in range(1,n+1) if n%d==0)
E4=np.array([1]+[240*sig(n,3) for n in range(1,N)],dtype=object)
E6=np.array([1]+[-504*sig(n,5) for n in range(1,N)],dtype=object)
def mul(a,b):
    r_=[0]*N
    for i in range(N):
        for j in range(N-i): r_[i+j]+=a[i]*b[j]
    return r_
E43=mul(mul(E4,E4),E4); E62=mul(E6,E6); D=[E43[i]-E62[i] for i in range(N)]
tau_=[0]*N; tau_[0]=1
for n in range(1,N):
    for _ in range(24):
        tau_=[tau_[i]-(tau_[i-n] if i>=n else 0) for i in range(N)]
Delta=[0]+tau_[:N-1]
P(f"   E₄³ − E₆² = {D[:5]} …: erster Koeffizient 3·240 + 2·504 = {3*240+2*504} = 12³ = {12**3}")
P(f"   (E₄³ − E₆²)/1728 = {[x//1728 for x in D[:7]]};  q∏(1 − qⁿ)²⁴ = {Delta[:7]}")
A_F=quad(lambda x: 1/np.sqrt(1-x*x),-0.5,0.5)[0]
P(f"   Fläche des Fundamentalbereichs ∫∫dx dy/y² = {A_F:.8f} = π/3 = 2π(1 − ½ − ⅓) = {2*pi*(1-0.5-1/3):.8f};  4π/Fläche = {4*pi/A_F:.6f}")
P("   Valenzformel: Nullstellen von f (Gewicht k) = k·Fläche/4π = k/12, Punkte i und ρ zählen ½ bzw. ⅓:")
P("   E₄: 4/12 = ⅓ (nur ρ), E₆: 6/12 = ½ (nur i), Δ: 12/12 = 1 (Spitze). 12 = 2/(1 − ½ − ⅓): aus den zwei symmetrischen Tori")
P("2h) Carter α_G ∝ α¹²(m_e/m_p)⁴: Beiträge zu α^{…} in α_G^{1/2}: Photosphäre T⁴ mit kT ∝ α² → 8, Zündung α² über T_c² → −4, Thomson α² → +2")
P(f"   12 = 2·(8 − 4 + 2) = {2*(8-4+2)}, 4 = 2·(4 − 2) = {2*(4-2)}: T⁴ (Stefan-Boltzmann, Bose-Integral π⁴/15), α² (Coulomb-Skalen), Verdopplung aus α_G^{{1/2}}; bei Kramers-Opazität 4 → 5")
Mm=1.0; rg=1e7  # r in GM/c²
Om=np.sqrt(1/rg**3); kap=Om*np.sqrt(1-6/rg)
P(f"2i) Perihel: Schwarzschild κ²/Ω² = 1 − 6GM/(rc²) = 1 − 18Δs/r; 2π(Ω/κ − 1) / (18πΔs/r) = {2*pi*(Om/kap-1)/(18*pi*(1/3)/rg):.6f}  (2π je Umlauf × 9)")
Nn=400; hh=1/Nn
ring=(2*np.eye(Nn)-np.roll(np.eye(Nn),1,0)-np.roll(np.eye(Nn),-1,0))/hh**2
Ni=Nn-1; hi=1/Nn; lin=(2*np.eye(Ni)-np.eye(Ni,k=1)-np.eye(Ni,k=-1))/hi**2
e_r=np.sort(np.linalg.eigvalsh(ring)); e_l=np.sort(np.linalg.eigvalsh(lin))
P(f"2j) Schließung: Lücke des Rings {e_r[1]:.4f} (4π² = {4*pi**2:.4f}), Rand {e_l[0]:.4f} (π² = {pi**2:.4f}), Verhältnis {e_r[1]/e_l[0]:.4f}")
alpha=1/137.035999
P(f"2k) Dirac: e·4πg = 2πħc·n (Kugelfluss = Schließung der Phase) ⇒ g/e = (2π/4π)/α = {0.5/alpha:.3f} = 1/(2α)")
Mb=1e12; T=hbar*c**3/(8*pi*G*Mb*kB); sig_SB=pi**2*kB**4/(60*hbar**3*c**2); A=16*pi*G**2*Mb**2/c**4
P(f"2l) Horizont: T = ħc/(24πkΔs) = ħκ/(2πkc) mit κ = c²/(12Δs); dM/dt = −σT⁴A/c² gegen −ħc⁴/(15360πG²M²): {(sig_SB*T**4*A/c**2)/(hbar*c**4/(15360*pi*G**2*Mb**2)):.6f}")
P(f"    15360π = 60·(8π)⁴/(16π·π²): 60 = Stefan-Boltzmann (Bose π⁴/15 und 4π²), 8π = 2π·4 (Schließung, κ), 16π = Horizontfläche; S = 4π(6Δs)²/4ℓ_P² = 36π(Δs/ℓ_P)²")
dip=quad(lambda th: np.sin(th)**3*2*pi,0,pi)[0]
P(f"2m) Eddington: ∫sin²θ dΩ = {dip:.6f} = 8π/3 (Thomson); L_Edd = 4π·3c²Δs·m_pc/((8π/3)r_e²) → Faktor {4*pi*3/(8*pi/3):.4f} = 9/2")
P(f"2n) reine Zahlen aus der Regel: ν(r*) = 1/(1 − e⁻¹) = {nu(1.0):.4f}; Rosette 360°/√2 = {360/np.sqrt(2):.2f}°; Hubble-Rand cH/2 = π·a₀")
P()

# =============== 3) Kreis als Unendlich-Eck ===============
P("3a) n-Eck (Umkreisradius 1): halber Umfang n·sin(π/n), umbeschrieben n·tan(π/n)")
P("   n       n·sin(π/n)     n·tan(π/n)     n²(π − n sin)   n²(n tan − π)")
for n in (6,12,24,96,1000,10**5):
    a_=n*np.sin(pi/n); b_=n*np.tan(pi/n)
    P(f"   {n:<7d} {a_:.10f}  {b_:.10f}  {n*n*(pi-a_):.5f}        {n*n*(b_-pi):.5f}")
P(f"   Grenzwerte π³/6 = {pi**3/6:.5f}, π³/3 = {pi**3/3:.5f}: Fehler ∝ 1/n²")
P("   Archimedes ab dem Sechseck: U_{2n} = 2U_nu_n/(U_n + u_n) (harmonisch), u_{2n} = √(U_{2n}u_n) (geometrisch)")
Uc,ui,n=2*np.sqrt(3),3.0,6; prev=None
for _ in range(5):
    Uc=2*Uc*ui/(Uc+ui); ui=np.sqrt(Uc*ui); n*=2; err=pi-ui
    P(f"   n = {n:3d}: {ui:.7f} < π < {Uc:.7f},  Fehlerverhältnis {(f'{prev/err:.4f}' if prev else '—')}"); prev=err
P(f"   n = 96: 223/71 = {223/71:.7f} < {ui:.7f},  {Uc:.7f} < 22/7 = {22/7:.7f}")
P("   Tangentenvieleck (Inkreisradius r = 1): Fläche/Umfang = n·tan(π/n)/(2n·tan(π/n)) = ½ für jedes n und für den Kreis (Kegelzerlegung, d: V = r·S/d)")
P()

P("3b) Umlaufphase der Bahn: Re τ_M = √(1 + χ(r*/r)), Milchstraße r* = %.2f kpc; Perizentren rücken je Radialperiode um 360°·q/p" % (np.sqrt(GM/a0)/kpc))
rsk=np.sqrt(GM/a0)/kpc
def s_of(x): return brentq(lambda s: chi(s)-(x*x-1), 1e-14, 745)
def r_of(x): return rsk/s_of(x)
fr=sorted({(p//gcd(p,q),q//gcd(p,q)) for q in range(1,8) for p in range(q+1,2*q) if 1<p/q<np.sqrt(2)}, key=lambda t: t[0]/t[1])
P("   p/q     Re τ      r [kpc]    Figur")
for p,q in fr:
    k_=p-q; k_=min(k_,p-k_); g_=gcd(p,k_)
    P(f"   {p}/{q}   {p/q:.5f}   {r_of(p/q):8.2f}   {{{p}/{k_}}}" + (" (Vieleck)" if k_==1 else " (Stern)"))
P("   irrational: Perizentren füllen den Ring dicht (Unendlich-Eck)")
P("   Häufung zur Kepler-Grenze 1: (n+1)/n")
for n in (5,10,30,100,1000):
    P(f"   n = {n:5d}: r = {r_of((n+1)/n):6.3f} kpc,  s = r*/r = {s_of((n+1)/n):.3f},  ln(n/2) + ln s = {np.log(n/2)+np.log(s_of((n+1)/n)):.3f}")
P("   Näherungsbrüche von √2 (Tiefgrenze):")
for p,q in ((7,5),(41,29),(239,169)):
    r=r_of(p/q); P(f"   {p}/{q}: √2 − p/q = {np.sqrt(2)-p/q:.3e},  r = {r:.4g} kpc = {r/rsk:.4g} r*,  r/(q²r*/2) = {r/(q*q*rsk/2):.4f}")
def mdist(xs,Q,lo=1.0,hi=1.5):
    fs=np.array(sorted({p/q for q in range(1,Q+1) for p in range(int(np.floor(lo*q)),int(np.ceil(hi*q))+1)}))
    i=np.clip(np.searchsorted(fs,xs),1,len(fs)-1)
    return np.minimum(abs(xs-fs[i-1]),abs(xs-fs[i])).mean()
xu=np.linspace(1,np.sqrt(2),400001)[1:-1]; Qs=(16,32,64,128,256,512)
md=np.array([mdist(xu,Q) for Q in Qs])
P("   mittlerer Abstand von Re τ ∈ [1, √2] zur nächsten Resonanz mit Nenner ≤ Q: " + ", ".join(f"Q={Q}: {m:.2e}" for Q,m in zip(Qs,md)))
P("   Abstand·Q²/ln Q = " + ", ".join(f"{m*Q*Q/np.log(Q):.3f}" for Q,m in zip(Qs,md)) + "  → fällt wie ln Q/Q²; beste Näherungen p/q einzeln: |x − p/q| < 1/q²")
cnt=len({(p,q) for q in range(1,51) for p in range(q+1,2*q) if gcd(p,q)==1 and 1<p/q<np.sqrt(2)})
P(f"   Zahl der Resonanzen mit q ≤ 50 in (1, √2): {cnt}; (√2 − 1)·3·50²/π² = {(np.sqrt(2)-1)*3*2500/pi**2:.1f}  (Dichte teilerfremder Paare 1/ζ(2) = 6/π²)")
xs_=np.sqrt(1+chi(rsk/8.2)); cf=[]; y=xs_
for _ in range(6): a_=int(np.floor(y)); cf.append(a_); y=1/(y-a_)
hp,hq=[1,cf[0]],[0,1]
for a_ in cf[1:]: hp.append(a_*hp[-1]+hp[-2]); hq.append(a_*hq[-1]+hq[-2])
P(f"   Sonnenort 8,2 kpc: Re τ = {xs_:.5f}, Kettenbruch {cf}; q²|x − p/q| = " + ", ".join(f"{hp[i]}/{hq[i]}: {hq[i]**2*abs(xs_-hp[i]/hq[i]):.3f}" for i in range(2,6)))
phi=(1+np.sqrt(5))/2; fib=[1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597]
P(f"   goldene Wicklung: n·‖nφ‖ bei Fibonacci-n → " + ", ".join(f"{n*abs(n*phi-round(n*phi)):.4f}" for n in fib[-4:]) + f";  1/√5 = {1/np.sqrt(5):.4f}")
P()

P("3c) Kristallographie: Drehung um 2π/n erhält ein Gitter nur, wenn ihre Spur 2cos(2π/n) ganz ist")
P("   erlaubt: n ∈ " + str([n for n in range(1,13) if abs(2*np.cos(2*pi/n)-round(2*np.cos(2*pi/n)))<1e-12]) + " → 360°, 180°, 120°, 90°, 60°")
def Ek(k,tau,Nt=60):
    q=mp.e**(2j*mp.pi*tau); Bk=mp.bernoulli(k)
    return 1-(2*k/Bk)*mp.fsum(sig(n,k-1)*q**n for n in range(1,Nt))
ti=mp.mpc(0,1); tr=mp.mpc(-0.5,mp.sqrt(3)/2)
P("   k:        " + "  ".join(f"{k:>9d}" for k in (4,6,8,10,12,14)))
P("   |E_k(i)|: " + "  ".join(f"{float(abs(Ek(k,ti))):9.2e}" for k in (4,6,8,10,12,14)) + "   (null genau für 4 ∤ k)")
P("   |E_k(ρ)|: " + "  ".join(f"{float(abs(Ek(k,tr))):9.2e}" for k in (4,6,8,10,12,14)) + "   (null genau für 6 ∤ k)")
jf=lambda t: 1728*Ek(4,t)**3/(Ek(4,t)**3-Ek(6,t)**2)
P(f"   j(i) = {float(mp.re(jf(ti))):.6f} (Quadrat, 90°),  |j(ρ)| = {float(abs(jf(tr))):.1e} (Sechseck, 60°/120°); 180° hat jedes Gitter (−1)")
P("   Drehung im Gitter: i·(1, i) = (i, −1), e^{iπ/3}·(1, ρ) = (1 + ρ, −1): ganzzahlig; nur diese zwei τ (mod SL(2,ℤ)) tragen 90° bzw. 60°/120°")
Ds=GM/(3*c**2); r8=8*kpc; ImT=np.log(r8/Ds)/(2*pi); ReT=np.sqrt(1+chi(rsk/8))
tm=mp.mpc(ReT-1,ImT)
P(f"   Bahn-Torus Milchstraße 8 kpc: Δs = {Ds/AU:.1f} AE, τ_M − 1 = {ReT-1:.3f} + {ImT:.3f}i, |j(τ_M)| = {float(abs(jf(tm))):.3e}, r/Δs = {r8/Ds:.3e}: an der Spitze, weit von i und ρ")
P()

# =============== 4) Wachsende Dimension ===============
P("4a) Einheitskugel V_n = π^{n/2}/Γ(n/2 + 1), Sphäre S_{n−1} = n·V_n, Würfel [−1,1]ⁿ")
V=lambda n: pi**(n/2)/gamma(n/2+1); S=lambda n: n*V(n)
P("   n    V_n       S_{n−1}    V_n/2ⁿ      Kac 2ⁿ/V_n   S/V Kugel  S/V Würfel(Inkreis 1)")
for n in (1,2,3,4,5,6,7,8,10,20):
    P(f"   {n:<3d} {V(n):8.4f}  {S(n):8.4f}   {V(n)/2**n:.3e}   {2**n/V(n):10.4g}   {S(n)/V(n):5.2f}      {2*n*2**(n-1)/2**n:5.2f}")
nV=minimize_scalar(lambda n: -V(n),bounds=(2,10),method='bounded',options={'xatol':1e-10}).x
nS=minimize_scalar(lambda n: -S(n),bounds=(2,12),method='bounded',options={'xatol':1e-10}).x
P(f"   Maximum V_n bei n = {nV:.4f}, S_{{n−1}} bei n = {nS:.4f}")
P(f"   V_n/V_(n−2) = 2π/n: max Abw. {max(abs(V(n)/V(n-2)-2*pi/n) for n in range(3,30)):.1e};  S_(n+1) = 2π·V_n: max Abw. {max(abs(S(n+2)-2*pi*V(n)) for n in range(1,30)):.1e}")
P(f"   Gauß: ∫e^{{−x²}}dx = {quad(lambda x: np.exp(-x*x),-np.inf,np.inf)[0]:.8f} = √π → π^{{n/2}}: ein π je zwei Dimensionen (eine Ebene, eine Schleife)")
P("4b) Bose-Integral in d Dimensionen (Zustandsdichte S_{d−1}x^{d−1}, Energie x): ∫x^d/(e^x − 1)dx = Γ(d+1)ζ(d+1)")
for d in (1,2,3,4,5):
    Id=quad(lambda x: x**d/np.expm1(x),0,np.inf)[0]
    form={1:"π²/6",2:"2ζ(3), keine π-Potenz",3:"π⁴/15",4:"24ζ(5), keine π-Potenz",5:"8π⁶/63"}[d]
    P(f"   d = {d}: {Id:.6f}  ({form})")
P(f"   8π⁶/63 = {8*pi**6/63:.6f}; π⁴/15 im Vakuumteil verlangt d = 3 (ungerades d gibt ζ(gerade) ∝ π^(d+1))")
P()

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"geometrische_muster.txt"),"w",encoding="utf-8").write(txt); print(txt)
