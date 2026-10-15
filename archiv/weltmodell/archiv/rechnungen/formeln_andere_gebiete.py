# Formeln anderer Gebiete gegen die Invarianten des Modells: Prüfung der Entsprechungen
# A Halbgerade gegen geschlossenes Ende (Thiele, Kabel, Warburg)   B Aufbauzeit (Morphogen) = Sattelzeit t*
# C Spektraldimension, Pólya, Watson   D kritische Verzweigung   E Entropie der Rücklaufzahl
# F Debye-Funktion   G Cauchy/Blanco-Mazzolo   H Akkumulationsverfahren an SPARC   I Matthiessen gegen Betrag
import numpy as np, sys, io, os, glob, math
from scipy.special import erfc, ellipk, zeta, gamma
from scipy.integrate import quad
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
nu=lambda s: 1/(1-np.exp(-s)); nB=lambda s: 1/np.expm1(s); chi=lambda s: s/np.expm1(s)

P("A) Dieselbe Gleichung D·G'' = Γ·G: offene Halbgerade gegen geschlossenes Ende bei Λ (Neumann, cosh-Form)")
P("   Punktmasse, s = r*/r, L = r*.  offen: q = e^{−s};  geschlossen: q = cosh((Λ/L − 1)s)/cosh((Λ/L)s)")
for lam in (10.,30.):
    qc=lambda s: np.cosh((lam-1)*s)/np.cosh(lam*s)
    f=lambda x: (1/(1-qc(1/x)))/(1/(1-np.exp(-1/x)))-1.10   # x = r/r*
    x10=brentq(f,1.,50*lam)
    P(f"   Λ/L = {lam:>4.0f}: Abweichung 10 % von e^(−s) bei r = {x10:.2f} r* = {x10/lam:.3f}·r*·Λ/L;  tief: 1 − q → s²(Λ/L − ½) (linear in Γ, a → konstant)")
s=1e-3; lam=10.; P(f"   Probe tief (Λ/L = 10, s = 10⁻³): (1 − q)/s² = {(1-np.cosh((lam-1)*s)/np.cosh(lam*s))/s**2:.4f}  (Λ/L − ½ = 9,5)")
P("   Warburg (Γ → iω): 1 − e^{−√(iωτ)} ≈ √(ωτ)·e^{iπ/4}, Phase 45°; geschlossen (coth-Form) unterhalb ω ≈ D/Λ² kapazitiv (90°).")
P()

P("B) Aufbauzeit der Rückkehr, falls sie ein Vorgang in t wäre (Synthese-Diffusion-Abbau, G(0,t) = 1 ab t = 0)")
P("   G(L,t) = ½[e^{−s}erfc(L/√(4Dt) − √(Γt)) + e^{s}erfc(L/√(4Dt) + √(Γt))],  τ_D = L²/D, Γ = s²/τ_D")
def G(s,t):   # t in Einheiten τ_D
    x=1/np.sqrt(4*t); y=s*np.sqrt(t)
    return 0.5*(np.exp(-s)*erfc(x-y)+np.exp(s)*erfc(x+y))
for s in (0.3,1.,3.):
    tau=quad(lambda t: 1-G(s,t)/np.exp(-s),0,np.inf,limit=400)[0]
    P(f"   s = {s}: mittlere Aufbauzeit ∫(1 − G/G_∞)dt = {tau:.5f} τ_D;  τ_D/(2s) = {1/(2*s):.5f};  Sattelzeit t* = ½√(τ_D/Γ) = {0.5/s:.5f}")
tD=2.998e8/1.2e-10/3.156e16   # Gyr
P(f"   τ_D = c/a₀ = {tD:.1f} Mrd. Jahre. Nach 13,8 Mrd. Jahren erreicht: s = 1: {G(1.,13.8/tD)/np.exp(-1):.2f}, s = 3: {G(3.,13.8/tD)/np.exp(-3):.2f}, s = 0,3: {G(0.3,13.8/tD)/np.exp(-0.3):.2f} des stationären q")
P()

P("C) Exponent ½ und Dimension der Rückkehrrichtung (Rückkehr eines Zufallswegs, z = e^{−Γ})")
P("   1D: 1 − F(z) = √(1 − z²);  2D (Quadratgitter): 1 − F = 1/P(z), P = (2/π)K(z);  3D: F(1) = 1 − 1/u₃ (Pólya)")
for Gm in (1e-2,1e-4,1e-6):
    z=np.exp(-Gm); one=np.sqrt(1-z*z); two=1/((2/np.pi)*ellipk(z*z))
    P(f"   Γ = {Gm:.0e}: 1D 1 − F = {one:.3e} (√(2Γ) = {np.sqrt(2*Gm):.3e});  2D 1 − F = {two:.4f} (π/ln(8/Γ) = {np.pi/np.log(8/Gm):.4f}, logarithmisch)")
u3=math.sqrt(6)/(32*math.pi**3)*math.gamma(1/24)*math.gamma(5/24)*math.gamma(7/24)*math.gamma(11/24)
P(f"   3D: Watson u₃ = {u3:.6f}, Rückkehrwahrscheinlichkeit 1 − 1/u₃ = {1-1/u3:.4f}; größte Verstärkung ν_max = u₃ = {u3:.4f}")
def slope(fun,G1,G2):  # d ln a/d ln a_N mit Γ ∝ a_N, a = a_N/(1 − F)
    return 1-(np.log(fun(G2))-np.log(fun(G1)))/(np.log(G2)-np.log(G1))
f1=lambda g: np.sqrt(1-np.exp(-2*g)); f2=lambda g: 1/((2/np.pi)*ellipk(np.exp(-2*g)))
P(f"   d ln a/d ln a_N tief: 1D {slope(f1,1e-6,1e-5):.4f} (½, Tully-Fisher), 2D {slope(f2,1e-6,1e-5):.4f} (→ 1, logarithmisch), 3D → 1 (ν beschränkt)")
P("   allgemein: Rückkehrdichte ∝ t^{d_s/2 − 2} ⇒ 1 − q ∝ Γ^{1 − d_s/2} ⇒ a ∝ a_N^{d_s/2}; flache Kurven verlangen d_s = 1.")
P()

P("D) Kritische Verzweigung (Poisson-Nachkommen, Mittel m): Gesamtzahl S, T(z) = E[z^S] = z·e^{m(T−1)}")
for m in (1.0,0.95):
    for Gm in (1e-3,1e-5):
        z=np.exp(-Gm); T=brentq(lambda T: T-z*np.exp(m*(T-1)),0,1-1e-15)
        P(f"   m = {m}: Γ = {Gm:.0e}: 1 − E[e^(−ΓS)] = {1-T:.4e};  √(2Γ) = {np.sqrt(2*Gm):.4e};  Γ/(1 − m) = {(f'{Gm/(1-m):.4e}' if m<1 else '–')}")
P("   m = 1: √-Gesetz (wie 1 − q ≈ s);  m < 1: linear in Γ (Boden, kein Tully-Fisher).")
P()

P("E) Rücklaufzahl n geometrisch verteilt, P(n) = (1 − q)qⁿ: Shannon-Entropie H = ln ν + χ")
for s in (0.3,1.,3.):
    q=np.exp(-s); n=np.arange(0,20000); p=(1-q)*q**n; p=p[p>0]
    P(f"   s = {s}: H = {-(p*np.log(p)).sum():.6f};  ln ν + χ = {np.log(nu(s))+chi(s):.6f};  Mittel Q = n_B = {nB(s):.4f}")
P()

P("F) Leiter bis zu einem Schnitt: ∫₀^x s^d n_B ds = (x^d/d)·D_d(x) (Debye-Funktion), x = s_max (Kompaktheit)")
for d in (1,2,3):
    row=[quad(lambda t: t**d/np.expm1(t),0,x)[0] for x in (0.3,1.,3.,30.)]
    P(f"   d = {d}: x = 0,3 / 1 / 3 / 30: " + " / ".join(f"{v:.4f}" for v in row) + f";  x^d/d bei 0,3: {0.3**d/d:.4f};  Γ(d+1)ζ(d+1) = {gamma(d+1)*zeta(d+1):.4f}")
P()

P("G) Mittlere Weglänge im Körper (Cauchy, Blanco-Mazzolo): ⟨ℓ⟩ = η_d·V/S, unabhängig von der Streuung")
for d in (1,2,3):
    P(f"   d = {d}: η_d = d·√π·Γ((d+1)/2)/Γ(d/2+1) = {d*math.sqrt(math.pi)*math.gamma((d+1)/2)/math.gamma(d/2+1):.5f}")
P("   d = 1, Halbgerade bei Λ reflektierend, Ein- und Austritt am Ort (V = Λ, S = 1): ⟨ℓ⟩ = 2Λ, ⟨T⟩ = 2Λ/c (Kac, masseunabhängig).")
# Monte Carlo: Zufallsweg auf 0..Λ, Start am Ort, Schritt ±h, reflektierend bei Λ, bis Rückkehr unter 0
rng=np.random.default_rng(1); Lam=1.0
for h in (0.1,0.05):
    N=int(round(Lam/h)); tot=0; runs=20000
    for _ in range(runs):
        x=0; k=0
        while True:
            x+= 1 if rng.random()<0.5 else -1; k+=1
            if x>N: x=N
            if x<0: break
        tot+=k
    P(f"   Zufallsweg, Schritt h = {h}: mittlere Weglänge bis zur Rückkehr = {tot/runs*h:.3f} Λ (Kac: 2(Λ + h))")
P()

P("H) Akkumulationsverfahren (Pharmakokinetik, Reaktor 1/M) an SPARC: s_beob = −ln(1 − g_bar/g_beob); Regel: s_beob² = g_bar/a₀")
base=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","daten","sparc")
tab={}
for line in open(os.path.join(base,"SPARC_Lelli2016c.mrt"),encoding="latin-1"):
    pa=line.split()
    if len(pa)>=18 and pa[1].isdigit() and pa[17].isdigit():
        tab[pa[0]]=(float(pa[5]),int(pa[17]))
gb=[];go=[]
for f in glob.glob(os.path.join(base,"Rotmod_LTG","*_rotmod.dat")):
    name=os.path.basename(f).replace("_rotmod.dat","")
    if name not in tab: continue
    inc,Q=tab[name]
    if Q>2 or inc<30: continue
    d=np.loadtxt(f,comments="#")
    R,V,eV,Vg,Vd,Vb=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]
    ok=(R>0)&(V>0)&(eV/V<0.1)
    k=1e6/3.0857e19
    gbar=(Vg*np.abs(Vg)+0.5*Vd**2+0.7*Vb**2)/R*k; gobs=V**2/R*k
    ok&=(gbar>0)&(gbar<gobs)
    gb+=list(gbar[ok]); go+=list(gobs[ok])
gb=np.array(gb); go=np.array(go); sb=-np.log(1-gb/go)
a0=np.median(gb/sb**2)
P(f"   Punkte: {len(gb)} (Q ≤ 2, Neigung ≥ 30°, δV/V < 0,1; Υ_disk = 0,5, Υ_bul = 0,7)")
P(f"   a₀ aus dem Median von g_bar/s_beob²: {a0:.3e} m/s²")
edges=np.logspace(-12,-8.5,8)
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(gb>=lo)&(gb<hi)
    if m.sum()<10: continue
    r=sb[m]**2*a0/gb[m]
    P(f"   g_bar {lo:.1e}–{hi:.1e}: n = {m.sum():4d}, Median s_beob²·a₀/g_bar = {np.median(r):.3f}  (16–84 %: {np.percentile(r,16):.2f}–{np.percentile(r,84):.2f})")
P("   Gerade durch den Ursprung heißt: die Exponentialform trägt; ein Anstieg bei kleinem g_bar wäre äußeres Feld (ν kleiner, s_beob größer).")
P("   Brauchbar nur für s ≲ 2 (g_bar ≲ 4a₀): für g_bar → g_beob divergiert −ln(1 − x), und die Auswahl g_bar < g_beob verzerrt nach unten.")
P()

P("I) Äußeres Feld: Betrag des Gesamtfelds (Modell) gegen unabhängige Verlustkanäle (Matthiessen: Raten addieren)")
GM=6.674e-11*6e10*1.989e30; r=50*3.0857e19; A0=1.2e-10
aN=GM/r**2; ae=0.03*A0
si=np.sqrt(aN/A0); se=np.sqrt(ae/A0)
P(f"   Milchstraße 50 kpc, g_e = 0,03 a₀: a_N = {aN/A0:.4f} a₀;  ν frei = {nu(si):.2f}")
P(f"   senkrecht s⁴ = s_i⁴ + s_e⁴: ν = {nu((si**4+se**4)**0.25):.2f};  gleichgerichtet s² = s_i² + s_e²: ν = {nu(np.sqrt(si**2+se**2)):.2f}")
P("   Matthiessen (Γ = Γ_i + Γ_e) gäbe für jede Richtung den gleichgerichteten Wert: Richtungsabhängigkeit trennt beide.")

txt="\n".join(out); print(txt)
open(os.path.splitext(os.path.abspath(__file__))[0]+".txt","w",encoding="utf-8").write(txt+"\n")
