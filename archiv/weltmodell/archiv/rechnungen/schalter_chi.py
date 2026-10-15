# Herleitung des Schalters χ = s/(e^s−1), s = √(a_N/a0): vier Ansätze rechnend geprüft
# (a) KAM/Chirikov  (b) Besetzung einer Mode  (c) Rückkehrstatistik  (d) Mittelung über Umläufe
import numpy as np, sys, io, os
from scipy import integrate
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; Gyr=3.156e16
M=6e10*Msun; GM=G*M; rs=np.sqrt(GM/a0); ds=GM/(3*c**2); vf=(GM*a0)**0.25
out=[]; P=lambda s="": out.append(s)
chi=lambda s: s/np.expm1(s)
nuB=lambda s: 1/(-np.expm1(-s))
rng=np.random.default_rng(20260928)

P("0) Zerlegung: a/a_N = ν = Σ_{n≥0} qⁿ = 1/(1−q),  χ = s·q/(1−q);  Tully-Fisher verlangt s·ν → 1 für s → 0 (1 − q ≈ s)")
P("   Tiefgrenzfall bei s = 0,01:")
for name,qf in [("q = e^{−s}",lambda s: np.exp(-s)),("q = 1/cosh s",lambda s: 1/np.cosh(s)),("q = e^{−s²}",lambda s: np.exp(-s*s))]:
    s=0.01; nu=1/(1-qf(s)); P(f"   {name:13s}  s·ν = {s*nu:9.4f}   a/a0 = s²ν = {s*s*nu:.4f}")
P()

# ---------------- (a) KAM / Chirikov ----------------
P("(a) Standardabbildung p' = p + K sin x, x' = x + p'  (Resonanz-Halbbreite 2√K, Abstand 2π; Überlappung bei K = π²/4)")
P("    Mit K ∝ ε = a0/a_N ist Breite/Abstand ∝ 1/s: die Wurzel entsteht aus der Pendelbreite ∝ √ε.")
def islands(K,n=3000,N=250):
    x=rng.uniform(0,2*np.pi,N*N); p=rng.uniform(-np.pi,np.pi,N*N)
    X=np.zeros_like(x); Y=np.zeros_like(x); mx=np.zeros_like(x); my=np.zeros_like(x)
    for _ in range(n):
        p=p+K*np.sin(x); p=(p+np.pi)%(2*np.pi)-np.pi; x=(x+p)%(2*np.pi)
        X+=p; Y+=np.mod(p,2*np.pi)-np.pi
        np.maximum(mx,np.abs(X),out=mx); np.maximum(my,np.abs(Y),out=my)
    return np.mean(mx<2*np.pi), np.mean(my<2*np.pi)
P("  a1) Anteil nicht auf rotierenden KAM-Tori (Inseln 0/1 und 1/2, Librationstest) bei kleinem K:")
P("      K      Insel 0/1   /√K      Pendel 16√K/4π²   Insel 1/2   /K")
for K in [0.01,0.04,0.1,0.3]:
    f1,f2=islands(K)
    P(f"   {K:6.2f}   {f1:.4f}    {f1/np.sqrt(K):.3f}    {16*np.sqrt(K)/(4*np.pi**2):.4f}          {f2:.4f}     {f2/K:.3f}")
P("      Anteil ≈ 0,405·√K ∝ 1/s: Zusatz χ·√(a_N a0) = Anteil·a0·s → konstant ≈ 0,64·a0 im ganzen Newton-Bereich (K = π²/4s²).")
def chaos(K,n=2500,N=250,thr=0.02):
    x=rng.uniform(0,2*np.pi,N*N); p=rng.uniform(-np.pi,np.pi,N*N)
    dx=np.ones_like(x); dp=np.zeros_like(x); ls=np.zeros_like(x)
    for _ in range(n):
        cK=K*np.cos(x); dp=dp+cK*dx; dx=dx+dp
        p=p+K*np.sin(x); x=(x+p)%(2*np.pi)
        nr=np.hypot(dx,dp); ls+=np.log(nr); dx/=nr; dp/=nr
    return np.mean(ls/n>thr)
P("  a2) Chaotischer Anteil (Lyapunov > 0,02 nach 2500 Schritten) gegen χ(s), s = π/(2√K) (s = 1 bei Überlappung):")
P("      K       s      chaotisch   χ(s)")
for K in [0.6,0.8,0.9716,1.2,1.5,2,3,5,8]:
    s=np.pi/(2*np.sqrt(K)); P(f"   {K:6.3f}  {s:5.3f}   {chaos(K):.4f}     {chi(s):.4f}")
P("      Grenzschicht für kleines K ∝ exp(−π²/√K) = exp(−2π·s); Koeffizient 1 hieße K(s=1) = π⁴ = 97,4 (dort längst voll chaotisch).")
P("  a3) Melnikov-Integrale des Pendels (Separatrix ẋ₀ = 2 sech t), numerisch gegen geschlossen:")
for W in [0.5,1,2,4]:
    I1=2*integrate.quad(lambda t: 2/np.cosh(t),0,np.inf,weight='cos',wvar=W)[0]
    I2=2*integrate.quad(lambda t: np.tanh(t)/np.cosh(t)**2,0,np.inf,weight='sin',wvar=W)[0]
    P(f"   Ω = {W}:  Antrieb {I1:.10f} = 2π·sech(πΩ/2) {2*np.pi/np.cosh(np.pi*W/2):.10f};   parametrisch {I2:.10f} = πΩ²/(2 sinh(πΩ/2)) {np.pi*W**2/(2*np.sinh(np.pi*W/2)):.10f}")
P("      Exponent immer πΩ/2 (Pole des sech-Pulses bei ±iπ/2), Vorfaktor je nach Störung verschieden.")
P()

# ---------------- (b) Besetzung einer Mode ----------------
P("(b) Mode der Frequenz ω bei einem in imaginärer Zeit periodischen Korrelator (Periode 2π, Singularität bei Koinzidenz):")
t=np.linspace(-80,80,1600001); eps=0.5
P("      ω      ∫e^{−iωt}·[−1/(4sinh²((t−iε)/2))]dt·e^{−ωε}   s/(e^s−1), s=2πω    FT[1/(4cosh²(t/2))]   Quotient   e^{−s/2}")
for w in [0.05,0.15,0.3,1.0]:
    I=np.trapezoid(np.exp(-1j*w*t)*(-1/(4*np.sinh((t-1j*eps)/2)**2)),t).real*np.exp(-w*eps)
    J=np.trapezoid(np.cos(w*t)/(4*np.cosh(t/2)**2),t); s=2*np.pi*w
    P(f"   {w:5.2f}    {I:.12f}                         {chi(s):.12f}      {J:.10f}          {chi(s)/J:.8f}   {np.exp(-s/2):.8f}")
P("      Singularität auf der reellen Achse (thermisch) → Bose s/(e^s−1); Separatrix-Puls (Pole bei halber Periode) → (s/2)/sinh(s/2).")
P("      χ = ⟨E⟩/kT (Planck-Mittelwert durch Gleichverteilungswert), a = a_N·(1 + n_B): Identitäten, keine Herleitung von s.")
P("  b1) Natürliches Verhältnis Energie/Temperatur bei Unruh- bzw. Horizonttemperatur ist linear: s' = a_N/a0.")
P("      Punktmasse 6·10¹⁰ M☉, v_f = (GMa0)^{1/4} = %.1f km/s:" % (vf/1e3))
P("      r [kpc]   v mit q = e^{−a_N/a0}   v mit q = e^{−√(a_N/a0)}")
for R in [10,30,100,300]:
    r=R*kpc; aN=GM/r**2
    vl=np.sqrt(aN/(-np.expm1(-aN/a0))*r); vb=np.sqrt(aN*nuB(np.sqrt(aN/a0))*r)
    P(f"   {R:6d}      {vl/1e3:8.1f} km/s            {vb/1e3:8.1f} km/s")
P("      Linear: a → a0 + a_N/2, v ∝ √r steigt, keine Tully-Fisher-Relation.")
P("  b2) Energie = ħ·Bahnfrequenz, Temperatur kT_GH = ħa0/c:")
for R in [8,50]:
    r=R*kpc; wK=np.sqrt(GM/r**3); s=np.sqrt(GM/r**2/a0)
    P(f"      r = {R} kpc: ħω_K/kT_GH = ω_K·c/a0 = {wK*c/a0:8.1f},   verlangt s = {s:.3f}")
wreq=np.sqrt(a0*a0)/c
P("      Verlangt ist ω = √(a_N a0)/c = √(a_N/ℓ0), ℓ0 = c²/a0: bei a_N = a0 Periode 2πc/a0 = %.0f Mrd. Jahre." % (2*np.pi/wreq/Gyr))
P()

# ---------------- (c) Rückkehrstatistik ----------------
P("(c) Rückkehr mit poissonschem Verlust (Rate Γ):")
Tl=1.0; esc=rng.exponential(1.0,2_000_000); nl=np.floor(esc/Tl)
P("  c1) Fester Umlauf T, ΓT = 1: mittlere Zahl vollständiger Umläufe (MC) %.5f, q/(1−q) = %.5f (Kac)." % (nl.mean(), np.exp(-1)/(-np.expm1(-1))))
P("      Mit Γ ∝ a_N und festem T folgt q = e^{−a_N/a0} (wie b1). Mit Γ = √(a_N/L), T = √(L/a0) folgt s exakt (gleiche Länge L gesetzt).")
P("      Verlust längs des Weges nach außen, Rate ∝ a_N: s_Weg = |Φ_N|/√(GMa0). Plummer b = 3 kpc, 6·10¹⁰ M☉:")
b=3*kpc
for R in [0.5,2,5,10,30]:
    r=R*kpc; aN=GM*r/(r*r+b*b)**1.5; Phi=GM/np.sqrt(r*r+b*b)
    sl=np.sqrt(aN/a0); sw=Phi/np.sqrt(GM*a0)
    P(f"      r = {R:4.1f} kpc:  s_lokal = {sl:.4f}  s_Weg = {sw:.4f}   ν_lokal = {nuB(sl):.4f}  ν_Weg = {nuB(sw):.4f}")
def fp(L,W,Tmax,two=False,B=512):
    pos=np.full(W,0 if two else L,dtype=np.int64); T=np.full(W,-1,dtype=np.int64); idx=np.arange(W); tt=0
    while tt<Tmax and idx.size:
        st=rng.integers(0,2,size=(idx.size,B),dtype=np.int8).astype(np.int64)*2-1
        path=pos[idx,None]+np.cumsum(st,axis=1)
        hit=(np.abs(path)>=L) if two else (path<=0)
        ah=hit.any(axis=1); fi=hit.argmax(axis=1)
        T[idx[ah]]=tt+fi[ah]+1; pos[idx]=path[:,-1]; idx=idx[~ah]; tt+=B
    return T
L=12; W=40000
T1=fp(L,W,60000); T2=fp(2*L,W,120000); T3=fp(L,W,20000,two=True)
P("  c2) Diffusive Rückkehr (Gitter-Irrfahrt, L = 12, %d Läufer), Verlust je Schritt 1−e^{−Γ}, s = L·√(2Γ):" % W)
P("      s     einseitig MC   Gitter z^L  e^{−s}     von 2L MC   e^{−2s}    geschlossener Kreis (±L) MC   1/cosh s")
for s in [0.3,1,2,4]:
    Gm=(s/L)**2/2; p=-np.expm1(-Gm); z=(1-np.sqrt(1-(1-p)**2))/(1-p)
    w=lambda T: np.where(T>0,(1-p)**np.maximum(T,0),0.0)
    P(f"   {s:4.1f}    {w(T1).mean():.4f}±{w(T1).std()/np.sqrt(W):.4f}  {z**L:.4f}      {np.exp(-s):.4f}    {w(T2).mean():.5f}    {np.exp(-2*s):.5f}           {w(T3).mean():.4f}±{w(T3).std()/np.sqrt(W):.4f}         {1/np.cosh(s):.4f}")
P("      Einseitig (nach außen offen): q = e^{−s}, n Rückkehren qⁿ, Zahl der Rückkehren ν = 1/(1−q). Geschlossener Kreis: 1/cosh s → a → 2a0 (Abschnitt 0).")
P("  c3) Laplace-Eindeutigkeit: E[e^{−λX}] = e^{−√λ} für alle λ nur für die Lévy-Dichte f(x) = x^{−3/2}e^{−1/(4x)}/(2√π) (erste Ankunft der Diffusion):")
for lam in [0.1,1,10]:
    I=integrate.quad(lambda x: np.exp(-lam*x)*x**-1.5*np.exp(-1/(4*x))/(2*np.sqrt(np.pi)),0,np.inf,limit=400)[0]
    P(f"      λ = {lam:5.1f}:  {I:.12f}   e^(−√λ) = {np.exp(-np.sqrt(lam)):.12f}")
P("      Verlust ∝ a_N^β, Rückkehrdauer stabil mit Index α: q = e^{−(a_N/a0)^{αβ}}; verlangt αβ = ½:")
P("      α = 1 (feste Dauer) ↔ β = ½ (Frequenz ∝ √a);  α = ½ (Diffusion) ↔ β = 1 (Rate ∝ a, wie Unruh/Horizont).")
P("  c4) Zeitskala bei Γ = a_N/c, τ_D = c/a0 = %.0f Mrd. Jahre: mittlere Dauer überlebender Rückkehren τ_D/(2s), bei s = 1: %.0f Mrd. Jahre." % (c/a0/Gyr, c/a0/2/Gyr))
P("      Als Echtzeit-Geschichte zu lang (mit a0 = cH0/2π ist c/(2a0) = π/H0); als stationärer Kern: q = e^{−κL} löst G'' = κ²G, κ²L² = a_N/a0.")
P()

# ---------------- (d) Mittelung über Umläufe ----------------
P("(d) Mittelung")
P("  d1) Kern-Nome Δs/r, zeitgemittelt über Kepler-Ellipsen mit großer Halbachse 8 kpc (Milchstraße 6·10¹⁰ M☉):")
Ea=np.linspace(0,2*np.pi,200001)
for e in [0,0.5,0.9]:
    r=8*kpc*(1-e*np.cos(Ea)); wt=(1-e*np.cos(Ea))
    P(f"      e = {e}:  ⟨Δs/r⟩ = {np.trapezoid(wt*ds/r,Ea)/np.trapezoid(wt,Ea):.4e}")
P("      verlangt q = e^{−s} bei 8 kpc: %.4f. Ein Potenzgesetz in Δs/r bleibt nach jeder Mittelung ein Potenzgesetz der Größe 10⁻⁷." % np.exp(-rs/(8*kpc)))
P("  d2) Summe über alle Harmonischen statt der einfachen Schleife, (s·ν)² = v⁴/(GMa0) bei der Punktmasse:")
P("      s        Bose (n=1)   Lambert 1+Σ1/(e^{ns}−1)   1+Σσ₁(n)e^{−ns}")
for s in [0.1,0.01,0.001]:
    n=np.arange(1,int(60/s)+1)
    nuL=1+np.sum(1/np.expm1(n*s)); nuS=1+np.sum(n*np.exp(-n*s)/(-np.expm1(-n*s)))
    P(f"   {s:6.3f}    {(s*nuB(s))**2:.4f}       {(s*nuL)**2:10.3f}                {(s*nuS)**2:.4e}")
P("      Nur die einfache Schleife hält v⁴ = GMa0; mit allen Harmonischen wächst v⁴ logarithmisch bzw. wie 1/s².")
P("  d3) Frequenzdarstellung: χ(s) = 1 − s/2 + 2Σ_k s²/(s² + (2πk)²) (Summe von Speicheranteilen mit Raten 2πk):")
k=np.arange(1,2_000_001)
for s in [0.5,2,5]:
    S=1-s/2+2*np.sum(s*s/(s*s+(2*np.pi*k)**2))+2*s*s/(4*np.pi**2*k[-1])
    P(f"      s = {s}:  Reihe {S:.10f}   χ = {chi(s):.10f}")
P()
P("Ergebnis: (a) liefert die Wurzel (Resonanzbreite ∝ √ε), nicht die Form. (b) liefert die Bose-Form, nicht die Wurzel.")
P("(c) einseitig diffusiv + Verlust ∝ a_N liefert q = e^{−√(a_N/a0)} exakt samt Tully-Fisher; ballistisch oder geschlossen nicht.")
P("(d) führt über die Laplace-Eindeutigkeit auf (c) zurück; Mittelung des Kern-Nomes oder aller Harmonischen trägt nicht.")
txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"schalter_chi.txt"),"w",encoding="utf-8").write(txt); print(txt)
