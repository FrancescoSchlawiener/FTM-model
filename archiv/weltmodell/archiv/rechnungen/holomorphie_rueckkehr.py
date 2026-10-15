# Holomorphie der Rückkehr: (a) Halbgerade streng aus (∇² + ∂_w²)δσ = 0, (b) warum nur der holomorphe Teil zählt.
# A1 Moden auf w ≥ 0 und auf der Schleife   A2 Halbgruppe, Positivität, Kontraktion   A3 Eindeutigkeit ohne Daten im Unendlichen
# A4 Schleife mit Quelle am Ort   A5 ortsabhängiger Operator A   A6 zeitartiges w
# B1 Szegő, Poisson, Schwarz   B2 Kern des Autors und Gittersumme   B3 Cauchy-Riemann-Paar   B4 Hardy-Projektion, Orientierung
# B5 Gewicht je Rücklauf (Erneuerung, Kanäle, Grenzwert χ → 1)   B6 Zahlen
import numpy as np, sys, io, os, sympy as sy
from scipy.linalg import expm, eigh
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19
out=[]; P=lambda s="": out.append(s)
nB=lambda s: 1/np.expm1(s); coth=lambda x: 1/np.tanh(x)

# ================= A1 Moden =================
P("A1) Mode f(w)·e^{ik·x} von (∇² + ∂_w²)u = 0:  f'' = k²f")
w,k,Lam=sy.symbols('w k Lambda',positive=True); f=sy.Function('f')
sol=sy.dsolve(sy.Eq(f(w).diff(w,2),k**2*f(w)))
P(f"    allgemein: {sol.rhs}   (k > 0)")
P("    beschränkt für w → ∞ ⇔ Koeffizient von e^{+kw} = 0 ⇒ f = f(0)·e^{−kw};  k = 0: f = A + Bw, beschränkt ⇔ B = 0")
Mono=sy.Matrix([[0,1],[k**2,0]])
ev=[sy.exp(Lam*e) for e in Mono.eigenvals().keys()]
P(f"    Schleife w ~ w + Λ: Monodromie exp(Λ·[[0,1],[k²,0]]) hat die Eigenwerte {ev}; = 1 nur für k = 0")
P(f"    k = 0: Monodromie [[1,Λ],[0,1]] ≠ I, periodisch nur f = konst.")
P()

# ================= A2 Halbgruppe =================
P("A2) Glättung P_w: Symbol e^{−|k|w} (x periodisch, 512 Punkte, Zufallsdaten)")
rng=np.random.default_rng(1); N=512; kk=np.fft.fftfreq(N,1/N)*2*np.pi/ (2*np.pi)   # ganze Wellenzahlen, Umfang 2π
f0=rng.normal(size=N); F0=np.fft.fft(f0)
Pw=lambda Fh,ww: Fh*np.exp(-np.abs(kk)*ww)
for a,b in ((0.1,0.3),(0.5,1.2),(2.0,3.0)):
    d=np.abs(np.fft.ifft(Pw(Pw(F0,a),b))-np.fft.ifft(Pw(F0,a+b))).max()
    P(f"    P_{a}·P_{b} − P_{a+b}: {d:.1e}")
xs=np.linspace(0,2*np.pi,N,endpoint=False)
for ww in (0.05,0.5,2.0):
    ker=(np.sinh(ww))/(np.cosh(ww)-np.cos(xs))/(2*np.pi)             # Poisson-Kern des Kreises, Masse 1
    fw=np.fft.ifft(Pw(F0,ww)).real
    P(f"    w = {ww}: Kern min {ker.min():.3e} ≥ 0, Masse {ker.sum()*2*np.pi/N:.12f};  max|P_w f| = {abs(fw).max():.3f} ≤ max|f| = {abs(f0).max():.3f}")
P("    Satz: Eine Familie T_w (w ≥ 0) mit T_0 = I, T_a T_b = T_{a+b}, verschiebungsinvariant in x, u = T_w f harmonisch in (x, w)")
P("    und max|T_w f| ≤ max|f| ist T_w = P_w: Multiplikator m_w(k) = e^{−wφ(k)} (Halbgruppe), φ² = k² (harmonisch), φ = +|k| (Kontraktion).")
P("    T_Λ = I verlangt e^{−|k|Λ} = 1 ⇒ k = 0. Fortsetzung auf w < 0 hätte |m| = e^{+|k||w|} > 1: keine Glättung, unbeschränkt.")
P()

# ================= A3 Eindeutigkeit =================
P("A3) Eindeutigkeit: Hälfte des Satzes ohne Halbgruppe. Beschränkt harmonisch auf w ≥ 0 mit u(x,0) = 0 ⇒ u ≡ 0")
P("    (ungerade Spiegelung w → −w gibt eine beschränkte harmonische Funktion auf dem ganzen Raum, Liouville ⇒ konstant = 0).")
P("    Ohne Beschränktheit gibt es u = w, u = w·x₁, u = sinh(kw)e^{ikx}: die Daten am Ort legen das Feld dann nicht fest.")
P("    Endliche Tiefe W mit Deckel (kL = s = 1, G bei w = L):")
for W in (2,5,20):
    gD=np.sinh(W-1)/np.sinh(W); gN=np.cosh(W-1)/np.cosh(W)
    P(f"      W = {W}L: Deckel u = 0 → {gD:.6f}, Deckel ∂u = 0 → {gN:.6f}, Halbgerade e^{{−1}} = {np.exp(-1):.6f}")
P()

# ================= A4 Schleife mit Quelle =================
P("A4) Schleife w ~ w + Λ mit dem Ort als Quelle (Sprung von ∂_w u am Ort): G(w) = cosh(k(w − Λ/2))/cosh(kΛ/2)")
for LamL in (2.0,4.0,10.0):
    kk_=1.0; Gl=lambda x: np.cosh(kk_*(x-LamL/2))/np.cosh(kk_*LamL/2)
    ws=np.linspace(0,LamL,401); A,B=np.meshgrid(ws,ws); msk=A+B<=LamL
    dS=np.abs(Gl(A)*Gl(B)-Gl(A+B))[msk].max()
    jump=2*kk_*np.tanh(kk_*LamL/2)
    P(f"    Λ = {LamL:4.1f}/k: max|G(a)G(b) − G(a+b)| = {dS:.3f};  G(Λ/2) = {Gl(LamL/2):.4f}, G(3Λ/4)/G(Λ/2) = {Gl(0.75*LamL)/Gl(LamL/2):.3f} > 1;"
      f"  Quellstärke −[∂_wG] = {jump:.4f} (Halbgerade: 1)")
P("    → Mit dem Ort als Quelle verträgt die Gleichung eine Schleife; dann ist G keine Halbgruppe und wächst jenseits Λ/2 wieder (Entglättung).")
P("      Mit dem Ort als Rand (Daten G(0) = 1, Glättungstiefe) gibt es keine Schleife. Das ist dieselbe Form wie der thermische Propagator im Kreis.")
P()

# ================= A5 ortsabhängiger Operator =================
P("A5) Ortsabhängige Kopplung: u_ww = A·u, A = −∂_x ν(x) ∂_x ≥ 0 (periodisch, ν > 0 zufällig, 200 Punkte)")
n=200; h=2*np.pi/n; xg=np.arange(n)*h
nuf=1.5+np.sin(xg)+0.3*np.cos(3*xg+1); nuh=0.5*(nuf+np.roll(nuf,-1))
Am=np.zeros((n,n))
for i in range(n):
    j=(i+1)%n; Am[i,i]+=nuh[i]/h**2; Am[j,j]+=nuh[i]/h**2; Am[i,j]-=nuh[i]/h**2; Am[j,i]-=nuh[i]/h**2
lam,V=eigh(Am); lam[np.abs(lam)<1e-10]=0.0; lam=np.clip(lam,0,None); Pa=lambda ww: (V*np.exp(-ww*np.sqrt(lam)))@V.T
P(f"    kleinste Eigenwerte von A: {lam[0]:.1e}, {lam[1]:.4f}, {lam[2]:.4f};  ‖P_0,3·P_0,7 − P_1‖ = {np.abs(Pa(0.3)@Pa(0.7)-Pa(1.0)).max():.1e}")
for LL in (1.0,10.0):
    P(f"    Λ = {LL}: Zahl der Eigenwerte von P_Λ gleich 1 (auf 10⁻⁹): {int(np.sum(np.abs(np.exp(-LL*np.sqrt(lam))-1)<1e-9))} (nur der Mittelwert)")
P("    → Homogenität ist nicht nötig: P_w = e^{−w√A} für jedes A ≥ 0; eine Schleife lässt nur den Kern von A (Mittelwert) zu.")
P()

# ================= A6 zeitartiges w =================
P("A6) Vorzeichen: wäre w zeitartig (h^{ww} = −1), gälte u_ww = −A·u mit Lösungen e^{±i√λ w}: beschränkt, keine Halbgruppe,")
P(f"    periodisch für √λ·Λ ∈ 2πℤ (z. B. λ = 1, Λ = 2π: max|M − I| = {np.abs(expm(2*np.pi*np.array([[0,1],[-1,0]]))-np.eye(2)).max():.1e}).")
P("    Die Halbgerade verlangt eine raumartige innere Richtung (elliptische Gleichung).")
P()

# ================= B1 Kerne =================
P("B1) Kerne auf dem Bahntorus, τ = x + iy, q = e^{2πiτ}, s = 2πy (y = Im τ_M = s/2π)")
P("    Szegő S = Σ_{n≥0}qⁿ = 1/(1 − q);  Poisson P = Σ_{n∈ℤ} e^{2πinx}e^{−2π|n|y} = (1 − |q|²)/|1 − q|²;  Schwarz H = (1 + q)/(1 − q)")
P("    s      x      Re S        ½(1 + P)     P           Re H        |S|²")
for s in (0.3,1.0,3.0):
    for x in (0.0,0.1,0.25):
        q=np.exp(2j*np.pi*x-s); S=1/(1-q); Pk=(1-abs(q)**2)/abs(1-q)**2; H=(1+q)/(1-q)
        P(f"    {s:3.1f}   {x:4.2f}   {S.real:10.6f}  {0.5*(1+Pk):10.6f}  {Pk:10.6f}  {H.real:10.6f}  {abs(S)**2:10.6f}")
P("    → Re S = ½(1 + P) und Re H = P überall; auf der Achse x = 0: S = 1 + n_B(s), H = P = coth(s/2) = 1 + 2n_B(s).")
P("      S und H sind beide holomorph und beide Reihen nur in q; sie unterscheiden sich im Gewicht (cₙ = 1 gegen cₙ = 2).")
P()

# ================= B2 Kern des Autors =================
P("B2) Kern des Autors Δ(τ) = q∏(1 − qⁿ)²⁴ = Σ_{n≥1} τ(n)qⁿ: Fourier-Koeffizienten längs Re τ bei Im τ = 0,3 (256 Punkte)")
M=256; xg2=np.arange(M)/M; y=0.3; qg=np.exp(2j*np.pi*(xg2+1j*y))
Dl=qg*np.prod([(1-qg**m)**24 for m in range(1,80)],axis=0)
cf=np.fft.fft(Dl)/M
tau_n={1:1,2:-24,3:252,4:-1472,5:4830}
P("    n:   " + "  ".join(f"{m}: {cf[m].real/np.exp(-2*np.pi*m*y):9.3f}" for m in range(1,6)) + "   (Soll τ(n) = 1, −24, 252, −1472, 4830)")
P(f"    Anteile mit n ≤ 0: max|c_n| = {max(abs(cf[0]),abs(cf[-40:]).max()):.1e};  normiert: c₁ = 1")
# Gittersumme G_4 über beide Richtungen der Rückkehr
tau0=0.13+0.8j; qs=np.exp(2j*np.pi*tau0); mm=np.arange(-3000,3001)
row=lambda nn: np.sum(1/(mm+nn*tau0)**4) if nn!=0 else np.sum(1/mm[mm!=0].astype(float)**4)
pos=sum(row(nn) for nn in range(1,40)); neg=sum(row(nn) for nn in range(-39,0))
sig3=lambda nn: sum(d**3 for d in range(1,nn+1) if nn%d==0)
qser=(2*np.pi*1j)**4/6*sum(sig3(nn)*qs**nn for nn in range(1,40))
P(f"    Gittersumme G₄: Zeilen n > 0: {pos:.10f};  Zeilen n < 0: {neg:.10f};  (2πi)⁴/3!·Σσ₃(n)qⁿ = {qser:.10f}")
P(f"    ganze Summe {row(0)+pos+neg:.10f} = 2ζ(4) + 2·(2πi)⁴/3!·Σσ₃qⁿ = {2*np.pi**4/90+2*qser:.10f}")
P("    → Rückkehr in beiden Richtungen (n und −n) gibt beide Male qⁿ, nie q̄ⁿ: Holomorphie wählt nicht die Richtung der Rückkehr,")
P("      die Gegenrichtung verdoppelt das Gewicht (Faktor 2 vor der q-Reihe).")
P()

# ================= B3 Cauchy-Riemann =================
P("B3) Paar σ + iσ̃ mit Cauchy-Riemann in (x, w): σ = Poisson-Fortsetzung des Orts, σ̃ = Hilbert-Transformierte (−i·sgn n)")
Mx=4096; ng=np.fft.fftfreq(Mx,1/Mx)
for s in (0.3,1.0,3.0):
    yy=s/(2*np.pi); damp=np.exp(-2*np.pi*np.abs(ng)*yy)
    sig=np.fft.ifft(damp).real*Mx; sigt=np.fft.ifft(-1j*np.sign(ng)*damp).real*Mx
    xg3=np.arange(Mx)/Mx; H=(1+np.exp(2j*np.pi*(xg3+1j*yy)))/(1-np.exp(2j*np.pi*(xg3+1j*yy)))
    P(f"    s = {s}: max|σ + iσ̃ − (1 + q)/(1 − q)| = {np.abs(sig+1j*sigt-H).max():.1e};  auf der Achse σ̃ = {sigt[0]:.1e}, σ + iσ̃ = {sig[0]:.6f} = coth(s/2) = {coth(s/2):.6f}")
P("    → Ein Cauchy-Riemann-Paar ist die holomorphe Ergänzung eines reellen harmonischen Feldes: positive Frequenzen mit Gewicht 2,")
P("      auf der Achse der Poisson-Wert 1 + 2n_B. Es gibt nicht 1 + n_B.")
P()

# ================= B4 Hardy, Orientierung =================
P("B4) Hardy-Projektion P₊ (nur n ≥ 0, Gewicht 1) des Orts: Fortsetzung = S; Re S = ½(1 + P).")
P("    Das ist dieselbe Aussage wie „jeder Rücklauf n ≥ 0 einmal“; eine eigene Begründung für n ≥ 0 liefert sie nicht.")
for s in (1.0,):
    q=np.exp(-s); S=1/(1-q); Sb=1/(1-np.conj(np.exp(2j*np.pi*0.1-s)))
    P(f"    Orientierung τ → −τ̄ (q → q̄): S(q̄) = conj S(q) (x = 0,1, s = 1: {1/(1-np.exp(2j*np.pi*0.1-s)):.5f} ↔ {Sb:.5f}); auf der Achse gleich ({S:.5f}).")
P("    → Im τ > 0 wählt q statt q̄, nicht das Gewicht; bei Re τ = 0 ist die Wahl unsichtbar.")
P()

# ================= B5 Gewicht je Rücklauf =================
P("B5) Gewicht je Rücklauf. Halbgerade ⇒ Erneuerung: G(mL) = G(L)^m = q^m (Halbgruppe). N gleichwertige Sinne, Rücklauflänge βL:")
P("    ν = 1 + N·Σ_{m≥1} e^{−βms} = 1 + N·n_B(βs),   χ(s) = s·(ν − 1) → N/β für s → 0,   Newton-Grenze ν − 1 → N·e^{−βs}")
P("    Fassung                       N  β   c₁(q = e^{−s})  χ(0)   v⁴/(GMa₀)   ν(a_N = a₀)   ν(0,1a₀)   ν(0,01a₀)")
for name,Nn,be,c1 in (("Szegő, ein Sinn",1,1,1),("Poisson, beide Sinne",2,1,2),("coth s (beide, Länge 2L)",2,2,0),("drei Sinne, Länge 3L",3,3,0)):
    nuv=lambda x: 1+Nn*nB(be*np.sqrt(x))
    P(f"    {name:28s} {Nn}  {be}   {c1:14d}  {Nn/be:4.1f}   {(Nn/be)**2:9.1f}   {nuv(1.0):10.4f}   {nuv(0.1):8.4f}   {nuv(0.01):8.3f}")
P("    Grenzwert des Autors χ → 1 (Tori zerfallen) ⇔ v⁴ = GMa₀ mit derselben Schwelle in s ⇔ N = β.")
P("    Rücklauflänge L = r* (β = 1) ⇒ N = 1: genau ein Sinn, ν = 1 + n_B(s). Dazu passt c₁ = 1, der erste Koeffizient des normierten Δ.")
P()

# ================= B6 Zahlen =================
P("B6) Zahlen Szegő gegen Poisson bei gleicher Schwelle, und gegen coth s (gleicher Tully-Fisher-Nullpunkt)")
P("    a_N/a₀    1 + n_B     1 + 2n_B    coth s     Poisson/Szegő   coth s/Szegő")
for x in (10,3,1,0.3,0.1,0.01,1e-4):
    s=np.sqrt(x)
    P(f"    {x:7g}   {1+nB(s):9.4f}   {1+2*nB(s):9.4f}   {coth(s):8.4f}   {(1+2*nB(s))/(1+nB(s)):10.4f}   {coth(s)/(1+nB(s)):10.4f}")
Mb=6e10*Msun; rs=np.sqrt(G*Mb/a0)
P(f"    Milchstraße als Punktmasse 6·10¹⁰ M☉, r* = {rs/kpc:.2f} kpc, v in km/s:")
for rr in (0.5,1,3,10):
    r=rr*rs; aN=G*Mb/r**2; s=np.sqrt(aN/a0)
    v=lambda nu_: np.sqrt(nu_*aN*r)/1e3
    P(f"      r = {rr:4.1f} r*: Szegő {v(1+nB(s)):6.1f}, Poisson {v(1+2*nB(s)):6.1f}, coth s {v(coth(s)):6.1f}")
P("    Newton-Grenze (Sonne bei 1 AE, s = 7031): ν − 1 = e^{−s} bzw. 2e^{−s}, beide unter 10⁻³⁰⁰⁰.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"holomorphie_rueckkehr.txt"),"w",encoding="utf-8").write(txt); print(txt)
