# Halbgerade als Spitze: Ist die offene innere Richtung (D·G'' = Γ·G, G(0) = 1, G beschränkt ⇒ q = e^{−s})
# die Richtung Im τ zur Spitze des Modulraums?
# 1 Zuordnung Mode λ = 0 ↔ Halbgerade (beide Wahlen von Im τ_M), Ort und Spitze, Stabilisator von ∞
# 2 Umlaufsinn und Gewicht (Szegő, Schwarz/Poisson, Bergman), Phase Re τ, Maass-Formen (λ ≠ 0), Spektrum von SL(2,ℤ) (Hejhal)
# 3 Geschlossene Schleife der Periode Λ in Im τ   4 prüfbare Folgen
import numpy as np, sys, io, os
from scipy.special import kv, gamma as Gam
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; AU=1.496e11; Rsun=6.957e8
out=[]; P=lambda s="": out.append(s)
nu_open=lambda s: 1/(-np.expm1(-s))
TP=2*np.pi

def Kir(R, x):
    """K_{iR}(x)·e^{πR/2} für x > 0 (Trapezregel des geraden Integranden, spektral genau)."""
    x=np.atleast_1d(np.asarray(x,float)); res=np.empty_like(x)
    for i,xi in enumerate(x):
        T=np.arccosh(max(60.0/xi,1.0))+1.0; h=0.01; t=np.arange(0,T+h,h)
        f=np.exp(-xi*np.cosh(t)+np.pi*R/2)*np.cos(R*t); res[i]=h*(f.sum()-0.5*f[0])
    return res

def d2(f, y, h=1e-4): return (f(y+h)-2*f(y)+f(y-h))/h**2

# ---------------- 1) Zuordnung ----------------
P("1a) Laplace-Beltrami y²(∂_x² + ∂_y²)f = λf, f = e^{2πinx}g(y)  ⇒  y²(g'' − 4π²n²g) = λg")
for n in (1,2,3):
    for sg in (+1,-1):
        g=lambda y,n=n,sg=sg: np.exp(sg*TP*n*y)
        y=0.37; P(f"    λ = 0, n = {n}, g = e^{{{'+' if sg>0 else '−'}2π{n}y}}: Residuum g'' − 4π²n²g = {abs(d2(g,y)-(TP*n)**2*g(y))/g(y):.1e}·g")
P("    beschränkt zur Spitze (y → ∞): nur e^{−2π|n|y}; n = 0: g = A + By ⇒ Konstante.")
x,y,h=0.13,0.29,1e-6
for n,name in ((1,"q = e^{2πiτ}"),(-1,"q̄ = e^{−2πiτ̄}")):
    f=lambda x,y,n=n: np.exp(2j*np.pi*n*x)*np.exp(-TP*abs(n)*y)
    fx=(f(x+h,y)-f(x-h,y))/(2*h); fy=(f(x,y+h)-f(x,y-h))/(2*h)
    lap=(f(x+1e-4,y)+f(x-1e-4,y)+f(x,y+1e-4)+f(x,y-1e-4)-4*f(x,y))/1e-8
    P(f"    n = {n:+d}: {name}: |Δf|/|f| = {abs(lap)/abs(f(x,y)):.1e}, |∂_τ̄ f| = {abs(0.5*(fx+1j*fy)):.1e}, |∂_τ f| = {abs(0.5*(fx-1j*fy)):.2e}, |f| = e^{{−2πy}} beschränkt")
P("    → Beschränktheit an der Spitze lässt {1, qⁿ, q̄ⁿ} zu; bei y = 0 ist |qⁿ| = |q̄ⁿ| = 1 für jedes n (G(0) = 1 von selbst).")
# Randwertaufgabe numerisch: g'' = (2π)²g auf [0, W], g(0) = 1, g(W) = 0
W=4.0; N=4000; yy=np.linspace(0,W,N+1); hh=yy[1]-yy[0]; k2=TP**2
A=np.diag(np.full(N-1,-2/hh**2-k2))+np.diag(np.full(N-2,1/hh**2),1)+np.diag(np.full(N-2,1/hh**2),-1)
b=np.zeros(N-1); b[0]=-1/hh**2; gg=np.linalg.solve(A,b)
for yv in (0.05,0.2,0.5):
    i=int(round(yv/hh)); P(f"    Randwertaufgabe g(0) = 1, g(4) = 0: g({yv}) = {gg[i-1]:.6f},  e^{{−2π·{yv}}} = {np.exp(-TP*yv):.6f}")
P()

P("1b) Zuordnung w = ℓ·Im τ: D·G'' = Γ·G wird g'' = (2π)²g für ℓ = 2π√(D/Γ); dann L√(Γ/D) = 2π·Im τ_M und G(L) = |q(τ_M)|")
P("    Wahl 1: Im τ_M = s/2π (= r*/2πr bei der Punktmasse), √(Γ/D) = 1/r ⇒ ℓ = 2πr, w = 2πr·Im τ; G(L) = e^{−s} = q (Rückkehranteil)")
P("      a_N/a₀     s      Im τ_M    q = e^{−2πIm τ}   im Spitzenbereich Im τ > 1?   S-Bild Im(−1/τ) = 2π/s")
for xN in (1e3, 4*np.pi**2, 10, 1, 0.1, 0.01):
    s=np.sqrt(xN); yM=s/TP
    P(f"      {xN:8.3g}  {s:6.3f}  {yM:8.4f}   {np.exp(-TP*yM):.4e}        {'ja' if yM>1 else 'nein'}                    {TP/s:8.3f}")
P(f"      Grenze Im τ_M = 1 ⇔ s = 2π ⇔ a_N = 4π²a₀ = {4*np.pi**2:.2f} a₀ = {4*np.pi**2*a0:.2e} m/s². Galaktische Punkte (s < 2π) liegen unter der Spitze.")
P("      Im τ = (w/L)·s/2π: Ort w = 0 ↔ Im τ = 0 (|q| = 1, innerer Zyklus 0); Spitze w → ∞ ↔ q → 0. Bei w = L dieselbe Gerade in a_N: q → 1 für a_N → 0, q → 0 für a_N → ∞ (Newton).")
P("    Wahl 2: Im τ_M = ln(r/Δs)/2π ⇒ G(L) = Δs/r = |q_Δs| (Kern-Nome, nicht der Rückkehranteil); L√(Γ/D) = ln(r/Δs) statt s")
def ds(M): return G*M/(3*c**2)
for name,M,r in (("Erde (Sonne, 1 AE)",Msun,AU),("Milchstraße 6·10¹⁰ M☉, 8 kpc",6e10*Msun,8*kpc),("Milchstraße, 100 kpc",6e10*Msun,100*kpc)):
    s=np.sqrt(G*M/r**2/a0); yM=np.log(r/ds(M))/TP
    P(f"      {name:30s}: Δs = {ds(M):.3e} m, Im τ_M = {yM:.3f}, Δs/r = {ds(M)/r:.2e}, s = {s:.4g}, ln(r/Δs)/s = {np.log(r/ds(M))/s:.3g}")
P(f"      Ort w = 0 ↔ r = Δs (|q_Δs| = 1); Spitze ↔ r → ∞; Spitzenbereich Im τ > 1 ⇔ r > e^{{2π}}Δs = {np.exp(TP):.0f} Δs.")
P()

P("1c) Modulgruppe und Im τ")
Sfix=[]; Stab=[]
R=range(-6,7)
for a in R:
    for b_ in R:
        for c_ in R:
            for d in R:
                if a*d-b_*c_!=1: continue
                if c_==0: Stab.append((a,b_,c_,d))
                if (b_==0 and c_==0) or (a==0 and d==0): Sfix.append((a,b_,c_,d))
P(f"    Stabilisator von ∞ (|Einträge| ≤ 6): {len(Stab)} Elemente, alle ±(1 b; 0 1): {all(abs(m[0])==1 and m[0]==m[3] for m in Stab)} → nur τ ↦ τ + b, keine Streckung in Im τ.")
P(f"    Elemente, die die Achse Re τ = 0 in sich abbilden: {sorted(set(Sfix))} → ±I und ±S (y ↦ 1/y): die Achse wird bei τ = i gefaltet, nicht geschlossen.")
P(f"    |q| ist nicht SL(2,ℤ)-invariant: τ = i/2 und −1/τ = 2i geben |q| = {np.exp(-np.pi):.4f} und {np.exp(-4*np.pi):.2e}, ν = {nu_open(np.pi):.4f} und {nu_open(4*np.pi):.7f}.")
P("    → ν lebt auf ℍ/⟨T⟩ (markierter Torus: Umlauf und Rückkehr unterschieden) = punktierte Kreisscheibe 0 < |q| < 1, Im τ ∈ (0, ∞).")
P()

# ---------------- 2) Umlaufsinn, Gewicht, Maass ----------------
P("2a) Holomorphe Funktionen mit Pol am Ort (q = 1) und ν(q = 0) = 1, auf der Achse Re τ = 0:")
P("      a_N/a₀   Szegő 1/(1−q)   Schwarz (1+q)/(1−q) = Poisson   Bergman 1/(1−q)²   (1+q²)/(1−q²) = coth s")
for xN in (3,1,0.3,0.1,0.01):
    s=np.sqrt(xN); q=np.exp(-s)
    P(f"      {xN:6g}   {1/(1-q):10.4f}   {(1+q)/(1-q):14.4f}              {1/(1-q)**2:12.4f}       {(1+q*q)/(1-q*q):10.4f}")
s=1e-4; q=np.exp(-s)
P(f"      s → 0: s·ν = {s/(1-q):.4f} | {s*(1+q)/(1-q):.4f} | s²ν = {s*s/(1-q)**2:.4f} (a → a₀) | {s*(1+q*q)/(1-q*q):.4f}")
P("      Koeffizienten cₙ von qⁿ: Szegő 1,1,1,… ; Schwarz 1,2,2,… ; Bergman n+1 ; coth s: 1,0,2,0,2,…")
P("    Randwerte bei |q| = 1: Szegő = Projektion von δ(Re τ) auf n ≥ 0 (Gegenläufe entfallen); Schwarz: Realteil = δ (Gegenläufe auf n > 0 geklappt).")
P("    G(0) = 1 für jede gezählte Mode heißt cₙ = 1: gezählt n ≥ 0 → Szegő = 1 + n_B; gezählt n ∈ ℤ (qⁿ und q̄ⁿ, beide beschränkt) → Poisson = 1 + 2n_B.")
P("    Bergman (cₙ = n+1) und coth s (nur gerade n, cₙ = 2) verletzen G(0) = 1 je Mode. Ob n < 0 mitzählt, entscheidet die Spitze nicht.")
P("    Die Funktion der Spitze hängt von Re τ ab; 1 + n_B verlangt Re τ = 0 (innerer Zyklus senkrecht zur Bahn):")
x=np.sqrt(2)-1
for s in (0.01,1.0):
    qq=np.exp(TP*1j*(x+1j*s/TP)); f=1/(1-qq)
    P(f"    Phase: Re τ = √2 − 1 (κ/Ω = √2), s = {s}: Re[1/(1−q)] = {f.real:.4f}, |1/(1−q)| = {abs(f):.4f}  (Re τ = 0: {nu_open(s):.4f})")
P()

P("2b) Maass-Moden λ ≠ 0: g = √y·K_ν(2πny), ν = √(¼ + λ); Prüfung der Gleichung y²(g'' − 4π²g) = λg:")
for lam in (-0.2,0.0,0.3,2.0):
    v=np.sqrt(0.25+lam); g=lambda y,v=v: np.sqrt(y)*kv(v,TP*y); y=0.4
    P(f"      λ = {lam:+.2f}: Residuum/|g| = {abs(y*y*(d2(g,y)-TP**2*g(y))-lam*g(y))/abs(g(y)):.1e}")
P("    Am Ort (y → 0) g ∝ y^{½−ν}(1 + … + C·y^{2ν}): endlicher Wert ≠ 0 nur für ν = ½ ⇔ λ = 0.")
P("    Normiert auf das führende Glied: q_λ(s) = (2/Γ(ν))(s/2)^ν K_ν(s)  (Wahl 1, s = 2π Im τ); λ = 0 gibt e^{−s}:")
P(f"      Prüfung ν = ½, s = 0,7: {(2/Gam(0.5))*(0.35)**0.5*kv(0.5,0.7):.6f} = e^{{−0,7}} = {np.exp(-0.7):.6f}")
P("      λ        ν      1 − q ∝ s^p (p bei s = 10⁻³)   d ln v/d ln r (Punktmasse, s = 10⁻³)   Grenzwert ν − ½")
def qlam(s,v): return (2/Gam(v))*(s/2)**v*kv(v,s)
for lam in (-0.2,-0.1,-0.0475,0.0,0.0525,0.1,0.5,2.0):
    v=np.sqrt(0.25+lam); s1,s2=1e-3,1.001e-3
    p=(np.log(1-qlam(s2,v))-np.log(1-qlam(s1,v)))/(np.log(s2)-np.log(s1))
    L1=np.log(s1/(1-qlam(s1,v))); L2=np.log(s2/(1-qlam(s2,v)))
    slope=-0.5*(L2-L1)/(np.log(s2)-np.log(s1))
    P(f"      {lam:+7.4f}  {v:6.4f}   {p:8.4f}                       {slope:+8.4f}                          {min(v,1)-0.5:+.4f}")
P("    a ∝ a_N^{1−ν}a₀^ν: flach und v⁴ ∝ M nur für ν = ½. |d ln v/d ln r| ≤ 0,05 ⇔ −0,0475 ≤ λ ≤ 0,0525; ≤ 0,02 ⇔ −0,0196 ≤ λ ≤ 0,0204.")
P()

P("2c) Spektrum von SL(2,ℤ): Maass-Spitzenformen, Hejhal-Verfahren (Hecke-Probe a₂a₃ = a₆)")
def pullback(x,y):
    while True:
        x=x-np.floor(x+0.5); r2=x*x+y*y
        if r2<1-1e-13: x,y=-x/r2,y/r2
        else: return x,y
def coeffs(Rr,par,M=20,Q=26,Y=0.42):
    tr=np.cos if par==0 else np.sin
    xm=(np.arange(1,Q+1)-0.5)/(2*Q); pb=np.array([pullback(x,Y) for x in xm]); xs,ys=pb[:,0],pb[:,1]
    n=np.arange(1,M+1); KY=np.array([Kir(Rr,TP*k*Y)[0] for k in n]); Ks=np.array([Kir(Rr,TP*l*ys) for l in n])
    Wm=(2/Q)*np.einsum('lm,km->kl',np.sqrt(ys)*Ks*tr(TP*np.outer(n,xs)),tr(TP*np.outer(n,xm)))
    V=np.diag(np.sqrt(Y)*KY)-Wm; a=np.linalg.solve(V[1:,1:],-V[1:,0]); return np.concatenate(([1.0],a))
def hh_(Rr,par): return coeffs(Rr,par,Y=0.40)[1]-coeffs(Rr,par,Y=0.46)[1]
found=[]
for par in (1,0):
    Rs=np.arange(8.0,14.5,0.1); hs=[hh_(Rr,par) for Rr in Rs]
    for i in range(len(Rs)-1):
        if np.sign(hs[i])!=np.sign(hs[i+1]):
            Rr=brentq(lambda t: hh_(t,par),Rs[i],Rs[i+1],xtol=1e-12); cf=coeffs(Rr,par)
            if abs(cf[1]*cf[2]-cf[5])<1e-6: found.append((Rr,par,cf[1]))
found.sort()
for Rr,par,a2 in found:
    P(f"      t = {Rr:.10f} ({'ungerade' if par else 'gerade'}), a₂ = {a2:+.6f}, λ = −(¼ + t²) = {-(0.25+Rr**2):.3f}")
t1=found[0][0]
P("    Außerdem λ = 0 (Konstante) und das kontinuierliche Spektrum λ ≤ −¼. Im Tully-Fisher-Band liegt nur λ = 0.")
P(f"    λ ≤ −¼: ν = it imaginär, g = √y·K_it(2πy) → 0 wie √y am Ort, schwingt für y < t/2π (t₁: y < {t1/TP:.3f}, Wahl 1: a_N < {t1**2:.1f} a₀).")
xs_=np.exp(np.linspace(np.log(1e-3),np.log(t1),60001)); Kv=Kir(t1,xs_); z=xs_[1:][np.sign(Kv[1:])!=np.sign(Kv[:-1])][::-1]
P(f"    Nullstellen von q ∝ √s·K_it₁(s) (Wahl 1): s = {', '.join(f'{v:.4f}' for v in z[:5])}, … ({len(z)} Nullstellen über s = 10⁻³); tief: Verhältnis {z[-2]/z[-1]:.4f} = e^{{π/t₁}} = {np.exp(np.pi/t1):.4f}")
P("    → q wechselt log-periodisch das Vorzeichen und fällt wie √s: a → a_N (Newton) im tiefen Bereich, keine flache Kurve.")
P()

# ---------------- 3) Geschlossene Schleife in Im τ ----------------
P("3a) Im τ geschlossen mit Periode Λ (flach, y ~ y + Λ; nur für λ = 0 zulässig, eine hyperbolische Streckung verträgt sich nicht mit Re τ mod 1)")
k=3; T_=np.array([[1,1],[0,1]]); D_=np.array([[k,0],[0,1]])
C_=np.linalg.inv(D_)@T_@D_
P(f"    Streckung τ ↦ {k}τ mit T: D⁻¹TD = τ ↦ τ + {C_[0,1]/C_[1,1]:.4f}; Verschiebungen um 1/kᵐ ⇒ Gruppe nicht diskret, Re τ mod 1 ginge verloren.")
P("    Mode n mit Wert 1 am Ort (y = 0 ≡ Λ): Gₙ(y) = cosh(2πn(y − Λ/2))/cosh(πnΛ)")
Lam=1.3; n=2; g=lambda y: np.cosh(TP*n*(y-Lam/2))/np.cosh(np.pi*n*Lam); y=0.41
P(f"      Prüfung Λ = {Lam}, n = {n}: Residuum {abs(d2(g,y)-(TP*n)**2*g(y))/g(y):.1e}, G(0) = {g(0):.6f}, G(Λ) = {g(Lam):.6f}, G(y) − G(Λ−y) = {g(y)-g(Lam-y):.1e}")
P(f"      gerade um den Ort ⇒ beide Umlaufsinne; qⁿ ist nicht periodisch (|q(y+Λ)/q(y)| = e^{{−2πΛ}} = {np.exp(-TP*Lam):.2e}); ohne Rand keine holomorphe Funktion außer Konstanten,")
P("      und eine Funktion mit nur einem einfachen Pol am Ort gibt es auf dem Torus nicht (Residuensumme 0). Bei Λ = 2·Im τ_M: G₁ = 1/cosh s.")
P("    Wahl 1: q_Λ(s) = G₁(s/2π) = cosh(s − πΛ)/cosh(πΛ); tief: 1 − q ≈ s·tanh(πΛ) ⇒ Tully-Fisher bleibt, v⁴ = GMa₀/tanh²(πΛ); q = 1 bei s = 2πΛ")
P("      Λ      v⁴-Faktor 1/tanh²(πΛ)   10 % Abweichung in a: q-Reihe ab s   Modensumme ΣGₙ ab s   r/r* (q-Reihe)   Pol bei s = 2πΛ")
def nuA(s,L): q=np.cosh(s-np.pi*L)/np.cosh(np.pi*L); return 1/(1-q)
def nuB(s,L):
    nn=np.arange(1,20000); A=np.abs(nn*(s-np.pi*L)); B=nn*np.pi*L
    return 1+np.sum(np.exp(A-B)*(1+np.exp(-2*A))/(1+np.exp(-2*B)))
for L in (0.5,1,2,5):
    ss=np.linspace(1e-3,TP*L-1e-3,20000)
    rA=np.array([nuA(s,L)/nu_open(s) for s in ss]); rB=np.array([nuB(s,L)/nu_open(s) for s in ss[::20]])
    iA=np.argmax(np.abs(rA-1)>0.1); iB=np.argmax(np.abs(rB-1)>0.1)
    sA=ss[iA] if np.abs(rA[iA]-1)>0.1 else np.nan; sB=ss[::20][iB] if np.abs(rB[iB]-1)>0.1 else np.nan
    P(f"      {L:4g}   {1/np.tanh(np.pi*L)**2:10.5f}                 {sA:8.3f}                          {sB:8.3f}            {1/sA:7.3f}          {TP*L:8.3f}")
sS=np.sqrt(G*Msun/Rsun**2/a0); sE=np.sqrt(G*Msun/AU**2/a0)
P(f"      Abweichung liegt innen (Newton-Bereich), nicht außen. Newton an der Sonnenoberfläche (s = {sS:.3g}) verlangt Λ > {sS/TP:.3g}; Erdbahn s = {sE:.0f} ⇒ Λ > {sE/TP:.0f}.")
P(f"      Dann ist die Schleife in Galaxien unsichtbar (Abweichung ~ e^{{−2πΛ}}).")
P("    Wahl 2: Kern-Nome Δs/r → G₁(ln(r/Δs)/2π) = (Δs/r)·(1 + (r/Δs)²e^{−2πΛ})/(1 + e^{−2πΛ}): wächst außen wieder an")
ydev=lambda L: brentq(lambda y: (1+np.exp(TP*2*y-TP*L))/(1+np.exp(-TP*L))-1.1,0,L)
for L in (5,6,8):
    P(f"      Λ = {L}: +10 % ab Im τ = {ydev(L):.4f} (= Λ/2 − {L/2-ydev(L):.4f}), r = {np.exp(TP*ydev(L)):.3e}·Δs")
yMW=np.log(100*kpc/ds(6e10*Msun))/TP
P(f"      Flach bis r_max verlangt Λ ≥ 2·Im τ(r_max) + {2*(5/2-ydev(5)):.3f}; Beispiel Milchstraße bis 100 kpc: Im τ = {yMW:.3f} ⇒ Λ ≥ {2*yMW+2*(5/2-ydev(5)):.2f}. Abweichungsradius ∝ Δs ∝ M.")
P("    Feste innere Länge Λ_w (Meter) ist keine Schleife in Im τ: dort wäre die Periode Λ_w/(2πr) ∝ 1/r, also von System und Radius abhängig.")
P()

# ---------------- 4) Prüfbare Folgen ----------------
P("4)  Offene Spitze (Λ = ∞, λ = 0): genau a = a_N/(1 − e^{−s}), keine neue Form.")
P("    Beimischung λ ≠ 0: im tiefen Bereich massenunabhängige Restneigung d ln v/d ln r = √(¼ + λ) − ½.")
P("    Geschlossene Schleife: Abweichung von der flachen Kurve bei Wahl 1 innen bei festem a_N ≈ (2πΛ)²a₀ (durch Newton ausgeschlossen),")
P("    bei Wahl 2 außen bei r ∝ M, bei fester innerer Länge bei r ∝ r* ∝ √M. Ein gefundener Knick trennt die Fassungen über seine Massenskalierung.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"halbgerade_spitze.txt"),"w",encoding="utf-8").write(txt); print(txt)
