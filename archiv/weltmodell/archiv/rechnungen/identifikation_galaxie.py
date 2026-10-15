# Welche Identifikation trägt 𝒫 im Zustandsraum einer Galaxie?
# A: reale Richtungen des Bahntorus (Winkel, Drehungen)   B: Rückkehr als Translation mit Gewicht |q| (Abel-Mittel)
# C: einseitig (kausal) gegen beidseitig (volle Gruppe)   D: Gittersummen = Eisenstein-Reihen   E: Zahlen Milchstraße
import numpy as np, sys, io, os
from math import factorial, pi
from scipy.special import zeta
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19
M=6e10*Msun; GM=G*M; rs=np.sqrt(GM/a0)
out=[]; P=lambda s="": out.append(s)
rng=np.random.default_rng(20260929)

# ---------------- A ----------------
P("A) Drehbilder einer Balkenverteilung (200 Punktmassen, |x| ≤ 3 kpc, σ_y = 0,5 kpc), Feldpunkt r = 20 kpc, 720 Winkel")
N=200; x=rng.uniform(-3,3,N)*kpc; y=rng.normal(0,0.5,N)*kpc; m=np.full(N,GM/N)
ph=np.linspace(0,2*np.pi,720,endpoint=False); R=20*kpc; X=R*np.cos(ph); Y=R*np.sin(ph)
def ar(xs,ys,ms):
    dx=X[:,None]-xs[None,:]; dy=Y[:,None]-ys[None,:]; d3=(dx*dx+dy*dy)**1.5
    ax=-(ms*dx/d3).sum(1); ay=-(ms*dy/d3).sum(1); return ax*np.cos(ph)+ay*np.sin(ph)
a_orig=ar(x,y,m); m2=lambda a: abs(np.fft.rfft(a)[2])/len(a)*2
P(f"   Original: ⟨a_r⟩ = {a_orig.mean():.6e} m/s²,  m=2-Anteil {m2(a_orig)/abs(a_orig.mean()):.3e}")
for n in (2,3,4,6):
    xs=np.concatenate([x*np.cos(2*np.pi*j/n)-y*np.sin(2*np.pi*j/n) for j in range(n)])
    ys=np.concatenate([x*np.sin(2*np.pi*j/n)+y*np.cos(2*np.pi*j/n) for j in range(n)])
    a=ar(xs,ys,np.tile(m/n,n))
    P(f"   n = {n}:  ⟨a_r⟩ relativ geändert {abs(a.mean()/a_orig.mean()-1):.1e},  m=2-Anteil {m2(a)/abs(a.mean()):.3e}")
P("   → Drehidentifikation löscht nur nichtaxiale Anteile (Fluchtmoden des Winkels); der bindende Monopol bleibt gleich.")
P()

# ---------------- B ----------------
P("B) Rückkehr als Translation T (quasiperiodisch: Drehung um den goldenen Winkel, Moden |k| ≤ 2000)")
P("   Abel-Mittel 𝒫_q = (1−q)Σ_{m≥0} qᵐTᵐ = (1−q)(I − qT)⁻¹;  𝒫 = Projektor auf T-invariante Funktionen (Satz 1)")
alpha=(np.sqrt(5)-1)/2; k=np.arange(1,2001); lam=np.exp(2j*np.pi*k*alpha)
P("       q      ‖𝒫_q − 𝒫‖    Verstärkung invariant 1/(1−q)   max Fluchtmoden |1/(1−qλ)|")
for q in (0.0,0.5,0.9,0.99,0.999,0.9999):
    P(f"   {q:8.4f}   {np.max(np.abs((1-q)/(1-q*lam))):9.3e}   {1/(1-q):14.1f}              {np.max(np.abs(1/(1-q*lam))):9.1f}")
P("   q → 0: 𝒫_q = I (offen, Newton);  q → 1: 𝒫_q → 𝒫 (voll geschlossen). Die Rückkopplungssumme (I − qT)⁻¹")
P("   verstärkt genau den wiederkehrenden Anteil um 1/(1−q); Fluchtmoden bleiben beschränkt (Kleinstabstand 1/min|1−λ|).")
P(f"   Kleinstes |1−λ_k| im Fenster: {np.min(np.abs(1-lam)):.3e}  → Schranke {1/np.min(np.abs(1-lam)):.1f}")
P()

# ---------------- C ----------------
P("C) Kausale Halbgruppe (m ≥ 0) gegen volle Gruppe (m ∈ ℤ), gleiche Tully-Fisher-Grenze")
P("   einseitig ν₁ = 1/(1 − e^{−s});   beidseitig Σ_{m∈ℤ} e^{−2s|m|} = coth s   (beide → 1/s für s → 0, → 1 für s → ∞)")
P("     a_N/a₀      s        ν₁        coth s    Abweichung")
for x_ in (100,10,3,1,0.3,0.1,0.01):
    s=np.sqrt(x_); n1=1/(-np.expm1(-s)); n2=1/np.tanh(s)
    P(f"   {x_:8.2f}  {s:7.3f}  {n1:9.4f}  {n2:9.4f}   {n2/n1-1:+7.1%}")
P("   → nur die kausale Hälfte gibt die gemessene Form; die volle Gruppe liegt bei a_N = a₀ um 17 % tiefer.")
P()

# ---------------- D ----------------
P("D) Gittersumme über ℤ + τℤ der Potenzkerne z^{−k} (Bildsumme von 𝒫) gegen q-Reihe (Eisenstein)")
def lattice(tau,k,Rmax=400.0):
    t=complex(tau); S=0j; Nn=int(Rmax/t.imag)+1
    for n in range(-Nn,Nn+1):
        yi=n*t.imag
        if abs(yi)>Rmax: continue
        w=np.sqrt(Rmax**2-yi**2); m0=int(np.floor(-w-n*t.real)); m1=int(np.ceil(w-n*t.real))
        mm=np.arange(m0,m1+1); z=mm+n*t; z=z[np.abs(z)<=Rmax]; z=z[z!=0]; S+=np.sum(z**(-k))
    return S
def qseries(tau,k,Nmax=400):
    q=np.exp(2j*np.pi*tau); n=np.arange(1,Nmax+1)
    sig=np.array([sum(d**(k-1) for d in range(1,j+1) if j%d==0) for j in n],float)
    return 2*zeta(k)+2*(2j*np.pi)**k/factorial(k-1)*np.sum(sig*q**n)
def lipschitz_direct(tau,k,n,Mmax=200000):
    mm=np.arange(-Mmax,Mmax+1); return np.sum((mm+n*tau)**(-k))
def lipschitz_q(tau,k,n,Dmax=400):
    d=np.arange(1,Dmax+1); return (-2j*np.pi)**k/factorial(k-1)*np.sum(d**(k-1.0)*np.exp(2j*np.pi*d*n*tau))
for tau,name in ((0.3+1.1j,"τ = 0,3 + 1,1i"),(1j/(2*pi),"τ = i/2π (s = 1, |q| = e⁻¹)")):
    P(f"   {name}")
    for kk in (4,6):
        L=lattice(tau,kk); Q=qseries(tau,kk)
        P(f"     G_{kk}: Gitter {L.real:+.8e}   q-Reihe {Q.real:+.8e}   rel. {abs(L/Q-1):.1e}")
tau=1j/(2*pi)
P("   Summe längs des Umlaufs (m) je Rückkehrschritt n bei τ = i/2π, k = 4: jede Zeile wird zur Reihe in qⁿ")
for n in (1,2,3):
    a=lipschitz_direct(tau,4,n); b=lipschitz_q(tau,4,n)
    P(f"     n = {n}:  Σ_m (m+nτ)⁻⁴ = {a.real:+.10e}   (2π)⁴/3!·Σ d³ q^{{nd}} = {b.real:+.10e}")
P("   Die niedrigste Reihe (d = 1) über alle Rückkehrschritte ist Σ_{n≥1} qⁿ = q/(1−q) = a/a_N − 1 der Kreisbahn.")
P()

# ---------------- E ----------------
P(f"E) Milchstraße baryonisch 6·10¹⁰ M☉, r* = {rs/kpc:.2f} kpc; Punktmasse s = r*/r, |q| = e^{{−s}}, Im τ = s/2π")
P("     r/r*     s       |q|     Im τ      ν = 1/(1−|q|)   r/r* + ½")
for f in (0.1,0.3,1,3,10):
    s=1/f; q=np.exp(-s); P(f"   {f:6.1f}  {s:6.2f}  {q:8.2e}  {s/(2*pi):6.3f}   {1/(1-q):10.4f}     {f+0.5:8.3f}")
P("   ν ist die mittlere Zahl der Durchläufe (geometrisch); im Fernfeld ν = r/r* + ½ + O(r*/r).")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"identifikation_galaxie.txt"),"w",encoding="utf-8").write(txt); print(txt)
