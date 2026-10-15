# Weitere Formeln Ramanujans, numerisch geprüft und am Modell gemessen
# A Ramanujan-Summation der Rückläufe   B Lambert-Reihe über alle Harmonischen   C bilaterale 1ψ1-Summe
# D Ramanujan-Summen als Perioden-Projektoren   E Basen 2, 3, 4, 6 und Bilderzahl der Zellgitter
import numpy as np, sys, io, os
from mpmath import mp, zeta, qp, mpf, nsum, inf
from math import gcd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
mp.dps = 30
out=[]; P=lambda s="": out.append(s)

# ---------------- A ----------------
P("A) ν(s) = Σ_{n≥0} e^{−ns} = 1/s + (1 + ζ(0)) + Σ_{k≥1} ζ(−k)(−s)^k/k!   (Ramanujan-Summen Σ1, Σn, Σn³ …)")
for s in (mpf('0.1'),mpf('0.5'),mpf(1)):
    exact=1/(1-mp.e**(-s))
    ser=1/s+1+zeta(0)+sum(zeta(-k)*(-s)**k/mp.factorial(k) for k in range(1,30))
    P(f"   s = {float(s):4.1f}:  exakt {float(exact):.15f}   Reihe {float(ser):.15f}")
P(f"   Glieder: 1 + ζ(0) = {float(1+zeta(0))},  −ζ(−1) = {float(-zeta(-1)):.6f} = 1/12,  ζ(−3)/3!·(−1) = {float(-zeta(-3)/6):.3e} = −1/720")
P("   M_schein/M = r/r* + ½ + r*/(12r) − (r*/r)³/720 + …: das erste Korrekturglied ist Ramanujans 1 + 2 + 3 + … = −1/12.")
P()

# ---------------- B ----------------
P("B) Lambert-Reihe über alle Harmonischen: L(s) = Σσ₁(n)e^{−ns} = (1 − E₂)/24 bei q = e^{−s}")
def sigma1(n): return sum(d for d in range(1,n+1) if n%d==0)
for s in (0.3,0.1,0.03):
    N=int(60/s); L=sum(sigma1(n)*np.exp(-n*s) for n in range(1,N))
    P(f"   s = {s:5.2f}:  L = {L:12.4f}   π²/(6s²) − 1/(2s) + 1/24 = {np.pi**2/(6*s*s)-1/(2*s)+1/24:12.4f}   a/a₀ = s²(1+L) = {s*s*(1+L):.4f}")
P(f"   → alle Harmonischen gleich gewichtet: a → (π²/6)·a₀ = {np.pi**2/6:.4f}·a₀, konstant, ohne Tully-Fisher; nur die Grundschwingung (Σqⁿ) gibt a → √(a_N a₀).")
P()

# ---------------- C ----------------
P("C) Ramanujans 1ψ1: Σ_{n∈ℤ} (a;q)_n/(b;q)_n zⁿ = (q, b/a, az, q/(az); q)_∞ / (b, q/a, z, b/(az); q)_∞")
def poch(x,q,n):
    if n>=0: return qp(x,q,n)
    return 1/qp(x*q**n,q,-n)
def invpoch(x,q,n):
    if n>=0: return 1/qp(x,q,n)
    return qp(x*q**n,q,-n)
def lhs(a,b,q,z,N=400):
    return sum(poch(a,q,n)*invpoch(b,q,n)*z**n for n in range(-N,N+1))
def rhs(a,b,q,z): return qp(q,q)*qp(b/a,q)*qp(a*z,q)*qp(q/(a*z),q)/(qp(b,q)*qp(q/a,q)*qp(z,q)*qp(b/(a*z),q))
a,b,q,z=mpf('0.9'),mpf('0.3'),mpf('0.6'),mpf('0.8')
P(f"   a = 0,9, b = 0,3, q = 0,6, z = 0,8 (|b/a| < |z| < 1):  Summe {float(lhs(a,b,q,z)):.12f}   Produkt {float(rhs(a,b,q,z)):.12f}")
S_neg=sum(poch(a,q,n)*invpoch(q,q,n)*z**n for n in range(-5,0))
P(f"   b = q: Glieder n < 0 haben den Faktor 1/(q;q)_n mit (1 − q⁰) = 0 im Nenner der Umkehrung → Summe der n < 0: {float(S_neg):.1e}")
P("   Mit a = q wird die einseitige Reihe zu Σ_{n≥0} zⁿ = 1/(1 − z); beidseitig verlangt |b/a| < |z| < 1 und endet für b = q von selbst bei n = 0.")
P()

# ---------------- D ----------------
P("D) Ramanujan-Summen c_q(n) = Σ_{(a,q)=1} e^{2πian/q}: Projektoren auf genaue Periode q")
Nn=12; T=np.roll(np.eye(Nn),1,axis=0)
def cq(q,n): return sum(np.exp(2j*np.pi*a*n/q) for a in range(1,q+1) if gcd(a,q)==1).real
E={}
for qq in [d for d in range(1,Nn+1) if Nn%d==0]:
    E[qq]=sum(cq(qq,n)*np.linalg.matrix_power(T,n) for n in range(Nn))/Nn
Ssum=sum(E.values())
P(f"   N = 12, Teiler q = {list(E)}:  ‖Σ E_q − I‖ = {np.linalg.norm(Ssum-np.eye(Nn)):.1e},  max ‖E_q² − E_q‖ = {max(np.linalg.norm(M@M-M) for M in E.values()):.1e},"
  f"  max ‖E_pE_q‖ (p≠q) = {max(np.linalg.norm(E[p]@E[r]) for p in E for r in E if p!=r):.1e}")
P(f"   Rang E_q = φ(q): {[int(round(np.trace(E[qq]))) for qq in E]}")
P("   E₁ = (1/N)ΣTⁿ ist der Schließungsprojektor 𝒫; die übrigen E_q trennen die geschlossenen Bahnen nach genauer Periode q (Resonanz p/q, Sternvieleck).")
sig=np.cos(2*np.pi*np.arange(60)*2/5)+0.5*np.cos(2*np.pi*np.arange(60)/3)+0.3*np.random.default_rng(1).normal(size=60)
T6=np.roll(np.eye(60),1,axis=0)
en={qq:np.linalg.norm(sum(cq(qq,n)*np.roll(sig,n) for n in range(60))/60)**2 for qq in (1,2,3,4,5,6,10,12,15)}
P("   Filterprobe (Signal mit Perioden 5 und 3 plus Rauschen, Länge 60), Energie je genaue Periode:")
P("   "+"  ".join(f"q={k}: {v:.1f}" for k,v in en.items()))
P()

# ---------------- E ----------------
P("E) Ramanujans Theorien zu den Basen r = 2, 3, 4, 6 und die Drehordnungen des Kaleidoskops")
P("   2cos(2π/r) ganzzahlig ⇔ r ∈ {1, 2, 3, 4, 6}:  " + ", ".join(f"r={r}: {2*np.cos(2*np.pi/r):+.3f}" for r in range(1,9)))
def a_c(qv,N=40):
    m=np.arange(-N,N+1); M,Nm=np.meshgrid(m,m); return np.sum(qv**(M*M+M*Nm+Nm*Nm))
def b_c(qv,N=40):
    w=np.exp(2j*np.pi/3); m=np.arange(-N,N+1); M,Nm=np.meshgrid(m,m); return np.sum(w**(M-Nm)*qv**(M*M+M*Nm+Nm*Nm)).real
def c_c(qv,N=40):
    m=np.arange(-N,N+1); M,Nm=np.meshgrid(m,m); return np.sum(qv**((M+1/3)**2+(M+1/3)*(Nm+1/3)+(Nm+1/3)**2))
for qv in (0.2,0.5):
    A,B,C=a_c(qv),b_c(qv),c_c(qv); P(f"   q = {qv}:  a³ = {A**3:.10f}   b³ + c³ = {B**3+C**3:.10f}   (kubische Identität)")
def count(form,Rmax):
    m=np.arange(-Rmax-2,Rmax+3); M,Nm=np.meshgrid(m,m); Q=form(M,Nm); return np.bincount(Q[Q<=Rmax*Rmax].ravel())
sq=count(lambda m,n: m*m+n*n,6); hx=count(lambda m,n: m*m+m*n+n*n,6)
P("   Bilder im Abstand² = k (Zellmaß 1): Quadratzelle r₂(k) = " + str([int(v) for v in sq[:14]]))
P("                                       Sechseckzelle     = " + str([int(v) for v in hx[:14]]))
P("   Die Thetareihen Σq^{m²+n²} und a(q) = Σq^{m²+mn+n²} zählen die Bilder je Abstand; 4 bzw. 6 nächste Bilder = Drehordnung.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"ramanujan_weitere.txt"),"w",encoding="utf-8").write(txt); print(txt)
