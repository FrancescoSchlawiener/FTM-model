# Formeln Ramanujans, deren Zweck unbekannt war oder ist: nachgerechnet und am Modell gemessen
# A Mock-Theta f(q) und Ramanujans Radialgrenze   B σ(q) aus dem verlorenen Notizbuch   C partielle Theta-Identität
# D Beinahe-Lösungen a³ + b³ = c³ ± 1 (Seite mit 1729)   E Q̃ als Summe unvollständiger Gammafunktionen
import sys, io, os
from mpmath import mp, mpf, qp, nsum, inf, quad, gammainc, polylog, exp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
mp.dps = 120
out=[]; P=lambda s="": out.append(s)

# ---------------- A ----------------
P("A) Mock-Theta f(q) = Σ q^{n²}/(−q;q)_n²  (letzter Brief 1920). Ramanujans Behauptung: an einer Einheitswurzel")
P("   gerader Ordnung 2k bleibt f(q) − (−1)^k·b(q) beschränkt, b(q) = (q;q²)_∞·(1 − 2q + 2q⁴ − 2q⁹ + …).")
def f_mock(q,N=4000):
    S=mpf(0); den=mpf(1)
    for n in range(N):
        if n>0: den*=(1+q**n)
        t=q**(n*n)/den**2; S+=t
        if n>20 and abs(t)<mpf(10)**(-30): break
    return S
def b_fun(q):
    th=1+2*sum((-1)**n*q**(n*n) for n in range(1,1200))
    pr=mpf(1); j=0
    while True:
        t=q**(2*j+1); pr*=(1-t); j+=1
        if abs(t)<mpf(10)**(-100): break
    return pr*th
P("   q = −r (radial gegen −1, Ordnung 2, k = 1):")
P("       r        f(q)              b(q)              f(q) + b(q)")
for r in (mpf('0.9'),mpf('0.97'),mpf('0.99'),mpf('0.995')):
    q=-r; fv=f_mock(q); bv=b_fun(q)
    P(f"   {float(r):6.3f}  {float(fv):16.6e}  {float(bv):16.6e}  {float(fv+bv):14.8f}")
P("   f und b wachsen ohne Grenze, ihre Summe bleibt beschränkt: die Behauptung hält.")
P()

# ---------------- B ----------------
P("B) σ(q) = Σ_{n≥0} q^{n(n+1)/2}/(−q;q)_n  (verlorenes Notizbuch): Koeffizienten")
N=1500
coef=[0]*(N+1)
# Potenzreihe: q^{n(n+1)/2} · Π_{j=1}^{n} 1/(1+q^j)
def mul_inv_1plus(a,j):   # a · 1/(1+q^j)
    b=a[:]
    for i in range(j,len(b)): b[i]-=b[i-j]
    return b
n=0
while n*(n+1)//2<=N:
    ser=[0]*(N+1); ser[n*(n+1)//2]=1
    for j in range(1,n+1): ser=mul_inv_1plus(ser,j)
    coef=[x+y for x,y in zip(coef,ser)]; n+=1
nz=sum(1 for c in coef[1:] if c==0)
P(f"   erste Koeffizienten: {coef[:16]}")
P(f"   bis q^{N}: {nz} von {N} Koeffizienten sind 0, größter Betrag {max(abs(c) for c in coef)}")
P()

# ---------------- C ----------------
P("C) Partielle Theta-Identität (verlorenes Notizbuch):")
P("   Σ qⁿ/((aq;q)_n(q/a;q)_n) = (1−a)Σ(−1)ⁿa^{3n}q^{n(3n+1)/2}(1 − a²q^{2n+1}) + Σ(−1)ⁿa^{2n+1}q^{n(n+1)/2} / ((aq;q)_∞(q/a;q)_∞)")
for a,q in ((mpf('0.4'),mpf('0.3')),(mpf('0.7'),mpf('0.55'))):
    L=nsum(lambda n: q**n/(qp(a*q,q,int(n))*qp(q/a,q,int(n))),[0,inf])
    R=(1-a)*nsum(lambda n: (-1)**n*a**(3*n)*q**(n*(3*n+1)/2)*(1-a**2*q**(2*n+1)),[0,inf]) \
      + nsum(lambda n: (-1)**n*a**(2*n+1)*q**(n*(n+1)/2),[0,inf])/(qp(a*q,q)*qp(q/a,q))
    P(f"   a = {float(a)}, q = {float(q)}:  links {float(L):.15f}   rechts {float(R):.15f}")
P("   Die rechte Seite enthält eine einseitige Thetareihe Σ_{n≥0}(−1)ⁿa^{2n+1}q^{n(n+1)/2}: nur n ≥ 0, nicht modular.")
P()

# ---------------- D ----------------
P("D) Seite mit 1729: Beinahe-Lösungen von a³ + b³ = c³")
for a,b,c in ((9,10,12),(64,94,103),(73,144,150),(135,138,172),(11161,11468,14258),(65601,67402,83802)):
    P(f"   {a}³ + {b}³ − {c}³ = {a**3+b**3-c**3:+d}")
P()

# ---------------- E ----------------
P("E) Die Strukturwirkung als Summe unvollständiger Gammafunktionen:")
P("   Q̃(s) = −4[s³Li₁ + 3s²Li₂ + 6sLi₃ + 6Li₄](e^{−s}) = −4∫_s^∞ x³/(eˣ−1)dx = −4·Σ_{n≥1} Γ(4, ns)/n⁴")
for s in (mpf('0.3'),mpf(1),mpf(3)):
    q=exp(-s)
    A=-4*(s**3*polylog(1,q)+3*s**2*polylog(2,q)+6*s*polylog(3,q)+6*polylog(4,q))
    B=-4*quad(lambda x: x**3/(exp(x)-1),[s,inf])
    C=-4*nsum(lambda n: gammainc(4,n*s)/n**4,[1,inf])
    P(f"   s = {float(s)}:  {float(A):.15f}  {float(B):.15f}  {float(C):.15f}")
P("   Γ(1−k, 4π|n|y) ist die Bauart des nichtholomorphen Teils harmonischer Maaß-Formen, deren holomorpher Teil die Mock-Theta sind.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"ramanujan_unbekannt.txt"),"w",encoding="utf-8").write(txt); print(txt)
