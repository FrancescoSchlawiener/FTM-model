# Welcher Nome trägt den Kern? Kern (Geschwindigkeitsform) gegen Rückkopplungssumme.
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30
out=[]; P=lambda s="": out.append(s)

def delta(q, N=400):
    n=np.arange(1,N+1); return q*np.prod((1-q**n)**24)

# 1) Kern-Nome q_Δs = Δs/r in realen Systemen
P("1) Kern-Nome q_Δs = Δs/r und Δ(τ)/q_Δs − 1")
for name,M,r in [("Sonne, 1 AE",Msun,1.496e11),("Sonne, r* (7030 AE)",Msun,np.sqrt(G*Msun/a0)),
                 ("Milchstraße 6e10, 8 kpc",6e10*Msun,8*3.086e19),("Milchstraße, 50 kpc",6e10*Msun,50*3.086e19)]:
    ds=G*M/(3*c**2); q=ds/r
    P(f"   {name:26s} q_Δs = {q:.3e}   Δ/q−1 = {delta(q)/q-1:+.3e}   r*/r = {np.sqrt(G*M/a0)/r:.3e}")

# 2) Geforderter Kern für die volle RAR, Geschwindigkeitsform a_K = (κ0/λ_eff)·M·v, κ0/λ_eff = 3c²/(r* v_rot)
M=6e10*Msun; GM=G*M; rs=np.sqrt(GM/a0); vrot=(GM*a0)**0.25; ds=GM/(3*c**2)
x=np.logspace(-2,2,9)   # r/r*
P(); P("2) Geforderter Kern M_req, damit a_K = (ν−1)·g_N mit ν = 1/(1−q), q = e^{−r*/r}")
P("   Behauptung: M_req = (Δs/r)·√s·q/√(1−q),  s = r*/r = −ln q  (→ Δs/r für r ≫ r*)")
P("     r/r*      M_req/(Δs/r)   √s·q/√(1−q)   Δ(τ_fb)/M_req")
for xi in x:
    r=xi*rs; s=1/xi; q=np.exp(-s); gN=GM/r**2; nu=1/(1-q)
    v=np.sqrt(r*gN*nu); aK=(nu-1)*gN
    Mreq=aK/(3*c**2/(rs*vrot)*v)
    F=np.sqrt(s)*q/np.sqrt(1-q)
    P(f"   {xi:8.3f}   {Mreq/(ds/r):12.6f}   {F:12.6f}   {delta(q)/Mreq:12.3e}")

# 3) Einheitliche Form über den Rückkopplungs-Nome allein
P(); P("3) Δs/r = (v_rot²/3c²)·s, also M_req = (v_rot²/3c²)·s^{3/2}·q/√(1−q)")
for xi in [0.1,1,10]:
    r=xi*rs; s=1/xi; q=np.exp(-s)
    P(f"   r/r* = {xi:5.1f}:  Δs/r = {ds/r:.6e}   (v_rot²/3c²)·s = {vrot**2/(3*c**2)*s:.6e}")

# 4) Kreisbahn aus der Bewegungsgleichung mit diesem Kern (exakt, λ_eff konstant) gegen RAR
P(); P("4) Kreisbahn aus v²/r = g_N + (κ0/λ_eff)·M_req(r)·v, gelöst für v, gegen v_RAR")
k=3*c**2/(rs*vrot)
for xi in x:
    r=xi*rs; s=1/xi; q=np.exp(-s); gN=GM/r**2
    Mk=vrot**2/(3*c**2)*s**1.5*q/np.sqrt(1-q)
    b=k*Mk*r; v=(b+np.sqrt(b*b+4*gN*r))/2          # v² − b v − gN r = 0
    vrar=np.sqrt(r*gN/(1-q))
    P(f"   r/r* = {xi:8.3f}   v = {v/1e3:9.3f} km/s   v_RAR = {vrar/1e3:9.3f} km/s   Δ = {v/vrar-1:+.1e}")

# 5) Andere Gedächtnisbereiche
P(); P("5) Langes Gedächtnis (λ_m ≪ ω, a_K = κ0·M·r): M_req ∝ s³·q/(1−q) → ∝ 1/r² für r ≫ r*,")
P("   also nicht Δs/r. Die Form Δs/r des Kerns gehört zur Geschwindigkeitsform (λ_eff fest).")
open("kern_nome.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))

# 6) Kern nur mit q_Δs (Summenform) gegen einheitlichen Kern
out=[]
P("6) Kern mit q_Δs allein (a_K = √(GM a0)/r) gegen einheitlichen Kern (a_K = (ν−1)·g_N)")
for xi in [0.3,1,3,10]:
    s=1/xi; q=np.exp(-s)
    v2sum=xi**-1+1; v2rar=xi**-1/(1-q)
    P(f"   r/r* = {xi:5.1f}:  v_Summe/v_RAR = {np.sqrt(v2sum/v2rar):.4f}")
GMs=G*Msun; rS=np.sqrt(GMs/a0); AU=1.496e11
for name,r in [("Erde (1 AE)",AU),("Saturn (9,5 AE)",9.54*AU)]:
    s=rS/r; gN=GMs/r**2
    P(f"   {name:16s} Zusatz q_Δs-Kern: {np.sqrt(GMs*a0)/r:.2e} m/s²   einheitlicher Kern: e^(−{s:.0f})·g_N ≈ 10^{(np.log10(gN)-s/np.log(10)):.0f} m/s²")
open("kern_nome.txt","a",encoding="utf-8").write("\n\n"+"\n".join(out)); print("\n".join(out))
