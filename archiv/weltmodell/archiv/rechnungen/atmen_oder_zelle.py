# B9: Atmen a0(z) ∝ H(z) gegen feste Zelle a0(z) ∝ 1+z
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from scipy.optimize import brentq
Om,Or=0.315,9.1e-5; OL=1-Om-Or; H0=67.4e3/3.0857e22
c=2.99792458e8; G=6.674e-11; hbar=1.054571817e-34; mp=1.67262192e-27
E=lambda z: np.sqrt(Om*(1+z)**3+Or*(1+z)**4+OL)
out=[]; P=lambda s="": out.append(s)

P("1) Schwelle, tiefer Bereich g_obs/g_bar ∝ √a0, Tully-Fisher v⁴ = G M a0")
P("   z     a0 Atmen  a0 Zelle   √-Verh. A/Z   Δlog v (A−Z)   Δlog M_bar bei festem v (A−Z)")
for z in [0.5,1,2,3,4,6,8,10]:
    A=E(z); Z=1+z
    P(f"  {z:4.1f}   {A:8.3f}  {Z:8.3f}    {np.sqrt(A/Z):8.4f}     {0.25*np.log10(A/Z):+.4f} dex     {-np.log10(A/Z):+.4f} dex")
zx=brentq(lambda z: E(z)-(1+z),0.5,5); P(f"   Gleichstand E(z) = 1+z bei z = {zx:.3f}")

P(); P("2) Heutige Drift der Schwelle d ln a0/dt")
q0=Om/2+Or-OL; Gyr=3.15576e16
P(f"   Atmen: Ḣ/H = −(1+q0)H0 = {-(1+q0):.4f} H0 = {-(1+q0)*H0*Gyr*100:.2f} %/Gyr;  v_flach: {-(1+q0)*H0*Gyr*25:.2f} %/Gyr")
P(f"   Zelle: −H0 = {-H0*Gyr*100:.2f} %/Gyr;  v_flach: {-H0*Gyr*25:.2f} %/Gyr")

P(); P("3) Proton-Kreuzung a0(z) = G m_p³c²/ħ²  (r*_p = λ̄_p, λ̄_p = ħ/(m_p c))")
a0x=G*mp**3*c**2/hbar**2; a00=c*H0/(2*np.pi); R=a0x/a00
zA=brentq(lambda z: E(z)-R,10,1e6)
P(f"   a0 nötig = {a0x:.4e} m/s²,  a0 heute (cH0/2π, H0 = 67,4) = {a00:.4e},  Verhältnis {R:.4e}")
P(f"   Atmen: z = {zA:.1f}     feste Zelle: z = {R-1:.4e}")
kB=8.617333e-5; T0=2.7255
P(f"   Plasmatemperatur: Atmen kT = {kB*T0*(1+zA):.4f} eV,  Zelle kT = {kB*T0*R:.3f} eV")

P(); P("4) Zelle = Horizont nur zu einer Epoche: L_phys = c²/a0")
P("   mitbewegter Umfang 2πc/(aH) in Einheiten des heutigen:  z: Wert")
for z in [0,0.5,1,2,5,1100]:
    P(f"   {z:7.1f}: {(1+z)/E(z):.4f}")
open("atmen_oder_zelle.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))

# 5) Massereiche Scheibe früh: Anteil der Zusatzbeschleunigung f = 1 − 1/ν innerhalb R
out=[]
P("5) Scheibe mit M_b = 1e11 M☉ (Punktmassen-Näherung), f_zus = 1 − 1/ν, a0(0) = 1,2e-10")
Msun=1.989e30; kpc=3.0857e19; a00=1.2e-10
nuf=lambda y: 1/(-np.expm1(-np.sqrt(y)))
for R in [5,10]:
    gN=G*1e11*Msun/(R*kpc)**2
    for z in [0,1,2]:
        A=E(z); Z=1+z
        fa=1-1/nuf(gN/(a00*A)); fz=1-1/nuf(gN/(a00*Z))
        P(f"   R = {R:2d} kpc, z = {z}: f_zus Atmen = {fa:.3f}, Zelle = {fz:.3f}")
open("atmen_oder_zelle.txt","a",encoding="utf-8").write("\n\n"+"\n".join(out)); print("\n".join(out))
