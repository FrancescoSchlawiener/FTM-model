# Netzwerk-Summe (03 §12) und Energiebilanz der Dunklen Energie (03 §13) des Autors, durchgerechnet
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from scipy.integrate import solve_ivp
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; Mpc=1e3*kpc
out=[]; P=lambda s="": out.append(s)
chi=lambda s: s*np.exp(-s)/(-np.expm1(-s))       # Schalter, a = a_N + χ(s)·√(a_N a0)
nu =lambda y: 1/(-np.expm1(-np.sqrt(y)))

# ---------- Teil A: Netzwerk-Summe ----------
P("A) Netzwerk-Summe: Galaxie (6e10 M☉, Hernquist a = 3 kpc) im Haufen")
GMg=6e10*Msun*G; ah=3*kpc; rsg=np.sqrt(GMg/a0)
S_eigen=GMg/(3*c**2); S_netz=G*1.7e14*Msun/(3*c**2)
P(f"   S_eigen = {S_eigen:.3e} m, S_Netz = {S_netz:.3e} m, Verhältnis {S_netz/S_eigen:.0f}")
r=15*kpc
aK=lambda S: (S/r)*(3*c**2/rsg)
P(f"   (1) linear eingesetzt in den Kern bei r = 15 kpc: a = (S/r)(3c²/r*) = {aK(S_eigen):.3e} (eigen) bzw. {aK(S_netz):.3e} m/s² (Netz)")
gobs_mw=(220e3)**2/(15*kpc)
P(f"       beobachtet (v ≈ 200–220 km/s) etwa {gobs_mw:.2e} m/s² → Netzsumme linear: Faktor {aK(S_netz)/gobs_mw:.0f} zu groß")
# (2) nichtlinear: äußeres Feld des Haufens (Baryonen 1e14 M☉ in 1 Mpc) am Ort der Galaxie
for Rcl in [0.5,1.0,2.0]:
    Mcl=1e14*Msun*min(1,(Rcl/1.0))    # grob: M_b(<R) ∝ R bis 1 Mpc
    geN=G*Mcl/(Rcl*Mpc)**2
    P(f"   (2) Galaxie bei R = {Rcl} Mpc im Haufen: g_e,N = {geN:.2e} m/s² = {geN/a0:.3f} a0")
    P("        r [kpc]   v isoliert [km/s]   v im Haufen [km/s]   Verhältnis")
    for rk in [5,10,20,40,80]:
        rr=rk*kpc; giN=GMg/(rr+ah)**2
        g_iso=giN*nu(giN/a0)
        g_efe=giN*nu(np.hypot(giN,geN)/a0)       # gemittelte Näherung für die innere Dynamik
        P(f"        {rk:5d}       {np.sqrt(g_iso*rr)/1e3:8.1f}            {np.sqrt(g_efe*rr)/1e3:8.1f}           {np.sqrt(g_efe/g_iso):.3f}")

# ---------- Teil B: Energiebilanz σ = ΣΔs_innen − ΔL_außen = 0 ----------
P(); P("B) Bilanz des Autors: ΔL_außen/R = ΣΔs/R über das beobachtbare Universum (Teilchenhorizont), H = d/dt(ΔL/R)")
Gly=9.4607e24; Gyr=3.15576e16
Ob=0.0493; H0=67.4e3/Mpc; rhoc=3*H0**2/(8*np.pi*G); rhob0=Ob*rhoc
K=4*np.pi*G*rhob0/(9*c**2)*Gly**2              # in Gly^-2
P(f"   ε = ΣΔs/r_p = (4πGρ_b/9c²)·r_p² = K·χ_p²/a,  K = {K:.3e} Gly⁻²;  heute mit χ_p = 46,5 Mrd Lj: ε = {K*46.5**2:.4f}")
# (i) im gemessenen Hintergrund (ΛCDM): dε/dt gegen H0
Om,Or,OL=0.315,9.1e-5,0.685
E=lambda z: np.sqrt(Om*(1+z)**3+Or*(1+z)**4+OL)
from scipy.integrate import quad
chip=quad(lambda z: 1/E(z),0,np.inf,limit=500)[0]*c/H0/Gly
H0g=H0*Gyr
deps=K*(2*chip*1.0-chip**2*H0g)                 # c = 1 Gly/Gyr, a = 1
P(f"   (i) im gemessenen Hintergrund: χ_p = {chip:.2f} Mrd Lj, ε = {K*chip**2:.4f}, dε/dt = {deps:.4f} /Gyr gegen H0 = {H0g:.4f} /Gyr → Verhältnis {deps/H0g:.3f}")
# (ii) selbstkonsistent: H = dε/dt ⇒ H = 2Kχ/(a²(1+ε))
def f(t,y):
    a,chi_=y; eps=K*chi_**2/a; H=2*K*chi_/(a**2*(1+eps))
    return [a*H, 1/a]
t0=1e-6; a_ini=1e-4
# Startwerte auf dem Potenzgesetz a ∝ t^(2/3): χ = 3 t/a, H = 2/(3t) ⇒ a aus H = 2Kχ/a² ⇒ a³ = 3K t² ... 
a_ini=(3*K)**(1/3)*t0**(2/3); chi_ini=3*t0/a_ini
sol=solve_ivp(f,[t0,1e4],[a_ini,chi_ini],rtol=1e-10,atol=1e-14,dense_output=True,
              events=lambda t,y: y[0]-1.0)
te=sol.t_events[0][0]; ae,che=sol.y_events[0][0]
epse=K*che**2/ae; He=2*K*che/(ae**2*(1+epse))
# Verzögerungsparameter numerisch
h=1e-3*te
aa=lambda t: sol.sol(t)[0]
add=(aa(te+h)-2*aa(te)+aa(te-h))/h**2; ad=(aa(te+h)-aa(te-h))/(2*h)
qd=-add*aa(te)/ad**2
P(f"   (ii) selbstkonsistent (H = dε/dt): Beginn wie a ∝ t^(2/3), solange ε klein ist.")
P(f"        Bei heutiger Baryonendichte (a = 1): t = {te:.1f} Gyr, H = {He/Gyr*Mpc/1e3:.2f} km/s/Mpc, ε = {epse:.4f}, q = {qd:.3f}")
P(f"        gemessen: t0 = 13,8 Gyr, H0 = 67–73 km/s/Mpc, q0 ≈ −0,55")
open("netz_und_bilanz.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))
