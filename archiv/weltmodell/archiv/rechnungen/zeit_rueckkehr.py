# Zeit der Rückkehr τ_D = L²/D. (1) Gleichwertigkeit L = r* ⇔ τ_D = c/a0; (2) Q̃ als thermische Größe:
# Zustandsdichte, Gewicht, Faktor 4, ob τ_D = β (KMS) aus der Wirkung folgt; (3) Kandidaten für τ_D;
# (4) a0(z) bei τ_D = β(z) = 2π/H(z).
import numpy as np, sys, io, os
import sympy as sp, mpmath as mp
from scipy import integrate
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; AU=1.496e11
Gyr=3.15576e16; Mpc=3.0857e22; lP=1.616255e-35
l0=c**2/a0; mp.mp.dps=30
out=[]; P=lambda s="": out.append(s)
nB=lambda s: 1/np.expm1(s)

# ---------------- 1) Gleichwertigkeit ----------------
P("1) Gleichwertigkeit L = r* ⇔ τ_D = c/a0")
Gs,M,r,cs,A,tau,Ds,g=sp.symbols('G M r c a_0 tau_D D g',positive=True)
aN=Gs*M/r**2; Gam=aN/cs
Dsol=sp.solve(sp.Eq(sp.sqrt(Gam/Ds),1/r),Ds)[0]
L=sp.sqrt(Dsol*tau); rst=sp.sqrt(Gs*M/A)
tsol=sp.solve(sp.Eq(L,rst),tau)
P(f"   Γ = a_N/c, √(Γ/D) = 1/r  ⇒  D = {Dsol}   (für alle r dieselbe)")
P(f"   L² = D·τ_D = {sp.expand(L**2)};   L = r* = √(GM/a0)  ⇒  τ_D = {tsol}")
P(f"   Probe L/r − √(a_N/a0) bei τ_D = c/a0: {sp.simplify((L/r).subs(tau,cs/A)-sp.sqrt(aN/A))}")
Dg=g*r**2/cs; Lg=sp.sqrt(Dg*tau)
P(f"   beliebige Verteilung, D = a_N r²/c:  L/r = {sp.simplify(Lg/r)};  = √(a_N/a0) genau für τ_D = c/a0")
for name,Mx in [("Sonne",Msun),("Milchstraße 6·10¹⁰ M☉",6e10*Msun)]:
    D=G*Mx/c; Lr=np.sqrt(G*Mx/a0)
    P(f"   {name:22s} D = {D:.3e} m²/s   L = r* = {Lr:.4e} m   τ_D = L²/D = {Lr**2/D/Gyr:.2f} Mrd. J.")
P()

# ---------------- 2) Q̃ als thermische Größe ----------------
P("2) Q̃ = −4[s³Li₁ + 3s²Li₂ + 6s·Li₃ + 6Li₄](e^{−s})")
x,n,s=sp.symbols('x n s',positive=True)
I=sp.integrate(x**3*sp.exp(-n*x),(x,s,sp.oo))
P(f"   ∫_s^∞ x³e^{{−nx}}dx = {sp.simplify(I*sp.exp(n*s))}·e^{{−ns}}  →  Σ_n: s³Li₁ + 3s²Li₂ + 6sLi₃ + 6Li₄, da 1/(eˣ−1) = Σe^{{−nx}}")
Li=lambda k,s: mp.polylog(k,mp.e**(-s))
Qp=lambda s: -4*(s**3*Li(1,s)+3*s**2*Li(2,s)+6*s*Li(3,s)+6*Li(4,s))
Qi=lambda s: -4*mp.quad(lambda t: t**3/mp.expm1(t),[s,s+20,mp.inf])
dev=max(abs(Qp(sv)/Qi(sv)-1) for sv in [1e-3,0.01,0.1,0.5,1,2,5,10,30])
P(f"   numerisch Q̃(Polylog) / (−4∫_s^∞ x³/(eˣ−1)dx) − 1, s = 10⁻³…30: max {float(dev):.1e};  Q̃(0) = −4π⁴/15 = {float(-4*mp.pi**4/15):.4f}")
# Feldgleichung
Y,l=sp.symbols('Y ell_0',positive=True)
sY=(l**2*Y)**sp.Rational(1,4)
dQds=4*s**3/(sp.exp(s)-1)
dLdY=sp.simplify((2/l**2)*dQds.subs(s,sY)*sp.diff(sY,Y))
P(f"   s⁴ = ℓ0²Y, dQ̃/ds = 4s³/(eˢ−1)  ⇒  ∂ℒ_𝒮/∂Y = {sp.simplify(dLdY.subs(Y,s**4/l**2))} = 2·n_B(s)")
P("   ∂ℒ/∂(∂σ) = −4B + 2·(∂ℒ/∂Y)·∇σ = 4[(1 + n_B)∇σ − a]  ⇒  ν = 1 + n_B(s),  q = n_B/(1 + n_B) = e^{−s}")
P("   Also Q̃ = −∫_s^∞ n_B(x) d(x⁴): Besetzung n_B, Modenmaß d(x⁴) = 4x³dx. Die 4 ist der Exponent von s⁴ = ℓ0²Y.")
P()
# Vergleich mit Einschleifen-Größen (Einheiten (kT)⁴/(2π²(ħc)³) je Freiheitsgrad)
f3=lambda s: mp.quad(lambda t: t**2*mp.log(1-mp.e**(-t)),[s,s+20,mp.inf])          # freie Energie 3D, scharf geschnitten
f3c=lambda s: -(s**2*Li(2,s)+2*s*Li(3,s)+2*Li(4,s))
f4=lambda s: mp.quad(lambda t: t**3*mp.log(1-mp.e**(-t)),[s,s+20,mp.inf])          # freie Energie 4D, geschnitten
f4c=lambda s: -(s**3*Li(2,s)+3*s**2*Li(3,s)+6*s*Li(4,s)+6*Li(5,s))
U3=lambda s: mp.quad(lambda t: t**3/mp.expm1(t),[s,s+20,mp.inf])                  # Energie 3D, geschnitten
pin=lambda s: s**3/3*mp.log(1-mp.e**(-s))+f3(s)                                    # 3D, Moden unter s auf s gehalten
fm=lambda s: mp.quad(lambda t: t**2*mp.log(1-mp.e**(-mp.sqrt(t*t+s*s))),[0,s+1,s+20,mp.inf])  # 3D mit Masse s/β
P(f"   Geschlossene Formen: f₃(s) = ∫_s^∞x²ln(1−e^{{−x}}) = −[s²Li₂ + 2sLi₃ + 2Li₄]  (Probe s = 1: {float(f3(1)):.6f} / {float(f3c(1)):.6f})")
P(f"                        f₄(s) = ∫_s^∞x³ln(1−e^{{−x}}) = −[s³Li₂ + 3s²Li₃ + 6sLi₄ + 6Li₅]  (Probe s = 1: {float(f4(1)):.6f} / {float(f4c(1)):.6f})")
P("      s      Q̃        −4·U₃     12·f₃(schnitt)  12·f₃(gehalten)  12·f₃(Masse s/β)")
for sv in [0.3,1,3]:
    P(f"   {sv:5.1f}  {float(Qp(sv)):8.4f}  {float(-4*U3(sv)):8.4f}   {float(12*f3(sv)):9.4f}      {float(12*pin(sv)):9.4f}        {float(12*fm(sv)):9.4f}")
P(f"   s = 0: alle −4π⁴/15 = {float(-4*mp.pi**4/15):.4f} (12·f₃(0) = −12π⁴/45).")
P("   Q̃ = −4·U₃ = 12·∫_0^∞ x² ln(1 − e^{−max(x,s)}) dx exakt; der scharfe Schnitt fehlt um den Randterm −4s³Li₁ = 4s³ln(1−e^{−s}).")
# Test: ν − 1 = dQ/d(s⁴) für jede Kandidatenform; Tully-Fisher verlangt s·ν → 1
def nu_from(Q,sv,h=1e-6):
    d=(Q(sv*(1+h))-Q(sv*(1-h)))/(2*sv*h); return 1+d/(4*sv**3)
P("   Test an der Feldgleichung: ν − 1 = dQ/d(s⁴), Tully-Fisher verlangt s·ν → 1 für s → 0:")
P("      s        Q̃ (−4U₃)   12·f₃ schnitt   12·f₃ Masse   f₄ schnitt (jeder Faktor)")
for sv in [0.1,0.01]:
    P(f"   {sv:6.2f}   {float(sv*nu_from(Qp,mp.mpf(sv))):9.4f}   {float(sv*nu_from(lambda t:12*f3c(t),mp.mpf(sv))):10.4f}   {float(sv*nu_from(lambda t:12*fm(t),mp.mpf(sv))):10.4f}   {float(sv*nu_from(lambda t:f4c(t),mp.mpf(sv))):10.4f}")
P("   geschnitten: ν − 1 = 3Li₁(e^{−s})/s (Logarithmus); Masse: ν − 1 → π²/(2s²), a → (π²/2)·a0; 4D: ν − 1 ∝ Li₁, a → a_N.")
P()
# KMS-Bildsumme in der Zeit gegen Windungen in w
P("   Bildsumme: KMS in imaginärer Zeit mit Periode β und Zerfallsrate ω: Σ_n e^{−nβω} = n_B(βω).")
P("   Mit der Rate der Gleichung (ω = Γ) und β = τ_D: βω = Γτ_D = s² → n_B(s²). Windungen in w: Σ_n e^{−n·L√(Γ/D)} = n_B(s).")
P("      s      ν = 1+n_B(s)   ν = 1+n_B(s²)")
for sv in [0.3,1,3]:
    P(f"   {sv:5.1f}   {1+nB(sv):9.4f}      {1+nB(sv*sv):9.4f}")
for sv in [1e-2,1e-3]:
    P(f"   s = {sv:g}: n_B(s): s·ν = {sv*(1+nB(sv)):.4f};   n_B(s²): s·ν = {sv*(1+nB(sv*sv)):.1f}  (a → a_N + a0, v² ∝ r)")
P("   β = τ_D gibt n_B(s) nur mit ω = s/τ_D = √(Γ/τ_D) = √(a_N·a0)/c; diese Rate kommt in D·G'' = Γ·G nicht vor.")
P()
# Lesarten von x
Mmw=6e10*Msun; rs_mw=np.sqrt(G*Mmw/a0); rs_sun=np.sqrt(G*Msun/a0)
vinf=lambda Mx:(G*Mx*a0)**0.25
P("   Lesarten von x (Milchstraße r* = %.2f kpc, v∞ = %.1f km/s; Sonne r* = %.0f AE, v∞ = %.0f m/s):"%(rs_mw/kpc,vinf(Mmw)/1e3,rs_sun/AU,vinf(Msun)))
P(f"   x = κ·L, κ = 1/r (Kreisbahnmode), L = r*: s = r*/r; Periode der Windungen L, Ausbreitungszeit L²/D = c/a0 = {c/a0/Gyr:.1f} Mrd. J. für jedes M")
P(f"   Welle in w mit c: Periode β = L/c: Sonne {rs_sun/c/86400:.1f} d, Milchstraße {rs_mw/c/3.156e7:.0f} a (hängt von M ab)")
P(f"   β = τ_D in der Wellen-Lesart verlangt die Geschwindigkeit L/τ_D = √(GMa0)/c = v∞²/c: Sonne {vinf(Msun)**2/c:.2e} m/s, Milchstraße {vinf(Mmw)**2/c:.3f} m/s")
P(f"   Horizontbad (Lichtmoden, kT = ħa0/c): x = k·ℓ0; die Kreisbahnmode bei r* hat x = ℓ0/r* = c²/v∞² = {l0/rs_mw:.2e} (Milchstraße), n_B = e^(−{l0/rs_mw:.1e})")
P()
# Normierung
P("   Normierung: ℒ_𝒮 trägt c⁴/(16πG)·(2/ℓ0²)·Q̃ (so gibt Q̃(0) genau Λ_s = 4π⁴/(15ℓ0²)); ein Freiheitsgrad bei kT = ħa0/c trägt")
P("   (kT)⁴/(2π²(ħc)³)·Q̃/12 = ħc·Q̃/(24π²ℓ0⁴). Verhältnis 3π·ℓ0²/ℓ_P² = %.2e. Die Wirkung legt s = βω und N·T⁴ fest, β allein nicht."%(3*np.pi*l0**2/lP**2))
P()

# ---------------- 3) Kandidaten für τ_D ----------------
P("3) Kandidaten für τ_D, a0 = c/τ_D   (flach, Ω_m = 0,315, Ω_Λ = 0,685 angenommen; Strahlung vernachlässigt)")
Om=0.315; OL=1-Om
E=lambda a: np.sqrt(Om/a**3+OL)
t0u=integrate.quad(lambda a: 1/(a*E(a)),0,1)[0]           # t0·H0
REu=integrate.quad(lambda a: 1/(a*a*E(a)),1,np.inf)[0]     # R_E·H0/c
P(f"   t0·H0 = {t0u:.4f},  R_E·H0/c = {REu:.4f},  H_Λ/H0 = √Ω_Λ = {np.sqrt(OL):.4f}")
P("   τ_D                              H0 = 67,4: τ_D [Mrd. J.]  a0 [10⁻¹⁰]  Abstand   |  H0 = 73: τ_D  a0  Abstand")
rows=[]
for lab,f in [("ℓ0/c (Wirkung, gemessenes a0)",None),("2π/H0 (KMS am Hubble-Radius)",lambda H:2*np.pi/H),
              ("2π/H_Λ (KMS de Sitter, exakt)",lambda H:2*np.pi/(H*np.sqrt(OL))),("2π·R_E/c (Ereignishorizont)",lambda H:2*np.pi*REu/H),
              ("1/H0",lambda H:1/H),("Weltalter t0",lambda H:t0u/H)]:
    line=f"   {lab:32s}"
    for h in [67.4,73]:
        H=h*1e3/Mpc
        td=c/a0 if f is None else f(H); aa=c/td
        line+=f"   {td/Gyr:8.1f}     {aa/1e-10:6.3f}   {aa/a0-1:+6.1%}  "
    P(line)
P(f"   Für a0 = 1,2·10⁻¹⁰ nötig: τ_D = {c/a0/Gyr:.1f} Mrd. J. = {c/a0*67.4e3/Mpc/(2*np.pi):.3f}·2π/H0 (67,4) bzw. {c/a0*73e3/Mpc/(2*np.pi):.3f}·2π/H0 (73);"
  f" H0 = {2*np.pi*a0/c*Mpc/1e3:.1f} km/s/Mpc")
P()

# ---------------- 4) τ_D = β(z) ----------------
P("4) τ_D = 2π/H(z)  ⇒  a0(z) = c·H(z)/2π,  a0(z)/a0(0) = H(z)/H0 (H(z) gemessen einsetzbar; hier ΛCDM wie oben)")
P("      z     H/H0    r*(z)/r*(0)   v∞ bei festem M_b   Δlog M_b bei festem v∞ [dex]")
for z in [0.5,1,2]:
    h=E(1/(1+z))
    P(f"   {z:5.1f}   {h:6.3f}     {h**-0.5:6.3f}        {h**0.25-1:+6.1%}             {-np.log10(h):+.3f}")
P("   Mit τ_D = ℓ0/c (Konstante) oder 2π/H_Λ: a0(z) = a0(0), keine Verschiebung.")
P("   Nullpunkt der baryonischen Tully-Fisher-Relation, v∞⁴ = G·M_b·a0(z), trennt beide.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"zeit_rueckkehr.txt"),"w",encoding="utf-8").write(txt); print(txt)
