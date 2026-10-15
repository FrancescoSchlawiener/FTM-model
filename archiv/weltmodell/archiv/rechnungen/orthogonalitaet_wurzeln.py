# Orthogonalität, Quadrate, Wurzeln, Vielecke im Modell.
# 1 Senkrechte Anteile: Quadratsumme, Minkowski-Form, komplexe Rate, Epizyklen-Summenregel, Außenfeld, Thaleskreis-Kriterium
# 2 Wurzeln als Höhen im Thaleskreis: Regel, Kosten-Minimum der Rückkehr, Kerr, Maxwell-Element, Irrfahrt, Leiter der Skalen
# 3 Kreis als Unendlich-Eck: Archimedes als Kathetensatz, Rückkehrsumme als Rechteckregel, Halbgerade als Stufenleiter
# 4 Umfang zu Fläche eine Stufe höher und tiefer: geschlossene Rückkehrrichtung (d = 1), Quellen (d = 0)
import numpy as np, sys, io, os, math
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
import warnings, mpmath as mpm
mpm.mp.dps=60
np.seterr(over='ignore'); warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; AU=1.495978707e11
hbar=1.054571817e-34; mp=1.67262192e-27; me=9.1093837e-31; ly=9.4607e15; yr=3.15576e7
lP=np.sqrt(hbar*G/c**3); l0=c**2/a0; pi=np.pi
out=[]; P=lambda s="": out.append(s)
chi=lambda s: s/np.expm1(s)
nu=lambda s: 1/(-np.expm1(-s))
ds3=lambda M: G*M/c**2          # 3Δs
rstar=lambda M: np.sqrt(G*M/a0)
lbar=lambda M: hbar/(M*c)

# ================= 1) Senkrechte Anteile =================
P("1a) Normalform θ̇ = ω, ρ̇ = −κρ in z = ρe^{iθ}: ż = (−κ + iω)z. Transversal = Realteil, tangential = Imaginärteil")
for kap,om in ((0.5,1.0),(1.0,0.2)):
    t=0.7; h=1e-6; z=lambda t: np.exp(-kap*t)*np.exp(1j*om*t)
    lam=(z(t+h)-z(t-h))/(2*h)/z(t)
    P(f"   κ = {kap}, ω = {om}: ż/z = {lam.real:+.6f}{lam.imag:+.6f}i,  |ż|/|z| = {abs(lam):.6f} = √(κ² + ω²) = {np.hypot(kap,om):.6f}")
P()
P("1b) Einrasten ψ̇ = Δω − 2K_d sin ψ: Schwebung = √(Δω² − (2K_d)²) (Kathete zur Hypotenuse Δω);")
P("    Amplitudenverhältnis r = (Δω − √(Δω² − 4K_d²))/2K_d = tan(θ/2), sin θ = 2K_d/Δω")
dw=1.0
for Kd in (0.1,0.3,0.45,0.499):
    T=quad(lambda p: 1/(dw-2*Kd*np.sin(p)),0,2*pi,limit=400)[0]
    Wb=np.sqrt(dw**2-4*Kd**2); th=np.arcsin(2*Kd/dw)
    P(f"   K_d = {Kd}: 2π/T = {2*pi/T:.6f}, √(Δω² − 4K_d²) = {Wb:.6f};  r = {(dw-Wb)/(2*Kd):.6f}, tan(θ/2) = {np.tan(th/2):.6f}, θ = {np.degrees(th):.2f}°")
P("   Einrasten bei θ = 90°: r = tan 45° = 1 (gleiche Amplituden). Kopplung K·e^{iφ}: K² = (K cos φ)² + (K sin φ)², nur K cos φ rastet ein.")
P()
P("1c) Pendelkanal ω² = ω_el² ± ω_g²; χ_g = ω_g²/(ω_el² + ω_g²) = sin²θ mit tan θ = ω_g/ω_el")
for deg in (10,30,45):
    th=np.radians(deg); wg=np.tan(th)
    P(f"   θ = {deg:2d}°: χ_g = {wg**2/(1+wg**2):.4f}, hängend ω²/ω_el² = {1+wg**2:.4f}, stehend {1-wg**2:+.4f}")
P("   Stehend knickt der Körper bei θ = 45° (χ_g = ½): Nullrichtung ω_g = ω_el, Diagonale des Quadrats.")
P()
P("1d) Nome als Umlauffaktor q = e^{2πiτ}: |q| = e^{−2π Im τ} (Zusammenziehen je Rücklauf), arg q = 2π Re τ (Drehung je Rücklauf)")
P("    Punkt 1 (Im τ = s/2π, |q| = e^{−s}), Re τ = √(1 + χ):")
P("   r/r*   |q|      Re τ     ν aus |q|   Re 1/(1−q)   |1/(1−q)|")
for x in (1,3,30):
    s=1/x; q=np.exp(-s)*np.exp(2j*pi*np.sqrt(1+chi(s)))
    P(f"   {x:4g}   {abs(q):.4f}   {np.sqrt(1+chi(s)):.4f}   {nu(s):8.4f}   {(1/(1-q)).real:8.4f}     {abs(1/(1-q)):.4f}")
s=2*pi; rt=np.sqrt(1+chi(s))
P(f"   Quadrat: Im τ = 1 bei s = 2π (a_N = 4π²a₀, r = r*/2π); dort Re τ = {rt:.5f} ≡ {rt-1:.5f} (mod 1): τ ≈ i bis auf Scherung {rt-1:.4f}.")
P(f"   Punkt 3: Im τ = ln(r/Δs)/2π ≥ 2,5 in Galaxien, kein Quadrat. Kepler: Re τ ≡ 0, Gitter rechtwinklig; Regel schert um √(1+χ) − 1 ≤ {np.sqrt(2)-1:.4f}.")
P()
P("1e) Epizyklen-Summenregel (Laplace): κ² + ν_z² − 2Ω² = 4πGρ. Leerer Raum: κ² + ν_z² = (√2·Ω)².")
P("    Regel, Punktmasse (kugelsymmetrisch, ν_z² = a/r), Einheiten GM = a₀ = 1:")
lna=lambda r: np.log(nu(1/r)/r**2)
h=1e-5
P("   r/r*   κ²/Ω²     ν_z²/Ω²   (κ²+ν_z²)/Ω²   2 + χ      4πGρ_X/Ω² (aus ∇·(a − a_N))   χ      2ζ/Ω = κ²/Ω²")
for x in (0.3,1,3,10,30):
    a=lambda r: np.exp(lna(r)); Om2=lambda r: a(r)/r
    k2=(Om2(x*(1+h))-Om2(x*(1-h)))/(2*h)+4*Om2(x)
    nz2=a(x)/x
    dX=lambda r: r**2*(a(r)-1/r**2)
    divX=(dX(x*(1+h))-dX(x*(1-h)))/(2*h*x)/x**2
    v=lambda r: np.sqrt(a(r)*r)
    zeta=((x*(1+h))*v(x*(1+h))-(x*(1-h))*v(x*(1-h)))/(2*h*x)/x
    Om=np.sqrt(Om2(x))
    P(f"   {x:4g}   {k2/Om2(x):.5f}   {nz2/Om2(x):.5f}   {(k2+nz2)/Om2(x):.5f}        {2+chi(1/x):.5f}    {divX/Om2(x):.5f}                       {chi(1/x):.5f}   {2*zeta/Om:.5f}")
P("   Der Fehlbetrag zur Diagonale √2·Ω ist χΩ²; genau diese Zahl liest der Beobachter als Dichte ρ_X = χΩ²/(4πG).")
P("   (Re τ)² = κ²/Ω² = 2ζ/Ω (ζ = Wirbelstärke): κ = √(2Ω·ζ) ist ein Produkt (Höhe), keine Quadratsumme.")
P()
P("1f) Y = h^{μν}∂σ∂σ = |∇σ|²: die Regel sieht den Betrag des Gesamtfelds. s² = |a_N,ges|/a₀.")
M=6e10*Msun; r=50*kpc; ai=G*M/r**2/a0; ae=0.03
si=np.sqrt(ai); sperp=(ai**2+ae**2)**0.25; spar=np.sqrt(ai+ae)
P(f"   Milchstraße 50 kpc: a_N/a₀ = {ai:.4f}, Außenfeld 0,03 a₀")
P(f"   frei: s = {si:.4f}, ν = {nu(si):.3f};  senkrecht (s⁴ = s_i⁴ + s_e⁴): s = {sperp:.4f}, ν = {nu(sperp):.3f};  parallel (s² = s_i² + s_e²): s = {spar:.4f}, ν = {nu(spar):.3f}")
P()
P("1g) Onsager-Casimir / GENERIC: L = L_s + L_a steht senkrecht in der Frobenius-Norm, xᵀL_a x = 0 (keine Arbeit)")
rng=np.random.default_rng(1); L=rng.normal(size=(4,4)); Ls=(L+L.T)/2; La=(L-L.T)/2; x=rng.normal(size=4)
P(f"   tr(L_sᵀL_a) = {np.trace(Ls.T@La):.1e};  ‖L‖² − ‖L_s‖² − ‖L_a‖² = {np.sum(L**2)-np.sum(Ls**2)-np.sum(La**2):.1e};  xᵀL_a x = {x@La@x:.1e}")
P()
P("1h) Thales-Kriterium: F = PF + (I−P)F ist Mode für Mode senkrecht ⇔ Re(p·(1−p)*) = 0 ⇔ |p − ½| = ½ (Eigenwert p sieht [0,1] unter 90°)")
th=np.linspace(0,2*pi,2001)
for q in (np.exp(-3),1/3,np.exp(-1),0.9,0.999):
    p=(1-q)/(1-q*np.exp(1j*th))
    rad=np.max(abs(abs(p-1/(1+q))-q/(1+q)))
    dfk=1-abs(p)**2-abs(1-p)**2
    P(f"   𝒫_q, q = {q:.4f}: Eigenwerte auf Kreis um 1/(1+q) mit Radius q/(1+q) (Abw. {rad:.0e}); größter Pythagoras-Fehlbetrag {dfk.max():.4f} = 4q(1−q)/(1+q)² = {4*q*(1-q)/(1+q)**2:.4f}")
P("   Größter Fehlbetrag ½ bei q = ⅓ (s = ln 3 = 1,099, r = 0,91 r*). Nur q → 0 (Newton) und q → 1 (geschlossen) sind senkrecht.")
xx=np.logspace(-2,2,9); w=1j*xx/(1+1j*xx)
P(f"   Maxwell-Element G*/γ = iωt/(1+iωt): max||w − ½| − ½| = {np.max(abs(abs(w-0.5)-0.5)):.1e}, max|Re(w(1−w)*)| = {np.max(abs((w*np.conj(1-w)).real)):.1e}: liegt ganz auf dem Thaleskreis.")
P("   Π⊥ (orthogonaler Projektor): Eigenwerte 0 und 1, die Endpunkte des Durchmessers.")
P()
P("1i) Minkowski-Form: Takt auf der Kreisbahn (Schwarzschild, Ω² = GM/r³ exakt)")
P("    1 = (dτ/dt)² + 6Δs/r + (v/c)²,  v/c = rΩ/c  ⇒  (dτ/dt)² = 1 − 9Δs/r")
for rr in (100,30,12,9):
    ut=1/np.sqrt(1-9/rr)   # u^t der Kreis-Geodäte, r in Δs
    P(f"   r = {rr:3d}Δs: 1/(u^t)² = {1/ut**2:.6f},  1 − 6Δs/r − (rΩ/c)² = {1-6/rr-3/rr:.6f}")
P("   Stillstand auf der Kreisbahn bei r = 9Δs (Photonensphäre), ruhend bei r₀ = 6Δs.")
P("   g_tt·g_rr = −1: √(−g_tt)·√(g_rr) = 1, Takt und Faltung sind Kehrwerte (geometrisches Mittel 1) → gleiche Hälften 6Δs/b der Ablenkung.")
P("   Schließungsregel (∂ξ)² = −1: (∂_tξ)² = 1 + |∇ξ|², γ² − (γv/c)² = 1 (Einheitshyperbel) → Staub.")
P()

# ================= 2) Wurzeln als Höhen =================
P("2a) Höhensatz h = √(p·q) im Thaleskreis über p + q; Kathetenverhältnis √(p/q)")
Mmw=6e10*Msun; rsmw=rstar(Mmw)
P("   Regel (Milchstraße): p = a_N, q = a₀, h = √(a_N a₀) (Tiefbereichsglied); Kathetenverhältnis = s, Winkel = arctan s")
for rk in (2,5,8.35,20,50):
    aN=G*Mmw/(rk*kpc)**2; s=np.sqrt(aN/a0)
    P(f"   r = {rk:5.2f} kpc: a_N/a₀ = {aN/a0:7.3f}, h/a₀ = {s:.3f}, h/a_N = {1/s:.3f}, Winkel {np.degrees(np.arctan(s)):5.1f}°")
P("   r*: p = q, Winkel 45°, h = a_N = a₀; innen dominiert p, außen h. Nur bei p = q ist AM = GM.")
P(f"   χ ist kein Mittel: a = a_N + χ·h, bei p = q χ(1) = {chi(1):.4f} (nicht ½), a = {1+chi(1):.4f} a₀")
P()
P("2b) Rückkehr als Kosten-Minimum: ψ(t) = L/√(4πDt³)·e^{−L²/4Dt} (erste Rückkehr), ∫ψ e^{−Γt}dt = e^{−s}")
P("    s = min_t [Γt + L²/(4Dt)] = min_t [t/τ_Γ + τ_D/(4t)], τ_Γ = c/a_N, τ_D = L²/D = c/a₀ (AM ≥ GM)")
for s in (0.3,1.0,3.0):
    L=1.0; D=1.0; Gam=s**2   # s² = Γ L²/D
    lap=quad(lambda t: L/np.sqrt(4*pi*D*t**3)*np.exp(-L**2/(4*D*t)-Gam*t),0,np.inf,limit=400)[0]
    f=lambda t: Gam*t+L**2/(4*D*t); m=minimize_scalar(f,bounds=(1e-6,1e3),method='bounded',options={'xatol':1e-12})
    tG=1/Gam; tD=L**2/D
    P(f"   s = {s}: ∫ψe^{{−Γt}} = {lap:.6f}, e^{{−s}} = {np.exp(-s):.6f};  Minimum {m.fun:.6f} bei t = {m.x:.5f} = ½√(τ_Γτ_D) = {0.5*np.sqrt(tG*tD):.5f}; Anteile {Gam*m.x:.4f} + {L**2/(4*D*m.x):.4f}")
P("   Beide Kosten sind am Minimum gleich (je s/2). r* ist der Ort τ_Γ = τ_D: das Dreieck der Zeiten ist gleichschenklig.")
P()
P("2c) Kerr: r_± = 3Δs ± √((3Δs)² − a_J²)  ⇒  r_+ + r_− = 6Δs = r₀ (Durchmesser), r_+·r_− = a_J² (Höhe)")
for j in (0,0.5,0.9,1.0):
    d3=1.0; aJ=j*d3; rp=d3+np.sqrt(d3**2-aJ**2); rm=d3-np.sqrt(d3**2-aJ**2)
    P(f"   j = {j}: r_+ = {rp:.4f}, r_− = {rm:.4f} (·3Δs), r_+ + r_− = {rp+rm:.4f}, r_+r_− − a_J² = {rp*rm-aJ**2:.1e}, Fläche 4π(r_+² + a_J²)/(4π r₀ r_+) = {(rp**2+aJ**2)/(2*d3*rp):.6f}")
P("   j = 1: p = q, Höhe = Radius (extrem). Ergosphäre r_E(θ): dieselbe Konstruktion mit Höhe a_J cos θ; Äquator Höhe 0 → r_E = r₀.")
P()
P("2d) Maxwell-Element: Speicher G' = γx²/(1+x²), Verlust G'' = γx/(1+x²), x = ωt_K:  G''² = G'·(γ − G')")
for x in (0.1,1,10):
    Gs=x**2/(1+x**2); Gl=x/(1+x**2)
    P(f"   x = {x:4}: G''² = {Gl**2:.6f}, G'(γ−G') = {Gs*(1-Gs):.6f}")
P("   Größter Verlust bei x = 1: G' = γ − G' (p = q), G'' = γ/2 (Radius).")
P()
P("2e) Wurzelwachstum = ungeordnete senkrechte Schritte: |Σv_i|² = Σ|v_i|² (Kreuzterme mitteln weg) → Diagonale des N-Würfels")
rng=np.random.default_rng(7)
for N,K in ((100,20000),(10000,1000)):
    R2=0.0
    for _ in range(K//1000 if N<=100 else K):
        if N<=100:
            u=rng.normal(size=(1000,N,3)); u/=np.linalg.norm(u,axis=2,keepdims=True); R2+=np.sum(np.sum(u,axis=1)**2)
        else:
            u=rng.normal(size=(N,3)); u/=np.linalg.norm(u,axis=1,keepdims=True); R2+=np.sum(u.sum(0)**2)
    R2/=K
    P(f"   N = {N:5d} Einheitsschritte zufälliger Richtung ({K} Läufe): ⟨R²⟩/N = {R2/N:.3f} ± {np.sqrt(2/3/K):.3f};  gleichgerichtet R = N, geschlossenes N-Eck R = 0")
P("   Modell: D = L²/τ_D = r*²a₀/c = GM/c = c·3Δs  ⇒  r*² = 3Δs·ℓ₀ = N·(3Δs)², N = ℓ₀/3Δs")
for nm,M in (("Sonne",Msun),("Milchstraße",Mmw)):
    d3=ds3(M); N=l0/d3; r0=2*d3
    P(f"   {nm:11s}: 3Δs = {d3:.4e} m, N = {N:.3e}, √N·3Δs = {np.sqrt(N)*d3:.4e} m = r* = {rstar(M):.4e} m;  Schritte der Länge r₀: ℓ₀/r₀ = {l0/r0:.3e} = 2(r*/r₀)² = {2*(rstar(M)/r0)**2:.3e}")
P()
P("2f) Zeit an der Schwelle: t* = √((r*/c)·τ_D), τ_D = c/a₀ (= 2π/H mit a₀ = cH/2π)")
for nm,M in (("Sonne",Msun),("Milchstraße",Mmw)):
    rs=rstar(M); ts=np.sqrt(rs/c*c/a0); Om=np.sqrt(G*M/rs**3)
    P(f"   {nm:11s}: r*/c = {rs/c/yr:.4g} a, τ_D = {c/a0/yr:.4g} a, t* = {ts/yr:.4g} a = 1/Ω_Kepler(r*) = {1/Om/yr:.4g} a;  f*_gw = 1/(πt*) = {1/(pi*ts):.3e} Hz")
P()
P("2g) Leiter der Mittel: 3Δs_p – ℓ_P – λ̄_p – 3Δs_* (M_* = m_Pl³/m_p²), jede Länge Höhe über ihren Nachbarn")
mPl=np.sqrt(hbar*c/G); Mst=mPl**3/mp**2
lad=[("3Δs_p",ds3(mp)),("ℓ_P",lP),("λ̄_p",lbar(mp)),("3Δs_*",ds3(Mst))]
aG=G*mp**2/(hbar*c); xs=-0.5*np.log(aG)
for i in range(3):
    P(f"   {lad[i+1][0]}/{lad[i][0]} = {lad[i+1][1]/lad[i][1]:.4e}")
P(f"   Schritt α_G^(−1/2) = {aG**-0.5:.4e}, x = ln = {xs:.3f};  M_* = {Mst/Msun:.3f} M☉")
P(f"   ℓ_P² = 3Δs·λ̄ (jede Masse): Proton {ds3(mp)*lbar(mp)/lP**2:.6f}, Elektron {ds3(me)*lbar(me)/lP**2:.6f};  λ̄_p² = ℓ_P·3Δs_*: {lbar(mp)**2/(lP*ds3(Mst)):.6f}")
P(f"   Spin ½: a_J = λ̄/2, a_J·3Δs = ℓ_P²/2 (Elektron {lbar(me)/2*ds3(me)/lP**2:.6f})")
P("   Lage in Einheiten x über ℓ_P:")
for nm,Lx in (("3Δs_p",ds3(mp)),("λ̄_p",lbar(mp)),("3Δs_*",ds3(Mst)),("r*_p",rstar(mp)),("r*☉",rstar(Msun)),("r*_MW",rsmw),("ℓ₀",l0)):
    P(f"   {nm:6s}: {np.log(Lx/lP)/xs:+.3f}")
P("   ℓ₀ liegt bei 3,23 x, nicht auf der Leiter: die Schwelle bringt eine eigene Länge.")
for nm,M in (("Proton",mp),("Sonne",Msun),("Milchstraße",Mmw)):
    P(f"   r*²·λ̄/(ℓ_P²·ℓ₀) {nm:11s} = {rstar(M)**2*lbar(M)/(lP**2*l0):.6f}  (massenunabhängig)")
P()
P("2h) Vakuumteil Λ_s = 4π⁴/(15ℓ₀²) gegen Tiefbereich: √(GMa₀)/r = Λ_s c² r/3  ⇒  r_Λ² = (45/4π⁴)·ℓ₀·r*")
Ls=4*pi**4/(15*l0**2)
for nm,M in (("Sonne",Msun),("Milchstraße",Mmw)):
    rs=rstar(M); rL=brentq(lambda r: np.sqrt(G*M*a0)/r-Ls*c**2*r/3,1e10,1e27)
    d3=ds3(M); fr=lambda X: np.log(X/d3)/np.log(l0/d3)
    P(f"   {nm:11s}: r_Λ = {rL:.4e} m = {rL/kpc/1e3:.3f} Mpc, √(45/4π⁴)·√(ℓ₀r*) = {np.sqrt(45/(4*pi**4))*np.sqrt(l0*rs):.4e} m;  Lage zwischen 3Δs (0) und ℓ₀ (1): r* {fr(rs):.3f}, r_Λ {fr(rL):.3f}")
P(f"   Faktor √(45/4π⁴) = {np.sqrt(45/(4*pi**4)):.4f}; mit voller Dunkler Energie (3·Ω_s) wäre r_Λ um √3 kleiner.")
P()
P("2i) p = q hängt bei Längen an der Faktorwahl: 3Δs·ℓ₀ = r₀·(ℓ₀/2) = r₀·πR_H (a₀ = cH/2π)")
M1=c**2*l0/G; M2=c**2*(l0/2)/(2*G)   # 3Δs = ℓ₀ ; 6Δs = ℓ₀/2
P(f"   3Δs = ℓ₀ bei M = c⁴/(Ga₀) = {M1/Msun:.3e} M☉;  r₀ = πR_H bei M = c⁴/(4Ga₀) = {M2/Msun:.3e} M☉ (Faktor 4). Bei Beschleunigungen (a_N gegen a₀) gibt es diese Freiheit nicht.")
P()

# ================= 3) Kreis als Unendlich-Eck =================
P("3a) Archimedes: s_{2n}² = 2·(1 − √(1 − s_n²/4)) = Durchmesser · Pfeilhöhe (Kathetensatz), Einheitskreis, s_6 = 1")
sn=mpm.mpf(1); sn_naiv=1.0; n=6
P("   n            n·s_n/2 (60 Stellen)   Fehler·n²/π³   naiv in Gleitkomma √(2 − √(4 − s_n²))")
for k in range(0,31):
    if k in (0,4,10,17,20,26,30):
        P(f"   {n:<12d} {mpm.nstr(n*sn/2,16):18s}  {mpm.nstr((mpm.pi-n*sn/2)*n**2/mpm.pi**3,7):10s}     {n*sn_naiv/2:.15f}")
    sn=sn/mpm.sqrt(2+mpm.sqrt(4-sn**2)); sn_naiv=np.sqrt(max(2-np.sqrt(4-sn_naiv**2),0)); n*=2
P("   Fehler → π³/(6n²) (Faktor 1/6); die Form mit Differenz fast gleicher Zahlen verliert in Gleitkomma ab n ≈ 10⁶ die Stellen.")
P()
P("3b) Rückkehrsumme als Rechteck-Vieleck: s·Σ_{n≥0}e^{−ns} = s·ν nähert ∫₀^∞e^{−x}dx = 1 mit Schrittweite s = r*/r")
P("   s       r/r*    sν − 1      (sν−1)/s    beidseitig (s/2)coth(s/2) − 1   /s²")
for s in (1,0.3,0.1,0.03,0.01):
    e1=s*nu(s)-1; e2=(s/2)/np.tanh(s/2)-1
    P(f"   {s:<6}  {1/s:6.1f}   {e1:.6f}    {e1/s:.5f}     {e2:.3e}                      {e2/s**2:.5f}")
P("   Einseitig (ein Umlaufsinn) Fehler s/2 + s²/12: daher M_schein = r/r* + ½ + r*/(12r). Beidseitig nur s²/12, wie das Vieleck ∝ 1/n².")
P()
P("3c) Halbgerade als Stufenleiter (Schritt h = r₀, N = r*/r₀ Stufen): D(G_{n+1} − 2G_n + G_{n−1})/h² = ΓG_n")
P("    exakt G_n = e^{−nκh}, κh = 2·arsinh(h/2λ); Abweichung von s: 1 − κλ ≈ (h/λ)²/24, h/λ = s·r₀/r*")
# Prüfung an endlicher Leiter
hl=0.3; Nn=200; A=np.zeros((Nn,Nn)); b=np.zeros(Nn)
for i in range(Nn):
    A[i,i]=-(2+hl**2)
    if i>0: A[i,i-1]=1
    if i<Nn-1: A[i,i+1]=1
b[0]=-1.0; Gv=np.linalg.solve(A,b)
P(f"   Prüfung h/λ = {hl}: G_10/G_9 = {Gv[10]/Gv[9]:.10f}, e^(−2 arsinh(h/2λ)) = {np.exp(-2*np.arcsinh(hl/2)):.10f}, Kontinuum e^(−h/λ) = {np.exp(-hl):.10f}")
for nm,M in (("Sonne",Msun),("Milchstraße",Mmw),("Haufen 1,7e14 M☉",1.7e14*Msun),("beob. Universum 7,5e22 M☉",7.5e22*Msun)):
    Ds=G*M/(3*c**2); N=c/np.sqrt(12*a0*Ds); hl=1/N
    corr=float(1-(2/mpm.mpf(hl))*mpm.asinh(mpm.mpf(hl)/2))
    P(f"   {nm:26s}: r*/r₀ = {N:.3e},  Abweichung von s bei s = 1: {corr:.2e}  ((r₀/r*)²/24 = {hl**2/24:.2e})")
P()

# ================= 4) Umfang zu Fläche: eine Stufe höher und tiefer =================
P("4a) Rückkehrrichtung bei Λ geschlossen (G'(Λ) = 0): q_Λ = cosh((Λ−L)/λ)/cosh(Λ/λ), L/λ = s = r*/r (Punktmasse, L = r*)")
def omq(s,LL):  # 1 − q_Λ
    return 1-np.cosh(s)+np.tanh(s*LL)*np.sinh(s)
LL=100.0
lnA=lambda r: np.log((1/r**2)/omq(1/r,LL))
P(f"   Λ/L = {LL:g}:")
P("   r/r*      a/a₀         d_eff     offen: d_eff")
for x in (0.3,1,10,30,100,300,1000,1e4):
    de=1-(lnA(x*(1+h))-lnA(x*(1-h)))/(np.log1p(h)-np.log1p(-h))
    P(f"   {x:<8g}  {np.exp(lnA(x)):.4e}   {de:.4f}    {3-chi(1/x):.4f}")
xr=brentq(lambda r: (1/omq(1/r,LL))/nu(1/r)-1.1,1,1e4)
P(f"   10 % über der offenen Regel ab r = {xr:.2f} r* = {xr/LL:.3f}·r*·Λ/L;  Boden a_∞ = a₀/(Λ/L − ½) = {1/(LL-0.5):.4e} a₀")
for LL2 in (4.6e5,4e7):
    xr2=brentq(lambda r: (1/omq(1/r,LL2))/nu(1/r)-1.1,1,1e10)
    P(f"   Λ/L = {LL2:.1e}: 10 % ab {xr2/LL2:.3f}·r*·Λ/L (Milchstraße {xr2*rsmw/kpc/1e6:.2f} Gpc), Boden {1/(LL2-0.5):.2e} a₀")
from scipy.special import gamma as Gf
P("   Flussfläche S_{d−1} = 2π^{d/2}/Γ(d/2): " + ", ".join(f"d = {d}: {2*pi**(d/2)/Gf(d/2):.4f}" for d in (3,2,1)) + "  (4π, 2π, 2: bei d = 1 kein π)")
P()
P("4b) Quellen im Volumen: a ∝ r, d_eff = 1 − 1 = 0 (kein Flussgesetz)")
P("   Homogene Kugel innen: a = GM r/R³;  Vakuumteil außen: a = +Λ_s c² r/3 (nach außen).")
Rp=0.84e-15
P(f"   Proton: g(R)/a₀ = {G*mp/Rp**2/a0:.0f} (kompakter Trichter, Newtonsch bis zur Glättungslänge); kein eigenes Umfanggesetz auf Teilchenebene.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"orthogonalitaet_wurzeln.txt"),"w",encoding="utf-8").write(txt); print(txt)
