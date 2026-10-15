# Umlaufsinn der Rückkehr: warum a = a_N/(1 − q) nur Rückläufe in einem Sinn zählt (Σ_{m≥0} qᵐ statt Σ_{m∈ℤ} q^{|m|}).
# 0 Vergleich   1 Kausalität (Operator, Gerade/Ring mit Verlust Γ, Kramers-Kronig, Aufbau in Echtzeit)
# 2 KMS (Bildsumme im Kreis imaginärer Zeit)   3 Orientierung des Torus (Reihe in q gegen q und q̄)
# 4 Drehsinn (Kern des Autors, Fall Zeit gegen Fall Scheibe, druckgestützte Systeme, Dissipation und Casimir-Anteil)
import numpy as np, sys, io, os
from scipy.integrate import quad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G=6.674e-11; c=2.99792458e8; a0=1.2e-10; Msun=1.989e30; kpc=3.0857e19; Gyr=3.156e16; tU=13.8*Gyr
out=[]; P=lambda s="": out.append(s)
nu=lambda s: 1/(-np.expm1(-s)); coth=lambda x: 1/np.tanh(x)

# ---------------- 0) Vergleich ----------------
P("0) ν = a/a_N, s = √(a_N/a₀), q = e^{−s}")
P("   a_N/a₀   Σ_{m≥0}qᵐ   Σ_{m∈ℤ}q^|m| = coth(s/2)   ½(1 + coth(s/2))   coth s (gleicher TF-Nullpunkt)   1/(1 − 1/cosh s)")
for x in (10,3,1,0.3,0.1,0.01):
    s=np.sqrt(x)
    P(f"   {x:6g}   {nu(s):9.4f}   {coth(s/2):9.4f}   {0.5*(1+coth(s/2)):9.4f}   {coth(s):9.4f}   {1/(1-1/np.cosh(s)):9.4f}")
P("   Tiefgrenze: 1/s + ½ (v⁴ = GMa₀) | 2/s (v⁴ = 4GMa₀) | = einseitig | 1/s | 2/s² (a → 2a₀)")
P("   Identität: 1/(1 − q) = ½ + ½·coth(s/2): ein Sinn mit vollem Gewicht = beide Sinne mit halbem Gewicht.")
P()

# ---------------- 1) Kausalität ----------------
s=1.0; q=np.exp(-s)
P("1a) Operator: T = Verschiebung um einen Rücklauf, Symbol z (|z| = 1 ↔ reelle Frequenz, z = e^{iωτ_R}). s = 1")
N=4096; z=np.exp(2j*np.pi*np.arange(N)/N)
R=1/(1-q*z); S=(1-q*q)/((1-q*z)*(1-q/z))
cR=np.fft.fft(R)/N; cS=np.fft.fft(S)/N; m=np.arange(1,21)
P(f"   (I − qT)⁻¹ = Σ_{{m≥0}}: Koeffizienten m ≥ 0: max|c_m − qᵐ| = {max(abs(cR[:21]-q**np.arange(21))):.1e};  m < 0: max|c_m| = {max(abs(cR[-m])):.1e}")
P(f"   Pol bei z = e^s = {np.exp(s):.3f} (außerhalb |z| ≤ 1): analytisch in der Scheibe ⇔ Träger m ≥ 0 ⇔ in ω analytisch für Im ω > 0")
P(f"   (1 − q²)(I − qT)⁻¹(I − qT⁻¹)⁻¹ = Σ_{{m∈ℤ}} q^|m|Tᵐ: max|c_m − q^|m|| = {max(max(abs(cS[:21]-q**np.arange(21))),max(abs(cS[-m]-q**m))):.1e};"
  f" Pole bei z = e^{{±s}}: beide Seiten")
P("   → einseitig = Green-Funktion des Operators erster Ordnung I − qT (kausaler Wiener-Hopf-Faktor);")
P("     beidseitig = Green-Funktion des symmetrischen Operators zweiter Ordnung (1+q²) − q(T + T⁻¹), diskret −D∂² + Γ")
P(f"   Die andere formale Inverse von I − qT, −Σ_{{m≥1}}(qT)^{{−m}}, divergiert: q^{{−50}} = {q**-50:.1e}")
P(f"   Gerade/ungerade: max|Re R(z) − ½(1 + S(z))| = {max(abs(R.real-0.5*(1+S.real))):.1e};  ungerader Teil Im R bei z = 1: {R[0].imag:.1e}")
P(f"   Statisch trägt nur der gerade (zeitumkehrsymmetrische) Teil: R(1) = {R[0].real:.4f} = ½(1 + S(1)) = {0.5*(1+S[0].real):.4f}")
P()

P("1b) Explizites Modell: Diffusion mit Verlust Γ auf der Geraden (Überlagerung des Rings). Einheiten 1/Γ, L²/D = s², s = 1")
P("    Rücklauf m in einem Sinn = erstes Erreichen von +mL: f_m(t) = m s/√(4πt³)·exp(−m²s²/(4t) − t), t ≥ 0")
fm=lambda t,mm,s: 0.0 if t<=0 else mm*s/np.sqrt(4*np.pi*t**3)*np.exp(-mm*mm*s*s/(4*t)-t)
P("    " + ",  ".join(f"∫f_{mm} = {quad(fm,0,np.inf,args=(mm,s))[0]:.6f} (e^{{−{mm}s}} = {np.exp(-mm*s):.6f})" for mm in (1,2,3)))
rt=lambda t: sum(fm(t,mm,s) for mm in range(1,60))
Rw=lambda w: 1/(1-np.exp(-s*np.sqrt(1-1j*np.asarray(w,dtype=complex))))
P(f"    Folgedichte r(t) = Σ_{{m≥1}} f_m, Träger t ≥ 0, ∫r = {quad(rt,0,np.inf,limit=200)[0]:.5f} = q/(1 − q) = {q/(1-q):.5f}")
for w in (0.5,2.0):
    re=quad(rt,0,np.inf,weight='cos',wvar=w)[0]; im=quad(rt,0,np.inf,weight='sin',wvar=w)[0]
    P(f"    ω = {w}Γ:  ∫r e^{{iωt}} = {re:.6f} + {im:.6f}i;  geschlossen R(ω) − 1 = {(Rw(w)-1).real:.6f} + {(Rw(w)-1).imag:.6f}i,  R(ω) = 1/(1 − e^{{−s√(1 − iω/Γ)}})")
X,Y=np.meshgrid(np.linspace(-200,200,801),np.linspace(0,20,201)); W=X+1j*Y
P(f"    Oberhalb der Achse (0 ≤ Im ω ≤ 20Γ, |Re ω| ≤ 200Γ): min|1 − q(ω)| = {abs(1-np.exp(-s*np.sqrt(1-1j*W))).min():.3f} > 0;"
  f" Singularitäten bei ω = −iΓ(1 + 4π²k²/s²), alle unten")
KK=lambda F,w0: quad(lambda w: F(w).imag,-3000,3000,weight='cauchy',wvar=w0,limit=2000)[0]/np.pi
R2=lambda w: 1+2*(Rw(w)-1)                      # beide Windungssinne auf dem Ring (Bildsumme), weiterhin t ≥ 0
Ssym=lambda w: 1+2*(Rw(w)-1).real+0j            # zeitsymmetrisch: r(|t|)
for name,F in (("ein Sinn, retardiert      ",Rw),("Ring, beide Sinne (Bilder)",R2),("zeitsymmetrisch r(|t|)    ",Ssym)):
    P(f"    {name}: statisch {F(0.0).real:.4f};  Kramers-Kronig Re F − 1 gegen (1/π)P∫Im F/(ω−ω₀):"
      f"  ω₀=0: {F(0.0).real-1:.5f} | {KK(F,0.0):.5f};  ω₀=Γ: {F(1.0).real-1:.5f} | {KK(F,1.0):.5f}")
for e in (0.99,0.9999):
    P(f"    ω = {e}iΓ:  |R| = {abs(Rw(1j*e)):.3f},  |zeitsymmetrisch| = |R(ω) + R(−ω) − 1| = {abs(Rw(1j*e)+Rw(-1j*e)-1):.1f}  (Pol bei +iΓ)")
P("    → Die Windungssumme beider Sinne auf dem Ring ist ebenso retardiert wie die eines Sinnes; Kausalität trennt die Sinne nicht.")
P("      Zeitsymmetrisch und Ring (beide Sinne) haben denselben statischen Wert coth(s/2).")
P()

P("1c) Aufbau in Echtzeit mit Γ = a_N/c, τ_D = L²/D = c/a₀ (Weltalter 13,8 Mrd. J. als Eingabe)")
tD=c/a0; P(f"    τ_D = {tD/Gyr:.1f} Mrd. J.;  t_U/τ_D = {tU/tD:.3f};  mittlere Laufzeit des m-ten Rücklaufs m·τ_D/(2s)")
for x in (1,0.1,0.01):
    s_=np.sqrt(x); T=tU/tD*s_*s_        # t_U in Einheiten 1/Γ
    arr=1+sum(quad(fm,0,T,args=(mm,s_),limit=200)[0] for mm in range(1,int(40/s_)+5))
    P(f"    a_N/a₀ = {x:5g}: bis t_U eingetroffen ν(t_U) = {arr:.4f}  gegen stationär ν = {nu(s_):.4f};  erster Rücklauf nach {tD/(2*s_)/Gyr:.0f} Mrd. J. im Mittel")
P("    → Als Vorgang in Echtzeit wäre die Summe heute nicht aufgebaut; sie muss stationär sein, und dann wählt kein Zeitpfeil einen Umlaufsinn.")
P("      Die Laufzeitdarstellung (D∂² − Γ)⁻¹ = ∫₀^∞ e^{t(D∂²−Γ)} dt ist stets einseitig in t (konvergiert für Γ > 0) und enthält doch beide Raumsinne.")
P()

P("1d) Stationär, finite Differenzen (h = L/400, s = 1): D·G'' = Γ·G auf der Halbgeraden ab dem Ort gegen Quelle in der vollen Geraden")
from scipy.linalg import solve_banded
def fdG(Wl,src_mid,s=1.0,nL=400):
    n=int(Wl*nL); h=1.0/nL; ab=np.zeros((3,n)); ab[0,1:]=1; ab[2,:-1]=1; ab[1,:]=-2-(s*h)**2; b=np.zeros(n)
    if src_mid: b[n//2]=-h                       # Punktquelle in der Mitte, Ränder G = 0 weit weg
    else: b[0]=-1                                # G(0) = 1 am Ort, G(W) = 0
    g=solve_banded((1,1),ab,b); return g,n
g,n=fdG(40,False); gm=[1.0]+[g[m*400-1] for m in range(1,6)]
g2,n2=fdG(80,True); c0=n2//2; gp=[g2[c0+m*400]/g2[c0] for m in range(0,6)]; gn=[g2[c0-m*400]/g2[c0] for m in range(1,6)]
P(f"    Halbgerade, G(0) = 1, beschränkt: Σ_{{m≥0}} G(mL) = {sum(gm)+sum(np.exp(-np.arange(6,60)*1.0)):.4f} (G(L) = {gm[1]:.5f}, e^{{−s}} = {np.exp(-1):.5f});"
  f" nur m ≥ 0 existiert, (∂ + √(Γ/D))G = 0")
P(f"    volle Gerade, Quelle am Ort: G(±L)/G(0) = {gp[1]:.5f} / {gn[0]:.5f};  Σ_{{m∈ℤ}} G(mL)/G(0) = {1+2*(sum(gp[1:])+sum(np.exp(-np.arange(6,60)*1.0))):.4f} = coth(½) = {coth(0.5):.4f}")
P("    → Der eine Sinn ist stationär gleichwertig mit einer Halbgeraden der inneren Richtung, die am Ort beginnt (Wiener-Hopf-Faktor).")
P()

# ---------------- 2) KMS ----------------
P("2) Bildsumme im Kreis imaginärer Zeit (Oszillator ω = 1, Umfang β = x). G_β(τ) = Σ_{n∈ℤ} e^{−|τ+nβ|}/2")
for x in (1.0,0.316):
    tau=0.3*x; nn=np.arange(-400,401)
    img=np.exp(-abs(tau+nn*x)).sum()/2; closed=np.cosh(x/2-tau)/(2*np.sinh(x/2))
    Ap=np.exp(-np.arange(0,400)*x).sum(); Am=np.exp(-np.arange(1,400)*x).sum()
    k=np.arange(0,4000); nB=(k*np.exp(-k*x)).sum()/np.exp(-k*x).sum()
    chiB=quad(lambda t: np.cosh(x/2-t)/(2*np.sinh(x/2)),0,x)[0]
    P(f"   x = {x}: Bildsumme {img:.6f} = geschlossen {closed:.6f} = [(Σ_{{n≥0}})e^{{−τ}} + (Σ_{{n≥1}})e^{{τ}}]/2 = {(Ap*np.exp(-tau)+Am*np.exp(tau))/2:.6f}")
    P(f"      Σ_{{n≥0}} = {Ap:.4f} = 1 + n_B,  Σ_{{n≥1}} = {Am:.4f} = n_B (Spur: {nB:.4f});  Verhältnis {Ap/Am:.4f} = e^x = {np.exp(x):.4f} (KMS)")
    P(f"      Propagator gleichzeitig 2G_β(0) = {np.exp(-abs(nn*x)).sum():.4f} = coth(x/2) = 1 + 2n_B;  Kommutator G^> − G^< = {Ap-Am:.4f};"
      f"  statische Suszeptibilität ∫₀^β G_β = {chiB:.4f} (unabhängig von β)")
P("   → Der Propagator ist beidseitig (coth), die Besetzung einseitig. + Windungen gehören zur positiven, − Windungen zur negativen Frequenz.")
P("     ν = 1/(1 − q) = 1 + n_B(s) = G^>/ρ (Emissionsgewicht einer Ordnung); Zusatz a − a_N = a_N·n_B(s), χ = s·n_B(s).")
P("     Die retardierte Antwort (Kommutator) trägt keinen Bose-Faktor; beide Ordnungen G^> + G^< geben coth(s/2); KMS: Gegensinn = q × Sinn.")
P()

# ---------------- 3) Orientierung des Torus ----------------
P("3) Torus: q = e^{2πiτ}, Im τ = s/2π. Holomorph f = Σ_{n≥0}qⁿ gegen reell-analytisch Σ_{n≥0}qⁿ + Σ_{n≥1}q̄ⁿ = (1 − |q|²)/|1 − q|²")
s=1.0; y=s/(2*np.pi); Nn=256; u=np.arange(Nn)/Nn; qq=np.exp(2j*np.pi*(u+1j*y))
f=1/(1-qq); Pk=(1-abs(qq)**2)/abs(1-qq)**2; cf=np.fft.fft(f)/Nn; cP=np.fft.fft(Pk)/Nn; n=np.arange(1,9)
kk=np.fft.fftfreq(Nn,1/Nn); hil=lambda g: np.fft.ifft(-1j*np.sign(kk)*np.fft.fft(g)).real
P(f"   Fourier in Re τ: f: n ≥ 0 max|c_n − e^{{−ns}}| = {max(abs(cf[n]-np.exp(-n*s))):.1e}, n < 0 max|c_n| = {max(abs(cf[-n])):.1e};"
  f"  reell-analytisch: c_{{±n}} = e^{{−|n|s}} (Fehler {max(abs(cP[-n]-np.exp(-n*s))):.1e})")
P(f"   Kramers-Kronig längs Re τ: max|Im f − H[Re f]| = {max(abs(f.imag-hil(f.real))):.1e};  reell-analytisch: Im = 0, max|H[P]| = {max(abs(hil(Pk))):.3f}")
P("   → Holomorph in Im τ > 0 ⇔ Träger n ≥ 0 ⇔ Kramers-Kronig (Titchmarsh): dieselbe Struktur wie Kausalität, mit der Windungszahl als Zeit.")
P("   Auf der Achse Re τ = 0 (Bindung nur über |q|):  a_N/a₀   f   Re f = ½(f + f̄)   Σq + Σq̄ = coth(s/2)   dasselbe mit q² = coth s   |f|²")
for x in (3,1,0.3,0.1,0.01):
    s_=np.sqrt(x); fq=nu(s_); P(f"      {x:6g}   {fq:8.4f}   {fq:8.4f}   {coth(s_/2):8.4f}   {coth(s_):8.4f}   {fq*fq:9.4f}")
qq2=np.exp(2j*np.pi*(0.1+1j*y))
P(f"   Umkehr der Orientierung τ → −τ̄ (q → q̄): |f(q)| = {abs(1/(1-qq2)):.6f} = |f(q̄)| = {abs(1/(1-np.conj(qq2))):.6f} bei Re τ = 0,1; auf der Achse gleich.")
P("   → Die Orientierung selbst ist in a unsichtbar; sichtbar ist nur, ob der zweite Sinn mit vollem Gewicht hinzukommt (coth).")
P()

# ---------------- 4) Drehsinn ----------------
P("4a) Kern des Autors auf der Kreisbahn, beide Drehsinne (κ₀ = r = 1, λ_m = 0,3): a = −κ₀∫₀^∞ e^{−λτ} Π v(t−τ) dτ")
lam=0.3
for om in (1.0,-1.0):
    ev=-quad(lambda t: np.exp(-lam*t)*om*np.sin(om*t),0,np.inf,limit=400)[0]          # Π = e_r⊗e_r
    od=-quad(lambda t: np.exp(-lam*t)*om*np.cos(om*t),0,np.inf,limit=400)[0]          # Π = e_r⊗e_φ − e_φ⊗e_r (radialer Teil)
    P(f"    ω = {om:+.0f}: e_r⊗e_r: a_r = {ev:.5f} (−ω²/(λ²+ω²) = {-om*om/(lam**2+om**2):.5f});  antisymmetrisch: a_r = {od:.5f} (−λω/(λ²+ω²))")
P("    → Der radiale Projektor gibt einen in ω geraden Term: gleich für mit- und gegenläufige Bahnen. Einen ungeraden Radialanteil")
P("      gäbe nur ein antisymmetrischer Projektor (Casimir-Typ, gekoppelt an den Drehsinn); er fehlt in der Gleichung.")
P()
P("4b) Zwei Fälle, Punktmasse (a_N = a₀(r*/r)²). Z: Sinn = Zeit (universell). D: Sinn = Drehsinn der Scheibe")
P("    D1: nur Rückläufe im Drehsinn der Scheibe zählen, Gegenläufer Newton;  D2: getragen vom Nettodrehimpuls, ν = 1 + λ·q/(1 − q), λ = |L_netto|/Σ|L|")
P("    r/r*    ν (Z, alle)   v_gegen/v_mit (D1) = v(λ=0)/v(λ=1) (D2)   v(λ=0,4)/v(λ=1) (D2, Scheiben 70:30)")
for rr in (0.5,1,2,3,5,10):
    s_=1/rr; nv=nu(s_); nl=1+0.4*(nv-1)
    P(f"    {rr:5g}   {nv:8.3f}   {1/np.sqrt(nv):8.3f}   {np.sqrt(nl/nv):8.3f}")
P("    Tiefgrenze D2: v⁴ = λ²·GMa₀ (Tully-Fisher-Nullpunkt hängt am Drehanteil); Z: v⁴ = GMa₀ für jedes λ.")
P()
P("4c) Druckgestützt (Plummer, isotrop): Virial ⟨v²⟩ = −(1/M)∫ρ r·a dV, σ₁ = √(⟨v²⟩/3). Z gegen D (λ = 0 → Newton)")
for Mb,b in ((1e7,0.3),(1e9,1.0),(1e11,3.0)):
    GM=G*Mb*Msun; bb=b*kpc
    aN=lambda x: GM/bb**2*x/(1+x*x)**1.5                       # x = r/b
    w=lambda x: 3*x*x*(1+x*x)**-2.5                            # Massengewicht dM/M pro dx
    vz=quad(lambda x: w(x)*bb*x*aN(x)*nu(np.sqrt(aN(x)/a0)),0,np.inf,limit=400)[0]
    vn=quad(lambda x: w(x)*bb*x*aN(x),0,np.inf,limit=400)[0]
    P(f"    M = {Mb:.0e} M☉, b = {b} kpc: σ₁(Z) = {np.sqrt(vz/3)/1e3:6.2f} km/s,  σ₁(Newton = D, λ = 0) = {np.sqrt(vn/3)/1e3:6.2f} km/s"
      f" (3π/32·GM/b: {np.sqrt(3*np.pi/32*GM/bb/3)/1e3:.2f}),  Verhältnis {np.sqrt(vz/vn):.2f};  Tiefgrenze Z √((2/9)√(GMa₀)) = {np.sqrt(2/9*np.sqrt(GM*a0))/1e3:.2f} km/s")
P()
P("4d) Resolvente des Modells 1b (s = 1): Re R gerade in ω (reaktiv), Im R ungerade (absorptiv, dissipativ)")
for w in (0.0,0.3,1.0,3.0):
    P(f"    ω = {w:3.1f}Γ: Re R = {Rw(w).real:.4f}, Im R(+ω) = {Rw(w).imag:+.4f}, Im R(−ω) = {Rw(-w).imag:+.4f}")
for g in (1.0,1e-2,1e-4):
    P(f"    Verlust Γ → {g:g}·Γ (s ∝ √Γ): q = {np.exp(-np.sqrt(g)):.4f}, ν = {nu(np.sqrt(g)):.1f}")
P(f"    Summenregel (Kramers-Kronig bei ω₀ = 0): R(0) − 1 = {Rw(0.0).real-1:.5f} = (1/π)P∫Im R/ω dω = {KK(Rw,0.0):.5f}; Ring: Im R₂ = 2 Im R → {KK(R2,0.0):.5f};"
  f" zeitsymmetrisch: Im = 0 → {KK(Ssym,0.0):.5f}")
P("    → Bei Kausalität ist die statische Verstärkung das Integral des absorptiven Spektrums; ein Sinn hat Im R, zwei Sinne 2·Im R.")
P("      Γ > 0 (dS/dt ≥ 0) macht die Reihe konvergent (Γ → 0: ν → ∞) und wählt einen Zeitsinn, wenn T eine Zeitverschiebung ist.")
P("      Einen Raumsinn kann Γ nicht wählen (Diffusion ist spiegelsymmetrisch); das könnte nur ein Casimir-Anteil mit dem Drehsinn Ω, der fehlt (4a).")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"umlaufsinn.txt"),"w",encoding="utf-8").write(txt); print(txt)
