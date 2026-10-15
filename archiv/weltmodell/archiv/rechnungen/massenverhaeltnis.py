# m_p/m_e und eine Zahl Im τ_X für α und α_G: nur strukturelle Wege des Modells, vorab eine Zufallsprobe.
# 1 Zufallsprobe: wie oft treffen einfache Ausdrücke aus Zahlen des Modells (π, kleine ganze Zahlen, Potenzen, Wurzeln) die Ziele zufällig?
# 2 Kopplungen als Längenverhältnisse, Δs·λ̄ = ℓ_P²/3, a_J·Δs = ℓ_P²/6: kürzt sich die Masse?
# 3 Carter-Kette (Exponenten 12 und 4 bzw. 5): was sie über m_p/m_e festlegt
# 4 Nome-Darstellung α_G = e^{−2π Im τ}: Werte von Im τ und ob eine Zahl beide Kopplungen trägt
import numpy as np, sys, io, os, itertools, sympy as sy
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out=[]; P=lambda s="": out.append(s)
rng=np.random.default_rng(7)

# Messwerte [G] (CODATA, aus dem Gedächtnis): relative Unsicherheit in Klammern
MU=1836.15267343; dMU=6e-11            # m_p/m_e
AINV=137.035999177; dA=1.5e-10         # 1/α
AG=5.906e-39                            # α_G = G m_p²/(ħc), Unsicherheit ~2·10⁻⁵ (aus G)
P(f"Ziele [G]: m_p/m_e = {MU} (±{dMU:.0e}), 1/α = {AINV} (±{dA:.1e}), α_G = {AG}")
P()

# ---------------- 1 Zufallsprobe ----------------
L2,L3,L5,Lp=np.log(2),np.log(3),np.log(5),np.log(np.pi)
def classA():                    # 2^a 3^b 5^c π^d, Exponenten in (1/n)·[−12,12], n = 1..4
    e=np.arange(-12,13); A,B,C,D=np.meshgrid(e,e,e,e,indexing="ij")
    base=(A*L2+B*L3+C*L5+D*Lp).ravel()
    return np.unique(np.round(np.concatenate([base/n for n in (1,2,3,4)]),12))
def classB():                    # Σ_{i=0}^{5} c_i π^i, c_i ∈ [−4,4]
    c=np.arange(-4,5); pw=np.pi**np.arange(6)
    vals=np.array(list(itertools.product(c,repeat=6)))@pw
    vals=vals[vals>0]; return np.unique(np.round(np.log(vals),12))
def classC():                    # (p/q)·π^k, p,q ≤ 64, k ∈ [−6,6]
    p=np.arange(1,65); Pq,Qq,K=np.meshgrid(p,p,np.arange(-6,7),indexing="ij")
    return np.unique(np.round((np.log(Pq)-np.log(Qq)+K*Lp).ravel(),12))
CL={"A: 2^a·3^b·5^c·π^d, Wurzeln bis 4":classA(),"B: Σ c_i·π^i (i ≤ 5, |c_i| ≤ 4)":classB(),"C: (p/q)·π^k (p,q ≤ 64)":classC()}

def best(lv, T):
    lt=np.log(T); j=np.searchsorted(lv,lt); cand=lv[max(j-1,0):j+1]
    k=np.argmin(abs(cand-lt)); return abs(np.expm1(cand[k]-lt)), cand[k]
def count(lv, T, d): lt=np.log(T); return int(np.sum(abs(lv-lt)<d))
def density(lv, T, w=0.2):
    lt=np.log(T); return np.sum(abs(lv-lt)<w)/(2*w)          # Werte je Einheit in ln x

P("1) Zufallsprobe. δ = |x/Ziel − 1|; erwartet N(δ) = Dichte·2δ; p = Anteil zufälliger Ziele (log-gleichverteilt, ×0,67…1,5) mit gleich gutem Treffer")
REC={}
for T,name in ((MU,"m_p/m_e"),(AINV,"1/α")):
    P(f"   Ziel {name} = {T}")
    P("   Klasse                                  Werte      Dichte/ln   bester δ    Ausdruck(Wert)     N(10⁻⁴) erw.   N(10⁻⁵) erw.   N(10⁻⁶) erw.   p(bester)")
    for cn,lv in CL.items():
        db,vb=best(lv,T); rho=density(lv,T)
        tg=np.exp(np.log(T)+rng.uniform(-np.log(1.5),np.log(1.5),4000))
        bd=np.array([best(lv,t)[0] for t in tg]); pval=np.mean(bd<=db)
        cells=[]
        for d in (1e-4,1e-5,1e-6):
            cells.append(f"{count(lv,T,d):3d} / {rho*2*d:6.2f}")
        P(f"   {cn:38s} {len(lv):9d}   {rho:9.0f}   {db:.1e}   {np.exp(vb):.6f}   " + "   ".join(cells) + f"   {pval:.2f}")
        REC[(name,cn)]=(db,pval,rho)
    P(f"   Auf die Messgenauigkeit (δ ≈ {dMU if name=='m_p/m_e' else dA:.0e}) erwartet: "
      + ", ".join(f"{density(lv,T)*2*(dMU if name=='m_p/m_e' else dA):.1e}" for lv in CL.values()) + " Treffer je Klasse; gefunden: "
      + ", ".join(str(count(lv,T,(dMU if name=='m_p/m_e' else dA))) for lv in CL.values()))
    P()
P(f"   Bekannte Beispiele: 6π⁵ = {6*np.pi**5:.4f} (δ = {abs(6*np.pi**5/MU-1):.1e});  4π³ + π² + π = {4*np.pi**3+np.pi**2+np.pi:.5f} (δ = {abs((4*np.pi**3+np.pi**2+np.pi)/AINV-1):.1e})")
pm=min(v[1] for v in REC.values()); P(f"   Kleinstes p aller 6 Prüfungen (3 Klassen × 2 Ziele): {pm:.2f}; Wahrscheinlichkeit, dass mindestens eine Prüfung so klein ausfällt: {1-(1-pm)**6:.2f}")
P("   → Beide liegen in ihrer Klasse im Bereich zufälliger Treffer (p in der Tabelle). Eine Beziehung zählt erst, wenn sie aus einer Gleichung")
P("     des Modells folgt und weit unter dem zufällig erreichbaren δ liegt; auf Messgenauigkeit ist in keiner Klasse ein Zufallstreffer zu erwarten.")
P()

# ---------------- 2 Längenverhältnisse ----------------
P("2) Kopplungen als Längenverhältnisse (symbolisch)")
m,G,hb,c,e0=sy.symbols('m G hbar c e0',positive=True)
lam=hb/(m*c); ds=G*m/(3*c**2); aJ=hb/(2*m*c); lP2=G*hb/c**3
P(f"   Δs·λ̄ = {sy.simplify(ds*lam)} = ℓ_P²/3: Masse kürzt sich ({sy.simplify(ds*lam/lP2)}·ℓ_P²)")
P(f"   a_J·Δs (Spin ½, a_J = ħ/2mc) = {sy.simplify(aJ*ds)} = {sy.simplify(aJ*ds/lP2)}·ℓ_P²: Masse kürzt sich")
mp,me=sy.symbols('m_p m_e',positive=True)
P(f"   α_G = 3Δs_p/λ̄_p = {sy.simplify(3*ds.subs(m,mp)/lam.subs(m,mp))};  m_p/m_e = λ̄_e/λ̄_p = {sy.simplify(lam.subs(m,me)/lam.subs(m,mp))}")
P("   → Die Beziehungen sind für jede Masse erfüllt; sie verbinden nur eine Masse mit ihrer eigenen Länge.")
P("     Zwischen α, α_G und m_p/m_e liefern sie keine Gleichung (drei Zahlen, null Beziehungen).")
P()

# ---------------- 3 Carter ----------------
P("3) Carter-Kette: α_G = 𝒞·α¹²·(m_e/m_p)⁴ (Thomson) bzw. 𝒞'·α¹²·(m_e/m_p)⁵ (Kramers); 𝒞, 𝒞' nicht hergeleitet")
a=1/AINV
Cc=AG/(a**12*MU**-4); Ck=AG/(a**12*MU**-5)
P(f"   Mit den Messwerten: 𝒞 = {Cc:.4f}, 𝒞' = {Ck:.0f}")
for CC in (1.0,Cc):
    P(f"   𝒞 = {CC:.3f}:  m_p/m_e = (𝒞α¹²/α_G)^{{1/4}} = {(CC*a**12/AG)**0.25:.2f};  1/α aus α_G und m_p/m_e = {1/((AG*MU**4/CC)**(1/12)):.3f}")
P(f"   Kramers, 𝒞' = 1: m_p/m_e = (α¹²/α_G)^{{1/5}} = {(a**12/AG)**0.2:.2f}")
P("   Empfindlichkeit: δμ/μ = ¼·δ𝒞/𝒞. Für m_p/m_e auf 10⁻¹⁰ müsste 𝒞 auf 4·10⁻¹⁰ bekannt sein; 𝒞 enthält Saha- und Gamow-Logarithmen")
P("   und die mittlere Opazität der Sterne. Die Kette legt Exponenten fest, keinen Zahlenwert von m_p/m_e.")
P()

# ---------------- 4 Nome ----------------
P("4) Nome-Darstellung |q| = e^{−2π Im τ}")
ImG=-np.log(AG)/(2*np.pi); ImA=np.log(AINV)/(2*np.pi)
lamp=1.32141e-15/(2*np.pi); ImP=np.log(3/AG)/(2*np.pi)
ImE=ImP+np.log(MU)/np.pi
P(f"   α_G = e^{{−2π Im τ_G}}: Im τ_G = {ImG:.4f};   α_G = |q|² (Kern-Nome der Kopplung): Im τ = {ImG/2:.4f}")
P(f"   α = e^{{−2π Im τ_α}}: Im τ_α = {ImA:.4f};   Kern-Nome Δs/λ̄ an der eigenen Quantenlänge: Proton Im τ_p = {ImP:.4f}, Elektron Im τ_e = {ImE:.4f}")
P(f"   m_p/m_e = e^{{π(Im τ_e − Im τ_p)}} = {np.exp(np.pi*(ImE-ImP)):.5f}: Umschreibung, keine Bestimmung")
P(f"   Carter in Nome-Form: Im τ_G = 12·Im τ_α + (2/π)·ln(m_p/m_e) − ln𝒞/2π = {12*ImA:.4f} + {2/np.pi*np.log(MU):.4f} − {np.log(Cc)/(2*np.pi):.4f} = {12*ImA+2/np.pi*np.log(MU)-np.log(Cc)/(2*np.pi):.4f}")
P("   Eine Zahl Im τ_X für α und α_G verlangte eine Gleichung α = f(Im τ_X) neben α_G = e^{−2π Im τ_X}; das Modell hat keine.")
P("   Selbst dann bliebe m_p/m_e über Carter an 𝒞 gebunden (μ = exp[(π/2)(Im τ_G − 12 Im τ_α)]·𝒞^{1/4}).")
# Nähe von Im τ_G zu einer ganzen Zahl: Zufallsprobe
fr=abs(ImG-round(ImG)); ratio=ImG/ImA
P(f"   Im τ_G liegt {fr:.4f} neben der ganzen Zahl {round(ImG)} (relativ {fr/ImG:.1e}). Zufall: Abstand ≤ {fr:.4f} zur nächsten ganzen Zahl mit"
  f" Wahrscheinlichkeit {2*fr:.3f}; geprüft wurden hier mindestens vier Im-τ-Werte ({ImG:.3f}, {ImG/2:.3f}, {ImP:.3f}, {ImE:.3f}), also ≈ {1-(1-2*fr)**4:.2f}.")
P(f"   ln α_G/ln α = Im τ_G/Im τ_α = {ratio:.4f} (keine Gleichung des Modells liefert diese Zahl)")
lvA=CL["A: 2^a·3^b·5^c·π^d, Wurzeln bis 4"]
for T,nm in ((ImG,"Im τ_G"),(ratio,"Im τ_G/Im τ_α"),(AINV/ImG,"(1/α)/Im τ_G")):
    db,vb=best(lvA,T); tg=np.exp(np.log(T)+rng.uniform(-np.log(1.5),np.log(1.5),3000))
    pv=np.mean(np.array([best(lvA,t)[0] for t in tg])<=db)
    P(f"   {nm} = {T:.5f}: bester Ausdruck der Klasse A {np.exp(vb):.5f} (δ = {db:.1e}), p = {pv:.2f}")
P()
P("Ergebnis: Keine Gleichung des Modells bestimmt m_p/m_e oder verbindet α und α_G über eine Zahl; alle Treffer liegen im Zufallsbereich.")

txt="\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"massenverhaeltnis.txt"),"w",encoding="utf-8").write(txt); print(txt)
