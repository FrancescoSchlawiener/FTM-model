# B1: Ω_b/Ω_m aus Kac für Flüsse: Anteil = Verweildauer/Zyklus = (1/H)/(2π/H)
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
c=2.99792458e8; Mpc=3.0857e22
out=[]; P=lambda s="": out.append(s)
nu=2*np.pi; q=1-1/nu; Q=nu-1
P("1) Kac für Flüsse: μ(A) = E[Verweildauer]/E[Zyklus]")
P(f"   Verweildauer 1/H, Zyklus 2π/H  →  Ω_b/Ω_m = 1/2π = {1/nu:.5f},  ν = {nu:.5f},  q = {q:.5f},  Q = {Q:.5f}")
P("   H kürzt sich: das Verhältnis ist zeitlich konstant, auch wenn H atmet.")
P()
P("2) Schleife je Relaxationszeit: Verlustanteil 1/2π, Summe Σ qⁿ·(1/H) = 2π/H = ein Zyklus")
P(f"   Σ qⁿ = {sum(q**n for n in range(20000)):.6f}")
P()
P("3) Stetige Form: Fluss erneuert sich mit der Rate Γ, eigene Zustände mit H")
P("   ν = H/Γ;  ν = 2π  ⟺  Γ = H/2π = a₀/c = θ_c  (kritische Rate der Schließungsregel)")
P()
P("4) Vergleich mit Messungen [G]")
for name,wb,sb,wc,sc in [("Planck 2018 TT,TE,EE+lowE+lensing",0.02237,0.00015,0.1200,0.0012),
                          ("Planck 2018 + BAO",0.02242,0.00014,0.11933,0.00091)]:
    nm=(wb+wc)/wb; snm=np.hypot(sc/wb, wc*sb/wb**2)
    P(f"   {name}: ν = {nm:.4f} ± {snm:.4f}  → Abstand zu 2π: {(nm-nu)/snm:+.2f} σ (ohne Korrelation)")
    P(f"      Vorhersage Ω_c h² = (2π−1)·Ω_b h² = {Q*wb:.5f} ± {Q*sb:.5f}")
wbbn,sbbn=0.02218,0.00055
P(f"   BBN (D/H) Ω_b h² = {wbbn} ± {sbbn}: Ω_c h² = {Q*wbbn:.4f} ± {Q*sbbn:.4f} gegen Planck 0,1200 ± 0,0012 → {(0.1200-Q*wbbn)/np.hypot(Q*sbbn,0.0012):+.2f} σ")
P()
P("5) Gleichwertige Form: a₀ = c·H₀·Ω_b/Ω_m")
for H in [67.4,70.0,73.0]:
    H0=H*1e3/Mpc
    P(f"   H₀ = {H}: c·H₀/ν(Planck 6,364) = {c*H0/6.364:.4e},  c·H₀/2π = {c*H0/nu:.4e} m/s²")
P()
P("6) Kosmischer Kern: Umlauffrequenz des Horizonts (euklidisch) ω = H, Relaxation λ = H → k = λ/ω = 1")
P("   Radiale Antwort ω/(λ²+ω²) bei festem λ maximal für ω = λ; λ_eff = ω + λ²/ω minimal (= 2ω) bei k = 1.")
open("rueckkehr_horizont.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))
