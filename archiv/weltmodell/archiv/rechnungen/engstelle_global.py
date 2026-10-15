# -*- coding: utf-8 -*-
"""
Globale Engstelle: der Horizont des Beobachters als Falte, Δs_H = r_H/6 = c/(6H).
Blick: r_H, Δs_H, S_H gelten je Beobachter und Epoche (Horizont, nicht G); Λ ist I (innere Richtung, Zelle).
1 dΔs_H/dt heute und über die Epoche, gegen Λ und die Mindestrate ε; woraus Δs_H besteht.
2 Kopplung Λ an Δs_H: (a) Λ ∝ Δs_H, (b) dΛ = ν_kosm·dΔs_H, (c) Mittel über alle Horizonte der Zelle
  mit Vakuumteil. Je Λ(a), w_s, w_DE (Rest Ω_x fest), a₀(z), Ω_s früh.
3 Beobachterfehler: Summe über die Horizonte der Zelle.   4 von Hand.
Konstanten nur aus ftm_kern.py.
"""
import math, os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, c, G, l_P, Gpc, yr, Gyr, OMEGA_M, OMEGA_R, OMEGA_B, OMEGA_DE, Q_KOSM, K_S,
                      H0_si, schliesslaenge, omega_vakuum, horizontentropie, kreislauf_mindestrate,
                      theta_, delta_N_eff)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def H_(t):
    P(); P("=" * 110); P(t); P("=" * 110)

H0 = H0_si()
Lam0 = schliesslaenge(A0)
Os0 = omega_vakuum(A0)
Ox = OMEGA_DE - Os0
nu_k = 1 + Q_KOSM
Th0 = theta_(A0, H0)
eps, S_Lam, H_s = kreislauf_mindestrate(Lam0)
SDOT = 7.1e86
def dsH(H): return c / (6 * H)

# Hintergrund des Kerns (Vakuumteil + Rest fest, w = −1): E², q
def E2_std(a): return OMEGA_M * a**-3 + OMEGA_R * a**-4 + OMEGA_DE
def q_std(a):
    E2 = E2_std(a)
    return (0.5 * OMEGA_M * a**-3 + OMEGA_R * a**-4 - OMEGA_DE) / E2

# ------------------------------------------------------------------------------------------------ 1
H_("1  Horizont als Falte: Δs_H = c/(6H), dΔs_H/dt = (c/6)(1 + q), d ln Δs_H/dt = (1 + q)·H")
q0 = q_std(1.0)
P(f"   heute: Δs_H = {dsH(H0)/Gpc:.4f} Gpc,  Λ = {Lam0/Gpc:.3f} Gpc,  Δs_H/Λ = {dsH(H0)/Lam0:.4f} = Θ/(6π) = {Th0/(6*math.pi):.4f}")
P(f"   q0 = {q0:.4f}:  dΔs_H/dt = {(1+q0)/6:.4f}·c = {c*(1+q0)/6:.3e} m/s;  d ln Δs_H/dt = {(1+q0)*H0:.3e} s⁻¹ = {(1+q0)*H0/eps:.2e}·ε")
S_H0 = horizontentropie(H0)
P(f"   S_H = 36π(Δs_H/ℓ_P)² = {36*math.pi*(dsH(H0)/l_P)**2:.4e} k (πc²/(ℓ_P²H²) = {S_H0:.4e});  dS_H/dt = 2(1+q)H·S_H = {2*(1+q0)*H0*S_H0:.3e} k/s = {2*(1+q0)*H0*S_H0/SDOT:.2e}·Ṡ_i")
P("   über die Epoche (Kern-Hintergrund, Λ fest):")
P("        z        1+q      dΔs_H/dt [c]   Δs_H [Gpc]   Δs_H/Λ     d ln Δs_H/dt / ε    Anteil Materie Ω_m(a)")
for z in (1e4, 1090, 10, 3, 1, 0.4, 0, -0.5, -0.9, -0.99):
    a = 1 / (1 + z); E = math.sqrt(E2_std(a)); H = H0 * E; qq = q_std(a)
    P(f"   {z:>8g}   {1+qq:8.4f}   {(1+qq)/6:10.4f}    {dsH(H)/Gpc:10.4e}  {dsH(H)/Lam0:9.3e}   {(1+qq)*H/eps:12.3e}       {OMEGA_M*a**-3/E**2:.4f}")
P(f"   Endzustand Λ fest: q → −1, dΔs_H/dt → 0, Δs_H → c/(6H_s) = {dsH(H_s)/Gpc:.3f} Gpc = {dsH(H_s)/Lam0:.4f}·Λ = (√45/6π²)·Λ = {math.sqrt(45)/(6*math.pi**2):.4f}·Λ")
P(f"   Endzustand Λ wächst mit ε: H = H_s ∝ 1/Λ ⇒ Δs_H ∝ Λ, d ln Δs_H/dt = ε genau (Folge von Λ, nicht Quelle)")
P("   Woraus Δs_H besteht (je Beobachter): Masse in der Hubble-Kugel M(<r_H) = Ω(a)·c³/(2GH) ⇒")
P("   Δs_H = Σ_Komponenten ΣΔs(<r_H) = [Ω_m(a) + Ω_r(a) + Ω_DE(a)]·Δs_H. Die Materiefaltung in der Kugel ändert sich wie")
P("   d ln(ρ_m r_H³)/dt = 3qH: sie fällt, sobald q < 0 (heute −1,59·H): die Kugel wächst nicht über neue Faltung, der Zuwachs von Δs_H ist")
for z in (3, 1, 0.4, 0):
    a = 1/(1+z); E2 = E2_std(a); qq = q_std(a)
    fm = OMEGA_M*a**-3/E2; fr = OMEGA_R*a**-4/E2; fd = OMEGA_DE/E2
    # d/d ln a von f_i·(1/E): Anteil am Zuwachs
    dlnE = -(1+qq)
    gm = (-3 - 3*dlnE) * fm; gr = (-4 - 3*dlnE) * fr; gd = (0 - 3*dlnE) * fd
    tot = gm + gr + gd
    P(f"      z = {z:>3}: d ln Δs_H/d ln a = {1+qq:.4f};  Beiträge Materie {gm:+.4f}, Strahlung {gr:+.5f}, Vakuumteil+Rest {gd:+.4f} (Summe {tot:.4f})")
P("   der Anteil, den Vakuumteil und Rest über das Volumen r_H³ liefern (mitbewegte Materiefaltung fest, Kugel läuft darüber).")
P("   ⇒ Das Wachsen des Horizonts ist keine neue Faltung: Δs_H ist das Budget der eigenen Hubble-Kugel, als Faltung gelesen.")

# ------------------------------------------------------------------------------------------------ 2
H_("2  Kopplung Λ an Δs_H: Λ(a), w_s = −1 + (2/3)·d ln Λ/d ln a, w_DE = −1 + (Ω_s/Ω_DE)(1 + w_s), Rest Ω_x fest")
P(f"   Ω_s0 = {Os0:.4f}, Ω_x = {Ox:.4f}, ν_kosm = {nu_k:.3f}, Θ0 = {Th0:.4f}")

def loese(Lam_fn):
    """E(a) selbstkonsistent: E² = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_x + Ω_s0·(Λ0/Λ(E))²."""
    def E_of(a):
        f = lambda E: E**2 - (OMEGA_M*a**-3 + OMEGA_R*a**-4 + Ox + Os0*(Lam0/Lam_fn(a, E))**2)
        return brentq(f, 1e-3, 1e12)
    return E_of

def tafel(name, Lam_fn):
    E_of = loese(Lam_fn)
    P(f"   {name}")
    P("        z       Λ/Λ0      a₀/a₀(0)   d lnΛ/d ln a    w_s        w_DE       Ω_s(a)/Ω_ges")
    rows = {}
    for z in (1090, 10, 3, 1, 0.5, 0.4, 0.2, 0.0, -0.5):
        a = 1/(1+z); h = 1e-4
        E = E_of(a); L = Lam_fn(a, E)
        Lp = Lam_fn(a*math.exp(h), E_of(a*math.exp(h))); Lm = Lam_fn(a*math.exp(-h), E_of(a*math.exp(-h)))
        dl = (math.log(Lp) - math.log(Lm)) / (2*h)
        ws = -1 + 2/3*dl
        Os = Os0*(Lam0/L)**2
        wde = -1 + Os/(Os+Ox)*(1+ws)
        rows[z] = (L/Lam0, wde, Os/E**2)
        P(f"   {z:>8g}   {L/Lam0:8.4f}   {Lam0/L:8.4f}   {dl:10.4f}   {ws:8.4f}   {wde:8.4f}    {Os/E**2:.4e}")
    # CPL aus w_DE(a) um a = 1
    def wde_a(a):
        h = 1e-4; E = E_of(a); L = Lam_fn(a, E)
        Lp = Lam_fn(a*math.exp(h), E_of(a*math.exp(h))); Lm = Lam_fn(a*math.exp(-h), E_of(a*math.exp(-h)))
        dl = (math.log(Lp) - math.log(Lm)) / (2*h); Os = Os0*(Lam0/L)**2
        return -1 + Os/(Os+Ox)*(2/3*dl)
    w0 = wde_a(1.0); wa = -(wde_a(1.0+1e-3) - wde_a(1.0-1e-3))/2e-3
    rde = [(z, Ox + Os0*(Lam0/Lam_fn(1/(1+z), E_of(1/(1+z))))**2) for z in np.linspace(0, 3, 31)]
    zmax = max(rde, key=lambda t: t[1])[0]
    P(f"      w0 = {w0:.4f}, w_a = {wa:+.4f} (Ableitung bei a = 1; gemessen w_a < 0);  ρ_DE(z) auf 0–3 maximal bei z = {zmax:.1f}" + (" = Rand: monoton, kein Maximum" if zmax >= 3 else ""))
    return rows

# (a) Λ ∝ Δs_H: Λ = Λ0·H0/H = Λ0/E
ra = tafel("(a) Λ = Λ0·Δs_H/Δs_H0 = Λ0/E (Λ atmet mit dem Horizont, Θ fest)", lambda a, E: Lam0/E)
fa = Os0 / 1.0
P(f"      Ω_s/Ω_ges fest = {Os0:.4f} in jeder Epoche ⇒ ΔN_eff = {delta_N_eff(Os0):.2f} (N_eff gemessen 2,99 ± 0,17);  a₀ ∝ H")
# (b) dΛ = ν·dΔs_H
rb = tafel("(b) dΛ = ν_kosm·dΔs_H: Λ = Λ0 + ν_kosm·(c/6H0)(1/E − 1)", lambda a, E: Lam0 + nu_k*dsH(H0)*(1/E - 1))
P(f"      früh (E → ∞): Λ → Λ0(1 − ν_kosm·Θ0/(6π)) = {1 - nu_k*Th0/(6*math.pi):.4f}·Λ0, a₀ früh = {1/(1 - nu_k*Th0/(6*math.pi)):.3f}·a₀(0); Modell: a₀ seit z ≈ 10 fest auf 0,3 %")
# (c) Mittel über alle Horizonte der Zelle, samt Vakuumteil: Λ = Λ_b + ΣΔs_s(Zelle) (ohne bzw. mit ν)
P("   (c) Summe über alle Horizonte der Zelle (Teil 3): Σ Δs_H = ΣΔs aller Anteile der Zelle. Nähme Λ den Vakuumteil der Zelle auf,")
P(f"       Λ = Λ_b + k·a³/Λ²·(Λ0²):  heute ΣΔs_s(Zelle)/Λ = Ω_s/Ω_m = {Os0/OMEGA_M:.3f} (ohne ν), Ω_s/Ω_b = {Os0/OMEGA_B:.2f} (mit ν)")
for lab, k in (("ohne ν", Os0/OMEGA_M), ("mit ν", Os0/OMEGA_B)):
    # Λ/Λ_b = x: x = 1 + kb·a³/x², kb aus heute: x0 = 1 + k (Λ0 = Λ_b(1+k)), Ω_s-Anteil ∝ (Λ0/Λ)²
    x0 = 1 + k; Lb = Lam0/x0
    kb = k*x0**2  # so dass x(1) = x0
    def xa(a): return brentq(lambda x: x - 1 - kb*a**3/x**2, 1, 1e9)
    ws_late = -1 + 2/3*(math.log(xa(1e3*1.001))-math.log(xa(1e3)))/math.log(1.001)
    ws_now = -1 + 2/3*(math.log(xa(1.001))-math.log(xa(0.999)))/(2*0.0010005)
    P(f"       {lab}: Λ_b = {Lb/Gpc:.2f} Gpc (heute Λ = 13,24 Gpc = ν·ΣΔs_b verlangt); w_s heute = {ws_now:.3f}, spät (a = 10³) = {ws_late:.4f} → −1/3 (Λ ∝ a, wie Krümmung 0,22)")
P("   Vorzeichen: Δs_H wächst immer (1 + q > 0 ohne Phantom) ⇒ Λ wächst ⇒ Ω_s ∝ Λ⁻² fällt ⇒ w_s > −1 und ρ_DE fällt monoton.")
P("   Die übrigen zwei Drittel (Ω_x) entstehen so nicht; ein Maximum von ρ_DE bei z ≈ 0,4 verlangte vorher schrumpfendes Λ.")

# ------------------------------------------------------------------------------------------------ 3
H_("3  Beobachterfehler: Horizont je Beobachter, Λ in der inneren Richtung der Zelle")
for a in (0.1, 1.0, 3.0):
    E = math.sqrt(E2_std(a)); H = H0*E
    L_cell = 27.6*Gpc                      # mitbewegte Zelle des Modells (a₀ = 1,1·10⁻¹⁰)
    V = (L_cell*a)**3
    N_H = V/(4*math.pi/3*(c/H)**3)
    sum_H = N_H*dsH(H)
    rho_tot = 3*H**2/(8*math.pi*G)
    fold_cell = G*rho_tot*V/(3*c**2)
    P(f"   a = {a}: N_H = {N_H:.3e} Hubble-Kugeln in der Zelle;  N_H·Δs_H = {sum_H:.4e} m,  G·ρ_ges·V/(3c²) = {fold_cell:.4e} m  (Verhältnis {sum_H/fold_cell:.6f})")
P("   ⇒ Σ_Horizonte Δs_H = ΣΔs_ges(Zelle) = (La)³H²/(8πc²): das Mittel über alle Beobachter ist die Ausdehnungsgleichung selbst,")
P("     als Faltung gelesen. Materieanteil mitbewegt fest (steckt schon in Λ = ν·ΣΔs_b), Strahlung fällt mit 1/a, der Vakuumteil")
P("     wäre Λ als Quelle seiner selbst (Kreis, Teil 2c). Neue Faltung liefert der Horizont auf keiner der beiden Lesarten.")
P(f"   Kopplung (a)/(b) setzt r_H (je Beobachter, Epoche) als Länge in Λ (I): Verwechslung Horizont statt innere Richtung,")
P(f"   dazu Doppelzählung (die Kugel läuft über vorhandene Faltung). Im Endzustand gilt Δs_H = (√45/6π²)Λ als Folge: H_s ∝ 1/Λ.")

# ------------------------------------------------------------------------------------------------ 4
H_("4  Von Hand")
P("   Δs_H = c/(6H0) = 2,998·10⁸/(6·2,184·10⁻¹⁸) = 2,288·10²⁵ m = 0,741 Gpc;  Λ = 13,24 Gpc ⇒ Δs_H/Λ = 0,0560")
P(f"   Kern: {dsH(H0)/Gpc:.4f} Gpc, {dsH(H0)/Lam0:.4f}")
P("   q0 = Ω_m/2 − Ω_DE = 0,157 − 0,687 = −0,530 ⇒ dΔs_H/dt = 0,470/6·c = 0,0783·c = 2,35·10⁷ m/s;  (1+q)H0 = 1,03·10⁻¹⁸ s⁻¹ = 2,7·10¹⁸·ε")
P("   (b) früh: 1 − ν·Δs_H/Λ = 1 − 6,364·0,0560 = 0,644 ⇒ a₀ früh = 1,55·a₀(0)")
P(f"   Kern: q0 = {q0:.4f}, {(1+q0)*H0:.3e} s⁻¹, {(1+q0)*H0/eps:.2e}·ε, 1 − νΔs_H/Λ = {1-nu_k*dsH(H0)/Lam0:.4f}")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "engstelle_global.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
