# -*- coding: utf-8 -*-
"""
Vorschläge (3) und (4) des Autors: Kontaktgeometrie für Atmen und w(z); Festkörper-Analoga.
Blick: Atmen a(t), H und w gehören dem Ganzen (G); ϖ ist die innere Richtung (I); die Regel a = ν g_N wirkt am Ort (L).
Eine Analogie im Labor prüft nur die Form einer Gleichung, nicht die Gravitation.

A Kontaktform auf T³ × (innere Richtung) × ℝ mit Reeb-Feld u; folgt w(a)?
B Anomale Geschwindigkeit ṗ × Ω gegen die Regel (ν − 1) g_N auf einer Kreisbahn.
C Analogon der Rückkehrgleichung mit Erneuerung: RC-Kette mit Rückführung (Bose-Form, Wurzel, versiegeltes Ende).
D Von Hand.
"""
import math, os, sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import A0, G, M_sun, kpc, nu, s_, OMEGA_M, OMEGA_R, omega_vakuum

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def H_(t):
    P(); P("=" * 100); P(t); P("=" * 100)

# ================================================================================================
H_("A1 Kontaktform mit Reeb-Feld ∂_t: explizite Konstruktion (sympy)")
t, x, y, z, w_ = sp.symbols('t x y z varpi', real=True)
k, ell = sp.symbols('k ell', positive=True)
X = [t, x, y, z, w_]
def dform(coef):
    """1-Form (Koeffizientenliste) -> 2-Form als antisymmetrische Matrix F_ij = ∂_i α_j − ∂_j α_i."""
    n = len(coef)
    return sp.Matrix(n, n, lambda i, j: sp.diff(coef[j], X[i]) - sp.diff(coef[i], X[j]))
def kontakt_volumen(coef):
    """Koeffizient von α∧(dα)² vor dt∧dx∧dy∧dz∧dϖ (bis auf Faktor 8): Σ ε α_i Pf-Anteile."""
    F = dform(coef)
    # α∧(dα)² = (1/4) ε^{ijklm} α_i F_jk F_lm · d⁵x /… ; direkt mit Levi-Civita:
    from sympy import LeviCivita
    tot = 0
    for i in range(5):
        for j in range(5):
            for kk in range(5):
                for l in range(5):
                    for m in range(5):
                        e = LeviCivita(i, j, kk, l, m)
                        if e != 0:
                            tot += e * coef[i] * F[j, kk] * F[l, m]
    return sp.simplify(tot / 8)
def reeb(coef):
    F = dform(coef)
    R = sp.symbols('R0:5')
    eqs = [sum(F[i, j] * R[i] for i in range(5)) for j in range(5)] + [sum(coef[i] * R[i] for i in range(5)) - 1]
    sol = sp.solve(eqs, R, dict=True)
    return [sp.simplify(sol[0][r]) for r in R] if sol else None

beta3 = [0, sp.cos(k*z), sp.sin(k*z), 0, 0]
alpha = [1, sp.exp(w_/ell)*sp.cos(k*z), sp.exp(w_/ell)*sp.sin(k*z), 0, 0]
P("α = dt + e^{ϖ/ℓ}(cos kz dx + sin kz dy)   auf T³(x,y,z) × [0,Λ]∋ϖ × ℝ∋t")
P(f"  α∧(dα)² ∝ {kontakt_volumen(alpha)}   (≠ 0 überall ⇒ Kontaktform)")
P(f"  Reeb-Feld R = {reeb(alpha)}   (= ∂_t: das Zeitfeld)")
alpha0 = [1, sp.cos(k*z), sp.sin(k*z), 0, 0]
P(f"  ohne ϖ-Abhängigkeit (e^{{ϖ/ℓ}} → 1): α∧(dα)² = {kontakt_volumen(alpha0)}  ⇒ keine Kontaktform")
f = sp.Function('f')(x, y, z)
alpha_f = [1, sp.cos(k*z), sp.sin(k*z), 0, f]
P(f"  ϖ-unabhängig mit dϖ-Anteil f(x,y,z)dϖ: α∧(dα)² ∝ {sp.simplify(kontakt_volumen(alpha_f))}")
P("    = (rot β₃·∇f)-Anteil; ∫_{T³} dβ₃∧df = ∫ d(β₃∧df) = 0 ⇒ wechselt das Vorzeichen, Nullstelle auf T³.")
P("  Satz (Stokes): Reeb = ∂_t mit α = dt + β (β zeitfrei, weil der Reeb-Fluss α erhält) verlangt dβ symplektisch")
P("  auf dem Schnitt Σ⁴. Σ⁴ geschlossen (T³ × S¹) ⇒ ∫(dβ)² = ∫d(β∧dβ) = 0, Widerspruch. Mit innerer Richtung als")
P("  Kreis S¹ gibt es keine Kontaktform mit Reeb-Feld u. Möglich nur mit Rand (ϖ ∈ [0,Λ]) UND ϖ-abhängigen Koeffizienten.")

# numerische Probe der Stokes-Schranke: zufälliges β auf T⁴, Pf(dβ) = F12F34 − F13F24 + F14F23
rng = np.random.default_rng(1)
n = 12
g1 = np.linspace(0, 2*np.pi, n, endpoint=False)
Th = np.meshgrid(g1, g1, g1, g1, indexing='ij')
def zufall_beta():
    comps = []
    for _ in range(4):
        fsum = np.zeros_like(Th[0]); dsum = [np.zeros_like(Th[0]) for _ in range(4)]
        for _ in range(6):
            kv = rng.integers(-2, 3, 4); a_, b_ = rng.normal(size=2)
            ph = sum(kv[i]*Th[i] for i in range(4))
            fsum += a_*np.cos(ph) + b_*np.sin(ph)
            for i in range(4):
                dsum[i] += kv[i]*(-a_*np.sin(ph) + b_*np.cos(ph))
        comps.append((fsum, dsum))
    return comps
vals = []
for trial in range(5):
    B = zufall_beta()
    F = [[B[j][1][i] - B[i][1][j] for j in range(4)] for i in range(4)]
    Pf = F[0][1]*F[2][3] - F[0][2]*F[1][3] + F[0][3]*F[1][2]
    vals.append((Pf.mean(), Pf.min(), Pf.max()))
P("  Probe, 5 zufällige β auf T⁴ (Gitter 12⁴): Mittel von Pf(dβ) / min / max:")
for v in vals:
    P(f"    {v[0]: .2e}  {v[1]: .2f}  {v[2]: .2f}")

# ------------------------------------------------------------------------------------------------
H_("A2 Atmen im Kontaktbild: Reeb-Fluss oder nicht?")
a = sp.Function('a')(t)
alpha_a = [1, a*sp.exp(w_/ell)*sp.cos(k*z), a*sp.exp(w_/ell)*sp.sin(k*z), 0, 0]
R_a = reeb(alpha_a)
P("α = dt + a(t)·e^{ϖ/ℓ}(cos kz dx + sin kz dy)  (Atmen in α eingesetzt)")
P(f"  Reeb-Feld: R_t = {sp.simplify(R_a[0])}")
P(f"             R_ϖ = {sp.simplify(R_a[4])}")
P("  ⇒ mit ȧ ≠ 0 ist ∂_t kein Reeb-Feld mehr (R erhält α: L_Rα = 0; Atmen ändert α). Atmen ist höchstens ein")
P("    konformer Kontaktfluss L_Xα = fα; die Dehnung e^{ϖ/ℓ} läuft längs ϖ (Liouville-Richtung), nicht längs u.")
P("  Metrische Lesart α = −u♭/c (g_tt = −c², g_ti = −c·β_i): i_u dα = Beschleunigung von u ⇒ Reeb ⇔ u geodätisch;")
P("  dα|_{u⊥} = 2ω (Wirbel von u); Kontakt ⇔ Pf(ω) ≠ 0 im 4-dim. u⊥ = T³ × ϖ: u muss in zwei Ebenen drehen,")
P("  eine davon mischt Raum und innere Richtung (Verwindung). Ein wirbelfreies, gleich atmendes u (u♭ = −c dt)")
P("  gibt dα = 0: Frobenius-integrierbar, das Gegenteil von Kontakt.")
P("  Der Reeb-Fluss trägt dα mit (L_R dα = 0): Koordinatenkomponenten fest ⇒ orthonormal ω_xy ∝ a⁻², ω_zϖ ∝ a⁻¹.")
P("  Kontakt-Torsion τ = L_u g = 2H·h (Raum): das Atmen ist genau die Torsion; H(t) selbst legt α nicht fest.")

# ------------------------------------------------------------------------------------------------
H_("A3 Folgt w(a)? Beiträge, die der Reeb-Fluss erzeugt: Potenzen a⁻ⁿ mit festem Exponenten")
P("Lie-getragene Wirbel: ω² ∝ a⁻⁴ (Raum-Raum) und a⁻² (Raum-ϖ) ⇒ Glieder mit w = n/3 − 1 = +1/3 bzw. −1/3.")
P("Dunkle Energie = Konstante + ε·a⁻ⁿ (Amplitude ε frei, Integrationskonstante der Wirbel):")
P("  n    ε/ρ_Λ     w0 = w(1)   wa = −dw/da|₁   Quadrant   (ε so, dass ρ_DE > 0 für a ≥ 0,5)")
def w_von_a(aa, terms):
    num = sum((n_/3)*e*aa**-n_ for n_, e in terms); den = 1 + sum(e*aa**-n_ for n_, e in terms)
    return -1 + num/den
def w0wa(terms, h=1e-6):
    return w_von_a(1.0, terms), -(w_von_a(1 + h, terms) - w_von_a(1 - h, terms))/(2*h)
def quadrant(w0, wa):
    return ("w0>−1, wa<0 (gemessen)" if wa < 0 else "w0>−1, wa>0") if w0 > -1 else ("w0<−1, wa<0" if wa < 0 else "w0<−1, wa>0")
for n_ in (2, 4):
    for eps in (0.05, 0.2, -0.02, -0.05):
        w0, wa = w0wa([(n_, eps)])
        P(f"  {n_}   {eps:+.2f}     {w0:+.4f}      {wa:+.4f}       {quadrant(w0, wa)}")
P("  Allgemein: w = −1 + (n/3)·f(a), f = εa⁻ⁿ/(1+εa⁻ⁿ); sign(1+w0) = sign(ε) = sign(w_a). Der gemessene Quadrant")
P("  (w0 > −1, w_a < 0) verlangt entgegengesetzte Vorzeichen: mit Konstante + einer Potenz unerreichbar.")
aa = np.linspace(0.5, 1.0, 50)
hits = []
for e2 in np.linspace(-0.3, 0.3, 121):
    for e4 in np.linspace(-0.3, 0.3, 121):
        if np.any(1 + e2*aa**-2 + e4*aa**-4 <= 0.05):
            continue
        w0, wa = w0wa([(2, e2), (4, e4)])
        if w0 > -1 and wa < 0:
            hits.append((e2, e4, w0, wa))
if hits:
    e2s = [h[0] for h in hits]; e4s = [h[1] for h in hits]
    P(f"  Beide Potenzen (a⁻², a⁻⁴), Vorzeichen frei: Quadrant erreicht für ε₂ ∈ [{min(e2s):+.3f}, {max(e2s):+.3f}], "
      f"ε₄ ∈ [{min(e4s):+.3f}, {max(e4s):+.3f}]")
    P(f"    nur mit ε₂ > 0 und ε₄ < 0 ({all(h[0] > 0 and h[1] < 0 for h in hits)}): zwei freie Amplituden mit Gegenvorzeichen,")
    P("    ε₄ < 0 heißt negativer Wirbelbeitrag; keine Vorhersage von w0, w_a.")
else:
    P("  Beide Potenzen: Quadrant im Raster nicht erreicht.")
P(f"  Zum Vergleich im Modell: Vakuumteil Ω_s = {omega_vakuum():.4f} bei festem Λ, w_s = −1.")

# ================================================================================================
H_("B Anomale Geschwindigkeit ṗ×Ω gegen die Regel (ν−1)g_N, Kreisbahn")
Mb = 5e10 * M_sun
P("Punktmasse 5·10¹⁰ M☉ (je Masseneinheit): Regel a = ν g_N (parallel zu g_N); Berry: ṙ = p − ṗ×Ω, ṗ = g_N(r).")
P("  2D, Ω = Ω_z ẑ fest (zeitumkehr-brechend): v_a = −g_N × Ω liegt in der Ebene, senkrecht zu g_N (tangential).")
P("  Stationäre Kreisbahn: |p| = g/φ̇, rφ̇ = g/φ̇ ± Ωg  ⇒  rφ̇² ∓ Ωgφ̇ − g = 0 (± = Umlaufsinn)")
def v_berry(g, r, Om, sinn):
    ph = (sinn*Om*g + math.sqrt(Om**2*g**2 + 4*r*g)) / (2*r)
    return r*ph
r_ref = 20*kpc
g_ref = G*Mb/r_ref**2
v_rule = math.sqrt(float(nu(g_ref))*g_ref*r_ref)
# Ω so, dass rechtläufig bei 20 kpc die Regel getroffen wird
from scipy.optimize import brentq
Om = brentq(lambda O: v_berry(g_ref, r_ref, O, +1) - v_rule, 0, 1e18)
P(f"  Ω so gewählt, dass die rechtläufige Bahn bei 20 kpc die Regel trifft: Ω = {Om:.3e} s")
P("   r[kpc]  a_N/a₀   v_N     v_Regel  v_Berry(+)  v_Berry(−)  [km/s]   g_N·v_a")
for rk in (2, 5, 10, 20, 40, 80):
    r = rk*kpc; g = G*Mb/r**2
    vN = math.sqrt(g*r); vR = math.sqrt(float(nu(g))*g*r)
    vp = v_berry(g, r, Om, +1); vm = v_berry(g, r, Om, -1)
    P(f"   {rk:5d}  {g/A0:7.3f}  {vN/1e3:6.1f}  {vR/1e3:7.1f}   {vp/1e3:8.1f}   {vm/1e3:8.1f}            0 (exakt, v_a ⊥ g_N)")
# unabhängige Probe: Integration der Bewegungsgleichung (rechtläufig, r = 20 kpc)
def rhs(tt, Y, sinn):
    rx, ry, px, py = Y
    r = math.hypot(rx, ry); gx, gy = -G*Mb*rx/r**3, -G*Mb*ry/r**3
    # ṙ = p − ṗ×Ω, ṗ×Ω (Ω ẑ) = (ṗ_y Ω, −ṗ_x Ω)
    return [px - gy*Om, py + gx*Om, gx, gy]      # Ω = +ẑ fest
for sinn in (+1, -1):
    vc = v_berry(g_ref, r_ref, Om, sinn); ph = vc/r_ref
    p0 = g_ref/ph
    Y0 = [r_ref, 0, 0, -sinn*p0]   # Sinn +1: Umlauf, in dem v_a zu p addiert
    T = 4*2*math.pi/ph
    sol = solve_ivp(rhs, (0, T), Y0, args=(sinn,), method='DOP853', rtol=1e-10, atol=1e-6*r_ref, dense_output=True)
    rr = np.hypot(sol.y[0], sol.y[1])
    ang = np.unwrap(np.arctan2(sol.y[1], sol.y[0]))
    P(f"  Integration, Sinn {sinn:+d}: r bleibt {rr.min()/kpc:.4f}–{rr.max()/kpc:.4f} kpc, v = r·φ̇ = "
      f"{abs(ang[-1]-ang[0])/T*r_ref/1e3:.2f} km/s (Formel {vc/1e3:.2f})")
P("  (Sinn −1: die Kreislösung existiert, ist in der Integration aber nicht stabil, die Bahn läuft nach außen.)")
P("  3D-Weyl-Knoten Ω = ±p̂/(2p²)·ħ-Einheiten: auf der Kreisbahn ist p̂ tangential, g_N radial ⇒ ṗ×Ω ∥ Bahnachse:")
P("    eine Drift aus der Ebene, kein radialer Anteil, v_c unverändert.")
P("  Ergebnis: v_a·g_N = 0 immer; die Regel ist (ν−1)g_N ∥ g_N. Berry hängt vom Umlaufsinn ab (Regel nicht, 07),")
P("    fällt ∝ g_N nach außen (Regel: v → (GMa₀)^{1/4} flach). Nicht identisch, auch nicht der Form nach.")

# ================================================================================================
H_("C Analogon der Rückkehr: RC-Kette (Kabelgleichung) mit Rückführung V(0) = V_ein + V(n_R)")
P("Kontinuum: ∂_t V + (g_m/c_m)V = (1/(r_a c_m))∂_x²V  ≡  ∂_t u + Γu = D∂_ϖ²u  (x ↔ ϖ, g_m ↔ Γ, n_R ↔ ϖ_R).")
P("Rückführung mit Verstärkung 1 (Summierverstärker) = Erneuerung u(0) = g_N + u(ϖ_R) ⇒ V(0)/V_ein = 1/(1 − e^{−s}).")
Rs, Cs, nR, N = 1e3, 100e-9, 20, 400
P(f"Bauteile: R = {Rs:.0f} Ω je Glied, C = {Cs*1e9:.0f} nF, Querleitwert G je Glied (einstellbar), Abgriff n_R = {nR},")
P(f"  {N} Glieder, Ende mit R_∞ = R/2 + √(R²/4 + R/G) (Halbgerade, beschränkt) bzw. offen (versiegelt, Schluss bei Λ).")
def kette(Gq, ende="inf", nL=None, Rtol=None, Gtol=None):
    """Knotenanalyse (DC). Knoten 0..N; Rückführung V0 = 1 + V[nR]."""
    NN = N if nL is None else nL
    Rv = np.full(NN, Rs) if Rtol is None else Rtol[:NN]
    Gv = np.full(NN + 1, Gq) if Gtol is None else Gtol[:NN + 1]
    A = np.zeros((NN + 1, NN + 1)); b = np.zeros(NN + 1)
    A[0, 0] = 1; A[0, nR] = -1; b[0] = 1.0
    for i in range(1, NN + 1):
        A[i, i] += Gv[i]
        A[i, i - 1] -= 1/Rv[i - 1]; A[i, i] += 1/Rv[i - 1]
        if i < NN:
            A[i, i + 1] -= 1/Rv[i]; A[i, i] += 1/Rv[i]
    if ende == "inf":
        Rinf = Rs/2 + math.sqrt(Rs**2/4 + Rs/Gq)
        A[NN, NN] += 1/Rinf  # Rest der Kette ab Knoten NN (Serien-R, dann Leiter) = R_∞
    V = np.linalg.solve(A, b)
    return V
P("   G [S]      RG        s=n_R√(RG)  s_eff       ν Kette    1/(1−e^{−s_eff})  1/(1−e^{−s})   1−V(n_R)/V(0)")
tab = []
for Gq in np.geomspace(1e-8, 1e-3, 11):
    V = kette(Gq)
    RG = Rs*Gq
    mu = 1 + RG/2 - math.sqrt(RG + RG**2/4)
    s_eff = -nR*math.log(mu); s = nR*math.sqrt(RG)
    tab.append((s, V[0]))
    P(f"   {Gq:.2e}  {RG:.2e}  {s:9.4f}   {s_eff:9.4f}  {V[0]:10.5f}  {1/(-math.expm1(-s_eff)):10.5f}      "
      f"{1/(-math.expm1(-s)):10.5f}    {1 - V[nR]/V[0]:.5f}")
ss = np.array([x[0] for x in tab]); nn = np.array([x[1] for x in tab])
sl = np.polyfit(np.log(ss[:3]**2), np.log(nn[:3] - 0.5), 1)[0]
P(f"  Wurzel: d ln(ν − ½)/d ln G im tiefen Bereich (s < 0,1) = {sl:.4f} (erwartet −½)")
P("  Die Kette gibt die Bose-Form exakt in s_eff = n_R·arcosh(1 + RG/2); s_eff = s(1 − RG/24 + …).")
# Versiegeltes Ende bei n_Λ (offener Abschluss): Schluss bei Λ
P("  Versiegeltes Ende (offen) bei n_Λ statt R_∞: Übertrag cosh((n_Λ−n_R)/λ)/cosh(n_Λ/λ), λ = 1/√(RG) Glieder")
for nL in (40, 200, 1000, 3000):
    Gq = (math.pi/(2*nL))**2/Rs     # Γ = kleinster Eigenwert Dπ²/(4Λ²): n_Λ/λ = π/2
    V1 = kette(Gq, ende="offen", nL=nL); V2 = kette(Gq)
    P(f"   n_Λ = {nL:4d} (n_Λ/n_R = {nL/nR:5.0f}), n_Λ/λ = π/2: ν_versiegelt/ν_Halbgerade − 1 = {V1[0]/V2[0]-1:.4f}")
P(f"   Grenzwert n_Λ ≫ n_R (07, Schluss bei Λ): 1/tanh(π/2) − 1 = {1/math.tanh(math.pi/2)-1:.4f}")
# Toleranzen
P("  Bauteiltoleranz 1 % (R und G, 200 Ziehungen), Streuung von ν:")
for Gq in (1e-7, 1e-5, 1e-4):
    vs = []
    for _ in range(200):
        Rt = Rs*(1 + 0.01*rng.standard_normal(N)); Gt = Gq*(1 + 0.01*rng.standard_normal(N + 1))
        vs.append(kette(Gq, Rtol=Rt, Gtol=Gt)[0])
    V0 = kette(Gq)[0]
    P(f"   G = {Gq:.0e}: ν = {V0:.4f} ± {np.std(vs):.4f} ({np.std(vs)/V0*100:.2f} %)")
tauD = nR**2*Rs*Cs
P(f"  Zeitskala τ_D = n_R²RC = {tauD*1e3:.1f} ms; Kern K̂(p) = e^{{−√((p+Γ)τ_D)}}: mit Wechselspannung ω ≈ 1/τ_D = "
  f"{1/(2*math.pi*tauD):.1f} Hz über Lock-in messbar (Warburg-Phase −45° bei Γ → 0).")
P("  Kennzeichnung: Analog-Prüfung der Form (Bose-Form, Wurzel, versiegeltes Ende), KEINE Prüfung der Gravitation.")
P("  Weyl-Halbmetall, ³He-B, Magnon-Kondensat: keines hat eine Halbgerade mit Erneuerung u(0) = g + u(ϖ_R);")
P("    ihre Bose-Formen sind Besetzungen von Quasiteilchen (Energie, Temperatur), die Regel zählt Rückläufe.")

# ================================================================================================
H_("D Von Hand")
P("  Kette: R = 1 kΩ, G = 10 µS ⇒ RG = 0,01; n_R = 20: s = 20·0,1 = 2; ν = 1/(1 − e⁻²) = 1/(1 − 0,13534) = 1,15652.")
V = kette(1e-5)
P(f"     Kette (Knotenanalyse): ν = {V[0]:.5f}; mit s_eff = 20·arcosh(1,005) = 20·0,099958 = 1,99917 → "
  f"{1/(1-math.exp(-1.99917)):.5f}")
P("  Berry: v_a = −ṗ×Ω ⇒ ṗ·v_a = −ṗ·(ṗ×Ω) = 0 (Spatprodukt mit zwei gleichen Vektoren).")
P("  Kontakt: T⁴ geschlossen, β global ⇒ ∫(dβ)² = ∫d(β∧dβ) = 0; ein symplektisches (dβ)² hat festes Vorzeichen ≠ 0.")
P("  Quadrant: ε = +0,1, n = 2: w0 = −1 + (2/3)·0,1/1,1 = −0,9394; bei a = 0,5: f = 0,4/1,4, w = −0,8095 > w0 ⇒ w_a > 0.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kontakt_atmen.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
