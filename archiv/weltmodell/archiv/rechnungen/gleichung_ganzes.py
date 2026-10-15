# -*- coding: utf-8 -*-
"""
Das Medium als eine Gleichung auf T³ × [0, Λ] (innere Richtung ϖ, geschlossen bei Λ).

Bulk:      (∇² + ∂_ϖ²) u = 0            für 0 < ϖ < Λ,   ∂_ϖu = 0 bei ϖ = Λ
Ort ϖ = 0: u(x, 0) = Φ_N(x) + u(x, ϖ_R)  (Erneuerung: Wert am Ort = Faltung + zurückgekehrter Wert aus der Tiefe ϖ_R)
Projektion auf ϖ = 0 je Mode k:  u(k,0) = Φ_N(k) · m(k),   m = 1/(1 − e^{−kϖ_R})  (offen),   1/(1 − cosh k(Λ−ϖ_R)/cosh kΛ)  (bei Λ geschlossen)
Zwei Lesarten des Operators A in ∂_ϖ²u = A·u:
  (O) A = −∇²   linear, nichtlokal (harmonische Fortsetzung, Caffarelli-Silvestre, DtN = |k|); ϖ_R aus der Gesamtmasse
  (R) A = Γ/D = a_N/(GM) lokal (Verlust Γ = a_N/c, D = GM/c), nichtlinear;  ϖ_R² = D·τ_D ⇒ s² = Γτ_D = a_N/a₀  (örtliche Regel)
Punktmasse, Kreisbahnmode k = 1/r: −∇² → 1/r² = Γ/D, beide gleich. Geprüft wird, wo sie sich trennen.
Teile: A Mode in ϖ mit finiten Differenzen · B Punktmasse (O) exakt · C Exponentialscheibe · D SPARC · E Ganzes (Nullmode, Kac) · F Handprobe
"""
import math
import numpy as np
from scipy import integrate, sparse
from scipy.sparse.linalg import spsolve
from scipy.special import j1, i0e, k0e, i1e, k1e
from ftm_kern import G, c, A0, kpc, Gpc, M_sun, Gyr, hbar, kB, nu, zelle_aus_a0
import sparc_flaechendichte as sf
import sparc_universal as su

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)


# ------------------------------------------------------------------ A  Mode in ϖ, finite Differenzen
def mode_fd(kappa_L, W_L, geschlossen, N=20000):
    """u'' = κ²u auf [0, W], Einheit L = ϖ_R. u(0) − u(L) = 1; bei W: u = 0 (offen, W groß) oder u' = 0 (geschlossen, W = Λ)."""
    h = W_L / N; iL = int(round(1.0 / h))
    k2 = kappa_L**2
    main = np.full(N + 1, -2.0 / h**2 - k2); off = np.full(N, 1.0 / h**2)
    Amat = sparse.diags([off, main, off], [-1, 0, 1], format='lil')
    b = np.zeros(N + 1)
    Amat[0, :] = 0; Amat[0, 0] = 1.0; Amat[0, iL] = -1.0; b[0] = 1.0
    if geschlossen:
        Amat[N, :] = 0; Amat[N, N] = 2.0 / h**2 * -1 - k2; Amat[N, N - 1] = 2.0 / h**2
    else:
        Amat[N, :] = 0; Amat[N, N] = 1.0
    u = spsolve(Amat.tocsr(), b)
    return u[0]

P("A  Mode in ϖ (finite Differenzen, Einheit ϖ_R = 1): u(0) gegen m(k)")
P("   κϖ_R |  FD offen  | 1/(1−e^{−κϖ_R}) |  FD Λ=20ϖ_R | cosh-Form")
for kl in (0.05, 0.3, 1.0, 3.0):
    fo = mode_fd(kl, 60.0 / min(kl, 1) if kl < 1 else 60.0, False, N=60000)
    fg = mode_fd(kl, 20.0, True, N=40000)
    ex = 1 / (-math.expm1(-kl)); eg = 1 / (1 - math.cosh(kl * 19) / math.cosh(kl * 20))
    P(f"   {kl:5.2f} | {fo:10.5f} | {ex:10.5f}      | {fg:10.5f}  | {eg:10.5f}")


# ------------------------------------------------------------------ B  Punktmasse, Lesart (O)
def extra_arctan(rho, Nmax=None):
    """Zusatz in GM/L²: Σ_{n≥1} (2/π)[arctan(ρ/n)/ρ² − n/(ρ(n²+ρ²))]  (Poisson-Glättung der Punktmasse in der Höhe nϖ_R)."""
    Nmax = Nmax or int(max(2e5, 2000 * rho))
    n = np.arange(1, Nmax + 1, dtype=float)
    s = np.sum((2 / math.pi) * (np.arctan(rho / n) / rho**2 - n / (rho * (n**2 + rho**2))))
    return s + (2 / math.pi) * (2 * rho / 3) / (2 * Nmax**2)   # Schwanz 2ρ/(3n³)


def extra_quad(rho):
    """Gleiche Größe als Fourier-Integral: (2/π)∫ x j₁(xρ) e^{−x}/(1−e^{−x}) dx, j₁ sphärisch."""
    def f(x):
        z = x * rho
        j = (math.sin(z) / z**2 - math.cos(z) / z) if z > 1e-4 else z / 3
        return x * j / math.expm1(x)
    v, _ = integrate.quad(f, 0, 80, limit=4000)
    return (2 / math.pi) * v


P("\nB  Punktmasse: Lesart (O) (harmonische Fortsetzung, ϖ_R = L) gegen örtliche Regel (R)")
P("   Einheit GM/L². Zusatz Reihe = Quadratur (Probe). Tiefgrenze (O): a → (2/π)·GM/(L·r)")
P("   r/L  | Zusatz Reihe | Zusatz Quadr. | (a−a_N)·r·L/GM | (2/π)")
for rho in (0.3, 1, 3, 10, 30, 100):
    e1 = extra_arctan(rho); e2 = extra_quad(rho)
    P(f"   {rho:5g} | {e1:12.6g} | {e2:12.6g}  | {e1 * rho:10.5f}     | {2 / math.pi:.5f}")
P("   ⇒ Tully-Fisher-Nullpunkt der Lesart (O) mit ϖ_R = r*: v⁴ = (2/π)²GMa₀ = 0,405·GMa₀. Gleicher Nullpunkt verlangt ϖ_R = (2/π)·r*.")

P("\n   Form bei gleichem Nullpunkt (ϖ_R = (2/π)r*): a_O/a_R, r in r*")
L_r = 2 / math.pi
for x in (0.1, 0.3, 0.5, 1, 2, 3, 10, 30):
    rho = x / L_r
    aN = 1 / x**2                                  # in GM/r*² = a₀
    aO = aN + extra_arctan(rho) / L_r**2
    aR = float(nu(aN * A0) * aN)
    P(f"   r/r* = {x:5g}: a_N/a₀ = {aN:8.4g}  a_O/a₀ = {aO:8.4f}  a_R/a₀ = {aR:8.4f}  Verhältnis {aO / aR:6.4f}")


# ------------------------------------------------------------------ C  Scheibe: Hankel-Form in der Ebene
def I_int(x):
    """I(x) = ∫₀^∞ n̄(√(x²+y²))/(x²+y²) dy, n̄(z) = 1/(e^z − 1); zerlegt in regulären Teil + 1/x² − π/(4x)."""
    def f(y):
        z = math.hypot(x, y)
        if z < 1e-3:
            return (1 / 12 - z**2 / 720) / z
        return (1 / math.expm1(z) - 1 / z + 0.5) / z**2
    v1, _ = integrate.quad(f, 0, 50, limit=500)
    v2, _ = integrate.quad(lambda y: (-1 / math.hypot(x, y) + 0.5) / (x * x + y * y), 50, np.inf, limit=200)
    return v1 + v2 + 1 / x**2 - math.pi / (4 * x)

XG = np.logspace(-4, math.log10(45), 500)
EG = np.array([(4 / (2 * math.pi)) * x * I_int(x) for x in XG])     # E(x) = m_eff − 1 (in der Ebene), x = kϖ_R
def E_eff(x):
    x = np.asarray(x, float); r = np.empty_like(x)
    lo = x < XG[0]; hi = x > XG[-1]; mid = ~(lo | hi)
    r[lo] = 2 / (math.pi * x[lo]); r[hi] = 0.0
    r[mid] = np.exp(np.interp(np.log(x[mid]), np.log(XG), np.log(np.maximum(EG, 1e-300))))
    return r

P("\nC0 Ebenen-Multiplikator E(x) = m_eff − 1, x = kϖ_R: Probe x·E → 2/π für x → 0")
for x in (1e-3, 1e-2, 0.1, 1, 3, 10):
    P(f"   x = {x:6g}: x·E = {x * E_eff(np.array([x]))[0]:.5f}")

XQ = np.linspace(0, 45, 450001)[1:]
EQ = E_eff(XQ); RQ = EQ - (2 / math.pi) * np.exp(-XQ) / XQ
def punkt_extra_ebene(rho):
    """∫E(x)J₁(xρ)x dx: analytischer Teil (2/π)(1 − 1/√(1+ρ²))/ρ + Rest numerisch."""
    ana = (2 / math.pi) * (1 - 1 / math.sqrt(1 + rho**2)) / rho
    return ana + np.trapezoid(RQ * j1(XQ * rho) * XQ, XQ)

P("   Probe Hankel in der Ebene gegen 3D-Reihe (Punktmasse): r/L, Hankel, Reihe")
for rho in (0.3, 1, 3, 10):
    P(f"   {rho:5g}: {punkt_extra_ebene(rho):.6f}  {extra_arctan(rho):.6f}")


def a_op_ebene(R, gN, Hk_fn, GM, L, kmax_L=45.0, nk=None):
    """a_O(R) = g_N + GM/L²·Punkt(R/L) + ∫(H(k) − GM)E(kL)J₁(kR)k dk."""
    Rmax = R.max()
    nk = nk or int(min(4e5, max(20000, 40 * kmax_L / L * Rmax)))
    k = np.linspace(0, kmax_L / L, nk + 1)[1:]
    Hk = Hk_fn(k)
    w = (Hk - GM) * E_eff(k * L) * k
    res = np.array([np.trapezoid(w * j1(k * Ri), k) for Ri in R])
    pm = np.array([punkt_extra_ebene(Ri / L) for Ri in R]) * GM / L**2
    return gN + pm + res


def expdisk(M, Rd, R):
    y = R / (2 * Rd); Sig0 = M / (2 * math.pi * Rd**2)
    V2 = 4 * math.pi * G * Sig0 * Rd * y**2 * (i0e(y) * k0e(y) - i1e(y) * k1e(y))
    return V2 / R

P("\nC  Exponentialscheibe (dünn), gleicher Nullpunkt ϖ_R = (2/π)r*: a_O/a_R in der Ebene (a_R = ν(g_R)g_R)")
Msc = 5e10 * M_sun; rs = math.sqrt(G * Msc / A0)
for fac in (0.3, 1.0, 3.0):
    Rd = fac * rs; L = (2 / math.pi) * rs
    R = np.array([0.25, 0.5, 1, 2, 3.5, 5, 8, 15]) * Rd
    gN = expdisk(Msc, Rd, R)
    aO = a_op_ebene(R, gN, lambda k: G * Msc / (1 + (k * Rd)**2)**1.5, G * Msc, L)
    aR = nu(gN) * gN
    P(f"   R_d = {fac} r*: " + "  ".join(f"{x:g}R_d:{o / r_:.3f}" for x, o, r_ in zip(R / Rd, aO, aR)))


# ------------------------------------------------------------------ D  SPARC
UPS_D, UPS_B = 0.5, 0.7
NAMEN = ['NGC3198', 'NGC2403', 'NGC6503', 'NGC2841', 'NGC7331', 'UGC00128', 'F568-3', 'NGC3109', 'DDO154', 'DDO168']

def galaxie(nm):
    import os
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat'), comments='#')
    R = a[:, 0] * kpc; Vo = a[:, 1] * 1e3; eV = a[:, 2] * 1e3
    Vg, Vd, Vb = a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3
    V2 = Vg * np.abs(Vg) + UPS_D * Vd**2 + UPS_B * Vb**2
    return R, Vo, eV, V2

def H_daten(R, V2):
    """H(k) = ∫g_N J₁(kR)R dR = ∫V²J₁(kR)dR bis R_l (linear interpoliert, V²(0) = 0) + GM_t·(1 − ∫₀^{kR_l}J₁(x)/x dx)."""
    Rl = R[-1]; GMt = V2[-1] * Rl
    Rg = np.linspace(0, Rl, 4001); V2g = np.interp(Rg, np.concatenate([[0], R]), np.concatenate([[0], V2]))
    xs = np.linspace(0, 3000, 3000001); f = np.where(xs > 0, j1(xs) / np.where(xs > 0, xs, 1), 0.5)
    T = integrate.cumulative_trapezoid(f, xs, initial=0)
    def Hk(k):
        out_ = np.empty_like(k)
        for i0 in range(0, len(k), 2000):
            kk = k[i0:i0 + 2000]
            out_[i0:i0 + 2000] = np.trapezoid(V2g[None, :] * j1(kk[:, None] * Rg[None, :]), Rg, axis=1)
        return out_ + GMt * (1 - np.interp(k * Rl, xs, T))
    return Hk, GMt

P("\nD  SPARC (Υ_d = 0,5, Υ_b = 0,7, a₀ = 1,1·10⁻¹⁰): χ²/N örtliche Regel (R) gegen Lesart (O)")
P("   (O1: ϖ_R = (2/π)·√(GM_t/a₀), O2: ϖ_R = √(GM_t/a₀); GM_t = V_bar²R am letzten Punkt)")
P("   Name       N  |  (R)    | (O1)    | (O2)    | a_O1/a_R innen, Mitte, außen")
summe = np.zeros(3); nsum = 0
for nm in NAMEN:
    R, Vo, eV, V2 = galaxie(nm)
    gN = V2 / R
    Hk, GMt = H_daten(R, V2)
    rs_ = math.sqrt(GMt / A0)
    aR = nu(np.abs(gN)) * gN
    aO1 = a_op_ebene(R, gN, Hk, GMt, (2 / math.pi) * rs_)
    aO2 = a_op_ebene(R, gN, Hk, GMt, rs_)
    chis = []
    for a in (aR, aO1, aO2):
        Vm = np.sqrt(np.maximum(a * R, 0))
        chis.append(np.sum(((Vm - Vo) / eV)**2))
    chis = np.array(chis); summe += chis; nsum += len(R)
    rr = aO1 / aR; i = [0, len(R) // 2, -1]
    P(f"   {nm:9s} {len(R):3d} | {chis[0] / len(R):7.2f} | {chis[1] / len(R):7.2f} | {chis[2] / len(R):7.2f} | "
      + ", ".join(f"{rr[j]:.3f}" for j in i))
P(f"   Summe      {nsum:3d} | {summe[0] / nsum:7.2f} | {summe[1] / nsum:7.2f} | {summe[2] / nsum:7.2f}")

P("\n   Υ_d frei (0,2 … 1,0, Schritt 0,1), je Galaxie bestes χ²/N und Υ_d: (R) gegen (O1)")
UPS_GRID = np.round(np.arange(0.2, 1.01, 0.1), 2)
sR = sO = 0.0
for nm in NAMEN:
    best = {}
    for u in UPS_GRID:
        UPS_D = u
        R, Vo, eV, V2 = galaxie(nm); gN = V2 / R
        Hk, GMt = H_daten(R, V2); rs_ = math.sqrt(GMt / A0)
        for key, a in (("R", nu(np.abs(gN)) * gN), ("O1", a_op_ebene(R, gN, Hk, GMt, (2 / math.pi) * rs_))):
            x2 = np.sum(((np.sqrt(np.maximum(a * R, 0)) - Vo) / eV)**2) / len(R)
            if key not in best or x2 < best[key][0]:
                best[key] = (x2, u)
    sR += best["R"][0] * len(R); sO += best["O1"][0] * len(R)
    P(f"   {nm:9s} | (R) {best['R'][0]:7.2f} bei Υ_d = {best['R'][1]:.1f} | (O1) {best['O1'][0]:7.2f} bei Υ_d = {best['O1'][1]:.1f}")
P(f"   Summe     | (R) {sR / nsum:7.2f}               | (O1) {sO / nsum:7.2f}")
UPS_D = 0.5


# ------------------------------------------------------------------ E  Das Ganze (Nullmode) und Kac
Lam = c**2 / (2 * A0); tauD = 2 * Lam / c
P("\nE  Projektion auf G")
P(f"   Λ = c²/(2a₀) = {Lam / Gpc:.2f} Gpc (Länge der inneren Richtung),  τ_D = 2Λ/c = {tauD / Gyr:.1f} Mrd. Jahre,  a₀ = c/τ_D = {c / tauD:.3e}")
P(f"   Periode 2Λ/c als Temperatur (keine Wärme): T_I = ħc/(2Λ k_B) = {hbar * c / (2 * Lam * kB):.3e} K")
Lz = zelle_aus_a0(A0); kmin = 2 * math.pi / Lz
P(f"   Raumzelle L = {Lz / Gpc:.1f} Gpc: kleinste Mode k = 2π/L, kΛ = {kmin * Lam:.2f}, DtN/k = tanh(kΛ) = {math.tanh(kmin * Lam):.4f}")
P("   Nullmode k = 0: u'' = 0, u(ϖ_R) = u(0) ⇒ m(0) = 1/(1 − 1): Erneuerung nur lösbar, wenn der Mittelwert keine Quelle ist (ρ − ρ̄).")
# Kac auf der Kette 0..N mit reflektierendem Deckel: mittlere Rückkehrzeit zum Ort
for N in (10, 100, 1000):
    h = Lam / N; deg = np.full(N + 1, 2.0); deg[0] = 1; deg[-1] = 1
    T = deg.sum() / deg[0] * h / c                     # Kac: ⟨T_0⟩ = Σdeg/deg(0) Schritte à h/c
    P(f"   Kac, Kette mit {N:4d} Ebenen: ⟨T_Ort⟩·c/(2Λ) = {T * c / (2 * Lam):.6f}")
# Abweichung der Lesart (O) durch den Abschluss bei Λ: Punktmasse Milchstraße
P("   Abschluss bei Λ in (O), Punktmasse 6·10¹⁰ M☉ (ϖ_R = (2/π)r*): Tiefgrenze bleibt bis r ~ Λ; Faktor m_Λ/m_offen bei k = 1/r:")
Mmw = 6e10 * M_sun; Lmw = (2 / math.pi) * math.sqrt(G * Mmw / A0)
for r in (1 * Gpc, 3 * Gpc, 8 * Gpc):
    k = 1 / r
    mo = 1 / (-math.expm1(-k * Lmw)); mg = 1 / (1 - math.cosh(k * (Lam - Lmw)) / math.cosh(k * Lam))
    P(f"     r = {r / Gpc:3.0f} Gpc: m_Λ/m_offen = {mg / mo:.4f}")


# ------------------------------------------------------------------ F  Handprobe
P("\nF  Handprobe (anderer Weg als der Code)")
P(f"   Erstes Glied n = 1 bei r = ϖ_R: (2/π)(arctan 1 − ½) = (2/π)(0,785398 − 0,5) = {(2 / math.pi) * (math.pi / 4 - 0.5):.5f}")
def glied1(rho):
    def f(x):
        z = x * rho; j = math.sin(z) / z**2 - math.cos(z) / z
        return x * j * math.exp(-x)
    return (2 / math.pi) * integrate.quad(f, 0, 80, limit=2000)[0]
P(f"   dasselbe als Fourier-Integral (2/π)∫x j₁(x)e^{{−x}}dx = {glied1(1.0):.5f}")
P("   Tiefgrenze: (2/π)∫₀^∞ j₁(x)dx = 2/π, da ∫j₁ = [−j₀]₀^∞ = 1  ⇒  a·r·L/GM → 0,63662")

open(__file__.replace('.py', '.txt'), 'w', encoding='utf-8').write("\n".join(out) + "\n")
