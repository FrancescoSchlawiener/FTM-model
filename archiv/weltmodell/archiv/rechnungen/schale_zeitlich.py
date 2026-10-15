# -*- coding: utf-8 -*-
"""
Die w-Schale zeitlich: Jede Lagrange-Schale durchläuft Ausdehnung, Umkehr und Kollaps; der Anteil der Ketten, der
momentan nicht zählt, wird mit λ = (1 − w)ω/N* übergeben. Füllt die Übergabe die Mulde der Momentaufnahme (A) bis heute?

Momentaufnahme (A) je Schale (eine Rückkehr, zwei Verlustkanäle Γ = a_N/c, θ̃):
    s_t² = a_N/a₀ + θ̃τ_D, q_t = e^{−s_t}, ρ = Γ/(Γ + θ̃), p_θ = (1 − ρ)(1 − q_t) (Mode), Rest R = (1 − ρ)q_t
    Ketten nach Endkanal: Γ-Ende (Wahrsch. ρ, je Kette ν_t − 1 Rückkehren), Mode (p_θ, 0 Rückkehren, zählt Q),
    Rest (R: geschlossen, dann über θ̃ gebrochen; zählt bis zur Übergabe nichts)
Zeitlich, P = Anteil der Ketten, die noch nicht übergeben sind (zählen nichts); übergebene Ketten zählen wie Γ-Ketten:
    a/a_N = 1 + p_θ·Q + (1 − p_θ − P)·(ν_t − 1)                                  (T)
    P = R: (A).  P = 0: alles Geschlossene zählt.  Nach der Umkehr θ̃ = 0 (Preisach): p_θ = 0, ν_t = ν(a_N), a/a_N = 1 + (1 − P)(ν − 1).
Zwei Lesarten des Gedächtnisses:
    M1 (ohne Gedächtnis): P = ξ·R_jetzt, ξ = e^{−I/N*}, I = ∫(1 − w)ω dt  → nach der Umkehr R = 0, sofort ν(a_N)
    M2 (gebrochen bleibt gebrochen): dP/dt = max(dR/dt, 0) − λP, λ = (1 − w(θ̃))ω/N*; nach der Umkehr zerfällt P mit ω/N*
Schalen: Top-Hat mit Λ (wie mischregel_schale), Kac-Lagrange-Masse (Mode Q·M_b bis R_ta gleichmäßig in r, außen Hintergrund),
    nach der Umkehr Newton-Einfall bis R_1/2, dann fest (wie naht_rq_umkehr); ω² = GM_L/R³; a₀ fest.
a_N = GM_b/r² am heutigen Ort (Punktmasse der Baryonen), gelesene Masse M_gel = r²a/G (gelesen, keine Masse am Ort).
Nur ftm_kern.
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

OUT = []
def P_(s=""):
    OUT.append(s); print(s, flush=True)

G, c, Ms, kpc, Gyr = K.G, K.c, K.M_sun, K.kpc, K.Gyr
pc = K.pc
a0, Q = K.A0, K.Q_KOSM
H0, Om = K.H0_si(), K.OMEGA_M
OL = 1 - Om
tauD = K.tau_D(a0)
rho_m = K.rho_m_mitbewegt()
NSTAR = (0.3, 1.0, 2.4)

def hdr(t):
    P_("\n" + "=" * 110); P_(t); P_("=" * 110)

Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_von_a = lambda a: quad(lambda y: 1 / (y * Hb(y)), 1e-8, a, limit=200)[0]
t0 = t_von_a(1.0)
AI = 1e-3; TI = t_von_a(AI); RHO_I = rho_m / AI**3

def rhs(t, y, M):
    return [y[1], -G * M / y[0]**2 + OL * H0**2 * y[0]]

def start(M, di):
    ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
    return ri, Hb(AI) * ri * (1 - di / 3)

def R_heute(M, di):
    """Radius heute: Ausdehnung, Umkehr, Einfall bis R_1/2, dann fest."""
    ri, vi = start(M, di)
    ev = lambda t, y, *a: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(rhs, (TI, t0), [ri, vi], args=(M,), events=ev, rtol=1e-11, atol=1e-8)
    if not len(s.t_events[0]):
        return s.y[0][-1]
    t1, R1 = s.t_events[0][0], s.y_events[0][0][0]
    ev2 = lambda t, y, *a: y[0] - R1 / 2; ev2.terminal = True; ev2.direction = -1
    s2 = solve_ivp(rhs, (t1, t0), [R1, -1e-9], args=(M,), events=ev2, rtol=1e-11, atol=1e-8)
    return R1 / 2 if len(s2.t_events[0]) else s2.y[0][-1]

def teile(xN, x):
    st = math.sqrt(xN + x); qt = math.exp(-st)
    rho = xN / (xN + x) if xN + x > 0 else 1.0
    return qt, rho, (1 - rho) * (1 - qt), (1 - rho) * qt, st

def nu_s(s):
    return 1 / (-math.expm1(-s))

def formel_T(xN, x, Pn):
    qt, rho, pT, R, st = teile(xN, x)
    return 1 + pT * Q + (1 - pT - Pn) * (nu_s(st) - 1)

def geschichte(M, di, Mb):
    """Volle Geschichte einer Schale: θ̃(t0)τ_D, I (N*-frei), P_M1/P_M2 je N*, t_Umkehr."""
    ri, vi = start(M, di)
    ev = lambda t, y, *a: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(rhs, (TI, t0), [ri, vi], args=(M,), events=ev, rtol=1e-11, atol=1e-8, dense_output=True)
    gedreht = len(s.t_events[0]) > 0
    t1 = s.t_events[0][0] if gedreht else t0
    # Gitter: logarithmisch in t und dicht vor t1
    tt = np.unique(np.concatenate([np.geomspace(TI, t1, 3000), t1 - np.geomspace(1e-9 * t1, 0.5 * (t1 - TI), 3000)]))
    tt = tt[(tt >= TI) & (tt <= t1)]
    Y = s.sol(tt); Rr, V = Y[0], Y[1]
    th = np.maximum(3 * V / Rr, 0.0)                     # θ fällt bis zur Umkehr monoton: θ̃ = θ
    x = th * tauD; xN = G * Mb / Rr**2 / a0
    st = np.sqrt(xN + x); qt = np.exp(-st); rho = xN / (xN + x)
    Rest = (1 - rho) * qt
    omw = np.exp(-np.sqrt(x))                            # 1 − w(θ̃) = e^{−√(θ̃τ_D)}
    om = np.sqrt(G * M / Rr**3)
    dt = np.diff(tt); dR = np.maximum(np.diff(Rest), 0.0)
    info = dict(gedreht=gedreht, t1=t1, x0=float(x[-1]) if not gedreht else 0.0, Rest0=float(Rest[-1]) if not gedreht else 0.0,
                Rmax=float(Rest.max()))
    P2 = {}
    for var, lam0 in (("a", omw * om), ("b", om)):      # λ·N*: a = (1 − w)ω (wie 03), b = ω (Rest ist schon der Anteil 1 − w)
        lm = 0.5 * (lam0[1:] + lam0[:-1])
        info["I_vor" + var] = float(np.sum(lm * dt))
        for Ns in NSTAR:
            e = np.exp(-lm * dt / Ns)
            # p_{k+1} = p_k e_k + dR_k √e_k  ⇒  p_end = Σ dR_k √e_k Π_{j>k} e_j
            cum = np.concatenate([np.cumsum(np.log(e)[::-1])[::-1][1:], [0.0]])
            P2[(var, Ns)] = float(np.sum(dR * np.sqrt(e) * np.exp(cum)))
    I = info["I_vora"]; info["I_vor"] = I
    I_nach = 0.0
    if gedreht:
        R1 = s.y_events[0][0][0]
        ev2 = lambda t, y, *a: y[0] - R1 / 2; ev2.terminal = True; ev2.direction = -1
        s2 = solve_ivp(rhs, (t1, t0), [R1, -1e-9], args=(M,), events=ev2, rtol=1e-11, atol=1e-8, dense_output=True)
        t2 = s2.t_events[0][0] if len(s2.t_events[0]) else t0
        u = np.linspace(t1, t2, 4000)
        I_nach = float(np.trapezoid(np.sqrt(G * M / s2.sol(u)[0]**3), u))
        if t2 < t0:
            I_nach += (t0 - t2) * math.sqrt(G * M / (R1 / 2)**3)
        info.update(R1=R1, t2=t2)
    info["I_nach"] = I_nach
    info["xi"] = {Ns: math.exp(-(I + I_nach) / Ns) for Ns in NSTAR}
    info["P1"] = {Ns: (math.exp(-I / Ns) * info["Rest0"] if not gedreht else 0.0) for Ns in NSTAR}
    info["P2"] = {Ns: P2[("a", Ns)] * math.exp(-I_nach / Ns) for Ns in NSTAR}
    info["P1b"] = {Ns: (math.exp(-info["I_vorb"] / Ns) * info["Rest0"] if not gedreht else 0.0) for Ns in NSTAR}
    info["P2b"] = {Ns: P2[("b", Ns)] * math.exp(-I_nach / Ns) for Ns in NSTAR}
    return info

def umkehr_heute(M):
    def fn(ld):
        ri, vi = start(M, 10**ld)
        ev = lambda t, y, *a: y[1]; ev.terminal = True; ev.direction = -1
        s = solve_ivp(rhs, (TI, 2 * t0), [ri, vi], args=(M,), events=ev, rtol=1e-11, atol=1e-8)
        return (s.t_events[0][0] if len(s.t_events[0]) else 2 * t0) - t0
    ld = brentq(fn, -4, 0, xtol=1e-13)
    ri, vi = start(M, 10**ld)
    ev = lambda t, y, *a: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(rhs, (TI, 2 * t0), [ri, vi], args=(M,), events=ev, rtol=1e-11, atol=1e-8)
    Rt = s.y_events[0][0][0]
    return Rt, 3 * M / (4 * math.pi * Rt**3) / rho_m, ld

def galaxie(Mb, U_in, U_out):
    M_mit = (1 + Q) * Mb
    R_ta, D_ta, ld_ta = umkehr_heute(M_mit)
    def M_L(r):
        if r <= R_ta:
            return Mb * (1 + Q * r / R_ta)
        return M_mit + 4 / 3 * math.pi * rho_m * (r**3 - R_ta**3)
    rows = []
    for u in list(-np.array(U_in)) + list(U_out):
        r = R_ta * (1 + u); M = M_L(r)
        fn = lambda ld: R_heute(M, 10**ld) - r
        if u > 0:
            ld = brentq(fn, -7.0, ld_ta + 1e-9, xtol=1e-13)
        else:
            ld = brentq(fn, ld_ta - 1e-9, 0.5, xtol=1e-13)
        g = geschichte(M, 10**ld, Mb)
        g.update(u=u, r=r, xN=G * Mb / r**2 / a0)
        rows.append(g)
    return dict(Mb=Mb, R_ta=R_ta, D_ta=D_ta, rows=rows)

def werte(g, Ns):
    """a/a_N je Schale: (A), M1, M2"""
    out = []
    for w in g["rows"]:
        xN, x = w["xN"], w["x0"]
        if w["gedreht"]:
            nuN = nu_s(math.sqrt(xN))
            A = nuN; m1 = nuN; m2 = 1 + (1 - w["P2"][Ns]) * (nuN - 1); m1b = nuN; m2b = 1 + (1 - w["P2b"][Ns]) * (nuN - 1)
        else:
            qt, rho, pT, R, st = teile(xN, x)
            A = formel_T(xN, x, R); m1 = formel_T(xN, x, w["P1"][Ns]); m2 = formel_T(xN, x, min(w["P2"][Ns], 1 - pT))
            m1b = formel_T(xN, x, w["P1b"][Ns]); m2b = formel_T(xN, x, min(w["P2b"][Ns], 1 - pT))
        out.append((w["r"], A, m1, m2, m1b, m2b))
    return np.array(out)

U_OUT = np.concatenate([np.geomspace(1e-5, 2e-2, 34), np.geomspace(2.5e-2, 3.0, 22)])
U_IN = np.concatenate([np.geomspace(1e-5, 2e-2, 14), np.geomspace(3e-2, 0.9, 22)])

GAL = {}
for mb in (1e10, 6e10, 1e12):
    GAL[mb] = galaxie(mb * Ms, U_IN, U_OUT)

# ================================================================== 1 Geschichte an der Schale
hdr("1  Geschichte der Schalen (Milchstraße 6·10¹⁰ M☉): Übergabe vor und nach der Umkehr")
g = GAL[6e10]
P_(f"   R_ta = {g['R_ta']/kpc:.2f} kpc, 1 + Δ_ta = {g['D_ta']:.3f}")
P_("   I_vor: λ·N* = (1 − w)ω (a) | ω (b);  P_M2 = noch nicht übergebener Anteil (a | b), N* = 0,3 / 1 / 2,4")
P_("   r/R_ta−1    r[kpc]    θ̃τ_D     Rest R   max R    I_vor a | b         I_nach   P_M2 a (0,3|1|2,4)      P_M2 b (0,3|1|2,4)     t_Umk.−t0 [Mrd. J.]")
for w in g["rows"]:
    if abs(w["u"]) in (1e-5,) or any(abs(abs(w["u"]) - v) / v < 0.08 for v in (3e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.9, 1.0, 3.0)):
        P_(f"   {w['u']:+9.1e}  {w['r']/kpc:8.2f}  {w['x0']:8.2e}  {w['Rest0']:.4f}  {w['Rmax']:.4f}  {w['I_vora']:.4f} | {w['I_vorb']:7.4f}  {w['I_nach']:7.3f}  "
           + " ".join(f"{w['P2'][n]:.3f}" for n in NSTAR) + "     " + " ".join(f"{w['P2b'][n]:.3f}" for n in NSTAR)
           + f"      {(w['t1']-t0)/Gyr:+.3f}")

# ================================================================== 2 Mulde, Tiefe und Breite
hdr("2  Mulde heute: a/a_N je Lesart; Tiefe (Minimum außen), Breite (a/a_N < 1 + Q/2) und Abfall der gelesenen Masse")
def auswerten(g, Ns, spalte):
    V = werte(g, Ns); r = V[:, 0]; f = V[:, spalte]
    Mb = g["Mb"]; R_ta = g["R_ta"]
    o = np.argsort(r); r, f = r[o], f[o]
    aus = r > R_ta
    k = np.argmin(np.where(aus, f, np.inf))
    tief = r[f < 1 + Q / 2]
    breite = (tief.min() / kpc, tief.max() / kpc) if tief.size else None
    Mg = f                                     # M_gel/M_b = a/a_N
    j = int(np.argmax(Mg)); jr = r[j]
    return dict(min=f[k], rmin=r[k], breite=breite, Mmax=Mg[j], rmax=jr, M_aus=Mg[-1],
                abfall=Mg[j] - Mg[j:].min(), r_ab_end=r[j + int(np.argmin(Mg[j:]))], r=r, f=f)
for mb, g in GAL.items():
    P_(f"\n   M_b = {mb:.0e} M☉: R_ta = {g['R_ta']/kpc:.1f} kpc, a_N/a₀(R_ta) = {G*g['Mb']/g['R_ta']**2/a0:.3e}, ν(R_ta) = {nu_s(math.sqrt(G*g['Mb']/g['R_ta']**2/a0)):.1f}")
    for lab, sp in (("(A) Momentaufnahme", 1), ("M1a ohne Gedächtnis", 2), ("M2a gebrochen bleibt", 3), ("M1b (λ = ω/N*)", 4), ("M2b (λ = ω/N*)", 5)):
        for Ns in (NSTAR if sp > 1 else (1.0,)):
            e = auswerten(g, Ns, sp)
            b = f"{e['breite'][0]:.1f}–{e['breite'][1]:.1f} kpc (Breite {e['breite'][1]-e['breite'][0]:.1f})" if e["breite"] else "keine"
            P_(f"      {lab:20s} N* = {Ns:3.1f}: Minimum außen {e['min']:6.2f} bei {e['rmin']/kpc:8.2f} kpc ({e['rmin']/g['R_ta']:.4f}·R_ta); "
               f"a/a_N < 1 + Q/2: {b}; M_gel max {e['Mmax']:6.2f}·M_b bei {e['rmax']/kpc:6.0f} kpc, fällt um {e['abfall']:6.2f}·M_b bis {e['r_ab_end']/kpc:6.0f} kpc; "
               f"bei {U_OUT[-1]+1:.0f}·R_ta {e['M_aus']:.3f}")

P_("\n   Profil Milchstraße a/a_N (= M_gel/M_b), N* = 1 (M1 | M2) und N* = 0,3 / 2,4 (M2):")
g = GAL[6e10]
V = {Ns: werte(g, Ns) for Ns in NSTAR}
P_("   r[kpc]    r/R_ta   ν(a_N)     (A)      M1a(1)    M2a(0,3)  M2a(1)   M2a(2,4)   M1b(1)   M2b(0,3)  M2b(1)   M2b(2,4)")
for i, w in enumerate(g["rows"]):
    if i % 3 == 0:
        P_(f"   {w['r']/kpc:8.2f}  {w['r']/g['R_ta']:.5f}  {nu_s(math.sqrt(w['xN'])):7.2f}  {V[1.0][i,1]:7.2f}  {V[1.0][i,2]:8.2f}  "
           f"{V[0.3][i,3]:8.2f}  {V[1.0][i,3]:8.2f}  {V[2.4][i,3]:8.2f}   {V[1.0][i,4]:7.2f}  {V[0.3][i,5]:8.2f}  {V[1.0][i,5]:8.2f}  {V[2.4][i,5]:8.2f}")

# ================================================================== 3 Linse
hdr("3  Linse 0,5–1 Mpc (Milchstraße, Kugel, Licht sieht die ganze gelesene Masse): ΔΣ = Σ̄(<R) − Σ(R) in M☉/pc²")
def M_gitter(g, Ns, sp):
    Mb, R_ta = g["Mb"], g["R_ta"]
    V = werte(g, Ns); r = V[:, 0]; f = V[:, sp]
    o = np.argsort(r); r, f = r[o], f[o]
    rin = np.geomspace(1 * kpc, r[0], 400)[:-1]
    fin = np.array([nu_s(math.sqrt(G * Mb / x**2 / a0)) for x in rin])
    rr = np.concatenate([rin, r]); ff = np.concatenate([fin, f])
    rg = np.unique(np.concatenate([rr, np.geomspace(rr[0], rr[-1], 6000)]))
    fg = np.interp(np.log(rg), np.log(rr), ff)
    return rg, Mb * fg

def zylinder(rg, M, Rp):
    rm = 0.5 * (rg[1:] + rg[:-1]); dm = np.diff(M)
    Mc = np.array([np.interp(R, rg, M) + np.sum(dm[rm > R] * (1 - np.sqrt(1 - (R / rm[rm > R])**2))) for R in Rp])
    Sig = np.gradient(Mc, Rp) / (2 * math.pi * Rp)
    return Mc / (math.pi * Rp**2) - Sig, Sig

Rp = np.geomspace(150 * kpc, 1500 * kpc, 461)
FAELLE = [("(A)", 1.0, 1), ("M1a N*=1", 1.0, 2), ("M2a N*=0,3", 0.3, 3), ("M2a N*=1", 1.0, 3), ("M2a N*=2,4", 2.4, 3),
          ("M1b N*=1", 1.0, 4), ("M2b N*=0,3", 0.3, 5), ("M2b N*=1", 1.0, 5), ("M2b N*=2,4", 2.4, 5)]
DS = {}
for lab, Ns, sp in FAELLE:
    rg, M = M_gitter(g, Ns, sp)
    DS[lab] = zylinder(rg, M, Rp)
pm = lambda R, m: m * g["Mb"] / (math.pi * R**2) * pc**2 / Ms
P_("   R[kpc]   " + "  ".join(f"{l:>11s}" for l, _, _ in FAELLE) + "   Punktm.(1+Q)M_b   Regel ν bis R, isotherm v_f²/4GR")
vf = K.v_flach(g["Mb"])
for Rk in (200, 300, 400, 500, 550, 570, 578, 585, 600, 650, 700, 800, 1000, 1200):
    j = int(np.argmin(abs(Rp - Rk * kpc)))
    P_(f"   {Rp[j]/kpc:7.1f}  " + "  ".join(f"{DS[l][0][j]*pc**2/Ms:11.3f}" for l, _, _ in FAELLE)
       + f"   {pm(Rp[j], 1+Q):10.3f}        {vf**2/(4*G*Rp[j])*pc**2/Ms:8.3f}")
P_("   Sprung von ΔΣ über R_ta (Wert bei 0,9·R_ta / 1,1·R_ta):")
for l, _, _ in FAELLE:
    a_ = np.interp(0.9 * g["R_ta"], Rp, DS[l][0]); b_ = np.interp(1.1 * g["R_ta"], Rp, DS[l][0])
    P_(f"      {l:11s}: {a_*pc**2/Ms:7.3f} → {b_*pc**2/Ms:7.3f} M☉/pc² (Faktor {a_/b_:5.1f}); ΔΣ(1 Mpc)/Punktmasse(1+Q)M_b = "
       f"{np.interp(1000*kpc, Rp, DS[l][0])/(pm(1000*kpc, 1+Q)*Ms/pc**2):.3f}")

# ================================================================== 4 von Hand
hdr("4  Von Hand (eingetippte Zahlen, ohne ftm_kern)")
H0h = 67.4e3 / 3.0856775814913673e22; tDh = 299792458.0 / 1.1e-10
w_ = GAL[6e10]["rows"]; iu = int(np.argmin([abs(r['u'] - 1e-5) if r['u'] > 0 else 9 for r in w_]))
om2 = 11.816 * 0.3134 / 2; ome2 = om2 - 0.6866
Ih = 2 * math.sqrt(om2) / (3 * ome2 * H0h * tDh)
P_(f"   I_vor an R_ta: nahe der Umkehr θ ≈ 3ω_e²Δt, ∫e^{{−√(kΔt)}}dΔt = 2/k ⇒ I ≈ 2ω/(3ω_e²τ_D); ω² = 11,816·0,3134/2·H₀² = {om2:.4f}H₀², "
   f"ω_e² = ω² − Ω_Λ H₀² = {ome2:.4f}H₀², H₀τ_D = {H0h*tDh:.3f}: I ≈ {Ih:.4f} (Skript an 1,00001·R_ta {w_[iu]['I_vora']:.4f}; lineares θ überschätzt weiter außen)")
xh, xNh, Ph = w_[iu]["x0"], w_[iu]["xN"], w_[iu]["P2"][1.0]
sth = math.sqrt(xNh + xh); qth = math.exp(-sth); rh = xNh / (xNh + xh); pth = (1 - rh) * (1 - qth)
P_(f"   Formel (T) an 1,00001·R_ta, N* = 1: x_N = {xNh:.4e}, θ̃τ_D = {xh:.4e}, s_t = {sth:.5f}, q_t = {qth:.5f}, ρ = {rh:.4f}, p_θ = {pth:.5f}, P = {Ph:.4f}")
P_(f"      1 + p_θ·5,3643 + (1 − p_θ − P)·q_t/(1 − q_t) = {1 + pth*5.3643 + (1-pth-Ph)*qth/(1-qth):.3f} (Skript M2 {werte(GAL[6e10], 1.0)[iu,3]:.3f});"
   f"  mit P = R = {(1-rh)*qth:.4f}: {1 + pth*5.3643 + rh*qth/(1-qth):.3f} (Skript (A) {werte(GAL[6e10], 1.0)[iu,1]:.3f})")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
