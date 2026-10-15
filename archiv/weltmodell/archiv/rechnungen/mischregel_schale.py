# -*- coding: utf-8 -*-
"""
Mischregel an der w-Schale (R_ta bis 1,16·R_ta) aus EINER Rückkehr mit zwei Verlustkanälen.

Kern (07): Dichte erster Rückkehrzeiten f(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/4T}, Laplace-Transformierte f̂(p) = e^{−√(pτ_D)}.
Zwei konstante Verlustraten auf demselben Rücklauf: Γ = a_N/c (Regel, Γτ_D = a_N/a₀) und θ̃ (Schließungsregel).
Je Rücklauf drei Ausgänge (konkurrierende Raten, jede Rate genau einmal):
    Rückkehr         q_t  = f̂(Γ + θ̃) = e^{−s_t},  s_t² = (Γ + θ̃)τ_D = a_N/a₀ + θ̃τ_D
    Ende über Γ      p_Γ  = ρ·(1 − q_t),          ρ = Γ/(Γ + θ̃)
    Ende über θ̃      p_θ  = (1 − ρ)·(1 − q_t)
(Beim Verlust zur Zeit T ist der Kanal unabhängig von T: P(Γ zuerst | Verlust) = ρ.)
Eine Kette von Rückläufen endet genau einmal; Zahl der Rückkehren N geometrisch, E[N] = q_t/(1 − q_t) = ν(s_t) − 1,
Endkanal unabhängig von N.
Zuordnung (aus den beiden Grenzfällen des Modells gelesen):
    Regel: Kette endet über Γ (Flucht in den offenen Bereich, Kac) → zählt alle ihre Rückkehren: E[N; Γ-Ende] = ρ(ν(s_t) − 1)
    Mode:  Rücklauf bricht schon beim ersten Mal über θ̃ ab (Schleife nie geschlossen, mitatmend) → 1 Anteil Q: p_θ·Q
    Rest:  Schleife geschlossen, dann über θ̃ gebrochen (Wahrscheinlichkeit (1 − ρ)q_t) → zählt jetzt nichts
           (im Hintergrund genau der nicht klumpende Anteil 1 − w, der mit λ = (1 − w)ω/N* übergeben wird)
    a/a_N = 1 + p_θ·Q + ρ·(ν(s_t) − 1)                                                     (A)
Grenzen: θ̃ = 0 → ρ = 1, p_θ = 0: ν(a_N).  Γ = 0 → ρ = 0, p_θ = 1 − e^{−√(θ̃τ_D)} = w: 1 + w·Q.
Zum Vergleich (Lesarten aus naht_rq_umkehr): (iv) 1 + wQ + (1 − w)(ν(s_t) − 1)  [θ̃ zweimal],
    (c) 1 + p_θQ + (ν(s_t) − 1) [alle Rückkehren zur Regel], (d) 1 + (1 − ρ)Q + ρ(ν(s_t) − 1) [ganze Kette nach Endkanal],
    (e) 1 + wQ + (1 − w)(ν(a_N) − 1) [w als Gewicht].
Danach gelesene Masse M_gel = r²a/G (Punktmasse, gelesen, keine Masse am Ort), gelesene Dichte, Zylindermasse und ΔΣ 0,5–1 Mpc.
Nur ftm_kern.
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize_scalar
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

G, c, Ms, kpc, Mpc, Gyr, pc = K.G, K.c, K.M_sun, K.kpc, K.Mpc, K.Gyr, K.kpc / 1e3
a0, Q = K.A0, K.Q_KOSM
H0, Om = K.H0_si(), K.OMEGA_M
OL = 1 - Om
tauD = c / a0
rho_m = Om * 3 * H0**2 / (8 * math.pi * G)

def hdr(t):
    print("\n" + "=" * 110 + "\n" + t + "\n" + "=" * 110)

# ------------------------------------------------------------------ Regel (A) und Vergleich
def teile(xN, x):
    """xN = a_N/a₀ = Γτ_D, x = θ̃τ_D → (q_t, ρ, p_Γ, p_θ, s_t)"""
    st = math.sqrt(xN + x)
    qt = math.exp(-st)
    rho = xN / (xN + x) if xN + x > 0 else 1.0
    return qt, rho, rho * (1 - qt), (1 - rho) * (1 - qt), st

def nu_s(s):
    return 1 / (-math.expm1(-s))

def f_A(xN, x):
    qt, rho, pG, pT, st = teile(xN, x)
    return 1 + pT * Q + rho * (nu_s(st) - 1)

def f_iv(xN, x):
    w = -math.expm1(-math.sqrt(x)); st = math.sqrt(xN + x)
    return 1 + w * Q + (1 - w) * (nu_s(st) - 1)

def f_c(xN, x):
    qt, rho, pG, pT, st = teile(xN, x)
    return 1 + pT * Q + (nu_s(st) - 1)

def f_d(xN, x):
    qt, rho, pG, pT, st = teile(xN, x)
    return 1 + (1 - rho) * Q + rho * (nu_s(st) - 1)

def f_e(xN, x):
    w = -math.expm1(-math.sqrt(x))
    return 1 + w * Q + (1 - w) * (nu_s(math.sqrt(xN)) - 1)

LES = [("(A) zwei Kanäle", f_A), ("(iv) θ̃ doppelt", f_iv), ("(c) alle Rückk.", f_c), ("(d) Kette/Endkanal", f_d), ("(e) w als Gewicht", f_e)]

# ------------------------------------------------------------------ Schalen (Top-Hat mit Λ, wie naht_rq_umkehr, hier eigener Code)
Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_von_a = lambda a: quad(lambda y: 1 / (y * Hb(y)), 1e-8, a, limit=200)[0]
t0 = t_von_a(1.0)
AI = 1e-3; TI = t_von_a(AI); RHO_I = rho_m / AI**3

def lauf(M, di, t_end, dense=False):
    ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
    vi = Hb(AI) * ri * (1 - di / 3)
    f = lambda t, y: [y[1], -G * M / y[0]**2 + OL * H0**2 * y[0]]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    return solve_ivp(f, (TI, t_end), [ri, vi], events=ev, rtol=1e-10, atol=1e-8, dense_output=dense)

def umkehr_heute(M):
    def fn(ld):
        s = lauf(M, 10**ld, 2 * t0)
        return (s.t_events[0][0] if len(s.t_events[0]) else 2 * t0) - t0
    ld = brentq(fn, -4, 0, xtol=1e-12)
    s = lauf(M, 10**ld, 2 * t0)
    Rt = s.y_events[0][0][0]
    return Rt, 3 * M / (4 * math.pi * Rt**3) / rho_m

def rechne_galaxie(Mb, n_innen=400, n_schale=90):
    M_mit = (1 + Q) * Mb
    R_ta, D_ta = umkehr_heute(M_mit)
    def M_L(r):                                   # Lagrange-Masse (Kac: Mode Q·M_b bis R_ta, außen Hintergrund)
        if r <= R_ta:
            return Mb * (1 + Q * r / R_ta)
        return M_mit + 4 / 3 * math.pi * rho_m * (r**3 - R_ta**3)
    def theta_til(r):
        """expandierende Schale, heute bei r: θ̃ = min θ der Geschichte (θ = 3V/R), 0 wenn umgekehrt"""
        M = M_L(r)
        fn = lambda ld: (lambda s: (s.y[0][-1] if s.status == 0 else -1e30) - r)(lauf(M, 10**ld, t0))
        ld = brentq(fn, -7.0, math.log10(0.6), xtol=1e-12)
        s = lauf(M, 10**ld, t0, dense=True)
        if s.status != 0:
            return 0.0
        tt = np.geomspace(TI * 10, t0, 600); Y = s.sol(tt)
        return max(float((3 * Y[1] / Y[0]).min()), 0.0)
    ua = np.concatenate([np.geomspace(1e-5, 0.16, n_schale), np.geomspace(0.17, 9.0, 40)])   # r/R_ta − 1
    TH = np.array([theta_til(R_ta * (1 + u)) * tauD for u in ua])
    return dict(Mb=Mb, R_ta=R_ta, D_ta=D_ta, ua=ua, TH=TH, M_L=M_L)

def x_von_r(g, r):
    """θ̃τ_D an r (0 innen; log-Interpolation in r/R_ta − 1)"""
    u = r / g["R_ta"] - 1
    if u <= 0:
        return 0.0
    if u < g["ua"][0]:
        return g["TH"][0] * u / g["ua"][0]            # linear an der Umkehr (θ̃ ∝ r − R_ta)
    return float(np.exp(np.interp(np.log(u), np.log(g["ua"]), np.log(np.maximum(g["TH"], 1e-300)))))

# ================================================================== 1 Herleitung und Grenzfälle
hdr("1  Herleitung: eine Rückkehr, zwei Verlustkanäle (Γ = a_N/c, θ̃), jede Rate einmal")
print("   E[e^{−(Γ+θ̃)T}] = e^{−s_t}; P(Verlust über θ̃ im Rücklauf) = θ̃/(Γ+θ̃)·(1 − e^{−s_t}) = p_θ;  P(über Γ) = ρ(1 − e^{−s_t})")
print("   Prüfung der Aufteilung durch Integration über f(T): Σ = q_t + p_Γ + p_θ, und p_θ gegen ∫ f(T)·θ̃/(Γ+θ̃)(1 − e^{−(Γ+θ̃)T}) dT")
for xN, x in ((2.26e-4, 0.0112), (0.03, 0.5), (1.0, 1.0)):
    qt, rho, pG, pT, st = teile(xN, x)
    # Integration in u = τ_D/T mit T in Einheiten τ_D
    fT = lambda T: math.sqrt(1 / (4 * math.pi)) * T**-1.5 * math.exp(-1 / (4 * T))
    num_q = quad(lambda T: fT(T) * math.exp(-(xN + x) * T), 0, np.inf, limit=400)[0]
    num_t = quad(lambda T: fT(T) * x / (xN + x) * (-math.expm1(-(xN + x) * T)), 0, np.inf, limit=400)[0]
    print(f"   a_N/a₀ = {xN:<8g} θ̃τ_D = {x:<6g}: q_t = {qt:.6f} (Integral {num_q:.6f}), p_θ = {pT:.6f} (Integral {num_t:.6f}), p_Γ = {pG:.6f}, Σ = {qt+pG+pT:.6f}")
print("   Grenzfälle von (A):")
for xN in (1.0, 0.03, 2.26e-4):
    print(f"      θ̃ = 0, a_N/a₀ = {xN:<8g}: (A) = {f_A(xN, 0.0):.6f}, ν(a_N) = {float(K.nu(xN*a0)):.6f}")
for x in (0.0112, 0.48045, 17.8):
    w = K.gewicht_w(x / tauD)
    print(f"      Γ = 0, θ̃τ_D = {x:<7g}: (A) = {f_A(0.0, x):.6f}, 1 + w·Q = {1 + w*Q:.6f} (w = {w:.4f})")
print("   Γ → 0 (a_N/a₀ = 10⁻¹²) in allen Lesarten (Soll 1 + wQ): " + ", ".join(f"{n} {f(1e-12, 17.8):.4f}" for n, f in LES) +
      f"  | Soll {1 + K.gewicht_w(17.8/tauD)*Q:.4f}")

# ================================================================== 2 Mulde, a_N → 0 analytisch, Milchstraße
hdr("2  Mulde: Minimum von a/a_N über θ̃ bei festem a_N (in der Schale wächst θ̃ von 0 an, a_N fast fest)")
print("   klein s_t: ρ(ν_t − 1) ≈ s_N²/s_t³, p_θQ ≈ s_t·Q  ⇒  s_t* = (3s_N²/Q)^{1/4}, Minimum ≈ 1 + (4/3)·Q·s_t* ∝ s_N^{1/2}; a_N → 0: Minimum → 1")
print("   a_N/a₀      s_N      " + "  ".join(f"{n:>18s}" for n, _ in LES) + "   Näherung (A)   θ̃τ_D bei Min (A)   w dort")
for xN in (1e-2, 1e-3, 2.2589e-4, 1e-4, 1e-5, 1e-6):
    mins = []
    for n, f in LES:
        rr = minimize_scalar(lambda lx: f(xN, 10**lx), bounds=(-9, 2), method="bounded", options={"xatol": 1e-9})
        mins.append((rr.fun, 10**rr.x))
    sN = math.sqrt(xN); ss = (3 * xN / Q)**0.25
    print(f"   {xN:8.2e}  {sN:7.4f}  " + "  ".join(f"{m:8.3f} @{x_:8.2e}" for m, x_ in mins) +
          f"   {1 + 4/3*Q*ss:8.3f}      {mins[0][1]:9.3e}        {K.gewicht_w(mins[0][1]/tauD):.4f}")

# Milchstraße
g = rechne_galaxie(6e10 * Ms)
R_ta, Mb = g["R_ta"], g["Mb"]
print(f"\n   Milchstraße M_b = 6·10¹⁰ M☉: R_ta = {R_ta/kpc:.2f} kpc (1 + Δ_ta = {g['D_ta']:.3f}), a_N/a₀(R_ta) = {G*Mb/R_ta**2/a0:.4e}")
print("   θ̃τ_D der Schale nahe R_ta (linear?): " + ", ".join(f"u = {u:.1e}: {t:.4g} ({t/u:.2f}·u)" for u, t in zip(g["ua"][:60:12], g["TH"][:60:12])))
print("   r[kpc]     r/R_ta−1   θ̃τ_D      w       ρ          q_t      p_θ     " + "  ".join(f"{n:>18s}" for n, _ in LES) + "   ν(a_N)")
tab_u = [-0.3, -1e-6, 1e-5, 3e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.16, 0.25, 0.5, 0.72, 1.0, 2.0, 4.0, 7.6]
for u in tab_u:
    r = R_ta * (1 + u); xN = G * Mb / r**2 / a0; x = x_von_r(g, r)
    qt, rho, pG, pT, st = teile(xN, x)
    print(f"   {r/kpc:8.2f}  {u:+9.1e}  {x:9.3e}  {K.gewicht_w(x/tauD):.4f}  {rho:9.3e}  {qt:.4f}  {pT:.4f}  " +
          "  ".join(f"{f(xN, x):18.3f}" for n, f in LES) + f"   {float(K.nu(xN*a0)):6.2f}")
RR = R_ta * (1 + np.geomspace(1e-6, 8, 6000))
FA = np.array([f_A(G * Mb / r**2 / a0, x_von_r(g, r)) for r in RR])
k = int(np.argmin(FA)); xk = x_von_r(g, RR[k])
print(f"   (A) Minimum a/a_N = {FA[k]:.3f} bei r = {RR[k]/kpc:.3f} kpc (r/R_ta − 1 = {RR[k]/R_ta-1:.2e}, θ̃τ_D = {xk:.4f}, w = {K.gewicht_w(xk/tauD):.3f})")
for lvl in (1 + Q / 2, 1 + 0.9 * Q, 1 + 0.95 * Q * K.gewicht_w(3 * H0)):
    j = k + int(np.argmax(FA[k:] >= lvl))
    print(f"      a/a_N erreicht {lvl:.3f} wieder bei r = {RR[j]/kpc:.1f} kpc ({RR[j]/R_ta:.3f}·R_ta)")
print(f"      weit außen → 1 + w(3H₀)·Q = {1 + K.gewicht_w(3*H0)*Q:.3f}")

# ================================================================== 3 gelesene Masse und Linse
hdr("3  Gelesene Masse M_gel(<r) = r²a/G = M_b·a/a_N (gelesen, keine Masse am Ort) und Linsenprofil (Kugel: Licht sieht die ganze Regel)")
def a_rel(r, form=f_A):
    xN = G * Mb / r**2 / a0
    return form(xN, x_von_r(g, r))
# Gitter: dicht an R_ta
rg = np.unique(np.concatenate([np.geomspace(0.5 * kpc, 0.999 * R_ta, 3000), R_ta * (1 + np.geomspace(1e-7, 1e-2, 3000)),
                               R_ta * (1 + np.geomspace(1.01e-2, 12, 3000)), [R_ta]]))
MG = {n: Mb * np.array([a_rel(r, f) for r in rg]) for n, f in LES}
MG["(i) Regel bis R_ta, dann 1+Q"] = Mb * np.where(rg <= R_ta, K.nu(G * Mb / rg**2), 1 + Q)
Mneg = MG["(A) zwei Kanäle"]
dM = np.diff(Mneg)
neg = dM < 0
print(f"   (A): M_gel(R_ta) = {Mneg[np.searchsorted(rg, R_ta)]/Mb:.2f}·M_b; gelesene Dichte negativ in r = {rg[:-1][neg].min()/kpc:.3f}–{rg[1:][neg].max()/kpc:.3f} kpc, "
      f"Summe der negativen gelesenen Masse {dM[neg].sum()/Mb:.2f}·M_b = {dM[neg].sum()/Ms:.2e} M☉ (gelesen)")
pos_out = dM[(~neg) & (rg[:-1] > R_ta)]
print(f"   danach positiv bis außen: +{pos_out.sum()/Mb:.2f}·M_b; M_gel(12·R_ta) = {Mneg[-1]/Mb:.3f}·M_b")
for rk in (300, 500, 570, 590, 620, 700, 800, 1000, 1500, 3000):
    print(f"      r = {rk:5d} kpc: M_gel/M_b = " + "  ".join(f"{n.split()[0]} {np.interp(rk*kpc, rg, MG[n])/Mb:7.3f}" for n in MG))

def zylinder(M):
    """M_cyl(R) = M(R) + ∫_R^∞ M'(r)(1 − √(1 − R²/r²)) dr  (Rand: M' = 0 jenseits des Gitters); Σ = dM_cyl/dR /(2πR)"""
    Rp = np.geomspace(100 * kpc, 1500 * kpc, 281)
    Mc = []
    rm = 0.5 * (rg[1:] + rg[:-1]); dm = np.diff(M)
    for R in Rp:
        sel = rm > R
        Mc.append(np.interp(R, rg, M) + np.sum(dm[sel] * (1 - np.sqrt(1 - (R / rm[sel])**2))))
    Mc = np.array(Mc)
    Sig = np.gradient(Mc, Rp) / (2 * math.pi * Rp)
    return Rp, Mc, Sig, Mc / (math.pi * Rp**2) - Sig
print("   Linse (M☉/pc²), Zylindermasse gelesen; ΔΣ = Σ̄(<R) − Σ(R):")
print("   R[kpc]   " + "  ".join(f"{n:>30s}" for n in ("(A) Σ | ΔΣ", "(i) Σ | ΔΣ", "(e) Σ | ΔΣ", "(iv) Σ | ΔΣ")))
ZY = {n: zylinder(MG[n]) for n in ("(A) zwei Kanäle", "(i) Regel bis R_ta, dann 1+Q", "(e) w als Gewicht", "(iv) θ̃ doppelt")}
Rp = ZY["(A) zwei Kanäle"][0]
for Rk in (200, 400, 500, 550, 570, 578, 590, 650, 750, 850, 1000, 1200):
    j = int(np.argmin(abs(Rp - Rk * kpc)))
    print(f"   {Rp[j]/kpc:7.1f}  " + "  ".join(f"{ZY[n][2][j]*pc**2/Ms:+13.3f} | {ZY[n][3][j]*pc**2/Ms:+13.3f}" for n in ZY))
S = ZY["(A) zwei Kanäle"]
band = (Rp >= 500 * kpc) & (Rp <= 1000 * kpc)
print(f"   (A) im Band 0,5–1 Mpc: Σ negativ für R = {Rp[band & (S[2] < 0)].min()/kpc:.0f}–{Rp[band & (S[2] < 0)].max()/kpc:.0f} kpc; "
      f"ΔΣ(500) = {np.interp(500*kpc, Rp, S[3])*pc**2/Ms:.2f}, ΔΣ(1000) = {np.interp(1000*kpc, Rp, S[3])*pc**2/Ms:.3f} M☉/pc²")
print(f"   Vergleich Punktmasse (1+Q)M_b: ΔΣ = M/(πR²) → 500 kpc {(1+Q)*Mb/(math.pi*(500*kpc)**2)*pc**2/Ms:.3f}, 1000 kpc {(1+Q)*Mb/(math.pi*(1000*kpc)**2)*pc**2/Ms:.3f} M☉/pc²")

# ================================================================== 4 andere Massen
hdr("4  Andere Massen: Tiefe und Lage der Mulde (A)")
for mb in (1e9, 1e10, 1e11, 1e12):
    gg = rechne_galaxie(mb * Ms, n_schale=50)
    RRg = gg["R_ta"] * (1 + np.geomspace(1e-7, 1, 3000))
    F = np.array([f_A(G * mb * Ms / r**2 / a0, x_von_r(gg, r)) for r in RRg])
    kk = int(np.argmin(F)); xNt = G * mb * Ms / gg["R_ta"]**2 / a0
    print(f"   M_b = {mb:.0e} M☉: R_ta = {gg['R_ta']/kpc:7.1f} kpc, a_N/a₀(R_ta) = {xNt:.3e}, ν(R_ta) = {float(K.nu(xNt*a0)):6.1f}; "
          f"Minimum {F[kk]:.3f} bei r/R_ta − 1 = {RRg[kk]/gg['R_ta']-1:.1e} (Näherung {1 + 4/3*Q*(3*xNt/Q)**0.25:.3f})")

# ================================================================== 5 von Hand
hdr("5  Prüfung von Hand (Zahlen eingetippt, ohne ftm_kern)")
Qh = 0.1200 / 0.02237
xN, x = 2.2589e-4, 0.186
st = math.sqrt(xN + x); qt = math.exp(-st); rho = xN / (xN + x)
print(f"   583 kpc (θ̃τ_D = 0,186 aus naht_rq_umkehr): s_t = {st:.4f}, q_t = {qt:.4f}, ρ = {rho:.2e}, ν_t − 1 = {qt/(1-qt):.3f}")
print(f"      (A) = 1 + (1 − ρ)(1 − q_t)Q + ρ(ν_t − 1) = 1 + {(1-rho)*(1-qt)*Qh:.4f} + {rho*qt/(1-qt):.5f} = {1 + (1-rho)*(1-qt)*Qh + rho*qt/(1-qt):.4f};  Skript {f_A(xN, x):.4f}")
ss = (3 * xN / Qh)**0.25
print(f"   Minimum bei a_N/a₀ = 2,2589·10⁻⁴: s_t* = (3·{xN:.4e}/{Qh:.4f})^¼ = {ss:.4f}, 1 + (4/3)Q s_t* = {1 + 4/3*Qh*ss:.4f}; "
      f"genauer mit ν − 1 = 1/s − ½ + s/12: {min(1 + (1-xN/s**2)*(1-math.exp(-s))*Qh + xN/s**2*(1/s-0.5+s/12) for s in np.linspace(0.05,0.3,2501)):.4f}")
print(f"   Γ = 0, θ̃τ_D = (ln 2)²: p_θ = 1 − e^(−ln 2) = 0,5 ⇒ (A) = 1 + Q/2 = {1 + Qh/2:.4f} (Skript {f_A(0, math.log(2)**2):.4f})")
vfh = (6e10 * 1.32712440018e20 * 1.1e-10)**0.25; pch = 3.0856775814913673e16; Msh = 1.98841e30
print(f"   Linse innen (Tiefbereich, isotherm M ∝ r): ΔΣ = v_f²/(4GR), v_f = (GMa₀)^¼ = {vfh/1e3:.1f} km/s: 200 kpc {vfh**2/(4*6.6743e-11*200e3*pch)*pch**2/Msh:.2f} M☉/pc² "
      f"(Skript mit Schale {np.interp(200*kpc, Rp, S[3])*pc**2/Ms:.2f}; Rest aus dem Kern r < r* und der negativen Schale)")
