# -*- coding: utf-8 -*-
"""
Kollaps nach der Umkehr: Ein Bereich atmet bis zum Umkehrpunkt mit (w ≈ 1, Staub mit 1 + Q); dort fällt w (Preisach, θ̃)
auf 0, die örtliche Regel greift, und die Mode ξ geht mit λ_m = ω_orb/N* in die Operatorform über.
Vergleich: derselbe Bereich, in dem die Mode bis zuletzt wie Staub mitgeht (Standardkugel, Δ_vir ≈ 18π²).

0 Spektrum (CAMB wie stoerungen_rekombination: A_s = 2,1·10⁻⁹, n_s = 0,9649 [G]); Umkehr-z der Peaks 1–3σ, δ_lin = 1,062
1 Top-Hat im Kern-Hintergrund bis zur Umkehr (Staub 1 + Q, beide Fälle gleich); w(t) mit a₀(z), Übergabe vor der Umkehr I_pre
2 Nach der Umkehr: Schalen (Baryonen, Mode), gleicher Drehimpuls in allen Fällen (Perizentrum 0,05–0,3 des Startradius unter Newton)
     a = a_N,b·[1 + (1 − R_ξ)(ν(a_N,b; a₀(z)) − 1)] + R_ξ·a_N,ξ ,   dR_ξ/dt = −ω_orb/N*·R_ξ  (w = 0 nach der Umkehr, Preisach)
     ω_orb² = (4πG/3)ρ̄, ρ̄ = gezählte Dichte (1 + Q)ρ_b innerhalb des Halbmassenradius r_h der Baryonen
   Fälle: Staub (R_ξ ≡ 1), N* = 0,3 / 1 / 2,4, R_ξ = 0 ab Umkehr, R_ξ = 0 mit a₀ fest (ohne a₀(z))
3 Zweiter Weg: Energiesatz + Virialsatz (homologe gleichförmige Kugel, R_ξ beim Kollaps fest)
4 Was entscheidet: a_N,b(R_ta/2)/a₀ gegen 0,0292 (dort ν = 1 + Q); Grenzmasse je z
5 Von Hand
Blick: w Zustand des Bereichs (G am Ort: θ des Bereichs gegen θ_c), Regel ν(a_N,b) je Schale (L), a₀(z) aus Λ (I).
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, pc, kpc, Mpc, yr, m_p, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B_H2, OMEGA_C_H2,
                      Q_KOSM, H0_KMS, H0_si, E_hintergrund, a0_frueh_relativ, nu, gewicht_w, rho_m_mitbewegt)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)

H0 = H0_si()
FB = 1 / (1 + Q_KOSM)
rho_m0 = rho_m_mitbewegt()
Myr = 1e6 * yr


def E(a):
    return E_hintergrund(1 / a - 1)


# Zeit-Tabelle (Einheiten 1/H0), auch über heute hinaus
AG = np.geomspace(1e-5, 4.0, 6000)
TG = np.concatenate([[0.0], np.cumsum([quad(lambda y: 1 / (y * E(y)), AG[i], AG[i + 1])[0] for i in range(len(AG) - 1)])])
TG += quad(lambda y: 1 / (y * E(y)), 1e-12, AG[0])[0]


def t_of_a(a): return float(np.interp(math.log(a), np.log(AG), TG))
def a_of_t(t): return float(np.exp(np.interp(t, TG, np.log(AG))))
def z_of_t(t): return 1 / a_of_t(t) - 1
def a0z(a): return A0 * a0_frueh_relativ(1 / a - 1)


T_HEUTE = t_of_a(1.0)

# ---------------------------------------------------------------- 0 Spektrum
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=0.06, tau=0.0544)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=200)
PD = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False,
                                        kmax=200, var1='delta_tot', var2='delta_tot', zmax=200)
KS = np.exp(np.linspace(math.log(1e-4), math.log(200.0), 4000))
P0 = PD.P(0.0, KS)


def sigM(M_msun):
    R = (3 * M_msun * M_sun / (4 * math.pi * rho_m0)) ** (1 / 3) / Mpc
    x = KS * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(np.trapezoid(KS**3 * P0 / (2 * math.pi**2) * W**2, np.log(KS)))


def D_rel(z): return math.sqrt(PD.P(z, 0.05) / PD.P(0.0, 0.05))


ZG = np.concatenate([np.linspace(0, 10, 401), np.geomspace(10.05, 150, 300)])
DG = np.array([D_rel(z) for z in ZG])
def z_of_D(d): return float(np.interp(d, DG[::-1], ZG[::-1])) if d <= 1 else None


MASSEN = (1e8, 1e10, 1e12, 1e14)
PEAKS = (1, 2, 3)
P("0  Spektrum (CAMB, A_s = 2,1·10⁻⁹, n_s = 0,9649 [G]); M = gezählte Masse der Kugel (1 + Q)M_b, M_b = M/(1 + Q), Q = %.3f" % Q_KOSM)
P("   Umkehr bei δ_lin = 1,062 (Top-Hat);  M [M☉]   σ(M, 0)   z_ta für 1σ | 2σ | 3σ")
ZTA = {}
for M in MASSEN:
    s0 = sigM(M); row = []
    for n in PEAKS:
        z = z_of_D(1.062 / (n * s0))
        ZTA[(M, n)] = z if (z is not None and z > 0.05) else None
        row.append("  nie " if ZTA[(M, n)] is None else "%6.2f" % z)
    P("   %.0e   %6.3f   %s" % (M, s0, " | ".join(row)))

# ---------------------------------------------------------------- 1 Top-Hat bis zur Umkehr
AI = 1e-3


def rhs_th(t, y):
    a, x, v = y
    return [a * E(a), v, -OMEGA_M / (2 * x * x) - OMEGA_R * a**-4 * x + OMEGA_DE * x]


def tophat(di, tmax, bis_kollaps):
    x = AI * (1 - di / 3)
    v = x * E(AI) * (1 - di / (3 * (1 + di)))
    def ev_ta(t, y): return y[2]
    ev_ta.direction = -1; ev_ta.terminal = not bis_kollaps
    def ev_c(t, y): return y[1] - 1e-3 * AI
    ev_c.direction = -1; ev_c.terminal = True
    return solve_ivp(rhs_th, (t_of_a(AI), tmax), [AI, x, v], events=[ev_ta, ev_c], rtol=1e-10, atol=1e-14,
                     dense_output=True, max_step=tmax / 2000)


def tophat_ta(z_ta):
    t_target = t_of_a(1 / (1 + z_ta))
    lo, hi = 1e-4, 0.5
    for _ in range(55):
        mid = math.sqrt(lo * hi)
        s = tophat(mid, 3 * t_target, False)
        tt = s.t_events[0][0] if len(s.t_events[0]) else 1e9
        lo, hi = (mid, hi) if tt > t_target else (lo, mid)
    s = tophat(hi, 4 * t_target, True)
    t_ta = s.t_events[0][0]; a_ta, x_ta, _ = s.y_events[0][0]
    t_c = s.t_events[1][0] if len(s.t_events[1]) else None
    ts = np.linspace(t_of_a(AI), t_ta, 100001)
    Y = s.sol(ts)
    th = 3 * Y[2] / Y[1] * H0                                   # θ = 3Ṙ/R [1/s], fällt monoton bis zur Umkehr
    thmin = np.minimum.accumulate(th)                            # Preisach: θ̃ = kleinstes θ der Geschichte
    w = np.array([gewicht_w(thmin[i], a0z(Y[0][i])) for i in range(len(ts))])
    om = np.sqrt(OMEGA_M / (2 * Y[1]**3))                       # ω_orb in H0, gezählte Dichte (1 + Q)ρ_b
    I_pre = np.trapezoid((1 - w) * om, ts)
    i9 = np.argmax(w < 0.9); i1 = np.argmax(w < 0.1) if (w < 0.1).any() else len(ts) - 1
    return dict(t_ta=t_ta, x_ta=x_ta, a_ta=a_ta, t_c=t_c, I_pre=I_pre, dt91=(ts[i1] - ts[i9]) / H0, s=s)


# ---------------------------------------------------------------- 2 Schalen nach der Umkehr
NS = 200
GOLD = (math.sqrt(5) - 1) / 2
FPS = 0.05 + 0.25 * ((np.arange(NS) * GOLD) % 1.0)             # Perizentrum/Startradius unter Newton, gestreut (Phasenmischung)


def schalen(ph, modus, Nstar=None, a0_fest=False, T_nach=8.0, nstep=80000, FP=None):
    """Einheiten: Länge R_ta, Masse M = (1+Q)M_b, G = 1, Zeit T0 = √(R_ta³/GM).
       Gibt Kollapszeit (Median erstes Perizentrum der Baryonen ab Umkehr), Zeitmittel über [4; 8] t_ta nach der Umkehr:
       r_h (Median Baryonen), g(r_h), M_b(<r_h) = ½; R_ξ beim Kollaps und am Ende."""
    FPv = FPS if FP is None else FP[0] + (FP[1] - FP[0]) * ((np.arange(NS) * GOLD) % 1.0)
    T0 = ph["T0"]; A_u = ph["A_u"]
    m = (np.arange(NS) + 0.5) / NS
    xb = m ** (1 / 3); xm = xb * (1 + 1e-9)
    L2 = 2 * m * xb * FPv / (1 + FPv)                            # Kepler, eingeschlossene gezählte Masse m
    L2all = np.concatenate([L2, L2])
    x = np.concatenate([xb, xm]); v = np.zeros(2 * NS)
    typ_b = np.concatenate([np.ones(NS), np.zeros(NS)])
    R = {"std": 1.0, "null": 0.0}.get(modus, math.exp(-ph["I_pre"] / Nstar) if Nstar else 1.0)
    tta_u = ph["t_ta"] / H0 / T0                                 # t_ta in T0
    dt = T_nach * tta_u / nstep
    t_u = ph["t_ta"] / H0 / T0
    tH0 = 1 / (H0 * T0)

    def acc(x, R, t_u):
        a = a_of_t(t_u / tH0)
        a0 = A0 if a0_fest else a0z(a)
        o = np.argsort(x)
        xs = x[o]; tb = typ_b[o]
        wb = tb * FB / NS; wm = (1 - tb) * (1 - FB) / NS
        cb = np.cumsum(wb) - 0.5 * wb; cm = np.cumsum(wm) - 0.5 * wm
        gb = cb / xs**2; gm = cm / xs**2
        n_ = nu(np.maximum(gb, 1e-300) * A_u, a0)
        g = gb * (1 + (1 - R) * (n_ - 1)) + R * gm
        gg = np.empty_like(g); gg[o] = g
        bg = (OMEGA_DE - OMEGA_R * a**-4) * (H0 * T0)**2
        return -gg + L2all / x**3 + bg * x, gg
    ac, gg = acc(x, R, t_u)
    tperi = np.full(NS, np.nan); R_koll = None
    rec = []
    for k in range(nstep):
        vold = v[:NS].copy()
        v += 0.5 * dt * ac; x += dt * v; t_u += dt
        if modus == "regel":
            rh = np.median(x[:NS])
            om = math.sqrt(1.0 / (2 * rh**3))                     # (4πG/3)ρ̄ mit ρ̄ = gezählte Masse ½M in r_h
            R *= math.exp(-om / Nstar * dt)
        ac, gg = acc(x, R, t_u)
        v += 0.5 * dt * ac
        neu = np.isnan(tperi) & (vold < 0) & (v[:NS] >= 0)
        tperi[neu] = t_u
        if R_koll is None and np.sum(~np.isnan(tperi)) >= NS / 2:
            R_koll = R
        tn = (t_u - tta_u) / tta_u
        if tn >= 4.0 and k % 20 == 0:
            ob = np.argsort(x[:NS]); i = ob[NS // 2]
            rr = x[i]; gi = gg[i]
            mb_in = FB * (NS // 2 + 0.5) / NS                       # Baryonen innerhalb (gezählte Einheit M)
            rec.append((rr, gi * rr**2 / mb_in, gi * rr, 0.5 * math.pi * math.sqrt(rr / (2 * gi))))
    rec = np.array(rec)
    tcol = (np.nanmedian(tperi) - tta_u) / tta_u
    return dict(tcol=tcol, rh=rec[:, 0].mean(), rh_sd=rec[:, 0].std() / rec[:, 0].mean(), Mgel=rec[:, 1].mean(),
                vc2=rec[:, 2].mean(), tff=rec[:, 3].mean(), Rk=R if R_koll is None else R_koll, Rend=R)


# ---------------------------------------------------------------- 3 Energiesatz (zweiter Weg)
def g_hom(r, Rk, Mb, a0, Rq):
    """gleichförmige Kugel Radius Rk, Feld bei r ≤ Rk, Rest der Mode Rq folgt den Baryonen (a_N,ξ = Q·a_N,b)."""
    aN = G * Mb * (r / Rk)**3 / r**2
    return aN * (1 + (1 - Rq) * (float(nu(aN, a0)) - 1)) + Rq * Q_KOSM * aN


def virial_energie(R_ta, Mb, a0, Rq, n=60):
    ms = (np.arange(n) + 0.5) / n; xs = ms ** (1 / 3)

    def f(Rv):
        lhs = rhs = 0.0
        for xx, mm in zip(xs, ms):
            aN = lambda r: G * Mb * mm / r**2
            gf = lambda r: aN(r) * (1 + (1 - Rq) * (float(nu(aN(r), a0)) - 1)) + Rq * Q_KOSM * aN(r)
            u = np.geomspace(xx * Rv, xx * R_ta, 200)
            lhs += np.trapezoid([gf(r) for r in u], u)
            rhs += 0.5 * xx * Rv * gf(xx * Rv)
        return lhs - rhs
    return brentq(f, 0.02 * R_ta, 0.999 * R_ta, xtol=1e-6 * R_ta)


# ---------------------------------------------------------------- Rechnung je Fall
def rho_bar(z): return rho_m0 * (1 + z)**3


FAELLE = (("Staub", "std", None, False), ("N*=0,3", "regel", 0.3, False), ("N*=1", "regel", 1.0, False),
          ("N*=2,4", "regel", 2.4, False), ("R_ξ=0", "null", None, False), ("R_ξ=0,a₀fest", "null", None, True))
P()
P("1+2  Kollaps je Masse und Peak.  t_col = Median des ersten Perizentrums der Baryonenschalen ab Umkehr (in t_ta, Staub normiert);")
P("   r_h = Halbmassenradius der Baryonen (Zeitmittel 4–8 t_ta nach Umkehr), R_v = (R_ta/2)·r_h/r_h,Staub;  ρ_b über ½M_b in R_v/(0,5)^{1/3} → n_H;")
P("   M_gel/M_b = g(r_h)r_h²/(G·½M_b);  v_c = √(g r_h);  t_ff = (π/2)√(r_h/2g) (Einfall mit der wirkenden Beschleunigung);  z_SF = z bei t_col + t_ff.")
P("   R_v(E) = Energiesatz, gleichförmige Kugel, R_ξ beim Kollaps fest (zweiter Weg).")
ERG = {}
for M in MASSEN:
    for n in PEAKS:
        z_ta = ZTA[(M, n)]
        if z_ta is None:
            continue
        th = tophat_ta(z_ta)
        RL = (3 * M * M_sun / (4 * math.pi * rho_m0)) ** (1 / 3)
        Mt = M * M_sun; Mb = FB * Mt
        R_ta = th["x_ta"] * RL
        T0 = math.sqrt(R_ta**3 / (G * Mt)); A_u = G * Mt / R_ta**2
        th.update(T0=T0, A_u=A_u)
        a0_ta = a0z(th["a_ta"])
        tta_s = th["t_ta"] / H0
        xN = G * Mb / R_ta**2 / a0_ta
        xNv = 4 * xN
        z_c = z_of_t(th["t_c"]) if th["t_c"] else float("nan")
        Dvir = Mt / (4 / 3 * math.pi * (0.5 * R_ta)**3) / rho_bar(z_c)
        P()
        P("   M = %.0e M☉, %dσ: z_ta = %.2f, t_ta = %.1f Mio. J., R_ta = %.3f kpc, a₀(z)/a₀ = %.4f, a_N,b/a₀(z) bei R_ta %.4f, bei R_ta/2 %.4f"
          " (ν = %.2f; 1 + Q = %.2f)" % (M, n, z_ta, tta_s / Myr, R_ta / kpc, a0_ta / A0, xN, xNv, float(nu(xNv, 1.0)), 1 + Q_KOSM))
        P("      vor der Umkehr: w 0,9 → 0,1 in %.3g Mio. J., I_pre = ∫(1−w)ω dt = %.3g;  Staub-Probe: t_c/t_ta = %.4f, z_c = %.2f, Δ_vir = %.1f (18π² = %.1f)"
          % (th["dt91"] / Myr, th["I_pre"], th["t_c"] / th["t_ta"], z_c, Dvir, 18 * math.pi**2))
        rows = {}
        for key, modus, Ns, fest in FAELLE:
            r = schalen(th, modus, Ns, fest)
            rows[key] = r
        S = rows["Staub"]
        P("      Staub: r_h/R_ta = %.3f (Schwankung im Zeitmittel %.1f %%), gleichförmige Kugel mit R_ta/2: 0,397" % (S["rh"], 100 * S["rh_sd"]))
        P("      Fall           t_col  z_col   R_ξ(Koll) R_ξ(Ende)  R_v[kpc] R_v/Std  R_v(E)/R_ta  n_H[cm⁻³] ρ_b/Std  M_gel/M_b  v_c[km/s] v_c/Std  t_ff[Mio.J] t_ff/Std  z_SF")
        out = {}
        for key, modus, Ns, fest in FAELLE:
            r = rows[key]
            f = r["rh"] / S["rh"]
            Rv = 0.5 * R_ta * f
            Mgel = r["Mgel"]
            vc = math.sqrt(r["vc2"] * A_u * R_ta) / 1e3
            tff = r["tff"] * T0
            tcol_abs = th["t_ta"] * (1 + r["tcol"] / S["tcol"])            # Staub normiert: t_c = 2 t_ta
            z_col = z_of_t(tcol_abs)
            z_sf = z_of_t(tcol_abs + tff * H0)
            nH = 0.76 * Mb / 2 / (4 / 3 * math.pi * (Rv * 0.5**(1 / 3))**3) / m_p / 1e6
            Rq = 1.0 if modus == "std" else r["Rk"]
            RvE = virial_energie(R_ta, Mb, A0 if fest else a0_ta, Rq)
            out[key] = dict(tcol=r["tcol"] / S["tcol"], z_col=z_col, Rv=Rv, f=f, nH=nH, Mgel=Mgel, vc=vc, tff=tff, z_sf=z_sf,
                            RvE=RvE / R_ta, Rk=r["Rk"], Rend=r["Rend"])
        for key in out:
            o = out[key]; s = out["Staub"]
            P("      %-13s  %5.3f  %6.2f   %6.3f   %6.3f    %8.3f  %5.3f    %6.4f     %9.3g  %6.2f    %6.2f     %7.1f   %5.2f    %8.2f    %5.2f   %s"
              % (key, o["tcol"], o["z_col"], o["Rk"], o["Rend"], o["Rv"] / kpc, o["f"], o["RvE"], o["nH"], 1 / o["f"]**3, o["Mgel"],
                 o["vc"], o["vc"] / s["vc"], o["tff"] / Myr, o["tff"] / s["tff"], ("%6.2f" % o["z_sf"]) if o["z_sf"] > -0.99 else " —"))
        ERG[(M, n)] = dict(out=out, z_ta=z_ta, xNv=xNv, th=th, R_ta=R_ta)

# ---------------------------------------------------------------- Prüfungen der Numerik
P()
P("   Numerik: 10¹⁰ M☉ 2σ mit doppelter Schrittzahl und mit anderer Streuung des Perizentrums (Verhältnisse zum Staub)")
th = ERG[(1e10, 2)]["th"]
for lab, kw in (("160000 Schritte", dict(nstep=160000)), ("Peri 0,02–0,15", dict(FP=(0.02, 0.15))), ("Peri 0,1–0,5", dict(FP=(0.1, 0.5)))):
    s = schalen(th, "std", None, False, **kw); r = schalen(th, "regel", 1.0, False, **kw); q0 = schalen(th, "null", None, False, **kw)
    P("      %-16s Staub t_col = %.3f t_ta, r_h/R_ta = %.3f;  N*=1: t_col/Std = %.3f, R_v/Std = %.3f;  R_ξ=0: t_col/Std = %.3f, R_v/Std = %.3f"
      % (lab, s["tcol"], s["rh"], r["tcol"] / s["tcol"], r["rh"] / s["rh"], q0["tcol"] / s["tcol"], q0["rh"] / s["rh"]))

# ---------------------------------------------------------------- 4 Was entscheidet
P()
P("4  Entscheidende Größe: x_v = a_N,b(R_ta/2)/a₀.  ν(x_v) = 1 + Q bei x_v = (ln(1 + 1/Q))² = %.4f." % (math.log(1 + 1 / Q_KOSM)**2))
P("   x_v ∝ M^{1/3}(1+z_ta)² (ρ_ta = (9π²/16)ρ̄_m(z_ta)).  Grenzmasse M_x (x_v = 0,0292): darunter bindet die übergegangene Regel stärker als Staub, darüber schwächer.")
xs_c = math.log(1 + 1 / Q_KOSM)**2


def xv(M, z):
    Rt = (3 * M * M_sun / (4 * math.pi * 9 * math.pi**2 / 16 * rho_bar(z))) ** (1 / 3)
    return G * FB * M * M_sun / (0.5 * Rt)**2 / (A0 * a0_frueh_relativ(z))


P("   z_ta     M_x [M☉]     ν(x_v) für 10⁸ / 10¹⁰ / 10¹² / 10¹⁴ M☉")
for z in (0.5, 1, 2, 3, 5, 8, 12, 20, 30):
    Mx = 10 ** brentq(lambda lm: math.log(xv(10**lm, z) / xs_c), -5, 30)
    P("   %5.1f   %10.3g      %s" % (z, Mx, "  ".join("%6.2f" % float(nu(xv(M, z), 1.0)) for M in MASSEN)))

# ---------------------------------------------------------------- 5 Handprüfung
P()
P("5  Von Hand (unabhängig vom Code)")
P("   Δ_vir Staub (EdS): ρ_ta/ρ̄(t_ta) = 9π²/16 = 5,5517; R_v = R_ta/2 → ×8; t_c = 2t_ta, ρ̄ ∝ t⁻² → ×4: 5,5517·32 = 177,65 = 18π²")
P("     (Code mit Strahlung und Vakuumteil bei z_c ≥ 4: 178,3–179,1; bei z_c ≤ 1,3 wächst er durch den Vakuumteil auf 191–227).")
P("   Umkehrzeit 10⁸ M☉ 2σ, z_ta = 13,21: t = 2/(3H₀√Ω_m)·(1+z)^{−3/2}·(1 − (3/2)Ω_r(1+z)/Ω_m)")
P("     = 2/(3·2,1843·10⁻¹⁸·0,55982) = 5,4520·10¹⁷ s; /14,21^{1,5} = /53,57 → 1,0177·10¹⁶ s = 322,5 Mio. J.; ·(1 − 1,5·4,14·10⁻³) = 320,5 Mio. J. (Code 320,7)")
P("   Tiefer Grenzfall, 10⁸ M☉ 1σ, R_ξ = 0, mittlere Baryonenschale (½M_b, r₀ = 0,7937·6,736 kpc = 1,6496·10²⁰ m):")
P("     v∞⁴ = G·½M_b·a₀(z) = 1,0427·10²⁷·1,0977·10⁻¹⁰ = 1,1446·10¹⁷ → v∞ = 18,39 km/s  (Code v_c(r_h) = 18,6)")
P("     radialer Einfall t = √(π/2)·r₀/v∞ = 1,2533·1,6496·10²⁰/1,8394·10⁴ = 1,124·10¹⁶ s = 356 Mio. J. = 0,390 t_ta")
P("     (Code: Perizentrum 0,323·1,239 t_ta = 0,400 t_ta; Staub radial 1,000 t_ta, mit Drehimpuls 1,239 t_ta)")
P("   Energiesatz tiefer Grenzfall: jede Schale v_i²ln(R_ta/R_v) = ½v_i² → R_v/R_ta = e^{−1/2} = 0,6065 (Code R_v(E)/R_ta = 0,6055 bei 10⁸ 1σ)")
P("   Kreisgeschwindigkeit übergegangen (tief): v_c(r_h) = (G·½M_b·a₀)^{1/4}, 10¹⁴ M☉: (7,857·10¹²·1,3271·10²⁰·1,1·10⁻¹⁰)^{1/4} = (1,147·10²³)^{1/4} = 582 km/s")
P("     (Code 601–611: dort x_v = 0,02–0,05, nicht ganz tief)")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
