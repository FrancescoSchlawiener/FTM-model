# -*- coding: utf-8 -*-
"""
Frühes Wachstum: Kollaps, Dichte und Zeit bis zur Virialisierung früher Halos, wenn nach dem Umkehrpunkt
die örtliche Regel greift (w → 0, Mode ξ geht mit λ = (1 − w)ω_orb/N* in die Operatorform über).

0 Spektrum (CAMB wie stoerungen_rekombination, A_s = 2,1·10⁻⁹, n_s = 0,9649 [G]): σ(M), D(z), Umkehr-z der Peaks 1–4σ
1 Kugel bis zum Umkehrpunkt: Staub mit 1 + Q im Kern-Hintergrund (w ≈ 1), R_ta, t_ta, Übergabe vor der Umkehr
2 Nach der Umkehr, Schalen (Baryonen + Rest der Mode, gleicher Drehimpuls):
     a = a_N,b·[1 + (1 − R_ξ)(ν(a_N,b; a₀(z)) − 1)] + R_ξ·a_N,ξ ,  dR_ξ/dt = −ω_orb/N*·R_ξ (w = 0, Preisach)
  Standard: R_ξ ≡ 1 (Mode als Staub bis zum Ende). N* = 0,3 und 2,4; Grenzfall R_ξ = 0 am Umkehrpunkt.
  Unabhängig daneben: Energiesatz (Virialradius der gleichförmigen Kugel) und Einfallzeit der Außenschale.
3 Kern nach dem Kollaps: gelesene Masse / Baryonen, Dichte der Baryonen, Scheibe nach dem Abkühlen (gleiches j)
4 Anzahl kollabierter Halos (Press-Schechter mit wirksamem δ_c) und frühe Sternmasse gegen JWST [G]
5 Handprüfung
Blick: w Zustand des Bereichs (G am Ort), Regel ν(a_N,b) lokal (L), a₀(z) = A0·a0_frueh_relativ(z) innen (I).
Kern: A0, G, c, M_sun, kpc, Mpc, Gyr, yr, OMEGA_*, Q_KOSM, H0_si, E_hintergrund, a0_frueh_relativ, nu, gewicht_w,
      rho_m_mitbewegt, v_flach.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)

H0 = H0_si()
tH = 1 / H0                                   # s
OL = OMEGA_DE
FB = 1 / (1 + Q_KOSM)                         # Baryonenanteil
rho_m0 = rho_m_mitbewegt()
NSTAR = (0.3, 2.4)


def E(a):
    return E_hintergrund(1 / a - 1)


def t_of_a(a):
    return quad(lambda y: 1 / (y * E(y)), 1e-10, a, limit=300, epsabs=1e-14)[0]   # H0⁻¹


def a_of_t(t):
    return brentq(lambda a: t_of_a(a) - t, 1e-6, 2.0, xtol=1e-12)


def a0z(a):
    return A0 * a0_frueh_relativ(1 / a - 1)


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


def D_rel(z):
    return math.sqrt(PD.P(z, 0.05) / PD.P(0.0, 0.05))


ZG = np.concatenate([np.linspace(0, 10, 401), np.geomspace(10.05, 150, 300)])
DG = np.array([D_rel(z) for z in ZG])


def z_of_D(d):
    return float(np.interp(d, DG[::-1], ZG[::-1])) if d <= 1 else None


P("0  Spektrum (CAMB, A_s, n_s [G]), Massen = Gesamtmasse der Kugel (1 + Q)M_b, M_b = M/(1 + Q), Q = %.3f" % Q_KOSM)
P("   Umkehr: δ_lin(z_ta) = 1,062 (Top-Hat), Kollaps Standard δ_lin = 1,686")
MASSEN = (1e8, 1e10, 1e12, 1e13)
PEAKS = (1, 2, 3, 4)
P("   M [M☉]    σ(M,0)   z_ta für 1σ | 2σ | 3σ | 4σ")
ZTA = {}
for M in MASSEN:
    s0 = sigM(M); row = []
    for n in PEAKS:
        z = z_of_D(1.062 / (n * s0))
        ZTA[(M, n)] = z
        row.append("  nie " if z is None or z < 0.3 else "%6.2f" % z)
    P("   %.0e   %6.3f   %s" % (M, s0, " | ".join(row)))
P("   (z_ta < 0,3 wird nicht gerechnet)")

# ---------------------------------------------------------------- 1 Kugel bis zur Umkehr (Standard, 1 + Q)
AI = 1e-3


def rhs_th(t, y):
    a, x, v = y
    return [a * E(a), v, -OMEGA_M / (2 * x * x) - OMEGA_R * a**-4 * x + OL * x]


def tophat(di, tmax, bis_kollaps=True):
    x = AI * (1 - di / 3)
    v = x * E(AI) * (1 - di / (3 * (1 + di)))
    def ev_ta(t, y): return y[2]
    ev_ta.direction = -1; ev_ta.terminal = not bis_kollaps
    def ev_c(t, y): return y[1] - 1e-3 * AI
    ev_c.direction = -1; ev_c.terminal = True
    s = solve_ivp(rhs_th, (t_of_a(AI), tmax), [AI, x, v], events=[ev_ta, ev_c], rtol=1e-10, atol=1e-14,
                  dense_output=True, max_step=(tmax) / 2000)
    return s


def tophat_ta(z_ta):
    """δ_i so, dass die Kugel bei z_ta umkehrt; gibt t_ta, x_ta, a_ta, t_c (x → 0), I_pre(N*=1)."""
    t_target = t_of_a(1 / (1 + z_ta))
    lo, hi = 1e-4, 0.3
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        s = tophat(mid, 3 * t_target, bis_kollaps=False)
        tt = s.t_events[0][0] if len(s.t_events[0]) else 1e9
        if tt > t_target: lo = mid
        else: hi = mid
    s = tophat(hi, 4 * t_target)
    t_ta = s.t_events[0][0]; a_ta, x_ta, _ = s.y_events[0][0]
    t_c = s.t_events[1][0] if len(s.t_events[1]) else None
    # Übergabe vor der Umkehr: I = ∫(1 − w)ω dt, w = gewicht_w(θ̃; a₀(z)), θ = 3ẋ/x (monoton fallend bis zur Umkehr)
    ts = np.linspace(t_of_a(AI), t_ta, 200001)
    Y = s.sol(ts); th = 3 * Y[2] / Y[1] * H0
    w = np.array([gewicht_w(th[i], a0z(Y[0][i])) for i in range(0, len(ts))])
    om = np.sqrt(OMEGA_M / (2 * Y[1]**3))
    I_pre = np.trapezoid((1 - w) * om, ts)
    j9 = np.argmax(w < 0.9)
    return dict(t_ta=t_ta, x_ta=x_ta, a_ta=a_ta, t_c=t_c, I_pre=I_pre, t9=ts[j9] if w[-1] < 0.9 else t_ta)




# ---------------------------------------------------------------- 2 Nach der Umkehr
NS = 200
FP = 0.1          # Perizentrum/Umkehrradius jeder Schale unter Newton (Drehimpuls aus der Zeit vor der Umkehr, in allen Fällen gleich)


def schalen(th, RL, modus, Nstar=None, T_nach=1.6, nstep=16000):
    """Baryonen- und Modeschalen ab der Ruhe am Umkehrpunkt.
       modus 'std' (R_ξ ≡ 1), 'regel' (R_ξ(t_ta) = e^{−I_pre/N*}, dR/dt = −ω_orb/N*·R, w = 0 nach der Umkehr), 'null' (R_ξ = 0).
       Einheiten: Länge R_L, Zeit 1/H0, G·M ↔ Ω_M/2.  Gibt Zeit bis zum ersten Perizentrum von 50 % / 90 % der Baryonen und R_ξ dort."""
    x_ta, t_ta = th["x_ta"], th["t_ta"]
    m = (np.arange(NS) + 0.5) / NS
    xb = x_ta * m ** (1 / 3); xm = xb * (1 + 1e-9)
    L2 = 2 * (OMEGA_M / 2) * m * (FP * xb) * (1 - FP)
    vb = np.zeros(NS); vm = np.zeros(NS)
    R = 1.0 if modus == "std" else (0.0 if modus == "null" else math.exp(-th["I_pre"] / Nstar))
    dt = T_nach * t_ta / nstep
    t = t_ta; a = th["a_ta"]
    agN = H0**2 * RL
    kb = OMEGA_M / 2 * FB / NS; km = OMEGA_M / 2 * (1 - FB) / NS
    typ = np.concatenate([np.ones(NS), np.zeros(NS)])

    def acc(xb, xm, R, a):
        allx = np.concatenate([xb, xm])
        w8 = np.concatenate([np.full(NS, kb), np.full(NS, km * R)])
        o = np.argsort(allx)
        xs = allx[o]; ws = w8[o]; ts_ = typ[o]
        cb = np.cumsum(ws * ts_) - 0.5 * ws * ts_
        cm = np.cumsum(ws * (1 - ts_)) - 0.5 * ws * (1 - ts_)
        gb = cb / xs**2; gm = cm / xs**2
        n_ = np.where(gb > 0, nu(np.maximum(gb, 1e-300) * agN, a0z(a)), 1.0)
        g = gb * (1 + (1 - R) * (n_ - 1)) + R * gm
        gg = np.empty_like(g); gg[o] = g
        bg = -OMEGA_R * a**-4 + OL
        return -gg[:NS] + L2 / xb**3 + bg * xb, -gg[NS:] + L2 / xm**3 + bg * xm
    ab_, am_ = acc(xb, xm, R, a)
    tperi = np.full(NS, np.nan); Rperi = np.full(NS, np.nan)
    for k in range(nstep):
        vold = vb.copy()
        vb += 0.5 * dt * ab_; vm += 0.5 * dt * am_
        xb += dt * vb; xm += dt * vm
        t += dt
        a += dt * a * E(a)
        if modus == "regel":
            R *= math.exp(-math.sqrt(OMEGA_M / (2 * max(np.quantile(xb, 0.9), 0.5 * x_ta)**3)) / Nstar * dt)
        ab_, am_ = acc(xb, xm, R, a)
        vb += 0.5 * dt * ab_; vm += 0.5 * dt * am_
        neu = np.isnan(tperi) & (vold < 0) & (vb >= 0)
        tperi[neu] = t; Rperi[neu] = R
        if not np.isnan(tperi).any():
            break
    ok = ~np.isnan(tperi)
    t50 = np.nanquantile(np.where(ok, tperi, np.inf), 0.5)
    t90 = np.nanquantile(np.where(ok, tperi, np.inf), 0.9)
    return (t50 - t_ta) / t_ta, (t90 - t_ta) / t_ta, float(np.nanmedian(Rperi)) if ok.any() else R, np.nanmin(xb) / x_ta


# unabhängig von den Schalen: Außenschale als Punktmasse, Energiesatz der gleichförmigen Kugel
def g_mix(r, Mi, a0, R):
    aN = G * Mi / r**2
    return aN * (1 + (1 - R) * (float(nu(aN, a0)) - 1)) + R * Q_KOSM * aN


def aussenschale(R_ta, Mb, a0, t_ta_s, R0=1.0, Nstar=None):
    """Außenschale ab der Ruhe bei R_ta bis r = 10⁻³R_ta; dR/dt = −ω/N*·R, ω² = G(1+Q)M_b/max(r, R_ta/2)³. Zeit in t_ta."""
    def f(t, y):
        r, v, R = y
        r = max(r, 1e-6 * R_ta)
        om = math.sqrt(G * (1 + Q_KOSM) * Mb / max(r, 0.5 * R_ta)**3)
        return [v, -g_mix(r, Mb, a0, min(max(R, 0), 1)), (-om / Nstar * R) if Nstar else 0.0]
    def ev(t, y): return y[0] - 1e-3 * R_ta
    ev.terminal = True
    s = solve_ivp(f, (0, 10 * t_ta_s), [R_ta, 0.0, R0], events=ev, rtol=1e-9, atol=[1e-6 * R_ta, 1e-3, 1e-12])
    return s.t_events[0][0] / t_ta_s, s.y_events[0][0][2]


def virial_energie(R_ta, Mb, a0, R, n=100):
    """R_v aus E(R_ta, Ruhe) = K + E_pot(R_v), 2K = Σ m r g; gleichförmige Kugel, Schalen ohne Kreuzung, R_ξ fest."""
    xs = (np.arange(n) + 0.5) / n
    def dphi(Mi, r1, r2):
        u = np.geomspace(r2, r1, 300)
        return np.trapezoid([g_mix(v, Mi, a0, R) for v in u], u)
    def f(Rv):
        lhs = sum(dphi(Mb * x**3, x * R_ta, x * Rv) for x in xs) / n
        rhs = 0.5 * sum(x * Rv * g_mix(x * Rv, Mb * x**3, a0, R) for x in xs) / n
        return lhs - rhs
    return brentq(f, 0.05 * R_ta, 0.999 * R_ta, xtol=1e-5 * R_ta)


def rho_bar(z):
    return rho_m0 * (1 + z)**3


P()
P("1+2  Kollaps je Masse und Peak (Standard: Mode Staub bis zum Ende; Regel: Übergang N* = 0,3 / 2,4; Grenzfall R_ξ = 0 ab Umkehr)")
P("   t_col = Zeit ab Umkehr (in t_ta) bis zum ersten Perizentrum von 50 % der Baryonenschalen (Schalen mit Kreuzung, Perizentrum 0,1 r_ta unter Newton,\n   daher Std 1,155 = Kepler-Halbperiode mit r_p/r_a = 0,1);")
P("   Außenschale = Punktmasse M_b ab Ruhe (unabhängig); R_v aus dem Energiesatz (gleichförmige Kugel, R_ξ beim Kollaps fest);")
P("   Δ = gelesene Dichte in R_v / ρ̄_m(z_col), z_col aus der Außenschale (Std: t_c = 2t_ta);  ρ_b/Std = (R_ta/2 / R_v)³ bei gleichem R_ta;  M_gel/M_b bei R_v.")
ERG = {}
for M in MASSEN:
    for n in PEAKS:
        z_ta = ZTA[(M, n)]
        if z_ta is None or z_ta < 0.3:
            continue
        th = tophat_ta(z_ta)
        RL = (3 * M * M_sun / (4 * math.pi * rho_m0)) ** (1 / 3)
        Mb = FB * M * M_sun
        R_ta = th["x_ta"] * RL
        a0_ta = a0z(th["a_ta"])
        aNb = G * Mb / R_ta**2
        tta_s = th["t_ta"] * tH
        z_c = 1 / a_of_t(th["t_c"]) - 1
        P()
        P("   M = %.0e M☉, %dσ: z_ta = %.2f, t_ta = %.0f Mio. J., R_ta = %.2f kpc, a_N,b(R_ta)/a₀(z) = %.2e, ν = %.1f (1 + Q = %.2f), a₀(z)/a₀ = %.4f, I_pre = %.4f"
          % (M, n, z_ta, tta_s / yr / 1e6, R_ta / kpc, aNb / a0_ta, float(nu(aNb, a0_ta)), 1 + Q_KOSM, a0_ta / A0, th["I_pre"]))
        rows = {}
        for key, modus, Ns in (("Standard", "std", None), ("N* = 0,3", "regel", 0.3), ("N* = 2,4", "regel", 2.4), ("R_ξ = 0", "null", None)):
            t50, t90, Rc, xmin = schalen(th, RL, modus, Ns)
            R0 = 1.0 if modus == "std" else (0.0 if modus == "null" else math.exp(-th["I_pre"] / Ns))
            ta_out, Rout = aussenschale(R_ta, Mb, a0_ta, tta_s, R0, Ns)
            Rc_use = 1.0 if modus == "std" else (0.0 if modus == "null" else Rc)
            Rv = 0.5 * R_ta if modus == "std" else virial_energie(R_ta, Mb, a0_ta, Rc_use)
            if modus == "std":
                Rv_chk = virial_energie(R_ta, Mb, a0_ta, 1.0)
            gv = g_mix(Rv, Mb, a0_ta, Rc_use)
            Mread = gv * Rv**2 / G
            t_col_abs = th["t_ta"] * (1 + ta_out)
            z_col = 1 / a_of_t(t_col_abs) - 1
            Dl = Mread / (4 / 3 * math.pi * Rv**3) / rho_bar(z_col)
            vc = math.sqrt(gv * Rv) / 1e3
            rows[key] = dict(t50=t50, t90=t90, Rc=Rc_use, ta_out=ta_out, Rv=Rv, Mread=Mread, z_col=z_col, D=Dl, vc=vc,
                             rhob=(0.5 * R_ta / Rv)**3, Mr=Mread / Mb)
        S = rows["Standard"]
        P("      Std-Probe: t_c/t_ta (Top-Hat) = %.3f, Δ_vir (Masse in R_ta/2 / ρ̄(z_c)) = %.1f, 18π² = %.1f;  Energiesatz mit R_ξ = 1: R_v/R_ta = %.4f"
          % (th["t_c"] / th["t_ta"], M * M_sun / (4 / 3 * math.pi * (0.5 * R_ta)**3) / rho_bar(z_c), 18 * math.pi**2, Rv_chk / R_ta))
        P("      Fall        t_col  /Std   t_90  Außenschale  R_ξ(Koll.)  z_col   R_v/R_ta  ρ_b/Std  M_gel/M_b   Δ      v_c(R_v)[km/s]  v_c/Std")
        for key, r in rows.items():
            P("      %-10s  %5.3f  %5.3f  %5.3f   %5.3f       %5.3f     %6.2f   %5.3f    %5.2f    %6.2f   %7.1f     %6.1f        %5.2f"
              % (key, r["t50"], r["t50"] / S["t50"], r["t90"], r["ta_out"], r["Rc"], r["z_col"], r["Rv"] / R_ta, r["rhob"], r["Mr"], r["D"], r["vc"], r["vc"] / S["vc"]))
        ERG[(M, n)] = dict(th=th, rows=rows, z_ta=z_ta, R_ta=R_ta, Mb=Mb, a0=a0_ta, M=M)

# ---------------------------------------------------------------- 3 Scheibe nach dem Abkühlen
LAM_SPIN = 0.035       # Spinparameter [G] (Bullock 2001), vor der Umkehr gleich in beiden Bildern
C_NFW = 4.0            # Konzentration früher Halos [G]
P()
P("3  Kern nach dem Abkühlen: Exponentialscheibe M_b mit j = √2·λ·V_vir·R_vir (λ = %.3f [G], aus der Zeit vor der Umkehr, in beiden Bildern gleich)" % LAM_SPIN)
P("   j = 2·r_d·v_c(2,2 r_d);  Std: v_c² = G[M_NFW(<r) (c = %.0f [G], ohne Kontraktion) + M_b(<r)]/r;  Regel: v_c² = r·ν(a_N,b)a_N,b, R_ξ = 0" % C_NFW)
P("   Baryonenmasse der Scheibe kugelgenähert M_b(<r) = M_b[1 − (1 + r/r_d)e^{−r/r_d}];  t_orb = 2π·2,2r_d/v_c;  Σ₀ = M_b/(2πr_d²)")
P("   M, Peak      z_c   r_d Std/Regel [kpc]   Σ₀ Std/Regel [M☉/pc²]   v_c Std/Regel [km/s]   M_gel/M_b(<2,2r_d) Std/Regel   t_orb Std/Regel [Mio. J.]")


def Mexp(r, Mb, rd):
    return Mb * (1 - (1 + r / rd) * math.exp(-r / rd))


def Mnfw(r, M, Rv, cc):
    mf = lambda x: math.log(1 + x) - x / (1 + x)
    return M * (1 - FB) * mf(cc * min(r, Rv) / Rv) / mf(cc)


SCH = {}
for key, e in ERG.items():
    M, Mb, a0_ = e["M"], e["Mb"], e["a0"]
    Rv = 0.5 * e["R_ta"]
    Vv = math.sqrt(G * M * M_sun / Rv)
    j = math.sqrt(2) * LAM_SPIN * Vv * Rv
    vS = lambda r, rd: math.sqrt(G * (Mnfw(r, M * M_sun, Rv, C_NFW) + Mexp(r, Mb, rd)) / r)
    vF = lambda r, rd: math.sqrt(r * float(nu(G * Mexp(r, Mb, rd) / r**2, a0_)) * G * Mexp(r, Mb, rd) / r**2)
    rdS = brentq(lambda rd: 2 * rd * vS(2.2 * rd, rd) - j, 1e-4 * Rv, 5 * Rv)
    rdF = brentq(lambda rd: 2 * rd * vF(2.2 * rd, rd) - j, 1e-4 * Rv, 5 * Rv)
    vs, vf = vS(2.2 * rdS, rdS), vF(2.2 * rdF, rdF)
    MgS = vs**2 * 2.2 * rdS / G / Mexp(2.2 * rdS, Mb, rdS)
    MgF = vf**2 * 2.2 * rdF / G / Mexp(2.2 * rdF, Mb, rdF)
    toS = 2 * math.pi * 2.2 * rdS / vs / yr / 1e6; toF = 2 * math.pi * 2.2 * rdF / vf / yr / 1e6
    Sg = lambda rd: Mb / (2 * math.pi * rd**2) / M_sun * pc**2
    SCH[key] = dict(toS=toS, toF=toF)
    P("   %.0e %dσ  %6.2f   %6.3f / %6.3f        %8.0f / %8.0f     %6.1f / %6.1f          %5.2f / %5.2f               %6.2f / %6.2f"
      % (M, key[1], e["rows"]["Standard"]["z_col"], rdS / kpc, rdF / kpc, Sg(rdS), Sg(rdF), vs / 1e3, vf / 1e3, MgS, MgF, toS, toF))

# ---------------------------------------------------------------- 4 Anzahl und frühe Sternmasse
P()
P("4  Anzahl kollabierter Halos bei z_obs (Press-Schechter, gleiche Linearphase; nur die Zeit Umkehr → Kollaps ändert sich)")
P("   Kollaps bei t_obs verlangt Umkehr bei t_ta mit t_ta(1 + f) = t_obs; f = Außenschale/t_ta (Std f ≈ 1);  δ_c,eff = 1,062·D(z_obs)/D(z_ta)")
P("   R_ta aus ρ_ta = (9π²/16)ρ̄_m(z_ta);  n_Regel/n_Std = (δ_R/δ_S)·exp(−(δ_R² − δ_S²)/2σ²) bei festem M")


def t_of_z(z): return t_of_a(1 / (1 + z)) * tH


def z_of_t(t): return 1 / a_of_t(t / tH) - 1


def f_regel(M, z_ta, Nstar, R0=1.0):
    Mb = FB * M * M_sun
    R_ta = (3 * M * M_sun / (4 * math.pi * 9 * math.pi**2 / 16 * rho_bar(z_ta))) ** (1 / 3)
    a0_ = A0 * a0_frueh_relativ(z_ta)
    if Nstar == "std":
        return aussenschale(R_ta, Mb, a0_, t_of_z(z_ta), 1.0, None)[0]
    if Nstar == 0:
        return aussenschale(R_ta, Mb, a0_, t_of_z(z_ta), 0.0, None)[0]
    return aussenschale(R_ta, Mb, a0_, t_of_z(z_ta), R0, Nstar)[0]


def z_ta_fuer(M, z_obs, Nstar):
    tobs = t_of_z(z_obs)
    def h(zt):
        return t_of_z(zt) * (1 + f_regel(M, zt, Nstar)) - tobs
    return brentq(h, z_obs + 1e-3, 3 * z_obs + 30, xtol=1e-4)


DC = {}
P("   M [M☉]  z_obs   σ(M,z_obs)  δ_c Std   δ_c N*=2,4  N*=0,3  R_ξ=0   n/n_Std N*=2,4   N*=0,3      R_ξ=0")
for M in (1e10, 1e11, 1e12):
    s0 = sigM(M)
    for zo in (7.5, 9.0, 11.0, 13.0):
        sg = s0 * D_rel(zo)
        dc = {}
        for Ns in ("std", 2.4, 0.3, 0):
            zt = z_ta_fuer(M, zo, Ns)
            dc[Ns] = 1.062 * D_rel(zo) / D_rel(zt)
        DC[(M, zo)] = (dc, sg)
        rat = lambda d: (d / dc["std"]) * math.exp(-(d**2 - dc["std"]**2) / (2 * sg**2))
        P("   %.0e  %5.1f    %6.3f      %5.3f     %5.3f     %5.3f   %5.3f    %9.3g     %9.3g   %9.3g"
          % (M, zo, sg, dc["std"], dc[2.4], dc[0.3], dc[0], rat(dc[2.4]), rat(dc[0.3]), rat(dc[0])))

P()
P("   Obergrenze der Sternmasse (alle Baryonen der Halos M_h = M*/f_b zu Sternen, ε = 1): ρ*(>M*) = f_b·ρ̄_m·erfc(δ_c/√2σ(M_h))")
P("   (PS kumulativ, δ_c des Halos bei M_h; Labbé 2023 [G]: ρ*(>10¹⁰ M☉) ≈ 10⁵–10⁶ M☉/Mpc³ bei z ≈ 7–9, an oder über dieser Grenze im Standard)")
from scipy.special import erfc
rhoM = rho_m0 * Mpc**3 / M_sun
P("   z_obs  M*[M☉]   ρ*_max Std   N*=2,4       N*=0,3       R_ξ=0     [M☉/Mpc³]")
for zo in (7.5, 9.0):
    for Ms in (1e10, 10**10.5):
        Mh = Ms / FB
        sg = sigM(Mh) * D_rel(zo)
        vals = []
        for Ns in ("std", 2.4, 0.3, 0):
            zt = z_ta_fuer(Mh, zo, Ns)
            d = 1.062 * D_rel(zo) / D_rel(zt)
            vals.append(FB * rhoM * erfc(d / (math.sqrt(2) * sg)))
        P("   %4.1f   %.1e   %9.3g    %9.3g    %9.3g    %9.3g" % (zo, Ms, *vals))

P()
P("   Sternentstehung (Annahme): Ṁ* = ε_ff·M_gas/t_ff, t_ff ∝ t_orb der Scheibe; Verhältnis der Raten Regel/Std = t_orb,Std/t_orb,Regel:")
for key, sv in SCH.items():
    P("   %.0e %dσ: %.2f" % (key[0], key[1], sv["toS"] / sv["toF"]))

# ---------------------------------------------------------------- 5 Handprüfung
P()
P("5  Von Hand (unabhängig vom Code oben)")
P("   Δ_vir Std: ρ_ta/ρ̄(t_ta) = 9π²/16 = 5,551; R_vir = R_ta/2 → ×8; t_c = 2t_ta, ρ̄ ∝ t⁻² → ×4: 5,551·32 = 177,65 = 18π²")
P("   Umkehr-z 2σ, 10⁸ M☉: D(z)/D(0) = 1,062/(2·5,898) = 0,0900; D(0)/a = g₀ = 2,5Ω_m/(Ω_m^{4/7} − Ω_Λ + (1 + Ω_m/2)(1 + Ω_Λ/70))")
P("     = 2,5·0,3134/(0,5153 − 0,6865 + 1,1567·1,0098) = 0,7835/0,9968 = 0,786; 1 + z = 1/(0,0900·0,786) = 14,1 → z_ta ≈ 13,1 (Code 13,21)")
P("   Tiefer Grenzfall, Außenschale 10⁸ M☉ 1σ: a = v∞²/r, v∞⁴ = GM_b·a₀(z) = 1,571·10⁷·1,327·10²⁰·1,0977·10⁻¹⁰ = 2,289·10¹⁷ → v∞ = 21,9 km/s;")
P("     t = ∫dr/√(2v∞² ln(R/r)) = √(π/2)·R/v∞ = 1,2533·2,080·10²⁰/2,187·10⁴ = 1,19·10¹⁶ s = 378 Mio. J. = 0,413 t_ta (Code mit vollem ν: 0,409)")
P("   Energiesatz tiefer Grenzfall: jede Schale v_i² ln(R_ta/R_v) = ½v_i² → R_v/R_ta = e^{−1/2} = 0,607 (Code 0,606)")
P("   Kepler mit r_p/r_a = 0,1: Halbperiode/radialer Einfall = ((1 + 0,1)/2)^{3/2}/(1/2)^{3/2} = 1,1^{1,5} = 1,154 (Schalen Std 1,155)")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
