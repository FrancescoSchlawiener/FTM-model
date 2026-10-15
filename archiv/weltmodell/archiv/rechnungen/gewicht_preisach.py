# -*- coding: utf-8 -*-
"""
Gewicht w der Schließungsregel als Preisach-Zustand.
  jetzige Form:  w = θ₊²/(θ₊² + θ_c²),  θ₊ = max(θ, 0)                       (Kern: gewicht_w)
  Preisach-Form: w = θ̃²/(θ̃² + θ_c²),   θ̃ = max(min_{t'≤t} θ(t'), 0)        (tiefster Kollaps der Geschichte)
1 Hintergrund, 2 Top-Hat-Kollaps (z_c = 3,4 / 2 / 1 / 0,5) und Wiederausdehnung,
3 frisch kollabierte Galaxien, 4 Haufenstoß, 5 Handprüfung.
Kern: A0, OMEGA_*, Q_KOSM, q_KOSM, E_hintergrund, gewicht_w, nu, hernquist_masse, uebergabezeit.
Nicht im Kern (hier nach der Rechnung zur Schließungsregel nachgebaut): Top-Hat, Hernquist-Galaxie
M_b = 6·10¹⁰, a = 1,5 kpc·(M_b/10¹⁰)^0,3, g/g_RAR mit Rest-Gedächtnis R_ξ, N* aus der RAR-Enge.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

OUT = []
def P(s=""):
    OUT.append(s); print(s)

H0 = H0_si()
tH = 1 / H0 / Gyr                     # 1/H0 in Gyr
thc = (A0 / c) / H0                   # θ_c in H0
OL = OMEGA_DE


def E(a):
    return E_hintergrund(1 / a - 1)


def w_jetzt(th):                      # th in H0-Einheiten, über den Kern
    return gewicht_w(th * H0)


def w_preisach(th_hist):              # Folge θ(t'), t' ≤ t
    tmin = np.minimum.accumulate(np.asarray(th_hist))
    return np.array([gewicht_w(max(x, 0.0) * H0) for x in tmin])


def t_of_a(a):
    return quad(lambda y: 1 / (y * E(y)), 1e-9, a, limit=200)[0]


t0 = t_of_a(1.0)

# ---------------------------------------------------------------- 1 Hintergrund
P("1  Hintergrund (A0 = %.2e, H0 = %.1f, θ_c/H0 = %.5f, θ_c = %.4f /Gyr)" % (A0, H0_KMS, thc, thc / tH))
zs = np.concatenate([np.geomspace(1100, 1e-3, 4000), [0.0]])
th_bg = np.array([3 * E_hintergrund(z) for z in zs])
dmono = np.max(np.diff(th_bg))
wj = np.array([w_jetzt(x) for x in th_bg]); wp = w_preisach(th_bg)
P("   θ_bg = 3H fällt monoton (größter Schritt %.1e ≤ 0): min der Geschichte = heutiger Wert" % dmono)
P("   1 − w heute: jetzige Form %.4e, Preisach %.4e, Differenz %.1e" % (1 - wj[-1], 1 - wp[-1], abs(wj[-1] - wp[-1])))
P("   1 − w bei z = 1090: %.2e (beide); mit a0(1090) = %.3f a0 (Kern a0_frueh_relativ): %.2e"
  % (1 - w_jetzt(3 * E_hintergrund(1090)), a0_frueh_relativ(1090), 1 - gewicht_w(3 * E_hintergrund(1090) * H0, A0 * a0_frueh_relativ(1090))))
P("   mit a0 = cH0/2π: 1 − w heute %.4e" % (1 - gewicht_w(3 * H0, a0_modell())))

# ---------------------------------------------------------------- 2 Top-Hat
ai = 1e-3


def rhs(t, y):
    a, x, v = y
    return [a * E(a), v, -OMEGA_M / (2 * x * x) - OMEGA_R * a**-4 * x + OL * x]


def ev_ta(t, y): return y[2]
ev_ta.terminal = False; ev_ta.direction = -1


def tophat(di, tmax):
    x = ai * (1 - di / 3)
    v = x * E(ai) * (1 - di / (3 * (1 + di)))
    # Kollaps: x < 1 % von x_ta; x_ta aus dem ersten Lauf
    s1 = solve_ivp(rhs, (t_of_a(ai), tmax), [ai, x, v], events=ev_ta, rtol=1e-10, atol=1e-13, max_step=1e-3)
    if len(s1.t_events[0]) == 0:
        return None
    t_ta = s1.t_events[0][0]; y_ta = s1.y_events[0][0]
    x_ta = y_ta[1]
    def ev_col(t, y): return y[1] - 0.01 * x_ta
    ev_col.terminal = True; ev_col.direction = -1
    s2 = solve_ivp(rhs, (t_of_a(ai), tmax), [ai, x, v], events=ev_col, rtol=1e-10, atol=1e-14,
                   max_step=1e-3, dense_output=True)
    t_c = s2.t_events[0][0] if len(s2.t_events[0]) else None
    return dict(di=di, t_ta=t_ta, x_ta=x_ta, t_c=t_c, sol=s2)


def collapse_at(zc):
    tc = t_of_a(1 / (1 + zc))
    lo, hi = 1e-4, 0.05
    for _ in range(45):
        mid = math.sqrt(lo * hi)
        r = tophat(mid, 1.6 * tc + 0.05)
        if r is None or r["t_c"] is None or r["t_c"] > tc:
            lo = mid
        else:
            hi = mid
    return tophat(hi, 1.6 * tc + 0.05)


def verlauf(r, n=400000):
    """θ(t), w_jetzt, w_Preisach, I_pre = ∫(1−w)ω dt bis t_c (Zone endet bei x_ta/2)."""
    ts = np.concatenate([np.linspace(t_of_a(ai), r["t_ta"] - 2e-3, n // 2),
                         np.linspace(r["t_ta"] - 2e-3, r["t_ta"] + 2e-3, n // 4),
                         np.linspace(r["t_ta"] + 2e-3, r["t_c"], n // 4)])
    ts = np.unique(ts)
    Y = r["sol"].sol(ts)
    x, v = Y[1], Y[2]
    th = 3 * v / x
    wj = np.array([w_jetzt(u) for u in th])
    wp = w_preisach(th)
    xe = np.where(ts < r["t_ta"], x, np.maximum(x, 0.5 * r["x_ta"]))
    om = np.sqrt(OMEGA_M / (2 * xe**3))
    Ij = np.trapezoid((1 - wj) * om, ts)
    Ip = np.trapezoid((1 - wp) * om, ts)
    return ts, x, th, wj, wp, Ij, Ip


def breite(ts, w):
    i9 = np.argmax(w < 0.9); i1 = np.argmax(w < 0.1)
    return ts[i9], ts[i1]


P()
P("2  Top-Hat im Kern-Hintergrund (Materie, Strahlung, Vakuumteil); θ = 3Ṙ/R")
Mb = 6e10 * M_sun
ab = 1.5 * kpc * (6e10 / 1e10) ** 0.3
Mtot = (1 + Q_KOSM) * Mb
rho_m0 = rho_m_mitbewegt()
RL = (3 * Mtot / (4 * math.pi * rho_m0)) ** (1 / 3)
res = {}
P("   z_c   z_ta   t_ta[Gyr]  r_vir[kpc]  ω_vir[1/Gyr]  w 0,9→0,1 [Myr]  vor Umkehr [Myr]  max|w_j−w_P| bis t_c  Hand Δt [Myr]")
for zc in (3.4, 2.0, 1.0, 0.5):
    r = collapse_at(zc)
    ts, x, th, wj, wp, Ij, Ip = verlauf(r)
    t9, t1 = breite(ts, wj)
    # Hand: θ ≈ θ̇_ta (t − t_ta), θ̇_ta = 3R̈/R;  w = 0,9 bei θ = 3θ_c, w = 0,1 bei θ = θ_c/3
    a_ta = r["sol"].sol(r["t_ta"])[0]
    thdot = 3 * (-OMEGA_M / (2 * r["x_ta"]**3) - OMEGA_R * a_ta**-4 + OL)
    dt_hand = (3 - 1 / 3) * thc / abs(thdot)
    ov = math.sqrt(OMEGA_M / (2 * (r["x_ta"] / 2)**3))
    res[zc] = dict(r=r, Ij=Ij, Ip=Ip, ov=ov, rv=r["x_ta"] / 2 * RL)
    P("   %4.1f  %5.2f   %7.3f    %7.1f      %6.3f         %6.2f           %6.2f           %.1e          %6.2f"
      % (zc, 1 / a_ta - 1, r["t_ta"] * tH, r["x_ta"] / 2 * RL / kpc, ov / tH, (t1 - t9) * tH * 1e3,
         (r["t_ta"] - t9) * tH * 1e3, np.max(np.abs(wj - wp)), dt_hand * tH * 1e3))
    if zc == 3.4:
        P("     Verlauf z_c = 3,4:  t[Gyr]   R/R_ta   θ/H0      w_jetzt   w_Preisach")
        for f in (0.2, 0.8, 0.99, 0.999, 1.0, 1.001, 1.5, 1.99):
            j = int(np.argmin(abs(ts - f * r["t_ta"])))
            P("                        %6.3f   %6.3f   %9.3f   %.4f    %.4f" % (ts[j] * tH, x[j] / r["x_ta"], th[j], wj[j], wp[j]))
P("   θ fällt bis zum Umkehrpunkt monoton, danach θ < 0 und virialisiert θ = 0: beide Formen gleich, I_pre gleich.")
P("   Handwert: Δt = (3 − 1/3)·θ_c/|θ̇_ta|, θ̇_ta = 3R̈/R am Umkehrpunkt.")

# Wiederausdehnung nach dem Kollaps
P()
P("   Wiederausdehnung eines virialisierten Bereichs (z_c = 1, Ereignis 0,5 Gyr nach t_c):")
P("   θ(t) = θ_e·e^{−(t−t_e)/τ}, Radiuszuwachs f = e^{∫θdt/3} − 1;  Verlust an Übergabe ΔI = ∫ w ω_vir dt")
ov1 = res[1.0]["ov"] / tH                       # 1/Gyr
thc_G = thc / tH
for f, tau in ((0.001, 0.1), (0.01, 0.1), (0.1, 0.1), (0.1, 1.0), (0.3, 1.0)):
    the = 3 * math.log(1 + f) / tau
    tt = np.linspace(0, 10 * tau, 20001)
    th = the * np.exp(-tt / tau)
    wj = th**2 / (th**2 + thc_G**2)
    dI = np.trapezoid(wj * ov1, tt)
    P("     f = %5.3f, τ = %.1f Gyr: θ_e = %.3f/Gyr = %.0f θ_c;  w_jetzt max %.5f, w_jetzt > 0,5 für %.3f Gyr;  ΔI = %.3f;  w_Preisach = 0"
      % (f, tau, the, the / thc_G, wj.max(), tt[wj > 0.5][-1] if np.any(wj > 0.5) else 0.0, dI))
P("   Atmende Schwingung um das Gleichgewicht R = R_v(1 + ε sin ω_vir t):  mittleres w der jetzigen Form (Preisach: 0)")
for eps in (1e-3, 1e-2, 0.1):
    ph = np.linspace(0, 2 * math.pi, 200001)
    th = 3 * eps * ov1 * np.cos(ph) / (1 + eps * np.sin(ph))
    wj = np.array([gewicht_w(u / Gyr) for u in th])
    P("     ε = %.0e: θ-Amplitude %.2f θ_c, ⟨w⟩ = %.3f, Übergaberate ⟨1−w⟩ω = %.3f der Rate bei w = 0" % (eps, 3 * eps * ov1 / thc_G, wj.mean(), 1 - wj.mean()))

# ---------------------------------------------------------------- 3 frisch kollabierte Galaxien
P()
P("3  Frisch kollabierte Galaxien (M_b = 6·10¹⁰ M☉, Hernquist a = %.2f kpc, Top-Hat (1+Q)M_b, Kern-Hintergrund)" % (ab / kpc))


def ratio_rar(rk, W, rv):
    rr = rk * kpc
    mb = hernquist_masse(Mb, ab, rr)
    gN = G * mb / rr**2
    n_ = float(nu(gN))
    g = gN + (1 - W) * (n_ - 1) * gN + W * G * Q_KOSM * Mb * min(rr / rv, 1.0) / rr**2
    return g / (n_ * gN)


RQ = brentq(lambda R: (float(nu(G * hernquist_masse(Mb, ab, R) / R**2)) - 1) * hernquist_masse(Mb, ab, R) - Q_KOSM * Mb, 1 * kpc, 500 * kpc)
zz = [3.4, 2.0, 1.0, 0.5]
lr = [math.log(res[z]["rv"] / RQ) for z in zz]
P("   R_Q = %.2f kpc; r_vir/R_Q: %s" % (RQ / kpc, ", ".join("z_c=%g: %.3f" % (z, res[z]["rv"] / RQ) for z in zz)))
zcQ = brentq(lambda zc: math.log(collapse_at(zc)["x_ta"] / 2 * RL / RQ), 3.0, 4.0, xtol=2e-3)
rQ = collapse_at(zcQ)
_, _, _, _, _, IjQ, IpQ = verlauf(rQ)
res["Q"] = dict(r=rQ, Ij=IjQ, Ip=IpQ, ov=math.sqrt(OMEGA_M / (2 * (rQ["x_ta"] / 2)**3)), rv=rQ["x_ta"] / 2 * RL)
P("   r_vir = R_Q bei z_c = %.3f" % zcQ)


def I_at(key, t, form):
    d = res[key]
    return (d["Ij"] if form == "j" else d["Ip"]) + d["ov"] * max(t - d["r"]["t_c"], 0.0)


cf05 = max(abs(ratio_rar(rk, 1.0, res[0.5]["rv"]) - 1) for rk in (10, 20, 30, 50))
NRAR = I_at(0.5, t0, "j") / math.log(cf05 / 0.05)
NRARp = I_at(0.5, t0, "p") / math.log(cf05 / 0.05)
P("   N* aus der RAR-Enge (z_c = 0,5, |g/g_RAR − 1| ≤ 5 %% bei 10–50 kpc heute): jetzige Form %.3f, Preisach %.3f" % (NRAR, NRARp))
P("   Kern uebergabezeit: 1/λ_m (z_c = 0,5, Δ aus Top-Hat, N* = %.2f) = %.2f Gyr; hier N*/ω_vir = %.2f Gyr"
  % (NRAR, uebergabezeit(OMEGA_M * (res[0.5]["r"]["x_ta"] / 2)**-3, 0.0, (NRAR, NRAR))[1] / Gyr,
     NRAR / (res[0.5]["ov"] / tH)))
P("   z_c     z_obs   R_ξ (jetzt)   R_ξ (Preisach)   g/g_RAR − 1 bei 10 / 20 / 50 kpc (Preisach)")
for key, zo in (("Q", 2.0), (3.4, 2.0), (2.0, 1.0), (1.0, 0.5), (0.5, 0.0)):
    tobs = t_of_a(1 / (1 + zo))
    Wj = math.exp(-I_at(key, tobs, "j") / NRAR)
    Wp = math.exp(-I_at(key, tobs, "p") / NRARp)
    rv = res[key]["rv"]
    zcs = zcQ if key == "Q" else key
    P("   %5.2f   %4.1f   %.4e    %.4e       %s" % (zcs, zo, Wj, Wp, "  ".join("%+.4f" % (ratio_rar(rk, Wp, rv) - 1) for rk in (10, 20, 50))))

# ---------------------------------------------------------------- 4 Haufenstoß
P()
P("4  Haufen 10¹⁵ M☉ (z = 0,3; Gas β = 2/3, r_c = 150 kpc, 10¹⁴ in 1 Mpc; BCG 1,5·10¹²; Galaxien 10¹³; NFW c = 4)")
zl = 0.3
rho_cz = 3 * (H0 * E_hintergrund(zl))**2 / (8 * math.pi * G)
rcut = 3000 * kpc
Mbeta = lambda r, rc: min(r, rcut) / rc - math.atan(min(r, rcut) / rc)
Mb_cl = lambda r: (1e14 * M_sun * Mbeta(r, 150 * kpc) / Mbeta(1000 * kpc, 150 * kpc)
                   + hernquist_masse(1.5e12 * M_sun, 20 * kpc, r) + hernquist_masse(1e13 * M_sun, 300 * kpc, r))
r200 = (3 * 1e15 * M_sun / (4 * math.pi * 200 * rho_cz)) ** (1 / 3); rs = r200 / 4
mf = lambda x: math.log(1 + x) - x / (1 + x)
Mref = lambda r: 1e15 * M_sun * mf(r / rs) / mf(4)
for rk in (100, 250, 1000, 1500):
    rr = rk * kpc; mb = Mb_cl(rr); n_ = float(nu(G * mb / rr**2))
    Wr = 1 - (n_ - 1) * mb / (q_KOSM * Mref(rr))
    P("   %5d kpc: ν = %.3f, stoßfrei R_ξ = %.3f, Operator %.3f" % (rk, n_, Wr, 1 - Wr))
P("   R_ξ = 1 − (ν−1)M_b/(q M_ref) enthält w nicht: vor dem Stoß (θ = 0) ist w in beiden Formen 0 → 80/20 unverändert.")
P("   Nach dem Stoß: auseinanderlaufende Bereiche, θ ≈ 3v/d (v, d angenommen):")
for v_ in (1000, 3000, 4500):
    for d_ in (300, 1000):
        th = 3 * v_ * 1e3 / (d_ * kpc)
        P("     v = %4d km/s, d = %4d kpc: θ = %.1f/Gyr = %.0f θ_c → w_jetzt = 1 − %.1e, w_Preisach = 0"
          % (v_, d_, th * Gyr, th / (A0 / c), 1 - gewicht_w(th)))
P("   Gasanteil der Übergabe nach dem Stoß 1 − (1 − e^{−x})/x, x = λ_m t: gilt nur mit Preisach (jetzige Form: λ = (1−w)ω/N* ≈ 0 solange θ > 0)")

# ---------------------------------------------------------------- 5 Handprüfung
P()
P("5  Handprüfung")
r_ = thc / 3
P("   1 − w_bg heute = r²/(1 + r²), r = (a0/c)/(3H0) = %.6f → %.4e" % (r_, r_**2 / (1 + r_**2)))
P("   von Hand: a0/c = 1,1e-10/2,998e8 = 3,669e-19 /s; 3H0 = 3·67,4/3,0857e19 = 6,553e-18 /s; r = 0,05599; r² = 3,135e-3; 1 − w = 3,125e-3")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
