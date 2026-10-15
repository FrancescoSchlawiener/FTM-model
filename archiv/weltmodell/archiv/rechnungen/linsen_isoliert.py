# -*- coding: utf-8 -*-
"""Linse um isolierte Galaxien bis 2,6 Mpc: Regel mit äußerem Feld und Rest der Mode ξ gegen KiDS-1000.

Messung (eingebettet, Quelle unten): Brouwer et al. 2021 (A&A 650, A113), isolierte Linsen, vier Massenbins,
g_obs = 4G·ΔΣ (Umrechnung wie in der Messung), g_bar = G·M*(1 + f_kalt)/R² (Punktmasse, ohne heißes Gas);
Mistele et al. 2024 (ApJL 969, L3), v_circ aus Entprojektion, Tabelle 1 (34–115 kpc).
Die Zahlen stammen aus dem öffentlichen Repositorium github.com/robertmonjo/lensing_RAR_KiDS
(Dateien Fig-9_RAR-KiDS-isolated_Massbin-*.txt, external_reference/mistele_gg_wl.txt); arxiv.org war nicht erreichbar,
gegen die Originalarbeit also nicht nachgeprüft.

Modell (nur ftm_kern):
  Regel      g = ν(|g_N|)·g_N, äußeres Feld g_e (Newton-Einheiten wie in 03): senkrecht |g_N| = √(g_i² + g_e²),
             gleichgerichtet g_i + g_e, Kugelmittel = Monopol der Radialkomponente ½∫ν(|g_i e_r + g_e e_z|)(g_i + g_e μ)dμ.
             Der Monopol ist genau das, was der gestapelte Tangentialschub sieht: Der wirbelbehaftete Teil hat keinen Fluss.
  Mode ξ     Rest ξ_R = exp(−∫(1 − w)ω dt / N*) je Lagrange-Schale (w aus ftm_kern.gewicht_w, ω² = GM/R³,
             nach Erreichen von r_ta/2 fest, wie in rest_mode_isolierte_galaxie.py); Schale so gewählt, dass sie zur
             Beobachtungszeit bei r liegt (gebunden, im Einfall oder noch expandierend).
             lokal (Ansatz aus 03):   g = g_N·[ν − ξ_R((ν − 1) − Q·min(r/r_vir, 1))]
             kumulativ (Schalen):     dM_L/dr = (1 − ξ_R)·dM_X/dr + ξ_R·Q·M_b/r_vir (r < r_vir),  M_X = (ν − 1)·M_b gelesen
  Linse      ΔΣ(R) = M_2D(<R)/(πR²) − Σ(R) aus M_L(<r) = g·r²/G, verglichen als 4G·ΔΣ.
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

G, c, Ms, kpc, Gyr = K.G, K.c, K.M_sun, K.kpc if hasattr(K, "kpc") else 1e3 * K.pc, K.Gyr
Mpc = 1e3 * kpc
a0, Q = K.A0, K.Q_KOSM
H0, Om = K.H0_si(), K.OMEGA_M
OL = 1 - Om

# ------------------------------------------------------------------ Messung
MSTAR = [1.5e10, 3.2e10, 4.6e10, 8.9e10]
f_kalt = lambda m: 10 ** (-0.69 * math.log10(m) + 6.63)       # wie in der Messung (Boselli 2014)
MB = [m * (1 + f_kalt(m)) for m in MSTAR]
R_MPC = np.array([3.5390e-02, 4.8108e-02, 6.5396e-02, 8.8896e-02, 1.2084e-01, 1.6427e-01, 2.2330e-01, 3.0354e-01,
                  4.1262e-01, 5.6090e-01, 7.6247e-01, 1.0365e+00, 1.4089e+00, 1.9152e+00, 2.6035e+00])
KIDS = [np.array(x).reshape(-1, 3) for x in (
 [8.9678e-12,6.1374e-12,1.9384e-11, 8.7469e-12,7.3406e-12,1.3834e-11, 8.0130e-12,6.7410e-12,1.0456e-11, 6.2018e-12,5.1260e-12,7.8323e-12,
  4.3587e-12,3.6570e-12,5.3527e-12, 2.9690e-12,2.5314e-12,3.5749e-12, 2.0745e-12,1.7913e-12,2.4488e-12, 1.5605e-12,1.3583e-12,1.8130e-12,
  1.2783e-12,1.1054e-12,1.4809e-12, 1.0607e-12,9.1520e-13,1.2117e-12, 7.8897e-13,6.7724e-13,9.0410e-13, 5.2832e-13,4.5092e-13,6.0996e-13,
  3.3855e-13,2.8545e-13,3.9559e-13, 2.1428e-13,1.7641e-13,2.5873e-13, 1.3393e-13,1.0657e-13,1.7144e-13],
 [1.5240e-11,9.8535e-12,3.3211e-11, 1.4124e-11,1.1616e-11,2.3146e-11, 1.2871e-11,1.1021e-11,1.7090e-11, 1.0679e-11,9.1144e-12,1.3259e-11,
  8.6016e-12,7.5652e-12,1.0177e-11, 7.1304e-12,6.3996e-12,8.1249e-12, 5.7500e-12,5.1995e-12,6.4133e-12, 4.3201e-12,3.9354e-12,4.7609e-12,
  3.0574e-12,2.7931e-12,3.3592e-12, 2.0942e-12,1.9101e-12,2.2974e-12, 1.3954e-12,1.2591e-12,1.5433e-12, 8.8727e-13,7.8691e-13,9.9376e-13,
  5.4361e-13,4.7481e-13,6.1821e-13, 3.2944e-13,2.8243e-13,3.8415e-13, 1.9842e-13,1.6559e-13,2.4103e-13],
 [1.6935e-11,1.1562e-11,3.5642e-11, 1.6440e-11,1.3986e-11,2.5662e-11, 1.5090e-11,1.3193e-11,1.9514e-11, 1.2426e-11,1.0795e-11,1.5236e-11,
  9.8112e-12,8.6715e-12,1.1521e-11, 7.5314e-12,6.7595e-12,8.5855e-12, 5.9140e-12,5.3595e-12,6.6054e-12, 4.6728e-12,4.2629e-12,5.1444e-12,
  3.5961e-12,3.3011e-12,3.9225e-12, 2.7037e-12,2.4906e-12,2.9365e-12, 1.9644e-12,1.8099e-12,2.1292e-12, 1.3641e-12,1.2494e-12,1.4839e-12,
  9.0829e-13,8.1959e-13,1.0011e-12, 5.8463e-13,5.1449e-13,6.6209e-13, 3.6568e-13,3.1151e-13,4.3290e-13],
 [1.9524e-11,1.4866e-11,3.2192e-11, 2.0661e-11,1.8503e-11,2.6481e-11, 2.0456e-11,1.8438e-11,2.3269e-11, 1.6525e-11,1.4850e-11,1.8635e-11,
  1.2646e-11,1.1571e-11,1.3953e-11, 1.0283e-11,9.5248e-12,1.1147e-11, 7.8173e-12,7.2633e-12,8.4194e-12, 5.6494e-12,5.2664e-12,6.0550e-12,
  4.2016e-12,3.9093e-12,4.5022e-12, 3.0784e-12,2.8632e-12,3.3035e-12, 2.1934e-12,2.0376e-12,2.3523e-12, 1.5082e-12,1.3923e-12,1.6282e-12,
  9.8307e-13,8.9601e-13,1.0749e-12, 6.2200e-13,5.5517e-13,6.9798e-13, 3.8775e-13,3.3420e-13,4.5779e-13])]
# Mistele et al. 2024, Tabelle 1: Bin, z, log M_b, log R/kpc, v, σ_stat, σ_sys (km/s)
MIST = np.array([
 [1,.24,10.10,1.53,117.6,24.7,2.8],[1,.24,10.10,1.66,115.2,23.9,1.2],[1,.24,10.10,1.80,129.5,20.5,5.0],[1,.24,10.10,1.93,134.2,19.3,2.2],[1,.24,10.10,2.06,137.7,18.6,1.9],
 [2,.30,10.66,1.53,114.3,33.1,18.1],[2,.30,10.66,1.66,148.8,23.9,0.7],[2,.30,10.66,1.80,161.4,21.0,4.8],[2,.30,10.66,1.93,164.4,19.9,2.7],[2,.30,10.66,2.06,196.0,16.3,0.7],
 [3,.35,10.96,1.53,167.3,19.7,5.8],[3,.35,10.96,1.66,207.9,14.9,6.2],[3,.35,10.96,1.80,204.2,14.4,6.9],[3,.35,10.96,1.93,202.3,14.1,4.9],[3,.35,10.96,2.06,211.5,13.2,0.2],
 [4,.39,11.29,1.53,234.6,13.7,0.6],[4,.39,11.29,1.66,219.0,13.6,0.2],[4,.39,11.29,1.80,266.0,10.6,2.8],[4,.39,11.29,1.93,267.8,10.2,1.7],[4,.39,11.29,2.06,242.3,11.0,4.5]])

# ------------------------------------------------------------------ Regel mit äußerem Feld
MU, WMU = np.polynomial.legendre.leggauss(64)

def g_regel(gi, ge, lage):
    """Radialer Monopol der Regel; ge in m/s² (Newton-Feld der Umgebung)."""
    gi = np.asarray(gi, float)
    if ge == 0:
        return K.nu(gi, a0) * gi
    if lage == "senkrecht":
        return K.nu(np.hypot(gi, ge), a0) * gi
    if lage == "gleich":
        return K.nu(gi + ge, a0) * gi
    out = np.empty_like(gi)
    for k, x in enumerate(gi.ravel()):
        gt = np.sqrt(np.maximum(x * x + ge * ge + 2 * x * ge * MU, 1e-60))
        out.ravel()[k] = 0.5 * np.sum(WMU * K.nu(gt, a0) * (x + ge * MU))
    return out

# ------------------------------------------------------------------ Rest der Mode ξ je Schale
Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_von_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
rho_c = 3 * H0**2 / (8 * math.pi * G)
AI = 1e-3; TI = t_von_a(AI); RHO_I = Om * rho_c / AI**3

def schale(M, di, t_obs):
    """Lagrange-Schale der Masse M, Anfangsüberdichte di bei a = 10⁻³. Gibt (R(t_obs), I = ∫(1 − w)ω dt)."""
    ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
    vi = Hb(AI) * ri * (1 - di / 3)
    def f(t, y):
        R, V, I = y
        w = K.gewicht_w(3 * V / R, a0)
        return [V, -G * M / R**2 + OL * H0**2 * R, (1 - w) * math.sqrt(G * M / R**3)]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (TI, t_obs), [ri, vi, 0.0], events=ev, rtol=1e-8, atol=1e-6)
    if not len(s.t_events[0]):
        return s.y[0][-1], s.y[2][-1]
    t1, (R1, V1, I1) = s.t_events[0][0], s.y_events[0][0]
    rta = R1
    ev2 = lambda t, y: y[0] - rta / 2; ev2.terminal = True; ev2.direction = -1
    s2 = solve_ivp(f, (t1, t_obs), [R1, -1e-9, I1], events=ev2, rtol=1e-8, atol=1e-6)
    if not len(s2.t_events[0]):
        return s2.y[0][-1], s2.y[2][-1]
    t2, I2 = s2.t_events[0][0], s2.y_events[0][0][2]
    return rta / 2, I2 + (t_obs - t2) * math.sqrt(G * M / (rta / 2)**3)

def I_schale(r, M, t_obs):
    """Schale, die zur Zeit t_obs bei r liegt."""
    fn = lambda ld: schale(M, 10**ld, t_obs)[0] - r
    lo, hi = -6.0, 0.0
    if fn(lo) < 0:
        return schale(M, 10**lo, t_obs)[1], "exp"
    ld = brentq(fn, lo, hi, xtol=1e-6)
    return schale(M, 10**ld, t_obs)[1], None

mu_nfw = lambda x: math.log(1 + x) - x / (1 + x)
def massen(Mb):
    Mvir = Mb * (1e12 / 6e10); cc = 10
    rvir = (3 * Mvir / (4 * math.pi * 200 * rho_c)) ** (1 / 3)
    vf = K.v_flach(Mb, a0)
    return {"flach": lambda r: vf**2 * r / G, "NFW": lambda r: Mvir * mu_nfw(r * cc / rvir) / mu_nfw(cc)}, rvir

R_XI = np.geomspace(30, 5000, 26) * kpc
_cache = {}
def I_profil(Mb, prof, z):
    key = (Mb, prof, z)
    if key not in _cache:
        Mf, _ = massen(Mb)
        t_obs = t_von_a(1 / (1 + z))
        _cache[key] = np.array([I_schale(r, Mf[prof](r), t_obs)[0] for r in R_XI])
    return _cache[key]

def xi_R(r, Mb, prof, z, Ns):
    if Ns is None:
        return np.zeros_like(r)
    I = I_profil(Mb, prof, z)
    lI = np.interp(np.log(r), np.log(R_XI), np.log(np.maximum(I, 1e-300)), left=np.log(I[0]), right=np.log(I[-1]))
    return np.exp(-np.exp(lI) / Ns)

# ------------------------------------------------------------------ Linse
RG = np.geomspace(0.05 * kpc, 3000 * Mpc, 2400)

def M_linse(Mb, ge=0.0, lage="senkrecht", Ns=None, prof="flach", z=0.25, form="kumulativ", Mb_xi=None):
    """Mb_xi: Masse, mit der der Rest ξ_R der Schalen gerechnet wird (Standard Mb)."""
    gN = G * Mb / RG**2
    g = g_regel(gN, ge, lage)
    _, rvir = massen(Mb)
    xi = xi_R(RG, Mb if Mb_xi is None else Mb_xi, prof, z, Ns)
    if form == "lokal":
        nu_ = g / gN
        return gN * (nu_ - xi * ((nu_ - 1) - Q * np.minimum(RG / rvir, 1))) * RG**2 / G
    MX = (g - gN) * RG**2 / G
    dMX = np.gradient(MX, RG)
    dM = (1 - xi) * dMX + xi * Q * Mb / rvir * (RG < rvir)
    return Mb + np.concatenate([[0], np.cumsum(0.5 * (dM[1:] + dM[:-1]) * np.diff(RG))])

def g_4GdS(Mb, R, **kw):
    """4G·ΔΣ(R) aus M_L(<r)."""
    ML = M_linse(Mb, **kw)
    Mx = ML - Mb
    dM = np.gradient(Mx, RG)
    out = []
    for Rp in np.atleast_1d(R):
        m = RG > Rp
        rr, dd = RG[m], dM[m]
        sq = np.sqrt(1 - (Rp / rr)**2)
        M2 = Mb + np.interp(Rp, RG, Mx) + np.trapezoid(dd * (1 - sq), rr)
        # Σ(R) = ∫ dM/dr /(2π r √(r² − R²)) dr, Substitution r = R cosh u
        u = np.linspace(1e-6, math.acosh(RG[-1] / Rp) * 0.999, 4000)
        ru = Rp * np.cosh(u)
        Sig = np.trapezoid(np.interp(ru, RG, dM) / (2 * math.pi * ru), u)
        out.append(4 * G * (M2 / (math.pi * Rp**2) - Sig))
    return np.array(out)

def g_3D(Mb, r, **kw):
    return np.interp(r, RG, M_linse(Mb, **kw)) * G / r**2

# ------------------------------------------------------------------ Auswertung
def chi2_bin(b, kw, rmax_mpc, fM=1.0):
    sel = R_MPC <= rmax_mpc + 1e-9
    d = KIDS[b][sel]
    mod = g_4GdS(MB[b] * fM, R_MPC[sel] * Mpc, Mb_xi=MB[b], **kw)
    sig = 0.5 * (np.log10(d[:, 2]) - np.log10(d[:, 1]))
    return np.sum(((np.log10(mod) - np.log10(d[:, 0])) / sig)**2), sel.sum()

def chi2_kids(bins_kw, rmax_mpc, fM=1.0):
    ch, n = 0.0, 0
    for b in range(4):
        x, m = chi2_bin(b, bins_kw, rmax_mpc, fM); ch += x; n += m
    return ch, n

FM_GITTER = 10 ** np.arange(-0.2, 0.61, 0.025)
def chi2_frei(kw, rmax_mpc):
    """Je Bin eigene Massenskala f_M (heißes Gas, M*-Systematik): gibt χ², n, f_M je Bin."""
    ch, n, fs = 0.0, 0, []
    for b in range(4):
        vals = [chi2_bin(b, kw, rmax_mpc, f)[0] for f in FM_GITTER]
        k = int(np.argmin(vals)); ch += vals[k]; n += (R_MPC <= rmax_mpc + 1e-9).sum(); fs.append(FM_GITTER[k])
    return ch, n, fs

if __name__ == "__main__":
    print(f"a₀ = {a0:.2e} m/s², Q = {Q:.3f}, t₀ = {t_von_a(1)/Gyr:.2f} Mrd. J.")
    print("Massen der Bins (M* → M_b = M*(1+f_kalt)), r*, v_flach:")
    for b in range(4):
        print(f"  Bin {b+1}: M* = {MSTAR[b]:.1e}, M_b = {MB[b]:.2e} M☉, r* = {K.r_stern(MB[b]*Ms)/kpc:.1f} kpc, "
              f"v_flach = {K.v_flach(MB[b]*Ms)/1e3:.0f} km/s, r_vir(Ansatz 03) = {massen(MB[b]*Ms)[1]/kpc:.0f} kpc")
    MB = [m * Ms for m in MB]

    # 0 Anschluss an 03: Milchstraße, lokaler Ansatz, z = 0
    print("\n0) Anschluss an 03 (Milchstraße 6·10¹⁰, z = 0, Ansatz lokal): δg/g = g/g_RAR − 1")
    Mmw = 6e10 * Ms
    for prof in ("NFW", "flach"):
        row = []
        for Ns in (0.3, 2.4):
            for rk in (100, 200, 300):
                r = np.array([rk * kpc])
                row.append(g_3D(Mmw, r, Ns=Ns, prof=prof, z=0.0, form="lokal")[0] / g_3D(Mmw, r)[0] - 1)
        print(f"  {prof:5s}: N* = 0,3: " + " ".join(f"{x:+.3f}" for x in row[:3]) + " | N* = 2,4: " + " ".join(f"{x:+.3f}" for x in row[3:]) + "  (100/200/300 kpc)")
    print("  kumulativ (Gaußfluss der Schalen), 3D, z = 0: δg/g bei 100/200/300/500/1000 kpc")
    for prof in ("flach", "NFW"):
        for Ns in (0.3, 2.4):
            r = np.array([100, 200, 300, 500, 1000]) * kpc
            x = g_3D(Mmw, r, Ns=Ns, prof=prof, z=0.0) / g_3D(Mmw, r) - 1
            print(f"    {prof:5s} N* = {Ns}: " + " ".join(f"{v:+.3f}" for v in x))
    print("  Rest ξ_R (z = 0, flach) je Schale, N* = 2,4:",
          " ".join(f"{rk}:{xi_R(np.array([rk*kpc]), Mmw, 'flach', 0.0, 2.4)[0]:.2g}" for rk in (100, 200, 300, 500, 1000, 2000)))

    # 1 Messpunkte: g_obs/g_bar gemessen und im Modell, Bin 3 (M_b ≈ Milchstraße)
    for b in (2, 0, 3):
        Mb = MB[b]
        gbar = G * Mb / (R_MPC * Mpc)**2
        d = KIDS[b]
        print(f"\n1) Bin {b+1} (M_b = {Mb/Ms:.2e}): g_obs/g_bar = 4GΔΣ/g_bar; Modell: Regel ohne/mit g_e (Kugelmittel), Mode kumulativ (flach, z = 0,25)")
        cols = [("Regel", {}),
                ("ge.003K", dict(ge=0.003 * a0, lage="kugel")),
                ("ge.01⊥", dict(ge=0.01 * a0, lage="senkrecht")),
                ("ge.01K", dict(ge=0.01 * a0, lage="kugel")),
                ("ge.03⊥", dict(ge=0.03 * a0, lage="senkrecht")),
                ("N*0.3", dict(Ns=0.3)), ("N*2.4", dict(Ns=2.4)),
                ("N*2.4lok", dict(Ns=2.4, form="lokal")), ("N*2.4NFW", dict(Ns=2.4, prof="NFW"))]
        mods = {k: g_4GdS(Mb, R_MPC * Mpc, **kw) for k, kw in cols}
        print("  R/kpc  g_bar      gemessen (−/+)          " + " ".join(f"{k:>8s}" for k, _ in cols))
        for i, R in enumerate(R_MPC):
            print(f"  {R*1e3:6.0f} {gbar[i]:.2e}  {d[i,0]/gbar[i]:7.1f} ({d[i,1]/gbar[i]:6.1f}–{d[i,2]/gbar[i]:6.1f})  "
                  + " ".join(f"{mods[k][i]/gbar[i]:8.1f}" for k, _ in cols))

    # 2 χ²
    print("\n2) χ² (log g, σ = halbe 16–84-Breite, Punkte als unabhängig gezählt) gegen KiDS, 4 Bins; Fenster R ≤ 0,30 / 1,04 / 2,60 Mpc")
    varianten = [("Regel, g_e = 0", {})]
    for ge in (0.0001, 0.0003, 0.001, 0.003, 0.01, 0.02, 0.03):
        for lage in ("senkrecht", "gleich", "kugel"):
            varianten.append((f"g_e = {ge:.4f}·a₀ {lage}", dict(ge=ge * a0, lage=lage)))
    for prof in ("flach", "NFW"):
        for form in ("kumulativ", "lokal"):
            for Ns in (0.3, 0.6, 1.2, 2.4):
                varianten.append((f"Mode N* = {Ns} {prof} {form}", dict(Ns=Ns, prof=prof, form=form)))
    varianten.append(("Mode N* = 2.4 flach kumulativ, z = 0", dict(Ns=2.4, z=0.0)))
    for name, kw in varianten:
        res = [chi2_kids(kw, rm) for rm in (0.3035, 1.0365, 2.6035)]
        print(f"  {name:42s} " + "  ".join(f"{c:7.1f}/{n}" for c, n in res))

    # 2b freie Massenskala (heißes Gas, M*-Systematik), Regel ohne g_e
    print("\n2b) Regel g_e = 0 mit gemeinsamer Massenskala f_M (M_b → f_M·M_b), χ² bis 0,30 / 1,04 / 2,60 Mpc")
    for fM in (0.7, 1.0, 1.3, 1.6, 2.0, 2.5):
        res = [chi2_kids({}, rm, fM) for rm in (0.3035, 1.0365, 2.6035)]
        print(f"  f_M = {fM:.1f}: " + "  ".join(f"{c:7.1f}/{n}" for c, n in res))

    # 2c je Bin freie Massenskala: Formtest
    print("\n2c) Formtest: je Bin eigene Massenskala f_M angepasst (Rest ξ_R mit der Nennmasse); χ²/n und f_M der Bins 1–4")
    for name, kw in [("Regel, g_e = 0", {}), ("g_e = 0,0003 Kugelmittel", dict(ge=0.0003 * a0, lage="kugel")),
                     ("g_e = 0,001 Kugelmittel", dict(ge=0.001 * a0, lage="kugel")),
                     ("g_e = 0,003 Kugelmittel", dict(ge=0.003 * a0, lage="kugel")),
                     ("g_e = 0,003 senkrecht", dict(ge=0.003 * a0, lage="senkrecht")),
                     ("g_e = 0,01 senkrecht", dict(ge=0.01 * a0, lage="senkrecht")),
                     ("Mode N* = 0,3 flach kumulativ", dict(Ns=0.3)), ("Mode N* = 0,6 flach kumulativ", dict(Ns=0.6)),
                     ("Mode N* = 1,2 flach kumulativ", dict(Ns=1.2)), ("Mode N* = 2,4 flach kumulativ", dict(Ns=2.4)),
                     ("Mode N* = 0,3 NFW kumulativ", dict(Ns=0.3, prof="NFW")), ("Mode N* = 2,4 NFW kumulativ", dict(Ns=2.4, prof="NFW"))]:
        out = []
        for rm in (0.3035, 1.0365, 2.6035):
            ch, n, fs = chi2_frei(kw, rm)
            out.append(f"{ch:6.1f}/{n} (" + " ".join(f"{f:.2f}" for f in fs) + ")")
        print(f"  {name:32s} " + "  ".join(out))

    # 3 Mistele 2024 (3D, v_circ), Regel und äußeres Feld
    print("\n3) Mistele 2024 (v_circ, 34–115 kpc, 20 Punkte): χ² in v")
    for name, kw in [("Regel", {}), ("g_e = 0,01 senkrecht", dict(ge=0.01 * a0, lage="senkrecht")),
                     ("g_e = 0,01 Kugelmittel", dict(ge=0.01 * a0, lage="kugel")),
                     ("g_e = 0,03 senkrecht", dict(ge=0.03 * a0, lage="senkrecht")),
                     ("Mode N* = 2,4 flach", dict(Ns=2.4))]:
        ch, lines = 0.0, []
        for row in MIST:
            Mb = 10**row[2] * Ms; r = 10**row[3] * kpc
            v = math.sqrt(g_3D(Mb, np.array([r]), z=row[1], **kw)[0] * r) / 1e3
            ch += ((v - row[4]) / math.hypot(row[5], row[6]))**2
            lines.append(v)
        print(f"  {name:24s} χ² = {ch:5.1f}/20   Modell v [km/s] je Bin bei 115 kpc: " + " ".join(f"{lines[i]:.0f}" for i in (4, 9, 14, 19)))
    print("  gemessen bei 115 kpc: " + " ".join(f"{MIST[i,4]:.0f}±{math.hypot(MIST[i,5],MIST[i,6]):.0f}" for i in (4, 9, 14, 19)))

    # 4 Einzelwerte für die Handprüfung
    Mb = MB[0]; R = 1.0365 * Mpc
    print("\n4) Einzelwerte für die Handprüfung (Bin 1, R = 1,0365 Mpc):")
    print(f"  g_bar = {G*Mb/R**2:.4e}, Regel 4GΔΣ = {g_4GdS(Mb, [R])[0]:.4e}, 3D g = {g_3D(Mb, np.array([R]))[0]:.4e}, √(g_bar a₀) = {math.sqrt(G*Mb/R**2*a0):.4e}")
    print(f"  g_e = 0,01 a₀ senkrecht: 3D g/g_bar = {g_3D(Mb, np.array([R]), ge=0.01*a0, lage='senkrecht')[0]/(G*Mb/R**2):.3f}")
    print(f"  ξ_R(1 Mpc, flach, z = 0,25, N* = 2,4) = {xi_R(np.array([R]), Mb, 'flach', 0.25, 2.4)[0]:.3f}, "
          f"I = ∫(1−w)ω dt = {math.exp(np.interp(math.log(R), np.log(R_XI), np.log(I_profil(Mb,'flach',0.25)))):.3f}")
    Mf, _ = massen(Mb)
    print(f"  M_flach(1 Mpc) = {Mf['flach'](R)/Ms:.3e} M☉, EdS-Umkehrzeit mit r_ta = 2R: {math.pi*math.sqrt((2*R)**3/(8*G*Mf['flach'](R)))/Gyr:.1f} Mrd. J.")

    # 5 Gegenweg für I = ∫(1 − w)ω dt: EdS-Zykloide ohne Λ (anderer Hintergrund, geschlossene Bahn)
    t_obs = t_von_a(1 / 1.25); M = Mf["flach"](R)
    def I_eds(eta_obs):
        B = t_obs / (eta_obs - math.sin(eta_obs)); A = (G * M * B**2) ** (1 / 3)
        et = np.linspace(1e-4, eta_obs, 200001)
        Rr = A * (1 - np.cos(et)); dt = B * (1 - np.cos(et))
        th = 3 * (np.sin(et) / (1 - np.cos(et))) / B
        thc = a0 / c
        w = np.where(th > 0, th**2 / (th**2 + thc**2), 0.0)
        return A * (1 - math.cos(eta_obs)), np.trapezoid((1 - w) * np.sqrt(G * M / Rr**3) * dt, et)
    eta = brentq(lambda e: I_eds(e)[0] - R, 0.1, 1.5 * math.pi)   # bis r_ta/2 (η = 3π/2)
    print(f"\n5) Gegenweg EdS-Zykloide (Bin 1, 1 Mpc, z = 0,25): η = {eta:.3f}, I = {I_eds(eta)[1]:.3f} (ODE mit Λ: siehe 4)")
