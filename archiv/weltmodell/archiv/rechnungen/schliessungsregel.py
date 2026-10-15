# -*- coding: utf-8 -*-
"""
FTM - Schliessungsregel (Vorschlag A13): Projektion z' = (I - w Pi_perp) F mit w = w(theta),
theta = nabla_mu u^mu (Expansionsrate des Zeitfelds).

Abschnitte der Ausgabe:
 0  Groessen, Wahl von w(theta)
 1  Wirkung/Feldgleichungen: Schallgeschwindigkeit, 0i-Probe (Impulserhaltung), 00-Probe (Hintergrund)
 2  Kosmos: w im Hintergrund, Spielzeug-Boltzmann mit w(a), lokales theta der Moden,
    Wachstum f sigma8(z), Einsatz der ersten Abweichung (k, z)
 3  Galaxie: spaerischer Kollaps, theta(t), w(t), Uebergabe W(t), Bilanz bei R_Q, g/g_RAR
 4  Galaxienhaufen und Bullet Cluster
 5  Fenster fuer die Uebergaberate, Vorhersagen

Nur numpy/math. Aufruf: python schliessungsregel.py  (schreibt schliessungsregel.txt daneben)
Spielzeug-Boltzmann: Ma & Bertschinger (1995), konforme Newton-Eichung, psi = phi, Neutrinos als
ideales Fluid, enge Photon-Baryon-Kopplung bis z* = 1090 (Spielzeug; verglichen werden Varianten
untereinander mit identischer Methode).
"""
import math
import os
import numpy as np

OUT = []


def P(s=""):
    OUT.append(str(s))


def hdr(s):
    P()
    P("=" * 100)
    P(s)
    P("=" * 100)


# ------------------------------------------------------------------ Eingaben
h = 0.6736
wb, wc = 0.02237, 0.1200
wg, wr = 2.4729e-5, 4.1835e-5
Ob, Oc, Or = wb / h**2, wc / h**2, wr / h**2
Om = Ob + Oc
OL = 1.0 - Om
fg = wg / wr
H0c = h / 2997.92458               # H0/c in 1/Mpc
zs = 1090.0
a_s = 1.0 / (1.0 + zs)
Q = wc / wb                        # 5,364: Gewicht des projizierten Flusses
q = Q / (1.0 + Q)                  # 0,8429
YHe = 0.245
sigT = 6.6524587e-29
mp = 1.67262192e-27
rhoc_h2 = 1.87847e-26
MpcM = 3.0856776e22
nH0 = (1 - YHe) * wb * rhoc_h2 / mp
tau0 = nH0 * sigT * MpcM
G = 6.674e-11
c = 2.99792458e8
Msun = 1.989e30
kpc = 3.0856776e19
a0 = float(__import__("os").environ.get("FTM_A0", "1.2e-10"))
H0s = 100 * h * 1e3 / MpcM          # 1/s
Gyr = 3.15576e16
tH_Gyr = 1.0 / H0s / Gyr
rho_c0 = 3 * H0s**2 / (8 * math.pi * G)
th_c = a0 / c                       # 1/s: Schwellrate theta_c
thc_H0 = th_c / H0s


def w_of_theta(th_H0):
    """w = theta_+^2/(theta_+^2 + theta_c^2), theta in Einheiten H0; theta <= 0 -> w = 0."""
    tp = max(th_H0, 0.0)
    return tp * tp / (tp * tp + thc_H0 * thc_H0)


def nu(y):
    return 1.0 / (1.0 - math.exp(-math.sqrt(y)))


def E(a):
    return math.sqrt(Om * a**-3 + OL)


# ==========================================================================================
hdr("0  Groessen und Wahl von w(theta)")
P(f"  Q = omega_c/omega_b = {Q:.4f},  q = Q/(1+Q) = {q:.5f},  nu_kos = 1/(1-q) = {1/(1-q):.4f},  y_Q = (ln q)^2 = {math.log(q)**2:.5f}")
P(f"  theta_c = a0/c = {th_c:.4e} 1/s = {thc_H0:.4f} H0   (a0 = 1,2e-10 m/s^2; a0/c ist die Atemrate des Torus pro Radiant, cH0/2pi -> {1/(2*math.pi):.4f} H0)")
P(f"  H0^-1 = {tH_Gyr:.3f} Gyr,  rho_c0 = {rho_c0:.4e} kg/m^3")
P("  w(theta) = theta_+^2/(theta_+^2+theta_c^2): C^1-glatt, w = 1 fuer theta >> theta_c, w = 0 fuer theta <= 0 (Umkehr, Einfall, gebunden)")

# ------------------------------------------------------------------ Spielzeug-Boltzmann


class BG:
    def __init__(self, Omx):
        self.Om = Omx
        self.OL = 1.0 - Omx - Or

    def E(self, a):
        return np.sqrt(Or * a**-4 + self.Om * a**-3 + self.OL)

    def Hc(self, a):
        return a * H0c * self.E(a)

    def frac(self, a):
        E2 = self.E(a) ** 2
        return Or * a**-4 / E2, Ob * a**-3 / E2, self.Om * a**-3 / E2, self.OL / E2

    def DM(self, a):
        x = np.linspace(math.log(a), 0.0, 6000)
        aa = np.exp(x)
        return float(np.trapezoid(1.0 / self.Hc(aa), x))

    def eta(self, a, n=4000):
        x = np.linspace(math.log(1e-10), math.log(a), n)
        aa = np.exp(x)
        return float(np.trapezoid(1.0 / self.Hc(aa), x) + 1e-10 / (H0c * math.sqrt(Or)))

    def kD(self, a):
        x = np.linspace(math.log(1e-8), math.log(a), 6000)
        aa = np.exp(x)
        R = 0.75 * wb / wg * aa
        tdot = tau0 / aa**2
        f = (R**2 / (1 + R) + 16.0 / 15.0) / (6 * (1 + R) * tdot) / self.Hc(aa)
        return 1.0 / math.sqrt(float(np.trapezoid(f, x)))


NV = 7   # 0 phi, 1 d_g, 2 th (Photon-Baryon), 3 d_nu, 4 th_nu, 5 d_P (proj. Fluss), 6 th_P


def wval(cfg, a):
    w = cfg.get("w", 1.0)
    return w(a) if callable(w) else w


def build_A(a, k, bg, cfg):
    Hc = bg.Hc(a)
    fr, fb, fm, fL = bg.frac(a)
    fgam, fnu = fr * fg, fr * (1 - fg)
    R = 0.75 * wb / wg * a
    nk = k.size
    A = np.zeros((nk, NV, NV))
    k2 = k * k
    w = wval(cfg, a)
    cons = cfg.get("cons", True)
    # Mitnahme: Anteil (1-w) des Flusses beschleunigt mit dem Photon-Baryon-Fluid; Rueckwirkung s = I/(I+L)
    eps = (1 - w) * Q * R / (1 + R)
    s = 1.0 / (1.0 + eps) if cons else 1.0
    S = np.zeros((nk, NV))
    S[:, 1] = fgam + 0.75 * fb
    S[:, 3] = fnu
    S[:, 5] = Q * fb
    rowphi = -0.5 * S
    rowphi[:, 0] += -1.0 - k2 / (3 * Hc**2)
    A[:, 0, :] = rowphi
    A[:, 1, :] = 4 * rowphi
    A[:, 1, 2] += -4.0 / 3.0 / Hc
    # th' = -th + k^2 phi/Hc + s f,  f = th/(1+R) + k^2 d_g/(4(1+R)Hc)  (nichtgravitative Beschleunigung)
    A[:, 2, 2] = -1.0 + s / (1 + R)
    A[:, 2, 1] = s * k2 / (4 * (1 + R) * Hc)
    A[:, 2, 0] = k2 / Hc
    A[:, 3, :] = 4 * rowphi
    A[:, 3, 4] += -4.0 / 3.0 / Hc
    A[:, 4, 3] = k2 / (4 * Hc)
    A[:, 4, 0] = k2 / Hc
    # projizierter Fluss: d_P' = -th_P/Hc + 3 phi';  th_P' = -th_P + k^2 phi/Hc + (1-w) s f
    A[:, 5, :] = 3 * rowphi
    A[:, 5, 6] += -1.0 / Hc
    A[:, 6, 6] = -1.0
    A[:, 6, 0] = k2 / Hc
    A[:, 6, 2] = (1 - w) * s / (1 + R)
    A[:, 6, 1] = (1 - w) * s * k2 / (4 * (1 + R) * Hc)
    return A


def resid_0i(a, k, y, bg, cfg):
    """0i-Gleichung: k^2 (dphi/dN + phi) = 1,5 Hc sum Omega_i (1+w_i) theta_i ; relative Abweichung."""
    A = build_A(a, k, bg, cfg)
    dphi = np.einsum("kj,kj->k", A[:, 0, :], y)
    Hc = bg.Hc(a)
    fr, fb, fm, fL = bg.frac(a)
    fgam, fnu = fr * fg, fr * (1 - fg)
    L = k**2 * (dphi + y[:, 0])
    terms = np.array([(4 / 3) * fgam * y[:, 2], fb * y[:, 2], (4 / 3) * fnu * y[:, 4], Q * fb * y[:, 6]])
    Rr = 1.5 * Hc * terms.sum(axis=0)
    return np.abs(L - Rr) / (1.5 * Hc * np.abs(terms).sum(axis=0))


def run(cfg, k, a_end=a_s, dN=0.0015, ai=1e-7, checks=()):
    bg = BG(Om)
    eta_i = bg.eta(ai)
    y = np.zeros((k.size, NV))
    y[:, 0] = 1.0
    y[:, 1] = -2.0
    y[:, 2] = 0.5 * k**2 * eta_i
    y[:, 3] = -2.0
    y[:, 4] = 0.5 * k**2 * eta_i
    y[:, 5] = -1.5
    y[:, 6] = 0.5 * k**2 * eta_i
    N = np.arange(math.log(ai), math.log(a_end), dN)
    N = np.append(N, math.log(a_end))
    I = np.eye(NV)[None, :, :]
    An = build_A(math.exp(N[0]), k, bg, cfg)
    res = {}
    for j in range(len(N) - 1):
        hN = N[j + 1] - N[j]
        An1 = build_A(math.exp(N[j + 1]), k, bg, cfg)
        rhs = y + 0.5 * hN * np.einsum("kij,kj->ki", An, y)
        y = np.linalg.solve(I - 0.5 * hN * An1, rhs[..., None])[..., 0]
        An = An1
        for ac in checks:
            if ac not in res and N[j + 1] >= math.log(ac) - 1e-12:
                res[ac] = resid_0i(math.exp(N[j + 1]), k, y, bg, cfg)
    return y, bg, res


KD_PLANCK = 0.1409


def spectrum(y, k, bg):
    Th = y[:, 1] / 4 + y[:, 0]
    Th1 = y[:, 2] / (3 * k)
    kD = KD_PLANCK
    return (Th**2 + Th1**2) * np.exp(-2 * (k / kD) ** 2)


def peaks(Pk, k, n=5, kmin=0.008):
    idx = [i for i in range(1, len(k) - 1) if k[i] > kmin and Pk[i] > Pk[i - 1] and Pk[i] >= Pk[i + 1]]
    return [(k[i], Pk[i]) for i in idx[:n]]


def ratios(pk):
    H = [p[1] for p in pk]
    while len(H) < 5:
        H.append(float("nan"))
    return H[1] / H[0], H[2] / H[0], H[2] / H[1], H[3] / H[2], H[4] / H[3]


kgrid = np.linspace(0.002, 0.32, 520)

# ==========================================================================================
hdr("1  Wirkung: Stabilitaet und Erhaltung (Mitnahme mit Rueckwirkung)")
Rst = 0.75 * wb / wg * a_s
P(f"  R(z*) = {Rst:.4f}.  Schallgeschwindigkeit des Photon-Baryon-Fluids mit mitgenommenem Anteil (1-w) des Flusses:")
P("  c_s^2 = 1/(3(1 + R + (1-w) Q R))   (w = 1: Standard; w = 0: Fluss voll mitgenommen)")
for w in (1.0, 0.99, 0.9, 0.5, 0.0):
    cs2 = 1 / (3 * (1 + Rst + (1 - w) * Q * Rst))
    P(f"   w = {w:4.2f}:  c_s^2(z*) = {cs2:.4f}  (>= 0)")
P("  projizierter Fluss (freier Anteil): c_s^2 = 0 exakt (Nebenbedingung (d chi)^2 = -1), kein Gradiententerm, Energiedichte mu >= 0")

checks = (1e-5, 1e-4, 5e-4, a_s)
sel = (kgrid > 0.01) & (kgrid < 0.3)
P()
P("  0i-Probe k^2(phi' + phi) = 1,5 Hc sum Omega_i (1+w_i) theta_i  (00-Gleichung treibt phi; 0i ist unabhaengige Kontrolle);"
  "  Abweichung normiert auf 1,5 Hc sum |Omega_i (1+w_i) theta_i|")
P("  Median / Maximum der relativen Abweichung ueber 0,01 < k < 0,3 /Mpc:")
runs = {}
for nm, cfg in (("w = 1", dict(w=1.0)),
                ("w = 0,9 mit Rueckwirkung", dict(w=0.9, cons=True)),
                ("w = 0,5 mit Rueckwirkung", dict(w=0.5, cons=True)),
                ("w = 0,9 ohne Rueckwirkung", dict(w=0.9, cons=False))):
    y, bg, res = run(cfg, kgrid, checks=checks)
    runs[nm] = (y, bg, cfg)
    P(f"   {nm:28s} " + "  ".join(f"a={ac:.0e}: {np.median(res[ac][sel]):.1e}/{np.max(res[ac][sel]):.1e}" for ac in checks))
P("  => Mit Rueckwirkung ist die 0i-Gleichung fuer jedes w auf Diskretisierungsniveau erfuellt (wie bei w = 1).")

# ==========================================================================================
hdr("2  Kosmos")
P("  2a  w im Hintergrund, theta = 3H:")
for z in (0, 0.5, 1, 2, 5, 10, 100, 1090):
    a = 1 / (1 + z)
    th = 3 * E(a)
    P(f"   z = {z:6g}:  theta/H0 = {th:10.3f}   1 - w = {1-w_of_theta(th):.3e}")

# Lokales theta der Moden bei z*: theta_lok/theta_bg = 1 - delta_P'/3
yref, bgref, _ = runs["w = 1"]
phip = (2 / 3) * math.sqrt(2.1e-9)
Aref = build_A(a_s, kgrid, bgref, dict(w=1.0))
dP = np.abs(np.einsum("kj,kj->k", Aref[:, 5, :], yref)) * phip
P(f"  2b  Moden bei z*: max |delta_P'|/3 (Delta pro ln k, 0,002-0,32/Mpc) = {dP.max()/3:.2e}; theta_lok = 3H(1 - delta'/3)")
th_s = 3 * E(a_s)
wl = w_of_theta(th_s * (1 - 3 * dP.max() / 3))
P(f"      auch bei 10 sigma: 1 - w(z*) = {1-w_of_theta(th_s*(1-10*dP.max()/3)):.2e}  (Hintergrund {1-w_of_theta(th_s):.2e})")

P()
P("  2c  Spielzeug-Boltzmann: Peakverhaeltnisse, Abweichung von H3/H2 gegen w = 1")
Pk_ref = spectrum(yref, kgrid, bgref)
rref = ratios(peaks(Pk_ref, kgrid))


def line(nm, y, bg):
    Pk = spectrum(y, kgrid, bg)
    rr = ratios(peaks(Pk, kgrid))
    P(f"   {nm:34s} H2/H1 {rr[0]:6.3f}  H3/H1 {rr[1]:6.3f}  H3/H2 {rr[2]:6.3f}  H4/H3 {rr[3]:6.3f}  H5/H4 {rr[4]:6.3f}  "
      f"|dH3/H2| {abs(rr[2]-rref[2])/rref[2]*100:6.3f}%  max|dH_i/H_i| {max(abs(rr[i]-rref[i])/rref[i] for i in range(5))*100:6.3f}%")


line("w = 1", yref, bgref)
wtheta = lambda a: w_of_theta(3 * math.sqrt(Or * a**-4 + Om * a**-3 + OL - Or))
y2, bg2, _ = run(dict(w=wtheta, cons=True), kgrid)
line("w(theta_bg), theta_c = a0/c", y2, bg2)
for w in (0.99, 0.97, 0.9):
    y3, bg3, _ = run(dict(w=w, cons=True), kgrid)
    line(f"w = {w} konstant (Empfindlichkeit)", y3, bg3)
y3, bg3, _ = runs["w = 0,5 mit Rueckwirkung"]
line("w = 0,5 konstant (Empfindlichkeit)", y3, bg3)
P(f"   Hintergrund: Omega_m = Omega_b(1+Q) = {Om:.4f}; Fluss und Uebergabe sind druckfrei -> a^3 rho_ges = const exakt")


# ------------------------------------------------------------------ Hintergrund-Zeit, lineares Wachstum (ohne Strahlung)
lnag = np.linspace(math.log(1e-4), 0.0, 20001)
ag = np.exp(lnag)
Eg = np.sqrt(Om * ag**-3 + OL)
tg = np.concatenate([[0.0], np.cumsum(0.5 * (1 / Eg[1:] + 1 / Eg[:-1]) * np.diff(lnag))]) + (2 / 3) * 1e-4**1.5 / math.sqrt(Om)
t0 = float(tg[-1])


def t_of_z(z):
    return float(np.interp(math.log(1 / (1 + z)), lnag, tg))


def z_of_t(t):
    return float(1 / math.exp(np.interp(t, tg, lnag)) - 1)


def growth(eps_fun=None):
    """D'' + (2 + dlnE/dlna) D' = 1,5 Om(a) (1 - eps) D, Start D = a bei a = 1e-3."""
    N = np.linspace(math.log(1e-3), 0.0, 6001)
    D, Dp = 1e-3, 1e-3
    Ds, fs = [D], [1.0]
    for j in range(len(N) - 1):
        hN = N[j + 1] - N[j]

        def rhs(n, s):
            a = math.exp(n)
            E2 = Om * a**-3 + OL
            Oma = Om * a**-3 / E2
            dlnE = -1.5 * Oma
            e = eps_fun(a) if eps_fun else 0.0
            return np.array([s[1], -(2 + dlnE) * s[1] + 1.5 * Oma * (1 - e) * s[0]])
        st = np.array([D, Dp])
        k1 = rhs(N[j], st); k2 = rhs(N[j] + hN / 2, st + hN / 2 * k1)
        k3 = rhs(N[j] + hN / 2, st + hN / 2 * k2); k4 = rhs(N[j] + hN, st + hN * k3)
        st = st + hN / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        D, Dp = st
        Ds.append(D); fs.append(Dp / D)
    return N, np.array(Ds), np.array(fs)


Ngr, Dl, fl = growth()


def Dlin(a):
    return float(np.interp(math.log(a), Ngr, Dl))


# Uebergabe im Hintergrund: d eps/dlna = (1-w_bg) sqrt(Om(a)/2) (1-eps) / N*, hier N* = 1 (obere Schranke, eps ~ 1/N*)
epsN = [0.0]
for j in range(len(lnag) - 1):
    a = ag[j]
    E2 = Om * a**-3 + OL
    rate = (1 - w_of_theta(3 * math.sqrt(E2))) * math.sqrt(Om * a**-3 / E2 / 2)
    epsN.append(epsN[-1] + rate * (1 - epsN[-1]) * (lnag[j + 1] - lnag[j]))
epsN = np.array(epsN)
eps_fun = lambda a: float(np.interp(math.log(a), lnag, epsN))
P()
P("  2d  Uebergabe im Hintergrund (Rate (1-w) Omega/N*, Omega = sqrt(G M/R^3) = H sqrt(Om(a)/2)); Werte fuer N* = 1, eps ~ 1/N*:")
for z in (10, 2, 1, 0.5, 0):
    P(f"      z = {z:4g}:  eps = {eps_fun(1/(1+z)):.3e}")
_, D2, f2 = growth(eps_fun)
s8 = 0.811
P("  2e  Wachstum f sigma8(z), sigma8 = 0,811 (Planck 2018, aus dem Gedaechtnis); Modell: Anteil eps als nicht klumpend (konservativ)")
P("      z      f_LCDM   fs8_LCDM   fs8_Modell(N*=1)   rel. Abw.")
for z in (0.0, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0):
    n = math.log(1 / (1 + z))
    fa = float(np.interp(n, Ngr, fl)); Da = float(np.interp(n, Ngr, Dl))
    fb_ = float(np.interp(n, Ngr, f2)); Db = float(np.interp(n, Ngr, D2))
    x1 = fa * s8 * Da / Dl[-1]
    x2 = fb_ * s8 * Db / Dl[-1]      # gleiche Normierung der Anfangsamplitude (CMB-normiert)
    P(f"     {z:4.1f}   {fa:.4f}   {x1:.4f}     {x2:.4f}           {(x2-x1)/x1:+.2e}")

# ------------------------------------------------------------------ spaerischer Kollaps (Top-Hat) im Hintergrund Om, OL
# x = R / R_L (R_L: mitbewegter Lagrange-Radius), t in 1/H0:  x'' = -Om/(2 x^2) + OL x


def tophat(di, ai=0.002, tstop=None, rec=False):
    x = ai * (1 - di / 3)
    v = x * E(ai) * (1 - di / (3 * (1 + di)))
    t = t_of_z(1 / ai - 1)
    tstop = t0 if tstop is None else tstop
    acc = lambda x: -Om / (2 * x * x) + OL * x
    Ipre = 0.0
    t_ta = x_ta = None
    hist = []
    while t < tstop:
        dt = 2e-3 * min(x**1.5 / math.sqrt(Om), 0.05)
        if t_ta is not None:
            dt = min(dt, 2e-3 * x / max(abs(v), 1e-30))
        # Uebergabe-Integrand (1-w) Omega vor dem Kollaps
        xe = x if t_ta is None else max(x, 0.5 * x_ta)     # Einfall: Zone endet bei r_vir = r_ta/2
        Ipre += (1 - w_of_theta(3 * v / x)) * math.sqrt(Om / (2 * xe**3)) * dt
        k1x, k1v = v, acc(x)
        k2x, k2v = v + dt / 2 * k1v, acc(x + dt / 2 * k1x)
        k3x, k3v = v + dt / 2 * k2v, acc(x + dt / 2 * k2x)
        k4x, k4v = v + dt * k3v, acc(x + dt * k3x)
        x += dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v += dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        t += dt
        if rec:
            hist.append((t, x, v))
        if t_ta is None and v < 0:
            t_ta, x_ta = t, x
        if t_ta is not None and x < 0.01 * x_ta:
            return dict(t_ta=t_ta, x_ta=x_ta, t_c=t, Ipre=Ipre, hist=hist, di=di)
    return dict(t_ta=t_ta, x_ta=x_ta, t_c=None, Ipre=Ipre, hist=hist, di=di)


def collapse_at(zc):
    tc = t_of_z(zc)
    lo, hi = 1e-3, 0.2
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        r = tophat(mid, tstop=tc * 1.5 + 0.05)
        if r["t_c"] is None or r["t_c"] > tc:
            lo = mid
        else:
            hi = mid
    return tophat(hi, tstop=tc * 1.5 + 0.05)


def RL(M):
    return (3 * M / (4 * math.pi * Om * rho_c0)) ** (1 / 3)


def handover_I(r, tobs):
    """Exponent der Uebergabe (ohne 1/N*) bis tobs; nach dem Kollaps theta = 0, w = 0, Omega = Omega_vir."""
    xv = r["x_ta"] / 2
    Ov = math.sqrt(Om / (2 * xv**3))
    return r["Ipre"] + Ov * max(tobs - r["t_c"], 0.0), Ov


# lineares delta bei Umkehr und Kollaps (Kontrolle gegen 1,06 / 1,686)
hdr("3  Galaxie: spaerischer Kollaps, theta_lok(t), w(t), Uebergabe W(t)")
for zc in (3.0, 0.0):
    r = collapse_at(zc)
    ai = 0.002
    dta = r["di"] * Dlin(math.exp(np.interp(r["t_ta"], tg, lnag))) / Dlin(ai)
    dc = r["di"] * Dlin(math.exp(np.interp(r["t_c"], tg, lnag))) / Dlin(ai)
    P(f"  Kontrolle Kollaps bei z_c = {zc}: lineares delta bei Umkehr {dta:.4f}, bei Kollaps {dc:.4f}")
    if zc == 3.0:
        delta_ta = dta

Mb = 6e10
Mtot_g = (1 + Q) * Mb
ab = 1.5 * kpc * (Mb / 1e10) ** 0.3


def hern(r, M, a):
    return M * r * r / (r + a) ** 2


rstar = math.sqrt(G * Mb * Msun / a0)
RQ_pt = rstar / (-math.log(q))
lo, hi = 1 * kpc, 1000 * kpc
for _ in range(80):
    mid = 0.5 * (lo + hi)
    mb = hern(mid, Mb * Msun, ab)
    f = (nu(G * mb / mid**2 / a0) - 1) * mb - Q * Mb * Msun
    lo, hi = (mid, hi) if f < 0 else (lo, mid)
RQ = 0.5 * (lo + hi)
P()
P(f"  M_b = {Mb:.0e} M_sun (Hernquist a = {ab/kpc:.2f} kpc), Fluss Q M_b = {Q*Mb:.3e}, Top-Hat-Masse (1+Q) M_b = {Mtot_g:.3e} M_sun")
P(f"  Bilanz: Operatorstruktur M_ph(<R) = (nu-1) M_b(<R) gleich Flussmasse Q M_b bei nu = 1+Q, d.h. q_lok = q_kos = {q:.5f}:")
P(f"   Punktmasse: R_Q = r*/(-ln q) = {rstar/kpc:.2f} kpc / {-math.log(q):.4f} = {RQ_pt/kpc:.2f} kpc;  Hernquist: R_Q = {RQ/kpc:.2f} kpc,"
  f" g_N(R_Q) = {G*hern(RQ,Mb*Msun,ab)/RQ**2:.3e} m/s^2")
RLg = RL(Mtot_g * Msun)
P(f"  Lagrange-Radius (mitbewegt) R_L = {RLg/kpc/1000:.3f} Mpc")
P("   z_c    z_ta   t_ta[Gyr] t_c[Gyr]  r_ta[kpc] r_vir[kpc]  r_vir/R_Q   Omega_vir[1/Gyr]  I_pre   I(z=0)")
gal = {}
for zc in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0):
    r = collapse_at(zc)
    I0, Ov = handover_I(r, t0)
    rv = r["x_ta"] / 2 * RLg
    gal[zc] = (r, I0, Ov, rv)
    P(f"   {zc:4.1f}  {z_of_t(r['t_ta']):5.2f}   {r['t_ta']*tH_Gyr:6.3f}   {r['t_c']*tH_Gyr:6.3f}   {r['x_ta']*RLg/kpc:8.1f}  {rv/kpc:8.1f}"
      f"    {rv/RQ:6.3f}      {Ov/tH_Gyr:8.3f}        {r['Ipre']:.4f}  {I0:7.3f}")
# z_c mit r_vir = R_Q
zz = sorted(gal)
lr = [math.log(gal[z][3] / RQ) for z in zz]
zcQ = float(np.interp(0.0, lr[::-1], zz[::-1]))
rQ = collapse_at(zcQ)
IQ, OvQ = handover_I(rQ, t0)
rvQ = rQ["x_ta"] / 2 * RLg
gal["Q"] = (rQ, IQ, OvQ, rvQ)
P(f"  r_vir = R_Q fuer z_c = {zcQ:.3f}: r_ta = {rQ['x_ta']*RLg/kpc:.1f} kpc, r_vir = {rvQ/kpc:.2f} kpc, t_c = {rQ['t_c']*tH_Gyr:.3f} Gyr, Omega_vir = {OvQ/tH_Gyr:.3f}/Gyr, I(z=0) = {IQ:.3f}")

# theta(t), w(t) fuer z_c = zcQ
r = tophat(rQ["di"], tstop=rQ["t_c"], rec=True)
P("  Verlauf fuer diesen Kollaps (theta_lok = 3 R'/R):")
P("     t[Gyr]    z      R/R_ta   theta_lok/H0   theta_bg/H0     w")
hist = r["hist"]
tt = np.array([hh[0] for hh in hist])
for frac in (0.2, 0.5, 0.8, 0.95, 0.99, 1.0, 1.3, 1.7, 1.99):
    tq = frac * r["t_ta"]
    j = int(np.argmin(abs(tt - tq)))
    t_, x_, v_ = hist[j]
    th = 3 * v_ / x_
    P(f"    {t_*tH_Gyr:7.3f}  {z_of_t(t_):6.2f}   {x_/r['x_ta']:6.3f}    {th:10.4f}    {3*E(1/(1+z_of_t(t_))):9.4f}   {w_of_theta(th):.4f}")
P("  Nach dem Kollaps (virialisiert, stationaer): theta = 0 exakt (Zeitfeld entlang des Killing-Vektors) -> w = 0.")
# Breite des w-Uebergangs
ws = np.array([w_of_theta(3 * hh[2] / hh[1]) for hh in hist])
i9 = int(np.argmax(ws < 0.9)); i1 = int(np.argmax(ws < 0.1))
P(f"  w faellt von 0,9 auf 0,1 zwischen t = {tt[i9]*tH_Gyr:.3f} und {tt[i1]*tH_Gyr:.3f} Gyr (Umkehr bei {r['t_ta']*tH_Gyr:.3f} Gyr)")

# g/g_RAR mit Rest-Fluss W (Fluss verteilt wie isotherme Kugel bis r_vir, Operator mit Gewicht 1-W)


def ratio_rar(rk, W, rv):
    rr = rk * kpc
    mb = hern(rr, Mb * Msun, ab)
    gN = G * mb / rr**2
    n_ = nu(gN / a0)
    g = gN + (1 - W) * (n_ - 1) * gN + W * G * Q * Mb * Msun * min(rr / rv, 1.0) / rr**2
    return g / (n_ * gN)


coef = max(abs(ratio_rar(rk, 1.0, rvQ) - 1) for rk in (5, 10, 20, 30, 50))
Wmax_gal = 0.05 / coef
P(f"  g/g_RAR mit Rest-Fluss W (Fluss als isotherme Kugel bis r_vir = {rvQ/kpc:.1f} kpc):  |g/g_RAR - 1| <= 5 % bei 5-50 kpc verlangt W <= {Wmax_gal:.4f}")

# ==========================================================================================
P()
P(f"  Form: Operatorstruktur M_ph(<r) = (nu-1) M_b(<r) gegen frei gefallenen Fluss (isotherme Kugel, Q M_b r/r_vir, r_vir = R_Q); Q*(-ln q) = {Q*(-math.log(q)):.4f}")
P("   r[kpc]   r/r*    M_ph/M_b   M_P,frei/M_b   Verhaeltnis")
for rk in (5, 10, 20, 30, 40, 50.64):
    rr = rk * kpc
    mb = hern(rr, Mb * Msun, ab)
    mph = (nu(G * mb / rr**2 / a0) - 1) * mb / (Mb * Msun)
    mpf = Q * min(rr / rvQ, 1.0)
    P(f"   {rk:6.2f}   {rr/rstar:5.2f}    {mph:7.4f}    {mpf:7.4f}       {mph/mpf:6.3f}")
P(f"  Umzuverteilende Masse (Integral |M_ph - M_P,frei| dr / R_Q) / (Q M_b):")
rrs = np.linspace(0.01, 1.0, 400) * rvQ
dif = [abs((nu(G*hern(x,Mb*Msun,ab)/x**2/a0)-1)*hern(x,Mb*Msun,ab)/(Mb*Msun) - Q*x/rvQ) for x in rrs]
P(f"   {float(np.trapezoid(dif, rrs/rvQ))/Q:.4f}")

hdr("4  Galaxienhaufen (10^15 M_sun) und Bullet Cluster")
zl = 0.3
rho_cz = rho_c0 * (Om * (1 + zl) ** 3 + OL)
rc_g, Mg1, rcut = 150 * kpc, 1e14 * Msun, 3000 * kpc


def Mbeta(r, rc):
    x = min(r, rcut) / rc
    return x - math.atan(x)


def Mb_cl(r):
    return (Mg1 * Mbeta(r, rc_g) / Mbeta(1000 * kpc, rc_g) + hern(r, 1.5e12 * Msun, 20 * kpc) + hern(r, 1e13 * Msun, 300 * kpc))


def nfw(M200, cc, rho):
    r200 = (3 * M200 / (4 * math.pi * 200 * rho)) ** (1 / 3)
    rs_ = r200 / cc
    m = lambda x: math.log(1 + x) - x / (1 + x)
    return (lambda r: M200 * m(r / rs_) / m(cc)), r200, rs_


Mref, r200c, rsc = nfw(1e15 * Msun, 4, rho_cz)
Mbz = Mb_cl(r200c)
P(f"  Baryonen: Gas beta = 2/3, r_c = 150 kpc, 1e14 in 1 Mpc (bis 3 Mpc); BCG 1,5e12 (a = 20 kpc); Galaxien 1e13 (a = 300 kpc); z = {zl}")
P(f"  Linsenmasse vertreten durch NFW M200 = 1e15, c = 4: r200 = {r200c/kpc:.0f} kpc; Baryonen in r200: {Mbz/Msun:.3e} (Anteil {Mbz/(1e15*Msun):.3f}, kosmisch {1/(1+Q):.3f})")
P(f"  g_N(r200)/a0 = {G*Mbz/r200c**2/a0:.4f}  gegen y_Q = {math.log(q)**2:.4f}: die Zone endet, bevor q_lok = q_kos erreicht ist.")
P("  Fluss der Zone M_P(<r): zwei Normierungen, NFW-Form (stossfrei): (i) Q M_b(<r200), (ii) q M_ref (kosmischer Anteil der Linsenmasse)")
P("  Energiebilanz: der Operator nimmt M_ph = (nu-1) M_b auf, der Rest W_rest = 1 - M_ph/M_P bleibt freier Fluss; Gesamtmasse M_b + M_P")
P("   r[kpc]  g_N/a0    nu    M_ref/(nu M_b)  | (i) M_P/M_b  W_rest  (M_b+M_P)/M_ref | (ii) M_P/M_b  W_rest  (M_b+M_P)/M_ref | reiner Operator nu M_b/M_ref")
Wrest = {}
for rk in (100, 250, 500, 1000, 1500):
    rr = rk * kpc
    mb = Mb_cl(rr)
    n_ = nu(G * mb / rr**2 / a0)
    mp1 = Q * Mbz * Mref(rr) / (1e15 * Msun)
    mp2 = q * Mref(rr)
    w1 = 1 - (n_ - 1) * mb / mp1
    w2 = 1 - (n_ - 1) * mb / mp2
    Wrest[rk] = (w1, w2)
    P(f"   {rk:5d}   {G*mb/rr**2/a0:6.3f}  {n_:6.3f}     {Mref(rr)/(n_*mb):6.3f}      | Operatoranteil 1-W_rest (i)/(ii): {1-w1:.3f}/{1-w2:.3f} |   {mp1/mb:6.3f}   {w1:6.3f}      {(mb+mp1)/Mref(rr):6.3f}      |   {mp2/mb:6.3f}   {w2:6.3f}      {(mb+mp2)/Mref(rr):6.3f}      |   {n_*mb/Mref(rr):6.3f}")

RLc = RL(1e15 * Msun)
tobs_cl = t_of_z(zl)
P(f"  Top-Hat 1e15 M_sun (R_L = {RLc/kpc/1000:.2f} Mpc), beobachtet bei z = {zl}:")
P("   z_c    z_ta   t_c[Gyr]  r_vir[kpc]  Omega_vir[1/Gyr]  (t_obs - t_c)[Gyr]")
clu = {}
for zc in (0.4, 0.5, 0.7, 1.0):
    r = collapse_at(zc)
    I_, Ov = handover_I(r, tobs_cl)
    clu[zc] = (r, I_, Ov)
    P(f"   {zc:4.1f}  {z_of_t(r['t_ta']):5.2f}   {r['t_c']*tH_Gyr:6.3f}   {r['x_ta']/2*RLc/kpc:7.0f}       {Ov/tH_Gyr:7.3f}        {(tobs_cl-r['t_c'])*tH_Gyr:6.3f}")

# Kapazitaetsgrenze: Zone r200 mit kosmischem Baryonenanteil; Operator fasst den ganzen Fluss, wenn g_N(r200) <= y_Q a0
P("  Kapazitaet: Operator fasst den ganzen Fluss der Zone, wenn g_N(r200) <= y_Q a0 (kosmischer Baryonenanteil, M_b = M200/(1+Q)):")
for zz_ in (0.0, 0.3, 1.0):
    Hz = H0s * E(1 / (1 + zz_))
    r200t = math.log(q) ** 2 * a0 * (1 + Q) / (100 * Hz**2)
    M200t = 100 * Hz**2 * r200t**3 / G
    P(f"   z = {zz_}: Grenzmasse M200 = {M200t/Msun:.3e} M_sun (r200 = {r200t/kpc:.0f} kpc); darueber bleibt Rest-Fluss")
P()
# Bullet Cluster
P()
zb = 0.296
rho_cb = rho_c0 * (Om * (1 + zb) ** 3 + OL)


def nfw2d(R, M200, cc, rho):
    r200 = (3 * M200 / (4 * math.pi * 200 * rho)) ** (1 / 3)
    rs_ = r200 / cc
    m = math.log(1 + cc) - cc / (1 + cc)
    x = R / rs_
    F = math.acosh(1 / x) / math.sqrt(1 - x * x) if x < 1 else math.acos(1 / x) / math.sqrt(x * x - 1)
    return M200 / m * (math.log(x / 2) + F)


P("  Bullet Cluster (z = 0,296): Hauptkomponente M200 = 1,5e15 M_sun, c = 4 (Groessenordnung aus dem Gedaechtnis);")
P("  Linsenmasse in 100 kpc um den Haupt-BCG ~2,8e13 M_sun (Clowe et al. 2006, aus dem Gedaechtnis)")
MP2 = q * nfw2d(100 * kpc, 1.5e15 * Msun, 4, rho_cb) / Msun
for rk in (100, 250):
    P(f"   stossfreier Anteil W_rest bei {rk} kpc (Haufenmodell oben, (i)/(ii)): {Wrest[rk][0]:.3f} / {Wrest[rk][1]:.3f}")
P(f"   Fluss (projiziert, 100 kpc): {MP2:.2e} M_sun; davon stossfrei (W_rest(100 kpc), (ii)): {Wrest[100][1]*MP2:.2e} M_sun")
P(f"   operatorgebundener Anteil (folgt den Baryonen, im Stoss dem Gas mit Nachlauf): {(1-Wrest[100][1])*MP2:.2e} M_sun")

hdr("5  Relaxation in die Operatorform: Rate lambda = Omega_vir/N* (pro Radiant der Zonenbahn)")
Nhi = IQ / math.log(1 / Wmax_gal)
P(f"  Galaxie (z_c = {zcQ:.2f}, r_vir = R_Q): Restabweichung W(z=0) <= {Wmax_gal:.3f}  ->  N* <= {Nhi:.2f}")
Ob_ = []
for zc in (0.5, 1.0):
    r = collapse_at(zc)
    Ob_.append(math.sqrt(Om / (2 * (r["x_ta"] / 2) ** 3)))
Nlo = max(0.26 / tH_Gyr * o for o in Ob_)
P(f"  Bullet: Nachlauf 1/lambda = N*/Omega_vir >= 0,26 Gyr mit Omega_vir = {Ob_[0]/tH_Gyr:.3f} bzw. {Ob_[1]/tH_Gyr:.3f} /Gyr (z_c = 0,5 bzw. 1)  ->  N* >= {Nlo:.3f}")
P(f"  Fenster: {Nlo:.3f} <= N* <= {Nhi:.2f}  ({Nlo/(2*math.pi):.3f} bis {Nhi/(2*math.pi):.2f} Umlaeufe der Zone)")
Nmid = Nhi
P("  Uebergabezeit 1/lambda = N*/Omega_vir an den Fenstergrenzen:")
P(f"   Galaxie (z_c = {zcQ:.2f}): {Nlo/OvQ*tH_Gyr:.3f} bis {Nhi/OvQ*tH_Gyr:.3f} Gyr")
for zc in (0.5, 1.0):
    P(f"   Haufen z_c = {zc}: {Nlo/clu[zc][2]*tH_Gyr:.3f} bis {Nhi/clu[zc][2]*tH_Gyr:.3f} Gyr")
P(f"   Uebergabe im Bullet waehrend 150 Myr seit dem Stoss (obere Grenze N* = {Nlo:.3f}): {1-math.exp(-Ob_[0]*0.15/tH_Gyr/Nlo):.3f}")

r05 = gal[0.5][0]
rv05 = r05["x_ta"] / 2 * RLg
cf05 = max(abs(ratio_rar(rk, 1.0, rv05) - 1) for rk in (10, 20, 30, 50))
NRAR = handover_I(r05, t0)[0] / math.log(cf05 / 0.05)
P(f"  RAR-Grenze N* = {NRAR:.2f} = {NRAR/(2*math.pi):.3f} Umlaeufe: 1/lambda Galaxie (z_c = {zcQ:.2f}) = {NRAR/OvQ*tH_Gyr:.3f} Gyr, Haufen z_c = 0,5: {NRAR/clu[0.5][2]*tH_Gyr:.3f} Gyr, z_c = 1: {NRAR/clu[1.0][2]*tH_Gyr:.3f} Gyr")
hdr("6  Vorhersagen (N* = obere Fenstergrenze: groesste Restabweichung; rechts N* aus der RAR-Enge bei z_c = 0,5)")
P("  6a  Galaxien M_b = 6e10: Rest-Fluss W und Abweichung g/g_RAR - 1 bei 10 / 20 / 50 kpc, nach Kollaps-Rotverschiebung und Beobachtungszeit")
P("   z_c     z_obs   W        10 kpc     20 kpc     50 kpc")
for zc in (0.5, 1.0, 2.0, 3.0, zcQ):
    r = gal[zc][0] if zc in gal else collapse_at(zc)
    rv = r["x_ta"] / 2 * RLg
    for zo in (2.0, 1.0, 0.5, 0.0):
        if zo >= zc:
            continue
        I_, Ov = handover_I(r, t_of_z(zo))
        W = math.exp(-I_ / Nmid)
        W2 = math.exp(-I_ / NRAR)
        P(f"   {zc:4.2f}   {zo:4.1f}   {W:.4f}   " + "   ".join(f"{ratio_rar(rk, W, rv)-1:+.4f}" for rk in (10, 20, 50))
          + f"   | N* = {NRAR:.2f}: W = {W2:.2e}  " + "  ".join(f"{ratio_rar(rk, W2, rv)-1:+.1e}" for rk in (10, 20, 50)))

P("   RAR-Enge als Schranke: |g/g_RAR - 1| <= 5 % bei 10-50 kpc zur Zeit z = 0 verlangt")
for zc in (0.5, 1.0, 2.0):
    r = gal[zc][0]
    rv = r["x_ta"] / 2 * RLg
    cf = max(abs(ratio_rar(rk, 1.0, rv) - 1) for rk in (10, 20, 30, 50))
    I_, Ov = handover_I(r, t0)
    P(f"    z_c = {zc}: W <= {0.05/cf:.3f}  ->  N* <= {I_/math.log(cf/0.05):.2f}")
P()
# Einsatz der ersten Abweichung: 2-sigma-Bereiche erreichen die Umkehr (BBKS-Transfer, aus dem Gedaechtnis)
P()
P("  6b  Erste Abweichung im Materiespektrum: 2-sigma-Bereiche erreichen die Umkehr (theta_lok = 0, w -> 0)")
P(f"      Bedingung 2 sigma(R_L, z) = delta_ta,lin = {delta_ta:.4f} (aus dem Top-Hat oben); Transfer BBKS mit Sugiyama-Form (aus dem Gedaechtnis), n_s = 0,965, sigma8 = 0,811")
Gam = Om * h * math.exp(-Ob - math.sqrt(2 * h) * Ob / Om)
kk = np.logspace(-4, 3, 6000)            # h/Mpc
qq = kk / Gam
T = np.log(1 + 2.34 * qq) / (2.34 * qq) * (1 + 3.89 * qq + (16.1 * qq) ** 2 + (5.46 * qq) ** 3 + (6.71 * qq) ** 4) ** -0.25
Pk0 = kk**0.965 * T**2


def sig(R):
    x = kk * R
    W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(float(np.trapezoid(kk**2 * Pk0 * W**2, kk)) / (2 * math.pi**2))


norm = 0.811 / sig(8.0)
P("      z     R_L [Mpc/h]   M [M_sun/h]    k_1 = 1/R_L [h/Mpc]   Anteil der Masse jenseits Umkehr (PS, alle Skalen >= M)")
for z in (0.0, 0.5, 1.0, 2.0, 3.0, 5.0):
    g_ = Dlin(1 / (1 + z)) / Dl[-1]
    lo, hi = 0.01, 50.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if 2 * norm * sig(mid) * g_ > delta_ta:
            lo = mid
        else:
            hi = mid
    R1 = math.sqrt(lo * hi)
    M1 = 4 / 3 * math.pi * R1**3 * Om * 2.775e11
    P(f"     {z:3.1f}    {R1:8.3f}     {M1:10.3e}        {1/R1:8.3f}                 {math.erfc(2/math.sqrt(2)):.4f}")
P("      (k_1 ist der Lagrange-Massstab; im Euler-Raum liegt die Struktur um (Dichte)^(1/3) ~ 5-6 bei hoeherem k)")


# Uebergabe vollzogen (W = 1/e): Kollaps-Rotverschiebung z_c* mit I(z_c*, z) = N*, Massstab der 2-sigma-Peaks mit Kollaps bei z_c*
dc_lin = None
r3 = collapse_at(3.0)
dc_lin = r3["di"] * Dlin(math.exp(np.interp(r3["t_c"], tg, lnag))) / Dlin(0.002)
zcg = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
tab = []
for zc in zcg:
    r = collapse_at(zc)
    xv = r["x_ta"] / 2
    tab.append((zc, r["t_c"], r["Ipre"], math.sqrt(Om / (2 * xv**3))))
P()
P(f"  6c  Uebergabe vollzogen (W = 1/e) fuer 2-sigma-Peaks: delta_c,lin = {dc_lin:.4f}; I(z_c, t) = I_pre + Omega_vir (t - t_c) (massenunabhaengig)")
P("      z_obs   z_c*    R_L [Mpc/h]   M [M_sun/h]    k_h = 1/R_L [h/Mpc]")
for zo in (0.0, 0.5, 1.0, 2.0):
    to = t_of_z(zo)
    Is = [(zc, Ip + Ov * (to - tc)) for zc, tc, Ip, Ov in tab if tc < to]
    zs_ = [a_[0] for a_ in Is]; Iv = [a_[1] for a_ in Is]
    zcs = float(np.interp(Nmid, Iv, zs_)) if Iv[-1] > Nmid else float("nan")
    g_ = Dlin(1 / (1 + zcs)) / Dl[-1]
    lo, hi = 0.005, 50.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if 2 * norm * sig(mid) * g_ > dc_lin:
            lo = mid
        else:
            hi = mid
    R1 = math.sqrt(lo * hi)
    P(f"      {zo:4.1f}   {zcs:5.2f}    {R1:8.3f}     {4/3*math.pi*R1**3*Om*2.775e11:10.3e}        {1/R1:7.3f}")
P("      Zwischen k_1 und k_h: Bereiche nach der Umkehr mit Rest-Fluss 1/e < W < 1 (Haufen, Gruppen) -> Uebergangsregime")

# f sigma8 bei N* = obere Fenstergrenze
_, D3, f3 = growth(lambda a: eps_fun(a) / Nmid)
P()
P(f"  6d  f sigma8: relative Abweichung vom LCDM-Wert bei N* = {Nmid:.2f} (Hintergrund-Uebergabe als nicht klumpend, konservativ):")
P("      " + "  ".join(f"z={z:.1f}: {(float(np.interp(math.log(1/(1+z)),Ngr,f3))*float(np.interp(math.log(1/(1+z)),Ngr,D3)))/(float(np.interp(math.log(1/(1+z)),Ngr,fl))*float(np.interp(math.log(1/(1+z)),Ngr,Dl)))-1:+.1e}" for z in (0.0, 0.5, 1.0, 2.0)))

# ------------------------------------------------------------------
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schliessungsregel.txt")
with open(path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n".join(OUT))
