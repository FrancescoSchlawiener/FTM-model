# -*- coding: utf-8 -*-
"""
FTM - Rueckstellung traegt die Zeitform, a0 folgt einem Horizont.

Rueckstellbilanz sigma = SUM Ds_innen - DL_aussen = 0.  SUM Ds aller Materie ist in einem mitbewegten Gebiet fest (Ds ~ M);
wachsen kann nur der Teil, der nicht mehr mitatmet.  Das Modell hat dafuer das Gewicht der Schliessungsregel
w = th+^2/(th+^2 + th_c^2), th_c = a0/c: mitatmend w ~ 1, gebunden w -> 0.  Innen zaehlt also die gebundene Faltung,
F = gebundener Anteil der Faltung = gebundener Massenanteil (Ds ~ M).

Ansatz:  Omega_R(z) = Omega_R0 F(z)/F(0)   (Beitrag der Rueckstellung zu H^2/H0^2, kein Stoff)
Vakuumteil Omega_s = ks Theta^2, Theta = 2 pi a0/(cH), ks = pi^2/45; a0 nach Hubble-Radius, Ereignishorizont,
Kodama-Hayward, de-Sitter-Ende.   E^2 = R/(1 - ks Theta^2),  R = Om a^-3 + Or a^-4 + Omega_R0 F/F0.
Einzige Normierung: E(0) = 1.  F aus Press-Schechter (Schwelle aus dem Gewicht w, Kollaps 1,686) und Sheth-Tormen,
mit dem linearen Wachstum des Modells (= Standardmodell), Mindestmasse 1e6, 1e8, 1e10 M_sun und M_Q(z).
Wachstum, F und H werden zusammen iteriert.

Abschnitte:
 0 Eingaben
 1 Welche Faltung innen zaehlt; F(z) je Verfahren und Mindestmasse (LCDM-Hintergrund)
 2 Gesamtbild selbstkonsistent je Horizont, Verfahren, Mindestmasse
 3 Fruehe Grenzen
 4 Stelle des Scheiterns: verlangte gegen gelieferte heutige Steigung von a0
 5 Pruefbare Groessen H(z), D_M, f sigma8
 6 Kurzfassung
Werte mit [G] aus dem Gedaechtnis, nicht geprueft.   Aufruf: python rueckstellung_zeitform.py
"""
import io
import math
import os
import sys

import numpy as np

np.seterr(divide="ignore", invalid="ignore", over="ignore", under="ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import erfc, gammaincc
from scipy.special import gamma as Gam

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []


def P(s=""):
    OUT.append(str(s))
    print(s, flush=True)


def hdr(s):
    P()
    P("=" * 110)
    P(s)
    P("=" * 110)


# ------------------------------------------------------------------ 0 Eingaben
G = 6.674e-11
c = 2.99792458e8
Msun = 1.989e30
Mpc = 3.0856776e22
h = 0.674
H0 = h * 1e5 / Mpc
wb = 0.02237
Ob = wb / h**2
Q = 5.364
Om = Ob * (1 + Q)
Or = 9.14e-5
ODE = 1 - Om - Or
ks = math.pi**2 / 45
a0_mod = c * H0 / (2 * math.pi)
a0_mess = 1.2e-10
lq = -math.log(Q / (1 + Q))
sig8 = 0.811            # [G] Planck 2018, gilt fuer LCDM mit denselben Om, H0 (linear = Standardmodell)
ns = 0.965              # [G]
Tcmb = 2.7255
d_c = 1.686
A_ST, a_ST, p_ST = 0.3222, 0.707, 0.3   # Sheth-Tormen [G]
DESI = [("Pantheon+", -0.838, -0.62), ("DESY5", -0.752, -0.86), ("Union3", -0.667, -1.09)]   # [G]

dx = 0.0025
xg = np.arange(-8300, 3201) * dx          # ln a von -20,75 (z = 1e9) bis 8
ag = np.exp(xg)
i0 = 8300
assert abs(xg[i0]) < 1e-12
XMAX = xg[-1]
x_grow0 = math.log(1e-2)


def xi(z):
    return int(round((math.log(1 / (1 + z)) - xg[0]) / dx))


hdr("0  Eingaben")
P(f"   H0 = 67,4 km/s/Mpc, Omega_b = {Ob:.5f}, Q = {Q}, Omega_m = Omega_b(1+Q) = {Om:.4f}, Omega_r = {Or:.2e}, Omega_DE = {ODE:.4f}")
P(f"   ks = pi^2/45 = {ks:.5f}, cH0/2pi = {a0_mod:.4e} m/s^2, sigma8(LCDM) = {sig8} [G], n_s = {ns} [G]")
P("   Einzige Normierung: E(0) = 1, d.h. Omega_R0 = Omega_DE - Omega_s,0.  sigma(M) am fruehen Wachstum festgemacht (D = a bei a = 0,01).")

# ------------------------------------------------------------------ Hilfen: Wachstum, sigma(M), F
E_L = np.sqrt(Om * ag**-3 + Or * ag**-4 + ODE)


def growth(E):
    lnE = np.log(E)
    dl = np.gradient(lnE, xg)
    m = xg >= x_grow0 - 1e-12
    xs = xg[m]

    def rhs(x, y):
        Ex = np.interp(x, xg, E)
        d = np.interp(x, xg, dl)
        return [y[1], -(2 + d) * y[1] + 1.5 * Om * math.exp(-3 * x) / Ex**2 * y[0]]

    s = solve_ivp(rhs, (xs[0], xs[-1]), [1e-2, 1e-2], t_eval=xs, rtol=1e-9, atol=1e-15, method="LSODA")
    D = np.empty_like(xg)
    f = np.empty_like(xg)
    D[m] = s.y[0]
    f[m] = s.y[1] / s.y[0]
    D[~m] = ag[~m]
    f[~m] = 1.0
    return D, f


D_L, f_L = growth(E_L)
D_L0 = D_L[i0]

kk = np.logspace(-5, 5, 20001)             # 1/Mpc
lnk = np.log(kk)


def T_EH(k):                                # Eisenstein-Hu 1998 ohne Schwingungen [G]
    om, ob, fb, th = Om * h * h, Ob * h * h, Ob / Om, Tcmb / 2.7
    s = 44.5 * np.log(9.83 / om) / np.sqrt(1 + 10 * ob**0.75)
    aG = 1 - 0.328 * np.log(431 * om) * fb + 0.38 * np.log(22.3 * om) * fb**2
    Ge = Om * h * (aG + (1 - aG) / (1 + (0.43 * k * s) ** 4))
    q = k / h * th**2 / Ge
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q * q)


def T_BBKS(k):                              # BBKS mit Sugiyama-Form [G]
    Gm = Om * h * np.exp(-Ob - math.sqrt(2 * h) * Ob / Om)
    q = k / h / Gm
    return np.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25


def Wth(y):
    return np.where(y < 1e-3, 1 - y * y / 10, 3 * (np.sin(y) - y * np.cos(y)) / np.maximum(y, 1e-3) ** 3)


rhom0 = Om * 2.775e11 * h * h                 # M_sun/Mpc^3
lMg = np.linspace(math.log(1e2), math.log(1e18), 321)


def sigma_table(Tf):
    Pk = kk**ns * Tf(kk) ** 2

    def sg(R):
        return math.sqrt(np.trapezoid(kk**3 * Pk * Wth(kk * R) ** 2, lnk) / (2 * math.pi**2))

    nrm = sig8 / sg(8 / h)
    Rg = (3 * np.exp(lMg) / (4 * math.pi * rhom0)) ** (1 / 3)
    return np.log(np.array([nrm * sg(R) for R in Rg]))


LSIG = {"EH": sigma_table(T_EH), "BBKS": sigma_table(T_BBKS)}


def sigma0(M, tf="EH"):
    return np.exp(np.interp(np.log(M), lMg, LSIG[tf]))


# Schwelle aus dem Gewicht der Schliessungsregel: w = 1/2 bei th = th_c, th = 3 H_lok, th_c = a0/c = Theta H/2pi
# => H_lok/H = Theta/(6 pi).  Kugel im EdS: H_lok/H = 3 sin(e)(e - sin e)/(2(1 - cos e)^2), delta_lin = (3/5)(3/4)^(2/3)(e - sin e)^(2/3)
def delta_w(Th=1.0):
    e = brentq(lambda e: 3 * math.sin(e) * (e - math.sin(e)) / (2 * (1 - math.cos(e)) ** 2) - Th / (6 * math.pi), 2.0, math.pi)
    return 0.6 * 0.75 ** (2 / 3) * (e - math.sin(e)) ** (2 / 3)


d_w = delta_w(1.0)
d_ta = 0.6 * 0.75 ** (2 / 3) * math.pi ** (2 / 3)


def F_of(sz, meth):
    sz = np.maximum(sz, 1e-300)
    if meth == "PS-w":
        return erfc(d_w / (math.sqrt(2) * sz))
    if meth == "PS-c":
        return erfc(d_c / (math.sqrt(2) * sz))
    nu = d_c / sz
    u = a_ST * nu * nu / 2
    return A_ST * (erfc(np.sqrt(u)) + 2**-p_ST / math.sqrt(math.pi) * Gam(0.5 - p_ST) * gammaincc(0.5 - p_ST, u))


def M_Q(Th, E):
    """R_Q = r*/(-ln q) = R_200:  GM = (1+Q)^3 a0^3 lq^6/(1e4 H^4), M Gesamtmasse (Omega_m), a0 = Theta c H/2pi."""
    Hs = H0 * E
    a0 = Th * c * Hs / (2 * math.pi)
    return (1 + Q) ** 3 * a0**3 * lq**6 / (1e4 * Hs**4) / G / Msun


def F_curve(D, Th, E, meth, Mmin, tf="EH"):
    Gz = D / D_L0
    M = M_Q(Th, E) if Mmin == "MQ" else Mmin * np.ones_like(D)
    return F_of(sigma0(M, tf) * Gz, meth)


# ------------------------------------------------------------------ 1 F(z)
hdr("1  Welche Faltung innen zaehlt; gebundener Anteil F(z) (LCDM-Hintergrund, Theta = 1 fuer M_Q)")
P("   Mitbewegt ist SUM Ds = G M/(3c^2) fest.  Traegt die Rueckstellung alle Materie, bleibt DL_aussen je mitbewegtem Gebiet fest")
P("   und DL/R faellt wie 1/a: keine wachsende Rueckstellung.  Wachsen kann nur der Anteil, der nicht mehr mitatmet.")
P("   Kriterium des Modells: Gewicht der Schliessungsregel w = 1/2 bei th = th_c = a0/c, d.h. H_lok/H = Theta/(6 pi) = "
  f"{1/(6*math.pi):.4f} (Theta = 1).")
P(f"   Kugelmodell: lineare Schwelle delta_w = {d_w:.4f} (Umkehr {d_ta:.4f}, Kollaps {d_c}).  Theta = 0,78 / 0,86: "
  f"delta_w = {delta_w(0.78):.4f} / {delta_w(0.86):.4f}.")
P(f"   Hintergrund: 1 - w = 1/(1 + 36 pi^2/Theta^2) = {1/(1+36*math.pi**2):.2e} (Theta = 1), vernachlaessigbarer Boden.")
P("   Verfahren: PS-w = Press-Schechter mit delta_w (Kriterium des Modells), PS-c = PS mit Kollaps 1,686, ST = Sheth-Tormen (1,686).")
P("   Transferfunktion Eisenstein-Hu ohne Schwingungen [G]; zum Vergleich BBKS.  Massen = Gesamtmasse (Omega_m); baryonisch M/(1+Q).")
P()
for M in (1e6, 1e8, 1e10, 1e12):
    P(f"   sigma(M = {M:.0e} M_sun, z = 0): EH {sigma0(M):.3f}, BBKS {sigma0(M,'BBKS'):.3f}")
P(f"   M_Q heute (Theta = 1) = {M_Q(1.0, 1.0):.3e} M_sun, sigma(M_Q) = {float(sigma0(M_Q(1.0,1.0))):.3f}")
P()
MMIN = [1e6, 1e8, 1e10, "MQ"]
METH = ["PS-w", "PS-c", "ST"]
ZT = [0.4, 1.0, 2.0, 3.0]


def mlab(M):
    return "M_Q(z)" if M == "MQ" else f"{M:.0e}"


P("   Verfahren  M_min    F(0)    F(z)/F(0) z=0,4 / 1 / 2 / 3          dlnF/dlna(0)  w_R(0)   w_R(1)   w_R(0) je phys. Vol.  "
  "CPL allein w0 / w_a   Maximum von F'(z)")
ones = np.ones_like(xg)
mc = (xg >= math.log(1 / 3.0) - 1e-12) & (xg <= 1e-12)
for meth in METH:
    for M in MMIN:
        F = F_curve(D_L, ones, E_L, meth, M)
        dl = np.gradient(np.log(np.maximum(F, 1e-300)), xg)
        wR = -1 - dl / 3
        w0R, waR = np.linalg.lstsq(np.vstack([np.ones(mc.sum()), 1 - ag[mc]]).T, wR[mc], rcond=None)[0]
        dF = np.gradient(F, xg)
        jj = xi(20.0) + int(np.argmax(dF[xi(20.0):i0 + 1]))
        P(f"   {meth:5s}    {mlab(M):7s}  {F[i0]:.4f}  " + " / ".join(f"{F[xi(z)]/F[i0]:.3f}" for z in ZT)
          + f"     {dl[i0]:+.3f}        {wR[i0]:+.3f}   {wR[xi(1.0)]:+.3f}   {-dl[i0]/3:+.3f}              "
          f"{w0R:+.3f} / {waR:+.3f}      z = {1/ag[jj]-1:.2f}")
P()
P("   Transferfunktion (PS-w): F(0) und F(1)/F(0) mit EH / BBKS:")
for M in (1e6, 1e8, 1e10):
    Fe = F_curve(D_L, ones, E_L, "PS-w", M, "EH")
    Fb = F_curve(D_L, ones, E_L, "PS-w", M, "BBKS")
    P(f"     M_min = {M:.0e}: F(0) = {Fe[i0]:.4f} / {Fb[i0]:.4f};  F(1)/F(0) = {Fe[xi(1)]/Fe[i0]:.4f} / {Fb[xi(1)]/Fb[i0]:.4f}")
P("   F waechst monoton (gebunden bleibt gebunden, Preisach): w_R <= -1 zu jeder Zeit.  Je physikalischem Volumen gelesen")
P("   (F a^-3) verduennt die Rueckstellung wie Materie (w_R = -dlnF/dlna/3 nahe 0) und kann keine Dunkle Energie sein.")


# ------------------------------------------------------------------ 2 Gesamtbild
def eps_R_arr(OR0, Fr, dFr):
    m = Om * ag**-3
    r = Or * ag**-4
    R = m + r + OR0 * Fr
    return R, (3 * m + 4 * r - OR0 * dFr) / (2 * R)


def Th_of_v(v):
    return np.sqrt(np.maximum(-np.expm1(-np.maximum(v, 0.0)), 0.0) / ks)


def horizon_event(Fr, dFr, x_end):
    def sol_for(OR0, xe):
        R, eR = eps_R_arr(OR0, Fr, dFr)

        def rhs(x, y):
            v = max(y[0], 0.0)
            Th = math.sqrt(max(-math.expm1(-v), 0.0) / ks)
            return [2 * ks * Th * Th * (Th - 1 + np.interp(x, xg, eR))]

        return solve_ivp(rhs, [XMAX, xe], [-math.log(1 - ks)], method="LSODA", rtol=1e-10, atol=1e-14, dense_output=True)

    def clo(OR0):
        Th0 = float(Th_of_v(sol_for(OR0, 0.0).y[0][-1]))
        return OR0 - ODE + ks * Th0**2

    try:
        OR0 = brentq(clo, 0.01, ODE - 1e-9, xtol=1e-11)
    except ValueError:                        # keine Loesung mit de-Sitter-Ende und 0 < Omega_R0 < Omega_DE
        return None, None
    return OR0, sol_for(OR0, x_end)


def horizon_kodama(Fr, dFr, x_end):
    def sol_for(OR0, xe):
        R, eR = eps_R_arr(OR0, Fr, dFr)

        def rhs(x, y):
            v = max(y[0], 0.0)
            Th = math.sqrt(max(-math.expm1(-v), 0.0) / ks)
            return [4 * (Th - 1 + np.interp(x, xg, eR) / 2)]

        ThR = 1 - eR[-1] / 2
        return solve_ivp(rhs, [XMAX, xe], [-math.log(1 - ks * ThR**2)], method="LSODA", rtol=1e-10, atol=1e-16, dense_output=True)

    def clo(OR0):
        Th0 = float(Th_of_v(sol_for(OR0, 0.0).y[0][-1]))
        return OR0 - ODE + ks * Th0**2

    try:
        OR0 = brentq(clo, 0.01, ODE - 1e-9, xtol=1e-11)
    except ValueError:                        # keine Loesung mit de-Sitter-Ende und 0 < Omega_R0 < Omega_DE
        return None, None
    return OR0, sol_for(OR0, x_end)


def solve_case(hor, meth, Mmin, tf="EH", it_max=60):
    E = E_L.copy()
    Th = np.ones_like(xg)
    res = 1.0
    OR0 = ODE - ks
    F = F_curve(growth(E)[0], Th, E, meth, Mmin, tf)
    Fr = F / F[i0]
    for it in range(it_max):
        D, fg = growth(E)
        F = F_curve(D, Th, E, meth, Mmin, tf)
        Fr = F / F[i0]
        dFr = np.gradient(Fr, xg)
        if hor == "Hubble-Radius":
            OR0 = ODE - ks
            R, _ = eps_R_arr(OR0, Fr, dFr)
            Enew = np.sqrt(R / (1 - ks))
            Thn = np.ones_like(xg)
        elif hor == "de-Sitter-Ende":
            Finf = Fr[-1]
            OR0 = ODE / (1 + ks * Finf / (1 - ks))
            Hinf2 = OR0 * Finf / (1 - ks)
            R, _ = eps_R_arr(OR0, Fr, dFr)
            Enew = np.sqrt(R + ks * Hinf2)
            Thn = math.sqrt(Hinf2) / Enew
        else:
            x_end = xg[0] if res < 1e-5 else math.log(1e-2) - 0.5
            OR0n, sol = (horizon_event if hor == "Ereignishorizont" else horizon_kodama)(Fr, dFr, x_end)
            if OR0n is None:
                res = float("inf")
                break
            OR0 = OR0n
            mm = xg >= x_end - 1e-12
            v = np.zeros_like(xg)
            v[mm] = sol.sol(xg[mm])[0]
            Thn = Th_of_v(v)
            R, _ = eps_R_arr(OR0, Fr, dFr)
            Enew = np.sqrt(R / (1 - ks * Thn**2))
            if x_end > xg[0]:
                Thn[~mm] = Thn[mm][0]
                Enew[~mm] = np.sqrt(R[~mm] / (1 - ks * Thn[~mm] ** 2))
        mgr = xg >= x_grow0
        if not np.all(np.isfinite(Enew[mgr])):
            res = float("inf")
            break
        res_new = float(np.max(np.abs(Enew[mgr] / E[mgr] - 1)))
        lam = 1.0 if it < 6 else 0.5                 # Daempfung, falls die Rueckkopplung ueber M_Q schwingt
        E, Th = (1 - lam) * E + lam * Enew, (1 - lam) * Th + lam * Thn
        if res_new < 1e-7 and res < 1e-5:          # Rauschboden der Wachstumsrechnung ~1e-8
            res = res_new
            break
        res = res_new
    D, fg = growth(E)
    return dict(E=E, Th=Th, OR0=OR0, Fr=Fr, D=D, fg=fg, F0=F[i0], res=res, it=it + 1)


def analyse(C):
    E, Th, OR0, Fr, D, fg = C["E"], C["Th"], C["OR0"], C["Fr"], C["D"], C["fg"]
    rs = ks * Th**2 * E**2
    rR = OR0 * Fr
    rde = rs + rR
    dl = np.gradient(np.log(rde), xg)
    wde = -1 - dl / 3
    A0 = Th * E
    sA = np.gradient(np.log(A0), xg)
    m = (xg >= math.log(1 / 3.0) - 1e-12) & (xg <= 1e-12)
    w0, wa = np.linalg.lstsq(np.vstack([np.ones(m.sum()), 1 - ag[m]]).T, wde[m], rcond=None)[0]
    ext = []                                     # Vorzeichenwechsel von rho_DE' auf Schritten 0,05 in ln a, |dln rho/dln a| > 2e-4
    last = 0
    for j in range(xi(10.0), i0 + 1, 20):
        sv = 0 if abs(dl[j]) < 2e-4 else int(np.sign(dl[j]))
        if sv == 0:
            continue
        if last != 0 and sv != last:
            ext.append(("Max" if last > 0 else "Min", float(1 / ag[j] - 1)))
        last = sv
    dM = np.concatenate([[0], np.cumsum(0.5 * (np.exp(-xg[1:]) / E[1:] + np.exp(-xg[:-1]) / E[:-1]) * dx)])
    DM = dM[i0] - dM
    dML = np.concatenate([[0], np.cumsum(0.5 * (np.exp(-xg[1:]) / E_L[1:] + np.exp(-xg[:-1]) / E_L[:-1]) * dx)])
    DML = dML[i0] - dML
    fs8 = fg * D / D_L0 * sig8
    fs8L = f_L * D_L / D_L0 * sig8
    C.update(rs=rs, rR=rR, rde=rde, wde=wde, w0=float(w0), wa=float(wa), ext=ext, A0=A0, sA=sA, DM=DM, DML=DML,
             fs8=fs8, fs8L=fs8L, dRR0=float(OR0 * np.gradient(Fr, xg)[i0]), drs0=float(np.gradient(rs, xg)[i0]))
    return C


HOR = ["Hubble-Radius", "Ereignishorizont", "Kodama-Hayward", "de-Sitter-Ende"]
hdr("2  Gesamtbild selbstkonsistent: H(z) mit Materie, Strahlung, Vakuumteil (Horizont) und Rueckstellung Omega_R0 F/F0")
P("   w_DE = w der gesamten Dunklen Energie (Vakuumteil + Rueckstellung); CPL-Fit ueber 0 < z < 2 wie bei den Horizonten allein.")
P("   Extrema von rho_DE in 0 < z < 10 (Maximum = Kreuzung w = -1 von unten nach oben mit der Zeit).")
RES = {}
for hor in HOR:
    P()
    P(f"   {hor}")
    P("   Verf.  M_min    Theta_0  a0(0)/1,2e-10  a0(1)/a0(0)  Omega_s,0  Omega_R0   w_s(0)   w_R(0)   w_DE(0)  w_DE(1)  CPL w0 / w_a      Extrema rho_DE")
    for meth in METH:
        for M in MMIN:
            C = analyse(solve_case(hor, meth, M))
            RES[(hor, meth, M)] = C
            Th0 = C["Th"][i0]
            wR0 = -1 - np.gradient(np.log(np.maximum(C["Fr"], 1e-300)), xg)[i0] / 3
            ws0 = -1 - 2 * C["sA"][i0] / 3
            ex = ", ".join(f"{t} z={z:.2f}" for t, z in C["ext"]) if C["ext"] else "keines"
            if not C["res"] < 1e-5:
                ex = "KEIN FIXPUNKT"
            P(f"   {meth:5s}  {mlab(M):7s}  {Th0:.4f}   {Th0*a0_mod/a0_mess:.3f}          {C['A0'][xi(1)]/C['A0'][i0]:.3f}        "
              f"{ks*Th0**2:.4f}     {C['OR0']:.4f}    {ws0:+.3f}   {wR0:+.3f}   {C['wde'][i0]:+.3f}   {C['wde'][xi(1)]:+.3f}   "
              f"{C['w0']:+.3f} / {C['wa']:+.3f}   {ex}   (Iter. {C['it']}, Rest {C['res']:.0e})")

P()
P("   Transferfunktion BBKS statt EH (PS-w, M_min = 1e8), CPL w0 / w_a:")
for hor in HOR:
    C = analyse(solve_case(hor, "PS-w", 1e8, "BBKS"))
    P(f"     {hor:17s}: {C['w0']:+.3f} / {C['wa']:+.3f}   (EH: {RES[(hor,'PS-w',1e8)]['w0']:+.3f} / {RES[(hor,'PS-w',1e8)]['wa']:+.3f})")
P("   DESI DR2 + CMB + SN [G]: w0 / w_a = -0,838 / -0,62 (Pantheon+), -0,752 / -0,86 (DESY5), -0,667 / -1,09 (Union3);")
P("   Maximum der Dunklen Energie bei z = 0,35 / 0,41 / 0,44.")

# ------------------------------------------------------------------ 3 fruehe Grenzen
hdr("3  Fruehe Grenzen")
P("   Omega_s und Omega_R bei z = 1e9 (Strahlungszeitalter), 1090, 3, 2; Anteil der gesamten DE an H^2 bei z = 2 und 3 gegen LCDM")
OL = ODE / E_L**2


def neff_equiv(Os):
    return (1 / (1 - Os) - 1) * (1 + 3.044 * 0.2271) / 0.2271


P("   Horizont           Verf.  M_min    Om_s(1e9)  dN_eff   Om_s(1090)  Om_R(1090)  Om_R(3)   Om_R(2)   Om_DE(3)/LCDM  Om_DE(2)/LCDM")
for hor in HOR:
    for meth in ("PS-w", "ST"):
        for M in (1e6, 1e10):
            C = RES[(hor, meth, M)]
            Os = ks * C["Th"] ** 2
            OmR = C["rR"] / C["E"] ** 2
            OmDE = C["rde"] / C["E"] ** 2
            P(f"   {hor:17s}  {meth:5s}  {mlab(M):7s}  {Os[0]:.2e}   {neff_equiv(Os[0]):.2f}     {Os[xi(1090)]:.2e}    {OmR[xi(1090)]:.1e}     "
              f"{OmR[xi(3)]:.4f}    {OmR[xi(2)]:.4f}    {OmDE[xi(3)]/OL[xi(3)]:.3f}          {OmDE[xi(2)]/OL[xi(2)]:.3f}")
P(f"   LCDM: Omega_Lambda(3) = {OL[xi(3)]:.4f}, Omega_Lambda(2) = {OL[xi(2)]:.4f}")

# ------------------------------------------------------------------ 4 Stelle des Scheiterns
hdr("4  Stelle des Scheiterns: heutige Aenderung beider Teile und verlangte Steigung von a0")
P("   rho_DE' = rho_s' + rho_R'  (' = d/dln a, Einheiten rho_krit,0).  rho_R' >= 0 immer (F monoton).  Ein Maximum heute oder frueher")
P("   verlangt rho_s' < -rho_R' heute und rho_s' > -rho_R' bei z > z_max.")
P("   Verlangt (CPL [G]): rho_s' = -3(1+w0) Omega_DE - rho_R'  =>  s0 = dln a0/dln a = rho_s'/(2 rho_s,0);  geliefert: s0 des Horizonts.")
P("   Horizont           Verf.  M_min    rho_R'(0)  rho_s'(0)  s0 geliefert  s0 verlangt P+/DESY5/U3      rho_R'(1)  rho_s'(1)")
for hor in HOR:
    for meth in METH:
        for M in MMIN:
            C = RES[(hor, meth, M)]
            if not C["res"] < 1e-5:
                P(f"   {hor:17s}  {meth:5s}  {mlab(M):7s}  kein Fixpunkt")
                continue
            rs0 = C["rs"][i0]
            req =[(-3 * (1 + w0) * ODE - C["dRR0"]) / (2 * rs0) for _, w0, _ in DESI]
            dRR1 = float(C["OR0"] * np.gradient(C["Fr"], xg)[xi(1)])
            drs1 = float(np.gradient(C["rs"], xg)[xi(1)])
            P(f"   {hor:17s}  {meth:5s}  {mlab(M):7s}  {C['dRR0']:+.4f}    {C['drs0']:+.4f}    {C['sA'][i0]:+.3f}        "
              + " / ".join(f"{r:+.2f}" for r in req) + f"          {dRR1:+.4f}    {drs1:+.4f}")
P()
P("   Verlangte Form des Vakuumteils bei gegebener Rueckstellung (Theta_0 = 1, Omega_R0 = Omega_DE - ks, LCDM-Wachstum):")
P("   rho_s,verl(z) = rho_DE,CPL(z) - Omega_R0 F(z)/F0;  a0(z)/a0(0) = sqrt(rho_s,verl/ks)")
for meth in ("PS-w", "ST"):
    for M in (1e6, 1e10):
        F = F_curve(D_L, ones, E_L, meth, M)
        Fr = F / F[i0]
        row = []
        for nm, w0, wa in DESI:
            rdeC = ODE * ag ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - ag))
            rsv = rdeC - (ODE - ks) * Fr
            r = np.sqrt(np.maximum(rsv, 0) / ks)
            mm = (xg >= math.log(1 / 4)) & (xg <= 0)
            j = np.argmax(np.where(mm, r, -1))
            zneg = [1 / ag[k] - 1 for k in range(xi(4), i0) if rsv[k] < 0]
            row.append(f"{nm}: Max z={1/ag[j]-1:.2f} ({r[j]:.2f}), a0(1)={r[xi(1)]:.2f}, a0(2)={r[xi(2)]:.2f}"
                       + (f", a0^2<0 ab z={min(zneg):.2f}" if zneg else ""))
        P(f"   {meth:5s} {mlab(M):6s}: " + ";  ".join(row))
P("   Zum Vergleich geliefert a0(1)/a0(0): Hubble 1,95, Ereignishorizont 1,21, Kodama-Hayward 0,93, de Sitter 1 (ohne Rueckstellung).")

# ------------------------------------------------------------------ 5 pruefbare Groessen
hdr("5  Pruefbare Groessen: H(z)/H_LCDM, D_M(z)/D_M,LCDM, f sigma8(z)/LCDM (gleiches Om, H0; sigma am fruehen Wachstum)")
P("   Horizont           Verf.  M_min    H: z=0,5 / 1 / 2          D_M: z=0,5 / 1 / 2 / 1090             fs8: z=0,3 / 0,5 / 1        sigma8(0)")
for hor in HOR:
    for meth in ("PS-w", "ST"):
        for M in (1e6, 1e10):
            C = RES[(hor, meth, M)]
            Hs = " / ".join(f"{C['E'][xi(z)]/E_L[xi(z)]:.4f}" for z in (0.5, 1, 2))
            Ds = " / ".join(f"{C['DM'][xi(z)]/C['DML'][xi(z)]:.4f}" for z in (0.5, 1, 2, 1090))
            Fs = " / ".join(f"{C['fs8'][xi(z)]/C['fs8L'][xi(z)]:.4f}" for z in (0.3, 0.5, 1))
            P(f"   {hor:17s}  {meth:5s}  {mlab(M):7s}  {Hs}    {Ds}    {Fs}    {sig8*C['D'][i0]/D_L0:.4f}")
P(f"   LCDM f sigma8: z=0,3 / 0,5 / 1 = " + " / ".join(f"{f_L[xi(z)]*D_L[xi(z)]/D_L0*sig8:.4f}" for z in (0.3, 0.5, 1)))

# ------------------------------------------------------------------ 6 Kurzfassung
hdr("6  Kurzfassung")
ok = {k: C for k, C in RES.items() if C["res"] < 1e-5}
cnt = sum(1 for C in ok.values() if any(t == "Max" for t, _ in C["ext"]))
P(f"   Faelle mit Fixpunkt: {len(ok)} von {len(RES)};  davon mit Maximum der Dunklen Energie in 0 < z < 10: {cnt}")
for hor in HOR:
    for lab, Ms in (("M_min fest", (1e6, 1e8, 1e10)), ("M_Q(z)    ", ("MQ",))):
        ks_ = [(m, M) for m in METH for M in Ms if (hor, m, M) in ok]
        if not ks_:
            P(f"   {hor:17s} {lab}: kein Fixpunkt")
            continue
        w0s = [ok[(hor, m, M)]["w0"] for m, M in ks_]
        was = [ok[(hor, m, M)]["wa"] for m, M in ks_]
        wd0 = [ok[(hor, m, M)]["wde"][i0] for m, M in ks_]
        P(f"   {hor:17s} {lab}: w_DE(0) {min(wd0):+.3f} ... {max(wd0):+.3f},  CPL w0 {min(w0s):+.3f} ... {max(w0s):+.3f},  "
          f"w_a {min(was):+.3f} ... {max(was):+.3f}")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rueckstellung_zeitform.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
