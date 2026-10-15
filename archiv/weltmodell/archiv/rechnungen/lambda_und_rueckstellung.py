# -*- coding: utf-8 -*-
# Λ als ganze innere Faltung der Zelle (a₀ fest seit dem Materiezeitalter) und Rückstellung als Bindungsrate dF/d ln a.
# Bausteine (Wachstum, σ(M), F nach PS mit δ_w und Sheth-Tormen) aus rueckstellung_zeitform.py übernommen.
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


OUT = []



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



(lambda *a: None)(f"   H0 = 67,4 km/s/Mpc, Omega_b = {Ob:.5f}, Q = {Q}, Omega_m = Omega_b(1+Q) = {Om:.4f}, Omega_r = {Or:.2e}, Omega_DE = {ODE:.4f}")
(lambda *a: None)(f"   ks = pi^2/45 = {ks:.5f}, cH0/2pi = {a0_mod:.4e} m/s^2, sigma8(LCDM) = {sig8} [G], n_s = {ns} [G]")
(lambda *a: None)("   Einzige Normierung: E(0) = 1, d.h. Omega_R0 = Omega_DE - Omega_s,0.  sigma(M) am fruehen Wachstum festgemacht (D = a bei a = 0,01).")

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



# ------------------------------------------------------------------ 1  Λ = ΣΔs_innen der Zelle
# ΣΔs = G E/(3c⁴) je mitbewegter Zelle, Materie (samt ξ) fest, Strahlung ∝ 1/a  ⇒  Λ(a) = Λ_m (1 + Ω_r/(Ω_m a)),  a₀ = c²/2Λ
# Wert: Λ_m = G ρ_m L³/(3c²)  ⇒  a₀ = 3c⁴/(2 G ρ_m L³),  ρ_m mitbewegt aus ω_m
rho_m = Om*3*H0**2/(8*math.pi*G)
def a0_of_L(Lgpc): L = Lgpc*1e3*Mpc; return 3*c**4/(2*G*rho_m*L**3)
def L_of_a0(a0): return (3*c**4/(2*G*rho_m*a0))**(1/3)/(1e3*Mpc)
# d_LSS = 2 × mitbewegter Abstand bis z = 1090 im ΛCDM-Hintergrund
xs = xg[xi(1090):i0+1]; dLSS = 2*c/H0*np.trapezoid(1/(ag[xi(1090):i0+1]*E_L[xi(1090):i0+1]), xs)/(1e3*Mpc)
print(f"Probe: a₀ = 3c⁴/(2Gρ_m L³) bei L = 27 Gpc: {a0_of_L(27):.4e} m/s² (Handwert 1,15·10⁻¹⁰ aus ΣΔs = 12,4 Gpc)")
print(f"d_LSS = {dLSS:.2f} Gpc;  Zelle L für a₀ = 1,04 / 1,10 / 1,19 / 1,2·10⁻¹⁰: " + " / ".join(f"{L_of_a0(v*1e-10):.2f}" for v in (1.04, 1.10, 1.19, 1.2)) + " Gpc")
print(f"L = 0,985·d_LSS (keine Kreispaare) ⇒ a₀ ≤ {a0_of_L(0.985*dLSS):.3e};  L = d_LSS ⇒ a₀ = {a0_of_L(dLSS):.3e} m/s²")
for z in (1e9, 1e4, 3400, 1090, 10, 1, 0):
    a = 1/(1+z); r = 1/(1 + Or/(Om*a)); Th = r*a0_of_L(dLSS)/(c*H0*np.interp(math.log(a), xg, E_L)/(2*math.pi))
    print(f"   z = {z:g}: a₀/a₀(heute) = {r:.4f},  Θ = 2πa₀/(cH) = {Th:.2e},  Ω_s = {ks*Th**2:.2e}")

# ------------------------------------------------------------------ 2  Rückstellung als Bindungsrate
# Vakuumteil fest (Θ0 = 1), ρ_R = Ω_R0 · (dF/d ln a)/(dF/d ln a)(0),  Ω_R0 = Ω_DE − ks;  selbstkonsistent mit dem Wachstum
def case_rate(meth, Mmin, Th0=1.0, it_max=80):
    E = E_L.copy()
    for it in range(it_max):
        D, f = growth(E)
        if Mmin == "MQ":
            M = M_Q(Th0/E, E)                                   # a₀ fest: Θ = Θ0/E
            F = F_of(sigma0(M)*D/D_L0, meth)
        else:
            F = F_of(sigma0(Mmin)*D/D_L0, meth)
        dF = np.gradient(F, xg); dF = np.maximum(dF, 0)
        rR = (ODE - ks*Th0**2)*dF/dF[i0]
        En = np.sqrt(Om*ag**-3 + Or*ag**-4 + ks*Th0**2 + rR)
        if np.max(np.abs(En/E - 1)) < 1e-9: break
        E = 0.5*(E + En)
    return dict(E=E, rde=ks*Th0**2 + rR, F=F, dF=dF)
def cpl(rde):
    m = (ag >= 1/3) & (ag <= 1.0); x = ag[m]; y = np.log(rde[m]/rde[i0]); best = None
    for w0 in np.arange(-1.6, -0.3, 0.005):
        A = -3*np.log(x) + 3*(x - 1); B = y + 3*(1 + w0)*np.log(x); wa = np.dot(A, B)/np.dot(A, A); r = np.sum((B - wa*A)**2)
        if best is None or r < best[0]: best = (r, w0, wa)
    return best[1], best[2]
def wz(rde, z): return -1 - np.gradient(np.log(rde), xg)[xi(z)]/3
print()
print("Rückstellung ∝ dF/d ln a, Vakuumteil fest | Verfahren | M_min | z des Maximums von ρ_DE | w_DE(0) | w_DE(1) | CPL w0 / wa | ρ_DE(z=2)/ρ_DE(0)")
for meth in ("PS-w", "ST"):
    for Mmin in (1e6, 1e8, 1e10, 1e11, 1e12, 1e13, "MQ"):
        C = case_rate(meth, Mmin); rde = C["rde"]
        m = (ag > 0.2) & (ag <= 1.0); jm = np.argmax(rde[m]); zmax = 1/ag[m][jm] - 1
        zs = f"{zmax:.2f}" if 0.01 < zmax < 3.9 else ("heute" if zmax <= 0.01 else "früh")
        w0, wa = cpl(rde); lab = "M_Q(z)" if Mmin == "MQ" else f"{Mmin:.0e}"
        print(f"   | {meth} | {lab} | {zs} | {wz(rde,0):+.3f} | {wz(rde,1):+.3f} | {w0:+.3f} / {wa:+.3f} | {rde[xi(2)]/rde[i0]:.3f}")
print("   Vergleich DESI DR2 + CMB + SN: w0 > −1, wa < 0; Pantheon+ −0,84…−0,86 / −0,58…−0,62; Maximum bei z ≈ 0,35–0,44 [G]")
