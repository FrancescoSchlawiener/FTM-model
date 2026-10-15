# -*- coding: utf-8 -*-
"""
FTM - Zeitgesetz der Schwelle a0(z), wenn a0 = c k T/hbar = kappa/2pi die Temperatur eines Horizonts ist.

Theta = 2 pi a0/(c H).  Vakuumteil: Omega_s = ks Theta^2, ks = pi^2/45.
Friedmann: E^2 = R(a)/(1 - ks Theta^2),  R = Om a^-3 + Or a^-4 + OLr  (OLr: Rest der Dunklen Energie, konstant,
aus E(z=0) = 1).  w_s = -1 - (2/3) dln a0/dln a.

Abschnitte der Ausgabe:
 0  Eingaben
 1  Selbstkonsistente Loesungen je Horizont
      1a Hubble-Radius c/H                       Theta = 1
      1b Ereignishorizont R_E = a int_t^inf c dt'/a   Theta = c/(H R_E);  Theta' = Theta (1-ks Theta^2)(Theta - 1 + eps_R)
      1c scheinbarer Horizont, Kodama-Hayward   Theta = 1 - eps/2;  v' = 4 (Theta(v) - Theta_R), v = -ln(1-ks Theta^2)
      1d de-Sitter-Endhorizont c/H_inf          Theta = H_inf/H
      1e Teilchenhorizont R_P = a int_0^t c dt'/a (zum Vergleich)   Theta' = Theta (1-ks Theta^2)(eps_R - 1 - Theta)
 2  Kennzahlen je Horizont: a0(0), a0(z)/a0(0), Extremum, w_s(z), w_DE, CPL-Fit, fruehe Epochen
 3  DESI-Bedarf, gerechnet mit dem eigenen Omega_s,0 jedes Horizonts
 4  Fruehe massereiche Scheiben: Zusatzanteil f = 1 - 1/nu = q = exp(-sqrt(g_N/a0))
 5  Ereignishorizont: Schar der Loesungen und Kreuzungssatz
 6  Kurzfassung

Nur numpy/scipy.  Aufruf: python zeitgesetz_a0.py  (schreibt zeitgesetz_a0.txt daneben)
"""
import io
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []


def P(s=""):
    OUT.append(str(s))
    print(s)


def hdr(s):
    P()
    P("=" * 100)
    P(s)
    P("=" * 100)


# ------------------------------------------------------------------ 0 Eingaben
G = 6.674e-11
c = 2.99792458e8
Msun = 1.989e30
Mpc = 3.0856776e22
kpc = Mpc / 1e3
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
tH = 1 / H0 / 3.15576e16            # Gyr
XMIN = math.log(1e-9)
XMAX = 8.0

hdr("0  Eingaben")
P(f"   H0 = 67,4 km/s/Mpc, Omega_b = {Ob:.5f}, Q = {Q}, Omega_m = {Om:.4f}, Omega_r = {Or:.2e}, Omega_DE = {ODE:.4f}")
P(f"   ks = pi^2/45 = {ks:.5f};  cH0/2pi = {a0_mod:.4e} m/s^2;  Messwert a0 = 1,2e-10 [G] verlangt Theta_0 = {a0_mess/a0_mod:.4f}")
P(f"   1/H0 = {tH:.3f} Gyr.  Kein Parameter angepasst: OLr folgt je Horizont aus E(0) = 1.")


def R_s(x, OLr):
    a = math.exp(x)
    return Om * a**-3 + Or * a**-4 + OLr


def epsR_s(x, OLr):
    a = math.exp(x)
    m, r = Om * a**-3, Or * a**-4
    return (3 * m + 4 * r) / (2 * (m + r + OLr))


def R_v(x, OLr):
    a = np.exp(x)
    return Om * a**-3 + Or * a**-4 + OLr


def E_L(x):
    return np.sqrt(R_v(x, ODE))


def Th_of_v(v):
    return np.sqrt(np.maximum(-np.expm1(-np.maximum(v, 0.0)), 0.0) / ks)


# ------------------------------------------------------------------ 1 Loesungen
hdr("1  Selbstkonsistente Loesungen (H haengt ueber Omega_s = ks Theta^2 von a0 ab)")

CASES = {}

# 1a Hubble-Radius
OLr_H = ODE - ks
CASES["Hubble-Radius"] = dict(OLr=OLr_H, Th=lambda x: np.ones_like(np.asarray(x, float)))
P(f"1a Hubble-Radius: Theta = 1, Omega_s = {ks:.4f} zu jeder Zeit, OLr = {OLr_H:.4f}.")
P(f"   E^2 = R/(1 - ks): wirkt wie G_kos = {1/(1-ks):.4f} G zu allen Zeiten; Form von H(z) wie LCDM mit Omega_m,eff = {Om/(1-ks):.4f}")


# 1b Ereignishorizont, Zweig mit de-Sitter-Ende (rueckwaerts aus der Zukunft, stabile Richtung)
def event_v_rhs(OLr):
    def rhs(x, y):
        v = max(y[0], 0.0)
        Th = math.sqrt(max(-math.expm1(-v), 0.0) / ks)
        return [2 * ks * Th * Th * (Th - 1 + epsR_s(x, OLr))]
    return rhs


def event_from_future(OLr):
    v_end = -math.log(1 - ks)
    return solve_ivp(event_v_rhs(OLr), [XMAX, XMIN], [v_end], method="LSODA",
                     rtol=1e-11, atol=1e-15, dense_output=True)


def closure_event(OLr):
    s = event_from_future(OLr)
    Th0 = float(Th_of_v(s.sol(0.0)[0]))
    return OLr - ODE + ks * Th0**2


OLr_E = brentq(closure_event, 0.02, ODE - 1e-9, xtol=1e-13)
solE = event_from_future(OLr_E)
CASES["Ereignishorizont"] = dict(OLr=OLr_E, Th=lambda x, s=solE: Th_of_v(s.sol(np.asarray(x, float))[0]))
Th0E = float(CASES["Ereignishorizont"]["Th"](0.0))
P()
P(f"1b Ereignishorizont (Loesung mit de-Sitter-Ende, Theta -> 1): Theta_0 = {Th0E:.4f}, OLr = {OLr_E:.4f}, Omega_s,0 = {ks*Th0E**2:.4f}")
# Probe ueber das Integral R_E = a int_a^inf da'/(a'^2 E)
xx = np.linspace(-3.0, XMAX, 200001)
ThE = CASES["Ereignishorizont"]["Th"](xx)
EE = np.sqrt(R_v(xx, OLr_E) / (1 - ks * ThE**2))
integ = np.exp(-xx) / EE
cum = np.concatenate([[0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(xx))])
chiE = (cum[-1] - cum) + math.exp(-XMAX) / EE[-1]
for z in (0.0, 1.0, 2.0):
    i = int(np.argmin(np.abs(xx - math.log(1 / (1 + z)))))
    P(f"   Probe z = {z}: Theta aus Differentialgleichung {ThE[i]:.5f}, aus c/(H a int da/(a^2 H)) {1/(EE[i]*math.exp(xx[i])*chiE[i]):.5f}")


# 1c Kodama-Hayward
def ThR_s(x, OLr):
    return 1 - epsR_s(x, OLr) / 2


def kodama_from_future(OLr):
    def rhs(x, y):
        v = max(y[0], 0.0)
        Th = math.sqrt(max(-math.expm1(-v), 0.0) / ks)
        return [4 * (Th - ThR_s(x, OLr))]
    v_end = -math.log(1 - ks * ThR_s(XMAX, OLr) ** 2)
    return solve_ivp(rhs, [XMAX, XMIN], [v_end], method="Radau", rtol=1e-10, atol=1e-18, dense_output=True)


def closure_kodama(OLr):
    s = kodama_from_future(OLr)
    Th0 = float(Th_of_v(s.sol(0.0)[0]))
    return OLr - ODE + ks * Th0**2


OLr_K = brentq(closure_kodama, 0.3, ODE - 1e-9, xtol=1e-12)
solK = kodama_from_future(OLr_K)
CASES["Kodama-Hayward"] = dict(OLr=OLr_K, Th=lambda x, s=solK: Th_of_v(s.sol(np.asarray(x, float))[0]))
Th0K = float(CASES["Kodama-Hayward"]["Th"](0.0))
P()
P(f"1c Kodama-Hayward, T = (hbar H/2pi k)(1 - eps/2) = hbar Ricci/(24 pi k H):")
P(f"   Vorwaerts instabil (Rate 2(1-ks Theta^2)/(ks Theta) je e-Faltung); genau eine Loesung endet in de Sitter (Theta -> 1).")
P(f"   Sie folgt rueckwaerts aus der Zukunft: Theta_0 = {Th0K:.4f}, OLr = {OLr_K:.4f}, Omega_s,0 = {ks*Th0K**2:.4f}")
xk = np.linspace(-1.2, 0.3, 30001)
ThK = CASES["Kodama-Hayward"]["Th"](xk)
EK2 = R_v(xk, OLr_K) / (1 - ks * ThK**2)
epsK = -0.5 * np.gradient(np.log(EK2), xk)
for z in (0.0, 1.0):
    i = int(np.argmin(np.abs(xk - math.log(1 / (1 + z)))))
    P(f"   Probe z = {z}: Theta = {ThK[i]:.5f}, 1 - eps/2 aus H(z) = {1 - epsK[i]/2:.5f}")

# 1d de-Sitter-Endhorizont
OLr_D = ODE * (1 - ks)
CASES["de-Sitter-Ende"] = dict(OLr=OLr_D, Th=lambda x: math.sqrt(ODE) / E_L(np.asarray(x, float)))
P()
P(f"1d de-Sitter-Endhorizont: H_inf^2 = H0^2 OLr/(1-ks) und E(0) = 1 geben H_inf = sqrt(Omega_DE) H0 = {math.sqrt(ODE):.4f} H0.")
P(f"   rho_s = ks rho_krit(H_inf) = konstant, H(z) exakt LCDM; OLr = {OLr_D:.4f}, Omega_s,0 = {ks*ODE:.4f}")


# 1e Teilchenhorizont (vorwaerts stabil, Start Theta = 1 im Strahlungszeitalter)
def particle_from_past(OLr):
    def rhs(x, y):
        Th = y[0]
        return [Th * (1 - ks * Th * Th) * (epsR_s(x, OLr) - 1 - Th)]
    return solve_ivp(rhs, [XMIN, 1.0], [1.0], method="LSODA", rtol=1e-11, atol=1e-14, dense_output=True)


def closure_particle(OLr):
    Th0 = float(particle_from_past(OLr).sol(0.0)[0])
    return OLr - ODE + ks * Th0**2


OLr_P = brentq(closure_particle, 0.3, ODE, xtol=1e-12)
solP = particle_from_past(OLr_P)
CASES["Teilchenhorizont"] = dict(OLr=OLr_P, Th=lambda x, s=solP: s.sol(np.asarray(x, float))[0])
Th0P = float(CASES["Teilchenhorizont"]["Th"](0.0))
P()
P(f"1e Teilchenhorizont: Theta = 1 im Strahlungszeitalter (R_P = c/H exakt), Theta_0 = {Th0P:.4f}, OLr = {OLr_P:.4f}")

# ------------------------------------------------------------------ 2 Kennzahlen
hdr("2  Kennzahlen je Horizont")

x = np.linspace(math.log(1 / 5.0), 0.0, 40001)
a = np.exp(x)
zgrid = 1 / a - 1
zs_tab = [0.25, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
mcpl = x >= math.log(1 / 3.0)
RES = {}


def at_z(arr, z):
    return float(np.interp(math.log(1 / (1 + z)), x, arr))


def neff_equiv(Th):
    return (1 / (1 - ks * Th**2) - 1) * (1 + 3.044 * 0.2271) / 0.2271


for name, C in CASES.items():
    OLr, Thf = C["OLr"], C["Th"]
    Th = Thf(x)
    E2 = R_v(x, OLr) / (1 - ks * Th**2)
    E = np.sqrt(E2)
    A0 = Th * E                                  # a0 / (cH0/2pi)
    s = np.gradient(np.log(A0), x)               # dln a0/dln a
    ws = -1 - 2 * s / 3
    rhoDE = OLr + ks * Th**2 * E2
    wDE = -1 - np.gradient(np.log(rhoDE), x) / 3
    A = np.vstack([np.ones(mcpl.sum()), 1 - a[mcpl]]).T
    w0f, waf = np.linalg.lstsq(A, wDE[mcpl], rcond=None)[0]
    eps0 = -0.5 * np.gradient(np.log(E2), x)[-1]
    # Extrema von a0 fuer 0 < z < 4
    ext = []
    sg = np.sign(s[:-50]) if np.max(np.abs(s)) > 1e-6 else np.zeros(1)
    for i in np.where(np.diff(sg) != 0)[0]:
        typ = "Maximum" if sg[i] > 0 else "Minimum"
        ext.append((typ, float(zgrid[i])))
    th1090 = float(Thf(math.log(1 / 1091)))
    th9 = float(Thf(math.log(1e-9 * 1.0001)))
    RES[name] = dict(Th0=float(A0[-1]), a00=A0[-1] * a0_mod, r={z: at_z(A0, z) / A0[-1] for z in zs_tab},
                     s0=float(s[-1]), ws={z: at_z(ws, z) for z in (0, 0.5, 1, 2)}, wDE0=float(wDE[-1]),
                     w0f=float(w0f), waf=float(waf), ext=ext, Os0=ks * float(Th[-1])**2, OLr=OLr,
                     Os1090=ks * th1090**2, neff=neff_equiv(th9), E1=at_z(E, 1.0) / float(E_L(math.log(0.5))),
                     eps0=float(eps0), A0=A0, E=E)
    R = RES[name]
    P()
    P(f"{name}:  Theta_0 = {R['Th0']:.4f},  a0(0) = {R['a00']:.3e} m/s^2 = {R['a00']/a0_mess:.3f} x 1,2e-10,  Omega_s,0 = {R['Os0']:.4f}, OLr = {OLr:.4f}, eps_0 = {R['eps0']:.3f}")
    P("   a0(z)/a0(0):  " + "  ".join(f"z={z}: {R['r'][z]:.3f}" for z in zs_tab))
    P("   Extrema von a0 (0<z<4): " + (", ".join(f"{t} bei z = {zz:.2f}" for t, zz in ext) if ext else "keines"))
    P(f"   dln a0/dln a (z=0) = {R['s0']:+.3f};  w_s(z=0/0,5/1/2) = " + "/".join(f"{R['ws'][z]:+.3f}" for z in (0, 0.5, 1, 2)))
    P(f"   gesamte Dunkle Energie (OLr + Vakuumteil): w(0) = {R['wDE0']:+.3f};  CPL-Fit 0<z<2: w0 = {w0f:+.3f}, w_a = {waf:+.3f}")
    P(f"   H(z=1)/H_LCDM(z=1) = {R['E1']:.4f};  Omega_s(z=1090) = {R['Os1090']:.2e};  Delta N_eff-Aequivalent bei z = 1e9: {R['neff']:.3f}")

# ------------------------------------------------------------------ 3 DESI
hdr("3  DESI-Bedarf [G: DESI DR2 + CMB + SN, w0/w_a aus dem Gedaechtnis], mit dem eigenen Omega_s,0 jedes Horizonts")
desi = [("Pantheon+", -0.838, -0.62), ("DESY5", -0.752, -0.86), ("Union3", -0.667, -1.09)]
P("   Traegt der Vakuumteil die ganze Zeitabhaengigkeit: rho_s(z) = rho_DE,CPL(z) - OLr, a0(z)/a0(0) = sqrt(rho_s(z)/Omega_s,0).")
P("   Das Maximum von a0 liegt dort, wo rho_DE,CPL maximal ist (w = -1): 1 - a = -(1+w0)/w_a, unabhaengig von Omega_s,0.")
for nm, w0, wa in desi:
    ac = 1 + (1 + w0) / wa
    P(f"   {nm:10s}: Maximum bei z = {1/ac-1:.3f}")
zz = np.linspace(0, 4, 40001)
aa = 1 / (1 + zz)
DREQ = {}
for name, R in RES.items():
    P()
    P(f"   {name} (Omega_s,0 = {R['Os0']:.4f}); eigene Werte: a0(0,4)/a0(0) = {R['r'][0.4]:.3f}, a0(1)/a0(0) = {R['r'][1.0]:.3f}, a0(2)/a0(0) = {R['r'][2.0]:.3f}, s0 = {R['s0']:+.3f}")
    for nm, w0, wa in desi:
        rde = aa ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - aa))
        r2 = 1 + (rde - 1) * ODE / R["Os0"]
        vals = []
        for z in (0.4, 1.0, 2.0):
            v = float(np.interp(z, zz, r2))
            vals.append(math.copysign(math.sqrt(abs(v)), v))
        zero = float(zz[np.argmax(r2 < 0)]) if np.any(r2 < 0) else float("nan")
        s0req = -1.5 * (1 + w0) * ODE / R["Os0"]
        DREQ[(name, nm)] = (vals, zero, s0req)
        P(f"     verlangt ({nm:9s}): a0(0,4)/a0(0) = {vals[0]:+.3f}, a0(1)/a0(0) = {vals[1]:+.3f}, a0(2)/a0(0) = {vals[2]:+.3f}, "
          f"a0^2 = 0 bei z = {zero:.2f}, s0 = {s0req:+.3f}")
P()
P("   Mittlerer Exponent p in a0 propto a^p (p > 0: a0 waechst mit der Zeit) in den Intervallen z = 2->1, 1->0,4, 0,4->0:")
P("   (DESI mit Omega_s,0 = ks, damit a0^2 > 0 bis z = 2)")


def pexp(r_lo, r_hi, z_lo, z_hi):
    return math.log(r_lo / r_hi) / math.log((1 + z_hi) / (1 + z_lo))


for nm, w0, wa in desi:
    rde = aa ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - aa))
    r = np.sqrt(np.maximum(1 + (rde - 1) * ODE / ks, 1e-30))
    rz = {z: float(np.interp(z, zz, r)) for z in (0.0, 0.4, 1.0, 2.0)}
    P(f"     DESI {nm:13s}: p = {pexp(rz[1.0], rz[2.0], 1, 2):+.2f} / {pexp(rz[0.4], rz[1.0], 0.4, 1):+.2f} / {pexp(rz[0.0], rz[0.4], 0, 0.4):+.2f}")
for name, R in RES.items():
    rz = {0.0: 1.0, 0.4: R["r"][0.4], 1.0: R["r"][1.0], 2.0: R["r"][2.0]}
    P(f"     {name:18s}: p = {pexp(rz[1.0], rz[2.0], 1, 2):+.2f} / {pexp(rz[0.4], rz[1.0], 0.4, 1):+.2f} / {pexp(rz[0.0], rz[0.4], 0, 0.4):+.2f}")

# ------------------------------------------------------------------ 4 fruehe Scheiben
hdr("4  Fruehe massereiche Scheiben: Zusatzanteil f = 1 - 1/nu = exp(-sqrt(g_N/a0)), M_b = 1e11 M_sun (Punktmasse)")
Mb = 1e11 * Msun
for Rk in (5.0, 10.0):
    gN = G * Mb / (Rk * kpc) ** 2
    P(f"   R = {Rk:.0f} kpc (g_N = {gN:.3e} m/s^2):")
    P(f"     a0 = 1,2e-10 fest:    f = {math.exp(-math.sqrt(gN/a0_mess)):.3f} bei jedem z")
    for name, R in RES.items():
        fz = []
        for z in (0.0, 1.0, 2.0):
            a0z = R["a00"] * (1.0 if z == 0 else R["r"][z])
            fz.append(math.exp(-math.sqrt(gN / a0z)))
        P(f"     {name:18s}: f(z=0/1/2) = {fz[0]:.3f} / {fz[1]:.3f} / {fz[2]:.3f}")
P("   Gemessen [G]: massereiche sternbildende Scheiben bei z = 0,9-2,4 innerhalb des Halblichtradius f <~ 0,2,")
P("   mit abnehmendem Anteil zu hoeherem z (Genzel et al. 2017, 2020).")

# ------------------------------------------------------------------ 5 Ereignishorizont: Schar und Kreuzungssatz
hdr("5  Ereignishorizont: Schar der Loesungen und Kreuzungssatz")
P("   Theta = c/(H R_E) mit a chi = 1/(E Theta) und chi = int_a^inf da/(a^2 E) + C.  Die Differentialgleichung gilt fuer jedes C;")
P("   der Ereignishorizont verlangt C = 0.  Vorwaerts: Theta = 1 ist instabil (Rate (1-ks) je e-Faltung).")
P("   Theta < Theta*_0 laeuft auf 0 (dann C > 0, kein Ereignishorizont); Theta > Theta*_0 laeuft auf 1/sqrt(ks) = "
  f"{1/math.sqrt(ks):.3f} (H -> unendlich).")
P("   Kreuzungssatz: bei Theta = 1 gilt Theta' = (1 - ks) eps_R > 0.  Theta kreuzt 1 nur aufwaerts, d.h. a0 (dln a0/dln a = Theta - 1)")
P("   kann nur ein Minimum haben, nie ein Maximum, solange eps_R > 0 (Materie, Strahlung, konstanter Rest).")
for z in (0.0, 0.4, 1.0):
    P(f"     Theta'(Theta=1) bei z = {z}: {(1-ks)*epsR_s(math.log(1/(1+z)), OLr_E):+.4f}")


def event_family(Th0):
    OLr = ODE - ks * Th0**2
    v0 = -math.log(1 - ks * Th0**2)
    back = solve_ivp(event_v_rhs(OLr), [0.0, math.log(1 / 6.0)], [v0], method="LSODA", rtol=1e-11, atol=1e-15, dense_output=True)
    fwd = solve_ivp(event_v_rhs(OLr), [0.0, 60.0], [v0], method="LSODA", rtol=1e-10, atol=1e-13, dense_output=True)
    xf = np.linspace(0, 60, 60001)
    vf = fwd.sol(xf)[0]
    Ef = np.sqrt(R_v(xf, OLr) * np.exp(vf))
    trip = float(np.trapezoid(1 / Ef, xf)) * tH
    xb = np.linspace(math.log(1 / 6.0), 0.0, 20001)
    vb = back.sol(xb)[0]
    Thb = Th_of_v(vb)
    E2b = R_v(xb, OLr) * np.exp(vb)
    A0 = Thb * np.sqrt(E2b)
    s = np.gradient(np.log(A0), xb)
    rhoDE = OLr + ks * Thb**2 * E2b
    wDE = -1 - np.gradient(np.log(rhoDE), xb) / 3
    ab = np.exp(xb)
    m = xb >= math.log(1 / 3.0)
    w0f, waf = np.linalg.lstsq(np.vstack([np.ones(m.sum()), 1 - ab[m]]).T, wDE[m], rcond=None)[0]
    ext = []
    sg = np.sign(s[:-20])
    for i in np.where(np.diff(sg) != 0)[0]:
        ext.append(("Maximum" if sg[i] > 0 else "Minimum", float(1 / ab[i] - 1)))
    r1 = float(np.interp(math.log(0.5), xb, A0)) / A0[-1]
    Thend = float(Th_of_v(vf[-1]))
    return OLr, trip, Thend, float(wDE[-1]), w0f, waf, ext, r1, float(s[-1])


P()
for Th0 in (0.80, 0.85):
    OLr_t = ODE - ks * Th0**2
    ft = solve_ivp(event_v_rhs(OLr_t), [0.0, 40.0], [-math.log(1 - ks * Th0**2)], method="LSODA", rtol=1e-10, atol=1e-15)
    P(f"   Probe Theta_0 = {Th0} < Theta*_0: Theta(ln a = 40) = {float(Th_of_v(ft.y[0][-1])):.2e}  (laeuft auf 0, C > 0)")
P()
P("   Schar (Theta_0 >= Theta*_0; OLr = Omega_DE - ks Theta_0^2):")
for Th0, lab in ((Th0E, "de-Sitter-Ende"), (1.0, "a0(0) = cH0/2pi"), (a0_mess / a0_mod, "a0(0) = 1,2e-10"),
                 (math.sqrt(ODE / ks), "OLr = 0 (Vakuumteil = ganze DE)")):
    if lab == "de-Sitter-Ende":
        R = RES["Ereignishorizont"]
        P(f"     Theta_0 = {Th0:.4f} ({lab:30s}): OLr = {R['OLr']:.4f}, Ende de Sitter; w(0) = {R['wDE0']:+.3f}, "
          f"CPL w0/w_a = {R['w0f']:+.3f}/{R['waf']:+.3f}, a0(1)/a0(0) = {R['r'][1.0]:.3f}, Extrema: {R['ext'] or 'keines'}")
        continue
    OLr, trip, Thend, wd0, w0f, waf, ext, r1, s0f = event_family(Th0)
    P(f"     Theta_0 = {Th0:.4f} ({lab:30s}): OLr = {OLr:.4f}, Theta(ferne Zukunft) = {Thend:.3f}, H -> unendlich nach {trip:.1f} Gyr; "
      f"w(0) = {wd0:+.3f}, CPL w0/w_a = {w0f:+.3f}/{waf:+.3f}, a0(1)/a0(0) = {r1:.3f}, a0(0) = {Th0*a0_mod:.3e}, s0 = {s0f:+.3f}, Extrema: "
      + (", ".join(f"{t} z = {zz_:.2f}" for t, zz_ in ext) if ext else "keines"))

# ------------------------------------------------------------------ 6 Kurzfassung
hdr("6  Kurzfassung")
P("   Horizont            Theta_0  a0(0)/1,2e-10  a0(0,4)/a0(0)  a0(1)/a0(0)  a0(2)/a0(0)  Extremum        w_s(0)   CPL w0 / w_a     Om_s(1090)  dN_eff")
for name, R in RES.items():
    ex = ", ".join(f"{t[:3]}. z={zz_:.2f}" for t, zz_ in R["ext"]) if R["ext"] else "keines"
    P(f"   {name:18s}  {R['Th0']:.3f}    {R['a00']/a0_mess:.3f}          {R['r'][0.4]:.3f}          {R['r'][1.0]:.3f}        {R['r'][2.0]:.3f}"
      f"        {ex:14s}  {R['ws'][0]:+.3f}   {R['w0f']:+.3f} / {R['waf']:+.3f}   {R['Os1090']:.1e}     {R['neff']:.2f}")
P("   DESI [G] verlangt: w0 > -1, w_a < 0, Maximum von a0 bei z = 0,35-0,44.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeitgesetz_a0.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
