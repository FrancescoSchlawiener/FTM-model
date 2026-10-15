# -*- coding: utf-8 -*-
"""
FTM - Atmen als Zustandsgleichung w(z) ohne gesetzte Massenskala (Kernfrage B5).

Abschnitte der Ausgabe:
 0  Eingaben
 1  Zustandsgleichung des Vakuumteils des Strukturterms: rho_s = (pi^2/45) Theta^2 rho_krit(z),
    w_s = -1 - (2/3) dln a0/dln a
 2  Bilanz des Strukturterms im Hintergrund:
    2a  nu-Teil / projizierter Fluss (Schliessungsregel): staubartig, Zusatzstrecken-Summe erhalten
    2b  Bindungsteil <Q_s> aus allen Trichtern (Halomodell, a0 = cH/2pi) und aus dem linearen Feld
 3  Rueckkopplung H <-> a0 <-> rho_s als Fixpunkt:
    3a  Kontraktion, Eindeutigkeit, Dimensionsargument (keine Massenskala)
    3b  Theta = 1 (de-Sitter-Temperatur)
    3c  Theta = 1 - eps/2 (Temperatur des dynamischen Horizonts), selbstkonsistent
    3d  a0 = c^2/L mit festem mitbewegtem L (Zelle), zum Vergleich
 4  Was die DESI-Zahlen von a0(z) verlangen (Inversion)
 5  Strukturabhaengige Rueckstellung: skalenfreie Wahl M = M_nl(z) und feste Skala M (Abtastung)
 6  Zusammenfassung

Nur numpy/scipy. Aufruf: python atmen_ohne_skala.py  (schreibt atmen_ohne_skala.txt daneben)
"""
import io
import math
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import erfc

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
h = 0.674
H0 = h * 1e5 / Mpc                      # 1/s
wb = 0.02237
Ob = wb / h**2
Q = 5.364                                # Omega_c/Omega_b = q/(1-q)
Om = Ob * (1 + Q)
Or = 9.14e-5                             # Photonen + Neutrinos (masselos)
ODE = 1 - Om - Or
Os = math.pi**2 / 45                     # Vakuumteil des Strukturterms, a0 = cH/2pi
a0_mod = c * H0 / (2 * math.pi)
a0_mess = 1.2e-10
fb = 1 / (1 + Q)
sig8, ns = 0.811, 0.965
dc = 1.686

hdr("0  Eingaben")
P(f"   H0 = {h*100:.1f} km/s/Mpc, Omega_b h^2 = {wb}, Q = {Q}, Omega_m = Omega_b(1+Q) = {Om:.4f}, Omega_r = {Or:.2e}")
P(f"   a0 = cH0/2pi = {a0_mod:.3e} m/s^2 (Modellwert, benutzt);  gemessen 1,2e-10 (nur Vergleich)")
P(f"   Omega_s = pi^2/45 = {Os:.4f}  (Anteil an Omega_DE = {ODE:.4f}: {Os/ODE:.4f})")
P(f"   Mit a0 = 1,2e-10 waere Omega_s = (4pi^4/15)(a0/c)^2/(3H0^2) = {4*math.pi**4/15*(a0_mess/c)**2/(3*H0**2):.4f}")
P(f"   sigma8 = {sig8}, n_s = {ns} (Planck 2018, aus dem Gedaechtnis), BBKS-Transfer mit Sugiyama-Form, delta_c = {dc}")

# ------------------------------------------------------------------ 1 Zustandsgleichung des Vakuumteils
hdr("1  Vakuumteil des Strukturterms als Zustandsgleichung")
P("   L_s enthaelt (2/l0^2) Q~(q), l0 = c^2/a0; bei y = 0 (q = 1) ist Q~ = -4pi^4/15.")
P("   Energiedichte: rho_s c^2 = (4pi^4/15) a0^2/(8 pi G).  Mit a0 = Theta * cH/2pi:")
P("      rho_s = (3H^2/8piG) * (pi^2/45) * Theta^2,   also Omega_s(z) = (pi^2/45) Theta(z)^2")
P("   Erhaltung (Bianchi, Materie getrennt erhalten): rho_s' + 3(1+w_s) rho_s = 0, rho_s propto a0^2:")
P("      w_s = -1 - (2/3) dln a0 / dln a")
chk = 4 * math.pi**4 / 15 * (a0_mod / c) ** 2 * c**2 / (8 * math.pi * G) / c**2 / (3 * H0**2 / (8 * math.pi * G))
P(f"   Probe: rho_s/rho_krit bei a0 = cH0/2pi = {chk:.6f} = pi^2/45 = {Os:.6f}")

# ------------------------------------------------------------------ Hilfsfunktionen: Hintergrund LCDM-artig, Wachstum, sigma(M)


def E_lcdm(a):
    return np.sqrt(Om * a**-3 + Or * a**-4 + ODE)


def growth():
    # D(a) aus der linearen Gleichung in ln a (LCDM-Hintergrund; der Vakuumteil aendert D nur schwach)
    def rhs(x, y):
        a = math.exp(x)
        E2 = Om * a**-3 + Or * a**-4 + ODE
        dlnE = (-3 * Om * a**-3 - 4 * Or * a**-4) / (2 * E2)
        Omz = Om * a**-3 / E2
        D, dD = y
        return [dD, -(2 + dlnE) * dD + 1.5 * Omz * D]
    x0 = math.log(1e-3)
    sol = solve_ivp(rhs, [x0, 0.0], [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-12)
    D1 = sol.sol(0.0)[0]
    return lambda a: sol.sol(np.log(a))[0] / D1, lambda a: sol.sol(np.log(a))[1] / sol.sol(np.log(a))[0]


Dg, fgrow = growth()

Gam = Om * h * math.exp(-Ob - math.sqrt(2 * h) * Ob / Om)
kk = np.logspace(-4, 3.5, 5000)          # h/Mpc
qq = kk / Gam
Tk = np.log(1 + 2.34 * qq) / (2.34 * qq) * (1 + 3.89 * qq + (16.1 * qq) ** 2 + (5.46 * qq) ** 3 + (6.71 * qq) ** 4) ** -0.25
Pk = kk**ns * Tk**2


def sigR(R):
    x = kk * R
    W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(float(np.trapezoid(kk**2 * Pk * W**2, kk)) / (2 * math.pi**2))


Anorm = sig8**2 / sigR(8.0) ** 2
Pk *= Anorm
rho_m_h = 2.775e11 * Om                 # M_sun/h pro (Mpc/h)^3, mitbewegt


def R_of_M(M):
    return (3 * M / (4 * math.pi * rho_m_h)) ** (1 / 3)


lnM_grid = np.linspace(math.log(1e5), math.log(1e17), 400)
sig_grid = np.array([sigR(R_of_M(math.exp(l))) for l in lnM_grid])


def sigM(M):
    return float(np.interp(math.log(M), lnM_grid, sig_grid))


def Mnl(z):
    s = dc / Dg(1 / (1 + z))
    # sigma faellt mit M: suche sigma(M) = s
    return math.exp(float(np.interp(-s, -sig_grid, lnM_grid)))


# ------------------------------------------------------------------ 2 Bilanz
hdr("2  Bilanz des Strukturterms im Hintergrund")
P("   2a  nu-Teil und projizierter Fluss (Schliessungsregel):")
P("       d/dt(rho_F) + 3H rho_F = -Gamma rho_F,  d/dt(rho_B) + 3H rho_B = +Gamma rho_F  =>  rho_b + rho_F + rho_B propto a^-3 exakt.")
P("       Die Faltungstiefe haengt nur an der Masse: Sum Delta s = G M_dyn/(3c^2). Innerhalb der Kapazitaetsgrenze R_Q gilt")
P("       M_dyn = (1+Q) M_b, gleich ob der Fluss frei faellt oder in Operatorform relaxiert ist. Also:")
P("       Sum Delta s pro mitbewegtem Volumen = G rho_m a^3 V_c /(3c^2) = const. Kein Faktor F_coll(>M), keine Massenskala.")
P("       Probe an der Grenze R_Q = r*/(-ln q): nu(R_Q) = 1/(1-q) = 1+Q, also Delta s_eff(R_Q) = (1+Q) G M_b/(3c^2), unabhaengig von a0:")
qk = Q / (1 + Q)
for M_b in (1e9, 6e10, 1e12, 1e14):
    for a0x in (a0_mod, 2 * a0_mod):
        rst = math.sqrt(G * M_b * Msun / a0x)
        RQ = rst / (-math.log(qk))
        nu = 1 / (1 - math.exp(-rst / RQ))
        P(f"       M_b = {M_b:8.1e} M_sun, a0 = {a0x:.2e}: r* = {rst/Mpc*1e3:8.2f} kpc, R_Q = {RQ/Mpc*1e3:8.2f} kpc, nu(R_Q) = {nu:.4f}")

P()
P("   2b  Bindungsteil: -(a0^2/8piG) <Q_s(y)>,  Q_s(y) = 4 int_0^s t^3/(e^t-1) dt,  s = sqrt(g/a0)")
P("       Verhaeltnis zum Vakuumteil: <Q_s>/(4pi^4/15). Jeder Trichter (Punktmasse M_b = f_b M) bis R_Q = r*/0,1709 bzw. bis R_200.")


def Qs(s):
    if s <= 0:
        return 0.0
    if s < 1e-3:
        return 4 * s**3 / 3
    s = min(s, 60.0)
    return 4 * quad(lambda t: t**3 / math.expm1(t), 0, s, limit=200)[0]


# J(u) = 4 pi int_0^u Q_s(1/x) x^2 dx  (r in Einheiten r*)
xg = np.logspace(-3, 3, 1500)
Qx = np.array([Qs(1 / x) for x in xg])
integrand = 4 * math.pi * Qx * xg**2
cumJ = np.concatenate([[4 * math.pi * (4 * math.pi**4 / 15) * xg[0] ** 3 / 3],
                       np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(xg))])
cumJ[1:] += cumJ[0]


def J(u):
    return float(np.interp(u, xg, cumJ))


P(f"       J(u) = int Q_s d^3r / r*^3:  J(1) = {J(1):.2f}, J(5,85) = {J(1/0.1709):.2f}, J(100) = {J(100):.2f}  (Schwanz (16pi/3) ln u)")
QV = 4 * math.pi**4 / 15


def ps_fnu(nu):
    return math.sqrt(2 / math.pi) * nu * math.exp(-nu**2 / 2)


def bind_fraction(z, rmode="RQ"):
    a = 1 / (1 + z)
    Hz = H0 * E_lcdm(a)
    a0z = c * Hz / (2 * math.pi)
    Dz = Dg(a)
    tot = 0.0
    lnMs = np.linspace(math.log(1e6), math.log(1e16), 300)
    dl = lnMs[1] - lnMs[0]
    for i, l in enumerate(lnMs):
        M = math.exp(l)                           # M_sun/h
        s = sigM(M) * Dz
        s2 = sigM(M * 1.01) * Dz
        dlns = (math.log(s2) - math.log(s)) / math.log(1.01)
        nu = dc / s
        dndlnM = rho_m_h / M * ps_fnu(nu) * abs(dlns)          # (h/Mpc)^3 mitbewegt
        n_phys = dndlnM * (1 + z) ** 3 * (h / Mpc) ** 3       # pro m^3 physikalisch
        Mb = fb * M / h * Msun
        rst = math.sqrt(G * Mb / a0z)
        if rmode == "RQ":
            u = 1 / 0.1709
        else:
            rho_c = 3 * Hz**2 / (8 * math.pi * G)
            R200 = (3 * (M / h * Msun) / (4 * math.pi * 200 * rho_c)) ** (1 / 3)
            u = R200 / rst
        tot += n_phys * rst**3 * J(u) * dl
    return tot / QV


P("       z     <Q_s>/(4pi^4/15) bis R_Q    bis R_200      M_nl(z) [M_sun/h]")
bf = {}
for z in (0.0, 0.5, 1.0, 2.0, 3.0):
    b1 = bind_fraction(z, "RQ")
    b2 = bind_fraction(z, "R200")
    bf[z] = b1
    P(f"      {z:3.1f}        {b1:10.3e}             {b2:10.3e}     {Mnl(z):10.3e}")

# lineares Feld
P()
P("       Lineares Feld zwischen den Trichtern: g_pec(k) = (3/2) Omega_m(z) H^2 a delta_k / k,  <Q_s> ~ (4/3)<(g/a0)^(3/2)>")
for z in (0.0, 1.0):
    a = 1 / (1 + z)
    Ez = E_lcdm(a)
    Hz = H0 * Ez
    Omz = Om * a**-3 / Ez**2
    D = Dg(a)
    intP = float(np.trapezoid(Pk * D**2, kk)) / (2 * math.pi**2)   # (Mpc/h)^2 * ... : P [(Mpc/h)^3] dk [h/Mpc] -> (Mpc/h)^2
    sig_disp = math.sqrt(intP) * Mpc / h * a                        # physikalische Laenge
    sg = 1.5 * Omz * Hz**2 * sig_disp
    a0z = c * Hz / (2 * math.pi)
    # Maxwell-Verteilung von |g| mit 3D-Streuung sg: <|g|^(3/2)> = sigma1^(3/2) * 2^(3/4) Gamma(9/4)/Gamma(3/2)
    s1 = sg / math.sqrt(3)
    mom = s1**1.5 * 2**0.75 * math.gamma(9 / 4) / math.gamma(1.5)
    P(f"       z = {z:.0f}: sigma_g = {sg:.2e} m/s^2 = {sg/a0z:.2e} a0(z),  <Q_s>/(4pi^4/15) = {(4/3)*mom/a0z**1.5/QV:.2e}")

# ------------------------------------------------------------------ 3 Rueckkopplung und Fixpunkt
hdr("3  Rueckkopplung H <-> a0 <-> rho_s als Fixpunkt")
P("   Friedmann mit dem Vakuumteil:  E^2 = R(a) + Omega_s Theta^2 E^2,   R(a) = Omega_m a^-3 + Omega_r a^-4 + Omega_Lr")
P("   Abbildung T(E^2) = R + Omega_s Theta^2 E^2 hat dT/dE^2 = Omega_s Theta^2 <= 0,219 < 1: Kontraktion,")
P("   also genau ein Fixpunkt  E^2 = R/(1 - Omega_s Theta^2), und er ist stabil (Iteration konvergiert mit Faktor 0,219 je Schritt).")
P("   Dimensionsargument: a0 = c k T/hbar mit T propto hbar H/k -> hbar und k fallen heraus. Aus G, c, H laesst sich")
P("   genau eine Masse bilden: M_H = c^3/(2GH). Heute:")
MH = c**3 / (2 * G * H0) / Msun
P(f"      M_H = {MH:.2e} M_sun,  sigma(M_H) = {sigM(MH*h):.2e},  F_coll(>M_H) = erfc(nu/sqrt2) = {erfc(dc/sigM(MH*h)/math.sqrt(2)):.1e}")
P("   Die Schwelle legt ebenfalls keine Masse fest: an der Umkehr (rho = 5,55 rho_m) gilt fuer die baryonische Neigung")
P("      g_N(R_ta)/a0 = (9pi^2/16)(Omega_b/2) H^2 R /(cH/2pi) = 17,4 Omega_b * H R/c, unabhaengig von der Masse bis auf R:")
for M in (1e10, 1e12, 1e14, 1e15):
    Rta = R_of_M(M) / (5.55) ** (1 / 3) / h * Mpc     # physikalisch bei z=0 (mitbewegt = physikalisch)
    ratio = 17.4 * Ob * H0 * Rta / c
    P(f"      M = {M:.0e} M_sun/h: R_ta = {Rta/Mpc:.3f} Mpc, g_N(R_ta)/a0 = {ratio:.2e}")
P("   -> jede gebundene Struktur liegt an ihrem Rand tief unter der Schwelle; die Schwelle schneidet keine Masse heraus.")


def w_from_rho(lna, lnrho):
    return -1 - np.gradient(lnrho, lna) / 3


def cpl_fit(a, w, amin=1 / 3.0):
    m = (a >= amin) & (a <= 1)
    A = np.vstack([np.ones(m.sum()), 1 - a[m]]).T
    sol, *_ = np.linalg.lstsq(A, w[m], rcond=None)
    return sol


def report(name, lna, rhoDE, extra=""):
    a = np.exp(lna)
    w = w_from_rho(lna, np.log(rhoDE))
    w0 = float(np.interp(0.0, lna, w))
    dw = np.gradient(w, a)
    wa_d = -float(np.interp(1.0, a, dw))
    w0f, waf = cpl_fit(a, w)
    zs = [0.0, 0.5, 1.0, 2.0]
    wz = [float(np.interp(math.log(1 / (1 + z)), lna, w)) for z in zs]
    P(f"   {name}: w(z=0/0,5/1/2) = " + "/".join(f"{x:+.3f}" for x in wz)
      + f";  w0 = {w0:+.3f}, w_a = {wa_d:+.3f} (Ableitung bei a=1);  CPL-Fit 0<z<2: w0 = {w0f:+.3f}, w_a = {waf:+.3f} {extra}")
    return w0f, waf


lna = np.linspace(math.log(1e-10), 0.0, 20001)
a = np.exp(lna)

# 3b Theta = 1
P()
P("   3b  Theta = 1 (T = hbar H/2pi k, a0 = cH/2pi zu jeder Zeit):")
OLr = 1 - Om - Or - Os
E2 = (Om * a**-3 + Or * a**-4 + OLr) / (1 - Os)
rhoDE = OLr + Os * E2
m = lna > math.log(1 / 4)
res_b = report("Theta=1", lna[m], rhoDE[m])
P(f"       Omega_s(z) = {Os:.4f} zu jeder Zeit. Friedmann: H^2 = (8piG/3) rho /(1-Omega_s): G_kos/G = {1/(1-Os):.4f}")
dNeff = (1 / (1 - Os) - 1) * (1 + 3.044 * 0.2271) / 0.2271
P(f"       Strahlungsaera: Delta N_eff-Aequivalent = {dNeff:.2f};  Rekombination: Anteil am Gesamtbudget {Os:.3f}")

# 3c Theta = 1 - eps/2, selbstkonsistent
P()
P("   3c  Theta = 1 - eps/2, eps = -dlnH/dln a (Temperatur des dynamischen Horizonts, T = (hbar H/2pi k)(1 - eps/2)).")
P("       Gleichbedeutend: Theta = (rho - 3p)/(4 rho) = Spur/(4 rho): Strahlung 0, Staub 1/4, Vakuum 1.")


lnaC = np.linspace(math.log(1e-10), 2.0, 20001)      # Gitter ueber heute hinaus, damit a = 1 innen liegt
aC = np.exp(lnaC)
iC = int(np.argmin(np.abs(lnaC)))


def theta_R(OLr_):
    Rr = Om * aC**-3 + Or * aC**-4 + OLr_
    epsR = (3 * Om * aC**-3 + 4 * Or * aC**-4) / (2 * Rr)
    return Rr, 1 - epsR / 2


def slow_solution(OLr_start=0.55, kmax=4):
    # Theta = Theta_R + Os Theta Theta' / (2 (1 - Os Theta^2)):  asymptotische Iteration (langsame Loesung)
    OLr_ = OLr_start
    for outer in range(80):
        Rr, ThR = theta_R(OLr_)
        Th_ = ThR.copy()
        hist = []
        for k in range(6):
            dTh = np.gradient(Th_, lnaC)
            Th_new = ThR + Os * Th_ * dTh / (2 * (1 - Os * Th_**2))
            hist.append(float(Th_new[iC]))
            if k < kmax:
                Th_ = Th_new
        OLr_new = 1 - Om - Or - Os * Th_[iC] ** 2
        if abs(OLr_new - OLr_) < 1e-12:
            break
        OLr_ = OLr_new
    E2_ = Rr / (1 - Os * Th_**2)
    return OLr_, E2_[:iC + 1], Th_[:iC + 1], hist


def family(Th0):
    # exakte Loesung der Differentialgleichung v' = 4(Theta(v) - Theta_R), rueckwaerts von heute (stabile Richtung)
    OLr_ = 1 - Om - Or - Os * Th0**2
    def rhs(x, y):
        aa = math.exp(x)
        Rr = Om * aa**-3 + Or * aa**-4 + OLr_
        ThR = 1 - (3 * Om * aa**-3 + 4 * Or * aa**-4) / (4 * Rr)
        v = max(y[0], 0.0)
        Thv = math.sqrt(max(-math.expm1(-v), 0.0) / Os)
        return [4 * (Thv - ThR)]
    v0 = -math.log(1 - Os * Th0**2)
    sol = solve_ivp(rhs, [0.0, math.log(1e-8)], [v0], dense_output=True, rtol=1e-10, atol=1e-14, method="LSODA")
    xx = np.linspace(math.log(1e-8), 0.0, 8001)
    v = sol.sol(xx)[0]
    aa = np.exp(xx)
    Rr = Om * aa**-3 + Or * aa**-4 + OLr_
    E2_ = Rr * np.exp(v)
    Th_ = np.sqrt(np.maximum(-np.expm1(-v), 0) / Os)
    return xx, OLr_, E2_, Th_


OLr2, E2c, Th, hist = slow_solution()
lna_save = lna
lna = lnaC[:iC + 1]
a = np.exp(lna)
rho_s = Os * Th**2 * E2c
rhoDE2 = OLr2 + rho_s
P("       Rueckkopplung jetzt differentiell: eps enthaelt dH/dt, also auch rho_s. Mit v = -ln(1 - Omega_s Theta^2):")
P("          v' = 4 (Theta(v) - Theta_R(a)),   Theta_R = 1 - eps_R/2 aus Materie, Strahlung, Omega_Lr allein.")
P("       Abweichungen von der langsamen Loesung wachsen vorwaerts mit der Rate dv'/dv = 2(1 - Omega_s Theta^2)/(Omega_s Theta) je e-Faltung.")
P("       Langsame Loesung: Theta = Theta_R + Omega_s Theta Theta'/(2(1 - Omega_s Theta^2)), iteriert; Theta_0 je Iteration: "
  + ", ".join(f"{x:.4f}" for x in hist[:6]) + "  (asymptotische Reihe; Streuung = Unsicherheit)")
P(f"       Fixpunkt: Theta_0 = {Th[-1]:.4f}, Omega_s,0 = {Os*Th[-1]**2:.4f}, Omega_Lr = {OLr2:.4f}")
for z in (0, 0.5, 1, 2, 5, 1090, 1e4, 1e9):
    ii = int(np.argmin(np.abs(lna - math.log(1 / (1 + z)))))
    rate = 2 * (1 - Os * Th[ii] ** 2) / (Os * max(Th[ii], 1e-30))
    P(f"       z = {z:>8g}: Theta = {Th[ii]:.4f}, Omega_s(z) = {Os*Th[ii]**2:.3e}, a0(z)/a0(0) = {Th[ii]*math.sqrt(E2c[ii])/Th[-1]:.3f}"
      f"  (Atmen a0 propto H: {math.sqrt(E2c[ii]):.3f});  Wachstumsrate der Abweichung {rate:.3g} je e-Faltung")
xx, OLf, E2f, Thf = family(float(Th[-1]))
P(f"       Probe: exakte Differentialgleichung rueckwaerts mit diesem Theta_0: Theta(z=1) = {np.interp(math.log(0.5), xx, Thf):.4f} "
  f"(langsam {np.interp(math.log(0.5), lna, Th):.4f}), Theta(z=1090) = {np.interp(math.log(1/1091), xx, Thf):.4f} (langsam {np.interp(math.log(1/1091), lna, Th):.4f})")
m = lna > math.log(1 / 4)
res_c = report("Theta=1-eps/2", lna[m], rhoDE2[m])
sel = lna > math.log(1 / 5)
ii = int(np.argmin(rho_s[sel])) + int(np.sum(~sel))
eps_all = -np.gradient(np.log(E2c), lna) / 2
iacc = int(np.argmin(np.abs(eps_all[sel] - 1))) + int(np.sum(~sel))
P(f"       rho_s hat ein Minimum bei z = {1/a[ii]-1:.3f}; Beginn der Beschleunigung (eps = 1) bei z = {1/a[iacc]-1:.3f}")
P(f"       a0 heute = Theta_0 cH0/2pi = {Th[-1]*a0_mod:.3e} m/s^2 (H0 = 67,4)")
P()
P("       Schar der exakten Loesungen, Parameter Theta_0 = 2 pi a0(0)/(c H0) (Integrationskonstante), rueckwaerts integriert:")
for Th0, lab in ((0.6, ""), (0.7, ""), (float(Th[-1]), "langsame Loesung"), (0.9, ""), (1.0, ""),
                 (1.2e-10 / (c * 70e3 / Mpc / (2 * math.pi)), "a0 = 1,2e-10, H0 = 70"), (a0_mess / a0_mod, "a0 = 1,2e-10, H0 = 67,4")):
    xx, OLf, E2f, Thf = family(Th0)
    rDE = OLf + Os * Thf**2 * E2f
    mm = xx > math.log(1 / 4)
    P(f"       Theta_0 = {Th0:.4f} {lab:24s} Omega_Lr = {OLf:.4f}, eps_0 = {2*(1-Th0):+.3f}:")
    report("         ", xx[mm], rDE[mm])


Th1090 = float(np.interp(math.log(1 / 1091), lna, Th))
Th_slow0 = float(Th[-1])
dlna0_c = float(np.gradient(np.log(Th * np.sqrt(E2c)), lna)[-1])
lna = lna_save
a = np.exp(lna)
m = lna > math.log(1 / 4)
# 3d Zelle
P()
P("   3d  a0 = c^2/L, L mitbewegt fest (a0 propto 1+z), Normierung a0(0) = cH0/2pi:")
rho_s3 = Os * a**-2
rhoDE3 = (ODE - Os) + rho_s3
res_d = report("Zelle", lna[m], rhoDE3[m], "(w_s = -1/3 exakt)")

# ------------------------------------------------------------------ 4 DESI-Inversion
hdr("4  Was die DESI-Zahlen von a0(z) verlangen (w0, w_a aus dem Gedaechtnis, DESI DR2 + CMB + SN)")
desi = [("Pantheon+", -0.838, -0.62), ("DESY5", -0.752, -0.86), ("Union3", -0.667, -1.09)]
P("   rho_s(z) = rho_DE,CPL(z) - (Omega_DE - Omega_s),  a0(z)^2/a0(0)^2 = rho_s(z)/Omega_s")
zlist = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
P("   Kombination      " + "  ".join(f"z={z:<4}" for z in zlist) + "   Nullstelle von a0^2   Maximum bei z")
for nm, w0, wa in desi:
    row = []
    zz = np.linspace(0, 5, 5001)
    aa = 1 / (1 + zz)
    rde = aa ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - aa))
    r2 = 1 + (rde - 1) * ODE / Os
    for z in zlist:
        v = float(np.interp(z, zz, r2))
        row.append(f"{math.copysign(math.sqrt(abs(v)), v):+.3f}")
    zero = zz[np.argmax(r2 < 0)] if np.any(r2 < 0) else float("nan")
    zmax = zz[np.argmax(r2)]
    P(f"   {nm:14s}   " + "  ".join(f"{x:6s}" for x in row) + f"      z = {zero:.2f}              {zmax:.2f}")
P("   (Eintraege: a0(z)/a0(0); negativ heisst a0^2 < 0, d.h. der Vakuumteil allein reicht nicht.)")
P("   Zum Vergleich Atmen a0 propto H:   " + "  ".join(f"{float(E_lcdm(1/(1+z))):6.3f}" for z in zlist))
P("   Mit a0 = 1,2e-10 (Omega_s = 0,292):")
Os_m = 4 * math.pi**4 / 15 * (a0_mess / c) ** 2 / (3 * H0**2)
for nm, w0, wa in desi:
    zz = np.linspace(0, 5, 5001)
    aa = 1 / (1 + zz)
    rde = aa ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - aa))
    r2 = 1 + (rde - 1) * ODE / Os_m
    row = [float(np.interp(z, zz, r2)) for z in zlist]
    zero = zz[np.argmax(r2 < 0)] if np.any(r2 < 0) else float("nan")
    P(f"   {nm:14s}   " + "  ".join(f"{math.copysign(math.sqrt(abs(v)), v):+6.3f}" for v in row) + f"      z = {zero:.2f}")
P()
P("   Benoetigte Steigung heute: dln a0/dln a = -3(1+w0) Omega_DE/(2 Omega_s):")
for nm, w0, wa in desi:
    P(f"   {nm:14s}: {-3*(1+w0)*ODE/(2*Os):+.3f}   (Atmen a0 propto H: {-1.5*Om:+.3f}; Zelle: -1; Theta=1-eps/2: {dlna0_c:+.3f})")

# ------------------------------------------------------------------ 5 strukturabhaengige Rueckstellung
hdr("5  Strukturabhaengige Rueckstellung rho_R = rho_Lr + Omega_s S(a)/S(1), Amplitude fest: beta = Omega_s/Omega_DE")
P(f"   beta = {Os/ODE:.4f} (nicht angepasst).  S1 = F_coll(>M,a) sqrt(H/H0)/a   (Delta s_eff propto sqrt(a0), R propto a)")
P("                                         S2 = F_coll(>M,a)/a            (Kapazitaetsgrenze R_Q: Delta s_eff unabhaengig von a0)")
lna5 = np.linspace(math.log(1 / 4), 0, 1201)
a5 = np.exp(lna5)
E5 = E_lcdm(a5)
D5 = np.array([Dg(x) for x in a5])


def run_S(Fc, label):
    out = []
    for nm, S in (("S1", Fc * np.sqrt(E5) / a5), ("S2", Fc / a5)):
        rho = (ODE - Os) + Os * S / S[-1]
        w = w_from_rho(lna5, np.log(rho))
        w0f, waf = cpl_fit(a5, w)
        im = int(np.argmax(S))
        zmax = 1 / a5[im] - 1 if 0 < im < len(S) - 1 else (99.0 if im == 0 else 0.0)
        out.append((nm, w0f, waf, zmax))
    P(f"   {label:34s} " + "   ".join(f"{nm}: w0 = {w0:+.3f}, w_a = {wa:+.3f}, z_max(S) = " + ("> 3 (faellt)" if zm > 90 else ("0 (steigt)" if zm == 0 else f"{zm:.2f}")) for nm, w0, wa, zm in out))
    return out


P("   a) skalenfrei, M = M_nl(z) (Kollapsbedingung sigma(M,z) = delta_c als Fixpunkt):")
Fnl = np.full_like(a5, erfc(1 / math.sqrt(2)))
P(f"      F_coll(>M_nl(z)) = erfc(1/sqrt2) = {erfc(1/math.sqrt(2)):.4f} fuer jedes z (Press-Schechter ist selbstaehnlich in nu = delta_c/sigma)")
run_S(Fnl, "M = M_nl(z):")
P("   b) skalenfrei, M = M_H(z) = c^3/(2GH): F_coll = 0, kein Beitrag.")
P("   b') M -> 0 (alle Masse traegt ihre Faltung; Schliessungsregel: der Fluss faltet von Anfang an mit):")
run_S(np.ones_like(a5), "F_coll = 1:")
P("   c) feste Masse M (Abtastung):")
resM = []
for lM in np.arange(10.0, 14.51, 0.25):
    M = 10**lM
    sM = sigM(M)
    Fc = erfc(dc / (sM * D5) / math.sqrt(2))
    o = run_S(Fc, f"M = 10^{lM:.2f} M_sun/h (nu0 = {dc/sM:.2f}):")
    resM.append((lM, dc / sM, o))
P()
P("   Naechste Punkte zu DESY5 (w0 = -0,752, w_a = -0,86; Fehler etwa 0,06 bzw. 0,22, aus dem Gedaechtnis), ohne Korrelation:")
for form in (0, 1):
    best = min(resM, key=lambda r: ((r[2][form][1] + 0.752) / 0.06) ** 2 + ((r[2][form][2] + 0.86) / 0.22) ** 2)
    d = math.sqrt(((best[2][form][1] + 0.752) / 0.06) ** 2 + ((best[2][form][2] + 0.86) / 0.22) ** 2)
    P(f"      {best[2][form][0]}: M = 10^{best[0]:.2f} M_sun/h, nu0 = {best[1]:.2f}, w0 = {best[2][form][1]:+.3f}, w_a = {best[2][form][2]:+.3f}, Abstand {d:.1f} (quadr. Summe der Einzelabstaende)")
    zeq = None
    for z in np.linspace(0, 3, 301):
        if Mnl(z) < 10 ** best[0]:
            zeq = z
            break
    P(f"          M_nl(z) = M bei z = {zeq:.2f}  (M_nl(0) = {Mnl(0):.2e} M_sun/h)")

# ------------------------------------------------------------------ 6 Zusammenfassung
hdr("6  Zusammenfassung")
P(f"   Bindungsteil <Q_s>/(4pi^4/15) heute: {bf[0.0]:.1e}  -> kein Hebel auf w.")
P(f"   Theta = 1       : CPL w0 = {res_b[0]:+.3f}, w_a = {res_b[1]:+.3f}; Delta N_eff-Aequivalent {dNeff:.2f}")
P(f"   Theta = 1-eps/2 : CPL w0 = {res_c[0]:+.3f}, w_a = {res_c[1]:+.3f}; Omega_s(z=1090) = {Os*Th1090**2:.4f}, BBN 0, Theta_0 = {Th_slow0:.4f}")
P(f"   Zelle a0 prop 1+z: CPL w0 = {res_d[0]:+.3f}, w_a = {res_d[1]:+.3f}")
P("   DESI verlangt a0(z)/a0(0) mit Maximum bei z ~ 0,3-0,5 und Abfall zu hohem z; Atmen (a0 propto H) steigt monoton.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "atmen_ohne_skala.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
