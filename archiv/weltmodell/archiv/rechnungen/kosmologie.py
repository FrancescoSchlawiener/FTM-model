# -*- coding: utf-8 -*-
"""
FTM - Fraktales Torus-Faltungs-Modell
Forschungsnotiz 08: Kosmologisches Leistungsspektrum (CMB, BAO), Wachstum, Torus-Topologie, Atmen.

Grundlagen des Modells (Notizen 01, 05; Dateien 03, 07):
    Struktur-Operator  nu(y) = 1/(1-q),  q = exp(-sqrt(y)),  y = a_N / a0
    Schwelle an das Atmen gekoppelt:  a0(z) = c * H(z) / (2 pi)

Alle Zahlen der Notiz 08 stammen aus diesem Skript.
Ausgabe: kosmologie.txt im selben Ordner.   Aufruf: python kosmologie.py
Nur numpy und math (kein scipy).
"""
import math
import os
import numpy as np

OUT = []


def P(s=""):
    OUT.append(s)


# ------------------------------------------------------------------ Konstanten
c = 2.99792458e8          # m/s
G = 6.674e-11             # m^3 kg^-1 s^-2
Mpc = 3.0856776e22        # m
Gpc = 1e3 * Mpc
LY = 9.4607e15            # m
ckm = c / 1e3             # km/s
arad = 7.5657e-16         # J m^-3 K^-4
Msun = 1.989e30

# Planck-2018-artige Eingabewerte (nur als Referenz und fuer die gemessene Geometrie)
h = 0.6736
wb = 0.02237              # omega_b = Omega_b h^2
wc = 0.1200               # omega_c (Staub-artige Komponente in LCDM)
Tcmb = 2.7255
Neff = 3.046
ns = 0.9649
sigma8 = 0.8111
zstar = 1089.92           # letzte Streuung
zdrag = 1059.94           # Ende der Baryonen-Mitnahme (drag epoch)

H0 = 100e3 * h / Mpc      # 1/s
rho_c100 = 3 * (100e3 / Mpc) ** 2 / (8 * math.pi * G)
wg = arad * Tcmb ** 4 / c ** 2 / rho_c100
wr = wg * (1 + 0.22711 * Neff)
Ob, Oc, Or = wb / h ** 2, wc / h ** 2, wr / h ** 2
OL = 1 - Ob - Oc - Or
wm = wb + wc
fb = wb / wm
a0_obs = 1.2e-10


def nu_s(y):
    """Struktur-Operator nu = 1/(1-exp(-sqrt(y)))."""
    return 1.0 / (1.0 - math.exp(-math.sqrt(y)))


def nu_v(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


def Lnu(y):
    """logarithmische Ableitung d ln nu / d ln y."""
    s = math.sqrt(y)
    e = math.exp(-s)
    return -s * e / (2 * (1 - e))


def a0_of_H(Hs):
    return c * Hs / (2 * math.pi)


def solve_y(nu_target):
    lo, hi = 1e-12, 1e3
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if nu_s(mid) > nu_target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# ------------------------------------------------------------------ Hintergrund-Modelle
class LCDM:
    """flach: Strahlung + Baryonen + staubartige Komponente + Rueckstellung (Lambda)."""
    def __init__(self, name, Om_dust, OLam):
        self.name, self.Od, self.OL = name, Om_dust, OLam

    def E2(self, z):
        z = np.asarray(z, dtype=float)
        return Or * (1 + z) ** 4 + (Ob + self.Od) * (1 + z) ** 3 + self.OL

    def rho_dust(self, z):   # Staub-Anteil in Einheiten rho_crit0
        return (Ob + self.Od) * (1 + np.asarray(z, float)) ** 3


class NuBackground:
    """Arbeitshypothese H-nu: der Operator wirkt im Hintergrund auf die Baryonen,
       rho_eff = nu(y_bg) * rho_b,   y_bg = pi * Omega_b(z) / X.
       Variante 'voll':  Omega_b(z) = rho_b / rho_crit(z) mit dem vollen H(z) (a0 ~ H, woertlich).
       Variante 'grav':  Omega_b(z) = rho_b / (rho_r + rho_eff) (nur gravitierende Dichte)."""
    def __init__(self, name, X, mode):
        self.name, self.X, self.mode = name, X, mode
        # heutiger Wert nu0 und Lambda aus Flachheit
        nu = 6.0
        for _ in range(200):
            if mode == "voll":
                y = math.pi * Ob / X
            else:
                y = math.pi * Ob / ((Or + nu * Ob) * X)
            nu = nu_s(y)
        self.nu0 = nu
        self.OL = 1 - Or - nu * Ob

    def nu(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        nu = np.full_like(z, 6.0)
        for _ in range(300):
            rb = Ob * (1 + z) ** 3
            if self.mode == "voll":
                E2 = Or * (1 + z) ** 4 + nu * rb + self.OL
            else:
                E2 = Or * (1 + z) ** 4 + nu * rb
            y = math.pi * rb / (E2 * self.X)
            nu = 0.5 * nu + 0.5 * nu_v(y)
        return nu

    def E2(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        return Or * (1 + z) ** 4 + self.nu(z) * Ob * (1 + z) ** 3 + self.OL

    def rho_dust(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        return self.nu(z) * Ob * (1 + z) ** 3


def Hkm(model, z):
    return 100 * h * np.sqrt(model.E2(z))       # km/s/Mpc


def Hsi(model, z):
    return H0 * np.sqrt(model.E2(z))


def rs(model, z, zmax=1e9, n=40000):
    x = np.linspace(math.log(1 + z), math.log(1 + zmax), n)
    zz = np.exp(x) - 1
    R = 3 * wb / (4 * wg) / (1 + zz)
    cs = ckm / np.sqrt(3 * (1 + R))
    integrand = cs / Hkm(model, zz) * (1 + zz)
    return np.trapezoid(integrand, x)          # Mpc (mitbewegt)


def DM(model, z, n=40000):
    if z <= 0:
        return 0.0
    x = np.linspace(0, math.log(1 + z), n)
    zz = np.exp(x) - 1
    return np.trapezoid(ckm / Hkm(model, zz) * (1 + zz), x)


def DV(model, z):
    return float((z * DM(model, z) ** 2 * ckm / np.ravel(Hkm(model, z))[0]) ** (1 / 3))


def zeq(model):
    lo, hi = 100.0, 1e5
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if float(np.ravel(model.rho_dust(mid))[0]) > Or * (1 + mid) ** 4:
            lo = mid
        else:
            hi = mid
    return mid


def H2H3(wb_, wm_, n_=ns):
    """Fitformeln Hu, Fukugita, Zaldarriaga & Tegmark 2001 (ApJ 549, 669):
       H2 = (dT_2/dT_1)^2, H3 = (dT_3/dT_1)^2."""
    H2 = 0.925 * wm_ ** 0.18 * 2.4 ** (n_ - 1) / (1 + (wb_ / 0.0164) ** (12 * wm_ ** 0.52)) ** 0.2
    H3 = 2.17 / (1 + (wb_ / 0.044) ** 2) * wm_ ** 0.59 * 3.6 ** (n_ - 1) / (1 + 1.63 * (1 - wb_ / 0.071) * wm_)
    return H2, H3


# Modelle
LC = LCDM("LCDM-aequivalent (Staub-Struktur + Lambda)", Oc, OL)
BO = LCDM("nur Baryonen + Lambda (ohne Staub-Term)", 0.0, 1 - Ob - Or)
nu_t = wm / wb
y_t = solve_y(nu_t)
X_cal = math.pi / (nu_t * y_t)
NV = NuBackground("H-nu 'voll' (a0 ~ H woertlich)", X_cal, "voll")
NG = NuBackground("H-nu 'grav' (nur gravitierende Dichte)", X_cal, "grav")

# ================================================================== 0 Eingaben
P("=" * 100)
P("FTM Notiz 08 - Kosmologie: CMB, BAO, Wachstum, Topologie, Atmen")
P("=" * 100)
P(f"h = {h}, omega_b = {wb}, omega_c(LCDM) = {wc}, T_CMB = {Tcmb} K, N_eff = {Neff}, n_s = {ns}, sigma8 = {sigma8}")
P(f"omega_gamma = {wg:.4e}, omega_r = {wr:.4e}")
P(f"Omega_b = {Ob:.4f}, Omega_c = {Oc:.4f}, Omega_r = {Or:.3e}, Omega_Lambda = {OL:.4f}, f_b = omega_b/omega_m = {fb:.4f}")
P(f"H0 = {100*h:.2f} km/s/Mpc = {H0:.4e} 1/s,  c/H0 = {ckm/(100*h):.1f} Mpc")
P(f"a0(z=0) = c H0/2pi = {a0_of_H(H0):.4e} m/s^2  (gemessen 1.2e-10; Verhaeltnis {a0_of_H(H0)/a0_obs:.3f})")

# ================================================================== 1 Schwelle frueh
P()
P("=" * 100)
P("1) H(z) und a0(z) = c H(z) / 2pi")
P("=" * 100)
P("Modell A: LCDM-aequivalenter Hintergrund (Strahlung + Baryonen + staubartige Komponente 0.2645 + Lambda)")
P("Modell B: nur Strahlung + Baryonen + Lambda-artige Rueckstellung (Omega_L = 1 - Omega_b - Omega_r)")
P(f"{'z':>7s} {'H_A[km/s/Mpc]':>14s} {'a0_A[m/s2]':>11s} {'a0_A/a0(0)':>10s} {'H_B[km/s/Mpc]':>14s} {'a0_B[m/s2]':>11s} {'a0_B/a0(0)':>10s}")
for z in [0, 0.5, 1, 2, 5, 10, 100, 1100]:
    HA, HB = float(Hkm(LC, z)), float(Hkm(BO, z))
    aA, aB = a0_of_H(float(Hsi(LC, z))), a0_of_H(float(Hsi(BO, z)))
    P(f"{z:7g} {HA:14.4g} {aA:11.3e} {aA/a0_of_H(H0):10.4g} {HB:14.4g} {aB:11.3e} {aB/a0_of_H(H0):10.4g}")
P("Die Schwelle lag bei der Rekombination um den Faktor ~2e4 (A) hoeher als heute.")

# ================================================================== 2 Stoerungen bei der Rekombination
P()
P("=" * 100)
P("2) Beschleunigung der Stoerungen bei z = 1100 im Vergleich mit a0(z)")
P("=" * 100)
z = 1100.0
Hz = float(Hsi(LC, z))
a0z = a0_of_H(Hz)
Obz = Ob * (1 + z) ** 3 / float(LC.E2(z))
Omz = (Ob + Oc) * (1 + z) ** 3 / float(LC.E2(z))
rho_crit_z = 3 * Hz ** 2 / (8 * math.pi * G)
P(f"H(1100) = {Hz:.4e} 1/s, a0(1100) = {a0z:.4e} m/s^2, c/H = {c/Hz/Mpc:.4f} Mpc (physikalisch)")
P(f"Omega_b(1100) = {Obz:.4f}, Omega_m(1100) [Modell A] = {Omz:.4f}, rho_b = {Obz*rho_crit_z:.3e} kg/m^3")
P("Formel (Kugel mit Radius R, Dichtekontrast delta):  g = (4 pi G/3) rho_x delta R = (1/2) Omega_x(z) H^2 R delta")
P("   y = g / a0(z) = pi * Omega_x(z) * delta * (R H / c)          [Befund: y ist durch R/R_Hubble und delta begrenzt]")
P("Potential-Form:  g ~ Phi / R  ->  y = 2 pi (Phi/c^2) (c/(R H))")
P()
P("a) Baryonen-Dichtekontrast delta_b (Quelle: nur Baryonen)")
P(f"{'R_mitb[Mpc]':>11s} {'R_phys[Mpc]':>11s} {'delta_b':>8s} {'g[m/s2]':>10s} {'y=g/a0':>10s} {'nu(y)':>8s} {'1/sqrt(y)':>9s} {'nu(1+L)':>8s}")
for Rc in [147.0, 50.0, 10.0, 1.0]:
    Rp = Rc / (1 + z) * Mpc
    for d in [1e-5, 1e-4, 1e-3]:
        g = 0.5 * Obz * Hz ** 2 * Rp * d
        y = g / a0z
        P(f"{Rc:11g} {Rp/Mpc:11.4g} {d:8.0e} {g:10.3e} {y:10.3e} {nu_s(y):8.1f} {1/math.sqrt(y):9.1f} {nu_s(y)*(1+Lnu(y)):8.1f}")
P()
P("b) Potential Phi/c^2 (gesamtes Potential, wie es der Sachs-Wolfe-Effekt misst: ~1e-5)")
P(f"{'R_mitb[Mpc]':>11s} {'Phi/c^2':>8s} {'g=Phi/R[m/s2]':>14s} {'y':>10s} {'nu(y)':>8s}")
for Rc in [147.0, 50.0, 10.0, 1.0]:
    Rp = Rc / (1 + z) * Mpc
    for ph in [1e-5, 1e-4, 1e-3]:
        g = ph * c ** 2 / Rp
        y = g / a0z
        P(f"{Rc:11g} {ph:8.0e} {g:14.3e} {y:10.3e} {nu_s(y):8.2f}")
P()
P("c) Kontrast, ab dem eine Stoerung die Schwelle erreicht (y = 1):  delta_1 = c / (pi Omega_x R H)")
for Rc in [147.0, 10.0, 1.0]:
    Rp = Rc / (1 + z) * Mpc
    P(f"   R = {Rc:6g} Mpc mitbewegt:  delta_1(Baryonen) = {c/(math.pi*Obz*Rp*Hz):.3e}   delta_1(gesamt, A) = {c/(math.pi*Omz*Rp*Hz):.3e}")
P()
P("d) Hubble-Kugel als Ganzes (F8, Notiz 05): g_H = G M_H / R_H^2 = c H / 2  ->  y_H = pi (jede Epoche, flach)")
P(f"   nu(pi) = {nu_s(math.pi):.4f},   d ln nu/d ln y bei y = pi: {Lnu(math.pi):.4f}")
P(f"   linearisierte Antwort parallel: nu(1+L) = {nu_s(math.pi)*(1+Lnu(math.pi)):.4f},  quer: {nu_s(math.pi):.4f},  Richtungsmittel nu(1+L/3) = {nu_s(math.pi)*(1+Lnu(math.pi)/3):.4f}")

# ================================================================== 3 Akustische Peaks
P()
P("=" * 100)
P("3) Akustische Peaks: Was muss der Operator fuer die CMB leisten?")
P("=" * 100)
Rstar = 3 * wb / (4 * wg) / (1 + zstar)
Rdrag = 3 * wb / (4 * wg) / (1 + zdrag)
P(f"Baryonen-Beladung R = 3 rho_b / (4 rho_gamma):  R(z*) = {Rstar:.4f},  R(z_drag) = {Rdrag:.4f}")
P("Einfaches Oszillator-Bild (konstantes Potential, Hu & Sugiyama): (Theta0+Psi) = (1/3+R) Psi cos(k r_s) - R Psi")
P(f"   Amplitude ungerade/gerade Peaks = 1 + 6R = {1+6*Rstar:.3f}  -> Leistungsverhaeltnis {(1+6*Rstar)**2:.2f}")
P("   (gilt nur, wenn Psi NICHT mit dem Fluid oszilliert -> genau dafuer braucht es die Staub-Rolle)")
P()
P("Fitformeln Hu et al. 2001 fuer die Peak-Hoehen (H2 = 2./1. Peak, H3 = 3./1. Peak, Temperatur^2):")
P(f"{'Fall':52s} {'omega_m':>8s} {'H2':>7s} {'H3':>7s} {'H3/H2':>7s}")
cases = [("LCDM (Planck-Werte)", wm),
         ("nur Baryonen (omega_m = omega_b), extrapoliert", wb),
         ("Staub-Rolle halb so stark", wb + wc / 2)]
for name, wmx in cases:
    a2, a3 = H2H3(wb, wmx)
    P(f"{name:52s} {wmx:8.4f} {a2:7.3f} {a3:7.3f} {a3/a2:7.3f}")
P("   Hinweis: ausserhalb von omega_m ~ 0.1-0.2 ist das eine Extrapolation der Fitformel (nur Tendenz).")
P()
P("Staub-Bedingung aus der Geometrie der Peaks:")
for m in [LC, BO]:
    rst, rd, dm = rs(m, zstar), rs(m, zdrag), DM(m, zstar)
    P(f"   {m.name:48s} z_eq = {zeq(m):7.0f}  r_s(z*) = {rst:7.2f} Mpc  r_d = {rd:7.2f} Mpc  D_M(z*) = {dm/1e3:6.3f} Gpc  100 theta* = {100*rst/dm:.4f}")
P("   Planck 2018 (Referenz): 100 theta* = 1.0411, r_d = 147.09 Mpc, z_eq = 3387  -> Modell A reproduziert die Referenz.")
P()
P("Arbeitshypothese H-nu: Operator wirkt im Hintergrund, rho_eff = nu(y_bg) rho_b, y_bg = pi Omega_b(z)/X")
P(f"   Ziel (Staub-Anteil wie Planck): nu_t = omega_m/omega_b = {nu_t:.4f}  ->  y_t = {y_t:.5f},  X = pi/(nu_t y_t) = {X_cal:.3f}")
P(f"   d.h. Bezugslaenge des Operators R_ref = R_Hubble / X = R_Hubble / {X_cal:.2f}  (a0 und R_ref beide an H gekoppelt)")
P(f"   Tiefes Regime: Fixpunkt nu -> X/pi = {X_cal/math.pi:.3f} (nur grobe Naeherung)")
P(f"   Variante 'voll': nu0 (heute) = {NV.nu0:.4f},  Omega_m,eff heute = {NV.nu0*Ob:.4f},  Omega_L = {NV.OL:.4f}")
P(f"   Variante 'grav': nu0 (heute) = {NG.nu0:.4f},  Omega_m,eff heute = {NG.nu0*Ob:.4f},  Omega_L = {NG.OL:.4f}")
P(f"   Zum Vergleich LCDM: Omega_m = {Ob+Oc:.4f}")
P()
P(f"{'z':>8s} {'nu voll':>9s} {'nu grav':>9s} {'w_S voll':>9s} {'w_S grav':>9s}")
for zz in [1e5, 1e4, 3400, 1100, 100, 10, 3, 1, 0.5, 0.0]:
    vals = []
    for m in [NV, NG]:
        e = 1e-3
        a1, a2_ = 1 / (1 + zz), 1 / (1 + zz) * (1 + e)
        z1, z2 = 1 / a1 - 1, 1 / a2_ - 1
        n1, n2 = float(m.nu(z1)[0]), float(m.nu(z2)[0])
        wS = -(1 / 3) * (math.log(n2 - 1) - math.log(n1 - 1)) / math.log(1 + e)
        vals += [n1, wS]
    P(f"{zz:8g} {vals[0]:9.4f} {vals[2]:9.4f} {vals[1]:9.4f} {vals[3]:9.4f}")
P("   w_S = -(1/3) d ln(nu-1)/d ln a  (w_S = 0 heisst: verhaelt sich wie Staub)")
P()
for m in [NV, NG]:
    rst, rd, dm = rs(m, zstar), rs(m, zdrag), DM(m, zstar)
    wmz = float(m.nu(zstar)[0]) * wb
    a2, a3 = H2H3(wb, wmz)
    P(f"   {m.name:40s} z_eq = {zeq(m):7.0f}  r_s(z*) = {rst:7.2f}  r_d = {rd:7.2f} Mpc  D_M(z*) = {dm/1e3:6.3f} Gpc  100 theta* = {100*rst/dm:.4f}  omega_m,eff(z*) = {wmz:.4f}  H2 = {a2:.3f} H3 = {a3:.3f}")

# ================================================================== 4 BAO und Wachstum
P()
P("=" * 100)
P("4) BAO-Skala und Wachstum f(k,z)")
P("=" * 100)
rd_LC = rs(LC, zdrag)
P("a) BAO-Abstaende D_V/r_d (DESI-Rotverschiebungen) und Abweichung der Hypothesen-Varianten von LCDM")
P(f"{'z_eff':>6s} {'LCDM D_V/r_d':>13s} {'voll':>9s} {'dev':>8s} {'grav':>9s} {'dev':>8s} {'nurBaryon':>10s} {'dev':>8s}")
rds = {m.name: rs(m, zdrag) for m in [NV, NG, BO]}
for zz in [0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330]:
    ref = DV(LC, zz) / rd_LC
    row = f"{zz:6.3f} {ref:13.4f}"
    for m in [NV, NG, BO]:
        v = DV(m, zz) / rds[m.name]
        row += f" {v:9.4f} {100*(v/ref-1):7.2f}%"
    P(row)

# --- lineares Leistungsspektrum (Eisenstein & Hu 1998, ohne Wiggles)
Om0 = Ob + Oc
theta2 = (Tcmb / 2.7) ** 2
s_EH = 44.5 * math.log(9.83 / wm) / math.sqrt(1 + 10 * wb ** 0.75)
alpha_G = 1 - 0.328 * math.log(431 * wm) * fb + 0.38 * math.log(22.3 * wm) * fb ** 2


def T_EH(k):   # k in h/Mpc
    Gam = Om0 * h * (alpha_G + (1 - alpha_G) / (1 + (0.43 * k * h * s_EH) ** 4))
    q = k * theta2 / Gam
    L0 = np.log(2 * math.e + 1.8 * q)
    C0 = 14.2 + 731 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q * q)


lnk = np.linspace(math.log(1e-5), math.log(1e3), 6000)
kk = np.exp(lnk)
Pk_shape = kk ** ns * T_EH(kk) ** 2


def Wth(x):
    x = np.asarray(x, float)
    out = np.where(x < 1e-3, 1 - x * x / 10, 3 * (np.sin(x) - x * np.cos(x)) / np.maximum(x, 1e-30) ** 3)
    return out


def sig_shape(R):
    return math.sqrt(np.trapezoid(kk ** 3 * Pk_shape / (2 * math.pi ** 2) * Wth(kk * R) ** 2, lnk))


Anorm = (sigma8 / sig_shape(8.0)) ** 2


def sigmaR(R):
    return math.sqrt(Anorm) * sig_shape(R)


def Delta_lin(k):   # dimensionslos, z = 0
    return math.sqrt(Anorm * k ** 3 * k ** ns * float(T_EH(k)) ** 2 / (2 * math.pi ** 2))


rho_m0 = Om0 * 3 * H0 ** 2 / (8 * math.pi * G)   # kg/m^3


def M_of_R(R_hMpc):
    R = R_hMpc / h * Mpc
    return 4 * math.pi / 3 * rho_m0 * R ** 3 / Msun   # M_sun


def R_of_M(M):   # M in M_sun -> R in Mpc/h
    R = (3 * M * Msun / (4 * math.pi * rho_m0)) ** (1 / 3)
    return R / Mpc * h


# --- Wachstum (ohne Strahlung, spaete Zeit), LCDM-Hintergrund
def E2_late(a):
    return Om0 * a ** -3 + (1 - Om0)


def Om_a(a):
    return Om0 * a ** -3 / E2_late(a)


def dlnH(a):
    return -1.5 * Om0 * a ** -3 / E2_late(a)


def growth(mu_func, a_ini=1e-3, n=6000):
    x = np.linspace(math.log(a_ini), 0.0, n)
    dx = x[1] - x[0]
    D, Dp = a_ini, a_ini
    Ds, fs = [D], [1.0]

    def rhs(xx, D, Dp):
        a = math.exp(xx)
        return Dp, -(2 + dlnH(a)) * Dp + 1.5 * Om_a(a) * mu_func(a) * D

    for i in range(n - 1):
        xx = x[i]
        k1 = rhs(xx, D, Dp)
        k2 = rhs(xx + dx / 2, D + dx / 2 * k1[0], Dp + dx / 2 * k1[1])
        k3 = rhs(xx + dx / 2, D + dx / 2 * k2[0], Dp + dx / 2 * k2[1])
        k4 = rhs(xx + dx, D + dx * k3[0], Dp + dx * k3[1])
        D += dx / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Dp += dx / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Ds.append(D)
        fs.append(Dp / D)
    return np.exp(x), np.array(Ds), np.array(fs)


aG, DG, fG = growth(lambda a: 1.0)
D0 = DG[-1]


def D_of_z(z):
    return float(np.interp(math.log(1 / (1 + z)), np.log(aG), DG)) / D0


def f_of_z(z, fs=fG):
    return float(np.interp(math.log(1 / (1 + z)), np.log(aG), fs))


P()
P(f"b) Referenz-Wachstum LCDM: sigma8 = {sigma8}, sigma(R=8 Mpc/h) Kontrolle = {sigmaR(8.0):.4f}")
P(f"{'z':>5s} {'D(z)':>7s} {'f(z)':>7s} {'Om(z)^0.55':>10s} {'f sigma8':>9s}")
for zz in [0.0, 0.3, 0.51, 0.71, 0.93, 1.32, 2.0]:
    P(f"{zz:5.2f} {D_of_z(zz):7.4f} {f_of_z(zz):7.4f} {Om_a(1/(1+zz))**0.55:10.4f} {f_of_z(zz)*sigma8*D_of_z(zz):9.4f}")

P()
P("c) Wirkt nu direkt auf die linearen Moden?  y(k,z) = pi Omega_m(z) Delta(k,z) (aH/(c k))")
P(f"{'k[h/Mpc]':>8s} {'z':>4s} {'Delta':>8s} {'y':>10s} {'nu(y)':>7s} {'Wachstumsexponent p(mu=nu)':>27s}")
for k in [0.01, 0.05, 0.1, 0.2, 0.5]:
    for zz in [0.0, 1.0, 2.0]:
        a = 1 / (1 + zz)
        Dl = Delta_lin(k) * D_of_z(zz)
        aH_c = a * 100 * math.sqrt(E2_late(a)) / ckm      # h/Mpc (H in h km/s/Mpc)
        y = math.pi * Om_a(a) * Dl * aH_c / k
        mu = nu_s(y)
        p = (-1 + math.sqrt(1 + 24 * mu)) / 4
        P(f"{k:8.2f} {zz:4.1f} {Dl:8.4f} {y:10.3e} {mu:7.2f} {p:27.3f}")
P("   p ist der Exponent D ~ a^p im Materie-Zeitalter bei konstantem mu (LCDM: p = 1).")

# --- Umschaltmodell: Operator wirkt nur in vom Atmen entkoppelten Regionen
delta_ta = 1.062          # linearer Kontrast bei Umkehr (EdS)
one_ta = 9 * math.pi ** 2 / 16   # 1 + Delta_ta = 5.55


def eps_switch(k, a):
    """eps(k,a) = F_dec * (f_b nu(y_ta) - 1)"""
    R = 1.0 / k                        # Mpc/h, Lagrange-Radius
    sig = sigmaR(R) * (float(np.interp(math.log(a), np.log(aG), DG)) / D0)
    F = math.erfc(delta_ta / (math.sqrt(2) * sig))
    r_ta = a * (R / h) * Mpc / one_ta ** (1 / 3)
    Hs = H0 * math.sqrt(E2_late(a))
    y = math.pi * Om_a(a) * fb * one_ta * r_ta * Hs / c
    return F * (fb * nu_s(y) - 1), F, y, nu_s(y)


P()
P("d) Umschaltmodell (Arbeitshypothese): eps(k,z) = A * F_dec(k,z) * [f_b nu(y_ta) - 1]")
P("   F_dec = erfc(delta_ta / (sqrt2 sigma(R=1/k, z))), delta_ta = 1.062; y_ta = pi Omega_m(z) f_b (1+Delta_ta) r_ta H/c, 1+Delta_ta = 5.55")
P(f"{'k[h/Mpc]':>8s} {'z':>4s} {'F_dec':>7s} {'y_ta':>9s} {'nu':>6s} {'f_b nu':>7s} {'eps(A=1)':>9s}")
for k in [0.02, 0.05, 0.1, 0.2]:
    for zz in [0.0, 0.5, 1.0, 2.0]:
        e, F, y, nv = eps_switch(k, 1 / (1 + zz))
        P(f"{k:8.2f} {zz:4.1f} {F:7.4f} {y:9.3e} {nv:6.2f} {fb*nv:7.3f} {e:9.4f}")

P()
P("   Wachstumsrate f(k,z) mit mu = 1 + A eps(k,a): relative Abweichung von LCDM")
zlist = [0.0, 0.5, 1.0, 1.5]
for A in [1.0, 0.1, 0.02]:
    P(f"   A = {A}")
    P(f"   {'k[h/Mpc]':>8s} " + " ".join([f"{'df/f z='+str(zz):>12s}" for zz in zlist]) + f" {'dD/D (z=0)':>15s}")
    for k in [0.02, 0.05, 0.1, 0.2]:
        aa, Dk, fk = growth(lambda a, k=k, A=A: 1 + A * eps_switch(k, a)[0])
        row = f"   {k:8.2f} "
        for zz in zlist:
            row += f" {100*(f_of_z(zz, fk)/f_of_z(zz)-1):11.3f}%"
        row += f" {100*(Dk[-1]/D0-1):14.3f}%"
        P(row)
P("   Befund: fuer kleines A ist die Abweichung proportional zu A; die Form in (k,z) ist vorgegeben.")

P()
P("e) Epochenabhaengige Schwelle direkt: RAR und Tully-Fisher bei z > 0 (a0 ~ H)")
P(f"{'z':>4s} {'a0(z)/a0(0)':>11s} {'tief: g_obs/g_bar rel.':>22s} {'dlog v (fest M) [dex]':>22s} {'dlog M (fest v) [dex]':>22s} {'nu(g=1e-11) z / z=0':>20s}")
for zz in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
    r = math.sqrt(E2_late(1 / (1 + zz)))
    a0z_ = a0_of_H(H0) * r
    g = 1e-11
    P(f"{zz:4.2f} {r:11.4f} {math.sqrt(r):22.4f} {0.25*math.log10(r):22.4f} {-math.log10(r):22.4f} {nu_s(g/a0z_)/nu_s(g/a0_of_H(H0)):20.4f}")

# ================================================================== 5 Topologie
P()
P("=" * 100)
P("5) Torus-Topologie und Kaleidoskop")
P("=" * 100)
chi = DM(LC, zstar)          # Mpc
d = 2 * chi
P(f"Abstand zur letzten Streuflaeche chi* = {chi/1e3:.3f} Gpc, Durchmesser d_LSS = {d/1e3:.3f} Gpc")
chi3 = DM(LC, 3.0)


def z_of_chi(ch):
    lo, hi = 0.0, 1100.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if DM(LC, mid, 4000) < ch:
            lo = mid
        else:
            hi = mid
    return mid


def count_clones(L, Rmax):
    n = int(Rmax / L) + 1
    cnt = 0
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                if (i, j, k) != (0, 0, 0) and L * math.sqrt(i * i + j * j + k * k) < Rmax:
                    cnt += 1
    return cnt


def jl(l, x):
    x = np.asarray(x, float)
    s, co = np.sin(x), np.cos(x)
    if l == 2:
        v = (3 / x ** 3 - 1 / x) * s - 3 * co / x ** 2
        ser = x ** 2 / 15 * (1 - x ** 2 / 14)
    elif l == 3:
        v = (15 / x ** 4 - 6 / x ** 2) * s - (15 / x ** 3 - 1 / x) * co
        ser = x ** 3 / 105 * (1 - x ** 2 / 18)
    else:
        v = (105 / x ** 5 - 45 / x ** 3 + 1 / x) * s - (105 / x ** 4 - 10 / x ** 2) * co
        ser = x ** 4 / 945 * (1 - x ** 2 / 22)
    return np.where(x < 0.3, ser, v)


def cut_fraction(l, xmin):
    x = np.linspace(1e-4, xmin, 20000)
    lost = np.trapezoid(jl(l, x) ** 2 / x, x)
    return 1 - lost * 2 * l * (l + 1)


P(f"Abstand bis z = 3 (Galaxien-Durchmusterungen): {chi3/1e3:.3f} Gpc")
P("Kubischer 3-Torus mit Kantenlaenge L:  l_min = pi d_LSS / L = 2 pi chi*/L;  Kreispaare wenn L < d_LSS, Winkelradius cos(alpha) = L/d_LSS")
P("C_l-Restanteil: kontinuierlicher Abschnitt bei k_min = 2 pi/L im Sachs-Wolfe-Integral (grobe Naeherung, ohne ISW und ohne Diskretheit)")
P(f"{'L/d':>6s} {'L[Gpc]':>7s} {'l_min':>6s} {'alpha[deg]':>10s} {'Klone<chi*':>10s} {'Klone<z=3':>9s} {'C2 Rest':>8s} {'C3 Rest':>8s} {'C4 Rest':>8s}")
U_author = 10e9 * LY / Mpc     # 10 Mrd Lj in Mpc
for Ld in [U_author / d, 0.25, 0.5, 0.75, 0.97, 0.985, 1.0, 1.2, 1.5, 2.0, 3.0]:
    L = Ld * d
    lmin = math.pi * d / L
    alpha = math.degrees(math.acos(L / d)) if L < d else float("nan")
    nc = count_clones(L, chi)
    n3 = count_clones(L, chi3)
    xm = 2 * math.pi * chi / L
    P(f"{Ld:6.3f} {L/1e3:7.2f} {lmin:6.2f} {alpha:10.1f} {nc:10d} {n3:9d} {cut_fraction(2, xm):8.3f} {cut_fraction(3, xm):8.3f} {cut_fraction(4, xm):8.3f}")
P(f"Kosmische Varianz des Quadrupols: sigma(C2)/C2 = sqrt(2/5) = {math.sqrt(0.4):.3f}")
P(f"Autor-Beispiel U = 10 Mrd Lj = {U_author/1e3:.3f} Gpc: erstes Eigenbild bei z = {z_of_chi(U_author):.3f}")
P(f"Grenze Planck 2015 (kubisch, R_i > 0.97 chi*): L > {0.97*2*chi/1e3:.2f} Gpc;  COMPACT 2023 (kuerzeste Schleife durch uns > 0.985 d): L > {0.985*d/1e3:.2f} Gpc")
P(f"COMPACT: anderswo kann die kuerzeste Schleife um Faktor 2 bis >=6 kuerzer sein: {0.985*d/2/1e3:.2f} Gpc bis {0.985*d/6/1e3:.2f} Gpc")

# ================================================================== 6 Atmen und w(z)
P()
P("=" * 100)
P("6) Atmen und Dunkle Energie: w(z)")
P("=" * 100)
P("DESI DR2 (arXiv:2503.14738), CPL w(a) = w0 + wa (1-a):")
desi = [("DESI+CMB+Pantheon+", -0.838, -0.62), ("DESI+CMB+Union3", -0.667, -1.09), ("DESI+CMB+DESY5", -0.752, -0.86)]
for n_, w0, wa in desi:
    ac = 1 + (1 + w0) / wa
    P(f"   {n_:20s} w0 = {w0:6.3f}, wa = {wa:6.2f}:  w = -1 bei z_c = {1/ac-1:.3f};  wa/(1+w0) = {wa/(1+w0):.2f};  w(z=1) = {w0+wa*0.5:.3f}, w(z=2) = {w0+wa*2/3:.3f}")
P("   w = -1 bei z_c bedeutet: die Dichte der Rueckstellung hat dort ein Maximum.")


def cpl_fit(zs, ws):
    a = 1 / (1 + np.asarray(zs))
    Amat = np.vstack([np.ones_like(a), 1 - a]).T
    sol, *_ = np.linalg.lstsq(Amat, np.asarray(ws), rcond=None)
    return sol


def w_from_rho(rho_func, zs, e=1e-4):
    out = []
    for zz in zs:
        a = 1 / (1 + zz)
        r1, r2 = rho_func(a * (1 - e)), rho_func(a * (1 + e))
        out.append(-1 - (1 / 3) * (math.log(r2) - math.log(r1)) / (math.log(1 + e) - math.log(1 - e)))
    return np.array(out)


zfit = np.linspace(0.0, 1.5, 61)

P()
P("a) Hypothese H-nu voll als scheinbare Dunkle Energie (Beobachter setzt Staub mit nu(z=10), Materie-Zeitalter, an):")
nuhi = float(NV.nu(10.0)[0])


def rho_DE_NV(a):
    zz = 1 / a - 1
    return NV.OL + (float(NV.nu(zz)[0]) - nuhi) * Ob * a ** -3


wv = w_from_rho(rho_DE_NV, zfit)
w0v, wav = cpl_fit(zfit, wv)
P(f"   nu(z=10) = {nuhi:.4f}, rho_DE,eff(z=0) = {rho_DE_NV(1.0):.4f}, w(z=0) = {wv[0]:.3f}, w(z=1) = {wv[40]:.3f}, CPL-Fit: w0 = {w0v:.3f}, wa = {wav:.3f}")

P()
P("b) Rueckstell-Ansatz: rho_DE(a) = rho_DE0 [ (1-beta) + beta S(a)/S(1) ],  S(a) = F_coll(>M,a) * sqrt(H(a)/H0) * a^(-m)")
P("   F_coll = erfc(1.686/(sqrt2 sigma(M) D(a)));  sqrt(H): tiefes Regime, Zusatzfaltung Delta s_eff ~ sqrt(G M a0) r ~ sqrt(a0)")
P(f"{'M[Msun/h]':>10s} {'m':>2s} {'sigma(M)':>8s} {'z_max(S)':>8s} {'beta':>5s} {'w0':>7s} {'wa':>7s} {'wa/(1+w0)':>9s}")
for M in [1e10, 1e11, 1e12, 3e12, 1e13, 1e14]:
    sM = sigmaR(R_of_M(M))
    for m_ in [0, 1]:
        def S(a, sM=sM, m_=m_):
            Da = float(np.interp(math.log(a), np.log(aG), DG)) / D0
            return math.erfc(1.686 / (math.sqrt(2) * sM * Da)) * math.sqrt(math.sqrt(E2_late(a))) * a ** (-m_)
        zg = np.linspace(0, 3, 301)
        Sv = [S(1 / (1 + zz)) for zz in zg]
        zmax = zg[int(np.argmax(Sv))]
        for beta in [0.1, 0.3]:
            def rho(a, beta=beta, S=S):
                return (1 - beta) + beta * S(a) / S(1.0)
            wz = w_from_rho(rho, zfit)
            w0_, wa_ = cpl_fit(zfit, wz)
            P(f"{M:10.0e} {m_:2d} {sM:8.3f} {zmax:8.2f} {beta:5.2f} {w0_:7.3f} {wa_:7.3f} {wa_/(1+w0_):9.2f}")

P()
P("c) beta so gewaehlt, dass w0 = -0.752 (DESY5-Kombination); welches wa folgt?")
for M in [1e11, 1e12, 3e12, 1e13, 1e14]:
    sM = sigmaR(R_of_M(M))
    for m_ in [0, 1]:
        def S(a, sM=sM, m_=m_):
            Da = float(np.interp(math.log(a), np.log(aG), DG)) / D0
            return math.erfc(1.686 / (math.sqrt(2) * sM * Da)) * math.sqrt(math.sqrt(E2_late(a))) * a ** (-m_)

        def fit(beta):
            wz = w_from_rho(lambda a: (1 - beta) + beta * S(a) / S(1.0), zfit)
            return cpl_fit(zfit, wz)
        lo, hi = 0.0, 1.0
        ok = fit(hi)[0] > -0.752
        if not ok:
            P(f"   M = {M:.0e}, m = {m_}: w0 = -0.752 nicht erreichbar (w0(beta=1) = {fit(hi)[0]:.3f})")
            continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if fit(mid)[0] > -0.752:
                hi = mid
            else:
                lo = mid
        w0_, wa_ = fit(hi)
        P(f"   M = {M:.0e}, m = {m_}: beta = {hi:.3f} -> w0 = {w0_:.3f}, wa = {wa_:.3f}   (DESY5: wa = -0.86 +0.23/-0.20)")

P()
P("Ende.")

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "kosmologie.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
print("\n".join(OUT))
