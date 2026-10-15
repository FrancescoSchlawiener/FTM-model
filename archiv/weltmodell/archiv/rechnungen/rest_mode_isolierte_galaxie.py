"""Rest der Mode ξ um eine isolierte Galaxie wie die Milchstraße, nach der Schließungsregel des Modells:
w = θ₊²/(θ₊² + θ_c²), Γ = (1 − w)·ω_orb/N*, ω_orb² = GM(<r)/r³. Zwei Wege:
 1) geschlossen: Umkehr im EdS-Kollaps t_ta = π√(r_ta³/8GM), danach w = 0, ω zwischen Umkehr- und Virialwert;
 2) numerisch: Schale im Hintergrund mit Dunkler Energie, w(θ) aus dem Kern, nach Erreichen von r_ta/2 fest."""
import sys, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, r"C:\DEV\FTM\weltmodell\rechnungen")
import ftm_kern as K

G, c, Ms, kpc, Gyr = K.G, K.c, K.M_sun, K.kpc, 3.15576e16
H0, a0 = K.H0_si(), K.A0
Om, OL = K.OMEGA_M, 1 - K.OMEGA_M
Q = K.Q_KOSM
Mb = 6e10 * Ms
nu = lambda y: 1 / (-math.expm1(-math.sqrt(y)))

# Hintergrund
Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_of_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a)[0]
t0 = t_of_a(1.0)

# Masse innerhalb der Schale vor der Übergabe: staubartig, NFW 10¹² M☉, c = 10 (Baryonen + Mode), oder Lesung der flachen Kurve
rho_c = 3 * H0**2 / (8 * math.pi * G)
Mvir, cc = 1e12 * Ms, 10
rvir = (3 * Mvir / (4 * math.pi * 200 * rho_c)) ** (1 / 3)
mu = lambda x: math.log(1 + x) - x / (1 + x)
M_nfw = lambda r: Mvir * mu(r * cc / rvir) / mu(cc)
vf = K.v_flach(Mb)
M_flach = lambda r: vf**2 * r / G

def weg1(r, M, Ns):
    rta = 2 * r
    tta = math.pi * math.sqrt(rta**3 / (8 * G * M))
    w_ta, w_vir = math.sqrt(G * M / rta**3), math.sqrt(G * M / r**3)
    dt = t0 - tta
    return tta, math.exp(-dt * w_ta / Ns), math.exp(-dt * w_vir / Ns)

def weg2(r, M, Ns):
    rta = 2 * r
    ai = 1e-3; ti = t_of_a(ai)
    rho_i = Om * rho_c / ai**3
    def lauf(di, voll=False):
        ri = (3 * M / (4 * math.pi * rho_i * (1 + di))) ** (1 / 3)
        vi = Hb(ai) * ri * (1 - di / 3)
        def f(t, y):
            R, V, L = y
            if R <= rta / 2 and V < 0:
                return [0, 0, -math.sqrt(G * M / (rta / 2)**3) / Ns]
            th = 3 * V / R
            w = K.gewicht_w(th)
            om = math.sqrt(G * M / R**3)
            return [V, -G * M / R**2 + OL * H0**2 * R, -(1 - w) * om / Ns]
        ev = lambda t, y: y[1]; ev.terminal = not voll; ev.direction = -1
        sol = solve_ivp(f, (ti, t0), [ri, vi, 0.0], events=ev, rtol=1e-8, atol=1e-3, max_step=0.02 * Gyr)
        return sol
    def rmax(di):
        s = lauf(di)
        return (s.y_events[0][0][0] if len(s.t_events[0]) else s.y[0].max()) - rta
    di = brentq(rmax, 1e-4, 0.5, xtol=1e-8)
    s = lauf(di, voll=True)
    # nach r_ta/2: ω fest (Virialwert), bis t0
    L = s.y[2][-1]
    return s.t_events[0][0], math.exp(L)

print(f"t0 = {t0/Gyr:.2f} Mrd. J., r_vir(10¹²) = {rvir/kpc:.0f} kpc, v_flach = {vf/1e3:.1f} km/s")
print("r [kpc] | Masse | t_Umkehr (1) | ξ_R Weg 1 (ω Umkehr…Virial) N* = 0,3 / 2,4 | t_Umkehr (2) | ξ_R Weg 2 N* = 0,3 / 2,4 | δg/g bei N* = 2,4 (Weg 2)")
for r_k in (100, 200, 300):
    r = r_k * kpc
    for lab, Mf in (("NFW", M_nfw), ("flach", M_flach)):
        M = Mf(r)
        t1, a1, b1 = weg1(r, M, 0.3); _, a2, b2 = weg1(r, M, 2.4)
        t2, x03 = weg2(r, M, 0.3); _, x24 = weg2(r, M, 2.4)
        y = G * Mb / r**2 / a0; n = nu(y)
        dg = -x24 * ((n - 1) - Q * min(r / rvir, 1)) / n
        print(f"{r_k} | {lab} {M/Ms:.1e} | {t1/Gyr:.1f} | {b1:.1e}…{a1:.1e} / {b2:.2f}…{a2:.2f} | {t2/Gyr:.1f} | {x03:.1e} / {x24:.2f} | {dg:+.2f} (ν = {n:.1f})")
