# -*- coding: utf-8 -*-
"""
Präzession geneigter Ringe unter der örtlichen Regel a = ν(|g_N|)·g_N (07, offen: Polarringe).
Scheibe: Miyamoto-Nagai, M = 6·10¹⁰ M☉, a = 3 kpc, b = 0,3 kpc (Annahme). Ring = Kreisbahn mit Neigung i.
Weg A: 3D-Bahnintegration (Newton und örtliche Regel), Knotenpräzession aus dem Drehimpulsvektor.
Weg B: kleine Neigung analytisch, ω_p = Ω − ν_z mit Frequenzen aus Differenzen des Felds (ν_z² = ν·ν_z,N², Ω² = ν·Ω_N²).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from ftm_kern import G, M_sun, pc, A0, nu, c

kpc = 1e3 * pc
GM = G * 6e10 * M_sun
aMN, bMN = 3 * kpc, 0.3 * kpc
KMSKPC = 1e3 / kpc


def gN(x):
    R2 = x[0]**2 + x[1]**2
    zb = math.sqrt(x[2]**2 + bMN**2)
    D = (R2 + (aMN + zb)**2) ** 1.5
    return np.array([-GM * x[0] / D, -GM * x[1] / D, -GM * x[2] * (aMN + zb) / (zb * D)])


def acc(x, regel):
    g = gN(x)
    return nu(np.linalg.norm(g)) * g if regel else g


def omega_kreis(R, regel):
    g = acc(np.array([R, 0, 0]), regel)
    return math.sqrt(-g[0] / R)


def nu_z(R, regel, h=1e-4 * kpc):
    gz = acc(np.array([R, 0, h]), regel)[2]
    return math.sqrt(-gz / h)


def praez(R, inc_deg, regel, n_orb=12):
    """Knotenpräzession [km/s/kpc] aus der Drehung des Drehimpulsvektors (gemittelt über n_orb Umläufe)."""
    inc = math.radians(inc_deg)
    Om = omega_kreis(R, regel)
    # Start: Ort auf der Knotenlinie, Geschwindigkeit geneigt; Kreisgeschwindigkeit aus dem Kraftbetrag an diesem Ort
    x0 = np.array([R, 0.0, 0.0])
    v = Om * R
    v0 = np.array([0.0, v * math.cos(inc), v * math.sin(inc)])
    T = 2 * math.pi / Om
    f = lambda t, y: np.concatenate([y[3:], acc(y[:3], regel)])
    sol = solve_ivp(f, (0, n_orb * T), np.concatenate([x0, v0]), method='DOP853', rtol=1e-10, atol=1e-6, dense_output=True)
    ts = np.linspace(0, n_orb * T, 40 * n_orb + 1)
    Y = sol.sol(ts)
    L = np.cross(Y[:3].T, Y[3:].T)
    phi = np.unwrap(np.arctan2(L[:, 1], L[:, 0]))       # Azimut der Knotenlinie + 90°
    k = np.polyfit(ts, phi, 1)[0]
    incl = np.degrees(np.arccos(L[:, 2] / np.linalg.norm(L, axis=1)))
    return k / KMSKPC, incl.mean(), Om / KMSKPC


if __name__ == "__main__":
    print("Miyamoto-Nagai M = 6·10¹⁰ M☉, a = 3 kpc, b = 0,3 kpc; a₀ = %.2e m/s²" % A0)
    print("\nWeg B, kleine Neigung: ω_p = Ω − ν_z [km/s/kpc]; Verhältnis örtlich/Newton gegen √ν")
    print("  R [kpc] |  ν(g_R) | ω_p Newton | ω_p örtlich | Verhältnis | √ν")
    for Rk in (8, 12, 20, 30, 40):
        R = Rk * kpc
        gR = -gN(np.array([R, 0, 0]))[0]
        n_ = nu(gR)
        pN = (omega_kreis(R, False) - nu_z(R, False)) / KMSKPC
        pL = (omega_kreis(R, True) - nu_z(R, True)) / KMSKPC
        print("  %6.0f  | %7.4f | %10.4f | %11.4f | %10.5f | %.5f" % (Rk, n_, pN, pL, pL / pN, math.sqrt(n_)))

    print("\nWeg A, Bahnintegration: Knotenpräzession ω_p [km/s/kpc] je Neigung (mittlere Neigung in Klammern)")
    print("  R [kpc] | i [°] | ω_p Newton          | ω_p örtlich         | örtlich/Newton | √ν(g_R) | Weg B (i → 0)")
    for Rk in (12, 20, 30):
        R = Rk * kpc
        gR = -gN(np.array([R, 0, 0]))[0]
        pB = (omega_kreis(R, True) - nu_z(R, True)) / KMSKPC
        for i in (2, 30, 60, 80, 89):
            pN, iN, _ = praez(R, i, False)
            pL, iL, _ = praez(R, i, True)
            print("  %6.0f  | %5.0f | %8.4f (%5.1f°)   | %8.4f (%5.1f°)   | %14.4f | %7.4f | %.4f" % (Rk, i, pN, iN, pL, iL, pL / pN, math.sqrt(nu(gR)), pB))
