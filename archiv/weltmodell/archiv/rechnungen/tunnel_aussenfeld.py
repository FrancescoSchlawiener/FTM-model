# -*- coding: utf-8 -*-
"""
Tunnelbild der Umgebungsregel. Barriere V = ħ|g_N|/c, s² = V/(M_min c²) = |g_N|/a₀.
Zusammensetzung zweier Felder (inneres g_i, äußeres g_e):
  Vektor (Regel des Modells): V = ħ|g_i + g_e|/c  → senkrecht s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e²
  Matthiessen (unabhängige Verluste, Raten addiert): s² = s_i² + s_e² in jeder Richtung
  Reihe (Barrieren hintereinander, Amplituden multipliziert): q = q_i·q_e, s = s_i + s_e
Signatur der Vektorform: liegt g_e in der Scheibenebene, hängt ν vom Azimut ab (einseitige Rotation, m = 1 und m = 2).
Im frei fallenden System der Galaxie: a_int(x) = ν(|g_i + g_e|)(g_i + g_e) − ν(g_e)g_e (Zentrum als Bezug).
Punktmasse 6·10¹⁰ M☉, g_e = 0,03·a₀ (Obergrenze des mittleren äußeren Felds im Modell).
"""
import math
import numpy as np
from ftm_kern import G, M_sun, pc, A0, nu

kpc = 1e3 * pc
M = 6e10 * M_sun


def nu_s(s):
    return 1 / (-np.expm1(-s))


if __name__ == "__main__":
    ge = 0.03 * A0
    print("1) ν bei 50 kpc (Milchstraße, g_e = 0,03·a₀) je Zusammensetzung:")
    gi = G * M / (50 * kpc)**2; si, se = math.sqrt(gi / A0), math.sqrt(ge / A0)
    print("   frei %.3f | Vektor senkrecht %.3f | Vektor gleichgerichtet %.3f | Matthiessen %.3f | Reihe %.3f"
          % (nu_s(si), nu_s((si**4 + se**4)**0.25), nu_s(math.sqrt(si**2 + se**2)), nu_s(math.sqrt(si**2 + se**2)), nu_s(si + se)))
    print("   (die Reihe senkt die Verstärkung am stärksten; nur die Vektorform hängt an der Richtung)")
    print("\n2) Signatur der Vektorform: g_e in der Scheibenebene; v_c²(φ) = −R·a_int·e_R; Fourier-Anteile von v_c(φ)")
    print("   R [kpc] | g_i/a₀ | v_c Mittel [km/s] | frei [km/s] | m = 1 Amplitude [%] | m = 2 Amplitude [%] | senkrechtes g_e: m = 1, m = 2")
    phi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    for Rk in (5, 10, 20, 30, 50):
        R = Rk * kpc; gi = G * M / R**2
        for lage in ("Ebene", "senkrecht"):
            eR = np.stack([np.cos(phi), np.sin(phi), 0 * phi])
            gie = -gi * eR
            gev = np.array([ge, 0, 0]) if lage == "Ebene" else np.array([0, 0, ge])
            gt = gie + gev[:, None]
            a = nu(np.linalg.norm(gt, axis=0)) * gt - (nu(ge) * gev)[:, None]
            vc2 = -R * np.sum(a * eR, axis=0)
            vc = np.sqrt(np.clip(vc2, 0, None)) / 1e3
            F = np.fft.rfft(vc) / len(vc)
            if lage == "Ebene":
                m1, m2, v0 = 2 * abs(F[1]) / F[0].real * 100, 2 * abs(F[2]) / F[0].real * 100, F[0].real
            else:
                p1, p2 = 2 * abs(F[1]) / F[0].real * 100, 2 * abs(F[2]) / F[0].real * 100
        vf = math.sqrt(nu(gi) * gi * R) / 1e3
        print("   %6.0f  | %6.3f | %17.1f | %11.1f | %19.2f | %19.2f | %.1e, %.1e" % (Rk, gi / A0, v0, vf, m1, m2, p1, p2))

    print("\n3) Bezug der Galaxie als Ganzes: B = ν(g_e)·g_e (Feldform, Impuls erhalten) gegen A = g_e (Masse im Newtonschen Zentrum, ν = 1)")
    print("   R [kpc] | m = 1 [%] Bezug B | m = 1 [%] Bezug A | Relativbeschleunigung außen/innen (A) in g_i·ν")
    for Rk in (5, 10, 20, 30, 50):
        R = Rk * kpc; gi = G * M / R**2
        eR = np.stack([np.cos(phi), np.sin(phi), 0 * phi]); gev = np.array([ge, 0, 0]); gt = -gi * eR + gev[:, None]
        res = []
        for ref in (nu(ge) * gev, gev):
            a = nu(np.linalg.norm(gt, axis=0)) * gt - ref[:, None]
            vc = np.sqrt(np.clip(-R * np.sum(a * eR, axis=0), 0, None))
            F = np.fft.rfft(vc) / len(vc); res.append(2 * abs(F[1]) / F[0].real * 100)
        rel = (nu(gi) - 1) * ge / (nu(gi) * gi)
        print("   %6.0f  | %17.1f | %17.1f | %.3f" % (Rk, res[0], res[1], rel))
