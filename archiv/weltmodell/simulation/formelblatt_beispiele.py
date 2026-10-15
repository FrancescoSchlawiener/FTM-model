# -*- coding: utf-8 -*-
"""Beispielzahlen des Formelblatts, live aus dem Kern. Aufruf: python formelblatt_beispiele.py"""
import math
import kern as K

MW = 6e10 * K.M_sun
LAM = K.schliesslaenge()
H0 = K.H0_si()
hbar = K.hbar


def zeile(text, wert):
    print(f"{text:58s} {wert}")


print("Rahmen (parametrisch)")
zeile("Δs Erde [mm]", f"{K.ds(K.M_earth) * 1e3:.4g}")
zeile("Δs Sonne [m]", f"{K.ds(K.M_sun):.4g}")
zeile("Δs Milchstraße [AE]", f"{K.ds(MW) / K.AU:.4g}")
zeile("Λ = c²/2a₀ [Gpc]", f"{LAM / K.Gpc:.3g}")
zeile("τ_D = 2Λ/c [Mrd. J.]", f"{K.tau_D() / K.Gyr:.3g}")

print("Rückkehr (transversal)")
for r_kpc in (2, 8.7, 30):
    aN = K.G * MW / (r_kpc * K.kpc)**2
    zeile(f"Milchstraße r = {r_kpc} kpc: s, q, ν", f"{float(K.s_(aN)):.3f}, {float(K.q_(aN)):.3f}, {float(K.nu(aN)):.3f}")
zeile("r* Sonne [AE]", f"{K.r_stern(K.M_sun) / K.AU:.4g}")
zeile("r* Milchstraße [kpc]", f"{K.r_stern(MW) / K.kpc:.3g}")
zeile("R_Q Milchstraße [kpc]", f"{K.r_Q(MW) / K.kpc:.3g}")
zeile("v_flach Milchstraße [km/s]", f"{K.v_flach(MW) / 1e3:.4g}")
zeile("RC-Kette ν (n_R = 20, RG = 0,01)", f"{K.rc_kette_nu(20, 0.01):.5g}")

print("Leiter")
zeile("r*² / (3Δs·ℓ₀), Sonne", f"{K.r_stern(K.M_sun)**2 / (3 * K.ds(K.M_sun) * K.ell0()):.6f}")
lp2 = hbar * K.G / K.c**3
lam_p = hbar / (K.m_p * K.c)
zeile("ℓ_P² / (3Δs·λ̄), Proton", f"{lp2 / (3 * K.ds(K.m_p) * lam_p):.6f}")

print("Ganzes")
zeile("Q = Ω_c/Ω_b, q_kosm", f"{K.Q_KOSM:.3f}, {K.q_KOSM:.4f}")
zeile("Zählung je Radiant 2π − 1", f"{2 * math.pi - 1:.3f}")
zeile("Θ = 2πa₀/(cH₀), Ω_s = (π²/45)Θ²", f"{K.theta_(K.A0, H0):.3f}, {K.omega_vakuum():.3f}")
