# -*- coding: utf-8 -*-
"""
Prüfstand: rechnet die Zahlen aus modell/01–09 mit dem einen Kern (ftm_kern.py) nach.
Aufruf: python ftm_pruefstand.py   (nach jeder Änderung am Modell)
Jede Zeile: Baustein | Aussage | Wert im Modell | Kern | Abweichung | OK/ABWEICHUNG
"""
import io, math, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import *

ARC = 180 / math.pi * 3600
MW = 6e10 * M_sun
CL = 1.7e14 * M_sun
a0L = A0
a0M = a0_modell()
H0 = H0_si()

checks = []
def P(b, text, modell, kern, tol=0.01):
    checks.append((b, text, modell, kern, tol))

# ---- 01 Grundlagen
P("01", "Δs Erde [mm]", 1.478, ds(M_earth) * 1e3, 0.002)
P("01", "Δs Sonne [m]", 492.2, ds(M_sun), 0.001)
P("01", "g Sonne [m/s²]", 274, neigung(M_sun, R_sun), 0.005)
P("01", "r₀ Erde [mm]", 8.87, r0(M_earth) * 1e3, 0.003)
P("01", "r₀ Sonne [km]", 2.95, r0(M_sun) / 1e3, 0.005)
P("01", "δ Erdoberfläche", 6.96e-10, delta(M_earth, R_earth), 0.005)
P("01", "Erdoberfläche [s/Jahr]", 0.022, delta(M_earth, R_earth) * yr, 0.02)
P("01", "δ Sonnenoberfläche", 2.12e-6, delta(M_sun, R_sun), 0.005)
P("01", "Sonnenoberfläche [s/Jahr]", 67, delta(M_sun, R_sun) * yr, 0.01)
# ---- 02 Teilchen bis Stern
P("02", "Δs Proton [m]", 4.14e-55, ds(m_p), 0.003)
P("02", "Δs/R Proton", 4.9e-40, ds(m_p) / R_p, 0.02)
P("02", "g Proton [m/s²]", 1.6e-7, neigung(m_p, R_p), 0.03)
P("02", "Δs/R Wasserstoff-Atom", 7.8e-45, ds(m_p) / a_B, 0.01)
P("02", "Δs Neutronenstern 1,4 M☉ [m]", 689, ds(1.4 * M_sun), 0.002)
P("02", "Δs/R Neutronenstern (12 km)", 0.057, ds(1.4 * M_sun) / 12e3, 0.01)
P("02", "g Neutronenstern [m/s²]", 1.3e12, neigung(1.4 * M_sun, 12e3), 0.02)
P("02", "Δs Roter Riese 18 M☉ [km]", 8.86, ds(18 * M_sun) / 1e3, 0.002)
P("02", "Δs/R Roter Riese (764 R☉)", 1.7e-8, ds(18 * M_sun) / (764 * R_sun), 0.03)
P("02", "g Roter Riese [m/s²]", 0.0084, neigung(18 * M_sun, 764 * R_sun), 0.01)
P("02", "α_G", 5.906e-39, alpha_G(), 0.001)
P("02", "Sternmasse M_* [M☉]", 1.853, sternmasse() / M_sun, 0.002)
P("02", "Δs_* [m]", 912, ds(sternmasse()), 0.002)
P("02", "Bahngeschwindigkeit Erde [km/s]", 29.79, math.sqrt(3 * c**2 * ds(M_sun) / AU) / 1e3, 0.001)
P("02", "r* Sonne [AE]", 7340, r_stern(M_sun) / AU, 0.001)
P("02", "Neutronenstern: Hülle bei r/r₀", 2.9, 12e3 / r0(1.4 * M_sun), 0.01)
P("02", "Neutronenstern: Takt langsamer um", 0.19, delta(1.4 * M_sun, 12e3), 0.01)
# ---- 03 Galaxie bis Kosmos
P("03", "Δs Milchstraße [AE]", 197.5, ds(MW) / AU, 0.002)
P("03", "Δs/R Milchstraße (15,3 kpc)", 6.2e-8, ds(MW) / (15.3 * kpc), 0.02)
P("03", "Δs Galaxienhaufen [Lj]", 8.8, ds(CL) / ly, 0.01)
P("03", "Δs/R Galaxienhaufen (2 Mpc)", 1.4e-6, ds(CL) / (2 * Mpc), 0.04)
P("03", "Δs Universum [Mrd Lj]", 3.9, ds(7.5e22 * M_sun) / ly / 1e9, 0.02)
P("03", "Δs/R Universum (46,5 Mrd Lj)", 0.084, ds(7.5e22 * M_sun) / (46.5e9 * ly), 0.01)
P("03", "r* Milchstraße [kpc]", 8.72, r_stern(MW) / kpc, 0.003)
for r, v in ((5, 250), (r_stern(MW) / kpc, 216), (15, 198), (50, 180)):
    P("03", f"v Milchstraße bei {r:.3g} kpc [km/s]", v, v_kreis_punkt(MW, r * kpc) / 1e3, 0.005)
P("03", "v_flach Milchstraße [km/s]", 172.0, v_flach(MW) / 1e3, 0.001)
P("03", "Grenzflächendichte a₀/(πG) [M☉/pc²]", 251, sigma_kollektiv(), 0.005)
ai = G * MW / (50 * kpc)**2
P("03", "ν frei (Milchstraße 50 kpc)", 6.25, float(nu(ai)), 0.002)
P("03", "ν senkrechtes Außenfeld 0,03 a₀", 5.36, nu_aussenfeld(ai, 0.03 * a0L, "senkrecht"), 0.002)
P("03", "ν gleichgerichtetes Außenfeld", 4.59, nu_aussenfeld(ai, 0.03 * a0L, "gleich"), 0.002)
lq = -math.log(q_KOSM)
P("03", "R_Q: a_N/a₀ beim kosmischen q", 0.029, lq**2, 0.01)
P("03", "R_Q Milchstraße [kpc]", 51, r_stern(MW) / lq / kpc, 0.02)
P("03", "1 − w heute (Standardwert a₀)", 1.46e-2, 1 - gewicht_w(3 * H0, a0L), 0.01)
P("03", "1 − w heute (a₀ = cH₀/2π)", 1.30e-2, 1 - gewicht_w(3 * H0, a0M), 0.01)
P("03", "ν an der Hubble-Kugel (g = πa₀)", 1.2, float(nu(math.pi * a0M, a0M)), 0.01)
P("03", "cH₀/2π bei 67,4 [10⁻¹⁰]", 1.04, a0M * 1e10, 0.005)
P("03", "cH₀/2π bei 70 [10⁻¹⁰]", 1.08, a0_modell(70) * 1e10, 0.005)
P("03", "H₀ für a₀ = 1,2·10⁻¹⁰ (Literatur) [km/s/Mpc]", 77.6, 2 * math.pi * A0_LIT / c * Mpc / 1e3, 0.002)
P("03", "Vakuumteil Ω_s bei Θ = 1", 0.22, K_S, 0.01)
P("03", "Vakuumteil mit dem Standardwert a₀", 0.24, omega_vakuum(a0L), 0.02)
P("03", "ΔN_eff bei Ω_s = 0,22 zu jeder Zeit", 2.09, delta_N_eff(K_S), 0.005)
dL = d_LSS()
P("03", "Durchmesser letzte Streufläche [Gpc]", 27.8, dL / Gpc, 0.003)
P("03", "Zelle für a₀ = 1,04·10⁻¹⁰ [Gpc]", 28.1, zelle_aus_a0(1.04e-10) / Gpc, 0.003)
P("03", "Zelle für a₀ = 1,19·10⁻¹⁰ [Gpc]", 26.9, zelle_aus_a0(1.19e-10) / Gpc, 0.003)
P("03", "a₀ ohne Kreispaare (L ≥ 0,985 d_LSS) [10⁻¹⁰]", 1.13, a0_aus_zelle(KREISPAAR_GRENZE * dL) * 1e10, 0.005)
P("03", "a₀(z = 10)/a₀(0)", 0.997, a0_frueh_relativ(10), 0.001)
P("03", "a₀(z = 1090)/a₀(0)", 0.76, a0_frueh_relativ(1090), 0.005)
a0c = a0_aus_zelle(dL) * a0_frueh_relativ(1090)
H1090 = H0 * E_hintergrund(1090)
P("03", "Vakuumteil bei z = 1090", 2.5e-10, K_S * theta_(a0c, H1090)**2 / E_hintergrund(1090)**0 / 1, 0.05)
P("03", "Schließlänge Λ für den Standardwert [Gpc]", 13.2, schliesslaenge(a0L) / Gpc, 0.01)
P("03", "Schließlänge Λ für a₀ = 1,04 [Gpc]", 14.0, schliesslaenge(1.04e-10) / Gpc, 0.01)
w_, lo, hi = uebergabezeit(500)
P("03", "1/λ_m heute, 500 ρ_krit, N* = 0,3 [Mrd J]", 0.28, lo / Gyr, 0.02)
w_, lo, hi = uebergabezeit(200)
P("03", "1/λ_m heute, 200 ρ_krit, N* = 2,4 [Mrd J]", 3.5, hi / Gyr, 0.02)
w_, lo, hi = uebergabezeit(500, z=2)
P("03", "1/λ_m Kollaps z = 2, 500 ρ_krit, N* = 0,3 [Mrd J]", 0.1, lo / Gyr, 0.12)
w_, lo, hi = uebergabezeit(200, z=2)
P("03", "1/λ_m Kollaps z = 2, 200 ρ_krit, N* = 2,4 [Mrd J]", 1.2, hi / Gyr, 0.05)
rate, S, Hend = kreislauf_mindestrate(schliesslaenge(a0L))
P("03", "Kreislauf: Mindestrate Λ̇/Λ [s⁻¹]", 3.8e-37, rate, 0.02)
P("03", "Kreislauf: je Jahr", 1.2e-29, rate * yr, 0.05)
P("03", "Kreislauf: w_s + 1", 2e-19, 2 / 3 * rate / Hend, 0.2)
P("03", "Kreislauf: S_Λ [k]", 9.3e122, S, 0.01)
P("03", "Kreislauf: H im Endzustand / H₀", 0.49, Hend / H0, 0.02)
# ---- 04
P("04", "Neutronenstern: Tage je Jahr", 69.6, delta(1.4 * M_sun, 12e3) * 365.25, 0.005)
P("04", "Grenze Rückgabe T_H = T_CMB [kg]", 4.5e22, grenzmasse_hawking(), 0.005)
# ---- 07
P("07", "ν bei r*", 1.582, float(nu(a0L)), 0.001)
P("07", "Bose-Stufen bei r* (s = 1)", 0.582, float(n_B(1.0)), 0.001)
s2, aN2, v2 = bose_stufe_max(2)
P("07", "Maximum d = 2: a_N/a₀", 2.54, aN2, 0.002)
P("07", "Maximum d = 2: Zusatz/a₀", 0.648, v2, 0.002)
P("07", "Maximum d = 3: a_N/a₀", 7.96, bose_stufe_max(3)[1], 0.002)
P("07", "beide Umlaufsinne: coth(1)", 1.313, 1 / math.tanh(1), 0.001)
P("07", "τ_D = c/a₀ [Mrd J]", 86, tau_D() / Gyr, 0.005)
P("07", "Ebenen bis r* Milchstraße", 1.5e6, r_stern(MW) / r0(MW), 0.05)
P("07", "Δs/r an der Erdbahn", 3.3e-9, ds(M_sun) / AU, 0.01)
P("07", "Δs/r Milchstraße bei 8 kpc", 1.2e-7, ds(MW) / (8 * kpc), 0.01)
P("07", "Im τ_M an der Erdbahn", 3.1, math.log(AU / ds(M_sun)) / (2 * math.pi), 0.01)
P("07", "Im τ_M Milchstraße 8 kpc", 2.5, math.log(8 * kpc / ds(MW)) / (2 * math.pi), 0.02)
P("07", "Rosette χ = 1: Winkel zwischen Perizentren [°]", 254.6, 360 / math.sqrt(2), 0.001)
P("07", "Erdbahn: s = r*/r", 7342, r_stern(M_sun) / AU, 0.001)
# ---- 08
P("08", "Ablenkung Sonnenrand [″]", 1.7516, ablenkung(M_sun, R_sun) * ARC, 0.001)
P("08", "r* Erde [AE]", 12.7, r_stern(M_earth) / AU, 0.005)
P("08", "Masse mit r* = λ̄ [kg]", 5.9e-29, masse_r_stern_gleich_compton(), 0.02)
P("08", "α_∞ Milchstraße [″]", 0.43, ablenkung_kollektiv(MW) * ARC, 0.02)
P("08", "α_∞ Galaxienhaufen [″]", 23, ablenkung_kollektiv(CL) * ARC, 0.02)
P("08", "Merkur-Perihel [″/Jh]", 43, perihel_pro_umlauf(M_sun, 5.7909e10, 0.2056) * ARC * 36525 / 87.969, 0.005)
aJ = 8.034e37 * 7.2921e-5 / (M_earth * c)
P("08", "Verwindung Erde a_J [m]", 3.27, aJ, 0.005)
P("08", "T_prograd − T_retrograd [ns]", 137, 4 * math.pi * aJ / c * 1e9, 0.005)
P("08", "Hulse-Taylor dP/dt", -2.40e-12, hulse_taylor(1.4398, 1.3886, 27906.98, 0.6171338), 0.005)
P("08", "Schatten M87* [μas]", 40, schatten_winkel(6.5e9 * M_sun, 16.8 * Mpc) * ARC * 1e6, 0.02)
P("08", "Erdbahn im Sonnentrichter [s/Jahr]", 0.31, 3 * ds(M_sun) / AU * yr, 0.01)
# ---- 09
P("09", "a_J Elektron [fm]", 193.08, compton_red(m_e) / 2 * 1e15, 0.0005)
P("09", "Bohr-Magneton [J/T]", 9.274e-24, e_el * hbar / (2 * m_e), 0.0005)
P("09", "Spinzahl Elektron", 2.9e44, 0.5 / alpha_G(m_e), 0.02)
P("09", "Spinzahl Proton", 8.5e37, 0.5 / alpha_G(m_p), 0.01)
P("09", "α/α_G", 1.2355e36, alpha / alpha_G(), 0.001)
P("09", "ω₁ = c/λ̄_p [rad/s]", 1.425e24, c / compton_red(m_p), 0.001)
P("09", "ω₂ = αc/a_B [rad/s]", 4.134e16, alpha * c / a_B, 0.001)
P("09", "Γ_G = α_G·ω₁ [1/s]", 8.42e-15, alpha_G() * c / compton_red(m_p), 0.002)


# ================================================================== Erweiterung: ganzes Modell
from scipy.integrate import quad
from scipy.optimize import brentq
# ---- 02
P("02", "δ bei 1 fm neben einem Proton", 1.24e-39, 3 * ds(m_p) / 1e-15, 0.01)
P("02", "Umlauffrequenz Wasserstoff-Hülle [Hz]", 6.58e15, alpha * c / (2 * math.pi * a_B), 0.002)
P("02", "Bindungssteilheit Wasserstoff [m/s²]", 9e22, (alpha * c)**2 / a_B, 0.02)
P("02", "Bindungssteilheit Wasserstoff [g⊕]", 9e21, (alpha * c)**2 / a_B / g_earth, 0.03)
P("02", "Faltung für diese Steilheit [kg]", 3.8e12, (alpha * c)**2 / a_B * a_B**2 / G, 0.01)
T_Fe, m_Fe = 1811.0, 55.845 * 1.66053906660e-27      # Schmelztemperatur Eisen [G]
P("02", "Erde: δ gegen 2(v_th/c)² (Eisen, Schmelzpunkt)", 116, delta(M_earth, R_earth) / (2 * kB * T_Fe / (m_Fe * c**2)), 0.01)
P("02", "Saturn: Faltung / Erde", 95, 3.7931187e16 / GM_earth, 0.005)
P("02", "Saturn: Kante / Erde (mittlerer Radius 58 232 km)", 1.14, (3.7931187e16 / 5.8232e7**2) / g_earth, 0.01)
M_pl = (1.89813e27 + 5.6832e26 + 8.6811e25 + 1.02409e26 + M_earth + 4.8675e24 + 6.4171e23 + 3.3011e23)
P("02", "Faltung des Sonnensystems [m]", 492.9, ds(M_sun + M_pl), 0.0005)
P("02", "Anteil der Sonne [%]", 99.87, 100 * M_sun / (M_sun + M_pl), 0.0002)
P("02", "Glättungslänge Sonne / r*", 26, 0.93 * pc / r_stern(M_sun), 0.02)
_wC = m_p * c**2 / hbar
_fs = kB * TEFF_sun / (alpha**2 * m_e * c**2); _LE = 4 * math.pi * c * G * M_sun * m_p / (8 * math.pi / 3 * r_e**2)
_sqC = (2 * math.pi**3 / 45) * _fs**4 * (M_sun / sternmasse()) * (alpha**2 / (3 * ds(M_sun) / R_sun))**2 * _LE / L_sun
P("02", "Carter-Identität √𝒞 (Sonne) gegen 𝒞 direkt", 2.944, _sqC**2, 1e-4)
P("02", "f_s = kT_eff/(α²m_ec²) Sonne", 0.0183, _fs, 0.002)
P("02", "L_Edd/L Sonne", 3.28e4, _LE / L_sun, 0.002)
P("09", "Faltungstakt Proton c/(3Δs_p) [s⁻¹]", 2.41e62, c / (3 * ds(m_p)), 0.005)
P("09", "ω_C² = Γ_G·c/(3Δs) (Proton)", 1.0, _wC**2 / (alpha_G() * _wC * c / (3 * ds(m_p))), 1e-9)
P("09", "Neutronenstern Umlauf an der Oberfläche [rad/s]", 1.04e4, math.sqrt(G * 1.4 * M_sun / 12e3**3), 0.005)
P("09", "Erdbahn [rad/s]", 1.99e-7, 2 * math.pi / yr, 0.005)
P("09", "Milchstraße Ω bei r* [rad/s]", 8.0e-16, v_kreis_punkt(6e10 * M_sun, r_stern(6e10 * M_sun, a0L), a0L) / r_stern(6e10 * M_sun, a0L), 0.01)
P("09", "Milchstraße Umlauf bei r* [Mio. Jahre]", 248, 2 * math.pi * r_stern(6e10 * M_sun, a0L) / v_kreis_punkt(6e10 * M_sun, r_stern(6e10 * M_sun, a0L), a0L) / yr / 1e6, 0.005)
P("03", "Θ heute (a₀ = 1,1)", 1.06, 2 * math.pi * a0L / (c * H0), 0.005)
import sparc_flaechendichte as _sfd
_t0 = _sfd.trend(_sfd.main()); _t1 = _sfd.trend(_sfd.main(ge=0.03))
_r4 = _sfd.trend(_sfd.main(0.4)); _r6 = _sfd.trend(_sfd.main(0.6))
P("03", "Σ_dyn(0) gemessen/Regel, Mittel [dex]", -0.05, _t0[0], 0.1)
P("03", "Σ_dyn(0) Rest bei Υ_d = 0,4 [dex]", -0.01, _r4[0], 0.5)
P("03", "Σ_dyn(0) Rest bei Υ_d = 0,6 [dex]", -0.09, _r6[0], 0.02)
P("03", "Σ_dyn(0) Trend je Dekade", 0.05, _t0[2], 0.05)
P("03", "Σ_dyn(0) Fehler des Trends", 0.025, _t0[3], 0.05)
P("03", "Σ_dyn(0) Trend mit Außenfeld 0,03·a₀", 0.02, _t1[2], 0.1)
P("08", "ε_L = ℓ_P/ℓ₀", 2.0e-62, l_P / (c**2 / a0L), 0.02)
P("08", "Dekaden 1/ε_L", 61.7, math.log10(c**2 / a0L / l_P), 0.001)
P("08", "M_min = ħa₀/c³ [kg]", 4.3e-70, hbar * a0L / c**3, 0.005)
import sparc_universal as _su
from scipy.stats import spearmanr as _sp
_o = _su.reste(); _r = np.array([_o[n]['rest'] for n in _o])
_d = np.array([_o[n]['innen'] - _o[n]['aussen'] for n in _o]); _d = _d[np.isfinite(_d)]
P("03", "Universalität: Galaxien", 149, len(_r), 1e-9)
P("03", "Universalität: Streuung des Rests [dex]", 0.16, _r.std(), 0.03)
P("03", "Universalität: innen − außen [dex]", -0.02, _d.mean(), 0.1)
_o8 = _su.reste(sig=8.0); _gd = [n for n in _o8 if _o8[n]['fgas'] > 0.5]
P("03", "Universalität: ρ(Rest, SB₀) gasreich, σ = 8 km/s", 0.36, _sp([_o8[n]['lSB'] for n in _gd], [_o8[n]['rest'] for n in _gd])[0], 0.01)
_nv = [n for n in _o if _o[n]['lVf'] == _o[n]['lVf']]
P("03", "Universalität: ρ(Rest, V_flach) ohne Druck", 0.16, _sp([_o[n]['lVf'] for n in _nv], [_o[n]['rest'] for n in _nv])[0], 0.03)
_o10 = _su.reste(sig=10.0); _nv = [n for n in _o10 if _o10[n]['lVf'] == _o10[n]['lVf']]
P("03", "Universalität: |ρ(Rest, V_flach)| mit σ = 10 km/s ≤ 0,01", 0.01, max(0.01, abs(_sp([_o10[n]['lVf'] for n in _nv], [_o10[n]['rest'] for n in _nv])[0])), 1e-9)
P("08", "Faltungstakt c/(3Δs) von M_U = a₀/c", 1.0, c / (3 * ds(c**4 / (G * a0L))) / (a0L / c), 1e-9)
P("09", "Spin ½: j = ½(m_Pl/M)² = 1 bei M = m_Pl/√2", 1.0, (hbar / 2) / ((m_Pl / math.sqrt(2)) * c) / (3 * ds(m_Pl / math.sqrt(2))), 1e-9)
import sparc_boseform as _bf
from scipy.optimize import minimize_scalar as _ms
_g, _gN = _bf.punkte(0.5); _B = _bf.bins(_gN, np.log10(_g / _gN), -2.0, math.log10(3), 14)
def _chi(f):
    r = _ms(lambda la: np.sum(((_B[:, 1] - np.log10(f(10**_B[:, 0] * A0, 10**la * 1e-10))) / _B[:, 3])**2), bounds=(-0.5, 0.5), method='bounded')
    return 10**r.x, r.fun
_aB, _cB = _chi(lambda x, a: nu(x, a)); _aS, _cS = _chi(_bf.nu_standard)
P("03", "a₀ an den SPARC-Punkten (g_N ≤ 3a₀) [10⁻¹⁰]", 1.08, _aB, 0.005)
P("07", "Standard-Form: χ²-Verhältnis zur Bose-Form (Υ_d = 0,5) im Bereich 3–8", 5.5, _cS / _cB, 0.02)
import sparc_streuung as _ss
_g0, _s0 = _ss.main(n_mc=200, dups=0.0); _g1, _s1 = _ss.main(n_mc=200, dups=0.1)
P("03", "Streuung des Rests mit Gasdruck [dex]", 0.14, _g0.std(), 0.03)
P("03", "aus Abstand und Neigung erwartet [dex]", 0.12, math.sqrt(np.mean(_s0**2)), 0.03)
P("03", "mit 0,1 dex Υ-Streuung erwartet [dex]", 0.13, math.sqrt(np.mean(_s1**2)), 0.03)
P("03", "Streuung der Regel selbst (Υ 0,1) [dex]", 0.05, math.sqrt(max(_g1.var() - np.mean(_s1**2), 0)), 0.03)
_ga, _sa = _ss.schwelle(n_mc=100, dups=0.0)
P("03", "a₀ je Galaxie direkt, Median [10⁻¹⁰]", 1.09, 10**np.median(_ga) * 1e10, 0.01)
P("03", "a₀ je Galaxie, robuste Streuung [dex]", 0.28, 1.4826 * np.median(np.abs(_ga - np.median(_ga))), 0.02)
_k, _ek, _ = _ss.anstieg(n_boot=500)
P("03", "Anstieg der Schwelle je Dekade Σ₀", 0.06, _k, 0.1)
P("03", "Fehler des Anstiegs", 0.04, _ek, 0.15)
P("08", "Λ_s·ℓ_P²", 1.0e-122, lambda_vakuum(schliesslaenge(a0L)) * l_P**2, 0.02)
_Mw = 6e10 * M_sun; _rw = r_stern(_Mw, a0L); _Ow = v_kreis_punkt(_Mw, _rw, a0L) / _rw
import tunnel_aussenfeld as _ta
def _m1(Rk, ref_nu=True):
    R = Rk * _ta.kpc; gi = G * _ta.M / R**2; ge = 0.03 * a0L; ph = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    eR = np.stack([np.cos(ph), np.sin(ph), 0 * ph]); gv = np.array([ge, 0, 0]); gt = -gi * eR + gv[:, None]
    ref = nu(ge, a0L) * gv if ref_nu else gv
    a = nu(np.linalg.norm(gt, axis=0), a0L) * gt - ref[:, None]; v = np.sqrt(np.clip(-R * np.sum(a * eR, axis=0), 0, None)); F = np.fft.rfft(v) / len(v)
    return 2 * abs(F[1]) / F[0].real * 100
P("03", "Einseitigkeit m = 1, 20 kpc, Bezug ν(g_e)g_e [%]", 13, _m1(20), 0.01)
P("03", "Einseitigkeit m = 1, 50 kpc, Bezug ν(g_e)g_e [%]", 20, _m1(50), 0.02)
P("03", "Einseitigkeit m = 1, 20 kpc, Bezug g_e [%]", 2, _m1(20, False), 0.02)
P("03", "Einseitigkeit m = 1, 50 kpc, Bezug g_e [%]", 36, _m1(50, False), 0.01)
# ---- Nachträge zu 02, 03, 06, 07, 08 (Zahlen der jüngsten Einträge)
P("08", "Proton: M/M_min = ω_C·τ_D", 3.9e42, (m_p * c**2 / hbar) * tau_D(a0L), 0.005)
P("08", "M_min = ħ/(c²τ_D) = ε_L·m_Pl", 1.0, (hbar / (c**2 * tau_D(a0L))) / (l_P / (c**2 / a0L) * m_Pl), 1e-9)
P("08", "M_U = c³τ_D/G = m_Pl/ε_L", 1.0, (c**3 * tau_D(a0L) / G) / (m_Pl / (l_P / (c**2 / a0L))), 1e-9)
P("07", "Abweichung bei r = (2/π)Λ: 1/tanh(π/2) − 1 [%]", 9.0, (1 / math.tanh(math.pi / 2) - 1) * 100, 0.005)
P("03", "Universalität: Fehler innen − außen [dex]", 0.013, np.std(_d) / math.sqrt(len(_d)), 0.05)
_gA, _sA = _ss.schwelle(n_mc=100, dups=0.0); _gB, _sB = _ss.schwelle(n_mc=100, dups=0.1)
P("03", "Schwelle je Galaxie: aus Fehlern erwartet, untere Grenze [dex]", 0.28, np.median(_sA), 0.02)
P("03", "Schwelle je Galaxie: aus Fehlern erwartet, obere Grenze [dex]", 0.31, np.median(_sB), 0.02)
_o3 = _su.reste(sig=8.0, ge=0.03); _g3 = [n for n in _o3 if _o3[n]['fgas'] > 0.5]
P("03", "Außenfeld 0,03·a₀ senkrecht: Rangkorrelation Rest–SB₀ gasreich", 0.23, _sp([_o3[n]["lSB"] for n in _g3], [_o3[n]["rest"] for n in _g3])[0], 0.01)
from scipy.special import i0 as _i0, i1 as _i1, k0 as _k0, k1 as _k1
from scipy.integrate import quad as _quad
_fy = lambda y: y * (_i0(y) * _k0(y) - _i1(y) * _k1(y)) if y < 300 else 1 / (8 * y * y)
def _sdyn(S0):
    gs = 2 * math.pi * G * S0 * M_sun / pc**2
    return S0 * (_quad(lambda y: float(nu(gs * _fy(y), a0L)) * _fy(y) / y, 0, 50, limit=400)[0] + _quad(lambda y: float(nu(gs / (8 * y * y), a0L)) / (8 * y**3), 50, 5e4, limit=200)[0])
_lel = lambda S: 1468 * (1 + S / 1571) ** 0.38 * (S / 1571) ** 0.62
P("03", "ideale Exponentialscheibe gegen Lelli-Verlauf, Σ₀ = 1000 [dex]", 0.08, math.log10(_sdyn(1000) / _lel(1000)), 0.05)
P("03", "ideale Exponentialscheibe gegen Lelli-Verlauf, Σ₀ = 10 [dex]", 0.22, math.log10(_sdyn(10) / _lel(10)), 0.05)
_kt = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kupfer_selbstpulsieren.txt"), encoding="utf-8").read()
_rows = [l.split("|") for l in _kt.split("3) Art der Hopf-Grenze")[1].split("4) Gegenprobe")[0].splitlines() if l.count("|") >= 6 and l.split("|")[0].strip().isdigit()]
_gap = [abs(float(r[2]) - float(r[3])) for r in _rows]; _ring = [float(r[4]) for r in _rows]; _qm = [float(r[0]) for r in _rows]
P("06", "Kupfer: Hopf vor der Falte, kleinster Abstand [Hz]", 0.03, min(_gap), 0.2)
P("06", "Kupfer: Hopf vor der Falte, größter Abstand [Hz]", 0.6, max(_gap), 0.1)
P("06", "Kupfer: Ringfrequenz untere Grenze [Hz]", 0.7, min(_ring), 0.1)
P("06", "Kupfer: Ringfrequenz obere Grenze [Hz]", 1.5, max(_ring), 0.01)
P("06", "Kupfer: kleinste Güte mit Hopf", 3000, min(_qm), 1e-9)
# ---- Atom als Ebene der Triade, zwei Klassen, Q je Radiant
_Tg = GM_sun / (0.93 * pc)**3
P("02", "Schwellenradius des Protons [fm]", 32, r_stern(m_p, a0L) * 1e15, 0.01)
P("02", "a_B / r*_Proton", 1660, a_B / r_stern(m_p, a0L), 0.002)
P("02", "Glättungslänge H-Atom in der galaktischen Gezeit [mm]", 2.7, (G * m_p / _Tg)**(1 / 3) * 1e3, 0.01)
P("08", "r*/R Proton", 38, r_stern(m_p, a0L) / R_p, 0.01)
P("08", "r*/R Erde", 3e5, r_stern(M_earth, a0L) / R_earth, 0.01)
P("08", "r*/R Sonne", 1.6e6, r_stern(M_sun, a0L) / R_sun, 0.02)
P("08", "r*/R Milchstraße (15 kpc)", 0.6, r_stern(6e10 * M_sun, a0L) / (15e3 * pc), 0.04)
P("08", "r*/R Haufen (1 Mpc)", 0.5, r_stern(1.7e14 * M_sun, a0L) / (1e6 * pc), 0.08)
_Qp = 0.1200 / 0.02237; _eQp = _Qp * math.hypot(0.0012 / 0.1200, 0.00015 / 0.02237)       # Planck 2018 [G]
_Qb = 0.11933 / 0.02242; _eQb = _Qb * math.hypot(0.00091 / 0.11933, 0.00014 / 0.02242)    # Planck 2018 + BAO [G]
P("03", "Q = 2π − 1", 5.283, 2 * math.pi - 1, 1e-4)
P("03", "Q Planck mit BAO [G]", 5.32, _Qb, 0.002)
P("03", "Fehler Q Planck mit BAO [G]", 0.05, _eQb, 0.06)
P("03", "Q Planck allein [G]", 5.36, _Qp, 0.002)
P("03", "Fehler Q Planck allein [G]", 0.07, _eQp, 0.08)
P("03", "Abstand 2π − 1 (Planck mit BAO) [σ]", 0.74, (_Qb - 2 * math.pi + 1) / _eQb, 0.07)
P("03", "Abstand 2π − 1 (Planck allein) [σ]", 1.3, (_Qp - 2 * math.pi + 1) / _eQp, 0.04)
P("03", "Ω_b/Ω_m = 1/2π", 0.159, 1 / (2 * math.pi), 0.002)
P("02", "m_p/m_e bei 𝒞 = 1", 1402, (m_p / m_e) / 2.944**0.25, 0.002)
_g4, _gN4 = _bf.punkte(0.4); _B4 = _bf.bins(_gN4, np.log10(_g4 / _gN4), -2.0, math.log10(3), 14)
_g6, _gN6 = _bf.punkte(0.6); _B6 = _bf.bins(_gN6, np.log10(_g6 / _gN6), -2.0, math.log10(3), 14)
def _chiB(B, f):
    return _ms(lambda la: np.sum(((B[:, 1] - np.log10(f(10**B[:, 0] * A0, 10**la * 1e-10))) / B[:, 3])**2), bounds=(-0.5, 0.5), method='bounded').fun
P("07", "Standardform: Fehlerquadratsumme-Verhältnis bei Υ_d = 0,4", 3.3, _chiB(_B4, _bf.nu_standard) / _chiB(_B4, lambda x, a: nu(x, a)), 0.03)
P("07", "Standardform: Fehlerquadratsumme-Verhältnis bei Υ_d = 0,6", 7.6, _chiB(_B6, _bf.nu_standard) / _chiB(_B6, lambda x, a: nu(x, a)), 0.03)
_o2 = _su.reste(sig=8.0, ge=0.02); _g2 = [n for n in _o2 if _o2[n]['fgas'] > 0.5]
P("03", "Außenfeld 0,02·a₀ senkrecht: Rangkorrelation Rest–SB₀ gasreich", 0.28, _sp([_o2[n]['lSB'] for n in _g2], [_o2[n]['rest'] for n in _g2])[0], 0.01)
_r3 = _sp([_o3[n]['lSB'] for n in _g3], [_o3[n]['rest'] for n in _g3])
P("03", "Außenfeld 0,03·a₀ senkrecht: p", 0.07, _r3[1], 0.08)
_ak = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "aussenfeld_kosmisch.txt"), encoding="utf-8").read()
import re as _re
P("03", "gemeinsames g_e frei, senkrecht [a₀]", 0.010, float(_re.search(r"senkrecht s⁴=s_i⁴\+s_e⁴\s*: A g_e/a₀ = ([\d.]+)", _ak).group(1)), 0.05)
P("03", "gemeinsames g_e frei, gleichgerichtet [a₀]", 0.003, float(_re.search(r"gleichgerichtet s²=s_i²\+s_e²\s*: A g_e/a₀ = ([\d.]+)", _ak).group(1)), 0.2)
_Lz = zelle_aus_a0(a0L); _Sb = G * rho_m_mitbewegt() / (1 + Q_KOSM) * _Lz**3 / (3 * c**2)
P("03", "Λ = ν_kosm·ΣΔs_b der Zelle, ν_kosm = 1 + Q", 1.0, schliesslaenge(a0L) / _Sb / (1 + Q_KOSM), 1e-9)
P("03", "Rückstellbilanz mit Rückläufen ν_kosm·ΣΔs/R", 0.53, ds(7.5e22 * M_sun) / (46.5e9 * ly) * (1 + Q_KOSM), 0.01)
_SH = math.pi * c**2 / (l_P**2 * H0**2); _SL = 45 / math.pi**3 * (schliesslaenge(a0L) / l_P)**2
P("03", "S_H/S_Λ heute", 0.24, _SH / _SL, 0.02)
P("03", "Horizontwachstum aus Materieverdünnung / Ṡ_i heute", 6.6e17, _SH * H0 * 3 * OMEGA_M / (7.1e86), 0.06)
from scipy.optimize import brentq as _bq
_OL = 1 - OMEGA_M - OMEGA_R
_af = _bq(lambda la: math.log(-(-1.5 * H0**2 * OMEGA_M * (10**la)**-3) / (l_P**2 * (H0 * math.sqrt(OMEGA_M * (10**la)**-3 + _OL))**3 * 7.1e86 / (2 * math.pi * c**2))), 3, 9)
_tf = 2 / (3 * H0 * math.sqrt(_OL)) * math.asinh(math.sqrt(_OL / OMEGA_M) * (10**_af)**1.5) / yr
P("04", "Bedingung −Ḣ verletzt ab [Jahre]", 2.5e11, _tf, 0.03)
_wp = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gewicht_rueckkehr.txt"), encoding="utf-8").read().split("\n")
_wr = {l.split()[0]: float(l.split()[4]) for l in _wp[101:106] if len(l.split()) >= 6 and l.split()[0] in ("3.4", "2.0", "1.0", "0.5")}
for _zc, _mw in (("3.4", 14), ("2.0", 45), ("1.0", 143), ("0.5", 309)):
    P("03", f"w von 0,9 auf 0,1, Kollaps bei z = {_zc} [Mio. Jahre]", _mw, _wr[_zc], 0.04)
P("03", "Preisach: 1 − w im Hintergrund gleich der jetzigen Form", 1.0, (1 - gewicht_w_preisach([6 * H0, 3 * H0], a0L)) / (1 - gewicht_w(3 * H0, a0L)), 1e-12)
_Om = 220e3 / (8 * kpc)
_sk = math.sqrt(G * 6e10 * M_sun / (8 * kpc) ** 2 / a0L)
_zk = complex(a0L * _sk ** 2 / c, _Om) * (c / a0L)
P("07", "Kern bei Bahnfrequenz |K̂(iΩ)|/K̂(0), Sonnenkreis", 2.1e-15, math.exp(-(_zk ** 0.5).real + _sk), 0.1)
import ringe_oertlich as _ro
for _Rk, _i, _w in ((12, 2, 1.42), (12, 89, 1.66), (20, 2, 1.69), (20, 89, 1.91), (30, 2, 2.00), (30, 89, 2.19)):
    P("07", f"Ring örtlich/Newton, {_Rk} kpc, {_i}°", _w, _ro.praez(_Rk * _ro.kpc, _i, True)[0] / _ro.praez(_Rk * _ro.kpc, _i, False)[0], 0.005)
P("07", "Ring i → 0: ω_p örtlich/Newton = √ν (20 kpc)", 1.0, (_ro.omega_kreis(20 * _ro.kpc, True) - _ro.nu_z(20 * _ro.kpc, True)) / (_ro.omega_kreis(20 * _ro.kpc, False) - _ro.nu_z(20 * _ro.kpc, False)) / math.sqrt(nu(-_ro.gN(np.array([20 * _ro.kpc, 0, 0]))[0])), 1e-4)
P("07", "Umlauf an der Glättungslänge der Sonne 2π/√T [Mio. Jahre]", 84, 2 * math.pi * math.sqrt((0.93 * pc)**3 / GM_sun) / yr / 1e6, 0.005)
T_gal = 4 * 15.3 * (15.3 + 11.9) * (1e3 / kpc)**2   # Oort A = 15,3, B = −11,9 km/s/kpc [G]
P("02", "Gezeitengrenze Sonne [pc]", 1.37, (GM_sun / T_gal)**(1 / 3) / pc, 0.01)
P("02", "Carter-Faktor 𝒞 = α_G/(α¹²(m_e/m_p)⁴)", 2.944, alpha_G() / (alpha**12 * (m_e / m_p)**4), 0.002)
P("02", "Eddington-Grenze Sonne [W]", 1.257e31, eddington(M_sun), 0.005)
P("02", "δ Sonnenoberfläche / α²", 0.04, delta(M_sun, R_sun) / alpha**2, 0.01)
P("02", "α = r_e/λ̄_e", 1 / 137.036, r_e / compton_red(m_e), 0.0001)
P("02", "Δs·λ̄/ℓ_P² (Proton)", 1 / 3, ds(m_p) * compton_red(m_p) / l_P**2, 1e-9)
P("02", "Im τ_G (α_G = e^{−2π Im τ})", 14.0096, -math.log(alpha_G()) / (2 * math.pi), 0.0001)
P("02", "Im τ_α", 0.7831, -math.log(alpha) / (2 * math.pi), 0.0002)
P("02", "magnetische/elektrische Elementarladung 1/(2α)", 68.5, 1 / (2 * alpha), 0.001)
# ---- 03
M_gal = 6e10 * M_sun
def v_im_haufen(r, D, M_cl=1e14 * M_sun, a=3 * kpc):
    gi = G * hernquist_masse(M_gal, a, r) / r**2; ge = G * M_cl * min(1.0, D / Mpc) / D**2   # Haufen: M(<R) ∝ R bis 1 Mpc
    return math.sqrt(float(nu(math.hypot(gi, ge))) * gi * r)
def v_frei(r, a=3 * kpc):
    gi = G * hernquist_masse(M_gal, a, r) / r**2
    return math.sqrt(float(nu(gi)) * gi * r)
P("03", "Galaxie im Haufen, 80 kpc, 0,5 Mpc: Anteil des freien Werts", 0.50, v_im_haufen(80 * kpc, 0.5 * Mpc) / v_frei(80 * kpc), 0.01)
P("03", "Galaxie im Haufen, 80 kpc, 2 Mpc: Anteil des freien Werts", 0.77, v_im_haufen(80 * kpc, 2 * Mpc) / v_frei(80 * kpc), 0.01)
a0_1090 = A0 * a0_frueh_relativ(1090)
P("03", "1 − w bei der Rekombination: Exponent (e^{−738})", 738, math.sqrt(3 * H0 * E_hintergrund(1090) / (a0_1090 / c)), 0.002)
def theta_zahl(n, hex_=False):
    cnt = 0
    for m in range(-5, 6):
        for k in range(-5, 6):
            if (m*m + m*k + k*k if hex_ else m*m + k*k) == n: cnt += 1
    return cnt
for n, v in zip(range(1, 6), (4, 4, 0, 4, 8)):
    P("03", f"Bilder quadratischer Torus, Abstand² = {n}", v if v else 1e-99, theta_zahl(n) if v else 1e-99, 1e-9)
for n, v in zip(range(1, 6), (6, 0, 6, 6, 0)):
    P("03", f"Bilder sechseckiger Torus, Abstand² = {n}", v if v else 1e-99, theta_zahl(n, True) if v else 1e-99, 1e-9)
P("03", "Zelle mindestens (0,985·d_LSS) [Gpc]", 27.4, KREISPAAR_GRENZE * d_LSS() / Gpc, 0.005)
# ---- 04
P("04", "Takt bei 100 r₀", 0.995, math.sqrt(1 - 1 / 100), 0.001)
P("04", "Takt bei 2 r₀", 0.707, math.sqrt(0.5), 0.001)
P("04", "S = 36πk(Δs/ℓ_P)² = 4πk·α_G (Sonne)", 1.0, 36 * math.pi * (ds(M_sun) / l_P)**2 / (4 * math.pi * alpha_G(M_sun)), 1e-9)
P("04", "Lebensdauer 3840(S/k)Δs/c = Hawking", 1.0, lebensdauer_falte(M_sun) / hawking_lebensdauer(M_sun), 1e-9)
P("04", "Temperatur global = ħH/(2πk) (Δs_H = c/6H)", 1.0, (hbar * c / (12 * math.pi * kB * c / (6 * H0))) / (hbar * H0 / (2 * math.pi * kB)), 1e-9)
L_sun_abs = 1.22e17                                  # von der Erde aufgenommene Leistung [G]
P("04", "Entropieexport der Erde [W/K]", 6.1e14, 4 / 3 * L_sun_abs * (1 / 255 - 1 / 5772), 0.01)
# ---- 05
P("05", "ν bei Im τ = ½", 1.045, 1 / (1 - math.exp(-math.pi)), 0.001)
P("05", "ν bei Im τ = 2", 1.0000035, 1 / (1 - math.exp(-4 * math.pi)), 1e-7)
for k_, v in ((2, 24), (4, 240), (6, 504)):
    from scipy.special import zeta, gamma
    P("05", f"(2π)^k/(Γ(k)ζ(k)), k = {k_}", v, (2 * math.pi)**k_ / (gamma(k_) * zeta(k_)), 1e-9)
P("05", "12 = 2/(1 − ½ − ⅓)", 12, 2 / (1 - 0.5 - 1 / 3), 1e-9)
P("05", "1728 = 3·240 + 2·504", 1728, 3 * 240 + 2 * 504, 1e-9)
x_um = brentq(lambda x: -1 / x + 0.005 * (x - 1)**2 - (0.5 * 1.6**2 - 1), 1.5, 50)
P("05", "Gedächtnis-Bindung: Umkehrpunkt (v = 1,6, γ = 0,01)", 9.75, x_um, 0.002)
# ---- 06
P("06", "Aufschaukeln im Dreieck: (Σt)(Σ1/t) − 1 ≥ 8", 8, 3 * 3 - 1, 1e-9)


def _tau_verh(Qm):
    d = dict(DRAHT); d['Qm'] = Qm
    p = draht(CU, d)
    Ql = p['m'] * p['w0'] / p['c_m'] / (1 + p['Gam']**2 / p['c_m'] / p['R0'])
    return (p['C_th'] / p['G_th']) / (2 * Ql / p['w0'])


_M = 6e10 * M_sun; _r = r_stern(_M, a0L); _v = v_kreis_punkt(_M, _r, a0L)
P("05", "x = 3H₀τ_D heute", 17.9, 3 * H0 * tau_D(a0L), 0.005)
P("06", "Kupferdraht τ_th/τ_m bei Q_m = 1000", 1.7, _tau_verh(1000), 0.01)
P("06", "Kupferdraht τ_th/τ_m bei Q_m = 3000", 0.74, _tau_verh(3000), 0.01)
# ---- 07
P("07", "Pólya: ν ≤ 1/(1 − 0,3405) in drei Dimensionen", 1.516, 1 / (1 - 0.340537), 0.001)
P("07", "beide Sinne voll: (1 + 2n_B)/(1 + n_B) bei s = 1", 1.37, (1 + 2 * float(n_B(1.0))) / (1 + float(n_B(1.0))), 0.002)
P("07", "Milchstraße bei r* mit beiden Sinnen [km/s]", 253, (G * MW * a0L)**0.25 * math.sqrt(float(nu(a0L, a0L))) / 1e3 * math.sqrt((1 + 2 * float(n_B(1.0))) / (1 + float(n_B(1.0)))), 0.003)
P("07", "Periode in Im τ: Sonnenoberfläche P >", 2.5e5, r_stern(M_sun) / R_sun / (2 * math.pi), 0.01)
P("07", "220 km/s gehören zu [M☉]", 1.60e11, (2.2e5)**4 / (G * A0) / M_sun, 0.003)
rs_MW = r_stern(MW)
for pq, rk in ((6 / 5, 5.9), (5 / 4, 8.2), (9 / 7, 10.9), (4 / 3, 18.0), (7 / 5, 108)):
    P("07", f"Resonanz Re τ_M = {pq:.4f} in der Milchstraße [kpc]", rk, radius_resonanz(pq, rs_MW) / kpc, 0.01)
P("07", "Re τ_M bei r* (Punktmasse)", 1.26, math.sqrt(realteil_quadrat_punkt(1.0)), 0.002)
for d_, v in ((1, math.pi**2 / 6), (2, 2 * 1.2020569), (3, math.pi**4 / 15)):
    P("07", f"Bose-Leiter: ∫ s^{d_} n_B ds", v, quad(lambda x: x**d_ / math.expm1(x), 0, 200)[0], 1e-6)
sm = bose_stufe_max(2)[0]
P("07", "größte Zusatzbeschleunigung Punktmasse bei r/r*", 0.628, 1 / sm, 0.002)
P("07", "größte Zusatzbeschleunigung Milchstraße [kpc]", 5.5, rs_MW / sm / kpc, 0.01)
x3 = 3.0
P("07", "M_schein/M Reihe bei r = 3r*", 1 / (-math.expm1(-1 / x3)), x3 + 0.5 + 1 / (12 * x3) - (1 / x3)**3 / 720, 1e-5)
# ---- 08
P("08", "r* Sonne [Lj]", 0.116, r_stern(M_sun) / ly, 0.005)
P("08", "δ(r*) = (v_flach/c)² (Milchstraße)", 1.0, (3 * ds(MW) / rs_MW) / (v_flach(MW) / c)**2, 1e-9)
P("08", "ISCO j = 0,9 mitläufig [3Δs]", 2.32, isco_kerr(0.9), 0.002)
P("08", "ISCO j = 0,9 gegenläufig [3Δs]", 8.72, isco_kerr(0.9, False), 0.002)
a_gpb = 7027.4e3; w_gpb = math.sqrt(GM_earth / a_gpb**3)
P("08", "Gravity Probe B geodätisch [mas/Jahr]", 6606, 4.5 * ds(M_earth) / a_gpb * w_gpb * yr * ARC * 1e3, 0.005)
P("08", "Gravity Probe B Mitführung [mas/Jahr]", 39.2, c * 3 * ds(M_earth) * aJ / (2 * a_gpb**3) * yr * ARC * 1e3, 0.05)
P("08", "Faltungsmesser GW170817 Δs_c [m]", 589, ds(1.1977 * M_sun), 0.002)
P("08", "Schwellenfrequenz Gravitationswellen (1 M☉) [Hz]", 1e-13, schwellenfrequenz_gw(M_sun), 0.15)
P("08", "ISS [s/Jahr]", 0.021, 3 * ds(M_earth) / (R_earth + 400e3) * yr, 0.02)


# ---- Nachtrag Gesamtdurchgang
lqk = -math.log(q_KOSM)
M_Q0 = (1 + Q_KOSM)**3 * A0**3 * lqk**6 / (1e4 * H0**4) / G / M_sun
P("03", "M_Q heute (R_Q = R₂₀₀) [M☉]", 2.8e14, M_Q0, 0.015)
_xw = 3 * H0 * c / a0L
_ww = gewicht_w(3 * H0, a0L)
P("07", "(1 − w)(ν − 1) zwischen Strukturen heute", 2e-4, (1 - _ww) / math.expm1(math.sqrt(_xw)), 0.1)
P("07", "Verluste unter einer Wurzel: Γτ_D = 1, θτ_D = 17,9", 1.30e-2, math.exp(-math.sqrt(1 + 3 * H0 * c / a0L)), 0.01)
P("07", "D = r₀c/2 = GM/c (Sonne) [m²/s]", G * M_sun / c, 6 * G * M_sun / (3 * c**2) * c / 2, 1e-12)
P("07", "linear nichtlokal: v⁴/(GMa₀) = (2/π)²", 0.405, (2 / math.pi)**2, 0.002)
P("04", "Endrate bei festem Λ: π²c/(√45·Λ) [km/s/Mpc]", 33.3, math.pi**2 * c / (math.sqrt(45) * c**2 / (2 * a0L)) * 3.0857e19, 0.005)
P("03", "r*(M_b,Zelle)/Λ = √(6/ν_kosm)", 0.97, math.sqrt(6 / (1 + Q_KOSM)), 0.005)
_n1 = 1 / (1 - math.exp(-1))
P("07", "Fisher-Metrik bei s = 1: ν(ν − 1) gegen 1/(4 sinh²(½))", 1 / (4 * math.sinh(0.5)**2), _n1 * (_n1 - 1), 1e-9)
P("03", "(1 + Q)·a₀-Grenze der Würfelzelle [10⁻¹⁰ m/s²]", 7.17, a0_aus_zelle(0.985 * d_LSS()) * (1 + Q_KOSM) * 1e10, 0.002)
P("03", "Würfel: a₀ ≤ bei Q = 5,36 [10⁻¹⁰]", 1.13, a0_aus_zelle(0.985 * d_LSS()) * 1e10, 0.005)
P("03", "Fluchtmaß ganze Tori gegen Rand: 16(3/2 − √2)", 1.3726, 16 * (1.5 - math.sqrt(2)), 1e-4)
P("03", "Zellgrenze im eigenen Atmen (1+Q)a₀ ≤ [10⁻¹⁰]", 6.74, a0_aus_zelle(0.985 * d_LSS()) * (1 + Q_KOSM) * 1e10 / 1.0209**3, 0.002)
P("03", "Würfel: Q ≤ bei a₀ = 1,1 (eigenes Atmen)", 5.13, a0_aus_zelle(0.985 * d_LSS()) * (1 + Q_KOSM) / 1.0209**3 / 1.1e-10 - 1, 0.003)
P("03", "M_b,Q heute (Baryonen, R_Q = R₂₀₀) [M☉]", 4.5e13, M_Q0 / (1 + Q_KOSM), 0.015)
P("09", "Dephasierungszeit ħ/kT bei 310 K [s]", 2.46e-14, hbar / (kB * 310), 0.002)
P("09", "Debye-Länge physiologische Salzlösung 0,15 M [nm]", 0.78, 0.304 / math.sqrt(0.15), 0.01)

# ---- Nachtrag Durchgang der Forschung (Notizen 01–75)
from scipy.integrate import quad as _quad
P("03", "ν bei g/a₀ = 25 (kompakter Kugelsternhaufen, g = 2,8·10⁻⁹)", 1.007, nu(2.76e-9, a0L), 0.001)
P("03", "ν bei g/a₀ = 0,13 (Zwerggalaxie, g = 1,44·10⁻¹¹)", 3.3, nu(1.44e-11, a0L), 0.01)
P("03", "σ Tiefgrenze 10⁷ M☉ [km/s]", 9.2, (4 / 81 * G * 1e7 * M_sun * a0L) ** 0.25 / 1e3, 0.005)
P("03", "σ Tiefgrenze 2·10⁸ M☉ [km/s]", 19.5, (4 / 81 * G * 2e8 * M_sun * a0L) ** 0.25 / 1e3, 0.005)
def _haufen(r):   # Notiz 28: β-Gas 10¹⁴ M☉ in 1 Mpc (r_c 150 kpc, bis 3 Mpc), BCG 1,5·10¹² (20 kpc), Galaxien 10¹³ (300 kpc); NFW 10¹⁵, c = 4 bei z = 0,3
    f = lambda x: x * x / (1 + (x / (150 * kpc))**2)
    Mg = 1e14 * M_sun * _quad(f, 0, min(r, 3 * Mpc))[0] / _quad(f, 0, Mpc)[0]
    Mb = Mg + hernquist_masse(1.5e12 * M_sun, 20 * kpc, r) + hernquist_masse(1e13 * M_sun, 300 * kpc, r)
    rc = 3 * (H0 * E_hintergrund(0.3))**2 / (8 * math.pi * G)
    r200 = (3 * 1e15 * M_sun / (4 * math.pi * 200 * rc)) ** (1 / 3)
    mu = lambda x: math.log(1 + x) - x / (1 + x)
    return 1e15 * M_sun * mu(4 * r / r200) / mu(4) / (nu(G * Mb / r**2, a0L) * Mb)
P("03", "Haufen: fehlender Faktor bei 100 kpc", 3.1, _haufen(100 * kpc), 0.02)
P("03", "Haufen: fehlender Faktor bei 1 Mpc", 1.6, _haufen(Mpc), 0.03)
P("03", "kosmische Streuung des Quadrupols", 0.63, math.sqrt(2 / 5), 0.005)
_rate, _S, _Hend = kreislauf_mindestrate(schliesslaenge(a0L))
P("03", "w_s + 1 heute aus der Mindestrate", 1e-19, 2 / 3 * _rate / H0, 0.2)
P("03", "Warp: 90° bei 10,86 km/s/kpc [Mrd. Jahre]", 0.14, math.pi / 2 / (10.86e3 / kpc) / (1e9 * yr), 0.03)
_tD = tau_D(a0L)
def _anteil(s, t):
    f = lambda u: (4 * math.pi * u**3) ** -0.5 * math.exp(-1 / (4 * u) - s * s * u)
    return _quad(f, 1e-6, t / _tD, limit=400)[0] / math.exp(-s)
P("07", "Rückkehr als Vorgang: Anteil nach 13,8 Mrd. J. (s = 1)", 0.19, _anteil(1.0, 13.8e9 * yr), 0.02)
P("07", "Rückkehr als Vorgang: Anteil nach 13,8 Mrd. J. (s = 0,3)", 0.103, _anteil(0.3, 13.8e9 * yr), 0.02)
P("07", "Rückkehrzeit Kugel r* im Hubble-Volumen, Sonne [Jahre]", 3e32, 4 / 3 * math.pi * (c / H0)**3 / (math.pi * r_stern(M_sun)**2 * c) / yr, 0.12)
P("07", "Periode in Im τ: Grenze Sonnenoberfläche", 2.5e5, math.sqrt(neigung(M_sun, R_sun) / a0L) / (2 * math.pi), 0.01)
P("08", "Knotenzahl an r* Milchstraße (j = 1329)", 5e-7, 2 * 1329 * (v_flach(MW) / c)**3, 0.08)
P("08", "Linse I(1)", 2.135, _quad(lambda t: math.cos(t) * nu(math.cos(t)**2 * a0L, a0L), 0, math.pi / 2)[0], 0.001)
P("02", "Radius des Ladungskreises 2ℓ_P/√α [ℓ_P]", 23.4, 2 / math.sqrt(alpha), 0.002)
P("02", "Umlauf um den Ladungskreis [GeV]", 5e17, hbar * c / (2 * l_P / math.sqrt(alpha)) / 1.602176634e-10, 0.05)
def _rueckgabe(M, OmDE):
    Hs = lambda a: H0 * math.sqrt(OMEGA_M / a**3 + OMEGA_R / a**4 + OmDE)
    an = T_CMB / (hbar * c**3 / (8 * math.pi * G * M * M_sun * kB))
    return _quad(lambda la: 1 / Hs(math.exp(la)), math.log(1e-8), math.log(an), limit=400)[0] / yr
P("04", "Einsatz der Rückgabe 10 M☉ (Dunkle Energie konstant) [Jahre]", 4e11, _rueckgabe(10, OMEGA_DE), 0.12)
P("04", "Einsatz der Rückgabe 10⁹ M☉ (nur Vakuumteil) [Jahre]", 1e12, _rueckgabe(1e9, K_S), 0.2)

def _wp(R, Om, wpN):   # örtliche Regel in der Ebene: Ω und ν_z mit √ν(g_R), also ω_p = √ν·ω_p,N (Newton-Frequenzen aus warp_praezession)
    return math.sqrt(float(nu((Om * 1e3 / kpc)**2 * R * kpc, a0L))) * wpN
P("03", "Warp: Aufwickeln 12–16 kpc, örtliche Regel [Mrd. Jahre]", 0.26, math.pi / (abs(_wp(12, 13.09, -32.40) - _wp(16, 8.64, -21.48)) * 1e3 / kpc * 1e9 * yr), 0.03)
P("03", "Warp: Aufwickeln 12–16 kpc, Newton [Mrd. Jahre]", 0.28, math.pi / ((32.40 - 21.48) * 1e3 / kpc * 1e9 * yr), 0.03)

_aN = G * MW / (8.2 * kpc)**2      # Pulsare am Sonnenort, Milchstraße als Punktmasse
P("03", "Pulsare: radial Ṗ/P, örtliche Regel [s⁻¹]", 6.3e-19, float(nu(_aN, a0L)) * _aN / c, 0.01)
P("03", "Pulsare: radial Ṗ/P, Newton [s⁻¹]", 4.1e-19, _aN / c, 0.02)

# ---- Verknüpfungen (Ordnung der Zahlen)
_Mc, _MU = masse_r_stern_gleich_compton(a0L), c**4 / (G * a0L)
P("08", "M_U = c⁴/(Ga₀) [kg]", 1.1e54, _MU, 0.01)
P("08", "M_c³·M_U/m_Pl⁴", 1.0, _Mc**3 * _MU / m_Pl**4, 1e-9)
P("08", "r*/λ̄ Proton = (m_p/M_c)^{3/2}", 151, r_stern(m_p, a0L) / compton_red(m_p), 0.005)
P("08", "r*/r₀ = 1/(2δ*) (Milchstraße)", 1.0, (r_stern(MW, a0L) / r0(MW)) * 2 * (v_flach(MW, a0L) / c)**2, 1e-9)
P("08", "f*_gw·π/Ω(r*)", 1.0, schwellenfrequenz_gw(M_sun, a0L) * math.pi / (v_flach(M_sun, a0L) / r_stern(M_sun, a0L)), 1e-9)
P("07", "a₀/c = Γ(r*) = 1/τ_D [s⁻¹]", 3.7e-19, (G * MW / r_stern(MW, a0L)**2) / c, 0.01)
P("07", "w bei θ = a₀/c = 1 − 1/e", 0.632, gewicht_w(a0L / c, a0L), 0.001)
P("07", "w = ½ bei θ = (ln 2)²·a₀/c", 0.5, gewicht_w(math.log(2)**2 * a0L / c, a0L), 1e-9)
_Th = 2 * math.pi * a0L / (c * H0)
P("03", "Θ = 2πa₀/(cH₀) heute", 1.06, _Th, 0.005)
P("03", "Ω_s/(1 − w) heute", 16.7, omega_vakuum(a0L) / (1 - gewicht_w(3 * H0, a0L)), 0.01)
P("03", "Zelle L·H₀/c = (8π²/(Ω_mΘ))^{1/3}", 6.2, zelle_aus_a0(a0L) * H0 / c, 0.005)
P("03", "Θ für L = d_LSS", 1.03, 8 * math.pi**2 / (OMEGA_M * (d_LSS() * H0 / c)**3), 0.005)
P("09", "Carter als Frequenzleiter Γ_G·ω₁³/(αω₂)⁴", 2.944, alpha_G() * (c / compton_red(m_p))**4 / (alpha * alpha * c / a_B)**4, 0.002)

# ---- Wicklung und Rückkehr (07): Drehung golden, Zelle μ = 0,01; Halbgerade 0 … N spiegelnd
_al = (math.sqrt(5) - 1) / 2; _ix = np.nonzero((np.arange(1, 400000) * _al) % 1 < 0.01)[0]; _Tw = np.diff(_ix)
P("07", "Drehung golden μ = 1 %: Werte 55, 89, 144 (Summe)", 288, float(np.unique(_Tw).sum()), 1e-9)
_Gw = lambda g: np.mean(np.exp(-g * _Tw))
P("07", "Drehung: Exponent von 1 − G bei Γ = 10⁻⁴ (→ 1)", 1.0, math.log((1 - _Gw(1.01e-4)) / (1 - _Gw(1e-4 / 1.01))) / (2 * math.log(1.01)), 0.01)
_N = 50; _Qh = np.diag(np.full(_N - 1, 0.5), 1) + np.diag(np.full(_N - 1, 0.5), -1); _Qh[-1, -2] = 1.0
P("07", "Kac Halbgerade: ⟨T₀⟩ = 2N (N = 50)", 100, 1 + np.linalg.solve(np.eye(_N) - _Qh, np.ones(_N))[0], 1e-9)

# ---- System als Ganzes im äußeren Feld (03): A = Massenmittel der örtlichen Regel
import re
_mu, _wmu = np.polynomial.legendre.leggauss(200)
def _Apl(X, ge):   # Plummer, GM/b² = X·a₀, g_e in a₀
    def f(r):
        gi = X * r / (1 + r * r)**1.5; gz = ge - gi * _mu; m = np.hypot(gz, gi * np.sqrt(1 - _mu**2))
        return 3 * r * r / (1 + r * r)**2.5 * 0.5 * np.sum(_wmu * nu(m * A0) * gz)
    return quad(f, 0, np.inf, limit=400)[0] / ge
P("03", "Plummer GM/b² = a₀, g_e = 0,1a₀: A/g_e", 2.37, _Apl(1, 0.1), 0.005)
P("03", "ν(0,1·a₀)", 3.69, float(nu(0.1 * A0)), 0.005)
_gg = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "galaxie_als_ganzes.txt"), encoding="utf-8").read()
P("03", "Exponentialscheibe MW, 0,03a₀: A/g_e", 1.27, float(re.search(r"Exponentialscheibe R_d = 3 kpc: A/g_e = ([\d.]+)", _gg).group(1)), 0.005)
_ex = _gg.split("Exponentialscheibe 3 kpc")[1]
_m1A = [float(x) for x in re.findall(r"R = \s*\d+ kpc.*\|\s*([\d.]+)\n", _ex)[:4]]
P("03", "m = 1 im Bezug A bis 30 kpc, höchstens [%]", 4.4, max(_m1A[:3]), 0.01)
P("03", "m = 1 im Bezug A bei 50 kpc [%]", 31, _m1A[3], 0.01)

# ---- Klassengrenzen und Spiegel des Trichters (07, 08); Gezeit am Sonnenort aus trichter_invarianten
_ti = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "trichter_invarianten.txt"), encoding="utf-8").read()
_Tg = float(re.search(r"T = √\(3/2\)\|E\^TF\| = ([\d.e+-]+)", _ti).group(1))
_Mx = a0L**3 / (G * _Tg**2)
P("07", "M× = a₀³/(GT²) am Sonnenort [M☉]", 3.3e8, _Mx / M_sun, 0.02)
P("07", "R× = a₀/T [pc]", 642, a0L / _Tg / pc, 0.01)
P("07", "ℓ/r* Sonne = (M×/M☉)^{1/6}", 26.2, (_Mx / M_sun)**(1 / 6), 0.005)

P("03", "S_Λ = S_H/Ω_s", 1.0, _SH / _SL / omega_vakuum(a0L), 1e-9)
P("07", "Milchstraße bei r* [km/s]", 216, (G * MW * a0L)**0.25 * math.sqrt(float(nu(a0L, a0L))) / 1e3, 0.003)

P("08", "Δs(M_U)/Λ", 2 / 3, G * (c**4 / (G * a0L)) / (3 * c**2) / schliesslaenge(a0L), 1e-9)

# ---- Schranken druckgestützter Systeme (03): σ_N und isolierter Wert aus gedaechtnis_umgebung
_gu = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gedaechtnis_umgebung.txt"), encoding="utf-8").read()
def _gsys(nm, key):
    return float(re.search(r"-- " + re.escape(nm) + r":.*?" + key + r" =? ?([\d.]+)", _gu).group(1))
for _nm, _v in (("Pal 14", 0.66), ("NGC 2419", 5.09), ("DF4", 7.15)):
    P("03", "σ_N " + _nm + " [km/s]", _v, _gsys(_nm, "σ_N"), 0.01)
P("03", "Draco isoliert [km/s]", 4.74, _gsys("Draco", "isoliert"), 0.01)

P("03", "Kreislauf: Λ̇ [m/s]", 1.6e-10, rate * schliesslaenge(a0L), 0.05)
_MH = OMEGA_B * c**3 / (2 * G * H0)          # Baryonen in der Hubble-Kugel; Faltung nur der Baryonen
P("04", "Σ S_Falten/S_H höchstens Ω_b² (36πk(Δs/ℓ_P)² mit ΣΔs_b)", OMEGA_B**2, 36 * math.pi * (G * _MH / (3 * c**2) / l_P)**2 / (math.pi * c**2 / (l_P**2 * H0**2)), 1e-9)

P("03", "χ bei R_Q = Q·(−ln q_kosm)", 0.917, float(chi(math.log(Q_KOSM / (1 + Q_KOSM))**2 * a0L, a0L)), 0.001)

# ---- Tiefe und Temperatur (02), innere Richtung als imaginäre Zeit (03, 07)
P("02", "T_c Sonne aus der Tiefe, n = 3, μ = 0,62 [K]", 1.2e7, 0.8543 * 3 * ds(M_sun) / 6.957e8 * 0.62 * 931.494e6 / 8.617333e-5, 0.02)
_ft = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fusion_trichter.txt"), encoding="utf-8").read()
P("02", "Grenzmasse aus Tiefe gegen Zählregel [M☉]", 1.456, float(re.search(r"n = 3 → ([\d.]+) M☉", _ft).group(1)), 0.001)
P("02", "Decke kT_max/(m_ec²) je μ/μ_e", 0.0964, float(re.search(r"\(μ/μ_e\)·([\d.]+)·m_ec² \(bei", _ft).group(1)), 0.001)
P("03", "T_I = ħa₀/(ck) [K]", 2.80e-30, hbar * a0L / (c * 1.380649e-23), 0.005)
P("03", "Θ künftig → √45/π", 2.135, math.sqrt(45) / math.pi, 0.001)

_ie = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "invarianten_ebenen.txt"), encoding="utf-8").read()
P("03", "L/d_LSS für Beobachter bei z = 1", 1.32, float(re.search(r"z_Beob =   1.0:.*?L/d_LSS = ([\d.]+)", _ie).group(1)), 0.005)
P("03", "L/d_LSS künftig", 0.73, float(re.search(r"L/d_LSS → ([\d.]+)", _ie).group(1)), 0.01)

P("07", "r*/Λ einer Masse mit ΣΔs_b = Λ/ν_kosm", 0.97, math.sqrt(6 / (1 + Q_KOSM)), 0.005)

P("03", "a_N/a₀(1090) bei 147 Mpc, Φ/c² = 10⁻⁵", 2.6, (1e-5 * c**2 / (147 * 1e6 * pc / 1091)) / (a0L * float(a0_frueh_relativ(1090))), 0.02)

_fg = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fruehe_galaxien_a0.txt"), encoding="utf-8").read()
P("03", "Preisach: Δlog M_b bei festem v (z_c = 3,4 → z = 2, 10 kpc) [dex]", 0.013, float(re.search(r"Δ = ([\d.]+) dex; N = ", _fg).group(1)), 0.03)
P("03", "Preisach: nötiger Massennullpunkt [dex]", 0.004, float(re.search(r"auf Δ/3 = ([\d.]+) dex", _fg).group(1)), 0.06)

P("03", "M_b tief bei kT = 1 keV, S = 2, μ = 0,6 [M☉]", 6.98e12, (2 * 1.602176634e-16 / (0.6 * m_p))**2 / (G * a0L) / M_sun, 0.005)
P("04", "δ Weißer Zwerg 0,6 M☉ (R = 8,696·10⁶ m)", 1.0e-4, 3 * ds(0.6 * M_sun) / 8.696e6, 0.03)

P("07", "ℓ Erde in galaktischer Gezeit [AE]", 2780, (G * M_earth / _Tg)**(1 / 3) / AU, 0.005)
P("07", "ℓ Mond in galaktischer Gezeit [AE]", 641, (G * 7.342e22 / _Tg)**(1 / 3) / AU, 0.005)

_tg = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "trichterkette_gewicht.txt"), encoding="utf-8").read()
P("03", "θ_c/(3H₀)", 0.056, a0L / (c * 3 * H0), 0.005)
P("03", "Umkehr heute: 1 + Δ_ta", 11.8, float(re.search(r"Umkehr heute \(θ = 0\): 1 \+ Δ_ta = ([\d.]+)", _tg).group(1)), 0.005)
P("03", "Einstein-de Sitter 1 + Δ_ta = 9π²/16", 5.55, 9 * math.pi**2 / 16, 0.001)

_sar = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sparc_aussenfeld_rand.txt"), encoding="utf-8").read()
P("03", "SPARC: Galaxien mit äußerstem a_N ≤ 0,05·a₀", 72, float(re.search(r"a_N ≤ 0,05·a₀: (\d+)", _sar).group(1)), 1e-9)
P("03", "SPARC: Reststeigung mit g_e = 0,03·a₀ senkrecht", -0.30, float(re.search(r"senkrecht 0.03 .*?Steigung (-?[\d.]+) ±", _sar).group(1)), 0.02)
P("03", "D564-8 letzter Punkt a_N/a₀ (von Hand 0,00756)", 0.00756, 78.82e6 / (3.07 * kpc) / a0L, 0.005)

_sd = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sparc_spektraldimension.txt"), encoding="utf-8").read()
P("07", "SPARC Spektraldimension d_s", 1.03, float(re.search(r"alle Punkte  d_s = ([\d.]+) ±", _sd).group(1)), 0.02)
P("03", "ν Milchstraße am Umkehrradius 0,58 Mpc", 67, float(nu(G * MW / (0.58e3 * kpc)**2, a0L)), 0.01)
P("03", "a₀-Grenze sechseckige Zelle", 1.30, 1.13 / (math.sqrt(3) / 2), 0.01)
P("03", "a₀-Grenze flächenzentrierte Zelle", 1.59, 1.127 * math.sqrt(2), 0.005)

P("03", "Filament r*_fil = 2Gμ/a₀ bei 10¹² M☉/Mpc [kpc]", 2.5, 2 * G * 1e12 * M_sun / (1e6 * pc) / a0L / kpc, 0.02)
P("03", "d_eff Kugel bei χ = 1", 2.0, 1 + 2 * (1 - 0.5), 1e-9)
P("01", "Leiter: (3Δs/r*)/(r*/ℓ₀) = 1 (Milchstraße)", 1.0, (3 * ds(MW) / r_stern(MW, a0L)) / (r_stern(MW, a0L) / (c**2 / a0L)), 1e-9)

P("07", "Glättungsrate D·k² = a_N/c (D = GM/c, k = 1/r; Sonne bei 1 AE)", 1.0, (G * M_sun / c) / AU**2 / (G * M_sun / AU**2 / c), 1e-12)
P("03", "ohne Gefälle ν = 1/w (θ̃ = 3H₀)", 1.0, gewicht_w(3 * H0, a0L) / (1 - math.exp(-math.sqrt(3 * H0 / (a0L / c)))), 1e-9)

_nr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "naht_rq_umkehr.txt"), encoding="utf-8").read()
P("03", "Minimum a/a_N an der w-Schale 2√(1+Q) − 1", 4.05, 2 * math.sqrt(1 + Q_KOSM) - 1, 0.002)
P("03", "R_Q/r* ≈ Q + ½ (tief)", 5.85, 1 / -math.log(Q_KOSM / (1 + Q_KOSM)), 0.005)
P("03", "SPARC jenseits R_Q: Rest ohne Kappe [dex]", -0.02, float(re.search(r"ohne Kappe, σ = 8 km/s  : Rest jenseits R_Q (-?[\d.]+)", _nr).group(1)), 0.2)

# ---- Ausgabe
bad = 0
print(f"{'Baustein':8s} | {'Aussage':52s} | {'Modell':>11s} | {'Kern':>11s} | {'Abw.':>7s} |")
for b, t, m, k, tol in checks:
    d = abs(k - m) / abs(m)
    ok = d <= tol
    bad += (not ok)
    print(f"{b:8s} | {t:52s} | {m:11.4g} | {k:11.4g} | {d*100:6.2f} % | {'OK' if ok else 'ABWEICHUNG'}")
print(f"\n{len(checks)} Zahlen geprüft, {bad} Abweichung(en).")
