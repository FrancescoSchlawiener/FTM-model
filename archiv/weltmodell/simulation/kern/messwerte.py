# -*- coding: utf-8 -*-
"""
Messwerte des Modells (Klasse und Regeln in messwert.py). Neue Messwerte kommen hierher, je Quelle eine Zeile;
konstanten.py liest die Kernwerte von hier.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
CODATA = "CODATA 2018"
PLANCK = "Planck Collaboration 2020, A&A 641, A6 (arXiv:1807.06209), Abstract und Tab. 2, TT,TE,EE+lowE+lensing"
_D03 = "03-von-der-galaxie-bis-zum-kosmos.md"

# ------------------------------------------------------------------ Naturkonstanten und Körper
_m("G", "G", 6.67430e-11, None, "m³/(kg·s²)", CODATA, "G")
_m("m_p", "m_p", 1.67262192369e-27, None, "kg", CODATA, "G")
_m("m_e", "m_e", 9.1093837015e-31, None, "kg", CODATA, "G")
_m("m_u", "m_u", 1.66053906660e-27, None, "kg", CODATA, "G")
_m("eps0", "ε₀", 8.8541878128e-12, None, "F/m", CODATA, "G")
_m("alpha", "α", 1 / 137.035999084, None, "1", CODATA + " (1/α = 137,035999084)", "G")
_m("r_e", "r_e", 2.8179403262e-15, None, "m", CODATA, "G")
_m("a_B", "a_B", 5.29177210903e-11, None, "m", CODATA, "G")
_m("R_p", "R_p", 0.8414e-15, None, "m", CODATA + ", Ladungsradius des Protons", "G")
_m("GM_sun", "GM☉", 1.32712440018e20, None, "m³/s²", "IAU", "G")
_m("R_sun", "R☉", 6.957e8, None, "m", "IAU 2015 B3 nominal", "G")
_m("L_sun", "P_L☉", 3.828e26, None, "W", "IAU 2015 B3 nominal", "G")
_m("TEFF_sun", "T_eff☉", 5772.0, None, "K", "IAU 2015 B3 nominal", "G")
_m("GM_earth", "GM⊕", 3.986004418e14, None, "m³/s²", "IAU/IERS", "G")
_m("R_earth", "R⊕", 6.371e6, None, "m", "mittlerer Erdradius", "G")
_m("BETA_PAGE", "β_Page", 1.4847, None, "1", "Page 1976: Parameter der Page-Kurve (Photonen, Neutrinos, Gravitonen)", "G")

# ------------------------------------------------------------------ Kosmos
_m("H0_KMS", "H₀", 67.4, None, "km/s/Mpc", PLANCK, "literatur")
_m("OMEGA_B_H2", "ω_b", 0.02237, 0.00015, "1", PLANCK + ", 1σ", "literatur")
_m("OMEGA_C_H2", "ω_c", 0.1200, 0.0012, "1", PLANCK + ", 1σ", "literatur")
_m("OMEGA_R", "Ω_r", 9.14e-5, None, "1", "Photonen + Neutrinos (Standardwert)", "G")
_m("T_CMB", "T_CMB", 2.72548, 0.00057, "K", "Fixsen 2009, ApJ 707, 916 (arXiv:0911.1955): FIRAS mit Literatur kombiniert", "literatur")
_m("N_EFF_STD", "N_eff", 3.046, None, "1", "Standardwert drei Neutrinoarten", "G")
_m("N_EFF_MESS", "N_eff", 2.99, 0.17, "1", "Planck 2018 + BAO (Abstract), 1σ", "G",
   modell_texte=((_D03, "gemessen N_eff = 2,99 ± 0,17"),))
_m("OMEGA_K_MESS", "Ω_K", 0.0007, 0.0037, "1", "Planck 2018 + BAO, Fehler 95 %", "G",
   modell_texte=((_D03, "gemessen Ω_K = 0,0007 ± 0,0037"),))
_m("KREISPAAR_GRENZE", "L_min", 0.985, None, "d_LSS",
   "COMPACT 2023 (JCAP 01/030): kürzeste Schleife durch unseren Ort > 0,985·d_LSS (Planck 2015, kubischer Torus: 0,97)", "G",
   modell_texte=((_D03, "länger als 0,985-mal"), (_D03, "L ≥ 0,985·d_LSS")))

# ------------------------------------------------------------------ Galaxien
_m("A0_LIT", "a₀", 1.2e-10, (0.02e-10, 0.24e-10), "m/s²",
   "McGaugh, Lelli, Schombert 2016: g† = 1,20 ± 0,02 (stat) ± 0,24 (sys), nur zum Vergleich", "G",
   modell_texte=((_D03, "Literaturwert 1,2·10⁻¹⁰"),))
_m("A0_SPARC", "a₀", (1.04e-10, 1.19e-10), None, "m/s²",
   "SPARC mit realistischen Fehlern von Abstand, Neigung und Masse-zu-Licht, Υ_d = 0,5 (eigene Rechnung, Notiz 74)",
   "daten", datei="sparc", modell_texte=((_D03, "getrennt geben 1,04–1,19"), (_D03, "a₀ = 1,04–1,19·10⁻¹⁰")))
_m("CRATER_II_SIGMA", "σ_v", 2.7, 0.3, "km/s", "Caldwell et al. 2017, Crater II, 1σ", "literatur",
   modell_texte=((_D03, "2,7 ± 0,3"),))

# ------------------------------------------------------------------ Sonnenort in der Milchstraße (Gezeit; Bausteine 02, 03, 06, 07)
_m("V0_SONNE", "v", 233.0, None, "km/s", "Kreisgeschwindigkeit am Sonnenort", "G")
_m("R0_SONNE", "r", 8.178, None, "kpc", "GRAVITY Collaboration 2019, A&A 625, L10 (arXiv:1904.05721): R₀ = 8178 ± 13 (stat) ± 22 (sys) pc", "literatur")
_m("OORT_A", "A_O", 15.3, None, "km/s/kpc", "Bovy 2017, MNRAS 468, L63 (arXiv:1610.07610): A = 15,3 ± 0,4 ± 0,5 km/s/kpc", "literatur")
_m("OORT_B", "B_O", -11.9, None, "km/s/kpc", "Bovy 2017, MNRAS 468, L63 (arXiv:1610.07610): B = −11,9 ± 0,4 km/s/kpc", "literatur")
_m("RHO_B_SONNE", "ρ", 0.084, None, "M☉/pc³", "Dichte der Baryonen am Sonnenort (McKee et al. 2015)", "G")

# ------------------------------------------------------------------ Stoffwerte Kupfer (Labor, 06)
_m("CU_RHO_EL", "ρ_el", 1.68e-8, None, "Ω·m", "Stoffwert Kupfer", "G")
_m("CU_AL_R", "α_R", 3.93e-3, None, "K⁻¹", "Stoffwert Kupfer", "G")
_m("CU_E", "E_Y", 128e9, None, "Pa", "Stoffwert Kupfer", "G")
_m("CU_AL_E", "α_E", -3.8e-4, None, "K⁻¹", "Stoffwert Kupfer", "G")
_m("CU_AL_TH", "α_th", 16.5e-6, None, "K⁻¹", "Stoffwert Kupfer", "G")
_m("CU_RHO", "ρ_Cu", 8960.0, None, "kg/m³", "Stoffwert Kupfer", "G")
_m("CU_CP", "c_p", 385.0, None, "J/(kg·K)", "Stoffwert Kupfer", "G")
_m("CU_K_TH", "k_th", 401.0, None, "W/(m·K)", "Stoffwert Kupfer", "G")

# ------------------------------------------------------------------ Handwerte der Gegenproben (Körper und Quellen)
_m("MERKUR_A", "a_b", 5.7909e10, None, "m", "Merkur, große Halbachse", "G")
_m("MERKUR_E", "e_B", 0.2056, None, "1", "Merkur, Exzentrizität", "G")
_m("MERKUR_P", "P_b", 87.969, None, "d", "Merkur, siderische Umlaufzeit", "G")
_m("ERDE_J", "J", 8.034e37 * 7.2921e-5, None, "J·s", "Erde, J = I·ω mit I = 8,034·10³⁷ kg·m², ω = 7,2921·10⁻⁵ s⁻¹", "G")
_m("M87_M", "M", 6.5e9, None, "M☉", "EHT Collaboration 2019, ApJL 875, L6 (arXiv:1906.11243): (6,5 ± 0,7)·10⁹ M☉", "literatur", modell_texte=(("08-licht-und-trichter.md", "6,5·10⁹ M☉"),))
_m("M87_D", "D_A", 16.8, None, "Mpc", "EHT Collaboration 2019, ApJL 875, L6 (arXiv:1906.11243): 16,8 (+0,8/−0,7) Mpc", "literatur", modell_texte=(("08-licht-und-trichter.md", "16,8 Mpc"),))
_m("GW150914_MF", "M", 62.0, None, "M☉", "LIGO, Virgo 2016, PRL 116, 061102: 62 M☉", "literatur")
_m("GW150914_F", "F_j", 0.5227, None, "1", "F(j = 0,67), Fit Berti, Cardoso, Will 2006", "G")
_m("PSR1913_M1", "m₁", 1.4398, None, "M☉", "PSR B1913+16, Pulsar", "G")
_m("PSR1913_M2", "m₂", 1.3886, None, "M☉", "PSR B1913+16, Begleiter", "G")
_m("PSR1913_PB", "P_b", 27906.98, None, "s", "PSR B1913+16, Umlaufzeit", "G")
_m("PSR1913_E", "e_B", 0.6171338, None, "1", "PSR B1913+16, Exzentrizität", "G", modell_texte=(("08-licht-und-trichter.md", "e = 0,617"),))
from . import messwerte_02  # noqa: E402,F401
from . import messwerte_03a  # noqa: E402,F401  (Messwerte 03, Teil a)
from . import messwerte_03b  # noqa: E402,F401  (Baustein 03, zweiter Teil)

from . import messwerte_05_06  # noqa: E402,F401
from . import messwerte_07  # noqa: E402,F401
from . import messwerte_08_09  # noqa: E402,F401
from . import messwerte_teilchen_licht  # noqa: E402,F401  (Gruppe Teilchen, Sterne und Licht)
from . import messwerte_galaxie  # noqa: E402,F401  (Gruppe Galaxie, Auftrag F)
from . import messwerte_kosmos  # noqa: E402,F401  (Gruppe Kosmos)
from . import messwerte_gesetze  # noqa: E402,F401  (Messwerte der Gesetze: Lelli 2016, 2019)
