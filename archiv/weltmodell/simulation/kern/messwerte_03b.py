# -*- coding: utf-8 -*-
"""
Messwerte zu Baustein 03, zweiter Teil (Expansion, Kaleidoskop, Dunkle Energie, Vorhersagen, Offen). Klasse und Regeln
in messwert.py. Alle aus dem Gedächtnis [G], Literatur im Netz gesperrt.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
_D = "03-von-der-galaxie-bis-zum-kosmos.md"
DESI = "DESI 2024/2025 mit Hintergrundstrahlung und Supernovae (Abstract), w₀w_a-Anpassung"
DRUCK = "Literaturwerte wie in gedaechtnis_umgebung (Plummer-Aufbau: M aus Leuchtkraft mal Υ, R_e, Abstand zum Wirt)"

# ------------------------------------------------------------------ Dunkle Energie, gemessene Form
_m("DE_MAX_Z", "z", 0.4, None, "1", DESI + ": Maximum der Dunklen Energie bei z ≈ 0,4", "G",
   modell_texte=((_D, "Maximum der Dunklen Energie bei z ≈ 0,4"),))
_m("DE_SIGNIFIKANZ", "w_Z", (2.8, 4.2), None, "σ", DESI + ": Abstand von w = −1 je nach Supernova-Satz", "G",
   modell_texte=((_D, "2,8–4,2σ gegen eine Konstante"),))

# ------------------------------------------------------------------ Entropie, Warp
_m("ENTROPIEPRODUKTION", "Ṡ_i", 7.1e86, None, "k_B/s",
   "heutige Entropieproduktion je Hubble-Horizont (Sterne, Staub, Schwarze Löcher; Schätzung der Literatur)", "G",
   modell_texte=((_D, "Ṡ_i = 7,1·10⁸⁶ k/s"),))
_m("WARP_PRAEZESSION", "ω_p", 10.86, None, "km/s/kpc", "Poggio et al. 2020: Warp der Milchstraße, rechtläufig", "G",
   modell_texte=((_D, "+10,9 km/s/kpc"),))

# ------------------------------------------------------------------ Druckgestützte Systeme (σ gemessen; Aufbau)
_m("PAL14_SIGMA", "σ_v", 0.38, 0.12, "km/s", "Jordi et al. 2009, A&A 504, 109 (arXiv:0903.4448): 0,38 ± 0,12 km/s ohne Ausreißer", "literatur", modell_texte=((_D, "0,38 ± 0,12"),))
_m("NGC2419_SIGMA", "σ_v", 4.1, 0.5, "km/s", "NGC 2419, Sichtlinienstreuung, 1σ", "G", modell_texte=((_D, "4,1 ± 0,5"),))
_m("FORNAX_SIGMA", "σ_v", 11.7, 0.9, "km/s", "Walker et al. 2009, AJ 137, 3100: 11,7 ± 0,9 km/s", "literatur", modell_texte=((_D, "11,7 ± 0,9"),))
_m("DF2_SIGMA", "σ_v", 8.5, 2.3, "km/s", "NGC 1052-DF2, Sichtlinienstreuung, 1σ", "G", modell_texte=((_D, "; 8,5;"),))
_m("DF4_SIGMA", "σ_v", 4.2, 2.2, "km/s", "NGC 1052-DF4, Sichtlinienstreuung, 1σ", "G", modell_texte=((_D, "; 4,2)"),))
_m("DRACO_SIGMA", "σ_v", 9.1, 1.2, "km/s", "Draco, Sichtlinienstreuung, 1σ", "G", modell_texte=((_D, "9,1 ± 1,2"),))
# Aufbau je System: baryonische Masse [M☉], Halblichtradius [kpc], heutiger Abstand zum Wirt [kpc]
for _k, _M, _Re, _d in (("PAL14", 2e4, 0.0247, 71.0), ("NGC2419", 9e5, 0.019, 90.0), ("CRATER2", 3.2e5, 1.066, 117.5),
                        ("DRACO", 5.8e5, 0.221, 76.0), ("FORNAX", 4.0e7, 0.71, 147.0), ("DF2", 2.0e8, 2.2, 80.0),
                        ("DF4", 1.5e8, 1.6, 165.0)):
    _m(f"{_k}_M", "M", _M, None, "M☉", DRUCK, "G")
    _m(f"{_k}_RE", "R_e", _Re, None, "kpc", DRUCK, "G")
    _m(f"{_k}_ABSTAND", "r", _d, None, "kpc", DRUCK, "G")

# ------------------------------------------------------------------ Verlauf der zentralen Flächendichte (Lelli)
_m("LELLI_SIGMA_DYN", "Σ_dyn", (1468.0, 1571.0, 0.38, 0.62), None, "M☉/pc²",
   "Lelli, McGaugh, Schombert, Pawlowski 2016 (ApJL 827, L19): Σ_dyn = Σ_M·(1 + Σ₀/Σ_C)^{α}·(Σ₀/Σ_C)^{β}, "
   "Werte (Σ_M, Σ_C, α, β)", "G")
_m("BAO_SKALA", "r", 147.0, None, "Mpc", "Schallhorizont der BAO (Planck 2018), Rechenpunkt linearer Störungen bei z = 1090",
   "G", modell_texte=((_D, "(147 Mpc:"),))
