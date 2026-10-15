# -*- coding: utf-8 -*-
"""
Messwerte der Gruppe Galaxie (Auftrag F): Baryonen am Sonnenort und das baryonische Massenmodell der Milchstraße
(McMillan 2017, ohne Halo). Klasse und Regeln in messwert.py. Alle aus dem Gedächtnis [G], Literatur im Netz gesperrt.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
MCM = "McMillan 2017 (MNRAS 465, 76), baryonische Teile ohne Halo"

_m("SIGMA_BAR_SONNE", "Σ", 47.1, None, "M☉/pc²", "Baryonen der Säule |z| < 1,1 kpc am Sonnenort (McKee et al. 2015)", "G")
# Bulge ρ = ρ₀/(1 + r′/r₀)^α·e^{−(r′/r_cut)²}, r′ = √(R² + (z/q)²): (ρ₀, r₀, r_cut, α, q)
_m("MCMILLAN_BULGE", "ρ", (98.4e9, 0.075, 2.1, 1.8, 0.5), None, "M☉/kpc³, kpc, kpc, 1, 1", MCM, "G")
# Sternscheiben ρ = Σ₀/(2z_d)·e^{−R/R_d − |z|/z_d}: (Σ₀, R_d, z_d)
_m("MCMILLAN_DUENN", "Σ₀", (896e6, 2.50, 0.300), None, "M☉/kpc², kpc, kpc", MCM, "G")
_m("MCMILLAN_DICK", "Σ₀", (183e6, 3.02, 0.900), None, "M☉/kpc², kpc, kpc", MCM, "G")
# Gasscheiben ρ = Σ₀/(4z_d)·e^{−R_m/R − R/R_d}·sech²(z/2z_d): (Σ₀, R_d, R_m, z_d)
_m("MCMILLAN_HI", "Σ₀", (53.1e6, 7.0, 4.0, 0.085), None, "M☉/kpc², kpc, kpc, kpc", MCM, "G")
_m("MCMILLAN_H2", "Σ₀", (2180e6, 1.5, 12.0, 0.045), None, "M☉/kpc², kpc, kpc, kpc", MCM, "G")
