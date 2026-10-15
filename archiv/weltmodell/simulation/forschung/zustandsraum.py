# -*- coding: utf-8 -*-
"""
Zustandsraum (Ebene Z), Rechnungen, noch nicht im Modell: offener Streifen eines eingebetteten Torus als
Kehrwert ν_T des offenen Anteils (q_form).
"""
from kern.konstanten import *
from kern.formeln import formel


@formel('torus_streifen_nu', 'ν_T = 2πR_T/(b_S·R_T ± 2r_T·sin(b_S/2))',
        ('ν_T',), ('R_T/r_T', 'b_S'), 'Z', 'forschung', 'Baustein 05', None)
def torus_streifen_nu(R_zu_r, breite=1.0, lage="aussen"):
    """Eingebetteter Torus (R, r), Flächenelement r(R + r cos θ)dθdφ. Offener Teil = Streifen der Breite
    'breite' (Radiant) im Rohrwinkel θ um θ = 0 (außen) bzw. θ = π (innen): Anteil (b·R ± 2r·sin(b/2))/(2πR);
    ν = 1/Anteil = 2πR/(b·R ± 2r·sin(b/2))."""
    sgn = 1.0 if lage == "aussen" else -1.0
    return 2 * math.pi * R_zu_r / (breite * R_zu_r + sgn * 2 * math.sin(breite / 2))
