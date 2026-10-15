# -*- coding: utf-8 -*-
"""
Rückkehr (Ebene I), Rechnungen, noch nicht im Modell: Maxwell-Form des Gewichts w (frühere Form, ersetzt
durch w = 1 − e^{−√(θ̃/θ_c)}).
"""
from kern.konstanten import *
from kern.formeln import formel


@formel('gewicht_w_maxwell', 'w = θ₊²/(θ₊² + θ_c²) (Maxwell-Speicheranteil)',
        ('w',), ('θ', 'a₀'), 'G', 'forschung', 'Baustein 03', None, alternative='gewicht_w')
def gewicht_w_maxwell(theta, a0=A0):
    """Frühere Form (Rechnung, nicht im Modell): Maxwell-Speicheranteil w = θ₊²/(θ₊² + θ_c²)."""
    th = max(theta, 0.0)
    thc = a0 / c
    return th**2 / (th**2 + thc**2)
