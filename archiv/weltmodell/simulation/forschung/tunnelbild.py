# -*- coding: utf-8 -*-
"""
Rückkehr als Tunneln (Schrödinger-Lesart D·G″ = Γ·G bei E = 0), Rechnung, nicht im Modell: ψ(r*) mit m = ħ/(2D),
D = r*²/τ_D und V = ħ·a_N/c gleich q = e^{−s} (Milchstraße, a_N bei 2r*). Das Bild ist eine Umformung derselben Regel
ohne neue Invariante (Thema rueckkehr_gleichung, gesperrt); der Modelltext führt es nicht. Ehemals Gesetz
M_rueckkehr_tunneln; jetzt Forschung, gerechnet im Test tests/test_tunnelbild.py.
"""
import numpy as np
from scipy.linalg import solve_banded

import kern as K
from kern.eingabe import EINGABEN

REGISTERMODUL = True


def tunnelbild_gegenprobe(a0=K.A0):
    """(ψ(r*) aus der Finite-Differenzen-Lösung, q = e^{−s} der Regel) für die Milchstraße; beide gleich auf 10⁻⁴."""
    Mt = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    rt = K.r_stern(Mt, a0)
    D = rt**2 / K.tau_D(a0)
    m = K.hbar / (2 * D)
    aN = K.G * Mt / (2 * rt)**2
    V = K.hbar * aN / K.c
    xx = np.linspace(0, 12 * rt, 60001)
    h = xx[1] - xx[0]
    kk = 2 * m * V / K.hbar**2
    ab = np.zeros((3, xx.size - 2))
    ab[0, 1:] = 1
    ab[1] = -(2 + kk * h**2)
    ab[2, :-1] = 1
    rr = np.zeros(xx.size - 2)
    rr[0] = -1
    psi = solve_banded((1, 1), ab, rr)[int(rt / h) - 1]
    return psi, float(K.q_(aN, a0))
