# -*- coding: utf-8 -*-
"""
Messwerte der Gruppe Kosmos (Klasse und Regeln in messwert.py). Alle aus dem Gedächtnis [G], Literatur im Netz gesperrt.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
_NFW = "Linsenreferenz eines massereichen Haufens aus dem Gedächtnis (NFW-Profil, Notiz 28)"

# ------------------------------------------------------------------ Haufen: gelesene Linse als Referenz
_m("NFW_M200", "M", 1e15, None, "M☉", _NFW + ", M₂₀₀", "G")
_m("NFW_C", "c_NFW", 4.0, None, "1", _NFW + ", Konzentration", "G")
_m("NFW_Z", "z", 0.3, None, "1", _NFW + ", Rotverschiebung der Linse", "G")

# ------------------------------------------------------------------ Anfangsspektrum für σ(M) im Lauf des eigenen Atmens
_PL = "Planck Collaboration 2020, A&A 641, A6, Tab. 2 (TT,TE,EE+lowE+lensing; A_s = e^3,044·10⁻¹⁰), Eingaben des Spektrums wie in de_fluchtmoden"
_m("A_S", "A_s", 2.1e-9, None, "1", _PL, "literatur")
_m("N_S", "n_s", 0.9649, None, "1", _PL, "literatur")
_m("TAU_RE", "τ_re", 0.0544, None, "1", _PL, "literatur")
_m("M_NU", "Σm_ν", 0.06, None, "eV", "kleinste Summe der Neutrinomassen, Eingabe des Spektrums wie in de_fluchtmoden", "G")
