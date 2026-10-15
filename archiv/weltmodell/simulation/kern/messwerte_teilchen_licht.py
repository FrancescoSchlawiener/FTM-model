# -*- coding: utf-8 -*-
"""
Messwerte der Gruppe Teilchen, Sterne und Licht (02, 04, 08, 09; Klasse und Regeln in messwert.py): Eingänge der
Formeln für starke Kopplung, Syntropie, Wellenquellen und Kohärenz.
"""
from .messwert import MESSWERTE

REGISTERMODUL = True
_m = MESSWERTE.neu
PDG = "PDG 2022"

# ------------------------------------------------------------------ Starke Kopplung (02)
_m("ALPHA_S_MZ", "α_s", 0.1180, 0.0009, "1", PDG + ", α_s(M_Z) (MS-bar), 1σ", "G")
_m("M_Z_GEV", "M_Z", 91.1876, 0.0021, "GeV", PDG + ", Masse des Z-Bosons", "literatur")
_m("M_CHARM", "m_q", 1.27, 0.02, "GeV", PDG + ", m_c(m_c) (MS-bar)", "G")
_m("M_BOTTOM", "m_q", 4.18, 0.03, "GeV", PDG + ", m_b(m_b) (MS-bar)", "G")
_m("M_TOP", "m_q", 172.7, 0.3, "GeV", PDG + ", Masse des Top-Quarks", "G")

# ------------------------------------------------------------------ Syntropie (04)
_m("ALBEDO_ERDE", "A_alb", 0.30, None, "1", "Albedo der Erde (Bond-Albedo, gerundet)", "G")
_m("S_FALTEN_HEUTE", "ΣS_F", 3.1e104, None, "1",
   "Egan, Lineweaver 2010: Entropie der Schwarzen Löcher im beobachtbaren Universum (überwiegend die großen), in k_B", "G")

# ------------------------------------------------------------------ Wellenquellen (08)
_m("GW150914_Q", "Q_RD", 3.1669, None, "1",
   "Q(j = 0,67) = q₁ + q₂(1 − j)^{q₃}, Fit Berti, Cardoso, Will 2006 (q = 0,7000; 1,4187; −0,4990)", "G")
_m("GW150914_E_RAD", "E_rad", 3.0, 0.5, "M☉c²", "GW150914, abgestrahlte Energie im Quellsystem (LIGO 2016)", "literatur")
_m("GW150914_DL", "D_L", 410.0, 170.0, "Mpc", "GW150914, Leuchtkraftabstand (LIGO 2016)", "literatur")

# ------------------------------------------------------------------ Kohärenz (09)
_m("EPS_R_WASSER_37", "ε_r", 74.2, None, "1", "Wasser bei 37 °C (Malmberg, Maryott 1956)", "G")
