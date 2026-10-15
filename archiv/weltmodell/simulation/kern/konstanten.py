# -*- coding: utf-8 -*-
"""
Konstanten des FTM-Kerns: SI-Definitionen und Einheiten als Literale; Messwerte aus messwerte.py, Setzungen aus
setzungen.py (einziger Ort ihrer Zahlen). Ebene: keine eigene; Zahlen, die alle Ebenen (L, I, Z, G) gemeinsam nutzen.
"""
import math
import numpy as np
from .formeln import formel
from .messwerte import MESSWERTE as _M
from .setzungen import SETZUNGEN as _S
from . import eingaben as _eingaben  # noqa: F401  (füllt das Eingabe-Register)

# ------------------------------------------------------------------ SI
G = _M["G"].wert
c = 299792458.0
hbar = 6.62607015e-34 / (2 * math.pi)   # h exakt (SI 2019)
kB = 1.380649e-23
e_el = 1.602176634e-19
N_A = 6.02214076e23        # exakt (SI 2019)
m_p = _M["m_p"].wert
m_e = _M["m_e"].wert
m_u = _M["m_u"].wert
eps0 = _M["eps0"].wert
alpha = _M["alpha"].wert
r_e = _M["r_e"].wert
a_B = _M["a_B"].wert
R_p = _M["R_p"].wert

GM_sun = _M["GM_sun"].wert
M_sun = GM_sun / G
R_sun = _M["R_sun"].wert
L_sun = _M["L_sun"].wert
TEFF_sun = _M["TEFF_sun"].wert
GM_earth = _M["GM_earth"].wert
M_earth = GM_earth / G
R_earth = _M["R_earth"].wert
g_earth = _S["G_ERDE"].wert

# ------------------------------------------------------------------ Einheiten
AU = 1.495978707e11
pc = 3.0856775814913673e16
kpc, Mpc, Gpc = 1e3 * pc, 1e6 * pc, 1e9 * pc
ly = 9.4607304725808e15
yr = 3.15576e7             # julianisches Jahr
Gyr = 1e9 * yr

l_P = math.sqrt(hbar * G / c**3)
m_Pl = math.sqrt(hbar * c / G)

# ------------------------------------------------------------------ Setzungen und Messwerte
A0 = _S["A0"].wert
A0_LIT = _M["A0_LIT"].wert
A0_SPARC = _M["A0_SPARC"].wert
H0_KMS = _M["H0_KMS"].wert
OMEGA_B_H2 = _M["OMEGA_B_H2"].wert
OMEGA_C_H2 = _M["OMEGA_C_H2"].wert
OMEGA_BC_H2_FEHLER = (_M["OMEGA_B_H2"].fehler, _M["OMEGA_C_H2"].fehler)
T_CMB = _M["T_CMB"].wert
N_EFF_STD = _M["N_EFF_STD"].wert
N_EFF_MESS = (_M["N_EFF_MESS"].wert, _M["N_EFF_MESS"].fehler)
OMEGA_K_MESS = (_M["OMEGA_K_MESS"].wert, _M["OMEGA_K_MESS"].fehler)
KREISPAAR_GRENZE = _M["KREISPAAR_GRENZE"].wert


@formel('H0_si', 'H = H₀·10³/Mpc',
        ('H',), ('H₀',), 'G', 'hilfe', 'Baustein 03', None)
def H0_si(H0=H0_KMS):
    return H0 * 1e3 / Mpc


h_ = H0_KMS / 100
Q_KOSM = OMEGA_C_H2 / OMEGA_B_H2                 # 5,364
q_KOSM = Q_KOSM / (1 + Q_KOSM)                   # 0,843
OMEGA_B = OMEGA_B_H2 / h_**2
OMEGA_M = OMEGA_B * (1 + Q_KOSM)                 # 0,3134
OMEGA_R = _M["OMEGA_R"].wert
K_S = math.pi**2 / 45                             # Vakuumteil bei Θ = 1
