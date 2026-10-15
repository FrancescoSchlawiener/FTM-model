# -*- coding: utf-8 -*-
"""Eingaben zu Baustein 02 (Vom Teilchen bis zum Stern); Klasse und Regeln in eingabe.py."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D = "02-vom-teilchen-bis-zum-stern.md"

# ------------------------------------------------------------------ Körper und Abstände
_e("ABSTAND_1FM", "r", 1e-15, "m", "Abstand neben dem Proton, Beispiel des Texts", (_D, "bei 1 fm neben einem Proton"))

# ------------------------------------------------------------------ Polytrope und Zusammensetzung
_e("N_KONVEKTIV", "n_P", 1.5, "1", "Polytrope der konvektiven Hülle", (_D, "konvektiv, n = 1,5"))
_e("N_STRAHLEND", "n_P", 3.0, "1", "Polytrope der strahlenden Hülle", (_D, "strahlend, n = 3"))
_e("MU_E", "μ_e", 2.0, "1", "Masse je Elektron in m_u für He-, C- und O-Kerne (Z/A = ½)")
_e("MU_HE", "μ_m", 4 / 3, "1", "voll ionisiertes Helium, A/(Z + 1)")
_e("MU_SI", "μ_m", 28 / 15, "1", "voll ionisiertes Silizium, A/(Z + 1)")
_e("MU_CO", "μ_m", 1.75, "1", "voll ionisierter Kohlenstoff-Sauerstoff-Kern (C/O 1:1 nach Masse)")
