# -*- coding: utf-8 -*-
"""Eingaben der Gruppe Teilchen, Sterne und Licht (Klasse und Regeln in eingabe.py): gewählte Rechenpunkte."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D2 = "02-vom-teilchen-bis-zum-stern.md"
_D8 = "08-licht-und-trichter.md"

# ------------------------------------------------------------------ Starke Kopplung (02)
_e("ALPHA_S_AENDERUNG", "α_s", 0.01, "1", "relative Änderung von α_s an der Planck-Skala im Beispiel der Empfindlichkeit",
   (_D2, "1 % mehr oder weniger in α_s"))
_e("ALPHA_S_STARK", "α_s", 0.5, "1", "Kopplung, deren Skala als Stellvertreter von Λ_QCD dient (m_p ∝ Λ bei fester "
   "Dynamik); jede Kopplung ≫ α_s(M_Z) gibt dieselbe Empfindlichkeit")

# ------------------------------------------------------------------ Zustand und Sterne (02)
_e("Z_EISEN", "Z_A", 26, "1", "Eisen im Beispiel der Zustandsschwelle der Erde", (_D2, "Die Erde liegt für Eisen"))
_e("BETA_GLEICHSTAND", "β_E", 0.5, "1", "Gasdruck gleich Strahlungsdruck", (_D2, "Strahlungsdruck und Gasdruck sind"))

# ------------------------------------------------------------------ Syntropie (04)
_e("Z_F_EINS", "z_f", 1.0, "1", "ohne Rotverschiebung: Vorfaktor von z_f² in der Mindestzahl der Rückverdichtung")

# ------------------------------------------------------------------ Kerne (06)
_e("TROEPFCHEN_A_MIN", "A_A", 16, "1", "kleinste Massenzahl der Anpassung des Tröpfchenmodells an AME2020 (leichte Kerne "
   "folgen Schalen, nicht dem Tröpfchen)")

# ------------------------------------------------------------------ Wellen (08)
_e("IOTA_KANTE", "ι", 90.0, "°", "Kantenlage der Quelle", (_D8, "der Wert gilt bei Kantenlage"))

# ------------------------------------------------------------------ Kohärenz (09)
_e("I_PHYSIOLOGISCH", "I_ion", 0.15, "mol/l", "Ionenstärke physiologischer Salzlösung (0,15 mol/l NaCl)")
