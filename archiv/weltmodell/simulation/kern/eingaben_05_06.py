# -*- coding: utf-8 -*-
"""Eingaben der Bausteine 05 (Zustandsraum) und 06 (Torus-Dynamik); Klasse und Regeln in eingabe.py."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D05 = "05-zustandsraum.md"
_D06 = "06-torus-dynamik.md"

# ------------------------------------------------------------------ 05: Resolvente der quasiperiodischen Rückkehr
_e("Q_RESOLVENTE", "q", 0.9999, "1", "Rückkehranteil nahe 1, an dem die Resolvente (I − qT)⁻¹ gerechnet ist "
   "(Rechnung identifikation_galaxie, q = 0,9 … 0,9999)", (_D05, "q = 0,9999"))

# ------------------------------------------------------------------ 05: Spielbeispiel des Gedächtnisses, zwei Kanäle
_e("SPIEL_V", "v_K", 1.6, "1", "Anfangsgeschwindigkeit des Spielbeispiels über der Fluchtgeschwindigkeit √2",
   (_D05, "v = 1,6 > v_esc"))
_e("SPIEL_GAMMA", "γ_K", 0.01, "1", "Kopplung des Maxwell-Spielkerns im Spielbeispiel", (_D05, "γ = 0,01"))
_e("KANAL_DH", "Δh", 1000.0, "m", "Höhenunterschied im Vergleich der Kanäle; das Verhältnis hängt nicht von ihm ab",
   (_D05, "Bei 1000 m Höhenunterschied"))

# ------------------------------------------------------------------ 06: Dissipation des durchgerechneten Torus
_e("TORUS_D", "D_Z", 0.5, "1", "Dissipation des durchgerechneten Torus T² × S", (_D06, "Dissipation D = 0,5"))

# ------------------------------------------------------------------ 06: RC-Kette als Nachbau der Erneuerung
_e("RC_N_R", "n_R", 20, "1", "Abgriff der Rückführung im Beispielaufbau (Rechnung kontakt_atmen)", (_D06, "n_R = 20"))
_e("RC_RG", "RG", 0.01, "1", "R = 1 kΩ, G = 10 µS je Glied: s = 2 im Beispielaufbau", (_D06, "RG = 0,01"))
_e("RC_RG_TIEF", "RG", (1e-5, 10 ** -4.5), "1",
   "die beiden tiefsten Punkte der Rechnung kontakt_atmen (s ≈ 0,06 und 0,11), an denen die Steigung der Wurzel "
   "d ln(ν − ½)/d ln RG gelesen ist")

# ------------------------------------------------------------------ 06: Stoffwerte der Elemente (PubChem)
_e("ELEMENTE_SONDERFAELLE", "Z_A", (2, 6, 33), "1",
   "aus T_Siede/T_Schmelz genommen: He (Schmelzpunkt nur unter Druck), C (sublimiert), As (sublimiert, Schmelzpunkt "
   "über Siedepunkt); Auswahl der Rechnung elemente_analyse")

# ------------------------------------------------------------------ 06: eingespannter Kupferdraht, Raster der Hopf-Suche
# (Q_m, V_D): je Güte zwei Antriebe V = V₁·(1, 2), V₁ = 10 mV·√(851/Q_l) gibt dieselbe Schleifenverstärkung wie
# Q_m = 1000 bei 10 mV (belastete Güte 851); Spannungen auf 1 µV wie in der Rechnung kupfer_selbstpulsieren.
_e("KUPFER_RASTER", "V_D", ((1e3, 10.000e-3), (1e3, 20.000e-3), (2e3, 7.579e-3), (2e3, 15.159e-3),
                            (3e3, 6.577e-3), (3e3, 13.155e-3), (1e4, 4.838e-3), (1e4, 9.676e-3),
                            (3e4, 4.211e-3), (3e4, 8.422e-3), (1e5, 3.968e-3), (1e5, 7.937e-3)), "V",
   "Rechenpunkte (Q_m, V) der Rechnung kupfer_selbstpulsieren, Teil 2: Güte 10³ … 10⁵, Antrieb gleicher "
   "Schleifenverstärkung und doppelt")
