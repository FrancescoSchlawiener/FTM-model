# -*- coding: utf-8 -*-
"""Eingaben der Formeln der Rückkehr, der Faltung und des Zustandsraums (Auftrag F, Gruppe Zustand): gewählte
Rechenpunkte der Beispiele; Klasse und Regeln in eingabe.py."""
from .eingabe import EINGABEN

REGISTERMODUL = True
_e = EINGABEN.neu
_D04 = "04-zeit-verschraenkung-syntropie-gedaechtnis.md"
_D05 = "05-zustandsraum.md"
_D06 = "06-torus-dynamik.md"
_D07 = "07-feldgleichung.md"

# ------------------------------------------------------------------ 05: selbsterregtes Minimalsystem (Rechnung wechselwirkungen, Einheiten L = C = 1)
_e("SELBST_G", "g_E", 0.30, "Ω", "Entdämpfung a − R_rad = 0,35 − 0,05 des Minimalsystems der Rechnung wechselwirkungen")
_e("SELBST_B", "b_I", 0.2, "Ω/A²", "kubische Sättigung des Minimalsystems der Rechnung wechselwirkungen")
_e("SELBST_DRDT", "dR/dT", 0.05, "Ω/K", "α_R·R₀ = 0,5·0,10 des Minimalsystems der Rechnung wechselwirkungen")
_e("SELBST_GTH", "G_th", 1.0, "W/K", "Wärmeleitwert des Minimalsystems der Rechnung wechselwirkungen")

# ------------------------------------------------------------------ 05: Resolvente, zwei Kanäle
_e("RESOLVENTE_MODEN", "K_M", 2000, "1", "Fluchtmoden je Umlaufsinn der Rechnung identifikation_galaxie (Wellenzahlen "
   "|k| ≤ 2000, zusammen 4000 Moden)")
_e("KANAL_FG", "f_g", 1.0, "1", "hängendes Pendel bzw. Draht: die Gravitation trägt die Rückstellung allein "
   "(Rechnung torus_wechselkraefte, Teil 5)")

# ------------------------------------------------------------------ 06: durchgerechneter Torus T² × S
_e("TORUS_K", "k_K", 0.05, "s⁻¹", "Kopplung des durchgerechneten Beispiels T² × S am Rastpunkt 1:1", wahl=True)
_e("TORUS_DOMEGA", "Δω", 0.02, "s⁻¹", "Verstimmung des durchgerechneten Beispiels T² × S am Rastpunkt 1:1", wahl=True)
_e("TORUS_PHI", "φ_K", 0.0, "rad", "rein dissipative Kopplung im durchgerechneten Beispiel T² × S")

_e("VO2_AUFBAU", "𝔄_VO₂", dict(T_0=300.0, T_C=340.0, breite=1.5, R_I=20.0, R_M=0.4, C_th=5e-3, G_th=5e-3, R_L=10.0,
                                C=0.02, kappa=-0.1, f_m=0.5, Q=1000.0, D_tr=0.2, lambda_B=0.063, Gamma=1e-4, x_D=0.3,
                                B=1.0, V_suche=(3.50, 3.60)),
   "µs, K, V, kΩ, mA, mW, nF, nJ/K, µm",
   "Annahmen der Rechnung bistabilitaet (Teil 1): Substrat 300 K, Übergang 340 K (Breite 1,5 K), 20 kΩ/0,4 kΩ, "
   "5 pJ/K, 5 µW/K, R_L = 10 kΩ, 20 pF, Brücke 0,5 MHz, Q = 1000, Dehnung 0,2 µm, κ = −0,1 K/µm, λ_B = 0,063 "
   "µm/(mA·T), Γ = 10⁻⁴ V·s/m je T, x_D = 0,3 µm, B = 1 T; Suchfenster der Hopf-Grenze 3,50–3,60 V")
_e("KURAMOTO_N", "N_K", (4, 8, 16, 32), "1", "Zahlen gekoppelter Zyklen der Rechnung torus_wechselkraefte (Teil 6b)")

# ------------------------------------------------------------------ 04, 07: Windungskorrelation, Gitter, Drehung
_e("CHSH_HARMONISCHE", "N_H", 2, "1", "Grundschwingung und erste Oberschwingung des Dreiecks",
   (_D04, "zwei Harmonische des Dreiecks"))
_e("WATSON_D", "d", 3, "1", "Rückkehr in drei Dimensionen (Pólya, Watson)", (_D07, "in drei ist ν ≤ 1,516"))
_e("DREHUNG_ZELLE", "μ(A)", 0.01, "1", "Bogen der goldenen Drehung, 1 % des Kreises", (_D07, "Zelle 1 %"))

# ------------------------------------------------------------------ 09: eingespannter Kupferdraht, Zerlegung der Schleife
_e("DRAHT_ARBEITSPUNKTE", "Δf_d", ((2e-3, -10.0), (2e-3, -2.0), (2e-3, 0.0), (2e-3, 2.0),
                                   (5e-3, -10.0), (5e-3, -7.0), (5e-3, -4.0), (5e-3, 0.0),
                                   (10e-3, -40.0), (10e-3, -30.0), (10e-3, -25.0), (10e-3, -20.0), (10e-3, -10.0)),
   "V, Hz", "Arbeitspunkte (V_D, Δf_d) der Tabelle der Rechnung torus_wechselkraefte (Teil 4a, Q_m = 1000), je Punkt "
   "die kälteste stationäre Lösung")
_e("DRAHT_Q_MERKLICH", "q", 0.1, "1", "Arbeitspunkte, an denen die Schleife merklich trägt (|q| ≥ 0,1); darunter ist "
   "die Wärme ohne Rückwirkung und das Verhältnis der Kanäle ohne Belang", wahl=True)
