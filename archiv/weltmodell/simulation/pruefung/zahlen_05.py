# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 05 (Zustandsraum). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "05_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Alle Prüfzahlen hier laufen im Parametersatz "standard"; nur 05_s_labor hängt an a₀ (Standardwert).
"""
import math

import kern as K
from kern.eingabe import EINGABEN
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "05"
D = datei(B)


# ------------------------------------------------------------------ Rückkopplungsschleifen
AUSSAGEN.neu("05_vorzeichen_schliessung", (D, "stabilisiert selbsterregt (q = −0,05)"), "Zahl")
PRUEFZAHLEN.neu("05_q_selbsterregt", B, "05_vorzeichen_schliessung", "−0,05",
                lambda p: K.schleife_selbsterregt(*(EINGABEN[k].wert for k in ("SELBST_G", "SELBST_B", "SELBST_DRDT",
                                                                               "SELBST_GTH"))), "1",
                "stabilisiert selbsterregt", haengt_an=("SELBST_G", "SELBST_B", "SELBST_DRDT", "SELBST_GTH"))

# ------------------------------------------------------------------ Schließung

# ------------------------------------------------------------------ Schließung in der Galaxie: Resolvente der Translation
AUSSAGEN.neu("05_resolvente_quasiperiodisch", (D, "invariante Mode 10⁴-fach bei q = 0,9999"), "Zahl")
PRUEFZAHLEN.neu("05_resolvente_invariant", B, "05_resolvente_quasiperiodisch", "10⁴",
                lambda p: K.schleifensumme(EINGABEN["Q_RESOLVENTE"].wert), "1", "invariante Mode",
                haengt_an=("Q_RESOLVENTE",), toleranz=1e-9,
                toleranz_grund="bloße Zehnerpotenz; 1/(1 − 0,9999) ist genau 10⁴")
PRUEFZAHLEN.neu("05_resolvente_flucht", B, "05_resolvente_quasiperiodisch", "567",
                lambda p: K.resolvente_flucht(EINGABEN["Q_RESOLVENTE"].wert, EINGABEN["RESOLVENTE_MODEN"].wert), "1",
                "Fluchtmoden höchstens", haengt_an=("Q_RESOLVENTE", "RESOLVENTE_MODEN"))

# ------------------------------------------------------------------ Gravitation: zwei Kanäle
AUSSAGEN.neu("05_zwei_kanaele", (D, "ist der differentielle Kanal 1,44·10⁹-mal stärker"), "Zahl")
# Verhältnis an der Erdoberfläche; beide Kanäle sind linear im Höhenunterschied Δh, er kürzt sich
PRUEFZAHLEN.neu("05_kanal_verhaeltnis", B, "05_zwei_kanaele", "1,44·10⁹",
                lambda p: K.kanal_verhaeltnis(K.M_earth, K.R_earth, EINGABEN["KANAL_FG"].wert), "1",
                "differentielle Kanal", haengt_an=("KANAL_FG", "R_earth", "GM_earth"))

# ------------------------------------------------------------------ Bindungsoperator: Spielbeispiel des Gedächtnisses
AUSSAGEN.neu("05_gedaechtnis_bindet", (D, "kehrt bei γ = 0,01 an x = 9,75 um"), "Zahl")
PRUEFZAHLEN.neu("05_gedaechtnis_umkehr", B, "05_gedaechtnis_bindet", "9,75",
                lambda p: K.spiel_umkehrpunkt(EINGABEN["SPIEL_V"].wert, EINGABEN["SPIEL_GAMMA"].wert), "1",
                "kehrt bei γ = 0,01", haengt_an=("SPIEL_V", "SPIEL_GAMMA"))

# ------------------------------------------------------------------ Modulformen
# ν = 1/(1 − |q|), |q| = e^{−2π·Im τ}; τ = i/2 und 2i stehen im Satz selbst (kein gewählter Rechenpunkt)
AUSSAGEN.neu("05_nu_bei_tau", (D, "τ = i/2 und 2i geben 1,045 und 1,0000035"), "Zahl")
PRUEFZAHLEN.neu("05_nu_tau_halb", B, "05_nu_bei_tau", "1,045", lambda p: 1 + K.n_B(2 * math.pi * 0.5), "1",
                "τ = i/2 und 2i")
PRUEFZAHLEN.neu("05_nu_tau_zwei", B, "05_nu_bei_tau", "1,0000035", lambda p: K.nu_aus_modul(2j), "1",
                "τ = i/2 und 2i")

AUSSAGEN.neu("05_diskriminante", (D, "q − 24q² + 252q³ − …"), "Zahl")
PRUEFZAHLEN.neu("05_delta_koeffizient_3", B, "05_diskriminante", "252", lambda p: K.delta_reihe(3)[2], "1",
                "q − 24q² + 252q³")

# Normierungen (2π)^k/(Γ(k)ζ(k)) = (2π)^k/ℐ_{k−1}
AUSSAGEN.neu("05_herkunft_zahlen", (D, "24, 240, 504 = (2π)^k/(Γ(k)ζ(k)) für k = 2, 4, 6"), "Zahl")
for _k, _z in ((2, "24"), (4, "240"), (6, "504")):
    PRUEFZAHLEN.neu(f"05_normierung_E{_k}", B, "05_herkunft_zahlen", _z,
                    lambda p, k=_k: (2 * math.pi)**k / K.bose_leiter_integral(k - 1), "1", "(2π)^k/(Γ(k)ζ(k))")
PRUEFZAHLEN.neu("05_zahl_1728", B, "05_herkunft_zahlen", "1728",
                lambda p: sum(f * (2 * math.pi)**k / K.bose_leiter_integral(k - 1) for f, k in ((3, 4), (2, 6))), "1",
                "1728 = 3·240 + 2·504")

AUSSAGEN.neu("05_vakuumteil_modul", (D, "Der Vakuumteil 4π⁴/15 = 4(2π)⁴/240"), "Zahl")
PRUEFZAHLEN.neu("05_vakuumteil_nenner", B, "05_vakuumteil_modul", "15", lambda p: math.pi**4 / K.bose_leiter_integral(3),
                "1", "Der Vakuumteil 4π⁴/15")
