# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 06 (Torus-Dynamik). So füllen: je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "06_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Alle Prüfzahlen hier laufen im Parametersatz "standard"; nur 06_s_labor hängt an a₀ (Standardwert).
"""
import math

import kern as K
from kern.eingabe import EINGABEN
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "06"
D = datei(B)

# ------------------------------------------------------------------ Kopplungsfenster: Schmelzen und Sieden (PubChem)
def _siede_zu_schmelz(p):
    """T_Siede/T_Schmelz aller Elemente mit beiden Werten, ohne die Sonderfälle der Auswahl."""
    aus = set(EINGABEN["ELEMENTE_SONDERFAELLE"].wert)
    werte = []
    for z in range(1, 119):
        if z in aus:
            continue
        e = K.pubchem(z)
        if math.isfinite(e["T_schm"]) and math.isfinite(e["T_sied"]):
            werte.append(e["T_sied"] / e["T_schm"])
    return werte


AUSSAGEN.neu("06_fenster_elemente", (D, "streut T_Siede/T_Schmelz von 1,02 bis 8,2"), "Zahl")
PRUEFZAHLEN.neu("06_siede_schmelz_min", B, "06_fenster_elemente", "1,02", lambda p: min(_siede_zu_schmelz(p)), "1",
                "T_Siede/T_Schmelz", haengt_an=("ELEMENTE_SONDERFAELLE",))
PRUEFZAHLEN.neu("06_siede_schmelz_max", B, "06_fenster_elemente", "8,2", lambda p: max(_siede_zu_schmelz(p)), "1",
                "T_Siede/T_Schmelz", haengt_an=("ELEMENTE_SONDERFAELLE",))

# ------------------------------------------------------------------ Stabilität: T² × S durchgerechnet
AUSSAGEN.neu("06_torus_durchgerechnet", (D, "Durchgerechnet (T² × S, Dissipation D = 0,5)"), "Zahl")
PRUEFZAHLEN.neu("06_rate_quer", B, "06_torus_durchgerechnet", "−0,5", lambda p: K.rate_quer(EINGABEN["TORUS_D"].wert),
                "1", "λ_⊥ = −0,5", haengt_an=("TORUS_D",))
PRUEFZAHLEN.neu("06_rate_rastpunkt", B, "06_torus_durchgerechnet", "−0,098",
                lambda p: K.rate_rastpunkt(EINGABEN["TORUS_K"].wert, EINGABEN["TORUS_PHI"].wert,
                                           EINGABEN["TORUS_DOMEGA"].wert), "1", "ein Exponent wird",
                haengt_an=("TORUS_K", "TORUS_PHI", "TORUS_DOMEGA"))

# Dreieck Energie–Struktur–Zustand: (Σt_i)(Σ1/t_i) − 1 ist nach der Mittelungleichung am kleinsten bei gleichen
# Zeitskalen t_i, dort gleich 3·3 − 1; geprüft wird dieser Kleinstwert.
GLEICHE_ZEITSKALEN = (1.0, 1.0, 1.0)
AUSSAGEN.neu("06_aufschaukeln_dreieck", (D, "oszillatorisch erst bei |q₃| > (Σt_i)(Σ1/t_i) − 1 ≥ 8"), "Zahl")
PRUEFZAHLEN.neu("06_dreieck_kleinster_wert", B, "06_aufschaukeln_dreieck", "8",
                lambda p: sum(GLEICHE_ZEITSKALEN) * sum(1 / t for t in GLEICHE_ZEITSKALEN) - 1, "1",
                "(Σt_i)(Σ1/t_i) − 1 ≥ 8")

# ------------------------------------------------------------------ Kippstelle im Kern, Geisterbahn
AUSSAGEN.neu("06_spaltbarkeit", (D, "die Kippstelle liegt bei x = 1 (größter gemessener Wert 0,93)"), "Zahl")
def _spaltbarkeit_max(p):
    """Größte Spaltbarkeit über alle gemessenen Kerne aus AME2020, Koeffizienten aus der Anpassung (A ≥ A_min)."""
    _, a_S, a_C, _, _ = K.troepfchen_anpassung(EINGABEN["TROEPFCHEN_A_MIN"].wert)
    x = 0.0
    for Z in range(1, 121):
        for A in range(Z, 301):
            try:
                k = K.ame2020(Z, A)
            except KeyError:
                continue
            if k["gemessen"]:
                x = max(x, K.spaltbarkeit(Z, A, a_C, a_S))
    return x


PRUEFZAHLEN.neu("06_spaltbarkeit_max", B, "06_spaltbarkeit", "0,93", _spaltbarkeit_max, "1", "größter gemessener Wert",
                haengt_an=("TROEPFCHEN_A_MIN",))

AUSSAGEN.neu("06_geisterbahn", (D, "um 0,0125–0,015 im gerechneten Beispiel"), "Zahl")
for _k, _z in (("min", "0,0125"), ("max", "0,015")):
    PRUEFZAHLEN.neu(f"06_geisterbahn_{_k}", B, "06_geisterbahn", _z,
                    None, "1", "Rastfenster", offen="fehlt_formel",
                beschreibung="Verschiebung des Rastfensters beim Durchfahren nahe dem verschwundenen "
                                 "Sattel-Knoten (Rampe durch die Adler-Gleichung)")

# ------------------------------------------------------------------ Laborprüfungen
AUSSAGEN.neu("06_labor_s", (D, "Im Labor ist s ≈ 3·10⁵"), "Zahl")
PRUEFZAHLEN.neu("06_s_labor", B, "06_labor_s", "3·10⁵", lambda p: K.s_(K.g_earth, p.a0), "1", "Im Labor ist s",
                haengt_an=("G_ERDE",))
PRUEFZAHLEN.neu("06_gezeit_labor", B, "06_labor_s", "10⁻²⁶",
                None, "m/s²", "Gezeit über 1 km",
                toleranz=0.7, toleranz_grund="bloße Zehnerpotenz: Größenordnung (0,3 … 1,7·10⁻²⁶)", offen="fehlt_formel",
                beschreibung="Änderung des galaktischen Zusatzes über 1 km, 𝒯·ℓ·(ν + dν/d ln g) am Sonnenort; "
                             "gezeit_sonnenort steht im Kern, ν + dν/d ln g nicht, Messwerte des Sonnenorts (Oort A, B, "
                             "v₀, R₀, ρ_b) kommen aus Teil 03")

AUSSAGEN.neu("06_rc_kette", (D, "R = 1 kΩ, n_R = 20, RG = 0,01: ν = 1,1567"), "Zahl")
PRUEFZAHLEN.neu("06_rc_nu", B, "06_rc_kette", "1,1567",
                lambda p: K.rc_kette_nu(EINGABEN["RC_N_R"].wert, EINGABEN["RC_RG"].wert), "1", "ν = 1,1567",
                haengt_an=("RC_N_R", "RC_RG"))


def _rc_steigung(p):
    """d ln(ν − ½)/d ln RG zwischen den beiden tiefsten Punkten (ν − ½ ≈ 1/s im Tiefbereich)."""
    n_R, (a, b) = EINGABEN["RC_N_R"].wert, EINGABEN["RC_RG_TIEF"].wert
    y = [math.log(K.rc_kette_nu(n_R, x) - 0.5) for x in (a, b)]
    return (y[1] - y[0]) / math.log(b / a)


PRUEFZAHLEN.neu("06_rc_wurzel", B, "06_rc_kette", "−0,499", _rc_steigung, "1", "(Steigung −0,499)",
                haengt_an=("RC_N_R", "RC_RG_TIEF"))
PRUEFZAHLEN.neu("06_rc_offenes_ende", B, "06_rc_kette", "9,0", lambda p: 100 * K.abweichung_schluss(math.pi / 2), "%",
                "1/tanh(π/2) − 1 = 9,0 %")

# ------------------------------------------------------------------ Tafel der Aufbauten
AUSSAGEN.neu("06_vo2_mikrobruecke", (D, "Hopf bei 3,55 V (3,26 MHz); bei 3,47 V"), "Zahl")
PRUEFZAHLEN.neu("06_vo2_hopf_spannung", B, "06_vo2_mikrobruecke", "3,55",
                lambda p: K.vo2_hopf(EINGABEN["VO2_AUFBAU"].wert)[0], "V", "VO₂ auf einer Mikrobrücke",
                haengt_an=("VO2_AUFBAU",))
PRUEFZAHLEN.neu("06_vo2_hopf_frequenz", B, "06_vo2_mikrobruecke", "3,26",
                lambda p: K.vo2_hopf(EINGABEN["VO2_AUFBAU"].wert)[1] / 1e6, "MHz", "VO₂ auf einer Mikrobrücke",
                haengt_an=("VO2_AUFBAU",))
PRUEFZAHLEN.neu("06_vo2_torus_spannung", B, "06_vo2_mikrobruecke", "3,47",
                None, "V", "VO₂ auf einer Mikrobrücke", offen="fehlt_formel",
                beschreibung="VO₂-Mikrobrücke im Magnetfeld (Leitfähigkeitssprung, Hall, Wärme; Hopf und "
                             "Neimark-Sacker)")

AUSSAGEN.neu("06_resonatoren", (D, "erstes Einrasten bei (N − 1)k ≈ 0,57–0,63 W"), "Zahl")
for _k, _z, _f in (("min", "0,57", min), ("max", "0,63", max)):
    PRUEFZAHLEN.neu(f"06_einrasten_{_k}", B, "06_resonatoren", _z,
                    lambda p, f=_f: f(K.kuramoto_einrasten(EINGABEN["KURAMOTO_N"].wert)), "1", "erstes Einrasten",
                    haengt_an=("KURAMOTO_N",))

# Kupferdraht: Raster (Q_m, V) als Eingabe, Hopf-Grenze und Falte des heißen Asts aus dem Kern
def _kupfer(p):
    """(Q_m, Abstand Hopf–Falte, Ringfrequenz) je Rasterpunkt mit Hopf-Grenze."""
    aus = []
    for Qm, V in EINGABEN["KUPFER_RASTER"].wert:
        abst = K.kupfer_hopf_abstand(Qm, V)
        if math.isfinite(abst):
            aus.append((Qm, abst, K.kupfer_ringfrequenz(Qm, V)))
    return aus


AUSSAGEN.neu("06_kupfer_hopf", (D, "endet der heiße Ast 0,0010–0,6 Hz vor seiner Falte"), "Zahl")
# Stoffwerte Kupfer [G] gehen über die Voreinstellung von draht() ein, nicht über einen Formelaufruf
_KH = ("KUPFER_RASTER", "CU_RHO_EL", "CU_AL_R", "CU_E", "CU_AL_E", "CU_AL_TH", "CU_RHO", "CU_CP", "CU_K_TH",
       "DRAHT_A", "DRAHT_ELL", "DRAHT_B", "DRAHT_SIG0")
PRUEFZAHLEN.neu("06_kupfer_zeitverhaeltnis", B, "06_kupfer_hopf", "0,74",
                lambda p: K.kupfer_zeitverhaeltnis(min(x[0] for x in _kupfer(p))), "1", "Wärmezeit/Ringzeit ≤ 0,74",
                haengt_an=_KH)
PRUEFZAHLEN.neu("06_kupfer_guete_min", B, "06_kupfer_hopf", "3000", lambda p: min(x[0] for x in _kupfer(p)), "1",
                "Güte ab 3000", haengt_an=_KH)
PRUEFZAHLEN.neu("06_kupfer_abstand_min", B, "06_kupfer_hopf", "0,0010", lambda p: min(x[1] for x in _kupfer(p)), "Hz",
                "vor seiner Falte", haengt_an=_KH)
PRUEFZAHLEN.neu("06_kupfer_abstand_max", B, "06_kupfer_hopf", "0,6", lambda p: max(x[1] for x in _kupfer(p)), "Hz",
                "vor seiner Falte", haengt_an=_KH)
PRUEFZAHLEN.neu("06_kupfer_ring_min", B, "06_kupfer_hopf", "0,20", lambda p: min(x[2] for x in _kupfer(p)), "Hz",
                "Ringfrequenz 0,20–1,5 Hz", haengt_an=_KH)
PRUEFZAHLEN.neu("06_kupfer_ring_max", B, "06_kupfer_hopf", "1,5", lambda p: max(x[2] for x in _kupfer(p)), "Hz",
                "Ringfrequenz 0,20–1,5 Hz", haengt_an=_KH)
