# -*- coding: utf-8 -*-
"""
Prüfzahlen Baustein 08 (Licht und Trichter): je Aussage des Modelltexts einmal AUSSAGEN.neu(kennung, (D, Stichwort),
"Zahl"), dann je gedruckte Zahl PRUEFZAHLEN.neu(...) mit den Feldern aus pruefzahl.py. Kennungen beginnen mit "08_".
Rechnung frisch aus dem Kern (K.…), Messwerte und Eingaben nur über ihre Register, genannt in haengt_an.
Winkel: 1 rad = (180/π)·3600″; Uhrgang je Jahr: julianisches Jahr K.yr.
"""
import math

import kern as K
from kern.eingabe import EINGABEN as E
from kern.messwert import MESSWERTE as M
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

B = "08"
D = datei(B)
_A = AUSSAGEN.neu
_P = PRUEFZAHLEN.neu
BOGENSEK = 180 / math.pi * 3600
M_MW = E["M_MILCHSTRASSE"].wert * K.M_sun
M_CL = E["M_HAUFEN"].wert * K.M_sun

# ------------------------------------------------------------------ Lichtablenkung
_A("08_ablenkung_sonnenrand", (D, "Sonnenrand (R☉ = 6,957·10⁸ m): 0,8756″ + 0,8756″ = 1,7512″"), "Zahl")
_P("08_ablenkung_sonnenrand", B, "08_ablenkung_sonnenrand", "1,7512", lambda p: K.ablenkung(K.M_sun, K.R_sun) * BOGENSEK,
   "″", "Sonnenrand", haengt_an=("GM_sun", "R_sun"))
_P("08_ablenkung_takt_haelfte", B, "08_ablenkung_sonnenrand", "0,8756",
   lambda p: K.ablenkung(K.M_sun, K.R_sun) / 2 * BOGENSEK, "″", "Sonnenrand", haengt_an=("GM_sun", "R_sun"))

# ------------------------------------------------------------------ Taktunterschied über die Zeit
_A("08_takt_jahr", (D, "ISS 0,021 s pro Jahr, Erdbahn im Sonnentrichter 0,31 s"), "Zahl")
_P("08_takt_iss", B, "08_takt_jahr", "0,021",
   lambda p: K.delta(K.M_earth, K.R_earth + E["ISS_HOEHE"].wert * 1e3) * K.yr, "s", "ISS",
   haengt_an=("GM_earth", "R_earth", "ISS_HOEHE"))
_P("08_takt_erdbahn", B, "08_takt_jahr", "0,31", lambda p: K.delta(K.M_sun, K.AU) * K.yr, "s",
   "Erdbahn im Sonnentrichter", haengt_an=("GM_sun",))

# ------------------------------------------------------------------ Schwelle zur Hülle: zwei Klassen
_A("08_schwelle_huelle", (D, "Nach der Lage von r* zur Hülle gibt es zwei Klassen"), "Zahl")
_P("08_rR_proton", B, "08_schwelle_huelle", "38", lambda p: K.schwelle_zu_huelle(K.m_p, K.R_p, p.a0), "1", "Proton",
   haengt_an=("m_p", "R_p"))
_P("08_rR_erde", B, "08_schwelle_huelle", "3·10⁵", lambda p: K.schwelle_zu_huelle(K.M_earth, K.R_earth, p.a0), "1", "Erde",
   haengt_an=("GM_earth", "R_earth"))
_P("08_rR_sonne", B, "08_schwelle_huelle", "1,6·10⁶", lambda p: K.schwelle_zu_huelle(K.M_sun, K.R_sun, p.a0), "1", "Sonne",
   haengt_an=("GM_sun", "R_sun"))
_P("08_rR_milchstrasse", B, "08_schwelle_huelle", "0,6",
   lambda p: K.schwelle_zu_huelle(M_MW, E["R_MILCHSTRASSE_HUELLE"].wert * K.kpc, p.a0), "1", "Galaxien und Haufen",
   haengt_an=("M_MILCHSTRASSE", "R_MILCHSTRASSE_HUELLE", "GM_sun"))
_P("08_rR_haufen", B, "08_schwelle_huelle", "0,5", lambda p: K.schwelle_zu_huelle(M_CL, E["R_HAUFEN_HUELLE"].wert * K.Mpc, p.a0),
   "1", "Galaxien und Haufen", haengt_an=("M_HAUFEN", "R_HAUFEN_HUELLE", "GM_sun"))

# ------------------------------------------------------------------ Leiter der Längen und Schnittmassen
_A("08_leiterzahl", (D, "ε_L = ℓ_P/ℓ₀ = a₀/a_P = 2,0·10⁻⁶²"), "Zahl")
_P("08_leiterzahl", B, "08_leiterzahl", "2,0·10⁻⁶²", lambda p: K.leiterzahl(p.a0), "1", "ε_L = ℓ_P/ℓ₀")
_A("08_schnittmasse_mitte", (D, "M_c = ε_L^{1/3}·m_Pl = 5,9·10⁻²⁹ kg"), "Zahl")
_P("08_schnittmasse_mitte", B, "08_schnittmasse_mitte", "5,9·10⁻²⁹", lambda p: K.masse_r_stern_gleich_compton(p.a0), "kg",
   "M_c = ε_L^{1/3}·m_Pl")
_A("08_schnittmasse_oben", (D, "M_U = m_Pl/ε_L = c³τ_D/G = 1,1·10⁵⁴ kg"), "Zahl")
_P("08_schnittmasse_oben", B, "08_schnittmasse_oben", "1,1·10⁵⁴", lambda p: K.schnittmasse_oben(p.a0), "kg",
   "M_U = m_Pl/ε_L")
_A("08_leiter_proton_vakuum", (D, "r*/λ̄ = (M/M_c)^{3/2} (Proton 151)"), "Zahl")
_P("08_rstern_compton_proton", B, "08_leiter_proton_vakuum", "151",
   lambda p: K.r_stern(K.m_p, p.a0) / K.compton_red(K.m_p), "1", "r*/λ̄ = (M/M_c)^{3/2}", haengt_an=("m_p",))
_P("08_vakuum_planck", B, "08_leiter_proton_vakuum", "1,0·10⁻¹²²",
   lambda p: K.lambda_vakuum(K.schliesslaenge(p.a0)) * K.l_P**2, "1", "Λ_s·ℓ_P²")

# ------------------------------------------------------------------ Tafel r*
_A("08_tafel_rstern", (D, "| System | r* |"), "Zahl")
_P("08_rstern_erde", B, "08_tafel_rstern", "12,7", lambda p: K.r_stern(K.M_earth, p.a0) / K.AU, "AE", "| Erde |",
   haengt_an=("GM_earth",))
_P("08_rstern_sonne_ae", B, "08_tafel_rstern", "7 340", lambda p: K.r_stern(K.M_sun, p.a0) / K.AU, "AE", "| Sonne |",
   haengt_an=("GM_sun",))
_P("08_rstern_sonne_lj", B, "08_tafel_rstern", "0,116", lambda p: K.r_stern(K.M_sun, p.a0) / K.ly, "Lj", "| Sonne |",
   haengt_an=("GM_sun",))
_P("08_rstern_milchstrasse", B, "08_tafel_rstern", "8,72", lambda p: K.r_stern(M_MW, p.a0) / K.kpc, "kpc",
   "| Milchstraße |", haengt_an=("M_MILCHSTRASSE", "GM_sun"))

# ------------------------------------------------------------------ Verstärkung an r* und Licht am kollektiven Trichter
_A("08_verstaerkung_rstern", (D, "a/a_N = 1/(1 − e⁻¹) = 1,582"), "Zahl")
_P("08_verstaerkung_rstern", B, "08_verstaerkung_rstern", "1,582", lambda p: K.nu(p.a0, p.a0), "1", "a/a_N = 1/(1 − e⁻¹)")
_A("08_linse_I1", (D, "I(1) = 2,135"), "Zahl")
_P("08_linse_I1", B, "08_linse_I1", "2,135",
   lambda p: K.ablenkung_linse(K.M_sun, K.r_stern(K.M_sun, p.a0), p.a0) / K.ablenkung(K.M_sun, K.r_stern(K.M_sun, p.a0)),
   "1", "I(1) =")
_A("08_ablenkung_grenzwert", (D, "Milchstraße 0,43″, Galaxienhaufen 23″"), "Zahl")
_P("08_alpha_inf_milchstrasse", B, "08_ablenkung_grenzwert", "0,43", lambda p: K.ablenkung_kollektiv(M_MW, p.a0) * BOGENSEK,
   "″", "Milchstraße", haengt_an=("M_MILCHSTRASSE", "GM_sun"))
_P("08_alpha_inf_haufen", B, "08_ablenkung_grenzwert", "23", lambda p: K.ablenkung_kollektiv(M_CL, p.a0) * BOGENSEK,
   "″", "Galaxienhaufen", haengt_an=("M_HAUFEN", "GM_sun"))

# ------------------------------------------------------------------ Rotation und Verwindung
_A("08_verwindung_erde", (D, "Erde: a_J = 3,27 m, Δs = 1,48 mm"), "Zahl")
_P("08_aJ_erde", B, "08_verwindung_erde", "3,27", lambda p: K.verwindung(M["ERDE_J"].wert, K.M_earth), "m", "Erde: a_J",
   haengt_an=("ERDE_J", "GM_earth"))
_P("08_ds_erde", B, "08_verwindung_erde", "1,48", lambda p: K.ds(K.M_earth) * 1e3, "mm", "Erde: a_J",
   haengt_an=("GM_earth",))
_A("08_umlauf_differenz", (D, "Erde: 137 ns pro Umlauf"), "Zahl")
_P("08_umlauf_differenz_erde", B, "08_umlauf_differenz", "137",
   lambda p: K.umlauf_differenz(K.verwindung(M["ERDE_J"].wert, K.M_earth)) * 1e9, "ns", "pro Umlauf",
   haengt_an=("ERDE_J", "GM_earth"))
_A("08_perihel_merkur", (D, "Merkur: 43″ pro Jahrhundert"), "Zahl")
_P("08_perihel_merkur", B, "08_perihel_merkur", "43",
   lambda p: (K.perihel_pro_umlauf(K.M_sun, M["MERKUR_A"].wert, M["MERKUR_E"].wert) * BOGENSEK
              * 36525 / M["MERKUR_P"].wert), "″/Jahrhundert", "Merkur:",
   haengt_an=("GM_sun", "MERKUR_A", "MERKUR_E", "MERKUR_P"))
_A("08_isco_kerr", (D, "Die innerste stabile Bahn liegt bei j = 0,9"), "Zahl")
_P("08_isco_mitlaeufig", B, "08_isco_kerr", "2,32", lambda p: K.isco_kerr(E["J_ISCO"].wert), "3Δs",
   "Die innerste stabile Bahn liegt bei j = 0,9", haengt_an=("J_ISCO",))
_P("08_isco_gegenlaeufig", B, "08_isco_kerr", "8,72", lambda p: K.isco_kerr(E["J_ISCO"].wert, mitlaeufig=False), "3Δs",
   "(gegenläufig)", haengt_an=("J_ISCO",))
_A("08_knotenzahl_schwelle", (D, "Knotenzahl ≈ 5·10⁻⁷ in der Milchstraße"), "Zahl")
_P("08_knotenzahl_milchstrasse", B, "08_knotenzahl_schwelle", "5·10⁻⁷",
   lambda p: K.knotendrehung(M_MW, K.r_stern(M_MW, p.a0), E["J_MILCHSTRASSE"].wert), "1", "Knotenzahl",
   haengt_an=("M_MILCHSTRASSE", "J_MILCHSTRASSE", "GM_sun"))

# ------------------------------------------------------------------ Gravity Probe B (Bahnmittel des Vektors längs der Erdachse)
_A("08_gravity_probe_b", (D, "Gravity Probe B (polarer Erdorbit, 642 km Höhe)"), "Zahl")
_P("08_gpb_geodaetisch", B, "08_gravity_probe_b", "6604",
   lambda p: K.geodaetisch(K.M_earth, E["GPB_A"].wert * 1e3) * K.yr * BOGENSEK * 1e3, "mas/Jahr", "geodätische Drehung",
   haengt_an=("GM_earth", "GPB_A"))
_P("08_gpb_mitfuehrung", B, "08_gravity_probe_b", "40,8",
   lambda p: (K.mitfuehrung_polarbahn(K.M_earth, K.verwindung(M["ERDE_J"].wert, K.M_earth), E["GPB_A"].wert * 1e3)
              * K.yr * BOGENSEK * 1e3), "mas/Jahr", "Mitführung nach ART", haengt_an=("GM_earth", "ERDE_J", "GPB_A"))

_P("08_gpb_leitstern", B, "08_gravity_probe_b", "39,1",
   lambda p: (K.leitstern_anteil(K.mitfuehrung_polarbahn(K.M_earth, K.verwindung(M["ERDE_J"].wert, K.M_earth), E["GPB_A"].wert * 1e3),
                                 math.radians(E["GPB_LEITSTERN_DEKL"].wert)) * K.yr * BOGENSEK * 1e3), "mas/Jahr",
   "Leitstern IM Peg", haengt_an=("GM_earth", "ERDE_J", "GPB_A", "GPB_LEITSTERN_DEKL"))

# ------------------------------------------------------------------ Gravitationswellen und Starkfeld (Formelblock)
_A("08_hulse_taylor", (D, "Hulse-Taylor: −2,40·10⁻¹²"), "Zahl")
_PSR = lambda e: K.hulse_taylor(M["PSR1913_M1"].wert, M["PSR1913_M2"].wert, M["PSR1913_PB"].wert, e)  # noqa: E731
_PSR_G = ("PSR1913_M1", "PSR1913_M2", "PSR1913_PB", "PSR1913_E", "GM_sun")
_P("08_hulse_taylor", B, "08_hulse_taylor", "−2,40·10⁻¹²", lambda p: _PSR(M["PSR1913_E"].wert), "1", "Hulse-Taylor:",
   haengt_an=_PSR_G)
_A("08_ringdown", (D, "Ringdown GW150914: 249,9 Hz"), "Zahl")
_P("08_ringdown_detektor", B, "08_ringdown", "249,9",
   lambda p: K.ringdown_frequenz(M["GW150914_MF"].wert * K.M_sun, M["GW150914_F"].wert) / (1 + M["GW150914_Z"].wert),
   "Hz", "Ringdown GW150914", haengt_an=("GW150914_MF", "GW150914_F", "GW150914_Z", "GM_sun"))
_A("08_schatten", (D, "Schatten M87*: Ø ≈ 40 μas"), "Zahl")
_P("08_schatten_stosszahl", B, "08_schatten", "15,6", lambda p: K.schatten_winkel(1.0, 1.0) / 2 / K.ds(1.0), "Δs",
   "b = 3√27·Δs")
_P("08_schatten_m87", B, "08_schatten", "40",
   lambda p: K.schatten_winkel(M["M87_M"].wert * K.M_sun, M["M87_D"].wert * K.Mpc) * BOGENSEK * 1e6, "μas",
   "Schatten M87*", haengt_an=("M87_M", "M87_D", "GM_sun"))
_A("08_faltungsmesser", (D, "GW150914 15 km, GW170817 590 m"), "Zahl")
_P("08_faltungsmesser_gw150914", B, "08_faltungsmesser", "15", lambda p: K.ds(M["GW150914_MC_DET"].wert * K.M_sun) / 1e3,
   "km", "GW150914", haengt_an=("GW150914_MC_DET", "GM_sun"))
_P("08_faltungsmesser_gw170817", B, "08_faltungsmesser", "590", lambda p: K.ds(M["GW170817_MC_DET"].wert * K.M_sun), "m",
   "GW170817", haengt_an=("GW170817_MC_DET", "GM_sun"))
_A("08_schwellenfrequenz", (D, "Ω(r*)/π ≲ 10⁻¹³ Hz"), "Zahl")
_P("08_schwellenfrequenz_sonne", B, "08_schwellenfrequenz", "10⁻¹³", lambda p: K.schwellenfrequenz_gw(K.M_sun, p.a0), "Hz",
   "f*_gw", haengt_an=("GM_sun",), toleranz=0.1,
   toleranz_grund="Obergrenze „≲“ als Größenordnung, gerechnet für die leichteste Quelle 1 M☉")

# ------------------------------------------------------------------ Doppelpulsar, Ringdown im Quellsystem, M87* (Text)
_A("08_doppelpulsar", (D, "Beim Doppelpulsar Hulse-Taylor (Gesamtmasse 2,83 M☉"), "Zahl")
_P("08_doppelpulsar_masse", B, "08_doppelpulsar", "2,83", lambda p: M["PSR1913_M1"].wert + M["PSR1913_M2"].wert, "M☉",
   "Gesamtmasse", haengt_an=("PSR1913_M1", "PSR1913_M2"))
_P("08_doppelpulsar_ds", B, "08_doppelpulsar", "1,39",
   lambda p: K.ds((M["PSR1913_M1"].wert + M["PSR1913_M2"].wert) * K.M_sun) / 1e3, "km", "Gesamtmasse",
   haengt_an=("PSR1913_M1", "PSR1913_M2", "GM_sun"))
_P("08_doppelpulsar_fe", B, "08_doppelpulsar", "11,86", lambda p: _PSR(M["PSR1913_E"].wert) / _PSR(0.0), "1",
   "f(e) =", haengt_an=_PSR_G)
_P("08_ringdown_abklingzeit", B, "08_ringdown", "4,03",
   lambda p: K.ringdown_abklingzeit(M["GW150914_MF"].wert * K.M_sun, M["GW150914_F"].wert, M["GW150914_Q"].wert)
   * (1 + M["GW150914_Z"].wert) * 1e3, "ms", "Ringdown GW150914",
   haengt_an=("GW150914_MF", "GW150914_F", "GW150914_Q", "GW150914_Z", "GM_sun"))
_A("08_gedaechtnis", (D, "Gedächtnis GW150914: 2,7·10⁻²²"), "Zahl")
_P("08_gedaechtnis_gw150914", B, "08_gedaechtnis", "2,7·10⁻²²",
   lambda p: K.gedaechtnis_welle(M["GW150914_E_RAD"].wert * K.M_sun * K.c**2 * (1 + M["GW150914_Z"].wert),
                                 M["GW150914_DL"].wert * K.Mpc,
                                 K.richtungsfaktor_gedaechtnis(math.radians(E["IOTA_KANTE"].wert))), "1",
   "Gedächtnis GW150914", haengt_an=("GW150914_E_RAD", "GW150914_DL", "GW150914_Z", "IOTA_KANTE", "GM_sun"))
_A("08_ringdown_quelle", (D, "Quellsystem 272,4 Hz"), "Zahl")
_P("08_ringdown_quelle", B, "08_ringdown_quelle", "272,4",
   lambda p: K.ringdown_frequenz(M["GW150914_MF"].wert * K.M_sun, M["GW150914_F"].wert), "Hz", "Quellsystem",
   haengt_an=("GW150914_MF", "GW150914_F", "GM_sun"))
