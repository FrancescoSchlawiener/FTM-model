# -*- coding: utf-8 -*-
"""
Offene Prüfzahlen: Zahlen des Modelltexts, deren Kernformel fehlt (offen = "fehlt_formel") oder die im Modell nur auf
einer Forschungsrechnung beruhen (offen = "nur_forschung"). Ohne Rechnung; die Beschreibung nennt die benötigte Formel
bzw. die Rechnung. Je Satz und Beschreibung eine Aussage, je Zahl eine Prüfzahl. Steht die Formel im Kern, zieht die
Prüfzahl mit Rechnung nach zahlen_0X.py um.
"""
from .aussage import AUSSAGEN
from .modelltext import datei
from .pruefzahl import PRUEFZAHLEN

FF, NF = "fehlt_formel", "nur_forschung"


def _offen(b, satz_nr, stichwort, zahlen, art, beschreibung):
    """Eine Aussage (Satz mit stichwort) und je Zahl eine offene Prüfzahl."""
    a = f"{b}_offen_{satz_nr:02d}"
    AUSSAGEN.neu(a, (datei(b), stichwort), "Zahl")
    for i, z in enumerate(zahlen, 1):
        PRUEFZAHLEN.neu(f"{a}_{i}", b, a, z, None, None, stichwort, offen=art, beschreibung=beschreibung)


# ------------------------------------------------------------------ 01
_offen("01", 1, "Scherung", ("10⁻⁶",), NF,
       "Scherung der Schleifen unter 10⁻⁶ der Atemrate (Rechnung vakuum_richtungen)")

# ------------------------------------------------------------------ 02
_offen("02", 1, "Größenordnungen zu klein", ("29", "36",), NF,
       "Verhältnis der eigenen Faltung zu Bindung und Phase je Element (Rechnung elemente_analyse, Notiz 75: trägt nicht)")
_offen("02", 8, "Jeans-Kriterium", ("10²¹",), NF,
       "Zustandsschwelle δ_Z = 2kT/(m_A c²) mit homogener Kugel, Median über PubChem (Notiz 75, Rechnung elemente_analyse)")
_offen("02", 12, "Hauptreihe", ("0,7", "1,3",), NF,
       "Temperaturskala gegen die Kerntemperaturen der Hauptreihe [G] (Rechnung fusion_trichter)")

# ------------------------------------------------------------------ 03
_offen("03", 1, "Neigung winzig", ("10⁻¹¹", "10⁻¹²",), FF,
       "Neigung der großen Skalen (Filamente, Wände, Leerräume)")
_offen("03", 2, "Straffung des Raumes", ("2·10⁻⁶",), FF,
       "Straffung S(<r)/r im Zentrum der Milchstraße (Massenprofil mit Kern)")
_offen("03", 13, "schwaches mittleres äußeres Feld", ("0,5",), FF,
       "a₀-Lesung an SPARC je Υ_d gegen die Grenze a₀ ≤ 1,13·10⁻¹⁰ der Zelle")
_offen("03", 14, "schwaches mittleres äußeres Feld", ("0,5",), FF,
       "a₀-Lesung an SPARC je Υ_d mit äußerem Feld (Senkung bei Υ_d = 0,5)")
_offen("03", 17, "Hintergrundstrahlung mit BAO", ("0,05",), FF,
       "σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an")
_offen("03", 18, "Hintergrundstrahlung mit BAO", ("0,7",), FF,
       "σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an (0,7σ, ± 0,05)")
_offen("03", 27, "unabhängig von N*", ("1,05",), FF,
       "Linse der Schale außerhalb R_ta ((1 + w·Q)·M_b ab 1,05·R_ta)")
_offen("03", 30, "seit der Umkehr", ("50", "100",), FF,
       "Radiant der Bahn seit der Umkehr bei 100 kpc (Rechnung haufen_nach_naht)")
_offen("03", 31, "wirksamen Tiefe", ("0,5",), FF,
       "f = 1/S_T je Ebene")
_offen("03", 32, "wirksamen Tiefe", ("0,85",), FF,
       "f = 1/S auf allen Ebenen (0,5–0,85)")
_offen("03", 34, "Dichtespitzen", ("20", "10⁸",), NF,
       "Umkehr-z je Masse (Rechnung kollaps_nach_umkehr, noch nicht im Modell)")
_offen("03", 37, "Winkelmoden der Hintergrundstrahlung", ("42",), NF,
       "Grenze der Prüfbarkeit beim 1,5-fachen Durchmesser der letzten Streufläche (Rechnung zelle_form, 1,5·d_LSS)")
_offen("03", 38, "kosmische Streuung", ("17",), NF,
       "Senkung des Quadrupols bei L = 0,985·d_LSS (Würfel) aus forschung.winkelmoden_zelle mit d_LSS im eigenen Atmen "
       "(forschung.lauf_eigenes_atmen S2); frisch 17,0 % (bei 1,5·d_LSS 2,4 %)")
_offen("03", 40, "Durchmessers", ("1,5",), NF,
       "Grenze der Prüfbarkeit beim 1,5-fachen Durchmesser der letzten Streufläche (Rechnung zelle_form)")
_offen("03", 41, "Quadrupolleistung", ("17", "30", "33", "76", "38",), NF,
       "Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im "
       "eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0")
_offen("03", 42, "Paare im rechten Winkel", ("1,01·10⁻¹⁰", "1,17",), NF,
       "Einsatz der Kreispaare zellgrenze(f³, L_min = d_LSS)/(1 + Q) mit d_LSS im eigenen Atmen (forschung.lauf_eigenes_atmen "
       "S2) und Q = 5,36: frisch Würfel 1,013, sechseckig 1,170·10⁻¹⁰")
_offen("03", 52, "früh kollabierten Galaxien", ("3·10⁻⁵",), FF,
       "Streuung der Beschleunigungsrelation bei früh kollabierten Galaxien")
_offen("03", 55, "Massennullpunkt", ("1,9", "0,013", "0,05", "0,004", "50",), FF,
       "Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)")
_offen("03", 56, "kollabiert ein Bereich früher", ("0,3", "0,96",), NF,
       "Kollapszeit der Mode gegen die Kugel (Rechnung kollaps_nach_umkehr, noch nicht im Modell)")
_offen("03", 57, "Endantwort", ("0,54", "0,59",), FF,
       "Endantwort der Mode auf eine Massenänderung δM/M")
_offen("03", 60, "Rest der Mode senkt die Linse", ("19", "39",), NF,
       "Linse um isolierte Galaxien mit Restmode (Rechnung linsen_isoliert, noch nicht im Modell)")
_offen("03", 61, "Umkehrskala", ("0,09", "0,15", "0,27",), FF,
       "Skala, unter der das Strukturwachstum von der Regel abweicht, je z (Umkehrskala)")
_offen("03", 70, "Zählung je Radiant", ("0,15", "0,5",), FF,
       "nötige Genauigkeit von Ω_c/Ω_b, um Q = 2π − 1 von der Zählung je Radiant zu trennen")
_offen("03", 71, "Zählung je Radiant", ("1,0",), FF,
       "σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an")
_offen("03", 77, "keplerisch", ("100", "200",), NF,
       "Linse um isolierte Galaxien bei kleinem äußerem Feld (0,001·a₀): Radius, ab dem sie keplerisch wird (Rechnung linsen_isoliert, noch nicht im Modell)")
_offen("03", 78, "Gravitationsmoment", ("0,35", "0,68", "0,025", "0,02",), FF,
       "Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017)")
_offen("03", 79, "Verwindung", ("10⁶",), FF,
       "Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017): Verwindung")
_offen("03", 80, "Knotenlinie", ("90",), NF,
       "Knotenlinie bei einfallendem Gas (Rechnung warp, noch nicht im Modell)")
_offen("03", 83, "Vorzeichen der Abweichung", ("0,012",), NF,
       "Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)")
_offen("03", 84, "stärkeres Außenfeld", ("2,35", "1,44", "10,2", "179",), NF,
       "Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)")

# ------------------------------------------------------------------ 04

# ------------------------------------------------------------------ 07
_offen("07", 1, "nichtlokale Fortsetzung", ("36",), FF,
       "χ²/N 36 statt 16 an SPARC (10 Galaxien, Υ_d frei) der linearen nichtlokalen Fortsetzung (Laplace, |k|, v⁴ = (2/π)²GMa₀); benötigt: Regel der nichtlokalen Fortsetzung an Rotationskurven und ihre Anpassung")
_offen("07", 2, "nichtlokale Fortsetzung", ("16",), FF,
       "χ²/N 16 der örtlichen Regel (10 Galaxien, Υ_d frei); benötigt: Anpassung je Galaxie mit freiem Υ_d. Zugleich ausnahme: Koeffizient in c⁴R/16πG; „16 Punkte“ unter s ≈ 0,1 ist Prüfzahl 07_punkte_unter_s_trenn")
_offen("07", 5, "fallende Kurven", ("0,13",), FF,
       "d ln v/d ln r = −0,13 bis −0,27 für andere Verlustraten (Bahnfrequenz, √(a_N a₀)/c, Neigung mit Rückkopplung), benötigt: Regel mit diesen Raten; zudem λ_m = 0,01–0,13 Ω am Sonnenkreis (benötigt: Umlauffrequenz des gebundenen Bereichs, Parametersatz fehlt) und f = +0,13 an SPARC (benötigt: Anpassung von f)")
_offen("07", 6, "fallende Kurven", ("0,27",), FF,
       "d ln v/d ln r = −0,27 für andere Verlustraten; benötigt: Regel mit diesen Raten")
_offen("07", 7, "Newton-Grenze", ("1,18",), FF,
       "ν ≥ 1,18 überall bei ϖ-abhängigen Koeffizienten; benötigt: Rückkehr mit ϖ-abhängigem Γ/D")
_offen("07", 10, "drei mit Δχ²", ("650",), FF,
       "Δχ² ≈ 650 für drei Dimensionen (Δχ² ≥ 7 für zwei); benötigt: Gitter-Rückkehr ν_d(a_N) in d = 2, 3 und χ² der tiefen Bins mit Bootstrap-Kovarianz")
_offen("07", 11, "vollen Bose-Form", ("0,1",), FF,
       "Wurzel von der vollen Bose-Form erst unter s ≈ 0,1 zu trennen; benötigt: Trennschärfe beider Formen an den SPARC-Punkten")
_offen("07", 12, "harmonische Mode", ("0,05",), FF,
       "Band |λ| ≲ 0,05 (d ln v/d ln r = √(¼ + λ) − ½) mit nur λ = 0 für SL(2,ℤ); benötigt: Steigungsband aus SPARC und Spektrum des Laplace auf SL(2,ℤ)")
_offen("07", 13, "Normierung", ("10¹²⁴",), FF,
       "Normierung der Rückkehr als freie Energie eines Feldes")
_offen("07", 20, "Wellenform", ("6·10⁻⁶",), NF,
       "Wellenform der Rückkehr gäbe w ∝ √M, Milchstraße 6·10⁻⁶ (Rechnung gleichung_zeit, Welle trägt nicht; keine Kernformel)")
_offen("07", 21, "Schrittlänge fällt heraus", ("10⁶",), FF,
       "Rückkehranteil q unabhängig von der Schrittlänge h (Irrfahrt bis 10⁶·r₀)")
_offen("07", 22, "Raumvolumen", ("3·10³²",), FF,
       "Rückkehrzeit zur Kugel r* im Hubble-Volumen (Sonne 3·10³² Jahre); benötigt: Rückkehrzeit im Raumvolumen ⟨T⟩ ∝ 1/M")
_offen("07", 24, "Masseverlust eines Zwergs", ("30",), FF,
       "30 % Masseverlust eines Zwergs (Bahnen beider Formen nach 3 Mrd. Jahren auf 5·10⁻⁶ gleich) und Umlaufarbeit durch die Schicht 11–30 %; benötigt: zeitabhängige Bahnen beider Wahlen des Imaginärteils, Feldform-Löser der Scheibe. „bei 30 kpc“ ist Eingabe RING_R30")
_offen("07", 25, "Masseverlust eines Zwergs", ("5·10⁻⁶",), FF,
       "Bahnen beider Formen nach einem Massensprung auf 5·10⁻⁶ gleich; benötigt: zeitabhängige Bahnen")
_offen("07", 26, "Zustand der Regel", ("0,01",), FF,
       "λ_m = 0,01–0,13 Ω am Sonnenkreis; benötigt: Umlauffrequenz ω_orb des gebundenen Bereichs der Milchstraße (Dichtefaktor fehlt im Parametersatz)")
_offen("07", 27, "Zustand der Regel", ("0,13",), FF,
       "d ln v/d ln r = −0,13 bis −0,27 für andere Verlustraten (Bahnfrequenz, √(a_N a₀)/c, Neigung mit Rückkopplung), benötigt: Regel mit diesen Raten; zudem λ_m = 0,01–0,13 Ω am Sonnenkreis (benötigt: Umlauffrequenz des gebundenen Bereichs, Parametersatz fehlt) und f = +0,13 an SPARC (benötigt: Anpassung von f)")
_offen("07", 28, "Spannung gelesen", ("10⁻⁹",), FF,
       "Spannung des Wirbelstroms der Feldform, als Energiedichte gelesen")
_offen("07", 29, "Arbeit eines Umlaufs", ("30",), FF,
       "30 % Masseverlust eines Zwergs (Bahnen beider Formen nach 3 Mrd. Jahren auf 5·10⁻⁶ gleich) und Umlaufarbeit durch die Schicht 11–30 %; benötigt: zeitabhängige Bahnen beider Wahlen des Imaginärteils, Feldform-Löser der Scheibe. „bei 30 kpc“ ist Eingabe RING_R30")
_offen("07", 30, "Bulge", ("5·10¹⁰",), FF,
       "Scheibe 5·10¹⁰ M☉ (Bulge 10¹⁰ M☉) der Rechnung des Energieaustauschs über Sternbahnen; benötigt: Wirbelteil der Feldform und Bahnrechnung")
_offen("07", 31, "Bulge", ("10¹⁰",), FF,
       "Wirbelteil der Feldform und Bahnrechnung damit (Scheibe 5·10¹⁰ M☉, Bulge 10¹⁰ M☉)")
_offen("07", 32, "σ_z²", ("0,00",), FF,
       "Energieänderung 0,00 ± 0,03 % von σ_z² und f = 0,00 ± 0,03 ohne IC 2574; benötigt: Bahnrechnung in der Feldform, Anpassung von f")
_offen("07", 33, "σ_z²", ("0,03",), FF,
       "Boden unter 0,03·a₀ (95 %) an SPARC, f = +0,03, f = 0,00 ± 0,03, Energieänderung 0,00 ± 0,03 % von σ_z²; benötigt: Anpassung a = ν(a_N)·a_N + B, Anpassung von f, Bahnrechnung in der Feldform")
_offen("07", 34, "σ_z²", ("0,6",), FF,
       "Υ_d = 0,6 derselben Vergleiche; zudem „je Stern um höchstens 0,6 %“ (Energieaustausch über Sternbahnen in der Feldform, benötigt: Wirbelteil der Feldform und Bahnrechnung damit)")
_offen("07", 35, "Polytrope n = 3", ("1,32",), FF,
       "Überschussradius der Polytrope n = 3: 1,32·GM/c²; benötigt: δR = ∫(√g_rr − 1)dr über ein Lane-Emden-Profil")
_offen("07", 36, "Wachstum ändert sich", ("10⁻⁵",), FF,
       "Wachstum fσ₈ mit der Summe der Verluste (1 − w)(ν − 1)")
_offen("07", 37, "Feld der Störung", ("10⁻⁴",), FF,
       "Rückkehr mit dem Feld der Störung allein (a_N ≈ 10⁻⁴·a₀)")

# ------------------------------------------------------------------ 09
_offen("09", 2, "Gitterschwingung", ("10¹³",), FF,
       "Frequenz ω₃ = v_s/r der Gitterschwingung (Eingaben v_s, r fehlen)")
