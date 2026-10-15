# -*- coding: utf-8 -*-
"""
Beschlüsse (Klasse und Regeln in beschluss.py). So füllen: BESCHLUESSE.neu(kennung, Aussage, (Datei, Stichwort),
Grundlage, "JJJJ-MM-TT", "hinein" | "heraus"[, ersetzt=…, neue_grundlage=…]). Datum = Commit, mit dem die Aussage ins
Modell kam (git log -S).

OFFEN          Kennung eines Beschlusses → (Art, Grund): die Grundlage trägt den Satz nicht; die Hauptsitzung entscheidet,
               der Modelltext bleibt bis dahin. Art ∈ ARTEN_OFFEN.
VORSCHLAEGE    offene Beschluss-Vorschläge zur Liste „Widerspruch“ (Modell gegen frische Rechnung):
               Vorschlag(kennung, modell_ort, modell, frisch, grund, satz); grund ∈ GRUENDE, satz = vorgeschlagener Text.
"""
from .beschluss import BESCHLUESSE

REGISTERMODUL = True

ARTEN_OFFEN = {"trägt nicht", "[G]", "noch nicht im Modell", "ohne Eintrag im Forschungsstand"}
GRUENDE = {"Rundung", "Blickrichtung", "veraltet", "Parametersatz"}
D00, D01, D02, D03 = ("00-bauplan.md", "01-grundlagen.md", "02-vom-teilchen-bis-zum-stern.md",
                      "03-von-der-galaxie-bis-zum-kosmos.md")
D04, D05, D06, D07 = ("04-zeit-verschraenkung-syntropie-gedaechtnis.md", "05-zustandsraum.md", "06-torus-dynamik.md",
                      "07-feldgleichung.md")
D08 = "08-licht-und-trichter.md"
D09 = "09-elektromagnetismus-und-frequenzleiter.md"

# ------------------------------------------------------------------ Beschlüsse
BESCHLUESSE.neu("a0_standardwert_1_1", "Das Modell rechnet mit a₀ = 1,1·10⁻¹⁰ m/s² statt 1,2·10⁻¹⁰",
                (D00, "das Modell rechnet mit 1,1·10⁻¹⁰"),
                "SPARC mit realistischen Fehlern 1,04–1,19·10⁻¹⁰ (Υ_d = 0,5); Zelle ohne Kreispaare ≤ 1,13·10⁻¹⁰",
                "2026-09-29", "hinein")
BESCHLUESSE.neu("gewicht_aus_rueckkehr", "Gewicht der Schließungsregel w = 1 − e^{−√(θ̃τ_D)} aus derselben Rückkehr",
                (D03, "1 − e^{−√"), "Rechnungen kac_schliessung und gewicht_rueckkehr (im Modell); ersetzt die "
                "Maxwell- und die Preisach-Form", "2026-09-29", "hinein")
BESCHLUESSE.neu("tunnelbild_heraus", "Rückkehr als Tunneln (Schrödinger-Lesart D·G″ = Γ·G) als eigenes Bild",
                (D07, "Rückkehr als Tunneln"), "Umformung derselben Regel ohne neue Invariante (Thema "
                "rueckkehr_gleichung, gesperrt); tunnel_aussenfeld nie eingetragen", "2026-09-29", "heraus")
BESCHLUESSE.neu("q_dynamisch", "Q ist dynamisch, nicht topologisch", (D03, "Q ist dynamisch"),
                "schluss_oben (zwingend nur Q ≥ 5); q_form: ausgezeichnete Formen geben nur ν = 1,00–1,26 bzw. ganze "
                "oder rationale Zahlen", "2026-09-29", "hinein")
BESCHLUESSE.neu("zaehlung_je_radiant", "Q = 2π − 1 bei Zählung der Rückkehr je Radiant, als offene Lesart",
                (D03, "Zählung je Radiant"), "Herleitung des Autors; kac_schliessung: 1 − q = 1/2π nur mit der Setzung "
                "„ein Radiant“; camb_q: trennbar erst bei Ω_c/Ω_b auf 0,15–0,5 %", "2026-09-29", "hinein")
BESCHLUESSE.neu("elemente_faltung_zu_klein", "Die eigene Faltung von Atom und Kern ist 29–36 Größenordnungen zu klein",
                (D02, "29–36 Größenordnungen"), "elemente_analyse an PubChem und AME2020 (Notiz 75): Der Satz gibt das "
                "negative Ergebnis wieder", "2026-09-29", "hinein")
BESCHLUESSE.neu("zustandsschwelle_masse", "Ein Trichter bestimmt den Zustand seines Stoffs ab etwa 10²¹ kg",
                (D02, "etwa 10²¹ kg"), "elemente_analyse, Median über PubChem (Notiz 75); Formel fehlt im Kern",
                "2026-09-29", "hinein")
BESCHLUESSE.neu("wellenform_kein_mitatmen", "Eine Wellenform gäbe kein Mitatmen (w ∝ √M, Milchstraße 6·10⁻⁶)",
                (D07, "Milchstraße 6·10⁻⁶"), "gleichung_zeit: Welle trägt nicht; der Satz gibt das negative Ergebnis "
                "wieder", "2026-09-29", "hinein")

# ------------------------------------------------------------------ Beschlüsse aus Forschung mit offener Grundlage
BESCHLUESSE.neu("fluchtmass_dunkle_energie", "Fluchtmaß mit verzögerter Ankunft: Ω_X = 0,32, w(0) = −0,99, im eigenen Atmen 0,36, mit dem "
                "Vakuumteil 0,56 von 0,687; übrige 18 %, w₀ = −1,01, w_a = −0,61, Maximum z ≈ 0,01, fσ₈ 9–15 %",
                (D03, "Geschachtelt gezählt und ab der Übergabe I_F"),
                "de_betrag: die verzögerte Ankunft gibt BAO-χ² 224 (sofortige Verschiebung 14), H₀ 94; die Form ist widerlegt",
                "2026-09-30", "heraus")
BESCHLUESSE.neu("kehrseite_bindung", "Dunkle Energie zum Teil als Kehrseite der Bindung (Ankunft über die Rückkehr)",
                (D03, "Kehrseite der Bindung"), "de_kehrseite_bindung: in der Lesart des Ganzen ≤ 0,004 mit falschem Vorzeichen; "
                "de_betrag: die verzögerte Ankunft trägt nicht", "2026-09-30", "heraus")
BESCHLUESSE.neu("aussenfeld_frei_angepasst", "Frei angepasst verlangen die Kurven ein gemeinsames Feld senkrecht "
                "0,010 ± 0,009·a₀, gleichgerichtet 0,003 ± 0,003·a₀", (D03, "Frei angepasst verlangen die Kurven selbst"),
                "aussenfeld_kosmisch trägt nicht; die gemeinsame Anpassung von g_e steht nicht im Kern", "2026-09-30", "heraus")
BESCHLUESSE.neu("df2_df4_aussenfeld", "DF2 und DF4 verlangten σ_v²/σ_v,N² ≤ 2,4 bzw. 1,5; erst ein 11- bzw. 206-fach "
                "stärkeres Außenfeld", (D03, "DF2 und DF4 verlangten"),
                "gedaechtnis_umgebung, zwerge_aussenfeld; frisch 2,35 / 1,44 und 10,2- / 179-fach (andere Konvention "
                "der Fehlerspanne)", "2026-09-29", "hinein")
BESCHLUESSE.neu("scherung_schleifen", "Scherung der Schleifen kleiner als 10⁻⁶ der Atemrate", (D01, "Scherung"),
                "3c ordnet die Zahl vakuum_richtungen zu; der Satz steht seit dem ersten Commit des Modells",
                "2026-09-29", "hinein")
BESCHLUESSE.neu("hauptreihe_faktor", "Die Temperaturskala trifft die Hauptreihe auf den Faktor 0,7–1,3",
                (D02, "Faktor 0,7–1,3"), "fusion_trichter gegen Kerntemperaturen der Hauptreihe", "2026-09-29",
                "hinein")
BESCHLUESSE.neu("alpha_s_planck", "1/α_s ≈ 53 an der Planck-Skala, |q̃| = 1,2·10⁻²⁰, m_p +41 % / −30 %",
                (D02, "1/α_s ≈ 53"), "Rechnung mikroskala mit α_s(M_Z) und Quarkmassen", "2026-09-29", "hinein")
BESCHLUESSE.neu("umkehr_z_je_masse", "Die örtliche Regel setzt zuerst bei z ≈ 20 (10⁸ M☉) ein", (D03, "z ≈ 20 (10⁸"),
                "kollaps_nach_umkehr: Umkehr-z je Masse aus dem Spektrum", "2026-09-29", "hinein")
BESCHLUESSE.neu("kollaps_frueher", "Nach der Umkehr kollabiert ein Bereich in 0,3–0,96 der Zeit einer Staubkugel",
                (D03, "0,3–0,96"), "kollaps_nach_umkehr", "2026-09-29", "hinein")
BESCHLUESSE.neu("linse_300_kpc", "Linse bei N* = 2,4 um 3–7 % (200 kpc) und 19–39 % (300 kpc) unter der Relation",
                (D03, "19–39 %"), "linsen_isoliert: kumulativ gerechnet −0,5 / −4 / −22 / −44 % bei 200 / 300 / 500 / "
                "1000 kpc; im Modell die Zahlen des lokalen Ansatzes; KiDS nur aus einem Repo", "2026-09-29", "hinein")
BESCHLUESSE.neu("kaleidoskop_fenster", "Zelle oberhalb des 1,5-fachen Durchmessers nicht mehr prüfbar (27–42 Gpc)",
                (D03, "1,5-fachen Durchmessers"), "Gesamtdurchgang des Modells (Kaleidoskop-Fenster); keine Rechnung "
                "im Forschungsstand", "2026-09-29", "hinein")
BESCHLUESSE.neu("warp_gas_90_grad", "Einfallendes Gas dreht die Knotenlinie nur bis 90°, bei der gemessenen Rate nach "
                "0,14 Mrd. Jahren", (D03, "Knotenlinie nur bis 90°"), "warp_drehmoment, warp_amplitude; Rate aus "
                "WARP_PRAEZESSION", "2026-09-29", "hinein")
BESCHLUESSE.neu("ds_milchstrasse", "Milchstraße: Δs = 197,4 AE", (D03, "197,4 AE"),
                "frische Rechnung, Rundung, frisch 197,41 AE", "2026-09-30", "hinein")
BESCHLUESSE.neu("ds_r_milchstrasse", "Δs/R Milchstraße 6,3·10⁻⁸", (D03, "6,3·10⁻⁸"),
                "frische Rechnung, Rundung, frisch 6,26·10⁻⁸", "2026-09-30", "hinein")
BESCHLUESSE.neu("rueckstellbilanz", "ΔL_Raum/R = ΣΔs/R = 0,083", (D03, "ΣΔs/R = 0,083"),
                "frische Rechnung, Rundung, frisch 0,0833", "2026-09-30", "hinein")
BESCHLUESSE.neu("m_b_q", "M_b,Q ≈ 4,45·10¹³ M☉ Baryonen (heute)", (D03, "4,45·10¹³ M☉"),
                "frische Rechnung, Rundung, frisch 4,45·10¹³", "2026-09-30", "hinein")
BESCHLUESSE.neu("umlauf_an_l", "Sonne: Umlauf an ℓ 84,5 Mio. Jahre", (D07, "Umlauf an ℓ 84,5 Mio"),
                "frische Rechnung, Rundung, frisch 84,4965", "2026-09-30", "hinein")
BESCHLUESSE.neu("pal14_newton", "Newton-Wert Pal 14 0,665 km/s", (D03, "(0,665; 5,09; 7,15 km/s)"),
                "frische Rechnung, Rundung, frisch 0,665", "2026-09-30", "hinein")
BESCHLUESSE.neu("ds_erdmasse", "Δs = GM/(3c²) = x · 1,4783 mm", (D01, "x · 1,4783 mm"),
                "frische Rechnung, Rundung, frisch 1,4783", "2026-09-30", "hinein")
BESCHLUESSE.neu("uhr_neutronenstern", "Neutronenstern: 69,5 Tage", (D04, "Neutronenstern 69,5 Tage"),
                "frische Rechnung, Rundung, frisch 69,5", "2026-09-30", "hinein")
BESCHLUESSE.neu("gw170817_laenge", "GW170817 590 m", (D08, "GW170817 590 m"),
                "frische Rechnung, Rundung, frisch 589,5", "2026-09-30", "hinein")
BESCHLUESSE.neu("v_erde", "Erde: 29,78 km/s", (D02, "Erde: 29,78 km/s"),
                "frische Rechnung, Rundung, frisch 29,785", "2026-09-30", "hinein")
BESCHLUESSE.neu("g_roter_riese", "Roter Riese: 0,0085", (D02, "| 0,0085 |"),
                "frische Rechnung, Rundung, frisch 0,00846", "2026-09-30", "hinein")
BESCHLUESSE.neu("einrasten_2_pi", "erstes Einrasten bei (N − 1)k ≈ 0,57–0,63 W", (D06, "0,57–0,63 W"),
                "frische Rechnung 0,569–0,630 (Kuramoto-Einrasten)", "2026-09-30", "hinein")
BESCHLUESSE.neu("sonnenrand", "Sonnenrand: 0,8756″ + 0,8756″ = 1,7512″ (R☉ = 6,957·10⁸ m)", (D08, "1,7512″"),
                "frische Rechnung, Parametersatz R☉, frisch 1,7512″", "2026-09-30", "hinein")

# ------------------------------------------------------------------ Beschlüsse dieses Laufs (frische Rechnung entscheidet)
BESCHLUESSE.neu("fluchtmoden_faser", "Fluchtmoden gehen als Operator je Raumpunkt in die innere Richtung über und tragen wie der "
                "Vakuumteil mit w = −1 bei (w_a < 0); Zahl der Zyklen je Raumpunkt offen", (D03, "je Raumpunkt in die innere Richtung über"),
                "de_betrag: sofortige Verschiebung gibt nahezu feste Form, w = −1,03 / −1,10 / −1,18 bei z = 0 / 1 / 2, BAO-χ² 14 "
                "(Λ 18); der Betrag ist an θ* angepasst und steht nicht im Modell", "2026-09-30", "hinein",
                ersetzt="fluchtmass_dunkle_energie", neue_grundlage="de_betrag: verzögerte Ankunft χ² 224, je Faltung χ² 669, "
                "sofortige Verschiebung 14")
BESCHLUESSE.neu("vakuumteil_rueckstellung", "Der Vakuumteil ist die Rückstellung des Ganzen; die Bindung macht Fluchtmoden inaktiv",
                (D03, "Der Vakuumteil ist die Rückstellung des Ganzen. Die Bindung"), "Aussage des Modells (Rückstellung, Bindungsoperator); "
                "ersetzt die Lesart als Kehrseite der Bindung", "2026-09-30", "hinein", ersetzt="kehrseite_bindung",
                neue_grundlage="de_betrag: verzögerte Ankunft trägt nicht")
BESCHLUESSE.neu("doppelungen_aufgeloest", "Jede Aussage steht nur an einer Stelle: Q-Herkunft, Zelle und Schwelle, Zahlen der Leiter, "
                "Erde und Sonne, Verstärkung an r*, Resonanzen", (D00, "Was noch fehlt"),
                "Doppelungsprüfung der Bausteine 00–09; 00 nennt je Baustein nur Inhalt und Lücke", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_00_09", "Frische Werte statt Rundung: g der Erde 9,82, α_G = 5,9·10⁻³⁹, Eddington 72,6 M_* = 134,6 M☉, Ω_s(1090) = 2,6·10⁻¹⁰, "
                "Streuung 0,27–0,30, σ(Q) 0,065 (1,3σ), fσ₈ 1,2·10⁻³, Kreislauf 2,6·10¹¹ Jahre, 567-fach", (D03, "0,27–0,30 dex"),
                "frische Läufe des Prüfstands", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_sparc", "SPARC: 1,06 ± 0,085 ± 0,11, Spanne bei Υ_b = 0,7 fest 0,96–1,19, Drittel 0,75 / 1,51 / 1,06 mit "
                "Anstieg +0,09 ± 0,08 (nicht gesichert), 74 Galaxien, Draco 1,2 km/s, Haufenstoß bis 1,8 Mrd. Jahre",
                (D03, "0,75 / 1,51 / 1,06"), "frische Läufe des Prüfstands (SPARC-Daten)", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_zwerge", "DF2 und DF4: σ_v²/σ_v,N² ≤ 2,35 bzw. 1,44, Außenfeld 10,2- bzw. 179-fach", (D03, "10,2- bzw. 179-fach"),
                "frische Rechnung gedaechtnis_umgebung, zwerge_aussenfeld (Konvention der Fehlerspanne)", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_universum", "Masse im eigenen Horizont 7,4·10²² M☉; Quadrupol des Sonnenorts 5,6·10⁻³⁰ s⁻²", (D03, "7,4·10²² M☉"),
                "frische Rechnung (Horizont 46,23 Mrd Lj; Gezeit des geglätteten Felds)", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_06", "Einrasten 0,57–0,63 W; Kupferdraht: 0,0010–0,6 Hz vor der Falte, Ringfrequenz 0,20–1,5 Hz (unterkritisch)",
                (D06, "0,0010–0,6 Hz"), "frische Rechnung (Kuramoto-Einrasten, Kupferdraht-Raster)", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_05", "Fluchtmoden höchstens 567-fach; Fluchtmoden gehen je Raumpunkt in die innere Richtung über",
                (D05, "je Raumpunkt in die innere Richtung"), "frische Rechnung; de_betrag", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_07", "Rand bei 0,66·r*·Λ/ϖ_R; Resonanz 5/4 bei 0,94 r*; Gleichstand 3,4 R_d; Ringfaktor 2,3–3,0",
                (D07, "0,66·r*·Λ/ϖ_R"), "frische Rechnung", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_08_09", "GP-B: 6604; Mitführung 40,8 mas/Jahr längs der Erdachse, zum Leitstern IM Peg 39,1; α/α_G = 1,23555·10³⁶",
                (D08, "vom Kreisel zum Leitstern IM Peg gesehen 39,1"), "frische Rechnung (Bahnmittel des Vektors, cos δ des Leitsterns)",
                "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_01_04", "Doppelte Tabelle Erde und Sonne in 01 gestrichen (steht in 02); Fluchtmaß trägt bei",
                (D01, "trägt als drittes, positives Glied bei"), "Doppelungsprüfung", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_04", "Kreislauf: mit konstanter Dunkler Energie ab 2,6·10¹¹ Jahren verletzt; Rate endet beim Vakuumteil",
                (D04, "ab etwa 2,6·10¹¹ Jahren"), "frische Rechnung", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_02", "g der Erde als GM/R² 9,82; Eddington 72,6 M_* = 134,6 M☉; α_G = 5,9·10⁻³⁹", (D02, "9,82"),
                "frische Rechnung", "2026-09-30", "hinein")
BESCHLUESSE.neu("frisch_09", "α/α_G = 1,23555·10³⁶; Umlauf bei r* ohne Zahl der Verstärkung (steht in 08)", (D09, "1,23555·10³⁶"),
                "frische Rechnung", "2026-09-30", "hinein")

OFFEN = {
    "df2_df4_aussenfeld": ("[G]", "Urteile der Zwerge offen, hängen an σ-Messwerten [G]"),
    "scherung_schleifen": ("trägt nicht", "einzige genannte Grundlage vakuum_richtungen trägt nicht"),
    "hauptreihe_faktor": ("[G]", "Vergleich mit Literaturwerten [G]; nach den Arbeitsregeln nur Forschungsstand"),
    "alpha_s_planck": ("[G]", "Eingänge α_s(M_Z) und Quarkmassen [G]; Lauf von α_s nicht im Kern"),
    "umkehr_z_je_masse": ("noch nicht im Modell", "Forschungsstand: Umkehr-z je Masse noch nicht im Modell, Spektrum [G]"),
    "kollaps_frueher": ("noch nicht im Modell", "Forschungsstand: kollaps_nach_umkehr noch nicht im Modell"),
    "linse_300_kpc": ("noch nicht im Modell", "linsen_isoliert noch nicht im Modell; kumulativ andere Zahlen"),
    "kaleidoskop_fenster": ("ohne Eintrag im Forschungsstand", "Zahl 1,5 und 42 Gpc ohne Rechnung im Forschungsstand"),
    "warp_gas_90_grad": ("[G]", "Rate hängt an WARP_PRAEZESSION [G]"),
}


# ------------------------------------------------------------------ Vorschläge zur Widerspruchsliste
class Vorschlag:
    """Offener Beschluss-Vorschlag: Modelltext gegen frische Rechnung, mit Grund und vorgeschlagenem Satz."""
    __slots__ = ("kennung", "modell_ort", "modell", "frisch", "grund", "satz")

    def __init__(self, kennung, modell_ort, modell, frisch, grund, satz):
        if grund not in GRUENDE:
            raise ValueError(f"{kennung}: Grund {grund!r}")
        if modell_ort is not None and (len(modell_ort) != 2 or not modell_ort[0].endswith(".md")):
            raise ValueError(f"{kennung}: modell_ort = (Datei, Stichwort) oder None (nur Register)")
        if not satz:
            raise ValueError(f"{kennung}: vorgeschlagener Satz fehlt")
        self.kennung, self.modell_ort, self.modell, self.frisch = kennung, modell_ort, modell, frisch
        self.grund, self.satz = grund, satz


VORSCHLAEGE = {}


def _v(kennung, ort, modell, frisch, grund, satz):
    if kennung in VORSCHLAEGE:
        raise ValueError(f"Vorschlag doppelt: {kennung}")
    VORSCHLAEGE[kennung] = Vorschlag(kennung, ort, modell, frisch, grund, satz)
