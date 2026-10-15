# -*- coding: utf-8 -*-
"""
Einordnung aller Archivdateien (Klasse und Regeln in katalog_eintrag.py): jedes Skript, dazu jede verwaiste Ausgabe
(.txt, .npy ohne gleichnamiges Skript) und jedes übersetzte Modul ohne Quelle. Je Datei eine Zeile, alphabetisch.

  übernommen   ziel = Kennung der Formel, Prüfzahl oder Gegenprobe, in der Logik oder Zahl jetzt steht; der Grund nennt
               weitere Kennungen. Ist die Zahl nur als offene Prüfzahl erfasst (Formel fehlt), ist das Ziel diese Prüfzahl.
  interessant  ziel = Kurztitel der Forschungszeile; die Rechnung ist ein Versuch, der noch nicht im Kern steht.
  nutzlos      Grund: ersetzt durch …, Vorstufe von …, Lesart aufgegeben, verwaiste Ausgabe.
"""
from katalog_eintrag import AUSGABE, KATALOG

_N = KATALOG.neu
Ü, I, N = "übernommen", "interessant", "nutzlos"
_VERGESSEN = "Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell"
_SCHWELLE = "Wert der Schwelle"

# ------------------------------------------------------------------ Skripte
_N("a0_lesungen", I, "a₀-Lesungen",
   "Lesungen im Modell (03) im Kern (a0_profiliert, profil_chi2), Prüfzahlen 03_sparc_a0_profil_wert, "
   "03_sparc_a0_profil_stat, 03_sparc_a0_varianten_unten (a0_je_galaxie und a0_an_punkten ohne Profil)")
_N("atmen_oder_zelle", N, None, "Vorstufe (Lesarten a₀ ∝ H gegen a₀ ∝ 1 + z); ersetzt durch a0_frueh_relativ, a₀ seit "
   "z ≈ 10 fest")
_N("atmen_ohne_skala", Ü, "zustand_vakuum",
   "Zustandsgleichung des Vakuumteils 1 + w_s = −⅔·d ln a₀/d ln a im Kern; der Fixpunkt H ↔ a₀ mit a₀ = cH/2π ist als "
   "Lesart aufgegeben")
_N("atmen_w", N, None, "Vorstufe (Λ_s = (π²/15)·H²/c²); ersetzt durch lambda_vakuum Λ_s = π⁴/(15Λ²)")
_N("ausgedehnte_massen", Ü, "zusatzpotential",
   "Zusatzpotential v_flach²·ln ν(s), Bose-Leiter (bose_leiter) und gelesene Dichte (dichte_gelesen) im Kern")
_N("aussenfeld_kosmisch", I,
   "Äußeres Feld als Einbettung in die kosmische Rückkehr, g_e = (ln q_kosm)²·a₀ = 0,029·a₀",
   "trägt nicht (Satz aus dem Modell genommen); Anstieg je Drittel im Kern (a0_je_drittel, 03_sparc_a0_drittel_anstieg)")
_N("aussenfeld_sparc", I, _SCHWELLE,
   "Senkung der gelesenen Schwelle durch das äußere Feld (7–14 % / 22–31 %) im Kern (a0_senkung_aussenfeld), "
   "Prüfzahlen 03_aussenfeld_senkrecht_gew bis 03_aussenfeld_gleich_ungew")
_N("axiom_mass", I, "Axiome und Maß", "Aussage im Modell (05, 07) als Beziehung; Rechnung (Jacobi-Determinante, "
   "Separatrixfläche) nicht im Kern")
_N("bistabilitaet", Ü, "preisach_ausgang",
   "Preisach-Zustand als Folge der nicht ausgelöschten Extremwerte im Kern; VO₂-Mikrobrücke als offene Prüfzahlen "
   "06_vo2_hopf_spannung, 06_vo2_hopf_frequenz, 06_vo2_torus_spannung")
_N("camb_q", N, None, "ohne Kern und ohne Ergebnis im Modell (CMB-Spektren Q = 5,364 gegen 2π − 1)")
_N("carter_herleitung", Ü, "carter_faktor", "Vorfaktor 2,944 als Prüfzahl 02_carter; Grenzmasse im Kern "
   "(grenzmasse_zaehlregel); Gleichstand von Gas- und Strahlungsdruck im Kern (eddington_quartik, "
   "02_gleichstand_druck_mstern, 02_gleichstand_druck_msonne)")
_N("compton_carter", N, None, "Vorstufe von carter_herleitung; ersetzt durch carter_faktor")
_N("compton_faltung", Ü, "faltungsrate",
   "α_G = 3Δs/λ̄ (alpha_G) und Γ_G = α_G·ω_C im Kern, Prüfzahl 09_faltungsrate_proton")
_N("de_bilanz_glieder", I, "Bilanz aller Glieder", "trägt nicht; kein Glied des Kerns")
_N("de_fluchtmoden", I, "DE als verdrängte Fluchtmoden",
   "trägt nicht in der Zeit (verzögerte Ankunft aus dem Modell genommen); Maß je Bereich im Kern (fluchtmass_bereich)")
_N("de_fluchtmoden_erneuerung", I, "Fluchtmaß mit Erneuerung", "trägt nicht; Erneuerung im Kern nur örtlich "
   "(erneuerung)")
_N("de_fluchtmoden_hierarchie", I, "Fluchtmaß hierarchisch", "trägt nicht")
_N("de_fluchtmoden_mass", Ü, "fluchtmass_bereich",
   "Maß je Bereich X_b = ln(μ_Regel/μ_Newton) (Milchstraße ln 284) im Kern, Prüfzahlen 03_fluchtmass_verhaeltnis, "
   "03_fluchtmass_ln; Übergang ins Atmen offen (atemrate_eff ohne Lauf Ẋ_an)")
_N("de_fluchtmoden_verzoegert", I, "Fluchtmaß verzögert über die Rückkehr",
   "trägt nicht in der Menge; Rate H_eff = H_b + ⅓·dX_an/dt im Kern (atemrate_eff), Lauf Ẋ_an fehlt")
_N("de_invariante", I, "Dunkle Energie als Invariante der Grundwerte", "Rechnung, noch nicht im Kern")
_N("de_kanal_atmen", I, "Kanal ins Atmen", "Aussage im Modell (03, 05); Irrfahrt auf [0, Λ] und Kanalwahl nicht im Kern")
_N("de_kehrseite_bindung", I, "DE als Kehrseite der Bindung", "trägt nicht")
_N("de_laufzeit", I, "Laufzeit des Übergangs", "Kandidaten K0–K4 gerechnet; Laufzeit fehlt im Kern")
_N("de_phasenseite", I, "Dunkle Energie über die Phasenseite des Ganzen", "trägt nicht")
_N("de_ratenschwelle", I, "Dunkle Energie über die Ratenschwelle", "trägt nicht (hängt an frei gewählter Masse)")
_N("de_rest_anstieg", Ü, "lauf_eigenes_atmen",
   "Schleife Atmen → Bildung → Maß als Forschungsformel (lauf_eigenes_atmen, fluchtmass_je_baryon, ankunft_rueckkehr); "
   "Zahlen nur Forschung, nicht im Modell (verzögerte Ankunft widerlegt)")
_N("de_uebergabe_anteil", I, "DE-Rest als Übergabe der Mode", "trägt nicht")
_N("de_zyklen", I, "Dunkle Energie je Zyklus der Zelle", "trägt nicht als Herleitung")
_N("dicke_alter_scheiben", Ü, "07_offen_31_1",
   "Energieaustausch über Sternbahnen (Scheibe 5·10¹⁰, Bulge 10¹⁰ M☉) als offene Prüfzahlen 07_offen_30_1, "
   "07_offen_31_1; Wirbelteil der Feldform und Bahnrechnung fehlen im Kern")
_N("dicke_alter_scheiben_gegenprobe", N, None, "Gegenprobe von dicke_alter_scheiben ohne eigene Zahl im Modell")
_N("ebenen_projektion", Ü, "glaettung_innen",
   "harmonische Fortsetzung P_ϖ = e^{−|k|ϖ} im Kern; Tiefen je Körper (Mond 482 AE) nicht im Modell")
_N("ebenengrenze_zwerge", I, "Ebenengrenze der Hierarchie als Umgebungsgröße der Zwerge", "trägt nicht")
_N("eine_struktur", Ü, "temperatur_innen",
   "Temperatur je Periode kT = ħ/(Periode 2Λ/c) im Kern; CHSH mit Oberschwingungen im Kern (chsh_oberschwingungen, 04_chsh_zwei_harmonische)")
_N("elemente_analyse", Ü, "06_siede_schmelz_min",
   "Siede- zu Schmelztemperatur aus PubChem (pubchem), Prüfzahlen 06_siede_schmelz_min, 06_siede_schmelz_max; "
   "Faltungsverhältnisse je Element 02_offen_01_1, 02_offen_08_1 (nur Forschung), Spaltbarkeit im Kern "
   "(troepfchen_anpassung, spaltbarkeit, 06_spaltbarkeit_max), Zustandsschwelle 02_zustand_eisen_erde")
_N("em_phasenseite", Ü, "moment_magnetisch", "μ = Q_el·c·a_J mit g = 2, Prüfzahl 09_bohr_magneton")
_N("engstelle_global", Ü, "faltung_horizont",
   "Δs_H = r_H/6 im Kern, Prüfzahlen 04_ds_H_heute, 04_ds_H_zu_Lambda; Kopplungen Λ ∝ Δs_H nicht im Modell")
_N("engstelle_rate", I, "Engstelle und Schwarze Löcher", "trägt nicht")
_N("feldgleichung", N, None, "frühe Fassung der Feldgleichung (Variante B, Maxwell-Kern κ₀e^{−λt}, δa = q·c²/L); "
   "ersetzt durch kern_zeit und a_regel")
_N("folgerungen", N, None, "Vorstufe (drei Radien, Einzeltrichter gegen kollektiv); ersetzt durch r_stern_leiter")
_N("formeln_andere_gebiete", Ü, "entropie_ruecklauf",
   "Shannon-Entropie der Rücklaufzahl ln ν + χ im Kern; Watson-Schranke 1,516 im Kern (watson_rueckkehr, 07_watson_drei)")
_N("forschung_uebergang", N, None, "Vorstufe von feldgleichung (Sprung bei r*, a₀ aus dem Atmen); Lesart aufgegeben")
_N("freiheitsgrad", N, None, "ohne Kern und ohne Ergebnis im Modell (Spielzeug-Boltzmann, Notiz 26)")
_N("fruehe_galaxien_a0", I, _VERGESSEN,
   "Preisach-Abweichung frisch kollabierter Galaxien als offene Prüfzahlen 03_offen_55_1 bis 03_offen_55_5; "
   "Antwort d ln a/d ln M = 1 − χ/2 im Kern (flussdimension)")
_N("fruehes_wachstum", N, None, "ohne Ausgabe; derselbe Ansatz ist in kollaps_nach_umkehr gerechnet")
_N("ftm_gesamtlauf", N, None, "ersetzt durch lauf.py (ein Lauf, Bericht bericht.md); Formeln alle im Kern")
_N("ftm_gesetze", Ü, "B_tully_fisher_grenzfall",
   "alle Gesetze im Gesetz-Register (pruefung/gesetze_A, _B, _C, _E, _M, _Z; A 24, B 12, C 16, Z 49, M 121, E 9); M22 und M23 als Forschungszeilen "
   "(Spiegel M ↔ m_Pl²/M, Rückkehrgleichung bei E = 0)")
_N("ftm_kern", N, None, "reiner Lader für simulation/kern (alter Inhalt: ftm_kern_alt)")
_N("ftm_kern_alt", Ü, "ds", "gesicherter alter Inhalt von ftm_kern (siehe ftm_kern)")
_N("ftm_pruefstand", Ü, "08_leiterzahl",
   "Alle 389 Prüfungen zugeordnet (Prüfzahl, Gesetz, offene Prüfzahl oder nicht mehr im Text). Entfallen, weil "
   "die Zahl nicht im Modelltext steht, alle aus ε_L (Prüfzahl 08_leiterzahl): 1/ε_L = 61,7 Dekaden, "
   "M_min = ħa₀/c³ = 4,3·10⁻⁷⁰ kg, Proton M/M_min = 3,9·10⁴²")
_N("fusion_trichter", Ü, "grenzmasse_zaehlregel",
   "Grenzmasse n = 3 und Decke kT_max im Kern (grenzmasse_zaehlregel, decke_faktor, tiefe_faktor), Prüfzahlen "
   "02_grenzmasse_zaehlregel, 02_decke_faktor; Temperaturskala der Hauptreihe 02_offen_12_1 (nur Forschung)")
_N("galaxie_als_ganzes", Ü, "massenmittel_A_scheibe",
   "System im äußeren Feld: A und m = 1 im Kern (massenmittel_A_plummer, einseitigkeit_m1), Prüfzahlen 03a_scheibe_A, "
   "03a_m1_scheibe_30, 03a_plummer_A")
_N("gedaechtnis_topologie", I, "Gedächtnis als Topologie",
   "Preisach-Form im Kern (gewicht_w_preisach); Zählung der Umkehrungen und N* als Radiant nicht im Kern")
_N("gedaechtnis_umgebung", Ü, "sigma_plummer_jeans",
   "Plummer mit Jeans im Kern, σ_N als Prüfzahlen 03_druck_newton_pal14, 03_druck_newton_ngc2419, "
   "03_druck_draco_isoliert; Gedächtnis-Lesart trägt nicht (03_offen_83_1, 03_offen_84_1 nur Forschung)")
_N("gegenlauf_dunkle_energie", N, None, "Lesart der Dunklen Energie vor dem Fluchtmaß (Bausteine aus "
   "rueckstellung_zeitform); ohne Ergebnis im Modell")
_N("gemeinsame_innere_richtung", I,
   "Eine oder viele innere Richtungen (Bündel je Falte gegen gemeinsamen Knoten bei Λ)",
   "je Falte eine Halbgerade im Modell (07) als Beziehung; Knotenrechnung nicht im Kern")
_N("geometrische_muster", Ü, "flussdimension",
   "d_eff = 3 − χ im Kern; übrige Muster (Kreis als Unendlich-Eck, Kugel und Würfel) ohne Ergebnis im Modell")
_N("gesamtbild", N, None, "Übersicht (Notiz 29) aus Kernwerten, ohne eigenes Ergebnis im Modell")
_N("gewicht_preisach", Ü, "gewicht_w_preisach", "θ̃ = kleinste Rate der Geschichte im Kern; Maxwell-Form daneben "
   "nur in forschung/ (gewicht_w_maxwell)")
_N("gewicht_rueckkehr", Ü, "schalen_dauer",
   "w = 1 − e^{−√x} im Kern (gewicht_w); Dauer w 0,9 → 0,1 als Prüfzahlen 03a_schale_z_c_34 bis 03a_schale_z_c_05; fσ₈ und Schale außerhalb der Naht im Kern (wachstum_uebergabe, gewicht_schale_aussen, 03_fsigma8_uebergabe, 03_w_aussen_133)")
_N("gewicht_scheibenfeld", Ü, "a_scheibe_ebene",
   "a = ν(√(g_R² + f·(2πGΣ)²))·g_R im Kern; Anpassung von f an SPARC im Kern (gewicht_f, gewicht_f_tafel, sigma_aus_rotation), "
   "Prüfzahlen 07_gewicht_f_unten, 07_gewicht_f_oben")
_N("gewicht_scheibenfeld_gegenprobe", N, None, "Gegenprobe von gewicht_scheibenfeld ohne eigene Zahl im Modell")
_N("glaettung_stetig", Ü, "gezeitendichte", "ρ_J = 3𝒯/(4πG) im Kern (Einsatz der Glättung über ρ − ρ_J); Quadrupol des Restfelds im Kern (quadrupol_glaettung, 07_quadrupol_glaettung)")
_N("gleichung_ganzes", Ü, "erneuerung",
   "Erneuerung u(0) = g_N + u(ϖ_R) im Kern; nichtlokale Lesart als offene Prüfzahlen 07_offen_01_1, 07_offen_02_1")
_N("gleichung_zeit", Ü, "rueckkehr",
   "û = e^{−√((p + Γ + θ)τ_D)} im Kern; Wellenform trägt nicht (07_offen_20_1 nur Forschung)")
_N("gravitationswellen", Ü, "hulse_taylor",
   "Hulse-Taylor, ISCO (isco_kerr), Ringdown (ringdown_frequenz), Schwellenfrequenz (schwellenfrequenz_gw) im Kern; "
   "Abklingzeit im Kern (ringdown_abklingzeit, 08_ringdown_abklingzeit), Gedächtnis 08_gedaechtnis_gw150914")
_N("halbgerade_hierarchie", I, _VERGESSEN, "Halbgerade der Hierarchie D·G″ = Γ·G mit Ebenengrenze; nicht im Kern")
_N("halbgerade_spitze", Ü, "nu_aus_modul", "ν = 1/(1 − |e^{2πiτ}|) im Kern, Prüfzahl 05_nu_tau_zwei; Maass-Formen "
   "ohne Ergebnis im Modell")
_N("harvey_nachlauf", I, _VERGESSEN, "Quelle nicht erreichbar; nicht im Kern")
_N("haufen_nach_naht", I, "Haufen nach der Naht", "Radiant seit der Umkehr als offene Prüfzahlen 03_offen_30_1, "
   "03_offen_30_2 (Bahnrechnung nicht portiert); Rest im Haufen im Kern (haufen_rest, 03_haufen_rest_1mpc)")
_N("herkunft_q", N, None, "Vorstufe von kac_schliessung (Herkunft des kosmischen q); ersetzt durch kac_zeit_schluss")
_N("hierarchie", N, None, "Varianten der Einzeltrichter (Notiz 13); ersetzt durch glaettungslaenge (Glättung über "
   "die Elterngezeit)")
_N("holomorphie_rueckkehr", Ü, "glaettung_innen", "Halbgruppe P_a·P_b = P_{a+b} der harmonischen Fortsetzung im Kern")
_N("hopf_knoten", I, "Hopf-Knoten", "Rechnung, noch nicht im Modell (Thema hopf_knoten)")
_N("horizont_kreislauf", Ü, "horizontentropie",
   "S_H und S_Λ im Kern (entropie_vakuum), Mindestrate (kreislauf_mindestrate); Bedingung an Ḣ im Kern (kreislauf_versiegen, weltzeit), "
   "Prüfzahl 04_kreislauf_versiegen")
_N("identifikation_galaxie", Ü, "05_resolvente_flucht", "Schranke 568 der Resolvente im Kern (resolvente_flucht, "
   "05_resolvente_flucht); Rückkehr als Abel-Mittel im Kern (schleifensumme)")
_N("info_geometrie", Ü, "fisher_metrik", "g = ν(ν − 1) und S_n = ln ν + χ (entropie_ruecklauf) im Kern")
_N("innere_richtung_zeit", Ü, "rueckkehr_geschlossen",
   "geschlossene Rückkehr als thermischer Propagator der Periode 2Λ/c im Kern (temperatur_innen); Bose-Form aus T_I "
   "und Unruh-Lesart tragen nicht")
_N("invarianten_ebenen", Ü, "d_LSS_beobachter", "L/d_LSS für einen Beobachter bei z = 1 als Prüfzahl "
   "03_zelle_beobachter_z1")
_N("kac_schliessung", Ü, "gewicht_w",
   "Form w = 1 − e^{−√(θ̃τ_D)} und ⟨T⟩ = 2Λ/c (kac_zeit_schluss) im Kern")
_N("kern_chi", N, None, "Vorstufe der Kopfgleichung (a = a_N + χ·(Δs/r)·3c²/r*); ersetzt durch kern_zeit")
_N("kern_gedaechtnis", I, "Gedächtniskern", "Maxwell-Kern der Kopfgleichung; Kern im Modell jetzt aus der Rückkehr "
   "(kern_zeit)")
_N("kern_nome", N, None, "Vorstufe (Nome des Kerns gegen Rückkopplungssumme); ersetzt durch kern_zeit")
_N("kern_schleife", N, None, "Vorstufe (Kern als Schleife auf der Kreisbahn); ersetzt durch kern_zeit")
_N("kollaps_nach_umkehr", I, "Kollaps nach der Umkehr",
   "Umkehr-z und Kollapszeit als offene Prüfzahlen 03_offen_34_1, 03_offen_56_1 (nur Forschung)")
_N("kontakt_atmen", Ü, "rc_kette_nu",
   "RC-Kette mit Rückführung im Kern, Prüfzahlen 06_rc_nu, 06_rc_wurzel; Kontaktform trägt nicht")
_N("kontraktion_atemrate", Ü, "hintergrund_modell",
   "Hintergrund mit Vakuumteil im Kern, H_s (H_vakuum) als Prüfzahl 04_endrate_H_s, Mindestrate 03_mindestrate")
_N("kopfgleichung_kern", Ü, "kern_zeit", "Zeitkern aus der Rückkehr und Dämpfung bei der Bahnfrequenz "
   "(kern_bahnfrequenz), Prüfzahl 07_kern_sonnenkreis")
_N("kosmologie", N, None, "Lesart a₀(z) = cH(z)/2π aufgegeben (a₀ seit z ≈ 10 fest, a0_frueh_relativ); CMB/BAO-"
   "Rechnung ohne Ergebnis im Modell")
_N("kovariante_e0", Ü, "gezeit_betrag", "Invariante der Gezeit 𝒯 = √(3/2)·|Ē^TF| und Glättungslänge "
   "(glaettungslaenge) im Kern")
_N("kupfer_selbstpulsieren", Ü, "kupfer_hopf_abstand",
   "Hopf-Grenze, Ringfrequenz, Zeitverhältnis im Kern (kupfer_hopf, kupfer_ringfrequenz, kupfer_zeitverhaeltnis), "
   "Prüfzahlen 06_kupfer_abstand_min bis 06_kupfer_ring_max")
_N("labor_vorhersage", Ü, "kipp_schleife",
   "Sattel-Knoten-Normalform 1 − q = 2(3a)^{1/4}δ^{1/2} im Kern; Sattelpunkt Sonne–Erde im Kern (sattelpunkt, a_N_aus_kreisbahn), "
   "Prüfzahlen 03_sattel_abstand, 03_sattel_feld_galaxie")
_N("lambda_bruecke", I, "Λ ohne ν festlegen (lambda_bruecke)", "trägt nicht; Λ = ν_kosm·ΣΔs_b bleibt")
_N("lambda_und_rueckstellung", N, None, "Lesart der Dunklen Energie vor dem Fluchtmaß (Bausteine aus "
   "rueckstellung_zeitform); ohne Ergebnis im Modell")
_N("licht_trichter", Ü, "delta", "Takt, δ = 1 − dτ/dt, r₀ und r* im Kern (takt, r0, r_stern)")
_N("linse_ohne_freiheit", Ü, "haufenstoss_trennzeit",
   "Trennzeit beim Haufenstoß t = I_F/ω_orb im Kern, Prüfzahlen 03_haufenstoss_unten, 03_haufenstoss_oben; I_F selbst "
   "im Kern (uebergabe_faltung), Prüfzahlen 03_I_F_keim_klein, 03_I_F_keim_gross, 03_I_F_flach")
_N("linsen", Ü, "ablenkung_linse", "α(b) = (12Δs/b)·I(b/r*) im Kern, Prüfzahl 08_linse_I1; Linsen-RAR nicht im Kern; "
   "Haufenrest im Kern (haufen_rest, linse_punktmasse, 03_haufen_rest_100kpc, 03_haufen_rest_1mpc)")
_N("linsen_isoliert", I, "Linsen um isolierte Galaxien (KiDS, Brouwer 2021; Mistele 2024)",
   "Linse mit Restmode als offene Prüfzahlen 03_offen_60_1, 03_offen_60_2 (nur Forschung)")
_N("massenverhaeltnis", Ü, "02_im_tau_G", "Im τ aus dem Nome (im_tau_aus_nome) für α_G und α (02_im_tau_alpha); "
   "Zufallsprobe für m_p/m_e ohne Ergebnis im Modell")
_N("maximum_schwelle", I, _VERGESSEN, "Kandidaten für ein Maximum von a₀(z) bei z ≈ 0,4; nicht im Kern")
_N("mikro_nome", I, "Mikroskala und Nome",
   "Nome-Schreibweise im Kern (im_tau_aus_nome); Triade der Teilchenebene und μ_p = 3μ_N nicht im Kern")
_N("mikroskala", Ü, "nome_qcd", "Lauf von α_s (alpha_s_lauf, skala_bei_kopplung) und Nome-Skala nome_qcd; Prüfzahlen "
   "02_alpha_s_planck, 02_nome_qcd_planck, 02_empfindlichkeit_plus, 02_empfindlichkeit_minus")
_N("mischregel_schale", I, "Aufteilung an der w-Schale aus einem Rücklauf mit zwei Verlusten",
   "Rechnung, noch nicht im Kern")
_N("mock_schatten", Ü, "Q_tilde_gammareihe", "Q̃(s) = −4Σ Γ(4, ns)/n⁴ als Gegenprobe von Q_tilde; Maaß-Lesart ohne "
   "Ergebnis im Modell")
_N("modell_h0", I, "H₀ im eigenen Atmen und Schallhorizont", "trägt nicht; Zahlen nur Forschung, nicht im Modell")
_N("modulraum_zahlen", I, "Modulraum-Vorschläge", "trägt nicht")
_N("naht_rq_umkehr", I, "Übergang an der w-Schale",
   "R_Q im Kern (r_Q); Rest jenseits R_Q im Kern (rest_jenseits), Prüfzahlen 03_rest_jenseits_rq_wert, 03_rest_jenseits_rq_fehler; Mischform nicht im Kern")
_N("netz_und_bilanz", N, None, "ohne Kern und ohne Ergebnis im Modell (Netzwerk-Summe, Energiebilanz des Autors "
   "durchgerechnet)")
_N("ordnung_pruefstand", N, None, "Exponenten der alten Prüfzahlen gegen Grundgrößen; ohne Ergebnis im Modell")
_N("ort_oder_system", I, "Ort oder System", "trägt nicht")
_N("orthogonalitaet_wurzeln", Ü, "epizyklen_summe", "(κ² + ω_z²)/Ω² = 2 + χ im Kern; übrige Muster ohne Ergebnis "
   "im Modell")
_N("parameter_trio", N, None, "Vorstufe (κ₀, λ, r₀ der Kopfgleichung); ersetzt durch kern_zeit")
_N("q_a0_kopplung", Ü, "zellgrenze", "(1 + Q)·a₀ ≤ K·f³ im Kern, Prüfzahl 03_zellgrenze_wuerfel; Formfaktor f³ je "
   "Zellform im Kern (formfaktor_zelle), Prüfzahlen 03_formfaktor_sechseckig bis 03_a0_max_flächenzentriert")
_N("q_form", Ü, "torus_streifen_nu", "Streifen eines Radiants als Forschungsformel; Q als Formverhältnis trägt nicht "
   "als Wert")
_N("radiale_beschraenkung", Ü, "dynamikmasse_ort", "Dynamikmasse am Sonnenort im Kern (dynamikmasse_ort, Prüfzahlen "
   "07_dynamikmasse_kern, 07_dynamikmasse_ort_schicht, 07_dynamikmasse_ort_ebene); Lesarten des Verlusts verglichen "
   "(a₀ = 1,2), entschieden für die örtliche Regel (a_regel), Außenfeld ersetzt durch nu_aussenfeld")
_N("ramanujan", Ü, "eisenstein_G", "Eisenstein-Reihen und Normierungen im Kern (bose_leiter_integral), Prüfzahlen "
   "05_normierung_E4, 05_zahl_1728; Δ(τ) als q-Reihe im Kern (delta_reihe, delta_modul), Prüfzahlen 05_delta_koeffizient_3, 07_delta_gegen_reihe")
_N("ramanujan_unbekannt", N, None, "ohne Ergebnis im Modell; Q̃ als Gammareihe ersetzt durch Q_tilde_gammareihe")
_N("ramanujan_weitere", Ü, "drehwinkel_zelle", "Drehordnungen aus der Bilderzahl der Zellgitter im Kern (drehwinkel_zelle), "
   "Prüfzahlen 03_drehung_moeglich_60 bis 03_drehung_sechseck_120")
_N("regelkreise", I, "Regelkreise je Ebene", "Schleifensumme und Kac-Budget im Kern (schleifensumme, kac_budget); "
   "Verstärkung je Ebene nicht im Kern")
_N("rest_mode_isolierte_galaxie", N, None, "alte Maxwell-Form von w; ersetzt durch gewicht_w")
_N("ringe_oertlich", Ü, "praezession_ring", "Präzession geneigter Ringe (knotenpraezession_klein, vertikalfrequenz), "
   "Prüfzahlen 07_ring_12_2 bis 07_ring_30_89")
_N("rotation", Ü, "horizont_kerr", "Kerr in Faltungsgrößen (kerr_omega_phi, kerr_omega_r, kerr_omega_theta, "
   "nullflaeche_takt) im Kern; Knotenpräzession geneigter Ringe im Kern (knotenpraezession_scheibe, knotenpraezession_frequenzen, 07_ringe_faktor_innen); Siegel-Periodenmatrix ohne Ergebnis im Modell")
_N("rueckkehr_2pi", N, None, "Vorstufe von kac_schliessung (ν = 2π); Zählung je Radiant als Prüfzahl 03a_Q_radiant")
_N("rueckkehr_horizont", N, None, "Vorstufe von kac_schliessung (Ω_b/Ω_m aus Kac); ersetzt durch kac_zeit_schluss")
_N("rueckkehr_zustandsraum", Ü, "kac_budget", "Kac-Budget Σ 1/⟨T_i⟩ = 1 im Kern; Fokker-Planck-Lesart ohne Ergebnis "
   "im Modell")
_N("rueckstellung_zeitform", N, None, "Lesart a₀ folgt einem Horizont und Maxwell-Gewicht aufgegeben; Bilanz ersetzt "
   "durch rueckstellbilanz")
_N("rueckweg_faltung", N, None, "Vorstufe von gleichung_ganzes (q als Green-Funktion von D∂² − Γ); ersetzt durch "
   "erneuerung")
_N("schale_zeitlich", Ü, "uebergabe_vor_umkehr", "Übergabe vor der Umkehr im Kern (uebergabe_vor_umkehr, "
   "rand_schale_aussen, gewicht_schale_aussen), Prüfzahlen 03_I_vor_umkehr, 03_schale_rand, 03_w_aussen_133; Linse der "
   "Schale je Lesart als offene Prüfzahl 03_offen_27_1 (Bahnrechnung fehlt)")
_N("schalter_chi", Ü, "chi", "χ = s/(e^s − 1) im Kern (Besetzung einer Mode, Gegenprobe chi_nome)")
_N("scheibe_feldform", Ü, "eta_scheibe", "Feldform gegen radiale Regel (0,53 / 0,79 / 0,92) im Kern (feldform_loeser, "
   "eta_scheibe, gleichstand_scheibe, feldform_naeherung), Prüfzahlen 07_feldform_viertel_tief, 07_feldform_rd_tief, "
   "07_feldform_rd_hoch")
_N("schliessungsregel", N, None, "Maxwell-Form von w (Vorschlag A13); ersetzt durch gewicht_w")
_N("schluss_oben", Ü, "q_fenster", "Q ≥ 5 aus r*(M_b,Zelle) ≤ Λ im Kern, Prüfzahl 03_rstern_zu_lambda")
_N("schwelle_anpassung", N, None, "Lesung von a₀ unter der Feldform bei a₀ = cH₀/2π; Lesart aufgegeben, ohne "
   "Ergebnis im Modell")
_N("schwelle_je_galaxie", I, _SCHWELLE,
   "Anstieg der Schwelle mit Σ₀ im Kern (a0_je_drittel), Prüfzahlen 03_sparc_a0_drittel_anstieg, "
   "03_sparc_a0_drittel_fehler; Lesung je Galaxie im Kern (a0_je_galaxie) ohne Profil der Störgrößen")
_N("schwelle_modell_mess", N, None, "Vorstufe von a0_lesungen (cH₀/2π gegen 1,2·10⁻¹⁰, Rolle von Υ)")
_N("schwelle_rueckkehrzeit", Ü, "kac_zeit_schluss", "⟨T_Ort⟩ = 2Λ/c für jede Schrittlänge im Kern; Lesart K2 ohne "
   "Ergebnis im Modell")
_N("setzungen_gleichung", Ü, "verlust_gamma",
   "Verlust je Schritt = Taktgefälle, Γ = a_N/c im Kern (Gegenprobe verlust_gamma_taktgefaelle)")
_N("skalen_gross", Ü, "03a_ds_milchstrasse", "Δs je Ebene (ds), Prüfzahlen 03a_ds_milchstrasse, 03a_ds_haufen, "
   "03a_ds_universum")
_N("skalen_klein", Ü, "ds", "Grundformeln A1–A3 im Kern (neigung_xy, takt), Prüfzahlen 02_ds_erde, 02_ds_proton")
_N("sparc_aussenfeld_rand", Ü, "03_aussen_zahl", "Auswahl der Außenpunkte als Prüfzahl; Reststeigung mit Bootstrap "
   "im Kern (reststeigung_aussen), Prüfzahlen 03_reststeigung_wert, 03_reststeigung_fehler")
_N("sparc_boseform", Ü, "a0_an_punkten", "Bose-Form an den SPARC-Punkten (fq_form, nu_ds) im Kern")
_N("sparc_flaechendichte", Ü, "toomre_sigma_dyn",
   "Toomre-Integral und ideale Scheibe (sigma_dyn_expscheibe) im Kern, Prüfzahlen 03_expscheibe_hoch, "
   "03_expscheibe_tief; Rest je Galaxie im Kern (sigma_dyn_rest), Prüfzahlen 03_sigma_dyn_rest_mittel bis 03_sigma_dyn_rest_anstieg_fehler")
_N("sparc_spektraldimension", Ü, "ds_anpassung", "d_s an SPARC als Prüfzahl 07_ds_sparc; Fehler über Bootstrap "
   "im Kern (ds_streuung, boden_anpassung, boden_grenze), Prüfzahlen 07_ds_fehler_ziehung, 07_boden_grenze")
_N("sparc_streuung", Ü, "streuung_messfehler", "Monte Carlo über die Messfehler und a₀ je Galaxie (a0_je_galaxie), "
   "Prüfzahlen 03a_streuung_rest, 03a_a0_je_galaxie")
_N("sparc_test", N, None, "Vorstufe von sparc_universal und sparc_streuung; ersetzt durch rest_galaxie")
_N("sparc_universal", Ü, "rest_galaxie", "Rest je Galaxie mit Gasdruck und Außenfeld, Prüfzahlen 03a_sparc_galaxien, "
   "03a_innen_aussen, 03a_rang_sb_gasreich")
_N("stern_trichter", N, None, "Vorstufe von carter_herleitung; ersetzt durch carter_faktor")
_N("stoerungen_rekombination", Ü, "03a_rest_w_rekombination",
   "w entscheidet bei z = 1090 (gewicht_hintergrund); a_N einer linearen Störung im Kern (neigung_stoerung), Prüfzahl 03_stoerung_aN_rekombination")
_N("straffung_rueckstellung", N, None, "Lesart der Dunklen Energie vor dem Fluchtmaß (Bausteine aus "
   "rueckstellung_zeitform); ohne Ergebnis im Modell")
_N("strittig_zweiter_weg", I, "Strittige Punkte, zweiter Weg", "zweiter Weg zu hopf_knoten, info_geometrie, "
   "modulraum_zahlen, kontakt_atmen; nicht im Kern")
_N("stufe_zwei", Ü, "bose_leiter_integral", "Stufe d = 2 mit Wert 2ζ(3) im Kern; Zusatzpotential (zusatzpotential)")
_N("syntropie", Ü, "grenzmasse_hawking",
   "T_H = T ⇒ M_T im Kern, Prüfzahl 04_grenzmasse_T_CMB; Entropieexport, Page-Zeit und Füllung des Horizonts "
   "im Kern (04_entropieexport_erde, 04_page_zeit, 04_fuellung_horizont)")
_N("tau_ebenen", Ü, "spiegel_radius", "Leiter symmetrisch um r* (r ↦ r*²/r) im Kern; Regel nicht spiegelsymmetrisch")
_N("tau_realteil", N, None, "Vorstufe (Re τ = κ/Ω); ersetzt durch realteil_quadrat_punkt")
_N("temperatur_tiefe_ebenen", Ü, "masse_aus_temperatur",
   "M_b = (S_T·kT/μm_p)²/(G·a₀) im Kern; f = 1/S_T je Ebene als offene Prüfzahlen 03_offen_31_1, 03_offen_32_1; "
   "Haufentemperaturen aus dem Modell genommen (Gedächtniswerte)")
_N("tensor", Ü, "Q_tilde", "Kinetikfunktion Q̃ der kovarianten Form im Kern (Gegenprobe Q_tilde_gammareihe)")
_N("theta_bezug", Ü, "gewicht_hintergrund", "θ̃ absolut gegen das feste Λ, 1 − w im Hintergrund, Prüfzahl "
   "03a_rest_w_heute")
_N("topologie_ganzes", Ü, "dirac_verhaeltnis", "Dirac 1/(2α) = 68,5 im Kern (dirac_verhaeltnis, 02_dirac_verhaeltnis); Spinzahl im Kern "
   "(spinzahl, 09_spinzahl_elektron)")
_N("torus_dynamik", Ü, "rastbreite", "Rastbedingung |Δω| ≤ 2k|cos φ| im Kern; Raten 06_rate_quer, 06_rate_rastpunkt "
   "im Kern (rate_quer, rate_rastpunkt)")
_N("torus_wechselkraefte", Ü, "06_einrasten_min",
   "Kuramoto-Schwelle (kuramoto_einrasten) 06_einrasten_min, 06_einrasten_max, Draht (draht_schleife) "
   "09_draht_kanaele_min")
_N("trenntests_b_dm", I, "Trenntests Bindung gegen Dunkle Materie",
   "η und Umkehrskala nicht im Kern; Skalen als offene Prüfzahlen 03_offen_61_1 bis 03_offen_61_3; Rest je Gasanteil im Kern (rest_gasanteil, 03_rest_gasanteil_unten)")
_N("triade_ebenen", Ü, "realteil_quadrat_punkt_triade",
   "Re τ_M(r*) und kT_I als Handwerte der Gegenproben (temperatur_innen_handwert)")
_N("trichter_invarianten", Ü, "schnittmasse_glaettung",
   "Gezeit am Sonnenort und M× = a₀³/(G𝒯²) im Kern (gezeit_sonnenort), Prüfzahlen 07_schnittmasse, 07_ell_erde")
_N("trichterkette_gewicht", Ü, "dichte_umkehr_heute", "1 + Δ_ta aus dem Top-Hat im Kern; Lokale Gruppe aus dem "
   "Modell genommen (Gedächtniswerte)")
_N("trichterkette_sonne", I,
   "Trichterkette der Sonne (Sgr A*, Milchstraße, Sgr-Zwerg, LMC, M31, Virgo, Großer Attraktor, Shapley)",
   "Rechnung, noch nicht im Modell; Eingaben aus dem Gedächtnis")
_N("tunnel_aussenfeld", Ü, "nu_aussenfeld", "Zusammensetzung senkrecht s⁴ = s_i⁴ + s_e⁴ im Kern, Einseitigkeit m = 1 "
   "(einseitigkeit_m1)")
_N("uebergabe_faltung", Ü, "uebergabe_faltung", "I_F im Kern (uebergabe_faltung), Prüfzahlen 03_I_F_keim_klein, "
   "03_I_F_keim_gross, 03_I_F_flach")
_N("umlaufsinn", Ü, "07_coth_schwelle", "ein Umlaufsinn, beide Sinne als Prüfzahlen 07_coth_schwelle, "
   "07_beide_sinne_faktor; gegenläufige Scheiben als Prüfzahl 03_gegenlaeufig_drehsinn")
_N("vakuum_richtungen", I, "Vakuumteil je Richtung", "trägt nicht; Scherung als offene Prüfzahl 01_offen_01_1 "
   "(nur Forschung)")
_N("wachstum_w_qs", I, "Gewicht vor Q_s im linearen Wachstum",
   "fσ₈ im Kern (wachstum_uebergabe), Prüfzahlen 03_fsigma8_uebergabe, 03_fsigma8_sofort; fσ₈ zwischen Strukturen als offene Prüfzahl 07_offen_36_1 (σ(R) aus einem Spektrum fehlt im Kern)")
_N("wahl_im_tau", Ü, "07_im_tau_erdbahn", "Wahl |q| = Δs/r mit Im τ aus dem Nome (im_tau_aus_nome), Prüfzahlen "
   "07_im_tau_erdbahn, 07_im_tau_milchstrasse")
_N("warp_amplitude", Ü, "warp_viertelzeit", "t_W = (π/2)/ω_p im Kern, Prüfzahl 03_warp_viertel")
_N("warp_drehmoment", Ü, "warp_aufwicklung", "Aufwicklung des Warps im Kern (milchstrasse_quelle, frequenzen_loesung, "
   "saeule_loesung, warp_aufwicklung), Prüfzahlen 03_warp_aufwicklung_regel, 03_warp_aufwicklung_newton; Drehmomente "
   "auf den geneigten Ring als offene Prüfzahl 03_offen_78_1 (lineare Antwort der Feldform fehlt)")
_N("warp_praezession", N, None, "Vorstufe von warp_drehmoment und warp_amplitude (a₀ = 1,2, QUMOND-Feldform)")
_N("wechselwirkungen", Ü, "schleife_thermisch", "q = (∂P/∂T)/G am thermoelektrischen Element im Kern; selbsterregtes "
   "Minimalsystem im Kern (schleife_selbsterregt, 05_q_selbsterregt)")
_N("wep_glaettung", I, "Gleichheit des Fallens und Glättungsschwelle mit der örtlichen Regel",
   "Δa = 0 am Schwerpunkt als Beziehung im Modell; Hill-Radius im Scheibenmodell als Prüfzahlen 07_hill_sonne_newton, 07_hill_sonne_regel; Mittel über die Glättungszone nicht im Kern")
_N("wicklung_rueckkehr", Ü, "kac_rueckkehrzeit", "Kac für die Kreisdrehung im Kern; Irrweg als Handwert der "
   "Gegenprobe rueckkehr_geschlossen_irrweg")
_N("wiederkehrende_muster", Ü, "bose_leiter", "s^d·n_B und Wien-Stellen je Dimension im Kern (bose_stufe_max)")
_N("wirkung_luecken", I, "Lücken der Wirkung", "Gegenspannung als offene Prüfzahl 07_offen_28_1; Wirbelanteil nicht "
   "im Kern")
_N("wirkung_oertliche_regel", N, None, "Vorstufe von wirkung_luecken (Wirkung der örtlichen Regel, 1/ν als Faktor)")
_N("zeit_rueckkehr", Ü, "tau_D", "τ_D = c/a₀ = ϖ_R²/D im Kern (Gegenprobe tau_D_rueckkehrtiefe)")
_N("zeitgesetz_a0", Ü, "omega_vakuum", "Ω_s = (π²/45)·Θ² und Θ (theta_) im Kern; Zeitgesetz a₀ ∝ H als Lesart "
   "aufgegeben")
_N("zelle_form", Ü, "formfaktor_zelle", "Formfaktoren im Kern (formfaktor_zelle, zellgrenze, drehwinkel_zelle, "
   "streuung_kosmisch, anteil_isotrop); Winkelmoden und d_LSS im eigenen Atmen als Forschungsformeln "
   "(winkelmoden_zelle, lauf_eigenes_atmen), Prüfzahlen 03_offen_38_1 bis 03_offen_42_2 nur Forschung")
_N("zelle_frueher_kosmos", N, None, "Vorstufe (Lesarten der Zelle); Lesart a₀ fest übernommen, ersetzt durch "
   "a0_frueh_relativ")
_N("zustandsraum", Ü, "05_delta_koeffizient_3", "q-Entwicklung mit τ(n) als offene Prüfzahl (Formel fehlt); "
   "Schleifensumme im Kern (schleifensumme)")
_N("zwerge_aussenfeld", Ü, "sigma_plummer_jeans", "Jeans mit gerichtetem äußerem Feld im Kern (sigma_plummer_jeans, Lage "
   "kugel/gleich/senkrecht), Prüfzahlen 03_druck_df2_gleichgerichtet, 03_druck_df2_ohne_feld")
_N("zwerge_tiefengedaechtnis", Ü, "02_helium_strahlend",
   "Decke kT_max für Weiße Zwerge (temperatur_max), Prüfzahlen 02_helium_zu_grenze, 02_helium_konvektiv")
_N("zwerge_uebergabe", I, "Zwerge vor der Übergabe", "trägt nicht")

# ------------------------------------------------------------------ verwaiste Ausgaben und Übersetzungen
# Namen aus Stamm und Endung (katalog_eintrag.AUSGABE); gelesen wird keine Datei.
for _s, _a in (("aussenfeld_sparc", "a012"), ("glaettung_stetig", "a011"), ("gravitationswellen", "a011"),
               ("radiale_beschraenkung", "a011"), ("rest_mode_isolierte_galaxie", "a012"), ("schliessungsregel", "a011"),
               ("warp_drehmoment", "a011"), ("warp_praezession", "a011"), ("zwerge_aussenfeld", "a011")):
    _N(f"{_s}_{_a}.{AUSGABE[0]}", N, None,
       f"verwaiste Ausgabe: Lauf von {_s} mit a₀ aus der Umgebung ({_a}); Einordnung beim Skript")
_N(f"ordnung_err.{AUSGABE[0]}", N, None, "verwaiste Ausgabe: Fehlerausgabe eines Laufs von ordnung_pruefstand")
_N(f"ordnung_pruefstand_E.{AUSGABE[1]}", N, None, "verwaiste Ausgabe: Exponentenmatrix von ordnung_pruefstand")
_N("__pycache__/faltungsknoten.cpython-313.pyc", N, None, "übersetztes Modul ohne Quelle; kein Skript faltungsknoten")
