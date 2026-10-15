# -*- coding: utf-8 -*-
"""
Themen der Forschung (Klassen und Regeln in thema.py). So füllen: THEMEN.neu(kennung, Frage, Invariante(…),
[Versuch(…), …], {gesperrter Weg: Grund}, Entscheider, Stand). Jede Forschungszeile mit "Rechnung" gehört als Versuch
zu einem Thema (Wächter, zählend); die Sammelzeile „Vergessene Funde …“ ist je Fund dem passenden Thema zugeordnet.
Datum = erster Commit des Skripts (git log); Versuche ohne Forschungszeile sind Skripte, die nie eingetragen wurden.
haengt_an nennt nur, worauf das Urteil selbst ruht; ein strukturelles „trägt nicht“ bleibt ohne Eintrag.
"""
from .thema import THEMEN, Invariante, Versuch

REGISTERMODUL = True

# Kurztitel der Forschungszeilen, die mehrfach vorkommen
VERGESSEN = "Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell"
PLANCK = ("OMEGA_B_H2", "OMEGA_C_H2")                     # Q, Ω_m, θ* und das Ziel 0,442 hängen daran
ZWERGE = ("DF2_SIGMA", "DF4_SIGMA", "FORNAX_SIGMA", "CRATER_II_SIGMA", "DRACO_SIGMA", "PAL14_SIGMA", "NGC2419_SIGMA")


def V(kennung, datum, weg, ergebnis, zeile=None, haengt_an=(), skripte=None):
    """Versuch; skripte = (kennung,) wenn nicht anders genannt."""
    return Versuch(kennung, datum, weg, ergebnis, haengt_an, zeile, (kennung,) if skripte is None else skripte)


# ------------------------------------------------------------------ Dunkle Energie
THEMEN.neu(
    "dunkle_energie",
    "Was ist der zeitfeste Teil der Dunklen Energie neben dem Vakuumteil Ω_s (k_DE/k_s ≈ 2,9–3,1) und woher kommt seine "
    "Form w(z)? Hinweis des Autors: Tafel 5 und 6 („DM gegen B“) meinen Dunkle Energie; Invariante der Frage: Dunkle Energie "
    "als Bindungs- bzw. Schließungsstruktur statt Quelle; Trenntests je Kanal.",
    Invariante("k_DE = ρ_DE·Gc²/a₀² (zeitfest, wie Q ein Verhältnis ohne Uhr); Vakuumteil k_s = π³/30 = 1,034. Größen "
               "mit H₀ oder t₀ gelten nur in einer Epoche", "GI", "Hintergrund"),
    [V("gegenlauf_dunkle_energie", "2026-09-29", "Gegenlauf: Rückkehr aus dem mitatmenden Anteil 1 − F, Rückstellung "
       "aus dem gebundenen Anteil F (Press-Schechter, Sheth-Tormen); nie in den Forschungsstand eingetragen",
       "offen", haengt_an=("frei:Mindestmasse", "DE_SIGNIFIKANZ")),
     V("lambda_und_rueckstellung", "2026-09-29", "Λ als ganze innere Faltung der Zelle, Rückstellung als Bindungsrate "
       "dF/d ln a; w₀ −0,64 bis −0,80 je nach Mindestmasse 10⁸–10¹³ M☉; nie eingetragen", "offen",
       haengt_an=("frei:Mindestmasse", "DE_SIGNIFIKANZ")),
     V("straffung_rueckstellung", "2026-09-29", "Rückstellung als Straffung ΣΔs_i/R_i gebundener Trichter, als Rate; "
       "nie eingetragen", "offen", haengt_an=("frei:Gewichtung", "DE_SIGNIFIKANZ")),
     V("atmen_w", "2026-09-29", "Strukturterm im Vakuumwert: Ω_s = π²/45 = 0,219 bei a₀ = cH₀/2π, w = −1; Syntropie-"
       "Spirale als Rest", "Teil im Modell"),
     V("horizont_kreislauf", "2026-09-29", "Horizontentropie in H-Form gegen Λ-Form; Mindestrate neuer Faltung "
       "Λ̇/Λ ≥ Ṡ_i/(2S_Λ) = 3,83·10⁻³⁷ s⁻¹; S_H/S_Λ = 0,244; nie eingetragen", "im Modell"),
     V("de_ratenschwelle", "2026-09-29", "Rückstellung über die Ratenschwelle w = ½ (fällt auf 1 % mit der Umkehr "
       "zusammen); gebundener Massenanteil F(z)", "trägt nicht", "Dunkle Energie über die Ratenschwelle",
       ("frei:Glättungsmasse",)),
     V("de_phasenseite", "2026-09-29", "Phase der Rückläufe in der inneren Richtung", "trägt nicht",
       "Dunkle Energie über die Phasenseite des Ganzen"),
     V("kontraktion_atemrate", "2026-09-29", "Kontraktion oder Atemrate: keine Umkehr, H_s = 33,3 km/s/Mpc, Mindestrate "
       "3,8·10⁻³⁷ s⁻¹; Faltungsverlust der Strahlung noch nicht im Modell", "Teil im Modell", "Kontraktion oder Atemrate"),
     V("engstelle_rate", "2026-09-29", "Einfall und Abstrahlung an Schwarzen Löchern als Quelle von Λ: Vorzeichen falsch "
       "(wachsendes Λ senkt Ω_s ∝ Λ⁻²); ρ_BH und S_obs [G] setzen nur den Betrag", "trägt nicht",
       "Engstelle und Schwarze Löcher"),
     V("engstelle_global", "2026-09-29", "Horizont als Engstelle; Kopplungen Λ ∝ Δs_H, dΛ = ν·dΔs_H, Vakuumteil als "
       "Quelle tragen nicht (w_a stets falsches Vorzeichen)", "Teil im Modell", "Globale Engstelle"),
     V("de_uebergabe_anteil", "2026-09-29", "Rest als übergebene Mode: zählt im Ganzen als Staub weiter; fiele sie weg, "
       "wüchse a₀ seit z = 10 um 1,1–3,2; Phantom-Vorzeichen", "trägt nicht", "DE-Rest als Übergabe der Mode"),
     V("de_bilanz_glieder", "2026-09-30", "Bilanz aller Glieder der Ausdehnung: kein positives Glied mit fester Dichte "
       "neben dem Vakuumteil (⟨ℒ_𝒮⟩ ≤ 1,7·10⁻⁵·Ω_s, Scherung, B² ≤ 10⁻⁹)", "trägt nicht", "Bilanz aller Glieder"),
     V("vakuum_richtungen", "2026-09-30", "Vakuumteil je Raumrichtung (3Ω_s): die 3 steckt schon in π⁴/15, zählte "
       "doppelt; Rückstellung ist die eine Volumenmode (b₃ = 1)", "trägt nicht", "Vakuumteil je Richtung"),
     V("de_kehrseite_bindung", "2026-09-30", "Überschuss der gezählten Faltung über 1 + Q je gebundenem Bereich; in der "
       "vom Modell verlangten Lesart des Ganzen ≤ 0,004 mit falschem Vorzeichen (die Summe 1,1–10,7 hing an der freien "
       "Mindestmasse)", "trägt nicht", "DE als Kehrseite der Bindung"),
     V("de_fluchtmoden", "2026-09-30", "inaktives Fluchtmaß je Bereich (Liouville je Baryon): Vorzeichen richtig, "
       "Ω_X = 0,18–0,33, Gipfel bei z ≈ 2–6", "trägt nicht", "DE als verdrängte Fluchtmoden", PLANCK + ("H0_KMS",)),
     V("de_fluchtmoden_hierarchie", "2026-09-30", "Fluchtmaß geschachtelt gezählt (Verschmelzungsbaum): w₀ ≈ 0, H bei "
       "z = 1–3 16–47 % zu hoch", "trägt nicht", "Fluchtmaß hierarchisch", PLANCK + ("H0_KMS",)),
     V("de_fluchtmoden_verzoegert", "2026-09-30", "Fluchtmaß verzögert über die erste Rückkehr (τ_D = 86 Gyr): Zeitpunkt "
       "richtig, Menge halb (Ω_X = 0,22 von 0,442)", "trägt nicht", "Fluchtmaß verzögert über die Rückkehr",
       PLANCK + ("H0_KMS",)),
     V("de_fluchtmoden_erneuerung", "2026-09-30", "Fluchtmaß mit Erneuerung: erste Rückkehr ist obere Schranke, "
       "Ω_X(0) = 0,223", "trägt nicht", "Fluchtmaß mit Erneuerung", PLANCK + ("H0_KMS",)),
     V("de_fluchtmoden_mass", "2026-09-30", "Maß erzwungen (ganze Tori innerhalb R_ta, ln B − ln A = 1,392): "
       "Ω_X(0) = 0,32, w₀ = −0,45, w_a = −4,1", "noch nicht im Modell", "Fluchtmaß, Definition des Maßes"),
     V("de_kanal_atmen", "2026-09-30", "Kanal ins Atmen: Maß wirkt auf die Schleifen, nicht auf Λ; kein zwingender "
       "Faktor 1,33 (6D-Monte-Carlo 1,3735)", "im Modell", "Kanal ins Atmen"),
     V("de_rest_anstieg", "2026-09-30", "Fluchtmaß selbstkonsistent (Atmen → Strukturbildung → Maß): Ω_X 0,328, zusammen "
       "83 % von 0,687 [G]; keine weitere zwingende Quelle", "offen", "Fluchtmaß selbstkonsistent"),
     V("de_laufzeit", "2026-09-30", "Laufzeit des Übergangs: K1 (Erstdurchgang bis ϖ_R) erzwungen, kein freier Hebel",
       "offen", "Laufzeit des Übergangs"),
     V("modell_h0", "2026-09-30", "geschlossenes H₀ im eigenen Atmen (Fluchtmaß K1): 63,45 km/s/Mpc, 100θ* 1,0198 gegen "
       "1,0411 [G], BAO χ² gegen das Standardmodell", "trägt nicht", "H₀ im eigenen Atmen und Schallhorizont",
       PLANCK + ("DE_SIGNIFIKANZ",)),
     V("de_invariante", "2026-09-30", "zeitfeste Größe k_DE = ρ_DE·Gc²/a₀² aus θ* [G]: 3,18 / 2,95 bei a₀ = 1,06 / 1,10; "
       "jede Lesart verlangt a₀ auf 0,2 %; Form fest wie der Vakuumteil", "noch nicht im Modell",
       "Dunkle Energie als Invariante der Grundwerte", skripte=("de_invariante", "modell_h0")),
     V("regelkreise", "2026-09-30", "Schleifenverstärkung je Ebene; Teilkreise derselben Rückkehr geben höchstens den "
       "Vakuumteil, nötig wäre eine zweite Rückkehr mit eigenem Budget", "noch nicht im Modell", "Regelkreise je Ebene"),
     V("de_zyklen", "2026-09-30", "Vakuumteil je Zyklus der Zelle (n = 1…8): als Herleitung nicht (x³ ist Flussgesetz "
       "einer Masse, b₁ kommt nicht vor); a₀ aus θ* je n; SPARC gibt n = 2,83 ± 1,19, nur n = 1 ausgeschlossen",
       "trägt nicht", "Dunkle Energie je Zyklus der Zelle"),
     V("maximum_schwelle", "2026-09-29", "Lage des Maximums der Dunklen Energie (Notiz 56): hängt an einer Konvention",
       "offen", VERGESSEN),
     V("lauf_eigenes_atmen", "2026-09-30", "Lauf des eigenen Atmens Ẋ_an(t) in forschung/kosmos.py (lauf_eigenes_atmen: "
       "Bildung im Vergleichshintergrund V0, selbstkonsistent S1, S2; Anfangsspektrum, Kugelumkehr, Fluchtmaß je Baryon): "
       "X_an/X heute 4,71 % (V0), 5,20 % (S2); d_LSS 2,09 % ferner, 28,37 Gpc, (1 + Q)a₀ ≤ 6,741·10⁻¹⁰·f³, Q ≤ 5,128 "
       "stimmen mit S2; Laufzeit S2 etwa 85 s über der Grenze des Laufs, darum nur Forschung", "Teil im Modell",
       haengt_an=PLANCK + ("A_S", "N_S", "TAU_RE", "M_NU"), skripte=()),
     V("de_adiabatisch", "2026-09-30", "Operator statt Quelle, adiabatische Rückwirkung (J = E/ω je schnellem Zyklus, "
       "Ebene I gegen das langsame Atmen der L_i): Die Form trägt (N fest je Zelle: w = 0; N ∝ Volumen: w = −1; "
       "ω ∝ 1/L: w = 1/3), der Betrag nicht (Kanalzählung 1 statt 2,9–3,1). Geometrischer Magnetismus ≈ 0 (nur "
       "Scherung, ≈ 10⁻¹⁰ rad). Reibung gibt w_a > 0 (gemessen w_a < 0). Lesart zur Prüfung: je Faltung gezählt = Mode ξ "
       "(w = 0), je Raumpunkt = Vakuumteil (w = −1)", "noch nicht im Modell", skripte=()),
     V("de_pql_energie", "2026-09-30", "𝒫_qL entzieht Energie: der symmetrische Teil von 𝒫_qL ist negativ semidefinit "
       "(Kette N = 64, Maximum der Spektralnorm bei q ≈ 0,49, null bei q = 0 und q → 1); die Energie geht ans Gedächtnis, "
       "nicht verloren. Rechnung, noch nicht im Modell", "noch nicht im Modell", skripte=()),
     V("de_fluchtmoden_faser", "2026-09-30", "Verschiebung der Fluchtmoden in die Faser (Betrag der Dunklen Energie): je Raumpunkt "
       "sofort (w = −1) gibt nahezu feste Form, w = −1,03 / −1,10 / −1,18 bei z = 0 / 1 / 2, BAO-χ² 14 (Λ 18), H₀ 68,95, "
       "Betrag an θ* angepasst; je Faltung (w = 0) χ² 669; verzögert χ² 224; Mindestmasse 10⁶–10¹² M☉ gibt k 3,00–3,68 [W]; "
       "offen: Umrechnung vom Fluchtmaß X in Zyklen je Raumpunkt", "Teil im Modell",
       haengt_an=PLANCK + ("H0_KMS", "frei:Mindestmasse"), skripte=()),
     V("de_trenntests_kanal", "2026-09-30", "Trenntests je Kanal (Lensing, Wachstum, Halo, Gedächtnis, Zustand): es fehlt Menge, "
       "keine Form je Kanal; Vorhersage: H_eff ±0,5 % je nach Umgebung", "noch nicht im Modell",
       haengt_an=PLANCK, skripte=()),
     V("de_torus_pythagoras", "2026-09-30", "Torus-Geometrie: Windungswinkel tan²α = κ²/Ω² = 1 + χ(s), cos²α = 1/(2 + χ) "
       "(Kepler 45°, tief 54,74°); kosmische Lesart Ω_DE/Ω_s = 2 + χ(s_kosm) = 2,917 liegt im Band, die Zuordnung "
       "(2Ω² ist in 07 der leere Wert) ist aber nicht begründet, aus 07 folgt höchstens 1 + χ = 1,917; Zufallsprobe "
       "1,7–6,1 % Treffer", "trägt nicht", haengt_an=PLANCK, skripte=()),
     V("de_topologie_ringe", "2026-09-30", "Torus mit Hauptkreis (parametrisch), Querringen (transversal), Windung "
       "(tangential): Takt je Ring ⇔ g_tt·g_rr = −1 exakt im Vakuum; Gauß-Bonnet ∫K = 2π(1 − T), Horizont bei 6Δs "
       "als Wechsel der Topologie; Ablenkung und Perihel je zur Hälfte Takt und Ringe; η − 1 = Wirbel-(Windungs-)anteil "
       "der Regel. Dunkle Energie: Bindung vergrößert Λ, Vakuumteil fällt, w ≥ −1 mit w₀ > −1, w_a < 0, Betrag aus dΛ "
       "heute 10⁻¹⁹ (trägt nicht); der Operator-Anteil (Fluchtmoden je Raumpunkt) bleibt offen. Offen: 11 Spannungen "
       "zwischen Modelltext und Bild (Rückkehrrichtung gegen Querring, Halbgerade gegen Ring, Zeit als δ gegen Phase, "
       "zwei Tori T³ und Bahn)", "noch nicht im Modell", skripte=()),
     V("de_nz_kac", "2026-09-30", "Betrag aus der Kac-Erneuerung: Fluchtanteil je Zyklus 1 − q = 1/ν_kosm, n_Z = X/ln ν = "
       "4,576/1,851 = 2,47, ρ_DE = ρ_s·[1 + n_Z·F(a)/F(1)], ohne Anpassung. a₀ = 1,1: k_DE/k_s = 3,47, 100θ* 1,0511, "
       "BAO-χ² 58; θ* verlangt a₀ = 1,025 (SPARC 1,106 ± 0,23), dort χ² 14, H₀ 69,05, w₀ −1,024. Setzungen: jeder "
       "Zyklus trägt k_s (nicht hergeleitet), X hängt an M_min [W]; andere Umrechnungen (q, X(1−q), X) liegen weit daneben",
       "noch nicht im Modell", haengt_an=PLANCK + ("frei:Mindestmasse",), skripte=()),
     V("de_kanal_s2", "2026-09-30", "Herleitung „ein Zyklus = ein Vakuumkanal“: folgt nicht. Q̃(0) ist die ganze Bose-Summe der "
       "Rückläufe einer inneren Richtung (n = 1 trägt 92,4 %); ein Rücklauf mehr ist schon darin (k = 1, a₀ = 1,91). Ein Kanal "
       "je Zyklus heißt eine eigene geschlossene innere Richtung je Zyklus (1/ℓ_eff² = (1 + n_Z F)/ℓ₀²), die das Modell nicht "
       "hat", "trägt nicht", skripte=()),
     V("de_gauss_bonnet_ganzes", "2026-09-30", "Winkeldefizite aller gebundenen Bereiche auf dem geschlossenen Ganzen: "
       "Ausgleich in den Leerräumen (L) als negative Krümmung, Winkelform der Rückstellbilanz; verhält sich wie Materie "
       "(w = −⅓·d ln F/d ln a), deckt 5–10 % der fehlenden Menge", "trägt nicht", skripte=()),
     V("de_henkel", "2026-09-30", "Bindung als Schnitt und Identifikation (Grundsatz des Autors): Identifikation mit der "
       "inneren Richtung ist eine Abbildung vom Grad m auf den vorhandenen Kreis, kein Henkel (b₁ bleibt); verwertbar: "
       "Wickelzahl n(M) = ln(B/A)/ln ν (Milchstraße 3,05; Mittel je Baryon 2,47)", "trägt nicht",
       skripte=()),
     V("de_buendel_verdrillung", "2026-09-30", "Innere Richtung als verdrilltes Kreisbündel: Σ Chern-Zahlen = 0 auf dem "
       "geschlossenen Ganzen; ganzzahlige Verdrillung lässt den Vakuumteil gleich, halbzahlige senkt ihn (π: −7/8)",
       "trägt nicht", skripte=()),
     V("de_lambda_zaehlung", "2026-09-30", "Λ als topologische Zählung: gebundene Wicklungen stecken schon in ν_kosm "
       "(bei z = 1090 wie heute); Λ fest bestätigt. Λ mit Bindung wachsend: a₀ wandert 4–55 %, w_a > 0; Differenz "
       "alle/mitatmend: Phantom, hängt an M_min", "trägt nicht", skripte=()),
     V("de_einrasten", "2026-09-30", "Einrasten T² → T¹: Maß bleibt (𝒫 doppelt stochastisch), nur die Stütze fällt; die "
       "Regel hebt die Kepler-Verriegelung außen auf, der verlorene Zyklus liegt innen (Anteil 10⁻⁴); das Fluchtmaß ist "
       "Bindung durch Rückkehr, nicht durch Einrasten", "trägt nicht", skripte=()),
     V("de_nichtadditiv", "2026-09-30", "Nicht additiv pro Teilchen (Autor, Tafel 5): globales Maß der Zelle X_glob = "
       "−ln[(1 − F) + ∫(−dF)·B⁻¹] ≈ −ln(1 − F), nur der identifizierte Anteil zählt; heute 2,02 (additiv 4,57), "
       "n_glob = 1,09, k_DE/k_s = 2,09 (Ziel 3,0 verlangt F = 0,975), M_min-abhängig −31 %; θ* verlangt a₀ = 1,32; "
       "Form wie M1", "trägt nicht", haengt_an=PLANCK + ("frei:Mindestmasse",), skripte=()),
     V("de_fundamentalgruppe", "2026-09-30", "Rang von π₁(T³) = 3 als Faktor: die drei Erzeugenden sind Raumschleifen "
       "der Zelle (Länge L), der Vakuumteil hängt an der einen inneren Richtung (Λ). Symmetrisch verteilt bleibt die "
       "Menge gleich; trägt jede Raumschleife einen Vakuumteil derselben Form, gibt das 1 + 3(Λ/L)² = 1,69; genau 3 "
       "nur mit einer eigenen inneren Schließung bei Λ je Erzeugender (zweiter Kanal); ganzzahlig und zeitfest: "
       "w = −1, w_a = 0", "trägt nicht", haengt_an=PLANCK, skripte=()),
     V("de_rollen_zaehlung", "2026-09-30", "Jede Rolle kehrt über dieselbe innere Richtung zurück: Q̃(0) = −4Γ(4)ζ(4) "
       "ist die volle Bose-Summe einer skalaren Mode (Faltung ξ, ein Sinn); im Grenzfall s → 0 ist auch die Faltung "
       "verlustfrei, die Konstante gilt also ebenso für die tangentiale Phase. Mit zwei Komponenten des Felds F: "
       "1 + 2 = 3·Ω_s (0,733 bei a₀ = 1,1; 0,681 bei 1,06), ohne zweite innere Richtung, Bilanz Λ = ν·ΣΔs unberührt. "
       "Offen: welche Felder zählen (Gravitation 1 oder 2, weitere Windungen), und der Teil ist zeitfest (w = −1, "
       "w_a = 0); als reiner Kamm ohne Zustandsdichte nur 1/12", "noch nicht im Modell", haengt_an=PLANCK, skripte=()),
     V("de_rollenfluss", "2026-09-30", "Energie zwischen den Rollen: Arbeit der Strahlung am Atmen Ω_r·z_i wächst "
       "mit der Anfangsepoche ohne Grenze und steckt schon als Faltungsverlust in Λ (w_s = −1,0002); Rollenwechsel "
       "der Bindung als Wärme kT·X je Baryon 10⁵-fach zu klein (𝒫 erhält Maß, nicht Energie). Vorzeichenregel: jeder "
       "monotone Übergang in die innere Richtung gibt w < −1; w₀ > −1 verlangt Rückfluss (Kreislauf, neue Faltung)",
       "trägt nicht", haengt_an=PLANCK, skripte=()),
     V("de_energieformen", "2026-09-30", "Alle Energieformen der Zelle als Massenäquivalent (Strahlung, Sternlicht, "
       "Fusion, Bewegung, Wärme, Magnetfelder): zusammen ≲ 10⁻⁴, 4000-fach unter dem Rest 0,44; als Maß zählt nur X "
       "(Umwandlungen erhalten das Maß); die Zeitformen der Bindungsenergie geben das Maximum bei z ≈ 0,4 nicht",
       "trägt nicht", haengt_an=PLANCK, skripte=()),
     V("de_muster_klein", "2026-09-30", "Selbstähnlich im Kleinen (Atom, Oszillator, Kern, Casimir): Nullpunkt ∝ Zahl "
       "der Moden, Bindung = Rahmen − Umlauf; Verhältnis hängt an der Potenzform (Coulomb 2, Oszillator 1, Kern 1,4), "
       "kein festes 1 : 3; π²/45 = Casimir auf geschlossener Länge = Strahlungsdruck (Koeffizient der Schließung)",
       "trägt nicht", skripte=()),
     V("de_muster_fluss", "2026-09-30", "Selbstähnlich im Großen (Stern: Virial, halbe Bindung geht als Leuchtkraft "
       "∝ Bindungsrate; Galaxie: Q·a_N abgeführt): der Bindungsanteil als Fluss ∝ dF/d ln a gibt w₀ −0,84 bis −0,99 "
       "mit w_a < 0, als Bestand ∝ F das falsche Vorzeichen; Lage des Maximums hängt an σ [W] (Spielzeug-F)",
       "offen", haengt_an=PLANCK, skripte=()),
     V("de_fluss_kern", "2026-09-30", "Bindungsrate dX/d ln a mit Kernzahlen (σ(M), Wachstum, Umkehr; M_min 10⁵–10⁸ M☉): "
       "w₀ −0,71 bis −0,81, w_a < 0, aber Maximum bei z ≈ 5–10; θ* verlangt H₀ ≈ 62, BAO-χ² 55–147 (ΛCDM 21); "
       "deckt sich mit F13 (Maximum wandert mit der kleinsten Masse, z ≈ 0,4 erst oberhalb 10¹³ M☉)",
       "trägt nicht", haengt_an=PLANCK + ("frei:Mindestmasse",), skripte=()),
     V("de_muster_schaukel", "2026-09-30", "Labor (Oberflächenspannung, Membran, Vakuum: w = −1 heißt keine Umlaufenergie, "
       "Energie nur durch den wandernden Parameter; Schaukel J = E/ω): Operatoranteil ∝ H gibt w₀ −0,89 aber w_a +0,39; "
       "w_a < 0 verlangt Anwachsen mit dem Parameter, dann Sättigung bei q = 1 (parametrische Resonanz), Raten fehlen",
       "offen", haengt_an=PLANCK, skripte=()),
     V("de_drei_schleifen", "2026-09-30", "Dehnung je Raumschleife: eine Gerade sieht nur die Projektion "
       "ε_xx = ε_rr·n_x², Mittel 1/3; die drei Schleifen teilen ΔL (Spur = ΔL der Bilanz), keine drei Schließungen "
       "bei Λ (die bräuchten 3·ΣΔs). „1 Schleife fest + 2 tragen die Bindung“ gibt 3·Ω_s mit w₀ −0,92, w_a −1,19, "
       "der Faktor 2 ist aber gesetzt (Ankunftsform liefert 1,31); Bindung ∝ F(a) gibt w < −1", "trägt nicht",
       haengt_an=PLANCK, skripte=()),
     V("de_figur_lb", "2026-09-30", "Zweite Länge der inneren Richtung ℓ_B ≈ 0,64–0,74·ℓ₀ aus einer Figur (Inversion, "
       "Höhensatz, Kappe, Schatten, Spirale): keine Figur erzeugt sie; ℓ_B bleibt Forderung", "trägt nicht", skripte=()),
     V("de_spirale_flaechen", "2026-09-30", "Flächen der Windungen der Rückkehrspirale: Vakuumteil ∝ 1/ℓ² gibt "
       "ρ_k ∝ q^{−2k} (divergent); Σq^{2k} = 1/(1 − q²) = 3,45 verlangt eine zweite innere Richtung (Doppelzählung, "
       "S2); zeitfest, w = −1, BAO-χ² 18 statt 14", "trägt nicht", haengt_an=PLANCK, skripte=()),
     V("de_treideln", "2026-09-30", "Treideln: Wirkung nach Winkel längs cos² / quer sin²; Sehwinkel B trennt die Ebenen "
       "(Sonne bis Haufen rein längs, Kosmos quer 0,078); Ω_s/cos²α = Ω_s(2 + χ) gibt 2,92 bei s_Q (+3,8 %), Grenzwert "
       "s → 0 genau 3, s nicht abgeleitet; ΣΔs = ΔL bleibt nur, wenn das Querstück in I geht: zurückgeholt netto null, "
       "nicht zurückgeholt ein zweiter Kanal in I; fester Winkel zeitfest (w_a = 0)", "trägt nicht",
       haengt_an=PLANCK, skripte=())],
    {"Ω_DE/Ω_s = 2 + χ aus dem Epizykel-Pythagoras": "2Ω² ist der Wert im leeren Raum, nur χΩ² wird gelesen; "
         "Zuordnung nach dem Ziel, Zufallsprobe ohne Beleg",
     "Reibung als Quelle von w_a": "Vorzeichen folgt aus der Struktur (Reibung ∝ H/ρ wächst mit der Zeit): w_a > 0, gemessen w_a < 0",
     "Zählung „drei Schleifen“ (Faktor 3 aus den Kanälen)": "isotropes Atmen: nur die Spur der Drücke zählt, jedes P_i = −ρ; "
         "Zählung 1 statt 2,9–3,1",
     "Summen über gebundene Bereiche (Kehrseite der Bindung, Fluchtmaß je Bereich, Ratenschwelle, Rückstellung als "
     "Bindungsrate)": "nicht zeitfest: folgen der Strukturbildung (Gipfel z ≈ 2–6); zeitfest ist nur k_DE; ≥ 12 Rechnungen",
     "wachsendes Λ als Quelle (Engstelle, Abstrahlung, Übergabe der Mode, Horizont)":
         "Ω_s ∝ Λ⁻² kann nur fallen: w_a mit falschem Vorzeichen",
     "Phase des Ganzen": "Rückläufe in der inneren Richtung tragen keine Phase (H¹ = 0)",
     "Faktoren n = 1, 2, 3, 4, 6, 8 vor dem Vakuumteil suchen":
         "Zahlensuche; als Herleitung trägt keiner (vakuum_richtungen, de_zyklen)",
     "weitere Stellschraube am Fluchtmaß (Hierarchie, Verzögerung, Erneuerung, Laufzeit, Selbstkonsistenz)":
         "Maß und Laufzeit sind erzwungen; keine weitere zwingende Quelle",
     "Urteil über BAO oder θ* gegen das Standardmodell": "Werte [G] und altes Paradigma; nur Hinweis",
    "verzögerte Ankunft des Fluchtmaßes": "BAO-χ² 224, H₀ 94; nur die sofortige Verschiebung trägt die Form",
    "Leiter-Summe (Bose-Leiter d = 1…3) als Betrag": "nur d = 3 hat eine Raumdichte im Vakuum; Summe 1,62 statt 2,9–3,1",
    "„eine Faser“ als Herleitung von k_s": "Zuordnung, kein Betrag; J und Zyklen je Raumpunkt sind nicht einzeln gesetzt"},
    "eine Herleitung einer zweiten Rückkehr mit eigenem Budget oder eines Terms der Zelle Σ|ω_i|²; als Messung a₀ auf "
    "1 % (σ(a₀)/a₀ ≈ 0,01) aus eigenen Daten, dann trennt k_DE/k_s die Lesarten; die Form w(z) ist mit unseren Daten nicht "
    "messbar",
    "offen. Vakuumteil Ω_s im Modell. Vakuum je Richtung und je Zyklus widersprechen sich nicht: beide finden, dass kein "
    "Faktor 3 hergeleitet ist (die 3 steckt in π⁴/15); de_zyklen liest nur die Zahl k_DE/k_s ≈ 3 aus θ* [G], und SPARC "
    "schließt allein n = 1 aus. Offen, hängt an ω_b, ω_c [G]: die Mengenurteile der Fluchtmaß-Kette und modell_h0.")

# ------------------------------------------------------------------ Herkunft von Q
THEMEN.neu(
    "herkunft_q",
    "Woher kommt Q = Ω_c/Ω_b ≈ 5,36 (ν_kosm = 1 + Q), das zeitfeste Gewicht der Mode ξ?",
    Invariante("Q ist ein Formverhältnis ohne Uhr (bei z = 1090 wie heute): q = e^{−2π·Im τ} Nome der Zelle, "
               "Im τ = 0,0272; Schranken ν = 2π (offener Teil ein Radiant) und ν ≥ 6 (Schluss r*(M_b,Zelle) ≤ Λ)", "GI",
               "Wirbel"),
    [V("herkunft_q", "2026-09-29", "q aus derselben Rückkehr wie a₀: Lesarten L1–L4 (geschlossene Halbgerade, Kac-Maß "
       "der atmenden Zelle, R_Q, Rand der Hubble-Kugel); atmende Zelle trägt nicht (ΔN_eff), Q = 2π − 1 nur als Zahlennähe",
       "offen", haengt_an=PLANCK),
     V("camb_q", "2026-09-29", "CAMB-Spektren mit Q = 5,364 gegen 2π − 1: wie genau Messungen die Zahlen trennen",
       "Teil im Modell", haengt_an=PLANCK),
     V("invarianten_ebenen", "2026-09-29", "Invarianten je Ebene L/I/G gegen Übereinstimmungen einer Epoche: S_G/S_I = "
       "Ω_s nur heute, ν_kosm = 6,364, H·τ_D = 5,95 heute; nie eingetragen", "noch nicht im Modell"),
     V("kac_schliessung", "2026-09-29", "Kac-Schließung der inneren Richtung: τ_D ∝ 2Λ/c, θ_c = 1/τ_D, w = 1 − e^{−√x}; "
       "1 − q = 1/2π nur mit der Setzung „ein Radiant“", "im Modell", "Kac-Schließung der inneren Richtung",
       skripte=("kac_schliessung", "gewicht_rueckkehr")),
     V("lambda_bruecke", "2026-09-29", "Λ ohne ν festlegen: drei vorab gewählte Kandidaten geben keine Lösung",
       "trägt nicht", "Λ ohne ν festlegen (lambda_bruecke)", ("frei:Kandidatenwahl",)),
     V("schluss_oben", "2026-09-29", "Schluss oben: Fixpunkte vorab festgelegt; zwingend nur ϖ_R(Zelle) ≤ Λ ⇔ Q ≥ 5",
       "Teil im Modell", "Schluss oben"),
     V("modulraum_zahlen", "2026-09-29", "Q über Phasenraum und Kapazität, m_p/m_e an CM-Punkten, Zahlensuche "
       "(1,47·10⁶ Ausdrücke)", "trägt nicht", "Modulraum-Vorschläge"),
     V("q_form", "2026-09-30", "Q als Formverhältnis: ausgezeichnete flache Tori und Einbettungen geben kein ν = 6,364; "
       "nur das Fenster 6 ≤ ν ≤ 2π…6,5", "trägt nicht", "Q als Formverhältnis", PLANCK)],
    {"Q als Einzelwert aus Formverhältnissen oder Zahlensuche": "q_form (835 Formverhältnisse) und modulraum_zahlen "
     "(1,47·10⁶ Ausdrücke) treffen wie Zufall",
     "Q = 2π − 1 über H₀τ_D": "Uhrgröße, gilt nur heute (Beobachterfehler)",
     "q aus einer Verlustrate im homogenen Hintergrund": "Γ = 0 gibt kein q",
     "Q als topologische Zahl": "die Form gibt nur ganze oder rationale Zahlen; Q ist dynamisch"},
    "Ω_c/Ω_b auf 0,15–0,5 % (heute 1,0–1,2 % [G]) trennt die Zählung je Radiant (2π − 1) von 5,36; oder eine Herleitung "
    "des offenen Teils je Umlauf aus der Form der Zelle",
    "offen. Im Modell: Fenster 5 ≤ Q aus dem Schluss, Zählung je Radiant als offene Lesart. Offen, hängt an ω_b, ω_c [G]: "
    "q_form und lambda_bruecke (Kandidatenwahl).")

# ------------------------------------------------------------------ a₀
THEMEN.neu(
    "a0_schwelle",
    "Welchen Wert hat die Schwelle a₀ = c²/2Λ, in welcher Lesart, und was bindet sie an die Zelle?",
    Invariante("Λ fest seit z ≈ 10 (keine Raumlänge): a₀ = c²/2Λ zeitfest; a₀ ≈ cH₀/2π gilt nur heute", "IG",
               "Hintergrund"),
    [V("zeitgesetz_a0", "2026-09-29", "Zeitgesetz der Schwelle aus der Horizonttemperatur, a₀(z) wandernd",
       "trägt nicht", haengt_an=("N_EFF_MESS",)),
     V("atmen_ohne_skala", "2026-09-29", "a₀ = Θ·cH/2π ohne gesetzte Skala, w(z) aus dem Atmen", "trägt nicht",
       haengt_an=("N_EFF_MESS",)),
     V("atmen_oder_zelle", "2026-09-29", "atmende Schwelle a₀(z) = cH(z)/2π gegen feste Zelle", "trägt nicht",
       haengt_an=("N_EFF_MESS",)),
     V("zelle_frueher_kosmos", "2026-09-29", "Zelle im frühen Kosmos: mit dem Horizont atmend ΔN_eff = 2,09, mitbewegt "
       "a₀(1)/a₀(0) = 2, fest Θ₀ = 1 ohne Folgen früh; Grundlage für „Λ keine Raumlänge“", "im Modell"),
     V("sparc_test", "2026-09-29", "SPARC-Scheibentest (153 Galaxien): mit errV allein 1,34, mit Fehlern von Abstand, "
       "Neigung, Υ 1,04–1,19·10⁻¹⁰", "im Modell"),
     V("schwelle_je_galaxie", "2026-09-29", "Wert der Schwelle: SPARC je Galaxie (Mittel 1,04, Streuung 0,21 dex, kein "
       "Anstieg mit Σ₀) und gemeinsam je Drittel (+0,18 ± 0,06); äußeres Feld senkt den gelesenen Wert "
       "(aussenfeld_sparc)", "im Modell", "Wert der Schwelle", skripte=("schwelle_je_galaxie", "aussenfeld_sparc")),
     V("sparc_streuung", "2026-09-29", "Streuung des Rests: Abstand und Neigung erzeugen 0,12 von 0,14 dex; a₀ je Galaxie "
       "Median 1,09", "im Modell", "SPARC, Streuung"),
     V("a0_lesungen", "2026-09-29", "a₀ profiliert mit Bootstrap 1,06 ± 0,09 ± 0,04; Tully-Fisher 1,53 erklärt, "
       "Log-Quadrate 0,95 durch Gewichtung", "im Modell", "a₀-Lesungen"),
     V("fruehe_galaxien_a0", "2026-09-29", "Genzel 2017/2020 und RC100: festes a₀ gibt f_DM(R_e) ≈ 0,1–0,3, a₀ ∝ H höbe "
       "es um 0,1–0,2 (Quellen nur als Suchtreffer)", "noch nicht im Modell", VERGESSEN)],
    {"a₀ = cH(z)/2π zu jeder Zeit (atmende Schwelle)": "Λ fest seit z ≈ 10; 87 [korr.] in 38 Notizen",
     "a₀ an SPARC mit anderer Gewichtung erneut lesen": "≥ 9 Lesungen (1,04 / 1,06 / 1,09 / 1,106); neu nur mit neuer "
     "Fehlerquelle oder neuen Daten"},
    "a₀ aus SPARC auf ±0,01·10⁻¹⁰ (nötig: Υ_d auf etwa 0,05 und Abstände); oder eine Herleitung von Λ aus ν_kosm·ΣΔs_b",
    "im Modell: Standardwert 1,1·10⁻¹⁰, SPARC 1,06 ± 0,09 ± 0,04, Zelle ohne Kreispaare ≤ 1,13. Offen, hängt an N_eff "
    "[G]: die Urteile gegen die atmende Schwelle (die Zelle im frühen Kosmos steht im Modell, ihr Vergleich mit N_eff ist "
    "ein Gedächtniswert).")

# ------------------------------------------------------------------ Zwerge und äußeres Feld
THEMEN.neu(
    "zwerge",
    "Warum weichen druckgestützte Systeme im äußeren Feld (DF2, DF4, Pal 14, NGC 2419, Fornax, Crater II, Draco) von der "
    "Regel ab?",
    Invariante("örtliche Regel ν(|g_N|) mit g_N = g_i + g_e am Ort; Vorzeichen der Abweichung ordnet sich nach g_i/a₀",
               "L", "Zentrum"),
    [V("zwerge_aussenfeld", "2026-09-29", "Feldform im homogenen äußeren Feld, Jeans-Streuung: DF2 17,7, DF4 19,2 km/s "
       "gegen 8,5 und 4,2 [G]", "im Modell", haengt_an=ZWERGE),
     V("tunnel_aussenfeld", "2026-09-29", "Rückkehr als Tunneln im äußeren Feld; Bezug der Galaxie A gegen B; Tunnelbild "
       "aus 07 wieder gestrichen, Rechnung nie eingetragen", "offen"),
     V("galaxie_als_ganzes", "2026-09-29", "Bezug der Galaxie = Massenmittel A der örtlichen Regel (zwischen g_e und "
       "ν(g_e)g_e)", "im Modell"),
     V("sparc_aussenfeld_rand", "2026-09-29", "äußerste SPARC-Punkte: Steigung der Außenkurven ohne äußeres Feld; "
       "s_kosm² = 0,0292", "im Modell"),
     V("ort_oder_system", "2026-09-29", "Zählung am Rand des Systems ohne Elternfeld: hebt alle Zwerge an (0 von 6 in 2σ); "
       "Rand als Raumradius von SPARC ausgeschlossen", "trägt nicht", "Ort oder System", ZWERGE),
     V("zwerge_uebergabe", "2026-09-29", "Zwerge vor der Übergabe: t_F = I_F/ω = 0,03–2,4 Gyr kürzer als ein Umlauf",
       "trägt nicht", "Zwerge vor der Übergabe", ZWERGE),
     V("aussenfeld_kosmisch", "2026-09-29", "äußeres Feld als Einbettung in die kosmische Rückkehr "
       "g_e = (ln q_kosm)²·a₀ = 0,029·a₀: gleichgerichtete Form, SPARC verlangt 0,003 ± 0,003·a₀", "trägt nicht",
       "Äußeres Feld als Einbettung in die kosmische Rückkehr, g_e = (ln q_kosm)²·a₀ = 0,029·a₀"),
     V("gedaechtnis_umgebung", "2026-09-29", "Gedächtnis der Rückkehr auf Bahnen im wechselnden Feld: Zeitmittel, σ "
       "höchstens 20 % verschoben, 0 von 7 getroffen", "trägt nicht",
       "Gedächtnis der Rückkehr im wechselnden äußeren Feld (Zwerge auf Bahnen)", ZWERGE),
     V("ebenengrenze_zwerge", "2026-09-29", "Ebenengrenze der Hierarchie als Umgebungsgröße: σ² höchstens 1,7 % gesenkt "
       "(DF2)", "trägt nicht", "Ebenengrenze der Hierarchie als Umgebungsgröße der Zwerge", ZWERGE),
     V("halbgerade_hierarchie", "2026-09-29", "Ebenengrenze R = (1 − √μ)/(1 + √μ) als Umgebungsgröße (Notiz 53)",
       "offen", VERGESSEN),
     V("zwerge_gezeitentensor", "2026-09-30", "Außenfeld als Gezeitentensor statt Vektor (g_e → 2g_e·R_e/D): wirkt "
       "wie ohne Außenfeld; Crater II 3,6 (gemessen 2,7 ± 0,3), Draco 4,7 (9,1), DF2 20,7 (8,5), DF4 19,4 (4,2), "
       "Pal 14 2,1 (0,38) km/s; 5 von 7 schlechter. Das Argument der Regel ist das ganze Taktgefälle am Ort, ein "
       "homogenes Feld lässt sich davon nicht trennen", "trägt nicht", haengt_an=ZWERGE, skripte=()),
     V("zwerge_fluchtmass_bahn", "2026-09-30", "Fluchtmaß längs der vergangenen Bahn als Argument der Regel: jede "
       "Gedächtnisform senkt höchstens bis Newton (ν ≥ 1); DF4, Pal 14 und NGC 2419 liegen gemessen unter Newton "
       "(4,2 < 7,15; 0,38 < 0,67; 4,1 < 5,1 km/s), Draco über dem Wert ohne Außenfeld; das Fluchtmaß folgt aus der "
       "Regel und kann sie nicht setzen (Kreis); Zeitmittel schon gerechnet (≤ 20 %)", "trägt nicht",
       haengt_an=ZWERGE, skripte=()),
     V("zwerge_geschichte", "2026-09-30", "Geschichte als Gedächtnis der Rückkehr: Aufbau nach Neustart dauert ~τ_D "
       "(bei 13,8 Mrd. J. ν ≈ 1,08 für alle), trennt Fornax nicht von DF2; die Regel ist stationär, kein Vorgang in "
       "der Zeit. Nachwirkung eines früher schwächeren Felds hebt Draco höchstens auf 4,7 km/s (gemessen 9,1), "
       "früher größere Masse senkt", "trägt nicht", haengt_an=ZWERGE, skripte=())],
    {"Einzelsysteme mit neuer Umgebungsgröße (Übergabe, Gedächtnis, Ebenengrenze, Tunnel, Rand des Systems)":
         "0 von 6 bzw. 7 getroffen; der Betrag ordnet sich nach keiner Größe",
     "Zählung am Rand des Systems (mitwandernd)": "gegen die Lokalität: Γ ist das ganze Taktgefälle",
     "kosmisches Außenfeld (ln q_kosm)²·a₀ gleichgerichtet": "SPARC verlangt 0,003 ± 0,003·a₀ (eigene Daten)"},
    "geprüfte Messwerte σ, Abstand, Masse und Umgebung der Zwerge (Netz gesperrt, heute [G]); oder eine Größe des Modells, "
    "die DF4 von Fornax trennt",
    "offen, hängt an den Messwerten der Zwerge [G]: alle Urteile „trägt nicht“ dieser Reihe. aussenfeld_kosmisch trägt "
    "nicht gegen SPARC (eigene Daten); seine frei angepasste Zahl steht im Modell (Beschluss offen).")

# ------------------------------------------------------------------ Warp
THEMEN.neu(
    "warp",
    "Was dreht die Knotenlinie des Warps der Milchstraße rechtläufig (+10,9 km/s/kpc [G])?",
    Invariante("Präzession eines geneigten Rings unter der örtlichen Regel (Kernform) gegen Feldform, Faktor 2,4–3,2",
               "L", "Wirbel"),
    [V("warp_praezession", "2026-09-29", "Warp-Präzession der Scheibe unter der Feldform: −1 bis −16 km/s/kpc "
       "(rückläufig)", "trägt nicht", haengt_an=("WARP_PRAEZESSION",)),
     V("warp_drehmoment", "2026-09-29", "Drehmoment von LMC, Sgr-Zwerg, M31, Lokaler Gruppe auf einen Ring bei 14 kpc: "
       "höchstens 0,68 km/s/kpc", "Teil im Modell", haengt_an=("WARP_PRAEZESSION",)),
     V("ringe_oertlich", "2026-09-29", "Präzession geneigter Ringe unter der örtlichen Regel gegen Feldform", "im Modell"),
     V("warp_amplitude", "2026-09-29", "Warp mit veränderlicher Amplitude: frei zerfallend −1 bis +1 km/s/kpc, beim "
       "Zerfall kurz bis +8; Messung trüge ohne Präzession bei kinematischer Amplitude 0,1–0,27 der geometrischen",
       "Teil im Modell", VERGESSEN),
     V("warp_komplexe_kopplung", "2026-09-30", "Gaseinfall als imaginärer Term in ω_p²: ein imaginärer Anteil gibt "
       "Wachsen oder Dämpfen, der Realteil der Präzession bleibt in erster Ordnung; das Vorzeichen geht im Quadrat "
       "verloren und kann aus ω_p² nicht kippen. Verlust gehört zu |q| (transversal), die Knotendrehung ist Phase "
       "(tangential); κ_Warp frei", "trägt nicht", haengt_an=("WARP_PRAEZESSION",), skripte=()),
     V("warp_mitnahme", "2026-09-30", "Laufender Einfall mit umlaufendem Muster (Treideln): w' = iω_f·w + κ(ψe^{iΩ_c t} − w) "
       "hat für jedes κ > 0 einen anziehenden Zustand mit Präzession genau Ω_c, rechtläufig, ohne Stopp bei 90°; "
       "der Warp wickelt nicht mehr auf (12 gegen 16 kpc 2,2° statt 180°); ohne Einfall rückläufig −37,2 km/s/kpc. "
       "Ω_c = 10,9 ist Eingabe (Umlauf der Scheibe bei 18,7 kpc); Amplitude 0,027–0,053 statt 0,1, verlangt "
       "Ṁ·sinψ ≈ 19 M☉/Jahr", "noch nicht im Modell", haengt_an=("WARP_PRAEZESSION",), skripte=())],
    {"Gravitationsmoment von LMC, Sgr, M31 und Lokaler Gruppe": "alle ≤ 0,68 km/s/kpc; viermal gerechnet",
     "Zählung am Rand des Systems für den Warp": "rückläufig in beiden Lesungen (ort_oder_system)"},
    "geprüfte Messung der Knotendrehung mit kinematischer Amplitude der Riesensterne (0,1–0,27 der geometrischen würde "
    "sie ohne Präzession erklären)",
    "offen, hängt an WARP_PRAEZESSION [G]; im Modell steht die Kinematik am Knoten und die rückläufigen Kippmoden.")

# ------------------------------------------------------------------ Naht und Übergabe
THEMEN.neu(
    "naht_uebergabe",
    "Wie geht die Mode ξ an der w-Schale (Umkehrradius R_ta) in die örtliche Regel über?",
    Invariante("Zählung der Umkehrungen ≥ 1 (Preisach θ̃ = min θ), Übergabe bei der ersten Faltung des Phasenraums "
               "I_F ≈ 0,9–1,26 rad; Gewicht w = 1 − e^{−√(θ̃τ_D)}", "LZ", "Wirbel"),
    [V("naht_rq_umkehr", "2026-09-29", "Mischform 1 + wQ + (1 − w)(ν(s_tot) − 1) an der w-Schale: Mulde 4,05, zählt θ̃ "
       "doppelt", "noch nicht im Modell", "Übergang an der w-Schale"),
     V("mischregel_schale", "2026-09-29", "Rücklauf mit zwei Verlusten, ρ = Γ/(Γ + θ̃)", "noch nicht im Modell",
       "Aufteilung an der w-Schale aus einem Rücklauf mit zwei Verlusten"),
     V("schale_zeitlich", "2026-09-29", "Schale zeitlich: Mulde bleibt, Linse ab 1,05·R_ta in allen Lesarten gleich",
       "Teil im Modell", "Schale zeitlich", skripte=("schale_zeitlich", "kollaps_nach_umkehr", "naht_rq_umkehr")),
     V("theta_bezug", "2026-09-29", "Bezug von θ̃: absolut gegen das feste Λ; relativ zum Ganzen (θ̃ − 3H) trägt nicht",
       "im Modell", "Bezug von θ̃"),
     V("wirkung_luecken", "2026-09-29", "Lücken der Wirkung: Wirbelanteil als Strom J, w vor Q_s", "Teil im Modell",
       "Lücken der Wirkung"),
     V("wachstum_w_qs", "2026-09-29", "Gewicht vor Q_s im linearen Wachstum: fσ₈ +9·10⁻⁶", "im Modell",
       "Gewicht vor Q_s im linearen Wachstum", skripte=("wachstum_w_qs", "gewicht_rueckkehr")),
     V("uebergabe_faltung", "2026-09-29", "Übergabe bei der ersten Faltung (Fillmore–Goldreich): I_F = 1,10–1,26 rad; "
       "Stufe statt e^{−I/N*} noch nicht im Modell", "Teil im Modell", "Übergabe bei der ersten Faltung"),
     V("kollaps_nach_umkehr", "2026-09-29", "Kollaps nach der Umkehr: Form im Modell; Umkehr-z je Masse aus dem Spektrum "
       "[G], Kernradius 0,58–0,61·R_ta", "Teil im Modell", "Kollaps nach der Umkehr"),
     V("stoerungen_rekombination", "2026-09-29", "lineare Störungen bei z = 1090: w entscheidet, nicht a_N/a₀",
       "im Modell", VERGESSEN)],
    {"Mischform 1 + wQ + (1 − w)(ν(s_tot) − 1)": "zählt θ̃ zweimal",
     "θ̃ relativ zum Ganzen (θ̃ − 3H)": "Überdichten ab z = 1090 gebunden, Peaks 6–8 % daneben",
     "statisches Budget je System (Aufteilung Operator/Rest)": "hängt am statischen Budget, trägt nicht (haufen_nach_naht)"},
    "Linse zwischen R_ta und 1,05·R_ta; Stoßnachlauf ab t_F = I_F/ω_orb an Haufenstößen",
    "Teil im Modell (I_F, Schale, Zuordnung der Ketten gesetzt); offen: ob die Übergabe λ = (1 − w)ω/N* die Mulde füllt "
    "und ob der Rest gebrochener Ketten ein Gedächtnis hat.")

# ------------------------------------------------------------------ Linse und Haufen
THEMEN.neu(
    "linse_haufen",
    "Trägt die Regel mit Übergabe die Linsen um isolierte Galaxien und das Profil der Haufen?",
    Invariante("gelesene Masse (1 + w·Q)·M_b jenseits R_ta, in R₂₀₀ Bilanz (1 + Q)·M_b; Linse = dynamische Masse für "
               "Kugeln", "LG", "Zentrum"),
    [V("haufen_nach_naht", "2026-09-29", "Haufen nach der Naht: Mode übergeben; gegen NFW-Referenz [G] Faktor 3,1 "
       "(100 kpc), in R₂₀₀ 9 %", "offen", "Haufen nach der Naht"),
     V("linse_ohne_freiheit", "2026-09-29", "Übergabe über die Bahnen einer Schale, N*_eff = 1,06–1,23; gegen KiDS "
       "(Repo-Daten) χ² 1094/48 gegen Regel 430/48", "Teil im Modell", "Linse ohne freien Parameter"),
     V("linsen_isoliert", "2026-09-29", "Linsen um isolierte Galaxien (KiDS aus einem Repo): bis 0,3 Mpc N* = 0,3–2,4, "
       "kumulativ −22/−44 % bei 500/1000 kpc", "noch nicht im Modell",
       "Linsen um isolierte Galaxien (KiDS, Brouwer 2021; Mistele 2024)"),
     V("harvey_nachlauf", "2026-09-29", "Harvey-Stapel der Haufenstöße gegen den Nachlauf: Fenster N* ≈ 0,3–1 "
       "(Quelle nicht erreichbar [G])", "offen", VERGESSEN)],
    {"NFW-Linse [G] als Maßstab des Haufenrests": "Gedächtniswert; offen ist das Profil, nicht die Bilanz"},
    "Linsen um isolierte Galaxien bei 200–1000 kpc auf etwa 5 % mit an den Originalarbeiten geprüften Daten; "
    "Haufenstöße nach t_F",
    "offen: Spannung zwischen der aus dem Modell bestimmten Übergabe und den Repo-Linsen jenseits 0,5 Mpc; Haufenprofil.")

# ------------------------------------------------------------------ Mikroskala
THEMEN.neu(
    "mikroskala",
    "Legt eine Gleichung des Modells m_p/m_e, α oder Im τ_X fest?",
    Invariante("Nome q = e^{2πiτ} auf der Teilchenebene: Betrag = Verlust, Phase = Windung (Ladung, Spin); Leiter "
               "λ̄ : ℓ_P : 3Δs mit Spiegel M ↔ m_Pl²/M", "LZ", "Hintergrund"),
    [V("mikro_nome", "2026-09-29", "Triade auf der Teilchenebene; drei Knoten je g = 2 geben das Quarkmodell mit gesetzten "
       "Drittelladungen; AME2020 ohne Leiter", "Teil im Modell", "Mikroskala und Nome"),
     V("hopf_knoten", "2026-09-29", "Windung plus Versatz i·a_J als inverse stereographische Abbildung ℝ³ → S³; "
       "Hopf-Zahl am Gitter und über die Gauß-Verlinkung", "noch nicht im Modell", "Hopf-Knoten"),
     V("em_phasenseite", "2026-09-29", "EM als Phasenseite: dF = 0, g = 2 aus dem imaginären Versatz, μ = μ_B",
       "im Modell", "EM als Phasenseite"),
     V("tau_ebenen", "2026-09-29", "τ je Ebene und Spiegelung: Leiter symmetrisch um r*, Regel nicht; keine SL(2,ℤ)-"
       "Ordnung", "Teil im Modell", "τ je Ebene und Spiegelung"),
     V("spiegel_m22", "2026-09-29", "Spiegel M ↔ m_Pl²/M als gespiegelte Quantendiffusion (Gesetz M22): ħ kürzt sich",
       "noch nicht im Modell", "Spiegel M ↔ m_Pl²/M (λ̄ ↔ 3Δs)", skripte=("ftm_gesetze",)),
     V("fusion_trichter", "2026-09-29", "Tiefe und Temperatur im Stern, Brennstufen als Leiter in δ; Zündtemperatur, "
       "Mindestmasse des Brennens und Häufigkeiten tragen nicht", "Teil im Modell",
       "Tiefe und Temperatur im Stern, Brennstufen als Leiter in δ")],
    {"m_p/m_e an CM-Punkten oder per Zahlensuche": "das Modell erzwingt keinen CM-Punkt; Treffer wie Zufall",
     "Gleichsetzung a_J = MHD-Helizität": "Helizität des eigenen Feldes ist null",
     "Quarkmodell-Drittel als Herleitung": "Drittelladungen gesetzt"},
    "eine Gleichung des Modells mit Viererterm, die Spin ½ bei H = 1 und eine Stabilitätsschranke trägt (Derrick), oder "
    "die Im τ_X festlegt",
    "offen: keine Gleichung des Modells legt m_p/m_e oder α fest.")

# ------------------------------------------------------------------ Zellform
THEMEN.neu(
    "zellform",
    "Welche Form hat die Zelle (Würfel, sechseckig, raum-, flächenzentriert)?",
    Invariante("(1 + Q)·a₀ ≤ 3c⁴/(2Gρ_b L_min³)·f³ mit Formfaktor f³ je Zelle; Drehsymmetrien flacher Tori", "G",
               "Hintergrund"),
    [V("zelle_form", "2026-09-30", "Zellform im eigenen Atmen: d_LSS ×1,0209, Würfel auf der Grenze; Ausrichtung der "
       "Quadrupolleistung trennt", "im Modell", "Zellform"),
     V("q_a0_kopplung", "2026-09-29", "Q und a₀ gekoppelt: (1 + Q)·a₀ ≤ 7,173·10⁻¹⁰·f³", "im Modell",
       "Q und a₀ gekoppelt"),
     V("winkelmoden_zelle", "2026-09-30", "Winkelmoden der Zelle in forschung/kosmos.py (winkelmoden_zelle): Quadrupol "
       "15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0; Fenster 17,0 % bei "
       "0,985·d_LSS und 2,4 % bei 1,5·d_LSS; Kreispaare ab 1,013 (Würfel) bzw. 1,170·10⁻¹⁰ (sechseckig)",
       "Teil im Modell", haengt_an=PLANCK, skripte=())],
    {"Zellform am Betrag der Quadrupolsenkung": "alle Formen 0,83–0,85; nur die Ausrichtung trennt"},
    "Kreispaare oder Ausrichtung der Winkelmoden im Himmel; a₀ auf ±0,01·10⁻¹⁰",
    "im Modell ohne festgelegte Form; hängt an Q [G] und dem Lauf des eigenen Atmens Ẋ_an(t).")

# ------------------------------------------------------------------ Zustandsraum und Labor
THEMEN.neu(
    "zustandsraum_labor",
    "Welche Aussagen des Zustandsraums sind prüfbar, und gibt es eine Laborvorhersage nur des Modells?",
    Invariante("drei Rollen (GENERIC), Schließung durch Translation, Rückkopplung als Resolvente 1/(1 − q); Maß "
               "erhalten (doppelt stochastisch)", "Z", "Wirbel"),
    [V("labor_vorhersage", "2026-09-29", "Sattelpunkt Sonne–Erde Newtonsch, Uhren: galaktische Gezeit 3·10⁻⁴¹ über "
       "1 km; Sattel-Knoten-Normalform", "im Modell", "Laborvorhersage"),
     V("kontakt_atmen", "2026-09-29", "Kontaktgeometrie und Festkörper: kein w(z), RC-Kette als exakte Nachbildung",
       "Teil im Modell", "Kontaktgeometrie und Festkörper"),
     V("axiom_mass", "2026-09-30", "zwei gekoppelte Tori: Maß erhalten, gebundenes Maß = Separatrixfläche", "im Modell",
       "Axiome und Maß"),
     V("gesetzespruefung", "2026-09-29", "Gesetzesprüfung (164 Gesetze): Grenzgesetze, Zustandsraum, Tafeln des Autors, "
       "Galaxiengesetze", "Teil im Modell", "Gesetzesprüfung", skripte=("ftm_gesetze",)),
     V("generic_praezision", "2026-09-30", "GENERIC-Präzisierung: Energie braucht M∇E = 0, der zweite Hauptsatz L∇S = 0; "
       "das parametrische Glied läuft über ∂_tH (dE/dt = ∂_pH·ṗ); Kupfer-Probe (1000 Zufallszustände: |L∇S| = 0, "
       "|M∇E| ≤ 1,6·10⁻¹⁷; ΔE/E₀ 3,2·10⁻¹⁴, ΔS = +1,5·10⁻⁶ J/K, ∫RI²dt = ∫T dS)", "Teil im Modell", skripte=()),
     V("noether_form", "2026-09-30", "Noether-Form μ̇ = −dH(ξ_M) (Abweichung 1,1·10⁻¹⁶); Erhaltung nur bei H, S, L, M "
       "invariant; die Form F = −∇⟨μ,ξ⟩ ist gesperrt (der Gradient ist die Symmetrierichtung, keine Kraft; Einheit "
       "Impuls/Länge)", "noch nicht im Modell", skripte=())],
    {"Impulsabbildung F = −∇⟨μ,ξ⟩ als Kraft": "Gradient von ⟨μ,ξ⟩ ist die Symmetrierichtung ξ_M, keine Kraft; richtig ist μ̇ = −dH(ξ_M)",
     "Laborvorhersage aus der Regel": "im Labor s ≈ 3·10⁵, galaktischer Anteil als Gezeit ≈ 10⁻²⁶ über 1 km",
     "Weyl, ³He-B, Magnonen als Analogon": "keine Halbgerade mit Erneuerung"},
    "ein Aufbau, der eine Aussage nur des Modells trägt (bisher keiner); sonst die Form (Sattel-Knoten, Wurzelgesetz, "
    "Kac-Budget) an Laborsystemen",
    "Form im Modell, keine eigene Laborvorhersage.")

# ------------------------------------------------------------------ Zeit und Gedächtnis
THEMEN.neu(
    "zeit_gedaechtnis",
    "Welche Zeitskala hat das Gedächtnis der Regel (Kern K, λ_K, λ_m = ω_orb/N*), und wie hängt die Zeit an der inneren "
    "Richtung?",
    Invariante("Kern aus derselben Rückkehr K̂(p) = e^{−√((p + Γ)τ_D)}; geschützt ist die Zählung der Umkehrungen; "
               "L und I bilden eine komplexe Zeit", "LI", "Wirbel"),
    [V("kern_gedaechtnis", "2026-09-29", "Maxwell-Form des Gedächtniskerns: alte Scheiben verlangen λ_K ≲ 3·10⁻⁴ Ω "
       "(Literatur); ersetzt durch den Kern aus der Rückkehr", "offen", "Gedächtniskern"),
     V("kopfgleichung_kern", "2026-09-29", "Kern aus derselben Rückkehr (Lévy-Dichte·e^{−ΓT}); am Sonnenkreis 2·10⁻¹⁵",
       "im Modell", "Kern aus derselben Rückkehr", skripte=("kopfgleichung_kern", "kern_gedaechtnis")),
     V("gedaechtnis_topologie", "2026-09-29", "Preisach θ̃ = min θ gleichwertig mit Zahl der Umkehrungen ≥ 1; Stufe bei "
       "z_c = 2 44,6 Mio. J.", "Teil im Modell", "Gedächtnis als Topologie"),
     V("innere_richtung_zeit", "2026-09-29", "innere Richtung als imaginäre Zeit; Bose-Form aus T_I und Unruh-Lesart "
       "tragen nicht", "Teil im Modell", "Innere Richtung als imaginäre Zeit")],
    {"Maxwell-Form, Preisach als Echtzeit, Mittag-Leffler für den Kern": "der Kern folgt aus der Rückkehr (Lévy); "
     "w = 1 − e^{−√x}",
     "N* = ½√(ωτ_D)": "trägt nicht"},
    "N* aus Haufenstößen und Linsen jenseits 200 kpc; λ_K aus der Radialstreuung alter Scheiben mit geprüften Daten",
    "Kern im Modell; N* offen (0,3–2,4).")

# ------------------------------------------------------------------ Feldform und Scheiben
THEMEN.neu(
    "feldform_scheiben",
    "Wirkt die Regel örtlich auf der Bahn (Kernform) oder als Feld (Feldform), und was zeigen Scheiben?",
    Invariante("örtliche Regel a = ν(|g_N|)·g_N am Ort der Bahn, in der Mittelebene f = 0; Metrik trägt den wirbelfreien "
               "Teil", "L", "Zentrum"),
    [V("sparc_boseform", "2026-09-29", "Bose-Form an den Punkten; größte Zusatzbeschleunigung bei 2,54a₀ an SPARC nicht "
       "trennbar (Υ_d, Bulge)", "offen", "Größte Zusatzbeschleunigung bei a_N = 2,54a₀ (0,648a₀)"),
     V("sparc_flaechendichte", "2026-09-29", "zentrale Flächendichte (Toomre) Galaxie für Galaxie: 0,05 dex unter der "
       "Regel, Trend +0,05 ± 0,025", "im Modell",
       "Zentrale Flächendichte (Toomre), Galaxie für Galaxie mit den eigenen Baryonen"),
     V("gleichung_ganzes", "2026-09-29", "Projektion auf den Ort: lokal exakt; nichtlokal |k| trägt nicht (χ²/N 36 gegen "
       "16)", "im Modell", "Projektion auf den Ort"),
     V("trenntests_b_dm", "2026-09-30", "Trenntests Bindung gegen Dunkle Materie: η = 1 für Kugeln, Scheiben innen "
       "1,10–1,91; Gasanteil ±0,013 dex", "Teil im Modell", "Trenntests Bindung gegen Dunkle Materie"),
     V("sparc_universal", "2026-09-29", "Rest je Galaxie 0,16 dex ohne Zusammenhang mit Eigenschaften; gemeinsames "
       "äußeres Feld 0,02–0,03·a₀ hebt den Anstieg auf", "im Modell", "SPARC, Universalität")],
    {"Feldform (f = 1) als Bahnregel": "SPARC schließt sie außerhalb 95 % aus",
     "nichtlokale Fortsetzung |k| (Caffarelli-Silvestre)": "v⁴ = 0,405·GMa₀, SPARC χ²/N 36 gegen 16",
     "ideale Exponentialscheibe als Maßstab der Flächendichte": "Galaxie für Galaxie trägt; die Idealscheibe liegt daneben"},
    "Gewicht f des senkrechten Felds an SPARC (±0,03); Linsen durch Spiralgalaxien innen (M_dyn/M_lens > 1); Gaia "
    "∂_z a_R ≠ ∂_R a_z",
    "im Modell (örtliche Regel); offen die Vorhersagen in C.")

# ------------------------------------------------------------------ Rückkehr als eine Gleichung
THEMEN.neu(
    "rueckkehr_gleichung",
    "Folgt die Regel ν = 1 + n_B(s) aus einer Gleichung mit wenigen Setzungen, und welche Invariante trägt sie?",
    Invariante("∂_t u + θu = D(∇² + ∂_ϖ²)u auf der Halbgeraden, eindimensional (d_s = 1), Kac ⟨T⟩ = 2Λ/c; ϖ in der "
               "inneren Richtung", "IZ", "Wirbel"),
    [V("gleichung_zeit", "2026-09-29", "Zeitform der einen Gleichung: parabolisch gibt q, 1 − w, K̂; Welle trägt nicht",
       "Teil im Modell", "Zeitform der einen Gleichung"),
     V("setzungen_gleichung", "2026-09-29", "Setzungen: Irrfahrt mit einem Takt gibt parabolische Form, eine Uhr, einen "
       "Sinn", "im Modell", "Setzungen der Gleichung"),
     V("eine_struktur", "2026-09-29", "eine Struktur: Nome mit Betrag = Verlust, Phase = Windung; CHSH 2√2 mit "
       "Grundschwingung", "im Modell", "Eine Struktur"),
     V("info_geometrie", "2026-09-29", "Fisher-Metrik g = ν(ν − 1); χ keine Geodäte; S = A/4ℓ_P²", "Teil im Modell",
       "Informationsgeometrie und Entropie"),
     V("strittig_zweiter_weg", "2026-09-29", "strittige Punkte auf zweitem Weg: Hopf-Zahl, Tori, χ, Clausius, "
       "Stabilisatoren, Mori-Kettenbruch", "offen", "Strittige Punkte, zweiter Weg"),
     V("schroedinger_m23", "2026-09-29", "Rückkehrgleichung als Schrödinger-Gleichung bei E = 0 (Gesetz M23): ħ kürzt "
       "sich; daraus 0,65·Λ als kleinster Eigenwert", "Teil im Modell",
       "Rückkehrgleichung als Schrödinger-Gleichung bei E = 0", skripte=("ftm_gesetze",)),
     V("tunnelbild_gesetz", "2026-09-30", "Tunnelbild ψ(r*) = q = e^{−s} (Milchstraße): rechnet auf 10⁻⁴, steht nicht im Modell; "
       "das Gesetz M_rueckkehr_tunneln ist Forschung (forschung/tunnelbild.py, Test tests/test_tunnelbild.py)",
       "noch nicht im Modell", "Rückkehrgleichung als Schrödinger-Gleichung bei E = 0"),
     V("gemeinsame_innere_richtung", "2026-09-29", "Bündel je Falte gegen gemeinsamen Knoten bei Λ: q = 1 − ΣΔs_b/Λ ist "
       "Λ = ν_kosm·ΣΔs_b umgeschrieben", "Teil im Modell",
       "Eine oder viele innere Richtungen (Bündel je Falte gegen gemeinsamen Knoten bei Λ)"),
     V("hannay_holonomie", "2026-09-30", "Hannay-Winkel π(cosh r − 1) einer geschlossenen Parameterschleife (Ebene Z, "
       "geschwindigkeitsunabhängig, null bei einem Parameter) und q als multiplikative Holonomie des Kopplungsgraphen; "
       "Forschung, ohne Modellfall mit zwei nichtvertauschenden Parametern", "noch nicht im Modell", skripte=())],
    {"ν = 1 + n_B(s) in neuer Sprache (Kac, Lévy, Irrfahrt, Tunnel, Schrödinger, imaginäre Zeit, Wicklung, Projektion, "
     "Topologie)": "> 20 Umformungen derselben Regel ohne neue Invariante",
     "Wellenform der Rückkehr": "w ∝ √M, keine Dämpfung bei Ω",
     "Rückkehrrichtung als Wicklung oder irrationale Drehung": "höchstens drei Rückkehrzeiten, keine Wurzel"},
    "eine neue Invariante oder eine Vorhersage, die eine Umformung vom Bestand trennt",
    "im Modell; weitere Umformungen gesperrt.")

# ------------------------------------------------------------------ Hierarchie der Trichter
THEMEN.neu(
    "trichter_hierarchie",
    "Wie glättet die Hierarchie der Trichter (ρ > ρ_J) von Körper zu Galaxie zu Haufen?",
    Invariante("Glättungslänge ℓ³ = GM/T mit der Gezeit des geglätteten Felds; Triade als Invarianten der Glättung "
               "(Fluss, harmonischer Elternteil, eigene Gezeit ∝ ϖ⁻³)", "LG", "Zentrum"),
    [V("ebenen_projektion", "2026-09-29", "Ebenen als Projektion: Tiefe je Körper (Mond 482 AE, Erde 2088 AE, Sonne "
       "0,70 pc)", "Teil im Modell", "Ebenen als Projektion"),
     V("trichterkette_gewicht", "2026-09-29", "Trichterkette der Sonne (Sgr A* bis Shapley); Lokale Gruppe linear "
       "mitatmend; Massen [G], aus dem Modell genommen", "noch nicht im Modell",
       "Trichterkette der Sonne (Sgr A*, Milchstraße, Sgr-Zwerg, LMC, M31, Virgo, Großer Attraktor, Shapley)",
       skripte=("trichterkette_gewicht", "trichterkette_sonne")),
     V("wep_glaettung", "2026-09-29", "Gleichheit des Fallens mit der örtlichen Regel: Δa = 0 am Schwerpunkt, gemittelt "
       "1,4–1,9·10⁻¹⁶ m/s²", "Teil im Modell",
       "Gleichheit des Fallens und Glättungsschwelle mit der örtlichen Regel"),
     V("trichter_invarianten", "2026-09-29", "Radien und Klassen des Trichters: Schnittpunkt M× = a₀³/(GT²)",
       "Teil im Modell", "Radien und Klassen des Trichters über alle Ebenen",
       skripte=("trichter_invarianten", "wep_glaettung"))],
    {"Hill-Radius als Glättungslänge": "trägt nicht (wep_glaettung)",
     "örtliche Regel auf Mpc-Skalen (Lokale Gruppe 1/r)": "Beobachterfehler; linear mitatmend"},
    "ℓ im Zentrum eines Wirts; Massen der Trichterkette aus eigenen Daten",
    "Teil im Modell; Trichterkette offen, hängt an Massen [G].")
