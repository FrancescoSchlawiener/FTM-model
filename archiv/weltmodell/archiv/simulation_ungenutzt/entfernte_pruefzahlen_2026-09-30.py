
# ---- aus pruefung/zahlen_09.py (Modelltext trägt die Zahl nicht mehr)
_P("09_milchstrasse_verstaerkung", B, "09_frequenzleiter", "1,582", lambda p: K.nu(p.a0, p.a0), "1",
   "Milchstraße (6·10¹⁰ M☉)")

# ---- aus pruefung/wiederholungen.py (Modelltext trägt die Zahl nicht mehr)
_w(datei("03"), "0,4", "Offen sind die übrigen 18", "DE_MAX_Z")
_w(datei("03"), "5,5", "Das Fenster 5 ≤ Q ≤ 5,5", "03_q_max_wuerfel")
_w(datei("03"), "5,36", "Die Zelle bindet beide Messgrößen", "03a_Q_kosmisch")
_w(datei("03"), "1,13·10⁻¹⁰", "Die Zelle bindet beide Messgrößen", "03_a0_max_wuerfel")

# ---- aus pruefung/zahlen_05.py (Modelltext trägt die Zahl nicht mehr)
AUSSAGEN.neu("05_atmen_angekommen", (D, "als Atmen der Raumschleifen an (heute 4,7 %)"), "Zahl")
PRUEFZAHLEN.neu("05_atmen_anteil_heute", B, "05_atmen_angekommen", "4,7",
                None, "%",
                "Atmen der Raumschleifen", offen="nur_forschung",
                beschreibung="angekommener Anteil X_an/X heute aus forschung.lauf_eigenes_atmen: V0 (Bildung im "
                             "Vergleichshintergrund) 4,71 %, selbstkonsistent S2 5,20 %")

# ---- aus pruefung/zahlen_07.py (Modelltext trägt die Zahl nicht mehr)
_P("07_nu_schwelle", B, "07_umlaufsinn", "1,582", lambda p: K.nu(p.a0, p.a0), "1", "bei a_N = a₀ 1,313 statt 1,582")
_P("07_rstern_zelle_zu_lambda", B, "07_netz_aller_falten", "0,97",
   lambda p: K.r_stern(K.schliesslaenge(p.a0) / (1 + K.Q_KOSM) / K.ds(1.0), p.a0) / K.schliesslaenge(p.a0), "Λ",
   "ϖ_R = r* = 0,97·Λ", haengt_an=("OMEGA_B_H2", "OMEGA_C_H2"))

# ---- aus pruefung/zahlen_offen.py (Modelltext trägt die Zahl nicht mehr)
_offen("03", 3, "gemeinsames Feld", ("0,010",), NF,
       "gemeinsames äußeres Feld frei angepasst, senkrecht (aussenfeld_kosmisch, Forschungsstand: trägt nicht); fehlt Formel")
_offen("03", 4, "gemeinsames Feld", ("0,009",), NF,
       "Fehler des frei angepassten senkrechten Felds")
_offen("03", 5, "gemeinsames Feld", ("0,003",), NF,
       "gemeinsames äußeres Feld frei angepasst, gleichgerichtet, mit Fehler")
_offen("03", 44, "ab der Übergabe I_F", ("0,32", "0,99", "0,36", "0,56",), NF,
       "Fluchtmaß der Dunklen Energie im Standardmodell gelesen (Rechnungen de_fluchtmoden, de_rest_anstieg, modell_h0: trägt nicht bzw. noch nicht im Modell)")
_offen("03", 45, "Maximum erst bei z", ("18", "1,01", "0,61", "0,01",), NF,
       "Fluchtmaß der Dunklen Energie im Standardmodell gelesen (Rechnungen de_fluchtmoden, de_rest_anstieg, modell_h0: trägt nicht bzw. noch nicht im Modell)")
_offen("03", 46, "letzte Streufläche", ("2,1", "28,4", "6,74·10⁻¹⁰",), NF,
       "forschung.lauf_eigenes_atmen (S2, selbstkonsistent): d_LSS 2,09 % ferner (28,37 Gpc), (1 + Q)·a₀ ≤ 6,741·10⁻¹⁰·f³ "
       "(V0: 2,23 %, 28,41 Gpc, 6,714)")
_offen("03", 47, "sechseckige Zelle liegt", ("5,13", "1,7", "3,2", "4,4",), NF,
       "forschung.lauf_eigenes_atmen (S2): Würfel Q ≤ 5,128 bei a₀ = 1,1; Abstände mit a₀ = 1,06 ± 0,09 ± 0,04 (SPARC-Lesung, "
       "Kernformel fehlt) und Q = 5,36 ± 0,07: sechseckig 1,65σ, raumzentriert 3,20σ, flächenzentriert 4,43σ")
_offen("03", 72, "Form w(z) des Fluchtmaßes", ("18",), NF,
       "Fluchtmaß der Dunklen Energie im Standardmodell gelesen (Rechnungen de_fluchtmoden, de_rest_anstieg, modell_h0: trägt nicht bzw. noch nicht im Modell)")

# ---- aus pruefung/zahlen_03b.py (Modelltext trägt die Zahl nicht mehr)
_A("03_de_anteil_standard", (D, "zusammen mit dem Vakuumteil 0,56 von 0,687"), "Zahl")
_P("03_omega_de_standardmodell", B, "03_de_anteil_standard", "0,687", lambda p: K.omega_de(), "1",
   "von 0,687", haengt_an=("OMEGA_R",) + KOSMOS)
_P("03_a0_anteil_grenze", B, "03_offen_herkunft_q", "98", lambda p: 100 * p.a0 / (K.zellgrenze(1.0) / (1 + K.Q_KOSM)),
   "%", "liegt heute bei 98 % seiner Grenze", haengt_an=("KREISPAAR_GRENZE",) + KOSMOS)

# ---- aus pruefung/zahlen_01.py (Modelltext trägt die Zahl nicht mehr)
AUSSAGEN.neu("01_tabelle_erde_sonne", (D, "| Δs | 1,478 mm | 492,2 m |"), "Zahl")
PRUEFZAHLEN.neu("01_ds_sonne", B, "01_tabelle_erde_sonne", "492,2", lambda p: K.ds(K.M_sun), "m", "| Δs |",
                haengt_an=("GM_sun",))
PRUEFZAHLEN.neu("01_ds_erde", B, "01_tabelle_erde_sonne", "1,478", lambda p: K.ds(K.M_earth) * 1e3, "mm", "| Δs |",
                haengt_an=("GM_earth",))
AUSSAGEN.neu("01_tabelle_g", (D, "| g | 9,81 m/s² | 274 m/s² |"), "Zahl")
PRUEFZAHLEN.neu("01_g_sonne", B, "01_tabelle_g", "274", lambda p: K.neigung(K.M_sun, K.R_sun), "m/s²", "| g |",
                haengt_an=("GM_sun", "R_sun"))

# ---- aus pruefung/zahlen_galaxie.py (Modelltext trägt die Zahl nicht mehr)
_p("03_sparc_a0_varianten_unten_skaliert", "03", "03_sparc_a0_profil", "1,00",
   lambda p: min(_a0_varianten(p)[0].values()) / 1e-10, "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0, satz="Upsilon_skaliert")
_p("03_sparc_a0_varianten_oben_skaliert", "03", "03_sparc_a0_profil", "1,13",
   lambda p: max(_a0_varianten(p)[0].values()) / 1e-10, "10⁻¹⁰ m/s²", _AS, haengt_an=_H_A0, satz="Upsilon_skaliert")
_w(D03, "0,53", "| In der Feldform liegen", "07_feldform_viertel_tief")
_w(D03, "1,06", "Die Zelle bindet beide Messgrößen", "03_sparc_a0_profil_wert")
_w(D03, "0,085", "Die Zelle bindet beide Messgrößen", "03_sparc_a0_profil_stat")
_w(D03, "0,11", "Die Zelle bindet beide Messgrößen", "03_sparc_a0_profil_syst")
_w(D03, "0,96", "Die Zelle bindet beide Messgrößen", "03_sparc_a0_varianten_unten")
_w(D03, "1,19", "Die Zelle bindet beide Messgrößen", "03_sparc_a0_varianten_oben")
_w(D03, "1,06", "- Dunkle Energie: ob Λ die", "03_sparc_a0_profil_wert")
_w(D03, "1,06", "mit a₀ = 1,06 (SPARC)", "03_sparc_a0_profil_wert")
