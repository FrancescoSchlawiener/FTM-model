# Bericht der FTM-Simulation

Commit f23e388 (mit nicht committeten Änderungen) · 2026-09-30T20:58:46 · Laufzeit 70.7 s

Aufgaben 91 in 4 parallelen Prozessen, davon 0 aus dem Zwischenspeicher (Schlüssel 70cdc81c3d272a7b); Rechenzeit aller Aufgaben 279 s; längste frisch gerechnete Aufgaben: Prüfzahlen 61 50 s, Prüfzahlen 55 29 s, Prüfzahlen 18 27 s, Gegenproben 4 22 s, pytest Gesetze 20 s.

Parametersätze:

- standard: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute
- horizont: a₀ horizont (1.04e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute (eigener Horizont)
- z1090: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 1090.0, gebunden, Beobachter hier, bei der Rekombination
- Upsilon_0.4: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.4, z = 0.0, gebunden, Beobachter hier, heute
- Upsilon_0.6: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.6, z = 0.0, gebunden, Beobachter hier, heute
- Upsilon_skaliert: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute, Υ_b mitskaliert [W]
- mitatmend_heute: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, mitatmend, Beobachter zwischen den Strukturen, heute
- Sonnenkreis_Punktmasse: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter Sonnenkreis, Milchstraße als Punktmasse, Sonnenkreis punktmasse
- Sonnenkreis_gemessen: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter Sonnenort, gemessen, Sonnenkreis gemessen

[G] = hängt an einem Messwert aus dem Gedächtnis (nicht nachgeprüft); Naturkonstanten (CODATA, IAU) tragen kein [G]. [W] = hängt an einer Wahl, die nicht im Modelltext steht (gewählte Eingabe, Lesart, Zuordnung, die der Nutzer noch bestätigen muss). Beide werden vererbt; die Spalte nennt die Wurzeln.

## Wächter

pytest: grün

| Wächter | Art | Funde | je Datei / erste Funde |
|---|---|---|---|
| symbole_eindeutig | streng | 0 |  |
| registrierung | streng | 0 |  |
| literale | streng | 0 |  |
| quellen | streng | 0 |  |
| tempsim | streng | 0 |  |
| modellsatz | streng | 0 |  |
| gegenproben | streng | 0 |  |
| doppelformeln | streng | 0 |  |
| katalog | streng | 0 |  |
| katalog_luecken | streng | 0 |  |
| themen | streng | 0 |  |
| epochen_eingabe | streng | 0 |  |
| zufall | streng | 0 |  |
| beschluesse | streng | 0 |  |
| daten | streng | 0 |  |
| gegenproben_abdeckung | zählend | 11 | ds_eff; ell0; gezeitendichte |
| rueckrichtung | zählend | 0 | 01: 0, 02: 0, 03: 0, 04: 0, 05: 0, 06: 0, 07: 0, 08: 0, 09: 0 |
| symbolkollision | zählend | 89 | 00-bauplan.md: τ_X kollidiert mit τ = Modul des Torus; q = e^{2πiτ}; 00-bauplan.md: α₁ kollidiert mit α = Feinstrukturkonstante r_e/λ̄_e; 00-bauplan.md: α₂ kollidiert mit α = Feinstrukturkonstante r_e |
| themen_luecken | zählend | 0 |  |
| uhr_muster | zählend | 3 | 03-von-der-galaxie-bis-zum-kosmos.md: „Λ = πc/H“ ohne „heute“; 03-von-der-galaxie-bis-zum-kosmos.md: „a₀ = cH“ ohne „heute“; 04-zeit-verschraenkung-syntropie-gedaechtnis.md: „H_s = π²c/(√45·Λ) = H₀“ o |

## Kern-Selbsttest

| Register | Stand |
|---|---|
| Symbole | 513 |
| Formeln | forschung 9, hilfe 37, modell 260 |
| Formeln mit [G] im Eingang | 18 |
| Messwerte | 142 (davon [G] 101; Konstanten CODATA und IAU zählen nicht) |
| Setzungen | 9 |
| Eingaben | 180 |
| Parametersätze | 9 |
| Formeln mit Standardargumenten ausgewertet | 40 |

## Zusammenfassung

| Baustein | Prüfzahlen | OK | ABWEICHUNG | OFFEN | davon fehlt Formel | davon nur Forschung |
|---|---|---|---|---|---|---|
| 01 | 4 | 3 | 0 | 1 | 0 | 1 |
| 02 | 76 | 71 | 0 | 5 | 0 | 5 |
| 03 | 274 | 219 | 0 | 55 | 31 | 24 |
| 04 | 20 | 20 | 0 | 0 | 0 | 0 |
| 05 | 13 | 13 | 0 | 0 | 0 | 0 |
| 06 | 24 | 20 | 0 | 4 | 4 | 0 |
| 07 | 112 | 86 | 0 | 26 | 25 | 1 |
| 08 | 45 | 45 | 0 | 0 | 0 | 0 |
| 09 | 21 | 20 | 0 | 1 | 1 | 0 |
| alle | 589 | 497 | 0 | 92 | 61 | 31 |

Rückrichtung je Fundstelle (erfasst = Prüfzahl, auch offen, Messwert, Eingabe, Setzung, Ausnahme oder Wiederholung; gedeckt = erfasst ohne offene Prüfzahlen):

| Datei | Zahlen | erfasst | gedeckt | offen | nicht erfasst |
|---|---|---|---|---|---|
| 01 | 5 | 5 | 4 | 1 | 0 |
| 02 | 89 | 89 | 84 | 5 | 0 |
| 03 | 407 | 407 | 352 | 55 | 0 |
| 04 | 27 | 27 | 27 | 0 | 0 |
| 05 | 23 | 23 | 23 | 0 | 0 |
| 06 | 26 | 26 | 22 | 4 | 0 |
| 07 | 146 | 146 | 120 | 26 | 0 |
| 08 | 69 | 69 | 69 | 0 | 0 |
| 09 | 23 | 23 | 22 | 1 | 0 |
| alle | 815 | 815 | 723 | 92 | 0 |

Gesetze: OK 230, SPANNUNG 2, WIDERSPRUCH 1, HINWEIS 1 (zusammen 234).

Herkunftsmarken ([G] nur aus Messwert-Einträgen, ohne Konstanten; [W] aus gewählten Eingaben, Lesarten und Zuordnungen; 101 Messwerte tragen [G]):

| Art | Zahl | mit [G] | mit [W] |
|---|---|---|---|
| Prüfzahlen | 589 | 148 | 12 |
| Gesetze | 234 | 25 | 10 |
| alle | 823 | 173 | 22 |

Wurzeln von [W]: DRAHT_Q_MERKLICH, GPB_A, ISS_HOEHE, J_MILCHSTRASSE, R_HAUFEN_HUELLE, R_MILCHSTRASSE_HUELLE, TORUS_DOMEGA, TORUS_K, T_KOERPER, Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt).

Zählende Wächter: gegenproben_abdeckung 11, rueckrichtung 0, symbolkollision 89, themen_luecken 0, uhr_muster 3.

## Prüfzahlen

OFFEN: die Kernformel fehlt (fehlt Formel) oder die Zahl beruht nur auf einer Forschungsrechnung (nur Forschung); weder OK noch ABWEICHUNG.

### Baustein 01

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 01_ds_je_erdmasse | 1,4783 | 1.47834 | mm | 4.26797e-05 | 5e-05 | standard | OK |  |
| 01_r0_erde | 8,87 | 8.87006 | mm | 5.60782e-05 | 0.005 | standard | OK |  |
| 01_r0_sonne | 2,95 | 2.95325 | km | 0.00325008 | 0.005 | standard | OK |  |
| 01_offen_01_1 | 10⁻⁶ |  |  |  |  | standard | OFFEN (nur Forschung: Scherung der Schleifen unter 10⁻⁶ der Atemrate (Rechnung vakuum_richtungen)) |  |

### Baustein 02

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 02_ds_proton | 4,14·10⁻⁵⁵ | 4.14039e-55 | m | 3.86158e-59 | 5e-58 | standard | OK |  |
| 02_ds_R_proton | 4,9·10⁻⁴⁰ | 4.92083e-40 | 1 | 2.08298e-42 | 5e-42 | standard | OK |  |
| 02_g_proton | 1,6·10⁻⁷ | 1.57688e-07 | m/s² | -2.31206e-09 | 5e-09 | standard | OK |  |
| 02_ds_R_wasserstoff | 7,8·10⁻⁴⁵ | 7.8242e-45 | 1 | 2.41959e-47 | 5e-47 | standard | OK |  |
| 02_bindung_wasserstoff | 9·10²² | 9.04422e+22 | m/s² | 4.42161e+20 | 5e+21 | standard | OK |  |
| 02_ds_erde | 1,478 | 1.47834 | mm | 0.00034268 | 0.0005 | standard | OK |  |
| 02_ds_R_erde | 2,3·10⁻¹⁰ | 2.32042e-10 | 1 | 2.04249e-12 | 5e-12 | standard | OK | [G: R_earth] |
| 02_g_erde | 9,82 | 9.82025 | m/s² | 0.000250487 | 0.005 | standard | OK | [G: R_earth] |
| 02_ds_sonne | 492,2 | 492.208 | m | 0.00834608 | 0.05 | standard | OK |  |
| 02_ds_R_sonne | 7,1·10⁻⁷ | 7.07501e-07 | 1 | -2.49914e-09 | 5e-09 | standard | OK |  |
| 02_g_sonne | 274 | 274.2 | m/s² | 0.200112 | 0.5 | standard | OK |  |
| 02_ds_neutronenstern | 689 | 689.092 | m | 0.0916845 | 0.5 | standard | OK |  |
| 02_ds_R_neutronenstern | 0,057 | 0.0574243 | 1 | 0.000424307 | 0.0005 | standard | OK |  |
| 02_g_neutronenstern | 1,3·10¹² | 1.29026e+12 | m/s² | -9.74017e+09 | 5e+10 | standard | OK |  |
| 02_ds_roter_riese | 8,86 | 8.85975 | km | -0.00024977 | 0.005 | standard | OK | [G: BETEIGEUZE_M] |
| 02_ds_R_roter_riese | 1,7·10⁻⁸ | 1.66689e-08 | 1 | -3.31132e-10 | 5e-10 | standard | OK | [G: BETEIGEUZE_M, BETEIGEUZE_R] |
| 02_g_roter_riese | 0,0085 | 0.00845577 | m/s² | -4.42251e-05 | 5e-05 | standard | OK | [G: BETEIGEUZE_M, BETEIGEUZE_R] |
| 02_delta_1fm | 10⁻³⁹ | 1.24212e-39 | 1 | 2.42116e-40 | 5e-40 | standard | OK |  |
| 02_alpha_em_laengen | 137,036 | 137.036 | 1 | -9.15652e-07 | 0.0005 | standard | OK |  |
| 02_umlauf_wasserstoff | 6,58·10¹⁵ | 6.57968e+15 | Hz | -3.1608e+11 | 5e+12 | standard | OK |  |
| 02_bindung_wasserstoff_g | 9·10²¹ | 9.21938e+21 | g⊕ | 2.19384e+20 | 5e+20 | standard | OK |  |
| 02_steilheit_masse | 3,8·10¹² | 3.79461e+12 | kg | -5.38714e+09 | 5e+10 | standard | OK |  |
| 02_r_stern_proton | 32 | 31.857 | fm | -0.142961 | 0.5 | standard | OK |  |
| 02_bohr_zu_r_stern | 1660 | 1661.1 | 1 | 1.09978 | 5 | standard | OK |  |
| 02_glaettung_wasserstoff | 2,7 | 2.71921 | mm | 0.019212 | 0.05 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 02_saturn_faltung | 95 | 95.1614 | 1 | 0.161414 | 0.5 | standard | OK | [G: SATURN_M] |
| 02_saturn_kante | 14 | 13.9076 | % | -0.0924043 | 0.5 | standard | OK | [G: R_earth, SATURN_M, SATURN_R] |
| 02_v_erde | 29,78 | 29.7847 | km/s | 0.00469183 | 0.005 | standard | OK |  |
| 02_ds_sonnensystem | 492,9 | 492.869 | m | -0.0313491 | 0.05 | standard | OK | [G: JUPITER_M, MARS_M, MERKUR_M, NEPTUN_M, SATURN_M, URANUS_M, VENUS_M] |
| 02_anteil_sonne | 99,87 | 99.866 | % | -0.00397177 | 0.005 | standard | OK | [G: JUPITER_M, MARS_M, MERKUR_M, NEPTUN_M, SATURN_M, URANUS_M, VENUS_M] |
| 02_a0_schwelle | 1,1·10⁻¹⁰ | 1.1e-10 | m/s² | -1.29247e-26 | 5e-12 | standard | OK |  |
| 02_glaettung_sonne_l | 0,93 | 0.933531 | pc | 0.00353062 | 0.005 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 02_hill_sonne | 1,37 | 1.37219 | pc | 0.00218791 | 0.005 | Sonnenkreis_gemessen | OK |  |
| 02_ds_stern | 912 | 912.188 | m | 0.187683 | 0.5 | standard | OK |  |
| 02_m_stern | 1,853 | 1.85326 | M☉ | 0.000255212 | 0.0005 | standard | OK |  |
| 02_carter_grenze_mstern | 0,7 | 0.701468 | M_* | 0.00146842 | 0.05 | standard | OK | [G: CARTER_GRENZE_BEOB] |
| 02_carter | 2,944 | 2.94412 | 1 | 0.000120845 | 0.0005 | standard | OK |  |
| 02_carter_eins | 1402 | 1401.75 | 1 | -0.252206 | 0.5 | standard | OK |  |
| 02_eddington_zu_L | 3,28·10⁴ | 32838.7 | 1 | 38.6933 | 50 | standard | OK |  |
| 02_delta_sonne_alpha | 0,04 | 0.0398582 | α² | -0.000141768 | 0.005 | standard | OK |  |
| 02_ns_huelle_r0 | 2,9 | 2.90237 | r₀ | 0.0023714 | 0.05 | standard | OK |  |
| 02_ns_takt | 19 | 19.0399 | % | 0.0398766 | 0.5 | standard | OK |  |
| 02_f_konvektiv | 0,54 | 0.538494 | 1 | -0.00150633 | 0.005 | standard | OK |  |
| 02_f_strahlend | 0,85 | 0.854316 | 1 | 0.00431644 | 0.005 | standard | OK |  |
| 02_T_c_sonne_homogen | 1,2 | 1.21525 | 10⁷ K | 0.0152508 | 0.05 | standard | OK | [G: MU_SONNE_HOMOGEN] |
| 02_T_c_sonne_kern | 1,7·10⁷ | 1.66607e+07 | K | -339304 | 500000 | standard | OK | [G: MU_SONNE_KERN] |
| 02_delta_pp | 7·10⁻⁷ | 6.98622e-07 | 1 | -1.3779e-09 | 5e-08 | standard | OK | [G: MU_SONNE_HOMOGEN, T_ZUEND_PP] |
| 02_delta_3a | 8·10⁻⁶ | 8.12148e-06 | 1 | 1.21482e-07 | 5e-07 | standard | OK | [G: T_ZUEND_3A] |
| 02_delta_si | 1,7·10⁻⁴ | 0.000174032 | 1 | 4.03175e-06 | 5e-06 | standard | OK | [G: T_ZUEND_SI] |
| 02_decke_faktor | 0,1 | 0.0963768 | (μ/μ_e)·m_ec² | -0.0036232 | 0.05 | standard | OK |  |
| 02_decke_kohlenstoff | 5·10⁸ | 5.00066e+08 | K | 66382.4 | 5e+07 | standard | OK |  |
| 02_grenzmasse_zaehlregel | 1,456 | 1.4563 | M☉ | 0.000297553 | 0.0005 | standard | OK |  |
| 02_helium_zu_grenze | 0,41 | 0.414339 | M_Grenze | 0.00433862 | 0.005 | standard | OK | [G: T_ZUEND_3A] |
| 02_helium_konvektiv | 0,40 | 0.400933 | M☉ | 0.000933447 | 0.005 | standard | OK | [G: T_ZUEND_3A] |
| 02_helium_strahlend | 0,60 | 0.6034 | M☉ | 0.00340032 | 0.005 | standard | OK | [G: T_ZUEND_3A] |
| 02_alpha_G_proton | 5,906·10⁻³⁹ | 5.90615e-39 | 1 | 1.49414e-43 | 5e-43 | standard | OK |  |
| 02_mp_me | 1836,15 | 1836.15 | 1 | 0.00267344 | 0.005 | standard | OK |  |
| 02_im_tau_G | 14,0096 | 14.0096 | 1 | -1.40588e-05 | 5e-05 | standard | OK |  |
| 02_im_tau_alpha | 0,7831 | 0.783081 | 1 | -1.89006e-05 | 5e-05 | standard | OK |  |
| 02_alpha_s_planck | 53 | 53.0196 | 1 | 0.0195885 | 0.5 | standard | OK | [G: ALPHA_S_MZ, M_BOTTOM, M_CHARM, M_TOP] |
| 02_nome_qcd_planck | 1,2·10⁻²⁰ | 1.16722e-20 | 1 | -3.27777e-22 | 5e-22 | standard | OK | [G: ALPHA_S_MZ, M_BOTTOM, M_CHARM, M_TOP] |
| 02_alpha_G_groessenordnung | 5,9·10⁻³⁹ | 5.90615e-39 | 1 | 6.14941e-42 | 2.95e-39 | standard | OK |  |
| 02_empfindlichkeit_plus | 41 | 41.1067 | % | 0.106651 | 0.5 | standard | OK | [G: ALPHA_S_MZ, M_BOTTOM, M_CHARM, M_TOP] |
| 02_empfindlichkeit_minus | −30 | -29.805 | % | 0.195013 | 5 | standard | OK | [G: ALPHA_S_MZ, M_BOTTOM, M_CHARM, M_TOP] |
| 02_dirac_verhaeltnis | 68,5 | 68.518 | 1 | 0.0179995 | 0.05 | standard | OK |  |
| 02_ladungskreis_radius | 23,4 | 23.4125 | ℓ_P | 0.0124752 | 0.05 | standard | OK |  |
| 02_ladungskreis_umlauf | 5·10¹⁷ | 5.2147e+17 | GeV | 2.14699e+16 | 5e+16 | standard | OK |  |
| 02_zustand_eisen_erde | 116 | 116.009 | 1 | 0.00933819 | 0.5 | standard | OK | [G: R_earth] |
| 02_f_s_sonne | 0,0183 | 0.0182788 | 1 | -2.11637e-05 | 5e-05 | standard | OK |  |
| 02_gleichstand_druck_mstern | 72,6 | 72.606 | M_* | 0.00597582 | 0.05 | standard | OK | [G: MU_SONNE_HOMOGEN] |
| 02_gleichstand_druck_msonne | 134,6 | 134.557 | M☉ | -0.0425969 | 0.05 | standard | OK | [G: MU_SONNE_HOMOGEN] |
| 02_offen_01_1 | 29 |  |  |  |  | standard | OFFEN (nur Forschung: Verhältnis der eigenen Faltung zu Bindung und Phase je Element (Rechnung elemente_analyse, Notiz 75: trägt nicht)) |  |
| 02_offen_01_2 | 36 |  |  |  |  | standard | OFFEN (nur Forschung: Verhältnis der eigenen Faltung zu Bindung und Phase je Element (Rechnung elemente_analyse, Notiz 75: trägt nicht)) |  |
| 02_offen_08_1 | 10²¹ |  |  |  |  | standard | OFFEN (nur Forschung: Zustandsschwelle δ_Z = 2kT/(m_A c²) mit homogener Kugel, Median über PubChem (Notiz 75, Rechnung elemente_analyse)) |  |
| 02_offen_12_1 | 0,7 |  |  |  |  | standard | OFFEN (nur Forschung: Temperaturskala gegen die Kerntemperaturen der Hauptreihe [G] (Rechnung fusion_trichter)) |  |
| 02_offen_12_2 | 1,3 |  |  |  |  | standard | OFFEN (nur Forschung: Temperaturskala gegen die Kerntemperaturen der Hauptreihe [G] (Rechnung fusion_trichter)) |  |

### Baustein 03

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 03_v_flach_milchstrasse | 172,0 | 172.034 | km/s | 0.0338436 | 0.05 | standard | OK |  |
| 03a_ds_milchstrasse | 197,4 | 197.413 | AE | 0.0125743 | 0.05 | standard | OK |  |
| 03a_dsR_milchstrasse | 6,3·10⁻⁸ | 6.25545e-08 | 1 | -4.45545e-10 | 5e-10 | standard | OK |  |
| 03a_ds_haufen | 8,8 | 8.8445 | Lj | 0.0444987 | 0.05 | standard | OK |  |
| 03a_dsR_haufen | 1,4·10⁻⁶ | 1.35587e-06 | 1 | -4.41322e-08 | 5e-08 | standard | OK |  |
| 03a_masse_universum | 7,4·10²² | 7.40535e+22 | M☉ | 5.35283e+19 | 5e+20 | standard | OK | [G: OMEGA_R, Ω_b] |
| 03a_ds_universum | 3,9 | 3.85274 | Mrd Lj | -0.0472568 | 0.05 | standard | OK | [G: OMEGA_R, Ω_b] |
| 03a_dsR_universum | 0,083 | 0.0833401 | 1 | 0.000340077 | 0.0005 | standard | OK | [G: OMEGA_R, Ω_b] |
| 03a_straffung_rand | 6·10⁻⁸ | 6.25545e-08 | 1 | 2.55446e-09 | 5e-09 | standard | OK |  |
| 03a_deff_kugel | 2 | 2.00002 | 1 | 1.58113e-05 | 0.5 | standard | OK |  |
| 03a_deff_filament | 1,5 | 1.50001 | 1 | 7.90565e-06 | 0.05 | standard | OK |  |
| 03a_r_stern_filament | 2,5 | 2.53424 | kpc | 0.0342407 | 0.05 | standard | OK |  |
| 03a_filament_querschnitt | 0,46 | 0.459646 | H₀ | -0.000353718 | 0.005 | standard | OK |  |
| 03a_v_milchstrasse_5kpc | 250 | 250.094 | km/s | 0.0936394 | 5 | standard | OK |  |
| 03a_v_milchstrasse_15kpc | 198 | 197.55 | km/s | -0.449543 | 0.5 | standard | OK |  |
| 03a_v_milchstrasse_50kpc | 180 | 179.586 | km/s | -0.413925 | 5 | standard | OK |  |
| 03a_v_milchstrasse_rstern | 216 | 216.378 | km/s | 0.378415 | 0.5 | standard | OK |  |
| 03a_sigma_1e7 | 9,2 | 9.21446 | km/s | 0.0144586 | 0.05 | standard | OK |  |
| 03a_sigma_2e8 | 19,5 | 19.4862 | km/s | -0.0137925 | 0.05 | standard | OK |  |
| 03a_sparc_galaxien | 149 | 149 | 1 | 0 | 0.5 | standard | OK |  |
| 03a_upsilon_d | 0,5 | 0.5 | M☉/L☉[3,6] | 0 | 0.05 | standard | OK |  |
| 03a_innen_aussen | −0,02 | -0.0210543 | dex | -0.00105433 | 0.005 | standard | OK |  |
| 03a_innen_aussen_fehler | 0,013 | 0.0132524 | dex | 0.000252397 | 0.0005 | standard | OK |  |
| 03a_rang_sb_gasreich | 0,36 | 0.363482 | 1 | 0.00348234 | 0.005 | standard | OK |  |
| 03a_rang_sb_ge002 | 0,28 | 0.279503 | 1 | -0.000496823 | 0.005 | standard | OK |  |
| 03a_rang_sb_ge003 | 0,23 | 0.232019 | 1 | 0.00201877 | 0.005 | standard | OK |  |
| 03a_rang_sb_ge003_p | 0,07 | 0.0650623 | 1 | -0.0049377 | 0.005 | standard | OK |  |
| 03a_streuung_rest | 0,14 | 0.137392 | dex | -0.00260768 | 0.005 | standard | OK |  |
| 03a_streuung_regel | 0,05 | 0.0500255 | dex | 2.55066e-05 | 0.005 | standard | OK |  |
| 03a_nu_kugelsternhaufen | 1,007 | 1.00678 | 1 | -0.000216345 | 0.0005 | standard | OK |  |
| 03a_nu_zwerg | 3,3 | 3.30348 | 1 | 0.00348234 | 0.05 | standard | OK |  |
| 03a_grenzflaechendichte | 251 | 251.207 | M☉/pc² | 0.207308 | 0.5 | standard | OK |  |
| 03a_haufen_nah | 50 | 50.2025 | % | 0.20249 | 5 | standard | OK |  |
| 03a_haufen_fern | 77 | 77.1918 | % | 0.191828 | 0.5 | standard | OK |  |
| 03a_nu_frei_50 | 6,25 | 6.24889 | 1 | -0.00110855 | 0.005 | standard | OK |  |
| 03a_nu_senkrecht_50 | 5,36 | 5.35553 | 1 | -0.00446555 | 0.005 | standard | OK |  |
| 03a_nu_gleich_50 | 4,59 | 4.58904 | 1 | -0.000962801 | 0.005 | standard | OK |  |
| 03a_plummer_A | 2,37 | 2.37337 | g_e | 0.00336584 | 0.005 | standard | OK |  |
| 03a_nu_ge01 | 3,69 | 3.68859 | 1 | -0.00141384 | 0.005 | standard | OK |  |
| 03a_scheibe_A | 1,27 | 1.27296 | g_e | 0.00296303 | 0.005 | standard | OK |  |
| 03a_m1_scheibe_30 | 4,4 | 4.37189 | % | -0.0281131 | 0.05 | standard | OK |  |
| 03a_m1_scheibe_50 | 31 | 30.8341 | % | -0.165947 | 0.5 | standard | OK |  |
| 03a_m1_punkt_20_ge | 2 | 1.9846 | % | -0.0154007 | 0.5 | standard | OK |  |
| 03a_m1_punkt_50_ge | 36 | 36.1113 | % | 0.111335 | 0.5 | standard | OK |  |
| 03a_m1_punkt_20_nu | 13 | 12.9472 | % | -0.0527879 | 0.5 | standard | OK |  |
| 03a_m1_punkt_50_nu | 20 | 19.6374 | % | -0.36259 | 5 | standard | OK |  |
| 03a_theta_heute | 1,06 | 1.05546 | 1 | -0.00453825 | 0.005 | standard | OK | [G: OMEGA_R] |
| 03a_z_theta | 0,3 | 0.299999 | 1 | -7.26668e-07 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03a_lambda_horizont | 14 | 13.9737 | Gpc | -0.0263237 | 0.5 | horizont | OK |  |
| 03a_lambda_literatur | 12 | 12.1361 | Gpc | 0.136113 | 0.5 | standard | OK | [G: A0_LIT] |
| 03a_vakuum_atmend | 22 | 21.9325 | % | -0.0675458 | 0.5 | horizont | OK |  |
| 03a_delta_neff | 2,09 | 2.09283 | 1 | 0.00283449 | 0.005 | horizont | OK | [G: N_EFF_MESS, N_EFF_STD] |
| 03a_vakuum_kruemmung | 0,22 | 0.219325 | 1 | -0.000675458 | 0.005 | horizont | OK |  |
| 03a_a0_horizont_674 | 1,04 | 1.0422 | 10⁻¹⁰ m/s² | 0.00219788 | 0.005 | horizont | OK |  |
| 03a_a0_horizont_70 | 1,08·10⁻¹⁰ | 1.0824e-10 | m/s² | 2.40136e-13 | 5e-13 | standard | OK |  |
| 03a_unter_literatur_70 | 10 | 9.79989 | % | -0.200113 | 5 | standard | OK | [G: A0_LIT] |
| 03a_unter_literatur_674 | 13 | 13.1502 | % | 0.150177 | 0.5 | horizont | OK | [G: A0_LIT] |
| 03a_h0_literatur | 77,6 | 77.6052 | km/s/Mpc | 0.00522398 | 0.05 | standard | OK | [G: A0_LIT] |
| 03a_tf_vflach | 1,5·10⁻¹⁰ | 1.52646e-10 | m/s² | 2.64583e-12 | 5e-12 | standard | OK |  |
| 03a_tf_s_median | 0,23 | 0.233574 | 1 | 0.00357414 | 0.005 | standard | OK |  |
| 03a_tf_regel | 0,97 | 0.966388 | 10⁻¹⁰ m/s² | -0.00361206 | 0.005 | standard | OK |  |
| 03a_a0_je_galaxie | 1,09·10⁻¹⁰ | 1.08782e-10 | m/s² | -2.18305e-13 | 5e-13 | standard | OK |  |
| 03a_a0_streuung_gemessen | 0,28 | 0.277449 | dex | -0.00255133 | 0.005 | standard | OK |  |
| 03a_a0_streuung_unten | 0,27 | 0.266809 | dex | -0.00319101 | 0.005 | standard | OK |  |
| 03a_a0_streuung_oben | 0,30 | 0.29986 | dex | -0.000140404 | 0.005 | standard | OK |  |
| 03a_anstieg_schwelle | 0,06 | 0.0639821 | dex/Dekade | 0.00398214 | 0.005 | standard | OK |  |
| 03a_anstieg_schwelle_fehler | 0,04 | 0.0395778 | dex/Dekade | -0.000422249 | 0.005 | standard | OK |  |
| 03a_a0_zelle | 1,13·10⁻¹⁰ | 1.12705e-10 | m/s² | -2.95417e-13 | 5e-13 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R] |
| 03a_a0_standard | 1,1·10⁻¹⁰ | 1.1e-10 | m/s² | -1.29247e-26 | 5e-12 | standard | OK |  |
| 03a_rest_w_heute | 1,5·10⁻² | 0.0146107 | 1 | -0.000389258 | 0.0005 | mitatmend_heute | OK | [G: OMEGA_R] |
| 03a_rest_w_horizont | 1,3·10⁻² | 0.0130156 | 1 | 1.55885e-05 | 0.0005 | horizont | OK | [G: OMEGA_R] |
| 03a_rest_w_rekombination | 738 | 738.42 | 1 | 0.420039 | 0.5 | z1090 | OK | [G: OMEGA_R] |
| 03a_Q_kosmisch | 5,36 | 5.36433 | 1 | 0.00432722 | 0.005 | standard | OK |  |
| 03a_q_kosmisch | 0,843 | 0.842874 | 1 | -0.000125799 | 0.0005 | standard | OK |  |
| 03a_Q_radiant | 5,283 | 5.28319 | 1 | 0.000185307 | 0.0005 | standard | OK |  |
| 03a_Q_bao | 5,32 | 5.32248 | 1 | 0.00247993 | 0.005 | standard | OK |  |
| 03a_Q_planck_fehler | 0,065 | 0.0645867 | 1 | -0.000413297 | 0.0005 | standard | OK |  |
| 03a_Q_planck_abstand | 1,3 | 1.25633 | σ | -0.0436746 | 0.05 | standard | OK |  |
| 03a_schale_z_c_34 | 14 | 14.3893 | Mio. Jahre | 0.389347 | 0.5 | standard | OK |  |
| 03a_schale_z_c_2 | 45 | 44.6416 | Mio. Jahre | -0.358366 | 0.5 | standard | OK |  |
| 03a_schale_z_c_1 | 143 | 142.765 | Mio. Jahre | -0.235197 | 0.5 | standard | OK |  |
| 03a_schale_z_c_05 | 309 | 309.456 | Mio. Jahre | 0.455972 | 0.5 | standard | OK |  |
| 03a_nachlauf_heute_min | 0,28 | 0.275257 | Mrd. Jahre | -0.00474328 | 0.005 | standard | OK | [G: OMEGA_R] |
| 03a_nachlauf_heute_max | 3,5 | 3.48175 | Mrd. Jahre | -0.0182473 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03a_nachlauf_z2_min | 0,1 | 0.0909687 | Mrd. Jahre | -0.0090313 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03a_nachlauf_z2_max | 1,2 | 1.15067 | Mrd. Jahre | -0.0493268 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03a_aN_RQ | 0,029 | 0.0292196 | a₀ | 0.000219649 | 0.0005 | standard | OK |  |
| 03a_RQ_milchstrasse | 51 | 51.009 | kpc | 0.00902876 | 0.5 | standard | OK |  |
| 03a_chi_RQ | 0,917 | 0.916965 | 1 | -3.49941e-05 | 0.0005 | standard | OK |  |
| 03a_Rta_milchstrasse | 0,58 | 0.580149 | Mpc | 0.000148581 | 0.005 | standard | OK |  |
| 03a_nu_Rta | 67 | 67.0369 | 1 | 0.0369447 | 0.5 | standard | OK |  |
| 03a_Rta_m_b_klein | 319 | 319.268 | kpc | 0.268068 | 0.5 | standard | OK |  |
| 03a_Rta_m_milchstrasse | 580 | 580.149 | kpc | 0.148581 | 5 | standard | OK |  |
| 03a_Rta_m_b_gross | 1482 | 1481.91 | kpc | -0.0889004 | 0.5 | standard | OK |  |
| 03a_w_mitbewegt | 0,985 | 0.985389 | 1 | 0.000389258 | 0.0005 | standard | OK | [G: OMEGA_R] |
| 03a_masse_bQ | 4,45·10¹³ | 4.45204e+13 | M☉ | 2.0396e+10 | 5e+10 | standard | OK | [G: OMEGA_R] |
| 03_zelle_mindestens | 27,4 | 27.3691 | Gpc | -0.0308891 | 0.05 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R] |
| 03_zelle_fenster_unten | 27 | 27.3691 | Gpc | 0.369111 | 0.5 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R] |
| 03_rueckstell_baryonen | 0,083 | 0.0833401 | 1 | 0.000340077 | 0.0005 | standard | OK | [G: OMEGA_R, Ω_b] |
| 03_rueckstell_ruecklaeufe | 0,53 | 0.530404 | 1 | 0.000403522 | 0.005 | standard | OK | [G: OMEGA_R, Ω_b] |
| 03_omega_s_theta1 | 0,22 | 0.219325 | 1 | -0.000675458 | 0.005 | horizont | OK |  |
| 03_omega_s_standard | 0,24 | 0.244327 | 1 | 0.00432743 | 0.005 | standard | OK |  |
| 03_omega_s_zu_rest | 16,7 | 16.7225 | 1 | 0.0224516 | 0.05 | standard | OK |  |
| 03_theta_heute | 1,06 | 1.05546 | 1 | -0.00453825 | 0.005 | standard | OK |  |
| 03_a0_fest_seit_z10 | 0,3 | 0.319778 | % | 0.0197784 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03_a0_rekombination | 0,76 | 0.758622 | 1 | -0.001378 | 0.005 | z1090 | OK | [G: OMEGA_R] |
| 03_omega_s_rekombination | 2,6·10⁻¹⁰ | 2.62105e-10 | 1 | 2.10517e-12 | 5e-12 | z1090 | OK | [G: OMEGA_R] |
| 03_zelle_a0_unten | 28,1 | 28.1123 | Gpc | 0.0123252 | 0.05 | standard | OK |  |
| 03_zelle_a0_oben | 26,9 | 26.8777 | Gpc | -0.0222922 | 0.05 | standard | OK |  |
| 03_zelle_a0_standard | 27,6 | 27.5916 | Gpc | -0.00839356 | 0.05 | standard | OK |  |
| 03_d_lss_heute | 27,8 | 27.7859 | Gpc | -0.0141006 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03_zelle_beobachter_z1 | 1,32 | 1.32454 | 1 | 0.00454335 | 0.005 | standard | OK | [G: OMEGA_R] |
| 03_a0_max_wuerfel | 1,13·10⁻¹⁰ | 1.12705e-10 | m/s² | -2.95417e-13 | 5e-13 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_fluchtmass_verhaeltnis | 284 | 284.094 | 1 | 0.0944065 | 0.5 | standard | OK |  |
| 03_fluchtmass_ln | 5,65 | 5.64931 | 1 | -0.0006934 | 0.005 | standard | OK |  |
| 03_entropie_verhaeltnis | 0,24 | 0.244327 | 1 | 0.00432743 | 0.005 | standard | OK |  |
| 03_schliesslaenge_standard | 13,2 | 13.2394 | Gpc | 0.0393962 | 0.05 | standard | OK |  |
| 03_mindestrate | 3,8·10⁻³⁷ | 3.82867e-37 | s⁻¹ | 2.86687e-39 | 5e-39 | standard | OK | [G: ENTROPIEPRODUKTION] |
| 03_mindestrate_jahr | 10⁻²⁹ | 1.20824e-29 | 1/Jahr | 2.08236e-30 | 5e-30 | standard | OK | [G: ENTROPIEPRODUKTION] |
| 03_lambda_zuwachs | 1,6·10⁻¹⁰ | 1.56411e-10 | m/s | -3.58928e-12 | 5e-12 | standard | OK | [G: ENTROPIEPRODUKTION] |
| 03_ws_heute | 1·10⁻¹⁹ | 1.16855e-19 | 1 | 1.6855e-20 | 5e-20 | standard | OK | [G: ENTROPIEPRODUKTION] |
| 03_pulsar_regel | 6,3·10⁻¹⁹ | 6.33683e-19 | s⁻¹ | 3.68327e-21 | 5e-21 | standard | OK |  |
| 03_pulsar_newton | 4,1·10⁻¹⁹ | 4.14871e-19 | s⁻¹ | 4.87129e-21 | 5e-21 | standard | OK |  |
| 03_rueckkopplung_zeit | 1 | 1.01906 | Mio. Jahre | 0.0190579 | 0.5 | standard | OK |  |
| 03_haufenstoss_unten | 0,9 | 0.917522 | Mrd. Jahre | 0.0175224 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03_haufenstoss_oben | 1,8 | 1.81341 | Mrd. Jahre | 0.0134129 | 0.05 | standard | OK | [G: OMEGA_R] |
| 03_expscheibe_hoch | 0,08 | 0.0829354 | dex | 0.0029354 | 0.005 | standard | OK | [G: LELLI_SIGMA_DYN] |
| 03_expscheibe_tief | 0,22 | 0.224632 | dex | 0.00463208 | 0.005 | standard | OK | [G: LELLI_SIGMA_DYN] |
| 03_sparc_auswahl | 153 | 153 | Galaxien | 0 | 0.5 | standard | OK |  |
| 03_q_abstand_fuenf | 5,6 | 5.6409 | σ | 0.0409014 | 0.05 | standard | OK |  |
| 03_q_max_wuerfel | 5,5 | 5.52081 | 1 | 0.0208077 | 0.05 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_rstern_zu_lambda | 0,97 | 0.970956 | 1 | 0.000955609 | 0.005 | standard | OK | [G: Ω_b] |
| 03_zellgrenze_wuerfel | 7,17·10⁻¹⁰ | 7.17289e-10 | m/s² | 2.88844e-13 | 5e-13 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_q_fehler_planck | 1,2 | 1.204 | % | 0.00400378 | 0.05 | standard | OK |  |
| 03_warp_viertel | 0,14 | 0.141428 | Mrd. Jahre | 0.0014284 | 0.005 | standard | OK | [G: WARP_PRAEZESSION] |
| 03_aussen_zahl | 74 | 74 | Galaxien | 0 | 0.5 | standard | OK |  |
| 03_aussen_kleinstes | 0,008 | 0.00756433 | a₀ | -0.000435669 | 0.0005 | standard | OK |  |
| 03_druck_regel_pal14 | 1,8 | 1.77681 | km/s | -0.0231876 | 0.05 | standard | OK | [G: PAL14_ABSTAND, PAL14_M, PAL14_RE] |
| 03_druck_regel_ngc2419 | 6,5 | 6.48635 | km/s | -0.0136519 | 0.05 | standard | OK | [G: NGC2419_ABSTAND, NGC2419_M, NGC2419_RE] |
| 03_druck_regel_fornax | 14,0 | 13.9648 | km/s | -0.0351984 | 0.05 | standard | OK | [G: FORNAX_ABSTAND, FORNAX_M, FORNAX_RE] |
| 03_druck_regel_df2 | 17,7 | 17.7023 | km/s | 0.00226007 | 0.05 | standard | OK | [G: DF2_ABSTAND, DF2_M, DF2_RE] |
| 03_druck_regel_df4 | 19,2 | 19.2128 | km/s | 0.012805 | 0.05 | standard | OK | [G: DF4_ABSTAND, DF4_M, DF4_RE] |
| 03_druck_regel_crater2 | 1,4 | 1.38699 | km/s | -0.0130093 | 0.05 | standard | OK | [G: CRATER2_ABSTAND, CRATER2_M, CRATER2_RE] |
| 03_druck_regel_draco | 3,3 | 3.33247 | km/s | 0.0324699 | 0.05 | standard | OK | [G: DRACO_ABSTAND, DRACO_M, DRACO_RE] |
| 03_druck_verstaerkung_df2 | 6,3 | 6.31319 | 1 | 0.0131916 | 0.05 | standard | OK | [G: DF2_ABSTAND, DF2_M, DF2_RE] |
| 03_druck_verstaerkung_df4 | 7,2 | 7.21123 | 1 | 0.0112263 | 0.05 | standard | OK | [G: DF4_ABSTAND, DF4_M, DF4_RE] |
| 03_druck_fornax_verlangt | 4,5 | 4.45009 | 1 | -0.0499086 | 0.05 | standard | OK | [G: FORNAX_ABSTAND, FORNAX_M, FORNAX_RE] |
| 03_druck_newton_geglaettet_fornax | 5,5 | 5.54627 | km/s | 0.046275 | 0.05 | standard | OK | [G: FORNAX_ABSTAND, FORNAX_M, FORNAX_RE] |
| 03_druck_newton_geglaettet_draco | 1,2 | 1.19707 | km/s | -0.00293307 | 0.05 | standard | OK | [G: DRACO_ABSTAND, DRACO_M, DRACO_RE] |
| 03_druck_newton_pal14 | 0,665 | 0.664916 | km/s | -8.35222e-05 | 0.0005 | standard | OK | [G: PAL14_ABSTAND, PAL14_M, PAL14_RE] |
| 03_druck_newton_ngc2419 | 5,09 | 5.08563 | km/s | -0.00436686 | 0.005 | standard | OK | [G: NGC2419_ABSTAND, NGC2419_M, NGC2419_RE] |
| 03_druck_newton_df4 | 7,15 | 7.15461 | km/s | 0.00461408 | 0.005 | standard | OK | [G: DF4_ABSTAND, DF4_M, DF4_RE] |
| 03_druck_draco_isoliert | 4,74 | 4.73786 | km/s | -0.00213699 | 0.005 | standard | OK | [G: DRACO_ABSTAND, DRACO_M, DRACO_RE] |
| 03_gezeit_quadrupol | 5,6·10⁻³⁰ | 5.55232e-30 | s⁻² | -4.76825e-32 | 5e-32 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 03_sparc_a0_profil_wert | 1,06 | 1.05828 | 10⁻¹⁰ m/s² | -0.00172107 | 0.005 | standard | OK |  |
| 03_sparc_a0_profil_stat | 0,085 | 0.0849193 | 10⁻¹⁰ m/s² | -8.07304e-05 | 0.0005 | standard | OK |  |
| 03_sparc_a0_profil_syst | 0,11 | 0.107476 | 10⁻¹⁰ m/s² | -0.0025238 | 0.005 | standard | OK |  |
| 03_sparc_a0_varianten_unten | 0,96 | 0.958157 | 10⁻¹⁰ m/s² | -0.00184333 | 0.005 | standard | OK |  |
| 03_sparc_a0_varianten_oben | 1,19 | 1.19267 | 10⁻¹⁰ m/s² | 0.00267115 | 0.005 | standard | OK |  |
| 03_sparc_a0_drittel_unten | 0,75 | 0.749921 | 10⁻¹⁰ m/s² | -7.86738e-05 | 0.005 | standard | OK |  |
| 03_sparc_a0_drittel_mitte | 1,51 | 1.51262 | 10⁻¹⁰ m/s² | 0.00262324 | 0.005 | standard | OK |  |
| 03_sparc_a0_drittel_oben | 1,06·10⁻¹⁰ | 1.0611e-10 | m/s² | 1.09926e-13 | 5e-13 | standard | OK |  |
| 03_sparc_a0_drittel_anstieg | 0,09 | 0.0904895 | 1 | 0.000489455 | 0.005 | standard | OK |  |
| 03_sparc_a0_drittel_fehler | 0,08 | 0.0782649 | 1 | -0.00173511 | 0.005 | standard | OK |  |
| 03_aussenfeld_senkrecht_gew | 7 | 7.40285 | % | 0.402851 | 0.5 | standard | OK |  |
| 03_aussenfeld_senkrecht_ungew | 14 | 14.3505 | % | 0.35048 | 0.5 | standard | OK |  |
| 03_aussenfeld_gleich_gew | 22 | 21.7189 | % | -0.28108 | 0.5 | standard | OK |  |
| 03_aussenfeld_gleich_ungew | 31 | 31.114 | % | 0.113974 | 0.5 | standard | OK |  |
| 03_sigma_dyn_rest_mittel | 0,05 | 0.0543969 | dex | 0.00439693 | 0.005 | standard | OK |  |
| 03_sigma_dyn_rest_ups04 | 0,01 | 0.0135651 | dex | 0.00356508 | 0.005 | Upsilon_0.4 | OK |  |
| 03_sigma_dyn_rest_ups06 | 0,09 | 0.0893901 | dex | -0.000609924 | 0.005 | Upsilon_0.6 | OK |  |
| 03_sigma_dyn_rest_anstieg | 0,05 | 0.0483444 | 1 | -0.00165564 | 0.005 | standard | OK |  |
| 03_sigma_dyn_rest_anstieg_fehler | 0,025 | 0.0249724 | 1 | -2.75998e-05 | 0.0005 | standard | OK |  |
| 03_rest_gasanteil_unten | −0,015 | -0.0154844 | dex | -0.000484397 | 0.0005 | Upsilon_0.4 | OK |  |
| 03_rest_gasanteil_oben | 0,018 | 0.0181817 | dex | 0.000181688 | 0.0005 | Upsilon_0.6 | OK |  |
| 03_rest_gasanteil_fehler | 0,013 | 0.0128427 | dex | -0.000157302 | 0.0005 | Upsilon_0.4 | OK |  |
| 03_rest_jenseits_rq_wert | −0,02 | -0.0169178 | dex | 0.00308219 | 0.005 | standard | OK |  |
| 03_rest_jenseits_rq_fehler | 0,05 | 0.0517572 | dex | 0.00175721 | 0.005 | standard | OK |  |
| 03_reststeigung_wert | −0,30 | -0.298252 | 1 | 0.001748 | 0.005 | standard | OK |  |
| 03_reststeigung_fehler | 0,07 | 0.0674292 | 1 | -0.00257078 | 0.005 | standard | OK |  |
| 03_sattel_abstand | 2,59·10⁸ | 2.58813e+08 | m | -186768 | 500000 | standard | OK |  |
| 03_sattel_feld_galaxie | 1,34 | 1.34162 | a₀ | 0.00162158 | 0.005 | Sonnenkreis_gemessen | OK | [G: V0_SONNE] |
| 03_sattel_ungeglaettet | 7·10⁻¹¹ | 7.12371e-11 | m/s² | 1.23713e-12 | 5e-12 | standard | OK |  |
| 03_warp_aufwicklung_regel | 0,26 | 0.255592 | Mrd. Jahre | -0.00440818 | 0.005 | standard | OK | [G: MCMILLAN_BULGE, MCMILLAN_DICK, MCMILLAN_DUENN, MCMILLAN_H2, MCMILLAN_HI] |
| 03_warp_aufwicklung_newton | 0,28 | 0.279937 | Mrd. Jahre | -6.29281e-05 | 0.005 | standard | OK | [G: MCMILLAN_BULGE, MCMILLAN_DICK, MCMILLAN_DUENN, MCMILLAN_H2, MCMILLAN_HI] |
| 03_druck_df2_gleichgerichtet | 17 | 16.9748 | km/s | -0.0251795 | 0.5 | standard | OK | [G: DF2_ABSTAND, DF2_M, DF2_RE] |
| 03_druck_df2_ohne_feld | 21 | 20.7056 | km/s | -0.294361 | 0.5 | standard | OK | [G: DF2_ABSTAND, DF2_M, DF2_RE] |
| 03_gegenlaeufig_drehsinn | 0,53 | 0.532418 | 1 | 0.00241778 | 0.005 | standard | OK |  |
| 03_drehung_moeglich_60 | 60 | 60 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_moeglich_90 | 90 | 90 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_moeglich_120 | 120 | 120 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_moeglich_180 | 180 | 180 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_jeder | 180 | 180 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_quadrat | 90 | 90 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_sechseck_60 | 60 | 60 | ° | 0 | 5 | standard | OK |  |
| 03_drehung_sechseck_120 | 120 | 120 | ° | 0 | 5 | standard | OK |  |
| 03_quadrupol_kosmische_streuung | 63 | 63.2456 | % | 0.245553 | 0.5 | standard | OK |  |
| 03_quadrupol_isotrop_achsen | 40 | 40 | % | 0 | 5 | standard | OK |  |
| 03_quadrupol_isotrop_m0 | 20 | 20 | % | 0 | 5 | standard | OK |  |
| 03_a0_max_sechseckig | 1,30 | 1.3014 | 10⁻¹⁰ m/s² | 0.00140042 | 0.005 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_a0_max_raumzentriert | 1,46 | 1.46408 | 10⁻¹⁰ m/s² | 0.00407548 | 0.005 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_a0_max_flächenzentriert | 1,59 | 1.59388 | 10⁻¹⁰ m/s² | 0.00388349 | 0.005 | standard | OK | [G: KREISPAAR_GRENZE, OMEGA_R, Ω_b] |
| 03_formfaktor_sechseckig | 1,15 | 1.1547 | 1 | 0.00470054 | 0.005 | standard | OK |  |
| 03_formfaktor_raumzentriert | 1,30 | 1.29904 | 1 | -0.000961894 | 0.005 | standard | OK |  |
| 03_formfaktor_flächenzentriert | 1,41 | 1.41421 | 1 | 0.00421356 | 0.005 | standard | OK |  |
| 03_I_F_keim_klein | 1,1 | 1.0955 | rad | -0.0045007 | 0.05 | standard | OK |  |
| 03_I_F_keim_gross | 1,25 | 1.25436 | rad | 0.00436497 | 0.005 | standard | OK |  |
| 03_I_F_flach | 1,75 | 1.749 | rad | -0.000998904 | 0.005 | standard | OK |  |
| 03_schale_rand | 1,16 | 1.15749 | R_ta | -0.00251286 | 0.005 | standard | OK | [G: OMEGA_R] |
| 03_I_vor_umkehr | 0,1 | 0.109766 | rad | 0.00976618 | 0.05 | standard | OK |  |
| 03_w_aussen_133 | 0,95 | 0.952113 | 1 | 0.00211295 | 0.005 | standard | OK |  |
| 03_fsigma8_uebergabe | 1,2·10⁻³ | 0.00118472 | 1 | -1.5276e-05 | 5e-05 | mitatmend_heute | OK | [G: OMEGA_R] |
| 03_fsigma8_sofort | 5·10⁻³ | 0.00492375 | 1 | -7.6248e-05 | 0.0005 | mitatmend_heute | OK | [G: OMEGA_R] |
| 03_stoerung_aN_rekombination | 2,6 | 2.59048 | a₀(z) | -0.00952051 | 0.05 | z1090 | OK | [G: BAO_SKALA, OMEGA_R] |
| 03_horizontwachstum_heute | 6,6·10¹⁷ | 6.5553e+17 | 1 | -4.46997e+15 | 5e+15 | mitatmend_heute | OK | [G: ENTROPIEPRODUKTION, OMEGA_R] |
| 03_ws_zerstrahlung | −1,0002 | -1.00019 | 1 | 5.62983e-06 | 5e-05 | mitatmend_heute | OK | [G: OMEGA_R] |
| 03_linse_milchstrasse_1mpc | 0,12 | 0.12155 | M☉/pc² | 0.0015497 | 0.005 | standard | OK |  |
| 03_haufen_rest_100kpc | 3,1 | 3.13155 | 1 | 0.0315484 | 0.05 | standard | OK | [G: NFW_C, NFW_M200, NFW_Z, OMEGA_R] |
| 03_haufen_rest_1mpc | 1,6 | 1.61354 | 1 | 0.0135415 | 0.05 | standard | OK | [G: NFW_C, NFW_M200, NFW_Z, OMEGA_R] |
| 03_offen_01_1 | 10⁻¹¹ |  |  |  |  | standard | OFFEN (fehlt Formel: Neigung der großen Skalen (Filamente, Wände, Leerräume)) |  |
| 03_offen_01_2 | 10⁻¹² |  |  |  |  | standard | OFFEN (fehlt Formel: Neigung der großen Skalen (Filamente, Wände, Leerräume)) |  |
| 03_offen_02_1 | 2·10⁻⁶ |  |  |  |  | standard | OFFEN (fehlt Formel: Straffung S(<r)/r im Zentrum der Milchstraße (Massenprofil mit Kern)) |  |
| 03_offen_13_1 | 0,5 |  |  |  |  | standard | OFFEN (fehlt Formel: a₀-Lesung an SPARC je Υ_d gegen die Grenze a₀ ≤ 1,13·10⁻¹⁰ der Zelle) |  |
| 03_offen_14_1 | 0,5 |  |  |  |  | standard | OFFEN (fehlt Formel: a₀-Lesung an SPARC je Υ_d mit äußerem Feld (Senkung bei Υ_d = 0,5)) |  |
| 03_offen_17_1 | 0,05 |  |  |  |  | standard | OFFEN (fehlt Formel: σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an) |  |
| 03_offen_18_1 | 0,7 |  |  |  |  | standard | OFFEN (fehlt Formel: σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an (0,7σ, ± 0,05)) |  |
| 03_offen_27_1 | 1,05 |  |  |  |  | standard | OFFEN (fehlt Formel: Linse der Schale außerhalb R_ta ((1 + w·Q)·M_b ab 1,05·R_ta)) |  |
| 03_offen_30_1 | 50 |  |  |  |  | standard | OFFEN (fehlt Formel: Radiant der Bahn seit der Umkehr bei 100 kpc (Rechnung haufen_nach_naht)) |  |
| 03_offen_30_2 | 100 |  |  |  |  | standard | OFFEN (fehlt Formel: Radiant der Bahn seit der Umkehr bei 100 kpc (Rechnung haufen_nach_naht)) |  |
| 03_offen_31_1 | 0,5 |  |  |  |  | standard | OFFEN (fehlt Formel: f = 1/S_T je Ebene) |  |
| 03_offen_32_1 | 0,85 |  |  |  |  | standard | OFFEN (fehlt Formel: f = 1/S auf allen Ebenen (0,5–0,85)) |  |
| 03_offen_34_1 | 20 |  |  |  |  | standard | OFFEN (nur Forschung: Umkehr-z je Masse (Rechnung kollaps_nach_umkehr, noch nicht im Modell)) |  |
| 03_offen_34_2 | 10⁸ |  |  |  |  | standard | OFFEN (nur Forschung: Umkehr-z je Masse (Rechnung kollaps_nach_umkehr, noch nicht im Modell)) |  |
| 03_offen_37_1 | 42 |  |  |  |  | standard | OFFEN (nur Forschung: Grenze der Prüfbarkeit beim 1,5-fachen Durchmesser der letzten Streufläche (Rechnung zelle_form, 1,5·d_LSS)) |  |
| 03_offen_38_1 | 17 |  |  |  |  | standard | OFFEN (nur Forschung: Senkung des Quadrupols bei L = 0,985·d_LSS (Würfel) aus forschung.winkelmoden_zelle mit d_LSS im eigenen Atmen (forschung.lauf_eigenes_atmen S2); frisch 17,0 % (bei 1,5·d_LSS 2,4 %)) |  |
| 03_offen_40_1 | 1,5 |  |  |  |  | standard | OFFEN (nur Forschung: Grenze der Prüfbarkeit beim 1,5-fachen Durchmesser der letzten Streufläche (Rechnung zelle_form)) |  |
| 03_offen_41_1 | 17 |  |  |  |  | standard | OFFEN (nur Forschung: Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0) |  |
| 03_offen_41_2 | 30 |  |  |  |  | standard | OFFEN (nur Forschung: Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0) |  |
| 03_offen_41_3 | 33 |  |  |  |  | standard | OFFEN (nur Forschung: Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0) |  |
| 03_offen_41_4 | 76 |  |  |  |  | standard | OFFEN (nur Forschung: Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0) |  |
| 03_offen_41_5 | 38 |  |  |  |  | standard | OFFEN (nur Forschung: Winkelmoden je Zellform aus forschung.winkelmoden_zelle (a₀ = 1,06 aus der SPARC-Lesung, Q = 5,36, χ = d_LSS/2 im eigenen Atmen): frisch Quadrupol 15,0–17,0 %, große Winkel 29,7–32,8 %, Würfel 75,8 % in E_g, sechseckig 38,2 % in m = 0) |  |
| 03_offen_42_1 | 1,01·10⁻¹⁰ |  |  |  |  | standard | OFFEN (nur Forschung: Einsatz der Kreispaare zellgrenze(f³, L_min = d_LSS)/(1 + Q) mit d_LSS im eigenen Atmen (forschung.lauf_eigenes_atmen S2) und Q = 5,36: frisch Würfel 1,013, sechseckig 1,170·10⁻¹⁰) |  |
| 03_offen_42_2 | 1,17 |  |  |  |  | standard | OFFEN (nur Forschung: Einsatz der Kreispaare zellgrenze(f³, L_min = d_LSS)/(1 + Q) mit d_LSS im eigenen Atmen (forschung.lauf_eigenes_atmen S2) und Q = 5,36: frisch Würfel 1,013, sechseckig 1,170·10⁻¹⁰) |  |
| 03_offen_52_1 | 3·10⁻⁵ |  |  |  |  | standard | OFFEN (fehlt Formel: Streuung der Beschleunigungsrelation bei früh kollabierten Galaxien) |  |
| 03_offen_55_1 | 1,9 |  |  |  |  | standard | OFFEN (fehlt Formel: Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)) |  |
| 03_offen_55_2 | 0,013 |  |  |  |  | standard | OFFEN (fehlt Formel: Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)) |  |
| 03_offen_55_3 | 0,05 |  |  |  |  | standard | OFFEN (fehlt Formel: Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)) |  |
| 03_offen_55_4 | 0,004 |  |  |  |  | standard | OFFEN (fehlt Formel: Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)) |  |
| 03_offen_55_5 | 50 |  |  |  |  | standard | OFFEN (fehlt Formel: Preisach-Verlauf frisch kollabierter Galaxien (Abweichung von der Relation je Kollaps- und Beobachtungsepoche)) |  |
| 03_offen_56_1 | 0,3 |  |  |  |  | standard | OFFEN (nur Forschung: Kollapszeit der Mode gegen die Kugel (Rechnung kollaps_nach_umkehr, noch nicht im Modell)) |  |
| 03_offen_56_2 | 0,96 |  |  |  |  | standard | OFFEN (nur Forschung: Kollapszeit der Mode gegen die Kugel (Rechnung kollaps_nach_umkehr, noch nicht im Modell)) |  |
| 03_offen_57_1 | 0,54 |  |  |  |  | standard | OFFEN (fehlt Formel: Endantwort der Mode auf eine Massenänderung δM/M) |  |
| 03_offen_57_2 | 0,59 |  |  |  |  | standard | OFFEN (fehlt Formel: Endantwort der Mode auf eine Massenänderung δM/M) |  |
| 03_offen_60_1 | 19 |  |  |  |  | standard | OFFEN (nur Forschung: Linse um isolierte Galaxien mit Restmode (Rechnung linsen_isoliert, noch nicht im Modell)) |  |
| 03_offen_60_2 | 39 |  |  |  |  | standard | OFFEN (nur Forschung: Linse um isolierte Galaxien mit Restmode (Rechnung linsen_isoliert, noch nicht im Modell)) |  |
| 03_offen_61_1 | 0,09 |  |  |  |  | standard | OFFEN (fehlt Formel: Skala, unter der das Strukturwachstum von der Regel abweicht, je z (Umkehrskala)) |  |
| 03_offen_61_2 | 0,15 |  |  |  |  | standard | OFFEN (fehlt Formel: Skala, unter der das Strukturwachstum von der Regel abweicht, je z (Umkehrskala)) |  |
| 03_offen_61_3 | 0,27 |  |  |  |  | standard | OFFEN (fehlt Formel: Skala, unter der das Strukturwachstum von der Regel abweicht, je z (Umkehrskala)) |  |
| 03_offen_70_1 | 0,15 |  |  |  |  | standard | OFFEN (fehlt Formel: nötige Genauigkeit von Ω_c/Ω_b, um Q = 2π − 1 von der Zählung je Radiant zu trennen) |  |
| 03_offen_70_2 | 0,5 |  |  |  |  | standard | OFFEN (fehlt Formel: nötige Genauigkeit von Ω_c/Ω_b, um Q = 2π − 1 von der Zählung je Radiant zu trennen) |  |
| 03_offen_71_1 | 1,0 |  |  |  |  | standard | OFFEN (fehlt Formel: σ_Q aus (ω_b, ω_c) mit BAO; Q_fehler nimmt keine Werte an) |  |
| 03_offen_77_1 | 100 |  |  |  |  | standard | OFFEN (nur Forschung: Linse um isolierte Galaxien bei kleinem äußerem Feld (0,001·a₀): Radius, ab dem sie keplerisch wird (Rechnung linsen_isoliert, noch nicht im Modell)) |  |
| 03_offen_77_2 | 200 |  |  |  |  | standard | OFFEN (nur Forschung: Linse um isolierte Galaxien bei kleinem äußerem Feld (0,001·a₀): Radius, ab dem sie keplerisch wird (Rechnung linsen_isoliert, noch nicht im Modell)) |  |
| 03_offen_78_1 | 0,35 |  |  |  |  | standard | OFFEN (fehlt Formel: Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017)) |  |
| 03_offen_78_2 | 0,68 |  |  |  |  | standard | OFFEN (fehlt Formel: Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017)) |  |
| 03_offen_78_3 | 0,025 |  |  |  |  | standard | OFFEN (fehlt Formel: Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017)) |  |
| 03_offen_78_4 | 0,02 |  |  |  |  | standard | OFFEN (fehlt Formel: Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017)) |  |
| 03_offen_79_1 | 10⁶ |  |  |  |  | standard | OFFEN (fehlt Formel: Drehmoment eines Störers auf einen geneigten Ring und Massenmodell der Milchstraße (McMillan 2017): Verwindung) |  |
| 03_offen_80_1 | 90 |  |  |  |  | standard | OFFEN (nur Forschung: Knotenlinie bei einfallendem Gas (Rechnung warp, noch nicht im Modell)) |  |
| 03_offen_83_1 | 0,012 |  |  |  |  | standard | OFFEN (nur Forschung: Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)) |  |
| 03_offen_84_1 | 2,35 |  |  |  |  | standard | OFFEN (nur Forschung: Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)) |  |
| 03_offen_84_2 | 1,44 |  |  |  |  | standard | OFFEN (nur Forschung: Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)) |  |
| 03_offen_84_3 | 10,2 |  |  |  |  | standard | OFFEN (nur Forschung: Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)) |  |
| 03_offen_84_4 | 179 |  |  |  |  | standard | OFFEN (nur Forschung: Gedächtnis der Umgebung und Außenfeld der Zwerge (Rechnungen gedaechtnis_umgebung, zwerge_aussenfeld: trägt nicht)) |  |

### Baustein 04

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 04_grenzmasse_T_CMB | 4,5·10²² | 4.50159e+22 | kg | 1.59484e+19 | 5e+20 | standard | OK |  |
| 04_takt_100r0 | 0,995 | 0.994987 | 1 | -1.25629e-05 | 0.0005 | standard | OK |  |
| 04_takt_2r0 | 0,707 | 0.707107 | 1 | 0.000106781 | 0.0005 | standard | OK |  |
| 04_uhr_erde | 0,022 | 0.0219681 | s | -3.18881e-05 | 0.0005 | standard | OK | [G: R_earth] |
| 04_uhr_sonne | 67 | 66.9812 | s | -0.0188418 | 0.5 | standard | OK |  |
| 04_uhr_neutronenstern | 69,5 | 69.5431 | d | 0.0431494 | 0.05 | standard | OK |  |
| 04_endrate_H_s | 33 | 33.3155 | km/s/Mpc | 0.315475 | 0.5 | standard | OK |  |
| 04_ds_H_heute | 0,74 | 0.741327 | Gpc | 0.00132655 | 0.005 | standard | OK |  |
| 04_ds_H_zu_Lambda | 0,056 | 0.055994 | 1 | -6.01505e-06 | 0.0005 | standard | OK |  |
| 04_delta_weisser_zwerg | 1,0·10⁻⁴ | 0.000101888 | 1 | 1.88824e-06 | 5e-06 | standard | OK |  |
| 04_entropieexport_erde | 6,1·10¹⁴ | 6.08258e+14 | W/K | -1.74181e+12 | 5e+12 | standard | OK | [G: ALBEDO_ERDE, R_earth] |
| 04_entropie_falte_36 | 36 | 36 | 1 | -1.42109e-14 | 0.5 | standard | OK |  |
| 04_page_zeit | 0,54 | 0.538101 | 1 | -0.00189911 | 0.005 | standard | OK | [G: BETA_PAGE] |
| 04_fuellung_horizont | 10⁻¹⁸ | 1.36839e-18 | 1 | 3.68387e-19 | 5e-19 | standard | OK | [G: S_FALTEN_HEUTE] |
| 04_rueckverdichtung_faktor | 1,5 | 1.4847 | 1 | -0.0153 | 0.05 | standard | OK | [G: BETA_PAGE] |
| 04_chsh_zwei_harmonische | 2,36 | 2.36256 | 1 | 0.00256023 | 0.005 | standard | OK |  |
| 04_rueckgabe_nenner | 15360 | 15360 | 1 | -1.81899e-12 | 5 | standard | OK |  |
| 04_rueckgabe_klein | 4·10¹¹ | 3.612e+11 | Jahr | -3.87996e+10 | 5e+10 | mitatmend_heute | OK | [G: OMEGA_R] |
| 04_rueckgabe_gross | 10¹² | 1.13599e+12 | Jahr | 1.35993e+11 | 5e+11 | mitatmend_heute | OK | [G: OMEGA_R] |
| 04_kreislauf_versiegen | 2,6·10¹¹ | 2.55384e+11 | Jahr | -4.61567e+09 | 5e+09 | mitatmend_heute | OK | [G: ENTROPIEPRODUKTION, OMEGA_R] |

### Baustein 05

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 05_q_selbsterregt | −0,05 | -0.05 | 1 | 1.38778e-17 | 0.005 | standard | OK |  |
| 05_resolvente_invariant | 10⁴ | 10000 | 1 | 1.10049e-09 | 1e-05 | standard | OK |  |
| 05_resolvente_flucht | 567 | 567.455 | 1 | 0.455009 | 0.5 | standard | OK |  |
| 05_kanal_verhaeltnis | 1,44·10⁹ | 1.43652e+09 | 1 | -3.48146e+06 | 5e+06 | standard | OK | [G: R_earth] |
| 05_gedaechtnis_umkehr | 9,75 | 9.74747 | 1 | -0.00253442 | 0.005 | standard | OK |  |
| 05_nu_tau_halb | 1,045 | 1.04517 | 1 | 0.000165705 | 0.0005 | standard | OK |  |
| 05_nu_tau_zwei | 1,0000035 | 1 | 1 | -1.26455e-08 | 5e-08 | standard | OK |  |
| 05_delta_koeffizient_3 | 252 | 252 | 1 | 0 | 0.5 | standard | OK |  |
| 05_normierung_E2 | 24 | 24 | 1 | 0 | 0.5 | standard | OK |  |
| 05_normierung_E4 | 240 | 240 | 1 | -5.68434e-14 | 5 | standard | OK |  |
| 05_normierung_E6 | 504 | 504 | 1 | -1.13687e-13 | 0.5 | standard | OK |  |
| 05_zahl_1728 | 1728 | 1728 | 1 | -4.54747e-13 | 0.5 | standard | OK |  |
| 05_vakuumteil_nenner | 15 | 15 | 1 | -3.55271e-15 | 0.5 | standard | OK |  |

### Baustein 06

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 06_siede_schmelz_min | 1,02 | 1.02274 | 1 | 0.00274417 | 0.005 | standard | OK |  |
| 06_siede_schmelz_max | 8,2 | 8.17735 | 1 | -0.0226536 | 0.05 | standard | OK |  |
| 06_rate_quer | −0,5 | -0.5 | 1 | 0 | 0.05 | standard | OK |  |
| 06_rate_rastpunkt | −0,098 | -0.0979796 | 1 | 2.04103e-05 | 0.0005 | standard | OK | [W: TORUS_DOMEGA, TORUS_K] |
| 06_dreieck_kleinster_wert | 8 | 8 | 1 | 0 | 0.5 | standard | OK |  |
| 06_spaltbarkeit_max | 0,93 | 0.92742 | 1 | -0.0025805 | 0.005 | standard | OK |  |
| 06_geisterbahn_min | 0,0125 |  | 1 |  |  | standard | OFFEN (fehlt Formel: Verschiebung des Rastfensters beim Durchfahren nahe dem verschwundenen Sattel-Knoten (Rampe durch die Adler-Gleichung)) |  |
| 06_geisterbahn_max | 0,015 |  | 1 |  |  | standard | OFFEN (fehlt Formel: Verschiebung des Rastfensters beim Durchfahren nahe dem verschwundenen Sattel-Knoten (Rampe durch die Adler-Gleichung)) |  |
| 06_s_labor | 3·10⁵ | 298633 | 1 | -1366.75 | 50000 | standard | OK |  |
| 06_gezeit_labor | 10⁻²⁶ |  | m/s² |  |  | standard | OFFEN (fehlt Formel: Änderung des galaktischen Zusatzes über 1 km, 𝒯·ℓ·(ν + dν/d ln g) am Sonnenort; gezeit_sonnenort steht im Kern, ν + dν/d ln g nicht, Messwerte des Sonnenorts (Oort A, B, v₀, R₀, ρ_b) kommen aus Teil 03) |  |
| 06_rc_nu | 1,1567 | 1.15667 | 1 | -3.15981e-05 | 5e-05 | standard | OK |  |
| 06_rc_wurzel | −0,499 | -0.499374 | 1 | -0.000373784 | 0.0005 | standard | OK |  |
| 06_rc_offenes_ende | 9,0 | 9.03314 | % | 0.0331411 | 0.05 | standard | OK |  |
| 06_vo2_hopf_spannung | 3,55 | 3.54822 | V | -0.00178261 | 0.005 | standard | OK |  |
| 06_vo2_hopf_frequenz | 3,26 | 3.25633 | MHz | -0.00367353 | 0.005 | standard | OK |  |
| 06_vo2_torus_spannung | 3,47 |  | V |  |  | standard | OFFEN (fehlt Formel: VO₂-Mikrobrücke im Magnetfeld (Leitfähigkeitssprung, Hall, Wärme; Hopf und Neimark-Sacker)) |  |
| 06_einrasten_min | 0,57 | 0.56875 | 1 | -0.00125 | 0.005 | standard | OK |  |
| 06_einrasten_max | 0,63 | 0.629688 | 1 | -0.0003125 | 0.005 | standard | OK |  |
| 06_kupfer_zeitverhaeltnis | 0,74 | 0.737672 | 1 | -0.00232817 | 0.005 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |
| 06_kupfer_guete_min | 3000 | 3000 | 1 | 0 | 500 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |
| 06_kupfer_abstand_min | 0,0010 | 0.00103247 | Hz | 3.24731e-05 | 5e-05 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |
| 06_kupfer_abstand_max | 0,6 | 0.623837 | Hz | 0.0238374 | 0.05 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |
| 06_kupfer_ring_min | 0,20 | 0.202767 | Hz | 0.0027675 | 0.005 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |
| 06_kupfer_ring_max | 1,5 | 1.50461 | Hz | 0.00460804 | 0.05 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] |

### Baustein 07

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 07_maximum_d3 | 7,96 | 7.96052 | a₀ | 0.000520131 | 0.005 | standard | OK |  |
| 07_maximum_d2_ort | 2,54 | 2.53964 | a₀ | -0.000361718 | 0.005 | standard | OK |  |
| 07_maximum_d2_wert | 0,648 | 0.64761 | a₀ | -0.000389762 | 0.0005 | standard | OK |  |
| 07_maximum_d2_punktmasse | 0,628 | 0.6275 | r* | -0.000499513 | 0.0005 | standard | OK |  |
| 07_maximum_d2_milchstrasse | 5,5 | 5.4714 | kpc | -0.028598 | 0.05 | standard | OK |  |
| 07_stufen_bei_rstern | 0,582 | 0.581977 | 1 | -2.32931e-05 | 0.0005 | standard | OK |  |
| 07_leiter_nenner_15 | 15 | 15 | 1 | -3.55271e-15 | 0.5 | standard | OK |  |
| 07_ds_r_erdbahn | 3,3·10⁻⁹ | 3.29021e-09 | 1 | -9.79043e-12 | 5e-11 | standard | OK |  |
| 07_ds_r_milchstrasse | 1,2·10⁻⁷ | 1.19635e-07 | 1 | -3.64604e-10 | 5e-09 | Sonnenkreis_Punktmasse | OK |  |
| 07_a0_standard | 1,1·10⁻¹⁰ | 1.1e-10 | m/s² | -1.29247e-26 | 5e-12 | standard | OK |  |
| 07_rate_a0_c | 3,7·10⁻¹⁹ | 3.66921e-19 | s⁻¹ | -3.0795e-21 | 5e-21 | standard | OK |  |
| 07_tau_D_mrd | 86 | 86.3623 | Mrd. Jahre | 0.36227 | 0.5 | standard | OK |  |
| 07_periode_sonnenoberflaeche | 2,5·10⁵ | 251280 | 1 | 1279.81 | 5000 | standard | OK |  |
| 07_ds_sparc | 1,03 | 1.03065 | 1 | 0.000648111 | 0.005 | standard | OK |  |
| 07_punkte_unter_s_trenn | 16 | 16 | Punkte | 0 | 0.5 | standard | OK |  |
| 07_coth_schwelle | 1,313 | 1.31304 | 1 | 3.52855e-05 | 0.0005 | standard | OK |  |
| 07_beide_sinne_faktor | 1,37 | 1.36788 | 1 | -0.00212056 | 0.005 | standard | OK |  |
| 07_milchstrasse_rstern_beide | 253 | 253.068 | km/s | 0.0683572 | 0.5 | standard | OK |  |
| 07_milchstrasse_rstern | 216 | 216.378 | km/s | 0.378415 | 0.5 | standard | OK |  |
| 07_vorgang_s1 | 19 | 18.6724 | % | -0.327625 | 0.5 | standard | OK | [G: WELTALTER] |
| 07_vorgang_s03 | 10 | 10.2761 | % | 0.276051 | 5 | standard | OK | [G: WELTALTER] |
| 07_theta_tau_D | 17,9 | 17.8591 | 1 | -0.0409389 | 0.05 | mitatmend_heute | OK |  |
| 07_wurzel_summe | 1,30·10⁻² | 0.0130014 | 1 | 1.35035e-06 | 5e-05 | mitatmend_heute | OK |  |
| 07_wurzel_produkt | 5,4·10⁻³ | 0.00537499 | 1 | -2.50083e-05 | 5e-05 | mitatmend_heute | OK |  |
| 07_q_milchstrasse_8kpc | 0,33624 | 0.336243 | 1 | 3.4381e-06 | 5e-06 | Sonnenkreis_Punktmasse | OK |  |
| 07_gewicht_rueckkehr_heute | 2·10⁻⁴ | 0.000216639 | 1 | 1.66391e-05 | 5e-05 | mitatmend_heute | OK |  |
| 07_rand_065 | 0,66 | 0.6569 | r*·Λ/ϖ_R | -0.00309973 | 0.005 | standard | OK |  |
| 07_flach_bis_N | 1,5 | 1.5223 | N | 0.0223011 | 0.05 | standard | OK |  |
| 07_abweichung_eigenwert | 9,0 | 9.03314 | % | 0.0331411 | 0.05 | standard | OK |  |
| 07_rand_0657 | 0,657 | 0.6569 | Λ | -9.97281e-05 | 0.0005 | standard | OK |  |
| 07_im_tau_erdbahn | 3,1 | 3.10866 | 1 | 0.00866441 | 0.05 | standard | OK |  |
| 07_im_tau_milchstrasse | 2,5 | 2.53674 | 1 | 0.0367415 | 0.05 | Sonnenkreis_Punktmasse | OK |  |
| 07_masse_220 | 1,60·10¹¹ | 1.60467e+11 | M☉ | 4.67248e+08 | 5e+08 | standard | OK |  |
| 07_rosette_winkel | 254,6 | 254.558 | ° | -0.0415587 | 0.05 | standard | OK |  |
| 07_resonanz_rstern_0_45 | 0,45 | 0.445039 | r* | -0.00496144 | 0.005 | standard | OK |  |
| 07_resonanz_rstern_0_68 | 0,68 | 0.682141 | r* | 0.00214062 | 0.005 | standard | OK |  |
| 07_resonanz_rstern_0_94 | 0,94 | 0.944949 | r* | 0.00494895 | 0.005 | standard | OK |  |
| 07_resonanz_rstern_2_1 | 2,1 | 2.0695 | r* | -0.0305021 | 0.05 | standard | OK |  |
| 07_resonanz_rstern_12 | 12 | 12.3311 | r* | 0.331069 | 0.5 | standard | OK |  |
| 07_resonanz_mw_5_9 | 5,9 | 5.94783 | kpc | 0.0478289 | 0.05 | standard | OK |  |
| 07_resonanz_mw_8_2 | 8,2 | 8.23935 | kpc | 0.0393491 | 0.05 | standard | OK |  |
| 07_resonanz_mw_10_9 | 10,9 | 10.9098 | kpc | 0.00983687 | 0.05 | standard | OK |  |
| 07_resonanz_mw_18_0 | 18,0 | 18.0447 | kpc | 0.0446948 | 0.05 | standard | OK |  |
| 07_resonanz_mw_108 | 108 | 107.519 | kpc | -0.480987 | 0.5 | standard | OK |  |
| 07_kern_sonnenkreis | 2·10⁻¹⁵ | 2.08697e-15 | 1 | 8.69678e-17 | 5e-16 | Sonnenkreis_Punktmasse | OK |  |
| 07_s_erdbahn | 7342 | 7342.33 | 1 | 0.332747 | 0.5 | standard | OK |  |
| 07_ring_12_2 | 1,42 | 1.41602 | 1 | -0.00398298 | 0.005 | standard | OK |  |
| 07_ring_12_89 | 1,66 | 1.66171 | 1 | 0.00171185 | 0.005 | standard | OK |  |
| 07_ring_20_2 | 1,69 | 1.69347 | 1 | 0.00346595 | 0.005 | standard | OK |  |
| 07_ring_20_89 | 1,91 | 1.91114 | 1 | 0.00113822 | 0.005 | standard | OK |  |
| 07_ring_30_2 | 2,00 | 1.9966 | 1 | -0.00339794 | 0.005 | standard | OK |  |
| 07_ring_30_89 | 2,19 | 2.18874 | 1 | -0.00126131 | 0.005 | standard | OK |  |
| 07_ell_erde | 2780 | 2778.19 | AE | -1.80698 | 5 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_ell_mond | 641 | 641.195 | AE | 0.194603 | 0.5 | Sonnenkreis_gemessen | OK | [G: M_MOND, RHO_B_SONNE, V0_SONNE] |
| 07_umlauf_an_ell | 84,5 | 84.4965 | Mio. Jahre | -0.00349318 | 0.05 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_schnittmasse | 3,3·10⁸ | 3.25325e+08 | M☉ | -4.6747e+06 | 5e+06 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_schnittradius | 642 | 642.049 | pc | 0.0485176 | 0.5 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_ell_zu_rstern_sonne | 26,2 | 26.2253 | 1 | 0.0252502 | 0.05 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_watson_drei | 1,516 | 1.51639 | 1 | 0.000386059 | 0.0005 | standard | OK |  |
| 07_drehung_golden_55 | 55 | 55 | 1 | 0 | 0.5 | standard | OK |  |
| 07_drehung_golden_89 | 89 | 89 | 1 | 0 | 0.5 | standard | OK |  |
| 07_drehung_golden_144 | 144 | 144 | 1 | 0 | 0.5 | standard | OK |  |
| 07_delta_gegen_reihe | 4,7·10⁻⁸ | 4.6558e-08 | 1 | -4.42045e-10 | 5e-10 | standard | OK |  |
| 07_feldform_viertel_hoch | 0,82 | 0.818343 | 1 | -0.0016573 | 0.005 | standard | OK |  |
| 07_feldform_viertel_tief | 0,53 | 0.529217 | 1 | -0.000782831 | 0.005 | standard | OK |  |
| 07_feldform_rd_tief | 0,79 | 0.792446 | 1 | 0.00244639 | 0.005 | standard | OK |  |
| 07_feldform_rd_hoch | 0,92 | 0.922096 | 1 | 0.00209624 | 0.005 | standard | OK |  |
| 07_feldform_gleichstand | 3,4 | 3.39641 | R_d | -0.0035886 | 0.05 | standard | OK |  |
| 07_feldform_naeherung_dex | 0,02 | 0.0194167 | dex | -0.000583262 | 0.005 | standard | OK |  |
| 07_gewicht_f_unten | 0,03 | 0.0254718 | 1 | -0.00452816 | 0.005 | Upsilon_0.4 | OK |  |
| 07_gewicht_f_oben | 0,13 | 0.132854 | 1 | 0.00285383 | 0.005 | Upsilon_0.6 | OK |  |
| 07_gewicht_f_ohne_ic2574 | 0,00 | 0.00400123 | 1 | 0.00400123 | 0.005 | standard | OK |  |
| 07_gewicht_f_ohne_fehler | 0,03 | 0.0305325 | 1 | 0.000532491 | 0.005 | standard | OK |  |
| 07_dynamikmasse_kern | 57 | 57.3285 | M☉/pc² | 0.328456 | 0.5 | standard | OK | [G: SIGMA_BAR_SONNE] |
| 07_dynamikmasse_ort_schicht | 69 | 69.1827 | M☉/pc² | 0.182665 | 0.5 | standard | OK | [G: SIGMA_BAR_SONNE] |
| 07_dynamikmasse_ort_ebene | 70 | 69.9924 | M☉/pc² | -0.00760034 | 5 | standard | OK | [G: SIGMA_BAR_SONNE] |
| 07_dynamikmasse_feldform | 92 | 92.3243 | M☉/pc² | 0.324301 | 0.5 | standard | OK | [G: MCMILLAN_BULGE, MCMILLAN_DICK, MCMILLAN_DUENN, MCMILLAN_H2, MCMILLAN_HI] |
| 07_ringe_faktor_innen | 2,3 | 2.32341 | 1 | 0.0234098 | 0.05 | standard | OK |  |
| 07_ringe_faktor_aussen | 3,0 | 3.03897 | 1 | 0.0389746 | 0.05 | standard | OK |  |
| 07_hill_sonne_newton | 188 | 188.03 | km/s | 0.0297874 | 0.5 | Sonnenkreis_gemessen | OK |  |
| 07_hill_sonne_regel | 229 | 228.608 | km/s | -0.39248 | 0.5 | Sonnenkreis_gemessen | OK |  |
| 07_quadrupol_glaettung | 8,6·10⁻³¹ | 8.58355e-31 | s⁻² | -1.64488e-33 | 5e-33 | Sonnenkreis_gemessen | OK | [G: RHO_B_SONNE, V0_SONNE] |
| 07_ds_fehler_ziehung | 0,07 | 0.074248 | 1 | 0.00424801 | 0.005 | standard | OK |  |
| 07_boden_grenze | 0,03 | 0.0315556 | a₀ | 0.00155562 | 0.005 | standard | OK |  |
| 07_boden_peclet | 0,016 | 0.0157778 | 1 | -0.000222192 | 0.0005 | standard | OK |  |
| 07_boden_bandbreite | 25 | 25.285 | ϖ_R | 0.285026 | 0.5 | standard | OK |  |
| 07_offen_01_1 | 36 |  |  |  |  | standard | OFFEN (fehlt Formel: χ²/N 36 statt 16 an SPARC (10 Galaxien, Υ_d frei) der linearen nichtlokalen Fortsetzung (Laplace, ¦k¦, v⁴ = (2/π)²GMa₀); benötigt: Regel der nichtlokalen Fortsetzung an Rotationskurven und ihre Anpassung) |  |
| 07_offen_02_1 | 16 |  |  |  |  | standard | OFFEN (fehlt Formel: χ²/N 16 der örtlichen Regel (10 Galaxien, Υ_d frei); benötigt: Anpassung je Galaxie mit freiem Υ_d. Zugleich ausnahme: Koeffizient in c⁴R/16πG; „16 Punkte“ unter s ≈ 0,1 ist Prüfzahl 07_punkte_unter_s_trenn) |  |
| 07_offen_05_1 | 0,13 |  |  |  |  | standard | OFFEN (fehlt Formel: d ln v/d ln r = −0,13 bis −0,27 für andere Verlustraten (Bahnfrequenz, √(a_N a₀)/c, Neigung mit Rückkopplung), benötigt: Regel mit diesen Raten; zudem λ_m = 0,01–0,13 Ω am Sonnenkreis (benötigt: Umlauffrequenz des gebundenen Bereichs, Parametersatz fehlt) und f = +0,13 an SPARC (benötigt: Anpassung von f)) |  |
| 07_offen_06_1 | 0,27 |  |  |  |  | standard | OFFEN (fehlt Formel: d ln v/d ln r = −0,27 für andere Verlustraten; benötigt: Regel mit diesen Raten) |  |
| 07_offen_07_1 | 1,18 |  |  |  |  | standard | OFFEN (fehlt Formel: ν ≥ 1,18 überall bei ϖ-abhängigen Koeffizienten; benötigt: Rückkehr mit ϖ-abhängigem Γ/D) |  |
| 07_offen_10_1 | 650 |  |  |  |  | standard | OFFEN (fehlt Formel: Δχ² ≈ 650 für drei Dimensionen (Δχ² ≥ 7 für zwei); benötigt: Gitter-Rückkehr ν_d(a_N) in d = 2, 3 und χ² der tiefen Bins mit Bootstrap-Kovarianz) |  |
| 07_offen_11_1 | 0,1 |  |  |  |  | standard | OFFEN (fehlt Formel: Wurzel von der vollen Bose-Form erst unter s ≈ 0,1 zu trennen; benötigt: Trennschärfe beider Formen an den SPARC-Punkten) |  |
| 07_offen_12_1 | 0,05 |  |  |  |  | standard | OFFEN (fehlt Formel: Band ¦λ¦ ≲ 0,05 (d ln v/d ln r = √(¼ + λ) − ½) mit nur λ = 0 für SL(2,ℤ); benötigt: Steigungsband aus SPARC und Spektrum des Laplace auf SL(2,ℤ)) |  |
| 07_offen_13_1 | 10¹²⁴ |  |  |  |  | standard | OFFEN (fehlt Formel: Normierung der Rückkehr als freie Energie eines Feldes) |  |
| 07_offen_20_1 | 6·10⁻⁶ |  |  |  |  | standard | OFFEN (nur Forschung: Wellenform der Rückkehr gäbe w ∝ √M, Milchstraße 6·10⁻⁶ (Rechnung gleichung_zeit, Welle trägt nicht; keine Kernformel)) |  |
| 07_offen_21_1 | 10⁶ |  |  |  |  | standard | OFFEN (fehlt Formel: Rückkehranteil q unabhängig von der Schrittlänge h (Irrfahrt bis 10⁶·r₀)) |  |
| 07_offen_22_1 | 3·10³² |  |  |  |  | standard | OFFEN (fehlt Formel: Rückkehrzeit zur Kugel r* im Hubble-Volumen (Sonne 3·10³² Jahre); benötigt: Rückkehrzeit im Raumvolumen ⟨T⟩ ∝ 1/M) |  |
| 07_offen_24_1 | 30 |  |  |  |  | standard | OFFEN (fehlt Formel: 30 % Masseverlust eines Zwergs (Bahnen beider Formen nach 3 Mrd. Jahren auf 5·10⁻⁶ gleich) und Umlaufarbeit durch die Schicht 11–30 %; benötigt: zeitabhängige Bahnen beider Wahlen des Imaginärteils, Feldform-Löser der Scheibe. „bei 30 kpc“ ist Eingabe RING_R30) |  |
| 07_offen_25_1 | 5·10⁻⁶ |  |  |  |  | standard | OFFEN (fehlt Formel: Bahnen beider Formen nach einem Massensprung auf 5·10⁻⁶ gleich; benötigt: zeitabhängige Bahnen) |  |
| 07_offen_26_1 | 0,01 |  |  |  |  | standard | OFFEN (fehlt Formel: λ_m = 0,01–0,13 Ω am Sonnenkreis; benötigt: Umlauffrequenz ω_orb des gebundenen Bereichs der Milchstraße (Dichtefaktor fehlt im Parametersatz)) |  |
| 07_offen_27_1 | 0,13 |  |  |  |  | standard | OFFEN (fehlt Formel: d ln v/d ln r = −0,13 bis −0,27 für andere Verlustraten (Bahnfrequenz, √(a_N a₀)/c, Neigung mit Rückkopplung), benötigt: Regel mit diesen Raten; zudem λ_m = 0,01–0,13 Ω am Sonnenkreis (benötigt: Umlauffrequenz des gebundenen Bereichs, Parametersatz fehlt) und f = +0,13 an SPARC (benötigt: Anpassung von f)) |  |
| 07_offen_28_1 | 10⁻⁹ |  |  |  |  | standard | OFFEN (fehlt Formel: Spannung des Wirbelstroms der Feldform, als Energiedichte gelesen) |  |
| 07_offen_29_1 | 30 |  |  |  |  | standard | OFFEN (fehlt Formel: 30 % Masseverlust eines Zwergs (Bahnen beider Formen nach 3 Mrd. Jahren auf 5·10⁻⁶ gleich) und Umlaufarbeit durch die Schicht 11–30 %; benötigt: zeitabhängige Bahnen beider Wahlen des Imaginärteils, Feldform-Löser der Scheibe. „bei 30 kpc“ ist Eingabe RING_R30) |  |
| 07_offen_30_1 | 5·10¹⁰ |  |  |  |  | standard | OFFEN (fehlt Formel: Scheibe 5·10¹⁰ M☉ (Bulge 10¹⁰ M☉) der Rechnung des Energieaustauschs über Sternbahnen; benötigt: Wirbelteil der Feldform und Bahnrechnung) |  |
| 07_offen_31_1 | 10¹⁰ |  |  |  |  | standard | OFFEN (fehlt Formel: Wirbelteil der Feldform und Bahnrechnung damit (Scheibe 5·10¹⁰ M☉, Bulge 10¹⁰ M☉)) |  |
| 07_offen_32_1 | 0,00 |  |  |  |  | standard | OFFEN (fehlt Formel: Energieänderung 0,00 ± 0,03 % von σ_z² und f = 0,00 ± 0,03 ohne IC 2574; benötigt: Bahnrechnung in der Feldform, Anpassung von f) |  |
| 07_offen_33_1 | 0,03 |  |  |  |  | standard | OFFEN (fehlt Formel: Boden unter 0,03·a₀ (95 %) an SPARC, f = +0,03, f = 0,00 ± 0,03, Energieänderung 0,00 ± 0,03 % von σ_z²; benötigt: Anpassung a = ν(a_N)·a_N + B, Anpassung von f, Bahnrechnung in der Feldform) |  |
| 07_offen_34_1 | 0,6 |  |  |  |  | standard | OFFEN (fehlt Formel: Υ_d = 0,6 derselben Vergleiche; zudem „je Stern um höchstens 0,6 %“ (Energieaustausch über Sternbahnen in der Feldform, benötigt: Wirbelteil der Feldform und Bahnrechnung damit)) |  |
| 07_offen_35_1 | 1,32 |  |  |  |  | standard | OFFEN (fehlt Formel: Überschussradius der Polytrope n = 3: 1,32·GM/c²; benötigt: δR = ∫(√g_rr − 1)dr über ein Lane-Emden-Profil) |  |
| 07_offen_36_1 | 10⁻⁵ |  |  |  |  | standard | OFFEN (fehlt Formel: Wachstum fσ₈ mit der Summe der Verluste (1 − w)(ν − 1)) |  |
| 07_offen_37_1 | 10⁻⁴ |  |  |  |  | standard | OFFEN (fehlt Formel: Rückkehr mit dem Feld der Störung allein (a_N ≈ 10⁻⁴·a₀)) |  |

### Baustein 08

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 08_ablenkung_sonnenrand | 1,7512 | 1.75119 | ″ | -9.6742e-06 | 5e-05 | standard | OK |  |
| 08_ablenkung_takt_haelfte | 0,8756 | 0.875595 | ″ | -4.8371e-06 | 5e-05 | standard | OK |  |
| 08_takt_iss | 0,021 | 0.0206703 | s | -0.000329665 | 0.0005 | standard | OK | [G: R_earth] [W: ISS_HOEHE] |
| 08_takt_erdbahn | 0,31 | 0.311493 | s | 0.00149335 | 0.005 | standard | OK |  |
| 08_rR_proton | 38 | 37.8619 | 1 | -0.138056 | 0.5 | standard | OK |  |
| 08_rR_erde | 3·10⁵ | 298789 | 1 | -1210.77 | 50000 | standard | OK | [G: R_earth] |
| 08_rR_sonne | 1,6·10⁶ | 1.57884e+06 | 1 | -21162.4 | 50000 | standard | OK |  |
| 08_rR_milchstrasse | 0,6 | 0.581291 | 1 | -0.0187094 | 0.05 | standard | OK | [W: R_MILCHSTRASSE_HUELLE] |
| 08_rR_haufen | 0,5 | 0.464123 | 1 | -0.0358767 | 0.05 | standard | OK | [W: R_HAUFEN_HUELLE] |
| 08_leiterzahl | 2,0·10⁻⁶² | 1.97816e-62 | 1 | -2.18411e-64 | 5e-64 | standard | OK |  |
| 08_schnittmasse_mitte | 5,9·10⁻²⁹ | 5.88617e-29 | kg | -1.38324e-31 | 5e-31 | standard | OK |  |
| 08_schnittmasse_oben | 1,1·10⁵⁴ | 1.10023e+54 | kg | 2.32331e+50 | 5e+52 | standard | OK |  |
| 08_rstern_compton_proton | 151 | 151.477 | 1 | 0.477364 | 0.5 | standard | OK |  |
| 08_vakuum_planck | 1,0·10⁻¹²² | 1.01646e-122 | 1 | 1.64606e-124 | 5e-124 | standard | OK |  |
| 08_rstern_erde | 12,7 | 12.7247 | AE | 0.0246877 | 0.05 | standard | OK |  |
| 08_rstern_sonne_ae | 7 340 | 7342.33 | AE | 2.33275 | 5 | standard | OK |  |
| 08_rstern_sonne_lj | 0,116 | 0.116101 | Lj | 0.000100691 | 0.0005 | standard | OK |  |
| 08_rstern_milchstrasse | 8,72 | 8.71936 | kpc | -0.000641094 | 0.005 | standard | OK |  |
| 08_verstaerkung_rstern | 1,582 | 1.58198 | 1 | -2.32931e-05 | 0.0005 | standard | OK |  |
| 08_linse_I1 | 2,135 | 2.13544 | 1 | 0.000443936 | 0.0005 | standard | OK |  |
| 08_alpha_inf_milchstrasse | 0,43 | 0.426768 | ″ | -0.00323246 | 0.005 | standard | OK |  |
| 08_alpha_inf_haufen | 23 | 22.7164 | ″ | -0.283563 | 0.5 | standard | OK |  |
| 08_aJ_erde | 3,27 | 3.27214 | m | 0.00213858 | 0.005 | standard | OK | [G: ERDE_J] |
| 08_ds_erde | 1,48 | 1.47834 | mm | -0.00165732 | 0.005 | standard | OK |  |
| 08_umlauf_differenz_erde | 137 | 137.158 | ns | 0.157907 | 0.5 | standard | OK | [G: ERDE_J] |
| 08_perihel_merkur | 43 | 42.9802 | ″/Jahrhundert | -0.0197919 | 0.5 | standard | OK | [G: MERKUR_A, MERKUR_E, MERKUR_P] |
| 08_isco_mitlaeufig | 2,32 | 2.32088 | 3Δs | 0.000883042 | 0.005 | standard | OK |  |
| 08_isco_gegenlaeufig | 8,72 | 8.71735 | 3Δs | -0.00264772 | 0.005 | standard | OK |  |
| 08_knotenzahl_milchstrasse | 5·10⁻⁷ | 5.02267e-07 | 1 | 2.26688e-09 | 5e-08 | standard | OK | [W: J_MILCHSTRASSE] |
| 08_gpb_geodaetisch | 6604 | 6603.88 | mas/Jahr | -0.123676 | 0.5 | standard | OK | [W: GPB_A] |
| 08_gpb_mitfuehrung | 40,8 | 40.8003 | mas/Jahr | 0.000345673 | 0.05 | standard | OK | [G: ERDE_J] [W: GPB_A] |
| 08_gpb_leitstern | 39,1 | 39.0507 | mas/Jahr | -0.0492759 | 0.05 | standard | OK | [G: ERDE_J] [W: GPB_A] |
| 08_hulse_taylor | −2,40·10⁻¹² | -2.40257e-12 | 1 | -2.56867e-15 | 5e-15 | standard | OK | [G: PSR1913_E, PSR1913_M1, PSR1913_M2, PSR1913_PB] |
| 08_ringdown_detektor | 249,9 | 249.922 | Hz | 0.0222482 | 0.05 | standard | OK | [G: GW150914_F] |
| 08_schatten_stosszahl | 15,6 | 15.5885 | Δs | -0.0115427 | 0.05 | standard | OK |  |
| 08_schatten_m87 | 40 | 39.6881 | μas | -0.311858 | 5 | standard | OK |  |
| 08_faltungsmesser_gw150914 | 15 | 14.8647 | km | -0.135308 | 0.5 | standard | OK | [G: GW150914_MC_DET] |
| 08_faltungsmesser_gw170817 | 590 | 589.518 | m | -0.482064 | 5 | standard | OK |  |
| 08_schwellenfrequenz_sonne | 10⁻¹³ | 1.00732e-13 | Hz | 7.31832e-16 | 1e-14 | standard | OK |  |
| 08_doppelpulsar_masse | 2,83 | 2.8284 | M☉ | -0.0016 | 0.005 | standard | OK | [G: PSR1913_M1, PSR1913_M2] |
| 08_doppelpulsar_ds | 1,39 | 1.39216 | km | 0.00216209 | 0.005 | standard | OK | [G: PSR1913_M1, PSR1913_M2] |
| 08_doppelpulsar_fe | 11,86 | 11.8568 | 1 | -0.00318427 | 0.005 | standard | OK | [G: PSR1913_E, PSR1913_M1, PSR1913_M2, PSR1913_PB] |
| 08_ringdown_abklingzeit | 4,03 | 4.03348 | ms | 0.00347675 | 0.005 | standard | OK | [G: GW150914_F, GW150914_Q] |
| 08_gedaechtnis_gw150914 | 2,7·10⁻²² | 2.70347e-22 | 1 | 3.46816e-25 | 5e-24 | standard | OK |  |
| 08_ringdown_quelle | 272,4 | 272.415 | Hz | 0.0152506 | 0.05 | standard | OK | [G: GW150914_F] |

### Baustein 09

| Prüfzahl | Modell | Rechnung | Einheit | Abweichung | Toleranz | Satz | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|---|---|
| 09_em_gegen_grav | 1,23555·10³⁶ | 1.23555e+36 | 1 | 1.63576e+30 | 5e+30 | standard | OK |  |
| 09_a_J_elektron | 193,08 | 193.08 | fm | -0.00036602 | 0.005 | standard | OK |  |
| 09_bohr_magneton | 9,274·10⁻²⁴ | 9.27401e-24 | J/T | 1.00784e-29 | 5e-28 | standard | OK |  |
| 09_spinzahl_elektron | 2,9·10⁴⁴ | 2.85419e+44 | 1 | -4.58082e+42 | 5e+42 | standard | OK |  |
| 09_spinzahl_proton | 8,5·10³⁷ | 8.46575e+37 | 1 | -3.42474e+35 | 5e+35 | standard | OK |  |
| 09_compton_takt_proton | 1,425·10²⁴ | 1.42549e+24 | rad/s | 4.86241e+20 | 5e+20 | standard | OK |  |
| 09_umlauf_wasserstoff | 4,134·10¹⁶ | 4.13414e+16 | rad/s | 1.37334e+12 | 5e+12 | standard | OK |  |
| 09_umlauf_wasserstoff_welle | 45,6 | 45.5634 | nm | -0.0366475 | 0.05 | standard | OK |  |
| 09_faltungsrate_proton | 8,42·10⁻¹⁵ | 8.41913e-15 | s⁻¹ | -8.65275e-19 | 5e-18 | standard | OK |  |
| 09_faltungstakt_proton | 2,41·10⁶² | 2.41356e+62 | s⁻¹ | 3.56278e+59 | 5e+59 | standard | OK |  |
| 09_neutronenstern_takt | 0,81 | 0.809601 | 1 | -0.000398766 | 0.005 | standard | OK |  |
| 09_neutronenstern_umlauf | 1,04·10⁴ | 10369.3 | rad/s | -30.7352 | 50 | standard | OK |  |
| 09_erdbahn | 1,99·10⁻⁷ | 1.99098e-07 | rad/s | 9.83675e-11 | 5e-10 | standard | OK |  |
| 09_milchstrasse_umlauf | 8,0·10⁻¹⁶ | 8.04227e-16 | rad/s | 4.22734e-18 | 5e-18 | standard | OK |  |
| 09_milchstrasse_umlaufzeit | 248 | 247.569 | Mio. Jahre | -0.430542 | 0.5 | standard | OK |  |
| 09_dephasierung_koerper | 2,46·10⁻¹⁴ | 2.46395e-14 | s | 3.94599e-17 | 5e-17 | standard | OK | [W: T_KOERPER] |
| 09_debye_physiologisch | 0,78 | 0.778668 | nm | -0.00133199 | 0.005 | standard | OK | [G: EPS_R_WASSER_37] [W: T_KOERPER] |
| 09_coulomb_feld_a0 | 110 | 110 | V/m | 0 | 5 | standard | OK |  |
| 09_draht_kanaele_min | 200 | 202.995 | 1 | 2.99544 | 50 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] [W: DRAHT_Q_MERKLICH] |
| 09_draht_kanaele_max | 1900 | 1940.01 | 1 | 40.012 | 50 | standard | OK | [G: CU_AL_E, CU_AL_R, CU_AL_TH, CU_CP, CU_E, CU_K_TH, CU_RHO, CU_RHO_EL] [W: DRAHT_Q_MERKLICH] |
| 09_offen_02_1 | 10¹³ |  |  |  |  | standard | OFFEN (fehlt Formel: Frequenz ω₃ = v_s/r der Gitterschwingung (Eingaben v_s, r fehlen)) |  |

## Gesetze

### Klasse A

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| A_schwarzschild_radius | Schwarzschild-Radius 2GM/c² = r₀ = 6Δs | 2953.25 | 2953.25 | 0.0e+00 | OK |  |
| A_zeitdilatation | Zeitdilatation √(1 − 2GM/rc²) = √(1 − 6Δs/r) | 1 | 1 | 0.0e+00 | OK |  |
| A_lichtablenkung | Lichtablenkung 4GM/(c²b) = 12Δs/b (Sonnenrand) | 8.49001e-06 | 8.49001e-06 | 2.2e-16 | OK |  |
| A_periheldrehung | Periheldrehung 6πGM/(c²a(1−e²)) = 18πΔs/(a(1−e²)) (Merkur) | 5.0186e-07 | 5.0186e-07 | 0.0e+00 | OK | [G: MERKUR_A, MERKUR_E] |
| A_photonensphaere | Photonensphäre 3GM/c² = 9Δs | 4429.88 | 4429.88 | 0.0e+00 | OK |  |
| A_schattenradius | Schattenradius 3√3·GM/c² = 3√27·Δs | 7672.77 | 7672.77 | 0.0e+00 | OK |  |
| A_isco_schwarzschild | ISCO Kerr (j = 0, mitläufig) als Minimum der Bahnenergie | 6 | 6 | 1.6e-08 | OK |  |
| A_isco_kerr_mitlaeufig | ISCO Kerr (j = 0,9, mitläufig) als Minimum der Bahnenergie | 2.32088 | 2.32088 | 1.3e-11 | OK |  |
| A_isco_kerr_gegenlaeufig | ISCO Kerr (j = 0,9, gegenläufig) als Minimum der Bahnenergie | 8.71735 | 8.71735 | 7.5e-09 | OK |  |
| A_hawking_temperatur | Hawking-Temperatur ħc³/(8πGMk) = ħc/(24πkΔs) | 6.17026e-08 | 6.17026e-08 | 2.2e-16 | OK |  |
| A_bekenstein_hawking | Bekenstein-Hawking S = 4πkGM²/(ħc) = 36πk(Δs/ℓ_P)² | 1.04889e+77 | 1.04889e+77 | 0.0e+00 | OK |  |
| A_hawking_lebensdauer | Hawking-Lebensdauer 5120πG²M³/(ħc⁴) = 3840(S/k)Δs/c | 6.61287e+74 | 6.61287e+74 | 3.3e-16 | OK |  |
| A_eddington | Eddington-Grenze 4πGMm_pc/σ_T = (9/2)m_pc³Δs/r_e² | 1.25707e+31 | 1.25707e+31 | 0.0e+00 | OK |  |
| A_gibbons_hawking | Gibbons-Hawking ħH/(2πk) = ħc/(12πkΔs_H), Δs_H = c/6H (heute) | 2.65535e-30 | 2.65535e-30 | 2.2e-16 | OK |  |
| A_kepler | Kepler T² = 4π²a³/(GM) = 4π²a³/(3c²Δs) (1 AE) | 3.15582e+07 | 3.15582e+07 | 0.0e+00 | OK |  |
| A_quadrupolleistung | Quadrupolleistung (32/5)G⁴(m₁m₂)²M/(c⁵r⁵) = (32/5)(c⁵/G)η_m²(3Δs/r)⁵ | 1.45793e+40 | 1.45793e+40 | 4.4e-16 | OK |  |
| A_hulse_taylor | Abnahme der Umlaufzeit Ṗ_b (Peters-Mathews) = Faltungsform (PSR B1913+16) | -2.40257e-12 | -2.40257e-12 | 8.9e-16 | OK | [G: PSR1913_E, PSR1913_M1, PSR1913_M2, PSR1913_PB] |
| A_geodaetisch | geodätische Präzession (3/2)(GM/c²r)ω = (9/2)(Δs/r)ω (Bahn von Gravity Probe B) | 1.01454e-12 | 1.01454e-12 | 2.2e-16 | OK | [W: GPB_A] |
| A_mitfuehrung | Mitführung im Polarorbit GJ/(2c²a³) = 3Δs·a_J·c/(2a³) | 6.26808e-15 | 6.26808e-15 | 1.1e-16 | OK | [G: ERDE_J] [W: GPB_A] |
| A_sternmasse | Sternmasse m_Pl³/m_p² (Landau-Chandrasekhar-Skala) | 3.68503e+30 | 3.68503e+30 | 0.0e+00 | OK |  |
| A_alpha_G | α_G(M) = 3Δs/λ̄ = GM²/(ħc) (M = 10⁻²⁰ kg) | 2.1111e-25 | 2.1111e-25 | 0.0e+00 | OK |  |
| A_ppn_gamma | γ = 1: Raumkrümmung je Masse wie in der ART (aus der Lichtablenkung der Faltungsform) | 1 | 1 | 0.0e+00 | OK |  |
| A_ppn_beta | β = 1: Nichtlinearität der Überlagerung wie in der ART (aus der Periheldrehung, γ aus der Ablenkung) | 1 | 1 | 0.0e+00 | OK | [G: MERKUR_A, MERKUR_E] |
| A_newton_erdbahn | Newton-Grenze: an der Erdbahn a = a_N (Zusatz e^{−s}, s = r*/r ≈ 7342, unter der Rechengenauigkeit) | 1 | 1 | 0.0e+00 | OK |  |

### Klasse B

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| B_beschleunigungsrelation | Beschleunigungsrelation: Form g/(1 − e^{−√(g/g†)}) = örtliche Regel ν(g_N)·g_N | 1.05382e-12 | 1.05382e-12 | 4.6e-15 | OK |  |
| B_tiefbereich_skaleninvariant | Tiefbereich skaleninvariant: a(λr)·λ = a(r) (a_N = 10⁻⁴·a₀, λ = 3) | 1 | 0.999625 | 3.8e-04 | OK |  |
| B_virialsatz_tief | Virialsatz im Tiefbereich σ⁴ = (4/81)GMa₀ (Plummer, örtliche Regel, jedes Profil) | 8.49062e+06 | 8.49168e+06 | 1.2e-04 | OK |  |
| B_tully_fisher_grenzfall | Tully-Fisher v⁴ = GMa₀ als Grenzfall der örtlichen Regel weit außen | 172034 | 172077 | 2.5e-04 | OK |  |
| B_btf_steigung | Baryonische Tully-Fisher: Steigung log M gegen log v aus der Regel (Grenzfall 4) gegen die Messung | 4 | 3.85 | +1.7σ | OK | [G: LELLI_BTF_STEIGUNG] |
| B_toomre_exponentialscheibe | Toomre: (1/2πG)∫V²/R² dR = Σ(0) für die Exponentialscheibe (Newton, abgetastete Kurve) | 0.159154 | 0.159155 | 4.4e-06 | OK |  |
| B_flaechendichte_beta | Zentrale Flächendichte, ideale Exponentialscheibe: Steigung β der Doppelpotenz (α = 1 − β) | 0.539347 | 0.62 | -2.0σ | SPANNUNG | [G: LELLI_BETA_FEHLER, LELLI_SIGMA_DYN] |
| B_flaechendichte_lelli_fit | Zentrale Flächendichte, ideale Exponentialscheibe gegen den Lelli-Fit bei Σ₀ = 10 / 100 / 1000 M☉/pc² [dex], innerhalb der Streuung | 0.224632 | 0 | +0.22 / +0.14 / +0.08 | SPANNUNG | [G: LELLI_SIGMA_DYN, LELLI_STREUUNG] |
| B_flaechendichte_sparc_rest | Zentrale Flächendichte Galaxie für Galaxie (eigene Baryonen): Rest log Σ_gem/Σ_Regel [dex] null | -0.0543969 | 0 | -1.3σ | OK |  |
| B_flaechendichte_sparc_trend | Zentrale Flächendichte Galaxie für Galaxie: Trend des Rests je Dekade Σ₀ null | 0.0483444 | 0 | +1.9σ | OK |  |
| B_crater2_schaetzer | Crater II im Außenfeld: Plummer-Virial gegen den Schätzer McGaugh (σ² ∝ M statt M/2) [km/s] | 1.33125 | 2.12437 | Faktor 2.55 in σ² | HINWEIS | [G: CRATER2_ABSTAND, CRATER2_M, CRATER2_RE] |
| B_aussenfeld_richtung | Äußeres Feld senkt ν, die Richtung zählt: frei > senkrecht > gleichgerichtet (Milchstraße, 50 kpc, 0,03·a₀) | 1 | 1 | ν 6.25 > 5.36 > 4.59 | OK |  |

### Klasse C

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| C_rueckkehr_verlust | Rückkehr mit Verlust: q = e^{−L√(Γ/D)} (Irrfahrt, D = ½, Γ = p, L = 12) | 0.341874 | 0.343225 | 4.0e-03 | OK |  |
| C_rueckkehr_diskret | diskret exakt: [(1 − √(1 − (1−p)²))/(1 − p)]^L gegen e^{−L√(2p)} | 0.341874 | 0.341262 | 1.8e-03 | OK |  |
| C_kac_halbgerade | Kac: mittlere Rückkehrzeit auf der bei N geschlossenen Halbgeraden = μ(X)/μ(0) = 2N (N = 15) | 30 | 30.3019 | 1.0e-02 | OK |  |
| C_watson | Pólya/Watson: Rückkehrsumme W₃ = ∫e^{−t}I₀(t/3)³dt = (√6/32π³)Γ(1/24)Γ(5/24)Γ(7/24)Γ(11/24) | 1.51639 | 1.51639 | 5.2e-09 | OK |  |
| C_bose_leiter_1 | Bose-Leiter ∫s^1n_B ds = Γ(1+1)ζ(1+1) | 1.64493 | 1.64493 | 1.1e-16 | OK |  |
| C_bose_leiter_2 | Bose-Leiter ∫s^2n_B ds = Γ(2+1)ζ(2+1) | 2.40411 | 2.40411 | 2.2e-16 | OK |  |
| C_bose_leiter_3 | Bose-Leiter ∫s^3n_B ds = Γ(3+1)ζ(3+1) | 6.49394 | 6.49394 | 1.3e-15 | OK |  |
| C_wien_verschiebung | Wien: Maximum von s³n_B bei s = 3 + W₀(−3e⁻³) | 2.82144 | 2.82144 | 0.0e+00 | OK |  |
| C_kippstelle_fenster | Kippstelle: Hysteresefenster ∝ b^{3/2} (Verhältnis b = 0,4 zu 0,1 = 8) | 8 | 8.00082 | 1.0e-04 | OK |  |
| C_fdt_daempfung_1 | Fluktuation-Dissipation ⟨x²⟩ = kT/k, Dämpfung γ = 1 (Langevin) | 0.5 | 0.503448 | 6.9e-03 | OK |  |
| C_fdt_daempfung_10 | Fluktuation-Dissipation ⟨x²⟩ = kT/k, Dämpfung γ = 10 (Langevin) | 0.5 | 0.50568 | 1.1e-02 | OK |  |
| C_onsager_casimir | Onsager-Casimir: Entropieerzeugung X·(L_s + L_a)·X = X·L_s·X (nur der symmetrische Teil) | 4.83589 | 4.83589 | 0.0e+00 | OK |  |
| C_ramanujan_tau | Ramanujan τ(n): Werte, Multiplikativität, Hecke, Kongruenz mod 691, Deligne-Schranke | 1 | 1 | τ(2) = -24 | OK |  |
| C_thetareihen_bilder | Bilder im Kaleidoskop (Zelle als Torus): r₂(n) = 4(d₁ − d₃), sechseckig 6(d₁ − d₂) mod 3, n ≤ 25 | 1 | 1 | n = 1…25; quadratisch 4, 4, 0, 4, 8, sechseckig 6, 0, 6, 6, 0 für n = 1…5 | OK |  |
| C_barkhausen | Barkhausen: Amplitude des Grenzzyklus A* = 2 (van der Pol, μ = 0,05; Schleifenverstärkung 1 im Mittel) | 2 | 1.99987 | 6.5e-05 | OK |  |
| C_kac_budget_ketten | Kac für jede Zerlegung: Σ 1/⟨T_i⟩ = 1 (drei geschlossene Ketten, 3, 5, 8 Zustände) | 1 | 1.00006 | 6.2e-05 | OK |  |

### Klasse Z

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| Z_poynting_joule | Poynting ∂u/∂t + ∇·S = −J·E: Feldverlust = Joule-Wärme | 15.2775 | 15.2775 | 8.3e-11 | OK |  |
| Z_maxwell_struktur | Maxwell als Strukturbedingung: ∇·B = 0 und ∇·D = ρ bleiben erhalten (Gitter 12³, Zufallsfelder) | 1 | 1 | 4e-15 | OK |  |
| Z_lorentz_ohne_arbeit | Lorentzkraft f = ρE + J×B tangential: ¦v¦ im Magnetfeld erhalten (Boris-Schritt) | 1.06301 | 1.06301 | 6.0e-15 | OK |  |
| Z_durchgehen | Thermisches Durchgehen bei Schleifenverstärkung q = 1: I²R₀α/G = 1/e (R = R₀e^{αT}) | 0.367879 | 0.367879 | 7.5e-09 | OK |  |
| Z_verstaerkung_verlangsamung | Verstärkung ist Verlangsamung: τ = τ₀/(1 − q) = τ₀·Σqⁿ (Wärmeschleife) | 1.34984 | 1.34984 | 2.3e-07 | OK |  |
| Z_schliessung_nullmode | Schließung: Ring hat einen Nullmode (Gesamtmenge bleibt), Strecke mit Rand keinen | 1 | 1 | -1e-14 / 9.87 | OK |  |
| Z_spektralluecke_kreis | Spektrallücke des geschlossenen Kreises 4π²D/L² (D = L = 1) | 39.4784 | 39.4776 | 2.1e-05 | OK |  |
| Z_fluchtrate_strecke | Fluchtrate der aufgeschnittenen Strecke π²D/L² | 9.8696 | 9.86955 | 5.1e-06 | OK |  |
| Z_parametrisch_wirkung | Parametrisch: Wirkung E/ω bleibt bei langsamer Änderung von g (g × 4 über ≈ 700 Perioden) | 0.5 | 0.5 | 3.4e-07 | OK |  |
| Z_vorspannung_resonanz | Gravitation: Vorspannung verschiebt Resonanz ω²(P)/ω²(0) = 1 − P/P_E (P = P_E/2; Knicken bei P_E = π²EI/L²) | 0.5 | 0.499995 | 9.3e-06 | OK |  |
| Z_selbsterregung | Strömung/Selbsterregung: Grenzzyklus genau bei q_self = a/(R + ¾bV²) = 1, V*² = 4(a − R)/(3b) | 0.258199 | 0.258202 | 1.3e-05 | OK |  |
| Z_abstrahlung | Akustik/Abstrahlung: reiner Verlust (q < 0), Beitrag −R_rad/2L zur Querrate λ_⊥ | 0.04 | 0.04 | 3.4e-11 | OK |  |
| Z_ermuedung | Ermüdung (Energie → Zustand → Struktur): Ḋ ∝ (1 − D)^{−m} bricht nach t_f = 1/((m + 1)Ḋ₀) (m = 2) | 0.333333 | 0.333333 | 2.2e-16 | OK |  |
| Z_alterung_arrhenius | Chemie/Alterung: Empfindlichkeit ∂ln Ȧ/∂T = E_a/T² (Arrhenius), einseitige Drift | 7.5 | 7.5 | 3.3e-11 | OK |  |
| Z_rauschen_banddicke | Rauschen: Banddicke ∝ √(D_⊥/¦λ_⊥¦), wächst für λ_⊥ → 0 (⟨x²⟩-Verhältnis bei λ = 1 und ¼) | 4 | 3.99588 | 1.0e-03 | OK |  |
| Z_kriechen | Kriechen: Spannung relaxiert σ = σ₀e^{−Et/η} (q < 0, stabilisierend, irreversibel) | 0.201897 | 0.201896 | 3.2e-06 | OK |  |
| Z_generic_rollen | GENERIC: 𝓛 = −𝓛ᵀ, 𝓛∇S = 0, 𝓜 = 𝓜ᵀ ≥ 0, 𝓜∇E = 0; E fest, S steigt (Oszillator mit Wärme) | 1 | 1 | ΔE 6e-13, ΔS +0.405 | OK |  |
| Z_einrasten | Einrasten genau für ¦Δω¦ ≤ 2k¦cos φ¦; bei φ = π/2 nie (156 Paare Δω, φ) | 1 | 1 | 40 von 156 gerastet | OK |  |
| Z_rastrate | Rate am Rastpunkt −√((2k cos φ)² − Δω²) (Δω = 0,6, k = 0,5, φ = 0; Abklingen einer Störung) | -0.8 | -0.799999 | 1.7e-06 | OK |  |
| Z_einrastzunge_02 | Einrastzunge 3:2: Rand bei ¦nω₂ − mω₁¦ = (m + n)k (k = 0.02) | 1 | 1.01 | 1.0e-02 | OK |  |
| Z_einrastzunge_04 | Einrastzunge 3:2: Rand bei ¦nω₂ − mω₁¦ = (m + n)k (k = 0.04, Breite linear) | 1 | 1.01 | 1.0e-02 | OK |  |
| Z_kuramoto | Kuramoto N = 40: volles Einrasten (N − 1)k/W, Simulation gegen Selbstkonsistenz | 0.635406 | 0.6375 | 3.3e-03 | OK |  |
| Z_torusdimension_grenze | Grenze N_max ≈ 1 + 2W/(πk): (N − 1)k_c/W → 2/π (N = 4000) | 0.636619 | 0.63662 | 1.9e-06 | OK |  |
| Z_koexistenz | Koexistenz zweier Zyklen genau für θ₁₂/β₂ < α₁/α₂ < β₁/θ₂₁ | 1 | 1 | 5 Verhältnisse | OK |  |
| Z_amplitudentod | Amplitudentod genau für K_d > μ und ¦Δω¦ > 2√(μ(2K_d − μ)) | 1 | 1 | 3 Punkte | OK |  |
| Z_preisach | Preisach: gleiche Folge der Extremwerte ⇒ gleicher Zustand; Rückkehrpunkt-Gedächtnis | 1 | 1 | Folge der Extremwerte | OK |  |
| Z_kippstelle_wurzel | Kippstelle: 1 − q = 2(3a)^{1/4}δ^{1/2} vor der Falte (a = 1, δ = 10⁻⁸) | 0.000263215 | 0.000263226 | 4.4e-05 | OK |  |
| Z_kippstelle_sprung | Kippstelle: Sprung an der Falte √(3a) | 1.73205 | 1.73205 | 3.3e-16 | OK |  |
| Z_dreieck_gleich | Oszillatorischer Einsatz im Dreieck Energie–Struktur–Zustand, Zeiten (1.0, 1.0, 1.0): ¦q₃¦ = (Σt)(Σ1/t) − 1 | 8 | 8 | 0.0e+00 | OK |  |
| Z_dreieck_getrennt | Oszillatorischer Einsatz im Dreieck Energie–Struktur–Zustand, Zeiten (0.1, 1.0, 30.0): ¦q₃¦ = (Σt)(Σ1/t) − 1 | 342.137 | 342.137 | 2.2e-16 | OK |  |
| Z_ruecklaeufe_netz | Rückläufe im Netz: (I − A)⁻¹ = ΣAⁿ (Spektralradius 0,9) | 59.2515 | 59.2515 | 2.2e-16 | OK |  |
| Z_gedaechtnis_bindet | Gedächtniskern bindet ohne Energieverlust: Umkehrpunkt aus der Erhaltung H = v²/2 − 1/x + γy²/2 (v > √2) | 9.74747 | 9.74747 | 1.7e-12 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_maxwell_speicher_0_5 | Maxwell-Element ωt_K = 0.5: Speicheranteil (ωt_K)²/(1 + (ωt_K)²) | 0.2 | 0.2 | 1.1e-16 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_maxwell_verlust_0_5 | Maxwell-Element ωt_K = 0.5: Verlustanteil ωt_K/(1 + (ωt_K)²) | 0.4 | 0.4 | 2.2e-16 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_maxwell_speicher_2_0 | Maxwell-Element ωt_K = 2: Speicheranteil (ωt_K)²/(1 + (ωt_K)²) | 0.8 | 0.8 | 2.2e-16 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_maxwell_verlust_2_0 | Maxwell-Element ωt_K = 2: Verlustanteil ωt_K/(1 + (ωt_K)²) | 0.4 | 0.4 | 2.2e-16 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_bindungsoperator | L_eff = 𝒫L: 𝒫_q vertauscht mit dem Generator, löscht keine Freiheitsgrade (voller Rang), q → 1 gibt den Projektor 𝒫 | 1 | 1 | Rang 64 von 64 | OK | [W: Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)] |
| Z_poisson_inversion | Poisson-Inversion: die gelesene Dichte ρ_X = ∇·(a − a_N)/(4πG) ist Ausgabe der Regel (außen isothermer Halo) | 5.41576e-25 | 5.41576e-25 | 2.6e-12 | OK |  |
| Z_faltung_additiv | Faltung additiv (ΣΔs), die Regel nicht additiv je Teilchen (Operator statt Quelle) | 1 | 1 | Faltung 0e+00, Regel -0.26 | OK |  |
| Z_linsentest_kugel_0_3 | Test 1: Kugel, Linsenmasse = dynamische Masse, η = 1 (b = 0.3 r*) | 5.11281e-06 | 5.11281e-06 | 6.3e-10 | OK |  |
| Z_linsentest_kugel_1_0 | Test 1: Kugel, Linsenmasse = dynamische Masse, η = 1 (b = 1 r*) | 2.81277e-06 | 2.81277e-06 | 9.7e-10 | OK |  |
| Z_linsentest_kugel_5_0 | Test 1: Kugel, Linsenmasse = dynamische Masse, η = 1 (b = 5 r*) | 2.20419e-06 | 2.20419e-06 | 1.2e-09 | OK |  |
| Z_ablenkung_weit_aussen | Ablenkung weit außen α → 2π(v_f/c)² + 6Δs/b (b = 200 r*) | 2.07232e-06 | 2.07232e-06 | 1.0e-06 | OK |  |
| Z_trenntest_wachstum | Test 2: Wachstum mit Skalen- und Epochenabhängigkeit (Quelle: nur über die Dichte) | 1 | 1 | Skala 3.53, Epoche 0.890 | OK | [G: OMEGA_R] |
| Z_trenntest_umfeld | Test 3: Form der gelesenen Verteilung an das Umfeld gebunden (äußeres Feld 0,03·a₀, 50 kpc) | 1 | 1 | ρ_X mit/ohne Außenfeld 0.43 | OK |  |
| Z_trenntest_gedaechtnis | Test 4: Pfadabhängigkeit über w (Gedächtnis statt Markow) | 1 | 1 | w 0.9854 gegen 0.8219 | OK |  |
| Z_trenntest_aequivalenz | Test 5: Äquivalenzprinzip: alle Materie folgt der Regel gleich, trotz verschiedener Glättungslänge (Erde, Mond) | 1 | 1 | ℓ Erde/Mond 4.33, a gleich | OK | [G: M_MOND, RHO_B_SONNE, V0_SONNE] |
| Z_trenntest_kanaele | Test 6: kanalspezifische Reste: in Scheiben innen η = M_dyn/M_lens > 1, außen 1 | 1 | 1 | Licht/Bahn innen 0.80, bei 5R_d 0.997 | OK |  |
| Z_triade_hodge | Triade als Form: jedes Feld auf dem geschlossenen Raum = Hintergrund (harmonisch, gleichförmig) + Zentrum (Gradient) + Wirbel (Rotation), eindeutig und orthogonal | 1 | 1 | Rest 1e-15, Überlapp 5e-18, div 1e-13, rot 6e-14, harmonisch 0e+00 | OK |  |

### Klasse M

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| M_neigung_A1 | (A1) g = x/y²·g⊕ gegen GM/R² (Sonne; g⊕ = 9,81 gesetzt) | 273.914 | 274.2 | 1.0e-03 | OK | [G: R_earth] |
| M_zeitdifferential_bahn | Zeitdifferential einer Bahn δ = 3Δs/r = v²/c² (Erdbahn, v aus dem Jahr) | 9.87063e-09 | 9.871e-09 | 3.8e-05 | OK |  |
| M_brechzahl_strahl | Brechzahl n = 1 + 6Δs/r als Strahlintegral: α = 12Δs/b (Takt und Faltung je die Hälfte) | 8.49001e-06 | 8.49001e-06 | 0.0e+00 | OK |  |
| M_schnitt_planck | Schnitt 3Δs = λ̄ bei m_Pl (Nullstelle) | 2.17643e-08 | 2.17643e-08 | 1.6e-15 | OK |  |
| M_schnitt_compton | Schnitt r* = λ̄ bei M_c = (ħ²a₀/(Gc²))^{1/3} (Nullstelle) | 5.88617e-29 | 5.88617e-29 | 1.3e-15 | OK |  |
| M_schnitt_oben | Schnitt r* = 3Δs bei M_U = c⁴/(Ga₀) = c³τ_D/G = m_Pl/ε_L (Nullstelle) | 1.10023e+54 | 1.10023e+54 | 1.2e-15 | OK |  |
| M_schnitt_unten | Schnitt λ̄ = ℓ₀ bei M_min = ħ/(c²τ_D) = ε_L·m_Pl (Nullstelle) | 4.30533e-70 | 4.30533e-70 | 2.0e-15 | OK |  |
| M_laengenleiter | Längenleiter M_c³·M_U = m_Pl⁴ (aus den Nullstellen) | 2.24379e-31 | 2.24379e-31 | 1.5e-14 | OK |  |
| M_produkt_fest | r*²·λ̄ = ℓ_P²·ℓ₀ für jede Masse (Stichprobe 10⁻²⁰ … 10⁴² kg) | 2.13436e-43 | 2.13436e-43 | 4.4e-16 | OK |  |
| M_produkt_faltung_compton | Produkt fest: 3Δs·λ̄ = ℓ_P² für jede Masse (Δs·λ̄/ℓ_P² = ⅓) | 8.7076e-71 | 8.7076e-71 | 2.2e-16 | OK |  |
| M_compton_mittel | Geometrisches Mittel λ̄_e² = r_e·a_B | 1.49119e-25 | 1.49119e-25 | 7.6e-12 | OK |  |
| M_sternskala | Sternskala: Δs_*²·Δs_p = (λ̄_p/3)³ ⇒ M_* = m_Pl³/m_p² | 3.68503e+30 | 3.68503e+30 | 0.0e+00 | OK |  |
| M_nullflaeche_takt | Nullfläche des Takts r_E(θ) = 3Δs + √(9Δs² − a_J²cos²θ) (Nullstelle von g_tt, θ = 0,7) | 2788.58 | 2788.58 | 3.1e-15 | OK |  |
| M_umlaufdifferenz_10 | T_mitläufig − T_gegenläufig = 4πa_J/c unabhängig von r (Metrik, r = 10 GM/c²) | 3.71373e-05 | 3.71373e-05 | 5.7e-15 | OK |  |
| M_umlaufdifferenz_1000 | T_mitläufig − T_gegenläufig = 4πa_J/c unabhängig von r (Metrik, r = 1000 GM/c²) | 3.71373e-05 | 3.71373e-05 | 8.4e-12 | OK |  |
| M_spin_halb | Spin ½: a_J·Δs = ℓ_P²/6 (Elektron); j = ½(m_Pl/M)² | 4.3538e-71 | 4.3538e-71 | 2.2e-16 | OK |  |
| M_magnetisches_moment | Magnetisches Moment μ = Q_el·c·a_J (g = 2) = μ_B = eħ/(2m_e) beim Elektron | 9.27401e-24 | 9.27401e-24 | 0.0e+00 | OK |  |
| M_rueckkehr_geschlossen | Rückkehr bei Λ geschlossen: G(x) = cosh((L − x)/λ)/cosh(L/λ) (L = 3λ, x = 1λ) | 0.373691 | 0.373691 | 4.4e-08 | OK |  |
| M_rueckkehr_offen | Rückkehr offen: q = e^{−s} (Differenzenrechnung, x = 1λ) | 0.367879 | 0.367879 | 4.2e-08 | OK |  |
| M_chi_elastizitaet_a0_0_3 | d ln a/d ln a₀ = χ/2 (s = 0.3) | 0.428744 | 0.428744 | 5.9e-11 | OK |  |
| M_chi_flussdimension_0_3 | d ln a/d ln a_N = 1 − χ/2, Flussdimension 3 − χ (s = 0.3) | 0.571256 | 0.571256 | 8.0e-11 | OK |  |
| M_chi_elastizitaet_a0_1_0 | d ln a/d ln a₀ = χ/2 (s = 1) | 0.290988 | 0.290988 | 1.1e-10 | OK |  |
| M_chi_flussdimension_1_0 | d ln a/d ln a_N = 1 − χ/2, Flussdimension 3 − χ (s = 1) | 0.709012 | 0.709012 | 5.6e-11 | OK |  |
| M_chi_elastizitaet_a0_3_0 | d ln a/d ln a₀ = χ/2 (s = 3) | 0.0785935 | 0.0785935 | 1.1e-09 | OK |  |
| M_chi_flussdimension_3_0 | d ln a/d ln a_N = 1 − χ/2, Flussdimension 3 − χ (s = 3) | 0.921406 | 0.921406 | 2.0e-11 | OK |  |
| M_spektraldimension_tief | Spektraldimension d_s = 2·d ln a/d ln a_N → 1 im Tiefbereich (s = 10⁻³) | 1.0005 | 1.0005 | 1.9e-10 | OK |  |
| M_zusatz_kern_nome | Zusatz der Rückkehr a − a_N = χ·(Δs/r)(3c²/r*) (Sonne, s = 0,7) | 5.31688e-11 | 5.31688e-11 | 0.0e+00 | OK |  |
| M_zusatzpotential_0_5 | Φ_X = v_f²·ln ν(s) als Integral von a − a_N (Punktmasse, s = 0.5) | 2.76054e+10 | 2.76054e+10 | 2.1e-15 | OK |  |
| M_zusatzpotential_2_0 | Φ_X = v_f²·ln ν(s) als Integral von a − a_N (Punktmasse, s = 2) | 4.3036e+09 | 4.3036e+09 | 1.2e-14 | OK |  |
| M_zusatzpotential_hernquist | Φ_X Hernquist = v_f²·ln[ν(s)/ν(s₀)] (b_H = 0,4 r*, r = 1,5 r*) | 2.39087e+10 | 2.39087e+10 | 2.2e-16 | OK |  |
| M_re_tau_punkt_0_5 | (Re τ_M)² = κ²/Ω² = 1 + χ, Punktmasse (r = 0.5 r*) | 1.31304 | 1.31304 | 2.5e-08 | OK |  |
| M_re_tau_punkt_2_0 | (Re τ_M)² = κ²/Ω² = 1 + χ, Punktmasse (r = 2 r*) | 1.77075 | 1.77075 | 9.8e-09 | OK |  |
| M_re_tau_hernquist | (Re τ_M)² = 3 + γ_N(1 − χ/2), Hernquist (r = 1 r*) | 2.06073 | 2.06073 | 3.6e-09 | OK |  |
| M_entropie_ruecklauf_0_3 | Entropie der Rücklaufzahl S_n = ln ν + χ (s = 0.3) | 2.20771 | 2.20771 | 0.0e+00 | OK |  |
| M_entropie_ruecklauf_1_0 | Entropie der Rücklaufzahl S_n = ln ν + χ (s = 1) | 1.04065 | 1.04065 | 2.2e-16 | OK |  |
| M_fisher_metrik | Fisher-Metrik der Rücklaufzahl g_F = −dν/ds = ν(ν − 1) = 1/(4 sinh²(s/2)) (s = 1) | 0.920674 | 0.920674 | 3.3e-16 | OK |  |
| M_bernoulli_reihe | M_schein/M = 1/(1 − e^{−r*/r}) = r/r* + ½ + r*/(12r) − (r*/r)³/720 + … (r = 3r*) | 3.52773 | 3.52773 | 3.8e-08 | OK |  |
| M_verzweigung_borel | Kritische Verzweigung (Borel): Summe über P(S = n) gegen T = z·e^{T−1} (Γ = 10⁻⁴) | 0.0140755 | 0.0140755 | 2.5e-13 | OK |  |
| M_verzweigung_wurzel | Kritische Verzweigung: 1 − E[e^{−ΓS}] = √(2Γ) (Γ = 10⁻⁸), der Exponent ½ der Rückkehr | 0.000141411 | 0.000141415 | 2.4e-05 | OK |  |
| M_verzweigung_potenz | Kritische Verzweigung: P(S = n)·n^{3/2}·√(2π) → 1 (n = 10⁷) | 1 | 1 | 2.2e-08 | OK |  |
| M_sehne_kreis | Mittlere Sehne im Kreis πA/P (η₂ = π; Monte Carlo) | 1.5708 | 1.57118 | 2.4e-04 | OK |  |
| M_sehne_kugel | Mittlere Sehne in der Kugel 4V/S (η₃ = 4; Monte Carlo) | 1.33333 | 1.33294 | 3.0e-04 | OK |  |
| M_dimension_drei | V/(S·R) = 1/d für die d-Kugel ⇒ Δs = GM/(3c²) für d = 3 | 0.333333 | 0.333333 | 1.1e-16 | OK |  |
| M_modulformen | Modulformen: q dE₂/dq = (E₂² − E₄)/12 usw., Δ = (E₄³ − E₆²)/1728 = q∏(1 − qⁿ)²⁴, D₁₂Δ = 0 | 1 | 1 | q = 0,05 | OK |  |
| M_herkunft_zahlen | Herkunft der Zahlen: 24, 240, 504 = (2π)^k/(Γ(k)ζ(k)), 1728 = 3·240 + 2·504, 12 = 2/(1 − ½ − ⅓), 4π⁴/15 = 4(2π)⁴/240 | 1 | 1 | 24, 240, 504, 1728, 12, 4π⁴/15 | OK |  |
| M_lambert_reihe | Lambert-Reihe Σσ₁(n)e^{−ns} = π²/(6s²) − 1/(2s) + 1/24 (Modultransformation, s = 0,1) | 159.535 | 159.535 | 4.4e-16 | OK |  |
| M_hecke | Hecke: Σσ₁(pm)q^m = (1 + p)ℓ(q) − p·ℓ(q^p) (Übertrag zwischen Skalenebenen) | 1 | 1 | p = 2, 3, 5 | OK |  |
| M_rogers_ramanujan | Rogers-Ramanujan Σq^{n²}/(q;q)_n = ∏1/((1 − q^{5n+1})(1 − q^{5n+4})) (q = ½) | 2.17267 | 2.17267 | 6.7e-16 | OK |  |
| M_ausschlussregel | Ausschlussregel Abstand ≥ 2: Zustände je Platz → φ (Fibonacci) | 1.61803 | 1.61803 | 0.0e+00 | OK |  |
| M_ramanujan_projektoren | Ramanujan-Summen: E_q Projektoren, ΣE_q = I, E_qE_q' = 0, E₁ = 𝒫 (ℤ₁₂) | 1 | 1 | ℤ₁₂ | OK |  |
| M_resolvente_kreis | 𝒫_q = (1 − q)(I − qT)⁻¹: Eigenwerte auf dem Kreis um 1/(1 + q) mit Radius q/(1 + q) | 0.270073 | 0.270073 | 8.0e-15 | OK |  |
| M_goldene_wicklung | Goldene Wicklung: n·‖nφ‖ → 1/√5 (n = F₂₀) | 0.447214 | 0.447214 | 1.4e-08 | OK |  |
| M_drei_abstaende | Drei-Abstände-Satz: Rückkehrzeiten nehmen ≤ 3 Werte an, der größte = Summe; Mittel = Kac | 1 | 1 | Kac 1.0000/1.0000/1.0001 | OK |  |
| M_kapitza | Kapitza-Pendel: zeitgemitteltes Einschließen oben genau für (Aω)² > 2gℓ ((Aω)²/gℓ = 3,24 / 1,21) | 1 | 1 | 3,24 / 1,21 | OK |  |
| M_aequivalenzschicht | Äquivalenzschicht: Schale und Vollkugel geben außen dasselbe Feld GM/r² (Dichte nicht eindeutig) | 1 | 1 | 1.0000000000 / 1.0000000000 | OK |  |
| M_kreuzterm | Kreuzterm beim Verschmelzen 8πG·M₁M₂k/(ħc) (aus der Horizontfläche) | 4.40534e+78 | 4.40534e+78 | 4.4e-16 | OK |  |
| M_kreuzterm_halb | Gleiche Massen: Kreuzterm = ½ der Endentropie | 0.5 | 0.5 | 0.0e+00 | OK |  |
| M_entropie_alpha_G | Entropie der Falte S = 36πk(Δs/ℓ_P)² = 4πk·α_G(M) (Sonne) | 1.04889e+77 | 1.04889e+77 | 2.2e-16 | OK |  |
| M_rueckgabe_rate | Rückgabe dM/dt = −ħc⁴/(15360πG²M²) (Ableitung von M(t) aus der Lebensdauer) | -3.96284e-09 | -3.96284e-09 | 3.3e-07 | OK |  |
| M_horizont_doppelt | Horizont bei gleichem Radius doppelt so heiß: (ħH/2πk)/(ħc/(4πk·r₀)), r₀ = c/H | 2 | 2 | 4.4e-16 | OK |  |
| M_nariai_drittel | Größte mögliche Falte (Nariai-Grenze, numerisch) hat ⅓ der Horizontentropie | 0.333333 | 0.333333 | 1.5e-10 | OK |  |
| M_falten_schranke | Alle Falten zusammen höchstens Ω_b²·S_H (eine Falte aus ΣΔs_b der Hubble-Kugel) | 0.00242489 | 0.00242489 | 5.6e-16 | OK | [G: Ω_b] |
| M_kein_endzustand | Kein Endzustand: −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k) aus Ṡ_H ≥ Ṡ_i | 4.82088e-51 | 4.82088e-51 | 0.0e+00 | OK |  |
| M_vakuumteil_theta | Λ_s = π⁴/(15Λ²) ⇒ Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH) | 0.244327 | 0.244327 | 4.4e-16 | OK |  |
| M_entropie_vakuum | S_Λ = S_H(H_s) = (45/π³)(Λ/ℓ_P)² bei H_s² = Λ_s c²/3 | 9.27215e+122 | 9.27215e+122 | 0.0e+00 | OK |  |
| M_entropie_vakuum_heute | S_Λ = S_H/Ω_s: die Entropie des Endzustands aus der heutigen Horizontentropie | 9.27215e+122 | 9.27215e+122 | 1.1e-16 | OK |  |
| M_endzustand_H | Endzustand H_s/H₀ = √Ω_s (nur der Vakuumteil) | 0.494295 | 0.494295 | 1.1e-16 | OK |  |
| M_theta_kuenftig | Θ im Endzustand (H → H_s) → √45/π | 2.13529 | 2.13529 | 2.2e-16 | OK |  |
| M_strukturwirkung_0_5 | Strukturwirkung: Q̃(s) + 4π⁴/15 = Q_s(s⁴), der Vakuumteil ist ihre Konstante (s = 0.5) | 0.137494 | 0.137494 | 3.6e-14 | OK |  |
| M_strukturwirkung_1_0 | Strukturwirkung: Q̃(s) + 4π⁴/15 = Q_s(s⁴), der Vakuumteil ist ihre Konstante (s = 1) | 0.899221 | 0.899221 | 1.4e-14 | OK |  |
| M_strukturwirkung_3_0 | Strukturwirkung: Q̃(s) + 4π⁴/15 = Q_s(s⁴), der Vakuumteil ist ihre Konstante (s = 3) | 10.2089 | 10.2089 | 3.1e-13 | OK |  |
| M_spaltbarkeit | Spaltbarkeit: Rückstellung ∝ (1 − x), Kippstelle bei x = E_C/(2E_S) = 1 (Ellipsoid numerisch) | 1 | 1 | 4.3e-08 | OK |  |
| M_bindung_viertel | Größte Bindung je Nukleon (Oberfläche + Coulomb) bei x = ¼, E_C = E_S/2 (für alle Koeffizienten) | 0.25 | 0.25 | 3.9e-08 | OK |  |
| M_aufbauzeit | Aufbauzeit der Rückkehr τ(x) = x/(2√(Dk)) = τ_D/(2s) (x = 3λ) | 1.5 | 1.49853 | 9.8e-04 | OK |  |
| M_leiterzahl | ε_L = ℓ_P/ℓ₀ = a₀/a_P (a_P = c²/ℓ_P) | 1.97816e-62 | 1.97816e-62 | 0.0e+00 | OK |  |
| M_spiegel_leiter | Spiegel der Leiter: M_min·M_U = m_Pl² (aus den Nullstellen) | 4.73687e-16 | 4.73687e-16 | 8.9e-16 | OK |  |
| M_min_planck_laenge | Bei M_min zugleich r* = ℓ_P (Nullstelle von λ̄ = ℓ₀ eingesetzt) | 1.61626e-35 | 1.61626e-35 | 1.0e-15 | OK |  |
| M_zelle_masse | M_U = ⅔·ρ_m·L³ (Zelle aus a₀ = 3c⁴/(2Gρ_mL³)) | 1.10023e+54 | 1.10023e+54 | 1.2e-14 | OK |  |
| M_entropie_leiter | S_max(Λ)/k = (45/π³)(Λ/ℓ_P)² = (45/4π³)·M_U/M_min | 9.27215e+122 | 9.27215e+122 | 3.2e-15 | OK |  |
| M_faltungstakt_MU | Faltungstakt der Zellmasse: c/(3Δs(M_U)) = a₀/c | 3.66921e-19 | 3.66921e-19 | 0.0e+00 | OK |  |
| M_ds_MU_lambda | Faltung der Zellmasse: Δs(M_U)/Λ = ⅔ | 0.666667 | 0.666667 | 1.3e-15 | OK |  |
| M_diffusion_spiegel | D·ħ/M = (ℓ_P c)² für jede Masse (D = r*²/τ_D aus der Rückkehr; Rechnung, nicht im Modell) | 1 | 1 | 2.2e-16 | OK |  |
| M_quantendiffusion_spiegel | Quantendiffusion der Spiegelmasse ħ/(2M′) = D/2, M′ = m_Pl²/M (Sonne) | 2.21341e+11 | 2.21341e+11 | 2.2e-16 | OK |  |
| M_faltungstakt_spiegel | Faltungstakt c/(3Δs) = Compton-Takt der Spiegelmasse (Proton) | 2.41356e+62 | 2.41356e+62 | 0.0e+00 | OK |  |
| M_compton_faltung | ω_C² = Γ_G·c/(3Δs) (Proton): Compton-Takt als Mittel von Faltungsrate und Faltungstakt | 2.03201e+48 | 2.03201e+48 | 2.2e-16 | OK |  |
| M_carter_frequenzen | Carter als Frequenzleiter: Γ_G·ω₁³/(αω₂)⁴ = 𝒞 | 2.94412 | 2.94412 | 2.0e-11 | OK |  |
| M_barriere_min | Barriere bei r* = Ruheenergie der untersten Schnittmasse: kT_I = ħa₀/c = M_min c² | 2.80262e-30 | 2.80262e-30 | 2.2e-15 | OK |  |
| M_schluss_abweichung | Geschlossene Rückkehr: Abweichung bei r = 2Λ/π gleich 1/tanh(π/2) − 1 (Γ = kleinster Eigenwert) | 0.0903314 | 0.0903701 | 4.3e-04 | OK |  |
| M_vakuum_planck | Λ_s·ℓ_P² = (4π⁴/15)·ε_L² | 1.01646e-122 | 1.01646e-122 | 0.0e+00 | OK |  |
| M_einseitigkeit | Einseitigkeit m = 1 bei 20 kpc, g_e = 0,03·a₀ in der Ebene: erste Ordnung gegen volle Rechnung | 0.129472 | 0.127913 | 1.2e-02 | OK |  |
| M_verlust_selbstkonsistent | Γ = a/c selbstkonsistent: d ln v/d ln r → −1/6 (fallende Kurven; die Regel nutzt Γ = a_N/c) | -0.166667 | -0.166717 | 3.0e-04 | OK |  |
| M_gewicht_maxwell | Frühere Form des Gewichts (Rechnung): Speicheranteil eines Maxwell-Elements mit t_K = τ_D (x = θτ_D; Zeitbereich) | 1 | 1 | 6.9e-11 | OK |  |
| M_gewicht_laplace | Gewicht der Schließungsregel aus der Rückkehr: 1 − w = Laplace-Transformierte der Dichte erster Rückkehrzeiten = e^{−√x} (x = 3) | 0.176921 | 0.176921 | 1.7e-13 | OK |  |
| M_maxwell_verlust_zeit | Verlustanteil des Maxwell-Elements x/(1 + x²) (Zeitbereich) | 1 | 1 | 2.1e-11 | OK |  |
| M_gedaechtnis_kreisbahn | Gedächtniskern auf der Kreisbahn: rω²/(λ_m² + ω²) = r·(Speicheranteil, x = ω/λ_m); numerisch | 1.83486 | 1.83486 | 1.3e-15 | OK |  |
| M_gedaechtnis_radial | Radialschwingung im Kern: Speicher κ²/(λ²+κ²), Verlust κλ/(λ²+κ²) (rundet Bahnen); numerisch | 1.16854 | 1.16854 | 3.0e-13 | OK |  |
| M_umlauf_takt | Ω·τ_D = (a/a₀)·(c/v) (Milchstraße bei r*) | 2191.83 | 2191.83 | 0.0e+00 | OK |  |
| M_hubble_kugel_nu | An der Hubble-Kugel ist g = cH/2 = πa₀ (heute, a₀ = cH₀/2π): ν(πa₀) = 1/(1 − e^{−√π}) | 1.2047 | 1.2047 | 0.0e+00 | OK |  |
| M_preisach_hintergrund | Preisach: 1 − w im mitatmenden Hintergrund gleich der jetzigen Form (θ fällt stetig) | 0.0146107 | 0.0146107 | 3.2e-15 | OK |  |
| M_nu_ohne_gefaelle | Ohne Gefälle ist die Rückkehr das Gewicht: ν = 1/w (θ̃ = 3H₀) | 1.01483 | 1.01483 | 0.0e+00 | OK |  |
| M_gewicht_eins | w bei θ = a₀/c: 1 − 1/e | 0.632121 | 0.632121 | 0.0e+00 | OK |  |
| M_gewicht_halb | w = ½ bei θ = (ln 2)²·a₀/c (Nullstelle) | 1.76288e-19 | 1.76288e-19 | 2.2e-16 | OK |  |
| M_ring_wurzel_nu | Ring bei Neigung → 0: ω_p örtlich/Newton = √ν (Miyamoto-Nagai, 20 kpc) | 1.69585 | 1.69585 | 3.2e-10 | OK |  |
| M_leiter_vflach | Leiterschritt δ* = 3Δs/r* = (v_flach/c)² (Milchstraße) | 3.29296e-07 | 3.29296e-07 | 1.1e-16 | OK |  |
| M_rstern_r0 | r*/r₀ = 1/(2δ*) (Milchstraße) | 1.51839e+06 | 1.51839e+06 | 1.1e-16 | OK |  |
| M_gw_schwelle | Schwellenfrequenz der Wellen f*_gw·π = Ω_N(r*) (Milchstraße) | 6.39409e-16 | 6.39409e-16 | 0.0e+00 | OK |  |
| M_diffusion_r0 | Diffusion der Rückkehr D = r*²/τ_D = r₀c/2 = GM/c (Sonne) | 4.42681e+11 | 4.42681e+11 | 0.0e+00 | OK |  |
| M_zelle_hubble | Zelle in Hubble-Längen L·H₀/c = (8π²/(ν_kosm Ω_b Θ))^{1/3} (L = d_LSS, a₀ aus der Zelle) | 6.24689 | 6.24689 | 1.1e-16 | OK | [G: OMEGA_R] |
| M_zelle_kuenftig | Durchmesser der letzten Streufläche im Endzustand (Beobachter bei a → ∞): L/d_LSS fällt auf seinen Grenzwert | 0.725998 | 0.725998 | 2.8e-10 | OK | [G: OMEGA_R] |
| M_leiter_ebenen | Leiter r*² = 3Δs·ℓ₀ auf jeder Ebene (Proton, Sonne, Milchstraße, Haufen): r* ist das geometrische Mittel von Faltung und Wirkungslänge | 1.01487e-27 | 1.01487e-27 | 2.2e-16 | OK |  |
| M_bruecke_LI_r0_lambda | Brücke L–I: r*² = r₀·Λ (Raum und innere Richtung) | 7.23886e+40 | 7.23886e+40 | 0.0e+00 | OK |  |
| M_bruecke_LI_gamma | Brücke L–I: Γ = a_N/c = D·k² (Glättungsrate mit D = GM/c, k = 1/r; Sonne bei 1 AE) | 1.97806e-11 | 1.97806e-11 | 0.0e+00 | OK |  |
| M_groessen_rollen | Jede Größe des Modells trägt höchstens ein Glied der Triade und höchstens eine Rolle im Zustandsraum; Zeichen stehen im Symbolregister | 1 | 1 | 0.0e+00 | OK |  |
| M_grenzlaengen_ebenen | Grenzlängen jeder Ebene 3Δs, r*, ℓ₀: Faltung und Schwellenradius im Raum (L), Wirkungslänge in der inneren Richtung (I); ihre Schritte sind gleich (δ* = 3Δs/r* = r*/ℓ₀) | 3.89903e-41 | 3.89903e-41 | 2.2e-16 | OK |  |
| M_bruecke_IG_lambda | Brücke I–G: Λ = ν_kosm·ΣΔs_b der Zelle, ν_kosm = 1 + Q | 4.08525e+26 | 4.08525e+26 | 1.1e-14 | OK | [G: Ω_b] |
| M_bruecke_LG_RQ | Brücke L–G: bei R_Q = r*/(−ln q_kosm) erreicht das örtliche q das kosmische | 0.842874 | 0.842874 | 0.0e+00 | OK |  |
| M_bruecke_ZL_regel | Brücke Z–L: Der Zustandsraum wirkt auf den Ort über die Regel, L_eff = 𝒫L (𝒫 vertauscht mit dem Erzeuger, löscht keine Freiheitsgrade) | 1 | 1 | ‖[𝒫, L]‖ = 8.4e-17 | OK |  |
| M_bruecke_ZG_kac | Brücke Z–G: Der Zustandsraum bindet sich ans Ganze über das Kac-Budget, Σ_i 1/⟨T_{A_i}⟩ = 1 für jede Zerlegung des geschlossenen Raums | 1 | 1 | 0.0e+00 | OK |  |
| M_budget_ladung | Kein Außen: auf dem geschlossenen Raum verlangt Gauß Σ Quelle = 0 (Σ Q_el = 0) | 1 | 1 | Rest 1e-15 / 0.10 | OK |  |
| M_budget_rueckstellung | Rückstellung σ = ΣΔs − ΔL = 0: Dehnung der Raumschleifen je Radius = Faltung der Baryonen je Radius | 0.0833401 | 0.0833401 | 0.0e+00 | OK | [G: OMEGA_R, Ω_b] |
| M_sattelpunkt_newton | Sattelpunkt Sonne–Erde ohne Zusatzbeschleunigung: er liegt in den Glättungszonen beider Körper | 1 | 1 | Sattel 2.59e+08 m, ℓ_Erde 2778 AE, ℓ_Sonne 192555 AE | OK | [G: RHO_B_SONNE, V0_SONNE] |
| M_glaettung_schwelle | Glättung: ℓ³ = GM/𝒯 und ℓ/r* = (M×/M)^{1/6} mit M× = a₀³/(G𝒯²) (Sonne im Feld der Galaxie) | 26.2253 | 26.2253 | 1.8e-15 | OK | [G: RHO_B_SONNE, V0_SONNE] |
| M_feldform_impuls | Feldform V: Impuls erhalten (Gesamtkraft null, ∂_i𝒮_ij = 0 im Newton-Grenzfall); die örtliche Regel allein nicht | 1 | 1 | Feldform 9e-05, örtliche Regel 2e-01 | OK |  |
| M_feldform_energie | Feldform V: Energie erhalten, die Kraft ist ein Gradient (Arbeit auf geschlossenem Weg null) | 1 | 1 | Umlaufarbeit / (¦F¦·Umfang) 2e-05 | OK |  |

### Klasse E

| Gesetz | Aussage | Modell | Vergleich | Maß | Ergebnis | [G] / [W] |
|---|---|---|---|---|---|---|
| E_a0_gegen_relation | Schwelle a₀ gegen die gemessene Relation g† (Fehler stat und sys) | 1.1e-10 | 1.2e-10 | -0.4σ | OK | [G: A0_LIT] |
| E_zelle_ohne_kreispaare | Zelle ohne Kreispaare: a₀ ≤ a₀(L = 0,985·d_LSS) (kubischer Torus) | 1.1e-10 | 1.12705e-10 | Abstand +2.4e-02 | OK | [G: KREISPAAR_GRENZE, OMEGA_R] |
| E_n_eff | N_eff mit dem Vakuumteil des Modells bei der Rekombination (a₀(z), H(z)): früh fast null | 3.046 | 2.99 | +0.3σ | OK | [G: N_EFF_MESS, N_EFF_STD, OMEGA_R] |
| E_omega_k | Krümmung Ω_K der flachen Zelle (Zelle atmet, Λ fest) gegen die Messung | 0 | 0.0007 | -0.4σ | OK | [G: OMEGA_K_MESS, OMEGA_R] |
| E_crater2 | Crater II σ_v, örtliche Regel (Jeans, Außenfeld der Milchstraße) | 1.38699 | 2.7 | -4.4σ | WIDERSPRUCH | [G: CRATER2_ABSTAND, CRATER2_M, CRATER2_RE] |
| E_gpb_geodaetisch | Gravity Probe B geodätische Drehung [mas/Jahr] | 6603.88 | 6602 | +0.1σ | OK | [W: GPB_A] |
| E_gpb_mitfuehrung | Gravity Probe B Mitführung [mas/Jahr] (Bahnmittel längs der Achse) | 40.8003 | 37.2 | +0.5σ | OK | [G: ERDE_J] [W: GPB_A] |
| E_hulse_taylor | Hulse-Taylor: gemessene zu erwarteter Abnahme der Umlaufzeit (Modell gleich ART) | 1 | 0.998 | +1.0σ | OK | [G: HT_VERHAELTNIS, PSR1913_E, PSR1913_M1, PSR1913_M2, PSR1913_PB] |
| E_schatten_m87 | Schatten M87* (Masse und Abstand EHT) [μas] | 39.6881 | 42 | -0.8σ | OK |  |

## Gegenproben

Wächter gegenproben_abdeckung: 11 von 260 Modellformeln ohne unabhängige Probe (davon 7 nur mit Umkehrprobe, 3 nur mit Proben, die eine Formel mit der geprüften teilen); Umkehrproben (art "umkehr") und teilende Proben zählen nicht.

| Gegenprobe | prüft | Art | größter Abstand | ruft sich selbst | Ergebnis |
|---|---|---|---|---|---|
| v_bar_quadrat_handwert | v_bar_quadrat | handwert | 0.0e+00 | nein | OK |
| g_gemessen_ohne_druck | g_gemessen | grenzfall | 0.0e+00 | nein | OK |
| g_gemessen_handwert | g_gemessen | handwert | 0.0e+00 | nein | OK |
| rest_galaxie_synthetisch | rest_galaxie | grenzfall | 1.3e-15 | nein | OK |
| a0_je_galaxie_synthetisch | a0_je_galaxie | grenzfall | 8.8e-13 | nein | OK |
| streuung_messfehler_abstand | streuung_messfehler | handwert | 4.1e-04 | nein | OK |
| toomre_flache_kurve | toomre_sigma_dyn | handwert | 7.9e-08 | nein | OK |
| g_expscheibe_fernfeld | g_expscheibe | grenzfall | 1.1e-04 | nein | OK |
| sigma_dyn_expscheibe_newton | sigma_dyn_expscheibe | grenzfall | 1.9e-07 | nein | OK |
| kuzmin_feld_potential | kuzmin_feld | zweiter_weg | 9.1e-08 | nein | OK |
| g_miyamoto_nagai_potential | g_miyamoto_nagai | zweiter_weg | 3.5e-09 | nein | OK |
| a_scheibe_ebene_radial | a_scheibe_ebene | grenzfall | 0.0e+00 | nein | OK |
| a_scheibe_ebene_kuzmin_tief | a_scheibe_ebene | grenzfall | 5.3e-04 | nein | OK |
| vertikalfrequenz_feld | vertikalfrequenz | zweiter_weg | 7.3e-12 | nein | OK |
| knotenpraezession_analytisch | knotenpraezession_klein | zweiter_weg | 3.7e-08 | nein | OK |
| praezession_ring_klein | praezession_ring | grenzfall | 9.4e-05 | nein | OK |
| einseitigkeit_m1_quadratur | einseitigkeit_m1 | zweiter_weg | 3.2e-08 | nein | OK |
| massenmittel_A_plummer_schwach | massenmittel_A_plummer | grenzfall | 2.4e-09 | nein | OK |
| massenmittel_A_scheibe_schwach | massenmittel_A_scheibe | grenzfall | 1.1e-12 | nein | OK |
| sigma_plummer_jeans_newton | sigma_plummer_jeans | zweiter_weg | 8.9e-05 | nein | OK |
| sigma_v_tief_handwert | sigma_v_tief | handwert | 1.6e-03 | nein | OK |
| gezeit_betrag_handwert | gezeit_betrag | handwert | 2.2e-16 | nein | OK |
| gezeit_sonnenort_newton | gezeit_sonnenort | grenzfall | 0.0e+00 | nein | OK |
| gezeit_hill_frequenzen | gezeit_hill | zweiter_weg | 2.2e-16 | nein | OK |
| glaettungslaenge_frequenz | glaettungslaenge | zweiter_weg | 2.0e-15 | nein | OK |
| gezeitendichte_mittel | gezeitendichte | umkehr | 6.7e-15 | nein | OK |
| schnittmasse_glaettung_schnitt | schnittmasse_glaettung | zweiter_weg | 1.4e-14 | nein | OK |
| schnittradius_glaettung_schnitt | schnittradius_glaettung | zweiter_weg | 7.0e-15 | nein | OK |
| fq_form_handwert | fq_form | handwert | 6.1e-15 | nein | OK |
| a0_an_punkten_synthetisch | a0_an_punkten | grenzfall | 5.2e-07 | nein | OK |
| nu_ds_bose | nu_ds | grenzfall | 0.0e+00 | nein | OK |
| ds_anpassung_synthetisch | ds_anpassung | grenzfall | 1.7e-10 | nein | OK |
| spektraldimension_steigung | spektraldimension_regel | zweiter_weg | 1.6e-10 | nein | OK |
| flussdimension_handwert | flussdimension | handwert | 2.2e-16 | nein | OK |
| dichte_gelesen_divergenz | dichte_gelesen | zweiter_weg | 4.0e-11 | nein | OK |
| r_stern_filament_schnitt | r_stern_filament | zweiter_weg | 2.2e-16 | nein | OK |
| r_stern_filament_handwert | r_stern_filament | handwert | 1.4e-02 | nein | OK |
| faltungsmesser_masse | faltungsmesser | zweiter_weg | 1.1e-16 | nein | OK |
| ds_eff_beschleunigung | ds_eff | zweiter_weg | 1.1e-16 | nein | OK |
| masse_aus_temperatur_kreis | masse_aus_temperatur | umkehr | 2.2e-16 | nein | OK |
| pulsar_rate_hand | pulsar_rate | handwert | 2.9e-07 | nein | OK |
| warp_viertelzeit_hand | warp_viertelzeit | handwert | 1.7e-08 | nein | OK |
| verstaerkung_druck_hand | verstaerkung_druck | handwert | 0.0e+00 | nein | OK |
| a_regel_hand | a_regel | handwert | 4.4e-07 | nein | OK |
| v_kreis_punkt_hand | v_kreis_punkt | zweiter_weg | 0.0e+00 | nein | OK |
| v_flach_schnitt | v_flach | zweiter_weg | 1.7e-15 | nein | OK |
| nu_aussenfeld_hand | nu_aussenfeld | zweiter_weg | 0.0e+00 | nein | OK |
| sigma_kollektiv_hand | sigma_kollektiv | handwert | 2.9e-05 | nein | OK |
| eta_scheibe_kuzmin | eta_scheibe | grenzfall | 3.8e-04 | nein | OK |
| feldform_naeherung_kuzmin | feldform_naeherung | grenzfall | 2.5e-04 | nein | OK |
| gleichstand_scheibe_nullstelle | gleichstand_scheibe | umkehr | 4.4e-16 | nein | OK |
| saeule_loesung_kuzmin | saeule_loesung | grenzfall | 7.4e-04 | nein | OK |
| a_N_aus_kreisbahn_vorwaerts | a_N_aus_kreisbahn | zweiter_weg | 2.2e-16 | nein | OK |
| sattelpunkt_nullstelle | sattelpunkt | zweiter_weg | 2.2e-16 | nein | OK |
| dynamikmasse_ort_hand | dynamikmasse_ort | zweiter_weg | 2.2e-16 | nein | OK |
| knotenpraezession_frequenzen_hand | knotenpraezession_frequenzen | handwert | 2.2e-16 | nein | OK |
| knotenpraezession_scheibe_bessel | knotenpraezession_scheibe | zweiter_weg | 6.7e-09 | nein | OK |
| warp_aufwicklung_hand | warp_aufwicklung | handwert | 0.0e+00 | nein | OK |
| a0_profiliert_synthetisch | a0_profiliert | grenzfall | 1.9e-05 | nein | OK |
| profil_chi2_least_squares | profil_chi2 | zweiter_weg | 3.1e-13 | nein | OK |
| a0_je_drittel_synthetisch | a0_je_drittel | grenzfall | 4.9e-04 | nein | OK |
| a0_senkung_ohne_feld | a0_senkung_aussenfeld | grenzfall | 1.3e-13 | nein | OK |
| sigma_dyn_rest_synthetisch | sigma_dyn_rest | grenzfall | 7.1e-15 | nein | OK |
| rest_gasanteil_synthetisch | rest_gasanteil | zweiter_weg | 3.1e-13 | nein | OK |
| rest_jenseits_synthetisch | rest_jenseits | grenzfall | 2.9e-15 | nein | OK |
| reststeigung_synthetisch | reststeigung_aussen | grenzfall | 4.4e-16 | nein | OK |
| ds_streuung_eine_ziehung | ds_streuung | grenzfall | 0.0e+00 | nein | OK |
| boden_grenze_ohne_fehler | boden_grenze | grenzfall | 1.3e-14 | nein | OK |
| boden_anpassung_synthetisch | boden_anpassung | grenzfall | 1.7e-14 | nein | OK |
| boden_drift_hand | boden_drift | handwert | 0.0e+00 | nein | OK |
| boden_quer_hand | boden_quer | handwert | 4.9e-08 | nein | OK |
| gewicht_f_synthetisch | gewicht_f | grenzfall | 1.0e-03 | nein | OK |
| quadrupol_glaettung_simpson | quadrupol_glaettung | zweiter_weg | 4.6e-06 | nein | OK |
| omega_vakuum_theta_eins | omega_vakuum | grenzfall | 5.6e-16 | nein | OK |
| lambda_vakuum_dichte_hand | lambda_vakuum | handwert | 9.1e-05 | nein | OK |
| lambda_vakuum_dichte | lambda_vakuum | zweiter_weg | 2.2e-16 | nein | OK |
| H_vakuum_H0_wurzel_omega | H_vakuum | zweiter_weg | 2.2e-16 | nein | OK |
| H_vakuum_33 | H_vakuum | handwert | 4.6e-04 | nein | OK |
| entropie_vakuum_hand | entropie_vakuum | zweiter_weg | 2.2e-16 | nein | OK |
| horizontentropie_faltung | horizontentropie | zweiter_weg | 2.2e-16 | nein | OK |
| zustand_vakuum_erhaltung | zustand_vakuum | zweiter_weg | 2.3e-12 | nein | OK |
| hintergrund_modell_de_sitter | hintergrund_modell | grenzfall | 0.0e+00 | nein | OK |
| hintergrund_modell_hand | hintergrund_modell | handwert | 8.6e-05 | nein | OK |
| atemrate_eff_ohne_flucht | atemrate_eff | grenzfall | 0.0e+00 | nein | OK |
| atemrate_eff_hand | atemrate_eff | handwert | 0.0e+00 | nein | OK |
| d_LSS_beobachter_gitter | d_LSS_beobachter | zweiter_weg | 8.1e-09 | nein | OK |
| teilchenhorizont_eds | teilchenhorizont | grenzfall | 2.1e-09 | nein | OK |
| rueckstellbilanz_faltung | rueckstellbilanz | zweiter_weg | 2.2e-16 | nein | OK |
| faltung_horizont_hubble_masse | faltung_horizont | zweiter_weg | 1.1e-16 | nein | OK |
| theta_temperaturen | theta_ | zweiter_weg | 2.2e-16 | nein | OK |
| a0_modell_theta | a0_modell | zweiter_weg | 2.2e-16 | nein | OK |
| gewicht_hintergrund_rueckkehr | gewicht_hintergrund | zweiter_weg | 0.0e+00 | nein | OK |
| zelle_hubble_zelle | zelle_hubble | zweiter_weg | 3.3e-15 | nein | OK |
| a0_aus_zelle_faltung | a0_aus_zelle | zweiter_weg | 0.0e+00 | nein | OK |
| zelle_aus_a0_umkehr | zelle_aus_a0 | umkehr | 3.3e-15 | nein | OK |
| zellgrenze_a0_zelle | zellgrenze | zweiter_weg | 1.1e-16 | nein | OK |
| zellgrenze_hand | zellgrenze | handwert | 3.9e-04 | nein | OK |
| a0_frueh_relativ_hand | a0_frueh_relativ | handwert | 5.8e-05 | nein | OK |
| delta_N_eff_hand | delta_N_eff | handwert | 1.4e-03 | nein | OK |
| r_Q_milchstrasse | r_Q | handwert | 1.8e-04 | nein | OK |
| masse_bQ_naht | masse_bQ | zweiter_weg | 2.6e-14 | nein | OK |
| dichte_umkehr_heute_eds | dichte_umkehr_heute | grenzfall | 1.3e-09 | nein | OK |
| umkehrradius_heute_trichterkette | umkehrradius_heute | zweiter_weg | 9.1e-06 | nein | OK |
| schalen_dauer_zykloide | schalen_dauer | grenzfall | 5.1e-03 | nein | OK |
| fluchtmass_bereich_wirkung | fluchtmass_bereich | handwert | 3.5e-06 | nein | OK |
| kreislauf_mindestrate_hand | kreislauf_mindestrate | zweiter_weg | 2.2e-16 | nein | OK |
| q_fenster_hand | q_fenster | zweiter_weg | 2.2e-16 | nein | OK |
| uebergabezeit_dichte | uebergabezeit | zweiter_weg | 2.2e-16 | nein | OK |
| haufenstoss_trennzeit_umlauf | haufenstoss_trennzeit | grenzfall | 0.0e+00 | nein | OK |
| omega_de_hand | omega_de | handwert | 0.0e+00 | nein | OK |
| weltzeit_materie_lambda | weltzeit | grenzfall | 2.2e-16 | nein | OK |
| weltzeit_materie_strahlung | weltzeit | grenzfall | 1.2e-15 | nein | OK |
| weltzeit_gegebene_atemrate | weltzeit | grenzfall | 4.4e-16 | nein | OK |
| horizontwachstum_ableitung | horizontwachstum | zweiter_weg | 1.5e-08 | nein | OK |
| kreislauf_versiegen_fern | kreislauf_versiegen | grenzfall | 2.2e-16 | nein | OK |
| formfaktor_zelle_hand | formfaktor_zelle | handwert | 4.4e-16 | nein | OK |
| drehwinkel_zelle_bilder | drehwinkel_zelle | zweiter_weg | 0.0e+00 | nein | OK |
| streuung_kosmisch_chi2 | streuung_kosmisch | zweiter_weg | 0.0e+00 | nein | OK |
| anteil_isotrop_hand | anteil_isotrop | handwert | 0.0e+00 | nein | OK |
| gewicht_schale_aussen_hintergrund | gewicht_schale_aussen | grenzfall | 1.0e-08 | nein | OK |
| gewicht_schale_aussen_tafel | gewicht_schale_aussen | handwert | 4.1e-06 | nein | OK |
| rand_schale_aussen_tafel | rand_schale_aussen | handwert | 4.2e-04 | nein | OK |
| uebergabe_vor_umkehr_bahn | uebergabe_vor_umkehr | handwert | 3.1e-04 | nein | OK |
| uebergabe_faltung_archiv | uebergabe_faltung | handwert | 9.3e-04 | nein | OK |
| wachstum_uebergabe_lineare_antwort | wachstum_uebergabe | grenzfall | 6.1e-05 | nein | OK |
| wachstum_uebergabe_archiv | wachstum_uebergabe | handwert | 4.0e-03 | nein | OK |
| gang_schwelle_strahlung_differenz | gang_schwelle_strahlung | zweiter_weg | 2.5e-09 | nein | OK |
| neigung_stoerung_hand | neigung_stoerung | handwert | 6.6e-06 | nein | OK |
| linse_punktmasse_plummer | linse_punktmasse | grenzfall | 2.0e-12 | nein | OK |
| haufen_rest_bose | haufen_rest | zweiter_weg | 2.2e-16 | nein | OK |
| haufen_baryonen_quad | haufen_baryonen | zweiter_weg | 6.7e-16 | nein | OK |
| masse_linse_nfw_quad | masse_linse_nfw | zweiter_weg | 1.8e-15 | nein | OK |
| verlust_gamma_taktgefaelle | verlust_gamma | zweiter_weg | 1.0e-08 | nein | OK |
| tau_D_rueckkehrtiefe | tau_D | zweiter_weg | 0.0e+00 | nein | OK |
| rueckkehr_exponent_minimum | rueckkehr_exponent | zweiter_weg | 2.2e-16 | nein | OK |
| rueckkehr_laplace | rueckkehr | zweiter_weg | 2.2e-16 | nein | OK |
| s_punktmasse | s_ | zweiter_weg | 2.2e-16 | nein | OK |
| q_bose | q_ | zweiter_weg | 0.0e+00 | nein | OK |
| nu_bose | nu | zweiter_weg | 2.2e-16 | nein | OK |
| nu_grenzfall_tief | nu | grenzfall | 5.0e-07 | nein | OK |
| nu_grenzfall_newton | nu | grenzfall | 0.0e+00 | nein | OK |
| chi_nome | chi | zweiter_weg | 2.2e-16 | nein | OK |
| chi_elastizitaet | chi | zweiter_weg | 3.4e-10 | nein | OK |
| n_B_rueckkehrsumme | n_B | zweiter_weg | 2.2e-16 | nein | OK |
| Q_tilde_vakuum | Q_tilde | handwert | 4.4e-16 | nein | OK |
| Q_tilde_gammareihe | Q_tilde | zweiter_weg | 2.2e-16 | nein | OK |
| schliesslaenge_kac | schliesslaenge | zweiter_weg | 2.2e-16 | nein | OK |
| ell0_leiter | ell0 | umkehr | 0.0e+00 | nein | OK |
| a0_aus_schliesslaenge_kac | a0_aus_schliesslaenge | zweiter_weg | 0.0e+00 | nein | OK |
| gewicht_w_laplace | gewicht_w | zweiter_weg | 2.2e-16 | nein | OK |
| gewicht_w_halb | gewicht_w | handwert | 0.0e+00 | nein | OK |
| gewicht_w_preisach_minimum | gewicht_w_preisach | zweiter_weg | 0.0e+00 | nein | OK |
| kern_zeit_ueberleben | kern_zeit | zweiter_weg | 7.9e-10 | nein | OK |
| kern_bahnfrequenz_wurzel | kern_bahnfrequenz | zweiter_weg | 1.3e-15 | nein | OK |
| glaettung_innen_randwert | glaettung_innen | zweiter_weg | 8.3e-06 | nein | OK |
| erneuerung_randwert | erneuerung | zweiter_weg | 3.1e-08 | nein | OK |
| eigenwert_schluss_matrix | eigenwert_schluss | zweiter_weg | 2.2e-08 | nein | OK |
| rueckkehr_geschlossen_matsubara | rueckkehr_geschlossen | zweiter_weg | 2.9e-12 | nein | OK |
| rueckkehr_geschlossen_irrweg | rueckkehr_geschlossen | handwert | 4.6e-07 | nein | OK |
| abweichung_schluss_grenzfall | abweichung_schluss | grenzfall | 6.6e-07 | nein | OK |
| abweichung_schluss_modelltext | abweichung_schluss | handwert | 3.7e-03 | nein | OK |
| abweichung_schluss_rc_kette | abweichung_schluss | zweiter_weg | 4.0e-04 | nein | OK |
| radius_schluss_eigenwert | radius_schluss | zweiter_weg | 2.2e-15 | nein | OK |
| temperatur_innen_periode | temperatur_innen | zweiter_weg | 0.0e+00 | nein | OK |
| temperatur_innen_handwert | temperatur_innen | handwert | 1.3e-04 | nein | OK |
| kac_zeit_schluss_markow | kac_zeit_schluss | zweiter_weg | 2.2e-12 | nein | OK |
| strukturteil_Q_tilde | strukturteil | zweiter_weg | 1.5e-13 | nein | OK |
| neigung_xy_newton | neigung_xy | zweiter_weg | 1.0e-03 | nein | OK |
| ds_leiter | ds | umkehr | 2.2e-16 | nein | OK |
| neigung_faltung | neigung | zweiter_weg | 0.0e+00 | nein | OK |
| takt_metrik | takt | zweiter_weg | 0.0e+00 | nein | OK |
| delta_metrik | delta | zweiter_weg | 2.2e-16 | nein | OK |
| delta_1fm_hand | delta | handwert | 1.3e-05 | nein | OK |
| r0_schwarzschild | r0 | zweiter_weg | 0.0e+00 | nein | OK |
| r_stern_leiterweg | r_stern | umkehr | 0.0e+00 | nein | OK |
| r_stern_leiter_direkt | r_stern_leiter | zweiter_weg | 0.0e+00 | nein | OK |
| leiter_schritt_innen | leiter_schritt | zweiter_weg | 2.2e-16 | nein | OK |
| leiter_schritt_flach | leiter_schritt | zweiter_weg | 1.1e-16 | nein | OK |
| gezeit_betrag_eigenwerte | gezeit_betrag | zweiter_weg | 2.2e-16 | nein | OK |
| glaettungslaenge_dichte | glaettungslaenge | umkehr | 0.0e+00 | nein | OK |
| gezeitendichte_wolke | gezeitendichte | umkehr | 6.4e-15 | nein | OK |
| schnittmasse_glaettung_gleichheit | schnittmasse_glaettung | zweiter_weg | 1.4e-14 | nein | OK |
| schnittradius_glaettung_gleichheit | schnittradius_glaettung | zweiter_weg | 7.2e-15 | nein | OK |
| glaettung_zu_schwelle_direkt | glaettung_zu_schwelle | zweiter_weg | 2.3e-15 | nein | OK |
| glaettung_zu_schwelle_modelltext | glaettung_zu_schwelle | handwert | 9.4e-04 | nein | OK |
| flussdimension_ableitung | flussdimension | zweiter_weg | 8.9e-12 | nein | OK |
| epizyklen_summe_bahn | epizyklen_summe | zweiter_weg | 5.7e-10 | nein | OK |
| realteil_quadrat_punkt_epizykel | realteil_quadrat_punkt | zweiter_weg | 9.1e-10 | nein | OK |
| realteil_quadrat_punkt_triade | realteil_quadrat_punkt | handwert | 5.5e-06 | nein | OK |
| radius_resonanz_bose | radius_resonanz | zweiter_weg | 2.9e-13 | nein | OK |
| radius_resonanz_modelltext | radius_resonanz | handwert | 1.1e-02 | nein | OK |
| radius_resonanz_modelltext_aussen | radius_resonanz | handwert | 2.7e-02 | nein | OK |
| im_tau_aus_nome_umkehr | im_tau_aus_nome | umkehr | 2.2e-16 | nein | OK |
| kac_rueckkehrzeit_drehung | kac_rueckkehrzeit | zweiter_weg | 1.1e-05 | nein | OK |
| kac_budget_markow | kac_budget | zweiter_weg | 1.0e-04 | nein | OK |
| kac_budget_wicklung_rueckwaerts | kac_budget | zweiter_weg | 1.4e-06 | nein | OK |
| bose_leiter_summe | bose_leiter | zweiter_weg | 4.4e-16 | nein | OK |
| bose_leiter_integral_quad | bose_leiter_integral | zweiter_weg | 2.2e-16 | nein | OK |
| bose_leiter_integral_flussmass | bose_leiter_integral | zweiter_weg | 5.7e-14 | nein | OK |
| fisher_metrik_ableitung | fisher_metrik | zweiter_weg | 2.0e-10 | nein | OK |
| fisher_metrik_sinh | fisher_metrik | zweiter_weg | 1.8e-15 | nein | OK |
| entropie_ruecklauf_shannon | entropie_ruecklauf | zweiter_weg | 1.3e-15 | nein | OK |
| schleifensumme_rueckläufe | schleifensumme | zweiter_weg | 6.7e-15 | nein | OK |
| kipp_schleife_normalform | kipp_schleife | grenzfall | 1.2e-06 | nein | OK |
| rastbreite_schwebung | rastbreite | zweiter_weg | 4.4e-10 | nein | OK |
| preisach_ausgang_volle_folge | preisach_ausgang | zweiter_weg | 0.0e+00 | nein | OK |
| rc_kette_nu_netz | rc_kette_nu | zweiter_weg | 7.9e-13 | nein | OK |
| rc_kette_nu_modelltext | rc_kette_nu | handwert | 2.7e-05 | nein | OK |
| eisenstein_G_gittersumme | eisenstein_G | zweiter_weg | 6.6e-07 | nein | OK |
| eisenstein_G_lemniskate | eisenstein_G | handwert | 8.9e-16 | nein | OK |
| nu_aus_modul_nome | nu_aus_modul | umkehr | 0.0e+00 | nein | OK |
| nu_aus_modul_modelltext | nu_aus_modul | handwert | 1.6e-04 | nein | OK |
| nu_aus_modul_modelltext_spitze | nu_aus_modul | handwert | 1.3e-08 | nein | OK |
| kupfer_hopf_abstand_routh | kupfer_hopf_abstand | zweiter_weg | 1.9e-06 | nein | OK |
| kupfer_ringfrequenz_routh | kupfer_ringfrequenz | zweiter_weg | 1.1e-06 | nein | OK |
| kupfer_ringfrequenz_alt | kupfer_ringfrequenz | handwert | 3.0e-04 | nein | OK |
| kupfer_zeitverhaeltnis_mode | kupfer_zeitverhaeltnis | zweiter_weg | 2.2e-16 | nein | OK |
| kupfer_zeitverhaeltnis_modelltext | kupfer_zeitverhaeltnis | handwert | 3.2e-03 | nein | OK |
| bose_stufe_max_lambert | bose_stufe_max | zweiter_weg | 3.3e-16 | nein | OK |
| kupfer_lyapunov_symbolisch | kupfer_lyapunov | zweiter_weg | 6.9e-07 | nein | OK |
| draht_schleife_symbolisch | draht_schleife | zweiter_weg | 4.6e-07 | nein | OK |
| delta_reihe_eisenstein | delta_reihe | zweiter_weg | 0.0e+00 | nein | OK |
| delta_modul_gespiegelt | delta_modul | zweiter_weg | 6.7e-15 | nein | OK |
| resolvente_flucht_fibonacci | resolvente_flucht | zweiter_weg | 1.8e-10 | nein | OK |
| spiel_umkehrpunkt_bahn | spiel_umkehrpunkt | zweiter_weg | 1.2e-12 | nein | OK |
| spiel_umkehrpunkt_modelltext | spiel_umkehrpunkt | handwert | 2.6e-04 | nein | OK |
| rate_rastpunkt_ableitung | rate_rastpunkt | zweiter_weg | 8.7e-12 | nein | OK |
| vo2_hopf_hurwitz | vo2_hopf | zweiter_weg | 3.0e-09 | nein | OK |
| schleife_selbsterregt_bahn | schleife_selbsterregt | zweiter_weg | 1.8e-02 | nein | OK |
| kuramoto_einrasten_klein | kuramoto_einrasten | zweiter_weg | 0.0e+00 | nein | OK |
| rate_quer_fluss | rate_quer | zweiter_weg | 2.1e-13 | nein | OK |
| watson_rueckkehr_gamma | watson_rueckkehr | zweiter_weg | 2.2e-16 | nein | OK |
| rueckkehrzeiten_drehung_je_start | rueckkehrzeiten_drehung | zweiter_weg | 0.0e+00 | nein | OK |
| chsh_oberschwingungen_frei | chsh_oberschwingungen | zweiter_weg | 2.2e-16 | nein | OK |
| chsh_oberschwingungen_grenzen | chsh_oberschwingungen | grenzfall | 2.2e-05 | nein | OK |
| kanal_verhaeltnis_ableitungen | kanal_verhaeltnis | zweiter_weg | 2.2e-10 | nein | OK |
| rueckgabe_rate_lebensdauer | rueckgabe_rate | zweiter_weg | 4.4e-16 | nein | OK |
| ds_sonne_handwert | ds | handwert | 1.7e-05 | nein | OK |
| r0_sonne_handwert | r0 | handwert | 1.1e-03 | nein | OK |
| compton_red_mittel | compton_red | zweiter_weg | 3.8e-12 | nein | OK |
| alpha_G_faltung | alpha_G | zweiter_weg | 2.2e-16 | nein | OK |
| sternmasse_faltung | sternmasse | zweiter_weg | 0.0e+00 | nein | OK |
| grenzmasse_hawking_temperatur | grenzmasse_hawking | zweiter_weg | 2.2e-16 | nein | OK |
| masse_r_stern_gleich_compton_nullstelle | masse_r_stern_gleich_compton | zweiter_weg | 7.0e-15 | nein | OK |
| eddington_thomson | eddington | zweiter_weg | 2.2e-16 | nein | OK |
| eddington_sonne_handwert | eddington | handwert | 5.2e-05 | nein | OK |
| lebensdauer_falte_page | lebensdauer_falte | zweiter_weg | 2.2e-16 | nein | OK |
| carter_faktor_sonne | carter_faktor | zweiter_weg | 1.0e-11 | nein | OK |
| carter_faktor_frequenzleiter | carter_faktor | zweiter_weg | 2.0e-11 | nein | OK |
| tiefe_faktor_potential | tiefe_faktor | zweiter_weg | 2.4e-11 | nein | OK |
| tiefe_faktor_analytisch | tiefe_faktor | grenzfall | 2.8e-10 | nein | OK |
| temperatur_kern_sonne_hand | temperatur_kern | handwert | 1.9e-05 | nein | OK |
| grenzmasse_zaehlregel_chandrasekhar | grenzmasse_zaehlregel | zweiter_weg | 5.3e-12 | nein | OK |
| grenzmasse_zaehlregel_handwert | grenzmasse_zaehlregel | handwert | 2.0e-04 | nein | OK |
| temperatur_max_langsam | temperatur_max | grenzfall | 1.0e-05 | nein | OK |
| decke_faktor_ableitung | decke_faktor | zweiter_weg | 5.6e-16 | nein | OK |
| leiterzahl_beschleunigung | leiterzahl | zweiter_weg | 0.0e+00 | nein | OK |
| leiterzahl_handwert | leiterzahl | handwert | 1.1e-02 | nein | OK |
| schnittmasse_unten_planck | schnittmasse_unten | zweiter_weg | 1.1e-16 | nein | OK |
| schnittmasse_oben_produkt | schnittmasse_oben | zweiter_weg | 2.2e-16 | nein | OK |
| schnittmasse_oben_handwert | schnittmasse_oben | handwert | 2.1e-04 | nein | OK |
| schwelle_zu_huelle_flaeche | schwelle_zu_huelle | zweiter_weg | 2.2e-16 | nein | OK |
| schwelle_zu_huelle_sonne_handwert | schwelle_zu_huelle | handwert | 1.3e-02 | nein | OK |
| spiegel_masse_faltung | spiegel_masse | zweiter_weg | 0.0e+00 | nein | OK |
| spiegel_radius_schliesslaenge | spiegel_radius | grenzfall | 0.0e+00 | nein | OK |
| verhaeltnis_em_grav_laengen | verhaeltnis_em_grav | zweiter_weg | 2.5e-12 | nein | OK |
| compton_takt_energie | compton_takt | zweiter_weg | 0.0e+00 | nein | OK |
| faltungsrate_mittel | faltungsrate | zweiter_weg | 3.3e-16 | nein | OK |
| faltungsrate_proton_handwert | faltungsrate | handwert | 1.0e-04 | nein | OK |
| faltungstakt_proton_handwert | faltungstakt | handwert | 1.5e-03 | nein | OK |
| faltungstakt_faltung | faltungstakt | zweiter_weg | 0.0e+00 | nein | OK |
| umlauf_wasserstoff_hartree | umlauf_wasserstoff | zweiter_weg | 5.0e-12 | nein | OK |
| umlauf_kepler_neutronenstern_handwert | umlauf_kepler | handwert | 3.0e-03 | nein | OK |
| umlauf_kepler_tempo | umlauf_kepler | zweiter_weg | 2.2e-16 | nein | OK |
| spinzahl_spin_halb | spinzahl | zweiter_weg | 2.2e-16 | nein | OK |
| spinzahl_elektron_handwert | spinzahl | handwert | 1.6e-02 | nein | OK |
| moment_magnetisch_bohr | moment_magnetisch | handwert | 1.1e-06 | nein | OK |
| einrast_phase_bewegung | einrast_phase | zweiter_weg | 1.8e-13 | nein | OK |
| dephasierzeit_koerper_handwert | dephasierzeit | handwert | 1.6e-03 | nein | OK |
| ablenkung_art | ablenkung | grenzfall | 6.3e-06 | nein | OK |
| ablenkung_sonnenrand_handwert | ablenkung | handwert | 6.8e-04 | nein | OK |
| ablenkung_kollektiv_rstern | ablenkung_kollektiv | zweiter_weg | 2.2e-16 | nein | OK |
| ablenkung_linse_sichtlinie | ablenkung_linse | zweiter_weg | 1.3e-15 | nein | OK |
| ablenkung_linse_I1_handwert | ablenkung_linse | handwert | 2.1e-04 | nein | OK |
| zusatzpotential_quadratur | zusatzpotential | zweiter_weg | 3.8e-14 | nein | OK |
| photonensphaere_bahn | photonensphaere | zweiter_weg | 2.2e-16 | nein | OK |
| schatten_winkel_bahn | schatten_winkel | zweiter_weg | 0.0e+00 | nein | OK |
| schatten_winkel_m87_handwert | schatten_winkel | handwert | 7.9e-03 | nein | OK |
| perihel_pro_umlauf_art | perihel_pro_umlauf | zweiter_weg | 2.2e-16 | nein | OK |
| perihel_merkur_handwert | perihel_pro_umlauf | handwert | 4.6e-04 | nein | OK |
| schwellenfrequenz_gw_umlauf | schwellenfrequenz_gw | zweiter_weg | 2.2e-16 | nein | OK |
| isco_kerr_nullstelle | isco_kerr | zweiter_weg | 5.6e-16 | nein | OK |
| hulse_taylor_peters | hulse_taylor | zweiter_weg | 8.9e-16 | nein | OK |
| hulse_taylor_handwert | hulse_taylor | handwert | 1.1e-03 | nein | OK |
| verwindung_flaeche | verwindung | zweiter_weg | 2.2e-16 | nein | OK |
| verwindung_erde_handwert | verwindung | handwert | 6.5e-04 | nein | OK |
| umlauf_differenz_frequenzen | umlauf_differenz | zweiter_weg | 1.1e-10 | nein | OK |
| umlauf_differenz_erde_handwert | umlauf_differenz | handwert | 5.0e-04 | nein | OK |
| kerr_omega_phi_kepler | kerr_omega_phi | grenzfall | 2.2e-16 | nein | OK |
| kerr_omega_r_schwarzschild | kerr_omega_r | grenzfall | 2.2e-16 | nein | OK |
| kerr_omega_theta_schwarzschild | kerr_omega_theta | grenzfall | 2.2e-16 | nein | OK |
| knotendrehung_frequenzen | knotendrehung | grenzfall | 6.7e-05 | nein | OK |
| nullflaeche_takt_aequator | nullflaeche_takt | grenzfall | 0.0e+00 | nein | OK |
| nullflaeche_takt_pol | nullflaeche_takt | grenzfall | 0.0e+00 | nein | OK |
| horizont_kerr_wurzeln | horizont_kerr | zweiter_weg | 1.8e-15 | nein | OK |
| mitfuehrung_dipol | mitfuehrung | zweiter_weg | 2.2e-16 | nein | OK |
| mitfuehrung_polarbahn_mittel | mitfuehrung_polarbahn | zweiter_weg | 4.4e-15 | nein | OK |
| geodaetisch_art | geodaetisch | zweiter_weg | 2.2e-16 | nein | OK |
| quadrupol_leistung_peters | quadrupol_leistung | zweiter_weg | 8.9e-16 | nein | OK |
| faltungsmesser_chirp_masse | faltungsmesser_chirp | zweiter_weg | 2.2e-16 | nein | OK |
| ringdown_frequenz_masse | ringdown_frequenz | zweiter_weg | 2.2e-16 | nein | OK |
| ringdown_quelle_handwert | ringdown_frequenz | handwert | 5.6e-05 | nein | OK |
| gedaechtnis_welle_energie | gedaechtnis_welle | zweiter_weg | 2.2e-16 | nein | OK |
| schleife_thermisch_sprung | schleife_thermisch | zweiter_weg | 5.6e-16 | nein | OK |
| ringdown_abklingzeit_daempfung | ringdown_abklingzeit | zweiter_weg | 2.2e-16 | nein | OK |
| richtungsfaktor_gedaechtnis_kosinus | richtungsfaktor_gedaechtnis | zweiter_weg | 6.7e-16 | nein | OK |
| richtungsfaktor_gedaechtnis_kante | richtungsfaktor_gedaechtnis | grenzfall | 0.0e+00 | nein | OK |
| entropie_falte_alpha_G | entropie_falte | zweiter_weg | 2.2e-16 | nein | OK |
| entropie_falte_smarr | entropie_falte | zweiter_weg | 4.4e-16 | nein | OK |
| page_zeit_gleichstand | page_zeit | zweiter_weg | 4.4e-16 | nein | OK |
| rueckverdichtung_teile_bilanz | rueckverdichtung_teile | zweiter_weg | 2.2e-16 | nein | OK |
| temperatur_abstrahlung_planck | temperatur_abstrahlung | zweiter_weg | 2.2e-16 | nein | OK |
| entropieexport_bose | entropieexport | zweiter_weg | 1.3e-15 | nein | OK |
| alpha_s_lauf_geschlossen | alpha_s_lauf | zweiter_weg | 3.4e-12 | nein | OK |
| skala_bei_kopplung_geschlossen | skala_bei_kopplung | zweiter_weg | 7.9e-11 | nein | OK |
| nome_qcd_pol | nome_qcd | zweiter_weg | 4.4e-16 | nein | OK |
| dirac_verhaeltnis_si | dirac_verhaeltnis | zweiter_weg | 7.9e-13 | nein | OK |
| radius_ladungskreis_kaluza_klein | radius_ladungskreis | zweiter_weg | 4.0e-13 | nein | OK |
| energie_ladungskreis_planck | energie_ladungskreis | zweiter_weg | 0.0e+00 | nein | OK |
| zustandsschwelle_maxwell | zustandsschwelle | zweiter_weg | 2.2e-16 | nein | OK |
| f_saha_hartree | f_saha | zweiter_weg | 1.0e-11 | nein | OK |
| eddington_quartik_druck | eddington_quartik | zweiter_weg | 5.3e-12 | nein | OK |
| eddington_quartik_chandrasekhar | eddington_quartik | handwert | 1.9e-03 | nein | OK |
| troepfchen_anpassung_normalgleichung | troepfchen_anpassung | zweiter_weg | 8.9e-13 | nein | OK |
| spaltbarkeit_energien | spaltbarkeit | zweiter_weg | 2.2e-16 | nein | OK |
| debye_laenge_bjerrum | debye_laenge | zweiter_weg | 2.2e-16 | nein | OK |
| debye_laenge_handwert | debye_laenge | handwert | 6.8e-04 | nein | OK |

## Themen

| Thema | Frage | Stand | Versuche (Urteil) |
|---|---|---|---|
| dunkle_energie | Was ist der zeitfeste Teil der Dunklen Energie neben dem Vakuumteil Ω_s (k_DE/k_s ≈ 2,9–3,1) und woher kommt seine Form w(z)? Hinweis des Autors: Tafel 5 und 6 („DM gegen B“) meinen Dunkle Energie; Invariante der Frage: Dunkle Energie als Bindungs- bzw. Schließungsstruktur statt Quelle; Trenntests je Kanal. | offen. Vakuumteil Ω_s im Modell. Vakuum je Richtung und je Zyklus widersprechen sich nicht: beide finden, dass kein Faktor 3 hergeleitet ist (die 3 steckt in π⁴/15); de_zyklen liest nur die Zahl k_DE/k_s ≈ 3 aus θ* [G], und SPARC schließt allein n = 1 aus. Offen, hängt an ω_b, ω_c [G]: die Mengenurteile der Fluchtmaß-Kette und modell_h0. | gegenlauf_dunkle_energie: offen; lambda_und_rueckstellung: offen; straffung_rueckstellung: offen; atmen_w: Teil im Modell; horizont_kreislauf: im Modell; de_ratenschwelle: offen, hängt an frei:Glättungsmasse; de_phasenseite: trägt nicht; kontraktion_atemrate: Teil im Modell; engstelle_rate: trägt nicht; engstelle_global: Teil im Modell; de_uebergabe_anteil: trägt nicht; de_bilanz_glieder: trägt nicht; vakuum_richtungen: trägt nicht; de_kehrseite_bindung: trägt nicht; de_fluchtmoden: trägt nicht; de_fluchtmoden_hierarchie: trägt nicht; de_fluchtmoden_verzoegert: trägt nicht; de_fluchtmoden_erneuerung: trägt nicht; de_fluchtmoden_mass: noch nicht im Modell; de_kanal_atmen: im Modell; de_rest_anstieg: offen; de_laufzeit: offen; modell_h0: offen, hängt an DE_SIGNIFIKANZ; de_invariante: noch nicht im Modell; regelkreise: noch nicht im Modell; de_zyklen: trägt nicht; maximum_schwelle: offen; lauf_eigenes_atmen: Teil im Modell; de_adiabatisch: noch nicht im Modell; de_pql_energie: noch nicht im Modell; de_fluchtmoden_faser: Teil im Modell; de_trenntests_kanal: noch nicht im Modell; de_torus_pythagoras: trägt nicht; de_topologie_ringe: noch nicht im Modell; de_nz_kac: noch nicht im Modell; de_kanal_s2: trägt nicht; de_gauss_bonnet_ganzes: trägt nicht; de_henkel: trägt nicht; de_buendel_verdrillung: trägt nicht; de_lambda_zaehlung: trägt nicht; de_einrasten: trägt nicht; de_nichtadditiv: offen, hängt an frei:Mindestmasse; de_fundamentalgruppe: trägt nicht; de_figur_lb: trägt nicht; de_spirale_flaechen: trägt nicht; de_treideln: trägt nicht |
| herkunft_q | Woher kommt Q = Ω_c/Ω_b ≈ 5,36 (ν_kosm = 1 + Q), das zeitfeste Gewicht der Mode ξ? | offen. Im Modell: Fenster 5 ≤ Q aus dem Schluss, Zählung je Radiant als offene Lesart. Offen, hängt an ω_b, ω_c [G]: q_form und lambda_bruecke (Kandidatenwahl). | herkunft_q: offen; camb_q: Teil im Modell; invarianten_ebenen: noch nicht im Modell; kac_schliessung: im Modell; lambda_bruecke: offen, hängt an frei:Kandidatenwahl; schluss_oben: Teil im Modell; modulraum_zahlen: trägt nicht; q_form: trägt nicht |
| a0_schwelle | Welchen Wert hat die Schwelle a₀ = c²/2Λ, in welcher Lesart, und was bindet sie an die Zelle? | im Modell: Standardwert 1,1·10⁻¹⁰, SPARC 1,06 ± 0,09 ± 0,04, Zelle ohne Kreispaare ≤ 1,13. Offen, hängt an N_eff [G]: die Urteile gegen die atmende Schwelle (die Zelle im frühen Kosmos steht im Modell, ihr Vergleich mit N_eff ist ein Gedächtniswert). | zeitgesetz_a0: offen, hängt an N_EFF_MESS; atmen_ohne_skala: offen, hängt an N_EFF_MESS; atmen_oder_zelle: offen, hängt an N_EFF_MESS; zelle_frueher_kosmos: im Modell; sparc_test: im Modell; schwelle_je_galaxie: im Modell; sparc_streuung: im Modell; a0_lesungen: im Modell; fruehe_galaxien_a0: noch nicht im Modell |
| zwerge | Warum weichen druckgestützte Systeme im äußeren Feld (DF2, DF4, Pal 14, NGC 2419, Fornax, Crater II, Draco) von der Regel ab? | offen, hängt an den Messwerten der Zwerge [G]: alle Urteile „trägt nicht“ dieser Reihe. aussenfeld_kosmisch trägt nicht gegen SPARC (eigene Daten); seine frei angepasste Zahl steht im Modell (Beschluss offen). | zwerge_aussenfeld: im Modell; tunnel_aussenfeld: offen; galaxie_als_ganzes: im Modell; sparc_aussenfeld_rand: im Modell; ort_oder_system: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA; zwerge_uebergabe: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA; aussenfeld_kosmisch: trägt nicht; gedaechtnis_umgebung: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA; ebenengrenze_zwerge: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA; halbgerade_hierarchie: offen; zwerge_gezeitentensor: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA; zwerge_fluchtmass_bahn: offen, hängt an DF2_SIGMA, DF4_SIGMA, DRACO_SIGMA, NGC2419_SIGMA |
| warp | Was dreht die Knotenlinie des Warps der Milchstraße rechtläufig (+10,9 km/s/kpc [G])? | offen, hängt an WARP_PRAEZESSION [G]; im Modell steht die Kinematik am Knoten und die rückläufigen Kippmoden. | warp_praezession: offen, hängt an WARP_PRAEZESSION; warp_drehmoment: Teil im Modell; ringe_oertlich: im Modell; warp_amplitude: Teil im Modell; warp_komplexe_kopplung: offen, hängt an WARP_PRAEZESSION |
| naht_uebergabe | Wie geht die Mode ξ an der w-Schale (Umkehrradius R_ta) in die örtliche Regel über? | Teil im Modell (I_F, Schale, Zuordnung der Ketten gesetzt); offen: ob die Übergabe λ = (1 − w)ω/N* die Mulde füllt und ob der Rest gebrochener Ketten ein Gedächtnis hat. | naht_rq_umkehr: noch nicht im Modell; mischregel_schale: noch nicht im Modell; schale_zeitlich: Teil im Modell; theta_bezug: im Modell; wirkung_luecken: Teil im Modell; wachstum_w_qs: im Modell; uebergabe_faltung: Teil im Modell; kollaps_nach_umkehr: Teil im Modell; stoerungen_rekombination: im Modell |
| linse_haufen | Trägt die Regel mit Übergabe die Linsen um isolierte Galaxien und das Profil der Haufen? | offen: Spannung zwischen der aus dem Modell bestimmten Übergabe und den Repo-Linsen jenseits 0,5 Mpc; Haufenprofil. | haufen_nach_naht: offen; linse_ohne_freiheit: Teil im Modell; linsen_isoliert: noch nicht im Modell; harvey_nachlauf: offen |
| mikroskala | Legt eine Gleichung des Modells m_p/m_e, α oder Im τ_X fest? | offen: keine Gleichung des Modells legt m_p/m_e oder α fest. | mikro_nome: Teil im Modell; hopf_knoten: noch nicht im Modell; em_phasenseite: im Modell; tau_ebenen: Teil im Modell; spiegel_m22: noch nicht im Modell; fusion_trichter: Teil im Modell |
| zellform | Welche Form hat die Zelle (Würfel, sechseckig, raum-, flächenzentriert)? | im Modell ohne festgelegte Form; hängt an Q [G] und dem Lauf des eigenen Atmens Ẋ_an(t). | zelle_form: im Modell; q_a0_kopplung: im Modell; winkelmoden_zelle: Teil im Modell |
| zustandsraum_labor | Welche Aussagen des Zustandsraums sind prüfbar, und gibt es eine Laborvorhersage nur des Modells? | Form im Modell, keine eigene Laborvorhersage. | labor_vorhersage: im Modell; kontakt_atmen: Teil im Modell; axiom_mass: im Modell; gesetzespruefung: Teil im Modell; generic_praezision: Teil im Modell; noether_form: noch nicht im Modell |
| zeit_gedaechtnis | Welche Zeitskala hat das Gedächtnis der Regel (Kern K, λ_K, λ_m = ω_orb/N*), und wie hängt die Zeit an der inneren Richtung? | Kern im Modell; N* offen (0,3–2,4). | kern_gedaechtnis: offen; kopfgleichung_kern: im Modell; gedaechtnis_topologie: Teil im Modell; innere_richtung_zeit: Teil im Modell |
| feldform_scheiben | Wirkt die Regel örtlich auf der Bahn (Kernform) oder als Feld (Feldform), und was zeigen Scheiben? | im Modell (örtliche Regel); offen die Vorhersagen in C. | sparc_boseform: offen; sparc_flaechendichte: im Modell; gleichung_ganzes: im Modell; trenntests_b_dm: Teil im Modell; sparc_universal: im Modell |
| rueckkehr_gleichung | Folgt die Regel ν = 1 + n_B(s) aus einer Gleichung mit wenigen Setzungen, und welche Invariante trägt sie? | im Modell; weitere Umformungen gesperrt. | gleichung_zeit: Teil im Modell; setzungen_gleichung: im Modell; eine_struktur: im Modell; info_geometrie: Teil im Modell; strittig_zweiter_weg: offen; schroedinger_m23: Teil im Modell; tunnelbild_gesetz: noch nicht im Modell; gemeinsame_innere_richtung: Teil im Modell; hannay_holonomie: noch nicht im Modell |
| trichter_hierarchie | Wie glättet die Hierarchie der Trichter (ρ > ρ_J) von Körper zu Galaxie zu Haufen? | Teil im Modell; Trichterkette offen, hängt an Massen [G]. | ebenen_projektion: Teil im Modell; trichterkette_gewicht: noch nicht im Modell; wep_glaettung: Teil im Modell; trichter_invarianten: Teil im Modell |

## Beschlüsse

| Beschluss | Aussage | Ort | Datum | Richtung |
|---|---|---|---|---|
| a0_standardwert_1_1 | Das Modell rechnet mit a₀ = 1,1·10⁻¹⁰ m/s² statt 1,2·10⁻¹⁰ | 00-bauplan.md: das Modell rechnet mit 1,1·10⁻¹⁰ | 2026-09-29 | hinein |
| gewicht_aus_rueckkehr | Gewicht der Schließungsregel w = 1 − e^{−√(θ̃τ_D)} aus derselben Rückkehr | 03-von-der-galaxie-bis-zum-kosmos.md: 1 − e^{−√ | 2026-09-29 | hinein |
| tunnelbild_heraus | Rückkehr als Tunneln (Schrödinger-Lesart D·G″ = Γ·G) als eigenes Bild | 07-feldgleichung.md: Rückkehr als Tunneln | 2026-09-29 | heraus |
| q_dynamisch | Q ist dynamisch, nicht topologisch | 03-von-der-galaxie-bis-zum-kosmos.md: Q ist dynamisch | 2026-09-29 | hinein |
| zaehlung_je_radiant | Q = 2π − 1 bei Zählung der Rückkehr je Radiant, als offene Lesart | 03-von-der-galaxie-bis-zum-kosmos.md: Zählung je Radiant | 2026-09-29 | hinein |
| elemente_faltung_zu_klein | Die eigene Faltung von Atom und Kern ist 29–36 Größenordnungen zu klein | 02-vom-teilchen-bis-zum-stern.md: 29–36 Größenordnungen | 2026-09-29 | hinein |
| zustandsschwelle_masse | Ein Trichter bestimmt den Zustand seines Stoffs ab etwa 10²¹ kg | 02-vom-teilchen-bis-zum-stern.md: etwa 10²¹ kg | 2026-09-29 | hinein |
| wellenform_kein_mitatmen | Eine Wellenform gäbe kein Mitatmen (w ∝ √M, Milchstraße 6·10⁻⁶) | 07-feldgleichung.md: Milchstraße 6·10⁻⁶ | 2026-09-29 | hinein |
| fluchtmass_dunkle_energie | Fluchtmaß mit verzögerter Ankunft: Ω_X = 0,32, w(0) = −0,99, im eigenen Atmen 0,36, mit dem Vakuumteil 0,56 von 0,687; übrige 18 %, w₀ = −1,01, w_a = −0,61, Maximum z ≈ 0,01, fσ₈ 9–15 % | 03-von-der-galaxie-bis-zum-kosmos.md: Geschachtelt gezählt und ab der Übergabe I_F | 2026-09-30 | heraus |
| kehrseite_bindung | Dunkle Energie zum Teil als Kehrseite der Bindung (Ankunft über die Rückkehr) | 03-von-der-galaxie-bis-zum-kosmos.md: Kehrseite der Bindung | 2026-09-30 | heraus |
| aussenfeld_frei_angepasst | Frei angepasst verlangen die Kurven ein gemeinsames Feld senkrecht 0,010 ± 0,009·a₀, gleichgerichtet 0,003 ± 0,003·a₀ | 03-von-der-galaxie-bis-zum-kosmos.md: Frei angepasst verlangen die Kurven selbst | 2026-09-30 | heraus |
| df2_df4_aussenfeld | DF2 und DF4 verlangten σ_v²/σ_v,N² ≤ 2,4 bzw. 1,5; erst ein 11- bzw. 206-fach stärkeres Außenfeld | 03-von-der-galaxie-bis-zum-kosmos.md: DF2 und DF4 verlangten | 2026-09-29 | hinein |
| scherung_schleifen | Scherung der Schleifen kleiner als 10⁻⁶ der Atemrate | 01-grundlagen.md: Scherung | 2026-09-29 | hinein |
| hauptreihe_faktor | Die Temperaturskala trifft die Hauptreihe auf den Faktor 0,7–1,3 | 02-vom-teilchen-bis-zum-stern.md: Faktor 0,7–1,3 | 2026-09-29 | hinein |
| alpha_s_planck | 1/α_s ≈ 53 an der Planck-Skala, ¦q̃¦ = 1,2·10⁻²⁰, m_p +41 % / −30 % | 02-vom-teilchen-bis-zum-stern.md: 1/α_s ≈ 53 | 2026-09-29 | hinein |
| umkehr_z_je_masse | Die örtliche Regel setzt zuerst bei z ≈ 20 (10⁸ M☉) ein | 03-von-der-galaxie-bis-zum-kosmos.md: z ≈ 20 (10⁸ | 2026-09-29 | hinein |
| kollaps_frueher | Nach der Umkehr kollabiert ein Bereich in 0,3–0,96 der Zeit einer Staubkugel | 03-von-der-galaxie-bis-zum-kosmos.md: 0,3–0,96 | 2026-09-29 | hinein |
| linse_300_kpc | Linse bei N* = 2,4 um 3–7 % (200 kpc) und 19–39 % (300 kpc) unter der Relation | 03-von-der-galaxie-bis-zum-kosmos.md: 19–39 % | 2026-09-29 | hinein |
| kaleidoskop_fenster | Zelle oberhalb des 1,5-fachen Durchmessers nicht mehr prüfbar (27–42 Gpc) | 03-von-der-galaxie-bis-zum-kosmos.md: 1,5-fachen Durchmessers | 2026-09-29 | hinein |
| warp_gas_90_grad | Einfallendes Gas dreht die Knotenlinie nur bis 90°, bei der gemessenen Rate nach 0,14 Mrd. Jahren | 03-von-der-galaxie-bis-zum-kosmos.md: Knotenlinie nur bis 90° | 2026-09-29 | hinein |
| ds_milchstrasse | Milchstraße: Δs = 197,4 AE | 03-von-der-galaxie-bis-zum-kosmos.md: 197,4 AE | 2026-09-30 | hinein |
| ds_r_milchstrasse | Δs/R Milchstraße 6,3·10⁻⁸ | 03-von-der-galaxie-bis-zum-kosmos.md: 6,3·10⁻⁸ | 2026-09-30 | hinein |
| rueckstellbilanz | ΔL_Raum/R = ΣΔs/R = 0,083 | 03-von-der-galaxie-bis-zum-kosmos.md: ΣΔs/R = 0,083 | 2026-09-30 | hinein |
| m_b_q | M_b,Q ≈ 4,45·10¹³ M☉ Baryonen (heute) | 03-von-der-galaxie-bis-zum-kosmos.md: 4,45·10¹³ M☉ | 2026-09-30 | hinein |
| umlauf_an_l | Sonne: Umlauf an ℓ 84,5 Mio. Jahre | 07-feldgleichung.md: Umlauf an ℓ 84,5 Mio | 2026-09-30 | hinein |
| pal14_newton | Newton-Wert Pal 14 0,665 km/s | 03-von-der-galaxie-bis-zum-kosmos.md: (0,665; 5,09; 7,15 km/s) | 2026-09-30 | hinein |
| ds_erdmasse | Δs = GM/(3c²) = x · 1,4783 mm | 01-grundlagen.md: x · 1,4783 mm | 2026-09-30 | hinein |
| uhr_neutronenstern | Neutronenstern: 69,5 Tage | 04-zeit-verschraenkung-syntropie-gedaechtnis.md: Neutronenstern 69,5 Tage | 2026-09-30 | hinein |
| gw170817_laenge | GW170817 590 m | 08-licht-und-trichter.md: GW170817 590 m | 2026-09-30 | hinein |
| v_erde | Erde: 29,78 km/s | 02-vom-teilchen-bis-zum-stern.md: Erde: 29,78 km/s | 2026-09-30 | hinein |
| g_roter_riese | Roter Riese: 0,0085 | 02-vom-teilchen-bis-zum-stern.md: ¦ 0,0085 ¦ | 2026-09-30 | hinein |
| einrasten_2_pi | erstes Einrasten bei (N − 1)k ≈ 0,57–0,63 W | 06-torus-dynamik.md: 0,57–0,63 W | 2026-09-30 | hinein |
| sonnenrand | Sonnenrand: 0,8756″ + 0,8756″ = 1,7512″ (R☉ = 6,957·10⁸ m) | 08-licht-und-trichter.md: 1,7512″ | 2026-09-30 | hinein |
| fluchtmoden_faser | Fluchtmoden gehen als Operator je Raumpunkt in die innere Richtung über und tragen wie der Vakuumteil mit w = −1 bei (w_a < 0); Zahl der Zyklen je Raumpunkt offen | 03-von-der-galaxie-bis-zum-kosmos.md: je Raumpunkt in die innere Richtung über | 2026-09-30 | hinein |
| vakuumteil_rueckstellung | Der Vakuumteil ist die Rückstellung des Ganzen; die Bindung macht Fluchtmoden inaktiv | 03-von-der-galaxie-bis-zum-kosmos.md: Der Vakuumteil ist die Rückstellung des Ganzen. Die Bindung | 2026-09-30 | hinein |
| doppelungen_aufgeloest | Jede Aussage steht nur an einer Stelle: Q-Herkunft, Zelle und Schwelle, Zahlen der Leiter, Erde und Sonne, Verstärkung an r*, Resonanzen | 00-bauplan.md: Was noch fehlt | 2026-09-30 | hinein |
| frisch_00_09 | Frische Werte statt Rundung: g der Erde 9,82, α_G = 5,9·10⁻³⁹, Eddington 72,6 M_* = 134,6 M☉, Ω_s(1090) = 2,6·10⁻¹⁰, Streuung 0,27–0,30, σ(Q) 0,065 (1,3σ), fσ₈ 1,2·10⁻³, Kreislauf 2,6·10¹¹ Jahre, 567-fach | 03-von-der-galaxie-bis-zum-kosmos.md: 0,27–0,30 dex | 2026-09-30 | hinein |
| frisch_sparc | SPARC: 1,06 ± 0,085 ± 0,11, Spanne bei Υ_b = 0,7 fest 0,96–1,19, Drittel 0,75 / 1,51 / 1,06 mit Anstieg +0,09 ± 0,08 (nicht gesichert), 74 Galaxien, Draco 1,2 km/s, Haufenstoß bis 1,8 Mrd. Jahre | 03-von-der-galaxie-bis-zum-kosmos.md: 0,75 / 1,51 / 1,06 | 2026-09-30 | hinein |
| frisch_zwerge | DF2 und DF4: σ_v²/σ_v,N² ≤ 2,35 bzw. 1,44, Außenfeld 10,2- bzw. 179-fach | 03-von-der-galaxie-bis-zum-kosmos.md: 10,2- bzw. 179-fach | 2026-09-30 | hinein |
| frisch_universum | Masse im eigenen Horizont 7,4·10²² M☉; Quadrupol des Sonnenorts 5,6·10⁻³⁰ s⁻² | 03-von-der-galaxie-bis-zum-kosmos.md: 7,4·10²² M☉ | 2026-09-30 | hinein |
| frisch_06 | Einrasten 0,57–0,63 W; Kupferdraht: 0,0010–0,6 Hz vor der Falte, Ringfrequenz 0,20–1,5 Hz (unterkritisch) | 06-torus-dynamik.md: 0,0010–0,6 Hz | 2026-09-30 | hinein |
| frisch_05 | Fluchtmoden höchstens 567-fach; Fluchtmoden gehen je Raumpunkt in die innere Richtung über | 05-zustandsraum.md: je Raumpunkt in die innere Richtung | 2026-09-30 | hinein |
| frisch_07 | Rand bei 0,66·r*·Λ/ϖ_R; Resonanz 5/4 bei 0,94 r*; Gleichstand 3,4 R_d; Ringfaktor 2,3–3,0 | 07-feldgleichung.md: 0,66·r*·Λ/ϖ_R | 2026-09-30 | hinein |
| frisch_08_09 | GP-B: 6604; Mitführung 40,8 mas/Jahr längs der Erdachse, zum Leitstern IM Peg 39,1; α/α_G = 1,23555·10³⁶ | 08-licht-und-trichter.md: vom Kreisel zum Leitstern IM Peg gesehen 39,1 | 2026-09-30 | hinein |
| frisch_01_04 | Doppelte Tabelle Erde und Sonne in 01 gestrichen (steht in 02); Fluchtmaß trägt bei | 01-grundlagen.md: trägt als drittes, positives Glied bei | 2026-09-30 | hinein |
| frisch_04 | Kreislauf: mit konstanter Dunkler Energie ab 2,6·10¹¹ Jahren verletzt; Rate endet beim Vakuumteil | 04-zeit-verschraenkung-syntropie-gedaechtnis.md: ab etwa 2,6·10¹¹ Jahren | 2026-09-30 | hinein |
| frisch_02 | g der Erde als GM/R² 9,82; Eddington 72,6 M_* = 134,6 M☉; α_G = 5,9·10⁻³⁹ | 02-vom-teilchen-bis-zum-stern.md: 9,82 | 2026-09-30 | hinein |
| frisch_09 | α/α_G = 1,23555·10³⁶; Umlauf bei r* ohne Zahl der Verstärkung (steht in 08) | 09-elektromagnetismus-und-frequenzleiter.md: 1,23555·10³⁶ | 2026-09-30 | hinein |

## Katalog

| Skript | Klasse | Ziel |
|---|---|---|
| a0_lesungen | interessant | a₀-Lesungen |
| atmen_oder_zelle | nutzlos |  |
| atmen_ohne_skala | übernommen | zustand_vakuum |
| atmen_w | nutzlos |  |
| ausgedehnte_massen | übernommen | zusatzpotential |
| aussenfeld_kosmisch | interessant | Äußeres Feld als Einbettung in die kosmische Rückkehr, g_e = (ln q_kosm)²·a₀ = 0,029·a₀ |
| aussenfeld_sparc | interessant | Wert der Schwelle |
| axiom_mass | interessant | Axiome und Maß |
| bistabilitaet | übernommen | preisach_ausgang |
| camb_q | nutzlos |  |
| carter_herleitung | übernommen | carter_faktor |
| compton_carter | nutzlos |  |
| compton_faltung | übernommen | faltungsrate |
| de_bilanz_glieder | interessant | Bilanz aller Glieder |
| de_fluchtmoden | interessant | DE als verdrängte Fluchtmoden |
| de_fluchtmoden_erneuerung | interessant | Fluchtmaß mit Erneuerung |
| de_fluchtmoden_hierarchie | interessant | Fluchtmaß hierarchisch |
| de_fluchtmoden_mass | übernommen | fluchtmass_bereich |
| de_fluchtmoden_verzoegert | interessant | Fluchtmaß verzögert über die Rückkehr |
| de_invariante | interessant | Dunkle Energie als Invariante der Grundwerte |
| de_kanal_atmen | interessant | Kanal ins Atmen |
| de_kehrseite_bindung | interessant | DE als Kehrseite der Bindung |
| de_laufzeit | interessant | Laufzeit des Übergangs |
| de_phasenseite | interessant | Dunkle Energie über die Phasenseite des Ganzen |
| de_ratenschwelle | interessant | Dunkle Energie über die Ratenschwelle |
| de_rest_anstieg | übernommen | lauf_eigenes_atmen |
| de_uebergabe_anteil | interessant | DE-Rest als Übergabe der Mode |
| de_zyklen | interessant | Dunkle Energie je Zyklus der Zelle |
| dicke_alter_scheiben | übernommen | 07_offen_31_1 |
| dicke_alter_scheiben_gegenprobe | nutzlos |  |
| ebenen_projektion | übernommen | glaettung_innen |
| ebenengrenze_zwerge | interessant | Ebenengrenze der Hierarchie als Umgebungsgröße der Zwerge |
| eine_struktur | übernommen | temperatur_innen |
| elemente_analyse | übernommen | 06_siede_schmelz_min |
| em_phasenseite | übernommen | moment_magnetisch |
| engstelle_global | übernommen | faltung_horizont |
| engstelle_rate | interessant | Engstelle und Schwarze Löcher |
| feldgleichung | nutzlos |  |
| folgerungen | nutzlos |  |
| formeln_andere_gebiete | übernommen | entropie_ruecklauf |
| forschung_uebergang | nutzlos |  |
| freiheitsgrad | nutzlos |  |
| fruehe_galaxien_a0 | interessant | Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell |
| fruehes_wachstum | nutzlos |  |
| ftm_gesamtlauf | nutzlos |  |
| ftm_gesetze | übernommen | B_tully_fisher_grenzfall |
| ftm_kern | nutzlos |  |
| ftm_kern_alt | übernommen | ds |
| ftm_pruefstand | übernommen | 08_leiterzahl |
| fusion_trichter | übernommen | grenzmasse_zaehlregel |
| galaxie_als_ganzes | übernommen | massenmittel_A_scheibe |
| gedaechtnis_topologie | interessant | Gedächtnis als Topologie |
| gedaechtnis_umgebung | übernommen | sigma_plummer_jeans |
| gegenlauf_dunkle_energie | nutzlos |  |
| gemeinsame_innere_richtung | interessant | Eine oder viele innere Richtungen (Bündel je Falte gegen gemeinsamen Knoten bei Λ) |
| geometrische_muster | übernommen | flussdimension |
| gesamtbild | nutzlos |  |
| gewicht_preisach | übernommen | gewicht_w_preisach |
| gewicht_rueckkehr | übernommen | schalen_dauer |
| gewicht_scheibenfeld | übernommen | a_scheibe_ebene |
| gewicht_scheibenfeld_gegenprobe | nutzlos |  |
| glaettung_stetig | übernommen | gezeitendichte |
| gleichung_ganzes | übernommen | erneuerung |
| gleichung_zeit | übernommen | rueckkehr |
| gravitationswellen | übernommen | hulse_taylor |
| halbgerade_hierarchie | interessant | Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell |
| halbgerade_spitze | übernommen | nu_aus_modul |
| harvey_nachlauf | interessant | Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell |
| haufen_nach_naht | interessant | Haufen nach der Naht |
| herkunft_q | nutzlos |  |
| hierarchie | nutzlos |  |
| holomorphie_rueckkehr | übernommen | glaettung_innen |
| hopf_knoten | interessant | Hopf-Knoten |
| horizont_kreislauf | übernommen | horizontentropie |
| identifikation_galaxie | übernommen | 05_resolvente_flucht |
| info_geometrie | übernommen | fisher_metrik |
| innere_richtung_zeit | übernommen | rueckkehr_geschlossen |
| invarianten_ebenen | übernommen | d_LSS_beobachter |
| kac_schliessung | übernommen | gewicht_w |
| kern_chi | nutzlos |  |
| kern_gedaechtnis | interessant | Gedächtniskern |
| kern_nome | nutzlos |  |
| kern_schleife | nutzlos |  |
| kollaps_nach_umkehr | interessant | Kollaps nach der Umkehr |
| kontakt_atmen | übernommen | rc_kette_nu |
| kontraktion_atemrate | übernommen | hintergrund_modell |
| kopfgleichung_kern | übernommen | kern_zeit |
| kosmologie | nutzlos |  |
| kovariante_e0 | übernommen | gezeit_betrag |
| kupfer_selbstpulsieren | übernommen | kupfer_hopf_abstand |
| labor_vorhersage | übernommen | kipp_schleife |
| lambda_bruecke | interessant | Λ ohne ν festlegen (lambda_bruecke) |
| lambda_und_rueckstellung | nutzlos |  |
| licht_trichter | übernommen | delta |
| linse_ohne_freiheit | übernommen | haufenstoss_trennzeit |
| linsen | übernommen | ablenkung_linse |
| linsen_isoliert | interessant | Linsen um isolierte Galaxien (KiDS, Brouwer 2021; Mistele 2024) |
| massenverhaeltnis | übernommen | 02_im_tau_G |
| maximum_schwelle | interessant | Vergessene Funde der Notizen (Durchsicht 01–75), noch nicht im Modell |
| mikro_nome | interessant | Mikroskala und Nome |
| mikroskala | übernommen | nome_qcd |
| mischregel_schale | interessant | Aufteilung an der w-Schale aus einem Rücklauf mit zwei Verlusten |
| mock_schatten | übernommen | Q_tilde_gammareihe |
| modell_h0 | interessant | H₀ im eigenen Atmen und Schallhorizont |
| modulraum_zahlen | interessant | Modulraum-Vorschläge |
| naht_rq_umkehr | interessant | Übergang an der w-Schale |
| netz_und_bilanz | nutzlos |  |
| ordnung_pruefstand | nutzlos |  |
| ort_oder_system | interessant | Ort oder System |
| orthogonalitaet_wurzeln | übernommen | epizyklen_summe |
| parameter_trio | nutzlos |  |
| q_a0_kopplung | übernommen | zellgrenze |
| q_form | übernommen | torus_streifen_nu |
| radiale_beschraenkung | übernommen | dynamikmasse_ort |
| ramanujan | übernommen | eisenstein_G |
| ramanujan_unbekannt | nutzlos |  |
| ramanujan_weitere | übernommen | drehwinkel_zelle |
| regelkreise | interessant | Regelkreise je Ebene |
| rest_mode_isolierte_galaxie | nutzlos |  |
| ringe_oertlich | übernommen | praezession_ring |
| rotation | übernommen | horizont_kerr |
| rueckkehr_2pi | nutzlos |  |
| rueckkehr_horizont | nutzlos |  |
| rueckkehr_zustandsraum | übernommen | kac_budget |
| rueckstellung_zeitform | nutzlos |  |
| rueckweg_faltung | nutzlos |  |
| schale_zeitlich | übernommen | uebergabe_vor_umkehr |
| schalter_chi | übernommen | chi |
| scheibe_feldform | übernommen | eta_scheibe |
| schliessungsregel | nutzlos |  |
| schluss_oben | übernommen | q_fenster |
| schwelle_anpassung | nutzlos |  |
| schwelle_je_galaxie | interessant | Wert der Schwelle |
| schwelle_modell_mess | nutzlos |  |
| schwelle_rueckkehrzeit | übernommen | kac_zeit_schluss |
| setzungen_gleichung | übernommen | verlust_gamma |
| skalen_gross | übernommen | 03a_ds_milchstrasse |
| skalen_klein | übernommen | ds |
| sparc_aussenfeld_rand | übernommen | 03_aussen_zahl |
| sparc_boseform | übernommen | a0_an_punkten |
| sparc_flaechendichte | übernommen | toomre_sigma_dyn |
| sparc_spektraldimension | übernommen | ds_anpassung |
| sparc_streuung | übernommen | streuung_messfehler |
| sparc_test | nutzlos |  |
| sparc_universal | übernommen | rest_galaxie |
| stern_trichter | nutzlos |  |
| stoerungen_rekombination | übernommen | 03a_rest_w_rekombination |
| straffung_rueckstellung | nutzlos |  |
| strittig_zweiter_weg | interessant | Strittige Punkte, zweiter Weg |
| stufe_zwei | übernommen | bose_leiter_integral |
| syntropie | übernommen | grenzmasse_hawking |
| tau_ebenen | übernommen | spiegel_radius |
| tau_realteil | nutzlos |  |
| temperatur_tiefe_ebenen | übernommen | masse_aus_temperatur |
| tensor | übernommen | Q_tilde |
| theta_bezug | übernommen | gewicht_hintergrund |
| topologie_ganzes | übernommen | dirac_verhaeltnis |
| torus_dynamik | übernommen | rastbreite |
| torus_wechselkraefte | übernommen | 06_einrasten_min |
| trenntests_b_dm | interessant | Trenntests Bindung gegen Dunkle Materie |
| triade_ebenen | übernommen | realteil_quadrat_punkt_triade |
| trichter_invarianten | übernommen | schnittmasse_glaettung |
| trichterkette_gewicht | übernommen | dichte_umkehr_heute |
| trichterkette_sonne | interessant | Trichterkette der Sonne (Sgr A*, Milchstraße, Sgr-Zwerg, LMC, M31, Virgo, Großer Attraktor, Shapley) |
| tunnel_aussenfeld | übernommen | nu_aussenfeld |
| uebergabe_faltung | übernommen | uebergabe_faltung |
| umlaufsinn | übernommen | 07_coth_schwelle |
| vakuum_richtungen | interessant | Vakuumteil je Richtung |
| wachstum_w_qs | interessant | Gewicht vor Q_s im linearen Wachstum |
| wahl_im_tau | übernommen | 07_im_tau_erdbahn |
| warp_amplitude | übernommen | warp_viertelzeit |
| warp_drehmoment | übernommen | warp_aufwicklung |
| warp_praezession | nutzlos |  |
| wechselwirkungen | übernommen | schleife_thermisch |
| wep_glaettung | interessant | Gleichheit des Fallens und Glättungsschwelle mit der örtlichen Regel |
| wicklung_rueckkehr | übernommen | kac_rueckkehrzeit |
| wiederkehrende_muster | übernommen | bose_leiter |
| wirkung_luecken | interessant | Lücken der Wirkung |
| wirkung_oertliche_regel | nutzlos |  |
| zeit_rueckkehr | übernommen | tau_D |
| zeitgesetz_a0 | übernommen | omega_vakuum |
| zelle_form | übernommen | formfaktor_zelle |
| zelle_frueher_kosmos | nutzlos |  |
| zustandsraum | übernommen | 05_delta_koeffizient_3 |
| zwerge_aussenfeld | übernommen | sigma_plummer_jeans |
| zwerge_tiefengedaechtnis | übernommen | 02_helium_strahlend |
| zwerge_uebergabe | interessant | Zwerge vor der Übergabe |
| aussenfeld_sparc_a012.txt | nutzlos |  |
| glaettung_stetig_a011.txt | nutzlos |  |
| gravitationswellen_a011.txt | nutzlos |  |
| radiale_beschraenkung_a011.txt | nutzlos |  |
| rest_mode_isolierte_galaxie_a012.txt | nutzlos |  |
| schliessungsregel_a011.txt | nutzlos |  |
| warp_drehmoment_a011.txt | nutzlos |  |
| warp_praezession_a011.txt | nutzlos |  |
| zwerge_aussenfeld_a011.txt | nutzlos |  |
| ordnung_err.txt | nutzlos |  |
| ordnung_pruefstand_E.npy | nutzlos |  |
| __pycache__/faltungsknoten.cpython-313.pyc | nutzlos |  |

## Änderungen seit dem letzten Lauf

Voriger Lauf: Commit 9068973, 2026-09-30T19:48:48.

### Werte

| Teil | Kennung | Art | Änderung |
|---|---|---|---|
| gesetze | Z_spektralluecke_kreis | geändert | vergleich: 39.4776 → 39.4776 |

### Modelltext (git diff modell/)

keine Änderung
