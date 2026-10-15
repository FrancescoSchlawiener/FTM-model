# Das FTM-Weltmodell und seine Simulation

## 1 Zusammenfassung

Die Simulation rechnet das Modell (modell/00–09) nach: ein Kern (kern/) hält alle Formeln, Konstanten, Messwerte, Eingaben und Setzungen je einmal; Register in pruefung/ und forschung/ binden jede Zahl, jedes Gesetz, jede Gegenprobe, jedes Thema und jeden Beschluss an den Modelltext. Ein Lauf (lauf.py) rechnet alles frisch und schreibt bericht.md, register.json und die erzeugten Abschnitte dieser Datei. Was das Modell sagt, steht nur in modell/.

## 2 Aufbau der Simulation

- **Kern.** Formeln (Dekorator, ein Eintrag je Größe), Symbolregister, Parametersätze, Messwerte mit Quelle, Eingaben mit Fundstelle, Setzungen.
- **Prüfzahlen.** Jede Zahl des Modelltexts einmal, aus dem Kern frisch gerechnet: OK, ABWEICHUNG oder OFFEN.
- **Gesetze und Gegenproben.** Gesetze je Klasse; jede Formel mit Gegenprobe auf einem zweiten Weg, soweit es einen gibt.
- **Wächter.** Streng (ein Fund stoppt den Lauf) und zählend (Lage im Bericht); pytest läuft im Lauf mit.
- **Forschung.** Themen mit Invariante, Versuchen und gesperrten Wegen; Beschlüsse für jede Änderung am Modelltext; offene Beschlüsse zur Entscheidung.
- **Archiv.** Unbenutztes und Entferntes liegt in weltmodell/archiv/, nichts wird gelöscht.

## 3 Blickrichtungen im Register

Jedes Symbol und jedes Gesetz nennt seine Ebene: L (am Ort), I (innere Richtung), Z (Zustandsraum), G (das Ganze), und den Zustand (fest, Epoche, je System); Formeln nennen zusätzlich die Rolle in der Triade. Die Definitionen der Ebenen und Rollen stehen im Modell (Bausteine 01, 05); die Tafeln unten sind daraus erzeugt.

## 4 Symboltafel

<!-- erzeugt: symbole -->
513 Zeichen, jedes einmal.

| Zeichen | Bedeutung | Einheit | Ebene | Zustand |
|---|---|---|---|---|
| `G` | Gravitationskonstante | m³/(kg·s²) | LIZG | fest |
| `c` | Lichtgeschwindigkeit, Skalierungsfaktor der Hülle | m/s | LIZG | fest |
| `ħ` | reduzierte Planck-Konstante | J·s | LIZG | fest |
| `k_B` | Boltzmann-Konstante | J/K | LIZG | fest |
| `e` | Elementarladung | C | LIZG | fest |
| `N_A` | Avogadro-Konstante | mol⁻¹ | LIZG | fest |
| `ε₀` | elektrische Feldkonstante | F/m | LIZG | fest |
| `m_u` | atomare Masseneinheit m(¹²C)/12 | kg | LIZG | fest |
| `m_p` | Protonmasse | kg | LIZG | fest |
| `m_e` | Elektronmasse | kg | LIZG | fest |
| `α` | Feinstrukturkonstante r_e/λ̄_e | 1 | LIZG | fest |
| `r_e` | klassischer Elektronenradius | m | LIZG | fest |
| `a_B` | Bohr-Radius | m | LIZG | fest |
| `R_p` | Ladungsradius des Protons | m | LIZG | fest |
| `M☉` | Sonnenmasse | kg | LIZG | fest |
| `GM☉` | Gravitationsparameter der Sonne | m³/s² | LIZG | fest |
| `R☉` | Sonnenradius | m | LIZG | fest |
| `P_L☉` | Leuchtkraft der Sonne | W | LIZG | fest |
| `T_eff☉` | Effektivtemperatur der Sonne | K | LIZG | fest |
| `M⊕` | Erdmasse | kg | LIZG | fest |
| `GM⊕` | Gravitationsparameter der Erde | m³/s² | LIZG | fest |
| `R⊕` | Erdradius | m | LIZG | fest |
| `g⊕` | Neigung an der Erdoberfläche, 9,81 m/s² | m/s² | L | fest |
| `AE` | Astronomische Einheit | m | LIZG | fest |
| `pc` | Parsec (kpc, Mpc, Gpc Vielfache) | m | LIZG | fest |
| `Lj` | Lichtjahr | m | LIZG | fest |
| `Jahr` | julianisches Jahr (Mrd. Jahre Vielfache) | s | LIZG | fest |
| `ℓ_P` | Planck-Länge | m | LIZG | fest |
| `m_Pl` | Planck-Masse | kg | LIZG | fest |
| `T_CMB` | Temperatur der Hintergrundstrahlung heute | K | G | Epoche |
| `M` | Masse (Faltung eines Körpers oder Systems) | kg | L | je System |
| `m` | Masse eines Teilchens | kg | L | je System |
| `R` | Radius der Hülle eines Körpers | m | L | je System |
| `r` | Abstand vom Zentrum im Raum | m | L | je System |
| `x` | Masse in Erdmassen, x = M/M⊕ (A1) | 1 | L | je System |
| `y` | Radius in Erdradien, y = R/R⊕ (A1) | 1 | L | je System |
| `g` | Neigung an der Kante des Trichters, GM/R² | m/s² | L | je System |
| `Δs` | Faltungstiefe GM/(3c²) | m | L | je System |
| `dτ/dt` | Abwicklungstakt √(1 − 6Δs/r) | 1 | L | je System |
| `t_E` | Eigenzeit am Ort (im Modelltext τ in dτ/dt) | s | L | je System |
| `δ` | Faltungsdifferential 1 − dτ/dt | 1 | L | je System |
| `r₀` | Grenzradius 6Δs (Takt steht) | m | L | je System |
| `r*` | Schwellenradius √(GM/a₀) = √(3Δs·ℓ₀) | m | L | je System |
| `λ̄` | reduzierte Compton-Länge ħ/(mc) | m | L | je System |
| `α_G` | gravitative Kopplung GM²/(ħc) = 3Δs/λ̄ | 1 | L | je System |
| `M_*` | Massenskala der Sterne m_Pl³/m_p² | kg | L | fest |
| `P_Edd` | Eddington-Leuchtkraft (9/2)·m_p·c³·Δs/r_e² | W | L | je System |
| `T` | Temperatur | K | L | je System |
| `T_H` | Temperatur einer Falte ħc/(24πkΔs) | K | L | je System |
| `M_T` | Masse der Falte mit T_H = T (die Epoche kommt über T, etwa T_CMB) | kg | L | je System |
| `t_F` | Lebensdauer der Falte 3840·(S/k)·Δs/c | s | L | je System |
| `ω_C` | Compton-Takt Mc²/ħ | s⁻¹ | L | je System |
| `a_J` | Verwindungslänge J/(Mc) | m | L | je System |
| `j` | Spinzahl a_J/(3Δs) | 1 | L | je System |
| `Q_el` | Ladung n·e (Windungszahl) | C | L | fest |
| `α_abl` | Ablenkwinkel des Lichts 12Δs/b | rad | L | je System |
| `α_∞` | kollektive Ablenkung 6πΔs/r* = 2π(v_flach/c)² | rad | L | je System |
| `b` | Stoßparameter des Lichts | m | L | je System |
| `Δφ_P` | Periheldrehung je Umlauf | rad | L | je System |
| `a_b` | große Halbachse einer Bahn | m | L | je System |
| `e_B` | Exzentrizität einer Bahn | 1 | L | je System |
| `α_S` | Winkeldurchmesser des Schattens 2b/D_A, b = 3√27·Δs | rad | L | je System |
| `D_A` | Abstand des Beobachters | m | L | je System |
| `r̃_ISCO` | innerste stabile Kreisbahn in Einheiten 3Δs | 1 | L | je System |
| `f*_gw` | Schwellenfrequenz der Gravitationswellen Ω(r*)/π | Hz | L | je System |
| `P_b` | Umlaufzeit eines Doppelsterns | s | L | je System |
| `Ṗ_b` | Abnahme der Umlaufzeit | 1 | L | je System |
| `m₁` | Masse der ersten Komponente (im Kern in M☉) | kg | L | je System |
| `m₂` | Masse der zweiten Komponente (im Kern in M☉) | kg | L | je System |
| `η_m` | Massenverhältnis m₁m₂/M² | 1 | L | je System |
| `η` | Massenverhältnis M_dyn/M_lens | 1 | L | je System |
| `ρ_J` | Gezeitendichte 3𝒯/(4πG) | kg/m³ | L | je System |
| `ℓ` | Glättungslänge, ℓ³ = GM/𝒯 | m | L | je System |
| `𝒯` | Gezeit des geglätteten Elternfelds √(3/2)·\|Ē^TF\| | s⁻² | L | je System |
| `σ` | Faltungspotential (kovariante Fassung) | 1 | L | je System |
| `σ_v` | Geschwindigkeitsstreuung | m/s | L | je System |
| `σ_el` | elektrische Leitfähigkeit | S/m | L | je System |
| `f_g` | Anteil der Schwere an der Steifigkeit | 1 | L | je System |
| `Λ_QCD` | Skala der starken Kopplung | J | L | fest |
| `a` | Beschleunigung nach der örtlichen Regel ν(a_N)·a_N | m/s² | L | je System |
| `a_N` | Newton-Beschleunigung \|g_N\| am Ort | m/s² | L | je System |
| `g_N` | Newton-Feld (Vektor) als Funktion des Orts | m/s² | L | je System |
| `a_i` | inneres Feld eines Systems | m/s² | L | je System |
| `a_e` | äußeres Feld am Ort eines Systems | m/s² | L | je System |
| `v` | Kreisgeschwindigkeit | m/s | L | je System |
| `v_flach` | flache Rotationsgeschwindigkeit (GMa₀)^{1/4} | m/s | L | je System |
| `Σ_a₀` | Grenzflächendichte a₀/(πG) | kg/m² | L | seit z≈10 |
| `M(<r)` | eingeschlossene Masse (Hernquist-Profil) | kg | L | je System |
| `b_H` | Skalenlänge des Hernquist-Profils | m | L | je System |
| `R_zyl` | Zylinderradius in der Scheibe | m | L | je System |
| `z_h` | Höhe über der Mittelebene | m | L | je System |
| `v_R` | Radialgeschwindigkeit | m/s | L | je System |
| `v_z` | Vertikalgeschwindigkeit | m/s | L | je System |
| `J_z` | Drehimpuls je Masse um die Achse | m²/s | L | je System |
| `t` | Zeit längs der Bahn | s | L | je System |
| `ω_orb` | Umlauffrequenz √((4πG/3)ρ̄) | s⁻¹ | L | je System |
| `λ_m` | Übergaberate der Mode ξ, ω_orb/N* | s⁻¹ | L | je System |
| `N*` | Radiant der Bahn bis zur Übergabe | rad | L | je System |
| `λ_K` | Zerfallsrate des Gedächtniskerns der Kopfgleichung | s⁻¹ | L | je System |
| `Δ_B` | Bindungsgrad eines Bereichs gegen ρ_krit, mit Rückläufen gezählt | 1 | L | je System |
| `θ` | Ausdehnungsrate eines Bereichs im Raum | s⁻¹ | L | je System |
| `θ̃` | kleinste Ausdehnungsrate der Geschichte, max(min θ, 0) | s⁻¹ | L | je System |
| `a₀` | Schwelle c²/(2Λ), Standardwert 1,1·10⁻¹⁰ | m/s² | I | seit z≈10 |
| `Λ` | Schließlänge der Rückkehr (innere Länge, keine Raumlänge) | m | I | seit z≈10 |
| `ℓ₀` | Länge der Wirkung c²/a₀ = 2Λ | m | I | seit z≈10 |
| `τ_D` | Rückkehrzeit ℓ₀/c = c/a₀ | s | I | seit z≈10 |
| `Γ` | Verlustrate a_N/c (Taktgefälle) | s⁻¹ | LI | je System |
| `θ_c` | kritische Ausdehnungsrate a₀/c = 1/τ_D | s⁻¹ | I | seit z≈10 |
| `ϖ` | Glättungstiefe der Rückkehr (Halbgerade) | m | I | je System |
| `ϖ_R` | Rückkehrtiefe, Punktmasse ϖ_R = r* | m | I | je System |
| `D` | Rückkehrdiffusion ϖ_R²/τ_D (Punktmasse GM/c) | m²/s | I | je System |
| `u` | Feld der Rückkehr in (x, ϖ) | m/s² | I | je System |
| `p` | Laplace-Variable des Kerns | s⁻¹ | I | je System |
| `s_p` | Exponent der einen Rückkehr √((p + Γ + θ)·τ_D) | 1 | I | je System |
| `û` | projizierte Rückkehr e^{−s_p}: q, 1 − w, K̂ | 1 | I | je System |
| `K̂` | Laplace-Bild des Kerns e^{−√((p+Γ)τ_D)} | 1 | I | je System |
| `x_θ` | Rate mal Zeitkonstante eines Gedächtnisses, θ₊τ_D | 1 | I | je System |
| `s` | Bose-Variable √(a_N/a₀), Punktmasse r*/r | 1 | LI | je System |
| `q` | Nome e^{2πiτ}; \|q\| = e^{−s} ist der Rückkehranteil, als Regelkreis die Schleifenverstärkung | 1 | LIZG | je System |
| `ν` | Verstärkung 1/(1 − q) = 1 + n_B(s) | 1 | LI | je System |
| `χ` | Schalter s/(e^s − 1) | 1 | LI | je System |
| `n_B` | Bose-Besetzung 1/(e^s − 1) | 1 | LI | je System |
| `Q̃` | Strukturwirkung −4∫_s^∞ x³n_B dx | 1 | LG | je System |
| `Q_s` | Strukturteil Q̃ + 4π⁴/15 | 1 | L | je System |
| `d` | Dimension der Stufe der Bose-Leiter | 1 | Z | fest |
| `s_d` | Maximum von s^d·n_B (Wien-Stelle) | 1 | Z | fest |
| `w` | Gewicht der Schließungsregel 1 − e^{−√(θ̃/θ_c)} | 1 | G | je System |
| `ξ` | Mode der Schließungsregel, (∂ξ)² = −1 | 1 | G | je System |
| `𝒫` | Schließungsoperator (1/N)ΣT_{mL} | 1 | Z | fest |
| `Π⊥` | Projektor auf Flucht- und Instabilitätsrichtungen | 1 | Z | je System |
| `L_eff` | Erzeuger der Bewegung mit Schließung, 𝒫L | s⁻¹ | Z | je System |
| `τ` | Modul des Torus; q = e^{2πiτ} | 1 | Z | je System |
| `τ_M` | Modul des Bahntorus, Re = ω_rad/ω_azim, Im = Bindung | 1 | Z | je System |
| `τ_Z` | Modul des flachen Torus der Zelle ℂ/(ℤ + τ_Zℤ) | 1 | Z | fest |
| `Re τ_M` | Umlaufphase ω_rad/ω_azim | 1 | Z | je System |
| `q_Δs` | Nome am Kern e^{2πiτ_M}, \|q_Δs\| = Δs/r | 1 | Z | je System |
| `p/k` | rationales Verhältnis einer Resonanz | 1 | Z | fest |
| `r_p/k` | Radius der Resonanz Re τ_M = p/k (Punktmasse) | m | L | je System |
| `Im τ_Z` | Imaginärteil des Moduls τ_Z | 1 | Z | fest |
| `R_T/r_T` | Radienverhältnis eines eingebetteten Torus | 1 | Z | fest |
| `b_S` | Breite des offenen Streifens im Rohrwinkel | rad | Z | fest |
| `ν_T` | Kehrwert des offenen Anteils der Torusfläche | 1 | Z | fest |
| `μ(X)` | Maß des ganzen Zustandsraums | 1 | Z | fest |
| `μ(A)` | Maß eines Bereichs A | 1 | Z | je System |
| `⟨T⟩` | mittlere Rückkehrzeit in A (Kac), in Schritten | 1 | Z | je System |
| `Σ1/⟨T_i⟩` | Kac-Budget einer Zerlegung, 1 bei Schließung | 1 | Z | fest |
| `ε_L` | Leiterzahl ℓ_P/ℓ₀ | 1 | Z | seit z≈10 |
| `M_min` | untere Schnittmasse ħ/(c²τ_D) | kg | Z | seit z≈10 |
| `M_c` | Schnittmasse r* = λ̄, ε_L^{1/3}·m_Pl | kg | Z | seit z≈10 |
| `M_U` | obere Schnittmasse c³τ_D/G | kg | Z | seit z≈10 |
| `ρ_el` | spezifischer Widerstand | Ω·m | L | je System |
| `α_R` | Temperaturbeiwert des Widerstands | K⁻¹ | L | je System |
| `E_Y` | Elastizitätsmodul | Pa | L | je System |
| `α_E` | Temperaturbeiwert des Elastizitätsmoduls | K⁻¹ | L | je System |
| `α_th` | Wärmeausdehnung | K⁻¹ | L | je System |
| `ρ_Cu` | Dichte des Drahts | kg/m³ | L | je System |
| `c_p` | spezifische Wärme | J/(kg·K) | L | je System |
| `k_th` | Wärmeleitfähigkeit | W/(m·K) | L | je System |
| `a_D` | Drahtradius | m | L | je System |
| `ℓ_D` | Drahtlänge | m | L | je System |
| `B` | Magnetfeld | T | L | je System |
| `Q_m` | mechanische Güte | 1 | L | je System |
| `σ_0` | mechanische Vorspannung | Pa | L | je System |
| `R_el,0` | elektrischer Widerstand bei ΔT = 0 | Ω | L | je System |
| `G_th` | Wärmeleitwert | W/K | L | je System |
| `C_th` | Wärmekapazität | J/K | L | je System |
| `m_eff` | effektive Masse der Mode | kg | L | je System |
| `Γ_B` | Kopplung Bℓ | T·m | L | je System |
| `k_st` | Steifigkeit der Mode k(ΔT) | N/m | L | je System |
| `ω₀` | Eigenfrequenz der Mode | s⁻¹ | L | je System |
| `γ_m` | mechanische Dämpfung m·ω₀/Q_m | kg/s | L | je System |
| `L` | Zelle, räumliche Schleife des Torus (mitbewegt; fest wie a₀ und ρ_m) | m | G | seit z≈10 |
| `L_min` | kürzeste Schleife ohne Kreispaare, 0,985·d_LSS | m | G | Epoche |
| `Q` | kosmisches Gewicht q/(1 − q) = Ω_c/Ω_b | 1 | G | fest |
| `q_kosm` | Rückkehranteil des Ganzen Q/(1 + Q) | 1 | G | fest |
| `ν_kosm` | mittlere Zahl der Rückläufe im Kosmos 1 + Q (Kac) | 1 | G | fest |
| `σ_Q` | Messfehler von Q aus σ(ω_b), σ(ω_c) | 1 | G | fest |
| `Q_min` | untere Grenze des Fensters von Q | 1 | G | fest |
| `Q_max` | obere Grenze des Fensters von Q (hängt an L_min, also am Alter des Beobachters) | 1 | G | Epoche |
| `H` | Atemrate ȧ/a | s⁻¹ | G | Epoche |
| `H₀` | Atemrate heute (im Kern in km/s/Mpc) | s⁻¹ | G | Epoche |
| `h_H` | H₀/(100 km/s/Mpc) | 1 | G | Epoche |
| `z` | Rotverschiebung | 1 | G | Epoche |
| `E(z)` | H(z)/H₀ im Hintergrund | 1 | G | Epoche |
| `d_C` | mitbewegter Abstand | m | G | Epoche |
| `d_LSS` | Durchmesser der letzten Streufläche (mitbewegt) | m | G | Epoche |
| `Θ` | 2πa₀/(cH) | 1 | G | Epoche |
| `Λ_s` | Vakuumteil π⁴/(15Λ²) | m⁻² | G | seit z≈10 |
| `Ω_s` | Vakuumteil als Anteil (π²/45)Θ² | 1 | G | Epoche |
| `K_s` | Vorfaktor π²/45 des Vakuumteils bei Θ = 1 | 1 | G | fest |
| `H_s` | Atemrate des Vakuumteils allein, H_s² = Λ_s c²/3 | s⁻¹ | G | seit z≈10 |
| `ρ_m` | mittlere Dichte der Materie, mitbewegt, ν_kosm·ρ_b | kg/m³ | G | fest |
| `ρ_b` | mittlere Dichte der Baryonen, mitbewegt | kg/m³ | G | fest |
| `Ω_b` | Anteil der Baryonen heute | 1 | G | Epoche |
| `Ω_m` | Anteil der Materie heute, Ω_b(1 + Q) | 1 | G | Epoche |
| `Ω_r` | Anteil der Strahlung heute | 1 | G | Epoche |
| `Ω_K` | Krümmungsanteil heute | 1 | G | Epoche |
| `Ω_DE` | Anteil der Dunklen Energie heute | 1 | G | Epoche |
| `H_DE` | Atemrate der Dunklen Energie | s⁻¹ | G | Epoche |
| `ω_b` | Ω_b h² | 1 | G | fest |
| `ω_c` | Ω_c h² | 1 | G | fest |
| `a₀(z)/a₀` | Schwelle früh gegen heute | 1 | G | Epoche |
| `f_N` | fester Anteil am Gesamtbudget im Strahlungszeitalter | 1 | G | Epoche |
| `N_eff` | effektive Zahl der Neutrinoarten | 1 | G | fest |
| `ΔN_eff` | Zusatz zu N_eff | 1 | G | Epoche |
| `S_H` | Horizontentropie πc²/(ℓ_P²H²) in Einheiten von k_B | 1 | G | Epoche |
| `Ṡ_i` | Entropieproduktion je Horizont in k_B/s | s⁻¹ | G | Epoche |
| `Λ̇/Λ` | Wachstumsrate der Schließlänge | s⁻¹ | G | Epoche |
| `w_Z` | Zustandsgleichung p/(ρc²) (w_Z0, w_Za, w_Zs) | 1 | G | Epoche |
| `Δw_Z` | Abstand der Zustandsgleichung von −1: 1 + w_Z | 1 | G | Epoche |
| `σ_B` | Faltungsbilanz ΣΔs − ΔL = 0 | m | G | fest |
| `R_Q` | Radius, an dem das örtliche q das kosmische erreicht | m | LG | je System |
| `M_b,Q` | Baryonenmasse, bei der R_Q die Grenze R₂₀₀ erreicht | kg | LG | Epoche |
| `Δ_ta` | Dichte gegen ρ_m einer heute umkehrenden Schale (1 + Δ_ta) | 1 | LG | Epoche |
| `R_ta` | Umkehrradius heute | m | LG | Epoche |
| `M_ges` | Maß des Ganzen (1 + Q)·M_b | kg | LG | je System |
| `a_GH` | Temperatur des Hubble-Horizonts als Beschleunigung cH/2π (heute nahe a₀, Θ = a₀/a_GH) | m/s² | G | Epoche |
| `d_H` | Teilchenhorizont eines Beobachters, mitbewegt | m | G | Epoche |
| `ΔL/R` | Dehnung der Raumschleifen je Radius = ΣΔs_b/R im eigenen Teilchenhorizont | 1 | G | Epoche |
| `Δs_H` | Faltung des eigenen Horizonts r_H/6 = c/(6H) (globale Engstelle) | m | G | Epoche |
| `S_Λ` | Horizontentropie des Vakuumteils allein (45/π³)(Λ/ℓ_P)² in Einheiten von k_B | 1 | IG | seit z≈10 |
| `d ln a₀/d ln a` | Gang der Schwelle mit dem Skalenfaktor (heute null) | 1 | IG | Epoche |
| `L·H/c` | Zelle in Hubble-Längen der Epoche | 1 | G | Epoche |
| `f³` | Formfaktor der Zelle ohne Kreispaare (Würfel 1, sechseckig 1,15, raumzentriert 1,30, flächenzentriert 1,41) | 1 | Z | fest |
| `(1+Q)a₀` | Grenze von (1 + Q)·a₀ aus der Zelle ohne Kreispaare, f³·3c⁴/(2Gρ_bL_min³) | m/s² | G | Epoche |
| `H_b` | Atemrate aus Bremsglied der Materie und Vakuumteil π⁴/(15Λ²), ohne Fluchtmaß | s⁻¹ | G | Epoche |
| `Ẋ_an` | Rate des über die Rückkehr angekommenen Fluchtmaßes X_an | s⁻¹ | ZG | Epoche |
| `H_eff` | Atemrate mit angekommenem Fluchtmaß H_b + ⅓Ẋ_an | s⁻¹ | G | Epoche |
| `ln a` | Atemmaß des Ganzen ln(L(t)/L(t₀)) = −ln(1 + z): die Epoche, auch in der Zukunft (z → −1) | 1 | G | Epoche |
| `t_A` | Weltzeit längs der Atemrate bis zur Epoche, ∫d ln a/H | s | G | Epoche |
| `Ṡ_H` | Wachstum der Horizontentropie −2S_H·Ḣ/H in k_B/s | s⁻¹ | G | Epoche |
| `Ω_X` | fester Rest der Atemrate im Vergleichshintergrund (Vakuumteil Ω_s, oder im Standardmodell gelesen 1 − Ω_m − Ω_r) | 1 | G | Epoche |
| `X` | inaktives Fluchtmaß je Baryon: ln-Maß der ganzen Tori im größten gebundenen Bereich, Regel gegen Newton (Zustandsraum, kein Stoff) | 1 | ZG | Epoche |
| `X_an` | über die Rückkehr angekommener Teil des Fluchtmaßes (Atmen der Raumschleifen δ ln L) | 1 | ZG | Epoche |
| `t_F/t_ta` | Übergabe bei der ersten Faltung in Einheiten der Umkehrzeit | 1 | LG | fest |
| `σ_M` | lineare Streuung der Dichte heute auf der Masse M (Kugelfilter) | 1 | G | Epoche |
| `D_+` | lineares Wachstum der Störungen D(a) | 1 | G | Epoche |
| `δ_ta` | lineare Schwelle der Umkehr: auf die Umkehrepoche fortgesetzter linearer Anfangskontrast | 1 | LG | Epoche |
| `A_s` | Amplitude des Anfangsspektrums | 1 | G | fest |
| `n_s` | Neigung des Anfangsspektrums | 1 | G | fest |
| `τ_re` | optische Tiefe der Reionisation | 1 | G | fest |
| `Σm_ν` | Summe der Neutrinomassen im Spektrum | eV | G | fest |
| `Δ(fσ₈)/fσ₈` | relative Änderung des linearen Wachstums fσ₈ durch die Schließungsregel gegen dasselbe Wachstum ohne Übergabe | 1 | G | Epoche |
| `Φ/c²` | Potentialtiefe einer linearen Störung in Einheiten von c² | 1 | G | Epoche |
| `Gitter` | Gitterform der Zelle (Würfel, sechseckig, raumzentriert, flächenzentriert) | 1 | Z | fest |
| `φ_D` | Drehwinkel einer Drehsymmetrie der Zelle (flacher Torus) | ° | Z | fest |
| `ℓ_W` | Ordnung einer Winkelmode am Himmel (2 Quadrupol, 3 Oktupol) | 1 | G | fest |
| `Anteil_ℓ` | Anteil der Leistung einer Winkelordnung in einer Auswahl ihrer Formen | 1 | G | fest |
| `C_ℓ/C_ℓ∞` | Leistung einer Winkelordnung in der geschlossenen Zelle gegen den unendlichen Raum, für den Beobachter heute | 1 | G | Epoche |
| `σ_ℓ` | kosmische Streuung: relative Streuung der Leistung einer Winkelordnung an einem Himmel | 1 | G | fest |
| `z_c` | Rotverschiebung, bei der ein Bereich kollabiert | 1 | LG | je System |
| `t_w` | Dauer, in der w vor der Umkehr einer Schale von 0,9 auf 0,1 fällt | s | LG | je System |
| `X_b` | Fluchtmaß je Bereich: ln des Maßes der ganzen Tori in der Naht, Regel gegen Newton am selben Rand | 1 | LZ | je System |
| `R/R_ta` | Radius einer Schale in Einheiten des Umkehrradius | 1 | LG | je System |
| `I_vor` | Übergabeintegral ∫(1 − w)ω dt einer heute umkehrenden Schale bis zur Umkehr | rad | LG | je System |
| `e_M` | Anfangsprofil eines Bereichs δM/M ∝ M^{−e} (Keim e = 1) | 1 | L | je System |
| `ε_j` | Drehimpuls der Schale ab der Umkehr, j² = ε_j·GM·r_ta | 1 | L | je System |
| `ΔΣ` | Linsensignal Σ̄(<R) − Σ(R) im Abstand R (gelesene Faltung, projiziert) | kg/m² | L | je System |
| `M_b(<r)` | eingeschlossene Baryonenmasse eines Aufbaus (Gas, Zentralgalaxie, Galaxien) | kg | L | je System |
| `M_L(<r)` | gelesene Linsenmasse innerhalb r nach einem Referenzprofil (gelesen, kein Stoff) | kg | L | je System |
| `F_H` | Rest im Haufen: gelesene Linsenmasse durch die Masse nach der örtlichen Regel, M_L/(ν·M_b) | 1 | L | je System |
| `c_NFW` | Konzentration des Referenzprofils der Linse (r₂₀₀/r_s) | 1 | L | je System |
| `n_P` | Index der Polytrope (1,5 konvektiv, 3 strahlend) | 1 | L | je System |
| `x₁_LE` | erste Nullstelle der Lane-Emden-Lösung θ(x₁) = 0 | 1 | Z | fest |
| `m₁_LE` | Massenzahl der Polytrope −x₁²θ′(x₁) | 1 | Z | fest |
| `ρ_c/ρ̄` | Zentraldichte gegen mittlere Dichte der Polytrope x₁³/(3m₁) | 1 | Z | fest |
| `f_T` | Anteil der Tiefe in kT_c = f·δ·m̄c², x₁/((n + 1)m₁) | 1 | Z | fest |
| `μ_m` | mittlere Teilchenmasse in Nukleonmassen (m̄ = μ_m·m_u; im Gas des Halos μ_m·m_p) | 1 | L | je System |
| `μ_e` | Masse je Elektron in m_u | 1 | L | je System |
| `T_c` | Kerntemperatur eines nicht entarteten Trichters, am Ort gemessen | K | L | je System |
| `M_Grenze` | Grenzmasse der Zählregel der Elektronen (Tiefe gegen hochrelativistische Zählregel) | kg | L | fest |
| `T_max` | höchste Kerntemperatur eines Trichters, den die Zählregel der Elektronen hält | K | L | je System |
| `g_T` | Decke kT_max/((μ_m/μ_e)m_ec²) an der Grenzmasse, max_x[x/4 − Π(x)] | 1 | Z | fest |
| `𝒞` | Carter-Faktor α_G/(α¹²(m_e/m_p)⁴) | 1 | Z | fest |
| `r*/R` | Lage der Schwelle zur Hülle √(g(R)/a₀); trennt kompakte und kollektive Klasse | 1 | L | je System |
| `M̃` | Spiegelmasse m_Pl²/M (vertauscht λ̄ und 3Δs) | kg | L | je System |
| `r_sp` | Spiegelort r*²/r (bildet r₀ auf Λ ab, den Raum auf die innere Richtung) | m | LI | je System |
| `α/α_G` | elektrische gegen gravitative Wechselwirkung zweier Teilchen | 1 | L | je System |
| `J` | Drehimpuls eines Körpers | J·s | L | je System |
| `μ_mag` | magnetisches Moment Q_el·c·a_J | J/T | L | je System |
| `Γ_G` | Faltungsrate α_G(M)·ω_C | s⁻¹ | L | je System |
| `ω_F` | Faltungstakt c/(3Δs) (Faltung im Zentrum, keine Raumlänge) | s⁻¹ | L | je System |
| `ω_H` | Umlauf des Elektrons im Wasserstoff αc/a_B | s⁻¹ | L | fest |
| `Δω` | Frequenzunterschied zweier Zyklen | s⁻¹ | Z | je System |
| `K_E` | Kopplung zweier Zyklen beim Einrasten | s⁻¹ | Z | je System |
| `Δθ` | Phasenlage zweier eingerasteter Zyklen | rad | Z | je System |
| `t_φ` | thermische Dephasierungszeit ħ/(k_BT) | s | L | je System |
| `r_ph` | Photonensphäre 9Δs | m | L | je System |
| `r̃` | Bahnradius in Einheiten 3Δs | 1 | L | je System |
| `ω_φ` | azimutale Grundfrequenz einer Kerr-Bahn | s⁻¹ | L | je System |
| `ω_r` | radiale Grundfrequenz einer Kerr-Bahn | s⁻¹ | L | je System |
| `ω_θ` | vertikale Grundfrequenz einer Kerr-Bahn | s⁻¹ | L | je System |
| `ζ` | Knotendrehung der Bahnebene (ω_φ − ω_θ)/ω_φ | 1 | L | je System |
| `ϑ` | Polwinkel gegen die Drehachse | rad | L | je System |
| `r_E` | Nullfläche des Takts ruhender Uhren 3Δs + √(9Δs² − a_J²cos²ϑ) | m | L | je System |
| `r_±` | Horizonte 3Δs ± √(9Δs² − a_J²) | m | L | je System |
| `ΔT_pr` | Umlaufzeit mitläufig gegen gegenläufig 4π·a_J/c | s | L | je System |
| `cos ψ` | Kosinus des Winkels zwischen Drehachse und Ort | 1 | L | je System |
| `Ω_LT` | Mitführung eines Kreisels (Frame-Dragging), Betrag | s⁻¹ | L | je System |
| `δ_L` | Winkel des Leitsterns zur Bahnebene (Deklination bei polarer Bahn) | rad | L | je System |
| `Ω_LS` | Mitführung vom Kreisel zum Leitstern gesehen, Ω_LT·cos δ_L | s⁻¹ | L | je System |
| `Ω_geo` | geodätische Drehung eines Kreisels auf Kreisbahn | s⁻¹ | L | je System |
| `P_GW` | Quadrupol-Leistung einer Kreisbahn | W | L | je System |
| `f_gw` | Frequenz einer Gravitationswelle im Detektor | Hz | L | je System |
| `ḟ_gw` | Anstieg der Frequenz (Chirp) | Hz/s | L | je System |
| `Δs_c` | Faltung der Chirp-Masse G·M_c/(3c²) | m | L | je System |
| `F_j` | Eigenfrequenz des Trichters in Einheiten c/(2π·3Δs), hängt nur am Endspin | 1 | L | je System |
| `f_RD` | Ringdown-Frequenz c/(2π·3Δs)·F(j) | Hz | L | je System |
| `E_rad` | abgestrahlte Energie einer Quelle | J | L | je System |
| `Δs_rad` | Faltung der abgestrahlten Energie G·E_rad/(3c⁴) | m | L | je System |
| `D_L` | Leuchtkraftabstand einer Quelle | m | G | Epoche |
| `κ_Θ` | Richtungsfaktor des Gedächtnis-Effekts | 1 | L | je System |
| `h_mem` | bleibende Verschiebung nach einer Gravitationswelle | 1 | L | je System |
| `Φ_X` | Zusatzpotential des kollektiven Trichters v_flach²·ln ν(s) | m²/s² | LI | je System |
| `x_F` | Fermi-Impuls der Elektronen in m_ec, λ̄_e·(3π²n_e)^{1/3} | 1 | L | je System |
| `Π_e` | Druck der Zählregel je Elektron in m_ec², P_e/(n_e m_ec²) | 1 | L | je System |
| `∂P/∂T` | Änderung der eingespeisten Heizleistung mit der Temperatur (Schleife Feld → Strom → Wärme → Widerstand) | W/K | L | je System |
| `α_s` | starke Kopplung bei der Skala μ_R (läuft mit der Skala) | 1 | L | je System |
| `μ_R` | Energieskala der laufenden Kopplung | GeV | L | je System |
| `M_Z` | Masse des Z-Bosons (Bezugsskala der starken Kopplung) | GeV | L | fest |
| `m_q` | Quarkmasse, Schwelle des Laufs (m_c, m_b, m_t) | GeV | L | fest |
| `n_f` | Zahl der leichten Quarksorten unter der Skala | 1 | L | je System |
| `q_QCD` | Nome der starken Kopplung Λ_QCD/M = e^{−2π/(β₀α_s)} (im Modelltext \|q̃\|) | 1 | L | je System |
| `g_D/e` | magnetische gegen elektrische Elementarladung 1/(2α) (Dirac-Bedingung) | 1 | L | fest |
| `r_Qel` | Radius des Ladungskreises 2ℓ_P/√α | m | L | fest |
| `E_Qel` | Energie eines Umlaufs um den Ladungskreis ħc/r_Qel | J | L | fest |
| `m_A` | Masse eines Atoms | kg | L | je System |
| `δ_Z` | Zustandsschwelle 2k_BT/(m_Ac²) = 2(v_th/c)² | 1 | L | je System |
| `f_s` | Oberflächentemperatur in Einheiten α²m_ec²: k_BT_eff/(α²m_ec²) | 1 | L | je System |
| `β_E` | Anteil des Gasdrucks am Druck (Eddington-Standardmodell, Polytrope n = 3) | 1 | L | je System |
| `M_β` | Masse eines Sterns mit Gasdruckanteil β_E im Eddington-Standardmodell | kg | L | je System |
| `S_F` | Entropie der Falte 36π(Δs/ℓ_P)² in Einheiten von k_B | 1 | L | je System |
| `ΣS_F` | Summe der Entropien aller Falten im Horizont, in Einheiten von k_B | 1 | G | Epoche |
| `t_P/t_F` | Page-Zeit als Anteil der Lebensdauer der Falte | 1 | Z | fest |
| `z_f` | Rotverschiebungsfaktor a₂/a₁ der Strahlung bis zur Rückverdichtung | 1 | G | je System |
| `N_RV` | Mindestzahl der Falten, deren Strahlung zu einer neuen Falte zurückverdichtet | 1 | Z | je System |
| `P_abs` | absorbierte Strahlungsleistung eines Körpers | W | L | je System |
| `A_alb` | Albedo (zurückgeworfener Anteil der einfallenden Strahlung) | 1 | L | je System |
| `T_A` | Abstrahlungstemperatur eines Körpers im Strahlungsgleichgewicht | K | L | je System |
| `J_S` | Entropieexport eines Körpers im Strahlungsstrom | W/K | L | je System |
| `a_Vol` | Volumenkoeffizient des Tröpfchenmodells | J | L | fest |
| `a_Obf` | Oberflächenkoeffizient des Tröpfchenmodells (E_S = a_Obf·A^{2/3}) | J | L | fest |
| `a_Cou` | Coulomb-Koeffizient des Tröpfchenmodells (E_C = a_Cou·Z²/A^{1/3}) | J | L | fest |
| `a_Asy` | Asymmetriekoeffizient des Tröpfchenmodells | J | L | fest |
| `a_Paar` | Paarungskoeffizient des Tröpfchenmodells | J | L | fest |
| `x_Sp` | Spaltbarkeit E_C/(2E_S) = a_Cou·Z²/(2a_Obf·A); Kippstelle bei 1 | 1 | L | je System |
| `Q_RD` | Güte des Ringdowns f_RD/(2·Dämpfungsrate), hängt nur am Endspin | 1 | L | je System |
| `τ_RD` | Abklingzeit des Ringdowns Q_RD/(πf_RD) | s | L | je System |
| `ι` | Neigung der Bahnachse einer Wellenquelle gegen die Sichtlinie (Kante: π/2) | rad | L | je System |
| `ε_r` | relative Permittivität eines Mediums | 1 | L | je System |
| `I_ion` | Ionenstärke einer Lösung ½Σc_iz_i² | mol/m³ | L | je System |
| `λ_D` | Debye-Länge √(ε₀ε_rk_BT/(2N_Ae²I_ion)) | m | L | je System |
| `T_R` | Rückkehrzeit in der inneren Richtung (Variable des Kerns, keine Zeit längs der Bahn) | s | I | je System |
| `K(T)` | Kern im Zeitbereich: Dichte erster Rückkehrzeiten mit dem Verlust Γ | s⁻¹ | I | je System |
| `Ω` | Winkelgeschwindigkeit der Kreisbahn, Ω² = a/r | s⁻¹ | L | je System |
| `\|K̂(iΩ)\|/K̂(0)` | Anteil des Kerns bei Bahnfrequenz gegen den stationären Kern | 1 | I | je System |
| `k_ϖ` | Wellenzahl einer Mode des Feldes am Ort | m⁻¹ | L | je System |
| `P_ϖ` | Glättung einer Mode über die Tiefe ϖ, e^{−\|k\|ϖ} (Halbgruppe P_a·P_b = P_{a+b}) | 1 | I | je System |
| `λ_Λ` | kleinster Eigenwert der bei Λ geschlossenen Rückkehr, Dπ²/(4Λ²) | s⁻¹ | I | je System |
| `Λ/ϖ_R` | Schließlänge gegen Rückkehrtiefe | 1 | I | je System |
| `q_Λ` | Rückkehranteil der bei Λ geschlossenen Rückkehr cosh(s(Λ/ϖ_R − 1))/cosh(sΛ/ϖ_R) | 1 | I | je System |
| `kΛ` | Schließlänge in Einheiten der Abschirmlänge, Λ·√(Γ/D) | 1 | I | je System |
| `ε_Λ` | Abweichung der Regel von Tully-Fisher durch den Schluss bei Λ, 1/tanh(kΛ) − 1 | 1 | I | je System |
| `r_Λ` | Radius, an dem Γ auf den kleinsten Eigenwert der geschlossenen Rückkehr fällt | m | LI | je System |
| `T_I` | Periode 2Λ/c der inneren Richtung als Temperatur geschrieben, kT_I = ħa₀/c (keine Wärme) | K | I | seit z≈10 |
| `x_V` | Argument des Strukturteils im Energiefunktional, \|∇Φ̄_N\|²/a₀² = s⁴ | 1 | L | je System |
| `Ē^TF` | spurfreie Hesse-Form des geglätteten Elternfelds (3×3), gemittelt über die eigene Zone | s⁻² | L | je System |
| `M×` | Masse, an der die Glättungsgrenze ℓ = R die Schwelle r* = R schneidet, a₀³/(G𝒯²) | kg | L | je System |
| `R×` | Radius dieses Schnitts, a₀/𝒯 | m | L | je System |
| `ℓ/r*` | Glättungslänge gegen Schwellenradius, (M×/M)^{1/6} | 1 | L | je System |
| `δ*` | Tiefe an der Schwelle (v_flach/c)², Schritt der Leiter 3Δs : r* : ℓ₀ | 1 | LI | je System |
| `d_eff` | Flussdimension einer Form der Dimension d unter der Regel 1 + (d − 1)(1 − χ/2); Kugel 1 − d ln a/d ln r = 3 − χ | 1 | L | je System |
| `κ_E` | Epizykelfrequenz der Kreisbahn | s⁻¹ | L | je System |
| `ω_z` | Vertikalfrequenz der Kreisbahn (im Modelltext ν_z) | s⁻¹ | L | je System |
| `(κ_E²+ω_z²)/Ω²` | Summe der Epizyklen im Maß der Umlauffrequenz, 2 + χ | 1 | L | je System |
| `ℬ_d` | Stufe d der Bose-Leiter s^d·n_B(s) | 1 | Z | je System |
| `ℐ_d` | Integral der Stufe d über s, Γ(d+1)ζ(d+1) | 1 | Z | fest |
| `g_F` | Fisher-Metrik der Rücklaufzahl −dν/ds = ν(ν − 1) | 1 | Z | je System |
| `S_n` | Entropie der Rücklaufzahl ln ν + χ (im Modelltext H) | 1 | Z | je System |
| `a_F` | Parameter der Falte in der Normalform ẋ = −(x³ − a_F·x − c) | 1 | Z | je System |
| `δ_F` | Abstand des Arbeitspunkts zur Falte | 1 | Z | je System |
| `k_K` | Stärke der Kopplung zweier Phasen | s⁻¹ | Z | je System |
| `φ_K` | Verzögerungswinkel der Kopplung | rad | Z | je System |
| `Δω_R` | größter Frequenzabstand, bei dem zwei Zyklen einrasten, 2k\|cos φ\| | s⁻¹ | Z | je System |
| `u_P` | Folge der Eingangswerte eines Preisach-Gedächtnisses | 1 | Z | je System |
| `α_P` | Aufwärtsschwelle je Hysteron | 1 | Z | je System |
| `β_P` | Abwärtsschwelle je Hysteron | 1 | Z | je System |
| `μ_P` | Gewicht je Hysteron (Preisach-Maß) | 1 | Z | je System |
| `y_P` | Ausgang Σ μ_P·Zustand des Preisach-Gedächtnisses | 1 | Z | je System |
| `n_R` | Glied der RC-Kette, dessen Spannung auf den Eingang zurückgeführt wird | 1 | Z | je System |
| `RG` | Längswiderstand mal Querleitwert je Glied der RC-Kette | 1 | Z | je System |
| `k_w` | Gewicht einer Modulform (gerade, ≥ 4) | 1 | Z | fest |
| `G_k` | Eisenstein-Reihe Σ′(m + nτ)^{−k}, 𝒫 auf dem Torus angewandt auf z^{−k} | 1 | Z | je System |
| `Δf_d` | Verstimmung der Anregung gegen die Eigenfrequenz | Hz | L | je System |
| `V_D` | Amplitude der Antriebsspannung | V | L | je System |
| `ΔT_D` | Erwärmung des heißen Asts im stationären Punkt | K | L | je System |
| `Δf_F` | Verstimmung, an der der heiße Ast endet (Falte) | Hz | L | je System |
| `Δf_H` | Abstand der Hopf-Grenze vor der Falte des heißen Asts | Hz | L | je System |
| `f_ring` | Ringfrequenz an der Hopf-Grenze | Hz | L | je System |
| `τ_th/τ_m` | Wärmezeit C_th/G_th gegen Ringzeit 2Q_l/ω₀ der Amplitude | 1 | L | je System |
| `β_Page` | Parameter der Page-Kurve der Abstrahlung einer Falte (Anteil der Emission) | 1 | L | fest |
| `γ_K` | Kopplung des Maxwell-Spielkerns (Feder γ_K, Dämpfer γ_K·t_K) | 1 | Z | je System |
| `v_K` | Anfangsgeschwindigkeit im Spielbeispiel (Einheiten GM = 1, x₀ = 1) | 1 | Z | je System |
| `Δh` | Höhenunterschied im Vergleich der beiden Kanäle | m | L | je System |
| `D_Z` | Dissipation D(Z) der Torus-Dynamik T² × S (transversale Rate −D_Z) | 1 | Z | je System |
| `ℓ₁` | erster Lyapunov-Koeffizient der Hüllkurve an der Hopf-Grenze (Koordinaten Re X/\|X₀\|, Im X/\|X₀\|, ΔT/K); > 0 unterkritisch, < 0 überkritisch | s⁻¹ | L | je System |
| `q_σ` | Anteil der Schleifenverstärkung am Draht über σ_el(T) allein | 1 | L | je System |
| `q_E` | Anteil der Schleifenverstärkung am Draht über E(T) allein | 1 | L | je System |
| `q_N` | Anteil der Schleifenverstärkung am Draht über die Wärmespannung N(T) allein | 1 | L | je System |
| `α_D` | Drehzahl einer irrationalen Drehung des Kreises (golden: (√5 − 1)/2) | 1 | Z | fest |
| `K_M` | Zahl der Fluchtmoden je Umlaufsinn (Wellenzahlen 1 … K_M) | 1 | Z | je System |
| `ν_F` | größte Verstärkung der Fluchtmoden 1/\|1 − q·e^{2πikα_D}\| der Resolvente (I − qT)⁻¹ | 1 | Z | je System |
| `x_U` | Umkehrpunkt des Spielbeispiels mit Maxwell-Spielkern (Einheiten GM = 1, x₀ = 1) | 1 | Z | je System |
| `c_Δ` | Koeffizienten der q-Reihe der Diskriminante Δ = Σ c_Δ(n)qⁿ (Ramanujan) | 1 | Z | fest |
| `Δ_τ` | Diskriminante Δ(τ) = q∏(1 − qⁿ)²⁴ des Torus ℂ/(ℤ + τℤ) | 1 | Z | je System |
| `λ_ψ` | Rate am Rastpunkt der Phasendifferenz, −√((2k cos φ)² − Δω²) | s⁻¹ | Z | je System |
| `N_K` | Zahl gekoppelter Zyklen an einem gemeinsamen Medium (Kuramoto) | 1 | Z | je System |
| `(N−1)k₁/W` | Kopplung beim ersten Einrasten im Maß der Bandbreite W, mal N − 1 | 1 | Z | je System |
| `λ_⊥` | transversale Rate des Torus T² × S, −D_Z | 1 | Z | je System |
| `G_d(1)` | größte Verstärkung der Gitter-Rückkehr in d Dimensionen, G_d(1) = 1/(1 − p_d) | 1 | I | fest |
| `T_μ` | Rückkehrzeiten (Schritte) einer irrationalen Drehung in einen Bogen vom Maß μ(A) | 1 | I | je System |
| `N_H` | Zahl der Harmonischen des Dreiecks in der Windungskorrelation | 1 | I | je System |
| `S_CHSH` | größte CHSH-Größe einer drehinvarianten Korrelation E(Δ) | 1 | I | je System |
| `Ṁ_F` | Rückgabe der Falte dM/dt | kg/s | L | je System |
| `R_K` | Verhältnis differentieller Kanal (Vorspannung → Resonanz) gegen Takt-Kanal bei radialer Verschiebung | 1 | L | je System |
| `𝔄_VO₂` | Satz der Annahmen des VO₂-Aufbaus (Übergang, Widerstände, Wärme, Kreis, Brücke, Feld; Einheiten µs, K, V, kΩ, mA, nF, nJ/K, µm) | gemischt | L | je System |
| `V_H` | Quellspannung an der Hopf-Grenze des VO₂-Elements | V | L | je System |
| `f_H` | Frequenz des geborenen Zyklus an der Hopf-Grenze | Hz | L | je System |
| `g_E` | Entdämpfung (negativer Widerstand) eines selbsterregten Umlaufs | Ω | L | je System |
| `b_I` | kubische Sättigung des Umlaufs (Glied b_I·I³) | Ω/A² | L | je System |
| `dR/dT` | Temperaturgang des Widerstands | Ω/K | L | je System |
| `V_obs` | gemessene Rotationsgeschwindigkeit | m/s | L | je System |
| `e_V` | Messfehler von V_obs | m/s | L | je System |
| `V_gas` | Beitrag des Gases (mit Helium) zur Rotation | m/s | L | je System |
| `V_disk` | Beitrag der Sternscheibe zur Rotation bei Υ = 1 | m/s | L | je System |
| `V_bul` | Beitrag des Bulges zur Rotation bei Υ = 1 | m/s | L | je System |
| `V_f` | flache Rotationsgeschwindigkeit, gemessen | m/s | L | je System |
| `Υ_d` | Masse-zu-Licht der Scheibe bei 3,6 μm | M☉/L☉[3,6] | L | je System |
| `Υ_b` | Masse-zu-Licht des Bulges bei 3,6 μm | M☉/L☉[3,6] | L | je System |
| `R_d` | Skalenlänge der Scheibe | m | L | je System |
| `SB₀` | zentrale Flächenhelligkeit der Scheibe bei 3,6 μm | L☉[3,6]/pc² | L | je System |
| `L_3,6` | Leuchtkraft bei 3,6 μm | L☉[3,6] | L | je System |
| `M_HI` | Masse des neutralen Wasserstoffs (im Datensatz in M☉) | kg | L | je System |
| `i_inc` | Neigung (Scheibe gegen die Himmelsebene, Ring gegen die Scheibe) | rad | L | je System |
| `Q_fl` | Güte einer SPARC-Kurve (1 hoch, 2 mittel, 3 niedrig) | 1 | L | je System |
| `Z_A` | Kernladungszahl | 1 | L | fest |
| `A_A` | Massenzahl | 1 | L | fest |
| `E_B/A` | Bindungsenergie je Nukleon | J | L | fest |
| `Δ_M` | Massenexzess | J | L | fest |
| `T_schm` | Schmelztemperatur bei Normaldruck | K | L | fest |
| `T_sied` | Siedetemperatur bei Normaldruck | K | L | fest |
| `E_ion` | erste Ionisierungsenergie | J | L | fest |
| `V_bar` | baryonische Rotation, V_bar² = V_gas\|V_gas\| + Υ_d·V_disk² + Υ_b·V_bul² | m/s | L | je System |
| `g_obs` | gemessene Zentripetalbeschleunigung (V_obs² + σ_v²·R/R_d)/R, mit dem Druck des Gases | m/s² | L | je System |
| `ε_R` | Rest je Galaxie: Median von log₁₀(g_obs/(ν(g_N)·g_N)) in dex | 1 | L | je System |
| `σ_ε` | Streuung einer Lesung aus den Messfehlern von Abstand, Neigung und Υ (dex) | 1 | L | je System |
| `a₀_gal` | Schwelle je Galaxie gelesen (Rest null) | m/s² | LI | je System |
| `a₀_prof` | Schwelle aus allen SPARC-Kurven zugleich, Abstand, Neigung und Υ je Galaxie als Störgrößen | m/s² | LI | je System |
| `a₀_pkt` | Schwelle an den SPARC-Punkten, Bose-Form an Bin-Medianen | m/s² | LI | je System |
| `Ȳ_bin` | Median einer Größe je Bin in log g_N (Rechenhilfe) | 1 | L | je System |
| `F_q` | Fehlerquadratsumme einer Anpassung an Bin-Mediane | 1 | L | je System |
| `d_s` | Spektraldimension der Rückkehr, 2·d ln a/d ln a_N | 1 | LI | je System |
| `Σ` | Flächendichte einer Scheibe am Ort | kg/m² | L | je System |
| `Σ₀` | zentrale Flächendichte einer Scheibe (SPARC: Υ_d·SB₀) | kg/m² | L | je System |
| `Σ_dyn` | zentrale dynamische Flächendichte (Toomre-Integral) (1/2πG)∫V²/R² dR | kg/m² | L | je System |
| `g_R` | Newton-Feld in der Scheibenebene, radial (Betrag) | m/s² | L | je System |
| `g_z` | Newton-Feld senkrecht zur Ebene (Betrag; an der Scheibe 2πGΣ) | m/s² | L | je System |
| `a_R` | radiale Beschleunigung in der Ebene nach der Regel | m/s² | L | je System |
| `f_Σ` | Gewicht des senkrechten Scheibenfelds in der Regel (Feldform 1, Bahnregel in der Ebene 0) | 1 | L | je System |
| `b_K` | Skala der Kuzmin-Scheibe | m | L | je System |
| `a_MN` | radiale Skala der Miyamoto-Nagai-Scheibe | m | L | je System |
| `b_MN` | vertikale Skala der Miyamoto-Nagai-Scheibe | m | L | je System |
| `ω_p` | Knotenpräzession eines geneigten Rings | s⁻¹ | L | je System |
| `Ω_N` | Kreisfrequenz nach Newton √(a_N/r), Punktmasse √(GM/r³), von der fernen Uhr gezählt | s⁻¹ | L | je System |
| `ω_z,N` | Vertikalfrequenz nach Newton | s⁻¹ | L | je System |
| `ρ` | Dichte der Baryonen am Ort | kg/m³ | L | je System |
| `ρ_X` | gelesene Dichte eines Beobachters, der Newton erwartet (Ausgabe seiner Rechnung, kein Stoff) | kg/m³ | L | je System |
| `A` | Massenmittel der Regel über das System ⟨ν(\|g_i + g_e\|)(g_i + g_e)⟩_M (Betrag, parallel zu g_e) | m/s² | L | je System |
| `A_m1` | Amplitude der Einseitigkeit m = 1 der Rotation (Anteil am Mittel) | 1 | L | je System |
| `R_e` | Halblichtradius, Skala der Plummer-Kugel | m | L | je System |
| `A_O` | Oort-Konstante A = −½R·dΩ/dR | s⁻¹ | L | je System |
| `B_O` | Oort-Konstante B = A − Ω | s⁻¹ | L | je System |
| `𝒯_H` | Gezeit des Hill-Radius, 4Ω² − κ² = 4A_O(A_O − B_O) (ungeglättet, andere Größe als 𝒯) | s⁻² | L | je System |
| `μ_fil` | Linienmasse eines Filaments | kg/m | L | je System |
| `r*_fil` | Schwelle eines Filaments 2Gμ/a₀ | m | L | je System |
| `Δs_eff` | gezählte Faltung ν·Δs (Faltung mal Rückkehrzahl; die Regel faltet nicht) | m | L | je System |
| `M_b` | baryonische Masse eines Systems, aus der Gastemperatur gelesen | kg | L | je System |
| `S_T` | logarithmisches Gefälle von Dichte und Temperatur, f = 1/S_T | 1 | L | je System |
| `f_He` | Gasmasse gegen M_HI (Helium), M_gas = f_He·M_HI | 1 | L | fest |
| `σ_Υ` | Streuung von Υ_d je Galaxie im Fehlermodell (dex) | 1 | L | fest |
| `Ṗ/P` | relative Änderung der Pulsperiode durch die Beschleunigung längs der Sichtlinie | s⁻¹ | L | je System |
| `t_W` | Zeit, bis die Knotenlinie des Warps bei gemessener Rate um 90° dreht | s | L | je System |
| `I_F` | Übergabeintegral ∫ω dt bis zur ersten Faltung des Phasenraumblatts | rad | L | je System |
| `t_T` | Zeit, bis sich Linsenmaximum und Gas beim Haufenstoß trennen | s | L | je System |
| `σ_v,N` | Geschwindigkeitsstreuung mit Newton allein (ohne Regel) | m/s | L | je System |
| `V_σ` | Verstärkung der Streuung (σ_v/σ_v,N)² druckgestützter Systeme | 1 | L | je System |
| `x_gl` | Lage des Gleichstands von Feldform und radialer Regel in der Scheibenebene, in R_d | 1 | L | je System |
| `ε_F` | gelöste Feldform gegen die örtliche Näherung ν(√(g_R² + (2πGΣ)²))·g_R in der Ebene | 1 | L | je System |
| `a₀_gem` | Schwelle gemeinsam je Drittel der Flächendichte angepasst (Υ_d fest) | m/s² | LI | je System |
| `a₀_e` | Schwelle, die man ohne äußeres Feld liest, wo die Regel mit äußerem Feld wirkt | m/s² | LI | je System |
| `k_a` | Anstieg von log₁₀ a₀ je Dekade der Flächendichte Σ₀ | 1 | L | je System |
| `ε_Σ` | Rest der zentralen Flächendichte log₁₀(Σ_dyn,gem/Σ_dyn,Regel) (dex) | 1 | L | je System |
| `k_Σ` | Anstieg des Rests ε_Σ je Dekade der Flächendichte Σ₀ | 1 | L | je System |
| `k_gas` | Änderung des Rests der Regel mit dem örtlichen Gasanteil innerhalb der Galaxien (dex) | 1 | L | je System |
| `ε_Q` | Rest der Regel an den Punkten jenseits R_Q (a_N < (ln q_kosm)²a₀), dex | 1 | L | je System |
| `k_aus` | Steigung des Rests gegen log a_N an den Außenpunkten (a_N < 0,1·a₀) | 1 | L | je System |
| `σ_ds` | Fehler der Spektraldimension aus Ziehungen der SPARC-Kurven | 1 | LI | je System |
| `B_a` | Boden der Regel a = ν(a_N)·a_N + B, in a₀ | 1 | LI | je System |
| `Pe` | Péclet-Zahl einer Drift der Rückkehr (Drift gegen Diffusion) | 1 | I | je System |
| `b_q` | Bandbreite einer kontrahierenden Querrichtung der Rückkehr | m | I | je System |
| `t_Wick` | Zeit, bis die Knotenlinien zweier Ringe um Δφ auseinanderlaufen (Aufwicklung des Warps) | s | L | je System |
| `𝒯_Q` | Quadrupol des Restfelds der Regel um einen geglätteten Körper im Feld der Elternebene | s⁻² | L | je System |
| `Σ_z` | gelesene Flächendichte der Säule \|z\| < z₁ am Ort, a_z(z₁)/(2πG) (Ausgabe einer Lesung, kein Stoff) | kg/m² | L | je System |
<!-- ende -->

## 5 Formeln je Ebene

Kennung, Gleichung und Herkunft der Formeln, die im Modelltext stehen (Formeln der Forschung und Rechenhilfen sind im Register geführt, aber nicht hier).

<!-- erzeugt: formeln -->
### Lokal (L), 171 Formeln

| Kennung | Gleichung | Herkunft |
|---|---|---|
| verlust_gamma | Γ = a_N/c | Baustein 07 |
| s_ | s = √(Γ·τ_D) = √(a_N/a₀) | Baustein 07 |
| q_ | q = e^{−√(a_N/a₀)} | Autor (FTM-PDF) |
| nu | ν = 1/(1 − e^{−s}) = 1 + n_B(s) | Autor (FTM-PDF) |
| chi | χ = s/(e^s − 1) | Autor (FTM-PDF) |
| n_B | n_B = 1/(e^s − 1) | Baustein 07 |
| Q_tilde | Q̃(s) = −4∫_s^∞ x³n_B(x)dx | Baustein 07 |
| radius_schluss | r_Λ = (2/π)·Λ·r*/ϖ_R | Baustein 07 |
| strukturteil | Q_s(x) = ∫₀^x (ν(√x′) − 1) dx′, x = s⁴ | Baustein 07 |
| neigung_xy | g = x/y²·g⊕ (A1) | Autor (FTM-PDF, A1) |
| ds | Δs = GM/(3c²) (A2) | Autor (FTM-PDF, A2) |
| neigung | g = GM/R² = 3c²Δs/R² | Autor (FTM-PDF) |
| takt | dτ/dt = √(1 − 6Δs/r) (A3) | Autor (FTM-PDF, A3) |
| delta | δ = 1 − dτ/dt | Autor (FTM-PDF) |
| r0 | r₀ = 6Δs | Autor (FTM-PDF) |
| r_stern | r* = √(GM/a₀) | Autor (FTM-PDF) |
| r_stern_leiter | r*² = 3Δs·ℓ₀ | Baustein 08 |
| leiter_schritt | δ* = 3Δs/r* = r*/ℓ₀ = (v_flach/c)² | Baustein 08 |
| gezeit_betrag | 𝒯 = √(3/2)·\|Ē^TF\| | Baustein 07 |
| glaettungslaenge | ℓ³ = GM/𝒯 | Baustein 07 |
| gezeitendichte | ρ_J = 3𝒯/(4πG) | Baustein 07 |
| schnittmasse_glaettung | M× = a₀³/(G𝒯²) | Baustein 07 |
| schnittradius_glaettung | R× = a₀/𝒯 | Baustein 07 |
| glaettung_zu_schwelle | ℓ/r* = (M×/M)^{1/6} | Baustein 07 |
| flussdimension | d_eff = 1 + (d − 1)(1 − χ/2); Kugel (d = 3): 1 − d ln a/d ln r = 3 − χ | Baustein 03, 07 |
| epizyklen_summe | (κ_E² + ω_z²)/Ω² = 2 + χ | Baustein 07 |
| kanal_verhaeltnis | R_K = f_g·(r − 6Δs)/(3Δs) | Baustein 05 |
| vo2_hopf | Re λ = 0 des elektrothermischen Paars der Jacobi-Matrix am Fixpunkt des VO₂-Elements | Baustein 06 |
| kupfer_hopf_abstand | Δf_H = Abstand der Hopf-Grenze vor der Falte des heißen Asts | Baustein 06 |
| kupfer_ringfrequenz | f_ring = \|Im λ\|/2π an der Hopf-Grenze | Baustein 06 |
| kupfer_zeitverhaeltnis | τ_th/τ_m = (C_th/G_th)/(2Q_l/ω₀), Q_l = Q_m/(1 + Γ_B²/(γ_m R_el,0)) | Baustein 06 |
| kupfer_lyapunov | ℓ₁ = Re[⟨p, C(q,q,q̄)⟩ − 2⟨p, B(q, J⁻¹B(q,q̄))⟩ + ⟨p, B(q̄, (2iω − J)⁻¹B(q,q))⟩]/2ω an der Hopf-Grenze der Hüllkurve (Re X, Im X in \|X₀\|, ΔT in K) | Baustein 06 |
| draht_schleife | q = (∂P/∂T)/G_th am stationären Punkt, zerlegt: q_σ (nur R(T)), q_E (nur E(T)), q_N (nur Wärmespannung N(T)), die übrigen Glieder am Arbeitspunkt festgehalten | Baustein 09 |
| uebergabezeit | ω_orb² = Δ_B·H²/2, 1/λ_m = N*/ω_orb; gibt (1/ω_orb, 1/λ_m min, max) | Baustein 03 |
| compton_red | λ̄ = ħ/(mc) | Baustein 02 |
| alpha_G | α_G = GM²/(ħc) = 3Δs/λ̄ | Baustein 02 |
| sternmasse | M_* = m_Pl³/m_p² | Baustein 02 |
| grenzmasse_hawking | T_H = T ⇒ M_T = ħc³/(8πGk_BT), T_H = ħc/(24πk_BΔs) | Baustein 04 |
| eddington | P_Edd = (9/2)·m_p·c³·Δs/r_e² | Baustein 02 |
| entropie_falte | S_F = A/(4ℓ_P²) = π·r₀²/ℓ_P² = 36π·(Δs/ℓ_P)² = 4π·α_G(M) (in k_B) | Baustein 04 |
| lebensdauer_falte | t_F = 3840·S_F·Δs/c | Baustein 04 |
| rueckgabe_rate | dM/dt = −ħc⁴/(15360π·G²·M²) | Baustein 04 |
| entropieexport | J_S = (4/3)·P_abs·(1/T_A − 1/T) | Baustein 04 |
| nome_qcd | q_QCD = Λ_QCD/M = e^{−2π/(β₀α_s(M))}, β₀ = 11 − 2n_f/3 | Baustein 02 |
| dirac_verhaeltnis | g_D/e = 1/(2α) (Dirac-Bedingung e·g = ħc/2, Gauß) | Baustein 02 |
| radius_ladungskreis | r_Qel = 2ℓ_P/√α (Ladung e als Windung um einen Kreis) | Baustein 02 |
| energie_ladungskreis | E_Qel = ħc/r_Qel = (√α/2)·m_Plc² | Baustein 02 |
| zustandsschwelle | δ_Z = 2k_BT/(m_Ac²) = 2(v_th/c)² | Baustein 02 |
| f_saha | f_s = k_BT_eff/(α²m_ec²) | Baustein 02 |
| eddington_quartik | M_β = C_E·(ħc/G)^{3/2}/(μ_m m_u)²·√(1 − β_E)/β_E², C_E = 4m₁(3)·√45/π^{3/2} | Baustein 02 |
| spaltbarkeit | x_Sp = E_C/(2E_S) = a_Cou·Z²/(2a_Obf·A) | Baustein 06 |
| debye_laenge | λ_D = √(ε₀ε_rk_BT/(2N_Ae²I_ion)) | Baustein 09 |
| temperatur_kern | T_c = f·δ·μ_m·m_u·c²/k_B, δ = 3Δs/R | Baustein 02 |
| grenzmasse_zaehlregel | M_Grenze = (K_ur·(3ρ_c/(4πρ̄))^{1/3}/(f·G))^{3/2}, K_ur = (3π²)^{1/3}ħc/(4(μ_e m_u)^{4/3}) | Baustein 02 |
| temperatur_max | kT_max = max_R μ_m m_u c²·[f·δ − (m_e/(μ_e m_u))·Π(x)], x = λ̄_e(3π²n_e)^{1/3} | Baustein 02 |
| schwelle_zu_huelle | r*/R = √(GM/(a₀R²)) = √(g(R)/a₀); r*/R = 1 ⇔ Σ̄ = a₀/(πG) | Baustein 08 |
| spiegel_masse | M̃ = m_Pl²/M: vertauscht λ̄ und 3Δs, bildet M_min auf M_U ab | Baustein 08 |
| spiegel_radius | r ↦ r*²/r: bildet r₀ = 6Δs auf Λ ab (Raum auf innere Richtung) | Autor (FTM-PDF) |
| verhaeltnis_em_grav | α/α_G(m) | Baustein 09 |
| compton_takt | ω_C = c/λ̄ = Mc²/ħ | Baustein 09 |
| faltungsrate | Γ_G = α_G(M)·ω_C = GM³c/ħ² | Baustein 09 |
| faltungstakt | ω_F = c/(3Δs) | Baustein 09 |
| umlauf_wasserstoff | ω₂ = αc/a_B | Baustein 09 |
| umlauf_kepler | ω = √(GM/r³) (ferne Uhr) | Baustein 09 |
| spinzahl | j = a_J/(3Δs) = (J/ħ)/α_G(M) | Baustein 09 |
| moment_magnetisch | μ = Q_el·c·a_J (g = 2 für jede Masse, Ladung, Drehung) | Baustein 09 |
| dephasierzeit | t_φ = ħ/(k_BT) | Baustein 09 |
| schleife_thermisch | q = (∂P/∂T)/G_th | Baustein 09 |
| a_regel | a = ν(a_N)·a_N = a_N/(1 − q) | Baustein 03 |
| v_kreis_punkt | v = √(ν(a_N)·a_N·r), a_N = GM/r² | Baustein 03 |
| v_flach | v_flach⁴ = GMa₀ | Baustein 03 |
| nu_aussenfeld | ν(\|g\|), senkrecht s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e² | Baustein 03 |
| sigma_kollektiv | Σ_a₀ = a₀/(πG), Rückgabe in M☉/pc² | Baustein 03 |
| v_bar_quadrat | V_bar² = V_gas\|V_gas\| + Υ_d·V_disk² + Υ_b·V_bul² | Baustein 03 |
| g_gemessen | g_obs = (V_obs² + σ_v²·R/R_d)/R | Baustein 03 |
| rest_galaxie | ε_R = Median log₁₀(g_obs/(ν(\|g_N + g_e\|)·g_N)), g_N = V_bar²/R, g_e senkrecht | Baustein 03 |
| a0_je_galaxie | a₀_gal: ε_R(a₀) = 0 je Galaxie | Baustein 03 |
| streuung_messfehler | σ_ε = Streuung einer Lesung über Ziehungen von Abstand, Neigung und Υ | Baustein 03 |
| toomre_sigma_dyn | Σ_dyn(0) = (1/2πG)∫₀^∞ V²/R² dR, innen V ∝ R, außen V = V_letzter | Baustein 03 |
| g_expscheibe | g_R = GM·R/(2R_d³)·[I₀K₀ − I₁K₁](R/2R_d) | Baustein 03 |
| sigma_dyn_expscheibe | Σ_dyn(0) = Σ₀·∫₀^∞ ν(2πGΣ₀·f(y))·f(y) dy/y, f(y) = y[I₀K₀ − I₁K₁](y), y = R/2R_d | Baustein 03 |
| kuzmin_feld | Kuzmin an der Scheibe (z = 0⁺): g_R = GMR/(R² + b²)^{3/2}, g_z = GMb/(R² + b²)^{3/2} = 2πGΣ | Baustein 07 |
| a_scheibe_ebene | a_R = ν(√(g_R² + f·(2πGΣ)²))·g_R | Baustein 07 |
| g_miyamoto_nagai | g_R = GMR/D, g_z = GMz(a + ζ)/(ζD), ζ = √(z² + b²), D = (R² + (a + ζ)²)^{3/2} | Baustein 07 |
| vertikalfrequenz | Kernform ω_z² = ω_z,N² + (ν − 1)Ω_N², örtliche Regel ω_z² = ν·ω_z,N² | Baustein 07 |
| knotenpraezession_klein | ω_p = Ω − ω_z (Neigung → 0), Frequenzen aus dem Feld am Ort | Baustein 07 |
| praezession_ring | ω_p = mittlere Drehrate des Drehimpulses eines geneigten Kreisrings (3D-Bahn, a = ν(\|g_N\|)·g_N) | Baustein 07 |
| einseitigkeit_m1 | A_m1 = 2\|F₁\|/F₀ von v_c(φ), v_c² = −R·(ν(\|g\|)g − B)·e_R, g = −g_i·e_R + g_e in der Ebene | Baustein 03 |
| massenmittel_A_plummer | A = ⟨ν(\|g_i + g_e\|)(g_i + g_e)⟩_M, Plummer-Kugel | Baustein 03 |
| massenmittel_A_scheibe | A = ⟨ν(\|g_i + g_e\|)(g_i + g_e)⟩_M, dünne Exponentialscheibe, g_e in der Ebene | Baustein 03 |
| sigma_plummer_jeans | σ_v² = ∫r²p·w dr/∫r²ρ·w dr, p = ∫_r^∞ρ·g dr, g = ⟨ν(\|g_i + g_e\|)(g_i + g_e)⟩ über die Kugel, Fenster R < R_e | Baustein 03 |
| sigma_v_tief | σ_v⁴ = (4/81)·GMa₀ | Baustein 03 |
| gezeit_sonnenort | 𝒯 aus diag(Φ_RR, Φ_φφ, Φ_zz): a(g_N) = v₀²/R₀, d ln g_N/d ln R = (2 d ln v/d ln R − 1)/(1 − χ/2), Φ_zz = 4πGρ − Φ_RR − Φ_φφ | Baustein 07 |
| gezeit_hill | 𝒯_H = 4Ω² − κ² = 4A_O(A_O − B_O) | Baustein 02 |
| fq_form | F_q = Σ((Ȳ_bin − log₁₀ν_Form(g_N; a₀))/σ)², Ȳ_bin Median von log₁₀(g_obs/g_N) | Baustein 07 |
| a0_an_punkten | a₀_pkt = argmin F_q(a₀), a₀ in 10^{±0,5}·10⁻¹⁰ | Baustein 07 |
| nu_ds | ν_{d_s} = 1/(1 − e^{−(a_N/a₀)^{1−d_s/2}}); d_s = 1 Bose-Form, tief a ∝ a_N^{d_s/2} | Baustein 07 |
| ds_anpassung | d_s aus ln g_obs = ln(ν_{d_s}(g_N; a₀)·g_N), robust (soft_l1, Skala 0,1), a₀ frei | Baustein 07 |
| spektraldimension_regel | d_s = 2·d ln a/d ln a_N = 2 − χ | Baustein 03 |
| dichte_gelesen | 4πGρ_X = χΩ² + 4πGρ·(ν(1 − χ/2) − 1), Ω² = νa_N/r | Baustein 07 |
| r_stern_filament | r*_fil = 2Gμ/a₀ | Baustein 03 |
| faltungsmesser | Δs = v_flach⁴/(3c²a₀) | Baustein 03 |
| ds_eff | Δs_eff = ν·Δs = a·r²/(3c²) | Baustein 03 |
| masse_aus_temperatur | M_b = (S_T·k_BT/(μ m_p))²/(G·a₀) | Baustein 03 |
| pulsar_rate | Ṗ/P = a_r/c (Doppler der Beschleunigung längs der Sichtlinie) | Baustein 03 |
| warp_viertelzeit | t_W = (π/2)/ω_p | Baustein 03 |
| verstaerkung_druck | V_σ = (σ_v/σ_v,N)² | Baustein 03 |
| eta_scheibe | η = ν(g_R)·g_R / a_R,Feldform in der Ebene einer dünnen Scheibe bei R = x·h | Baustein 07 |
| gleichstand_scheibe | x_gl: η(x_gl·R_d) = 1 zwischen R_d und 8R_d (Feldform gleich radialer Regel) | Baustein 07 |
| feldform_naeherung | ε_F = a_R,Feldform / (ν(√(g_R² + (2πGΣ)²))·g_R) in der Ebene bei R = x·h | Baustein 07 |
| a_N_aus_kreisbahn | a_N: ν(a_N)·a_N = v²/r (die Regel rückwärts) | Baustein 03 |
| sattelpunkt | d_S = D·√M₂/(√M₁ + √M₂): GM₂/d_S² = GM₁/(D − d_S)² | Baustein 03 |
| dynamikmasse_ort | Σ_z = a_z(z₁)/(2πG): Kern ν_z² = ν_z,N² + (ν − 1)Ω_N² → Σ_b + (ν(g_R) − 1)g_R·z₁/(2πGR); örtliche Regel a_z = ν(\|g_N\|)·g_z → ν(√(g_R² + (2πGΣ_b)²))·Σ_b (Schicht) bis ν(g_R)·Σ_b (Ebene) | Baustein 07 |
| knotenpraezession_scheibe | ω_p = Ω − ω_z in der Ebene einer Exponentialscheibe mit Schicht der Höhe h: ω_z,N² = 4πGΣ/(2h) − (dv_N²/dR)/R, Ω² = ν(g_R)·Ω_N², ω_z nach Kern- oder Feldform | Baustein 07 |
| knotenpraezession_frequenzen | ω_p = √ν·Ω_N − ω_z(ω_z,N, Ω_N, ν; Form) | Baustein 07 |
| saeule_loesung | Σ_z = a_z(R, z₁)/(2πG) aus einer Lösung (Feldform oder Newton) | Baustein 07 |
| warp_aufwicklung | t_Wick = Δφ/\|ω_p(R₁) − ω_p(R₂)\| | Baustein 03 |
| quadrupol_glaettung | Q₂ = −(3/2)(g_e/r_e)·J(β, λ), β = g_e/a₀, λ = ℓ/r_e, r_e = √(GM/g_e); J = ½∫∫(ν(β\|h e_r + e_z\|) − 1)·k(ρ, θ)·ρ sin θ dρ dθ + (4/15)(−χ/2)ν/λ³, h = ρ/(ρ² + λ²)^{3/2} | Baustein 07 |
| profil_chi2 | χ²_prof = min über (u_D, δ_i, u_Υ) von Σ[(ln g_obs − u_D + 2 ln(sin i/sin(i + δ_i)) − ln ν(g_N)g_N)/σ]² + (u_D/e_D)² + (δ_i/e_i)² + (u_Υ/σ_Υ)², σ² = (2e_V/V)² + σ_int², Υ = Υ_d·e^{u_Υ} | Baustein 03 |
| a0_profiliert | a₀_prof = argmin Σ_Galaxien χ²_prof(a₀) (Parabel im Gitter ln a₀), Fehler: Bootstrap über Galaxien | Baustein 03 |
| a0_je_drittel | a₀ je Drittel von Σ₀ = Υ_d·SB₀: argmin der Summe der χ²_prof (Abstand, Neigung als Störgrößen, Υ_d fest, σ = 2e_V/V); Anstieg d log a₀/d log Σ₀ über die drei Drittel, Fehler Bootstrap über Galaxien | Baustein 03 |
| a0_senkung_aussenfeld | a₀_e = argmin Σ w·[ln(ν(\|g_N + g_e\|)·g_N) − ln(ν(g_N; a₀_e)·g_N)]² an den SPARC-Punkten | Baustein 03 |
| sigma_dyn_rest | ε_Σ = ⟨log₁₀(Σ_dyn,gem/Σ_dyn,Regel)⟩ je Galaxie (Toomre-Integral), Anstieg k_Σ gegen log₁₀Σ₀ | Baustein 03 |
| rest_gasanteil | k_gas = d(Rest − Median der Galaxie)/d f_gas über alle Punkte, f_gas = V_gas\|V_gas\|/V_bar² | Baustein 03 |
| rest_jenseits | ε_Q = Median über Galaxien des Medians von log₁₀(g_obs/(ν(g_N)g_N)) an Punkten a_N < a_Q; Fehler: Bootstrap über Galaxien ⊕ Monte Carlo (Abstand, Neigung, Υ, e_V) | Baustein 03 |
| reststeigung_aussen | k_aus = Σ XY/Σ X², X, Y: log₁₀(a_N/a₀) und Rest je Galaxie zentriert, Punkte a_N < 0,1·a₀ der Galaxien mit äußerstem a_N ≤ a_rand; Regel mit äußerem Feld ν(\|g_N + g_e\|) | Baustein 03 |
| ds_streuung | σ(d_s) = Streuung von d_s (ds_anpassung) über Ziehungen von Galaxien, Abstand, Neigung und Υ | Baustein 07 |
| boden_anpassung | B_a: ln g_obs = ln(ν(g_N)·g_N + B_a·a₀) an den Punkten s < s_max, robust (soft_l1, 0,1) | Baustein 07 |
| boden_grenze | B_95 = 95-%-Quantil von B_a über Ziehungen (mindestens 0) | Baustein 07 |
| gewicht_f | f̂ = argmin_f min_{a₀} Σ_Galaxien χ²(f, a₀) (Parabel im Gitter), Fehler: Jackknife über Galaxien | Baustein 07 |
| r_Q | R_Q = r*/(−ln q_kosm) | Baustein 03 |
| masse_bQ | M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴) | Baustein 03 |
| dichte_umkehr_heute | 1 + Δ_ta aus Top-Hat im flachen Hintergrund mit Vakuumteil 1 − Ω_m | Baustein 03 |
| umkehrradius_heute | R_ta = (3M_ges/(4πρ_m(1 + Δ_ta)))^{1/3} | Baustein 03 |
| schalen_dauer | Dauer von w = 0,9 bis w = 0,1 vor der Umkehr: Ṙ² = 2μ(1/R − 1/R_ta) − Ω_Λ(R_ta² − R²), θ = 3Ṙ/R, w = 1 − e^{−√(θ/θ_c)}, Kollaps (R = 0,01·R_ta) bei z_c | Baustein 03 |
| gewicht_schale_aussen | w = 1 − e^{−√(θ/θ_c)}, θ = 3Ṙ/R heute aus Ṙ² = 2μ/R + Ω_Λ R² − K mit Alter t₀, 2μ/Ω_m = (1 + Δ_ta)(R_ta/R)³ (Punktmasse, Einheiten H₀ = 1, R = 1 heute) | Baustein 03 |
| rand_schale_aussen | R/R_ta, an dem gewicht_schale_aussen den Wert w erreicht (Umkehrung) | Baustein 03 |
| uebergabe_vor_umkehr | I_vor = ∫₀^{R_ta}(1 − w)·ω dR/Ṙ, ω² = GM/R³, θ = 3Ṙ/R (Schale, die heute umkehrt) | Baustein 03 |
| fluchtmass_bereich | X_b = ln(μ_Regel/μ_Newton), μ = ∫4πr²·(4π/3)v_R³/(1 − r²/R_ta²)dr, v_R² = 2∫_r^{R_ta} a dr′ | Baustein 03 |
| haufenstoss_trennzeit | t_T = I_F/ω_orb (Übergabe beim Gas erst nach I_F Radiant seit dem Stoß) | Baustein 03 |
| uebergabe_faltung | I_F = ∫₁^{τ_F} ω dτ, ω² = (π²/8)/Λ³ (Masse innen fest bis zur Faltung); Λ″ = −(π²/8)τ^p𝓜(Λτ^{−s})/Λ² + (π²/8)ε_j/Λ³, s = 2/(9e) + 2/3, p = 2/(3e), 𝓜 aus der Bahn selbst; τ_F: Λτ^{−s} = λ_c (äußerster Rücklauf) | Baustein 03 |
| linse_punktmasse | ΔΣ(R) = M/(πR²) (Linsensignal einer Punktmasse im Abstand R) | Baustein 03 |
| haufen_rest | F_H = M_L(<r)/(ν(GM_b/r²)·M_b(<r)) (gelesene Linse gegen die örtliche Regel allein) | Baustein 03 |
| ablenkung | α_abl = 12Δs/b | Baustein 08 |
| ablenkung_kollektiv | α_∞ = 6πΔs/r* = 2π(v_flach/c)² | Baustein 08 |
| perihel_pro_umlauf | Δφ_P = 18πΔs/(a_b(1 − e_B²)) | Baustein 08 |
| schatten_winkel | α_S = 2b/D_A, b = 3√27·Δs = 3√3·GM/c² | Baustein 08 |
| isco_kerr | r̃_ISCO = 3 + z₂ ∓ √((3 − z₁)(3 + z₁ + 2z₂)) (Kerr, Bardeen), Einheit 3Δs | Baustein 08 |
| schwellenfrequenz_gw | f*_gw = a₀^{3/4}/(π(3c²Δs)^{1/4}) | Baustein 08 |
| hulse_taylor | Ṗ_b = −(192π/5)·η_m·(2π·3Δs/(cP_b))^{5/3}·f(e_B) | Baustein 08 |
| ablenkung_linse | α(b) = (12Δs/b)·I(b/r*), I(x) = ∫₀^{π/2} cos θ·ν(cos²θ/x²) dθ (ν über a_N/a₀) | Baustein 08 |
| zusatzpotential | Φ_X = v_flach²·ln ν(s) = −v_flach²·ln(1 − e^{−s}), s = r*/r (vom Zentrum gezählt) | Baustein 08 |
| photonensphaere | r_ph = 9Δs = 3GM/c² | Baustein 08 |
| verwindung | a_J = J/(Mc) | Baustein 08 |
| umlauf_differenz | T_prograd − T_retrograd = 4π·a_J/c (unabhängig von r) | Baustein 08 |
| kerr_omega_phi | ω_φ = c³/(GM(r̃^{3/2} ± j)), r̃ = r/3Δs | Baustein 08 |
| kerr_omega_r | ω_r² = ω_φ²(1 − 6/r̃ + 8j·r̃^{−3/2} − 3j²/r̃²) | Baustein 08 |
| kerr_omega_theta | ω_θ² = ω_φ²(1 − 4j·r̃^{−3/2} + 3j²/r̃²) | Baustein 08 |
| knotendrehung | ζ = (ω_φ − ω_θ)/ω_φ ≈ 2j·(3Δs/r)^{3/2} (führende Ordnung) | Baustein 08 |
| nullflaeche_takt | r_E(ϑ) = 3Δs + √(9Δs² − a_J²cos²ϑ) | Baustein 08 |
| horizont_kerr | r_± = 3Δs ± √(9Δs² − a_J²); a_J² = r₊r₋, r₊ + r₋ = r₀ = 6Δs | Baustein 08 |
| mitfuehrung | \|Ω_Mitführung\| = (3Δs·a_J·c/r³)·\|3(ĵ·r̂)r̂ − ĵ\| = (3Δs·a_J·c/r³)·√(1 + 3cos²ψ) | Baustein 08 |
| mitfuehrung_polarbahn | Ω̄_Mitführung = ½·3Δs·a_J·c/r³ längs ĵ (Mittel über eine polare Kreisbahn) | Baustein 08 |
| leitstern_anteil | Ω_LS = Ω_LT·cos δ_L | Baustein 08 |
| geodaetisch | Ω_geodätisch = (9/2)·(Δs/r)·ω_Bahn, ω_Bahn = √(GM/r³) | Baustein 08 |
| quadrupol_leistung | P = (32/5)·(c⁵/G)·η_m²·(3Δs/r)⁵ | Baustein 08 |
| faltungsmesser_chirp | ḟ = (96/5)·π^{8/3}·(3Δs_c/c)^{5/3}·f^{11/3} ⇒ Δs_c = (c/3)·(ḟ/((96/5)π^{8/3}f^{11/3}))^{3/5} | Baustein 08 |
| ringdown_frequenz | f = c/(2π·3Δs)·F(j) | Baustein 08 |
| ringdown_abklingzeit | τ_RD = Q_RD/(πf_RD) = 6Δs·Q_RD/(c·F(j)) | Baustein 08 |
| richtungsfaktor_gedaechtnis | κ(ι) = sin²ι·(17 + cos²ι)/96 (Kante ι = π/2: 17/96) | Baustein 08 |
| gedaechtnis_welle | h_mem = 12·κ(Θ)·Δs_rad/D_L, Δs_rad = G·E_rad/(3c⁴) | Baustein 08 |

### Innere Richtung (I), 25 Formeln

| Kennung | Gleichung | Herkunft |
|---|---|---|
| tau_D | τ_D = c/a₀ | Baustein 07 |
| rueckkehr_exponent | s_p = √((p + Γ + θ)·τ_D) | Baustein 07 |
| rueckkehr | û = e^{−√((p + Γ + θ)·τ_D)} | Baustein 07 |
| schliesslaenge | Λ = c²/(2a₀) | Autor (FTM-PDF) |
| ell0 | ℓ₀ = 2Λ = c²/a₀ | Baustein 08 |
| a0_aus_schliesslaenge | a₀ = c²/(2Λ) | Autor (FTM-PDF) |
| kern_zeit | K(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/4T − ΓT} | Baustein 07 |
| kern_bahnfrequenz | \|K̂(iΩ)\|/K̂(0) = e^{−Re√((Γ + iΩ)τ_D) + √(Γτ_D)} | Baustein 07 |
| glaettung_innen | P_ϖ = e^{−\|k\|ϖ}, P_a·P_b = P_{a+b} | Baustein 07 |
| erneuerung | u(ϖ) = u(0)·e^{−ϖ√(Γ/D)}, u(0) = g_N + u(ϖ_R) = g_N/(1 − q) | Baustein 07 |
| eigenwert_schluss | λ_Λ = Dπ²/(4Λ²) | Baustein 07 |
| rueckkehr_geschlossen | q_Λ = cosh(s(Λ/ϖ_R − 1))/cosh(s·Λ/ϖ_R) | Baustein 07 |
| abweichung_schluss | ε_Λ = 1/tanh(kΛ) − 1 | Baustein 07 |
| temperatur_innen | kT_I = ħa₀/c = ħ/(Periode 2Λ/c) | Baustein 07 |
| kac_zeit_schluss | ⟨T⟩ = 2Λ/c = τ_D | Baustein 07 |
| watson_rueckkehr | G_d(1) = ∫₀^∞ e^{−t}·I₀(t/d)^d dt = 1/(1 − p_d) | Baustein 07 |
| rueckkehrzeiten_drehung | T_μ = {n ≥ 1: erste Rückkehr von x + nα_D (mod 1) in [0, μ)} | Baustein 07 |
| chsh_oberschwingungen | S = max_θ \|3E(θ) − E(3θ)\|, E(Δ) = −Σ_{m<N_H} w_m cos((2m+1)Δ), w_m ∝ (2m+1)⁻² | Baustein 04 |
| zustand_vakuum | 1 + w_s = −⅔·d ln a₀/d ln a | Baustein 03 |
| a0_aus_zelle | a₀ = 3c⁴/(2Gρ_mL³), ρ_m = ν_kosm·ρ_b | Baustein 03 |
| zelle_aus_a0 | L = (3c⁴/(2Gρ_m a₀))^{1/3} | Baustein 03 |
| entropie_vakuum | S_Λ = S_H(H_s) = (45/π³)(Λ/ℓ_P)² (in k_B) | Baustein 03 |
| gang_schwelle_strahlung | d ln a₀/d ln a = x/(1 + x), x = Ω_r(1 + z)/Ω_m (Strahlung verliert Faltung durch Rotverschiebung) | Baustein 03 |
| boden_drift | B = 2Pe·a₀ (Drift zum Torus hin) | Baustein 07 |
| boden_quer | B = √(2/π)·(ϖ_R/b)·a₀ (kontrahierende Querrichtung der Bandbreite b) | Baustein 07 |

### Zustandsraum (Z), 37 Formeln

| Kennung | Gleichung | Herkunft |
|---|---|---|
| bose_stufe_max | d(1 − e^{−s}) = s: Maximum von s^d·n_B; gibt (s_d, a_N/a₀, Wert) | Baustein 07 |
| realteil_quadrat_punkt | (Re τ_M)² = 1 + χ(s), Punktmasse | Baustein 07 |
| radius_resonanz | 1 + χ(s) = (p/k)², r_p/k = r*/s | Baustein 07 |
| im_tau_aus_nome | Im τ_Z = −ln\|q\|/2π, q = e^{2πiτ_Z} | Baustein 05 |
| kac_rueckkehrzeit | ⟨T⟩ = μ(X)/μ(A) | Baustein 06 |
| kac_budget | Σ_i 1/⟨T_i⟩ = Σ_i μ(A_i)/μ(X) = 1 für jede Zerlegung | Baustein 06 |
| bose_leiter | ℬ_d = s^d·n_B(s) | Baustein 07 |
| bose_leiter_integral | ℐ_d = ∫₀^∞ s^d·n_B ds = Γ(d+1)ζ(d+1) | Baustein 07 |
| fisher_metrik | g_F = −dν/ds = ν(ν − 1) | Baustein 07 |
| entropie_ruecklauf | S_n = ln ν + χ | Baustein 07 |
| schleifensumme | Σ_{n≥0} qⁿ = 1/(1 − q) | Baustein 05 |
| schleife_selbsterregt | q = (∂P/∂T)/G_th = −g_E·2/(3b_I)·(dR/dT)/G_th, P = g_E⟨I²⟩, ⟨I²⟩ = 2(g_E − R)/(3b_I) | Baustein 05 |
| kipp_schleife | q = 1 − 2(3a)^{1/4}δ^{1/2} vor der Falte | Baustein 06 |
| rastbreite | \|Δω\| ≤ 2k\|cos φ\| | Baustein 06 |
| preisach_ausgang | y_P = Σ μ_P·Zustand(α_P, β_P) nach der Folge der nicht ausgelöschten Extremwerte | Baustein 06 |
| rc_kette_nu | ν = 1/(1 − e^{−s}), s = n_R·arcosh(1 + RG/2) | Baustein 06 |
| eisenstein_G | G_k(τ) = Σ′(m + nτ)^{−k} = 2ζ(k) + 2(2πi)^k/(k−1)!·Σ σ_{k−1}(n)qⁿ | Baustein 05 |
| nu_aus_modul | ν = 1/(1 − \|e^{2πiτ}\|) | Baustein 05 |
| delta_reihe | Δ = q∏_{n≥1}(1 − qⁿ)²⁴ = Σ c_Δ(n)qⁿ | Baustein 05 |
| delta_modul | Δ(τ) = q∏(1 − qⁿ)²⁴, q = e^{2πiτ} | Baustein 05 |
| resolvente_flucht | ν_F = max_{1≤k≤K_M} 1/\|1 − q·e^{2πikα_D}\| | Baustein 05 |
| spiel_umkehrpunkt | v_K²/2 − 1 = −1/x_U + γ_K(x_U − 1)²/2 (H = v²/2 − 1/x + γ_K·y²/2, y = x − 1, t_K → ∞) | Baustein 05 |
| rate_rastpunkt | λ_ψ = −√((2k cos φ)² − Δω²) | Baustein 06 |
| kuramoto_einrasten | θ̇_i = ω_i + k_K·Σ_j sin(θ_j − θ_i), ω_i gleichabständig im Band W; (N − 1)k₁/W beim ersten Einrasten | Baustein 06 |
| rate_quer | λ_⊥ = −D_Z | Baustein 06 |
| formfaktor_zelle | f³ = ℓ_min³/V (kürzeste geschlossene Schleife hoch drei durch das Volumen der Zelle) | Baustein 03 |
| drehwinkel_zelle | φ_D = 360°/n mit e^{2πi/n}·(ℤ + τℤ) = ℤ + τℤ (Drehsymmetrien des flachen Torus ℂ/(ℤ + τℤ)) | Baustein 03 |
| masse_r_stern_gleich_compton | r* = λ̄ bei M_c = (ħ²a₀/(Gc²))^{1/3} = ε_L^{1/3}·m_Pl | Baustein 08 |
| page_zeit | t_P/t_F = 1 − (β/(1 + β))^{3/2} | Baustein 04 |
| rueckverdichtung_teile | N_RV = β·z_f²: (N/z_f)²S₀ ≥ N·β·S₀ | Baustein 04 |
| carter_faktor | 𝒞 = α_G(m_p)/(α¹²(m_e/m_p)⁴) | Baustein 02 |
| tiefe_faktor | f = θ_n = x₁/((n + 1)·m₁), kT_c = f·δ·m̄c² | Baustein 02 |
| decke_faktor | kT_max → (μ_m/μ_e)·m_ec²·max_x[x/4 − Π(x)] an der Grenzmasse | Baustein 02 |
| leiterzahl | ε_L = ℓ_P/ℓ₀ = a₀/a_P, a_P = c²/ℓ_P | Baustein 08 |
| schnittmasse_unten | M_min = ε_L·m_Pl = ħ/(c²τ_D): λ̄ = ℓ₀ und r* = ℓ_P | Baustein 08 |
| schnittmasse_oben | M_U = m_Pl/ε_L = c³τ_D/G: r* = 3Δs = ℓ₀ | Baustein 08 |
| einrast_phase | dΔθ/dt = Δω − K sin Δθ = 0 ⇒ Δθ = arcsin(Δω/K) für K ≥ \|Δω\| | Baustein 09 |

### Das Ganze (G), 27 Formeln

| Kennung | Gleichung | Herkunft |
|---|---|---|
| gewicht_w | w = 1 − e^{−√(θ̃/θ_c)}, θ_c = a₀/c | Baustein 03 |
| gewicht_w_preisach | w(θ̃), θ̃ = max(min θ der Geschichte, 0) | Baustein 03 |
| a0_modell | a_GH = cH/2π (Horizonttemperatur als Beschleunigung; heute nahe a₀) | Baustein 03 |
| weltzeit | t_A(a) = ∫₀^a da′/(a′H) = (1/H₀)∫_{−∞}^{ln a} d ln a′/E(a′) (Anfang bis a = 10⁻⁸ analytisch) | Baustein 04 |
| horizontwachstum | Ṡ_H = −2S_H·Ḣ/H, −Ḣ = H₀²(³⁄₂Ω_m a⁻³ + 2Ω_r a⁻⁴) (Verdünnung der Materie und Strahlung) | Baustein 03 |
| kreislauf_versiegen | ln a mit Ṡ_H(a) = Ṡ_i, d. h. −Ḣ = ℓ_P²H³Ṡ_i/(2πc²) (ab dort muss neue Faltung den Horizont tragen) | Baustein 04 |
| hintergrund_modell | H_b² = (8πG/3)(ρ_m(1+z)³ + ρ_r(1+z)⁴) + Λ_s c²/3, ρ = ω·3(100 km/s/Mpc)²/(8πG) | Baustein 04 |
| d_LSS_beobachter | d_LSS(z_B) = 2(c/H₀)∫_{z_B}^{1090} dz′/E(z′) (mitbewegt, heutige Einheiten) | Baustein 03 |
| teilchenhorizont | d_H = (c/H₀)∫_z^∞ dz′/E(z′) (mitbewegt) | Baustein 03 |
| rueckstellbilanz | ΔL/R = ΣΔs_b/R = 4πGρ_b(z)R²/(9c²), R = d_H/(1 + z) | Baustein 03 |
| faltung_horizont | Δs_H = r_H/6 = c/(6H) | Baustein 04 |
| theta_ | Θ = 2πa₀/(cH) | Baustein 03 |
| omega_vakuum | Ω_s = (π²/45)·Θ² | Baustein 03 |
| lambda_vakuum | Λ_s = π⁴/(15Λ²) | Baustein 03 |
| H_vakuum | H_s² = Λ_s·c²/3 | Baustein 03 |
| gewicht_hintergrund | 1 − w = e^{−√(6π/Θ)} (θ̃ = 3H) | Baustein 03 |
| zelle_hubble | L·H/c = (8π²/(ν_kosm·Ω_b·Θ))^{1/3} | Baustein 03 |
| zellgrenze | (1 + Q)·a₀ ≤ f³·3c⁴/(2Gρ_b L_min³), L_min = 0,985·d_LSS | Baustein 03 |
| a0_frueh_relativ | a₀(z)/a₀ = 1/(1 + Ω_r/(Ω_m a)), a = 1/(1 + z) | Baustein 03 |
| delta_N_eff | ΔN_eff = f_N/(1 − f_N)·(1 + 0,2271·N_eff)/0,2271 | Baustein 03 |
| horizontentropie | S_H = πc²/(ℓ_P²H²) (in k_B) | Baustein 03 |
| kreislauf_mindestrate | Λ̇/Λ ≥ Ṡ_i/(2S_Λ); gibt (Λ̇/Λ, S_Λ, H_s) | Baustein 03 |
| streuung_kosmisch | σ_ℓ = √(2/(2ℓ + 1)) (ein Himmel: 2ℓ + 1 Moden je Winkelordnung) | Baustein 03 |
| anteil_isotrop | Anteil = n/(2ℓ + 1) (n von 2ℓ + 1 Formen einer Winkelordnung im unendlichen, isotropen Raum) | Baustein 03 |
| wachstum_uebergabe | D″ + (2 + d ln E/d ln a)D′ = ³⁄₂Ω_m(a)(1 − ε)D, dε/d ln a = (1 − w)√(Ω_m(a)/2)(1 − ε)/N*, 1 − w = e^{−√(3H/θ_c(a))}; sofort: Quelle ³⁄₂Ω_m(a)(1 − (1 − w)Q/(1 + Q))D; gibt Δ(fσ₈)/fσ₈ = ΔD′/D′ bei z | Baustein 03 |
| neigung_stoerung | a_N = (Φ/c²)·c²·(1 + z)/d (Neigung einer linearen Störung der mitbewegten Skala d zur Epoche z) | Baustein 03 |
| q_fenster | Q ≥ 5 aus r*(M_b,Zelle)/Λ = √(6/ν_kosm) ≤ 1; Q ≤ (1 + Q)a₀-Grenze/a₀ − 1 | Baustein 03 |
<!-- ende -->

## 6 Verfahren

**Aufbau.** `kern/` enthält die Register: Symbole (jedes Zeichen einmal), Formeln (jede Funktion trägt genau einen Eintrag mit Gleichung, Ebene, Triadenrolle und Herkunft), Messwerte (Art G nur bei Werten aus dem Gedächtnis), Setzungen, Eingaben und Parametersätze. `pruefung/` enthält Prüfzahlen (Zahlen des Modelltexts), Gesetze und Gegenproben; `forschung/` Rechnungen, Themen und Beschlüsse; `waechter.py` die Wächter. Formeln und Konstanten stehen nur im Kern, kein Skript rechnet eigene.

**Ablauf.** `lauf.py` lädt alle Module, führt die Aufgaben in parallelen Prozessen aus (mit Zwischenspeicher), lässt pytest und die Wächter laufen, schreibt `bericht.md` und `register.json` und ruft am Ende `doku_erzeugen.py` auf. Bei einem strengen Fund stoppt der Lauf, ohne die Ergebnisse zu schreiben.

**Parametersätze.** Jede Prüfung rechnet in einem benannten Satz aus Lesart von a₀, H₀, Υ_d, Rotverschiebung, Zustand (gebunden oder mitatmend) und Beobachter:

<!-- erzeugt: parametersaetze -->
- **standard**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute
- **horizont**: a₀ horizont (1.04e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute (eigener Horizont)
- **z1090**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 1090.0, gebunden, Beobachter hier, bei der Rekombination
- **Upsilon_0.4**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.4, z = 0.0, gebunden, Beobachter hier, heute
- **Upsilon_0.6**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.6, z = 0.0, gebunden, Beobachter hier, heute
- **Upsilon_skaliert**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter hier, heute, Υ_b mitskaliert [W]
- **mitatmend_heute**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, mitatmend, Beobachter zwischen den Strukturen, heute
- **Sonnenkreis_Punktmasse**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter Sonnenkreis, Milchstraße als Punktmasse, Sonnenkreis punktmasse
- **Sonnenkreis_gemessen**: a₀ standard (1.1e-10 m/s²), H₀ = 67.4 km/s/Mpc, Υ_d = 0.5, z = 0.0, gebunden, Beobachter Sonnenort, gemessen, Sonnenkreis gemessen
<!-- ende -->

## 7 Prüfung

- **Prüfzahlen.** Jede Zahl im Modelltext hat eine Prüfzahl: Wert im Text gegen Rechnung aus dem Kern. Toleranz ist die halbe letzte gedruckte Stelle, absolut; eine eigene Toleranz braucht einen Grund. OFFEN heißt: die Formel fehlt im Kern oder die Zahl beruht nur auf einer Forschungsrechnung.
- **Gesetze.** Klasse A Grenzgesetze der Physik, B Galaxiengesetze der Literatur, C Zustandsraum-Gesetze, Z Tafeln des Autors, M eigene Gesetze aus 01–09, E Messwerte mit Quelle (Abstand in σ). Ergebnis OK, SPANNUNG, WIDERSPRUCH oder HINWEIS.
- **Gegenproben.** Eine zweite Rechnung auf anderem Weg oder ein Handwert prüft die Formel; eine Prüfung mit demselben Code zählt nicht.
- **Wächter.** Streng (ein Fund stoppt den Lauf): Eindeutigkeit der Symbole, Registrierung, Literale, Quellen, Modellsatz, Doppelformeln, Katalog, Themen, Beschlüsse, Daten u. a. Zählend (Lücken werden gezählt): Rückrichtung, Symbolkollision, Themenlücken, Abdeckung der Gegenproben, Uhr-Muster.
- **[G] und [W].** [G] hängt an einem Messwert aus dem Gedächtnis, nicht nachgeprüft; Naturkonstanten tragen kein [G]. [W] hängt an einer Wahl, die nicht im Modelltext steht und die der Autor bestätigt. Beide werden auf jede abgeleitete Zahl vererbt.
- **Regel für den Modelltext.** Zahlen kommen erst ins Modell, wenn sie unabhängig nachgerechnet sind; jede Änderung an `modell/` braucht einen Beschluss mit Datum.

<!-- erzeugt: pruefstand -->
Wächter im letzten Lauf: 20 Prüfungen, Funde je Wächter: symbole_eindeutig 0, registrierung 0, literale 0, quellen 0, tempsim 0, modellsatz 0, gegenproben 0, doppelformeln 0, katalog 0, katalog_luecken 0, themen 0, epochen_eingabe 0, zufall 0, beschluesse 0, daten 0, gegenproben_abdeckung 11, rueckrichtung 0, symbolkollision 89, themen_luecken 0, uhr_muster 3.

Messwerte 142 (davon mit [G] 115), Setzungen 9, Eingaben 180, Gegenproben 337.
<!-- ende -->

## 8 Ergebnisse des letzten Laufs

<!-- erzeugt: ergebnisse -->
Lauf f23e388, 2026-09-30T20:58:46.

**Prüfzahlen je Baustein**

| Baustein | Zahlen | OK | ABWEICHUNG | OFFEN | mit [G] | mit [W] |
|---|---|---|---|---|---|---|
| 01 | 4 | 3 | 0 | 1 | 0 | 0 |
| 02 | 76 | 71 | 0 | 5 | 27 | 0 |
| 03 | 274 | 219 | 0 | 55 | 75 | 0 |
| 04 | 20 | 20 | 0 | 0 | 8 | 0 |
| 05 | 13 | 13 | 0 | 0 | 1 | 0 |
| 06 | 24 | 20 | 0 | 4 | 6 | 1 |
| 07 | 112 | 86 | 0 | 26 | 13 | 0 |
| 08 | 45 | 45 | 0 | 0 | 15 | 7 |
| 09 | 21 | 20 | 0 | 1 | 3 | 4 |
| alle | 589 | 497 | 0 | 92 | 148 | 12 |

**Gesetze je Klasse**

| Klasse | Gesetze | OK | SPANNUNG | WIDERSPRUCH | HINWEIS |
|---|---|---|---|---|---|
| A | 24 | 24 | 0 | 0 | 0 |
| B | 12 | 9 | 2 | 0 | 1 |
| C | 16 | 16 | 0 | 0 | 0 |
| E | 9 | 8 | 0 | 1 | 0 |
| M | 124 | 124 | 0 | 0 | 0 |
| Z | 49 | 49 | 0 | 0 | 0 |

**Gesetze ohne OK**

| Gesetz | Aussage | Ergebnis |
|---|---|---|
| B_flaechendichte_beta | Zentrale Flächendichte, ideale Exponentialscheibe: Steigung β der Doppelpotenz (α = 1 − β) | SPANNUNG |
| B_flaechendichte_lelli_fit | Zentrale Flächendichte, ideale Exponentialscheibe gegen den Lelli-Fit bei Σ₀ = 10 / 100 / 1000 M☉/pc² [dex], innerhalb der Streuung | SPANNUNG |
| B_crater2_schaetzer | Crater II im Außenfeld: Plummer-Virial gegen den Schätzer McGaugh (σ² ∝ M statt M/2) [km/s] | HINWEIS |
| E_crater2 | Crater II σ_v, örtliche Regel (Jeans, Außenfeld der Milchstraße) | WIDERSPRUCH |

**Abweichungen der Prüfzahlen** (0)

| Prüfzahl | Baustein | Modelltext | Rechnung | Soll |
|---|---|---|---|---|
<!-- ende -->

## 9 Offene Punkte

<!-- erzeugt: offen -->
**Offene Beschlüsse** (Entscheidung der Hauptsitzung; der Modelltext bleibt bis dahin)

| Kennung | Frage | vorgeschlagener Satz |
|---|---|---|
| df2_df4_aussenfeld | DF2 und DF4 verlangten σ_v²/σ_v,N² ≤ 2,4 bzw. 1,5; erst ein 11- bzw. 206-fach stärkeres Außenfeld | [G]: Urteile der Zwerge offen, hängen an σ-Messwerten [G] |
| scherung_schleifen | Scherung der Schleifen kleiner als 10⁻⁶ der Atemrate | trägt nicht: einzige genannte Grundlage vakuum_richtungen trägt nicht |
| hauptreihe_faktor | Die Temperaturskala trifft die Hauptreihe auf den Faktor 0,7–1,3 | [G]: Vergleich mit Literaturwerten [G]; nach den Arbeitsregeln nur Forschungsstand |
| alpha_s_planck | 1/α_s ≈ 53 an der Planck-Skala, \|q̃\| = 1,2·10⁻²⁰, m_p +41 % / −30 % | [G]: Eingänge α_s(M_Z) und Quarkmassen [G]; Lauf von α_s nicht im Kern |
| umkehr_z_je_masse | Die örtliche Regel setzt zuerst bei z ≈ 20 (10⁸ M☉) ein | noch nicht im Modell: Forschungsstand: Umkehr-z je Masse noch nicht im Modell, Spektrum [G] |
| kollaps_frueher | Nach der Umkehr kollabiert ein Bereich in 0,3–0,96 der Zeit einer Staubkugel | noch nicht im Modell: Forschungsstand: kollaps_nach_umkehr noch nicht im Modell |
| linse_300_kpc | Linse bei N* = 2,4 um 3–7 % (200 kpc) und 19–39 % (300 kpc) unter der Relation | noch nicht im Modell: linsen_isoliert noch nicht im Modell; kumulativ andere Zahlen |
| kaleidoskop_fenster | Zelle oberhalb des 1,5-fachen Durchmessers nicht mehr prüfbar (27–42 Gpc) | ohne Eintrag im Forschungsstand: Zahl 1,5 und 42 Gpc ohne Rechnung im Forschungsstand |
| warp_gas_90_grad | Einfallendes Gas dreht die Knotenlinie nur bis 90°, bei der gemessenen Rate nach 0,14 Mrd. Jahren | [G]: Rate hängt an WARP_PRAEZESSION [G] |

**Themen**

| Thema | Frage | Stand |
|---|---|---|
| dunkle_energie | Was ist der zeitfeste Teil der Dunklen Energie neben dem Vakuumteil Ω_s (k_DE/k_s ≈ 2,9–3,1) und woher kommt seine Form w(z)? Hinweis des Autors: Tafel 5 und 6 („DM gegen B“) meinen Dunkle Energie; Invariante der Frage: Dunkle Energie als Bindungs- bzw. Schließungsstruktur statt Quelle; Trenntests je Kanal. | offen. Vakuumteil Ω_s im Modell. Vakuum je Richtung und je Zyklus widersprechen sich nicht: beide finden, dass kein Faktor 3 hergeleitet ist (die 3 steckt in π⁴/15); de_zyklen liest nur die Zahl k_DE/k_s ≈ 3 aus θ* [G], und SPARC schließt allein n = 1 aus. Offen, hängt an ω_b, ω_c [G]: die Mengenurteile der Fluchtmaß-Kette und modell_h0. |
| herkunft_q | Woher kommt Q = Ω_c/Ω_b ≈ 5,36 (ν_kosm = 1 + Q), das zeitfeste Gewicht der Mode ξ? | offen. Im Modell: Fenster 5 ≤ Q aus dem Schluss, Zählung je Radiant als offene Lesart. Offen, hängt an ω_b, ω_c [G]: q_form und lambda_bruecke (Kandidatenwahl). |
| a0_schwelle | Welchen Wert hat die Schwelle a₀ = c²/2Λ, in welcher Lesart, und was bindet sie an die Zelle? | im Modell: Standardwert 1,1·10⁻¹⁰, SPARC 1,06 ± 0,09 ± 0,04, Zelle ohne Kreispaare ≤ 1,13. Offen, hängt an N_eff [G]: die Urteile gegen die atmende Schwelle (die Zelle im frühen Kosmos steht im Modell, ihr Vergleich mit N_eff ist ein Gedächtniswert). |
| zwerge | Warum weichen druckgestützte Systeme im äußeren Feld (DF2, DF4, Pal 14, NGC 2419, Fornax, Crater II, Draco) von der Regel ab? | offen, hängt an den Messwerten der Zwerge [G]: alle Urteile „trägt nicht“ dieser Reihe. aussenfeld_kosmisch trägt nicht gegen SPARC (eigene Daten); seine frei angepasste Zahl steht im Modell (Beschluss offen). |
| warp | Was dreht die Knotenlinie des Warps der Milchstraße rechtläufig (+10,9 km/s/kpc [G])? | offen, hängt an WARP_PRAEZESSION [G]; im Modell steht die Kinematik am Knoten und die rückläufigen Kippmoden. |
| naht_uebergabe | Wie geht die Mode ξ an der w-Schale (Umkehrradius R_ta) in die örtliche Regel über? | Teil im Modell (I_F, Schale, Zuordnung der Ketten gesetzt); offen: ob die Übergabe λ = (1 − w)ω/N* die Mulde füllt und ob der Rest gebrochener Ketten ein Gedächtnis hat. |
| linse_haufen | Trägt die Regel mit Übergabe die Linsen um isolierte Galaxien und das Profil der Haufen? | offen: Spannung zwischen der aus dem Modell bestimmten Übergabe und den Repo-Linsen jenseits 0,5 Mpc; Haufenprofil. |
| mikroskala | Legt eine Gleichung des Modells m_p/m_e, α oder Im τ_X fest? | offen: keine Gleichung des Modells legt m_p/m_e oder α fest. |
| zellform | Welche Form hat die Zelle (Würfel, sechseckig, raum-, flächenzentriert)? | im Modell ohne festgelegte Form; hängt an Q [G] und dem Lauf des eigenen Atmens Ẋ_an(t). |
| zustandsraum_labor | Welche Aussagen des Zustandsraums sind prüfbar, und gibt es eine Laborvorhersage nur des Modells? | Form im Modell, keine eigene Laborvorhersage. |
| zeit_gedaechtnis | Welche Zeitskala hat das Gedächtnis der Regel (Kern K, λ_K, λ_m = ω_orb/N*), und wie hängt die Zeit an der inneren Richtung? | Kern im Modell; N* offen (0,3–2,4). |
| feldform_scheiben | Wirkt die Regel örtlich auf der Bahn (Kernform) oder als Feld (Feldform), und was zeigen Scheiben? | im Modell (örtliche Regel); offen die Vorhersagen in C. |
| rueckkehr_gleichung | Folgt die Regel ν = 1 + n_B(s) aus einer Gleichung mit wenigen Setzungen, und welche Invariante trägt sie? | im Modell; weitere Umformungen gesperrt. |
| trichter_hierarchie | Wie glättet die Hierarchie der Trichter (ρ > ρ_J) von Körper zu Galaxie zu Haufen? | Teil im Modell; Trichterkette offen, hängt an Massen [G]. |

**Beschlüsse**

48 Beschlüsse; ins Modell: 44, heraus: 4. Zur Bestätigung durch den Autor stehen die Wahlen [W] (10):

- DRAHT_Q_MERKLICH
- GPB_A
- ISS_HOEHE
- J_MILCHSTRASSE
- R_HAUFEN_HUELLE
- R_MILCHSTRASSE_HUELLE
- TORUS_DOMEGA
- TORUS_K
- T_KOERPER
- Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)
<!-- ende -->
