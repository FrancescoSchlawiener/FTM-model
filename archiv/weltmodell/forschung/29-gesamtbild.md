# Das Gesamtbild: alle Größen zusammen

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Die Notiz rechnet durchgehend mit a₀ = cH₀/2π = 1,0416·10⁻¹⁰ m/s² (im Text nicht umgerechnet); mit dem Standardwert a₀ = 1,2·10⁻¹⁰ sind r* Sonne 7 030 statt 7 547 AE, r* Milchstraße 8,35 statt 8,96 kpc, v_flach 175,8 statt 169,7 km/s, M_c 6,1 statt 5,78·10⁻²⁹ kg, r*_p/λ̄_p 145 statt 155,7. Räumliche Zelle L = c²/a₀ = 2πR_H (Kurzfassung 6 und 7, Variablenkarte, 3.2, 5, 6.1, 6.2, 6.3, G13, P2, P3) → trägt nicht: c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung; die Zelle folgt aus a₀ = 3c⁴/(2Gρ_mL³), L(1,2·10⁻¹⁰) = 26,8 Gpc statt 24,27 Gpc, ohne Kreispaare a₀ ≤ 1,13·10⁻¹⁰ statt 1,066·10⁻¹⁰; die Schranken H₀ ≤ 69,3 und Ω_m ≥ 0,297 entfallen. a₀ ∝ H(z) bzw. ∝ 1 + z (5.4, P1: +33,8 % bzw. +41,4 % bei z = 1) und „das Atmen setzt die Schwelle“ (6.1, 6.2) → trägt nicht (ΔN_eff = 2,09 bzw. Krümmung 0,22); a₀ ist seit z ≈ 10 fest. Proton-Kreuzung bei z = 1119 und M_c im Endzustand 30,4 MeV → setzen a₀ = cH(z)/2π voraus; bei festem Λ bleibt r*_p > λ̄_p, heute 145·λ̄_p; Zahlennähe markiert (Rekombination, Pionmasse, 2π-Faktoren in Ω_s). Hypothese 1/λ_m = v_flach/a₀ (4.8, P4, 6.3) → trägt nicht (ein Maxwell-Element gäbe 1 − q ∝ Γ ohne Wurzel); λ_m = ω_orb/N*, für eine bei z = 2 kollabierte Galaxie 1/λ_m = 0,1–1,2 Mrd. Jahre. q als Anteil „pro Umlauf“, Umläufe tragen die Bindung nach außen → Rückkehr als stationäre Schleife am Ort. Übergabe bei r_J ≈ 40 r* → an der Glättungslänge (Sonne 0,93 pc ≈ 27 r*). Stärke des Kerns κ₀ = GM/r³, κ̃₀/λ_m = q → κ₀/λ_eff = 3c²/(r*·v_rot). Im τ_M = r*/(2πr) als Festlegung → Wahl des Autors offen, nur Schreibweise. [G] an a₀ = 1,20 ± 0,02 ± 0,24, Ω_m und Ω_m h², Oort A und B, Masse und Abstand von M31, Q = 5,364 ± 0,065. Umbenannt (ohne Marke): χ (Spinzahl) → j; χ* (Abstand zur letzten Streufläche) → d_*; y = a_N/a₀ → s² (s Bose-Variable; y bleibt R/R⊕); Ω_Λ → Ω_DE.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/gesamtbild.py`

**Worum es geht.** Der Autor fragt, was sich ergibt, wenn man alle Variablen und das ganze Bild zugleich ansieht. Diese Notiz ordnet alle Größen in eine Karte, baut für sieben Körper die vollständige Längenleiter, klärt die Rolle der Zahl 2π, prüft sieben Querverbindungen mit Zufallsprüfung und nimmt die Deutung von a₀ als Grundfrequenz einer Torus-Schleife gesondert auf (sie trägt nicht: c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge) [korr.]. Am Ende stehen das Gesamtbild, der minimale Parametersatz und die neuen Beziehungen.

Zahlenwert der Schwelle: Diese Notiz rechnet durchgehend mit dem Modellwert a₀ = cH₀/2π bei H₀ = 67,36 km/s/Mpc, also a₀ = 1,0416·10⁻¹⁰ m/s², und mit R_H = c/H₀ = 4,4506 Gpc. Der gemessene Wert a₀ = 1,20·10⁻¹⁰ m/s² liegt um den Faktor 1,152 höher. Er steht nur dort, wo es ausdrücklich gesagt wird.

---

## Kurzfassung

1. **Zwei exakte Beziehungen tragen die ganze Längenleiter.** Für jeden Körper gilt
   - **Δs·λ̄ = ℓ_P²/3**, unabhängig von der Masse,
   - **r*² = 6π·Δs·R_H**.

   Mit μ = M/m_Pl und N_H = R_H/ℓ_P = 8,50·10⁶⁰ liegen dann alle Längen fest: λ̄ = ℓ_P/μ, Δs = ℓ_Pμ/3, r₀ = 2ℓ_Pμ und r* = ℓ_P√(2πμN_H). Masseunabhängig sind nur ℓ_P und R_H.
2. **Die Entropie der Falte ist die Kopplung.** Es gilt S/k = 4π·α_G(M) mit α_G(M) = 3Δs/λ̄ = (M/m_Pl)². Außerdem kT_H = Mc²/(8π·α_G(M)) und T_H/T_GH = R_H/(2r₀).
3. **Die Proton-Kreuzung ist eine Längengleichheit.** H/2π = α_G·ω_p ist dasselbe wie **r*_p = λ̄_p**: Der Schwellenradius des Protons schrumpft auf seine Compton-Länge. Heute gilt r*_p = 155,7·λ̄_p. Das setzt a₀(z) = cH(z)/2π voraus, das nicht trägt; bei festem Λ bleibt r*_p > λ̄_p immer, heute 145·λ̄_p (a₀ = 1,2·10⁻¹⁰). [korr.]
4. **Die 2π-Regel.** 2π rechnet zwischen „pro Radiant" und „pro Umlauf" um. Das 2π in a₀ ist dasselbe wie in der Gibbons-Hawking-Temperatur des Hubble-Horizonts: **a₀ = κ_dS/2π = c·kT_GH/ħ**, und a₀·(2π/H₀) = c. Das 2π des Nomes gehört zu einem anderen Kreis, dem azimutalen Umlauf des Bahntorus. Ω_s = π²/45 ist (2π)²/180: Aus der Summe über alle Rückläufe kommen vier Faktoren 2π, zwei davon nimmt a₀² wieder weg. [korr.] Zahlennähe ohne Herleitung, trägt nicht.
5. **Die Rückkopplungssumme ist eine Rückkehrstatistik.** ν = 1/(1−q) ist exakt die Kac-Rückkehrzeit in eine Zelle vom Maß 1−q. Q = q/(1−q) ist die mittlere Zahl der Rückläufe vor der Flucht. Kosmisch ist 1−q = Ω_b/Ω_m. Der Zahlenwert q = 0,843 selbst bleibt ohne Herleitung; die Zufallsprüfung findet so viele „schöne" Treffer, wie der Zufall erwarten lässt.
6. **Topologie-Zelle.** Mit a₀ = c²/L folgt L = 2πR_H = 27,96 Gpc, und L/d_LSS = π/I(Ω_m) = 1,008. H₀ kürzt sich heraus. Dass L so nahe am Durchmesser der letzten Streufläche liegt, ist eine Eigenschaft der Epoche und von Ω_m, keine Folge des Modells. Nähme man die Deutung ernst, ergäben sich a₀ ≤ 1,066·10⁻¹⁰ m/s², H₀ ≤ 69,3 km/s/Mpc und Ω_m ≥ 0,297; der Quadrupol wäre nach der diskreten Modensumme nur um etwa 15 % unterdrückt. [korr.] Die Deutung trägt nicht: c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge. Die Zelle folgt aus a₀ = 3c⁴/(2Gρ_mL³): Bei a₀ = 1,2·10⁻¹⁰ ist L = 26,8 Gpc. Ohne Kreispaare ist L ≥ 0,985·d_LSS (COMPACT 2023), also a₀ ≤ 1,13·10⁻¹⁰; Schranken auf H₀ und Ω_m folgen nicht. [korr.]
7. **Freie Parameter, ehrlich gezählt:** Auf Galaxienskala keiner. Im Kosmos mindestens vier: Q, die Form von w(θ_u), λ_m und der Rest der Dunklen Energie; die Zellgröße L kommt als fünfter hinzu, solange die Deutung L = 2πR_H nicht übernommen ist. Diese Deutung trägt nicht; gelesen als ganze innere Faltung bindet a₀ = 3c⁴/(2Gρ_mL³) die Zelle an die Schwelle. [korr.]

---

## 1. Variablenkarte

### 1.1 Alle Größen

Rolle im Zustandsraum: **par** = parametrisch (langsam, formt den Torus), **tan** = tangential (Umlauf, erhält Energie und Entropie), **trans** = transversal (Kontraktion, Dissipation, Entropieerzeugung).

**A. Faltung und Trichter**

| Symbol | Bedeutung | Formel / Definition | Einheit | Rolle | Skala |
|---|---|---|---|---|---|
| M, x | Masse, Masse in Erdmassen | x = M/M⊕ = Δs/Δs⊕ | kg, 1 | par | 9,1·10⁻³¹ kg (Elektron) … 1,5·10⁵³ kg (Universum) |
| y | Radius in Erdradien | y = R/R⊕ | 1 | par | |
| g | Hüll-Beschleunigung, Neigung der Kante | g = GM/r² = x/y²·g⊕ = 3c²Δs/r² | m/s² | par (Vorspannung) | g(R)/a₀ von 0,057 (Haufen) bis 2,6·10¹² (Sonne) |
| Δs | Faltungstiefe, Zusatzstrecke im Zentrum | Δs = GM/(3c²) | m | par | 2,3·10⁻⁵⁸ m (Elektron) … 1203 Mpc (Universum) |
| dτ/dt | Gangrate, Abwicklungstakt | √(1 − 6Δs/r) | 1 | par | |
| δ | Faltungs-Differential | 1 − dτ/dt ≈ 3Δs/r | 1 | par | |
| ε | Kurzform | ε = 3Δs/r = gr/c² | 1 | par | |
| r₀ | Grenzradius, Nullstelle des Takts, Engstelle | r₀ = 6Δs = 2GM/c² | m | trans (Engstelle, Entropie) | Sonne 2954 m |
| R | Hülle, Oberfläche | gemessen | m | par | |
| r* | Schwellenradius, g(r*) = a₀ | r* = √(GM/a₀) = √(r₀·πR_H) | m | Übergang par → tan | Sonne 7547 AE, Milchstraße 8,96 kpc |
| r_J | Gezeitengrenze im Elternfeld | (GM/\|T\|)^{1/3} | m | par (Elternfeld) | Sonne 1,372 pc |
| r_eq | Radius gleicher Neigung wie das Elternfeld | √(GM/g_Eltern) | m | par | |
| v_flach | flache Rotationsgeschwindigkeit | v⁴ = GMa₀ | m/s | tan | |
| Δs_eff | wirksame Faltung jenseits r* | Δs·ν(r) | m | par | |

**B. Schwelle, Atmen, Rückkopplung**

| Symbol | Bedeutung | Formel / Definition | Einheit | Rolle | Skala |
|---|---|---|---|---|---|
| H, H₀ | Atemrate des Torus | H = L̇/L | 1/s | par (global) | H₀ = 2,183·10⁻¹⁸ s⁻¹ |
| R_H | Hubble-Radius | c/H₀ | m | par | 4,4506 Gpc |
| a₀ | Schwelle | cH₀/2π = c²/(2πR_H) = κ_dS/2π | m/s² | par | gemessen 1,2·10⁻¹⁰; Modellwert 1,0416·10⁻¹⁰ (H₀ = 67,36) |
| s | Bose-Variable; s² = Feldstärke in Schwelleneinheiten | s² = a_N/a₀ | 1 | par | |
| τ_M | Modulparameter des Bahntorus | Re τ_M = ω_rad/ω_azim, Im τ_M = r*/(2πr) (Punkt 1 des Autors; das Modell schreibt \|q_Δs\| = Δs/r, Punkt 3; die Wahl ist offen, in der stationären Regel sind beide auf allen Bahnen gleich, χ(s) = s·q/(1 − q), sie betrifft nur die Schreibweise) [korr.] | 1 | Re: tan, Im: par | |
| q | Rückkehranteil am Ort, je Rücklauf der stationären Schleife (kein Anteil pro Bahnumlauf); als Nome des Bahntorus in Punkt 1 für die Punktmasse der Rückkehranteil, der Kern-Nome heißt q_Δs [korr.] | \|q\| = e^{−2π Im τ_M} = e^{−r*/r} = e^{−s} | 1 | par (Lage), Wirkung Bindung | q(r*) = 1/e für jeden Körper |
| ν | Rückkopplungssumme | Σqⁿ = 1/(1−q) = a/a_N | 1 | Bindung durch Rückläufe | ν(r*) = 1,582 für jeden Körper |
| a_eff | Zusatzfeld im Fernfeld | √(GMa₀)/r | m/s² | | |
| r*/3 | Kernlänge der Gedächtnisregel | r*/3 = GM/(3v_flach²) | m | par | |

**C. Gedächtnis, Operator, Schließung**

| Symbol | Bedeutung | Formel / Definition | Einheit | Rolle | Skala |
|---|---|---|---|---|---|
| 𝖪 | Gedächtniskern | κ₀·M_k(τ_M)·e_r⊗e_r·e^{−λ_mΔt} | s⁻² | par/trans | |
| κ₀ | Stärke des Kerns (Geschwindigkeitsform) | festgelegt ist nur κ₀/λ_eff = 3c²/(r*·v_rot), λ_eff = ω + λ_m²/ω [korr.] | s⁻² | par | |
| κ̃₀ | Stärke des Kerns (Beschleunigungsform) | κ̃₀/λ_m = q trägt nicht; gebunden ist κ₀/λ_eff = 3c²/(r*·v_rot) [korr.] | s⁻¹ | par | |
| λ_m | Zerfallsrate des Gedächtnisses | Nachlaufzeit 1/λ_m | s⁻¹ | trans (Relaxation) | Bedarf aus Haufenstößen: 1/λ_m ≥ 170–260 Mio a |
| M_k, Δ | Spitzenform, Nahfeldform | Δ = qΠ(1−qⁿ)²⁴ | 1 | | |
| E₂ | Zusammenhang der parametrischen Richtung | 1 − 24Σσ(n)qⁿ | 1 | par | |
| T_n | Hecke-Operator, Skalenwechsel | Mittel über Untergitter vom Index n; ℓ^{(p)}(s) = (1+p)ℓ(s) − pℓ(p·s) | – | Skalen | |
| c_eff | Stärke einer Ausschlussregel | 1 − 3/(2k+1) | 1 | trans | |
| 𝒫, Π⊥ | Schließungsprojektor, Fluchtmoden | 𝒫 = (1/M)ΣT_{mL}, Π⊥ = I − 𝒫 | – | tan (Identifikation) | |
| w(θ_u) | Projektionsgewicht der Fluchtrichtung | ż = (I − w·Π⊥)F | 1 | trans | w = 1 Hintergrund, w → 0 gebunden |
| θ_u | Expansionsrate des Zeitfelds | ∇_μu^μ | 1/s | par | |
| Q | Verstärkung der Schließungsregel (ν − 1), von einem Beobachter als Ω_c/Ω_b gelesen | Q = q/(1−q) = Ω_c/Ω_b | 1 | | 5,364 ± 0,065 [G] |
| Kac-Budget | Rückkehrzeiten eines geschlossenen Raums | Σ_i 1/⟨T_i⟩ = 1, ⟨T_A⟩ = μ(X)/μ(A) | 1 | tan (Schließung) | |
| κ_gap | Spektrallücke der Schließung | 4π²D/L² = D/R² (R = L/2π) | 1/s | tan (Topologie) | |
| q_ij, \|q₃\| | Schleifenverstärkungen im Dreieck Energie–Struktur–Zustand | Hopf-Schwelle \|q₃\|_krit = (Σt_i)(Σ1/t_i) − 1 ≥ 8 | 1 | alle drei | |
| W | Wicklungszahl | Frequenzverhältnis zweier Umläufe | 1 | tan | irrational: kein Einrasten |
| K_d, K_r | dissipative und reaktive Kopplung | K cos φ, K sin φ | 1/s | trans, tan | Einrasten bei \|Δω\| ≤ 2K_d |
| N_max | Zahl freier Zyklen | 1 + 2W_Band/(πk) | 1 | | |

**D. Teilchen, Kopplungen, Frequenzen**

| Symbol | Bedeutung | Formel / Definition | Einheit | Rolle | Skala |
|---|---|---|---|---|---|
| λ̄ | reduzierte Compton-Länge | ħ/(Mc) = ℓ_P²/(3Δs) | m | tan (Compton-Takt) | Proton 2,103·10⁻¹⁶ m |
| α_G(M) | Faltung in Compton-Einheiten | 3Δs/λ̄ = GM²/(ħc) = (M/m_Pl)² | 1 | par/tan | Proton 5,906·10⁻³⁹ |
| r_e | elektromagnetischer Radius des Elektrons | α·λ̄_e | m | tan | 2,818·10⁻¹⁵ m |
| α | Feinstruktur | r_e/λ̄_e | 1 | tan | 1/137,036 |
| a_J | Verwindung | J/(Mc) = (J/ħ)·λ̄ | m | tan (pro Radiant) | Elektron 1,931·10⁻¹³ m |
| j | Spinzahl | a_J/(3Δs) = (J/ħ)/α_G(M) | 1 | | Proton 8,47·10³⁷ |
| ω₁ | Compton-Takt Proton | c/λ̄_p | rad/s | tan | 1,4255·10²⁴ |
| ω₂ | Umlauf Elektron im H-Atom | α²m_ec²/ħ | rad/s | tan | 4,1341·10¹⁶ |
| ω₃ | akustische Schwingung | v_s/r | rad/s | tan | ≈ 10¹³ (Annahme) |
| Γ_G | Faltungsrate Proton | α_G·ω_p = Gm_p³c/ħ² | 1/s | par | 8,4191·10⁻¹⁵ |

**E. Thermodynamik der Falte und des Horizonts**

| Symbol | Bedeutung | Formel / Definition | Einheit | Rolle | Skala |
|---|---|---|---|---|---|
| S | Entropie der Falte | 36πk(Δs/ℓ_P)² = 4πk·α_G(M) | J/K | trans | Sonne 1,05·10⁷⁷ k |
| T_H | Hawking-Temperatur | kT_H = ħc/(24πΔs) = Mc²/(8πα_G(M)) | K | trans | |
| T_GH | Temperatur des Hubble-Horizonts | ħH₀/(2πk) | K | par | 2,654·10⁻³⁰ K |
| Ω_s, Λ_s | Vakuumwert des Strukturterms | Λ_s = 24ζ(4)(a₀/c²)², Ω_s = π²/45 | 1 | par | 0,2193 = 0,320·Ω_DE |
| L | Umfang der Torus-Zelle | a₀ = 3c⁴/(2Gρ_mL³); die Deutung L = c²/a₀ = 2πR_H trägt nicht (c²/a₀ = ℓ₀ = 2Λ ist keine Raumlänge) [korr.] | m | tan (Topologie) | 26,8 Gpc bei a₀ = 1,2·10⁻¹⁰ [korr.] |

---

## 2. Die Längenleiter eines Trichters

### 2.1 Alle natürlichen Längen

| Körper | λ̄ | Δs | r₀ = 6Δs | R | r* | r_J | R_H |
|---|---|---|---|---|---|---|---|
| Elektron | 3,862·10⁻¹³ m | 2,255·10⁻⁵⁸ m | 1,353·10⁻⁵⁷ m | r_e = 2,818·10⁻¹⁵ m | 7,640·10⁻¹⁶ m | – | 4,45 Gpc |
| Proton | 2,103·10⁻¹⁶ m | 4,140·10⁻⁵⁵ m | 2,484·10⁻⁵⁴ m | 8,409·10⁻¹⁶ m | 3,274·10⁻¹⁴ m | – | 4,45 Gpc |
| Erde | 5,890·10⁻⁶⁸ m | 1,478 mm | 8,87 mm | 6371 km | 1,956·10⁹ km | 1,496·10⁶ km | 4,45 Gpc |
| Sonne | 1,769·10⁻⁷³ m | 492,4 m | 2954 m | 6,957·10⁵ km | 7547 AE | 1,372 pc | 4,45 Gpc |
| Milchstraße | 2,949·10⁻⁸⁴ m | 197,4 AE | 1184 AE | 15,33 kpc | 8,96 kpc | 391 kpc | 4,45 Gpc |
| Haufen | 1,040·10⁻⁸⁷ m | 2,712 pc | 16,27 pc | 2 Mpc | 477 kpc | – | 4,45 Gpc |
| Universum | 2,345·10⁻⁹⁶ m | 1203 Mpc | 7220 Mpc | 14 260 Mpc | 10 050 Mpc | – | 4,45 Gpc |

Eingaben: Die Gezeitengrenze der Erde ist ihr Hill-Radius bei 1 AE. Die der Sonne folgt aus den Oort-Konstanten A = 15,3 und B = −11,9 km/s/kpc [G]. Die der Milchstraße ist der Hill-Radius gegen M31, mit M_bar(M31) = 1,5·10¹¹ M☉ [G] und D = 765 kpc [G]. Für Haufen und Universum gibt es kein Elternfeld mit Gezeitengrenze. Für Teilchen ist r_J bedeutungslos.

### 2.2 Universell oder massenabhängig

**Universell**, also für alle sieben Körper auf zehn Stellen gleich:

| Beziehung | Quotient links/rechts |
|---|---|
| Δs·λ̄ = ℓ_P²/3 | 1,0000000000 |
| r*² = r₀·πR_H = 6πΔs·R_H | 1,0000000000 |
| 3Δs/λ̄ = α_G(M) = (M/m_Pl)² | 1,0000000000 |
| GM/v_flach² = r* (Übergangsradius Newton/Operator) | 1,000000 |
| q(r*) = 1/e, ν(r*) = 1/(1 − 1/e) | 0,36788 bzw. 1,5820 |
| δ(r*) = (v_flach/c)² in erster Ordnung | 1,0000000000 |

Die letzte Zeile gilt exakt nur in erster Ordnung in Δs/r*. Mit dem vollen Takt δ = 1 − √(1 − 6Δs/r*) ist der Quotient für den Haufen 1,000009 und für das Universum 1,3068, denn dort ist r₀/r* = 0,72.

**Massenabhängig** (μ = M/m_Pl):

| Körper | μ | r*/λ̄ | R/r₀ | r*/R | r_J/r* | R_H/r* | g(R)/a₀ | q(R) |
|---|---|---|---|---|---|---|---|---|
| Elektron | 4,19·10⁻²³ | 1,98·10⁻³ | 2,1·10⁴² | 0,271 | – | 1,8·10⁴¹ | 0,074 | 0,76 |
| Proton | 7,69·10⁻²⁰ | 155,7 | 3,4·10³⁸ | 38,9 | – | 4,2·10³⁹ | 1516 | 1,2·10⁻¹⁷ |
| Erde | 2,74·10³² | 3,3·10⁷⁹ | 7,2·10⁸ | 3,1·10⁵ | 0,00077 | 7,0·10¹³ | 9,4·10¹⁰ | ≈ 0 |
| Sonne | 9,14·10³⁷ | 6,4·10⁸⁷ | 2,4·10⁵ | 1,6·10⁶ | 37,5 | 1,2·10¹¹ | 2,6·10¹² | ≈ 0 |
| Milchstraße | 5,48·10⁴⁸ | 9,4·10¹⁰³ | 2,7·10⁶ | 0,58 | 43,6 | 5,0·10⁵ | 0,34 | 0,56 |
| Haufen | 1,55·10⁵² | 1,4·10¹⁰⁹ | 1,2·10⁵ | 0,24 | – | 9,3·10³ | 0,057 | 0,79 |
| Universum | 6,89·10⁶⁰ | 1,3·10¹²² | 1,98 | 0,70 | – | 0,44 | 0,50 | 0,49 |

Zwei Befunde fallen auf:
- Beim Elektron liegt r* innerhalb der Compton-Länge (r*/λ̄ = 0,002) und sogar innerhalb von r_e. Beim Proton liegt r* heute 156 Compton-Längen außen.
- r_J/r* ist für Sonne und Milchstraße fast gleich, 37,5 und 43,6. Die Übergabe liegt aber nicht bei r_J: Ein kompakter Trichter bleibt bis zu seiner Glättungslänge ℓ Newtonsch, ℓ³ = GM/T mit T = √(3/2)·|E^TF| (voller spurfreier Gezeitentensor); für die Sonne ist ℓ ≈ 0,93 pc ≈ 27 r*, 1,37 pc ist der klassische Jacobi-Radius. [korr.]

### 2.3 Das Formelnetz

Es gibt vier massenabhängige Längen aus Konstanten: λ̄ ∝ M⁻¹, Δs ∝ M, r₀ ∝ M und r* ∝ M^{1/2}. Masseunabhängig sind ℓ_P und R_H. Masseunabhängige Produkte gibt es genau zwei, und sie legen alles fest:

```
                   ℓ_P² /3                       6π·R_H
      λ̄  ─────────────────────  Δs  ──────────────────────────  r*²
   (ħ/Mc)      Δs·λ̄ = ℓ_P²/3    (GM/3c²)    r*² = 6π·Δs·R_H     (GM/a₀)
      │                          │                               │
      │ 3Δs/λ̄ = α_G(M)           │ r₀ = 6Δs                       │ g(r*) = a₀ = c²/(2πR_H)
      │ S/k = 4π·α_G(M)          │ T_H/T_GH = R_H/(2r₀)           │ q(r*) = 1/e
      │ kT_H = Mc²/(8πα_G)       │ ω_H = kT_H/ħ = c/(4πr₀)        │ δ(r*) ≈ (v/c)²,  GM/v² = r*
      ▼                          ▼                               ▼
   Planck-Einheiten:   λ̄ = ℓ_P/μ,   Δs = ℓ_P·μ/3,   r₀ = 2ℓ_P·μ,   r* = ℓ_P·√(2π·μ·N_H),   R_H = N_H·ℓ_P
                       μ = M/m_Pl,  N_H = R_H/ℓ_P = 8,50·10⁶⁰
```

Die Probe r*/(ℓ_P√(2πμN_H)) = 1,0000000000 gilt für alle sieben Körper. Auf der Leiter log₁₀(Länge/ℓ_P):

| Körper | λ̄ | Δs | r* |
|---|---|---|---|
| Elektron | 22,38 | −22,86 | 19,68 |
| Proton | 19,11 | −19,59 | 21,31 |
| Erde | −32,44 | 31,96 | 47,08 |
| Sonne | −37,96 | 37,48 | 49,84 |
| Milchstraße | −48,74 | 48,26 | 55,23 |
| Haufen | −52,19 | 51,71 | 56,96 |
| Universum | −60,84 | 60,36 | 61,28 |

R_H liegt bei 60,93. λ̄ und Δs spiegeln sich am Punkt ℓ_P/√3, denn ihr geometrisches Mittel ist für jede Masse dasselbe. r* wächst nur mit halber Steigung, denn es ist das geometrische Mittel aus r₀ und πR_H. Bei Teilchen liegt Δs unter ℓ_P und λ̄ darüber, bei makroskopischen Körpern umgekehrt.

**Für rotierende Teilchen** mit J = ħ/2 folgt ebenso exakt a_J·Δs = ℓ_P²/6 und j = a_J/(3Δs) = (J/ħ)/α_G(M). Das Elektron hat j = 2,85·10⁴⁴, das Proton 8,47·10³⁷. Teilchen sind im Faltungsbild extrem überdrehte Trichter.

### 2.4 Besondere Massen

| Bedingung | Masse | Bemerkung |
|---|---|---|
| λ̄ = r₀ | μ = 1/√2, 1,539·10⁻⁸ kg | Planck-Skala |
| λ̄ = r* | M_c = (ħ²a₀/(Gc²))^{1/3} = 5,78·10⁻²⁹ kg = 32,4 MeV/c² = 63,5 m_e | im de-Sitter-Endzustand 30,4 MeV/c²; das setzt a₀ = cH/2π voraus, das nicht trägt; bei festem Λ bleibt M_c [korr.] |
| r* = R_H | r₀ = R_H/π, M = 2,94·10⁵² kg | |
| r₀ = R_H | Hubble-Masse M_H = c³/(2GH₀) = 9,25·10⁵² kg | dort gilt r* = √π·R_H (Probe 1,772454) |
| r* = r₀ | r₀ = πR_H, M = 2,90·10⁵³ kg | |

Das beobachtbare Universum hat r₀/R_H = 1,622 und liegt damit zwischen der Hubble-Masse und der Masse mit r* = r₀.

**Die Proton-Kreuzung als Längengleichheit.** [korr.] Diese Rechnung setzt a₀(z) = cH(z)/2π voraus; das trägt nicht: Atmete a₀ mit H, bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute. Bei festem Λ bleibt r*_p > λ̄_p immer, heute r*_p = 145·λ̄_p (a₀ = 1,2·10⁻¹⁰). Aus r*² = GM/a₀ und a₀ = cH/2π folgt:

```
r*_p = λ̄_p   ⟺   a₀(z) = G m_p³ c²/ħ²   ⟺   H(z)/2π = α_G·ω_p
```

Die Probe 2π·α_G·ω_p = 5,2899·10⁻¹⁴ s⁻¹ stimmt auf alle Stellen. Die Kreuzung wäre unter a₀ ∝ H der Moment, in dem der Schwellenradius des Protons auf seine eigene Compton-Länge schrumpft; bei festem Λ tritt sie nicht ein. [korr.] Teilchen mit M < M_c = 32,4 MeV/c² haben heute r* < λ̄. Das Elektron gehört dazu und bleibt es bei festem a₀. [korr.]

**Nahe Zahl:** M_c hat die Größenordnung der Pionmasse. Die bekannte Form ohne 2π, (ħ²H₀/(Gc))^{1/3} = 59,8 MeV/c², liegt einen Faktor 2,3 unter m_π. Das ist eine bekannte Koinzidenz großer Zahlen ohne Mechanismus im Modell. [korr.] Zahlennähe ohne Herleitung, trägt nicht.

### 2.5 Die Skalenleiter mit dem Lineal λ̄_p

Misst man jede Faltung in der Compton-Länge des Protons, so gilt 3Δs/λ̄_p = (M/m_p)·α_G. Die Exponenten zur Basis α_G:

| Stufe | 3Δs/λ̄_p | Exponent |
|---|---|---|
| Proton | 5,906·10⁻³⁹ | 1,0000 |
| Planck-Masse | 7,685·10⁻²⁰ | 0,5000 |
| Sternskala M_* = m_Pl³/m_p² (1,853 M☉, Δs_* = 912,2 m) | 1,301·10¹⁹ | −0,5000 |
| Sonne | 7,023·10¹⁸ | −0,4930 |
| Milchstraße | 4,213·10²⁹ | −0,7749 |
| Haufen | 1,194·10³³ | −0,8652 |
| Hubble-Masse | 3,265·10⁴¹ | −1,0859 |
| Universum | 5,297·10⁴¹ | −1,0914 |

Exakt in Stufen von α_G^{1/2} sind nur Proton, Planck-Masse und Sternskala; das ist die Beziehung Δs_*²·Δs_p = (λ̄_p/3)³ (Probe 1,0000000000). Galaxie, Haufen und Kosmos liegen **nicht** auf halbzahligen Stufen. Der Kosmos läge bei −1, wenn R_H/λ̄_p = 1/α_G wäre. Tatsächlich ist R_H/λ̄_p = 6,53·10⁴¹ gegen 1/α_G = 1,69·10³⁸, also Faktor 3857 daneben. Das ist Diracs Koinzidenz großer Zahlen. Sie hängt an der Epoche und ist kein Teil der Leiter.

---

## 3. Die Zahl 2π und die Schließung

### 3.1 Wo 2π auftritt

| Stelle | Form | welcher Kreis | Befund |
|---|---|---|---|
| Schwelle | a₀ = cH₀/2π = c²/(2πR_H) | Hubble-Schleife: Umfang 2πR_H, Periode 2π/H₀ | a₀·(2π/H₀) = c exakt. Die Schwelle baut in einer Periode die Geschwindigkeit c auf. |
| Horizont-Temperatur | kT_GH = ħH₀/2π, κ_dS = c²/R_H | euklidische Zeit des Hubble-Horizonts, Periode 2π/H₀ | a₀ = κ_dS/2π = c·kT_GH/ħ exakt (1,0416·10⁻¹⁰ beidseitig) |
| Grundmode einer Lichtschleife | f₁ = c/L, L = 2πR_H | dieselbe Hubble-Schleife | a₀ = c·f₁, der Wellenoperator liefert ω₁ = 2πc/L = H₀ |
| Nome | q = e^{2πiτ_M} | **azimutaler Umlauf des Bahntorus** | Phase e^{2πi Re τ_M} = radialer Phasenüberschuss pro Umlauf: echte Schließung. Betrag e^{−r*/r}: kein 2π in der Physik, das 2π in Im τ_M = r*/(2πr) ist Konvention. |
| Spektrallücke | 4π²D/L² = D(2π/L)² | Kreis der Länge L | zwei Faktoren, weil der Operator zweiter Ordnung ist; mit dem Radius R = L/2π lautet sie D/R² ohne 2π |
| Proton-Kreuzung | ħH/2π = α_G·m_p c² = G m_p²/λ̄_p | Horizonttemperatur gegen Gravitations-Selbstenergie | Epoche hängt an der Quantenlänge: ħ/(m_p c) → z = 1119, h/(m_p c) → z = 104 (Abschnitt 3.3); als Kreuzung r*_p = λ̄_p setzt das a₀ ∝ H voraus, das nicht trägt [korr.] |
| Strukturterm | Ω_s = π²/45 = (2π)²/180 | Summe über alle Rückläufe n und Hubble-Schleife | Λ_s = 24ζ(4)·(a₀/c²)², ζ(4) = (2π)⁴/1440. Vier Faktoren 2π aus der Rücklaufsumme, zwei davon hebt a₀² auf. [korr.] Zahlennähe ohne Herleitung, trägt nicht. |
| Hawking | kT_H = ħκ/(2πc) | euklidische Zeit der Falte | ω_H = c/(4πr₀) |
| Verwindung, Compton | a_J = J/(Mc), λ̄ = ħ/(Mc) | Drehung, Compton-Phase | Längen pro Radiant; pro Umlauf wären es 2πa_J und λ = 2πλ̄ |

Die Rechnungen im Einzelnen:
- 24ζ(4) = ∫₀^∞ 4t³·Σ_n e^{−nt} dt = 25,9757576091 = 4π⁴/15. Die Integrandensumme Σ_n e^{−nt} ist genau die Rückkopplungssumme ν − 1 = q/(1−q) mit t = s.
- Ω_s = 0,2193245422 = π²/45 = (2π)²/180.
- T_GH(H₀) = 2,654·10⁻³⁰ K, κ_dS = 6,544·10⁻¹⁰ m/s², κ_dS/2π = a₀.
- Die Unruh-Temperatur der Schwelle ist T_U(a₀) = T_GH/2π (0,159155).

### 3.2 Die gemeinsame Regel

> **Regel der Schließung.** 2π tritt genau dort auf, wo eine Winkelkoordinate θ ≡ θ + 2π geschlossen ist, und rechnet zwischen „pro Radiant" (R, ω, ħ, λ̄, a_J) und „pro Umlauf" (L, f, h, λ, 2πa_J) um. Eine Gleichung zwischen Raten oder Längen ist nur dann eine physikalische Aussage, wenn beide Seiten in derselben Zählung stehen.

Ist es jedes Mal dieselbe Schließung? **Nein, es ist dieselbe Operation an verschiedenen Kreisen.**
- **Ein und derselbe Kreis** sind die Hubble-Schleife (Umfang 2πR_H), die Periode des Atmens (2π/H₀) und die euklidische Zeitperiode des Hubble-Horizonts. Hier fallen drei Lesarten von a₀ exakt zusammen: Lichtumlauf, Atemfrequenz und Horizonttemperatur. Das ist die stärkste Stütze für A3. Sie ist allerdings eine Umformung und keine Herleitung: Alle drei Formen enthalten dieselbe Setzung „H₀ ist eine Kreisfrequenz". [korr.] Als Stütze trägt das nicht: c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine räumliche Schleife, und Λ atmet nicht mit H (sonst ΔN_eff = 2,09); a₀ = c²/2Λ ≈ cH₀/2π gilt nur heute.
- **Ein anderer Kreis** ist der azimutale Umlauf des Bahntorus, zu dem das 2π des Nomes gehört. Ein dritter ist die Compton-Phase des Teilchens.
- Ω_s verbindet zwei davon: die Rücklaufsumme der Bahntori (ζ(4)) und die Hubble-Schleife (a₀²).

### 3.3 Folge für die Proton-Kreuzung

[korr.] Die Gleichsetzung mit r*_p = λ̄_p setzt a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe dann zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), und r*_p > λ̄_p gilt immer. Die Gleichheit der Energien kT_GH = α_G·m_p c² bei z = 1119 bleibt als Rechnung; ihre Nähe zur Rekombination ist Zahlennähe.

Die drei Zählungen geben drei verschiedene Epochen:

| Zählung | Bedingung | z | kT |
|---|---|---|---|
| gemischt | H/2π = α_G·ω_p, gleichwertig r*_p = λ̄_p | 1119,2 | 0,2631 eV |
| einheitlich pro Radiant | H = α_G·ω_p, gleichwertig G m_p R_H/c² = λ̄_p² | 348,9 | 0,0822 eV |
| einheitlich pro Umlauf | H = α_G·ω_p/2π (a₀ pro Umlauf gegen λ = 2πλ̄) | 104,1 | 0,0247 eV |

Vergleich: Die letzte Streuung liegt bei z* = 1089,9 (0,2562 eV), die Gleichheit von Materie und Strahlung bei z_eq = 3448,7 (0,8102 eV). Nur die Form in der ersten Zeile trifft die Rekombination. [korr.] Zahlennähe ohne Herleitung, trägt nicht.

**Als Gleichheit zweier physikalischer Energien** hat diese Form eine eindeutige Lesart:

```
k·T_GH = ħH/2π  =  G·m_p²/λ̄_p  =  α_G·m_p c²
```

- Links steht die Temperatur, die ein Beobachter am Hubble-Horizont misst (Gibbons-Hawking). Das 2π ist dort die euklidische Periode 2π/H, also physikalisch und nicht Konvention.
- Rechts steht die Gravitationsenergie zweier Protonen im Abstand ħ/(m_p c). Das ist die Reichweite eines Feldes der Masse m_p (Yukawa-Länge), also die natürliche Quantenlänge des Teilchens.

Beide Seiten sind messbare Energien. Die Bedingung lautet damit: **Die Horizonttemperatur fällt auf die Gravitations-Selbstenergie des Protons auf seiner Yukawa-Reichweite.** Das geschieht bei z = 1119. [korr.] Zahlennähe ohne Herleitung, trägt nicht; eine Kreuzung r*_p = λ̄_p ist es nur unter a₀ ∝ H. Unberührt bleibt die Beziehung zwischen Konstanten m_p c²·α_G^{1/4} = 0,2601 eV gegen α³m_ec² = 0,1986 eV (Carter), denn sie enthält kein H.

---

## 4. Querverbindungen

Bewertung: **exakt** = folgt aus den Definitionen; **Umformung** = exakt, aber ohne neuen Inhalt; **Nähe** = Größenordnung oder Zahlennähe ohne Mechanismus; **nichts**.

### 4.1 Spektrallücke der Schließung gegen a₀

Die Lücke eines geschlossenen Kreises mit Umfang L ist κ_gap = D(2π/L)². Soll a₀ = c·κ_gap mit L = 2πR_H gelten, braucht es D = c·R_H/2π = 6,55·10³³ m²/s. Dieses D ist aus denselben Größen gebaut, die es erklären soll. Die natürliche Form braucht kein D: Der Wellenoperator −c²∂² auf dem Kreis L = 2πR_H hat die Grundfrequenz ω₁ = 2πc/L = H₀. Dann ist a₀ = c·ω₁/2π = c·f₁.

**Bewertung: Umformung.** a₀ als Grundmode einer Lichtschleife ist gleichwertig mit der Setzung L = 2πR_H. Diffusiv folgt a₀ nicht, weil D frei bliebe.

### 4.2 Kac-Budget gegen Rückkopplungssumme

Drei Aussagen sind exakt dieselbe Struktur:

```
ν = Σ qⁿ = 1/(1−q)           Rückkopplungssumme
u = 1/(1−F)                  Pólya: erwartete Besuchszahl bei Rückkehrwahrscheinlichkeit F
⟨T_A⟩ = 1/μ(A)               Kac: mittlere Rückkehrzeit in eine Zelle vom Maß μ(A)
```

Mit μ(A) = 1−q ist die Kac-Zeit genau ν. Das Budget mit zwei Zellen, (1−q) + q = 1, ist erfüllt. **q ist also ein Rückkehranteil je Rücklauf.** Die Rückläufe sind keine Bahnumläufe: Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter), kein Vorgang längs der Bahn. [korr.] Q = q/(1−q) ist die mittlere Zahl der Rückläufe, bevor der Zustand entweicht.

Kosmisch gilt Q = Ω_c/Ω_b = 5,364 ± 0,065 [G]. Daraus folgen q = Ω_c/Ω_m = 0,84287 ± 0,00159 und 1−q = Ω_b/Ω_m = 0,15713. Die Kac-Zeit in die baryonische Zelle ist ⟨T_b⟩ = 6,364, in die Strukturzelle ⟨T_S⟩ = 1,186.

In Galaxien ist q(r) = e^{−r*/r} < 1: Der Zustand ist transient und kehrt nicht sicher zurück. Positive Rekurrenz, also volle Schließung, entspricht q → 1, dem tiefen Fernfeld. Die offenen Gitter haben Rückkehrwahrscheinlichkeiten von 0,34 (sc), 0,28 (bcc) und 0,26 (fcc), keine davon nahe q.

**Bewertung: exakte Umformung.** Sie gibt q eine Bedeutung (Rückkehranteil je Rücklauf der stationären Schleife am Ort, 1−q = Leck je Rücklauf) [korr.], legt aber den Wert nicht fest.

### 4.3 Schleifenverstärkung im Dreieck gegen q = 0,843

Die Hopf-Schwelle des Dreiecks ist |q₃|_krit = 8,0000 bei gleichen Zeitskalen und 122,21 bei (1, 10, 100). Das kosmische q ist eine einzelne, statische Schleife mit q < 1 und der Reserve 1−q = 0,157. Bis zur oszillatorischen Schwelle fehlt mindestens der Faktor 9,5.

**Bewertung: nichts.** Es sind verschiedene Größen. Verträglich ist beides: Die kosmische Rückkopplung ist unterkritisch, sonst wäre ν unendlich.

### 4.4 Entropie der Falte, Informationserhalt, Gedächtnis

Exakt gelten:
- S/k = 4π·α_G(M): 7,42·10⁻³⁸ für das Proton, 1,05·10⁷⁷ für die Sonne, 5,97·10¹²² für das beobachtbare Universum.
- S = 1·k bei M = m_Pl/√(4π) = 6,14·10⁻⁹ kg.
- T_H/T_GH = R_H/(2r₀). Eine Falte ist genau dann heißer als der Hubble-Horizont, wenn r₀ < R_H/2. Für das beobachtbare Universum ist das Verhältnis 0,308.
- **Vereinigung ist transversal.** Δs ist additiv (parametrisch), S ist quadratisch in Δs. Beim Verschmelzen entsteht darum der Kreuzterm S(M₁+M₂) − S(M₁) − S(M₂) = 8πGM₁M₂/(ħc)·k, bei gleichen Massen 50,0 % der Endentropie (2·1 M☉: 2,10·10⁷⁷ k). Die Addition der Faltung ist also immer mit Entropieerzeugung verbunden.

Gedächtnis: Die tangentiale Dynamik erhält die Entropie (Liouville, Kac). Der Hecke-Übertrag zwischen den Skalen ist unitär beschränkt. Entropie entsteht nur transversal: beim Verdichten zur Falte und bei der Vereinigung.

**Bewertung:** Die drei Formeln sind exakte Umformungen der Hawking-Größen in Faltungsgrößen. Neu ist die Lesart S = 4πk·α_G(M): Die Entropie einer Falte ist ihre Gravitationskopplung. Eine quantitative Brücke zum Gedächtnis gibt es noch nicht.

### 4.5 Hecke-Leiter gegen Skalenleiter

Der Hecke-Operator T_n schickt q → qⁿ, also s → n·s. Das ist exakt r* → n·r* bei festem r, also **M → n²·M**. Die Stufen der physikalischen Leiter verlangen:

| Stufe | Massenfaktor | nötiger Index n |
|---|---|---|
| Proton → Planck | 1,30·10¹⁹ | 3,6·10⁹ |
| Planck → Sternskala | 1,69·10³⁸ | 1,3·10¹⁹ |
| Sonne → Milchstraße | 6,0·10¹⁰ | 2,4·10⁵ |
| Milchstraße → Haufen | 2,8·10³ | 53 |

**Bewertung: nur Analogie.** Die Übersetzung T_n ⟺ M → n²M ist exakt. Keine Stufe zeichnet ein kleines n aus, die Hecke-Struktur legt also keine physikalische Leiter fest.

### 4.6 Die eV-Ebene gegen die Frequenzleiter

Zwei Leitern von Konstanten:

| Stufe | E | ω = E/ħ [rad/s] |
|---|---|---|
| m_ec² | 5,110·10⁵ eV | 7,76·10²⁰ |
| α·m_ec² | 3729 eV | 5,67·10¹⁸ |
| α²·m_ec² (= ħω₂) | 27,21 eV | 4,134·10¹⁶ |
| **α³·m_ec²** | **0,1986 eV** | 3,02·10¹⁴ |
| α⁴·m_ec² | 1,449 meV | 2,20·10¹² |
| m_pc² (= ħω₁) | 9,383·10⁸ eV | 1,4255·10²⁴ |
| **m_pc²·α_G^{1/4}** | **0,2601 eV** | 3,95·10¹⁴ |
| m_pc²·α_G^{1/2} | 7,21·10⁻¹¹ eV | 1,10·10⁵ |
| m_pc²·α_G (= ħΓ_G) | 5,54·10⁻³⁰ eV | 8,42·10⁻¹⁵ |
| kT_rek (z*) | 0,2562 eV | 3,89·10¹⁴ |
| ħω₃ | 6,58 meV | 10¹³ |
| ħH₀ | 1,44·10⁻³³ eV | 2,18·10⁻¹⁸ |

Die eV-Ebene ist der Kreuzungspunkt der elektromagnetischen Leiter (Stufe α³) und der gravitativen Leiter (Stufe α_G^{1/4}). Das ist Carters Beziehung α_G ≈ α¹²(m_e/m_p)⁴, und sie hat einen Mechanismus (Sternhüllen). ω₃ gehört zu keiner der beiden Leitern. Es ist eine Materialgröße der Born-Oppenheimer-Stufen (α²m_ec²·(m_e/m_p)^{1/2} = 0,635 eV, ·(m_e/m_p) = 14,8 meV).

**Zufallsprüfung:** Paare (a, b) mit a ∈ [0, 24] und b ∈ [−12, 12], für die α^a·(m_e/m_p)^b innerhalb eines Faktors 3 an α_G liegt, gibt es 8: (1,11), (4,9), (7,7), (12,4), (15,2), (18,0), (21,−2), (24,−4). Die Exponenten 12 und 4 sind also nicht durch die Zahlen, sondern nur durch den Sternmechanismus ausgezeichnet.

**Bewertung:** Carter ist eine bekannte Nähe mit Mechanismus. Die Nähe von kT_rek hängt an der 2π-Zählung (Abschnitt 3.3). Zwischen ω₃ und dieser Ebene besteht keine Verbindung.

### 4.7 Zufallsprüfung für q = 0,8429

742 einfache Kandidaten wurden geprüft: Brüche a/b mit b ≤ 20, sowie Potenzen, Exponentialformen und 1 − 1/(mx) aus π, e, φ, √2, √3, ln 2, γ, ζ(3) und 2π. Im 1σ-Fenster [0,8413; 0,8445] liegen 4 Treffer, erwartet sind bei gleicher Dichte 4,2. Im 2σ-Fenster liegen 11. Die Treffer: erf(1) = 0,84270, 1 − 1/(11γ), exp(−ζ(3)/7) und 16/19. Ausgeschlossen sind e^{−1/6} (+2,3σ) und schwach 1 − 1/2π (−1,3σ).

**Bewertung: nichts.** q = 0,843 ist Ω_c/Ω_m und bleibt eine Eingabe.

### 4.8 Nachlaufzeit als Dynamikzeit an der Schwelle (Hypothese)

[korr.] Diese Hypothese setzt einen Maxwell-Kern voraus; das trägt nicht: Als Gedächtniskern ist die Rückkehr die Dichte erster Rückkehrzeiten (∝ Δt^{−3/2}); ein Maxwell-Element gäbe 1 − q ∝ Γ ohne Wurzel. Die Rückkehr ist eine stationäre Schleife am Ort, und das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m). Im Modell gilt λ_m = ω_orb/N* (N* ≈ 0,3–2,4); für eine bei z = 2 kollabierte Galaxie ist 1/λ_m = 0,1–1,2 Mrd. Jahre.

Der Toy-Kern ist ein Maxwell-Element, und ein Maxwell-Element relaxiert auf der Zeitskala seiner eigenen Dynamik. Nahe liegt deshalb 1/λ_m = r*/v_flach = v_flach/a₀ ∝ M^{1/4}:

| System | v_flach | r* | v/a₀ | 2π·v/a₀ |
|---|---|---|---|---|
| Zwerggalaxie 10⁸ M☉ | 34,3 km/s | 0,37 kpc | 10,4 Mio a | 0,07 Mrd a |
| Milchstraße 6·10¹⁰ M☉ | 169,7 km/s | 8,96 kpc | 51,6 Mio a | 0,32 Mrd a |
| Gruppe 10¹³ M☉ | 609,8 km/s | 115,7 kpc | 185,5 Mio a | 1,17 Mrd a |
| Haufen 1,7·10¹⁴ M☉ | 1238 km/s | 477 kpc | 376,7 Mio a | 2,37 Mrd a |
| Haufen 10¹⁵ M☉ | 1928 km/s | 1157 kpc | 586,7 Mio a | 3,69 Mrd a |

Unter dieser Annahme erfüllten Haufen den Bedarf von 1/λ_m ≥ 170–260 Mio. Jahren [korr.]. Der Faktor 2π (pro Radiant oder pro Umlauf) ist nach der Regel aus 3.2 offen.

**Bewertung: trägt nicht** (Maxwell-Kern); λ_m = ω_orb/N* bleibt mit N* als offener Zahl. [korr.]

---

## 5. a₀ als Grundfrequenz einer Torus-Schleife

[korr.] Dieser Abschnitt setzt L = c²/a₀ als räumliche Zelle voraus; das trägt nicht: c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge. Im Modell folgt die Zelle aus a₀ = 3c⁴/(2Gρ_mL³): Bei a₀ = 1,2·10⁻¹⁰ ist L = 26,8 Gpc. Ohne Kreispaare ist L ≥ 0,985·d_LSS (COMPACT 2023), also a₀ ≤ 1,13·10⁻¹⁰. Schranken auf H₀ und Ω_m folgen nicht.

### 5.1 Die Deutung

Die Deutung lautet: a₀ = c·f₁ mit der Grundfrequenz f₁ = c/L einer geschlossenen Lichtschleife vom Umfang L. Dann ist L = c²/a₀, und mit a₀ = cH₀/2π folgt L = 2πR_H = 27,96 Gpc. Das ist dieselbe Hubble-Schleife wie in 3.1, jetzt als räumliche Zelle gelesen.

Der Durchmesser der letzten Streufläche ist d_LSS = 2d_* = 2·I(Ω_m)·R_H mit

```
I(Ω_m) = d_*/R_H = ∫₀^{z*} dz / E(z),     E(z) = √(Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_DE)
```

Mit Ω_m = 0,3153 [G] ist I = 3,11592, d_* = 13,868 Gpc und d_LSS = 27,735 Gpc. Damit gilt

```
L / d_LSS = 2πR_H / (2·I·R_H) = π / I(Ω_m) = 1,00824
```

**H₀ kürzt sich heraus.** Das Verhältnis hängt nur von Ω_m und von Ω_r h² ab.

### 5.2 Frage 1: Folgt eine Zelle genau an der Grenze des Messbaren?

Numerisch ja: L = 1,008·d_LSS. Das folgt aber nicht aus dem Modell. Es folgt aus I(Ω_m) ≈ π, und das ist eine Eigenschaft von ΛCDM zur heutigen Epoche:

| Ω_m | 0,20 | 0,25 | 0,28 | 0,30 | 0,3153 | 0,33 | 0,35 | 0,40 | 1,0 |
|---|---|---|---|---|---|---|---|---|---|
| I = d_*/R_H | 3,727 | 3,416 | 3,266 | 3,178 | 3,116 | 3,060 | 2,989 | 2,832 | 1,937 |
| L/d_LSS | 0,843 | 0,920 | 0,962 | 0,989 | 1,008 | 1,027 | 1,051 | 1,109 | 1,622 |

Die Probe gegen Einstein-de-Sitter, 2(1 − 1/√(1+z*)) = 1,9394, stimmt (1,9366 mit Strahlung).

**Epoche.** Für einen Beobachter zu einer anderen Zeit, mit derselben Entkopplung:

| Skalenfaktor a_o | 0,25 | 0,5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| 2πR_H(t)/d_LSS(t) (Schleife folgt dem Atmen) | 1,662 | 1,492 | 1,008 | 0,504 | 0,239 |
| 2πR_H,0/d_LSS(t) (feste Zelle) | 1,899 | 1,336 | 1,008 | 0,857 | 0,793 |

In beiden Lesarten gilt L ≈ d_LSS nur heute. Es ist dieselbe Art Koinzidenz wie „warum jetzt" bei der Dunklen Energie.

### 5.3 Frage 2: Liegt d_* ≈ π·R_H etwas zugrunde?

- I = π gilt bei Ω_m = 0,3089. Planck liefert 0,3153 ± 0,0073 [G], der Abstand ist 0,89σ.
- Bei festem Ω_m h² = 0,1430 [G] (vom CMB festgelegt) und variablem H₀:

| H₀ [km/s/Mpc] | Ω_m | I | L/d_LSS |
|---|---|---|---|
| 64,00 | 0,3491 | 2,990 | 1,051 |
| 67,36 | 0,3152 | 3,117 | 1,008 |
| 70,00 | 0,2918 | 3,214 | 0,977 |
| 73,00 | 0,2683 | 3,324 | 0,945 |

- **Zufallsprüfung:** Ω_m ∈ [0,2; 0,4] überstreicht I ∈ [2,832; 3,727]. In diesem Intervall liegen fünf einfache Zahlen: 3, π, √10, e + ½ und 2φ.

**Bewertung: Parameter-Koinzidenz von ΛCDM.** Sie ist nicht unwahrscheinlich und hängt an Ω_m und an der Epoche. Das Modell erklärt sie nicht.

### 5.4 Frage 3: Was folgt, wenn man die Deutung ernst nimmt?

Die Deutung wird prüfbar, sobald man die Messgrenze für die Zelle durch unseren Ort ansetzt, L ≥ 0,985·d_LSS. Die Folgen gelten nur unter der Deutung L = c²/a₀, die nicht trägt; im Modell gibt die Zelle aus a₀ = 3c⁴/(2Gρ_mL³) die Grenze a₀ ≤ 1,13·10⁻¹⁰, und die Schranken auf H₀ und Ω_m entfallen. [korr.]

| Folge | Wert |
|---|---|
| Obergrenze der Schwelle | a₀ ≤ c²/(0,985·d_LSS) = 1,066·10⁻¹⁰ m/s² (nur unter der Deutung; im Modell a₀ ≤ 1,13·10⁻¹⁰ [korr.]). Der gemessene Wert 1,20 ± 0,02 (stat) ± 0,24 (sys) [G] liegt darüber, aber innerhalb seines systematischen Fehlers. |
| Hubble-Konstante bei festem Ω_m h² | H₀ ≤ 69,32 km/s/Mpc (nur unter der Deutung; entfällt [korr.]). |
| Materiedichte bei festem h | Ω_m ≥ 0,2974 (nur unter der Deutung; entfällt [korr.]). |
| Kreispaare | Keine, solange L > d_LSS, also für Ω_m > 0,3089. Darunter entstehen kleine Kreise (L = 0,985·d_LSS: α = 9,9°). |
| Richtung | Die Deutung legt keine Orientierung der Zelle fest. Eine Vorzugsrichtung ist keine Vorhersage. |

**Quadrupol**, gerechnet als diskrete Sachs-Wolfe-Modensumme im kubischen 3-Torus (P ∝ k⁻³, ohne ISW; Kontrolle bei L = 20·d_LSS: 0,9999):

| Fall | L [Gpc] | L/d_LSS | Kreisradius α | C₂-Rest | C₃-Rest | C₄-Rest |
|---|---|---|---|---|---|---|
| a₀ = 1,20·10⁻¹⁰ | 24,27 | 0,875 | 28,9° | 0,666 | 1,072 | 1,000 |
| a₀ = 1,10·10⁻¹⁰ | 26,48 | 0,955 | 17,3° | 0,795 | 1,042 | 0,981 |
| Messgrenze 0,985·d | 27,32 | 0,985 | 9,9° | 0,830 | 1,028 | 0,989 |
| **a₀ = cH₀/2π, L = 2πR_H** | **27,96** | **1,008** | – | **0,851** | 1,021 | 0,994 |
| a₀ = 1,00·10⁻¹⁰ | 29,13 | 1,050 | – | 0,879 | 1,014 | 0,997 |
| 1,2·d | 33,28 | 1,200 | – | 0,935 | 1,004 | 0,999 |
| 1,5·d | 41,60 | 1,500 | – | 0,976 | 1,001 | 1,000 |

Bei L = 2πR_H ist der Quadrupol um 15 % unterdrückt. Für a₀ von 1,0 bis 1,1·10⁻¹⁰ sind es 12 bis 20 %. Die Zuordnung der Zeilen zu a₀ setzt L = c²/a₀ voraus und trägt nicht; im Modell ist L(1,2·10⁻¹⁰) = 26,8 Gpc statt 24,27 Gpc. Die Werte je L bleiben: Für L = 27–42 Gpc ist der Quadrupol um 2–17 % gesenkt. [korr.] Die kosmische Varianz des Quadrupols beträgt 63 %. **Über das Leistungsspektrum allein ist die Zelle nicht nachweisbar.** Nachweisbar ist sie nur über die Korrelationen der a_ℓm zwischen verschiedenen ℓ und m.

**Unterscheidung vom Atmen.** [korr.] Beide Lesarten setzen eine wandernde Schwelle voraus; das trägt nicht: Atmete a₀ mit H (Λ mit dem Hubble-Radius), bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Mit den Schleifen atmend (Λ ∝ a, hier a₀ ∝ 1 + z) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 bis auf 0,3 % fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute. Die Tabelle gilt nur unter diesen Annahmen. Ist a₀ die Grundfrequenz einer festen Zelle, dann gilt a₀(z) = c²(1+z)/L. Folgt a₀ dagegen dem Atmen, dann gilt a₀(z) ∝ H(z). Beide Lesarten sind unterscheidbar; im tiefen Bereich gilt g_obs/g_bar ∝ √a₀:

| z | 0,5 | 1 | 2 | 3 |
|---|---|---|---|---|
| Atmen, a₀ ∝ H(z) | +15,0 % | +33,8 % | +74,2 % | +113,8 % |
| feste Zelle, a₀ ∝ 1+z | +22,5 % | +41,4 % | +73,2 % | +100,0 % |

Den größten Unterschied gibt es bei z ≈ 0,5–1 und bei z ≳ 3. Bei z = 2 fallen beide zusammen.

**Gesamtbewertung: Deutung.** Die Zahlennähe L ≈ d_LSS ist keine Folge des Modells. Die Folgen der Deutung, die drei Schranken (a₀, H₀, Ω_m) und die Form von a₀(z), gelten nur unter ihr; sie trägt nicht. Im Modell folgt aus der Zelle ohne Kreispaare a₀ ≤ 1,13·10⁻¹⁰, und a₀ ist seit z ≈ 10 fest. [korr.]

---

## 6. Das Gesamtbild

### 6.1 In einem Satz

> **Ein Trichter ist durch seine Masse und die Schwelle a₀ festgelegt: Die Faltung Δs = ℓ_P²/(3λ̄) setzt parametrisch die Vorspannung, am Ort der Bahn kehrt der Zustand mit dem Anteil q = e^{−r*/r} in sich zurück (eine stationäre Schleife am Ort, kein Vorgang längs der Bahn), und die Schließlänge Λ der Rückkehr, eine Länge der inneren Richtung (heute ≈ πc/H₀, fest), setzt über a₀ = c²/2Λ die Schwelle, an der der Trichter transversal ins Netzwerk übergeht. Zwei exakte Beziehungen, Δs·λ̄ = ℓ_P²/3 und r*² = 3Δs·ℓ₀ (heute ≈ 6π·Δs·R_H), spannen die ganze Leiter von der Compton-Länge bis zur Länge der Wirkung auf.** [korr.]

### 6.2 Diagramm

Das Diagramm zeigt einen makroskopischen Trichter von innen nach außen. Bei Teilchen kehrt sich der untere Teil um (Δs < ℓ_P < λ̄), und r* liegt nahe λ̄.

```
 log(Länge)
   │
 R_H ─┼──── Hubble-Radius; SCHLIESSLÄNGE Λ ≈ πR_H (innere Richtung, keine Raumlänge) [korr.] ──
   │      a₀ = c²/2Λ, heute ≈ κ_dS/2π           [parametrisch global: Λ setzt die Schwelle, heute fest] [korr.]
   │                          ▲
   │                          │  r*² = r₀·πR_H   (geometrisches Mittel aus Engstelle und Ganzem)
   │                          │
 ℓ   ─┼──── Übergabe ans Elternfeld: Glättungslänge (Sonne 0,93 pc ≈ 27 r*) [parametrisch: Elternfeld] [korr.]
   │                          │
 r*  ─┼──── Schwelle: g = a₀, q = 1/e, ν = 1,582, δ ≈ (v/c)²            [stationäre Rückkehr am Ort
   │        jenseits: q → 1, ν = 1/(1−q) wächst (Kac-Rückkehr)           trägt die Bindung] [korr.]
   │                          │
 R   ─┼──── Hülle: g = 3c²Δs/R²                                          [parametrisch: Vorspannung]
   │        Umläufe ω(s), Re τ_M = Phase pro Umlauf                     [tangential, frei, Newton]
   │                          │
 r₀  ─┼──── Engstelle 6Δs: dτ/dt = 0, S = 4πk·α_G(M), T_H                [transversal: Entropie,
 Δs  ─┼──── Faltungstiefe, additiv                                        Verdichtung, Rückgabe]
   │                          │  Δs·λ̄ = ℓ_P²/3
 ℓ_P ─┼──── Planck-Länge (Spiegelpunkt ℓ_P/√3 von Δs und λ̄)
   │                          │
 λ̄   ─┼──── Compton-Länge ħ/(Mc): für Körper winzig; für Teilchen der Takt ω = c/λ̄
   │         (Proton heute r* = 156 λ̄; die Kreuzung ist r* = λ̄)
   │
```

Die Rollen greifen so ineinander:
- **Parametrisch** wirken Δs und die Schließlänge Λ der Rückkehr (heute fest, a₀ = c²/2Λ). [korr.] Beide sind additive bzw. globale Größen und setzen die Landschaft: die Vorspannung g innen, die Schwelle a₀ außen.
- **Tangential** wirken die Umläufe. Innen (bis r*) laufen sie frei und Newtonsch, außen trägt die Rückkehr am Ort die Bindung (ν), eine stationäre Schleife, kein Vorgang längs der Bahn. [korr.] Die Phase Re τ_M ist ihr Umlaufanteil.
- **Transversal** wirken die Engstelle bei r₀ (Entropie, Hawking-Rückgabe), die Vereinigung von Falten (Kreuzterm der Entropie) und die Projektion der Fluchtrichtungen (w, Π⊥).
- Die Schließung erscheint doppelt: tangential als Identifikation (Kac, Spektrallücke) und global als Schließlänge Λ der Rückkehr, die a₀ setzt (heute ≈ πc/H₀). [korr.]

### 6.3 Minimaler Parametersatz

**Für die Längenleiter reichen vier Konstanten und eine Masse:** G, c, ħ, H₀ und M. Daraus folgen λ̄, Δs, r₀, r*, S, T_H, v_flach, δ(r*) und q(r) ohne weitere Eingabe. R, J und r_J brauchen den Körper bzw. seine Umgebung.

**Freie Parameter, ehrlich gezählt:**

| Bereich | frei | Bemerkung |
|---|---|---|
| Galaxien, stationär | 0 | a₀ = cH₀/2π, q = e^{−√(a_N/a₀)}. Offen sind nur die 15 % zwischen dem Modellwert cH₀/2π und dem gemessenen Wert 1,2·10⁻¹⁰ m/s². |
| Hierarchie E0 | 0 | Glättungslänge ℓ aus dem Elternfeld, ℓ³ = GM/T [korr.] |
| Kosmos: Gewicht | 1 | Q = 5,364 = Ω_c/Ω_b; es ersetzt Ω_c h² von ΛCDM |
| Kosmos: Projektionsgewicht | ≥ 1 | Form des Übergangs w(θ_u) |
| Gedächtnis | 1 | N* in λ_m = ω_orb/N*; die Hypothese 4.8 trägt nicht [korr.] |
| Dunkle Energie | ≥ 1 | 68 % bleiben nach Ω_s offen; für DESI zwei (w₀, w_a) |
| Topologie | 1 (+3) | L und Orientierung; die Deutung L = 2πR_H trägt nicht, gelesen als ganze innere Faltung bindet a₀ = 3c⁴/(2Gρ_mL³) die Zelle an die Schwelle [korr.] |
| Standard-Kosmologie | 4 | Ω_b h², A_s, n_s, τ_reio wie in ΛCDM |

Im Kosmos hat das Modell also heute nicht weniger Parameter als ΛCDM, sondern dieselbe Zahl in anderer Gestalt. Der Gewinn liegt auf Galaxienskala: Dort fallen alle Halo-Parameter weg.

---

## 7. Was sich neu ergibt

### 7.1 Neue, geprüfte Beziehungen

| # | Beziehung | Art | Probe |
|---|---|---|---|
| G1 | Δs·λ̄ = ℓ_P²/3 für jede Masse | exakt | 1,0000000000 (7 Körper) |
| G2 | Leiter in Planck-Einheiten: λ̄ = ℓ_P/μ, Δs = ℓ_Pμ/3, r* = ℓ_P√(2πμN_H) | exakt, aus G1 und r*² = 6πΔsR_H | 1,0000000000 |
| G3 | S/k = 4π·α_G(M), kT_H = Mc²/(8πα_G(M)), T_H/T_GH = R_H/(2r₀) | exakte Umformung | 1,0000000000 |
| G4 | Spin ½: a_J·Δs = ℓ_P²/6, j = (J/ħ)/α_G(M) | exakt | 1,000000 |
| G5 | Übergangsradius des Toy-Ansatzes GM/v_flach² = r* | exakt | 1,000000 |
| G6 | Proton-Kreuzung ⟺ r*_p = λ̄_p; heute r*_p = 155,7·λ̄_p; Grenzmasse M_c = 32,4 MeV/c² | exakt nur unter a₀(z) = cH(z)/2π, das nicht trägt; bei festem Λ bleibt r*_p > λ̄_p, heute 145·λ̄_p (a₀ = 1,2·10⁻¹⁰) [korr.] | 5,2899·10⁻¹⁴ beidseitig |
| G7 | a₀ = κ_dS/2π = c·kT_GH/ħ, a₀·(2π/H₀) = c, T_U(a₀) = T_GH/2π | exakte Umformung | 1,0416·10⁻¹⁰ beidseitig |
| G8 | Ω_s = (2π)²/180 aus 24ζ(4) = ∫4t³(ν−1)dt | exakt | 0,2193245422 |
| G9 | ν = Kac-Rückkehrzeit in die Zelle vom Maß 1−q; Q = Rückläufe vor der Flucht; kosmisch 1−q = Ω_b/Ω_m | exakte Umformung | Σ1/⟨T⟩ = 1,000000 |
| G10 | Vereinigung von Falten ist transversal: Kreuzterm 8πGM₁M₂k/(ħc) | exakt | 50,0 % bei gleichen Massen |
| G11 | Hecke T_n ⟺ M → n²M | exakte Übersetzung | – |
| G12 | Stufen α_G^{±1/2} nur für Proton, Planck-Masse, Sternskala; Galaxie −0,775, Haufen −0,865, Kosmos −1,09 | Befund | – |
| G13 | L/d_LSS = π/I(Ω_m), unabhängig von H₀ | exakt nur unter der Deutung L = 2πR_H, die nicht trägt (c²/a₀ = ℓ₀ = 2Λ ist keine Raumlänge); im Modell L = 26,8 Gpc bei a₀ = 1,2·10⁻¹⁰ [korr.] | 1,00824 |

### 7.2 Neue prüfbare Vorhersagen

| # | Vorhersage | Messung |
|---|---|---|
| P1 | Die Zeitabhängigkeit der Schwelle unterscheidet Atmen und Zelle: im tiefen Bereich +33,8 % (Atmen) gegen +41,4 % (Zelle) bei z = 1 und +15,0 % gegen +22,5 % bei z = 0,5. Beide Lesarten tragen nicht (Atmen mit H: ΔN_eff = 2,09; mit den Schleifen: Krümmung 0,22). Im Modell ist a₀ seit z ≈ 10 fest: Die Beschleunigungsrelation bei z = 0,5–1 ist die von heute. [korr.] | Rotationskurven und Linsen bei z = 0,5–1 (Euclid, JWST) |
| P2 | Mit der Zelle aus a₀ = 3c⁴/(2Gρ_mL³) und ohne Kreispaare (L ≥ 0,985·d_LSS, COMPACT 2023): a₀ ≤ 1,13·10⁻¹⁰ m/s². Die Schranken H₀ ≤ 69,3 km/s/Mpc und Ω_m ≥ 0,297 setzten die Raumlänge L = c²/a₀ voraus und entfallen. [korr.] | RAR-Kalibrierung, Suche nach Kreispaaren (Planck/LiteBIRD) [korr.] |
| P3 | Die Zelle unterdrückt den Quadrupol nur um etwa 15 % (12–20 % für a₀ = 1,0–1,1·10⁻¹⁰), das liegt unter der kosmischen Varianz. Die Zuordnung zu a₀ setzt L = c²/a₀ voraus und trägt nicht; für L = 27–42 Gpc sind es 2–17 %. [korr.] Nachweis nur über a_ℓm-Korrelationen. | Topologie-Likelihood mit voller Kovarianz |
| P4 | Hypothese 1/λ_m = v_flach/a₀ ∝ M^{1/4}: Haufen (1,7·10¹⁴ M☉) 377 Mio a, Gruppen 186 Mio a, Galaxien 52 Mio a. Die Hypothese trägt nicht (ein Maxwell-Kern gäbe 1 − q ∝ Γ ohne Wurzel); im Modell ist λ_m = ω_orb/N*, für eine bei z = 2 kollabierte Galaxie 1/λ_m = 0,1–1,2 Mrd. Jahre, und in Galaxien gibt es keinen messbaren Nachlauf. [korr.] | Haufen- und Gruppenstöße unterschiedlicher Masse |

Die Notation ist im Forschungsstand vereinheitlicht.
