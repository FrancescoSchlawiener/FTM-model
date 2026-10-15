# Orthogonalität und Wurzeln: Quadrat, Thaleskreis, Vieleck

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Milchstraße bei geschlossener Rückkehr „Λ/L = 4,6·10⁵–4·10⁷ bei 2,5–220 Gpc, Boden 2·10⁻⁶–2,5·10⁻⁸ a₀“ → mit dem Kern Λ = c²/(2a₀) = 12,1–14,0 Gpc, L = r*: Λ/L = 1,45–1,56·10⁶, Knick r₁₀ = 8,0–9,2 Gpc, Boden 6,4–6,9·10⁻⁷ a₀ (der alte obere Wert 4·10⁷ gehört zu einem Zwerg von 10⁸ M☉, ein Knick bei 220 Gpc verlangte Λ = 335 Gpc); Thales-Strecken p, q → p, t (q ist der Rückkehranteil); Hebel χ_g → f_g; Transportmatrix L → L_T; Schrittlänge ℓ → Δx; Vieleckseite s_n → k_n (s ist die Bose-Variable). Abschnitt 2 („Grenze“, Suchregel 5) setzt an einer Stelle ℓ₀ = 2πR_H, also a₀ = cH₀/2π (nur heute).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/orthogonalitaet_wurzeln.py` (Ausgabe `orthogonalitaet_wurzeln.txt`)
> Aufgenommen ist nur Hergeleitetes und Nachgerechnetes. s = √(a_N/a₀) = r*/r (Punktmasse), χ = s/(e^s − 1), ν = 1/(1 − e^{−s}), ℓ₀ = c²/a₀.

## 1  Senkrechte Anteile

Senkrecht heißt: Der Kreuzterm verschwindet. Dann addieren sich Quadrate, |v|² = Σvᵢ²; gleichgerichtete Anteile addieren sich linear; ein Produkt zweier Anteile ist eine Höhe (Abschnitt 2).

| Stelle | senkrechte Zerlegung | Hypotenuse |
|---|---|---|
| Normalform θ̇ = ω, ρ̇ = −κρ | ż = (−κ + iω)z: transversal reell, tangential imaginär | \|ż\|/\|z\| = √(κ² + ω²) |
| Kopplung K·e^{iφ} | K cos φ (rastet ein) und K sin φ (verschiebt nur Frequenzen) | K |
| Einrasten | Schwebung √(Δω² − 4K_d²); Amplitudenverhältnis r = tan(θ/2), sin θ = 2K_d/Δω | Δω; Einrasten bei θ = 90°, r = 1 |
| Nome q = e^{2πiτ} | ln q = −2π Im τ + 2πi Re τ: Zusammenziehen und Drehung je Rücklauf | 2π\|τ\|, von der Regel nicht genutzt |
| Pendelkanal | ω² = ω_el² ± ω_g², f_g = sin²θ mit tan θ = ω_g/ω_el | hängend √(ω_el² + ω_g²) |
| Kernform vertikal | ν_z² = ν_z,N² + (ν − 1)ω_N² (Steifigkeiten addieren) | ν_z |
| Faltungsgleichung | Y = h^{μν}∂σ∂σ = \|∇σ\|², h projiziert senkrecht zu u | \|∇σ\| = a_N/c² |
| Onsager-Casimir, GENERIC | L_T = L_s + L_a, tr(L_sᵀL_a) = 0, xᵀL_a x = 0 (keine Arbeit) | ‖L_T‖ |

**Epizyklen.** Im leeren Raum gilt κ² + ν_z² = 2Ω²: Radial- und Vertikalfrequenz sind Katheten, die Hypotenuse ist √2·Ω, die Diagonale des Quadrats mit Seite Ω. Allgemein ist κ² + ν_z² − 2Ω² = ∇²Φ = 4πGρ. Mit der Regel (Punktmasse, κ² = (1 + χ)Ω², ν_z² = Ω²):

| r/r* | 0,3 | 1 | 3 | 10 | 30 |
|---|---|---|---|---|---|
| (κ² + ν_z²)/Ω² = 2 + χ | 2,123 | 2,582 | 2,843 | 2,951 | 2,983 |

Der Fehlbetrag χΩ² ist genau ∇·(a − a_N). Wer den leeren Pythagoras erwartet, liest ihn als Dichte ρ_X = χΩ²/(4πG), außen v²/(4πGr²), den isothermen Halo. Das Dreieck öffnet die geänderte Regel, kein Stoff.

**Minkowski-Form.** Wo die Zeit eine Kathete ist, trägt sie ein Minus. Kreisbahn (Ω² = GM/r³ exakt): 1 = (dτ/dt)² + 6Δs/r + v²/c², also (dτ/dt)² = 1 − 9Δs/r; Stillstand auf der Bahn an der Photonensphäre 9Δs, ruhend bei r₀ = 6Δs. Schließungsregel (∂ξ)² = −1: (∂_tξ)² = 1 + |∇ξ|², die Einheitshyperbel γ² − (γv/c)² = 1; der Gradient ist eine Vierergeschwindigkeit, daher fällt die Mode wie Staub ohne Druck. Stehender Körper: ω² = ω_el² − ω_g² = 0 bei θ = 45°, f_g = ½, er knickt auf der Diagonale des Quadrats. Kerr √((3Δs)² − a_J²) und Schwebung √(Δω² − 4K_d²): Die Kathete verschwindet im Extremfall bzw. beim Einrasten.

**Kein Pythagoras.** (Re τ_M)² = 3 + γ_N(1 − χ/2) ist die Kettenregel (d ln a/d ln r = γ_N·d ln a/d ln a_N). Mit der Wirbelstärke ζ gilt κ² = 2Ω·ζ, also Re τ_M = √(2ζ/Ω): ein Produkt (Kepler ζ = Ω/2, Tiefbereich ζ = Ω). Die Feldform ν_z² = ν·ν_z,N² ist ein Faktor. g_tt·g_rr = −1 macht Takt und Faltung zu Kehrwerten (geometrisches Mittel 1, Ψ = Φ); ihre gleichen Hälften 6Δs/b der Ablenkung addieren sich linear, ebenso die Faltungen ΣΔs.

**τ als Rechteck und Quadrat.** Kepler: Re τ ≡ 0 (mod 1), das Gitter des Bahntorus ist rechtwinklig; die Regel schert es um √(1 + χ) − 1 ≤ 0,414. Mit Punkt 1 (Im τ = s/2π) ist er bei s = 2π (a_N = 4π²a₀, r = r*/2π) ein Quadrat τ ≈ i bis auf 0,006 Scherung, mit Punkt 3 (Im τ ≥ 2,5) nie. Die Regel nutzt nur die Kathete Im τ; mit der Phase fiele die Verstärkung (Re 1/(1 − q) = 0,869 bei r*, 0,509 bei 30 r*).

**Außenfeld.** Die Regel sieht |∇σ|: senkrecht s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e². Milchstraße bei 50 kpc (a_N = 0,028 a₀), Außenfeld 0,03 a₀: ν = 6,50 frei, 5,46 senkrecht, 4,68 gleichgerichtet.

**Thales-Kriterium.** F = PF + (I − P)F ist Mode für Mode senkrecht genau dann, wenn jeder Eigenwert p die Strecke [0, 1] unter 90° sieht, |p − ½| = ½. Π⊥ hat die Endpunkte 0 und 1. Das Maxwell-Element G*/γ = iωt_K/(1 + iωt_K) liegt ganz auf dem Kreis, weil der Dämpferanteil das Zeitintegral des Federanteils ist (Faktor 1/iω, 90°). Der Rückkehroperator 𝒫_q hat seine Eigenwerte auf dem Kreis um 1/(1 + q) mit Radius q/(1 + q), innen; der Pythagoras-Fehlbetrag ist bis zu 4q(1 − q)/(1 + q)², am größten (½) bei q = ⅓ (r = 0,91 r*), bei r* 0,497. Senkrecht zerlegen 𝒫_q und I − wΠ⊥ nur an den Enden (q bzw. w = 0 oder 1).

## 2  Wurzeln als Höhen im Thaleskreis

Im Thaleskreis über p + t (zwei Strecken; t hier keine Zeit) ist die Höhe h = √(p·t), das Verhältnis der Katheten √(p/t). AM ≥ GM, gleich nur bei p = t (gleichschenklig, 45°, h = Radius).

| Stelle | p | t | Höhe | p = t |
|---|---|---|---|---|
| Regel | a_N | a₀ | √(a_N a₀), das Tiefbereichsglied | r* (Kathetenverhältnis s) |
| Rückkehr | τ_Γ = c/a_N | τ_D = c/a₀ | 2·t_opt | r* |
| Schwelle | 3Δs | ℓ₀ | r* | nur Produkt fest |
| Planck | 3Δs | λ̄ | ℓ_P | Planck-Masse |
| Stern | ℓ_P | 3Δs_* | λ̄_p | – |
| Kerr | r_+ | r_−, r_+ + r_− = r₀ = 6Δs | a_J; Horizontfläche 4π(r_+² + a_J²) = 4π·r₀r_+ | j = 1 |
| Maxwell-Element | G' | γ − G' | G'' | ωt_K = 1, größter Verlust |
| Zeit an der Schwelle | r*/c | τ_D | t* = 1/Ω_Kepler(r*) | – |
| Vakuumteil | ℓ₀ | (45/4π⁴)·r* | r_Λ | – |

**Warum der Übergang bei r* liegt.** Die Rückkehr ist die Dichte erster Rückkehrzeiten ψ = L/√(4πDt³)·e^{−L²/4Dt}, gewichtet mit dem Verlust e^{−Γt}; das Integral ist exakt e^{−s}. Der Exponent ist das Minimum zweier Kosten:

```
s = min_t [ t/τ_Γ + τ_D/(4t) ] = √(τ_D/τ_Γ),     t_opt = ½·√(τ_Γ·τ_D),     am Minimum je s/2
```

Die Wurzel ist AM ≥ GM. Bei r* ist τ_Γ = τ_D, das Dreieck der Zeiten gleichschenklig, q = e^{−1}; innen verliert die Bahn schneller, als sie zurückkehrt (q → 0), außen langsamer (q → 1).

**Wurzelwachstum = Irrfahrt = Diagonale des N-Würfels.** Bei N ungeordneten Schritten mitteln die Kreuzterme weg, R² = N·Δx² (Schrittlänge Δx; gerechnet 1,004 ± 0,006 bei N = 100); gleichgerichtete wachsen wie N, das geschlossene Vieleck gibt 0. Im Modell ist D = L²/τ_D = GM/c = c·3Δs, also r*² = 3Δs·ℓ₀ = N·(3Δs)² mit N = ℓ₀/3Δs (Sonne 5,1·10²³, Milchstraße 8,5·10¹²), gleichwertig ℓ₀/r₀ = 2(r*/r₀)² Schritte der Länge r₀. r* ist die Diagonale eines N-Würfels der Kante 3Δs; ebenso ℓ_P² = (1/α_G)·(3Δs)².

**Kette der Skalen.** Planck, Teilchen und Stern liegen auf einer Leiter mit dem Schritt α_G^{−1/2} = 1,301·10¹⁹ (Logarithmus x = 44,01): 3Δs_p (−1), ℓ_P (0), λ̄_p (+1), 3Δs_* (+2); das ist Δs_*²·Δs_p = (λ̄_p/3)³. ℓ₀ liegt bei 3,23 x, nicht auf der Leiter: Die Schwelle bringt eine eigene Länge. Zwischen 3Δs (0) und ℓ₀ (1) teilen Schnitte von Potenzgesetzen den Logarithmus: r* bei 0,500 (a_N ∝ r⁻² gegen a₀), r_Λ bei 0,71–0,73 (Tiefbereich ∝ r⁻¹ gegen Vakuumteil ∝ r).

```
r*² = 3Δs·ℓ₀,     r_Λ² = (45/4π⁴)·ℓ₀·r*   (√(GMa₀)/r = Λ_s c² r/3),     r*²·λ̄ = ℓ_P²·ℓ₀   für jede Masse
```

Milchstraße: r_Λ = 4,84 Mpc (nur Vakuumteil; mit der ganzen Dunklen Energie √3 kleiner), t* = 46 Mio. Jahre; Sonne t* = 94 000 Jahre. f*_gw = 1/(πt*).

**Grenze.** Der Höhensatz legt nur das Produkt fest. Bei Längen verschiebt ein Faktor den Punkt p = t: 3Δs·ℓ₀ = r₀·πR_H (a₀ = cH/2π) gibt p = t bei 5,1·10²³ bzw. 1,3·10²³ M☉. Bedeutung hat p = t nur, wo p und t einzeln festliegen (a_N gegen a₀, τ_Γ gegen τ_D, r_+ gegen r_−, G' gegen γ − G'). Der Schalter ist kein Mittel: bei p = t ist χ = 0,582.

## 3  Kreis als Unendlich-Eck

- **Archimedes** verdoppelt mit dem Kathetensatz (Seite k_n): k_{2n}² = Durchmesser · Pfeilhöhe = 2(1 − √(1 − k_n²/4)). (π − n·k_n/2)·n²/π³ → 1/6, gerechnet bis n = 6,4·10⁹ mit 60 Stellen; die Form √(2 − √(4 − k_n²)) verliert in Gleitkomma ab n ≈ 10⁶ die Stellen.
- **Rückkehrsumme als Rechteck-Vieleck.** s·ν = s·Σ_{n≥0}e^{−ns} ist die Rechteckregel für ∫₀^∞e^{−x}dx = 1 mit der Schrittweite s = r*/r. Weit außen wird das Vieleck zum Kreis, ν → r/r* (Umfanggesetz). Einseitig gezählt (ein Umlaufsinn) ist der Fehler s/2 + s²/12, daher M_schein = r/r* + ½ + r*/(12r); beidseitig ist (s/2)·coth(s/2) − 1 ≈ s²/12, von zweiter Ordnung wie beim Vieleck (s = 0,01: (sν − 1)/s = 0,501, beidseitig /s² = 0,0833).
- **Halbgerade als Stufenleiter.** Mit dem Schritt r₀ und N = r*/r₀ Stufen fällt G je Stufe exakt um e^{−2 arsinh(h/2λ)}; s weicht um (s·r₀/r*)²/24 ab: Sonne 3·10⁻²⁵, Milchstraße 2·10⁻¹⁴, Haufen 6·10⁻¹¹, beobachtbares Universum 2·10⁻² (r*/r₀ = 1,3). Grob wird die Leiter erst bei der Masse des Hubble-Volumens.

## 4  Umfang zu Fläche eine Stufe höher und tiefer

**Höher: d = 1.** Schließt die Rückkehrrichtung bei Λ (G'(Λ) = 0), wird für kleine s 1 − q ≈ s²(Λ/L − ½), also a → a₀/(Λ/L − ½) = konstant: d_eff geht 3 → 2 → 1, die Flussfläche S_{d−1} von 4π über 2π auf 2. Bei Λ/L = 100 ist d_eff = 2,418 bei r*, 2,049 bei 10 r*, 1,45 bei 100 r*, 1,007 bei 1000 r*. 10 % Abweichung von Tully-Fisher ab 0,655·r*·Λ/L; Milchstraße mit Λ = c²/(2a₀) = 12,1–14,0 Gpc und L = r*: Λ/L = 1,45–1,56·10⁶, Knick r₁₀ = 8,0–9,2 Gpc, Boden 6,4–6,9·10⁻⁷ a₀ [korr.]. Jede geschlossene Richtung nimmt eine Flussdimension; bei d = 1 ist die Flussfläche ein Punktpaar, das π verschwindet. Gebundene Systeme erreichen diese Stufe nicht.

**Volumen: d_eff = 0.** Wo eine Quelle das Volumen füllt, gilt a ∝ r: innen in einer homogenen Kugel und beim Vakuumteil (a = +Λ_s c²r/3, nach außen). Das ist kein Flussgesetz; der Vakuumteil trifft den Tiefbereich bei r_Λ.

**Tiefer: kein Umfanggesetz.** Teilchen sind kompakte Trichter (Proton g(R)/a₀ = 1318), Newtonsch bis zu ihrer Glättungslänge; auf Teilchenebene gibt es im Modell kein Flussgesetz mit Umfang.

## 5  Suchregeln

1. **Kreuzterm prüfen.** Senkrechte Anteile addieren sich im Quadrat, gleichgerichtete linear. Beispiel: Steifigkeiten (Pendel, Kernform ν_z²), Außenfeld senkrecht s⁴ = s_i⁴ + s_e⁴, Onsager-Casimir. Gilt nicht für ΣΔs und 6Δs/b + 6Δs/b (gleichgerichtet) und die Feldform ν·ν_z,N² (Faktor).
2. **Leerer Raum ist der Pythagoras der Epizyklen.** κ² + ν_z² = 2Ω²; eine geänderte Regel öffnet das Dreieck, und der Fehlbetrag wird als Dichte gelesen (4πGρ_X = χΩ²). Gilt nicht in Materie: dort trägt 4πGρ_b mit.
3. **Die Zeit trägt das Minus.** 1 = (dτ/dt)² + 6Δs/r + v²/c², (∂ξ)² = −1, stehender Körper: Die Nullrichtung (45°) heißt Stillstand, Staub oder Knicken. Gilt nicht für räumliche Anteile (Y: plus).
4. **Tangential imaginär, transversal reell.** Rate −κ + iω, Kopplung K cos φ + iK sin φ, ln q = −2π Im τ + 2πi Re τ. Einrasten dreht die Rate; die Regel nutzt nur eine Kathete (|q|). Gilt nicht für die parametrische Rolle: Sie verschiebt die Ebene selbst.
5. **Wurzel = Höhe, Übergang bei p = t.** Jede Wurzel √(p·t) als Höhe über zwei Grenzen lesen; der Übergang liegt beim gleichschenkligen Dreieck, wo zwei Kosten gleich sind (s = min[t/τ_Γ + τ_D/4t], Kerr j = 1, Maxwell ωt_K = 1). Gilt nicht, wo nur das Produkt festliegt (3Δs·ℓ₀ = r₀·πR_H), und nicht für χ selbst.
6. **√N heißt ungeordnet und senkrecht.** √N gehört zu unkorrelierten Schritten (Diagonale des N-Würfels), N zu gleichgerichteten, 0 zum geschlossenen Vieleck. Beispiel r*² = N(3Δs)², ℓ_P² = (1/α_G)(3Δs)². Gilt nicht für Massen: M_* = N_*·m_p addiert linear.
7. **Vieleck gegen Kreis nach Schrittordnung.** Beidseitige Zählung gibt Fehler ∝ Schritt² (Archimedes π³/6n², Stufenleiter (r₀/r*)²/24, coth), einseitige ∝ Schritt (½ in M_schein). Jede geschlossene Richtung nimmt eine Flussdimension (3 → 2 → 1, π nur bis d = 2). Gilt nicht in Quellen (a ∝ r, d_eff = 0).
