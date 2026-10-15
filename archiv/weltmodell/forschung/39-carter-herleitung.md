# Carters Beziehung aus dem Sterntrichter: woher 12 und 4 kommen

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: optische Tiefe τ = Nσ_T/R² → τ = 3Nσ_T/(4πR²) (homogene Kugel; nur dazu passt τ_T = 3,9·10¹⁰, mit Nσ_T/R² wären es 1,63·10¹¹ und Faktor 340 statt 1400), in Faltungslängen (Δs/Δs_p)·(8π/3)r_e²/R² → (Δs/Δs_p)·2r_e²/R²; Wert 3,9·10¹⁰ und „Faktor 1400“ bleiben. [G] an M_div 1,2–1,4 M☉, T_eff 6200–7000 K, Opazitäten 1,2 und 20 cm²/g, H⁻-Bindung 0,754 eV, Hayashi 3000–4500 K, 0,075 M☉. Umbenannt: Faltungsdifferential ε = 3Δs/R → δ (ε_*, ε☉ → δ_*, δ☉); Lane-Emden ξ₁ → x₁; Ionisationsenergie χ = α²m_ec²/2 → E_ion; Saha-Logarithmus Λ_s = 1/(2f_s) → ℓ_S.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/carter_herleitung.py` (Ausgabe `carter_herleitung.txt`)

## Frage

Warum gilt

```
3Δs_p/λ̄_p  ≈  (r_e/λ̄_e)¹² · (λ̄_p/λ̄_e)⁴          ⟺          α_G ≈ α¹²·(m_e/m_p)⁴ ,     𝒞 := α_G/(α¹²(m_e/m_p)⁴) = 2,944 ?
```

Gesucht ist eine physikalische Gleichgewichtsbedingung, aus der die Exponenten 12 und 4 folgen. Sie soll in der Sprache des Modells formuliert sein: ein Stern als Trichter, in dem die Faltung (parametrisch, α_G) den Energiefluss antreibt und elektromagnetische Raten ihn abführen (transversal: Strahlungsdiffusion; tangential: atomare und nukleare Schwellen, Konvektionsumläufe). Danach wird geprüft, ob das Einrasten gekoppelter Zyklen, dΔθ/dt = Δω − K sin Δθ, eine gleichwertige Bedingung liefert.

Die Beziehung hat vier Längen: die Faltungstiefe des Protons Δs_p = Gm_p/(3c²), den EM-Radius des Elektrons r_e = αλ̄_e und die beiden Compton-Längen λ̄_p, λ̄_e. In vierter Wurzel lautet sie (3Δs_p/λ̄_p)^{1/4} = 𝒞^{1/4}·r_e³λ̄_p/λ̄_e⁴ mit 𝒞^{1/4} = 1,310.

## 1. Die Größen des Trichters

| Größe im Stern | Formel | Länge im Modell | Rolle |
|---|---|---|---|
| Masse | M | Faltungstiefe Δs = GM/(3c²) | parametrisch |
| Teilchenzahl | N = M/m_p | Δs/Δs_p | parametrisch |
| Sternskala | M_* = m_Pl³/m_p² = 1,853 M☉ | 3Δs_*/λ̄_p = α_G^{−1/2}, Δs_* = 912 m, N_* = α_G^{−3/2} | parametrisch |
| Hydrostatik | kT_c = θ·GMm_p/R | kT_c = θ·m_pc²·δ, **δ = 3Δs/R** [korr.] | parametrisch |
| Zündung p+p | kT_c = f_nuc·α²m_pc², E_G = π²α²m_pc² = 493 keV | E_G = π²ħc·r_e²/(λ̄_e²λ̄_p) | tangential (Coulomb-Barriere) |
| Thomson-Streuung | σ_T = (8π/3)r_e² | EM-Radius | transversal (Photonendiffusion) |
| Eddington-Leuchtkraft | L_Edd = 4πcGMm_p/σ_T | **L_Edd = (9/2)·m_pc³·Δs/r_e²** | Faltung durch EM-Fläche |
| Strahlungsdruckanteil | 1 − β ≈ (μ_M μ_mol²/C_E)², μ_M = M/M_* | (Δs/Δs_*)² | parametrisch |
| Ionisationsschwelle | kT = f_s·α²m_ec² | α²m_ec² = ħc·r_e²/λ̄_e³, a_B = λ̄_e²/r_e | tangential (Elektronumlauf) |
| Konvektion | ∇_rad > ∇_ad | neuer Umlauf (Zirkulation) | tangential, entsteht |

Die Buchstaben f_nuc, f_s, θ, μ_mol, C_E und κ̄/κ_T sind reine Zahlen. Sie werden in Abschnitt 4 mit Sonnenzahlen belegt.

## 2. Die Kette Schritt für Schritt

Einheiten ħ = c = k = 1, also G = m_Pl^{−2} und α_G = (m_p/m_Pl)².

**2.1 Hydrostatik (parametrisch).** Das Virial liefert kT_c = θ·GMm_p/R. Im Modell ist GMm_p/R = m_pc²·3Δs/R, also

```
kT_c = θ · m_p c² · δ ,        δ = 3Δs/R   (Faltungsdifferential am Trichterrand)   [korr.]
```

Die Zentraltemperatur ist das Faltungsdifferential am Rand, gemessen in der Ruheenergie des Protons. Bei fester Teilchenzahl und Dichte gilt kT_c = θ·α_G·N^{2/3}·ħc·n^{1/3}, also T_c ∝ α_G. Auf der Sternskala N = μ_M α_G^{−3/2} wird daraus kT_c = θμ_M·α_G^{−1/2}·ħc/R: Die Faltung legt nur fest, wie groß der Trichter bei gegebener Temperatur ist.

**2.2 Zündung (tangential, Coulomb-Barriere).** Die Temperatur stellt sich dort ein, wo Protonen die Coulomb-Barriere durchtunneln: kT_c = E_G/(4(Λ_G/3)³) mit dem Gamow-Logarithmus Λ_G = 13,5 (Sonne), also f_nuc = 0,0271. Zusammen mit 2.1:

```
δ_* = 3Δs_*/R_* = (f_nuc/θ) · (r_e/λ̄_e)²          (Sonne: 2,1226·10⁻⁶ = 0,040·α²,  Rotverschiebung cδ = 636 m/s)   [korr.]
```

Das Faltungsdifferential an der Oberfläche eines brennenden Sterns ist durch die EM-Kopplung festgelegt, nicht durch α_G. Der Radius folgt als R_* = θμ_M α_G^{−1/2} λ̄_p/(f_nuc α²).

**2.3 Transport (transversal, Thomson).** Die Strahlung diffundiert durch Streuung an Elektronen. Exakt gilt L_r = (4πcGM_r/κ)·dP_rad/dP. Im Eddington-Standardmodell (Polytrope n = 3) ist dP_rad/dP = 1 − β konstant, und

```
L = λ · L_Edd ,     λ = (κ_T/κ̄)(1 − β) ,     L_Edd = (9/2) m_p c³ Δs/r_e² ,
M/M_* = C_E · μ_mol⁻² · √(1−β)/β² ,     C_E = 4·x₁²|θ'(x₁)|·√45/π^{3/2} = 9,7255   [korr.]
```

Die Lane-Emden-Lösung gibt x₁ = 6,89685, x₁²|θ'| = 2,01824 [korr.], ρ_c/ρ̄ = 54,18. Für kleine 1 − β gilt 1 − β = (μ_M μ_mol²/C_E)² ∝ (Δs/Δs_*)². Die Sternskala M_* ist damit die Faltung, bei der der Strahlungsdruck (reine Photonenzahl) gegen den Gasdruck (Protonenzahl) ins Spiel kommt; Gleichstand β = 1/2 erst bei 72 M_* = 133 M☉.

**2.4 Photosphäre.** Mit L = 4πR²σ_SB T_eff⁴, σ_SB = π²/60:

```
T_eff⁴ = (45/2π³) · [μ_M μ_mol⁴ f_nuc² / (C_E² θ² (κ̄/κ_T))] · α² m_e² m_p³ / m_Pl
```

Die Masse geht als T_eff⁴ ∝ μ_M ein.

**2.5 Grenzbedingung.** Eine Sternhülle bleibt radiativ, wenn die Leuchtkraft, die der Transport liefert, nur von einer Photosphäre oberhalb der Wasserstoffschwelle abgestrahlt werden kann. Liegt T_eff darunter, liegt die Ionisationszone von Wasserstoff direkt unter der Oberfläche: Die Opazität steigt dort steil an (gebunden-frei, H⁻), und der adiabatische Gradient sinkt von 0,4 auf 0,1 (Abschnitt 6). Dann gilt ∇_rad > ∇_ad, und die Hülle konvektiert. Die Grenze ist

```
kT_eff = f_s · α² m_e c²
```

Einsetzen und Quadrieren (α_G tritt als m_p/m_Pl = α_G^{1/2} auf):

```
α_G  =  (2π³/45)² · C_E⁴ (κ̄/κ_T)² θ⁴ f_s⁸ / (μ_M² μ_mol⁸ f_nuc⁴)  ·  α¹² (m_e/m_p)⁴
```

Die Monom-Rechnung des Skripts gibt diese Zeile exakt, mit dem Zahlenfaktor (2π³/45)² = 1,899.

**2.6 Herkunft der Exponenten.** Beitrag jedes Kettenglieds zum Exponenten in α_G:

| Glied | α | m_e | m_p | μ_M |
|---|---:|---:|---:|---:|
| Photosphäre kT_eff = f_s α² m_e (T_eff⁴ ∝ α⁸ m_e⁴) | +16 | +8 | 0 | 0 |
| Zündung kT_c = f_nuc α² m_p (über R⁻², T_c²) | −8 | 0 | −4 | 0 |
| Thomson σ_T ∝ α²/m_e² | +4 | −4 | 0 | 0 |
| Gravitation (GMm_p)⁻¹ = (μ_M m_Pl/m_p)⁻¹ | 0 | 0 | −2 | +2 |
| Strahlungsdruck 1 − β ∝ μ_M² | 0 | 0 | 0 | −4 |
| α_G = m_p²/m_Pl² | 0 | 0 | +2 | 0 |
| **Summe** | **12** | **4** | **−4** | **−2** |

Gelesen: In α_G^{1/2} steht α⁸ (Photosphäre) · α⁻⁴ (Zündung) · α² (Thomson) = α⁶ und m_e⁴ · m_e⁻² = m_e². Die Verdopplung auf 12 und 4 kommt daher, dass die Gravitation nur als Wurzel m_p/m_Pl = α_G^{1/2} in den Stern eingeht: G M_* m_p = ħc·m_Pl/m_p. Die 4 im Exponenten von T_eff (Stefan-Boltzmann) und die beiden Potenzen 2 der EM-Skalen (α² in Ionisation und Zündung) sind die eigentlichen Zahlenquellen.

In Längen: Die Photosphärenschwelle ist ħc·r_e²/λ̄_e³, die Zündskala ħc·r_e²/(λ̄_e²λ̄_p), die Streufläche r_e². Die Faltung tritt nur als Δs/Δs_* und als δ = 3Δs/R [korr.] auf.

**2.7 Die Grenzmasse.** Aus T_eff⁴ ∝ μ_M folgt M_div/M_* ∝ 𝒞^{−1/2}. Carters Beziehung sagt also: **Der Stern, dessen Hülle gerade von radiativ zu konvektiv wechselt, hat die Faltung Δs_* der Sternskala.** Beobachtet ist M_div ≈ 1,2–1,4 M☉ = 0,65–0,75 M_* (T_eff ≈ 6200–7000 K) [G].

## 3. Die einfachste Form: optische Tiefe gleich Temperaturverhältnis hoch vier

Ohne Strukturzahlen lässt sich dieselbe Bedingung als Gleichheit zweier Abflussraten desselben Strahlungsvorrats aT_c⁴R³ schreiben:

```
Γ_diff = c/(τR),  τ = 3Nσ_T/(4πR²)     (transversal: Diffusion durch Thomson-Streuung; homogene Kugel)   [korr.]
Γ_surf = (c/R)·(T_ion/T_c)⁴            (Abstrahlung einer Oberfläche an der Wasserstoffschwelle)

Γ_diff = Γ_surf   ⟺   τ = (T_c/T_ion)⁴   ⟺   α_G ∝ α¹²(m_e/m_p)⁴
```

In Faltungslängen: τ = (Δs/Δs_p)·2r_e²/R² [korr.] und T_c/T_ion = (f_nuc/f_s)·λ̄_e/λ̄_p. Die Exponenten sind dieselben (Skript Teil 3). **Die Carter-Bedingung ist: Die optische Tiefe des Sterntrichters ist die vierte Potenz des Verhältnisses von Kern- zu Atomtemperatur.** Bei der Sonne ist τ_T = 3,9·10¹⁰ (mit τ = 3Nσ_T/(4πR²); Nσ_T/R² gäbe 1,63·10¹¹ und den Faktor 340) [korr.] und (T_c/T_eff)⁴ = 5,5·10¹³. Den Faktor 1400 dazwischen tragen Temperaturgradient, Zentralkonzentration und die reale Opazität; sie sind in λ gebündelt.

## 4. Der Vorfaktor 2,94 mit Sonnenzahlen

**4.1 Exakte Identität.** L = 4πR²σ_SB T_eff⁴ lässt sich ohne jede Näherung umschreiben:

```
√𝒞 = (2π³/45) · f_s⁴ · (Δs/Δs_*) · [ (r_e/λ̄_e)² / (3Δs/R) ]² · (L_Edd/L)
```

Die Zentraltemperatur fällt heraus (θ/f_nuc = α²/δ). Mit den Sonnenzahlen f_s = 0,01828, Δs☉/Δs_* = 0,5396, δ☉ = 2,1226·10⁻⁶ [korr.], L☉/L_Edd = 3,045·10⁻⁵ ergibt sich 𝒞 = 2,9441, gleich dem direkten Wert. Die Identität gilt für jeden Stern mit seinen eigenen Zahlen. Sie prüft nichts, sie zeigt aber, woraus 2,94 besteht.

**4.2 Zerlegung mit dem Eddington-Modell.** Mit λ = (κ_T/κ̄)(1 − β) und der Eddington-Quartik zerfällt √𝒞 = 1,7158 in sechs Gruppen:

| Gruppe | Herkunft | Sonne (μ_mol = 0,60 / 0,62) |
|---|---|---|
| A = 2π³/45 | Stefan-Boltzmann und Thomson-Geometrie | 1,378 |
| B = C_E² | Polytrope n = 3 | 94,6 |
| C = κ̄/κ_T | Opazität über der Thomson-Opazität | 13,1 / 14,9 |
| D = θ²/μ_mol⁴ | Virialzahl und Molekulargewicht | 3,56 / 3,12 |
| E = f_s⁴/f_nuc² | Saha- und Gamow-Logarithmus | 1,52·10⁻⁴ |
| F = 1/μ_M | Masse in Einheiten M_* | 1,853 |

Das Produkt ist 1,713; mit dem exakten statt des genäherten 1 − β wird es 1,7158.

Die Gruppe C ist die einzige, die nicht aus einer Formel folgt. Die Identität verlangt eine mittlere Opazität κ̄ = 5,2–5,9 cm²/g. Die Sonne hat im Zentrum rund 1,2 cm²/g und an der Unterkante der Konvektionszone rund 20 cm²/g [G]. Der geforderte Wert liegt in diesem Bereich. Das ist die eigentliche numerische Prüfung: Mit realen Opazitäten ist der Vorfaktor verträglich.

**4.3 Was die 2,94 bedeutet.** Die Nähe von 𝒞 zu 1 ist eine Kompensation: A·B·D·F ≈ 8·10² und C ≈ 14 stehen gegen E ≈ 1,5·10⁻⁴. E besteht aus Logarithmen. ℓ_S = 1/(2f_s) = 27 [korr.] gehört zur Photosphäre, Λ_G = 13,5 zur Zündung. Beide hängen nur logarithmisch von den Konstanten ab und bleiben beim Potenzzählen unsichtbar. Eine Herleitung der Zahl 2,94 gibt es deshalb nicht. Hergeleitet sind die Exponenten und die Aussage M_div ~ M_*.

**4.4 Empfindlichkeit.** Wegen 𝒞 ∝ f_s⁸ entspricht der Faktor 2,94 nur 14 % in der Grenztemperatur. In der Grenzmasse sind es Faktor 1,72 (Thomson) bzw. 1,20 (Kramers). Setzt man alle Strukturzahlen auf 1, ergibt die reine Skala T_eff(M_*) = 28,7 eV = 1,06·α²m_ec².

**4.5 Grenzmasse aus der Sonne.** Mit festen Strukturzahlen (Homologie) gilt T_eff ∝ M^{1/4} bei Thomson- und ∝ M^{3/4} bei Kramers-Opazität. Für T_crit = 6200 / 6500 / 7000 K folgt M_div = 1,33 / 1,61 / 2,16 M☉ (Thomson) und 1,10 / 1,17 / 1,29 M☉ (Kramers). Beobachtet sind 1,2–1,4 M☉ [G].

## 5. Opazitätsgesetz, Hayashi-Grenze, Massenfenster

**5.1 Die 12 hängt nicht am Opazitätsgesetz, die 4 schon.** Sterne um 1 M☉ werden innen von Kramers-Opazität (frei-frei, gebunden-frei) beherrscht, nicht von Thomson-Streuung. Mit κ̄/κ_T = C_K·α·n λ̄_e³·(m_ec²/kT_c)^{7/2} und n ∝ N/R³ ergibt die Monom-Rechnung

```
κ̄/κ_T ∝ (m_e/m_p)^{1/2} f_nuc^{-1/2} θ^{-3} μ_M^{-2}          (α kürzt sich, weil kT_c ∝ α² m_p)
α_G ∝ α¹² (m_e/m_p)⁵ · μ_M⁻⁶ · …                                 (M_div ∝ 𝒞'^{-1/6})
```

Der Exponent 12 von α bleibt erhalten. Der Exponent von m_e/m_p wird 5 statt 4, und α_G/(α¹²(m_e/m_p)⁵) = 5406. Carters Form mit der 4 ist die Thomson-Idealisierung. Weil m_e/m_p in unserer Welt fest ist, lassen sich 4 und 5 an einem einzelnen Stern nicht unterscheiden. Der Unterschied steckt dann ganz im Vorfaktor, und der hängt bei Kramers zusätzlich über die gebunden-freie Absorption an der Metallizität.

**5.2 Hayashi-Grenze.** Ein voll konvektiver Stern (Polytrope n = 1,5) mit H⁻-Opazität κ ∝ ρ^{1/2}T⁹ an der Photosphäre P = (2/3)g/κ hat T_eff ∝ M^{7/51}L^{1/102}. Das Skript löst das lineare System und erhält genau diese Exponenten. Die Temperatur ist praktisch fest, ihre Skala ist die H⁻-Bindung 0,754 eV [G] = 0,0277·α²m_ec². Die Hayashi-Linie bei 3000–4500 K [G] entspricht f = 0,010–0,014. Die Carter-Grenze liegt bei f_s ≈ 0,020, die halbe Ionisation bei Photosphärendichte bei f = 0,044. Alle drei sind Bruchteile derselben EM-Skala α²m_ec²; die Bruchteile sind Logarithmen und Tiefenstruktur.

**5.3 Massenfenster.** Die untere Grenze (Elektronenentartung vor der Zündung, kT_max ~ α_G² N^{4/3} m_ec²) hat die Skala (f_nuc α² m_p/m_e)^{3/4} M_* = 0,012 M_*; beobachtet sind 0,075 M☉ [G] = 0,040 M_*. Die obere Grenze (β = 1/2) liegt bei 72 M_*. Die schwache Form der Beziehung, „die Grenzmasse liegt innerhalb des Sternfensters", verlangt bei Thomson-Skalierung nur 2·10⁻⁴ < 𝒞 < 9·10². Das ist eine Breite von 6,6 Dekaden in α_G. Die scharfe Form „M_div ≈ M_*" ist eine Aussage über Faktoren von etwa 2.

## 6. Die Konvektionsschwelle als Entstehen eines Umlaufs

Für reinen Wasserstoff ist

```
∇_ad = [2 + x(1−x)(5/2 + E_ion/kT)] / [5 + x(1−x)(5/2 + E_ion/kT)²] ,       E_ion = α²m_ec²/2   [korr.]
```

Bei Photosphärendichte n_H = 10²³ m⁻³ (Saha):

| T | f = kT/(α²m_ec²) | Ionisationsgrad x | ∇_ad |
|---:|---:|---:|---:|
| 5000 K | 0,016 | 1,3·10⁻⁵ | 0,399 |
| 8000 K | 0,025 | 6,8·10⁻³ | 0,258 |
| 10000 K | 0,032 | 0,056 | 0,131 |
| 12000 K | 0,038 | 0,22 | 0,100 |
| 14000 K | 0,044 | 0,50 | 0,104 |
| 20000 K | 0,063 | 0,96 | 0,269 |

Der radiative Gradient lautet in denselben Größen

```
∇_rad = (λ/4) · (κ/κ_T) · (P_gas/P_rad) ,        P_gas/P_rad ∝ (Δs_*/Δs)²
```

Er trägt die Faltung (λ und das Druckverhältnis) und die EM-Opazität. ∇_ad trägt allein die atomare Bindung E_ion [korr.]. Konvektion setzt ein, wenn ∇_rad > ∇_ad. Das Quadrat der Auftriebsfrequenz ist N² ∝ ∇_ad − ∇. Im radiativen Bereich ist N reell: Ein verschobenes Gaselement schwingt, das ist ein tangentialer Umlauf (Schwerewelle). An der Schwelle geht N durch null. Jenseits davon wächst die Verschiebung an und schließt sich zu einer Zirkulation, einem neuen Umlauf.

In der Zustandsraum-Sprache ist das eine Änderung der Zahl zugänglicher Umläufe n(p). Der Parameterpunkt p = (Δs, R, T_c) wird von der Faltung gesetzt. Die Schwellfläche setzen die EM-Größen: σ_T bzw. κ im Transport, E_ion [korr.] im adiabatischen Gradienten. Carters Beziehung sagt, dass die Faltung Δs_* der Sternskala auf diese Schwellfläche fällt.

## 7. Liefert das Einrasten eine gleichwertige Bedingung?

**7.1 Ratenform.** Mit Γ_G = α_G ω_p = 8,42·10⁻¹⁵ s⁻¹, ω_p = m_pc²/ħ, ω_e = m_ec²/ħ und der Umlauffrequenz ω_a = α²ω_e = 4,13·10¹⁶ s⁻¹ lautet Carters Beziehung exakt

```
Γ_G · ω_p³  =  𝒞 · (α ω_a)⁴            (2,44·10⁵⁸ = 2,944 × 8,28·10⁵⁷)
```

Die Frequenz (Γ_G ω_p³)^{1/4} = α_G^{1/4}ω_p = 3,95·10¹⁴ s⁻¹ ist kein Zyklus eines Wasserstoffsystems. Sie entsteht nur als Kombination in der Sternkette.

**7.2 Einzelnes Wasserstoffsystem.** Das Einrasten dΔθ/dt = Δω − K sin Δθ verlangt K ≥ |Δω|, einen linearen Vergleich zweier Raten. Setzt man K = Γ_G gegen die Verstimmung Δω ~ ω_a, ist K/Δω = 2,0·10⁻³¹. Ein Einrasten der Faltungsrate an einen EM-Zyklus eines Atoms gibt es nicht. Allgemeiner: Die Raten des Wasserstoffs haben die Form α^n(m_e/m_p)^k ω_e mit Exponenten, die aus der Atomphysik folgen:

| Rate | n | k | Wert [s⁻¹] |
|---|---:|---:|---:|
| Umlauf (Rydberg) α²ω_e/2 | 2 | 0 | 2,07·10¹⁶ |
| Feinstruktur 2p α⁴ω_e/32 | 4 | 0 | 6,9·10¹⁰ |
| Lyman-α-Zerfall (2/3)⁸α⁵ω_e | 5 | 0 | 6,27·10⁸ |
| 2s→1s, zwei Photonen | 8 | 0 | 8,23 |
| Hyperfeinaufspaltung (4/3)g_pα⁴(m_e/m_p)ω_e | 4 | 1 | 8,9·10⁹ |
| 21-cm-Zerfall (64/81)g_p³α¹³(m_e/m_p)³ω_e | 13 | 3 | 2,87·10⁻¹⁵ |

Die Gleichung Γ_G = Γ_X ergibt Carter genau dann, wenn Γ_X ∝ α¹²(m_e/m_p)³ω_e ist, also (n, k) = (12, 3). Diese Stelle ist in der Tabelle nicht besetzt.

**7.3 Ein naher Treffer, der keiner ist (Look-elsewhere).** Die 21-cm-Linie liegt mit Γ_G/A₁₀ = 2,935 zahlenmäßig fast bei 𝒞 = 2,944. Die Zerlegung zeigt den Grund:

```
Γ_G/A₁₀ = 𝒞 · 81/(64·g_p³·α) ,        (64/81)·g_p³·α = 1,005
```

Der g-Faktor des Protons (eine QCD-Zahl) liegt in der dritten Potenz zufällig bei 81/(64α), und der Exponent von α ist 13 statt 12. Bei sechs geprüften Raten ist ein Treffer auf 0,3 % in einer Größe, die zwei unabhängige Konstanten zufällig verknüpft, ohne physikalischen Gehalt. Die 21-cm-Linie wird hier nicht verwendet.

**7.4 Was stattdessen gilt.** α_G erscheint in der Sternkette nur über G M m_p = μ_M·α_G^{−1/2}·ħc, also über die kollektive Teilchenzahl N_* = α_G^{−3/2}. Die Faltung einzelner Protonen koppelt an keinen einzelnen EM-Zyklus. Die Faltung des ganzen Trichters setzt den Parameterpunkt. Die Bedingung, die Carter ergibt, hat deshalb zwei Teile:

1. **Ratengleichheit** zweier Abflüsse desselben Vorrats (Abschnitt 3): Γ_diff (transversal, σ_T ∝ r_e²) = Γ_surf bei der Wasserstoffschwelle (tangential gesetzt, α²m_ec²).
2. **Schwelle** eines tangentialen Umlaufs (Abschnitt 6): N² ∝ ∇_ad − ∇_rad geht durch null.

Die Schwelle hat dieselbe mathematische Klasse wie die Einrastschwelle. Beim Einrasten verschwindet die Schlupffrequenz Ω = √(Δω² − K²) bei K = |Δω|. Bei der Konvektion verschwindet die Auftriebsfrequenz bei ∇_rad = ∇_ad. Die Rollen sind aber andere. Beim Einrasten fallen zwei bestehende Umläufe zu einem zusammen (T² → T¹). Bei der Konvektion wird eine Schwingung durch eine Zirkulation ersetzt, und die Zahl der Umläufe ändert sich über den Parameter. Eine Einrastbedingung „Faltungsrate gegen EM-Rate" ergibt Carters Exponenten nicht.

## 8. Beobachterfehler-Selbstcheck

Die Beziehung α_G ≈ α¹²(m_e/m_p)⁴ enthält nur Eigenschaften von Proton und Elektron. Man könnte sie deshalb als lokale Eigenschaft des Protons lesen, etwa als „Gleichgewicht im einzelnen Proton". Die Rechnung zeigt etwas anderes: Die Faltung geht nur als kollektive Größe N_* = α_G^{−3/2} bzw. Δs_* ein, und die Bedingung ist eine Eigenschaft des ganzen Trichters. Sie sagt, wo seine Transportregel (Diffusion oder Zirkulation) umschlägt. Die Deutung als Eigenschaft eines einzelnen Objekts wäre die Verwechslung einer globalen Regeländerung mit einem lokalen Träger.

## 9. Was hergeleitet ist und was nicht

| Aussage | Art |
|---|---|
| Exponenten 12 und 4 aus Hydrostatik + Zündung (α²m_p) + Thomson-Transport (r_e²) + Stefan-Boltzmann + Photosphärenschwelle (α²m_e) + Grenzmasse ~ M_* | Herleitung (Potenzrechnung, im Skript als Monom-Rechnung) |
| Herkunft: α-Exponent 12 = 2·(8 − 4 + 2), m_e/m_p-Exponent 4 = 2·(4 − 2); die 2 aus α_G^{1/2} = m_p/m_Pl | Herleitung |
| δ_* = (f_nuc/θ)(r_e/λ̄_e)² [korr.]: Das Faltungsdifferential brennender Sterne ist EM-bestimmt | Herleitung (Hydrostatik + Gamow), Sonne 2,12·10⁻⁶ |
| L_Edd = (9/2)m_pc³Δs/r_e², 1 − β ∝ (Δs/Δs_*)², τ = (Δs/Δs_p)·2r_e²/R² [korr.] | exakte Umformungen |
| Carter ⟺ τ_Trichter = (T_c/T_ion)⁴ ⟺ M_div ~ M_* | Herleitung bis auf Strukturzahlen |
| α-Exponent 12 auch bei Kramers-Opazität; m_e/m_p-Exponent dann 5 | Herleitung |
| 𝒞 = 2,94 als Produkt von Strukturzahlen der Sonne; verlangt κ̄ ≈ 5–6 cm²/g, im solaren Bereich 1,2–20 cm²/g [G] | exakte Identität + Abschätzung |
| Zahlenwert 2,94 | nicht hergeleitet; Kompensation von Logarithmen (≈ 10⁻⁴) gegen Struktur- und Opazitätszahlen (≈ 10⁴) |
| Grenzmasse aus Homologie: 1,1–2,2 M☉ gegen beobachtet 1,2–1,4 M☉ [G] | Abschätzung |
| Einrasten Faltungsrate gegen EM-Rate eines Wasserstoffsystems | ergibt Carter nicht (K/Δω = 2·10⁻³¹; Gitterpunkt (12, 3) unbesetzt) |
| Konvektionsschwelle als N² = 0, gleiche Klasse wie die Einrastschwelle | strukturelle Entsprechung, keine Herleitung aus der Torus-Gleichung |
| Warum unser 𝒞 nahe 1 liegt | offen; das Modell hat dafür bisher keine Gleichung |

## Literatur

- B. Carter, *Large number coincidences and the anthropic principle in cosmology*, in: M. S. Longair (Hrsg.), Confrontation of Cosmological Theories with Observational Data, IAU Symp. 63, Reidel 1974, S. 291–298.
- B. J. Carr, M. J. Rees, *The anthropic principle and the structure of the physical world*, Nature 278 (1979) 605–612.
- W. H. Press, A. P. Lightman, *Dependence of macrophysical phenomena on the values of the fundamental constants*, Phil. Trans. R. Soc. A 310 (1983) 323–336.
- J. D. Barrow, F. J. Tipler, *The Anthropic Cosmological Principle*, Oxford University Press 1986, Kap. 5.
- A. S. Eddington, *The Internal Constitution of the Stars*, Cambridge University Press 1926.
- C. Hayashi, *Stellar evolution in early phases of gravitational contraction*, PASJ 13 (1961) 450–452.
- R. Kippenhahn, A. Weigert, A. Weiss, *Stellar Structure and Evolution*, 2. Aufl., Springer 2012.
- G. B. Rybicki, A. P. Lightman, *Radiative Processes in Astrophysics*, Wiley 1979.
- J. N. Bahcall, A. M. Serenelli, S. Basu, *New solar opacities, abundances, helioseismology, and neutrino fluxes*, ApJ 621 (2005) L85.
- R. Adler, *A study of locking phenomena in oscillators*, Proc. IRE 34 (1946) 351–357.
