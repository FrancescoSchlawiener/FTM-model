# Umlaufsinn der Rückkehr

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „der eine Sinn ist gleichwertig mit Halbgerade, Holomorphie in Im τ > 0 und positiver Frequenz“ → festgelegt durch die Erneuerung G(mL) = q^m zusammen mit χ → 1; die Halbgerade allein legt ihn nicht fest, Kausalität und Holomorphie erzwingen ihn nicht, τ_D ist keine Periode in imaginärer Zeit; offene Frage „warum Halbgerade“ → beantwortet (eine Schleife der Periode P verlangte e^{−|k|P} = 1, nur k = 0); Drehimpuls L → J.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/umlaufsinn.py` (Ausgabe `umlaufsinn.txt`)
> s = √(a_N/a₀), q = e^{−s}, ν = a/a_N = Σ_{m≥0} qᵐ. Offen war, warum nur Rückläufe in einem Sinn zählen.

## 0  Statisch zählt das Gewicht, nicht die Richtung

```
1/(1 − q) = ½ + ½·coth(s/2)        ein Sinn mit vollem Gewicht = beide Sinne mit halbem Gewicht
```

| a_N/a₀ | 10 | 3 | 1 | 0,3 | 0,1 | 0,01 |
|---|---|---|---|---|---|---|
| Σ_{m≥0} qᵐ | 1,044 | 1,215 | 1,582 | 2,371 | 3,689 | 10,51 |
| Σ_{m∈ℤ} q^{\|m\|} = coth(s/2) (v⁴ = 4GMa₀) | 1,088 | 1,430 | 2,164 | 3,742 | 6,377 | 20,02 |
| beide Sinne, gleicher Nullpunkt: coth s | 1,004 | 1,065 | 1,313 | 2,005 | 3,267 | 10,03 |

Zerlegt man die Summe in den geraden (m ↔ −m) und den ungeraden Teil, trägt statisch nur der gerade: Re R = ½(1 + coth(s/2)) auf 4·10⁻¹⁶, der ungerade Teil, der den Sinn trägt, ist bei ω = 0 null. Festgelegt ist das Gesamtgewicht q je Rückkehrschritt; coth entsteht nur, wenn der zweite Sinn mit vollem Gewicht hinzukommt.

## 1  Kausalität

(I − qT)⁻¹ = Σ_{m≥0} qᵐTᵐ ist die einzige beschränkte Inverse des Operators erster Ordnung I − qT (die andere Reihe wächst wie q^{−m}, bei m = 50 auf 5·10²¹). Symbol 1/(1 − qz), Pol bei z = e^s außerhalb der Einheitsscheibe: Träger m ≥ 0, mit z = e^{iωτ_R} analytisch für Im ω > 0. Die beidseitige Summe ist die Green-Funktion des symmetrischen Operators zweiter Ordnung, die einseitige dessen kausaler Wiener-Hopf-Faktor:

```
Σ_{m∈ℤ} q^{|m|}Tᵐ = (1 − q²)·[(1 + q²) − q(T + T⁻¹)]⁻¹        diskret (−D∂² + Γ)⁻¹, Pole bei z = e^{±s}
```

**Modell:** Diffusion mit Verlust Γ auf der Geraden (Überlagerung des Rings), Einheit 1/Γ, L²/D = s². Rücklauf m in einem Sinn = erstes Erreichen von +mL, f_m(t) = ms/√(4πt³)·e^{−m²s²/4t − t}, ∫f_m = e^{−ms}. Die Folgedichte Σf_m hat Träger t ≥ 0 und die Transformierte

```
R(ω) = 1/(1 − e^{−s√(1 − iω/Γ)})        Singularitäten bei ω = −iΓ(1 + 4π²k²/s²), alle unten
```

| s = 1 | statisch | Kramers-Kronig ω₀ = 0 | ω₀ = Γ |
|---|---|---|---|
| ein Sinn | 1,582 | 0,58198 = 0,58198 | 0,36752 = 0,36752 |
| Ring, beide Sinne | 2,164 | 1,16395 = 1,16395 | 0,73504 = 0,73504 |
| zeitsymmetrisch r(\|t\|) | 2,164 | 1,16395 ≠ 0 | 0,73504 ≠ 0 (Pol bei +iΓ) |

Die einseitige Summe ist also die retardierte Resolvente. Die Windungssumme beider Sinne ist aber ebenso retardiert: Beide Umlaufsinne einer räumlichen Schleife liegen in der Zukunft, und statisch hat sie denselben Wert coth(s/2) wie die zeitsymmetrische Fassung. Kausalität erzwingt einen Sinn nur, wenn T die Geschichte verschiebt.

**Echtzeit:** Mit Γ = a_N/c und τ_D = c/a₀ = 79 Mrd. J. dauert der m-te Rücklauf im Mittel m·τ_D/(2s), bei a_N = a₀ also 40 Mrd. J. Bis 13,8 Mrd. J. wären ν = 1,08 (a_N = a₀), 1,09 (0,1 und 0,01 a₀) aufgebaut statt 1,58, 3,69, 10,5. Als Vorgang in Echtzeit trägt die Kausalität nicht; die Summe muss stationär sein. Stationär ist (D∂² − Γ)⁻¹ = ∫₀^∞ e^{t(D∂² − Γ)}dt in t stets einseitig, im Raum aber beidseitig. Einseitig wird es nur auf einer Halbgeraden ab dem Ort (G(0) = 1, G beschränkt, (∂ + √(Γ/D))G = 0): Σ_{m≥0} G(mL) = 1,582, bei einer Quelle in der vollen Geraden Σ_{m∈ℤ} = 2,164 (finite Differenzen).

## 2  Positive Frequenz (KMS)

Oszillator im Kreis imaginärer Zeit, Umfang β, x = βω (auf 6 Stellen, n_B aus der Spur):

```
G_β(τ) = Σ_{n∈ℤ} e^{−ω|τ+nβ|}/2ω = [(Σ_{n≥0} e^{−nx})·e^{−ωτ} + (Σ_{n≥1} e^{−nx})·e^{+ωτ}]/2ω = [(1 + n_B)e^{−ωτ} + n_B e^{ωτ}]/2ω
```

Windungen in + Richtung gehören zur positiven Frequenz, in − Richtung zur negativen; KMS: Gegensinn = e^{−x} × Sinn. Bei x = 1: Besetzung n_B = Σ_{n≥1} = 0,582; Emissionsgewicht G^>/ρ = Σ_{n≥0} = 1,582; Propagator G^K/ρ = 1 + 2n_B = Σ_{n∈ℤ} = 2,164; Kommutator G^> − G^< = 1; statische Suszeptibilität ∫₀^β G_β = 1/ω² für jedes β.

Der Propagator ist also beidseitig, die Besetzung einseitig. Mit x = s stellt die Rückkopplung das Emissionsgewicht einer Ordnung dar: ν = 1 + n_B(s), Zusatz a − a_N = a_N·n_B(s), χ = s·n_B(s). Die retardierte Antwort trägt keinen Bose-Faktor. Der eine Sinn ist hier ein Vorzeichen der Frequenz; der Gegensinn ist nicht null, sondern q-mal so stark, und mitgezählt gibt er coth(s/2).

## 3  Orientierung des Torus

q = e^{2πiτ}, Im τ = s/2π. f = Σ_{n≥0} qⁿ hat längs Re τ nur Fourier-Anteile n ≥ 0 (die negativen ≤ 2·10⁻¹⁷), und Im f ist die Hilbert-Transformierte von Re f (6·10⁻¹⁶): Holomorphie in Im τ > 0 ist dieselbe Struktur wie Kausalität, mit der Windungszahl als Zeit. Den einen Sinn legt Holomorphie aber nicht fest: Sie wählt nur q statt q̄; der Schwarz-Kern (1 + q)/(1 − q) ist ebenfalls holomorph und gibt 1 + 2n_B. [korr.] Die reell-analytische Reihe Σqⁿ + Σq̄ⁿ = (1 − |q|²)/|1 − q|² ist beidseitig und reell.

| a_N/a₀ (Re τ = 0) | 3 | 1 | 0,3 | 0,1 | 0,01 |
|---|---|---|---|---|---|
| f = ½(f + f̄) | 1,215 | 1,582 | 2,371 | 3,689 | 10,51 |
| Σqⁿ + Σq̄ⁿ = coth(s/2) | 1,430 | 2,164 | 3,742 | 6,377 | 20,02 |
| \|f\|² | 1,476 | 2,503 | 5,622 | 13,61 | 110,4 |

Die Orientierung selbst ist in a unsichtbar: τ → −τ̄ (q → q̄) lässt |q| gleich, und die symmetrisierte Fassung ½(f + f̄) ist auf der Achse gleich f. Die reell-analytische Fassung mit vollem Gewicht gibt v⁴ = 4GMa₀, mit gleichem Nullpunkt coth s (bei a_N = a₀ 17 % tiefer); das Produkt |f|² gibt a → a₀ ohne Tully-Fisher.

## 4  Drehsinn der Scheibe

Der Kern des Autors mit dem Radialprojektor e_r⊗e_r gibt auf der Kreisbahn a_r = −κ₀rω²/(λ_m² + ω²), gerade in ω (λ_m = 0,3: −0,91743 für ω = ±1). Einen ungeraden Anteil gäbe nur ein antisymmetrischer Projektor e_r⊗e_φ − e_φ⊗e_r (Casimir-Typ, ∓0,27523); er fehlt.

**Fall Z (Sinn = Zeit):** ν = 1/(1 − q) für jede Komponente und jedes System. **Fall D (Sinn = Drehsinn):** D1: nur Rückläufe im Drehsinn der Scheibe zählen, Gegenläufer laufen Newton; D2: getragen vom Nettodrehimpuls, ν = 1 + λ·q/(1 − q), λ = |J_netto|/Σ|J|, Tiefgrenze v⁴ = λ²GMa₀.

| r/r* (Punktmasse) | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| ν (Z, alle Komponenten) | 1,582 | 2,541 | 3,528 | 5,517 | 10,51 |
| D1 v_gegen/v_mit = D2 v(λ=0)/v(λ=1) | 0,795 | 0,627 | 0,532 | 0,426 | 0,308 |
| D2, zwei Scheiben 70:30 (λ = 0,4) | 0,883 | 0,798 | 0,755 | 0,713 | 0,676 |

Druckgestützt (Plummer, isotrop, Virialsatz), σ₁ in km/s, Z gegen D mit λ = 0 (Newton): 10⁷ M☉, b = 0,3 kpc: 9,79 gegen 3,75 (Faktor 2,61); 10⁹ M☉, 1 kpc: 33,4 gegen 20,6; 10¹¹ M☉, 3 kpc: 133,7 gegen 118,7. Tiefgrenze Z: σ₁⁴ = (4/81)GMa₀.

**Trennbeobachtung:** gegenläufige Sternscheiben (Z: gleiche Kreisgeschwindigkeit beider Komponenten nach Korrektur der asymmetrischen Drift; D1: 0,53 : 1 bei 3r*), Tully-Fisher-Nullpunkt gegen Drehanteil (Z: unabhängig; D2: v ∝ √λ), druckgestützte Zwerge im Tiefbereich (Z: σ₁⁴ = (4/81)GMa₀; D: Newton).

**Dissipation und Casimir-Anteil:** R = Re R (gerade, reaktiv) + i·Im R (ungerade, absorptiv). Summenregel R(0) − 1 = (1/π)P∫Im R/ω dω = 0,58198 (Ring 1,16395, zeitsymmetrisch 0): Bei Kausalität ist die statische Verstärkung das Integral des absorptiven Spektrums. Γ > 0 (dS/dt ≥ 0) macht die Reihe konvergent (Γ → 0: ν → ∞) und wählt einen Zeitsinn, wenn T die Geschichte verschiebt; einen Raumsinn wählt Γ nicht, weil die Diffusion spiegelsymmetrisch ist. Das könnte nur ein antisymmetrischer Casimir-Anteil, der an den Drehsinn Ω koppelt; er fehlt im Kern.

## Ergebnis

- Statisch gilt 1/(1 − q) = ½ + ½coth(s/2): Die Form legt das Gesamtgewicht q je Rückkehrschritt fest, nicht die Richtung.
- Die einseitige Summe ist die retardierte Resolvente von I − qT (Träger t ≥ 0, analytisch oben, Kramers-Kronig exakt). Kausalität erzwingt sie nur, wenn T die Geschichte verschiebt; auf einer räumlichen Schleife ist auch die Summe beider Sinne retardiert (coth(s/2)).
- In Echtzeit wäre die Summe heute nur zu ν ≈ 1,08 aufgebaut; sie ist stationär. Stationär ist der eine Sinn verträglich mit einer Halbgeraden der inneren Richtung ab dem Ort, mit Holomorphie in Im τ > 0 und mit positiver Frequenz (ν = 1 + n_B(s)); festgelegt wird er durch keine davon allein, sondern durch die Erneuerung G(mL) = q^m aus der Halbgeraden zusammen mit dem Grenzwert χ → 1 des Autors (beide Sinne gäben im tiefen Bereich χ → 2). Kausalität erzwingt ihn nicht (beide Windungssinne sind retardiert), Holomorphie auch nicht (der Schwarz-Kern (1 + q)/(1 − q) ist holomorph und gibt 1 + 2n_B); τ_D ist die Laufzeit einer Windung, keine Periode in imaginärer Zeit. [korr.]
- Der Kern ist gerade im Drehsinn, die Orientierung des Torus geht über |q| nicht ein: Der Sinn ist universell (Fall Z). Fall D trennt sich an gegenläufigen Scheiben, am Tully-Fisher-Nullpunkt gegen λ und an druckgestützten Zwergen.
- **Beantwortet:** warum die Rückkehr auf einer Halbgeraden läuft und nicht auf einer geschlossenen Schleife: Die innere Richtung ϖ ist eine Glättungstiefe, Glättungen setzen sich zusammen (e^{−|k|a}·e^{−|k|b} = e^{−|k|(a+b)}), und eine Schleife der Periode P verlangte e^{−|k|P} = 1, was nur für k = 0 (den Mittelwert) geht; ϖ ≥ 0 ist daher eine Halbgerade mit dem Ort als unterster Ebene. [korr.]
