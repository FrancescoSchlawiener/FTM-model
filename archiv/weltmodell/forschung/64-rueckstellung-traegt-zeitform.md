# Rückstellung aus gebundener Faltung und Zeitform der Dunklen Energie

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: M_Q = 2,41·10¹⁴ M☉ als Wert für a₀ = cH₀/2π gekennzeichnet, mit a₀ = 1,2·10⁻¹⁰ 3,7·10¹⁴ M☉ (M_Q ∝ a₀³); Umbenennungen s₀ = d ln a₀/d ln a → m₀, Ω_Λ → Ω_DE. Die Notiz rechnet mit a₀ nach einem Horizont (Hubble-Radius: a₀ = cH₀/2π = 1,04·10⁻¹⁰); 1 − w = 2,8·10⁻³ und M_Q = 2,41·10¹⁴ M☉ gelten für cH₀/2π.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/rueckstellung_zeitform.py` (Ausgabe `rueckstellung_zeitform.txt`)
> Eingaben wie im Modell: H₀ = 67,4 km/s/Mpc, Ω_m = Ω_b(1 + Q) = 0,3134, Ω_r = 9,14·10⁻⁵, k_s = π²/45; lineares Wachstum = Standardmodell mit σ₈ = 0,811, n_s = 0,965, Transferfunktion Eisenstein-Hu, Sheth-Tormen-Parameter 0,3222 / 0,707 / 0,3 **[G]**. Einzige Normierung E(0) = 1. **[G]** = aus dem Gedächtnis, nicht geprüft.

## Welche Faltung innen zählt

In einem mitbewegten Gebiet ist ΣΔs = GM/(3c²) fest. Zählte alle Materie, bliebe ΔL_außen je Gebiet fest und ΔL/R fiele wie 1/a. Wachsen kann nur der Anteil, der nicht mehr mitatmet. Das Kriterium dafür steckt im Gewicht der Schließungsregel, w = θ₊²/(θ₊² + θ_c²), θ_c = a₀/c: w = ½ bei θ = 3H_lok = θ_c, also H_lok/H = Θ/6π = 0,053. Im Kugelmodell ist das die lineare Schwelle δ_w = 1,041 (Umkehr 1,062, Kollaps 1,686). Der Hintergrund trägt 1 − w = 2,8·10⁻³ und bleibt vernachlässigbar.

Innen zählt also die gebundene Faltung, F = ΣΔs_gebunden/ΣΔs = gebundener Massenanteil (Δs ∝ M). Ansatz für den Beitrag der Rückstellung zu H²:

```
Ω_R(z) = Ω_R0·F(z)/F₀,     E² = R/(1 − k_sΘ²),     R = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_R0·F/F₀,     Ω_R0 = Ω_DE − k_sΘ₀²
```

Je physikalischem Volumen gelesen (F·a⁻³) verdünnt die Rückstellung wie Materie (w = −0,02 … −0,09) und scheidet aus. F wächst monoton, denn gebunden bleibt gebunden; damit gilt w_R ≤ −1 zu jeder Zeit.

F aus Press-Schechter mit δ_w (PS-w) und mit 1,686 (PS-c) sowie Sheth-Tormen (ST), ΛCDM-Hintergrund, Massen als Gesamtmasse (baryonisch M/6,36); σ(M) = 8,09 / 5,85 / 3,86 bei 10⁶ / 10⁸ / 10¹⁰ M☉, M_Q = 2,41·10¹⁴ M☉ mit σ = 0,84 (gerechnet mit a₀ = cH₀/2π; mit a₀ = 1,2·10⁻¹⁰ ist M_Q = 3,7·10¹⁴ M☉, M_Q ∝ a₀³) [korr.]:

| Verfahren | M_min | F(0) | F(1)/F₀ | F(2)/F₀ | w_R(0) | w_R(1) | Rückstellung allein: CPL w₀ / w_a |
|---|---|---|---|---|---|---|---|
| PS-w | 10⁶ | 0,90 | 0,93 | 0,85 | −1,020 | −1,058 | −1,007 / −0,11 |
| PS-w | 10⁸ | 0,86 | 0,90 | 0,78 | −1,029 | −1,085 | −1,009 / −0,17 |
| PS-w | 10¹⁰ | 0,79 | 0,83 | 0,66 | −1,046 | −1,142 | −1,010 / −0,30 |
| PS-w | M_Q(z) | 0,21 | 0,33 | 0,09 | −1,30 | −1,83 | −1,11 / −1,59 |
| ST | 10⁶ | 0,64 | 0,85 | 0,71 | −1,047 | −1,120 | −1,026 / −0,21 |
| ST | 10¹⁰ | 0,48 | 0,71 | 0,47 | −1,091 | −1,259 | −1,030 / −0,51 |
| ST | M_Q(z) | 0,05 | 0,15 | 0,01 | −1,52 | −2,46 | −1,19 / −2,85 |

**Mindestmasse:** Sie setzt nur die Stärke. |1 + w_R(0)| = 0,02–0,05 (10⁶ M☉), 0,03–0,06 (10⁸), 0,05–0,09 (10¹⁰), bei M_Q 0,3–0,7; PS-c liegt zwischen PS-w und ST. Eine andere Transferfunktion (BBKS) ändert F(1)/F₀ um 0,5–0,7 %, w_a des Gesamtbilds um höchstens 0,01. Die Rückstellung allein liegt stets bei w₀ < −1, w_a < 0.

## Gesamtbild, selbstkonsistent

Wachstum, F und H(z) iteriert; der Vakuumteil folgt dem Horizont wie ohne Rückstellung. Spannen über PS-w, PS-c, ST und M_min = 10⁶–10¹⁰ M☉; w_DE ist die gesamte Dunkle Energie (Vakuumteil + Rückstellung), CPL über 0 < z < 2.

| Horizont | Θ₀ | a₀(0) / 1,2·10⁻¹⁰ | w_DE(0) | w_DE(1) | CPL w₀ / w_a | Extremum von ρ_DE |
|---|---|---|---|---|---|---|
| Hubble-Radius | 1 | 0,87 | −0,90 … −1,01 | −0,47 | −1,02 … −1,20 / +1,15 … +1,52 | keines |
| Ereignishorizont | 0,871–0,897 | 0,76–0,78 | −0,99 … −1,05 | −0,95 … −1,06 | −1,02 … −1,08 / +0,07 … +0,16 | Minimum z = 0,05–1,7 (PS-w, 10⁶: keines) |
| Kodama-Hayward | 0,787–0,818 | 0,68–0,71 | −1,04 … −1,10 | −1,03 … −1,21 | −1,09 … −1,15 / −0,09 … +0,17 | Minimum z = 1,2–2,2 |
| de-Sitter-Ende | 0,839–0,875 | 0,73–0,76 | −1,02 … −1,07 | −1,04 … −1,19 | −1,01 … −1,04 / −0,08 … −0,37 | keines |

Mit M_Q(z) wächst die Rückstellung so schnell, dass nur Hubble-Radius (w₀ = −1,8 … −2,5), Ereignishorizont mit PS-w (−1,11 / +0,25) und de-Sitter-Ende mit PS-w (−1,05 / +0,05) einen selbstkonsistenten Zustand haben; die übrigen sieben Fälle haben keinen. In keinem der 41 selbstkonsistenten Fälle hat die Dunkle Energie ein Maximum. DESI DR2 + CMB + SN **[G]**: w₀ / w_a = −0,84 / −0,62, −0,75 / −0,86, −0,67 / −1,09, Maximum bei z = 0,35–0,44.

Alle Fälle enden in de Sitter (F sättigt, sobald das Wachstum steht); die Bedingung des Kreislaufs, fallendes H, ist damit wie bei konstanter Dunkler Energie verletzt.

## Frühe Grenzen

- Hubble-Radius: Ω_s = 0,22 zu jeder Zeit (ΔN_eff = 2,09), Dunkle Energie bei z = 3 fast siebenmal so groß wie bei ΛCDM, D_M(1090) 8–9 % kürzer, σ₈ fällt (Wachstum ab a = 0,01) auf 0,50–0,51.
- Ereignishorizont, Kodama-Hayward, de-Sitter-Ende: Ω_s(1090) = 2,7·10⁻⁵ / 8·10⁻³ / 3·10⁻¹⁰, ΔN_eff = 0.
- Rückstellung: Ω_R(1090) ≈ 0, Ω_R(3) = 0,007–0,020, Ω_R(2) = 0,026–0,051 (Ω_DE bei ΛCDM 0,033 bzw. 0,075) [korr.]. Die gesamte Dunkle Energie liegt bei z = 3 bei 0,5–1,3 des ΛCDM-Werts. Früh ist beides klein.

## Wo es nicht trägt

Mit ρ_DE′ = ρ_s′ + ρ_R′ (′ = d/d ln a, Einheiten der kritischen Dichte heute) und ρ_R′ ≥ 0 verlangt ein Maximum, dass der Vakuumteil heute schneller fällt, als die Rückstellung wächst, früher aber langsamer.

| | Hubble-Radius | Ereignishorizont | Kodama-Hayward | de-Sitter-Ende |
|---|---|---|---|---|
| ρ_R′ heute / bei z = 1 | +0,05 … +0,22 / +0,11 … +0,37 | +0,03 … +0,14 / +0,08 … +0,30 | +0,03 … +0,15 / +0,09 … +0,33 | +0,03 … +0,14 / +0,09 … +0,31 |
| ρ_s′ heute / bei z = 1 | −0,25 … −0,20 / −2,0 … −2,1 | −0,04 / −0,20 | +0,05 … +0,07 / −0,04 … −0,01 | 0 / 0 |
| m₀ = d ln a₀/d ln a geliefert [korr.] | −0,57 … −0,46 | −0,13 … −0,10 | +0,20 … +0,23 | 0 |
| m₀ verlangt (drei DESI-Kombinationen) [korr.] | −0,9 … −2,1 | −1,1 … −2,3 | −1,4 … −2,8 | −1,2 … −2,5 |

- Die Rückstellung hat den richtigen Gang: Ihr Zuwachs ist früh groß und heute klein; sie allein gibt w_a < 0.
- Das Vorzeichen von 1 + w₀ kann nur der Vakuumteil liefern. Bei allen Horizonten fällt er früher schneller als heute (Hubble-Radius, Ereignishorizont) oder wächst heute (Kodama-Hayward) oder steht (de Sitter). Die Summe hat deshalb höchstens ein Minimum, das Spiegelbild der gemessenen Form.
- Verlangte Form des Vakuumteils bei gegebener Rückstellung (Θ₀ = 1): a₀ muss weiterhin ein Maximum haben, bei z = 0,42–0,81, 1,13–1,38-mal der heutige Wert, a₀(1)/a₀(0) = 1,04–1,35. Die Rückstellung verschiebt das verlangte Maximum von a₀ nur später und höher.

## Ergebnis

1. Innen zählt die gebundene Faltung; das Gewicht der Schließungsregel legt die Schwelle fest (δ_w = 1,04, nahe der Umkehr). F wächst monoton, die Rückstellung ist daher zu jeder Zeit w_R ≤ −1.
2. Die Rückstellung allein trägt das Vorzeichen w_a < 0 der Messung, aber w₀ < −1 (−1,01 … −1,03 bei M_min = 10⁶–10¹⁰ M☉, bis −1,2 bei M_Q). Die Mindestmasse ändert nur die Stärke.
3. Mit dem Vakuumteil nach einem der vier Horizonte hat die Dunkle Energie in keinem Fall ein Maximum. Die drei früh verträglichen Horizonte liegen bei w₀ < −1; der Hubble-Radius bei w_a > 0 und ist früh ausgeschlossen.
4. Die Stelle: Die Umkehr bei z ≈ 0,4 verlangt ein Glied, das heute fällt. Die Rückstellung kann nur wachsen; der Vakuumteil müsste a₀ mit Maximum bei z = 0,4–0,8 und heutiger Steigung −0,9 … −2,8 haben, geliefert werden −0,57 … +0,23. Der Ansatz „Zeitform aus der Rückstellung, a₀ nach einem Horizont“ trägt damit nicht.
5. Prüfbare Folgen der Zeitform entstehen daher nicht. Im Ansatz selbst (früh verträgliche Horizonte) läge H(z = 0,5–1) zwischen −3 % und +0,6 % gegen ΛCDM, fσ₈ zwischen −2 % und +5 %, D_M(1090) zwischen −0,3 % und +0,9 %.
