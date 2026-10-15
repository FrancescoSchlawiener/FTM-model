# Halbgerade als Spitze

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Knick einer festen inneren Länge „bei r ∝ r* ∝ √M“ → massenunabhängig bei r ≈ 0,657·Λ (mit dem Kern gleich für 10⁸, 6·10¹⁰ und 10¹⁴ M☉; heute 8,0–9,2 Gpc); Periode dieser Länge in Im τ „abhängig von System und Radius“ → nur vom Radius; q als Nome des Bahntorus gekennzeichnet; innere Richtung w → ϖ, ihre Periode Λ_w → P_ϖ; dimensionslose Periode in Im τ Λ → P (q_Λ → q_P); Bessel-Index ν → m_B; Längenmaß ℓ → ℓ_τ.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/halbgerade_spitze.py` (Ausgabe `halbgerade_spitze.txt`)
> τ = x + iy, q = e^{2πiτ} (Nome des Bahntorus: bei der ersten Wahl, Punkt 1, für die Punktmasse der Rückkehranteil, bei der zweiten, Punkt 3, der Kern-Nome q_Δs [korr.]), s = √(a_N/a₀). Frage: Ist die offene innere Richtung (D·G'' = Γ·G, G(0) = 1, G beschränkt ⇒ q = e^{−s}) die Richtung Im τ zur Spitze des Modulraums?

## 1. Mode λ = 0 und Halbgerade

```
y²(∂_x² + ∂_y²)f = λf,   f = e^{2πinx}g(y)   ⇒   y²(g'' − 4π²n²g) = λg
λ = 0:   g'' = (2πn)²g,   g = e^{±2πny},   beschränkt für y → ∞:   g = e^{−2π|n|y} = |qⁿ|
ϖ = ℓ_τ·Im τ,  ℓ_τ = 2π√(D/Γ):   D·G'' = Γ·G  ist die Mode n = 1,   L√(Γ/D) = 2π·Im τ_M,   G(L) = |q(τ_M)|
```

- **G(0) = 1** ist von selbst erfüllt: Bei y = 0 ist |qⁿ| = 1 für jedes n. Der Ort ist der Einheitskreis |q| = 1 (Torus mit verschwindendem innerem Zyklus), die Spitze ist q = 0. Für λ = 0 zählt nur die konforme Struktur: Der Halbzylinder ist die punktierte Scheibe 0 < |q| < 1, der Ort ihr Rand, die Spitze ihr Mittelpunkt.
- **Wahl Im τ_M = s/2π:** ℓ_τ = 2πr, ϖ = 2πr·Im τ, G(L) = e^{−s} = q. Wegen Im τ = (ϖ/L)·s/2π sind die Halbgerade in ϖ bei festem r und der Verlauf in a_N bei ϖ = L dieselbe Gerade: Ort ↔ a_N → 0 (q → 1), Spitze ↔ a_N → ∞ (Newton, q → 0).
- **Wahl Im τ_M = ln(r/Δs)/2π:** G(L) = Δs/r, der Kern-Nome, nicht der Rückkehranteil. Dazu müsste L√(Γ/D) = ln(r/Δs) statt s sein (Milchstraße 6·10¹⁰ M☉ bei 8 kpc: 15,9 statt 1,04). Ort ↔ r = Δs, Spitze ↔ r → ∞; e^{−s} steckt hier in χ, nicht in der Mode.
- **Modulgruppe:** Der Stabilisator von ∞ in SL(2,ℤ) enthält nur τ ↦ τ + b. Re τ ist geschlossen, Im τ hat keine Identifikation; die Achse Re τ = 0 wird nur durch S (y ↦ 1/y) bei τ = i gefaltet.
- **Wo die Punkte liegen:** Der Spitzenbereich von ℍ/SL(2,ℤ) ist y > 1, bei der ersten Wahl s > 2π, a_N > 4π²a₀ = 39,5 a₀. Galaktische Punkte liegen darunter und sind unter S Punkten tief in der Spitze gleich (Im(−1/τ) = 2π/s). Auf dem Modulraum fallen Ort und Newton-Grenze also in dasselbe Ende; |q| ist dort nicht erklärt (τ = i/2 und 2i geben ν = 1,045 und 1,0000035). ν lebt auf ℍ/⟨T⟩, dem markierten Torus mit unterschiedenem Umlauf und Rückkehr, Im τ ∈ (0, ∞). Bei der zweiten Wahl liegen alle realen Systeme im Spitzenbereich (Erde 3,11, Milchstraße 2,54; y > 1 ⇔ r > 535 Δs).

## 2. Umlaufsinn, Gewicht, Maass-Formen

- **Umlaufsinn:** Beschränkt an der Spitze mit λ = 0 sind 1, qⁿ und q̄ⁿ; q̄ⁿ ist ebenso harmonisch und beschränkt (∂_τ q̄ⁿ = 0). Die Spitze wählt keinen Sinn; Holomorphie ist eine zusätzliche Bedingung.
- **Gewicht:** G(0) = 1 für jede gezählte Mode heißt Koeffizient cₙ = 1. Gezählt n ≥ 0: Σqⁿ = 1/(1 − q) = 1 + n_B (Szegő-Kern der Scheibe, δ am Ort auf n ≥ 0 projiziert). Gezählt n ∈ ℤ: Poisson-Kern (1 − |q|²)/|1 − q|², auf der Achse coth(s/2) = 1 + 2n_B. Bergman 1/(1 − q)² (cₙ = n + 1) und coth s (nur gerade n) verletzen G(0) = 1 je Mode.

| a_N/a₀ | 1 | 0,1 | 0,01 | s → 0 |
|---|---|---|---|---|
| Szegő 1/(1 − q) | 1,582 | 3,689 | 10,51 | 1/s, v⁴ = GMa₀ |
| Poisson, n ∈ ℤ | 2,164 | 6,377 | 20,02 | 2/s, v⁴ = 4GMa₀ |
| Bergman 1/(1 − q)² | 2,503 | 13,61 | 110,4 | 1/s², a → a₀ |

- **Phase:** Die Funktion der Spitze hängt von Re τ ab. 1 + n_B verlangt Re τ = 0 (innerer Zyklus senkrecht zur Bahn); mit Re τ = √2 − 1 (κ/Ω = √2) fällt Re[1/(1 − q)] im tiefen Bereich auf ½ (0,503 bei s = 0,01).

```
λ ≠ 0:   g = √y·K_{m_B}(2πny),   m_B = √(¼ + λ),   am Ort g ∝ y^{½−m_B}
auf das führende Glied normiert:   q_λ(s) = (2/Γ(m_B))·(s/2)^{m_B}·K_{m_B}(s),   1 − q_λ ∝ s^{2m_B} (m_B < 1),   a ∝ a_N^{1−m_B}·a₀^{m_B}
```

- Einen endlichen Wert G(0) ≠ 0 hat nur λ = 0. Im tiefen Bereich ist d ln v/d ln r = m_B − ½ (für m_B ≥ 1: +½, Boden a → konst.): λ = −0,2 / −0,0475 / 0 / 0,0525 / 0,5 gibt −0,276 / −0,050 / 0 / +0,050 / +0,346. Flach auf ±0,05 heißt −0,0475 ≤ λ ≤ 0,0525.
- **Spektrum von SL(2,ℤ)** (Hejhal-Verfahren, Hecke-Probe a₂a₃ = a₆): Spitzenformen bei t = 9,5337 (ungerade), 12,1730, 13,7798 (gerade), 14,3585, also λ = −(¼ + t²) = −91,14, −148,4, …; dazu das kontinuierliche Spektrum λ ≤ −¼ und λ = 0 (Konstante). Im Tully-Fisher-Band liegt nur λ = 0.
- Für λ ≤ −¼ geht g am Ort wie √y gegen 0 und schwingt für y < t/2π (t₁: a_N < 91 a₀ bei der ersten Wahl). q wechselt bei s = 6,08; 4,13; 2,90; … das Vorzeichen (Verhältnis → e^{π/t₁} = 1,390) und fällt wie √s: a → a_N, keine flache Kurve.

## 3. Geschlossene Schleife in Im τ

```
y ~ y + P:   Gₙ(y) = cosh(2πn(y − P/2)) / cosh(πnP)          erste Wahl:   q_P(s) = cosh(s − πP) / cosh(πP)
tief:   1 − q ≈ s·tanh(πP)   ⇒   v⁴ = GMa₀ / tanh²(πP);        q = 1 bei s = 2πP
```

- Nur die flache Identifikation ist möglich (λ = 0). Eine hyperbolische Streckung erzeugt mit τ ↦ τ + 1 Verschiebungen um 1/kᵐ; die Gruppe ist nicht diskret.
- Gₙ ist gerade um den Ort: beide Sinne. Auf dem Torus gibt es keine holomorphe Funktion außer Konstanten und keine mit nur einem einfachen Pol am Ort (Residuensumme 0). Die einseitige Reihe Σqⁿ braucht das offene Ende.
- 1/cosh s entsteht nur, wenn der Ort gegenüber liegt (P = 2·Im τ_M). Bei fester Periode bleibt der tiefe Bereich Tully-Fisher (P = 1: v⁴ um 0,75 % höher); die Abweichung liegt innen, im Newton-Bereich: 10 % ab s = 3,87 / 10,2 / 29,0 für P = 1 / 2 / 5, Pol bei s = 2πP (die Modensumme ΣGₙ gibt 3,91 / 10,18 / 29,03). Newton an der Sonnenoberfläche (s = 1,5·10⁶) verlangt P > 2,4·10⁵; dann ist die Schleife in Galaxien unsichtbar.
- Zweite Wahl: Der Nome wächst außen wieder an, +10 % ab Im τ = P/2 − 0,183, also bei r ∝ Δs ∝ M. Flach bis r_max verlangt P ≥ 2·Im τ(r_max) + 0,37 (Milchstraße bis 100 kpc: P ≥ 6,24).
- Eine feste innere Länge in Metern ist keine Schleife in Im τ: Ihre Periode wäre P_ϖ/(2πr), abhängig vom Radius, nicht von der Masse. [korr.]

## 4. Prüfbare Folgen

- Offene Spitze mit λ = 0: genau a = a_N/(1 − e^{−s}), keine neue Form.
- Beimischung λ ≠ 0: im tiefen Bereich eine massenunabhängige Restneigung d ln v/d ln r = √(¼ + λ) − ½.
- Knick der flachen Kurve: Schleife in Im τ bei der ersten Wahl innen bei festem a_N ≈ (2πP)²a₀ (durch Newton ausgeschlossen), bei der zweiten außen bei r ∝ M; feste innere Länge (Abschluss der Rückkehr bei Λ, gleichwertig einer Schleife P_ϖ = 2Λ) massenunabhängig bei r ≈ 0,657·Λ, mit dem Kern gleich für 10⁸, 6·10¹⁰ und 10¹⁴ M☉, heute 8,0–9,2 Gpc. [korr.] Ein gefundener Knick trennt die Fassungen über seine Massenskalierung.

## Ergebnis

Es trägt: Die Mode λ = 0 auf dem Halbzylinder ist genau D·G'' = Γ·G; G(0) = 1 ist |qⁿ| = 1 auf dem Einheitskreis, das Gewicht 1 je Mode folgt daraus. λ = 0 ist der einzige Wert mit endlichem G(0), mit Tully-Fisher und im Spektrum von SL(2,ℤ) innerhalb des Bandes. Im τ hat in der Modulgruppe keine Identifikation; eine Schleife fester Periode verschont Tully-Fisher und zerstört die Newton-Grenze: Bei der ersten Wahl ist die Halbgerade gleichwertig mit der Newton-Grenze q → 0 an der Spitze.

Es trägt nicht: Die Spitze wählt keinen Umlaufsinn (q̄ⁿ ist ebenso beschränkt, Poisson 1 + 2n_B ebenso zulässig). Der Spitzenbereich von ℍ/SL(2,ℤ) beginnt bei Im τ = 1; galaktische Punkte liegen darunter, und auf dem Modulraum sind Ort und Newton-Grenze dasselbe Ende. Getragen wird die Halbgerade von ℍ/⟨T⟩ (markierter Torus, Re τ = 0), nicht vom Modulraum. Bei der zweiten Wahl trägt die Spitze den Nome Δs/r, nicht e^{−s}.
