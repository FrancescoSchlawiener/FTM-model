# Kern-Nome Δs/r und der glatte Schalter χ

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Abschnitt 2, „s/(e^s − 1) hat die Form der Bose-Besetzung einer Mode“ → die Besetzung ist n_B = 1/(e^s − 1), χ = s·n_B ist die Stufe d = 1 der Bose-Leiter; Ergebnis, „Offen: warum χ die Bose-Form hat … Lesart Anteil zerfallener KAM-Tori … Herleitung steht aus“ → hergeleitet aus der Rückkehr (D·G″ = Γ·G) und der Wirkung (ν = 1 + n_B); KAM-Zerfall allein gibt nur die Wurzel, eine Horizonttemperatur allein nur die Bose-Form.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/kern_chi.py`, `rechnungen/kern_schleife.py` (Abschnitte 2 und 5), Ausgaben `.txt`
> Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Beispielsystem: baryonische Milchstraße, 6·10¹⁰ M☉.

## Ausgangspunkt: die Feldgleichung des Autors

```
r̈ = −∇Φ_bar − χ(a_N) · κ₀ ∫ Δ(τ_M) (e_r⊗e_r) e^{−λ_m(t−t′)} ṙ(t′) dt′

q_Δs = e^{2πiτ_M},   |q_Δs| = Δs/r,   Im τ_M = ln(r/Δs)/2π,   Re τ_M = ω_rad/ω_azim
χ = 0 für a_N > a₀ (gekoppelt, KAM-Tori stabil),   χ = 1 für a_N < a₀ (kritisch, 𝒫 aktiv)
```

Im Fernfeld gibt das a = a_N + χ·(Δs/r)·(3c²/r*). Der Nome des Kerns ist Δs/r, die „Straffung“ des Raumes. Offen war nur der Übergang zwischen den beiden Werten von χ.

## 1. Der Kern-Nome ist in allen Systemen klein

| System | q_Δs = Δs/r | Δ(τ)/q_Δs − 1 |
|---|---|---|
| Sonne, 1 AE | 3,29·10⁻⁹ | −7,9·10⁻⁸ |
| Sonne, bei r* = 7030 AE | 4,68·10⁻¹³ | −1,1·10⁻¹¹ |
| Milchstraße, 8 kpc | 1,20·10⁻⁷ | −2,9·10⁻⁶ |
| Milchstraße, 50 kpc | 1,91·10⁻⁸ | −4,6·10⁻⁷ |

Mit diesem Nome ist Δ(τ_M) auf sechs Stellen gleich q_Δs. Der Kern ist im Fernfeld linear in der Straffung Δs/r.

## 2. Der glatte Schalter

Mit s = √(a_N/a₀) (für die Punktmasse s = r*/r) wird χ zu

```
χ(s) = s/(e^s − 1)
```

- a_N ≫ a₀: χ ≈ s·e^{−s} → 0. Das ist der gekoppelte Bereich des Autors. An der Erde (s = 7031) ist der Zusatz 10⁻³⁰⁵⁵ m/s².
- a_N ≪ a₀: χ ≈ 1 − s/2 → 1. Das ist der kritische Bereich des Autors.

Eingesetzt in die Fernfeldform des Autors folgt exakt die Beschleunigungsrelation:

```
a = a_N + χ(s) · (Δs/r)(3c²/r*) = a_N + χ(s)·√(a_N a₀) = a_N · e^s/(e^s − 1) = a_N/(1 − e^{−s})
```

| r/r* | χ | a/a_N |
|---|---|---|
| 0,1 | 0,00045 | 1,000045 |
| 0,3 | 0,1233 | 1,0370 |
| 1 | 0,5820 | 1,5820 |
| 3 | 0,8426 | 3,5277 |
| 10 | 0,9508 | 10,508 |
| 100 | 0,9950 | 100,50 |

Für eine ausgedehnte Verteilung (Hernquist-Kugel, a = 3 kpc) stimmt a = a_N + χ(s)·√(a_N a₀) von 1 bis 60 kpc auf zehn Stellen mit der Beschleunigungsrelation überein.

**Damit bleiben beide Größen des Autors, wie er sie gesetzt hat.** Der Nome des Kerns ist Δs/r mit Im τ_M = ln(r/Δs)/2π. Der Faktor e^{−s} steht im Schalter χ, nicht im Nome. Die Bose-Besetzung einer Mode mit s als Verhältnis von Energie zu Temperatur ist n_B = 1/(e^s − 1); χ = s/(e^s − 1) = s·n_B ist die Stufe d = 1 der Bose-Leiter s^d·n_B. [korr.]

## 3. Der Realteil folgt aus χ

Re τ_M = ω_rad/ω_azim = κ/Ω = √(2(1 + d ln v/d ln r)). Mit γ_N = d ln a_N/d ln r gilt allgemein

```
(Re τ_M)² = 3 + γ_N · (1 − χ/2)
```

- Newton (χ = 0): (Re τ_M)² = 3 + γ_N, die gewöhnliche Epizykelbeziehung.
- Punktmasse (γ_N = −2): **(Re τ_M)² = 1 + χ**.
- Kontrolle an einer Hernquist-Kugel von 1 bis 60 kpc: Formel und numerische Ableitung der Kurve stimmen auf acht Stellen überein.

Der Schalter χ des Autors, also der Zustand der KAM-Tori, legt damit das Frequenzverhältnis der Bahn fest. χ = 0 gibt Re τ_M = 1 (geschlossene Keplerbahn), χ = 1 gibt Re τ_M = √2 (Rosette).

| r/r* | Re τ_M | Abstand der Perizentren | Präzession je Radialperiode |
|---|---|---|---|
| 0,1 | 1,00023 | 359,92° | −0,08° |
| 0,3 | 1,05986 | 339,67° | −20,33° |
| 1 | 1,25777 | 286,22° | −73,78° |
| 3 | 1,35742 | 265,21° | −94,79° |
| 10 | 1,39672 | 257,75° | −102,25° |
| 100 | 1,41245 | 254,88° | −105,12° |

Die Resonanzradien des Autors, an denen Re τ_M = m/n rational wird, folgen aus 1 + χ(s) = (m/n)².

## 4. Ergebnis

- Die Feldgleichung des Autors mit Kern-Nome Δs/r gibt die Beschleunigungsrelation exakt, wenn χ als χ = s/(e^s−1) mit s = √(a_N/a₀) glatt wird. Beide Grenzwerte von χ bleiben die des Autors.
- Es gibt nur einen Nome, Δs/r. Die Exponentialfunktion gehört zum Schalter.
- Real- und Imaginärteil von τ_M hängen über χ zusammen: (Re τ_M)² = 3 + γ_N(1 − χ/2), für die Punktmasse 1 + χ.
- Die Bose-Form von χ ist hergeleitet: Der Rückweg ist eine abgeschirmte Ausbreitung D·G″ = Γ·G mit Γ = a_N/c und τ_D = c/a₀, das gibt q = e^{−s}; die Strukturwirkung gibt ν = 1 + n_B(s), also χ = s·n_B. Der Zerfall der KAM-Tori allein gibt nur die Wurzel, eine Horizonttemperatur allein nur die Bose-Form ohne Wurzel. [korr.]
