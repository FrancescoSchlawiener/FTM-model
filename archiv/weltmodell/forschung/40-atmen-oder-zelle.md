# Atmende Schwelle oder feste Zelle

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Die Notiz rechnet mit a₀ = cH₀/2π = 1,042·10⁻¹⁰ m/s² (die Tabelle in Abschnitt 3 mit 1,2·10⁻¹⁰). „Feste Zelle“ a₀ = c²(1 + z)/L₀ und „Soll die Zelle bestehen, muss sie selbst mit dem Horizont atmen“ → im Modell ist c²/a₀ = ℓ₀ = 2Λ keine Raumlänge; die Zelle atmet mit (L = L̄·a) und ist über die innere Faltung gebunden, a₀ = 3c⁴/(2Gρ_mL³) ist fest; beide Lesarten der Notiz tragen nicht. Genzel 2017/2020 [G] stützt die feste Schwelle. Tabelle Abschnitt 3: a₀ vermerkt (f(z = 0) = 0,116 / 0,340 bei 1,2·10⁻¹⁰; mit 1,042·10⁻¹⁰ 0,099 / 0,315). Umbenannt: L = c²/a₀ → ℓ₀ (L₀ → ℓ₀(0)); Bremsparameter q₀ → q_B.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/atmen_oder_zelle.py` (Ausgabe `atmen_oder_zelle.txt`)
> Kosmologie: Ω_m = 0,315, Ω_r = 9,1·10⁻⁵, H₀ = 67,4 km/s/Mpc. Schwelle heute a₀ = cH₀/2π = 1,042·10⁻¹⁰ m/s².

## Die zwei Lesarten

| Lesart | Schwelle | Länge ℓ₀ = c²/a₀ [korr.] |
|---|---|---|
| Atmen | a₀(z) = cH(z)/2π | ℓ₀ = 2πc/H(z), der Umfang des Hubble-Horizonts |
| feste Zelle | a₀(z) = c²(1+z)/ℓ₀(0) [korr.] | ℓ₀ ∝ Skalenfaktor, fester mitbewegter Umfang |

[korr.] Beide Lesarten tragen nicht. c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge. Atmete a₀ mit H, bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); atmete Λ mit den Schleifen (a₀ ∝ 1 + z), wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell atmet die Zelle mit (L = L̄·a) und ist über die innere Faltung an die Schwelle gebunden: a₀ = 3c⁴/(2Gρ_mL³) ist fest (seit z ≈ 10, bei der Rekombination 0,76·a₀); bei a₀ = 1,2·10⁻¹⁰ ist L = 26,8 Gpc, ohne Kreispaare ist L ≥ 0,985·d_LSS (COMPACT 2023), also a₀ ≤ 1,13·10⁻¹⁰. Die Rechnungen unten gelten nur unter den beiden Lesarten.

Heute sind beide gleich, das ist die Festlegung ℓ₀(0) = 2πc/H₀. [korr.] Der mitbewegte Umfang des Horizonts, 2πc/(aH), ist aber nicht fest:

| z | 0 | 0,5 | 1 | 2 | 5 | 1100 |
|---|---|---|---|---|---|---|
| 2πc/(aH) relativ zu heute | 1 | 1,134 | 1,117 | 0,989 | 0,723 | 0,047 |

Eine mitbewegte Zelle fällt also nur zu einer Epoche mit dem Horizont zusammen. Beides zugleich gilt nur heute, bei z = 1,917 (Gleichstand H(z)/H₀ = 1+z) und nirgends sonst.

## 1. Rotationskurven und Tully-Fisher

Im tiefen Bereich gilt g_obs/g_bar ∝ √a₀ und v⁴ = G M a₀.

| z | a₀/a₀(0) Atmen | a₀/a₀(0) Zelle | √(Atmen/Zelle) | Δlog v bei festem M | Δlog M bei festem v |
|---|---|---|---|---|---|
| 0,5 | 1,322 | 1,500 | 0,939 | −0,014 dex | +0,055 dex |
| 1 | 1,791 | 2,000 | 0,946 | −0,012 dex | +0,048 dex |
| 2 | 3,033 | 3,000 | 1,005 | +0,001 dex | −0,005 dex |
| 3 | 4,568 | 4,000 | 1,069 | +0,014 dex | −0,058 dex |
| 4 | 6,334 | 5,000 | 1,126 | +0,026 dex | −0,103 dex |
| 6 | 10,44 | 7,000 | 1,221 | +0,043 dex | −0,174 dex |
| 10 | 20,53 | 11,00 | 1,366 | +0,068 dex | −0,271 dex |

- Bei z = 0,5–1 liegt die Zelle 5–6 % über dem Atmen, im Nullpunkt der Tully-Fisher-Relation sind das 0,05 dex in der Masse. Das ist kleiner als die heutige Streuung der Massenbestimmung bei diesen Rotverschiebungen.
- Die Trennschärfe wächst oberhalb von z ≈ 3. Bei z = 6 liegt das Atmen um 22 % in g_obs/g_bar über der Zelle, bei z = 10 um 37 %. Rotierende Scheiben bei z = 4–7 aus ALMA ([CII]) und JWST sind damit der schärfste Test.

## 2. Heutige Drift der Schwelle

| Lesart | d ln a₀/dt | d ln v_flach/dt |
|---|---|---|
| Atmen: Ḣ/H = −(1+q_B)H₀ = −0,473 H₀ [korr.] | −3,26 %/Gyr | −0,81 %/Gyr |
| feste Zelle: −H₀ | −6,89 %/Gyr | −1,72 %/Gyr |

Direkt messbar ist das nicht. Bei fester baryonischer Masse folgt v_flach = (G M a₀)^{1/4} der Schwelle. Mit der heutigen Rate hochgerechnet sinkt v_flach über 14 Gyr um rund 11 % (Atmen) bzw. 24 % (Zelle).

## 3. Massereiche Scheiben bei z = 1–2

Beide Lesarten lassen a₀ mit z wachsen. Für dieselbe Scheibe wird der Zusatzanteil der Beschleunigung, f = 1 − 1/ν, zu früheren Zeiten größer (M_b = 10¹¹ M☉, Punktmassen-Näherung):

| R | z = 0 | z = 1 (Atmen / Zelle) | z = 2 (Atmen / Zelle) |
|---|---|---|---|
| 5 kpc | 0,116 | 0,200 / 0,218 | 0,290 / 0,288 |
| 10 kpc | 0,340 | 0,447 / 0,467 | 0,539 / 0,537 |

[korr.] Die Tabelle steht auf a₀(0) = 1,2·10⁻¹⁰ m/s², nicht auf dem Wert der Kopfzeile; mit a₀(0) = 1,042·10⁻¹⁰ gibt der Kern bei z = 0 f = 0,099 (5 kpc) und 0,315 (10 kpc).

Gemessen wurde an massereichen sternbildenden Scheiben bei z ≈ 0,9–2,4 das Gegenteil: außen fallende Rotationskurven und innerhalb des Halblichtradius kleine Zusatzanteile, f ≲ 0,2 (Genzel et al. 2017, 2020 [G]). Das spricht gegen jedes a₀, das mit z wächst, also gegen beide Lesarten gleichermaßen. Druckstützung des Gases und die Massenmodelle der Scheiben sind dabei die Unsicherheiten. Ein a₀, das zu früheren Zeiten nicht wächst, würde auch den DESI-Bedarf (a₀(1)/a₀(0) ≈ 1,03) erfüllen. Die Messung stützt damit die feste Schwelle des Modells (a₀ seit z ≈ 10 fest). [korr.]

## Ergebnis

- Heute sind beide Lesarten gleich. Bei z = 1,92 sind sie noch einmal gleich, dazwischen weichen sie um höchstens 6 % in √a₀ ab.
- Zwischen den beiden Lesarten trennt erst z ≳ 4: 13 % (z = 4) bis 37 % (z = 10) in g_obs/g_bar.
- Massereiche Scheiben bei z = 1–2 zeigen kleinere Zusatzanteile, als beide Lesarten vorhersagen [G]. Das ist ein offener Widerspruch zu jedem a₀, das mit z wächst, und stützt die feste Schwelle des Modells. [korr.]
- Setzt man die Zelle gleich c²/a₀, ist eine Zelle mit fester Topologie mit dem Atmen nicht vereinbar. [korr.] Das trägt nicht: c²/a₀ = ℓ₀ = 2Λ ist keine Raumlänge. Im Modell atmet die Zelle mit (L = L̄·a) und ist über die innere Faltung gebunden, a₀ = 3c⁴/(2Gρ_mL³) ist fest; mit dem Horizont atmen muss sie nicht.
