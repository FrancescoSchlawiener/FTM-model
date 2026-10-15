# Massenverhältnis und eine Zahl für α und α_G

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: nur Bezeichnungen, Nome der Kopplungen q → q̃ (q ist der Rückkehranteil); inhaltlich keine Korrektur nötig, [G]-Marken vollständig.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/massenverhaeltnis.py` (Ausgabe `massenverhaeltnis.txt`)
> Ziele [G] (CODATA, aus dem Gedächtnis): m_p/m_e = 1836,15267343 (±6·10⁻¹¹ relativ), 1/α = 137,035999177 (±1,5·10⁻¹⁰), α_G = 5,906·10⁻³⁹.

## 1 Zufallsprobe

Klassen einfacher Ausdrücke aus Zahlen, die im Modell vorkommen:

- **A** 2^a·3^b·5^c·π^d, Exponenten in ℤ/n mit |·| ≤ 12/n, n ≤ 4 (1,5 Mio. Werte; enthält 3, 6, 12, 24, 45, 1728, Wurzeln);
- **B** Σ_{i≤5} c_iπ^i mit |c_i| ≤ 4 (0,27 Mio.);
- **C** (p/q)·π^k mit p, q ≤ 64, |k| ≤ 6 (33 000).

δ = |x/Ziel − 1|; erwartet ist N(δ) = Dichte·2δ; p ist der Anteil von 4000 Zufallszielen (Faktor 0,67–1,5 um das Ziel), deren bester Treffer gleich gut oder besser ist.

| Ziel | Klasse | Werte je Einheit ln x | bester δ | N(10⁻⁴) gefunden / erwartet | N(10⁻⁵) | p |
|---|---|---|---|---|---|---|
| m_p/m_e | A | 39 910 | 9,7·10⁻⁶ | 9 / 8,0 | 1 / 0,8 | 0,59 |
| | B | 26 825 | 2,3·10⁻² | 0 / 5,4 | 0 / 0,5 | 0,49 |
| | C | 978 | 1,9·10⁻⁵ (6π⁵) | 1 / 0,2 | 0 / 0,02 | 0,04 |
| 1/α | A | 59 715 | 1,2·10⁻⁵ | 12 / 11,9 | 0 / 1,2 | 0,85 |
| | B | 26 390 | 2,2·10⁻⁶ (4π³ + π² + π) | 6 / 5,3 | 2 / 0,5 | 0,11 |
| | C | 2 075 | 4,0·10⁻⁴ | 0 / 0,4 | 0 / 0,04 | 0,87 |

- Die gefundenen Anzahlen folgen der Zufallserwartung. Das kleinste p (0,04, 6π⁵ in C) tritt bei sechs Prüfungen mit Wahrscheinlichkeit 0,21 mindestens einmal auf. Beide bekannten Ausdrücke liegen zudem 3·10⁵ bzw. 1,5·10⁴ Messfehler neben dem Wert.
- Auf Messgenauigkeit erwartet jede Klasse 10⁻⁷–10⁻⁵ Zufallstreffer; gefunden wird keiner.
- **Schwelle:** Eine Beziehung zählt nur, wenn sie aus einer Gleichung des Modells folgt und weit unter δ ≈ 10⁻⁵ trifft.

## 2 Strukturelle Wege des Modells

**Längenverhältnisse.** α = r_e/λ̄_e, α_G = 3Δs_p/λ̄_p, m_p/m_e = λ̄_e/λ̄_p. Symbolisch gilt Δs·λ̄ = Għ/(3c³) = ℓ_P²/3 und a_J·Δs = ℓ_P²/6 (Spin ½) für jede Masse: m kürzt sich. Die Beziehungen verbinden eine Masse mit ihrer eigenen Länge; zwischen α, α_G und m_p/m_e geben sie keine Gleichung (drei Zahlen, keine Beziehung).

**Carter-Kette.** α_G = 𝒞·α¹²·(m_e/m_p)⁴ mit 𝒞 = 2,944 aus den Messwerten (Kramers-Form mit Exponent 5: 𝒞′ = 5406). Mit 𝒞 = 1 folgte m_p/m_e = 1402 (−24 %) bzw. 1/α = 125,2; Kramers mit 𝒞′ = 1 gäbe 329. Wegen δμ/μ = ¼·δ𝒞/𝒞 müsste 𝒞 für m_p/m_e auf 10⁻¹⁰ selbst auf 4·10⁻¹⁰ bekannt sein; 𝒞 enthält Saha- und Gamow-Logarithmen und die mittlere Opazität. Die Kette legt Exponenten fest, keinen Wert.

**Nome-Darstellung** |q̃| = e^{−2π Im τ}:

| Größe | Im τ |
|---|---|
| α_G = e^{−2π Im τ_G} | 14,0096 (bei α_G = \|q̃\|²: 7,0048) |
| α = e^{−2π Im τ_α} | 0,7831 |
| Kern-Nome Δs/λ̄ an der eigenen Quantenlänge: Proton / Elektron | 14,1844 / 16,5767 |

- m_p/m_e = e^{π(Im τ_e − Im τ_p)}: Umschreibung, keine Bestimmung.
- Carter in Nome-Form: Im τ_G = 12·Im τ_α + (2/π)·ln(m_p/m_e) − ln 𝒞/2π = 9,3970 + 4,7845 − 0,1719.
- Eine Zahl Im τ_X für α und α_G verlangt neben α_G = e^{−2π Im τ_X} eine zweite Gleichung α = f(Im τ_X). Das Modell hat keine; der Weg über den Lauf der Kopplungen braucht Teilcheninhalt, Vereinigungsskala und die QCD-Zahl C, die keine Größen des Modells sind. Selbst mit ihr bliebe m_p/m_e an 𝒞 gebunden: μ = e^{(π/2)(Im τ_G − 12 Im τ_α)}·𝒞^{1/4}.
- Nähen im Zufallsbereich: Im τ_G liegt relativ 6,8·10⁻⁴ neben 14 (Wahrscheinlichkeit 0,019 für einen Wert, 0,07 bei vier geprüften Im-τ-Werten); ln α_G/ln α = 17,890 (bester Ausdruck der Klasse A, p = 0,06); (1/α)/Im τ_G = 9,782 (p = 0,37).

## Ergebnis

- **Nicht hergeleitet.** Keine Gleichung des Modells bestimmt m_p/m_e, keine verbindet α und α_G über eine gemeinsame Zahl Im τ_X. Die Längenbeziehungen gelten für jede Masse, die Carter-Kette hat einen offenen Vorfaktor, die Nome-Form ist eine Umschreibung.
- Alle Zahlennähen (6π⁵, 4π³ + π² + π, Im τ_G ≈ 14) liegen im Bereich zufälliger Treffer und folgen aus keiner Gleichung.
- **Gegen das Modell geprüft:** kein Widerspruch. „m_p/m_e offen“, „eine Zahl Im τ_X offen“ und „Zahlenfaktor ≈ 2,9 in Carters Beziehung nicht hergeleitet“ bleiben so stehen; 𝒞 = 2,944 stimmt mit ≈ 2,9 überein.
