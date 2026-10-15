# Zeit der Rückkehr

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „dieselbe offene Größe“ und „ℓ₀/c = 2π/H bleibt Setzung“ → τ_D = 2Λ/c (Kac), Λ die ganze innere Faltung der Zelle, ≈ 2π/H₀ nur heute; Abschnitt 4 (τ_D = 2π/H(z), a₀ ∝ H) als nicht tragende Annahme gekennzeichnet; H₀ = 73 mit [G]; innere Richtung w → ϖ; allgemeine Wirkungsform Q → F; H_Λ → H_DE, Ω_Λ → Ω_DE.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/zeit_rueckkehr.py` (Ausgabe `zeit_rueckkehr.txt`)
> a₀ = 1,2·10⁻¹⁰ m/s², ℓ₀ = c²/a₀, n_B(x) = 1/(eˣ − 1). Beispiel: Milchstraße 6·10¹⁰ M☉, r* = 8,35 kpc.

q = e^{−s} folgt aus D·G'' = Γ·G mit Γ = a_N/c und τ_D = L²/D = c/a₀. Gefragt ist, woher τ_D kommt und ob es die Periode β = 2π/H eines Horizontgleichgewichts ist.

## 1. Eine Größe, nicht zwei

```
Γ = a_N/c,  √(Γ/D) = 1/r   ⇒   D = a_N·r²/c        (Punktmasse: D = GM/c, für alle r dieselbe)
L² = D·τ_D   ⇒   L/r = √(a_N·τ_D/c) = √(a_N/a₀)   ⇔   τ_D = c/a₀;     Punktmasse: L = r* ⇔ τ_D = c/a₀
```

Sonne D = 4,43·10¹¹ m²/s, L = 7031 AE; Milchstraße D = 2,66·10²² m²/s, L = 8,35 kpc; beide τ_D = 79,2 Mrd. Jahre. Innere Länge und Rückkehrzeit sind dieselbe Größe; festgelegt ist sie über Kac: τ_D = 2Λ/c ist die mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden, Λ die ganze innere Faltung der Zelle (heute Λ ≈ πc/H₀, also τ_D ≈ 2π/H₀ nur heute). [korr.] Für ausgedehnte Massen gilt das mit D = a_N r²/c an jedem Ort.

## 2. Q̃ als thermische Größe

```
Q̃ = −4∫_s^∞ x³ n_B(x) dx = −∫_s^∞ n_B(x) d(x⁴),        ∫_s^∞ x³e^{−nx}dx = e^{−ns}(s³/n + 3s²/n² + 6s/n³ + 6/n⁴)
∂ℒ_𝒮/∂Y = 2·n_B(s)   (s⁴ = ℓ₀²Y)    ⇒    ν = 1 + n_B(s),     q = n_B/(1 + n_B) = e^{−s}
```

Die Summe über n gibt die Polylog-Form; numerisch gleich auf 10⁻¹⁷ (s = 10⁻³ … 30). Die Faltungsgleichung trägt den Bose-Faktor 1 + n; q ist das Verhältnis n/(1 + n), das in einem Gleichgewicht der Periode β gleich e^{−βω} ist. Festgelegt ist damit nur das Produkt βω = s.

Vergleich mit Einschleifen-Größen (je Freiheitsgrad, Einheit (kT)⁴/(2π²(ħc)³), masselos, drei Raumrichtungen). Aus der Feldgleichung folgt für jede Form F: ν − 1 = dF/d(s⁴); Tully-Fisher verlangt s·ν → 1:

| Form | F bei s = 0,3 / 1 / 3 | s·ν bei s = 0,01 |
|---|---|---|
| Q̃ | −25,94 / −25,08 / −15,77 | 1,005 |
| freie Energie, Moden unter s entfernt: −12[s²Li₂ + 2sLi₃ + 2Li₄] | −25,80 / −23,24 / −10,25 | 13,8 |
| freie Energie mit Masse s/β | −25,23 / −20,36 / −6,69 | 489 |
| freie Energie, vier Raumrichtungen, geschnitten | – | 0,02 |

- Die freie Energie eines geschnittenen Feldes ist Q̃ in keiner Dimension. Ihre Ableitung nach dem Schnitt ist ln(1 − e^{−s}), nicht n_B(s): ν − 1 = 3Li₁(e^{−s})/s wächst logarithmisch. Mit Masse folgt a → (π²/2)·a₀, mit vier Richtungen a → a_N.
- Exakt gilt Q̃ = −4·U₃, U₃ = ∫_s^∞x³n_B dx die Energie der Moden oberhalb des Schnitts. Gleichwertig ist Q̃ = 12·∫_0^∞x² ln(1 − e^{−max(x,s)})dx: die freie Energie, wenn die Moden unter dem Schnitt nicht entfernt, sondern auf ω_s gehalten werden. Vom scharfen Schnitt trennt sie der Randterm 4s³·ln(1 − e^{−s}).
- x³ ist die Zustandsdichte x² (drei Richtungen, masselos) mal die Energie x. Die 4 ist der Exponent in s⁴ = ℓ₀²Y: Q̃ = −∫n_B d(x⁴).
- Normierung: Q̃(0) = −4π⁴/15 gibt Λ_s = 4π⁴/(15ℓ₀²) mit dem Vorfaktor c⁴/(16πG)·2/ℓ₀². Ein Freiheitsgrad bei kT = ħa₀/c trägt ħc·Q̃/(24π²ℓ₀⁴); das Verhältnis ist 3πℓ₀²/ℓ_P² = 2·10¹²⁴. Die Wirkung legt nur s = βω und N·T⁴ fest, β allein nicht.

| Lesart von x | x | Zeit |
|---|---|---|
| Windungen in ϖ (Diffusion) | κL, κ = 1/r, L = r* | Ausbreitungszeit L²/D = c/a₀, für jedes M gleich |
| Bildsumme in imaginärer Zeit mit der Rate der Gleichung, ω = Γ | Γβ = s² bei β = τ_D | ν = 1 + n_B(s²), a → a_N + a₀, v² ∝ r |
| Welle in ϖ mit c | κL | β = L/c: Sonne 40,6 d, Milchstraße 27 000 a |
| Horizontgleichgewicht mit Lichtmoden, kT = ħa₀/c | k·ℓ₀ | Kreisbahnmode bei r*: x = c²/v∞² = 2,9·10⁶ |

β = τ_D gibt n_B(s) nur mit ω = s/τ_D = √(Γ/τ_D) = √(a_N·a₀)/c; diese Rate kommt in D·G'' = Γ·G nicht vor. Als Welle in ϖ verlangt β = τ_D die Geschwindigkeit L/τ_D = v∞²/c (Sonne 4·10⁻⁴ m/s, Milchstraße 103 m/s), also eine von M abhängige. Für alle Massen gleich ist τ_D nur als Ausbreitungszeit L²/D, weil L² und D beide ∝ M sind.

x ist also κL; die Bose-Form kommt aus den Windungen in ϖ, nicht aus einer periodischen Zeit.

## 3. Kandidaten für τ_D

a₀ = c/τ_D; flach, Ω_m = 0,315 für beide H₀ angenommen, ohne Rückwirkung des Vakuumteils auf H.

| τ_D | H₀ = 67,4: τ_D [Mrd. J.] · a₀ [10⁻¹⁰] · Abstand | H₀ = 73 [G]: τ_D · a₀ · Abstand |
|---|---|---|
| ℓ₀/c (Wirkung, a₀ gemessen) | 79,2 · 1,200 · 0 | 79,2 · 1,200 · 0 |
| 2π/H₀ (Hubble-Radius) | 91,2 · 1,042 · −13 % | 84,2 · 1,129 · −6 % |
| 2π/H_DE (de Sitter, H_DE = H₀√Ω_DE) | 110,1 · 0,863 · −28 % | 101,7 · 0,934 · −22 % |
| 2π·R_E/c (Ereignishorizont, R_E = 1,150 c/H₀) | 104,8 · 0,907 · −24 % | 96,8 · 0,982 · −18 % |
| 1/H₀ | 14,5 · 6,55 · +446 % | 13,4 · 7,09 · +491 % |
| Weltalter t₀ = 0,951/H₀ | 13,8 · 6,89 · +474 % | 12,7 · 7,46 · +522 % |

Ein exaktes Gleichgewicht der Periode 2π/H gibt es nur in de Sitter; dort liegt a₀ am weitesten unter dem Messwert. 1/H₀ und t₀ ohne 2π liegen um das Fünf- bis Sechsfache darüber. Für a₀ = 1,2·10⁻¹⁰ ist τ_D = 0,868·2π/H₀ (67,4) bzw. 0,941·2π/H₀ (73) nötig, als 2π/H₀ allein H₀ = 77,6 km/s/Mpc.

## 4. Folge von τ_D = β(z)

[korr.] Diese Rechnung setzt τ_D = 2π/H(z), also a₀ ∝ H, voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute.

```
τ_D = 2π/H(z)   ⇒   a₀(z)/a₀(0) = H(z)/H₀,     v∞⁴ = G·M_b·a₀(z),     r* ∝ a₀^{−1/2}
```

| z | H/H₀ | r*(z)/r*(0) | v∞ bei festem M_b | log M_b bei festem v∞ |
|---|---|---|---|---|
| 0,5 | 1,322 | 0,870 | +7,2 % | −0,12 dex |
| 1 | 1,790 | 0,747 | +15,7 % | −0,25 dex |
| 2 | 3,032 | 0,574 | +32,0 % | −0,48 dex |

H(z) ist messbar einzusetzen; die Zahlen gelten für den obigen ΛCDM-Verlauf. Mit τ_D = ℓ₀/c oder 2π/H_DE bleibt a₀(z) = a₀(0).

## Ergebnis

| Teil | trägt | Stelle |
|---|---|---|
| L = r* ⇔ τ_D = c/a₀ | ja, eine Größe | D = a_N r²/c aus Γ und der Abschirmung 1/r |
| Q̃ = −4∫_s^∞x³n_B dx | ja, exakt | gliedweise Summe über n |
| Bose-Faktor in der Gleichung | ja | ν = 1 + n_B(s), q = n_B/(1 + n_B) |
| Q̃ als freie Energie eines geschnittenen Feldes | nein | Ableitung ln(1 − e^{−s}) statt n_B; kein Tully-Fisher |
| Q̃ als Energie oberhalb des Schnitts | ja, der Form nach | x³ = x²·x; 4 aus s⁴ = ℓ₀²Y |
| Normierung als Einschleifen-Größe | nein | Faktor 3πℓ₀²/ℓ_P² = 2·10¹²⁴ |
| τ_D = β aus der Wirkung | nein | nur βω = s festgelegt; KMS mit ω = Γ gibt n_B(s²); β = τ_D verlangt ω = √(a_N a₀)/c |
| τ_D = ℓ₀/c | ja | einzige Länge der Wirkung; ihr Wert ist die Konstante a₀; nach Kac ℓ₀/c = 2Λ/c [korr.] |

τ_D ist die Ausbreitungszeit einer Windung, L²/D, nicht die Periode eines Gleichgewichts. Aus der Wirkung folgt τ_D = ℓ₀/c; nach Kac ist τ_D = 2Λ/c, die mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden, mit Λ, der ganzen inneren Faltung der Zelle; ≈ 2π/H₀ gilt nur heute (6–13 % unter dem Messwert). [korr.] Prüfbar ist die Zeitabhängigkeit: Folgte τ_D dem Hubble-Radius, sänke der Tully-Fisher-Nullpunkt bei z = 1 um 0,25 dex in M_b; im Modell ist Λ heute fest und a₀ seit z ≈ 10 fest, der Nullpunkt bleibt. [korr.]
