# Wahl des Imaginärteils von τ_M: wo (A) und (B) sich trennen

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. (A) entspricht Punkt 3 des Autors (|q_Δs| = Δs/r, Schalter χ), (B) Punkt 1 (|q| = e^{−r*/r}). Geändert: „Trennende Messungen“ über das Gedächtnisintegral (Abschnitt 6 und Ergebnis: Sprung von a_N, (A) sofort, (B) mit 1/λ_m; Instabilitätsband 0,07–0,55 r* in (B); Radialstreuung alter Scheiben) → das Gedächtnis wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; stationär sind beide Wahlen auf allen Bahnen gleich, nach einer Massenänderung gibt es keine messbare Verzögerung; die Wahl betrifft die Schreibweise. Merkur „10⁻⁷⁸⁹¹ m/s²“ → s = 18160, Zusatz 10⁻⁷⁸⁸⁸ m/s². Umbenannt: Nome von (A) q_A → q_Δs (Kern-Nome); Nome von (B) q_B → q (für die Punktmasse der Rückkehranteil).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/wahl_im_tau.py` (Ausgabe `wahl_im_tau.txt`)
> Schwelle a₀ = 1,2·10⁻¹⁰ m/s². Beispiel: baryonische Milchstraße 6·10¹⁰ M☉, r* = 8,35 kpc, v_f = (GMa₀)^{1/4} = 175,8 km/s.

## Die beiden Wahlen

```
(A)  |q_Δs| = Δs/r,       Im τ_M = ln(r/Δs)/2π,   Zusatz = χ(s)·Δ(τ_M)·3c²/r*,      χ = s/(e^s − 1)
(B)  |q| = e^{−s},        Im τ_M = s/2π,          Zusatz = a_N·Σ|q|ⁿ = a_N·|q|/(1 − |q|)
s = √(a_N/a₀)  (Punktmasse: s = r*/r),   Re τ_M = ω_rad/ω_azim in beiden gleich
```

[korr.] q_A heißt jetzt q_Δs (Kern-Nome), q_B heißt q; für die Punktmasse ist |q| = e^{−r*/r} der Rückkehranteil.

## 1. Auf Kreisbahnen sind beide dieselbe Größe

χ(s)·(Δs/r)·(3c²/r*) = a_N·|q|/(1 − |q|), also χ(s) = s·|q|/(1 − |q|). [korr.] Die Gleichung gilt exakt (numerisch auf 6·10⁻¹⁵ für 10⁸–10¹⁴ M☉ und r/r* = 0,01–100). Der Schalter von (A) ist der Nome von (B), mal s. Unterscheiden können sich beide nur über (a) die höheren Glieder der q-Reihe, (b) die Phase des Nomes, (c) die Stelle, an der der Schalter im Gedächtnisintegral steht.

## 2. Höhere Fourier-Glieder

**(A):** Δ/q_Δs = 1 − 24q_Δs + 252q_Δs² − …, mit q_Δs = (v_f²/3c²)·s [korr.]. Die Beschleunigung ändert sich um a_A/a_B − 1 = (Δ/q_Δs − 1)(1 − 1/ν), am stärksten bei r = r*: max |a_A/a_B − 1| = 24·(v_f²/3c²)/e = 2,94·(v_f/c)².

| M [M☉] | v_f [km/s] | Im τ_A bei r* | max \|a_A/a_B − 1\| |
|---|---|---|---|
| 10⁸ | 35,5 | 3,05 | 4,1·10⁻⁸ |
| 6·10¹⁰ | 175,8 | 2,54 | 1,0·10⁻⁶ |
| 10¹² | 355 | 2,32 | 4,1·10⁻⁶ |
| 10¹⁴ | 1123 | 1,95 | 4,1·10⁻⁵ |
| 10¹⁵ | 1998 | 1,77 | 1,3·10⁻⁴ |

In (A) hängt Im τ_M bei r* von der Masse ab, die Beschleunigungsrelation bleibt universell bis auf 2,94(v_f/c)². Die Koeffizienten τ_R(n) = −24, 252, −1472, … sind in keinem System messbar.

**(B):** Im τ_B = s/2π hängt nicht von der Masse ab; die Relation ist exakt universell. Die Bindung ist die Reihe mit allen Koeffizienten 1. Setzt man Δ(τ_B) selbst als Kern ein, bleibt nur das erste Glied gleich:

| r/r* | 0,1 | 0,2 | 0,304 | 1 | 3 |
|---|---|---|---|---|---|
| Δ(τ_B) / Σ\|q\|ⁿ [korr.] | 0,999 | 0,84 | 0,37 | 4,7·10⁻⁸ | 3·10⁻³⁷ |

Δ(τ_B) hat sein Maximum bei 0,304 r* und fällt danach mit Δ(i/t) = t¹²Δ(it) ab; ab r = r*/2π liegt τ_B jenseits des selbstdualen Punkts τ = i. In (B) ist Δ deshalb nur die Nahfeldform (Unterschied 25|q|² für r ≪ r*) [korr.], die Bindung selbst trägt die geometrische Reihe.

## 3. Phase des Nomes

Mit Re τ = √(1 + χ) (Punktmasse) ändert die Phase in (A) nur die höheren Glieder: |Δ(τ)|/Δ(i Im τ) − 1 = 2,9·10⁻⁶ bei r*, 4,9·10⁻⁷ bei 10 r*. In (B) wirkt sie auf die ganze Bindung: statt ν = 1/(1 − |q|) = 1,582 (r*) bzw. 10,51 (10 r*) gäbe Re 1/(1 − q) 0,869 bzw. 0,528 und |1/(1 − q)| 0,924 bzw. 0,554. (B) verlangt, dass die Phase in keiner Ordnung eingeht; (A) verlangt es nur für das erste Glied.

## 4. Resonanzradien

Re τ_M ist in beiden gleich, also auch die Radien mit Re τ = m/n. Die Δ-Glieder von (A) verschieben sie um höchstens 5·10⁻⁶ relativ. Verschieden ist die Größe des n-ten Glieds der Reihe, das dort phasengleich wird:

| m/n | r [kpc] | \|τ_R(n)\|·\|q_Δs\|ⁿ [korr.] | \|q\|ⁿ [korr.] |
|---|---|---|---|
| 6/5 | 5,70 | 7·10⁻³¹ | 6,6·10⁻⁴ |
| 5/4 | 7,89 | 3·10⁻²⁵ | 1,5·10⁻² |
| 4/3 | 17,28 | 4·10⁻²⁰ | 0,24 |
| 7/5 | 102,95 | 3·10⁻³⁷ | 0,67 |

Eine Feinstruktur der Rotationskurven an diesen Radien in dieser Größe gäbe es nur in (B), und nur wenn die Phase eingeht.

## 5. Sonnensystem

Beide geben denselben Zusatz: 10⁻⁷⁸⁸⁸ m/s² an Merkur (s = 18160) [korr.], 10⁻³⁰⁵⁵ an der Erde, 10⁻¹⁰⁷ an Neptun. Dort ist Im τ_A = 2,96 / 3,11 / 3,65, Im τ_B = 2890 [korr.] / 1119 / 37. In (A) ist der Nome q_Δs = GM/(3c²r) [korr.], genau ein Drittel des post-Newtonschen Parameters GM/(c²r), und τ_M liegt bei Im τ ≈ 3; der Kern verschwindet allein über χ. Ginge τ_M in (A) an einer weiteren Stelle ohne χ ein, entstünden Glieder der Größe q_Δs ≈ 10⁻⁹–10⁻⁸ [korr.], also von der Ordnung der post-Newtonschen Effekte. In (B) steht τ_M an der Spitze; jede Funktion von τ_M ist dort e^{−s}-klein.

## 6. Gedächtnisintegral: wo der Schalter steht

[korr.] Diese Rechnung setzt ein Gedächtnis der einzelnen Bahn voraus (Geschwindigkeitskern mit Nachlauf 1/λ_m auf Galaxienbahnen); das trägt nicht: Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m = ω_orb/N*, N* ≈ 0,3–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Geschwindigkeitskern würde exzentrische Bahnen runden, alte Scheiben verlangen dafür λ_m ≪ 0,01·Ω. In Galaxien gibt es keinen messbaren Nachlauf; nach einer Massenänderung folgt die Regel ohne messbare Verzögerung. Stationär sind beide Wahlen auf allen Bahnen gleich (χ(s) = s·q/(1 − q)); die Wahl betrifft die Schreibweise. Die Trennungen dieses Abschnitts gelten nur unter der Annahme.

In der Feldgleichung steht χ vor dem Integral, Δ(τ_M) darin. In (A) wirkt der Schalter also zur Zeit t, der Nome längs der Bahn; in (B) ist der Schalter Teil des Nomes und läuft mit der Bahngeschichte. In Einheiten GM = r* = a₀ = 1:

```
r̈ = −e_r/r² − C·κ(r)·(e_r·m)·e_r,     ṁ = K·v − λ_m·m,     κ(r) = (λ_m² + ω_c²)/(r·ω_c²)
(A): C = χ(s(t)),  K = s             (B): C = 1,  K = s·χ(s) = s²|q|/(1 − |q|)      [korr.]
```
κ(r) ist an der Kreisbahn am Ort geeicht; auf jeder Kreisbahn geben beide exakt a = a_N/(1 − e^{−s}).

**Sprung des Baryonenfelds.** Ändert sich a_N an einem Ort plötzlich (a_N → f·a_N), springt in (A) der Schalter sofort, in (B) nichts. Unmittelbar danach gilt Zusatz_A/Zusatz_B = χ(s₂)/χ(s₁):

| s₁ | f = 0,5 | f = 2 | Endwert (beide), f = 0,5 / 2 |
|---|---|---|---|
| 0,3 | 1,05 | 0,94 | 0,74 / 1,32 |
| 1 | 1,18 | 0,78 | 0,84 / 1,10 |
| 3 | 1,84 | 0,39 | 1,30 / 0,56 |

Im Übergang läuft die Sofortantwort von (A) dem Endwert entgegen (s₁ ≤ 1) oder schießt über ihn hinaus (s₁ = 3); (B) läuft nur mit 1/λ_m nach. Das gilt nur mit einem Gedächtnis der einzelnen Bahn; im Modell folgt die Regel einer Massenänderung in beiden Wahlen ohne messbare Verzögerung. [korr.]

**Kreisbahnen.** Mit Gedächtnis zerfällt die radiale Schwingung in beiden Wahlen in zwei Moden (bei r*, λ_m/Ω = 0,01: (A) 1,50 und 0,84 Ω, (B) 1,38 und 0,91 Ω; ohne Gedächtnis 1,26 Ω). In (A) sind alle Moden gedämpft (λ_m/Ω = 10⁻³–10², r/r* = 0,01–100). In (B) wachsen sie in einem Band:

| λ_m/Ω | 0,001 | 0,01 | 0,1 | 0,3 | 1 | 2 |
|---|---|---|---|---|---|---|
| größte Rate [Ω] | 0,114 | 0,109 | 0,073 | 0,031 | 0,0027 | 0,0002 |
| Band (Rate > 10⁻³Ω) [r*] | 0,07–0,55 | 0,08–0,56 | 0,10–0,55 | 0,12–0,46 | 0,14–0,24 | – |

Mit κ fest an der Ausgangsbahn statt am Ort bleibt das Band (0,08–0,50 r*, 0,104 Ω). Nichtlinear wächst in (B) bei 0,35 r* eine Auslenkung von 10⁻⁴ in 15 Umläufen auf 1,2; in (A) klingt sie auf 7·10⁻⁵ ab. In der Milchstraße ist das Band 0,6–4,6 kpc, 1/Ω(0,35 r*) = 9,3 Myr; jede Nachlaufzeit 1/λ_m ≥ 0,1 Gyr gibt λ_m/Ω ≤ 0,09 und e-Faltung ≤ 0,13 Gyr.

**Exzentrische Bahnen** (Gedächtnis auf der Bahn im lokalen Feld eingeschwungen, dann volle Dynamik). Hub r_a − r_p nach 20 Radialperioden relativ zum Anfang:

| r_p–r_a [r*] | λ_m/Ω = 0,01: (A) | (B) | λ_m/Ω = 0,1: (A) | (B) |
|---|---|---|---|---|
| 0,2–1 | 0,45 | 2,5 (nach 18; r_a × 2,2) | 0,004 | 1,29 |
| 0,3–1,5 | 0,57 | 0,33 | 0,005 | 0,27 |
| 1–5 | 0,52 | 0,36 | 0,007 | 0,007 |
| 3–15 | 0,46 | 0,41 | 0,002 | 0,008 |

Außerhalb des Bands wird der Hub in beiden kleiner; (A) und (B) unterscheiden sich dort nur um Faktoren bis 1,6 in der Abklingzeit. Bahnen, die im Band laufen, bleiben in (B) exzentrisch oder werden weiter, in (A) werden sie rund.

## Ergebnis

| Stelle | (A) | (B) | trennende Messung |
|---|---|---|---|
| höhere Glieder, Masse | ≤ 2,94(v_f/c)²: 10⁻⁶ Galaxien, 10⁻⁴ Haufen | exakt universell | keine erreichbar |
| Phase, Resonanzradien | ≤ 3·10⁻⁶ | Größe 1; \|q\|ⁿ = 7·10⁻⁴–0,7 an m/n [korr.] | Feinstruktur an Re τ = m/n (nur falls Phase eingeht) |
| Sonnensystem | 10⁻³⁰⁵⁵ m/s² | gleich | keine |
| Sprung von a_N | sofort χ(s₂)/χ(s₁), z. B. 1,84 bzw. 0,39 (s₁ = 3, f = ½ bzw. 2) | nur Nachlauf 1/λ_m | nur mit Bahngedächtnis; im Modell keine messbare Verzögerung, keine Trennung [korr.] |
| Kreisbahnen 0,07–0,55 r* | gedämpft | wachsend bis 0,11 Ω, falls λ_m ≲ Ω | nur mit Bahngedächtnis; das Gedächtnis wirkt auf den Zustand der Regel, keine Trennung [korr.] |

Stationär sind beide Wahlen bis auf 2,94(v_f/c)² gleich. Zeitabhängig trennten sie sich nur über die Stelle des Schalters in einem Gedächtnisintegral der einzelnen Bahn. Das Gedächtnis wirkt aber auf den Zustand der Regel, nicht auf die einzelne Bahn: In der stationären Regel sind beide Wahlen auf allen Bahnen gleich, nach einer Massenänderung gibt es keine messbare Verzögerung, und die Wahl betrifft die Schreibweise. [korr.]
