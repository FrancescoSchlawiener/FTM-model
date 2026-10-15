# Radiale Beschränkung der Rückkehr

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: [G] an R₀ = 8,2 kpc und am Literaturwert 71 ± 6 M☉/pc² (Kuijken & Gilmore 1991); „fehlende Scheinmasse“ → gelesene Masse. Umbenennungen: w (Glättungstiefe) → ϖ; χ² → X².

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/radiale_beschraenkung.py` (Ausgabe `radiale_beschraenkung.txt`)
> Daten: SPARC (Lelli, McGaugh, Schombert 2016), 153 Galaxien, 3166 Punkte, Υ_d = 0,5, Υ_b = 0,7 [A]. Sonnenort: Σ_bar(|z| < 1,1 kpc) = 47,1 M☉/pc², g_N,R = 1,374·10⁻¹⁰ m/s² (Punktmasse 1,244), R₀ = 8,2 kpc [G]; Literaturwert 71 ± 6 M☉/pc² (Kuijken & Gilmore 1991) [G]. [korr.]

## Zwei Wahlen in jeder Lesart

Jede Lesart legt fest, **welches Gefälle** den Verlust Γ setzt (Argument von ν), und **auf welche Richtung** der Bewegung die Verstärkung wirkt (Anwendung). Rotationskurven prüfen das Argument; die Dynamikmasse am Sonnenort prüft die Anwendung, denn dort ist x = 2πGΣ/g_R = 0,30 und das Argument ändert ν nur um 1 %.

```
F  Feldform:     ∇²Φ = ∇·[ν(|∇Φ_N|)∇Φ_N]                         aus der Wirkung, wirbelfrei
a  Kernform:     a = g_N + (ν(g_R) − 1)(e_r⊗e_r)g_N,   e_r vom Zentrum
b  Tensor:       a_i = ν(|g_i|)·g_i                                je Torusrichtung eine Rückkehr
c  Torusmittel:  a = ν(⟨|g_N|⟩_Bahn)·g_N                           Rückkehr als Vorgang an der Bahn
d  Bahnregel:    a = g_N + (ν − 1)(n⊗n)g_N = ν(|g_N(x)|)·g_N(x),   n = Dσ/|Dσ| am Ort des Sterns
```

- **Die eine Rückkehr gehört zu keiner Torusrichtung.** Sie läuft in der Glättungstiefe ϖ [korr.] und wird an dem Ort gezählt, an dem die Bahn gerade ist; Γ ist der Betrag des Gefälles dort. Rückkehren je Torusrichtung, gemeinsam gezählt, gäben d_s = 2 oder 3 und keine Wurzel; getrennt gezählt geben sie b.
- **In der Mittelebene ist das Gefälle radial.** Aus der Spiegelung folgt g_z(0) = 0, also n = e_R und n⊗n = e_R⊗e_R. Dort braucht die radiale Beschränkung keine Regel; sie ist die Symmetrie des Orts.
- **Die Tiefe der Feldform kommt aus der Feldgleichung.** Die Wirkung gibt nicht ν∇Φ_N am Ort, sondern dessen wirbelfreien Teil; der Unterschied ist der Wirbelanteil von ν∇Φ_N (∇ × (ν∇Φ_N) = ∇ν × ∇Φ_N). Er trägt das über und in der Schicht abgesenkte ν (|g_N|² = g_R² + g_z²) in die Ebene: Die Feldform liest die Rückkehr über die Umgebung des Orts, die Bahngleichung am Ort des Sterns.
- **Schichtmittel:** In einer selbstgravitierenden Schicht ist |g_z|/2πGΣ unter der Dichte gleichverteilt auf [0, 1], für jedes Profil: ⟨|g_z|⟩ = πGΣ, ⟨g_z²⟩ = (2πGΣ)²/3. Für einen Tracer mit der Verteilung der Massenschicht gilt in d
  a_R/[ν(g_R)g_R] = ∫₀¹ ν(g_R√(1 + x²u²)) du / ν(g_R); ein Tracer in der Mittelebene gibt genau die radiale Form.

## Dicke Scheibe

Exponentialscheibe, sech² mit z₀ = 0,15 R_d, R_d = r*, Löser in (R, z); Verhältnis zur radialen Form ν(g_R)·g_R:

| R/R_d | 0,25 | 0,5 | 1 | 2 |
|---|---|---|---|---|
| Feldform (Löser) | 0,726 | 0,787 | 0,859 | 0,940 |
| örtliche Regel in der Mittelebene | 1,000 | 1,000 | 1,000 | 1,000 |
| d, Tracer = Massenschicht | 0,757 | 0,850 | 0,919 | 0,967 |
| Feldform dünn, ν(√(g_R² + (2πGΣ)²)) | 0,557 | 0,681 | 0,803 | 0,915 |

Die Feldform senkt, obwohl das Gefälle in der Ebene radial ist. Innen liegt sie nahe beim Schichtmittel d, außen zwischen d und der dünnen Scheibe (R_d = 0,3 und 3 r* ebenso; bei R_d/4 Feldform 0,854 und 0,683, d 0,862 und 0,721).

## SPARC

Lesarten in der Schichtnäherung mit g_R = g_bar; c und d mit Tracer = Massenschicht (obere Grenze, ein Tracer in der Mittelebene gibt 0). β = Anteil der Tiefe der dünnen Feldform (Regression mit quadratischem Verlauf in log g_bar), δ = Abweichung bei gleichem g_bar.

| | Daten | F | a, b | c | d |
|---|---|---|---|---|---|
| β, alle Punkte | 0,39 ± 0,14 | 1 | 0 | 0,66 | 0,55 |
| β, ohne die drei innersten je Galaxie | 0,08 ± 0,18 | 1 | 0 | 0,55 | 0,47 |
| δ bei x = 1–3 (729 Punkte) | −0,025 ± 0,018 | −0,102 | 0 | −0,050 | −0,044 |
| δ bei x ≥ 3 (174 Punkte) | −0,143 ± 0,039 | −0,282 | 0 | −0,189 | −0,157 |
| bestes a₀ [10⁻¹⁰ m/s²] | | 1,47 | 1,34 | 1,39 | 1,39 |
| X² − X²(radial), X²/N ≈ 50 [korr.] | | +6669 | 0 | +865 | +435 |

## Sonnenort, Kugeln, äußeres Feld

| Lesart | Ebene (β) | Kugeln | ν_z² in der Ebene | Σ_dyn [M☉/pc²] | ν bei 50 kpc: senkrecht / in der Ebene |
|---|---|---|---|---|---|
| F | 1 | gleich | ν·ν_z,N² | 70,8 (72,6) | 5,46 / 4,68 |
| a, e_r vom Zentrum | 0 | gleich | ν_z,N² + (ν−1)ω_N² | 58,1 (57,9) | 6,50 / 4,68 |
| b | 0 | nicht kugelförmig | a_z ∝ √z, ν_z divergiert | 106 | 6,50 / 4,68 |
| c | 0 … 0,66 | Kreisbahnen gleich, druckgestützte anders | ν·ν_z,N² | 70,8–71,7 | 5,46 / 5,72 |
| d | 0 … 0,55 | gleich | ν·ν_z,N² | 70,8 (72,6) | 5,46 / 4,68 |

Klammern: Punktmasse. Milchstraße 6·10¹⁰ M☉, g_e = 0,03 a₀, frei ν = 6,50; in der Ebene bei F, a, b, d am gleichgerichteten Punkt, bei c gemittelt über die Kreisbahn. b schräg zur Achse: |a| 10–18 % größer, Richtung bis 7° neben radial. Σ_dyn der Feldform ist hier örtlich genähert (a_z = ν(|g_N|)·g_z) und fällt so mit d zusammen; bei Anwendung längs g_N ist Σ_dyn ≈ ν·Σ_N,dyn, die Spanne 71–95 der Feldform kommt aus dem Baryonenmodell.

## Ergebnis

1. **Nichts beschränkt die Rückkehr auf eine Richtung.** Sie ist eine, am Ort der Bahn gezählt, mit dem Betrag des Gefälles dort; in der Mittelebene ist dieser Betrag g_R. Die Tiefe der Feldform entsteht erst durch die wirbelfreie Ergänzung der Feldgleichung.
2. **Lesart d trägt beides:** radiale Rotationskurve für Tracer in der Ebene, mit Tracer = Massenschicht 0,55 der Tiefe (Daten 0,39 ± 0,14; bei x ≥ 3 −0,157 gegen −0,143 ± 0,039), Dynamikmasse 70,8 (Kernform mit e_r 58,1), Kugeln und äußeres Feld wie im Modell.
3. **Herkunft:** d ist die Bahngleichung mit n⊗n, stationär am Ort des Sterns. e_r⊗e_r und n⊗n stimmen in der Ebene und bei Kugeln überein; über der Ebene braucht die Dynamikmasse n⊗n. Aus der Wirkung folgt d nicht: Deren Feldgleichung gibt den wirbelfreien Teil, also F. Für Scheiben ist die Bahngleichung damit keine Näherung der Wirkung; d als Grundgleichung ist eine eigene Setzung (die Regel wirkt auf Bahnen, nicht als Potential).
4. **Wo es nicht trägt:** Ohne die drei innersten Punkte zeigen die Daten keine Tiefe (0,08 ± 0,18), d mit Tracer = Schicht 0,47; das passt nur, wenn die Tracer dort dünner liegen als die Massenschicht (Gas in der Ebene). Die Dicke der Tracer steht nicht in den Daten; ohne freie Zahl gibt d nur das Band 0 … 0,55. Global bevorzugt X² die radiale Form (+435) [korr.]. d hat keinen Wirkungsterm: Die Kraft enthält den Wirbelanteil ∇ν × ∇Φ_N, Energie- und Impulsbilanz sind nicht durch eine Wirkung gesichert; prüfbar ist ∂_z a_R ≠ ∂_R a_z, was keine Massenverteilung erzeugt.
5. **Die anderen Lesarten:** b hängt an Achsen, lässt ν_z divergieren und gibt 106 M☉/pc²; a mit e_r vom Zentrum gibt 58 M☉/pc²; c ist tiefer als d (0,66), ändert druckgestützte Kugeln und das äußere Feld in der Ebene (5,72 statt 4,68) und verlangt die Rückkehr als Vorgang längs der Bahn.
6. **Beobachterfehler:** Die Beschränkung ist der Ort, an dem die Rückkehr gezählt wird, eine Regel. Die Dicke des Tracers ist die Vertikalwirkung J_z seiner Bahnen, ein Zustand des Torus. Wer die Tiefe der Feldform als fehlende Masse in der Scheibe oder den Überschuss in Σ_dyn als dunkle Scheibe liest (gelesene Masse), macht aus der Regel einen Stoff [korr.]; den Wirbelanteil von d kann keine Dichte darstellen.
