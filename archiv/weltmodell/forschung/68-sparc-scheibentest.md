# SPARC-Scheibentest: algebraische Relation und Feldform an gemessenen Rotationskurven

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: a₀ = 1,34 und die starke Trennung der Abstandsmethoden (1,10 gegen 1,52) als Ergebnis nur mit errV gewichtet (X²/N ≈ 50) gekennzeichnet; mit den Fehlern von Abstand, Neigung und Masse-zu-Licht a₀ = 1,04–1,19·10⁻¹⁰, Trennung nur 3–13 %; Vorhersage B (Steigung +0,117 ± 0,042) unverändert, mit Hinweis: nur Geschwindigkeitsfehler, mit vollen Fehlern nicht wiederholt, Σ₀ und a₀ hängen beide an Υ_d. Umbenennungen: χ² → X²; Q (Qualitätsstufe) → „Qualitätsstufe“.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/sparc_test.py` (Ausgabe `sparc_test.txt`)
> Daten: SPARC (Lelli, McGaugh, Schombert 2016; astroweb.cwru.edu/SPARC), `daten/sparc/Rotmod_LTG/*.dat` und `SPARC_Lelli2016c.mrt`. Auswahl Qualitätsstufe ≤ 2 [korr.], Neigung ≥ 30°: 153 Galaxien, 3166 Punkte. [A] Annahme.

## Rechnung

```
g_obs = V_obs²/R,   g_bar = (V_gas|V_gas| + Υ_d·V_disk|V_disk| + Υ_b·V_bul²)/R,   Υ_d = 0,5, Υ_b = 0,7 [A], fest
algebraisch:  g = ν(g_bar/a₀)·g_bar,   ν = 1/(1 − e^{−s}),  s = √(g/a₀)        a_N = 0,873 a₀ → ν = 1,6470, a = 1,438 a₀
Feldform:     ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N],   g_F = ν(g_bar/a₀)·g_bar · (a_Feld/a_alg)_Löser
x = 2πGΣ/g_bar,   Σ = Υ_d·SBdisk + Σ_gas,   Fehler σ(ln g_obs) = 2·errV/V_obs
```

- Feldform je Galaxie aus ihrer Flächendichte gelöst (Finite Volumen in (R, z), 440 × 440 Zellen bis 3000 kpc): Sterne aus SBdisk an den Messradien (innen und außen exponentiell mit R_d der Tabelle), Gas als dünne Scheibe aus V_gas zurückgerechnet (Σ ≥ 0, glatt, Gesamtmasse an 1,33·M_HI [A] gebunden), Bulge als Kugel aus V_bul. Sternscheibe dünn oder sech² mit z₀ = 0,15 R_d.
- Löser: Plummer-Kugel Feldform = algebraisch auf 8·10⁻⁵ (bei a_N = 0,878 a₀: 1,4432 a₀, von Hand 1,4432); Kuzmin-Scheibe exakt auf 0,5 %; Exponentialscheibe R_d = r*: 0,606 / 0,717 / 0,829 bei 0,25 / 0,5 / 1 R_d.
- Massenmodell gegen die Daten: Newton-Feld/g_bar = 1,033 bei dünner Sternscheibe (16–84 %: 1,00–1,12), 1,002 bei sech² mit z₀ = 0,15 R_d (0,98–1,03); V_disk der Daten gehört zu einer Scheibe dieser Dicke. V_gas wird auf 0,07 km/s wiedergegeben (Median, 90 % unter 0,5). Gas trägt bei x ≥ 1 im Median 27 % von Σ.
- Näherung in der Ebene: a_Feld/a_alg = ν(√(g_R² + (2πGΣ)²))/ν(g_R) trifft den Löser im Mittel auf −0,001 dex, 90 % der Punkte innerhalb 0,022 dex.
- Einzelwert NGC 2403, R = 0,16 kpc (nur Sterne, a₀ = 1,2·10⁻¹⁰): g_bar = 0,455 a₀, x = 8,9; algebraisch 0,927 a₀, Feldform (Näherung) 0,52 a₀, gemessen 1,013 a₀.

## Vorhersage A: Residuen gegen x

Δ = log g_obs − log[ν(g_bar/a₀)·g_bar]; δ = Δ bei gleichem g_bar (minus einem in log(g_bar/a₀) quadratischen Verlauf, angepasst an die Punkte mit x < 0,3). Fehler durch Ziehen ganzer Galaxien. a₀ = 1,2·10⁻¹⁰; mit cH₀/2π ändern sich alle Werte um höchstens 0,005 dex.

| x | Punkte | δ Daten | δ Feldform dünn / dick |
|---|---|---|---|
| < 0,3 | 898 | 0,000 ± 0,013 | 0,000 / 0,000 |
| 0,3–1 | 1365 | +0,007 ± 0,013 | −0,019 / −0,019 |
| 1–3 | 729 | −0,025 ± 0,018 | −0,095 / −0,086 |
| ≥ 3 | 174 | −0,143 ± 0,039 | −0,240 / −0,198 |

Nach der Vorhersage selbst geordnet (Δ_F dünn): bei −0,02…−0,05 Daten −0,002 gegen −0,034, bei −0,05…−0,1 −0,010 gegen −0,075, bei −0,1…−0,2 −0,040 gegen −0,140, unter −0,2 (119 Punkte, 53 Galaxien, Median x = 4,8) −0,170 ± 0,050 gegen −0,287.

Regression Δ = c₀ + c₁l + c₂l² + β·Δ_F, l = log(g_bar/a₀): Feldform β = 1, algebraisch β = 0 (Probe an Feldform-Werten ohne Rauschen: 1,00).

| Auswahl | β (dünn) | Punkte mit Δ_F < −0,05 |
|---|---|---|
| alle, a₀ = 1,2 / cH₀/2π | 0,43 ± 0,15 / 0,44 ± 0,15 (dick 0,38 ± 0,17) | 871 |
| Υ_d = 0,4 / 0,6 | 0,26 ± 0,14 / 0,55 ± 0,16 | 933 / 821 |
| nur Qualitätsstufe 1 [korr.] | 0,16 ± 0,15 | 584 |
| ohne die drei innersten Punkte je Galaxie | −0,03 ± 0,18 | 548 |
| nur 2·errV/V_obs < 0,1 | −0,29 ± 0,32 | 241 |

Bei großem x fallen die Daten in die Richtung der Feldform, aber nur mit etwa 40 % der vorhergesagten Tiefe. β = 1 liegt 3,8 σ daneben (Qualitätsstufe 1 [korr.]: 5,6 σ; ohne die innersten Punkte: 5,7 σ). β = 0 liegt 2,9 σ daneben, getragen allein von den innersten drei Punkten je Galaxie.

## Vorhersage B: a₀ je Galaxie gegen Σ₀

[korr.] Diese Werte sind nur mit den Geschwindigkeitsfehlern (errV) gewichtet; mit den vollen Fehlern (Abstand, Neigung, Masse-zu-Licht) ist die Steigung nicht wiederholt. Σ₀ = Υ_d·SBdisk und das angepasste a₀ hängen beide an Υ_d.

149 Galaxien mit ≥ 5 Punkten, a₀ je Galaxie bei festem Υ angepasst, Σ₀ = Υ_d·SBdisk(zentral) der Tabelle. Median a₀ in 10⁻¹⁰ m/s²; Vorhersage als a₀,angepasst/a₀,wahr, aus der Feldform an denselben Radien mit denselben Gewichten.

| Σ₀ [M☉/pc²] | < 30 | 30–100 | 100–300 | 300–1000 | > 1000 | d log a₀/d log Σ₀ |
|---|---|---|---|---|---|---|
| Galaxien | 10 | 49 | 26 | 36 | 28 | |
| Daten, algebraisch angepasst | 0,58 | 0,87 | 1,01 | 1,22 | 1,26 | +0,117 ± 0,042 |
| Daten, Feldform angepasst | 0,68 | 1,12 | 1,13 | 1,38 | 1,28 | +0,069 ± 0,044 |
| Vorhersage Feldform (÷ a₀) | 0,77 | 0,74 | 0,83 | 0,89 | 0,98 | +0,051 ± 0,007 |

- Die algebraische Regel sagt überall das wahre a₀ und Steigung 0. Das angepasste a₀ fällt mit Σ₀ in der Richtung der Feldform, doppelt so steil (Daten minus Feldform +0,066 ± 0,043; gegen die algebraische Regel 2,8 σ). Mit der Feldform angepasst bleibt eine Steigung von +0,069 ± 0,044. Streuung je Galaxie 0,35 dex; Abstände aus dem Hubble-Fluss und aus TRGB, Cepheiden, UMa, SN geben dieselbe Steigung (+0,13 / +0,14). Für die Vorhersage geben a₀ = 1,2·10⁻¹⁰ und cH₀/2π sowie dünn und dick dieselben Werte auf 0,01.

## Gesamt: ein a₀ für alle Punkte

| Form | bestes a₀ [10⁻¹⁰] | X²/N beim besten [korr.] | bei 1,2 | bei cH₀/2π (H₀ = 67,4 / 73) |
|---|---|---|---|---|
| algebraisch | 1,339 | 49,7 | 51,3 | 57,8 / 53,5 |
| Feldform dünn | 1,446 | 52,6 | 57,1 | 66,1 / 60,4 |
| Feldform dick | 1,448 | 52,3 | 56,8 | 66,0 / 60,2 |

- Die Feldform verlangt ein um 8 % größeres a₀ und passt schlechter (ΔX² = 8960 dünn, 8010 dick) [korr.]. Am stärksten trennen sich beide bei x = 1–3 (ΔX² = 6610 für die algebraische Relation) [korr.]; bei x ≥ 3 liegt die Feldform knapp vorn (100). Je Galaxie bevorzugen 61 die Feldform, 92 die algebraische Relation; größte Beiträge IC 2574 (+3470) und UGC 2953 (+3140) für die Feldform, UGC 128 (−2790) und NGC 2403 (−1930) für die algebraische. X²/N ≈ 50 [korr.]: errV allein trägt die Unsicherheit in Υ, Abstand und Neigung nicht; die Differenzen hängen an wenigen Galaxien mit kleinem errV.

## Grenzen, gerechnet

| Variante | a₀ algebraisch | a₀ Feldform dünn |
|---|---|---|
| Υ_d = 0,4 / 0,5 / 0,6 (Υ_b fest) | 1,69 / 1,34 / 1,08 | 1,84 / 1,45 / 1,16 |
| Neigungsfehler 2·δi/tan i im Fehler | 1,20 | 1,33 |
| nur Qualitätsstufe 1 (93 Galaxien) [korr.] | 1,20 | 1,29 |
| Abstände aus TRGB, Cepheiden, UMa, SN (70) | 1,10 | 1,21 |
| Abstände aus dem Hubble-Fluss (83; H₀ = 73 in den Daten) | 1,52 | 1,62 |
| alle, Hubble-Fluss-Abstände auf H₀ = 67,4 umgerechnet (g_obs ∝ 1/D) | 1,17 | 1,26 |

- **Υ:** fest angenommen; 0,1 in Υ_d ändert a₀ um 25 % und β um 0,15. Ein kleineres Υ hebt die Residuen der sternbeherrschten inneren Punkte, also genau dort, wo die Feldform senkt.
- **Abstand, Neigung:** Nur mit errV gewichtet (X²/N ≈ 50) trennen die beiden Abstandsmethoden a₀ stärker (1,10 gegen 1,52) als die beiden Formen; mit den Fehlern von Abstand, Neigung und Masse-zu-Licht trennen sie sich nur um 3–13 %, und a₀ liegt bei 1,04–1,19·10⁻¹⁰ (Υ_d = 0,5). [korr.] g_bar und x hängen nicht vom Abstand ab. Die Neigung geht nur als Fehler ein; in den Gewichten senkt sie a₀ um 10 %.
- **Gas:** Aus V_gas allein ist die Gasmasse schlecht bestimmt (ohne Bindung 0,82–2,6 × 1,33·M_HI; einzelne V_gas-Kurven sind unstetig oder null). Die Bindung an M_HI setzt die Masse, die Verteilung bleibt eine Rückrechnung. Ohne Gas in Σ hätten 562 statt 903 Punkte x ≥ 1.
- **Dicke:** Die Dicke von V_disk senkt die Feldform-Vorhersage bei x ≥ 3 von −0,24 auf −0,20 dex, sonst kaum. Gas dünn, Bulge kugelförmig angenommen.
- **Innere Punkte, Umgebung:** Der Rest in β sitzt in den innersten drei Punkten je Galaxie, wo Auflösung und nicht kreisförmige Bewegung die Messung am stärksten betreffen (nicht gerechnet). Kein äußeres Feld.

## Ergebnis

- **A:** Bei gleichem g_bar liegen Punkte mit großem x tiefer, mit etwa 40 % der Tiefe der Feldform (β = 0,43 ± 0,15); die volle Tiefe liegt 3,8–5,7 σ neben den Daten. Ohne die innersten drei Punkte je Galaxie verschwindet die Abhängigkeit (β = −0,03 ± 0,18).
- **B:** Das je Galaxie angepasste a₀ fällt mit Σ₀ (Steigung +0,12 ± 0,04, nur mit errV gewichtet [korr.]), steiler als die Feldform voraussagt (+0,05); die algebraische Regel sagt 0.
- **Gesamt:** algebraisch a₀ = 1,34, Feldform 1,45 (Υ fest, Gewichte nur aus errV, X²/N ≈ 50); mit den Fehlern von Abstand, Neigung und Masse-zu-Licht liegt a₀ bei 1,04–1,19·10⁻¹⁰ (Υ_d = 0,5). [korr.] Die algebraische Relation passt besser, am deutlichsten bei x = 1–3.
- **Schwelle:** Das angepasste a₀ hängt stärker an Υ, Gewichten und Abstandsmethode (algebraisch 1,08–1,69) als an der Form; mit Abständen aus TRGB, Cepheiden, UMa, SN ist es 1,10 (algebraisch) bzw. 1,21 (Feldform) gegen cH₀/2π = 1,04–1,13.
- **Näherung:** In der Ebene realer Scheiben gibt ν(√(g_R² + (2πGΣ)²))·g_R die Feldform auf 0,02 dex (90 % der Punkte).
