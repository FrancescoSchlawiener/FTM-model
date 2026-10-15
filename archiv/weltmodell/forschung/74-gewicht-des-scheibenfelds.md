# Gewicht des Scheibenfelds in der Regel

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Hinweis bei den Ergebnissen: Gegenprobe (Jackknife über Galaxien) gibt f = +0,03 bis +0,13, ohne IC 2574 0,00 ± 0,03, a₀ bei f = 0 und Υ_d = 0,5 1,08–1,19·10⁻¹⁰ (diese Werte gelten); Gegenprobe N (Υ, D, i je Galaxie frei) liest eine Regeländerung als Masse und entfällt für a₀. Umbenennungen: χ² → X²; Q (Qualitätsstufe) → „Qualitätsstufe“.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/gewicht_scheibenfeld.py` (Ausgabe `gewicht_scheibenfeld.txt`)
> Daten: SPARC (Lelli, McGaugh, Schombert 2016), `daten/sparc/Rotmod_LTG/*.dat` und `SPARC_Lelli2016c.mrt`; Qualitätsstufe ≤ 2 [korr.], Neigung ≥ 30°: 153 Galaxien, 3166 Punkte; Primärabstand (f_D = 2–5: TRGB, Cepheiden, UMa, SN) 70 Galaxien, Hubble-Fluss (f_D = 1) 83. [A] Annahme, [G] aus dem Gedächtnis.

## Familie

```
a = ν(√(g_R² + f·(2πGΣ)²)/a₀)·g_R,   ν = 1/(1 − e^{−s}),  s = √(|g|/a₀)       a_N = 0,873 a₀ → ν = 1,6470, a = 1,4378 a₀
f = 0: radiale Form;  f = 1: Feldform der dünnen Scheibe (in der Ebene auf 0,02 dex)
g_R = (V_gas|V_gas| + Υ_d·V_disk|V_disk| + Υ_b·V_bul²)/R,   Σ = Υ_d·SBdisk + Σ_gas,   Υ_b = 0,7 [A]
Σ_gas aus V_gas als dünne Scheibe zurückgerechnet (V_gas auf 0,07 km/s im Median)
f < 0: Regel in ln ν gespiegelt (die eigene Flächendichte höbe ν), misst nur die Richtung
```

## Fehler

- Abstand: g_obs ∝ 1/D, g_R und Σ hängen nicht von D ab, also σ = e_D/D (Median 0,08 Primär, 0,30 Hubble-Fluss).
- Neigung: V_obs ∝ 1/sin i, also σ = 2·e_i·cot i (e_i Median 3°; die Neigung wirkt nur auf V_obs [A]).
- Υ: 0,1 dex je Galaxie [A], also σ = (∂ln a/∂ln Υ)·0,23; groß nur, wo Sterne tragen, also gerade bei großem x = 2πGΣ/g_R.
- **P (Hauptrechnung):** die drei Anteile quadratisch zu 2·errV/V addiert; X²/FG 0,9–2,8 statt 50 mit errV allein. Weil sie je Galaxie gemeinsam wirken, kommen die Grenzen aus dem Ziehen ganzer Galaxien (1000-mal), nicht aus ΔX². [korr.]
- **N (Gegenprobe):** Abstand, Neigung und Υ je Galaxie mit denselben Vorgaben mitgeschätzt (je Galaxie gegen scipy auf 0,001 in X²). X²/FG bleibt 2,7–6,1 [korr.], und 35–36 Galaxien stoßen an die Grenze |δΥ| = 0,3 dex [A]: Bei kleinem errV setzen wenige Galaxien das Ergebnis.

## 1. Verlangtes f (P; bestes f, 16–84 % / 2,5–97,5 % beim Ziehen)

| Υ_d | Punkte | Primärabstand | Hubble-Fluss | alle |
|---|---|---|---|---|
| 0,4 | alle | +0,12 (0,00…0,40 / −0,04…0,78) | −0,04 (−0,17…0,00 / −0,40…0,05) | +0,03 (−0,01…0,16 / −0,04…0,39) |
| 0,5 | alle | +0,21 (0,01…0,56 / −0,03…1,06) | −0,03 (−0,13…0,01 / −0,36…0,08) | +0,07 (0,00…0,25 / −0,03…0,53) |
| 0,6 | alle | +0,30 (0,03…0,71 / −0,02…1,31) | −0,01 (−0,11…0,03 / −0,34…0,12) | +0,12 (0,01…0,35 / −0,03…0,74) |
| 0,4 / 0,5 / 0,6 | ohne innerste drei | +0,15 / +0,21 / +0,27 | −0,22 / −0,19 / −0,15 | +0,05 / +0,09 / +0,13 |

- Primärabstände verlangen ein kleines positives Gewicht, Hubble-Fluss-Abstände keines oder ein negatives. Ohne die innersten drei Punkte bleibt die Primärgruppe gleich; die Hubble-Gruppe fällt unter 0 (obere 68-%-Grenze −0,07 bis +0,01).
- Das positive f der Primärgruppe trägt eine Galaxie: ohne IC 2574 f = +0,02 (−0,01…0,14 / −0,05…0,44), alle ohne IC 2574 +0,01 (−0,02…0,05 / −0,07…0,12). Je Galaxie streut das beste f über den ganzen Bereich (Υ_d = 0,5: 71 Galaxien unter 0, 79 über 0).
- N: Primär +0,05 / +0,09 / +0,14 (97,5 % bei 0,13 / 0,22 / 0,39), Hubble-Fluss am Rand −1, getragen von UGC 128, NGC 6015, UGC 5764; alle +0,01 bis +0,03, ohne die innersten drei unter 0. errV allein, Υ fest 0,5: Primär +0,25, Hubble-Fluss −0,19, alle +0,06.

## 2. Begründete Werte ohne Anpassung

In einer symmetrischen Schicht ist g_z(z) = 2πGΣ·F(z), F der Anteil der Säule innerhalb |z|. Über die Dichte gemittelt ist F gleichverteilt auf [0, 1], für jedes Profil (gleichförmig, sech², exponentiell nachgerechnet):

```
⟨|g_z|⟩ = πGΣ                 → f = 1/4    (gleichförmige Schicht und sech² gleich)
⟨g_z²⟩ = (2πGΣ)²/3            → f = 1/3
⟨ν(√(g_R² + g_z²))⟩ gemittelt  → wirksam f = 0,26 (Median, Punkte mit x ≥ 1; 5–95 %: 0,15–0,29)
Mittelebene g_z = 0            → f = 0;   Gas in eigener sech²-Schicht, h = z₀/2 bzw. z₀/4 [A]: f = 0,09–0,14 bzw. 0,03–0,05
f = ¼ gibt 33 % (x = 1), 54–57 % (x = 3), 65–70 % (x = 5) der Tiefe von f = 1
```

| Wert | Primärabstand | Hubble-Fluss | alle |
|---|---|---|---|
| f = 0 | innerhalb 95 %, knapp außerhalb 68 % (f ≤ 0 in 7–16 %) | innerhalb 68 % | innerhalb 68 % (Υ_d = 0,6: 95 %) |
| f = ¼, ⅓, ⟨ν⟩ | innerhalb 68 %; ohne IC 2574 innerhalb 95 % | außerhalb 95 % | innerhalb 95 %; ohne IC 2574 außerhalb |
| f = 1 | außerhalb 95 % bei Υ_d = 0,4, innerhalb bei 0,5–0,6 | außerhalb 95 % | außerhalb 95 % |

In der Primärgruppe passt ⟨ν⟩ so gut wie das beste feste f (ΔX²/(X²/FG) von −7 bis +1) [korr.], in der Hubble-Gruppe um 58–113 schlechter.

## 3. f und a₀ (a₀ in 10⁻¹⁰ m/s² bei f = 0 / ¼ / 1; Υ_d = 0,5, alle Punkte)

| Behandlung | Primärabstand | Hubble-Fluss | alle |
|---|---|---|---|
| errV allein, Υ fest | 1,10 / – / 1,24 | 1,52 / – / 1,64 | 1,34 / – / 1,47 |
| P: Υ fest, Anteile als Fehler | 1,04 / 1,15 / 1,32 | 1,18 / 1,29 / 1,45 | 1,08 / 1,18 / 1,35 |
| N: Υ, D, i je Galaxie frei (entfällt für a₀, liest eine Regeländerung als Masse) [korr.] | 1,20 / 1,15 / 1,07 | 1,24 / 1,10 / 0,92 | 1,22 / 1,12 / 0,97 |

- Bei festem Υ steigt a₀ mit f: ln[a₀(1)/a₀(0)] = +0,15 bis +0,24 in allen Varianten, bei f = ¼ um 6–10 %. Ist Υ je Galaxie frei, fällt a₀ mit f (−0,07 bis −0,14 Primär, −0,27 bis −0,31 Hubble-Fluss): Ein größeres f senkt die inneren, sternbeherrschten Punkte, die Anpassung hebt Υ, außen wird weniger Verstärkung verlangt. Die Kopplung läuft über Υ.
- Beim jeweils besten f (P): Primär 1,33 / 1,13 / 0,98, Hubble-Fluss 1,46 / 1,17 / 0,96 für Υ_d = 0,4 / 0,5 / 0,6. Υ verschiebt a₀ stärker als f und als die Abstandsmethode.
- Mit realistischen Fehlern liegen die Methoden nah beieinander: bei f = 0 1,04 gegen 1,18 (P) und 1,20 gegen 1,24 (N), beim jeweils besten f 1,13 gegen 1,17 (P), statt 1,10 gegen 1,52 mit errV allein; dort setzen wenige Galaxien den Wert. cH₀/2π = 1,04 für H₀ = 67,4 [G], 1,13 für H₀ = 73 (Annahme der Hubble-Fluss-Abstände in den Daten).

## 4. Beobachterfehler

- f ist ein Gewicht in der Regel: wie stark das eigene senkrechte Feld der Scheibe in das Argument von ν eingeht. Es ändert keine Masse. f < 1 heißt nicht, die Scheibe habe weniger, dickere oder verborgene Flächendichte (Σ_wirksam = √f·Σ wäre die Lesart als Stoff); f < 0 ist kein negativer Stoff, nur die Richtung der Abweichung.
- Die Verwechslung sitzt in der Anpassung selbst: Ist Υ je Galaxie frei, nimmt die Sternmasse einen Teil des Gewichts auf (N: kleineres f, umgekehrte Richtung von a₀). Eine Regeländerung wird dann als Masse gelesen. Darum trägt die Hauptrechnung Υ als Fehler, nicht als freie Masse.

## Ergebnis

[korr.] Die Gegenprobe (Rechnung `rechnungen/gewicht_scheibenfeld_gegenprobe.py`, Jackknife über Galaxien) gibt zusammen f = +0,03 bis +0,13 (Υ_d = 0,4–0,6), ohne IC 2574 0,00 ± 0,03 (Υ_d = 0,5), und a₀ bei f = 0 und Υ_d = 0,5 1,08–1,19·10⁻¹⁰; diese Werte gelten. Die Gegenprobe N (Υ, D, i je Galaxie frei) liest eine Regeländerung als Masse (Abschnitt 4) und entfällt deshalb für a₀.

- Die Daten verlangen kein einheitliches Gewicht: mit Primärabständen f = +0,1 bis +0,3, mit Hubble-Fluss-Abständen f ≈ 0 oder darunter, zusammen +0,03 bis +0,12 (Υ_d = 0,4–0,6). Die volle Feldform f = 1 liegt zusammen außerhalb 95 %.
- Das positive f hängt an IC 2574; ohne sie ist f zusammen +0,01 (95 %: −0,07…0,12), die radiale Form.
- Die über die Sternschicht gemittelten Werte f = ¼, ⅓ und ⟨ν⟩ (wirksam 0,26) hängen nicht vom Profil ab; sie passen die Primärgruppe, nicht die Hubble-Gruppe, zusammen nur mit IC 2574. Mittelebene (f = 0) und eine dünnere Gasschicht (f = 0,03–0,14 [A]) liegen näher an den Daten.
- a₀ steigt bei festem Υ mit f (von f = 0 bis 1 um 16–27 %) und fällt bei freiem Υ; Υ verschiebt a₀ stärker (0,96–1,46). Mit realistischen Fehlern trennen sich Primär- und Hubble-Fluss-Abstände in a₀ nur um 3–13 % statt 39 %.
