# Schwelle aus der Anpassung der Feldform

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Feldform als wahre Bahnregel als nicht tragende Voraussetzung gekennzeichnet (Bahnen folgen der örtlichen Regel; SPARC f = +0,03 bis +0,13, Feldform außerhalb 95 %); Ergebnis und prüfbare Folgen (−4 % je Dekade Σ₀, Residuen bis −0,12 dex) auf die Feldform beschränkt, gemessen bei x = 1–3 −0,025 ± 0,018 statt −0,095; Umrechnung „H₀ = 70, a₀/(2πG) = 118 M☉/pc²“ → H₀ = 67,4, 119,0 M☉/pc² (mit H₀ = 70: 123,6); die Beispiele R_d/r* gehören zu H₀ = 67,4. Die Notiz rechnet mit a₀ = cH₀/2π (1,04·10⁻¹⁰ bei H₀ = 67,4) als wahrem Wert.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/schwelle_anpassung.py` (Ausgabe `schwelle_anpassung.txt`)
> Feldform ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N], ν = 1/(1 − e^{−s}), s = √(|g_N|/a₀), wahres a₀ = cH₀/2π; Finite Volumen in (R, z), 420 × 420 Zellen bis 80 r*. Angepasst wird mit a = ν(k·g_R/a₀)·k·g_R, g_R radiales Newton-Feld in der Ebene, k = Υ_angenommen/Υ_wahr. [A] Annahme, [G] Messwert aus dem Gedächtnis.
> [korr.] Diese Rechnung setzt die Feldform als Bahnregel voraus; das trägt nicht: Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N, die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Für Kugeln sind beide gleich, die Zahlen bleiben dort.

## Löser

Plummer-Kugel (Feldform = algebraische Relation exakt): Abweichung ≤ 10⁻⁴ bis 25 r*. Kuzmin-Scheibe (exakt a = ν(|g_N|)·g_R): 0,05–0,2 % bei h = 3 r*, 0,1–0,9 % bei h = 0,3 r*.

## Die Scheibe in der Feldform

Dünne Exponentialscheibe, Verhältnis Feldform zu algebraischer Relation ν(g_R)·g_R in der Ebene:

| R_d/r* | 2πGΣ₀/a₀ | R/R_d = 0,25 | 0,5 | 1 | 2 | 3 | 5 | 8 |
|---|---|---|---|---|---|---|---|---|
| 0,2 | 25 | 0,907 | 0,940 | 0,960 | 0,979 | 0,995 | 1,013 | 1,012 |
| 0,3 | 11 | 0,820 | 0,878 | 0,922 | 0,966 | 0,993 | 1,016 | 1,013 |
| 0,7 | 2,0 | 0,653 | 0,753 | 0,850 | 0,943 | 0,989 | 1,020 | 1,015 |
| 1 | 1,0 | 0,606 | 0,717 | 0,829 | 0,936 | 0,988 | 1,021 | 1,016 |
| 3 | 0,11 | 0,529 | 0,656 | 0,792 | 0,925 | 0,987 | 1,023 | 1,017 |

Innen liegt die Feldform darunter (bei R_d/4 0,91 bis 0,53), Gleichstand bei 3,5 R_d, außen (5–8 R_d) 1,2–2,3 % darüber. Im Tiefbereich gilt für die Kuzmin-Scheibe a_Feld/a_alg = √(cos θ). Eine dicke Scheibe (sech², z₀ = 0,1 R_d) wird in der Mittelebene innen etwas weniger gesenkt (R_d = 3 r*, R_d/2: 0,72 statt 0,66), außen gleich.

Umrechnung [A]: mit H₀ = 67,4 ist a₀/(2πG) = 119,0 M☉/pc² (mit H₀ = 70: 123,6) [korr.], also Σ₀ = 119·(r*/R_d)² M☉/pc² (die Tabellen rechnen mit 118). Beispiele (H₀ = 67,4) [korr.]: 6·10¹⁰ M☉ mit R_d = 2,6 kpc → R_d/r* = 0,29; 10⁹ M☉ mit 1,5 kpc → 1,3; 10⁸ M☉ mit 1 kpc → 2,7.

## Anpassung ohne äußeres Feld

Stichprobe [A]: 8 dünne Scheiben R_d/r* = 0,2 … 3 (logarithmisch gleich verteilt, Σ₀ = 3000 … 13 M☉/pc²), dazu 5 mit Plummer-Bulge (b = 0,15 R_d, B/T = 0,2–0,4) bei R_d = 0,2–0,45 r*. Punkte bei 0,5–8 R_d im Abstand 0,5 R_d, Fehler 5 % in v. „Fest“: Υ wahr, nur a₀ frei. „Vorgabe“: Υ und a₀ je Galaxie frei, Υ mit Vorgabe 0,1 dex wie in echten Anpassungen.

| R_d/r* | Σ₀ [M☉/pc²] | a₀ fest | Vorgabe: k | Vorgabe: a₀ |
|---|---|---|---|---|
| 0,2 | 2950 | 1,014 | 0,95 | 1,17 |
| 0,45 | 580 | 0,965 | 0,79 | 1,41 |
| 1 | 118 | 0,932 | 0,75 | 1,33 |
| 3 | 13 | 0,912 | 0,85 | 1,09 |
| mit Bulge | – | 0,997–1,021 | 0,91–0,99 | 1,05–1,17 |

| alle Punkte gemeinsam | Υ fest | Υ mit Vorgabe |
|---|---|---|
| 13 Galaxien, 0,5–8 R_d | 0,955 | 1,17 |
| nur Punkte 0,5–5 / 1–8 / 2–8 R_d | 0,883 / 0,984 / 1,015 | 1,11 / 1,15 / 1,09 |
| 8 Scheiben dünn / dick | 0,937 / 0,945 | 1,26 / 1,22 |
| 8 Scheiben mit Gasscheibe (R_g = 2,5 R_d, Gasanteil 0,1–0,7 [A]) | 0,863 | 1,14 |

Ohne Vorgabe für Υ ist die Anpassung bei Scheiben niedriger Flächendichte entartet (Υ → 0, nur k·a₀ bestimmt); global gibt sie 1,28.

Mit festem Υ ergibt die Feldform in m/s²: 0,995 / 1,033 / 1,078·10⁻¹⁰ für H₀ = 67,4 / 70 / 73. Soll das angepasste a₀ bei festem Υ = 1,20·10⁻¹⁰ sein, muss Υ_wahr = 1,10–1,13·Υ_angenommen sein (H₀ = 67,4–70; algebraisch 1,07–1,09), mit Gasscheibe 1,27–1,32 (algebraisch 1,10–1,15). Mit Υ_angenommen[3,6 µm] = 0,50 [G] heißt das Υ_wahr ≈ 0,55–0,56, mit Gas 0,64–0,66.

## Äußeres Feld

Angepasst ohne äußeres Feld, 13 Galaxien, global (senkrecht / gleichgerichtet):

| g_e/a₀ | Feldform, Υ fest | Feldform, Υ mit Vorgabe | algebraisch, Υ fest |
|---|---|---|---|
| 0,01 | 0,81 / 0,71 | 1,04 / 0,93 | 0,86 / 0,74 |
| 0,03 | 0,62 / 0,51 | 0,87 / 0,71 | 0,65 / 0,54 |
| 0,05 | 0,50 / 0,41 | 0,75 / 0,58 | 0,53 / 0,43 |

Scheiben hoher Flächendichte (R_d = 0,2 r*) bleiben bei 0,95–1,01, Scheiben niedriger Flächendichte (R_d = 3 r*) fallen auf 0,11–0,45: Ihre Außenpunkte liegen tief unter g_e.

## Ergebnis

[korr.] Alle Punkte gelten nur unter der Annahme, die Feldform sei die Bahnregel und a₀ = cH₀/2π der wahre Wert; die Bahnen folgen aber der örtlichen Regel.

- **Υ fest** (so ist der Messwert 1,20·10⁻¹⁰ gewonnen [G]): Die Feldform schließt den Abstand nicht, sie vergrößert ihn. Das angepasste a₀ liegt 4–6 % unter dem wahren (mit Gas 14 %); aus 10–13 % Abstand werden 14–17 % (mit Gas 22–25 %). Nur eine Anpassung allein der Außenpunkte (2–8 R_d) hebt a₀ um 1,5 %.
- **Υ frei mit Vorgabe:** Das angepasste a₀ liegt 14–26 % über dem wahren und schließt den Abstand ganz (1,22–1,27·10⁻¹⁰ bei H₀ = 67,4–70), aber nur, weil das angepasste Υ bis 25 % (mit Gas bis 44 %) unter dem wahren liegt. Mit äußerem Feld 0,01 a₀ bleiben +4 % (senkrecht) bzw. −7 % (gleichgerichtet), also teilweise; ab 0,03 a₀ liegt es 13–42 % darunter.
- **Wer trägt:** Bei festem Υ ziehen die Scheiben niedriger Flächendichte (Σ₀ ≲ 250 M☉/pc², R_d ≳ 0,7 r*) a₀ über ihre inneren Punkte nach unten; Scheiben hoher Flächendichte mit Bulge liegen bei 1,00–1,02. Bei freiem Υ heben die mittleren Scheiben (Σ₀ ≈ 100–1300 M☉/pc²) a₀ am stärksten (1,30–1,41), weil ihr Übergangsbereich in der Scheibe liegt. Das äußere Feld wirkt fast nur über die Scheiben niedriger Flächendichte.
- **Masse-zu-Licht-Verhältnis:** Bei festem Υ ist cH₀/2π in der Feldform mit Υ_wahr = 1,10–1,13·Υ_angenommen verträglich (algebraisch 1,07–1,09), mit Gasscheibe 1,27–1,32 (algebraisch 1,10–1,15), jeweils ohne äußeres Feld.

## Prüfbare Folgen

[korr.] Diese Folgen hängen an der Feldform als Bahnregel; sie gelten nur für sie.

- Bei festem Υ fiele in der Feldform das je Galaxie angepasste a₀ mit der zentralen Flächendichte: 1,01 (Σ₀ ≈ 3000 M☉/pc²), 0,97 (580), 0,93 (120), 0,91 (13), etwa 4 % je Dekade. [korr.]
- In der Feldform hinge die Abweichung von der algebraischen Relation bei gleichem g_bar an 2πGΣ/g_R: +0,007 dex (unter 0,3), −0,016 (0,3–1), −0,07 (1–2), −0,12 dex (2–4). [korr.] Betroffen sind die inneren Punkte der Scheiben niedriger Flächendichte. Die Residuen der gemessenen Relation hängen nach [G] nicht von Flächenhelligkeit oder Radius ab; an SPARC nachgeprüft liegen die Daten bei 2πGΣ/g_R = 1–3 bei −0,025 ± 0,018 statt −0,095 der Feldform. [korr.]
- Bei freiem Υ liegt das angepasste Υ unter dem Wert der Sternpopulationen, am stärksten bei Σ₀ ≈ 100–600 M☉/pc² (k ≈ 0,75–0,8).
- Das angepasste a₀ hängt am inneren Rand der Kurve: 0,88 mit Punkten ab 0,5 R_d bis 5 R_d, 1,015 nur mit 2–8 R_d.
