# Schwelle: Modellwert und Messwert

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Das Modell setzt a₀ = cH₀/2π“ → im Modell ist a₀ = c²/(2Λ), heute ≈ cH₀/2π; „a₀ = cH₀/2π ist damit nicht bestätigt“ → an SPARC mit realistischen Fehlern 1,04–1,19·10⁻¹⁰ m/s², im Bereich von cH₀/2π; ohne Kreispaare a₀ ≤ 1,13·10⁻¹⁰. Die Notiz vergleicht cH₀/2π = 1,04–1,13·10⁻¹⁰ mit dem Literaturwert 1,2·10⁻¹⁰.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/schwelle_modell_mess.py` (Ausgabe `schwelle_modell_mess.txt`)
> Beobachtungswerte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.

## Ausgangslage

Im Modell ist a₀ = c²/(2Λ) mit der Schließlänge Λ der Rückkehr, heute ≈ cH₀/2π = c·kT_GH/ħ. [korr.] Gemessen ist a₀ = 1,20 ± 0,02 (stat) ± 0,24 (sys) ·10⁻¹⁰ m/s² [G]. Die Messung passt die Rotationskurven von rund 150 Scheibengalaxien mit **derselben** Funktion ν = 1/(1−e^{−√(g_bar/a₀)}) an, die im Modell die Rückkopplungssumme ist. Eine andere Übergangsfunktion ist also nicht die Ursache des Unterschieds. Die systematische Unsicherheit stammt fast ganz aus dem Masse-zu-Licht-Verhältnis der Sterne, Υ[3,6 µm] = 0,5 M☉/L☉ [G].

## 1. Abstand

| H₀ [km/s/Mpc] | cH₀/2π [10⁻¹⁰ m/s²] | Abstand in σ_sys |
|---|---|---|
| 67,4 | 1,042 | 0,66 |
| 70,0 | 1,082 | 0,49 |
| 73,0 | 1,129 | 0,30 |

a₀ = 1,20·10⁻¹⁰ m/s² verlangt H₀ = 77,6 km/s/Mpc. Der Unterschied liegt innerhalb der systematischen Unsicherheit der Messung.

## 2. Welches Masse-zu-Licht-Verhältnis schließt die Lücke?

Die Rechnung skaliert den Sternanteil der baryonischen Beschleunigung mit einem Faktor k, das Gas bleibt fest. Dann wird mit k = 1 angepasst, über g_bar = 10⁻¹² bis 10⁻⁹ m/s². Gesucht ist das k, bei dem das Modell-a₀ als 1,20·10⁻¹⁰ erscheint.

| H₀ | Sternanteil 100 % | 70 % | 50 % |
|---|---|---|---|
| 67,4 | k = 1,090 → Υ = 0,545 | 1,129 → 0,565 | 1,180 → 0,590 |
| 70,0 | 1,066 → 0,533 | 1,094 → 0,547 | 1,131 → 0,566 |

Die Anpassung skaliert nicht wie a₀ ∝ 1/Υ. Das gilt nur im tiefen Bereich. Im Newton-Bereich bestimmt Υ die Kurve selbst, und die Anpassung verteilt den Unterschied.

**Ergebnis:** Der Modellwert ist gleichbedeutend mit Υ[3,6 µm] ≈ 0,53–0,59 statt 0,50. Das liegt im Bereich der Sternpopulationsmodelle (etwa 0,4–0,7 [G]).

## 3. Das äußere Feld wirkt entgegen

Scheibengalaxien stehen im Feld ihrer Umgebung, und der Schalter χ hängt am Gesamtfeld. Die Rechnung erzeugt Kurven mit äußerem Feld, g = ν(√(g_N² + g_e²)/a₀)·g_N, und passt sie ohne äußeres Feld an (g_N = 10⁻¹² bis 10⁻⁹ m/s²):

| g_e | angepasstes a₀ / wahres a₀ |
|---|---|
| 0 | 1,000 |
| 0,01 a₀ | 0,94 |
| 0,03 a₀ | 0,76–0,77 |
| 0,1 a₀ | 0,43–0,44 |

Eine Anpassung ohne äußeres Feld liefert ein zu **kleines** a₀. Bei typischen Umgebungsfeldern von einigen Hundertstel a₀ [G] liegt das wahre a₀ also 5–30 % **über** dem angepassten Wert 1,2·10⁻¹⁰ m/s². Der Abstand zum Modellwert cH₀/2π (1,04–1,08·10⁻¹⁰) wächst damit von 10–13 % auf bis zu etwa 40 %. Die Näherung gewichtet alle Radien gleich; eine Anpassung an einen Katalog mit dem Umgebungsfeld jeder einzelnen Galaxie steht aus.

## 4. Ergebnis

- Ohne äußeres Feld liegt cH₀/2π 10–13 % unter dem gemessenen a₀. Das entspricht Υ[3,6 µm] ≈ 0,53–0,59 statt 0,50.
- Das äußere Feld vergrößert den Abstand. Ein höheres Υ und das äußere Feld wirken gegeneinander, und ihre Summe ist nicht gerechnet.
- a₀ = 1,20·10⁻¹⁰ verlangt H₀ = 77,6 km/s/Mpc; mit äußerem Feld läge der verlangte Wert noch höher.
- Gegen den Literaturwert 1,20·10⁻¹⁰ allein wäre a₀ ≈ cH₀/2π nicht bestätigt. An den SPARC-Rotationskurven mit den Fehlern von Abstand, Neigung und Masse-zu-Licht liegt der gelesene Wert aber bei 1,04–1,19·10⁻¹⁰ m/s² (Υ_d = 0,5), im Bereich von cH₀/2π; ohne Kreispaare (L ≥ 0,985·d_LSS) ist a₀ ≤ 1,13·10⁻¹⁰. Geprüft wird der Wert weiter an Υ und am Umgebungsfeld der Galaxien. [korr.]
