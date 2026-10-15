# Drehmoment am Warp der Milchstraße

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: nur Bezeichnungen, Drehimpuls des Rings L → J_R, spezifischer Drehimpuls j_⊥, j_c → h_⊥, h_c (j ist die Spinzahl); inhaltlich keine Korrektur nötig, [G]-Marken vollständig.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/warp_drehmoment.py` (Ausgabe `warp_drehmoment.txt`)
> Messwert: rechtläufige Knotendrehung +10,9 km/s/kpc bei R = 12–16 kpc [G]. Werte mit [G] sind aus dem Gedächtnis; Massen, Abstände, Einfallraten und die Neigung θ_w sind Annahmen.

## Ansatz

- **Feld:** baryonische Milchstraße (Massenmodell McMillan 2017 [G]), M_b = 6,68·10¹⁰ M☉, kugelgemittelt, Feldform ν = 1/(1 − e^{−√y}), a₀ = 1,2·10⁻¹⁰ m/s². Bei 14 kpc: g_N/a₀ = 0,36, ν = 2,23, v = 202 km/s (Newton 136), Ω = 14,5 km/s/kpc.
- **Störer:** Punktmasse; die Feldgleichung wird um den Kugelhintergrund (mit äußerem Feld 0,02 a₀) linearisiert, radial verstärkt mit ν − 1 + yν′, tangential mit ν − 1, alle geraden Multipole bis zur Konvergenz. Proben: konstantes ν − 1 gibt u = (ν − 1)f auf 10⁻³; die Multipolsumme trifft das exakte Newton-Drehmoment.
- **Größe:** größte Knotenrate ω = τ_Knoten/(J_R sin θ_w) eines Rings bei 14 kpc, θ_w = 0,1 (Warp-Höhe ≈ 1,4 kpc), über alle Azimute des Störers. Ihr Vorzeichen hängt von der Lage des Störers zur Knotenlinie ab.
- **Hierarchie:** Glättungslängen im Feld der Milchstraße ℓ = d(M/3M_b)^{1/3}: LMC 12,6 kpc, Sgr 2,4 kpc, M31 619 kpc. Der Ring liegt außerhalb dieser Zonen, also im kollektiven, verstärkten Feld; die Hierarchie ändert die Antwort nicht.

## Zahlen

Knotenrate in km/s/kpc, größter rechtläufiger / rückläufiger Wert:

| Quelle | Annahme [G] | Newton | Feld × ν(Scheibe) | volle Feldform | Feldverstärkung |
|---|---|---|---|---|---|
| Sagittarius-Zwerg | 5·10⁸ M☉ baryonisch, 18 kpc, Breite −20° | +0,41 / −0,32 | +0,61 / −0,48 | +0,69 / −0,50 | ×2,6 |
| Sagittarius-Zwerg | 1·10⁸ M☉ | +0,08 / −0,06 | +0,12 / −0,10 | +0,14 / −0,10 | ×2,6 |
| Große Magellansche Wolke | 3,2·10⁹ M☉, 50 kpc, −33° | +0,075 / −0,084 | +0,11 / −0,13 | +0,33 / −0,36 | ×6,4 |
| M31 | 1·10¹¹ M☉, 780 kpc, −22° | +4·10⁻⁴ / −6·10⁻⁴ | – | +0,017 / −0,026 (obere Schranke) | ×66 |
| äußeres Feld der Lokalen Gruppe | g_e = 0,01 / 0,02 / 0,03 a₀ | 0 | – | 0,002 / 0,009 / 0,020 | ∝ g_e² |
| Verwindung a_J = 4,26 pc (j = 1332) | – | +3,0·10⁻⁶ | +6,6·10⁻⁶ | +1,5·10⁻⁵ (Strom und Feld × ν) | rechtläufig |

Die Verstärkung entsteht, weil LMC und M31 im schwachen Feld sitzen (am Ort der LMC g_N ≈ 0,03 a₀); am Ort der Scheibe allein (× ν = 2,2) wäre sie kleiner. Für +10,9 km/s/kpc bräuchte die LMC baryonisch 1,1·10¹¹ M☉ (1,6 M_b), der Sagittarius-Zwerg 7,9·10⁹ M☉.

**Schräg einfallendes Gas.** Drehimpulszufuhr dJ_R/dt = Ṁ·h_⊥; mit Gas, das Kreisbahn-Drehimpuls h_c unter dem Winkel ψ bringt und ganz in der Außenscheibe landet (R > 12 kpc, M_o = 9,2·10⁹ M☉):

```
ω = (Ṁ/M_o)·sin ψ / sin θ_w
```

| Ṁ [M☉/Jahr, Annahme] | θ_w = 0,05 | θ_w = 0,1 |
|---|---|---|
| 1 (ψ = 90°) | 2,1 | 1,1 |
| 3 (ψ = 90°) | 6,4 | 3,2 |

+10,9 km/s/kpc verlangt Ṁ·sin ψ = 5,1 (θ_w = 0,05) bzw. 10,2 M☉/Jahr (θ_w = 0,1). Das Vorzeichen setzt die Richtung von h_⊥ gegen die Knotenlinie; die Drehung endet nach höchstens 90°, wenn die Scheibe zu h_⊥ gekippt ist (bei 5–11 km/s/kpc nach 0,14–0,3 Gyr). Das Modell ändert hier nur h_c (× 1,49): Stammt h_⊥ aus Kreisbahnen, kürzt sich das; ist h_⊥ fest, sinkt ω um diesen Faktor.

## Ergebnis

- Kein Gravitationsdrehmoment der Feldform erreicht den Messwert. Es fehlen bei der LMC ein Faktor 30, beim Sagittarius-Zwerg (5·10⁸ M☉) ein Faktor 16, bei M31 ein Faktor 400, beim äußeren Feld ein Faktor 500. Die Verstärkung der Feldform hebt die Gezeiten von LMC (× 6,4), Sgr (× 2,6) und M31 (× 66), schließt die Lücke aber nicht. Das Vorzeichen dieser Momente setzt die Lage des Störers, nicht das Modell.
- Ein festes rechtläufiges Vorzeichen hat nur die Verwindung a_J; sie ist verstärkt wie unverstärkt 7·10⁵- bis 4·10⁶-fach zu klein.
- Die Größenordnung erreicht nur die Drehimpulszufuhr durch schräg einfallendes Gas, bei Ṁ·sin ψ ≈ 5–10 M☉/Jahr in die Außenscheibe. Ihr Vorzeichen setzt die Richtung des Einfalls; aus dem Modell folgt weder diese Rate noch diese Richtung.
- **Gegen das Modell geprüft:** kein Widerspruch. Die Aussage „äußeres Feld, Magellansche Wolke und Verwindung reichen nicht; es fehlt ein äußeres Drehmoment“ gilt, erweitert um Sagittarius-Zwerg und M31. Das fehlende Moment kann nach dieser Rechnung nur die Zufuhr von Drehimpuls sein, nicht eine Kraft des verstärkten Feldes.
