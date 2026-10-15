# Der Wardenclyffe-Turm im Faltungsmodell

**Theoretischer Aufbau nach Teslas Bau, mit den Formen des Modells gerechnet**

Francesco Silvano Schauer

---

Nikola Tesla baute den Turm von 1901 an in Shoreham auf Long Island, 1917 wurde er abgerissen. Er sollte Signale und Energie durch die Erde senden: Die Erde ist der Schwinger, der Turm sein Sender, Teslas verstärkender Sender schaukelt die Schwingung auf. Diese Ausarbeitung legt den Turm als vollständigen Aufbau aus, von der Erdung bis zur Kuppel, und rechnet jede Größe mit den Formen des Faltungsmodells (FTM). Drei Formen tragen ihn:

```
ν = Σqⁿ = 1/(1 − q)                 Resonanz ist die Summe der Rückkehren
f_n = √(n(n + 1))·c/(2πR)           die Erdschale schließt in Kugelfunktionen
E = Q_el/(4πε₀r²)                   das Feld ist Windung, nicht Faltung: 1/r², ohne Schwelle
```

Alle Zahlen kommen aus simulation/turm.py mit den Konstanten des Kerns.

**Eingaben**

- Teslas Baumaße: Turm 57 m hoch (187 Fuß), Kuppel 20,7 m Durchmesser (68 Fuß), Schacht 36,6 m (120 Fuß), darunter 16 Eisenrohre von 91,4 m (300 Fuß)
- Auslegung (gewählt): Einspeisung 200 kW, Wendel von 6 m Durchmesser aus Kupferrohr von 5 cm, Ganghöhe zwei Rohrdurchmesser, Primärkreis mit Güte 50 und 50 kV
- Stoffe: Kupfer 5,96·10⁷ S/m, Boden 100 Ωm, Durchschlag der Luft 3 MV/m (trocken, auf Meereshöhe)

**Aufbau von unten nach oben**

- Erdung: 16 Eisenrohre im Schacht, 91,4 m in den Boden
- Primärkreis: Kondensatorbank und flache Spule am Fuß des Turms, vom Kraftwerk gespeist
- Wendel: eine einlagige Spule aus Kupferrohr, vom Fuß bis zur Kuppel im Turm
- Kuppel: die Endkapazität oben, glatt und mit großem Radius

## Modul A: Die Kuppel

```
C = 4πε₀r,     V_max = E_Luft·r,     Q = C·V,     W = ½·C·V²
```

Als Kugel gerechnet hat die Kuppel 1,15 nF. Bis zum Durchschlag der Luft trüge sie 31 MV, die Ladung 0,036 C und die Energie 0,56 MJ. Ihr Feld ist Windung, nicht Faltung: Coulomb bleibt 1/r², eine Schwelle gibt es nicht. Der einzige Hebel ist der große glatte Radius. Er senkt das Feld an der Oberfläche bei gleicher Spannung, und genau das war Teslas Form.

## Modul B: Die Wendel im Turm

```
L = μ₀πr²N²·K_N/h,     f = 1/(2π·√(L·(C_K + C_W))),     f_¼ = c/(4·NπD)
```

Ohne Wendel schwingen Turm und Schacht (93,6 m) als Viertelwelle bei 801 kHz. Die Wendel senkt das. Sie füllt die Turmhöhe mit 570 Windungen, so dicht, wie die Ganghöhe erlaubt. K_N ist die Nagaoka-Zahl der Stromlage, exakt aus den elliptischen Integralen. Wheelers Näherung 1/(1 + 0,9r/h) weicht davon um 0,2 % ab. Die Wendel hat 0,194 H, ihre Eigenkapazität ist 742 pF, mit der Kuppel zusammen 1,89 nF. Ihr Draht ist 10,7 km lang. Mit der Kuppel schwingt sie bei 8,31 kHz, als Viertelwelle des Drahts schlösse die stehende Welle bei 6,98 kHz. Beide liegen nahe beieinander: Die Welle schließt auf dem Draht, Fuß an der Erde, Spitze offen an der Kuppel.

## Modul C: Primärkreis und Einrasten

Der Primärkreis ist auf die Wendel abgestimmt, Δω = 0: Beide Zyklen rasten ein, die Energie schwingt zwischen ihnen hin und her. Die Kopplung muss im Fenster liegen. Zu schwach, und die Schwingung zerfällt vor der Übergabe, zu stark, und die Moden spalten auf. Die Grenze ist die kritische Kopplung k = 1/√(Q_p·Q_s) = 0,0018, die Übergabe dauert dann 1/(2kf) = 33 ms. Im Primärkreis laufen 192 J um, bei 50 kV verlangt das 153 nF und 2,39 mH.

## Modul D: Verluste und Güte

Alle Verluste liegen im transversalen Kanal oder gehen als Abgabe nach außen, die Windung selbst verliert nichts:

- Joule-Wärme im Kupferrohr: Der Strom dringt nur 0,72 mm tief ein, die Hauttiefe. Das gibt 1,60 Ω, ohne die Nähewirkung der Windungen.
- Abstrahlung des oben belasteten Turms, 160π²(h/λ)²: 3,9 mΩ.
- Erdung, 16 Rohre ρ/(2πℓ)·(ln(4ℓ/a) − 1): 86 mΩ als untere Grenze.

Die Güte der Wendel ist Q_s = ωL/R = 5969, eine obere Grenze. Die Rückkehr hebt die Amplitude damit um ν = Q_s/π = 1900: Das ist Teslas verstärkender Sender als Summe der Rückkehren, dieselbe Form wie die Regel und die Kette des SFFR.

## Modul E: Erdung und Gegenplatte

Die Erde hat als Kugel 709 µF. Auf dem geschlossenen Raum ist ΣQ_el = 0: Was die Kuppel abgibt, nimmt die Erde auf und gibt es im nächsten Halbtakt zurück. Im Betrieb trägt die Kuppel 5,7 mC, das verschiebt das Potential der Erde um 8,0 V. Der Turm pumpt keine Ladung in die Erde, er schiebt sie hin und her, mit 344 A durch die Rohre.

## Modul F: Der Betriebspunkt

Bei 200 kW fließen in der Wendel 344 A. Die Kuppel schwingt mit 4,91 MV Scheitelspannung, zwei Wege geben dasselbe: über den Strom, V = √2·I/(ωC), und über die gespeicherte Energie, W = Q_s·P/ω = 22,9 kJ, V = √(2W/C). Das Feld an der Kuppel ist 475 kV/m, 16 % des Durchschlags. Von den 200 kW werden 189 kW im Kupfer warm, 10,1 kW in der Erdung, 465 W strahlt der Turm ab.

| Windungen | Frequenz | Güte | Kuppel | abgestrahlt |
|---|---|---|---|---|
| 100 | 47,4 kHz | 2003 | 1,19 MV | 28,9 kW |
| 200 | 23,7 kHz | 3322 | 2,17 MV | 6,0 kW |
| 400 | 11,8 kHz | 4935 | 3,74 MV | 1,11 kW |
| 570 | 8,31 kHz | 5969 | 4,91 MV | 0,47 kW |

Je dichter die Wendel, desto tiefer die Frequenz, desto höher die Kuppel und desto weniger strahlt der Turm ab. Teslas Ziel war der Strom durch die Erde, nicht die Abstrahlung: Darum gilt die dichteste Wendel. Die Bodenwelle hat dann in 100 km 2,0 mV/m, in 1000 km 0,20 mV/m.

## Modul G: Die Erdschale

Eine Welle läuft in 0,1335 s einmal um die Erde. Je Umlauf geschlossen gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz. Als stehende Wellen auf der Kugelschale schließen die Moden in Kugelfunktionen: 10,59, 18,34, 25,94 und 33,49 Hz, mit 2n + 1 Moden je Grad. Das ist dieselbe Form wie die Schwellen N(N + 1) der Atomhülle und die Schließung je Umlauf des Elektrons, 2π·a_B = λ: Die Erdschale zählt wie eine Atomhülle. Verluste in Boden und Ionosphäre senken die Kugelmoden zur Zählung je Umlauf hin.

Zwei Zyklen rasten ein, wenn |Δω| ≤ 2k·|cos φ|. Der Turm schwingt 784-mal schneller als die tiefste Kugelmode, die Kopplung müsste halb so schnell sein wie seine eigene Schwingung. Damit die Kuppel selbst bei 10,59 Hz schwänge, bräuchte sie eine Spule von 1,96·10⁵ H. Der Turm rastet an der Erdschale nicht ein, er treibt sie nur.

Die geschlossene Schale behält, was sie bekommt, sie vermehrt es nicht: Die Rückkehr hebt die umlaufende Amplitude, abgeben kann die Schale im Gleichgewicht je Umlauf nur, was der Turm einspeist (fester Haushalt der Wiederkehr, Kac). Bei 200 kW kommen auf der Erdoberfläche 3,9·10⁻¹⁰ W/m² an. Für 1 W/m² überall müsste der Turm 5,1·10¹⁴ W einspeisen.

## Modul H: Die Regel am Turm

An der Erdoberfläche ist g = 9,82 m/s², also s = √(g/a₀) = 3,04·10⁵ und der Rückkehranteil der Faltung q = e^{−s}, lg q = −1,3·10⁵. Die Faltung kehrt nur weit unter a₀ zurück, und das Feld des Turms trägt keine Regel, es ist Phase, nicht Faltung. Aus dem Medium holt der Turm nichts, seine Energie ist die seines Kraftwerks.

## Schluss

Teslas Bild trifft die Form: Die Erde mit ihrer Schale ist ein geschlossener Schwinger, ihre Moden zählen wie die einer Atomhülle, und Resonanz ist die Summe der Rückkehren. Als Aufbau trägt der Turm: Mit 200 kW schwingt die Kuppel bei 8,31 kHz mit 4,91 MV, weit unter dem Durchschlag. Er scheitert am Haushalt: Die Rückkehr hebt Amplituden, nicht Leistung, und der Turm schwingt 784-mal zu schnell, um an der Erdschale einzurasten. Als Sender von Signalen trägt er, als Kraftwerk für alle nicht.

## Optimierung der Bauart

**Fragestellung.** Die Module A bis H rechnen Teslas Maße. Hier ist die Form frei, Material und Leistung bleiben gleich: dasselbe Kupferrohr (10,7 km, 5 cm, Ganghöhe zwei Rohrdurchmesser) und 200 kW. Zielgröße ist die Güte Q_s der Wendel. An ihr hängen alle Leistungen des Baus: die Rückkehr ν = Q_s/π, die gespeicherte Energie W = Q_s·P/ω und das Kopplungsfenster U = k·Q_s einer Übertragung.

**Methode.** Die Wendel hat den Durchmesser D und die Höhe h. Bei gleichem Rohr und gleicher Ganghöhe ist die Drahtlänge πDh/(g·d), das Produkt D·h bleibt fest bei 342 m². Die Kuppel ist so weit wie die Wendel und schirmt deren obere Windungen. Die Induktivität folgt exakt aus der Nagaoka-Zahl. Die Kopplung zweier gleicher Bauten Seite an Seite folgt aus dem Neumann-Integral über beide Stromlagen:

```
M = μ₀/(4π)·∮∮ dl₁·dl₂/|r₁ − r₂|,     k = M/L,     U = k·Q_s,     η = U²/(1 + √(1 + U²))²
```

Ein zweiter Weg rechnet den Fluss des Feldes jedes Rings, geschlossen aus den elliptischen Integralen, durch die Windungen des anderen Baus. Beide Wege stimmen auf 12 Stellen. Fern gibt der Dipol Seite an Seite, k = r²(h + 0,9r)/(4d³), in 10 km dasselbe auf 0,2 %. Nahe überschätzt er die Kopplung.

**Güte des Kupfers.** Mit N = h/(g·d) wird die Güte des Leiters

```
Q_D = ωL/R_D = π·K_N·D/(2g·δ)
```

unabhängig von der Zahl der Windungen. Eine weitere Wendel hebt D, senkt aber die Frequenz, die Hauttiefe δ wächst, und K_N fällt. Dazwischen liegt ein flaches Maximum. Die geschlossene Form stimmt mit der vollen Rechnung auf 0,7 %.

| D/h | D | h | Windungen | Frequenz | Güte | 50 % bis |
|---|---|---|---|---|---|---|
| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |
| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |
| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |
| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |
| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |

Das Maximum liegt bei D/h = 1,01, eine Wendel so weit wie hoch. Die Reichweite einer Übertragung wächst mit flacherer Form weiter, weil k mit r² geht, die Güte nicht.

**Ergebnis.** Drei Bauten aus derselben Leiterlänge:

| | Teslas Maße | Kupfer, D = h | supraleitend |
|---|---|---|---|
| Durchmesser der Wendel | 6 m | 18,6 m | 13,6 m |
| Höhe der Wendel | 57 m | 18,4 m | 25,1 m |
| Windungen | 570 | 184 | 251 |
| Kuppel | 20,7 m | 18,6 m | 13,6 m |
| Erdrohre | 16 | 16 | 64 |
| Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz |
| Güte Q_s | 5969 | 10765 | 7,33·10⁵ |
| Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ |
| Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag |
| Kuppel bei 200 kW | 4,91 MV | 8,06 MV | über dem Durchschlag |
| Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW |
| Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ |
| 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m |

Mit demselben Kupfer gibt die Form D = h das 1,80-fache der Güte und das 2,69-fache der gespeicherten Energie. Die Kuppel schwingt bei 200 kW mit 8,06 MV statt 4,91 MV, abgestrahlt werden nur noch 26 W.

**Supraleitende Fassung.** Supraleitend fällt der Widerstand des Leiters weg: Das Einrasten ist vollständig, der transversale Kanal geschlossen. Es bleiben Abstrahlung und Erdung, die Erdung überwiegt. 64 Erdrohre geben 21 mΩ als untere Grenze ohne gegenseitige Abschirmung. Die Güte ist darum eine obere Grenze, ohne die Wechselverluste des Supraleiters und ohne Kühlung. Bei gleicher Leiterlänge liegt ihr Maximum bei D/h = 0,54, sie ist 68-mal größer als in der Form aus Kupfer.

Die Grenze ist nun die Kuppel. Sie erreicht den Durchschlag mit 17,6 kW, dann laufen 895 A in der Wendel und 0,30 MJ schwingen. Die Decke ist W_max = ½·C·(E_Luft·r)²: Eine doppelt so weite Kuppel trägt 1,82 MJ und erreicht den Durchschlag bei 107 kW. Die Amplitude läuft in 34 s auf ihren Endwert zu, die Resonanz ist 9,5 mHz breit (Kupfer, D = h: 0,61 s und 0,52 Hz). Das hält nur ein Antrieb aus Halbleitern, der in der Phase der Wendel einrastet, ψ̇ = Δω − K·sin ψ.

**Größe.** Wachsen alle Längen außer dem Rohr mit dem Faktor s, so gilt L ∝ s³, C ∝ s, f ∝ s⁻² und δ ∝ s. Die Güte des Kupfers Q_D = π·K_N·D/(2gδ) bleibt gleich. Die gespeicherte Energie bei gleicher Leistung wächst mit s², die Reichweite mit s. Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m.

**Grenzen.** Keine Form rastet an der Erdschale ein: Die Wendel aus Kupfer schwingt 526-mal, die supraleitende 655-mal schneller als die tiefste Kugelmode. Der Haushalt bleibt, die Schale gibt ab, was eingespeist wird. Für Signale ist die Wendel der falsche Bau. Als Viertelwelle ohne Wendel steht der Strahlungswiderstand von 36,5 Ω gegen 86 mΩ der Erdung, 99,8 % der Leistung gehen als Welle hinaus. In 1000 km kommen 4,2 mV/m an, das 21-fache der Wendel in Teslas Maßen.

**Folgerung.** Bei gleichem Material entscheidet die Form der Wendel: so weit wie hoch. Sie hebt die Rückkehr um das 1,80-fache und die gespeicherte Energie um das 2,69-fache. Supraleitend wächst die Güte um das 68-fache, dann begrenzt die Kuppel. Die Übertragung an einen gleichen, eingerasteten Bau reicht von 54 m über 130 m bis 454 m und wächst mit der Größe. Energie für alle bleibt aus: Die Rückkehr hebt Amplituden, nicht Leistung.
