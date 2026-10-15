# Bericht 5

2026-10-01T03:19:07 · Commit 50370bb · a₀ = 1,1·10⁻¹⁰ m/s², H₀ = 67,4 km/s/Mpc

## Atom: Wasserstoff

Elektron und Proton auf der Kreisbahn um den gemeinsamen Schwerpunkt, drei Umläufe, nur Coulomb: die Neigung der Faltung ist dort 2,3·10³⁹-mal schwächer als die Windung.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| Umlauffrequenz der Hülle | 6,5797·10¹⁵ | 6,5815·10¹⁵ | Hz | OK |
| Tempo der Hülle v/c = α | 0,0072974 | 0,0072993 |  | OK |
| Coulomb gegen Neigung | 2,27·10³⁹ | 2,2687·10³⁹ |  | OK |

## Sonnensystem

Fünf Planeten auf Kreisbahnen im Trichter der Sonne, zwölf Jahre; die Sonne ist ein Einzeltrichter (Newton-Bereich), die Eigenzeit der Planeten läuft mit dem Takt √(1 − 6Δs☉/r) mit.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| Erde: Umlaufzeit | 365,26 | 365,26 | Tage | OK |
| Erde: Tempo v = c·√(3Δs☉/r) | 29,785 | 29,785 | km/s | OK |
| Erde: Uhr im Sonnentrichter je Jahr | 0,31149 | 0,31152 | s | OK |
| Merkur: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉) | 0,24085 | 0,24085 | Jahre | OK |
| Venus: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉) | 0,6152 | 0,6152 | Jahre | OK |
| Erde: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉) | 1 | 1 | Jahre | OK |
| Mars: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉) | 1,8808 | 1,8808 | Jahre | OK |
| Jupiter: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉) | 11,867 | 11,867 | Jahre | OK |

## Merkur: Drehung des Perihels

Merkur auf seiner Ellipse, zehn Umläufe; Takt und Faltungsstrecke als Glied erster Ordnung. Gegenprobe ohne das Glied: das Perihel steht.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| Perihel mit Takt und Faltung je Jahrhundert | 42,981 | 42,981 | ″ | OK |
| Gegenprobe ohne Glied, Betrag unter | 0,01 | 3,9006·10⁻⁶ | ″ | OK |

## Galaxie: Rotationskurve

Vier Teststerne auf Kreisbahnen um eine Milchstraße als Punktmasse (6·10¹⁰ M☉), zwei Milliarden Jahre, mit der Rückkehrregel a = a_N/(1 − q); Gegenprobe nach Newton.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| Rotation bei 5 kpc (Regel) | 250 | 250,09 | km/s | OK |
| Bahnradius bleibt 5 kpc | 5 | 5 | kpc | OK |
| Rotation bei 8,72 kpc (Regel) | 216 | 216,38 | km/s | OK |
| Bahnradius bleibt 8,72 kpc | 8,72 | 8,72 | kpc | OK |
| Rotation bei 15 kpc (Regel) | 198 | 197,55 | km/s | OK |
| Bahnradius bleibt 15 kpc | 15 | 15 | kpc | OK |
| Rotation bei 50 kpc (Regel) | 180 | 179,59 | km/s | OK |
| Bahnradius bleibt 50 kpc | 50 | 50 | kpc | OK |
| Newton bei 50 kpc (ohne Regel, fällt wie Kepler) | 71,841 | 71,841 | km/s | OK |
| tiefer Wert v⁴ = GMa₀ | 172 | 172,03 | km/s | OK |

## Zahlen im Modelltext

| Baustein | Stelle | Kern | Urteil |
|---|---|---|---|
| 01 | Erde Δs = 1,478 mm | 1,4783 | OK |
| 01 | Sonne Δs = 492,2 m | 492,21 | OK |
| 01 | Erde 8,87 mm | 8,8701 | OK |
| 01 | Sonne 2,95 km | 2,9533 | OK |
| 02 | 5,906·10⁻³⁹ | 5,9061 | OK |
| 02 | 1/137,036 | 137,04 | OK |
| 02 | 0,529·10⁻¹⁰ m | 0,52918 | OK |
| 02 | 6,58·10¹⁵ Hz | 6,5797 | OK |
| 02 | Erde 29,78 km/s | 29,785 | OK |
| 02 | 1,853 M☉ | 1,8533 | OK |
| 03 | Δs = 197,4 AE | 197,41 | OK |
| 03 | ν = 1,582 | 1,582 | OK |
| 03 | r* = 8,72 kpc | 8,7194 | OK |
| 03 | 172,0 km/s | 172,03 | OK |
| 03 | Λ = 13,2 Gpc | 13,239 | OK |
| 03 | 86,4 Mrd. Jahre | 86,362 | OK |
| 03 | L = 27,6 Gpc | 27,592 | OK |
| 03 | = 5,36 | 5,3643 | OK |
| 03 | q = 0,843 | 0,84287 | OK |
| 03 | Θ = 2πa₀/(cH₀) = 1,055 | 1,0555 | OK |
| 03 | Ω_s = 0,244 | 0,24433 | OK |
| 03 | H₀·√Ω_s = 33,3 km/s/Mpc | 33,315 | OK |
| 03 | k_s = π³/30 = 1,034 | 1,0335 | OK |
| 03 | Ω_DE/Ω_s = 2,81 | 2,8102 | OK |
| 03 | a₀ = 1,06·10⁻¹⁰ | 1,0646 | OK |
| 04 | 0,022 s | 0,021968 | OK |
| 04 | Sonnenoberfläche 67 s | 66,981 | OK |
| 04 | auf der Erdbahn sind es 0,31 s | 0,31149 | OK |
| 04 | Takt 0,81 | 0,8096 | OK |
| 08 | 1,7512″ | 1,7512 | OK |
| 08 | 0,8756″ | 0,8756 | OK |
| 08 | 42,98″ | 42,981 | OK |
| 08 | r* = 7342 AE | 7342,3 | OK |

**Ergebnis:** alles OK · Bahndaten in `berichte/bericht_0005/` · Laufzeit 53 s
