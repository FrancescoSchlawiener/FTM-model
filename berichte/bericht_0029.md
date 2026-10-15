# Bericht 29

2026-10-01T07:27:53 · Commit ed693ea · a₀ = 1,06·10⁻¹⁰ m/s², H₀ = 67,4 km/s/Mpc

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
| Rotation bei 5 kpc (Regel) | 250 | 249,25 | km/s | ABWEICHUNG |
| Bahnradius bleibt 5 kpc | 5 | 5 | kpc | OK |
| Rotation bei 8,72 kpc (Regel) | 216 | 215,22 | km/s | ABWEICHUNG |
| Bahnradius bleibt 8,72 kpc | 8,72 | 8,72 | kpc | OK |
| Rotation bei 15 kpc (Regel) | 198 | 196,21 | km/s | ABWEICHUNG |
| Bahnradius bleibt 15 kpc | 15 | 15 | kpc | OK |
| Rotation bei 50 kpc (Regel) | 180 | 178,07 | km/s | ABWEICHUNG |
| Bahnradius bleibt 50 kpc | 50 | 50 | kpc | OK |
| Newton bei 50 kpc (ohne Regel, fällt wie Kepler) | 71,841 | 71,841 | km/s | OK |
| tiefer Wert v⁴ = GMa₀ | 172 | 170,45 | km/s | ABWEICHUNG |

## Labor: RC-Kette mit Rückführung

Leiterkette (R = 1 kΩ, R·G = 0,01) mit Rückführung des Abgriffs 20 auf den Eingang, als Knotennetz gelöst; die Spannung am Eingang ist die Verstärkung ν. Die Daten tragen statt der Zeit die Knotennummer.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| Verstärkung ν = 1/(1 − e^{−s}), s = n·arcosh(1 + RG/2) | 1,1567 | 1,1567 |  | OK |
| Abfall je Glied e^{−arcosh(1 + RG/2)} | 0,90488 | 0,90488 |  | OK |

## Kosmos: Atmen mit Rückstellung und Fluchtmaß

Dunkle Energie als Rückstellung mal Rückkehrzahl des Ganzen, ρ_DE = ν_G·ρ_s. Lesart L (r*/R = 0,770 heute): ν_G = 2,787, verlangt sind nach Planck 3,026; H₀ = 65,55 km/s/Mpc, Weltalter 14,49 Mrd. J., Ω_DE = 0,669, w = −1,34 heute und −1,23 bei z = 1 (die Hülle wächst bei fester Faltung, ν_G steigt). Lesart I (r*/Λ = 0,971, zeitfest): ν_G = 2,218, H₀ = 60,91 km/s/Mpc, Weltalter 14,42 Mrd. J., Ω_DE = 0,616, w = −1. Gesehen, wie die Messung liest (Abstand zur letzten Streufläche, jede Schale in ihrer früheren Epoche): Lesart L wirkt wie ein fester Anteil 2,016 statt 3,026, der Abstand ist 2,5 % länger als nach Planck; Lesart I 2,218, 1,9 %. Triade des Ganzen (b₁ = 3 Windungsklassen, fester Anteil 3): H₀ = 67,20 km/s/Mpc, Weltalter 13,83 Mrd. J., Ω_DE = 0,685, Abstand zur letzten Streufläche 0,1 % gegen Planck. Die Daten tragen a(t) der Lesarten L, I, T als x-Koordinate.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| ν_G Lesart L: Quadratur gegen Monte-Carlo | 2,7618 | 2,7873 |  | OK |
| ν_G Lesart I: Quadratur gegen Monte-Carlo | 2,2128 | 2,2182 |  | OK |
| Weltalter L: Integral gegen Lauf | 14,488 | 14,488 | Mrd. J. | OK |
| Weltalter I: Integral gegen Lauf | 14,422 | 14,422 | Mrd. J. | OK |
| Gesehen: fester Anteil von Lesart I gibt sich selbst zurück | 2,2182 | 2,2182 |  | OK |
| H₀ Triade aus dem Kern gegen den Lauf | 67,202 | 67,202 | km/s/Mpc | OK |
| H₀ L von Hand 100·√(ω_m + ω_r + ν_G·ω_s) | 65,55 | 65,55 | km/s/Mpc | OK |

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
| 03 | r* = 8,72 kpc | 8,8824 | ABWEICHUNG |
| 03 | 172,0 km/s | 170,45 | ABWEICHUNG |
| 03 | Λ = 13,2 Gpc | 13,739 | ABWEICHUNG |
| 03 | 86,4 Mrd. Jahre | 89,621 | ABWEICHUNG |
| 03 | L = 27,6 Gpc | 27,934 | ABWEICHUNG |
| 03 | = 5,36 | 5,3643 | OK |
| 03 | q = 0,843 | 0,84287 | OK |
| 03 | Θ = 2πa₀/(cH₀) = 1,055 | 1,0171 | ABWEICHUNG |
| 03 | Ω_s = 0,244 | 0,22688 | ABWEICHUNG |
| 03 | H₀·√Ω_s = 33,3 km/s/Mpc | 32,104 | ABWEICHUNG |
| 03 | k_s = π³/30 = 1,034 | 1,0335 | OK |
| 03 | Ω_DE/Ω_s = 2,81 | 3,0263 | ABWEICHUNG |
| 03 | a₀ = 1,06·10⁻¹⁰ | 1,0646 | OK |
| 02 | Grenzmasse 1,456 M☉ | 1,4563 | OK |
| 03 | ν = 6,25 frei | 6,1439 | ABWEICHUNG |
| 03 | 5,36 senkrecht | 5,3097 | ABWEICHUNG |
| 03 | 4,59 gleichgerichtet | 4,5511 | ABWEICHUNG |
| 03 | R_Q = 51,0 kpc | 51,963 | ABWEICHUNG |
| 03 | a_N = 0,029·a₀ | 0,02922 | OK |
| 03 | χ = Q·(−ln q) = 0,917 | 0,91697 | OK |
| 03 | r*/Λ = 0,971 | 0,97096 | OK |
| 06 | ν = 1,1567 | 1,1567 | OK |
| 06 | s ≈ 3·10⁵ | 3,0422 | OK |
| 08 | Sonne 4,43 und 8,86 km | 4,4299 | OK |
| 08 | 8,86 km | 8,8598 | OK |
| 03 | nach Planck a₀ = 1,065·10⁻¹⁰ | 1,0646 | OK |
| 03 | SPARC-Mittel 0,9 % von ω_DE | 0,85424 | OK |
| 03 | etwa ±17 % | 16,976 | OK |
| 03 | (SPARC) H₀ = 67,2 km/s/Mpc | 67,202 | OK |
| 03 | und Ω_DE = 0,685, mit | 0,68466 | OK |
| 03 | 1,1·10⁻¹⁰ H₀ = 68,9 | 68,948 | OK |
| 03 | 68,9 und 0,700 | 0,70043 | OK |
| 03 | erst bei z = −0,08 ein | -0,10209 | ABWEICHUNG |
| 03 | = 0,47·v_flach | 0,4714 | OK |
| 03 | unter v_flach = 17–21 | 16,971 | OK |
| 03 | 17–21 km/s. In SPARC | 21,213 | OK |
| 03 | bei z_eq = 3402 | 3402 | OK |
| 03 | bei 3–4·10¹² | 3,1173 | OK |
| 03 | 3–4·10¹² M☉ (mit | 3,8339 | OK |
| 03 | Bereiche ab 1,4·10¹¹ | 1,4395 | OK |
| 03 | bis 1,1·10¹³ M☉ | 1,166 | ABWEICHUNG |
| 03 | M_min = 0,6–1,4·10⁷ M☉ | 0,58961 | OK |
| 03 | 1,4·10⁷ M☉; kaltes | 1,4395 | OK |
| 03 | zum Endwert 1,23 | 1,2328 | OK |
| 03 | drei Teile tun es bei z = 0,33 | 0,29501 | ABWEICHUNG |
| 03 | mit drei Teilen in 57,7 km/s/Mpc | 55,606 | ABWEICHUNG |
| 03 | 9,2 km/s bei 10⁷ M☉ | 9,1295 | ABWEICHUNG |
| 03 | 19,5 km/s bei 2·10⁸ M☉ | 19,307 | ABWEICHUNG |
| 03 | 10¹² M☉/Mpc: 2,5 kpc | 2,6299 | ABWEICHUNG |
| 03 | (ν = 1,007) | 1,0068 | OK |
| 03 | (ν = 3,3) | 3,3035 | OK |
| 03 | a₀/(πG) = 251 M☉/pc² | 242,07 | ABWEICHUNG |
| 03 | a_N = 2,54·a₀ | 2,5396 | OK |
| 03 | beträgt 0,648·a₀ | 0,64761 | OK |
| 07 | Rosetten mit 254,6° | 254,56 | OK |
| 07 | e^{−7342} | 7479,6 | ABWEICHUNG |
| 08 | bei der Milchstraße 0,43″ | 0,41894 | ABWEICHUNG |
| 08 | ε_L = ℓ_P/ℓ₀ = 2,0·10⁻⁶² | 1,9062 | OK |
| 08 | ε_L^{1/3}·m_Pl = 5,9·10⁻²⁹ kg | 5,8139 | ABWEICHUNG |
| 08 | c³τ_D/G = 1,1·10⁵⁴ kg | 1,1418 | OK |
| 08 | Δs(M_U) = ⅔·Λ | 2 | OK |
| 08 | (4π⁴/15)·ε_L² = 1,0·10⁻¹²² | 0,94388 | OK |
| 02 | r* = 32 fm | 32,453 | OK |
| 02 | a_B = 1661·r* | 1630,6 | ABWEICHUNG |
| 04 | Sonne 1,05·10⁷⁷ k | 1,0489 | OK |
| 04 | 0,022 s | 0,021968 | OK |
| 04 | Sonnenoberfläche 67 s | 66,981 | OK |
| 04 | auf der Erdbahn sind es 0,31 s | 0,31149 | OK |
| 04 | Takt 0,81 | 0,8096 | OK |
| 08 | 1,7512″ | 1,7512 | OK |
| 08 | 0,8756″ | 0,8756 | OK |
| 08 | 42,98″ | 42,981 | OK |
| 08 | r* = 7342 AE | 7479,6 | ABWEICHUNG |

**Ergebnis:** 31 Abweichung(en) · Bahndaten in `berichte/bericht_0029/` · Laufzeit 60 s
