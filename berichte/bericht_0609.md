# Bericht 609

2026-10-25T05:01:29 · Commit a081658 · a₀ = 1,06·10⁻¹⁰ m/s², H₀ = 67,4 km/s/Mpc

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
| Rotation bei 5 kpc (Regel) | 249,25 | 249,25 | km/s | OK |
| Bahnradius bleibt 5 kpc | 5 | 5 | kpc | OK |
| Rotation bei 8,88 kpc (Regel) | 214,41 | 214,4 | km/s | OK |
| Bahnradius bleibt 8,88 kpc | 8,88 | 8,88 | kpc | OK |
| Rotation bei 15 kpc (Regel) | 196,21 | 196,21 | km/s | OK |
| Bahnradius bleibt 15 kpc | 15 | 15 | kpc | OK |
| Rotation bei 50 kpc (Regel) | 178,07 | 178,07 | km/s | OK |
| Bahnradius bleibt 50 kpc | 50 | 50 | kpc | OK |
| Newton bei 50 kpc (ohne Regel, fällt wie Kepler) | 71,841 | 71,841 | km/s | OK |
| tiefer Wert v⁴ = GMa₀ | 170,45 | 170,45 | km/s | OK |

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

## Daten: SPARC und AME2020

SPARC, 141 Galaxien, 2688 Punkte, Abstand, Neigung und Masse-zu-Licht je Galaxie mitgeschätzt, Masse-zu-Licht der Scheiben so, dass Stern und Gas derselben Regel folgen (Υ_d = 0,67): a₀ = 1,11·10⁻¹⁰ m/s² (16–84 %: 1,00–1,21); χ² Bose 2917,6, einfach 2970,7, Standard 2958,4; Dimension der Rückkehr 0,92 ± 0,05; Rest gegen Gasanteil innerhalb der Galaxien 0,00 ± 0,04 dex; je Galaxie gegen Masse ρ = −0,06, gegen Gasanteil ρ = 0,07, gegen Hubble-Typ ρ = −0,03, gegen Flächenhelligkeit ρ = 0,03. AME2020: kleinste Masse je Nukleon A = 56, größte Bindung A = 62, Spaltbarkeit am Maximum der Tröpfchenformel 0,248. Die Daten tragen log g_bar, log g_obs je Punkt.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| a₀ aus SPARC gegen a₀ aus dem Atmen (Toleranz 1,5 σ) | 1,0646 | 1,1063 | 10⁻¹⁰ m/s² | OK |
| Dimension der Rückkehr (Toleranz 2 σ) | 1 | 0,91682 |  | OK |
| Bose-Form vor Standard-Form | 1 | 1 |  | OK |
| kleinste Masse je Nukleon bei A | 56 | 56 |  | OK |
| größte Bindung je Nukleon bei A | 62 | 62 |  | OK |
| Spaltbarkeit am Maximum der Bindung | 0,25 | 0,24768 |  | OK |
| Bindung des Wasserstoffs ist die halbe Spreizung α²/2 (PubChem) | 2,6611·10⁻⁵ | 2,6611·10⁻⁵ |  | OK |
| Rest je Galaxie unabhängig von Masse (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Gasanteil (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Typ (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Flächenhelligkeit (p > 0,05) | 1 | 1 |  | OK |

## Zahlen im Modelltext

| Baustein | Stelle | Kern | Urteil |
|---|---|---|---|
| A:wardenclyffe_de | Wheelers Näherung 1/(1 + 0,9r/h) weicht davon um 0,2 % ab | 0,20258 | OK |
| A:wardenclyffe_de | Kupferrohr (10,7 km, 5 cm | 10,743 | OK |
| A:wardenclyffe_de | bleibt fest bei 342 m² | 342 | OK |
| A:wardenclyffe_de | Beide Wege stimmen auf 12 Stellen | 0 | OK |
| A:wardenclyffe_de | in 10 km dasselbe auf 0,2 % | 0,18496 | OK |
| A:wardenclyffe_de | Die geschlossene Form stimmt mit der vollen Rechnung auf 0,7 % | 0,74421 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 10,129 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 33,764 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 337 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 8,4642 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 9395,7 | OK |
| A:wardenclyffe_de | | 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m | | 90,54 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 13,077 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 26,153 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 261 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 7,1538 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 10280 | OK |
| A:wardenclyffe_de | | 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m | | 106,19 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 18,493 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 18,493 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 184 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 5,6214 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 10737 | OK |
| A:wardenclyffe_de | | 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m | | 129,49 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 26,153 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 13,077 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 130 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 4,434 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 10208 | OK |
| A:wardenclyffe_de | | 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m | | 157,13 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 32,031 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 10,677 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 106 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 3,8919 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 9508,1 | OK |
| A:wardenclyffe_de | | 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m | | 175,97 | OK |
| A:wardenclyffe_de | Das Maximum liegt bei D/h = 1,01 | 1,0099 | OK |
| A:wardenclyffe_de | | Durchmesser der Wendel | 6 m | 18,6 m | 13,6 m | | 18,585 | OK |
| A:wardenclyffe_de | | Durchmesser der Wendel | 6 m | 18,6 m | 13,6 m | | 13,623 | OK |
| A:wardenclyffe_de | | Höhe der Wendel | 57 m | 18,4 m | 25,1 m | | 18,402 | OK |
| A:wardenclyffe_de | | Höhe der Wendel | 57 m | 18,4 m | 25,1 m | | 25,105 | OK |
| A:wardenclyffe_de | | Windungen | 570 | 184 | 251 | | 184 | OK |
| A:wardenclyffe_de | | Windungen | 570 | 184 | 251 | | 251 | OK |
| A:wardenclyffe_de | | Kuppel | 20,7 m | 18,6 m | 13,6 m | | 18,585 | OK |
| A:wardenclyffe_de | | Kuppel | 20,7 m | 18,6 m | 13,6 m | | 13,623 | OK |
| A:wardenclyffe_de | | Erdrohre | 16 | 16 | 64 | | 64 | OK |
| A:wardenclyffe_de | | Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz | | 8,3082 | OK |
| A:wardenclyffe_de | | Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz | | 5,5742 | OK |
| A:wardenclyffe_de | | Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz | | 6,9409 | OK |
| A:wardenclyffe_de | | Güte Q_s | 5969 | 10765 | 7,33·10⁵ | | 5969 | OK |
| A:wardenclyffe_de | | Güte Q_s | 5969 | 10765 | 7,33·10⁵ | | 10765 | OK |
| A:wardenclyffe_de | | Güte Q_s | 5969 | 10765 | 7,33·10⁵ | | 7,3307 | OK |
| A:wardenclyffe_de | | Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ | | 1900 | OK |
| A:wardenclyffe_de | | Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ | | 3426,6 | OK |
| A:wardenclyffe_de | | Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ | | 2,3334 | OK |
| A:wardenclyffe_de | | Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag | | 22,869 | OK |
| A:wardenclyffe_de | | Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag | | 61,472 | OK |
| A:wardenclyffe_de | | Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag | | 1 | OK |
| A:wardenclyffe_de | | Kuppel bei 200 kW | 4,91 MV | 8,06 MV | über dem Durchschlag | | 4,9142 | OK |
| A:wardenclyffe_de | | Kuppel bei 200 kW | 4,91 MV | 8,06 MV | über dem Durchschlag | | 8,0567 | OK |
| A:wardenclyffe_de | | Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW | | 7,9845 | OK |
| A:wardenclyffe_de | | Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW | | 2,3945 | OK |
| A:wardenclyffe_de | | Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW | | 17,644 | OK |
| A:wardenclyffe_de | | Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ | | 0,91299 | OK |
| A:wardenclyffe_de | | Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ | | 0,73597 | OK |
| A:wardenclyffe_de | | Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ | | 0,29657 | OK |
| A:wardenclyffe_de | | 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m | | 54,264 | OK |
| A:wardenclyffe_de | | 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m | | 129,96 | OK |
| A:wardenclyffe_de | | 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m | | 453,96 | OK |
| A:wardenclyffe_de | das 1,80-fache der Güte und das 2,69-fache der gespeicherten Energie | 1,8035 | OK |
| A:wardenclyffe_de | das 1,80-fache der Güte und das 2,69-fache der gespeicherten Energie | 2,688 | OK |
| A:wardenclyffe_de | abgestrahlt werden nur noch 26 W | 26,406 | OK |
| A:wardenclyffe_de | 64 Erdrohre geben 21 mΩ | 21,487 | OK |
| A:wardenclyffe_de | liegt ihr Maximum bei D/h = 0,54 | 0,54265 | OK |
| A:wardenclyffe_de | sie ist 68-mal größer als in der Form aus Kupfer | 68,097 | OK |
| A:wardenclyffe_de | dann laufen 895 A in der Wendel | 895,12 | OK |
| A:wardenclyffe_de | Eine doppelt so weite Kuppel trägt 1,82 MJ und erreicht den Durchschlag bei 107 kW | 1,8192 | OK |
| A:wardenclyffe_de | Eine doppelt so weite Kuppel trägt 1,82 MJ und erreicht den Durchschlag bei 107 kW | 107,32 | OK |
| A:wardenclyffe_de | Die Amplitude läuft in 34 s auf ihren Endwert zu, die Resonanz ist 9,5 mHz breit | 33,618 | OK |
| A:wardenclyffe_de | Die Amplitude läuft in 34 s auf ihren Endwert zu, die Resonanz ist 9,5 mHz breit | 9,4683 | OK |
| A:wardenclyffe_de | (Kupfer, D = h: 0,61 s und 0,52 Hz) | 0,61472 | OK |
| A:wardenclyffe_de | (Kupfer, D = h: 0,61 s und 0,52 Hz) | 0,51781 | OK |
| A:wardenclyffe_de | Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m | 253,7 | OK |
| A:wardenclyffe_de | Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m | 262,66 | OK |
| A:wardenclyffe_de | Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m | 1,0312 | OK |
| A:wardenclyffe_de | Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m | 528,13 | OK |
| A:wardenclyffe_de | Die Güte des Kupfers Q_D = π·K_N·D/(2gδ) bleibt gleich | 1 | OK |
| A:wardenclyffe_de | Die Wendel aus Kupfer schwingt 526-mal, die supraleitende 655-mal schneller | 526,31 | OK |
| A:wardenclyffe_de | Die Wendel aus Kupfer schwingt 526-mal, die supraleitende 655-mal schneller | 655,34 | OK |
| A:wardenclyffe_de | gegen 86 mΩ der Erdung, 99,8 % | 85,947 | OK |
| A:wardenclyffe_de | gegen 86 mΩ der Erdung, 99,8 % | 99,765 | OK |
| A:wardenclyffe_de | In 1000 km kommen 4,2 mV/m an, das 21-fache der Wendel | 4,2377 | OK |
| A:wardenclyffe_de | In 1000 km kommen 4,2 mV/m an, das 21-fache der Wendel | 20,713 | OK |
| A:wardenclyffe_de | Supraleitend wächst die Güte um das 68-fache | 68,097 | OK |
| A:wardenclyffe_de | reicht von 54 m über 130 m bis 454 m | 453,96 | OK |
| A:wardenclyffe_de | hat die Kuppel 1,15 nF | 1,1516 | OK |
| A:wardenclyffe_de | trüge sie 31 MV | 31,05 | OK |
| A:wardenclyffe_de | die Ladung 0,036 C | 0,035757 | OK |
| A:wardenclyffe_de | die Energie 0,56 MJ | 0,55513 | OK |
| A:wardenclyffe_de | Turm und Schacht (93,6 m) | 93,6 | OK |
| A:wardenclyffe_de | als Viertelwelle bei 801 kHz | 800,73 | OK |
| A:wardenclyffe_de | mit 570 Windungen | 570 | OK |
| A:wardenclyffe_de | Die Wendel hat 0,194 H | 0,19376 | OK |
| A:wardenclyffe_de | ihre Eigenkapazität ist 742 pF | 742,38 | OK |
| A:wardenclyffe_de | mit der Kuppel zusammen 1,89 nF | 1,894 | OK |
| A:wardenclyffe_de | Ihr Draht ist 10,7 km lang | 10,744 | OK |
| A:wardenclyffe_de | schwingt sie bei 8,31 kHz | 8,3082 | OK |
| A:wardenclyffe_de | schlösse die stehende Welle bei 6,98 kHz | 6,9757 | OK |
| A:wardenclyffe_de | k = 1/√(Q_p·Q_s) = 0,0018 | 0,0018305 | OK |
| A:wardenclyffe_de | 1/(2kf) = 33 ms | 32,878 | OK |
| A:wardenclyffe_de | Im Primärkreis laufen 192 J um | 191,56 | OK |
| A:wardenclyffe_de | 153 nF und 2,39 mH | 153,25 | OK |
| A:wardenclyffe_de | 153 nF und 2,39 mH | 2,3946 | OK |
| A:wardenclyffe_de | nur 0,72 mm tief ein | 0,71523 | OK |
| A:wardenclyffe_de | Das gibt 1,60 Ω | 1,6046 | OK |
| A:wardenclyffe_de | 160π²(h/λ)²: 3,9 mΩ | 3,9404 | OK |
| A:wardenclyffe_de | 86 mΩ als untere Grenze | 85,947 | OK |
| A:wardenclyffe_de | Q_s = ωL/R = 5969 | 5969 | OK |
| A:wardenclyffe_de | ν = Q_s/π = 1900 | 1900 | OK |
| A:wardenclyffe_de | Die Erde hat als Kugel 709 µF | 708,87 | OK |
| A:wardenclyffe_de | trägt die Kuppel 5,7 mC | 5,6592 | OK |
| A:wardenclyffe_de | um 8,0 V | 7,9834 | OK |
| A:wardenclyffe_de | mit 344 A durch die Rohre | 343,56 | OK |
| A:wardenclyffe_de | fließen in der Wendel 344 A | 343,56 | OK |
| A:wardenclyffe_de | mit 4,91 MV Scheitelspannung | 4,9142 | OK |
| A:wardenclyffe_de | mit 4,91 MV Scheitelspannung | 4,9142 | OK |
| A:wardenclyffe_de | W = Q_s·P/ω = 22,9 kJ | 22,869 | OK |
| A:wardenclyffe_de | Das Feld an der Kuppel ist 475 kV/m, 16 % | 474,8 | OK |
| A:wardenclyffe_de | Das Feld an der Kuppel ist 475 kV/m, 16 % | 15,827 | OK |
| A:wardenclyffe_de | 189 kW im Kupfer warm | 189,39 | OK |
| A:wardenclyffe_de | 10,1 kW in der Erdung | 10,144 | OK |
| A:wardenclyffe_de | 465 W strahlt der Turm ab | 465,08 | OK |
| A:wardenclyffe_de | | 100 | 47,4 kHz |  | 47,357 | OK |
| A:wardenclyffe_de | | 100 | 47,4 kHz |  | 2002,6 | OK |
| A:wardenclyffe_de | | 100 | 47,4 kHz |  | 1,1922 | OK |
| A:wardenclyffe_de | | 100 | 47,4 kHz |  | 28,897 | OK |
| A:wardenclyffe_de | | 200 | 23,7 kHz |  | 23,678 | OK |
| A:wardenclyffe_de | | 200 | 23,7 kHz |  | 3321,6 | OK |
| A:wardenclyffe_de | | 200 | 23,7 kHz |  | 2,1715 | OK |
| A:wardenclyffe_de | | 200 | 23,7 kHz |  | 5,9912 | OK |
| A:wardenclyffe_de | | 400 | 11,8 kHz |  | 11,839 | OK |
| A:wardenclyffe_de | | 400 | 11,8 kHz |  | 4935,5 | OK |
| A:wardenclyffe_de | | 400 | 11,8 kHz |  | 3,7433 | OK |
| A:wardenclyffe_de | | 400 | 11,8 kHz |  | 1,1128 | OK |
| A:wardenclyffe_de | | 570 | 8,31 kHz |  | 8,3082 | OK |
| A:wardenclyffe_de | | 570 | 8,31 kHz |  | 5969 | OK |
| A:wardenclyffe_de | | 570 | 8,31 kHz |  | 4,9142 | OK |
| A:wardenclyffe_de | | 570 | 8,31 kHz |  | 0,46508 | OK |
| A:wardenclyffe_de | in 100 km 2,0 mV/m | 2,0459 | OK |
| A:wardenclyffe_de | in 1000 km 0,20 mV/m | 0,20459 | OK |
| A:wardenclyffe_de | in 0,1335 s einmal um die Erde | 0,13353 | OK |
| A:wardenclyffe_de | gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz | 7,4892 | OK |
| A:wardenclyffe_de | gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz | 14,978 | OK |
| A:wardenclyffe_de | gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz | 22,467 | OK |
| A:wardenclyffe_de | 10,59, 18,34, 25,94 und 33,49 Hz | 10,591 | OK |
| A:wardenclyffe_de | 10,59, 18,34, 25,94 und 33,49 Hz | 18,345 | OK |
| A:wardenclyffe_de | 10,59, 18,34, 25,94 und 33,49 Hz | 25,943 | OK |
| A:wardenclyffe_de | 10,59, 18,34, 25,94 und 33,49 Hz | 33,493 | OK |
| A:wardenclyffe_de | Der Turm schwingt 784-mal schneller | 784,43 | OK |
| A:wardenclyffe_de | die Kopplung müsste halb so schnell sein | 0,99873 | OK |
| A:wardenclyffe_de | bei 10,59 Hz schwänge, bräuchte sie eine Spule von 1,96·10⁵ H | 1,9609 | OK |
| A:wardenclyffe_de | Bei 200 kW kommen auf der Erdoberfläche 3,9·10⁻¹⁰ W/m² an | 3,9211 | OK |
| A:wardenclyffe_de | müsste der Turm 5,1·10¹⁴ W einspeisen | 5,1006 | OK |
| A:wardenclyffe_de | ist g = 9,82 m/s² | 9,8203 | OK |
| A:wardenclyffe_de | s = √(g/a₀) = 3,04·10⁵ | 3,0437 | OK |
| A:wardenclyffe_de | lg q = −1,3·10⁵ | -1,3219 | OK |
| A:wardenclyffe_de | schwingt die Kuppel bei 8,31 kHz mit 4,91 MV | 8,3082 | OK |
| A:wardenclyffe_de | schwingt die Kuppel bei 8,31 kHz mit 4,91 MV | 4,9142 | OK |
| A:wardenclyffe_de | der Turm schwingt 784-mal zu schnell | 784,43 | OK |
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
| 03 | r* = 8,88 kpc | 8,8824 | OK |
| 03 | 170,4 km/s | 170,45 | OK |
| 03 | Λ = 13,7 Gpc | 13,739 | OK |
| 03 | 89,6 Mrd. Jahre | 89,621 | OK |
| 03 | L = 27,9 Gpc | 27,934 | OK |
| 03 | = 5,36 | 5,3643 | OK |
| 03 | q = 0,843 | 0,84287 | OK |
| 03 | Θ = 2πa₀/(cH₀) = 1,020 | 1,0201 | OK |
| 03 | Ω_s = 0,228 | 0,22822 | OK |
| 03 | H₀·√Ω_s = 32,1 km/s/Mpc | 32,104 | OK |
| 03 | k_s = π³/30 = 1,034 | 1,0335 | OK |
| 03 | Ω_DE/Ω_s = 3,03 | 3,0263 | OK |
| 03 | a₀ = 1,06·10⁻¹⁰ | 1,0646 | OK |
| 02 | ω₃ = 2,018 der Beiwert der Form n = 3 | 2,0182 | OK |
| 02 | m_Pl³/(μ_e·m_u)² = 1,456 M☉ | 1,4563 | OK |
| 02 | im Tal nach der Tröpfchenformel Z = 26,7, zwischen ⁵⁶Fe und ⁶²Ni | 1 | OK |
| 02 | höchstens kT = 0,096·(μ/μ_e)·m_ec² | 0,096377 | OK |
| 02 | für Kohlenstoff und Sauerstoff 5,0·10⁸ K | 5,0007 | OK |
| 03 | an den blauen, die noch Sterne bilden, χ² = 21 in 15 Bins | 20,615 | OK |
| F:lesungen | Sechs gasreiche (Mancera Piña | 6 | OK |
| F:lesungen | äußerer Punkt bei a_N = 0,018 bis 0,039·a₀ | 0,017693 | OK |
| F:lesungen | äußerer Punkt bei a_N = 0,018 bis 0,039·a₀ | 0,039216 | OK |
| F:lesungen | drehen wie Newton allein, im Mittel +0,01 dex | 0,0078052 | OK |
| F:lesungen | Die Regel verlangt 0,40 dex mehr Geschwindigkeit | 0,40217 | OK |
| F:lesungen | einen 2,5-mal kleineren Sinus der Neigung | 2,5245 | OK |
| 03 | die Zelle erreicht die letzte Streufläche bei Θ = 1,020, für jedes a₀ | 1,0199 | OK |
| 03 | Mit der Strahlung ist es Θ = 1,024 bis 1,030 (a₀ = 0,9 bis 1,2) | 1,0243 | OK |
| 03 | Mit der Strahlung ist es Θ = 1,024 bis 1,030 (a₀ = 0,9 bis 1,2) | 1,0298 | OK |
| 03 | das Mittel aus Neigung (parametrisch) und Schwelle (transversal), a = √(a_N·a₀) | 1,0005 | OK |
| 03 | Ihr Produkt ist fest, v_flach², gleich sind beide bei r* | 1 | OK |
| 03 | Ihr Produkt ist fest, v_flach², gleich sind beide bei r* | 1 | OK |
| 03 | 1,012 (Zelle und innen, L/ℓ₀ = Θ^{2/3}/(π·Ω_m)^{1/3}) | 1,0121 | OK |
| 03 | 1,012 (Zelle und innen, L/ℓ₀ = Θ^{2/3}/(π·Ω_m)^{1/3}) | 1 | OK |
| 03 | und 1,023 (Zelle und Horizont, L·H/2πc = (π·Θ·Ω_m)^{−1/3}) | 1,0232 | OK |
| 03 | und 1,023 (Zelle und Horizont, L·H/2πc = (π·Θ·Ω_m)^{−1/3}) | 1 | OK |
| 03 | liegt das Ganze in diesem Fenster von 2 % | 2,3248 | OK |
| 03 | dauert das Fenster 0,68 Mrd. Jahre | 0,67735 | OK |
| 03 | Es öffnete sich vor 0,58 Mrd. Jahren (z = 0,041) | 0,58108 | OK |
| 03 | Es öffnete sich vor 0,58 Mrd. Jahren (z = 0,041) | 0,041154 | OK |
| 03 | und schließt sich in 0,10 Mrd. Jahren | 0,09627 | OK |
| 03 | Im besten Punkt ist L/D_S = 0,9984 ± 0,0007 | 0,99841 | OK |
| 03 | Im besten Punkt ist L/D_S = 0,9984 ± 0,0007 | 6,7837·10⁻⁴ | OK |
| 03 | also α = 3,2° (2,4° bis 3,9°) | 3,2281 | OK |
| 03 | also α = 3,2° (2,4° bis 3,9°) | 2,4424 | OK |
| 03 | also α = 3,2° (2,4° bis 3,9°) | 3,8572 | OK |
| 03 | (vier Zufallsfelder, M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉ | 0,9785 | OK |
| 03 | (vier Zufallsfelder, M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉ | 3,492 | OK |
| 03 | je über 31 000 Körper in R₂₀₀), liest das Profil | 31000 | OK |
| 03 | liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2) | 3,974 | OK |
| 03 | liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2) | 3,3021 | OK |
| 03 | liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2) | 5,1672 | OK |
| 03 | Mit ihm ist der Rest an X-COP im Kern −0,02 dex | -0,019045 | OK |
| 03 | mit der Konzentration aus dem Einfall Q = 5,6 ± 0,2 | 5,6275 | OK |
| 03 | mit der Konzentration aus dem Einfall Q = 5,6 ± 0,2 | 0,20065 | OK |
| F:lesungen | (SH0ES 2022, 3492 Gleichungen) liest H₀ = 73,0 ± 1,0 | 73,043 | OK |
| F:lesungen | (SH0ES 2022, 3492 Gleichungen) liest H₀ = 73,0 ± 1,0 | 1,0071 | OK |
| 03 | kreuzt die anderen zuerst in der Mitte, bei ω·t = π/(2√2) = 1,11 | 1,1107 | OK |
| 03 | kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11 | 1,1045 | OK |
| 03 | kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11 | 1,0783 | OK |
| 03 | kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11 | 1,1075 | OK |
| F:lesungen | liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36° | 31,387 | OK |
| F:lesungen | liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36° | 28,358 | OK |
| F:lesungen | liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36° | 36,3 | OK |
| 03 | HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5 | 9,7622 | OK |
| 03 | HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5 | 6,5401 | OK |
| 03 | X-COP liest c = 4,2 | 4,2194 | OK |
| F:lesungen | An KiDS-1000 (ξ₊ und ξ₋, 225 Werte | 225 | OK |
| F:lesungen | liegen Modell und ΛCDM gleich, χ² = 307,8 | 307,8 | OK |
| F:lesungen | liegen Modell und ΛCDM gleich, χ² = 307,8 | 307,79 | OK |
| F:lesungen | frei S8 = 0,73 statt 0,82 (Δχ² = 38) | 0,73147 | OK |
| F:lesungen | frei S8 = 0,73 statt 0,82 (Δχ² = 38) | 0,82106 | OK |
| F:lesungen | frei S8 = 0,73 statt 0,82 (Δχ² = 38) | 37,675 | OK |
| F:lesungen | Bis M₂₀₀ = 2,5·10¹⁴ M☉ liegt R₂₀₀ jenseits der Naht | 2,5239 | OK |
| 03 | liest HIFLUGCS Q = 5,5 ± 0,3 | 5,5168 | OK |
| 03 | liest HIFLUGCS Q = 5,5 ± 0,3 | 0,27882 | OK |
| 03 | (Rest im Kern −0,04 dex, in der Form des Gases +0,41) | -0,042279 | OK |
| 03 | (Rest im Kern −0,04 dex, in der Form des Gases +0,41) | 0,41074 | OK |
| 03 | HIFLUGCS c = 5,5 (χ² = 447 statt 601 | 5,4966 | OK |
| 03 | HIFLUGCS c = 5,5 (χ² = 447 statt 601 | 446,78 | OK |
| 03 | HIFLUGCS c = 5,5 (χ² = 447 statt 601 | 600,66 | OK |
| 03 | innerhalb 300 kpc −0,04 statt +0,21 dex | -0,043283 | OK |
| 03 | innerhalb 300 kpc −0,04 statt +0,21 dex | 0,21215 | OK |
| 03 | im Kern das 6,8-fache des Gases tragen | 6,8017 | OK |
| F:lesungen | das Maß des Ganzen 1 + Q, das dort dasselbe gibt, 0,39 dex | 0,39128 | OK |
| F:lesungen | (c = 4,0, R₂₀₀ um 45 kpc) | 3,974 | OK |
| F:lesungen | (c = 4,0, R₂₀₀ um 45 kpc) | 45,372 | OK |
| F:lesungen | und sie verlangt dort 0,13 dex mehr | 0,13217 | OK |
| F:lesungen | Die gasreichen isolierten Zwerge liegen dagegen 0,37 dex über dieser Form | 0,37049 | OK |
| 03 | χ² = 21 in 15 Bins | 15 | OK |
| 03 | mit 0,003·a₀ im Betrag schon 534 | 533,69 | OK |
| 03 | χ² = 27483 statt 503 | 27483 | OK |
| 03 | Naht (30 bis 66 kpc) | 30,379 | OK |
| 03 | Naht (30 bis 66 kpc) | 66,01 | OK |
| 03 | Die 22 klassischen Satelliten | 22 | OK |
| 03 | liegen ohne Wirtsfeld +0,15 ± 0,04 dex über der Regel | 0,14523 | OK |
| 03 | liegen ohne Wirtsfeld +0,15 ± 0,04 dex über der Regel | 0,0437 | OK |
| 03 | mit ihm im Betrag +0,46 ± 0,05 | 0,45901 | OK |
| 03 | mit ihm im Betrag +0,46 ± 0,05 | 0,048514 | OK |
| F:galaxien | gasreichen isolierten 0,00 ± 0,04 dex | -0,0011822 | OK |
| F:galaxien | gasreichen isolierten 0,00 ± 0,04 dex | 0,040615 | OK |
| F:galaxien | die gasarmen isolierten +0,13 ± 0,07 | 0,12991 | OK |
| F:galaxien | die gasarmen isolierten +0,13 ± 0,07 | 0,071752 | OK |
| F:galaxien | die 36 klassischen Satelliten +0,11 ± 0,03 | 36 | OK |
| F:galaxien | die 36 klassischen Satelliten +0,11 ± 0,03 | 0,1096 | OK |
| F:galaxien | die 36 klassischen Satelliten +0,11 ± 0,03 | 0,030478 | OK |
| F:galaxien | die ultraschwachen +0,37 dex | 0,3668 | OK |
| F:galaxien | liest Υ_V = 1,96 für die gasreichen | 1,9647 | OK |
| F:galaxien | 5,2 für die klassischen Satelliten | 5,1942 | OK |
| F:galaxien | weder das Verhältnis der Dichten zum Wirt | 0,024285 | OK |
| F:galaxien | noch das Perizentrum der Bahn | -0,0097672 | OK |
| F:galaxien | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 2,0906 | OK |
| F:galaxien | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 0,28449 | OK |
| F:galaxien | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 188 | OK |
| F:galaxien | liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49 | 0,011471 | OK |
| F:galaxien | liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49 | 0,028374 | OK |
| F:galaxien | liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49 | 2,4853 | OK |
| F:galaxien | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | 7 | OK |
| F:galaxien | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | -1,0305 | OK |
| F:galaxien | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | 0,19643 | OK |
| F:galaxien | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03 | 27 | OK |
| F:galaxien | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03 | -0,026844 | OK |
| F:galaxien | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03 | 0,026277 | OK |
| F:galaxien | die 27 ultraschwachen +0,07 ± 0,04 | 27 | OK |
| F:galaxien | die 27 ultraschwachen +0,07 ± 0,04 | 0,066168 | OK |
| F:galaxien | die 27 ultraschwachen +0,07 ± 0,04 | 0,037226 | OK |
| F:galaxien | die gasarmen isolierten +0,05 ± 0,08 | 0,048077 | OK |
| F:galaxien | die gasarmen isolierten +0,05 ± 0,08 | 0,07922 | OK |
| F:galaxien | ebenso −0,03 und +0,06 | -0,030422 | OK |
| F:galaxien | ebenso −0,03 und +0,06 | 0,057619 | OK |
| F:galaxien | (36 Satelliten der Milchstraße mit Eigenbewegung) | 36 | OK |
| 03 | X-COP (12 Haufen | 12 | OK |
| 03 | a_N = 0,036·a₀ (0,031 bis 0,065) | 0,036204 | OK |
| 03 | a_N = 0,036·a₀ (0,031 bis 0,065) | 0,031225 | OK |
| 03 | a_N = 0,036·a₀ (0,031 bis 0,065) | 0,06472 | OK |
| 03 | gegen 0,029·a₀. Dort | 0,02922 | OK |
| 03 | M/M_Gas = 7,9 (4,7 bis 10,5) | 7,8963 | OK |
| 03 | M/M_Gas = 7,9 (4,7 bis 10,5) | 4,7366 | OK |
| 03 | M/M_Gas = 7,9 (4,7 bis 10,5) | 10,467 | OK |
| 03 | gegen 1 + Q = 6,28, die Bilanz | 6,2832 | OK |
| 03 | Sternen von einem Viertel des Gases | 0,25673 | OK |
| 03 | in R₂₀₀ das 1,7-fache des Gases | 1,7423 | OK |
| 03 | (0,1·R₅₀₀) das 12-fache | 12,404 | OK |
| 03 | über alle 108 Punkte +0,03 dex auf den Daten (χ² = 601) | 0,033248 | OK |
| 03 | über alle 108 Punkte +0,03 dex auf den Daten (χ² = 601) | 600,66 | OK |
| 03 | über alle 108 Punkte +0,03 dex | 108 | OK |
| 03 | die örtliche Regel 0,17 dex darunter (χ² = 4173) | 0,16872 | OK |
| 03 | die örtliche Regel 0,17 dex darunter (χ² = 4173) | 4173,4 | OK |
| 03 | das Maß 1 + Q dort nur das 2,9-fache | 2,9059 | OK |
| 03 | R_Q = 52,0 kpc | 51,963 | OK |
| 03 | a_N = 0,029·a₀ | 0,02922 | OK |
| 03 | χ = Q·(−ln q) = 0,917 | 0,91697 | OK |
| 03 | r*/Λ = 0,971 | 0,97096 | OK |
| 09 | μ_B = 9,274·10⁻²⁴ J/T | 9,274 | OK |
| 09 | m_p/m_e = 1836,2 | 1836,2 | OK |
| 09 | beim Elektron 2,9·10⁴⁴ | 2,8542 | OK |
| 09 | m_Pl/√2 = 1,54·10⁻⁸ kg | 1,539 | OK |
| 09 | a_J = 3,27 m | 3,273 | OK |
| 09 | 3Δs = 4,4 mm | 4,435 | OK |
| 09 | also j = 738 | 737,99 | OK |
| F:galaxien | Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43) | 65 | OK |
| F:galaxien | Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43) | -0,10235 | OK |
| F:galaxien | Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43) | 0,36 | OK |
| F:galaxien | Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43) | 0,435 | OK |
| 07 | gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern) | 4,2603 | OK |
| 07 | gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern) | 28,976 | OK |
| 07 | gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern) | 7 | OK |
| 07 | in 3 Gyr auf das 2- bis 3-fache ihres Halbmassenradius | 1,9763 | OK |
| 07 | in 3 Gyr auf das 2- bis 3-fache ihres Halbmassenradius | 3,2937 | OK |
| 07 | ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton | 1,1738 | OK |
| 07 | ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton | 0,95915 | OK |
| 07 | ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton | 0,20831 | OK |
| 03 | Über die ganze Schicht gemittelt liest die Bahn das senkrechte Feld für jedes Profil gleichverteilt zwischen 0 und 2πGΣ | 0 | OK |
| 03 | innerhalb 2R_d im Median um 0,03 dex der Beschleunigung, 4 % des Tempos | 0,03161 | OK |
| 03 | innerhalb 2R_d im Median um 0,03 dex der Beschleunigung, 4 % des Tempos | 3,5738 | OK |
| 03 | 4 % des Tempos (flächenschwach 0,05, flächenhell 0,01) | 0,049872 | OK |
| 03 | 4 % des Tempos (flächenschwach 0,05, flächenhell 0,01) | 0,0092151 | OK |
| F:galaxien | zählte die Massenfunktion bei W50 = 50 bis 63 km/s das 218- und 508-fache der Galaxien | 218,46 | OK |
| F:galaxien | zählte die Massenfunktion bei W50 = 50 bis 63 km/s das 218- und 508-fache der Galaxien | 507,7 | OK |
| F:galaxien | bei 355 bis 447 km/s das 2,0- und 4,4-fache | 1,9975 | OK |
| F:galaxien | bei 355 bis 447 km/s das 2,0- und 4,4-fache | 4,4391 | OK |
| F:galaxien | der Linie einer Scheibe an der Schwelle mit der gelesenen Streuung von 12 km/s | 12,362 | OK |
| F:galaxien | unter W50 = 100 km/s auf das 0,4- bis 4,6-fache | 0,43286 | OK |
| F:galaxien | unter W50 = 100 km/s auf das 0,4- bis 4,6-fache | 4,5715 | OK |
| F:galaxien | bei 100 bis 250 km/s bleibt das 7- bis 21-fache | 6,8075 | OK |
| F:galaxien | bei 100 bis 250 km/s bleibt das 7- bis 21-fache | 20,878 | OK |
| F:galaxien | verlangte ALFALFA ε = 0,02 bis 0,05 | 0,017589 | OK |
| F:galaxien | verlangte ALFALFA ε = 0,02 bis 0,05 | 0,047717 | OK |
| F:galaxien | das 20- bis 56-fache ihrer Baryonen als nie abgekühltes Gas | 19,957 | OK |
| F:galaxien | das 20- bis 56-fache ihrer Baryonen als nie abgekühltes Gas | 55,853 | OK |
| F:galaxien | hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32 | 21,001 | OK |
| F:galaxien | hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32 | 340,4 | OK |
| F:galaxien | hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32 | 701,66 | OK |
| F:galaxien | hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32 | 31,62 | OK |
| F:galaxien | liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache | 177,17 | OK |
| F:galaxien | liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache | 100,51 | OK |
| F:galaxien | liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache | 304,41 | OK |
| 03 | In Bereichen schwerer als M liegt der Anteil erfc(δ_c/√(2S)) | 0 | OK |
| 03 | mit 0,90 bis 0,92 des Tempos in der Ebene | 0,90194 | OK |
| 03 | mit 0,90 bis 0,92 des Tempos in der Ebene | 0,91835 | OK |
| 03 | nach Newton mit den Baryonen allein mit 0,86 bis 0,88 | 0,85874 | OK |
| 03 | nach Newton mit den Baryonen allein mit 0,86 bis 0,88 | 0,88137 | OK |
| 03 | Bei 8·R_d sind es 0,96 bis 0,97 | 0,96386 | OK |
| 03 | Bei 8·R_d sind es 0,96 bis 0,97 | 0,96753 | OK |
| 03 | Polarringe laufen langsamer als die Scheibe ihrer Galaxie | 0 | OK |
| F:galaxien | fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11} | -0,31796 | OK |
| F:galaxien | fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11} | 0,11468 | OK |
| F:galaxien | fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11} | 0 | OK |
| F:galaxien | Ohne Gang läge es 2,8 σ daneben | 2,7727 | OK |
| F:galaxien | R_d ∝ R₂₀₀ ∝ H^{−2/3}, 3,0 σ | 3,0408 | OK |
| F:lokale_gruppe | ruht nach 13,8 Gyr die Schale bei 0,89 Mpc | 0,88632 | OK |
| F:lokale_gruppe | die 53 Zwerge zwischen 0,6 und 3 Mpc | 53 | OK |
| F:lokale_gruppe | liegen im Median 10 km/s darunter (4 bis 17) | 9,6391 | OK |
| F:lokale_gruppe | liegen im Median 10 km/s darunter (4 bis 17) | 4,4508 | OK |
| F:lokale_gruppe | liegen im Median 10 km/s darunter (4 bis 17) | 17,228 | OK |
| F:lokale_gruppe | ruhte erst die Schale bei 2,7 Mpc | 2,7437 | OK |
| F:lokale_gruppe | ruht nach 13,8 Gyr die Schale bei 0,89 Mpc | 13,828 | OK |
| F:lokale_gruppe | um isolierte Galaxien bei 1,9 Mpc 0,03 dex unter ihr | 0,033402 | OK |
| F:lokale_gruppe | die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc | 0,83264 | OK |
| F:lokale_gruppe | die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc | 1,2345 | OK |
| F:lokale_gruppe | die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc | 1,2345 | OK |
| F:galaxien | an SPARC λ' ∝ M_b^{−0,10 ± 0,02} | -0,099693 | OK |
| F:galaxien | an SPARC λ' ∝ M_b^{−0,10 ± 0,02} | 0,020814 | OK |
| F:galaxien | wie die Schwankung σ(M) des Spektrums über dieselben Bereiche (−0,11) | -0,11297 | OK |
| F:galaxien | der Halbmassenradius geht von 0,59 auf 0,46·r* | 0,58689 | OK |
| F:galaxien | der Halbmassenradius geht von 0,59 auf 0,46·r* | 0,45692 | OK |
| F:galaxien | an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r* | 1,1738 | OK |
| F:galaxien | an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r* | 0,95915 | OK |
| F:galaxien | an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r* | 2,3476 | OK |
| F:galaxien | an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r* | 2,0923 | OK |
| 07 | für die Kuzmin-Scheibe exakt, über der das Feld das eines Punktes im Spiegelbild unter ihr ist: a = ν(√(g_R² + (2πGΣ)²)/a₀)·g_R | 0 | OK |
| 07 | In einer dünnen Scheibe liest die Feldform ν am vollen Betrag der Neigung an der Oberfläche der Schicht | 0 | OK |
| 07 | An SPARC (141 Galaxien, Gas aus V_gas | 141 | OK |
| 07 | liegt die Feldform als Bahnregel Δχ² = 368 über ν·g_N | 368,21 | OK |
| 07 | mit der Steigung 0,05 ± 0,12, wo sie 1 verlangt | 0,045768 | OK |
| 07 | mit der Steigung 0,05 ± 0,12, wo sie 1 verlangt | 0,12279 | OK |
| 03 | innerhalb 2R_d im Median um 0,05 dex | 0,05303 | OK |
| 03 | (flächenschwach 0,08, flächenhell 0,02 | 0,07928 | OK |
| 03 | (flächenschwach 0,08, flächenhell 0,02 | 0,021817 | OK |
| 03 | außerhalb 4R_d und bei Kugeln sind beide gleich | 0 | OK |
| 03 | über einer Kuzmin-Scheibe genau den Wert, den das Licht in der Ebene liest | 0 | OK |
| 03 | Innerhalb 2R_d der SPARC-Scheiben sind das im Median 0,07 dex der Beschleunigung | 0,074106 | OK |
| 03 | Innerhalb 2R_d der SPARC-Scheiben sind das im Median 0,07 dex der Beschleunigung, 8 % des Tempos | 8,1779 | OK |
| 03 | (flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01 | 0,11192 | OK |
| 03 | (flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01 | 0,023957 | OK |
| 03 | (flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01 | 0,0084751 | OK |
| F:galaxien | Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s | 62,524 | OK |
| F:galaxien | Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s | 131,58 | OK |
| F:galaxien | Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s | 29,324 | OK |
| F:galaxien | Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s | 63,362 | OK |
| F:galaxien | In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124 | 3,751 | OK |
| F:galaxien | In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124 | 5,3785 | OK |
| F:galaxien | In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124 | 0,016025 | OK |
| F:galaxien | In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124 | 0,12395 | OK |
| F:galaxien | am Ende trägt die gemischte Mode nur 0,013 bis 0,043 | 0,013408 | OK |
| F:galaxien | am Ende trägt die gemischte Mode nur 0,013 bis 0,043 | 0,04316 | OK |
| F:galaxien | an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135) | 0,072764 | OK |
| F:galaxien | an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135) | 0,046687 | OK |
| F:galaxien | an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135) | 0,13467 | OK |
| F:galaxien | aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157) | 0,097702 | OK |
| F:galaxien | aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157) | 0,057298 | OK |
| F:galaxien | aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157) | 0,15675 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 0,97967 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 3,4919 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 0,95268 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 0,99241 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 0,77952 | OK |
| F:galaxien | An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20 | 1,1995 | OK |
| F:galaxien | bis zur Umkehr gerechnet ist der Anteil der Drehung λ' = J/(√2·M₂₀₀·V₂₀₀·R₂₀₀) = 0,022 bis 0,053 | 0,022276 | OK |
| F:galaxien | bis zur Umkehr gerechnet ist der Anteil der Drehung λ' = J/(√2·M₂₀₀·V₂₀₀·R₂₀₀) = 0,022 bis 0,053 | 0,053413 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 0,013162 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 0,0315 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 2,6857 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 9,3897 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 68,14 | OK |
| F:galaxien | Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe | 83,886 | OK |
| 03 | In einer Scheibe an der Schwelle liegt Gas an der Grenze Q_T = 1 innerhalb 4,2·R_d über der Gezeitendichte | 4,2236 | OK |
| 03 | In SPARC liegen dort vier Zehntel des Gases, aus V_gas im Median 0,37, nach Masse 0,40 | 0,36808 | OK |
| 03 | In SPARC liegen dort vier Zehntel des Gases, aus V_gas im Median 0,37, nach Masse 0,40 | 0,40479 | OK |
| 03 | unabhängig aus M_HI und R_HI 0,36 und 0,40 | 0,35541 | OK |
| 03 | unabhängig aus M_HI und R_HI 0,36 und 0,40 | 0,39825 | OK |
| 03 | Sechs Zehntel verlangen Gas wenig über der Grenze, Q_T = 1,2 bis 1,4 | 1,2034 | OK |
| 03 | Sechs Zehntel verlangen Gas wenig über der Grenze, Q_T = 1,2 bis 1,4 | 1,4049 | OK |
| 03 | Q_T = 1,2 bis 1,4, aus der Flächendichte 1,28 | 1,2773 | OK |
| 03 | aus der Flächendichte, die V_gas als dünne Scheibe verlangt, ebenso 0,37 und 0,40 | 0,37384 | OK |
| 03 | aus der Flächendichte, die V_gas als dünne Scheibe verlangt, ebenso 0,37 und 0,40 | 0,40298 | OK |
| 03 | Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ | 1,5 | OK |
| 03 | Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ | 1,5 | OK |
| 03 | Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ | 1,6832 | OK |
| 03 | Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ | 1,6832 | OK |
| 03 | im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s | 19 | OK |
| 03 | im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s | 28 | OK |
| 03 | im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s | 26,4 | OK |
| 03 | im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s | 30,15 | OK |
| 03 | Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s | 10,694 | OK |
| 03 | Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s | 12,362 | OK |
| 03 | Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s | 19,308 | OK |
| 03 | Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s | 21,619 | OK |
| 03 | aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5 | 0,58745 | OK |
| 03 | aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5 | 0,55463 | OK |
| 03 | aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5 | 0,51527 | OK |
| 03 | aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5 | 0,48284 | OK |
| 03 | aus V_gas und aus M_HI und R_HI 0,5 bis 0,6 | 0,46123 | OK |
| 03 | aus V_gas und aus M_HI und R_HI 0,5 bis 0,6 | 0,64998 | OK |
| 03 | zwischen einem Drittel und dem Dreifachen der gelesenen sind es 0,3 bis 0,8 | 0,27532 | OK |
| 03 | zwischen einem Drittel und dem Dreifachen der gelesenen sind es 0,3 bis 0,8 | 0,84788 | OK |
| 03 | wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen | 0,48422 | OK |
| 03 | wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen | 0,637 | OK |
| 03 | wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen | 0,5862 | OK |
| 03 | wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen | 0,77116 | OK |
| 03 | Mit diesem F₀ gehalten lesen alle Präzisionsdaten χ² = 2025,2 bis 2026,3 | 2025,2 | OK |
| 03 | Mit diesem F₀ gehalten lesen alle Präzisionsdaten χ² = 2025,2 bis 2026,3 | 2026,3 | OK |
| 03 | Δχ² = 1,5 bis 2,6 unter ΛCDM, mit einer freien Zahl weniger | 1,4791 | OK |
| 03 | Δχ² = 1,5 bis 2,6 unter ΛCDM, mit einer freien Zahl weniger | 2,5897 | OK |
| 03 | Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003) | 1,0582 | OK |
| 03 | Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003) | 1,0637 | OK |
| 03 | Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003) | 0,0032183 | OK |
| F:galaxien | An den 36 Satelliten der Milchstraße mit Eigenbewegung | 36 | OK |
| F:galaxien | Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9) | 0,44288 | OK |
| F:galaxien | Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9) | 1,3415 | OK |
| F:galaxien | Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9) | 0,69363 | OK |
| F:galaxien | Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9) | 1,8619 | OK |
| F:galaxien | bei allen übrigen liegt er über 2,6·R_e | 2,622 | OK |
| F:galaxien | im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit | 1 | OK |
| F:galaxien | im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit | 2 | OK |
| F:galaxien | im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit | 6 | OK |
| F:galaxien | im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit | 0,28011 | OK |
| F:galaxien | im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit | 0,28011 | OK |
| F:galaxien | GC5 den Kern in 6,7 bis 21 Gyr | 6,7336 | OK |
| F:galaxien | GC5 den Kern in 6,7 bis 21 Gyr | 20,653 | OK |
| F:galaxien | GC1 in 27 bis 84 Gyr | 27,156 | OK |
| F:galaxien | GC1 in 27 bis 84 Gyr | 83,778 | OK |
| F:galaxien | in 1,2 bis 4,8 Gyr | 1,2256 | OK |
| F:galaxien | in 1,2 bis 4,8 Gyr | 4,7622 | OK |
| F:galaxien | GC5 den Kern in 6,7 bis 21 Gyr | 6,7336 | OK |
| F:galaxien | GC5 den Kern in 6,7 bis 21 Gyr | 20,653 | OK |
| F:galaxien | GC1 in 27 bis 84 Gyr | 27,156 | OK |
| F:galaxien | GC1 in 27 bis 84 Gyr | 83,778 | OK |
| F:galaxien | in 1,2 bis 4,8 Gyr | 1,2256 | OK |
| F:galaxien | in 1,2 bis 4,8 Gyr | 4,7622 | OK |
| 09 | Ein Stern der Form n = 3 (k = 0,075) erreicht höchstens j = 52 (Sonne) | 0,075358 | OK |
| 09 | Ein Stern der Form n = 3 (k = 0,075) erreicht höchstens j = 52 (Sonne) | 51,725 | OK |
| 09 | der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1 | 1,4166 | OK |
| 09 | der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1 | 1,0293 | OK |
| 09 | der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1 | 0,094888 | OK |
| 09 | der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1 | 2,1048 | OK |
| 09 | der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51 | 0,71018 | OK |
| 09 | der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51 | 9,1615 | OK |
| 09 | der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51 | 0,17369 | OK |
| 09 | der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51 | 0,51337 | OK |
| 09 | 1 + α/(2π) = 1,001161 | 1,0012 | OK |
| 09 | Gemessen ist 1,001160 | 1,0012 | OK |
| F:ganzes | +(α/2π)² = +1,35·10⁻⁶ | 1,3489 | OK |
| F:ganzes | gezählt +0,5, −0,25, +0,125, −0,0625, +0,031 | 5 | OK |
| F:ganzes | +0,5, −0,328, +1,181, −1,91, +9,2 | 5 | OK |
| F:ganzes | alle fünf Vorzeichen wie die halbe Drehung | 5 | OK |
| F:ganzes | Rückkehr −1,76·10⁻⁶ | -1,7576 | OK |
| F:ganzes | also −0,33·(α/π)² | -0,32574 | OK |
| 09 | B = M²c²/(eħ) = 4,41·10⁹ T | 4,414 | OK |
| 09 | E = B·c = 1,32·10¹⁸ V/m | 1,3233 | OK |
| 09 | = −2,99·10⁻¹¹ m³/mol | -2,9858 | OK |
| 09 | das sind −2,38·10⁻⁶ cm³/mol | -2,376 | OK |
| 09 | ν(s) = ½ + 1/s + ½·L(s/2) | 1 | OK |
| 09 | m = tanh(μB/kT) = 1 − 2f(2μB/kT) | 1 | OK |
| 09 | m ≈ √(3(1 − T/T_C)) | 0,9996 | OK |
| 09 | h_c ≈ (2/3)(1 − T/T_C)^{3/2} | 1,0002 | OK |
| 09 | eingerastet (ohne Spannung) für |I| ≤ I_c | 1 | OK |
| 09 | Φ₀ = h/(2e) = 2,068·10⁻¹⁵ Wb | 2,0678 | OK |
| 09 | sind das 17 nm | 16,805 | OK |
| 05 | M = m_p·(α/α_G)^{3/2} = 2,30·10²⁷ kg | 2,2971 | OK |
| 05 | Tiefe dort δ = α²·m_e/m_p = 2,9·10⁻⁸ | 2,9002 | OK |
| 05 | (Erde δ = 7,0·10⁻¹⁰) | 6,9613 | OK |
| 05 | (Sonne δ = 2,1·10⁻⁶) | 2,1225 | OK |
| 05 | Verhältnis α^{3/2} = 6,2·10⁻⁴ | 6,2337 | OK |
| 05 | a₀/(πG) = 242 M☉/pc² | 242,07 | OK |
| 05 | Wasserstoff: (8π/3)·α³ = 3,3·10⁻⁶ | 3,2555 | OK |
| 05 | Erde um die Sonne: 2,3·10⁻²⁴ | 2,3382 | OK |
| 05 | in 1,6·10⁻¹¹ s auf den Kern | 1,5562 | OK |
| 05 | trägt die Windung 7/8 der Faltung | 0,875 | OK |
| F:ganzes | ρ_DE = ρ_s·(3 + 16,8·F) | 16,846 | OK |
| F:ganzes | das Maß 17,2 inaktiv | 17,196 | OK |
| F:lesungen | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 0,0071084 | OK |
| F:lesungen | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 0,02398 | OK |
| F:lesungen | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 20,615 | OK |
| F:lesungen | Über alle 15 Bins | 15 | OK |
| F:lesungen | die späten +0,02 ± 0,02 | 0,023309 | OK |
| F:lesungen | die späten +0,02 ± 0,02 | 0,022461 | OK |
| F:lesungen | die roten +0,23 ± 0,01 darüber (χ² 387) | 0,23069 | OK |
| F:lesungen | die roten +0,23 ± 0,01 darüber (χ² 387) | 0,012377 | OK |
| F:lesungen | die roten +0,23 ± 0,01 darüber (χ² 387) | 387,06 | OK |
| F:lesungen | die frühen +0,21 ± 0,01 | 0,20705 | OK |
| F:lesungen | die frühen +0,21 ± 0,01 | 0,013496 | OK |
| F:lesungen | liegen die roten +0,13 ± 0,02 | 0,13063 | OK |
| F:lesungen | liegen die roten +0,13 ± 0,02 | 0,024146 | OK |
| F:lesungen | die frühen +0,14 ± 0,02 über der Regel | 0,1356 | OK |
| F:lesungen | die frühen +0,14 ± 0,02 über der Regel | 0,024879 | OK |
| F:lesungen | wächst es auf +0,29 und +0,26 | 0,29199 | OK |
| F:lesungen | wächst es auf +0,29 und +0,26 | 0,26448 | OK |
| F:lesungen | am Gasanteil (Rangkorrelation −0,04) | -0,036389 | OK |
| F:lesungen | noch am Typ (−0,10) | -0,095828 | OK |
| F:lesungen | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | -0,31212 | OK |
| F:lesungen | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | 168,74 | OK |
| F:lesungen | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | 20,615 | OK |
| F:lesungen | und die späten um −0,26 (157 gegen 14) | -0,2575 | OK |
| F:lesungen | und die späten um −0,26 (157 gegen 14) | 157,15 | OK |
| F:lesungen | und die späten um −0,26 (157 gegen 14) | 14,488 | OK |
| F:lesungen | trifft die roten im Mittel (χ² 146 gegen 387) | 146,32 | OK |
| F:lesungen | trifft die roten im Mittel (χ² 146 gegen 387) | 387,06 | OK |
| F:lesungen | und die frühen (85 gegen 268) | 85,364 | OK |
| F:lesungen | und die frühen (85 gegen 268) | 267,75 | OK |
| F:lesungen | mit χ² 34 bei 20 Werten | 34,369 | OK |
| 03 | bei R₅₀₀ 3,4 bis 7,0 keV | 3,3822 | OK |
| 03 | bei R₅₀₀ 3,4 bis 7,0 keV | 6,9556 | OK |
| F:ganzes | im Median nur 0,4 % des Gases innerhalb R₅₀₀ | 0,37249 | OK |
| 03 | höchstens 3,7 % des Gases abkühlbar | 3,6704 | OK |
| 03 | bei 0,1·R₅₀₀ 4,8 bis 9,6 keV | 4,7804 | OK |
| 03 | bei 0,1·R₅₀₀ 4,8 bis 9,6 keV | 9,6451 | OK |
| 03 | Abkühlzeit von 13,6 bis 38 Gyr | 13,603 | OK |
| 03 | Abkühlzeit von 13,6 bis 38 Gyr | 37,828 | OK |
| F:ganzes | in 13,8 Gyr im Median nur | 13,828 | OK |
| 03 | in 13,8 Gyr höchstens 3,7 % | 13,828 | OK |
| F:ganzes | liegt die Dichte des neutralen Wasserstoffs in sternbildenden Galaxien bei (6,5 ± 1,0 ± 1,2)·10⁷ M☉/Mpc³ | 6,5 | OK |
| F:ganzes | das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14) | 1,1028 | OK |
| F:ganzes | das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14) | 0,26501 | OK |
| F:ganzes | das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14) | 1,1393 | OK |
| F:ganzes | Der Zustandsteil verlangt am Maximum das 1,74-fache des heutigen gebundenen Gases (Sättigung beim Standardwert) oder das 2,2-fache (eigene Form des Modells) | 1,7371 | OK |
| F:ganzes | Der Zustandsteil verlangt am Maximum das 1,74-fache des heutigen gebundenen Gases (Sättigung beim Standardwert) oder das 2,2-fache (eigene Form des Modells) | 2,1855 | OK |
| F:ganzes | das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien | 3,4596 | OK |
| F:ganzes | ihr Verhältnis ist r_e/(3Δs_p) = 2,269·10³⁹ | 2,2687 | OK |
| F:ganzes | ihr Verhältnis ist r_e/(3Δs_p) = 2,269·10³⁹ | 2,2687 | OK |
| F:ganzes | Das Elektron liegt von der Schließung ℓ₀/λ̄_e = 2,196·10³⁹ entfernt, 3 % weniger | 2,1957 | OK |
| F:ganzes | Das Elektron liegt von der Schließung ℓ₀/λ̄_e = 2,196·10³⁹ entfernt, 3 % weniger | 3,2169 | OK |
| F:ganzes | Gemessen ist r*_p/√(λ̄_e·r_e) = 0,984 | 0,98378 | OK |
| F:ganzes | a₀ = G·m_p/(α³·a_B²) = 1,026·10⁻¹⁰ m/s² | 1,0259 | OK |
| F:ganzes | a₀ = G·m_p/(α³·a_B²) = 1,026·10⁻¹⁰ m/s² | 1,0259 | OK |
| F:ganzes | N_b = (3/(2(1 + Q)))·α·(m_Pl²/(m_e·m_p))² = 1,68·10⁸⁰ | 1,6838 | OK |
| F:ganzes | Ist der Formteil 1,026, lesen alle Präzisionsdaten χ² = 2020,5 (mit 1,06: 2019,7) | 2020,5 | OK |
| F:ganzes | Ist der Formteil 1,026, lesen alle Präzisionsdaten χ² = 2020,5 (mit 1,06: 2019,7) | 2019,7 | OK |
| F:ganzes | der Zustandsteil trägt heute 13,3 ± 0,9 % der Dunklen Energie | 13,281 | OK |
| F:ganzes | der Zustandsteil trägt heute 13,3 ± 0,9 % der Dunklen Energie | 0,94814 | OK |
| F:ganzes | Das verlangt F = 2,37 ± 0,17 % der Baryonen als dünnes gebundenes Gas | 2,3651 | OK |
| F:ganzes | Das verlangt F = 2,37 ± 0,17 % der Baryonen als dünnes gebundenes Gas | 0,16885 | OK |
| F:ganzes | verlangt diese Grenze vom Formteil a₀ ≥ 1,054, dort ist F = 1,21 ± 0,13 % | 1,2078 | OK |
| F:ganzes | verlangt diese Grenze vom Formteil a₀ ≥ 1,054, dort ist F = 1,21 ± 0,13 % | 0,12617 | OK |
| F:ganzes | das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien | 0,50226 | OK |
| F:ganzes | das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien | 11419 | OK |
| F:ganzes | der Zustandsteil verlangt in der eigenen Form des Modells dort das 0,34-fache | 0,34395 | OK |
| F:ganzes | bei z = 1,04 ist ihre Flächendichte das 1,48-fache der Schwelle | 1,04 | OK |
| F:ganzes | bei z = 1,04 ist ihre Flächendichte das 1,48-fache der Schwelle | 1,481 | OK |
| F:ganzes | dann 3 % weniger des Gases unter der Gezeitendichte | 3 | OK |
| F:ganzes | dann 3 % weniger des Gases unter der Gezeitendichte | 3 | OK |
| F:ganzes | dann 3 % weniger des Gases unter der Gezeitendichte | 3 | OK |
| F:ganzes | bei z = 0,34 1 % weniger | 1 | OK |
| F:ganzes | bei der dreifachen auf das 1,4-fache | 1,4242 | OK |
| F:ganzes | mit dem gezählten F₀ = 0,67 % beim Standardwert | 0,6741 | OK |
| F:ganzes | nur auf das 1,58-fache, 3,7 Fehler unter der Messung | 1,5831 | OK |
| F:ganzes | nur auf das 1,58-fache, 3,7 Fehler unter der Messung | 3,7362 | OK |
| F:ganzes | Δχ² = 4,6 über dem Gang mit dem Maximum an der Sättigung, 2,0 über ΛCDM | 4,5571 | OK |
| F:ganzes | Δχ² = 4,6 über dem Gang mit dem Maximum an der Sättigung, 2,0 über ΛCDM | 1,9675 | OK |
| F:ganzes | gäbe bei z = 1,04 das 0,67-fache, 5,6 Fehler darunter | 0,67203 | OK |
| F:ganzes | gäbe bei z = 1,04 das 0,67-fache, 5,6 Fehler darunter | 5,5501 | OK |
| F:ganzes | (das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber) | 2,6718 | OK |
| F:ganzes | (das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber) | 1,5685 | OK |
| F:ganzes | (das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber) | 1,6999 | OK |
| F:ganzes | (das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber) | 2,1507 | OK |
| F:ganzes | nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter | 0,38077 | OK |
| F:ganzes | nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter | 0,071 | OK |
| F:ganzes | nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter | 4 | OK |
| F:ganzes | und liegt dann 0,3 über ΛCDM | 0,29969 | OK |
| F:ganzes | = −0,328, ein rationaler Teil 1,368 | -0,32848 | OK |
| F:ganzes | = −0,328, ein rationaler Teil 1,368 | -0,32848 | OK |
| F:ganzes | = −0,328, ein rationaler Teil 1,368 | 1,3681 | OK |
| F:ganzes | ζ(n) = Σ1/kⁿ, mit +1,724 | 1,724 | OK |
| F:ganzes | ζ(n) = Σ1/kⁿ, mit +1,724 | 1,724 | OK |
| F:ganzes | ln 2 = Σ(−1)^{k+1}/k, mit −3,421 | -3,4205 | OK |
| F:ganzes | ln 2 = Σ(−1)^{k+1}/k, mit −3,421 | -3,4205 | OK |
| F:ganzes | Δ₀/kT_c = π·e^{−γ_E} = 1,764 | 1,7639 | OK |
| F:ganzes | Δ₀/kT_c = π·e^{−γ_E} = 1,764 | 1,7639 | OK |
| F:ganzes | das Maximum liegt jenseits z = 1 | 1 | OK |
| F:ganzes | χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein | 1414,1 | OK |
| F:ganzes | χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein | 1416,9 | OK |
| F:ganzes | χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein | 1432,6 | OK |
| F:ganzes | den Preis zahlen die Spektren der Hintergrundstrahlung, Δχ² = 7,2 | 7,1772 | OK |
| F:ganzes | Ist der Formteil frei, liest er a₀ = 1,046 | 1,0461 | OK |
| F:ganzes | das Gas folgt mit 1,9 Fehlern an beiden Punkten | 1,9006 | OK |
| F:ganzes | das Gas folgt mit 1,9 Fehlern an beiden Punkten | 1,8719 | OK |
| F:ganzes | und die Präzisionsdaten liegen 1,2 über ΛCDM | 1,1548 | OK |
| F:ganzes | Dann liegen die Abstände mit dem Gas 2,5 höher als ungekoppelt | 2,5225 | OK |
| F:ganzes | Die Messung liegt 1,3 und 0,9 Fehler darunter | 1,2524 | OK |
| F:ganzes | Die Messung liegt 1,3 und 0,9 Fehler darunter | 0,87168 | OK |
| 07 | Mit D = GM/c ist die Tiefe der Rückkehr ϖ_R = √(D·τ_D) = √(GM/a₀) = r* | 1 | OK |
| F:ganzes | bei 30 kAU (1,5 M☉) ν = 3,8 | 3,8003 | OK |
| F:ganzes | die 1,95-fache Geschwindigkeit | 1,9494 | OK |
| F:ganzes | (Gaia DR3, bis 130 pc, 1 bis 30 kAU) | 130 | OK |
| F:ganzes | der Zuwachs von 20 % nach MOND fehlt | 20 | OK |
| 03 | die Glättungslänge ℓ = 349 kAU | 349,19 | OK |
| 03 | eines Paares von 1 M☉ (277 kAU) | 277,15 | OK |
| F:lesungen | mit χ² 34 bei 20 Werten | 20 | OK |
| F:lesungen | im Mittel je Bin −0,03 bis +0,03 dex | -0,027279 | OK |
| F:lesungen | im Mittel je Bin −0,03 bis +0,03 dex | 0,032415 | OK |
| F:ganzes | verlangt F = 1,8 % | 1,7808 | OK |
| F:ganzes | mit a₀ = 1,06 10,2 ± 1,8 % | 10,202 | OK |
| F:ganzes | mit a₀ = 1,06 10,2 ± 1,8 % | 1,8185 | OK |
| F:ganzes | also F = 1,82 ± 0,32 % | 1,8168 | OK |
| F:ganzes | also F = 1,82 ± 0,32 % | 0,32385 | OK |
| F:ganzes | mit freiem a₀ 7,9 ± 3,2 % | 7,9227 | OK |
| F:ganzes | mit freiem a₀ 7,9 ± 3,2 % | 3,2196 | OK |
| F:ganzes | F = 1,4 ± 0,6 % | 1,4109 | OK |
| F:ganzes | F = 1,4 ± 0,6 % | 0,57337 | OK |
| F:ganzes | ist die ganze Dunkle Energie heute Formteil (a₀ = 1,085) | 1,0854 | OK |
| F:ganzes | ist die Höhe 4,6 ± 1,8 % | 4,5595 | OK |
| F:ganzes | ist die Höhe 4,6 ± 1,8 % | 1,7511 | OK |
| F:ganzes | F = 0,81 ± 0,31 % | 0,81197 | OK |
| F:ganzes | F = 0,81 ± 0,31 % | 0,31184 | OK |
| F:ganzes | Baryonen um z = 0,4 als Gas | 0,39831 | OK |
| 03 | gibt er a₀ = 1,08 | 1,0825 | OK |
| 03 | (H₀ = 68,0) | 68,047 | OK |
| 03 | Δχ² = 5,9 neben ihm | 5,863 | OK |
| 03 | statt 21,4 mit dem gemessenen ω_c | 21,438 | OK |
| 03 | zusammen lesen a₀ = 1,088 | 1,0883 | OK |
| 03 | (1,084 bis 1,092) | 1,0842 | OK |
| 03 | 1,084 bis 1,092) | 1,0925 | OK |
| 03 | bei H₀ = 68,2 | 68,153 | OK |
| 03 | und Ω_m = 0,30 | 0,30196 | OK |
| 03 | Abstände von 79 der 141 Galaxien | 79 | OK |
| 03 | 79 der 141 Galaxien | 141 | OK |
| 03 | Hubble-Fluss mit H₀ = 73 | 73 | OK |
| 03 | H₀ = 73: a₀ = 1,085 | 1,0854 | OK |
| 03 | (0,98 bis 1,19) | 0,98421 | OK |
| 03 | 0,98 bis 1,19), bei | 1,19 | OK |
| 03 | bei Υ_d = 0,66 | 0,663 | OK |
| 03 | w₀ = −0,85 | -0,85297 | OK |
| 03 | w_a = −0,49 | -0,49109 | OK |
| 03 | Δχ² = 8,5 gegen den Formteil allein | 8,5189 | OK |
| 03 | liegt dann bei z = 0,43 | 0,42734 | OK |
| 03 | (0,34 bis 0,55) | 0,34063 | OK |
| 03 | 0,34 bis 0,55), zwischen | 0,54745 | OK |
| 03 | Faltung bei z = 0,32 | 0,32218 | OK |
| 03 | beschleunigten Dehnung bei z = 0,67 | 0,66584 | OK |
| 03 | Der Standardwert 1,06 | 1,06 | OK |
| 09 | g_p = 5,586 ist Inhalt | 5,5857 | OK |
| 09 | Ohne Rückkehr gibt sie 1418,95 MHz | 1419 | OK |
| 09 | mit der ersten 1420,60 MHz | 1420,6 | OK |
| 09 | gemessen sind 1420,41 MHz | 1420,4 | OK |
| 09 | vom g-Faktor zu 0,88·α/2π | 0,88127 | OK |
| 09 | Die restlichen 1,4·10⁻⁴ | 1,3775 | OK |
| F:ganzes | schließt nach 2, 8, 20, 40, 70 (κ = 0,1 | 1 | OK |
| F:ganzes | Mit μ = 0,02 statt 0,03 erscheint 126, dann fehlt 28 | 1 | OK |
| F:ganzes | 356 Spins (42 %), 716 Paritäten (84 %) | 3,5672·10⁵ | OK |
| F:ganzes | 325 Spins und 558 Paritäten (54 %, 92 %) | 3,2556·10⁵ | OK |
| 06 | ν = 1,1567 | 1,1567 | OK |
| 06 | s ≈ 3·10⁵ | 3,0422 | OK |
| 08 | Sonne 4,43 und 8,86 km | 4,4299 | OK |
| 08 | 8,86 km | 8,8598 | OK |
| 02 | Über 93 Elemente | 93 | OK |
| 02 | von 0,81 (As | 0,81376 | OK |
| 02 | bis 8,2 (Ga) | 8,1773 | OK |
| 02 | Ionisierungsenergie −0,48 | -0,48024 | OK |
| 03 | im Median 0,71·r* | 0,70523 | OK |
| 03 | (16 bis 84 %: 0,37 bis 1,28) | 0,36537 | OK |
| 03 | 0,37 bis 1,28) über vier | 1,2765 | OK |
| 03 | M_b^{0,31 ± 0,02} | 0,31436 | OK |
| 03 | 0,31 ± 0,02} | 0,017354 | OK |
| 01 | für Wasserstoff 2,661·10⁻⁵ | 2,6611 | OK |
| 03 | Q = 2π − 1 = 5,283 | 5,2832 | OK |
| 03 | gemessen 5,36 ± 0,06 | 0,064587 | OK |
| 03 | ± 0,06, 1,3 σ | 1,2563 | OK |
| 01 | ⁵⁶Fe trägt 0,94 % weniger | 0,94368 | OK |
| 03 | Materie: Υ_d = 0,67 | 0,668 | OK |
| 03 | liest Pantheon+ dagegen H₀ = 73,5 | 73,544 | OK |
| 03 | H₀ = 73,5 ± 1,0 | 1,0287 | OK |
| 03 | (496 Supernovae, 0,023 < z < 0,15) | 496 | OK |
| 03 | zwischen 73,0 und 74,5 | 73,004 | OK |
| 03 | 73,0 und 74,5. Der | 74,488 | OK |
| 03 | Dipol ist ΔH/H = 0,9 % | 0,88081 | OK |
| 03 | 0,9 % (0,7 σ) | 0,66484 | OK |
| 03 | erst bei z > 0,1 | 1 | OK |
| 03 | Teilkatalog auf 1,8 σ | 1,811 | OK |
| 03 | allein, ohne freie Zahl (Ω_m = 0,315) | 0,31525 | OK |
| 03 | SPARC, 141 Galaxien | 141 | OK |
| 03 | 141 Galaxien, 2688 Punkte | 2688 | OK |
| 03 | In SPARC (175 Galaxien) | 175 | OK |
| 03 | v = 249 km/s bei 5 kpc | 249,25 | OK |
| 03 | 214 bei r* | 214,38 | OK |
| 03 | 196 bei 15 kpc | 196,21 | OK |
| 03 | 178 bei 50 kpc | 178,07 | OK |
| 03 | kleinste gemessene Rotation 17,8 km/s | 17,8 | OK |
| 03 | die kleinste flache 18,6 km/s | 18,6 | OK |
| 03 | Pantheon+, 1590 Supernovae | 1590 | OK |
| 03 | nur Δχ² = 2,0 über | 1,9631 | OK |
| 03 | ΛCDM (Ω_m = 0,304) | 0,30393 | OK |
| 03 | Ω_m = 0,304), 1,4 σ | 1,4011 | OK |
| 03 | Der Formteil trifft es mit χ² = 6,8 | 6,7556 | OK |
| 08 | (10⁶ M☉: s² = 1,4·10¹⁷ | 1,4355 | OK |
| 08 | am Horizont, r* = 36 pc) | 36,262 | OK |
| 08 | = 1 − ⅔√2 = 5,7 % | 5,7191 | OK |
| 08 | sechsmal so viel wie Kernbindung bis zum Eisen | 6,0604 | OK |
| 08 | (j = 1) sind es 1 − 1/√3 = 42 % | 42,265 | OK |
| 08 | M_irr = M_rot = M/√2: Das ist j = 1 | 1 | OK |
| 08 | M_irr·M_rot = Jc/(2G) | 1 | OK |
| 08 | 4π·(r₊² + a_J²) = 4π·r₀(M_irr)² | 1 | OK |
| 08 | am Treffpunkt j = 1 sind es 1 − 1/√2 = 29 % | 29,289 | OK |
| 08 | die Entropie ist dort halb so groß | 0,5 | OK |
| 02 | P/(ρc²) = (dτ/dt|_R)/(dτ/dt|_0) − 1 | 1 | OK |
| 02 | schwach ihr Unterschied ½·δ, die halbe Tiefe | 0,5 | OK |
| 02 | (3/2)·√(1 − r₀/R) − ½. Bei gleichem R | 1 | OK |
| 02 | von 1,4 M☉ und 12 km höchstens mit 0,71 | 0,7144 | OK |
| 02 | Hülle jedes solchen Körpers schneller als ⅓ | 0,33333 | OK |
| 02 | höchstens um z = 2 gerötet an | 2 | OK |
| 02 | außerhalb von (9/8)·r₀ = 27Δs/4 | 6,75 | OK |
| 02 | also 3Δs/R unter 4/9 gegen ½ | 0,44444 | OK |
| 02 | bei 1,4 M☉ über 4,65 km (Buchdahl) | 4,6514 | OK |
| 02 | 3Δs/r* = ½·√(6/(1 + Q)) = 0,49 | 0,4886 | OK |
| 08 | die n-te mit 1/n² | 1 | OK |
| 08 | wächst die Ablenkung so um (π²/6)·(b/r*)² für b ≪ r* | 1 | OK |
| 08 | Bei 4D_E = 5,1 Gpc wächst die Ablenkung dort um 32 % | 32,195 | OK |
| 08 | Bei 4D_E = 5,1 Gpc wächst die Ablenkung dort um 32 % | 32,195 | OK |
| 08 | und der Ring liegt bei 1,20·R_E | 1,198 | OK |
| 08 | Am Trichter vorbei verlängert sie sich um Δt = (r₀/c)·ln(4r₁r₂/b²) | 1 | OK |
| 08 | c·|∂Δt/∂b| = 2r₀/b = 12Δs/b | 1 | OK |
| 08 | R_E²/(2c·D_E) = r₀/c | 1 | OK |
| 08 | bei der Sonne 9,85 µs, am Ort der Linse gezählt | 9,851 | OK |
| 08 | Ihre Summe ist die Lage der Quelle, ihr Produkt −1 | -1 | OK |
| 08 | y = √(b₊/|b₋|) − √(|b₋|/b₊) | 1 | OK |
| 08 | (r₀/c)·(y·√(y² + 4) + 2·ln((√(y² + 4) + y)/(√(y² + 4) − y))) | 1 | OK |
| 05 | für Masse und Ladung des Protons 1,1·10¹⁸ | 1,1116 | OK |
| 05 | für Masse und Ladung des Protons 1,1·10¹⁸ | 1,1116 | OK |
| 05 | ℓ_J/ℓ_D = √(α/α_G), bei jeder Dichte und Temperatur | 1 | OK |
| 05 | N^{2/3}·α_G = α, dieselbe Bedingung | 1 | OK |
| 05 | ρ·ℓ_J³ = M·n·ℓ_D³ mit M = m·(α/α_G)^{3/2} | 1 | OK |
| 04 | also C = −(3/2)·Nk | -1,5 | OK |
| 04 | also C = −(3/2)·Nk | -1,5073 | OK |
| 04 | an der Grenzmasse bei ħk_F = 1,06·m_e·c | 1,0576 | OK |
| 04 | an der Grenzmasse bei ħk_F = 1,06·m_e·c | 1,0576 | OK |
| 04 | weil M·m̃ = m_Pl²/2 fest ist | 0,5 | OK |
| 04 | C = −Mc²/T = −4πk·M/m̃ = −2S | -2 | OK |
| 04 | C = −Mc²/T = −4πk·M/m̃ = −2S | -2 | OK |
| 07 | ⟨Δn²⟩ = n_B·ν = n_B + n_B² = 1/(4·sinh²(s/2)) | 9,992·10⁻⁸ | OK |
| 07 | ⟨Δn²⟩ = −dν/ds = −d(½·coth(s/2))/ds | 0,0050209 | OK |
| 07 | d ln ν/d ln s = −s·⟨Δn²⟩/ν = −χ | 8,2688·10⁻⁵ | OK |
| 07 | q = ½, bei a_N = (ln 2)²·a₀ = 0,48·a₀ | 0,48045 | OK |
| 07 | q = ½, bei a_N = (ln 2)²·a₀ = 0,48·a₀ | 0,48045 | OK |
| 07 | so viel wie Newton, ν = 2 | 2 | OK |
| 07 | n_F·(1 − n_F) = −dn_F/ds ≤ ¼ | 0,25 | OK |
| 07 | n_F·(1 − n_F) = −dn_F/ds ≤ ¼ | 1,3924·10⁻⁵ | OK |
| 03 | = π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³ | 1,7163 | OK |
| 03 | = π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³ | 1,7163 | OK |
| 03 | = π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³ | 1,7163 | OK |
| 03 | 2b₁·k_s = π³/5 = 6,20 | 6,2013 | OK |
| 03 | 2b₁·k_s = π³/5 = 6,20 | 6,2013 | OK |
| 03 | 3H_∞²/(8πG) = b₁·ρ_s | 1 | OK |
| 03 | 2π·max_x (1 − e^{−x}I₀(x))/√x = 3,36 statt π | 3,3583 | OK |
| 03 | 2π·max_x (1 − e^{−x}I₀(x))/√x = 3,36 statt π | 3,3583 | OK |
| 03 | für eine Scheibe aus Gas bei R_d 0,259·v | 0,25937 | OK |
| 03 | für eine Scheibe aus Gas bei R_d 0,259·v | 0,25937 | OK |
| 03 | bei 2R_d 0,144·v (dort | 0,14418 | OK |
| 03 | bei 2R_d 0,144·v (dort | 0,14418 | OK |
| 03 | (dort ν = 2,16, längs | 2,1615 | OK |
| 03 | längs ν(1 − χ/2) = 1,382) | 1,3819 | OK |
| 03 | Streut das Gas mit mehr als 0,224·v_flach | 0,22432 | OK |
| 03 | Streut das Gas mit mehr als 0,224·v_flach | 0,22432 | OK |
| F:galaxien | liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD | 26 | OK |
| F:galaxien | liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD | 26 | OK |
| F:galaxien | liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD | 20 | OK |
| F:galaxien | liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD | 20 | OK |
| F:galaxien | am nächsten UGC 1281 (Sdm) mit 0,22·v_flach, von Hand 0,19·v_flach | 0,21968 | OK |
| F:galaxien | am nächsten UGC 1281 (Sdm) mit 0,22·v_flach, von Hand 0,19·v_flach | 0,18903 | OK |
| 03 | bei 2R_d mindestens 26,7 km/s | 26,687 | OK |
| 03 | bei 2R_d mindestens 26,7 km/s | 26,687 | OK |
| 03 | Er rechnet bei 2R_d 0,104·v | 0,10434 | OK |
| 03 | (19,3 km/s) und liest Q_T | 19,312 | OK |
| 03 | für die Milchstraße ohne Nachbarn 2,2 Mpc | 2,1675 | OK |
| 03 | für die Milchstraße ohne Nachbarn 2,2 Mpc | 2,1686 | OK |
| 03 | Um eine Galaxie, deren Mode übergeben ist (d_eff = 2) | 2,002 | OK |
| 03 | Ω = H_∞, einmal in | 1,001 | OK |
| 03 | Ω = H_∞, einmal in | 1 | OK |
| 03 | 2π/H_∞ = Θ_∞·τ_D = 110 Mrd. Jahren | 110,49 | OK |
| 03 | 2π/H_∞ = Θ_∞·τ_D = 110 Mrd. Jahren | 110,49 | OK |
| 03 | bei M₂₀₀ = 10¹⁵ M☉ 7,0 Mpc, mit | 7,0321 | OK |
| 03 | bei M₂₀₀ = 10¹⁵ M☉ 7,0 Mpc, mit | 7,0321 | OK |
| 03 | mit Ω = √3·H_∞. Zwischen | 1,7321 | OK |
| 02 | Das Mitfalten senkt die Grenzmasse nur auf 1,417 M☉ | 1,4166 | OK |
| 02 | wenn der Halter in der Mitte schon mit 0,999·c läuft | 0,99904 | OK |
| 02 | Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉ | 0,38426 | OK |
| 02 | Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉ | 0,71018 | OK |
| 02 | Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉ | 0,71018 | OK |
| 02 | mit R = 9,16 km und 3Δs/R = 0,114 | 9,1615 | OK |
| 02 | mit R = 9,16 km und 3Δs/R = 0,114 | 9,1615 | OK |
| 02 | mit R = 9,16 km und 3Δs/R = 0,114 | 0,11447 | OK |
| 02 | mit R = 9,16 km und 3Δs/R = 0,114 | 0,11447 | OK |
| 02 | der Halter läuft in der Mitte mit 0,64·c | 0,64042 | OK |
| 02 | der Halter läuft in der Mitte mit 0,64·c | 0,64042 | OK |
| 02 | ω₃·(√(3π)/2)·m_Pl³/m_n² = 5,73 M☉ | 5,7255 | OK |
| 02 | ω₃·(√(3π)/2)·m_Pl³/m_n² = 5,73 M☉ | 5,7252 | OK |
| 02 | das Mitfalten schneidet sie um den Faktor 8 | 8,0621 | OK |
| 02 | ist 0,384 eine reine Formzahl, gleich für jede Masse des Trägers | 0,38426 | OK |
| 02 | ist 0,384 eine reine Formzahl, gleich für jede Masse des Trägers | 0,38426 | OK |
| 02 | M = 9,726·√(1 − β)/(μ²β²)·m_Pl³/m_H² | 9,7255 | OK |
| 02 | M = 9,726·√(1 − β)/(μ²β²)·m_Pl³/m_H² | 9,7256 | OK |
| 02 | Die 45 ist die 3 des Drittels mal der 15 | 45 | OK |
| 02 | Die 45 ist die 3 des Drittels mal der 15 | 45 | OK |
| 02 | Licht und Gas halten gleich, β = ½, bei 146 M☉ | 146,11 | OK |
| 02 | Licht und Gas halten gleich, β = ½, bei 146 M☉ | 146,12 | OK |
| 02 | (Y_p = 0,2449, μ | 0,2449 | OK |
| 02 | μ = 1/(2X + 3Y_p/4) = 0,590 | 0,59036 | OK |
| 02 | Die Zählregel der Elektronen zählt nach Fermi, b = (3π²)^{1/3}/4 | 1 | OK |
| 08 | um einen Haufen mit (1 + Q)·M_b = 10¹⁵ M☉ 6,4·10⁻⁶ unter eins | 6,4304 | OK |
| 08 | um einen Haufen mit (1 + Q)·M_b = 10¹⁵ M☉ 6,4·10⁻⁶ unter eins | 6,4304 | OK |
| 08 | √(1 − 9Δs/R_max) = √(1 − (M/M_N)^{2/3}) | 2,4425·10⁻⁶ | OK |
| 08 | Die Bahnuhr bleibt √(1 − 9Δs/r), die Photonensphäre | 2,2574·10⁻⁴ | OK |
| 08 | Ihre Summe ist null und ihre Quadratsumme 2 | 2 | OK |
| 08 | Umkreisradius 2/√3 | 1 | OK |
| 08 | mit cos 3θ = −M/M_N, die Projektionen | 2,276·10⁻⁶ | OK |
| 08 | von θ = 30° ohne Masse | 30 | OK |
| 08 | bis 60° bei M_N dreht | 59,973 | OK |
| 08 | Ohne Masse sind sie u = 0 und ±1 | 0 | OK |
| 08 | r_F² + r_F·r_H + r_H² = (c/H_∞)² | 2,2204·10⁻⁷ | OK |
| 08 | r_F² + r_F·r_H + r_H² = (c/H_∞)² | 1,2212·10⁻⁶ | OK |
| 08 | M_N = c³/(3√3·G·H_∞) = (√5/6π²)·M_U | 6,9389·10⁻⁶ | OK |
| 08 | = 2,2·10²² M☉,     r_N | 2,1682 | OK |
| 08 | r_N = c/(√3·H_∞) = 3,11 Gpc | 3,1127 | OK |
| 08 | fallen Horizont der Falte, Photonensphäre, R_max und Horizont des Ganzen in r_N zusammen | 0,0081649 | OK |
| 08 | H_∞·τ_D = 2π/Θ_∞ = 2π²/√15 | 2,2204·10⁻⁷ | OK |
| 08 | bleibt die Wärmelänge 2π·r_N | 1 | OK |
| 08 | bleibt die Wärmelänge 2π·r_N | 0,99998 | OK |
| 08 | Gegen t werden dort beide kalt | 5,1641·10⁻⁵ | OK |
| 08 | Der Radius r_N liegt bei 0,58·r* | 0,58296 | OK |
| 08 | um ν = 1,22 verstärkte | 1,2194 | OK |
| 08 | in der Ecke bei M_U zusammen | 2,2204·10⁻⁴ | OK |
| 08 | M_N liegt davor, bei x = 0,98 | 0,97694 | OK |
| 08 | M_N liegt davor, bei x = 0,98 | 0,97694 | OK |
| 08 | liegt 6,6-mal darüber | 6,6207 | OK |
| 08 | liegt 6,6-mal darüber | 6,6207 | OK |
| 03 | ε_b = ρ_b·c²/n_γ = 0,574 eV | 0,57388 | OK |
| 03 | ε_b = ρ_b·c²/n_γ = 0,574 eV | 0,57388 | OK |
| 03 | je Photon (π⁴/15)/(2ζ(3))·kT = 2,70·kT | 2,7012 | OK |
| 03 | je Photon (π⁴/15)/(2ζ(3))·kT = 2,70·kT | 2,7012 | OK |
| 03 | N_eff = 3,046 | 3,046 | OK |
| 03 | je 7/8 bei (4/11)^{1/3} der Temperatur | 0,875 | OK |
| 03 | je 7/8 bei (4/11)^{1/3} der Temperatur | 0,71377 | OK |
| 03 | das 1,69-fache | 1,6918 | OK |
| 03 | das 1,69-fache | 1,6918 | OK |
| 03 | bei kT_eq = 0,80 eV | 0,79924 | OK |
| 03 | bei kT_eq = 0,80 eV | 0,79924 | OK |
| 03 | ℓ = 2,0·10³·λ_H | 2,0187 | OK |
| 03 | ℓ = 2,0·10³·λ_H | 2,0187 | OK |
| 03 | μ = 1/(X + Y_p/4) = 1,225 | 1,225 | OK |
| 03 | μ = 1/(X + Y_p/4) = 1,225 | 1,225 | OK |
| 03 | mit 1 + Q gelesen M_J = 1,0·10³ M☉ | 0,99566 | OK |
| 03 | mit 1 + Q gelesen M_J = 1,0·10³ M☉ | 0,99566 | OK |
| 03 | und die Jeans-Masse ist zeitfest | 1 | OK |
| 03 | an der letzten Streufläche ist M_J = 2,8·10¹⁶ M☉ | 2,7895 | OK |
| 03 | an der letzten Streufläche ist M_J = 2,8·10¹⁶ M☉ | 2,7895 | OK |
| 03 | ℓ_J = 0,43·c/H | 0,42537 | OK |
| 03 | ℓ_J = 0,43·c/H | 0,42537 | OK |
| 03 | = 2,8·10¹³ | 2,8017 | OK |
| 03 | = 2,8·10¹³ | 2,8017 | OK |
| 08 | 398 M☉/pc² | 398,19 | OK |
| 08 | 398 M☉/pc² | 398,19 | OK |
| 08 | Linse bei z = 0,5 | 0,5 | OK |
| 08 | Quelle bei z = 2 | 2 | OK |
| 08 | ist κ = 0,20 | 0,1976 | OK |
| 08 | ist κ = 0,20 | 0,1976 | OK |
| 08 | die Bilder lesen 25 % mehr Faltung | 24,627 | OK |
| 08 | die Bilder lesen 25 % mehr Faltung | 24,627 | OK |
| 08 | −½c² und dem Drehimpuls J/M je Masse | -0,5 | OK |
| 08 | Perizentrum und Apozentrum sind r₋ und r₊ | -0,5 | OK |
| 08 | ihre Exzentrizität e = √(1 − j²) = (r₊ − r₋)/r₀ | 1 | OK |
| 08 | wo das Tempo der Flucht c erreicht | 1 | OK |
| 08 | J/M = 3Δs·c | 1 | OK |
| 08 | mit c durchlaufen | 1 | OK |
| 08 | −½c²/j² | -0,5 | OK |
| 08 | Abbildung der Scheitel auf die Wendepunkte, r₀(M_irr)² = r₀·r₊ | 1 | OK |
| 08 | a_J·r₀ = r₀(M_irr)·r₀(M_rot) | 1 | OK |
| 08 | j der Drehimpuls gegen den des Kreises gleicher Energie | 0,6 | OK |
| 08 | j der Drehimpuls gegen den des Kreises gleicher Energie | 0,6 | OK |
| 08 | e = 0,8 | 0,8 | OK |
| 08 | T = 0,89·T₀ | 0,88889 | OK |
| 08 | S = 0,9·S₀ | 0,9 | OK |
| 08 | −GM/r − ½H_∞²r² | -0,5 | OK |
| 08 | −GM/r − ½H_∞²r² | -0,5 | OK |
| 03 | an der Schwelle 2,14-mal | 2,1385 | OK |
| 03 | an der Schwelle 2,14-mal | 2,1385 | OK |
| 03 | bei a_N = 0,1·a₀ 10,7-mal | 10,714 | OK |
| 03 | bei a_N = 0,1·a₀ 10,7-mal | 10,714 | OK |
| 03 | bei a_N = 0,1·a₀ 10,7-mal | 10,714 | OK |
| 03 | 0 für schnelle Körper und 3/5 für langsame | 0,0016667 | OK |
| 03 | 0 für schnelle Körper und 3/5 für langsame | 0,0016667 | OK |
| 03 | 0 für schnelle Körper und 3/5 für langsame | 0,59999 | OK |
| 03 | 0 für schnelle Körper und 3/5 für langsame | 0,59999 | OK |
| 03 | ν² = Σ(n + 1)·qⁿ | 1 | OK |
| 03 | in einer Kugel ν(1 − χ/2)·ρ_b + χΩ²/(4πG) | 1 | OK |
| 08 | weiter die Faltung der Linse selbst | 1,0292 | OK |
| 03 | ϑ = 3,49 | 3,4875 | OK |
| 03 | linear δ_v = −2,72 | -2,7172 | OK |
| 03 | linear δ_v = −2,72 | -2,7164 | OK |
| 03 | 1 + δ = 0,20 | 0,20472 | OK |
| 03 | 1 + δ = 0,20 | 0,20476 | OK |
| 03 | 1,70-fachen | 1,6967 | OK |
| 03 | 1,70-fachen | 1,6966 | OK |
| 03 | H_v = (4/3)·H | 4 | OK |
| 03 | H_v = (4/3)·H | 3,9999 | OK |
| 03 | bei denselben δ_v und 1 + δ | -2,7168 | OK |
| 03 | bei denselben δ_v und 1 + δ | -2,7167 | OK |
| 03 | bei denselben δ_v und 1 + δ | 0,20486 | OK |
| 03 | bei denselben δ_v und 1 + δ | 0,20487 | OK |
| 03 | auf 0,1 % | 0,10243 | OK |
| 03 | auf 0,1 % | 0,10269 | OK |
| 03 | umkehrt (Θ = 1,020) | 1,0201 | OK |
| 03 | 11,8-fache | 11,763 | OK |
| 03 | 11,8-fache | 11,763 | OK |
| 03 | (3π/4)² = 5,55 | 5,5517 | OK |
| 03 | (3π/4)² = 5,55 | 5,5517 | OK |
| 03 | 0,72·R_max | 0,71735 | OK |
| 03 | 0,72·R_max | 0,71735 | OK |
| 03 | bei 2R_d das 2,08-fache | 2,0813 | OK |
| 03 | 0,888·σ²/√(2πGΣa₀) | 0,8879 | OK |
| 03 | 0,549·σ²/(πGΣ) | 0,54931 | OK |
| 03 | = ρ_s·c² = 1,74·10⁻¹⁰ Pa | 1,7399 | OK |
| 03 | höchstens π³a₀²/(30G) | 1,7399 | OK |
| 03 | bei R_d 2,16 | 2,1649 | OK |
| 03 | bei R_d 2,16 | 2,1649 | OK |
| 03 | bei 2R_d 3,45 | 3,4487 | OK |
| 03 | bei 2R_d 3,45 | 3,4487 | OK |
| 03 | (dort d_eff = 1,58) | 1,5801 | OK |
| 03 | gäbe dort 2,98 | 2,9823 | OK |
| 03 | bei 2R_d 7,45 | 7,4545 | OK |
| 03 | bei 2R_d 7,45 | 7,4545 | OK |
| 03 | x₀ = 3,830 | 3,83 | OK |
| 03 | bei 217,5 GHz | 217,51 | OK |
| 03 | um den Anteil 4y | 4 | OK |
| 03 | die Stufe d = 2, bleibt | 3,6944·10⁻¹⁶ | OK |
| 03 | ΔT/T = −2y | 2 | OK |
| 03 | ∫s^d·⟨Δn²⟩ds = d·Γ(d)ζ(d) | 1 | OK |
| 03 | s·coth(s/2) = s·(2ν − 1) = d | 1 | OK |
| 03 | = 1,09·√(GM_b·a₀) | 1,0879 | OK |
| 03 | heute 2,3·10³ Gyr | 2,2554 | OK |
| 03 | noch 1,4·10² Gyr | 1,4096 | OK |
| 03 | heute kreuzt (Θ = 1,020) | 1,0201 | OK |
| 08 | mit der Determinante 1 und der Spur 2·cos(2π·Re τ_M) | 1 | OK |
| 08 | mit der Determinante 1 und der Spur 2·cos(2π·Re τ_M) | 2 | OK |
| 08 | λ² − Spur·λ + 1 = 0 | 1 | OK |
| 08 | liegen die Eigenwerte auf dem Einheitskreis | 1 | OK |
| 08 | Bei 18Δs fallen sie in eins zusammen | 18 | OK |
| 08 | (δu, δu′) → (δu + 2π·δu′, δu′) | 2 | OK |
| 08 | e^{±2π·√(6δ − 1)} | 1 | OK |
| 08 | an der Photonensphäre mit der Spur 2·cosh 2π = 535,5 | 535,49 | OK |
| 08 | Bei 24Δs = 8d ist Re τ_M = ½ | 24 | OK |
| 08 | Bei 24Δs = 8d ist Re τ_M = ½ | 8 | OK |
| 08 | Bei 24Δs = 8d ist Re τ_M = ½ | 0,5 | OK |
| 08 | kehrt nach einem Umlauf mit umgekehrtem Vorzeichen zurück | -1 | OK |
| 08 | schließt nach zwei Umläufen | 2 | OK |
| 08 | Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität | 1 | OK |
| 08 | Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität | 2 | OK |
| 08 | Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität | 0 | OK |
| 08 | Die Regel allein senkt die Spur von 2 höchstens bis 2·cos(2π√2) = −1,72 | 2 | OK |
| 08 | Die Regel allein senkt die Spur von 2 höchstens bis 2·cos(2π√2) = −1,72 | -1,7164 | OK |
| 08 | bringt um eine Punktmasse keine Kreisbahn aus der Drehung | 1 | OK |
| 08 | An der letzten stabilen Kreisbahn um einen gebundenen Bereich liegt die zweite Scherung | 0 | OK |
| 08 | An der letzten stabilen Kreisbahn um einen gebundenen Bereich liegt die zweite Scherung | 0 | OK |
| 08 | und bis R_max dehnt die Abbildung | 1 | OK |
| 08 | mit der Identität bei R_max/√3 | 3 | OK |
| 08 | der halben Drehung bei √(7/15)·R_max | 7 | OK |
| 08 | liegt die halbe Drehung bei R_max/∛5 | 5 | OK |
| 01 | hat ihre Abbildung die Determinante 1 | 1 | OK |
| 01 | und die Spur 2·cos(√K·λ) | 1 | OK |
| 01 | für K < 0 ist sie 2·cosh(√−K·λ) | 1 | OK |
| 01 | für K = 0 schert sie (Spur 2) | 2 | OK |
| 01 | q ist am größten bei |k| = κ/σ, dort gleich 1/Q_T | 1 | OK |
| 01 | q ist am größten bei |k| = κ/σ, dort gleich 1/Q_T | 0,99998 | OK |
| 01 | die Jeans-Länge, Q_T = 1 | 1 | OK |
| 01 | Ihre Antwort Σqⁿ = 1/(1 − q) ist ω_frei²/ω² | 1 | OK |
| 01 | Ihre Antwort Σqⁿ = 1/(1 − q) ist ω_frei²/ω² | 1 | OK |
| 01 | (dz/dλ)² − (E/2)·z² = GM/2 | 0,5 | OK |
| 01 | K = −E/2 | 1 | OK |
| 01 | der rollende Kreis mit θ = 2√K·λ | 1 | OK |
| 01 | steht z mit umgekehrtem Vorzeichen | -1 | OK |
| 01 | r³ = (9/2)·GM·t² | 4,5 | OK |
| 01 | r³ = (9/2)·GM·t² | 4,5 | OK |
| 03 | 2b₁·ρ_s ist ihre Gezeitendichte 3𝒯/(4πG) | 1 | OK |
| 03 | R_max³ = GM_dyn/H_∞² | 1 | OK |
| 03 | Θ = Θ_∞/√3 | 1 | OK |
| 03 | Θ = 1,020 | 1,0201 | OK |
| 03 | 1,63·R_max | 1,6315 | OK |
| 03 | 4,34-fache | 4,3423 | OK |
| 02 | A² = 3c²/(8πGρ) = R³/r₀ | 1 | OK |
| 02 | für die Erde 1,14 AE | 1,1414 | OK |
| 02 | sin²ψ_R = 2δ | 1 | OK |
| 02 | an der Hülle ist das der Takt | 1 | OK |
| 02 | (3/2)·cos ψ_R − ½·cos ψ | 1 | OK |
| 02 | In der Mitte steht sie bei cos ψ_R = ⅓ | 0,33333 | OK |
| 02 | 70,5° vom Pol | 70,529 | OK |
| 02 | bis 36° vom Pol | 35,943 | OK |
| 02 | R³/(6A²) = Δs | 1 | OK |
| 02 | Bei 90° läge die Hülle am Äquator, R = A = r₀ | 90 | OK |
| 02 | c²/A² = (8πG/3)·ρ + H_∞² | 1 | OK |
| 02 | die halbe Kugel vom Radius c/H_∞ | 1 | OK |
| 02 | ihr Äquator der Horizont des Ganzen | 90 | OK |
| 02 | im Augenblick seiner Umkehr | 1 | OK |
| 08 | = (Θ_∞/2π)²·ℓ₀ = 1,06 Gpc | 1,0578 | OK |
| 08 | = 8,5·10²⁰ M☉, x = 0,95 | 8,51 | OK |
| 08 | = 8,5·10²⁰ M☉, x = 0,95 | 0,95416 | OK |
| 08 | bei v_flach = a₀/H_∞ = 0,196·c | 0,19621 | OK |
| 08 | (6·10¹⁰ M☉) nach der Regel 345 | 345,1 | OK |
| 08 | ·τ_D = 0,45 Mrd. Jahre | 0,45049 | OK |
| 08 | ungedreht 3/(2√2) − 1 = 0,061 | 0,06066 | OK |
| 08 | 0,061 mit t_e = 27 Mio. Jahre | 27,327 | OK |
| 08 | (j = 1) √3 − 1 = 0,73 | 0,73205 | OK |
| 08 | 0,73 mit t_e = 330 Mio. Jahre | 329,78 | OK |
| 08 | Θ = 0,1, nach 0,95 Mrd. Jahren | 0,953 | OK |
| 08 | um den Faktor e³⁵ | 34,874 | OK |
| 08 | bei j = 1 um den Faktor 18 | 17,99 | OK |
| 05 | Bis zum Fall vergeht t = E/((l + 2)·P) | 1 | OK |
| 05 | Bis zum Fall vergeht t = E/((l + 2)·P) | 1 | OK |
| 05 | Die Erde fiele so in 1,1·10²³ Jahren auf die Sonne | 1,0692 | OK |
| 05 | Die Erde fiele so in 1,1·10²³ Jahren auf die Sonne | 1,0692 | OK |
| 05 | t = (5/512)·(3Δs/c)·δ⁻⁴ | 1 | OK |
| 05 | N = δ^{−5/2}/(64√2·π) | 1 | OK |
| 05 | von δ = 0,001 an 1,1·10⁵ Umläufe | 1,1121 | OK |
| 05 | von δ = 0,001 an 1,1·10⁵ Umläufe | 1,1121 | OK |
| 05 | bei je einer Sonnenmasse in 13 Stunden | 13,361 | OK |
| 05 | bei je einer Sonnenmasse in 13 Stunden | 13,361 | OK |
| 05 | n = 3 über die Windung | 3 | OK |
| 05 | n = 5 über eine Verformung | 5 | OK |
| 05 | seit höchstens t = Ω/((n − 1)·|Ω̇|) = E/(l·P) | 1 | OK |
| 05 | seit höchstens t = Ω/((n − 1)·|Ω̇|) = E/(l·P) | 1 | OK |
| 03 | r* = 0,36 kpc | 0,36262 | OK |
| 03 | bis 6,8 kpc | 6,7847 | OK |
| 03 | mit 36,0 km/s | 36,012 | OK |
| 03 | ν = 6,03 | 6,0305 | OK |
| 03 | Newton, 14,7 km/s | 14,664 | OK |
| 03 | = 120,5 km/s | 120,53 | OK |
| 03 | als 20,7 km/s | 20,739 | OK |
| 03 | 7,2·10⁶ M☉ | 7,1864 | OK |
| 03 | 5,0·10⁸ M☉ | 5,0305 | OK |
| F:lokale_gruppe | laufen die Schwerpunkte mit 42,1 km/s | 42,114 | OK |
| F:lokale_gruppe | laufen die Schwerpunkte mit 42,1 km/s | 42,114 | OK |
| F:lokale_gruppe | je Richtung auf ±48 und ±44 km/s | 47,837 | OK |
| F:lokale_gruppe | je Richtung auf ±48 und ±44 km/s | 44,157 | OK |
| F:lokale_gruppe | vor 5,3 Gyr bei 54 kpc | 5,2558 | OK |
| F:lokale_gruppe | vor 5,3 Gyr bei 54 kpc | 5,2558 | OK |
| F:lokale_gruppe | vor 5,3 Gyr bei 54 kpc | 53,613 | OK |
| F:lokale_gruppe | vor 5,3 Gyr bei 54 kpc | 53,613 | OK |
| F:lokale_gruppe | nur 0,2 % der Baryonen der Milchstraße | 0,19929 | OK |
| F:lokale_gruppe | 2,5 % der von M31 | 2,4652 | OK |
| F:lokale_gruppe | In 3,0 Gyr kommen sich beide wieder so nahe | 2,9635 | OK |
| F:lokale_gruppe | In 3,0 Gyr kommen sich beide wieder so nahe | 2,9635 | OK |
| F:lokale_gruppe | stehen ihre 12 hellsten Satelliten | 12 | OK |
| F:lokale_gruppe | bei 0,44 ± 0,08 | 0,44127 | OK |
| F:lokale_gruppe | bei 0,44 ± 0,08 | 0,44127 | OK |
| F:lokale_gruppe | bei 0,44 ± 0,08 | 0,083333 | OK |
| F:lokale_gruppe | frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰) | 7,936 | OK |
| F:lokale_gruppe | frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰) | 5,2254 | OK |
| F:lokale_gruppe | frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰) | 9,9246 | OK |
| F:lokale_gruppe | frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰) | 0,50001 | OK |
| F:lokale_gruppe | frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰) | 0,50008 | OK |
| F:lokale_gruppe | gegen 6,6·10¹⁰ M☉ ihres Bereichs | 6,6051 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 92,941 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 16,744 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 10,225 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,38702 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 1,1577 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,58847 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 92,941 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 16,744 | OK |
| F:lokale_gruppe | stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82 | 10,224 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,38707 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 1,1577 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,58847 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,57735 | OK |
| F:lokale_gruppe | liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58 | 0,57735 | OK |
| F:lokale_gruppe | läuft die LMC zwischen 45 und 182 kpc | 44,717 | OK |
| F:lokale_gruppe | läuft die LMC zwischen 45 und 182 kpc | 44,717 | OK |
| F:lokale_gruppe | läuft die LMC zwischen 45 und 182 kpc | 181,53 | OK |
| F:lokale_gruppe | läuft die LMC zwischen 45 und 182 kpc | 181,53 | OK |
| F:lokale_gruppe | liegt er bei 19 kpc und nur Pictor II darin | 18,525 | OK |
| F:lokale_gruppe | liegt er bei 19 kpc und nur Pictor II darin | 18,525 | OK |
| F:lokale_gruppe | Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc | 1,472 | OK |
| F:lokale_gruppe | Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc | 1,4716 | OK |
| F:lokale_gruppe | Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc | 68,51 | OK |
| F:lokale_gruppe | Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc | 68,51 | OK |
| F:lokale_gruppe | und 7 der 8 Zwerge | 7 | OK |
| F:lokale_gruppe | und 7 der 8 Zwerge | 7 | OK |
| F:lokale_gruppe | und 7 der 8 Zwerge | 8 | OK |
| F:lokale_gruppe | liegt er bei 19 kpc und nur Pictor II darin | 1 | OK |
| F:lokale_gruppe | liegt er bei 19 kpc und nur Pictor II darin | 1 | OK |
| F:lesungen | bei 14 kpc das 2,2-fache von Newton | 2,2394 | OK |
| F:lesungen | bei 14 kpc das 2,2-fache von Newton | 2,2394 | OK |
| 03 | voll ionisiert in der Urmischung 1,93 | 1,9302 | OK |
| 03 | kT_ab = 1,4 meV | 1,4052 | OK |
| 03 | kT_ab = 1,4 meV | 1,4018 | OK |
| 03 | (16 K) | 16,306 | OK |
| 03 | z_ab = 5,0 | 4,9828 | OK |
| 03 | z_ab = 5,0 | 4,9687 | OK |
| 03 | noch 0,68 der Temperatur des Lichts | 0,67781 | OK |
| 03 | 2,8·10⁻¹⁰ der Energie des Lichts | 2,7711 | OK |
| 03 | 0,025/(1 + z) | 0,025011 | OK |
| 08 | wo γ² = 2 | 2 | OK |
| 08 | 𝒯 = 5GM/r³ | 5 | OK |
| 08 | M_V = √5·(¼/δ_m)^{3/2}·m | 1 | OK |
| 08 | M_V = 9,0·10⁷ M☉ | 90,39 | OK |
| 08 | bei 0,6 AE | 0,59579 | OK |
| 08 | gleich 4GM/r³ | 4 | OK |
| 08 | M_V = 2·(⅙/δ_m)^{3/2}·m | 1 | OK |
| 08 | für die Sonne 4,4·10⁷ M☉ | 44,008 | OK |
| 03 | tief ¾ | 0,75 | OK |
| 03 | zwei Drittel nah | 0,66667 | OK |
| 03 | tief 2·ln(R/b_min)/N | 2 | OK |
| 03 | ein Viertel des Werts nach Newton | 0,25 | OK |
| 03 | (ν = 3,3) 3,5-mal | 3,3035 | OK |
| 03 | (ν = 3,3) 3,5-mal | 3,5064 | OK |
| 03 | (ν = 10,5) 29-mal | 10,508 | OK |
| 03 | (ν = 10,5) 29-mal | 28,582 | OK |
| 03 | dort 1,1 und 2,9 | 1,0881 | OK |
| 03 | dort 1,1 und 2,9 | 2,9008 | OK |
| 03 | (X² = 6) | 6 | OK |
| 03 | sind das 0,0074 | 0,0073832 | OK |
| 03 | (γ = 9/2) | 4,5 | OK |
| 03 | je Relaxation 1,0·10⁻³ | 0,99904 | OK |
| 03 | von der Schwelle aus 1,2·10⁻¹⁵ | 1,2203 | OK |
| 03 | √5/(3π) = 0,080 | 79,676 | OK |
| 03 | √5/(3π) = 0,080 | 79,676 | OK |
| 03 | zwei bei 0,420 | 420,32 | OK |
| 03 | einer bei 0,420 | 420,32 | OK |
| 03 | keiner bei 0,080 | 79,676 | OK |
| 03 | mindestens zwei positive Eigenwerte | 1 | OK |
| 03 | 3σ/√(10π) = 0,535·σ | 535,24 | OK |
| 03 | 3σ/√(10π) = 0,535·σ | 535,24 | OK |
| 03 | trägt 1/3 von Σλ² | 1 | OK |
| 03 | die Scherung 2/3 | 1 | OK |
| 03 | (σ²/15) | 15 | OK |
| 03 | je ein Achtel | 1 | OK |
| 03 | fällt sie nur mit λ > 0,71 | 0,70988 | OK |
| 03 | d_eff = 3 − χ, tief 2 | 2 | OK |
| 03 | in der Tiefe 1, 2, 2 | 2 | OK |
| 03 | v² = 2(1 + Q)·Gμ_b | 2 | OK |
| 01 | (29 + 6√6)/(29 − 6√6) = 3,055 | 3,0551 | OK |
| 01 | (29 + 6√6)/(29 − 6√6) = 3,055 | 3,0551 | OK |
| 01 | im Mittel 6,11 Fäden | 6,1102 | OK |
| 01 | im Mittel 6,11 Fäden | 6,1102 | OK |
| 01 | 2b₁ = 6 | 6 | OK |
| 03 | nicht 1, 1,5, 2 | 1,5 | OK |
| 10 | Je hundertfach kleineren Verlust wächst die Verweilzeit im Kanal zehnfach | 9,9551 | OK |
| 03 | im Feld der Sterne c/v_flach ≤ 0,58·c/σ_Gas | 0,57735 | OK |
| 03 | im Feld der Sterne c/v_flach ≤ 0,58·c/σ_Gas | 0,57735 | OK |
| 02 | Carter, 5,9·10⁻³⁹ | 5,9061 | OK |
| 02 | gegen 2,0·10⁻³⁹ | 20,061 | OK |
| 01 | auf jede Stelle gleich | 0 | OK |
| 01 | kT_ab ∝ α^{−4/5}·m_e^{6/5} | -0,8 | OK |
| 01 | kT_ab ∝ α^{−4/5}·m_e^{6/5} | 1,2 | OK |
| 01 | bis auf den Faktor 3 | 2,9441 | OK |
| 01 | 3 % auseinander | 3,2169 | OK |
| 03 | (D/H und Y_p, ω_b = 0,0219 ± 0,0002) | 0,021929 | OK |
| 03 | (D/H und Y_p, ω_b = 0,0219 ± 0,0002) | 2,1727·10⁻⁴ | OK |
| 03 | ist (1 + Q)·ω_b = 0,1384 ± 0,0009 | 0,13842 | OK |
| 03 | ist (1 + Q)·ω_b = 0,1384 ± 0,0009 | 8,5701·10⁻⁴ | OK |
| 03 | und Q = 5,32 ± 0,06 | 5,3186 | OK |
| 03 | und Q = 5,32 ± 0,06 | 0,060881 | OK |
| 03 | 2π − 1 liegt dann Δχ² = 0,3 daneben, das gemessene Verhältnis Δχ² = 0,6 | 0,33938 | OK |
| 03 | 2π − 1 liegt dann Δχ² = 0,3 daneben, das gemessene Verhältnis Δχ² = 0,6 | 0,55753 | OK |
| 03 | ist Q = 5,40 ± 0,08 | 539,67 | OK |
| 03 | ist Q = 5,40 ± 0,08 | 0,078183 | OK |
| 03 | Der Zustandsteil hebt Q um 0,08 | 0,078141 | OK |
| 03 | Der Zustandsteil hebt Q um 0,08 | 0,081142 | OK |
| 03 | Gegen das Atmen aller Präzisionsdaten (H₀ = 68,25) | 68,249 | OK |
| 03 | linear H_v = (1 − f·δ_L/3)·H mit f = 0,52 heute | 0,51615 | OK |
| 03 | linear mindestens δ_L = −0,41 tief | -0,40824 | OK |
| 03 | 0,59 Gpc im Radius messen | 0,59038 | OK |
| 03 | Diese Tiefe ist das 43-fache der Schwankung | 42,879 | OK |
| 03 | für jeden Rand bis z = 0,3 mindestens Δχ² = 41 daneben | 41,074 | OK |
| 03 | (r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2 | 147,41 | OK |
| 03 | (r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2 | 68,293 | OK |
| 03 | (r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2 | 0,23352 | OK |
| 03 | die Cepheiden 73,2 ± 1,1 | 73,201 | OK |
| 03 | die Cepheiden 73,2 ± 1,1 | 1,0642 | OK |
| 03 | um 5·lg(73,2/68,3) = 0,15 mag dunkler | 0,15069 | OK |
| F:lesungen | 68,25 ± 0,10, sind es 4,7 σ | 0,097269 | OK |
| F:lesungen | 68,25 ± 0,10, sind es 4,7 σ | 4,7376 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 70,39 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 1,9358 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 24 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 68,81 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 2,2241 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 67,8 | OK |
| F:lesungen | liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB) | 2,72 | OK |
| F:lesungen | 1,1, 0,3 und 0,2 Fehler neben dem Atmen | 1,1045 | OK |
| F:lesungen | 1,1, 0,3 und 0,2 Fehler neben dem Atmen | 0,25193 | OK |
| F:lesungen | 1,1, 0,3 und 0,2 Fehler neben dem Atmen | 0,16502 | OK |
| F:lesungen | das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4 | 2,6957 | OK |
| F:lesungen | das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4 | 11,1 | OK |
| F:lesungen | das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4 | 7,5 | OK |
| F:lesungen | das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4 | 11,3 | OK |
| F:lesungen | das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4 | 11,4 | OK |
| 03 | ν_G(−ln q) − 1 = 16,8 bei Q = 2π − 1 | 16,846 | OK |
| 03 | heute 1,21 % der Baryonen als neutrales Gas in Galaxien (ALFALFA α.100) | 1,2106 | OK |
| 03 | unabhängig davon 1,17 % (HIPASS) | 1,1718 | OK |
| 03 | bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019) | 0,34 | OK |
| 03 | bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019) | 1,1207 | OK |
| 03 | bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019) | 0,26932 | OK |
| 03 | verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 % | 1,06 | OK |
| 03 | verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 % | 0,72897 | OK |
| 03 | verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 % | 0,11229 | OK |
| 03 | sechs Zehntel des gezählten Gases | 0,61195 | OK |
| 03 | Δχ² = 2,6 unter ΛCDM | 2,646 | OK |
| 03 | den Formteil a₀ = 1,043 ± 0,003 | 1,0431 | OK |
| 03 | den Formteil a₀ = 1,043 ± 0,003 | 0,003277 | OK |
| 03 | (DESI DR2 mit r_d aus ω_b, Pantheon+, D/H) 1,05 ± 0,01 | 1,0523 | OK |
| 03 | (DESI DR2 mit r_d aus ω_b, Pantheon+, D/H) 1,05 ± 0,01 | 0,0081115 | OK |
| 03 | alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM | 1,0465 | OK |
| 03 | alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM | 0,0033879 | OK |
| 03 | alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM | 4,1097 | OK |
| 03 | Der Standardwert mit allem gezählten Gas liegt Δχ² = 17 über ΛCDM | 16,752 | OK |
| 03 | je Baryon für v_flach = 200 bis 20 km/s das 27- bis 266-fache | 26,588 | OK |
| 03 | je Baryon für v_flach = 200 bis 20 km/s das 27- bis 266-fache | 265,88 | OK |
| 03 | nur 0,003 bis 0,03 % der Baryonen als dünnes Gas | 0,0027417 | OK |
| 03 | nur 0,003 bis 0,03 % der Baryonen als dünnes Gas | 0,027417 | OK |
| F:lesungen | bis z = 0,34 einen Anstieg auf das 2,2-fache | 0,34 | OK |
| F:lesungen | bis z = 0,34 einen Anstieg auf das 2,2-fache | 2,1632 | OK |
| F:lesungen | mit dem gezählten Gang liegen sie Δχ² = 1,6 höher | 1,6364 | OK |
| F:lesungen | den linearen Teil bei 0,83 ± 0,07 | 0,82687 | OK |
| F:lesungen | den linearen Teil bei 0,83 ± 0,07 | 0,070524 | OK |
| F:lesungen | den Rest bei 0,68 ± 0,06 (χ² = 269) | 0,68101 | OK |
| F:lesungen | den Rest bei 0,68 ± 0,06 (χ² = 269) | 0,05598 | OK |
| F:lesungen | den Rest bei 0,68 ± 0,06 (χ² = 269) | 269 | OK |
| F:lesungen | bei k = 1/Mpc 0,99, bei 5/Mpc 0,87 | 0,98836 | OK |
| F:lesungen | bei k = 1/Mpc 0,99, bei 5/Mpc 0,87 | 0,86907 | OK |
| F:lesungen | 10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio | 10,86 | OK |
| F:lokale_gruppe | (Local Volume Database, Σm_i/M = 0,08) | 0,081846 | OK |
| F:lokale_gruppe | sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°) | 22,477 | OK |
| F:lokale_gruppe | sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°) | 22,477 | OK |
| F:lokale_gruppe | sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°) | 82,033 | OK |
| F:lokale_gruppe | sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°) | -43,483 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | 32 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | 4 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | 56 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | 9 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | -34 | OK |
| F:lokale_gruppe | mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°) | 10 | OK |
| F:lokale_gruppe | trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben | 70,24 | OK |
| F:lokale_gruppe | trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben | 2,3808 | OK |
| F:lokale_gruppe | trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben | 22,286 | OK |
| F:lesungen | 10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio | 0,03 | OK |
| F:lesungen | 10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio | 3,2 | OK |
| F:lesungen | an 12 Millionen Riesensternen | 12 | OK |
| F:lesungen | 6,5 bis 6,7 Fehler neben der Messung | 6,4927 | OK |
| F:lesungen | 6,5 bis 6,7 Fehler neben der Messung | 6,6912 | OK |
| F:lesungen | rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung | 9,4081 | OK |
| F:lesungen | rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung | 9,4413 | OK |
| F:lesungen | rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung | 0,44334 | OK |
| F:lesungen | rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung | 0,45373 | OK |
| F:lesungen | kippen die Satelliten den Rand um 0,06° | 0,056633 | OK |
| F:lesungen | kippen die Satelliten den Rand um 0,06° | 0,056772 | OK |
| F:lesungen | dreht die Knotenlinie rückläufig mit 13,9 km/s/kpc | 13,861 | OK |
| F:lesungen | mit dem 25-fachen Antrieb der LMC rückläufig mit 21,8 km/s/kpc, und der Rand kippt dann um 0,49° | 21,762 | OK |
| F:lesungen | mit dem 25-fachen Antrieb der LMC rückläufig mit 21,8 km/s/kpc, und der Rand kippt dann um 0,49° | 0,48818 | OK |
| F:lesungen | nur um 0,04° gegen die innere Scheibe | 0,039309 | OK |
| F:lesungen | nur um 0,04° gegen die innere Scheibe | 0,039792 | OK |
| F:lesungen | rückläufig, mit 10 bis 11 km/s/kpc | 9,9166 | OK |
| F:lesungen | rückläufig, mit 10 bis 11 km/s/kpc | 10,552 | OK |
| F:lesungen | (v/c)² = 4·10⁻⁷ des Umlaufs | 4,4123 | OK |
| F:lesungen | 27 liegen dort −0,03 ± 0,02 dex auf der Regel | 27 | OK |
| F:lesungen | 27 liegen dort −0,03 ± 0,02 dex auf der Regel | -0,030704 | OK |
| F:lesungen | 27 liegen dort −0,03 ± 0,02 dex auf der Regel | 0,016596 | OK |
| F:lesungen | die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27 | 7 | OK |
| F:lesungen | die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27 | -0,0065699 | OK |
| F:lesungen | die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27 | 0,02881 | OK |
| F:lesungen | die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27 | -0,26654 | OK |
| F:lesungen | die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°) | 7,4008 | OK |
| F:lesungen | die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°) | 32 | OK |
| F:lesungen | die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°) | 3 | OK |
| F:lesungen | über alle 15 Bins um einen festen Faktor höher | 15 | OK |
| F:lesungen | +0,22 ± 0,03 dex (χ² = 17) | 0,22175 | OK |
| F:lesungen | +0,22 ± 0,03 dex (χ² = 17) | 0,027604 | OK |
| F:lesungen | +0,22 ± 0,03 dex (χ² = 17) | 16,892 | OK |
| F:lesungen | frei angepasstem f = 3,1 ± 0,2 | 3,1478 | OK |
| F:lesungen | frei angepasstem f = 3,1 ± 0,2 | 0,19078 | OK |
| F:lesungen | (χ² = 27, die Regel allein 387) | 27,053 | OK |
| F:lesungen | (χ² = 27, die Regel allein 387) | 387,06 | OK |
| F:lesungen | bei den blauen f = 1,05 ± 0,12 | 1,0513 | OK |
| F:lesungen | bei den blauen f = 1,05 ± 0,12 | 0,12386 | OK |
| 03 | beide aus dem Auslassen je eines Haufens auf ± 0,8 und ± 2,5 genau | 0,76781 | OK |
| 03 | beide aus dem Auslassen je eines Haufens auf ± 0,8 und ± 2,5 genau | 2,5211 | OK |
| 03 | Um die Spitze ist die mittlere lineare Faltung innen flach | 1 | OK |
| 03 | und die inneren Schalen kehren zugleich um | 1 | OK |
| 03 | Fallen sie radial, wird das Profil innen steiler als NFW | 1 | OK |
| 03 | wird das Profil innen flacher als NFW | 1 | OK |
| 07 | An SPARC (2688 Punkte) | 2688 | OK |
| 07 | sind das höchstens 2·10⁻⁵ dex | 1,9094 | OK |
| 07 | Δχ² = 0,01 an SPARC | 0,010607 | OK |
| 07 | 0,06 an HIFLUGCS | 0,060328 | OK |
| 07 | als Faltung 1·10⁻⁴ dex | 1,0937 | OK |
| 07 | unter 1/6 | 0,16667 | OK |
| 07 | Sonne 4,5 | 4,5283 | OK |
| 07 | Milchstraße 2,6 | 2,5521 | OK |
| 08 | Steigung 1 hält die Spreizung GM/(Rc²) | 1 | OK |
| 08 | ½ die Neigung GM/R² | 1 | OK |
| 08 | ⅓ die Gezeit GM/R³ und mit ihr die mittlere Dichte | 1 | OK |
| 08 | im Tiefbereich ¼ den Umlauf, Ω⁴ = GMa₀/R⁴ | 1 | OK |
| 08 | im Tiefbereich ¼ den Umlauf, Ω⁴ = GMa₀/R⁴ | 1 | OK |
| 08 | Steigung 0 hält eine Zeit, ℓ₀ = c·τ_D | 1 | OK |
| 08 | Auf −1 hält λ̄ die Wirkung M·λ̄·c = ħ | 1 | OK |
| 08 | auf −½ hält r̃ die Zeit der eigenen Quantenbewegung, M·r̃²/ħ = τ_D | 1 | OK |
| 08 | liegt bei der Masse ihres Treffpunkts mal Z^{1/(k₁ − k₂)} | 1 | OK |
| 08 | mit dem Exponenten 2 zwischen Spreizung und Neigung | 1 | OK |
| 08 | 3/2 zwischen Spreizung und Gezeit | 1 | OK |
| 08 | 6 zwischen Neigung und Gezeit | 1 | OK |
| 08 | 4 zwischen Neigung und Umlauf im Tiefbereich | 1 | OK |
| 08 | in m_p·(α/α_G)^{3/2} bei δ = α²·m_e/m_p | 1 | OK |
| 08 | darum liegt ihr Treffpunkt für jeden Halter bei M_* | 1 | OK |
| 08 | M_V = √5·m_p·(λ_H/(12Δs_p))^{3/2} | 1 | OK |
| 08 | M_V = √5·m_p·(λ_H/(12Δs_p))^{3/2} | 1 | OK |
| 08 | fällt dagegen die Masse des Sterns heraus | 1 | OK |
| 08 | auf der Kreisbahn 2·m_p·(λ_H/(18Δs_p))^{3/2} | 1 | OK |
| 08 | parallel zu den Kanten der Reihe nach M_U | 1 | OK |
| 08 | der Reihe nach M_U, M_min | 1 | OK |
| 08 | M_min und m_Pl | 1 | OK |
| 08 | parallel zur Seitenhalbierenden von λ̄ nach r* M_c | 1 | OK |
| 08 | zu der von 3Δs nach r̃ m_Pl²/M_c | 1 | OK |
| 08 | M_U = c⁴/(G·a₀) ohne ħ | 1 | OK |
| 08 | M_min = ħa₀/c³ ohne G | 1 | OK |
| 08 | m_Pl = √(ħc/G) ohne a₀ | 1 | OK |
| 08 | GMa₀/c⁴ = M/M_U | 1 | OK |
| 08 | nur um einen Faktor der Form | 1 | OK |
| 08 | eine Dekade der Masse ist in x 0,016 | 0,016202 | OK |
| 04 | g ist darin die halbe Steigung von c²·(dτ/dt)² über r | 1 | OK |
| 04 | Sie spürt die Neigung als g/(dτ/dt) | 1 | OK |
| 04 | Im gleichmäßig gehaltenen Körper ist T/T_R = 1 + P/(ρc²) | 1 | OK |
| 04 | mindestens 13 % wärmer als an der Hülle | 13,326 | OK |
| 04 | mindestens 13 % wärmer als an der Hülle | 13,326 | OK |
| 01 | Laufen beide parallel, d = n + 2 | 2 | OK |
| 01 | Laufen beide parallel, d = n + 2 | 2 | OK |
| 01 | ihr Halter ist der Umlauf, L²/r³ mit n = 2 | 2 | OK |
| 01 | Die 4 in κ² = (4 − d_eff)·Ω² ist n + 2 | 4 | OK |
| 01 | Die 4 in κ² = (4 − d_eff)·Ω² ist n + 2 | 4 | OK |
| 01 | darum zählt die Rückstellung, d = 0, vierfach | 4 | OK |
| 01 | darum zählt die Rückstellung, d = 0, vierfach | 4 | OK |
| 01 | Nach Newton, d = 3, ist K = (n − 1)·⟨r·g⟩/⟨r²⟩ | 0,5 | OK |
| 01 | Nach Newton, d = 3, ist K = (n − 1)·⟨r·g⟩/⟨r²⟩ | 0,5 | OK |
| 01 | die Bindung je Masse ist −((n − 1)/n)·⟨r·g⟩ | -0,33333 | OK |
| 01 | Gas, n = 2, dreht | 2 | OK |
| 01 | bindet mit der Hälfte | 0,5 | OK |
| 01 | wie die Epizyklen an seiner Hülle, K = GM/R³ | 1 | OK |
| 01 | Der Halter bei c, n = 1 | 1 | OK |
| 01 | läuft parallel zu Newton | 0,011102 | OK |
| 01 | seine Bindung ist null | 0 | OK |
| 01 | 5/18 bei β = ½ | 0,27778 | OK |
| 01 | 5/18 bei β = ½ | 0,27778 | OK |
| 01 | β/2 bei kleinem β | 0,50001 | OK |
| 01 | β/2 bei kleinem β | 0,50006 | OK |
| 01 | n − 1 = β·(4 − 3β)/(8 − 7β) | 1 | OK |
| 01 | d über 3 | 1 | OK |
| 01 | und dort ist K = 0 | -1,8713·10⁻⁵ | OK |
| 01 | Mit der Regel ist d = 3 − χ̄ | 3 | OK |
| 01 | Im Tiefbereich, d = 2 | 2 | OK |
| 01 | hängt r·g nicht am Radius | 0,99998 | OK |
| 01 | bei R_max/∛4 | 0,62996 | OK |
| 01 | wie die letzte stabile Kreisbahn um einen Haufen | 1 | OK |
| 01 | das Vierfache der Gezeitendichte des Ganzen | 4 | OK |
| 08 | Jedes Glied F_i ∝ u^{d_i − 3} | 3 | OK |
| 08 | Σ w_i = 1 | 1 | OK |
| 08 | eine Abweichung schwingt mit (Re τ_M)² = 4 − d_eff | 4 | OK |
| 08 | Newton d = 3 | 3 | OK |
| 08 | die Regel 3 − χ | 3 | OK |
| 08 | die Rückstellung d = 0 mit negativem Anteil | -2,6645·10⁻¹⁵ | OK |
| 08 | die Rückstellung d = 0 mit negativem Anteil | -1 | OK |
| 08 | Die Faltung fügt das Glied 3GMu²/c² hinzu, mit d = 5 | 5 | OK |
| 08 | darum trägt Newton nur 1 − 3δ, das Quadrat der umlaufenden Uhr | 1 | OK |
| 08 | und die Faltung den Rest 3δ | 3 | OK |
| 08 | So ist d_eff = 3 + 6δ | 6 | OK |
| 08 | der Punkt beim Bruchteil 3δ auf der Strecke von 3 nach 5 | 3 | OK |
| 08 | und (Re τ_M)² = 1 − 6δ | 1 | OK |
| 08 | und (Re τ_M)² = 1 − 6δ | 6 | OK |
| 08 | unter 4 die Drehung, bei 4 die Scherung, darüber die Dehnung | 4 | OK |
| 08 | Die Identität liegt bei 3 und 0 | 0 | OK |
| 08 | hält Re τ_M·Re τ_M′ = 2 fest | 2 | OK |
| 08 | sein Spiegelpunkt ist der tiefe Bereich bei 2 | 2 | OK |
| 08 | sein Spiegelpunkt ist der tiefe Bereich bei 2 | 2 | OK |
| 08 | Die halbe Drehung liegt bei 15/4 | 3,75 | OK |
| 08 | Die halbe Drehung liegt bei 15/4 | -2 | OK |
| 08 | beide über 15/4 hinweg auf 4 | 4 | OK |
| 08 | beide über 15/4 hinweg auf 4 | 4 | OK |
| 08 | beide über 15/4 hinweg auf 4 | 3,75 | OK |
| 08 | Innen liegt die Scherung in der Mitte der Strecke | 18 | OK |
| 08 | bei 18Δs tragen Newton und Faltung je die Hälfte | 0,5 | OK |
| 08 | bei 18Δs tragen Newton und Faltung je die Hälfte | 0 | OK |
| 08 | der Photonensphäre 9Δs | 9 | OK |
| 08 | ist d_eff = 5 und (Re τ_M)² = −1 | 5 | OK |
| 08 | ist d_eff = 5 und (Re τ_M)² = −1 | -1 | OK |
| 08 | Newton trägt nichts | 0 | OK |
| 05 | t = E/(|k − 1|·P) auch für k unter eins | 1 | OK |
| 05 | Steigung 1 hält die Zeit E/P fest | 1 | OK |
| 05 | P ∝ M⁻², k = −2 | -2 | OK |
| 05 | und lebt E/(3P) | 3 | OK |
| 05 | An der Grenze wächst sie mit k = 1 | 1 | OK |
| 05 | um den Faktor e je Zeit E/P | 1 | OK |
| 05 | Ẽ^{2 − k} | 2 | OK |
| 05 | kehrt k − 1 um | -1 | OK |
| 05 | und hält die Zeit fest | 1 | OK |
| 05 | Fest bleibt nur das Wachsen an der Grenze | 1 | OK |
| 05 | nimmt ihre Spiegelmasse in derselben Zeit um den Faktor e ab | 1 | OK |
| 05 | wird aus k = −2 gespiegelt 2 − k = 4 | 4 | OK |
| 05 | steigt mit Ṫ ∝ T⁴ | 4 | OK |
| 05 | nach derselben Zeit E/(3P) | 3 | OK |
| 08 | nur bei der Spur 0 oder ±1 | 1 | OK |
| 08 | bei Re τ_M = ¼ und ¾, die Spur null | 4,8986·10⁻¹⁶ | OK |
| 08 | jede Richtung kehrt nach zwei Umläufen wieder | 1,2859·10⁻¹⁵ | OK |
| 08 | die Bahn nahe am Kreis schließt nach vier | 4 | OK |
| 08 | bei ⅙ und ⅓, ⅔ und ⅚, ist die Spur ±1 | 1,7764·10⁻¹⁵ | OK |
| 08 | jede Richtung kehrt nach drei Umläufen wieder | 8,063·10⁻¹⁵ | OK |
| 08 | die Bahn schließt nach drei oder sechs | 1 | OK |
| 08 | zwei solche Folgen | 2 | OK |
| 08 | zwei solche Folgen | 2 | OK |
| 08 | bis zur halben Drehung bei 24Δs | 24 | OK |
| 08 | Re τ_M von 1 bis √2 | 1 | OK |
| 08 | Re τ_M von 1 bis √2 | 2 | OK |
| 08 | durchläuft eine Folge | 1 | OK |
| 08 | erreicht die Identität bei 3/2 nicht | 1,5 | OK |
| 08 | die Dreierdrehung bei a_N = 3,17·a₀ | 3,1677 | OK |
| 08 | den Spiegel bei 1,12·a₀ | 1,1199 | OK |
| 08 | und die Dreierdrehung bei 0,23·a₀ | 0,23349 | OK |
| 08 | d²u/dφ² + u = K_el·W/(p_φ·c)² + α_φ²·u | 1 | OK |
| 08 | Auf dem Kreis ist α_φ genau sein Tempo gegen c, bei jedem Tempo | 1 | OK |
| 08 | das Quadrat der umlaufenden Uhr, 1 − α_φ² | 1 | OK |
| 08 | das Quadrat der umlaufenden Uhr, 1 − α_φ² | 1 | OK |
| 08 | den Rest α_φ² mit d = 4 | 4 | OK |
| 08 | (Re τ_M)² = 1 − α_φ² | 1 | OK |
| 08 | dreht als Rosette um denselben Winkel | 1 | OK |
| 08 | ein Drittel der 3δ der Faltung | 0,33333 | OK |
| 08 | die den Takt 2δ mitträgt | 2 | OK |
| 08 | von 3 nach 4 | 3 | OK |
| 08 | von 3 nach 4 | 4 | OK |
| 08 | ist halb so lang wie die von 3 nach 5 | 0,5 | OK |
| 08 | schwach um πα_φ² vor | 1,0001 | OK |
| 08 | ein Sechstel der Faltung | 0,16667 | OK |
| 08 | Die Scherung liegt bei K_el = p_φ·c | 1 | OK |
| 08 | Die Scherung liegt bei K_el = p_φ·c | 1 | OK |
| 08 | mit dem Drehimpuls ħ bei |Q₁Q₂| = q_P² | 1 | OK |
| 08 | im Treffpunkt der Windung, Zα = 1 | 1 | OK |
| 08 | Darüber dehnt die Abbildung | 1 | OK |
| 08 | So erreicht die Bahn um eine Windung c an ihrer Scherung | 1 | OK |
| 08 | die um eine Faltung erst am Ende ihrer Strecke | 1 | OK |
| 08 | die um eine Faltung erst am Ende ihrer Strecke | 5 | OK |
| 08 | (dτ/dt)² = 1 − 6Δs/r + r_Q·3Δs/r² | 0 | OK |
| 08 | r_Q·3Δs = (ℓ_P·Q_el/q_P)² | 0 | OK |
| 08 | r² − 6Δs·r + a_J² + r_Q·3Δs = 0 | 0 | OK |
| 08 | Die Höhe im Thaleskreis wird die Hypotenuse aus a_J und √(r_Q·3Δs) | 0 | OK |
| 08 | beide Horizonte treffen sich bei j² + r_Q/(3Δs) = 1 | 0 | OK |
| 08 | beide Horizonte treffen sich bei j² + r_Q/(3Δs) = 1 | 0 | OK |
| 08 | r_Q/(3Δs) = (Q_el/q_P)²/α_G | 0 | OK |
| 08 | Ohne Drehung liegt der Treffpunkt bei Q_el/q_P = M/m_Pl | 0 | OK |
| 08 | r_e/(3Δs_e) = α/α_G(m_e) = 4,2·10⁴² | 4,1656 | OK |
| 08 | r_e/(3Δs_e) = α/α_G(m_e) = 4,2·10⁴² | 4,1656 | OK |
| 08 | Ohne Drehung ist M_irr die Faltung innerhalb des Horizonts | 0 | OK |
| 08 | M = M_irr + Q_el²/(8πε₀c²r₊) | 0 | OK |
| 08 | M² = (M_irr + M_el)² + M_rot² | 0 | OK |
| 08 | M_irr·M_el = (m_Pl·Q_el/(2q_P))² | 0 | OK |
| 08 | r₀(M_irr)·r₀(M_el) = r_Q·3Δs | 0 | OK |
| 08 | M_irr = M_el = M/2 | 0 | OK |
| 08 | Dort lässt sich die Hälfte der Ruheenergie herausholen | 0,5 | OK |
| 08 | die Entropie ist ein Viertel der ungeladenen Falte gleicher Masse | 0,25 | OK |
| 08 | ħc·(r₊ − r₋)/(4πk·r₀(M_irr)²) | 0 | OK |
| 08 | ħc·(r₊ − r₋)/(4πk·r₀(M_irr)²) | 0 | OK |
| 08 | ohne Windung mit r₀(M_irr)² = r₀·r₊ | 0 | OK |
| 09 | das Mittel aus Schwelle und Tiefe, a_J² = r*·3Δs | 1 | OK |
| 09 | für die Milchstraße (6·10¹⁰ M☉) 1759 | 1758,8 | OK |
| 09 | ihre Horizonte stehen 89,97° zur Tiefe | 89,967 | OK |
| 09 | für die Milchstraße (6·10¹⁰ M☉) 1759 | 1758,8 | OK |
| 09 | eine Scheibe an der Schwelle mit J/M = 2R_d·v_flach hat √2·c/v_flach | 1 | OK |
| 09 | Ihre Spinzahl fällt mit der Masse wie M^{−1/4} | -0,25 | OK |
| 09 | erreicht eins erst am Treffpunkt M_U | 1 | OK |
| 09 | ist sie 2π/Θ_∞ = 5,1 | 5,0966 | OK |
| 09 | ist sie 2π/Θ_∞ = 5,1 | 5,0966 | OK |
| 03 | im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc | 49,592 | OK |
| 03 | im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc | 49,592 | OK |
| 03 | im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc | 3,0653 | OK |
| 03 | verliert ihren Tiefbereich innerhalb von 0,33 Mpc | 0,32723 | OK |
| 03 | verliert ihren Tiefbereich innerhalb von 0,33 Mpc | 0,32723 | OK |
| 05 | t = (4 − d)·E_kin/((l·(d − 2) + 2)·P) | 1 | OK |
| 05 | t = (4 − d)·E_kin/((l·(d − 2) + 2)·P) | 1 | OK |
| 05 | Nach Newton, d = 3, ist das E/((l + 2)·P) | 1 | OK |
| 05 | Nach Newton, d = 3, ist das E/((l + 2)·P) | 1 | OK |
| 05 | je Faktor m·v_flach² | 1,0002 | OK |
| 05 | jede Bahn fällt in t = E_kin/P | 1 | OK |
| 05 | jede Bahn fällt in t = E_kin/P | 1 | OK |
| F:lesungen | bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache | 3,468 | OK |
| F:lesungen | bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache | 2,3885 | OK |
| F:lesungen | bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache | 1,6534 | OK |
| F:lesungen | sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06 | 5,7964 | OK |
| F:lesungen | sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06 | 0,69364 | OK |
| F:lesungen | sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06 | 0,057905 | OK |
| 03 | a₀ = 1,10·10⁻¹⁰ m/s² (1,09 bis | 1,092 | OK |
| 03 | hier (1 + Q)² = 40,5-mal Chandrasekhar | 40,505 | OK |
| 03 | hier (1 + Q)² = 40,5-mal Chandrasekhar | 40,505 | OK |
| 03 | reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc | 49,592 | OK |
| 03 | reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc | 49,592 | OK |
| 03 | reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc | 51,963 | OK |
| 03 | reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc | 51,963 | OK |
| 03 | und das Verhältnis ist 1,04 | 1,0437 | OK |
| 03 | und das Verhältnis ist 1,04 | 1,0437 | OK |
| 03 | ist es (1 + Q)·(1 − 1/e) = 4,0 | 4,023 | OK |
| 03 | ist es (1 + Q)·(1 − 1/e) = 4,0 | 4,023 | OK |
| 03 | mit γ = 5/3 und Sternen von einem Viertel des Gases, ist q = 0,024 | 0,024317 | OK |
| 03 | mit γ = 5/3 und Sternen von einem Viertel des Gases, ist q = 0,024 | 0,024317 | OK |
| 03 | und das Gas verstärkt sie um 2,5 % | 2,4923 | OK |
| 03 | und das Gas verstärkt sie um 2,5 % | 2,4923 | OK |
| 03 | schaukelt sie erst bei γ_P = 4π²·γ/3 = 22 auf | 21,932 | OK |
| 03 | schaukelt sie erst bei γ_P = 4π²·γ/3 = 22 auf | 21,932 | OK |
| 02 | die Erde bis 1,50·10⁶ km | 1,4966 | OK |
| 02 | die Erde bis 1,50·10⁶ km | 1,4966 | OK |
| 02 | bei der mittleren Dichte der Erde 2,4 Stunden | 2,4349 | OK |
| 02 | bei der mittleren Dichte der Erde 2,4 Stunden | 2,4349 | OK |
| 02 | die kürzer läuft als 4,8 Stunden | 4,8152 | OK |
| 02 | die kürzer läuft als 4,8 Stunden | 4,8152 | OK |
| 02 | um jeden Planeten 3·√(π/(Gρ)) | 1 | OK |
| 03 | am Ort der Sonne (8,25 kpc) 1,7 Wasserstoffatome je cm³ | 1,7066 | OK |
| 03 | am Ort der Sonne (8,25 kpc) 1,7 Wasserstoffatome je cm³ | 1,7066 | OK |
| 03 | bei 2R_d 0,57, bei 20 kpc 0,18 | 0,57155 | OK |
| 03 | bei 2R_d 0,57, bei 20 kpc 0,18 | 0,57155 | OK |
| 03 | bei 2R_d 0,57, bei 20 kpc 0,18 | 0,18393 | OK |
| 03 | bei 2R_d 0,57, bei 20 kpc 0,18 | 0,18393 | OK |
| F:lesungen | Bei 10° ist das Achsenverhältnis der Sterne, 0,86, keine Neigung | 0,86 | OK |
| F:lesungen | mindestens um 13 % oval | 12,779 | OK |
| F:lesungen | mindestens um 13 % oval | 12,795 | OK |
| F:lesungen | das 2,1-fache ihrer heutigen Baryonen einmal gebunden | 2,1478 | OK |
| F:lesungen | für jede Masse höchstens das 0,3-fache ihrer Baryonen dazu | 0,29414 | OK |
| F:lesungen | für jede Masse höchstens das 0,3-fache ihrer Baryonen dazu | 0,29414 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 2,9455 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 2,9455 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 3,8459 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 3,8437 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 4,4257 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 4,4011 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 5,2621 | OK |
| 03 | (ν = 2,9 und 3,8) c = 4,4 und 5,3 | 5,2654 | OK |
| F:lesungen | (1 − 1/2π)² = 0,71 | 0,70702 | OK |
| F:lesungen | (1 − 1/2π)² = 0,71 | 0,70702 | OK |
| F:lesungen | ist χ² = 275 statt 307,8 | 275,15 | OK |
| F:lesungen | ist χ² = 275 statt 307,8 | 275,13 | OK |
| F:lesungen | ist χ² = 275 statt 307,8 | 307,78 | OK |
| F:lesungen | Frei liest der Rest dann 0,67 ± 0,06 | 0,67088 | OK |
| F:lesungen | Frei liest der Rest dann 0,67 ± 0,06 | 0,66582 | OK |
| F:lesungen | Frei liest der Rest dann 0,67 ± 0,06 | 0,055936 | OK |
| F:lesungen | Frei liest der Rest dann 0,67 ± 0,06 | 0,055638 | OK |
| F:lesungen | q² liegt Δχ² = 0,4 daneben | 0,41739 | OK |
| F:lesungen | q² liegt Δχ² = 0,4 daneben | 0,43646 | OK |
| F:galaxien | (Rangkorrelation −0,75) | -0,7549 | OK |
| F:galaxien | (Rangkorrelation −0,75) | -0,7549 | OK |
| F:galaxien | an den dichten mitgeeicht (Υ_V = 0,8) | 0,80225 | OK |
| F:galaxien | an den dichten mitgeeicht (Υ_V = 0,8) | 0,79326 | OK |
| F:galaxien | die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33) | -8,6469·10⁻⁴ | OK |
| F:galaxien | die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33) | 0,03336 | OK |
| F:galaxien | die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33) | 32,713 | OK |
| F:galaxien | Newton bleibt damit bei +0,12 und verlangt rund 2 km/s | 0,12273 | OK |
| F:galaxien | Newton bleibt damit bei +0,12 und verlangt rund 2 km/s | 2,1808 | OK |
| F:galaxien | Newton bleibt damit bei +0,12 und verlangt rund 2 km/s | 2,2772 | OK |
| 09 | Ω = 2GJ/(c²r³) = 6Δs·c·a_J/r³ | 1 | OK |
| 09 | 3Δs·c·a_J/(2r³), einem Viertel davon | 0,25 | OK |
| 08 | Die Temperatur ħc·(r₊ − r₋)/(4πk·r₀·r₊) | 1 | OK |
| 08 | Erst bei M_U/4 = 1,4·10²³ M☉ | 1,4355 | OK |
| 08 | Am Horizont ist s² = M_U/(4M) | 1 | OK |
| 08 | kT = m̃c²/(4π). Strahlt | 1 | OK |
| 03 | (r_d = 147,1 Mpc aus ω_b und ω_c) | 147,08 | OK |
| 03 | liest DESI dann a₀ = 1,10 | 1,0997 | OK |
| 03 | m/s² (1,09 bis 1,11 | 1,1073 | OK |
| 03 | mit H₀ = 68,9 bei χ² = 10,4 | 68,933 | OK |
| 03 | bei χ² = 10,4 für 13 Werte | 10,375 | OK |
| 03 | für 13 Werte | 13 | OK |
| 03 | an sechs Rotverschiebungen | 6 | OK |
| 01 | (PubChem, 13,598 eV) | 13,598 | OK |
| 03 | allein (ω_b, ω_c nach Planck), ist a₀ = 1,09 | 1,0885 | OK |
| 03 | (1,07 bis 1,11 bei Δχ² = 1) | 1,0681 | OK |
| 03 | 1,07 bis 1,11 bei Δχ² | 1,1089 | OK |
| 03 | mit H₀ = 68,4: Rotationskurven | 68,442 | OK |
| 02 | nach Z = 2, 10, 18, 36, 54, 86 | 1 | OK |
| 02 | bei N = 16, 20, 28, 50, 82, 126 | 1 | OK |
| 02 | Z = 14, 20, 28, 38, 50, 82 (AME | 1 | OK |
| 02 | Faktor 4,6 bei He | 4,5599 | OK |
| 02 | Halterlänge λ_H = 1,25 fm | 1,2493 | OK |
| 02 | v/c ≈ λ̄_p/λ_H = 0,17 | 0,16835 | OK |
| 02 | Formel zu λ_H = 1,33 fm | 1,3324 | OK |
| 02 | unabhängig davon 0,14 | 0,13741 | OK |
| 02 | angepasst an die 2509 gemessenen Kerne | 2509 | OK |
| 02 | AME2020: x = 0,25 | 0,24768 | OK |
| 02 | (Z = 118, A = 294: x = 0,98) | 0,98153 | OK |
| 02 | (Z = 118, A = 294 | 1 | OK |
| 02 | Die Hülle hält bis Z = 170 bis 171 | 170,06 | OK |
| 02 | Hülle hält bis Z = 170 bis 171: | 170,93 | OK |
| 02 | (36 Paare, AME2020) | 36 | OK |
| 02 | Zα = 1 (Z = 137) gilt für einen punktförmigen Kern | 137 | OK |
| 02 | Mit dem Edelgas Og schließt die siebte Periode bei Z = 118 (PubChem) | 1 | OK |
| 02 | Die Perioden sind 2, 8, 8, 18, 18, 32, 32 lang | 1 | OK |
| 02 | bindet an der Schwelle E = 0 genau bei den Kopplungen N(N + 1) | 1 | OK |
| 02 | Sein Feld folgt der Spiegelform innen auf 0,6 % | 0,55266 | OK |
| 02 | Seine Höhe v = 3,44·Z^{2/3} liegt von Ne bis Og im Fenster | 1 | OK |
| 02 | Seine Zustände an der Schwelle sind Kugelfunktionen, u = √r·P_N^{2l+1}((r − R)/(r + R)) | 1 | OK |
| 01 | Spiegelfeld der Hülle (Kugelfunktionen auf S², Grad n + l, Ordnung 2l + 1) | 1 | OK |
| 02 | erst bei Zα = 0,86 | 0,86109 | OK |
| 02 | (2n² je Schale: 2, 10, 28, 60) | 1 | OK |
| 02 | nach 2, 8, 20, 40, 70, 112 | 1 | OK |
| 02 | also nach 28, 50, 82, 126 (k = 2 bis 5) | 1 | OK |
| 02 | Die nächste Schließung wäre 184 | 184 | OK |
| 02 | schließen nach 2, 10, 18, 36, 54, 86 (PubChem) | 1 | OK |
| 02 | 2,8 bei Rn | 2,7551 | OK |
| 03 | 40 km/s 72 % der | 72,082 | OK |
| 03 | 150 km/s 12 % | 12,363 | OK |
| 03 | Schwelle zu a₀ = 1,11·10⁻¹⁰ | 1,1063 | OK |
| 03 | (16 bis 84 %: 1,00 bis 1,21) | 1,0023 | OK |
| 03 | 1,00 bis 1,21). Die Bose | 1,2125 | OK |
| 03 | die einfache (Δχ² = 53) | 53,128 | OK |
| 03 | Standard-Form (Δχ² = 41) | 40,747 | OK |
| 07 | Dimension 0,92 ± 0,05 | 0,91682 | OK |
| 03 | nach Planck a₀ = 1,065·10⁻¹⁰ | 1,0646 | OK |
| 03 | m/s² H₀ = 67,2 km/s/Mpc | 67,202 | OK |
| 03 | und Ω_DE = 0,685, mit | 0,68466 | OK |
| 03 | 1,1·10⁻¹⁰ H₀ = 68,9 | 68,948 | OK |
| 03 | 68,9 und 0,700 | 0,70043 | OK |
| 03 | erst bei z = −0,10 ein | -0,10209 | OK |
| 03 | = 0,47·v_flach | 0,4714 | OK |
| 03 | unter v_flach = 17 bis 21 | 16,971 | OK |
| 03 | 17 bis 21 km/s. In SPARC | 21,213 | OK |
| 03 | bei z_eq = 3402 | 3402 | OK |
| 03 | bei 3 bis 4·10¹² | 3,1173 | OK |
| 03 | 3 bis 4·10¹² M☉ (mit | 3,8339 | OK |
| 03 | Bereiche ab 6,4·10¹⁰ | 6,3977 | OK |
| 03 | Bereiche ab 6,4·10¹⁰ | 6,4002 | OK |
| 03 | bis 5,2·10¹² M☉ | 5,1822 | OK |
| 03 | bei M > 9·σ⁴/(Ga₀), trägt es sein Feld selbst, bei M > (81/4)·σ⁴/(Ga₀) | 1 | OK |
| 03 | M_Gas = 0,6 bis 1,4·10⁷ M☉ | 0,58961 | OK |
| 03 | 1,4·10⁷ M☉. Kaltes | 1,4395 | OK |
| 03 | zum Endwert 1,23 | 1,2328 | OK |
| 03 | drei Teile tun es bei z = 0,30 | 0,29501 | OK |
| 03 | Mit drei Teilen in H_∞ = c·√Λ_s = 55,6 km/s/Mpc | 55,606 | OK |
| 03 | 9,1 km/s bei 10⁷ M☉ | 9,1295 | OK |
| 03 | 19,3 km/s bei 2·10⁸ M☉ | 19,307 | OK |
| 03 | 10¹² M☉/Mpc: 2,6 kpc | 2,6299 | OK |
| 03 | (ν = 1,007) | 1,0068 | OK |
| 03 | (ν = 3,3) | 3,3035 | OK |
| F:galaxien | Mit Υ_V = 0,91, geeicht | 0,91476 | OK |
| F:galaxien | an den dichten (g_N > 3·a₀, 28) | 28 | OK |
| F:galaxien | dünnen (g_N < a₀, 17) | 17 | OK |
| F:galaxien | mit der Regel +0,10 ± 0,03 dex über | 0,099792 | OK |
| F:galaxien | mit der Regel +0,10 ± 0,03 dex über | 0,03336 | OK |
| F:galaxien | ohne Halo, +0,41 ( | 0,41462 | OK |
| F:galaxien | (Υ_V = 1,08, χ² 54 | 1,081 | OK |
| F:galaxien | (Υ_V = 1,08, χ² 54 | 53,732 | OK |
| F:galaxien | χ² 54 gegen 235) | 235 | OK |
| F:galaxien | Die gut gemessenen (11) | 11 | OK |
| F:galaxien | liegen +0,03 gegen +0,34 | 0,033561 | OK |
| F:galaxien | liegen +0,03 gegen +0,34 | 0,33698 | OK |
| F:galaxien | 10 kpc vom Zentrum (9) +0,06 | 9 | OK |
| F:galaxien | 10 kpc vom Zentrum (9) +0,06 | 0,062422 | OK |
| F:galaxien | (9) +0,06 gegen +0,43 | 0,43249 | OK |
| F:galaxien | Database (52 mit gemessener | 52 | OK |
| F:galaxien | v_esc(8,3 kpc) = 528 (+24, −25) km/s | 528 | OK |
| F:galaxien | v_esc(8,3 kpc) = 528 (+24, −25) km/s | 24 | OK |
| F:galaxien | v_esc(8,3 kpc) = 528 (+24, −25) km/s | 25 | OK |
| F:galaxien | v_esc(8,3 kpc) = 528 (+24, −25) km/s | 8,3 | OK |
| F:galaxien | bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter | 0,4231 | OK |
| F:galaxien | bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter | 510,98 | OK |
| F:galaxien | bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter | 510,98 | OK |
| F:galaxien | bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter | 0,68062 | OK |
| F:galaxien | Bis zur größten Umkehr gezählt sind es 590 km/s | 590,08 | OK |
| F:galaxien | für ihre Masse von 10¹² M☉ | 1 | OK |
| F:galaxien | mit γ = 9/2, der Flussdimension der Kugel | 4,5 | OK |
| F:galaxien | in jede Richtung, (σ/c)² = (2/9)·(v_flach/c)² | 0,22222 | OK |
| F:galaxien | in jede Richtung, (σ/c)² = (2/9)·(v_flach/c)² | 0,22223 | OK |
| F:galaxien | Sternhaufen im Median bei R_h = 0,38·r* | 0,37906 | OK |
| F:galaxien | Sternhaufen im Median bei R_h = 0,38·r* | 0,37906 | OK |
| F:galaxien | 38 der 52 innerhalb | 38 | OK |
| F:galaxien | 38 der 52 innerhalb | 52 | OK |
| F:galaxien | die Zwerge bei 16·r*, 76 der 77 außerhalb | 15,825 | OK |
| F:galaxien | die Zwerge bei 16·r*, 76 der 77 außerhalb | 15,825 | OK |
| F:galaxien | die Zwerge bei 16·r*, 76 der 77 außerhalb | 76 | OK |
| F:galaxien | die Zwerge bei 16·r*, 76 der 77 außerhalb | 77 | OK |
| F:galaxien | (Rangkorrelation −0,51 an 14 Haufen | -0,51209 | OK |
| F:galaxien | (Rangkorrelation −0,51 an 14 Haufen | 14 | OK |
| F:galaxien | am heutigen Abstand zu 0,13 bis 0,64 füllt | 0,13104 | OK |
| F:galaxien | am heutigen Abstand zu 0,13 bis 0,64 füllt | 0,63692 | OK |
| F:galaxien | zu 0,13 bis 0,64 füllt (+0,31) | 0,30549 | OK |
| 03 | a₀/(πG) = 242 M☉/pc² | 242,07 | OK |
| 03 | a_N = 2,54·a₀ | 2,5396 | OK |
| 03 | beträgt 0,648·a₀ | 0,64761 | OK |
| 07 | Rosetten mit 254,6° | 254,56 | OK |
| 07 | e^{−7480} | 7479,6 | OK |
| 08 | bei der Milchstraße 0,42″ | 0,41894 | OK |
| 08 | Für die eingeschlossene Masse im Bild ist am Ring deshalb G·M/R_E² = c²/(4D_E), für jede Masse und jede runde Verteilung | 1 | OK |
| 08 | Für die eingeschlossene Masse im Bild ist am Ring deshalb G·M/R_E² = c²/(4D_E), für jede Masse und jede runde Verteilung | 1 | OK |
| 08 | Für einen Punkt ist im Newton-Bereich R_E² = 2r₀·D_E wie r*² = r₀·Λ, also (R_E/r*)² = 2D_E/Λ | 1 | OK |
| 08 | In Flächendichte gelesen sind Ring und Schwelle dieselbe Form, g = πG·Σ̄ | 1 | OK |
| 08 | In Flächendichte gelesen sind Ring und Schwelle dieselbe Form, g = πG·Σ̄ | 1 | OK |
| 08 | Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06) | 5,1476 | OK |
| 08 | Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06) | 5,1476 | OK |
| 08 | Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06) | 0,98886 | OK |
| 08 | Jeder Ring liegt innerhalb der Schwelle der Masse, die er einschließt, bei G·M/R_E² ≥ 5,3·a₀ | 5,338 | OK |
| 08 | Jeder Ring liegt innerhalb der Schwelle der Masse, die er einschließt, bei G·M/R_E² ≥ 5,3·a₀ | 5,338 | OK |
| 08 | ε_L = ℓ_P/ℓ₀ = 1,9·10⁻⁶² | 1,9062 | OK |
| 08 | ε_L^{1/3}·m_Pl = 5,8·10⁻²⁹ kg | 5,8139 | OK |
| 08 | r̃ ist die Schwelle der Spiegelmasse, r̃ = √2·r*(m̃) | 1 | OK |
| 08 | und r*/r̃ = M/m_Pl | 1 | OK |
| 08 | Beide kreuzen sich bei m_Pl in der Länge des Vakuumteils | 1 | OK |
| 08 | m_Pl²/M_c = ε_L^{−1/3}·m_Pl = 8,1·10¹² kg | 8,1474 | OK |
| 08 | m_Pl²/M_c = ε_L^{−1/3}·m_Pl = 8,1·10¹² kg | 8,1474 | OK |
| 08 | mit 3Δs = r̃, in derselben Länge ℓ_P^{2/3}·ℓ₀^{1/3} | 1 | OK |
| 08 | mit 3Δs = r̃, in derselben Länge ℓ_P^{2/3}·ℓ₀^{1/3} | 1 | OK |
| 08 | liegen Teilchen und Atome unter der Linie r̃ (Wasserstoff R/r̃ = 10⁻¹⁶) | -16 | OK |
| 08 | Erde und Galaxie darüber (10²⁷ und 10⁴⁹) | 27 | OK |
| 08 | Erde und Galaxie darüber (10²⁷ und 10⁴⁹) | 49 | OK |
| 08 | kreuzt r̃ bei 2,0·10⁻⁸ kg und r* bei 3,1·10⁻⁸ kg | 2,0221 | OK |
| 08 | kreuzt r̃ bei 2,0·10⁻⁸ kg und r* bei 3,1·10⁻⁸ kg | 3,1442 | OK |
| 08 | Leichtere Körner liegen wie die Atome unter r̃ und über r* | 1 | OK |
| 08 | schwerere wie Erde und Stern über r̃ und unter r* | 1 | OK |
| 08 | ihre Summe ist die Tiefe, ihr Produkt die Verwindung | 1 | OK |
| 08 | ihre Summe ist die Tiefe, ihr Produkt die Verwindung | 1 | OK |
| 08 | und liegen auf dem Kreis |r±| = a_J um den Knoten, im Winkel cos φ = 1/j zur Tiefe | 1 | OK |
| 08 | und liegen auf dem Kreis |r±| = a_J um den Knoten, im Winkel cos φ = 1/j zur Tiefe | 1 | OK |
| 08 | (Erde φ = 89,92°) | 89,922 | OK |
| 09 | die Spinzahl ist diese Zahl durch die Kopplung, j = (J/ħ)/α_G | 1 | OK |
| 09 | Die Drehung erreicht c, wenn die Zahl ihrer Quanten α_G erreicht, J = GM²/c | 1 | OK |
| 01 | denn bei λ̄/α_G schließt die Quantenlänge nach einem Umlauf um die eigene Faltung | 1 | OK |
| 01 | denn bei λ̄/α_G schließt die Quantenlänge nach einem Umlauf um die eigene Faltung | 1 | OK |
| 01 | (beim Proton 1,15 Mpc) | 1,154 | OK |
| 08 | Zugleich ist r̃² = ħ·τ_D/M | 1 | OK |
| 08 | r̃ = √(2D·τ_D), wie seine Faltung mit D = GM/c bis r* = √(D·τ_D) | 1 | OK |
| 07 | Eine Ausdehnung mit der Rate H gibt H·τ_D = 2π·T_H/T_I = 2π/Θ | 1 | OK |
| 07 | Eine Ausdehnung mit der Rate H gibt H·τ_D = 2π·T_H/T_I = 2π/Θ | 1 | OK |
| 07 | mit der Temperatur des Horizonts kT_H = ħH/(2π) | 1 | OK |
| 08 | Beim Elektron ist λ̄_e = 510·r* | 509,88 | OK |
| 08 | sein Feld liegt außerhalb der Quantenlänge ganz unter a₀ | 1 | OK |
| 08 | beim Proton r* = 154·λ̄_p | 154,31 | OK |
| 08 | Unter m_Pl²/M_c reicht die eigene Quantenbewegung weiter als die Tiefe, r̃ > 3Δs | 1 | OK |
| 08 | Ihr eigenes Spiegelbild, m̃ = M, ist die Falte bei m_Pl/√2, wo r₀ = λ̄ ist | 1 | OK |
| 08 | Ihr eigenes Spiegelbild, m̃ = M, ist die Falte bei m_Pl/√2, wo r₀ = λ̄ ist | 1 | OK |
| 08 | bei der ein Teilchen mit Spin ½ die Spinzahl j = 1 erreicht | 1 | OK |
| 08 | seine Baryonen tragen (3/4π)·M_U, mit r̃ = 2,05·ℓ_P | 2,0467 | OK |
| 08 | Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1 | 1 | OK |
| 08 | Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1 | 1 | OK |
| 08 | Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1 | 1 | OK |
| 08 | r* (½, ½) und r̃ (−½, ½) | 1 | OK |
| 08 | r* (½, ½) und r̃ (−½, ½) | 1 | OK |
| 08 | der Schwerpunkt ist ℓ_P^{2/3}·ℓ₀^{1/3} (0, ⅓) | 1 | OK |
| 08 | darum teilt M_c die Leiter im Drittel | 1 | OK |
| 08 | darum teilt M_c die Leiter im Drittel | 1 | OK |
| 08 | Die Treffmassen sind x = ±1, wo λ̄ und 3Δs die Ecke ℓ₀ erreichen | 1 | OK |
| 08 | Die Treffmassen sind x = ±1, wo λ̄ und 3Δs die Ecke ℓ₀ erreichen | 1 | OK |
| 08 | λ̄·3Δs = ℓ_P² und ebenso r*·r̃ = ℓ_P·ℓ₀ | 1 | OK |
| 08 | λ̄·3Δs = ℓ_P² und ebenso r*·r̃ = ℓ_P·ℓ₀ | 1 | OK |
| 08 | um λ̄ die Tiefe auf λ̄/α_G (−3, 0) | 1 | OK |
| 08 | Tully-Fisher (1, −1) und α_I (−1, −1) | 1 | OK |
| 08 | Tully-Fisher (1, −1) und α_I (−1, −1) | 1 | OK |
| 08 | r_e, λ̄_e und a_B sind n = 1, 0 und −1 | 1 | OK |
| 08 | α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)² | 1 | OK |
| 08 | α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)² | 1 | OK |
| 08 | α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)² | 1 | OK |
| 08 | Zwei Spiegel nacheinander strecken jede Länge um das Quadrat des Abstands ihrer Mitten | 1 | OK |
| 01 | Am Ort gelesen, gegen den Abstand r, wird es zum Nome: s = r*/r | 1 | OK |
| 08 | Das Elektron liegt bei x = −0,36, links von M_c, das Proton bei −0,31 | -0,36258 | OK |
| 08 | Das Elektron liegt bei x = −0,36, links von M_c, das Proton bei −0,31 | -0,3097 | OK |
| 08 | Erde und Sonne bei 0,53 und 0,62 | 0,52557 | OK |
| 08 | Erde und Sonne bei 0,53 und 0,62 | 0,61505 | OK |
| 08 | die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99 | 0,78968 | OK |
| 08 | die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99 | 0,85808 | OK |
| 08 | die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99 | 0,98992 | OK |
| 08 | r* und r̃ sind seine Diagonalen und kreuzen sich in der Mitte | 1 | OK |
| 01 | Um den Ort gespiegelt wird die Tiefe zur Länge der Neigung, g = c²/(r²/3Δs) | 1 | OK |
| 01 | die Schwelle ist der Ort, an dem dieses Spiegelbild die Schließung ℓ₀ erreicht | 1 | OK |
| 08 | Aus M, ħ, G, c und a₀ lassen sich nur zwei reine Zahlen bilden, M/m_Pl und ε_L | 2 | OK |
| 08 | der Horizont 2·3Δs, die Wärmelänge der Falte 8π·3Δs und c mal ihrer Lebensdauer 5120π·α_G·3Δs | 1 | OK |
| 08 | der Horizont 2·3Δs, die Wärmelänge der Falte 8π·3Δs und c mal ihrer Lebensdauer 5120π·α_G·3Δs | 1 | OK |
| 01 | Kreuzung in der Mitte nach dem ganzen, ω·t = π/(2√2) nach der Umkehr | 1 | OK |
| 01 | der gebundene Bereich als rollender Kreis, r = (r_U/2)·(1 − cos θ) | 1 | OK |
| 03 | Mit g = (1 + Q)·g_bar ist am Rand R_Δ jedes Haufens a_N = Δ·H²·R/(2(1 + Q)) | 1 | OK |
| 03 | R = 2(1 + Q)·(ln q)²·a₀/(Δ·H²), mit 1 + Q = 2π also R = 2Θ·(ln q)²·(c/H)/Δ, für Δ = 200 heute 1,3 Mpc | 1 | OK |
| 03 | R = 2(1 + Q)·(ln q)²·a₀/(Δ·H²), mit 1 + Q = 2π also R = 2Θ·(ln q)²·(c/H)/Δ, für Δ = 200 heute 1,3 Mpc | 1,3258 | OK |
| 03 | Dort liegt die Masse M = 8·(ln q)⁶·Θ³·M_H/Δ² = 2,5·10¹⁴ M☉ | 2,5239 | OK |
| 03 | Dort liegt die Masse M = 8·(ln q)⁶·Θ³·M_H/Δ² = 2,5·10¹⁴ M☉ | 2,5239 | OK |
| 01 | ohne Rückstellung linear δ_c = (3/20)·(12π)^{2/3} = 1,686 | 1,6865 | OK |
| 01 | ohne Rückstellung linear δ_c = (3/20)·(12π)^{2/3} = 1,686 | 1,6865 | OK |
| 01 | und auf das 18π²-fache des Mittels verdichtet | 1 | OK |
| 01 | Druck des Lichts p = u/3 | 1 | OK |
| 03 | R_max = v_flach/H_∞, mit H_∞ = c·√Λ_s = 55,6 km/s/Mpc | 55,606 | OK |
| 03 | R_max = v_flach/H_∞, mit H_∞ = c·√Λ_s = 55,6 km/s/Mpc | 55,606 | OK |
| 03 | für die Milchstraße (6·10¹⁰ M☉) sind es 3,1 Mpc | 3,0653 | OK |
| 03 | für die Milchstraße (6·10¹⁰ M☉) sind es 3,1 Mpc | 3,0675 | OK |
| 03 | R_max³ = (1 + Q)·GM_b/H_∞², bei M₂₀₀ = 10¹⁵ M☉ 11 Mpc | 11,163 | OK |
| 03 | R_max³ = (1 + Q)·GM_b/H_∞², bei M₂₀₀ = 10¹⁵ M☉ 11 Mpc | 11,163 | OK |
| 03 | denn −d(ln n_B)/ds = ν und damit dΦ/dr = ν·a_N | 1 | OK |
| 03 | denn −d(ln n_B)/ds = ν und damit dΦ/dr = ν·a_N | 1 | OK |
| 03 | Newton allein ist v_flach²·ln q, der Logarithmus des Rückkehranteils | 1 | OK |
| 03 | Die Rückkehr fügt den Logarithmus der Summe der Rückläufe hinzu, q·ν = n_B | 1 | OK |
| 03 | Die Bindung zwischen zwei Orten ist so v_flach² mal dem Logarithmus des Verhältnisses ihrer Rückkehrzahlen | 1 | OK |
| 03 | mit n_B(r*/R_max) ≈ R_max/r* = (√15/2π²)·c/v_flach | 0,99855 | OK |
| 03 | an der Schwelle 3,43·v_flach = 585 km/s | 3,4304 | OK |
| 03 | an der Schwelle 3,43·v_flach = 585 km/s | 3,4308 | OK |
| 03 | an der Schwelle 3,43·v_flach = 585 km/s | 584,7 | OK |
| 03 | am Ort der Sonne (8,25 kpc) 591 km/s | 590,58 | OK |
| 03 | am Ort der Sonne (8,25 kpc) 591 km/s | 8,249 | OK |
| 03 | √2·v_flach = 241 km/s | 241,05 | OK |
| 03 | √2·v_flach = 241 km/s | 241,05 | OK |
| 03 | Sie teilt 4Ω² zwischen Gezeit mit Umlauf und Epizyklen | 4 | OK |
| 03 | Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc | 3,7934 | OK |
| 03 | Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc | 3,7911 | OK |
| 03 | Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc | 1,9079 | OK |
| 03 | Verliert er neun Zehntel davon als Gas, rückt der Rand auf 3,7 kpc, ohne dieses Gedächtnis auf 2,1 kpc | 3,7419 | OK |
| 03 | Verliert er neun Zehntel davon als Gas, rückt der Rand auf 3,7 kpc, ohne dieses Gedächtnis auf 2,1 kpc | 2,1262 | OK |
| 03 | also r_J = v_flach/(√d_eff·Ω), linear im Tempo | 0,99249 | OK |
| 03 | gibt die Jeans-Gleichung σ² = v_flach²/γ an jedem Radius | 1,0005 | OK |
| 03 | hat nur mit σ < v_flach/√3 = 0,58·v_flach eine endliche Masse | 0,57736 | OK |
| 03 | Trägt sich eine isotherme Kugel ganz im Tiefbereich selbst, ist σ⁴ = v_flach⁴/γ² = (4/81)·GMa₀, also γ = 9/2 | 4,5 | OK |
| 03 | Sie fällt außen wie r^{−9/2}, Virialsatz und Profil lesen dieselbe Zahl | 4,5 | OK |
| 03 | Sie fällt außen wie r^{−9/2}, Virialsatz und Profil lesen dieselbe Zahl | 1 | OK |
| 03 | Darum ist 9/2 die Obergrenze, mit Newtonschem Kern fällt die Kugel flacher | 1 | OK |
| 03 | Sein Spiegelbild um r* ist c_s²/a₀ = ℓ₀·(c_s/c)² | 1 | OK |
| 03 | Das Bahntempo der Regel fällt nie unter v_flach | 1 | OK |
| 03 | f(S) = δ_c/√(2π)·S^{−3/2}·e^{−δ_c²/(2S)} hat die Form T^{−3/2}·e^{−τ_D/(4T)}, mit S/(2δ_c²) statt T/τ_D | 1 | OK |
| 03 | Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c, der typisch bindenden Masse M_z, und dort sind es 0,317 | 0,31731 | OK |
| 03 | Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c, der typisch bindenden Masse M_z, und dort sind es 0,317 | 0,31731 | OK |
| 03 | Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c | 1 | OK |
| 03 | Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c | 1 | OK |
| 03 | Ebenso hat die Rückkehr je logarithmischem Schritt der Zeit ihren Gipfel bei T = τ_D/2 | 2 | OK |
| 03 | Der Faktor 2 von Press und Schechter ist der Spiegel an der Schwelle | 2 | OK |
| 03 | Ohne ihn fände über alle Skalen nur die Hälfte der Baryonen einen Bereich | 0,5 | OK |
| F:lokale_gruppe | liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database) | 783,8 | OK |
| F:lokale_gruppe | liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database) | 90,18 | OK |
| F:lokale_gruppe | liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database) | 783,8 | OK |
| F:lokale_gruppe | liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database) | 90,18 | OK |
| F:lokale_gruppe | Newton verlangt für 90 km/s das 19-fache der Baryonen | 18,573 | OK |
| F:lokale_gruppe | aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉) | 6,6051 | OK |
| F:lokale_gruppe | aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉) | 6,595 | OK |
| F:lokale_gruppe | aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉) | 1,4921 | OK |
| F:lokale_gruppe | aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉) | 1,4913 | OK |
| F:lokale_gruppe | v₂ = 214 km/s aus den Baryonen der Bereiche | 214,04 | OK |
| F:lokale_gruppe | v₂ = 214 km/s aus den Baryonen der Bereiche | 214,04 | OK |
| F:lokale_gruppe | ist das Paar radial seit 13,1 Gyr auf der Bahn | 13,116 | OK |
| F:lokale_gruppe | ist das Paar radial seit 13,1 Gyr auf der Bahn | 13,116 | OK |
| F:lokale_gruppe | mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr | 13,475 | OK |
| F:lokale_gruppe | mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr | 13,475 | OK |
| F:lokale_gruppe | mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr | 13,81 | OK |
| 03 | also H = H_∞·coth(3H_∞t/2) und Θ = Θ_∞·tanh(3H_∞t/2) mit Θ_∞ = √15/π | 1 | OK |
| 03 | Mit H_∞·τ_D = 2π/Θ_∞ ist die Zeit | 1 | OK |
| 03 | bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95 | 1,0201 | OK |
| 03 | bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95 | 13,833 | OK |
| 03 | bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95 | 13,828 | OK |
| 03 | bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95 | 0,9507 | OK |
| 03 | bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95 | 0,95039 | OK |
| 03 | im Formteil ½ − (3/2)·(Θ/Θ_∞)², heute −0,53 | -0,52699 | OK |
| 03 | im Formteil ½ − (3/2)·(Θ/Θ_∞)², heute −0,53 | -0,52694 | OK |
| 03 | bei Θ = Θ_∞/√3 = 0,712 | 0,71176 | OK |
| 03 | bei Θ = Θ_∞/√3 = 0,712 | 0,71188 | OK |
| 03 | Formteil und Faltung sind gleich bei Θ_∞/√2 = 0,872 | 0,87173 | OK |
| 03 | Formteil und Faltung sind gleich bei Θ_∞/√2 = 0,872 | 0,87164 | OK |
| 03 | Beide Stufen hängen nur an Θ, nicht an a₀ | 1 | OK |
| 03 | dauert das Fenster 0,68 Mrd. Jahre | 0,67735 | OK |
| 03 | im Formteil heute 0,51·L und 0,18·L | 0,50765 | OK |
| 03 | im Formteil heute 0,51·L und 0,18·L | 0,50765 | OK |
| 03 | im Formteil heute 0,51·L und 0,18·L | 0,18363 | OK |
| 03 | im Formteil heute 0,51·L und 0,18·L | 0,18363 | OK |
| 03 | c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70 | 0,70092 | OK |
| 03 | c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70 | 0,70092 | OK |
| 03 | c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70 | 0,70092 | OK |
| 03 | Mit der Strahlung ist es 0,69 | 0,69128 | OK |
| 03 | Mit der Strahlung ist es 0,69 | 0,69128 | OK |
| 03 | Der Radius der letzten Streufläche wächst dabei auf 0,68·L | 0,68124 | OK |
| 03 | Der Radius der letzten Streufläche wächst dabei auf 0,68·L | 0,68125 | OK |
| 03 | die drei Kreispaare erscheinen und wachsen bis höchstens α = 43° | 42,781 | OK |
| 03 | die drei Kreispaare erscheinen und wachsen bis höchstens α = 43° | 42,815 | OK |
| 03 | Ihr Durchmesser, 1,36·L, bleibt unter dem Abstand √2·L der Abbilder über die Diagonalen der Flächen | 1,3625 | OK |
| 03 | Ihr Durchmesser, 1,36·L, bleibt unter dem Abstand √2·L der Abbilder über die Diagonalen der Flächen | 1 | OK |
| 03 | = (5/6)·B(5/6, 2/3)·a_× = 1,44·a_× | 1,4373 | OK |
| 03 | = (5/6)·B(5/6, 2/3)·a_× = 1,44·a_× | 1,4373 | OK |
| 03 | Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 % | 60,703 | OK |
| 03 | Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 % | 60,703 | OK |
| 03 | Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 % | 70,994 | OK |
| 03 | Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 % | 70,997 | OK |
| 03 | Linear wächst Struktur von heute an im Formteil nur noch um 41 % | 40,858 | OK |
| 03 | Linear wächst Struktur von heute an im Formteil nur noch um 41 % | 40,85 | OK |
| 03 | denn (a_×/a)³ = (Θ_∞/Θ)² − 1 mit Θ_∞ = √15/π | 1 | OK |
| 08 | Bei 6Δs steht die ruhende Uhr | 6 | OK |
| 08 | Bei 9Δs steht die umlaufende | 9 | OK |
| 08 | dort läuft die Bahn für die ruhende Uhr am Ort mit c | 1 | OK |
| 08 | Bei 18Δs steht die radiale Schwingung | 18 | OK |
| 08 | die Bahn läuft dort für die ruhende Uhr mit c/2 | 2 | OK |
| 08 | am kleinsten bei 18Δs und bei 12Δs = 4d genau eins | 18 | OK |
| 08 | am kleinsten bei 18Δs und bei 12Δs = 4d genau eins | 12 | OK |
| 08 | und läuft gegen die Fernen je Umlauf um 2π·(1 − √(1 − 3δ)) vor, schwach 3πδ | 3 | OK |
| 08 | und läuft gegen die Fernen je Umlauf um 2π·(1 − √(1 − 3δ)) vor, schwach 3πδ | 1 | OK |
| 08 | In erster Ordnung ist Re τ_M = 1 − 3δ: Um 3δ verfehlt die Bahn den vollen Umlauf, die Drehung des Perihels | 1 | OK |
| 08 | An der Photonensphäre ist (Re τ_M)² = −1 | -1 | OK |
| 08 | mit b_c = 9√3·Δs, dem Rand des Schattens | 9 | OK |
| 08 | sein Umlauf c/(2π·b_c) bei einer Sonnenmasse 6,2 kHz | 6,2185 | OK |
| 08 | sein Umlauf c/(2π·b_c) bei einer Sonnenmasse 6,2 kHz | 6,2185 | OK |
| 01 | Wo es 2r₀ erreicht, liegt der Horizont der Falte | 1 | OK |
| 01 | wo es 4D_E erreicht, der Einsteinring | 1 | OK |
| 05 | q = e^{−r*/R} ≈ 0, Sonne r*/R = 1,6·10⁶ | 1,6084 | OK |
| 05 | δ(R₂₀₀) = 2,4·10⁻⁵ (X-COP) | 2,3812 | OK |
| 08 | bis zur Naht, 5,8·r*, ist sie 0,025 der Höhe breit | 5,7687 | OK |
| 08 | bis zur Naht, 5,8·r*, ist sie 0,025 der Höhe breit | 0,025264 | OK |
| F:lesungen | Mit θ̃ = H ist θ̃·τ_D = 2π/Θ = 6,16 | 6,1595 | OK |
| F:lesungen | Mit θ̃ = H ist θ̃·τ_D = 2π/Θ = 6,16 | 6,1595 | OK |
| F:lesungen | also ν ≤ 1,09 an jedem Ort | 1,0912 | OK |
| F:lesungen | lägen dann −0,01 dex auf der Rechnung | -0,011039 | OK |
| F:lesungen | haben fünf der sechs weniger Nachbarn | 5 | OK |
| F:lesungen | Ihr Rang liegt im Mittel bei 0,27 statt 0,5, zwei Fehler darunter | 0,26528 | OK |
| F:lesungen | Ihr Rang liegt im Mittel bei 0,27 statt 0,5, zwei Fehler darunter | 1,9917 | OK |
| F:lesungen | in 1 Mpc bei 0,26 | 0,25694 | OK |
| F:lesungen | Nur AGC 219533 liegt dicht (0,87) | 0,87333 | OK |
| F:lesungen | ihr W50 liegt bei gleicher HI-Masse +0,03 ± 0,02 dex wie das der übrigen | 0,02711 | OK |
| F:lesungen | ihr W50 liegt bei gleicher HI-Masse +0,03 ± 0,02 dex wie das der übrigen | 0,021106 | OK |
| 03 | In den Größen des Modells ist es a₀/(√κ·H·v_flach), mit R_ta³ = GM/(κ·H²) | 1 | OK |
| 03 | die Umlaufrate an der Schwelle, v_flach/r* = a₀/v_flach | 1 | OK |
| 03 | also (c/v_flach)·Θ/(2π·√κ), tangential gegen parametrisch | 1 | OK |
| 03 | also (c/v_flach)·Θ/(2π·√κ), tangential gegen parametrisch | 1 | OK |
| 03 | Kehrt er heute um, ist mit der Rückstellung κ = 0,3 | 0,29519 | OK |
| 03 | Kehrt er heute um, ist mit der Rückstellung κ = 0,3 | 0,29519 | OK |
| 03 | mit Δ = (3π/4)² nach dem halben Umlauf, κ = 0,44 | 0,44179 | OK |
| 03 | mit Δ = (3π/4)² nach dem halben Umlauf, κ = 0,44 | 0,44179 | OK |
| 03 | mit Δ = (3π/4)² nach dem halben Umlauf | 1 | OK |
| 08 | bei M_min λ̄ = r̃ = ℓ₀ und r* = ℓ_P | 1 | OK |
| 08 | bei M_min λ̄ = r̃ = ℓ₀ und r* = ℓ_P | 1 | OK |
| 08 | bei M_U 3Δs = r* = ℓ₀ und r̃ = ℓ_P | 1 | OK |
| 08 | Das Dreieck schließt, α_G·α_I = (v_flach/c)⁴ | 1 | OK |
| 08 | α_I = λ̄/ℓ₀ = kT_I/(Mc²) = M_min/M | 1 | OK |
| 08 | α_I = λ̄/ℓ₀ = kT_I/(Mc²) = M_min/M | 1 | OK |
| 08 | ist ihr Mittel √(α_G·α_I) | 1 | OK |
| 08 | Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U | 1 | OK |
| 08 | Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U | 1 | OK |
| 08 | Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U | 1 | OK |
| 08 | bei M_c sind beide Quantenkopplungen gleich, α_G = α_I = ε_L^{2/3} | 1 | OK |
| 08 | bei M_c sind beide Quantenkopplungen gleich, α_G = α_I = ε_L^{2/3} | 1 | OK |
| 08 | kehrt α_G um und tauscht die beiden anderen, α_I(m_Pl²/M) = (v_flach/c)⁴ | 1 | OK |
| 08 | kehrt α_G um und tauscht die beiden anderen, α_I(m_Pl²/M) = (v_flach/c)⁴ | 1 | OK |
| 08 | Beim Elektron ist α_I = 4,6·10⁻⁴⁰ | 4,5544 | OK |
| 04 | kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π) | 1 | OK |
| 04 | kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π) | 1 | OK |
| 04 | kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π) | 1 | OK |
| 04 | Ihre Ruheenergie erreicht die Falte bei M = m_Pl/√(8π) = 4,3·10⁻⁹ kg | 1 | OK |
| 04 | Ihre Ruheenergie erreicht die Falte bei M = m_Pl/√(8π) = 4,3·10⁻⁹ kg | 4,3414 | OK |
| 04 | Am Horizont ist das 4π·r₀: Die Wärme der Falte ist die Wärme ihrer Neigung am Horizont | 1 | OK |
| 04 | Am Horizont ist das 4π·r₀: Die Wärme der Falte ist die Wärme ihrer Neigung am Horizont | 1 | OK |
| 04 | Die Falte erreicht diese Grenze, mit R = r₀ ist sie 4πk·α_G | 1 | OK |
| 04 | Die Falte erreicht diese Grenze, mit R = r₀ ist sie 4πk·α_G | 1 | OK |
| 08 | Zwei gleiche Massen im Abstand a strahlen P = (64/5)·(c⁵/G)·δ⁵ mit δ = 3Δs/a | 1 | OK |
| 08 | ihr Treffpunkt c⁵/G, die Ruheenergie je Lichtzeit der Tiefe, hängt an keiner Masse | 1 | OK |
| 08 | Die Falte selbst strahlt als schwarzer Körper (c⁵/G)/(15360π·α_G) | 1 | OK |
| 08 | lebt sie t = 5120π·α_G·3Δs/c, eine Sonnenmasse 2,1·10⁶⁷ Jahre | 1 | OK |
| 08 | lebt sie t = 5120π·α_G·3Δs/c, eine Sonnenmasse 2,1·10⁶⁷ Jahre | 2,0955 | OK |
| 08 | Die Länge c·t ist der Punkt (3, 0) der Leiter, um ℓ_P gespiegelt zu λ̄/α_G (−3, 0) | 1 | OK |
| 08 | Die Länge c·t ist der Punkt (3, 0) der Leiter, um ℓ_P gespiegelt zu λ̄/α_G (−3, 0) | 1 | OK |
| 03 | S_Λ·Λ_s·ℓ_P² = 3πk wie bei de Sitter | 1 | OK |
| 03 | S_Λ = (45/4π³)·k/ε_L² = 1,0·10¹²³ k | 0,99851 | OK |
| 03 | S_Λ = (45/4π³)·k/ε_L² = 1,0·10¹²³ k | 0,99851 | OK |
| 09 | Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c | 1 | OK |
| 09 | Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c | 1 | OK |
| 09 | Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c | 1 | OK |
| 08 | c³τ_D/G = 1,1·10⁵⁴ kg | 1,1418 | OK |
| 08 | Δs(M_U) = ⅔·Λ | 2 | OK |
| 08 | (4π⁴/15)·ε_L² = 0,94·10⁻¹²² | 0,94388 | OK |
| 02 | r* = 32 fm | 32,453 | OK |
| 02 | a_B = 1631·r* | 1630,6 | OK |
| 04 | Friedmanns H² = 8πG·ρ_krit/3 sagt, dass diese Kugel ihre eigene Falte ist, 6Δs_H = c/H | 1 | OK |
| 04 | Friedmanns H² = 8πG·ρ_krit/3 sagt, dass diese Kugel ihre eigene Falte ist, 6Δs_H = c/H | 1 | OK |
| 04 | Sonne 1,05·10⁷⁷ k | 1,0489 | OK |
| 01 | 0,022 s | 0,021968 | OK |
| 01 | Sonnenoberfläche 67 s | 66,981 | OK |
| 01 | auf der Erdbahn sind es 0,31 s | 0,31149 | OK |
| 01 | Takt 0,81 | 0,8096 | OK |
| 08 | 1,7512″ | 1,7512 | OK |
| F:ganzes | γ − 1 = (2,1 ± 2,3)·10⁻⁵ | 2,1 | OK |
| F:ganzes | (2,1 ± 2,3)·10⁻⁵ | 2,3 | OK |
| F:ganzes | (gemessen 43,11 ± 0,45″) | 43,11 | OK |
| F:ganzes | 43,11 ± 0,45″) | 0,45 | OK |
| F:ganzes | geht 4,46·10⁻¹⁰ vor | 4,4567 | OK |
| F:ganzes | 38,5 µs am Tag | 38,506 | OK |
| F:ganzes | | 4,465·10⁻¹⁰ (Werkseinstellung | 4,4647 | OK |
| F:ganzes | | 6605 mas/a | | 6604,8 | OK |
| F:ganzes | 6602 ± 18 mas/a | 6601,8 | OK |
| F:ganzes | 6602 ± 18 | 18,3 | OK |
| F:ganzes | | 40,8 mas/a | | 40,818 | OK |
| F:ganzes | 37,2 ± 7,2 mas/a | 37,2 | OK |
| F:ganzes | 37,2 ± 7,2 | 7,2 | OK |
| F:ganzes | Ṗ_b = −2,40·10⁻¹² | -2,4021 | OK |
| F:ganzes | 0,9983 ± 0,0016 davon | 0,9983 | OK |
| F:ganzes | 0,9983 ± 0,0016 | 0,0016 | OK |
| F:ganzes | | 39,7 µas | | 39,688 | OK |
| F:ganzes | 42 ± 3 µas | 42 | OK |
| F:ganzes | gleich c auf 10⁻¹⁵ | 1 | OK |
| 07 | ½, ¾ und ⅞ für d = 1, 2, 3 | 3 | OK |
| 07 | ∫n_F ds = ln 2 = 0,693 | 0,69315 | OK |
| 02 | alle 884 gerade-geraden Kerne haben 0⁺ | 884 | OK |
| 02 | sind 849 von 850 Spins halbzahlig | 849 | OK |
| 02 | 849 von 850 Spins | 850 | OK |
| 02 | bei den 352 ungerade-ungeraden alle ganzzahlig | 352 | OK |
| F:ganzes | (11/4)^{1/3} = 1,401 | 1,401 | OK |
| F:ganzes | T_ν = 1,945 K | 1,9454 | OK |
| 07 | ∏(1 + qⁿ), tief ln Z = π²/(12s) − ½·ln 2 + s/24 | 1 | OK |
| 02 | positiv, 10,7·A^{−1/2} MeV | 10,672 | OK |
| 02 | Stabil sind nach NUBASE2020 253 Kerne | 253 | OK |
| 02 | 253 Kerne, 147 gerade-gerade | 147 | OK |
| 02 | 102 mit ungeradem A | 102 | OK |
| 02 | nur 4 ungerade-ungerade | 4 | OK |
| 07 | (Euler, bei q = ½ beide 2,384) | 2,3842 | OK |
| 07 | bei q = ½ beide 2,384) | 2,3842 | OK |
| 01 | Abgabe der Windung quer zur Schwingung 2/3 | 0,66667 | OK |
| F:ganzes | denselben vollen Umlauf geben, ν = 2π | 1 | OK |
| F:ganzes | d/π = 0,955 | 0,95493 | OK |
| F:ganzes | √(d/π) = 0,977 | 0,97721 | OK |
| F:ganzes | δ = 2·(ln q)⁴·Θ²/Δ = 9,4·10⁻⁶ heute | 9,3959 | OK |
| F:ganzes | δ = 2·(ln q)⁴·Θ²/Δ = 9,4·10⁻⁶ heute | 9,3959 | OK |
| F:ganzes | Θ_∞/2π = √15/(2π²) = 0,196 | 0,19621 | OK |
| F:ganzes | Θ_∞/2π = √15/(2π²) = 0,196 | 0,19621 | OK |
| F:ganzes | für die Milchstraße 1759·0,196 = 345 | 1758,8 | OK |
| F:ganzes | für die Milchstraße 1759·0,196 = 345 | 345,1 | OK |
| F:ganzes | für die Milchstraße 1759·0,196 = 345 | 345,1 | OK |
| F:ganzes | Spinzahl weit über 2π/Θ_∞ = 5,1 | 5,0966 | OK |
| F:ganzes | Bei 1 + Q = 6·(9/8)² = 7,59 ist sie 4/9 | 7,5938 | OK |
| F:ganzes | Bei 1 + Q = 6·(9/8)² = 7,59 ist sie 4/9 | 7,5937 | OK |
| F:ganzes | Im Fenster 5 < Q < 6,59 dazwischen | 6,5938 | OK |
| F:ganzes | das gemessene Verhältnis, 1 + Q = 6,36 | 6,3643 | OK |
| F:ganzes | ganz bleibt nur der Knoten, beim Anteil 0,080 der Masse | 0,079676 | OK |
| F:galaxien | heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11 | 0,063843 | OK |
| F:galaxien | heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11 | 0,063843 | OK |
| F:galaxien | heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11 | 0,10917 | OK |
| F:galaxien | heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11 | 0,10917 | OK |
| F:galaxien | wächst R_d wie M^{5/12 − 0,11} = M^{0,30} | 0,3037 | OK |
| F:galaxien | wächst R_d wie M^{5/12 − 0,11} = M^{0,30} | -0,11297 | OK |
| F:galaxien | an SPARC wachsen die Scheiben wie M^{0,31 ± 0,02} | 0,31436 | OK |
| F:galaxien | an SPARC wachsen die Scheiben wie M^{0,31 ± 0,02} | 0,017354 | OK |
| F:ganzes | atmet die Sonne in 2π·√(R³/GM) = 2,8 h | 2,7801 | OK |
| F:ganzes | atmet die Sonne in 2π·√(R³/GM) = 2,8 h | 2,7801 | OK |
| F:ganzes | mit κ = 1,26·Ω in 0,20 Gyr | 1,2578 | OK |
| F:ganzes | mit κ = 1,26·Ω in 0,20 Gyr | 1,2578 | OK |
| F:ganzes | mit κ = 1,26·Ω in 0,20 Gyr | 0,20238 | OK |
| F:ganzes | mit κ = 1,26·Ω in 0,20 Gyr | 0,20238 | OK |
| F:ganzes | K = 100·H₀² − 4·H_∞², in 9,3 Gyr | 9,2699 | OK |
| F:ganzes | K = 100·H₀² − 4·H_∞², in 9,3 Gyr | 9,2699 | OK |
| F:ganzes | = −3·H_∞² für jeden Träger | -3 | OK |
| F:ganzes | = −3·H_∞² für jeden Träger | -3 | OK |
| F:ganzes | = −3·H_∞² für jeden Träger | -3 | OK |
| F:ganzes | (1/(√3·H_∞) = 10,2 Gyr) | 10,152 | OK |
| F:ganzes | (1/(√3·H_∞) = 10,2 Gyr) | 10,152 | OK |
| F:ganzes | H_∞ = 2π²a₀/(√15·c) = 55,6 km/s/Mpc | 55,606 | OK |
| F:ganzes | H_∞ = 2π²a₀/(√15·c) = 55,6 km/s/Mpc | 55,606 | OK |
| F:ganzes | 1/H_∞ = 17,6 Gyr | 17,584 | OK |
| F:ganzes | 1/H_∞ = 17,6 Gyr | 17,584 | OK |
| F:ganzes | Heute atmet der Formteil mit H₀ = 67,2 | 67,202 | OK |
| F:ganzes | Heute atmet der Formteil mit H₀ = 67,2 | 67,202 | OK |
| F:ganzes | τ_D = c/a₀ = 89,6 Gyr für jede Masse | 89,621 | OK |
| F:ganzes | τ_D = c/a₀ = 89,6 Gyr für jede Masse | 89,621 | OK |
| F:ganzes | H_∞·τ_D = 2π/Θ_∞ = 5,10 | 5,0966 | OK |
| F:ganzes | H_∞·τ_D = 2π/Θ_∞ = 5,10 | 5,0966 | OK |
| F:ganzes | Eine Rückkehr dauert 5,1 Zeiten des Atmens | 5,0966 | OK |
| F:ganzes | um den Faktor 163 | 163,47 | OK |
| F:ganzes | um den Faktor 163 | 163,47 | OK |
| F:ganzes | gleichwertig Θ_∞/2π = 0,196 | 0,19621 | OK |
| F:ganzes | gleichwertig Θ_∞/2π = 0,196 | 0,19621 | OK |
| F:ganzes | hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22 | 1,7639 | OK |
| F:ganzes | hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22 | 1,7795 | OK |
| F:ganzes | hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22 | 0,21948 | OK |
| F:ganzes | heute um 1,1 %, und H fällt dabei stets | 1,0604 | OK |
| F:ganzes | heute um 1,1 %, und H fällt dabei stets | 1 | OK |
| F:ganzes | ρ_T = 3H_∞²/(4πG) = 1,2·10⁻²⁶ kg/m³ | 1,1616 | OK |
| F:ganzes | ρ_T = 3H_∞²/(4πG) = 1,2·10⁻²⁶ kg/m³ | 1,1616 | OK |
| F:ganzes | sieben Wasserstoffatome je Kubikmeter | 6,9408 | OK |
| F:ganzes | sieben Wasserstoffatome je Kubikmeter | 6,9408 | OK |
| F:ganzes | die Sonne mit ρ̄/ρ_T = 1,2·10²⁹ | 1,2137 | OK |
| F:ganzes | die Sonne mit ρ̄/ρ_T = 1,2·10²⁹ | 1,2137 | OK |
| F:ganzes | ein Haufen in R₂₀₀ heute mit 146 | 146,06 | OK |
| F:ganzes | ein Haufen in R₂₀₀ heute mit 146 | 146,06 | OK |
| F:ganzes | im Formteil bei z = 0,63 | 0,63161 | OK |
| F:ganzes | im Formteil bei z = 0,63 | 0,63109 | OK |
| F:ganzes | kT = ħH/(2π), heute 2,6·10⁻³⁰ K | 2,6476 | OK |
| F:ganzes | kT = ħH/(2π), heute 2,6·10⁻³⁰ K | 2,6476 | OK |
| F:ganzes | im Endzustand 2,2·10⁻³⁰ K | 2,1907 | OK |
| F:ganzes | im Endzustand 2,2·10⁻³⁰ K | 2,1907 | OK |
| F:ganzes | 1/(30·S_Λ/k) = 3,3·10⁻¹²⁵ | 3,3383 | OK |
| F:ganzes | 1/(30·S_Λ/k) = 3,3·10⁻¹²⁵ | 3,3383 | OK |
| F:ganzes | kT_I = ħa₀/c = 2,7·10⁻³⁰ K | 2,7007 | OK |
| F:ganzes | kT_I = ħa₀/c = 2,7·10⁻³⁰ K | 2,7007 | OK |
| F:ganzes | Heute ist das Atmen 2 % kühler als die Rückkehr | 1,9682 | OK |
| F:ganzes | Heute ist das Atmen 2 % kühler als die Rückkehr | 1,9682 | OK |
| F:ganzes | im Endzustand um den Faktor Θ_∞ = 1,23 | 1,2328 | OK |
| F:ganzes | im Endzustand um den Faktor Θ_∞ = 1,23 | 1,2328 | OK |
| F:ganzes | Heute liegt das Mittel bei 0,23·ρ_T | 0,23022 | OK |
| F:ganzes | Heute liegt das Mittel bei 0,23·ρ_T | 0,23022 | OK |
| F:ganzes | erst ab dem 17-fachen des Mittels | 17,374 | OK |
| F:ganzes | erst ab dem 17-fachen des Mittels | 17,374 | OK |
| F:lesungen | im besten Punkt α = 3,2° | 3,2 | OK |
| F:lesungen | im besten Punkt α = 3,2° | 3,2281 | OK |
| F:lesungen | gibt das beste Paar S = 0,78 | 0,78195 | OK |
| F:lesungen | die beste Dreiheit senkrechter Paare S₃ = 0,58 | 0,58374 | OK |
| F:lesungen | gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63 | 0,76468 | OK |
| F:lesungen | gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63 | 0,81657 | OK |
| F:lesungen | gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63 | 0,57657 | OK |
| F:lesungen | gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63 | 0,63419 | OK |
| F:lesungen | Ein Kreis von 20° Umfang | 20,096 | OK |
| F:lesungen | Streuung σ = 0,23 je Paar nur 1/σ² = 18 unabhängige Flecken | 0,23399 | OK |
| F:lesungen | Streuung σ = 0,23 je Paar nur 1/σ² = 18 unabhängige Flecken | 18,265 | OK |
| F:lesungen | am ersten Gipfel (ℓ = 220) 11 % von C_ℓ, im Tal (ℓ = 400) 47 % | 10,732 | OK |
| F:lesungen | am ersten Gipfel (ℓ = 220) 11 % von C_ℓ, im Tal (ℓ = 400) 47 % | 46,875 | OK |
| F:lesungen | die Zelle erwartet S = 0,25 bis 0,37 | 0,24626 | OK |
| F:lesungen | die Zelle erwartet S = 0,25 bis 0,37 | 0,37274 | OK |
| F:lesungen | liest die Messung an den Achsen S = 0,25 ± 0,06 | 0,25219 | OK |
| F:lesungen | liest die Messung an den Achsen S = 0,25 ± 0,06 | 0,063578 | OK |
| F:lesungen | die Zerlegung dort 0,28 bis 0,37 | 0,28055 | OK |
| F:lesungen | die Zerlegung dort 0,28 bis 0,37 | 0,36575 | OK |
| F:ganzes | (CAMB, x_e = 2,9·10⁻⁴) | 2,8774 | OK |
| F:ganzes | trägt das Gas die Wärme des Lichts bis z_t = 123 | 123,38 | OK |
| F:ganzes | in der geschlossenen Form 122 | 122,49 | OK |
| F:ganzes | z = 17 auf 7,1 K gegen 49,1 K des Lichts | 7,0998 | OK |
| F:ganzes | z = 17 auf 7,1 K gegen 49,1 K des Lichts | 49,059 | OK |
| F:ganzes | z = 17 auf 7,1 K gegen 49,1 K des Lichts | 49,059 | OK |
| F:ganzes | mit dem Übergang y 7,3 K, in der vollen Rechnung von CAMB 6,9 K | 7,3286 | OK |
| F:ganzes | mit dem Übergang y 7,3 K, in der vollen Rechnung von CAMB 6,9 K | 6,8832 | OK |
| F:ganzes | A₁₀ = μ₀ω³μ_B²/(3πħc³) = 2,87·10⁻¹⁵ s⁻¹ | 2,8689 | OK |
| F:ganzes | A₁₀ = μ₀ω³μ_B²/(3πħc³) = 2,87·10⁻¹⁵ s⁻¹ | 2,8689 | OK |
| F:ganzes | bei z = 17 (79 MHz) τ = 0,09 | 78,911 | OK |
| F:ganzes | bei z = 17 (79 MHz) τ = 0,09 | 0,09123 | OK |
| F:ganzes | das Signal gegen das Licht −0,20 K (−196 bis −210 mK) | -0,20325 | OK |
| F:ganzes | das Signal gegen das Licht −0,20 K (−196 bis −210 mK) | -196,11 | OK |
| F:ganzes | das Signal gegen das Licht −0,20 K (−196 bis −210 mK) | -210,43 | OK |
| F:ganzes | EDGES liest bei 78 MHz eine Aufnahme von 0,5 K | 78 | OK |
| F:ganzes | EDGES liest bei 78 MHz eine Aufnahme von 0,5 K | 0,5 | OK |
| F:ganzes | das 2,4- bis 2,5-fache dessen, was das Modell dort zulässt | 2,376 | OK |
| F:ganzes | das 2,4- bis 2,5-fache dessen, was das Modell dort zulässt | 2,5496 | OK |
| F:ganzes | SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 % | 55 | OK |
| F:ganzes | SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 % | 85 | OK |
| F:ganzes | SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 % | 95,3 | OK |
| F:ganzes | bei x = 0,25000 | 0,25 | OK |
| F:ganzes | bei x = 0,25000 | 0,25 | OK |
| F:ganzes | senkt es um 2Z/(2Z + 1) auf 0,245 | 0,24535 | OK |
| F:ganzes | senkt es um 2Z/(2Z + 1) auf 0,245 | 0,24535 | OK |
| F:ganzes | Bei z = 10 fehlen 0,29 % | 0,29291 | OK |
| F:ganzes | Bei z = 10 fehlen 0,29 % | 0,29291 | OK |
| F:ganzes | ein Prozent erst bei z = 34 | 34,383 | OK |
| F:ganzes | ein Prozent erst bei z = 34 | 34,383 | OK |
| F:ganzes | an der letzten Streufläche 24 % | 24,253 | OK |
| F:ganzes | an der letzten Streufläche 24 % | 24,253 | OK |
| F:ganzes | dort trug der Vakuumteil das 0,57-fache von heute | 0,57375 | OK |
| F:ganzes | dort trug der Vakuumteil das 0,57-fache von heute | 0,57375 | OK |
| F:ganzes | Seit z ≈ 10 ist die Schwelle so auf 0,3 % fest | 0,29291 | OK |
| F:ganzes | noch um a_eq = 0,03 % | 0,029386 | OK |
| F:ganzes | noch um a_eq = 0,03 % | 0,029386 | OK |
| F:ganzes | erst bei v_flach = π²c/√15 = 2,55·c | 2,5483 | OK |
| F:ganzes | erst bei v_flach = π²c/√15 = 2,55·c | 2,5483 | OK |
| F:ganzes | trägt die mittlere Achse 19/30 von ⟨J²⟩ | 19 | OK |
| F:ganzes | trägt die mittlere Achse 19/30 von ⟨J²⟩ | 18,99 | OK |
| F:ganzes | jede äußere 11/60 | 11 | OK |
| F:ganzes | jede äußere 11/60 | 11,044 | OK |
| F:ganzes | jede äußere 11/60 | 10,976 | OK |
| F:ganzes | ⟨(λ₁ − λ₃)²⟩ = (19/15)·σ² | 19 | OK |
| F:ganzes | ⟨(λ₁ − λ₃)²⟩ = (19/15)·σ² | 19,003 | OK |
| F:ganzes | gegen je (11/30)·σ² | 11 | OK |
| F:ganzes | die Schwelle des Bereichs, δ_c = 1,686 | 1,6865 | OK |
| F:ganzes | die des Leerraums, δ_v = −2,72 | -2,7172 | OK |
| F:ganzes | mit δ_c/(δ_c + |δ_v|) = 0,38 in einem Leerraum | 0,38297 | OK |
| F:ganzes | mit δ_c/(δ_c + |δ_v|) = 0,38 in einem Leerraum | 0,38426 | OK |
| F:ganzes | und mit 0,62 in einem Bereich | 0,61703 | OK |
| F:ganzes | und mit 0,62 in einem Bereich | 0,61574 | OK |
| F:ganzes | Mit 1 + δ = 0,20 am Rand | 0,20472 | OK |
| F:ganzes | verlangten diese Leerräume das 1,9-fache des Raums | 1,8707 | OK |
| F:ganzes | verlangten diese Leerräume das 1,9-fache des Raums | 1,877 | OK |
| F:ganzes | Knoten und Leerraum tragen je 0,080 der Masse | 0,079676 | OK |
| F:ganzes | Knoten und Leerraum tragen je 0,080 der Masse | 0,079676 | OK |
| F:ganzes | Snider, 22,5 m) | 22,5 | OK |
| F:ganzes | | 2,46·10⁻¹⁵ | | 2,4585 | OK |
| F:ganzes | gemessen auf 1 % genau | 1 | OK |
| F:ganzes | Masse aus der Bahn 4,30·10⁶ M☉ | 4,296 | OK |
| F:ganzes | | 12,2′ je Umlauf | | 12,22 | OK |
| F:ganzes | 1,10 ± 0,19 davon | 1,1 | OK |
| F:ganzes | 1,10 ± 0,19 davon (GRAVITY) | 0,19 | OK |
| F:ganzes | | 53,4 µas | | 53,422 | OK |
| F:ganzes | 51,8 ± 2,3 µas | 51,8 | OK |
| F:ganzes | ± 2,3 µas (Event | 2,3 | OK |
| F:ganzes | Ṗ_b = −1,248·10⁻¹² | -1,2478 | OK |
| F:ganzes | | 0,75 ± 0,13 ± 0,16 g | | 0,75 | OK |
| F:ganzes | 0,75 ± 0,13 ± 0,16 g | 0,13 | OK |
| F:ganzes | ± 0,13 ± 0,16 g | | 0,16 | OK |
| F:ganzes | bestätigt auf 1,3·10⁻⁴ | 1,3 | OK |
| 03 | 615 Bins, τ | 615 | OK |
| 03 | τ = 0,054 ± 0,007 und Eichung | 0,054 | OK |
| 03 | 0,054 ± 0,007 und Eichung | 0,007 | OK |
| 03 | allein liegt das Modell Δχ² = 1,0 über | 0,9953 | OK |
| 03 | an allen Daten 0,01. | 0,014855 | OK |
| 03 | an allen Daten 0,01. | 0,014855 | OK |
| 03 | DESI mit χ² = 18,6 | 18,608 | OK |
| 03 | 18,6 (13 Werte) | 13 | OK |
| 03 | ΛCDM mit 31,4 | 31,438 | OK |
| 03 | gewinnt das Modell Δχ² = 9,6 | 9,6195 | OK |
| 03 | w₀ = −0,84 ± 0,05 | -0,83943 | OK |
| 03 | w₀ = −0,84 ± 0,05 | 0,052894 | OK |
| 03 | w_a = −0,56 ± 0,18 | -0,56368 | OK |
| 03 | w_a = −0,56 ± 0,18 | 0,18311 | OK |
| 03 | bei z = 0,40 ± 0,07 | 0,39831 | OK |
| 03 | bei z = 0,40 ± 0,07 | 0,073489 | OK |
| 03 | die Schwelle liest 1,070 ± 0,013 | 1,0704 | OK |
| 03 | die Schwelle liest 1,070 ± 0,013 | 0,013404 | OK |
| 03 | der Standardwert nur Δχ² = 0,8 daneben | 0,83463 | OK |
| F:lesungen | nur auf 1,3 % genau | 1,2522 | OK |
| 03 | die Dunkle Energie heute als a₀ = 1,085 ± 0,003 | 1,0854 | OK |
| 03 | die Dunkle Energie heute als a₀ = 1,085 ± 0,003 | 0,0025843 | OK |
| 03 | Das Maximum liegt dann bei z = 0,32 | 0,31808 | OK |
| 03 | wie ΛCDM und liegt Δχ² = 8,7 darunter | 8,7035 | OK |
| 03 | der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils | 5,4431 | OK |
| 03 | der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils | 0,66616 | OK |
| F:ganzes | F = 0,97 ± 0,12 % | 0,96933 | OK |
| F:ganzes | F = 0,97 ± 0,12 % | 0,11863 | OK |
| F:ganzes | am Maximum (z = 0,30) F = 1,68 ± 0,36 % | 1,6838 | OK |
| F:ganzes | am Maximum (z = 0,30) F = 1,68 ± 0,36 % | 0,36412 | OK |
| F:ganzes | am Maximum (z = 0,30) F = 1,68 | 0,29741 | OK |
| F:ganzes | bei z = 2 läge F bei −3,0 % | -2,9718 | OK |
| F:ganzes | lesen alle Präzisionsdaten χ² = 2023,9 | 2023,9 | OK |
| F:ganzes | 3,9 unter ΛCDM bei gleich vielen Zahlen, 4,2 über dem Gang | 3,8955 | OK |
| F:ganzes | 3,9 unter ΛCDM bei gleich vielen Zahlen, 4,2 über dem Gang | 4,2478 | OK |
| 10 | die Faltung das 8,1·10⁻³⁷-fache der Windung | 8,0936 | OK |
| 10 | numerisch auf 10⁻⁵ gleich | 1 | OK |
| 10 | auf das 1,582-fache | 1,582 | OK |
| 10 | Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache | 1 | OK |
| 10 | Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache | 9 | OK |
| 10 | Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache | 99 | OK |
| 10 | von 90 auf 99 % der Rate elffach | 11 | OK |
| 10 | Aus dem Gleichgewicht numerisch gelöst gleich auf 10⁻¹⁰ | 1 | OK |
| 10 | bei s = 0,1 um 0,508 und 0,492 Gaben | 0,50833 | OK |
| 10 | bei s = 0,1 um 0,508 und 0,492 Gaben | 0,49167 | OK |
| 10 | Gabe für Gabe summiert gleich auf 10⁻¹² | 1 | OK |
| 10 | in einer gewürfelten Folge von Binden und Lösen gezählt gleich auf 0,2 % | 1 | OK |
| 10 | W = S/(S + k_Z) = τ·[L]/(K_d + (1 + τ)·[L]) | 1 | OK |
| 10 | in einer gewürfelten Folge mit Grundtakt gezählt gleich auf 0,2 % | 1 | OK |
| 10 | Wirkstoffs mit Δω = +0,3·K auf das 0,39-fache | 0,39439 | OK |
| 10 | Δω = −0,3·K, auf das 2,79-fache | 2,7945 | OK |
| 10 | Ein genau passender Wirkstoff verliert immer, auf das 0,30-fache | 0,29773 | OK |
| 10 | aus dem abgetasteten Potential gleich auf 10⁻⁹ | 1 | OK |
| 10 | ΔV ≈ 2K − π·|Δω| | 1 | OK |
| 10 | bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11 | 2 | OK |
| 10 | bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11 | 5 | OK |
| 10 | bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11 | 11 | OK |
| 10 | aus der selbstkonsistenten Besetzung abgeleitet gleich auf 10⁻⁹ | 1 | OK |
| 10 | zwischen 0,99·T_K und 0,999·T_K mit dem Exponenten 1,505 | 1,5047 | OK |
| 10 | über die Form der eingerasteten Windungen gleich | 1 | OK |
| 10 | bei T = 2·T_N zwei Drittel | 0,66667 | OK |
| 10 | Gruppe für Gruppe im selben Zustand wie nach −7, −1, zu 96 % besetzt | 1 | OK |
| 10 | Gruppe für Gruppe im selben Zustand wie nach −7, −1, zu 96 % besetzt | 96,333 | OK |
| 10 | wieder im Zustand nach −7, −1, −5, zu 13 % besetzt | 1 | OK |
| 10 | wieder im Zustand nach −7, −1, −5, zu 13 % besetzt | 13 | OK |
| 10 | dazwischen bei −3 waren es 40 % | 39,667 | OK |
| 10 | sobald der Zustand an jedem Extremwert einschwingt | 1 | OK |
| 10 | aus finiten Differenzen gleich auf 10⁻⁸ | 1 | OK |
| 10 | bei x = 1 halb so groß | 0,5 | OK |
| 10 | bei x = 10 elfmal so lang | 11 | OK |
| 10 | im Kanal zehnfach, 10,5, 100,5, 1000,5 | 10,508 | OK |
| 10 | im Kanal zehnfach, 10,5, 100,5, 1000,5 | 100,5 | OK |
| 10 | im Kanal zehnfach, 10,5, 100,5, 1000,5 | 1000,5 | OK |
| 10 | auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1 | 3,599 | OK |
| 10 | auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1 | 6,8149 | OK |
| 10 | auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1 | 10,133 | OK |
| 10 | in freier Lösung kaum, 1,83 und 1,98 | 1,8262 | OK |
| 10 | in freier Lösung kaum, 1,83 und 1,98 | 1,9803 | OK |
| 10 | Aus finiten Differenzen gleich auf 10⁻⁴ | 1 | OK |
| 10 | bei s = 1, also E = 1 − 1/e = 0,63 | 0,63212 | OK |
| 10 | hält der Körper den Wirkstoff ν = 1,582 Durchgänge | 1,582 | OK |
| 10 | dieselbe Verstärkung wie die Regel an r* | 1 | OK |
| 10 | beide mit der Rate (1 − q)/τ, numerisch integriert gleich auf 10⁻¹¹ | 1 | OK |
| 10 | bei |q| = 1 je die Hälfte | 0,5 | OK |
| 10 | bei |q| = 1 je die Hälfte | 0,5 | OK |
| 10 | bei |q| = 3 bleibt ein Viertel | 0,25 | OK |
| 10 | Bei q = 0,9 wächst die Wirkung auf das Zehnfache | 10 | OK |
| 10 | läuft nach dem Absetzen mit dem Neunfachen nach | 9 | OK |
| 10 | klingt mit einem Zehntel der Rate ab | 0,1 | OK |
| 10 | für m = 1 genau K/π = 0,064 breit | 0,063662 | OK |
| 10 | für m = 2 noch 0,0032, für m = 3 noch 0,00028 | 0,0031726 | OK |
| 10 | für m = 2 noch 0,0032, für m = 3 noch 0,00028 | 2,7917·10⁻⁴ | OK |
| 10 | mit den Exponenten 1,000, 1,996 und 2,994 | 1 | OK |
| 10 | mit den Exponenten 1,000, 1,996 und 2,994 | 1,9964 | OK |
| 10 | mit den Exponenten 1,000, 1,996 und 2,994 | 2,994 | OK |
| 10 | Aus der Umlaufzahl in 10⁶ Schritten gleich auf 10⁻⁵ | 1 | OK |
| 10 | Bei E = 0,1 sind es 10 Durchgänge gegen 9,49 | 10 | OK |
| 10 | Bei E = 0,1 sind es 10 Durchgänge gegen 9,49 | 9,4912 | OK |
| 10 | bei E = 0,01 ist der Unterschied 0,501 | 0,50084 | OK |
| 10 | Gliedweise summiert gleich auf 10⁻¹² | 1 | OK |
| 10 | Gliedweise summiert sind beide gleich auf 10⁻¹⁵ | 1 | OK |
| 10 | bei n_H = 2 mit der halben Rate, numerisch integriert gleich auf 10⁻⁵ | 0,5 | OK |
| 10 | bei n_H = 2 mit der halben Rate, numerisch integriert gleich auf 10⁻⁵ | 1 | OK |
| 10 | bei J = 4,01 das 0,333-fache von (J − 4)^{3/2} | 0,33308 | OK |
| 10 | Durch Zählen der Lösungen gleich auf 10⁻⁸ | 1 | OK |
| 10 | Vor jeder Falte fällt die Rate der Rückkehr mit der Wurzel des Abstands | 1 | OK |
| 10 | und ihr Verhältnis K_d bleibt | 1 | OK |
| 10 | geben stationäre Verteilung und mittlere Übergangszeiten dasselbe, numerisch gleich auf 10⁻¹² | 1 | OK |
| 10 | Bei c_eff = 10·K_d hält der Wirkstoff 11-mal so lange | 11 | OK |
| 10 | (bei s = 0,2 das 1,90-fache) | 1,9048 | OK |
| 10 | J = 2·kT gibt n_H = 2 | 2 | OK |
| 10 | numerisch abgeleitet gleich auf 10⁻⁶ | 1 | OK |
| 10 | bei W = 10·k sind es 7 | 7 | OK |
| F:andocken | An 68 Kinasehemmern gegen 442 Kinasen | 68 | OK |
| F:andocken | An 68 Kinasehemmern gegen 442 Kinasen | 442 | OK |
| F:andocken | haben die Rangkorrelation 0,43, mit 1 µM als Schwelle 0,41 | 0,4316 | OK |
| F:andocken | haben die Rangkorrelation 0,43, mit 1 µM als Schwelle 0,41 | 0,41225 | OK |
| F:andocken | ihr oberes Zehntel beginnt bei 0,69 der stärksten statt bei 0,85 | 0,69075 | OK |
| F:andocken | im ähnlichsten Viertel 38 %, im fernsten 24 % | 37,645 | OK |
| F:andocken | im ähnlichsten Viertel 38 %, im fernsten 24 % | 23,69 | OK |
| F:andocken | mit Gruppen aus zwei Aminosäuren 35 % und 25 % | 35,109 | OK |
| F:andocken | mit Gruppen aus zwei Aminosäuren 35 % und 25 % | 24,911 | OK |
| F:andocken | ihr oberes Zehntel beginnt bei 0,69 der stärksten statt bei 0,85 | 0,84792 | OK |
| 10 | θ(h) + θ(−h) = 1 | 1 | OK |
| 10 | Bei Δω = 1,01·K dauert ein Schlupf 44/K, bei Δω = 2·K nur 3,6/K | 44,318 | OK |
| 10 | Bei Δω = 1,01·K dauert ein Schlupf 44/K, bei Δω = 2·K nur 3,6/K | 3,6276 | OK |
| 10 | bei voller Passung 2K, bei Δω = K/2 nur noch 34 % davon | 2 | OK |
| 10 | bei voller Passung 2K, bei Δω = K/2 nur noch 34 % davon | 34,243 | OK |
| 10 | Bei [L] = K_d braucht er das Doppelte seiner eigenen Bindungskonstante | 2 | OK |
| 10 | bei [L] = 3·K_d das Vierfache | 4 | OK |
| 10 | Aus der stationären Verteilung der Übergänge numerisch gleich auf 10⁻¹² | 1 | OK |
| 10 | drei Punkte legt immer eine Drehung auf ihr Spiegelbild, es bleibt kein Abstand | 1 | OK |
| 10 | ein mittlerer Abstand von 0,71 der Kantenlänge, über alle Drehungen abgetastet ebenso | 0,70711 | OK |
| 10 | ein mittlerer Abstand von 0,71 der Kantenlänge, über alle Drehungen abgetastet ebenso | 0,70711 | OK |
| 10 | Geschlossene Formen und numerische Integration stimmen auf 10⁻⁵ | 1 | OK |
| F:andocken | ist am größten bei N_e = (3c/a)^{3/2} | 1 | OK |
| F:andocken | gleich bei 68 von 68 | 68 | OK |
| F:andocken | die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37 | 0,33307 | OK |
| F:andocken | die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37 | 0,0055136 | OK |
| F:andocken | die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37 | 0,36965 | OK |
| F:andocken | Größe und stärkste Bindung 0,25 | 0,25072 | OK |
| F:andocken | wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8} | 1,9124 | OK |
| F:andocken | wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8} | 0,68481 | OK |
| F:andocken | wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8} | 2,4102 | OK |
| F:andocken | wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8} | 0,76841 | OK |
| F:andocken | An den K_i und K_d von ChEMBL (5985 Spiegelpaare | 5985 | OK |
| F:andocken | über vertauschtes @ und @@ in den SMILES zu 95 % bestätigt, 10 246 Paare | 95,455 | OK |
| F:andocken | über vertauschtes @ und @@ in den SMILES zu 95 % bestätigt, 10 246 Paare | 10246 | OK |
| F:andocken | im Median um 0,51 in pK | 0,51 | OK |
| F:andocken | die Rangkorrelation mit der Zahl der Atome ist −0,07 | -0,065065 | OK |
| F:andocken | bis 15 Atome 14 % der Paare, bei 31 bis 40 Atomen 22 % | 13,8 | OK |
| F:andocken | bis 15 Atome 14 % der Paare, bei 31 bis 40 Atomen 22 % | 21,569 | OK |
| F:andocken | An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare | 643 | OK |
| F:andocken | An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare | 30 | OK |
| F:andocken | An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare | 15 | OK |
| F:andocken | alle Antagonisten | 1 | OK |
| F:andocken | führt nur 22 k_off mit Ziel und Wirkstoff, keines trifft eine Wirkweise | 22 | OK |
| F:andocken | führt nur 22 k_off mit Ziel und Wirkstoff, keines trifft eine Wirkweise | 0 | OK |
| F:andocken | 235 599 Moleküle | 2,356·10⁵ | OK |
| F:andocken | an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1 | 2436 | OK |
| F:andocken | an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1 | 95,525 | OK |
| F:andocken | an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1 | 99,918 | OK |
| F:andocken | zwischen 12 und 13 Atomen von 9,5 auf 10,4 | 9,47 | OK |
| F:andocken | zwischen 12 und 13 Atomen von 9,5 auf 10,4 | 10,369 | OK |
| F:andocken | bis 66 Atome zwischen 9,9 und 10,9 | 66 | OK |
| F:andocken | bis 66 Atome zwischen 9,9 und 10,9 | 9,92 | OK |
| F:andocken | bis 66 Atome zwischen 9,9 und 10,9 | 10,92 | OK |
| F:ganzes | Dann ist F heute 0,55 ± 0,27 % | 0,55314 | OK |
| F:ganzes | Dann ist F heute 0,55 ± 0,27 % | 0,26779 | OK |
| F:ganzes | schmale Spitze bei z = 0,61, χ² = 2017,9 | 0,61322 | OK |
| F:ganzes | schmale Spitze bei z = 0,61, χ² = 2017,9 | 2017,9 | OK |
| F:ganzes | (F_s = 5 ± 5 %) | 5,1988 | OK |
| F:ganzes | (F_s = 5 ± 5 %) | 5,2529 | OK |
| F:ganzes | am Maximum (z = 0,30) 1,21 ± 0,33 % | 1,2089 | OK |
| F:ganzes | am Maximum (z = 0,30) 1,21 ± 0,33 % | 0,3287 | OK |
| F:ganzes | am Maximum (z = 0,30) 1,21 ± 0,33 % | 0,29858 | OK |
| 03 | lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006 | 5,3583 | OK |
| 03 | lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006 | 0,048776 | OK |
| 03 | lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006 | 1,076 | OK |
| 03 | lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006 | 0,0064264 | OK |
| 03 | 2π − 1 liegt Δχ² = 2,4 daneben. | 2,4206 | OK |
| 03 | ist L das 0,998-fache des Durchmessers | 0,9984 | OK |
| 03 | den Quadrupol der Hintergrundstrahlung um 13 % | 12,692 | OK |
| 03 | hebt den Oktupol um 2 % | 2,4427 | OK |
| 03 | Sachs-Wolfe: 16 % und 2 % | 15,761 | OK |
| 03 | Sachs-Wolfe: 16 % und 2 % | 2,3101 | OK |
| 03 | ändert das χ² nur um 0,08 | 0,083492 | OK |
| 03 | gegen 1434 mit 282 | 1434,4 | OK |
| 03 | 1434 mit 282 freien Zahlen | 282 | OK |
| 03 | χ² = 2918 gegen | 2917,6 | OK |
| 03 | liegt die Regel 736 vorn | 735,67 | OK |
| 03 | vorn, 2926 gegen | 2925,5 | OK |
| 03 | 2926 gegen 3661 | 3661,2 | OK |
| 03 | höchstens 0,02 dex | 0,018186 | OK |
| 03 | An SPARC liegen 131 Punkte jenseits der Naht | 131 | OK |
| 03 | ist +0,013 ± 0,006 dex | 0,012707 | OK |
| 03 | +0,013 ± 0,006 dex | 0,005793 | OK |
| 03 | lägen sie 0,094 dex darüber | 0,093703 | OK |
| 03 | bei 31 der 141 Galaxien am Rand | 31 | OK |
| 03 | um 0,58 dex: Der Halo | 0,57936 | OK |
| 02 | (κ = 0,12, μ = 0,03: Schließzahlen | 0,12 | OK |
| 02 | μ = 0,03: Schließzahlen | 0,03 | OK |
| 02 | Schließzahlen 2, 8, 20, 28, 50, 82), gegen | 1 | OK |
| 02 | A^{−1/3} = 38,0 MeV | 38,049 | OK |
| 02 | liest unabhängig davon 39,3 | 39,341 | OK |
| 02 | 39,3 (3 %) | 3,3977 | OK |
| 02 | also λ_H = 1,23 fm | 1,2286 | OK |
| 02 | streut die Bindung um 1,83 MeV | 1,8275 | OK |
| 02 | um 1,26 statt 1,93 MeV | 1,2617 | OK |
| 02 | 1,26 statt 1,93 MeV | 1,9274 | OK |
| F:ganzes | (82 gerade-gerade Strahler | 82 | OK |
| F:ganzes | streut um 0,52 Dekaden | 0,51739 | OK |
| F:ganzes | drei Zahlen (0,51) | 0,51272 | OK |
| F:ganzes | liegt fest 1,24 Dekaden zu schnell | 1,2412 | OK |
| F:ganzes | (N = 118 bis 126, 13 Strahler) | 13 | OK |
| F:ganzes | 1,69 ± 0,07 Dekaden zu schnell | 1,6856 | OK |
| F:ganzes | 1,69 ± 0,07 Dekaden | 0,072575 | OK |
| F:ganzes | (N = 128 bis 136, 23 Strahler) | 23 | OK |
| F:ganzes | nur 0,94 ± 0,03 | 0,94489 | OK |
| F:ganzes | 0,94 ± 0,03. | 0,026246 | OK |
| F:ganzes | schließt sich 5,5-mal leichter | 5,5038 | OK |
| F:ganzes | in einem von 17 Anläufen | 17,427 | OK |
| 08 | 0,8756″ | 0,8756 | OK |
| 08 | 42,98″ | 42,981 | OK |
| 08 | r* = 7480 AE | 7479,6 | OK |
| F:galaxien | An RC100 (100 sternbildende Galaxien | 100 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 1,0035 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 0,93533 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 1,0755 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 0,86685 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 1,2099 | OK |
| F:galaxien | liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01 | 1,0129 | OK |
| F:galaxien | Das g_N der Tabelle ist das Feld der Baryonen aus Masse und Radius: V_bar²·R_e/(GM_b) liegt im Median bei 0,60, die dünne Scheibe gibt bei R_e 0,63 | 0,59539 | OK |
| F:galaxien | Das g_N der Tabelle ist das Feld der Baryonen aus Masse und Radius: V_bar²·R_e/(GM_b) liegt im Median bei 0,60, die dünne Scheibe gibt bei R_e 0,63 | 0,62515 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 1,2759 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 1,1705 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 1,3904 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 109,71 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 100 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 1,3852 | OK |
| F:galaxien | Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06 | 59 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | 0,14102 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | -0,077966 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | 0,35882 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | -0,29287 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | -0,57568 | OK |
| F:galaxien | Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01) | -0,010529 | OK |
| F:galaxien | Wüchse a₀ mit H(z), läge χ² um 15 und 20 höher | 15,238 | OK |
| F:galaxien | Wüchse a₀ mit H(z), läge χ² um 15 und 20 höher | 19,98 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 26 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 0,95697 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 0,81338 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 1,1147 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 1,1233 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 1,0087 | OK |
| F:galaxien | In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23) | 1,2334 | OK |
| F:galaxien | liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02 | -0,029196 | OK |
| F:galaxien | liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02 | 0,050071 | OK |
| F:galaxien | liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02 | -0,039479 | OK |
| F:galaxien | liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02 | 0,017308 | OK |
| 07 | über die Verteilung integriert gleich auf 10⁻¹⁵ | 1 | OK |
| 08 | r*(m_Pl) = √(ℓ_P·ℓ₀) = 0,117 mm | 0,11706 | OK |
| 08 | r*(m_Pl) = √(ℓ_P·ℓ₀) = 0,117 mm | 1 | OK |
| 08 | ρ_s·c² = (π³/30)·ħc/r*(m_Pl)⁴ | 1 | OK |
| 08 | mit 1 + Q = 2π also ε_L·N_b = (3/4π)·m_Pl/m_p | 1 | OK |
| 08 | Die Zelle trägt N_b = 1,6·10⁸⁰ Baryonen, aus ρ_b·L³ gezählt ebenso | 1,6296 | OK |
| 08 | Die Zelle trägt N_b = 1,6·10⁸⁰ Baryonen, aus ρ_b·L³ gezählt ebenso | 1,6088 | OK |
| 03 | Θ ist das Verhältnis zweier Schleifen, der Wärmelänge des Horizonts 2πc/H₀ zur Schleife der inneren Richtung ℓ₀ = 2Λ | 1 | OK |
| 03 | Gleich lang sind sie, wenn drei Vakuumteile Ω_DE = π²/15 = 0,658 tragen | 0,65797 | OK |
| 03 | Gleich lang sind sie, wenn drei Vakuumteile Ω_DE = π²/15 = 0,658 tragen | 0,65797 | OK |
| 01 | Beide hängen über das Ganze zusammen: δ·(R/r*)² = R/ℓ₀ mit δ = 3Δs/R | 1 | OK |
| 01 | Bei der Verwindung ist es die Spinzahl, j = k·c·v/v_K² mit k = I/(MR²) | 1 | OK |
| 08 | M_c teilt die Leiter im Drittel: Dort ist λ̄ = r* = ℓ_P^{2/3}·ℓ₀^{1/3} | 1 | OK |
| 08 | M_c teilt die Leiter im Drittel: Dort ist λ̄ = r* = ℓ_P^{2/3}·ℓ₀^{1/3} | 1 | OK |
| F:ganzes | Die Schließung muss das Sechseck umfassen, 1 + Q ≥ 2d, und 1 + Q = 2π ist der Kreis durch seine Ecken, r₀/Λ = d/π | 1 | OK |
| F:ganzes | Die Schließung muss das Sechseck umfassen, 1 + Q ≥ 2d, und 1 + Q = 2π ist der Kreis durch seine Ecken, r₀/Λ = d/π | 1 | OK |
| F:ganzes | Die Schwelle des Ganzen, r* = √(r₀·Λ), ist der Spiegel zwischen Sechseck und Kreis | 1 | OK |
| 02 | An ihrer Schwelle ist δ = 3Δs/r* = r*/ℓ₀ = (v_flach/c)² | 1 | OK |
| 02 | an derselben Stelle, r*/Λ = √(6/(1 + Q)) = 0,977, und die Menge fällt heraus | 0,97721 | OK |
| 02 | an derselben Stelle, r*/Λ = √(6/(1 + Q)) = 0,977, und die Menge fällt heraus | 1 | OK |
| 03 | In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt | 4 | OK |
| 03 | In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt | 3,7525 | OK |
| 03 | In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt | 13,755 | OK |
| 03 | In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2 | 16,3 | OK |
| 03 | In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2 | 4,3 | OK |
| 03 | In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2 | 10,566 | OK |
| 03 | In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2 | 2,2189 | OK |
| 03 | Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05 | 0,05 | OK |
| 03 | Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05 | 0,02 | OK |
| 03 | Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05 | 0,36 | OK |
| 03 | Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05 | 0,2 | OK |
| 03 | Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7 | 7 | OK |
| 03 | Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7 | 0,93897 | OK |
| 03 | Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7 | 1 | OK |
| 03 | Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7 | 6,7315 | OK |
| 06 | gilt Σ1/⟨T⟩ = 1 | 1 | OK |
| 07 | in drei kehrte die Irrfahrt nur zu 0,341 zurück (Pólya) | 0,34054 | OK |
| 07 | in drei kehrte die Irrfahrt nur zu 0,341 zurück (Pólya) | 1 | OK |
| 07 | die Summe bliebe bei 1/(1 − 0,341) = 1,516 | 1,5164 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 141 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 2,1483 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 1,557 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 2,6459 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 0,7513 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 0,64339 | OK |
| 03 | An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05) | 1,0525 | OK |
| 03 | Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87 | 5,6569 | OK |
| 03 | Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87 | 0,5651 | OK |
| 03 | Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87 | 0,57795 | OK |
| 03 | Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87 | 0,86709 | OK |
| 03 | Die Standard-Form läge mit ihrem Gipfel bei 0,49·a₀ | 0,48587 | OK |
| F:ganzes | Im Katalog α.100 von ALFALFA (Code 1, bis cz = 15 000 km/s, 23 956 Quellen über der Grenze, an der die Hälfte erfasst ist) | 23956 | OK |
| F:ganzes | Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 % | 1,2106 | OK |
| F:ganzes | Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 % | 1,1217 | OK |
| F:ganzes | Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 % | 1,4659 | OK |
| F:ganzes | Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 % | 0,2449 | OK |
| F:ganzes | (gezählt in Zellen von 2° oder 1°: 6801 oder 6095 Quadratgrad) | 6800,9 | OK |
| F:ganzes | (gezählt in Zellen von 2° oder 1°: 6801 oder 6095 Quadratgrad) | 6094,8 | OK |
| F:ganzes | Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 % | 1171 | OK |
| F:ganzes | Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 % | 1,1718 | OK |
| F:ganzes | Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 % | 1,1718 | OK |
| F:ganzes | Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 % | 1,2201 | OK |
| AB | abhandlung_de.md: alle 395 Zahlen stehen im Modell | 395 | OK |
| AB | abhandlung_en.md: alle 395 Zahlen stehen im Modell | 395 | OK |
| RB | rechnen_de.md: alle 158 gerechneten Zahlen stehen | 158 | OK |
| RB | rechnen_en.md: alle 158 gerechneten Zahlen stehen | 158 | OK |
| EL | elemente_de.md entspricht der Rechnung | 1 | OK |
| EL | elemente_en.md entspricht der Rechnung | 1 | OK |

**Ergebnis:** alles OK · Bahndaten in `berichte/bericht_0609/` · Laufzeit 4 s
