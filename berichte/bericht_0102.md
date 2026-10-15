# Bericht 102

2026-10-02T13:06:15 · Commit 2514d47 · a₀ = 1,06·10⁻¹⁰ m/s², H₀ = 67,4 km/s/Mpc

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

SPARC, 141 Galaxien, 2688 Punkte, Abstand, Neigung und Masse-zu-Licht je Galaxie mitgeschätzt, Masse-zu-Licht der Scheiben so, dass Stern und Gas derselben Regel folgen (Υ_d = 0,67): a₀ = 1,11·10⁻¹⁰ m/s² (16–84 %: 1,00–1,21); χ² Bose 2917,6, einfach 2970,7, Standard 2958,4; Dimension der Rückkehr 0,92 ± 0,05; Rest gegen Gasanteil innerhalb der Galaxien 0,00 ± 0,04 dex; je Galaxie gegen Masse ρ = −0,06, gegen Gasanteil ρ = 0,07, gegen Hubble-Typ ρ = −0,03, gegen Flächenhelligkeit ρ = 0,03. AME2020: kleinste Masse je Nukleon A = 56, größte Bindung A = 62, Spaltbarkeit am Maximum der Tröpfchenformel 0,247. Die Daten tragen log g_bar, log g_obs je Punkt.

| Aussage | Modell | Simulation | Einheit | Urteil |
|---|---|---|---|---|
| a₀ aus SPARC gegen a₀ aus dem Atmen (Toleranz 1,5 σ) | 1,0646 | 1,1063 | 10⁻¹⁰ m/s² | OK |
| Dimension der Rückkehr (Toleranz 2 σ) | 1 | 0,91682 |  | OK |
| Bose-Form vor Standard-Form | 1 | 1 |  | OK |
| kleinste Masse je Nukleon bei A | 56 | 56 |  | OK |
| größte Bindung je Nukleon bei A | 62 | 62 |  | OK |
| Spaltbarkeit am Maximum der Bindung | 0,25 | 0,24697 |  | OK |
| Bindung des Wasserstoffs ist die halbe Spreizung α²/2 (PubChem) | 2,6611·10⁻⁵ | 2,6611·10⁻⁵ |  | OK |
| Rest je Galaxie unabhängig von Masse (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Gasanteil (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Typ (p > 0,05) | 1 | 1 |  | OK |
| Rest je Galaxie unabhängig von Flächenhelligkeit (p > 0,05) | 1 | 1 |  | OK |

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
| 03 | r* = 8,88 kpc | 8,8824 | OK |
| 03 | 170,4 km/s | 170,45 | OK |
| 03 | Λ = 13,7 Gpc | 13,739 | OK |
| 03 | 89,6 Mrd. Jahre | 89,621 | OK |
| 03 | L = 27,9 Gpc | 27,934 | OK |
| 03 | = 5,36 | 5,3643 | OK |
| 03 | q = 0,843 | 0,84287 | OK |
| 03 | Θ = 2πa₀/(cH₀) = 1,017 | 1,0171 | OK |
| 03 | Ω_s = 0,227 | 0,22688 | OK |
| 03 | H₀·√Ω_s = 32,1 km/s/Mpc | 32,104 | OK |
| 03 | k_s = π³/30 = 1,034 | 1,0335 | OK |
| 03 | Ω_DE/Ω_s = 3,03 | 3,0263 | OK |
| 03 | a₀ = 1,06·10⁻¹⁰ | 1,0646 | OK |
| 02 | Grenzmasse 1,456 M☉ | 1,4563 | OK |
| 03 | ν = 6,14 frei | 6,1439 | OK |
| 03 | 5,31 senkrecht | 5,3097 | OK |
| 03 | 4,55 gleichgerichtet | 4,5511 | OK |
| 03 | R_Q = 52,0 kpc | 51,963 | OK |
| 03 | a_N = 0,029·a₀ | 0,02922 | OK |
| 03 | χ = Q·(−ln q) = 0,917 | 0,91697 | OK |
| 03 | r*/Λ = 0,971 | 0,97096 | OK |
| 06 | ν = 1,1567 | 1,1567 | OK |
| 06 | s ≈ 3·10⁵ | 3,0422 | OK |
| 08 | Sonne 4,43 und 8,86 km | 4,4299 | OK |
| 08 | 8,86 km | 8,8598 | OK |
| 02 | Über 93 Elemente | 93 | OK |
| 02 | von 0,81 (As | 0,81376 | OK |
| 02 | bis 8,2 (Ga) | 8,1773 | OK |
| 02 | Ionisierungsenergie −0,48 | -0,48024 | OK |
| 03 | im Median 0,71·r* | 0,70523 | OK |
| 03 | (16–84 %: 0,37–1,28) | 0,36537 | OK |
| 03 | 0,37–1,28) über vier | 1,2765 | OK |
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
| 03 | = 0,20 mag heller | 0,19583 | OK |
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
| 00 | 1,3 σ neben der Messung | 1,2563 | OK |
| 03 | Pantheon+, 1590 Supernovae | 1590 | OK |
| 03 | nur Δχ² = 2,0 über | 1,9631 | OK |
| 03 | ΛCDM (Ω_m = 0,304) | 0,30393 | OK |
| 03 | Ω_m = 0,304), 1,4 σ | 1,4011 | OK |
| 03 | (w₀ = −0,89, | -0,89032 | OK |
| 03 | w_a = −0,21) | -0,21286 | OK |
| 03 | Δχ² = 6,9 unter | 6,871 | OK |
| 03 | dem Formteil, 1,8 σ | 1,7736 | OK |
| 03 | r_d = 149,0 Mpc | 149,05 | OK |
| 03 | allein (ω_b, ω_c nach Planck), ist a₀ = 1,09 | 1,0885 | OK |
| 03 | (1,07–1,11 bei Δχ² = 1) | 1,0681 | OK |
| 03 | 1,07–1,11 bei Δχ² | 1,1089 | OK |
| 03 | mit H₀ = 68,4: Rotationskurven | 68,442 | OK |
| 03 | lesen die Abstände a₀ = 1,08 | 1,0815 | OK |
| 03 | 1,08 (1,06–1,10) | 1,0613 | OK |
| 03 | (1,06–1,10), H₀ | 1,1019 | OK |
| 03 | H₀ = 68,0, bei gleichem | 68,005 | OK |
| 02 | nach Z = 2, 10, 18, 36, 54, 86 | 1 | OK |
| 02 | bei N = 16, 20, 28, 50, 82, 126 | 1 | OK |
| 02 | Z = 8, 16, 20, 28, 50, 82 (AME | 1 | OK |
| 02 | Faktor 4,6 bei He | 4,5599 | OK |
| 02 | Halterlänge λ_H = 1,25 fm | 1,2517 | OK |
| 02 | v/c ≈ λ̄_p/λ_H = 0,17 | 0,16802 | OK |
| 02 | Formel zu λ_H = 1,32 fm | 1,3197 | OK |
| 02 | (48 Paare, AME2020) | 48 | OK |
| 02 | unabhängig davon 0,14 | 0,13741 | OK |
| 02 | (Z = 118, A = 294 | 1 | OK |
| 02 | A = 294: x = 1,00 | 1,0036 | OK |
| 02 | erst bei Zα = 0,86 | 0,86109 | OK |
| 02 | (2n² je Schale: 2, 10, 28, 60) | 1 | OK |
| 02 | nach 2, 8, 20, 40, 70, 112 | 1 | OK |
| 02 | also nach 28, 50, 82, 126 (k = 2 bis 5) | 1 | OK |
| 02 | die nächste Schließung wäre 184 | 184 | OK |
| 02 | schließen nach 2, 10, 18, 36, 54, 86 (PubChem) | 1 | OK |
| 02 | 2,8 bei Rn | 2,7551 | OK |
| 03 | 40 km/s 72 % der | 72,082 | OK |
| 03 | 150 km/s 12 % | 12,363 | OK |
| 03 | Schwelle zu a₀ = 1,11·10⁻¹⁰ | 1,1063 | OK |
| 03 | (16–84 %: 1,00–1,21) | 1,0023 | OK |
| 03 | 1,00–1,21). Die Bose | 1,2125 | OK |
| 03 | die einfache (Δχ² = 53) | 53,128 | OK |
| 03 | Standard-Form (Δχ² = 41) | 40,747 | OK |
| 03 | Rückkehr ist 0,92 ± 0,05 | 0,91682 | OK |
| 03 | 0,92 ± 0,05; der Rest | 0,0526 | OK |
| 07 | Dimension 0,92 ± 0,05 | 0,91682 | OK |
| 02 | Tröpfchenformel: x = 0,25 | 0,24697 | OK |
| 03 | nach Planck a₀ = 1,065·10⁻¹⁰ | 1,0646 | OK |
| 03 | m/s² H₀ = 67,2 km/s/Mpc | 67,202 | OK |
| 03 | und Ω_DE = 0,685, mit | 0,68466 | OK |
| 03 | 1,1·10⁻¹⁰ H₀ = 68,9 | 68,948 | OK |
| 03 | 68,9 und 0,700 | 0,70043 | OK |
| 03 | erst bei z = −0,10 ein | -0,10209 | OK |
| 03 | = 0,47·v_flach | 0,4714 | OK |
| 03 | unter v_flach = 17–21 | 16,971 | OK |
| 03 | 17–21 km/s. In SPARC | 21,213 | OK |
| 03 | bei z_eq = 3402 | 3402 | OK |
| 03 | bei 3–4·10¹² | 3,1173 | OK |
| 03 | 3–4·10¹² M☉ (mit | 3,8339 | OK |
| 03 | Bereiche ab 1,4·10¹¹ | 1,4395 | OK |
| 03 | bis 1,2·10¹³ M☉ | 1,166 | OK |
| 03 | M_min = 0,6–1,4·10⁷ M☉ | 0,58961 | OK |
| 03 | 1,4·10⁷ M☉; kaltes | 1,4395 | OK |
| 03 | zum Endwert 1,23 | 1,2328 | OK |
| 03 | drei Teile tun es bei z = 0,30 | 0,29501 | OK |
| 03 | mit drei Teilen in 55,6 km/s/Mpc | 55,606 | OK |
| 03 | 9,1 km/s bei 10⁷ M☉ | 9,1295 | OK |
| 03 | 19,3 km/s bei 2·10⁸ M☉ | 19,307 | OK |
| 03 | 10¹² M☉/Mpc: 2,6 kpc | 2,6299 | OK |
| 03 | (ν = 1,007) | 1,0068 | OK |
| 03 | (ν = 3,3) | 3,3035 | OK |
| 03 | a₀/(πG) = 242 M☉/pc² | 242,07 | OK |
| 03 | a_N = 2,54·a₀ | 2,5396 | OK |
| 03 | beträgt 0,648·a₀ | 0,64761 | OK |
| 07 | Rosetten mit 254,6° | 254,56 | OK |
| 07 | e^{−7480} | 7479,6 | OK |
| 08 | bei der Milchstraße 0,42″ | 0,41894 | OK |
| 08 | ε_L = ℓ_P/ℓ₀ = 2,0·10⁻⁶² | 1,9062 | OK |
| 08 | ε_L^{1/3}·m_Pl = 5,8·10⁻²⁹ kg | 5,8139 | OK |
| 08 | c³τ_D/G = 1,1·10⁵⁴ kg | 1,1418 | OK |
| 08 | Δs(M_U) = ⅔·Λ | 2 | OK |
| 08 | (4π⁴/15)·ε_L² = 1,0·10⁻¹²² | 0,94388 | OK |
| 02 | r* = 32 fm | 32,453 | OK |
| 02 | a_B = 1631·r* | 1630,6 | OK |
| 04 | Sonne 1,05·10⁷⁷ k | 1,0489 | OK |
| 01 | 0,022 s | 0,021968 | OK |
| 01 | Sonnenoberfläche 67 s | 66,981 | OK |
| 01 | auf der Erdbahn sind es 0,31 s | 0,31149 | OK |
| 01 | Takt 0,81 | 0,8096 | OK |
| 08 | 1,7512″ | 1,7512 | OK |
| 08 | 0,8756″ | 0,8756 | OK |
| 08 | 42,98″ | 42,981 | OK |
| 08 | r* = 7480 AE | 7479,6 | OK |
| AB | abhandlung_de.md: alle 101 Zahlen stehen im Modell | 101 | OK |
| AB | abhandlung_en.md: alle 101 Zahlen stehen im Modell | 101 | OK |

**Ergebnis:** alles OK · Bahndaten in `berichte/bericht_0102/` · Laufzeit 71 s
