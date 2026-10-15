# Bericht 252

2026-10-13T23:26:13 · Commit bb69471 · a₀ = 1,06·10⁻¹⁰ m/s², H₀ = 67,4 km/s/Mpc

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
| 03 | (KiDS-1000, 60 Werte bis 2,6 Mpc) | 60 | OK |
| 03 | (KiDS-1000, 60 Werte bis 2,6 Mpc) | 2,6035 | OK |
| 03 | ohne äußeres Feld, χ² = 503 | 503,4 | OK |
| 03 | 10⁻⁴·a₀ im Betrag kostet Δχ² = 783 | 783,39 | OK |
| 03 | bei 35 kpc 0,25 dex darunter | 0,25311 | OK |
| 03 | bei 35 kpc 0,25 dex darunter | 35,39 | OK |
| 03 | bei 561 kpc 0,15 darüber | 0,14542 | OK |
| 03 | bei 561 kpc 0,15 darüber | 560,9 | OK |
| 03 | bei 2,6 Mpc 0,11 darunter | 0,10814 | OK |
| 03 | (15 Bins in g_bar, Brouwer 2021): ohne Feld χ² = 21 | 20,615 | OK |
| 03 | die blauen Linsen +0,01 ± 0,02 dex auf der Regel | 0,0071084 | OK |
| 03 | Sechs gasreiche (Mancera Piña | 6 | OK |
| 03 | äußerer Punkt bei a_N = 0,018 bis 0,039·a₀ | 0,017693 | OK |
| 03 | äußerer Punkt bei a_N = 0,018 bis 0,039·a₀ | 0,039216 | OK |
| 03 | drehen wie Newton allein, im Mittel +0,01 dex | 0,0078052 | OK |
| 03 | Die Regel verlangt 0,40 dex mehr Geschwindigkeit | 0,40217 | OK |
| 03 | die blauen Linsen +0,01 ± 0,02 dex auf der Regel | 0,02398 | OK |
| 03 | die roten innen +0,13 ± 0,02 darüber | 0,13063 | OK |
| 03 | die roten innen +0,13 ± 0,02 darüber | 0,024146 | OK |
| 03 | (15 Bins in g_bar, Brouwer 2021) | 15 | OK |
| 03 | mit 0,003·a₀ 534 | 533,69 | OK |
| 03 | mit 0,01·a₀ 996 | 995,97 | OK |
| 03 | χ² = 27483 statt 503 | 27483 | OK |
| 03 | Naht (30–66 kpc) | 30,379 | OK |
| 03 | Naht (30–66 kpc) | 66,01 | OK |
| 03 | isolierten +0,04 ± 0,04 dex | 0,038145 | OK |
| 03 | isolierten +0,04 ± 0,04 dex | 0,040404 | OK |
| 03 | die 22 klassischen Satelliten | 22 | OK |
| 03 | +0,15 ± 0,04 dex ohne Wirtsfeld | 0,14523 | OK |
| 03 | +0,15 ± 0,04 dex ohne Wirtsfeld | 0,0437 | OK |
| 03 | +0,46 ± 0,05 mit ihm | 0,45901 | OK |
| 03 | +0,46 ± 0,05 mit ihm | 0,048514 | OK |
| 03 | (χ² = 76 gegen 341) | 76,497 | OK |
| 03 | (χ² = 76 gegen 341) | 340,98 | OK |
| 03 | gasreichen isolierten 0,00 ± 0,04 dex | -0,0011822 | OK |
| 03 | gasreichen isolierten 0,00 ± 0,04 dex | 0,040615 | OK |
| 03 | die gasarmen isolierten +0,13 ± 0,07 | 0,12991 | OK |
| 03 | die gasarmen isolierten +0,13 ± 0,07 | 0,071752 | OK |
| 03 | die 36 klassischen Satelliten +0,11 ± 0,03 | 36 | OK |
| 03 | die 36 klassischen Satelliten +0,11 ± 0,03 | 0,1096 | OK |
| 03 | die 36 klassischen Satelliten +0,11 ± 0,03 | 0,030478 | OK |
| 03 | die ultraschwachen +0,37 dex | 0,3668 | OK |
| 03 | liest Υ_V = 1,96 für die gasreichen | 1,9647 | OK |
| 03 | 5,2 für die klassischen Satelliten | 5,1942 | OK |
| 03 | weder das Verhältnis der Dichten zum Wirt | 0,024285 | OK |
| 03 | noch das Perizentrum der Bahn | -0,0097672 | OK |
| 03 | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 2,0906 | OK |
| 03 | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 0,28449 | OK |
| 03 | lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge) | 188 | OK |
| 03 | liegen dann −0,002 ± 0,028 dex auf der Regel und lesen Υ_V = 1,93 | -0,0016533 | OK |
| 03 | liegen dann −0,002 ± 0,028 dex auf der Regel und lesen Υ_V = 1,93 | 0,028347 | OK |
| 03 | liegen dann −0,002 ± 0,028 dex auf der Regel und lesen Υ_V = 1,93 | 1,9292 | OK |
| 03 | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | 7 | OK |
| 03 | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | -1,0305 | OK |
| 03 | die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20) | 0,19643 | OK |
| 03 | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,04 ± 0,03 | 27 | OK |
| 03 | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,04 ± 0,03 | -0,041832 | OK |
| 03 | Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,04 ± 0,03 | 0,025948 | OK |
| 03 | die 27 ultraschwachen +0,05 ± 0,04 | 27 | OK |
| 03 | die 27 ultraschwachen +0,05 ± 0,04 | 0,047817 | OK |
| 03 | die 27 ultraschwachen +0,05 ± 0,04 | 0,037299 | OK |
| 03 | die gasarmen isolierten +0,04 ± 0,08 | 0,035615 | OK |
| 03 | die gasarmen isolierten +0,04 ± 0,08 | 0,08079 | OK |
| 03 | ebenso −0,03 und +0,06 | -0,030422 | OK |
| 03 | ebenso −0,03 und +0,06 | 0,057619 | OK |
| 03 | (36 Satelliten der Milchstraße mit Eigenbewegung) | 36 | OK |
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
| F:ganzes | ν_G(−ln q) − 1 = 16,8 bei Q = 2π − 1 | 16,846 | OK |
| F:ganzes | 17,2 mit dem gemessenen Q | 17,196 | OK |
| F:ganzes | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 0,0071084 | OK |
| F:ganzes | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 0,02398 | OK |
| F:ganzes | die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21) | 20,615 | OK |
| F:ganzes | Über alle 15 Bins | 15 | OK |
| F:ganzes | die späten +0,02 ± 0,02 | 0,023309 | OK |
| F:ganzes | die späten +0,02 ± 0,02 | 0,022461 | OK |
| F:ganzes | die roten +0,23 ± 0,01 darüber (χ² 387) | 0,23069 | OK |
| F:ganzes | die roten +0,23 ± 0,01 darüber (χ² 387) | 0,012377 | OK |
| F:ganzes | die roten +0,23 ± 0,01 darüber (χ² 387) | 387,06 | OK |
| F:ganzes | die frühen +0,21 ± 0,01 | 0,20705 | OK |
| F:ganzes | die frühen +0,21 ± 0,01 | 0,013496 | OK |
| F:ganzes | liegen die roten +0,13 ± 0,02 | 0,13063 | OK |
| F:ganzes | liegen die roten +0,13 ± 0,02 | 0,024146 | OK |
| F:ganzes | die frühen +0,14 ± 0,02 über der Regel | 0,1356 | OK |
| F:ganzes | die frühen +0,14 ± 0,02 über der Regel | 0,024879 | OK |
| F:ganzes | wächst es auf +0,29 und +0,26 | 0,29199 | OK |
| F:ganzes | wächst es auf +0,29 und +0,26 | 0,26448 | OK |
| F:ganzes | am Gasanteil (Rangkorrelation −0,04) | -0,036389 | OK |
| F:ganzes | noch am Typ (−0,10) | -0,095828 | OK |
| F:ganzes | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | -0,31212 | OK |
| F:ganzes | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | 168,74 | OK |
| F:ganzes | verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21) | 20,615 | OK |
| F:ganzes | und die späten um −0,26 (157 gegen 14) | -0,2575 | OK |
| F:ganzes | und die späten um −0,26 (157 gegen 14) | 157,15 | OK |
| F:ganzes | und die späten um −0,26 (157 gegen 14) | 14,488 | OK |
| F:ganzes | trifft die roten im Mittel (χ² 146 gegen 387) | 146,32 | OK |
| F:ganzes | trifft die roten im Mittel (χ² 146 gegen 387) | 387,06 | OK |
| F:ganzes | und die frühen (85 gegen 268) | 85,364 | OK |
| F:ganzes | und die frühen (85 gegen 268) | 267,75 | OK |
| F:ganzes | mit χ² 34 bei 20 Werten | 34,369 | OK |
| F:ganzes | (bei R₅₀₀ 3,4 bis 7,0 keV) | 3,3822 | OK |
| F:ganzes | (bei R₅₀₀ 3,4 bis 7,0 keV) | 6,9556 | OK |
| F:ganzes | nur 0,4 % des Gases innerhalb R₅₀₀ | 0,37249 | OK |
| F:ganzes | höchstens 3,7 % | 3,6704 | OK |
| F:ganzes | in 13,8 Gyr nur | 13,8 | OK |
| F:ganzes | die Glättungslänge ℓ = 349 kAU | 349,19 | OK |
| F:ganzes | (ℓ/2^{1/3} = 277 kAU) | 277,15 | OK |
| F:ganzes | bei 30 kAU (1,5 M☉) ν = 3,8 | 3,8003 | OK |
| F:ganzes | die 1,95-fache Geschwindigkeit | 1,9494 | OK |
| F:ganzes | (Gaia DR3, bis 130 pc, 1 bis 30 kAU) | 130 | OK |
| F:ganzes | der Zuwachs von 20 % nach MOND fehlt | 20 | OK |
| 03 | die Glättungslänge ℓ = 349 kAU | 349,19 | OK |
| 03 | jedes gebundenen Paares (277 kAU) | 277,15 | OK |
| F:ganzes | mit χ² 34 bei 20 Werten | 20 | OK |
| F:ganzes | im Mittel je Bin −0,03 bis +0,03 dex | -0,027279 | OK |
| F:ganzes | im Mittel je Bin −0,03 bis +0,03 dex | 0,032415 | OK |
| F:ganzes | verlangt F = 1,8 % | 1,7808 | OK |
| F:ganzes | mit a₀ = 1,06 10,2 ± 1,8 % | 10,202 | OK |
| F:ganzes | mit a₀ = 1,06 10,2 ± 1,8 % | 1,8185 | OK |
| F:ganzes | also F = 1,82 ± 0,32 % | 1,8168 | OK |
| F:ganzes | also F = 1,82 ± 0,32 % | 0,32385 | OK |
| F:ganzes | mit freiem a₀ 7,9 ± 3,2 % | 7,9227 | OK |
| F:ganzes | mit freiem a₀ 7,9 ± 3,2 % | 3,2196 | OK |
| F:ganzes | F = 1,4 ± 0,6 % | 1,4109 | OK |
| F:ganzes | F = 1,4 ± 0,6 % | 0,57337 | OK |
| F:ganzes | Baryonen um z = 0,4 als Gas | 0,39831 | OK |
| 03 | gibt er a₀ = 1,08 | 1,0825 | OK |
| 03 | (H₀ = 68,0) | 68,047 | OK |
| 03 | Δχ² = 5,9 neben ihm | 5,863 | OK |
| 03 | statt 21,4 mit dem gemessenen ω_c | 21,438 | OK |
| 03 | liest Q = 5,22 | 5,2194 | OK |
| 03 | Q = 5,22 ± 0,05 | 0,04843 | OK |
| 03 | 2π − 1 liegt Δχ² = 1,7 daneben | 1,7106 | OK |
| 03 | Verhältnis Δχ² = 8,7 | 8,6694 | OK |
| 03 | ist Q = 5,26 | 5,2627 | OK |
| 03 | Q = 5,26 ± 0,04 | 0,044922 | OK |
| 03 | 2π − 1 liegt Δχ² = 0,2 daneben | 0,20759 | OK |
| 03 | das gemessene Verhältnis 5,0 | 4,9963 | OK |
| 03 | zusammen lesen a₀ = 1,088 | 1,0883 | OK |
| 03 | (1,084–1,092) | 1,0842 | OK |
| 03 | 1,084–1,092) | 1,0925 | OK |
| 03 | mit den Spitzen 1,085 | 1,0854 | OK |
| 03 | (1,082–1,089) | 1,0815 | OK |
| 03 | 1,082–1,089) | 1,0892 | OK |
| 03 | bei H₀ = 68,2 | 68,153 | OK |
| 03 | und Ω_m = 0,30 | 0,30196 | OK |
| 03 | Abstände von 79 der 141 Galaxien | 79 | OK |
| 03 | 79 der 141 Galaxien | 141 | OK |
| 03 | Hubble-Fluss mit H₀ = 73 | 73 | OK |
| 03 | H₀ = 73: a₀ = 1,085 | 1,0854 | OK |
| 03 | (0,98–1,19) | 0,98421 | OK |
| 03 | 0,98–1,19), bei | 1,19 | OK |
| 03 | bei Υ_d = 0,66 | 0,663 | OK |
| 03 | w₀ = −0,85 | -0,85297 | OK |
| 03 | w_a = −0,49 | -0,49108 | OK |
| 03 | Δχ² = 8,5 gegen den Formteil allein | 8,5189 | OK |
| 03 | gemessenen Q Δχ² = 12,7 | 12,652 | OK |
| 03 | liegt dann bei z = 0,43 | 0,42734 | OK |
| 03 | (0,34–0,55) | 0,34063 | OK |
| 03 | 0,34–0,55), zwischen | 0,54745 | OK |
| 03 | Faltung bei z = 0,32 | 0,32218 | OK |
| 03 | beschleunigten Dehnung bei z = 0,67 | 0,66584 | OK |
| 03 | Der Standardwert 1,06 | 1,06 | OK |
| 09 | g_p = 5,586 ist Inhalt | 5,5857 | OK |
| 09 | Ohne Rückkehr gibt sie 1418,95 MHz | 1419 | OK |
| 09 | mit der ersten 1420,60 MHz | 1420,6 | OK |
| 09 | gemessen sind 1420,41 MHz | 1420,4 | OK |
| 09 | vom g-Faktor zu 0,88·α/2π | 0,88127 | OK |
| 09 | die restlichen 1,4·10⁻⁴ | 1,3775 | OK |
| F:ganzes | schließt nach 2, 8, 20, 40, 70 (κ = 0,1 | 1 | OK |
| F:ganzes | Mit μ = 0,02 statt 0,03 erscheint 126, dann fehlt 28 | 1 | OK |
| 02 | 850 ungerade Kerne mit eindeutig gemessenem Spin | 850 | OK |
| 02 | 321 Spins mit Parität (38 %), 641 Paritäten (75 %) | 3,2164·10⁵ | OK |
| 02 | 606 Kerne) 299 Spins und 520 Paritäten (49 %, 86 %) | 2,9952·10⁸ | OK |
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
| 03 | (H₀ = 68,9) sind es 0,14 mag | 0,14059 | OK |
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
| 00 | 1,3 σ neben Plancks Wert | 1,2563 | OK |
| 03 | Pantheon+, 1590 Supernovae | 1590 | OK |
| 03 | nur Δχ² = 2,0 über | 1,9631 | OK |
| 03 | ΛCDM (Ω_m = 0,304) | 0,30393 | OK |
| 03 | Ω_m = 0,304), 1,4 σ | 1,4011 | OK |
| 03 | (w₀ = −0,89, | -0,89032 | OK |
| 03 | w_a = −0,21) | -0,21286 | OK |
| 03 | Δχ² = 6,9 unter | 6,871 | OK |
| 03 | dem Formteil, 1,8 σ | 1,7736 | OK |
| 03 | r_d = 149,0 Mpc | 149,05 | OK |
| 03 | der Formteil trifft es mit χ² = 6,8 | 6,7556 | OK |
| 08 | (10⁶ M☉: s² = 1,4·10¹⁷ | 1,4355 | OK |
| 08 | am Horizont, r* = 36 pc) | 36,262 | OK |
| 08 | = 1 − ⅔√2 = 5,7 % | 5,7191 | OK |
| 08 | sechsmal so viel wie Kernbindung bis zum Eisen | 6,0604 | OK |
| 08 | (j = 1) sind es 1 − 1/√3 = 42 % | 42,265 | OK |
| 08 | ist T/T_I = M_U/(8πM) | 1 | OK |
| 08 | r₀ = Λ/2π, dieselbe Zählung | 1 | OK |
| 08 | eine 2,3·10²² mal kürzere Periode | 2,2847 | OK |
| 08 | Erst bei M_U/4 = 1,4·10²³ M☉ | 1,4355 | OK |
| 08 | Am Horizont ist s² = M_U/(4M) | 1 | OK |
| 08 | kT = m̃c²/(4π). Gegen | 1 | OK |
| 03 | (r_d = 147,1 Mpc aus ω_b und ω_c) | 147,08 | OK |
| 03 | liest DESI dann a₀ = 1,10 | 1,0997 | OK |
| 03 | a₀ = 1,10·10⁻¹⁰ m/s² (1,09– | 1,092 | OK |
| 03 | m/s² (1,09–1,11 | 1,1073 | OK |
| 03 | mit H₀ = 68,9 bei χ² = 10,4 | 68,933 | OK |
| 03 | bei χ² = 10,4 für 13 Werte | 10,375 | OK |
| 03 | für 13 Werte | 13 | OK |
| 03 | der Standardwert 1,06 liegt Δχ² = 27 | 27,42 | OK |
| 03 | Mit Q = 2π − 1 statt des gemessenen ebenso 1,10 | 1,1008 | OK |
| 03 | (r_d = 147,6 Mpc, H₀ = 68,9 | 147,55 | OK |
| 03 | r_d = 147,6 Mpc, H₀ = 68,9 | 68,851 | OK |
| 03 | H₀ = 68,9, χ² = 10,3) | 10,289 | OK |
| 03 | ist Q = 5,5 (5,3–5,7) | 5,4754 | OK |
| 03 | Q = 5,5 (5,3– | 5,2772 | OK |
| 03 | (5,3–5,7); 2π − 1 | 5,6736 | OK |
| 03 | 2π − 1 liegt Δχ² = 0,9 daneben | 0,94025 | OK |
| 03 | an sechs Rotverschiebungen | 6 | OK |
| 01 | (PubChem, 13,598 eV) | 13,598 | OK |
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
| 02 | auf 0,32 Protonen genau | 0,31507 | OK |
| 02 | 219 Massenzahlen bis A = 238 | 219 | OK |
| 02 | erst bei Z = 134 (A = 370) | 133,66 | OK |
| 02 | Z = 134 (A = 370) | 370,27 | OK |
| 02 | erst bei Z = 134 | 134 | OK |
| 02 | Schließzahlen des Modells bei Z = 181 bis 191 | 181,25 | OK |
| 02 | bei Z = 181 bis 191 | 191,12 | OK |
| 02 | sinkt von 3,39 auf 1,72 MeV | 3,3921 | OK |
| 02 | auf 1,72 MeV, das leichteste | 1,7213 | OK |
| 02 | trifft sie auf 0,26 statt 0,32 | 0,26484 | OK |
| 02 | Massen nur bis Z = 110 | 110 | OK |
| 02 | Die Hülle hält bis Z = 170 bis 171 | 170,06 | OK |
| 02 | Hülle hält bis Z = 170 bis 171: | 170,93 | OK |
| 02 | sind die Elemente bis Z = 133 | 133 | OK |
| 02 | mit N ≈ 202–235 | 1 | OK |
| 02 | reichen sie bis zur Grenze der Hülle, Z = 170 | 170 | OK |
| 02 | (x = 1,06 ohne | 1,0614 | OK |
| 02 | 0,84 mit Schalenglied) | 0,83917 | OK |
| 02 | gegen den Alpha-Zerfall Sekunden bis Stunden | 1 | OK |
| 02 | (36 Paare, AME2020) | 36 | OK |
| 02 | Zα = 1 (Z = 137) gilt für einen punktförmigen Kern | 137 | OK |
| 02 | an 82 gerade-geraden Alpha-Strahlern | 82 | OK |
| 02 | trifft das auf 0,51 Dekaden | 0,51272 | OK |
| 02 | Mit dem Edelgas Og schließt die siebte Periode bei Z = 118 (PubChem) | 1 | OK |
| 02 | die Perioden sind 2, 8, 8, 18, 18, 32, 32 lang | 1 | OK |
| 02 | 99 der 118 Elemente füllen genau nach n + l | 99 | OK |
| 02 | die 19 Abweichungen | 19 | OK |
| 02 | Acht schließen dabei eine halbe oder volle Schale | 8 | OK |
| 02 | sieben beginnen die f-Schale mit einem d-Elektron | 7 | OK |
| 02 | vier geben ein s-Elektron an die d-Schale | 4 | OK |
| 02 | mit 8s (Z = 119, 120) und der ersten g-Schale | 1 | OK |
| 02 | und schlösse bei 168 | 168 | OK |
| 02 | bindet an der Schwelle E = 0 genau bei den Kopplungen N(N + 1) | 1 | OK |
| 02 | Sein Feld folgt der Spiegelform innen auf 0,6 % | 0,55266 | OK |
| 02 | seine Höhe v = 3,44·Z^{2/3} liegt von Ne bis Og im Fenster | 1 | OK |
| 02 | (104,8 zwischen 90 und 110) | 1 | OK |
| 02 | Seine Zustände an der Schwelle sind Kugelfunktionen, u = √r·P_N^{2l+1}((r − R)/(r + R)) | 1 | OK |
| 01 | Spiegelfeld der Hülle (Kugelfunktionen auf S², Grad n + l, Ordnung 2l + 1) | 1 | OK |
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
| 03 | M_Gas = 0,6–1,4·10⁷ M☉ | 0,58961 | OK |
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
| 08 | ε_L = ℓ_P/ℓ₀ = 1,9·10⁻⁶² | 1,9062 | OK |
| 08 | ε_L^{1/3}·m_Pl = 5,8·10⁻²⁹ kg | 5,8139 | OK |
| 08 | c³τ_D/G = 1,1·10⁵⁴ kg | 1,1418 | OK |
| 08 | Δs(M_U) = ⅔·Λ | 2 | OK |
| 08 | (4π⁴/15)·ε_L² = 0,94·10⁻¹²² | 0,94388 | OK |
| 02 | r* = 32 fm | 32,453 | OK |
| 02 | a_B = 1631·r* | 1630,6 | OK |
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
| 07 | (Euler; bei q = ½ beide 2,384) | 2,3842 | OK |
| 07 | bei q = ½ beide 2,384) | 2,3842 | OK |
| 01 | Abgabe der Windung quer zur Schwingung 2/3 | 0,66667 | OK |
| F:ganzes | denselben vollen Umlauf geben, ν = 2π | 1 | OK |
| F:ganzes | d/π = 0,955 | 0,95493 | OK |
| F:ganzes | √(d/π) = 0,977 | 0,97721 | OK |
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
| 03 | (ω_b = 0,02238 | 0,022379 | OK |
| 03 | ω_c = 0,1199, H₀ | 0,11991 | OK |
| 03 | H₀ = 67,37). An | 67,367 | OK |
| 03 | trifft Plancks Werte | 0,058587 | OK |
| 03 | trifft Plancks Werte | 0,078874 | OK |
| 03 | allein liegt das Modell Δχ² = 1,0 über | 0,9953 | OK |
| 03 | an allen Daten 0,01; | 0,014855 | OK |
| 03 | an allen Daten 0,01; | 0,014855 | OK |
| 03 | DESI mit χ² = 18,6 | 18,608 | OK |
| 03 | 18,6 (13 Werte) | 13 | OK |
| 03 | ΛCDM mit 31,4 | 31,438 | OK |
| 03 | liest die Schwelle 1,081 ± 0,003 | 1,0815 | OK |
| 03 | liest die Schwelle 1,081 ± 0,003 | 0,002952 | OK |
| 03 | der Standardwert liegt Δχ² = 52,6 | 52,591 | OK |
| 03 | gewinnt das Modell Δχ² = 9,6 | 9,6195 | OK |
| 00 | w_a < 0, Δχ² = 9,6) | 9,6195 | OK |
| 03 | w₀ = −0,84 ± 0,05 | -0,83943 | OK |
| 03 | w₀ = −0,84 ± 0,05 | 0,052894 | OK |
| 03 | w_a = −0,56 ± 0,18 | -0,56368 | OK |
| 03 | w_a = −0,56 ± 0,18 | 0,18311 | OK |
| 03 | bei z = 0,40 ± 0,07 | 0,39831 | OK |
| 03 | bei z = 0,40 ± 0,07 | 0,073489 | OK |
| 03 | beschleunigten Dehnung (z = 0,65) | 0,64557 | OK |
| 03 | Faltung (z = 0,31) | 0,30609 | OK |
| 03 | DESI liegt dann bei χ² = 9,3 | 9,3031 | OK |
| 03 | die Schwelle liest 1,070 ± 0,013 | 1,0704 | OK |
| 03 | die Schwelle liest 1,070 ± 0,013 | 0,013404 | OK |
| 03 | der Standardwert nur Δχ² = 0,8 daneben | 0,83463 | OK |
| 03 | liegt Δχ² = 8,8 darunter | 8,77 | OK |
| 03 | liegt noch 2,5 tiefer | 2,4842 | OK |
| 03 | nur auf 1,3 % genau | 1,2522 | OK |
| 03 | liest ΛCDM selbst Q = 5,278 | 5,2778 | OK |
| 03 | an den Spektren allein 5,367 | 5,3672 | OK |
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
| 02 | bei 126 (2,03 statt | 2,0267 | OK |
| 02 | 2,03 statt 1,62 MeV | 1,6224 | OK |
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
| AB | abhandlung_de.md: alle 181 Zahlen stehen im Modell | 181 | OK |
| AB | abhandlung_en.md: alle 181 Zahlen stehen im Modell | 181 | OK |
| RB | rechnen_de.md: alle 156 gerechneten Zahlen stehen | 156 | OK |
| RB | rechnen_en.md: alle 156 gerechneten Zahlen stehen | 156 | OK |
| EL | elemente_de.md entspricht der Rechnung | 1 | OK |
| EL | elemente_en.md entspricht der Rechnung | 1 | OK |

**Ergebnis:** alles OK · Bahndaten in `berichte/bericht_0252/` · Laufzeit 133 s
