# 07 Feldgleichung

Die klassischen Theorien sind Markow-Prozesse: Der nächste Zustand hängt nur vom jetzigen ab. Im Modell hängt die Bewegung zusätzlich von der Rückkehr in der inneren Richtung ab. Diese Rückkehr ändert die Bewegungsregel. Sie fügt keine Quelle hinzu.

## Die eine Rückkehrgleichung

Die Rückkehr ist eine Irrfahrt mit einem Takt in der inneren Richtung ϖ (Glättungstiefe, eine Halbgerade vom Ort aus). Je Schritt geht der Anteil verloren, den das Taktgefälle am Ort vorgibt:

```
∂_t u + (Γ + θ)·u = D·∂_ϖ²u,     Γ = a_N/c,     D = GM/c
```

Γ ist der Verlust aus dem Gefälle am Ort, θ der aus der Ausdehnung des Bereichs. Stationär, ohne Ausdehnung, bleibt D·G″ = Γ·G mit G(0) = 1, und die Rückkehr an die Tiefe ϖ_R ist

```
q = G(ϖ_R) = e^{−ϖ_R·√(Γ/D)} = e^{−s},     s² = Γ·τ_D = a_N/a₀,     τ_D = ϖ_R²/D = c/a₀
```

Mit D = GM/c ist die Tiefe der Rückkehr ϖ_R = √(D·τ_D) = √(GM/a₀) = r*: Die innere Richtung reicht so tief, wie die Schwelle im Raum weit ist, und s = r*/r.

Am Ort gezählt, jeder Rücklauf ein Ereignis in einem Sinn, ergibt die Erneuerung u(0) = g_N + u(ϖ_R) die Regel:

```
a = ν(|g_N|)·g_N,     ν = Σ_{n≥0} qⁿ = 1/(1 − q) = 1 + n_B(s),     n_B = 1/(e^s − 1)
```

Die Verstärkung hat die Form eins plus einer Bose-Besetzung. Gezählt werden Rückläufe, keine Teilchen. Tief (s → 0) gibt sie v⁴ = GMa₀, hoch (s → ∞) exakt Newton. Der Schalter χ = s/(e^s − 1) läuft von 0 (Newton) nach 1 (tiefer Bereich). Die Wurzel verlangt eine eindimensionale Rückkehr. Genauer ist sie der erste Durchgang einer Irrfahrt in einer Dimension: Die Zeit T, in der sie die Strecke der Rückkehr zurücklegt, ist verteilt wie T^{−3/2}, und ihr Mittel ⟨e^{−x·T}⟩ ist e^{−√x} mit x = a_N/a₀, über die Verteilung integriert gleich auf 10⁻¹⁵. Nur eine Dimension trägt den Tiefbereich: In zwei Dimensionen wüchse die Summe der Rückläufe ohne Verlust nur mit dem Logarithmus, in drei kehrte die Irrfahrt nur zu 0,341 zurück (Pólya), die Summe bliebe bei 1/(1 − 0,341) = 1,516, und flache Kurven gäbe es nicht. An den SPARC-Rotationskurven gelesen ist ihre Dimension 0,92 ± 0,05.

## Die Schwelle aus der Schließung

Eine exakt offene Richtung hätte unendliche mittlere Rückkehrzeit. Die innere Richtung ist bei Λ geschlossen. Nach Kac ist die mittlere Rückkehrzeit zum Ort dann ⟨T⟩ = 2Λ/c, unabhängig von der Masse. Mit τ_D = ⟨T⟩ folgt

```
a₀ = c²/(2Λ),     r*² = 6Δs·Λ = 3Δs·ℓ₀,     ℓ₀ = 2Λ
```

Bei r* sind Verlustzeit 1/Γ und Rückkehrzeit τ_D gleich. Die Rückkehr ist kein Vorgang längs der Bahn, sondern eine stationäre Schleife am festen Ort, denn τ_D ist länger als das Weltalter. Die Bahn läuft ohne Nachlauf.

## Die Rückkehr als Spiegelbild

Liest man die Zeit der Rückkehr als imaginäre Zeit, ist die Rückkehrgleichung die Schrödinger-Gleichung in der inneren Richtung ϖ, mit D = ħ/(2m̃) und der Neigung als Potential, V = ħΓ = ħa_N/c. Der Trichter M kehrt zurück wie die Quantenbewegung seiner Spiegelmasse m̃ = ħc/(2GM) = m_Pl²/(2M), und deren Quantenlänge ist sein Grenzradius, λ̄(m̃) = 2GM/c² = r₀. Der Rückkehranteil q = e^{−s} ist der Tunnelfaktor dieses Spiegelbilds durch die Neigung, dieselbe Form wie e^{−√(E_G/E)} beim Durchgang zweier Kerne durch ihre Coulomb-Schwelle. Darum ist die Periode der Rückkehr, τ_D, eine Temperatur, kT_I = ħ/τ_D, und die Regel hat die Bose-Form: Wärme ist Periodizität in imaginärer Zeit. Gravitation und Quant sind so zwei Ansichten desselben Spiegels, M in der Raumrichtung, m̃ = m_Pl²/(2M) in der inneren.

## Die Formen der Regel

Dieselbe Funktion bildet eine Leiter, jede Stufe die vorige mal s: s^d·n_B(s) mit d = 0 (Summe der Rückläufe), 1 (Schalter χ), 2 (Zusatzbeschleunigung), 3 (Strukturwirkung). Integriert gibt jede Stufe Γ(d + 1)ζ(d + 1): π²/6, 2ζ(3), π⁴/15. Bei r* (s = 1) sind alle Stufen gleich. Ihre Maxima sind Wiens Verschiebungskonstanten in d Dimensionen, d(1 − e^{−s}) = s. Der Tiefbereich ist die Rayleigh-Jeans-Grenze, die Newton-Grenze der Wien-Schwanz. Mit der halben Drehung der Windung, q → −q, wird jede Stufe der Leiter s^d·n(s) mit 1 − 2^{−d} gewichtet: ½, ¾ und ⅞ für d = 1, 2, 3. Die Stufe d = 0, die Summe aller Rückläufe, bleibt endlich, ∫n_F ds = ln 2 = 0,693: ein Bit, die zwei Sinne einer Windung. Bei der Faltung wächst sie ohne Grenze, das ist der Tiefbereich. Die Rückläufe einer Windung summieren sich höchstens zu ½: Auch mit Rückkehr hätte die Windung keinen Tiefbereich.

**Gleichgewicht, Nullpunkt und Spiegel der Schleifen:** Die Regel ist ein Gleichgewicht nach Einstein: Abgabe ∝ 1 + n_B, Aufnahme ∝ n_B, ihr Verhältnis ist der Rückkehranteil, q = n_B/ν = e^{−s}. Newton ist der spontane Teil, die Rückkehr der angeregte. Kehrt man den Sinn der inneren Richtung um, tauschen beide die Rollen, ν(s) + ν(−s) = 1. Zerlegt ist ν = ½ + ½·coth(s/2): eine feste Hälfte und der Schwankungsteil eines Schwingers. Die feste Hälfte ist das ½ in der Umrechnung in Masse. Die Regel zählt Rückläufe einer Schleife. Zählt man alle Vielfachen der Schleife zusammen, Z = ∏ 1/(1 − qⁿ), ist das bis auf q^{1/24} der Kehrwert der η-Funktion des Torus, und deren Spiegelung, η(−1/τ) = √(−iτ)·η(τ), bildet den Newton-Bereich auf den tiefen ab, s ↔ 4π²/s, mit dem Spiegelpunkt s = 2π, einem vollen Umlauf. Mit s = −ln|q| = 2π·Im τ liegt der Treffpunkt der Leiter, s = 1 bei r*, bei Im τ = 1/2π, einem Radiant. Der Spiegelpunkt bei Im τ = 1, dem vollen Umlauf: Die Schwelle zählt je Radiant, der Spiegel je Umlauf, wie beim Ganzen, dessen offener Teil ein Radiant eines Umlaufs ist. Tief ist ln Z = π²/(6s) + ½·ln(s/2π) − s/24: der führende Teil ist die erste Stufe der Leiter, π²/6. Bei der Windung, die sperrt, zählen alle Vielfachen ∏(1 + qⁿ), tief ln Z = π²/(12s) − ½·ln 2 + s/24, beide bis auf Glieder in e^{−4π²/s} bzw. e^{−2π²/s}. Die halbe Drehung halbiert den führenden Teil und kehrt den Nullpunkt um, +s/24 gegen −s/24: Wie Ladungen sich aufheben und Massen sich addieren, heben sich die Nullpunkte von Windung und Faltung auf. Der Spiegel bildet die Sperre auf die ungeraden Moden ab, ∏(1 + qⁿ) = (q·q̃)^{−1/24}·∏(1 − q̃^{2n−1})/√2 mit q̃ = e^{−2π²/s}. Der Spiegelpunkt der Windung liegt bei s = π√2, der der Faltung bei 2π. Windungen, die einander sperren, besetzen jede Mode höchstens einmal, ∏(1 + qⁿ). Rückläufe, die nur in ungeraden Moden kommen, zählen ∏1/(1 − q^{2n−1}). Beide sind gleich (Euler, bei q = ½ beide 2,384): Die Sperre der halben Drehung ist ein voller Umlauf auf dem doppelten Weg, der nach einem Umlauf das Vorzeichen wechselt.

**Entropie und Temperatur:** Die geometrische Verteilung der Rücklaufzahl hat bei festem Mittel die größte Entropie, H = ln ν + χ. Ihre Schwankung ist ⟨Δn²⟩ = n_B·ν = n_B + n_B² = 1/(4·sinh²(s/2)), Einsteins Form der Schwankung der Strahlung: n_B zählt die Rückläufe einzeln wie unabhängige Ereignisse, n_B² ihr Bündeln. Die Verstärkung ist das Verhältnis von Schwankung zu Mittel, ν = ⟨Δn²⟩/n_B, und Newton (ν = 1) ist die Grenze unabhängiger Rückläufe. Die Schwankung ist die Antwort der Regel auf die Tiefe, ⟨Δn²⟩ = −dν/ds = −d(½·coth(s/2))/ds: Nur der Teil des Schwingers antwortet, die feste Hälfte nicht. So ist d ln ν/d ln s = −s·⟨Δn²⟩/ν = −χ: Der Schalter, der die Steigung der Relation stellt, ist die Schwankung der Rückläufe. Sie gehört der Verteilung in der inneren Richtung, am Ort wirkt nur ihr Mittel ν. Beide Teile sind gleich, wo die Hälfte zurückkehrt, q = ½, bei a_N = (ln 2)²·a₀ = 0,48·a₀. Dort trägt der angeregte Teil so viel wie der spontane, die Rückkehr so viel wie Newton, ν = 2. Bei größerem a_N überwiegen die einzelnen Rückläufe, bei kleinerem das Bündeln. Mit der Sperre der halben Drehung, jede Mode höchstens einmal besetzt, ist sie n_F·(1 − n_F) = −dn_F/ds ≤ ¼: Die Windung schwankt begrenzt, die Faltung ohne Grenze. Die Schwelle vergleicht die Unruh-Temperatur des ruhenden Beobachters am Ort mit der Periode der inneren Richtung: s² = 2π·T_U/T_I, kT_I = ħa₀/c. Eine Periode, keine Wärme. Ebenso ist der Verlust aus der Ausdehnung eine Temperatur gegen dieselbe Periode: Eine Ausdehnung mit der Rate H gibt H·τ_D = 2π·T_H/T_I = 2π/Θ, mit der Temperatur des Horizonts kT_H = ħH/(2π).

**Kern-Nome und Bahntorus:** Der Bahntorus hat den Modul τ_M mit Re τ_M = ω_rad/ω_azim (Umlaufphase) und Im τ_M (Bindung). Für die Punktmasse ist (Re τ_M)² = 1 + χ: χ = 0 gibt geschlossene Keplerbahnen, χ = 1 Rosetten mit 254,6° zwischen den Perizentren. Die Kräfte r^a und r^A sind Spiegelpartner, wenn (a + 3)(A + 3) = 4: Die Abbildung z → z^{(a+3)/2} führt die Bahnen der einen in die der anderen über. Newton (r⁻²) und der gehaltene Körper (r¹) sind ein solches Paar, z → z². Der Spiegelpunkt der gebundenen Bahnen ist a = −1, die Kraft des tiefen Bereichs (der zweite, a = −5, führt in den Fall ins Zentrum): Die Regel läuft von Newton zum Treffpunkt der Kraftformen. Die Regel hebt die Kepler-Entartung κ = Ω auf und ersetzt sie durch örtliche Resonanzen bei rationalem κ/Ω. Der Kern-Nome ist |q_Δs| = Δs/r, die Spitzenform Δ = q_Δs·∏(1 − q_Δsⁿ)²⁴ verschwindet für q_Δs → 0. Die Bindung läuft nur über den Betrag des Nomes. Ginge die Phase mit ein, fiele die Verstärkung außen unter 1.

**Zwei Nome:** Der Autor lässt offen, welcher Nome den Imaginärteil trägt, |q| = e^{−r*/r} oder der Kern-Nome |q_Δs| = Δs/r. Mit dem ersten ist Im τ = r*/(2πr) die Tiefe der Rückkehr gegen den Umfang am Ort, der Modul des Torus aus Umlauf und innerer Richtung. Mit dem zweiten ist |q_Δs| ein Drittel der Tiefe δ = 3Δs/r des Trichters, Im τ = ln(r/Δs)/2π, und der Zusatz ist χ·Δ(q_Δs)·3c²/r* mit der Spitzenform Δ(q_Δs) = q_Δs·∏(1 − q_Δsⁿ)²⁴. Am Ort sind beide eins zu eins verbunden, (r*/r)² = 3·(Δs/r)·ℓ₀/r, und der Schalter des zweiten trägt den Nome des ersten, χ/s = Σ_{n≥1} e^{−ns}. Beide geben so dieselbe Regel bis auf den Faktor ∏(1 − q_Δsⁿ)²⁴ ≈ 1 − 8δ vor dem Zusatz. Auf die Regel wirkt er als 1 − 8(v/c)²·q(1 − q), mit der Geschwindigkeit v der Kreisbahn und dem Nome q des Orts oder, in Haufen, des Ganzen, also um höchstens 2(v/c)². An SPARC (2688 Punkte), an der Linse um isolierte Galaxien (KiDS-1000), an HIFLUGCS und X-COP sind das höchstens 2·10⁻⁵ dex, Δχ² = 0,01 an SPARC und 0,06 an HIFLUGCS, mit der ganzen Masse des Feldes als Faltung 1·10⁻⁴ dex: An den Daten sind beide nicht unterscheidbar. Getrennt sind sie in der Form. Mit e^{−r*/r} liegt die Schwelle für jede Masse bei Im τ = 1/2π, das Ganze hat den Nome Q/(1 + Q), und die halbe Drehung q → −q gibt die Fermi-Form der Windung ohne Tiefbereich. Mit Δs/r bleibt |q_Δs| außerhalb von r₀ = 6Δs unter 1/6 und erreicht den Nome des Ganzen nie, die Schwelle liegt je nach Masse bei Im τ = ln(3c²/v_flach²)/2π (Sonne 4,5, Milchstraße 2,6), und die halbe Drehung kehrt nur das Vorzeichen des Zusatzes um, ν ≈ 1 − n_B: Die Windung hätte einen Tiefbereich mit umgekehrtem Vorzeichen, durch null bei s = ln 2.

## Struktur statt Masse

Der Beobachter rechnet den Zusatz in eine Masse um, für die Punktmasse

```
M_schein/M = 1/(1 − e^{−s}) = r/r* + ½ + r*/(12r) − (r*/r)³/720 + …,     Φ_X = v_flach²·ln ν(s)
```

Die Korrekturen ½ und 1/12 sind Bernoulli-Zahlen: der Fehler der einzelnen Rückläufe gegen das Kontinuum, wie beim Vieleck gegen den Kreis. Das Zusatzpotential ist der Logarithmus der Verstärkung und geht jenseits von r* in v_flach²·ln(r/r*) über. In der Masse gelesen wechselt der Anteil der Baryonen an der gelesenen Dichte bei der größten Zusatzbeschleunigung das Vorzeichen: Dünne Materie verstärkt sie, dichte mindert sie.

## Strukturwirkung und Vakuumterm

Die Wirkung der Metrik trägt dieselbe Funktion wie die Regel:

```
Q̃(s) = −4∫_s^∞ x³·n_B(x) dx,     Q̃(0) = −4π⁴/15
```

Die Konstante −4π⁴/15 ist ein globaler Vakuumterm: Er hängt weder an der Glättung noch an einzelnen Trichtern. Er ist der Vakuumteil der Dunklen Energie. Die x³ ist die Zustandsdichte der drei Raumdimensionen.

**Ganz und örtlich gespiegelt:** Dieselbe Funktion hat zwei Lesarten. Über alle s integriert, für das Ganze, gibt jede Stufe Γ(k)ζ(k) mit k = d + 1. Am Ort, nahe dem tiefen Bereich entwickelt, hat sie die Beiwerte ζ(0) = −½, ζ(−1) = −1/12, ζ(−3) = 1/120: n_B(s) = 1/s + Σ ζ(−n)·(−s)ⁿ/n!. Beide verbindet der Spiegel ζ(1 − k) = 2·(2π)^{−k}·cos(πk/2)·Γ(k)ζ(k), mit (2π)^{−k}, der Zählung je Radiant. So ist die Konstante des Vakuumteils, Γ(4)ζ(4) = π⁴/15, das Spiegelbild von ζ(−3) = 1/120, das in der Umrechnung in Masse als −(r*/r)³/720 erscheint, und π²/6 das von −1/12: Was das Ganze trägt und was der Beobachter am Ort als Korrektur liest, ist dieselbe Größe von zwei Seiten.

## Metrik und Regel

Die Regel setzt am Erzeuger der Bewegung an (L_eff = 𝒫L), über die Rückkehrschritte entlang der Bahn. Auf einer Kreisbahn gibt das a = ν·g_N mit der Neigung g_N in der Ebene der Bahn. Ein Erzeuger, der nicht von Zeit und Richtung abhängt, erhält Energie und Drehimpuls. Die Metrik trägt den wirbelfreien Teil der Regel:

```
ν·g_N = −∇Φ + ∇×A,     ∇²Φ = ∇·[ν·∇Φ_N],     Ψ = Φ
```

Licht folgt −∇Φ. Der Wirbel kommt aus ∇×(ν·g_N) = ∇ν × g_N: Er verschwindet, wo ν auf den Flächen gleicher Newton-Tiefe konstant ist, bei Kugel, Faden und Wand, und dort lesen Bahnen und Licht dieselbe Neigung. In einer dünnen Scheibe liest die Feldform ν am vollen Betrag der Neigung an der Oberfläche der Schicht, mit deren senkrechtem Feld 2πGΣ, für die Kuzmin-Scheibe exakt, über der das Feld das eines Punktes im Spiegelbild unter ihr ist: a = ν(√(g_R² + (2πGΣ)²)/a₀)·g_R. Die Bahnen lesen nur die Neigung an ihrem Ort: An SPARC (141 Galaxien, Gas aus V_gas, Sterne aus der Flächenhelligkeit der Scheibe) liegt die Feldform als Bahnregel Δχ² = 368 über ν·g_N, und die Reste von ν·g_N gehen mit ihrer Vorhersage mit der Steigung 0,05 ± 0,12, wo sie 1 verlangt. Darum ist in Scheiben innen die Linsenmasse kleiner als die dynamische. In der Triade ist −∇Φ das Zentrum der Regel und ∇×A ihr Wirbel: Die Bahn liest beide, das Licht nur das Zentrum.

Als Kraft am Ort, ohne Erzeuger, hielte sich keine kalte Scheibe: Mit ν·g_N als Kraft je Körper dehnt sich eine Scheibe an der Schwelle (Streuung 0,15·v_c), sobald sie Arme und Balken bildet, in 3 Gyr auf das 2- bis 3-fache ihres Halbmassenradius und gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern). Mit einem Erzeuger, hier der Feldform, bleibt sie an der Schwelle, ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton.

Der Strukturterm ist eine Reihe in q mit |q| = e^{−s}, die für q → 0 in jeder Ordnung verschwindet. Er trägt zu keinem PPN-Koeffizienten bei, γ = β = 1, und Gravitationswellen laufen mit c. An der Erdbahn ist der Zusatz e^{−7480}.

## Hierarchie der Trichter

Kompakte Trichter, dichter als die Gezeitendichte ρ_J = 3𝒯/(4πG) ihrer Umgebung (𝒯 Gezeit des geglätteten Elternfelds), bleiben bis zu ihrer Glättungslänge ℓ, ℓ³ = GM/𝒯, Newtonsch. Verstärkt wird nur das kollektive, über diese Zonen geglättete Feld. Weil die Regel am Ort auf das geglättete Feld wirkt, in dem das eigene Feld am Schwerpunkt null ist, fallen Körper am selben Ort gleich.
