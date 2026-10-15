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

Am Ort gezählt, jeder Rücklauf ein Ereignis in einem Sinn, ergibt die Erneuerung u(0) = g_N + u(ϖ_R) die Regel:

```
a = ν(|g_N|)·g_N,     ν = Σ_{n≥0} qⁿ = 1/(1 − q) = 1 + n_B(s),     n_B = 1/(e^s − 1)
```

Die Verstärkung hat die Form eins plus einer Bose-Besetzung. Gezählt werden Rückläufe, keine Teilchen. Tief (s → 0) gibt sie v⁴ = GMa₀, hoch (s → ∞) exakt Newton. Der Schalter χ = s/(e^s − 1) läuft von 0 (Newton) nach 1 (tiefer Bereich). Die Wurzel verlangt eine eindimensionale Rückkehr. An den SPARC-Rotationskurven gelesen ist ihre Dimension 0,92 ± 0,05.

## Die Schwelle aus der Schließung

Eine exakt offene Richtung hätte unendliche mittlere Rückkehrzeit. Die innere Richtung ist bei Λ geschlossen. Nach Kac ist die mittlere Rückkehrzeit zum Ort dann ⟨T⟩ = 2Λ/c, unabhängig von der Masse. Mit τ_D = ⟨T⟩ folgt

```
a₀ = c²/(2Λ),     r*² = 6Δs·Λ = 3Δs·ℓ₀,     ℓ₀ = 2Λ
```

Bei r* sind Verlustzeit 1/Γ und Rückkehrzeit τ_D gleich. Die Rückkehr ist kein Vorgang längs der Bahn, sondern eine stationäre Schleife am festen Ort, denn τ_D ist länger als das Weltalter. Die Bahn läuft ohne Nachlauf.

## Die Rückkehr als Spiegelbild

Liest man die Zeit der Rückkehr als imaginäre Zeit, ist die Rückkehrgleichung die Schrödinger-Gleichung in der inneren Richtung ϖ, mit D = ħ/(2m̃) und der Neigung als Potential, V = ħΓ = ħa_N/c. Der Trichter M kehrt zurück wie die Quantenbewegung seiner Spiegelmasse m̃ = ħc/(2GM) = m_Pl²/(2M), und deren Quantenlänge ist sein Grenzradius, λ̄(m̃) = 2GM/c² = r₀. Der Rückkehranteil q = e^{−s} ist der Tunnelfaktor dieses Spiegelbilds durch die Neigung, dieselbe Form wie e^{−√(E_G/E)} beim Durchgang zweier Kerne durch ihre Coulomb-Schwelle. Darum ist die Periode der Rückkehr, τ_D, eine Temperatur, kT_I = ħ/τ_D, und die Regel hat die Bose-Form: Wärme ist Periodizität in imaginärer Zeit. Gravitation und Quant sind so zwei Ansichten desselben Spiegels, M in der Raumrichtung, m_Pl²/M in der inneren.

## Die Formen der Regel

Dieselbe Funktion bildet eine Leiter, jede Stufe die vorige mal s: s^d·n_B(s) mit d = 0 (Summe der Rückläufe), 1 (Schalter χ), 2 (Zusatzbeschleunigung), 3 (Strukturwirkung). Integriert gibt jede Stufe Γ(d + 1)ζ(d + 1): π²/6, 2ζ(3), π⁴/15. Bei r* (s = 1) sind alle Stufen gleich. Ihre Maxima sind Wiens Verschiebungskonstanten in d Dimensionen, d(1 − e^{−s}) = s. Der Tiefbereich ist die Rayleigh-Jeans-Grenze, die Newton-Grenze der Wien-Schwanz. Mit der halben Drehung der Windung, q → −q, wird jede Stufe der Leiter s^d·n(s) mit 1 − 2^{−d} gewichtet: ½, ¾ und ⅞ für d = 1, 2, 3. Die Stufe d = 0, die Summe aller Rückläufe, bleibt endlich, ∫n_F ds = ln 2 = 0,693: ein Bit, die zwei Sinne einer Windung. Bei der Faltung wächst sie ohne Grenze, das ist der Tiefbereich. Die Rückläufe einer Windung summieren sich höchstens zu ½: Auch mit Rückkehr hätte die Windung keinen Tiefbereich.

**Gleichgewicht, Nullpunkt und Spiegel der Schleifen:** Die Regel ist ein Gleichgewicht nach Einstein: Abgabe ∝ 1 + n_B, Aufnahme ∝ n_B, ihr Verhältnis ist der Rückkehranteil, q = n_B/ν = e^{−s}. Newton ist der spontane Teil, die Rückkehr der angeregte. Kehrt man den Sinn der inneren Richtung um, tauschen beide die Rollen, ν(s) + ν(−s) = 1. Zerlegt ist ν = ½ + ½·coth(s/2): eine feste Hälfte und der Schwankungsteil eines Schwingers. Die feste Hälfte ist das ½ in der Umrechnung in Masse. Die Regel zählt Rückläufe einer Schleife. Zählt man alle Vielfachen der Schleife zusammen, Z = ∏ 1/(1 − qⁿ), ist das bis auf q^{1/24} der Kehrwert der η-Funktion des Torus, und deren Spiegelung, η(−1/τ) = √(−iτ)·η(τ), bildet den Newton-Bereich auf den tiefen ab, s ↔ 4π²/s, mit dem Spiegelpunkt s = 2π, einem vollen Umlauf. Mit s = −ln|q| = 2π·Im τ liegt der Treffpunkt der Leiter, s = 1 bei r*, bei Im τ = 1/2π, einem Radiant. Der Spiegelpunkt bei Im τ = 1, dem vollen Umlauf: Die Schwelle zählt je Radiant, der Spiegel je Umlauf, wie beim Ganzen, dessen offener Teil ein Radiant eines Umlaufs ist. Tief ist ln Z = π²/(6s) + ½·ln(s/2π) − s/24: der führende Teil ist die erste Stufe der Leiter, π²/6. Bei der Windung, die sperrt, zählen alle Vielfachen ∏(1 + qⁿ), tief ln Z = π²/(12s) − ½·ln 2 + s/24, beide bis auf Glieder in e^{−4π²/s} bzw. e^{−2π²/s}. Die halbe Drehung halbiert den führenden Teil und kehrt den Nullpunkt um, +s/24 gegen −s/24: Wie Ladungen sich aufheben und Massen sich addieren, heben sich die Nullpunkte von Windung und Faltung auf. Der Spiegel bildet die Sperre auf die ungeraden Moden ab, ∏(1 + qⁿ) = (q·q̃)^{−1/24}·∏(1 − q̃^{2n−1})/√2 mit q̃ = e^{−2π²/s}. Der Spiegelpunkt der Windung liegt bei s = π√2, der der Faltung bei 2π. Windungen, die einander sperren, besetzen jede Mode höchstens einmal, ∏(1 + qⁿ). Rückläufe, die nur in ungeraden Moden kommen, zählen ∏1/(1 − q^{2n−1}). Beide sind gleich (Euler, bei q = ½ beide 2,384): Die Sperre der halben Drehung ist ein voller Umlauf auf dem doppelten Weg, der nach einem Umlauf das Vorzeichen wechselt.

**Entropie und Temperatur:** Die geometrische Verteilung der Rücklaufzahl hat bei festem Mittel die größte Entropie, H = ln ν + χ. Die Schwelle vergleicht die Unruh-Temperatur des ruhenden Beobachters am Ort mit der Periode der inneren Richtung: s² = 2π·T_U/T_I, kT_I = ħa₀/c. Eine Periode, keine Wärme.

**Kern-Nome und Bahntorus:** Der Bahntorus hat den Modul τ_M mit Re τ_M = ω_rad/ω_azim (Umlaufphase) und Im τ_M (Bindung). Für die Punktmasse ist (Re τ_M)² = 1 + χ: χ = 0 gibt geschlossene Keplerbahnen, χ = 1 Rosetten mit 254,6° zwischen den Perizentren. Die Kräfte r^a und r^A sind Spiegelpartner, wenn (a + 3)(A + 3) = 4: Die Abbildung w = z^{(a+3)/2} führt die Bahnen der einen in die der anderen über. Newton (r⁻²) und der gehaltene Körper (r¹) sind ein solches Paar, w = z². Der Spiegelpunkt der gebundenen Bahnen ist a = −1, die Kraft des tiefen Bereichs (der zweite, a = −5, führt in den Fall ins Zentrum): Die Regel läuft von Newton zum Treffpunkt der Kraftformen. Die Regel hebt die Kepler-Entartung κ = Ω auf und ersetzt sie durch örtliche Resonanzen bei rationalem κ/Ω. Der Kern-Nome ist |q_Δs| = Δs/r, die Spitzenform Δ = q_Δs·∏(1 − q_Δsⁿ)²⁴ verschwindet für q_Δs → 0. Die Bindung läuft nur über den Betrag des Nomes. Ginge die Phase mit ein, fiele die Verstärkung außen unter 1.

## Struktur statt Masse

Der Beobachter rechnet den Zusatz in eine Masse um, für die Punktmasse

```
M_schein/M = 1/(1 − e^{−s}) = r/r* + ½ + r*/(12r) − (r*/r)³/720 + …,     Φ_X = v_flach²·ln ν(s)
```

Die Korrekturen ½ und 1/12 sind Bernoulli-Zahlen: der Fehler der einzelnen Rückläufe gegen das Kontinuum, wie beim Vieleck gegen den Kreis. Das Zusatzpotential ist der Logarithmus der mittleren Zahl der Durchläufe und geht jenseits von r* in v_flach²·ln(r/r*) über. In der Masse gelesen wechselt der Anteil der Baryonen an der gelesenen Dichte bei der größten Zusatzbeschleunigung das Vorzeichen: Dünne Materie verstärkt sie, dichte mindert sie.

## Strukturwirkung und Vakuumterm

Die Wirkung der Metrik trägt dieselbe Funktion wie die Regel:

```
Q̃(s) = −4∫_s^∞ x³·n_B(x) dx,     Q̃(0) = −4π⁴/15
```

Die Konstante −4π⁴/15 ist ein globaler Vakuumterm: Er hängt weder an der Glättung noch an einzelnen Trichtern. Er ist der Vakuumteil der Dunklen Energie. Die x³ ist die Zustandsdichte der drei Raumdimensionen. Auf dem geschlossenen Torus steht die Konstante je Schleife einmal, also b₁ = 3-fach. Der örtliche Teil bleibt richtungsfrei.

**Ganz und örtlich gespiegelt:** Dieselbe Funktion hat zwei Lesarten. Über alle s integriert, für das Ganze, gibt jede Stufe Γ(k)ζ(k) mit k = d + 1. Am Ort, nahe dem tiefen Bereich entwickelt, hat sie die Beiwerte ζ(0) = −½, ζ(−1) = −1/12, ζ(−3) = 1/120: n_B(s) = 1/s + Σ ζ(−n)·(−s)ⁿ/n!. Beide verbindet der Spiegel ζ(1 − k) = 2·(2π)^{−k}·cos(πk/2)·Γ(k)ζ(k), mit (2π)^{−k}, der Zählung je Radiant. So ist die Konstante des Vakuumteils, Γ(4)ζ(4) = π⁴/15, das Spiegelbild von ζ(−3) = 1/120, das in der Umrechnung in Masse als −(r*/r)³/720 erscheint, und π²/6 das von −1/12: Was das Ganze trägt und was der Beobachter am Ort als Korrektur liest, ist dieselbe Größe von zwei Seiten.

## Metrik und Regel

Die Regel wirkt am Ort der Bahn auf alle Materie, auch ruhende. Die Metrik trägt nur ihren wirbelfreien Teil:

```
ν·g_N = −∇Φ + ∇×A,     ∇²Φ = ∇·[ν·∇Φ_N],     Ψ = Φ
```

Licht folgt der Metrik und sieht nur −∇Φ. Bahnen sehen auch den Wirbel ∇×A. Für Kugeln ist ∇×A = 0, dort sind Bahnen und Licht gleich. In Scheiben innen ist die dynamische Masse größer als die Linsenmasse. Prüfbar ist der Wirbel an ∂_z a_R ≠ ∂_R a_z.

Die Regel hat keine Wirkung (außer bei Kugeln): Ihre Kraft −ν∇Φ_N ist nur dann ein Gradient, wenn ν auf den Newton-Äquipotentialflächen konstant ist. Sie setzt am Erzeuger der Bewegung an (L_eff = 𝒫L), nicht an der Energie. Auf Sternbahnen mittelt sich der Arbeitsaustausch weg. Die einfachste Wirkung, die Feldform ∇·[ν∇Φ_N], setzt in einer Scheibe das senkrechte Scheibenfeld wie ein äußeres Feld zusammen. In einer früheren Anpassung an die SPARC-Rotationskurven (Υ_d = 0,5, ohne IC 2574) ging die Feldform nicht ein: Die Regel wirkt am Ort der Bahn.

Der Strukturterm ist eine Reihe in q mit |q| = e^{−s}, die für q → 0 in jeder Ordnung verschwindet. Er trägt zu keinem PPN-Koeffizienten bei, γ = β = 1, und Gravitationswellen laufen mit c. An der Erdbahn ist der Zusatz e^{−7480}.

## Hierarchie der Trichter

Kompakte Trichter, dichter als die Gezeitendichte ρ_J = 3𝒯/(4πG) ihrer Umgebung (𝒯 Gezeit des geglätteten Elternfelds), bleiben bis zu ihrer Glättungslänge ℓ, ℓ³ = GM/𝒯, Newtonsch. Verstärkt wird nur das kollektive, über diese Zonen geglättete Feld. Weil die Regel am Ort auf das geglättete Feld wirkt, in dem das eigene Feld am Schwerpunkt null ist, fallen Körper am selben Ort gleich.
