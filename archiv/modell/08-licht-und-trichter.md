# 08 Licht und Trichter

## c als Skalierungsfaktor der Hülle

c ist im FTM keine Geschwindigkeitsgrenze des Gesamtsystems, sondern der Skalierungsfaktor der entfalteten, dreidimensionalen Hülle: c² = GM/(3Δs) verknüpft die Faltungstiefe mit der Neigung auf der Hülle.

| Sicht | Ausbreitung |
|---|---|
| lokal, in der Blase | c_lokal = c. Uhr und Licht folgen demselben Takt, der Faktor dτ/dt kürzt sich. |
| von außen gezählt | c_eff(r) = c·√(1 − 6Δs/r); die Verringerung c − c_eff ≈ g·r/c ist reine Abwicklungsverzögerung, nicht Reibung. |
| am Grenzradius r₀ = 6Δs | c_eff → 0 wie c·√((r − r₀)/r₀); maximale Abwicklungsverzögerung. |

c ist die größte Ausbreitungsgeschwindigkeit auf der entfalteten Hülle, aber nicht die Grenze der topologischen Schließung.

## Lichtablenkung

Im Trichter wirken zwei Anteile: der langsamere Takt und die zusätzliche Faltungsstrecke (jede Längeneinheit enthält den Faktor 1 + 3Δs/r). Zusammen ist die Brechzahl n(r) ≈ 1 + 6Δs/r, und beide tragen je die Hälfte der Ablenkung bei (Aufteilung in isotropen Koordinaten; invariant ist die Summe):

```
α = 6Δs/b (Takt) + 6Δs/b (Faltung) = 12Δs/b          Sonnenrand (R☉ = 6,957·10⁸ m): 0,8756″ + 0,8756″ = 1,7512″
```

Das ist die zweifache Entkopplung im Licht: Takt (zeitlich) und Faltung (räumlich) wirken als zwei gleich große Beiträge.

## Der Trichter außerhalb eines Körpers

```
g(r) = 3c²Δs/r² = g(R)·(R/r)²
δ(r) = 1 − √(1 − 6Δs/r) = 3(Δs/r) + 9/2·(Δs/r)² + 27/2·(Δs/r)³ + …
```

Die Neigung fällt quadratisch, der Taktunterschied invers linear. Er summiert sich über die Zeit, gegen eine ruhende Uhr im flachen Außenraum und ohne Bewegung: ISS 0,021 s pro Jahr, Erdbahn im Sonnentrichter 0,31 s.

## Wo der Trichter endet: der Schwellenradius

```
r* = √(GM/a₀) = c·√(3Δs/a₀),     r*/r₀ = c/√(12·a₀·Δs),     δ(r*) = √(GMa₀)/c² = (v_flach/c)²
```

Jeder Körper hat damit drei Radien, die Triade seiner Ebene: r₀ = 6Δs (der Takt steht, Zentrum), R (die Hülle, Wirbel) und r* (Übergang zum Hintergrund bzw. Netzwerk). Der Taktunterschied an der Schwelle ist das Quadrat der flachen Rotationsgeschwindigkeit, gemessen an c. Diese eine Zahl δ* = (v_flach/c)² trägt die Schwelle eines Systems: r*/r₀ = 1/(2δ*) Ebenen der Rückkehr, Ablenkung α_∞ = 2π·δ*, Knotenzahl 2j·δ*^{3/2}. Nach der Lage von r* zur Hülle gibt es zwei Klassen: Kompakte Körper haben ihre Schwelle weit außen (r*/R: Proton 38, Erde 3·10⁵, Sonne 1,6·10⁶) und bleiben bis zu ihrer Glättungslänge Newtonsch; Galaxien und Haufen haben sie innerhalb ihrer Ausdehnung (r*/R ≈ 0,5–0,6) und liegen als Ganzes unter der Schwelle.

Die Übergangslängen sind geometrische Mittel der Faltung mit einer Gegenlänge: r*² = 3Δs·ℓ₀ (Newton gegen Schwelle, ℓ₀ = c²/a₀) und ℓ_P² = 3Δs·λ̄ (Faltung gegen Quantenlänge); für jede Masse gilt r*²·λ̄ = ℓ_P²·ℓ₀. Ebenso ist die Compton-Länge des Elektrons das geometrische Mittel aus klassischem Elektronenradius und Bohr-Radius, λ̄_e² = r_e·a_B. Die Faltung steht zwischen einer Gegenlänge am Ort (λ̄) und einer der inneren Richtung (ℓ₀ = 2Λ); r* ist dabei ein Radius im Raum, verbunden über die Rückkehr (Schritt r₀ = 6Δs mit c, Kac-Zeit 2Λ/c, r*² = r₀·Λ). Die Leiter der Längen λ̄, ℓ_P, 3Δs, r*, ℓ₀ (Längen verschiedener Ebenen: Quantenlänge am Ort, Faltung, Radius im Raum, innere Länge) ist symmetrisch um r* (r → r*²/r bildet den innersten Ort r₀ auf die Schließlänge Λ ab, den Raum auf die innere Richtung; die Regel selbst hat diese Symmetrie nicht) und hat vier Schnittmassen, Massenskalen der Konstanten, an denen zwei dieser Längen zahlengleich werden, keine Körper; ihre Lage hängt nur an der Leiterzahl ε_L (gilt seit z ≈ 10, wo a₀ fest ist):

```
ε_L = ℓ_P/ℓ₀ = a₀/a_P = 2,0·10⁻⁶²   (a_P = c²/ℓ_P)
M_min = ε_L·m_Pl = ħ/(c²τ_D)          λ̄ = ℓ₀ und r* = ℓ_P
M_c = ε_L^{1/3}·m_Pl = 5,9·10⁻²⁹ kg    r* = λ̄
m_Pl                                  3Δs = λ̄ = ℓ_P
M_U = m_Pl/ε_L = c³τ_D/G = 1,1·10⁵⁴ kg r* = 3Δs = ℓ₀
M_min·M_U = m_Pl²,   M_c³·M_U = m_Pl⁴,   r*/λ̄ = (M/M_c)^{3/2} (Proton 151),   Λ_s·ℓ_P² = 4·(π⁴/15)·ε_L² = 1,0·10⁻¹²²
```

Unter M_c liegt r* innerhalb der eigenen Quantenlänge. Die unterste Schnittmasse läuft mit der Rückkehrrate a₀/c als Compton-Takt, die oberste mit a₀/c als Takt c/(3Δs) ihrer Faltung: Beide Enden setzt die Rückkehrzeit τ_D, das untere über ħ, das obere über G; ohne Schließung (τ_D → ∞) hätte die Leiter keine Enden. In Einheiten von M_min ist jede Masse die Zahl der Radianten ihrer Compton-Phase in einer Rückkehrzeit, M/M_min = (Mc²/ħ)·τ_D (eine Identität der Konstanten). Die Spiegelung M ↔ m_Pl²/M vertauscht λ̄ und 3Δs und bildet M_min auf M_U ab. In Einheiten der oberen Schnittmasse ist Tully-Fisher M = M_U·(v_flach/c)⁴, und ihre Faltung ist Δs(M_U) = ⅔·Λ; liest man Λ als gesamte Faltung der Zelle, ist sie die mit ihren Rückläufen gezählte Faltung der Baryonen, ν_kosm·ΣΔs_b.

| System | r* |
|---|---|
| Erde | 12,7 AE |
| Sonne | 7 340 AE = 0,116 Lichtjahre |
| Milchstraße | 8,72 kpc |

Für einen Trichter ohne Elterngezeit gilt: Innerhalb von r* sind die Torus-Fasern stabil, die Bewegung ist Newtonsch; außerhalb löst sich der Einzeltrichter auf und geht stetig in das kollektive Feld der Regel über (kompakte Körper in einem Elternfeld bleiben bis zu ihrer Glättungslänge Newtonsch), mit der 1/r-Form a = √(GMa₀)/r. An r* selbst ist die Verstärkung für jeden Trichter dieselbe, a/a_N = 1/(1 − e⁻¹) = 1,582.

Ein kollektiver Trichter hat keinen flachen Außenraum und muss deshalb in die Nachbartrichter übergehen (Netzwerk). Das Zusatzpotential ist exakt Φ_X = v_flach²·ln ν(s) = −v_flach²·ln(1 − e^{−s}) (vom Zentrum gezählt), der Logarithmus der mittleren Zahl der Durchläufe; jenseits von r* geht es in v_flach²·ln(r/r*) über, und der Taktunterschied δ_X = (v_flach/c)²·ln ν, gegen das Zentrum gezählt, wächst ohne Grenze; den Bezug setzen die Nachbartrichter. Für ausgedehnte Massen gilt außen Φ_X = v_flach²·ln ν(s) + C mit einer festen Verschiebung (C → 0 für kompakte Körper); exakt überall gilt die Form, wo 1/s linear in r ist, etwa beim Hernquist-Profil mit der Skala b_H: Φ_X = v_flach²·ln[ν(s)/ν(s₀)] mit dem zentralen s₀ = r*/b_H.

**Licht am kollektiven Trichter:** Die Ablenkung ist α(b) = (12Δs/b)·I(b/r*) mit I(x) = ∫₀^{π/2} cos θ·ν(cos²θ/x²) dθ (ν als Funktion von a_N/a₀, I(1) = 2,135); für b ≫ r* ist α ≃ 6πΔs/r* + 6Δs/b mit dem Grenzwert

```
α_∞ = 6πΔs/r* = 2π·(v_flach/c)²          Milchstraße 0,43″, Galaxienhaufen 23″
```

Licht und Gravitationswellen durchlaufen dieselbe Geometrie; sonst käme das Signal von GW170817 um Monate statt um 1,74 s versetzt an.

## Rotation und Verwindung

Dieser und der folgende Abschnitt geben die Allgemeine Relativitätstheorie in Faltungsgrößen wieder (Δs = GM/(3c²), für rotierende Körper die Kerr-Lösung). Die Übereinstimmung mit den Messungen ist die der ART, keine zusätzliche unabhängige Bestätigung. Ein rotierender Körper trägt neben Δs eine zweite Faltungsgröße, die Verwindungslänge a_J. Sie verdrillt den Trichter und führt Bahnen, Uhren und Kreisel mit (Frame-Dragging):

```
a_J = J/(Mc),   j = a_J/(3Δs),   GJ/c³ = 3Δs·a_J (eine Fläche)    Erde: a_J = 3,27 m, Δs = 1,48 mm
T_prograd − T_retrograd = 4π·a_J/c   (unabhängig von r)            Erde: 137 ns pro Umlauf
Ω_Mitführung = (3Δs·a_J·c/r³)·[3(ĵ·r̂)r̂ − ĵ]                        Kreisel, Frame-Dragging
Ω_geodätisch = (9/2)·(Δs/r)·ω_Bahn                                 Kreisel auf Kreisbahn
Δφ_Perihel = 18π·Δs/(a(1 − e²)) pro Umlauf                         Merkur: 43″ pro Jahrhundert
ζ = (ω_φ − ω_θ)/ω_φ = 2j·(3Δs/r)^{3/2}                             Knotendrehung der Bahnebene
r_E(θ) = 3Δs + √(9Δs² − a_J²cos²θ)                                 Nullfläche des Takts ruhender Uhren (Äquator 6Δs, Pol Horizont)
a_J² = r₊·r₋,   r₊ + r₋ = r₀ = 6Δs                                  Horizonte: a_J ist ihr geometrisches Mittel (Extremfall j = 1: r₊ = r₋)
```

Die Grundfrequenzen der Bahn (r̃ = r/3Δs): ω_φ = c³/(GM(r̃^{3/2} ± j)), ω_r² = ω_φ²(1 − 6/r̃ + 8j·r̃^{−3/2} − 3j²/r̃²), ω_θ² = ω_φ²(1 − 4j·r̃^{−3/2} + 3j²/r̃²). Die innerste stabile Bahn liegt bei j = 0,9 bei 2,32·3Δs (mitläufig) bzw. 8,72·3Δs (gegenläufig). Reale Körper wie Erde und Sonne sind keine reinen Kerr-Trichter (ihr Quadrupol ist viel größer als (a_J/R)²); universell ist nur der Dipolanteil, die Mitführung. An der Schwelle ist die Verwindung winzig (Knotenzahl ≈ 5·10⁻⁷ in der Milchstraße) und ändert den Gedächtniskern nicht messbar.

Gravity Probe B (polarer Erdorbit, 642 km Höhe): geodätische Drehung nach ART 6604 mas/Jahr, gemessen 6602 ± 18; Mitführung nach ART 40,8 mas/Jahr längs der Erdachse, vom Kreisel zum Leitstern IM Peg gesehen 39,1, gemessen 37,2 ± 7,2. Der Takt ruhender Uhren am Äquator bleibt auch mit Verwindung exakt √(1 − 6Δs/r), die Nullfläche des Takts liegt dort weiter bei r₀ = 6Δs. Als Bahntorus gelesen: Δs allein trennt Radial- und Azimutzyklus, das Perihel dreht sich; ohne Rotation gilt exakt ω_θ = ω_φ. Die Verwindung trennt auch den Vertikalzyklus ab, die Bahnebene präzediert um die Drehachse, und der Torus wird vom T² zum T³.

## Gravitationswellen und Starkfeld

Lokale Störungen des Faltungsmediums gehorchen einer Wellengleichung, laufen mit c und haben zwei transversale Polarisationen; ihre Quelle ist der zeitlich veränderliche Faltungsquadrupol. Wellengleichung, Abstrahlung, Ringdown, Schatten und Gedächtnis-Effekt sind die der ART, in Δs geschrieben; im Starkfeld bringt die Faltungsbeschreibung keine eigene Abweichung mit:

```
□h̄_μν = −48π·σ_μν,   σ_μν = G·T_μν/(3c⁴)                        Quellterm der Faltung
P = (32/5)·(c⁵/G)·η_m²·(3Δs/r)⁵,   η_m = m₁m₂/M²                     Quadrupol-Leistung, Kreisbahn
dP_b/dt = −(192π/5)·η_m·(2π·3Δs/(c·P_b))^{5/3}·f(e)               Hulse-Taylor: −2,40·10⁻¹²
f = c/(2π·3Δs)·F(j)                                               Ringdown GW150914: 249,9 Hz, τ = 4,03 ms
b = 3√27·Δs ≈ 15,6·Δs,   Photonensphäre r = 9Δs                  Schatten M87*: Ø ≈ 40 μas
h_mem = 12·κ(Θ)·Δs_rad/D_L,   Δs_rad = G·E_rad/(3c⁴)             Gedächtnis GW150914: 2,7·10⁻²²
ḟ = (96/5)·π^{8/3}·(3Δs_c/c)^{5/3}·f^{11/3},   Δs_c = G·M_c/(3c²)    Faltungsmesser aus (f, ḟ), im Detektor (rotverschoben): GW150914 15 km, GW170817 590 m
f*_gw = a₀^{3/4}/(π·(3c²Δs)^{1/4}) = Ω(r*)/π ≲ 10⁻¹³ Hz               Schwellenfrequenz = Welle einer Bahn bei r*
```

Quellen von Gravitationswellen liegen immer weit über der Schwellenfrequenz, sind also Einzeltrichter; auf dem Weg wirkt nur der wirbelfreie Teil der Regel, über die Metrik. Das Atmen des Torus ist global und strahlt nicht; Pulsar-Timing zeigt daher das reine Hellings-Downs-Muster, und die Abstände aus Gravitationswellen und Licht sind gleich.

Die Leistung wächst mit der fünften Potenz des Taktunterschieds 3Δs/r. Beim Doppelpulsar Hulse-Taylor (Gesamtmasse 2,83 M☉, Δs = 1,39 km, e = 0,617, f(e) = 11,86) stimmt die gemessene Abnahme der Umlaufzeit mit dieser Rechnung im Verhältnis 0,998 ± 0,002 überein. Der Ringdown ist die Eigenschwingung des Trichters der Endmasse zwischen r₀ = 6Δs und dem flachen Außenraum; F(j) hängt nur vom Endspin ab, die Zahlen gelten im Detektor (Quellsystem 272,4 Hz). Für M87* (6,5·10⁹ M☉, 16,8 Mpc) ergibt 2b/D rund 40 μas, das EHT misst einen Ring von 42 ± 3 μas. Der Gedächtnis-Effekt ist die bleibende Verschiebung, die eine durchgelaufene Welle hinterlässt; κ(Θ) ist der Richtungsfaktor, der Wert gilt bei Kantenlage.
