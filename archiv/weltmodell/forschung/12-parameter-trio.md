# Das Parameter-Trio (κ̃₀, λ_m, r*/3) und die Frage A gegen B

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: a₀ = cH₀/2π als „Atmen pro Umlaufwinkel“, 2π als Periode der Hubble-Schleife bzw. euklidische Periode 2π/H₀, r*/3 „aus H₀“ → a₀ = c²/(2Λ), ≈ cH₀/2π nur heute; τ_D ist die Laufzeit einer Windung, keine Periode in imaginärer Zeit; stationäre Regel als Exponentialkern κ̃₀e^{−λ_m t′} mit κ̃₀/λ_m = q → der Kern ist die Dichte erster Rückkehrzeiten ∝ Δt^{−3/2}, λ_m die Übergaberate von ξ; Nachlauf des Strukturfelds mit 1/λ_m → kein messbarer Nachlauf in Galaxien; λ_m-Kandidaten (H₀ usw.) → λ_m = ω_orb/N*; „nur eine Größe offen“ → offen sind N* und die Herkunft von Q, Λ ist heute fest; Umbenennung Integrationszeit s → t′; Im τ_M = r*/(2πr) als Punkt 1 des Autors gekennzeichnet. Abschnitt 2 rechnet mit a₀ = cH₀/2π = 1,082·10⁻¹⁰ (H₀ = 70), die Abschnitte 1 und 3 mit 1,2·10⁻¹⁰.

> Arbeitsstand, **noch nicht im Gesamtmodell**.
> Rechnung: `rechnungen/parameter_trio.py`

## 1. A und B sind dieselbe Formel

| | Formel | Milchstraße (6·10¹⁰ M☉) |
|---|---|---|
| A | Referenzlänge GM/(3v₀²), mit v₀ = 220 km/s | 1,778 kpc |
| B | Referenzlänge r*/3 = ⅓·√(GM/a₀), mit dem gemessenen a₀ = 1,2·10⁻¹⁰ m/s² | 2,783 kpc |

Aus der Feldgleichung folgt in Variante B v⁴ = G·M·a₀. Setzt man das in A ein:

```
GM/(3v⁴ᐟ²) = GM/(3√(GM·a₀)) = ⅓·√(GM/a₀) = r*/3
```

**A und B sind identisch.** Die eine universelle Formel lautet:

```
r*/3 = GM/(3·v_flach²)
```

Der Zahlenunterschied kommt allein vom Datenpaar (M, v₀), nicht von der Theorie:
- Zu M = 6·10¹⁰ M☉ gehört v_flach = 175,8 km/s.
- Zu v₀ = 220 km/s gehört M_bar = 1,47·10¹¹ M☉.

Welches Paar für die Milchstraße gilt, ist eine Frage der Messung ihrer baryonischen Masse und ihrer äußeren Rotationskurve. Mit dem Faltungsmesser Δs = v⁴/(3c²a₀) lässt sich das direkt prüfen.

## 2. Die Referenzlänge r*/3 aus Naturkonstanten und dem Atmen des Torus

[korr.] Dieser Abschnitt setzt die Schwelle als Atmen des Torus voraus, a₀ = cH/2π mit dem Hubble-Radius; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell ist a₀ = c²/(2Λ) mit der Schließlänge Λ der Rückkehr; Λ ist die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute. τ_D = ℓ₀/c = c/a₀ ist die Laufzeit einer Windung, keine Periode in imaginärer Zeit.

Die Schwelle ist die Beschleunigung des atmenden Ganzen pro Umlaufwinkel:

```
a₀ = c²/(2π·R_H) = c·H₀/2π,      R_H = c/H₀
```

Das 2π ist die Periode der Hubble-Schleife: a₀ = κ_dS/2π = c·f₁ mit κ_dS = c²/R_H und der Grundfrequenz f₁ = c/(2πR_H) = H₀/2π einer Lichtschleife vom Umfang 2πR_H (euklidische Periode 2π/H₀). [korr.] Das gilt nur heute (a₀ ≈ cH₀/2π); die Rückkehrzeit τ_D ist die Laufzeit einer Windung, keine Periode in imaginärer Zeit. Im Modell folgt a₀ = c²/(2Λ) aus der Kac-Rückkehrzeit auf der bei Λ geschlossenen Halbgeraden, heute Λ ≈ πc/H₀. Damit gilt heute:

```
r*/3 = ⅓·√(6π · Δs · R_H)
```

**r*/3 folgt aus G, c, a₀ = c²/(2Λ) (heute ≈ cH₀/2π) und der Masse, über Δs: r*² = 3Δs·ℓ₀ mit ℓ₀ = c²/a₀ = 2Λ. Es ist kein freier Parameter.** [korr.] ħ tritt auf galaktischen Skalen nicht auf. Die Skala kommt aus der globalen Topologie, nicht aus der Quantenwelt.

| M [M☉] | Δs [m] | r*/3 [kpc] | v_flach [km/s] |
|---|---|---|---|
| 10⁹ | 4,9·10¹¹ | 0,38 | 62 |
| 10¹⁰ | 4,9·10¹² | 1,20 | 110 |
| 6·10¹⁰ | 3,0·10¹³ | 2,93 | 171 |
| 1,47·10¹¹ | 7,2·10¹³ | 4,59 | 214 |
| 10¹² | 4,9·10¹⁴ | 11,96 | 346 |

Die Werte gelten mit a₀ = cH₀/2π = 1,082·10⁻¹⁰ m/s² bei H₀ = 70 km/s/Mpc.

## 3. κ̃₀ und λ_m: Nur ihr Verhältnis ist festgelegt

[korr.] Dieser Abschnitt setzt für die stationäre Regel einen Exponentialkern κ̃₀e^{−λ_m t′} voraus; das trägt nicht: Als Gedächtniskern ist die stationäre Regel die Dichte erster Rückkehrzeiten (∝ Δt^{−3/2}), der Zwanzig-Mori-Kern einer offenen Diffusionskette; ein Maxwell-Element (Exponentialkern) oder endlich viele Stufen geben 1 − q ∝ Γ ohne Wurzel und damit kein Tully-Fisher. λ_m ist die Übergaberate der Mode ξ in die stationäre Regel (λ_m = ω_orb/N*), nicht die Zerfallsrate der stationären Regel.

Der Gedächtniskern hat hier die Beschleunigungsform κ̃₀·e^{−λ_m t′} mit der Kernstärke κ̃₀ in s⁻¹. Sie ist von der Kernstärke κ₀ (in s⁻²) der Geschwindigkeitsform zu unterscheiden. Als stationäres Feldgedächtnis gälte:

```
g = g_N + ∫₀^∞ κ̃₀ e^{−λ_m t′} g dt′   ⟹   g = g_N / (1 − κ̃₀/λ_m)
```

Der Vergleich mit der Rückkopplungssumme g = g_N/(1−q) ergibt:

```
κ̃₀/λ_m = q(r) = exp(−r*/r),      Im τ_M = r*/(2πr)
```

[korr.] Im τ_M = r*/(2πr) ist Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise.

Das Verhältnis κ̃₀/λ_m wäre unter dem Exponentialkern also **kein Parameter, sondern die Topologie selbst**, nämlich der Nome des Torus an der Stelle r (für die Punktmasse der Rückkehranteil q). [korr.] Für die Milchstraße:

| r [kpc] | 2 | 5 | 8,35 (r*) | 20 | 50 |
|---|---|---|---|---|---|
| κ̃₀/λ_m = q | 0,015 | 0,188 | 0,368 | 0,659 | 0,846 |

**λ_m allein** tritt in allen stationären Systemen nicht auf. Es wirkt nur, wenn sich die Massenverteilung ändert. [korr.] In Galaxien läuft dabei nichts nach: Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Geschwindigkeitskern würde exzentrische Bahnen runden, alte Scheiben verlangen dafür λ_m ≪ 0,01·Ω. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche). Frühere Kandidaten:

| λ_m | Nachlaufzeit |
|---|---|
| H₀ (Atmen) | 14 Mrd. Jahre |
| 1/(1 Mrd. a) | 1 Mrd. Jahre |
| 1/(100 Mio. a) | 0,1 Mrd. Jahre |

[korr.] Keiner davon gilt: Im Modell ist λ_m = ω_orb/N* mit ω_orb² = (4πG/3)ρ̄ und N* ≈ 0,3–2,4 (aus dem Bullet Cluster); bei der 200–500-fachen kritischen Dichte des gebundenen Bereichs ist 1/λ_m heute 0,28–3,5 Mrd. Jahre.

**Wie λ_m bestimmt wird:** Verschmelzende Galaxienhaufen sind das natürliche Messinstrument. Der Versatz zwischen Linsenmaximum und Gas hängt von der Zeit seit dem Stoß im Verhältnis zu 1/λ_m ab. Ein Beispiel ist der Bullet Cluster.

## Ergebnis

| Parameter | Herleitung | frei? |
|---|---|---|
| r*/3 | ⅓√(3Δs·ℓ₀), ℓ₀ = c²/a₀ = 2Λ, heute ⅓√(6π·Δs·R_H) [korr.] | nein |
| a₀ | c²/(2Λ), Schließlänge Λ der Rückkehr, heute ≈ c·H₀/2π [korr.] | nein, Λ heute fest |
| κ̃₀/λ_m | q = exp(−r*/r), nur unter einem Exponentialkern; die stationäre Regel ist die Dichte erster Rückkehrzeiten [korr.] | nein |
| λ_m | Übergaberate der Mode ξ, λ_m = ω_orb/N*; Nachlauf nur, wo ξ nicht ganz übergegangen ist [korr.] | N* messbar an Haufenstößen und Linsen um isolierte Galaxien |

[korr.] Offen sind N* (damit 1/λ_m) und die Herkunft von Q = 5,36; Λ ist heute fest.
