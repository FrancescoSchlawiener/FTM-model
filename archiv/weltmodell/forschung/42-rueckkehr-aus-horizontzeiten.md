# Eine Deutung von Ω_b/Ω_m über die zwei Zeiten des Horizonts

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: ν = 2π, „H kürzt sich heraus“, H/2π = a₀/c (Abschnitte 2–3) → trägt nicht: Zahlennähe ohne Herleitung, Q = 5,364 ist aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0); bei festem Λ wäre ν = 2ΛH/c zeitabhängig; die eigene Folge (c) liegt 11–14σ neben der Messung (5,46 gegen 6,364). λ_m = H (Abschnitt 4) → λ_m = ω_orb/N*, nur in gebundenen Bereichen. Auswahl von 2π aus 21 einfachen Zahlen als Zahlennähe gekennzeichnet. Die Notiz setzt a₀ = cH/2π. Umbenannt: (∂χ)² → (∂ξ)².

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/rueckkehr_horizont.py` (Ausgabe `rueckkehr_horizont.txt`)
> Messwerte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.

## Ausgangspunkt

Die Rückkopplungssumme ν = 1/(1−q) ist eine Kac-Rückkehrzeit. Kosmisch gilt 1 − q = Ω_b/Ω_m. Diese Notiz prüft eine Deutung der Hypothese ν = 2π über Zeiten des Modells.

**Vorab.** Im 2σ-Band der Messung liegen 21 einfache Zahlen (19/3, 25/4, 32/5, 3π−3, e²−1, …). 2π ist unter ihnen nur durch eine Deutung ausgezeichnet, und die folgende Deutung wurde gesucht, nachdem 2π als Hypothese feststand. Sie ist deshalb kein Beleg für 2π, sondern ein Angebot, das an eigenen Folgen geprüft werden muss (Abschnitt 3c). [korr.] Zahlennähe ohne Herleitung, trägt nicht.

**Lesart.** Es gibt nur baryonische Materie. Die Einschränkung der Freiheitsgrade (L_eff = 𝒫L, Projektor I − Π⊥) ändert die Regel ihrer Bewegung. ν ist die **Verstärkung** dieser Regel: Die Gravitation der Baryonen wirkt ν-fach. Ein Beobachter, der ungeänderte Regeln annimmt, liest die Verstärkung als zweite Dichte, Ω_c = (ν − 1)·Ω_b. „Ω_c“ ist in dieser Notiz nur der Name der Messgröße, keine Komponente.

## 1. Kac für Flüsse

Für einen stationären ergodischen Fluss im Phasenraum und eine Menge A gilt (Kac, Satz von Birkhoff):

```
μ(A) = E[Verweildauer in A] / E[Zyklus zwischen zwei Eintritten] ,        mittlere Rückkehrzeit = 1/μ(A)
```

Im Modell ist A der Teil des Zyklus, in dem die Freiheitsgrade der Materie **uneingeschränkt** sind. Außerhalb von A wirkt die Einschränkung 𝒫. Die Verstärkung ist die mittlere Rückkehrzeit nach A, gemessen in Verweildauern:

```
ν = Ω_m/Ω_b = 1/μ(A) = T_Zyklus / τ_frei
```

Das ist eine Aussage über die Zeitanteile einer Regel. Weder ein Objekt noch eine Dichte kommt dazu.

## 2. Die zwei Zeiten des Horizonts

[korr.] Diese Rechnung setzt ν = 2π und a₀ = cH/2π zu jeder Zeit voraus; das trägt nicht: Q = 2π − 1 bzw. ν = 2π ist Zahlennähe ohne Herleitung; Q = 5,364 ist aus der Messung eingesetzt, und im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0). Im Modell ist a₀ = c²/(2Λ) mit heute festem Λ; ≈ cH₀/2π gilt nur heute. Dann wäre ν = cH/a₀ = 2ΛH/c zeitabhängig. Die eigene Folge (c) liegt 11–14σ neben der Messung (5,46 gegen 6,364).

Der Hubble-Horizont hat zwei eigene Zeiten. Beide sind durch H festgelegt:

| Zeit | Wert | Herkunft |
|---|---|---|
| Zyklus | T = 2π/H | Euklidische Periode des Horizonts. Dieselbe Periode ergibt die Temperatur kT_GH = ħH/2π und damit die Schwelle a₀ = c·kT_GH/ħ = cH/2π. |
| Verweildauer | τ = 1/H | Relaxationszeit des Horizonts, **gewählt**: Die langsamste Abklingrate eines konform gekoppelten Feldes im statischen de-Sitter-Gebiet ist H [G]. Andere Felder klingen mit anderen Raten ab (Vielfache von H [G]); mit 2H oder 3H gäbe dieselbe Rechnung ν = 4π bzw. 6π. |

Damit folgt:

```
Ω_b/Ω_m = (1/H)/(2π/H) = 1/2π
ν = 2π = 6,28319,     q = 1 − 1/2π = 0,84085,     ν − 1 = 2π − 1 = 5,28319
```

**H kürzt sich heraus**, wenn beide Zeiten zu jeder Zeit an H hängen. Die Verstärkung bliebe dann konstant, während der Kosmos atmet; das verlangen die CMB und das Wachstum: Das gemessene Verhältnis Ω_c/Ω_b ist seit vor der Rekombination fest. Mit heute festem Λ kürzt sich H nicht heraus (ν = 2ΛH/c). [korr.]

## 3. Drei gleichwertige Formen

**(a) Schleife je Verweildauer.** Pro Relaxationszeit 1/H fällt der Anteil 1/2π des Zyklus in den freien Bereich, der Anteil q = 1 − 1/2π in den eingeschränkten. Die Rückkopplung Σ qⁿ · (1/H) = 2π/H summiert sich zu genau einem Zyklus.

**(b) Zyklusrate und Schließung.** Die Zyklusrate des Horizonts ist 1/T = H/2π. Sie ist gleich der kritischen Rate der Schließungsregel:

```
1/T = H/2π = a₀/c = θ_c
```

[korr.] H/2π = a₀/c gilt nur heute, im Rahmen a₀ ≈ cH₀/2π; im Modell ist a₀ = c²/(2Λ) mit festem Λ. Mit θ_c wird das Gewicht w = θ²/(θ²+θ_c²) gebildet. Unter der Annahme a₀ = cH/2π wären die Schwelle der Galaxien, die kritische Rate der Schließung und die Zyklusrate des Horizonts **eine** Größe und ν das Verhältnis aus Relaxationsrate und Zyklusrate: ν = H/θ_c. [korr.]

**(c) Schwelle und Baryonenanteil.**

```
a₀ = c·H·Ω_b/Ω_m
```

Das ist die schärfste eigene Folge der Deutung, und sie steht gegen die Messung. Mit dem gemessenen a₀ = 1,2·10⁻¹⁰ m/s² verlangt sie ν = cH₀/a₀ = 5,46 (H₀ = 67,4) bzw. 5,67 (H₀ = 70). Planck misst 6,364 ± 0,065, das sind 11–14 σ. Umgekehrt verlangt das Planck-Verhältnis a₀ = 1,03·10⁻¹⁰ (H₀ = 67,4) bzw. 1,07·10⁻¹⁰ m/s² (H₀ = 70). Das äußere Feld der Galaxien vergrößert den Abstand noch (bis etwa 40 %). Die Deutung steht und fällt damit mit a₀ = cH₀/2π.

## 4. Der kosmische Kern sitzt in der Resonanz

Auf der Kreisbahn wirkt der Gedächtniskern radial mit ω/(λ_m² + ω²). Bei fester Gedächtnisrate λ_m ist diese Antwort für ω = λ_m am größten, und die effektive Rate λ_eff = ω + λ_m²/ω ist dort am kleinsten (= 2ω). Für den Horizont gilt:
- euklidische Kreisfrequenz ω = H (Periode 2π/H),
- Relaxation λ_m = H.

[korr.] λ_m ist im Modell die Übergaberate der Mode ξ in die stationäre Regel, λ_m = ω_orb/N* mit ω_orb² = (4πG/3)ρ̄ und N* ≈ 0,3–2,4, und gilt nur in gebundenen Bereichen; λ_m = H für den Horizont ist nicht festgelegt. Unter dieser Setzung wäre das Verhältnis k = λ_m/ω = 1, das Maximum der Kernantwort. Das folgt aus der Wahl der Verweildauer 1/H und ist keine unabhängige Stütze.

## 5. Vergleich mit Messungen

Die Messgröße ist ν = (Ω_b h² + Ω_c h²)/Ω_b h².

| Datensatz [G] | ν | Abstand zu 2π | vorhergesagtes Ω_c h² = (2π−1)·Ω_b h² |
|---|---|---|---|
| Planck 2018 (TT,TE,EE+lowE+lensing) | 6,364 ± 0,065 | +1,26 σ | 0,1182 ± 0,0008 |
| Planck 2018 + BAO | 6,323 ± 0,053 | +0,75 σ | 0,1185 ± 0,0007 |
| BBN (D/H), Ω_b h² = 0,0222 ± 0,0006 | – | +0,90 σ gegen Planck-Ω_c h² | 0,1172 ± 0,0029 |

Die Abstände sind ohne die Korrelation zwischen Ω_b h² und Ω_c h² gerechnet. Die Rechnung mit dem Boltzmann-Code bei fester Winkelgröße des Schallhorizonts ergibt 1,0–1,5 σ.

## 6. Was hergeleitet ist und was gesetzt

- **Hergeleitet:** Nach Kacs Lemma ist die Verstärkung ν eine Rückkehrzeit, also Zyklus durch Verweildauer. Zwei Zeiten, die beide proportional zu 1/H sind, ergeben eine zeitlich konstante Verstärkung. Im Modell hängt die Schwelle aber an Λ, das heute fest ist, nicht an H; so wäre ν = 2ΛH/c zeitabhängig. [korr.]
- **Gesetzt:** Der Zyklus ist die euklidische Periode 2π/H, die Verweildauer ist 1/H. Die Verweildauer ist aus mehreren möglichen Abklingraten ausgewählt, und zwar die, die 2π ergibt.
- **Gegen die Messung:** Zusammen mit a₀ = cH₀/2π verlangt die Deutung ν = cH₀/a₀. Das sind 5,5–5,7 gegen 6,364 ± 0,065 gemessen.
- **Offen:** Die Nebenbedingung (∂ξ)² = −1 [korr.] mit der euklidischen Periodizität verknüpfen und die Verweildauer aus der Wirkung ableiten statt sie zu wählen.
- **Test von 2π selbst:** Die nächste Generation der CMB-Messungen (σ(ν)/ν ≈ 0,2–0,3 % [G]) trennt 2π von 6,36 mit 3–5 σ, sofern der heutige Mittelwert bleibt.

## Literatur

- Kac, M. (1947): On the notion of recurrence in discrete stochastic processes. Bull. AMS 53, 1002.
- Gibbons, G. W., Hawking, S. W. (1977): Cosmological event horizons, thermodynamics, and particle creation. Phys. Rev. D 15, 2738.
- López-Ortega, A. (2006): Quasinormal modes of D-dimensional de Sitter spacetime. Gen. Rel. Grav. 38, 1565.
