# Herleitung des Schalters χ = s/(e^s − 1)

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Schluss, „Offen ist, den diffusiven Rückweg aus der Wirkung abzuleiten“ → abgeleitet: Halbgerade aus (∇² + ∂_ϖ²)δσ = 0, Bose-Form aus der Wirkung (∂ℒ_𝒮/∂Y = 2n_B).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/schalter_chi.py` (Ausgabe `schalter_chi.txt`)

## Frage

Die Grenzwerte stammen vom Autor: χ = 0 für stabile KAM-Tori, χ = 1 für zerfallene. Dazwischen steht

```
a = a_N + χ·√(a_N a₀) = a_N·Σ_{n≥0} qⁿ = a_N/(1 − q),     q = e^{−s},     s = √(a_N/a₀),     χ = s·q/(1 − q)
```

Die Form hat drei Teile: die Summe über Rückläufe Σqⁿ, ein exponentielles Gesetz je Rücklauf und die Wurzel in s. Der Tiefgrenzfall legt fest, was verlangt ist: v⁴ = GMa₀ gilt genau dann, wenn 1 − q ≈ s für s → 0.

| Rückkehrgesetz | s·ν bei s = 0,01 | Tiefgrenzfall |
|---|---|---|
| q = e^{−s} | 1,005 | a → √(a_N a₀), Tully-Fisher |
| q = 1/cosh s | 200 | a → 2a₀ |
| q = e^{−s²} (linear in a_N) | 100 | a → a₀ |

## (a) Zerfall der KAM-Tori, Resonanzüberlappung

Modell ist die Standardabbildung p′ = p + K sin x, x′ = x + p′. Eine Resonanz hat die Halbbreite 2√K, der Abstand der Resonanzen ist 2π. Mit der Störstärke K ∝ ε = a₀/a_N (Schwelle gegen Newton-Bindung) wird Abstand/Breite ∝ √(a_N/a₀). **Die Wurzel entsteht hier von selbst**, aus der Pendelbreite ∝ √ε. Die Form entsteht nicht:

- **Anteil außerhalb rotierender KAM-Tori** (Inseln, Librationstest): 0,40·√K für K = 0,01…0,3, wie die Pendelfläche 16√K/4π² = 0,405·√K. Das ist ∝ 1/s, ein Potenzgesetz. Der Zusatz χ·√(a_N a₀) wäre im ganzen Newton-Bereich konstant, etwa 0,64·a₀ (Überlappung bei s = 1).
- **Chaotischer Anteil** (Lyapunov-Exponent), s = π/(2√K):

| s | 2,03 | 1,59 | 1,11 | 0,70 | 0,55 |
|---|---|---|---|---|---|
| chaotisch | 0,031 | 0,401 | 0,754 | 0,982 | 1,000 |
| χ(s) | 0,307 | 0,406 | 0,545 | 0,690 | 0,748 |

  Der Übergang ist viel steiler und bei endlichem s voll; die Werte treffen sich nur in einem Punkt. Für kleines K ist die chaotische Schicht ∝ exp(−π²/√K) = e^{−2πs}. Der Koeffizient 1 verlangte K(s = 1) = π⁴ ≈ 97, weit im voll chaotischen Bereich.
- **Aufspaltung der Separatrix** (Melnikov, auf 10 Stellen): 2π·sech(πΩ/2) bei Antrieb, πΩ²/(2 sinh(πΩ/2)) bei parametrischer Störung. Der Exponent ist stets πΩ/2, der Vorfaktor hängt von der Störung ab.

**(a) trägt nicht.** Es liefert die Wurzel, aber weder den Koeffizienten 1 noch die Bose-Form.

## (b) Besetzung einer Mode bei der Temperatur a₀

χ ist genau der Planck-Mittelwert einer Mode, geteilt durch den Gleichverteilungswert: χ = ⟨E⟩/kT mit s = ħω/kT. Entsprechend ist a = a_N(1 + n_B) der Bose-Verstärkungsfaktor. χ = 1 heißt Gleichverteilung (Tori zerfallen), χ = 0 eine eingefrorene Mode.

**Woher die Bose-Form kommt** (auf 12 Stellen): Ist ein Korrelator in imaginärer Zeit periodisch (Periode β) und liegt seine Singularität bei Koinzidenz, so antwortet eine Mode der Frequenz ω genau mit s/(e^s − 1), s = ωβ. Jede weitere Singularität im Abstand nβ ist ein Rücklauf mit dem Gewicht qⁿ. Beim Separatrix-Puls aus (a) liegen die Pole um einen halben Abstand neben der reellen Achse. Er gibt (s/2)/sinh(s/2); der Quotient zur Bose-Form ist genau e^{−s/2}.

**Die Wurzel entsteht nicht:**
- Unruh- und Horizonttemperatur (kT_GH = ħa₀/c) geben das lineare Verhältnis s′ = a_N/a₀. Dann gilt a → a₀ + a_N/2. Die Punktmasse 6·10¹⁰ M☉ liefe mit 340 km/s bei 30 kpc und 1054 km/s bei 300 kpc, statt flach bei 188 und 177 km/s.
- Energie = ħ·Bahnfrequenz, Temperatur kT_GH: s = 1818 bei 8 kpc, verlangt ist 1,04.
- s = √(a_N/a₀) verlangt eine Mode ω = √(a_N/ℓ₀), ℓ₀ = c²/a₀, mit einer Periode von 497 Mrd. Jahren bei a_N = a₀. Allgemein braucht es ω(a) = √(a/L) für Mode und Temperatur mit derselben Länge L. L kürzt sich heraus; dass beide dieselbe Länge haben, ist gesetzt.

**(b) trägt die Bose-Form und die Summe über Rückläufe, nicht die Wurzel.**

## (c) Rückkehr mit poissonschem Verlust

q ist die Wahrscheinlichkeit, einen Rücklauf ohne Flucht zu überstehen; der Verlust geschieht mit der Rate Γ.

- **Fester Umlauf T:** q = e^{−ΓT}. Die mittlere Zahl vollständiger Umläufe ist q/(1 − q) (Kac; Monte-Carlo 0,5831 gegen 0,5820). Mit Γ ∝ a_N folgt q = e^{−a_N/a₀}, der lineare Fall aus (b). Die Wurzel verlangt Γ = √(a_N/L) und T = √(L/a₀), wieder mit gesetzter gemeinsamer Länge.
- **Verlust längs des Weges nach außen** (Rate ∝ a_N): s = |Φ_N|/√(GMa₀). Für die Punktmasse ist das r*/r, aber nicht lokal: Im Plummer-Körper (b = 3 kpc) ist bei 0,5 kpc s = 2,75 statt √(a_N/a₀) = 1,11.
- **Diffusive Rückkehr, nach außen offen:** Ein Zustand im Abstand L von der Schleife diffundiert zurück. Vor dem Verlust kommt er mit q = E[e^{−ΓT}] = e^{−L√(Γ/D)} an. Die Wurzel entsteht von selbst. Ist die Verlustrate linear, Γ ∝ a_N (wie in (b)), und die Rückkehrzeit τ_D = L²/D ∝ 1/a₀, dann folgt

```
s = √(Γ·τ_D) = √(a_N/a₀),     q = e^{−s},     n Rückläufe: qⁿ,     ν = Σqⁿ = 1/(1 − q),     χ = s/(e^s − 1)
```

  a₀ ist die Beschleunigung, bei der Verlustrate und Rückkehrrate gleich sind. Es bleibt kein Zahlenfaktor frei. Irrfahrt mit 40 000 Läufern: q = 0,7388 / 0,3674 / 0,1349 bei s = 0,3 / 1 / 2, gegen e^{−s} = 0,7408 / 0,3679 / 0,1353. Vom doppelten Abstand aus ergibt sich q².
- **Geschlossener Kreis** (Rückkehr über ±L): q = 1/cosh s (0,6474 gegen 0,6481 bei s = 1). Dann 1 − q ≈ s²/2 und a → 2a₀, ohne Tully-Fisher. Der Rückweg muss einseitig offen sein; das passt zum rein radialen Kern e_r⊗e_r.

**Eindeutigkeit.** Der Verlust sei poissonsch mit Rate ∝ a_N^β, die Rückkehrdauer zufällig. Dann gilt q = e^{−√(a_N/a₀)} für alle a_N genau dann, wenn die Dauer einseitig stabil mit dem Index α = 1/(2β) verteilt ist (Laplace-Eindeutigkeit; Lévy-Transformierte auf 12 Stellen). Zwei Wege haben einen gewöhnlichen Mechanismus: feste Dauer mit Rate ∝ √a (eine Frequenz) oder diffusive Rückkehr mit Rate ∝ a.

**Zeitskala.** Festgelegt ist nur das Produkt Γ·τ_D. Mit Γ = a_N/c und τ_D = c/a₀ (79 Mrd. Jahre) dauern die überlebenden Rückkehren im Mittel τ_D/(2s), bei s = 1 also 40 Mrd. Jahre (π/H₀ für a₀ = cH₀/2π). Als Geschichte in Echtzeit ist das zu lang. Die Rechnung trägt als stationärer Ausbreitungskern: q = e^{−κL} löst G″ = κ²G mit κ²L² = a_N/a₀, eine abgeschirmte Ausbreitung zweiter Ordnung.

**(c) trägt in der diffusiven, einseitigen Form exakt, in der ballistischen und der geschlossenen nicht.**

## (d) Mittelung über Umläufe

- **Kern-Nome:** Δs/r, zeitgemittelt über Kepler-Ellipsen (Halbachse 8 kpc, e = 0; 0,5; 0,9), ist genau Δs/a = 1,20·10⁻⁷; verlangt ist q = 0,352. Ein Potenzgesetz bleibt nach jeder Mittelung ein Potenzgesetz dieser Größe. Mit Δ(τ_M) ≈ q_Δs gilt dasselbe.
- **Alle Harmonischen statt der einfachen Schleife:** Mit ν = 1 + Σ_n 1/(e^{ns} − 1) wird v⁴/(GMa₀) = 9,0 / 27 / 56 bei s = 0,1 / 0,01 / 0,001, mit Σσ₁(n)e^{−ns} wächst es wie 1/s². Nur die einfache Schleife (1,10 / 1,01 / 1,001) hält v⁴ = GMa₀.
- **Mittelung eines Verlustfaktors** e^{−(a_N/a₀)X} über Umlaufdauern X gibt e^{−√(a_N/a₀)} nur für die Lévy-Verteilung, also genau (c).
- **Frequenzdarstellung (exakt):** χ(s) = 1 − s/2 + 2Σ_{k≥1} s²/(s² + (2πk)²). Jeder Summand hat die Form des speichernden Anteils eines Maxwell-Elements, (ωt)²/(1 + (ωt)²) mit ωt = s/2πk. Das ordnet χ dem Gedächtniskern zu, legt s aber nicht fest.

**(d) trägt nicht.** Soweit eine Mittelung trägt, ist sie (c).

## Ergebnis

| Ansatz | Σqⁿ | Gesetz je Rücklauf | Wurzel | Tully-Fisher | trägt |
|---|---|---|---|---|---|
| (a) KAM, Chirikov | – | Inseln ∝ 1/s, Schicht e^{−2πs} | ja (Breite ∝ √ε) | ja (χ → 1) | nein |
| (b) Bose-Besetzung | ja | e^{−s} | nein (linear) | nein | nein |
| (c) fester Umlauf | ja | e^{−s} | nur mit gesetzter Länge | nur dann | nein |
| (c) geschlossener Kreis | ja | 1/cosh s | ja | nein | nein |
| (c) diffusiv, nach außen offen | ja | e^{−s} | ja | ja | ja, mit einer Annahme |
| (d) Mittelung | – | Potenzgesetz | nur über Lévy = (c) | nein | nein |

(b) liefert die Bose-Form mit dem linearen Verhältnis a_N/a₀ als Produkt aus Verlustrate und Rückkehrzeit. (c) macht daraus über die diffusive Rückkehr die Wurzel.

**Die fehlende Annahme:** Die Rückkehr in die Schleife ist ein diffusiver, nach außen offener Weg (zweiter Ordnung), dessen Zeitskala nicht von a_N abhängt; der Verlust geschieht mit einer zu a_N proportionalen Rate. a₀ ist die Beschleunigung, bei der beide Raten gleich sind. Mit dieser Annahme folgt χ = s/(e^s − 1) ohne freie Zahl; ohne sie ist die Bose-Form gesetzt. Der zweite exakte Weg (feste Dauer, Rate ∝ √a) braucht statt dessen eine gesetzte gemeinsame Länge für Mode und Schwelle. Der diffusive Rückweg ist aus der Wirkung abgeleitet: Die um die innere Richtung ϖ erweiterte Faltungsgleichung (∇² + ∂_ϖ²)δσ = 0 setzt das Feld harmonisch fort (e^{−|k|ϖ}), ϖ ≥ 0 ist eine Halbgerade mit dem Ort als unterster Ebene, und mit L = r* folgt genau e^{−s}; die Bose-Form folgt aus der Wirkung, ∂ℒ_𝒮/∂Y = 2n_B(s), also ν = 1 + n_B. [korr.]
