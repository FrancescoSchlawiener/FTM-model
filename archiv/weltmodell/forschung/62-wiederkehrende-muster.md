# Wiederkehrende Muster und die Rolle der Dimension

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Phantomdichte“ → gelesene Dichte, 4πGρ_X = χΩ² außerhalb der Masse, in der Masse + 4πGρ_b·∂(a − a_N)/∂a_N (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt; kein Stoff); „scheinbare Masse“ → gelesene Masse; „Nullpunkt einer Schleife (1/12 = −ζ(−1)) eindimensional“ als Zahlennähe gekennzeichnet; Skala der Profile L → R_s.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/wiederkehrende_muster.py`, `rechnungen/stufe_zwei.py`, `rechnungen/ausgedehnte_massen.py` (Ausgaben `.txt`)
> s = √(a_N/a₀) (Punktmasse r*/r), n_B = 1/(e^s − 1). Aufgenommen ist nur Exaktes und Nachgerechnetes.

## 1  Die Bose-Leiter

Die Hauptgrößen der Regel sind eine Leiter derselben Funktion, jede Stufe die vorige mal s:

| d | s^d·n_B(s) | Größe im Modell | ∫₀^∞ = Γ(d+1)ζ(d+1) | Größe im Modell |
|---|---|---|---|---|
| 0 | n_B | ν − 1 = M_schein/M − 1 (Summe der Rückläufe) | divergiert (1/s) | Tiefbereich |
| 1 | s·n_B | χ, Verlust an Flussdimension (d_eff = 3 − χ), gelesene Dichte 4πGρ_X = χΩ² (außerhalb der Masse) [korr.] | π²/6 = 1,6449 | ∫(a − a_N) d ln r = (π²/6)a₀ je Punktmasse; Boden (π²/6)a₀ über alle Harmonischen |
| 2 | s²·n_B | (a − a_N)/a₀, Zusatzbeschleunigung | 2ζ(3) = 2,4041 | ∫(Ω² − Ω_N²) d ln r = 2ζ(3)·Ω_N(r*)² je Punktmasse |
| 3 | s³·n_B | ¼·dQ̃/ds, Strukturwirkung | π⁴/15 = 6,4939 | Vakuumteil −Q̃(0) = 4π⁴/15 |

- Bei r* (s = 1) sind alle Stufen gleich: ν − 1 = χ = (a − a_N)/a₀ = ¼·dQ̃/ds = 0,582.
- Die Form s^d·n_B ist das Planck-Spektrum in d Dimensionen (Zustandsdichte s^{d−1}, Energie s). Die Wien-Stellen liegen bei s → 0 (d = 1), s = 1,594 (d = 2, a_N = 2,54a₀, Wert 0,648) und s = 2,821 (d = 3, a_N = 7,96a₀).
- **Maß der Stufe:** Jede Stufe gehört zum Maß ihres Flussgesetzes, ∫(a − a_N)·(r*/r)^{d−1} d ln r = a₀·Γ(d+1)ζ(d+1) (auf 10 Stellen): d = 1 mit d ln r (Skala), d = 2 mit d(1/r) (Umfang, Krümmung der Bahnkreise), d = 3 mit d(1/r²) (Fläche, Newton-Feld).
- **Die Größe der Stufe 2:** Mit (a − a_N)·r*/r = r*(Ω² − Ω_N²) wird daraus ∫(Ω² − Ω_N²) d ln r = 2ζ(3)·Ω_N(r*)²: Die zusätzliche Umlauffrequenz im Quadrat, über den Logarithmus des Radius summiert, ist 2ζ(3) mal das Quadrat der Kepler-Frequenz am Schwellenradius (Milchstraße: Ω_N(r*) = 6,8·10⁻¹⁶ s⁻¹, Umlauf 292 Mio. Jahre).
- **Zusatzpotential:** Φ_X(r) = ∫₀^r(a − a_N)dr = v_f²·ln ν(s) = −v_f²·ln(1 − e^{−s}) exakt (vom Zentrum gezählt; auf 10 Stellen bei 0,3 bis 30 r*). ln ν ist der Logarithmus der Zustandssumme einer Bose-Mode; außen v_f²·ln(r/r*), Takt δ_X = (v_f/c)²·ln ν. Dieselbe Zahl als Potential: ∫Φ_X d(a_N/a₀) = 2ζ(3)·v_f².

## 2  Dimension in den Zahlen

**Die 3 in Δs ist die Raumdimension.** Der Überschussradius der homogenen d-Kugel ist ∫Gm(r)/(rc²)dr = GM/(d·c²) (nachgerechnet d = 1 … 5), und 1/d = V/(S·R):

```
Δs = (GM/c²) · V/(S·R)         Faltungstiefe = Schwarzschild-Länge mal Volumen durch (Oberfläche × Radius)
```

Damit sind die ganzen Zahlen des Modells Vielfache der Dimension, die der ART einfache Zahlen: Takt 3Δs/r (ART 1), Grenzradius 6Δs = 2d (ART 2), Photonensphäre 9Δs = 3d (ART 3), Ablenkung 12Δs/b = 4d (ART 4).

**Flussdimension.** Newton ist Fluss durch die Kugelfläche (d = 3), der Tiefbereich Fluss durch einen Umfang (d = 2), eine bei Λ geschlossene Rückkehrrichtung gäbe d = 1 (konstante Beschleunigung). Die Planck-Leiter folgt demselben Schritt: Zusatzbeschleunigung zweidimensional, Strukturwirkung und Vakuumteil dreidimensional, Summe über alle Harmonischen und Nullpunkt einer Schleife (1/12 = −ζ(−1)) eindimensional. [korr.] Zahlennähe ohne Herleitung, trägt nicht (für den Nullpunkt einer Schleife).

## 3  Übergangslängen als geometrische Mittel

```
ℓ_P² = 3Δs·λ̄  (Faltung = Quantenlänge)        r*² = 3Δs·ℓ₀  (Newton = Schwelle)        r*²·λ̄ = ℓ_P²·ℓ₀  für jede Masse
```

Nachgerechnet für Proton, Sonne und Milchstraße: ℓ_P = 1,6162·10⁻³⁵ m fest, r* ∝ √M. Beide Übergänge haben dieselbe Form, eine Faltung mal eine Gegenlänge; klein ist die Gegenlänge die Quantenlänge, groß die Schwellenlänge.

## 4  Ausgedehnte Massen

Gerechnet für Plummer (Kern), Hernquist (Spitze), homogene Kugel und dünne Exponentialscheibe (Mittelebene, Regel algebraisch), Skala R_s = 0,03 bis 3 r*.

**Ortsmaß.** Die Leiter ∫(a − a_N)(r*/r)^{d−1} d ln r ist exakt nur für die Punktmasse. Für d = 1 hängt sie vom Profil ab (Kern: Faktor 0,5 bis 3,2) und divergiert logarithmisch, wo a − a_N im Zentrum endlich bleibt (Hernquist). Für d = 2 und 3 divergiert sie bei jedem ausgedehnten Profil, weil das Gewicht 1/r bzw. 1/r² im Zentrum wächst und a − a_N dort nicht exponentiell verschwindet.

**Maß der Bose-Variablen.** Längs des Radius liefert jeder monotone Ast von s(r) genau ∫₀^{s_max} s^d n_B ds, für jedes Profil (auf 5 Stellen: Plummer, Hernquist, Kugel, Scheibe). Ein Profil mit Kern (a_N(0) = 0) durchläuft s zweimal, hinein und hinaus; eine Spitze einmal; die Punktmasse einmal bis s_max = ∞ und gibt so die vollen Werte Γ(d+1)ζ(d+1). Die Leiter ist also eine Eigenschaft des Zustandsraums (der Variablen s), nicht des Ortes; s_max = √(a_N,max/a₀) misst die Kompaktheit.

**Zusatzpotential.** Außerhalb der Masse gilt Φ_X = v_f²·ln ν(s) + C mit einer festen Verschiebung C (Plummer +0,145 / −0,21 / −1,01 v_f² bei R_s = 0,3 / 1 / 3 r*; C → 0 für kompakte Profile). Exakt gilt die geschlossene Form überall, wo 1/s linear in r ist; das trifft das Hernquist-Profil (a_N = GM/(r + R_s)², eine verschobene Punktmasse): Φ_X = v_f²·ln[ν(s)/ν(s₀)], s₀ = r*/R_s, nachgerechnet C = −ln ν(s₀) = −0,0363 / −0,4587 / −1,2607.

**Gelesene Dichte in der Masse** [korr.] (Lesung eines Beobachters, der ungeänderte Regeln annimmt; kugelsymmetrisch, exakt, nachgerechnet am Plummer-Profil auf 6 Stellen):

```
4πGρ_X = χΩ² + 4πGρ_b·∂(a − a_N)/∂a_N,     ∂(a − a_N)/∂a_N = n_B(1 − sν/2)
```

Außerhalb der Masse bleibt χΩ². In der Masse kommt ein Anteil proportional zur Baryonendichte hinzu, der genau bei sν = 2, also 2(1 − e^{−s}) = s, das Vorzeichen wechselt: bei a_N = 2,54·a₀, der Wien-Stelle der Stufe 2 (größte Zusatzbeschleunigung). Darunter verstärkt dünne Materie die gelesene Dichte, darüber mindert dichte Materie sie.

## Ergebnis

Drei Muster kehren überall wieder: die Bose-Leiter s^d·n_B mit den Konstanten Γ(d+1)ζ(d+1), die Dimension als Faktor (die 3 in Δs, die Flussdimension 3 → 2 → 1, die Planck-Dimension jeder Größe) und die Übergangslänge als geometrisches Mittel einer Faltung mit einer Gegenlänge. χ ist dabei die eine Größe, die an allen Stellen den Wechsel misst: Schalter, verlorene Flussdimension, Anstieg der gelesenen Masse und gelesene Dichte χΩ² (außerhalb der Masse) [korr.].
