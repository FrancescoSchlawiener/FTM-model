# Welche Identifikation trägt 𝒫 in einer Galaxie?

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Gewicht |q| = e^{−2π Im τ_M} als Tatsache → als Punkt 1 des Autors gekennzeichnet (Modell: Punkt 3, |q_Δs| = Δs/r; Wahl offen, betrifft nur die Schreibweise); der eine Sinn „aus Kausalität, kausale Halbgruppe, aus der Vergangenheit“ → festgelegt durch die Erneuerung G(mL) = q^m zusammen mit χ → 1, Kausalität erzwingt ihn nicht; „offen: warum das Gewicht genau e^{−s} ist“ → im Modell hergeleitet (D·G'' = Γ·G, Γ = a_N/c, τ_D = 2Λ/c).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/identifikation_galaxie.py` (Ausgabe `identifikation_galaxie.txt`)
> Beispiel: baryonische Milchstraße 6·10¹⁰ M☉, r* = 8,35 kpc.

**Frage.** Schließung bindet, wenn eine Translationsgruppe den Raum identifiziert: 𝒫 = (1/M)Σ T_{mL} vertauscht mit dem Generator, lässt nur die invarianten Moden durch und öffnet eine Spektrallücke. Welche Identifikation leistet das im Zustandsraum einer Galaxie? Der Autor bildet den Phasenraum auf einen parametrischen Torus mit zwei Zyklen ab, Umlauf und Rückkehr, τ_M = ω_rad/ω_azim + i·Im τ_M.

## 1  Die realen Richtungen binden nicht

Die Winkel des Bahntorus sind auf jeder gebundenen Bahn schon identifiziert (θ ~ θ + 2π); das ist Newton. Auf einer offenen Bahn (E ≥ 0) gibt es keinen Torus. Die Winkelschließung folgt also aus der Bindung, sie erzeugt keine.

Drehbilder einer Balkenverteilung (200 Massen, Feldpunkt 20 kpc): n-zählige Drehidentifikation löscht ab n = 3 den nichtaxialen Anteil (m = 2 von 1,4 % auf 10⁻¹⁷), der Monopol ⟨a_r⟩ bleibt exakt gleich. Drehungen und Umlaufphase (Re τ_M) entfernen Fluchtmoden des Winkels, nicht die radiale Flucht.

## 2  Die Rückkehr bindet

Bindend ist allein die zweite Richtung: die Rückkehr mit Gewicht |q| = e^{−2π·Im τ_M}. [korr.] Das ist Punkt 1 des Autors (Im τ_M = r*/(2πr), |q| = e^{−r*/r}, für die Punktmasse der Rückkehranteil); das Modell schreibt den Kern-Nome |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise. Als Operator ist das eine Translation T (Verschiebung der Geschichte um eine Rückkehr), gewichtet und kausal summiert:

```
𝒫_q = (1 − q)·Σ_{m≥0} qᵐTᵐ = (1 − q)(I − qT)⁻¹          Abel-Mittel der Rückkehr
q → 0:  𝒫_q = I   (offen, Newton, Spitze)                q → 1:  𝒫_q → 𝒫   (voll geschlossen)
(I − qT)⁻¹ = Rückkopplungssumme:  invariante Mode × 1/(1 − q),  Fluchtmoden beschränkt
```

Gerechnet für quasiperiodische Rückkehr (Drehung um den goldenen Winkel, 4000 Moden): Die invariante Mode wird um 10, 100, 1000, 10⁴ verstärkt (q = 0,9 … 0,9999), die Fluchtmoden höchstens um 568 = 1/min|1 − λ|. ‖𝒫_q − 𝒫‖ fällt erst nahe q = 1 (0,49 bei 0,999; 0,057 bei 0,9999); wie schnell, bestimmen die kleinsten Abstände der Rückkehrfrequenzen von der Resonanz.

Damit ist die Galaxie ein Fall von Satz 1: Die Identifikation ist eine Translation, nämlich die Rückkehr, und die Beschleunigungsrelation a = a_N/(1 − q) ist ihre Resolvente. Nahfeld heißt offen, Fernfeld heißt geschlossen. Die Verstärkung ν = 1/(1 − q) ist die mittlere Zahl der Durchläufe; im Fernfeld gilt ν = r/r* + ½:

| r/r* | s | \|q\| | Im τ_M (Punkt 1) [korr.] | ν | r/r* + ½ |
|---|---|---|---|---|---|
| 0,3 | 3,33 | 0,036 | 0,53 | 1,037 | 0,8 |
| 1 | 1 | 0,368 | 0,16 | 1,582 | 1,5 |
| 3 | 0,33 | 0,717 | 0,053 | 3,528 | 3,5 |
| 10 | 0,1 | 0,905 | 0,016 | 10,508 | 10,5 |

## 3  Nur die Halbgruppe [korr.]

Summiert man über die volle Gruppe (m ∈ ℤ, Rückkehr in beiden Richtungen) statt über die Halbgruppe m ≥ 0, entsteht bei gleicher Tully-Fisher-Grenze coth s statt 1/(1 − e^{−s}):

| a_N/a₀ | 10 | 3 | 1 | 0,3 | 0,1 | 0,01 |
|---|---|---|---|---|---|---|
| 1/(1 − e^{−s}) | 1,044 | 1,215 | 1,582 | 2,371 | 3,689 | 10,51 |
| coth s | 1,004 | 1,065 | 1,313 | 2,005 | 3,267 | 10,03 |

Die volle Gruppe liegt bei a_N = a₀ 17 % tiefer. Die Galaxie trägt also die Halbgruppe der Rückkehr (m ≥ 0). Festgelegt wird dieser eine Sinn durch die Erneuerung G(mL) = q^m zusammen mit dem Grenzwert χ → 1 des Autors; Kausalität erzwingt ihn nicht, denn beide Windungssinne sind retardiert, und Holomorphie auch nicht (der Schwarz-Kern (1 + q)/(1 − q) ist ebenfalls holomorph und gibt 1 + 2n_B). [korr.]

## 4  Modulformen als Bildsummen

Summiert man einen Potenzkern z^{−k} über das Gitter ℤ + τℤ (das ist 𝒫 auf dem Torus), entsteht die Eisenstein-Reihe G_k(τ). Gerechnet bei τ = 0,3 + 1,1i und τ = i/2π (s = 1): Gitter und q-Reihe stimmen auf 10⁻¹⁰ bis 10⁻¹⁵. Die Summe längs des Umlaufs macht aus jedem Rückkehrschritt n eine Reihe in qⁿ:

```
Σ_m (m + nτ)^{−k} = (−2πi)^k/(k−1)! · Σ_{d≥1} d^{k−1} q^{nd}        (n = 1, 2, 3 auf 10 Stellen)
```

Die Modulformen des Autors sind damit die Bildsummen der Identifikation, und die Modulgruppe ist die freie Wahl, welcher Gittervektor Umlauf und welcher Rückkehr heißt. Die niedrigste Reihe (d = 1) über alle Rückkehrschritte ist Σ qⁿ = q/(1 − q), der Zusatz der Kreisbahn.

## Ergebnis

- 𝒫 trägt in der Galaxie die **Rückkehr als Translation**, kausal summiert mit Gewicht |q| = e^{−2π Im τ_M} (Punkt 1 des Autors; Wahl des Imaginärteils offen, nur Schreibweise [korr.]). Die realen Richtungen (Winkel, Drehungen, Umlaufphase) sind schon Newton und binden nicht zusätzlich.
- Die Beschleunigungsrelation ist die Resolvente dieser Translation; Nahfeld offen (𝒫_q = I), Fernfeld geschlossen (𝒫_q → 𝒫).
- Nur die Halbgruppe (ein Sinn) gibt die gemessene Form; die volle Gruppe (coth s) liegt 17 % tiefer. Festgelegt wird der Sinn durch die Erneuerung zusammen mit χ → 1, nicht durch Kausalität. [korr.]
- Die Modulformen folgen als Gittersummen aus der Identifikation.

Das Gewicht je Rückkehr e^{−s} ist im Modell hergeleitet: Der Rückweg ist eine abgeschirmte Ausbreitung D·G'' = Γ·G mit Γ = a_N/c und τ_D = L²/D = 2Λ/c (Kac), also q = G(L) = e^{−s}, Punktmasse L = r*. [korr.] **Offen:** wie die Rückkehr auf exzentrische Bahnen wirkt. [korr.] Beantwortet im Modell: Die Rückkehr ist eine stationäre Schleife am festen Ort; exzentrische Bahnen werden nicht gerundet, es gibt keinen Nachlauf längs der Bahn.
