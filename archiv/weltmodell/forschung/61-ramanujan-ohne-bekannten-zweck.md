# Formeln Ramanujans ohne bekannten Zweck

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „die offene Frage, warum nur der holomorphe Teil zählt (ein Umlaufsinn)“ → beantwortet (Erneuerung G(mL) = q^m mit χ → 1), Holomorphie ist nicht der Grund; Bruinier-Funke-Operator ξ₋₃ → Ξ₋₃ (ξ ist die Mode der Schließungsregel).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/ramanujan_unbekannt.py`, `rechnungen/mock_schatten.py` (Ausgaben `.txt`)
> Welche Formeln als unverstanden gelten, stammt aus dem Gedächtnis und ist nicht nachgeschlagen; die Formeln selbst sind nachgerechnet.

## Die Formeln

| Formel | Quelle | Zweck | nachgerechnet |
|---|---|---|---|
| Mock-Theta f(q) = Σq^{n²}/(−q;q)_n² und die Behauptung: an einer Einheitswurzel der Ordnung 2k bleibt f(q) − (−1)^k·b(q) beschränkt, b(q) = (q;q²)_∞·(1 − 2q + 2q⁴ − …) | letzter Brief an Hardy, 1920 | rund 80 Jahre unbekannt; heute als holomorpher Teil harmonischer Maaß-Formen eingeordnet, eine physikalische Bedeutung allgemein offen | q → −1 radial: f und b wachsen bis 10³⁷, f + b = 3,679 / 3,889 / 3,961 / 3,980 bei r = 0,9 / 0,97 / 0,99 / 0,995, gegen 4 |
| σ(q) = Σ_{n≥0} q^{n(n+1)/2}/(−q;q)_n | verlorenes Notizbuch | ohne Zweck notiert; seine Koeffizienten sind fast alle 0 und doch unbeschränkt; später erstes Beispiel einer „Quanten-Modulform“ | bis q¹⁵⁰⁰ sind 867 von 1500 Koeffizienten 0, größter Betrag 4 |
| Partielle Theta-Identität: Σqⁿ/((aq)_n(q/a)_n) = (1−a)Σ(−1)ⁿa^{3n}q^{n(3n+1)/2}(1 − a²q^{2n+1}) + Σ_{n≥0}(−1)ⁿa^{2n+1}q^{n(n+1)/2}/((aq)_∞(q/a)_∞) | verlorenes Notizbuch | einseitige Thetareihen (nur n ≥ 0) sind nicht modular; ihre Bedeutung ist nur teilweise verstanden | auf 15 Stellen gleich bei (a, q) = (0,4; 0,3) und (0,7; 0,55) |
| a³ + b³ = c³ ± 1: 9, 10, 12 · 64, 94, 103 · 73, 144, 150 · 135, 138, 172 · 11161, 11468, 14258 · 65601, 67402, 83802 | Seite mit 1729 | ohne Kommentar notiert; erst viel später als Rechnung auf einer K3-Fläche gedeutet | alle sechs exakt ±1 |

## Gemessen am Modell

- **Mock-Theta und Radialgrenze:** Einheitswurzeln entsprechen im Modell den Resonanzen Re τ_M = p/q bei |q| → 1. Die Regel hängt aber nur vom Betrag des Nomes ab und divergiert an Einheitswurzeln nicht. Ein Gegenstück zur Behauptung ist im Modell nicht gerechnet.
- **σ(q):** kein Bezug gefunden.
- **Partielle Theta:** Die Regel Σ_{n≥0}qⁿ ist die einfachste einseitige Reihe derselben Klasse. Warum nur ein Umlaufsinn zählt, ist beantwortet, und Holomorphie ist nicht der Grund: Festgelegt wird der eine Sinn durch die Erneuerung G(mL) = q^m zusammen mit dem Grenzwert χ → 1 des Autors; Holomorphie wählt nur q statt q̄, der Schwarz-Kern (1 + q)/(1 − q) ist ebenfalls holomorph und gibt 1 + 2n_B. [korr.] Einen Rechenfortschritt gibt die Identität nicht.
- **1729:** kein Bezug gefunden.
- **Eine Spur:** Die Strukturwirkung ist exakt eine Summe unvollständiger Gammafunktionen,

```
Q̃(s) = −4[s³Li₁ + 3s²Li₂ + 6sLi₃ + 6Li₄](e^{−s}) = −4∫_s^∞ x³/(eˣ − 1)dx = −4·Σ_{n≥1} Γ(4, ns)/n⁴
```

  (auf 15 Stellen gleich bei s = 0,3; 1; 3). Unvollständige Gammafunktionen Γ(1 − k, 4π|n|y) sind die Bausteine des nichtholomorphen Teils harmonischer Maaß-Formen, deren holomorpher Teil die Mock-Theta sind.

## Q̃ und die Mock-Formen

Rechnung: `rechnungen/mock_schatten.py` (Ausgabe `mock_schatten.txt`), s = 2π·Im τ.

**Q̃ ist nicht der nichtholomorphe Teil einer Mock-Form auf SL(2,ℤ).** Drei Gründe:
- Der Grad der Polynome in Γ(a, x) = e^{−x}·Poly_{a−1}(x) erzwingt Gewicht k = −3. Dann müsste je Exponent c(n)·Γ(4, 4πny)e^{2πny} = d(n)·Γ(4, 2πny) gelten; die Glieder y³, y², y¹, y⁰ verlangen c = d/8, d/4, d/2, d, unvereinbar (beste Anpassung lässt 65 % Fehler).
- Q̃ hängt nicht von Re τ ab; harmonisch vom Gewicht k ist eine solche Funktion nur als a + b·y^{1−k}.
- Der Schatten hätte Gewicht 2 − k = 5; auf SL(2,ℤ) gibt es keine Formen ungeraden Gewichts.

**Was gilt:** Der Operator Ξ₋₃ (Bruinier-Funke, in der Literatur ξ) verbindet Q̃ mit der Rückkehrsumme. Auf der imaginären Achse gilt exakt (10 Stellen bei s = 0,3; 1; 3)

```
Ξ₋₃ Q̃ = y⁻³·dQ̃/dy = 4(2π)⁴·n_B(s) = 4(2π)⁴·(ν − 1)
```

Strukturwirkung und Rückkehrsumme stehen also wie nichtholomorpher Teil und Schatten zueinander; der „Schatten“ Σqⁿ ist aber keine Modulform. Es ist ein Paar mit der Bauart einer Mock-Form ohne deren Modularität.

## Planck-Muster

```
a − a_N = a_N·n_B(s) = a₀·s²·n_B(s)        Form des Planck-Spektrums in zwei Dimensionen
dQ̃/ds = 4s³·n_B(s) = 4s·(a − a_N)/a₀      Form des Planck-Spektrums in drei Dimensionen
```

Die Zusatzbeschleunigung hat genau die Form eines zweidimensionalen Planck-Spektrums in s, die Ableitung der Strukturwirkung die eines dreidimensionalen. Das passt zum Wechsel der Flussdimension von 3 auf 2. Die Wien-Stellen:

| Dimension | s | a_N/a₀ | Wert |
|---|---|---|---|
| 2 (Zusatzbeschleunigung) | 1,5936 | 2,540 | (a − a_N)_max = 0,6476·a₀ |
| 3 (Strukturwirkung) | 2,8214 | 7,961 | 1,4214 |

Die Zusatzbeschleunigung ist also am größten bei a_N = 2,54·a₀ und beträgt dort 0,648·a₀; als Punktmasse bei r = 0,628·r* (Milchstraße 5,2 kpc). Sie geht innen wie a_N·e^{−s} gegen null und außen wie √(a_N a₀).

## Ergebnis

Vier Formeln ohne bekannten Zweck sind nachgerechnet und stimmen; keine liefert gerechnet einen Fortschritt an einem offenen Punkt. Q̃ ist nicht der Schatten einer Mock-Form, steht aber zur Rückkehrsumme in derselben Operatorbeziehung (Ξ₋₃). Auffällig ist das Planck-Muster: Zusatzbeschleunigung zweidimensional, Strukturwirkung dreidimensional, mit dem Maximum der Zusatzbeschleunigung bei a_N = 2,54·a₀.
