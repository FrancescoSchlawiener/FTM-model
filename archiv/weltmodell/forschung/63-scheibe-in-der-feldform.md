# Die Scheibe in der Feldform

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Die Scheibe wirkt auf die eigene Radialbewegung wie ein senkrechtes Außenfeld“ → gilt nur in der Feldform; nach der örtlichen Regel ist g_z in der Mittelebene null, die eigene Flächendichte senkt ν dort nicht (SPARC: f = +0,03 bis +0,13, Feldform außerhalb 95 %). Folgen: neben Kern- und Feldform die örtliche Regel a = ν(|g_N|)·g_N ergänzt; sie gilt.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/scheibe_feldform.py` (Ausgabe `scheibe_feldform.txt`)
> Feldform ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N], achsensymmetrisch in (R, z), Finite Volumen, 300 × 300 Zellen bis 80 r*. Einheiten G = M = a₀ = 1.

## Prüfung des Lösers

- **Plummer-Kugel** (Feldform = algebraische Regel exakt): Abweichung höchstens 10⁻⁴ von 0,1 bis 10 r*.
- **Kuzmin-Scheibe:** Oberhalb der Ebene ist ihr Newton-Feld das einer Punktmasse unter der Ebene, ν·g_N ist dort wirbelfrei, also exakt g = ν(|g_N|)·g_N. Der Löser trifft das auf 0,1–1,8 %.

## Das Ergebnis

In der Ebene einer dünnen Scheibe setzt sich der Betrag des Newton-Felds senkrecht zusammen, aus dem radialen Anteil g_R und dem vertikalen 2πGΣ:

```
Feldform:      a_R = ν(|g_N|)·g_R,   |g_N|² = g_R² + (2πGΣ)²     (s⁴ = s_R⁴ + s_z⁴)
algebraisch:   a_R = ν(g_R)·g_R
Kuzmin, tief:  a_Feld/a_alg = √(g_R/|g_N|) = √(cos θ),   θ Neigung des Newton-Felds gegen die Ebene
```

In der Feldform wirkt die Scheibe auf die eigene Radialbewegung wie ein senkrechtes Außenfeld: Ihre Flächendichte senkt ν. Nach der örtlichen Regel a = ν(|g_N|)·g_N gilt das nicht: In der Mittelebene ist g_z aus der Spiegelung null, die eigene Flächendichte senkt ν dort nicht. Die Bahnen folgen der örtlichen Regel, die Metrik (Feldform) trägt nur deren wirbelfreien Teil; an SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03, die Feldform (f = 1) liegt außerhalb 95 %. [korr.] Exponentialscheibe, Verhältnis Feldform zu algebraischer Regel in der Ebene:

| R/R_d | 0,25 | 0,5 | 1 | 2 | 3 | 5 | 8 |
|---|---|---|---|---|---|---|---|
| R_d = 0,3 r* (hohe Flächendichte) | 0,820 | 0,878 | 0,922 | 0,966 | 0,992 | 1,016 | 1,013 |
| R_d = 1 r* | 0,605 | 0,715 | 0,829 | 0,936 | 0,988 | 1,021 | 1,016 |
| R_d = 3 r* (niedrige Flächendichte) | 0,528 | 0,658 | 0,793 | 0,924 | 0,986 | 1,023 | 1,017 |

Innen liegt die Feldform unter der algebraischen Regel, umso stärker, je niedriger die Flächendichte: bei R_d/4 0,82 (R_d = 0,3 r*) bis 0,53 (R_d = 3 r*), bei R_d 0,79–0,92. Bei etwa 3,5 R_d kreuzen beide, außen (5–8 R_d) liegt die Feldform 1,3–2,3 % darüber. Weit draußen, wo das Feld kugelförmig wird, fallen beide zusammen.

## Folgen

- Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N (Projektor n⊗n, am Ort der Bahn gezählt); sie gilt. In der Mittelebene ist g_z = 0, sie gibt dort genau die algebraische Form ν(g_R)·g_R; für Kugeln ist sie gleich Kern- und Feldform. Die Feldform ist ihr wirbelfreier Teil, den die Metrik trägt. [korr.]
- Kernform (radialer Kern mit dem Radialfeld) und Feldform geben für Kugeln dieselbe Rotationskurve, für Scheiben nicht. Sie sind damit nicht nur an der Vertikaldynamik trennbar, sondern an den inneren Rotationskurven von Scheiben niedriger Flächendichte.
- In der Feldform ist die Beschleunigungsrelation für Scheiben keine einzige Kurve in g_R allein; sie hängt zusätzlich vom Verhältnis 2πGΣ/g_R ab. Prüfbar: innere Punkte von Galaxien niedriger Flächendichte müssten systematisch unter der algebraischen Relation liegen, um so mehr, je größer 2πGΣ/g_R.
- Das Zusammensetzen im Quadrat ist dasselbe wie beim äußeren Feld (s⁴ = s_i⁴ + s_e⁴) und folgt aus dem Argument der Wirkung, s⁴ = ℓ₀²·|∇σ|².
