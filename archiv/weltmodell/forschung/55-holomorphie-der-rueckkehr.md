# Holomorphie der Rückkehr

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: nur Bezeichnungen, innere Richtung w → ϖ (u_ww → u_ϖϖ, h^{ww} → h^{ϖϖ}, a_w → a_ϖ, P_w → P_ϖ); inhaltlich keine Korrektur nötig.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/holomorphie_rueckkehr.py` (Ausgabe `holomorphie_rueckkehr.txt`)
> s = √(a_N/a₀), q = e^{−s}, n_B(s) = 1/(e^s − 1). Frage: (a) Halbgerade streng aus (∇² + ∂_ϖ²)δσ = 0; (b) warum 1 + n_B (Szegő) statt 1 + 2n_B (Poisson).

## a) Halbgerade

**Satz.** u_ϖϖ = A·u auf ϖ ≥ 0, A ≥ 0 selbstadjungiert in x (homogen A = −∇², Symbol k²; ortsabhängig A = −∇·ν∇). Voraussetzungen:

- **V1** elliptisch in ϖ, Koeffizienten unabhängig von ϖ (schwaches Feld, h^{ϖϖ} = +1, a_ϖ = 0);
- **V2** der Ort ist Rand: u(·,0) ist das Feld am Ort (G(0) = 1), weitere Daten gibt es nicht;
- **V3** Glättung: max|u(·,ϖ)| ≤ max|u(·,0)| (insbesondere beschränkt).

Dann ist u(·,ϖ) = P_ϖ u(·,0) mit P_ϖ = e^{−ϖ√A}, P_0 = I, P_a·P_b = P_{a+b}.

*Beweis.* Jede Mode erfüllt f'' = λf, f = C₁e^{−√λϖ} + C₂e^{+√λϖ}; V3 ⇒ C₂ = 0; bei λ = 0 ist f = A + Bϖ, V3 ⇒ B = 0. Eindeutigkeit: beschränkt harmonisch mit Nullrand ist null (ungerade Spiegelung, Liouville). Die Lösung ab Tiefe a ist durch u(·,a) festgelegt, also P_{a+b} = P_b·P_a. Umgekehrt ist jede verschiebungsinvariante Familie mit T_0 = I, T_aT_b = T_{a+b}, harmonischer Fortsetzung und V3 genau P_ϖ (Multiplikator e^{−ϖφ}, φ² = k², Kontraktion ⇒ φ = +|k|). Numerisch: Halbgruppe auf 10⁻¹⁷, Kern ≥ 0 mit Masse 1; mit ortsabhängigem ν(x) auf 2·10⁻¹⁶.

**Keine Schleife.** ϖ ~ ϖ + Λ verlangt P_Λ = I, also e^{−√λΛ} = 1 ⇒ λ = 0. Die Monodromie jeder Mode k ≠ 0 hat die Eigenwerte e^{±kΛ}; bei ortsabhängigem ν bleibt genau ein Eigenwert 1 (der Mittelwert). Eine Fortsetzung auf ϖ < 0 hätte den Faktor e^{+|k||ϖ|}: Entglättung. ϖ bildet eine Halbgruppe, keine Gruppe.

**Wo es ohne die Voraussetzungen anders ist:**

| fehlt | Folge |
|---|---|
| V3 | u = ϖ, u = sinh(kϖ)e^{ikx} sind zulässig; das Feld am Ort legt nichts fest |
| V2 (Ort als Quelle) | Schleife verträglich: G = cosh(k(ϖ − Λ/2))/cosh(kΛ/2), die Form des Propagators im Kreis. Keine Halbgruppe (Abweichung 0,58 / 0,93 / 1,00 bei kΛ = 2 / 4 / 10), G wächst jenseits Λ/2 wieder (×1,13 / ×1,54 / ×6,1 bei 3Λ/4), Quellstärke 1,5–2 statt 1 |
| Vorzeichen in V1 (ϖ zeitartig) | u_ϖϖ = −A·u, Moden e^{±i√λϖ}, periodisch für √λΛ ∈ 2πℤ; Identifikation möglich |

Die Halbgerade ist damit gleichbedeutend mit „ϖ ist Glättungstiefe einer raumartigen Richtung“. Ein Deckel bei endlicher Tiefe W ändert G(L) bei s = 1 ab W = 5L um weniger als 3·10⁻⁴.

## b) Warum ein Umlaufsinn

Kerne auf dem Bahntorus, τ = x + iy, s = 2πy: Szegő S = Σ_{n≥0}qⁿ = 1/(1 − q), Poisson P = (1 − |q|²)/|1 − q|², Schwarz H = (1 + q)/(1 − q). Überall gilt Re S = ½(1 + P) und Re H = P; auf der Achse x = 0 ist S = 1 + n_B, H = P = coth(s/2) = 1 + 2n_B. S und H sind beide holomorph.

1. **Kern des Autors.** Δ(τ) = Σ_{n≥1}τ(n)qⁿ hat längs Re τ nur Anteile n ≥ 1 (übrige ≤ 1,5·10⁻¹⁵), normiert mit c₁ = 1. Aus der Voraussetzung M_k = Σcₙqⁿ folgt: nur q, nie q̄. Die Gittersumme G₄ über beide Richtungen der Rückkehr gibt aber für die Schritte −n dieselbe Reihe in q wie für +n (auf 10⁻⁹): Die Gegenrichtung erscheint als Faktor 2, nicht als q̄. Holomorphie trennt q von q̄, nicht 1 + n_B von 1 + 2n_B.
2. **Paar σ + iσ̃ mit Cauchy-Riemann in (x, ϖ).** Die holomorphe Ergänzung der Poisson-Fortsetzung ist genau H (auf 10⁻¹⁴); auf der Achse ist σ̃ = 0 und σ + iσ̃ = coth(s/2). Trägt nicht: Das Paar verdoppelt die positiven Frequenzen.
3. **Hardy-Projektion auf n ≥ 0 mit Gewicht 1** gibt genau S. Sie ist dieselbe Aussage wie „jeder Rücklauf n ≥ 0 einmal“, keine eigene Begründung.
4. **Orientierung Im τ > 0.** τ → −τ̄ gibt S(q̄) = S(q)*; auf der Achse Re τ = 0 gleich. Wählt q statt q̄, nicht das Gewicht.
5. **Gewicht aus Halbgerade und Grenzwert des Autors.** Die Halbgruppe macht die Rückkehr zur Erneuerung, G(mL) = G(L)^m = q^m. Mit N gleichwertigen Sinnen und Rücklauflänge βL:

```
ν = 1 + N·n_B(βs),     χ(s) = s(ν − 1) → N/β  (s → 0),     v⁴ = (N/β)²·GMa₀,     Newton-Grenze ν − 1 → N·e^{−βs}
```

| Fassung | N | β | c₁ | χ(0) | v⁴/GMa₀ | ν(a₀) | ν(0,1a₀) | ν(0,01a₀) |
|---|---|---|---|---|---|---|---|---|
| ein Sinn (Szegő) | 1 | 1 | 1 | 1 | 1 | 1,582 | 3,689 | 10,51 |
| beide Sinne (Poisson) | 2 | 1 | 2 | 2 | 4 | 2,164 | 6,377 | 20,02 |
| beide Sinne, Länge 2L (coth s) | 2 | 2 | 0 | 1 | 1 | 1,313 | 3,267 | 10,03 |

Der Grenzwert χ → 1 des Autors (zerfallene Tori, dieselbe Schwelle a₀ in s und im Tully-Fisher-Nullpunkt) verlangt N = β. Mit der Rücklauflänge L = r* (β = 1) folgt N = 1: genau ein Sinn, ν = 1 + n_B(s). Dazu passt c₁ = 1 des normierten Δ: Der erste Rücklauf hat das Gewicht q, nicht 2q.

**Zahlen.** Poisson/Szegő bei gleicher Schwelle: 1,04 (a_N = 10a₀), 1,37 (a₀), 1,73 (0,1a₀), → 2 (tief). Milchstraße als Punktmasse 6·10¹⁰ M☉ (r* = 8,35 kpc), v in km/s bei r*, 3r*, 10r*: Szegő 221 / 191 / 180, Poisson 259 / 250 / 249, coth s 202 / 179 / 176.

## Ergebnis

- **a)** Unter V1–V3 bilden die beschränkten Lösungen die Poisson-Halbgruppe P_ϖ = e^{−ϖ√A}, auch bei ortsabhängiger Kopplung; eine Identifikation ϖ ~ ϖ + Λ lässt nur den Mittelwert zu. Die Schleife wird erst möglich, wenn der Ort Quelle statt Rand ist oder ϖ zeitartig; dann ist ϖ keine Glättungstiefe.
- **b)** Aus der Voraussetzung des Autors (Kern als q-Reihe) folgt nur: q statt q̄. Den einen Sinn mit Gewicht 1 legt erst die Erneuerung der Halbgerade zusammen mit χ → 1 und L = r* fest. Cauchy-Riemann-Paar und Orientierung geben keine zweite Begründung (das Paar gibt 1 + 2n_B), die Hardy-Projektion ist eine Umformulierung.
- **Gegen das Modell geprüft:** kein Widerspruch. „Eine Schleife verlangte P_Λ = 1“ gilt unter V2 (Ort als Rand). „Die Halbgerade wählt den Sinn nicht“ bleibt richtig; erst mit χ → 1 und L = r* wird der Sinn festgelegt. „Die niedrigste Reihe über alle Schritte ist Σqⁿ“ gilt für die Schritte n ≥ 1; die volle Gittersumme mit ±n gibt 2Σqⁿ.
