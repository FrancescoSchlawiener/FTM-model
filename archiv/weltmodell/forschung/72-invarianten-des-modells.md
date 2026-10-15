# Invarianten des Modells

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „a₀ = c²/U, ℓ₀ = U (Umfang der atmenden Zelle), r*² = 3Δs·U“ → a₀ = c²/(2Λ), ℓ₀ = 2Λ, Λ innere Länge und heute fest, r*² = 3Δs·ℓ₀; (r*/ℓ_P)² = 4,10·10⁴² / 4,88·10⁹⁹ / 2,93·10¹¹⁰ → 3,56·10⁴² / 4,23·10⁹⁹ / 2,54·10¹¹⁰ (ℓ₀ = c²/a₀, a₀ = 1,2·10⁻¹⁰); M(r* = λ̄) = 5,8·10⁻²⁹ → 6,1·10⁻²⁹ kg; „es fehlt eine Invariante für a₀(z), die bei z ≈ 0,4 umkehrt“ → a₀ ist fest, die Form der Dunklen Energie wird in der Rückstellung gesucht; „Q nur Kac-Lesart 2π − 1“ → Zahlennähe ohne Herleitung, Q aus der Messung; Flächendichte der Scheibe im Quadrat nur außerhalb der Mittelebene. Umbenennungen: U → ℓ₀; Λ (QCD-Skala) → Λ_QCD. Die Zahlen rechneten mit ℓ₀ = 2πc/H₀, also a₀ = cH₀/2π = 1,04·10⁻¹⁰.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Analyse der Struktur; Zahlen mit kurzer Rechnung geprüft (Proton, Sonne, Milchstraße, H₀ = 67,4 km/s/Mpc).

## Fünf Formen, die überall wiederkehren

| Form | Stellen im Modell | Was sie festlegt |
|---|---|---|
| **Produkt fest, Übergang bei Gleichheit** (geometrisches Mittel) | g·R² = 3c²Δs · Δs·λ̄ = ℓ_P²/3 · r*² = 3Δs·ℓ₀ · s² = Γ·τ_D · a_J² = r₊r₋ | Der Übergang liegt dort, wo beide Faktoren gleich sind: Planck-Masse (3Δs = λ̄), r* (a_N = a₀), Verlustzeit = Rückkehrzeit, extremer Kerr-Trichter (r₊ = r₋) |
| **Summe der Rückläufe** (Schließung) | Σqⁿ = 1/(1 − q) · Kac ⟨T⟩ = μ(X)/μ(A), Σ1/⟨T⟩ = 1 · ΣΔs additiv · Rückstellbilanz ΣΔs_innen = ΔL_außen | Bindung und Verstärkung ohne Quelle; das Budget der Wiederkehr ist fest |
| **Dimension als Faktor** | Δs = (GM/c²)·V/(S·R) · d_eff = 3 − χ · s^d·n_B mit Γ(d+1)ζ(d+1) · Drehordnungen 1, 2, 3, 4, 6 | die ganzen Zahlen 3, 6, 9, 12; Fläche gegen Umfang; die Konstanten π²/6, 2ζ(3), π⁴/15 |
| **Quadrat bei Senkrechtem** | κ² + ν_z² = (2 + χ)Ω² · s⁴ = s_i⁴ + s_e⁴ · (Re τ_M)² = 1 + χ · tangential ⊥ transversal | äußeres Feld, Flächendichte der Scheibe (nur außerhalb der Mittelebene, dort ist ihr senkrechtes Feld null) [korr.], gelesene Dichte als offener Pythagoras |
| **Exponent eines Verhältnisses** (Nome) | q = e^{2πiτ} · e^{−s} · α_G als Nome · Λ_QCD/M = e^{−2π/(β₀α_s)} [korr.] | große Unterschiede aus mäßigen Zahlen (10⁻³⁹ aus 1/α_s ≈ 53) |

## Die Übergangsleiter

Nach Kac ist a₀ = c²/(2Λ): Die Länge der Wirkung ist ℓ₀ = c²/a₀ = 2Λ, die doppelte Schließlänge der Rückkehr, eine Länge der inneren Richtung, keine Raumlänge; Λ ist heute fest. [korr.] Damit haben beide Enden des Modells dieselbe Form:

```
klein:  ℓ_P² = 3Δs·λ̄        Faltung gegen Quantenlänge (λ̄ = ħ/Mc, lokal, ∝ 1/M)
groß:   r*²  = 3Δs·ℓ₀       Faltung gegen Länge der Wirkung (ℓ₀ = 2Λ, global)
        (r*/ℓ_P)² = ℓ₀/λ̄    für jede Masse
```

Geprüft mit ℓ₀ = c²/a₀ (a₀ = 1,2·10⁻¹⁰): (r*/ℓ_P)² = ℓ₀/λ̄ = 3,56·10⁴² (Proton), 4,23·10⁹⁹ (Sonne), 2,54·10¹¹⁰ (Milchstraße). [korr.] Die Faltung Δs steht zwischen zwei Gegenlängen, einer lokalen und einer globalen. Das ist dieselbe Paarung wie bei der Engstelle, lokal r₀ = 6Δs, global r₀,H = 6Δs_H = c/H; heute ist ℓ₀ ≈ 2π·r₀,H. [korr.]

Beide Übergänge fallen zusammen, wo der Schwellenradius die Quantenlänge erreicht: r* = λ̄ bei M = (ħ²a₀/(Gc²))^{1/3} = 6,1·10⁻²⁹ kg [korr.]. Darunter liegt r* innerhalb der eigenen Quantenlänge; die galaktische Regel hat dort keinen Bahnsinn mehr.

## Wo das Modell noch keine Invariante hat

- **Form der Dunklen Energie:** Alle gerechneten Größen (Horizonte, 16 weitere, Kac, gebundene Faltung) wachsen oder fallen einseitig. a₀ ist fest (seit z ≈ 10), Λ heute fest; der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie. Die Form der übrigen zwei Drittel (Maximum bei z ≈ 0,4) wird in der Rückstellung gesucht, als Bindungsrate dF/d ln a. [korr.]
- **Q = 5,36:** keine Gleichung; Q = 5,364 ist aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0). Die Kac-Lesart 2π − 1 ist Zahlennähe ohne Herleitung. [korr.]
- **m_p/m_e, Im τ_X:** In allen Produkt-Invarianten kürzt sich die Masse; eine Gleichung, die sie festlegt, fehlt.
- **DF2/DF4 gegen Crater II:** Es fehlt eine Größe, die im starken Gezeitenfeld die Verstärkung senkt, im schwachen nicht.

## Beobachterfehler

Alle fünf Formen sind Beziehungen zwischen Längen, Zeiten und Regeln. Keine führt einen Stoff ein: Die Bose-Form zählt Rückläufe, die Dichte χΩ²/4πG ist eine gelesene, der Pythagoras der Epizyklen öffnet sich durch die Regel, nicht durch Masse.
