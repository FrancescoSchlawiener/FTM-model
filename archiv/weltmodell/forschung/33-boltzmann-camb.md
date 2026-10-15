# Boltzmann-Rechnung mit CAMB: die Schließungsregel im CMB

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: 1 − w bei der Rekombination 9·10⁻¹² (ohne Strahlung in H) → mit Strahlung in H(z) 7,0·10⁻¹² bei festem a₀ = 1,2·10⁻¹⁰, 4·10⁻¹² mit a₀(1090) = 0,76·a₀ (Modellwert); heute 3,7·10⁻³ bleibt; Q = 2π − 1 als Hypothese → Zahlennähe ohne Herleitung, Q = 5,364 aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0); „Vorhersage“, Ergebnis 2 und nächster Schritt nur noch als Vergleich bzw. unter der Annahme. Umbenennung (ohne Marke): Δχ² (Chi-Quadrat) → ΔX². [G] ergänzt: Planck-Werte.

> Arbeitsstand, **noch nicht im Gesamtmodell**.
> Rechnung: `rechnungen/camb_q.py` (CAMB 2.0.4)

## Warum CAMB die Schließungsregel direkt rechnet

Im mitatmenden Kosmos ist das Projektionsgewicht praktisch exakt w = 1: 1−w = 4·10⁻¹² bei der Rekombination (mit Strahlung in H(z) und a₀(1090) = 0,76·a₀; bei festem a₀ = 1,2·10⁻¹⁰ 7,0·10⁻¹², 9·10⁻¹² nur ohne Strahlung in H) und 3,7·10⁻³ heute. [korr.] Mit w = 1 sind die Gleichungen des projizierten Flusses linear identisch mit denen eines druckfreien Flusses der Dichte Q·ρ_b. Ein Standard-Boltzmann-Code rechnet die Schließungsregel deshalb ohne Umbau, mit ω_c = Q·ω_b.

**Linear testet der Code den Wert von Q.** Die Abweichungen der Schließungsregel vom Standardmodell entstehen erst nichtlinear beim Kollaps (w → 0, Wellenzahlen oberhalb k₁). Dafür braucht es eine N-Körper-Rechnung.

## Rechnung: Q = 5,364 gegen Q = 2π − 1

[korr.] Der Vergleichswert Q = 2π − 1 ist eine Zahlennähe ohne Herleitung und trägt nicht; Q = 5,364 ist aus der Messung eingesetzt, und im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0). Die Rechnung zeigt nur, wie gut Messungen die beiden Zahlen trennen.

Planck-Parameter [G] (ω_b = 0,02237, n_s = 0,9649, A_s = 2,1·10⁻⁹, τ_reio = 0,0544, Σm_ν = 0,06 eV), mit festgehaltenem Winkel des Schallhorizonts θ* = 0,0104092:

| | ω_c | H₀ | z* | r_drag [Mpc] | z_eq |
|---|---|---|---|---|---|
| Q = 5,364 | 0,11999 | 67,34 | 1089,91 | 147,105 | 3402,7 |
| Q = 2π − 1 = 5,283 | 0,11818 | 67,98 | 1089,75 | 147,585 | 3359,5 |

Relative Änderung der Temperaturleistung:

| ℓ | 2 | 30 | 220 | 540 | 810 | 1130 | 2000 |
|---|---|---|---|---|---|---|---|
| ΔC_ℓ/C_ℓ | +1,59 % | +0,37 % | +0,78 % | +0,71 % | +0,22 % | +0,34 % | −0,06 % |

## Unterscheidbarkeit

| Messung | ΔX², nur ω_c geändert (obere Schranke) |
|---|---|
| Planck-ähnlich (f_sky 0,6, 33 µK', 7,3') | 15,2 (3,9 σ) |
| CMB-S4-ähnlich (f_sky 0,4, 1 µK', 1,4') | 25,6 (5,1 σ) |
| nur kosmische Varianz (f_sky 0,6) | 38,7 (6,2 σ) |

Werden A_s, n_s und H₀ mitangepasst, nimmt ein Teil des Unterschieds andere Parameter auf. Maßgeblich ist dann die marginalisierte Unsicherheit:

| Datensatz | ω_c gemessen | Vergleichswert (2π−1)·ω_b [korr.] | Abstand |
|---|---|---|---|
| Planck 2018 | 0,1200 ± 0,0012 [G] | 0,1182 | 1,5 σ |
| Planck 2018 + BAO | 0,11933 ± 0,00091 [G] | 0,11845 (ω_b = 0,02242) | 1,0 σ |

## Ergebnis

1. Die Schließungsregel reproduziert das gemessene CMB-Spektrum mit Q = 5,364 exakt, weil sie linear mit Staub der Dichte Q·ρ_b übereinstimmt.
2. Die Hypothese „Rückkehrzeit = ein Umlauf" (Q = 2π−1) liegt 1,0–1,5 σ neben den heutigen Messungen und ist damit verträglich. [korr.] Zahlennähe ohne Herleitung, trägt nicht; die Verträglichkeit ist kein Beleg.
   - Sie verschiebt H₀ bei festem θ* um +0,64 km/s/Mpc.
   - Die Leistung bei ℓ ≈ 200–550 steigt um etwa 0,7–0,8 %.
3. Eine Messung mit CMB-S4-Qualität unterscheidet beide Werte ohne Mitanpassung auf etwa 5 σ. Mit Mitanpassung muss σ(ω_c) unter etwa 0,0004 fallen.

## Nächster Schritt

Nur falls Q = 2π − 1 eine Herleitung erhält: eine vollständige Anpassung mit der Planck-Likelihood, z. B. mit Cobaya, in der Q = 2π−1 festgesetzt wird. [korr.] Dazu die Frage, ob der Anstieg von H₀ und die Verschiebung von r_drag zu den BAO-Daten von DESI passen.
