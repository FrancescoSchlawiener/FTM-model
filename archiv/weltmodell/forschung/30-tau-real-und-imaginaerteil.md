# τ_M: Realteil und Imaginärteil haben verschiedene Rollen

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Im τ_M = r*/(2πr), damit ist die Frage nach der Definition von τ_M beantwortet“ → das ist Punkt 1 des Autors, das Modell schreibt |q_Δs| = Δs/r (Punkt 3); die Wahl ist offen und betrifft in der stationären Regel nur die Schreibweise. Bezeichnung getrennt: q ist hier der Nome in Punkt 1 (für die Punktmasse der Rückkehranteil), der Kern-Nome heißt q_Δs. Tabelle (κ/ω_orb, ν) bleibt.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/tau_realteil.py`

## Ausgangspunkt

Der Autor definiert den Torus-Parameter als

```
τ_M = ω_rad/ω_azim + i·Im τ_M,        q = e^{2πiτ_M} = e^{−2π·Im τ_M} · e^{2πi·ω_rad/ω_azim}
```

Der Nome zerfällt in einen **Betrag** |q| = e^{−2π·Im τ_M} und eine **Phase** e^{2πi·Re τ_M}.

In dieser Notiz ist Im τ_M = r*/(2πr) gesetzt (Punkt 1 des Autors); dann ist |q| = e^{−r*/r} für die Punktmasse der Rückkehranteil. Der Kern-Nome des Modells heißt q_Δs, mit |q_Δs| = Δs/r (Punkt 3). [korr.]

## Modul-Symmetrie T

q ist invariant unter τ_M → τ_M + 1. Es zählt also nur Re τ_M mod 1.

- **Kepler-Bereich (Newton, g ≫ a₀):** ω_rad/ω_azim = 1 exakt (geschlossene Ellipsen). Wegen τ_M → τ_M+1 ist das gleichwertig mit Re τ_M = 0. **Die Definition des Autors und die rein imaginäre Form τ_M = i·Im τ_M stimmen im Newton-Bereich exakt überein.**
- **Flacher Bereich (g ≪ a₀):** κ/ω_orb → √2, also Re τ_M mod 1 → 0,414. Die Phase dreht q aus der reellen Achse.

## Rechnung (Milchstraße, 6·10¹⁰ M☉, r* = 8,35 kpc)

Für die Rotationskurve ist ω_rad/ω_azim = κ/ω_orb = √(2(1 + d ln v/d ln r)), berechnet aus der Kurve selbst.

| r [kpc] | g_N/a₀ | κ/ω_orb | Re τ_M mod 1 | ν aus \|q\| | Re ν mit Phase |
|---|---|---|---|---|---|
| 0,5 | 279 | 1,0000 | 0,000 | 1,000 | 1,000 |
| 2 | 17,4 | 1,032 | 0,032 | 1,016 | 1,015 |
| 5 | 2,79 | 1,178 | 0,178 | 1,232 | 1,054 |
| 8,35 | 1,00 | 1,258 | 0,258 | 1,582 | 0,869 |
| 30 | 0,077 | 1,367 | 0,367 | 4,116 | 0,583 |
| 120 | 0,0048 | 1,402 | 0,402 | 14,88 | 0,519 |

Die gemessene Radial-Acceleration-Relation verlangt die Spalte „ν aus |q|". Mit der Phase fiele ν im Außenbereich unter 1, und die Bindung würde schwächer statt stärker.

## Ergebnis: die Rollen der beiden Teile

| Teil von τ_M | Größe | Rolle im Zustandsraum | Wirkung |
|---|---|---|---|
| **Im τ_M** | Betrag \|q\| = e^{−2π·Im τ_M} = e^{−r*/r} (Punkt 1; mit Punkt 3 \|q_Δs\| = Δs/r) [korr.] | parametrisch/transversal | **Bindung**: Rückkopplungssumme ν = 1/(1−\|q\|) |
| **Re τ_M = ω_rad/ω_azim** | Phase e^{2πi·Re τ_M} | tangential | **Umlaufphase**: Rosettenform und Präzession der Bahn, keine Bindung |

```
Bindung:     ν = 1/(1 − |q|),     |q| = e^{−2π Im τ_M} = e^{−r*/r},     Im τ_M = r*/(2πr)
Umlauf:      Phase φ = 2π·(ω_rad/ω_azim mod 1)   → Periheldrehung pro Umlauf
```

Die Festlegung Im τ_M = r*/(2πr) ist Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise. [korr.] Beide Teile der Definition des Autors bleiben und haben getrennte Aufgaben:
- Der Imaginärteil trägt die Bindung.
- Der Realteil trägt die Phase der Umläufe.

Im Newton-Bereich ist die Phase wegen der T-Symmetrie exakt null. Das ist der geschlossene Kepler-Umlauf.

**Prüfbare Folge:** Die Phase Re τ_M mod 1 = κ/ω_orb − 1 ist der Überschuss der radialen Phase pro azimutalem Umlauf. Im flachen Bereich einer Galaxie sind das 0,414·360° ≈ 149° radiale Phase pro Umlauf. Gleichwertig liegen aufeinanderfolgende Perizentren um 2π·ω_orb/κ = 360°/√2 ≈ 254,6° auseinander. Das ist die bekannte Rosettenbewegung von Sternbahnen. Sie folgt hier aus demselben τ_M, das die Bindung trägt.
