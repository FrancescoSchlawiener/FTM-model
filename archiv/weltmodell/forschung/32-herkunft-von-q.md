# Herkunft von q = 0,843: die kosmische Rückkehrzeit

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Hypothese ν = 2π, Q = 2π − 1 → Zahlennähe ohne Herleitung (21 Kandidaten im 2σ-Band); als Kac-Zeit ν = 2ΛH/c hinge die 2π an Λ = πc/H, das nicht trägt, bei festem Λ liefe ν mit H, obwohl Ω_m/Ω_b fest ist; Q = 5,364 ist aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0); Vorhersage und Offen nur unter der Hypothese; „freie, nicht projizierte Materie“ → Baryonen, die Mode ξ ist eine Einschränkung der Freiheitsgrade, keine Materie. Umbenennung (ohne Marke): (∂χ)² = −1 → (∂ξ)² = −1. [G] ergänzt: Planck-Werte, Zielwert künftiger CMB-Messungen.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/rueckkehr_2pi.py`

## Umformulierung

Die Rückkopplungssumme ν = 1/(1−q) ist exakt die Kac-Rückkehrzeit. q ist der Rückkehranteil pro Umlauf. Kosmisch gilt mit Q = q/(1−q) = Ω_c/Ω_b:

```
1 − q = Ω_b/Ω_m          (Anteil der Baryonen; der Rest ist die Mode ξ, eine Einschränkung der Freiheitsgrade, keine Materie) [korr.]
ν_kosmisch = 1/(1−q) = Ω_m/Ω_b
```

Die Frage nach der Herkunft von 5,36 ist damit die Frage: **Wie lang ist die kosmische Rückkehrzeit der Baryonen?** [korr.]

## Hypothese: Die Rückkehrzeit ist ein voller Umlauf

[korr.] Diese Hypothese setzt voraus, dass die kosmische Kac-Rückkehrzeit genau ein Umlauf ist; das trägt nicht: Als Kac-Zeit ist ν = 2ΛH/c, und die 2π hinge an Λ = πc/H. Mit dem Hubble-Radius atmend bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); bei festem Λ liefe ν mit H, obwohl Ω_m/Ω_b fest ist. Im Modell ist Q = 5,364 aus der Messung eingesetzt; im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0).

Im geschlossenen Torus kehrt eine Schleife nach einem vollen Umlauf zurück, also nach 2π Radiant. Setzt man die mittlere Rückkehrzeit der Baryonen (Kac: 1/Anteil) gleich einem Umlauf, folgt: [korr.]

```
ν_kosmisch = Ω_m/Ω_b = 2π
q = 1 − 1/2π = 0,84085,      Ω_c/Ω_b = 2π − 1 = 5,2832,      Ω_b/Ω_m = 1/2π = 0,15915
```

| Messung | ν = Ω_m/Ω_b | Abstand zu 2π = 6,2832 |
|---|---|---|
| Planck 2018 (TT,TE,EE+lowE+lensing) | 6,364 ± 0,065 [G] | +1,26 σ |
| Planck 2018 + BAO | 6,322 ± 0,053 [G] | +0,75 σ |

Nur unter dieser Hypothese folgte: Ω_c h² = (2π−1)·Ω_b h² = 0,1182 bei Ω_b h² = 0,02237. [korr.]

## Einordnung

- **Zahlennähe:** Im 2σ-Band liegen 21 „natürliche" Kandidaten, zum Beispiel 19/3, 25/4, 32/5, 3π−3, e²−1 und φ³+2. 2π unterscheidet sich von ihnen nur durch die Deutung. Der Treffer allein ist kein Beleg. [korr.] Zahlennähe ohne Herleitung, trägt nicht.
- **Deutung:** Kac-Rückkehrzeit der Baryonen = ein Umlauf des Torus. [korr.] Das ist als Prinzip formuliert, aber nicht aus der Wirkung hergeleitet.
- **Was sie leisten würde:** Mit ν = 2π wäre das Modell im Kosmos **ohne freie Zahl für die Dunkle Materie**. Q wäre dann eine Folge der Schließung.

## Der Test

Der Unterschied zwischen 2π und dem Messwert beträgt 0,6–1,3 %. Eine Prüfung auf 5σ verlangt σ(Ω_m/Ω_b)/(Ω_m/Ω_b) ≤ 0,12–0,25 %. Heute sind es etwa 0,8 %. Die nächste Generation von CMB-Messungen zielt auf etwa 0,2–0,3 % je Dichteparameter [G]. Damit wird die Hypothese entscheidbar.

## Offen

Offen ist die Herkunft von Q = 5,364 selbst. Unter der Hypothese oben wäre aus der Wirkung der Schließungsregel herzuleiten, warum die Rückkehrzeit der Baryonen genau einen Umlauf beträgt; dafür müsste die Nebenbedingung (∂ξ)² = −1 des projizierten Flusses mit der Periodizität des kosmischen Torus verknüpft werden. Die Hypothese trägt aber nicht. [korr.]
