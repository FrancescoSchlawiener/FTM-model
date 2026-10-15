# Zwerge DF2 und DF4 im äußeren Feld

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Regel des Modells (Feldform …)“ → örtliche Regel a = ν(|g_N|)·g_N, Zahlen gleich (Monopol beider gleich); „Kriterium g(R) > a₀“ → geglättet wird nur der Überschuss (ρ − ρ_J)₊ der Quelldichte, a₀ kommt darin nicht vor; Umbenennung L = d ln ν/d ln g_N → k_ν.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/zwerge_aussenfeld.py` (Ausgabe `zwerge_aussenfeld.txt`)
> Feldform im homogenen äußeren Feld, achsensymmetrisch in (R, z) über die volle z-Achse, Finite Volumen 220 × 440 Zellen bis 500 kpc; Streuung aus der isotropen Jeans-Gleichung, σ_los im Fenster R < R_e. [A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.

## Annahmen

| | M (Υ = 2) [G] | R_e [G] | Nachbar [A] | Abstand, projiziert [G] | g_Ne/a₀ | σ gemessen [G] |
|---|---|---|---|---|---|---|
| DF2 | 2,0·10⁸ M☉ | 2,2 kpc | 10¹¹ M☉ | 80 kpc | 0,018 | 8,5 +2,3/−3,1 km/s |
| DF4 | 1,5·10⁸ M☉ | 1,6 kpc | 10¹¹ M☉ | 165 kpc | 0,0043 | 4,2 +4,4/−2,2 km/s |

Plummer-Kugel mit b = R_e [A]; Nachbar als Punktmasse, sein Feld über den Zwerg homogen.

## Woher die 18–20 km/s kommen

- Schätzer σ² = r_½·g(r_½)/3 mit r_½ = 4/3·R_e und M(<r_½) = M/2.
- Außenfeld: Newton-Feld des Nachbarn, projizierter Abstand als 3D-Abstand.
- Glättung: keine. Der Zwerg ist kollektiv (g_N(r_½) = 0,014 bzw. 0,020 a₀), geglättet werden nur seine Sterne (1 M☉: ℓ = 2–3 pc).
- Zusammensetzung: vektoriell, ν(|g_i + g_e|) über die Kugel gemittelt. Nach dem Satz von Gauß ist das genau der Monopol der Feldform, ⟨g_r⟩ = g_i + ⟨(ν(|g_N|) − 1)(g_i + g_e cos θ)⟩. Die 18,0 und 19,7 km/s sind also schon die Feldform; mit eigener Verstärkung bis zur Gezeitengrenze entstehen 20–21 km/s.

## Feldform im homogenen äußeren Feld

```
Φ = −ν_e·g_e·z + φ,     ∇²φ = ∇·[ν(|∇Φ_N|)∇Φ_N − ν_e∇Φ_N,e],     Φ_N = Φ_N,innen + g_e·z
fern:  φ = −ν_e·GM/r·[1 + (k_ν/2)·sin²θ],     k_ν = d ln ν/d ln g_N = −(ν_e − 1)·s_e/2     (DF2: ν_e = 7,93, k_ν = −0,47)
```

θ zählt gegen die Feldrichtung. Außen ist die Anziehung längs des Felds ν_e·GM/r², quer um den Faktor 1 + k_ν/2 schwächer, im Mittel ν_e(1 + k_ν/3)·GM/r². [korr.] Probe: isolierte Kugel gegen ν(g_N)·g_N auf 0,07 %, mit Außenfeld der Monopol auf dem Gitter gegen das Flussmittel auf 0,07 %.

## Vergleich im projizierten Abstand (km/s)

| | Newton | isoliert | gleichgerichtet s² = s_i² + s_e² | senkrecht s⁴ = s_i⁴ + s_e⁴ | Feldform | Feldform, Blick längs / quer zum Feld |
|---|---|---|---|---|---|---|
| DF2 | 7,0 | 21,1 | 17,3 | 18,7 | 18,1 | 17,0 / 18,6 |
| DF4 | 7,2 | 19,8 | 18,8 | 19,6 | 19,6 | 19,2 / 19,8 |

Die Feldform liegt zwischen den algebraischen Formen, bei schwachem Außenfeld auf der senkrechten. Schätzer und Jeans-Rechnung unterscheiden sich um weniger als 1 %. Die Blickrichtung folgt aus dem Tensor-Virialsatz, σ_zz²/σ_xx² = 0,84 (DF2) und 0,94 (DF4): längs des Felds ist die Streuung kleiner.

## Abstand zum Nachbarn

Der 3D-Abstand ist mindestens der projizierte; jeder größere Abstand schwächt das Außenfeld und hebt σ (Feldform, km/s):

| D/D_proj | 1 | 1,5 | 2 | 3 |
|---|---|---|---|---|
| DF2 | 18,1 | 20,2 | 20,8 | 21,0 |
| DF4 | 19,6 | 19,7 | 19,8 | 19,8 |

Der projizierte Abstand ist der günstigste Fall. Selbst ein Feld vom Betrag des verstärkten Nachbarfelds (ν·g_Ne = 0,144 a₀) im Argument gäbe für DF2 erst 11,6 km/s.

## Glättungsradius

Wird die Masse des ganzen Zwergs über ℓ geglättet (im Argument, in der Quelle und in der gespürten Struktur, wie bei einem kompakten Trichter), sinkt σ zum Newton-Wert (Feldform, km/s):

| ℓ | 0 | ½ r_½ | r_½ | ℓ_J | 2ℓ_J | ∞ |
|---|---|---|---|---|---|---|
| DF2 (r_½ = 2,87 kpc, ℓ_J = 3,97 kpc) | 18,1 | 12,3 | 9,4 | 8,4 | 7,3 | 7,0 |
| DF4 (r_½ = 2,09 kpc, ℓ_J = 5,94 kpc) | 19,6 | 14,8 | 11,4 | 7,7 | 7,3 | 7,2 |

ℓ_J = [GM/(√(3/2)·|E^TF|)]^{1/3} im verstärkten Feld des Nachbarn.

## Ergebnis

Mit der Regel des Modells (örtliche Regel a = ν(|g_N|)·g_N, Zwerg kollektiv, ℓ = 0; im Monopol gleich der Feldform, die Zahlen bleiben) [korr.] liegt die Streuung bei 18,1 (DF2) und 19,6 km/s (DF4), auch im günstigsten Abstand und in jeder Blickrichtung mindestens 17,0 km/s, also nicht im gemessenen Bereich. Für den oberen Rand der Messung (10,8 bzw. 8,6 km/s) fehlt eines davon:

1. **Außenfeld:** g_Ne ≥ 0,21 a₀ im Argument für DF2, das 12-Fache des Newton-Felds des Nachbarn (mit 10¹¹ M☉ ein Abstand ≤ 23 kpc bei 80 kpc projiziert, oder bei 80 kpc ≥ 1,2·10¹² M☉ baryonisch); für DF4 ≥ 0,96 a₀, das 224-Fache.
2. **Masse:** M ≤ 7·10⁷ M☉ (Υ ≤ 0,7) für DF2, ≤ 1,7·10⁷ M☉ (Υ ≤ 0,22) für DF4. Der Mittelwert von DF4 (4,2 km/s) liegt schon unter Newton mit Υ = 2 (7,2 km/s) und verlangt dort Υ ≤ 0,7.
3. **Glättung des ganzen Zwergs:** ℓ ≥ 2,0 kpc (DF2) bzw. 4,0 kpc (DF4), für den Mittelwert von DF2 3,8 kpc. Mit ℓ = ℓ_J aus der Formel der Hierarchie lägen beide im gemessenen Bereich (8,4 und 7,7 km/s). Dieselbe Glättung senkt Crater II an der Milchstraße [G: 3,2·10⁵ M☉, R_e = 1,07 kpc, 117,5 kpc, σ = 2,7 ± 0,3] von 1,42 auf 0,93 km/s, also weiter vom Messwert weg. Nach der Regel der Hierarchie wird aber nur der Überschuss (ρ − ρ_J)₊ der Quelldichte geglättet, ρ_J = 3T/(4πG), a₀ kommt darin nicht vor; Sterne sind Einzeltrichter, Sternsysteme kollektiv, die Zwerge werden also nicht als Ganzes geglättet. [korr.]

Crater II verlangt ein schwächeres, DF2 und DF4 verlangen ein stärkeres Außenfeld, als die Feldform gibt; keine Änderung allein der Stärke des Außenfelds trifft beide. Es fehlt eine Größe des Modells, die bei DF2 (g_i ≈ g_e) und DF4 (g_i > g_e) die wirksame Verstärkung σ²/σ²_Newton auf höchstens 2,4 bzw. 1,5 senkt (Feldform: 6,6 bzw. 7,5) und bei Crater II (g_i ≪ g_e) nicht. Unsicher bleiben die Eingaben Masse (Υ), 3D-Abstand und Gleichgewicht der Zwerge.
