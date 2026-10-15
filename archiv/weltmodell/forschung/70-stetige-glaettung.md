# Stetige Glättung

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: nur Umbenennungen, keine inhaltliche Korrektur: Poisson-Tiefe w, w_J, w_b → ϖ, ϖ_J, ϖ_b; Stufe θ → 𝟙; Stufe χ_i → 𝟙_i.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/glaettung_stetig.py` (Ausgabe `glaettung_stetig.txt`)
> a₀ = 1,2·10⁻¹⁰ m/s². [A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft. Sternsysteme: Plummer mit b = R_e [A], Masse mit Υ = 2 [A], isotrope Jeans-Gleichung, σ im Fenster R < R_e, Feldform im homogenen Feld des Wirts (Punktmasse; Milchstraße 6·10¹⁰, NGC 1052 10¹¹ M☉ [A]).

## 1 Stetiger Einsatz

Die Plummer-Wolke der Hierarchie hat die Zentraldichte 3M/(4πℓ³). Mit ℓ³ = GM/T ist das genau die Gezeitendichte ρ_J = 3T/(4πG), T = √(3/2)|⟨Ē^TF⟩|: die Elterngezeit in Dichteeinheiten geschrieben, kein Stoff. Die Glättung senkt also jeden Körper auf die Dichte ab, bei der Eigenfeld und Elterngezeit gleich sind. Ein Körper, dessen größte Dichte schon darunter liegt, hat nichts abzugeben. Daraus folgt die stetige Regel:

```
ρ̄ = min(ρ, ρ_J) + Σ_i M_i^c·K_{ℓ_i^c},     M_i^c = ∫_{V_i}(ρ − ρ_J) dV,     (ℓ_i^c)³ = G·M_i^c/T_i,     V_i Komponente von {ρ > ρ_J}
```

- Statt g(R) > a₀ gilt ρ > ρ_J, also Eigenfeld gegen Elterngezeit; a₀ kommt im Kriterium nicht vor. Die Regel hat keine freie Zahl.
- Punktmasse: M^c = M, also die bisherige Hierarchie. Für die Sonne ist M^c/M − 1 = −1,4·10⁻²³.
- Die Glättung setzt bei ρ_max = ρ_J mit M^c = 0 ein. Überschussanteil eines Plummer-Körpers (Skala b): 0,85 / 0,51 / 0,24 / 0,07 / 0,005 / 0 bei b/ℓ = 0,1 / 0,3 / 0,5 / 0,7 / 0,9 / 1, am Einsatz ≈ 1,6·(1 − b/ℓ)^{5/2}.
- Energie einer Wolke von 1 M☉ am Sonnenort: linear im äußeren Feld ist V_s = −c_e·W mit W = (1/8πG)∫|∇φ̄|² und c_e = ν_e(1 + K/3) − 1 = 0,367; die volle Rechnung stimmt auf 0,1 m²/s². Die scharfe Regel springt beim Einsatz um ΔV/M = 517 m²/s² (homogene Wolke). Die stetige Regel gibt 0 bei ρ₀ = ρ_J und −0,4 / −3 / +4 / +170 m²/s² bei ρ₀/ρ_J = 1,001 / 1,01 / 1,1 / 1,5.
- Poisson-Form (Tiefe ϖ, e^{−|k|ϖ}, P_a·P_b = P_{a+b}): Die Punktmasse in Tiefe ϖ hat die Zentraldichte M/(π²ϖ³). Geglättet wird bis ϖ_J = (4/3π)^{1/3}·ℓ = 0,75·ℓ; ein Körper der eigenen Tiefe ϖ_b erhält (ϖ_J − ϖ_b)_+. [korr.] Der Einsatzpunkt ist derselbe.
- Das Modell legt fest, wo die Glättung einsetzt und wo sie voll ist. Nicht festgelegt ist der Verlauf dazwischen (Überschuss, umschlossener Anteil (1 − b²/ℓ²)^{3/2}, Poisson-Tiefe) und die Kernform (Plummer in 𝒢, Poisson auf der Halbgeraden). Das ist eine freie Funktion, keine freie Zahl.

## 2 Welche Dichte

Die Regel vergleicht eine Dichte mit ρ_J. Dafür gibt es zwei Lesarten:
- **(i) Dichte der Quelle:** Sterne sind Punkte, zwischen ihnen ist ρ = 0. Körper sind nur Sterne und dichte Wolken; jedes Sternsystem bleibt kollektiv. Das ist die bisherige Regel mit stetigem Einsatz.
- **(ii) Dichte der Stufe:** Die Sterne werden zuerst geglättet, dann wird ihr System mit der Gezeitendichte seines Wirts verglichen. Körper ist dann jedes System, das dichter ist als seine Gezeitendichte.

## 3 Prüfung

**Sonne** (beide Lesarten): ℓ = 0,934 pc = 27,4 r*, Quadrupol aus der Glättung 9,2·10⁻³¹ s⁻², γ − 1 = 5·10⁻⁵ bei 10 kAE und 1,3·10⁻³ bei 30 kAE.
**Weite Doppelsterne** (2 × 0,75 M☉ [A]): Der über seine Zone geglättete Begleiter trägt zur Gezeit nur ∝ (d/ℓ)² bei. Der Fixpunkt ist eindeutig, ℓ = 0,847–0,849 pc bei 5–30 kAE (ohne Begleiter 0,848 pc), und γ − 1 ≤ 9·10⁻⁴ bis 30 kAE: Newton.

| System | g_i/a₀ | g_e/a₀ | ρ̄/ρ_J | M^c/M | Newton | (i) Feldform | (ii) stetig | gemessen [G] |
|---|---|---|---|---|---|---|---|---|
| kompakter Haufen [A] (2·10⁵ M☉, 3 pc, 10 kpc) | 7,6 | 0,70 | 9·10³ | 0,96 | 6,03 | 6,25 | 6,03 | – |
| Pal 14 (2·10⁴ M☉, 25 pc, 71 kpc) | 0,011 | 0,014 | 594 | 0,90 | 0,66 | 1,81 | 0,66 | 0,38 ± 0,12 |
| NGC 2419 (9·10⁵ M☉, 19 pc, 90 kpc) | 0,85 | 0,009 | 1·10⁵ | 0,99 | 5,09 | 6,57 | 5,09 | 4,1 ± 0,5 |
| Crater II | 1·10⁻⁴ | 0,005 | 0,5 | 0,05 | 0,40 | 1,42 | 1,13 | 2,7 ± 0,3 |
| Draco | 0,004 | 0,012 | 30 | 0,68 | 1,20 | 3,40 | 1,03 | 9,1 ± 1,2 |
| Fornax | 0,027 | 0,003 | 444 | 0,88 | 5,55 | 14,25 | 5,46 | 11,7 ± 0,9 |
| DF2 | 0,014 | 0,018 | 7,2 | 0,49 | 7,04 | 18,06 | 3,71 | 8,5 +2,3/−3,1 |
| DF4 | 0,020 | 0,004 | 123 | 0,81 | 7,15 | 19,60 | 6,77 | 4,2 +4,4/−2,2 |

σ in km/s. ρ̄/ρ_J ist die mittlere Dichte in r_½ gegen die Gezeitendichte der baryonischen Elterngezeit (Hesse-Form von σ̄, wie im Modell).

- **(i)** ändert kein Sternsystem. Pal 14 liegt bei 1,81 statt 0,38, DF2 und DF4 bei 18 und 20 km/s, der kompakte Haufen 4 % über Newton.
- **(ii)** bringt Pal 14, NGC 2419 und DF4 auf Newton. Crater II füllt seine Gezeitenzone über (ρ̄/ρ_J = 0,5) und wird nur wenig geglättet (M^c/M = 0,05). Zugleich werden aber Fornax und Draco Newtonsch, 5,5 statt 11,7 und 1,0 statt 9,1 km/s. Draco und Crater II liegen schon in der Feldform zu tief.
- **Kraftgesetz:** Aus der Energie folgt 𝒢†f = (1 − 𝟙)f + 𝟙·[⟨f⟩_K + (ℓ/3)∂_ℓ⟨f⟩_K], mit 𝟙 als Stufe {ρ > ρ_J}. [korr.] Am Rand r_c von V springt das gespürte Potential um Δ = Φ_s(r_c) − Φ̃. Bei Sternen ist das belanglos. Bei Zwergen ist der Sprung so groß wie die ganze Struktur (DF2 −373, DF4 −362, Fornax −163 (km/s)²) und zieht nach außen: DF2 fällt auf 3,7 km/s, unter Newton (ohne Randterm wären es 8,3). Die glatte Form (ρ − ρ_J)²_+/ρ mit dem Kraftgewicht 1 − (ρ_J/ρ)² hat keine freie Zahl, verteilt den Term aber nur; DF2 hat dann kein Gleichgewicht. Die Werte 9,1 / 7,3 / 1,36 km/s (DF2 / DF4 / Crater II) entstehen nur, wenn jeder Anteil seine Kraft getrennt spürt; das folgt nicht aus der Energie.
- **Gezeit des verstärkten Wirtsfelds** statt der baryonischen: Die Überschüsse werden kleiner (DF2 0,18, Crater II 0), DF2 kommt auf 4,8 km/s, Fornax auf 4,2 km/s; DF4 und Draco haben kein Gleichgewicht.
- **Konsequent angewandt** sind auch die Wirte Körper. Die Milchstraße (Plummer 3 kpc [A]) im Feld von M31 (10¹¹ M☉ [A], 785 kpc [G]) hat ℓ = 459 kpc und M^c/M = 0,994. Ihre Kreisbahn wird Newtonsch: 112 statt 193 km/s bei 20 kpc, 72 statt 183 km/s bei 50 kpc. Alle Satelliten liegen in dieser Zone.

## 4 Warum keine Größe DF4 von Fornax trennt

Fornax ist in jeder Kenngröße des Modells kompakter und stärker gegen die Elterngezeit gebunden als DF4: g_i/a₀ 0,027 gegen 0,020, ρ̄/ρ_J 444 gegen 123, r½/r* 4,3 gegen 5,0; das äußere Feld ist ähnlich (g_e/a₀ 0,003 gegen 0,004). Verlangt wird aber σ²/σ²_N = 4,5 bei Fornax und höchstens 1,45 bei DF4 (oberer Rand); die Feldform gibt 6,6 und 7,5. Jede Regel, die mit Kompaktheit oder Bindung wächst, glättet Fornax mindestens so stark wie DF4.

Ordnet man alle Systeme, so ist das Vorzeichen in g_i/a₀ (und in g_i/g_e) geordnet. Crater II und Draco (g_i < g_e/2) verlangen mehr Verstärkung als die Feldform ((σ_gem/σ_F)² = 3,6 und 7,2), alle Systeme ab g_i ≈ 0,011 a₀ weniger (Pal 14 0,04, DF4 0,05, DF2 0,22, NGC 2419 0,39, Fornax 0,67). Der Betrag ist in keiner Kenngröße geordnet, und die Schwelle für das Vorzeichen läge zwischen 0,004 und 0,011 a₀; das wäre eine freie Zahl.

## 5 Kovariante Fassung

```
𝒢[σ] = D²σ − (4πG/c²)·[(ρ_σ − ρ_J)_+ − Σ_i γ_i M_i^c 𝒦_{ℓ_i^c}(ρ_i)],     M_i^c = ∫_{V_i}(ρ_σ − ρ_J) dV,     ρ_J = 3√(3/2)|Ē^TF|/(4πG)
```

- Die Regel passt zur Nebenbedingung D²σ̄ = 𝒢[σ]. Es ändert sich nur der Quellterm: 𝟙_i·D²σ wird zum Überschuss. [korr.] M_i^c ist ein Volumenintegral skalarer Größen auf dem Ruheschnitt, und ρ_J hängt wie ℓ über Ē von σ̄ ab (Fixpunkt wie bisher). Der Newton-Grenzfall ist Abschnitt 1, der Einsatz ist stetig.
- Unstetig bleibt die Zerlegung in Komponenten. Verschmelzen zwei Kerne (je 1 M☉, Sattel erreicht ρ_J), springt die Glättung von zwei Plummer-Wolken (ℓ = 0,93 pc) auf eine (1,18 pc): ΔV/M = 206 m²/s², gegen 517 m²/s² beim scharfen Einsatz.
- Stetig in der Energie heißt nicht stetig in der Kraft; es bleibt der Randterm aus Abschnitt 3.
- Lesart (ii) verlangt 𝒢 auf der schon geglätteten Dichte, also D²σ̄ = 𝒢[σ̄] als Fixpunkt. Damit wird jedes gebundene System Körper.

## Ergebnis

Ein stetiger Einsatz folgt aus dem Modell ohne freie Zahl: Die Glättungswolke hat die Gezeitendichte als Zentraldichte, also wird nur der Überschuss über ρ_J geglättet. Der Energiesprung beim Einsatz verschwindet (517 → 0 m²/s²), Sonne, Quadrupol und weite Doppelsterne bleiben wie bisher, und die Regel passt zur kovarianten Nebenbedingung.

## Wo es nicht trägt

1. Mit der Dichte der Quelle ändert sich an keinem Sternsystem etwas; DF2 und DF4 bleiben bei 18 und 20 km/s, Pal 14 bei 1,8 km/s.
2. Mit der Dichte der Stufe trennt die Regel DF4, Pal 14 und NGC 2419 von Crater II. Dasselbe Kriterium (dichter als die Gezeitendichte) erfüllt aber jedes gebundene System: Fornax, Draco und, konsequent, die Kreisbahn der Milchstraße werden Newtonsch.
3. Das Kraftgesetz aus der Energie hat am Rand der geglätteten Zone einen Sprung in der Größe der ganzen Struktur. DF2 kommt auf 3,7 km/s, unter Newton; die glatte Form verschiebt den Term nur.
4. An dieser Stelle fehlt die Größe: Fornax ist in g_i/a₀, ρ̄/ρ_J und r½/r* kompakter als DF4, verlangt aber mindestens das Dreifache der Verstärkung. Keine Größe des Modells trennt beide. Der Unterschied muss in den Eingaben liegen (Υ, Abstand, Gleichgewicht) oder in einer Größe, die das Modell nicht enthält.
