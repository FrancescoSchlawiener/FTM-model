# Freiheitsgrad-Einschränkung statt Zusatzquelle

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: w = w(θ_u) als Momentanwert (w → 0 bei θ_u = 0) → das Gewicht hängt am tiefsten Kollaps der Geschichte (Preisach) und bleibt; die Übergabe der Mode ξ in die Operatorform läuft mit λ_m = ω_orb/N*.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/freiheitsgrad.py`

---

## Kurzfassung

1. **Was im CMB gebraucht wird, in Freiheitsgraden gesagt.** Ein druckfreier, adiabatisch an die Baryonen gebundener **Impuls-Freiheitsgrad** mit dem Gewicht 5,36·ρ_b.
   - Seine Druck- und Streurichtung ist inaktiv. Er schwingt deshalb nicht mit und verschiebt die Nulllage der akustischen Schwingung.
   - Seine Fallrichtung ist aktiv. Er fällt in die Mulden und vertieft sie.
2. **Die Projektion des Autors trägt das CMB.** ż = (I − w·Π⊥)F, wobei Π⊥ die Fluchtrichtung aus der Mulde ist (Druck bzw. Streuung). Bei w = 1 ergibt sich das gemessene Peakmuster exakt. Die Zahl 5,36 erscheint als **Gewicht des projizierten Flusses**, nicht als Stoffmenge.
3. **In Galaxien** liefert der reine Operator mit stationärem Gedächtnis die Radial-Acceleration-Relation exakt: g = g_N/(1−q), ohne jeden Überschuss.
4. **Die Einheit beider Regime** verlangt ein zustandsabhängiges Projektionsgewicht:
   - w = 1 im atmenden Hintergrund (CMB, großräumiges Wachstum),
   - w → 0 in gebundenen, nicht mehr mitatmenden Systemen (Galaxien).

---

## 1. Die Rolle im CMB als Aussage über Freiheitsgrade

Das Photon-Baryon-Fluid ist ein Oszillator: Druck treibt die Baryonen aus der Mulde, Gravitation zieht sie zurück. Gebraucht werden zwei Freiheitsgrade pro Fourier-Mode, (δ_s, θ_s), mit drei Eigenschaften:

| | Eigenschaft | Wirkung im CMB |
|---|---|---|
| F1 | Druck- und Streurichtung **inaktiv** | Der Potentialanteil schwingt nicht mit; die Nulllage der Schwingung verschiebt sich (Baryon-Loading). |
| F2 | eigene Geschwindigkeit, **fällt in die Mulde** | δ_s wächst (logarithmisch im Strahlungs-, linear im Materie-Zeitalter); die Mulde wird tiefer und bleibt stehen. |
| F3 | adiabatisch an die Baryonen gebunden | keine Isokrümmung |

**Analytisch (fester Ψ):** Θ₀+Ψ = (1/3+R)Ψ·cos(k r_s) − RΨ. Mit R(z*) = 0,6219 ist das Amplitudenverhältnis ungerader zu gerader Peaks 1 + 6R = 4,731. Nur ein Potentialanteil, der nicht zur akustischen Mode gehört, verschiebt die Nulllage.

**Gebundenes gegen freies Potential bei z*** (Poisson-Quelle in Einheiten Ω_b·δ):

| k [1/Mpc] | 0,02 | 0,05 | 0,10 | 0,20 |
|---|---|---|---|---|
| gebunden Q·δ_s | −45,96 | −132,13 | −230,07 | −349,95 |
| Hüllkurve \|δ_b\| | 4,68 | 3,85 | 4,11 | 4,05 |
| Verhältnis | 9,83 | 34,36 | 55,99 | 86,32 |

5,36 ist das Verhältnis der Hintergrunddichten, also das Gewicht des Freiheitsgrads. Seine Wirkung bei z* entsteht durch das Hineinfallen (F2).

## 2. Die Projektion (I − w·Π⊥) auf den Baryonenfluss

Im Sinn des Autors ist Π⊥ die Projektion auf die **Fluchtrichtung**. Für ein Baryon in der Mulde ist das die Druck- und Streukraft. Der projizierte Fluss ż = (I − w·Π⊥)F behält Gravitation und Trägheit, die Druckrichtung ist inaktiv.
- Er hat **keine eigenen Anfangsdaten**: Er startet im Zustand der Baryonen und ist damit von selbst adiabatisch (F3).
- In der Poisson-Gleichung zählt er mit dem Gewicht Q = 5,364.
- Der Hintergrund ist Ω_eff = Ω_b(1+Q) ∝ a⁻³, also von selbst staubartig.

Verhältnis 3. zu 2. Peak (Spielzeug-Boltzmann, relative Änderung gegen den Fall w = 1):

| w | 0,75 | 0,9 | 0,97 | 1 |
|---|---|---|---|---|
| Abweichung | 18,7 % | 8,2 % | 2,6 % | 0 |

Für ≤ 10 % Abweichung genügt w ≥ 0,9, für ≤ 3 % w ≥ 0,97.

**Formale Einordnung:** Eliminiert man den projizierten Fluss, bleibt für die Baryonen ein Gedächtnisterm. Dessen Kern ist die Green-Funktion des druckfreien Hineinfallens, und sie **wächst** mit der Zeit (∝ ln a bzw. ∝ a). Projektion und Gedächtnis sind zwei Seiten desselben Formalismus (Zwanzig 1961, Mori 1965). Der Kern, der das CMB trägt, ist ein wachsender Kern, kein abklingender.

**Trennbarkeit:** Linear sind die Gleichungen des projizierten Flusses mit w = 1 identisch mit denen einer druckfreien Staubkomponente der Dichte Q·ρ_b. Nach der Trennlogik des Autors (Diagramm 6) trennen sich die beiden Lesarten erst dort, wo w vom Zustand abhängt, also in gebundenen Systemen (Trennklassen 1 und 3).

## 3. Welche Rückkopplungsstärke entspricht 5,36?

Gebundener zu freiem Anteil ist Σ_{n≥1} qⁿ = q/(1−q) = 5,364. Daraus folgen

```
q = 0,84287,    ν = 1/(1−q) = 6,3643,    y = (ln q)² = 0,02922
```

Offene Frage: Welche Größe des Modells legt bei z ≈ 1100 genau diesen Wert fest?

## 4. Galaxien

Baryonen als Hernquist-Profil, a₀ = 1,2·10⁻¹⁰ m/s². Verhältnis g/g_RAR für den reinen Operator (stationäres Gedächtnis, q lokal):

| r [kpc] | 5 | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|
| M_b = 6·10¹⁰ M☉ | 1,0000 | 1,0000 | 1,0000 | 1,0000 | 1,0000 |
| M_b = 10⁹ M☉ | 1,0000 | 1,0000 | 1,0000 | 1,0000 | 1,0000 |

Das stationäre Gedächtnis konvergiert gegen g = g_N/(1−q) = ν·g_N, per Iteration geprüft. q hängt dabei lokal von y = g_N/a₀ ab.

## 5. Die Schließungsregel für beide Regime

Die Projektion trägt den Kosmos, der reine Operator die Galaxien. Eine einzige Regel vereint beide, wenn das Projektionsgewicht vom Zustand abhängt:

```
w = w(θ_u),   θ_u = ∇_μ u^μ  (Expansionsrate des Zeitfelds, das Atmen)
w = 1       im mitatmenden Hintergrund       → CMB und großräumiges Wachstum wie gemessen
w → 0       in gebundenen Systemen (θ_u = 0)   → reiner Operator, g = g_N/(1−q)
```

[korr.] Das Gewicht ist kein Momentanwert von θ_u: Im Modell ist w = θ₊²/(θ₊² + θ_c²) mit θ₊ = max(θ, 0) und θ_c = a₀/c, und es hängt am tiefsten Kollaps der Geschichte (Preisach), nicht am heutigen Zustand; einmal gefallen, bleibt es. Die Mode ξ geht danach mit der Übergaberate λ_m = ω_orb/N* (N* ≈ 0,3–2,4) in die Operatorform über.

Lesart im Modell: **Solange ein Bereich mit dem Torus atmet, ist seine Fluchtrichtung vollständig geschlossen. Sobald er sich vom Atmen entkoppelt und gebunden ist, bleibt nur die stationäre Rückkopplung.**

## 6. Nächste Rechnungen

1. **w(θ_u) als Wirkung** formulieren. Prüfen: Energie-Impuls-Erhaltung (00- und 0i-Gleichung zugleich) und Stabilität.
2. **Boltzmann-Code (CLASS/hi_class)** mit Projektionsmodul w(k, a). Der Grenzfall w = 1 muss das gemessene Spektrum exakt ergeben. Dann fσ₈(k, z) und das Linsenpotential mit dem Übergang w → 0 bei der Entkopplung rechnen.
3. **N-Körper-Test** des Übergangs vom projizierten Fluss zum Operator ν in einer kollabierenden Galaxie. Dabei prüfen, ob Galaxienhaufen und der Bullet Cluster ihre Linsenmasse behalten.

---

## Literatur

- Hu W., Sugiyama N. 1995, ApJ 444, 489.
- Ma C.-P., Bertschinger E. 1995, ApJ 455, 7.
- Hu W., Fukugita M., Zaldarriaga M., Tegmark M. 2001, ApJ 549, 669.
- Planck Collaboration 2020, A&A 641, A6.
- Zwanzig R. 1961, Phys. Rev. 124, 983; Mori H. 1965, Prog. Theor. Phys. 33, 423.
- Autor: „Zustandsraum Term, Toy model" und Diagramme 5 und 6 des Zustandsraum-Projekts.
