# Halbgerade der Hierarchie

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Abschnitt 2 setzt Λ als an den Hubble-Radius gebundene Raumlänge (πc/H = U/2) voraus → als nicht tragend gekennzeichnet, Λ ist die ganze innere Faltung der Zelle, heute fest, keine Raumlänge; Zeile „πc/H = U/2 = 12,1 Gpc“ → c²/(2a₀) = ℓ₀/2; Zeile „13,5 Gpc (U ≥ 27 Gpc)“ → 13,5 Gpc (ℓ₀ = 27 Gpc, keine Raumzelle; die Schranke aus den Kreispaaren gilt für die Raumzelle L); „Hierarchie schließt erst an der Hubble-Kugel (L_H ≈ Λ)“ → schließt bei Λ, L_H ≈ Λ nur heute; Zelle U → ℓ₀ (c²/a₀) bzw. L (Raumzelle); innere Richtung w → ϖ. Der Kopf setzt H = 2πa₀/c, also a₀ = cH/2π (mit a₀ = 1,2·10⁻¹⁰ ist das H = 77,6 km/s/Mpc).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/halbgerade_hierarchie.py` (Ausgabe `halbgerade_hierarchie.txt`)
> a₀ = 1,2·10⁻¹⁰ m/s², H = 2πa₀/c. Beispiele: Milchstraße 6·10¹⁰ M☉ (r* = 8,35 kpc), Zwerg 10⁸ M☉ (r* = 0,34 kpc).

## Frage

q = e^{−s} ist G(L) von D·G'' = Γ·G auf ϖ ≥ 0 (G(0) = 1, G beschränkt). Ist ϖ die Richtung der Hierarchie der Faltungen, und folgt daraus die Halbgerade statt einer Schleife?

## 1 Irrfahrt auf den Ebenen

Kette n = 0, 1, 2, … (0 = Ort), Sprungrate κ je Richtung, Verlust Γ je Zeit, G_n = E_n[e^{−ΓT₀}]:

```
κ(G_{n+1} − 2G_n + G_{n−1}) = Γ·G_n  ⇒  G_n = λⁿ,   q_m = G_m = exp(−2m·arsinh(s/2m)) → e^{−s}
τ_D = m²/κ,   D = κh²,   L = m·h
```

- **Halbgerade:** Sie entsteht aus drei Bedingungen: Der Ort ist das untere Ende, Auf- und Abwärtsschritte sind gleich häufig, die oberste Ebene ist nicht mit der untersten verbunden. Ob der Boden absorbiert (Ankunft = Rückkehr) oder reflektiert (Quelle am Ort), ist gleich; das normierte Profil stimmt auf 10⁻¹² überein.
- **Wenige Ebenen:** m = 1 / 3 / 10 Ebenen treffen ν schon auf 4,7 / 0,63 / 0,06 %. Der Newton-Schwanz ist aber (m/s)^{2m} statt e^{−s} (bei s = 10 ×216 für m = 1); e^{−s} bis s verlangt m ≫ s.
- **„Nur aufwärts“ als Bewegung** trägt nicht: Zurück führt dann nur eine geschlossene Kette, 1 − q → s², a → a₀ (Boden, v² = a₀r). Drift aufwärts begrenzt ν (Pe = 0,05: ν ≤ 10,5, außen Kepler), abwärts gibt den Boden 2Pe·a₀. „Nur aufwärts“ ist die Ordnung der Ebenen, der Ort liegt unten; die Bewegung ist symmetrisch.
- **Keine Schleife:** Als Glättungslänge ist ϖ additiv. Die harmonische Fortsetzung (∇² + ∂_ϖ²)u = 0 ist die Poisson-Glättung mit Symbol e^{−|k|ϖ}: P_a·P_b = P_{a+b}, P_0 = 1, beschränkt nur für ϖ ≥ 0 (ϖ < 0 wäre Entglättung, e^{+|k|ϖ}). Eine Schleife verlangte P_Λ = 1 für ein Λ > 0; das erfüllt nur k = 0, der Mittelwert.
- **Abschirmzahl:** Für die Punktmasse ist D = L²/τ_D = GM/c = 3Δs·c und Γ = a_N/c, also √(Γ/D) = 1/r genau, die Wellenzahl der Kreisbahnmode.
- **L = r\*:** L ist die Glättungslänge, bei der die geglättete Punktmasse die Schwelle erreicht. Mit der Oberflächenneigung GM/L² des Kriteriums (gleichmäßige Kugel) ist L = r* genau. Die Spitze der Plummer- bzw. Poisson-Glättung (α = 0,385 bzw. 0,186) gäbe L = 0,62 bzw. 0,43·r* und verschöbe Tully-Fisher um ×2,6 bzw. ×5,4.
- **Zahl der Ebenen:** Läuft der Schritt mit c, ist D = hc/2, also h = 2GM/c² = 6Δs = r₀ und m = r*/r₀: Milchstraße 1,5·10⁶ (r₀ = 1185 AE), Zwerg 3,6·10⁷. Die Kette ist im Grenzfall (Diskretheit s³/24m² ≈ 10⁻¹⁴), der Weg in τ_D ist c·τ_D = ℓ₀. Galaxie, Gruppe und Haufen sind keine Schritte, sondern Stellen, an denen D wechseln kann (3b).

## 2 Wo die Hierarchie schließt

[korr.] Diese Rechnung setzt Λ als Raumlänge voraus, an den Hubble-Radius gebunden (Λ = πc/H als halbe Raumzelle, c²/a₀ als die Raumzelle); das trägt nicht: Mit dem Hubble-Radius atmend bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge; die räumliche Zelle folgt aus a₀ = 3c⁴/(2Gρ_mL³), bei a₀ = 1,2·10⁻¹⁰ ist L = 26,8 Gpc. Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute. Die Werte für r₁₀ = 0,657·Λ gelten für jedes feste Λ.

Oben liegt der geschlossene Raum; über ihm gibt es keine Ebene, der Deckel reflektiert. Bei durchlässigen Ebenen (D gleich):

```
q = cosh(k(Λ − L))/cosh(kΛ),  k = s/L;     tief: 1 − q ≈ s·tanh(Λ/r)  ⇒  10 % bei r₁₀ = Λ/artanh(1/1,1) = 0,657·Λ
```

| Λ | Λ/L Milchstraße | Λ/L Zwerg | r₁₀ (beide) | a_N/a₀ bei r₁₀ (MS) |
|---|---|---|---|---|
| c/H = 3,86 Gpc | 4,6·10⁵ | 1,1·10⁷ | 2,5 Gpc | 1·10⁻¹¹ |
| c²/(2a₀) = ℓ₀/2 = 12,1 Gpc [korr.] | 1,5·10⁶ | 3,6·10⁷ | 8,0 Gpc | 1·10⁻¹² |
| 13,5 Gpc (ℓ₀ = 27 Gpc, keine Raumzelle) [korr.] | 1,6·10⁶ | 4,0·10⁷ | 8,9 Gpc | 9·10⁻¹³ |

- Weil L = r*, hängt r₁₀ nicht von der Masse ab. Die Abweichung liegt bei Gpc, weit hinter jeder Rotationskurve: mit flachen Kurven verträglich. Die mittlere Rückkehrzeit ist endlich (Milchstraße, Λ = c/H: 3,7·10¹⁶ Jahre).
- **Spitze:** Die Hubble-Kugel (GM_H = c³/2H, Randneigung cH/2) hat L_H = √π·c/H, also L_H/Λ = 0,56–1,8. Für das Ganze ist der Rückweg geschlossen (ν_H = 1,220 statt 1,205 bei Λ = πc/H), für jedes Teilsystem ist L/Λ ≤ 1,2·10⁻⁴ (Haufen). Die Schließung liegt für jedes Teilsystem 10⁴–10⁷·L höher; das ist die Halbgerade. [korr.] L_H ≈ Λ gilt nur heute: Λ ist keine Raumlänge und folgt dem Hubble-Radius nicht, sondern ist die ganze innere Faltung der Zelle und heute fest; die Hubble-Kugel ist kein Deckel der Hierarchie.

## 3 Umgebung

**(a) Im Verlust.** Die Regel liest Γ = a_N,ges/c am Ort. Einbettung hebt Γ; L = r* (eigene Masse) und Λ (Spitze) bleiben. Das ist der Effekt des äußeren Feldes.

**(b) An der Ebenengrenze.** Sind die Faltungen je Ebene additiv, D_Ebene = 3c·ΣΔs = GM_Ebene/c, springt D oberhalb der eigenen Glättungslänge ℓ = d·(M/3M_p)^{1/3} um μ = M_p/M. Bei Schritten mit c ist der Fluss D·G′ stetig:

```
R = (1 − √μ)/(1 + √μ),   q = e^{−kL}·(1 + R·e^{−2k(ℓ−L)})/(1 + R·e^{−2kℓ});   tief: ν_b/ν = (1 + R·e^{−2ℓ/r})/(1 − R·e^{−2ℓ/r})
```

Für μ ≫ 1 ist R → −1: Die Grenze wirkt fast absorbierend, die Kurve fällt dort. 10 % bei r₁₀ = 0,67–0,69·ℓ; weit außen (r ≫ ℓ√μ) folgt eine tiefere Ebene v_f·μ^{−1/4}. Beispiele (d angenommen): Milchstraße mit g_e = 0,03a₀ bei d = 1 Mpc (μ = 430, ℓ = 92 kpc = 11 r*), Zwerg bei 100 kpc von 6·10¹⁰ M☉ (μ = 600, ℓ = 8,2 kpc = 24 r*). a/a_frei:

| | r/ℓ | äußeres Feld ⊥ / ∥ | Ebenengrenze | beide (⊥) |
|---|---|---|---|---|
| Milchstraße | 0,5 / 0,7 / 1 | 0,87 / 0,48 · 0,72 / 0,38 · 0,54 / 0,28 | 0,96 · 0,90 · 0,77 | 0,86 · 0,70 · 0,52 |
| Zwerg | 0,5 / 0,7 / 1 | 0,85 / 0,44 · 0,68 / 0,34 · 0,50 / 0,25 | 0,97 · 0,90 · 0,78 | 0,83 · 0,66 · 0,48 |

- Die Grenze wirkt schwächer als das äußere Feld und käme hinzu (1–2,5 Prozentpunkte).
- Liegt die Ebene des Sterns (ℓ = 0,93 pc) unter der der Galaxie, sperrt dieselbe Grenze die Rückkehr (bei 8 kpc q = 10⁻¹⁸ statt 0,35). Der Ort ist die Ebene des kollektiven, über die Zonen geglätteten Feldes.

## 4 Prüfbare Folge

Nur bei (b). Milchstraße bei 50 kpc, gleiches g_e = 0,03a₀: d = 0,5 / 1 / 2 Mpc gibt ℓ = 73 / 92 / 116 kpc und aus der Grenze a/a_frei = 0,91 / 0,95 / 0,98; das äußere Feld gibt in allen drei Fällen 0,84 (⊥). Bei gleichem äußeren Feld hinge die Außenkurve am Abstand zur Elternmasse. Bei durchlässigen Ebenen gibt es keine prüfbare Folge.

## Ergebnis

Die Hierarchie gibt die Halbgerade: Der Ort ist das untere Ende, die Glättungslänge ist additiv und nur aufwärts beschränkt, eine Schleife ist ausgeschlossen. Die Bewegung auf den Ebenen muss symmetrisch sein; mit c-Schritten ist h = r₀ und m ~ 10⁶–10⁷, das Kontinuum D·G'' = Γ·G ist genau. Die Rückkehr schließt erst bei Λ, der ganzen inneren Faltung der Zelle (keine Raumlänge; L_H ≈ Λ gilt nur heute) [korr.]; jedes Teilsystem hat L/Λ ≤ 10⁻⁴, die Regel gilt bis 0,66·Λ ≈ Gpc.

## Wo es nicht trägt

- L = r* folgt aus dem Schwellenkriterium (gleichmäßige Kugel), nicht aus der Poisson-Glättung, die ϖ additiv macht.
- Ob D über die Ebenen gleich bleibt (2) oder je Ebene additiv wächst (3b), legt der Ansatz nicht fest. Nur im zweiten Fall gibt es einen Effekt bei 0,68·ℓ und die Folge in 4.
- Γ gleich auf allen Ebenen ist Annahme des Ansatzes; mit Γ = D/r² je Ebene ist R = (1 − μ)/(1 + μ), r₁₀ bleibt 0,67·ℓ.
