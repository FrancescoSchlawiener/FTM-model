# Kosmologie: Leistungsspektrum (CMB, BAO), Wachstum und Torus-Topologie

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Grundannahme a₀(z) = cH(z)/2π (a₀(1100) = 2,45·10⁻⁶, RAR +34 % bei z = 1) → trägt nicht, Abschnitte 1, 2.1–2.3, 4.3, 4.4 und 6.2 gelten nur unter dieser Annahme, im Modell ist a₀ seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); „Störungen grundsätzlich im Tief-a₀-Regime“ → mit fester Schwelle bei der Rekombination an der Schwelle oder darüber (Φ/c² = 10⁻⁵: a_N/a₀ ≈ 2,4 bei 147 Mpc, ≈ 350 bei 1 Mpc); y_H = π „in jeder Epoche“ → nur heute; offene Kernfrage des nicht mitschwingenden Potentials → beantwortet durch die Schließungsregel (∂ξ)² = −1; A_L = 0 bzw. „Wirkt ν auf Licht?“ → die Metrik trägt den wirbelfreien Teil der Regel, für Kugeln η = 1; konstanter Anteil 1 − β „unerklärt“ → Vakuumteil π²/45; „Faktor 2π anpassen?“ → SPARC 1,04–1,19·10⁻¹⁰ im Bereich cH₀/2π; „oberhalb 2·d_LSS unsichtbar“ → oberhalb des 1,5-fachen Durchmessers nicht prüfbar; Tabelle 5.2 C₂-Rest bei L = d_LSS 0,489 (grober Abschnitt) → diskrete Modensumme ≈ 0,84; „Halo-Profile“ → Profile der gelesenen Masse; Umbenennung χ* → d_*, Λ → Λ_kosm, Ω_Λ → Ω_DE, ν(1+L) → ν(1 + d ln ν/d ln y). Die Notiz rechnet mit a₀ = cH₀/2π = 1,04·10⁻¹⁰ (H₀ = 67,36).

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/kosmologie.py`
> Größen: Struktur-Operator ν = 1/(1−q), q = exp(−√y), y = a_N/a₀, und die Kopplung der Schwelle an die Atemrate a₀ = c·H/2π. [korr.] Diese Kopplung trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute.
> Zahlenwert: Diese Notiz rechnet mit dem Modellwert a₀(0) = c·H₀/2π = 1,0416·10⁻¹⁰ m/s² (H₀ = 67,36 km/s/Mpc). Der gemessene Wert ist 1,2·10⁻¹⁰ m/s².

---

## Kurzfassung

1. **Die Schwelle lag früher viel höher.** Mit a₀(z) = c·H(z)/2π war sie bei z = 1100 um den Faktor 2,35·10⁴ größer als heute: a₀(1100) = 2,45·10⁻⁶ m/s². [korr.] Das gilt nur unter der Annahme a₀ ∝ H, die nicht trägt; im Modell ist a₀ seit z ≈ 10 fest, bei der Rekombination 0,76·a₀.
2. **Alle kosmologischen Störungen liegen tief im a₀-Regime.** Es gilt y = π·Ω_x(z)·δ·(R·H/c). Weil δ ≪ 1 ist und R kleiner als der Hubble-Radius, folgt y ≪ 1. Bei der Rekombination ist für realistische Amplituden (δ_b = 10⁻⁵…10⁻³, Φ/c² ~ 10⁻⁵) y = 2·10⁻⁸ bis 10⁻², die Verstärkung ν liegt bei etwa 9 bis 7500. [korr.] Das gilt nur unter a₀ ∝ H. Mit fester Schwelle (a₀(1090) = 0,76·1,2·10⁻¹⁰) liegen die Störungen bei der Rekombination an der Schwelle oder darüber: Φ/c² = 10⁻⁵ gibt bei 147 Mpc (mitbewegt) a_N/a₀ ≈ 2,4 (ν = 1,27), bei 1 Mpc ≈ 350.
3. **Offene Kernfrage: das nicht mitschwingende Potential.** Das CMB verlangt ein Potential, das nicht mit dem Photon-Baryon-Fluid mitschwingt. Ohne ein solches Potential sänke das Verhältnis 3. zu 2. Peak von gemessen 0,98 auf 0,46, und der Winkel des Schallhorizonts verschöbe sich um 32 %. Das Modell soll dieses Potential als Einschränkung der Freiheitsgrade bzw. Gedächtnis des Operators liefern, nicht als zusätzliche Substanz. Wie? [korr.] Beantwortet durch die Schließungsregel (∂ξ)² = −1: Die Druckrichtung der konformen Mode ist geschlossen; die Mode ξ atmet mit (w ≈ 1), wirkt wie druckfreier Staub und ist linear gleich dem Standardmodell bei Q = Ω_c/Ω_b = 5,36.
4. **Wachstum mit Umschaltregel.** Der Operator wirkt nur in Regionen, die sich vom Atmen entkoppelt haben. Daraus folgt eine prüfbare Form von Δf/f(k,z): null für k ≲ 0,05 h/Mpc, anwachsend zu k ≳ 0,1 h/Mpc und zu z → 0. [korr.] Die Zahlen setzen a₀ ∝ H voraus und gelten nur unter dieser Annahme; im Modell leistet die Umschaltung das Gewicht w der Schließungsregel.
5. **Die sauberste Trennbeobachtung** ist die Entwicklung der Radialbeschleunigungs-Relation (RAR) mit z. Im tiefen Regime gilt g_obs/g_bar ∝ √(H(z)/H₀). Das sind +15 % bei z = 0,5, +34 % bei z = 1 und +74 % bei z = 2. Bei fester Rotationsgeschwindigkeit ist die baryonische Masse bei z = 2 um 0,48 dex kleiner. [korr.] Das gilt nur unter a₀ ∝ H; im Modell ist a₀ seit z ≈ 10 fest, die Relation verschiebt sich nicht.
6. **Das Kaleidoskop als große Torus-Zelle.** Mit den Messungen verträglich ist eine Grundzelle mit L ≳ 27,4 Gpc (durch unseren Ort). Im Fenster L ≈ 1,0–1,5·d_LSS ist der Quadrupol um 2–17 % unterdrückt, unter der kosmischen Varianz. Die Signatur sind Korrelationen der a_ℓm zwischen verschiedenen ℓ und m.
7. **Das Atmen als w(z).** Die Lesart ΔL_außen/R = ΣΔs_eff/R mit Δs_eff ∝ √a₀ ergibt ein Maximum der Rückstelldichte. Damit landet w₀ > −1, w_a < 0 im DESI-Quadranten. Die DESI-Form w_a/(1+w₀) ≈ −3,5 verlangt eine Massenskala von ~3·10¹² M☉/h. Das ist noch eine Wahl, keine Vorhersage. [korr.] Der Ansatz setzt Δs_eff ∝ √a₀ ∝ √H voraus (a₀ ∝ H), das trägt nicht; im Modell zählt die innere Faltung je mitbewegter Zelle (1 + Q)·ΣΔs_b fest, Λ ist heute fest, und der konstante Anteil ist der Vakuumteil π²/45, rund ein Drittel der Dunklen Energie.

---

## 0. Eingaben

Planck-2018-artige Werte dienen als Referenz der gemessenen Geometrie (Planck 2018 VI):

| Größe | Wert |
|---|---|
| h, ω_b, ω_c (ΛCDM) | 0,6736; 0,02237; 0,1200 |
| T_CMB, N_eff, n_s, σ₈ | 2,7255 K; 3,046; 0,9649; 0,8111 |
| ω_γ, ω_r | 2,4729·10⁻⁵; 4,1835·10⁻⁵ |
| Ω_b, Ω_c, Ω_r, Ω_DE | 0,0493; 0,2645; 9,220·10⁻⁵; 0,6861 |
| f_b = ω_b/ω_m | 0,1571 |
| c/H₀ | 4450,6 Mpc |
| a₀(0) = c·H₀/2π | 1,0416·10⁻¹⁰ m/s² (0,868 × 1,2·10⁻¹⁰) |

**Kontrolle des Codes.** Das ΛCDM-äquivalente Modell ergibt r_d = 147,10 Mpc, 100·θ* = 1,0396 und z_eq = 3402. Planck-Referenz: 147,09 Mpc; 1,0411; 3387. Die Abweichung in θ* von 0,14 % stammt aus der vereinfachten Neutrino-Behandlung (masselos).

---

## 1. Die Schwelle in der frühen Zeit

Flaches Modell: H²(z) = H₀²·[Ω_r(1+z)⁴ + Ω_b(1+z)³ + Ω_eff(z)]. Gerechnet wird mit dem Hintergrund, den die Messungen verlangen (**Hintergrund A**): nicht mitschwingender Anteil Ω_c = 0,2645 wie in ΛCDM, dazu die Rückstellung Λ_kosm.

[korr.] Die Spalten a₀,A setzen a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest, bei der Rekombination 0,76·a₀; ≈ cH₀/2π gilt nur heute.

| z | H_A [km/s/Mpc] | a₀,A [m/s²] | a₀,A/a₀(0) |
|---|---|---|---|
| 0 | 67,36 | 1,042·10⁻¹⁰ | 1 |
| 0,5 | 89,0 | 1,376·10⁻¹⁰ | 1,321 |
| 1 | 120,5 | 1,863·10⁻¹⁰ | 1,788 |
| 2 | 203,9 | 3,153·10⁻¹⁰ | 3,027 |
| 5 | 557,8 | 8,626·10⁻¹⁰ | 8,281 |
| 10 | 1380 | 2,134·10⁻⁹ | 20,49 |
| 100 | 3,886·10⁴ | 6,009·10⁻⁸ | 577 |
| 1100 | 1,586·10⁶ | 2,452·10⁻⁶ | 2,354·10⁴ |

---

## 2. Zustand der Störungen bei der Rekombination (z = 1100)

### 2.1 Eine allgemeine Formel

[korr.] Diese Rechnung setzt a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute.

Für eine Kugel vom Radius R mit Dichtekontrast δ einer Komponente x gilt g = (4πG/3)·ρ_x·δ·R = ½·Ω_x(z)·H²·R·δ. Geteilt durch a₀ = c·H/2π:

```
y = g / a₀(z) = π · Ω_x(z) · δ · (R·H/c)          (Dichte-Form)
y = 2π · (Φ/c²) · (c/(R·H))                          (Potential-Form, g ≈ Φ/R)
```

**Befund:** Weil a₀ an H hängt, fällt H aus der Dichte-Form bis auf das Verhältnis R/R_Hubble heraus. Eine lineare (δ < 1) Störung innerhalb des Horizonts (R < c/H) hat immer y < π·Ω_x. **Kosmologische Störungen lägen damit grundsätzlich im Tief-a₀-Regime.** [korr.] Das gilt nur unter a₀ ∝ H. Mit fester Schwelle (a₀(1090) = 0,76·1,2·10⁻¹⁰) liegt eine Störung mit Φ/c² = 10⁻⁵ bei 147 Mpc (mitbewegt) bei a_N/a₀ ≈ 2,4 (ν = 1,27), bei 1 Mpc bei ≈ 350 (a_N ≈ Φ/λ_phys): an der Schwelle oder darüber, nicht tief darunter. Unter a₀ ∝ H erreicht eine Störung die Schwelle y = 1 erst bei: [korr.]

| R (mitbewegt) | δ₁ (nur Baryonen) | δ₁ (gesamte Materie, Hintergrund A) |
|---|---|---|
| 147 Mpc | 3,80 | 0,60 |
| 10 Mpc | 55,8 | 8,77 |
| 1 Mpc | 558 | 87,7 |

### 2.2 Zahlen bei z = 1100

[korr.] Die Spalten y und ν setzen a₀(1100) = cH(1100)/2π voraus; das trägt nicht (der Vakuumteil bliebe zu jeder Zeit bei 22 %, wie ΔN_eff = 2,09). Im Modell ist a₀ bei der Rekombination 0,76·1,2·10⁻¹⁰; damit läge in Tabelle b die Zeile 147 Mpc bei a_N/a₀ ≈ 2,4 (ν = 1,27) statt 8,90·10⁻⁵; eine Störung mit Φ/c² = 10⁻⁵ bei 1 Mpc läge bei ≈ 350.

Werte bei z = 1100: H = 5,139·10⁻¹⁴ s⁻¹, a₀ = 2,452·10⁻⁶ m/s², Ω_b(1100) = 0,1187, Ω_m(1100) = 0,7556 (Hintergrund A), ρ_b = 5,61·10⁻¹⁹ kg/m³.

**a) Quelle: nur Baryonen** (Auszug)

| R mitb. | δ_b | g [m/s²] | y | ν(y) | lineare Antwort ν(1 + d ln ν/d ln y) |
|---|---|---|---|---|---|
| 147 Mpc | 10⁻⁵ | 6,46·10⁻¹² | 2,63·10⁻⁶ | 617 | 309 |
| 147 Mpc | 10⁻³ | 6,46·10⁻¹⁰ | 2,63·10⁻⁴ | 62,1 | 31,3 |
| 10 Mpc | 10⁻⁴ | 4,39·10⁻¹² | 1,79·10⁻⁶ | 748 | 374 |
| 1 Mpc | 10⁻⁵ | 4,39·10⁻¹⁴ | 1,79·10⁻⁸ | 7471 | 3736 |

**b) Gesamtes Potential** (Sachs-Wolfe-Niveau Φ/c² ~ 10⁻⁵)

| R mitb. | Φ/c² | g [m/s²] | y | ν(y) |
|---|---|---|---|---|
| 147 Mpc | 10⁻⁵ | 2,18·10⁻¹⁰ | 8,90·10⁻⁵ | 106,5 |
| 50 Mpc | 10⁻⁵ | 6,41·10⁻¹⁰ | 2,62·10⁻⁴ | 62,3 |
| 10 Mpc | 10⁻⁵ | 3,21·10⁻⁹ | 1,31·10⁻³ | 28,2 |
| 1 Mpc | 10⁻⁴ | 3,21·10⁻⁷ | 0,131 | 3,30 |

Die Tabellen sind Auszüge (ausgewählte R und δ).

### 2.3 Folgen für ν

[korr.] Diese Folgerungen setzen a₀(z) = cH(z)/2π voraus; das trägt nicht (der Vakuumteil bliebe zu jeder Zeit bei 22 %, wie ΔN_eff = 2,09). Im Modell ist a₀ seit z ≈ 10 fest, bei der Rekombination 0,76·a₀.

- Im Tief-Regime gilt ν ≈ 1/√y. Die Tabelle bestätigt das: ν = 617 gegen 1/√y = 616. Die Verstärkung wächst also mit **kleinerem** δ. Der Operator ist **nichtlinear in der Quelle**, denn g ∝ √(g_N·a₀) ∝ √δ.
- Für eine kleine Zusatzstörung auf einem Hintergrundfeld wirkt längs des Feldes nur die linearisierte Antwort ν·(1 + d ln ν/d ln y). Das ist ν/2.
- Zum Vergleich die Hubble-Kugel als Ganzes. Mit g_H = c·H/2 ist y_H = π heute (a₀ ≈ cH₀/2π), nicht in jeder Epoche [korr.], und daraus folgt ν(π) = 1,205. Der globale Faltungszustand selbst liegt knapp über der Schwelle, jede Einzelstörung liegt tief darunter.

**Befund 2.** [korr.] Mit fester Schwelle (a₀(1090) = 0,76·1,2·10⁻¹⁰) liegen die Störungen bei z = 1100 an der Schwelle oder darüber, nicht tief darunter: Φ/c² = 10⁻⁵ gibt bei 147 Mpc (mitbewegt) a_N/a₀ ≈ 2,4 (ν = 1,27), bei 1 Mpc ≈ 350 (a_N ≈ Φ/λ_phys). Die Faktoren 10 bis 7000 gelten nur unter a₀ ∝ H. Linear wirkt die Schließungsregel mitatmend (w ≈ 1) wie das Standardmodell bei Q; in die Form des Operators geht die Mode ξ erst in gebundenen Bereichen über (Abschnitt 4.3). Was die Peaks verlangen, zeigt Abschnitt 3.

---

## 3. Was die akustischen Peaks verlangen

### 3.1 Das nicht mitschwingende Potential (halbquantitativ)

Das Photon-Baryon-Fluid schwingt im Potential Ψ. Im einfachsten Oszillatorbild mit konstantem Ψ (Hu & Sugiyama 1995) gilt:

```
(Θ₀ + Ψ)(k) = (1/3 + R)·Ψ·cos(k r_s) − R·Ψ ,     R = 3ρ_b/(4ρ_γ)
```

- Mit R(z*) = 0,622 ist die Amplitude der ungeraden (Kompressions-) zu den geraden Peaks gleich 1 + 6R = 4,73. Das gilt **nur, wenn Ψ nicht mitschwingt**. Die Nulllage −RΨ braucht ein festes Potential.
- Wird Ψ vom schwingenden Fluid selbst getragen, zerfällt es beim Eintritt in den Horizont (Strahlungsantrieb). Die Nulllage verschwindet, die Ungerade/Gerade-Modulation bricht ein, und die höheren Peaks werden gedämpft.
- In ΛCDM hält CDM (Druck null, keine Kopplung an die Photonen) Ψ nach der Gleichheitsepoche fest.

Fitformeln für die Peak-Höhen (Hu, Fukugita, Zaldarriaga & Tegmark 2001):

| Fall | ω_m | H₂ = 2./1. | H₃ = 3./1. | H₃/H₂ |
|---|---|---|---|---|
| ΛCDM (Planck-Werte) | 0,1424 | 0,460 | 0,450 | 0,978 |
| nicht mitschwingender Anteil halb so stark | 0,0824 | 0,439 | 0,346 | 0,788 |
| ohne nicht mitschwingenden Anteil, nur Baryonen (extrapoliert) | 0,0224 | 0,372 | 0,171 | 0,460 |

Ohne nicht mitschwingendes Potential wäre der 3. Peak nur halb so hoch wie der 2. Gemessen sind beide etwa gleich hoch. Das entspricht der ΛCDM-Zeile, die an Planck angepasst ist. Für ω_m = 0,0224 liegt die Fitformel außerhalb ihres Eichbereichs, die Zeile zeigt also nur die Tendenz.

Die Geometrie der Peaks verlangt denselben Anteil im Hintergrund:

| Hintergrund | z_eq | r_s(z*) | r_d | D_M(z*) | 100·θ* |
|---|---|---|---|---|---|
| A: gemessener Hintergrund (ω_m = 0,1424) + Λ_kosm | 3402 | 144,44 Mpc | 147,10 Mpc | 13,894 Gpc | 1,0396 |
| B: nur Baryonen + Λ_kosm | 534 | 195,12 Mpc | 199,55 Mpc | 27,745 Gpc | 0,7033 |

### 3.2 Die offene Kernfrage

Wirkungen sind im Modell Geometrie bzw. Struktur, keine Substanz (Beobachter-Regel). Das CMB verlangt ein Potential, das nicht mit dem Photon-Baryon-Fluid mitschwingt: keine Kopplung an den Photonendruck, Störungen, die dem Potential folgen und nicht schwingen, und im Hintergrund vor z ~ 1000 der gemessene Anteil ω_m − ω_b = 0,1200. Das ist das 5,36-Fache der Baryonen (ω_m/ω_b = 6,364).

**Frage:** Wie liefert das Modell dieses Potential als Einschränkung der Freiheitsgrade bzw. Gedächtnis des Operators, nicht als zusätzliche Substanz? In gebundenen, vom Atmen entkoppelten Systemen übernimmt der Operator ν diese Rolle (Abschnitt 4.3).

[korr.] Beantwortet durch die Schließungsregel (∂ξ)² = −1: eine eingeschränkte konforme Mode der Metrik, deren Druckrichtung geschlossen ist (keine Kopplung an den Photonendruck, kein eigener Druck). Mitatmend (w ≈ 1) wirkt sie wie druckfreier Staub; linear ist sie damit gleich dem Standardmodell bei Q = Ω_c/Ω_b = 5,36 (aus der Messung eingesetzt). In gebundenen Bereichen (w → 0) geht die Mode mit λ_m = ω_orb/N* in die Form des Operators über. Ein Beobachter, der ungeänderte Regeln annimmt, liest ihre Wirkung als Dichte Q·ρ_b; ein Stoff ist das nicht.

### 3.3 Die nötige Rechnung: Boltzmann-Code mit FTM-Sektor

Grundlage sind CLASS (Blas, Lesgourgues & Tram 2011) oder CAMB (Lewis, Challinor & Lasenby 2000). Für den späten Operator eignen sich MGCAMB (Zhao et al. 2009) oder hi_class (Zumalacárregui et al. 2017). Konforme Newton-Eichung (Ma & Bertschinger 1995):

```
ds² = a²(τ_c)·[ −(1+2ψ) dτ_c² + (1−2φ) dx² ] ,   ℋ = a'/a,   ' = d/dτ_c (konforme Zeit τ_c)
```

**Hintergrund**

```
ℋ² = (8πG/3)·a²·(ρ_γ + ρ_ν + ρ_b + ρ_R) + Beitrag des nicht mitschwingenden Potentials (Mode ξ der Schließungsregel, wie druckfreier Staub mit Q·ρ_b gelesen, 3.2) [korr.]
ρ_R' = −3ℋ(1 + w_R(a))ρ_R              (Rückstellung, Abschnitt 6)
```

**Einstein-Gleichungen (Störungen)**

```
k²φ + 3ℋ(φ' + ℋψ) = −4πG a² Σ_i ρ_i δ_i
k²(φ' + ℋψ)        =  4πG a² Σ_i (ρ_i + P_i) θ_i
k²(φ − ψ)          = 12πG a² Σ_i (ρ_i + P_i) σ_i
```

**Baryonen** (unverändert, Kopplung nur an die Photonen)

```
δ_b' = −θ_b + 3φ'
θ_b' = −ℋθ_b + c_b² k² δ_b + (4ρ_γ/3ρ_b)·a n_e σ_T (θ_γ − θ_b) + k²ψ
```

**Später Operator** (quasistatisch, ab der Entkopplung vom Atmen, Abschnitt 4)

```
k²ψ       = −4πG a² μ(k,a) ρ̄_m Δ_m ,        μ = 1 + A·ε(k,a)
k²(φ + ψ) = −8πG a² Σ(k,a) ρ̄_m Δ_m ,        Σ = 1 + A_L·ε(k,a)
```

Wirkt der Operator nur auf Bahnen und nicht auf Lichtablenkung (Trennklasse I), dann ist A_L = 0. Das ergibt einen messbaren Unterschied zwischen dynamischem und Linsen-Potential. [korr.] Im Modell ist A_L ≠ 0: Die Metrik trägt den wirbelfreien Teil der Regel (Ψ = Φ); für Kugeln und außerhalb von Scheiben ist η = 1 (A_L = A), nur in Scheiben innen sieht Licht weniger.

**Nichtlinear** (N-Körper, Poisson-Löser mit epochenabhängiger Schwelle; mitbewegte Koordinaten x, ∇_r = ∇_x/a)

```
∇_x²Φ = ∇_x·[ ν(|∇_xΦ_N| / (a·a₀(z))) ∇_xΦ_N ] ,   a₀(z) = c H(z)/2π ,   ∇_x²Φ_N = 4πG a² (ρ_b − ρ̄_b)
```

[korr.] Die epochenabhängige Schwelle a₀(z) = cH(z)/2π trägt nicht (der Vakuumteil bliebe zu jeder Zeit bei 22 %, wie ΔN_eff = 2,09); im Modell ist a₀ seit z ≈ 10 fest.

**Prüfschritte.**
1. Mit dem gemessenen Hintergrund (ΛCDM-Anteil ω_m) und A = 0 muss exakt ΛCDM herauskommen (Code-Kontrolle).
2. Die Grenze an A und ε aus TT/TE/EE, Linsenpotential, BAO und fσ₈ bestimmen.

---

## 4. BAO und Wachstum der Strukturen

### 4.1 Schallhorizont und BAO-Skala

Die BAO-Skala ist r_d = ∫ c_s dz/H von z_drag = 1059,94 bis ∞. Im ΛCDM-äquivalenten Hintergrund ist r_d = 147,10 Mpc. Jede Veränderung des nicht mitschwingenden Anteils vor z ~ 1000 verschiebt r_d und damit alle BAO-Abstände. Zum Vergleich der Hintergrund ohne diesen Anteil:

| z_eff | ΛCDM D_V/r_d | nur Baryonen |
|---|---|---|
| 0,295 | 8,060 | −19,83 % |
| 0,510 | 12,842 | −14,54 % |
| 0,706 | 16,476 | −9,57 % |
| 0,934 | 19,968 | −3,87 % |
| 1,321 | 24,499 | +5,06 % |
| 1,484 | 26,005 | +8,45 % |
| 2,330 | 31,398 | +22,24 % |

DESI misst diese Größen auf etwa ein Prozent genau. Der Hintergrund muss deshalb vor z ~ 1000 den gemessenen Anteil tragen; dann ist die BAO-Skala die von ΛCDM.

### 4.2 Wachstum: lineare Moden

Mit dem linearen Spektrum (Eisenstein & Hu 1998, ohne Wiggles, σ₈ = 0,8111) ist y(k,z) = π·Ω_m(z)·Δ(k,z)·aH/(ck):

| k [h/Mpc] | z | Δ(k,z) | y | ν(y) |
|---|---|---|---|---|
| 0,01 | 0 | 0,033 | 1,08·10⁻³ | 30,9 |
| 0,05 | 0 | 0,285 | 1,87·10⁻³ | 23,6 |
| 0,10 | 0 | 0,540 | 1,78·10⁻³ | 24,2 |
| 0,20 | 1 | 0,532 | 1,96·10⁻³ | 23,1 |
| 0,50 | 0 | 1,392 | 9,16·10⁻⁴ | 33,6 |

Lineare Moden liegen also tief im a₀-Regime. Die gemessenen Wachstumsraten verlangen auf linearen Skalen D ∝ a (Referenzwerte ΛCDM: fσ₈ = 0,473 (z = 0,3), 0,474 (0,51), 0,461 (0,71), 0,439 (0,93), 0,395 (1,32)). Der Operator wirkt deshalb nicht auf den mitatmenden linearen Fluss, sondern nur in entkoppelten Regionen (4.3).

### 4.3 Die Umschaltregel des Modells und die prüfbare Formel

[korr.] Diese Rechnung setzt die Schwelle a₀ = cH/2π als Beschleunigung des atmenden Mediums voraus (y_ta ∝ H(z)); das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest; ≈ cH₀/2π gilt nur heute. Die Umschaltung leistet im Modell das Gewicht w der Schließungsregel: Vor dem Umkehrpunkt fällt w innerhalb von 7 Mio. Jahren von 0,9 auf 0,1, und die Mode ξ geht mit λ_m = ω_orb/N* in die Form des Operators über.

Das FTM liefert eine natürliche Regel: **Der Operator wirkt dort, wo eine Region nicht mehr mitatmet.** Die Schwelle a₀ = c·H/2π ist selbst die Beschleunigung des atmenden Mediums. Eine Region, die noch mit dem Hubble-Fluss expandiert, ist kein Trichter. Die Entkopplung beginnt mit der Umkehr (Turnaround, linearer Kontrast δ_ta = 1,062, nichtlinear 1 + Δ_ta = 9π²/16 = 5,55). Dort übernimmt der Operator die Rolle des nicht mitschwingenden Potentials:

```
μ(k,z) = 1 + A · ε(k,z)
ε(k,z) = F_dec(k,z) · [ f_b · ν(y_ta(k,z)) − 1 ]
F_dec  = erfc( δ_ta / (√2 · σ(R = 1/k, z)) )                         (Anteil entkoppelter Masse)
y_ta   = π · Ω_m(z) · f_b · (1+Δ_ta) · r_ta · H(z)/c ,   r_ta = a·R/(1+Δ_ta)^{1/3}
```

A ist der Umschaltgrad (0 = ΛCDM, 1 = das nicht mitschwingende Potential wird in entkoppelten Regionen vollständig von ν getragen). Die Wachstumsgleichung

```
D'' + (2 + d ln H/d ln a) D' − (3/2) Ω_m(a) μ(k,a) D = 0      (' = d/d ln a),   f = d ln D/d ln a
```

ergibt:

| k [h/Mpc] | z | F_dec | ν(y_ta) | f_b·ν | ε (A = 1) |
|---|---|---|---|---|---|
| 0,05 | 0 | 0,0071 | 18,1 | 2,84 | 0,013 |
| 0,10 | 0 | 0,125 | 25,4 | 3,98 | 0,373 |
| 0,10 | 1 | 0,0116 | 17,1 | 2,69 | 0,020 |
| 0,20 | 0 | 0,333 | 35,6 | 5,60 | 1,530 |
| 0,20 | 1 | 0,111 | 24,0 | 3,77 | 0,306 |

| A | k [h/Mpc] | Δf/f z = 0 | z = 0,5 | z = 1 | z = 1,5 | ΔD/D (z = 0) |
|---|---|---|---|---|---|---|
| 1 | 0,05 | 0,16 % | 0,005 % | 0 | 0 | 0,013 % |
| 1 | 0,10 | 7,72 % | 1,63 % | 0,26 % | 0,03 % | 1,27 % |
| 1 | 0,20 | 40,99 % | 17,03 % | 6,49 % | 2,22 % | 11,60 % |
| 0,1 | 0,10 | 0,78 % | 0,16 % | 0,03 % | 0,003 % | 0,13 % |
| 0,1 | 0,20 | 4,26 % | 1,73 % | 0,65 % | 0,22 % | 1,13 % |

**Die prüfbare Aussage für DESI und Euclid** (nur unter a₀ ∝ H und dieser Umschaltregel [korr.]):
1. **Skalenabhängigkeit.** Δf/f ≈ 0 für k ≲ 0,05 h/Mpc. Der Effekt setzt bei k ≈ 0,1 h/Mpc ein und wächst steil mit k. ΛCDM erwartet dort ∂f/∂k ≈ 0.
2. **Epochenabhängigkeit.** Δf/f wächst zu z → 0, bei k = 0,2 h/Mpc um den Faktor 6,3–6,5 zwischen z = 1 und z = 0.
3. **Form fest, Amplitude A frei.** Für kleines A ist Δf/f ∝ A (Tabelle: A = 0,1 gibt ein Zehntel von A = 1). Das Verhältnis Δf/f(0,2)/Δf/f(0,1) bei z = 0 ist für kleines A unabhängig von A (5,49 bei A = 0,1; 5,5 bei A = 0,02; 5,3 bei A = 1). Gemessen wird also die **Form**.
4. **Grenze.** A = 1 gäbe 17 % bei z = 0,5 und k = 0,2 h/Mpc. Das steht mit den heutigen fσ₈-Messungen in Spannung und muss mit dem vollen Code (3.3) geprüft werden. Realistisch ist A ≲ 0,1.
5. Ein Vorbehalt: Die Formel nutzt die lineare Wachstumsgleichung auf quasilinearen Skalen und die Press-Schechter-Statistik (Press & Schechter 1974). Sie ist halbquantitativ und gibt die Form, nicht die letzte Stelle.

### 4.4 Die schärfste Trennbeobachtung: die Schwelle wandert mit H

[korr.] Diese Rechnung setzt a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest; ≈ cH₀/2π gilt nur heute. Die Schwelle wandert nicht.

Unabhängig von der Umschaltregel folgt aus a₀(z) = c·H(z)/2π für jede gebundene Galaxie g_obs = g_bar·ν(g_bar/a₀(z)). Im tiefen Regime gilt g_obs/g_bar ∝ √(a₀(z)/a₀(0)):

| z | a₀(z)/a₀(0) | g_obs/g_bar (tief, relativ) | Δlog v bei festem M | Δlog M_bar bei festem v | ν(g = 10⁻¹¹) relativ |
|---|---|---|---|---|---|
| 0,25 | 1,140 | 1,068 | +0,014 dex | −0,057 dex | 1,058 |
| 0,5 | 1,321 | 1,149 | +0,030 dex | −0,121 dex | 1,128 |
| 1,0 | 1,788 | 1,337 | +0,063 dex | −0,252 dex | 1,288 |
| 1,5 | 2,364 | 1,538 | +0,093 dex | −0,374 dex | 1,460 |
| 2,0 | 3,026 | 1,740 | +0,120 dex | −0,481 dex | 1,633 |
| 3,0 | 4,557 | 2,135 | +0,165 dex | −0,659 dex | 1,972 |

Messbar ist das mit Galaxie-Galaxie-Linsen als RAR in Rotverschiebungsbins (KiDS-1000: Brouwer et al. 2021, heute bei kleinem z; Euclid bis z ~ 1) und mit Kinematik bei z ≈ 1–2 (z. B. Genzel et al. 2017). [korr.] Unter der Kopplung a₀ ∝ H ergäbe sich eine Verschiebung der RAR um den Faktor 1,34 bei z = 1. Das Modell sagt dagegen eine zeitlich feste Relation voraus: a₀ ist seit z ≈ 10 fest, eine Verdopplung bis z = 1 ist ausgeschlossen.

---

## 5. Torus-Topologie und Kaleidoskop, quantitativ

### 5.1 Grundgrößen

Der Abstand zur letzten Streufläche ist d_* = 13,894 Gpc, ihr Durchmesser d_LSS = 27,788 Gpc. Galaxien-Durchmusterungen reichen bis z = 3, das sind 6,516 Gpc.

**Messgrenzen:**
- Planck 2015 XVIII, kubischer Torus: Der Inkugelradius erfüllt R_i > 0,97·d_*, also L > 26,95 Gpc.
- COMPACT (Petersen et al. 2023): Die kürzeste nicht zusammenziehbare Schleife **durch unseren Ort** ist länger als 0,985·d_LSS, also L > 27,37 Gpc. Anderswo kann die kürzeste Schleife um den Faktor 2 bis ≥ 6 kürzer sein, also 13,69 bis 4,56 Gpc.

### 5.2 Kubischer 3-Torus mit Kantenlänge L

ℓ_min = π·d_LSS/L. Kreispaare („circles in the sky“, Cornish, Spergel & Starkman 1998) gibt es für L < d_LSS, mit Winkelradius cos α = L/d_LSS. Der C_ℓ-Rest ist ein grober Sachs-Wolfe-Abschnitt bei k_min = 2π/L, ohne integrierten Sachs-Wolfe-Effekt und ohne die Diskretheit der Moden.

| L/d_LSS | L [Gpc] | ℓ_min | α | Eigenbilder < d_* | Eigenbilder < z = 3 | C₂ Rest | C₃ Rest | C₄ Rest |
|---|---|---|---|---|---|---|---|---|
| 0,25 | 6,95 | 12,6 | 75,5° | 26 | 0 | 0,020 | 0,036 | 0,071 |
| 0,50 | 13,89 | 6,28 | 60,0° | 0 | 0 | 0,088 | 0,134 | 0,311 |
| 0,75 | 20,84 | 4,19 | 41,4° | 0 | 0 | 0,189 | 0,555 | 0,856 |
| 0,97 | 26,95 | 3,24 | 14,1° | 0 | 0 | 0,455 | 0,827 | 0,969 |
| 0,985 | 27,37 | 3,19 | 9,9° | 0 | 0 | 0,472 | 0,838 | 0,972 |
| 1,0 | 27,79 | 3,14 | – | 0 | 0 | 0,489 (diskret ≈ 0,84) [korr.] | 0,848 | 0,974 |
| 1,2 | 33,35 | 2,62 | – | 0 | 0 | 0,673 | 0,934 | 0,993 |
| 1,5 | 41,68 | 2,09 | – | 0 | 0 | 0,831 | 0,979 | 0,998 |
| 2,0 | 55,58 | 1,57 | – | 0 | 0 | 0,936 | 0,996 | 1,000 |
| 3,0 | 83,36 | 1,05 | – | 0 | 0 | 0,986 | 1,000 | 1,000 |

[korr.] Die Spalten C₂–C₄ sind der grobe Abschnitt bei k_min; die diskrete Modensumme gibt bei L = d_LSS einen C₂-Rest von ≈ 0,84 statt 0,489. Im Fenster 1,0–1,5·d_LSS ist der Quadrupol damit um 2–17 % gesenkt, weniger als seine kosmische Streuung (63 %).

Die kosmische Varianz des Quadrupols beträgt σ(C₂)/C₂ = 0,632. Der gemessene Quadrupol ist niedrig, und über 60° fehlen Korrelationen (Copi, Huterer, Schwarz & Starkman 2015). Beides ist auffällig, aber einzeln nicht zwingend.

### 5.3 Das Kaleidoskop als große Torus-Zelle

Das Kaleidoskop hat im Modell drei Teile:

- **K1 (topologisch).** Die Grundzelle ist durch unseren Ort mindestens 27,4 Gpc groß (≥ 0,985·d_LSS). Innerhalb des Galaxien-Volumens gibt es **kein** Eigenbild. Das Kaleidoskop des Lichts reicht höchstens bis an die letzte Streufläche.
- **K2 (Verdrehung).** Die „gegeneinander verdrehten Pfade“ des Autors sind in kompakten flachen Räumen nur mit Drehwinkeln von 60°, 90°, 120° und 180° möglich. Eine irrationale Verdrehung (Wicklungszahl W des Torus) geht nur in einer Richtung, in der der Raum nicht geschlossen ist (Schrauben-Identifikation). In verdrehten Topologien darf die kürzeste Schleife anderswo 4,6–13,7 Gpc sein (COMPACT). Dort bleibt Raum für eine lokale Mehrfachabbildung, die wir von unserem Ort aus nicht direkt sehen.
- **K3 (fraktal statt Kopie).** Was man als „Wiederholung“ sieht, ist im FTM die **Selbstähnlichkeit** der Ebenen (g(s²x, sy) = g(x, y)) und der Blick in die Zeit. Es sind keine identischen Kopien derselben Galaxien. Die Zahl der realen Grundsysteme ist dann die Zahl der beobachteten Galaxien.

**Prüfbare Vorhersage (FTM-T).** Soll der geschlossene Torus des Modells messbar sein, muss L im Fenster 0,985–1,5·d_LSS (27–42 Gpc) liegen. Dort gibt es keine Kreispaare (oder nur sehr kleine, α < 10°), aber:
- einen um 2–17 % unterdrückten Quadrupol (C₂-Rest 0,83–0,98 nach der diskreten Modensumme) bei fast unverändertem ℓ = 4. Das liegt unter der kosmischen Varianz; allein über das Leistungsspektrum ist die Zelle nicht nachweisbar;
- Korrelationen zwischen a_ℓm mit verschiedenem ℓ, m, wie sie COMPACT IIa (Eskilt et al. 2024) für die Nachweisbarkeit berechnet. Das ist der eigentliche Test;
- eine Vorzugsrichtung der Grundzelle, die in Temperatur **und** Polarisation gleich sein muss.

Oberhalb des 1,5-fachen Durchmessers (L > 1,5·d_LSS) ist der Torus nicht mehr nachweisbar; das Modell wäre dann in dieser Hinsicht nicht prüfbar. [korr.]

---

## 6. Atmen und Dunkle Energie: w(z)

### 6.1 Was DESI zeigt

DESI DR2 (arXiv:2503.14738) mit der CPL-Form w(a) = w₀ + w_a(1−a) (Chevallier & Polarski 2001; Linder 2003):

| Kombination | w₀ | w_a | w = −1 bei z_c | w_a/(1+w₀) | w(z = 1) | w(z = 2) |
|---|---|---|---|---|---|---|
| DESI+CMB+Pantheon+ | −0,838 | −0,62 | 0,354 | −3,83 | −1,148 | −1,251 |
| DESI+CMB+Union3 | −0,667 | −1,09 | 0,440 | −3,27 | −1,212 | −1,394 |
| DESI+CMB+DESY5 | −0,752 | −0,86 | 0,405 | −3,47 | −1,182 | −1,325 |

**Übersetzung:** Die Dichte der Rückstellung hatte bei z_c ≈ 0,35–0,44 ein **Maximum** und nimmt seither ab.

### 6.2 Die Form, die das Modell nahelegt

[korr.] Dieser Rückstell-Ansatz setzt Δs_eff ∝ √(GMa₀) ∝ √H voraus, also a₀ ∝ H; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell zählt die innere Faltung je mitbewegter Zelle (1 + Q)·ΣΔs_b fest, Λ ist heute fest, a₀ wandert nicht. Der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie; die übrigen zwei Drittel werden in der Rückstellung als Bindungsrate dF/d ln a gesucht.

Das Modell liefert zwei Zutaten. Erstens die Bilanz ΔL_außen/R = ΣΔs_innen/R. Zweitens gilt für jede kollektive Struktur Δs_eff = Δs·ν, im tiefen Regime also Δs_eff ∝ √(G M a₀) ∝ √H. Mit dem Anteil F_coll(>M, a) der Masse in gebundenen Strukturen und der physikalischen Länge R ∝ a ergibt sich der **Rückstell-Ansatz**:

```
ρ_R(a) = ρ_R0 · [ (1−β) + β · S(a)/S(1) ] ,     S(a) = F_coll(>M, a) · √(H(a)/H₀) · a^(−m)
w(a)   = −1 − (1/3) · [β S / ((1−β) S(1) + β S)] · [ d ln F_coll/d ln a + ½ d ln H/d ln a − m ]
```

Mit physikalischer Länge R ist m = 1.

**Die Physik der Form.** Drei Tendenzen treten gegeneinander an: die Faltung wächst mit der Strukturbildung (F_coll ↑), die Schwelle sinkt mit dem Atmen (√H ↓), und der Außenraum dehnt sich (1/a ↓). Solange die Faltung schneller wächst als die beiden anderen Effekte, steigt die Rückstelldichte (w < −1). Danach fällt sie (w > −1). Der Wechsel liegt bei d ln F_coll/d ln a = m − ½·d ln H/d ln a.

| M [M☉/h] | m | z_max von S | w₀ (β = 0,3) | w_a (β = 0,3) | w_a/(1+w₀) |
|---|---|---|---|---|---|
| 10¹⁰ | 1 | 2,16 | −0,881 | −0,022 | −0,19 |
| 10¹¹ | 1 | 1,43 | −0,874 | −0,161 | −1,28 |
| 10¹² | 1 | 0,78 | −0,875 | −0,337 | −2,70 |
| 3·10¹² | 1 | 0,50 | −0,889 | −0,407 | −3,66 |
| 10¹³ | 1 | 0,22 | −0,929 | −0,414 | −5,80 |
| 10¹⁴ | 1 | 0,00 | −1,127 | +0,018 | −0,14 |

Mit β so gewählt, dass w₀ = −0,752 (DESY5):

| M [M☉/h] (m = 1) | β | w_a | DESY5 |
|---|---|---|---|
| 10¹¹ | 0,618 | −0,333 | −0,86 (+0,23/−0,20) |
| 10¹² | 0,618 | −0,669 | |
| 3·10¹² | 0,652 | −0,913 | |
| 10¹³ | 0,760 | −1,360 | |

**Befund 6** (nur unter der Annahme a₀ ∝ H [korr.]).
1. Die Lesart mit physikalischer Länge (m = 1) legt unter dieser Annahme ein **Maximum der Rückstellung** und damit den DESI-Quadranten w₀ > −1, w_a < 0 nahe. [korr.]
2. Die Form w_a/(1+w₀) ≈ −3,3 bis −3,8 trifft unter dieser Annahme die Massenskala 3·10¹² M☉/h (Gruppen), mit w_a = −0,91 gegen −0,86 gemessen. [korr.] Diese Skala ist gewählt und nicht hergeleitet. Setzt man sie bei ~10¹⁰ M☉ an, ab der Galaxien kollektiv sind, liegt das Verhältnis bei −0,2 bis −1,3.
3. Ein konstanter Anteil 1 − β ≈ 0,35 bleibt. [korr.] Im Modell ist er der Vakuumteil Ω_s = π²/45, rund ein Drittel der Dunklen Energie.

---

## 7. Offene Punkte als nächste Rechenschritte

1. **Boltzmann-Lauf (CLASS/hi_class).** Die Gleichungen aus 3.3 einbauen und zwei Kontrollen rechnen:
   - (a) gemessener Hintergrund, A = 0 muss ΛCDM exakt reproduzieren;
   - (b) A und A_L frei, Grenzen aus Planck + DESI + fσ₈.
2. **Das nicht mitschwingende Potential (Kernfrage).** Das CMB verlangt ein Potential, das nicht mit dem Photon-Baryon-Fluid mitschwingt, in der Stärke des 5,36-Fachen der Baryonen. Wie liefert das Modell es als Einschränkung der Freiheitsgrade bzw. Gedächtnis des Operators, nicht als zusätzliche Substanz? [korr.] Beantwortet durch die Schließungsregel (∂ξ)² = −1: Die Mode ξ atmet mit (w ≈ 1), wirkt wie druckfreier Staub und ist linear gleich dem Standardmodell bei Q = 5,36.
3. **Übergang in entkoppelten Regionen in N-Körper-Simulationen.** Den Poisson-Löser aus 3.3 mit a₀(z) ∝ H und der Umschaltregel „Umkehr = Entkopplung“ bauen. [korr.] a₀ ∝ H trägt nicht (Vakuumteil wie ΔN_eff = 2,09); im Modell ist a₀ seit z ≈ 10 fest, die Umschaltung leistet das Gewicht w der Schließungsregel. Zu messen sind die Profile der gelesenen Masse [korr.], die RAR bei z = 0–2 und fσ₈(k). Ziel ist es, das heuristische ε(k,z) aus 4.3 durch Simulationsergebnisse zu ersetzen.
4. **RAR-Evolution gegen Daten.** Die Tabelle 4.4 ist mit Linsen-RAR in z-Bins (KiDS, Euclid) und mit Kinematik bei z ≈ 1–2 zu vergleichen. Das ist der schnellste Test von a₀ ∝ H. [korr.] Das Modell sagt eine zeitlich feste Relation voraus (a₀ seit z ≈ 10 fest).
5. **Linsen gegen Dynamik.** A_L gegen A festlegen. Wirkt ν auch auf Licht? [korr.] Ja: Die Metrik trägt den wirbelfreien Teil der Regel (Ψ = Φ); für Kugeln und außerhalb von Scheiben ist η = 1, in Scheiben innen η > 1. Die Größe E_G bzw. der Quotient der Potentiale φ/ψ ist eine direkte Prüfung von Trennklasse I.
6. **Topologie-Likelihood.** Für L = 0,985–1,5·d_LSS die volle Kovarianz der a_ℓm (Eigenmoden des 3-Torus, wie in COMPACT IIa) rechnen, mit ISW und diskreten Moden. Die grobe C₂-Rest-Tabelle aus 5.2 ist dadurch zu ersetzen.
7. **Rückstell-Ansatz festlegen.** Die Massenskala M aus dem Modell herleiten statt sie zu wählen. Kandidaten sind die Skala, ab der Strukturen im Mittel die Schwelle unterschreiten, oder die Gewichtung aller kollektiven Strukturen mit ihrem ν. Danach die volle Distanz-Likelihood (DESI + SN) statt des CPL-Fits auf w(a) rechnen. Der konstante Rest 1 − β ist im Modell der Vakuumteil Ω_s = π²/45, rund ein Drittel der Dunklen Energie. [korr.]
8. **Schwelle a₀(0).** c·H₀/2π = 1,04·10⁻¹⁰ liegt 13 % unter 1,2·10⁻¹⁰. Zu klären ist, ob der Faktor 2π oder die Definition von H im Operator (lokal gegenüber global) anzupassen ist. [korr.] Nicht nötig: An SPARC liegt der gelesene Wert mit realistischen Fehlern von Abstand, Neigung und Masse-zu-Licht bei 1,04–1,19·10⁻¹⁰ (Υ_d = 0,5), im Bereich von cH₀/2π. Das ändert alle Tabellen in 4.4 nur in der Normierung, nicht im Verlauf.

---

## Literatur

- Planck Collaboration 2020, *Planck 2018 results VI. Cosmological parameters*, A&A 641, A6.
- Planck Collaboration 2016, *Planck 2015 results XVIII. Background geometry and topology of the Universe*, A&A 594, A18.
- COMPACT: Petersen P. et al. 2023, *Cosmic topology. Part I. Limits on orientable Euclidean manifolds from circle searches*, JCAP 01, 030.
- COMPACT: Eskilt J. R. et al. 2024, *Cosmic topology. Part IIa. Eigenmodes, correlation matrices, and detectability of orientable Euclidean manifolds*, JCAP 03, 036.
- Cornish N. J., Spergel D. N., Starkman G. D. 1998, *Circles in the sky: finding topology with the microwave background radiation*, Class. Quantum Grav. 15, 2657.
- Copi C. J., Huterer D., Schwarz D. J., Starkman G. D. 2015, *Lack of large-angle TT correlations persists in WMAP and Planck*, MNRAS 451, 2978.
- DESI Collaboration 2025, *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints*, arXiv:2503.14738.
- Hu W., Sugiyama N. 1995, *Anisotropies in the cosmic microwave background: an analytic approach*, ApJ 444, 489.
- Hu W., Fukugita M., Zaldarriaga M., Tegmark M. 2001, *CMB observables and their cosmological implications*, ApJ 549, 669.
- Eisenstein D. J., Hu W. 1998, *Baryonic features in the matter transfer function*, ApJ 496, 605.
- Ma C.-P., Bertschinger E. 1995, *Cosmological perturbation theory in the synchronous and conformal Newtonian gauges*, ApJ 455, 7.
- Blas D., Lesgourgues J., Tram T. 2011, *The Cosmic Linear Anisotropy Solving System (CLASS) II*, JCAP 07, 034.
- Lewis A., Challinor A., Lasenby A. 2000, *Efficient computation of CMB anisotropies in closed FRW models*, ApJ 538, 473.
- Zhao G.-B., Pogosian L., Silvestri A., Zylberberg J. 2009, *Searching for modified growth patterns with tomographic surveys*, Phys. Rev. D 79, 083513.
- Zumalacárregui M., Bellini E., Sawicki I., Lesgourgues J., Ferreira P. G. 2017, *hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System*, JCAP 08, 019.
- Press W. H., Schechter P. 1974, ApJ 187, 425.
- Chevallier M., Polarski D. 2001, IJMPD 10, 213; Linder E. V. 2003, Phys. Rev. Lett. 90, 091301.
- Brouwer M. M. et al. 2021, *The weak lensing radial acceleration relation*, A&A 650, A113.
- Genzel R. et al. 2017, *Strongly baryon-dominated disk galaxies at the peak of galaxy formation ten billion years ago*, Nature 543, 397.
