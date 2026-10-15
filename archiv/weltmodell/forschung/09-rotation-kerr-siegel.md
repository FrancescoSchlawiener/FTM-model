# Kopplung an rotierende Massen: Drehfaltung, T³ und Siegelsche Periodenmatrix

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Umbenennung Spinzahl χ → j (χ₁₀, χ₁₈ bleiben die Igusa-Formen), Schwerepotentiale U_Äq, U_Pol → Φ_Äq, Φ_Pol, Sp(4)-Matrizen U → S, Knotenfrequenz ν_Knoten → ω_Knoten; Sonnenort „Feldform trifft 71 ± 6, Kernform 2σ darunter; die Wahl K/F ist die Modellentscheidung“ → die örtliche Regel ergibt ebenfalls 71, die Feldform liegt an SPARC außerhalb 95 %; Realteil im Kern Δ „offen“ → die Bindung läuft nur über den Betrag des Nomes; ζ* ≤ 5·10⁻⁷ → ≈ 5·10⁻⁷ (2j(v_flach/c)³ = 5,36·10⁻⁷); Im τ_M = r*/(2πr) als Punkt 1 des Autors gekennzeichnet, „Kern ∝ e^{−2π Im τ_M} ≡ 0“ nur in Punkt 1 (in Punkt 3 schaltet χ ab, am ISCO |q_Δs| ≈ 0,14); [G] für J₂, J_Erde, Merkur 42,98″ und die Warp-Präzession. Die Notiz rechnet mit a₀ = 1,2·10⁻¹⁰.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors (Punkt I.2).
> Rechnung: `rechnungen/rotation.py` (121 Prüfungen, alle bestanden; jede Zahl unten stammt aus dieser Rechnung).
> Größen: Δs = GM/(3c²), dτ/dt = √(1−6Δs/r), r₀ = 6Δs, r* = √(3c²Δs/a₀), q = exp(−r*/r) (Rückkehranteil der Punktmasse), τ_M = ω_rad/ω_azim + i·Im τ_M mit Im τ_M = r*/(2πr). [korr.] Das ist Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise.

---

## 0 · Ergebnis in Kürze

1. **Drehfaltung.** Die zweite Faltungsgröße ist die Verwindungslänge a_J = J/(Mc) mit j = a_J/(3Δs). Die Größe GJ/c³ = 3Δs·a_J ist eine *Fläche*.
   - Die Verwindung zeigt sich exakt als Uhrzeit: Prograd läuft ein Umlauf um 4πa_J/c länger als retrograd, und zwar unabhängig vom Radius. Für die Erde sind das 137 ns.
2. **Der FTM-Takt überlebt die Rotation am Äquator unverändert.** Für ruhende Uhren gilt in der Äquatorebene exakt dτ/dt = √(1−6Δs/r), für jedes a_J.
   - Die Nullfläche des Takts ist die Ergofläche. Sie liegt am Äquator immer bei r₀ = 6Δs und schrumpft zum Pol hin auf den Horizont r₊.
3. **T³ statt T².** Ohne Rotation ist ω_θ = ω_φ exakt. Der Vertikalzyklus ist dann an den Azimutzyklus gekoppelt, und die Bahn bleibt effektiv auf T². Rotation trennt beide Zyklen durch die Knotenpräzession, und die Bahn wird zum echten T³.
4. **Siegel-Matrix.** Die natürliche Kopplungsgröße ist die *Knotenzahl* ζ = (ω_φ−ω_θ)/ω_φ. Im schwachen Feld gilt ζ = 2j(3Δs/r)^{3/2}, an der Schwelle ζ* = 2j(v_flach/c)³.
5. **χ₁₀ verschwindet am Diagonalrand, Δ·Δ ist der Koeffizient.** Genauer gilt χ₁₀ = π²ε²·Δ(τ₁)Δ(τ₂)·[1 − (π²/3)E₂(τ₁)E₂(τ₂)ε² + …]. Die Faktorisierung steckt also im Koeffizienten von ε², nicht im Wert selbst.
   - Nahe der Spitze gilt χ₁₀ = q₁q₂·sin²(πε).
   - Rotation geht damit **quadratisch** und unabhängig vom Drehsinn in den Kern ein.
6. **Größenordnung.** Die relativistische Verwindung ändert den Gedächtniskern nicht messbar: An der Schwelle ist ζ* ≈ 5·10⁻⁷ (Milchstraße). [korr.]
   - Die Abplattung von Scheiben liefert dagegen ζ = O(1). Die Siegel-Struktur ist also eine Aussage über **Scheibendicke**, nicht über Frame-Dragging.
7. **Prüfbar.** Kernform (Projektor e_r⊗e_r) und Feldform (∇·[ν∇Φ_N]) geben dieselbe Rotationskurve, aber verschiedene Vertikaldynamik.
   - Am Sonnenort ergibt sich eine Dynamikmasse von 58 bzw. 71–73 M☉/pc². Gemessen sind 71 ± 6. [korr.] Die örtliche Regel a = ν(|g_N|)·g_N ergibt ebenfalls 71; die Feldform liegt an SPARC außerhalb 95 % (senkrechtes Scheibenfeld mit f = +0,03 bis +0,13, ohne IC 2574 0,00 ± 0,03).
8. **Δ im Fernfeld.** Δ verschwindet im Fernfeld auf der imaginären Achse (Re τ_M = 0). Mit Re τ_M = κ/ω_orb = √2 wächst |Δ| wie (Im τ_M)⁻⁶.

---

## 1 · Drehfaltung: wie J den Trichter verwindet

### 1.1 Die zweite Faltungsgröße

| Größe | Formel | Dimension | Bedeutung |
|---|---|---|---|
| Faltungstiefe | Δs = GM/(3c²) | Länge | radiale Zusatzstrecke |
| **Verwindung** | **a_J = J/(Mc)** | Länge | azimutale Zusatzstrecke pro Radiant (siehe 1.4) |
| Spinzahl | j = a_J/(GM/c²) = a_J/(3Δs) | 1 | a_J/Δs = 3j |
| Drehfläche | GJ/c³ = 3Δs·a_J | **Fläche** | geometrisches Mittel √(GJ/c³) ist eine Länge |

Der extreme Kerr-Fall j = 1 bedeutet a_J = 3Δs. Dort fallen beide Horizonte bei r = 3Δs = r₀/2 zusammen.

### 1.2 Kerr-Metrik in FTM-Größen (Boyer-Lindquist, M → 3c²Δs/G)

```
ds² = −(1 − 6Δs·r/Σ) c²dt² − (12Δs·a_J·r·sin²θ/Σ) c dt dφ + (Σ/Δ_K) dr² + Σ dθ²
      + (r² + a_J² + 6Δs·a_J²·r·sin²θ/Σ) sin²θ dφ²

Σ = r² + a_J² cos²θ,     Δ_K = r² − 6Δs·r + a_J²      (Δ_K ≠ Δs ≠ Δ(τ_M))
```

| Struktur | FTM-Form | Prüfung |
|---|---|---|
| Horizonte | r± = 3Δs ± √(9Δs² − a_J²) | Δ_K(r±) = 0 für j = 0,3 / 0,9 / 0,998 ✓ |
| Ergofläche | r_E(θ) = 3Δs + √(9Δs² − a_J² cos²θ) | g_tt(r_E) = 0 bei 30°, 60°, 90° ✓ |
| **Äquator** | **r_E(90°) = 6Δs = r₀ für jedes a_J** | ✓ exakt |
| Pol | r_E(0°) = r₊ | ✓ |
| Mitnahme (exakt) | ω = −g_tφ/g_φφ = 6Δs·a_J·c·r / [(r²+a_J²)² − a_J²Δ_K sin²θ] | – |
| Mitnahme (schwach) | ω = 2GJ/(c²r³) = 6cΔs·a_J/r³ | Verhältnis 0,9904 / 0,99992 / 0,9999992 bei r = 10 / 100 / 1000 · 3Δs ✓ |

### 1.3 Der Takt: Äquator gegen Pol

Eine ruhende Uhr hat den Takt dτ/dt = √(1 − 6Δs·r/Σ). In der Äquatorebene ist Σ = r², also gilt **exakt die FTM-Formel √(1 − 6Δs/r), unabhängig von a_J**. Die Rechnung prüft das bei j = 0,5 und 0,9 auf alle Stellen. Am Pol gilt dagegen √(1 − 6Δs·r/(r² + a_J²)).

Lesart: Die Verwindung ändert die Takttiefe *entlang der Drehachse*, nicht in der Drehebene. Die Nullfläche des Takts wandert von r₀ (Äquator) zu r₊ (Pol).

j = 0,9, gleiche Boyer-Lindquist-Koordinate r (Einheiten 3Δs; gilt für Sgr A* und M87* gleich):

| r/3Δs | ruhend, Äquator | ruhend, Pol | ZAMO, Äquator | ZAMO, Pol | Pol/Äquator |
|---|---|---|---|---|---|
| 1,5 | – (Ergoregion) | 0,1400 | 0,1204 | 0,1400 | – |
| 2,0 | 0 (= r₀) | 0,4104 | 0,3796 | 0,4104 | – |
| 2,5 | 0,4472 | 0,5402 | 0,5170 | 0,5402 | 1,208 |
| 3 | 0,5774 | 0,6232 | 0,6067 | 0,6232 | 1,079 |
| 4 | 0,7071 | 0,7239 | 0,7154 | 0,7239 | 1,024 |
| 6 | 0,8165 | 0,8210 | 0,8180 | 0,8210 | 1,0055 |
| 10 | 0,8944 | 0,8953 | 0,8946 | 0,8953 | 1,0010 |
| 20 | 0,9487 | 0,9488 | 0,9487 | 0,9488 | 1,0001 |

Die Spalte „ZAMO" gibt den Takt eines Beobachters ohne Drehimpuls an (α = √(ΣΔ_K/A)); er existiert auch innerhalb der Ergoregion. Auf der Achse gilt g_tφ = 0, dort stimmen ZAMO-Takt und statischer Takt überein ✓.

**Erde, mitrotierende Uhren an der Oberfläche:** Die Einzelbeiträge liegen bei 10⁻¹²:
- Fliehkraft: 1,20·10⁻¹²
- Abplattung: 2,34·10⁻¹²

Sie heben sich jedoch auf. Die Schwerepotentiale Φ_Äq = 62 636 796,0 m²/s² und Φ_Pol = 62 636 701,4 m²/s² ergeben eine Taktdifferenz von 1,05·10⁻¹⁵. Dieser Rest stammt aus dem vernachlässigten J₄ bzw. der Geoidhöhe. Das Geoid ist eine Äquipotentialfläche, die Uhren gehen also gleich.

Die Kerr-eigenen Anteile sind:
- Mitnahme: 6Δs·a_J·ω_orb/(cR) = 1,1·10⁻²¹
- a_J²-Term: 3Δs·a_J²/R³ = 1,8·10⁻²²

Die Verwindung der Erde ist mit Uhren an der Oberfläche nicht zugänglich.

### 1.4 Die Verwindung als Uhrzeit (exakt)

Für äquatoriale Kreisbahnen gilt ω_φ = c³/(GM(r̃^{3/2} ± j)) (Bardeen, Press & Teukolsky 1972). Daraus folgt:

```
T_φ = 2π·r^{3/2}/√(GM)  ±  2π·a_J/c           ⇒     T_pro − T_retro = 4πa_J/c = 4πJ/(Mc²)
```

Die Differenz hängt nicht vom Radius ab; die Rechnung prüft sie bei r̃ = 6, 20 und 10⁴. Pro Radiant Umlauf addiert der verwundene Trichter also die Lichtlaufstrecke a_J. Das ist das azimutale Gegenstück zur radialen Zusatzstrecke Δs.

In der Literatur heißt das der **gravitomagnetische Uhreffekt** (Cohen & Mashhoon 1993). Er ist noch nicht gemessen.

### 1.5 Zahlen

| Körper | Δs | a_J = J/(Mc) | j | 4πa_J/c |
|---|---|---|---|---|
| Erde (J = 5,86·10³³ [G]) | 1,478 mm | 3,273 m | 738 | 137 ns |
| Sonne (J = 1,92·10⁴¹) | 492,3 m | 322,0 m | 0,218 | 13,5 µs |
| Sgr A* (4,3·10⁶ M☉, j = 0,9 angenommen) | 2,117·10⁹ m | 5,716·10⁹ m | 0,9 | 240 s |
| M87* (6,5·10⁹ M☉, j = 0,9 angenommen) | 3,200·10¹² m | 8,640·10¹² m | 0,9 | 3,6·10⁵ s |
| Milchstraße (J ≈ 2M·R_d·v, 220 km/s) | 2,954·10¹³ m | 1,18·10¹⁷ m (3,8 pc) | 1 329 | 156 a |

Erde, Sonne und Milchstraße sind **keine Kerr-Körper**. Das zeigt der Quadrupol:
- Kerr verlangt J₂ = (a_J/R)² = 2,6·10⁻¹³ (Erde) bzw. 2,1·10⁻¹³ (Sonne).
- Gemessen sind 1,08·10⁻³ bzw. ≈ 2,2·10⁻⁷ [G].

Universell ist nur der Dipolanteil der Verwindung in Ordnung a_J¹, also Lense-Thirring. Höhere Formmomente gehören zum Material, nicht zur Faltung.

**Präzession im schwachen Feld:**

| Größe | Formel (FTM) | Rechnung | Messung |
|---|---|---|---|
| GP-B geodätisch (r = 7027,4 km) | (3/2)(GM)^{3/2}/(c²r^{5/2}) = (9/2)(Δs/r)·ω_K | 6 603 mas/a | 6 601,8 ± 18,3 (Everitt et al. 2011) |
| GP-B Mitnahme | GJ/(2c²r³) = 3cΔs·a_J/(2r³) (Bahnmittel, polar) | 40,8 mas/a; × cos 16,84° (Leitstern) = 39,1 | 37,2 ± 7,2 (Vorhersage 39,2) |
| LAGEOS-Knoten | 2GJ/(c²a³(1−e²)^{3/2}) | 30,67 mas/a | Ciufolini & Pavlis 2004: Summe beider Satelliten, ≈ 99 % der Vorhersage |
| LAGEOS 2 | ebenso | 31,49 mas/a | ebenso |
| Merkur, Perihel (Schwarzschild) | 18πΔs/(a(1−e²)) pro Umlauf | 42,99″/Jhdt | 42,98″/Jhdt [G] |
| Merkur, Lense-Thirring | Knoten bzw. Perihellänge | +1,02 bzw. −2,04 mas/Jhdt | – |
| Sonnenoberfläche | ω = 2GJ/(c²R³) | 8,47·10⁻¹³ rad/s | – |
| Milchstraße bei 15 kpc | 2GJ/(c²r³) | 1,9·10⁻⁶ km/s/kpc | Warp-Präzession ≈ 10,9 km/s/kpc (Poggio et al. 2020) [G] |

Für einen Kreisel gilt die allgemeine Form ω_LT = (G/c²r³)·[3(J·r̂)r̂ − J] (Schiff 1960). Der Präzessionswinkel einer Bahn pro Umlauf ist **2πζ** mit der Knotenzahl ζ aus Abschnitt 2, im schwachen Feld 4πj(3Δs/r)^{3/2}.

**Schwarze Löcher (j = 0,9):**

| | Sgr A* | M87* |
|---|---|---|
| r₊ | 9,12·10⁹ m (1,4359·3Δs) | 1,379·10¹³ m |
| r₋ | 3,58·10⁹ m | 5,42·10¹² m |
| Ergofläche Äquator = r₀ = 6Δs | 1,270·10¹⁰ m | 1,920·10¹³ m |
| ω_H = c·a_J/(r₊² + a_J²) | 1,479·10⁻² rad/s (7,08 min) | 9,79·10⁻⁶ rad/s (7,43 d) |
| ISCO prograd / retrograd | 2,3209 / 8,7174 · 3Δs | ebenso |
| Umlaufzeit am ISCO | 9,84 min | 10,33 d |

Am ISCO verschwindet κ² exakt ✓. Der Spin 0,9 ist eine **Annahme**; für Sgr A* ist er nicht sicher bestimmt.

---

## 2 · Präzedierende Eigenfrequenzen und der Bahntorus T³

### 2.1 Die drei Grundfrequenzen (äquatoriale Kreisbahn, prograd, r̃ = r/(3Δs))

```
ω_φ = c³ / (GM (r̃^{3/2} + j))
ω_r² = ω_φ² (1 − 6/r̃ + 8j r̃^{−3/2} − 3j²/r̃²)          (radiale Epizyklenfrequenz κ)
ω_θ² = ω_φ² (1 − 4j r̃^{−3/2} + 3j²/r̃²)                (vertikale Frequenz)
ω_Knoten = ω_φ − ω_θ,   ω_Periapsis = ω_φ − ω_r
```

Die Formeln stammen aus Okazaki, Kato & Fukue 1987, Kato 1990 und Wagoner 1999. Retrograde Bahnen erhält man mit j → −j.

**Unabhängige numerische Bestätigung**, jede gegen die Formel:
- **A:** ω_φ aus der Kreisbahnbedingung ∂_r g_tt + 2ω_φ∂_r g_tφ + ω_φ²∂_r g_φφ = 0. Abweichung ≤ 6·10⁻¹¹.
- **B:** ω_r und ω_θ aus der Linearisierung der Hamiltonschen Geodätengleichung (H = ½g^{μν}p_μp_ν, Ableitungen numerisch, geteilt durch u^t). Abweichung ≤ 2·10⁻⁷ für (j, r̃) = (0; 8), (0,5; 7), (0,9; 5), (0,9; 12), (−0,7; 12).
- **C:** volle nichtlineare Geodäte mit RK4 bei j = 0,9, r̃ = 8, Frequenzen aus Nulldurchgängen. Abweichung ≤ 3·10⁻⁶. Ergebnis: (ω_φ, ω_r, ω_θ) = (0,042504; 0,030950; 0,039846)·c³/GM.

**Grenzfälle:**
- j = 0: **ω_θ = ω_φ exakt**, bei jedem Radius ✓.
- Schwaches Feld: Knoten = 2GJ/(c²r³)·[1 − (3/4)j r̃^{−1/2}], Periapsis = 3(GM)^{3/2}/(c²r^{5/2})·[1 − (4/3)j r̃^{−1/2}]. Beides geprüft bei r̃ = 10⁴ bis 10⁶. Der führende Term ist Lense-Thirring bzw. die Einstein-Periheldrehung.

### 2.2 Vom T² zum T³

Eine Bahn um eine nichtrotierende Masse hat zwar drei Winkel, aber ω_θ = ω_φ. Die Vertikalschwingung ist nur eine gekippte Bahnebene, der θ-Zyklus ist an den φ-Zyklus phasenverriegelt. Nach der Regel „zwei phasenverriegelte Frequenzen bilden einen einzigen Zyklus" ist der invariante Torus effektiv ein **T²**.

Rotation hebt die Verriegelung auf, weil der Knoten mit ω_φ − ω_θ ≠ 0 wandert. Die Bahn liegt dann auf einem echten **T³** mit drei rational unabhängigen Frequenzen.

Die Kerr-Geodäten sind integrabel (Carter 1968). Die drei Wirkungen bleiben erhalten, und es gibt kein Chaos. Resonanzen treten nur auf, wenn die Frequenzverhältnisse rational sind. Beispiel: Die 3:2-Resonanz ω_θ : ω_r liegt bei j = 0,9 bei r̃ = 4,96263 (vgl. Abramowicz & Kluźniak 2001):
- Sgr A*: f_θ = 5,53·10⁻⁴ Hz, f_r = 3,68·10⁻⁴ Hz, f_Knoten = 7,6·10⁻⁵ Hz
- M87*: f_θ = 3,66·10⁻⁷ Hz, f_r = 2,44·10⁻⁷ Hz

Frequenzverhältnisse bei j = 0,9 (Knotenzahl ζ = 1 − ω_θ/ω_φ):

| r̃ | ω_r/ω_φ | ω_θ/ω_φ | ζ | 2j r̃^{−3/2} | ω_r/ω_φ (j = 0) |
|---|---|---|---|---|---|
| 2,321 (ISCO) | 0 | 0,6580 | 0,3420 | 0,5091 | – |
| 4 | 0,4981 | 0,8378 | 0,1622 | 0,2250 | – |
| 6 | 0,6499 | 0,9069 | 0,0931 | 0,1225 | – (ISCO) |
| 10 | 0,7768 | 0,9542 | 0,04582 | 0,05692 | 0,6325 |
| 20 | 0,8800 | 0,9828 | 0,01724 | 0,02013 | 0,8367 |
| 100 | 0,9731 | 0,9983 | 0,001680 | 0,001800 | 0,9695 |

### 2.3 Konstruktion der Periodenmatrix Ω ∈ ℍ₂ (Setzung)

In Analogie zu τ_M = ω_rad/ω_azim + i·Im τ_M wird φ als Takt gewählt. Die beiden übrigen Zyklen liefern je eine Rotationszahl, dazu kommt die Knotenzahl als Kopplung:

```
        ⎛ ω_r/ω_φ + it        ζ          ⎞
Ω  =    ⎜                                ⎟ ,     ζ = (ω_φ − ω_θ)/ω_φ ,    t = Im τ_M = r*/(2πr)
        ⎝      ζ          ω_θ/ω_φ + it   ⎠
```

[korr.] t = r*/(2πr) ist Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3, t = ln(r/Δs)/2π). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise.

Eigenschaften, alle numerisch geprüft:
- (i) Ω = Ωᵀ und Im Ω = t·𝟙 > 0.
- (ii) j → 0 ⇒ ζ → 0. Ω wird diagonal, liegt also am **Rand der Zerfällung** (Produkt zweier elliptischer Tori).
- (iii) a_J → −a_J ⇒ ζ → −ζ. Das ist die Sp(4,ℤ)-Transformation S = diag(1, −1), und χ₁₀ ist darunter invariant ✓. Der Drehsinn ist also keine Strukturgröße.
- (iv) Im schwachen Feld gilt ζ = 2j(3Δs/r)^{3/2}. An der Schwelle r* folgt **ζ* = 2j(v_flach/c)³ = 2j·δ(r*)^{3/2}**, mit δ(r*) = 3Δs/r* = (v_flach/c)². Die Identität ist für alle Körper auf 10⁻¹⁵ geprüft.

**Vorbehalt:** Die Kerr-Bewegung ist separabel (Carter). Die linearen r- und θ-Epizyklen koppeln nicht. Die Nebendiagonale ζ ist deshalb **eine Setzung**. Sie ist die einzige Kopplungsgröße erster Ordnung im Problem: der gemeinsame φ-Zyklus, den Knotenlauf und Vertikalzyklus teilen (γ_θ = γ_φ − γ_Knoten).

ℍ₂ statt ℍ₃: Neben dem Takt φ bleiben zwei Rotationszahlen. Für ℍ₃ bräuchte man eine dritte unabhängige Zahl, siehe offene Punkte.

---

## 3 · Δ(τ_M) für präzedierende Tori: die Igusa-Form χ₁₀

### 3.1 Definition und Modularität

Die Thetakonstanten mit Charakteristik lauten θ[a;b](Ω) = Σ_{n∈ℤ²} exp(πi(n+a/2)ᵀΩ(n+a/2) + πi(n+a/2)ᵀb). Davon sind genau 10 gerade (a·b ≡ 0 mod 2). Nach Igusa (1962) ist

```
χ₁₀(Ω) = −2⁻¹⁴ · ∏_{10 gerade m} θ_m(Ω)²          (Gewicht 10, Spitzenform für Sp(4,ℤ))
```

Geprüft am Punkt Ω₀ = ((0,3+1,1i, 0,2+0,3i), (0,2+0,3i, −0,1+0,9i)), jeweils auf 10⁻¹⁵:
- χ₁₀(−Ω⁻¹) = det(Ω)¹⁰χ₁₀(Ω)
- Translationen in τ₁ und in der Nebendiagonale
- SΩSᵀ für S = ((1,0),(1,1)) und S = diag(1,−1)

Grundlagen aus Geschlecht 1, ebenfalls geprüft (η ist hier die Dedekind-Funktion):
- θ₂θ₃θ₄ = 2η³
- θ₁′(0) = 2πη³
- Δ(−1/τ_M) = τ_M¹²Δ(τ_M)
- d log Δ/dτ_M = 2πi·E₂

### 3.2 Verhalten am Rand: χ₁₀ verschwindet, Δ(τ₁)Δ(τ₂) ist der Koeffizient

**Diagonal (nichtrotierend, ζ = 0):** θ[11;11] ist ungerade × ungerade und verschwindet exakt (numerisch 1,4·10⁻¹⁶). Die übrigen neun Faktoren ergeben ∏θ² = 2¹²·η(τ₁)¹⁸η(τ₂)¹⁸ ✓. Also gilt

```
χ₁₀(diag(τ₁, τ₂)) = 0 .
```

**χ₁₀ verschwindet auf dem Diagonalrand quadratisch, und der Koeffizient ist das Produkt:**

```
χ₁₀( τ₁  ε ; ε  τ₂ ) = π² ε² Δ(τ₁) Δ(τ₂) · [ 1 − (π²/3) E₂(τ₁) E₂(τ₂) ε² + O(ε⁴) ]
```

Herleitung der führenden Ordnung: θ[11;11] ≈ 2πiε·(iη₁³)(iη₂³) = −2πiε η₁³η₂³. Mit dem Rest 2¹²η₁¹⁸η₂¹⁸ und −2⁻¹⁴ ergibt sich π²ε²η₁²⁴η₂²⁴.

Die nächste Ordnung folgt aus dem Gewicht (14, 14) des ε⁴-Koeffizienten: Der einzige quasimodulare Faktor vom Gewicht (2, 2) ist E₂⊗E₂. Die Konstante −π²/3 liefert der Fourier-Jacobi-Grenzfall.

Numerisch bei τ₁ = 0,1+0,8i, τ₂ = −0,2+1,0i:

| ε | Q_S/ε² mit Q_S = χ₁₀/(π²Δ₁Δ₂) |
|---|---|
| 10⁻¹ | 0,97170 + 0,00178i |
| 10⁻² | 0,999716 + 0,0000187i |
| 10⁻⁴ | 0,99999997 |

Der Koeffizient der nächsten Ordnung ergibt nach Richardson 2,8419304 − 0,1873539i. Der Vergleichswert (π²/3)E₂E₂ ist 2,8419301 − 0,1873540i ✓.

**Weitere Randformen (geprüft):**

| Rand | Form | Prüfgenauigkeit |
|---|---|---|
| Spitze q₁, q₂ → 0 | χ₁₀ = q₁q₂·sin²(πε) + O(q³) | 2·10⁻⁷ bei Im τ_M = 3 |
| Ein Zyklus fällt weg (q₂ → 0, Fourier-Jacobi) | χ₁₀ → ¼·q₂·η(τ₁)¹⁸·θ₁(τ₁, z)² | 3·10⁻¹⁰ |
| Symmetrie | χ₁₀(ε) = χ₁₀(−ε) = χ₁₀(ε + 1) | Rechengenauigkeit |
| ε ∈ ℤ | χ₁₀ = 0 | – |

η¹⁸θ₁² ist die Jacobi-Spitzenform vom Gewicht 10 und Index 1 (Eichler & Zagier 1985).

### 3.3 Was daraus für den Gedächtniskern folgt

**Neue Größe:** Die Kopplung der Präzessionsfrequenz misst die **Kopplungszahl**

```
Q_S(Ω) := χ₁₀(Ω) / (π² Δ(τ_r) Δ(τ_θ))  =  ζ² + O(ζ⁴),        im Nahfeld (q → 0) exakt:  Q_S = sin²(πζ)/π²
```

mit der Knotenzahl ζ = ω_Knoten/ω_φ. Q_S ist der Siegel-Ausdruck dafür, wie stark der Radialzyklus und der Vertikalzyklus „verklebt" sind.

Der rotierende Kern lautet damit (Vorschlag):

```
K_rot = κ₀ · e^{−λ_mΔt} · [ M_k(τ_r) e_r⊗e_r  +  π² Q_S(Ω) · M_k(τ_r) Δ(τ_θ) · 𝒯 ]
```

Dabei ist π²Q_SΔ_rΔ_θ = χ₁₀, und 𝒯 ist eine noch festzulegende Tensorstruktur in der (e_r, e_z)-Ebene (offener Punkt 1).

Diese Folgerungen hängen nicht von 𝒯 ab:

1. **Quadratisch und ohne Drehsinn.** Rotation ändert den Kern erst in Ordnung ζ², und nur über den Betrag der Präzession. Das folgt aus der Sp(4,ℤ)-Invarianz unter ζ → −ζ.
2. **Resonanzen schalten die Kopplung aus.** Bei ζ ∈ ℤ, also ω_θ/ω_φ ∈ ℤ, zerfällt der Torus wieder, und die Kopplung verschwindet. Das ist eine Folge der Setzung und damit prüfbar, sofern sie gilt.
3. **Größenordnung Kerr.** Für Schwarze Löcher ist Im τ_M = r*/(2πr) ~ 10⁷ (Sgr A* am ISCO: 2,4·10⁷). Dann gilt Q_S = sin²(πζ)/π² exakt, der Kern selbst ist aber ∝ e^{−2π Im τ_M} ≡ 0. [korr.] Das gilt in Punkt 1 des Autors (Im τ_M = r*/(2πr)). In Punkt 3 (|q_Δs| = Δs/r) ist der Nome am ISCO nicht klein (|q_Δs| ≈ 0,14); dort schaltet der Schalter χ(s) = s/(e^s − 1) mit s = r*/r ~ 10⁸ den Kern ab. Gedächtniskern und Kerr-Präzession leben in **getrennten Regimen**. An der Schwelle, wo der Kern wirkt, ist ζ* = 2j(v_flach/c)³:

| System | ζ* |
|---|---|
| Erde | 1,8·10⁻¹⁹ |
| Sonne | 7,3·10⁻¹⁹ |
| Sgr A* | 2,8·10⁻¹³ |
| M87* | 6,9·10⁻¹¹ |
| Milchstraße | 5,4·10⁻⁷ |

   Die relativistische Verwindung verändert den Gedächtniskern also nirgends messbar; Q_S ≈ ζ*² ≤ 3·10⁻¹³.

4. **Fernfeld wird nichtperturbativ.** Die Entwicklung gilt nur für ε² ≪ 3/(π²|E₂E₂|):
   - Bei Re τ_M = 1,41 / 1,0 ist (π²/3)|E₂E₂| = 54 (Im τ_M = 0,3), 2 658 (Im τ_M = 0,1) und 25 940 (Im τ_M = 0,05). Bei Im τ_M = 0,05 und ε = 10⁻³ ist der ε⁴-Term bereits sichtbar (numerisch).
   - In Scheiben (Abschnitt 4) springt Q_S von Radius zu Radius um Größenordnungen, zwischen 3·10⁻¹⁷ und 3·10⁻².
   - Siegel- und Modulformen sind deshalb **Nahfeldformen**. Das Fernfeld trägt die Rückkopplungssumme ν = 1/(1−q).

### 3.4 Δ im Fernfeld bei reellem Anteil von τ_M

Δ ist an beiden Enden null für **Re τ_M = 0**. Mit τ_M = ω_rad/ω_azim + it, t = Im τ_M, ist bei flacher Rotationskurve Re τ_M = κ/ω_orb = √2:

| t | \|Δ(it)\| | \|Δ(√2 + it)\| | t⁶·\|Δ(√2 + it)\| |
|---|---|---|---|
| 1 | 1,8·10⁻³ | 1,9·10⁻³ | 1,9·10⁻³ |
| 0,3 | 1,5·10⁻³ | 2,35 | 1,7·10⁻³ |
| 0,1 | 5,2·10⁻¹⁶ | 1,05·10³ | 1,0·10⁻³ |
| 0,03 | 2,1·10⁻⁷³ | 2,4·10⁶ | 1,7·10⁻³ |

Die Größe (Im τ_M)⁶|Δ(τ_M)| ist SL(2,ℤ)-invariant (geprüft) und bleibt beschränkt. Bei irrationalem Realteil **wächst |Δ| im Fernfeld wie t⁻⁶**. [korr.] Für den Kern ist das nicht offen: Die Bindung läuft nur über den Betrag des Nomes; ginge die Phase (der Realteil) mit ein, fiele die Verstärkung außen unter 1. Der Realteil Re τ_M = ω_rad/ω_azim ist die Umlaufphase.

---

## 4 · Folgerungen für rotierende Scheiben

Modell: exponentielle Scheibe mit M = 6·10¹⁰ M☉, R_d = 2,6 kpc (Bland-Hawthorn & Gerhard 2016) und Höhe h. Die Rotationskurve folgt der Freeman-Formel mit Bessel-Funktionen; diese sind numerisch gegen Tabellenwerte auf 10⁻¹⁰ geprüft. Übergang: ν = 1/(1 − q) mit q = exp(−√(a_N/a₀)).

### 4.1 Hängt der Operator von der Scheibendicke ab? — Das hängt von seiner Lesart ab

Das Modell enthält zwei Formen des Operators:

| | Kernform K (Projektor e_r⊗e_r, Kugel-Radialrichtung) | Feldform F (∇²Φ = ∇·[ν∇Φ_N]) |
|---|---|---|
| Rotationskurve ω_orb², κ² | ν·ω_N² | ν·ω_N² (identisch) |
| Vertikalfrequenz | ν_z² = ν_z,N² + (ν−1)·ω_N² | **ν_z² = ν·ν_z,N²** (exakt in der Mittelebene bei g_R = ν g_N,R; geprüft auf 10⁻⁷) |
| erhalten bleibt | ν_z² − ω_orb² (Newton) ✓ | ν_z/ω_orb (Newton) ✓ |
| Abhängigkeit von der Dicke | Zusatzterm ist dickenunabhängig, Abplattung wird verdünnt | volle Abplattung bleibt erhalten |
| Siegel-Kopplung ζ = 1 − ν_z/ω_orb | → 0 außen (der Kern macht den Torus kugelig und zerfällt selbst) | bleibt auf Newton-Wert |

Frequenzen bei h = 0,3 kpc (ω_orb und ω_p = ω_orb − ν_z in km/s/kpc):

| R [kpc] | ν | v [km/s] | ω_orb | κ/ω_orb | ν_z/ω_orb N = F | ν_z/ω_orb K | ω_p N | ω_p K | ω_p F |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 1,267 | 213 | 53,3 | 1,589 | 3,434 | 3,085 | −115,2 | −111,1 | −129,7 |
| 8,2 | 1,522 | 230 | 28,1 | 1,378 | 3,322 | 2,756 | −52,8 | −49,3 | −65,1 |
| 12 | 1,874 | 220 | 18,3 | 1,289 | 2,833 | 2,179 | −24,6 | −21,6 | −33,6 |
| 16 | 2,325 | 209 | 13,0 | 1,274 | 2,212 | 1,635 | −10,4 | −8,3 | −15,8 |
| 20 | 2,812 | 201 | 10,0 | 1,293 | 1,683 | 1,285 | −4,1 | −2,9 | −6,8 |
| 25 | 3,431 | 194 | 7,8 | 1,323 | 1,277 | 1,088 | −1,2 | −0,7 | −2,2 |
| 30 | 4,045 | 190 | 6,4 | 1,345 | 1,105 | 1,027 | −0,3 | −0,2 | −0,7 |

Dickenabhängigkeit bei 8,2 kpc, ν_z/ω_orb:

| h [kpc] | N = F | K |
|---|---|---|
| 0,15 | 4,64 | 3,81 |
| 0,3 | 3,32 | 2,76 |
| 0,6 | 2,40 | 2,04 |

Die Knotenzahl der Scheibe ζ = 1 − ν_z/ω_orb liegt zwischen −2,4 und −0,3. Sie ist damit um **sechs (Milchstraße) bis neunzehn (Erde) Größenordnungen größer** als die relativistische Knotenzahl. Wenn die Siegel-Struktur irgendwo physikalisch wirkt, dann über die Abplattung.

### 4.2 Prüfbarer Effekt 1: Vertikale Dynamikmasse am Sonnenort

Gesucht ist die Dynamikmasse in |z| < 1,1 kpc bei Σ_bar = 47,1 M☉/pc² (McKee, Parravano & Hollenbach 2015):

| g_N,R am Sonnenort | Kernform K | Feldform F |
|---|---|---|
| Scheibenmodell 1,374·10⁻¹⁰ m/s² | 58,1 M☉/pc² | 70,8 M☉/pc² |
| Punktmasse 1,244·10⁻¹⁰ m/s² | 57,9 M☉/pc² | 72,6 M☉/pc² |

Gemessen sind 71 ± 6 M☉/pc² (Kuijken & Gilmore 1991). **Die Feldform trifft den Wert. Die Kernform liegt etwa 2σ darunter.** [korr.] Die örtliche Regel a = ν(|g_N|)·g_N (Projektor längs n am Ort der Bahn) ergibt ebenfalls 71 M☉/pc². Die Feldform scheidet aus: An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %.

Die Rechnung ist eine Abschätzung: Die Feldform ist algebraisch genähert, und das Modell enthält keinen Bulge. [korr.] Die Wahl zwischen K und F ist nicht die Modellentscheidung: Die Bahnen folgen der örtlichen Regel; die Metrik (Feldform) trägt nur deren wirbelfreien Teil. Für Kugeln sind alle drei Formen gleich.

### 4.3 Prüfbarer Effekt 2: Knotenpräzession und Warps im Tief-a₀-Regime

- **Vorzeichen.** Ein gekippter Ring präzediert frei mit ω_p = ω_orb − ν_z < 0, also rückläufig, und zwar in allen drei Varianten. Die mit Gaia gemessene Warp-Präzession der Milchstraße ist dagegen rechtläufig, +10,9 km/s/kpc (Poggio et al. 2020, Interpretation umstritten) [G]. Weder freie Präzession noch relativistische Mitnahme erklärt sie; die Mitnahme liefert 1,9·10⁻⁶ km/s/kpc, also einen Faktor 6·10⁶ zu wenig. Das Modell braucht dafür ein äußeres Drehmoment, wie jede andere Theorie auch.
- **Aufwickeln.** Die Differentialpräzession zwischen 12 und 20 kpc wickelt einen Warp in folgenden Zeiten auf:

| Variante | Aufwickelzeit |
|---|---|
| Kernform K | 0,33 Gyr |
| Newton N | 0,30 Gyr |
| Feldform F | 0,23 Gyr |

  Das FTM löst das klassische Aufwickelproblem also *nicht* von selbst. Die Kernform wickelt langsamer, um den Faktor 1,43 gegenüber der Feldform.
- **Unterscheidbar ist das Verhältnis der Präzessionsraten.** Es gilt ω_p^F/ω_p^K = 2,4 bei 20 kpc und 3,2 bei 25 kpc. Asymptotisch ist es ≈ ν(y), das folgt aus den erhaltenen Größen in 4.1. Äußere gekippte Ringe mit messbarer Präzession, etwa Polarringe oder verbogene HI-Scheiben, trennen beide Lesarten. Tief im a₀-Regime wird der Unterschied umso größer, je kleiner a_N/a₀ ist.

### 4.4 Siegel-Größen der Scheibe (h = 0,3 kpc, r* = 8,349 kpc)

| R [kpc] | t = Im τ_M | ζ_K | Q_S,K | ζ_F | Q_S,F |
|---|---|---|---|---|---|
| 4 | 0,332 | −2,085 | 8,1·10⁻³ | −2,434 | 8,5·10⁻⁴ |
| 8,2 | 0,162 | −1,756 | 1,5·10⁻² | −2,322 | 4,2·10⁻⁶ |
| 12 | 0,111 | −1,179 | 1,6·10⁻⁵ | −1,833 | 3,3·10⁻² |
| 16 | 0,083 | −0,635 | 4,9·10⁻⁴ | −1,212 | 2,1·10⁻¹⁰ |
| 20 | 0,066 | −0,285 | 3,0·10⁻¹⁷ | −0,683 | 2,8·10⁻⁴ |

Q_S ist bei Im τ_M ≲ 0,3 keine glatte Funktion des Radius, sondern springt um Größenordnungen (siehe 3.3, Punkt 4). **Als Fernfeldkern für Scheiben ist χ₁₀ ungeeignet.** Er bleibt eine Nahfeld- und Strukturaussage.

---

## 5 · Lesart nach der Beobachter-Regel

- Die Verwindung a_J ist eine Länge, keine Substanz. Sie wirkt als Zusatzlaufzeit 2πa_J/c pro Umlauf und als Nebendiagonale g_tφ. Frame-Dragging ist deshalb keine Kraft eines Objekts, sondern die Schrägstellung der Zeitachse gegen den Azimut.
- ζ und Q_S sind Zahlen der Periodenmatrix, also reine Struktur. Die Masse ρ_bar bleibt in allen Varianten die einzige Quelle.
- Kern- und Feldform unterscheiden sich allein darin, **welche Richtung** der Operator als „radial" liest: die Kugelrichtung e_r oder die Feldrichtung ∇Φ_N. Das ist eine geometrische Wahl, und sie ist messbar (4.2 und 4.3). [korr.] Im Modell gilt keine der beiden als Bahnregel, sondern die örtliche Regel a = ν(|g_N|)·g_N: Sie liest die Feldrichtung am Ort der Bahn, ohne die wirbelfreie Ergänzung der Feldform; am Sonnenort gibt sie 71 M☉/pc², die Feldform liegt an SPARC außerhalb 95 %.

---

## 6 · Offene Punkte, nächste Rechenschritte

1. **Tensorform 𝒯 des Rotationskerns.** Kandidat ist der Gradient ∂χ₁₀/∂Ω_ij als symmetrischer Tensor ↔ e_i⊗e_j. Er ist nur auf dem Divisor χ₁₀ = 0 kovariant. Dort ist die gemischte Komponente ∝ ζ, also linear, und koppelt Radialgedächtnis an Vertikalbewegung. Zu prüfen ist, ob daraus ein Kippmoment (Warp-Antrieb) entsteht und ob es mit der Sp(4,ℤ)-Invarianz verträglich ist.
2. **Fernfeld-Analogon der Rückkopplungssumme** für zwei Zyklen, z. B. ν₂ = 1/((1−q_r)(1−q_θ)) mit einem Kopplungsterm in ζ. Bedingung: Es muss sich bei ζ = 0 auf ν = 1/(1−q) reduzieren und die RAR erhalten.
3. **Realteil von τ_M im Kern Δ festlegen:** i·r*/(2πr) oder κ/ω_orb + i·Im τ_M. Davon hängt ab, ob Δ im Fernfeld verschwindet oder wie (Im τ_M)⁻⁶ wächst (3.4). [korr.] Festgelegt: Die Bindung läuft nur über den Betrag des Nomes; ginge die Phase mit ein, fiele die Verstärkung außen unter 1.
4. **Geschlecht 3 (ℍ₃).** Man nimmt den Knotenzyklus bzw. den Eigenspin der Zentralmasse als dritte Rotationszahl hinzu. Die Spitzenform ist dann χ₁₈, das Produkt der 36 geraden Thetakonstanten; sie verschwindet auf dem hyperelliptischen Ort (Igusa 1967). Zu prüfen ist das Randverhalten gegen χ₁₀·Δ.
5. **Nichtäquatoriale Kerr-Bahnen** (Carter-Konstante): exakte Frequenzen für beliebige Neigung nach Schmidt (2002). Dazu ein Test der Setzung ζ bei endlicher Neigung.
6. **Scheiben voll rechnen:** QUMOND-Poisson-Lösung statt algebraischer Näherung, Eigengravitation des Warps (modifizierte Kippmoden, Sparke & Casertano 1988). Danach Vergleich mit Polarring-Galaxien und mit σ_z-Messungen externer Scheiben, um K gegen F zu entscheiden. [korr.] Im Modell ist das entschieden: Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N, die Feldform liegt an SPARC außerhalb 95 %; für die örtliche Regel sind die Präzessionsraten geneigter Ringe noch nicht gerechnet.
7. **KAM-Bedingung.** Die Kolmogorov-Nichtentartung det(∂²H/∂J_i∂J_j) ≠ 0 ist bei j = 0 verletzt, weil H nur von J_θ + J_φ abhängt. Die Hesse-Matrix aus den Kerr-Wirkungen ist numerisch zu berechnen. Zu prüfen ist, ob ihre Nebendiagonale eine *dynamisch begründete* Alternative zur Setzung ζ liefert.
8. **Direkte Messung der Verwindung:** Der gravitomagnetische Uhreffekt (137 ns pro Erdumlauf) verlangt die Bahnbestimmung prograder und retrograder Satelliten auf Millimeterbruchteile. Er wäre der sauberste Test von a_J als Faltungsgröße.

---

## Literatur

- R. P. Kerr, Phys. Rev. Lett. 11, 237 (1963); R. H. Boyer, R. W. Lindquist, J. Math. Phys. 8, 265 (1967)
- B. Carter, Phys. Rev. 174, 1559 (1968)
- J. M. Bardeen, W. H. Press, S. A. Teukolsky, ApJ 178, 347 (1972)
- J. Lense, H. Thirring, Phys. Z. 19, 156 (1918); L. I. Schiff, Phys. Rev. Lett. 4, 215 (1960)
- C. W. F. Everitt et al., Phys. Rev. Lett. 106, 221101 (2011) — Gravity Probe B
- I. Ciufolini, E. C. Pavlis, Nature 431, 958 (2004) — LAGEOS
- J. M. Cohen, B. Mashhoon, Phys. Lett. A 181, 353 (1993) — gravitomagnetischer Uhreffekt
- A. T. Okazaki, S. Kato, J. Fukue, PASJ 39, 457 (1987); S. Kato, PASJ 42, 99 (1990); R. V. Wagoner, Phys. Rep. 311, 259 (1999)
- M. A. Abramowicz, W. Kluźniak, A&A 374, L19 (2001) — 3:2-Resonanz
- W. Schmidt, Class. Quantum Grav. 19, 2743 (2002) — Frequenzen beliebiger Kerr-Bahnen
- J. Igusa, Amer. J. Math. 84, 175 (1962); Amer. J. Math. 89, 817 (1967)
- M. Eichler, D. Zagier, *The Theory of Jacobi Forms*, Birkhäuser (1985)
- V. A. Gritsenko, V. V. Nikulin, Amer. J. Math. 119, 181 (1997) — χ₁₀ als Borcherds-Produkt
- M. Milgrom, MNRAS 403, 886 (2010) — QUMOND (Feldform)
- K. Kuijken, G. Gilmore, ApJ 367, L9 (1991); C. F. McKee, A. Parravano, D. J. Hollenbach, ApJ 814, 13 (2015)
- E. Poggio et al., Nature Astronomy 4, 590 (2020) — Warp-Präzession (Gaia)
- L. S. Sparke, S. Casertano, MNRAS 234, 873 (1988); K. C. Freeman, ApJ 160, 811 (1970)
- J. Bland-Hawthorn, O. Gerhard, ARA&A 54, 529 (2016)
- F. P. Pijpers, MNRAS 297, L76 (1998) — Drehimpuls der Sonne
- Event Horizon Telescope Collaboration, ApJL 875, L6 (2019) — M87*; GRAVITY Collaboration, A&A 657, L12 (2022) — Sgr A*
