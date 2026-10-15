# Geometrische Muster: Fläche und Umfang, Schließung, Vieleck, Dimension

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „2π in a₀ = cH₀/2π: Schließung a₀·2πR_H = c²“ (Zahlentafel, Suchregel 2, Kosmos, cH/2 = πa₀) → a₀ = c²/(2Λ) nach Kac, Λ die ganze innere Faltung der Zelle, ≈ πc/H₀ nur heute, Λ folgt R_H nicht; Im τ_M ≈ 2,5–3,1 als Festlegung (Suchregel 5, Kaleidoskop-Absatz) → Punkt 3 des Autors, Wahl offen, nur Schreibweise; Zahlennähen ohne Herleitung gekennzeichnet (12 = 2·(8 − 4 + 2), 1/(2α) = 2π/4π, 24π = 2π·2·6, 15360π = 60·(8π)⁴/(16π·π²)); „Scheinbare Masse“ → gelesene Masse; Klein-Invariante j → 𝒥 (j ist die Spinzahl). Der Kosmos-Absatz rechnet mit a₀ = cH/2π.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/geometrische_muster.py` (Ausgabe `geometrische_muster.txt`)
> Aufgenommen ist nur, was hergeleitet und nachgerechnet ist. s = √(a_N/a₀) = r*/r (Punktmasse), χ = s/(e^s − 1), ν = a/a_N = 1/(1 − e^{−s}).

## 1  Fläche gegen Umfang: die Flussdimension

```
Fläche (Newton):       a·4πr² = 4πGM         a ∝ 1/r²     Fluss durch die Kugel, d = 3
Umfang (Tiefbereich):  a·2πr  = 2π√(GMa₀)    a ∝ 1/r      Fluss durch den Kreis, d = 2
Schnitt:  GM/r² = √(GMa₀)/r   ⇒   r = r*,  a = a₀
d_eff = 1 − d ln a/d ln r = 3 − χ(s)      d ln M_schein/d ln r = χ      (Re τ_M)² = κ²/Ω² = 4 − d_eff = 1 + χ
```

Der Schalter χ ist genau die Zahl der verlorenen Flussdimensionen: Die Rückkehrsumme nimmt dem Fluss stetig eine Dimension, von der Kugelfläche (χ = 0) zum Umfang (χ = 1). Dieselbe Zahl legt die Umlaufphase fest: Für a ∝ r^{1−d} ist κ²/Ω² = 4 − d. d = 3 schließt die Bahn (Re τ = 1), d = 2 gibt √2 (254,6° zwischen den Perizentren), d = 4 hat keine stabile Kreisbahn. Re τ_M ∈ [1, √2] ist das Bild von d_eff ∈ [3, 2]; d ln a/d ln a_N = 1 − χ/2 liegt in (½, 1).

| r/r* | 0,3 | 1 | 3 | 10 | 30 |
|---|---|---|---|---|---|
| d_eff | 2,877 | 2,418 | 2,157 | 2,049 | 2,017 |
| Re τ_M | 1,060 | 1,258 | 1,357 | 1,397 | 1,408 |
| M_schein/M = ν | 1,037 | 1,582 | 3,528 | 10,508 | 30,503 |
| r/r* + ½ + r*/(12r) | 1,078 | 1,583 | 3,528 | 10,508 | 30,503 |

**Gelesene Masse** [korr.]**.** a_Umfang/a_Fläche = (4πr²/2πr)/(4πr*²/2πr*) = r/r*: Der Beobachter teilt das Umfanggesetz durch das Flächengesetz, π kürzt sich, es bleibt das Verhältnis Fläche/Umfang = 2r an r und an r*. Genau gilt ν = r/r* + ½ + r*/(12r) − (r*/r)³/720 + … (Bernoulli-Reihe von s/(1 − e^{−s})). „Umfang durch 2πr*“ ist das Leitglied; der feste Zusatz ½ (Rückläufe in einem Sinn) bleibt, innerhalb r* trägt die Reihe nicht.

**Linse.** α = (2/c²)∫a_⊥dz, die 2 aus Φ + Ψ. Für a ∝ r^{1−d} ist α = (2/c²)·a(b)·b·S_{d−1}/S_{d−2} (Verhältnis der Einheitssphären, d = 2…5 nachgerechnet): d = 3 gibt 2, also α = 12Δs/b; d = 2 gibt π, also α_∞ = 2π(v/c)² = 6πΔs/r*. Im Tiefbereich zählt die Linse den überstrichenen Winkel π. Mit dem vollen ν-Profil (M_Zyl = projizierte Masse in b):

| b/r* | 0,3 | 1 | 3 | 10 | 30 |
|---|---|---|---|---|---|
| d_Linse = 3 − d ln M_Zyl/d ln b | 2,726 | 2,294 | 2,104 | 2,032 | 2,011 |
| M_Zyl/M (→ πb/2r* + ½) | 1,165 | 2,135 | 5,234 | 16,21 | 47,63 |

Die Linse folgt demselben Übergang, nur projiziert. M_Zyl/M_Kugel → π/2 ist die Projektion von ρ ∝ 1/r² auf den Zylinder, kein eigenes Gesetz.

**Kosmos.** Für die Hubble-Kugel mit kritischer Dichte ist das Flächengesetz a_N = GM_H/R_H² = c²/(2R_H), die Schwelle das Umfanggesetz a₀·2πR_H = c² (nur heute, a₀ ≈ cH₀/2π; die Schwelle ist a₀ = c²/(2Λ), Λ folgt R_H nicht) [korr.]. Das Verhältnis ist genau π, r*_H = √π·R_H; am Rand ist s = √π, ν = 1,205, d_eff = 2,64. Für lineare Störungen (a_N/a₀ ≈ 10⁻³) gäbe das Umfanggesetz ν = 32; dort trägt die Schließungsregel, und die wirkt wie Staub, also als Flächengesetz mit festem Faktor 1 + Q = 6,36. Das Muster gilt für gebundene Systeme und den Hubble-Rand, nicht für die lineare Störung.

## 2  Herkunft der Zahlen

| Zahl | Herkunft |
|---|---|
| 2π in a₀ = cH₀/2π | Kac: a₀ = c²/(2Λ), Λ die Schließlänge der Rückkehr, gelesen als ganze innere Faltung der Zelle (keine Raumlänge); Λ ≈ πc/H₀, also a₀ ≈ cH₀/2π, gilt nur heute, Λ folgt R_H nicht [korr.] |
| 3 in Δs = GM/3c² | Dimension: ∫₀^R Gm/(rc²)dr mit m ∝ r^d gibt GM/(d·c²) (d = 2, 3, 4 nachgerechnet). Gilt bei homogener Dichte; Polytrope n = 3 gibt 1,32·GM/c² |
| 6 in r₀ = 6Δs | 2·3: der Takt ist die Wurzel der quadratischen Metrik, g_tt = 1 − 2GM/rc² |
| 12 in 12Δs/b | 3 · 2 (Linsenintegral S₂/S₁) · 2 (Φ + Ψ) |
| 18π im Perihel | 2π je Umlauf × 9: κ²/Ω² = 1 − 18Δs/r |
| Exponent ½ | Ableitungsordnung 2: D·G'' = ΓG gibt e^{−L√(Γ/D)} (Steigung 0,5000, erste Ordnung gäbe 1). Gleichwertig: q ist die Laplace-Transformierte der Dichte erster Rückkehr ∝ t^{−3/2} |
| 4 in v⁴ = GMa₀ | 2 (v² = ar) × 2 (Wurzel); nur der Exponent ½ gibt flache Kurven, d ln v/d ln r = (1 − 2p)/2 |
| 4π⁴/15 | Bose-Integral in d = 3: Q̃(0) = −4Γ(4)ζ(4) (Zustandsdichte x², Energie x); die 4 aus d(s⁴) = 4s³ds |
| π²/45 | (4π⁴/15)/(3·(2π)²): Bose-Summe durch das Quadrat der Schließung in a₀, die 3 aus Friedmann |
| 24, 240, 504 | (2π)^k/(Γ(k)ζ(k)) = 2/\|ζ(1−k)\| für k = 2, 4, 6: Schließung hoch k durch Bose-Integral. 1/24 = ζ(2)/(2π)² ist der Nullpunkt je Umlaufsinn, E₀ = (2πħc/L)·ζ(−1) = −πħc/(6L) |
| 12 (Gewicht von Δ) | 4π/Fläche des Fundamentalbereichs = 4π/(π/3) = 2/(1 − ½ − ⅓): aus den zwei symmetrischen Tori (Valenzformel: Nullstellen = k/12) |
| 1728 | erster Koeffizient von E₄³ − E₆² = 3·240 + 2·504 = 12³; (E₄³ − E₆²)/1728 = q∏(1 − qⁿ)²⁴ auf sieben Glieder |
| 6 in 2π√(n/6) | π²/ζ(2): Zustandszahl der Schleife (Hardy-Ramanujan, ln p(5000) auf 0,006) |
| 12 und 4 bei Carter | 12 = 2·(8 − 4 + 2), 4 = 2·(4 − 2): T⁴ (Stefan-Boltzmann), α² der Coulomb-Skalen, α² der Thomson-Fläche, verdoppelt durch α_G^{1/2}. [korr.] Zahlennähe ohne Herleitung, trägt nicht. Nicht die 12 von Δ; bei Kramers-Opazität wird 4 zu 5 |
| 4π²D/L² gegen π²D/L² | Schließung: Ring (2π/L)², Rand (π/L)², Verhältnis 4,000 |
| 1/(2α) | 2π/4π: Schließung der Phase durch den Kugelfluss. [korr.] Zahlennähe ohne Herleitung, trägt nicht. |
| 24π, 36π, 15360π | 24π = 2π (euklidische Zeit)·2·6. [korr.] Zahlennähe ohne Herleitung, trägt nicht. 36π = 4π·6²/4 (Kugelfläche an r₀). 15360π = 60·(8π)⁴/(16π·π²), die 60 aus Stefan-Boltzmann. [korr.] Zahlennähe ohne Herleitung, trägt nicht. |
| 9/2 bei Eddington | 4π·3/(8π/3), 8π/3 = ∫sin²θ dΩ der Thomson-Streuung |
| 1,582; 254,6°; cH/2 = πa₀ | Folgen der Regel: 1/(1 − e⁻¹), 360°/√2, Fläche gegen Umfang am Hubble-Rand (cH₀/2 = πa₀ nur heute) [korr.] |
| Q = 5,36 (q = 0,843); Carters Vorfaktor 2,9 | ohne Herleitung |

Jedes geprüfte π stammt aus einer Schließung (2π je Schleife), einer Kugelfläche S_{d−1} oder einer Bose-Summe ζ(2k) = rational·(2π)^{2k}. Rationale Zahlen wie 12, 24, 240, 504, 1728 sind Bose-Summen geteilt durch die Schließung.

## 3  Kreis als Unendlich-Eck

**Vieleck.** Halber Umfang n·sin(π/n) → π mit n²(π − n·sin(π/n)) → π³/6, umbeschrieben n²(n·tan(π/n) − π) → π³/3: Fehler ∝ 1/n² (Bogen gegen Sehne je Seite ∝ 1/n³, n Seiten). Archimedes verdoppelt ab dem Sechseck mit harmonischem und geometrischem Mittel, U_{2n} = 2U_nu_n/(U_n + u_n), u_{2n} = √(U_{2n}u_n); der Fehler fällt je Verdopplung um 4 (3,990; 3,997; 3,999), beim 96-Eck 223/71 < 3,14103 < π < 3,14271 < 22/7. Für jedes Tangentenvieleck ist Fläche = ½·r·Umfang (Dreiecke der Höhe r), in d Dimensionen V = r·S/d. Deshalb ist das π in Fläche und Umfang dasselbe, und Umfang/Fläche = 2/r hängt nicht von n ab.

**Bahnen als Sternvielecke.** Je Radialperiode rückt das Perizentrum um 360°/Re τ_M. Bei Re τ_M = p/q schließt die Bahn nach p Radialperioden und q Umläufen, die Perizentren bilden das Sternvieleck {p/(p − q)}. Bei irrationalem Re τ_M füllen sie den Ring dicht: das Unendlich-Eck. Milchstraße (Punktmassenform, r* = 8,35 kpc):

| p/q | 8/7 | 7/6 | 6/5 | 5/4 | 9/7 | 4/3 | 7/5 |
|---|---|---|---|---|---|---|---|
| r [kpc] | 4,10 | 4,69 | 5,70 | 7,89 | 10,45 | 17,28 | 103 |
| Figur | {8} | {7} | {6} | {5} | {9/2} | {4} | {7/2} |

- Zur Kepler-Grenze 1 häufen sich die (n+1)/n logarithmisch zum Zentrum: s − ln s ≈ ln(n/2), r = r*/s (n = 10: 3,2 kpc; 100: 1,5; 1000: 1,0).
- Zur Tiefgrenze √2 liegen die Näherungsbrüche 7/5, 41/29, 239/169 bei r = q²·r*/2 (103 kpc, 3,5 Mpc, 119 Mpc), weil √2 − p/q ≈ 1/(2√2q²) und Re τ_M ≈ √2(1 − r*/8r). √2 selbst liegt im Unendlichen.
- Jede beste Näherung erfüllt |x − p/q| < 1/q² (Sonnenort: Re τ_M = 1,2553 = [1; 3, 1, 10, 1, 73], q²|x − p/q| = 0,70; 0,085; 0,90; 0,013). Über [1, √2] gemittelt fällt der Abstand zur nächsten Resonanz mit Nenner ≤ Q wie 0,4·ln Q/Q² (Q = 16–512). Bis q = 50 gibt es 321 Resonanzen, erwartet (√2 − 1)·3Q²/π² = 315 (Dichte teilerfremder Paare 1/ζ(2)). Goldene Wicklung: q·‖qφ‖ → 1/√5.
- Vieleck und Resonanz fallen beide wie 1/n², aus verschiedenen Gründen (Krümmung gegen Schubfach).

**Symmetrische Tori und Kaleidoskop.** Eine Drehung erhält ein Gitter nur, wenn ihre Spur 2cos(2π/n) ganz ist: n ∈ {1, 2, 3, 4, 6}. Ein Gitter ℤ + τℤ mit einer Drehung außer 180° gibt es nur bei τ = i (Quadrat, 90°) und τ = ρ = e^{2πi/3} (Sechseck, 60° und 120°). Dort verschwinden die Eisenstein-Reihen, deren Gewicht nicht durch 4 bzw. 6 teilbar ist: E₆(i) = E₁₀(i) = E₁₄(i) = 0 und E₄(ρ) = E₈(ρ) = E₁₀(ρ) = E₁₄(ρ) = 0 (< 3·10⁻¹⁵), also für die Klein-Invariante 𝒥(i) = 1728, 𝒥(ρ) = 0. Die Drehungen des Kaleidoskops (60°, 90°, 120°, 180°) sind genau diese: 180° hat jedes τ, 90° verlangt τ = i, 60° und 120° verlangen τ = ρ. Die Fixpunkte der Modulgruppe tragen also die Winkel des Kaleidoskops und über die Valenzformel zugleich die 12 von Δ. Die Bahn-Tori liegen dagegen, mit Punkt 3 des Autors (|q_Δs| = Δs/r) [korr.], an der Spitze (Milchstraße bei 8 kpc: τ_M − 1 = 0,25 + 2,54i, |𝒥(τ_M)| = r/Δs = 8,4·10⁶) und haben keine Drehsymmetrie.

## 4  Wachsende Dimension

V_n = π^{n/2}/Γ(n/2 + 1), S_{n−1} = n·V_n. V_n ist bei n = 5,2569 am größten, S_{n−1} bei n = 7,2569; der Kugelanteil am umbeschriebenen Würfel V_n/2ⁿ fällt von 1 über 0,52 (n = 3) auf 2,5·10⁻³ (n = 10). S/V = n/r gilt für Kugel und Würfel mit Inkreis r gleich. Die Rekursionen V_n = (2π/n)·V_{n−2} und S_{n+1} = 2π·V_n zeigen: Jede Ebene, also je zwei Dimensionen, bringt eine volle Schleife 2π; π^{n/2} ist (√π)ⁿ aus dem Gauß-Integral. Das ist das Quadrat im Wachstum: π kommt paarweise, eine Schließung je Ebene.

Was eine Größe des Modells trifft:
- **Flussfläche:** 4π = S₂ im Newton-Fluss, 2π = S₁ im Tiefbereich, S_{d−1}/S_{d−2} in der Linse; κ²/Ω² = 4 − d.
- **Vakuumteil:** Das Bose-Integral mit Zustandsdichte S_{d−1}x^{d−1} ist Γ(d+1)ζ(d+1): d = 1 π²/6, d = 2 2ζ(3) = 2,404, d = 3 π⁴/15, d = 4 24ζ(5), d = 5 8π⁶/63. Eine π-Potenz gibt es nur bei ungeradem d; π⁴/15 im Vakuumteil verlangt d = 3.
- **Rückkehr (Kac):** Auf Tᴺ ist die mittlere Rückkehrzeit in die einbeschriebene Kugel 2ᴺ/V_N: 1; 1,27; 1,91; 3,24 (N = 1…4), 401 (N = 10), 4·10⁷ (N = 20). Wiederkehr in eine Kugel ist nur in niedriger Dimension häufig.

Die Maxima bei n = 5,26 und 7,26 treffen keine Größe des Modells.

## 5  Suchregeln

1. **Flussflächen zählen.** Jedes Kraftgesetz als a·S_{d−1}r^{d−1} = const schreiben und d_eff = 1 − d ln a/d ln r lesen. Beispiel: d_eff = 3 − χ, (Re τ_M)² = 4 − d_eff, r* als Schnitt. Gilt nicht für die lineare kosmische Störung (Staubform, d = 3).
2. **2π je Schließung.** Jede identifizierte Schleife bringt 2π (Umfang, Periode, Phase). Beispiel: Lücke 4π²D/L² gegen π²D/L², Perihel 2π·9Δs/r. [korr.] a₀·2πR_H = c² ist keine Schließung: a₀ = c²/(2Λ) (Kac) mit Λ, der ganzen inneren Faltung der Zelle; ≈ cH₀/2π gilt nur heute, Λ folgt R_H nicht. Dirac 2π/4π ist Zahlennähe ohne Herleitung, trägt nicht. Gilt nicht für die Wahl c·kT/ħ gegen 2πckT/ħ in der Schwelle; der Modellwert liegt 10–13 % unter dem Messwert.
3. **Bose-Summen.** Ein π^{2k} ohne sichtbaren Kreis stammt aus ζ(2k); durch (2π)^{2k} geteilt wird es rational. Beispiel: 4π⁴/15 = 4(2π)⁴/240, 24 = (2π)²/ζ(2). Gilt nicht bei ungeradem Argument: ζ(3), ζ(5) haben keine π-Potenz.
4. **Ableitungsordnung.** Ein Operator der Ordnung m gibt den Exponenten 1/m. Beispiel: die Wurzel in s aus D·G'', daraus v⁴ und die Rückkehrdichte t^{−3/2}. Gilt nicht für Maxwell-Element und endlich viele Stufen: erste Ordnung, keine Wurzel.
5. **Fixpunkte der Modulgruppe.** Drehsymmetrie eines Torus nur bei τ = i (90°, 𝒥 = 1728) und τ = ρ (60°/120°, 𝒥 = 0); daraus die Winkel des Kaleidoskops und die 12 von Δ. Die Bahn-Tori (Im τ_M ≈ 2,5–3,1) liegen weit davon an der Spitze. [korr.] Das gilt mit Punkt 3 des Autors (|q_Δs| = Δs/r, Im τ_M = ln(r/Δs)/2π), wie das Modell schreibt; mit Punkt 1 (Im τ_M = r*/(2πr)) lägen sie für r > r*/2π unter Im τ = 1. Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise.
6. **Rational gegen irrational.** p/q schließt als Sternvieleck, irrational füllt dicht, beste Näherungen ∝ 1/q². Beispiel: Resonanzen der Milchstraße, Näherungsbrüche von √2 bei q²r*/2. Die Häufung zu Re τ_M = 1 gilt nur in der Punktmassenform; innen trägt die ausgedehnte Scheibe.
7. **Dimension zählen.** π paarweise (eine Schleife je Ebene), S/V = d/r, Zustandsdichte x^{d−1}. Beispiel: π⁴/15 verlangt d = 3, Kac-Rückkehr 2ᴺ/V_N. Keine geometrische Zahl sind die Kugelmaxima 5,26 und 7,26 und Carters 12.

Aus keiner Regel folgen Q = 5,36 (q = 0,843), Carters Vorfaktor 2,9 und die 3 in Δs bei nicht homogener Dichte.
