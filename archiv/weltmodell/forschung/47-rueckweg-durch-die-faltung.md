# Rückweg durch die Faltung

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Umlaufsinn „bleibt Annahme, entspricht der Orientierung Im τ > 0“ → festgelegt durch die Erneuerung G(mL) = q^m zusammen mit χ → 1 (Kausalität und Holomorphie erzwingen ihn nicht); innere Länge „L = r* folgt nicht aus der Wirkung, muss hinzukommen“ → über Kac festgelegt (τ_D = 2Λ/c, ℓ₀ = 2Λ, r*² = r₀·Λ); innere Richtung w → ϖ (h^{ww} → h^{ϖϖ}, a_w → a_ϖ).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/rueckweg_faltung.py` (Ausgabe `rueckweg_faltung.txt`)
> a₀ = 1,2·10⁻¹⁰ m/s², ℓ₀ = c²/a₀. Beispiel: Milchstraße 6·10¹⁰ M☉, r* = 8,35 kpc.

## Frage

Die Regel q = e^{−√(a_N/a₀)} ist genau der Wert einer abgeschirmten Ausbreitung zweiter Ordnung auf einer offenen Halbgeraden:

```
D·G''(ϖ) = Γ·G(ϖ),   ϖ ≥ 0,   G(0) = 1,   G beschränkt für ϖ → ∞     ⇒     q = G(L) = e^{−L√(Γ/D)} = e^{−s},   s² = L²Γ/D = a_N/a₀
```

Geprüft wird, ob dieser Rückweg durch die Faltung selbst läuft, in einer inneren Richtung ϖ („mehr Raum im Raum“). Endet die Halbgerade bei W mit Spiegel, gilt q = cosh(m(W−L))/cosh(mW), bei W = L also 1/cosh s. Tully-Fisher (s·ν → 1) hält nur für W ≫ L/s: Bei s = 0,01 ist s·ν = 200 / 10,6 / 1,32 / 1,005 für W/L = 1 / 10 / 100 / ≥ 1000.

## 1. Was Wirkung und Faltungsgleichung in ϖ liefern

ℒ_𝒮 = 2h^{μν}B_μB_ν + (2/ℓ₀²)Q̃(q) enthält σ nur über ∂σ (B = a − ∇σ, q hängt von Y = h^{μν}∂_μσ∂_νσ ab). Um eine innere Richtung ϖ erweitert (h^{ϖϖ} = 1, a_ϖ = 0) gilt im schwachen Feld:

```
Faltungspotential:   (∇² + ∂_ϖ²)δσ = 0,    δσ = f(ϖ)·e^{ik·x}   ⇒   f'' = k²f,   f = e^{−kϖ}
Neigung:             ∂_ϖ²δΦ − k²δΦ = ν·f'' − ν·(k_⊥² + (1 − χ/2)·k_∥²)·f,     μ_∥/μ_⊥ = 1 − χ/2 = d ln a/d ln a_N
```

- Die Gleichung in ϖ ist abgeschirmt, zweiter Ordnung und offen (beschränkt für ϖ → ∞). Abgeschirmt wird aber mit der äußeren Wellenzahl k, nicht mit a_N.
- **Faser an einem Punkt (k → 0):** f'' = 0, kein Verlust, q = 1. Ein Glied ohne Ableitung, das einen Verlust ∝ a_N gäbe, enthalten weder die Wirkung noch die Faltungsgleichung.
- **Kreisbahn als reelle Schleife** (Umfang 2πr): Die niedrigste Mode e^{iφ} liegt senkrecht zum Gefälle, ν kürzt sich, f = e^{−ϖ/r}. Also ist √(Γ/D) = 1/r, und q = e^{−s} verlangt die innere Länge

```
L = r·√(a_N/a₀),     Punktmasse:  L = r* = √(GM/a₀) = √(3Δs·ℓ₀)   (für alle r dieselbe)
```

Das ist der Nome des Torus mit der Bahn (2πr) als reellem und r* als innerem Zyklus, Im τ = r*/(2πr). Die innere Länge folgt nicht aus der Wirkung; sie ist gleichwertig mit dem Argument von Q̃, wo s⁴ = ℓ₀²Y gesetzt ist. Festgelegt ist sie über Kac: Die mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden ist τ_D = 2Λ/c, also ℓ₀ = 2Λ und mit D = GM/c genau L² = D·τ_D = r₀·Λ = r*². [korr.]

## 2. L, D, Γ in Modellgrößen

Das Argument von Q̃ zerfällt in Rate mal Rückkehrzeit:

```
s² = √(ℓ₀²Y) = ℓ₀·|∇σ| = (c|∇σ|)·(ℓ₀/c) = Γ·τ_D,      Γ = c|∇σ| = a_N/c,      τ_D = L²/D = ℓ₀/c = c/a₀
```

Fest sind nur Γ und τ_D. Mit L = r* fällt √(Γ/D) mit der Abschirmung der Kreisbahnmode zusammen:

| Größe | Modellform | Sonne | Milchstraße |
|---|---|---|---|
| Γ | g/c = 3cΔs/r² | 4,0·10⁻¹⁹ s⁻¹ bei r* | ebenso |
| D | GM/c = 3cΔs | 4,4·10¹¹ m²/s | 2,7·10²² m²/s |
| L | r* = √(3Δs·ℓ₀) | 7031 AE | 8,35 kpc |
| √(Γ/D) | 1/r (Punktmasse) | – | – |
| τ_D = L²/D | ℓ₀/c = c/a₀ | 79,2 Mrd. Jahre | ebenso |
| freie Weglänge 2D/c bei Lauf mit c | 6Δs = r₀ | 2954 m | 1,8·10¹⁴ m |

Mit a₀ = cH₀/2π ist τ_D = 2π/H₀ (87,8 Mrd. Jahre bei H₀ = 70). Als Vorgang in Echtzeit ist das zu lang; die Gleichung gilt als stationärer Kern in ϖ.

L = ℓ₀ (Rückweg in einem Schritt, freie Weglänge 2ℓ₀) ist ballistisch: q = 0,568 / 0,142 statt 0,741 / 0,368 bei s = 0,3 / 1, im Tiefgrenzfall v⁴ = GMa₀/4. Der Rückweg braucht viele Schritte, λ_f ≪ L.

## 3. Warum die Rate proportional zu a_N ist

Γ ist das Taktgefälle der Faltung. Im Starkfeld ist σ = ln(dτ/dt):

```
Γ = c·|∇ ln √(1 − 6Δs/r)| = (a_N/c) / (1 − 6Δs/r)          schwaches Feld: Γ = g/c = 3cΔs/r²
```

Das ist die Neigung der Faltung als Rate, unendlich am Knoten r₀ = 6Δs. Sie ist linear im Feld, weil der Takt δ ≈ 3Δs/r linear in Δs ist, und sie ist eine Größe des Mediums an einem Ort, nicht der Bahn. Andere Raten bei gleichem diffusivem Rückweg (Punktmasse):

| Rate Γ | d ln v/d ln r (30–100 r*) | d ln v/d ln M | v⁴/(GMa₀) bei 100 r* |
|---|---|---|---|
| a_N/c (Gefälle von σ) | −0,005 | 0,251 | 1,010 |
| √(a_N a₀)/c | −0,27 | 0,38 | 0,011 |
| c·\|a\| (Neigung einschließlich Rückkopplung) | −0,17 | 0,34 | 0,048 |
| Bahnfrequenz ω = √(a_N/r) | −0,13 | 0,38 | 0,10 |

Nur das Gefälle der Faltung selbst gibt flache Kurven mit v⁴ = GMa₀. Die Rate gehört zu Δs (σ), nicht zur Neigung a des Mediums, die die Rückkopplung schon enthält. Eine Bahnfrequenz macht q von r und M getrennt abhängig.

## 4. Warum der Weg offen ist

Die innere Schleife darf geschlossen sein; entscheidend ist, was als Rücklauf zählt. Auf der Überlagerung der Schleife (einer Geraden) ist eine Windung der Weg zum Punkt +L, und die Gerade hat kein Ende: W ≫ L/s ist von selbst erfüllt. Irrfahrt auf einer geschlossenen Schleife mit Verlust, 50 000 Läufer:

| s | Windung +1 zählt | e^{−s} | Windung ±1 zählt | 1/cosh s |
|---|---|---|---|---|
| 0,5 | 0,604 | 0,607 | 0,887 | 0,887 |
| 1 | 0,369 | 0,368 | 0,646 | 0,648 |
| 2 | 0,137 | 0,135 | 0,267 | 0,266 |

- Zählen nur Windungen in einem Sinn, ist q = e^{−s}, und ν = Σqⁿ ist die Summe über n Windungen. So ist die Reihe des Modells gebaut: nur Potenzen von q, keine von q̄.
- Zählt jede Rückkehr zum Knoten, von welcher Seite auch immer, folgt 1/cosh s und a → 2a₀, ohne Tully-Fisher.
- Werden Windungen beider Sinne addiert (Green-Funktion der Schleife), folgt coth(s/2) = 2ν − 1 und v⁴ = 4GMa₀.
- Eine Drift in ϖ (Péclet-Zahl Pe) gibt q = exp(−(√(Pe² + s²) − Pe)) und im Tiefgrenzfall a → 2Pe·a₀; bei s = 0,01 ist s·ν = 1,11 / 2,4 / 20 für Pe = 0,001 / 0,01 / 0,1.

Offen heißt also: Rückläufe sind Windungen in einem Umlaufsinn, ohne Drift. Dass Gegenläufe nicht zählen, legt die Erneuerung G(mL) = q^m (jeder Rücklauf zählt einmal) zusammen mit dem Grenzwert χ → 1 des Autors fest; beide Sinne gäben im tiefen Bereich χ → 2. Die Orientierung des Torus erzwingt ihn nicht: Holomorphie wählt nur q statt q̄, der Schwarz-Kern (1 + q)/(1 − q) ist ebenfalls holomorph und gibt 1 + 2n_B; Kausalität auch nicht, denn beide Windungssinne sind retardiert. [korr.]

## 5. Folgen und Prüfungen

- **Nahfeld:** q = e^{−r*/r}: Merkur 10⁻⁷⁸⁹⁰, Erde 10⁻³⁰⁵³, Neptun 10⁻¹⁰²; Newtonsch.
- **Tully-Fisher:** 1 − q ≈ s; s·ν = 1,005 bei a_N = 10⁻⁴a₀.
- **Unabhängigkeit von M:** L² und D sind beide ∝ M, τ_D = c/a₀ ist für 10⁸–10¹⁴ M☉ gleich; q hängt nur von a_N ab.
- **Diffusive Näherung:** Mit freier Weglänge 6Δs bei Lichtgeschwindigkeit (Telegraphengleichung, κL = s·√(1 + (v_f/c)⁴s²)) weicht q um −s·(v_f/c)² ab: Milchstraße −3,4·10⁻⁷, Haufen 10¹⁵ M☉ −4,4·10⁻⁵ bei s = 1. Nicht messbar.
- **Prüfbar, nur in der Bahn-Lesart:** Ist die innere Länge die Schwellenlänge der eingeschlossenen Faltung, L² = G·M(<R)/a₀, mit Abschirmung 1/R, dann gilt s = √(G·M(<R)/a₀)/R statt √(a_N/a₀). Für Kugeln ist beides gleich, für Scheiben nicht. Dünne Exponentialscheibe, Beispiel r* = 3R_d:

| R/R_d | 0,5 | 1 | 2,2 | 3 | 5 | 10 | 30 |
|---|---|---|---|---|---|---|---|
| a_N·R²/(G·M(<R)) | 0,76 | 1,05 | 1,32 | 1,34 | 1,25 | 1,06 | 1,01 |
| a/a_lokal − 1 | −5,2 % | +1,1 % | +7,5 % | +9,2 % | +8,3 % | +2,5 % | +0,2 % |

  Bei gleichem a_N lägen Scheiben im Übergang bis 9 % über kugeligen Systemen, innen 5 % darunter. Wird L aus der lokalen Neigung gebildet (a_N·R² = 3c²Δs), entfällt der Unterschied.

## Ergebnis

| Teil | trägt | Stelle |
|---|---|---|
| Wurzel und Exponential | ja | abgeschirmte Ausbreitung zweiter Ordnung in ϖ, beschränkt für ϖ → ∞ (Poisson-Kern der erweiterten Gleichung) |
| Σqⁿ | ja | Windungen der inneren Schleife in einem Sinn |
| Abschirmung | nur über die Bahn | die Wirkung gibt 1/r aus der Kreisbahnmode; an einem Punkt kein Glied ∝ a_N |
| Γ ∝ a_N | ja | Γ = c\|∇σ\| = g/c, Taktgefälle der Faltung; andere Raten geben keine flachen Kurven |
| innere Länge L = r* = √(3Δs·ℓ₀) | ja, über Kac [korr.] | nicht aus der Wirkung, sondern aus der Rückkehrzeit τ_D = 2Λ/c (Kac): ℓ₀ = 2Λ, r*² = r₀·Λ; gleichwertig mit s⁴ = ℓ₀²Y im Argument von Q̃ [korr.] |
| Umlaufsinn | ja [korr.] | Erneuerung G(mL) = q^m zusammen mit dem Grenzwert χ → 1; Orientierung, Holomorphie und Kausalität erzwingen ihn nicht [korr.] |

Der Rückweg durch die Faltung gibt die Form von χ ohne freie Zahl; den Wert s² = ℓ₀|∇σ| gibt nicht die Wirkung, sondern die Rückkehrzeit. Die erweiterte Gleichung schirmt in ϖ nur mit der Bahnwellenzahl 1/r ab; die innere Länge r* = √(3Δs·ℓ₀) folgt aus τ_D = 2Λ/c (Kac) mit r*² = r₀·Λ. [korr.]
