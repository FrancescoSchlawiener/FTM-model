# Schwelle als Rückkehrzeit im geschlossenen Raum

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Deckel (b) atmende Zelle (a₀ ∝ H) und (c) feste mitbewegte Zelle (Λ ∝ a) als nicht tragend gekennzeichnet (Vakuumteil zu jeder Zeit 22 % bzw. wie Krümmung 0,22); Abschnitt 3 als Rechnung unter einem Zeitgesetz a₀(z) gekennzeichnet; „nur die atmende Zelle zulässig“ → keiner der gerechneten Deckel; „jede Länge, bei der der Raum schließt, wächst; a₀ fällt monoton“ → a₀ früh kleiner (0,76·a₀ bei der Rekombination, 0 bei der Nukleosynthese), seit z ≈ 10 fest; feste Zelle U₀ = 24,3 bzw. 10,8 Gpc < 27 Gpc → Zelle aus a₀ = 3c⁴/(2Gρ_mL³), L = 26,9–28,1 Gpc; Ergebniszeilen zu a₀(z) auf die Annahme beschränkt, fehlende Umkehr bei z ≈ 0,4 → Form der Dunklen Energie in der Rückstellung. Umbenennungen: w (Ebenen der Halbgeraden) → ϖ; U (atmende Zelle, c²/a₀) → ℓ₀; U₀ (feste Raumzelle) → L₀; Ω_Λr → Ω_DE,r. Die Notiz rechnet in K1b und K2 mit a₀ = cH/2π (heute 1,04·10⁻¹⁰).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/schwelle_rueckkehrzeit.py` (Ausgabe `schwelle_rueckkehrzeit.txt`)
> Eingaben: H₀ = 67,4 km/s/Mpc, Ω_m = Ω_b(1+Q) = 0,3134 (Q = 5,364), Ω_r = 9,14·10⁻⁵, k_s = π²/45. Werte mit [G] aus dem Gedächtnis (σ₈ = 0,811, n_s = 0,965, CPL-Werte DESI DR2).

## Frage

a₀ = c²/ℓ₀ = c/τ_D. Ist τ_D die mittlere Rückkehrzeit nach Kac, ⟨T_A⟩ = t·μ(X)/μ(A), und gibt das ein a₀(z) mit Maximum bei z ≈ 0,4? Jede Lesart selbstkonsistent mit E² = (R + Ω_DE,r)/(1 − k_sΘ²) [korr.], Θ = 2πa₀/(cH), w_s = −1 − (2/3)·d ln a₀/d ln a.

## 1 Kac auf der Halbgeraden der Hierarchie

Der Rückweg ist eine symmetrische Kette auf den Ebenen ϖ [korr.]: Schritte r₀ = 6Δs mit c (D = r₀c/2), der Ort unten, der Deckel bei Λ reflektierend.

```
⟨T_Ort⟩ = Σ_j deg(j)·h_j/c / deg(0) = 2Λ/c          für jede Schrittfolge h_j
τ_D = ⟨T_Ort⟩   ⇔   ℓ₀ = 2Λ   ⇔   r*² = 3Δs·ℓ₀ = r₀·Λ
```

- Exakt aus den mittleren Trefferzeiten (N = 10 … 1000 Ebenen, gleiche und zufällige Schritte, also auch springendes D): ⟨T⟩·c/2Λ = 1 auf 12 Stellen.
- Die mittlere Rückkehrzeit einer geschlossenen Kette hängt nur von ihrer Länge ab, nicht vom Schritt r₀ ∝ M. Das ist der Grund, warum τ_D für jede Masse gleich ist. Die Tiefe, die in dieser Zeit erreicht wird, ist √(D·⟨T⟩) = √(r₀Λ) = r*: Der Schwellenradius ist das geometrische Mittel aus Grenzradius und Deckel.
- Zahlen: Mit Λ = c²/2a₀ = 12,14 Gpc folgt r* = 7031 AE (Sonne, r₀ = 2954 m) und 8,35 kpc (Milchstraße, r₀ = 1185 AE); mit Λ = πc/H₀ = 13,97 Gpc 7544 AE und 8,96 kpc.
- Vakuumteil: Θ = π·(c/H)/Λ. Ω_s < 1 verlangt Λ > π√k_s·c/H = 1,47·c/H, Ω_s < Ω_DE heute Λ > 1,78·c/H₀.

## 2 Lesarten

| Lesart | X | A | τ_D | a₀ |
|---|---|---|---|---|
| K1 Ort | Halbgerade [0, Λ] | unterste Ebene | 2Λ/c | c²/2Λ |
| K2 gebunden | Orte × Ebenen (Faltung × Länge in ϖ) [korr.] | Ort der gebundenen Orte, Anteil F | 2Λ/(cF) | (c²/2Λ)·F |
| K3 Volumen | Volumen V, Bewegung isotrop mit c | Kugel r* bzw. gebundene Bereiche | 4V/(c·S_A) | c/⟨T⟩ |

Deckel in K1: (a) eigener Schwellenradius der Hubble-Kugel, Λ² = GM_H/a₀ mit ℓ₀ = 2Λ ⇒ Λ = 2GM_H/c² = c/H, a₀ = cH/2, Θ = π, Ω_s = π⁴/45 = 2,17 > 1: kein H(z). Ebenso jeder Horizont, der in de Sitter gegen c/H geht. Mit Λ = √π·c/H festgehalten: a₀ = 1,85·10⁻¹⁰, Ω_s = 0,689 > Ω_DE. (b) Atmende Zelle ℓ₀ = 2πc/H, Deckel am Gegenpunkt Λ = ℓ₀/2 (gleichwertig: Schleife vom Umfang ℓ₀, ⟨T⟩ = ℓ₀/c): a₀ = c²/ℓ₀ = cH/2π. (c) Feste mitbewegte Zelle L = L₀a, Λ = L₀a/2: a₀ = c²/(L₀a). [korr.] (b) und (c) tragen nicht: Mit Λ = πc/H bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); mit den Schleifen atmend (Λ ∝ a, Fall c) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). c²/a₀ = ℓ₀ = 2Λ ist eine Länge der inneren Richtung, keine Raumlänge. Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute; die Zelle folgt aus a₀ = 3c⁴/(2Gρ_mL³).

F in K2 ist der Anteil der Faltung in Bereichen jenseits der Umkehr (δ = 1,062) oder des Kollapses (1,686), F = erfc(δ/(√2·σ(M_min)·D)), lineares Feld mit Ω_m = Ω_b(1+Q). Die kleinste gebundene Masse M_min legt das Modell nicht fest; gerechnet wird die Schar. Mit Λ = ℓ₀/2 = πc/H ist a₀ = (cH/2π)·F und Θ = F [korr.].

## 3 Selbstkonsistent

[korr.] Diese Rechnung setzt ein Zeitgesetz a₀(z) nach dem Deckel voraus (a₀ ∝ H, ∝ 1/a oder ∝ H·F); das trägt nicht (Deckel b und c oben). Im Modell ist a₀ seit z ≈ 10 fest; die Form der Dunklen Energie wird in der Rückstellung gesucht, als Bindungsrate dF/d ln a.

| Lesart | a₀(0) [10⁻¹⁰] | a₀(1)/a₀(0) | z_max | max/heute | a₀(1)/max | a₀(2)/a₀(0) | d ln a₀/d ln a | w_s(0) / w_s(1) | Ω_s(1090) |
|---|---|---|---|---|---|---|---|---|---|
| K1b atmende Zelle | 1,042 | 1,953 | – | – | – | – | −0,60 | −0,60 / −0,16 | 0,22 |
| K1c feste Zelle, L₀ = 27 Gpc [korr.] | 1,079 | 2,000 | – | – | – | – | −1,00 | −0,33 / −0,33 | 5·10⁻⁴ |
| K1c, L₀ = c²/a₀ = 24,3 Gpc [korr.] | 1,200 | 2,000 | – | – | – | – | −1,00 | −0,33 / −0,33 | 6·10⁻⁴ |
| K2, M_min = 10^12,77 M☉/h (Umkehr) | 0,523 | 0,957 | 0,40 | 1,009 | 0,949 | 0,649 | −0,02 | −0,98 / −1,26 | ≈ 0 |
| K2c, feste Zelle 27 Gpc, 10^13,09 | 0,482 | 0,943 | 0,40 | 1,089 | 0,865 | 0,460 | −0,47 | −0,69 / −1,60 | ≈ 0 |
| verlangt [G] | 1,2 | – | 0,35–0,44 | 1,11–1,26 | 0,87 | 0,39–0,51 | −0,76 bis −1,56 | w₀ > −1, w_a < 0 | klein |

- **K2-Schar:** kleines M_min gibt ein frühes Maximum (z = 8,5 bei 10⁶, 3,3 bei 10¹⁰ M☉/h), großes ein noch wachsendes a₀ (ab 10¹³). Im Fenster 0,35–0,44 liegt das Maximum bei log M_min = 12,74–12,79 (Umkehr) bzw. 11,34–11,41 (Kollaps), beide mit F(0) = 0,50: Die Lage legt nur F(0) fest, nicht die Konvention. Gesamte Dunkle Energie dort w₀/w_a = −0,98/−0,08.
- **K3:** (a) Kugel r* im Hubble-Volumen: ⟨T⟩ = V/(πr*²c) = 3·10³² (Sonne), 5·10²¹ (Milchstraße), 2·10¹⁸ Jahre (Haufen) gegen τ_D = 7,9·10¹⁰; ⟨T⟩ ∝ 1/M gibt a₀ ∝ √M und v⁴ ∝ M^{3/2}. (b) Gebundene Bereiche: mittlere freie Weglänge 0,007–0,25 Gpc (z = 0, M_min = 10⁸–10¹⁴) gegen ℓ₀ = 24,3 Gpc.

## 4 Frühe Grenze

K1b hält Ω_s = k_s zu jeder Zeit, im Strahlungszeitalter wie ΔN_eff = 2,09. K1c (Ω_s ∝ a²/(aH)²) und K2, K2c (F → 0) schalten sich früh ab.

## 5 Wo es nicht trägt

1. **Wert:** Kac legt ℓ₀ = 2Λ fest, der Wert hängt am Deckel. Von den gerechneten Deckeln ist keiner zulässig: Die atmende Zelle (a₀ = cH/2π, heute 13 % unter dem Messwert, der Wert der Horizonttemperatur, hier ohne Temperatur) hält den Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). [korr.] Ein Deckel am Schwellenradius der Hubble-Kugel oder an einem Horizont (→ c/H) ist mit dem Vakuumteil unvereinbar. Eine feste Zelle verlangte mit Λ = L₀a/2 für den Messwert L₀ = 24,3 Gpc < 27 Gpc; c²/a₀ = ℓ₀ = 2Λ ist aber keine Raumlänge: Im Modell folgt die Zelle aus a₀ = 3c⁴/(2Gρ_mL³), L = 26,9–28,1 Gpc für a₀ = 1,19–1,04·10⁻¹⁰; ohne Kreispaare ist L ≥ 0,985·d_LSS, also a₀ ≤ 1,13·10⁻¹⁰. [korr.]
2. **Zeitgesetz in K1:** a₀ ∝ 1/Λ. Im Modell ist Λ die ganze innere Faltung der Zelle, früh trug die Strahlung mehr bei: a₀ ist früh kleiner (0,76·a₀ bei der Rekombination, 0 bei der Nukleosynthese) und seit z ≈ 10 fest; es fällt nicht monoton. [korr.]
3. **Maximum in K2:** Ein inneres Maximum entsteht, sobald μ(A)/μ(X) wächst. Es bricht an drei Stellen: a₀(0) ≤ (cH₀/2π)·F(0), der Messwert verlangte F(0) = 1,15 > 1, beim Maximum im Fenster ist a₀(0) = 0,52·10⁻¹⁰; M_min ist frei; das Maximum ist 1 % hoch statt 11–26 %.
4. **Steigung heute:** Verlangt ist ein a₀, das schneller fällt als H: −0,76 … −1,56 gegen −ε₀ = −0,64 … −0,81 im DESI-Hintergrund [G]. Mit atmender Zelle und nicht schrumpfender Menge A gilt d ln a₀/d ln a ≥ −ε (heute −0,60 in K1b, −0,47 in K2). Es bräuchte eine Menge A, deren Anteil seit z ≈ 0,4 schrumpft, um 0,12–0,75 je e-Faltung; gebundene Bereiche wachsen nur.
5. **K2c** kommt der Form am nächsten (Breite 0,865, a₀(2)/a₀(0) = 0,46), bleibt aber zu flach (1,09, Steigung −0,47); der Messwert verlangte L₀ = 10,8 Gpc < 27 Gpc, die Zelle atmet nicht, L₀ und M_min sind frei; im Modell ist die Zelle L = 26,9–28,1 Gpc. [korr.]
6. **K3:** Im Volumen hängt die Rückkehrzeit an M oder ist um Größenordnungen zu kurz. Nur die eindimensionale Halbgerade gibt eine massenunabhängige Rückkehrzeit.

## Ergebnis

- Kac auf der geschlossenen Halbgeraden gibt τ_D = 2Λ/c für jede Masse und r*² = r₀·Λ; das erklärt, warum τ_D nicht von M abhängt.
- Mit der atmenden Zelle (Λ = ℓ₀/2) folgt a₀ = c²/ℓ₀ = cH/2π ohne Horizonttemperatur, mit deren Zeitgesetz a₀ ∝ H (monoton, früh ΔN_eff = 2,09); das trägt nicht (gemessen N_eff = 2,99 ± 0,17). Im Modell ist Λ heute fest, a₀ seit z ≈ 10 fest, ≈ cH₀/2π nur heute. [korr.]
- Unter der Annahme eines Zeitgesetzes a₀(z) nach dem Deckel [korr.] gibt ein Maximum bei z ≈ 0,4 nur ein wachsender gebundener Anteil F, a₀ = (cH/2π)·F; dann ist F(0) = 0,5, a₀(0) = 0,52·10⁻¹⁰, M_min frei und das Maximum 1 % hoch.
- Die verlangte heutige Steigung fällt steiler als H; in a₀ ∝ H·μ(A)/μ(X) verlangte sie eine Menge A, die seit z ≈ 0,4 schrumpft.
- Das Zeitgesetz a₀(z) folgt aus Kac nicht. Im Modell wandert a₀ nicht (seit z ≈ 10 fest); die Form der Dunklen Energie mit der Umkehr bei z ≈ 0,4 wird in der Rückstellung gesucht, als Bindungsrate dF/d ln a. [korr.]
