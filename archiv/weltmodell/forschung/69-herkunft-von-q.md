# Herkunft des kosmischen q aus der Rückkehr

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: mit H atmende Zelle (Λ = πc/H) als nicht tragend gekennzeichnet (ΔN_eff = 2,09; Λ im Modell heute fest); L2 (Q = 2π − 1) als Zahlennähe ohne Herleitung gekennzeichnet, Q = 5,364 aus der Messung; lineares Feld heute a_N/a₀ = 0,010 → gilt mit a₀ = cH₀/2π, mit a₀ = 1,2·10⁻¹⁰ 0,008–0,009. Umbenennungen: U (atmende Zelle, c²/a₀) → ℓ₀; Würfel U³ und feste Zelle U₀ ohne Zeichen (Kantenlänge 2πc/H, 27 Gpc). Das Ergebnis „Herkunft von Q offen“ stimmt mit dem Modell. Die Notiz rechnet mit a₀(z) = cH(z)/2π.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/herkunft_q.py` (Ausgabe `herkunft_q.txt`)
> Messwerte [G] aus dem Gedächtnis, ungeprüft: Planck 2018 (TT,TE,EE+lowE+lensing) ω_b = 0,02237 ± 0,00015, ω_c = 0,1200 ± 0,0012 ⇒ Q = 5,364 ± 0,065 (1,2 %), q = 0,8429 ± 0,0016; mit BAO ω_b = 0,02242 ± 0,00014, ω_c = 0,11933 ± 0,00091 ⇒ Q = 5,323 ± 0,053. Dazu H₀ = 67,4, σ₈ = 0,811, n_s = 0,965 [G].

## Frage

Linear ist die Schließungsregel das Standardmodell bei Q = q/(1 − q); der Beobachter liest die Wirkung der Mode als Dichte Q·ρ_b. Folgt q aus derselben Rückkehr, die a₀ = c²/ℓ₀ gibt: Halbgerade [0, Λ], Λ = ℓ₀/2 = πc/H, Schritte r₀ = 6Δs mit c, ⟨T_Ort⟩ = 2Λ/c? [korr.] Mit H atmend (Λ = πc/H zu jeder Zeit) trägt das nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); Λ ≈ πc/H₀ gilt nur heute.

## 1 Lesarten

```
D·G'' = Γ·G,  G(0) = 1,  G'(Λ) = 0   ⇒   q = G(L) = cosh((Λ−L)s/L)/cosh(Λs/L),   s² = Γ·τ_D = a_N/a₀
kleines s:  1 − q = (2Λ/L − 1)·s²/2                    (diskrete Kette mit Verlust: gleich auf 2·10⁻⁴)
```

**L1 Rückkehr auf der geschlossenen Halbgeraden.** Im homogenen Hintergrund ist a_N = 0, also Γ = 0 und q = 1 für jedes Λ: ν = ∞, kein endliches Q. Im linearen Feld (Effektivwert über alle k) ist a_N/a₀ = 0,010 heute (mit a₀ = cH₀/2π; mit a₀ = 1,2·10⁻¹⁰ 0,008–0,009, BBKS- bzw. Eisenstein-Hu-Spektrum, σ₈ = 0,811 [G]) [korr.] und 0,0015 bei z = 1090, also Q = 9,6 bzw. 25. Der Deckel wirkt erst bei Λk ≈ 1; für k ≥ 0,01 h/Mpc ist Λk ≥ 94.

**L2 Kac-Maß in der atmenden Zelle.** [korr.] Diese Lesart setzt die mit H atmende Zelle voraus (Λ = πc/H); das trägt nicht (ΔN_eff = 2,09). Q = 2π − 1 ist Zahlennähe ohne Herleitung; Q = 5,364 ist aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0). ν = μ(X)/μ(A), gezählt in Schritten des Ortes. Der Ort des Kosmos ist die globale Engstelle mit Δs_H = c/6H, also r₀ = c/H. Die Kette auf [0, πc/H] gibt für jede Schrittfolge ⟨T_Ort⟩·H = 2π (auf 10 Stellen), auf der Schleife vom Umfang ℓ₀ ebenso [korr.]:

```
1 − q = μ(Ort)/μ(X) = r₀/2Λ = 1/2π   ⇒   ν = 2π,   q = 0,8408,   Q = 2π − 1 = 5,283
Ω_b/Ω_m = μ(Ort)/μ(X),   Q = μ(Rückweg)/μ(Ort)
```

Die Hubble-Kugel erfüllt dabei dieselbe Beziehung wie jede Masse: r*_H² = GM_H/a₀ = π(c/H)² = r₀·Λ. Andere Maße derselben Zelle geben anderes: Ort 2c/H Q = 2,14, Ort r*_H = √π·c/H Q = 2,54, Hubble-Kugel im Würfel der Kantenlänge 2πc/H [korr.] Q = 58,2, Zählung in Umläufen ℓ₀/c [korr.] Q = 0.

**L3 Kriterium R_Q.** Lokales q = kosmisches q heißt a_N = (ln q)²·a₀ = 0,029 a₀, R_Q = 5,85 r*; das schreibt q nur um. Legt man die Zone auf den Virialradius des Kollapses (18π²·ρ̄_b(z_c), a₀ = cH₀/2π), folgt Q = 1/(e^s − 1) mit s = r*/r_vir ∝ M_b^{1/6}(1 + z_c): Q = 1,9–25 für M_b = 10⁹–10¹² M☉ und z_c = 1–6 (6·10¹⁰ M☉ bei z_c = 3: 6,1). Kein fester Wert.

**L4 Rand der Hubble-Kugel.** a_N = cH/2 = π·a₀, s = √π; geschlossen mit L = r*_H, Λ/L = √π: q = cosh(π − √π)/cosh π = 0,181, Q = 0,22 (offen 0,170). Nur mit Ω_m bzw. Ω_b als Faltung [G]: Q = 0,60 bzw. 2,09. Verlangt wäre am Rand 0,029 a₀, 108-mal weniger.

| Lesart | q | Q | zeitlich fest |
|---|---|---|---|
| L1 homogen | 1 | ∞ | – |
| L1 lineares Feld | 0,906 (heute), 0,962 (z = 1090) | 9,6 / 25 | nein |
| L2 Engstelle r₀ = c/H | 0,8408 | 5,283 | ja |
| L3 Zone des Kollapses | – | 1,9–25 | hängt an M_b, z_c |
| L4 Rand der Hubble-Kugel | 0,181 | 0,22 | ja |

## 2 Zeitabhängigkeit

Fest ist nur L2 mit atmender Zelle: Λ = πc/H und r₀ = c/H wachsen gleich, ν = 2π zu jeder Zeit. [korr.] Die mit H atmende Zelle trägt nicht (Vakuumteil zu jeder Zeit 22 %, ΔN_eff = 2,09); im Modell ist Λ heute fest, a₀ seit z ≈ 10 fest. Eine feste Zelle von 27 Gpc [korr.] gäbe Q = 5,07 heute und 145 bei z = 1090; L1 9,6 → 25; die Kac-Zeit durch das Weltalter, 2π/(H·t) − 1, 6,1 → 9,6.

Ein zeitabhängiges Q widerspricht dem exakt staubartigen Hintergrund der Schließungsregel (a³ mal Baryonen- plus gelesene Dichte fest). Die Strahlung liest Q bei z ≈ 1090 als 5,36 ± 0,06 [G]; L1 gäbe dort 25, die feste Zelle 145, das Peakmuster wäre ein ganz anderes.

L2 verlangt ω_c = (2π − 1)·ω_b = 0,1182: gegen Planck −1,5σ in ω_c, −1,3σ in Q; mit BAO −1,0σ bzw. −0,75σ [G]. Eine Prüfung auf 5σ verlangt σ(Q)/Q ≤ 0,3 %.

## 3 Zufallsprobe

Treffer im 1σ-Band der Messung (gefunden / erwartet = Dichte·2δ):

| Klasse | Werte | Q (1,2 %) | q (0,19 %) |
|---|---|---|---|
| A 2^a·3^b·5^c·π^d, Exponenten in ℤ/n, n ≤ 4 | 1,5 Mio. | 1889 / 1906 | 300 / 311 |
| C (p/q)·π^k, p, q ≤ 64, \|k\| ≤ 6 | 32 755 | 53 / 51 | 10 / 8,5 |
| D x∘y aus 20 Zahlen des Modells (1, 2, 3, 4, 6, 9, 12, 24, π, 2π, √π, π², e, φ, √2, ζ(3), π²/6, π⁴/15, π²/45, 1/(1 − e⁻¹)), ∘ ∈ +, −, ×, ÷ | 805 | 3 / 5,3 | 0 / 0,4 |

- 2π − 1 liegt 1,5 % neben Q; gleich nah oder näher liegen in D 6, in C 69 Ausdrücke. Als Zahl zählt der Treffer nichts.
- Die eigenen Lesarten mit festem, endlichem Wert (L2, L4) sind 8 Werte zwischen ν = 1 und 59. Dass einer zufällig auf 1,3 % an ν = 6,364 liegt, hat p = 0,04–0,05. Eine einzige zwingende Lesart hätte p = 0,006; p ≤ 10⁻³ verlangte Treffer und Messung auf 0,2 %.

## Wo es nicht trägt

1. **Mechanismus:** In der Galaxie entsteht 1 − q aus dem Verlust Γ = a_N/c auf dem Rückweg (q = ⟨e^{−ΓT}⟩). Im homogenen Hintergrund ist Γ = 0, dieselbe Gleichung gibt q = 1. L2 braucht eine zweite Zuordnung, 1 − q = μ(Ort)/μ(X), die aus der Rückkehrgleichung nicht folgt; für eine einzelne Masse gäbe sie ν = ℓ₀/6Δs [korr.] (Milchstraße 5·10¹², Sonne 3·10²³) statt 1–30.
2. **Wahl:** Maß des Ortes (c/H, nicht 2c/H oder r*_H) und Zähleinheit (Schritt 1/H, nicht Umlauf 2π/H) sind gesetzt; andere Maße derselben Zelle geben Q = 0 bis 58.
3. **Wert:** 1,3σ unter Planck, 0,75σ unter Planck + BAO [G]; die Zufallsprobe schlägt die Schwelle nicht (p ≈ 0,05).
4. L1, L3 und L4 geben kein festes Q: L1 ist homogen unendlich und im linearen Feld zeitabhängig, L3 hängt an Masse und Kollapszeit, L4 gibt Q = 0,22 statt 5,36.

## Ergebnis

- Aus der Rückkehrgleichung selbst (Verlust auf der bei Λ = ℓ₀/2 geschlossenen Halbgeraden) [korr.] folgt kein kosmisches q.
- Aus dem Kac-Maß derselben Halbgeraden folgt mit dem Ort r₀ = 6Δs_H = c/H, nur bei mit H atmender Zelle [korr.]: ν = 2π, Q = 2π − 1 = 5,283, zeitlich fest, 0,75–1,3σ unter der Messung [G]. [korr.] Zahlennähe ohne Herleitung, trägt nicht.
- Ohne freie Zahl trägt das nicht: Die Zuordnung 1 − q = μ(Ort)/μ(X) ist eine eigene Setzung, Maß und Zähleinheit sind gewählt, und p ≈ 0,05 liegt nicht deutlich unter der Zufallsschwelle. Die Herkunft von Q bleibt offen; Q = 2π − 1 ist bei σ(Q)/Q ≤ 0,3 % entscheidbar.
