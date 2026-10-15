# Der Übergang am Schwellenradius r*

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: q = e^{−r*/r} als Nome des Bahntorus in Punkt 1 des Autors gekennzeichnet (für die Punktmasse der Rückkehranteil; Kern-Nome q_Δs = Δs/r); Im τ_M = r*/(2πr) als Festlegung → offene Wahl zwischen Punkt 1 und Punkt 3 (nur Schreibweise); RAR-Tabelle: „gemessene Kurve“ → „angepasste Relation“, „identisch“ → gleich nach Konstruktion, keine Bestätigung; „153 Galaxien“ [G]; Sonne Newtonsch bis zur Gezeitengrenze 1,37 pc → bis zur Glättungslänge ℓ ≈ 0,93 pc; Befund 4 und Vorschlag 4: Schwelle aus dem Atmen, a₀(z) = cH(z)/2π (bei z = 2 dreifach) → a₀ = c²/(2Λ), seit z ≈ 10 fest, ≈ cH₀/2π nur heute; Vorschlag 5: Feldform als Feldgleichung der Bewegung → Gleichung der Metrik (wirbelfreier Teil der örtlichen Regel). Die Zahlen rechnen mit a₀ = 1,2·10⁻¹⁰.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/forschung_uebergang.py`

## Ausgangslage

Die Feldgleichung hat die Form

```
r̈ = −∇Φ_bar − χ(a_N)·κ₀·∫ Δ(τ_M(r))·(e_r⊗e_r)·e^{−λ_m(t−t′)}·ṙ(t′) dt′
```

Der Faktor χ(a_N) ist die Kurzschreibweise für den Übergang zwischen zwei Grenzfällen: χ → 0 für a_N ≫ a₀ (Newton) und χ → 1 für a_N ≪ a₀ (kollektives Feld). Mit der Referenzlänge r*/3 wird die Zusatzbeschleunigung im Fernfeld a_eff = √(GM·a₀)/r.

Die Form des Übergangs dazwischen liefert die Rückkopplungssumme a = a_N/(1−q) mit q = e^{−r*/r}. Sie ist stetig und glatt.

Als Zahlenwert der Schwelle gilt a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Der Modellwert c·H₀/2π (Befund 4) ist 1,04·10⁻¹⁰ m/s² bei H₀ = 67,4 bzw. 1,08·10⁻¹⁰ m/s² bei H₀ = 70 km/s/Mpc.

## Befund 1 — Die Fassung mit der Referenzlänge r*/3 liefert die Tully-Fisher-Relation exakt

Für jede Masse gilt v_rot⁴ = G·M·a₀. Die Rechnung ergibt das Verhältnis 1,000000 für Massen von 10⁹ bis 2·10¹¹ M☉.

| M [M☉] | r* [kpc] | v_rot [km/s] |
|---|---|---|
| 10⁹ | 1,08 | 63 |
| 10¹⁰ | 3,41 | 112 |
| 6·10¹⁰ | 8,35 | 176 |
| 2·10¹¹ | 15,24 | 238 |

Diese Fassung kommt damit ohne freie Parameter aus und reproduziert die gemessene Masse-Geschwindigkeits-Beziehung der Galaxien.

## Befund 2 — Die Rückkopplungsschleifen liefern den stetigen Übergang

Der Operator 𝒫 identifiziert die Fluchtmoden „mit den inneren Rückkopplungsschleifen". Nimmt man das wörtlich, entsteht folgendes Bild:
- Jeder Umlauf durch die geschlossene Schleife führt einen Anteil q der Beschleunigung ins System zurück.
- Selbstkonsistent gilt dann:

```
a = a_N + q·a   ⟹   a = a_N · Σ_{n≥0} qⁿ = a_N / (1 − q)
```

Die q-Reihe ist also die Summe über alle Rückläufe. Das ist die Struktur einer q-Entwicklung M_k(τ_M) = Σ c_n qⁿ, hier mit allen Koeffizienten c_n = 1.

Der Nome koppelt an den Schwellenradius r* = √(3c²Δs/a₀); q ist hier der Nome des Bahntorus in Punkt 1 des Autors und für die Punktmasse der Rückkehranteil, der Kern-Nome des Modells heißt q_Δs = Δs/r [korr.]. Dann gilt:

```
q(r) = exp(−r*/r)          (Punktmasse; allgemein q = exp(−√(a_N/a₀)))
Im τ_M = r*/(2π r)        ⟹   |q| = |e^{2πiτ_M}| = e^{−r*/r}
```

Im τ_M = r*/(2πr) ist keine Festlegung, sondern Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise. [korr.]

Im starken Potential (r → 0) gilt Im τ_M = r*/(2πr) → ∞, also q → 0 und damit exakt Newton. Der Imaginärteil trägt die Bindung. Der Realteil Re τ_M = ω_rad/ω_azim trägt die Umlaufphase und geht in den Betrag |q| nicht ein.

q hängt über r* an der Faltungstiefe Δs. Ein freier Parameter tritt nicht auf.

Vergleich der Übergangsfunktionen ν = a/a_N:

| y = a_N/a₀ | angepasste Relation (RAR) [korr.] | Rückkopplung Σqⁿ |
|---|---|---|
| 100 | 1,000 | 1,000 |
| 10 | 1,044 | 1,044 |
| 3 | 1,215 | 1,215 |
| 1 | 1,582 | 1,582 |
| 0,3 | 2,371 | 2,371 |
| 0,1 | 3,689 | 3,689 |
| 0,01 | 10,51 | 10,51 |

Die Rückkopplungsform ist dieselbe Funktion wie die an die Rotationsdaten von 153 Galaxien [G] angepasste Relation (Radial Acceleration Relation), a_N/(1 − e^{−√y}); die Gleichheit der Spalten gilt nach Konstruktion und ist keine Bestätigung. [korr.] Sie ist glatt und hat das tiefe Grenzverhalten √(a_N·a₀). Damit liefert auch sie die Tully-Fisher-Relation aus Befund 1.

Im Sonnensystem ergibt sie:

| Ort | q | Zusatzbeschleunigung |
|---|---|---|
| Erdbahn | e^{−7030} ≈ 0 | 0 |
| Saturn | 2·10⁻³¹⁹ | 10⁻³²³ m/s² |
| Neptun | 4·10⁻¹⁰² | 2·10⁻¹⁰⁷ m/s² |
| Sonne bei r* | 0,37 | 7·10⁻¹¹ m/s² |

Die letzte Zeile gibt den Wert des isolierten Trichters an. Im Feld der Milchstraße bleibt die Sonne als kompakter Trichter (g(R) > a₀) bis zu ihrer Glättungslänge ℓ ≈ 0,93 pc (27 r*) Newtonsch, aus ℓ³ = GM/T mit T = √(3/2)·|E^TF| (voller spurfreier Gezeitentensor); 1,37 pc ist der klassische Jacobi-Radius. Erst dort geht sie in das kollektive Feld der Galaxie über. [korr.]

## Befund 3 — Rolle von Δ(τ_M)

Auf der imaginären Achse ist Δ an **beiden** Enden null:
- bei q → 0 (Nahfeld),
- bei q → 1 (tiefes Fernfeld), wegen Δ(i/t) = t¹²Δ(it).

| t | q | Δ(it) |
|---|---|---|
| 3 | 6,5·10⁻⁹ | 6,5·10⁻⁹ |
| 1 | 1,9·10⁻³ | 1,8·10⁻³ |
| 0,5 | 4,3·10⁻² | 1,4·10⁻² |
| 0,2 | 0,28 | 5,5·10⁻⁶ |
| 0,1 | 0,53 | 5·10⁻¹⁶ |

Δ ist deshalb die Form für das **Nahfeld** (verschwindet an der Spitze). Im Fernfeld (q → 1), wo die flache Rotationskurve den vollen Effekt braucht, tragen alle Rückläufe bei, und die Summe mit den Koeffizienten 1 trifft die Daten.

## Befund 4 — Was die Schwelle a₀ setzt

[korr.] Dieser Befund setzt voraus, dass die Schwelle mit dem Atmen mitläuft (a₀ ∝ H); das trägt nicht: Der Vakuumteil bliebe dann zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037). Im Modell gilt a₀ = c²/(2Λ) aus der Kac-Rückkehrzeit; Λ ist die ganze innere Faltung der Zelle und heute fest, a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute.

Heute trifft die Schwelle den Wert aus dem **Atmen des Torus** [korr.]:

```
a₀ ≈ c·H₀ / 2π
```

| H₀ [km/s/Mpc] | c·H₀/2π [m/s²] |
|---|---|
| 67,4 | 1,04·10⁻¹⁰ |
| 70,0 | 1,08·10⁻¹⁰ |
| 73,0 | 1,13·10⁻¹⁰ |

Das liegt innerhalb von 6–13 % der gemessenen Schwelle 1,2·10⁻¹⁰ m/s².

Lesart im FTM (nur heute, mit a₀ ≈ cH₀/2π) [korr.]: Der Einzeltrichter bleibt isoliert, solange seine lokale Neigung g größer ist als die Beschleunigung, mit der das atmende Medium an ihm zieht. Fällt g darunter, geht der Trichter in das kollektive Faltungsfeld über. Genau das ist die Aussage des Autors zu r*.

Das KAM-Bild beschreibt, was jenseits von r* geschieht: Die Tori verlieren ihre Isolation.

**Prüfbare Folgerung:** Nur wenn a₀ aus dem Atmen stammte, gälte a₀(z) = c·H(z)/2π, bei z = 2 etwa dreimal so hoch. Im Modell ist a₀ seit z ≈ 10 fest; eine Verdopplung bis z = 1 ist ausgeschlossen, frühe Galaxien haben dieselbe Schwelle. [korr.]

## Vorschlag für das Modell (wartet auf Freigabe)

1. Referenzlänge r*/3, a_eff = √(GM·a₀)/r. Der Grenzradius ist r₀ = 6Δs.
2. Der Übergang χ(a_N) ist die Rückkopplungssumme:
   ```
   a = a_N / (1 − q),   q = exp(−√(a_N/a₀)) = exp(−r*/r),   Im τ_M = r*/(2πr)
   ```
   Im τ_M = r*/(2πr) ist Punkt 1 des Autors, q dort der Nome des Bahntorus (für die Punktmasse der Rückkehranteil); die Wahl gegen Punkt 3 (|q_Δs| = Δs/r) ist offen und betrifft in der stationären Regel nur die Schreibweise. [korr.]
3. Δ(τ_M) ist die Nahfeld-Form.
4. Die Schwelle ist a₀ = c²/(2Λ) mit der Schließlänge Λ der Rückkehr (Kac), heute fest; a₀ ≈ c·H₀/2π gilt nur heute, eine Kopplung an das Atmen trägt nicht (N_eff, Krümmung). [korr.] Das KAM-Bild beschreibt den Zustand jenseits von r*.
5. Feldform (Poisson) als Gleichung der Metrik, mit derselben Übergangsfunktion. [korr.] Als Feldgleichung der Bewegung trägt sie nicht: Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N; die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Für Kugeln sind beide gleich:
   ```
   ∇²Φ_N = 4πG·ρ_bar
   ∇²Φ   = ∇·[ ν(|∇Φ_N|/a₀) · ∇Φ_N ],   ν = Σ qⁿ = 1/(1 − q)
   ```
   Die Masse ρ_bar bleibt unverändert. Die Zusatzwirkung steckt allein im Operator ν, also in der Struktur und nicht in der Quelle. Das ist genau das Beobachterprinzip des Autors.
