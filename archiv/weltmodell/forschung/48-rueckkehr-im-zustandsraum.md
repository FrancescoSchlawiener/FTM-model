# Rückkehr im Zustandsraum

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „τ_D = c/a₀ = 2π/H₀ = 79 Gyr“ → c/a₀ = 79,2 Mrd. Jahre, 2π/H₀ = 91,2 Mrd. Jahre (gleich nur bei H₀ = 77,6); Herkunft von τ_D „offen“ → τ_D = 2Λ/c (Kac); „statische Schleife am Ort, also in der Feldform ∇·[ν∇Φ_N]“ → örtliche Regel ν(|g_N|) mit dem Projektor n⊗n, die Feldform ist nur ihr wirbelfreier Teil; „warum nur eine Umlaufrichtung“ → beantwortet (Erneuerung mit χ → 1); Rate des Lévy-Kerns λ_m = a_N/c → Verlustrate Γ = a_N/c (λ_m ist die Übergaberate ω_orb/N*); Rate w je Richtung → k_s. Der Echtzeit-Absatz in Abschnitt 1 (q_t mit tH₀) setzt τ_D = 2π/H₀, also a₀ = cH₀/2π = 1,04·10⁻¹⁰.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/rueckkehr_zustandsraum.py` (Ausgabe `rueckkehr_zustandsraum.txt`)
> a₀ = 1,2·10⁻¹⁰ m/s². Beispiel: Punktmasse 6·10¹⁰ M☉, r* = 8,35 kpc, v_f = 175,8 km/s.

## Frage

q = e^{−s}, s = √(a_N/a₀), ist exakt die Wahrscheinlichkeit, dass ein Zustand aus dem Abstand L diffusiv zur Schleife zurückkehrt, bevor ein poissonscher Verlust Γ ihn nimmt: q = E[e^{−ΓT}] = e^{−L√(Γ/D)}, T Lévy-verteilt. Geprüft wird, ob eine Richtung im Zustandsraum oder der Gedächtniskern diese Rückkehr liefert.

## 1 Fokker-Planck quer zum Torus

```
D q'' + b(ρ) q' − Γ q = 0,   q(0) = 1,   s² = ΓL²/D,   ν = 1/(1 − q)        Tully-Fisher: s·ν → 1 für s → 0
```

**Drift** b = −v: −ln q = √(Pe² + s²) − Pe mit Pe = vL/2D.

| Pe | s·ν bei s = 0,1 | 0,01 | 0,001 | Tiefgrenze |
|---|---|---|---|---|
| −0,1 (vom Torus weg) | 0,466 | 0,055 | 0,0055 | a → 5,52·a_N |
| 0 | 1,051 | 1,005 | 1,0005 | a → √(a_N a₀) |
| 0,001 | 1,061 | 1,110 | 2,415 | Boden a → 2Pe·a₀ |
| 0,1 | 2,465 | 20,05 | 200,0 | Boden a → 0,2·a₀ |

Die Wurzel verlangt |Pe| ≪ s bis zum kleinsten s, also Drift genau null.

**Kontraktion** (Normalform ρ̇ = −κρ + Rauschen, κ = |λ_⊥|, Bandbreite σ² = D/κ): q = e^{y²/4}·D_{−Γ/κ}(y)/D_{−Γ/κ}(0), y = L/σ (parabolische Zylinderfunktion; Differenzenverfahren auf 6 Stellen gleich).

| L/σ | s·ν bei s = 1 | 0,1 | 0,01 | 0,001 | Boden a/a₀ |
|---|---|---|---|---|---|
| 0,01 | 1,582 | 1,054 | 1,261 | 8,07 | 0,0080 |
| 0,1 | 1,587 | 1,333 | 8,36 | 83,0 | 0,083 |
| 1 | 2,097 | 11,19 | 110,9 | 1109 | 1,109 |

Die Wurzel gilt nur für s ≫ L/σ. Darunter ist die mittlere Rückkehrzeit endlich, 1 − q ≈ Γ⟨T⟩ wird linear, und a läuft gegen den Boden ≈ √(2/π)·(L/σ)·a₀. Die kontrahierende Querrichtung trägt nicht; verlangt ist λ_⊥ → 0, ein unbeschränkt dickes Band.

**Wirkung I.** Zerfallen die Tori, diffundiert die Wirkung mit D = (K²/4)(1 − 2J₂(K) + 2J₂(K)²) (Standardabbildung, gerechnet 0,65/1,35/0,99/0,87·K²/4 bei K = 10/20/40/80). Mit K ∝ a₀/a_N hängt D selbst von a_N ab. Allgemein gilt für D ∝ a_N^α, Γ ∝ a_N^β: s ∝ a_N^{(β−α)/2}. Die Wurzel verlangt β − α = 1; Wirkungsdiffusion mit Γ ∝ a_N gibt den Exponenten 3/2.

**Verlust ∝ a_N.** λ_⊥ ist eine Drift, kein Verlust. Die Rückstellrate beim Einrasten, √(4K_d² − Δω²), ist ∝ a_N nur, wenn die Kopplung es ist. Die einzige lokale, universelle Rate proportional zu a_N ist a_N/c:

```
Γ = a_N/c = 3·|q_Δs|·c/r,     τ_D = L²/D = c/a₀ = 79,2 Mrd. Jahre,     s² = Γ·τ_D = a_N/a₀     [korr.]
```

[korr.] c/a₀ und 2π/H₀ sind verschieden: c/a₀ = 79,2 Mrd. Jahre, 2π/H₀ = 91,2 Mrd. Jahre (H₀ = 67,4); gleich wären sie nur bei H₀ = 77,6 km/s/Mpc. Die Form s² = 6π·(Δs/r)·(R_H/r) gilt nur mit τ_D = 2π/H₀.

L√(Γ/D) = √(a_N/a₀) kommt also genau dann heraus, wenn die Rückkehr τ_D = c/a₀ dauert. Keine der untersuchten Richtungen liefert τ_D; es folgt aus Kac: Die mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden ist τ_D = 2Λ/c, für jede Masse gleich, also a₀ = c²/(2Λ) (heute Λ ≈ πc/H₀); mit D = GM/c folgt L = r*. [korr.]

**Echtzeit.** Als Vorgang in der Zeit scheitert die Rückkehr zweifach:
- Seit der Zeit t sind nur Rückkehren mit T < t abgeschlossen: q_t(s → 0) = erfc(√(π/2tH₀)) = 0,012 / 0,076 / 0,21 bei tH₀ = 0,5 / 1 / 2. Bei s = 0,1 und tH₀ = 1 ist s·ν_t = 0,108 statt 1,051; auf 10 % genau verlangte es tH₀ = 420.
- Wirkt der Kern auf die Geschwindigkeit, sieht die Kreisbahn K̂(−iω) statt K̂(0) = e^{−s}, mit ωτ_D = (c/v)·s²ν:

| r/r* | 0,3 | 1 | 3 | 10 | 100 |
|---|---|---|---|---|---|
| ωτ_D | 10 567 | 2 145 | 616 | 175 | 17 |
| lg(\|K̂(−iω)\|/K̂(0)) | −30,1 | −13,8 | −7,5 | −4,0 | −1,3 |

Die Unterdrückung ändert sich über die Scheibe um 29 Größenordnungen; eine feste Stärke κ₀ gleicht das nicht aus.

**Ergebnis 1.** Die Form trägt, die Zeit nicht. Die Regel ist exakt die Rückkehr längs einer driftfreien, offenen Richtung, deren Diffusion nicht von a_N abhängt, mit Verlust a_N/c. Die Querrichtung des Torus (kontrahierend) und die Wirkungsrichtung zerfallener Tori (D hängt von a_N ab) sind das nicht. Der Ansatz bricht an der Zeitskala: τ_D = c/a₀ = 79,2 Mrd. Jahre [korr.] ist länger als das Alter und als jede Bahnperiode. Die Rückkehr kann kein Vorgang längs der Bahn sein, nur eine stationäre Ausbreitung (die Resolvente bei Γ).

## 2 Gedächtniskern

**Lévy-Kern.** Dichte erster Rückkehrzeiten mal Überlebenswahrscheinlichkeit statt e^{−λ_mΔt}:

```
K(Δt) = √(τ_D/4π)·Δt^{−3/2}·e^{−τ_D/4Δt}·e^{−ΓΔt},     ∫K dΔt = e^{−√(Γτ_D)}      (auf 10 Stellen)
```

Die Schleifenverstärkung ist genau e^{−s}, wenn Γτ_D = a_N/a₀, also mit der Verlustrate Γ = a_N/c; die Übergaberate λ_m = ω_orb/N* der Mode ξ ist eine andere Größe. [korr.] Die Zeit 1/Γ = c/a_N ist bei a_N = a₀ 79,2 Mrd. Jahre, als Nachlaufzeit gelesen weit über der dynamischen Zeit. [korr.]

**Maxwell-Element und Zwanzig-Mori.** Eliminiert man eine Kette gleicher Koordinaten ohne Rückstellung (Rate k_s je Richtung [korr.]), die an der Schleife endet und im Abstand n₀ gespeist wird, ist der Zwanzig-Mori-Kern K(u) = k_s·[e^{Au}]_{1,n₀} die Lévy-Dichte (n₀ = 50: Verhältnis 1,00000 für u ≥ τ_D; K̂(1/τ_D) = 0,367886 gegen e^{−1} = 0,367879).

| Kern | s·ν bei s = 1 | 0,1 | 0,01 |
|---|---|---|---|
| Maxwell-Element (eine Stufe) | 2,000 | 10,10 | 100,0 |
| Kette in Reihe, gerichtet (50 Stufen) | 1,591 | 10,05 | 100,0 |
| Kette mit Hin- und Rücklauf, offen | 1,582 | 1,051 | 1,005 |

Endlich viele Stufen (Prony-Reihe) geben ein rationales K̂ und damit 1 − q ∝ Γ. Die Wurzel verlangt einen Verzweigungspunkt von K̂ bei Γ = 0, also ein lückenloses Kontinuum eliminierter Moden: eine halbunendliche Kette von Maxwell-Stufen mit Hin- und Rücklauf (Diffusionsleitung), abgegriffen im Abstand L. Der Kern des druckfreien Einfalls, G(t,t′) = (3/5)t′[(t/t′)^{2/3} − t′/t] (numerisch auf 6 Stellen), wächst und hat kein endliches Integral; aus ihm entsteht keine Verstärkung q < 1.

**Kreisbahnen und Rundung.** Der Lévy-Kern hat bei Δt → 0 kein Gewicht. Es fehlt genau der Speicheranteil, der in der Exponentialform die Kreisbahn bindet; bei Δt ~ τ_D mittelt die Bahnphase den Kern weg. Auf der Geschwindigkeit bindet er daher keine Kreisbahn und rundet keine exzentrische: Speicher- und Verlustanteil sind bei der Radialfrequenz gleich stark unterdrückt. Trägt die Rückkehr, dann als statische Schleife am festen Ort, a = a_N/(1 − q(a_N)). Am Ort der Bahn gezählt sieht sie den Betrag des Gefälles dort: Das ist die örtliche Regel a = ν(|g_N|)·g_N mit dem Projektor n⊗n; die Feldform ∇·[ν∇Φ_N] trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Für Kugeln sind beide gleich. [korr.] Exzentrische Bahnen werden dann nicht gerundet, die Radialstreuung alter Scheiben bleibt, und ein Nachlauf nach Feldänderungen folgt aus der Rückkehr nicht.

## 3 Kac-Budget

Exakte Lévy-Rückkehr heißt ⟨T⟩ = ∞: Der Kac-Anteil der Schleife ist null, μ(X) = ∞. Auf einem geschlossenen Raum ist ⟨T⟩ endlich und 1 − q = Γ⟨T⟩ + … linear. Beides verträgt sich auf zwei Wegen.

**Offene Unterrichtung, geschlossen bei Λ.** q = cosh((Λ−L)s/L)/cosh(Λs/L), ⟨T⟩ = L(2Λ−L)/2D (Kette: beide exakt).

| Λ/L | s·ν bei s = 0,1 | 0,01 | 0,001 | Boden a/a₀ = 2/(2Λ/L − 1) | 10 % Abweichung bei s = |
|---|---|---|---|---|---|
| 1 | 20,08 | 200,0 | 2000 | 2 | 2,47 |
| 10 | 1,403 | 10,56 | 105,3 | 0,105 | 0,157 |
| 100 | 1,051 | 1,322 | 10,08 | 0,010 | 0,015 |
| 1000 | 1,051 | 1,005 | 1,314 | 0,001 | 0,0015 |
| ∞ | 1,051 | 1,005 | 1,0005 | 0 | – |

Die Regel gilt für s ≫ L/Λ (10 % bei s ≈ 1,5·L/Λ); Λ = L ist der geschlossene Rückweg mit 1/cosh s. Die Kontraktion quer zum Torus wirkt ebenso, mit σ statt Λ.

**Schließung der Ausbreitungsrichtung.** 𝒫 = (1/M)ΣT_{mL}, angewandt auf die abgeschirmte Ausbreitung −D∂² + Γ, summiert über Windungen. Am Ausgangspunkt gilt (Ring mit 400 Plätzen, 8 Stellen):

```
G_Ring(0)/G_Gerade(0) = Σ_{n∈ℤ} e^{−|n|s} = coth(s/2) = 1 + 2Σ_{n≥1} qⁿ = 2ν − 1
```

Der Ring ist geschlossen (Nullmode mit Gewicht 1/L, Kac), die Gerade nullrekurrent; ihr Verhältnis trägt die Wurzel für alle s, ohne Boden. Es zählt aber beide Umlaufrichtungen. Bei gleichem Tiefgrenzfall entsteht coth s statt ν(s): 2,164 / 1,313 / 1,037 gegen 2,541 / 1,582 / 1,157 bei s = 0,5 / 1 / 2. Die Regel ist die Summe über eine Umlaufrichtung (n ≥ 0); eine Ausbreitung zweiter Ordnung auf dem Ring gibt beide.

**Ergebnis 3.** Das Kac-Budget passt nicht zu einer exakt offenen Richtung, wohl aber zu einem geschlossenen Raum mit offener Unterrichtung, deren Schließungslänge Λ weit über L liegt. Die Regel gilt dann bis s ≈ L/Λ.

## 4 Prüfbare Aussagen

1. **Ende des flachen Bereichs.** Ist die Rückkehrrichtung bei Λ geschlossen (oder durch λ_⊥ auf die Bandbreite σ begrenzt), weicht die Beschleunigung ab r ≈ 0,65·r*·Λ/L um mehr als 10 % von Tully-Fisher ab und läuft gegen den Boden a₀·L/Λ; dort steigt v² ∝ r. Flache Rotationskurven bis N·r* verlangen Λ/L ≳ 1,5·N.
2. **Form am Übergang.** Eine Umlaufrichtung gibt ν(s), beide geben coth s: bei a_N = a₀ 1,582 gegen 1,313 bei gleichem Tully-Fisher-Nullpunkt, trennbar in der Beschleunigungsrelation.
3. **Örtliche Regel.** Mit Γ = a_N/c ist die Rückkehr kein Vorgang längs der Bahn: keine Rundung exzentrischer Bahnen, Vertikaldynamik in der Mittelebene ν_z² = ν·ν_z,N² (dort ist g_N,z = 0, örtliche Regel und Feldform geben dasselbe), kein Nachlauf aus der Rückkehr. Die Bahnen folgen der örtlichen Regel ν(|g_N|); die Feldform ist nur ihr wirbelfreier Teil. [korr.]

## Offen

- Welche Richtung des Zustandsraums offen, driftfrei und von a_N unabhängig ist; Querkontraktion, Wirkung und Umlauf erfüllen nicht alle drei Bedingungen.
- Woher τ_D kommt, ist beantwortet: τ_D = 2Λ/c (Kac), mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden; ≈ 2π/H₀ gilt nur heute. [korr.]
- Warum nur eine Umlaufrichtung zählt, ist beantwortet: Die Erneuerung G(mL) = q^m (jeder Rücklauf zählt einmal) legt sie zusammen mit dem Grenzwert χ → 1 des Autors fest; beide Richtungen gäben im tiefen Bereich χ → 2. Kausalität erzwingt sie nicht (beide Windungssinne sind retardiert), Holomorphie auch nicht (der Schwarz-Kern (1 + q)/(1 − q) ist ebenfalls holomorph und gibt 1 + 2n_B). [korr.]
