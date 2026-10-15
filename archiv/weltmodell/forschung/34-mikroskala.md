# Die Mikroskala: Ladung, Massenverhältnis, Feinstrukturkonstante

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: gemeinsame Zahl Im τ_X für α und α_G über 1/α_em(M_X) = (8/3)/α_X (Ergebnis B, 3.2, 3.3, 4) → diese Beziehung ist nicht Teil des Modells; die zweite Gleichung α = f(Im τ_X) fehlt; im Standardmodell kreuzen sich die Kopplungen nicht, die MSSM-Rechnung (einschleifig) gibt 1/α_em(M_Z) = 126,68 statt gemessen 127,95 [G]. [G] an m_u, m_d, v, α_s(M_Z), α(0), α(M_Z). Hinweis zum Bahn-Nome |q| = e^{−r*/r} nach den Bezeichnungsregeln. Umbenannt (ohne Marke): q und Q als Ladung → Q_el (q_n → Q_el,n); w (Windungszahl) → n_w; w (Nahtbreite) → b_N; ν (Verdrehung) → ϑ; χ (Eichfunktion) → λ_g; g (Polstärke) → g_m; Eichkopplungen α₁, α₂, α₃ → α_Y, α_W, α_s (die Hyperladungskopplung ohne GUT-Normierung α_Y → α′); Λ (QCD-Skala) → Λ_QCD.

> Rechnung: `rechnungen/mikroskala.py`

**Worum es geht.** Der Autor möchte die Mikroskala geometrisch schließen, und zwar mit drei Herleitungen:
- (C) Ladung als Windungszahl eines Faltungsumlaufs,
- (A) m_p/m_e ≈ 1836,15 als Verhältnis geodätischer Eigenwerte zweier Tori,
- (B) α ≈ 1/137,036 als stabiler Resonanzpunkt der Modul-Symmetrie.

Diese Notiz prüft alle drei streng. Beobachter-Regel: Ladung und Masse werden als Topologie bzw. Geometrie des Mediums behandelt, nicht als Substanz.

## Ergebnis vorab

| Teil | Aussage | Kategorie |
|---|---|---|
| C | Ladung = Windungszahl n ∈ ℤ der Phase um den Faltungsumlauf; Q_el = n·e; Dirac-Bedingung; Chern-Zahl; nicht-ganzzahlige Windung kostet eine Schnittlinie mit Energie ∝ Länge | **hergeleitet** (mathematischer Satz) |
| C | Ladung = Impuls um den kompakten Kreis, R_KK = 2ℓ_P/√α = 23,41 ℓ_P | **hergeleitet** (Kaluza-Klein), mit ehrlicher Grenze beim Elektron |
| C | Drittelladungen: Q_el + t/3 ∈ ℤ über das Z₃-Zentrum der Farbe; frei sind nur Zustände mit ganzer Gesamtwindung | **hergeleitet** als Konsistenzregel; die Drittel selbst sind Eingabe |
| A | Dimensionstransmutation Λ_QCD = M·e^{−2π/(β₀α_s)} hat exakt die Form eines Nomes; α_G = C²\|q̃\|² | **Brücke zu etablierter Physik** |
| A | Gravitation ist für Protonen schwach, weil 1/α_s(M_Pl) = 53 im Exponenten steht | **Brücke**; die Eingabe α_s(M_Pl) bleibt offen |
| A | m_p/m_e | **offen**: verbindet den QCD-Nome (Proton) mit der Yukawa-Kopplung (Elektron) |
| B | Nach Erklärung von α_G bleibt für α die Vereinheitlichungsrelation 1/α_em = (8/3)/α_X bei M_X. Diese Beziehung ist nicht Teil des Modells; für eine gemeinsame Zahl Im τ_X fehlt die zweite Gleichung α = f(Im τ_X); im Standardmodell kreuzen sich die Kopplungen nicht, die MSSM-Rechnung (einschleifig) gibt 1/α_em(M_Z) = 126,68 statt gemessen 127,95 [G]. [korr.] | **Brücke** (hypothetische Vereinheitlichung) |

---

## 1 Ladung als Windungszahl (C)

### 1.1 Satz 1: Windung auf dem Faltungsumlauf ist ganzzahlig

**Satz.** Sei ψ = |ψ|·e^{iθ} ein stetiges komplexes Feld auf einem geschlossenen Umlauf S¹ (Umfang L) mit ψ ≠ 0. Dann ist

```
n_w = (1/2π) ∮ dθ = (1/2πi) ∮ dψ/ψ   ∈ ℤ.
```

**Beweis.** θ lässt sich längs des Umlaufs stetig fortsetzen (Hochhebung über e^{iθ}). Eindeutigkeit des Feldes verlangt ψ(L) = ψ(0), also e^{iθ(L)} = e^{iθ(0)} und damit θ(L) − θ(0) ∈ 2πℤ. ∎

**Zusatz (Stabilität).** n_w ist unter stetiger Verformung konstant, denn eine stetige ganzzahlige Größe kann nicht springen. n_w ändert sich nur, wenn ψ irgendwo durch null geht. Das ist ein Kern, an dem die Faltungsamplitude verschwindet.

Probe im Skript: Ein Feld mit Phasenlauf 2πϑ hat die Randauslenkung |ψ(L) − ψ(0)| = 2|sin πϑ|. Für ϑ = 0, 1, −2 ist sie null, für ϑ = 1/2 ist sie 2,0000 und für ϑ = 1/3 ist sie 1,7321.

### 1.2 Satz 2: Torus, Chern-Zahl und Dirac-Bedingung

Auf einem Faltungstorus T² mit Fläche A und konstantem Feld B wählt man A_y = Bx. Die Verschiebung um eine Periode L₁ in x ist nur bis auf eine Eichtransformation dieselbe:

```
ψ(x+L₁, y) = e^{i Q_el B L₁ y/ħ} ψ(x, y),        ψ(x, y+L₂) = ψ(x, y).
```

Die beiden Übergänge sind nur verträglich, wenn der Umlauf um die Ecke des Torus die Phase 1 liefert (Kozykel-Bedingung):

```
e^{i Q_el B L₁ L₂/ħ} = 1     ⟺     Q_el·Φ = 2πnħ     ⟺     c₁ = Q_el·Φ/h = n ∈ ℤ.
```

- c₁ ist die erste Chern-Zahl des U(1)-Bündels über dem Torus.
- Das Flussquant ist h/e = 4,135668·10⁻¹⁵ Wb.
- Im Skript ergibt die Kozykel-Phase für Φ = h/e und 2h/e den Wert +1, für Φ = h/(2e) aber −1. Der halbe Fluss ist also nicht erlaubt.

**Dirac-Bedingung.** Um einen Monopol mit Fluss Φ_m legt man zwei Eichpotentiale A_N, A_S auf die Nord- und die Südhalbkugel. Auf dem Äquator unterscheiden sie sich um einen Gradienten ∇λ_g, und die Übergangsfunktion e^{iQ_el λ_g/ħ} muss auf dem Äquatorkreis eindeutig sein. Das ist Satz 1, angewandt auf diesen Kreis:

```
Q_el ∮(A_N − A_S)·dl = Q_el·Φ_m = 2πnħ    (SI, Φ_m in Wb)
e·g_m = nħc/2                             (Gauß, Polstärke g_m),   für n = 1: g_m/e = 1/(2α) = 68,518
```

**Im Modell:** Die Ladung eines Faltungsknotens ist die Windungszahl seiner Phase um den Umlauf. Ein Knoten mit n ∈ ℤ schließt ohne Randauslenkung und ist topologisch stabil. Er kann nur zerfallen, wenn die Amplitude irgendwo auf null geht.

### 1.3 Nicht-ganzzahlige Windung kostet eine Schnittlinie mit Energie ∝ Länge

Man zwingt dem Torus eine Verdrehung ϑ auf, also einen Phasenlauf 2πϑ pro Umlauf in x, mit der Naht ψ(x+L₁) = e^{2πiϑ}ψ(x). Die Energie ist die Ginzburg-Landau-Form

```
E = K ∫ |∇ψ|² + (λ/4) ∫ (|ψ|² − 1)².
```

**Satz 3.** Es gilt E_min = L₂·e₁(ϑ). Dabei ist L₂ die Länge der Naht (Schnittlinie), und e₁(ϑ) ist die minimale Energie einer einzelnen x-Linie. Es gilt e₁(ϑ) = 0 genau für ϑ ∈ ℤ.

**Beweis.**
- Der y-Gradient ist nicht negativ. Also ist E größer oder gleich der Summe der Linienenergien, und das ist mindestens L₂·e₁(ϑ).
- Die in y konstante Lösung erreicht diese Schranke.
- Eine Linie mit Energie null muss konstant sein und |ψ| = 1 haben. Das ist mit der Naht nur verträglich, wenn e^{2πiϑ} = 1 gilt, also ϑ ∈ ℤ.
- Wegen der Kompaktheit wird das Minimum angenommen, also ist e₁(ϑ) > 0 für ϑ ∉ ℤ. ∎

Numerisch (Gitter 24 × N_y, K = λ = 1, Gradientenfluss aus Zufallsstart):

| ϑ | N_y = 8 | N_y = 16 | N_y = 32 | E/N_y | analytisch (\|ψ\| = 1) |
|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 | 0 |
| 1/4 | 0,8187 | 1,6373 | 3,2746 | 0,1023 | 0,1028 |
| 1/3 | 1,4501 | 2,9002 | 5,8005 | 0,1813 | 0,1827 |
| 1/2 | 3,2290 | 6,4579 | 12,9159 | 0,4036 | 0,4106 |
| 1, 2 | 0 | 0 | 0 | 0 | 0 |

Die Energie verdoppelt sich mit der Länge der Schnittlinie und verschwindet genau für ganzzahliges ϑ. Die Werte liegen etwas unter der Formel mit starrem |ψ| = 1, weil die Amplitude ein wenig nachgibt.

Wird die Verdrehung nicht über den Umlauf verteilt, sondern in einer Naht der Breite b_N gesammelt, so ist die Spannung σ = K|δ|²/b_N mit |δ| = 2|sin πϑ|. Für ϑ = 1/2 ergeben sich 4, 8 bzw. 40 K bei b_N = 1, 0,5 bzw. 0,1. Für b_N → 0 divergiert σ, ein scharfer Schnitt ist also unendlich teuer.

Physikalisch ist die Verdrehung ϑ dasselbe wie ein Aharonov-Bohm-Fluss ϑ·h/e durch den Ring. Die Energie ist periodisch in ϑ mit Minima bei ganzen Zahlen (Byers-Yang, aus dem Gedächtnis).

### 1.4 Kaluza-Klein: Ladung als Impuls um den kompakten Kreis

Die 5D-Metrik ist ds² = g_μν dx^μ dx^ν + (dy + κA_μ dx^μ)², mit y ∈ [0, 2πR). Man integriert die Einstein-Hilbert-Wirkung (c³/16πG₅)∫R₅ über y und setzt G = G₅/(2πR). Dann entsteht (c³/16πG)∫(R₄ − κ²F²/4). Gleichsetzen mit der Maxwell-Wirkung −(1/4μ₀c)∫F² liefert

```
κ² = 16πG ε₀ / c².
```

Ein Feld auf dem Kreis ist nur eindeutig als Fourier-Reihe Σ φ_n(x) e^{iny/R} mit **n ∈ ℤ**, nach Satz 1. Sein Impuls ist p₅ = nħ/R. Die kovariante Ableitung koppelt p₅ an A_μ mit der Ladung

```
Q_el,n = κ·p₅ = n · √(16πGε₀) · ħ/(c·R)       (SI)
Q_el,n = 2n · √G · ħ/(c·R)                     (Gauß)
```

Für Q_el,1 = e folgt

```
R_KK = √(16πGε₀)·ħ/(c·e) = 2ℓ_P/√α = 3,7841·10⁻³⁴ m = 23,412 ℓ_P,        Umfang 2πR_KK = 2,3776·10⁻³³ m.
```

**Umformung in der Sprache der Längenleiter.** Für die n = 1-Mode gilt λ̄₁ = R_KK, und aus Δs·λ̄ = ℓ_P²/3 folgt Δs₁ = 2,3011·10⁻³⁷ m. Dann ist

```
α = (2ℓ_P/R_KK)² = 12·Δs₁/λ̄₁ = 0,0072973526,            im Vergleich zu   α_G = 3·Δs_p/λ̄_p.
```

Das ist eine exakte Umformung und keine Herleitung. Die Frage „warum α" wird damit zur Frage „warum ist der Ladungskreis 23,41 Planck-Längen groß". Streng gilt die Relation mit der Kopplung an der Kompaktifizierungsskala. Hier ist mit α(0) gerechnet.

**Ehrliche Grenze.** Die Mode n = 1 hat die Masse ħ/(cR) = (√α/2)·m_Pl = 5,2147·10¹⁷ GeV. Ihr Verhältnis von elektrischer Abstoßung zu Gravitation ist e²/(4πε₀Gm₁²) = 4. Das ist die Extremalitätsbedingung Q_el² = 4GM² der Kaluza-Klein-Theorie (aus dem Gedächtnis). Beim Elektron ist dieses Verhältnis 4,1656·10⁴². Sein e/m ist also 1,02·10²¹-mal größer als das der Kaluza-Klein-Mode. Damit gilt:
- Die **Ladung** des Elektrons kann eine Windungszahl sein.
- Seine **Masse** kann nicht aus dem Impuls um denselben Kreis stammen.

Ladung (Topologie) und Masse (Geometrie bzw. Dynamik) haben verschiedene Ursprünge.

### 1.5 Chiralität = Vorzeichen der Windung

Kehrt man die Orientierung des Umlaufs um (y → −y), so wird n zu −n. Das ist in Kaluza-Klein genau die Ladungskonjugation (zusammen mit A → −A). Die „Chiralität des Wirbeleinstroms" ist also das Vorzeichen von n: Links- und rechtsdrehende Knoten sind Teilchen und Antiteilchen. n = 0 ist neutral.

### 1.6 Drittelladungen der Quarks

Quarks haben Q_el = 2/3 bzw. −1/3. Eine Windung von 1/3 um den U(1)-Umlauf hätte die Randauslenkung e^{2πi/3} − 1 ≠ 0. Sie ist trotzdem konsistent, weil das Quark zugleich im Farbraum sitzt. Das Zentrum Z₃ von SU(3) wirkt auf ein Triplett als ω = e^{2πi/3}. Die volle Drehung in U(1) zusammen mit dem Zentrumselement ist die Identität, wenn gilt:

```
e^{2πiQ_el}·ω^t = 1     ⟺     Q_el + t/3 ∈ ℤ,         t = (Quarks − Antiquarks) mod 3   (Trialität)
```

Das Skript prüft die Regel für u, d, ū, d̄, e⁻, ν, W⁺, Gluon, p, n, π⁺, Δ⁺⁺, Ω⁻ und das Diquark uu. Sie gilt für alle, und genau die farbneutralen Zustände (t ≡ 0) haben ganzzahliges Q_el.

**Im Modell:** Die Randauslenkung der U(1)-Windung 1/3 wird durch eine Drittel-Verdrehung im Farbraum aufgehoben. Die **Gesamtwindung** des kombinierten Bündels, (U(1)×SU(3))/Z₃, ist geschlossen. Ein einzelnes Quark lässt aber eine unkompensierte Farbnaht zurück.

Genau so verhält sich der Einschluss: Zwischen Farbladungen bildet sich ein Flussschlauch mit Energie ∝ Länge. Die Saitenspannung ist σ ≈ 0,18 GeV² = 0,912 GeV/fm (Gitter-QCD, aus dem Gedächtnis). Das ist dieselbe Struktur wie in Satz 3: eine Schnittlinie mit Energie proportional zur Länge.

**Ehrlich:**
- Die Regel Q_el + t/3 ∈ ℤ ist hergeleitet.
- Dass die Drittel überhaupt auftreten, ist eine Eingabe aus dem Standardmodell (Hyperladungen).
- Der Einschluss selbst ist in der QCD nicht bewiesen, aber numerisch gesichert.
- Die Lesart des Flussschlauchs als Farbnaht ist eine Deutung im Modell.

---

## 2 Massenverhältnis (A)

### 2.1 Stand der Physik

- **Proton:** m_p = 938,272 MeV. Die Ruhemassen der Valenzquarks (2m_u + m_d = 9,02 MeV [G] im MSbar-Schema bei 2 GeV) machen **0,96 %** davon aus. Der Rest ist Feldenergie der Gluonen und Bewegungsenergie der Quarks („Masse ohne Masse", Wilczek). Es gilt m_p = C·Λ_QCD, wobei C eine reine Zahl aus der QCD-Dynamik ist, die Gitter-QCD berechnet. Feinere Zerlegungen (σ-Term, Anteil der s-Quarks) verschieben den Quarkanteil um einige Prozent. Das Bild bleibt gleich.
- **Elektron:** m_e = y_e·v/√2 mit v = 246,22 GeV [G]. Das ergibt y_e = 2,935·10⁻⁶, eine freie Yukawa-Kopplung.

Die beiden Massen haben also **verschiedene Ursprünge**: m_p ist dynamisch (Einschluss), m_e ist parametrisch (Higgs).

### 2.2 Die Dimensionstransmutation ist ein Nome

Auf Einschleifen-Niveau gilt μ·dα_s/dμ = −(β₀/2π)·α_s² mit β₀ = 11 − 2n_f/3. Daraus folgt

```
1/α_s(μ) = 1/α_s(M) + (β₀/2π)·ln(μ/M),         Λ_QCD = M·exp(−2π/(β₀·α_s(M)))     (RG-invariant)
```

In der Standardform der QCD-Kopplung τ = θ/2π + 4πi/g² = θ/2π + i/α_s ist das Instanton-Gewicht exakt ein Nome:

```
q = e^{2πiτ} = e^{iθ}·e^{−2π/α_s},        (Λ_QCD/M)^{β₀} = |q(M)|.
```

Mit der Umskalierung τ̃ = τ/β₀ wird daraus genau die verlangte Form:

```
Λ_QCD/M = |q̃| = e^{−2π·Im τ̃},        Im τ̃ = 1/(β₀·α_s(M)).
```

- **Skalenwechsel ist eine Verschiebung von Im τ̃:** Im τ̃(μ) = Im τ̃(M) + ln(μ/M)/2π.
- **Re τ = θ/2π ist eine Umlaufphase:** θ ist physikalisch 2π-periodisch, also ist τ → τ + 1 eine exakte Symmetrie (T). Wie beim Bahn-τ trägt der Imaginärteil die Bindung (hier: die Einschlussskala) und der Realteil eine Phase, die keine Bindung trägt. Gemessen ist θ verschwindend klein (Neutronen-Dipolmoment), also Re τ ≈ 0.
- **Ehrliche Grenze:** Die Inversion S: τ → −1/τ ist **keine** Symmetrie der QCD. Die volle Modulgruppe wirkt nur in speziellen Theorien (N = 4 Super-Yang-Mills, Montonen-Olive). Geteilt wird die Nome-Form und T, nicht die ganze Modulgruppe.
- **Formale Entsprechung:** Beim Bahn-τ gilt |q| = e^{−r*/r} (Punkt 1 des Autors; für die Punktmasse ist dieses q der Rückkehranteil, der Kern-Nome des Modells heißt q_Δs = Δs/r) [korr.], hier |q̃| = Λ_QCD/M. Die Rolle von r*/r übernimmt ln(M/Λ_QCD). Das ist eine formale Parallele, keine Identität.

### 2.3 Warum Gravitation schwach ist, in der Sprache des Nomes

```
α_G = (m_p/m_Pl)² = C²·(Λ_QCD/m_Pl)² = C²·|q̃(m_Pl)|² ≈ exp(−4π/(β_eff·α_s(M_Pl))).
```

Lauf von α_s(M_Z) = 0,1180 [G] nach oben, mit den Schwellen m_c, m_b, m_t:

| | 1 Schleife | 2 Schleifen |
|---|---|---|
| α_s(M_GUT = 2·10¹⁶ GeV) | 0,02206 (1/45,33) | 0,02182 (1/45,83) |
| α_s(M_Pl) | 0,01906 (1/52,48) | **0,01886 (1/53,02)** |

- Rein einschleifig mit Schwellen ist Λ_QCD,3 = 142,51 MeV. Man erhält denselben Wert, ob man von α_s(m_c) oder von α_s(M_Pl) ausgeht; das ist die RG-Invarianz. Daraus folgt C₁ = m_p/Λ_QCD,3 = 6,584. Der Zahlenwert von C hängt von der Konvention für Λ_QCD ab: Mit der Einschleifen-Form an m_c, aber zweischleifig gelaufenem α_s, ergibt sich 196,7 MeV bzw. C = 4,770.
- Der Nome am Planck-Punkt ist |q̃| = Λ_QCD,3/M_Pl = 1,167·10⁻²⁰, mit Im τ̃_eff = 7,305. Damit gilt α_G = C₁²|q̃|² = 5,906·10⁻³⁹, wie gemessen. Das ist per Konstruktion so, denn C ist aus m_p genommen.
- Die Schwellen zählen: exp(−4π/(β₀α_s(M_Pl))) ergibt 1,5·10⁻³² mit β₀ = 9 und 1,2·10⁻⁴¹ mit β₀ = 7. Das effektive β₀ ist **7,49** (1 Schleife) bzw. **7,57** (2 Schleifen) und liegt zwischen den Werten für sechs und drei Quarksorten.
- **Empfindlichkeit** (2 Schleifen, C fest): Ändert man α_s(M_Pl) um ±1 %, ändert sich m_p um +41,1 % bzw. −29,8 %. Bei ±5 % wird m_p zu 5710 MeV bzw. 146 MeV.

**Antwort auf die Frage:** Das gemessene m_p verlangt α_s(M_Pl) = 0,01886 bzw. α_s(M_GUT) = 0,02182. Die winzige Zahl 10⁻³⁸ ist die Exponentialfunktion einer mäßigen Zahl:

```
ln(1/α_G) = 2·ln(m_Pl/m_p) = 2 · 44,01,        4π/(β_eff·α_s(M_Pl)) mit 1/α_s(M_Pl) = 53.
```

**Lesart im Modell.** Gravitation ist nicht „schwach", denn G ist keine Kopplung. α_G ist klein, weil der Proton-Torus durch einen Nome tief in der Spitze (Im τ̃ ≈ 7,3) festgelegt ist. Die Masse ist ein Bindungsmaß, keine Substanz.

**Offen bleibt:** Warum ist 1/α_s(M_Pl) ≈ 53? Die Hierarchie ist damit auf eine Zahl der Größenordnung 50 zurückgeführt, aber nicht erklärt.

### 2.4 Das Massenverhältnis m_p/m_e

Auf einem Umlauf der Länge L hat der Laplace-Operator die Eigenwerte k = 2π|n|/L, also m = ħk/c und λ̄ = L/(2π|n|). Für zwei Tori folgt

```
m_p/m_e = (|n_p|/|n_e|)·(L_e/L_p)
```

Das Massenverhältnis ist damit ein Verhältnis der Umlauflängen und Quantenzahlen. Festgelegt wird es erst, wenn beide Umläufe aus dem Modell bestimmt sind. m_p und m_e stammen aus verschiedenen Mechanismen: das Proton aus dem QCD-Nome, das Elektron aus der Yukawa-Kopplung. Ein gemeinsames Eigenwertproblem müsste y_e mit e^{−2π/(β₀α_s)} verknüpfen. Das ist die offene Aufgabe.

### 2.5 Kategorien für (A)

| Aussage | Kategorie |
|---|---|
| m_p ≈ C·Λ_QCD, Quarkmassen ≈ 1 % | etablierte Physik |
| Λ_QCD/M = \|q̃\|, Im τ̃ = 1/(β₀α_s); T-Symmetrie = θ-Periodizität | **Brücke** (exakte formale Identität) |
| α_G = C²\|q̃(m_Pl)\|², 1/α_s(M_Pl) = 53 erklärt 10⁻³⁸ | **Brücke**; Eingabe α_s(M_Pl) offen |
| m_p/m_e als Verhältnis zweier Umläufe | **offen**; Umlauflängen und Quantenzahlen müssen aus dem Modell folgen |

---

## 3 Feinstrukturkonstante (B)

### 3.1 α läuft

α(0) = 1/137,035999177 [G] mit einer relativen Unsicherheit von 1,5·10⁻¹⁰. Bei der Z-Masse gilt α(M_Z) = 1/127,95 [G] (MSbar), der Unterschied ist Δ(1/α) = 9,09.
- Die Leptonen allein (einschleifig) schieben 1/α auf 132,73.
- Der Rest stammt aus Hadronen und schweren Teilchen.

α(0) enthält also die Massenverhältnisse m_e, m_μ, m_τ und die hadronischen Schwellen. Eine rein geometrische Formel, die α(0) auf 10⁻¹⁰ trifft, müsste diese Massen mitenthalten. Eine geometrische Aussage über α gilt deshalb an einer hohen Skala und muss diese Skala nennen.

### 3.2 Konstruktiv: Was bleibt für α, wenn α_G erklärt ist?

**Carter.** Die Beziehung α_G ≈ α¹²(m_e/m_p)⁴ nach α aufgelöst:

```
α ≈ α_G^{1/12}·(m_p/m_e)^{1/3} = 1/125,24          (gemessen 1/137,036; Faktor 0,914; in α_G Faktor 2,944)
```

Mit α_G = C²|q̃|² aus Teil 2 wird das zu α ≈ (C|q̃|)^{1/6}·(m_p/m_e)^{1/3}. α wäre damit an den QCD-Nome und an die Yukawa-Kopplung des Elektrons gebunden, aber nur auf etwa 9 %. Carter ist eine Größenordnungsbeziehung (Sternhüllen) und kann α nicht auf 10⁻¹⁰ festlegen.

**Vereinheitlichung (die schärfere Beziehung).** Mit der Hyperladungskopplung α′ und α_Y = (5/3)α′ (GUT-Normierung) gilt 1/α_em = 1/α′ + 1/α_W. Bei Vereinigung (α_Y = α_W = α_X) auf der Skala M_X folgt

```
1/α_em(M_X) = (8/3)/α_X.
```

Einschleifig, ab M_Z = 91,19 GeV:

| | Kreuzung α_Y = α_W | Kreuzung α_W = α_s | 1/α_X |
|---|---|---|---|
| Standardmodell | 10^13,01 GeV | 10^16,99 GeV | keine Vereinigung (≈ 47) |
| MSSM (hypothetisch, Supersymmetrie nicht nachgewiesen) | 10^16,30 GeV | 10^16,36 GeV | 24,31 |

Umgekehrt sagt die MSSM-Vereinigung aus α_s(M_Z) und sin²θ_W den Wert 1/α_em(M_Z) = **126,68** voraus, gemessen ist 127,95 [G]. Die Rechnung ist einschleifig und ohne Schwellen.

**Folgerung:** Wenn es eine Vereinigung gibt, bestimmt **eine** Zahl, Im τ_X = 1/α_X, zwei Dinge:
- α(M_Z) logarithmisch, über den Lauf,
- α_G exponentiell, über den QCD-Nome.

Dazu kommen das geladene Spektrum und y_e. Die Frage nach α ist dann die Frage nach dem Imaginärteil von τ an der Faltungsskala. In der Kaluza-Klein-Form ist das die Frage nach dem Radius des Ladungskreises in Planck-Längen, α = (2ℓ_P/R_KK)².

[korr.] Diese Beziehung ist nicht Teil des Modells: Für eine gemeinsame Zahl Im τ_X fehlt die zweite Gleichung α = f(Im τ_X). Im Standardmodell kreuzen sich die drei Kopplungen nicht (Tabelle oben), und die MSSM-Rechnung (einschleifig) gibt 1/α_em(M_Z) = 126,68 statt gemessen 127,95 [G]. Die Folgerung gilt nur unter einer Vereinigung, die das Modell nicht enthält.

### 3.3 Kategorien für (B)

| Aussage | Kategorie |
|---|---|
| α läuft; α(0) enthält Leptonen- und Hadronenschwellen | etablierte Physik |
| α = (2ℓ_P/R_KK)² = 12Δs₁/λ̄₁ | exakte Umformung (Kaluza-Klein) |
| 1/α_em(M_X) = (8/3)/α_X; ein Im τ_X für α und α_G | **Brücke** zu (hypothetischer) Vereinheitlichung; nicht Teil des Modells, die zweite Gleichung α = f(Im τ_X) fehlt, im Standardmodell keine Kreuzung, MSSM (einschleifig) 126,68 statt 127,95 [G] [korr.] |

---

## 4 Gesamtergebnis

**Streng hergeleitet (C).** Die Ladung ist im Modell topologisch:
- Die Phase eines Faltungsknotens um seinen Umlauf hat eine ganzzahlige Windung n. Die Ladung ist Q_el = n·e, und das Vorzeichen ist die Chiralität des Einstroms.
- Auf dem Torus ist n die erste Chern-Zahl. Die Dirac-Bedingung Q_el·Φ_m = 2πnħ ist derselbe Satz auf dem Äquatorkreis.
- Eine nicht-ganzzahlige Windung erzwingt eine Schnittlinie mit Energie ∝ Länge.
- Drittelladungen sind nur zusammen mit einer Farbverdrehung erlaubt (Q_el + t/3 ∈ ℤ). Frei sind nur farbneutrale, ganzzahlig geladene Zustände.
- In Kaluza-Klein-Form ist der Ladungskreis R_KK = 2ℓ_P/√α = 23,41 ℓ_P groß.

**Echte Brücke (A, B).** Die Dimensionstransmutation ist exakt ein Nome, |q̃| = Λ_QCD/M = e^{−2π Im τ̃} mit Im τ̃ = 1/(β₀α_s). Der Realteil ist die θ-Phase mit τ → τ + 1.
- Das erklärt in der Sprache des Modells, warum Gravitation für Nukleonen schwach ist: α_G = C²|q̃(m_Pl)|² mit 1/α_s(M_Pl) = 53.
- Nur über eine Vereinheitlichung hinge α an derselben Zahl Im τ_X. Diese Beziehung ist nicht Teil des Modells; die zweite Gleichung α = f(Im τ_X) fehlt, im Standardmodell kreuzen sich die Kopplungen nicht, und die MSSM-Rechnung (einschleifig) gibt 126,68 statt gemessen 127,95 [G]. [korr.]

**Offen.** m_p/m_e verknüpft zwei Mechanismen (Nome der QCD und Yukawa-Kopplung). α ist eine Größe an einer hohen Skala und hängt über den Lauf an einer Zahl Im τ_X.

**Nächster sinnvoller Schritt.** Im Modell sollte der Faltungsumlauf benannt werden, der die U(1)-Windung trägt, und sein Radius sollte aus der Längenleiter bestimmt werden. Gesucht ist **eine** Größe: Im τ_X, gleichwertig R_KK/ℓ_P an der Faltungsskala. Jeder geometrische Kandidat hat dann einen harten, zweifachen Test. Aus dieser einen Eingabe müssen zugleich
- 1/α(M_Z) = 127,95 über den logarithmischen Lauf und
- α_G = 5,906·10⁻³⁹ über den exponentiellen QCD-Nome

folgen. Ein Kandidat muss beide Werte ohne Anpassung treffen. Die Yukawa-Kopplung y_e bleibt dabei eine eigene, offene Eingabe.

---

## Literatur

Alle Angaben aus dem Gedächtnis; Seitenzahlen vor einer Veröffentlichung prüfen.

- Kaluza, Th. (1921): Zum Unitätsproblem der Physik. Sitzungsber. Preuß. Akad. Wiss. Berlin, 966–972.
- Klein, O. (1926): Quantentheorie und fünfdimensionale Relativitätstheorie. Z. Phys. 37, 895–906.
- Dirac, P. A. M. (1931): Quantised singularities in the electromagnetic field. Proc. R. Soc. Lond. A 133, 60–72.
- Byers, N.; Yang, C. N. (1961): Theoretical considerations concerning quantized magnetic flux in superconducting cylinders. Phys. Rev. Lett. 7, 46.
- Carter, B. (1974): Large number coincidences and the anthropic principle in cosmology. IAU Symp. 63, 291–298.
- Wu, T. T.; Yang, C. N. (1975): Concept of nonintegrable phase factors and global formulation of gauge fields. Phys. Rev. D 12, 3845.
- Montonen, C.; Olive, D. (1977): Magnetic monopoles as gauge particles? Phys. Lett. B 72, 117.
- 't Hooft, G. (1979): A property of electric and magnetic flux in non-abelian gauge theories. Nucl. Phys. B 153, 141.
- Amaldi, U.; de Boer, W.; Fürstenau, H. (1991): Comparison of grand unified theories with electroweak and strong coupling constants measured at LEP. Phys. Lett. B 260, 447.
- Wilczek, F. (1999/2000): Mass without mass I und II. Physics Today 52(11), 11 und 53(1), 13.
