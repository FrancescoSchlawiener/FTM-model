# 05 Zustandsraum

Der Zustandsraum beschreibt das Modell von der Seite der Dynamik: nicht, was im Raum ist, sondern welche Bewegungen ein System ausführen kann, wie sie gebunden werden und wie sich ihre Möglichkeiten langsam verschieben.

## Der Grundraum

```
X = X_EM × X_Mech × X_Wärme × X_Mat / ~          ż = F(z; p)
```

Ein Punkt z ist ein vollständiger Systemzustand zu einem Zeitpunkt (alle Felder, Temperaturen, Spannungen, Parameter), kein Ort und kein Teilchen. Das Zeichen ~ steht für die Identifikation der Phasen, θ ~ θ + 2π: Aus Geraden werden Kreise, aus Produkten von Kreisen Tori.

Das Vektorfeld F steckt in fünf Blöcken:

| Block | Gleichung | Rolle |
|---|---|---|
| Maxwell | ∇·D = ρ, ∇·B = 0, ∇×E = −∂_tB, ∇×H = J + ∂_tD | Struktur: wo Feldlinien enden, geschlossene Felder, Zeit erzeugt Umlauf |
| Stoffgesetze | D = ε(E,T,…)E, B = μ(H,T,…)H, J = σ_el(T,E,…)E | hier sitzt die Kopplung an Temperatur, Mechanik, Vorgeschichte |
| Lorentz, Poynting | f = ρE + J×B, ∂_tu + ∇·S = −J·E | Feld → Bewegung, Feld → Wärme |
| Wärme | ρc ∂_tT = ∇·(k∇T) + P_V | langsames Gedächtnis |
| Mechanik | ρü = ∇·σ + f, σ = C:ε | Formänderung der Kopplung selbst |

Daraus folgt die Trennung in schnelle Variablen und langsame Parameter: ż_f = F(z_f; p), ṗ = εG(z_f, p), ε ≪ 1. Die Multiphysik ist keine Summe von Kräften, sondern eine Organisation von Flussrichtungen.

**Geschlossener Wirkkreis:** Energie, Struktur und Zustand wirken in allen sechs Richtungen aufeinander (Energie verformt Struktur, Struktur lenkt Energie, Struktur altert zu Zustand, Zustand verändert Form, Energie hinterlässt Spuren, Zustand verändert die Antwort). Der Kreis hat keinen Anfang und keinen äußeren Auslöser.

**Freiheitsgrad = nicht kontrahierende Umlaufrichtung.** n unabhängige Umläufe spannen einen Torus Tⁿ auf. Mit einer elektromagnetischen und einer mechanischen Phase ist der natürliche Attraktor ein Torus T², kein Kreis.

## Drei Rollen

Tangential, transversal und parametrisch bezeichnen nicht die Art einer Kraft, sondern ihre Rolle im Zustandsraum.

| Rolle | Gleichung | Wirkung | Träger |
|---|---|---|---|
| tangential | θ̇ = ω(p) | Umläufe und Phasen | Elektromagnetik, Eigenmoden, Bahnen |
| transversal | ρ̇ = −κ(p)ρ | zieht auf den Torus, stabilisiert | Dissipation, Sättigung, Strahlung |
| parametrisch | ṗ = εG | bewegt den Torus selbst (Faserbündel π: M → P) | Temperatur, Geometrie, Material, Gravitation |

**Normalform:** θ̇ = ω(p), ρ̇ = −κ(p)ρ, ṗ = εG(θ,ρ,p). Die Dynamik lebt auf einem Torus, der durch Dissipation stabilisiert und durch Struktur langsam verformt wird.

**Die Rollen als Gleichung (GENERIC):**

```
ż = 𝓛(z)∇E  +  𝓜(z)∇S + F_Tor  +  ε·F_par
𝓛 = −𝓛ᵀ, 𝓛∇S = 0 (tangential: erhält E und S, divergenzfrei)      𝓜 = 𝓜ᵀ ≥ 0, 𝓜∇E = 0 (transversal: erzeugt S)
```

Symmetrische Kreuzkopplungen (Onsager, z. B. Seebeck/Peltier mit Π = S_b·T) sind dissipativ und transversal; antisymmetrische (Casimir: Lorentz/Induktion, Thermoelastik) verschieben verlustfrei und sind tangential. Weil der tangentiale Teil divergenzfrei ist, gilt Σλ_k = ⟨div F_trans⟩: Eine tangentiale Kopplung allein kann nie einrasten; Einrasten braucht Dissipation. Netto-Energie fließt zwischen zwei Umläufen erst, wenn ihre Phasen gebunden sind.

**Rückkopplungsschleifen:** Eine Schleife im Kopplungsgraphen hat die Schleifenverstärkung q (Produkt der Übertragungsfaktoren); alle Rückläufe summieren sich zu Σ qⁿ = 1/(1 − q), im Netz (1 − A)⁻¹ für Spektralradius < 1. q < 0 dämpft, 0 < q < 1 verstärkt konvergent, q ≥ 1 schaukelt auf (statisch Durchgehen, dynamisch Geburt eines Zyklus). Transversale Kontraktion ist die Sättigung der Schleifenverstärkung am Punkt q = 1: Die Kontraktionsrate ist λ_⊥ ∝ ∂_A ln q < 0 bei der Amplitude A*, an der q = 1 wird. Das Vorzeichen gehört der ganzen Schließung, nicht dem einzelnen Glied: Bei steigendem Widerstand R(T) geht dieselbe Kette stromgetrieben durch (q > 0) und stabilisiert selbsterregt (q = −0,05); bei fallendem R(T) kehrt sich das um. Die Verstärkung ist zugleich Verlangsamung: δz = F_p·δp/(1 − q), τ_rel = τ₀/(1 − q). Parametrisch heißt: offene Kette, auf schnellen Skalen q ≈ 0 (Gravitation → Vorspannung → Steifigkeit hat keinen Rückweg). In der Energiebilanz der gekoppelten Felder steht jeder innere Kopplungsterm zweimal mit entgegengesetztem Vorzeichen; Energie zuführen können nur die Gravitation und von außen getriebene Parameter. Dieselbe Summe beschreibt die Anhäufung bei wiederholter Gabe, 1/(1 − e^{−kτ}), und einen selbsterregten Ereignisprozess mit Verzweigungsgrad q: Die Zusatzwirkung ist die Wiederholung des Ausgangs, keine zweite Sorte.

**Schließung:** Jeder Schnitt einer periodischen Richtung erzeugt einen Rand, über den Bahnen entkommen können (Fluchtmoden). Die Identifikation θ ~ θ + 2π beseitigt den Rand und damit die Fluchtmoden. Als Operator: 𝒫 = (1/M)·Σ_m T_{mL} vertauscht mit dem Generator (𝒫𝓛 = 𝓛𝒫). Auf dem geschlossenen Kreis der Länge L bleibt ein Nullmode mit Spektrallücke 4π²D/L², ein Rand hat keinen Nullmode und die Fluchtrate π²D/L²; die Bindung skaliert wie κ_eff = c·D/L². Längs des Torus bindet die Topologie, quer die Dissipation. 𝒫 erhält das Maß (doppelt stochastisch): Die inaktiven Fluchtrichtungen werden nicht entfernt, sondern gehen in die Rückkehr über: Sie gehen als Operator je Raumpunkt in die innere Richtung über und tragen dort wie der Vakuumteil mit w = −1 bei. Λ ändern sie nicht; Λ setzt allein die Faltung mit ihren Rückläufen. Das gebundene Maß wächst nur, solange ein Parameter wandert (Bindungsstärke oder Schleifenlänge): Der Parameter bewirkt die Bindung, nicht umgekehrt.

**Schließung in der Galaxie:** Winkel, Drehungen und Umlaufphase sind auf jeder gebundenen Bahn schon identifiziert; sie entfernen Fluchtmoden des Winkels, nicht die radiale Flucht, und binden nicht zusätzlich (Drehbilder eines Balkens lassen den Monopol exakt gleich). Bindend ist die Rückkehr: die Verschiebung T um einen Rücklauf in der inneren Richtung, kausal summiert mit dem Rückkehranteil q = e^{−√(a_N/a₀)}.

```
𝒫_q = (1 − q)·Σ_{m≥0} qᵐTᵐ = (1 − q)(I − qT)⁻¹        q → 0: 𝒫_q = I (offen, Newton)     q → 1: 𝒫_q → 𝒫 (geschlossen)
```

Die Rückkopplungssumme (I − qT)⁻¹ ist die Resolvente dieser Translation: Sie verstärkt genau den wiederkehrenden Anteil um 1/(1 − q), Fluchtmoden bleiben beschränkt (quasiperiodische Rückkehr: invariante Mode 10⁴-fach bei q = 0,9999, Fluchtmoden höchstens 567-fach). Damit ist die Galaxie ein Fall der Schließung durch Translation; Nahfeld heißt Rückkehr offen, Fernfeld Rückkehr geschlossen (nicht die Topologie des Raumes). Ein senkrechter Projektor ist 𝒫_q nur an den Enden: Dazwischen liegen seine Eigenwerte auf dem Kreis um 1/(1 + q) mit Radius q/(1 + q), am weitesten von senkrecht bei q = ⅓.

**Gravitation** ist keine Umlaufrichtung, sondern eine langsame Verschiebung der Einbettung: Sie bewegt das System von Torus zu Torus. Sie legt fest, welche Bahnen existieren, welche Freiheitsgrade unabhängig sind und welche Bahnen geschlossen bleiben (bei rationalem Frequenzverhältnis T² → T¹). Sie ist nicht abschirmbar; ihre Regel ist nicht additiv pro Teilchen, die Faltung ist es. Das gilt, solange ihr Anteil an der Rückstellung klein ist (f_g ≪ 1); trägt sie die Rückstellung allein (f_g → 1, Pendel), wird sie tangential, und bei f_g = ½ kann ein stehender Körper knicken.

Sie wirkt über zwei Kanäle:

```
ω_m,beob(r) = √(k_el/m ± 3c²Δs/(ℓ·r²)) · √(1 − 6Δs/r)
```

Der erste Faktor ist der differentielle Kanal (Vorspannung → Steifigkeit → Resonanz); er ändert die Topologie. Der zweite ist die Zeitdilatation; sie verschiebt alle Frequenzen gleich und parametrisiert nur die Zeit um. Bei 1000 m Höhenunterschied ist der differentielle Kanal 1,44·10⁹-mal stärker; erst am Neutronenstern werden beide vergleichbar. Eine langsame Änderung der Gravitation leistet die parametrische Arbeit P = ½·m·ġ·x²/ℓ und hält die Wirkung J = E/ω fest: Der Torus rückt auf die Nachbarfaser gleicher Wirkung.

T² → T¹ geht auf drei Wegen: Einrasten (beide Zyklen schwingen weiter, gebunden), Knicken und Verlust eines Zyklus (inverse Hopf-Verzweigung). Vor dem Einrasten wird der Torus „dünn“: Die Phasenordnung |⟨e^{iψ}⟩| steigt.

## Der Bindungsoperator

Bindung ohne zusätzliche Quelle: Der Operator schränkt die zugänglichen Bahnen ein, statt Masse hinzuzufügen.

```
L_eff = 𝒫L,    ż = (I − Π⊥)F(z)         Π⊥ projiziert auf Flucht- und Instabilitätsrichtungen
Φ_eff = Φ_bar·(1 + f(ℓ/λ_S))            multiplikativ; ℓ Bahnskala, λ_S Schließungsskala aus Geschichte und Dissipation
ẍ = −GM/x² − γ ∫ K(t − t′) ẋ(t′) dt′     Gedächtniskern
```

- 𝒫 löscht keine Freiheitsgrade; der Zustand behält seine volle Dimension. Bestimmte Richtungen werden dynamisch inaktiv, der zugängliche Bereich schrumpft.
- Spielbeispiel, nicht der Kern der Regel: Mit dem Maxwell-Kern K(Δt) = e^{−Δt/t_K} wird die Gleichung lokal in einer zusätzlichen Koordinate y: ẋ = v, v̇ = −GM/x² − γy, ẏ = v − y/t_K. Die Vergangenheit steckt in y, einer Koordinate des jetzigen Zustands. Für H = v²/2 − GM/x + γy²/2 gilt Ḣ = −γy²/t_K.
- Im Spielbeispiel bindet unendliches Gedächtnis (t_K → ∞) ohne Energieverlust: Ein Teilchen mit v = 1,6 > v_esc = √2 kehrt bei γ = 0,01 an x = 9,75 um, statt zu entkommen. Der Effekt skaliert mit der Geschichte, nicht mit der Masse.
- Der Spielkern ist ein Maxwell-Element (Feder γ, Dämpfer γt_K): Der Speicheranteil γ(ωt_K)²/(1 + (ωt_K)²) bindet, der Verlustanteil γωt_K/(1 + (ωt_K)²) reibt. Bindendes Gedächtnis ist Struktur, alterndes ist Zustand.
- Mit Gewicht w (1 in mitatmenden Bereichen, → 0 in gebundenen): Eliminiert man die projizierte Richtung, bleibt für die Baryonen ein Gedächtniskern (Zwanzig-Mori), die Green-Funktion des druckfreien Einfalls, in der Zeit der Bahn; er wächst und klingt nicht ab (Kern der Mode, nicht der Regel).

| | Quelle (Dunkle Materie) | Operator 𝒫 |
|---|---|---|
| Stelle | ρ → ρ + ρ_DM | L → 𝒫L |
| Potential | Φ_bar + Φ_DM (additiv) | Φ_bar(1 + f) (multiplikativ) |
| lokalisierbar | ja | nein |
| Freiheitsgrade | zusätzliche | vorhandene werden versklavt |
| Zeitstruktur | Markow | Gedächtnis |

**Poisson-Inversion:** Jede gemessene Zusatzbeschleunigung lässt sich formal in eine Dichte umrechnen, ρ_X = ∇·(a_obs − a_bar)/(4πG); für a = v₀²/r ist das ein isothermer Halo. Diese Dichte ist Ausgabe der Rechnung, nicht Eingabe der Dynamik: Man beobachtet globale Bindung und projiziert sie auf lokale Inventarparameter. Eindeutig ist diese Dichte nicht einmal: Nach dem Satz der Äquivalenzschicht erzeugen verschiedene Verteilungen dasselbe äußere Feld. Ein Laborbild ist die Paul-Falle, deren zeitgemitteltes Einschließen wie eine Ladung gelesen würde, die nicht da ist.

## Modulformen als Werkzeug

Ein flacher Torus ist ℂ/(ℤ + τℤ). Modulformen beschreiben, wie er verformt werden kann, ohne dass seine Topologie reißt. Die Spitze q = e^{2πiτ} → 0 bedeutet, dass ein Umlauf gegen den anderen verschwindet (T² → T¹); dort verschwindet eine Spitzenform.

```
D_k = q d/dq − (k/12)E₂                      bildet Modulformen vom Gewicht k auf Gewicht k+2 ab
q dE₂/dq = (E₂² − E₄)/12,   q dE₄/dq = (E₂E₄ − E₆)/3,   q dE₆/dq = (E₂E₆ − E₄²)/2
Δ = (E₄³ − E₆²)/1728 = q∏(1 − qⁿ)²⁴ = q − 24q² + 252q³ − …,     D₁₂Δ = 0
```

Weitere Entsprechungen:

- **Gittersumme:** 𝒫 auf dem Torus, angewandt auf einen Potenzkern z^{−k}, ist die Eisenstein-Reihe G_k(τ) = Σ′(m + nτ)^{−k}. Summiert man zuerst längs des Umlaufs, wird jeder Rückkehrschritt n zur Reihe (−2πi)^k/(k−1)!·Σ_d d^{k−1}q^{nd}. Die niedrigste Reihe (d = 1) über alle Schritte einer Richtung (n ≥ 1) ist Σqⁿ = q/(1 − q), der Zusatz der Kreisbahn; die volle Gittersumme mit ±n gibt das Doppelte. Die Modulgruppe vertauscht Umlauf und Rückkehr, die Verstärkung ν aber nicht (τ = i/2 und 2i geben 1,045 und 1,0000035): Getragen wird ν vom markierten Torus ℍ/⟨τ ↦ τ + 1⟩, nicht vom Modulraum, dessen Spitzenbereich erst bei a_N = 4π²a₀ beginnt.
- **Lambert-Reihe:** Σσ₁(n)qⁿ = (1 − E₂)/24 ist die Rückkopplung über alle Harmonischen. Für die Beschleunigung trägt nur die Grundschwingung: Bei q = e^{−s} wächst die Reihe wie π²/(6s²), und alle Harmonischen gleich gewichtet gäben a → (π²/6)·a₀ ohne Tully-Fisher.
- **Ramanujan-Summen:** c_k(n) = Σ_{(a,k)=1}e^{2πian/k} bilden die Projektoren Π_k = (1/N)Σc_k(n)Tⁿ auf genaue Periode k, mit ΣΠ_k = I und Π₁ = 𝒫. Eine Bahn mit Re τ_M = p/k liegt ganz in Π_k; so lassen sich Resonanzen nach ihrer Ordnung trennen.
- **Herkunft der Zahlen:** 24, 240, 504 = (2π)^k/(Γ(k)ζ(k)) für k = 2, 4, 6 sind die Normierungen von E₂, E₄, E₆; 12 = 2/(1 − ½ − ⅓) kommt aus den beiden Fixpunkten der Modulgruppe (Ordnung 2 bei τ = i, Ordnung 3 bei τ = e^{2πi/3}); 1728 = 3·240 + 2·504. Der Vakuumteil 4π⁴/15 = 4(2π)⁴/240 verlangt drei Raumdimensionen (Zustandsdichte x²).
- **Bilaterale Summe 1ψ1:** Σ_{n∈ℤ}(a;q)_n/(b;q)_n zⁿ = (q, b/a, az, q/(az); q)_∞/(b, q/a, z, b/(az); q)_∞ ist die q-Form der Wahl zwischen beiden Umlaufsinnen und einem. Für b = q enthält jedes Glied mit n < 0 den Faktor 1 − q⁰ = 0; die Summe endet von selbst beim Ort, und mit a = q bleibt Σ_{n≥0}zⁿ. Welche Größe b physikalisch ist, ist offen.
- **Rogers-Ramanujan:** Der Kettenbruch führt für q → 1 auf die goldene Wicklung. Die Ausschlussregel „Abstand ≥ 2“ senkt die Zahl der Zustände von 2 auf φ je Platz; allgemein S ≃ 2π√(c_eff·n/6), c_eff = 1 − 3/(2k+1).
- **Mock-Theta mit Schatten:** geschlossene Dynamik unter permanenter Dissipation (tangential + transversal).
- **Hecke-Operatoren:** wechseln zwischen Skalenebenen, Σσ₁(pm)q^m = (1+p)ℓ(q) − p·ℓ(q^p). Der Übertrag zwischen den Skalen ist unitär (|U_k| ≤ k+1): er wächst nicht und verlischt nicht.
- **Grundverschiebung q^{−1/24}:** die Nullpunktenergie einer Schleife, E₀ = −πħc/(6ℓ_S) bei der Schleifenlänge ℓ_S; 24 ist die Multiplikatorordnung in Δ = η_D²⁴ (Dedekind).
- **Master-Theorem:** löst das Gedächtnisintegral geschlossen.
- **Siegel-Form χ₁₀ (Sp(4,ℤ)):** Ein präzedierender Bahntorus braucht zwei Modulparameter; die Rotation geht quadratisch und ohne Drehsinn in den Kern ein und verschwindet bei ganzzahliger Knotenzahl.

Die Ableitung führt nicht aus dem System heraus („Rechnen ohne Rechnen“). Die Zerlegung q df/dq = D_k f (erlaubt) + (k/12)E₂f (symmetriebrechend) entspricht F = (I − Π⊥)F + Π⊥F: Die Modulformen sind die Eigenmoden von 𝒫.

## Trennbeobachtungen

| Test | Quelle erwartet | Operator erwartet |
|---|---|---|
| dynamische Masse gegen Linsenmasse, η = M_dyn/M_lens | η konstant | im Gleichgewicht η = 1 für Kugeln und außerhalb von Scheiben, in Scheiben innen η > 1 (Wirbelanteil der Regel); sonst Abweichung nur als Nachlauf der Mode ξ (Haufenstoß, frisch kollabierte Bereiche) |
| Wachstum von Struktur | nur über die Materiedichte | zusätzliche Skalen- und Epochenabhängigkeit |
| Form der scheinbaren Verteilung | universelle Profile | stärker an Umfeld und Dissipationsgeschichte gebunden |
| Dissipation und Gedächtnis in der Dynamik | keines (Markow) | Pfadabhängigkeit (w), Nachlauf nur über die Mode ξ |
| Äquivalenzprinzip | universelle Kopplung | alle Materie folgt der Regel gleich (Sterne, ruhendes Gas); Licht sieht nur ihren wirbelfreien Teil |
| Konsistenz über alle Kanäle | ein Zusatz bedient alle Kanäle | kanalspezifische Reste möglich |

Laufen alle Tests in den Daten gleich, sind beide Formen nur Umparametrisierungen voneinander; läuft eine robust auseinander, sind sie unterscheidbar.
