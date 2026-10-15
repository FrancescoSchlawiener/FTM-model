# 06 Torus-Dynamik

Ein Torus ist die Schließungslogik mehrerer gekoppelter Zyklen. Eine Kraft ist ein Vektorfeld auf dem Zustandsraum, eine Kopplung eine Rate mit Phasenlage.

## Grundgleichungen

```
θ̇ = ω(Z) + K(θ, Z),     Ż = G(θ, Z) − D(Z)·Z
```

θ sind die Zyklen, Z die langsame Struktur, ω die Eigenfrequenzen, K die Querkopplung, G die Einspeisung, D die Dissipation. Jede Wechselwirkung erscheint als Beitrag zu ω, K, G oder D.

| Klasse | Frage | Rolle |
|---|---|---|
| A Energie | Wie viel fließt? | tangential in Schwingung, transversal als Verlust |
| B Phase | Wer moduliert wessen Zeit? | tangential |
| C Struktur | Wer ändert die Randbedingungen? | parametrisch |
| D Dissipation | Was wird gelöscht? | transversal |

## Einrasten

Für zwei Phasen mit ψ = θ₂ − θ₁ gilt ψ̇ = Δω − 2k·cos φ·sin ψ. Sie rasten ein, wenn |Δω| ≤ 2k|cos φ|. Dann zieht sich der Torus auf einen Kreis zusammen (T² → T¹). Nur der dissipative Anteil k·cos φ rastet ein. Der reaktive verschiebt nur Frequenzen. Vor dem Einrasten wird der Torus dünn: Die Phasenordnung steigt.

T² → T¹ geht auf drei Wegen: Einrasten (beide Zyklen schwingen weiter, gebunden), Knicken, Verlust eines Zyklus.

**Kopplungsfenster:** k_min < k < k_Rast. Zu starke Kopplung verriegelt die Phasen, zu schwache lässt den Torus zerfasern. Ein gemeinsames Medium mit Paarkopplung k trägt im Frequenzband W höchstens N ≈ 1 + 2W/(πk) freie Zyklen. Darüber rasten sie ein.

## Stabilität

Die Zahl unabhängiger Zyklen ist die Zahl der neutralen Lyapunov-Exponenten. Quer zum Torus sind sie negativ. Stabilität ist Bewegung ohne Wiederholungspflicht: Die Rückwirkung muss greifen, bevor sich eine Störung aufschaukelt. Verzweigungen sind Rollenwechsel: Die Geburt eines Zyklus macht eine transversale Richtung tangential, Einrasten eine tangentiale transversal.

**Kippstelle:** Bei q = 1 divergiert die Rückkopplungssumme Σqⁿ. Vor der Falte gilt 1 − q ∝ δ^{1/2}. Das Hysteresefenster öffnet mit dem Abstand hoch 3/2, die Rückkehrrate fällt mit dessen Wurzel (kritische Verlangsamung).

**Hysterese ohne Speicher:** Ein Preisach-Zustand ist nur die Folge der nicht ausgelöschten Extremwerte. Rauschen und Dissipation stammen aus derselben Dämpfung. Die Dicke des Torusbandes ist D_⊥/|λ_⊥|.

## Kac-Budget

Auf einem geschlossenen Zustandsraum gilt für jede Zerlegung in Bereiche

```
Σᵢ 1/⟨T_{Aᵢ}⟩ = 1,     ⟨T_A⟩ = μ(X)/μ(A)
```

Schließung verteilt ein festes Budget an Wiederkehr. Ein offener Raum setzt es auf null. Für die Verweilzeiten eines geschlossenen, stationären Systems zwischen zwei Zuständen gilt Σ1/⟨τ⟩ = 1: ein direkter Test auf Geschlossenheit.

## Labor

Die Form der Rückkehr lässt sich nachbauen: Eine RC-Kette mit Rückführung des Abgriffs n_R auf den Eingang, V(0) = V_ein + V(n_R), ist die Erneuerung u(0) = g_N + u(ϖ_R). Sie gibt die Bose-Form exakt, ν = 1/(1 − e^{−s}) mit s = n_R·arcosh(1 + RG/2). Bei R·G = 0,01 und n_R = 20 ist ν = 1,1567. Das prüft die Gleichung, nicht die Gravitation: Im Labor ist s ≈ 3·10⁵, eine Laborvorhersage, die nur das Modell macht, gibt es nicht.
