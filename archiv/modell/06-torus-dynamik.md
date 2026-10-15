# 06 Torus-Dynamik

Ein Torus ist die Schließungslogik mehrerer gekoppelter Zyklen ohne gemeinsame Zählvariable. Eine Kraft ist ein Vektorfeld auf dem Zustandsraum, eine Kopplung eine Rate mit Phasenlage, Dichte eine Kontaktfrequenz.

## Grundgleichungen

```
𝓜 ≅ Tⁿ × S,     θ̇ = ω(Z) + K(θ,Z),     Ż = G(θ,Z) − D(Z)·Z
```

θ sind die Zyklen (jeder mit eigener Zeitkonstante und Phase), Z die langsame Struktur (Geometrie, Temperatur, Material, Dichte), ω die Eigenfrequenzen, K die Querkopplung, G die Einspeisung, D die Dissipation. Jede Wechselwirkung erscheint als Beitrag zu ω, K, G oder D. Der konservative Kern hat fast erhaltene Aktionen I_i (θ̇ = ∂H/∂I, İ = −∂H/∂θ).

| Kanal | Wirkung | Formel |
|---|---|---|
| A Phase → Phase | ein Zyklus moduliert die Geschwindigkeit eines anderen | θ̇_i = ω_i + Σ_j k_ij sin(θ_j − θ_i − φ_ij) |
| B Medium → Phase | die Struktur taktet die Zyklen | ω_i(Z), k_ij(Z) |
| C Phase → Medium | schnelle Zyklen pumpen die Struktur | Ż = G(θ,Z) − D(Z)Z |
| D Dissipation | löscht Richtungen, bildet Attraktoren | stabilisiert nur quer zum Torus |

## Einrasten

Für zwei Phasen mit ψ = θ₂ − θ₁ gilt ψ̇ = Δω − 2k cos φ sin ψ.

- Die Zyklen rasten genau dann ein, wenn |Δω| ≤ 2k|cos φ|. Dann laufen beide mit einer Frequenz, und der Torus zieht sich auf einen Kreis zusammen (T² → S¹).
- Wächst der Verzögerungswinkel über |φ| = π/2, wird aus anziehender Rückkopplung abstoßende, ohne dass sich ihr Betrag ändert.
- An einer Resonanz m:n ist die Einrastzunge |nω₂ − mω₁| ≤ (m + n)k_mn; ihre Breite wächst linear mit der Kopplung.

**Kopplungsfenster:** k_min < k_ij(ρ,T) < k_Rast. Dichte taktet die Begegnungen (steigt und sättigt), Temperatur verschmiert die Phasen (Maximum bei mittlerem T). Zu starke Kopplung verriegelt die Phasen, der Torus kollabiert; zu schwache lässt ihn zerfasern. Bei Stoffen ist fest verriegelt (das Gitter als Identifikation), flüssig im Fenster, gasförmig zerfasert: Schmelzen gibt die Identifikation des Gitters auf, nicht den Zusammenhalt. Die beiden Ränder setzen verschiedene Kopplungen; über alle Elemente streut T_Siede/T_Schmelz von 1,02 bis 8,2 und hängt an der Ionisierungsenergie.

## Stabilität

```
λ_∥ ≈ 0 (entlang der Phasen neutral),   λ_⊥ < 0 (quer anziehend),   m·ω̄ ≠ 0 für alle m ≠ 0 (kein Einrasten)
```

Die Zahl unabhängiger Zyklen ist die Zahl der neutralen Lyapunov-Exponenten: keiner (Fixpunkt), einer (Grenzzyklus), n (n-Torus); ein positiver bedeutet Chaos. Stabilität ist Bewegung ohne Wiederholungspflicht: Die Rückwirkung muss greifen, bevor sich eine Störung aufschaukelt (τ_θ < τ_s ≤ τ_Auf).

Durchgerechnet (T² × S, Dissipation D = 0,5): Bei inkommensurablen Frequenzen (ω₂/ω₁ = φ) gibt es zwei neutrale Exponenten, λ_⊥ = −0,5, und der Umlaufschnitt füllt den Kreis dicht. Nahe 1:1 oder 3:2 rastet das Paar ein; ein Exponent wird −0,098, genau wie die Formel vorhersagt. Bei zu großer Kopplung rasten auch inkommensurable Frequenzen ein.

**Drei Kippmechanismen** bei unveränderten Kopplungen: (1) die Verzögerung überschreitet die Aufschaukelzeit (lange Integrationszeiten wie Wärme oder Trägheit), (2) Zwischenspeicher drehen die Phase, bis aus Gegenkopplung Mitkopplung wird, (3) die Rückkopplung wirkt auf einer anderen Skala als die Störung (lokal stabil, global instabil). Das System verliert dabei keinen Attraktor, es verliert die Fähigkeit, ihn rechtzeitig zu erreichen; auch langsame Dissipation kann destabilisieren. Im Dreieck Energie–Struktur–Zustand setzt ein Aufschaukeln statisch bei Σq ≥ 1 ein, oszillatorisch erst bei |q₃| > (Σt_i)(Σ1/t_i) − 1 ≥ 8: Zeitskalentrennung schützt.

**Zeitverengung:** Steigen die Kopplungsraten, fallen viele Relaxationszeiten auf eine zusammen, {τ₁, τ₂, …} → τ*. Der Torus wird schmal (verliert Dicke), nicht flach (die Dimension bleibt).

**Attraktorlandschaft:** Mehrere Tori bilden eine getrennte Vereinigung T²₁ ∪ T²₂ ∪ … mit Einzugsgebieten und Separatrizen; gekoppelte Tori bilden T² × T² = T⁴. Wechsel gibt es nur über Parameteränderung, Rauschen oder Verzweigung (Torus → Chaos oder Fixpunkt → Torus); bei fraktalen Grenzen genügen kleinste Änderungen. Die Attraktorzeit ist die Verweildauer nahe derselben Stabilitätsform.

## Kopplungsklassen

| Klasse | Frage | Richtung |
|---|---|---|
| A Energie | Wie viel fließt? | tangential in Schwingung, transversal als Verlust |
| B Phase | Wer moduliert wessen Zeit? | tangential; entscheidet Torus oder Kreis |
| C Struktur | Wer ändert die Randbedingungen? | parametrisch; formt den Torus |
| D Dissipation | Was wird gelöscht? | transversal |

Elektromagnetik: A + B, sie treibt Phasen. Mechanik: A + B + C, begrenzt und koppelt. Wärme: B + C + D, integriert über die Zeit und löscht. Strömung: A + B + D, koppelt mit. Akustik und Rauschen: B + D, Verlustpfad bzw. Phasengrenze. Chemie: C, irreversible Verformung des Torus. Gravitation: C, langsame Parametrisierung. Gedächtnis: C + D, Verengung des Zustandsraums. Abstrahlung ist der Grenzfall zwischen tangential und transversal. Viele Frequenzlinien bedeuten nicht viele Zyklen: Zwei verriegelte Frequenzen sind ein Zyklus.

## Dichte, Elektromagnetik und Gravitation

Dichte taktet nicht die Kraft, sondern die Rate, mit der Kopplungen wirksam werden. Eine Wechselwirkung setzt sich durch, wenn ihre Reaktionszeit kürzer ist als die Zeit zwischen Begegnungen. Elektromagnetik ist schnell, stark, abschirmbar und neutralisiert sich (|Q_el| ≤ Σ|Q_el,a|). Gravitation ist langsam, nicht abschirmbar und addiert sich zwangsweise (M = Σ m_a). Sie gewinnt nicht durch Stärke, sondern durch Unvermeidlichkeit, und versklavt die schnellen Freiheitsgrade: y ≃ y*(Z). Freiheitsgrade verschwinden nicht, sie verlieren ihre Unabhängigkeit.

## Zeit und Gedächtnis

- **Zeit** = Persistenz + Verzögerung + Rückkopplung. Die Zahl der Zeitrichtungen ist die Zahl der neutralen Phasen. Mehrere nicht synchronisierte Ebenen lassen sich nicht auf eine Achse projizieren: Zeit wird gefaltet, nicht addiert. Gegenwart ist der Moment, in dem die Rückkopplung rechtzeitig schließt.
- **Gedächtnis** ist die Topologie verlorener Möglichkeiten: Gedächtnis_[t₁,t₂] = 𝓡_{t₁} \ 𝓡_{t₂}. Hysterese ändert die Regel selbst, f_{t+1} = 𝒯(f_t, x_t).
- **Gedächtniskern:** Aus Kanal C und D folgt x_s(t) = ∫ 𝒦(t,t′) G dt′ mit 𝒦 = exp(−∫D du). Über ω(x_s) hängt die Phasengeschwindigkeit von der ganzen durchlaufenen Geschichte ab.

## Der Kosmos als Feldsystem

Zyklus: Differenz → Integration → Auflösung → Re-Integration. Kompression wirkt wie Masse, Metrikverzerrung wie Gravitation, Diffusion wie Entropie.

```
K_{t+1}(i,j) = K_t(i,j)·exp(r_C(C_i + C_j) − r_D(D_i + D_j)),     g = −log K,     δ log K ≈ (r_C − r_D)(C_i + C_j)
```

Nutzung verdichtet, Nähe steigt, Wege stabilisieren sich: Der Raum deformiert sich durch Nutzung. Überverdichtung wird durch Diffusion zurückgenommen. Metastabil ist das Band r_D ≈ r_C (Verdichtungsrate r_C, Diffusionsrate r_D); darunter läuft die Krümmung in die Singularität, darüber wird sie geglättet. Der Ω-Operator nimmt einer Region das Flussmonopol, wenn Persistenz und Zentrierung gemeinsam wachsen. Gedächtnis ist eingefrorene Dissipation: Das System kehrt leichter dorthin zurück, wo es sich schon einmal halten konnte.

## Wechselkräfte, Bistabilität, Hysterese

- **Kac-Budget:** Auf einem geschlossenen Zustandsraum gilt für jede Zerlegung in Bereiche Σ_i 1/⟨T_{A_i}⟩ = 1 (mittlere Rückkehrzeit ⟨T_A⟩ = μ(X)/μ(A)). Schließung verteilt ein festes Budget an Wiederkehr; ein offener Raum setzt es auf null.
- **Nur der dissipative Anteil rastet ein:** Bei Kopplung K mit Phase φ rastet K_d = K cos φ ein und kann Amplituden auslöschen (K_d > μ und |Δω| > 2√(μ(2K_d − μ))); der reaktive Anteil verschiebt nur Frequenzen.
- **Form des Torus:** Amplitude A*² = [μ − K_d(1 − ⟨cos ψ⟩)]/β, Streuung Var A = D/(4βA*²). Vor dem Einrasten wird er dünn: Das Amplitudenverhältnis r = (Δω − √(Δω² − 4K_d²))/2K_d geht gegen 1, die Schwebung √(Δω² − 4K_d²) gegen null. Er zerfällt auch, wenn die Verstärkung μ gegen null geht.
- **Welche Zyklen nebeneinander bestehen:** Zwei Zyklen an gemeinsamer Verstärkung (İ_i = 2I_i(α_i − β_iI_i − θ_ijI_j)) koexistieren genau für θ₁₂/β₂ < α₁/α₂ < β₁/θ₂₁. Parametrische Größen verschieben α und θ und schalten so Zyklen ein und aus.
- **Die Torusdimension ist begrenzt:** Ein gemeinsames Medium mit Paarkopplung k trägt höchstens N_max ≈ 1 + 2W/(πk) freie Zyklen im Frequenzband W; darüber rasten sie ein.
- **Verzweigungen sind Rollenwechsel:** Hopf und Neimark-Sacker machen eine transversale Richtung tangential (Geburt eines Zyklus bzw. Torus), Einrasten macht eine tangentiale Richtung wieder transversal, Chaos macht eine transversale dehnend.
- **Kippstelle = Divergenz der Rückkopplungssumme:** Bei q = 1 bzw. det(1 − G) = 0 divergiert Σqⁿ; vor der Falte gilt 1 − q = 2(3a)^{1/4}δ^{1/2}, der Sprung wächst mit a^{1/2}. Dort öffnet das Hysteresefenster mit dem Abstand hoch 3/2, und die Rückkehrrate fällt mit dessen Wurzel (kritische Verlangsamung). Im Kern ist die Spaltbarkeit x = E_C/(2E_S) diese Verstärkung: Die Rückstellung einer Verformung ist ∝ 1 − x, die Antwort ∝ 1/(1 − x), die Kippstelle liegt bei x = 1 (größter gemessener Wert 0,93); die Bindung je Nukleon ist am größten bei x = ¼ (E_C = E_S/2).
- **Hysterese ohne Speicher:** Ein Preisach-Zustand ist nur die Folge der nicht ausgelöschten Extremwerte; wenige Zahlen genügen für eine lange Geschichte. Der Gedächtniskern Σg_i e^{−h_iΔt} ist sein geglätteter Kleinsignal-Grenzfall: Er gilt für aktive Barrieren und kleine Amplituden (unter kT/|∂E_b|), Preisach für eingefrorene Barrieren (E_b > kT·ln(f₀t)); Preisach ist ratenunabhängig, der Kern nicht.
- **Hysterese ohne Bistabilität:** Nahe einem verschwundenen Sattel-Knoten läuft der Fluss langsam („Geisterbahn“). Beim Durchfahren verschiebt sich das Rastfenster je nach Richtung (um 0,0125–0,015 im gerechneten Beispiel), obwohl statisch kein zweiter Attraktor existiert: Pfadabhängigkeit allein aus der Langsamkeit des Flusses.
- **Rauschen und Dissipation sind eine Größe:** Beide stammen aus derselben Dämpfung γ. Die Dicke des Torus-Bandes ist ⟨δρ_⊥²⟩ ≃ D_⊥/|λ_⊥| mit D_⊥ ∝ γk_BT, die Phasendiffusion D_θ ∝ γk_BT/E*. Für λ_⊥ → 0 wird das Band beliebig dick.
- **Der kritische Punkt als Attraktor der Gegenkopplung:** Hat eine Regel am kritischen Punkt unendliche Steigung (Exponent 2β < 1), strebt q → −∞; der Arbeitspunkt rastet dort ein, und Schwankungen sinken um 1 − q (σ_T² = kT²/(C(1 − q))).
- **Rauschen zwischen zwei Zuständen:** Die Verweilzeiten folgen Kramers; auf einem geschlossenen, stationären System gilt Σ1/⟨τ⟩ = 1. Das ist ein direkter Test auf Geschlossenheit.

**Laborprüfungen:** Sie prüfen die Form (Sattel-Knoten, Wurzelgesetz, Kac-Budget), die das Modell mit der gewöhnlichen nichtlinearen Dynamik teilt; eine Laborvorhersage, die nur das Modell macht, gibt es nicht: Im Labor ist s ≈ 3·10⁵, und der galaktische Anteil ist über das Labor gleichförmig (Gezeit über 1 km 10⁻²⁶ m/s²). Die Form der Rückkehr lässt sich nachbauen: Eine RC-Kette (Kabelgleichung) mit Rückführung des Abgriffs n_R auf den Eingang, V(0) = V_ein + V(n_R), ist die Erneuerung u(0) = g_N + u(ϖ_R). Sie gibt die Bose-Form exakt, ν = 1/(1 − e^{−s}) mit s = n_R·arcosh(1 + RG/2) (R = 1 kΩ, n_R = 20, RG = 0,01: ν = 1,1567), im Tiefbereich die Wurzel (Steigung −0,499), und ein offenes Kettenende bei n_Λ = (π/2)·λ die Abweichung 1/tanh(π/2) − 1 = 9,0 %; Bauteile auf 1 % genügen, der Kern ist bei etwa 4 Hz messbar. Das prüft die Gleichung, nicht die Gravitation.

| Aufbau | Vorhersage |
|---|---|
| VO₂-Element mit Lastwiderstand | bistabil genau für R_L < max(−dV/dI); Fensterbreite ∝ ε^{3/2}, Sprung ∝ ε^{1/2}, Rückkehrrate ∝ Abstand^{1/2} |
| VO₂ auf einer Mikrobrücke im Magnetfeld | Hopf bei 3,55 V (3,26 MHz); bei 3,47 V: ohne Feld ein Zyklus, bei 1–2 T ein Torus, bei 4 T eingerastet |
| zwei selbsterregte Oszillatoren mit Phasenschieber | Rastbereich \|Δω\| ≤ 2K cos φ, bei φ = π/2 kein Einrasten; Amplitudentod bei K_d > μ |
| N Resonatoren auf gemeinsamem Substrat | erstes Einrasten bei (N − 1)k ≈ 0,57–0,63 W; nie mehr als 1 + 2W/(πk) freie Zyklen |
| eingespannter Kupferdraht im Magnetfeld | thermomechanische Bistabilität ab 5–6 mV, Vorzeichen je nach Resonanzflanke. Ist die Ringzeit der Schwingung länger als die Wärmezeit (Wärmezeit/Ringzeit ≤ 0,74, Güte ab 3000) und der Antrieb stark genug, endet der heiße Ast 0,0010–0,6 Hz vor seiner Falte in einer unterkritischen Hopf-Verzweigung (Ringfrequenz 0,20–1,5 Hz) und springt auf den kalten Ast; dauerhaftes Selbstpulsieren gibt es nicht |
| Ferritkern im Induktionsfeld nahe der Curie-Temperatur | der Arbeitspunkt rastet knapp unter T_C ein, T_C − T* ∝ P₀^{−1/(2β)} |
| thermisch bistabiler Resonator mit langsamer Alterung | die Erholungsrate nach kleinen Stößen fällt vor dem Sprung wie √(t_F − t); heilt die Alterung schneller aus, bleibt der Sprung aus |
