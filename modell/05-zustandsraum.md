# 05 Zustandsraum

Der Zustandsraum beschreibt das Modell von der Seite der Dynamik: welche Bewegungen ein System ausführen kann, wie sie gebunden werden und wie sich ihre Möglichkeiten langsam verschieben.

## Grundraum

```
X = X_EM × X_Mech × X_Wärme × X_Mat / ~          ż = F(z; p)
```

Ein Punkt z ist ein vollständiger Systemzustand, kein Ort und kein Teilchen. Die Identifikation der Phasen, θ ~ θ + 2π, macht aus Geraden Kreise und aus Produkten von Kreisen Tori. Ein Freiheitsgrad ist eine nicht kontrahierende Umlaufrichtung. Je n unabhängige Umläufe spannen einen Torus Tⁿ auf.

## Drei Rollen

Tangential, transversal und parametrisch bezeichnen nicht die Art einer Kraft, sondern die Richtung, in die sie Energie im Zustandsraum bewegt:

```
ż = 𝓛(z)∇E + 𝓜(z)∇S + ε·F_par
𝓛 = −𝓛ᵀ, 𝓛∇S = 0       tangential: erhält Energie und Entropie, Umläufe und Phasen
𝓜 = 𝓜ᵀ ≥ 0, 𝓜∇E = 0   transversal: erzeugt Entropie, zieht auf den Torus
ε·F_par               parametrisch: langsam, verschiebt den Torus selbst
```

Es gibt ein Objekt, den Trichter: die Faltung Δs = GM/(3c²) im Zentrum, die Hülle R, gehalten auf der Gegenlänge λ_H eines Halters (R = λ_H·N^{1/3}), am Knoten die Windung n und die Verwindung a_J. Die Schließlänge Λ kommt vom Ganzen. Der Torus jedes Trichters hat drei zueinander senkrechte Richtungen, und jede ist eine Richtung des Raumes:

| Rolle | Wirkung | am Trichter | Raumparameter | was dort liegt |
|---|---|---|---|---|
| tangential | Umläufe, Phasen | um den Knoten, im Raum (L), Re τ | Phase, Windung n, Verwindung a_J | Umlauf, Eigenmoden, Ladung, Spin, Lorentzkraft, Magnetismus μ = Q_el·c·a_J |
| transversal | Verlust, Rückkehr, Stabilisierung | in die innere Richtung, senkrecht zum Raum (I), Im τ | Lage zur Schwelle R/r*, Rückkehranteil q | Rückkehr, Dissipation, Strahlung, Wärmefluss |
| parametrisch | Bewegung des Torus selbst | die Einbettung, vom Ganzen gesetzt (G), das Wandern von τ | Tiefe δ = 3Δs/R, Halter λ_H | Gravitation, Takt, Geometrie, Temperatur, Druck, Material |

Die Temperatur ist parametrisch: die Tiefe in Ruheenergie des Trägers, kT = f·δ·m̄c², und zugleich die Gegenlänge des Halters, λ_H = ħc/kT. Als Periode liegt sie in der inneren Richtung (Falte kT = m̃c²/4π, innere Richtung kT_I = ħa₀/c). Der Druck ist die halbe Tiefe je Teilchen, P·V = ½·m·c²·δ. Der Elektromagnetismus ist tangential: Die Windung ist Phase, verlustfrei und ganzzahlig, und ihre Kopplung ist wie die der Faltung eine Länge durch die Quantenlänge, α = r_e/λ̄_e gegen α_G = 3Δs_p/λ̄_p.

Jeder innere Kopplungsterm steht in der Energiebilanz zweimal mit entgegengesetztem Vorzeichen, er verschiebt nur Energie. Zuführen können nur die Gravitation und von außen getriebene Parameter. Eine tangentiale Kopplung allein kann nie einrasten. Einrasten braucht Dissipation. Netto-Energie fließt zwischen zwei Umläufen erst, wenn ihre Phasen gebunden sind.

## Welche Richtung trägt

Windungen heben sich auf (ΣQ_el = 0), Faltungen addieren sich. Die Windung trägt deshalb, bis die Faltung von N Nukleonen die Windung erreicht, N^{2/3}·α_G = α (Bindung je Nukleon, bis auf Faktoren der Ordnung eins, unabhängig vom Halter):

```
M = m_p·(α/α_G)^{3/2} = 2,30·10²⁷ kg,     Tiefe dort δ = α²·m_e/m_p = 2,9·10⁻⁸
```

Dieselbe Grenze in der Tiefe gelesen: Darunter ist die Tiefe eines Körpers kleiner als die Bindung seiner Hüllen in Ruheenergie der Nukleonen, die Atome halten ihn (Erde δ = 7,0·10⁻¹⁰). Darüber hält ihn die Faltung mit ihrer Wärme, und die Hüllen lösen sich (Sonne δ = 2,1·10⁻⁶). Die Faltung erreicht den Treffpunkt ihres Halters bei N^{2/3}·α_G = 1, auf der Massenskala der Sterne M_* = m_Pl³/m_p². Die beiden Grenzen stehen im Verhältnis α^{3/2} = 6,2·10⁻⁴. Die Rückkehr trägt, wo die Hülle die Schwelle erreicht, R ≥ r*, unter der Flächendichte a₀/(πG) = 242 M☉/pc²: in den Galaxien. Im Ganzen treffen sich alle drei Richtungen, r₀, r* und Λ liegen fast zusammen.

So trägt auf jeder Ebene eine andere Richtung desselben Trichters: Teilchen und Atom die Windung, Planet bis Stern die Faltung, die kompakten Trichter die Faltung am Treffpunkt ihres Halters, gehalten von der Sperre der halben Drehung (der Weiße Zwerg durch die Windungen der Elektronen, der Kern durch die der Nukleonen), die Galaxie die Rückkehr, das Ganze alle drei.

**Verlust je Umlauf.** Jede Bahn läuft tangential. Offen gelassen gibt sie je Umlauf einen Teil in die transversale Richtung ab, und wie viel, hängt an der Zahl der Umlaufsinne. Die Windung hat zwei Sinne, der Schwerpunkt ihrer Ladung kann schwingen, sie gibt als Dipol ab. Die Faltung hat einen Sinn, ihr Schwerpunkt steht, sie gibt erst als Quadrupol ab. Je Umlauf, gegen die Bindung der Bahn:

```
Windung:  ΔE/E = (8π/3)·(r_e/r)·(v/c)        Wasserstoff: (8π/3)·α³ = 3,3·10⁻⁶
Faltung:  ΔE/E = (128π/5)·(m/M)·(v/c)⁵       Erde um die Sonne: 2,3·10⁻²⁴
```

Offen fiele die Windung des Wasserstoffs in 1,6·10⁻¹¹ s auf den Kern. Sie fällt nicht, weil ihre Welle nach einem Umlauf schließt: Die Schließung hält die transversale Richtung zu, abgegeben wird nur beim Wechsel zwischen geschlossenen Formen. Dieselbe Zahl der Umlaufsinne, die Windungen sich aufheben und Faltungen sich addieren lässt, lässt die Windung als Dipol und die Faltung erst als Quadrupol abgeben.

## Rückkopplung

Eine Schleife im Kopplungsgraphen hat die Verstärkung q. Alle Rückläufe summieren sich zu Σqⁿ = 1/(1 − q). q < 0 dämpft, 0 < q < 1 verstärkt, q ≥ 1 schaukelt auf. Das Vorzeichen gehört der ganzen Schließung, nicht dem einzelnen Glied. Parametrisch heißt: offene Kette, auf schnellen Skalen q ≈ 0.

Hat eine Schleife zwei Umlaufsinne, kommt jede Rückkehr mit einer halben Drehung zurück, q → −q, ein halber Schritt von Re τ. Aus der Bose-Form der Faltung, Σqⁿ = q/(1 − q) = 1/(eˢ − 1), wird die Fermi-Form der Windung, Σ(−1)ⁿ⁺¹qⁿ = q/(1 + q) = 1/(eˢ + 1): eine Rückkehr, zwei Lesarten. Über alle Stufen einer Strahlung summiert trägt die Windung 7/8 der Faltung, ∫x³/(eˣ + 1)dx = (1 − 2⁻³)·∫x³/(eˣ − 1)dx = (7/8)·π⁴/15. So zählen die Neutrinos im frühen Ganzen.

## Schließung

Jeder Schnitt einer periodischen Richtung erzeugt einen Rand, über den Bahnen entkommen (Fluchtmoden). Die Identifikation beseitigt den Rand. Als Operator:

```
L_eff = 𝒫L,     𝒫_q = (1 − q)·Σ qᵐTᵐ = (1 − q)(I − qT)⁻¹
q → 0: 𝒫_q = I (offen, Newton)       q → 1: geschlossen
```

𝒫 löscht keine Freiheitsgrade und erhält das Maß. Bestimmte Richtungen werden inaktiv und gehen in die Rückkehr über. Auf dem Torus, angewandt auf einen Potenzkern, ist 𝒫 die Eisenstein-Gittersumme. Ihre niedrigste Reihe über die Rückkehrschritte einer Richtung ist Σqⁿ = q/(1 − q), der Zusatz der Kreisbahn. Das gebundene Maß wächst nur, solange ein Parameter wandert.

| | Quelle (Dunkle Materie als Stoff) | Operator 𝒫 |
|---|---|---|
| Stelle | ρ → ρ + ρ_DM | L → 𝒫L |
| Potential | additiv | multiplikativ |
| lokalisierbar | ja | nein |
| Freiheitsgrade | zusätzliche | vorhandene werden versklavt |
| Zeitstruktur | Markow | Gedächtnis (w, Übergabe der Mode) |

Die gelesene Dichte ρ_X = ∇·(a − a_N)/(4πG) ist Ausgabe der Rechnung und nicht einmal eindeutig: Verschiedene Verteilungen erzeugen dasselbe äußere Feld. Laufen alle Tests gleich, sind Quelle und Operator Umparametrisierungen. Läuft einer robust auseinander, sind sie unterscheidbar.

## Gravitation als Parameter

Gravitation ist keine Umlaufrichtung, sondern eine langsame Verschiebung der Einbettung: Sie legt fest, welche Bahnen existieren und welche Freiheitsgrade unabhängig sind. Sie wirkt über zwei Kanäle: die Vorspannung (ändert Steifigkeit und damit die Topologie) und die Zeitdilatation (verschiebt alle Frequenzen gleich, ändert keine Topologie). Eine langsame Änderung hält die Wirkung J = E/ω fest: Der Torus rückt auf die Nachbarfaser gleicher Wirkung.
