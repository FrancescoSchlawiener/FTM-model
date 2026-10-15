# Was logisch aus dem Modell folgt

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. F13 ist aktuell und bis auf zwei Stellen unverändert. Geändert: a₀ = cH₀/2π als Modellformel → a₀ = c²/(2Λ), heute ≈ cH₀/2π; F2 und Aussage 1: Grenze Einzeltrichter/kollektiv bei g(R) = a₀ → Gezeitendichte ρ_J, Sternsysteme kollektiv (kompakter Kugelsternhaufen ν = 1,008, sehr dünne kollektiv mit vollem Zusatz); F3: r* als geometrisches Mittel aus r₀ und der Länge des atmenden Ganzen → aus r₀ und der Schließlänge Λ (innere Länge, heute fest, Λ ≈ πR_H nur heute); F8: „globale Neigung setzt die Regel“ → Gleichheit nur heute, g_U 5,2·10⁻¹¹ → 5,11·10⁻¹¹ m/s², r*_U/R_U 0,66 → 0,65; F10: Lesart Q = 2π − 1 als Zahlennähe ohne Herleitung markiert; [G] an Freeman-Grenze, Beobachtungsaussagen in F2 und DESI-Werten in F13. Umbenennungen: ω_vir → ω_orb; ν als Peak-Höhe (Press-Schechter) → p_h. F3 und F8 rechnen mit a₀ = cH₀/2π = 1,08·10⁻¹⁰ (H₀ = 70), die korrigierten F8-Zahlen und die übrigen Abschnitte mit a₀ = 1,2·10⁻¹⁰. F1: r* des Universums 9,4 → 9,3 Gpc (Kern, a₀ = 1,2·10⁻¹⁰).

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Verwendet werden ausschließlich die Formeln des Modells: g = GM/r² = x/y²·g⊕, GM = 3c²Δs, dτ/dt = √(1−6Δs/r), r₀ = 6Δs, r* = √(3c²Δs/a₀), a₀ = c²/(2Λ) (heute ≈ cH₀/2π) [korr.]. Dazu kommen die Feldgleichung (Referenzlänge r*/3, im Fernfeld v⁴ = GM·a₀) und die Rückkopplungssumme a = a_N/(1−q), q = exp(−√(a_N/a₀)).
> Zahlenwert der Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen) in F1, F2 und F4–F6. F3 und F8 rechnen mit dem Modellwert c·H₀/2π = 1,08·10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc).
> Rechnung: `rechnungen/folgerungen.py`

---

## F1 — Jeder Körper hat drei Radien

| Radius | Formel | Bedeutung |
|---|---|---|
| r₀ | 6Δs | Takt steht still (dτ/dt = 0), tiefster Punkt der Falte |
| R | Oberfläche | Hülle, hier wirkt g(x,y) |
| r* | √(3c²Δs/a₀) | Schwelle: bis hier gilt der Einzeltrichter, jenseits das kollektive Feld |

Das ist die Triade des Modells als Längenmaß:
- Zentrum ↔ r₀
- Wirbel bzw. Hülle ↔ R
- Übergang zum Hintergrund bzw. Netzwerk ↔ r*

| Körper | r₀ = 6Δs | R | r* |
|---|---|---|---|
| Erde | 8,9 mm | 6 371 km | 12,2 AE |
| Sonne | 2,95 km | 696 000 km | 0,111 Lj |
| Neutronenstern | 4,1 km | 12 km | 0,13 Lj |
| Milchstraße | 1 180 AE | 15,3 kpc | 8,35 kpc |
| Galaxienhaufen | 53 Lj | 2 Mpc | 0,44 Mpc |
| Universum | 7,2 Gpc | 14,3 Gpc | 9,3 Gpc [korr.] |

## F2 — Einzeltrichter oder kollektiv: ein Kriterium aus g(x,y)

[korr.] Diese Einteilung setzt voraus, dass g(R) > a₀ einen Einzeltrichter kennzeichnet; das trägt nicht: Geglättet wird nur der Überschuss (ρ − ρ_J)₊ der Quelldichte über die Gezeitendichte ρ_J = 3T/(4πG); a₀ kommt darin nicht vor. Sterne sind Einzeltrichter, Sternsysteme kollektiv. Das Kriterium unten sagt nur, ob ein System als Ganzes unter der Schwelle liegt: Ein kompakter Kugelsternhaufen (g/a₀ ≈ 23) ist kollektiv, zeigt aber kaum einen Zusatz (ν = 1,008); sehr dünne Kugelsternhaufen sind kollektiv und zeigen ihn voll.

Ein Körper liegt als Ganzes **unter der Schwelle**, wenn seine eigene Hüll-Neigung unter a₀ liegt [korr.]:

```
g(R) < a₀   ⟺   x/y² · g⊕ < a₀   ⟺   y > √(x·g⊕/a₀)
```

Gleichwertig ist die Bedingung über die Flächendichte: g = πGΣ < a₀ ⟺ Σ < a₀/(πG) ≈ 274 M☉/pc².

| Körper | g(R)/a₀ | Regime |
|---|---|---|
| Erde, Jupiter, Sonne, Beteigeuze, Weißer Zwerg, Neutronenstern | 10⁷ bis 10²² | Einzeltrichter |
| **Kugelsternhaufen** (5·10⁵ M☉, 5 pc) | **23** | **kollektiv, ν = 1,008** [korr.] |
| **Zwerggalaxie** (10⁸ M☉, 1 kpc) | **0,12** | **kollektiv, ν = 3,4** [korr.] |
| Milchstraße | 0,30 | kollektiv |
| Galaxienhaufen | 0,05 | kollektiv |
| Universum | 0,43 | kollektiv |

**Folgerung:** Der sichtbare Zusatz setzt zwischen kompakten Kugelsternhaufen und Zwerggalaxien ein, weil ν von g/a₀ abhängt (1,008 gegen 3,4), nicht weil die einen Einzeltrichter wären [korr.]. So sieht es die Beobachtung [G]:
- Kugelsternhaufen zeigen keinen Dunkle-Materie-Effekt.
- Zwerggalaxien erscheinen am stärksten von „Dunkler Materie" dominiert.

Im Modell folgt die Stärke des Zusatzes aus g/a₀ über ν, die Einteilung in Einzeltrichter und kollektiv aus der Gezeitendichte ρ_J, ohne dass irgendwo Masse hinzugefügt wird. [korr.]

Die Grenzdichte a₀/(2πG) ≈ 137 M☉/pc² hat die Größenordnung der beobachteten Obergrenze für die zentrale Flächenhelligkeit von Scheibengalaxien (Freeman-Grenze, ~140 L☉/pc² [G]). Scheibengalaxien liegen also genau an der Schwelle.

## F3 — Die Schwelle verbindet den lokalen Horizont mit dem Ganzen

[korr.] Die Lesart unten setzt voraus, dass c²/(2a₀) die Länge des atmenden Ganzen (des Hubble-Radius) ist; das trägt nicht: Atmete sie mit dem Hubble-Radius, bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); mit den Schleifen atmend (Λ ∝ a) wirkte er wie eine Krümmung von 0,22 (gemessen Ω_K = 0,0007 ± 0,0037), so auch F10. Im Modell ist c²/(2a₀) = Λ die Schließlänge der Rückkehr, eine Länge der inneren Richtung und keine Raumlänge, als ganze innere Faltung der Zelle heute fest; Λ ≈ πR_H gilt nur heute.

Aus r*² = GM/a₀ und GM = c²r₀/2 folgt:

```
r*² = r₀ · c²/(2a₀)
```

Mit c²/(2a₀) = Λ ist r*² = r₀·Λ. Setzt man den heutigen Wert a₀ ≈ c·H₀/2π ein, dann ist Λ ≈ π·R_H mit dem Hubble-Radius R_H = c/H₀ [korr.]:

```
r* = √(r₀ · π·R_H) = √(6π · Δs · R_H)
```

**Der Schwellenradius ist das geometrische Mittel aus dem lokalen Horizont r₀ und der Schließlänge Λ der Rückkehr.** Nur heute ist Λ ≈ πR_H; in diesem Sinn trägt jeder Einzeltrichter die Größe des Hubble-Radius in seiner Reichweite. [korr.] Werte mit dem Modellwert a₀ = cH₀/2π = 1,08·10⁻¹⁰ m/s² (H₀ = 70): Sonne 0,117 Lj, Milchstraße 8,79 kpc. Das stimmt mit √(GM/a₀) überein, wenn a₀ = cH₀/2π gesetzt wird.

## F4 — Der Takt an der Schwelle ist die flache Geschwindigkeit

```
δ(r*) = 3Δs/r* = √(GM·a₀)/c² = (v_flach/c)²
```

Die Zeitabweichung an der Schwelle ist in erster Ordnung in Δs/r* exakt das Quadrat der flachen Rotationsgeschwindigkeit in Einheiten von c. Numerisch gilt die Gleichheit in dieser Ordnung für alle Körper auf alle Stellen:

| Körper | δ(r*) | v_flach |
|---|---|---|
| Zwerggalaxie | 1,40·10⁻⁸ | 35,5 km/s |
| Milchstraße | 3,44·10⁻⁷ | 175,8 km/s |
| Galaxienhaufen | 1,83·10⁻⁵ | 1 283 km/s |

## F5 — Faltungsmesser: die Faltung aus der Rotationskurve ablesen

Aus GM = 3c²Δs und v⁴ = GM·a₀ folgt:

```
Δs = v_flach⁴ / (3c² · a₀)
```

Aus der gemessenen flachen Geschwindigkeit einer Galaxie lässt sich ihre Faltungstiefe direkt bestimmen:

| v_flach [km/s] | 60 | 110 | 176 | 220 | 300 |
|---|---|---|---|---|---|
| Δs [AE] | 2,7 | 30 | 198 | 484 | 1 673 |

## F6 — Die Zusatzstrecken der Dunklen Materie

Die FTM-Aussage lautet: „Dunkle Materie ist der Summeneffekt der im Raum liegenden Zusatzstrecken Δs." Sie verbindet sich mit dem Operator der Zustandsraum-Theorie. Die wirksame Faltung jenseits von r* ist:

```
Δs_eff(r) = a(r)·r²/(3c²) = Δs · ν(r),     ν = 1/(1−q),  q = exp(−√(a_N/a₀))
```

| r [kpc] | 5 | 8,35 (r*) | 15,3 | 30 | 50 | 100 |
|---|---|---|---|---|---|---|
| ν | 1,23 | 1,58 | 2,38 | 4,12 | 6,50 | 12,5 |
| Zusatzstrecke Δs_eff − Δs [AE] | 46 | 115 | 273 | 615 | 1 087 | 2 268 |

Tief jenseits von r* wächst die Zusatzstrecke linear mit r: Δs_eff ≈ √(GM·a₀)·r/(3c²). **Die Zusatzstrecken wachsen, ohne dass Masse hinzukommt.** Das ist die Brücke zwischen dem FTM-Satz über die Summe der Zusatzstrecken und dem Operator 𝒫: Die Zusatzstrecke entsteht durch die Struktur, nicht durch eine Quelle.

## F7 — Ein kollektiver Trichter hat keinen flachen Außenraum

Jenseits von r* ist das Potential logarithmisch, Φ ≈ v_flach²·ln r. Das Faltungs-Differential zwischen r* und k·r* ist Δδ = (v_flach/c)²·ln k. Es wächst ohne Grenze:

| k | 2 | 10 | 100 | 1000 |
|---|---|---|---|---|
| Δδ (Milchstraße) | 2,4·10⁻⁷ | 7,9·10⁻⁷ | 1,6·10⁻⁶ | 2,4·10⁻⁶ |

**Folgerung:** Ein Einzeltrichter (Stern) klingt zum flachen Außenraum hin ab. Ein kollektiver Trichter (Galaxie) kann das nicht. Er **muss** in die Trichter seiner Nachbarn übergehen. Das Netzwerk aus Galaxien, Haufen und Filamenten ist im Modell damit eine logische Notwendigkeit und keine Zusatzannahme. Es entspricht der FTM-Aussage vom „vernetzten Gesamtgewebe".

## F8 — Das Ganze liegt auf seiner eigenen Schwelle

[korr.] Die Folgerung unten setzt voraus, dass die Neigung des Ganzen die Schwelle setzt (a₀ ∝ H); das trägt nicht (N_eff, Krümmung, so auch F10). Im Modell ist a₀ = c²/(2Λ) mit Λ, der ganzen inneren Faltung der Zelle, heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute.

- Die Hüll-Neigung des beobachtbaren Universums ist g_U = GM_U/R_U² = 5,11·10⁻¹¹ m/s² = 0,43·a₀ (M_U = 7,5·10²² M☉, R_U = 14,3 Gpc, a₀ = 1,2·10⁻¹⁰). [korr.]
- Die Hubble-Kugel bei kritischer Dichte hat g_H = c·H₀/2 = 3,4·10⁻¹⁰ m/s². Daraus folgt g_H/π = c·H₀/2π = 1,08·10⁻¹⁰ m/s², heute also a₀ [korr.].
- Der Schwellenradius des Universums beträgt r*_U = 0,65·R_U. [korr.]

**Folgerung:** Heute ist die Schwelle a₀ bis auf den Faktor π die Randneigung der Hubble-Kugel, und das Universum liegt auf dieser Schwelle. Diese Gleichheit gilt nur heute; die globale Neigung setzt die Regel nicht, sie trifft sie heute. [korr.]

## F9 — Die Regel wirkt auf die Bewegung, Licht sieht nur ihren wirbelfreien Teil

Drei Aussagen des Modells:

1. Eine statische Metrik gibt einem langsamen Körper die Beschleunigung −∇Φ. Sie ist immer wirbelfrei.
2. Die örtliche Bahnregel a = ν(|g_N|)·g_N ist in Scheiben nicht wirbelfrei, ∇×a = ∇ν × g_N ≠ 0. Die Rotationskurven stützen sie: Das senkrechte Scheibenfeld geht an SPARC mit f = +0,03 bis +0,13 ein (ohne IC 2574 0,00 ± 0,03), die örtliche Regel gibt in der Mittelebene f = 0, die Feldform f = 1 liegt außerhalb 95 %.
3. Licht folgt nur der Metrik.

Zerlegt man die örtliche Regel in einen wirbelfreien und einen quellfreien Teil, ist der wirbelfreie Teil genau die Feldform:

```
ν·g_N = −∇Φ_Feld + ∇×A,     ∇²Φ_Feld = ∇·[ν∇Φ_N]
```

**Folgerung:** Die Regel kann nicht allein über die Metrik wirken. Die Metrik trägt höchstens −∇Φ_Feld; den quellfreien Rest ∇×A erhalten nur massive Körper. Die Regel wirkt also auf die Bewegung, wie es der Satz des Autors sagt: Der Gedächtnisterm ändert die Bewegungsregel. Für Kugeln ist ∇×A = 0, dort sind Bahnen und Licht gleich.

| | Kugel, Punktmasse, weit außen | Scheibe innen (R ≲ 3,5 R_d) |
|---|---|---|
| Bahnen | ν(\|g_N\|)·g_N | in der Ebene ν(g_R)·g_R |
| Licht | dasselbe | Feldform, in der Ebene bei R_d/4 das 0,53–0,82-Fache der Bahnbeschleunigung, bei R_d das 0,79–0,92-Fache |
| M_dyn/M_lens | 1 | > 1 |

**Zweite Folge, Ruhe:** Die Gedächtnisform (Term ∝ ṙ) verschwindet für ruhende Körper; die stationäre örtliche Regel wirkt auch in Ruhe. Heißes Gas im hydrostatischen Gleichgewicht ruht im Mittel. In der Gedächtnisform sähe es nur Newton, in der stationären Lesart dieselbe Verstärkung wie die Sterne. Das trennt beide Lesarten in ein und derselben Galaxie. Gemessen: Die hydrostatischen Röntgenhalos isolierter elliptischer Galaxien (NGC 720, NGC 1521) zeigen die Verstärkung über 0,1–10 a₀. Die Regel wirkt stationär am Ort.

## F10 — Die Schließlänge der Rückkehr ist keine Raumlänge

Nach Kac ist a₀ = c²/(2Λ), Λ die Schließlänge der Halbgeraden der Rückkehr. Der Vakuumteil der Strukturwirkung hängt an derselben Länge:

```
Ω_s = (π²/45)·Θ²,   Θ = 2πa₀/(cH),   w_s = −1 − (2/3)·d ln a₀/d ln a
```

Drei Lesarten, jede mit Θ = 1 heute (Ω_m = 0,315, Ω_r = 9,1·10⁻⁵, H₀ = 67,4 km/s/Mpc):

| Λ atmet | a₀(t) | Ω_s bei der Nukleosynthese | Ω_s bei z = 1090 | w_s heute | a₀(1)/a₀(0) | Abstand zur letzten Streufläche gegen ΛCDM |
|---|---|---|---|---|---|---|
| mit dem Hubble-Radius, Λ = πc/H | ∝ H | 0,219 (wie ΔN_eff = 2,09) | 0,219 | −0,60 | 1,96 | 0,904 |
| mit den Schleifen des Torus, Λ ∝ a | ∝ 1 + z | 2·10⁻¹⁵ | 5·10⁻⁴ | −⅓ | 2,00 | 0,952 |
| nicht, Λ fest | konstant | 2·10⁻³³ | 4·10⁻¹⁰ | −1 | 1 | 1 |

Gemessen (Planck 2018 mit BAO): N_eff = 2,99 ± 0,17, Krümmung Ω_K = 0,0007 ± 0,0037.

**Folgerung:**
- Λ atmet weder mit dem Horizont noch mit den Schleifen des Torus. Mit dem Horizont bliebe der Vakuumteil früh bei 22 %; mit den Schleifen wirkte er im Ausdehnungsgesetz wie eine Krümmung von 0,22. Beide verdoppelten a₀ bis z = 1.
- Λ ist also keine Raumlänge, sondern eine Länge der inneren Richtung der Rückkehr. Die räumliche Zelle des Kaleidoskops atmet mit, L(t) = L̄·a(t); über ihre ganze innere Faltung setzt sie Λ (F13).
- Früh muss a₀/H gegen null gehen. Dann fehlt der Vakuumteil bei Nukleosynthese und Rekombination, wie es die Messung verlangt.
- Ein fester Vakuumteil hielte die Atemrate dauerhaft über einem Mindestwert; der Kreislauf verlangt eine fallende Atemrate. Ob dafür Λ wächst oder die Rückstellung, ist offen. Heute ist Λ ≈ πc/H₀ ≈ 14 Gpc.
- Was Λ = πc/H als Raumschleife voraussetzte, entfällt: die Grenzen für Ω_m und H₀ und die Lesart Q = 2π − 1 (bei festem Λ wäre ν = 2ΛH/c zeitabhängig). [korr.] Q = 2π − 1 (ν = 2π) ist Zahlennähe ohne Herleitung und trägt nicht; Q = 5,364 ist aus der Messung eingesetzt, im homogenen Hintergrund gibt die Rückkehr kein q (Γ = 0). Über die innere Faltung bindet die Zelle a₀ neu (F13).

## F11 — Das Gedächtnis ist eins

Im Modell stehen vier Stücke: die bleibende Knickspur, die Mode ξ mit dem Gewicht w (tiefster Kollaps, Preisach), ihr Übergang in die Operatorform mit Γ = (1 − w)·ω_orb/N* und der Nachlauf λ_m = ω_orb/N* [korr.], mit demselben N* ≈ 0,3–2,4.

- ξ startet im Zustand der Baryonen, fällt druckfrei und trägt ihre Bahngeschichte weiter. Wird Gas abgebremst, läuft ξ auf der alten Bahn weiter. Die stationäre Regel sitzt um die heutigen Baryonen. Geht ξ mit der Rate Γ in sie über, wandert die Struktur mit derselben Rate zum Gas. Das ist der Nachlauf: λ_m = Γ(w → 0) = ω_orb/N*.
- Bleibend ist w: Ein einmal gebundener Bereich bleibt gebunden. Das ist die Knickspur als Zustand der Regel. Die Rillen sind die Mulden von ξ, in die nach der Rekombination die Baryonen fallen.

| System | 1/ω_orb | 1/λ_m = N*/ω_orb |
|---|---|---|
| Galaxie, Kollaps bei z = 2 (200-fache kritische Dichte damals) | 0,48 Mrd. Jahre | 0,14–1,15 Mrd. Jahre |
| Haufen, 500-fache kritische Dichte | 0,92 Mrd. Jahre | 0,28–2,2 Mrd. Jahre |
| Haufen, 200-fache kritische Dichte | 1,45 Mrd. Jahre | 0,44–3,5 Mrd. Jahre |

Der Bullet Cluster verlangt 1/λ_m ≥ 0,17–0,26 Mrd. Jahre, die Haufenwerte erfüllen das.

**Folgerung:**
- In alten Galaxien ist die Übergabe abgeschlossen (0,1–1,2 Mrd. Jahre nach dem Kollaps). Die Regel ist dort stationär, ohne eigene Streuung. Gezeitenzwerge und Galaxien nach Ausbrüchen folgen ihr sofort.
- In Haufen bleibt die Übergabe unvollständig. Der stoßfreie Rest von ξ folgt den Galaxien (Bullet Cluster 79 %), der Operatorteil läuft mit 1/λ_m zum Gas.
- Als Geschwindigkeitskern auf jeder Bahn wäre λ_m am Sonnenkreis 0,01–0,13 Ω. Exzentrische Bahnen würden dann in wenigen Umläufen rund; alte Scheiben verlangten λ_m ≪ 0,01 Ω. Das Gedächtnis wirkt also auf den Zustand der Regel, nicht auf die einzelne Bahn.

## F12 — Die örtliche Regel hat keine Wirkung

Die örtliche Regel ist F = ν(|g_N|)·g_N = −ν∇Φ_N.

- F steht überall senkrecht auf den Newton-Äquipotentialflächen, F·(∇×F) = 0. 1/ν ist der integrierende Faktor, die Arbeit längs eines Weges ist −∫ν dΦ_N.
- Eine Wirkung gibt ruhender Materie immer einen Gradienten (Gleichgewicht aus einem Energieminimum). F ist genau dann ein Gradient, wenn ν auf den Newton-Äquipotentialflächen konstant ist: bei Kugeln und über der Kuzmin-Scheibe.
- Geschwindigkeitsabhängige Glieder einer Wirkung helfen nicht. Ungerade in v wären sie vom Drehsinn abhängig; die Regel ist aber gerade in ω. Eine gerade Trägheit ändert ruhende Materie nicht. Die hydrostatischen Röntgenhalos verlangen aber die Verstärkung in Ruhe, also in der Kraft.

**Folgerung:** Außerhalb der Kugel hat die örtliche Regel keine Wirkung. Sie setzt am Erzeuger der Bewegung an, L_eff = 𝒫L, nicht an der Energie. Auf geschlossenen Wegen tauscht sie Arbeit mit dem Medium. Die einfachste Wirkung, die Feldform, liegt an SPARC außerhalb 95 %.

Wie stark die Regel nicht wirbelfrei ist, gemessen als Arbeit eines geschlossenen Umlaufs gegen die Arbeit auf dem Hinweg (G = M = a₀ = 1):

| System | innen | bei r* bzw. R_d–3R_d | außen |
|---|---|---|---|
| abgeplattete Ellipse, Miyamoto-Nagai a = 0,25, b = 0,05–0,25 | 4–5 % (R = 0,3) | 1–2 % (R = 1) | 0,1–0,3 % (R = 3) |
| Exponentialscheibe, Umlauf durch die Schicht, R_d = 0,3 / 1 / 3 r* | 11 / 25 / 30 % (R_d/4 bis R_d) | 5 / 9 / 11 % | 1 % (3–6 R_d) |

In der Scheibe sitzt der Wirbel in der Schicht selbst: In der Ebene gilt ν(g_R), knapp darüber ν(√(g_R² + (2πGΣ)²)).

Ruhendes Gas bleibt im Gleichgewicht. Aus ∇P = ρF folgt P = P(Φ_N) und ρν = −dP/dΦ_N. Der Druck ist also auf den Newton-Äquipotentialflächen konstant, und die Temperatur wächst längs einer Isobare mit ν. In den Ellipsen oben sind das innen höchstens 6 %, bei r* 1–2 %.

Auf Sternbahnen mittelt sich der Austausch weg, weil sich die Bahnen in (R, z) nicht schließen. In einer Scheibe wie der Milchstraße (Scheibe 5·10¹⁰ M☉, R_d = 2,6 kpc, Dicke h = 0,15–0,3 kpc; Bulge 10¹⁰ M☉), gerechnet über 10 Mrd. Jahre gegen eine wirbelfreie Kraft mit derselben Rotationskurve (224 km/s bei 8 kpc):

| Scheibe | ΔE Mittel / Streuung | Anteil an σ_z² | Dicke: Regel / wirbelfrei |
|---|---|---|---|
| h = 0,3 kpc, σ_z = 18 km/s | +1,8 / 10 (km/s)² (3D: +0,8 / 7) | 0,3–0,6 % | +0,5 / +0,5 % (3D: +0,4 / +0,4 %) |
| h = 0,3 kpc, σ_z = 30 km/s | +7,2 / 22 (km/s)² | 0,8 % | +0,7 / +0,3 % |
| h = 0,15 kpc, σ_z = 12 km/s | +2,2 / 11 (km/s)² | 1,5 % | +0,6 / +0,5 % |

Die Dicke alter Scheiben bleibt vom Arbeitsaustausch unberührt.

## F13 — Die Dunkle Energie im Ganzen

**Rückstellung und Schwelle sind eine Länge.** Die Rückkehr läuft durch die innere Richtung der Faltung und schließt bei Λ. Deren ganze Länge ist die ganze innere Faltung der geschlossenen Zelle; nach der Bilanz des Autors (σ = ΣΔs_innen − ΔL_außen = 0) ist sie die äußere Dehnung. Der Vakuumteil der Strukturwirkung hängt nur an Λ, Λ_s = π⁴/(15Λ²): Er ist die Energie der Rückstellung.

```
Λ = ΣΔs_innen(Zelle) = ΔL_außen = G·ρ_m·L³/(3c²),     a₀ = c²/(2Λ) = 3c⁴/(2G·ρ_m·L³)
```

| Größe | Wert (Ω_m = 0,313, H₀ = 67,4 km/s/Mpc) |
|---|---|
| a₀ bei L = 27 Gpc | 1,174·10⁻¹⁰ m/s² (von Hand nachgerechnet) |
| L für a₀ = 1,04 / 1,10 / 1,19 / 1,2·10⁻¹⁰ | 28,11 / 27,59 / 26,88 / 26,80 Gpc |
| Durchmesser der letzten Streufläche d_LSS | 27,79 Gpc; bei L = d_LSS ist a₀ = 1,077·10⁻¹⁰ |
| ohne Kreispaare, kürzeste Schleife durch unseren Ort L ≥ 0,985·d_LSS (COMPACT 2023; Planck 2015, kubischer Torus: 0,97) | a₀ ≤ 1,13·10⁻¹⁰ |
| a₀ früh (Strahlung erhöht die Faltung) | z = 10: 0,997; z = 1090: 0,76; z = 3400: 0,50; Nukleosynthese: 0 |
| Vakuumteil | heute 0,234 (Θ = 1,03), z = 1: 0,073, z = 1090: 2,5·10⁻¹⁰ |

- Die innere Faltung je mitbewegter Zelle ist heute fest ((1 + Q)·ΣΔs_b); Λ ist seit dem Materiezeitalter praktisch fest (Kreislauf unten), der Vakuumteil eine Konstante, früh null. Das erfüllt N_eff und die Hintergrundstrahlung von selbst.
- Die Schwelle wird an die Zelle gebunden: Für den gemessenen Bereich von a₀ liegt die Zelle am Durchmesser der letzten Streufläche, prüfbar an Quadrupol, Winkelkorrelationen und Kreispaaren.

**Die übrigen zwei Drittel.** Mit fester Faltung gibt die Bilanz keine Zeitform. Der Autor sagt aber auch: Die Atemrate ist die zeitliche Änderung der Rückstellung. Als Rate gelesen folgt ein Teil der Rückstellung der Bindung; gerechnet mit ρ_R ∝ dF/d ln a, Vakuumteil fest, selbstkonsistent (F nach Press-Schechter mit δ_w und nach Sheth-Tormen):

| Verfahren, kleinste gebundene Masse | Maximum von ρ_DE bei z | w_DE(0) | CPL w₀ / w_a |
|---|---|---|---|
| PS-w 10¹⁰ / 10¹² / 10¹³ M☉ | 3,9 / 2,0 / 1,2 | −0,70 / −0,72 / −0,74 | −0,68 / −0,08 · −0,66 / −0,32 · −0,64 / −0,65 |
| ST 10¹⁰ / 10¹² / 10¹³ M☉ | 2,2 / 1,1 / 0,6 | −0,74 / −0,77 / −0,80 | −0,71 / −0,26 · −0,67 / −0,66 · −0,64 / −1,17 |

Das ist der gemessene Quadrant (Pantheon+ mit DESI DR2 und Hintergrundstrahlung: −0,84…−0,86 / −0,58…−0,62; Maximum bei z ≈ 0,35–0,44 [G]). Die Lage des Maximums hängt an der kleinsten gebundenen Masse; für z ≈ 0,4 läge sie über 10¹³ M☉ (mit der Kollapsschwelle 1,69 statt δ_w bei 10¹³ M☉), und w₀ bleibt dort 0,05–0,2 zu hoch. Das Modell legt diese Masse nicht fest.

- **M_Q legt sie nicht fest.** Bei festem a₀ ist M_Q ∝ H⁻⁴. Mit der Peak-Höhe p_h = δ/(σ(M_Q)D) und σ ∝ M^{−α} (α ≈ 0,27) gilt d ln p_h/d ln a = 4αε − f [korr.]: im Materiezeitalter +0,6, heute ≈ 0. Der Anteil oberhalb M_Q fällt also früh und steht heute; seine Rate wäre negativ.
- **Die Straffung braucht keine Mindestmasse**, trifft die Form aber nicht. Rückstellung ∝ Σ(Δs_i/R_i) der gebundenen Trichter (Gewicht M^{2/3}), R aus 200-facher Dichte:

| R aus | gelesen als | Maximum bei z | CPL w₀ / w_a (PS-w · ST) |
|---|---|---|---|
| mittlerer Materiedichte | Rate je ln a | 0,90 · 0,66 | −0,37 / −1,40 · −0,47 / −1,46 |
| mittlerer Materiedichte | Rate je Zeit | 1,97 · 1,36 | −0,31 / −0,76 · −0,34 / −1,07 |
| mittlerer Materiedichte | Bestand | wächst | −1,13 / −0,34 · −1,19 / −0,26 |
| kritischer Dichte | Rate | Rate wird negativ | – |
| kritischer Dichte | Bestand | wächst | −1,34 / +0,13 · −1,41 / +0,25 |

**Kreislauf.** Nach dem Autor verdichtet sich das Medium an der Engstelle wieder zu Masse; es entsteht neue innere Faltung, Λ wächst, der Vakuumteil fällt. Trägt der Vakuumteil die Atemrate, H² = c²Λ_s/3 = π⁴c²/(45Λ²), ist die Horizontentropie die der ganzen inneren Faltung:

```
S_max = πkc²/(ℓ_P²H²) = (45/π³)·k·(Λ/ℓ_P)²      (lokale Engstelle: S = 36πk(Δs/ℓ_P)²)
dS_max/dt ≥ Ṡ_i   ⇔   Λ̇/Λ ≥ Ṡ_i/(2S_max) = ℓ_P²H²Ṡ_i/(2πkc²)
```

| Größe | Wert (Λ = 13,5 Gpc, Ṡ_i = 7,1·10⁸⁶ k/s je Horizont) |
|---|---|
| S_max im Endzustand | 9,6·10¹²² k |
| H im Endzustand | 0,49·H₀ |
| Mindestrate Λ̇/Λ | 3,7·10⁻³⁷ s⁻¹ = 1,2·10⁻²⁹ je Jahr |
| neue Faltung in der ganzen Zelle (27,8 Gpc) | 6·10¹⁷ kg/s ≈ 10⁻⁵ M☉ je Jahr |
| Wirkung auf den Vakuumteil | w_s = −1 + 2·10⁻¹⁹ im Endzustand, heute −1 + 1·10⁻¹⁹ |
| Verlust durch Rotverschiebung der Strahlung heute | d ln Λ/d ln a = −Ω_r/Ω_m = −2,9·10⁻⁴ (w_s = −1,0002), verschwindet mit 1/a |

Das Axiom „kein Endzustand“ verlangt also eine neue innere Faltung, deren Mindestrate die Entropieproduktion festlegt; auf die Dunkle Energie heute wirkt sie nicht. Wo und wie das Medium sich an der Engstelle verdichtet, legt das Modell noch nicht fest.

**Offen:** eine Größe, die die kleinste gebundene Masse der Rückstellung festlegt; Ort und Weg der neuen inneren Faltung an der Engstelle.

---

## Neue prüfbare Aussagen aus F1–F13

1. Die Grenze Einzeltrichter/kollektiv setzt die Gezeitendichte ρ_J, nicht a₀: Sterne sind Einzeltrichter, Sternsysteme kollektiv. Die Stärke des Zusatzes folgt aus g/a₀: Kompakte Sternsysteme (g/a₀ ≈ 23) zeigen kaum einen (ν = 1,008), diffuse mit g(R) < a₀ den vollen; auch sehr dünne Kugelsternhaufen sind kollektiv. Prüfbar an ausgedehnten, lichtschwachen Kugelsternhaufen und ultrakompakten Zwergen. [korr.]
2. Der Faltungsmesser Δs = v_flach⁴/(3c²a₀) liefert für jede Galaxie die baryonische Masse aus der Rotationskurve.
3. Die Schwelle liegt heute bei cH₀/2π; in der Zeit folgt sie H nicht.
4. Scheibengalaxien: innen liegt die Linsenmasse unter der Bahnmasse, außerhalb von etwa 3,5 R_d sind beide gleich.
5. Heißes Gas im hydrostatischen Gleichgewicht zeigt dieselbe Zusatzbeschleunigung wie die Sterne derselben Galaxie; an isolierten elliptischen Galaxien gemessen.
6. Die Beschleunigungsrelation früher Scheiben (z = 1–2) zeigt keine Verdopplung der Schwelle.
7. Junge Gezeitenzwerggalaxien liegen schon auf der Beschleunigungsrelation; einen Nachlauf gibt es nur, wo ξ nicht übergegangen ist (Haufen, frisch kollabierte Bereiche).
8. Heißes Gas in abgeplatteten Systemen: Isobaren folgen den Newton-Äquipotentialflächen der Baryonen, die Temperatur längs einer Isobare wächst mit ν (Ellipsen innen bis 6 %).
9. Die Beschleunigungsrelation bei z = 0,5–1 legt die Zustandsgleichung des Vakuumteils fest: w_s = −1 − ⅔·d ln a₀/d ln a; mit Λ als ganzer innerer Faltung wandert sie nicht.
10. Die Zelle liegt am Durchmesser der letzten Streufläche (26,9–28,1 Gpc); ohne Kreispaare (L ≥ 0,985·d_LSS) ist a₀ ≤ 1,13·10⁻¹⁰.
