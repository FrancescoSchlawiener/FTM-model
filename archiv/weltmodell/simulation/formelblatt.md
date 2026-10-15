# Formelblatt des FTM-Modells

Die wenigen geometrischen Formeln, die das Ganze tragen, geordnet nach den drei Rollen im Zustandsraum. Je Formel ein einfaches Beispiel. Die Zahlen rechnet `python formelblatt_beispiele.py` aus dem Kern; die Prüfung aller übrigen Zahlen läuft im Hintergrund (`lauf.py`).

**Ein Medium, geschlossen, ohne Außen.** Jede Erscheinung ist die Rückkehr eines geschlossenen Weges, gewichtet mit einem Nome q = e^{2πiτ}:
- **Betrag |q|:** Verlust und Rückkehr (transversal).
- **Phase:** Windung (tangential).
- **Rahmen:** Er setzt, wie weit zurückgekehrt wird (parametrisch).

## 1 Rahmen (parametrisch): Faltung und Schließung

| Formel | Bedeutung | Beispiel |
|---|---|---|
| Δs = GM/(3c²) | Masse ist Faltung des Mediums (die 3 ist die Raumdimension) | Erde 1,478 mm · Sonne 492,2 m · Milchstraße 197,4 AE |
| a₀ = c²/(2Λ) | Schwelle aus der Schließlänge Λ der inneren Richtung (keine Raumlänge) | Λ = 13,2 Gpc bei a₀ = 1,1·10⁻¹⁰ m/s² |
| τ_D = 2Λ/c | Kac-Zeit: mittlere Rückkehr auf der bei Λ geschlossenen Richtung, für jede Masse gleich | 86,4 Mrd. Jahre |

## 2 Rückkehr (transversal): Verlust, Verstärkung, Regel

| Formel | Bedeutung | Beispiel (Milchstraße, 6·10¹⁰ M☉) |
|---|---|---|
| s = √(a_N/a₀), q = e^{−s} | Anteil, der je Umlauf zurückkehrt | 2 kpc: q = 0,013 · 8,7 kpc: 0,367 · 30 kpc: 0,748 |
| ν = 1/(1 − q), a = ν·a_N | Summe der Rückläufe = Verstärkung (Schleife mit Verstärkung q) | ν = 1,01 · 1,58 · 3,97 |
| r* = √(GM/a₀) | Radius, an dem die Rückkehr übernimmt (s = 1) | Sonne 7342 AE · Milchstraße 8,72 kpc |
| v⁴ = GM·a₀ | tiefer Bereich: flache Kurve (Tully-Fisher) | 172 km/s |
| û = e^{−√((p + Γ + θ)·τ_D)} | die eine Rückkehrgleichung, am Ort gelesen: p = θ = 0 gibt q, Γ = 0 gibt das Gewicht w, p ≠ 0 den Gedächtniskern | – |
| Labor: RC-Kette mit Rückführung | dieselbe Form, nachgebaut | ν = 1,1567 (20 Glieder, RG = 0,01) |

**Die Regel als Winkel der Windung.** Eine Bahn ist ein Torus: Umlauf Ω auf dem Hauptkreis, radiale Schwingung κ auf dem Querring. Die Windung κ/Ω kippt mit der Rückkehr:

tan²α = κ²/Ω² = 1 + χ(s),   cos²α = 1/(2 + χ)   (Pythagoras der Anteile parametrisch : transversal = 1 : 1 + χ)

| Stelle | α | parametrisch : transversal |
|---|---|---|
| Kepler und Coulomb (Sonnensystem, Atom) | 45,0° | 0,50 : 0,50 |
| Milchstraße bei r* (8,7 kpc) | 51,5° | 0,39 : 0,61 |
| Milchstraße bei R_Q (51 kpc) | 54,2° | 0,34 : 0,66 |
| tiefer Bereich, flache Kurve | 54,74° | 1/3 : 2/3 |

## 3 Windung (tangential): Phase ohne Verlust

| Formel | Bedeutung |
|---|---|
| arg q = 2π·Re τ | Windung: Ladung, Spin, Umlauf, Verschränkung |
| ΣQ_el = 0 | geschlossenes Budget der Windung |
| Windung ändert den Takt nicht | Das Argument der Regel ist das Taktgefälle g_N, nie qE/m; darum hat der Elektromagnetismus keine Schwelle. Joule-Wärme J·E ist der Kanal von der Windung in die transversale Rolle (dS/dt ≥ 0), keine Phase |

## 4 Leiter: dieselbe Form auf jeder Ebene

| Formel | Beispiel |
|---|---|
| r*² = 3Δs·ℓ₀ (ℓ₀ = 2Λ, groß) | Sonne: r*²/(3Δs·ℓ₀) = 1,000000 |
| ℓ_P² = 3Δs·λ̄ (λ̄ = ħ/Mc, klein) | Proton: 1,000000 |

Die Faltung Δs steht zwischen einer Gegenlänge am Ort (λ̄) und einer der inneren Richtung (ℓ₀).

## 5 Das Ganze

| Formel | Bedeutung | Zahl |
|---|---|---|
| Q = q/(1 − q) = Ω_c/Ω_b | mittlere Zahl der Rückläufe vor der Flucht; zeitfest, eine Form wie π | 5,364 gemessen; Zählung je Radiant 2π − 1 = 5,283 |
| Σ1/⟨T_i⟩ = 1 | Kac-Budget: Die Schließung verteilt ein festes Budget an Rückkehr | – |
| Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH) | Vakuumteil: Rückstellung des Ganzen, als Energie gelesen | Θ = 1,055, Ω_s = 0,244 |
| Bindung als Operator | Die Bindung macht Fluchtmoden inaktiv; sie gehen in die innere Richtung je Raumpunkt über und wirken wie der Vakuumteil (w = −1, w_a < 0). Dunkle Energie ist die Kehrseite der Bindung, keine Quelle | offen: wie viele Zyklen je Raumpunkt; zusammen verlangt θ* etwa das Dreifache des Vakuumteils |

## 6 Topologie des Trichters (symbolisch bewiesen)

| Aussage | Formel | Beispiel |
|---|---|---|
| Takt je Querring gleich ⇔ keine Quelle am Ort | n·T = 1 ⇔ g_tt·g_rr = −1 ⇔ ρc² + p_r = 0; im Vakuum-Trichter exakt: Ringdichte = 1/Takt | gemessen als γ = 1 |
| Masse = Überschuss an Querringen | Δs = ∫(√g_rr − 1)dr = GM/(3c²) (homogene Kugel; die 3 aus ∫r²/R³) | Sonne 492,2 m |
| Der Trichter schließt sich | Gauß-Bonnet: ∫K dA = 2π·(1 − T). Bei 6Δs ist T = 0 und die Krümmung voll 2π: Der Horizont ist ein Wechsel der Topologie | – |
| Stufen | 9Δs: Ringdichte √3, T/r extremal (Photonensphäre) · 18Δs: κ = 0, die Windung schwingt nicht mehr (innerste stabile Bahn) | – |
| Die Windung spürt Takt und Ringe gleich | Ablenkung 6Δs/b (Takt) + 6Δs/b (Ringe); Perihel (12π + 12π − 6π)·Δs/p (Takt + Ringe + Nichtlinearität) | Merkur 28,65″ + 28,65″ − 14,33″ = 42,98″/Jh |
| Nur die Bahn sieht die Windung der Regel | νg_N = −∇Φ + ∇×A; Metrik und Linse tragen den wirbelfreien Teil, die Bahn auch den Wirbel: η = M_dyn/M_lens = 1 + Wirbel/wirbelfrei | Kugel η ≤ 1,05; Scheibe innen bis 1,2–1,9, Vorzeichenwechsel bei 3,4R_d |

## Die drei Rollen

- **tangential:** Zyklen, Phase, Windung. Beispiele: Umlauf, Schwingung, Ladung.
- **transversal:** Verlust, Dämpfung, Rückkehr. Beispiele: q, ν, Wärme, Gedächtnis.
- **parametrisch:** langsame Verformung des Rahmens. Beispiele: Faltung Δs, Schließlänge Λ, Atmen, Gravitation.

Dunkle Materie und Dunkle Energie sind beide keine Quelle, sondern Änderungen der zugänglichen Bewegung:
- Die Regel q, ν ist die Rückkehr je Faltung.
- Die Dunkle Energie ist die inaktive Flucht je Raumpunkt.
