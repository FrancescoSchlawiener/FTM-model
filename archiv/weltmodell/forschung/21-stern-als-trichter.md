# Der Stern als Faltungstrichter: eine Temperatur verbindet Kosmos und Sterne

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „kT_× ist die Temperatur der Proton-Kreuzung bei der Rekombination, dort wird das Proton kosmisch zum Einzeltrichter“ → setzt a₀ ∝ H voraus, trägt nicht; a₀ seit z ≈ 10 fest, r*_p/λ̄_p bleibt über 145; Merksatz und Zusammenfassung „dieselbe Größe auf drei Skalen“ als Zahlennähe markiert; Ausblick „Selbstkonsistenz mit globalem Atmen (a₀ ∝ H)“ → a₀ ∝ H trägt nicht. Abschnitte 1 und 3 (Δs_* = 912 m, Carter-Exponenten) bleiben. [G] ergänzt: Sonnenwerte (Zentraltemperatur, Leuchtkraft, Oberflächentemperatur).

> Rechnung: `rechnungen/stern_trichter.py`

## 1. Warum Sterne etwa eine Sonnenmasse haben — in Faltungsgrößen

Die Massenskala der Sterne ist M_* = m_Pl³/m_p² (Landau/Chandrasekhar) = 1,85 M☉. Ihre Faltungstiefe:

```
Δs_*² · Δs_p = (λ̄_p / 3)³            ⟺            3Δs_*/λ̄_p = α_G^{−1/2}
```

| Trichter | Δs | 3Δs/λ̄_p |
|---|---|---|
| Proton | 4,1·10⁻⁵⁵ m | α_G = 5,9·10⁻³⁹ |
| Stern-Skala M_* | **912 m** | α_G^{−1/2} = 1,3·10¹⁹ |
| Chandrasekhar-Grenze (1,44 M☉) | 709 m | |
| Sonne | 492 m | |

**Lesart:** Auf der logarithmischen Leiter der Faltung liegt der Stern so weit über der Compton-Länge (Faktor α_G^{−1/2}), wie die Faltung des Protons unter ihr liegt (Faktor α_G, halbe Schrittweite). Die typische Sternfaltung von einigen hundert Metern folgt allein aus Δs_p und λ̄_p.

## 2. Die Photosphäre in Faltungsgrößen

Ein Stern ist ein Trichter, in dem sich zwei Rollen des Zustandsraums die Waage halten:
- **parametrisch (Faltung):** Hydrostatik, kT_c ~ G·M·m_p/R = 3c²Δs·m_p/R
- **tangential/transversal (Elektromagnetik):** Strahlungstransport mit dem Thomson-Querschnitt σ_T = (8π/3)·r_e². Dabei ist r_e = α·λ̄_e der elektromagnetische Radius des Elektrons (α = r_e/λ̄_e).

Mit M ~ M_* und T_c ~ α²m_p c² folgt für die Oberflächentemperatur:

```
kT_s  ~  √( α·m_e c² · kT_× ),        kT_× = m_p c²·α_G^{1/4} = 0,2601 eV
```

**kT_× läge nur unter a₀ ∝ H auf der Temperatur einer Proton-Kreuzung.** [korr.] Dort fiele die Horizonttemperatur des atmenden Kosmos (ħH/2π) auf die Gravitations-Selbstenergie des Protons (α_G·m_p c²), und das geschähe bei der Rekombination. [korr.] Diese Gleichheit ist r*_p = λ̄_p in gemischter Zählung (H/2π pro Umlauf gegen α_G·ω_p pro Radiant); die Nähe zur Rekombination hängt an dieser Zählung, Carters Beziehung nicht. [korr.] Als Kreuzung setzt das a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute; r*_p/λ̄_p bleibt über 145 (bei der Rekombination 166), eine kosmische Kreuzung gibt es nicht. kT_× bleibt als Energieskala m_p c²·α_G^{1/4}.

> **Ein und dieselbe Temperatur kT_× legt fest, wann das Proton kosmisch zum Einzeltrichter wird und auf welcher Temperaturskala die Oberflächen der Sterne liegen.** [korr.] Der erste Teil gilt nur unter a₀ ∝ H und trägt nicht; es bleibt die Temperaturskala der Sternoberflächen.

## 3. Carter folgt als Balancebedingung

Die Photosphäre liegt an der Ionisationsgrenze von Wasserstoff, also am Übergang von Strahlungs- zu Konvektionshülle, wenn gilt:

```
kT_s ~ α² m_e c²     ⟺     kT_× ~ α³ m_e c²     ⟺     α_G ~ α¹²·(m_e/m_p)⁴      (Carter)
```

Im Trichter-Bild ist Carters Koinzidenz damit die Bedingung: **Die parametrische Faltung (Hydrostatik) und der tangentiale Strahlungstransport (σ_T ∝ r_e²) treffen sich genau an der Bindungsgrenze des Wasserstoffs.** Die Exponenten 12 und 4 entstehen aus den Potenzen in Hydrostatik, Opazität und Stefan-Boltzmann:
- T_s⁴ ~ α²·m_p³·m_e²/m_Pl
- gleichgesetzt mit (α²m_e)⁴
- ergibt α⁶m_e² = m_p³/m_Pl, quadriert α¹²m_e⁴ = m_p⁶/m_Pl², also α¹²(m_e/m_p)⁴ = α_G.

## 4. Vorfaktor-Prüfung

| Größe | Grobrechnung (Thomson, Virial) | Sonne gemessen |
|---|---|---|
| T_c | 1,16·10⁷ K | 1,57·10⁷ K (Modell) [G] |
| L | 2,8·10²⁸ W | 3,83·10²⁶ W [G] |
| T_s | ~17 000 K | 5 772 K [G] |

- Die reine Skala √(α m_e c²·kT_×) = 31 eV liegt ohne Vorfaktoren weit über den 0,5 eV der Sonne.
- Die Grobrechnung mit Virial und Thomson-Opazität trifft T_c auf Faktor 1,4, L auf Faktor 70 und T_s auf Faktor 3.
- Die Skalen und die **Exponenten** stimmen. Die Vorfaktoren brauchen die reale Opazität (Kramers, H⁻), die Fusion und die Temperaturgradienten.

## 5. Einordnung

| Aussage | Stand |
|---|---|
| Δs_*² · Δs_p = (λ̄_p/3)³ (Sternfaltung aus Protonfaltung) | exakte Umformung der Landau-Chandrasekhar-Skala, neue FTM-Lesart |
| kT_s ~ √(α m_e c² · kT_×): Sternoberflächen an der Proton-Kreuzungstemperatur | Größenordnungsrechnung, Exponenten exakt |
| Carter = Balance Faltung/Strahlungstransport an der Wasserstoff-Bindungsgrenze | Standard-Herleitung (Carr & Rees-Typ, aus dem Gedächtnis), im FTM neu gelesen |
| **Warum** unser Universum auf dieser Balance sitzt | offen; bisher anthropisch begründet |

**Zusammengefasst:** Die Spur „Compton-Länge des Protons zur Rekombinationszeit" führt zu einem Knoten, an dem drei Dinge zusammenhängen:
1. das kosmische Atmen,
2. die Faltung des einzelnen Protons,
3. die Oberflächentemperatur der Sterne.

Im FTM ist das dieselbe Größe, kT_× = m_p c²·α_G^{1/4}, gelesen auf drei Skalen. Das ist die fraktale Selbstähnlichkeit des Modells an einer konkreten Zahl. [korr.] Zahlennähe ohne Herleitung, trägt nicht; ohne die Kreuzung (a₀ ist seit z ≈ 10 fest) fehlt der Bezug zum kosmischen Atmen, es bleiben Proton und Sternoberflächen.

**Offen und nächster Schritt:** Die Frage „warum sitzt die Welt auf dieser Balance?" lässt sich im FTM so stellen: Stellt sich der Kopplungswert α_G durch die Geschlossenheit des Torus selbst ein? Eine Selbstkonsistenz über ein mit H atmendes a₀ (a₀ ∝ H) trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 % (ΔN_eff = 2,09 gegen gemessen N_eff = 2,99 ± 0,17); im Modell ist Λ die ganze innere Faltung der Zelle und heute fest. [korr.] Das ist rein spekulativ, bis eine Gleichung dafür vorliegt.
