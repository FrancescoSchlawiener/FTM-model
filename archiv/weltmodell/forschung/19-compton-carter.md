# Was hinter der Proton-Kreuzung steckt: Carters Koinzidenz

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Kreuzung r*_p = λ̄_p bei z ≈ 1100 (setzt a₀(z) = cH(z)/2π voraus) → trägt nicht, r*_p/λ̄_p = 145 (a₀ = 1,2·10⁻¹⁰), 156 (cH₀/2π), 166 bei z = 1090 (0,76·a₀), nie 1; r* = λ̄ nur bei M = 6,1·10⁻²⁹ kg; Ergebnisse, Status und die Lesart „Balance tangential/parametrisch“ nur noch unter dieser Annahme bzw. entfallen; „Carter nicht hergeleitet, Herleitung über gleiche Kopplungsraten gesucht“ → Exponenten 12 und 4 im Modell aus der Sternkette hergeleitet, ein Einrasten einzelner Zyklen erzeugt Carter nicht, offen nur 𝒞 = 2,94; Zahlennähe markiert (Saha-Logarithmus nahe 1/(2α), „zwei Abweichungen heben sich auf“). Umbenennung (ohne Marke): L (Saha-Logarithmus) → l_S. [G] ergänzt: Carter 1974. Carters Beziehung selbst (𝒞 = 2,94) bleibt.

> Rechnung: `rechnungen/compton_carter.py`

## Der Befund

[korr.] Die Kreuzung „Neigung des Protons an seiner Compton-Grenze = Schwelle a₀(z)“ setzt a₀(z) = cH(z)/2π voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute. Damit ist r*_p/λ̄_p = 145 (a₀ = 1,2·10⁻¹⁰), 156 (cH₀/2π) und 166 bei z = 1090, nie 1; r* = λ̄ gilt nur bei M = 6,1·10⁻²⁹ kg. Die Energieskalen unten und Carters Beziehung (ohne H) bleiben davon unberührt.

Drei Energieskalen fallen bei z ≈ 1100 zusammen:

| Skala | Herkunft | Wert |
|---|---|---|
| m_p c²·α_G^{1/4} | Gravitation + Nukleonmasse | 0,2601 eV |
| kT bei z = 1100 | gemessene Rekombination | 0,2586 eV |
| α³·m_e c² | Elektromagnetismus + Elektronmasse | 0,1986 eV |

Nur unter der Annahme a₀(z) = cH(z)/2π fände die Kreuzung „Neigung des Protons an seiner Compton-Grenze = Schwelle a₀(z)" bei **kT ≈ m_p c²·α_G^{1/4}** statt. [korr.] Der Grund: Nahe der Gleichheit von Materie und Strahlung gilt H ~ T²/m_Pl. Aus H = G·m_p³/ħ² folgt dann T² ~ m_p³/m_Pl, also T ~ m_p·(m_p/m_Pl)^{1/2} = m_p·α_G^{1/4}. In Längenform lautet die Kreuzung r*_p = λ̄_p. Die Bedingung H/2π = α_G·ω_p vergleicht Umläufe pro Sekunde (links) mit einer Rate pro Radiant (rechts); bei einheitlicher Zählung läge sie bei z = 349 (pro Radiant) bzw. z = 104 (pro Umlauf). Die Nähe zur Rekombination hängt also an der Zählung. Carters Beziehung enthält kein H und bleibt davon unberührt.

Die Rekombination findet bei **kT ≈ α³·m_e c²** statt. Nach Saha ist kT_rec = 13,6 eV/l_S, und der Logarithmus l_S ≈ 43–53 liegt nahe 1/(2α) = 69. [korr.] Zahlennähe ohne Herleitung, trägt nicht.

**Beide Skalen sind gleich, wenn gilt:**

```
m_p·α_G^{1/4} ≈ α³·m_e      ⟺      α_G ≈ α¹²·(m_e/m_p)⁴
```

Das ist **Carters Koinzidenz** (Carter 1974) [G]. Sie ist eine der bekanntesten Beziehungen zwischen den Naturkonstanten: Genau durch sie liegen typische Sterne an der Grenze zwischen strahlungs- und konvektionsgetriebenen Hüllen.

Die Zahlen:
- α_G = 5,906·10⁻³⁹
- α¹²(m_e/m_p)⁴ = 2,006·10⁻³⁹
- Verhältnis 2,94, in der Energie also Faktor 2,94^{1/4} = 1,31

## Was das bedeutet

1. **Der Treffer beim Proton ist kein isolierter Zufall.** Er ist Carters Koinzidenz in anderer Gestalt, eine Balance zwischen Gravitation (α_G) und Elektromagnetismus (α) über das Massenverhältnis Elektron/Proton. Das gilt für die Energieskala m_p c²·α_G^{1/4} ≈ kT_rec; als Kreuzung r*_p = λ̄_p trägt der Treffer nicht. [korr.]

2. **Unter der Annahme a₀ ∝ H wäre die Kreuzung eine Eigenschaft des Wasserstoffs.** [korr.] Carters Beziehung betrifft allein Elektron und Proton. Für Helium-Kerne greift die Bedingung nicht: Helium rekombiniert bei Z²·13,6 eV, und die Kreuzung wandert mit m^{3/2}.

3. **Warum der Treffer so genau ist (1 %):** Zwei Abweichungen von je rund 30 % heben sich fast auf. [korr.] Zahlennähe ohne Herleitung, trägt nicht.
   - m_p α_G^{1/4} liegt um den Faktor 1,310 über α³m_e (Carter).
   - kT_rec liegt um den Faktor 1,302 über α³m_e (Saha-Logarithmus gegen 1/2α).

   Das ist teils Glück. Die Größenordnung ist aber durch Carter festgelegt.

## Die Lesart im FTM

[korr.] Diese Lesart setzt die Kreuzung r*_p = λ̄_p bei der Rekombination voraus, also a₀(z) = cH(z)/2π; das trägt nicht (r*_p/λ̄_p ≥ 145, bei der Rekombination 166, nie 1). Die Lesart entfällt damit.

Im Zustandsraum des Modells ist Elektromagnetismus **tangential** (schnelle Zyklen, Phase) und Gravitation **parametrisch** (langsame Modulation). Nach der Torus-Dynamik („Dichte als Taktgeber") dominiert der Elektromagnetismus, wo Struktur fein, dicht und schnell ist, die Gravitation, wo sie grob, dünn und langsam ist. Carters Koinzidenz beschreibt genau den **Balancepunkt beider Rollen**:
- in Sternen: Strahlung gegen Konvektion,
- kosmisch: die Epoche, in der der Wasserstoff aus dem Strahlungsfluid entlassen wird **und** der Faltungstrichter des Protons an seine Quantengrenze stößt.

> **In FTM-Sprache (nur unter der Annahme a₀ ∝ H; entfällt) [korr.]:** Die Rekombination wäre der Moment, in dem sich tangentiale (EM) und parametrische (Faltung) Rolle auf der Skala des einzelnen Protons die Waage halten.

## Status

| Aussage | Stand |
|---|---|
| Kreuzung ħH/2π = G m_p²/λ̄_p bei z = 1119 (Quantenlänge ħ/(m_p c); mit h/(m_p c) bei z = 104) ⟺ m_p α_G^{1/4} ≈ kT_rec | gerechnet unter a₀(z) = cH(z)/2π; das trägt nicht, r*_p/λ̄_p ≥ 145 zu jeder Zeit [korr.] |
| ⟺ Carters Koinzidenz α_G ≈ α¹²(m_e/m_p)⁴ (Faktor 2,9) | bekannte Beziehung; im Modell sind die Exponenten 12 und 4 aus der Sternkette hergeleitet, offen ist 𝒞 = 2,94 [korr.] |
| Deutung als Balance tangential/parametrisch | setzt die Kreuzung voraus, entfällt [korr.] |

## Nächster Schritt

Die Exponenten 12 und 4 sind im Modell hergeleitet: Druckgleichgewicht, Zündung der Kernfusion, Strahlungstransport und Abstrahlung an der Wasserstoffschwelle ergeben zusammen α_G ∝ α¹²(m_e/m_p)⁴. [korr.] Ein Einrasten einzelner Zyklen erzeugt Carters Beziehung nicht; die Bedingung „Kopplungsraten gleich" aus der Torus-Dynamik (Klasse A/B gegen Klasse C) liefert sie also nicht. [korr.] Offen ist nur der Zahlenfaktor 𝒞 = 2,94. [korr.]
