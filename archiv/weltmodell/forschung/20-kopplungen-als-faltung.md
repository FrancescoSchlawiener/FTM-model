# Die Kopplungskonstanten als Faltungsverhältnisse

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Kreuzung r*_p = λ̄_p bei z = 1119 und a₀ = c·kT_GH/ħ = cH/2π zu jeder Zeit → setzt a₀ ∝ H voraus, trägt nicht; a₀ ≈ cH₀/2π nur heute; Kreuzung, Merksatz und Kette nur noch unter dieser Annahme; r*_p = 155,7·λ̄_p → 145·λ̄_p (a₀ = 1,2·10⁻¹⁰); M < 32,4 MeV/c² → 6,1·10⁻²⁹ kg = 34 MeV/c² (32,4 gehört zu cH₀/2π); Begründung „bleibt es wegen H_Λ“ → a₀ heute fest; „Carter nicht hergeleitet, Bedingung aus der Torus-Dynamik gesucht“ → Exponenten 12 und 4 aus der Sternkette hergeleitet, ein Einrasten einzelner Zyklen erzeugt Carter nicht, offen nur 𝒞 = 2,94; Zahlennähe markiert (Wahl ħ oder h verschiebt z von 1119 auf 104). Umbenennung (ohne Marke): H_Λ → H_DE. Abschnitt 3 rechnete mit a₀ = cH₀/2π = 1,04·10⁻¹⁰.

> Rechnung: `rechnungen/compton_faltung.py`

## 1. Gravitation und Elektromagnetismus in derselben Sprache

| Kopplung | als Längenverhältnis | Wert |
|---|---|---|
| Gravitation | α_G = **3·Δs_p / λ̄_p** | 5,906149·10⁻³⁹ |
| Elektromagnetismus | α = **r_e / λ̄_e** | 0,0072973526 = 1/137,036 |

- **α_G** ist die Faltungstiefe des Protons, gemessen in seiner eigenen Compton-Länge (mal 3). Δs_p = 4,14·10⁻⁵⁵ m, λ̄_p = 2,10·10⁻¹⁶ m.
- **α** ist der elektromagnetische Radius des Elektrons (klassischer Radius r_e = 2,82·10⁻¹⁵ m), gemessen in seiner Compton-Länge λ̄_e = 3,86·10⁻¹³ m.

Beide Naturkonstanten sind damit dasselbe Muster: **Wirkungslänge durch Quantenlänge.** Im FTM ist das die Brücke zwischen Faltung (parametrisch) und Elektromagnetik (tangential).

## 2. Carters Koinzidenz als reine Längenbeziehung

```
α_G ≈ α¹²·(m_e/m_p)⁴      ⟺      3·Δs_p · λ̄_e¹⁶  ≈  r_e¹² · λ̄_p⁵
```

Numerisch: 3,04·10⁻²⁵³ m¹⁷ gegen 1,03·10⁻²⁵³ m¹⁷, Verhältnis 2,94. In der Sprache des Modells verbindet die Beziehung die **Faltungstiefe des Protons** mit dem **elektromagnetischen Radius** und den **Compton-Längen** von Elektron und Proton.

## 3. Die Kreuzung als Gleichheit zweier Energien

[korr.] Diese Rechnung setzt a₀(z) = cH(z)/2π zu jeder Zeit voraus; das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest (bei der Rekombination 0,76·a₀); ≈ cH₀/2π gilt nur heute. Damit ist r*_p/λ̄_p = 145 (a₀ = 1,2·10⁻¹⁰), 156 (cH₀/2π), 166 bei z = 1090, nie 1.

Nur unter dieser Annahme erreicht der Schwellenradius des Protons r*_p = √(G m_p/a₀) seine Compton-Länge λ̄_p, und zwar wenn [korr.]

```
k·T_GH(z) = ħH(z)/2π  =  α_G·m_p c²  =  G·m_p²/λ̄_p        erfüllt bei z = 1119  (Plasma-Temperatur kT = 0,263 eV)
```

- links steht die Temperatur des Hubble-Horizonts (Gibbons-Hawking),
- rechts die Gravitations-Selbstenergie des Protons auf seiner eigenen Compton-Länge.

In Längenform lautet die Kreuzung r*_p = λ̄_p. Beide Seiten der Energieform sind physikalisch: Das 2π der Horizonttemperatur ist die euklidische Periode 2π/H, und λ̄_p = ħ/(m_p c) ist die Reichweite eines Feldes der Masse m_p. Die Epoche hängt an dieser Wahl der Quantenlänge: Mit der Reichweite ħ/(m_p c) liegt sie bei z = 1119, mit der vollen Compton-Wellenlänge h/(m_p c) bei z = 104. [korr.] Zahlennähe ohne Herleitung, trägt nicht. Carters Beziehung enthält kein H und hängt von dieser Wahl nicht ab.

Heute entspricht die Schwelle der Horizonttemperatur, a₀ ≈ c·kT_GH/ħ = cH₀/2π; zu anderen Zeiten gilt das nicht, a₀ ist seit z ≈ 10 fest. [korr.]

> **Die Horizonttemperatur des atmenden Kosmos fällt zur Zeit der Rekombination auf die Gravitations-Selbstenergie des Protons.** [korr.] Zahlennähe ohne Herleitung, trägt nicht; als Kreuzung r*_p = λ̄_p gilt der Satz nur unter a₀ ∝ H.

Heute ist r*_p = 145·λ̄_p (a₀ = 1,2·10⁻¹⁰; 155,7 mit cH₀/2π). [korr.] Teilchen mit M < 6,1·10⁻²⁹ kg = 34 MeV/c² haben r* < λ̄ (32,4 MeV/c² mit cH₀/2π). [korr.] Das Elektron gehört dazu und bleibt es, denn a₀ ist heute fest; die Begründung über die Horizonttemperatur des Endzustands (H_DE) setzt a₀ ∝ H voraus. [korr.] Nur unter dieser Annahme träte die Gleichheit beim Proton ein; bei festem a₀ gilt r* = λ̄ nur bei M = 6,1·10⁻²⁹ kg. [korr.]

## 4. Die Kette

| Stufe | Energie |
|---|---|
| Kreuzung kT_GH = α_G·m_p c² (Plasma-Temperatur; nur unter a₀ ∝ H) [korr.] | 0,2631 eV |
| m_p c²·α_G^{1/4} | 0,2601 eV |
| letzte Streuung (z = 1090) | 0,2562 eV |
| α³·m_e c² | 0,1986 eV |

Unter der Annahme a₀ ∝ H läge die Kreuzung bei der Rekombination, weil [korr.]
- die Atemrate nahe der Materie-Strahlungs-Gleichheit H ~ T²/m_Pl ist, woraus T ~ m_p α_G^{1/4} folgt, und
- m_p α_G^{1/4} ≈ α³m_e gilt. Das ist Carters Koinzidenz.

## 5. Offene Aufgabe

Im FTM lautet die Frage nun: **Warum ist 3Δs_p/λ̄_p ≈ (r_e/λ̄_e)¹²·(λ̄_p/λ̄_e)⁴?**

Die Exponenten 12 und 4 sind im Modell aus der Sternkette hergeleitet: Druckgleichgewicht, Zündung der Kernfusion, Strahlungstransport und Abstrahlung an der Wasserstoffschwelle ergeben zusammen α_G ∝ α¹²(m_e/m_p)⁴. [korr.] Ein Gleichgewicht der Faltungsrate α_G·ω_p (Klasse C) mit elektromagnetischen Raten α^n·ω_e (Klasse A/B), also ein Einrasten einzelner Zyklen, erzeugt Carters Beziehung nicht. [korr.]

Einordnung: Die Identitäten in Abschnitt 1 sind algebraisch exakt, also Umformungen. Neu ist die Lesart. Carters Exponenten sind im Modell hergeleitet; offen ist nur der Zahlenfaktor 𝒞 = 2,94. [korr.]
