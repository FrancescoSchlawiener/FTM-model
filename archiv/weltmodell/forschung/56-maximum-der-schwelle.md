# Maximum der Schwelle bei z ≈ 0,4

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Die Frage „a₀(z) mit Maximum bei z ≈ 0,4“ ist durch das feste Λ überholt (a₀ wandert nicht, seit z ≈ 10 fest); das Ergebnis „keine Größe trifft“ passt zum Modell. M_Q = 2,4·10¹⁴ M☉ gilt für a₀ = cH₀/2π → mit a₀ = 1,2·10⁻¹⁰ 3,7·10¹⁴ M☉ (M_Q ∝ a₀³). „Kapazitätsgrenze R_Q“ → Radius R_Q, an dem das lokale q das kosmische erreicht. Steigung s₀ → γ₀; Ω_Λr → Ω_DE,r. Die Notiz rechnet mit a₀(0) = cH₀/2π = 1,04·10⁻¹⁰ (Θ₀ = 1).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/maximum_schwelle.py` (Ausgabe `maximum_schwelle.txt`)
> Eingaben wie im Modell: H₀ = 67,4 km/s/Mpc, Ω_m = Ω_b(1 + Q) = 0,3134 (Q = 5,364), Ω_r = 9,14·10⁻⁵, k_s = π²/45, a₀(0) = cH₀/2π (Θ₀ = 1). Werte mit [G] aus dem Gedächtnis.

## Verlangte Form

[korr.] Diese Rechnung setzt eine mit der Zeit wandernde Schwelle a₀(z) voraus, die die ganze Zeitabhängigkeit der Dunklen Energie im Vakuumteil trägt; das trägt nicht: Λ ist heute fest, a₀ wandert nicht (seit z ≈ 10 fest, bei der Rekombination 0,76·a₀). Der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie und konstant; die übrigen zwei Drittel werden in der Rückstellung als Bindungsrate dF/d ln a gesucht. Die neue innere Faltung des Kreislaufs (Λ̇/Λ ≈ 10⁻²⁹ je Jahr) gibt heute w_s = −1 + 1·10⁻¹⁹, unmessbar.

Trägt der Vakuumteil die ganze Zeitabhängigkeit, ist a₀(z)/a₀(0) = √(ρ_s/k_s) mit ρ_s = ρ_DE,CPL − (Ω_DE − k_s), CPL-Werte DESI DR2 + CMB + SN [G]:

| | z_max | max/heute | a₀(1)/max | a₀(2)/a₀(0) | γ₀ = d ln a₀/d ln a |
|---|---|---|---|---|---|
| Pantheon+ | 0,35 | 1,11 | 0,87 | 0,51 | −0,76 |
| DESY5 | 0,40 | 1,18 | 0,87 | 0,44 | −1,16 |
| Union3 | 0,44 | 1,26 | 0,87 | 0,39 | −1,56 |

Die Werte bei z = 2 sind Fortsetzung der CPL-Form.

## Rechnung

Jede Größe X(z) wird als a₀(z) = a₀(0)·X(z)/X(0) eingesetzt; mit ρ_s = k_s(X/X₀)² und E² = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_DE,r + k_s(X/X₀)² werden H(z), das Wachstum D(z) und X zusammen iteriert. Die Lage des Maximums hängt nicht davon ab, mit welcher Potenz X in a₀ eingeht. Größen, die Raten sind, koppeln über Ḣ zurück; dort findet die Iteration keine Lösung mit a₀ > 0, angegeben ist ein Schritt im ΛCDM-Hintergrund (*). Epochen: ρ_m = ρ_DE bei z = 0,30, ä = 0 bei z = 0,64, ρ_m = ρ_s (Vakuumteil fest) bei z = −0,11.

| Größe | z_max | max/heute | X(1)/max | γ₀ | w₀ / w_a | Befund |
|---|---|---|---|---|---|---|
| Gewicht 1 − w (θ = 3H, θ_c = a₀/c), a₀ = cH/2π | – | – | – | 0 | – | konstant 1/(1 + 36π²) = 2,8·10⁻³ |
| 1 − w bei festem a₀ | Zukunft | – | – | +0,65 | −1,23 / +0,30 | wächst monoton |
| Rate d(1 − w)/dt, a₀ fest | 0,64 | 1,48 | 0,95 | −1,59 | −0,55 / −1,15 * | Maximum bei ä = 0 |
| Anteil jenseits Umkehr F(>M_Q) | Zukunft | – | – | +1,05 | −1,28 / +0,41 | wächst monoton |
| Rate dF/d ln a | 0,35 | 1,11 | 0,80 | −0,68 | −0,86 / −0,61 * | Lage hängt an Konventionen |
| Rate dF/dt | 1,01 | 1,59 | 1,00 | −1,15 | −0,54 / −1,09 * | zu früh |
| lineares Wachstum dD/d ln a = fD | 0,44 | 1,10 | 0,92 | −0,52 | −0,88 / −0,43 | Lage passt; zu breit, zu flach |
| lineares Wachstum dD/dt | früh | – | – | −1,09 | −0,75 / +0,11 | kein spätes Maximum |
| Ω_m·Ω_DE (Gleichstand) | 0,35 | 1,16 | 0,61 | −0,90 | −0,88 / −0,58 | Lage, Steigung passen; zu schmal |
| Ω_m·Ω_s (Vakuumteil fest) | 0,32 | 1,09 | 0,75 | −0,59 | −0,89 / −0,53 | Lage passt; Form nicht |
| mitbewegter Hubble-Radius c/(aH) | 0,68 | 1,12 | 0,99 | −0,43 | −0,89 / −0,30 | Maximum bei ä = 0 |
| d(H⁻²)/dt ∝ ε/H | 0,64 | 1,48 | 0,95 | −1,59 | −0,56 / −1,14 * | Maximum bei ä = 0 |
| Zerfall des Potentials (1 − f)D/a | Zukunft | – | – | +1,13 | −1,32 / +0,51 | wächst noch |
| Rückstellung ΣΔs(<c/H)/(c/H) = Ω_b(z)/6 | 16 | – | – | −1,43 | −0,63 / −0,13 | fällt seit dem Materiezeitalter |
| Rückstellung ΣΔs(<R_P)/R_P (heute 0,084) | 1,97 | 1,17 | 0,99 | −0,37 | −0,92 / −0,13 | zu früh, flach |
| ΣΔs(<c/H) absolut ∝ (aH)⁻³ | 0,76 | 1,33 | 0,98 | −0,96 | −0,70 / −0,78 | zu früh, zu breit |
| Bindungsteil ⟨Q_s⟩, lineares Feld | 0,90 | 1,58 | 1,00 | – | – | Betrag 3–5·10⁻⁶, negativ |
| Bindungsteil ⟨Q_s⟩, Trichter bis R_Q | 0,33 | 1,15 | 0,73 | – | – | Betrag 5·10⁻⁷, negativ |

- **M_Q** ist die Masse, bei der der Radius R_Q = r*/(−ln q), an dem das lokale q das kosmische erreicht, den Radius R₂₀₀ erreicht [korr.] (heute 2,4·10¹⁴ M☉ bei a₀ = cH₀/2π; mit a₀ = 1,2·10⁻¹⁰ 3,7·10¹⁴ M☉, M_Q ∝ a₀³ [korr.]). Die Lage des Maximums von dF/d ln a wandert mit den Konventionen: Schwellwert Umkehr oder Kollaps (1,062 bzw. 1,686), Überdichte 100 / 200 / 500, Schwelle atmend oder fest geben z_max zwischen −0,35 und 5,9. Der Treffer bei 0,35 gehört zu einer Wahl.
- **Bindungsteil:** Q_s geht mit negativem Vorzeichen in ρ_s ein; an seinem Maximum ist ρ_s am kleinsten, nicht am größten.

## Ergebnis

- Keine Größe des Modells trifft die verlangte Form. Die Lage z ≈ 0,35–0,44 fällt zwischen den Gleichstand ρ_m = ρ_DE (z = 0,30) und den Übergang ä = 0 (z = 0,64); was mit dem Materieanteil wächst und mit der Beschleunigung einfriert, kehrt dort um. Drei konventionsfreie Größen liegen im Fenster (fD bei 0,44, Ω_m·Ω_DE bei 0,35, Ω_m·Ω_s bei 0,32); keine erfüllt zugleich Breite (a₀(1)/max = 0,87) und heutige Steigung (−0,76 bis −1,56).
- Gewicht der Schließungsregel und Rückstellbilanz haben im Fenster kein Maximum; das Gewicht ist bei atmender Schwelle konstant, ΣΔs/R fällt seit dem Materiezeitalter.
- Für keine der Größen gibt es eine Gleichung des Modells, die sie an a₀ bindet. Das Maximum ist damit nicht hergeleitet.
- **Gegen das Modell geprüft:** kein Widerspruch. „Aus einer Horizonttemperatur allein folgt das nicht; es fehlt eine Größe, die im späten Materiezeitalter wächst und bei z ≈ 0,4 umkehrt“ gilt, erweitert um Gewicht, Strukturanteil, Wachstum, Gleichstand und Rückstellung. Der Wert 1 − w = 3,7·10⁻³ heute gilt mit a₀ = 1,2·10⁻¹⁰; mit a₀ = cH₀/2π ist er 2,8·10⁻³.
