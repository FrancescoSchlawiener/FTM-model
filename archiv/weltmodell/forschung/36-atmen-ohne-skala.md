# Das Atmen als Zustandsgleichung w(z) ohne gesetzte Skala

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Die Notiz rechnet mit a₀ = Θ·cH/2π, bei Θ = 1 heute a₀ = cH₀/2π = 1,04·10⁻¹⁰ m/s². „Es fehlt das Zeitgesetz der Schwelle … muss aus dem Temperaturgesetz des Horizonts kommen“ (Kurzfassung, Ergebnis 4–5) → trägt nicht: Λ ist heute fest, a₀ wandert nicht, Horizonte als Deckel tragen nicht (eigene Zahlen: ΔN_eff = 2,09; a₀(0) = 0,81·10⁻¹⁰); die übrigen zwei Drittel der Dunklen Energie werden als Bindungsrate dF/d ln a gesucht. „F_coll = 1 … eine Massenskala kommt nicht vor“ (Kurzfassung, Abschnitt 5, Ergebnis 2) → gilt nur für den Bestand; als Rate dF/d ln a hängt die Rückstellung an der kleinsten gebundenen Masse (für ein Maximum bei z ≈ 0,4 bei 10¹³ M☉ oder darüber). Atmen a₀ ∝ H in Abschnitt 4 ohne Rückkopplung → selbstkonsistent wie Abschnitt 3.2: Steigung −0,47 → −0,60, a₀(1)/a₀(0) 1,79 → 1,95, Tabellenzeile „Atmen a₀ ∝ H“ mit dem Kern neu (z = 0,25: 1,14 → 1,18). [G] an σ₈, n_s, Ω_b h². ρ_F, ρ_B als gelesene Dichten gekennzeichnet. Umbenannt: (∂χ)² → (∂ξ)²; Ω_Λr → Ω_DE,r, ρ_Λr → ρ_DE,r; Steigung s₀ → m₀; ν = δ_c/σ → p_h; a₀ = c²/L mit L als Umfang → a₀ = c²/ℓ₀.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: rechnungen/atmen_ohne_skala.py
> Schwelle: Die Formeln benutzen a₀ = Θ·cH/2π mit dem Modellwert bei Θ = 1, H₀ = 67,4 km/s/Mpc: a₀ = 1,04·10⁻¹⁰ m/s². Der gemessene Wert 1,2·10⁻¹⁰ m/s² erscheint nur im Vergleich; mit ihm wäre Ω_s = 0,291 statt π²/45 = 0,219.
> Kosmische Eingaben: Ω_b h² = 0,02237 [G], Q = 5,364, also Ω_m = Ω_b(1+Q) = 0,3134, Ω_r = 9,14·10⁻⁵, σ₈ = 0,811 [G], n_s = 0,965 [G]. Werte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.

---

## Kurzfassung

| Frage | Ergebnis |
|---|---|
| Was vom Strukturterm wirkt im Hintergrund? | Nur sein Vakuumteil: ρ_s = (π²/45)·Θ²·ρ_krit(z) mit Θ = 2πa₀/(cH). Der ν-Teil und die Gedächtnismode der Schließungsregel wirken zusammen exakt staubartig. Der Bindungsteil ⟨Q_s⟩ beträgt 5·10⁻⁷ (Trichter) bzw. 5·10⁻⁵ (lineares Feld) des Vakuumteils. |
| Zustandsgleichung | w_s = −1 − (2/3)·d ln a₀/d ln a. Das ganze w(z) des Strukturterms ist das Zeitgesetz der Schwelle. |
| Rückkopplung H ↔ a₀ ↔ ρ_s | H² = H_R²/(1 − Ω_sΘ²). Die Abbildung ist eine Kontraktion (Faktor ≤ 0,219), der Fixpunkt ist eindeutig. Aus G, c und H lässt sich nur die Hubble-Masse c³/(2GH) = 4,7·10²² M☉ bilden. Der Fixpunkt legt Ω_s(z) fest, **keine Massenskala**. |
| Schließungsregel und Skala | Die Gedächtnismode wirkt mit ihrem Gewicht von Anfang an. Die wirksame Summe pro mitbewegtem Volumen ist daher erhalten; der Anteil F_coll(>M) und mit ihm die Massenskala fallen heraus. Das gilt nur für den Bestand; als Rate dF/d ln a gelesen hängt die Rückstellung an der kleinsten gebundenen Masse (für ein Maximum bei z ≈ 0,4 bei 10¹³ M☉ oder darüber). [korr.] |
| Θ = 1 (Temperatur ħH/2πk) | Ω_s = 0,219 zu jeder Zeit, w_s = w_ges; w₀ = −0,87, w_a > 0; im Strahlungszeitalter wie ΔN_eff = 2,09. |
| Θ = 1 − ε/2 (Temperatur des dynamischen Horizonts) | Eindeutige Lösung: Θ₀ = 0,777, Ω_s,0 = 0,132, a₀(0) = 0,81·10⁻¹⁰ m/s². Im Strahlungszeitalter null, bei z = 1090 0,8 %. w₀ = −1,02; CPL (0 < z < 2): w₀ = −1,07, w_a = +0,24. Minimum von ρ_s bei z = 0,73. |
| a₀ = c²/ℓ₀ mit festem mitbewegtem ℓ₀ [korr.] | w_s = −1/3 exakt; w₀ = −0,79, CPL w₀ = −0,82, w_a = +0,51. |
| DESI-Quadrant (w₀ > −1, w_a < 0) | Keine der skalenfreien Formen liegt darin. |
| Was DESI von a₀(z) verlangt | a₀(z)/a₀(0) mit Maximum 1,10–1,26 bei z = 0,35–0,44, bei z = 1 gleich 0,97–1,09, bei z = 2 gleich 0,39–0,51; heutige Steigung d ln a₀/d ln a = −0,76 bis −1,56. Das Atmen a₀ ∝ H liefert selbstkonsistent (Abschnitt 3.2) −0,60 und bei z = 1 den Wert 1,95. [korr.] |
| Die fehlende Zahl | Keine Masse im Fixpunkt. Trüge a₀(z) die Form, fehlten die heutige Steigung m₀ = d ln a₀/d ln a (für w₀) und für w_a < 0 ein Maximum von a₀(z) bei z_c ≈ 0,4. [korr.] Das trägt nicht: Λ ist heute fest, a₀ wandert nicht (seit z ≈ 10 fest, bei der Rekombination 0,76·a₀); Horizonte als Deckel tragen nicht, wie die Zahlen dieser Notiz zeigen (Θ = 1: ΔN_eff = 2,09, gemessen N_eff = 2,99 ± 0,17; dynamischer Horizont: a₀(0) = 0,81·10⁻¹⁰ m/s²). Der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie; die übrigen zwei Drittel werden in der Rückstellung als Bindungsrate dF/d ln a gesucht. Die Beschleunigungsrelation bei z = 0,5–1 prüft, dass a₀ nicht wandert. [korr.] |

---

## 1. Zustandsgleichung des Vakuumteils

Die Wirkung enthält im Strukturterm (2/ℓ₀²)·Q̃(q) mit ℓ₀ = c²/a₀. Im feldfreien Bereich (y = 0, q = 1) ist Q̃ = −4π⁴/15. Nach der Variation nach der Metrik steht dort ½g_αβL_s, und das wirkt wie −Λ_s g_αβ:

$$
\rho_s c^2=\frac{4\pi^4}{15}\,\frac{a_0^2}{8\pi G}.
$$

Schreibt man a₀ = Θ·cH/2π, dann folgt

$$
\frac{8\pi G}{3}\rho_s=\frac{\pi^2}{45}\,\Theta^2H^2,\qquad \Omega_s(z)=\frac{\pi^2}{45}\,\Theta(z)^2 .
$$

Bei Θ = 1 ist das der Wert π²/45 = 0,2193 (Probe: 0,219325). Der Vakuumteil ist eine Eigenschaft der Kinetikfunktion des Strukturterms, also der Geometrie, und hängt allein an a₀. Materie ist getrennt erhalten, die Bianchi-Identität verlangt daher ρ̇_s + 3H(1 + w_s)ρ_s = 0. Mit ρ_s ∝ a₀² folgt

$$
\boxed{\,w_s=-1-\frac{2}{3}\,\frac{d\ln a_0}{d\ln a}\,}
$$

Fällt a₀, ist w_s > −1. Steigt a₀, ist w_s < −1. Die Zustandsgleichung des Strukturterms ist also das Zeitgesetz der Schwelle, sonst nichts. Abschnitt 2 zeigt, dass die übrigen Teile des Strukturterms daran nichts ändern.

---

## 2. Bilanz des Strukturterms im Hintergrund

### 2.1 ν-Teil und Gedächtnismode

Gelesen von einem Beobachter, der ungeänderte Regeln annimmt (Poisson-Inversion), gilt nach der Schließungsregel ρ̇_F + 3Hρ_F = −Γρ_F und ρ̇_B + 3Hρ_B = +Γρ_F, also ρ_b + ρ_F + ρ_B ∝ a⁻³ exakt. ρ_F und ρ_B sind gelesene Dichten, kein Stoff: die Umrechnung der Mode ξ (frei) und der Regel (gebunden) in Dichten; die Bilanz ist eine Buchführung dieser Lesung. [korr.] Die freie Gedächtnismode ist ein eingeschränkter Freiheitsgrad der Metrik ((∂ξ)² = −1, druckfrei) [korr.], der gebundene Anteil ist die Strukturenergie des Operators. Beide Formen gehen mit demselben Gewicht (1+Q)M_b in die Friedmann- und die Poisson-Gleichung ein.

**Folge für die Zusatzstrecken.** Die Faltungstiefe hängt nur an der Masse, ΣΔs = G·M_dyn/(3c²). An der Kapazitätsgrenze R_Q = r*/(−ln q) gilt ν(R_Q) = 1 + Q. Das ist unabhängig von a₀, denn R_Q ∝ r* ∝ a₀^{−1/2} wandert mit:

| M_b [M☉] | a₀ | r* [kpc] | R_Q [kpc] | ν(R_Q) |
|---|---|---|---|---|
| 6·10¹⁰ | 1,04·10⁻¹⁰ | 8,96 | 52,4 | 6,364 |
| 6·10¹⁰ | 2,08·10⁻¹⁰ | 6,33 | 37,1 | 6,364 |
| 10¹⁴ | 1,04·10⁻¹⁰ | 366 | 2140 | 6,364 |
| 10¹⁴ | 2,08·10⁻¹⁰ | 259 | 1513 | 6,364 |

Damit ist die wirksame Summe pro mitbewegtem Volumen, G ρ_m a³V_c/(3c²) mit ρ_m = (1+Q)ρ_b, konstant. Die Faltungstiefen ΣΔs selbst stammen nur von den Baryonen; der Faktor 1+Q ist die Verstärkung der Regel. Es gibt keinen Faktor „Anteil der kollabierten Masse“. Die Gedächtnismode wirkt schon vor dem Kollaps mit, beim Kollaps wechselt nur die Form, nicht die wirksame Summe.

### 2.2 Bindungsteil

Zum Vakuumteil kommt in jedem Feld der Beitrag −(a₀²/8πG)·Q_s(y) mit Q_s = 4∫₀^s t³/(eᵗ−1) dt, s = √(g/a₀). Gerechnet wird mit a₀(z) = cH(z)/2π, damit r* = √(GM_b/a₀) mit dem Atmen wächst. Jeder Trichter zählt als Punktmasse M_b = M/(1+Q) bis R_Q bzw. R₂₀₀; die Häufigkeit folgt aus Press-Schechter mit BBKS-Transfer. Das Volumenintegral eines Trichters ist r*³·J(R/r*) mit J(1) = 11,5, J(5,85) = 36,3 und J(100) = 82,8 (Schwanz (16π/3) ln u):

| z | ⟨Q_s⟩/(4π⁴/15) bis R_Q | bis R₂₀₀ | M_nl(z) [M☉/h] |
|---|---|---|---|
| 0 | 5,2·10⁻⁷ | 5,6·10⁻⁷ | 3,9·10¹² |
| 0,5 | 5,8·10⁻⁷ | 6,6·10⁻⁷ | 6,6·10¹¹ |
| 1 | 4,4·10⁻⁷ | 5,2·10⁻⁷ | 1,0·10¹¹ |
| 2 | 1,9·10⁻⁷ | 2,5·10⁻⁷ | 3,0·10⁹ |

Hinzu kommt das lineare Feld zwischen den Trichtern: σ_g = 1,0·10⁻¹² m/s² = 0,0098·a₀ heute, also ⟨Q_s⟩/(4π⁴/15) = 4,7·10⁻⁵ (z = 0) und 7,4·10⁻⁵ (z = 1).

Der Bindungsteil verschiebt ρ_s um weniger als 10⁻⁴ und w entsprechend. Der Grund ist, dass die Strukturenergie eines Trichters von der Ordnung M·v²/c² ist.

**Ergebnis von Abschnitt 2.** Im Hintergrund wirkt vom Strukturterm allein der Vakuumteil. w(z) folgt vollständig aus a₀(z).

---

## 3. Die Rückkopplung H ↔ a₀ ↔ ρ_s

### 3.1 Fixpunkt

Friedmann mit dem Vakuumteil, E = H/H₀:

$$
E^2=R(a)+\Omega_s\Theta^2E^2,\qquad R(a)=\Omega_m a^{-3}+\Omega_r a^{-4}+\Omega_{DE,r}.
$$

Die Abbildung T(E²) = R + Ω_sΘ²E² hat die Ableitung Ω_sΘ² ≤ 0,219 < 1. Sie ist eine Kontraktion, der Fixpunkt E² = R/(1 − Ω_sΘ²) ist eindeutig und stabil. Ω_DE,r [korr.] ist der Rest der Dunklen Energie, den der Strukturterm nicht trägt.

**Keine Skala aus dem Fixpunkt.** In a₀ = c·kT/ħ mit T ∝ ħH/k fallen ħ und k heraus. Aus G, c und H lässt sich genau eine Masse bilden:

M_H = c³/(2GH) = 4,65·10²² M☉, σ(M_H) = 0,11, F_coll(>M_H) = 8·10⁻⁵³.

Auch die Schwelle schneidet keine Masse heraus. Am Umkehrradius (ρ = 5,55 ρ̄_m) gilt für die baryonische Neigung g_N/a₀ = (9π²/16)(Ω_b/2)·H²R/(cH/2π) = 17,4·Ω_b·HR/c. Das ergibt 4,9·10⁻⁵ bei 10¹⁰ M☉/h und 2,3·10⁻³ bei 10¹⁵ M☉/h: Jede gebundene Struktur liegt an ihrem Rand tief unter der Schwelle.

Die Selbstkonsistenz legt also Ω_s(z) fest. Eine Massenskala folgt aus ihr nicht, weil in der Schleife keine vorkommt.

### 3.2 Θ = 1: Temperatur ħH/2πk zu jeder Zeit

Ω_s(z) = 0,2193 ist dann konstant. Der Vakuumteil läuft mit der Gesamtdichte mit (w_s = w_ges) und wirkt in der Friedmann-Gleichung wie G_kos/G = 1/(1 − Ω_s) = 1,281.

| | w(0) | w(0,5) | w(1) | w(2) | w_a = −dw/da bei a = 1 | CPL-Fit 0 < z < 2 |
|---|---|---|---|---|---|---|
| Θ = 1 | −0,872 | −0,668 | −0,459 | −0,201 | +0,17 | w₀ = −0,98, w_a = +1,07 |

Im Strahlungszeitalter entspricht der konstante Anteil ΔN_eff = 2,09, bei der Rekombination 21,9 % des Gesamtbudgets. Der Quadrant ist w₀ > −1 mit w_a > 0.

### 3.3 Θ = 1 − ε/2: Temperatur des dynamischen Horizonts

Die Horizonttemperatur ħH/2πk gilt exakt für den de-Sitter-Horizont. Für einen expandierenden Raum mit ε = −Ḣ/H² hat der scheinbare Horizont die Temperatur

$$
T=\frac{\hbar H}{2\pi k}\Big(1-\frac{\varepsilon}{2}\Big)=\frac{\hbar\,\mathcal R}{24\pi k\,H},\qquad \mathcal R=6(\dot H+2H^2),
$$

also

$$
a_0=\frac{c\,\mathcal R}{24\pi H},\qquad \Theta=1-\frac{\varepsilon}{2}=\frac{\rho-3p}{4\rho}.
$$

Θ ist die Spur des Energie-Impuls-Tensors pro 4ρ: 0 für Strahlung, ¼ für Staub, 1 für Vakuum. Im de-Sitter-Grenzfall ist das wieder a₀ = cH/2π.

**Die Rückkopplung wird differentiell**, denn ε enthält Ḣ und damit ρ_s selbst. Mit v = −ln(1 − Ω_sΘ²) und Θ_R = 1 − ε_R/2 (aus Materie, Strahlung und Ω_DE,r [korr.] allein) gilt

$$
v'=4\big(\Theta(v)-\Theta_R(a)\big),\qquad \Theta(v)=\sqrt{(1-e^{-v})/\Omega_s},\qquad '=d/d\ln a .
$$

Abweichungen von der langsamen Lösung wachsen vorwärts mit der Rate λ = 2(1 − Ω_sΘ²)/(Ω_sΘ) je e-Faltung: 10,2 heute, 21,8 bei z = 1, 47,5 bei z = 1090, 10⁷ bei z = 10⁹. Im Strahlungszeitalter ist Θ = 0 erzwungen. Jede andere Lösung bräuchte dort eine Anfangsabweichung der Größe e^{−∫λ dN}. Übrig bleibt eine Lösung, die langsame:

$$
\Theta=\Theta_R+\frac{\Omega_s\,\Theta\,\Theta'}{2(1-\Omega_s\Theta^2)} .
$$

Die Iteration ergibt für Θ₀ die Folge 0,7762 / 0,7804 / 0,7784 / 0,7770 / 0,7771, also **Θ₀ = 0,777 ± 0,002**. Mit der Schließung E(1) = 1 folgt Ω_DE,r = 0,554 [korr.] und Ω_s,0 = 0,132. Probe: Die exakte Differentialgleichung, von heute aus rückwärts integriert (stabile Richtung), trifft Θ(z=1) = 0,4033 und Θ(z=1090) = 0,1906 wie die langsame Lösung.

| z | 0 | 0,5 | 1 | 2 | 5 | 1090 | 10⁴ | 10⁹ |
|---|---|---|---|---|---|---|---|---|
| Θ | 0,777 | 0,545 | 0,403 | 0,301 | 0,256 | 0,191 | 0,064 | 0,000 |
| Ω_s(z) | 0,132 | 0,065 | 0,036 | 0,020 | 0,014 | 0,0080 | 9·10⁻⁴ | 2·10⁻¹³ |
| a₀(z)/a₀(0) | 1 | 0,92 | 0,93 | 1,17 | 2,75 | 5700 | 9·10⁴ | – |
| zum Vergleich H(z)/H₀ | 1 | 1,31 | 1,78 | 3,03 | 8,32 | 23 233 | 1,1·10⁶ | – |

| | w(0) | w(0,5) | w(1) | w(2) | w_a = −dw/da bei a = 1 | CPL-Fit 0 < z < 2 |
|---|---|---|---|---|---|---|
| Θ = 1 − ε/2 | −1,023 | −1,017 | −0,974 | −0,849 | −0,01 | w₀ = −1,07, w_a = +0,24 |

- Der Vakuumteil schaltet sich im Strahlungszeitalter von selbst ab. Bei der Nukleosynthese ist er null, bei der Rekombination beträgt er 0,8 %.
- ρ_s hat ein **Minimum** bei z = 0,73, nahe am Beginn der Beschleunigung (ε = 1 bei z = 0,62). Davor fällt ρ_s wie verdünnter Staub, danach steigt es. In der Zeitrichtung kreuzt w bei z ≈ 0,8 den Wert −1 von oben nach unten.
- Der Quadrant ist w₀ < −1 mit w_a > 0, also das Spiegelbild des DESI-Befunds.
- Die Schwelle heute beträgt a₀ = 0,777·cH₀/2π = 0,81·10⁻¹⁰ m/s². Bis z = 1 bleibt sie nahezu konstant (0,92–0,93).

**Die Schar der übrigen Lösungen.** Θ₀ ist die Integrationskonstante der Differentialgleichung. Rückwärts laufen alle Lösungen innerhalb von ~0,1 e-Faltungen auf die langsame zusammen:

| Θ₀ | Ω_DE,r [korr.] | ε₀ | w(0) | w(0,5) | CPL w₀ / w_a |
|---|---|---|---|---|---|
| 0,60 | 0,608 | +0,80 | −0,685 | −1,015 | −0,94 / −0,03 |
| 0,70 | 0,579 | +0,60 | −0,876 | −1,016 | −1,01 / +0,11 |
| 0,777 (langsam) | 0,554 | +0,45 | −1,024 | −1,017 | −1,07 / +0,24 |
| 1,00 | 0,467 | 0,00 | −1,454 | −1,021 | −1,30 / +0,71 |
| 1,151 (a₀ = 1,2·10⁻¹⁰, H₀ = 67,4) | 0,396 | −0,30 | −1,748 | −1,025 | −1,51 / +1,15 |

Die Schar ändert w nur für z ≲ 0,1. Bei z = 0,5 liegen alle Lösungen bei w = −1,02. Der gemessene Wert a₀ = 1,2·10⁻¹⁰ m/s² verlangt in dieser Lesart ε₀ < 0, also ein heute wachsendes H.

### 3.4 a₀ = c²/ℓ₀ mit festem mitbewegtem ℓ₀ [korr.]

Diese Form wächst mit (1+z) und ist unabhängig von H. Es gilt ρ_s ∝ a⁻², w_s = −1/3 exakt:

| | w(0) | w(0,5) | w(1) | w(2) | CPL-Fit 0 < z < 2 |
|---|---|---|---|---|---|
| Zelle | −0,787 | −0,658 | −0,565 | −0,461 | w₀ = −0,82, w_a = +0,51 |

Der Wert w₀ liegt im DESI-Bereich, das Vorzeichen von w_a nicht.

---

## 4. Was die DESI-Zahlen von a₀(z) verlangen

Verwendet werden die CPL-Werte aus DESI DR2 + CMB + SN **[G]**:

| Kombination | w₀ | w_a |
|---|---|---|
| + Pantheon+ | −0,838 | −0,62 |
| + DESY5 | −0,752 | −0,86 |
| + Union3 | −0,667 | −1,09 |

Trägt der Vakuumteil die ganze Zeitabhängigkeit, dann gilt ρ_s(z) = ρ_DE,CPL(z) − (Ω_DE − Ω_s), und daraus folgt a₀(z)/a₀(0) = √(ρ_s(z)/Ω_s). Negative Einträge bedeuten a₀² < 0:

| a₀(z)/a₀(0) | z = 0,25 | 0,5 | 0,75 | 1 | 1,5 | 2 | 2,5 | Maximum bei | a₀² = 0 bei |
|---|---|---|---|---|---|---|---|---|---|
| Pantheon+ | 1,101 | 1,097 | 1,044 | 0,965 | 0,764 | 0,510 | −0,13 | z = 0,35 | z = 2,47 |
| DESY5 | 1,162 | 1,176 | 1,120 | 1,027 | 0,779 | 0,442 | −0,39 | z = 0,41 | z = 2,27 |
| Union3 | 1,223 | 1,256 | 1,199 | 1,094 | 0,806 | 0,393 | −0,50 | z = 0,44 | z = 2,18 |
| Atmen a₀ ∝ H, selbstkonsistent (Abschnitt 3.2) [korr.] | 1,176 | 1,398 | 1,659 | 1,952 | 2,622 | 3,383 | 4,222 | – | – |
| dynamischer Horizont | – | 0,92 | – | 0,93 | – | 1,17 | – | Minimum z ≈ 0,7 | – |

Mit a₀ = 1,2·10⁻¹⁰ m/s² (Ω_s = 0,291) liegt a₀² = 0 bei z = 2,6–3,2, sonst bleibt das Bild gleich.

Für w₀ allein reicht eine Zahl, die heutige Steigung:

$$
w_0=-1-\frac{2}{3}\,\frac{\Omega_s}{\Omega_{DE}}\,m_0,\qquad m_0=\frac{d\ln a_0}{d\ln a}\Big|_{z=0}.
$$

| | m₀ [korr.] |
|---|---|
| verlangt: Pantheon+ / DESY5 / Union3 | −0,76 / −1,16 / −1,56 |
| Atmen a₀ ∝ H, selbstkonsistent (m₀ = −ε₀) | −0,60 [korr.] |
| Zelle a₀ ∝ 1+z | −1 |
| dynamischer Horizont | +0,18 |

Ergebnis:
1. Bis z ≈ 0,25 stimmen die DESI-Lesart und das Atmen überein (1,16 gegen 1,18). Ab z ≈ 0,5 trennen sie sich, bei z = 1 stehen 1,03 gegen 1,95. [korr.]
2. Die DESI-Form verlangt ein Maximum der Schwelle bei z_c = 0,35–0,44. Keine der drei skalenfreien Formen hat dort ein Maximum. Das Atmen steigt monoton, die Zelle ebenfalls, der dynamische Horizont hat an dieser Stelle ein Minimum.
3. Oberhalb von z ≈ 2,2–2,5 kann der Vakuumteil mit der Amplitude π²/45 die CPL-Extrapolation allein nicht tragen. Dort liegen allerdings kaum noch DESI-Daten; die CPL-Form ist dort eine Extrapolation.

---

## 5. Strukturabhängige Rückstellung und die Rolle der Masse

Die Rückstellung des geschlossenen Torus lautet ΔL_außen = ΣΔs_innen, relativ ΣΔs/R. Als Ansatz für die Dichte wird ρ_R = ρ_DE,r + Ω_s·S(a)/S(1) gesetzt. [korr.] Die Amplitude ist fest, **β = Ω_s/Ω_DE = 0,3195**, und wird nicht angepasst. Zwei Lesarten der Zusatzstrecke einer kollektiven Struktur:
- S1 = F_coll(>M, a)·√(H/H₀)/a: Δs_eff ∝ √(GMa₀) bei fester Länge, R ∝ a.
- S2 = F_coll(>M, a)/a: Δs_eff an der Kapazitätsgrenze R_Q, dort unabhängig von a₀.

**Skalenfreie Wahlen**

| Wahl | F_coll | S1: w₀ / w_a | S2: w₀ / w_a | Maximum von S |
|---|---|---|---|---|
| M = M_nl(z), Kollapsbedingung σ(M, z) = δ_c | erfc(1/√2) = 0,3173 für jedes z | −0,905 / +0,42 | −0,903 / +0,13 | keines (fällt) |
| M → 0: die Gedächtnismode wirkt von Anfang an mit | 1 (Bestand; als Rate dF/d ln a hängt die Rückstellung an der kleinsten gebundenen Masse) [korr.] | −0,905 / +0,42 | −0,903 / +0,13 | keines |
| M = M_H(z) = c³/(2GH) | 0 | – | – | – |

Press-Schechter ist selbstähnlich in p_h = δ_c/σ. [korr.] Jede Skala, die sich selbstkonsistent mit der Kollapsbedingung bewegt, gibt deshalb ein konstantes F_coll und kein Maximum. Mit der Schließungsregel gilt F_coll = 1 (Abschnitt 2.1), und das Ergebnis ist dasselbe. Das gilt für den Bestand; als Rate dF/d ln a gelesen hängt die Rückstellung an der kleinsten gebundenen Masse, für ein Maximum bei z ≈ 0,4 bei 10¹³ M☉ oder darüber. [korr.]

**Feste Masse (Abtastung 10¹⁰–10¹⁴·⁵ M☉/h, Auszug)**

| M [M☉/h] | p_h = δ_c/σ(M) [korr.] | S1: w₀ / w_a / z_max | S2: w₀ / w_a / z_max |
|---|---|---|---|
| 10¹⁰ | 0,47 | −0,859 / −0,12 / 1,94 | −0,885 / −0,24 / 1,12 |
| 10¹¹ | 0,60 | −0,851 / −0,28 / 1,26 | −0,890 / −0,33 / 0,70 |
| 10¹² | 0,81 | −0,865 / −0,44 / 0,65 | −0,922 / −0,37 / 0,32 |
| 10¹³ | 1,17 | −0,972 / −0,33 / 0,13 | −1,030 / −0,20 / 0 |
| 10¹⁴ | 1,81 | −1,204 / +0,22 / 0 | −1,231 / +0,29 / 0 |

Am nächsten an DESY5 liegt S1 bei M = 10¹² M☉/h mit w₀ = −0,865 und w_a = −0,44. Diese Masse ist M_nl bei z = 0,39, also die Masse, die zur Zeit des Maximums gerade kollabiert. Das ist eine Umschreibung von z_c, keine Herleitung. Mit der parameterfreien Amplitude gilt für jede Masse |1 + w₀| ≤ 0,15. Der DESI-Wert w₀ = −0,75 ist damit auch mit gewählter Skala nicht erreichbar.

---

## 6. Ergebnis

1. **Bilanz.** Vom Strukturterm wirkt im Hintergrund allein der Vakuumteil ρ_s = (π²/45)Θ²ρ_krit. Der ν-Teil wirkt mit der Gedächtnismode exakt staubartig. Der Bindungsteil ist mit ≤ 10⁻⁴ vernachlässigbar. Es gilt w_s = −1 − (2/3)·d ln a₀/d ln a. Das Atmen als Zustandsgleichung ist also das Zeitgesetz der Schwelle.
2. **Fixpunkt ohne Skala.** Die Rückkopplung H ↔ a₀ ↔ ρ_s hat einen eindeutigen, stabilen Fixpunkt (Kontraktion mit Faktor ≤ 0,219). Er legt Ω_s(z) fest. In der Schleife kommt keine Masse außer c³/(2GH) vor. Mit der Schließungsregel fällt auch der Anteil der kollabierten Masse heraus, weil die Gedächtnismode von Anfang an mitwirkt. Das gilt für den Bestand: Im Fixpunkt und im Bestand kommt keine Massenskala vor. Als Rate dF/d ln a gelesen hängt die Rückstellung dagegen an der kleinsten gebundenen Masse, für ein Maximum bei z ≈ 0,4 bei 10¹³ M☉ oder darüber. [korr.]
3. **Die drei skalenfreien Zeitgesetze der Schwelle**
   - Θ = 1: w₀ = −0,87, w_a > 0, mit ΔN_eff-Äquivalent 2,09 im Strahlungszeitalter.
   - Dynamischer Horizont, a₀ = cℛ/(24πH): eindeutige Lösung ohne freie Zahl. Θ₀ = 0,777, Ω_s,0 = 0,132, a₀(0) = 0,81·10⁻¹⁰ m/s², null im Strahlungszeitalter, 0,8 % bei der Rekombination; w₀ = −1,02, CPL w₀ = −1,07, w_a = +0,24.
   - Zelle a₀ ∝ 1+z: w₀ = −0,79, w_a = +0,51.

   Den DESI-Quadranten erreicht keines der drei.
4. **Was ein Zeitgesetz der Schwelle verlangte.** Trüge a₀(z) allein die Form der Dunklen Energie, wäre für w₀ die heutige Steigung m₀ = d ln a₀/d ln a nötig: DESY5 verlangt m₀ = −1,16 (Pantheon+ −0,76, Union3 −1,56), das Atmen liefert selbstkonsistent −0,60. Für w_a < 0 bräuchte es zusätzlich ein Maximum von a₀(z) bei z_c = 0,35–0,44 mit a₀,max/a₀(0) = 1,10–1,26. [korr.]
5. **Das trägt nicht.** [korr.] De-Sitter-Temperatur und dynamischer Horizont geben m₀ = −ε₀ = −0,60 bzw. +0,18 und kein Maximum; Θ = 1 gäbe im Strahlungszeitalter ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17), der dynamische Horizont a₀(0) = 0,81·10⁻¹⁰ m/s². Horizonte als Deckel tragen also nicht. Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ wandert nicht (seit z ≈ 10 fest, bei der Rekombination 0,76·a₀). Der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie; die übrigen zwei Drittel werden in der Rückstellung als Bindungsrate dF/d ln a gesucht. Messbar bleibt a₀(z = 1)/a₀(0) in der Beschleunigungsrelation früher Galaxien: 1 bei fester Schwelle, 1,95 beim Atmen, 1,03 bei DESY5, 0,93 beim dynamischen Horizont. [korr.]

---

## Literatur

- G. W. Gibbons, S. W. Hawking, *Cosmological event horizons, thermodynamics, and particle creation*, Phys. Rev. D 15, 2738 (1977) [G]
- S. A. Hayward, *Unified first law of black-hole dynamics and relativistic thermodynamics*, Class. Quantum Grav. 15, 3147 (1998) [G]
- R.-G. Cai, S. P. Kim, *First law of thermodynamics and Friedmann equations of Friedmann-Robertson-Walker universe*, JHEP 02, 050 (2005) [G]
- R.-G. Cai, L.-M. Cao, Y.-P. Hu, *Hawking radiation of apparent horizon in a FRW universe*, Class. Quantum Grav. 26, 155018 (2009) [G]
- M. Chevallier, D. Polarski, Int. J. Mod. Phys. D 10, 213 (2001); E. V. Linder, Phys. Rev. Lett. 90, 091301 (2003) [G]
- DESI Collaboration, *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints*, arXiv:2503.14738 (2025) [G]
- W. H. Press, P. Schechter, ApJ 187, 425 (1974); J. M. Bardeen, J. R. Bond, N. Kaiser, A. S. Szalay, ApJ 304, 15 (1986); N. Sugiyama, ApJS 100, 281 (1995) [G]
- Planck Collaboration, *Planck 2018 results VI*, A&A 641, A6 (2020) [G]
