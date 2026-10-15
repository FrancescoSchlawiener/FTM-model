# Zeitgesetz der Schwelle aus der Horizonttemperatur

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Voraussetzung, dass a₀(z) als Horizonttemperatur wandert und die Form der Dunklen Energie trägt (Ansatz, Ergebnis 5) → trägt nicht: Λ ist heute fest, a₀ wandert nicht; der Vakuumteil ist ein Drittel der Dunklen Energie, die übrigen zwei Drittel werden als Bindungsrate dF/d ln a gesucht; das Ergebnis der Notiz (kein Horizont trägt) passt zum Modell. Die Notiz rechnet a₀ aus Horizonttemperaturen (Hubble-Radius: a₀ = cH₀/2π = 1,04·10⁻¹⁰ m/s²). Umbenannt: ρ_Λr, Ω_Λr → ρ_DE,r, Ω_DE,r; holographischer Parameter c² → c_h².

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/zeitgesetz_a0.py` (Ausgabe `zeitgesetz_a0.txt`)
> Eingaben: H₀ = 67,4 km/s/Mpc, Ω_m = Ω_b(1+Q) = 0,3134 (Q = 5,364), Ω_r = 9,14·10⁻⁵, Ω_DE = 0,6865, k_s = π²/45 = 0,2193. Kein Wert angepasst. Werte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.

## Ansatz

[korr.] Diese Rechnung setzt voraus, dass a₀(z) als Horizonttemperatur wandert und die Form der Dunklen Energie trägt; das trägt nicht: Λ ist heute fest, a₀ wandert nicht (seit z ≈ 10 fest, bei der Rekombination 0,76·a₀). Der Vakuumteil Ω_s = π²/45 ist ein Drittel der Dunklen Energie; die übrigen zwei Drittel werden in der Rückstellung als Bindungsrate dF/d ln a gesucht. Der Kreislauf verlangt neue innere Faltung mit der Mindestrate Λ̇/Λ ≥ Ṡ_i/(2S_max) ≈ 10⁻²⁹ je Jahr; das gibt heute w_s = −1 + 1·10⁻¹⁹ (im Endzustand 2·10⁻¹⁹), also unmessbar. Das Ergebnis dieser Notiz, dass kein Horizont das gemessene Zeitgesetz trägt, passt dazu.

Die Schwelle ist die Temperatur eines Horizonts als Beschleunigung, a₀ = c·kT/ħ = κ/2π (κ Oberflächenbeschleunigung). Mit Θ = 2πa₀/(cH) und ε = −Ḣ/H²:

```
Ω_s = k_s·Θ²,     H² = (8πG/3)(ρ_m + ρ_r + ρ_DE,r)/(1 − k_sΘ²),     w_s = −1 − (2/3)·d ln a₀/d ln a      [korr.]
```

ρ_DE,r [korr.] ist der Teil der Dunklen Energie, den der Vakuumteil nicht trägt (konstant); er folgt aus H(0) = H₀. Θ hängt über ε oder über Integrale von 1/(aH) wieder von H ab; jedes Zeitgesetz ist daher eine Rückkopplung und wird selbstkonsistent gelöst.

| Horizont | Radius, Temperatur | Θ | d ln a₀/d ln a |
|---|---|---|---|
| Hubble-Radius | c/H | 1 | −ε |
| Ereignishorizont | R_E = a∫_t^∞ c dt′/a | c/(H·R_E) | Θ − 1 |
| scheinbarer Horizont, Kodama-Hayward | c/H, κ = cH(1 − ε/2) | 1 − ε/2 = ℛ/(12H²) | −ε − ε′/(2 − ε) |
| de-Sitter-Endhorizont | c/H_∞ | H_∞/H | 0 |
| Teilchenhorizont (zum Vergleich) | R_P = a∫_0^t c dt′/a | c/(H·R_P) | −(1 + Θ) |

ℛ ist der Ricci-Skalar, ′ = d/d ln a.

## Selbstkonsistenz

Mit ε_R aus Materie, Strahlung und Rest allein und v = −ln(1 − k_sΘ²):

```
Ereignishorizont:   Θ′ = Θ(1 − k_sΘ²)(Θ − 1 + ε_R)
Kodama-Hayward:     v′ = 4(Θ − 1 + ε_R/2)
Teilchenhorizont:   Θ′ = Θ(1 − k_sΘ²)(ε_R − 1 − Θ)
```

- **Hubble-Radius:** Ω_s = 0,219 zu jeder Zeit. Die Rückkopplung wirkt in jeder Epoche wie G_kos = 1,28·G: im Strahlungszeitalter wie ΔN_eff = 2,09, spät wie Ω_m,eff = 0,40 in der Form von H(z). Selbstkonsistent ist a₀(1)/a₀(0) = 1,95 (ohne Rückkopplung 1,79).
- **Ereignishorizont:** Die Gleichung gilt für jede Integrationskonstante C des Integrals; der Ereignishorizont verlangt C = 0. Vorwärts ist Θ = 1 instabil. Unterhalb eines Grenzwerts Θ*₀ läuft Θ auf 0 (dann C > 0, kein Ereignishorizont), oberhalb auf 1/√k_s = 2,14 (H → ∞ in endlicher Zeit). Genau eine Lösung endet in de Sitter: Θ*₀ = 0,863. Das Integral selbst gibt dieselben Θ auf fünf Stellen. Die Form ist die einer holographischen Dunklen Energie mit Parameter c_h² = π²/45 [korr.] neben einem konstanten Rest [G].
- **Kodama-Hayward:** vorwärts instabil mit der Rate 2(1 − k_sΘ²)/(k_sΘ) je e-Faltung. Genau eine Lösung endet in de Sitter; rückwärts aus der Zukunft gerechnet ist Θ₀ = 0,777. Jedes andere Θ₀ wäre eine gesetzte Zahl.
- **de-Sitter-Ende:** H_∞² = H₀²Ω_DE,r/(1 − k_s) [korr.] zusammen mit H(0) = H₀ gibt H_∞ = √Ω_DE·H₀ = 0,829 H₀. ρ_s ist konstant, H(z) ist exakt ΛCDM.
- **Teilchenhorizont:** Im Strahlungszeitalter ist R_P = c/H exakt, also Θ = 1.

## Ergebnisse

„Gesamte DE“ ist Rest plus Vakuumteil; das ist die Größe, die DESI misst. CPL-Fit über 0 < z < 2.

| | Hubble-Radius | Ereignishorizont | Kodama-Hayward | de-Sitter-Ende | Teilchenhorizont |
|---|---|---|---|---|---|
| Θ₀ | 1 | 0,863 | 0,777 | 0,829 | 0,323 |
| a₀(0) [10⁻¹⁰ m/s²] | 1,04 | 0,90 | 0,81 | 0,86 | 0,34 |
| a₀(0) / 1,2·10⁻¹⁰ | 0,87 | 0,75 | 0,68 | 0,72 | 0,28 |
| a₀(z)/a₀(0), z = 0,4 / 1 / 2 | 1,30 / 1,95 / 3,38 | 1,07 / 1,21 / 1,48 | 0,93 / 0,93 / 1,17 | 1 / 1 / 1 | 1,58 / 2,65 / 4,84 |
| Extremum von a₀ | keines | keines | Minimum z = 0,73 | konstant | keines |
| d ln a₀/d ln a heute | −0,60 | −0,14 | +0,19 | 0 | −1,32 |
| w_s bei z = 0 / 1 | −0,60 / −0,16 | −0,91 / −0,72 | −1,13 / −0,84 | −1 / −1 | −0,12 / −0,02 |
| gesamte DE: w(0) | −0,87 | −0,98 | −1,02 | −1 | −0,97 |
| CPL w₀ / w_a | −0,98 / +1,07 | −1,00 / +0,20 | −1,07 / +0,24 | −1 / 0 | −1,06 / +0,58 |
| Ω_s bei z = 1090 | 0,22 | 3·10⁻⁵ | 0,008 | 3·10⁻¹⁰ | 0,13 |
| ΔN_eff-Äquivalent (z = 10⁹) | 2,09 | 0 | 0 | 0 | 2,09 |

## Vergleich mit DESI

DESI DR2 + CMB + SN **[G]**: w₀/w_a = −0,838/−0,62 (Pantheon+), −0,752/−0,86 (DESY5), −0,667/−1,09 (Union3). Trägt der Vakuumteil die ganze Zeitabhängigkeit, liegt das Maximum von a₀ dort, wo ρ_DE maximal ist, bei 1 − a = −(1 + w₀)/w_a, also z = 0,35 / 0,41 / 0,44, unabhängig von Ω_s,0.

Mittlerer Exponent p in a₀ ∝ a^p (p > 0: a₀ wächst mit der Zeit):

| | z = 2 → 1 | 1 → 0,4 | 0,4 → 0 |
|---|---|---|---|
| DESI verlangt (Ω_s,0 = k_s) | +1,6 bis +2,5 | +0,4 | −0,3 bis −0,7 |
| Hubble-Radius | −1,36 | −1,13 | −0,79 |
| Ereignishorizont | −0,51 | −0,34 | −0,20 |
| Kodama-Hayward | −0,59 | +0,02 | +0,21 |
| de-Sitter-Ende | 0 | 0 | 0 |
| Teilchenhorizont | −1,49 | −1,44 | −1,37 |

- **Heutige Steigung:** Mit dem eigenen Ω_s,0 verlangt DESI d ln a₀/d ln a = −0,76 bis −1,56 (Hubble-Radius), −1,02 bis −2,10 (Ereignishorizont), −1,26 bis −2,59 (Kodama-Hayward), −1,11 bis −2,28 (de Sitter). Geliefert werden −0,60, −0,14, +0,19 und 0.
- **Reichweite:** Mit Ω_s,0 < k_s verlangt die CPL-Form schon bei z = 1,75–2,0 ein a₀² = 0 (bei Ω_s,0 = k_s erst bei z = 2,2–2,5).
- **Kreuzungssatz (Ereignishorizont):** Bei Θ = 1 gilt Θ′ = (1 − k_s)ε_R > 0 (0,44 heute, 0,73 bei z = 0,4). Θ kreuzt 1 nur aufwärts; wegen d ln a₀/d ln a = Θ − 1 hat a₀ höchstens ein Minimum, nie ein Maximum, solange ε_R > 0. Hubble-Radius und Teilchenhorizont haben ein festes Vorzeichen (−ε bzw. −(1 + Θ)), de Sitter ist konstant. Kodama-Hayward zeigt das Spiegelbild: ein Minimum bei z = 0,73 statt eines Maximums bei z ≈ 0,4.
- **Schar des Ereignishorizonts:** Die Lösungen oberhalb Θ*₀ enden mit H → ∞. Θ₀ = 1,151 trifft a₀(0) = 1,2·10⁻¹⁰ m/s², mit Ende nach 32 Gyr, CPL −1,07/+0,34 und einem Minimum von a₀ bei z = 0,25. Trägt der Vakuumteil die ganze Dunkle Energie (Ω_DE,r = 0 [korr.]), ist Θ₀ = 1,77, a₀(0) = 1,84·10⁻¹⁰ m/s², w₀ = −1,54, Ende nach 16 Gyr. Alle liegen bei w₀ < −1, w_a > 0.

## Frühe massereiche Scheiben

Zusatzanteil f = 1 − 1/ν = q = e^{−√(g_N/a₀)} für M_b = 10¹¹ M☉ als Punktmasse bei R = 5 kpc (in Klammern 10 kpc), mit dem absoluten a₀(z) jedes Horizonts:

| | z = 0 | z = 1 | z = 2 |
|---|---|---|---|
| Hubble-Radius | 0,10 (0,32) | 0,19 (0,44) | 0,28 (0,53) |
| Ereignishorizont | 0,08 (0,29) | 0,10 (0,32) | 0,13 (0,36) |
| Kodama-Hayward | 0,07 (0,27) | 0,07 (0,26) | 0,09 (0,30) |
| de-Sitter-Ende | 0,08 (0,28) | 0,08 (0,28) | 0,08 (0,28) |

Gemessen **[G]**: Massereiche sternbildende Scheiben bei z ≈ 0,9–2,4 zeigen innerhalb des Halblichtradius f ≲ 0,2, zu höherem z abnehmend (Genzel et al. 2017, 2020). Bei fester Scheibe lässt kein Horizont f mit z fallen; am nächsten kommen de-Sitter-Ende (fest) und Kodama-Hayward (0,07 → 0,09), dann der Ereignishorizont (0,08 → 0,13). Der Hubble-Radius hebt f bei 5 kpc bis z = 1 fast auf das Doppelte, bis z = 2 fast auf das Dreifache. Die höhere Flächendichte früher Scheiben ist nicht gerechnet.

## Ergebnis

1. **Welcher Horizont:** Hubble-Radius und Teilchenhorizont haben im Strahlungszeitalter Θ = 1 und damit ΔN_eff = 2,09; erlaubt sind einige Zehntel **[G]**. Der Hubble-Radius wirkt zusätzlich in jeder Epoche wie G_kos = 1,28·G. Selbstkonsistent bleiben Ereignishorizont, Kodama-Hayward und de-Sitter-Ende; alle drei schalten sich früh ab (Ω_s(1090) ≤ 0,8 %).
2. **Schwelle heute:** Diese drei geben Θ₀ < 1, a₀(0) = 0,81–0,90·10⁻¹⁰ m/s², 25–33 % unter 1,2·10⁻¹⁰. Der Hubble-Radius liegt 13 % darunter. Den Messwert erreicht nur ein Zweig des Ereignishorizonts, der mit H → ∞ endet.
3. **Frühe Scheiben:** a₀(1)/a₀(0) = 1,21 (Ereignishorizont), 0,93 (Kodama-Hayward), 1 (de Sitter) gegen 1,95 beim Hubble-Radius. Die drei selbstkonsistenten Formen vertragen den kleinen Zusatz bei z = 1–2.
4. **DESI:** Keiner liegt im Quadranten w₀ > −1, w_a < 0. Der Ereignishorizont gibt heute w_s = −0,91, aber w_a > 0 und kein Maximum; Kodama-Hayward gibt das Spiegelbild; de Sitter gibt w = −1.
5. **Was fehlt:** DESI verlangt ein a₀, das zwischen z = 2 und z = 0,4 mit der Zeit wächst (p ≈ +2 bzw. +0,4) und danach fällt. Bei allen Horizonten fällt a₀ zwischen z = 2 und z = 1 oder bleibt fest. Beim Ereignishorizont entstünde ein Maximum nur, wenn der übrige Inhalt beim Durchgang Θ = 1 mit der Zeit wächst (ε_R < 0); dann steckte das Zeitgesetz schon im Rest. Aus einer Horizonttemperatur allein folgt das DESI-Zeitgesetz nicht. Es fehlt eine Größe, die im späten Materiezeitalter mit der Zeit wächst und bei z ≈ 0,4 umkehrt. [korr.] Im Modell liegt sie nicht in a₀: Λ ist heute fest, a₀ wandert nicht; die übrigen zwei Drittel der Dunklen Energie werden in der Rückstellung als Bindungsrate dF/d ln a gesucht.

## Literatur

- G. W. Gibbons, S. W. Hawking, Phys. Rev. D 15, 2738 (1977) [G]
- S. A. Hayward, Class. Quantum Grav. 15, 3147 (1998); R.-G. Cai, S. P. Kim, JHEP 02, 050 (2005); R.-G. Cai, L.-M. Cao, Y.-P. Hu, Class. Quantum Grav. 26, 155018 (2009) [G]
- M. Li, *A model of holographic dark energy*, Phys. Lett. B 603, 1 (2004) [G]
- DESI Collaboration, DESI DR2 Results II, arXiv:2503.14738 (2025) [G]
- R. Genzel et al., Nature 543, 397 (2017); R. Genzel et al., ApJ 902, 98 (2020) [G]
- Planck Collaboration, Planck 2018 results VI, A&A 641, A6 (2020) [G]
