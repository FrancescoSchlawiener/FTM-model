# Gravitationswellen und das Starkfeld-Regime

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Umbenennung Massenverhältnis η → η_m, Eichvektor ξ_ν → ε_ν; Ringdown ohne Spin: 194,7 Hz als Wert im Quellsystem gekennzeichnet, im Detektor 178,7 Hz (mit dem Fit 176,1 Hz); zwei Gedächtnisarten (Bahngedächtnis e^{−λ_mΔt} auf gebundenen Bahnen) → ein Gedächtnis (ξ, w, λ_m), kein Nachlauf in Galaxien; Pulsartest in der Feldform (Σ_eff ≈ ν·Σ_b = 78, ν ≈ 1,57) → gilt nur unter der Feldform, die örtliche Regel gibt am Sonnenort 71 M☉/pc²; Scheindichte/Scheinmasse → gelesene Dichte/Masse; Starkfeld „übernommen“ → dτ/dt = √(1 − 6Δs/r) im Modell aus G^t_t = G^r_r hergeleitet; η = 1 allgemein → nur für Kugeln und außerhalb von Scheiben; Haufenrest: Kandidat a₀(z) = cH(z)/2π trägt nicht, im Modell schließt ihn der stoßfreie Rest der Mode ξ; R₀ = 8,2 kpc und v_c = 233 km/s mit [G]. Die Notiz rechnet mit a₀ = 1,2·10⁻¹⁰.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/gravitationswellen.py` (alle Zahlen dieser Notiz).
> Formeln des Modells: Δs = GM/(3c²), g = 3c²Δs/r², dτ/dt = √(1−6Δs/r), r₀ = 6Δs, r* = √(3c²Δs/a₀), ν = 1/(1−q), q = exp(−√(a_N/a₀)) = exp(−r*/r), a₀ ≈ cH₀/2π.
> Zahlenwert der Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Beim Atmen-Takt H₀/2π = a₀/c steht der Modellwert c·H₀/2π = 1,08·10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc).
> Konstanten: GM☉ = 1,32712·10²⁰ m³/s² (IAU), damit Δs☉ = 492,208 m, GM☉/c² = 3Δs☉ = 1 476,6 m, GM☉/c³ = 3Δs☉/c = 4,9255 µs.

**Schreibweise.** In dieser Notiz steht r* für den Schwellenradius, r_T für die Takt-Koordinate (Schildkröten-Koordinate) aus Abschnitt 0.5, j = cJ/(GM²) für den dimensionslosen Spin und a(t) = L(t)/L(t₀) für den Skalenfaktor des Atmens. a_N ist die Newton-Beschleunigung.

---

## Kurzfassung

1. Lokale Störungen des Faltungsmediums gehorchen einer Wellengleichung □h̄_μν = −48π·σ_μν mit der **Faltungsdichte** σ_μν = G·T_μν/(3c⁴). Sie laufen mit c und haben genau zwei transversale Polarisationen. Ihre Quelle ist der zeitlich veränderliche **Faltungsquadrupol**. Die Abstrahlung ist P = (32/5)(c⁵/G)·η_m²·(3Δs/r)⁵, also die fünfte Potenz des Taktunterschieds.
2. Im Starkfeld ist der Operator wirkungslos, und zwar in keiner endlichen Ordnung. Am Horizont von GW150914 ist q ≈ 10^(−1,96·10¹⁰). Wellenform und Ringdown sind die der ART. Das FTM bringt dort keine eigene Vorhersage mit und erbt die ART.
3. Der Ringdown ist die Eigenschwingung des Faltungstrichters zwischen r₀ = 6Δs und dem flachen Außenraum. Für GW150914 folgt f = c/(2π·3Δs)·F(j) = 272,4 Hz im Quellsystem bzw. **249,9 Hz** im Detektor. Die Dämpfungszeit beträgt **4,03 ms** im Detektor.
4. Das „Gedächtnis des Raumes" hat ein messbares Gegenstück, den nichtlinearen Memory-Effekt. In Faltungsgrößen gilt h_mem = 12·κ(Θ)·Δs_rad/D_L. Für GW150914 ergibt das 2,7·10⁻²² bei Kantenlage.
5. Abweichungen sind nur im Tief-a₀-Regime möglich. Gravitationswellen-Quellen liegen dort nie, weil die Schwellenfrequenz f*_gw ≲ 10⁻¹³ Hz beträgt. Prüfbar bleiben die Ausbreitung durch kollektive Trichter, Pulsar-Beschleunigungen und die Frage nach einer „Atmungs"-Polarisation.

---

## 0. Wellengleichung für lokale Störungen des Faltungsmediums

### 0.1 Was „Faltung" als Feld ist

Die Faltung Δs ist eine **Zusatzstrecke**. Sie ändert Abstände zwischen Punkten. Die natürliche Feldgröße ist deshalb die Änderung des Abstandsquadrats:

```
ds² = (η_μν + h_μν) dx^μ dx^ν ,     |h_μν| ≪ 1
```

h_μν ist ein symmetrischer Tensor mit 10 Komponenten. Das ist keine Zusatzannahme, sondern die Lesart „Faltung = Streckenänderung" in allen Richtungen, einschließlich der Zeitrichtung (Takt).

**Statischer Grenzfall des Modells.** Um eine Masse muss h den Takt des Autors liefern:

```
h_00 = 6Δs/r   ⇒   dτ/dt = √(1 − h_00) = √(1 − 6Δs/r)
```

Probe (Skript, Abschnitt 0):

| Ort | h_00 | 1−√(1−h_00) | 3Δs/r |
|---|---|---|---|
| Sonne, 1 AE | 1,974·10⁻⁸ | 9,8706·10⁻⁹ | 9,8706·10⁻⁹ |
| Sonnenoberfläche | 4,245·10⁻⁶ | 2,1225·10⁻⁶ | 2,1225·10⁻⁶ |
| Neutronenstern | 0,3445 | 0,1904 | 0,1723 |

### 0.2 Die Feldgleichung

Drei Regeln des Modells legen die lineare Gleichung fest:

- **R1, Hüll-Regel:** c_lokal = c. Jede lokale Umordnung der Hülle breitet sich mit dem Skalierungsfaktor c aus. Der Laplace-Operator des statischen Falls wird deshalb zum Wellenoperator □ = −∂²_t/c² + ∇².
- **R2, Faltungserhaltung:** Δs_ges = const. Die Quelle ist erhalten, ∂_μT^μν = 0. Die Feldgleichung muss daher mit einer Eichfreiheit h_μν → h_μν + ∂_με_ν + ∂_νε_μ verträglich sein.
- **R3, statischer Grenzfall:** h_00 = 6Δs/r.

Mit der spurumgekehrten Größe h̄_μν = h_μν − ½η_μν h und der Lorenz-Eichung ∂^μh̄_μν = 0 gibt es genau eine lineare Gleichung, die R1 bis R3 erfüllt:

```
□ h̄_μν = −48π · σ_μν ,        σ_μν = G·T_μν / (3c⁴)        (Faltungsdichte, Einheit m⁻²)
```

Der Faktor 48π = 3·16π ist der bekannte Faktor 16πG/c⁴, umgerechnet in Faltungseinheiten. Für ruhende Materie gilt σ_00 = Gρ/(3c²), also Δs pro Volumen.

**Probe statisch.** h̄_00 = 12Δs/r erfüllt ∇²h̄_00 = −48π·Δs·δ³(x). Der Gauß-Fluss ergibt −7,4223·10⁴ m, genau −48π·Δs☉. Daraus folgt:

```
h_00 = h̄_00/2 = 6Δs/r        h_ij = (6Δs/r)·δ_ij
```

Die Faltung wirkt also **gleich stark auf Takt und Strecke**. Das hat eine messbare Folge: Die Lichtablenkung am Sonnenrand beträgt α = 12Δs/b = **1,7512″**, der gemessene Wert. Der Takt (h_00) allein trägt davon die Hälfte bei, 6Δs/b = 0,8756″; die gemessene Ablenkung verlangt Takt und Strecke gemeinsam. Das Faltungsfeld ist ein Tensor.

**Status.** Die Gleichung ist die linearisierte Einstein-Gleichung. Das ist kein Zufall: Ein masseloses symmetrisches Tensorfeld, das an eine erhaltene Quelle koppelt, führt bei linearer Ordnung eindeutig auf diese Form (Fierz & Pauli 1939). Die nichtlineare Vervollständigung durch Selbstkopplung führt auf die Einstein-Gleichungen (Deser 1970). Das FTM liefert hier die **Lesart**, nicht eine neue Gleichung.

### 0.3 (a) Ausbreitung mit c, zwei transversale Polarisationen

Im Vakuum gilt □h̄_μν = 0, also ω = c·|k|. Die Störungen laufen mit c, ohne Dispersion.

Zählung der Freiheitsgrade:

```
10 (symmetrischer Tensor) − 4 (Lorenz-Eichung) − 4 (Resteichung) = 2
```

In der transversal-spurfreien Eichung (TT) mit Ausbreitung entlang z bleiben:

```
h_xx = −h_yy = h₊ ,   h_xy = h_× ,   k^i h_ij = 0 ,   h^i_i = 0
```

Beide Moden stehen **quer** zur Ausbreitung (tangential zur Wellenfront) und sind **spurfrei**. Ein Ring freier Teilchen wird elliptisch verformt, seine Fläche bleibt erhalten. Die relative Längenänderung beträgt δL/L = h/2.

Die Spur, also die isotrope „Atmungs"-Verformung, ist **kein** propagierender Freiheitsgrad. Sie ist über Zwangsgleichungen instantan an die Quelle gebunden, so wie h_00 = 6Δs/r. In Worten: **Die Faltungstiefe Δs selbst wandert nicht als Welle, nur die Umordnung ihrer Form.**

### 0.4 (b) Quelle: der Faltungsquadrupol

- **Monopol:** Nach R2 ist Δs_ges erhalten. Ein Monopol kann deshalb nicht strahlen. Ein kugelsymmetrisch pulsierender Stern ist stumm (Birkhoff).
- **Dipol:** Der Schwerpunkt der Faltung bewegt sich gleichförmig (Impulserhaltung). Auch der Dipol strahlt nicht.
- **Niedrigste Ordnung:** der Quadrupol.

Faltungsquadrupol (Einheit m³):

```
𝒬_ij = Σ_k Δs_k (x_i x_j − δ_ij r²/3)          (Δs_k = G m_k/(3c²))
```

Fernfeld und Leistung:

```
h_ij^TT = (6/(c² D)) · 𝒬̈_ij^TT(t − D/c)
P       = (9/5) · (c⁵/G) · ⟨(d³𝒬_ij/dt³ / c³)²⟩ ,        c⁵/G = 3,6283·10⁵² W
```

Das ist die Quadrupolformel von Einstein, nur mit GM = 3c²Δs umgeschrieben.

**Kreisförmiges Doppelsystem** mit Gesamtfaltung Δs, η_m = m₁m₂/M² und Abstand r:

```
P = (32/5) · (c⁵/G) · η_m² · (3Δs/r)⁵
```

3Δs/r ist in erster Ordnung das Faltungs-Differential δ. **Die abgestrahlte Leistung ist die fünfte Potenz des Taktunterschieds am Bahnradius.** Werte für GW150914 (η_m = 0,2471):

| r | 3Δs/r | P |
|---|---|---|
| 18Δs (ISCO) | 0,1667 | 1,82·10⁴⁸ W |
| 12Δs | 0,25 | 1,39·10⁴⁹ W |
| 9Δs | 0,333 | 5,84·10⁴⁹ W |

Die gemessene Spitzenleistung von 3,6·10⁴⁹ W liegt zwischen 12Δs und 9Δs, wie erwartet für die Verschmelzung.

**Chirp-Faltung.** Mit der Chirp-Masse M_c = (m₁m₂)^{3/5}/M^{1/5} ist die Chirp-Faltung

```
Δs_c = G·M_c/(3c²)
```

Die Frequenzentwicklung (führende Ordnung) lautet in Faltungsgrößen:

```
ḟ      = (96/5) · π^{8/3} · (3Δs_c/c)^{5/3} · f^{11/3}
Δs_c   = (c/3) · [ (5/96) · π^{−8/3} · f^{−11/3} · ḟ ]^{3/5}              (Faltungsmesser aus dem Chirp)
t_c(f) = (5/256) · (3Δs_c/c)^{−5/3} · (πf)^{−8/3}
h₀     = (4/D_L) · (3Δs_c)^{5/3} · (πf/c)^{2/3}                          (optimale Orientierung)
```

| | M_c (Quelle) | Δs_c Quelle | Δs_c Detektor ((1+z)M_c) |
|---|---|---|---|
| GW150914 (36 + 29 M☉, z = 0,09) | 28,10 M☉ | 13 829 m | 15 074 m |
| GW170817 (1,46 + 1,27 M☉, z = 0,0099) | 1,185 M☉ | 583,2 m | 589,0 m |

Nachrechnung, jeweils mit der Detektor-Chirp-Faltung:

| Ereignis | f [Hz] | ḟ [Hz/s] | Zeit bis Verschmelzung | Zyklen | h₀ (optimal, D_L) |
|---|---|---|---|---|---|
| GW150914, 410 Mpc | 35 | 79,7 | 0,165 s | 9 | 9,3·10⁻²² |
| GW150914 | 150 | 1,66·10⁴ | 3,4 ms | 1 | 2,45·10⁻²¹ |
| GW170817, 40 Mpc | 24 | 0,0899 | **100,1 s** | 3 844 | 3,3·10⁻²³ |
| GW170817 | 100 | 16,8 | 2,23 s | 356 | 8,6·10⁻²³ |

- Die Dauer von rund 0,2 s im Band (GW150914) und von rund 100 s ab 24 Hz (GW170817) sind die beobachteten Größenordnungen.
- Die Rückrechnung von Δs_c aus (f, ḟ) liefert den Eingabewert exakt, der Faltungsmesser ist also konsistent.

**Unabhängige Probe, Hulse-Taylor-Pulsar PSR B1913+16:**
- Δs_c = 605,78 m
- Ṗ_b = −(192π/5)·(2π·3Δs_c/(c·P_b))^{5/3}·f(e) = **−2,4021·10⁻¹²**
- ART-Literaturwert: −2,40263·10⁻¹² (Weisberg & Huang 2016)

**Faltungsbilanz GW150914.** Es gilt Δs(65 M☉) − Δs(62 M☉) = 1 476,6 m = Δs(3 M☉). Die Faltung, die dem Endloch fehlt, ist als Wellenenergie abgewandert.

### 0.5 (c) Ringdown als Eigenschwingung des Faltungstrichters

Nichtlinear hat der statische Trichter die Form g_tt = −(dτ/dt)², g_rr = (dτ/dt)⁻² (Schwarzschild). Störungen um diesen Trichter gehorchen der Regge-Wheeler-Gleichung. In Faltungsgrößen lautet sie:

```
d²ψ/dr_T² + [ω²/c² − V(r)] ψ = 0
V(r)  = (1 − 6Δs/r) · [ ℓ(ℓ+1)/r² − 18Δs/r³ ]                (axial, Spin 2)
r_T   = r + 6Δs · ln(r/(6Δs) − 1)        (Takt-Koordinate, dr_T = dr / (dτ/dt)²)
```

Die Randbedingungen lesen sich im Modell anschaulich:
- **Bei r → r₀ = 6Δs**, dem Grenzradius mit Takt = 0, ist r_T → −∞. Die Welle läuft nur hinein, sie „verschwindet in der Falte".
- **Bei r → ∞**, dem flachen Außenraum, läuft sie nur hinaus.

Der Trichter ist also ein offener Resonator zwischen r₀ und dem Außenraum. Sein Potentialwall liegt bei ℓ = 2 an r = 9,842Δs, also nahe der Photonensphäre 9Δs.

**Drei Rechnungen für ℓ = m = 2, n = 0, j = 0 (Einheiten GM/c³ = 3Δs/c):**

| Methode | M·ω |
|---|---|
| Leaver-Kettenbruch (exakt, im Skript gelöst) | **0,373672 − 0,088962 i** (Literatur: 0,373672 − 0,088962 i) |
| WKB 1. Ordnung (Schutz & Will 1985), eigene Auswertung | 0,3988 − 0,0883 i |
| Eikonal, Photonensphäre | (ℓ+½)/√27 − i·½/√27 = 0,4811 − 0,0962 i |

Das Eikonal-Bild ist die anschaulichste Lesart (Cardoso et al. 2009):
- Die Ringdown-Frequenz ist (ℓ+½)-mal die Umlauffrequenz des Lichts auf der Photonensphäre, ω_ph = c/(3√27·Δs).
- Die Dämpfung ist die Rate, mit der diese instabile Lichtbahn zerfällt.

Bei ℓ = 2 ist dieses Bild auf etwa 30 % genau, der Kettenbruch ist exakt.

**Die Verbindung zum Horizont r₀ = 6Δs** ist also eine Randbedingung und kein Ort der Schwingung. Das Frequenzmuster setzt der Wall bei ≈ 9–10Δs.

### 0.6 (d) Atmen a(t) gegen lokale Wellen

| | Atmen a(t) = L(t)/L(t₀) | lokale Welle h_ij^TT(t, x) |
|---|---|---|
| Art | globaler Modus des ganzen Mediums | lokale Störung auf dem Medium |
| Wellenzahl | k = 0 (homogen) | k ≠ 0 |
| Symmetrie | ℓ = 0, Spur (isotrop) | spurfrei, transversal, Spin 2 |
| Gleichung | Friedmann-Zwang, an die mittlere Faltungsdichte gebunden | Wellengleichung ḧ + 3Hḣ + (ck/a)²h = 0 |
| Zeitskala | 1/H₀ = 13,97 Gyr; H₀/2π = a₀/c = 3,61·10⁻¹⁹ Hz | f_LIGO/(H₀/2π) = 2,8·10²⁰, f_PTA/(H₀/2π) = 2,8·10¹⁰ |
| strahlt? | nein (Monopol, siehe 0.4) | ja |

Die Wellen spüren das Atmen **nur parametrisch**:
- Die Frequenz wird gestreckt: f_beob = f_emit · L(t_e)/L(t₀) = f_emit/(1+z).
- Die Amplitude fällt wie 1/(a·r). Daraus folgt d_L^GW = d_L^EM, also Ξ = 1.

Deshalb stehen in 0.4 Quell- und Detektorwerte getrennt: Detektor = (1+z) × Quelle, bei Massen bzw. Faltungen.

**Wichtige Grenze.** Wenn die Schleifen **ungleich** atmen (L₁, L₂, L₃ mit verschiedenen Raten), ist das eine homogene Scherung, also der Grenzfall k → 0 einer Tensorwelle. Sie ist vom CMB begrenzt (Saadeh et al. 2016):
- (σ/H)₀ < 10⁻⁶ für die schwächste, reguläre Tensormode, daraus Ω_σ < 3,3·10⁻¹³;
- (σ/H)₀ < 4,7·10⁻¹¹ für die Vektormode.

**Die drei Torus-Schleifen atmen also gleich, bis auf etwa 10⁻⁶.**

---

## 1. Der Operator im Starkfeld

Mit a = g = 3c²Δs/r² gilt √(a/a₀) = r*/r, also q = exp(−r*/r). Die Zusatzbeschleunigung relativ zu Newton beträgt ν − 1 = q/(1−q) ≈ q.

**GW150914** (Quellsystem: M = 65 M☉, M_f = 62 M☉, j = 0,67). Es gilt Δs(M_f) = 30 517 m, r₀ = 183,10 km und r*(M_f) = 55 352 AE = 0,268 pc.

| Ort | r | a_N [m/s²] | a_N/a₀ | r*/r | q = ν − 1 |
|---|---|---|---|---|---|
| Horizont Schwarzschild r₀ = 6Δs | 183,1 km | 2,45·10¹¹ | 2,0·10²¹ | 4,52·10¹⁰ | 10^(−1,96·10¹⁰) |
| Horizont Kerr r₊ = 1,7424 GM/c² | 159,5 km | 3,23·10¹¹ | 2,7·10²¹ | 5,19·10¹⁰ | 10^(−2,25·10¹⁰) |
| Photonensphäre 9Δs | 274,7 km | 1,09·10¹¹ | 9,1·10²⁰ | 3,01·10¹⁰ | 10^(−1,31·10¹⁰) |
| ISCO Kerr (j = 0,67): 3,528 GM/c² = 10,58Δs | 323,0 km | 7,89·10¹⁰ | 6,6·10²⁰ | 2,56·10¹⁰ | 10^(−1,11·10¹⁰) |
| ISCO Schwarzschild 18Δs(M) | 575,9 km | 2,60·10¹⁰ | 2,2·10²⁰ | 1,47·10¹⁰ | 10^(−6,39·10⁹) |
| Inspiral f_gw = 150 Hz | 338,7 km | 7,52·10¹⁰ | 6,3·10²⁰ | 2,50·10¹⁰ | 10^(−1,09·10¹⁰) |
| Inspiral f_gw = 35 Hz | 893,6 km | 1,08·10¹⁰ | 9,0·10¹⁹ | 9,49·10⁹ | 10^(−4,12·10⁹) |
| Inspiral f_gw = 20 Hz | 1 298 km | 5,12·10⁹ | 4,3·10¹⁹ | 6,53·10⁹ | 10^(−2,84·10⁹) |

Die ISCO-Frequenz beträgt f_gw(ISCO) = 67,6 Hz im Quellsystem.

**GW170817** (1,46 + 1,27 M☉):

| Ort | r | a_N/a₀ | q |
|---|---|---|---|
| NS-Oberfläche | 12 km | 1,1·10²² | 10^(−4,60·10¹⁰) |
| Kontakt | 24 km | 5,2·10²¹ | 10^(−3,14·10¹⁰) |
| Inspiral 24 Hz | 399 km | 1,9·10¹⁹ | 10^(−1,89·10⁹) |

Die ISCO-Frequenz beträgt 1 611 Hz. Zum Vergleich: Beim Hulse-Taylor-Pulsar ist q = 10^(−3,94·10⁵), am Detektor im Erdfeld ist q = 10^(−1,24·10⁵).

Die Tabelle verwendet die Newton-Beschleunigung GM/r². Die Eigenbeschleunigung eines ruhenden Beobachters ist im Starkfeld größer, q wäre damit noch kleiner. Die Werte oben sind also obere Schranken.

**In welcher Ordnung weicht das FTM ab?** Mit ε_a = a₀/a lautet der Operator q = exp(−ε_a^{−1/2}). Alle Ableitungen nach ε_a verschwinden bei ε_a → 0. Die Korrektur ist deshalb **in keiner endlichen post-Newtonschen Ordnung** vorhanden, sie liegt „jenseits aller Ordnungen". Vergleich im Skript:

| ε_a | q | ε_a¹⁰ |
|---|---|---|
| 10⁻² | 4,5·10⁻⁵ | 10⁻²⁰ |
| 10⁻⁶ | 10^(−434) | 10⁻⁶⁰ |

Bei GW150914 ist ε_a ≈ 10⁻²¹.

**Folgerung:**
- Inspiral, Verschmelzung und Ringdown sind im FTM **identisch mit der ART**. Die Abweichung ist ≲ 10^(−10⁹) in der Amplitude.
- Die Tests der ART mit LIGO/Virgo erlauben Abweichungen der PN-Koeffizienten von etwa 10 bis 100 %. Sie sind für das FTM daher keine Gefahr, aber auch keine Bestätigung.
- Die Starkfeld-Dynamik übernimmt das FTM als nichtlineare Vervollständigung (Abschnitt 0.2). Ihre Herleitung aus eigenen Regeln ist offen (Abschnitt 7).

**Der Operator wirkt parametrisch.** Die Welle selbst erzeugt an den LIGO-Spiegeln eine Gezeitenbeschleunigung von 1,78·10⁻¹² m/s² = 0,0148·a₀. Sie ist nicht das Argument von ν. Das Argument ist das Newton-Feld |∇Φ_N| des Trichters, in dem der Detektor liegt, wie in der Feldform ∇·[ν(|∇Φ_N|/a₀)∇Φ_N]. Am Detektor ist das die Erde mit 9,81 m/s², also q = 10^(−1,24·10⁵).

---

## 2. Ringdown in FTM-Größen

Die Zeitskala ist GM/c³ = 3Δs/c. Für die Grundmode ℓ = m = 2, n = 0 von Kerr verwendet die Notiz den Fit von Berti, Cardoso & Will (2006, PRD 73, 064030):

```
M·ω_R = F(j) = f₁ + f₂(1−j)^{f₃},        f₁ = 1,5251,  f₂ = −1,1568,  f₃ = 0,1292
Q_QNM(j) = b₁ + b₂(1−j)^{b₃},           b₁ = 0,7000,  b₂ = 1,4187,  b₃ = −0,4990
τ_QNM = Q_QNM/(π f)
```

In Faltungsgrößen:

```
f_QNM = c/(2π·3Δs) · F(j)          τ_QNM = (3Δs/c) · 2Q_QNM(j)/F(j)
```

Probe bei j = 0: F = 0,3683 gegen 0,3737 aus dem Kettenbruch, Q_QNM = 2,119 gegen 2,100. Der Fit ist dort auf etwa 1,5 % genau.

| j | F(j) | Q_QNM(j) | M·ω_I |
|---|---|---|---|
| 0 | 0,3683 | 2,119 | 0,0869 |
| 0,5 | 0,4674 | 2,705 | 0,0864 |
| **0,67** | **0,5227** | **3,167** | **0,0825** |
| 0,9 | 0,6660 | 5,176 | 0,0643 |
| 0,99 | 0,8870 | 14,82 | 0,0299 |

**GW150914** (M_f = 62 M☉, j = 0,67, z = 0,09):

| Größe | Wert |
|---|---|
| Δs(M_f) | 30 517 m, also 3Δs = 91,55 km |
| 3Δs/c | 0,30538 ms |
| c/(2π·3Δs) | 521,17 Hz |
| f_QNM (Quelle) | 521,17 Hz × 0,5227 = **272,4 Hz** |
| τ_QNM (Quelle) | **3,70 ms** |
| f_QNM (Detektor, ÷(1+z)) | **249,9 Hz** |
| τ_QNM (Detektor, ×(1+z)) | **4,03 ms** |
| Wellenlänge c/f | 1 101 km = 36,1Δs = 6,0·r₀ |
| τ_QNM in Takt-Einheiten 3Δs/c | 12,1 |

Die Masse ist mit ±4 M☉ unsicher. Im Detektor ergeben sich damit:

| M_f | f | τ_QNM |
|---|---|---|
| 58 M☉ | 267,1 Hz | 3,77 ms |
| 62 M☉ | 249,9 Hz | 4,03 ms |
| 66 M☉ | 234,8 Hz | 4,29 ms |

Das trifft die Referenzwerte ~250 Hz und ~4 ms. Ohne Spin (j = 0) läge die Frequenz im Quellsystem bei 194,7 Hz (521,17 Hz · 0,37367), im Detektor bei 178,7 Hz (mit dem Fit 176,1 Hz). [korr.] Die Messung verlangt also den Spin, also die Rotation des Trichters.

**Lesart:** Die Ringdown-Wellenlänge beträgt rund sechs Grenzradien. Das Endloch klingt in etwa 12 Takt-Einheiten 3Δs/c aus. Das FTM fügt hier nichts hinzu (Abschnitt 1). Es drückt die ART-Zahlen in der Faltungstiefe des Endlochs aus.

---

## 3. Ausbreitung: GW170817 / GRB 170817A

Der Gammablitz kam 1,74 s nach dem Wellensignal an, bei D ≥ 26 Mpc:

| Größe | Wert |
|---|---|
| c·Δt/D (obere Grenze für c_gw > c) | 6,5·10⁻¹⁶ |
| mit bis zu 10 s Eigenverzögerung des Blitzes | 3,1·10⁻¹⁵ |
| Literatur (Abbott et al. 2017, ApJL 848, L13) | −3·10⁻¹⁵ ≤ Δv/c ≤ +7·10⁻¹⁶ |

**Was das FTM dazu sagt, in drei Sätzen:**

- **P1, Ausbreitung:** Transversale Störungen der Hülle laufen mit c_lokal = c (0.3). Der Operator ν tritt in der Wellengleichung nicht auf. Er ist weder Masse- noch Brechungsterm, die Dispersion lautet ω = ck.
- **P2, Parametrik:** ν hängt nur vom Hintergrundfeld |∇Φ_N|/a₀ ab, nicht von Amplitude oder Frequenz der Welle (Abschnitt 1).
- **P3, Universalität:** Im Tief-a₀-Regime ändert ν die **Geometrie des Hintergrunds**, durch die Licht und Welle gemeinsam laufen. Seine Wirkung auf Wellen ist deshalb **transversal**: zusätzliche Ablenkung (Linse) und zusätzliche Shapiro-Verzögerung, beide für Licht und Gravitationswelle gleich. Auf die Ausbreitungsgeschwindigkeit wirkt ν nicht.

**Warum P3 notwendig ist.** Sähe nur das Licht den logarithmischen Trichter jenseits r*, aber nicht die Welle (oder umgekehrt), entstünde ein Laufzeitunterschied. Für den Weg von 8,2 auf 200 kpc durch die Milchstraße:

```
Δt = 2(v_flach/c)² · [r_e − R₀ − R₀ ln(r_e/R₀)] / c
```

| v_flach | Δt |
|---|---|
| 175,8 km/s | 136 Tage |
| 230 km/s | 232 Tage |

Beobachtet sind 1,74 s. Dieselbe Abschätzung für Theorien, in denen der Ersatz für Dunkle Materie nur an Materie und Licht koppelt, ergibt etwa 400 Tage (Boran, Desai, Kahya & Woodard 2018, PRD 97, 041501). **Das FTM besteht diesen Test, weil die Zusatzwirkung Geometrie ist (Beobachter-Regel), kein Feld, das an einen Teil der Welt koppelt.** Das ist eine Festlegung des Modells.

**Kann das Atmen eine Dispersion erzeugen?** Selbst wenn das Atmen der Welle eine Masse der Größe ħ·H₀/2π gäbe (= 2,4·10⁻³⁴ eV), wäre das weit unter den Grenzen (1,2·10⁻²² eV aus GW150914, 1,76·10⁻²³ eV aus GWTC-3). Die Laufzeitabweichung ½(ω_a/ω)² beträgt:

| Band | Frequenz | Abweichung |
|---|---|---|
| LIGO | 100 Hz | 1,7·10⁻⁴³ |
| LISA | 1 mHz | 1,7·10⁻³³ |
| PTA | 10 nHz | 1,7·10⁻²³ |

Das Atmen kann die Ausbreitungsgeschwindigkeit also nicht messbar verändern.

---

## 4. Gedächtnis des Raumes = Gravitationswellen-Memory

Die ART sagt eine **bleibende** Verschiebung freier Testmassen nach dem Durchgang einer Welle voraus:
- **linear:** Zel'dovich & Polnarev 1974; Braginsky & Thorne 1987;
- **nichtlinear**, also von der Wellenenergie selbst erzeugt: Christodoulou 1991 (PRL 67, 1486); Wiseman & Will 1991; Thorne 1992; Blanchet & Damour 1992.

Nach dem Signal kehrt der Abstand nicht zurück. Das ist die **Knickspur** des FTM: „Gefalteter Raum vergisst nicht", und zwar als messbare Geometrie ohne Speicher-Objekt.

**Formel in Faltungsgrößen.** Das führende Memory eines kreisförmigen Doppelsystems ist h_mem = (2η_mMx/D)·sin²Θ(17+cos²Θ)/96 (Wiseman & Will 1991; Favata 2009). Die abgestrahlte Energie ist E_rad = η_mMx/2 = 3c⁴·Δs_rad/G. Daraus folgt:

```
h_mem = 12 · κ(Θ) · Δs_rad / D_L ,       κ(Θ) = sin²Θ (17 + cos²Θ)/96
κ(90°) = 0,1771  (Kante: h_mem = 2,125·Δs_rad/D_L),   κ(150°) = 0,0462,   κ(0°) = 0
```

**Die bleibende Knickspur ist die Faltungstiefe der abgestrahlten Energie, geteilt durch die Entfernung.** Δs_rad ist genau die Faltung, die dem Endloch fehlt (0.4).

**GW150914** (E_rad = 3,0 M☉c², D_L = 410 Mpc, (1+z)·Δs_rad = 1 609,5 m, Spitzenamplitude 1,0·10⁻²¹):

| Orientierung | h_mem | Anteil an h_peak = 1,0·10⁻²¹ |
|---|---|---|
| Kante (90°) | **2,7·10⁻²²** | 27 % |
| nahezu frontal (Beispiel 150°) | 7,1·10⁻²³ | 7 % |

- Bleibende Armverschiebung von LIGO bei Kantenlage: ΔL = ½h·L = 5,4·10⁻¹⁹ m.
- Einschränkung: Das ist eine Abschätzung in führender Ordnung, kalibriert über E_rad. Mit numerischer Relativität kalibrierte Rechnungen (Favata 2009) können davon um einen Faktor der Ordnung 1 abweichen. Die Größenordnung 10⁻²², also 5–25 % der Spitzenamplitude, ist robust.

**Nachweisaussichten:**

| Detektor | Aussicht |
|---|---|
| LIGO, Einzelereignis | Das Memory ist eine langsame Stufe, LIGO sieht nur ihren Anstieg. Für GW150914-artige Ereignisse ist das SNR ≈ 0,5 pro Ereignis. **Stacking:** SNR 3 mit ~35, SNR 5 mit ~90 Ereignissen bei Design-Empfindlichkeit (Lasky et al. 2016, PRL 117, 061102). Erste Katalogsuchen: Hübner et al. 2020 (PRD 101, 023011). Ausblick: Grant & Nichols 2023 (PRD 107, 064056). |
| GW170817 | h_mem ≈ 4,8·10⁻²³ (Kante, nur Inspiral), nicht nachweisbar |
| **LISA** | Schwarzes-Loch-Paar mit 10⁶ M☉, z = 1 (D_L = 6,79 Gpc): h_mem ≈ 4,8·10⁻¹⁹ bei einer Anstiegszeit ~10³ s, also mitten im LISA-Band. Nachweisbar bis z ≈ 2 (Favata 2009). Bei 10⁷ M☉ und z = 3 ist h_mem ≈ 2,5·10⁻¹⁸. |
| **PTA** | Paar mit 10¹⁰ M☉ bei z = 0,2: h_mem ≈ 1,9·10⁻¹⁴. Das ergibt einen linear wachsenden Timing-Rest von rund 4,6 µs über 15 Jahre (10⁹ M☉ bei z = 0,05: 1,8 µs). Solche Ereignisse sind nachweisbar, aber selten. Suchen: Arzoumanian et al. 2015 (ApJ 810, 150), Aggarwal et al. 2020 (ApJ 889, 38). Theorie: van Haasteren & Levin 2010, Seto 2009, Pshirkov, Baskaran & Postnov 2010. |

**Gedächtnis im FTM: eins, keine zweite Art auf gebundenen Bahnen.** [korr.] Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m = ω_orb/N*, N* ≈ 0,3–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Geschwindigkeitskern e^{−λ_mΔt} würde exzentrische Bahnen runden, alte Scheiben verlangen dafür λ_m ≪ 0,01·Ω. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche). Im Kleinen zeigt sich dieselbe Knickspur als Memory-Effekt:

| | Knickspur im Kosmos: Zustand der Schließungsregel [korr.] | Knickspur im Kleinen (Raumgedächtnis) |
|---|---|---|
| wirkt auf | den Zustand der Regel in gebundenen Bereichen, nicht auf die einzelne Bahn [korr.] | Abstände freier Massen |
| Regime | Übergabe von ξ in die Regel, wo die Mode nicht ganz übergegangen ist [korr.] | alle Regime, auch Starkfeld |
| Dauer | Gewicht w bleibend, Übergabe in 1/λ_m = N*/ω_orb [korr.] | bleibend [korr.] |
| Beleg | Haufenstöße (N* aus dem Bullet Cluster) [korr.] | Memory-Effekt der ART |

Der Memory-Effekt bestätigt die FTM-Aussage vom bleibenden Gedächtnis. Er unterscheidet das FTM aber nicht von der ART, weil beide ihn gleich vorhersagen.

---

## 5. Wo das FTM abweichen könnte: das Tief-a₀-Regime

### 5.1 Warum Gravitationswellen-Quellen nie im Tief-a₀-Regime liegen

Ein Doppelsystem erreicht die Bahnbeschleunigung a₀ erst beim Abstand r*. Die Gravitationswellen-Frequenz dort ist die **Schwellenfrequenz**:

```
f*_gw = (1/π)·√(a₀/r*) = a₀^{3/4} / (π·(3c²Δs)^{1/4})
```

| M [M☉] | r* | f*_gw | Umlaufzeit |
|---|---|---|---|
| 1 | 0,034 pc | 1,1·10⁻¹³ Hz | 5,9·10⁵ a |
| 60 | 0,264 pc | 3,9·10⁻¹⁴ Hz | 1,6·10⁶ a |
| 10⁶ | 34 pc | 3,4·10⁻¹⁵ Hz | 1,9·10⁷ a |
| 10⁹ | 1,08 kpc | 6,0·10⁻¹⁶ Hz | 1,0·10⁸ a |
| 10¹⁰ | 3,41 kpc | 3,4·10⁻¹⁶ Hz | 1,9·10⁸ a |

Alle Detektoren liegen weit oberhalb von f*: PTA (10⁻⁹–10⁻⁷ Hz), LISA und LIGO. Selbst die Quellen im PTA-Band liegen tief im Newton-Regime. Ein Paar mit 10⁹ M☉ bei f_gw = 10 nHz hat den Abstand 3 425 AE, a = 4,2·10⁹·a₀ und q = 10^(−2,8·10⁴).

**Die Quellen sind immer Einzeltrichter.** Das FTM kann also nur **auf dem Weg** (5.2, 5.3) oder **im Hintergrund** (5.4) abweichen, nicht an der Quelle. Nur die CMB-B-Moden reichen in den Bereich um f* (etwa 10⁻¹⁸ bis 10⁻¹⁶ Hz), sie stammen aber nicht von gebundenen Paaren.

### 5.2 Pulsar-Timing: das galaktische Feld am Ort der Pulsare

Die Sonne liegt (R₀ = 8,2 kpc [G]) im Feld der baryonischen Milchstraße bei a_N = 1,244·10⁻¹⁰ m/s² = 1,036·a₀. Daraus folgen:

| Größe | Wert |
|---|---|
| q | 0,361 |
| ν | 1,566 |
| a = ν·a_N | 1,947·10⁻¹⁰ m/s² |
| gemessen, v_c = 233 km/s [G] | 2,146·10⁻¹⁰ m/s² (das grobe Punktmassen-Modell liegt 9 % darunter) |

**Prüfbare Größe:** Pulsare messen ihre Beschleunigung in der Galaxis direkt über Ṗ/P = a_los/c. Newton gibt dafür 4,15·10⁻¹⁹ s⁻¹, das FTM 6,50·10⁻¹⁹ s⁻¹ (Messprinzip: Chakrabarti et al. 2021, ApJL 907, L26).

[korr.] Die folgende Rechnung setzt die Feldform ∇²Φ = ∇·[ν∇Φ_N] als Bahnregel voraus; das trägt nicht: Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N; die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Im Modell gibt die örtliche Regel am Sonnenort die gelesene Flächendichte 71 M☉/pc².

Der Unterschied zu einem Halo aus Dunkler Materie liegt nicht im Betrag in der Ebene (beide sind an die Rotationskurve angepasst), sondern in der **Vertikalstruktur**. In der Feldform ergibt die Umrechnung von ν durch einen Beobachter, der ungeänderte Regeln annimmt, eine scheibenförmige gelesene Dichte (kein Stoff): [korr.]

```
ρ_X = ∇·[(ν−1)∇Φ_N] / (4πG)          (Lesung)
```

Ein Halo ist dagegen fast kugelförmig. In der Feldform gilt in führender Ordnung nahe der Scheibe Σ_eff ≈ ν·Σ_b; für eine angenommene baryonische Scheibe von Σ_b = 50 M☉/pc² ergäbe das 78 M☉/pc², mit der örtlichen Regel sind es am Sonnenort 71 M☉/pc². [korr.] **Test:** Das Profil a_z(z) der Pulsar-Beschleunigungen muss die gelesene Dichte als Scheibe zeigen, festgelegt durch die Baryonen; der Faktor ν ≈ 1,57 gilt nur in der Feldform. [korr.]

Auf die Gravitationswellen selbst, also Hellings-Downs-Korrelation und Ausbreitung, wirkt das nach P1–P3 nicht.

### 5.3 Wellen durch kollektive Trichter (Galaxien, Haufen)

Jenseits von r* ist der Trichter logarithmisch. Licht und Wellen erfahren dieselbe konstante Ablenkung:

```
α_∞ = 2π·(v_flach/c)² = 2π·√(GM·a₀)/c² = 6π·Δs/r*
```

| Linse | v_flach | α_∞ |
|---|---|---|
| Milchstraße 6·10¹⁰ M☉ | 175,8 km/s | 0,446″ |
| Galaxie 2·10¹¹ M☉ | 237,6 km/s | 0,814″ |
| Haufen 1,7·10¹⁴ M☉ | 1 283 km/s | 23,7″ |

Zwei prüfbare Aussagen:

1. **Gleiche Bilder, gleiche Verzögerungen.** Stark gelinste Gravitationswellen (in O4/O5 erwartet) müssen dieselben Bild-Zeitverzögerungen und Vergrößerungsverhältnisse zeigen wie das Licht der Wirtsgalaxie. Die Linsenmasse ist durch die **Baryonen** über ν festgelegt: M_eff(<b) = b·v_flach²/G. Ein Halo aus Dunkler Materie hat diese Bindung nicht, das FTM hat sie. Das ist ein Linsen-RAR-Test mit zwei Botenarten.
2. **Wellenoptik.** Unterhalb von f ≈ c³/(8πG·M_eff) beugt sich die Welle um die Linse, und die Verstärkung sinkt (Takahashi & Nakamura 2003). Mit M_eff(<100 kpc) ergibt sich:
   - Milchstraße: 7,2·10¹¹ M☉, f ≈ 1,1·10⁻⁸ Hz;
   - Haufen: 3,8·10¹³ M☉, f ≈ 2·10⁻¹⁰ Hz.

   Das liegt im PTA-Band. Die frequenzabhängige Verstärkung prüft dort das Profil der gelesenen Masse (∝ b) statt nur ihres Integrals. [korr.]

**Einschränkung:** Eine Übergangsfunktion vom MOND-Typ erklärt die Masse von Galaxienhaufen nur etwa zur Hälfte (Sanders 2003; Pointecouteau & Silk 2005; Bullet Cluster: Clowe et al. 2006). Im Modell schließt diesen Rest der stoßfreie Rest der Mode ξ der Schließungsregel: Die Operatorform fasst 21–71 % der Mode, der Rest trifft die Linsenmasse auf 1–5 %. [korr.]

### 5.4 Ein Beitrag des Atmens zu einem Hintergrund

- **Isotropes Atmen strahlt nicht** (0.4, 0.6). Es liefert keinen Wellenhintergrund, nur die Rotverschiebung.
- **Ungleiches Atmen** der Schleifen ist Scherung. Sie ist begrenzt auf (σ/H)₀ < 10⁻⁶, also Ω_σ < 3,3·10⁻¹³ (Saadeh et al. 2016).
- **Beim Atmen-Takt** f_b = H₀/2π = a₀/c = 3,6·10⁻¹⁹ Hz gilt allgemein Ω_GW(f_b) = (2π²/3H₀²)·f_b²·h_c² = **h_c²/6**. Eine Schwankung der Schleifenlängen mit der relativen Amplitude h_c trägt also Ω ≈ h_c²/6 bei.
- **Polarisationstest (PTA).** Hätten lokale Wellen eine Atmungs-Polarisation (Spur-Mode), zeigten Pulsarpaare das Korrelationsmuster (3 + cos ζ)/8 (Lee, Jenet & Price 2008) statt Hellings-Downs:

| ζ | 0° | 30° | 60° | 90° | 120° | 180° |
|---|---|---|---|---|---|---|
| Hellings-Downs (Tensor) | 0,500 | 0,212 | −0,082 | −0,145 | −0,011 | 0,250 |
| Atmen-Mode (3 + cos ζ)/8 | 0,500 | 0,483 | 0,438 | 0,375 | 0,313 | 0,250 |

**FTM-Vorhersage:** Das Atmen ist global (k = 0) und nicht lokal. Die Korrelation ist **reines Hellings-Downs**, ohne Atmungs-Anteil. Die NANOGrav-15-Jahres-Daten zeigen Hinweise auf Hellings-Downs (Agazie et al. 2023, ApJL 951, L8). Ein signifikanter Anteil der Atmungs-Mode stünde im Gegensatz zu dieser Vorhersage.

- **Standard-Sirenen:** Ξ = d_L^GW/d_L^EM = 1 (0.6). Theorien mit modifizierter Reibung sagen Ξ ≠ 1 voraus (Belgacem et al. 2018). LISA und das Einstein-Teleskop messen Ξ auf Prozent-Niveau.

---

## 6. Schatten eines Schwarzen Lochs (EHT)

Ein Photon hat den Stoßparameter b(r) = r/(dτ/dt) = r/√(1−6Δs/r). Das Minimum liegt bei

```
r_ph = 9Δs  (= 3GM/c² = 1,5·r₀),   dτ/dt(r_ph) = √(1/3) = 0,5774,   b_c = 9Δs/√(1/3) = 3√27·Δs = 15,588·Δs
```

Die kritische Bahn hängt in einer statischen Kugelgeometrie **nur vom Takt g_tt ab**, nicht vom Raumanteil. Der Schattenrand ist also rein durch dτ/dt bestimmt: **Der Schatten ist das Bild des Taktes.**

| | M87* | Sgr A* (4,0·10⁶, 8,2 kpc) | Sgr A* (GRAVITY: 4,297·10⁶, 8,277 kpc) |
|---|---|---|---|
| Δs | 3,20·10¹² m = 21,4 AE | 1,97·10⁹ m = 0,013 AE | 2,12·10⁹ m |
| r₀ = 6Δs | 128,3 AE | 0,079 AE | 0,085 AE |
| Photonensphäre 9Δs | 192,5 AE | 0,118 AE | 0,127 AE |
| θ_g = 3Δs/D | 3,819 µas | 4,815 µas | 5,124 µas |
| Schatten-Ø 2·3√27·Δs/D | **39,7 µas** | **50,0 µas** | 53,3 µas |
| Ring-Ø ≈ 11·θ_g (EHT-Kalibrierung) | 42,0 µas | 53,0 µas | 56,4 µas |
| gemessen | 42 ± 3 µas (EHT 2019, ApJL 875, L1) | 51,8 ± 2,3 µas (EHT 2022, ApJL 930, L12) | |
| q an 9Δs | 10^(−1,28·10⁶) | 10^(−5,2·10⁷) | |

- Der leuchtende Ring liegt etwas außerhalb der kritischen Kurve. Deshalb wird er über α ≈ 11·θ_g kalibriert.
- Der Spin verändert den Durchmesser nur um wenige Prozent (Johannsen & Psaltis 2010).
- Der Operator ist am Schattenrand wirkungslos. Das FTM sagt also den ART-Schatten voraus, **ausgedrückt durch Δs und den Takt**.

---

## 7. Offene Punkte – nächste Schritte

1. **Nichtlineare Vervollständigung aus eigenen Regeln.** Abschnitt 0 leitet die lineare Gleichung aus R1 bis R3 her. Im Modell ist der statische Trichter nicht übernommen: dτ/dt = √(1 − 6Δs/r) folgt im Vakuum aus G^t_t = G^r_r (g_tt·g_rr = −1) und ist damit hergeleitet, keine Setzung. Kerr ist die Wiedergabe der ART; offen bleiben dynamische Starkfeldsysteme (Verschmelzung). [korr.]
2. **Den Operator ν in die Metrik schreiben.** P3 verlangt, dass ν Geometrie ist, die Licht und Wellen gleich sehen. Das muss konkret werden: ein Φ = Φ_N + Φ_ν in h_00 **und** h_ij, mit demselben Faktor (γ = 1 auch im Tief-a₀-Regime). Dann gilt für Kugeln und außerhalb von Scheiben M_dyn = M_lens (η = 1); in Scheiben innen (R ≲ 3,5 R_d) ist η > 1, weil die Metrik nur den wirbelfreien Teil der Regel trägt. [korr.] Prüfen an der schwachen Linsenwirkung von Galaxien (Linsen-RAR).
3. **Bahn- gegen Raumgedächtnis.** [korr.] Geklärt: Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m = ω_orb/N*) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Kern e^{−λ_mΔt} auf gebundenen Bahnen ist keine zweite Gedächtnisart, λ_m keine Mischungsrate von Knickspuren. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche).
4. **Galaxienhaufen.** Eine Übergangsfunktion ν allein reicht in Haufen nicht, es fehlt etwa ein Faktor 2. Kandidaten innerhalb des Modells:
   - Summenwirkung kollektiver Trichter statt eines einzelnen ν(a_N);
   - höheres a₀ wegen a₀(z) = cH(z)/2π bei der Entstehung. [korr.] Das trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist a₀ seit z ≈ 10 fest; ≈ cH₀/2π gilt nur heute.

   [korr.] Im Modell schließt den Faktor der stoßfreie Rest der Mode ξ der Schließungsregel: Die Operatorform fasst 21–71 % der Mode, der Rest trifft die Linsenmasse auf 1–5 %.
5. **Linsen mit zwei Botenarten.** Für erwartete gelinste Gravitationswellen-Ereignisse die FTM-Vorhersage vorab festlegen: Zeitverzögerung und Vergrößerung aus den Baryonen über ν. Im PTA-Band zusätzlich die frequenzabhängige Verstärkung für Haufenlinsen (f ≲ 10⁻⁹ Hz).
6. **Pulsar-Beschleunigungen.** Das vertikale Profil a_z(z) aus Ṗ_b von Binärpulsaren mit der örtlichen Regel a = ν(|g_N|)·g_N (nicht der Feldform ∇·[ν∇Φ_N], die an SPARC außerhalb 95 % liegt) und einem realistischen Baryonenmodell der Scheibe rechnen, statt mit der Punktmasse. [korr.] Das ist ein unabhängiger Test der Regel am Ort der Sonne.
7. **Atmungs-Polarisation.** Die FTM-Vorhersage „reines Hellings-Downs, keine (3 + cos ζ)/8-Komponente" in die PTA-Analysen einbringen und an den kommenden IPTA-Daten prüfen.
8. **Memory als Bestätigung des Raumgedächtnisses.** Stacking mit LIGO/Virgo/KAGRA (etwa 35 bzw. 90 Ereignisse) und LISA-Einzelereignisse prüfen die Aussage „gefalteter Raum vergisst nicht" quantitativ. Mit der FTM-Formel h_mem = 12κΔs_rad/D_L eine eigene Vorhersageliste für die Katalogereignisse erstellen.

---

## Literatur

- Abbott et al. (LIGO/Virgo) 2016, PRL 116, 061102 (GW150914); PRL 116, 221101 (ART-Tests, Gravitonmasse).
- Abbott et al. 2017, PRL 119, 161101 (GW170817); ApJL 848, L13 (GW170817/GRB 170817A).
- Agazie et al. (NANOGrav) 2023, ApJL 951, L8.
- Aggarwal et al. 2020, ApJ 889, 38; Arzoumanian et al. 2015, ApJ 810, 150.
- Bardeen, Press & Teukolsky 1972, ApJ 178, 347 (ISCO).
- Belgacem, Dirian, Foffa & Maggiore 2018, PRD 98, 023510 (modifizierte GW-Ausbreitung, Ξ).
- Berti, Cardoso & Will 2006, PRD 73, 064030.
- Blanchet & Damour 1992, PRD 46, 4304; Christodoulou 1991, PRL 67, 1486; Thorne 1992, PRD 45, 520; Wiseman & Will 1991, PRD 44, R2945.
- Boran, Desai, Kahya & Woodard 2018, PRD 97, 041501.
- Braginsky & Thorne 1987, Nature 327, 123; Zel'dovich & Polnarev 1974, Sov. Astron. 18, 17.
- Cardoso, Miranda, Berti, Witek & Zanchin 2009, PRD 79, 064016.
- Chakrabarti et al. 2021, ApJL 907, L26.
- Clowe et al. 2006, ApJL 648, L109; Pointecouteau & Silk 2005, MNRAS 364, 654; Sanders 2003, MNRAS 342, 901.
- Deser 1970, Gen. Rel. Grav. 1, 9; Fierz & Pauli 1939, Proc. R. Soc. A 173, 211.
- EHT Collaboration 2019, ApJL 875, L1; 2022, ApJL 930, L12.
- Favata 2009, ApJL 696, L159.
- Grant & Nichols 2023, PRD 107, 064056; Hübner, Talbot, Lasky & Thrane 2020, PRD 101, 023011.
- Hellings & Downs 1983, ApJL 265, L39; Lee, Jenet & Price 2008, ApJ 685, 1304.
- Johannsen & Psaltis 2010, ApJ 718, 446.
- Lasky, Thrane, Levin, Blackman & Chen 2016, PRL 117, 061102.
- Leaver 1985, Proc. R. Soc. A 402, 285; Regge & Wheeler 1957, Phys. Rev. 108, 1063; Schutz & Will 1985, ApJL 291, L33.
- Peters & Mathews 1963, Phys. Rev. 131, 435; Weisberg & Huang 2016, ApJ 829, 55.
- Pshirkov, Baskaran & Postnov 2010, MNRAS 402, 417; Seto 2009, MNRAS 400, L38; van Haasteren & Levin 2010, MNRAS 401, 2372.
- Saadeh, Feeney, Pontzen, Peiris & McEwen 2016, PRL 117, 131302.
- Takahashi & Nakamura 2003, ApJ 595, 1039.
