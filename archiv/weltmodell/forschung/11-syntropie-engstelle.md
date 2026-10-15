# Syntropie und die Engstelle: lokal und global

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: 4.4 Folgerung „prüfbare Abweichung von ΛCDM, verträglich mit a₀ = cH/2π ∝ H, DESI zeigt Hinweise“ → Mindestrate neuer innerer Faltung Λ̇/Λ ≥ Ṡ_i/(2S_max), heute w_s = −1 + 1·10⁻¹⁹ (Endzustand 2·10⁻¹⁹), unmessbar, DESI erklärt sie nicht; a₀ ∝ H trägt nicht (N_eff); Abschnitt 6 Punkt 3: Verlust der Strahlung durch Rotverschiebung als Rückstellung σ → senkt Λ (heute w_s = −1,0002), nicht die Rückstellung; Abschnitt 6 Punkt 4: Mindestabnahme von H „aus a₀ = cH/2π und (S3)“ → aus (S3) allein, als Mindestrate neuer innerer Faltung, unmessbar. Umbenennungen (ohne Marke): χ (Zentrierung des FGMS) → c_z; im Gain G^Ω α → A_Ω, τ_D → 𝒟_c; Ω_Λ → Ω_DE, H_Λ → H_DE, Λ (kosmologische Konstante) → Λ_kosm. [G] ergänzt: Entropiedichte der Hintergrundstrahlung, Leistung und Temperaturen im Beispiel Erde, Werte von Egan und Lineweaver, DESI. Rechnet mit H₀ = 70 km/s/Mpc statt 67,4; a₀ kommt nicht als Zahl vor. Entropieproduktion σ_Σ, σ_ges → Ṡ_{i,Σ}, Ṡ_i (Notation).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: rechnungen/syntropie.py (alle Zahlen dieser Datei stammen aus ihrer Ausgabe).
> Eingaben: H₀ = 70 km/s/Mpc, Ω_m = 0,3, Ω_DE = 0,7, T_CMB = 2,7255 K, M_U = 1,5·10⁵³ kg, R_obs = 4,4·10²⁶ m.

## Fragestellung

Das Modell sagt: Außen wächst die Entropie. An der Engstelle des Torus, dem zentralen Rückfluss-Wirbel, verdichtet sich das Medium und ordnet sich neu. Dabei geht keine Information verloren.

Offen ist, was die „Engstelle" physikalisch ist. Zwei Lesarten stehen nebeneinander:

- **(L) lokal:** der Innenraum jedes dichten Trichters, also r ≤ r₀ = 6Δs.
- **(G) global:** das Loch des Torus, also die topologische Schließungsregion.

**Ergebnis vorweg:** Beide Lesarten gelten. Sie werden durch dieselbe Formel beschrieben, S = 36πk·(Δs/ℓ_P)², und unterscheiden sich nur darin, welches Δs man einsetzt. Thermodynamisch gilt dabei: Die Verdichtung zur Falte **erhöht** die Entropie. Syntropie ist Ordnungsbildung bei wachsender Gesamtentropie und nur mit Entropieexport möglich. Der Kreislauf schließt sich in der Form, nicht im Zustand. Er ist eine Spirale.

---

## 1 — Die lokale Engstelle in FTM-Größen

### 1.1 Formeln

Mit GM = 3c²Δs, r₀ = 6Δs = 2GM/c² und der Planck-Länge ℓ_P = √(ħG/c³) = 1,6162·10⁻³⁵ m:

```
Entropie        S   = k·A/(4ℓ_P²),  A = 4π(6Δs)²   ⇒   S = 36π·k·(Δs/ℓ_P)²
Temperatur      T_H = ħc³/(8πGMk)                  ⇒   T_H = ħc/(24π·k·Δs)
Verdampfung     t_ev = 5120π·G²M³/(ħc⁴)            ⇒   t_ev = 138240π·(Δs/ℓ_P)²·Δs/c = 3840·(S/k)·(Δs/c)
Smarr           M c² = 2·T_H·S
```

Die dritte Zeile zeigt: **Die Lebensdauer einer Falte ist 3840-mal ihre Entropie (in k), gemessen in der Lichtlaufzeit über Δs.** Das Skript prüft alle Umformungen gegen die Standardformen. Die Quotienten sind 1,000000000000.

Zu t_ev: Das ist Hawkings Abschätzung für masselose Strahlung ohne Grauekörperfaktoren. Page (1976) hat mit den realen Teilchenarten gerechnet. Das ändert den Vorfaktor um einen Faktor der Ordnung eins, die Skalierung mit M³ bleibt.

### 1.2 Zahlen

| Falte | Δs | r₀ = 6Δs | S/k | T_H [K] | t_ev [a] |
|---|---|---|---|---|---|
| stellar, 10 M☉ | 4,92 km | 29,5 km | 1,05·10⁷⁹ | 6,17·10⁻⁹ | 2,10·10⁷⁰ |
| 1 M☉ (hypothetisch) | 492 m | 2,95 km | 1,05·10⁷⁷ | 6,17·10⁻⁸ | 2,10·10⁶⁷ |
| Sgr A* (4,15·10⁶ M☉) | 0,0137 AE | 0,0819 AE | 1,81·10⁹⁰ | 1,49·10⁻¹⁴ | 1,50·10⁸⁷ |
| M87* (6,5·10⁹ M☉) | 21,4 AE | 128 AE | 4,43·10⁹⁶ | 9,49·10⁻¹⁸ | 5,76·10⁹⁶ |
| Milchstraße (6·10¹⁰ M☉) | 197 AE | 1 180 AE | 3,78·10⁹⁸ | 1,03·10⁻¹⁸ | 4,53·10⁹⁹ |
| Universum (1,5·10⁵³ kg) | 1 200 Mpc | 7 220 Mpc | 5,97·10¹²² | 8,18·10⁻³¹ | (9·10¹³⁵, physikalisch nicht erreichbar, siehe 3.3) |

### 1.3 Vergleich mit der Materie-Entropie

| System | S/k |
|---|---|
| Sonne als Gaskugel (Sackur-Tetrode, Virialtemperatur 6,95·10⁶ K, μ = 0,605) | 2,56·10⁵⁸ (21,4 k pro Baryon) |
| dieselbe Masse als Falte | 1,05·10⁷⁷ |
| Verhältnis | **4,1·10¹⁸** |
| CMB-Photonen im beobachtbaren Universum (1 479 k/cm³ [G] · 3,57·10⁸⁰ m³) | 5,28·10⁸⁹ |
| Reliktneutrinos | 5,04·10⁸⁹ |
| eine einzige Falte von 10⁶ M☉ | 1,05·10⁸⁹ |

Eine Falte von 2,2·10⁶ M☉, also etwa Sgr A* in halber Größe, trägt so viel Entropie wie die gesamte Hintergrundstrahlung des beobachtbaren Universums.

### 1.4 Folgerung: Verdichtung ist Entropiezunahme (Penrose)

Für nicht gravitierende Freiheitsgrade ist Gleichverteilung der Zustand höchster Entropie. Für gravitierende gilt das Umgekehrte: Klumpung erhöht die Entropie, und die Falte ist ihr Maximum (Penrose 1979). Der Kollaps der Sonne zur Falte würde die Entropie um 18,6 Größenordnungen erhöhen.

Damit ist der Satz „an der Engstelle … ordnet sich das Medium neu" thermodynamisch so zu lesen:

- **Mikroskopisch** ist die Engstelle der Ort **maximaler** Entropie pro Masse.
- **Makroskopisch** ist sie der Ort **maximaler Einfachheit**. Nach dem Keine-Haare-Satz wird eine stationäre Falte durch drei Zahlen vollständig beschrieben: M, J und die Ladung Q_el. Das ist die geometrisch einfachste Form, die das Medium annehmen kann.

Das „ordnet sich" trifft die makroskopische Seite: Form wird auf Geometrie reduziert. Das entspricht der Beobachter-Regel, denn an der Engstelle bleibt kein Objekt übrig, nur Geometrie. Thermodynamische Ordnung im Sinne niedriger Entropie ist sie nicht.

### 1.5 Syntropie sauber definiert: dissipative Strukturen

Für ein Teilsystem Σ (Stern, Planet, Galaxie, Lebewesen) gilt die Bilanz nach Prigogine:

```
dS_Σ/dt = Ṡ_{i,Σ} − J,          Ṡ_{i,Σ} ≥ 0   (innere Produktion),   J = Entropieexport über den Rand
```

**Definition (Syntropie).** Syn_Σ := S_Σ^eq − S_Σ ≥ 0. Das ist der Abstand zum Gleichgewichtszustand des Teilsystems bei gleicher Energie und gleichem Volumen, also die Negentropie im Sinne Schrödingers.

```
dSyn_Σ/dt = J − Ṡ_{i,Σ}        Syntropie wächst  ⇔  J > Ṡ_{i,Σ}
dS_ges/dt = Ṡ_{i,Σ} + dS_Umgebung/dt ≥ Ṡ_{i,Σ} ≥ 0,     mit dS_Umgebung/dt ≥ J
⇒   dEnt_außen/dt ≥ dSyn_Σ/dt        (Ungleichung, nicht Gleichung)
```

Die Kreislaufbeziehung zwischen Entropie außen und Syntropie innen ist also eine **Ungleichung**. Ordnung innen ist immer durch mindestens ebenso viel Entropie außen bezahlt.

**Beispiel Erde.** Absorbiert werden 1,22·10¹⁷ W [G] bei 5 772 K [G], abgestrahlt wird bei 255 K [G]. Der Export beträgt J = (4/3)·P·(1/T_aus − 1/T_ein) = 6,07·10¹⁴ W/K, das sind 4,40·10³⁷ k/s. Diese Ordnungsbildung ist möglich, weil eine **Falte** (die Sonne) als heiße Quelle in einem **kalten Außenraum** steht. Die Verdichtung erzeugt den Temperaturgradienten, von dem jede Syntropie lebt.

**Rollen im Zustandsraum.** Die transversale Kontraktion ρ̇ = −κρ verkleinert das Phasenraumvolumen. In thermostatisierten Modellsystemen ist die Kontraktionsrate des Phasenraums gleich der exportierten Entropieproduktion, σ = k·⟨Σκ⟩ (Evans, Cohen und Morriss 1993; Gallavotti und Cohen 1995). **Syntropie ist damit die transversale Rolle:** Anziehung auf den Torus-Attraktor, bezahlt durch Export. Tangentiale Freiheit, also Umlauf und Phase, bleibt erhalten.

---

## 2 — Informationserhalt: Unitarität und Page-Kurve

Die Aussage des Modells lautet: „Der Kosmos verliert keine Information." In der heutigen Physik entspricht das der **Unitarität** der Zeitentwicklung. Ein reiner Gesamtzustand bleibt rein, die feinkörnige (von-Neumann-)Entropie des Ganzen ist konstant. Das steht nicht im Widerspruch zum zweiten Hauptsatz, denn dieser betrifft die **grobkörnige** Entropie. Die Information bleibt erhalten, sie wird nur unzugänglich verteilt.

**Page-Kurve.** Im einfachen Modell gilt S_BH(t) = S₀(1 − t/t_ev)^{2/3} und S_rad = β(S₀ − S_BH). Bei unitärer Entwicklung ist die feinkörnige Strahlungsentropie S_vN = min(S_rad, S_BH): Sie steigt, erreicht ein Maximum und fällt am Ende auf null.

| β | Page-Zeit t_P/t_ev | S_BH(t_P)/S₀ | M(t_P)/M₀ |
|---|---|---|---|
| 1 (Idealfall) | 0,6464 | 0,5000 | 0,7071 |
| 1,4847 (Photonen + Gravitonen, Page 2013) | **0,5381** | 0,5975 | 0,7730 |

Die Page-Zeit liegt also bei etwa der Hälfte der Verdampfungszeit. Das Modell reproduziert Pages Wert von 0,5381 t_ev. Beispiele: 1 M☉ ergibt t_Page = 1,13·10⁶⁷ a, M87* ergibt 3,10·10⁹⁶ a.

**Stand der Forschung.**
- Hawking (1976) folgerte aus der semiklassischen Rechnung einen Informationsverlust.
- Page (1993) zeigte, wie die Entropie bei Unitarität verlaufen *müsste*.
- Über AdS/CFT (Maldacena 2003) ist die Unitarität in Anti-de-Sitter-Räumen gut gestützt.
- Seit 2019 reproduzieren Inselformel und Replika-Wurmlöcher die Page-Kurve in semiklassischer Gravitation (Penington 2020; Almheiri, Engelhardt, Marolf und Maxfield 2019; Übersicht bei Almheiri et al. 2021).
- **Offen** bleibt der mikroskopische Mechanismus für reale, asymptotisch flache oder de-Sitter-artige Schwarze Löcher. „Keine Information geht verloren" ist die heute mehrheitlich gestützte Arbeitshypothese, aber kein bewiesener Satz.

**FTM-Lesart.** Die Falte ist kein Tresor, der Information verschluckt. Sie ist eine Umordnung, in der Information in Korrelationen zwischen Falte und Strahlung steht. Nach der Page-Zeit wird diese Information schrittweise an die Strahlung zurückgegeben. Das ist die quantenmechanische Fassung des Ω-Schritts (Abschnitt 4).

---

## 3 — Die globale Engstelle: das Loch des Torus als Horizont

### 3.1 Horizontentropien

| Horizont | Radius | S/k |
|---|---|---|
| Hubble-Radius heute, c/H₀ | 4,283 Gpc | 2,10·10¹²² |
| Ereignishorizont heute, R_CEH = c∫₁^∞ da/(a²H) | 15,93 Mrd Lj | 2,73·10¹²² |
| de-Sitter-Endzustand, c/H_DE | 5,119 Gpc | **3,00·10¹²²** (Maximum) |

Egan und Lineweaver (2010) geben S_CEH = 2,6·10¹²² k [G] an. Die eigene Rechnung liefert 2,73·10¹²² k und stimmt damit überein. Die Gibbons-Hawking-Temperatur des Endzustands beträgt T_GH = ħH_DE/(2πk) = 2,31·10⁻³⁰ K. Tiefer kann der Außenraum nicht abkühlen.

### 3.2 Die Brückenidentität: Das Ganze ist eine Falte mit r₀ = c/H

In einem flachen Universum bei kritischer Dichte ist die Masse in der Hubble-Kugel

```
M_H = (4π/3)·ρ_c·(c/H)³ = c³/(2GH) = 8,90·10⁵² kg
Δs_H = G·M_H/(3c²) = c/(6H)      ⇒      r₀(M_H) = 6Δs_H = c/H
S_Falte(M_H) = 36π·k·(Δs_H/ℓ_P)² = π·k·c²/(ℓ_P²H²) = S_Hubble        (Quotient 1,000000000000)
```

**Die Hubble-Kugel ist ihr eigener Grenzradius.** Die globale Horizontentropie ist exakt die Entropieformel der lokalen Falte, eingesetzt mit Δs_H = c/(6H). Das ist die Aussage „Das Ganze liegt auf seiner eigenen Schwelle", denn G·M_H/(c/H)² = cH/2 = g_H.

Ein Unterschied bleibt, und er ist physikalisch: die Oberflächengravitation.

| | Grenzradius | Temperatur | Oberflächengravitation |
|---|---|---|---|
| lokal | r₀ = 6Δs | T = ħc/(24πkΔs) | κ = c²/(2r₀) |
| global | r_H = 6Δs_H | T = ħc/(12πkΔs_H) | κ = c²/r_H |

Bei gleichem Radius ist der globale Horizont doppelt so heiß: 2,76·10⁻³⁰ K statt 1,38·10⁻³⁰ K. Die Falte ist nach innen geschlossen, der Horizont nach außen. Das ist die Umkehrung von Innen und Außen, die der Torus topologisch verlangt.

### 3.3 Das Universum ist keine Einzelfalte

Das Modell-Universum hat M_U/M_H = 1,686, also r₀(M_U) = 1,686·c/H₀ und S(M_U) = 1,99·S_dS. Das übersteigt das Maximum. In einem de-Sitter-Raum ist die größte mögliche Falte die Nariai-Lösung:

- Masse M_N = c³/(3√3·G·H_DE) = 4,09·10⁵² kg,
- Entropie S_dS/3 = 1,00·10¹²² k,
- beide Horizonte zusammen 2/3·S_dS.

Außerdem gilt T_H(M_U) < T_GH (Verhältnis 0,35). Eine Falte dieser Größe gäbe nie etwas zurück. **Das Loch des Torus ist deshalb nicht eine riesige Falte, sondern eine Struktur eigener Art:** eine Horizontfläche, die kein Objekt einschließt. Sie hängt vom Beobachter ab, denn jeder Beobachter hat seinen eigenen de-Sitter-Horizont. Das passt zum Grundsatz des Modells: kein ausgezeichnetes Zentrum, jeder Wirbel hat sein eigenes Zentrum. Nach der Beobachter-Regel ist die globale Engstelle Geometrie der Schließung und kein Ding.

### 3.4 Haushalt: Die Schließung setzt die Obergrenze, die Falten füllen sie

| Größe | S/k |
|---|---|
| Sonne (Materie) | 2,6·10⁵⁸ |
| CMB | 5,3·10⁸⁹ |
| S_obs heute, dominiert von supermassereichen Falten (Egan und Lineweaver 2010) | 3,1·10¹⁰⁴ [G] |
| Obergrenze Σ S_Falten, wenn alle Materie der Hubble-Kugel in einer Falte läge, Ω_m²·S_H | 1,89·10¹²¹ |
| S_CEH | 2,7·10¹²² |
| S_dS (Maximum) | 3,0·10¹²² |

- **Füllgrad heute:** S_obs/S_dS = 1,03·10⁻¹⁸.
- Die 3,1·10¹⁰⁴ k entsprechen 10¹¹ Falten mit √⟨M²⟩ = 1,7·10⁸ M☉. Das zeigt, dass der Haushalt von den größten Falten getragen wird. Eine einzelne M87*-Falte liefert 1,4·10⁻⁸ von S_obs.
- Wegen S ∝ M² ist Σ S_i ≤ 4πG(Σ M_i)²/(ħc). Verschmelzen erhöht die Entropie. Die Falten allein erreichen das globale Maximum nie: höchstens S_dS/3 (Nariai), realistisch höchstens Ω_m²·S_H.

**Brücke (präzise).** Σ S_Falten → S_global gilt **nicht** durch Verdichtung allein. Das Maximum wird erst erreicht, wenn die Falten über die Rückgabe 𝓡 zerstrahlen und die Strahlung über den Horizont abfließt. Nach dem verallgemeinerten zweiten Hauptsatz (Bekenstein 1974, Gibbons und Hawking 1977) wächst S_Materie + S_Falten + S_Horizont monoton, und der Endzustand ist S → S_dS. Die lokalen Falten sind **Zwischenspeicher auf dem Weg zur globalen Schließung**.

---

## 4 — Zusammenführung mit dem Ω-Operator

### 4.1 Der Ω-Operator physikalisch: Hawking-Rückgabe 𝓡 als Abbau des Flussmonopols

Im FGMS ist die Zentrierung c_z = maxᵢ Σⱼ F_ij das Flussmonopol einer Region. Ω zieht es in Richtung Gleichverteilung U. Das physikalische Gegenstück ist die Rückgabe 𝓡 der in der Falte konzentrierten Energie Mc² als isotrope Strahlung:

```
dM/dt = −A/M²,     A = ħc⁴/(15360π·G²) = 3,96·10¹⁵ kg³/s
relative Rate  |Ṁ/M| = A/M³ = 1/(3·t_ev)
```

Im Vergleich mit dem FGMS ergeben sich zwei Befunde.

**Befund 𝓡-1 (Persistenz).** Persistenz und Zentrierung wachsen korreliert. Physikalisch gilt t_ev ∝ M³, die Persistenz wächst also streng mit der Zentrierung. Der Fixierungsindex 𝒟 des Ω-Operators ist für große Falten tatsächlich groß.

**Befund 𝓡-2 (Rate).** Die Rate von 𝓡 **fällt** mit der Zentrierung wie M⁻³:

| Falte | \|Ṁ/M\| [1/a] | L_H [W] |
|---|---|---|
| 1 M☉ | 1,59·10⁻⁶⁸ | 9,0·10⁻²⁹ |
| Sgr A* | 2,22·10⁻⁸⁸ | 5,2·10⁻⁴² |
| M87* | 5,79·10⁻⁹⁸ | 2,1·10⁻⁴⁸ |

Die Natur baut kleine Monopole zuerst ab, große zuletzt. Der Gain des FGMS, G^Ω = A_Ω·sigmoid((𝒟 − 𝒟_c)/τ_s), ist daran anzupassen (Abschnitt 6, Punkt 1).

**Auslöser.** 𝓡 arbeitet nur, wenn der Außenraum kälter ist als die Falte:

```
𝓡 aktiv   ⇔   T_H(M) > T_außen(a) = T_CMB/a     (Boden: T_GH)
```

Heute gilt T_H = T_CMB bei M = 4,50·10²² kg (0,61 Mondmassen). **Jede astrophysikalische Falte ist heute kälter als der CMB und wächst netto.**

| Falte | Start der Nettoverdampfung bei a | t_start [a] |
|---|---|---|
| 10 M☉ | 4,4·10⁸ | 3,4·10¹¹ |
| Sgr A* | 1,8·10¹⁴ | 5,6·10¹¹ |
| M87* | 2,9·10¹⁷ | 6,8·10¹¹ |
| Milchstraße | 2,7·10¹⁸ | 7,2·10¹¹ |

Nach wenigen 10¹¹ Jahren hat die Expansion den Außenraum so weit gekühlt, dass die Rückgabe startet. Verglichen mit t_ev ist das verschwindend früh. **Das Atmen hält den Außenraum kalt und schaltet damit 𝓡 ein.**

**Zur Zuordnung „Entropie ↔ Diffusion" im FGMS.** Die Zuordnung gilt für nicht gravitierende Freiheitsgrade. Für gravitierende kehrt sich das Vorzeichen um. Die thermodynamische Entropie wächst in **beiden** Richtungen des Kreislaufs:

- bei der Verdichtung (c_z steigt): Stern → Falte, Faktor 4,1·10¹⁸,
- bei 𝓡 (c_z fällt): S_rad = β·ΔS_BH mit β ≈ 1,48 > 1 (Zurek 1982; Page 2013).

Ob verdichtet oder entzentriert wird, entscheidet nicht die Entropie, sondern die **Temperaturordnung** T_H gegen T_außen.

### 4.2 Der geschlossene Kreislauf und seine Zeitskalen

```
 (1) Verdichtung        Wolke → Stern → Falte           t_ff = 3,4·10⁵ a,  t_KH = 3,1·10⁷ a
         │              Akkretion (Salpeter, ε = 0,1)   t_Sal = 5,0·10⁷ a
         │              S steigt um ~10¹⁹ (Penrose), solange T_H < T_außen: Falte wächst
         ▼
 (2) Engstelle r₀ = 6Δs S = 36πk(Δs/ℓ_P)² maximal pro Masse; Form = (M, J, Q_el)
         │              Information in Korrelationen, nicht vernichtet
         ▼
 (3) Rückgabe 𝓡         Start bei T_außen < T_H:  t_start ≈ 3·10¹¹ – 7·10¹¹ a
         │              Dauer t_ev = 2,1·10⁶⁷ a (1 M☉) … 4,5·10⁹⁹ a (6·10¹⁰ M☉) … 2,1·10¹⁰⁰ a (10¹¹ M☉)
         │              Page-Zeit 0,5381·t_ev: Information kehrt in die Strahlung zurück
         ▼
 (4) Verteilung         Strahlung im atmenden Außenraum: T ∝ 1/a, Energie ∝ 1/a, Entropie bleibt
         │              Abfluss über den Horizont (globale Engstelle), Boden T_GH = 2,3·10⁻³⁰ K
         ▼
 (5) erneute Verdichtung  nur erlaubt, wenn  N ≥ β·z_f²   (siehe 4.3)
```

Die Zeitskalen liegen um 60 bis 92 Größenordnungen auseinander. Verdichtung ist schnell, Rückgabe extrem langsam. Das Monopol wird also nicht „abgebaut", sobald es entsteht. Es besteht über Zeiträume, gegen die das heutige Weltalter nur den Bruchteil 6,6·10⁻⁵⁸ (1 M☉) bis 3,0·10⁻⁹⁰ (Milchstraße) ausmacht (t_ev/t₀ = 1,52·10⁵⁷ bzw. 3,28·10⁸⁹).

### 4.3 Verträglichkeit mit dem zweiten Hauptsatz

**(a) Rückverdichtung.** Angenommen, N Falten der Masse M zerstrahlen. Die Strahlung trägt die Entropie N·β·S₀ und die Energie NMc². Bis zur neuen Verdichtung wird sie um z_f = a₂/a₁ rotverschoben. Eine neue Falte aus dieser Energie hätte S = (N/z_f)²·S₀. Der zweite Hauptsatz erlaubt das nur, wenn

```
N ≥ β·z_f²          z_f = 1: N ≥ 1,485;   z_f = 10: N ≥ 148,5;   z_f = 10³: N ≥ 1,49·10⁶
```

Die Rückverdichtung ist nur **hierarchisch aufwärts** möglich: mehr Masse in weniger, größeren Falten. Das fügt sich in das fraktale Bild. Expansion erschwert sie mit z_f². Nur eine Kontraktionsphase (z_f < 1) erleichtert sie. Das ist genau das Tolman-Problem zyklischer Modelle: Jeder Zyklus beginnt mit mehr Entropie (Tolman 1934).

**(b) Ohne Atmen gibt es keinen Kreislauf.** In einem festen Kasten ist eine Falte mit ihrer eigenen Strahlung stabil im Gleichgewicht, wenn E_rad < Mc²/4 (Hawking 1976). Für 1 M☉ gilt das bis V_max = 4,1·10⁹⁰ m³, das ist weit mehr als das beobachtbare Universum. Dieser Zustand ist maximale Entropie und stationär, aber tot. Erst wenn V wächst und T_außen fällt, wird das Gleichgewicht instabil. Die Falte hat die negative Wärmekapazität C = −8πGM²k/(ħc) = −2,9·10⁵⁴ J/K für 1 M☉, und die Rückgabe setzt ein. **Das Atmen ist die thermodynamische Voraussetzung dafür, dass der Kreislauf überhaupt läuft.**

**(c) Der Kreislauf ist eine Spirale.** In jeder Runde wächst S_ges. Stationär ist nur die **Form**: Verdichtung und Verteilung wechseln sich ab, und Ordnung besteht nur durch Durchfluss (FGMS: „Ordnung existiert nur durch Durchfluss"). Stationär ist nicht der **Zustand**. Der Kreislauf ist ein **offener Durchfluss-Kreislauf mit globaler Senke**.

### 4.4 Bedingung für einen stationären Syntropie-Kreislauf

Ein Teilsystem Σ (lokale Falten und die Strukturen in ihrem Strahlungsstrom) in einem kausalen Gebiet mit Horizont ist stationär syntrop, wenn gleichzeitig gilt:

```
(S1) Temperaturordnung   T_Falte  >  T_Struktur  >  T_außen(a)  ≥  T_GH
(S2) Export              J = Ṡ_{i,Σ} ≥ 0                  (Prigogine: d_eS = −d_iS)
(S3) Senkenkapazität     dS_max/dt ≥ Ṡ_i,   S_max = π·k·c²/(ℓ_P²H²)
                         ⇔   −Ḣ  ≥  ℓ_P²·H³·Ṡ_i / (2π·c²·k)
(S4) Rückverdichtung     N ≥ β·z_f²
```

(S3) ist die mathematische Fassung von „das Atmen hält den Außenraum kalt": **Die globale Engstelle, also die Horizontfläche, muss mindestens so schnell wachsen, wie lokal Entropie erzeugt wird.**

Zahlen:
- Heutige Produktion: Ṡ_i ≈ S_obs/t₀ = 7,1·10⁸⁶ k/s.
- Benötigt wird −Ḣ ≥ 3,8·10⁻⁵⁴ s⁻². ΛCDM liefert heute −Ḣ = 1,5·H₀²·Ω_m = 2,3·10⁻³⁶ s⁻². Die Reserve beträgt 6·10¹⁷.
- In ΛCDM geht Ḣ gegen null. Bei konstantem σ verletzt der Kreislauf (S3) ab a ≈ 1,0·10⁶, das heißt bei t ≈ 2,4·10¹¹ a. Danach zehrt er von der Restkapazität S_dS − S. Bei heutiger Rate wäre sie in 1,3·10²⁸ a erschöpft. Tatsächlich sinkt σ, und der Endzustand ist das de-Sitter-Gleichgewicht.

**Folgerung.** Mit Λ_kosm = konstant läuft jeder Syntropie-Kreislauf in den Wärmetod des de-Sitter-Gleichgewichts. Die Aussage des Modells „ein atmender Kosmos läuft nicht in einen Endzustand" verlangt, dass **H weiter fällt**, dass die Horizontfläche also unbegrenzt wächst: w ≠ −1. Die Forderung besteht, ist aber unmessbar: Im Modell ist sie die Mindestrate neuer innerer Faltung, Λ̇/Λ ≥ Ṡ_i/(2S_max) (Ṡ_i die Entropieproduktion Ṡ_i), heute 3,7·10⁻³⁷ s⁻¹, rund 10⁻²⁹ je Jahr; das gibt heute w_s = −1 + 1·10⁻¹⁹, im Endzustand (S_max = 9,6·10¹²² k, H_end = 0,49·H₀) 2·10⁻¹⁹. [korr.] Die Hinweise auf veränderliche Dunkle Energie in den BAO-Daten von DESI (2024/2025) [G] erklärt sie nicht. [korr.] Eine an H gebundene Schwelle a₀ = cH/2π ∝ H trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17). Im Modell ist Λ die ganze innere Faltung der Zelle und heute fest; a₀ ist seit z ≈ 10 fest (bei der Rekombination 0,76·a₀), ≈ cH₀/2π gilt nur heute. [korr.]

---

## 5 — Ergebnis: Die Engstelle ist lokal UND global

| | lokal: Falte | global: Loch des Torus |
|---|---|---|
| Ort | r ≤ r₀ = 6Δs, pro Wirbel | Horizont r_H = c/H = 6Δs_H, pro Beobachter |
| Faltungstiefe | Δs = GM/(3c²) | Δs_H = c/(6H) |
| Entropie | 36πk(Δs/ℓ_P)² | 36πk(Δs_H/ℓ_P)² = πkc²/(ℓ_P²H²) |
| Temperatur | ħc/(24πkΔs) | ħc/(12πkΔs_H) |
| Oberflächengravitation | c²/(2r₀) (nach innen geschlossen) | c²/r_H (nach außen geschlossen) |
| Zahlen | 10⁷⁷ (1 M☉) bis 3,8·10⁹⁸ k (Milchstraße) | S_dS = 3,0·10¹²² k |
| Rolle | Zwischenspeicher, Quelle des Temperaturgradienten | Obergrenze und endgültige Senke |
| Zeitskala | t_ev = 3840·(S/k)·Δs/c | Füllzeit (S_dS − S)/σ |
| im FGMS | hohes c_z, lange Persistenz | Gleichverteilung U, Viabilitätsgrenze |

**Die Brücke** hat drei Teile:
1. **Formal:** dieselbe Entropieformel, Δs → Δs_H. Die Hubble-Kugel bei kritischer Dichte ist eine Falte mit r₀ = c/H (Quotient exakt 1).
2. **Quantitativ:** Σ S_Falten ≤ Ω_m²·S_H < S_dS/3 < S_dS. Die Falten füllen den Haushalt heute zu 10⁻¹⁸ und können ihn allein nie füllen.
3. **Dynamisch:** Der Weg zu S_global führt **durch** die lokale Engstelle. Verdichtung konzentriert, 𝓡 gibt zurück, der Horizont nimmt auf.

Kurzfassung: **Die globale Schließung setzt das Maß, die lokalen Falten sind die Stationen auf dem Weg dorthin.**

---

## 6 — Offene Punkte (konstruktiv)

1. **Ω-Gain an die Physik anpassen.** Vorschlag für die physikalische Lesart des FGMS: G^Ω(i) ∝ Θ(T_H(cᵢ) − T_außen)·cᵢ⁻³, also eine Rate, die der Temperatur folgt. Prüfen, ob die Kopplung Persistenz × Zentrum als *Auslöser* erhalten bleiben kann, während die *Rate* der Temperatur folgt.
2. **Zwei Entropien im FGMS trennen.** Die Verteilungsentropie des Flusses H(F), maximal bei U, und die gravitative Entropie, maximal bei maximaler Zentrierung, sollten zwei getrennte Größen werden. Die Zuordnung „Entropie ↔ Diffusion" gilt nur für die erste.
3. **Faltungsbilanz.** Δs_ges = konstant gilt für Materie. Für Strahlung nach 𝓡 gilt d(Δs_rad)/dt = −H·Δs_rad, weil die Energie mit 1/a fällt und die ART keine globale Energieerhaltung kennt. Dieser Term ist nicht die Rückstellung: Die Bilanz des Autors ist σ = ΣΔs_innen − ΔL_außen = 0. Der Verlust der Strahlung durch Rotverschiebung senkt Λ (heute w_s = −1,0002); dieser Anteil verschwindet mit 1/a. [korr.]
4. **Stationarität (S3) quantitativ.** Welches H(t) fordert das FTM? Aus (S3) ergibt sich eine Mindestabnahme von H; im Modell ist sie die Mindestrate neuer innerer Faltung, Λ̇/Λ ≥ Ṡ_i/(2S_max) (Ṡ_i = Ṡ_i), mit heute w_s = −1 + 1·10⁻¹⁹, also nicht an DESI- oder Supernova-Daten prüfbar. [korr.] Eine an H gebundene Schwelle (a₀ = cH/2π) trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie ΔN_eff = 2,09 (gemessen N_eff = 2,99 ± 0,17); a₀ ist seit z ≈ 10 fest. [korr.]
5. **Kontraktionsphase.** Falls „Atmen" echte Kontraktion einschließt, trifft das Modell auf das Tolman-Problem: Die Horizontfläche schrumpft, gegen den verallgemeinerten zweiten Hauptsatz. Zu klären ist, ob das FTM-Atmen nur eine Modulation der Rate (H > 0, Ḣ schwankend) ist. Das wäre verträglich. Eine echte Kontraktion wäre es nicht ohne zusätzlichen Mechanismus. Penroses konforme zyklische Kosmologie (2010) ist ein bekannter, umstrittener Versuch.
6. **Information am globalen Horizont.** Ob die Information, die über den de-Sitter-Horizont abfließt, erhalten bleibt (de-Sitter-Komplementarität), ist offen. Hinzu kommen die Probleme von Poincaré-Wiederkehr und Fluktuationen im de-Sitter-Gleichgewicht (Dyson, Kleban und Susskind 2002).
7. **Verschränkung und Page-Kurve.** Der gemeinsame Knoten bei r ≤ r₀ mit d_innen = 0 und die Inselformel, in der das Innere der Falte zur Strahlung „gehört", sind strukturell verwandt. Eine formale Verbindung wäre ein eigenes Forschungsthema.

---

## Literatur

- Almheiri, A., Engelhardt, N., Marolf, D., Maxfield, H. (2019): The entropy of bulk quantum fields and the entanglement wedge of an evaporating black hole. JHEP 12, 063.
- Almheiri, A., Hartman, T., Maldacena, J., Shaghoulian, E., Tajdini, A. (2021): The entropy of Hawking radiation. Rev. Mod. Phys. 93, 035002.
- Bekenstein, J. D. (1973): Black holes and entropy. Phys. Rev. D 7, 2333.
- Bekenstein, J. D. (1974): Generalized second law of thermodynamics in black-hole physics. Phys. Rev. D 9, 3292.
- Dyson, L., Kleban, M., Susskind, L. (2002): Disturbing implications of a cosmological constant. JHEP 10, 011.
- Egan, C. A., Lineweaver, C. H. (2010): A larger estimate of the entropy of the universe. ApJ 710, 1825.
- Evans, D. J., Cohen, E. G. D., Morriss, G. P. (1993): Probability of second law violations in shearing steady states. Phys. Rev. Lett. 71, 2401.
- Gallavotti, G., Cohen, E. G. D. (1995): Dynamical ensembles in nonequilibrium statistical mechanics. Phys. Rev. Lett. 74, 2694.
- Gibbons, G. W., Hawking, S. W. (1977): Cosmological event horizons, thermodynamics, and particle creation. Phys. Rev. D 15, 2738.
- Hawking, S. W. (1975): Particle creation by black holes. Commun. Math. Phys. 43, 199.
- Hawking, S. W. (1976): Black holes and thermodynamics. Phys. Rev. D 13, 191.
- Hawking, S. W. (1976): Breakdown of predictability in gravitational collapse. Phys. Rev. D 14, 2460.
- Maldacena, J. (2003): Eternal black holes in anti-de Sitter. JHEP 04, 021.
- Nicolis, G., Prigogine, I. (1977): Self-Organization in Nonequilibrium Systems. Wiley.
- Page, D. N. (1976): Particle emission rates from a black hole. Phys. Rev. D 13, 198.
- Page, D. N. (1993): Information in black hole radiation. Phys. Rev. Lett. 71, 3743.
- Page, D. N. (2013): Time dependence of Hawking radiation entropy. JCAP 09, 028.
- Penington, G. (2020): Entanglement wedge reconstruction and the information paradox. JHEP 09, 002.
- Penrose, R. (1979): Singularities and time-asymmetry. In: Hawking, S. W., Israel, W. (Hrsg.), General Relativity: An Einstein Centenary Survey. Cambridge University Press.
- Penrose, R. (2010): Cycles of Time. Bodley Head.
- Schrödinger, E. (1944): What is Life? Cambridge University Press.
- Tolman, R. C. (1934): Relativity, Thermodynamics and Cosmology. Oxford University Press.
- Zurek, W. H. (1982): Entropy evaporated by a black hole. Phys. Rev. Lett. 49, 1683.
