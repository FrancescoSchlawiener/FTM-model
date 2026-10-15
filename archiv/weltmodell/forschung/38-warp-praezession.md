# Warp-Präzession der Milchstraßenscheibe

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Dynamik aus dem Energiefunktional (Kraft −∇Φ der Feldform, 3.2) → die Bahnen folgen der örtlichen Regel, die Metrik trägt deren wirbelfreien Teil; Warp-Zahlen −1 bis −16 km/s/kpc als Feldform-Zahlen markiert; Sonnenort „FTM (volle Feldform) 94,6 M☉/pc²“ → örtliche Regel 71 M☉/pc². Gedächtnis-Nachlauf 1/λ_m = 170–260 Myr in der Milchstraße (Abschnitt 4, Offen 3) → in Galaxien kein messbarer Nachlauf, das Gedächtnis wirkt auf den Zustand der Regel. Äußeres Feld „≤ 4·10⁻³ km/s/kpc“ → 0,02 km/s/kpc, Magellansche Wolke bei 14 kpc 0,32 → 0,36 km/s/kpc (spätere Rechnung mit verstärkten Feldern); dazu: einfallendes Gas dreht die Knotenlinie höchstens um 90°. Rückkopplungsanteil ρ_ν, ρ₂ als gelesene Dichte gekennzeichnet. Im τ_M = r*/(2πR) als Punkt 1 des Autors gekennzeichnet. Umbenannt: Spinzahl χ = 1362 → j = 1361; Abplattung q → ε_a; Wachstumsvariable s (Abschnitt 4) → p.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: rechnungen/warp_praezession.py
> Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Vergleich mit dem Modellwert cH₀/2π = 1,082·10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc) in 3.4.
> Werte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft. Vorzeichen: Präzession > 0 heißt rechtläufig (im Drehsinn der Scheibe).

---

## Kurzfassung

| Frage | Ergebnis bei R = 12–16 kpc |
|---|---|
| Messwert (Gaia, Poggio et al. 2020 [G]) | **+10,9 ± 0,1 (stat) ± 3,2 (sys) km/s/kpc**, rechtläufig |
| Freie Präzession im Feld der baryonischen Milchstraße, Außenscheibe starr gekippt | **−1,0 bis −1,2 km/s/kpc** (Newton mit denselben Baryonen: −1,1 bis −2,2); Feldform-Zahl [korr.] |
| dasselbe, einzelner Ring von 2 kpc Breite gekippt | **−12,4 bis −16,1 km/s/kpc**; Feldform-Zahl [korr.] |
| Anteil der ruhenden Feldgeometrie allein | −17 bis −31 km/s/kpc; die mitgeführte Feldgeometrie hebt davon 50–95 % auf |
| Gedächtnis-Nachlauf 1/λ_m = 170–260 Myr | langsame Mode: Präzession ≈ 0 (−0,04 bis −0,3), Kippung klingt ab; schnelle Mode: −19,6 bzw. −29 km/s/kpc, klingt in 0,2–0,5 Gyr ab. Gilt nur unter der Annahme eines Nachlaufs in der Galaxie; dort gibt es keinen messbaren. [korr.] |
| Äußeres Feld der Lokalen Gruppe (g_e = 0,01–0,03 a₀) | Drehmoment erst in Ordnung g_e²; **≤ 0,02 km/s/kpc** (spätere Rechnung mit verstärkten Feldern; hier ≤ 4·10⁻³) [korr.] |
| Große Magellansche Wolke (3,2·10⁹ M☉ baryonisch [G]) | Feldgleichung verstärkt die Gezeit um den Faktor 6,4; Präzession ≤ 0,04, Knotenverschiebung ≤ 0,37 km/s/kpc (spätere Rechnung mit verstärkten Feldern bei 14 kpc: 0,36) [korr.] |
| Verwindung a_J = J/(Mc) = 4,35 pc (Lense-Thirring) | 3·10⁻⁶ km/s/kpc (mit Faktor ν: 6·10⁻⁶) |
| Knotenzahl ζ = ω_p/Ω | FTM (Feldform) [korr.]: −0,06 bis −0,09 (starr), −0,89 bis −0,94 (Ring); Gaia verlangt +0,61 bis +0,83 |

**Ergebnis.** Im FTM präzediert ein geneigter Teil der baryonischen Scheibe in allen gerechneten Lesarten rückläufig oder gar nicht; die Beträge sind in der Feldform gerechnet. [korr.] Äußeres Feld, Magellansche Wolke und Verwindung liefern zusammen weniger als 0,4 km/s/kpc. Schräg einfallendes Gas dreht die Knotenlinie höchstens um 90°, bei 10,9 km/s/kpc nach 0,14 Mrd. Jahren. [korr.] Zum Gaia-Wert fehlen damit mindestens 10,5 km/s/kpc, und das Vorzeichen der inneren Dynamik ist entgegengesetzt. Der Grund ist strukturell: Die Feldgeometrie des Operators folgt der abgeplatteten Baryonenverteilung. Die ausgerichtete Lage ist deshalb ein Minimum der Energie, und daraus folgt rückläufige Präzession.

---

## 1. Messgröße

Die Scheibe der Milchstraße ist jenseits von etwa 10–12 kpc verbogen [G]. Poggio et al. (2020) haben aus Gaia-Eigenbewegungen von Riesensternen ein kinematisches Modell angepasst. Darin ist die Höhe der Scheibe

$$
h(R,\varphi,t)=h_w(R)\,\sin\!\big(\varphi-\varphi_w(t)\big),\qquad \dot\varphi_w=\omega_p ,
$$

und die Vertikalgeschwindigkeit der Sterne ist v_z = ∂h/∂t + Ω ∂h/∂φ = h_w(Ω − ω_p)cos(φ − φ_w) (bei fester Amplitude). Das Ergebnis war ω_p = +10,9 km/s/kpc. Mit dem FTM-Ω bei 14 kpc (15,2 km/s/kpc) ist das ω_p/Ω = +0,72.

Die Warp-Präzession ist die Knotenbewegung eines gekippten Scheibenteils. Gesucht ist also das Drehmoment, das die innere Scheibe, das äußere Feld und die Rotation auf einen gekippten Außenteil ausüben, im nichtlinearen Feld des FTM.

---

## 2. Das Feld der baryonischen Milchstraße

### 2.1 Massenmodell (McMillan 2017 [G], nur Baryonen)

| Teil | Dichte | Parameter | M [M☉] |
|---|---|---|---|
| Bulge | ρ₀(1+r′/r₀)^{−α} e^{−(r′/r_c)²}, r′² = R² + (z/ε_a)² [korr.] | r₀ = 0,075, r_c = 2,1 kpc, α = 1,8, ε_a = 0,5 | 9,23·10⁹ |
| dünne Scheibe | Σ₀e^{−R/R_d} e^{−\|z\|/z_d}/(2z_d) | Σ₀ = 896 M☉/pc², R_d = 2,50, z_d = 0,30 kpc | 3,52·10¹⁰ |
| dicke Scheibe | ebenso | Σ₀ = 183, R_d = 3,02, z_d = 0,90 | 1,05·10¹⁰ |
| HI | Σ₀e^{−R_m/R−R/R_d} sech²(z/2z_d)/(4z_d) | Σ₀ = 53,1, R_d = 7, R_m = 4, z_d = 0,085 | 1,07·10¹⁰ |
| H₂ | ebenso | Σ₀ = 2180, R_d = 1,5, R_m = 12, z_d = 0,045 | 1,23·10⁹ |
| **Summe** | | | **M_b = 6,68·10¹⁰** |

Damit ist r* = √(GM_b/a₀) = 8,81 kpc.

### 2.2 Feldgleichung und Lösungsweg

Im schwachen Feld gilt

$$
\nabla^2\Phi=\nabla^2\Phi_N+\nabla\!\cdot\!\big[(\nu(y)-1)\nabla\Phi_N\big],\qquad y=\frac{|\nabla\Phi_N|}{a_0},\qquad \nu=\frac1{1-q},\quad q=e^{-\sqrt y}.
$$

Die Sterne der Scheibe sind über ihre Gezeitengrenzen (Größenordnung pc) geglättet. Auf der Skala der Galaxie ist Φ̄_N daher das Newton-Potential der glatten Verteilung.

Gelöst wird in zwei Schritten auf einem axialsymmetrischen Finite-Volumen-Gitter (R, z ≥ 0; 435 × 240 Zellen; ΔR = 0,1 kpc bis 30 kpc, Δz = 0,01 kpc bis 0,4 kpc, danach geometrisch wachsend bis 1,5 Mpc):

1. **Newton-Teil.** ∇²Φ_N = 4πGρ mit Monopol und Quadrupol als Randwert.
2. **Rückkopplungsteil** Φ_ν = Φ − Φ_N. Die Quelle ist die diskrete Divergenz des Flusses F = (ν−1)∇Φ_N an den Zellflächen. Am Außenrand wird derselbe Fluss als Neumann-Bedingung gesetzt. Dadurch fallen die Randflächen exakt heraus, und das System ist von selbst verträglich. Der rotationsfreie Teil von ∇Φ_ν ist F; der Wirbelanteil ist in der Lösung enthalten.

**Kontrollen:** Homogene Kugel: dΦ_N/dR auf 2·10⁻⁴ bis 5·10⁻⁴. Feldform außerhalb der Kugel: dΦ/dR = ν·g_N auf 1–4·10⁻⁴. Gittermassen der Scheiben gegen 2πΣ₀R_d² auf 10⁻⁴.

### 2.3 Rotationskurve und vertikale Dynamikmasse

| R [kpc] | v_Newton [km/s] | v_FTM [km/s] | g_N/a₀ | ν = v²/v_N² | beobachtet [G] |
|---|---|---|---|---|---|
| 4 | 185,7 | 205,4 | 2,33 | 1,22 | – |
| 8,2 | 177,4 | 220,0 | 1,04 | 1,54 | 229 ± 1 (Eilers et al. 2019) |
| 12 | 157,1 | 216,2 | 0,56 | 1,89 | ≈ 222 |
| 14 | 147,1 | 213,0 | 0,42 | 2,10 | ≈ 218 |
| 16 | 138,2 | 210,1 | 0,32 | 2,31 | ≈ 214 |
| 20 | 124,0 | 205,3 | 0,21 | 2,74 | ≈ 209 |

Die Rotationskurve des FTM mit diesen Baryonen liegt 2–4 % unter den Messwerten.

Die vertikale Dynamikmasse Σ(|z| < 1,1 kpc) am Sonnenort beträgt:
- baryonisch 54,7 M☉/pc²,
- Newton 60,7 M☉/pc²,
- FTM (volle Feldform) 94,6 M☉/pc²; die Bahnen folgen aber der örtlichen Regel a = ν(|g_N|)·g_N, und die gibt am Sonnenort 71 M☉/pc² [korr.],
- gemessen 71 ± 6 M☉/pc² (Kuijken & Gilmore 1991 [G]).

Eine Variante mit um 20 % leichteren Sternscheiben ergibt 84,2 M☉/pc² bei v(8,2) = 206,9 km/s. Die Präzessionsraten unten ändern sich zwischen beiden Modellen um weniger als 1 km/s/kpc (3.4).

---

## 3. Präzession eines geneigten Scheibenteils

### 3.1 Einzelring im vollen Feld

Für einen dünnen, leicht gekippten Ring in einem festen axialsymmetrischen Potential gilt

$$
\omega_p=\Omega-\nu_z\;\approx\;-\frac{\nu_z^2-\Omega^2}{2\Omega},\qquad \Omega^2=\frac{\Phi_R}{R},\quad \nu_z^2=\Phi_{zz}\big|_{z=0},\quad \kappa^2=\Phi_{RR}+3\Omega^2 .
$$

Werte in km/s/kpc:

| Form | R | Ω | κ/Ω | ν_z/Ω | Ω − ν_z | −(ν_z²−Ω²)/2Ω |
|---|---|---|---|---|---|---|
| Newton | 12 | 13,09 | 1,092 | 3,475 | −32,4 | −72,5 |
| Newton | 14 | 10,50 | 1,048 | 3,474 | −26,0 | −58,1 |
| Newton | 16 | 8,64 | 1,024 | 3,487 | −21,5 | −48,2 |
| FTM | 12 | 18,02 | 1,351 | 3,479 | −44,7 | −100,1 |
| FTM | 14 | 15,22 | 1,339 | 3,474 | −37,6 | −84,2 |
| FTM | 16 | 13,13 | 1,337 | 3,484 | −32,6 | −73,1 |

Zwei Befunde:

1. **ν_z/Ω ist in Newton und FTM gleich** (3,47–3,49). Die Feldform erhält das Verhältnis von Vertikal- und Umlauffrequenz in der Mittelebene; ν_z² und Ω² werden beide etwa mit ν multipliziert. Die volle numerische Lösung bestätigt das auf 0,1 %.
2. **Einzelne Sternbahnen bilden keinen gekippten Ring.** Bei ν_z ≈ 3,5 Ω schwingt ein Stern pro Umlauf 3,5-mal durch die Scheibe. ν_z² wird von der örtlichen Scheibendichte beherrscht. Ein Warp ist deshalb eine *gemeinsame* Verbiegung, bei der die Nachbarringe mitkippen und die örtliche Dichte nicht zurückstellt. Die Tabelle ist eine formale Obergrenze des Betrags. Die maßgebliche Rechnung folgt in 3.2.

### 3.2 Kippsteifigkeit aus dem Energiefunktional

[korr.] Diese Rechnung setzt die Feldform ∇²Φ = ∇·[ν∇Φ_N] als Bahnregel voraus (Dynamik aus dem Energiefunktional, Kraft −∇Φ); das trägt nicht: Die Bahnen folgen der örtlichen Regel a = ν(|g_N|)·g_N, die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Für Kugeln sind beide gleich. Die Präzessionsraten unten (−1 bis −16 km/s/kpc) sind Feldform-Zahlen; für die örtliche Regel sind sie nicht gerechnet.

Die potentielle Energie des FTM im schwachen Feld ist ein Funktional der Newton-Felder allein:

$$
V[\rho]=-\frac1{8\pi G}\int\Big[|\nabla\Phi_N|^2+a_0^2\,Q_s\!\Big(\frac{|\nabla\Phi_N|^2}{a_0^2}\Big)\Big]d^3x,\qquad Q_s'=\nu-1 .
$$

Ihre Variation ist δV = ∫(Φ_N + Φ_ν)δρ d³x. In der Feldform ist die Kraft auf die Masse also −∇Φ mit dem vollen Potential der Feldgleichung. [korr.]

**Kippung.** Die Masse wird zerlegt in einen Innenteil ρ_i (Bulge und Scheibe innerhalb R_w, bleibt eben) und einen Außenteil ρ_o, der starr um die x-Achse um θ gekippt wird. Der Erzeuger der Drehung ist D = −z∂_y + y∂_z. Damit gilt ρ_o^θ = e^{−θD}ρ_o und Φ_{N,o}^θ = e^{−θD}Φ_{N,o}.

Für kleine θ ist V = V₀ + ½Kθ². Das Drehmoment −Kθ steht längs der Knotenlinie. Aus dL/dt = τ folgt für den waagrechten Anteil von L_o

$$
\boxed{\;\omega_p=-\frac{K}{L_o}\;},\qquad L_o=\int\rho_o\,R\,v_\varphi\,d^3x .
$$

K > 0 (ausgerichtete Lage ist ein Energieminimum) bedeutet rückläufige Präzession.

**Zweite Ableitung.** Mit ψ = DΦ_{N,o} (Newton-Antwort auf die Kippung: Φ_N^{(1)} = −ψ) und dem linearisierten Fluss δ[(ν−1)∇Φ_N] ergibt sich

$$
K=K_\text{ruh}+K_\text{mit},
$$
$$
K_\text{ruh}=\int\rho_o\,\big\langle D^2(\Phi-\Phi_{N,o})\big\rangle\,d^3x ,
\qquad
K_\text{mit}=-\frac1{4\pi G}\int\Big[(\nu-1)\,\langle|\nabla\psi|^2\rangle+y\,\nu'(y)\,\langle(\hat n\cdot\nabla\psi)^2\rangle\Big]d^3x .
$$

Dabei ist n̂ = ∇Φ_N/|∇Φ_N|, und ⟨⟩ bezeichnet das Mittel über φ. Für axialsymmetrisches Ψ(R, z) und h = RΨ_z − zΨ_R gilt

$$
\langle D^2\Psi\rangle=-\frac{zh}{R}-\frac{zR}{2}\,\partial_R\!\Big(\frac hR\Big)+\frac R2\,\partial_z h ,
\qquad
\langle|\nabla\psi|^2\rangle=\tfrac12\Big(h_R^2+h_z^2+\frac{h^2}{R^2}\Big),\qquad
\langle(\hat n\cdot\nabla\psi)^2\rangle=\tfrac12\big(\hat n_Rh_R+\hat n_zh_z\big)^2
$$

(im letzten Ausdruck mit h aus Φ_{N,o}).

**Bedeutung der beiden Teile:**
- **K_ruh:** Der Außenteil dreht sich in der *ruhenden* Feldgeometrie des ebenen Zustands. Sein eigenes Newton-Feld dreht mit und ist abgezogen.
- **K_mit:** Die Feldgeometrie des Rückkopplungsteils wird mit dem Außenteil mitgedreht. Das ist die Antwort der nichtlinearen Feldgleichung auf die neue Lage der Baryonen. Es ist keine zusätzliche Quelle; die Baryonen bleiben die einzige Quelle.

**Grenzfälle** (analytisch und numerisch geprüft):
- Dünner Ring, äußeres Ψ: Es folgt ⟨D²Ψ⟩ = (R²/2)(ν_z² − Ω²) und damit ω_p = −(ν_z² − Ω²)/(2Ω).
- Kugelsymmetrisches Ψ: h ≡ 0, also K = 0. Numerisch ergibt sich K/K_typ = −4·10⁻⁵.
- Newton (ν ≡ 1): K = ∫ρ_i⟨D²Φ_{N,o}⟩, die gewöhnliche Wechselwirkung von Innen- und Außenteil.
- ν − 1 = c konstant: K = (1 + c)·K_Newton. Das beruht auf der Identität (1/4πG)∫⟨|∇ψ|²⟩ = ∫ρ_o⟨D²Φ_{N,o}⟩, die auf dem Gitter auf 4·10⁻⁴ (Außenteil) bzw. 3·10⁻³ (Ring) erfüllt ist.

### 3.3 Zahlen

ω in km/s/kpc; < 0 heißt rückläufig.

| gekippter Teil | M_o [M☉] | Newton | ruhende Geometrie | mitgeführte Geometrie | **FTM (Feldform)** [korr.] |
|---|---|---|---|---|---|
| Außenteil R > 10 | 1,25·10¹⁰ | −2,23 | −20,86 | +19,86 | **−1,00** |
| Außenteil R > 12 | 9,19·10⁹ | −1,76 | −19,69 | +18,50 | **−1,19** |
| Außenteil R > 14 | 6,94·10⁹ | −1,40 | −18,57 | +17,36 | **−1,21** |
| Außenteil R > 16 | 5,36·10⁹ | −1,14 | −17,43 | +16,25 | **−1,17** |
| Ring 11–13 | 2,73·10⁹ | −11,11 | −31,32 | +15,26 | **−16,06** |
| Ring 13–15 | 1,87·10⁹ | −9,04 | −29,23 | +15,24 | **−13,98** |
| Ring 15–17 | 1,34·10⁹ | −7,52 | −27,84 | +15,47 | **−12,37** |

Die beiden Grenzfälle klammern einen echten Warp ein. Dessen Neigung wächst stetig mit R:
- **Starr gekippter Außenteil:** Die Eigengravitation des Warps ist vollständig.
- **Einzelner Ring:** Die Nachbarn bleiben eben.

### 3.4 Empfindlichkeit

| Variante | Außenteil R > 12 | Ring 13–15 |
|---|---|---|
| a₀ = 1,2·10⁻¹⁰, McMillan | −1,19 | −13,98 |
| a₀ = cH₀/2π = 1,082·10⁻¹⁰ (H₀ = 70) | −1,16 | −13,47 |
| Sternscheiben × 0,8 (M_b = 5,77·10¹⁰) | −0,97 | −14,52 |

Das Ergebnis hängt weder merklich von der Schwelle noch vom Massenmodell ab.

### 3.5 Lesart

Die ruhende Feldgeometrie hält den gekippten Außenteil stark fest, mit −17 bis −31 km/s/kpc. Der Rückkopplungsteil des Feldes sitzt aber dort, wo die Baryonen sitzen: Seine Quelle ∇·[(ν−1)∇Φ_N] enthält den Term (ν−1)·4πGρ_b. Kippt der Außenteil, kippt diese Geometrie mit und hebt beim starren Außenteil 94 %, beim einzelnen Ring rund die Hälfte der Rückstellung auf.

Übrig bleibt eine Kippsteifigkeit mit K > 0 in jedem gerechneten Fall. Die ausgerichtete Lage ist ein Minimum von V, und das Vorzeichen der Präzession ist rückläufig. Das ist eine Eigenschaft der Regel, nicht eines zusätzlichen Körpers: Der Operator formt das Feld nach der abgeplatteten Baryonenverteilung. Eine vom Baryonenbild unabhängige Form, etwa eine gekippte oder gestreckte, hat er nicht.

---

## 4. Gedächtnis-Nachlauf

[korr.] Diese Rechnung setzt einen Nachlauf der Feldgeometrie mit 1/λ_m = 170–260 Myr in der Milchstraße voraus; das trägt nicht: Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m = ω_orb/N*, N* ≈ 0,3–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Geschwindigkeitskern würde exzentrische Bahnen runden, alte Scheiben verlangen dafür λ_m ≪ 0,01·Ω. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche). Die Zahlen dieses Abschnitts gelten nur unter der Annahme.

Folgt die Feldgeometrie mit der Nachlaufzeit 1/λ_m (Kern λ_m e^{−λ_mΔt}), dann folgt der mitgeführte Anteil mit dem Faktor λ_m/(λ_m + p). Mit der komplexen Kippung Θ ∝ e^{pt} gilt [korr.]:

$$
L_o\,p=-i\Big[K_\text{ruh}+K_\text{mit}\,\frac{\lambda_m}{\lambda_m+p}\Big]
\quad\Longrightarrow\quad
p^2+(\lambda_m+ia)\,p+i\lambda_m(a+b)=0,\qquad a=\frac{K_\text{ruh}}{L_o},\; b=\frac{K_\text{mit}}{L_o}.
$$

Im p ist die Präzession, Re p das Anwachsen. Die Grenzfälle: λ_m → ∞ liefert p = −i(a+b), die Präzession aus 3.3; λ_m → 0 liefert p = −ia, die ruhende Geometrie.

| gekippter Teil | 1/λ_m [Myr] | Mode 1: Im p [km/s/kpc] | Re p [1/Gyr] | Mode 2: Im p | Re p [1/Gyr] |
|---|---|---|---|---|---|
| Außenteil R > 12 | 0 | −1,19 | 0 | – | – |
| | 100 | −0,23 | −0,50 | −19,45 | −9,50 |
| | 170 | −0,09 | −0,33 | −19,60 | −5,55 |
| | 260 | −0,04 | −0,23 | −19,65 | −3,62 |
| Ring 13–15 | 0 | −13,98 | 0 | – | – |
| | 100 | −0,84 | −4,77 | −28,39 | −5,23 |
| | 170 | −0,29 | −2,81 | −28,94 | −3,07 |
| | 260 | −0,12 | −1,84 | −29,11 | −2,01 |

Bei einer Nachlaufzeit ≥ 170 Myr zerfällt die Bewegung in zwei Moden:
- **Langsame Mode:** Scheibe und Feldgeometrie kippen gemeinsam. Sie präzedieren kaum und richten sich in 0,4–4 Gyr wieder aus.
- **Schnelle Mode:** Der Außenteil läuft rückläufig mit −19,6 bis −29 km/s/kpc um die ruhende Geometrie und klingt in 0,2–0,5 Gyr ab.

Der Nachlauf ist ein transversaler (dissipativer) Beitrag. Er ändert die Amplitude, erzeugt aber keine rechtläufige Präzession: Im p ≤ 0 in allen Fällen. In der Milchstraße gibt es diesen Nachlauf nicht messbar; es gelten die Zahlen ohne Nachlauf (3.3). [korr.]

---

## 5. Äußeres Feld der Lokalen Gruppe (EFE)

### 5.1 Symmetrie

Ein homogenes äußeres Feld g_e ê übt auf eine Masse in Newton kein Drehmoment aus. In der Feldgleichung geht es über ν(|∇Φ_N + g_{N,e}|) in die innere Feldgeometrie ein.

Das Drehmoment auf eine Scheibe mit Achse n ist ein axialer Vektor. Aus n (axial) und ê (polar) lässt sich in erster Ordnung in ê keiner bilden: (n·ê)n und n×ê sind polar. Der erste axiale Term ist (n·ê)(n×ê), also **quadratisch in g_e**. Daran ändert der Drehsinn der Scheibe nichts, denn die Massenverteilung ist spiegelsymmetrisch zur Mittelebene.

### 5.2 Rechnung

Die Milchstraße wird als kugelsymmetrisches Newton-Feld g_N(r) aus der Gittermasse angesetzt, dazu kommt das homogene Newton-Außenfeld g_{N,e} mit ν(g_{N,e})g_{N,e} = g_e. Der Rückkopplungsanteil der Quelle, geschrieben als gelesene Dichte ρ_ν (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt; kein Stoff) [korr.], ist

$$
4\pi G\rho_\nu=-\nabla\!\cdot\!\big[(\nu-1)\,\mathbf g\big]=(\nu-1)\,4\pi G\rho_b-\frac{\nu'(y)}{a_0}\,\mathbf g\cdot\nabla|\mathbf g|,\qquad \mathbf g=-g_N(r)\,\hat r+g_{N,e}\,\hat e .
$$

Projiziert wird auf P₂(r̂·ê) mit 400 Gauß-Legendre-Knoten. Der Beitrag zum Potential ist Φ₂ = C(r)r²P₂(r̂·ê) mit

$$
C(R)=-\frac{4\pi G}{5}\Big[R^{-5}\!\int_0^R\rho_2r^4dr+\int_R^\infty\frac{\rho_2}{r}dr\Big].
$$

Ein Ring in diesem Feld präzediert um ê mit |Ω_EFE| = (3/2)|C cos β|/Ω. Die Tabelle gibt den Höchstwert cos β = 1.

| g_e | g_{N,e}/a₀ | r_e (g_N = g_{N,e}) | C(14 kpc) [(km/s/kpc)²] | Ω_EFE bei 12 / 14 / 16 kpc [km/s/kpc] |
|---|---|---|---|---|
| 0,01 a₀ | 9,9·10⁻⁵ | 885 kpc | −4,4·10⁻³ | 3,6 / 4,3 / 5,0 ·10⁻⁴ |
| 0,02 a₀ | 3,9·10⁻⁴ | 445 kpc | −1,7·10⁻² | 1,4 / 1,7 / 2,0 ·10⁻³ |
| 0,03 a₀ | 8,7·10⁻⁴ | 298 kpc | −3,9·10⁻² | 3,2 / 3,8 / 4,4 ·10⁻³ |

Kontrolle: Für g_e → 0 verschwindet die gelesene Dichte ρ₂ (max. 4·10⁻⁴ M☉/kpc³). [korr.]

Die Anisotropie entsteht dort, wo das innere Newton-Feld auf das äußere abgefallen ist (r_e ≈ 300–900 kpc). Am Ort des Warps wirkt sie als schwache Quadrupol-Gezeit mit |C|/Ω² ≈ 7·10⁻⁵ (bei 0,02 a₀). **Das äußere Feld der Lokalen Gruppe trägt zur Präzession höchstens 4·10⁻³ km/s/kpc bei.** Eine spätere Rechnung mit verstärkten Feldern gibt höchstens 0,02 km/s/kpc; dieser Wert gilt. [korr.]

---

## 6. Große Magellansche Wolke

Masse (baryonisch) M_L = 3,2·10⁹ M☉ [G], Abstand d = 50 kpc [G], galaktische Breite b = −33° [G]. Damit ist |cos β| = 0,545 zwischen Scheibenachse und Richtung zur Wolke.

**Lineare Antwort der Feldgleichung.** Auf dem kugelsymmetrischen Hintergrund (ν aus |g_N| und g_{N,e} = 3,9·10⁻⁴ a₀) entkoppeln die Multipole. Für den Newton-Anteil δΦ_N = f(r)P_l gilt

$$
f=-\frac{GM_L r^l}{d^{\,l+1}}\ (r<d),\qquad f=-\frac{GM_L d^{\,l}}{r^{\,l+1}}\ (r>d).
$$

Der Rückkopplungsanteil u(r)P_l löst

$$
(r^2u')'-l(l+1)\,u=(r^2\alpha f')'-l(l+1)(\nu-1)f,\qquad \alpha=\nu-1+y\,\nu'(y),
$$

weil die Radialkomponente des Flusses mit α, die Tangentialkomponente mit ν − 1 verstärkt wird. Kontrolle: Bei konstantem ν − 1 = 1,5 und ν′ = 0 folgt u/f = 1,5016.

| R [kpc] | C_Newton [(km/s/kpc)²] | C_FTM | Faktor | Ω_prec Newton | Ω_prec FTM | Knotenrate bei Neigung θ_w = 0,1 |
|---|---|---|---|---|---|---|
| 12 | −0,110 | −0,719 | 6,53 | 6,9·10⁻³ | 3,3·10⁻² | 0,27 |
| 14 | −0,110 | −0,709 | 6,44 | 8,6·10⁻³ | 3,8·10⁻² | 0,32 (gilt: 0,36) [korr.] |
| 16 | −0,110 | −0,700 | 6,35 | 1,0·10⁻² | 4,4·10⁻² | 0,37 |

(Präzession und Knotenrate in km/s/kpc.) [korr.] Eine spätere Rechnung mit verstärkten Feldern gibt bei 14 kpc 0,36 km/s/kpc; dieser Wert gilt.

Die Wolke sitzt im schwachen Feld der Milchstraße (g_N/a₀ ≈ 0,03 bei 50 kpc). Die Feldgleichung verstärkt ihre Gezeit am Ort des Warps deshalb um den Faktor 6,4. Das ist die Rückkopplungssumme am Ort der Quelle, keine zusätzliche Masse.

Die Knotenrate ist die größte Geschwindigkeit, mit der das Drehmoment (3/2)|C|·cos β·sin β/Ω die Knotenlinie eines um θ_w = 0,1 geneigten Rings verschieben kann. Um so 10,9 km/s/kpc zu erreichen, müsste die Wolke baryonisch 1,1·10¹¹ M☉ haben, also das 1,6-Fache der Milchstraße.

Die Bahn der Wolke ist nahezu polar [G]. Ihre m = 1-Anregung der Scheibe hat deshalb eine feste Knotenrichtung und eine zeitlich veränderliche Stärke. Sie treibt die Neigung an (Ȧ/A ≲ 0,3 km/s/kpc bei θ_w = 0,1), aber keinen Umlauf der Knotenlinie.

---

## 7. Verwindung und Umlaufphase

### 7.1 Verwindung a_J

Aus dem Modell folgt J = ∫ρRv_φ d³x = 8,72·10¹³ M☉ kpc km/s. Daraus ergeben sich:
- a_J = J/(M_bc) = 4,35 pc,
- Δs = GM_b/(3c²) = 3,29·10¹³ m,
- Spinzahl j = a_J/(3Δs) = 1361. [korr.]

Die Knotenpräzession durch die Verwindung ist ω_LT = 2GJ/(c²R³):

| R [kpc] | ω_LT [km/s/kpc] | mit Faktor ν(R) |
|---|---|---|
| 12 | 4,8·10⁻⁶ | 9,2·10⁻⁶ |
| 14 | 3,0·10⁻⁶ | 6,4·10⁻⁶ |
| 16 | 2,0·10⁻⁶ | 4,7·10⁻⁶ |

Die Verwindung wirkt rechtläufig, liegt aber sechs Größenordnungen unter dem Messwert.

In der Siegel-Form des präzedierenden Torus geht die Knotenzahl nur als ζ² ein (Kopplungszahl 𝒦_S = ζ² + O(ζ⁴), invariant unter a_J → −a_J). Die Drehfaltung legt im Kern also keinen Drehsinn fest. Einen Drehsinn trägt nur der lineare Verwindungsterm, eben ω_LT.

### 7.2 Re τ_M und vertikale Phase

Im Bahntorus trägt Re τ_M = ω_rad/ω_azim = κ/Ω die Umlaufphase in der Ebene. Sie ist die Rosette mit der Apsidendrehung Ω − κ. Die Knotenbewegung gehört zum Vertikalzyklus. Deshalb wird analog Re τ_z := ν_z/Ω definiert. Die Knotenzahl ist dann

$$
\zeta=1-\mathrm{Re}\,\tau_z=\frac{\omega_p}{\Omega}.
$$

Für den gekippten Scheibenteil steht an Stelle von ν_z die kollektive Vertikalfrequenz aus 3.2, also ζ = −K/(L_oΩ).

| R | Im τ_M = r*/(2πR) (Punkt 1) [korr.] | Re τ_M = κ/Ω | Ω − κ [km/s/kpc] | ν_z/Ω (Stern) | ζ starr | ζ Ring | ζ Gaia |
|---|---|---|---|---|---|---|---|
| 12 | 0,117 | 1,351 | −6,3 | 3,48 | −0,056 | −0,89 | +0,61 |
| 14 | 0,100 | 1,339 | −5,2 | 3,47 | −0,078 | −0,92 | +0,72 |
| 16 | 0,088 | 1,337 | −4,4 | 3,48 | −0,092 | −0,94 | +0,83 |

[korr.] Im τ_M = r*/(2πR) ist Punkt 1 des Autors; das Modell schreibt |q_Δs| = Δs/r (Punkt 3). Die Wahl ist offen; in der stationären Regel sind beide auf allen Bahnen gleich (χ(s) = s·q/(1 − q)), sie betrifft nur die Schreibweise.

Im axialsymmetrischen Feld sind Radial- und Vertikalzyklus in linearer Ordnung getrennt. Die Nebendiagonale der Periodenmatrix ist dann null. Re τ_M bestimmt daher die Phase der Rosette, treibt aber die Knotenlinie nicht.

Der Gaia-Wert verlangt Re τ_z,eff = 1 − ζ = 0,28 bei 14 kpc. Dafür müsste die vertikale Rückstellung des gekippten Teils weit unter der Umlauffrequenz liegen. Die Feldgeometrie müsste also längs der Scheibenachse gestreckt sein. Die Feldgeometrie des FTM folgt den Baryonen und ist abgeplattet; deshalb ist ζ < 0 bei allen Radien.

---

## 8. Bilanz bei R = 12–16 kpc

| Anteil | Präzession [km/s/kpc] | Sinn |
|---|---|---|
| freie Präzession im nichtlinearen Feld (Feldform) [korr.], starrer Außenteil | −1,0 bis −1,2 | rückläufig |
| freie Präzession (Feldform) [korr.], einzelner Ring | −12,4 bis −16,1 | rückläufig |
| mit Gedächtnis-Nachlauf ≥ 170 Myr (nur unter dieser Annahme; in Galaxien kein messbarer Nachlauf) [korr.] | −0,3 bis 0 (langsame Mode), −19,6 bis −29 (schnelle, abklingend) | rückläufig bzw. keiner |
| äußeres Feld der Lokalen Gruppe | ≤ 0,02 (spätere Rechnung mit verstärkten Feldern) [korr.] | je nach Richtung |
| Große Magellansche Wolke (Präzession / größte Knotenrate) | ≤ 0,04 / ≤ 0,37 (bei 14 kpc 0,36 nach der späteren Rechnung) [korr.] | je nach Lage |
| Verwindung | 3·10⁻⁶ bis 1·10⁻⁵ | rechtläufig |
| **Gaia** | **+10,9 ± 3,2** | **rechtläufig** |

Mit dem Ansatz „baryonische Milchstraße + Feldgleichung + äußeres Feld + Magellansche Wolke + Verwindung“ liegt die Präzession bei R = 12–16 kpc zwischen −16 und +0,4 km/s/kpc (innere Dynamik in der Feldform gerechnet). [korr.] Der Abstand zum Gaia-Wert beträgt 10,5 bis 27 km/s/kpc, das sind 3,3 bis 8,4 systematische Standardabweichungen. Schräg einfallendes Gas dreht die Knotenlinie höchstens um 90°, bei 10,9 km/s/kpc nach 0,14 Mrd. Jahren. [korr.]

---

## 9. Offen

1. **Eigenmode des Warps.** Gerechnet sind nur die beiden Grenzfälle „starr“ und „einzelner Ring“. Ein echter Warp mit stetig wachsender Neigung θ(R) verlangt die volle quadratische Form V = ½Σ K_kl θ_kθ_l zwischen Ringen und ihre Eigenwerte, also das Analogon der Kippmoden nach Sparke & Casertano im Funktional V. In allen gerechneten Fällen ist K > 0. Ist die Form positiv definit, bleibt jede Mode rückläufig; zu bestimmen bleibt der Betrag, der zwischen den Grenzfällen −1 und −16 km/s/kpc liegen sollte (Feldform-Zahlen; für die örtliche Regel nicht gerechnet). [korr.]
2. **Lesart der Messgröße.** Der Gaia-Wert stammt aus einem Modell mit fester Amplitude und gemeinsamer Knotenlinie von Höhe und Geschwindigkeit. Eine zeitlich veränderliche Amplitude liefert einen zusätzlichen Term Ȧ sin(φ − φ_w) in v_z und verschiebt die angepasste Rate. Im FTM ist die von außen angeregte Änderung klein (Ȧ/A ≲ 0,3 km/s/kpc). Wie groß die interne Änderung aus der abklingenden schnellen Mode (Abschnitt 4, nur unter der Annahme eines Nachlaufs) [korr.] ist, hängt vom Anfangszustand ab und ist nicht gerechnet.
3. **Gedächtnis als Orientierung.** Hat sich die Orientierung der Baryonen in den letzten 1/λ_m geändert, etwa durch Einfall von Gas, dann ist die ruhende Feldgeometrie gegen die heutige Scheibe gekippt. Sie wirkt dann auf den Außenteil wie eine gekippte Rückstellachse. Präzession um diese Achse bleibt rückläufig; die sichtbare Knotenbewegung gegen die innere Scheibe hängt aber von der Vorgeschichte ab. Nicht gerechnet. [korr.] Das setzt einen Nachlauf in der Galaxie voraus; das trägt nicht: Das Gedächtnis wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn, und in Galaxien folgt die Regel einer Änderung der Baryonen ohne messbaren Nachlauf. Einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche).
4. **Dreidimensionale Rechnung** für die Magellansche Wolke (nichtlinear innerhalb ~11 kpc um die Wolke) und für die Richtung des äußeren Feldes.
5. **Vertikale Dynamikmasse am Sonnenort.** Mit dem hier benutzten Massenmodell ergibt die volle Feldform 94,6 M☉/pc², mit leichteren Scheiben 84,2; gemessen sind 71 ± 6 [G]. Die örtliche Regel, der die Bahnen folgen, gibt am Sonnenort 71 M☉/pc². [korr.] Die Präzessionsraten ändern sich zwischen beiden Modellen um weniger als 1 km/s/kpc.

---

## Literatur

- E. Poggio et al., Nature Astronomy 4, 590 (2020) — Warp-Präzession aus Gaia DR2 [G]
- P. J. McMillan, MNRAS 465, 76 (2017) — Massenmodell der Milchstraße [G]
- A.-C. Eilers, D. W. Hogg, H.-W. Rix, M. K. Ness, ApJ 871, 120 (2019) — Rotationskurve [G]
- K. Kuijken, G. Gilmore, ApJ 367, L9 (1991) — vertikale Dynamikmasse [G]
- M. Milgrom, MNRAS 403, 886 (2010) — quasilineare Formulierung (Feldform)
- L. S. Sparke, S. Casertano, MNRAS 234, 873 (1988) — Kippmoden von Scheiben
- J. Binney, ARA&A 30, 51 (1992) — Warps
- J. Lense, H. Thirring, Phys. Z. 19, 156 (1918)
