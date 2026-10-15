# Kovariante Fassung der Einzeltrichter-Hierarchie E0

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „Licht und Materie koppeln nur an die Metrik“, „Die Materie spürt Φ = Φ_N + 𝒢†Φ_s“ (Kurzfassung, 2.1, 3.1) → Bahnen und Gas folgen der örtlichen Regel a = ν(|g_N|)·g_N, die Metrik (Feldform) trägt nur deren wirbelfreien Teil; für Kugeln gleich. Ganze Körper ersetzt, Sprung ΔV/M „nicht gelungen“ (Kurzfassung, 2.2) → im Modell wird nur der Überschuss (ρ − ρ_J)₊ geglättet, der Sprung beim Einsatz verschwindet. Kern M(q̄) = (Δs/r*)·s̄^{3/2}·q̄/√(1−q̄) (2.5) → nicht der Kern des Modells (dort Δ(τ_M) mit |q_Δs| = Δs/r und χ). ρ_J an der Sonne „8,7 H-Atome/cm³“ → 11,9 m_p/cm³ (8,7 nur mit 1,36 m_p je H-Atom); Punktwert der Gezeit „ℓ = 25,7 AE = 0,0037 r*“ → 36,6 AE = 0,0052 r*. [G] an σ_Cassini, σ_EPM und „Cassini verlangt ℓ ≥ 1,30 r*“; Zahlen im Planetensystem als mit der unstetigen Glättung gerechnet markiert. Umbenannt: χ_i (charakteristische Funktion) → 𝟙_i; Multiplikator ζ → μ_L; geodätische Länge ρ_i → d_i (Kernargument ρ → d); η als Massenverschiebung (3.2) → δρ_ε.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: rechnungen/kovariante_e0.py
> Zahlenwert der Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen), durchgehend. Der Modellwert c·H₀/2π = 1,08·10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc) erscheint nur in einer Vergleichsspalte (Abschnitt 7.3).
> Werte mit **[G]** sind Literatur- und Beobachtungswerte aus dem Gedächtnis, nicht nachgeprüft.

---

## Kurzfassung

| Frage | Ergebnis |
|---|---|
| Wirkung | FTM-AS mit einem geglätteten Faltungspotential σ̄ als eigenem Feld, festgelegt durch eine Nebenbedingung D²σ̄ = 𝒢[σ] mit Multiplikator μ_L [korr.]. Das Argument von q ist q̄ = exp[−(ℓ₀²Ȳ)^{1/4}], Ȳ = h^{μν}∇_μσ̄∇_νσ̄. |
| Glättungsoperator | 𝒢 ersetzt im Quellterm von σ jeden kompakten Körper durch eine Plummer-Wolke gleicher Komar-Masse, gebildet im Ruhesystem des Körpers mit geodätischem Abstand. Körper = beschränkte Zusammenhangskomponente von {ρ_σ > ρ_J}. Im Modell wird nur der Überschuss (ρ_σ − ρ_J)₊ ersetzt, durch eine Plummer-Verteilung der Masse ∫(ρ_σ − ρ_J)dV. [korr.] |
| Glättungsradius | ℓ_i³ = GM_i / (√(3/2)·\|⟨Ē^TF⟩_i\|): spurfreie Hesse-Form von c²σ̄ auf dem Ruheschnitt, **gemittelt über die eigene Zone** (Fixpunkt). Das Zonenmittel entfernt das eigene Feld und die Gezeiten der eigenen Begleiter; übrig bleibt die Elterngezeit. |
| Newton-Grenzfall | exakt das Energiefunktional E0, mit \|T_i\| := √(3/2)·\|⟨∂∂Φ̄_N^TF⟩_{ℓ_i}\|. Körper spüren das über ihre Zone gemittelte Strukturpotential. Diskrete Probe: Wirkung = −V_E0 auf 6·10⁻¹⁴, Kraftgesetz auf 3·10⁻¹⁰. |
| Linsen | Ψ = Φ exakt; Licht spürt Φ_s, Sterne ⟨Φ_s⟩_ℓ. η − 1 ≈ 4·10⁻⁸. Gas und Bahnen folgen der örtlichen Regel a = ν(\|g_N\|)·g_N; die Metrik (Feldform) trägt nur deren wirbelfreien Teil, für Kugeln sind beide gleich. [korr.] |
| Erhaltung | Die Wirkung ist aus g, u, σ, σ̄, μ_L [korr.] gebaut, ohne Materie. Noether-Identität ∇^μ𝒮_μν = 0 auf der Lösung, mit Bianchi ∇_μT^μν = 0. |
| c_T | = c. Die Gezeit kommt aus dem Skalar σ̄, nicht aus der Weyl-Projektion der Metrik. Mit der vollen Riemann-Projektion würden die Zonen an Gravitationswellen koppeln: \|δc_T/c\| ≈ 3·10⁻³ bei 1 nHz. |
| Sonne | ℓ = 0,934 pc = 27,4 r*. Q₂ = 9,2·10⁻³¹ s⁻² (3·10⁻⁴ σ_Cassini [G]). Zusatz im Planetensystem: isotroper Hooke-Term κ = 2,0·10⁻³⁰ s⁻², an Saturn δa/g_N = 4,5·10⁻¹⁴; alle Periheldrehungen ≤ 1,1·10⁻³ σ_EPM [G]. γ − 1 ≤ 9·10⁻⁴ bis 30 kAE. Die Zahlen im Planetensystem sind noch mit der unstetigen Glättung (ganze Körper ersetzt) gerechnet. [korr.] |
| Kovariante Korrekturen | Vorzugssystem (v/c)² = 1,5·10⁻⁶ relativ, also ~10⁻³⁶ s⁻² absolut; Retardierung ebenso; Komar- statt σ-Fluss 2·10⁻⁶ (Sonne), 24 % (Neutronenstern). |
| Nicht gelungen | Masse und Schwerpunkt dynamischer Starkfeldsysteme; kosmologischer Grenzfall der Körperdefinition. Der Sprung beim Einsetzen der Glättung (ΔV/M = 5·10² m²/s² für eine Wolke von 1 M☉ an der Schwelle) entsteht nur, weil hier ganze Körper ersetzt werden; im Modell wird nur der Überschuss (ρ − ρ_J)₊ der Quelldichte geglättet, und der Sprung beim Einsatz verschwindet. [korr.] |

---

## 1. Was kovariant zu machen ist

In FTM-AS bildet der Strukturterm den Nome aus dem ganzen Gradienten des Faltungspotentials,

$$
q=\exp\!\big[-(\ell_0^2\mathcal Y)^{1/4}\big],\qquad \mathcal Y=h^{\mu\nu}\nabla_\mu\sigma\nabla_\nu\sigma,\qquad \ell_0=c^2/a_0 ,
$$

und im schwachen Feld gilt c²σ = Φ_N. E0 verlangt stattdessen das über die Gezeitengrenzen geglättete Feld:

$$
V=V_N-\frac{a_0^2}{8\pi G}\int Q_s\!\Big(\frac{|\nabla\bar\Phi_N|^2}{a_0^2}\Big)d^3x,\qquad
\nabla^2\bar\Phi_N=4\pi G\Big[\sum_i m_iK_{\ell_i}(\mathbf x-\mathbf x_i)+\rho_g\Big],\qquad
\ell_i=\Big(\frac{Gm_i}{|\mathsf T_i|}\Big)^{1/3},
$$

mit Q_s(y) = ∫₀^y (ν(√y′) − 1) dy′ und dem Plummer-Kern K_ℓ(x) = (3/4πℓ³)(1 + |x|²/ℓ²)^{−5/2}.

Vier Bausteine brauchen eine kovariante Definition:

1. **Körper:** Welche Teile der Quelle sind kompakte Körper, welche diffus?
2. **Masse und Ort:** m_i und x_i ohne Bezug auf Newtonsche Größen.
3. **Glättungskern:** Abstand und Dichte der Plummer-Wolke.
4. **Gezeit:** Welcher Tensor, welche Invariante, an welchem Ort? Dazu gehört, dass die Gezeit des Körpers selbst und die seiner eigenen Begleiter nicht mitzählen dürfen.

Alle vier werden aus g, u und σ gebaut. Die Materie tritt im Strukturteil nicht auf.

---

## 2. Die Wirkung FTM-AS-E0

### 2.1 Felder

| Feld | Rolle |
|---|---|
| g_μν | Metrik; Licht koppelt nur an sie. Die Materie koppelt an sie, ihre Bahnen folgen aber der örtlichen Regel a = ν(\|g_N\|)·g_N; die Metrik trägt nur deren wirbelfreien Teil (für Kugeln gleich). [korr.] |
| u^μ, u·u = −1 | Zeitfeld des Mediums, h^{μν} = g^{μν} + u^μu^ν, D_μ = h_μ^ν∇_ν |
| σ | Faltungspotential, B_μ = a_μ − ∇_μσ wie in FTM-AS |
| σ̄ | geglättetes Faltungspotential (Hilfsfeld, keine eigene Dynamik) |
| μ_L [korr.] | Multiplikator der Glättungs-Nebenbedingung; auf der Lösung ist c²μ_L das Strukturpotential Φ_s |

### 2.2 Körper, Masse, Weltlinie

Geometrische Dichte auf dem Ruheschnitt:

$$
\rho_\sigma=\frac{c^2}{4\pi G}\,e^{-\sigma}D_\mu\big(e^{\sigma}D^\mu\sigma\big)\;\xrightarrow{\text{schwaches Feld}}\;\rho .
$$

**Körper.** Ein Körper V_i ist eine beschränkte Zusammenhangskomponente der Menge

$$
\{x:\ \rho_\sigma(x)>\rho_J(x)\},\qquad \rho_J=\frac{3\sqrt{3/2}\,|\bar E^\text{TF}|}{4\pi G}.
$$

ρ_J ist die Gezeitendichte: Eine Kugel mittlerer Dichte ρ_J hat genau den Radius ℓ aus 2.4. Sterne, Planeten und Kometen liegen um viele Größenordnungen darüber. Das diffuse interstellare Gas liegt darunter. An der Sonne ist ρ_J = 0,293 M☉/pc³, das entspricht 11,9 m_p/cm³ (8,7 Wasserstoffatome je cm³ nur mit 1,36 m_p je H-Atom, also mit Helium). [korr.] Kugelsternhaufen und Zwerggalaxien zerfallen in ihre Sterne, weil zwischen den Sternen ρ_σ = 0 ist; sie bleiben aufgelöst wie in E0.

[korr.] Hier wird jeweils der ganze Körper ersetzt; das gibt beim Einsetzen der Glättung einen Sprung. Im Modell wird nur der Überschuss (ρ_σ − ρ_J)₊ der Quelldichte geglättet, durch eine Plummer-Verteilung gleicher Masse M^c = ∫(ρ_σ − ρ_J)dV (Komar-Masse des Überschusses); der Sprung beim Einsatz verschwindet.

**Masse als Komar-Fluss.** Über eine geschlossene Fläche um V_i im Ruheschnitt:

$$
M_i=\frac{c^2}{4\pi G}\oint e^{\sigma}\,n^\mu D_\mu\sigma\,dA .
$$

Im Innern kompakter Körper ist q = 0 und σ = ln N, also e^σDσ = DN. Das ist die Komar-Masse. Im statischen Vakuum hängt sie nicht von der Fläche ab. Der Fluss ohne e^σ gäbe M/N(R), das sind 2,1·10⁻⁶ zu viel an der Sonnenoberfläche und 24 % bei einem Neutronenstern (1,4 M☉, 12 km).

**Weltlinie.** z_i ist der Schwerpunkt der Komponente im Ruhesystem ihres eigenen Impulses (Massendipol null), U_i die Vierergeschwindigkeit dieser Linie.

### 2.3 Glättungskern

Für einen Punkt x sei d_i(x) [korr.] die Länge der Geodäte von z_i nach x, die senkrecht auf U_i steht. Die Plummer-Wolke des Körpers ist der Strom

$$
\bar J_i^\mu=M_i\,\mathcal K_{\ell_i}\big(d_i(x)\big)\,U_i^\mu,\qquad \mathcal K_\ell(d)=\frac{3}{4\pi\ell^3}\Big(1+\frac{d^2}{\ell^2}\Big)^{-5/2}.
$$

Die Glättung ist damit im Ruhesystem des Körpers erklärt, ohne Bezug auf u. Das Medium sieht die Dichte −u_μJ̄_i^μ = γ_i M_i𝒦, mit γ_i = −u·U_i.

Der geglättete Quellterm:

$$
\boxed{\;\mathcal G[\sigma]=D^2\sigma-\sum_i\Big[\mathbb{1}_i\,D^2\sigma-\frac{4\pi G}{c^2}\,\gamma_iM_i\,\mathcal K_{\ell_i}(d_i)\Big]\;}
$$

Hier ist 𝟙_i die charakteristische Funktion von V_i. [korr.] Jeder Körper wird aus der Quelle herausgenommen und durch seine Plummer-Wolke gleicher Masse ersetzt. Das diffuse Gas bleibt unverändert.

### 2.4 Gezeit und Glättungsradius

**Gezeitentensor.** Die spurfreie Hesse-Form des geglätteten Faltungspotentials auf dem Ruheschnitt:

$$
\bar E_{\mu\nu}=c^2\Big(D_\mu D_\nu\bar\sigma-\tfrac13 h_{\mu\nu}D^2\bar\sigma\Big).
$$

In statischen Lagen ist das der elektrische Teil der Krümmung des Skalarsektors, im schwachen Feld die spurfreie Newtonsche Gezeit ∂_i∂_jΦ̄_N^TF.

**Zonenmittel.** Die Gezeit wird über die eigene Zone gemittelt und auf z_i parallel verschoben:

$$
\langle\bar E\rangle_i=\int_\Sigma \mathcal K_{\ell_i}(d_i)\,P_{z_i\leftarrow x}\,\bar E(x)\,dV .
$$

**Glättungsradius als Fixpunkt:**

$$
\boxed{\;\ell_i^3\,\sqrt{3/2}\,\big|\langle\bar E\rangle_i\big|=GM_i\;},\qquad |E|=\sqrt{E_{\mu\nu}E^{\mu\nu}} .
$$

**Normierung.** Für eine Punktmasse als Eltern auf der Kreisbahn ist Ē^TF = (GM/R³)·diag(−2, 1, 1), also √(3/2)|Ē^TF| = 3GM/R³. Das ist genau der Jacobi-Wert 4A(A−B) = 3Ω² der Kepler-Bahn. Der Faktor √(3/2) macht die Invariante zur Jacobi-Gezeit, ohne den Bahnrahmen zu benutzen.

**Warum das Zonenmittel.** Für eine Quelle im Abstand d vom Schwerpunkt gilt für jede Kugelschale um z_i:
- r > d: Das Schalenmittel der Hesse-Form verschwindet. Jede Ableitung einer irregulären Kugelfunktion vom Grad l ist eine irreguläre Kugelfunktion vom Grad l + 1, das Mittel über die Kugel ist null.
- r < d: Das Feld ist harmonisch, das Mittel ist der Zentralwert.

Daraus folgt exakt

$$
\langle\partial\partial\Phi^\text{TF}\rangle_\ell=\partial\partial\Phi^\text{TF}(z_i)\cdot F(d/\ell),\qquad F(x)=\frac{x^3}{(1+x^2)^{3/2}} .
$$

Numerische Probe (Quadratur über Schalen): 1,757443 gegen 1,757479 bei d/ℓ = 0,3; 0,707099 gegen 0,707107 bei d/ℓ = 1; 0,063245 gegen 0,063246 bei d/ℓ = 3.

Folgen:
- Das eigene, kugelsymmetrische Feld trägt nichts bei (d = 0).
- Ein Begleiter im Abstand d ≪ ℓ trägt ≈ 3GM_p/ℓ³ bei, unabhängig von d, also nur im Verhältnis seiner Masse.
- Äußere Quellen (d ≫ ℓ) zählen voll.

Das ist die Hierarchie „Elterngezeit“ ohne Baum: Welche Quelle Eltern und welche Kind ist, folgt aus dem Mittel über die Zone.

Mit dem Punktwert statt des Zonenmittels wäre an der Sonne Jupiters Gezeit (8,1·10⁻¹⁹ s⁻²) um den Faktor 10¹¹ größer als die galaktische, und ℓ läge nach ℓ³ = GM☉/T bei 36,6 AE = 0,0052 r*. [korr.] Cassini verlangt ℓ ≥ 1,30 r* [G].

### 2.5 Die Wirkung

$$
\boxed{\;
S=\frac{c^3}{16\pi G}\int d^4x\sqrt{-g}\Big[R+2h^{\mu\nu}B_\mu B_\nu+\frac{2}{\ell_0^2}Q_s\big(\ell_0^2\bar{\mathcal Y}\big)-2\Lambda_s-\lambda_u(u_\mu u^\mu+1)+4\mu_L\big(D_\mu D^\mu\bar\sigma-\mathcal G[\sigma]\big)\Big]+S_m[g,\psi]\;}
$$

$$
\bar{\mathcal Y}=h^{\mu\nu}\nabla_\mu\bar\sigma\nabla_\nu\bar\sigma,\qquad \bar q=\exp\!\big[-(\ell_0^2\bar{\mathcal Y})^{1/4}\big],\qquad \bar\nu=\frac1{1-\bar q},\qquad \Lambda_s=\frac{4\pi^4}{15\,\ell_0^2}.
$$

**Aufteilung der Kinetikfunktion.** Die Kinetikfunktion von FTM-AS, Q̃ = −4[s³Li₁(q) + 3s²Li₂(q) + 6sLi₃(q) + 6Li₄(q)], zerfällt exakt in Q̃ = Q_s − 4π⁴/15. Q_s ist der Strukturteil mit Q_s(0) = 0, der Rest ist die globale Konstante −2Λ_s. Die Glättung wirkt nur auf das Argument von Q_s; die Konstante hängt nicht von σ̄ und nicht von den ℓ_i ab. Eine Vakuumenergie je Zone entsteht deshalb nicht. Ob Λ_s physikalisch ist, bleibt eine eigene Frage.

**Kern der Bahngleichung.** [korr.] Die Form M(q̄) = (Δs/r*)·s̄^{3/2}·q̄/√(1−q̄) mit s̄ = −ln q̄ stammt aus einer älteren Rechnung und ist nicht der Kern des Modells. Dort steht Δ(τ_M) mit dem Kern-Nome |q_Δs| = Δs/r und dem Schalter χ = s/(e^s − 1); das geglättete Feld geht über s̄ in χ ein.

### 2.6 Feldgleichungen

| Variation | Gleichung |
|---|---|
| δμ_L [korr.] | D²σ̄ = 𝒢[σ] (Glättung) |
| δσ̄ | D_μ[(ν̄ − 1)D^μσ̄] = D²μ_L + (Terme aus ∂ℓ_i/∂σ̄ auf den Weltlinien) [korr.] |
| δσ | Gleichung von FTM-AS, der Strukturteil geht über den adjungierten Operator: ∇_μ[h^{μν}(∇_νσ − a_ν)] + 𝒢†μ_L = 0 [korr.] |
| δu, δg | wie FTM-AS, dazu die Beiträge von μ_L(D²σ̄ − 𝒢) über h und über die Weltliniendaten [korr.] |

Der adjungierte Operator wirkt auf eine Funktion f wie

$$
\mathcal G^\dagger f=\Big(1-\sum_i\mathbb{1}_i\Big)f+\sum_i\mathbb{1}_i\,\langle f\rangle_{\mathcal K_{\ell_i}},
$$

also: außerhalb der Körper unverändert, im Körper i ersetzt durch das Plummer-Mittel über seine Zone.

---

## 3. Newton-Grenzfall

### 3.1 Herleitung

Mit ds² = −(1 + 2Φ/c²)c²dt² + (1 − 2Ψ/c²)dx², statischem u, c²σ = Φ_N, c²σ̄ = Φ̄_N und c²μ_L = φ [korr.] wird die Wirkung zu

$$
L_\text{NR}=\frac1{8\pi G}\Big[(\nabla\Psi)^2-2\nabla\Phi\cdot\nabla\Psi+(\nabla\Phi-\nabla\Phi_N)^2+a_0^2Q_s\Big(\frac{|\nabla\bar\Phi_N|^2}{a_0^2}\Big)+2\varphi\big(\nabla^2\bar\Phi_N-\mathcal G[\nabla^2\Phi_N]\big)\Big]-\rho\Phi .
$$

Der Faktor folgt aus (c⁴/16πG)·4·(φ/c²)·(∇²Φ̄_N − 𝒢)/c² = (1/8πG)·2φ(…).

| Variation | Ergebnis |
|---|---|
| δΨ | Ψ = Φ |
| δΦ | ∇²Φ_N = 4πGρ |
| δφ | ∇²Φ̄_N = 𝒢[4πGρ] = 4πG[Σ_i m_iK_{ℓ_i} + ρ_g], die Glättungsgleichung von E0 |
| δΦ̄_N | ∇²φ = ∇·[(ν̄ − 1)∇Φ̄_N], also φ = Φ_s, das Strukturpotential von E0 |
| δΦ_N | Φ = Φ_N + 𝒢†Φ_s |

Kein Feld hat im Newton-Grenzfall eine Zeitableitung. Sie lassen sich daher eliminieren. Auf der Lösung heben sich die Terme −(1/4πG)∇Φ·∇Φ_N und −ρΦ gegenseitig auf, die Nebenbedingung verschwindet, und es bleibt

$$
-L_\text{NR}\big|_\text{Lösung}=-\frac1{8\pi G}\int|\nabla\Phi_N|^2-\frac{a_0^2}{8\pi G}\int Q_s\Big(\frac{|\nabla\bar\Phi_N|^2}{a_0^2}\Big)=V_\text{E0}.
$$

Die Bewegung der Materie folgt aus L = T − V_E0. Das ist E0, einschließlich der Kraft aus der Ortsabhängigkeit der ℓ_i, mit der Festlegung

$$
|\mathsf T_i|:=\sqrt{3/2}\;\big|\langle\partial_i\partial_j\bar\Phi_N^\text{TF}\rangle_{\ell_i}\big| .
$$

In der Feldform spürt die Materie Φ = Φ_N + 𝒢†Φ_s, Licht Φ_s am Ort. [korr.] Für die Bahnen gilt das nur bei Kugeln: Bahnen und Gas folgen der örtlichen Regel a = ν(|g_N|)·g_N, die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Ein Körper spürt ⟨Φ_s⟩_{K_{ℓ_i}}, konstant über sein Gebiet: Er bewegt sich als Einheit, ohne innere Strukturkraft. Das ist die Regel „Körper i spürt Φ_N + K_{ℓ_i}∗Φ_s“ von E0. Hier folgt sie aus der Adjungierten der Glättung.

### 3.2 Numerische Probe (3D-Gitter)

Das Gitter hat 14³ Knoten, G = a₀ = 1. Ein Körper der Masse 20 liegt auf 7 Knoten, er wird über ℓ = 2,5 geglättet. Dazu kommt Gas der Masse 8. Es gilt |∇Φ̄_N|/a₀ = 0,03 … 1,6, der Übergangsbereich ist also voll erfasst. Die Felder werden aus den diskreten Gleichungen gelöst, dann wird die diskrete Wirkung geprüft:

| Probe | Ergebnis |
|---|---|
| Richtungsableitungen δS/δ(Ψ, Φ, Φ_N, Φ̄_N, φ) in der Lösung | ≤ 3·10⁻¹¹, Vergleichsgröße 0,02–0,23 |
| Wirkung in der Lösung gegen −V_E0 | 277,9674557766 gegen 277,9674557766, Differenz 6·10⁻¹⁴ |
| Strukturpotential im Körpergebiet | konstant −2,4856964828 = K·φ auf allen 7 Knoten |
| Kraftgesetz: dV_E0/dε bei Verschiebung von Masse (neu gelöst) | −13,0475168348 |
| Σ Φ·δρ_ε [korr.] mit Φ = Φ_N + 𝒢†φ, δρ_ε = dρ/dε | −13,0475168344 |
| dasselbe ohne Zonenmittel im Körper, Σ(Φ_N + φ)·δρ_ε [korr.] | −13,2309343419, also 1,4 % daneben |

Das Kraftgesetz von E0 gilt also genau dann, wenn der Körper das Zonenmittel spürt.
