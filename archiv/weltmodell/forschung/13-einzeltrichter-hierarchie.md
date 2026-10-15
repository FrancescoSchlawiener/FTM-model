# Einzeltrichter und Hierarchie – Trichter in Trichtern

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Glättungslänge der Sonne ℓ = r_J = 1,372 pc = 40 r* (mit T = 4A(A − B)) → ℓ³ = GM/T mit T = √(3/2)·|E^TF| (voller spurfreier Gezeitentensor): aus der Norm |T|_F = 6,00·10⁻³⁰ s⁻² dieser Notiz 0,91 pc, im Modell 0,93 pc = 27 r*; 1,37 pc ist der klassische Jacobi-Radius; Quadrupol aus der Glättung 2,9·10⁻³¹ (für 1,37 pc) → 9,3·10⁻³¹ s⁻² (für 0,93 pc); kovariante Glättungslänge aus dem elektrischen Teil der Riemann-Krümmung → aus der spurfreien Hesse-Form von c²σ̄ mit √(3/2), nicht aus der Weyl-Projektion; „η = 1 gilt unverändert“ → η = 1 für Kugeln und außerhalb von Scheiben, in Scheiben innen η > 1; „Körper spüren die QUMOND-Feldform“ → die Bahnen folgen der örtlichen Regel, die Feldform trägt nur deren wirbelfreien Teil (für Kugeln und Zwerge gleich, die Zahlen bleiben); Glättungsgewichte aus g(R)/a₀ → Überschuss (ρ − ρ_J)₊ über die Gezeitendichte, a₀ kommt nicht vor; Kugelsternhaufen als Einzeltrichter → als Sternsysteme kollektiv, kompakte mit ν = 1,008; „lokale Strukturdichte“, „Strukturpotential“ → gelesene Dichte der Regel, Zusatzpotential der Regel. Umbenennungen (ohne Marke): χ_i (charakteristische Funktion) → 𝟙_i; η (Slip) → γ_s; Δϖ̇ (Periheldrehung) → Δω̇_p; s (Abstand der Doppelsterne) → d. Rechnet mit a₀ = 1,2·10⁻¹⁰, in einzelnen Spalten mit cH₀/2π = 1,08·10⁻¹⁰ (H₀ = 70).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: rechnungen/hierarchie.py
> Literatur- und Beobachtungswerte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.
> Zahlenwert der Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Wo der Modellwert c·H₀/2π = 1,0824·10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc) benutzt wird, steht es dabei.

---

## Kurzfassung

| Frage | Ergebnis |
|---|---|
| Lässt sich „Trichter in Trichtern“ konsistent formulieren? | Ja, als **Energiefunktional** mit geglättetem Elternfeld (1.3): Jeder kompakte Körper wird über seine Gezeitengrenze geglättet, verstärkt wird nur das kollektive Feld. Energie, Impuls, Drehimpuls und Galilei-Invarianz gelten exakt. |
| Was kostet das? | **Lokalität.** Ein homogenes äußeres Feld lässt sich im Vakuum durch keine lokale Mittelung vom eigenen Feld trennen (Mittelwertsatz, 1.5). Dazu kommt eine kleine Zusatzkraft aus der Umgebungsabhängigkeit der Zonen: 4·10⁻⁹ g_obs. |
| Kovariant? | Ja, als FTM-AS mit nichtlokalem Glättungsoperator auf den Ruheschnitten des Mediums (2). |
| Cassini | Die Glättung liefert Q₂ = 9,3·10⁻³¹ s⁻² (bei ℓ = 0,93 pc), dazu die galaktische Gezeit 4,3·10⁻³⁰ s⁻². Cassini verlangt eine Übergabe außerhalb von **ℓ ≥ 1,30 r\*** (2σ); die Glättungslänge liegt bei 27 r\*. [korr.] |
| Weite Doppelsterne (1,5 M☉) | **γ = 1,000** bei 5–30 kAE. |
| Oortsche Wolke | Newtonsch bis 10⁵ AE; die Oort-Spitze wird mit 1/a_app = +4,25·10⁻⁵ AE⁻¹ wiedergegeben. |
| Zwerggalaxien | Aufgelöst über ihre Sterne; es wirkt das volle kollektive Feld mit äußerem Feld (EFE). |
| Folgerung | Kompakte Einzeltrichter (g(R) > a₀) bleiben bis zu ihrer Gezeitengrenze im Elternfeld Newtonsch; verstärkt wird nur das kollektive, über diese Zonen geglättete Feld. Das ist verträglich mit Cassini, EPM, Oort-Spitze und Banik et al. **Schärfster Test: weite Doppelsterne bei 5–30 kAE**, Vorhersage γ = 1,000. |

---

## 0. Ausgangsgrößen

**Sonnenumgebung [G]:** v₀ = 233 km/s, R₀ = 8,178 kpc. Daraus folgt g_obs = v₀²/R₀ = 2,151·10⁻¹⁰ m/s², mit v₀ = 229 km/s sind es 2,078·10⁻¹⁰. Das entspricht 1,8–2,3·10⁻¹⁰ m/s² ≈ 1,5–2 a₀.

**Wichtige Präzisierung:** Im QUMOND-Argument von ν steht das **Newtonsche** Feld, nicht das beobachtete. Die Umkehr der RAR, g_obs = ν(g_N/a₀)·g_N, ergibt:

| a₀ | g_obs/a₀ | g_Ne | ν_e | K = dln ν/dln g |
|---|---|---|---|---|
| cH₀/2π = 1,0824·10⁻¹⁰ | 1,988 | 1,484·10⁻¹⁰ = 1,371 a₀ | 1,4494 | −0,2631 |
| 1,2·10⁻¹⁰ | 1,793 | 1,429·10⁻¹⁰ = 1,191 a₀ | 1,5056 | −0,2759 |

**Galaktische Gezeiten [G: A = 15,3, B = −11,9 km/s/kpc, ρ₀ = 0,10 M☉/pc³]:**
- Komponenten: Φ_RR = −9,71·10⁻³¹, Φ_φφ = 7,77·10⁻³¹, Φ_zz = 5,87·10⁻³⁰ s⁻².
- 4A(A−B) = 1,748·10⁻³⁰ s⁻², Norm |T|_F = 6,00·10⁻³⁰ s⁻².
- Das ist das **effektive** Gezeitenfeld.

**Hierarchie.** Jeder Trichter k hat drei Radien:
- r\*_k = √(GM_k/a₀),
- r_eq,k = √(GM_k/g_Eltern): hier ist die eigene Newtonsche Neigung gleich der des Elternfeldes,
- r_J,k: die Gezeitengrenze.

| Trichter (in Eltern) | r\* | r_eq | r_J | r_J/r\* | g_Eltern/a₀ |
|---|---|---|---|---|---|
| Erde (Sonne) | 12,18 AE | 0,0017 AE | 0,01 AE | 0,0008 | 4,9·10⁷ |
| Jupiter (Sonne) | 217,2 AE | 0,161 AE | 0,355 AE | 0,0016 | 1,8·10⁶ |
| Sonne (Milchstraße) | 7030 AE | 6442 AE | 2,83·10⁵ AE = 1,372 pc (klassischer Jacobi-Radius; Glättungslänge ℓ = 0,93 pc = 27 r\*) [korr.] | 40,3 | 1,19 |
| Doppelstern 1,5 M☉ (Milchstraße) | 8610 AE | 7890 AE | 3,24·10⁵ AE | 37,6 | 1,19 |
| Milchstraße (Lokale Gruppe) [G] | 8,35 kpc | 493 kpc | 399 kpc | 47,7 | 2,9·10⁻⁴ |

Planeten erreichen ihr eigenes r\* nie, denn die Gezeitengrenze liegt 600- bis 1000-mal weiter innen. Sie sind reine Newton-Trichter. Die eigentliche Frage betrifft den Bereich zwischen r\* und r_J der Sonne, also 7 000 bis 283 000 AE, und die entsprechende Zone jedes Sterns. Für Sterne in Scheibengalaxien gilt r_eq ≈ r\* (hier 0,916 r\*), weil die Scheiben selbst an der Schwelle liegen (Freeman-Grenze).

---

## 1. Fraktale Feldzerlegung

### 1.1 Die drei Rollen des äußeren Feldes

Ein Trichter k hat Schwerpunkt X_k, Masse M_k und Zone ℓ_k. Das Feld der Umgebung (Eltern, Geschwister) wird an X_k entwickelt:

$$
\mathbf g_\text{ext}(\mathbf x)=\underbrace{\mathbf A_k}_{\text{homogen}}+\underbrace{\mathsf T_k\,(\mathbf x-\mathbf X_k)}_{\text{Gezeiten}}+O(|\mathbf x-\mathbf X_k|^2),\qquad
\ell_k=\ell\big(M_k,\mathsf T_k\big)\ \ \text{(Parameter)} .
$$

1. **Homogen, freier Fall.** A_k beschleunigt alle Glieder des Trichters gleich. In Relativkoordinaten x′ = x − X_k fällt dieser Anteil exakt heraus (Äquivalenzprinzip). A_k steht deshalb **nicht** im Argument von ν.
2. **Gezeiten, additiv.** Das Gezeitenfeld T_k wirkt linear auf die Glieder, und zwar mit dem **effektiven** Elternfeld (Newton plus Elternstruktur). Das ist die gewöhnliche galaktische Gezeit, wie sie die Oortschen Konstanten messen.
3. **Parameter.** Das Elternfeld legt die Ausdehnung der Zone fest: ℓ_k ist die Glättungslänge ℓ_k = (GM_k/T_k)^{1/3} mit T_k = √(3/2)·|E^TF| aus dem vollen spurfreien Gezeitentensor, nicht mit 4A(A−B). Schon mit der Norm |T|_F = 6,00·10⁻³⁰ s⁻² aus Abschnitt 0 folgt für die Sonne (GM☉/|T|_F)^{1/3} = 0,91 pc, also der Modellwert ℓ ≈ 0,93 pc = 27 r\*; 1,37 pc aus 4A(A−B) ist der klassische Jacobi-Radius. [korr.] Das ist die Rolle eines langsamen Parameters im Zustandsraum, mit schnellen Variablen z_f und langsamen Parametern p:

$$
\dot z_\text f=F(z_\text f;\,p),\qquad \dot p=\epsilon\,G(z_\text f,p),\qquad p=(\mathbf A_k,\mathsf T_k,\ell_k).
$$

Die schnellen Variablen z_f sind die Bahnen im Trichter. Das Elternfeld ist kein Teil des eigenen Wirbels, sondern stellt dessen Rahmen.

### 1.2 Bausteine

Die Rechnung verwendet folgende Bausteine:
- **Plummer-Kern.** K_ℓ(x) = (3/4πℓ³)(1 + |x|²/ℓ²)^{−5/2}. Das Feld einer über ℓ geglätteten Masse ist GM r/(r² + ℓ²)^{3/2}.
- **Strukturenergie-Dichte.**
  $$Q_s(y)=\int_0^y\big(\nu(\sqrt{y'})-1\big)\,dy'=4\int_0^{s}\frac{t^3\,dt}{e^t-1},\qquad s=\sqrt{|\mathbf g|/a_0}=-\ln q .$$
  Die Normierung ist Q_s(0) = 0. Das ist die Strukturfunktion Q̃ des Modells ohne deren Konstante −4π⁴/15. Diese Normierung ist **erzwungen**. Mit der Konstanten stünde in jeder Zone eine Vakuumenergie ∝ ℓ³, die für das Plummer-Fenster sogar divergiert. Der Vakuumterm darf also höchstens global auftreten, nicht je Trichter.
- **Übergabe.** W_ℓ(r) = 1 − r³/(r² + ℓ²)^{3/2} ist der Anteil der Masse, den das Elternfeld bei r noch **nicht** sieht. Eigenes Feld und Elternfeld teilen sich die Masse wie in einer Ewald-Zerlegung: Der kurzreichweitige Teil g − g_K bleibt beim Körper und ist Newtonsch, das kollektive Feld übernimmt den geglätteten Teil g_K.

### 1.3 Das Energiefunktional (konsistente Theorie, E0)

Alle Kräfte kommen aus **einer** translationsinvarianten Energie; dadurch gilt actio = reactio exakt. Punktmassen m_i bei x_i, dazu diffuses Gas ρ_g. Die Newtonsche Energie ist exakt:

$$
V_N=-\frac1{8\pi G}\int|\nabla\Phi_N|^2d^3x,\qquad \nabla^2\Phi_N=4\pi G\rho .
$$

**E0: Verstärkt wird nur das kollektive Feld.**

$$
\boxed{\;V=V_N-\frac{a_0^2}{8\pi G}\int Q_s\!\Big(\frac{|\nabla\bar\Phi_N|^2}{a_0^2}\Big)d^3x\;},\qquad
\nabla^2\bar\Phi_N=4\pi G\Big[\sum_i m_i\,K_{\ell_i}(\mathbf x-\mathbf x_i)+\rho_g\Big],\qquad \ell_i=\Big(\frac{Gm_i}{|\mathsf T(\mathbf x_i)|}\Big)^{1/3}.
$$

Jeder kompakte Körper (Planet, Stern, Komet) wird über seine Gezeitengrenze im kollektiven Feld geglättet. Das kollektive Feld ist QUMOND auf dieser geglätteten Dichte. Die Bahnen folgen dabei der örtlichen Regel a = ν(|g_N|)·g_N; QUMOND ist die Feldform und trägt nur deren wirbelfreien Teil (Feldgleichung unten). [korr.] Es gibt keinen Baum und keine eigenen Strukturterme. Die Hierarchie entsteht allein aus den ℓ_i:
- Sterne einer Zwerggalaxie haben ℓ_i ≪ Größe der Galaxie. Die Galaxie ist aufgelöst und kollektiv.
- Ein Doppelstern bei 20 kAE hat ℓ ≈ 3·10⁵ AE ≫ d. Er ist nicht aufgelöst und innen Newtonsch.

**Feldgleichung, schwaches Feld:**

$$
\nabla^2\Phi=\nabla^2\Phi_N+\nabla\!\cdot\!\big[(\nu(|\nabla\bar\Phi_N|/a_0)-1)\nabla\bar\Phi_N\big],\qquad \text{Körper } i \text{ spürt } \Phi_N+K_{\ell_i}*\Phi_s .
$$

[korr.] Das ist die Feldform, also die Metrik. Die Körper selbst folgen der örtlichen Regel a = ν(|ḡ_N|)·ḡ_N im geglätteten Feld; die Feldform trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit f = +0,03 bis +0,13 ein, ohne IC 2574 mit 0,00 ± 0,03; die Feldform (f = 1) liegt außerhalb 95 %. Für Kugeln und kugelnahe Systeme (Sonne mit Zone, Doppelsterne, Oortsche Wolke, Zwerge) sind beide gleich, die Zahlen dieser Notiz bleiben.

### 1.4 Erhaltungssätze

| Eigenschaft | Begründung | Status |
|---|---|---|
| Impuls, actio = reactio | V hängt nur von Differenzen x_i − x_j ab, auch über ℓ_i. Daraus folgt Σ F_i = 0. Für ein isoliertes Zweikörpersystem gilt F₁₂ = −F₂₁ exakt. Beispiel Sonne und ein Feldstern: Sie wechselwirken über V_N und über den kollektiven Term, beides symmetrisch. | exakt |
| Energie | L = T − V hängt nicht explizit von der Zeit ab. | exakt, solange ℓ_i als Funktion der Orte genommen wird. Läuft ℓ als langsamer Parameter nach (ṗ = εG), tauscht das System Energie der Ordnung ε mit p (adiabatisch). |
| Drehimpuls | K ist kugelsymmetrisch. | exakt |
| Galilei | V hängt nicht von den Geschwindigkeiten ab. ℓ_i hängt über T nur von Ortsdifferenzen ab. Ein gleichförmig bewegter Trichter verhält sich identisch. | exakt (nichtrelativistisch) |
| WEP innerhalb eines Trichters | Alle Glieder spüren das Elternfeld mit demselben ℓ. | exakt |

### 1.5 Was die Form kostet

**(a) Lokalität, unvermeidbar.** Im Vakuum ist Φ_N harmonisch. Der Mittelwert eines harmonischen Feldes über jede Kugel um x ist gleich seinem Wert in x. Deshalb kann **keine** lokale Mittelung das homogene äußere Feld vom eigenen Feld trennen. An einem Punkt sind beide ununterscheidbar, das ist das Äquivalenzprinzip. Die Trennung braucht die Kenntnis, wo die Quellen liegen. In E0 steckt diese in der Glättung um die Körper, also in einer Greenschen Funktion über die Zone.

**(b) ℓ-Kraft.** Die Strukturenergie hängt über ℓ(p) vom Ort des Trichters im Elternfeld ab. Das ergibt eine Kraft auf den Trichter als Ganzes, deren Gegenkraft auf die Eltern wirkt, sodass der Impuls erhalten bleibt:

| | dV/dℓ | Kraft je Masse (ℓ ∝ R^{2/3}) | relativ zu g_obs |
|---|---|---|---|
| E0 (Glättung im kollektiven Feld) | 1,59·10¹⁶ J/m | 8,96·10⁻¹⁹ m/s² | 4,2·10⁻⁹ |

Der Wert liegt weit unter jeder Messbarkeit.

**(c) WEP zwischen Körpern mit verschiedenem ℓ (Abschätzung).** Die Rechnung nimmt für die gelesene Dichte der Regel am Sonnenort (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt; kein Stoff) 0,01 M☉/pc³ an [G]. [korr.] Das Zusatzpotential der Regel in der Galaxie, über verschiedene ℓ gemittelt, ergibt für Sonne gegen Jupiter, Komet oder einen Stern mit 0,1 M☉ den Unterschied Δa ≈ 1,6–2,0·10⁻¹⁸ m/s², also Δa/g_obs ≈ 7–9·10⁻⁹. [korr.] Innerhalb des Sonnensystems gibt es diesen Effekt nicht, weil alle Glieder das Elternfeld mit ℓ_Sonne sehen.

**(d) Zerlegung in Körper.** E0 braucht keinen Baum, nur die Angabe, was ein „kompakter Körper“ ist, der geglättet wird. Topologische Änderungen wie Verschmelzen oder Zerfallen laufen über glatte Gewichte, damit V stetig bleibt.

**So wird der Preis klein gehalten:**
- glatte Plummer-Übergabe,
- ℓ nicht von Hand gesetzt, sondern aus der Gezeit (Parameter p),
- ein gemeinsames ℓ für alle Glieder eines Trichters,
- Q_s(0) = 0.

---

## 2. Kovariante Form: FTM-AS mit Glättungsoperator (FTM-AS-H)

In der kovarianten Wirkung FTM-AS ändert sich nur das Argument der Kinetikfunktion:

$$
S=\frac{c^3}{16\pi G}\!\int\! d^4x\sqrt{-g}\Big[R+2h^{\mu\nu}B_\mu B_\nu+\frac{2}{\ell_0^2}\,\tilde Q_s(\bar q)-\lambda_u(u_\mu u^\mu+1)\Big]+S_m,\qquad
\bar q=\exp\!\big[-(\ell_0^2\,\bar{\mathcal Y})^{1/4}\big],\quad \bar{\mathcal Y}=h^{\mu\nu}D_\mu\bar\sigma\,D_\nu\bar\sigma ,
$$

$$
\bar\sigma=\mathcal S[\sigma]:=D^{-2}\Big[\sum_i \mathcal K_{\ell_i}\big(\mathbb 1_i\,D^2\sigma\big)\Big]+D^{-2}\big[(1-\textstyle\sum_i\mathbb 1_i)\,D^2\sigma\big] .
$$

Die Zutaten:
- D² = D_μD^μ ist der räumliche Laplace-Operator auf dem Ruheschnitt des Mediums (Normale u^μ).
- 𝟙_i ist die charakteristische Funktion des i-ten kompakten Quellgebiets von D²σ, also des Gebiets, in dem der Trichter i seine eigene Faltung hat.
- 𝒦_ℓ ist der Plummer-Kern mit geodätischem Abstand auf dem Schnitt.
- Den Glättungsradius liefert die spurfreie Hesse-Form Ē^TF von c²σ̄, gemittelt über die eigene Zone, sodass nur die Elterngezeit bleibt: ℓ_i³ = GM_i/(√(3/2)·|⟨Ē^TF⟩_i|), und M_i ist der Gauß-Fluss von Dσ über das Quellgebiet. Die Gezeit kommt aus σ̄, nicht aus der Weyl-Projektion der Metrik (elektrischer Teil der Krümmung); so bleiben Ψ = Φ und c_T = c. [korr.]

Alles ist aus g, u und σ gebaut, ohne Kopplung an die Materie. S ist deshalb diffeomorphismeninvariant, und die Noether-Identität ∇^μ𝒮_μν = 0 gilt auf der Lösung. Der Strukturterm bleibt reine Geometrie, das entspricht der Beobachter-Regel.

Feldgleichung für σ:

$$
\nabla_\mu\big[h^{\mu\nu}(\nabla_\nu\sigma-a_\nu)\big]+\mathcal S^\dagger\Big\{\nabla_\mu\big[h^{\mu\nu}(\nu(\bar q)-1)\nabla_\nu\bar\sigma\big]\Big\}=0 .
$$

Im schwachen Feld ergibt sich exakt das Funktional aus 1.3, mit Ψ = Φ. Der Strukturteil koppelt nicht an die räumliche Metrik, also gilt für den Slip γ_s = Ψ/Φ = 1.

Einordnung:
- 𝒮 ist nichtlokal auf dem Ruheschnitt. Das passt zum ohnehin instantanen σ in FTM-AS.
- Die Glättung geschieht im Ruhesystem des Mediums. Ein gegen das Medium bewegter Trichter hat einen um Terme der Ordnung (v/c)² verformten Kern. Das ist ein Vorzugssystem-Effekt derselben Art wie die Parameter α₁, α₂ und nicht gerechnet.
- Die Zerlegung 𝟙_i von D²σ in kompakte Quellgebiete ist wohldefiniert für getrennte Körper, nicht aber für Gas und Flüssigkeiten. Für diffuse Materie gilt 𝟙 = 0 (keine Glättung).

---

## 3. Cassini, Quadrupol und Planeten

### 3.1 Q₂ als Funktion des Übergaberadius

Die Elternstruktur sieht die Sonne über ℓ = λ_ℓ r_eq geglättet. Der eigene Anteil der Sonne ist kugelsymmetrisch und trägt deshalb nichts bei. Mit g_e und r_eq aus Abschnitt 0 gilt

$$
Q_2=-\tfrac32\,\frac{g_e}{r_{eq}}\,J(\beta,\lambda_\ell),\qquad
J=\tfrac12\!\iint\!(\nu(\beta|\hat g|)-1)\Big[-\frac{3(3c^2-1)h}{\rho^2}+\frac{9c-15c^3}{\rho^2}\Big]d\rho\,dc+\frac{4}{15}\frac{\beta\nu'(\beta)}{\lambda_\ell^3},\qquad h=\frac{\rho}{(\rho^2+\lambda_\ell^2)^{3/2}} .
$$

Der Oberflächenterm entsteht, weil die geglättete Masse im Zentrum ein reguläres, anisotrop verstärktes Feld hat. Er ist die lineare EFE-Antwort (4/15)(dν/dln g)·GM/ℓ³ und erklärt den Abfall ∝ λ_ℓ⁻³.

**Kontrollen:** Für λ_ℓ = 0,5 / 1 / 3 stimmt die Quadratur mit einem unabhängigen Gitter-Hauptwert überein: −2,64982 gegen −2,64908, −0,14625 gegen −0,14612, −0,004557 gegen −0,004552.

**Cassini-Bedingung.** Cassini verlangt ℓ ≥ 1,62 / 1,42 / 1,30 r_eq (1σ / 2σ / 3σ), also ℓ ≥ 10 422 / 9 159 / 8 374 AE = 1,48 / 1,30 / 1,19 r\*.

| λ_ℓ = ℓ/r_eq | ℓ [AE] | Q₂ [s⁻²] (a₀ = 1,2·10⁻¹⁰) | Abstand zu Cassini [σ] | Q₂ (a₀ = cH₀/2π) |
|---|---|---|---|---|
| 43,9 (ℓ = r_J, klassischer Jacobi-Radius 1,37 pc) [korr.] | 283 034 | 2,90·10⁻³¹ | −1,0 | 2,82·10⁻³¹ |

Die Glättungslänge ist ℓ ≈ 0,93 pc = 27 r\* (Abschnitt 1.1); dort ist der Quadrupol aus der Glättung 9,3·10⁻³¹ s⁻². [korr.] E0 erfüllt die Cassini-Bedingung mit großem Abstand.

**Galaktische Gezeit.** Sie ist auch Newtonsch vorhanden. Der spurfreie Tensor ist diag(−2,86; −1,12; 3,98)·10⁻³⁰ s⁻². Auf die Hees-Form projiziert ergibt das Q₂ = 4,3·10⁻³⁰ s⁻², drei Größenordnungen unter der Messgenauigkeit.

### 3.2 Periheldrehungen Merkur bis Neptun

Die Periheldrehungen folgen aus den Gauß-Gleichungen, gemittelt über die Bahn, mit der Richtung zum galaktischen Zentrum in Ekliptik-Koordinaten (ekliptikale Länge 266,84°, Breite −5,54°). Probe mit konstanter Radialbeschleunigung: numerisch 4,126208·10⁻¹⁷, analytisch 4,126208·10⁻¹⁷ rad/s. Die Bahnelemente (J2000) und die EPM2011-Werte (Pitjeva & Pitjev 2013) sind [G].

| Planet | Δω̇_p je 10⁻²⁶ s⁻² | E0 + Gezeit | EPM2011 [G] | Q₂-Grenze (1σ) |
|---|---|---|---|---|
| Merkur | 1,10·10⁻⁴ | 5·10⁻⁸ | −0,020 ± 0,030 | 2,7·10⁻²⁴ |
| Venus | 5,43·10⁻⁵ | 2·10⁻⁸ | 0,026 ± 0,016 | 2,9·10⁻²⁴ |
| Erde | 4,22·10⁻⁴ | 2·10⁻⁷ | 0,0019 ± 0,0019 | 4,5·10⁻²⁶ |
| Mars | −4,16·10⁻⁴ | −2·10⁻⁷ | −0,00020 ± 0,00037 | 8,9·10⁻²⁷ |
| Jupiter | −3,00·10⁻³ | −1·10⁻⁶ | 0,587 ± 0,283 | 9,4·10⁻²⁵ |
| Saturn | 1,39·10⁻² | 6·10⁻⁶ | −0,0032 ± 0,0047 | 3,4·10⁻²⁷ |
| Uranus | −2,67·10⁻² | −1·10⁻⁵ | −3,89 ± 3,90 | 1,5·10⁻²⁴ |
| Neptun | 2,71·10⁻² | 1·10⁻⁵ | −4,44 ± 5,40 | 2,0·10⁻²⁴ |

Alle Werte in mas/Jahr, die letzte Spalte in s⁻².

Die Periheldaten allein ergeben aus Saturn eine Grenze von 3,4·10⁻²⁷ s⁻² für Q₂. Das passt in der Größenordnung zur Cassini-Grenze von Hees et al. (2014) und ist eine unabhängige Plausibilitätsprobe der erinnerten EPM-Werte. E0 liegt rund drei Größenordnungen oder mehr unter jeder Grenze, bei Mars 3,3 und bei Saturn 2,9 Größenordnungen. INPOP-Auswertungen liefern Schranken derselben Größenordnung [G, Werte nicht verifiziert].

**Zusatzbeschleunigungen.** Der Quadrupol erzeugt |Q₂|·a, in E0 von 2,7·10⁻¹⁹ m/s² (Merkur) bis 2,1·10⁻¹⁷ m/s² (Neptun).

Der eigene Strukturterm (ν − 1)g_N der Sonne ist bei Neptun 10^(−106,6) m/s² und bei Merkur 10^(−7888) m/s², also praktisch null.

---

## 4. Vorhersagen

### 4.1 Weite Doppelsterne

Gerechnet ist die kugelgemittelte Einwärtsbeschleunigung. Sie ist über die Flussidentität ∮∇Φ·dA = ∮ν∇Φ_N·dA exakt. Das Ergebnis ist das Verhältnis γ = ⟨g⟩/g_N, bei M_tot = 1,5 M☉ und a₀ = 1,2·10⁻¹⁰.

| d [kAE] | g_N [m/s²] | **E0** |
|---|---|---|
| 2 | 2,22·10⁻⁹ | 1,000 |
| 3 | 9,88·10⁻¹⁰ | 1,000 |
| 5 | 3,56·10⁻¹⁰ | 1,000 |
| 7 | 1,82·10⁻¹⁰ | 1,000 |
| 10 | 8,90·10⁻¹¹ | 1,000 |
| 15 | 3,95·10⁻¹¹ | 1,000 |
| 20 | 2,22·10⁻¹¹ | 1,000 |
| 30 | 9,88·10⁻¹² | 1,000 |

Das Geschwindigkeitsverhältnis ist v/v_N = √γ = 1,000 bei 10, 20 und 30 kAE.

Zusätzliche Kenngrößen:
- **Isotropie.** In E0 ist γ innerhalb der Zone isotrop; eine Richtungsabhängigkeit relativ zum galaktischen Zentrum tritt nicht auf.
- **Galaktische Gezeit.** Bei 30 kAE beträgt sie 2,7·10⁻³ von g_N.

**Datenlage [G]:**
- *Chae (2023, ApJ 952, 128)* findet γ_g = 1,43 ± 0,06 für g_N ≲ 10⁻¹⁰ m/s². Chae (2024) und Hernandez (2023, MNRAS) finden Ähnliches.
- *Banik et al. (2024, MNRAS 527, 4573)* finden Newton bevorzugt, mit hoher Signifikanz (16–19 σ) gegenüber MOND mit EFE.
- Umstritten sind die Reinheit der Stichproben (unaufgelöste Dreifachsysteme), Projektionseffekte und Geschwindigkeitsfehler.

E0 sagt γ = 1,000 voraus; das entspricht dem Befund von Banik et al. Die Messlage ist offen, die Entscheidung wird Gaia DR4 bringen (5.3).

### 4.2 Äußere Oortsche Wolke, 10⁴ bis 10⁵ AE

**Kreisbahnen, γ = ⟨g⟩/g_N:**

| r [AE] | E0 | Gezeit/g (Newton) |
|---|---|---|
| 10 000 | 1,000 | 1,5·10⁻⁴ |
| 50 000 | 1,002 | 1,9·10⁻² |
| 100 000 | 1,014 | 0,151 |

**Scheinbare Bahnenergie.** Die „ursprüngliche“ Energie 1/a_orig eines langperiodischen Kometen wird bei etwa 250 AE bestimmt, wo das Strukturpotential konstant ist. Dort gilt Newton. Ein Komet mit wahrem Aphel r_A hat deshalb die Newtonsche Energie

$$
E_N=-\frac{GM}{r_A}+\big[\Phi_s(r_A)-\Phi_s(0)\big],\qquad \frac1{a_\text{app}}=-\frac{2E_N}{GM}.
$$

| wahres Aphel r_A [AE] | Newton | E0 |
|---|---|---|
| 10 000 | 20,0 | 20,0 |
| 20 000 | 10,0 | 10,0 |
| 50 000 | 4,0 | 4,0 |
| 100 000 | 2,0 | 1,99 |

Alle Werte für 1/a_app in 10⁻⁵ AE⁻¹. Die Oort-Spitze bei 1/a = 3·10⁻⁵ AE⁻¹ entspricht E_N = −1,33·10⁴ m²/s².

**Gezeiten-Einschuss über die Jupiter-Saturn-Barriere.** Das Maß ist die Änderung der Wurzel der Periheldistanz r_P, Δ√r_P ≥ √15 − √5 = 1,64 √AE je Umlauf, mit Δ√r_P ≈ ½ν_z²⟨r²⟩P/√(2GM):

| | r_A,krit [AE] | treibend | 1/a_app(r_A,krit) [10⁻⁵ AE⁻¹] |
|---|---|---|---|
| Newton | 47 005 | Gezeit | +4,25 (Probe: a_krit ≈ 2,4·10⁴ AE, Literatur ~2–3·10⁴ [G]) |
| E0 | 47 014 | Gezeit | +4,25 |

**Befund.** E0 ist bis 10⁵ AE Newtonsch (γ = 1,014 bei 10⁵ AE) und reproduziert Newton. Neu eingeschossene Kometen erscheinen gebunden mit 1/a_app ≈ +4,25·10⁻⁵ AE⁻¹, in Übereinstimmung mit der beobachteten Spitze bei +2…5·10⁻⁵ AE⁻¹ [G].

*Vorbehalte:* Die Bahnen sind radial genähert, und das Potential ist richtungsgemittelt. Die Aussage über die scheinbare Energie ist dagegen robust, denn sie folgt nur aus ΔΦ_s(r_A).

### 4.3 Zwerggalaxien im Feld großer Galaxien

**Crater II ist ein Satellit der Milchstraße**, DF2 und DF4 gehören zur NGC-1052-Gruppe. Als Andromeda-Beispiel dient And XIX, dessen Eingangswerte unsicher sind [G].

Der Schätzer ist σ² = r_½·g(r_½)/3 (Wolf et al. 2010 [G]) mit M(<r_½) = M/2. Die Beschleunigung g folgt aus derselben Flussformel wie in 4.1.

| System | r_½ [kpc] | g_Ne/a₀ | ℓ_J [kpc] | Newton | isoliert | E0 | σ_obs [G] |
|---|---|---|---|---|---|---|---|
| Crater II (MW) | 1,42 | 0,0050 | 0,67 | 0,40 | 4,11 | 1,41 | 2,7 ± 0,3 |
| And XIX (M31), unsicher | 2,13 | 0,0132 | 0,76 | 0,52 | 5,17 | 1,45 | 4,7 +1,6/−1,4 |
| NGC 1052-DF2 | 2,93 | 0,0181 | 4,01 | 6,99 | 21,1 | 18,0 | 8,5 +2,3/−3,1 |
| NGC 1052-DF4 | 2,13 | 0,0043 | 5,97 | 7,10 | 19,8 | 19,7 | 4,2 +4,4/−2,2 |

Alle σ in km/s.

- **E0 für Zwerge.** Ihre Sterne haben kleine ℓ, die Zwerggalaxie ist aufgelöst, und es wirkt die volle Regel im Gesamtfeld mit äußerem Feld. Die Bahnen folgen der örtlichen Regel; für kugelnahe Zwerge ist sie gleich der Feldform (QUMOND mit EFE), die Zahlen der Tabelle bleiben. [korr.]
- **Normierung.** Der Absolutwert des Schätzers ist im EFE-dominierten Fall um etwa √2 unsicher. McGaugh (2016) sagt mit seiner Normierung 2,1 km/s für Crater II voraus [G]. Die Verhältnisse zwischen den Spalten sind belastbar.
- **DF2 und DF4** liegen in E0 bei 18,0 bzw. 19,7 km/s, gemessen sind 8,5 bzw. 4,2. Famaey, McGaugh & Milgrom (2018) erhalten mit EFE 13–14 km/s für DF2 [G]. Der Unterschied liegt bei Wirtsmasse, 3D-Abstand, Schätzer und Interpolationsfunktion. Die Spannung betrifft jede QUMOND-artige Fassung mit EFE gleichermaßen und ist keine Frage der Hierarchie; die Diskussion um Entfernung, Neigung und Gezeitenzustand ist offen [G].

---

## 5. Folgerung

### 5.1 Verträglichkeit

- **Cassini und EPM:** erfüllt, Q₂ drei Größenordnungen und mehr unter jeder Grenze.
- **Oort-Spitze:** erfüllt, Newtonsch bis 10⁵ AE.
- **Weite Doppelsterne:** γ = 1,000, im Einklang mit Banik et al.; Chae und Hernandez finden größere Werte. Die Messlage ist offen.
- **Zwerggalaxien:** volles kollektives Feld mit EFE; Crater II im Rahmen der Normierungsunsicherheit, DF2 und DF4 offen für jede QUMOND-artige Fassung.

### 5.2 E0 als die im Modell angelegte Einzeltrichter-Hierarchie

E0 formuliert drei Aussagen des Modells:
1. Ein Körper mit g(R) > a₀ ist ein Einzeltrichter und zeigt keinen „Dunkle-Materie“-Effekt. Beispiele sind die Sonne und Sterne. Kugelsternhaufen sind als Sternsysteme kollektiv, kompakte (g/a₀ ≈ 23) zeigen aber kaum einen Zusatz (ν = 1,008). [korr.] Diffuse Systeme mit g(R) < a₀ sind kollektiv. In E0 gilt das nicht nur an der Oberfläche, sondern für die ganze Zone des Körpers bis zu seiner Gezeitengrenze.
2. „Dunkle Materie ist der Summeneffekt der im Raum liegenden Zusatzstrecken.“ Verstärkt wird in E0 genau dieser Summeneffekt, also das über die Einzeltrichter geglättete kollektive Feld, nicht die Einzelfalte.
3. „Einzeltrichter bis zur Schwelle, danach kollektives Netzwerk-Faltungsfeld.“ Für kompakte Körper liegt die Übergabe an der Glättungslänge ℓ (nahe der Gezeitengrenze), nicht bei r\*. Für die Sonne ist ℓ ≈ 0,93 pc = 27 r\*; der klassische Jacobi-Radius ist 1,37 pc. [korr.]

   Für Galaxien bleibt r\* die Schwelle des kollektiven Feldes, und RAR, Tully-Fisher und Faltungsmesser gelten unverändert. Für Kugeln und außerhalb von Scheiben ist η = M_dyn/M_lens = 1; in Scheiben innen (R ≲ 3,5 R_d) ist η > 1. [korr.]

Fraktal ist E0 in dem Sinn, dass jede Stufe dieselbe Regel hat: Ein Trichter ist bis zur Gezeitengrenze im Feld seiner Eltern Newtonsch. Die Eltern sehen ihn geglättet über diese Grenze. Verstärkt wird nur, was auf der jeweiligen Stufe aufgelöst verteilt ist. Das äußere Feld ist langsamer Parameter p und bestimmt ℓ_i, und das ist seine einzige Wirkung auf das Innere, neben den Gezeiten.

### 5.3 Prüfung von E0

**Weite Doppelsterne bei 5–30 kAE.** E0 sagt γ = 1,000 voraus (Rest der Glättung ≤ 3·10⁻⁴, Gezeit 2,7·10⁻³ bei 30 kAE), isotrop.

Ein sauberer Nachweis von γ ≳ 1,1 bei d ≳ 7 kAE ist mit E0 unvereinbar. Voraussetzung ist, dass Dreifachsysteme kontrolliert sind und die Stichproben der konkurrierenden Gruppen denselben Befund zeigen. Gaia DR4 ist für Dezember 2026 angekündigt [G].

**Quadrupol im Sonnensystem.** E0 erlaubt nur die galaktische Gezeit von Q₂ = 4,3·10⁻³⁰ s⁻² (dazu 9,3·10⁻³¹ s⁻² aus der Glättung bei ℓ = 0,93 pc). [korr.] Ein anomaler Q₂ ≳ 10⁻²⁷ s⁻² in Richtung des galaktischen Zentrums, etwa aus INPOP/EPM mit BepiColombo, wäre mit E0 ebenso unvereinbar.

---

## 6. Offene Punkte

1. **Körperzerlegung kovariant.** Die Aufteilung von D²σ in kompakte Quellgebiete 𝟙_i ist für Gas und für verschmelzende Systeme nicht eindeutig. Im Modell wird nur der Überschuss (ρ − ρ_J)₊ der Quelldichte über die Gezeitendichte ρ_J = 3T/(4πG) geglättet; a₀ kommt darin nicht vor, der Einsatz ist stetig. Sterne sind Einzeltrichter, Sternsysteme kollektiv; ein kompakter Kugelsternhaufen (g/a₀ ≈ 23) hat ν = 1,008. Offen ist die Form des Übergangs zwischen Einsatz und voller Glättung. [korr.]
2. **Definition von ℓ.** Gewählt ist hier die Jacobi-Grenze mit Newtonscher Eigenmasse und effektiver Elterngezeit 4A(A−B) (1,37 pc); im Modell gilt ℓ³ = GM/T mit T = √(3/2)·|E^TF| (0,93 pc). [korr.] Eine Alternative ist die effektive Eigenneigung. Sie ist für E0 unkritisch, weil Zwerge über ihre Sterne aufgelöst werden.
3. **Aufgelöste, aber dünne Sternsysteme.** Dünne Kugelsternhaufen wie Pal 14, Pal 4 und NGC 2419 [G] sind in E0 aufgelöst und werden kollektiv verstärkt, mit EFE. Die gemessenen, eher Newtonschen Dispersionen [G] sind ein eigener Test für E0.
4. **Vakuumterm.** Die Hierarchie erzwingt Q_s(0) = 0 je Trichter (1.2). Ob Λ_s = (π²/15)H₀²/c² als globale Größe physikalisch ist, bleibt offen.
5. **Galaxienhaufen.** Der Massenfehlbetrag in Galaxienhaufen bleibt bestehen, weil E0 auf Galaxien- und Haufenskala das volle kollektive Feld ist.
6. **DF2 und DF4.** Unabhängig von der Hierarchie zu klären: Entfernung, Neigung, Gezeitenzustand.
7. **Rückwirkung auf FTM-AS.** Der Glättungsoperator 𝒮 ersetzt das Argument von q. Die klassischen Tests (γ = β = 1, Lichtablenkung, Merkur) sind davon unberührt, weil q im Sonnensystem 10⁻¹⁰⁰ oder kleiner ist.

---

## Quellen

- Hees, A. et al. (2014), *Phys. Rev. D* 89, 102002 (Cassini, Q₂)
- Milgrom, M. (2010), *MNRAS* 403, 886 (QUMOND)

Aus dem Gedächtnis zitiert [G], nicht nachgeprüft:
- Pitjeva, E. V. & Pitjev, N. P. (2013), *MNRAS* 432, 3431 (EPM2011, Periheldrehungen)
- Chae, K.-H. (2023), *ApJ* 952, 128; Chae, K.-H. (2024), *ApJ* 960, 114 (weite Doppelsterne)
- Hernandez, X. (2023), *MNRAS* 525, 1401 (weite Doppelsterne)
- Banik, I. et al. (2024), *MNRAS* 527, 4573 (weite Doppelsterne, Newton bevorzugt)
- Heisler, J. & Tremaine, S. (1986), *Icarus* 65, 13; Duncan, M., Quinn, T. & Tremaine, S. (1987), *AJ* 94, 1330 (Gezeiten-Einschuss, Oort-Spitze)
- Bovy, J. (2017), *MNRAS* 468, L63 (Oortsche Konstanten); GRAVITY Collaboration (2019), *A&A* 625, L10 (R₀)
- Torrealba, G. et al. (2016), *MNRAS* 459, 2370; Caldwell, N. et al. (2017), *ApJ* 839, 20; McGaugh, S. (2016), *ApJL* 832, L8 (Crater II)
- Collins, M. et al. (2013), *ApJ* 768, 172 (And XIX)
- van Dokkum, P. et al. (2018), *Nature* 555, 629; Danieli, S. et al. (2019), *ApJL* 874, L12; van Dokkum, P. et al. (2019), *ApJL* 874, L5 (DF2, DF4)
- Famaey, B., McGaugh, S. & Milgrom, M. (2018), *MNRAS* 480, 473; Kroupa, P. et al. (2018), *Nature* 561, E4; Haslbauer, M. et al. (2019), *A&A* 626, A47 (DF2/DF4 und EFE)
- McGaugh, S. & Milgrom, M. (2013), *ApJ* 766, 22 (EFE-Formel für Zwerge)
- Wolf, J. et al. (2010), *MNRAS* 406, 1220 (Massenschätzer)
