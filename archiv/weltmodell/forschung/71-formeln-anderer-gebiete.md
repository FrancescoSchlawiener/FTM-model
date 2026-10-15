# Formeln anderer Gebiete an den Invarianten des Modells

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Paul-Falle „κ² + ν_z² − 2Ω² − 4πGρ_b = χΩ²“ → gilt nur außerhalb der Masse, in der Masse 4πGρ_X = χΩ² + 4πGρ_b·∂(a − a_N)/∂a_N; Verweilzeit „a₀ = c²/U“ → a₀ = c²/(2Λ), Λ keine Raumlänge; SPARC „a₀ = 1,22·10⁻¹⁰“ als Wert ohne Fehlermodell gekennzeichnet, mit den Fehlern von Abstand, Neigung und Masse-zu-Licht 1,04–1,19·10⁻¹⁰, ohne Kreispaare ≤ 1,13·10⁻¹⁰. Umbenennungen: w (Glättungstiefe) → ϖ; U → ℓ₀; Thiele-Wirkungsgrad η → η_Th; Weglängenfaktor η → η_d; kinematische Zähigkeit ν → ν_kin; Ladung q (Paul-Falle) → Q_el; Förder- und Durchflussrate Q → V̇; mittlere Zahl L in Littles Gesetz → N̄.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/formeln_andere_gebiete.py` (Ausgabe `.txt`). [G] = Messwert aus dem Gedächtnis, nicht nachgeprüft.
> Grad: **gleich** = dieselbe Gleichung; **Struktur** = gleiche Form, andere Herkunft; **ähnlich** = nur Anklang. Übertragen wird nur Regel oder Form, nie ein Stoff oder Träger.

## 1  Abgeschirmte Rückkehr: D·G'' = Γ·G, q = e^{−L√(Γ/D)}, s² = Γ·τ_D

- **Thiele-Modul** (Verfahrenstechnik, Katalysatorkorn) · D c'' = k c, φ = L√(k/D); dichte Mitte: cosh(φ(1 − x/L))/cosh φ, Wirkungsgrad η_Th = tanh φ/φ [korr.] · **gleich** (s = φ).
  Gewinn: Mit der Aris-Länge L = V/S werden Platte, Zylinder und Kugel im abgeschirmten Bereich formgleich; dieselbe Größe V/(S·R) = 1/d steht in Δs, ein Formbeiwert für Scheibe gegen Kugel. Im Diffusionsbereich halbieren sich scheinbare Ordnung und Aktivierungsenergie (Rate ∝ √(kD)), wie d ln a/d ln a_N → ½.
  Bruch: Thiele erhält die Wurzel im Fluss bei großem φ, das Modell in der Summe 1/(1 − q) bei kleinem s: gleiche Wurzel, andere Größe (Struktur, nicht gleich).
- **Kabelgleichung** (Neurobiologie) · λ²V'' = V, λ = √(r_m/r_i); unendlich e^{−x/λ}, versiegeltes Ende cosh((L − x)/λ)/cosh(L/λ), Eingangsleitwert G_∞·tanh(L/λ) · **gleich** (s = elektrotonische Länge).
  Gewinn: Das versiegelte Ende ist die bei Λ geschlossene Rückkehr. Gerechnet: 10 % Abweichung von e^{−s} ab r = 0,64–0,65·r*·Λ/L, tief 1 − q → s²(Λ/L − ½), a konstant (Flussdimension 1). Ralls Regel d^{3/2} = Σd_i^{3/2} sagt, wann ein verzweigter Rückweg wie eine einzige Halbgerade wirkt.
  Bruch: Die Membran leckt über Kanäle; im Modell ist der Verlust Γ = a_N/c ein Taktgefälle, kein Leck.
- **Warburg-Impedanz** (Elektrochemie) · halbunendlich Z ∝ 1/√(iω) (45°); endlich reflektierend Z ∝ coth√(iωτ_d)/√(iωτ_d), durchlässig tanh-Form; τ_d = L²/D · **gleich** mit Γ → iω, s² ↔ iωτ_d.
  Gewinn: fertiges Messverfahren, das offene, reflektierend und absorbierend geschlossene Halbgerade am Tieffrequenzende trennt (45° bleibt, → 90°, → reell); der Knick 45° → 90° ist der Übergang Tully-Fisher → konstante Beschleunigung bei Kac-Schließung. Als Frequenzantwort hat die Rückkehr den Verlustwinkel 45° (Speicher = Verlust, wie das kritische Gel nach Winter–Chambon mit n = ½).
  Bruch: Γ ist im Modell reell.
- **Morphogen-Gradient** (Entwicklungsbiologie) · ∂_t c = Dc'' − kc, λ = √(D/k); mittlere Aufbauzeit bei fester Quelle τ(x) = x/(2√(Dk)) · **gleich**.
  Gewinn (gerechnet): τ(L) = τ_D/(2s) ist exakt die Sattelzeit t* = ½√(τ_D/Γ) der Kostenminimierung s = min_t[Γt + τ_D/4t]. Wäre die Rückkehr ein Vorgang in t, stünde nach 13,8 Mrd. Jahren (τ_D = 79,2 Mrd.) erst 22 % des stationären q bei s = 1, 12 % bei s = 0,3. Randbedingung: ϖ ist keine Zeitrichtung, Stationarität ist verlangt. [korr.]
  Bruch: Im Gewebe läuft echte Zeit.
- **Skin-Effekt, Stokes-Schicht, Womersley** (Elektrotechnik, Strömung, Blutfluss) · e^{−(1+i)y/δ}, δ = √(2D/ω); α_W = R√(ω/ν_kin) [korr.] · **Struktur** (Γ komplex).
  Gewinn: zeigt, was eine Phase täte: |q| und arg q wären gleich, die Umlaufphase Re τ_M hinge an s. Im Modell ist Re τ_M = ω_rad/ω_azim ∈ [1, √2] frei von s; das bestätigt ein reelles Γ.
  Bruch: äußerer periodischer Antrieb.
- **Leckender Grundwasserleiter, Debye-Hückel, London, Brinkman** (Hydrogeologie, Elektrolyte, Supraleitung, poröse Medien) · (V̇/2πT)·K₀(r/B) [korr.], B = √(T·b'/K'); e^{−r/λ_D}/r; e^{−x/λ_L}; √K · **Struktur**, Ausbreitung in 2 bzw. 3 Dimensionen.
  Gewinn: Randbedingung: K₀ und e^{−r}/r divergieren am Ort, nur die eindimensionale (ebene) Ausbreitung erlaubt G(0) = 1. Die London-Tiefe ist die Starrheit der Phase des geordneten Mediums, kein Stoff: ein fertiges Beispiel für „Regel statt Stoff“.
- **Fortsetzung nach oben** (Geophysik, Schwere- und Magnetfeld) · F̂(k, h) = F̂(k, 0)·e^{−|k|h}, P_a·P_b = P_{a+b}, nach unten instabil; Tiefe aus ln P(k) = const − 2|k|h (Spector–Grant) · **gleich**.
  Gewinn: ϖ ist die Fortsetzungshöhe [korr.], q = e^{−r*/r} die Kreismode k = 1/r in der Höhe r*. Messverfahren: ln q gegen 1/r ist eine Gerade mit Steigung −r*. Satz der Äquivalenzschicht: Jedes Außenfeld entsteht aus unendlich vielen Dichteverteilungen; die gelesene Dichte ρ_X ist nicht eindeutig und darum kein Befund über Inhalt.
  Bruch: Dort ist h eine Raumhöhe über echten Quellen.

## 2  Rückkopplungssumme Σqⁿ = 1/(1 − q) = 1 + n_B, Kippstelle q = 1

- **Akkumulation bei wiederholter Gabe** (Pharmakokinetik) · R = 1/(1 − e^{−kτ}); Spitze C₀R, Tal C₀e^{−kτ}R; kτ = −ln(1 − 1/R) · **gleich** (s = kτ).
  ν = 1 + n_B heißt: Spitze = frische Gabe + Rest aller früheren (Tal = n_B). Newton: zwischen den Gaben ganz abgebaut; Tully-Fisher: Gabe viel schneller als Abbau, R ≈ 1/s.
  Gewinn (gerechnet an SPARC, 2573 Punkte, Υ = 0,5/0,7): s_beob = −ln(1 − g_bar/g_beob) liest s direkt; s_beob²·a₀/g_bar liegt für g_bar < 10⁻¹⁰ m/s² im Median bei 1,01–1,14 mit a₀ = 1,22·10⁻¹⁰ (ohne Fehlermodell) [korr.]. Linearität in g_bar prüft die Exponentialform, ein Anstieg bei den kleinsten g_bar das äußere Feld. Brauchbar nur für s ≲ 2.
- **Quellvervielfachung** (Reaktorphysik) · M = 1/(1 − k_eff), Annäherung an k = 1 über die Gerade 1/M; Periode ℓ/(k − 1) · **gleich**.
  Gewinn: dieselbe Linearisierung; „Verstärkung ist Verlangsamung“ (τ₀/(1 − q)) ist dort Messpraxis. Bruch: Der Reaktor hat einen überkritischen Ast, das Modell mit q = e^{−s} nie.
- **Keynes-Multiplikator, Leontief-Inverse** (Wirtschaft) · 1/(1 − c); x = (I − A)⁻¹y = ΣAⁿy, produktiv genau für Spektralradius < 1 (Hawkins–Simon); mittlere Fortpflanzungslänge ΣnAⁿ/ΣAⁿ · **gleich**.
  Gewinn: Zerlegung in direkte und indirekte Anteile im Netz; die mittlere Rundenzahl ist Q = ν − 1. Das Sparparadoxon ist der Beobachterfehler in Reinform: Einzelregel gegen Kreislaufregel. Bruch: c ist Verhalten, q folgt aus der Geometrie.
- **Selbsterregte Punktprozesse** (Hawkes; Nachbeben ETAS; Märkte) · λ̄ = μ/(1 − n*), Anteil endogener Ereignisse n* · **gleich** im Mittel.
  Gewinn: Lesart des kosmischen q = 0,843: Die Baryonen sind die exogenen Ereignisse (1 − q = 0,157 = Ω_b/Ω_m), die Zusatzwirkung ist die Wiederholung derselben Ereignisse, keine zweite Sorte. Märkte n* ≈ 0,7–0,9 [G].
  Bruch: Hawkes ist zufällig (Fano-Faktor 1/(1 − n*)²); das Modell ist deterministisch, dieses Rauschen darf nicht übertragen werden.
- **Fabry–Pérot** (Optik) · Σ(Re^{iδ})ⁿ, T = 1/(1 + F sin²(δ/2)), F = 4R/(1 − R)²; Lebensdauer ≈ τ_Umlauf/(1 − R) · **gleich** mit komplexem q.
  Gewinn: Mit Phase läge |1/(1 − q)| zwischen 1/(1 + |q|) und 1/(1 − |q|), also auch unter 1; die Verzögerung ist ν mal Umlaufzeit (im Modell ν·r/c). Bruch: Im Modell trägt nur |q|.
- **Kleinverstärkungssatz, May-Stabilität** (Regelungstechnik, Ökologie) · Kreis stabil für Produkt der Verstärkungen < 1; zufälliges Netz aus N Arten mit Verknüpfung C und Stärke σ stabil für σ√(NC) < 1 · **Struktur**.
  Gewinn: eine Komplexitätsgrenze für das Netz der Rückkopplungsschleifen, NC·σ² < 1, neben der Grenze N_max ≈ 1 + 2W/(πk) für Zyklen. Bruch: May mittelt über Zufallsmatrizen, der Kopplungsgraph des Modells ist gebaut.

## 3  Der Exponent ½: Dimension und Kritikalität der Rückkehr

- **Pólya, Spektraldimension** (Statistische Physik, Netzwerke) · Rückkehrdichte ∝ t^{d_s/2 − 2} ⇒ 1 − q ∝ Γ^{1 − d_s/2} ⇒ a ∝ a_N^{d_s/2} · **gleich** mit d_s = 1.
  Gerechnet: 1D 1 − F = √(2Γ) (Steigung 0,500); 2D 1 − F = π/ln(8/Γ) (Steigung 0,93 bei Γ = 10⁻⁶, → 1); 3D Rückkehr 0,3405, Verstärkung höchstens u₃ = 1,516 (Watson).
  Gewinn: Messverfahren d_s = 2·d ln a/d ln a_N im Tiefbereich; Tully-Fisher verlangt d_s = 1 genau. Randbedingung: Eine räumliche Rückkehr trägt keine flachen Kurven. Sichere Rückkehr heißt unendlicher effektiver Widerstand ins Unendliche (Doyle–Snell).
- **Kritische Verzweigung, neuronale Lawinen** (Neurobiologie) · Galton–Watson mit Mittel m; bei m = 1 Größen ∝ S^{−3/2}, 1 − E[e^{−ΓS}] ≈ √(2Γ); bei m < 1 linear Γ/(1 − m) (gerechnet); Kortex: Exponent ≈ 3/2, Verzweigung ≈ 1 [G] · **Struktur**.
  Gewinn: Tully-Fisher heißt: Die Kaskade der Rückkehr ist genau kritisch (Schritte auf und ab gleich häufig). Jede Drift gibt einen Boden; die Lawinenstatistik (Größen-, Dauerexponent, Skalenfunktion) beschreibt die Form des Übergangs.
  Bruch: Im Modell gibt es keine einzelnen Rückläufe als Ereignisse, nur ihre Summe.
- **Kleiber, West–Brown–Enquist** (Biologie) · B ∝ M^{d/(d+1)} = M^{3/4}, Zeiten ∝ M^{1/(d+1)} · **ähnlich**.
  Gleiche Viertelpotenz wie v ∝ M^{1/4}, andere Herkunft: Im Modell ist 4 = 2·2 aus der Wurzel der Rückkehr und v² = a·r; in d Raumdimensionen gäbe es überhaupt nur für d = 3 flache Kurven (v² ∝ √M·r^{(3−d)/2}), keine Reihe M^{1/(d+1)}.
  Gewinn: Das Netz hat zwei Flussregime (pulsierend und Poiseuille, Umschalten bei α_W ≈ 1) und daher eine gekrümmte Skalenlinie [G], wie χ zwischen Flussdimension 3 und 2. Bruch: WBE braucht ein stoffliches Netz; 1/4 = 1/(3 + 1) ist kein Beleg.

## 4  Bose-Leiter s^d·n_B und Γ(d+1)ζ(d+1)

- **Planck, Stefan–Boltzmann, Wien, Rayleigh–Jeans** (Strahlung) · u ∝ x³n_B, ∫ = π⁴/15, Wien-Stelle 2,821; in d Dimensionen T^{d+1} · **gleich** (Funktion).
  Gewinn: Grenzfälle lesbar: Tiefbereich n_B ≈ 1/s ist Gleichverteilung (Rayleigh–Jeans), Newton der Wien-Schwanz e^{−s} (die Rückläufe frieren aus). Ohne den Schwanz divergierte die Leiter wie die Ultraviolett-Summe; die Exponentialform sichert die Newton-Grenze. Bruch: keine Temperatur.
- **Debye-Funktion** (Festkörper) · D_d(x) = (d/x^d)∫₀^x t^d/(e^t − 1)dt, x = Θ_D/T · **gleich**: Jeder monotone Ast von s(r) liefert (s_max^d/d)·D_d(s_max).
  Gewinn: Kompaktheit ablesen wie Θ_D aus der Wärmekapazität; diffus (s_max ≪ 1) ≈ s_max^d/d, kompakt → Γ(d+1)ζ(d+1); gerechnet: d = 3 bei s_max = 3 erst 39 % des Grenzwerts. Bruch: Debye schneidet bei fester Modenzahl, das Modell bei der größten Beschleunigung.
- **Geordnete Korrelation, KMS** (Quantenoptik, Statistische Physik) · Emission ∝ 1 + n, Absorption ∝ n, G^>/G^< = e^s; symmetrisch (Fluktuation–Dissipation) coth(s/2) = 1 + 2n_B · **Struktur**.
  Gewinn: ν = 1 + n_B ist die geordnete Funktion (ein Sinn), coth die symmetrische (beide Sinne), q = n_B/(1 + n_B) das KMS-Verhältnis. Die Wahl des Umlaufsinns ist dieselbe Wahl; für Übergänge eines Detektors (Unruh–DeWitt) ist die geordnete Funktion die richtige. Bruch: keine Wärme; das gibt eine Sprache, keinen Beweis.
- **Maximale Entropie der Rücklaufzahl** (Informationstheorie) · P(n) = (1 − q)qⁿ ist bei festem Mittel Q die Verteilung größter Entropie (Jaynes) · **gleich**.
  Gewinn (gerechnet): H = ln ν + χ = Φ_X/v_f² + χ (s = 1: 1,0407 nats; tief ≈ ln(1/s) + 1). Die Bose-Form folgt ohne Wärme aus der Zählung; Zusatzpotential und Schalter sind freie Energie und Energie einer Zählverteilung. Bruch: MaxEnt braucht eine Nebenbedingung mit Multiplikator s; im Modell kommt s aus der Geometrie.

## 5  Flussdimension, Epizyklen, gelesene Dichte

- **Paul-Falle, Kapitza-Pendel** (Ionenfallen) · ladungsfrei verlangt Earnshaw Σω_i² = 0; im Wechselfeld Pseudopotential Ψ = Q_el²|∇φ₀|²/(4mΩ²) [korr.], Σω_i² > 0 · **Struktur**.
  Wer Earnshaws Summe erwartet, liest eine Ladung, die nicht da ist, wie κ² + ν_z² = (2 + χ)Ω² → 4πGρ_X = χΩ². Die Wirkung hängt am Quadrat des Gradienten, wie die Strukturwirkung an s⁴ = ℓ₀²|∇σ|².
  Gewinn: Laborfall des Beobachterfehlers mit messbaren Frequenzen. Nach dem Invarianzsatz ω_c² = ω₊² + ω_z² + ω₋² (Brown–Gabrielse) misst κ² + ν_z² − 2Ω² = χΩ² außerhalb der Masse den Schalter am Ort aus drei Frequenzen; in der Masse (kugelsymmetrisch) ist κ² + ν_z² − 2Ω² − 4πGρ_b = 4πGρ_X = χΩ² + 4πGρ_b·∂(a − a_N)/∂a_N. [korr.] Bruch: Die Falle braucht einen äußeren Antrieb.
- **Effektive und mitgeführte Masse** (Festkörper, Strömung) · m* = ħ²/(∂²E/∂k²); Kugel in Flüssigkeit m + ½ρV · **Struktur**: Masse als Ausgabe der Bewegungsregel (M_schein/M = ν).
  Gewinn: Beide sind Tensoren; zu erwarten ist eine richtungs- und formabhängige scheinbare Masse (Scheibe gegen Kugel, √(cos θ)). Bruch: Dort ändert die Periodizität bzw. die Verdrängung die Regel, im Modell die Rückkehr.

## 6  Zusammensetzung im Betrag: s⁴ = s_i⁴ + s_e⁴

- **Matthiessen-Regel, Debye-Abschirmung mehrerer Sorten** (Festkörper, Elektrolyte) · unabhängige Kanäle: Raten addieren, κ² = Σκ_i² · **Gegenbild**.
  Gewinn (gerechnet): Wäre Γ eine Summe unabhängiger Kanäle, gälte s² = s_i² + s_e² in jeder Richtung (Milchstraße 50 kpc, 0,03a₀: ν = 4,68 statt 5,46 senkrecht, frei 6,50). Die Richtungsabhängigkeit des äußeren Felds prüft, dass Γ der Betrag eines Vektors ist.

## 7  Übergangslängen als geometrische Mittel

- **Zufallsweg, Diffusionslänge, Grenzschicht** (Sternaufbau, Reaktorphysik, Strömung) · Photonenweg: Weite √(ℓ·Weg); Neutronen L² = λ_tr·λ_a/3; Blasius δ ≈ 5√(ν_kin·x/U) [korr.] · **gleich** für r*.
  r*² = r₀·Λ ist die Weite eines Zufallswegs mit Schritt r₀ = 6Δs über den Weg Λ; das gibt r* ∝ √M (Schritt ∝ M, Λ fest) und τ_D = 2Λ/c mit D = r₀c/2.
  Gewinn: Übergangslänge = Diffusionslänge, eine Lesart für alle Längen der Form √(Schritt·Weg). Bruch: Für ℓ_P² = 3Δs·λ̄ gibt es im Modell keinen Weg mit diesen Schritten: dort nur ähnlich.

## 8  Kac, Rückkehrzeit, Budget

- **Mittlere Weglänge** (Cauchy, Blanco–Mazzolo; Optik trüber Medien) · ⟨ℓ⟩ = η_d·V/S, η_d = 2, π, 4 [korr.] für d = 1, 2, 3, unabhängig von der Streulänge; in trüben Medien gemessen [G] · **gleich**.
  d = 1, am Ort ein- und austretend, bei Λ reflektierend: ⟨ℓ⟩ = 2Λ, ⟨T⟩ = 2Λ/c (Zufallsweg gerechnet: 2(Λ + h)). Gewinn: Die Masseunabhängigkeit von τ_D ist die Streuunabhängigkeit dieses Satzes; V/S steht zugleich in Δs und im Thiele-Korn. Bruch: gleichmäßiger Eintritt über den Rand vorausgesetzt.
- **Verweilzeit, Littles Gesetz, Jährlichkeit** (Verfahrenstechnik, Warteschlangen, Hydrologie) · τ̄ = V/V̇ für jede Durchmischung; N̄ = λW [korr.]; T = 1/p · **Struktur** (⟨T⟩ = μ(X)/μ(A)).
  Gewinn: a₀ = c²/ℓ₀ = c²/(2Λ) als Umsatzgleichung: Rückkehrzeit 2Λ/c = Länge der geschlossenen Halbgeraden hin und zurück durch Fluss c; Λ ist eine Länge der inneren Richtung, keine Raumlänge. [korr.] Bruch: verlangt stationären Durchfluss.
- **Drei-Abstände-Satz, Phyllotaxis** (Zahlentheorie, Botanik) · Für jede irrationale Drehung nehmen die Rückkehrzeiten in ein Intervall höchstens drei Werte an, der größte die Summe der beiden anderen (Slater); ihr Mittel ist Kac. Goldener Winkel 137,5° aus Wachstum mit Abstoßung (Douady–Couder) [G] · **gleich**.
  Gewinn: Vorhersage für Rosettenbahnen mit irrationalem Re τ_M: Die Perizentren kehren in einen Sektor nur in drei Abständen zurück; die goldene Wicklung wird durch eine Regel ohne Feinabstimmung ausgewählt. Bruch: Die Pflanze wählt über einen Hemmstoff; übertragbar ist nur die Regel.

## 9  Einrasten, Kippen, Gedächtnis, Rollen

- **Kuramoto** (Netzwerke, Biologie, Stromnetze) · K_c = 2/(πg(0)), Ordnung r = √(1 − K_c/K) · **gleich**: N_max ≈ 1 + 2W/(πk) ist K_c mit g(0) = 1/W (2/π = 0,637).
  Gewinn: Wurzelgesetz der Ordnung jenseits der Schwelle als Messgröße für die Phasenordnung vor dem Einrasten.
- **Spitzen-Katastrophe, Landau-Hysterese, Frühwarnzeichen** (Magnetismus, Ökologie, Klima) · Koerzitivfeld ∝ |t|^{3/2}, Sprung ∝ |t|^{1/2}, Engpasszeit ∝ ε^{−1/2}; AR(1)-Varianz σ²/(1 − ρ²) · **gleich** (Falten-Normalform).
  Gewinn: Die Exponenten der VO₂-Vorhersage sind universell, nicht stoffgebunden. Bruch: Galaxien erreichen q = 1 nie.
- **Kaiser-Effekt, Vorbelastung, Rückkehrpunkt-Gedächtnis** (Geologie, Bodenmechanik, Gummi, Suspensionen) · Schallemission setzt erst jenseits der früheren Höchstspannung ein (Messung der Paläospannung); die Vorbelastung erscheint als Knick der Verdichtungskurve (Casagrande); unter Rauschen bleibt nur die größte Scheramplitude [G] · **Struktur** (Preisach: nur Extremwerte).
  Gewinn: Messprinzip für „Gewicht am tiefsten Kollaps“: Nur ein Kollaps über den bisherigen Höchstwert ändert die Bindung, schwächere liegen auf der Wiederbelastungslinie. Die Knickspur hat im Boden eine wörtliche Form. Bruch: Dort tragen Mikrorisse und Körner; im Modell kein Träger, nur die Regel je Bereich.
- **Port-Hamilton, Landschaft und Fluss** (Regelungstechnik, Systembiologie) · ẋ = (J − R)∇H + Bu, J = −Jᵀ, R ≥ 0; F = −D∇U + J_ss/P · **gleich** den drei Rollen (J tangential, R transversal, u parametrisch).
  Gewinn: Passivität Ḣ ≤ yᵀu als Stabilitätsnachweis gekoppelter Schleifen; Zyklen brauchen einen Rotationsfluss, der das detaillierte Gleichgewicht bricht.

## 10  Kaleidoskop und Windung

- **Kristallographische Beschränkung, Quasikristalle** (Kristallographie) · gitterverträglich nur 2cos θ ∈ ℤ: 60°, 90°, 120°, 180°; fünfzählig nur durch Schnitt und Projektion mit irrationaler Steigung (Fibonacci-Kette = Gerade goldener Steigung auf T²) · **gleich**.
  Gewinn: Die irrationale Wicklung des Modells ist die Konstruktion eines Quasikristalls: Ordnung ohne Periode mit scharfem Punktspektrum. Zu erwarten sind scharfe Linien, kein Rauschen, im Spektrum einer quasiperiodischen Wicklung.
- **Halbzahlige Defekte** (Flüssigkristalle) · nematisch (n ≡ −n) ist Windung ½ erlaubt; bei polarer Ordnung hängt an halber Windung eine Wand mit Energie ∝ Länge · **Struktur** wie Q_el + t/3 ∈ ℤ und die Schnittlinie. Gewinn: Energie und Bewegung solcher Wände sind dort gerechnet und gemessen.

## Beobachterfehler in den Quellgebieten

Nicht übertragbar sind die Träger: Membranleck, Katalysatorstoff, Wirkstoff, Hemmstoff, Mikrorisse, zufällige Einzelereignisse, Temperatur. Übertragbar ist, was die Gebiete selbst schon als Regel erkannt haben: Photonenmasse im Supraleiter als Starrheit, effektive Masse als Dispersion, Scheinladung der Paul-Falle, nicht eindeutige Äquivalenzschicht, endogene Ereignisse bei Hawkes, Sparparadoxon.

## Ergebnis: die fünf ergiebigsten Entsprechungen

1. **Spektraldimension und kritische Verzweigung:** Der Exponent ½ heißt d_s = 1 und Kritikalität der Rückkehr. Messgröße d_s = 2·d ln a/d ln a_N; Randbedingung: Rückkehr eindimensional und driftfrei.
2. **Akkumulation (Pharmakokinetik, Reaktor 1/M):** s_beob = −ln(1 − g_bar/g_beob) linearisiert die Regel; an SPARC ohne Fehlermodell a₀ = 1,22·10⁻¹⁰, Median 1,01–1,14 im Tiefbereich; mit den Fehlern von Abstand, Neigung und Masse-zu-Licht liegt a₀ bei 1,04–1,19·10⁻¹⁰ (Υ_d = 0,5); ohne Kreispaare ist a₀ ≤ 1,13·10⁻¹⁰. [korr.] Test der Exponentialform und des äußeren Felds.
3. **Kabel, Thiele, Warburg:** dieselbe Gleichung mit versiegeltem Ende. Die Impedanzdiagnostik (45° → 90°) unterscheidet offene und bei Λ geschlossene Rückkehr; die Aufbauzeit τ_D/(2s) = t* verlangt Stationarität.
4. **Fortsetzung nach oben (Geophysik):** ϖ ist die Fortsetzungshöhe [korr.]; ln q gegen 1/r gibt r*; die Äquivalenzschicht macht die gelesene Dichte nicht eindeutig.
5. **Entropie der Rücklaufzahl:** Die Bose-Form ist die geometrische Verteilung größter Entropie; H = ln ν + χ verbindet Zusatzpotential und Schalter ohne Wärme.
