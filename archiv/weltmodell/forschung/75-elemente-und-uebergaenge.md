# Elemente und Übergänge: was den nächsten Zustand der Materie bestimmt

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: nur Umbenennungen, keine inhaltliche Korrektur: Suszeptibilität χ → κ_m; Gitter Λ → 𝓖; Spearman-Koeffizient ρ → r_S.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/elemente_analyse.py` (Ausgabe `elemente_analyse.txt`)
> Daten: PubChem-Periodensystem, 118 Elemente (`daten/elemente/pubchem_periodensystem.csv`; Radius = Van-der-Waals-Radius, Dichte der Gase = Gasdichte); AME2020 (Huang, Wang, Kondev, Audi, Naimi 2021; `daten/elemente/mass_1.mas20.txt`), gemessen = ohne #: 2550 von 3558 Kernen. [A] Annahme, [G] aus dem Gedächtnis, nicht in den Daten.

Mit Daten: Schmelzpunkt 103, Siedepunkt 93, Ionisierungsenergie (IE) 102, Elektronenaffinität (EA) 57, Elektronegativität (EN) 95, Radius 99, Dichte fest oder flüssig 85 Elemente. He (Schmelzpunkt unter Druck), C und As (sublimieren) sind aus T_b/T_m genommen. Zufallsschutz: Spearman r_S [korr.] mit 95-%-Bootstrap, Permutations-p, Bonferroni.

## 1 Phasen

- **Das einzelne Atom legt den Schmelzpunkt nicht fest.** ln T_m gegen IE, EN, EA, Radius: |r_S| ≤ 0,20, nicht signifikant; Vorhersage aus IE, EN, r (Leave-one-out) R² = −0,12. Es trägt, was erst im Festkörper entsteht: Atomvolumen r_S = −0,56 [−0,68; −0,41] [korr.], R² = 0,45; mit der d-Füllung 0,53 (Zufall ≤ 0,01). T_b ebenso (Atomvolumen −0,47, Dichte +0,53, IE nur −0,35).
- **d-Block:** T_m über die Gruppe ist je Reihe eine Parabel, R² = 0,78 / 0,96 / 0,97 (3d / 4d / 5d), Scheitel bei Gruppe 5,6 / 6,5 / 6,9, im Mittel 6,35 (Permutation p < 2·10⁻⁴): halbe Füllung der 12 s+d-Plätze (Friedel, bekannt). Bindung am stärksten, wenn alle bindenden Zustände besetzt sind und kein lockernder.
- **Sprünge zum Nachbarn:** C→N 3823→63 K, Ne→Na 25→371 K, Au→Hg 1337→234 K, Si→P 1687→317 K, Ar→K 84→337 K. Sie liegen, wo ein durchgehendes Bindungsnetz (Metall, kovalentes Gitter) in abgeschlossene Einheiten (Atom, [G] Moleküle N₂, P₄) übergeht. Der Übergang zählt die Bindungen zwischen den Einheiten, nicht die in ihnen. Waagrechter und senkrechter Nachbar sind für T_m, T_b gleich ähnlich (p > 0,1).
- **Fenster T_b/T_m:** 90 Elemente, Median 1,82, von 1,023 (Xe) bis 8,18 (Ga). Nach Block: Edelgase 1,04, Halogene 1,32, Erdalkali 1,58, Nichtmetalle 1,70, Übergangsmetalle 1,77, Lanthanoide 1,90, Halbmetalle 2,08, Actinoide 2,62, Alkali 3,12, Post-Übergangsmetalle 3,37 (Kruskal-Wallis p = 3·10⁻⁸). Es hängt an der Bindung des äußeren Elektrons: IE gegen T_b/T_m r_S = −0,50, in Nachbardifferenzen −0,48 [korr.]. Je loser das Elektron, desto breiter das flüssige Fenster.
- [G] Richards (ΔS_m ≈ k je Atom) und Trouton (ΔS_b ≈ 10,5 k) geben ΔH_b/ΔH_m ≈ 10,5·T_b/T_m ≈ 19: Schmelzen kostet nur etwa 5 % der Bindung.

## 2 Elektronen

- **Schwelle am Schalenschluss:** Die sechs größten Abfälle der IE von Z nach Z+1 (von 101) sind genau die sechs Schritte Edelgas → Alkalimetall (24,59 → 5,39 eV bei He → Li usw.); Zufall 1/C(101,6) = 8·10⁻¹⁰. Knicke innerhalb der Periode (Gruppe 2 > 13, 15 > 16) in den Perioden 2–4, nicht in 5 und 6.
- **Sättigung:** Für Z ≤ 56 fehlt die EA (kein oder kaum gebundenes Anion) bei He, Be, N, Ne, Mg, Ar, Ca, Mn, Zn, Kr, Sr, Cd, Xe, Ba. Alle 14 haben nur volle oder halbvolle (p, d) Valenzunterschalen; solche gibt es 19 (hypergeometrisch p = 2·10⁻⁹). Die übrigen fünf (P, As, Sb, Tc, Pd) binden trotzdem (0,55–1,07 eV). Die Schale ist eine Zählregel, 2(2l+1) Plätze je Unterschale; ist sie voll, hat das nächste Elektron keinen Zustand.
- **Bindung als Längenverhältnis:** IE/(m_e c²) = r_e/(2r_IE), r_IE = e²/(8πε₀·IE), beim Wasserstoff r_IE = a_B. Über 99 Elemente liegt IE·r_vdW im Median bei 1,02·e²/(4πε₀) (16–84 %: 0,89–1,31), also IE/(m_e c²) ≈ r_e/r_vdW. Die Steigung ln IE gegen ln r ist aber −1,39 [−1,69; −1,12], nicht −1: Größenordnung, kein Gesetz.

## 3 Kern

- **Maximum:** größte Bindung je Nukleon ⁶²Ni 8794,56, ⁵⁸Fe 8792,25, ⁵⁶Fe 8790,36 keV. Die kleinste Masse je Nukleon, also die kleinste Faltung Δs/A, hat ⁵⁶Fe (vor ⁶⁰Ni, ⁶²Ni), weil ⁶²Ni mehr der schwereren Neutronen trägt. Jeder freiwillige Kernübergang (A bleibt) senkt die Faltung je Nukleon; ihr tiefster Punkt ist ⁵⁶Fe.
- **Tröpfchenmodell** (bekannt) an 2484 gemessenen Kernen (A ≥ 16): a_V = 15,41, a_S = 16,86, a_C = 0,695, a_A = 22,50, a_P = 12,03 MeV, rms 3,3 MeV. Hergeleitet: Bei festem Z/A ist d(B/A)/dA = 0 genau bei E_C = E_S/2, also Spaltbarkeit x = E_C/(2E_S) = ¼. Gerechnet: x = 0,240 am Maximum des Modells (A = 58), 0,249 (⁵⁶Fe), 0,261 (⁶²Ni).
- **Kippstelle:** Verformungsenergie (2/5)·E_S·α²·(1 − x); die Antwort auf eine Quadrupolstörung wächst wie 1/(1 − x), bei x = 1 verschwindet die Rückstellung (Bohr-Wheeler). Größtes gemessenes x = 0,927 (²⁶⁹Ds), 1/(1 − x) = 14; so schwere Kerne bestehen nur durch Schalen.
- **Magische Zahlen ohne Vorgabe:** Schalenlücke S_2n(N) − S_2n(N+2) minus Tröpfchengang, Median über alle Z. Die vier größten N-Lücken sind 82, 50, 126, 28 (3,2 / 3,2 / 3,0 / 2,4 MeV, dann 152 mit 0,6), 20 auf Rang 7 (p = 2·10⁻⁴); bei Z sind die drei größten 50, 28, 82, 20 auf Rang 7 (p = 5·10⁻³). 8 wird nicht gefunden, weil das Tröpfchenmodell leichte Kerne schlecht trägt.
- **Tropflinie:** S_n < 0 bei 17 von 2393, S_p < 0 bei 76 von 2343 gemessenen Kernen (4H bis 16B; Protonenseite bis ²¹⁹Np). Der Kern gibt ein Nukleon ab, sobald die Separationsenergie das Vorzeichen wechselt.
- **Mehrfachstabilität aus der Paarung:** Lokale Massenminima längs der Isobaren (stabil gegen einfachen β-Zerfall): ungerades A 127 von 127 Ketten genau eines; gerades A 74 Ketten zwei, 5 drei, 48 eines. Gerade-gerade und ungerade-ungerade Kerne liegen auf zwei Parabeln im Abstand 2a_P/√A; zwischen den Minima führt nur der doppelte β-Zerfall.

## 4 Trichter je Element und Kern

| | H | Fe | U |
|---|---|---|---|
| Δs = Gm/(3c²) | 4,1·10⁻⁵⁵ m | 2,3·10⁻⁵³ m | 9,8·10⁻⁵³ m |
| Δs/r_vdW | 3,5·10⁻⁴⁵ | 1,2·10⁻⁴³ | 4,1·10⁻⁴³ |
| Takt am Kern 3Δs/R, R = 1,2 A^{1/3} fm [A] | 1,0·10⁻³⁹ | 1,5·10⁻³⁸ | 3,9·10⁻³⁸ |
| α_G = 3Δs/λ̄ = (m/m_Pl)² | 5,9·10⁻³⁹ | 1,8·10⁻³⁵ | 3,3·10⁻³⁴ |
| Kernbindung B/A/(m_u c²) | 3,0·10⁻³ | 9,4·10⁻³ | 8,2·10⁻³ |
| Hüllenbindung IE/(m_e c²) | 2,7·10⁻⁵ | 1,5·10⁻⁵ | 1,2·10⁻⁵ |
| Phase kT_m/(m c²) | 1,3·10⁻¹² | 3,0·10⁻¹² | 5,5·10⁻¹³ |

- **Die eigene Faltung wirkt nicht:** Gravitation zweier Nachbaratome G m²/(2r_vdW) ist 1,5·10⁻³³ bis 7·10⁻³⁰ von kT_m und 4·10⁻³⁷ bis 2·10⁻³² von der IE; im Kern ist der Takt 1,6·10⁻³⁶ der Bindung (⁵⁶Fe). Die Rückkehrregel ruht: Jedes Labor liegt im Erdfeld, q = e^{−√(g⊕/a₀)} = e^{−2,9·10⁵}. Die Kernbindung ändert die Faltung eines Atoms um höchstens 0,94 %.
- **Die kollektive Faltung wirkt ab einer Körpergröße:** Für eine homogene Kugel ist Zentraldruck mal Atomvolumen P_c·V_A = ½·m_A c²·δ mit dem Takt δ = 3Δs/R. Der Trichter verschiebt die Phasengrenzen seines Stoffs, sobald
  ```
  δ ≳ δ_krit = 2kT/(m_A c²) = 2(v_th/c)²
  ```
  Derselbe Quotient entscheidet, ob ein Körper sein Gas hält (Jeans, bekannt). Über 85 feste und flüssige Elemente ist δ_krit(T_m) = 2·10⁻¹³ bis 6·10⁻¹¹, das sind Kugeln von im Median 7·10²⁰ kg und 270 km Radius (16–84 %: 170–590 km); Eisen 4·10²¹ kg, 495 km. Die Erde (δ = 7,0·10⁻¹⁰) liegt 116-fach über δ_krit(Fe). Größenordnung; [G] Kleinkörper werden ab einigen hundert Kilometern rund.

## 5 Regeln des Modells an den Daten

- **Rückkopplungssumme und Kippstelle q = 1:** Im Kern ist q = x: Rückstellung ∝ 1 − x, Antwort ∝ 1/(1 − x), Kippstelle x = 1; gemessen 0,26 (⁶²Ni) bis 0,93. In der Hülle ist q = (4π/3)α′/V_A die Rückkopplung der lokalen Felder, ε = (1 + 2q)/(1 − q), Metall ab q = 1 (Herzfeld). Mit geschätzter Polarisierbarkeit [A] liegt die trennende Schwelle bei 0,25 [0,20; 0,56] (α′ aus IE) bzw. 0,72 [0,27; 1,08] (α′ = r_vdW³); IE allein trennt Metalle und Nichtmetalle besser (0,99 richtig). Mit diesen Daten nicht entscheidbar. [G] Magnetismus hat dieselbe Form (Stoner, κ_m = κ_m,0/(1 − I·N(E_F))) [korr.], keine Daten.
- **Schmelzen und Sieden sind keine Kippstellen.** Als Übergänge erster Art liegen T_m und T_b zwischen zwei Stellen q = 1 (Grenzen von Unterkühlung und Überhitzung), am Gleichstand der freien Energien im Hysteresefenster. Hysterese und kritisches Verlangsamen stehen nicht in den Daten (keine Unterkühlung, keine Halbwertszeiten, keine Schwingungsfrequenzen).
- **Kopplungsfenster:** fest = verriegelt (Gitterphasen gebunden), flüssig = gekoppelt, nicht verriegelt, Gas = zerfasert. Schmelzen ist der Eintritt ins Fenster von der verriegelten Seite, Sieden der Austritt nach unten. Setzte eine Kopplung beide Ränder, wäre T_b/T_m fest; gemessen 1,02–8,18. Am schmalsten ist das Fenster bei reiner Paarkopplung (Edelgase 1,02–1,10); bei loser Hüllenbindung überlebt der Zusammenhalt den Verlust der Ordnung.
- **Bistabilität:** Die Paarung macht gerade Isobaren mehrfach stabil (79 von 127), ungerade nie. Übergang über eine Schwelle höherer Masse.
- **Bose-Leiter, Planck-Form:** Keine Wärmekapazität, kein Spektrum, keine Debye-Temperatur in den Daten; nicht prüfbar. [G] Lindemann: ⟨u²⟩ = (ħ/2mω)·coth(ħω/2kT) = (ħ/2mω)(1 + 2n_B), Schmelzen im Rayleigh-Jeans-Bereich.
- **Geometrische Mittel:** In der Hülle ist λ̄_e² = r_e·a_B exakt (r_e = αλ̄_e = α²a_B): Die Compton-Länge steht zwischen klassischem Radius und Bahnradius wie ℓ_P² = 3Δs·λ̄ zwischen Faltung und Quantenlänge. Der Schwellenradius r* eines Atoms (31 fm bei H bis 0,52 pm bei Og) liegt zwischen Kern und Hülle, ohne Wirkung.
- **Drehordnungen:** Ein Kristall ist ein Zustand auf dem Torus ℝ³/𝓖 [korr.]. Eine Gitterdrehung hat ganzzahlige Spur, 2cos(2π/n) ∈ ℤ, also n ∈ {1, 2, 3, 4, 6}: derselbe Satz wie für die Drehungen flacher Tori. Kristallstrukturen stehen nicht in den Daten.

## 6 Musterscan

- 9 Größen (IE, EN, EA, r, Atomvolumen, ln T_m, ln T_b, T_b/T_m, größtes B/A je Element), 36 Paare, als Werte und als Nachbardifferenzen (gegen den glatten Gang mit Z), Bonferroni ×72: 12 bzw. 6 Paare signifikant, alle bekannte Chemie (IE–EN +0,90, T_m–T_b +0,89, IE–r −0,73, Atomvolumen–T_m −0,56, IE–T_b/T_m −0,50).
- **Kern und Hülle hängen nur über dieselbe Zählregel zusammen.** An den kernmagischen Z = 8, 20, 28, 50, 82 zeigen T_m, T_b und IE nichts (p = 0,68–0,76). IE und größtes B/A laufen in Nachbardifferenzen mit (r_S = +0,38) [korr.]: Alle Edelgase haben gerades Z, und gerades Z bindet im Kern stärker. Ohne die Umgebung der Edelgase bleibt r_S = +0,18 (p = 0,09) [korr.]. Beides ist die Regel „zwei Zustände je Bahn“, kein Einfluss des Kerns auf die Hülle. Die Schalenzahlen beider Ebenen sind verschieden (Hülle 2, 10, 18, 36, 54, 86; Kern 2, 8, 20, 28, 50, 82, 126).

## Beobachterfehler

Kein Übergang bringt einen Stoff hinzu. Schmelzen hebt die Identifikation x ~ x + 𝓖 [korr.] auf (das Gitter als Torus) und gibt die Verschiebung frei, die Schersteifigkeit verschwindet, bei etwa 5 % Bindungsverlust; Sieden hebt die Kopplung auf; Ionisation gibt einem Elektron freie Bewegung; ein Kernübergang besetzt um, hin zu kleinerer Faltung je Nukleon. Die latente Wärme als „Wärmestoff“ war genau dieser Fehler [G]. Schalen und magische Zahlen sind Zählregeln, keine Gegenstände.

## Ergebnis

- **Trägt:** Auf jeder Ebene tritt der nächste Zustand ein, wenn eine Energiebilanz das Vorzeichen wechselt (T = ΔE/ΔS, S_n = 0, Q_β = 0), begrenzt durch eine Zählregel (Schalen, Paare, halbe Füllung) und ein Verhältnis zweier Terme (x = E_C/2E_S, Maximum bei ¼, Kippstelle bei 1). Bei den Phasen entscheidet das Netz, nicht das Atom. Trichter bestimmen den Zustand erst kollektiv, ab δ ≈ 2kT/(m_A c²), also ab Körpern von etwa 10²¹ kg.
- **Trägt nicht:** die eigene Faltung von Atom und Kern (29 bis 36 Größenordnungen zu klein), die Rückkehrregel (q = 0 im Labor), eine Wirkung des Kerns auf Hülle oder Phasen. Nicht entscheidbar mit diesen Daten: Herzfeld-Kippstelle, Hysterese und Verlangsamen, Bose-Form, Drehordnungen je Element.
