# Ramanujan-Werkzeuge für das FTM

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. y = a_N/a₀ bleibt als Abkürzung. Geändert: q = e^{−√y} als Nome gekennzeichnet (Punkt 1 des Autors, für die Punktmasse der Rückkehranteil; Kern-Nome q_Δs = Δs/r); „a₀ ist die einzige Integrationskonstante“ → falsch: die Scharkonstante c multipliziert q (s → s − ln c), eine andere Schwelle potenziert q, mit c ≠ 1 fiele Tully-Fisher weg; a₀ = c²/(2Λ) aus der Kac-Rückkehrzeit; „Lage des Δ-Nahfelds vorhergesagt“ → keine messbare Vorhersage; S-Dualität als „exakte Fassung“ von Entropie außen und Syntropie innen (Umschlag bei y = 4π²) → die Modulgruppe erhält ν nicht; vervollständigter, exakt modularer Gedächtniskern → nur unter Annahme, Q̃ ist keine Mock-Form auf SL(2,ℤ); „RAR hergeleitet“ aus dem Geschwindigkeitskern mit κ₀ = GM/r³ → Wiedergabe, Herleitung im Modell aus D·G″ = Γ·G, Geschwindigkeitskern ausgeschlossen; 2π in a₀ ≈ cH₀/2π und 3/π aus E₂ → Zahlennähe ohne Herleitung. Umbenennungen: s := ln q → ln q = −s mit s = √y; Integrations- und Mellin-Variable s → x (Kernvariable, Zeitabstand), σ (Mellin-Exponent) bzw. t = Im τ_M; Λ(σ) als L-Funktion → 𝔏(σ); ξ_k als Operator → 𝔛_k. Die Zahlen rechnen mit a₀ = 1,2·10⁻¹⁰.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/ramanujan.py` (alle Zahlen dieser Notiz).

## 0 Konventionen und Leitfrage

- $\tau_M$ liegt in der oberen Halbebene, der Nome ist $q = e^{2\pi i\tau_M}$, und $D := q\,\frac{d}{dq} = \frac{1}{2\pi i}\frac{d}{d\tau_M}$.
- **FTM-Kopplung:** $q = e^{-\sqrt y}$ mit $y = a_N/a_0$. Das ist gleichwertig mit $\tau_M = i t$ und
  $$t = \operatorname{Im}\tau_M = \frac{\sqrt y}{2\pi} = \frac{r^*}{2\pi r}.$$
  Die Spitze $q\to0$ ist das Newton-Regime, der Rand $q\to1$ das tiefe Regime. Das ist Punkt 1 des Autors: $q$ ist hier der Nome des Bahntorus und für die Punktmasse zugleich der Rückkehranteil $q=e^{-s}$, $s=\sqrt y$; der Kern-Nome des Modells heißt $q_{\Delta s}=\Delta s/r$ (Punkt 3). Die Wahl ist offen und betrifft in der stationären Regel nur die Schreibweise. [korr.]
- **Selbstdualer Punkt.** Unter der Inversion $S:\tau_M\mapsto-1/\tau_M$ gilt $t\mapsto1/t$. Der Fixpunkt $\tau_M=i$ liegt bei
  $$y = 4\pi^2 = 39{,}478,\qquad r = r^*/(2\pi).$$
  Er liegt also **nicht** bei der Schwelle $y=1$. Der galaktische Übergangsbereich $y\in[0{,}01;\,100]$ entspricht $t\in[0{,}016;\,1{,}59]$ und liegt damit überwiegend auf der Seite $q\to1$. Die Inversion $S$ erhält die Verstärkung $\nu$ nicht ($\tau_M=i/2$: $q=e^{-\pi}$, $\nu=1{,}045$; $\tau_M=2i$: $q=e^{-4\pi}$, $\nu=1{,}0000035$): Die Regel hängt am markierten Torus (nur $\tau_M\mapsto\tau_M+1$), nicht am Modulraum; $\tau_M=i$ ist für sie kein ausgezeichneter Punkt. [korr.]
- **Rollen.** tangential = Zyklen auf $T^n$, transversal = Kontraktion bzw. Dissipation, parametrisch = Bewegung des Punkts $\tau_M$ in der Basis $P$ des Bündels $\pi:M\to P$.
- **Leitfrage** jeder Abbildung: Liegt dieselbe Mathematik in derselben Rolle vor? Zahlengleichheit allein zählt nicht.

---

## 1 Serre-Ramanujan-Ableitung: $E_2$ als Zusammenhang der parametrischen Richtung

**Formeln.**
$$
D_k = q\frac{d}{dq} - \frac{k}{12}E_2,\qquad
E_2\!\Big(\frac{a\tau_M+b}{c\tau_M+d}\Big) = (c\tau_M+d)^2E_2(\tau_M) + \frac{12\,c\,(c\tau_M+d)}{2\pi i}.
$$
Den inhomogenen Term kann man umschreiben:
$$\frac{12c(c\tau_M+d)}{2\pi i} = \frac{12}{2\pi i}\,(c\tau_M+d)^2\,\partial_\tau_M\log(c\tau_M+d).$$
Das ist genau das Transformationsgesetz eines Zusammenhangs, $A\mapsto gAg^{-1}+g\,dg^{-1}$, mit der Übergangsfunktion $g=(c\tau_M+d)$.

**Numerische Prüfung** bei $\tau_M = 0{,}3+1{,}1i$:
- Transformationsgesetz von $E_2$ für vier Matrizen (bis $\operatorname{Im}(g\tau_M)=0{,}0129$): relativer Fehler ≤ $4{,}1\cdot10^{-15}$.
- Die gewöhnliche Ableitung ist **nicht** kovariant: $DE_4(g\tau_M)/((c\tau_M+d)^6DE_4(\tau_M)) = -0{,}0038+1{,}4851i$.
- Der Störterm $4c(c\tau_M+d)^5E_4/(2\pi i)$ stimmt auf $5\cdot10^{-14}$.
- Die Serre-Ableitung ist kovariant: $D_4E_4(g\tau_M)/((c\tau_M+d)^6D_4E_4) = 1{,}000000000000$.
- $D_4E_4=-E_6/3$ gilt auf $1{,}3\cdot10^{-16}$.
- $D_{12}\Delta = 0$ (per Differenzenquotient): $7\cdot10^{-10}$.

**Abbildung aufs Modell.** Ein Kern $M_k(\tau_M)$ vom Gewicht $k$ ist ein Schnitt eines Geradenbündels über dem Modulraum der Bahntori. Das Gewicht gibt an, wie sich $M_k$ bei einem Wechsel der Zyklenbasis ($SL_2(\mathbb Z)$) umskaliert. Die parametrische Wirkung $\dot p=\epsilon G$ bewegt $\tau_M$ in der Basis. Die Änderung einer Faser-Größe entlang der Basis ist nur mit einem Zusammenhang basisunabhängig definiert, und dieser Zusammenhang ist $E_2$.

**Beobachter-Lesart.** Der Term $\frac{k}{12}E_2f$ entsteht allein aus der Wahl, welcher Zyklus „radial" und welcher „azimutal" heißt. Er ist ein Koordinateneffekt wie ein Christoffel-Symbol, keine Substanz.

**Folgerungen, gerechnet.**
1. **Steigungsgesetz des Δ-Kerns.** Aus $D_{12}\Delta=0$ und $t=\sqrt y/(2\pi)$ folgt
   $$\boxed{\frac{d\ln\Delta}{d\ln r} = \sqrt{y}\;E_2(\tau_M(r))}$$
   Das Profil folgt also aus einem **Wert** und nicht aus einer Ableitung („Rechnen ohne Rechnen"). Numerisch stimmt die Formel bei $y=1;\,10;\,50$ auf alle 8 Stellen (−27,47841760; −0,48303740; +6,92656479).
2. **Maximum des Δ-Kerns entlang $r$.** Es liegt an der Nullstelle von $E_2$ auf der imaginären Achse, $t_0 = 0{,}5235217000$. Daraus folgt
   $$y_0 = (2\pi t_0)^2 = 10{,}82,\quad a_N = 1{,}30\cdot10^{-9}\ \text{m/s}^2,\quad r_{\max} = 0{,}304\,r^*.$$
   Für die Milchstraße ($r^*=8{,}35$ kpc) liegt es bei 2,54 kpc. Die Rechnung stimmt, eine messbare Vorhersage ist das nicht: Die stationäre Regel ist bei beiden Wahlen des Imaginärteils gleich; bei Punkt 1 trägt die Reihe $\sum q^n$, Δ ist bei $r^*$ nur $4{,}7\cdot10^{-8}$ davon. [korr.]
3. **Am selbstdualen Radius** $r=r^*/2\pi$ gilt wegen $E_2(i)=3/\pi$ exakt $d\ln\Delta/d\ln r = 6$, und $E_6(i)=0$. Das invariante $|\Delta|(\operatorname{Im}\tau_M)^6$ ist dort maximal ($t=1{,}000$).

---

## 2 Ramanujans Differentialsystem und Chazy: die Flussgleichung der Torusform

**Formeln.** Mit $P,Q,R=E_2,E_4,E_6$:
$$
DP = \frac{P^2-Q}{12},\qquad DQ = \frac{PQ-R}{3},\qquad DR = \frac{PR-Q^2}{2}.
$$
Für $y=\pi iE_2$ und $'=d/d\tau_M$ gilt die Chazy-Gleichung
$$y'''=2yy''-3y'^2 .$$

**Numerische Prüfung.**
- Die drei Gleichungen gelten bei $q=0{,}05;\ 0{,}3;\ 0{,}7;\ 0{,}05+0{,}2i$ auf ≤ $7{,}9\cdot10^{-15}$. Die Ableitungen sind dabei gliedweise aus den Reihen gebildet, also unabhängig vom System.
- Das Chazy-Residuum an drei Punkten ist ≤ $1{,}4\cdot10^{-15}$.

**Abbildung aufs Modell.** In der FTM-Kopplung gilt $\ln q=-s$ mit $s=\sqrt y=r^*/r$ [korr.]. Das System ist deshalb ein **autonomer polynomialer Fluss der Torusform entlang des Radius**. Die parametrische Dynamik braucht keinen äußeren Input: Jede Änderung der Form ist eine algebraische Kombination der vorhandenen Formen.

**Folgerung (eigenes Ergebnis): Newton ist ein Fixpunkt; die eine Integrationskonstante ist nicht $a_0$ [korr.].**
- Der Punkt $(P,Q,R)=(1,1,1)$, also $q=0$ bzw. das Newton-Regime, ist ein Fixpunkt des Flusses.
- Die Linearisierung $J$ hat die Eigenwerte $\{1,0,0\}$, mit $\operatorname{Rang}J=2$ und $\operatorname{Rang}J^2=1$. Zum Eigenwert 0 gehört also ein Jordan-Block.
- Der Eigenvektor zu 1 ist **exakt $(-24,240,-504)$**, das sind die Eisenstein-Koeffizienten.

Daraus folgt:
- Die einzige an der Spitze reguläre Abweichung von Newton ist $\propto q = e^{-r^*/r}$.
- Die Jordan-Richtungen erzeugen Abweichungen $\propto\ln q = -r^*/r$, also Potenzgesetze im Newton-Regime. Das Sonnensystem zeigt keine solchen Abweichungen.
- Die regulären Lösungen bilden genau die Schar $E(cq)$. Die Scharkonstante $c$ multipliziert $q$ und verschiebt $s$ um $\ln c$ ($s\to s-\ln c$). Eine andere Schwelle $a_0'$ potenziert dagegen $q$ ($q\to q^{\sqrt{a_0/a_0'}}$). Mit $c\ne1$ ginge der Rückkehranteil $cq$ im Tiefbereich gegen $c$ statt gegen 1, und Tully-Fisher fiele weg. $c$ ist also nicht die Wahl von $a_0$; die Kopplung verlangt $c=1$. [korr.]

Numerisch: Ein RK4-Fluss ab $(1,1,1)+c\,q_0(-24,240,-504)$ mit $q_0=10^{-9}$ reproduziert $E_{2,4,6}(cq)$ auf ≤ $2{,}7\cdot10^{-7}$. Der Lauf mit $c=0{,}5$ bis $q=0{,}6$ ergibt dieselben Werte wie der mit $c=1$ bis $q=0{,}3$ ($P=-17{,}26793$).

**Satz.** Der parametrische Fluss lässt bei regulärem Newton-Grenzfall genau eine Integrationskonstante zu, die Scharkonstante $c$; die FTM-Kopplung legt $c=1$ fest. Die Schwelle $a_0$ ist davon unabhängig: Im Modell folgt $a_0=c^2/(2\Lambda)$ ($c$ hier die Lichtgeschwindigkeit) aus der Kac-Rückkehrzeit auf der bei $\Lambda$ geschlossenen Halbgeraden der Rückkehr. [korr.]

Profil entlang einer Galaxie (Auszug):

| $y$ | $q$ | $E_2$ | $\Delta$ | $d\ln\Delta/d\ln r$ |
|---:|---:|---:|---:|---:|
| 100 | $4{,}5\cdot10^{-5}$ | 0,99891 | $4{,}5\cdot10^{-5}$ | +9,989 |
| 39,48 ($\tau_M=i$) | $1{,}9\cdot10^{-3}$ | 0,95493 | $1{,}79\cdot10^{-3}$ | +6,000 |
| 10,82 | 0,0373 | 0,00000 | $1{,}45\cdot10^{-2}$ (Max.) | 0 |
| 1 | 0,368 | −27,48 | $2{,}7\cdot10^{-8}$ | −27,48 |
| 0,1 | 0,729 | −356,8 | $2{,}3\cdot10^{-39}$ | −112,8 |

Die natürliche Grenze der Chazy-Lösungen auf $|q|=1$ ist aus der Literatur bekannt und hier nicht geprüft. In der Kopplung liegt das tiefe Regime $y\to0$ genau auf diesem Rand.

---

## 3 Hecke-Operatoren: die Mathematik der Skalenselbstähnlichkeit

**Formeln.** Für Gewicht $k$ ist $T_n$ ein Mittel über alle Untergitter vom Index $n$:
$$
T_nf(\tau_M) = n^{k-1}\sum_{ad=n}\sum_{b\bmod d} d^{-k}f\!\Big(\frac{a\tau_M+b}{d}\Big),\qquad T_n\Delta=\tau_R(n)\,\Delta .
$$
Dazu gelten
- $\tau_R(mn)=\tau_R(m)\tau_R(n)$ für $\gcd(m,n)=1$,
- $\tau_R(p^{k+1})=\tau_R(p)\tau_R(p^k)-p^{11}\tau_R(p^{k-1})$,
- $|\tau_R(p)|\le2p^{11/2}$ (Deligne).

**Numerische Prüfung** ($\tau_R(n)$ exakt bis $n=1000$ aus $q\prod(1-q^n)^{24}$):
- Multiplikativität an 2988 Paaren, Rekursion an 25 Fällen: keine Verletzung.
- $\max|\tau_R(p)|/(2p^{5,5}) = 0{,}9594$ für die 168 Primzahlen ≤ 1000.
- $\max|\tau_R(n)|/(d(n)n^{5,5}) = 1$ (Gleichheit nur bei $n=1$).
- Sato-Tate-Momente von $x_p=\tau_R(p)/p^{5,5}$: $\langle x^2\rangle=0{,}937$ (erwartet 1), $\langle x^4\rangle=1{,}826$ (erwartet 2).
- Gitter-Mittel direkt ausgewertet: $T_n\Delta/\Delta = -24;\ 252;\ -1472;\ 4830;\ -6048;\ -16744$ für $n=2\ldots7$, auf ≤ $10^{-10}$.
- Die Zahl der Untergitter vom Index $n$ in $\mathbb Z^2$ (Brute Force) ist $\sigma_1(n)$ für $n=1\ldots8$.

**Abbildung aufs Modell.**
- Ein Untergitter vom Index $n$ ist ein Torus, der den ursprünglichen $n$-fach überdeckt. Das ist eine Ebene mit $n$-fach feinerer Wicklung.
- $T_n$ mittelt den Kern über alle Arten, diese feinere Ebene einzubetten. $T_n$ ist also der **Skalenwechsel**.
- „Δ ist Hecke-Eigenform" heißt: Der Kern auf der feineren Ebene ist derselbe Kern mal einer Zahl. Das ist die exakte Form von „dieselbe Torus-Gleichung auf allen Ebenen, nur die Besetzung ändert sich".
- Multiplikativität heißt: Unabhängige (teilerfremde) Skalensprünge setzen sich zusammen, ohne sich zu stören.

**Folgerungen.**
1. **Beschränkte Leiter.** Mit $\tau_R(p)=2p^{11/2}\cos\theta_p$ ist der normierte Kern nach $k$ Stufen gleich $\tau_R(p^k)/p^{11k/2}=U_k(\cos\theta_p)$ (Tschebyschow 2. Art). Er bleibt durch $k+1$ beschränkt. Numerisch für $p=2$, $k=1\ldots9$ exakt; die Werte liegen zwischen −1,04 und +0,91.
   - Die Satake-Parameter $p^{11/2}e^{\pm i\theta_p}$ liegen auf einem Kreis. Der Skalenübertrag wirkt also **unitär (tangential)**: Über viele Ebenen gibt es weder exponentielles Anwachsen noch Absterben.
   - Das ist die Hecke-Fassung von „kein Informationsverlust über die Skalen".
2. **Konvergenz des Gedächtniskerns.** Aus $|\tau_R(n)|\le d(n)n^{11/2}$ folgt: Die Koeffizienten wachsen nur polynomial.
   - $\sum\tau_R(n)q^n$ konvergiert auf der ganzen Scheibe $|q|<1$ (Radius exakt 1).
   - $L(\Delta,\sigma)$ konvergiert absolut für $\operatorname{Re}\sigma>13/2$ [korr.].
   - Das Mellin-Integral des Kerns ist sogar ganz (Punkt 9).
3. **Vorzeichen.** 81 der 168 Multiplikatoren $\tau_R(p)$ mit $p\le1000$ sind negativ. Eine Bindung, die sich über Skalen aufsummiert, trägt deshalb der Eisenstein-Teil (Punkt 12a).

---

## 4 Rogers-Ramanujan-Kettenbruch: Schleife in der Schleife

**Formeln.**
$$
R(q)=\cfrac{q^{1/5}}{1+\cfrac{q}{1+\cfrac{q^2}{1+\dots}}}
= q^{1/5}\prod_{n\ge1}\frac{(1-q^{5n-1})(1-q^{5n-4})}{(1-q^{5n-2})(1-q^{5n-3})}.
$$
Weitere Beziehungen:
- $R(e^{-2\pi})=\sqrt{\tfrac{5+\sqrt5}{2}}-\varphi$,
- $R^{-5}-11-R^5=\big(\eta(\tau_M)/\eta(5\tau_M)\big)^6$,
- $R(e^{-2\pi t})=\dfrac{1-\varphi R(e^{-2\pi/t})}{\varphi+R(e^{-2\pi/t})}$.

**Numerische Prüfung.**
- Kettenbruch = Produkt auf 15 Stellen.
- $R(e^{-2\pi})=0{,}284079043840412$, die geschlossene Form stimmt auf alle Stellen.
- Mit reeller Fünfter-Wurzel-Konvention ist $R(-e^{-\pi})=-0{,}557536515835051 = -\big(\sqrt{(5-\sqrt5)/2}-1/\varphi\big)$ auf $10^{-16}$.
- Die Ikosaeder-Identität gilt auf $10^{-15}$, die S-Transformation auf alle Stellen.
- Grenzwert $q\to1$: $R\to1/\varphi=0{,}618033988750$. Die Annäherung ist exponentiell,
  $$1/\varphi-R\simeq(1+\varphi^{-2})e^{-4\pi^2/(5\varepsilon)},\qquad q=e^{-\varepsilon}.$$
  Bei $\varepsilon=1$ ist das $5{,}146\cdot10^{-4}$ gegen gemessen $5{,}145\cdot10^{-4}$.
- Der Randwert ist Nullstelle der Ikosaeder-Gleichung: $\varphi^5-11-\varphi^{-5}=2\cdot10^{-15}$.

**Abbildung aufs Modell.** Mit $x_0/x_1=K(q)=q^{1/5}/R(q)$ ist der Kettenbruch das Ebenenverhältnis der **gewichteten Fibonacci-Rekursion**
$$x_n = x_{n+1}+q^{n+1}x_{n+2}.$$
Jede Ebene $n$ speist sich aus der nächsten und, mit Gewicht $q^{n+1}$, aus der übernächsten. Das ist die verschachtelte Rückkopplung. Die geometrische Reihe $\sum q^n$ ist dagegen die einfache Schleife.

Bei voller Kopplung ($q\to1$) geht das Verhältnis gegen $\varphi$: $K(0{,}999)=1{,}61771$.

**Folgerung 1 (KAM).** $\varphi=[1;1,1,\dots]$ ist die am schlechtesten rational approximierbare Zahl. Für die besten Näherungen gilt $n^2|\varphi-m/n|\to1/\sqrt5=0{,}447214$ (gemessen 0,447214). Die verschachtelte Rückkopplung wählt im tiefen Regime also von selbst die KAM-stabilste Wicklung.

Für $\pi$ gilt dagegen:

| $m/n$ | $n^2\lvert\pi-m/n\rvert$ |
|---|---:|
| 22/7 | 0,062 |
| 355/113 | **0,0034** |
| 103993/33102 | 0,63 |

Der Teilnenner 292 macht $355/113$ zu einer Beinahe-Resonanz. Siehe Punkt 11.

**Folgerung 2: Rolle im Modell.** Ramanujans Formel
$$\cfrac{1}{1+\cfrac{aq}{1+\cfrac{aq^2}{1+\dots}}}=\frac{G(aq)}{G(a)},\qquad G(a)=\sum_n\frac{a^nq^{n^2}}{(q;q)_n},$$
ist numerisch geprüft (14 Stellen). Bei Verschachtelungstiefe 1 reduziert sich der Kettenbruch exakt auf die einfache Schleife $\nu_\text{RAR}=1/(1-q)$.

**Befund.** Die gemessene Kurve verlangt den Pol der **einfachen** Schleife bei $q=1$; die Übergangsfunktion ist die einfache Schleife. Der Kettenbruch gehört in die **Wicklungsstruktur**: Er liefert die goldene Wicklung im tiefen Regime.

---

## 5 Rogers-Ramanujan-Identitäten: Nicht-Möglichkeit regelt

**Formeln.**
$$
\sum_{n\ge0}\frac{q^{n^2}}{(q;q)_n}=\prod_{n\ge0}\frac{1}{(1-q^{5n+1})(1-q^{5n+4})},\qquad
\sum_{n\ge0}\frac{q^{n^2+n}}{(q;q)_n}=\prod_{n\ge0}\frac{1}{(1-q^{5n+2})(1-q^{5n+3})}.
$$
Kombinatorisch: Die Zahl der Partitionen mit Teilabstand ≥ 2 ist gleich der Zahl der Partitionen in Teile $\equiv\pm1\bmod5$.

**Numerische Prüfung.**
- Die Koeffizienten beider Identitäten stimmen bis $q^{60}$ exakt.
- Direkte Zählung mit Abstand ≥ 2 gegen Produktseite: gleich für $n\le25$.
- Werte bei $q=0{,}5$: 2,17266875084966, bei $q=0{,}9$: 437,66943537610 (Summe = Produkt).
- Andrews-Gordon $k=3$ ($f_1\le2$, $f_j+f_{j+1}\le2$ ⇔ Teile $\not\equiv0,\pm3\bmod7$): gleich für $n\le25$.

**Abbildung aufs Modell.** Die Teile einer Partition sind besetzte tangentiale Moden ($n$ = $n$-te Harmonische).
- „Abstand ≥ 2" ist eine **Ausschlussregel im Modenraum**: Modus $k$ und $k+1$ dürfen nicht zugleich besetzt sein.
- Das ist die Harte-Kerne-Bedingung und die eindimensionale Form von Baxters Hard-Hexagon-Modell. Dessen kritische Aktivität ist $z_c=(11+5\sqrt5)/2=11{,}09017=\varphi^5$. Das ist dieselbe Zahl wie der Randwert aus Punkt 4, $R^{-5}=\varphi^5$ bei $q\to1$.
- Die Regel ist transversal: $\Pi_\perp$ entfernt Konfigurationen, ohne die Dimension des Raums zu ändern.

**Folgerungen (quantitativ).**
1. **Harte Kerne auf $N$ Plätzen:** $F_{N+2}$ Zustände (geprüft $N\le16$).
   - Die Entropie pro Platz ist $\ln\varphi=0{,}481212$ statt $\ln2=0{,}693147$.
   - Die Ausschlussregel vernichtet also **0,2119 nats pro Freiheitsgrad (30,6 %)**.
   - Das ist die Zahl zu „Gedächtnis = verlorene Möglichkeiten".
2. **Moden-Entropie:** $\varepsilon\ln G(e^{-\varepsilon})\to\pi^2/15$ (gemessen 0,657812 bei $\varepsilon=10^{-3}$, Grenzwert 0,657974). Frei gilt $\pi^2/6$.
   - Der Ausschluss senkt die effektive Modenzahl von $c=1$ auf $c_\text{eff}=2/5$.
   - Der Faktor $q^{-1/60}=q^{-c_\text{eff}/24}$ macht $G$ modular.
   - Allgemeine Formel in 12c.

---

## 6 Mock-Theta-Funktionen: geschlossene Dynamik plus Schatten

**Formeln.**
$$f(q)=\sum_{n\ge0}\frac{q^{n^2}}{(-q;q)_n^2}=1+q-2q^2+3q^3-3q^4+3q^5-5q^6+7q^7-\dots$$
Zwegers' Vervollständigung:
$$H=F-R,\qquad F=\big(q^{-1/24}f(q),\ 2q^{1/3}\omega(q^{1/2}),\ 2q^{1/3}\omega(-q^{1/2})\big).$$
$R$ ist das nicht-holomorphe Integral des Schattens (Thetafunktionen $g_0,g_1,g_2$ vom Gewicht 3/2):
$$R(\tau_M)=2i\sqrt3\int_{-\bar\tau_M}^{i\infty}\frac{(g_1,g_0,-g_2)(z)}{\sqrt{-i(z+\tau_M)}}\,dz .$$
Gliedweise ergibt jeder Summand einen Beitrag $\frac{i}{|a|}\operatorname{erfc}\!\big(|a|\sqrt{2\pi\operatorname{Im}\tau_M}\big)e^{-\pi ia^2\tau_M}$.

$H$ transformiert als vektorwertige Form vom Gewicht 1/2:
$$H(-1/\tau_M)=\sqrt{-i\tau_M}\,M H(\tau_M),\qquad H(\tau_M+1)=\operatorname{diag-Block}(\zeta_{24}^{-1};\ \zeta_3\text{-Tausch})H(\tau_M).$$

**Numerische Prüfung.**
- Koeffizienten wie oben. Die Asymptotik $a(n)\sim(-1)^{n-1}e^{\pi\sqrt{n/6-1/144}}/(2\sqrt{n-1/24})$ ist bei $n=50;100;200;400$ auf 1,0035; 0,99931; 0,99994; 0,999999 erfüllt.
- **Zwegers bei drei Punkten:** S- und T-Gesetz für $H$ gelten auf ≤ $5\cdot10^{-16}$. Der holomorphe Teil $F$ allein verfehlt das S-Gesetz um 77–91 %.
- **Ramanujans Grenzwert** $q\to-1$ radial ($q=-e^{-\varepsilon}$): $f$ und $b=(q;q)_\infty/(-q;q)_\infty^2$ wachsen bis $\pm1{,}28\cdot10^{19}$ ($\varepsilon=0{,}01$), aber $f+b$ bleibt endlich: 3,0751; 3,4815; 3,6917; 3,8274; 3,9250; 3,9613. Die Extrapolation ergibt 3,9976, also den Grenzwert 4.
- Einfachstes Beispiel derselben Struktur: $E_2^*=E_2-3/(\pi\operatorname{Im}\tau_M)$ ist exakt modular vom Gewicht 2 (Verhältnis 1,000000000000), $E_2$ allein nicht ($-0{,}593-0{,}473i$).

**Abbildung aufs Modell (Rollen).**

| Baustein | Mathematik | Rolle |
|---|---|---|
| holomorpher Teil $F$ | $q$-Reihe aus den Zyklen, nur fast modular | tangential: geschlossene Torus-Dynamik |
| Schatten $g$ und Integral $R$ | erfc-Kerne in $\operatorname{Im}\tau_M$ = Lösungen der Wärmeleitungsgleichung, gemessen durch $\mathfrak X_k=2i(\operatorname{Im}\tau_M)^k\bar\partial$ [korr.] | transversal: Diffusion bzw. Dissipation |
| Vervollständigung $H$ | exakt modular | symmetrisches Gesamtsystem |

Keiner der beiden Teile ist allein symmetrisch, die Summe ist es exakt.

Beobachter-Lesart: Die „Nicht-Modularität" des Beobachteten (nur $F$ gemessen) ist kein Defekt des Mediums. Sie ist der weggelassene transversale Teil.

**Folgerung: Anteil des Schattens.** $|R_0/F_0|$ entlang $\tau_M=it$ (FTM-Kopplung $y=(2\pi t)^2$):

| $t$ | $y$ | $\lvert R_0/F_0\rvert$ |
|---:|---:|---:|
| 0,15 | 0,89 | 0,999 |
| 0,5 | 9,87 | 0,890 |
| 1 | 39,5 | 0,611 |
| 2 | 158 | 0,296 |
| 4 | 632 | 0,081 |
| 8 | 2527 | 0,0076 |

Der Anteil fällt asymptotisch wie $e^{-\pi t/6}/\sqrt t$. Der transversale Anteil des Gedächtnisses ist im tiefen Regime also gleich groß wie der tangentiale und verschwindet zur Spitze hin exponentiell. Auswahlregel dazu in 12b.

---

## 7 Dedekind-η, Partitionen, $1+2+3+\dots\mapsto-1/12$

**Formeln.**
$$
\frac{1}{\prod_{n\ge1}(1-q^n)}=\sum_n p(n)q^n=\frac{q^{1/24}}{\eta(\tau_M)},\qquad
\eta(-1/\tau_M)=\sqrt{-i\tau_M}\,\eta(\tau_M),\qquad
\eta(\tau_M+1)=e^{2\pi i/24}\eta(\tau_M).
$$
Zur Summe: $\sum ne^{-n\varepsilon}=1/\varepsilon^2-1/12+O(\varepsilon^2)$. Die halbe Summe gibt $\tfrac12\sum n\mapsto-1/24$.

**Numerische Prüfung.**
- $\eta(-1/\tau_M)/(\sqrt{-i\tau_M}\eta)=1{,}000000000000$ an drei Punkten.
- Der Faktor $q^{1/24}$ ist nötig: Ohne ihn ergibt das Verhältnis 1,197+0,094i; 1,157+0,419i; 0,799−0,060i.
- $\eta(\tau_M+1)/\eta(\tau_M)=e^{2\pi i/24}$, und $\eta^{24}$ ist $T$-invariant (1,000000000000).
- $\sum ne^{-n\varepsilon}-1/\varepsilon^2=-0{,}0833329$ bei $\varepsilon=0{,}01$.
- $p(0\ldots15)=1,1,2,3,5,7,11,15,22,30,42,56,77,101,135,176$.
- $\eta(i)=\Gamma(1/4)/(2\pi^{3/4})=0{,}768225422326$.

**Abbildung aufs Modell.**
- $1/\prod(1-q^n)$ zählt die Besetzungen der tangentialen Harmonischen $n=1,2,3,\dots$ eines Zyklus. Jede ist ein Oszillator.
- Die Grundverschiebung $q^{-1/24}$ ist die regularisierte Nullpunktssumme. Sie ist **erzwungen**: Nur mit ihr ist die Zählung unabhängig von der Wahl der Zyklenbasis.
- Das gehört auf die Ebene des Quanten-Vakuum-Felds ($\Delta s=0$). Der ungefaltete Hintergrund hat trotzdem eine Grundverschiebung.
- Beobachter-Regel: Die Verschiebung ist eine Eigenschaft der geschlossenen Schleife (Geometrie), kein Inhalt. Für eine Schleife der Länge $L$ (masseloses Feld, beide Umlaufrichtungen) ist $E_0=-\pi\hbar c/(6L)$, bei $L=1\,\mu$m also $-1{,}66\cdot10^{-20}$ J. Sie hängt nur von $L$ ab.

**Folgerung.** Die 24 in $\Delta=\eta^{24}$ ist die **Ordnung des Multiplikators** von η. Sie ist die kleinste Potenz, die die Grundverschiebung 1/24 ganzzahlig macht.

---

## 8 Hardy-Ramanujan: Zustandszahl und Entropie

**Formel.**
$$p(n)\sim\frac{e^{\pi\sqrt{2n/3}}}{4n\sqrt3},\qquad S(n)=\ln p(n)\simeq\pi\sqrt{2n/3}-\ln(4\sqrt3\,n).$$

**Numerische Prüfung** ($p(n)$ exakt bis $10^4$):

| $n$ | $p(n)/\text{HR}$ | $S=\ln p$ | $\pi\sqrt{2n/3}-\ln(4\sqrt3n)$ |
|---:|---:|---:|---:|
| 10 | 0,8731 | 3,738 | 3,873 |
| 100 | 0,9563 | 19,066 | 19,110 |
| 1000 | 0,9860 | 72,258 | 72,272 |
| 10000 | 0,9956 | 245,360 | 245,364 |

Die Rademacher-Reihe mit 12 Termen trifft $p(100)=190\,569\,292$ und $p(200)=3\,972\,999\,029\,388$ exakt.

**Abbildung aufs Modell.** $n$ ist die gesamte Anregung in Einheiten des Grundzyklus, $p(n)$ die Zahl der Verteilungen auf die tangentialen Harmonischen, also die Zahl der Mikrozustände. Die Entropie wächst wie $\sqrt n$, nicht wie $n$.

**Folgerung (Syntropie ↔ Entropie).** Die Cardy-Form ist
$$S=2\pi\sqrt{c\,n/6},$$
mit $c=1$ gilt $2\pi\sqrt{1/6}=\pi\sqrt{2/3}=2{,}5650997$ (geprüft). Der Entropie-Koeffizient ist durch die **Grundverschiebung** $c/24$ aus Punkt 7 vollständig festgelegt. Mechanismus: Die S-Dualität $t\leftrightarrow1/t$ verbindet in der Zustandssumme
- den geordneten Grenzfall $q\to0$ (Grundzustand, „Zentrum") mit
- dem ungeordneten Grenzfall $q\to1$ (viele Zustände, „äußere Bögen").

Für die Zustandssumme ist die Zahl außen durch die Ordnung innen bestimmt. Eine exakte Fassung der Autoraussage, dass Entropie außen und Syntropie im Zentrum ein Kreislauf sind, ist das nicht: Die Modulgruppe erhält die Verstärkung $\nu$ nicht ($\tau_M=i/2$: $q=e^{-\pi}$, $\nu=1{,}045$; $\tau_M=2i$: $q=e^{-4\pi}$, $\nu=1{,}0000035$), $\tau_M\mapsto-1/\tau_M$ bildet also keinen Zustand der Regel auf einen gleichwertigen ab, und einen Umschlagpunkt bei $y=4\pi^2$ hat die Regel nicht. [korr.]

---

## 9 Ramanujans Master-Theorem und das Gedächtnisintegral

**Formel.**
$$\int_0^\infty x^{\sigma-1}\sum_k\frac{\varphi(k)(-x)^k}{k!}\,dx=\Gamma(\sigma)\,\varphi(-\sigma).$$
Mellin-Exponent $\sigma$, Kernvariable $x$ [korr.].

**Numerische Prüfung** (Quadratur gegen Theorem):

| Kern $K(x)$ [korr.] | $\varphi(k)$ | Ergebnis | Übereinstimmung |
|---|---|---|---|
| $(1+x)^{-\alpha}$ (algebraisches Gedächtnis) | $(\alpha)_k$ | $B(\sigma,\alpha-\sigma)$ | 12 Stellen bei $\alpha=2{,}5$, $\sigma=0{,}5;\,1;\,1{,}7$ |
| $(1-e^{-x})/x$ | $1/(k+1)$ | $\Gamma(\sigma)/(1-\sigma)$ | 8–9 Stellen |
| **Rückkopplungskern** $\sum_nq^ne^{-n\lambda_m x}$ | $\mathrm{Li}_{-k}(q)$ | $\Gamma(\sigma)\lambda_m^{-\sigma}\mathrm{Li}_\sigma(q)$ | 12 Stellen |

**Anwendung aufs Kreisbahn-Integral.** Man entwickelt $\sin(\omega_\text{orb} x)$ und setzt die Momente $\int x^mK=\Gamma(m+1)\varphi(-m-1)/\lambda_m^{m+1}$ ein [korr.]:
$$
I(\omega_\text{orb})=\int_0^\infty K(x)\sin(\omega_\text{orb} x)\,dx=\frac1\lambda_m\sum_{m\ge0}(-1)^m\Big(\frac\omega_\text{orb}\lambda_m\Big)^{2m+1}\varphi(-2m-2).
$$
- Für $\varphi\equiv1$ ergibt das $\omega_\text{orb}/(\lambda_m^2+\omega_\text{orb}^2)$.
- Für den Rückkopplungskern ergibt es geschlossen
  $$\boxed{I(\omega_\text{orb})=\sum_{n\ge1}\frac{q^n\,\omega_\text{orb}}{n^2\lambda_m^2+\omega_\text{orb}^2}}$$
  Jeder Rücklauf $n$ wirkt als eigene Lorentz-Kurve mit Rate $n\lambda_m$.
- Numerisch sind Quadratur, Reihe und geschlossene Form gleich: 0,166508894018 ($q=e^{-1}$, $\omega_\text{orb}/\lambda_m=0{,}5$) bzw. 0,312304513604.

**Folgerung (eigenes Ergebnis): Die RAR als Grenzfall des Summen-über-Rückläufe-Kerns bei langem Gedächtnis (Wiedergabe, keine Herleitung) [korr.].**

[korr.] Diese Rechnung setzt einen Geschwindigkeitskern auf jeder Bahn mit $\kappa_0=GM/r^3$ voraus; das trägt nicht. Hergeleitet wird die Relation damit nicht, denn $q$ steckt schon im Kern. Ein Geschwindigkeitskern würde exzentrische Bahnen runden, alte Scheiben verlangen dafür $\lambda_m\ll0{,}01\,\Omega$. Die Rückkehr ist eine stationäre Schleife am Ort ($\tau_D=79$ Mrd. Jahre, länger als das Weltalter); das Gedächtnis ist eins (Mode ξ, Gewicht $w$, Übergaberate $\lambda_m=\omega_\text{orb}/N^*$, $N^*\approx0{,}3$–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn. Im Modell folgt die Form aus der Rückkehr $D\,G''=\Gamma\,G$ (hier $D$ Diffusion, $\Gamma=a_N/c$ Verlustrate; $q=e^{-s}$), und $\kappa_0/\lambda_\text{eff}=3c^2/(r^*v_\text{rot})$.
- Für $\lambda_m\ll\omega_\text{orb}$ gilt $\omega_\text{orb} I\to q/(1-q)=\nu_\text{RAR}-1$ (gemessen 0,5820 → 0,58198 bei $\lambda_m/\omega_\text{orb}=10^{-3}$).
- Auf der Kreisbahn (Radius $r$, Winkelgeschwindigkeit $\omega_\text{orb}$) ist die radiale Gedächtnisbeschleunigung $a_\mathsf K=\kappa_0\,r\,\omega_\text{orb}\,I(\omega_\text{orb})$. Damit folgt $a_\mathsf K\to\kappa_0\,r\,q/(1-q)$.
- Mit der natürlichen Kopplung **$\kappa_0=GM/r^3=\omega_K^2$** (Kepler-Frequenz², Einheit s⁻²) wird daraus exakt
  $$a=a_N\,\nu_\text{RAR}.$$

Unter dieser Annahme ist die angepasste Übergangsfunktion der Grenzfall unendlichen Gedächtnisses des Rückkopplungskerns, der Grenzfall „Bindung ohne Energieverlust“; als Bahnkern auf jeder Bahn ist er ausgeschlossen. [korr.]

Weitere Kennwerte:
- $K(0)=\mathrm{Li}_0(q)=q/(1-q)$.
- Das integrierte Gedächtnis ist $\int K=-\ln(1-q)/\lambda_m$ (0,458675 bei $q=e^{-1}$).
- Beim algebraischen Kern $(1+x)^{-3}$ [korr.] existiert nur das erste Moment ($\varphi(-\sigma)$ hat einen Pol bei $\sigma=3$). Es gilt $I/\omega_\text{orb}\to0{,}5$ (0,4354 bei $\omega_\text{orb}=0{,}1$, 0,4924 bei $\omega_\text{orb}=0{,}01$), aber keine Potenzreihe darüber hinaus. Ein langreichweitiges Gedächtnis erzeugt also nicht-analytische Bahnkorrekturen.

**Mellin-Momente des Δ-Kerns.**
$$\mathfrak L(\sigma)=\int_0^\infty t^{\sigma-1}\Delta(it)\,dt=(2\pi)^{-\sigma}\Gamma(\sigma)L(\Delta,\sigma).$$
Vervollständigte L-Funktion $\mathfrak L$, $t=\operatorname{Im}\tau_M$ [korr.].
- Geprüft bei $\sigma=8;10;12$ (z. B. $\mathfrak L(12)=1{,}0486273129\cdot10^{-2}$ aus Integral und Reihe).
- Weil Δ an **beiden** Spitzen verschwindet, konvergiert das Integral für **jedes** $\sigma$.
- Die Funktionalgleichung $\mathfrak L(\sigma)=\mathfrak L(12-\sigma)$ gilt auf 12 Stellen ($\sigma=0,2,4,6$). Der Zentralwert ist $\mathfrak L(6)=1{,}544879\cdot10^{-3}$. [korr.]

Der Δ-Kern ist also ein Gedächtnis mit endlichen Momenten jeder Ordnung, das mit seinem Spiegelbild übereinstimmt.

---

## 10 Lambert-Reihe: Rückkopplung über alle Harmonischen

**Formeln.**
$$
\sum_{n\ge1}\sigma_1(n)q^n=\sum_{n\ge1}\frac{n\,q^n}{1-q^n}=\frac{1-E_2}{24},
\qquad
\sum_{n\ge1}\frac{q^n}{1-q^n}=\sum_{n\ge1}d(n)q^n .
$$
Numerisch bei $q=0{,}3$: Die drei Ausdrücke der ersten Identität ergeben 0,76116395560637; $\sum q^n/(1-q^n)=0{,}56686583469627$.

**Abbildung aufs Modell.** Die Harmonische $n$ ist eine Schleife mit $q^n=e^{-\sqrt{n^2y}}$. Der Term $n=1$ von $\sum q^n/(1-q^n)$ ist exakt $\nu_\text{RAR}-1$, also die einfache Schleife. Die Lambert-Reihe ist die Rückkopplung über **alle** Harmonischen.

**Befund.** Die Übergangsfunktion ist die einfache Schleife $\nu_\text{RAR}=1/(1-q)$. $E_2$ ist das Objekt für die **Steigung** (Punkt 1) und die **Skalenleiter** (12a).

---

## 11 Ramanujans 1/π-Reihe: π als modulare Größe

**Formel.**
$$\frac1\pi=\frac{2\sqrt2}{9801}\sum_{k\ge0}\frac{(4k)!\,(1103+26390k)}{(k!)^4\,396^{4k}} .$$

**Numerische Prüfung** (110 Stellen): Jeder Term liefert ≈ 8 Stellen (7,1 · 15,2 · 23,2 · … · 63,3 Stellen nach 1…8 Termen). Der singuläre Modul dahinter:
$$e^{\pi\sqrt{58}}=24\,591\,257\,751{,}9999998222=396^4-104-1{,}78\cdot10^{-7}.$$
Außerdem gilt $E_2(i)=0{,}954929658551372=3/\pi$ exakt, also $E_2^*(i)=0$.

**Abbildung aufs Modell.** π tritt in die modulare Welt **durch die Vervollständigung** ein:
$$E_2^*=E_2-\frac{3}{\pi\operatorname{Im}\tau_M}.$$
Am quadratischen Torus $\tau_M=i$, im FTM am selbstdualen Radius $r^*/2\pi$, heben sich der tangentiale Teil $E_2$ und der transversale Teil $3/(\pi\operatorname{Im}\tau_M)$ exakt auf. Ramanujans Reihen sind Auswertungen dieser Struktur an CM-Punkten ($\tau_M=i\sqrt{58}$). π ist dort eine Konstante der Metrik bzw. Krümmung des Modulraums, keine Wicklungszahl.

**Folgerung.** π und die Wicklung sind zwei getrennte Größen:
- **π** ist die Konstante der Kreis- und Metrikstruktur (Umfang, $E_2^*$, $a_0\approx cH_0/2\pi$). [korr.] Die Verbindung des 2π in $a_0\approx cH_0/2\pi$ mit $E_2(i)=3/\pi$ ist Zahlennähe ohne Herleitung, trägt nicht: Im Modell kommt das 2π aus $a_0=c^2/(2\Lambda)$ mit $\Lambda\approx\pi c/H_0$ (nur heute).
- **Die KAM-stabile Wicklung** ist nobel bzw. golden (Punkt 4). π $=[3;7,15,1,292,\dots]$ hat die Beinahe-Resonanz 355/113 mit $n^2|\pi-m/n|=0{,}0034$, $\varphi$ dagegen 0,447.

Ein Torus mit Wicklung π wäre unter Störung deutlich weniger stabil als einer mit goldener Wicklung.

---

## 12 Eigene Formeln

### 12a Hecke-Skalenleiter

**Formel.** $E_2$ (und damit die Lambert-Rückkopplung $\ell=\sum\sigma(n)q^n$) ist Hecke-Eigenfunktion vom Gewicht 2 mit Eigenwert $\sigma_1(p)=1+p$, also der Zahl der Untergitter vom Index $p$. In FTM-Form ($q=e^{-\sqrt y}$), mit $\ell^{(p)}(y):=\sum_m\sigma(pm)e^{-m\sqrt y}$, der Rückkopplung des $p$-fach gröberen Torus:
$$
\boxed{\;\ell^{(p)}(y)=(1+p)\,\ell(y)-p\,\ell(p^2y)\;}
$$
Für den Δ-Kern gilt analog
$$\sum_m\tau_R(pm)q^m=\tau_R(p)\Delta(q)-p^{11}\Delta(q^p).$$

**Prüfung.**
- $T_pE_2/E_2=3;4;6;8$ für $p=2,3,5,7$ auf $10^{-15}$.
- Die Leiter gilt auf 12 Stellen, z. B. $y=1$, $p=2$: 3,154001886595 beidseitig.
- Die Δ-Leiter gilt auf 12 Stellen, z. B. $y=4$, $p=3$: −413,0282220449.

**Lesart.** Die Rückkopplung einer Ebene ist exakt aus der nächstfeineren Ebene ($\ell(y)$) und der $p^2$-fach näher liegenden ($\ell(p^2y)$) bestimmt. Die Ebenen sind nicht unabhängig. Das ist die rechnerische Form der Selbstähnlichkeit.

**Auswahlregel.** Eisenstein-Anteile haben positive Multiplikatoren ($\sigma_{k-1}(p)>0$), die Spitzenform Δ wechselnde (81 von 168 negativ). Eine Bindung, die sich über Skalen **aufsummiert**, kann deshalb nur vom Eisenstein-Teil getragen werden. Δ bleibt Nahfeld-Form.

### 12b Vervollständigter Gedächtniskern: tangential + transversal

[korr.] Dieser Abschnitt setzt voraus, dass der Gedächtniskern des Modells eine exakt modulare Vervollständigung hat; das trägt nicht: Die Strukturwirkung $\tilde Q$ und die Rückkehrsumme stehen zwar wie nichtholomorpher Teil und Schatten zueinander, $\tilde Q$ ist aber keine Mock-Form auf $SL(2,\mathbb Z)$, denn die Gammapolynome passen nicht, und Formen vom Gewicht 5 gibt es dort nicht. Die Selektionsregel gilt daher nur für einen angenommenen Kern dieser Form.

**Formel.**
$$
\hat K(\tau_M)=\underbrace{K_\text{hol}(\tau_M)}_{\text{tangential, }q\text{-Reihe}}\;+\;\underbrace{\int_{-\bar\tau_M}^{i\infty}\frac{g(z)}{(-i(z+\tau_M))^{k}}\,dz}_{\text{transversal, Schatten }g\text{ vom Gewicht }2-k},\qquad
\hat K\ \text{exakt modular vom Gewicht }k.
$$

**Selektionsregel (eigenes Ergebnis).** Der Rest im Newton-Regime (Spitze) hängt an der Art des Schattens:
- **Eisenstein-Schatten** (Beispiel $E_2^*$, Schatten konstant): Rest $3/(\pi\operatorname{Im}\tau_M)$, also ein **Potenzgesetz**. Die Planetenbahnen zeigen keinen solchen Rest.
- **Spitzenform-Schatten** (Beispiel $f$, Schatten Theta vom Gewicht 3/2): erfc-Faktor, Rest $\propto e^{-2\pi\kappa\operatorname{Im}\tau_M}$, **exponentiell**. Für $f$ gilt $\kappa=1/12$ (gemessen $|R_0/F_0|\propto e^{-\pi t/6}$).
- Mit der Annahme $|\nu-1|<10^{-10}$ bei Saturn (Größenordnung, im Skript als Annahme gesetzt) folgt $\kappa>0{,}031$ bei FTM-Kopplung ($f$ erfüllt das: Rest $2{,}8\cdot10^{-27}$).

**Satz (nur unter dieser Annahme) [korr.].** Ein transversaler (dissipativer) Anteil des Gedächtnisses ist mit dem Sonnensystem nur verträglich, wenn sein Schatten eine Spitzenform ist. Dann ist er im geordneten Regime exponentiell unterdrückt und im tiefen Regime gleich groß wie der tangentiale Anteil (Punkt 6).

### 12c Ausschluss-Entropie

**Formel.** Für die Andrews-Gordon-Klasse $k$ (höchstens $k-1$ Besetzungen in je zwei benachbarten Moden; $k\to\infty$: kein Ausschluss) gilt
$$
\boxed{\;c_\text{eff}=1-\frac{3}{2k+1},\qquad S(n)\simeq2\pi\sqrt{\frac{c_\text{eff}\,n}{6}}\;}
$$
Mit Meinardus-Vorfaktor $\prod_a\Gamma(a/m)/\sqrt{2\pi}\cdot A^{1/4}/(2\sqrt\pi n^{3/4})$ und $A=\pi^2c_\text{eff}/6$.

**Prüfung** (Zählung bis $n=4000$, Verhältnis Zahl/Vorhersage bei $n=500/2000/4000$):

| Regel | $c_\text{eff}$ | $n=500$ | $n=2000$ | $n=4000$ |
|---|---:|---:|---:|---:|
| frei (Bosonen) | 1 | 0,980 | 0,990 | 0,993 |
| Pauli (verschiedene Teile) | 1/2 | 0,992 | 0,996 | 0,997 |
| RR ($k=2$) | 2/5 | 0,989 | 0,995 | 0,996 |
| AG $k=3$ | 4/7 | 0,990 | 0,995 | 0,997 |
| AG $k=4$ | 2/3 | 0,991 | 0,995 | 0,997 |

**Lesart.** Die Stärke einer Ausschlussregel ist eine Zahl $c_\text{eff}$, und der Entropieverlust ist $\pi\sqrt{2n/3}\,(1-\sqrt{c_\text{eff}})$. Für RR bei $n=1000$ sind das 29,8 nats. Das ist das Maß für „Gedächtnis als Topologie verlorener Möglichkeiten". Zugleich ist $c_\text{eff}/24$ die Grundverschiebung der jeweiligen Zählung, so dass Ausschluss (transversal) und Vakuumverschiebung (Punkt 7) dieselbe Zahl sind.

### 12d Zusammenfassung der weiteren eigenen Ergebnisse

- **Parametrischer Fluss** (Punkt 2): Newton ist Fixpunkt mit Eigenwerten $\{1,0,0\}$. Die an der Spitze reguläre Integrationskonstante ist die Scharkonstante $c$ (Kopplung: $c=1$), nicht $a_0$. [korr.]
- **Steigungsgesetz** (Punkt 1): $d\ln\Delta/d\ln r=\sqrt y\,E_2$. Δ hat sein Maximum bei $y_0=10{,}82$ ($r=0{,}304\,r^*$); das ist keine messbare Vorhersage, die stationäre Regel ist bei beiden Wahlen des Imaginärteils gleich. [korr.]
- **RAR aus dem Gedächtnis** (Punkt 9): Rückkopplungskern $\sum q^ne^{-n\lambda_m x}$, langes Gedächtnis und $\kappa_0=GM/r^3$ geben $a=a_N/(1-q)$ wieder, $q$ steckt im Kern; hergeleitet ist die Form im Modell aus der Rückkehr $D\,G''=\Gamma\,G$, ein Geschwindigkeitskern auf jeder Bahn ist ausgeschlossen. [korr.]

---

## 13 Übersicht der Abbildungen

| # | Ramanujan-Struktur | FTM-Struktur | Rolle | Kernzahl |
|---|---|---|---|---|
| 1 | $E_2$ als Zusammenhang, $D_k$ | Ableitung entlang der Basis $P$ | parametrisch | $y_0=10{,}82$, $r_{\max}=0{,}304r^*$ |
| 2 | $P,Q,R$-System, Chazy | Fluss der Torusform entlang $r$ | parametrisch | Eigenwerte $\{1,0,0\}$, Vektor $(-24,240,-504)$ |
| 3 | Hecke $T_n$, Δ-Eigenform | Skalenwechsel zwischen Ebenen | Skalen | $\lvert U_k\rvert\le k+1$, RP 0,959 |
| 4 | Rogers-Ramanujan-Kettenbruch | verschachtelte Schleife, goldene Wicklung | tangential | $R\to1/\varphi$, $K(q\to1)\to\varphi$ |
| 5 | RR-Identitäten | Ausschlussregel im Modenraum | transversal | $c_\text{eff}=2/5$, −0,212 nats/FG |
| 6 | Mock-Theta + Schatten | geschlossene Dynamik + Dissipation | tangential + transversal | Zwegers $10^{-16}$, Grenzwert 4 |
| 7 | η, $q^{1/24}$ | Grundverschiebung des Vakuum-Hintergrunds | Hintergrund | $-1/24$; 24 = Multiplikatorordnung |
| 8 | Hardy-Ramanujan, Cardy | Entropie außen ↔ Ordnung innen (nur Zustandssumme; ν nicht S-invariant) [korr.] | S-Dualität | $\pi\sqrt{2/3}=2\pi\sqrt{1/6}$ |
| 9 | Master-Theorem, Mellin | Gedächtnisintegral in geschlossener Form | Gedächtnis | $\omega_\text{orb} I\to\nu_\text{RAR}-1$ |
| 10 | Lambert-Reihe, $E_2$ | Rückkopplung über alle Harmonischen | tangential | $\sum\sigma_1(n)q^n=(1-E_2)/24$ |
| 11 | 1/π-Reihe, $E_2(i)=3/\pi$ | π als Vervollständigungskonstante | Metrik | 8 Stellen/Term |
| 12 | eigene Formeln | Leiter, Selektionsregel, Ausschluss-Entropie | – | siehe 12a–d |

## 14 Vorschläge für das Gesamtmodell (warten auf Freigabe)

1. Die Übergangsfunktion ist die einfache Schleife $\nu=1/(1-q)$. [korr.] Hergeleitet ist sie im Modell aus der Rückkehr $D\,G''=\Gamma\,G$; der Rückkopplungskern mit $\kappa_0=GM/r^3$ gibt sie nur wieder ($q$ steckt im Kern), als Geschwindigkeitskern auf jeder Bahn ist er ausgeschlossen (Punkt 9).
2. $E_2$ wird als Zusammenhang der parametrischen Richtung ins Modell aufgenommen, mit Steigungsgesetz $d\ln\Delta/d\ln r=\sqrt y\,E_2$. Das Δ-Maximum bei $y_0=10{,}82$ ist keine prüfbare Vorhersage: Die stationäre Regel ist bei beiden Wahlen des Imaginärteils gleich, bei Punkt 1 trägt $\sum q^n$. [korr.]
3. Entfällt: Die reguläre Integrationskonstante des Ramanujan-Flusses ist die Scharkonstante $c$, nicht $a_0$ (sie multipliziert $q$, eine andere Schwelle potenziert $q$); $a_0=c^2/(2\Lambda)$ ($c$ hier die Lichtgeschwindigkeit) folgt aus der Kac-Rückkehrzeit (Punkt 2). [korr.]
4. Wicklung: π ist die Metrik-Konstante, die KAM-stabile Wicklung ist golden (Punkte 4, 11).
5. Dissipative Gedächtnisanteile nur mit Spitzenform-Schatten zulassen (12b), nur falls der Kern eine modulare Vervollständigung hat; $\tilde Q$ hat keine auf $SL(2,\mathbb Z)$. [korr.]
