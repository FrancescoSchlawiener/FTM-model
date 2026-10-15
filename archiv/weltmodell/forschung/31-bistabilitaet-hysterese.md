# Bistabilität, Hysterese und die Geburt des Torus

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Nachlauf in Galaxien (2.5, N4, K2, Offen) mit τ_R = r/v_c, τ_rel = 359 Mio. Jahre bei 20 kpc und 1,9 Mrd. Jahre bei 50 kpc, Außenpunkte neben der RAR → trägt nicht: Die Rückkehr ist eine stationäre Schleife am Ort, die Rückkopplung ist stationär, höchstens ν·r/c ≈ 1 Mio. Jahre; kein messbarer Nachlauf in Galaxien, Endantwort 0,54–0,59·δM/M bleibt. §3.5, N6, K5, Offen: gleichverteilte Schwellen (w = 1 − θ_min/θ_vir), Einfallsschalen w = 0,48–0,80, δg/g = Qw/ν (K5: 0,21–0,35 bei 100 kpc, 0,07–0,12 bei 300 kpc, 0,02–0,04 bei 1 Mpc) → im Modell w = θ₊²/(θ₊² + θ_c²), also w = 0 für jedes θ ≤ 0; die Mode geht dort mit λ_m = ω_orb/N* über; um eine Galaxie wie die Milchstraße folgt die Linse bis 100 kpc der Relation, bei 200 kpc liegt sie mit N* = 2,4 um 3–7 %, bei 300 kpc um 17–44 % darunter, mit N* = 0,3 nicht messbar; Vorzeichen negativ. [G] an T_C(VO₂) und 2β. Umbenannt (ohne Marke): w (Übergangsbreite) → Δ_T; s(T) (Schaltfunktion) → h(T); Γ = Bℓ (Lorentz-Kopplung) → G_y.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/bistabilitaet.py`

**Worum es geht.** Ein System, das sich über eine Schleife selbst beeinflusst, hat drei Möglichkeiten. Es kann einen Arbeitspunkt halten, es kann kippen und es kann umlaufen. Diese Notiz rechnet alle drei an einem realen Bauteil durch und überträgt die Ergebnisse auf das Modell. Fünf Fragen:

1. Wie entsteht ein Torus im Labor, Schritt für Schritt von Fixpunkt über Grenzzyklus zu $T^2$ und weiter zu Einrasten oder Chaos? Welche Lyapunov-Signatur hat jede Stufe, und wie verteilen sich die drei Rollen?
2. Warum liegt die Kippstelle genau dort, wo die Rückkopplungssumme $\sum q^n=1/(1-q)$ divergiert? Was folgt daraus für Galaxien, deren $q=e^{-\sqrt y}$ nie 1 erreicht?
3. Wie speichert Hysterese Vergangenheit, ohne sie zu speichern? Wie verhält sich das Preisach-Modell zum Gedächtniskern, und ist die Schließungsregel $w(\theta_u)$ ein Hysteron?
4. Wie teilt ein Magnetfeld den thermoelektrischen Transport in einen tangentialen und einen transversalen Anteil?
5. Welche Laborversuche prüfen das mit Zahlen, und was folgt neu für das Modell?

Die Beobachter-Regel gilt durchgehend. Gedächtnis ist hier fortwirkende Einschränkung, also der gegenwärtige Zustand. Gespeichert wird nichts. Wirkungen sind Terme, Schwellen und Projektionen, keine Substanzen.

---

## 1  Torus-Geburt im Labor

### 1.1 Das Bauteil

Ein VO₂-Element sitzt auf einer beidseitig eingespannten Mikrobrücke. Die Brücke führt den Strom und liegt in einem Magnetfeld $B$. VO₂ hat bei $T_C\approx340$ K [G] einen Metall-Isolator-Übergang (Morin 1959): Der Widerstand fällt über wenige Kelvin um mehr als eine Größenordnung. Unter Stromfluss ergibt sich daraus ein negativer differentieller Widerstand.

Vier Zustandsgrößen: Spannung $V$ über dem Element (parallel liegt die Kapazität $C$), Temperatur $T$, Auslenkung $x$ und Geschwindigkeit $\dot x$ der Brücke.

$$
\begin{aligned}
C\dot V&=\frac{V_s-V}{R_L}-I,\qquad I=\frac{V-G_y\dot x}{R(T,x)},\\
C_{\rm th}\dot T&=I^2R-G\,(T-T_0),\\
\ddot x&=-\frac{\omega_m}{Q_m}\dot x-\omega_m^2\,y\Big(1+\frac{y^2}{x_D^2}\Big),\qquad y=x-D_{\rm tr}\big(1-h(T)\big)-\frac{G_y}{k}I,\\
R(T,x)&=R_m\Big(\frac{R_i}{R_m}\Big)^{h(T-\kappa x)},\qquad h(T)=\frac{1}{1+e^{(T-T_C)/\Delta_T}} .
\end{aligned}
$$

Die Kopplungen im Einzelnen:
- **Joule und Kühlung:** $I^2R$ heizt, $G(T-T_0)$ kühlt.
- **Wärme ändert die elektrische Regel:** $R(T)$ fällt im Übergang steil ab.
- **Übergangsdehnung:** Der strukturelle Übergang biegt die Brücke um $D_{\rm tr}$ (Bimorph-Wirkung).
- **Dehnung verschiebt den Übergang:** $T_C\to T_C+\kappa x$.
- **Lorentzkraft und Bewegungs-EMK:** $G_yI$ bzw. $G_y\dot x$ mit $G_y=B\ell$.
- **Geometrische Nichtlinearität:** Die Brücke versteift sich mit der Amplitude (Duffing-Länge $x_D$).

| Größe | Wert | Herkunft |
|---|---|---|
| $T_C$ | 340 K [G] | Literatur (Morin 1959) |
| Übergangsbreite $\Delta_T$, $R_i$, $R_m$ | 1,5 K; 20 kΩ; 0,4 kΩ | Annahme (glatter Übergang ohne eigene Hysterese) |
| $C_{\rm th}$, $G$, $\tau_{\rm th}$ | 5 pJ/K; 5 µW/K; 1,0 µs | Annahme (Mikroelement) |
| $R_L$, $C$ | 10 kΩ; 20 pF | Annahme |
| $f_m$, Güte $Q_m$ | 0,50 MHz; 1000 (Vakuum) | Annahme |
| $D_{\rm tr}$, $\kappa$, $x_D$ | 0,20 µm; −0,10 K/µm; 0,30 µm | Annahme; das Vorzeichen von $\kappa$ hängt von der Einbaulage ab |
| $G_y/k$ pro Tesla | 0,063 µm/(mA·T), mit $\ell=100$ µm und $k=1{,}58$ N/m | Annahme |

Die eigene Hysterese von VO₂ ist bewusst weggelassen. Jede Hysterese, die unten auftaucht, entsteht aus der Rückkopplung.

### 1.2 Die Verzweigungsfolge

Parameter ist die Vorspannung $V_s$, zusätzlich das Feld $B$. Die Lyapunov-Spektren sind nach Benettin et al. (1980) gerechnet (RK4, $dt=0{,}01$ µs, 2500 µs Einschwingen, 3000 µs Mittelung). Ein Exponent mit $|\lambda|<0{,}01$/µs gilt als null.

**Fixpunkt und Hopf.** Am Fixpunkt ist das elektrothermische Eigenwertpaar für $V_s>V_H$ gedämpft. Bei

$$
V_H=3{,}5482\ \text{V},\qquad f_H=3{,}256\ \text{MHz}
$$

geht sein Realteil durch null (Arbeitspunkt $T^*=340{,}18$ K, $R=2{,}489$ kΩ, $I=0{,}284$ mA). Die Brücke ist dort mit rund $-0{,}00046$/µs gedämpft, die reine mechanische Dämpfung beträgt $-\omega_m/2Q=-0{,}00157$/µs.

**Grenzzyklus und Neimark-Sacker.** Unterhalb von $V_H$ läuft die elektrothermische Schwingung. Ob die Brücke mitschwingt, entscheidet ihre Floquet-Rate auf diesem Zyklus. Gemessen wird sie am Abklingen einer Anfangsauslenkung:

| $V_s$ | 3,52 V | 3,50 V | 3,48 V |
|---|---|---|---|
| Floquet-Rate der Brücke | −0,00034/µs | −0,00021/µs | −0,00005/µs |

Die Rate geht bei $V_{NS}\approx3{,}472$ V ($B=1$ T) durch null. Dort wird die zweite Frequenz geboren.

**Die ganze Folge** (bei $B=1$ T, außer wo angegeben):

| $V_s$ | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ | $\lambda_4$ [1/µs] | $\Delta x$ [µm] | $f_x$, $f_V$ [MHz] | Zustand |
|---|---|---|---|---|---|---|---|
| 3,60 V | −0,0016 | −0,0025 | −0,4823 | −0,4827 | 0,030 (abklingend) | – | Fixpunkt |
| 3,47 V, $B=0$ | +0,0005 | −0,0004 | −0,0013 | −1,5111 | 0,086 (abklingend) | 0,5086; 2,9511 | Grenzzyklus |
| 3,47 V | +0,0009 | −0,0012 | −0,0009 | −1,5076 | 0,130 | 0,5148; 2,9514; Verhältnis 0,1744 | Torus $T^2$ |
| 3,46 V | +0,0015 | −0,0003 | −0,0016 | −1,7110 | 0,263 | 0,5248; 2,9121; Verhältnis 0,1802 | Torus $T^2$ |
| 3,44 V | +0,0088 | −1,0916 | −1,0950 | −2,8848 | 4,04 | 2,7031 = 2,7031 | eingerastet 1:1 |
| 3,20 V | +0,0011 | −0,4253 | −0,4238 | −8,8293 | 2,72 | 1,8256 = 1,8256 | eingerastet 1:1 |
| 2,87 V | +0,0007 | −0,0328 | −0,0313 | −5,1076 | 0,99 | 0,7745 = 0,7745 | eingerastet 1:1 |
| 2,80 V | +0,0004 | −0,0029 | −0,0024 | −4,1881 | 0,40 | 0,5508; 0,4406; Verhältnis 1,2500 | eingerastet 5:4 |
| 2,70 V | −0,0033 | −0,0040 | −0,8437 | −7,0720 | 0,002 | – | Fixpunkt |

**Das Feld als Schalter.** Bei festem $V_s=3{,}47$ V:

| $B$ | 0 | 1 T | 2 T | 4 T |
|---|---|---|---|---|
| Zustand | Grenzzyklus | $T^2$ | $T^2$ (Brücke 0,289 µm) | eingerastet 1:1 (Brücke 4,42 µm) |

**Chaos** tritt bei einer steiferen Brücke auf ($f_m=2$ MHz, $x_D=0{,}1$ µm, $\kappa=-0{,}05$ K/µm): $\lambda_1=+0{,}334$/µs bei $V_s=2{,}90$ V und $\lambda_1=+0{,}944$/µs bei 2,80 V, jeweils mit negativer Summe der Exponenten. Der Weg über den Torus ins Chaos entspricht dem Szenario von Ruelle, Takens und Newhouse (Newhouse et al. 1978).

**Lyapunov-Signaturen** (Vorzeichen von $\lambda_1\ge\lambda_2\ge\lambda_3\ge\lambda_4$):

| Fixpunkt | Grenzzyklus | Torus $T^2$ | eingerastet | Chaos |
|---|---|---|---|---|
| $(-,-,-,-)$ | $(0,-,-,-)$ | $(0,0,-,-)$ | $(0,-,-,-)$, rationales $f_x/f_V$ | $(+,0,-,-)$ |

Die Zahl der verschwindenden Exponenten ist die Dimension des Torus. Die Summe aller Exponenten ist in jeder Stufe negativ: Das Phasenvolumen schrumpft, der Attraktor wird von Dissipation getragen.

### 1.3 Die drei Rollen je Stufe

| Stufe | tangential (Umlauf) | transversal (Dämpfung) | parametrisch |
|---|---|---|---|
| Fixpunkt | keine Richtung | alle vier Richtungen kontrahieren: Kühlung $G/C_{\rm th}$, Umladung über $R\,C$, Reibung $\omega_m/2Q_m$ | $V_s$, $B$, $T_0$ legen den Arbeitspunkt fest |
| Hopf → Grenzzyklus | eine Richtung: die Phase der elektrothermischen Schwingung | drei; die Amplitude sättigt über die Form von $R(T)$ | $V_s$ schiebt den Arbeitspunkt in den steilen Teil von $R(T)$ |
| Neimark-Sacker → $T^2$ | zwei Phasen, Rotationszahl $f_x/f_V$ irrational | zwei | $V_s$ und $B$ schieben die Floquet-Rate der Brücke durch null |
| Einrasten | eine Phase; die Differenzphase wird kontrahiert | drei, eine davon aus der Phasenkopplung | $B$ stellt die Kopplungsstärke ein, $V_s$ das Frequenzverhältnis |
| Chaos | eine Phase | eine Richtung dehnt ($\lambda_1>0$); die Sättigung faltet zurück | Steifigkeit und Nichtlinearität der Brücke |

Jede Verzweigung ist ein **Rollenwechsel**:
- Hopf macht eine transversale Richtung tangential.
- Neimark-Sacker macht eine zweite transversale Richtung tangential.
- Einrasten macht eine tangentiale Richtung wieder transversal. Die Kontraktion der Differenzphase stammt aus dem dissipativen Teil der Kopplung.
- Chaos macht eine transversale Richtung dehnend.

Der parametrische Anteil verschiebt dabei nur die Vorzeichen der Exponenten und erzeugt selbst keine Richtung.

### 1.4 Befunde

1. **Die elektrothermische Schleife pulsiert.** Der Übergang zum Umlauf braucht zwei Raten gleicher Ordnung: die Umladerate des Kreises und die elektrothermische Rate im steilen Teil von $R(T)$. Letztere liegt weit über $1/\tau_{\rm th}$, entsprechend liegt $f_H=3{,}256$ MHz weit über $1/(2\pi\tau_{\rm th})$. Die Kapazität teilt die Schleife nach Frequenzen auf:
   - Auf kurzen Zeiten hält sie die Spannung fest. Dort koppelt die Wärme mit, denn steigendes $T$ senkt $R$ und hebt $V^2/R$.
   - Auf langen Zeiten ist das Element stromgesteuert. Dort koppelt sie gegen, denn $I^2R$ fällt mit $R$.

   Das ergibt die verzögerte Gegenkopplung, die den Hopf-Punkt erzeugt.
2. **Die zweite Frequenz ist eine eigene Mode, keine erzwungene Antwort.** Eine lineare, gedämpfte Brücke, die vom Zyklus angetrieben wird, schwingt nur mit dessen Frequenz. $T^2$ entsteht erst, wenn die Schleife aus Dehnung, Übergang und Wärme die Brücke im Mittel über den Zyklus entdämpft.
3. **Das Magnetfeld ist parametrisch.** Es liefert keine eigene Richtung, entscheidet aber, ob $T^2$ existiert. Bei 3,47 V macht 1 T aus dem Grenzzyklus einen Torus, 4 T rastet ihn ein.
4. **Einrasten erfolgt sprunghaft.** Beim Einrasten springt die Brückenamplitude (Spitze-Spitze) von 0,26 auf 4,04 µm. Die Duffing-Versteifung zieht die mechanische Frequenz auf die elektrothermische. Die Amplitude wächst also nicht stetig an, sondern springt.
5. Am unteren Ende des Wegs (2,80 V) liegt ein 5:4-Fenster, dann ein Fixpunkt. Ohne Brücke ist der Kreis bei $R_L=10$ kΩ stationär bistabil zwischen 2,7374 und 2,7509 V ($R_L$ liegt knapp unter $R_c=10{,}967$ kΩ, §2.4).

---

## 2  Kippstelle = Divergenz der Rückkopplungssumme

### 2.1 Allgemein

Ein Arbeitspunkt, der sich über eine Schleife selbst erzeugt, erfüllt $z^*=F(z^*;p)$. Die stationäre Schleifenverstärkung ist $q=\partial_zF$. Eine kleine Änderung des Parameters läuft beliebig oft durch die Schleife:

$$
\delta z=F_p\,\delta p+q\,\delta z\quad\Rightarrow\quad
\boxed{\;\delta z=\frac{F_p\,\delta p}{1-q}=F_p\,\delta p\sum_{n\ge0}q^n\;}
$$

Jeder Summand ist ein Umlauf. Nach dem Satz über implizite Funktionen setzt sich der Ast $z^*(p)$ fort, solange $q\ne1$ ist. Bei $q=1$ gilt $dp/dz=0$: Der Ast faltet sich um. Das ist der Sattel-Knoten. Zwischen zwei Falten liegen drei Äste. Der mittlere hat $q>1$ und ist instabil, die äußeren haben $q<1$. Dort liegen Bistabilität und Hysterese.

Für die Dynamik $\tau_0\dot z=-z+F(z;p)$ ist die Relaxationsrate

$$
\lambda=-\frac{1-q}{\tau_0},\qquad \tau_{\rm rel}=\frac{\tau_0}{1-q}=\nu\,\tau_0 .
$$

**Dieselbe Zahl $\nu=1/(1-q)$ ist statische Verstärkung und dynamische Verlangsamung.** Die Kippstelle ist der Punkt, an dem beide divergieren. Bei mehreren Größen tritt $\det(I-G)$ an die Stelle von $1-q$, mit der Matrix $G$ der Schleifenverstärkungen; die Falte liegt bei $\det(I-G)=0$.

### 2.2 Spitze und Breite der Hysterese

Dort, wo zwei Falten zusammenlaufen, gilt zusätzlich $\partial_z q=0$ (Spitze, cusp). In der Umgebung reduziert sich die Stationärbedingung auf die Normalform

$$
x^3-a\,x-b=0,\qquad 1-q\propto3x^2-a .
$$

- **Falten** bei $3x_f^2=a$, also $b_f=\mp2(a/3)^{3/2}$.
- **Breite** des Hysteresefensters:
$$
\boxed{\;\Delta b=4\Big(\frac a3\Big)^{3/2}\propto a^{3/2}\;}
$$
- **Sprunghöhe** $\Delta x=3\sqrt{a/3}=\sqrt{3a}\propto a^{1/2}$.

Zahlen der Normalform: $a=0{,}1$ ergibt $\Delta b=0{,}02434$ und $\Delta x=0{,}54772$. $a=1$ ergibt 0,76980 und 1,73205, $a=4$ ergibt 6,15840 und 3,46410.

### 2.3 Kritische Verlangsamung

Entwickelt man am stabilen Ast vor der Falte, $b=b_f-\delta$, dann ist

$$
\boxed{\;1-q=2\,(3a)^{1/4}\,\delta^{1/2}\;}
$$

Die Relaxationsrate verschwindet also wie die Wurzel des Abstands, und die Relaxationszeit divergiert wie $\delta^{-1/2}$. Direkt an der Spitze ($a=0$) gilt $1-q=3\,b^{2/3}$. Numerisch bei $a=1$: $1-q=0{,}026437$ gegen 0,026321 aus der Formel bei $\delta=10^{-4}$, und 0,084377 gegen 0,083236 bei $10^{-3}$.

### 2.4 Nachrechnung am VO₂-Element

Das Element ohne Kapazität, mit Lastwiderstand: $G(T-T_0)=P(T)=V_s^2R/(R+R_L)^2$, $q=P'(T)/G$.

**Identität** (isolierender Ast, $R_L=0{,}8\,R_c$):

| Abstand zur Falte | $q$ | $dT^*/dV_s$ numerisch | $(\partial_VP/G)/(1-q)$ | Teilsummen $N=0,1,5,50,500$ |
|---|---|---|---|---|
| 30 mV | 0,3083 | 36,805 K/V | 36,806 K/V | 25,5; 33,3; 36,8; 36,8; 36,8 |
| 3 mV | 0,7286 | 96,654 K/V | 96,670 K/V | 26,2; 45,4; 82,2; 96,7; 96,7 |
| 0,3 mV | 0,9095 | 292,941 K/V | 292,908 K/V | 26,5; 50,6; 127,1; 290,6; 292,9 |

Je näher die Falte, desto mehr Umläufe braucht die Summe.

**Spitze.** Sie liegt bei $R_c=10{,}96694$ kΩ. Das ist auf alle angegebenen Stellen der größte negative differentielle Widerstand des Elements, $\max(-dV/dI)$. Bistabil ist das Element also genau dann, wenn die Lastgerade flacher liegt als der steilste NDR-Abschnitt. Bei $R_L=0{,}8R_c$ liegen die Falten bei 2,62683 V und 2,57738 V, mit $q=1{,}00001$ bzw. 0,99977 (Gitterauflösung).

**Breite und Sprung** gegen den Abstand $\varepsilon=1-R_L/R_c$:

| $\varepsilon$ | 0,001 | 0,00316 | 0,01 | 0,0316 |
|---|---|---|---|---|
| $\Delta V_s$ | 0,0156 mV | 0,0881 mV | 0,497 mV | 2,83 mV |
| Temperatursprung | 0,273 K | 0,481 K | 0,844 K | 1,472 K |

Die angepassten Exponenten sind **1,504** für die Breite und **0,488** für den Sprung. Die Normalform verlangt 3/2 und 1/2.

**Kritische Verlangsamung** vor der Falte bei 2,62683 V. Die lineare Rate $G(1-q)/C_{\rm th}$ und die simulierte Abklingrate stimmen auf 1 % überein:

| $V_{s,f}-V_s$ | 100 mV | 10 mV | 1 mV | 0,3 mV | 0,03 mV |
|---|---|---|---|---|---|
| $q$ | 0,0749 | 0,5403 | 0,8379 | 0,9095 | 0,9710 |
| $\lambda$ (Simulation) | 0,930/µs | 0,462/µs | 0,163/µs | 0,0910/µs | 0,0293/µs |

Der Exponent für Abstände ≤ 1 mV ist **0,490**. Die Relaxationszeit wächst von der Größenordnung $\tau_{\rm th}=1{,}0$ µs auf 11,0 µs bei 0,3 mV und 34,1 µs bei 0,03 mV Abstand. Das ist das Frühwarnsignal vor dem Kippen (Scheffer et al. 2009).

### 2.5 Galaxien: nahe an der Kippstelle, ohne sie zu erreichen

Im Modell gilt $g=g_N/(1-q)$ mit $q=e^{-\sqrt y}$ und $y=g_N/a_0$. Für jedes $y>0$ ist $q<1$. Daraus folgt:

1. **Keine Falte.** Es gilt
$$
\frac{d\ln g}{d\ln g_N}=1-\frac{q\sqrt y}{2(1-q)}\in\big(\tfrac12,\,1\big),
$$
numerisch 0,5025 bis 0,9998. $g(g_N)$ ist streng monoton. Die Radial-Acceleration-Relation hat genau einen Ast, keine Bistabilität und keine Hysterese.
2. **Der Rand nähert sich der Kippstelle nur asymptotisch.** Für kleine $y$ gilt $1-q\to\sqrt y$ (0,00995 gegen 0,01 bei $y=10^{-4}$). Die Verstärkung $\nu$ wächst nach außen ohne Grenze, und mit ihr die absolute Empfindlichkeit $dg/dg_N=\nu\,d\ln g/d\ln g_N\to1/(2\sqrt y)$. Die relative Antwort strebt dagegen gegen 1/2.
3. **Kritische Verlangsamung ohne Übergang.** Braucht eine Runde der Rückkopplung die Zeit $\tau_R$, dann folgt aus $\tau_R\dot g=-g+g_N+q\,g$
$$
\tau_{\rm rel}=\nu\,\tau_R .
$$
Zwei Grenzfälle: Die Signallaufzeit $\tau_R=r/c$ oder die dynamische Zeit $\tau_R=r/v_c$. [korr.] Das Modell legt $\tau_R$ fest: Die Rückkehr ist eine stationäre Schleife am Ort ($\tau_D=79$ Mrd. Jahre, länger als das Weltalter), kein Vorgang längs der Bahn; die Rückkopplung ist stationär, die Verlangsamung höchstens $\nu\,r/c\approx1$ Mio. Jahre (50 kpc). Die dynamische Zeit $r/v_c$ trägt nicht; ein Geschwindigkeitskern würde exzentrische Bahnen runden, alte Scheiben verlangen dafür $\lambda_m\ll0{,}01\,\Omega$. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit $1/\lambda_m$ gibt es nur, wo die Mode ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche).

**Milchstraße** ($M_b=6\cdot10^{10}\,M_\odot$, Hernquist mit $a=3$ kpc; beides Annahmen, die Größenordnung aus dem Gedächtnis):

| $r$ | $y$ | $q$ | $\nu$ | $d\ln g/d\ln g_N$ | $dg/dg_N$ | $v_c$ | $\tau_{\rm rel}$ für $\tau_R=r/v_c$ | $\tau_{\rm rel}$ für $\tau_R=r/c$ |
|---|---|---|---|---|---|---|---|---|
| 10 kpc | 0,4125 | 0,5261 | 2,110 | 0,6435 | 1,358 | 179,5 km/s | 115 Mio. J. | 0,07 Mio. J. |
| 20 kpc | 0,1318 | 0,6956 | 3,285 | 0,5853 | 1,923 | 179,0 km/s | 359 Mio. J. | 0,21 Mio. J. |
| 30 kpc | 0,0640 | 0,7765 | 4,474 | 0,5606 | 2,508 | 178,4 km/s | 736 Mio. J. | 0,44 Mio. J. |
| 40 kpc | 0,0377 | 0,8235 | 5,666 | 0,5470 | 3,099 | 177,9 km/s | 1246 Mio. J. | 0,74 Mio. J. |
| 50 kpc | 0,0248 | 0,8543 | 6,861 | 0,5383 | 3,694 | 177,5 km/s | 1889 Mio. J. | 1,12 Mio. J. |

Bei flacher Rotationskurve gilt $\nu\propto r$ und $t_{\rm dyn}\propto r$, also $\tau_{\rm rel}\propto r^2$. Die vorletzte Spalte gilt nur unter der Annahme $\tau_R=r/v_c$, die nicht trägt; im Modell ist die letzte Spalte die Obergrenze. [korr.]

**Prüfbare Aussage.** Ändert sich die baryonische Masse um $\delta M/M$ (Gaszufuhr, Ausfluss nach einem Sternentstehungsausbruch, Verschmelzung), dann folgt die äußere Beschleunigung mit

$$
\frac{\delta g}{g}(t)=\frac{d\ln g}{d\ln g_N}\,\frac{\delta M}{M}\,\big(1-e^{-t/\tau_{\rm rel}}\big).
$$

Für 1 % Massenänderung beträgt die Endantwort 0,59 % bei 20 kpc und 0,54 % bei 50 kpc. Wäre $\tau_R$ die dynamische Zeit, dann wäre nach einer dynamischen Zeit erst ein Teil erreicht; es fehlten noch 74 % (20 kpc) bzw. 86 % (50 kpc). [korr.]
- Nur unter dieser Annahme lägen Außenpunkte von Galaxien mit jüngerer Massenänderung ($t<\tau_{\rm rel}$) gegen die RAR verschoben, in Richtung des alten Zustands, und die Verschiebung wüchse nach außen wie $r^2$. [korr.]
- Ist $\tau_R=r/c$, dann ist die Verzögerung mit höchstens 1,1 Mio. Jahren unmessbar klein.

Im Modell ist die Rückkopplung stationär (Verlangsamung höchstens $\nu\,r/c$): Die Außenbeschleunigung folgt einer Massenänderung ohne messbare Verzögerung, mit der Endantwort 0,54–0,59·δM/M. [korr.]

---

## 3  Hysterese als Gedächtnis ohne Speicher

### 3.1 Das Preisach-Modell

Ein Hysteron $\hat\gamma_{\alpha\beta}$ mit $\alpha\ge\beta$ schaltet auf $+1$, wenn der Eingang $u$ den Wert $\alpha$ erreicht, und auf $-1$ bei $u\le\beta$. Dazwischen behält es seinen Wert. Das Preisach-Modell (Preisach 1935) ist die gewichtete Überlagerung

$$
f(t)=\iint_{\alpha\ge\beta}\mu(\alpha,\beta)\,\hat\gamma_{\alpha\beta}[u](t)\,d\alpha\,d\beta .
$$

Im Dreieck $\alpha\ge\beta$ trennt eine Treppenlinie die Hysteronen auf $+1$ von denen auf $-1$. Ihre Ecken sind die **dominanten Extremwerte** der Vergangenheit, eine alternierende Folge $M_1>M_2>\dots$ und $m_1<m_2<\dots$ Zwei Eigenschaften folgen daraus:
- **Auslöschung (wiping-out):** Übersteigt $u$ ein früheres Maximum $M_k$, dann werden alle späteren Paare $(M_j,m_j)$ mit $j\ge k$ gelöscht, ebenso gespiegelt für Minima.
- **Rückkehrpunkt-Gedächtnis:** Schließt der Eingang eine Nebenschleife und kehrt zum Anfangswert zurück, dann kehrt auch der Zustand exakt zurück.

Nach Mayergoyz (1986) ist eine Hysterese genau dann durch Preisach darstellbar, wenn sie Auslöschung und Kongruenz der Nebenschleifen zeigt.

### 3.2 Preisach gegen Gedächtniskern

| | Gedächtniskern $y=\int K(t-t')u(t')\,dt'$ | Preisach |
|---|---|---|
| Linearität | linear, Superposition gilt | nichtlinear |
| Zeit | zeitinvariant, ratenabhängig (außer $K=\delta$) | ratenunabhängig: invariant unter jeder monotonen Umparametrisierung der Zeit |
| was wirkt fort | ein gewichtetes Mittel über die ganze Vergangenheit, mit der Zeit schwindend | nur die nicht ausgelöschten Extremwerte, nicht schwindend |
| Ableitung | überall linearisierbar | an Umkehrpunkten verschiedene Steigung auf- und abwärts, also keine Linearisierung |

**Wann geht das eine ins andere über?** Man gebe jedem Hysteron eine Barriere $E_b(u)$ und thermisch aktivierte Raten $k=\nu_0e^{-E_b/kT}$. Dann gilt:
- Hysteronen mit $E_b>kT\ln(\nu_0t_{\rm obs})$ sind **eingefroren**. Sie bilden den Preisach-Teil.
- Die übrigen **relaxieren**. Für kleine Amplituden ist jedes ein Debye-Relaxator, und zusammen ergeben sie einen Prony-Kern
$$
K(\Delta t)=\sum_i2\,w_i\,s_i\,e^{-K_i\Delta t}.
$$
Er ist eine Summe von Termen der Form $\nu e^{-h\Delta t}$.
- Linear ist die Antwort nur für Amplituden unterhalb der **Glättungsbreite** $kT/|\partial_uE_b|$. Über diese Breite verschmiert die Temperatur die Schaltkante.

Der Kern ist also der thermisch geglättete Kleinsignal-Grenzfall eines Ensembles bistabiler Elemente. Preisach ist sein Grenzfall bei eingefrorener Aktivierung und großer Amplitude.

### 3.3 Nachrechnung

7260 Hysteronen auf einem Gitter, $\mu\propto\exp[-c^2/(2\cdot0{,}3^2)-(h-0{,}25)^2/(2\cdot0{,}15^2)]$, mit Mitte $c$ und halber Breite $h$.

- **Auslöschung.** Geschichte A mit den Extremwerten (−1; 0,8; −0,6; 0,5; −0,3; 0,2) und Geschichte B mit (−1; 0,8; −0,6; 0,4; −0,2; 0,35; −0,1; 0,5; −0,3; 0,2) haben denselben Stapel (−1; 0,8; −0,6; 0,5; −0,3; 0,2). Alle 7260 Hysteronen stehen gleich, der Ausgang ist in beiden Fällen 0,0987617341. Bei gleicher Zukunft (30 zufällige Umkehrpunkte) sind die Ausgänge exakt gleich (Abweichung 0).
- **Rückkehrpunkt.** $f(0{,}5)=0{,}5468901531$ vor und nach der Nebenschleife 0,5 → −0,2 → 0,5. Auch der Zustand ist identisch.
- **Gedächtnis ohne Speicher.** Ein Zufallsweg mit 200 000 Schritten hat 100 054 Umkehrpunkte. Am Ende stehen 6 Werte im Stapel, unterwegs höchstens 17. Der Zustand aus diesen 6 Werten ist identisch mit dem Zustand aus allen 200 000 Eingangswerten.
- **Ratenunabhängigkeit.** Geschichte A mit 40 oder 400 Schritten je Ast ergibt in Preisach beide Male 0,098762. Ein Kern mit $\tau_0=50$ liefert dagegen 0,007383 bzw. 0,138177. Für die Geschichten A und B liefert der Kern 0,007383 bzw. 0,038318, obwohl deren Preisach-Zustand gleich ist.

**Thermisch aktivierte Hysteronen** (576 Elemente, $\nu_0=10^6$, $E_b=40\,h(1\mp(u-c)/h)^2$):

| $kT$ | eingefroren | Remanenz $f(0\,|\,{+1})-f(0\,|\,{-1})$ | Rückkehrpunkt-Fehler | Ratenabhängigkeit (10× langsamer) |
|---|---|---|---|---|
| 0,2 | 83,3 % | 0,6385 | $3{,}0\cdot10^{-3}$ | $9{,}6\cdot10^{-3}$ |
| 1,0 | 16,7 % | 0,0355 | $1{,}1\cdot10^{-3}$ | $2{,}8\cdot10^{-3}$ |
| 3,0 | 0 | 0 | 0 | 0 |

Kalt ist das Ensemble nahezu ein Preisach-System: Es hat große Remanenz, und Rückkehrpunkt und Ratenunabhängigkeit gelten bis auf $10^{-3}$ bis $10^{-2}$. Heiß ist das Gedächtnis verschwunden, und der Ausgang ist eine Funktion des momentanen Eingangs.

Abweichung zwischen der vollen Antwort und dem Prony-Kern:

| | Glättungsbreite $kT/2E_0$ | Amplitude 0,003 | 0,03 | 0,3 |
|---|---|---|---|---|
| $kT=3$, Gleichgewicht | 0,0375 | 0,0020 | 0,0965 | 0,1403 |
| $kT=1$, nach Vorgeschichte | 0,0125 | 0,0553 | 0,9159 | 0,9067 |

Der Kern trägt, solange die Amplitude unter der Glättungsbreite bleibt. Darüber schalten Hysteronen, und die Antwort wird Preisach-artig.

### 3.4 Die Schließungsregel $w(\theta_u)$ als Hysteron

Das Modell verlangt ein Projektionsgewicht $w$: $w=1$ im mitatmenden Hintergrund, $w\to0$ in gebundenen Systemen. Dabei ist $\theta_u=\nabla_\mu u^\mu$ die Expansionsrate des Zeitfelds. Ein Bereich, der beim Kollaps vom Atmen entkoppelt, bleibt entkoppelt, auch wenn $\theta_u$ lokal später wieder schwankt. Das ist genau das Verhalten eines Hysterons, und die Regel lautet als Preisach-Operator:

$$
\boxed{\;w[\theta_u](\mathbf x,t)=\iint_{\alpha\ge\beta}\mu(\alpha,\beta)\,\hat\gamma_{\alpha\beta}\big[\theta_u(\mathbf x,\cdot)\big](t)\,d\alpha\,d\beta,\qquad\hat\gamma\in\{0,1\}\;}
$$

- $\hat\gamma=1$ heißt mitatmend.
- **Entkopplung** bei $\theta_u\le\beta$ (Kollaps).
- **Wiedereinkopplung** erst bei $\theta_u\ge\alpha$. Liegt $\alpha$ beim Hubble-Fluss, also $\theta_u/3H\approx1$, dann geschieht das für einen gebundenen Bereich praktisch nie.

Damit wird $w$ eine Funktion des **tiefsten** Punkts der Geschichte allein:

$$
w=F(\theta_{\min})=\int\mu(\beta)\,\Theta(\beta<\theta_{\min})\,d\beta .
$$

Der Bereich trägt seine Vergangenheit als gegenwärtige Einschränkung, nicht als Aufzeichnung. Spätere, flachere Einbrüche ändern nichts (Auslöschung). Nur ein tieferer Kollaps senkt $w$ weiter.

**Kugelkollaps** ($R=A(1-\cos\eta)$, $t=B(\eta-\sin\eta)$, Hintergrund Einstein-de Sitter; Gunn & Gott 1972):

| $\eta/\pi$ | 0,5 | 1,0 (Umkehr) | 1,1 | 1,25 | 1,4 | 1,5 (Virialisierung) |
|---|---|---|---|---|---|---|
| $R/R_{\rm ta}$ | 0,500 | 1,000 | 0,976 | 0,854 | 0,655 | 0,500 |
| $\theta_u/3H_{\rm bg}$ | +0,856 | 0 | −0,458 | −1,687 | −4,454 | −8,569 |

**Bei der Virialisierung ist $\theta_{\min}/3H_{\rm bg}=-8{,}5686$, unabhängig von der Masse.** Ein vollständiger Kugelkollaps hinterlässt also bei jeder Masse denselben Zustand.

### 3.5 Größe des Effekts: die Bindung hängt von der Kollapsgeschichte ab

[korr.] Diese Abschätzung setzt gleichverteilte Schwellen und die Interpolation $g=\nu g_N+wQg_N$ voraus; das trägt nicht: Im Modell gilt $w=\theta_+^2/(\theta_+^2+\theta_c^2)$, also $w=0$ für jedes $\theta\le0$ (innerhalb der Umkehr); die Mode geht dort mit $\lambda_m=\omega_{\rm orb}/N^*$ über. Mit dieser Regel gerechnet (Rechnung rest_mode_isolierte_galaxie, zwei Wege): Um eine Galaxie wie die Milchstraße ist die Mode bis 100 kpc ganz übergegangen, die Linse folgt der Relation; bei 200 kpc liegt sie mit $N^*=2{,}4$ um 3–7 %, bei 300 kpc um 17–44 % unter der Relation, mit $N^*=0{,}3$ nicht messbar. Das Vorzeichen ist negativ: Der Rest der Mode senkt die Linse. Die Zahlen unten gelten nur unter den Annahmen der Notiz.

Annahmen für die Abschätzung:
- Die Schwellen $\beta$ sind gleichverteilt in $[\theta_{\rm vir},0]$. Dann gilt $w=1-\theta_{\min}/\theta_{\rm vir}$.
- Der projizierte Fluss trägt linear interpoliert bei: $g=\nu g_N+w\,Q\,g_N$. Dabei ist $Q=5{,}364$ sein Gewicht im mitatmenden Hintergrund (das Verhältnis 5,36 aus dem CMB).

Damit gilt

$$
\frac{\delta g}{g}=\frac{Q\,w}{\nu}.
$$

- **Virialisierte Innenbereiche:** $w=0$, reiner Operator, RAR exakt. Die gemessene intrinsische Streuung der RAR (≲ 0,057 dex, aus dem Gedächtnis; Lelli et al. 2017) begrenzt dort $w<0{,}026$ bei $\nu=1$, $<0{,}052$ bei $\nu=2$, $<0{,}078$ bei $\nu=3$ und $<0{,}157$ bei $\nu=6$. Mit $\theta_{\min}=\theta_{\rm vir}$ ist das von selbst erfüllt.
- **Auslöschung in Zahlen:** Ein früher, ungestörter Kollaps mit späteren kleinen Schwankungen und derselbe Kollaps mit einer späten Verschmelzung bis $\theta_u=-4$ ergeben beide das Stapel-Minimum −8,569 und damit $w=0$.
- **Abgebrochener Kollaps:** Hält ein Gezeitenfeld den Kollaps bei $\theta_{\min}=-3$ an, dann ist unter der Annahme $w=0{,}650$; bei $\nu=3$ folgte $\delta g/g=1{,}162$. Im Modell ist dort $w=0$. [korr.]
- **Einfallsbereich:** Schalen zwischen $R_{\rm ta}/2$ und $R_{\rm ta}$ hätten unter der Annahme $0<w<1$; im Modell ist dort $w=0$. [korr.] Für ein System mit $M_b=6\cdot10^{10}M_\odot$ (Punktmasse), nur unter der Annahme [korr.]:

| $r$ | $\nu$ | $\delta g/g$ bei $w=0{,}803$ ($R=0{,}854R_{\rm ta}$) | $\delta g/g$ bei $w=0{,}480$ ($R=0{,}655R_{\rm ta}$) |
|---|---|---|---|
| 100 kpc | 12,48 | 0,345 | 0,206 |
| 300 kpc | 36,43 | 0,118 | 0,071 |
| 1000 kpc | 120,27 | 0,036 | 0,021 |

**Trennsignal.** Der Effekt ist keine Funktion des heutigen Zustands, sondern des tiefsten Kollapses:
1. Material, das schon einmal voll kollabiert war, ist entkoppelt und bleibt es. Gezeitenzwerggalaxien, die aus der Scheibe einer Muttergalaxie entstehen, erben deren $w=0$ und liegen auf der RAR.
   - In einem Bild mit Teilchen-Dunkelmaterie enthalten sie kaum Dunkelmaterie und sollten Newtonsch sein.
   - Gemessen wurde bei den Gezeitenzwergen um NGC 5291 fehlende Masse (Bournaud et al. 2007), und ihre Rotation ließ sich mit der Beschleunigungsrelation beschreiben (Gentile et al. 2007). Beides ist aus dem Gedächtnis zitiert.
2. Unvollständig kollabierte Bereiche (Gezeitenabbruch, Einfallszone) trügen unter der Annahme eine Zusatzbeschleunigung $Q\,w/\nu$. Im Modell ist dort $w=0$; eine Abweichung trägt nur der noch nicht übergegangene Rest der Mode, sie hängt an der Zeit seit der Umkehr und an $N^*$ und senkt die Linse (um die Milchstraße bei 200 kpc 3–7 %, bei 300 kpc 17–44 % mit $N^*=2{,}4$). [korr.]
3. Spätere, schwächere Verschmelzungen hinterlassen in $w$ keine Spur.

---

## 4  Magnetfeld teilt den Transport in tangential und transversal

### 4.1 Zerlegung nach Onsager-Casimir

Die Flüsse sind Ladungsstrom $\mathbf J$ und Wärmestrom $\mathbf J_Q$ in der Ebene senkrecht zu $\mathbf B$:

$$
\mathbf J=\sigma\mathbf E-\alpha\nabla T,\qquad \mathbf J_Q=T\tilde\alpha\,\mathbf E-\bar\kappa\nabla T .
$$

Mit den Kräften $X=(\mathbf E/T,\ \nabla(1/T))$ gilt $J=LX$ mit

$$
L(B)=\begin{pmatrix}T\sigma(B)&T^2\alpha(B)\\T^2\tilde\alpha(B)&T^2\bar\kappa(B)\end{pmatrix},\qquad
\dot s=X^{\!\top}LX .
$$

Die Onsager-Casimir-Relation $L_{ij}(B)=L_{ji}(-B)$ (Onsager 1931, Casimir 1945) lautet in Matrixform $L^{\!\top}(B)=L(-B)$. Daraus folgt die Zerlegung

$$
L_s=\tfrac12\big(L+L^{\!\top}\big)=\tfrac12\big(L(B)+L(-B)\big),\qquad
L_a=\tfrac12\big(L-L^{\!\top}\big)=\tfrac12\big(L(B)-L(-B)\big).
$$

- $L_s$ ist **gerade in $B$**. Er enthält die Leitfähigkeit mit magnetischem Widerstand, Seebeck und Peltier sowie die Wärmeleitung.
- $L_a$ ist **ungerade in $B$**. Er enthält Hall ($\sigma_{xy}$), Nernst-Ettingshausen ($\alpha_{xy}$) und Righi-Leduc ($\kappa_{xy}$).

**Der antisymmetrische Teil erzeugt keine Entropie.** Für jede Kraft gilt $X^{\!\top}L_aX=0$, also $\dot s=X^{\!\top}L_sX\ge0$. Der Fluss $J_a=L_aX$ steht senkrecht auf der Kraft. Er dreht den Fluss um die Kraft, ohne Entropie zu erzeugen, so wie die Lorentzkraft die Geschwindigkeit dreht, ohne Arbeit zu leisten.

**Zuordnung in GENERIC** (Grmela & Öttinger 1997; $\dot z=\mathcal L\nabla E+\mathcal M\nabla S$, $\mathcal L=-\mathcal L^{\!\top}$, $\mathcal M=\mathcal M^{\!\top}\ge0$):
- $L_a$ gehört zur Rolle **tangential**, zu $\mathcal L$: antisymmetrisch, ungerade in $B$, entropieneutral.
- $L_s$ gehört zur Rolle **transversal**, zu $\mathcal M$: symmetrisch, positiv, entropieerzeugend.

Das Magnetfeld selbst ist **parametrisch**. Es erzeugt keinen Fluss, es legt fest, welcher Anteil des Transports umläuft und welcher dissipiert.

### 4.2 Bismut bei 300 K

Zweibandmodell (Elektronen und Löcher, kompensiert):
- $n=p=2{,}7\cdot10^{24}$ m⁻³, $\mu_e=1{,}2$ und $\mu_h=0{,}6$ m²/Vs. Das ist die Größenordnung für Bismut bei Raumtemperatur, aus dem Gedächtnis.
- Band-Seebeck $-130$ und $+60$ µV/K sowie Gitterwärmeleitung 2 W/mK: Annahmen.

Ergebnis ohne Feld: $\rho=1{,}284\cdot10^{-6}$ Ωm und $S=-66{,}7$ µV/K. Beide passen zu den Tabellenwerten für Bismut (etwa $1{,}3\cdot10^{-6}$ Ωm und −60 bis −70 µV/K, aus dem Gedächtnis).

| $B$ | 0,1 T | 1 T | 3 T |
|---|---|---|---|
| $\rho_{xx}$ [Ωm] | $1{,}294\cdot10^{-6}$ | $2{,}209\cdot10^{-6}$ | $9{,}607\cdot10^{-6}$ |
| Magnetowiderstand $\Delta\rho/\rho$ | 0,007 | 0,720 | 6,480 |
| Hall-Winkel $\tan\theta_H$ | −0,060 | −0,349 | −0,241 |
| Nernst $E_y/\partial_xT$ | −7,6 µV/K | −76,0 µV/K | −228 µV/K |
| Righi-Leduc $\kappa_{yx}/\kappa_{xx}$ | 0,035 | 0,151 | 0,038 |
| $\kappa_{xx}$ ($J=0$) [W/mK] | 9,51 | 6,83 | 4,59 |
| $\lvert J_a\rvert/\lvert J_s\rvert$ bei Kraft längs $x$ | 0,060 | 0,350 | 0,252 |
| $\max\lvert X\!\cdot\!L_aX\rvert/\max X\!\cdot\!L_sX$ | $4\cdot10^{-18}$ | $3\cdot10^{-17}$ | $2\cdot10^{-17}$ |

- Der Nernst-Koeffizient des Modells ist $-76$ µV/(K·T), in diesem Modell unabhängig von $B$. Für Bismut wird bei Raumtemperatur die Größenordnung $10^{-4}$ V/(K·T) berichtet, bei tiefen Temperaturen ein um Größenordnungen größerer Wert (Behnia et al. 2007; beides aus dem Gedächtnis).
- Die Entropieproduktion des antisymmetrischen Teils verschwindet auf Maschinengenauigkeit.
- Der tangentiale Anteil erreicht bei $\mu B\approx1$ ein Drittel des transversalen und fällt bei starkem Feld wieder, weil der magnetische Widerstand den symmetrischen Teil umformt.
- **Bridgman-Beziehung** als Folge von Onsager-Casimir (Callen 1948): Der Ettingshausen-Koeffizient $P$ und $TN/\kappa$ stimmen bei $B=0{,}01$ T überein: $-2{,}3819\cdot10^{-3}$ gegen $-2{,}3819\cdot10^{-3}$ K·m/(A·T), Verhältnis 1,00000.

---

## 5  Anwendungen

### 5.1 Curie-Punkt als stabile Randlage

Ein Ferrit wird induktiv geheizt. Die Heizleistung folgt dem Quadrat der Magnetisierung, $P=P_0(1-T/T_C)^{2\beta}$ mit $2\beta=0{,}72$ [G] (3D-Heisenberg). Weitere Werte: $T_C=493$ K, $T_0=293$ K, $C=0{,}75$ J/K (1 g Ferrit) und $G=0{,}01$ W/K, also $\tau_0=75$ s (Annahmen). Die Schleifenverstärkung ist

$$
q=-\frac{P_0}{G}\,\frac{2\beta}{T_C}\Big(1-\frac{T^*}{T_C}\Big)^{2\beta-1}\ \longrightarrow\ -\infty\quad(T^*\to T_C),
$$

weil $2\beta<1$ ist. Die Regel, nach der EM heizt, verliert am kritischen Punkt ihre Wirkung, und zwar mit unendlicher Steigung. Jede Abweichung wird deshalb mit beliebig großer Verstärkung zurückgeholt.

- $\nu=1/(1-q)\to0$: Die Rückkopplungssumme wird unterdrückt, nicht verstärkt.
- Der Restabstand fällt wie $T_C-T^*=T_C\big[(T_C-T_0)G/P_0\big]^{1/2\beta}$.
- Die Kippstelle ($q\to+1$) und die Randlage ($q\to-\infty$) sind die beiden Enden derselben Größe.

| $P_0/G$ | $T_C-T^*$ | $q$ | $1/(1-q)$ | $dT^*/d\ln P_0$ | $\tau_{\rm rel}=\tau_0/(1-q)$ | $\sigma_T$ thermisch | $\delta T$ bei 1 % Leistungsrauschen |
|---|---|---|---|---|---|---|---|
| 100 K | 156,27 K | −0,20 | 0,832 | 36,39 K | 62,4 s | $1{,}3\cdot10^{-9}$ K | 0,364 K |
| 1000 K | 39,01 K | −2,97 | 0,252 | 40,54 K | 18,9 s | $9{,}8\cdot10^{-10}$ K | 0,405 K |
| 10⁴ K | 2,122 K | −67,14 | 0,0147 | 2,904 K | 1,10 s | $2{,}6\cdot10^{-10}$ K | 0,029 K |
| 10⁵ K | 0,088 K | −1637 | 0,0006 | 0,122 K | 0,046 s | $5{,}2\cdot10^{-11}$ K | 0,0012 K |

Die thermische Schwankung ist $\sigma_T^2=k_BT^2/\big(C(1-q)\big)$. Die Gegenkopplung unterdrückt sie um $1-q$, wie in Übergangskanten-Sensoren mit elektrothermischer Gegenkopplung (Irwin 1995, aus dem Gedächtnis). Praktisch bestimmt das Rauschen der Versorgung die Schwankung: 29 mK bei 1 % Leistungsrauschen und $P_0/G=10^4$ K. Die Relaxationszeit sinkt dort von 75 s auf 1,1 s.

### 5.2 Laborversuche mit Zahlen

Alle Zahlen gelten für das VO₂-Element aus §1.1 mit Lastwiderstand ($R_c=10{,}967$ kΩ).

**(a) Breite gegen Abstand zur Spitze.** $R_c$ wird aus der gemessenen Kennlinie bestimmt, als steilster negativer Abschnitt. Vorhersage für $\varepsilon=1-R_L/R_c=0{,}001$, 0,01 und 0,0316: Fensterbreite 0,0156, 0,497 und 2,83 mV, Exponent 3/2 (gerechnet 1,504). Der Temperatursprung skaliert mit Exponent 1/2 (gerechnet 0,488).

**(b) Kritische Verlangsamung.** Bei $R_L=0{,}8R_c$ endet der isolierende Ast bei 2,62683 V. Die Abklingrate nach einem kleinen Stoß fällt von 0,930/µs (100 mV vor der Falte) auf 0,0910/µs (0,3 mV) und 0,0293/µs (0,03 mV); die Relaxationszeit wächst auf 11,0 µs bzw. 34,1 µs. Vorhersage: Rate $\propto\sqrt{V_{s,f}-V_s}$ (gerechnet 0,490). Messbar mit einem Pump-Probe-Puls auf die Versorgung.

**(c) Verweilzeiten im Rauschen und Kac-Budget.** Im Fenster 2,57738 bis 2,62683 V (49,45 mV breit) hat das Element zwei Zustände, kalt und heiß. Die Langevin-Gleichung lautet $\dot T=f(T)+\sqrt{2D}\,\xi$, mit dem Potential $U=-\int f\,dT$.
- Das Eigenrauschen $D_{\rm th}=k_BT^2G/C_{\rm th}^2=3{,}07\cdot10^{-7}$ K²/µs ist viel zu klein: $\Delta U/D_{\rm th}=3{,}5\cdot10^6$. Spontane Sprünge gibt es praktisch nie.
- Eingespeist wird deshalb Leistungsrauschen. Für $D=\Delta U/4$ braucht es $\sqrt{S_P}=3{,}65\cdot10^{-9}$ W/√Hz.
- Kramers: $\tau_i=2\pi\big/\sqrt{U''_i\,|U''_s|}\;e^{\Delta U_i/D}$ (Kramers 1940; Hänggi, Talkner, Borkovec 1990).

| $V_s$ | $\Delta U_1$, $\Delta U_2$ [K²/µs] | $D$ | Verweilzeit kalt: Sim / Kramers | heiß: Sim / Kramers | Sprünge |
|---|---|---|---|---|---|
| 2,60629 V | 1,0635; 1,0635 | $\Delta U/3$ | 211 / 185 µs | 110 / 88 µs | 6006 / 6034 |
| 2,60629 V | 1,0635; 1,0635 | $\Delta U/4$ | 530 / 502 µs | 282 / 238 µs | 2244 / 2271 |
| 2,61143 V | 0,6872; 1,3813 | $\Delta U_1/3$ | 246 / 206 µs | 1555 / 1823 µs | 861 / 784 |

Bei gleich hohen Barrieren bleibt das Verweilzeitverhältnis nahe dem Krümmungsverhältnis: simuliert 1,92 und 1,88, nach Kramers 2,11. Es ist unabhängig von $D$. Das trennt Barriere und Krümmung im Experiment.

**Kac-Test.** Zerlegt man in die Zellen kalt und heiß, dann gilt im geschlossenen, stationären System nach Kac (1947) $\langle\tau_{\mathrm{ret},i}\rangle=1/p_i$ und damit
$$
\sum_i\frac1{\langle\tau_{\mathrm{ret},i}\rangle}=1 .
$$
Gemessen: 1,0125, 1,0324 und 1,0386 für die drei Läufe, bei Anteilen $p=(0{,}668;\,0{,}332)$, $(0{,}666;\,0{,}334)$ und $(0{,}113;\,0{,}887)$. Die Abweichung von 1 stammt aus den endlichen Laufzeiten, also aus angeschnittenen Aufenthalten.

Im Versuch prüft die Summe, ob das System geschlossen und stationär ist. Ein dritter Zustand (Durchbrennen), eine Drift durch Alterung oder ein Leck verschieben sie gegen den Anteil der Zeit, der den beiden Zellen verloren geht.

**(d) Torus-Geburt.** Die Brücke aus §1 bei 1 T:
- Hopf bei 3,548 V mit 3,26 MHz,
- zweite Frequenz ab etwa 3,472 V bei 0,51 MHz,
- Einrasten 1:1 ab 3,44 V mit einem Amplitudensprung der Brücke von 0,26 auf 4,0 µm.

Bei 3,47 V schaltet das Feld zwischen den Zuständen: 0 T Zyklus, 1–2 T Torus, 4 T eingerastet. Messgröße ist das Spektrum von Spannung und Auslenkung (optischer Hebel).

---

## 6  Ergebnis

### Neu für das Modell

**N1 · Verzweigungen sind Rollenwechsel.** Hopf und Neimark-Sacker machen transversale Richtungen tangential. Einrasten macht eine tangentiale Richtung über den dissipativen Kopplungsanteil wieder transversal, Chaos macht eine transversale Richtung dehnend. Die Zahl der verschwindenden Lyapunov-Exponenten ist die Dimension des Torus. Der Weg Fixpunkt → Grenzzyklus → Torus ist damit eine Folge von Vorzeichenwechseln im Spektrum, gesteuert von parametrischen Größen.

**N2 · Das Magnetfeld ist ein parametrischer Schalter für Tori.** Am selben Arbeitspunkt erzeugt es den Torus (0 → 1 T) und rastet ihn ein (4 T). Im Transport legt es fest, welcher Anteil umläuft ($L_a$, ungerade in $B$, entropieneutral) und welcher dissipiert ($L_s$).

**N3 · $\nu=1/(1-q)$ ist Verstärkung und Verlangsamung zugleich.** Die Kippstelle liegt bei $q=1$ bzw. $\det(I-G)=0$. Dort divergiert die Rückkopplungssumme, das Hysteresefenster öffnet mit $a^{3/2}$ und die Relaxationsrate schließt mit $\delta^{1/2}$.

**N4 · Galaxien sind monostabil, aber nahe kritisch am Rand.** Mit $q=e^{-\sqrt y}<1$ gibt es keine Falte und keine Hysterese in der RAR. Nach außen wachsen Verstärkung und Relaxationszeit, $\tau_{\rm rel}=\nu\tau_R\propto r^2$ bei flacher Rotationskurve. Das ist kritische Verlangsamung ohne Übergang. Weil die Rückkopplung stationär ist, bleibt sie unter $\nu\,r/c\approx1$ Mio. Jahre (50 kpc), ohne messbaren Nachlauf. [korr.]

**N5 · „Zustand ist Gedächtnis“ ohne Speicherung.** Ein Preisach-Zustand ist die Folge der nicht ausgelöschten Extremwerte. Nach 200 000 Schritten reichen 6 Zahlen. Der Gedächtniskern $\sum\nu_ie^{-h_i\Delta t}$ ist der thermisch geglättete Kleinsignal-Grenzfall desselben Ensembles. Kriterium: Kern für Amplituden unter $kT/|\partial_uE_b|$ und aktive Barrieren, Preisach für eingefrorene Barrieren.

**N6 · Die Schließungsregel ist ein Preisach-Operator.** Es gilt $w=F(\theta_{\min})$. Die Bindung eines Bereichs hängt vom tiefsten Kollaps seiner Geschichte ab, nicht vom heutigen Zustand. Der Kugelkollaps setzt $\theta_{\min}/3H=-8{,}57$ unabhängig von der Masse. Damit sind virialisierte Innenbereiche einheitlich entkoppelt. Im Modell ist $w=0$ schon für jedes $\theta\le0$, also ab der Umkehr; die Geschichte zeigt sich nur dort, wo die Mode ξ noch nicht ganz übergegangen ist (Übergaberate $\lambda_m=\omega_{\rm orb}/N^*$), nicht über gestufte Schwellen. [korr.]

**N7 · Die kritische Stelle als Attraktor der Gegenkopplung.** Hat die Regel am kritischen Punkt unendliche Steigung ($2\beta<1$), dann strebt $q\to-\infty$. Der Arbeitspunkt rastet dort ein, und Schwankungen sinken um $1-q$.

### Prüfbare Aussagen

**Labor**

| Nr | Aufbau | Vorhersage |
|---|---|---|
| B1 | VO₂-Element mit einstellbarem Lastwiderstand | bistabil genau für $R_L<\max|dV/dI|$; Fensterbreite $\propto\varepsilon^{3/2}$ (0,497 mV bei $\varepsilon=0{,}01$), Sprung $\propto\varepsilon^{1/2}$ |
| B2 | dasselbe, kleiner Stoß vor der Falte | Relaxationsrate $\propto(V_{s,f}-V_s)^{1/2}$: 0,930 → 0,091 → 0,029/µs bei 100 → 0,3 → 0,03 mV |
| B3 | dasselbe, eingespeistes Rauschen | Verweilzeiten nach Kramers (Abweichung 5–25 % bei $\Delta U/D=3$–4); Verhältnis bei gleichen Barrieren ≈ 2, unabhängig von $D$; $\sum1/\langle\tau_{\rm ret}\rangle=1$ genau dann, wenn das System geschlossen und stationär ist |
| B4 | VO₂ auf Mikrobrücke im Magnetfeld | Hopf bei 3,548 V (3,26 MHz), zweite Frequenz ab ≈ 3,47 V, Einrasten mit Amplitudensprung; Feld schaltet bei 3,47 V zwischen Zyklus, Torus und Einrasten; bei steiferer Brücke Chaos mit $\lambda_1=0{,}3$–0,9/µs. Die Zahlen hängen von den Annahmen zu $\kappa$, $D_{\rm tr}$, $Q$ ab, die Reihenfolge nicht. |
| B5 | Ferritkern im Induktionsfeld | $T_C-T^*\propto P_0^{-1/2\beta}$; 1 % Leistungsrauschen gibt 29 mK bei $P_0/G=10^4$ K; Relaxation $1/(1-q)$-fach schneller |
| B6 | Bismut, 300 K, 1 T | Nernst-, Hall- und Righi-Leduc-Anteile verschwinden in der Entropiebilanz; $P\kappa=TN$ (Bridgman). Kontrolle mit bekannter Physik. |

**Kosmos**

| Nr | Aussage | Messung |
|---|---|---|
| K1 | Die RAR hat genau einen Ast: keine Bistabilität, kein Hysteresezweig, $d\ln g/d\ln g_N\in(1/2,1)$. | Streuung und Form der RAR |
| K2 | Äußere Beschleunigungen folgen Massenänderungen ohne messbare Verzögerung: Die Rückkopplung ist stationär, $\tau_{\rm rel}\le\nu\,r/c\approx1$ Mio. Jahre (50 kpc); Endantwort 0,54 bis 0,59·$\delta M/M$ zwischen 50 und 20 kpc. Einen Nachlauf mit $1/\lambda_m$ gibt es nur, wo die Mode ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche). [korr.] | Rotationskurven von Galaxien mit datierten Ausbrüchen, Ausflüssen, Gaszufuhr |
| K3 | Material, das einmal voll kollabiert war, bleibt entkoppelt: Gezeitenzwerggalaxien liegen auf der RAR. | Kinematik von Gezeitenzwergen |
| K4 | Späte, schwächere Verschmelzungen hinterlassen keine Spur in der Bindung (Auslöschung). | RAR-Residuen gegen Verschmelzungsgeschichte |
| K5 | Um isolierte Galaxien wie die Milchstraße ist die Mode ξ bis 100 kpc ganz übergegangen, die Linse folgt dort der Relation. Weiter außen senkt der Rest der Mode die Linse unter die Relation: bei $N^*=2{,}4$ um 3–7 % (200 kpc) und 17–44 % (300 kpc), bei $N^*=0{,}3$ nicht messbar. Das misst $N^*$. Die Abschätzung $\delta g/g=Qw/\nu$ (0,21–0,35 bei 100 kpc, 0,07–0,12 bei 300 kpc, 0,02–0,04 bei 1 Mpc) setzte gleichverteilte Schwellen voraus und trägt nicht. [korr.] | Linsensignal um isolierte Galaxien bis zum Umkehrradius |

---

## Literatur

- Benettin G., Galgani L., Giorgilli A., Strelcyn J.-M. 1980, Meccanica 15, 9.
- Bournaud F. et al. 2007, Science 316, 1166 (aus dem Gedächtnis).
- Behnia K., Méasson M.-A., Kopelevich Y. 2007, Phys. Rev. Lett. 98, 166602 (aus dem Gedächtnis).
- Callen H. B. 1948, Phys. Rev. 73, 1349 (Onsager-Relationen für thermomagnetische Effekte, Bridgman-Beziehung; aus dem Gedächtnis).
- Casimir H. B. G. 1945, Rev. Mod. Phys. 17, 343.
- Gentile G., Famaey B., Combes F., Kroupa P., Zhao H. S., Tiret O. 2007, A&A 472, L25 (aus dem Gedächtnis).
- Grmela M., Öttinger H. C. 1997, Phys. Rev. E 56, 6620; Öttinger H. C., Grmela M. 1997, Phys. Rev. E 56, 6633.
- Gunn J. E., Gott J. R. 1972, ApJ 176, 1.
- Hänggi P., Talkner P., Borkovec M. 1990, Rev. Mod. Phys. 62, 251.
- Irwin K. D. 1995, Appl. Phys. Lett. 66, 1998 (aus dem Gedächtnis).
- Kac M. 1947, Bull. Amer. Math. Soc. 53, 1002.
- Kramers H. A. 1940, Physica 7, 284.
- Lelli F., McGaugh S. S., Schombert J. M., Pawlowski M. S. 2017, ApJ 836, 152 (Streuung der RAR; aus dem Gedächtnis).
- Mayergoyz I. D. 1986, Phys. Rev. Lett. 56, 1518.
- Morin F. J. 1959, Phys. Rev. Lett. 3, 34.
- Newhouse S., Ruelle D., Takens F. 1978, Commun. Math. Phys. 64, 35.
- Onsager L. 1931, Phys. Rev. 37, 405; 38, 2265.
- Preisach F. 1935, Z. Phys. 94, 277.
- Scheffer M. et al. 2009, Nature 461, 53 (Frühwarnsignale, kritische Verlangsamung).

---

## Selbstcheck

| Prüffrage | Antwort |
|---|---|
| Wird Gedächtnis als Speicherobjekt eingeführt? | Nein. Der Preisach-Zustand ist die gegenwärtige Treppenlinie im Schwellendreieck. Der Kern ist ein Faltungsgewicht, $w$ eine Funktion des tiefsten Kollapses. Nichts wird aufgezeichnet. |
| Werden Wirkungen als Substanz eingeführt? | Nein. Kopplungen sind Terme, Schwellen und Matrixanteile, das Magnetfeld ist ein Parameter. |
| Stammen alle Zahlen aus dem Skript? | Ja. Stoffwerte, die aus der Literatur stammen, sind gekennzeichnet; Bauteil-, Ferrit- und Galaxienparameter sind Annahmen. |
| Mystik oder Numerologie? | Keine. Verwendet werden Hopf- und Neimark-Sacker-Verzweigung, Lyapunov-Spektren, Sattel-Knoten- und Spitzennormalform, Preisach, Debye/Prony, Onsager-Casimir, GENERIC, Kramers, Kac und der Kugelkollaps. |

**Offen:**
- $\tau_R$, die Dauer einer Runde der Rückkopplung in Galaxien. Im Modell ist die Rückkopplung stationär (höchstens $\nu\,r/c$); K2 ist keine messbare Verzögerung. [korr.]
- Die Verteilung $\mu(\alpha,\beta)$ der Schwellen von $w$. Hier ist sie gleichverteilt angenommen. Im Modell gilt $w=\theta_+^2/(\theta_+^2+\theta_c^2)$, $w=0$ für $\theta\le0$; die Größe von K5 legt $N^*$ fest. [korr.]
- Die Interpolation $g=\nu g_N+wQg_N$ ist eine Annahme für die Abschätzung, nicht aus einer Wirkung abgeleitet.
