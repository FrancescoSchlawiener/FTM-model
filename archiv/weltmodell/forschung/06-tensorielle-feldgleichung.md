# Die tensorielle Feldgleichung

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Die Notiz rechnet mit a₀ = cH₀/2π = 1,08·10⁻¹⁰ (H₀ = 70), wo nicht 1,2·10⁻¹⁰ steht. Geändert: FTM-AS/QUMOND als Grundgleichung der Bewegung mit η = 1 überall → Gleichung der Metrik (wirbelfreier Teil der örtlichen Regel), η = 1 für Kugeln und außerhalb von Scheiben, in Scheiben innen η > 1; Sonne Newtonsch bis zur Gezeitengrenze, homogenes äußeres Feld nicht in ν, Q₂ = 4,7·10⁻³¹ s⁻² → Glättungslänge ℓ ≈ 0,93 pc, äußeres Feld geht in ν ein, Quadrupol aus der Glättung 9·10⁻³¹ s⁻²; ℓ₀ = 2πR_H und „a₀ dynamisch, a₀(z) ∝ H(z)“ → ℓ₀ = 2Λ (innere Länge), a₀ seit z ≈ 10 fest, ≈ cH₀/2π nur heute; kurzes Gedächtnis (λ_m ≫ ω, κ₀ = λ_m²) → Geschwindigkeitskern auf jeder Bahn ausgeschlossen, κ₀/λ_eff = 3c²/(r*·v_rot); Offene Punkte: 1 gelöst (stetige Glättung ρ > ρ_J), 4 und 7 durch die Schließungsregel geschlossen, 5 entfällt, 6 beantwortet (Vakuumteil = Energie der Rückstellung), 2 und 3 offen; [G] an Drehimpuls der Sonne, VLBI-Genauigkeit, Euclid-Faktor und Faktor 2 in Haufen. Umbenennungen: Λ_obs und nackte Konstante Λ → Λ_kosm, Ω_Λ → Ω_DE; a_K → a_J; Diffeomorphismus ξ → ε^μ; PPN-Parameter ξ → κ_W; ρ_ph → gelesene Dichte ρ_X; q = e^{2πiτ_M} als Nome in Punkt 1 gekennzeichnet (für die Punktmasse der Rückkehranteil). 𝒬 = u·∇σ steht schon als \mathcal Q.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/tensor.py` (alle Zahlen dieser Notiz).

---

## Kurzfassung

| Frage | Ergebnis |
|---|---|
| Kovariante Form | Einstein-Aether-Skalar-Theorie „FTM-AS“: Metrik $g_{\mu\nu}$, Einheits-Zeitfeld $u^\mu$ (Abwicklungsrichtung des Mediums), Faltungspotential $\sigma$. Materie koppelt an $g_{\mu\nu}$; der Wirbelanteil der örtlichen Regel wirkt zusätzlich auf ihre Bewegung [korr.]. |
| Strukturterm | $\mathcal S_{\mu\nu}$ entsteht allein aus Geometrie ($g$, $u$, $\sigma$), ohne Kopplung an Materie. Er ist eine Summe über Potenzen von $q=e^{2\pi i\tau_M}$ und verschwindet für $q\to0$ in jeder Ordnung. Wegen der Diffeomorphismen-Invarianz gilt $\nabla^\mu\mathcal S_{\mu\nu}=0$. |
| Schwacher Grenzfall | ergibt exakt $\nabla^2\Phi_N=4\pi G\rho$, $\nabla^2\Phi=\nabla\cdot[\nu\nabla\Phi_N]$ und dazu $\Psi=\Phi$. Numerisch geprüft. |
| Faltung und Takt | $\Delta s$ ist der Excess-Radius aus $g_{rr}$ (Feynman), $d\tau/dt=\sqrt{-g_{tt}}$. Die Beziehung $d\tau/dt=\sqrt{1-6\Delta s/r}$ **folgt** aus $G^t{}_t=G^r{}_r$ im Vakuum. |
| Klassische Tests | Lichtablenkung 1,7512″, Merkur 42,981″/Jh, GP-B 6603,9 bzw. 39,05 mas/Jahr, LAGEOS 30,66 mas/Jahr. Der Zusatzterm der isolierten Quelle liegt bei $10^{-102}$ bis $10^{-656489}$. PPN: $\gamma=\beta=1$ in jeder Ordnung. |
| Sonnensystem im Feld der Milchstraße | Die Sonne ist ein Einzeltrichter und bleibt bis zu ihrer Glättungslänge $\ell\approx0{,}93$ pc ($27\,r^*$) im Feld der Milchstraße Newtonsch; 1,37 pc ist der klassische Jacobi-Radius. Die Regel sieht den Betrag des Gesamtfelds, das äußere Feld geht in $\nu$ ein. Der Quadrupol aus der Glättung ist $9\cdot10^{-31}$ s⁻², weit unter der Cassini-Schranke $(3\pm3)\cdot10^{-27}$ s⁻². [korr.] |
| Linsen- gegen dynamische Masse | Für Kugeln und außerhalb von Scheiben ist $\eta=M_\text{dyn}/M_\text{lens}=1$, in Übereinstimmung mit KiDS-1000; in Scheiben innen ($R\lesssim3{,}5\,R_d$) ist $\eta>1$, weil die Bahnen der örtlichen Regel folgen und die Metrik nur deren wirbelfreien Teil trägt [korr.]. Prüfbar mit Euclid über Linsen-RAR und Linsen-Tully-Fisher-Nullpunkt $v^4=GMa_0$. |

---

## 1. Kovariante Form

### 1.1 Theorieklasse

[korr.] Die Fassung unten ist die Gleichung der Metrik, nicht die der Bahnen: Die Bahnen folgen der örtlichen Regel $a=\nu(|g_N|)\,g_N$; die Metrik (Feldform) trägt nur deren wirbelfreien Teil. An SPARC geht das senkrechte Scheibenfeld mit $f=+0{,}03$ bis $+0{,}13$ ein, ohne IC 2574 mit $0{,}00\pm0{,}03$; die Feldform ($f=1$) liegt außerhalb 95 %. Für Kugeln sind beide gleich, die Zahlen bleiben dort.

Die nichtrelativistische Feldform der Metrik [korr.] ist dem Aufbau nach **QUMOND** (Milgrom 2010). Es gibt ein Hilfspotential $\Phi_N$, das die Poisson-Gleichung der Baryonen erfüllt, und das physikalische Potential $\Phi$ folgt aus einer zweiten, linearen Gleichung mit dem Operator $\nu$. Die kovariante Fassung gehört zur Klasse der **Vektor-Skalar-Theorien mit nichtkanonischer Kinetik**, in der $\Phi=\Psi$ möglich ist. Vorbilder: verallgemeinerter Einstein-Aether (Zlosnik, Ferreira & Starkman 2007), TeVeS (Bekenstein 2004), **AeST** (Skordis & Złośnik 2021).

Der Grund für den Vektor: Eine MOND-artige Wirkung muss **in erster Ordnung** der Potentiale wirken, also im selben Maß wie $\rho$. Die Feldenergie eines Skalars ist dagegen quadratisch und damit um $v^2/c^2\sim10^{-6}$ zu klein. In AeST sorgt ein bilinearer Mischterm $J^\mu\nabla_\mu\phi$ dafür: Im statischen Fall ist $J^\mu$ die Beschleunigung des Einheitsvektors, also $\nabla\Phi$. Die FTM-Fassung übernimmt genau diesen Mechanismus in der einfachsten Form.

### 1.2 Die Wirkung (FTM-AS)

Felder und ihre FTM-Lesart:

| Feld | Bedeutung im FTM |
|---|---|
| $g_{\mu\nu}$ | Geometrie des Faltungsmediums; Licht und Materie laufen auf ihr |
| $u^\mu$, $u_\mu u^\mu=-1$ | Abwicklungsrichtung des Mediums (Takt). $h^{\mu\nu}=g^{\mu\nu}+u^\mu u^\nu$ projiziert auf den Raum, den $u$ sieht |
| $a_\mu=u^\nu\nabla_\nu u_\mu$ | Beschleunigung des Mediums. Statisch gilt $c^2 a_i=\partial_i\Phi$, das ist die **Neigung** $g$ (numerisch geprüft: Abweichung $4{,}7\cdot10^{-10}$) |
| $\sigma$ (dimensionslos) | Faltungspotential des Einzeltrichters. Im schwachen Feld gilt $c^2\sigma=\Phi_N$ |

$$
\boxed{\;
S=\frac{c^3}{16\pi G}\int d^4x\,\sqrt{-g}\,\Big[R\;+\;\underbrace{2\,h^{\mu\nu}B_\mu B_\nu+\frac{2}{\ell_0^2}\,\tilde Q(q)}_{L_s\ \text{(Struktur)}}\;-\;\lambda_u\,(u_\mu u^\mu+1)\Big]\;+\;S_m[g_{\mu\nu},\psi]\;}
$$

$$
B_\mu=a_\mu-\nabla_\mu\sigma,\qquad \mathcal Y=h^{\mu\nu}\nabla_\mu\sigma\nabla_\nu\sigma,\qquad
\ell_0=\frac{c^2}{a_0}=2\Lambda\ \ (\text{heute}\approx2\pi c/H_0=2\pi R_H),
$$

$$
q=e^{2\pi i\tau_M},\qquad 2\pi i\tau_M=-\big(\ell_0^2\,\mathcal Y\big)^{1/4}=-\sqrt{\frac{|\nabla\Phi_N|}{a_0}}\;\stackrel{\text{Punktmasse}}{=}\;-\frac{r^*}{r},\qquad r^*=\sqrt{\frac{3c^2\Delta s}{a_0}} .
$$

[korr.] $\ell_0=2\Lambda$ ist eine Länge der inneren Richtung der Rückkehr, keine Raumlänge; $\Lambda$ ist heute fest, $\approx\pi c/H_0$ gilt nur heute. $q$ ist hier der Nome des Bahntorus in Punkt 1 des Autors ($2\pi i\tau_M=-s$), für die Punktmasse zugleich der Rückkehranteil $q=e^{-s}$; der Kern-Nome des Modells heißt $q_{\Delta s}=\Delta s/r$ (Punkt 3). Die Wahl ist offen und betrifft in der stationären Regel nur die Schreibweise.

$\lambda_u$ ist der Lagrange-Multiplikator der Normierung $u_\mu u^\mu=-1$. Die Kinetikfunktion ist die aufintegrierte Rückkopplungssumme mit der Normierung $\tilde Q\to0$ für $q\to0$:

$$
\frac{d\tilde Q}{d\mathcal Y'}=\nu-1=\sum_{n\ge1}q^n,\qquad
\tilde Q=-4\sum_{n\ge1}q^n\Big(\frac{s^3}{n}+\frac{3s^2}{n^2}+\frac{6s}{n^3}+\frac{6}{n^4}\Big)
=-4\big[s^3\mathrm{Li}_1(q)+3s^2\mathrm{Li}_2(q)+6s\,\mathrm{Li}_3(q)+6\,\mathrm{Li}_4(q)\big],
$$

mit $s=-\ln q$ und $\mathcal Y'=\ell_0^2\mathcal Y$. Numerische Probe: Reihe, Polylog-Form und direktes Integral stimmen auf 9 bis 10 Stellen überein. Die numerische Ableitung $d\tilde Q/d\mathcal Y'$ trifft $\nu-1$ bei $s=0{,}05$ bis $12$ auf alle ausgegebenen Stellen. Grenzwerte:

- $q\to0$: $\tilde Q\approx-4q(s^3+3s^2+6s+6)\to0$. Bei $s=12$ ist $\tilde Q=-0{,}055$.
- $q\to1$: $\tilde Q(0)=-4\pi^4/15=-25{,}976$.

Mathematische Randbemerkung ohne physikalische Deutung: $\nu-1=1/(e^s-1)$ hat dieselbe Form wie das Bose-Integral, daher $\int_0^\infty 4s^3/(e^s-1)\,ds=4\pi^4/15$.

In FTM-Größen stehen im Strukturterm also:
- $\Delta s$ über $r^*$,
- $a_0=c^2/(2\Lambda)$ über $\ell_0=2\Lambda$, heute $\approx cH_0/2\pi$ bzw. $\ell_0\approx2\pi R_H$ [korr.] ($\ell_0=8{,}30\cdot10^{26}$ m $=26\,909$ Mpc; Modellwert $a_0=1{,}08\cdot10^{-10}$ m/s² bei $H_0=70$ km/s/Mpc; die Zahlen dieser Notiz nutzen diesen Wert, außer wo der gemessene Wert $1{,}2\cdot10^{-10}$ m/s² ausdrücklich steht),
- $q=e^{2\pi i\tau_M}$ als einzige Variable von $\tilde Q$.

Probe für die Sonne bei 5000 AE: $(\ell_0^2\mathcal Y)^{1/4}=r^*/r=1{,}4804$, $\tau_M=0{,}2356\,i$, $q=0{,}2276$.

### 1.3 Feldgleichungen

**Variation nach $\sigma$** (Faltungsgleichung):

$$
\boxed{\;\nabla_\mu\Big[h^{\mu\nu}\big(\nu(q)\,\nabla_\nu\sigma-a_\nu\big)\Big]=0\;}
$$

Das ist die FTM-Feldgleichung in kovarianter Form. Die räumliche Divergenz der **Neigung des Mediums** $a$ ist gleich der Divergenz des mit $\nu$ gewichteten Einzeltrichter-Feldes.

**Variation nach $u^\alpha$** (Aether-Gleichung, $\mathcal Q=u^\mu\nabla_\mu\sigma$, $\hat B^\mu=h^{\mu\nu}B_\nu$):

$$
-4\mathcal Q B_\alpha+4\hat B^\mu\nabla_\alpha u_\mu-4\nabla_\nu(\hat B_\alpha u^\nu)+4(\nu-1)\mathcal Q\nabla_\alpha\sigma-2\lambda_u u_\alpha=0 .
$$

Statisch und entlang des Killing-Vektors gilt $\mathcal Q=0$, $\nabla_\alpha u_\mu=-u_\alpha a_\mu$, $\theta_u=0$. Die Gleichung reduziert sich dann auf $-(8\hat B\cdot a+2\lambda_u)u_\alpha=0$. Ihre räumlichen Komponenten verschwinden identisch, und die Zeitkomponente legt $\lambda_u=-4\hat B\cdot a$ fest. Das statische Medium ist also eine konsistente Lösung. Dasselbe folgt aus der Zeitumkehr: $S$ ist gerade in $u$, die Gleichung daher ungerade.

**Variation nach $g^{\alpha\beta}$** (bei festem $u^\mu$):

$$
\boxed{\;G_{\alpha\beta}=\frac{8\pi G}{c^4}\,T_{\alpha\beta}+\mathcal S_{\alpha\beta}\;}
$$

$$
\mathcal S_{\alpha\beta}=-2B_\alpha B_\beta-2(\nu-1)\nabla_\alpha\sigma\nabla_\beta\sigma+\tfrac12 g_{\alpha\beta}L_s-\lambda_u u_\alpha u_\beta
-4\nabla_\lambda\big(u^\lambda\hat B_{(\alpha}u_{\beta)}\big)+2\nabla_\lambda\big(\hat B^\lambda u_\alpha u_\beta\big).
$$

**Energieerhaltung.** $S_s=\int\sqrt{-g}L_s$ ist ein Skalar aus $g$, $u$ und $\sigma$. Unter einem Diffeomorphismus $\varepsilon^\mu$ [korr.] gilt $\delta S_s=0$. Daraus folgt die Noether-Identität

$$
\nabla^\alpha\mathcal S_{\alpha\beta}\;\propto\;E_\sigma\nabla_\beta\sigma+(\text{Terme}\propto E_u)\;=\;0\quad\text{auf der Lösung}.
$$

Materie koppelt minimal an $g$, also gilt $\nabla^\mu T_{\mu\nu}=0$. Beides passt zur Bianchi-Identität $\nabla^\mu G_{\mu\nu}\equiv0$.

**Struktur, keine Quelle.** $\sigma$ und $u$ koppeln nicht an Materie. $\sigma$ wird allein durch die Geometrie festgelegt: Im Starkfeld gilt $\sigma=\ln N$, im schwachen Feld wird $c^2\sigma$ zu $\Phi_N$ über die Metrikgleichung. $\mathcal S_{\mu\nu}$ ist damit ein Term der Geometrie. Wer ihn in $T_{\mu\nu}$ verschiebt und die gelesene Dichte $\rho_X=\nabla\cdot[(\nu-1)\nabla\Phi_N]/4\pi G$ (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt) als „Dunkle-Materie-Dichte“ einführt, deutet Geometrie als Stoff (Fehlschluss des Beobachters). [korr.] Im quasistatischen Limes sieht diese Dichte wie druckloser Staub aus (1.5). Genau deshalb wird sie für einen Teilchenstoff gehalten.

### 1.4 Grenzfall $q\to0$: Einstein exakt

Für $q\to0$ gilt $\nu-1=O(q)$. Die $\sigma$-Gleichung wird im statischen Fall zu $D_i\big(N\,D^i(\sigma-\ln N)\big)=0$ und hat die Lösung $\sigma=\ln N$, also $B_\mu=0$. Damit werden auch $\lambda_u=0$, $L_s=O(q)$ und $\mathcal S_{\mu\nu}=O(q)$. Es gilt $G_{\mu\nu}=8\pi G\,T_{\mu\nu}/c^4$. Dasselbe gilt für jede stationäre Lage mit $u$ entlang des Killing-Vektors, denn dort ist $a_\mu=\nabla_\mu\ln N$ immer ein Gradient.

$q=\exp(-\sqrt{GM/(a_0r^2)})$ fällt schneller ab als jede Potenz von $GM/(rc^2)$, denn in PN-Zählung gilt $q=\exp(-c\sqrt{U/(a_0 r)})$. **Deshalb trägt der Strukturterm zu keinem Koeffizienten der post-Newtonschen Entwicklung bei** (2.3).

Tensormoden: $u$ und $\sigma$ koppeln auf flachem Hintergrund linear nicht an transversal-spurfreie Störungen. Daher gilt $c_T=c$, wie von GW170817 verlangt.

### 1.5 Schwacher quasistatischer Grenzfall

Für $ds^2=-(1+2\Phi/c^2)c^2dt^2+(1-2\Psi/c^2)d\vec x^2$ wird die Einstein-Hilbert-Wirkung zu $(1/8\pi G)[(\nabla\Psi)^2-2\nabla\Phi\cdot\nabla\Psi]$. Numerisch aus $\sqrt{-g}R$ geprüft: Koeffizient $-0{,}663700$ gegen $-0{,}663698$, relative Abweichung $2{,}3\cdot10^{-6}$.

Mit $c^2a_i=\partial_i\Phi$ und $c^2\sigma=\Phi_N$:

$$
L_\text{NR}=\frac1{8\pi G}\Big[(\nabla\Psi)^2-2\nabla\Phi\cdot\nabla\Psi+(\nabla\Phi-\nabla\Phi_N)^2+a_0^2\,\tilde Q\Big]-\rho\Phi .
$$

| Variation | Gleichung |
|---|---|
| $\delta\Psi$ | $\nabla^2(\Phi-\Psi)=0\Rightarrow\boxed{\Psi=\Phi}$ |
| $\delta\Phi$ | $\nabla^2\Psi-\nabla^2(\Phi-\Phi_N)=4\pi G\rho\Rightarrow\boxed{\nabla^2\Phi_N=4\pi G\rho_\text{bar}}$ |
| $\delta\Phi_N$ | $\boxed{\nabla^2\Phi=\nabla\cdot[\nu(|\nabla\Phi_N|/a_0)\nabla\Phi_N]}$, $\nu=\sum_{n\ge0}q^n$ |

Das ist die Feldgleichung der Metrik [korr.] $\nabla^2\Phi=\nabla\cdot[\nu\nabla\Phi_N]$, zusätzlich mit $\Psi=\Phi$. Die diskrete Variationsprobe an einer Plummer-Kugel bestätigt das: Die Konfiguration mit $\Psi=\Phi$ macht die Wirkung stationär (Residuen $\le6\cdot10^{-7}$).

Kugelsymmetrisch lautet die 00-Komponente in reinen FTM-Größen

$$
G_{00}=\frac{6}{r^2}\,\frac{d}{dr}\big[\Delta s_\text{eff}(r)\big],\qquad \Delta s_\text{eff}=\nu\,\Delta s(r),\qquad
\mathcal S_{00}=\frac6{r^2}\frac{d}{dr}\big[\Delta s_\text{eff}-\Delta s\big].
$$

Der Strukturterm ist die radiale Ableitung der **Zusatzfaltung** $\Delta s_\text{eff}-\Delta s=(\nu-1)\,\Delta s$. Für $q\to0$ bleibt $G_{00}=(6/r^2)\,\Delta s(r)'=8\pi G\rho/c^2$.

### 1.6 Faltung und Takt aus einer einzigen Metrik (Zusatzanforderung des Autors)

**(i) $\Delta s$ ist der räumliche Excess-Radius aus $g_{rr}$.** Die $G^t{}_t$-Gleichung allein liefert $g_{rr}=1/(1-2Gm(r)/c^2r)$. Für eine homogene Kugel ist der Überschuss des gemessenen Radius über den Flächenradius

$$
\delta R=\int_0^R\big(\sqrt{g_{rr}}-1\big)dr=R\Big[\frac{\arcsin\sqrt x}{\sqrt x}-1\Big]\;\xrightarrow{\;x=2GM/c^2R\ \ll 1\;}\;\frac{GM}{3c^2}=\Delta s .
$$

Die FTM-Faltungstiefe ist genau Feynmans „excess radius“ (Feynman Lectures II, Kap. 42). Werte:

| Körper | $\Delta s$ | $\delta R$ exakt | $\delta R/\Delta s$ |
|---|---|---|---|
| Erde | 1,478343 mm | 1,478343 mm | 1,000001 |
| Sonne | 492,2083 m | 492,2093 m | 1,000002 |
| Sirius B | 501,08 m | 501,20 m | 1,000233 |
| Neutronenstern (1,4 M☉, 12 km) | 689,1 m | 824,8 m | 1,196901 |

Damit ist auch die Aussage bestätigt, dass $\Delta s=GM/3c^2$ „nur von der Masse abhängt“: In erster Ordnung gilt das für homogene Körper. Für andere Profile ist $\delta R=(G/c^2)\int m(r)/r\,dr$; $GM/3c^2$ ist dann die homogene Normierung. Im Starkfeld wächst der Excess, beim Neutronenstern um 20 %.

**(ii) $d\tau/dt=\sqrt{-g_{tt}}$ und die Kopplung an $\Delta s$ folgt.** Mit $ds^2=-e^{2\alpha}c^2dt^2+e^{2\beta}dr^2+r^2d\Omega^2$ gilt $G^t{}_t-G^r{}_r=-(2/r)e^{-2\beta}(\alpha'+\beta')$. Im Vakuum ist die rechte Seite null. Daraus folgt $\alpha=-\beta$, also $g_{tt}g_{rr}=-1$. Numerische Probe bei $r=10\,GM/c^2$: Für Schwarzschild ist $g_{tt}g_{rr}=-1{,}000000$ und $G^r{}_r=-3\cdot10^{-9}$ (numerisch null). Eine Metrik mit $g_{tt}g_{rr}=-0{,}8$ (etwa gekrümmter Takt bei $g_{rr}=1$) hat $G^r{}_r=+2{,}5\cdot10^{-3}$ und ist keine Vakuumlösung.

Takt und Faltung sind also zwei Komponenten derselben Lösung. Außerhalb des Körpers gilt $g_{rr}=1/(1-6\Delta s/r)$ und damit $\sqrt{-g_{tt}}=\sqrt{1-6\Delta s/r}$. Numerisch identisch auf 12 Stellen. Die Formel des Autors wird also hergeleitet und nicht gesetzt.

**(iii) Mit Strukturterm.** Weil $\mathcal S_{\mu\nu}$ im quasistatischen Limes keine anisotrope Spannung erster Ordnung hat, bestimmt dieselbe Zusatzfaltung $\Delta s_\text{eff}=\nu\Delta s$ beide Komponenten: $g_{rr}=1/(1-6\Delta s_\text{eff}/r)$ und $\Phi'=3c^2\Delta s_\text{eff}/r^2$. Diese Einheit von Faltung und Takt ist genau $\Psi=\Phi$; für Kugeln und außerhalb von Scheiben folgt daraus $\eta=1$ (Abschnitt 3) [korr.].

### 1.7 Kovarianter Projektor und die Integrodifferentialgleichung (Zusatzanforderung)

Der räumliche Projektor $\mathbf e_r\otimes\mathbf e_r$ wird ersetzt durch

$$
\Pi^\mu{}_\nu=n^\mu n_\nu,\qquad n_\mu=\frac{D_\mu\sigma}{|D\sigma|},\qquad D_\mu=h_\mu{}^\nu\nabla_\nu .
$$

$n$ ist die Normale der Niveauflächen des Faltungspotentials im Ruheraum des Mediums; es gilt $n\cdot u=0$. Für eine Punktmasse ist $D\sigma\propto\nabla\Phi_N\propto\mathbf e_r$, also $n=\mathbf e_r$. Kovariante Fassung der Bahngleichung mit Gedächtnis für eine Weltlinie $x(\tau_p)$ mit Vierergeschwindigkeit $U$:

$$
U^\nu\nabla_\nu U^\mu=-\chi\,\kappa_0\,\mathcal P_U{}^\mu{}_\rho\,\Delta\big(\tau_M(x)\big)\,\Pi^\rho{}_\nu(x)\int_{-\infty}^{\tau_p}e^{-\lambda_m(\tau_p-\tau')}\,g^\nu{}_{\nu'}(x,x')\,h^{\nu'}{}_{\alpha'}(x')\,U^{\alpha'}(\tau')\,d\tau' .
$$

Hier ist $g^\nu{}_{\nu'}$ der Paralleltransport entlang der Bahn und $\mathcal P_U=\delta+UU/c^2$ der Projektor, der $U\cdot U=-c^2$ erhält. $h\,U$ ist die Geschwindigkeit relativ zum Medium.

**Grenzfall.** Im schwachen Feld bei langsamer Bewegung gilt:
- $\tau_p\to t$, $U\to(c,\dot{\mathbf r})$, $hU\to\dot{\mathbf r}$,
- $g^\nu{}_{\nu'}\to\delta$, $\Pi\to\mathbf e_r\otimes\mathbf e_r$ **am Zeitpunkt $t$**,
- der Geodätenteil wird zu $\ddot{\mathbf r}+\nabla\Phi$.

Damit folgt exakt die Gleichung des Autors:

$$
\ddot{\mathbf r}=-\nabla\Phi_\text{bar}-\chi\kappa_0\int_{-\infty}^t\Delta(\tau_M)\,(\mathbf e_r\otimes\mathbf e_r)\,e^{-\lambda_m(t-t')}\dot{\mathbf r}(t')\,dt' .
$$

Die Auswertung von $\mathbf e_r$ zur Zeit $t$ folgt dabei von selbst, weil $\Pi$ am Endpunkt steht.

**Verbindung zur Metrik-Theorie.** [korr.] Diese Rechnung setzt einen Geschwindigkeitskern auf jeder Bahn mit kurzem Gedächtnis ($\lambda_m\gg\omega_\text{orb}$, $\kappa_0=\lambda_m^2$) voraus; das trägt nicht: Die Rückkehr ist eine stationäre Schleife am Ort ($\tau_D=79$ Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht $w$, Übergaberate $\lambda_m=\omega_\text{orb}/N^*$, $N^*\approx0{,}3$–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn; ein Geschwindigkeitskern würde exzentrische Bahnen runden (hier ist das Glied $\propto n\,(n\cdot v)$ dissipativ), alte Scheiben verlangen dafür $\lambda_m\ll0{,}01\,\Omega$. In Galaxien gibt es keinen messbaren Nachlauf. Im Modell sind $\kappa_0$ und $\lambda_m$ über $\kappa_0/\lambda_\text{eff}=3c^2/(r^*v_\text{rot})$, $\lambda_\text{eff}=\omega+\lambda_m^2/\omega$, an das Fernfeld gebunden. Für kurzes Gedächtnis ($\lambda_m\gg\omega_\text{orb}$) wird der Kern zu $-(\kappa_0q/\lambda_m)\,n\,(n\cdot v)+(\kappa_0q/\lambda_m^2)\,n\,(n\cdot\dot v)+\dots$. Auf der Kreisbahn mit $M_k=q$ und $\kappa_0=\lambda_m^2$ ergibt sich unter dieser Annahme selbstkonsistent $a=a_N/(1-q)=\nu a_N$, also die Rückkopplungssumme [korr.]. Numerisch für die Milchstraße:

| $r/r^*$ | $\nu$ | $a/a_N$ bei $\lambda_m/\omega_K=1$ | 10 | 1000 |
|---|---|---|---|---|
| 0,5 | 1,15652 | 1,06995 | 1,15446 | 1,15652 |
| 1 | 1,58198 | 1,20072 | 1,56789 | 1,58198 |
| 5 | 5,51666 | 1,48991 | 4,60213 | 5,51652 |

Einordnung: Die Integrodifferentialgleichung wirkt als Kraft auf Bahnen und ist eine Bahngleichung, keine metrische Wirkung. Der Term erster Ordnung $\propto n\,(n\cdot v)$ ist dissipativ. Aus einer Wirkung folgt nur die Metrik (FTM-AS: Licht, Feldform); die Bahnen folgen der örtlichen Regel $a=\nu(|g_N|)\,g_N$, die außer bei Kugeln keine Wirkung hat, und die Metrik trägt nur deren wirbelfreien Teil (an SPARC $f=+0{,}03$ bis $+0{,}13$, ohne IC 2574 $0{,}00\pm0{,}03$; die Feldform $f=1$ liegt außerhalb 95 %) [korr.]. Auf Kreisbahnen gibt die Integrodifferentialgleichung dasselbe wie FTM-AS; bei radialer Bewegung unterscheiden sich beide. **Vorschlag:** FTM-AS ist die Gleichung der Metrik, die örtliche Regel die der Bahnen; als Geschwindigkeitskern auf jeder Bahn ist die Integrodifferentialgleichung ausgeschlossen. [korr.]

### 1.8 Nebenbefund: ein Vakuumterm

Wegen $\tilde Q(q\to0)=0$ ist das Starkfeld frei von einer kosmologischen Konstanten. In feldfreien Gebieten ($\mathcal Y=0$, $q=1$) wirkt $\tfrac12 g_{\alpha\beta}L_s$ dagegen wie $-\Lambda_s g_{\alpha\beta}$ mit

$$
\Lambda_s=\frac{4\pi^4}{15}\Big(\frac{a_0}{c^2}\Big)^2\stackrel{a_0=cH_0/2\pi}{=}\frac{\pi^2}{15}\frac{H_0^2}{c^2}=3{,}77\cdot10^{-53}\ \text{m}^{-2},\qquad \frac{\Lambda_s}{\Lambda_\text{kosm}}=\frac{\pi^2}{45\,\Omega_\text{DE}}=0{,}320 .
$$

Beobachtete kosmologische Konstante $\Lambda_\text{kosm}$, Anteil $\Omega_\text{DE}$ [korr.]. Das Vorzeichen ist positiv, die Größenordnung ist die der Dunklen Energie. Es passt zur FTM-Lesart „Rückstellung der Geometrie, kein Stoff“. Vorsicht ist geboten, weil eine additive Konstante in $\tilde Q$ Konvention ist. Der Befund ist erst dann physikalisch, wenn keine nackte Konstante $\Lambda_\text{kosm}$ auftritt [korr.]. Das ist Teil der offenen Punkte.

---

## 2. Klassische Tests in FTM-Größen

### 2.1 Formeln und Werte

$\Delta s_\odot=492{,}2083$ m, $\Delta s_\oplus=1{,}4783$ mm. Wegen $\gamma=\beta=1$ (2.4) gelten die ART-Formeln, umgeschrieben mit $GM=3c^2\Delta s$:

| Test | FTM-Formel | Wert | Messung |
|---|---|---|---|
| Lichtablenkung am Sonnenrand | $\alpha=12\Delta s/b$ | 1,7512″ | VLBI: $\gamma$ auf ~$10^{-4}$ [G] |
| Merkur-Perihel | $\Delta\varphi=18\pi\Delta s/(a(1-e^2))$ | $5{,}0187\cdot10^{-7}$ rad/Umlauf = 42,981″/Jh | ≈ 43″/Jh |
| Shapiro (Cassini-Geometrie $b=1{,}6R_\odot$, $r_C=7{,}4$ AE) | $\Delta t=\frac{12\Delta s}{c}\frac{1+\gamma}{2}\ln\frac{4r_\oplus r_C}{b^2}$ | 259,9 μs (hin und zurück) | $\gamma-1=(2{,}1\pm2{,}3)\cdot10^{-5}$ |
| Geodätische Präzession (GP-B) | $\Omega=\tfrac92(\Delta s_\oplus/a)\,\omega_\text{orb}$ | 6603,9 mas/Jahr | $6601{,}8\pm18{,}3$ |

### 2.2 Drehfaltung und Frame-Dragging

Analog zu $GM=3c^2\Delta s$ wird für rotierende Falten definiert:

$$
\boxed{\;GJ=3c^3\,\Delta s_J\qquad\Longleftrightarrow\qquad \Delta s_J=\frac{GJ}{3c^3}=\Delta s\cdot a_J,\quad a_J=\frac{J}{Mc}\;}
$$

Die Drehfaltung ist eine **Fläche**: Faltungstiefe mal Kerr-Hebelarm (Verwindung $a_J$ [korr.]). Das gravitomagnetische Potential lautet $g_{0i}=-6\,\Delta s_J\,(\mathbf n\times\mathbf x)_i/r^3$, und die Präzession einer Kreiselachse ist $\boldsymbol\Omega_{LT}=(3c\,\Delta s_J/r^3)\,[3(\mathbf n\cdot\hat r)\hat r-\mathbf n]$.

| | $J$ [kg m²/s] | $a_J$ [korr.] | $\Delta s_J$ |
|---|---|---|---|
| Erde | $5{,}8585\cdot10^{33}$ | 3,2721 m | $4{,}8374\cdot10^{-3}$ m² |
| Sonne | $1{,}92\cdot10^{41}$ [G] | 322,08 m | $1{,}5853\cdot10^{5}$ m² |

| Test | FTM-Formel | Wert | Messung |
|---|---|---|---|
| GP-B, polar, Bahnmittel | $\langle\Omega\rangle=\tfrac32c\Delta s_J/a^3$ | 40,80 mas/Jahr, projiziert auf IM Peg ($\delta=16{,}84°$): **39,05** | $37{,}2\pm7{,}2$ (Everitt et al. 2011) |
| LAGEOS, Knoten | $\dot\Omega=6c\Delta s_J/(a^3(1-e^2)^{3/2})$ | 30,66 mas/Jahr | ART-Wert auf wenige % (Ciufolini et al.) |
| LAGEOS 2, Knoten | ebenso | 31,49 mas/Jahr | ebenso |

### 2.3 Wert von $q$ und Größe des Zusatzterms je Test

Relativer Zusatz ist $\nu-1\approx q=e^{-s}$ mit $s=\sqrt{a_N/a_0}$. Die Spalte „rel. Zusatz“ ist konservativ mit dem größeren $a_0=1{,}2\cdot10^{-10}$ gerechnet.

| Test | $a_N$ [m/s²] | $\log_{10}q$ ($cH_0/2\pi$) | $\log_{10}q$ ($1{,}2\cdot10^{-10}$) | rel. Zusatz | Messgenauigkeit |
|---|---|---|---|---|---|
| Lichtablenkung, $b=R_\odot$ | 274,2 | −691 232 | −656 489 | $10^{-656489}$ | $10^{-4}$ |
| Merkur, $r=a$ | 0,03957 | −8304 | −7887 | $10^{-7887}$ | $10^{-4}$ |
| Cassini, Strahl bis 7,4 AE | $1{,}083\cdot10^{-4}$ | −434 | −413 | $10^{-413}$ | $2{,}3\cdot10^{-5}$ |
| Erdbahn | $5{,}93\cdot10^{-3}$ | −3215 | −3053 | $10^{-3053}$ | – |
| Neptun | $6{,}56\cdot10^{-6}$ | −107 | −102 | $10^{-102}$ | – |
| GP-B | 8,071 | −118 594 | −112 634 | $10^{-112634}$ | 0,19 |
| LAGEOS | 2,648 | −67 923 | −64 509 | $10^{-64509}$ | ~0,05 |
| zum Vergleich: Sonne bei $r^*$ | $1{,}2\cdot10^{-10}$ | −0,5 | −0,4 | 0,582 | – |

Für eine **isolierte** Quelle liegt der Zusatzterm also um Hunderte bis Hunderttausende Größenordnungen unter jeder Messgenauigkeit.

### 2.4 PPN-Parameter

Im PPN-Ansatz $g_{00}=-1+2U-2\beta U^2$, $g_{ij}=(1+2\gamma U)\delta_{ij}$ trägt $\mathcal S_{\mu\nu}=O(q)$ zu keinem Koeffizienten bei (1.4). Daher gilt

$$
\gamma=1,\qquad \beta=1,\qquad \eta_\text{Nordtvedt}=4\beta-\gamma-3=0,\qquad \kappa_W=\zeta_i=0 .
$$

Whitehead-Parameter $\kappa_W$ [korr.]. Selbst im MOND-Regime bleibt $\Psi=\Phi$ (für Kugeln und außerhalb von Scheiben $\eta=1$ [korr.]), und damit $\gamma=1$ auf allen Skalen. **Nicht** gezeigt sind die Vorzugssystem-Parameter $\alpha_1,\alpha_2$: Das Sonnensystem bewegt sich relativ zum kosmischen Medium, und $u$ ist dann nicht zum Killing-Vektor der Sonne parallel (Abschnitt 4).

### 2.5 Das äußere Feld der Milchstraße und der Cassini-Quadrupol

[korr.] Dieser Abschnitt setzt voraus, dass die eigene Neigung $g(R)>a_0$ den Einzeltrichter kennzeichnet, dass er bis zur Gezeitengrenze Newtonsch bleibt und dass das homogene äußere Feld nicht in $\nu$ eingeht; das trägt nicht. Geglättet wird nur der Überschuss $(\rho-\rho_J)_+$ der Quelldichte, $\rho_J=3T/(4\pi G)$; $a_0$ kommt darin nicht vor. Die Glättungslänge der Sonne ist $\ell\approx0{,}93$ pc ($27\,r^*$) aus $\ell^3=GM/T$ mit $T=\sqrt{3/2}\,|E^{TF}|$ (voller spurfreier Gezeitentensor); 1,37 pc ist der klassische Jacobi-Radius. Der Quadrupol aus der Glättung ist $9\cdot10^{-31}$ s⁻². Die Regel sieht den Betrag des Gesamtfelds ($\Gamma=a_N/c$ enthält das Gesamtgefälle): Milchstraße bei 50 kpc im Feld $0{,}03\,a_0$: $\nu=6{,}50$ frei, 5,46 senkrecht, 4,68 gleichgerichtet.

Das Sonnensystem liegt im Newtonschen Feld der Milchstraße mit $g_e\approx1{,}2\cdot10^{-10}$ m/s², also bei $g_e\approx a_0$. Das Modell denkt die Sonne als **Einzeltrichter**: Das Kriterium ist die eigene Neigung des Körpers. Ein kompakter Trichter mit $g(R)>a_0$ bleibt bis zu seiner Gezeitengrenze im Elternfeld Newtonsch; erst dort geht er ins kollektive Feld der Milchstraße über. In $\nu$ gehen das eigene Feld und die Gezeiten der Umgebung ein, das homogene äußere Feld nicht (starkes Äquivalenzprinzip); das gilt nur unter der Annahme oben, im Modell geht das äußere Feld ein [korr.]. Der Quadrupol im Sonnensystem ist dann das galaktische Gezeitenfeld

$$
Q_2\sim\frac{g_e}{R_0}=4{,}7\cdot10^{-31}\ \text{s}^{-2},\qquad R_0\approx8{,}2\ \text{kpc (galaktozentrischer Abstand der Sonne)}.
$$

Das ist nur die Größenordnung des galaktischen Gezeitenfelds; aus der Glättung folgt $9\cdot10^{-31}$ s⁻² [korr.]. Die Cassini-Bahnverfolgung ergibt $Q_2=(3\pm3)\cdot10^{-27}$ s⁻² (Hees et al. 2014). Der Wert liegt weit unter dieser Schranke.

Nur unter der Annahme oben wäre diese Fassung nicht mehr QUMOND. Im Modell geht das äußere Feld richtungsabhängig ein (senkrecht $s^4=s_i^4+s_e^4$, gleichgerichtet $s^2=s_i^2+s_e^2$); prüfbar an Satelliten-Zwerggalaxien im Feld ihrer Muttergalaxie (Abschnitt 4). [korr.]

---

## 3. Linsenmasse gegen dynamische Masse

[korr.] Dieser Abschnitt setzt voraus, dass die Bahnen der Feldform folgen, also Dynamik $\Phi$ misst; das trägt nicht: Die Bahnen folgen der örtlichen Regel $a=\nu(|g_N|)\,g_N$, die Metrik trägt nur deren wirbelfreien Teil (an SPARC $f=+0{,}03$ bis $+0{,}13$, ohne IC 2574 $0{,}00\pm0{,}03$; die Feldform $f=1$ liegt außerhalb 95 %). Für Kugeln und außerhalb von Scheiben ist $\eta=M_\text{dyn}/M_\text{lens}=1$; in Scheiben innen ($R\lesssim3{,}5\,R_d$) ist $\eta>1$.

Dynamik misst $\Phi$, Lichtablenkung misst $(\Phi+\Psi)/2$. In FTM-AS wirkt der Operator auf $\Phi$ und $\Psi$ gleich, weil der Strukturteil nicht an die räumliche Metrik koppelt; $\delta\Psi$ erzwingt $\Psi=\Phi$ (1.5). Damit gilt für Kugeln und außerhalb von Scheiben [korr.]

$$
\eta=\frac{M_\text{dyn}}{M_\text{lens}}=1 .
$$

Im tiefen Regime gilt für den Nullpunkt der Linsen-Tully-Fisher-Relation $v_\text{lens}^4=GMa_0$, derselbe Wert wie in der Dynamik.

**Stand der Daten.**
- Brouwer et al. (2021, KiDS-1000) haben mit schwachen Linsen um isolierte Galaxien aus dem KiDS-bright-Katalog (~10⁶ Galaxien) die RAR etwa zwei Größenordnungen tiefer in $g_\text{bar}$ verfolgt. Die Linsen-RAR der spättypischen Galaxien folgt der RAR aus Rotationskurven. Frühtypische Galaxien liegen höher; wie weit, hängt von der unsicheren Masse des heißen Gases ab.
- Mistele, McGaugh, Lelli, Schombert & Li (2024) leiten aus denselben Daten Kreisgeschwindigkeiten ab, die bis ~1 Mpc flach bleiben. Ihre Linsen-Tully-Fisher-Relation stimmt mit der kinematischen überein.

**Vergleich:** Beide Befunde stimmen mit $\eta=1$ außerhalb der Scheiben überein, wo schwache Linsen messen [korr.]: Die Linsen-RAR der spättypischen Galaxien und die Linsen-Tully-Fisher-Relation folgen den kinematischen. $\Psi=\Phi$ ist außerdem die Fassung, in der Faltung und Takt aus derselben Metrik kommen (1.6) [korr.].

**Entscheidende Euclid-Messung.** Gemessen werden sollen die Linsen-RAR und die Linsen-Tully-Fisher-Relation isolierter Galaxien in Bins von $g_\text{bar}=10^{-11}$ bis $10^{-13}$ m/s², mit Baryonenmassen einschließlich des heißen Gases. Dazu kommt Dynamik derselben Stapel aus Satellitenkinematik bzw. HI-Rotationskurven, damit $\eta$ direkt bestimmt wird. Die Vorhersage lautet $\eta=1$ in jedem Bin, soweit die Linsen außerhalb der Scheiben messen; in Scheiben innen ($R\lesssim3{,}5\,R_d$) $\eta>1$ [korr.]. Die statistischen Fehler sinken gegenüber KiDS-1000 grob um das 8,2-Fache (Fläche × Quellendichte) [G]. Die Genauigkeit der $\eta$-Bestimmung wird damit durch die Systematik der Baryonenmassen begrenzt, nicht durch die Statistik.

---

## 4. Offene Punkte und nächste Schritte

1. **Einzeltrichter-Kriterium kovariant (vorrangig).** Nächster Schritt: eine kovariante Fassung von FTM-AS, in der $\nu$ vom eigenen Feld des Trichters abhängt und das homogene äußere Feld nicht eingeht (2.5): Kompakte Trichter mit $g(R)>a_0$ bleiben bis zur Gezeitengrenze im Elternfeld Newtonsch, erst dort gehen sie ins kollektive Feld über. [korr.] Gelöst durch die stetige Glättung: Geglättet wird nur der Überschuss $(\rho-\rho_J)_+$ über die Gezeitendichte $\rho_J=3T/(4\pi G)$, kovariant über ein geglättetes Faltungspotential $\bar\sigma$ mit der Gezeit aus der spurfreien Hesse-Form von $c^2\bar\sigma$; die Glättungslänge der Sonne ist $\ell\approx0{,}93$ pc. Das äußere Feld geht in $\nu$ ein.
   - Kandidat: Das Argument von $q$ wird relativ zum frei fallenden Bezugssystem des Systems gebildet.
   - Testfälle: $Q_2$ im Sonnensystem, Rotationskurven mit äußerem Feld (Zwerggalaxien nahe großer Galaxien), weite Doppelsterne bei $\sim r^*\approx7000$ AE. Die Datenlage bei den weiten Doppelsternen ist umstritten.
2. **Aether-Sektor wohldefiniert machen (offen [korr.]).** Für stabile Aether-Dynamik braucht es Einstein-Aether-Terme $c_1\ldots c_3$.
   - $c_{13}=0$ wegen GW170817.
   - $\alpha_1,\alpha_2$ im bewegten Sonnensystem berechnen; Schranken nach Foster & Jacobson 2006 und Oost, Mukohyama & Wang 2018.
   - Zu prüfen ist, ob die Kombination $B=a-D\sigma$ die Vorzugssystem-Effekte im Starkfeld unterdrückt.
3. **Instantanes $\sigma$ (offen [korr.]).** $\sigma$ hat im Ruhesystem des Mediums keine Zeitableitung, die Gleichung ist elliptisch wie in QUMOND. Ein Term $\mathcal K(\mathcal Q)$ mit $\mathcal Q=u^\mu\nabla_\mu\sigma$ (wie in AeST) gäbe $\sigma$ Dynamik. Für die lineare Kosmologie ist das nicht mehr Voraussetzung, die trägt die Schließungsregel (Punkt 4) [korr.].
4. **Kosmologie.** CMB-Spitzen und Materie-Leistungsspektrum sind nicht behandelt. Offene Kernfrage: Das CMB verlangt ein Potential, das nicht mit dem Photon-Baryon-Fluid mitschwingt. Kann das Modell es als Einschränkung der Freiheitsgrade bzw. als Gedächtnis des Operators liefern, nicht als zusätzliche Substanz? [korr.] Im Modell geschlossen durch die Schließungsregel $(\partial\xi)^2=-1$, eine Einschränkung der Freiheitsgrade (die Druckrichtung der konformen Mode ist geschlossen, kein Stoff): Linear ist sie gleich dem Standardmodell bei $\Omega_c/\Omega_b=Q=5{,}36$ (aus der Messung eingesetzt); Spitzen und Strukturwachstum folgen.
5. **$a_0=cH_0/2\pi$ dynamisch.** [korr.] Entfällt: $a_0\propto H$ trägt nicht. Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie $\Delta N_\text{eff}=2{,}09$ (gemessen $N_\text{eff}=2{,}99\pm0{,}17$); mit den Schleifen atmend ($\Lambda\propto a$) wirkte er wie eine Krümmung von 0,22 (gemessen $\Omega_K=0{,}0007\pm0{,}0037$). Im Modell ist $a_0=c^2/(2\Lambda)$, $\Lambda$ die ganze innere Faltung der Zelle und heute fest; $a_0$ ist seit $z\approx10$ fest (bei der Rekombination $0{,}76\,a_0$), $\approx cH_0/2\pi$ gilt nur heute. Die Folge $a_0(z)\propto H(z)$ gibt es nicht.
6. **Vakuumterm.** Klären, ob $\Lambda_s=(\pi^2/15)H_0^2/c^2$ (0,320 $\Lambda_\text{kosm}$) physikalisch ist, also keine nackte Konstante $\Lambda_\text{kosm}$ auftritt, oder nur Normierung. Im ersten Fall ist die räumliche Verteilung dieses Terms (voll in Leerräumen, null in Galaxien) eine eigene Vorhersage. [korr.] Im Modell ist $\Lambda_s=\pi^4/(15\Lambda^2)$ der Vakuumteil der Strukturwirkung, die Energie der Rückstellung, heute ein Drittel der Dunklen Energie ($\Omega_s=\pi^2/45$ bei $\Theta=1$); er hängt nur an $\Lambda$ und tritt nur global auf, nicht je Trichter; eine räumliche Verteilung (voll in Leerräumen, null in Galaxien) hat er nicht.
7. **Galaxienhaufen.** QUMOND-artige Theorien behalten in Haufen einen Massenfehlbetrag von etwa Faktor 2 [G]. FTM-AS erbt das; offen ist, ob ein kollektiver Haufenterm ihn schließt. [korr.] Im Modell geschlossen durch die Schließungsregel: Im Haufen fasst die Operatorform nur einen Teil der Mode ξ, der stoßfreie Rest schließt den fehlenden Faktor und trifft die Linsenmasse auf 1–5 %.
8. **Integrodifferentialgleichung.** Klären, ob sie nur als Bahn-Näherung von FTM-AS behalten wird (1.7). Wenn ja, dann mit $M_k=q$ und $\kappa_0=\lambda_m^2$ im Kurzgedächtnis-Grenzfall; $\Delta(\tau_M)$ bliebe Nahfeld-Form. [korr.] Als Geschwindigkeitskern auf jeder Bahn ist sie ausgeschlossen (er würde exzentrische Bahnen runden; das Gedächtnis wirkt auf den Zustand der Regel); im Modell ist die Rückkehr stationär und $\kappa_0/\lambda_\text{eff}=3c^2/(r^*v_\text{rot})$.

---

## Quellen (Vorbilder und Daten)

- Milgrom, M. (2009), *Phys. Rev. D* 80, 123536 (BIMOND); Milgrom, M. (2010), *MNRAS* 403, 886 (QUMOND)
- Bekenstein, J. & Milgrom, M. (1984), *ApJ* 286, 7; Bekenstein, J. (2004), *Phys. Rev. D* 70, 083509 (TeVeS)
- Zlosnik, T., Ferreira, P. & Starkman, G. (2007), *Phys. Rev. D* 75, 044017
- Skordis, C. & Złośnik, T. (2021), *Phys. Rev. Lett.* 127, 161302 (AeST)
- Foster, B. & Jacobson, T. (2006), *Phys. Rev. D* 73, 064015; Oost, J., Mukohyama, S. & Wang, A. (2018), *Phys. Rev. D* 97, 124023
- Bertotti, B., Iess, L. & Tortora, P. (2003), *Nature* 425, 374 (Cassini)
- Everitt, C. W. F. et al. (2011), *Phys. Rev. Lett.* 106, 221101 (Gravity Probe B)
- Hees, A. et al. (2014), *Phys. Rev. D* 89, 102002; Hees, A. et al. (2016), *MNRAS* 455, 449; Desmond, H., Hees, A. & Famaey, B. (2024), *MNRAS*
- Brouwer, M. et al. (2021), *A&A* 650, A113 (KiDS-1000, Linsen-RAR); Mistele, T. et al. (2024), *ApJL* 969, L3
- Feynman, R. P., *Lectures on Physics* II, Kap. 42 (excess radius)
