# Die Schließungsregel: Projektion der Fluchtrichtung mit Gewicht w(θ_u)

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: a₀/c als „Atemfrequenz H₀/2π des Torus“ (Kopf, 1.1) → trägt nicht (a₀ ∝ H gäbe ΔN_eff = 2,09, gemessen N_eff = 2,99 ± 0,17); a₀ = c²/2Λ ist seit z ≈ 10 fest, ≈ cH₀/2π gilt nur heute. N* bis 12,5 bzw. 20,3 als „obere Grenze“ → Fenster N* ≈ 0,3–2,4 (spät kollabierte Galaxien verlangen N* ≤ 2,38); Werte mit N* = 20 bzw. 20,3 (fσ₈, k_h = 0,62/1,35/2,87/11,5 h/Mpc, +21 % in T5) als außerhalb des Fensters markiert. 1 − w bei der Rekombination 9,2·10⁻¹² (bei 10σ 9,3·10⁻¹²) → 7,0·10⁻¹² (bei 10σ 7,0·10⁻¹²) bei festem a₀ = 1,2·10⁻¹⁰, 4·10⁻¹² mit a₀(1090) = 0,76·a₀; bei z = 100 1,2·10⁻⁸ → 1,12·10⁻⁸; Ursache: H(z) ohne Strahlung gerechnet. Bullet Cluster: in 150 Myr höchstens 33 % → 44 % (1/λ_m ≥ 0,26 Gyr). Kurzfassung: R_Q = 50,9 kpc → 49 kpc (Punktmasse, a₀ = 1,2·10⁻¹⁰; das Hernquist-Profil gibt mit dem Kern 50,9 kpc); λ_m = ω_vir/N* → ω_orb/N*; N* ≤ 2,4 bis 12,5 → N* ≈ 0,3–2,4. Wortlaut: ρ_ph, M_ph → gelesene Dichte ρ_X bzw. gelesene Masse M_X der Operatorform; ρ_B, T_B^{μν} als Lesung des übergegangenen Anteils, kein eigener Sektor; „Kapazität ρ_B ≤ ρ_ph“ → die Operatorform fasst (ν − 1)·M_b; Halo-Profil, Halo-Überschuss, Halo-Streuung und „isotherme Kugel“ der Mode umformuliert. [G] an Δ_ζ² und n_s. Umbenannt (ohne Marke): χ (Gedächtnismode) → ξ, (∂χ)² → (∂ξ)²; ξ_R, ξ_rest (Restanteil) → R_ξ; s (Rückwirkungsfaktor) → k_R; Ω_Λ → Ω_DE; ρ_ph → ρ_X, M_ph → M_X.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors (Vorschlag A13).
> Rechnung: rechnungen/schliessungsregel.py
> Literatur- und Beobachtungswerte mit **[G]** sind aus dem Gedächtnis zitiert und nicht nachgeprüft.
> Zahlenwert der Schwelle: a₀ = 1,2·10⁻¹⁰ m/s² (gemessen), also θ_c = a₀/c = 4,00·10⁻¹⁹ s⁻¹. a₀/c ist keine Atemfrequenz: a₀ = c²/2Λ ist seit z ≈ 10 fest, a₀ ≈ cH₀/2π gilt nur heute. [korr.]

---

## Lesart

Es kommt kein Stoff hinzu. ξ ist die **konforme Mode der Metrik**, die durch die Nebenbedingung (∂ξ)² = −1 eingeschränkt wird. Die Druckrichtung ist geschlossen, übrig bleibt eine Mode ohne eigene Dynamik quer zur Bahn. Diese Mode hat keine eigenen Anfangsdaten. Sie startet im Zustand der Baryonen und trägt deren Bewegungsgeschichte weiter. Sie ist damit die **Gedächtniskoordinate** des Bindungsoperators: Die Vergangenheit der Baryonenbewegung steckt in einer Koordinate des gegenwärtigen Zustands.

Der Lagrange-Multiplikator μ der Einschränkung wirkt im Energie-Impuls-Tensor wie druckfreier Staub. Ein Beobachter, der ungeänderte Regeln annimmt, rechnet diese Wirkung per Poisson-Inversion in eine Dichte um und erhält Q·ρ_b. Diese Dichte ist Ausgabe der Inversion, nicht Eingabe der Dynamik. „Scheinbare Masse“ und „scheinbare Dichte“ bezeichnen in dieser Notiz genau diese Umrechnung.

Die Gedächtnismode trägt die Geschichte der Baryonenbahnen. Ihre Wirkung liegt deshalb dort, wo die Baryonen **waren**, nicht dort, wo sie jetzt sind. Das ist der pfadabhängige Gedächtniseffekt: Im Bullet Cluster folgt sie den Bahnen vor dem Stoß (stoßfrei), während das Gas durch den Stoß abgebremst ist.

**Offen.** Q = Ω_c/Ω_b ist hier eingesetzt (Planck), nicht hergeleitet. Die Peaklage der CMB ist deshalb keine unabhängige Bestätigung: Im linearen Bereich ist die Schließungsregel mit ω_c = Q·ω_b gleich ΛCDM.

---

## Kurzfassung

| Frage | Ergebnis |
|---|---|
| Wirkung | Die Gedächtniskoordinate ist ein Potential ξ mit der Nebenbedingung $g^{\mu\nu}\partial_\mu\xi\partial_\nu\xi=-1$. Diese Einschränkung schließt die Druckrichtung ($c_s^2=0$ exakt). Das Gewicht $w(\theta_u)$ steuert zwei Dinge: die **Mitnahme** durch die Druckkraft (Anteil $1-w$) und die **Relaxation in die Operatorform**, sobald ein Bereich nicht mehr mitatmet. |
| Wahl von w | $w=\theta_+^2/(\theta_+^2+\theta_c^2)$ mit $\theta_c=a_0/c=0{,}183\,H_0$. Im Hintergrund gilt $1-w\le3{,}7\cdot10^{-3}$ (heute) und $7{,}0\cdot10^{-12}$ bei $z_*$ (Strahlung in $H(z)$; mit $a_0(1090)=0{,}76\,a_0$ sind es $4\cdot10^{-12}$). [korr.] Beim Kollaps fällt $w$ innerhalb von 7 Myr vor dem Umkehrpunkt von 0,9 auf 0,1. |
| Erhaltung | $\nabla_\mu T^{\mu\nu}_\text{ges}=0$ gilt nach Bauart. Die 0i-Gleichung ist im Boltzmann-Spielzeug für jedes $w$ auf $\le5\cdot10^{-6}$ (Median) erfüllt. Nötig ist dafür die Rückwirkung der Mitnahme auf das Photon-Baryon-Fluid. |
| Stabilität | $c_s^2=1/[3(1+R+(1-w)QR)]\ge0$ für alle $w\in[0,1]$. Die freie Gedächtnismode verhält sich druckfrei mit $\mu\ge0$, die Relaxationsrate ist $\ge0$. |
| CMB | $w(\theta_u)$ ergibt dasselbe Peakmuster wie $w=1$ (Abweichung 0,000 %). |
| Wachstum | $f\sigma_8$ weicht um höchstens $5\cdot10^{-4}$ ab (N* = 1) bzw. $2\cdot10^{-5}$ (N* = 20, außerhalb des Fensters N* ≈ 0,3–2,4 [korr.]). Der Hintergrund bleibt exakt staubartig. |
| Galaxie | Die Energie der Gedächtnismode bleibt als **gebundene Energie** erhalten: Sie geht in die Strukturenergie des Operators über. Die Bilanz $QM_b=M_X(<R_Q)$ (gelesene Masse der Operatorform) ist exakt am Radius $R_Q=r^*/(-\ln q)$, an dem das lokale $q$ gleich dem kosmischen ist. Für $M_b=6\cdot10^{10}\,M_\odot$ ist $R_Q=49$ kpc ($r^*=8{,}35$ kpc, $q=0{,}843$; als Hernquist-Profil 50,9 kpc). [korr.] Frei fallende Gedächtnismode und Operatorform unterscheiden sich im Mittel nur um 4,8 %. Danach gilt $g/g_\text{RAR}\to1$. |
| Haufen | Die Operatorform fasst in einem Haufen nur 21–71 % der Gedächtnismode (100–1500 kpc). Der Rest bleibt als freie, stoßfreie Gedächtnismode, und genau sie schließt den fehlenden Faktor 1,2–3,0: $(M_b+M_P)/M_\text{ref}=0{,}99$–$1{,}05$. Über $M_{200}\approx3{,}7\cdot10^{14}\,M_\odot$ (heute) bleibt Rest-Gedächtnis, darunter nicht. |
| Bullet Cluster | In 100 kpc um den Haupt-BCG sind 79 % der scheinbaren Masse stoßfrei und liegen bei den Galaxien. Der Rest (21 %) folgt den Baryonen mit Nachlauf. |
| Nachlaufzeit | Relaxation mit $\lambda_m=\omega_\text{orb}/N^*$. [korr.] Aus der Bullet-Schranke folgt $N^*\ge0{,}28$, aus der Enge der RAR auch für spät kollabierte Galaxien $N^*\le2{,}38$: Fenster $N^*\approx0{,}3$–2,4. [korr.] |

---

## 1. Wirkung

### 1.1 Felder

Zu den Feldern von FTM-AS (Metrik $g_{\mu\nu}$, Zeitfeld $u^\mu$ mit $u_\mu u^\mu=-1$, Faltungspotential σ, Strukturterm mit $\nu=1/(1-q)$) treten die Größen der Einschränkung. Es sind keine neuen Stoffe, sondern eine eingeschränkte Mode der Metrik, ihr Multiplikator und das Atmen:

| Feld | Bedeutung |
|---|---|
| ξ | Gedächtniskoordinate: die durch (∂ξ)² = −1 eingeschränkte konforme Mode der Metrik. Sie hat keine eigenen Anfangsdaten, sondern startet im Zustand der Baryonen und trägt deren Bewegungsgeschichte. |
| μ | Lagrange-Multiplikator der Einschränkung. Ein Beobachter, der ungeänderte Regeln annimmt, liest ihn als Dichte $\rho_F$ (Poisson-Inversion). |
| $\rho_B$ | Lesung des in die Operatorform übergegangenen Anteils der Mode (Strukturenergie des Operators, in Ruhe im Zeitfeld $u$), gelesen als Dichte; kein eigener Sektor [korr.] |
| $\theta_u=\nabla_\mu u^\mu$ | Expansionsrate des Zeitfelds (das Atmen) |

Das Gewicht:

$$
w(\theta_u)=\frac{\theta_+^2}{\theta_+^2+\theta_c^2},\qquad \theta_+=\max(\theta_u,0),\qquad \theta_c=\frac{a_0}{c}=4{,}00\cdot10^{-19}\ \text{s}^{-1}=0{,}183\,H_0 .
$$

**Warum diese Form.**
- $a_0/c$ ist die einzige Rate, die das Modell selbst mitbringt. Eine Atemfrequenz ist sie nicht: Atmete $a_0$ mit $H$ ($a_0=cH/2\pi$), bliebe der Vakuumteil zu jeder Zeit bei 22 %, im Strahlungszeitalter wie $\Delta N_\text{eff}=2{,}09$ (gemessen $N_\text{eff}=2{,}99\pm0{,}17$). Im Modell ist $a_0=c^2/2\Lambda$ mit der Schließlänge Λ als ganzer innerer Faltung der Zelle heute fest; $a_0$ ist seit $z\approx10$ fest (bei der Rekombination $0{,}76\,a_0$), $a_0\approx cH_0/2\pi$ gilt nur heute. [korr.]
- Jeder Bereich, der mit dem Hintergrund expandiert, hat $\theta_u\ge3H\gg\theta_c$ und damit $w\approx1$, zu jeder Epoche bis heute.
- $w$ fällt erst, wenn die lokale Expansion zum Stillstand kommt (Umkehrpunkt). Kontraktion ($\theta_u<0$) ist kein Mitatmen, daher $w=0$.
- Ein stationäres gebundenes System hat das Zeitfeld entlang des Killing-Vektors. Dort gilt $\theta_u=0$ exakt, also $w=0$.
- Die Form ist $C^1$-glatt. Das Ergebnis hängt praktisch nicht vom genauen Wert von $\theta_c$ ab, solange $\theta_c\ll H$ gilt: Der Übergang wird dann zur Stufe am Umkehrpunkt.

### 1.2 Reversibler Teil: freie Gedächtnismode als eingeschränkter Freiheitsgrad

$$
S=S_\text{FTM-AS}[g,u,\sigma]+S_m[g,\psi]-\int d^4x\,\sqrt{-g}\;\frac{\mu c^2}{2}\big(g^{\mu\nu}\partial_\mu\xi\,\partial_\nu\xi+1\big).
$$

Feldgleichungen:

| Variation | Gleichung | Bedeutung |
|---|---|---|
| δμ | $g^{\mu\nu}\partial_\mu\xi\partial_\nu\xi=-1$ | $U^\mu=-\nabla^\mu\xi$ ist ein Einheitsvektor. Es folgt $U^\nu\nabla_\nu U_\mu=U^\nu\nabla_\mu U_\nu=\tfrac12\nabla_\mu(U\cdot U)=0$: Die Mode fällt geodätisch. |
| δξ | $\nabla_\mu(\mu U^\mu)=0$ | Kontinuität |
| δg | $T_F^{\mu\nu}=\mu\,U^\mu U^\nu$ | druckfrei, $c_s^2=0$ exakt |

Die Nebenbedingung ist die **Schließung der Fluchtrichtung**. Eine Mode mit Einheitsgradient hat keinen Gradiententerm, kann also keinen Druck tragen. Gravitation und Trägheit bleiben erhalten (Form: mimetischer Staub, Chamseddine & Mukhanov 2013 [G]; Staub aus Nebenbedingung, Lim, Sawicki & Vikman 2010 [G]).

### 1.3 Mitnahme bei w < 1

Bei $w<1$ wird der Anteil $1-w$ der Gedächtnismode von der nichtgravitativen Kraft des Photon-Baryon-Fluids (Druck, Thomson-Kopplung) mitgenommen. Die Kraft wirkt dann auf die gemeinsame Trägheit:

$$
\nabla_\mu T_F^{\mu\nu}=J_p^\nu,\qquad \nabla_\mu T_{\gamma b}^{\mu\nu}=-J_p^\nu,\qquad
J_p^\nu=(1-w)\,\rho_F\,\frac{h^\nu{}_\alpha F^\alpha_\text{ng}}{I_{\gamma b}+(1-w)\rho_F},\quad I_{\gamma b}=\rho_b+\tfrac43\rho_\gamma .
$$

Linear in Newton-Eichung ($'=d/d\ln a$, $\mathcal H$ konformes Hubble-Maß) folgt:

$$
\theta_{\gamma b}'=-\theta_{\gamma b}+\frac{k^2\phi}{\mathcal H}+k_R\,f,\qquad
\theta_F'=-\theta_F+\frac{k^2\phi}{\mathcal H}+(1-w)\,k_R\,f,\qquad
\delta_F'=-\frac{\theta_F}{\mathcal H}+3\phi',
$$

$$
f=\frac{\theta_{\gamma b}}{1+R}+\frac{k^2\delta_\gamma}{4(1+R)\mathcal H},\qquad k_R=\Big[1+\frac{(1-w)QR}{1+R}\Big]^{-1}.
$$

Linear ist das gleichbedeutend damit, dass der Anteil $1-w$ mit dem Fluid mitläuft und der Anteil $w$ frei fällt.

### 1.4 Irreversibler Teil: Relaxation in die Operatorform

Wo $w<1$ ist, verliert die Gedächtnismode seinen freien Fallfreiheitsgrad. Er relaxiert in die Form, die der Operator verlangt. Gelesen als Dichte (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt; kein Stoff) ist das [korr.]:

$$
\rho_X=\frac{c^2}{4\pi G}\,D_\mu\big[(\nu-1)\,D^\mu\sigma\big]\;\xrightarrow{\text{schwaches Feld}}\;\frac{\nabla\cdot[(\nu-1)\nabla\Phi_N]}{4\pi G}.
$$

Die Relaxation geht nur so weit, wie die Operatorform die Mode fassen kann: Sie fasst $(\nu-1)\,M_b$, gelesen $\rho_B\le\rho_X$. In gelesenen Dichten geschrieben [korr.]:

$$
\nabla_\mu(\rho_F U^\mu)=-\Gamma\rho_F,\qquad \nabla_\mu T_B^{\mu\nu}=+\Gamma\rho_F U^\nu,\qquad
\Gamma=\frac{\omega_\text{orb}}{N^*}\,\big(1-w(\theta_u)\big),\qquad \omega_\text{orb}^2=\tfrac13 R_{\mu\nu}u^\mu u^\nu=\tfrac{4\pi G}{3}\bar\rho .
$$

$\omega_\text{orb}$ ist die Bahnfrequenz der Zone, gemittelt über die Zone. $N^*$ ist die Zahl der Radiant Bahnphase pro e-Faltung der Relaxation und der einzige neue Parameter. $\Gamma\ge0$ ist der entropieerzeugende (transversale) Teil einer GENERIC-Struktur (Öttinger & Grmela 1997 [G]). Der gebundene Anteil ist kein zweiter Stoff neben dem Operator, sondern dessen Strukturenergie. $\rho_F$ und $\rho_B$ (mit $T_B^{\mu\nu}$) sind gelesene Dichten, kein eigener Sektor: die Umrechnung der Mode ξ (frei) und der Regel (übergegangen) in Dichten; der Austausch $\Gamma\rho_F$ ist die Übergabe der Mode in die Regel, als Bilanz dieser Lesung geschrieben. [korr.] Im stationären Zustand ist $T_B^{00}=\rho_X c^2$ genau der Strukturterm von FTM-AS.

### 1.5 Die Anforderungen (a) bis (e)

**(a) Feldgleichungen:** 1.2 bis 1.4.

**(b) Erhaltung.** Jeder Austausch ($J_p$, $\Gamma\rho_F U$) tritt mit entgegengesetztem Vorzeichen in zwei Bilanzen auf. [korr.] Daher gilt $\nabla_\mu T^{\mu\nu}_\text{ges}=0$, verträglich mit der Bianchi-Identität.
- **00:** $\dot\rho_F+3H\rho_F=-\Gamma\rho_F$ und $\dot\rho_B+3H\rho_B=+\Gamma\rho_F$, also $\rho_b+\rho_F+\rho_B\propto a^{-3}$ exakt; das ist die Bilanz der gelesenen Dichten: Die Mode geht in die Regel über, kein Stoff wandert zwischen Sektoren. [korr.]
- **0i:** Im Spielzeug-Boltzmann treibt die 00-Gleichung φ. Die 0i-Gleichung $k^2(\phi'+\phi)=\tfrac32\mathcal H\sum_i\Omega_i(1+w_i)\theta_i$ dient als unabhängige Probe. Relative Abweichung (Median/Maximum über $0{,}01<k<0{,}3$/Mpc, bei $a=10^{-4}$ bis $a_*$):

| | Median | Maximum |
|---|---|---|
| w = 1 | 1,0–4,7·10⁻⁶ | ≤ 3,9·10⁻⁵ |
| w = 0,9 mit Rückwirkung | 1,0–3,7·10⁻⁶ | ≤ 3,7·10⁻⁵ |
| w = 0,5 mit Rückwirkung | 1,0–3,5·10⁻⁶ | ≤ 1,2·10⁻⁴ |
| w = 0,9 ohne Rückwirkung auf das Fluid | bis 5,2·10⁻³ | bis 5,0·10⁻² |

Die Rückwirkung (Faktor $k_R$ in der Fluidgleichung) ist also Teil der Regel.

**(c) Stabilität.**
- Schallgeschwindigkeit bei $z_*$ ($R=0{,}622$): $c_s^2=0{,}2055$ (w = 1), 0,2014 (0,99), 0,1705 (0,9), 0,1013 (0,5), 0,0672 (0). Sie ist nie negativ, es gibt also keine Gradienteninstabilität.
- Die freie Gedächtnismode ist druckfrei: kein Gradiententerm, $c_s^2=0$, Energiedichte $\mu\ge0$, also kein Geist.
- Die Relaxation hat $\Gamma\ge0$ und hält beide Dichten positiv.

**(d) Grenzfall w = 1.** Hier ist $J_p=0$ und $\Gamma=0$. Die Gleichungen sind Zeile für Zeile die einer druckfreien Mode mit der scheinbaren Dichte $Q\rho_b$ mit $Q=5{,}364$. Das Peakmuster ist identisch (Abschnitt 2).

**(e) Wohin geht die Energie bei w → 0?** Sie bleibt **als gebundene Energie am Ort**: Die Gedächtnismode geht in die Strukturenergie des Operators über. Zwei Zahlen tragen das:

1. **Integrale Bilanz.** Kugelsymmetrisch gilt $M_\text{dyn}(<r)=\nu M_b(<r)$, also ist die gelesene Masse der Operatorform $M_X(<r)=(\nu-1)M_b(<r)$. [korr.] Die scheinbare Masse der Gedächtnismode der Zone ist $QM_b$. Gleichheit gilt genau bei $\nu=1+Q$, also $q_\text{lok}=q_\text{kos}=0{,}84287$. Das ist der Radius
$$
R_Q=\frac{r^*}{-\ln q}=\frac{r^*}{0{,}1709},\qquad g_N(R_Q)=(\ln q)^2a_0=0{,}0292\,a_0\quad\text{(Punktmasse)}.
$$
Für $M_b=6\cdot10^{10}\,M_\odot$ ($r^*=8{,}35$ kpc) ist $R_Q=48{,}8$ kpc als Punktmasse und $50{,}9$ kpc als Hernquist-Profil ($g_N=2{,}9\cdot10^{-12}$ m/s²).
2. **Form.** Die frei fallende Gedächtnismode (angesetzt mit gelesener Masse ∝ r bis $r_\text{vir}=R_Q$, der Form einer isothermen Kugel; kein Stoff) [korr.] und die Operatorform liegen nahe beieinander, weil $Q\cdot(-\ln q)=0{,}917$ ist (Galaxie $6\cdot10^{10}\,M_\odot$):

| r [kpc] | r/r* | $M_X/M_b$ | $M_{P,\text{frei}}/M_b$ | Verhältnis |
|---|---|---|---|---|
| 5 | 0,60 | 0,217 | 0,530 | 0,41 |
| 10 | 1,20 | 0,671 | 1,059 | 0,63 |
| 20 | 2,40 | 1,754 | 2,119 | 0,83 |
| 30 | 3,59 | 2,904 | 3,178 | 0,91 |
| 40 | 4,79 | 4,075 | 4,237 | 0,96 |
| 50,6 | 6,07 | 5,332 | 5,364 | 0,99 |

Umzuverteilen ist im Mittel $\int_0^{R_Q}|M_X-M_{P,\text{frei}}|\,dr/R_Q=4{,}8\,\%$ der scheinbare Masse der Gedächtnismode, vor allem innerhalb von $2r^*$. Lokal fließt die Differenz über das Strukturfeld: σ ist auf dem Ruheschnitt elliptisch, also nichtlokal. Die Bilanz schließt über die Zone.

In Haufen fasst die Operatorform mit $(\nu-1)M_b$ nicht die ganze Mode (Abschnitt 4). [korr.] Der Rest bleibt freie Gedächtnismode, auch das ist energieerhaltend.

---

## 2. Kosmos

### 2.1 Das Gewicht im Hintergrund, θ = 3H

| z | 0 | 0,5 | 1 | 2 | 10 | 100 | 1090 |
|---|---|---|---|---|---|---|---|
| 1 − w | 3,7·10⁻³ | 2,1·10⁻³ | 1,2·10⁻³ | 4,1·10⁻⁴ | 8,9·10⁻⁶ | 1,12·10⁻⁸ [korr.] | 7,0·10⁻¹² [korr.] |

$H(z)$ mit Strahlung, bei festem $a_0=1{,}2\cdot10^{-10}$; mit dem Modellwert $a_0(1090)=0{,}76\,a_0$ ist $1-w$ bei der Rekombination $4\cdot10^{-12}$. [korr.]

Die Moden bei $z_*$ atmen mit: $\theta_{u,\text{lok}}=3H(1-\delta'/3)$ mit $\max|\delta'|/3=6{,}2\cdot10^{-4}$ (pro ln k, normiert auf $\Delta_\zeta^2=2{,}1\cdot10^{-9}$ [G]). Auch bei 10σ ist $1-w(z_*)=7{,}0\cdot10^{-12}$. [korr.]

### 2.2 Peakmuster (Spielzeug-Boltzmann, Verhältnisse der Peakhöhen im k-Raum)

| Fall | H2/H1 | H3/H1 | H3/H2 | H4/H3 | H5/H4 | größte Abweichung |
|---|---|---|---|---|---|---|
| w = 1 | 0,708 | 1,151 | 1,625 | 0,597 | 0,856 | – |
| **w(θ), θ_c = a₀/c** | 0,708 | 1,151 | 1,625 | 0,597 | 0,856 | **0,000 %** |
| w = 0,99 konstant | 0,679 | 1,117 | 1,644 | 0,588 | 0,857 | 4,1 % |
| w = 0,97 konstant | 0,626 | 1,051 | 1,680 | 0,570 | 0,859 | 11,7 % |

Das Muster verlangt $w\ge0{,}99$ bei der Rekombination für eine Abweichung von höchstens wenigen Prozent. $w(\theta_u)$ liegt dort um $4$–$7\cdot10^{-12}$ unter 1. [korr.]

### 2.3 Hintergrund und lineares Wachstum

- **Staubartig.** Gedächtnismode (frei und in die Regel übergegangen, gelesen als $\rho_F$ und $\rho_B$) und Baryonen sind druckfrei, und die Übergabe erhält die Summe. [korr.] Daher gilt $a^3(\rho_b+\rho_F+\rho_B)=\text{const}$ exakt, mit $\Omega_m=\Omega_b(1+Q)=0{,}3138$.
- **Relaxierter Anteil im Hintergrund** (Rate $(1-w)\omega_\text{orb}/N^*$): $\varepsilon=1{,}1\cdot10^{-3}/N^*$ heute und $2{,}9\cdot10^{-4}/N^*$ bei $z=1$.
- **fσ₈.** Obere Schranke: Der relaxierte Anteil wird so behandelt, als klumpe er nicht; $\sigma_8=0{,}811$ [G].

| z | 0 | 0,3 | 0,5 | 0,8 | 1,0 | 1,5 | 2,0 |
|---|---|---|---|---|---|---|---|
| fσ₈ (ΛCDM) | 0,4265 | 0,4729 | 0,4742 | 0,4527 | 0,4314 | 0,3745 | 0,3245 |
| rel. Abweichung, N* = 1 | −4,7·10⁻⁴ | −3,0·10⁻⁴ | −2,2·10⁻⁴ | −1,4·10⁻⁴ | −1,1·10⁻⁴ | −6·10⁻⁵ | −4·10⁻⁵ |
| rel. Abweichung, N* = 20,3 (außerhalb des Fensters N* ≈ 0,3–2,4 [korr.]) | −2,3·10⁻⁵ | | −1,1·10⁻⁵ | | −5,4·10⁻⁶ | | −1,8·10⁻⁶ |

### 2.4 Wo die erste Abweichung einsetzt (prüfbar)

Die Regel greift erst, wenn ein Bereich den Umkehrpunkt erreicht ($\theta_{u,\text{lok}}=0$). Als Maß dient: 2σ-Bereiche erreichen die Umkehr, also $2\sigma(R_L,z)=\delta_\text{ta,lin}=1{,}076$ (aus dem Top-Hat). Das betrifft 4,6 % der Masse (Press-Schechter). Transfer nach BBKS mit Sugiyama-Form [G], $n_s=0{,}965$ [G].

| z | 0 | 0,5 | 1 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|
| $k_1=1/R_L$ [h/Mpc] | 0,072 | 0,101 | 0,143 | 0,267 | 0,462 | 1,20 |
| M [M☉/h] | 9,8·10¹⁴ | 3,5·10¹⁴ | 1,3·10¹⁴ | 1,9·10¹³ | 3,7·10¹² | 2,1·10¹¹ |

Die Relaxation ist vollzogen ($R_\xi=1/e$ für 2σ-Peaks, $\delta_c=1{,}703$, $N^*=20{,}3$) ab $k_h=0{,}62$ h/Mpc (z = 0), 1,35 (0,5), 2,87 (1) und 11,5 (2). Diese Werte gelten für $N^*=20{,}3$ und liegen außerhalb des Fensters $N^*\approx0{,}3$–2,4 (spät kollabierte Galaxien verlangen $N^*\le2{,}38$); für das Fenster ist $k_h$ nicht gerechnet. [korr.] Zwischen $k_1$ und $k_h$ liegt das Übergangsregime (Gruppen, Haufen). Beide Maßstäbe sind Lagrange-Maßstäbe; im Euler-Raum erscheint die Struktur bei etwa 5- bis 6-fach höherem k. Oberhalb von $k_1(z)$ ist das lineare Wachstum wie im Standardmodell.

---

## 3. Galaxie: der Übergang w → 0 beim Kollaps

**Ansatz.** Baryonen $M_b=6\cdot10^{10}\,M_\odot$ (Hernquist, $a=2{,}57$ kpc), scheinbare Masse der Gedächtnismode $QM_b=3{,}22\cdot10^{11}$, Top-Hat mit $(1+Q)M_b=3{,}82\cdot10^{11}\,M_\odot$ im Hintergrund $\Omega_m=0{,}314$, $\Omega_\text{DE}$. Kontrolle: Lineares δ bei der Umkehr ist 1,076 und beim Kollaps 1,703 ($z_c=3$). $\theta_{u,\text{lok}}=3\dot R/R$; nach dem Kollaps gilt $r_\text{vir}=r_\text{ta}/2$ und θ = 0.

| $z_c$ | $z_\text{ta}$ | $r_\text{ta}$ [kpc] | $r_\text{vir}$ [kpc] | $r_\text{vir}/R_Q$ | $\omega_\text{vir}$ [1/Gyr] |
|---|---|---|---|---|---|
| 0,5 | 1,49 | 292 | 146 | 2,87 | 0,74 |
| 1 | 2,24 | 228 | 114 | 2,23 | 1,08 |
| 2 | 3,79 | 155 | 77,6 | 1,52 | 1,92 |
| 3 | 5,36 | 117 | 58,5 | 1,15 | 2,93 |
| **3,63** | 6,36 | 101 | **50,6** | **1** | **3,64** |
| 5 | 8,52 | 78 | 39,1 | 0,77 | 5,36 |

Die Zone, deren Gedächtnismode die Operatorstruktur genau füllt, entspricht einem Kollaps bei $z_c=3{,}63$ ($t_c=1{,}73$ Gyr).

**Verlauf von θ und w** (Kollaps bei $z_c=3{,}63$):

| t [Gyr] | z | R/R_ta | θ_lok/H₀ | θ_bg/H₀ | w |
|---|---|---|---|---|---|
| 0,173 | 20,5 | 0,53 | 141,5 | 167,8 | 1,0000 |
| 0,692 | 7,54 | 0,975 | 12,9 | 42,0 | 0,9998 |
| 0,822 | 6,61 | 0,999 | 3,07 | 35,4 | 0,9964 |
| 0,856 | 6,41 | 1,000 | 0,62 | 34,0 | 0,918 |
| 0,864 | 6,36 | 1,000 | −0,02 | 33,6 | 0 |
| 1,47 | 4,16 | 0,66 | −87 | 19,9 | 0 |

$w$ fällt zwischen 0,857 und 0,864 Gyr von 0,9 auf 0,1, also innerhalb von 7 Myr vor dem Umkehrpunkt. Danach bleibt $w=0$: im Einfall, weil θ < 0 ist, und nach der Virialisierung, weil θ = 0 ist.

**Dynamik in der Galaxie.** Das Rest-Gedächtnis $R_\xi$ klingt mit $\exp[-\int(1-w)\omega_\text{orb}\,dt/N^*]$ ab. Die Abweichung ist linear in $R_\xi$:

$$
\frac{g}{g_\text{RAR}}=1-R_\xi\,\frac{\nu-1}{\nu}+R_\xi\,\frac{Q\,M_b\min(r/r_\text{vir},1)}{\nu\,M_b(<r)} .
$$

Für $|g/g_\text{RAR}-1|\le5\,\%$ bei 5–50 kpc genügt $R_\xi\le0{,}104$. Der Exponent erreicht heute 45,8 (für $z_c=3{,}63$). Im Fenster $N^*\approx0{,}3$–2,4 ist die Bedingung erfüllt (bis $N^*=20{,}3$ auch, das liegt aber außerhalb des Fensters) [korr.], bei $N^*=2{,}38$ ist $R_\xi=4{,}5\cdot10^{-9}$. **Die Galaxie folgt dem reinen Operator, $g/g_\text{RAR}\to1$, ohne Überschuss über die Relation.** [korr.]

**Was mit der Gedächtnismode geschieht** (Zahlen aus 1.5 (e)): Sie bleibt am Ort und geht in die Strukturenergie des Operators über. Bis $R_Q=50{,}9$ kpc (Hernquist-Profil; als Punktmasse 49 kpc) ist die Bilanz exakt ($3{,}22\cdot10^{11}\,M_\odot$). [korr.] Umverteilt werden im Mittel 4,8 %, vor allem innerhalb von $2r^*$ ($r^*=8{,}35$ kpc). Sie verschwindet nicht und wird auch nicht nach außen abgegeben.

**Nachlaufzeit als Übergangszeit.** Zwei Zeiten sind zu unterscheiden:
- Das Umschalten von $w$ dauert 7 Myr.
- Die Relaxation in die Operatorform dauert $1/\lambda_m=N^*/\omega_\text{vir}$.

Die Bullet-Schranke $1/\lambda_m\ge0{,}26$ Gyr bei $\omega_\text{vir}=0{,}74$–$1{,}08$/Gyr ergibt $N^*\ge0{,}28$. Die Enge der RAR auch für spät kollabierte Galaxien verlangt $N^*\le2{,}38$ ($z_c=0{,}5$), $\le4{,}66$ ($z_c=1$) und $\le12{,}5$ ($z_c=2$). Im Fenster $0{,}28\le N^*\le2{,}38$, also höchstens 0,38 Umläufe der Zone, dauert die Relaxation in der Galaxie mit $z_c=3{,}63$ 0,08–0,65 Gyr und im Haufen ($z_c=0{,}5$) 0,38–3,2 Gyr. Die Nachlaufzeit passt also als Übergangszeit.

---

## 4. Galaxienhaufen und Bullet Cluster

### 4.1 Haufen von 10¹⁵ M☉ (z = 0,3)

**Baryonen.** Gas als β-Modell ($\beta=2/3$, $r_c=150$ kpc, $10^{14}\,M_\odot$ in 1 Mpc, bis 3 Mpc), BCG $1{,}5\cdot10^{12}$ ($a=20$ kpc), Galaxien $10^{13}$ ($a=300$ kpc). Die Linsenmasse vertritt NFW mit $M_{200}=10^{15}$, $c=4$, also $r_{200}=1903$ kpc.

**Was die Operatorform fasst.** [korr.] Am Zonenrand ist $g_N(r_{200})/a_0=0{,}071>y_Q=0{,}029$. Die Zone endet, bevor $q_\text{lok}=q_\text{kos}$ erreicht ist. Die Operatorform fasst daher nur $(\nu-1)M_b$. Der Rest $R_\xi=1-(\nu-1)M_b/M_P$ bleibt freie Gedächtnismode, und die scheinbare Gesamtmasse ist $M_b+M_P$. Die Gedächtnismode der Zone wird wie folgt angesetzt: (ii) kosmischer Anteil $q\,M_\text{ref}$, verteilt wie ein NFW-Profil (stoßfrei).

| r [kpc] | $g_N/a_0$ | ν | fehlender Faktor $M_\text{ref}/(\nu M_b)$ | $R_\xi$ | $(M_b+M_P)/M_\text{ref}$ | reiner Operator $\nu M_b/M_\text{ref}$ |
|---|---|---|---|---|---|---|
| 100 | 0,368 | 2,20 | 3,03 | 0,79 | 0,99 | 0,33 |
| 250 | 0,288 | 2,41 | 2,58 | 0,73 | 1,00 | 0,39 |
| 500 | 0,207 | 2,74 | 2,09 | 0,64 | 1,02 | 0,48 |
| 1000 | 0,125 | 3,36 | 1,56 | 0,46 | 1,03 | 0,64 |
| 1500 | 0,088 | 3,89 | 1,23 | 0,29 | 1,05 | 0,81 |

**Das Rest-Gedächtnis mit $R_\xi$ zwischen 0 und 1 liefert den fehlenden Faktor 1,2–3,0.** Das Modell trifft die Linsenmasse auf 1–5 %. Setzt man die Gedächtnismode statt über (ii) direkt als $Q$ mal die Baryonen dieses Modells an, liegt das Ergebnis um 1,34–1,40 höher. Das ist genau das Verhältnis des Baryonenanteils dieses Gasmodells (0,222) zum kosmischen (0,157). Nötig wäre bei Zeitabhängigkeit allein $R_\xi=0{,}32$–$0{,}66$; die Grenze $(\nu-1)M_b$ der Operatorform liefert es ohne zusätzlichen Parameter. [korr.]

**Grenzmasse.** Bei kosmischem Baryonenanteil fasst der Operator die ganze Gedächtnismode, wenn $g_N(r_{200})\le y_Qa_0$ ist. Daraus folgt:

| z | 0 | 0,3 | 1 |
|---|---|---|---|
| $M_{200}$-Grenze [M☉] | 3,7·10¹⁴ | 1,9·10¹⁴ | 3,6·10¹³ |

Darunter (Galaxien, Gruppen) wirkt der reine Operator, darüber bleibt Rest-Gedächtnis. Haufen sind damit **Übergangssysteme**, und das nicht wegen ihres Alters, sondern weil die Operatorform ihrer Zone, $(\nu-1)M_b$, nicht die ganze Mode fasst. [korr.]

### 4.2 Bullet Cluster

- **Das Rest-Gedächtnis ist stoßfrei.** Sein freier Anteil fällt geodätisch und geht mit den Galaxien durch den Stoß.
- **Abschätzung.** Hauptkomponente $M_{200}=1{,}5\cdot10^{15}\,M_\odot$, $c=4$ [G]. Scheinbare Masse der Gedächtnismode in 100 kpc um den BCG: $5{,}1\cdot10^{13}\,M_\odot$. Davon sind $4{,}0\cdot10^{13}$ stoßfrei ($R_\xi=0{,}79$) und liegen bei den Galaxien; die Linsenmasse dort beträgt ~$2{,}8\cdot10^{13}$ (Clowe et al. 2006 [G]). Die Größenordnung stimmt, die genaue Bilanz verlangt die gemessenen Karten.
- **Operatorgebundener Anteil.** $1{,}1\cdot10^{13}\,M_\odot$ (21 %) folgt den Baryonen, im Stoß also dem Gas, mit dem Nachlauf $1/\lambda_m\ge0{,}26$ Gyr. In den 150 Myr seit dem Stoß relaxieren höchstens $1-e^{-0{,}15/0{,}26}=44\,\%$ davon neu. [korr.]

Damit ist die Trennklasse „Linsenmasse gegen Gas" eine Zahl: der stoßfreie Anteil der scheinbaren Masse, $R_\xi(r)$.

---

## 5. Trennsignale gegen das Standardmodell

| # | Signal | Schließungsregel | Standardmodell |
|---|---|---|---|
| T1 | **Kleinskaliges Wachstum nach Kollaps** | Linear gleich bis $k_1(z)$: 0,072 (z = 0), 0,143 (z = 1), 0,267 h/Mpc (z = 2) (Lagrange). Darüber relaxiert die Gedächtnismode in die Operatorform: Innerhalb von $R_Q$ geht die gelesene Masse in $M_X(<r)=(\nu-1)M_b(<r)$ über, statt in ein NFW-Profil. [korr.] | NFW auf allen Skalen |
| T2 | **RAR ohne eigene Streuung** [korr.] | Für Galaxien mit $z_c\ge2$ ist heute $\lvert g/g_\text{RAR}-1\rvert\le3\cdot10^{-5}$ ($N^*=2{,}38$). Eine eigene Streuung aus einem Zusatzanteil gibt es nicht; es bleibt nur die Streuung der Baryonenmassen. [korr.] | Streuung aus Halo-Konzentration und Baryonenanteil |
| T3 | **Haufen als Übergangssysteme** | Stoßfreier Anteil $R_\xi=0{,}79$ (100 kpc) bis 0,29 (1,5 Mpc). Grenzmasse $3{,}7\cdot10^{14}$ (z = 0), $1{,}9\cdot10^{14}$ (z = 0,3), $3{,}6\cdot10^{13}\,M_\odot$ (z = 1). Unterhalb gilt der reine Operator. | stoßfreier Anteil 1 überall |
| T4 | **Versatz Linse–Gas in Verschmelzungen** | Der Anteil $1-R_\xi$ (21 % bei 100 kpc) folgt dem Gas mit Nachlauf. Das Linsenmaximum liegt nahe den Galaxien, aber nicht vollständig. | Linse ganz bei den Galaxien |
| T5 | **Kollapsgeschichte** | Galaxien, die kurz nach ihrem Kollaps beobachtet werden, tragen Rest-Gedächtnis. Beispiel $z_c=3{,}63$, gesehen bei $z=2$: $R_\xi=0{,}044$, $g/g_\text{RAR}-1=+1{,}3\,\%$ (10 kpc) bei $N^*=2{,}38$; bei $N^*=20{,}3$, außerhalb des Fensters $N^*\approx0{,}3$–2,4, wären es $+21\,\%$. [korr.] Spät kollabierte Galaxien ($z_c=1$, gesehen bei $z=0{,}5$): $-2\,\%$ bis $-6\,\%$ (10–50 kpc, $N^*=2{,}38$). Das Vorzeichen wechselt mit dem Radius: innen mehr, außen weniger als die RAR. | keine Kopplung der RAR an $z_c$ |
| T6 | **Nachlauf-Parameter** | Ein einziges $N^*$ im Fenster $0{,}28\le N^*\le2{,}4$ legt Bullet-Nachlauf, Galaxien-Relaxation und T5 zugleich fest. | – |

Messbar sind T2 und T5 mit Rotationskurven früher Galaxien (JWST, ALMA) und mit Linsen-RAR in Stapeln nach Umgebung bzw. Entstehungszeit (Euclid). T3 und T4 lassen sich an Haufenstichproben mit Gasprofilen und an Verschmelzungen mit geschätzter Stoßzeit prüfen.

---

## 6. Ergebnis und offene Punkte

**Ergebnis.** A13 trägt als eine Regel:
- **CMB:** unverändert gegenüber $w=1$ (0,000 %).
- **Wachstum:** linear wie im Standardmodell bis $k_1(z)$, $f\sigma_8$ auf $\le5\cdot10^{-4}$.
- **Hintergrund:** exakt staubartig.
- **Galaxien:** reiner Operator, $g/g_\text{RAR}\to1$. Die Energie der Gedächtnismode bleibt als Strukturenergie erhalten, exakt bis $R_Q$.
- **Haufen:** Weil die Operatorform nur $(\nu-1)M_b$ fasst, bleibt ein stoßfreier Rest der Mode stehen, und dieser schließt den Faktor ~2. [korr.]

Neu gegenüber der reinen Projektion ist nur der Parameter $N^*$ der Relaxation, dessen Fenster eng ist. In der Lesart des Autors heißt das: **Bindung ist die Einschränkung der Gedächtnismode auf die Operatorform.** Solange ein Bereich atmet, fällt die Gedächtnismode frei. Hört er auf zu atmen, verliert sie diesen Freiheitsgrad und bleibt als Struktur des Operators bestehen, so weit die Zone sie fassen kann.

**Offen.**
1. **Relaxation aus einem Potential.** Der dissipative Teil ist auf Fluidebene angesetzt. Offen sind die GENERIC-Klammer mit Entropiefunktional und die Frage, ob $N^*$ aus dem Modell folgt.
2. **Grenze der Operatorform kovariant.** Die Bedingung, dass die Operatorform höchstens $(\nu-1)M_b$ fasst (gelesen $\rho_B\le\rho_X$), und die Zone [korr.] (Mittelung von Ω) brauchen die Zonenzerlegung der Einzeltrichter-Hierarchie in kovarianter Form.
3. **Operatorform für lineare Störungen.** Der relaxierte Anteil im Hintergrund ($\varepsilon\sim10^{-3}/N^*$) ist hier als nicht klumpend behandelt. Wie wirkt die Operatorform auf lineare Störungen, bei denen $g_N\to0$ und ν groß ist?
4. **Schalenweise Relaxation.** Innere Schalen kollabieren früher. Rechnet man Schale für Schale statt mit einer Zone, sinkt T5 bei kleinen Radien. Dazu gehört ein N-Körper-Test des Übergangs.
5. **Boltzmann-Code.** Die Regel ist in CLASS/hi_class mit $w(k,a)$ und der Mitnahme einzubauen, und $k_1$, $k_h$ sind im nichtlinearen Spektrum zu bestimmen.
6. **Haufen mit echten Profilen.** Zu prüfen ist $R_\xi(r)$ an CLASH/LoCuSS-Haufen [G] und das Bullet-Modell mit gemessenen Karten.
7. **Herkunft von q = 0,843.** Die Regel benutzt das Gewicht; welche Größe des Modells es festlegt, bleibt offen.

---

## Literatur

- Ma C.-P., Bertschinger E. 1995, ApJ 455, 7.
- Hu W., Sugiyama N. 1995, ApJ 444, 489.
- Planck Collaboration 2020, A&A 641, A6 (Parameter, $\sigma_8=0{,}811$ [G]).
- Zwanzig R. 1961, Phys. Rev. 124, 983; Mori H. 1965, Prog. Theor. Phys. 33, 423.
- Skordis C., Złośnik T. 2021, Phys. Rev. Lett. 127, 161302 (AeST).
- Milgrom M. 2010, MNRAS 403, 886 (QUMOND).

Aus dem Gedächtnis [G]:
- Chamseddine A. H., Mukhanov V. 2013, JHEP 11, 135 (mimetische Materie).
- Lim E. A., Sawicki I., Vikman A. 2010, JCAP 05, 012 (Staub aus Nebenbedingung).
- Öttinger H. C., Grmela M. 1997, Phys. Rev. E 56, 6620 und 6633 (GENERIC).
- Gunn J. E., Gott J. R. 1972, ApJ 176, 1 (sphärischer Kollaps).
- Press W. H., Schechter P. 1974, ApJ 187, 425.
- Bardeen J. M., Bond J. R., Kaiser N., Szalay A. S. 1986, ApJ 304, 15 (BBKS); Sugiyama N. 1995, ApJS 100, 281.
- Clowe D. et al. 2006, ApJ 648, L109 (Bullet Cluster).
