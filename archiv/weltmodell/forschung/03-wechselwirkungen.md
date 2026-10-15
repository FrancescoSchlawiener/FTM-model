# Die Wechselwirkungen im Zustandsraum

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert (nur Bezeichnungen und [G], überall in der Notiz): χ als Temperaturleitfähigkeit → D_T; χ als Suszeptibilität → κ_m; Λ(p) als transversaler Jacobi-Block → J_⊥(p); Gyrator-Konstante Γ = Bℓ_w → G_y = Bℓ_g; η als Viskosität → μ_v; ν₀ als Arrhenius-Vorfaktor → f₀; ξ(t) als Rauschen → n(t); ξ als Schwingungsamplitude der Abstrahlung → û; w als Strömungsgeschwindigkeit → u_s; σ als Leitfähigkeit → σ_el (σ₀ → σ_el,0); σ als mechanische Spannung → T_s (σ₀ → T_s,0, σ_v → T_s,v, σ_a, σ_f′ → T_s,a, T′_s,f); σ als Entropiekoordinate C_th ln T → S_th; s als Entropiedichte → s_th; Q̂ als Ladungsamplitude → Q̂_el; τ_D als Zeitskala der Schädigung → τ_Sch; Durchbiegung δ → δ_b; Masse und Radius von Sirius B [G]. L (Induktivität) und q (Schleifenverstärkung) bleiben, wie im Modell.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/wechselwirkungen.py`. Jede Zahl in dieser Notiz stammt aus deren Ausgabe.

**Worum es geht.** Die Tafeln 1–4 des Autors (*Übersicht & Wirkfelder*, *Generatoren von Flussrichtungen*, *geschlossene Wirklogik*, *Topologie*) benennen die Wirkfelder, ihre Rollen (tangential, transversal, parametrisch) und ihre Kopplungsklassen (A Energie, B Phase, C Struktur, D Dissipation). Diese Notiz schreibt jede Kopplung als Term eines einzigen Vektorfelds $\dot z = F(z)$ und zeigt, dass sich die drei Rollen des Autors **exakt** aus der Thermodynamik ergeben:

- *tangential* ist der reversible Teil: Er erhält Energie und Entropie.
- *transversal* ist der dissipative Teil zusammen mit den offenen Toren: Er erhält die Energie, erzeugt aber Entropie, oder er tauscht mit der Umgebung aus.
- *parametrisch* ist der langsame Teil der Ordnung $\varepsilon$.

Danach werden die Rückkopplungsschleifen der Tafeln als Zyklen im Kopplungsgraphen berechnet, die Gravitation wird als langsamer Parameter bis zum Kollaps $T^2\to T^1$ verfolgt, und ein Minimalsystem wird durchgerechnet.

---

## 1  Der Zustandsvektor

### 1.1 Feldform

$$
z \;=\; \big(\,\underbrace{\mathbf E,\mathbf B}_{z_{\rm EM}},\;\underbrace{\mathbf u,\dot{\mathbf u}}_{z_{\rm mech}},\;\underbrace{T}_{z_{\rm th}},\;\underbrace{\mathbf C,\sigma_{\rm el},\varepsilon,\mu,\alpha,S_b;\;D,A,\boldsymbol\varepsilon_p,\mathbf M_{\rm irr}}_{z_{\rm mat}}\;;\;\underbrace{\rho,\mathcal G,\mathbf T_{s,0},\Phi}_{p}\big)
$$

| Block | Variablen | Bedeutung |
|---|---|---|
| $z_{\rm EM}$ | $\mathbf E(\mathbf x),\mathbf B(\mathbf x)$, abgeleitet $\mathbf D,\mathbf H,\rho_e,\mathbf J$ | Feldzustand |
| $z_{\rm mech}$ | Auslenkung $\mathbf u(\mathbf x)$, Geschwindigkeit $\dot{\mathbf u}$, Dehnung $\boldsymbol\varepsilon=\operatorname{sym}\nabla\mathbf u$ | Konfiguration und Bewegung |
| $z_{\rm th}$ | Temperatur $T(\mathbf x)$ bzw. Entropiedichte $s_{\rm th}$ [korr.] | thermischer Zustand |
| $z_{\rm mat}$ | Steifigkeit $\mathbf C$, Leitfähigkeit $\sigma_{\rm el}$, Permittivität $\varepsilon$, Permeabilität $\mu$, Ausdehnungskoeffizient $\alpha$, Seebeck-Koeffizient $S_b$ | Stoffparameter |
| | innere Variablen: Schädigung $D\in[0,1]$, Alterungsgrad $A\in[0,1]$, plastische Dehnung $\boldsymbol\varepsilon_p$, irreversible Magnetisierung $\mathbf M_{\rm irr}$ | Vorgeschichte als Koordinate |
| $p$ | Dichte $\rho$, Geometrie $\mathcal G$, Vorspannung $\mathbf T_{s,0}$ [korr.], Gravitationspotential $\Phi$ mit $\mathbf g=-\nabla\Phi$, $d\tau/dt=\sqrt{1+2\Phi/c^2}$ | langsame Parameter |

Die inneren Variablen $D, A, \boldsymbol\varepsilon_p, \mathbf M_{\rm irr}$ sind die Stelle, an der **Gedächtnis** in den Zustandsraum eingeht. Sie sind keine Speicherobjekte, sondern zusätzliche Koordinaten des aktuellen Punkts, und ihre Dynamik ist einseitig ($\dot D\ge0$, $\dot A\ge0$). Das entspricht genau der Gedächtniskoordinate $y$ des Memory-Kernels im Zustandsraum.

### 1.2 Modalform und Torus-Koordinaten

Projiziert man auf je eine Eigenmode (Galerkin), dann wird jeder schnelle Block zu einem Oszillator mit Amplitude und Phase:

$$
Q_{\rm el} = \hat Q_{\rm el}\cos\theta_{\rm EM},\quad I = -\hat Q_{\rm el}\,\omega_{\rm EM}\sin\theta_{\rm EM};\qquad
x = X\cos\theta_{\rm mech},\quad v = -X\,\omega_m\sin\theta_{\rm mech}.
$$

$Q_{\rm el}$ ist die modale Ladung, $\hat Q_{\rm el}$ ihre Amplitude [korr.], $I$ der Strom, $x$ die modale Auslenkung und $v$ die Geschwindigkeit. Die Aktionen sind $J_{\rm EM}=E_{\rm EM}/\omega_{\rm EM}$ und $J_{\rm mech}=E_{\rm mech}/\omega_m$. In der Torus-Normalform heißen die Koordinaten dann

$$
\theta=(\theta_{\rm EM},\theta_{\rm mech})\in T^2,\qquad
\rho_\perp=(\hat Q_{\rm el}-\hat Q_{\rm el}^*,\;X-X^*,\;T-T^*),\qquad
p=(D,A,g,\dots).
$$

---

## 2  Die Kopplungsmatrix

### 2.1 Kurzform

Rolle: **t** = tangential, **⊥** = transversal, **p** = parametrisch. Klasse: A Energie, B Phase, C Struktur, D Dissipation. $\tau_{\rm EM}=2\pi/\omega_{\rm EM}$, $\tau_m=2\pi/\omega_m$, $\tau_{\rm th}=C_{\rm th}/G$ (konzentriert) bzw. $\ell^2/D_T$ (Diffusion, Temperaturleitfähigkeit $D_T$) [korr.], $\tau_\sigma=\varepsilon/\sigma_{\rm el}$ (dielektrische Relaxation).

| Nr | Paar | Gleichung (Kern) | Energiefluss | Rolle | Klasse | Zeitskala |
|---|---|---|---|---|---|---|
| K1 | EM → Mech | $\mathbf f=\rho_e\mathbf E+\mathbf J\times\mathbf B=\nabla\!\cdot\mathbf T_M-\partial_t\mathbf g_{\rm EM}$ | Feld → Bewegung, Leistung $(\mathbf J\times\mathbf B)\cdot\dot{\mathbf u}$ | t (Ströme), ⊥ (Maxwell-Druck) | A (+B) | $\tau_{\rm EM}$ |
| K2 | Mech → EM | $\mathbf E'=\mathbf E+\dot{\mathbf u}\times\mathbf B$ (Bewegungs-EMK), $\omega_{\rm EM}(\mathbf u)=1/\sqrt{L(\mathbf u)C(\mathbf u)}$ | Bewegung → Feld | t (Induktion), p (Geometrie) | A, B/C | $\tau_m$ |
| K3 | EM → Wärme | $\mathbf J'\cdot\mathbf E'=\sigma_{\rm el}\lvert\mathbf E'\rvert^2$; Wirbelstrom, Hysterese $f\oint H\,dB$ | Feld → Wärme, irreversibel | ⊥ | A/D | $\tau_\sigma$ … $\tau_{\rm EM}$ |
| K4 | Mech → Wärme | $\mathbf T_{s,v}:\dot{\boldsymbol\varepsilon}=\boldsymbol\mu_v\,\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}$ [korr.], Reibung $F_R v$, Zener-Dämpfung | Bewegung → Wärme | ⊥ | A/D | $\tau_m$ |
| K5 | Wärme → Mech/Geometrie | $\boldsymbol\varepsilon_{\rm th}=\alpha\,\Delta T\,\mathbf 1$, $\mathbf T_s=\mathbf C:(\boldsymbol\varepsilon-\alpha\Delta T\mathbf 1)$ [korr.] → $k(T)$ → $\omega_m(T)$ | reversibler Kelvin-Anteil $\pm T\alpha\,\mathbf C:\dot{\boldsymbol\varepsilon}$ | p (+t) | C | $\tau_{\rm th}$ |
| K6 | Wärme → EM | $\sigma_{\rm el}(T),\mu(T),\varepsilon(T)$; $\mathbf J=\sigma_{\rm el}(\mathbf E'-S_b\nabla T)$, $\mathbf j_{\rm th}=\Pi\mathbf J-k\nabla T$, $\Pi=S_bT$ | Wärme ↔ Strom (Onsager) | p; ⊥ (Tor) | C; A | $\tau_{\rm th}$ |
| K7 | EM ↔ Mech über Material | Piezo $\mathbf S=\mathbf s^E\mathbf T+\mathbf d^{\!\top}\mathbf E$, $\mathbf D=\mathbf d\,\mathbf T+\boldsymbol\varepsilon^T\mathbf E$; Magnetostriktion $\varepsilon=\tfrac32\lambda_s(M/M_s)^2$ | reversibel, beide Richtungen | t | B (+A) | $\tau_m$ |
| K8 | Mech → Akustik → Abstrahlung | $P_{\rm ac}=\tfrac12R_{\rm rad}\,\omega^2\lvert\hat u\rvert^2$ [korr.], $R_{\rm rad}=\rho_0c_0S\,\sigma_{\rm rad}(ka)$ | Bewegung → Schall → außen | ⊥ (+t: Mitschwingmasse) | D (+B) | $2m/R_{\rm rad}$ |
| K9 | EM → Abstrahlung | $P=\dfrac{Q_{\rm el}^2a^2}{6\pi\varepsilon_0c^3}$ (Larmor), $R_{\rm rad}=80\pi^2(\ell/\lambda)^2\,\Omega$ | Feld → außen | ⊥ (Grenzfall zu t) | D | $2L/R_{\rm rad}$ |
| K10 | Strömung ↔ Mech/Wärme | $\rho(\partial_t\mathbf u_s+\mathbf u_s\!\cdot\!\nabla\mathbf u_s)=-\nabla p+\mu_v\nabla^2\mathbf u_s$ [korr.], $F_y=\tfrac12\rho U^2 d\,[a_1\dot y/U-a_3(\dot y/U)^3]$, $h(U)$ | Strömung → Bewegung (Selbsterregung), Wärme → Strömung (Kühlung) | t + ⊥ | A + B + D | $d/U$ |
| K11 | Chemie/Alterung | $\dot A=f_0\,e^{-E_a/k_BT}(1-A)\ge0$ [korr.]; $\sigma_{\rm el}(A),\,k(A)$ | irreversibel, Energie → Zustand | p (irreversibel) | C/D | $e^{E_a/k_BT}/f_0$ [korr.] |
| K12 | Ermüdung | $\dot D=f_D(Y,D)\ge0$, $Y=-\partial\Psi/\partial D$; $k=k_0(1-D)$ | Belastung → Zustand, Wärme $Y\dot D$ | p (+⊥) | C + D | $1/\dot D$ |
| K13 | Rauschen/Fluktuation | $\langle n(t)n(t')\rangle=2\gamma k_BT\,\delta(t-t')$, $S_V=4k_BTR$, $S_x(\omega)=\tfrac{2k_BT}{\omega}\operatorname{Im}\kappa_m(\omega)$ [korr.] | Wärme ↔ Phase/Amplitude, im Mittel null | ⊥ (Band), t (Phasendiffusion) | B + D | alle |
| K14 | Gravitation | $k_{\rm eff}=k_{\rm el}\pm mg/\ell$, $\omega_i^{\rm obs}=\omega_i\,d\tau/dt$ | parametrische Arbeit $\tfrac12\dot k\,x^2$ | p | C | $g/\dot g$ |
| K15 | Zustand → Energie | $R(A),\;k(D),\;\mu(\mathbf M_{\rm irr})$ (Preisach) | ändert die Antwort, nicht den Inhalt | p | C | $\tau_{\rm Sch},\tau_A$ [korr.] |
| K16 | Zustand → Struktur | Kriechen $\dot\varepsilon_{\rm cr}=A_NT_s^n e^{-E_{\rm cr}/RT}$ [korr.], Setzen, $\boldsymbol\varepsilon_p$ | Spannung → bleibende Form | p | C | $\tau_{\rm cr}$ |

Die sechs Pfeile der Tafel *Übersicht & Wirkfelder* verteilen sich so auf die Zeilen: Energie → Struktur (K1 Maxwell-Druck, K5, K7), Struktur → Energie (K2 Geometrie, K6 $\sigma_{\rm el}(T)$ [korr.], K14), Energie → Zustand (K3, K11, K12), Zustand → Energie (K15, K13), Struktur → Zustand (K12 über Spannungskonzentration), Zustand → Struktur (K16).

### 2.2 Die Paare im Einzelnen

**K1 · EM → Mechanik (Lorentzkraft, Maxwell-Spannung).**
$$
\mathbf f=\rho_e\mathbf E+\mathbf J\times\mathbf B=\nabla\cdot\mathbf T_M-\varepsilon_0\mu_0\,\partial_t\mathbf S,\qquad
T_{M,ij}=\varepsilon_0\big(E_iE_j-\tfrac12\delta_{ij}E^2\big)+\tfrac1{\mu_0}\big(B_iB_j-\tfrac12\delta_{ij}B^2\big).
$$
Die Leistungsdichte zerfällt exakt mit $\mathbf E'=\mathbf E+\dot{\mathbf u}\times\mathbf B$:
$$
\mathbf J\cdot\mathbf E=\underbrace{(\mathbf J\times\mathbf B)\cdot\dot{\mathbf u}}_{\text{Arbeit an der Materie}}+\underbrace{\mathbf J\cdot\mathbf E'}_{\text{Umsatz im Ruhesystem}} .
$$
Der erste Term *verschiebt* Energie ins Mechanische, und er ist tangential. Der Maxwell-Druck $p_M=B^2/2\mu_0$ an Grenzflächen wirkt dagegen als statische Querkraft bzw. Vorspannung, also transversal und in Tafel 3 („Maxwell-Spannungskräfte – transversal“) genau so eingeordnet. Konzentriert man das auf eine Mode (Tauchspule, Wandler), dann gilt $F=G_yI$ mit $G_y=B\ell_g$ (Leiterlänge $\ell_g$). [korr.]

**K2 · Mechanik → EM (Induktion, Geometrie).** Zu K1 gehört das Gegenstück $U_{\rm ind}=G_yv$. Es ist dieselbe Kopplungskonstante mit umgekehrtem Vorzeichen in der Energiebilanz: $-G_yvI$ im Kreis, $+G_yIv$ in der Mechanik. [korr.] Das ist ein **Gyrator**, antisymmetrisch und damit ein Onsager-Casimir-Paar, das ungerade im Magnetfeld ist. Über die Geometrie wirkt die Mechanik außerdem parametrisch: $\omega_{\rm EM}(x)=1/\sqrt{L(x)C(x)}$. Die Kräfte $F=\tfrac12I^2\partial_xL+\tfrac12V^2\partial_xC$ schließen dabei die Schleife S5 (Pull-in, §4).

**K3 · EM → Wärme.** Ohmsch gilt $\mathbf J'\cdot\mathbf E'=\sigma_{\rm el}|\mathbf E'|^2\ge0$. Für Wirbelströme in Blechen der Dicke $d$ gilt $p_{\rm eddy}=\pi^2d^2f^2B_p^2/(6\rho_{\rm el})$ (klassischer Steinmetz-Anteil). Die Hysterese liefert $p_h=f\oint H\,dB$ pro Volumen. Energie fließt vom Feld in den Wärmeknoten, und dabei entsteht Entropie $\sigma_{\rm el}|\mathbf E'|^2/T$. Bei Tafel 3 („EM → Wärme ist keine Ergänzung, sie ist strukturell eingebaut“) ist das wörtlich zu nehmen: Derselbe Term $\mathbf J\cdot\mathbf E$ steht mit entgegengesetztem Vorzeichen im Poynting-Satz und in der Wärmegleichung.

**K4 · Mechanik → Wärme.** Viskose innere Dämpfung $\boldsymbol\mu_v\,\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}$ [korr.], Coulomb-Reibung $\mu_RN|v|$ und thermoelastische Zener-Dämpfung
$$
Q_Z^{-1}=\Delta_E\,\frac{\omega\tau_Z}{1+\omega^2\tau_Z^2},\qquad \Delta_E=\frac{E\alpha^2T}{\rho c_p}.
$$
Die Zener-Dämpfung ist die dissipative Spur des *reversiblen* Kelvin-Terms aus K5 und entsteht erst, wenn innerhalb der Struktur Temperaturgradienten mit der Relaxationszeit $\tau_Z$ auftreten.

**K5 · Wärme → Mechanik/Geometrie.** Es gilt $\mathbf T_s=\mathbf C:(\boldsymbol\varepsilon-\alpha\Delta T\mathbf 1)$. In einer eingespannten Struktur wird die Ausdehnung zur Vorspannung $T_{s,0}=-E\alpha\Delta T$ und damit zur Steifigkeitsänderung. [korr.] Konzentriert ergibt das $k(T)=k_0-\kappa_T(T-T_0)$ und
$$
\frac{\partial\omega_m}{\partial T}=-\frac{\kappa_T}{2m\,\omega_m}.
$$
Thermodynamisch konsistent ist das nur mit dem **Kelvin-Term** in der Wärmegleichung, $C_{\rm th}\dot T\ni T\,\partial_Tk\;x\dot x=-\kappa_TT\,x v$. Er ist reversibel, also tangential, und im Minimalsystem enthalten.

**K6 · Wärme → EM.** Parametrisch wirken $\sigma_{\rm el}(T)=\sigma_{\rm el,0}/(1+\alpha_R\Delta T)$ [korr.], $\mu(T)$ (Curie) und $\varepsilon(T)$. Als Transport kommen Seebeck und Peltier hinzu:
$$
\mathbf J=\sigma_{\rm el}(\mathbf E'-S_b\nabla T),\qquad \mathbf j_{\rm th}=\Pi\,\mathbf J-k\nabla T,\qquad \Pi=S_bT\ \ (\text{zweite Kelvin-Relation}).
$$
Diese zweite Kelvin-Relation ist die Onsager-Symmetrie $L_{12}=L_{21}$ (§3.5).

**K7 · EM ↔ Mechanik über das Material.** Beim Piezoeffekt folgen direkter und inverser Effekt aus einer einzigen freien Enthalpie $G=-\tfrac12\mathbf s\mathbf T^2-\mathbf d\mathbf T\mathbf E-\tfrac12\varepsilon\mathbf E^2$. Deshalb ist $\mathbf d$ in beide Richtungen dasselbe (Maxwell-Reziprozität). Die Kopplung ist reversibel und verschiebt die Frequenz:
$$
k^D=\frac{k^E}{1-k_c^2},\qquad k_c^2=\frac{d^2}{s^E\varepsilon^T}.
$$
Das ist Klasse B im Sinn des Autors: *wer moduliert wessen Zeit*. Bei der Magnetostriktion gilt $\varepsilon=\tfrac32\lambda_s(M/M_s)^2$, der inverse Effekt ist der Villari-Effekt, und der $\Delta E$-Effekt verschiebt ebenfalls $\omega_m$.

**K8 · Mechanik → Akustik → Abstrahlung.** Die Strahlungsimpedanz $Z_{\rm rad}=R_{\rm rad}+iX_{\rm rad}$ hat einen Realteil, der Energie nach außen trägt ($P_{\rm ac}=\tfrac12R_{\rm rad}\omega^2|\hat u|^2$ mit der Schwingungsamplitude $\hat u$ [korr.], transversal) und einen Imaginärteil, die mitschwingende Masse, die $\omega_m$ verschiebt (Klasse B). Für eine atmende Kugel gilt $R_{\rm rad}=\rho_0c_0\,4\pi a^2\,(ka)^2/(1+(ka)^2)$.

**K9 · EM-Abstrahlung.** Die Larmor-Leistung bzw. der Strahlungswiderstand $R_{\rm rad}$ bilden den äußeren Abfluss. Auch hier hat die Rückwirkung (Abraham-Lorentz) einen Blindanteil. Deshalb nennt der Autor die Abstrahlung den „Grenzfall zwischen transversal und tangential“.

**K10 · Strömung als Mitkoppler.** Strömung wirkt dreifach. Erstens liefert sie eigene Zyklen (Wirbelablösung mit der Strouhal-Frequenz $f=St\,U/d$) und damit tangentiale Richtungen. Zweitens kann sie negative Dämpfung erzeugen (Galloping, Flattern): $F_y=\tfrac12\rho U^2d\,[a_1\dot y/U-a_3(\dot y/U)^3]$. Das ist genau die van-der-Pol-Form $(a_mv-b_mv^3)$ des Minimalsystems, Energie fließt also aus der Strömung in die Schwingung. Drittens transportiert sie Wärme über $h(U)\propto Nu(U)\,k/d$ und ist damit transversal. „Mitkoppler“ heißt: Strömung schließt Schleifen zwischen Mechanik, Wärme und Struktur, ohne selbst Auslöser zu sein.

**K11 · Chemie/Alterung.**
$$
\dot A=f_0\,e^{-E_a/k_BT}\,(1-A)\ \ge 0 .
$$
Der Vorfaktor heißt $f_0$ [korr.]. Der Alterungsgrad $A$ verändert die Materialparameter, etwa $R=R_0(1+\alpha_R\Delta T)(1+\beta_AA)$. Die Temperaturempfindlichkeit beträgt $\partial\ln\dot A/\partial T=E_a/(k_BT^2)$. Für $E_a=0{,}7$ eV beschleunigt eine Erwärmung von 300 K auf 310 K die Alterung um den Faktor **2,395**. Die Dynamik ist irreversibel ($\dot A\ge0$) und parametrisch ($\varepsilon$-langsam).

**K12 · Ermüdung als gespeicherte Belastung.** Die Schädigung $D$ ist eine innere Variable der freien Energie, $\Psi=\tfrac12k_0(1-D)x^2+\dots$. Die zu $D$ konjugierte Kraft ist die Energiefreisetzungsrate
$$
Y=-\frac{\partial\Psi}{\partial D}=\tfrac12k_0x^2\ \ge0,\qquad \dot D=f_D(Y,D)\ \ge0,\qquad \text{Dissipation } Y\dot D\ \ge0 .
$$
Beispiele für $f_D$ sind Kachanov $\dot D=C\,(T_s/(1-D))^m$ und Miner/Basquin $N_f=(T_{s,a}/T_{s,f}')^{1/b}$ [korr.]. Die freigesetzte elastische Energie $Y\dot D$ geht in Wärme über. Die Steifigkeit sinkt, und $\omega_m$ wandert: Das ist der Satz des Autors „Ermüdung macht Systeme pfadabhängig“.

**K13 · Rauschen (Fluktuations-Dissipations-Beziehung).** Jeder dissipative Kanal ist zugleich eine Rauschquelle:
$$
\langle n(t)n(t')\rangle=2\gamma k_BT\,\delta(t-t'),\qquad S_V=4k_BTR,\qquad S_x(\omega)=\frac{2k_BT}{\omega}\operatorname{Im}\kappa_m(\omega)
$$
(Einstein, Johnson-Nyquist, Callen-Welton; Rauschkraft $n(t)$, Suszeptibilität $\kappa_m$ [korr.]). Transversal verbreitert das Rauschen den Torus zu einem Band der Dicke $\langle\delta\rho_\perp^2\rangle\simeq D_\perp/|\lambda_\perp|$ mit $D_\perp\propto\gamma k_BT$. Tangential erzeugt es Phasendiffusion $D_\theta\propto\gamma k_BT/E^*$. Das Rauschen ist damit „eine Eigenschaft des Zustands selbst“, wie der Autor sagt, denn Stärke und Dissipation sind ein und dieselbe Größe $\gamma$.

**K14 · Gravitation als parametrische Modulation** wird in §5 ausführlich behandelt: Sie wirkt über die Vorspannung $k_{\rm eff}=k_{\rm el}\pm mg/\ell$, über die Zeitdilatation $\omega_i^{\rm obs}=\omega_i\,d\tau/dt$ und über die Geometrie (Durchbiegung → Spalt → $G_y$, $C$) [korr.].

**K15/K16 · Zustand → Energie, Zustand → Struktur.** Das sind die Rückkopplungen des Gedächtnisses. $R(A)$, $k(D)$ und $\mu(\mathbf M_{\rm irr})$ ändern die *Antwort*, also Parameter des Vektorfelds. Kriechen (Norton) und plastische Dehnung ändern die *Form*. Formal ist Gedächtnis hier $\partial F/\partial(D,A,\boldsymbol\varepsilon_p,\mathbf M_{\rm irr})\ne0$: eine Abhängigkeit des Vektorfelds von Koordinaten, die einseitig laufen. Ein gespeichertes Objekt gibt es dabei nicht.

---

## 3  Das Gesamtsystem

### 3.1 Feldform

$$
\begin{aligned}
&\partial_t\mathbf B=-\nabla\times\mathbf E,\qquad \partial_t\mathbf D=\nabla\times\mathbf H-\mathbf J,\qquad \nabla\cdot\mathbf D=\rho_e,\ \ \nabla\cdot\mathbf B=0,\\
&\mathbf D=\varepsilon(T,\boldsymbol\varepsilon,A)\mathbf E+\mathbf d:\mathbf T_s,\qquad \mathbf B=\mu(T,D,\mathbf M_{\rm irr})\mathbf H,\qquad \mathbf J=\sigma_{\rm el}(T,A)\big(\mathbf E+\dot{\mathbf u}\times\mathbf B-S_b\nabla T\big),\\
&\rho\,\ddot{\mathbf u}=\nabla\cdot\big(\mathbf T_s+\mathbf T_M\big)+\rho\,\mathbf g,\qquad
\mathbf T_s=\mathbf C(D):(\boldsymbol\varepsilon-\alpha\Delta T\mathbf 1-\boldsymbol\varepsilon_p)-\mathbf d^{\!\top}\mathbf E+\boldsymbol\mu_v:\dot{\boldsymbol\varepsilon},\\
&\rho c\,\partial_tT=\nabla\cdot(k\nabla T)+\mathbf J\cdot\mathbf E'+\boldsymbol\mu_v\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}+Y\dot D-T\alpha\,\mathbf 1:\mathbf C:\dot{\boldsymbol\varepsilon}-T\,\mathbf J\cdot\nabla S_b,\\
&\dot D=f_D(Y,D)\ge0,\qquad \dot A=f_0e^{-E_a/k_BT}(1-A)\ge0,\qquad \dot{\boldsymbol\varepsilon}_p=f_p(\mathbf T_s,T),\qquad \dot p=\varepsilon\,G_p(z,p).
\end{aligned}
$$

Spannung $\mathbf T_s$, Viskosität $\boldsymbol\mu_v$, Vorfaktor $f_0$. [korr.]

### 3.2 Energiebilanz: welche Terme verschieben und welche dissipieren

Die Energiedichten sind $u_{\rm EM}=\tfrac12(\mathbf E\cdot\mathbf D+\mathbf B\cdot\mathbf H)$, $\tfrac12\rho|\dot{\mathbf u}|^2$ (kinetisch), $W_{\rm el}=\tfrac12\boldsymbol\varepsilon:\mathbf C:\boldsymbol\varepsilon$ (elastisch) und $\rho e_{\rm th}$ (thermisch). Die vier Teilbilanzen lauten:

$$
\begin{aligned}
\text{Poynting:}\quad &\partial_tu_{\rm EM}+\nabla\cdot\mathbf S=-\mathbf J\cdot\mathbf E=-(\mathbf J\times\mathbf B)\cdot\dot{\mathbf u}-\mathbf J\cdot\mathbf E',\\
\text{Mechanik:}\quad &\partial_t\big(\tfrac12\rho\dot{\mathbf u}^2+W_{\rm el}\big)-\nabla\cdot(\mathbf T_s\dot{\mathbf u})=+(\mathbf J\times\mathbf B)\cdot\dot{\mathbf u}-\boldsymbol\mu_v\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}-Y\dot D+T\alpha\,\mathbf 1:\mathbf C:\dot{\boldsymbol\varepsilon}+\rho\mathbf g\cdot\dot{\mathbf u}+\tfrac12\boldsymbol\varepsilon:\dot{\mathbf C}_{\rm ext}:\boldsymbol\varepsilon,\\
\text{Wärme:}\quad &\partial_t(\rho e_{\rm th})+\nabla\cdot\mathbf j_{\rm th}=+\mathbf J\cdot\mathbf E'+\boldsymbol\mu_v\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}+Y\dot D-T\alpha\,\mathbf 1:\mathbf C:\dot{\boldsymbol\varepsilon}.
\end{aligned}
$$

Summiert ergibt sich
$$
\boxed{\;\partial_t\big(u_{\rm EM}+\tfrac12\rho\dot{\mathbf u}^2+W_{\rm el}+\rho e_{\rm th}\big)+\nabla\cdot\big(\mathbf S+\mathbf j_{\rm th}-\mathbf T_s\dot{\mathbf u}\big)=\rho\,\mathbf g\cdot\dot{\mathbf u}+\tfrac12\boldsymbol\varepsilon:\dot{\mathbf C}_{\rm ext}:\boldsymbol\varepsilon\;}
$$

(Spannung $\mathbf T_s$, Viskosität $\boldsymbol\mu_v$ [korr.].) Rechts stehen nur noch die äußeren bzw. parametrischen Einspeisungen, das sind Gravitation und extern getriebene Parameteränderungen. Jeder innere Kopplungsterm taucht **zweimal mit entgegengesetztem Vorzeichen** auf. Er verschiebt also nur Energie:

| Term | von → nach | verschiebt | erzeugt Entropie |
|---|---|:-:|:-:|
| $(\mathbf J\times\mathbf B)\cdot\dot{\mathbf u}$ | EM → Mechanik | ● | – |
| $T\alpha\,\mathbf 1:\mathbf C:\dot{\boldsymbol\varepsilon}$ (Kelvin) | Mechanik ↔ Wärme | ● | – |
| Seebeck/Peltier ($\Pi=S_bT$) | Wärme ↔ Strom | ● | – |
| Piezo, Magnetostriktion | EM ↔ Mechanik | ● | – |
| $\mathbf J\cdot\mathbf E'$ | EM → Wärme | ● | ● $\sigma_{\rm el}|\mathbf E'|^2/T$ [korr.] |
| $\boldsymbol\mu_v\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}$ [korr.] | Mechanik → Wärme | ● | ● |
| $Y\dot D$ | Elastik → Wärme | ● | ● |
| $\nabla\cdot\mathbf j_{\rm th}$, Abstrahlung | System → Umgebung | Fluss über den Rand | ● (in der Umgebung) |

**Verschieben** heißt: Der Term steht in zwei Teilbilanzen. **Dissipieren** heißt: Der Term steht zusätzlich in der Entropieproduktion. Auch Dissipation vernichtet keine Energie. Sie verschiebt Energie in den Wärmeknoten, und zwar *einseitig*.

### 3.3 Entropieproduktion

$$
\Pi_s=\frac{\mathbf J\cdot\mathbf E'_{\rm eff}}{T}+\frac{\boldsymbol\mu_v\dot{\boldsymbol\varepsilon}:\dot{\boldsymbol\varepsilon}}{T}+\frac{Y\dot D}{T}+\frac{\mathbf H\cdot\dot{\mathbf M}_{\rm irr}}{T}+\frac{k|\nabla T|^2}{T^2}=\sum_\alpha J_\alpha X_\alpha=\mathbf X^{\!\top}\mathbf L\,\mathbf X\ \ge\ 0 .
$$

Jeder Summand ist ein Produkt aus Fluss und konjugierter Kraft (Viskosität $\boldsymbol\mu_v$ [korr.]). Im linearen Bereich ist $J_\alpha=\sum_\beta L_{\alpha\beta}X_\beta$ mit positiv semidefinitem $\mathbf L$. Die Terme für Ermüdung und Alterung sind **einseitig** ($\dot D,\dot A\ge0$) und tragen deshalb nicht-negativ bei.

### 3.4 Die exakte Zerlegung in die drei Rollen: GENERIC

Die Rollen des Autors entsprechen genau der Struktur des Formalismus von Grmela und Öttinger (GENERIC, *General Equation for the Non-Equilibrium Reversible-Irreversible Coupling*):

$$
\boxed{\;\dot z=\underbrace{\mathcal L(z)\,\nabla E}_{F_{\rm tan}}\;+\;\underbrace{\mathcal M(z)\,\nabla S\;+\;F_{\rm port}}_{F_{\rm trans}}\;+\;\varepsilon\,F_{\rm par}(z)\;}
$$

$$
\mathcal L=-\mathcal L^{\!\top},\quad \mathcal L\nabla S=0;\qquad \mathcal M=\mathcal M^{\!\top}\succeq0,\quad \mathcal M\nabla E=0 .
$$

Aus diesen vier Bedingungen folgt sofort:

| Teil | $dE/dt$ | $dS/dt$ | Phasenraumvolumen | Rolle beim Autor |
|---|---|---|---|---|
| $F_{\rm tan}=\mathcal L\nabla E$ | $\nabla E\cdot\mathcal L\nabla E=0$ | $\nabla S\cdot\mathcal L\nabla E=-\nabla E\cdot\mathcal L\nabla S=0$ | erhalten (Liouville) | Rotation, Zyklen, Phasen |
| $\mathcal M\nabla S$ | $\nabla E\cdot\mathcal M\nabla S=0$ | $\nabla S\cdot\mathcal M\nabla S\ge0$ | kontrahiert | Kontraktion, Verlust |
| $F_{\rm port}$ | $P_{\rm in}-P_{\rm out}$ | Entropiefluss | — | Quellen, Abflüsse, Kühlung |
| $\varepsilon F_{\rm par}$ | parametrische Arbeit | $\ge0$ (einseitige Variablen) | $O(\varepsilon)$ | Bewegung des Torus selbst |

Die Onsager-Symmetrie ist $\mathcal M=\mathcal M^{\!\top}$. Die Casimir-Antisymmetrie der magnetfeldungeraden Kopplungen (Lorentz/Induktion) sitzt in $\mathcal L=-\mathcal L^{\!\top}$. **Die Tafel „Generatoren von Flussrichtungen“ ist damit die Aussage, dass jedes thermodynamisch zulässige Vektorfeld diese Form hat.** Tangential heißt „liegt in $\mathcal L$“, transversal heißt „liegt in $\mathcal M$ oder an einem Tor“, parametrisch heißt „ist langsam“.

**Folgesatz (Liouville).** Weil $F_{\rm tan}$ divergenzfrei ist (in den Koordinaten $(Q_{\rm el},I,x,v,S_{\rm th})$ mit $S_{\rm th}=C_{\rm th}\ln T$ [korr.]), gilt für die Lyapunov-Exponenten
$$
\sum_k\lambda_k=\langle\operatorname{div}F\rangle=\langle\operatorname{div}F_{\rm trans}\rangle+O(\varepsilon).
$$
**Folglich kann eine tangentiale Kopplung allein niemals verriegeln.** Beim Einrasten $T^2\to T^1$ wird eine Phasenrichtung kontrahierend, und diese Kontraktion muss aus dem transversalen Teil stammen. Verriegelung ist immer ein Zusammenspiel von Phase (B) und Dissipation bzw. Sättigung (D).

### 3.5 Onsager und Casimir, konkret

*Thermoelektrisches Element* mit den Kräften $X_e=\Delta V/T$, $X_q=\Delta T/T^2$ und den Flüssen $I$, $J_q$:
$$
\mathbf L=\begin{pmatrix}T/R & S_bT^2/R\\ S_bT^2/R & S_b^2T^3/R+KT^2\end{pmatrix},\qquad L_{12}=L_{21}\iff \Pi=S_bT .
$$
Die Rechnung (Teil 4) liefert am Arbeitspunkt $L=\begin{pmatrix}10{,}95218&0{,}66287\\0{,}66287&1{,}50537\end{pmatrix}$ mit $L_{12}-L_{21}=0$ und den Eigenwerten 1,45909 und 10,99846, beide $>0$. Über 5000 Zufallskräfte ist das Minimum von $\mathbf X^{\!\top}\mathbf L\mathbf X$ gleich $2{,}398\cdot10^{-3}\ge0$. In der Flussdarstellung heben sich die Kreuzterme exakt auf: $\Pi_s=RI^2/T+K\Delta T^2/T^2$. Im dynamischen Modell ist der Seebeck/Peltier-Beitrag zur Gesamtentropieproduktion $|\Delta\Pi|=0$.

*Reaktive Kopplungen* sind antisymmetrisch, und zwar in der Metrik, die aus $\nabla^2E$ bzw. $-\nabla^2S$ folgt:

| Kopplung | Matrixelemente | Charakter |
|---|---|---|
| Gyrator (Lorentz/Induktion) | $J[\tilde I,\tilde v]=-0{,}03000$, $J[\tilde v,\tilde I]=+0{,}03000$ | Casimir, ungerade im Magnetfeld |
| Thermoelastik (Metrik $C_{\rm th}/T$) | $J[\tilde v,\tilde T]=+0{,}00861$, $J[\tilde T,\tilde v]=-0{,}00861$ | reversibel (Kelvin) |

**Regel:** Dissipative Kreuzkopplungen sind symmetrisch (Onsager), reaktive antisymmetrisch (Casimir). Die symmetrischen erzeugen Entropie und gehören zur Rolle *transversal*. Die antisymmetrischen verschieben verlustfrei und gehören zur Rolle *tangential*.

### 3.6 Das Minimalsystem

Die Gleichungen, genau so gerechnet, mit $L=C=m=\ell=1$:

$$
\begin{aligned}
\dot Q_{\rm el}&=I,\\
L\dot I&=\underbrace{-Q_{\rm el}/C-G_yv}_{\rm tan}\;\underbrace{-R(T,A)I-bI^3}_{\rm irr}\;\underbrace{+(a-R_{\rm rad})I+S_b(T-T_a)}_{\rm port},\\
\dot x&=v,\\
m\dot v&=\underbrace{-k(T,D,g)\,x+G_yI}_{\rm tan}\;\underbrace{-cv-b_mv^3}_{\rm irr}\;\underbrace{+(a_m-c_{\rm ac})v}_{\rm port},\\
C_{\rm th}\dot T&=\underbrace{-\kappa_TTxv}_{\rm tan}\;\underbrace{+RI^2+bI^4+cv^2+b_mv^4}_{\rm irr}\;\underbrace{-G(T-T_a)-S_bTI}_{\rm port}\;\underbrace{+Y\dot D}_{\rm par},\\
\dot D&=\varepsilon_D\,(Y/Y_{\rm ref})^{m_D},\qquad \dot A=f_0e^{-E_a/T}(1-A),\qquad \dot g=\varepsilon_g .
\end{aligned}
$$

Dabei ist ($G_y$ Gyrator-Konstante, $f_0$ Arrhenius-Vorfaktor [korr.])
$$
R=R_0(1+\alpha_R(T-T_0))(1+\beta_AA),\qquad k=k_{\rm el}(1-D)-\kappa_T(T-T_0)+mg/\ell,\qquad Y=\tfrac12k_{\rm el}x^2 .
$$

Energie und Entropie sind
$$
E=\frac{Q_{\rm el}^2}{2C}+\frac{LI^2}{2}+\frac{mv^2}{2}+\frac{K(D,g)\,x^2}{2}+C_{\rm th}T,\qquad K=k-T\partial_Tk,\qquad S=C_{\rm th}\ln T+\tfrac12\kappa_Tx^2 .
$$

Die Parameter sind $a=0{,}35$, $b=0{,}2$, $R_0=0{,}1$, $\alpha_R=0{,}5$, $R_{\rm rad}=0{,}05$, $\beta_A=3$, $k_{\rm el}=0{,}2$, $\kappa_T=0{,}05$, $c=0{,}02$, $c_{\rm ac}=0{,}01$, $a_m=0{,}1$, $b_m=0{,}2$, $G_y=0{,}03$ [korr.], $C_{\rm th}=20$, $G=1$, $T_a=T_0=1$, $S_b=0{,}05$ und $E_a=10$.

**Identitäten der Zerlegung.** Geprüft an 2000 Zufallszuständen (Teil 1):

| Identität | Ergebnis |
|---|---|
| $\max\lvert\nabla E\cdot F_{\rm tan}\rvert$ | $1{,}4\cdot10^{-15}$ |
| $\max\lvert\nabla S\cdot F_{\rm tan}\rvert$ | $8{,}3\cdot10^{-17}$ |
| $\max\lvert\operatorname{div}F_{\rm tan}\rvert$ in $(Q_{\rm el},I,x,v,S_{\rm th})$ [korr.] | $5{,}6\cdot10^{-11}$ (Differenzenquotient) |
| $\max\lvert\nabla E\cdot F_{\rm irr}\rvert$ | $7{,}1\cdot10^{-15}$ |
| $\min\ \nabla S\cdot F_{\rm irr}$ | $7{,}0\cdot10^{-6}\ \ge0$ |
| $\max\lvert\nabla E\cdot F_{\rm port}-(P_{\rm in}-P_{\rm grav}-P_{\rm out})\rvert$ | $1{,}3\cdot10^{-15}$ |
| $\max\lvert\nabla E\cdot F_{\rm par}-P_{\rm grav}\rvert$ | $1{,}5\cdot10^{-17}$ (die Schädigung ist energieneutral: $-Y\dot D$ elastisch, $+Y\dot D$ Wärme) |
| $\min\ \nabla S\cdot F_{\rm par}$ | $\ge0$ |

### 3.7 Die Struktur der Jacobimatrix

**Rotierender Block.** $J_{\rm tan}$ im Ruhepunkt, energienormiert ($Q_{\rm el}/\sqrt C,\sqrt LI,\sqrt kx,\sqrt mv$), bei $g=0{,}5$:
$$
J_{\rm tan}=\begin{pmatrix}0&1&0&0\\-1&0&0&-0{,}03\\0&0&0&0{,}8303\\0&0{,}03&-0{,}8303&0\end{pmatrix},\qquad \lVert J+J^{\!\top}\rVert=0 .
$$
Die Matrix ist exakt antisymmetrisch, und ihre Eigenwerte $\pm i\,0{,}82915$ und $\pm i\,1{,}00144$ sind die Normalmoden. Die Diagonalblöcke $\begin{pmatrix}0&\omega\\-\omega&0\end{pmatrix}$ sind die beiden Umläufe. Das antisymmetrische Nebendiagonalelement $\pm G_y/\sqrt{Lm}$ ist die Casimir-Kopplung. Auf dem Torus selbst ist $J_{\rm tan}$ zustandsabhängig, bleibt aber spurfrei: $\operatorname{tr}J_{\rm tan}=4\cdot10^{-12}$ in $S_{\rm th}$-Koordinaten. [korr.]

**Kontrahierender Block.** Momentan ist $\operatorname{tr}J_{\rm trans}=-0{,}2972$. Über den Attraktor gemittelt ist die Kontraktion quer zum Torus negativ definit, mit $\operatorname{diag}(-\mu_e,-\mu_m,-G_{\rm eff}/C_{\rm th})=(-0{,}1895,\,-0{,}0700,\,-0{,}0525)$. Das ist die Vorhersage aus den Schleifenverstärkungen (§4). Gemessen werden $-0{,}18255$, $-0{,}06648$ und $-0{,}05361$ (§6).

**Langsamer Block.** Die Blocknormen der vollen Jacobimatrix bei $g=0{,}5$ mit aktiven $\varepsilon$-Termen:

| Zeile \ Spalte | EM $(Q_{\rm el},I)$ | Mech $(x,v)$ | $T$ | langsam $(D,A,g)$ |
|---|---|---|---|---|
| EM | 1,414 | 0,030 | 0,0763 | 0,174 |
| Mech | 0,030 | 1,245 | 0,00402 | 0,0820 |
| $T$ | 0,0146 | 0,0189 | 0,0482 | 0,00458 |
| langsam | 0 | $2{,}1\cdot10^{-8}$ | $9{,}0\cdot10^{-5}$ | $1{,}3\cdot10^{-5}$ |

Die langsame *Zeile* ist um vier bis acht Größenordnungen kleiner. Das ist $\varepsilon$: Die schnellen Variablen verändern die Parameter kaum. Die langsame *Spalte* ist dagegen von der Ordnung 0,1, denn die Parameter formen die schnellen Blöcke stark um. Genau diese Asymmetrie macht eine Wirkung *parametrisch*: **starker Einfluss auf den Fluss bei vernachlässigbarer Rückwirkung.**

**Normalform.** Um den Torus gilt
$$
\dot\theta=\omega(p)+O(\rho_\perp),\qquad \dot\rho_\perp=J_\perp(p)\,\rho_\perp+O(\rho_\perp^2),\qquad \dot p=\varepsilon G(\theta,\rho_\perp,p)
$$
mit dem transversalen Jacobi-Block $J_\perp\prec0$ [korr.]. Das ist die Torus-Normalform, hier mit expliziten Einträgen.

---

## 4  Rückkopplungsschleifen

### 4.1 Der Kopplungsgraph und die Reihe $\sum q^n$

Die Knoten sind die Blöcke von $z$, die Kanten die Einträge $a_{ij}=\partial\dot z_i/\partial z_j$, bei statischen Schleifen die stationären Übertragungsfaktoren. Eine **Schleife** ist ein Zyklus $i_1\to i_2\to\dots\to i_1$. Ihre Verstärkung ist das Produkt der Kantenfaktoren:
$$
q=\prod_{\rm Kanten}a_{i_{k+1}i_k}.
$$
Ein Signal, das die Schleife $n$-mal durchläuft, trägt $q^n$ bei. Die Summe über alle Rückläufe ist
$$
\boxed{\;\sum_{n=0}^\infty q^n=\frac1{1-q}\;}\qquad\text{bzw. für das Netz:}\qquad \sum_{n=0}^\infty\mathbf A^n=(\mathbf 1-\mathbf A)^{-1}\ \ \text{konvergent}\iff\rho(\mathbf A)<1 .
$$
Das ist die Neumann-Reihe. Die Mason-Regel für Signalflussgraphen ist ihre Zyklenzerlegung. Für dynamische Schleifen ist $q=q(i\omega)$, und dann gelten drei Fälle:

| Bedingung | Regime |
|---|---|
| $q<0$ (statisch) bzw. Nyquist-Ortskurve fern von 1 | **stabilisierend**: Antwort um $1/(1+\lvert q\rvert)$ gedämpft |
| $0<q<1$ | verstärkend, aber konvergent |
| $q\ge1$ bzw. $q(i\omega)=1$ (Barkhausen) | Reihe divergiert: **Aufschaukeln**, statisch Durchgehen, dynamisch Hopf-Einsatz eines Zyklus |

Eine Laufzeit $t_L$ in der Schleife dreht die Phase, $q\to q\,e^{-i\omega t_L}$. Dadurch kann eine gegenkoppelnde Schleife mitkoppelnd werden. Das ist die Formel zur „Phasenverdrehung durch Zwischenebenen“.

### 4.2 Die Schleifen der Tafeln

| Schleife | Zyklus (Tafel) | Verstärkung $q$ | Befund (Skript) |
|---|---|---|---|
| **S1** Lorentz/Induktion | EM-Feld → Strom → Mechanik → EM (Tafel 3, III) | $q_{LI}(\omega)=\dfrac{G_y^2\omega^2}{Lm(\omega^2-\omega_e^2)(\omega^2-\omega_m^2)}$ [korr.] | $q_{LI}(\omega_\pm)=1{,}000000000$: Die Normalmoden $\omega_\pm$ sind genau die Barkhausen-Punkte. $\omega=0{,}95$: $q=-0{,}03911$, $\sum q^n=0{,}96236$; $\omega=0{,}5$: $q=+0{,}00068$ |
| **S2** Joule-Material | EM-Feld → Wirbelstrom → Wärme → $\sigma_{\rm el}(T)$ [korr.] → EM (Tafel 3, IV/V) | $q=\dfrac{\partial P/\partial T}{G}$; stromgetrieben $+\alpha RI_0^2/G$, spannungsgetrieben $-\alpha P/G$, selbsterregt $-\dfrac{(a-R_{\rm rad})\,2\alpha_RR_0}{3bG}$ | stromgetrieben $q=0{,}5$: $\Delta T=1{,}00000=\Delta T_0/(1-q)$, Teilsummen $0{,}5;0{,}75;0{,}875;0{,}9844;1{,}0000$; $q=0{,}9$: Faktor 10; $q=1{,}2$: $50\,\Delta T_0$ nach $t=11{,}17\,C/G$ (Durchgehen). Selbsterregt: $q=-0{,}0500$, $\Delta T=0{,}22100/1{,}05=0{,}21048$ |
| **S3** Selbsterregung | Bewegung → Strömung → Kraft → Bewegung (Tafel 1, Strömung); aktives Element ↔ Kreis | $q_{\rm self}(A)=\dfrac{a}{R_{\rm tot}+\tfrac34bA^2}$, am Grenzzyklus $q=1$ | $q_{\rm self}(A^*)=1{,}000000$ für EM und Mechanik; $\lambda_\perp=\tfrac{aA^*}{2L}\,\partial_A\ln q=-\mu$: $-0{,}18948$ bzw. $-0{,}07000$ |
| **S4** thermo-mechanisch | Bewegung → Wärme → Ausdehnung → Steifigkeit → Resonanz (Tafel 1, Wärme/Mechanik; Tafel 3, V) | $q_{TM}=\dfrac{\partial\lvert X\rvert^2}{\partial\omega_m}\,\dfrac{\partial\omega_m}{\partial T}\,\dfrac1{G_{\rm eff}}\,\dfrac{\partial P}{\partial\lvert X\rvert^2}$ | Beim selbsterregten Oszillator ist $\partial\lvert X\rvert^2/\partial\omega_m\approx0$: offene Kette, $T$ verschiebt nur $\omega_m$. Bei fremderregter Resonanz: thermische Bistabilität für $q_{TM}\ge1$ |
| **S5** Maxwell-Druck | Feld → Grenzfläche → Mechanik → Spalt → Feld (Tafel 3, I) | $q_{ES}=\dfrac{\varepsilon_0AV^2}{k(d-x)^3}$, $x=x_0/(1-q)$ | am Pull-in-Punkt $x=d/3$, $V_{PI}=\sqrt{8kd^3/27\varepsilon_0A}$: $q_{ES}=1{,}000000$ |
| **S6** Abstrahlung | Schwingung → Schall/Strahlung → Energieverlust → Schwingung (Tafel 1, Akustik, Abstrahlung) | $q<0$ (reiner Verlust) | liefert den Anteil $-R_{\rm rad}/2L$ bzw. $-c_{\rm ac}/2m$ zu $\lambda_\perp$ |
| **S7** Ermüdung | Energie → Zustand ($D$) → Struktur ($k$) → Energie (Tafel 1: „veränderte Reaktion beschleunigt weitere Schädigung“) | $\dot D\propto(1-D)^{-m}$ (bei fester Geschwindigkeitsamplitude: $Y\propto x^2\propto1/k$) | integrierende Mitkopplung: endliche Lebensdauer $t_f=(1-D_0)^{m+1}/((m+1)\dot D_0)=0{,}3333/\dot D_0$ für $m=2$ |
| **S8** Alterung | EM → Wärme → Chemie ($A$) → $R$ → EM (Tafel 1, Chemie) | Empfindlichkeit $\partial\ln\dot A/\partial T=E_a/T^2=6{,}8248$ (Modell) | gegenkoppelnd in der Amplitude, aber einseitig: Drift bis $\mu_e=0$, Zyklustod bei $A_c=0{,}6563$ (§6) |
| **S9** Rauschen | Zustand ($T$) → Rauschen → Phase/Amplitude → Zustand (Tafel 1, Rauschen) | $D_\perp\propto\gamma k_BT$, Banddicke $\propto\sqrt{D_\perp/\lvert\lambda_\perp\rvert}$ | stabil, solange $\lambda_\perp<0$; die Banddicke wächst mit $T$ und divergiert bei $\lambda_\perp\to0$ |
| **S10** Gravitation | Gravitation → Vorspannung → Steifigkeit → Resonanz (Tafel 1) | Rückweg zu $g$: $q=0$ | **offene Kette**, deshalb parametrisch |
| **S11** Kriechen | Spannung → Kriechdehnung → Spannungsumlagerung → Spannung (Tafel 1, Struktur ⇄ Zustand) | $q<0$ (Relaxation) | stabilisierend, irreversibel |

**Netzform von S2** (Teil 5). Die Knoten sind $(\Delta T,\Delta R,\Delta P)$ mit $\mathbf A=\begin{pmatrix}0&0&1/G\\\alpha_RR_0&0&0\\0&\partial P/\partial R&0\end{pmatrix}$. Der Spektralradius ist $\rho(\mathbf A)=0{,}36840=|q|^{1/3}$. Die Differenz $\lVert(\mathbf 1-\mathbf A)^{-1}-\sum_{n<60}\mathbf A^n\rVert$ beträgt $2\cdot10^{-16}$, und die Gesamtverstärkung von $\Delta P_0$ nach $\Delta T$ ist $0{,}95238=(1/G)/(1-q)$.

**Das Vorzeichen gehört zur Schließung, nicht zum Glied.** S2 ist in allen drei Fällen dieselbe physikalische Kette $R(T)$. Wird der Strom von außen festgehalten, dann ist $q>0$, und oberhalb von 1 folgt thermisches Durchgehen. Setzt der Kreis seine Amplitude selbst (Selbsterregung), dann ist $q=-0{,}05$ und die Schleife stabilisiert. Ob eine Schleife aufschaukelt, ist also eine Eigenschaft der **gesamten Topologie des Kopplungsgraphen**, nicht des Widerstands. Das ist die Beobachter-Regel in Formelform.

### 4.3 Wann ein Torus-Attraktor entsteht

In Schleifensprache lautet die Bedingung „tangentiale Freiheit + transversale Kontraktion + langsame Parametrik, konsistent gekoppelt“ so:

1. **Jeder Zyklus $i$ hat eine Selbstschleife am Barkhausen-Punkt mit fallender Steigung:**
   $$q_{ii}(A_i^*)=1,\qquad \partial_A\ln q_{ii}\big|_{A_i^*}<0\ \Rightarrow\ \lambda_{\perp,i}=\frac{a_iA_i^*}{2L_i}\,\partial_A\ln q_{ii}<0 .$$
   Transversale Kontraktion ist nichts anderes als die *Sättigung der Schleifenverstärkung*.
2. **Die Kreuzschleifen zwischen den Zyklen bleiben unter der Rastschwelle:** $|\Delta\omega(p)|>\Delta\omega_R(G_y,\mu,A^*)$ (Rastbreite $\Delta\omega_R$) [korr.]. Dann bleiben beide Phasen frei, $\lambda_{\parallel,1}=\lambda_{\parallel,2}=0$.
3. **Die langsamen Schleifen konvergieren:** $\rho(\mathbf A_{\rm slow})<1$, bei Zeitskalentrennung $\varepsilon\ll|\lambda_\perp|\ll\omega$.
4. **Liouville:** $F_{\rm tan}$ ist divergenzfrei, also gilt $\sum\lambda=\langle\operatorname{div}F_{\rm trans}\rangle<0$.

Ist (2) verletzt, dann rastet der Torus ein: $T^2\to T^1$, eine Phasenrichtung wird kontrahierend. Ist (1) für einen Zyklus verletzt ($q_{ii}<1$ für alle $A$), dann stirbt dieser Zyklus (inverse Hopf-Bifurkation), und es folgt ebenfalls $T^2\to T^1$, aber durch Verlust statt durch Bindung. Ist (3) verletzt, dann bricht die Zeitskalentrennung, und der Torus wird von seinen eigenen Parametern fortgetragen (Durchgehen, Ermüdungsbruch).

---

## 5  Gravitation genau

### 5.1 Die FTM-Größen

$$
g=\frac{GM}{r^2}=\frac{3c^2\,\Delta s}{r^2},\qquad \Delta s=\frac{GM}{3c^2},\qquad \frac{d\tau}{dt}=\sqrt{1-\frac{6\Delta s}{r}}\quad(6\Delta s=2GM/c^2=r_s).
$$

| Körper | $g$ [m/s²] | $\Delta s$ [m] | $d\tau/dt$ | $1-d\tau/dt$ |
|---|---|---|---|---|
| Erde | 9,8203 | $1{,}47834\cdot10^{-3}$ | 0,999999999304 | $6{,}96\cdot10^{-10}$ |
| Sonne | 274,20 | 492,208 | 0,999997877495 | $2{,}12\cdot10^{-6}$ |
| Sirius B (Masse und Radius [G]) | $4{,}016\cdot10^{6}$ | 501,068 | 0,999740793253 | $2{,}59\cdot10^{-4}$ |
| Neutronenstern 1,4 M☉, 12 km | $1{,}290\cdot10^{12}$ | 689,092 | 0,809601233782 | 0,190 |

Kontrolle: $\Delta s_\oplus=1{,}4783$ mm und $r_0=6\Delta s_\oplus=8{,}870$ mm.

### 5.2 Kanal 1: gemeinsame Frequenzverschiebung (Zeitdilatation)

Alle tangentialen Zyklen am Ort $r$ laufen, von fern gesehen, langsamer:
$$
\omega_i^{\rm obs}=\omega_i^{\rm lok}\,\frac{d\tau}{dt}\qquad\text{für alle }i .
$$
Das ist eine parametrische Frequenzverschiebung, aber eine **gemeinsame**, und sie ist eine Umparametrisierung der Zeit, $t\mapsto\tau$. Daraus folgt:
$$
W^{\rm obs}=\frac{\omega_m^{\rm obs}}{\omega_e^{\rm obs}}=\frac{\omega_m^{\rm lok}}{\omega_e^{\rm lok}}=W^{\rm lok}.
$$
Die Rotationszahl ist invariant. Selbst beim Neutronenstern ($1-d\tau/dt=0{,}19$) bleibt $W^{\rm obs}-W^{\rm lok}=-1{,}1\cdot10^{-16}$, also Rundungsniveau. **Zeitdilatation allein kann den Torus weder einrasten lassen noch zerstören.** Sie verschiebt alle Umläufe gleich und lässt die Topologie $T^2$ unberührt. Im Grenzfall $r\to6\Delta s$ geht $d\tau/dt\to0$: Von außen erscheinen alle Zyklen eingefroren, lokal ist der Torus unverändert. Das ist die Beobachter-Regel an der Gravitation. Die Frequenzverschiebung ist eine Relation zwischen Orten, keine lokale Veränderung der Zyklen.

Der Gradient ergibt die Pound-Rebka-Form $\Delta\omega/\omega=g\,h/c^2$: Für $h=22{,}5$ m sind das $2{,}4585\cdot10^{-15}$.

### 5.3 Kanal 2: differentielle Modulation (Vorspannung → Steifigkeit → Resonanz)

Die Gravitation erzeugt eine Vorspannung $N=Mg$. Für eine Zugspannung (Pendel, gespannte Saite) gilt $k_{\rm eff}=k_{\rm el}+mg/\ell$, für eine Druckspannung (invertiertes Pendel, Stütze) $k_{\rm eff}=k_{\rm el}-mg/\ell$. Damit ist
$$
\omega_m(g)=\sqrt{\frac{k_{\rm el}(T,D)\pm mg/\ell}{m}},\qquad \frac{\partial\ln\omega_m}{\partial\ln g}=\pm\frac{mg}{2\ell\,k_{\rm eff}},\qquad \frac{\partial\omega_{\rm EM}}{\partial g}=0 .
$$
Für die Saite gilt $\omega_n=\tfrac{n\pi}{\ell}\sqrt{N/\mu_\ell}$ und damit $\partial\ln\omega_n/\partial\ln g=Mg/2N$. Die Rotationszahl $W(g)=\omega_m(g)/\omega_{\rm EM}$ ändert sich also. Mit FTM-Größen:
$$
\boxed{\;\omega_m^{\rm obs}(r)=\sqrt{\frac{k_{\rm el}}{m}\pm\frac{3c^2\Delta s}{\ell\,r^2}}\;\sqrt{1-\frac{6\Delta s}{r}},\qquad
W^{\rm obs}(r)=\frac{1}{\omega_{\rm EM}}\sqrt{\frac{k_{\rm el}}{m}\pm\frac{3c^2\Delta s}{\ell\,r^2}}\;}
$$
Der zweite Faktor (Kanal 1) kürzt sich heraus. Nur der erste (Kanal 2) wirkt auf die Topologie.

**Größenordnung** (Teil 8): Hebt man das System um $h=1000$ m an, dann gilt $\delta g/g=-3{,}14\cdot10^{-4}$. Beim Pendel folgt $\delta\omega_m/\omega_m=-1{,}57\cdot10^{-4}$, beim EM-Zyklus 0. Der gemeinsame Zeitdilatationsfaktor ändert sich nur um $1{,}09\cdot10^{-13}$. Im Labor ist der differentielle Kanal also **$1{,}44\cdot10^9$-mal stärker**. Im Minimalsystem beträgt $\partial\ln\omega_m/\partial\ln g=0{,}3626$ bei $g=0{,}5$.

### 5.4 Wie Gravitation Kopplungen und zugängliche Zyklen verändert

- **Frequenzen:** über $\omega_m(g)$, wie oben.
- **Kopplungen:** Die Durchbiegung $\delta_b=Mg/k$ verändert Spalte und Überdeckungen, also $G_y=B\ell_g(\delta_b)$, $C(\delta_b)$ und $L(\delta_b)$. Die Rastbreite $\Delta\omega_R$ skaliert mit $G_y$ (im Minimalsystem in $g$ gemessen: Fenster $[0{,}76;0{,}86]$ bei $G_y=0{,}03$). [korr.] Außerdem ändert die Vorspannung die Arbeitstemperatur über Reibung und damit $\mu_e(T)$.
- **Zugängliche Zyklen $n(p)$:** Es gibt drei Wege.
  1. *Einrasten.* $W(g)$ wandert in die Rastzunge um $W=1$, und es folgt $T^2\to T^1$ durch Bindung. Im Minimalsystem liegt die Resonanzmitte bei $g=0{,}8105$, das Rastfenster bei $g\in[0{,}76;0{,}86]$.
  2. *Knicken.* Bei Druckvorspannung geht $k_{\rm eff}\to0$ bei $g_c=k_{\rm el}\ell/m-\kappa_T\Delta T\,\ell/m=0{,}1895$. Dann gilt $\omega_m\to0$, der mechanische Umlauf verschwindet, und es folgt $T^2\to T^1$ durch Verlust.
  3. *Horizont.* $r\to6\Delta s$: gemeinsames Einfrieren (Kanal 1), keine Topologieänderung.
- **Langsamkeit und adiabatische Invariante:** Die parametrische Arbeit der Gravitation ist $P_{\rm grav}=\tfrac12m\dot g\,x^2/\ell$ (im Lauf von Teil 2 im Mittel $1{,}47\cdot10^{-5}$). Ohne Dissipation ist sie genau die Arbeit, die die Wirkung $J=E_m/\omega_m$ konstant hält. Eine langsame Gravitationsänderung verschiebt den Torus also auf die Nachbarfaser gleicher Wirkung. Das ist das *Deformieren, ohne lokal zu zerstören* des Autors, als Adiabatensatz.

---

## 6  Durchrechnung des Minimalsystems

### 6.1 Energiebilanz und Entropieproduktion (Teil 2)

Gerechnet wurde mit RK4, $dt=0{,}01$, $t=1000$. Alle langsamen Terme waren aktiv: $\varepsilon_D=10^{-5}$, $f_0=0{,}05$ [korr.], $\dot g=10^{-4}$, und $g$ lief von 0,5 bis 0,6.

| Größe | Wert |
|---|---|
| $\max\lvert E(t)-E(0)-\int(P_{\rm in}-P_{\rm out})dt\rvert$ | $5{,}2\cdot10^{-8}$ |
| Durchsatz $\int(P_{\rm in}+P_{\rm out})dt$ | 471,853 |
| **relativer Bilanzfehler** | $\mathbf{1{,}1\cdot10^{-10}}$ (Forderung $<10^{-6}$) |
| $\max\lvert S(t)-S(0)-\int(\Pi-\Phi_S)dt\rvert$ | $2{,}7\cdot10^{-10}$ |
| $\min_t\Pi(t)$ | $3{,}89\cdot10^{-6}\ \ge0$ |
| Endzustand | $T=1{,}1948$, $D=1{,}371\cdot10^{-3}$, $A=1{,}214\cdot10^{-2}$ |

**Leistungsflüsse**, Zeitmittel über $t=0\dots1000$:

| Kanal | Mittel |
|---|---|
| DC-Quelle → EM | 0,2164 |
| Strömung → Mechanik | 0,02183 |
| EM → Wärme: Joule / Sättigung | 0,06928 / 0,1154 |
| Mechanik → Wärme: innere Dämpfung / Sättigung | 0,004367 / 0,01541 |
| EM → außen (Abstrahlung) | 0,03091 |
| Mechanik → außen (Schall) | 0,002183 |
| Wärme → Umgebung | 0,2005 |
| EM → Mechanik (Lorentz $G_y\langle Iv\rangle$ [korr.]) | $3{,}49\cdot10^{-4}$ |
| Gravitation → Mechanik (parametrisch) | $1{,}47\cdot10^{-5}$ |
| Seebeck → EM / Peltier aus dem Wärmeknoten | $5{,}5\cdot10^{-6}$ / $3{,}0\cdot10^{-5}$ |
| Mechanik → Wärme reversibel ($\kappa_TTxv$) | $2{,}6\cdot10^{-6}$ |
| Schädigung $Y\dot D$ | $7{,}4\cdot10^{-8}$ |

Die großen Flüsse laufen quer: Quelle → Zyklus → Wärme → Umgebung. Die tangentialen Flüsse (Lorentz, Kelvin) mitteln auf dem freien Torus nahezu zu null.

### 6.2 Torus-Attraktor und Lyapunov-Signatur (Teil 6)

Fünfdimensionaler schneller Raum $(Q_{\rm el},I,x,v,T)$, $D,A,g$ eingefroren ($\varepsilon=0$), $t=3000$:

| | $g=0{,}5$ (frei) | $g=0{,}80$ (eingerastet) |
|---|---|---|
| Lyapunov-Spektrum | $+0{,}00004,\ -0{,}00001,\ -0{,}05361,\ -0{,}06648,\ -0{,}18255$ | $-0{,}00000,\ -0{,}00904,\ -0{,}05537,\ -0{,}12176,\ -0{,}14887$ |
| Signatur | $(0,0,-,-,-)$: $T^2$ | $(0,-,-,-,-)$: $T^1$ |
| $\sum\lambda$ / $\langle\operatorname{div}F\rangle$ | $-0{,}30261$ / $-0{,}30261$ | $-0{,}33505$ / $-0{,}33505$ |
| $\bar\omega_{\rm EM}$, $\bar\omega_{\rm mech}$ | 0,999111, 0,829245 | 1,004428, 1,004435 |
| Rotationszahl $W$ | 0,829983 | 1,000007 |
| $\lvert\langle e^{i\psi}\rangle\rvert$ (Phasendifferenz) | 0,0434 | 0,9999 |
| Poincaré-Schnitt $\theta_{\rm EM}=0$: belegte $\theta_{\rm mech}$-Bins | 100,0 % (397 Durchgänge) | 1,4 % (399 Durchgänge) |
| Lorentz-Nettofluss $G_y\langle Iv\rangle$ [korr.] | $+0{,}00019$ | $+0{,}01130$ |
| $T^*$ | 1,2017 | 1,2077 |

Die **Vorhersage aus den Schleifen** (§4, S2/S3) lautete $\lambda_{\perp}=-0{,}1895,\,-0{,}0700,\,-0{,}0525$ und $T^*=1{,}21048$. Gemessen wurden $-0{,}18255,\,-0{,}06648,\,-0{,}05361$ und $T^*=1{,}2017$. Die Abweichungen von wenigen Prozent kommen von der Gyratorkopplung und den höheren Harmonischen.

Im freien Fall gilt $\lambda_\parallel\approx0$ auf zwei Richtungen, auf allen anderen $\lambda_\perp<0$. Der Poincaré-Schnitt füllt einen geschlossenen Kreis: Das ist der $T^2$-Attraktor. Beim Einrasten wird eine Phasenrichtung kontrahierend mit $\lambda=-0{,}00904$. Der Schnitt zieht sich auf einen Punkt zusammen, und die Summe der Exponenten stimmt in beiden Fällen auf fünf Stellen mit $\langle\operatorname{div}F\rangle$ überein (Liouville, §3.4).

### 6.3 Kollaps $T^2\to T^1$ unter dem Gravitationsparameter (Teil 7)

**Statischer Scan** ($g$ fest, Start außerhalb der Rast):

| $g$ | $\omega_{m,0}$ | $\rho$ | $\lvert\langle e^{i\psi}\rangle\rvert$ |
|---|---|---|---|
| 0,20 | 0,6239 | 0,62406 | 0,021 |
| 0,50 | 0,8305 | 0,83022 | 0,056 |
| 0,70 | 0,9433 | 0,94378 | 0,198 |
| 0,74 | 0,9640 | 0,96611 | 0,446 |
| 0,75 | 0,9694 | 1,00056 | 0,907 |
| 0,76 | 0,9744 | 1,00003 | 0,9997 |
| 0,80 | 0,9949 | 1,00001 | 0,9999 |
| 0,84 | 1,0147 | 1,00000 | 0,9998 |
| 0,86 | 1,0244 | 1,00007 | 0,965 |
| 0,87 | 1,0294 | 1,03091 | 0,511 |
| 0,90 | 1,0436 | 1,04673 | 0,235 |
| 1,00 | 1,0903 | 1,09462 | 0,108 |

Das Rastfenster ist $g\in[0{,}76;\,0{,}86]$, entsprechend $\omega_{m,0}\in[0{,}974;\,1{,}024]$. Vor dem Einrasten wächst $|\langle e^{i\psi}\rangle|$ stetig von 0,02 auf 0,45: Der Torus wird **dünn**, bevor er kollabiert. Die Bahn verweilt immer länger nahe der künftigen geschlossenen Kurve. Das ist das *$T^2$ dünn machen* des Autors als messbare Größe. Beide Oszillatoren schwingen im Rastfenster weiter ($\max|I|\ge1{,}065$, $\max|v|\ge0{,}527$). Der Kollaps ist also eine Bindung, kein Zyklustod.

**Langsame Rampe** ($\dot g=\pm5\cdot10^{-5}$, Fenster $\Delta t=250$):

| $g$ auf | $\rho_{\rm auf}$ | $g$ ab | $\rho_{\rm ab}$ |
|---|---|---|---|
| 0,7375 | 0,95683 | 0,8625 | 1,02785 |
| 0,7500 | 0,97063 | 0,8500 | 0,99945 |
| 0,7625 | 1,00001 | 0,8375 | 0,99873 |
| 0,8000 | 1,00190 | 0,8000 | 0,99709 |
| 0,8500 | 1,00221 | 0,7500 | 0,99943 |
| 0,8625 | 1,00082 | 0,7375 | 0,99930 |
| 0,8750 | 1,00113 | 0,7250 | 0,97778 |
| 0,8875 | 1,01922 | 0,7125 | 0,95312 |

Als Rastkriterium in der Rampe gilt $|\rho-1|<5\cdot10^{-3}$. Im Rastzustand wandert nur die Rastphase $\psi^*(g)$ adiabatisch mit, was die kleinen Abweichungen $\pm(1\dots3)\cdot10^{-3}$ mit dem Vorzeichen der Rampe erklärt. Ein einziger Phasensprung pro Fenster ergäbe dagegen $\ge0{,}025$.

| Richtung | Rastfenster |
|---|---|
| Rampe aufwärts | $g\in[0{,}7625;\,0{,}8750]$ |
| Rampe abwärts | $g\in[0{,}7375;\,0{,}8500]$ |
| statisch | $[0{,}76;\,0{,}86]$ |

In beiden Richtungen bleibt der Rastzustand etwa 0,0125 bis 0,015 über den statischen Rand hinaus erhalten. Der Weg durch den Parameterraum bestimmt also, wo $T^1$ endet: Das Fenster verschiebt sich mit der Fahrtrichtung.

**Bistabilitätstest.** Statische Läufe, die im Rastzustand von $g=0{,}80$ starten, ergeben bei $g=0{,}72;\,0{,}74;\,0{,}87;\,0{,}88;\,0{,}89$ die Werte $\rho=0{,}95480;\,0{,}96606;\,1{,}03194;\,1{,}03711;\,1{,}04145$. Außerhalb des statischen Fensters gibt es **keinen** koexistierenden Rastattraktor. Die Rampen-Hysterese ist daher ein **Nachlauf** beim langsamen Durchfahren der Sattel-Knoten-Bifurkation (Geisterbahn der verschwundenen Rastkurve). Sie ist keine Bistabilität. Das Gedächtnis steckt hier in der Lage des Zustands relativ zum „Geist“ des verschwundenen Attraktors, also in der Geometrie des Flusses, und es klingt mit der Zeit ab.

**Alterungsweg** ($g=0{,}5$, $A$ fest): Der EM-Umlauf stirbt, sobald $\mu_e(A)=a-R_{\rm rad}-R(T^*,A)$ das Vorzeichen wechselt.

| $A$ | $\mu_e$ | $\max\lvert I\rvert$ | $\max\lvert v\rvert$ |
|---|---|---|---|
| 0,30 | +0,0991 | 0,854 | 0,746 |
| 0,50 | +0,0419 | 0,576 | 0,726 |
| 0,60 | +0,0150 | 0,357 | 0,709 |
| 0,65 | +0,0019 | 0,136 | 0,690 |
| 0,70 | −0,0132 | 0,056 (Restanregung über $G_y$ [korr.]) | 0,683 |
| 0,80 | −0,0435 | 0,056 | 0,682 |

Vorhergesagt war $A_c=0{,}6563$, gemessen liegt der Übergang zwischen 0,65 und 0,70. Das ist der zweite Weg $T^2\to T^1$: kein Einrasten, sondern der Verlust eines Umlaufs (inverse Hopf-Bifurkation). Oberhalb von $A_c$ läuft der EM-Kreis nur noch als getriebener Mitläufer des mechanischen Zyklus.

### 6.4 Fluktuations-Dissipations-Beziehung (Teil 9)

Langevin-Oszillator mit $\langle nn\rangle=2ck_BT\delta$ [korr.], 4000 Realisierungen:
$$
\langle\tfrac12mv^2\rangle/(\tfrac12k_BT)=1{,}0062,\qquad \langle\tfrac12kx^2\rangle/(\tfrac12k_BT)=1{,}0070 .
$$
Die Gleichverteilung wird innerhalb der Statistik erreicht: Die Rauschstärke ist durch die Dämpfung festgelegt. Real ergibt Johnson-Nyquist bei 1 kΩ und 300 K $S_V=1{,}6568\cdot10^{-17}$ V²/Hz, also 4,070 nV/√Hz.

---

## 7  Befunde

1. **Die drei Rollen des Autors sind die GENERIC-Struktur.** Tangential ist die reversible Klammer $\mathcal L\nabla E$ (energie- und entropieneutral, divergenzfrei). Transversal ist die dissipative Klammer $\mathcal M\nabla S$ zusammen mit den Toren (energieneutral, entropieerzeugend bzw. austauschend). Parametrisch ist der $\varepsilon$-langsame Rest. Im Minimalsystem ist diese Zerlegung exakt und auf $10^{-15}$ geprüft.
2. **Onsager und Casimir teilen die Kreuzkopplungen auf die Rollen auf.** Symmetrische Kreuzkopplungen ($L_{ij}=L_{ji}$) sind dissipativ und transversal. Antisymmetrische (Gyrator, Kelvin-Thermoelastik) sind reaktiv und tangential.
3. **Verriegelung braucht Dissipation.** $\operatorname{div}F_{\rm tan}=0$ erzwingt, dass die beim Einrasten auftretende tangentiale Kontraktion ($-0{,}00904$) aus $F_{\rm trans}$ stammt. Mit reinen Hamilton-Kopplungen bleibt $T^2$ immer $T^2$.
4. **Netto-Energiefluss zwischen Zyklen setzt Phasenbindung voraus.** Auf dem freien Torus pendelt die Lorentzkopplung Energie hin und her, im Mittel fließen $1{,}9\cdot10^{-4}$. Im Rastzustand fließen $1{,}13\cdot10^{-2}$, rund 60-mal mehr. Klasse A (Energie) zwischen zwei Umläufen setzt Klasse B (Phase) voraus.
5. **Transversale Kontraktion ist die Steigung der Schleifenverstärkung am Barkhausen-Punkt**, $\lambda_\perp=\tfrac{aA^*}{2L}\partial_A\ln q$. Die vorhergesagten Werte stimmen auf wenige Prozent mit den gemessenen Exponenten überein.
6. **Das Vorzeichen einer Schleife gehört der Schließung, nicht dem Glied.** Dieselbe Kette $R(T)$ ist mit $q>0$ durchgehend (stromgetrieben) oder mit $q=-0{,}05$ stabilisierend (selbsterregt).
7. **Parametrisch heißt: offene Kette.** Alle Schleifen durch $p$ haben auf schnellen Skalen $q\approx0$ (langsame Zeile der Jacobimatrix $\le9\cdot10^{-5}$), während $p$ die schnellen Blöcke stark umformt (langsame Spalte $\sim0{,}1$).
8. **Zeitdilatation ändert keine Topologie.** Sie ist ein gemeinsamer Faktor aller tangentialen Frequenzen und damit eine Umparametrisierung der Zeit. Topologisch wirksam ist nur der differentielle Kanal (Vorspannung → Steifigkeit), im Labor $1{,}44\cdot10^9$-mal stärker. $T^2\to T^1$ entsteht durch Einrasten ($g\in[0{,}76;0{,}86]$) oder Knicken ($g_c=0{,}1895$).
9. **Gedächtnis ist Koordinate und Flussgeometrie.** Die Vorgeschichte steckt in $D$, $A$ und $T$, also in einseitig laufenden Koordinaten des aktuellen Punkts. Die Rampe zeigt eine zweite Form: Das Rastfenster verschiebt sich mit der Fahrtrichtung um 0,0125 bis 0,015 (Nachlauf an der Geisterbahn des Sattel-Knotens), obwohl statisch kein zweiter Attraktor existiert. Pfadabhängigkeit entsteht hier ohne Speicher, allein aus der Langsamkeit des Flusses nahe einer verschwundenen Invariante.

---

## 8  Selbstcheck: Beobachter-Regel

| Prüffrage | Antwort |
|---|---|
| Wird irgendwo eine Wirkung als lokales Objekt oder Substanz eingeführt? | Nein. Jede Kopplung ist ein Term von $F$, eine Kante im Kopplungsgraphen oder ein Eintrag von $\mathcal L$ bzw. $\mathcal M$. |
| Ist Gedächtnis ein Speicher? | Nein. $D$, $A$ und $\boldsymbol\varepsilon_p$ sind Koordinaten des gegenwärtigen Zustands mit einseitiger Dynamik. Ihr Einfluss ist $\partial F/\partial D\ne0$, eine Verformung des Vektorfelds. Die Rampen-Hysterese ist ein Nachlauf an der Geisterbahn des verschwundenen Rastattraktors. Einen koexistierenden Attraktor gibt es statisch nicht (Bistabilitätstest). |
| Wird Rauschen als äußere Störung verdinglicht? | Nein. Es ist über die FDT dieselbe Größe $\gamma$ wie die Dissipation, also eine Eigenschaft des Zustands. |
| Ist die Gravitation eine Kraft auf die Zyklen? | Nur als Vorspannung. Topologisch wirkt sie als Parameter $p$, und ihre Zeitdilatation ist eine Relation zwischen Beobachtungsorten ($\rho^{\rm obs}=\rho^{\rm lok}$). |
| Wird ein Aufschaukeln einem Bauteil zugeschrieben? | Nein. Das Vorzeichen von $q$ hängt von der Schließung ab (Befund 6). |
| Mystik? | Keine. Verwendet werden Maxwell, Lorentz, Poynting, Fourier, Kelvin, Onsager/Casimir, GENERIC, Einstein/Johnson/Callen-Welton, Arrhenius, Kachanov/Miner, Barkhausen/Nyquist und die Neumann-Reihe. |
| Stammen alle Zahlen aus dem Skript? | Ja, aus der Ausgabe der Rechnung. |

**Offen für den Autor:**
- Die Rastbreite $W(G_y,\mu)$ [korr.] ist bisher nur numerisch bestimmt. Eine geschlossene Formel für reaktiv gekoppelte van-der-Pol-Oszillatoren fehlt noch.
- Die thermische Bistabilität (S4) ist nur als Bedingung formuliert und noch nicht gerechnet.
- Die Feldform (§3.1) ist aufgestellt, gerechnet wurde aber nur das modale Minimalsystem.
