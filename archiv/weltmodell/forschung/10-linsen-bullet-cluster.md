# Gravitationslinsen quantitativ: Galaxien, Haufen, Bullet Cluster

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Umbenennung Index A des Operatorfelds → Op (A ist im Modell das Vektorpotential), Phantomdichte ρ_ph → gelesene Dichte ρ_X, Phantommasse/Phantomfeld → gelesene Masse/Zusatzanteil, ρ_mem → ρ_bar,mem (Baryonendichte an der Gedächtnislage), Faktor s → ζ_r; Linse „trifft die RAR bis 1 Mpc auf 2 % (KiDS)“ → Werte [G], im Modell folgt die Linse um isolierte Galaxien bis 100 kpc der Relation, weiter außen senkt der Rest der Mode ξ sie; η = 1 „in stationären Systemen“ → für Kugeln und außerhalb von Scheiben, in Scheiben innen η > 1; a₀(z) = cH(z)/2π bzw. „höheres a₀“ → trägt nicht; „Massenrest offen, Faktor ~2“ → im Modell schließt der stoßfreie Rest der Mode ξ den Faktor; „es fehlt eine kovariante Form mit η = 1“ → im Modell Ψ = Φ, γ = β = 1, c_GW = c; junge Gezeitenzwerge mit verminderter Massendiskrepanz → liegen im Modell auf der Relation; „unsicher/prüfen“ → [G]. Die Notiz rechnet mit a₀ = 1,2·10⁻¹⁰.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors (Punkt II.2).
> Rechnung: `rechnungen/linsen.py` (jede Zahl unten stammt aus dieser Rechnung).
> Größen: $\Delta s = GM/(3c^2)$, $r^* = \sqrt{GM/a_0}$, Operator
> $\nabla^2\Phi = \nabla\cdot[\nu(|\nabla\Phi_N|/a_0)\,\nabla\Phi_N]$, $\nu = 1/(1-q)$, $q = e^{-\sqrt{a_N/a_0}}$, $a_0 = 1{,}2\cdot10^{-10}$ m/s² (gemessener Wert; der Modellwert $cH_0/2\pi = 1{,}08\cdot10^{-10}$ m/s² bei $H_0 = 70$ km/s/Mpc erscheint nur in Abschnitt 3 als Vergleich).
> Leitregel: Die Zusatzwirkung ist Struktur, keine Quelle. In der Rechnung erscheint sie nur als $\nabla\cdot[(\nu-1)\nabla\Phi_N]$, das heißt als Umformung des baryonischen Feldes. Eine zweite Dichte kommt nicht vor.

**Ergebnis in drei Sätzen.**
1. Die Galaxie-Galaxie-Linsen verlangen $\eta = M_\text{dyn}/M_\text{lens} = 1$: Der Operator wirkt auf beide Potentiale. Dieses Operatorfeld trägt im Folgenden den Index Op.
2. Haufen behalten damit denselben Rest wie MOND. Er liegt beim Faktor 1,5–3 in der Masse und beim Faktor 2,4–4,3 im Einstein-Radius-Bereich. [korr.] Im Modell schließt diesen Rest der stoßfreie Rest der Mode ξ der Schließungsregel: Die Operatorform fasst 21–71 % der Mode, der Rest trifft die Linsenmasse auf 1–5 %.
3. Im Bullet Cluster folgt das Konvergenzmaximum des momentanen Operatorfeldes dem Gas des Subhaufens. Der Gedächtnis-Nachlauf mit $1/\lambda_m \gtrsim 170$–$260$ Myr zieht es zu den Galaxien; offen bleibt ein Restversatz von 15–35 %, den Massenrest schließt im Modell der stoßfreie Rest der Mode ξ. [korr.]

---

## 1 Linsengleichung im Modell

### 1.1 Ablenkung

Im schwachen Feld gilt $ds^2 = -(1+2\Phi/c^2)c^2dt^2 + (1-2\Psi/c^2)\,d\mathbf x^2$. Die Ablenkung ist

$$
\hat\alpha(b) = \frac{1}{c^2}\int\nabla_\perp(\Phi+\Psi)\,dl .
$$

Mit $\Phi=\Psi=-GM/r$ ergibt das $4GM/(c^2b)$. In der Sprache der Faltung heißt das:

$$
\alpha_N(b) = \frac{12\,\Delta s}{b}.
$$

Kontrolle an der Sonne: $12\Delta s_\odot/R_\odot = 1{,}7512''$ (GR: $1{,}751''$).

**Takt und Zusatzstrecke.** Die Brechzahl aus der Taktverzögerung $c_\text{eff} = c\sqrt{1-6\Delta s/r}$ liefert $6\Delta s/R_\odot = 0{,}8756''$, also die Hälfte. Die andere Hälfte kommt aus dem räumlichen Anteil $\Psi$. Im FTM ist das die Rolle von $\Delta s$ als *Zusatzstrecke*. Takt (Φ) und Zusatzstrecke (Ψ) tragen je $6\Delta s/b$ bei.

### 1.2 Der Operator wirkt auf beide Potentiale

Der Operator wirkt auf $\Phi$ und $\Psi$. Für eine Kugel gilt

$$
\hat\alpha_\text{Op} = \frac{2}{c^2}\int g_\perp\,dz \quad\text{mit}\quad g=\nu g_N,\qquad \eta(r) = 1 .
$$

### 1.3 Formeln in $\Delta s$ und $r^*$ (Punktmasse)

Setzt man $z = b\tan\theta$ und $x = b/r^*$, dann gilt $a_N/a_0 = (r^*/r)^2 = \cos^2\theta/x^2$ und damit

$$
\hat\alpha_\text{Op}(b) = \frac{12\,\Delta s}{b}\,I\!\Big(\frac{b}{r^*}\Big),\qquad
I(x) = \int_0^{\pi/2}\cos\theta\;\nu\!\Big(\frac{\cos^2\theta}{x^2}\Big)d\theta .
$$

Für $b\gg r^*$ geht die Ablenkung gegen einen konstanten Wert:

$$
\hat\alpha_\text{Op} \simeq \frac{6\pi\,\Delta s}{r^*} + \frac{6\,\Delta s}{b},\qquad
\frac{6\pi\Delta s}{r^*} = \frac{2\pi\sqrt{GMa_0}}{c^2} = 2\pi\frac{v_\text{flach}^2}{c^2}.
$$

Das ist die Ablenkung einer isothermen Kugel.

Zahlen für $M = 10^{11}\,M_\odot$ ($\Delta s = 329$ AE, $r^* = 10{,}78$ kpc, $\hat\alpha_\infty = 0{,}5755''$); $\alpha_N = 12\Delta s/b$ ist die Ablenkung der Baryonen allein:

| $b/r^*$ | $I$ | $\alpha_N$ ["] | $\alpha_\text{Op}$ ["] |
|---|---|---|---|
| 0,3 | 1,165 | 1,221 | 1,422 |
| 1 | 2,135 | 0,366 | 0,782 |
| 10 | 16,21 | 0,0366 | 0,594 |
| 100 | 157,6 | 0,00366 | 0,577 |

---

## 2 Galaxie-Galaxie-Linsen: $\Delta\Sigma(R)$ von 10 kpc bis 1 Mpc

**Modell.** Die Baryonen einer isolierten Galaxie sind als Hernquist-Kugel angesetzt, mit $a = 1{,}5\,\text{kpc}\,(M/10^{10}M_\odot)^{0,3}$ (Annahme). Heißes Gas ist nicht enthalten. Aus der projizierten Linsenmasse folgen $\Delta\Sigma = \bar\Sigma(<R)-\Sigma(R)$ und $g_\text{obs}\approx 4G\Delta\Sigma$ in der SIS-Umrechnung. Verglichen wird mit der Kurve $g_\text{RAR} = \nu(g_\text{bar}/a_0)\,g_\text{bar}$.

$\Delta\Sigma$ in M☉/pc²:

| $R$ [kpc] | $10^{10}$: Baryonen (Newton) | Op | $10^{11}$: Baryonen (Newton) | Op |
|---|---|---|---|---|
| 10 | 22,9 | 73,6 | 172,3 | 266,7 |
| 50 | 1,19 | 14,7 | 11,1 | 48,9 |
| 100 | 0,307 | 7,35 | 2,97 | 23,9 |
| 300 | 0,035 | 2,45 | 0,346 | 7,81 |
| 1000 | 0,0032 | 0,734 | 0,0316 | 2,33 |

Verhältnis $4G\Delta\Sigma/g_\text{RAR}$, also Linse gegen die RAR-Kurve, über alle drei Massen $10^{10}$, $3\cdot10^{10}$ und $10^{11}\,M_\odot$:

| | 10 kpc | 100 kpc | 1 Mpc |
|---|---|---|---|
| Op ($\eta = 1$) | 1,00–1,02 | 1,00–1,01 | 1,00 |
| Baryonen allein (Newton) | 0,31–0,66 | 0,04–0,13 | 0,004–0,014 |

**Vergleich mit der Beobachtung.** Brouwer et al. 2021 (A&A 650, A113; KiDS-1000, isolierte Galaxien) finden, dass die Linsen-RAR der Kurve $g_\text{obs}(g_\text{bar})$ aus den Rotationskurven folgt, weit über die Rotationskurven hinaus bis zu Radien von etwa 1 Mpc [G].
- Zwei Einschränkungen von dort [G]: Frühe Typen liegen höher als späte, und bei den kleinsten $g_\text{bar}$ spielt der unbekannte heiße Gashalo eine Rolle.
- Auch Mistele et al. 2024 berichten flache Kreisgeschwindigkeiten aus KiDS-Linsen bis etwa 1 Mpc [G].

**Befund.**
- Das Operatorfeld mit $\eta = 1$ reproduziert die RAR in der Linse auf 2 %. [korr.] Das gilt für das stationäre Operatorfeld ohne die Mode ξ. Nach der Schließungsregel des Modells ist um isolierte Galaxien wie die Milchstraße die Mode ξ bis 100 kpc ganz übergegangen, die Linse folgt dort der Relation; weiter außen senkt der Rest der Mode die Linse unter die Relation (bei N* = 2,4 um 3–7 % bei 200 kpc, 17–44 % bei 300 kpc; bei N* = 0,3 nicht messbar). Ein Treffer auf 2 % bis 1 Mpc verlangte kleines N*.
- Auch bei Haufen stimmen hydrostatische Röntgenmassen und Linsenmassen auf etwa 10–30 % überein (hydrostatischer Bias, z. B. „Weighing the Giants“) [G].

**Folge für die Trennbeobachtung I des Autors.** Innerhalb des Operators gilt für Kugeln und außerhalb von Scheiben $\eta = 1$; in Scheiben innen ($R \lesssim 3{,}5\,R_d$) ist $\eta > 1$. [korr.] Die Unterscheidung „Struktur statt Quelle“ zeigt sich deshalb in der Zeit- und Pfadabhängigkeit (Abschnitt 5).

---

## 3 Massereicher Galaxienhaufen: Einstein-Radius und Rest

**Baryonisches Modell (Annahme, realistische Größenordnung).**
- Gas als β-Modell mit $\beta = 2/3$, $r_c = 150$ kpc und $M_\text{gas}(<1\,\text{Mpc}) = 10^{14}\,M_\odot$, abgeschnitten bei 3 Mpc.
- BCG mit ICL: $1{,}5\cdot10^{12}\,M_\odot$, $a = 20$ kpc.
- Mitgliedsgalaxien: $10^{13}\,M_\odot$, $a = 300$ kpc.
- Lage $z_l = 0{,}3$.
- Die „beobachtete“ Linsenmasse vertreten NFW-Profile typischer massereicher Linsenhaufen ($M_{200}=10^{15}$, $c=4$–6; $2\cdot10^{15}$, $c=5$).

| $r$ bzw. $R$ [kpc] | $g_N/a_0$ | $\nu$ | $M_\text{dyn,Op}/M_\text{bar}$ (3D) | Rest 3D (NFW/Op) | $M_{2D,\text{Op}}/M_{2D,\text{bar}}$ | Rest 2D (NFW 1e15,c4 / Op) |
|---|---|---|---|---|---|---|
| 50 | 0,553 | 1,91 | 1,91 | 2,74 | 2,50 | 2,42 |
| 100 | 0,368 | 2,20 | 2,20 | 3,15 | 2,77 | 2,29 |
| 250 | 0,288 | 2,41 | 2,41 | 2,66 | 3,13 | 1,87 |
| 500 | 0,207 | 2,74 | 2,74 | 2,15 | 3,67 | 1,49 |
| 1000 | 0,125 | 3,36 | 3,36 | 1,59 | 4,58 | 1,09 |

**Linsenmasse in 500 kpc.**
- Baryonen allein (Newton): $7{,}5\cdot10^{13}\,M_\odot$.
- **Op: $2{,}76\cdot10^{14}\,M_\odot$.**
- NFW-Referenzen: $4{,}1$–$6{,}9\cdot10^{14}$.

Das Operatorfeld Op erreicht damit die beobachtete Größenordnung ($10^{14}$–$10^{15}$) in 500 kpc. Es bleibt der Faktor 1,5 bis 2,5.

**Einstein-Radius ($z_s = 2$).** Op ergibt $\theta_E = 1{,}57''$ ($R_E = 7{,}0$ kpc), und dieser Ring kommt praktisch nur vom BCG. Die NFW-Referenzen ergeben $9{,}6''$ bis $29{,}3''$. Beobachtete massereiche Haufen haben typisch 10–50″ [G].

Nötige Masse für einen gegebenen Einstein-Radius:

| $\theta_E$ | $R_E$ [kpc] | nötig [M☉] | Op | Rest | nötig/Baryonen |
|---|---|---|---|---|---|
| 10″ | 44,5 | $1{,}47\cdot10^{13}$ | $6{,}0\cdot10^{12}$ | 2,44 | 5,98 |
| 20″ | 89,1 | $5{,}87\cdot10^{13}$ | $1{,}84\cdot10^{13}$ | 3,19 | 8,69 |
| 30″ | 133,6 | $1{,}32\cdot10^{14}$ | $3{,}58\cdot10^{13}$ | 3,68 | 10,5 |
| 45″ | 200,4 | $2{,}97\cdot10^{14}$ | $6{,}89\cdot10^{13}$ | 4,31 | 13,0 |

**Hat der Operator denselben Rest wie MOND?** Ja, und zwar quantitativ.
- Bei $M_{2D}(<500\,\text{kpc})$ liegt der Rest mit $\nu_\text{FTM}$ bei 1,49, mit $\nu_\text{simple}$ bei 1,48 und mit $\nu_\text{standard}$ bei 1,69.
- Bei 100 kpc sind es 2,29, 2,28 und 2,69.
- $\nu = 1/(1-q)$ ist im Haufenbereich ($g_N/a_0\approx0{,}1$–0,5) fast gleich der „simple“-Funktion. Deshalb erbt das Modell den bekannten Befund: Sanders 2003 (MNRAS 342, 901), Pointecouteau & Silk 2005 und Angus, Famaey & Buote 2008 (MNRAS 387, 1470) finden den Faktor ~2 [G].
- Der Rest ist **zentral konzentriert**. Er beträgt 3D 3,2 bei 100 kpc und 1,6 bei 1 Mpc, im Einstein-Radius-Bereich 2,4–4,3. [korr.] Im Modell schließt ihn der stoßfreie Rest der Mode ξ der Schließungsregel: Die Operatorform fasst 21–71 % der Mode, der stoßfreie Rest (79 % bei 100 kpc bis 29 % bei 1,5 Mpc) trifft die Linsenmasse auf 1–5 %.
- Setzt man $a_0 = cH_0/2\pi$ ($H_0 = 70$ km/s/Mpc), wird der Rest bei 500 kpc etwas größer (1,55). Mit $a_0(z) = cH(z)/2\pi$ ($=1{,}26\cdot10^{-10}$ bei $z=0{,}3$) wird er etwas kleiner (1,45). Die Änderung ist klein. [korr.] $a_0(z) = cH(z)/2\pi$ trägt nicht: Der Vakuumteil bliebe zu jeder Zeit bei 22 %, im Strahlungszeitalter wie $\Delta N_\text{eff} = 2{,}09$ (gemessen $N_\text{eff} = 2{,}99 \pm 0{,}17$). Im Modell ist $a_0$ seit $z \approx 10$ fest; $\approx cH_0/2\pi$ gilt nur heute.

---

## 4 Bullet Cluster 1E 0657-56

### 4.1 Rechnung

**Verfahren.**
- Der volle 3D-Operator wird auf einem Gitter gerechnet: 192×192×96 Zellen zu 20×20×40 kpc, isolierte Randbedingung über Nullauffüllung und FFT-Faltung mit $1/r$.
- Die gelesene Dichte (Umrechnung der Regel durch einen Beobachter, der ungeänderte Regeln annimmt; kein Stoff) ist $\rho_X = -\nabla\cdot[(\nu-1)\,\mathbf g_N]/(4\pi G)$. Sie wird entlang der Sichtlinie aufsummiert. [korr.]
- Die Konvergenz ist $\kappa = \Sigma/\Sigma_\text{crit}$ mit $z_l = 0{,}296$ und $z_s = 1$, also $\Sigma_\text{crit} = 2{,}84\cdot10^{9}\,M_\odot/\text{kpc}^2$.
- **Validierung** an einer Kugel gegen die 1D-Rechnung: 2,6 % Abweichung bei 100 kpc, 0,1 % bei 250 und 500 kpc.
- Es entsteht keine negative gelesene Masse. [korr.]

**Baryonen (Modellannahme).** Geometrie gerundet nach Clowe et al. 2006 (ApJ 648, L109), Abb. 1, mit ±50 kpc Unsicherheit. Massen in der Größenordnung von Clowe et al. 2006 und Bradač et al. 2006 (ApJ 652, 937), nicht exakt übernommen.

| Komponente | Lage $x$ [kpc] | Masse [M☉] | Form |
|---|---|---|---|
| Haupt-Galaxien + BCG | 0 | $1{,}0\cdot10^{13}$ + $1{,}0\cdot10^{12}$ | Hernquist 250 / 25 kpc |
| Hauptgas | 200 | $1{,}5\cdot10^{14}$ | β = 2/3, $r_c$ = 200 kpc |
| Bullet-Gas | 560 | $2{,}0\cdot10^{13}$ | β = 2/3, $r_c$ = 60 kpc |
| Sub-Galaxien + BCG | 720 | $5\cdot10^{12}$ + $6\cdot10^{11}$ | Hernquist 150 / 20 kpc |

Gas zu Sternen ist 10,2. Das Feld ist tief im Übergangsbereich: $\nu = 2{,}59$ (Haupt-BCG), 5,73 (Hauptgas), 8,76 (Mitte), 3,75 (Bullet-Gas), 2,22 (Sub-BCG).

### 4.2 Folgen die Maxima dem Gas oder den Galaxien? (momentanes Feld, $1/\lambda_m \to 0$)

Konvergenzmaxima auf der Stoßachse. Die Karten sind wie eine Linsenrekonstruktion geglättet (Gauß-σ):

| Karte | σ = 50 kpc: Haupt / Sub | σ = 100 kpc: Haupt / Sub |
|---|---|---|
| Baryonen allein | 15 / 560 | 134 / 544 |
| **Op (momentan)** | **19 / 547** | **203 / 520** |
| Galaxien bzw. Gas liegen bei | 0 bzw. 200 / 720 bzw. 560 | |

- **Subhaufen:** Das Maximum des momentanen Feldes Op liegt **beim Bullet-Gas** (547 bzw. 520 kpc), nicht bei den Sub-Galaxien (720). Beobachtet liegt das Linsenmaximum bei den Galaxien (Clowe et al. 2006: Versatz Linse–Plasma mit 8σ). Dorthin führt der Gedächtnis-Nachlauf (Abschnitt 5).
- **Haupthaufen:** Das Ergebnis ist unentschieden und hängt am Sternmodell.
  - Mit Sternen ×1 ($2{,}6\cdot10^{12}$ in 100 kpc) liegt das Maximum bei σ = 50 kpc bei den Galaxien (19 kpc), bei σ = 100 kpc beim Gas (203 kpc).
  - Mit Sternen ×0,25 liegt es auch bei σ = 50 kpc beim Gas (259 kpc). Diese Sternmasse ($6{,}4\cdot10^{11}$ in 100 kpc) liegt nahe an den ~$5\cdot10^{11}$, die Clowe et al. angeben [G].
  - Erst Sterne ×2 ziehen auch den Sub-Peak zu den Galaxien (697).
- Der Grund: Der Operator verstärkt die Baryonen mit dem Faktor $\nu$. Dominiert das Gas die Baryonen, dominiert es auch $\kappa$. Der Zusatzanteil (ν − 1) allein, als gelesene Masse, hat sein Sub-Maximum bei 537 kpc, also ebenfalls am Gas. [korr.]

### 4.3 Wie hoch ist der Rest?

Aperturmassen in 100 kpc um die BCGs:
- Op: $1{,}51\cdot10^{13}$ (Haupt-BCG) und $1{,}34\cdot10^{13}$ (Sub-BCG).
- Die Linsenmassen sind ~$2{,}8$ und $2{,}0\cdot10^{13}$ (Clowe et al. 2006, Tab. 2) [G].
- Damit ergibt sich ein Rest von **1,9 (Haupt) und 1,5 (Sub)**. Am Bullet-Gas ist er 1,1, weil das momentane Feld Op dort zu viel Masse hat.
- Das passt zum MOND-Befund von Angus et al. 2007 (ApJ 654, L13) [G]: Es fehlt ein Faktor ~2, und diese fehlende Masse muss **stoßfrei** sein, also bei den Galaxien sitzen.

---

## 5 Gedächtnis-Nachlauf: Das Operatorfeld läuft der Zustandsgeschichte nach

### 5.1 Formulierung

Der Gedächtniskern $e^{-\lambda_m(t-t')}$ wird auf die Schwerpunktbewegung jeder kohärenten Komponente angewandt:

$$
\mathbf u(t) = \lambda_m\int_{-\infty}^{t}e^{-\lambda_m(t-t')}\,\mathbf v(t')\,dt',\qquad \dot{\mathbf X}_\text{mem} = \mathbf u,
$$

$$
\nabla^2\Phi = 4\pi G\rho_\text{bar} + \nabla\cdot\Big[\big(\nu(|\nabla\Phi_N[\rho_\text{bar,mem}]|/a_0)-1\big)\,\nabla\Phi_N[\rho_\text{bar,mem}]\Big].
$$

Die Quelle ist unverändert $\rho_\text{bar}$. Die **Struktur** wird aus der Gedächtnislage gebildet: $\rho_\text{bar,mem}$ ist die Baryonendichte derselben Baryonen, an die Gedächtnislage $\mathbf X_\text{mem}$ ihrer Bahngeschichte verschoben, keine eigene Dichte. [korr.] Eine neue Substanz tritt nicht auf. Der Nachlauf hat drei Eigenschaften:
- **Stationär** gilt $\rho_\text{bar,mem}=\rho_\text{bar}$. Galaxien, Rotationskurven, Sonnensystem und die Abschnitte 2–3 bleiben unverändert.
- **Galaxien** sind stoßfrei. Ihre Geschwindigkeit bleibt konstant, also gilt $\mathbf X_\text{mem} = \mathbf X_\text{gal}$.
- **Gas** wird beim Stoß abgebremst. Seine Struktur bleibt bei der alten Bahn, also bei den Galaxien:

$$
\frac{X_\text{mem,gas}-X_\text{gal}}{X_\text{gas}-X_\text{gal}} = \varphi(\lambda_m t),\qquad
\varphi_\text{sofort}(x) = 1-\frac{1-e^{-x}}{x},\qquad
\varphi_\text{konst}(x) = 1-\frac2x+\frac{2(1-e^{-x})}{x^2}.
$$

### 5.2 Ergebnis für den Bullet Cluster

Grundlage: $t = 150$ Myr und 4500 km/s, also 690 kpc Durchgangsweg. Gemessen wird $\beta$, der Versatz κ–Galaxien geteilt durch den Versatz Gas–Galaxien. $\beta = 0$ heißt „bei den Galaxien“.

| $\varphi$ | $\beta_\text{Haupt}$ | $\beta_\text{Sub}$ | Sterne ×0,25: $\beta_\text{Haupt}$ / $\beta_\text{Sub}$ |
|---|---|---|---|
| 0 (unendliches Gedächtnis) | 0,05 | 0,15 | 0,13 / 0,35 |
| 0,2 | 0,07 | 0,23 | 0,24 / 0,60 |
| 0,5 | 0,10 | 0,52 | 0,81 / 0,83 |
| 1 (kein Nachlauf, momentanes Feld Op) | 0,09 | 1,08 | 1,29 / 1,15 |

**Nötige Nachlaufzeit.** Das Kriterium ist $\beta\le 0{,}25$ bei σ = 50 kpc. Es ergibt $\varphi\le0{,}243$, also:

$$
1/\lambda_m \;\ge\; 256\ \text{Myr (sofortige Abbremsung)}\quad\text{bzw.}\quad 167\ \text{Myr (konstante Abbremsung)}.
$$

Das ist etwa die Kreuzungszeit (156 Myr) bis zum 1,7-Fachen davon. Folgen:
- Der Nachlauf zieht das Maximum zu den Galaxien. Mit $1/\lambda_m\to\infty$ bleibt aber ein **Restversatz** von $\beta \approx 0{,}15$ (Sterne ×1) bzw. 0,35 (Sterne ×0,25). Das Gas liegt als Newton-Quelle weiterhin am wahren Ort.
- Bei σ = 100 kpc wird $\beta\le0{,}25$ auch mit unendlichem Gedächtnis nicht erreicht.
- Mit Sternen ×0,25 wird $\beta\le0{,}4$ erst bei $1/\lambda_m\ge1{,}9$ Gyr erreicht.
- **Massenrest.** Die Aperturmassen bei $\varphi=0$ ergeben einen Rest von 1,6 (Haupt) und 1,3 (Sub). [korr.] Im Modell schließt ihn der stoßfreie Rest der Mode ξ der Schließungsregel: Die Operatorform fasst 21–71 % der Mode, der Rest trifft die Linsenmasse auf 1–5 %; 79 % der gelesenen Masse folgen den Galaxien, 21 % dem abgebremsten Gas mit Nachlauf.
- Behandelt man den Rest ebenfalls als Gedächtnisstruktur an den stoßfreien Komponenten (Faktor $\zeta_r$ auf den Zusatzanteil), sinkt $\beta_\text{Sub}$ bei Sternen ×1 von 0,15 auf 0,12 ($\zeta_r=2$). [korr.] Die Aperturmassen erreichen dann $2{,}9$ und $2{,}6\cdot10^{13}$.

### 5.3 Prüfbare Vorhersage: $\beta$ hängt von der Zeit seit dem Stoß ab

Die Tabelle gilt für die Bullet-Geometrie mit Sternen ×1, σ = 50 kpc und sofortiger Abbremsung:

| $t$ [Gyr] | $1/\lambda_m$ = 300 Myr | 600 Myr | 1000 Myr | 3000 Myr | ohne Nachlauf | ΛCDM (stoßfreie DM) |
|---|---|---|---|---|---|---|
| 0,15 | 0,23 | 0,19 | 0,17 | 0,15 | 1,08 | 0 |
| 0,5 | 0,55 | 0,30 | 0,23 | 0,17 | 1,08 | 0 |
| 1,0 | 0,84 | 0,55 | 0,34 | 0,20 | 1,08 | 0 |
| 2,0 | 0,97 | 0,84 | 0,66 | 0,27 | 1,08 | 0 |

**Vorhersage.** Der Versatz κ–Gas ist $(1-\beta(t))\cdot d_\text{Gas-Gal}$. Mit Gedächtnis-Nachlauf **wandert das Linsenmaximum mit der Zeit seit dem Stoß vom Galaxienort zum Gas**. In ΛCDM bleibt es bei den Galaxien. Soll $\beta\le0{,}25$ auch in alten Verschmelzungen gelten, braucht man
- $1/\lambda_m \ge 854$ Myr bei $t = 0{,}5$ Gyr,
- $1/\lambda_m \ge 1195$ Myr bei $t = 0{,}7$ Gyr,
- $1/\lambda_m \ge 1708$ Myr bei $t = 1$ Gyr.

Kandidaten, alle Zeitmarken [G]:
- **Musket Ball**, $t\approx0{,}7$ Gyr (Dawson 2013): $\beta = 0{,}73$ für $1/\lambda_m = 300$ Myr, $0{,}27$ für 1 Gyr.
- **El Gordo**, $t\approx0{,}5$ Gyr.
- **Stapelanalyse:** Harvey et al. 2015 (Science 347, 1462; 72 Kollisionen) finden die Linsenmasse im Mittel bei den Sternen. Ihr Wert für $\beta$ ist mit null verträglich [G]. Schon der Restversatz 0,15–0,35 des Nachlaufs ist damit in Spannung, sofern die Messgrößen vergleichbar sind.

**Zweite Folgerung des Nachlaufs.** [korr.] Liefe der Nachlauf auch auf Galaxienbahnen, zeigten junge Gezeiten-Zwerggalaxien, die jünger als $1/\lambda_m$ sind, eine verminderte Massendiskrepanz, weil ihre Gedächtnisstruktur noch auf der Bahn des Mutter-Gases läge. Das trägt nicht: Die Rückkehr ist eine stationäre Schleife am Ort (τ_D = 79 Mrd. Jahre, länger als das Weltalter). Das Gedächtnis ist eins (Mode ξ, Gewicht w, Übergaberate λ_m = ω_orb/N*, N* ≈ 0,3–2,4) und wirkt auf den Zustand der Regel, nicht auf die einzelne Bahn. In Galaxien gibt es keinen messbaren Nachlauf; einen Nachlauf mit 1/λ_m gibt es nur, wo ξ nicht ganz übergegangen ist (Haufen, frisch kollabierte Bereiche). Im Modell liegen junge Gezeitenzwerggalaxien schon auf der Beschleunigungsrelation; prüfbar an NGC 5291 (Bournaud et al. 2007, Lelli et al. 2015; Alter ~0,36 Gyr [G]).

---

## Übersicht

| Frage | Befund |
|---|---|
| Linsen | $\eta=1$: Der Operator wirkt auf Φ und Ψ. |
| Galaxie-Galaxie-Linsen bis 1 Mpc | Op trifft die RAR auf 2 %, in Übereinstimmung mit KiDS-1000 [G]. [korr.] Mit der Mode ξ folgt die Linse um isolierte Galaxien bis 100 kpc der Relation, weiter außen liegt sie bei N* = 2,4 darunter (3–7 % bei 200 kpc, 17–44 % bei 300 kpc). |
| Haufen | Linsenmasse in 500 kpc: $2{,}8\cdot10^{14}$, also im richtigen Bereich, aber Rest 1,5–3,2. Einstein-Radien um den Faktor 2,4–4,3 in der Masse zu klein. Derselbe Rest wie bei MOND. [korr.] Im Modell schließt ihn der stoßfreie Rest der Mode ξ (Linsenmasse auf 1–5 %). |
| Bullet Cluster, momentanes Feld | Maximum des Subhaufens liegt beim Gas. Rest ≈ 1,5–1,9 (Linsenwerte [G]). |
| Gedächtnis-Nachlauf | Mit $1/\lambda_m\gtrsim$ 170–260 Myr Maximum nahe den Galaxien. Restversatz 15–35 % bleibt offen; den Massenrest schließt im Modell der stoßfreie Rest der Mode ξ. [korr.] Klare Zeitvorhersage $\beta(t)$. |

Die Zusatzwirkung bleibt in allen Rechnungen Struktur: Sie ist die Divergenz eines umgeformten baryonischen Feldes. Der Bullet Cluster zeigt aber, dass die **Lage** dieser Struktur der Geschichte der stoßfreien Komponenten folgen muss und dass ihr **Betrag** mit dem Operatorfeld allein in Haufen um ~2 zu klein ist. [korr.] Im Modell trägt diese Geschichte die Mode ξ der Schließungsregel: Ihr stoßfreier Rest schließt den Faktor und trifft die Linsenmasse auf 1–5 %.

---

## 6 Offene Punkte als nächste Rechenschritte

1. **Metrische Form des Operators mit $\eta = 1$.** [korr.] Im Modell liegt sie vor: Die kovariante Wirkung gibt $\Psi = \Phi$, $\gamma = \beta = 1$ und $c_\text{GW} = c$; die Metrik trägt den wirbelfreien Teil der Regel, für Kugeln und außerhalb von Scheiben ist $\eta = 1$, in Scheiben innen $\eta > 1$. GW170817 und die Sonnensystem-Grenzen sind damit erfüllt. Die Lichtablenkung setzt sich aus Takt- und Zusatzstrecken-Anteil zusammen ($6\Delta s/b$ je Anteil).
2. **Echte Bullet-Daten.** Die Modellkomponenten sind durch die Chandra-Gaskarte und die publizierte κ-Karte (Clowe 2006, Bradač 2006) zu ersetzen, mit Pixelvergleich und Sichtlinien-Tiefe. Zu verbessern sind außerdem die Gitterauflösung ($h_z$ = 20 kpc) und die Box-Größe (Konvergenztest). Die unsicheren Tabellenwerte (Aperturmassen, Lagen, Sternmassen) sind nachzuschlagen.
3. **Haufenstichprobe.** Zu rechnen ist $\text{Rest}(r) = M_\text{lens}/M_\text{Op}$ für CLASH- bzw. LoCuSS-Haufen mit gemessenen Gasprofilen. Frage: Folgt der Rest den Galaxien (stoßfrei) oder dem Gas? Und hängt er von $z$ ab, wie über $a_0(z)=cH(z)/2\pi$ vorhergesagt (1,49 → 1,45 bei $z = 0{,}3$)? [korr.] $a_0(z) = cH(z)/2\pi$ trägt nicht (der Vakuumteil bliebe zu jeder Zeit bei 22 %, wie $\Delta N_\text{eff} = 2{,}09$); im Modell ist $a_0$ seit $z \approx 10$ fest, eine $z$-Abhängigkeit über $a_0$ gibt es nicht.
4. **Galaxienlinsen im Detail.** Zu rechnen ist $\Delta\Sigma$ mit heißem Gashalo, Zwei-Halo-Term und einem möglichen Außenfeld-Effekt des Operators, getrennt nach frühen und späten Typen. Der Vergleich mit den KiDS-Datenpunkten soll Punkt für Punkt erfolgen, nicht nur qualitativ.
5. **Gedächtnis-Nachlauf als Feldtheorie.** Die Transportgleichung für $\rho_\text{bar,mem}$ ist sauber zu formulieren (Kontinuität mit $\mathbf u$ statt $\mathbf v$, Trennung Bulk-Fluss/innere Bahnen). Danach stehen drei Prüfungen an: Stationarität exakt nachweisen, Gezeiten-Zwerge rechnen, $\beta(t)$ gegen die Stichprobe von Harvey et al. 2015 mit geschätzten Stoßzeiten fitten.
6. **Woher der Faktor ~2 in Haufen?** Im Modell gibt es zwei Kandidaten, beide ohne neue Substanz. Erstens ein höheres $a_0$ bzw. eine andere Summenform $\nu$ bei $g_N/a_0\sim0{,}1$–0,5. Das darf die Galaxien-RAR nicht verletzen, dort ist der Spielraum < 10 %. Zweitens eine Strukturkomponente, die an stoßfreie Freiheitsgrade koppelt (Test 5 des Autors: „kalt vs. heiß“). [korr.] Das höhere $a_0$ trägt nicht ($a_0$ ist seit $z \approx 10$ fest; an $H$ gebunden bliebe der Vakuumteil zu jeder Zeit bei 22 %, wie $\Delta N_\text{eff} = 2{,}09$). Den Faktor schließt im Modell die zweite Möglichkeit: der stoßfreie Rest der Mode ξ der Schließungsregel. Die Operatorform fasst 21–71 % der Mode, der Rest trifft die Linsenmasse auf 1–5 %.
