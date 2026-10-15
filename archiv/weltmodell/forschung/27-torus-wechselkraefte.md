# Torus-Wechselkräfte

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: „für die Galaxie fehlt die konkrete Identifikation“ (§1.4, Selbstcheck, Offen) → im Modell Schließung durch Translation: Rückkehr T im Zustandsraum, Rückkopplungssumme als Resolvente (I − qT)⁻¹, die Rückkehr ist eine stationäre Schleife am Ort. Umbenennungen (ohne Marke): Γ (Gruppe) → 𝔊; H₀ (Energie im Toy Model) → E₀; ξ_j (Rauschen) → n_j; η (relative Dicke) → h_r; η (Dämpfer) → d_D; s (Strukturkoordinate) → x_s, damit 𝓛(s) → 𝓛(x_s); ν, ν_S, ν_k (Alterungsraten) → k_A, k_{A,S}, k_{A,k}; χ_g (Hebel) → f_g; zusätzlich Γ = Bℓ (Wandlerfaktor) → Bℓ ausgeschrieben. [G] ergänzt: Stoffwerte Kupfer, T_C = 493 K.

> Arbeitsstand, **noch nicht im Gesamtmodell**. Wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/torus_wechselkraefte.py`

**Worum es geht.** Die drei Rollen des Autors (tangential, transversal, parametrisch) entsprechen der GENERIC-Zerlegung ż = 𝓛∇E + 𝓜∇S + ε·F_par. Diese Notiz fragt, wie die Wirkungen **als Kräfte am Torus** gegeneinander stehen: was den Torus zusammenhält, was seine Form festlegt und wann er dünn wird, einrastet oder zerfällt. Sechs Fragen:

1. Wie wird der Satz „Jeder Schnitt erzeugt Rand. Rand erzeugt Fluchtmoden. Schließung eliminiert Randbedingungen durch Identifikation.“ zu einem Satz über Spektren?
2. Welches Kräftegleichgewicht aus tangentialen, transversalen und parametrischen Anteilen legt Radius, Dicke und Frequenzverhältnis des Torus fest?
3. Wie sehen die sechs Übergänge des Dreiecks Energie–Struktur–Zustand als Kopplungsoperatoren aus, und wann stabilisiert sich der Kreis selbst?
4. Wie groß ist die Schleifenverstärkung von „Wärme verändert die Regeln, nach denen EM wirkt“ in einem realen Material?
5. Wie genau wirkt Gravitation als Vorspannung, im Labor und an kompakten Sternen?
6. Was folgt daraus Neues?

Die Beobachterfehler-Regel des Autors gilt durchgehend: Keine Wirkung wird als Objekt oder Substanz eingeführt. Jede Wirkung ist ein Generator, ein Projektor, eine Kante im Kopplungsgraphen oder ein Parameter des Flusses.

---

## 1  Schließung eliminiert Fluchtmoden

### 1.1 Der Satz als Operatoraussage

Sei $\tilde X$ ein offener Raum, etwa die Gerade $\mathbb R$. Die Identifikation $x\sim x+L$ ist die Wirkung der Gruppe $\mathfrak G=L\mathbb Z$ durch Translationen $T_{mL}$, und der geschlossene Raum ist $X=\tilde X/\mathfrak G=S^1_L$. Der Generator der Dynamik $\mathcal L$ (Diffusion $D\partial_x^2$, Liouville $-v\partial_x$ oder freier Schrödinger-Operator) ist translationsinvariant: $[\mathcal L,T_{mL}]=0$.

Die Rechnung geht über eine endliche Überlagerung, den Kreis der Länge $ML$, der für $M\to\infty$ in die Gerade übergeht. Auf ihm ist

$$
\boxed{\;\mathcal P=\frac1M\sum_{m=0}^{M-1}T_{mL},\qquad \Pi_\perp=I-\mathcal P\;}
$$

**Satz 1 (Schließung als Projektion).**
1. $\mathcal P$ ist ein orthogonaler Projektor ($\mathcal P^2=\mathcal P=\mathcal P^\ast$), und er vertauscht mit dem Generator: $[\mathcal P,\mathcal L]=0$.
2. Das Bild von $\mathcal P$ sind genau die $L$-periodischen Funktionen, also die Funktionen auf dem geschlossenen Kreis $S^1_L$. Die Dynamik auf $S^1_L$ ist die Einschränkung von $\mathcal L$ auf $\operatorname{Bild}\mathcal P$:
$$
\mathcal L_{\rm eff}=\mathcal P\mathcal L=\mathcal L\mathcal P=\mathcal P\mathcal L\mathcal P .
$$
3. In Fourier-Darstellung gilt $T_{mL}e^{ikx}=e^{ikmL}e^{ikx}$. $\mathcal P$ lässt daher nur die Wellenzahlen des Gitters $k\in(2\pi/L)\mathbb Z$ durch. $\Pi_\perp$ ist die Spektralprojektion auf alle übrigen $k$, insbesondere auf das langwellige Band $0<|k|<2\pi/L$.
4. Das Spektrum von $-D\partial_x^2$ geht dabei vom kontinuierlichen $[0,\infty)$ ohne normierbare Eigenfunktionen in das diskrete $\{D(2\pi n/L)^2\}$ über. Es behält den Nullmode (die Konstante, also die erhaltene Gesamtwahrscheinlichkeit) und öffnet eine Lücke $4\pi^2D/L^2$.

Das ist die Formel des Autors $\dot z=(I-\Pi_\perp)F(z)$ bzw. $L_{\rm eff}=\mathcal PL$ mit einem expliziten Projektor. Entscheidend ist Punkt 2: Weil $\mathcal P$ mit $\mathcal L$ vertauscht, kommt **kein Term zum Generator hinzu**. Das Vektorfeld bleibt, wie es ist, nur der Raum der zulässigen Zustände wird kleiner. Bindung ohne zusätzliche Quelle heißt also: Der Generator ist unverändert, und die Fluchtmoden liegen im Kern der Projektion. Eine Quelle würde dagegen den Generator ändern, $\mathcal L\to\mathcal L+\nabla\Phi_{\rm DM}\cdot\partial_v$. Das ist die Trennung „DM = Quellenänderung, B = Operatoränderung“ aus Tafel 6, jetzt als Rechenoperation.

**Rand gegen Schließung.** Das Intervall $[0,L]$ mit absorbierendem Rand hat ebenfalls ein diskretes Spektrum, aber ohne Nullmode: Alle Moden klingen ab, die langsamste mit $\lambda_1=\pi^2D/L^2$. Der Rand wirkt als Senke, über die alles entweicht. Genau das meint „Rand erzeugt Fluchtmoden“.

Rechnung, Teil 1a/1b (diskreter Laplace-Operator, $D=L=1$):

| Raum | kleinster Eigenwert | Moden unter der Kreislücke $4\pi^2D/L^2$ |
|---|---|---|
| Intervall, Rand absorbierend | $\lambda_1=9{,}7714$ (analytisch $\pi^2=9{,}8696$) | kein Nullmode |
| Kreis (identifiziert) | $\lambda_0=0$, dann $\lambda_{1,2}=39{,}4752$ ($4\pi^2=39{,}4784$) | nur die Konstante |
| offen, $\Lambda=10L$ / $100L$ / $1000L$ | $9{,}86\cdot10^{-2}$ / $9{,}87\cdot10^{-4}$ / $9{,}87\cdot10^{-6}$ | 20 / 200 / 2000 (Weyl: $2\Lambda/L$) |

Der Projektor auf dem Überlagerungskreis $M=8$ (32 Punkte pro $L$): $\lVert\mathcal P^2-\mathcal P\rVert=0$ und $\lVert[\mathcal P,\Delta]\rVert=0$. Rang $\mathcal P=32$, Rang $\Pi_\perp=224$. Das Spektrum von $\mathcal P\Delta\mathcal P$ auf $\operatorname{Bild}\mathcal P$ stimmt bis auf $3{,}9\cdot10^{-12}$ mit dem Kreisspektrum überein. Unter der halben Lücke liegt im Bild von $\mathcal P$ nur die Konstante, im Bild von $\Pi_\perp$ liegen 10 Moden. Die langsamste davon hat den Eigenwert $0{,}6168=(2\pi/ML)^2$: Das sind die langwelligen Fluchtmoden, und die Identifikation schneidet sie heraus.

### 1.2 Verbleib, Rückkehr, Wiederkehr

**Diffusion** (Teil 1c; Start in der Fenstermitte, Zeit in $L^2/D$):

| $t$ | offen: Anteil im Fenster $L$ | Intervall mit Rand: Überleben | Kreis | Rückkehrdichte offen | Rückkehrdichte Kreis |
|---|---|---|---|---|---|
| 0,1 | 0,736448 | 0,4744875 | 1 | 0,89206 | 1,03859 |
| 1 | 0,276326 | $6{,}59\cdot10^{-5}$ | 1 | 0,28209 | 1,00000 |
| 10 | 0,089021 | $1{,}74\cdot10^{-43}$ | 1 | 0,08921 | 1,00000 |
| 100 | 0,028204 | 0 | 1 | 0,02821 | 1,00000 |

Offen fällt die Rückkehrdichte wie $1/\sqrt{4\pi Dt}$ gegen null, auf dem Kreis strebt sie gegen $1/L$. Der **Aufenthaltsanteil** (Zeitmittel) im Fenster $a=0{,}1L$ um den Start beträgt auf dem Kreis 0,10071, 0,10001 und 0,10000 bei $t=10$, $10^3$ und $10^5$, strebt also gegen $a/L$. Offen sind es 0,01772, 0,00178 und 0,00018, gemäß $a/\sqrt{\pi Dt}$ gegen null.

**Liouville-Fluss** $\dot x=v$ (Teil 1d, Gauß-Ensemble $\sigma_v=1$): Offen liegen im Fenster $L$ nach $t=1$, 10 und 100 noch 0,38292, 0,03988 und 0,00399 des Ensembles. Auf dem Kreis bleibt alles im Raum, und die Dichte mischt sich zur Gleichverteilung ($\max|\rho L-1|=5{,}35\cdot10^{-9}$ bei $t=1$). Das ist Phasenmischung ohne Verlust. Das Koopman-Spektrum ist auf dem Kreis rein punktförmig ($iv\,2\pi n/L$), offen kontinuierlich.

**Poincaré und Kac.** Ein maßerhaltender Fluss auf einem Raum endlichen Volumens kehrt fast überall beliebig nahe zurück (Poincaré). Die mittlere Rückkehrzeit in eine Menge $A$ ist bei Ergodizität $\langle T_A\rangle=\mu(X)/\mu(A)$ (Kac). Gemessen an der Kreisdrehung um den goldenen Winkel: $|A|=0{,}01$ ergibt 100,001, $|A|=0{,}001$ ergibt 999,807. Die offene Translation kehrt nie zurück.

**Wellenpaket** (Teil 1e, $\hbar=m=1$): Auf dem Ring ist das Spektrum $E_n=2\pi^2n^2/L^2$ kommensurabel. Das Paket kehrt deshalb exakt zurück, bei $T_{\rm rev}=mL^2/\pi\hbar=0{,}31831$ mit $|\langle\psi_0|\psi(t)\rangle|^2=1{,}000000$. Bei $T_{\rm rev}/4$ gibt es eine Teilwiederkehr (0,500000). Auf der Geraden (Ring der Länge $64L$) bleibt der Überlapp bei $2{,}3\cdot10^{-19}$, und schon der ruhende Anteil zerfließt auf $1{,}57\cdot10^{-2}$. Allgemein gilt: Schließung macht das Spektrum diskret, damit ist die Dynamik fastperiodisch und kehrt beliebig nahe zurück. Exakte Wiederkehr braucht zusätzlich ein kommensurables Spektrum.

**Pólya** (Teil 1f): Die Irrfahrt auf $\mathbb Z^3$ ist transient. Das Watson-Integral ergibt $u=1{,}516386$, die Rückkehrwahrscheinlichkeit ist $1-1/u=0{,}340537$. Auf dem Torus $\mathbb Z_N^3$ kehren 100 % zurück, im Mittel nach 63,6 bzw. 216,8 Schritten für $N=4$ bzw. 6 (Kac: $N^3=64$ bzw. 216). In einer und zwei Dimensionen ist schon die offene Irrfahrt rekurrent, aber nur nullrekurrent (mittlere Rückkehrzeit unendlich). Die Schließung macht daraus **positive Rekurrenz** mit endlicher mittlerer Rückkehrzeit. Das ist die genaue Bedeutung von „gebunden“.

### 1.3 Bindung als Funktion der Schließungslänge

Wie stark die Schließung bindet, hängt nur von $L$ ab (Teil 1g, $D=1$):

| $L$ | $\kappa_{\rm var}=12D/L^2$ | $\kappa_{\rm gap}=4\pi^2D/L^2$ | Fluchtrate mit Rand $\pi^2D/L^2$ | Kac-Rückkehr (Gitter $h=0{,}01$) |
|---|---|---|---|---|
| 1 | 12,000 | 39,478 | 9,870 | 100 Schritte |
| 2 | 3,000 | 9,870 | 2,467 | 200 |
| 10 | 0,120 | 0,395 | 0,0987 | 1000 |
| 100 | 0,0012 | 0,00395 | 0,00099 | 10000 |

Dabei ist $\kappa_{\rm var}$ die Steifigkeit eines harmonischen Potentials, das dieselbe stationäre Varianz erzeugt wie der Kreis ($L^2/12$), und $\kappa_{\rm gap}$ die Rückstellrate, die dieselbe Relaxationszeit liefert. Das **effektive Potential** der Schließung ist also

$$
U_{\rm eff}(x)=\tfrac12\kappa_{\rm eff}x^2,\qquad \kappa_{\rm eff}=c\,\frac{D}{L^2},\quad c=12\ \text{(Varianz)}\ \text{bzw.}\ 4\pi^2\ \text{(Relaxation)} .
$$

Das Toy Model des Autors mit unendlichem Gedächtnis ($K\equiv1$) gehorcht derselben Form. Mit $GM=1$, $x_0=1$, $v_s=1{,}6$ und $E_0=0{,}28$ liegt der Umkehrpunkt bei $x_{\max}=5{,}3256$ ($\gamma=0{,}05$), 9,7475 ($\gamma=10^{-2}$), 26,2243, 77,5372 und 751,1085 ($\gamma=10^{-6}$). Das Verhältnis $\gamma\ell_B^2/2E_0$ mit $\ell_B=x_{\max}-x_0$ geht dabei über 1,67061 und 1,13619 gegen 1,00475:

$$
\boxed{\;\text{Bindungsstärke}\;\approx\;\frac{\text{Energieskala}}{(\text{Schließungslänge})^2}:\qquad \kappa_{\rm eff}\simeq c\,\frac{D}{L^2},\qquad \gamma\simeq\frac{2E_0}{\ell_B^2}\;}
$$

Eine engere Schließung bindet quadratisch stärker. Die Quelle bleibt dabei dieselbe.

### 1.4 Einordnung in den Torus und Grenze der Aussage

Am Torus-Bündel des Autors gibt es zwei Sorten Fluchtrichtungen, und sie werden auf verschiedene Weise stillgelegt:

- **tangential**: durch Identifikation ($\theta\sim\theta+2\pi$). Die Umlaufrichtungen sind kompakt, eine Bahn kann entlang des Torus nicht entweichen, die Rückkehrzeiten sind endlich (Kac). Das ist Satz 1.
- **transversal**: durch Kontraktion ($\lambda_\perp<0$, Abschnitt 2). Hier entfernt nicht die Topologie, sondern die Dissipation die Fluchtrichtungen.

„Lokale Freiheit entlang der Bahn, globale Gebundenheit durch Schließung“ heißt damit genau: Entlang des Torus bindet die Topologie, quer zu ihm die Dissipation. $\Pi_\perp$ des Autors umfasst beides.

**Grenze.** Satz 1 zeigt, *wie* eine Identifikation bindet und dass sie dabei den Generator nicht ändert. Für die Galaxie ist die Identifikation im Modell die Rückkehr: Schließung durch Translation, die Verschiebung $T$ der Geschichte um einen Rücklauf im Zustandsraum, kausal summiert mit dem Rückkehranteil $q=e^{-\sqrt{a_N/a_0}}$. Die Rückkopplungssumme $(I-qT)^{-1}$ ist die Resolvente dieser Translation, und $\mathcal P_q=(1-q)(I-qT)^{-1}$ geht für $q\to1$ in $\mathcal P$ über. Die Rückkehr ist eine stationäre Schleife am Ort. [korr.]

---

## 2  Die drei Rollen als Kräfte am Torus

### 2.1 Das Modell

Für Schritt 1–4 des Autors werden zwei selbsterregte Zyklen (EM-Phase und mechanische Phase) in Normalform geschrieben (Stuart-Landau), gekoppelt über ein gemeinsames Medium:

$$
\dot z_j=(\mu+i\omega_j)z_j-\beta|z_j|^2z_j+Ke^{i\varphi}(z_k-z_j)+\sqrt{2D}\,n_j,\qquad z_j=A_je^{i\theta_j}.
$$

Die Kopplung zerfällt selbst in zwei Rollen: $K_d=K\cos\varphi$ ist **dissipativ** (symmetrisch, Onsager) und $K_r=K\sin\varphi$ **reaktiv** (antisymmetrisch, Casimir). In Polarform mit $\psi=\theta_2-\theta_1$:

$$
\begin{aligned}
\dot A_1&=\underbrace{\mu A_1-\beta A_1^3}_{\text{Verstärkung/Sättigung}}\;\underbrace{-K_d(A_1-A_2\cos\psi)-K_rA_2\sin\psi}_{\text{Querkraft der Kopplung}}&&\text{(transversal)}\\
\dot\theta_1&=\omega_1(p)+K_d\frac{A_2}{A_1}\sin\psi+K_r\Big(\frac{A_2}{A_1}\cos\psi-1\Big)&&\text{(tangential)}\\
\dot p&=\varepsilon\,G(\theta,A,p):\ \mu(p),\ \omega_j(p),\ K(p),\ D(p)&&\text{(parametrisch)}
\end{aligned}
$$

Bei gleichen Amplituden fällt $K_r$ aus der Phasendifferenz heraus:

$$
\dot\psi=\Delta\omega-2K_d\sin\psi .
$$

### 2.2 Das Kräftegleichgewicht, das die Torusform festlegt

Drei Gleichgewichte bestimmen die Form des Torus:

| Gleichgewicht | Kräfte | Ergebnis |
|---|---|---|
| **Radius** (transversal) | Verstärkung $\mu$ gegen Sättigung $\beta A^2$ und dissipative Querkraft $K_d(1-\langle\cos\psi\rangle)$ | $A^{*2}=\dfrac{\mu-K_d\,(1-\langle\cos\psi\rangle)}{\beta}$; frei: $\langle\cos\psi\rangle=0$; eingerastet: $\cos\psi^*=\sqrt{1-(\Delta\omega/2K_d)^2}$ |
| **Dicke** (transversal gegen Rauschen) | Rauschen $D$ gegen Rückstellrate $\lvert\lambda_\perp\rvert=2\beta A^{*2}$ | $\operatorname{Var}A=\dfrac{D/2}{2\beta A^{*2}}$, relative Dicke $h_r=\sigma_A/A^*$ |
| **Frequenzverhältnis** (tangential) | Verstimmung $\Delta\omega(p)$ gegen dissipative Phasenkopplung $2K_d$ | frei: Schwebung $\bar\omega_2-\bar\omega_1=\sqrt{\Delta\omega^2-4K_d^2}$; eingerastet: $\bar\omega_2=\bar\omega_1$, Phasendicke $\operatorname{Var}\psi=\dfrac{D}{A^{*2}\,2K_d\cos\psi^*}$ |

Die Parameter $p$ greifen in alle drei Gleichgewichte ein, über $\mu(p)$, $\omega_j(p)$, $K(p)$ und $D(p)$. Das ist die Rechenform von „Multiphysik ist keine Summe von Kräften, sondern eine Organisation von Flussrichtungen“: Die Torusform ist ein Gleichgewicht zwischen Rollen, nicht zwischen Kräften gleicher Art.

**Dünn werden.** Auf dem freien Torus ist die invariante Dichte der Phasendifferenz $\rho(\psi)\propto1/(\Delta\omega-2K_d\sin\psi)$. Ihre Konzentration ist

$$
r=\lvert\langle e^{i\psi}\rangle\rvert=\frac{\Delta\omega-\sqrt{\Delta\omega^2-4K_d^2}}{2K_d}\;\xrightarrow{\ |\Delta\omega|\to2K_d\ }\;1 .
$$

Der Torus wird in Phasenrichtung dünn, bevor er einrastet, und die Schwebung verschwindet mit einem Wurzelgesetz. Das ist das „$T^2$ dünn machen“ des Autors in geschlossener Form.

**Einrasten** $T^2\to T^1$: $|\Delta\omega|\le2K_d$.

**Zerfallen** hat drei Wege:
1. *Amplitudentod durch Kopplung* ($T^2\to$ Fixpunkt). Die Linearisierung im Ursprung hat die Eigenwerte $\mu-K_d+i\bar\omega\pm\sqrt{K_d^2-\Delta\omega^2/4}$. Der Ursprung wird stabil, beide Zyklen erlöschen, genau dann wenn
$$K_d>\mu\quad\text{und}\quad|\Delta\omega|>2\sqrt{\mu\,(2K_d-\mu)} .$$
2. *Hopf-Tod eines Zyklus*: $\mu_j(p)\to0$. Das ist der Alterungsweg aus 03 §6.3.
3. *Zerfasern*: Die relative Dicke $h_r$ divergiert mit $(\mu-K_d)^{-1}$. Der Torus bleibt dem Typ nach erhalten, ist aber vom Rauschen nicht mehr aufgelöst.

**Satz 2 (Verriegelung braucht den dissipativen Kopplungsanteil).** Einrasten und Amplitudentod hängen nur von $K_d=K\cos\varphi$ ab. Eine rein reaktive Kopplung ($\varphi=\pi/2$) rastet bei gleichen Amplituden nie ein, wie groß $K$ auch ist. Das folgt auch aus dem Satz von Liouville: Ohne Dissipation rastet kein Umlauf ein.

### 2.3 Durchgerechnet

Parameter $\mu=\beta=1$, $\bar\omega=10$ (Teil 2a, RK4, $t=600$):

| Fall | $K$ | $\Delta\omega$ | $\varphi$ | Vorhersage | $A$ vorhergesagt | $A$ simuliert | Schwebung vorh. / sim. | $r$ vorh. / sim. |
|---|---|---|---|---|---|---|---|---|
| A | 0,1 | 1,00 | 0 | freier $T^2$ | 0,9487 | 0,9479 | 0,9798 / 0,9799 | 0,101 / 0,101 |
| B | 0,1 | 0,21 | 0 | freier $T^2$, dünn | 0,9487 | 0,9476 | 0,0640 / 0,0642 | 0,730 / 0,729 |
| C | 0,3 | 0,50 | 0 | eingerastet $T^1$ | 0,9305 | 0,9305 | 0 / 0 | 1 / 1,000 |
| D | 1,5 | 4,00 | 0 | Amplitudentod ($\Delta\omega_c=2{,}8284$) | 0 | 0,0000 | – | – |
| E | 1,5 | 2,00 | 0 | eingerastet trotz $K>\mu$ | 0,7862 | 0,7862 | 0 / 0 | 1 / 1,000 |
| F | 0,3 | 0,50 | $\pi/2$ | frei (rein reaktiv) | 1,0000 | 0,9688 | 0,5000 / 0,5213 | 0 / 0,005 |

Fall F rastet nicht ein, obwohl $|\Delta\omega|<2K$ gilt (Satz 2). Die Abweichung des Radius in F kommt von der Amplitudenmodulation durch $K_r$ (0,783 bis 1,103).

**Dicke** (Teil 2b, stochastisches Heun, 400 Realisierungen, $D=10^{-3}$):

| | Rauschanteil $\operatorname{Var}A$ | vorhergesagt | relative Dicke |
|---|---|---|---|
| A, frei | $2{,}810\cdot10^{-4}$ | $2{,}779\cdot10^{-4}$ | 0,0177 |
| C, eingerastet | $3{,}211\cdot10^{-4}$ | $2{,}892\cdot10^{-4}$ | 0,0193 |
| C, Phase $\operatorname{Var}\psi$ | $3{,}609\cdot10^{-3}$ | $3{,}487\cdot10^{-3}$ | |

Im Fall A ist die deterministische Schwebung ($1{,}052\cdot10^{-3}$) herausgerechnet. Im eingerasteten Fall liegt die Amplitudenvarianz 11 % über der Einmoden-Formel, weil bei $\sin\psi^*\ne0$ Phasenschwankungen über $K_dA_2\cos\psi$ in die Amplitude übersprechen.

### 2.4 Parametrischer Pfad: Einrasten, Dünnwerden, Zerfall

Die Temperatur als Parameter (Teil 2c): $\mu(\vartheta)=1-0{,}1\vartheta$ (Verluste steigen), $\Delta\omega(\vartheta)=0{,}6+0{,}05\vartheta$ (differentielle Frequenzdrift), $K_d=0{,}4$, $D(\vartheta)=10^{-3}(1+0{,}1\vartheta)$.

| $\vartheta$ | $\mu$ | $\Delta\omega$ | Regime | $A$ | Schwebung | $r$ | relative Dicke |
|---|---|---|---|---|---|---|---|
| 0 | 1,00 | 0,60 | eingerastet $T^1$ | 0,9298 | 0 | 1,000 | – |
| 3,8 | 0,62 | 0,79 | eingerastet, Radius schrumpft | 0,5320 | 0 | 1,000 | – |
| 4,2 | 0,58 | 0,81 | frei $T^2$, dünn | 0,4077 | 0,1360 | 0,843 | 0,1047 |
| 5,0 | 0,50 | 0,85 | frei $T^2$ | 0,2941 | 0,2921 | 0,697 | 0,1936 |
| 5,8 | 0,42 | 0,89 | frei, zerfasert | 0,1311 | 0,3957 | 0,618 | 0,9937 |
| 6,2 | 0,38 | 0,91 | Amplitudentod | 0,0000 | – | – | – |

Die Schwellen folgen analytisch: Entriegeln bei $\vartheta=4{,}0$ ($\Delta\omega=2K_d$), Zerfall bei $\vartheta=6{,}0$. An dieser Stelle ist $\Delta\omega=0{,}9>2\sqrt{\mu(2K_d-\mu)}=0{,}8$. Kurz vor dem Tod erreicht die relative Dicke 0,99: Der Torus zerfasert, bevor er erlischt.

**Zahlenbeispiel** (Annahmen, Teil 2d): zwei selbsterregte MEMS-Oszillatoren bei $f_0=1$ MHz, Nettoverstärkung $\mu_0=2\pi\cdot50$ Hz, dissipative Kopplung über die gemeinsame Einspannung $K_d=2\pi\cdot20$ Hz, eine $\vartheta$-Einheit entspricht 10 K. Der Rastbereich ist $|\Delta f|\le40$ Hz. Kalt ($\Delta f=30$ Hz) ist das Paar eingerastet. Bei einer differentiellen Drift von 0,25 Hz/K und einem Verstärkungsabfall von 0,5 Hz/K entriegelt es bei $\Delta T=40$ K (die Schwebung setzt mit dem Wurzelgesetz ein), und bei $\Delta T=60$ K erlöschen beide Zyklen.

---

## 3  Die sechs Übergänge des Dreiecks als Wechselkräfte

### 3.1 Die Übergänge als Kopplungsoperatoren

Die drei Ecken werden zu drei Koordinaten des Zustandsraums:

- $E$: Energie in den schnellen Zyklen (bzw. deren Wirkungen $J_i$),
- $x_s$: reversible Strukturkoordinate (Wärmedehnung, Vorspannung, Geometrie, Domänenordnung),
- $Z$: einseitige Zustandskoordinate, also Gedächtnis ($D$, $A$, $\varepsilon_p$ aus 03).

| Übergang | Tafeltext | Kopplungsoperator | Beispiel |
|---|---|---|---|
| Energie → Struktur | verformt, erwärmt, ordnet | $t_S\dot x_s\ni aE$ | Wärmedehnung $\alpha\,\Delta T$, Maxwell-Druck $B^2/2\mu_0$ |
| Struktur → Energie | lenkt, verstärkt, dämpft | $\dot E\ni2P\,\mathcal L(x_s)-2\gamma(x_s)E$, $\mathcal L(x_s)=\dfrac1{1+(\Delta_0-x_s)^2}$ | Resonanzlage, Verluste, Modenauswahl |
| Struktur → Zustand | altert unterschiedlich | $\dot Z\ni k_{A,S}\,x_s^2E$ | Spannungskonzentration, Kachanov $(\sigma/(1-D))^m$ |
| Zustand → Struktur | verändert die Form | $t_S\dot x_s\ni d\,Z$ | Kriechen, Setzen, bleibende Verformung |
| Energie → Zustand | hinterlässt Spuren | $\dot Z\ni k_AE$ | Erwärmung, Ermüdung (Miner), Arrhenius-Chemie |
| Zustand → Energie | verändert die Antwort | $\gamma(Z)=\gamma_0(1+bZ)$ | Widerstand, Permeabilität, Rauschen $D(Z)$ |

Hinzu kommt das Ausheilen $-hZ$. Ohne Ausheilen ($h=0$) ist $Z$ ein reiner Integrator.

### 3.2 Der geschlossene Kreis als 3×3-System

Linearisiert um einen Arbeitspunkt, mit dimensionslosen Kantenverstärkungen $g_{ij}$ (von $j$ nach $i$) und den Zeitskalen $t_E\llt_S\llt_Z$:

$$
\dot{\mathbf x}=J\mathbf x,\qquad J=\operatorname{diag}\!\Big(\frac1{t_E},\frac1{t_S},\frac1{t_Z}\Big)(G-I).
$$

Das charakteristische Polynom $\lambda^3+a_1\lambda^2+a_2\lambda+a_3$ hat die Koeffizienten

$$
a_1=\sum_i\frac1{t_i},\qquad a_2=\sum_{i<j}\frac{1-q_{ij}}{t_it_j},\qquad
\boxed{\;a_3=\frac{\det(I-G)}{t_Et_St_Z}=\frac{1-q_{ES}-q_{EZ}-q_{SZ}-q_{ESZ}-q_{EZS}}{t_Et_St_Z}\;}
$$

Dabei sind $q_{ij}=g_{ij}g_{ji}$ die drei 2-Zyklen und $q_{ESZ}$, $q_{EZS}$ die beiden Umlaufrichtungen des 3-Zyklus. Nach Hurwitz ist der Kreis genau dann stabil, wenn $a_1>0$, $a_3>0$ und $a_1a_2>a_3$ gilt. Teil 3a prüft das an 20 000 Zufallssystemen ($g\sim\mathcal N(0,1)$, $t_i$ zwischen 0,1 und 1000): Die Übereinstimmung mit den Eigenwerten beträgt 100,000 %. 14 295 Systeme sind stabil, 5233 instabil über $a_3\le0$ und 472 instabil über $a_1a_2<a_3$.

Daraus folgen **zwei Arten, sich aufzuschaukeln**:
- **Durchgehen** (statisch, reeller Eigenwert durch null): $\sum q\ge1$. Die Neumann-Reihe aus 03 §4 divergiert. Beispiele sind thermisches Durchgehen, tertiäres Kriechen (die Schleife Struktur ⇄ Zustand ist positiv) und Ermüdungsbruch.
- **Oszillatorisches Aufschaukeln** (Hopf): $a_1a_2<a_3$. Das braucht einen verzögerten negativen 3-Zyklus. Für den reinen 3-Zyklus $E\to S\to Z\to E$ mit den Vorzeichen $(+,+,-)$ liegt die Schwelle bei

$$
\boxed{\;|q_3|_{\rm krit}=\Big(\sum_it_i\Big)\Big(\sum_i\frac1{t_i}\Big)-1\;\ge\;8\;}
$$

(Teil 3b, numerisch = analytisch):

| $(t_E,t_S,t_Z)$ | $\lvert q_3\rvert_{\rm krit}$ |
|---|---|
| (1, 1, 1) | 8,0000 |
| (1, 3, 10) | 19,0667 |
| (1, 10, 100) | 122,2100 |
| (1, 100, $10^4$) | 10 202,0201 |

Gleiche Zeitskalen sind am gefährlichsten; das ist die Sekantenbedingung $\sec^3(\pi/3)=8$. **Zeitskalentrennung schützt quadratisch.** Das ist die formale Fassung von „Rückwirkung muss schneller sein als Aufschaukelung“ (`06` §3.5).

### 3.3 Selbststabilisierung, Alterung, Kippen

Im nichtlinearen Modell (Teil 3c) ist $\mathcal L(x_s)$ die Resonanzlage, die die Struktur verschiebt, mit $aP=5$, $\Delta_0=3$ und $t_S=10$. Ohne Zustand hat das System drei Fixpunkte, $x_s^*=1$, 1,382 und 3,618, mit den Schleifenverstärkungen 0,8, 1,236 und $-3{,}2361$. Der heiße Ast ($x_s^*=3{,}618$) ist durch eine **negative** Schleife stabilisiert, der Mittelast hat $q>1$. Bistabil ist das System für $2{,}9113<A<5{,}0887$ mit $A=aP/\gamma$; an den Rändern liegen die Falten mit $q=1$.

| Szenario | $b$ | $d$ | $k_A$ | $h$ | Ende $E$, $x_s$, $Z$ | Kippen |
|---|---|---|---|---|---|---|
| (i) Ausheilen | 0,5 | 0 | $10^{-3}$ | $5\cdot10^{-3}$ | 0,7117; 3,5584; 0,1423 | nie: **Attraktor** |
| (ii) Alterung | 0,5 | 0 | $10^{-3}$ | 0 | 0,0655; 0,3275; 1,7506 | $t=2289$ |
| (iii) Alterung ×4 | 0,5 | 0 | $4\cdot10^{-3}$ | 0 | 0,0483; 0,2414; 2,8138 | $t=616$ |
| (iv) $Z$ schiebt $x_s$ weg von der Resonanz | 0 | +0,5 | $10^{-3}$ | 0 | 0,4650; 4,0727; 3,4969 | nie: selbstbegrenzend |
| (v) $Z$ schiebt $x_s$ zur Falte | 0 | −0,5 | $10^{-3}$ | 0 | 0,0397; −1,9185; 4,2343 | $t=4818$ |

- **Ausheilschwelle.** Ein Attraktor auf dem heißen Ast existiert, solange $h>h_c\approx bk_AE^*/(\gamma_c-1)=5{,}04\cdot10^{-4}$ ist ($\gamma_c=aP/A_f=1{,}7174$). Heilung, die schneller ist als Schädigung, macht aus dem Kreis einen Attraktor.
- **Kippzeit.** Bei $Z_c=(\gamma_c-1)/b=1{,}4348$ verschwindet der heiße Ast. Die Abschätzung $Z_c/(k_A\langle E\rangle)=2197{,}5$ trifft den Faltendurchgang bei $t_F=2176$. Der Sprung folgt $113$ Zeiteinheiten später, das ist die Verzögerung an der Geisterbahn des Sattel-Knotens (03 §6.3). Viermal schnellere Alterung verkürzt die Kippzeit um den Faktor 3,716.
- **Vorzeichen gehört der Schließung.** (iv) und (v) sind dieselbe Kante $Z\to S$ mit entgegengesetztem Vorzeichen. Einmal begrenzt sich die Alterung selbst, einmal kippt das System. Das ist Befund 6 aus 03, jetzt für die Gedächtniskante.
- **Kritisches Verlangsamen.** Der langsamste Eigenwert des schnellen $(E,x_s)$-Blocks vor der Falte:

| $\gamma_c-\gamma$ | $t_F-t$ | $q$ | $\lambda_s$ | $\lambda_s/\sqrt{\gamma_c-\gamma}$ |
|---|---|---|---|---|
| 0,3 | 963 | −1,684 | −0,2874 | −0,5247 |
| 0,03 | 102 | 0,249 | −0,0745 | −0,4302 |
| 0,003 | 10 | 0,777 | −0,0218 | −0,3982 |
| 0,001 | 3 | 0,873 | −0,0124 | −0,3917 |

Es gilt $\lambda_s\propto\sqrt{\gamma_c-\gamma}\propto\sqrt{t_F-t}$. Das Kippen durch Alterung kündigt sich also an: Die Erholung nach kleinen Störungen wird mit einem Wurzelgesetz langsamer.

### 3.4 „Zustand ist Gedächtnis, wie heute noch wirkt“: der Gedächtniskern

Eliminiert man $Z$ aus dem linearen Kreis, so ergibt sich

$$
Z(t)=\int_{-\infty}^{t}e^{-h(t-t')}\big[k_AE(t')+k_{A,S}x_s^2E(t')\big]\,dt',\qquad
\dot E(t)\ni-2b\,E(t)\int_{-\infty}^{t}\underbrace{k_A\,e^{-h(t-t')}}_{\mathcal K_Z(t-t')}E(t')\,dt' .
$$

**Formalisierung des Satzes.** Der Zustand $Z$ ist das Faltungsintegral der Vergangenheit mit dem Kern $\mathcal K_Z(\Delta t)=k_Ae^{-h\Delta t}$. „Wie heute noch wirkt“ ist das Gewicht $\mathcal K_Z(t-t')$, mit dem die Energie zur Zeit $t'$ jetzt über die Kante $Z\to E$ die Antwort verändert. Mit $h=0$ wird nichts vergessen (Integrator, Miner-Regel), mit $h>0$ schwindet das Gedächtnis. Mehrere innere Variablen mit verschiedenen $h_k$ ergeben eine Prony-Reihe $\sum_kk_{A,k}e^{-h_k\Delta t}$ und im Grenzfall algebraische Kerne $\Delta t^{-\alpha}$ (Kriechnachgiebigkeit, $1/f$-Rauschen).

**Verbindung zum Toy Model.** Der Memory-Term des Autors,
$$
\ddot x=-\partial_x\Phi-\gamma y,\qquad \dot y=\dot x-y/t_K,
$$
ist exakt ein Teilchen an einem **Maxwell-Element**: einer Feder $\gamma$ in Reihe mit einem Dämpfer $d_D=\gamma t_K$. Der Kern $K(t)=e^{-t/t_K}$ ist dessen Relaxationsmodul, $G(t)=\gamma e^{-t/t_K}$. Die Auslenkung teilt sich in $x=y+z$:

- $y$ ist die Dehnung der Feder. Sie ist reversibel, speichert $\tfrac12\gamma y^2$ und gehört zur Ecke **Struktur**.
- $z$ ist der Fluss durch den Dämpfer, $\dot z=y/t_K$. Er ist einseitig dissipativ, $\dot H=-(\gamma/t_K)y^2$ (wie in `05` §5.2), und gehört zur Ecke **Zustand**.

Damit ist der Toy-Kern das Dreieck selbst. Energie (Bewegung) spannt die Struktur (Feder), die Struktur wirkt auf die Energie zurück ($-\gamma y$), und der Zustand (Dämpfer) entscheidet, wie schnell die Struktur vergisst. Bei harmonischer Bewegung (Teil 3d, $\gamma=0{,}1$) ist der Speicheranteil $\gamma(\omega t_K)^2/(1+(\omega t_K)^2)$ und bindet, der Verlustanteil $\gamma t_K/(1+(\omega t_K)^2)$ reibt:

| $\omegat_K$ | Speicher (Formel) | simuliert $(\omega^2-\omega_0^2)/\gamma$ | Abklingrate simuliert | Formel mit tatsächlichem $\omega$ |
|---|---|---|---|---|
| 0,1 | 0,00990 | 0,00992 | 0,00496 | 0,00495 |
| 1 | 0,50000 | 0,52559 | 0,02497 | 0,02437 |
| 10 | 0,99010 | 0,99179 | 0,00451 | 0,00451 |
| 100 | 0,99990 | 0,99992 | 0,00045 | 0,00045 |

**Befund.** Die Bindung ohne Energieverlust ($K\equiv1$) steckt im Toy Model in der **reversiblen Strukturkoordinate**, nicht im einseitigen Zustand. Der Zustand legt nur fest, wie lange die Bindung hält: Für $\omega t_K\gg1$ bindet der Kern, für $\omega t_K\ll1$ reibt er. Das präzisiert den Satz des Autors: Gedächtnis, das *bindet*, ist Struktur mit langer Relaxationszeit. Gedächtnis, das *altert*, ist Zustand.

---

## 4  Wärme als langsames Gedächtnis: Rückwirkung auf die EM-Regeln

### 4.1 Die Schleifenverstärkung der Kette

Die Kette der Tafel 3 lautet: Wärme → $\sigma(T),\mu(T),\varepsilon(T)$ → EM → Kraft → Verformung → Feld → Wärme. Mit der Heizleistung $P(T)$ und dem Wärmeleitwert $G$ gilt stationär $G\Delta T=P(\Delta T)$, und die Schleifenverstärkung ist

$$
q=\frac1G\frac{dP}{dT}=\underbrace{\frac1G\frac{\partial P}{\partial\sigma}\frac{d\sigma}{dT}}_{q_\sigma}+\underbrace{\frac1G\frac{\partial P}{\partial\omega_m}\frac{\partial\omega_m}{\partial E}\frac{dE}{dT}}_{q_E}+\underbrace{\frac1G\frac{\partial P}{\partial\omega_m}\frac{\partial\omega_m}{\partial N}\frac{dN}{dT}}_{q_N},\qquad
\delta T=\frac{\delta T_0}{1-q}.
$$

Die Temperatur tritt in $P$ **nur über Koeffizienten** auf: $\sigma$, den Elastizitätsmodul $E$ und die Wärmespannung $N$. Sie liefert keinen Kraftterm. Das macht die Rückkopplung **parametrisch**: Wärme ändert die Regeln, nach denen EM wirkt, und das ist Kanal B aus `06` (Medium → Phase) mit einer Schleifenverstärkung. Für $q\ge1$ divergiert $1/(1-q)$; es gibt dann mehrere stationäre Temperaturen, Sprünge und Hysterese.

### 4.2 Kupfer: eingespannter Draht im Magnetfeld

Das Beispiel ist ein magnetomotorischer Resonator: ein Kupferdraht, beidseitig eingespannt, im Feld $B$, mit einer Wechselspannung $V$ bei $f_d$ betrieben. Die Lorentzkraft $B\ell I$ treibt die Grundmode, die Bewegungs-EMK $B\ell\dot x$ wirkt als Impedanz $Z_m=(B\ell)^2/(c+i(m\omega-k/\omega))$ in den Kreis zurück. Geheizt wird durch $\tfrac12|I|^2(R+\operatorname{Re}Z_m)$.

Stoffwerte Kupfer [G]: $\rho_{\rm el}=1{,}68\cdot10^{-8}$ Ωm, $\alpha_R=3{,}93\cdot10^{-3}$/K, $E=128$ GPa, $(1/E)\,dE/dT\approx-3{,}8\cdot10^{-4}$/K, $\alpha=1{,}65\cdot10^{-5}$/K, $k=401$ W/mK. Annahmen: Radius 50 µm, Länge 20 mm, Vorspannung 20 MPa, $B=0{,}02$ T, $Q_m=1000$ (Vakuum). Galerkin-Mode $1-\cos(2\pi x/\ell)$.

| Größe | Wert |
|---|---|
| $R_0$, $G_{\rm th}=12kA/\ell$, $\tau_{\rm th}$ | 42,781 mΩ; 1,8897 mW/K; 0,2867 s |
| $f_{m0}$, $R_m/R_0$, belastete Güte $Q_l$ | 1610,73 Hz; 0,175; 851,0 |
| $\partial\ln\omega_m/\partial T$ über Wärmespannung / über $E(T)$ | $-3{,}7855\cdot10^{-2}$/K / $-5{,}3778\cdot10^{-5}$/K |
| Verschiebung pro Kelvin | 32,22 Linienbreiten |
| Knicken durch Wärmespannung | $\Delta T_b=13{,}208$ K |

Stationäre Lösungen und zerlegte Schleifenverstärkung (Auswahl):

| $V$ | $f_d-f_{m0}$ | $\Delta T$ | $q$ | $q_\sigma$ | $q_E$ | $q_N$ |
|---|---|---|---|---|---|---|
| 2 mV | −2 Hz | 0,0222 K | −0,151 | $-7{,}8\cdot10^{-5}$ | $-2{,}2\cdot10^{-4}$ | −0,1512 |
| 2 mV | 0 | 0,0236 K | +0,066 | $-8{,}8\cdot10^{-5}$ | $+9{,}3\cdot10^{-5}$ | +0,0655 |
| 5 mV | −10 Hz | 0,1450 K | −0,612 | $-5{,}3\cdot10^{-4}$ | $-8{,}7\cdot10^{-4}$ | −0,6104 |
| 5 mV | −7 Hz | 0,1511 K | +0,160 | $-5{,}8\cdot10^{-4}$ | $+2{,}3\cdot10^{-4}$ | +0,1600 |
| 10 mV | −40 Hz | 0,6060 K | −0,476 | $-2{,}3\cdot10^{-3}$ | $-6{,}7\cdot10^{-4}$ | −0,4726 |

Größte Schleifenverstärkung über alle Lösungen und Bistabilitätsfenster (Scan in 0,25-Hz-Schritten):

| $V$ | $\max q$ | Fenster mit drei stationären Temperaturen |
|---|---|---|
| 2 mV | +0,153 | keins |
| 5 mV | +0,603 | keins |
| 6 mV | +1,232 | bei −11,5 Hz |
| 7 mV | +1,898 | −16,25 … −15,75 Hz |
| 10 mV | +3,924 | −35,0 … −32,5 Hz |

Ein quasistatischer Frequenzsweep bei 10 mV zeigt die Hysterese im selben Fenster, mit einem Temperatursprung von 88,5 mK (525,9 gegen 614,3 mK bei −32,5 Hz).

**Befunde.**
1. Die Schleife schließt sich tatsächlich über den ganzen Weg des Autors. Heizen erzeugt Wärmespannung, die Wärmespannung verschiebt die Resonanz, die Resonanz bestimmt über die Bewegungs-EMK die elektrische Impedanz und damit Strom und Joule-Wärme.
2. **Die Regeländerung läuft überwiegend über die Struktur.** Der strukturelle Kanal übertrifft den Leitwertkanal $\sigma(T)$ um das Verhältnis $Q_l\,|\partial\ln\omega/\partial T|/\alpha_R=8197$. In der Tabelle ist $q_N/q_\sigma$ je nach Arbeitspunkt 200 bis 1900. $\sigma(T)$ allein bleibt bei $|q_\sigma|\le2{,}4\cdot10^{-3}$.
3. Das Vorzeichen wechselt mit der Seite der Resonanz. Nähert die Erwärmung die Resonanz der Anregung an, dann wächst der Bewegungswiderstand, Strom und Heizung sinken, und die Schleife gegenkoppelt ($q<0$, Selbstverriegelung). Auf der anderen Flanke koppelt sie mit. Die Bistabilität setzt zwischen 5 und 6 mV ein. Die Zeitskalen sind $\tau_{\rm th}=0{,}287$ s und $\tau_m=2Q_l/\omega=0{,}168$ s; sie liegen damit im gefährlichen Bereich von §3.2, und eine dynamische Rechnung (Selbstpulsieren) steht noch aus.

### 4.3 Ferrit: Selbstregelung am Curie-Punkt

Modell: Heizleistung $P=P_0(1-T/T_C)^{0{,}72}$ (Mean-Field-Skalierung $M_s^2$), $T_C=493$ K (MnZn typisch) [G], $T_0=293$ K. Das ist die Umkehrung des Kupferfalls: Die Regel, nach der EM heizt, erlischt am Curie-Punkt.

| $P_0/G$ | $T^*$ | $T_C-T^*$ | $q$ | $1/(1-q)$ | $dT^*/d\ln P_0$ | ohne Regelung |
|---|---|---|---|---|---|---|
| 100 K | 336,727 K | 156,27 K | −0,201 | 0,832 | 36,39 K | 100 K |
| 1000 K | 453,990 K | 39,01 K | −2,971 | 0,252 | 40,54 K | 1000 K |
| 10 000 K | 490,878 K | 2,12 K | −67,143 | 0,0147 | 2,90 K | 10 000 K |

Die Temperatur rastet am Curie-Punkt ein. Eine Verzehnfachung der Leistung von 1000 auf 10 000 K verschiebt $T^*$ nur um 36,9 K. Das ist Wärme als parametrischer **Gegenkoppler**: Der Arbeitspunkt ist der Punkt, an dem die EM-Regel ihre Wirkung verliert (Prinzip der Curie-Punkt-Lötspitze).

---

## 5  Gravitation als permanenter Hintergrund

### 5.1 Vorspannung → Steifigkeit → Resonanz

Eine konstante Kraft verschiebt in linearer Elastizität nur die Ruhelage. Die Steifigkeit ändert sie erst über die **geometrische Steifigkeit** der Vorspannung. Für einen senkrechten Kragbalken unter Eigengewicht (Galerkin-Mode $1-\cos(\pi x/2\ell)$) gilt

$$
k_{\rm el}=\frac{\pi^4}{32}\frac{EI}{\ell^3}=3{,}0440\,\frac{EI}{\ell^3},\qquad k_{\rm geo}=\Big(\frac{\pi^2}{16}-\frac14\Big)\rho Ag=0{,}3669\,\rho Ag,\qquad
\omega^2=\frac{k_{\rm el}\pm k_{\rm geo}}{m_{\rm eff}}
$$

mit + hängend (Zug) und − stehend (Druck). Kontrolle: Stehend knickt der Balken bei $\ell_c^3=8{,}298\,EI/(\rho Ag)$; exakt (Greenhill) ist $7{,}837$, die Galerkin-Schranke liegt 5,9 % darüber. Mit dem **Gravitationsanteil der Steifigkeit** $f_g=k_{\rm geo}/(k_{\rm el}+k_{\rm geo})$ lautet das Hebelgesetz

$$
\boxed{\;\frac{\partial\ln\omega}{\partial\ln g}=\tfrac12f_g\ \ (\text{hängend}),\qquad -\frac{f_g}{2(1-2f_g)}\ \ (\text{stehend})\;}
$$

In FTM-Größen gilt mit $g=3c^2\Delta s/r^2$ und dem Takt $\sqrt{1-6\Delta s/r}$:

$$
\boxed{\;\omega^{\rm obs}(r)=\sqrt{\frac{k_{\rm el}\pm\kappa_g\,\rho A\,3c^2\Delta s/r^2}{m_{\rm eff}}}\;\sqrt{1-\frac{6\Delta s}{r}}\;},\qquad \kappa_g=0{,}3669 .
$$

### 5.2 Labor gegen Stern und Neutronenstern

Teil 5, Stahl mit $\sigma_y=250$ MPa:

| Körper | $g$ [m/s²] | $\Delta s$ [m] | Takt | $h_{\max}=\sigma_y/\rho g$ | $\ell_c$ (Stab 10 mm) | $f_g$ Si-Balken 100×1 µm | $f_g$ Stahlstab 1 m, 10 mm | $f_g$ Stahldraht 1 m, 0,1 mm |
|---|---|---|---|---|---|---|---|---|
| Erde | 9,820 | $1{,}4782\cdot10^{-3}$ | 0,999999999304 | 3243 m | 5,03 m | $1{,}96\cdot10^{-7}$ | $7{,}38\cdot10^{-3}$ | 0,9867 |
| Sonne | 274,3 | 492,31 | 0,999997877 | 116 m | 1,66 m | $5{,}47\cdot10^{-6}$ | 0,172 | 0,9995 |
| Sirius B | $4{,}017\cdot10^6$ | 501,17 | 0,999740741 | 7,93 mm | 67,7 mm | 0,0742 | 0,9997 | 1,0000 |
| Neutronenstern | $1{,}291\cdot10^{12}$ | 689,23 | 0,809558 | 24,7 nm | 0,99 mm | 1,000 | 1,000 | 1,000 |

(Sonne und Sirius B zeigen hier nur die Skala von $g$, keine realen Bauteile.) Im Labor ist Gravitation für steife Mikrostrukturen ein Effekt von $10^{-7}$ und für schlanke Stäbe von $10^{-2}$. Der Frequenzunterschied zwischen hängend und stehend beträgt $1{,}958\cdot10^{-7}$ (Si-Balken) bzw. $7{,}460\cdot10^{-3}$ (Stahlstab) und ist messbar. Auf einem Neutronenstern gibt es makroskopische Strukturen gar nicht mehr: Stahl trägt sich nur bis 25 nm Höhe.

### 5.3 Der Rollenwechsel der Gravitation

- Für $f_g\ll1$ ist Gravitation **parametrisch**: Sie verschiebt Resonanzen um $\tfrac12f_g\,\delta\ln g$ und formt den Torus, ohne auf ihm zu laufen.
- Für $f_g\to1$ ist Gravitation selbst die Rückstellkraft (Pendel, hängender Draht, g-Moden von Sternen). Dann **erzeugt** sie eine Umlaufrichtung und ist tangential.
- Stehend und bei $f_g=\tfrac12$ verschwindet die Steifigkeit: Der Zyklus erlischt (Knicken). Das ist der Weg $T^2\to T^1$ durch Verlust (03 §5.4).

Welche Rolle Gravitation spielt, legt also $f_g$ fest und nicht die Gravitation selbst. Das ist der Satz des Autors „tangential meint nicht die Art der Kraft, sondern ihre Rolle“ mit einer Kennzahl. Der Satz „Gravitation ist keine Umlaufrichtung“ gilt im Bereich $f_g\ll1$, in dem technische Systeme fast immer liegen.

### 5.4 Differentieller Kanal gegen Takt-Kanal

Bei einer radialen Verschiebung ändert sich die Frequenz über zwei Kanäle: differentiell über $\tfrac12f_g\,d\ln g$ und gemeinsam über den Takt. Ihr Verhältnis ist in FTM-Form

$$
R=\frac{\tfrac12f_g\,|d\ln g/dr|}{d\ln\sqrt{1-6\Delta s/r}/dr}=f_g\,\frac{r-6\Delta s}{3\Delta s}.
$$

$R/f_g$ beträgt $1{,}4367\cdot10^9$ auf der Erde (wie in 03: Pendel um 1 m angehoben, $\delta\omega/\omega=-1{,}5696\cdot10^{-7}$ gegen $1{,}0925\cdot10^{-16}$ über den Takt), $4{,}7105\cdot10^5$ an der Sonne, $3{,}8556\cdot10^3$ an Sirius B und **3,8036 am Neutronenstern**. Erst am Neutronenstern werden beide Kanäle vergleichbar. Für $f_g=1$ sind sie gleich groß bei $r=9\Delta s=3GM/c^2$. Das ist zahlenmäßig der Radius der Photonensphäre. Die Übereinstimmung ist allerdings nicht robust: Mit der Eigenbeschleunigung $a=g/\sqrt{1-6\Delta s/r}$ des statischen Beobachters wird $R'=f_g[(r-6\Delta s)/3\Delta s+\tfrac12]$, und die Gleichheit rückt nach $r=7{,}5\Delta s$.

---

## 6  Neue Ergebnisse

Neu aus den Wechselkräften:

**N1 · Erhaltungssatz: das Kac-Budget.** Auf einem geschlossenen, ergodischen Zustandsraum gilt für jede Zerlegung in Zellen $A_i$: $\langle T_{A_i}\rangle=\mu(X)/\mu(A_i)$, also
$$
\boxed{\;\sum_i\frac1{\langle T_{A_i}\rangle}=1\;}
$$
Das Aufenthaltsbudget ist erhalten und endlich. Auf einem offenen Raum sind alle mittleren Rückkehrzeiten unendlich, die Summe ist 0. Rechnung (Teil 6c, 10 Zellen): Rückkehrzeiten 9,9999 bis 10,0001, $\sum1/\langle T_i\rangle=1{,}000000$. Das ist die erhaltene Größe hinter „globale Gebundenheit durch Schließung“: Schließung verteilt ein festes Budget an Wiederkehr, Öffnung setzt es auf null.

**N2 · Schließung ändert den Generator nicht** (Satz 1). $[\mathcal P,\mathcal L]=0$, und die Fluchtmoden sind $\operatorname{Bild}\Pi_\perp$. Die Bindung skaliert wie Energieskala/$L^2$ (§1.3).

**N3 · Verriegelungssatz** (Satz 2). Nur der dissipative Kopplungsanteil $K_d=K\cos\varphi$ rastet ein und tötet Amplituden; der reaktive Anteil verschiebt nur Frequenzen.

**N4 · Die drei Enden des Torus mit Schwellen:** dünn ($r\to1$, Schwebung mit Wurzelgesetz), eingerastet ($|\Delta\omega|\le2K_d$), zerfallen ($K_d>\mu$ und $|\Delta\omega|>2\sqrt{\mu(2K_d-\mu)}$; oder $\mu\to0$; oder $h_r\to1$).

**N5 · Auswahlregel: welche Zyklen gleichzeitig existieren.** Zwei Zyklen an einer gemeinsamen Verstärkung ($\dot I_i=2I_i(\alpha_i-\beta_iI_i-\theta_{ij}I_j)$, Lamb) existieren genau dann stabil nebeneinander, wenn
$$
\frac{\theta_{12}}{\beta_2}<\frac{\alpha_1}{\alpha_2}<\frac{\beta_1}{\theta_{21}}
$$
gilt. Daraus folgt $C=\theta_{12}\theta_{21}/\beta_1\beta_2<1$, aber nicht umgekehrt (Teil 6a):

| $\theta_{12},\theta_{21}$ | $C$ | Fenster erfüllt | Ergebnis |
|---|---|---|---|
| 0,5; 0,5 | 0,25 | ja | Koexistenz, $I=(0{,}70;0{,}60)$ |
| 0,9; 0,9 | 0,81 | ja | Koexistenz, $I=(0{,}763;0{,}263)$ |
| 1,5; 1,5 | 2,25 | nein | ein Zyklus verdrängt |
| 2,0; 0,3 | 0,60 | nein | ein Zyklus verdrängt, **obwohl** $C<1$ |

Parametrische Wirkungen verschieben $\alpha_i(p)$ und $\theta_{ij}(p)$ und schalten so Zyklen ein und aus. Das ist Klasse C „entscheidet, welche Phasenrichtungen existieren dürfen“ als Ungleichung.

**N6 · Stabilitätsgrenze für die Zahl freier Zyklen.** $N$ Zyklen mit Frequenzen in einem Band $W$, gekoppelt über ein gemeinsames Medium mit der Paarkopplung $k$ ($\dot\theta_i=\omega_i+k\sum_j\sin(\theta_j-\theta_i)$). Das erste Einrasten liegt bei (Teil 6b)

| $N$ | $k_1$ | $(N-1)k_1$ | $N k_s$ (vollständig synchron) |
|---|---|---|---|
| 4 | 0,1925 | 0,5775 | 0,800 |
| 8 | 0,0813 | 0,5687 | 0,720 |
| 16 | 0,0413 | 0,6187 | 0,680 |
| 32 | 0,0203 | 0,6297 | 0,660 |

mit dem Grenzwert von Kuramoto $2W/\pi=0{,}6366$. Also:
$$
\boxed{\;N_{\max}\simeq1+\frac{2W}{\pi k}\;}
$$
Mehr unabhängige Umläufe trägt ein gemeinsames Medium nicht. Die Grenze wächst mit der Bandbreite und fällt mit der Kopplung. Die Torusdimension ist damit eine beschränkte Ressource.

**N7 · Aufschaukeln des Dreiecks.** Statisch bei $\sum q\ge1$, oszillatorisch erst bei $|q_3|>(\sum t_i)(\sum1/t_i)-1\ge8$. Zeitskalentrennung schützt.

**N8 · Kippen durch Alterung ist vorhersagbar.** Ausheilschwelle $h_c=bk_AE^*/(\gamma_c-1)$, Kippzeit $\propto1/k_A$, Frühwarnung $\lambda_s\propto\sqrt{t_F-t}$, danach eine Geisterverzögerung von der Ordnung $10\,t_S$.

**N9 · Der Toy-Kern ist ein Maxwell-Element.** Bindendes Gedächtnis ist Struktur mit langer Relaxationszeit, alterndes Gedächtnis ist Zustand.

**N10 · Wärme ändert die EM-Regeln vor allem über die Struktur.** Im Kupferbeispiel ist der Wärmespannungskanal an den gerechneten Arbeitspunkten 200- bis 1900-mal wirksamer als $\sigma(T)$. Das Vorzeichen hängt von der Resonanzflanke ab.

**N11 · Die Rolle der Gravitation legt $f_g$ fest.** Differentieller Kanal und Takt-Kanal werden erst am Neutronenstern vergleichbar ($R/f_g=3{,}80$).

### Prüfbare Aussagen (Labor)

| Nr | Aufbau | Vorhersage |
|---|---|---|
| L1 | Kupferdraht wie in §4.2, Spannungsantrieb, Frequenzsweep | thermomechanische Bistabilität ab 5 bis 6 mV; bei 10 mV Hysteresefenster −35 … −32,5 Hz und Sprung von rund 90 mK. Gegenkopplung auf der einen, Mitkopplung auf der anderen Resonanzflanke. Die Zahlen hängen von den Annahmen $Q_m$, $B$ und Vorspannung ab; Skalierung über $q\propto V^2\,Q_l\,|\partial\ln\omega/\partial T|\,R_m/R$. |
| L2 | zwei selbsterregte Oszillatoren, Kopplungsphase $\varphi$ über einen Phasenschieber einstellbar | Rastbereich $|\Delta\omega|\le2K\cos\varphi$; bei $\varphi=\pi/2$ kein Einrasten. Vor dem Einrasten $r=(\Delta\omega-\sqrt{\Delta\omega^2-4K_d^2})/2K_d$ und Schwebung $\propto\sqrt{\Delta\omega^2-4K_d^2}$. Amplitudentod bei $K_d>\mu$ und $|\Delta\omega|>2\sqrt{\mu(2K_d-\mu)}$. |
| L3 | $N$ Resonatoren auf gemeinsamem Substrat (MEMS-Array) | erstes Einrasten bei $(N-1)k\approx0{,}58$ bis $0{,}63\,W$; mehr als $1+2W/\pi k$ freie Zyklen gibt es nicht. |
| L4 | thermisch bistabiler Resonator mit langsamer Alterung (Oxidation, Ermüdung) | Erholungsrate nach kleinen Stößen fällt wie $\sqrt{t_F-t}$ vor dem Sprung; mit Ausheilen oberhalb $h_c$ kein Sprung. |
| L5 | Stahlstab 1 m, Ø 10 mm, hängend gegen stehend | relative Frequenzdifferenz $7{,}46\cdot10^{-3}$; im Si-Mikrobalken $1{,}96\cdot10^{-7}$. |
| L6 | Ferritkern im Induktionsfeld (Kontrolle mit bekannter Technik) | Einrasten nahe $T_C$, $dT^*/d\ln P_0\to$ wenige K bei hoher Leistung. |

L2, L5 und L6 prüfen nur etablierte Physik, sie bestätigen die Rollenzerlegung. L1, L3 und L4 prüfen die Kombination der Aussagen dieser Notiz.

---

## 7  Selbstcheck

| Prüffrage | Antwort |
|---|---|
| Wird eine Wirkung als Objekt oder Substanz eingeführt? | Nein. Schließung ist ein Projektor, der mit dem Generator vertauscht, Bindung eine Eigenschaft des zulässigen Zustandsraums. Die Kopplungen sind Terme und Kanten, Gravitation ist ein Parameter mit der Kennzahl $f_g$. |
| Ist Gedächtnis ein Speicher? | Nein. Es ist ein Faltungskern, $\mathcal K_Z(\Delta t)=k_Ae^{-h\Delta t}$, bzw. der Relaxationsmodul eines Maxwell-Elements. Gespeichert wird nichts; die Vergangenheit steckt in der aktuellen Koordinate. |
| Was zeigt §1 für Galaxien? | §1 zeigt, wie Identifikation ohne Quelle bindet. In der Galaxie ist sie im Modell die Rückkehr, eine Translation im Zustandsraum mit der Resolvente $(I-qT)^{-1}$ (§1.4). [korr.] |
| Mystik oder literarische Zitate? | Keine. Verwendet werden Spektraltheorie, Poincaré, Kac, Pólya/Watson, Stuart-Landau, Adler, Hurwitz, Sekantenbedingung, Sattel-Knoten, Maxwell-Viskoelastizität, magnetomotorische Transduktion, Galerkin/Greenhill, Lamb/Lotka-Volterra und Kuramoto. |
| Stammen alle Zahlen aus dem Skript? | Ja, aus `rechnungen/torus_wechselkraefte.txt`. Stoffwerte für Kupfer sind Literaturwerte ($dE/dT$ gerundet); Draht, $Q_m$, $B$, MEMS- und Ferrit-Parameter sind als Annahmen gekennzeichnet. |

**Offen:**
- Satz 1 ist für Translationsgruppen bewiesen; die Galaxie ist im Modell ein solcher Fall, Schließung durch die Rückkehr als Translation (§1.4). [korr.]
- Der Kupferresonator ist quasistatisch gerechnet. Weil $\tau_m\approx\tau_{\rm th}$ ist, kann nach §3.2 ein thermomechanisches Selbstpulsieren auftreten; die dynamische Rechnung fehlt.
- Die Amplitudenvarianz im eingerasteten Fall liegt 11 % über der Einmoden-Formel. Die Kopplung zwischen Amplitude und Phase ist dort noch nicht geschlossen gerechnet.
- $N_{\max}$ gilt für gleichabständige Frequenzen und globale Kopplung. Zufällige Frequenzen und lokale Kopplung (Nachbarn auf dem Substrat) sind offen.
