# 09 Elektromagnetismus und Frequenzleiter

## Rolle im Zustandsraum

Die Lorentzkraft f = ρE + J×B wirkt **tangential**: An bewegter Ladung leistet sie keine Arbeit ((v×B)·v = 0), ihr Block ist antisymmetrisch und erhält Energie und Entropie. Sie erzeugt Umläufe und Phasen; auf Dipole wirkt das Drehmoment N = μ×B. Joule-Wärme und Hysterese sind **transversal**, für sie gilt dS/dt ≥ 0.

**Das Magnetfeld teilt den Transport:** Die Transportmatrix zerfällt nach Onsager-Casimir in einen symmetrischen Teil L_s (gerade im Feld: Leitfähigkeit, Seebeck, Peltier, Wärmeleitung; erzeugt Entropie, transversal) und einen antisymmetrischen Teil L_a (ungerade im Feld: Hall, Nernst, Righi-Leduc; entropieneutral, tangential). Das Feld legt also fest, welcher Anteil umläuft und welcher dissipiert. Am selben Arbeitspunkt kann es einen Torus erzeugen und wieder einrasten: Es wirkt als parametrischer Schalter.

**Thermische Rückwirkung:** Wärme ändert die elektromagnetischen Regeln über σ_el(T), μ(T), ε(T) und vor allem über Wärmespannung in der Struktur. Die Schleife Feld → Strom → Wärme → Widerstand → Feld hat die Schleifenverstärkung q = (∂P/∂T)/G. Ihr Vorzeichen hängt am Antrieb und an der Richtung von R(T): Bei steigendem Widerstand ist die stromgetriebene Schleife mitkoppelnd (bei q ≥ 1 thermisches Durchgehen, Bistabilität mit Hysterese) und die selbsterregte gegenkoppelnd; bei fallendem Widerstand (etwa VO₂ am Metall-Isolator-Übergang) ist es umgekehrt. In mechanisch eingespannten Leitern wirkt die Wärme über die Wärmespannung in der Struktur 200- bis 1900-mal stärker als über σ_el(T); das Vorzeichen wechselt dann mit der Flanke der Resonanz.

Die elektrische Wechselwirkung zweier Protonen ist α/α_G = 1,23555·10³⁶-mal stärker als ihre gravitative.

## Verwindung und magnetisches Moment

Ein rotierender Trichter ist verwunden, a_J = J/(Mc) ist diese Verwindung als Länge. Für das Elektron (J = ħ/2) ist a_J = λ̄_e/2 = r_e/(2α) = 193,08 fm. Faltung und Windung sitzen am selben Knoten; die Verwindung ist ein imaginärer Versatz i·a_J dieses Orts. Beide Multipolreihen sind M(ia_J)^l und Q_el(ia_J)^l: J = Mc·a_J und μ = Q_el·c·a_J, also g = 2 für jede Masse, Ladung und Drehung (eine starr drehende Ladung im Raum gäbe g = 1); beim Elektron ist das das Bohrsche Magneton μ_B = eħ/(2m_e) = 9,274·10⁻²⁴ J/T. Die Abweichung g − 2 folgt daraus nicht. Für jeden Spin ½ gilt a_J·Δs = ℓ_P²/6: Die Spinzahl a_J/(3Δs) = (J/ħ)/α_G(M) ist für Teilchen riesig (Elektron 2,9·10⁴⁴, Proton 8,5·10³⁷): Ihre Verwindung übertrifft ihre Faltung bei weitem. Für Spin ½ ist j = ½·(m_Pl/M)²: Extremal (j = 1) wird ein solcher Trichter bei M = m_Pl/√2, nahe der Planck-Masse; alle leichteren sind überextremal.

## Frequenzleiter

| Ebene (Blickrichtung) | Frequenz | Wert |
|---|---|---|
| Compton-Takt des Protons (Eigenzeit am Ort) | ω₁ = c/λ̄_p | 1,425·10²⁴ rad/s |
| Elektron im Wasserstoff | ω₂ = αc/a_B | 4,134·10¹⁶ rad/s (45,6 nm) |
| Molekül- und Gitterschwingung | ω₃ = v_s/r | ≈ 10¹³ rad/s |
| Faltungsrate des Protons (Konstantenskala) | Γ_G = α_G·ω₁ = Gm_p³c/ħ² | 8,42·10⁻¹⁵ s⁻¹ |
| Faltungstakt des Protons (Faltung im Zentrum, keine Raumlänge) | c/(3Δs_p) | 2,41·10⁶² s⁻¹ |
| Neutronenstern (1,4 M☉, 12 km), Umlauf an der Oberfläche (ferne Uhr; am Ort Takt 0,81) | √(GM/R³) | 1,04·10⁴ rad/s |
| Erdbahn | 2π/Jahr | 1,99·10⁻⁷ rad/s |
| Milchstraße (6·10¹⁰ M☉), Umlauf bei r* (v mit der Regel) | v/r* | 8,0·10⁻¹⁶ rad/s (248 Mio. Jahre) |

Carters Beziehung ist eine Aussage über diese Leiter: Γ_G·ω₁³ = 𝒞·(α·ω₂)⁴.

**Drei Raten jeder Masse:** Faltungstakt c/(3Δs), Compton-Takt ω_C = c/λ̄ und Faltungsrate Γ_G = α_G(M)·ω_C. Der Compton-Takt ist ihr geometrisches Mittel, ω_C² = Γ_G·c/(3Δs). Unter der Planck-Masse ist der Faltungstakt am schnellsten und die Faltungsrate am langsamsten, darüber umgekehrt, bei m_Pl sind alle drei gleich.

## Einrasten und Grenzen der Kohärenz

Zwei Zyklen mit Frequenzunterschied Δω und Kopplung K folgen dΔθ/dt = Δω − K sin Δθ. Sie rasten ein, wenn K ≥ |Δω| (dann Δθ = arcsin(Δω/K)); sonst laufen sie quasi-periodisch auf dem Torus weiter.

In warmer, feuchter Umgebung ist Kohärenz eng begrenzt: Die thermische Dephasierungszeit bei Körpertemperatur ist ħ/kT = 2,46·10⁻¹⁴ s, die Debye-Länge in physiologischer Salzlösung 0,78 nm. Quantenkohärenz zerfällt dort in Femto- bis Pikosekunden, elektrische Felder reichen nur über Nanometer.

**Phasenseite:** Die Windung hat keinen Verlust und zählt in beiden Sinnen; ihre Rücklaufsumme ist der Kamm ganzer Windungen, ohne Skala. Die Phase ändert den Takt nicht: Das Argument der Regel ist das Taktgefälle g_N, nie qE/m. Coulomb bleibt bei jedem Feld 1/r², auch wo qE/m unter a₀ liegt (1 kg mit 1 pC bei 110 V/m); einen elektromagnetischen Schwelleneffekt gibt es nicht. Auf dem geschlossenen Raum folgen dF = 0 (die Phase ist eine Windung) und ΣQ_el = 0; die Quellengleichung selbst braucht eine Wirkung für F.
