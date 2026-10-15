# Systemarchitektur des Schauer-Faltungs-Fusionsreaktors (SFFR)

**Ausgearbeitet und nachgerechnet**

Francesco Silvano Schauer

---

Der Schauer-Faltungs-Fusionsreaktor (SFFR) verlässt die thermische Fusion, wie sie ITER und NIF verfolgen. Er heizt kein Plasma auf Hunderte Millionen Grad, sondern baut Kerne durch Neutroneneinfang auf: Das Neutron trägt keine Ladung und verschmilzt ohne Coulombwall mit dem Trägerkern. Im Zentrum eines geschlossenen Torus-Resonators stößt ein Injektor eine Kette von Neutronen an. Jedes Neutron kehrt mit dem Anteil q als Neutron der nächsten Generation zurück, und die Kette verstärkt die Quelle um

```
ν = 1/(1 − q) = 1 + q + q² + q³ + …
```

Das ist dieselbe Form wie die Rückkehr im Faltungsmodell (FTM), ν = 1/(1 − q) mit der Summe der Rückkehren Σqⁿ. Die Schwelle des Reaktors ist q = 1. Der SFFR bleibt darunter, bei q = 0,97: Die Kette lebt vom Injektor und endet mit ihm. Fusion und Spaltung stehen dabei im Gleichgewicht. Der Einfang baut aus ²³²Th das ²³³U auf und aus ⁶Li das ⁴He, das ²³³U spaltet und liefert die Neutronen der Kette.

Alle Zahlen kommen aus simulation/reaktor.py: Kernenergien aus AME2020, Halbwertszeiten aus NUBASE2020, Querschnitte bei 2200 m/s nach Sears (1992) und ENDF, Neutronen je Spaltung nach den IAEA Standards 2017, die Spallation nach der Messung am PSI (2003), der Wandler nach Backhaus und Swift (1999).

**Aufbau von innen nach außen**

- Zentrum: der Phasen-Injektor, ein Protonenstrahl auf ein Ziel aus Blei und Bismut
- Kavität: ein Ring aus Graphit, durchströmt von flüssigem Brennstoff mit ²³²Th und ²³³U
- Innenwand: ein Reflektor aus Graphit
- Mantel: ⁶Li im Vakuum der bifilaren Torusspule
- außen: die Resonanz-Abnehmer

## Modul A: Der topologische Schwellen-Resonator (bifilare T²-Faltungsgeometrie)

**1. Geometrische Anordnung der Wicklung**

- **Bifilare Wicklung.** Eine einfache Torusspule schreitet mit ihren Windungen einmal um den Ring fort und wirkt außen wie eine einzige Schleife mit dem Moment I·πR². Bei I = 5 kA und R = 1 m gibt sie in 10 m Abstand noch 3,1·10⁻⁶ T. Die zweite Lage schreitet gegensinnig fort und hebt diese Schleife auf. Dann ist das Feld innen B = μ₀NI/(2πr), 1,0 T bei N = 1000 Windungen und r = 1 m, und außen null. Die Kapazität zwischen den Lagen bildet mit der Induktivität einen Schwingkreis: Bei seiner Eigenfrequenz heben sich kapazitive und induktive Phase genau auf. Die Spule hat keinen Eisenkern und damit keine Hystereseverluste.
- **Kac-Wickelschema.** Die Wicklung erhält eine Verdrillung mit dem Umlaufverhältnis ι = (√5 − 1)/2, dem goldenen Schnitt. Eine Feldlinie schließt dann nie genau, sie kehrt aber jeder Stelle beliebig nahe zurück, zu einem Bereich vom Anteil μ der Fläche im Mittel nach 1/μ Umläufen (Kac). Die Flächen mit dem goldenen Verhältnis halten Störungen am längsten stand (KAM). Geladene Teilchen bleiben auf ihren Flächen, der Verlust quer zum Feld bleibt klein.

**2. Der Schwellenradius r* und das Einrasten**

Ein Neutron entsteht mit rund 2 MeV und rastet im Graphit ein, wenn es auf die Wärmebewegung des Graphits gebremst ist (0,0253 eV). Dazu braucht es 115 Stöße und 143 µs. Danach wandert es im Gleichgewicht mit dem Graphit, bis es eingefangen wird. Der Radius, den es dabei überstreicht, ist die Diffusionslänge, im Graphit der Dichte 1,7 g/cm³ L = 54 cm. Das ist der Schwellenradius der Kavität: Der Ring aus Graphit ist einige L dick, damit jedes Neutron innerhalb seines Radius eingefangen wird. Im Graphit allein lebt ein eingerastetes Neutron 15 ms.

Die Schwelle selbst ist q = 1. Unterhalb trägt der Injektor die Kette, oberhalb trüge sie sich selbst und wüchse. Der Abstand 1 − q = 0,03 ist die Reserve des SFFR.

## Modul B: Der Wardenclyffe-Phasen-Injektor (E0-Kopplung und Stehwelle)

**1. Die geschlossene Reaktorkavität als T²-Torus**

- **Die Innenwand.** Ein dicker Reflektor aus Graphit wirft von 100 Neutronen, die ihn treffen, 94 zurück, die Albedo ist (1 − 2D/L)/(1 + 2D/L) = 0,94 mit der Diffusionskonstante D. Die Wand ist der Spiegel der Kavität.
- **Der zentrale E0-Injektor.** Ein Linearbeschleuniger aus supraleitenden Kavitäten führt Protonen auf 0,59 GeV und lenkt sie auf das Ziel aus Blei und Bismut im Zentrum. Jedes Proton zerschlägt dort Kerne und setzt 9,6 Neutronen frei (gemessen am PSI). Der Injektor strahlt nicht ab: Seine stehende Welle gibt ihre Energie an den Strahl.
- **Stehwellen-Ausbildung.** Im Kern bildet der Neutronenfluss selbst eine stehende Welle. Er erfüllt die Helmholtz-Gleichung ∇²φ + B²φ = 0 wie ein Feld im Hohlraum, die Grundmode reicht von Wand zu Wand und hat ihr Maximum am Injektor.

**2. Eliminierung des Skin-Effekts und der Kernerwärmung**

In den Kavitäten des Injektors schwingt die Welle als Stehwelle in supraleitendem Niob. Der Oberflächenwiderstand ist so klein, dass fast die ganze eingespeiste Leistung in den Strahl geht. Im Zeitmittel sind elektrische und magnetische Feldenergie gleich, E_el = E_mag, wie in jedem Resonator. Die Torusspule ohne Eisenkern hat keine Hysterese.

## Modul C: Phasen-getaktetes Injektionssystem (sanfter Einfang)

**1. Das Wirkprinzip des Neutroneneinfangs**

- **Keine Coulomb-Sperre.** Das Neutron trägt keine Windung der Ladung und erreicht den Kern bei jeder Energie.
- **Faltungstiefe und Halterlänge.** Der gehaltene Trichter des ²³³U hat den Radius λ_H·A^{1/3} mit λ_H = 1,25 fm und die Fläche 1,86 b. Eingefangen wird ein thermisches Neutron aber mit 576 b, dem 310-fachen. Seine Welle ħ/(m_n·v) ist 3700-mal größer als der Kern und überdeckt den Trichter ganz. Je langsamer das Neutron, desto größer die Fläche, mit der es einrastet: Der Einfang ist sanft, weil das Neutron langsam ist. Beim ⁶Li ist der Querschnitt 940 b, das 5800-fache seiner Kernfläche.
- **Kollektive Einrastung im Paar.** Im ²³³U schließt das Neutron ein Paar und macht 6,85 MeV frei. Im ²³²Th öffnet es ein Paar und macht nur 4,79 MeV frei (beides AME2020). ²³³U spaltet, ²³²Th bleibt und baut sich über zwei β-Zerfälle zum ²³³U auf: ²³³Th lebt 21,8 min, ²³³Pa 27,0 Tage (NUBASE2020). So baut der Einfang den Brennstoff, den die Spaltung verbraucht.
- **Gleichgewicht.** Von jedem im ²³³U absorbierten Neutron kommen η = ν·σ_f/(σ_f + σ_γ) = 2,29 zurück: eines für die nächste Spaltung, eines für den Aufbau des nächsten ²³³U aus ²³²Th, 0,29 für ⁶Li, Wand und die Reserve unter q = 1. Bei q = 0,97 nimmt das ²³³U 42,3 % der Neutronen auf, das ²³²Th ebenso viel, 15,3 % bleiben für den Rest.

**2. Mechanik der Injektionstaktung**

- Der Injektor liefert die Protonen in Paketen im Takt seiner Hochfrequenz. Die Kette glättet sie: Eine Generation dauert 0,90 ms (Bremsung und Leben bis zur Absorption, wenn der Graphit 5 % der Neutronen aufnimmt), der Kern sieht eine gleichmäßige Quelle.
- **Bindungsenergie direkt ins Feld.** Der Einfang im ⁶Li gibt 4,78 MeV ganz an geladene Teilchen, α mit 2,05 MeV und Triton mit 2,73 MeV. Im Vakuum der Torusspule laufen sie gegen einen Sammler auf Hochspannung und geben ihre Energie unmittelbar als Feldenergie ½CU² ab, ohne Umweg über Wärme. Ihre Bahnradien bei 1 T sind 0,21 m und 0,41 m, größer als der Spalt, sie erreichen den Sammler. Die langsamen Sekundärelektronen hält das Feld fest. Das sind 117 MeV je Proton.

## Modul D: Der passive FTM-Regelkreis und direkte Energieauskopplung

**1. Thermodynamische Selbstlimitierung**

Steigt die Leistung, dehnt sich der flüssige Brennstoff aus und tritt zum Teil aus dem Kern. Dann nehmen ²³³U und ²³²Th weniger Neutronen auf, q fällt und mit ihm die Verstärkung:

| Brennstoffdichte | q | ν = 1/(1 − q) |
|---|---|---|
| Betrieb | 0,970 | 33,3 |
| 1 % weniger | 0,968 | 31,7 |
| 10 % weniger | 0,954 | 21,6 |

Der Regelkreis des Autors gilt damit wörtlich: Überleistung, Ausdehnung, Rückkehranteil q fällt, die Einrastung bricht ab, die Reaktion geht zurück. Da q < 1 bleibt, gibt es kein Durchgehen. Wird der Strahl abgeschaltet, klingt die Kette mit l/(1 − q) = 30 ms ab, l die Generationszeit. Die verzögerten Neutronen der Spaltprodukte tragen danach nur einen kleinen Rest, der in Minuten vergeht. Die Nachwärme führt ein passiver Kreislauf ab.

**2. Direkte elektromagnetische Phasen-Auskopplung (Verzicht auf Turbinen)**

- **Resonanz-Abnehmer.** Die Wärme des Kerns treibt thermoakustische Wandler: In einem Resonator mit Helium erzeugt der Temperaturunterschied eine Schallwelle, ohne bewegte Teile, mit dem Wirkungsgrad 0,30 von Wärme zu Schall (Backhaus und Swift 1999). Wo der Druck am stärksten schwingt, sitzen lineare Generatoren, die Wechselstrom in der Phase der Welle abgeben. Aus 100 MW Wärme werden 30 MW Schall, den die Generatoren in Strom wandeln.
- **Direkte Wandlung.** Die geladenen Produkte aus dem ⁶Li geben ihre Energie unmittelbar als Spannung ab (Modul C).

## Steuerungs-Matrix für den Reaktorbetrieb

| Parameter | Betriebswert | Funktion im SFFR |
|---|---|---|
| Rückkehranteil q | 0,97 | Schwelle q < 1, die Kette lebt vom Injektor |
| Verstärkung ν = 1/(1 − q) | 33,3 | Neutronen je Quellneutron |
| Neutronen je Proton | 9,6 bei 0,59 GeV | Anstoß der Kette |
| Spaltungen je Proton | 125 | Wärme 22,5 GeV je Proton, das 38-fache der Strahlenergie |
| Neutronen je Absorption im ²³³U | η = 2,29 | Spaltung und Aufbau im Gleichgewicht |
| Einfangfläche | 310-fache Kernfläche | sanfter Einfang durch die Welle |
| Einrasten im Graphit | 115 Stöße, 143 µs | thermisches Gleichgewicht |
| Schwellenradius | L = 54 cm | Reichweite des eingerasteten Neutrons |
| Innenwand | Albedo 0,94 | Rückkehr der Neutronen |
| Feldverhältnis | E_el = E_mag | supraleitende Kavitäten |
| Ausdehnung um 10 % | q = 0,954, ν = 21,6 | Selbstbegrenzung |
| Strahl aus | Abklingzeit 30 ms | Abschaltung |
| Strahl für 100 MW Wärme | 4,45 mA, 2,6 MW | Antrieb |

## Fazit der Auslegung

Der SFFR baut Kerne durch sanften Neutroneneinfang auf, ohne Plasma und ohne Zündung. Der Torus-Resonator hält die Neutronen als stehende Welle im Graphit und wirft sie an seiner Wand zurück, der Injektor stößt die Kette an, Fusion durch Einfang und Spaltung stehen im Gleichgewicht. Jedes Proton von 0,59 GeV setzt 22,5 GeV Wärme frei. Die Selbstbegrenzung folgt aus q < 1 und aus der Ausdehnung des Brennstoffs, die Energie verlässt den Reaktor über stehende Schallwellen und geladene Teilchen, ohne Turbine. Die Größenordnung entspricht der Anlage MYRRHA, die in Belgien für 600 MeV, 4 mA und 100 MW Wärme entsteht.

## 1. Das Start- und Initialisierungsprotokoll

- **Ruhezustand.** Der Kern ist leer, q ist nahe null.
- **Frequenz-Sweep.** Jede supraleitende Kavität wird mit einer Frequenzrampe auf die gemeinsame Resonanz gestimmt, dann läuft der Strahl mit kleinem Strom.
- **Annäherung an die Schwelle.** Der Brennstoff wird stufenweise eingefüllt. Nach jeder Stufe wird ν als Verhältnis von Neutronen zu Strahlstrom gemessen und q = 1 − 1/ν bestimmt. Der Betrieb beginnt bei q = 0,97, die Schwelle q = 1 wird nie erreicht.
- **Kein externer Anlasser nötig.** Im Betrieb trägt der Reaktor seinen Injektor selbst. Gerechnet mit 30 % von Wärme zu Strom braucht er dafür einen Beschleuniger mit mindestens 8,8 % Wirkungsgrad. Bei einem Drittel nimmt der Beschleuniger 7,9 MW der 30 MW.

## 2. Die Echtzeit-Phasenschnittstelle

- **Kuramoto-Regelung.** Eine FPGA-Regelung misst Phase und Amplitude jeder Kavität und korrigiert sie über GaN- und SiC-Verstärker. Die Kavitäten sind gekoppelte Oszillatoren, ihre Ordnung r = |⟨e^{iθ}⟩| ist bei gaußisch gestreuter Phase e^{−σ²/2}, und dieser Anteil des Energiegewinns erreicht den Strahl. Bei 1° Streuung sind es 0,99985, bei 10° noch 0,985.
- **Phasen-Kopplung statt Amplituden-Gewalt.** Geregelt wird die Phase jeder Kavität gegen den Strahl, θ₁ − θ₂ = 0. Die Leistung des Reaktors stellt allein der Strahlstrom.

## 3. Messtechnik und In-situ-Diagnostik

Der Zustand im Kern wird ohne Sonde von außen gelesen:

- **Verstärkung.** ν ist das Verhältnis der Neutronen an Detektoren außerhalb der Wand zum Strahlstrom, daraus q = 1 − 1/ν.
- **Strahlpause.** Eine Pause von einigen Millisekunden zeigt das Abklingen der Kette mit α = (1 − q)/l = 33 je Sekunde. Aus α folgt q, unabhängig von der Eichung der Detektoren. MYRRHA unterbricht seinen Strahl dafür einmal in der Sekunde.
- **Impedanz.** Der Strahl belastet die Kavitäten und verschiebt Phase und Amplitude ihrer Welle. Daraus liest die Regelung den Strahlstrom, ohne den Strahl zu berühren.

## 4. Isotopen-Management und Transmutation

- **Modus A, Energie.** Fusion durch Einfang und Spaltung im Gleichgewicht: ²³²Th + n baut ²³³U auf, ²³³U spaltet, ⁶Li + n gibt ⁴He und Tritium. Je Gigawattjahr Wärme werden 423 kg ²³³U gespalten und 459 kg aufgenommen, ebenso viel baut der Einfang aus ²³²Th nach. Verbraucht wird Thorium.
- **Modus B, Transmutation.** Langlebige Aktiniden aus abgebrannten Brennelementen werden in den Kern gegeben und in der Kette gespalten. Der Einfang baut dabei gezielt Isotope auf, ²³³U aus Thorium und Tritium aus ⁶Li.

## 5. Die Netz-Ankopplung

Ein Stromnetz ist ein Netz gekoppelter Oszillatoren: Generatoren und Umrichter laufen synchron wie im Modell von Kuramoto. Der SFFR speist seinen Wechselstrom über einen netzbildenden Umrichter ein, der sich auf die Phase des Netzes einrastet. Über weite Strecken trägt Gleichstrom hoher Spannung die Leistung, so verbinden sich die Netze der Kontinente.
