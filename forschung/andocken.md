# Andocken

Rechnung, noch nicht im Modell.

**Kopplungsfenster an Davis u. a. 2011.** An 68 Kinasehemmern gegen 442 Kinasen (Davis u. a. 2011) binden die stärkeren mehr: Die stärkste Bindung jedes Hemmers und die Zahl seiner Bindungen unter 10 µM haben die Rangkorrelation 0,43, mit 1 µM als Schwelle 0,41. Das trennt aber nicht von einer bloßen Verschiebung des ganzen Profils. Die Verstimmung lässt sich an der Sequenz ablesen: Je ähnlicher eine Kinase der stärksten Partnerin eines Hemmers ist (Gruppen aus drei Aminosäuren, ohne nahezu gleiche Mutanten), desto eher bindet er sie, im ähnlichsten Viertel 38 %, im fernsten 24 %, mit Gruppen aus zwei Aminosäuren 35 % und 25 %. Eine scharfe Grenze des Fensters zeigt dieses grobe Maß nicht. Die Form der Bindungen eines Hemmers, mit gleichmäßig verteilter Verstimmung und der Nachweisgrenze als Boden ohne freie Zahl gerechnet, trägt nicht: Die Daten haben weniger starke Nebenbindungen, ihr oberes Zehntel beginnt bei 0,69 der stärksten statt bei 0,85.

**Wirkstoff als gehaltener Trichter.** Ein Wirkstoff aus N schweren Atomen, gehalten auf dem Abstand seiner Bindungen, ist ein gehaltener Trichter wie Kern und Stern, R = λ_H·N^{1/3}. Er bindet nur über seine Hülle, die mit N^{2/3} wächst. Das ist der Kern im Spiegel: Der Kern bindet im Innern, ∝ A, und verliert an der Hülle, ∝ A^{2/3}. Der Wirkstoff bindet an der Hülle, sein Inneres berührt nichts. Die obere Grenze der gelesenen Bindung ist darum

```
pK = a·N^{2/3} − c
```

bis zum Treffpunkt, an dem seine Hülle die der Tasche erreicht: Die Tasche ist gefüllt, der Kontakt wächst nicht mehr. c trägt den Nullpunkt der Lesung, die Bezugskonzentration von 1 mol/l, eine Wahl des Beobachters. Die Bindung je Atom, (a·N^{2/3} − c)/N, ist am größten bei N_e = (3c/a)^{3/2}, wo die Bindung über die Hülle dreimal der Nullpunkt ist. Wo ein Wirkstoff je Atom am besten bindet, hängt so an der gewählten Bezugskonzentration und ist keine Eigenschaft des Moleküls.

Mit der Hülle wächst die Kopplung K und mit ihr das Kopplungsfenster: Größere Wirkstoffe rasten bei mehr Partnern ein. An Davis u. a. 2011 (schwere Atome aus den SMILES gezählt und von PubChem gelesen, gleich bei 68 von 68) haben Größe und Zahl der Partner unter 10 µM die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37, Größe und stärkste Bindung 0,25. Die Zahl der Partner wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8}, steiler als die Hülle: Gezählt wird nur, was über der Schwelle des Nachweises bindet, und wächst das Fenster, wächst auch der Teil seiner Breite über der Schwelle.

An ChEMBL (alle K_i und K_d mit pChEMBL, 235 599 Moleküle, N aus der Ligandeneffizienz, an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1) springt die stärkste Bindung, das oberste Prozent je Größe, zwischen 12 und 13 Atomen von 9,5 auf 10,4 und bleibt von dort bis 66 Atome zwischen 9,9 und 10,9. Ob dort die Tasche gefüllt ist oder die Messung nicht tiefer liest, trennen diese Daten nicht. Unter 13 Atomen sind die Moleküle zu wenige, um den Anstieg mit N^{2/3} zu prüfen.

**Spiegelpaare an ChEMBL.** An den K_i und K_d von ChEMBL (5985 Spiegelpaare, deren InChI bis auf die Schicht /m0 gegen /m1 gleich ist, über vertauschtes @ und @@ in den SMILES zu 95 % bestätigt, 10 246 Paare mit gemeinsamem Ziel) unterscheiden sich Spiegelbilder am selben Ziel im Median um 0,51 in pK, bei jeder Größe ähnlich, die Rangkorrelation mit der Zahl der Atome ist −0,07. Die Hülle zählt die Kontakte um das Zentrum nicht: Ohne Unterschied (unter 0,1) sind bis 15 Atome 14 % der Paare, bei 31 bis 40 Atomen 22 %.

**Leiter der Takte an den Daten.** An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare aus Ziel und Wirkstoff eine Wirkweise im Guide to Pharmacology (Fassung 2026.3.1), alle Antagonisten. BindingDB (Teilmenge aus Artikeln, Stand 10.2026) führt nur 22 k_off mit Ziel und Wirkstoff, keines trifft eine Wirkweise. Die Reihe der Wirkweisen nach k_off lässt sich an beiden nicht prüfen.

Gemessen sind hier Zahl und Größe der Partner an Davis u. a. 2011, die stärkste Bindung je Größe, die Spiegelpaare und die Kinetik an ChEMBL.
