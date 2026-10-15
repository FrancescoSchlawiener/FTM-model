# Weitere Formeln Ramanujans

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Zahlennähen ohne Herleitung gekennzeichnet („dieselbe Zahl −1/12 wie Grundverschiebung und Nullpunktenergie“; „Kaleidoskop-Winkel sind die Signaturen von Ramanujans alternativen Basen“, zwei Zahlenmengen gleichgesetzt); historische Angabe zu den Signaturen mit [G]; „scheinbare Masse“ → gelesene Masse; Lambert-Reihe L(s) → Σ₁(s); η (Dedekind) als Wort ausgeschrieben (η ist das Massenverhältnis).

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/ramanujan_weitere.py` (Ausgabe `ramanujan_weitere.txt`)

Geprüft wurden Formeln, die über Serre-Ableitung, Differentialsystem, Hecke, Rogers-Ramanujan, Mock-Theta, Dedekinds Eta-Funktion, Hardy-Ramanujan, Master-Theorem, Lambert-Reihe und 1/π-Reihe hinausgehen. Jede ist numerisch nachgerechnet; aufgenommen ist, was eine Größe des Modells trifft.

## 1  Ramanujan-Summation der Rückläufe

```
ν(s) = Σ_{n≥0} e^{−ns} = 1/s + (1 + ζ(0)) + Σ_{k≥1} ζ(−k)(−s)^k/k! = 1/s + ½ + s/12 − s³/720 + …
```

Auf 15 Stellen gleich (s = 0,1; 0,5; 1). Die Glieder sind die regularisierten Summen der Rückläufe: ½ = 1 + ζ(0) (Σ1), 1/12 = −ζ(−1) (Σn = 1 + 2 + 3 + … = −1/12), −1/720 = −ζ(−3)/3! (Σn³). Für die gelesene Masse [korr.] heißt das M_schein/M = r/r* + ½ + r*/(12r) − (r*/r)³/720 + …: Das erste Korrekturglied zum Fernfeld ist dieselbe Zahl −1/12, die auch die Grundverschiebung q^{−1/24} und die Nullpunktenergie einer Schleife trägt. [korr.] Zahlennähe ohne Herleitung, trägt nicht.

## 2  Lambert-Reihe über alle Harmonischen

```
Σ₁(s) = Σσ₁(n)e^{−ns} = (1 − E₂)/24 = π²/(6s²) − 1/(2s) + 1/24 + O(e^{−4π²/s})
```

Auf 4 Stellen gleich bei s = 0,3 / 0,1 / 0,03. Zählen alle Harmonischen einer Bahn gleich, wächst die Verstärkung wie 1/s² statt 1/s: a/a₀ = s²(1 + Σ₁) → π²/6 = 1,645, eine konstante Beschleunigung ohne Tully-Fisher. Die Beschleunigungsrelation trägt nur die Grundschwingung (Σqⁿ, d = 1).

## 3  Bilaterale Summe 1ψ1

```
Σ_{n∈ℤ} (a;q)_n/(b;q)_n · zⁿ = (q, b/a, az, q/(az); q)_∞ / (b, q/a, z, b/(az); q)_∞,     |b/a| < |z| < 1
```

Nachgerechnet bei a = 0,9, b = 0,3, q = 0,6, z = 0,8: beide Seiten 0,252956080128. Für b = q enthält jedes Glied mit n < 0 den Faktor (1 − q⁰) = 0: Die beidseitige Summe endet von selbst bei n = 0 (Rest 10⁻³¹), und mit a = q bleibt Σ_{n≥0} zⁿ = 1/(1 − z). Die Formel ist damit die q-Form der Wahl zwischen einem und beiden Umlaufsinnen: Der eine Sinn entsteht, wenn der Nennerparameter gleich q ist, also durch eine Nullstelle beim Schritt vom Ort rückwärts. Welche physikalische Größe b entspricht, ist offen.

## 4  Ramanujan-Summen als Perioden-Projektoren

```
c_q(n) = Σ_{(a,q)=1} e^{2πian/q},     E_q = (1/N)·Σ_n c_q(n)·Tⁿ,     Σ_{q|N} E_q = I,     E_pE_q = δ_pq E_q,     Rang E_q = φ(q)
```

Nachgerechnet für N = 12 (Fehler 10⁻¹⁵). E₁ = (1/N)ΣTⁿ ist der Schließungsprojektor 𝒫; die übrigen E_q zerlegen seinen Rest nach genauer Periode q. Eine Bahn mit Re τ_M = p/q (Sternvieleck) liegt ganz in E_q. Filterprobe (Perioden 5 und 3 mit Rauschen, Länge 60): Energie 30,8 bei q = 5, 7,5 bei q = 3, sonst ≤ 1,1. Werkzeug, um Resonanzen in Bahn- oder Kartendaten nach ihrer Ordnung zu suchen.

## 5  Basen 2, 3, 4, 6 und Kaleidoskop

Ramanujan hat die Theorie der elliptischen Funktionen außer zur klassischen Basis zu den Signaturen 3, 4 und 6 entwickelt [G]. Gemeinsame Bedingung ist 2cos(2π/r) ∈ ℤ, erfüllt genau für r = 1, 2, 3, 4, 6; das sind die Drehordnungen, die ein kompakter flacher Raum zulässt. [korr.] Zahlennähe ohne Herleitung, trägt nicht. Für die Sechseckzelle gilt seine kubische Identität a³ = b³ + c³ (auf 10 Stellen bei q = 0,2 und 0,5), und a(q) = Σq^{m²+mn+n²} zählt die Bilder je Abstand: 6, 0, 6, 6, 0, 0, 12 für Abstand² = 1 … 7; die Quadratzelle mit Σq^{m²+n²}: 4, 4, 0, 4, 8. Die Zahl nächster Bilder ist die Drehordnung.

## Ergebnis

- Die Korrekturen zum Fernfeld sind die Ramanujan-Summen der Rückläufe; das erste Glied ist −ζ(−1) = 1/12.
- Nur die Grundschwingung gibt Tully-Fisher; alle Harmonischen gäben a → (π²/6)a₀.
- 1ψ1 ist die q-Form von einem gegen beide Umlaufsinne; der eine Sinn entsteht durch eine Nullstelle beim Rückschritt.
- Ramanujan-Summen zerlegen 𝒫 nach genauer Periode und sind ein Suchfilter für Resonanzen.
- Die Kaleidoskop-Winkel sind die Signaturen von Ramanujans alternativen Basen; seine kubischen Thetareihen zählen die Bilder der Sechseckzelle. [korr.] Zahlennähe ohne Herleitung, trägt nicht (gilt für die Gleichsetzung von Winkeln und Signaturen; das Zählen der Bilder durch die Thetareihen ist nachgerechnet).
