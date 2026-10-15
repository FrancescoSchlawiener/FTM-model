# Netzwerk-Summe und Rückstell-Bilanz, durchgerechnet

> **Korrigiert am 29.09.2026** gegen das Modell und den Rechenkern ftm_kern.py; Bezeichnungen und Konstanten nach der Notation in 00-forschungsstand.md. Korrigierte Stellen tragen [korr.]. Geändert: Straffung ΣΔs_i/r_i als offener Kandidat → gerechnet: braucht keine Mindestmasse, trifft die Form nicht (CPL w₀ = −0,31 … −0,47) und hängt an der Definition der Trichtergröße. [G] an Ω_b, „beobachtet 1·10⁻¹⁰“, R = 46,5 Mrd Lj, Alter 13,8 Gyr, H = 67–73 km/s/Mpc, q ≈ −0,55. Umbenannt: mitbewegter Teilchenhorizont χ_p → d_p; Bremsparameter q → q_B.

> Arbeitsstand, wartet auf die Freigabe des Autors.
> Rechnung: `rechnungen/netz_und_bilanz.py` (Ausgabe `netz_und_bilanz.txt`)
> Schwelle a₀ = 1,2·10⁻¹⁰ m/s² (gemessen). Kosmologie: Ω_b = 0,0493 [G], H₀ = 67,4 km/s/Mpc.

## A. Die Netzwerk-Summe

Aussage des Autors: Eine Galaxie im Haufen liegt im Trichter-Netzwerk aller anderen. An ihrem Außenrand liegt die Netzsumme S_Netz = Σ Δs_j, das 2833-Fache ihrer eigenen Zusatzstrecke. Der Raum ist dort „strammer gezogen“.

### A1. Linear in den Kern eingesetzt

Die Feldgleichung des Autors hat im Fernfeld a = (S/r)·(3c²/r*). Mit der eigenen Summe S_eigen = 2,95·10¹³ m gibt sie bei 15 kpc a = 6,7·10⁻¹¹ m/s², mit S_Netz = 8,37·10¹⁶ m dagegen 1,9·10⁻⁷ m/s². Beobachtet ist dort etwa 1·10⁻¹⁰ m/s² [G]. Die Netzsumme lässt sich also nicht linear in die Innendynamik einer Galaxie einsetzen, sie wäre 1800-mal zu groß.

### A2. Über den Operator

Die Bindung des Autors ist nicht additiv, sondern multiplikativ, und sie wirkt auf das Gesamtfeld. Das homogene Feld des Haufens fällt aus der Innendynamik heraus (Äquivalenzprinzip). Es bleibt aber im Argument des Schalters, denn χ hängt am **gesamten** Feld |g_innen + g_Haufen|. Das Haufenfeld hebt s an, χ sinkt, und die innere Verstärkung wird kleiner.

Die Rechnung setzt eine Galaxie von 6·10¹⁰ M☉ (Hernquist, a = 3 kpc) in einen Haufen mit 10¹⁴ M☉ Baryonen in 1 Mpc. Die innere Beschleunigung ist genähert durch g = ν(√(g_i,N² + g_e,N²)/a₀)·g_i,N.

| r [kpc] | isoliert [km/s] | bei 0,5 Mpc | bei 1 Mpc | bei 2 Mpc |
|---|---|---|---|---|
| 5 | 176,4 | 175,9 | 176,3 | 176,4 |
| 10 | 179,5 | 175,2 | 178,3 | 179,4 |
| 20 | 179,0 | 155,5 | 168,8 | 178,2 |
| 40 | 177,9 | 120,5 | 137,6 | 168,8 |
| 80 | 177,0 | 88,5 | 101,8 | 136,5 |

Das äußere Feld des Haufens (0,03–0,23 a₀) lässt die Rotationskurve außen fallen: bei 80 kpc auf 50–77 % des isolierten Werts. Innen bleibt sie unverändert.

**Folge.** Das Netzwerk wirkt auf die Galaxie **als Ganzes**: Sie bewegt sich in der Straffung des Haufens (dessen eigene Rotations- bzw. Geschwindigkeitsskala). Auf ihre Sterne wirkt es nur über die Gezeiten und über den Schalter. Außenbahnen in dichter Umgebung laufen dadurch **langsamer**, nicht schneller. Das ist eine prüfbare Umgebungsabhängigkeit: Galaxien gleicher Masse haben in Haufen und Gruppen außen fallende, im Feld flache Rotationskurven. In Rotationskurven-Katalogen ist ein solcher Effekt berichtet (Chae et al. 2020 [G]).

## B. Die Rückstell-Bilanz der Dunklen Energie

Aussage des Autors: σ = ΣΔs_innen − ΔL_außen = 0, also ΔL_außen/R = ΣΔs/R. Für das beobachtbare Universum (R = 46,5 Mrd Lj [G], sichtbare Materie) ist das 0,0844. Die Atemrate H ist die zeitliche Änderung dieser Rückstellung.

### B1. Die Rückstellung als Funktion der Zeit

Die Summe läuft über den Teilchenhorizont r_p = a·d_p [korr.] mit der sichtbaren Materie ρ_b = ρ_b,0/a³:

```
ε(t) = ΣΔs/r_p = (4πG ρ_b / 9c²) · r_p² = K · d_p²/a ,        K = 3,90·10⁻⁵ (Mrd Lj)⁻²      [korr.]
```

Heute gilt ε = 0,0831 (d_p = 46,1 Mrd Lj) [korr.]. Das ist die Zahl des Autors.

Weil Δs nur von der Masse abhängt, ändert sich ΣΔs innerhalb eines festen mitbewegten Gebiets nie. ε wächst nur, weil der Horizont neue Masse einschließt, und fällt, weil r_p mitwächst.

### B2. Im gemessenen Expansionsverlauf

Mit dem gemessenen Verlauf H(z) gilt heute dε/dt = −0,0021 pro Gyr gegen H₀ = 0,0689 pro Gyr. Die Rückstellung **nimmt heute ab**. Ihre Änderung ist 3 % von H₀ und hat das umgekehrte Vorzeichen.

### B3. Selbstkonsistent mit H = dε/dt

Setzt man die Regel des Autors als Bewegungsgleichung an, so folgt mit d(d_p)/dt = c/a [korr.]:

```
H = 2Kc·d_p / (a²(1 + ε))      [korr.]
```

Solange ε klein ist, beginnt die Lösung wie a ∝ t^{2/3}, also wie ein materiedominiertes Universum. Bei heutiger Baryonendichte (a = 1) ergibt sich:

| Größe | Bilanz des Autors | gemessen |
|---|---|---|
| Alter | 217 Gyr | 13,8 Gyr [G] |
| H | 2,9 km/s/Mpc | 67–73 km/s/Mpc [G] |
| Verzögerungsparameter q_B [korr.] | +0,53 (gebremst) | ≈ −0,55 (beschleunigt) [G] |

## Ergebnis

- **Netzwerk-Summe:** Linear addiert ist sie 1800-mal zu groß. Über den nicht-additiven Operator wirkt das Netzwerk als äußeres Feld. Es verlangsamt die Außenbahnen einer Galaxie im Haufen (bei 80 kpc auf 50–77 %), statt sie zu beschleunigen. Das ist eine eigenständige Vorhersage zur Umgebungsabhängigkeit.
- **Rückstell-Bilanz:** Der Wert ΔL/R = 0,084 stimmt. Als Zeitgesetz gibt die Bilanz aber keine beschleunigte Expansion: Die Rückstellung nimmt heute ab, und die selbstkonsistente Lösung ist gebremst mit einem H, das 23-mal zu klein ist.
- **Offen:** Welche Größe statt ΣΔs/r_p die Rückstellung tragen soll. Die Summe der Faltungstiefen ist fest, weil Δs nur an der Masse hängt. Eine wachsende Rückstellung braucht eine Größe, die sich mit der Strukturbildung ändert. Die Straffung ΣΔs_i/r_i der einzelnen Trichter, die beim Kollaps wächst, ist gerechnet: Sie braucht keine Mindestmasse, trifft die gemessene Form aber nicht (CPL w₀ = −0,31 … −0,47) und hängt an der Definition der Trichtergröße. [korr.]
