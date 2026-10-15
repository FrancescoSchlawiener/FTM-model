"""Die Elemente: Messdaten (PubChem, AME2020), Rechnung des Modells, mögliche neue Elemente.
Schreibt abhandlung/elemente_de.md und elemente_en.md; der Lauf prüft, dass beide genau dieser Rechnung entsprechen."""
import csv
import math
import os

import numpy as np

import kern as K
import daten as D
import kerne as KE

WURZEL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZIEL = {"de": os.path.join(WURZEL, "abhandlung", "elemente_de.md"), "en": os.path.join(WURZEL, "abhandlung", "elemente_en.md")}
KOMMA = [","]
L = "spdfg"


def k(x, n):
    """Zahl mit n Nachkommastellen; Dezimalzeichen nach Sprache (Komma deutsch, Punkt englisch)."""
    return f"{x:.{n}f}".replace(".", KOMMA[0]).replace("-", "−")


def konfig_nl(Z):
    """Konfiguration nach n + l als Text, ab dem letzten Edelgas."""
    edel = {2: "He", 10: "Ne", 18: "Ar", 36: "Kr", 54: "Xe", 86: "Rn", 118: "Og"}
    rumpf = max([e for e in edel if e < Z], default=0)
    teile = [f"{n}{L[l]}{min(Z, b) - a + 1}" for n, l, a, b in K.fuellung(8) if rumpf < a <= Z]
    return (f"[{edel[rumpf]}] " if rumpf else "") + " ".join(teile)


def konfig_text(Z, b):
    """Gerechnete Konfiguration als Text, ab dem letzten Edelgas nach n + l (bei gleichem n + l zuerst kleineres n)."""
    edel = {2: "He", 10: "Ne", 18: "Ar", 36: "Kr", 54: "Xe", 86: "Rn", 118: "Og"}
    rumpf = max([e for e in edel if e < Z], default=0)
    kern = {(n, l): b - a + 1 for n, l, a, b in K.fuellung(9) if b <= rumpf}
    teile = [f"{n}{L[l]}{v}" for (n, l), v in sorted(b.items(), key=lambda x: (x[0][0] + x[0][1], x[0][0])) if kern.get((n, l)) != v]
    return (f"[{edel[rumpf]}] " if rumpf and all(b.get(k) == v for k, v in kern.items()) else "") + " ".join(teile)


def entstehung():
    """Konfigurationen aus der Etappe entstehung (Hülle beim Fallen von δ): {Z: {(n, l): Zahl}}."""
    from huelle_ablage import lesen
    return {int(z): {(n, l): m for n, l, m in d["konf"]} for z, d in lesen().get("_entstehung", {}).items()}


def zaehlweisen():
    """Konfigurationen der Suche nach den zwei Zählweisen, von n + l aus ohne Messung (Etappen suche, ueber_118):
    {Z: {(n, l): Zahl}}."""
    from huelle_ablage import lesen
    return {int(z): {(n, l): m for n, l, m in d if m} for z, d in lesen().get("_suche", {}).items()}


def term_text(g):
    """Grundterm ^{2S+1}L_J als Text und Moment g·J [μ_B]."""
    from fractions import Fraction
    if g["J"] == 0 and g["S"] == 0 and g["L"] == 0:
        return "¹S₀ · 0"
    J = Fraction(g["J"]).limit_denominator(2)
    hoch = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    tief = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return f"{str(int(round(2 * g['S'] + 1))).translate(hoch)}{'SPDFGHIKLMNOQRTUVWXYZ'[int(round(g['L']))]}{str(J).translate(tief)} · {k(g['mu'], 2)}"


def systemname(Z):
    """Systematischer Name und Symbol nach IUPAC aus den Ziffern (nil, un, bi, tri, quad, pent, hex, sept, oct, enn)."""
    w = ["nil", "un", "bi", "tri", "quad", "pent", "hex", "sept", "oct", "enn"]
    name = "".join(w[int(d)] for d in str(Z)) + "ium"
    name = name.replace("iium", "ium").replace("nnn", "nn")
    return name.capitalize(), "".join(w[int(d)][0] for d in str(Z)).capitalize()


def letzte_schale(Z):
    n, l, a, b = [x for x in K.fuellung(8) if x[2] <= Z <= x[3]][0]
    return f"{n}{L[l]}"


def periode(Z):
    edel = (0, 2, 10, 18, 36, 54, 86, 118, 168)
    p = max(i for i in range(len(edel)) if edel[i] < Z) + 1
    return p, Z - edel[p - 1], edel[p] - edel[p - 1]


TR, SC = KE.TROPFEN, KE.SCHALE


def dauer(lg, sprache):
    """Zeit 10^lg s in der passenden Einheit."""
    t = 10**lg
    for grenze, faktor, de, en in ((86400, 86400, "Tage", "days"), (3600, 3600, "Stunden", "hours"), (60, 60, "Minuten", "minutes"), (0, 1, "Sekunden", "seconds")):
        if t >= grenze:
            return f"{k(t / faktor, 0) if t / faktor >= 10 else k(t / faktor, 1)} {de if sprache == 'de' else en}"


def tropfen():
    p = KE.anpassung(TR)[0]
    return dict(aV=p[0], aS=p[1], aC=p[2], aA=p[3], aP=p[4])


def tal_n(Z):
    return KE.n_tal(Z, TR)


def x_kipp(Z, A):
    return KE.x(Z, A, TR)


def x_schale(Z):
    """x mit Schalenglied im Tal dieser Formel."""
    return KE.x(Z, KE.a_tal(Z, SC), SC)


def b_je_a(Z, A):
    """Bindung je Nukleon nach der Tröpfchenformel [keV]."""
    return float(KE.bindung(Z, A, TR)[0]) / A * 1e3


def kerne_je_z():
    """Je Z aus AME2020: Zahl gemessen, Zahl gesamt, kleinstes und größtes gemessenes A, festester gemessener Kern (A, B/A [keV])."""
    aus = {}
    for Z, N, A, B, gem in D.ame2020_markiert():
        if Z < 1:
            continue
        e = aus.setdefault(Z, dict(gem=0, ges=0, amin=None, amax=None, fest=None))
        e["ges"] += 1
        if gem:
            e["gem"] += 1
            e["amin"] = A if e["amin"] is None else min(e["amin"], A)
            e["amax"] = A if e["amax"] is None else max(e["amax"], A)
            if e["fest"] is None or B > e["fest"][1]:
                e["fest"] = (A, B)
    return aus


def zahlen():
    """Kennzahlen für Text und Prüfung, alle aus der Simulation der Kerne."""
    rows = list(KE._pubchem_zeilen())
    edel, abw, n_el = KE.konfigurationen()
    ame = D.ame2020_markiert()
    return dict(rows=rows, edel=edel, abw=list(abw), n_el=n_el, t=tropfen(), kz=kerne_je_z(), lam_H=KE.halterlaenge(KE.a_C(TR)),
                n_kerne=sum(1 for x in ame if x[0] >= 1), n_gem=sum(1 for x in ame if x[0] >= 1 and x[4]),
                z_gem_max=KE.bekannte()["z_gem_max"], eins_durch_alpha=1 / K.alpha)


def text(sprache="de"):
    KOMMA[0] = "," if sprache == "de" else "."
    T = lambda de, en: de if sprache == "de" else en
    Zd = zahlen()
    rows, t, kz = Zd["rows"], Zd["t"], Zd["kz"]
    sym = {int(r["AtomicNumber"]): r["Symbol"] for r in rows}
    gem = KE.belegung_gemessen()
    nl = __import__("huelle_1_loeser").kerne_belegung_nl
    abw = dict(Zd["abw"])
    art = {"schliesst": T("schließt halbe oder volle Schale", "closes a half or full shell"),
           "d_vor_f": T("d vor f", "d before f"), "s_an_d": T("s an d", "s to d")}
    gr = KE.grenzen()
    at = gr["tropfen"][1]
    zsicher = gr["z_sicher"]                                        # Kippstelle beider Formeln
    ks = gr["schale"]                                               # Kippstelle mit Schalenglied
    zkr = gr["huelle"]                                              # Grenze der Hülle
    zmax = gr["z_max"]                                              # letztes Element unter der Grenze der Hülle
    neu = range(119, zmax + 1)
    von = T("von", "of")
    o = []
    w = o.append
    w(T("# Die Elemente: gemessen, gerechnet, möglich\n", "# The Elements: Measured, Computed, Possible\n"))
    w(T("Diese Datei stellt alle bekannten Elemente mit ihren Messdaten neben die reine Rechnung des Modells und leitet daraus ab, "
        "welche Elemente über die bekannten hinaus möglich sind. Jede Zahl ist aus den Daten gelesen oder mit dem Kern der Simulation "
        "gerechnet; die Datei wird von der Simulation geschrieben und bei jedem Lauf mit ihr verglichen.\n",
        "This file places all known elements with their measured data next to the pure calculation of the model and derives from it "
        "which elements are possible beyond the known ones. Every number is read from the data or computed with the core of the simulation; "
        "the file is written by the simulation and compared with it on every run.\n"))
    w(T(f"**Daten.** PubChem: {Zd['n_el']} Elemente mit Masse, Elektronenkonfiguration, Ionisierungsenergie, Atomradius, "
        f"Elektronegativität, Schmelz- und Siedepunkt, Dichte. AME2020: {Zd['n_kerne']} Kerne, davon {Zd['n_gem']} mit gemessener "
        f"Masse; gemessen bis Z = {Zd['z_gem_max']}, darüber sind die Massen geschätzt.\n",
        f"**Data.** PubChem: {Zd['n_el']} elements with mass, electron configuration, ionization energy, atomic radius, "
        f"electronegativity, melting and boiling point, density. AME2020: {Zd['n_kerne']} nuclei, {Zd['n_gem']} of them with measured "
        f"mass; measured up to Z = {Zd['z_gem_max']}, above that the masses are estimated.\n"))
    w("---\n")
    w(T("## 1 Was gerechnet wird\n", "## 1 What Is Computed\n"))
    w(T("Das Modell rechnet die Elemente auf zwei Ebenen, Hülle und Kern, mit je einer schließenden Form und einer Kopplung, die sie öffnet.\n",
        "The model computes the elements on two levels, shell and nucleus, each with a closing form and a coupling that opens it.\n"))
    w(T("**Hülle.** Die reine 1/r²-Neigung schließt nach 2n² je Schale. Die inneren Schalen öffnen diese Schließung: "
        "Bahnen mit kleiner Windung tauchen tief ein und binden fester, die Schalen füllen nach n + l, bei gleichem n + l zuerst die mit kleinerem n. "
        "Die Kopplung der Hülle ist Zα: das Tempo des innersten Elektrons gegen c. Bei Zα = 1 erreicht es c, der Treffpunkt der Hülle, "
        f"bei Z = 1/α = {k(Zd['eins_durch_alpha'], 3)}, für einen punktförmigen Kern.\n",
        "**Shell.** The pure 1/r² slope closes after 2n² per shell. The inner shells open this closure: "
        "orbits with small winding dip deep and bind more tightly, the shells fill by n + l, for equal n + l the one with smaller n first. "
        "The coupling of the shell is Zα: the speed of the innermost electron relative to c. At Zα = 1 it reaches c, the meeting point of the shell, "
        f"at Z = 1/α = {k(Zd['eins_durch_alpha'], 3)}, for a point nucleus.\n"))
    kt, kg, tvt, tvs = KE.kippstelle(TR), KE.kippstelle_ganz(TR), KE.tal_vergleich(TR), KE.tal_vergleich(SC)
    ks6, ks6g = KE.kippstelle(SC), KE.kippstelle_ganz(SC)
    w(T("**Kern.** Der Kern ist ein gehaltener Trichter. Seine Bindung folgt der Tröpfchenformel, einmal angepasst an alle "
        f"{KE.anpassung(TR)[2]} gemessenen Kerne von AME2020 ab A = 12 (geschätzte Massen nicht)\n",
        "**Nucleus.** The nucleus is a held funnel. Its binding follows the liquid-drop formula, fitted once to all "
        f"{KE.anpassung(TR)[2]} measured nuclei of AME2020 from A = 12 (estimated masses excluded)\n"))
    w("```\nB = a_V·A − a_S·A^{2/3} − a_C·Z(Z − 1)/A^{1/3} − a_A·(N − Z)²/A ± a_P/√A\n```\n")
    w(T(f"mit a_V = {k(t['aV'], 2)}, a_S = {k(t['aS'], 2)}, a_C = {k(t['aC'], 4)}, a_A = {k(t['aA'], 2)}, a_P = {k(t['aP'], 2)} MeV, "
        f"Streuung {k(KE.anpassung(TR)[1], 2)} MeV. Das Coulomb-Glied, die Windung der Protonen, gibt die Halterlänge λ_H = (3/5)·αħc/a_C = "
        f"{k(Zd['lam_H'], 2)} fm. Die Spaltbarkeit x = E_C/(2E_S) = a_C·Z²/(2a_S·A) vergleicht Windung und Haltung; bei x = 1 verschwindet "
        "die Rückstellung einer Verformung, die Kippstelle des Kerns.\n",
        f"with a_V = {k(t['aV'], 2)}, a_S = {k(t['aS'], 2)}, a_C = {k(t['aC'], 4)}, a_A = {k(t['aA'], 2)}, a_P = {k(t['aP'], 2)} MeV, "
        f"scatter {k(KE.anpassung(TR)[1], 2)} MeV. The Coulomb term, the winding of the protons, gives the holder length λ_H = (3/5)·αħc/a_C = "
        f"{k(Zd['lam_H'], 2)} fm. The fissility x = E_C/(2E_S) = a_C·Z²/(2a_S·A) compares winding and holding; at x = 1 the restoring "
        "of a deformation vanishes, the tipping point of the nucleus.\n"))
    w(T("**Tal der Stabilität.** Bei festem A ist der Kern mit der kleinsten Faltungsenergie Z·E_H + N·E_n − B der stabilste, mit "
        "E = 3c⁴·Δs/G für das Wasserstoffatom und das Neutron. Gegen das leichteste gemessene Isobar je A (A = 20 bis 238, "
        f"{tvt[2]} Massenzahlen) trifft das Tal der Tröpfchenformel, je A das ganze Z kleinster Faltungsenergie, im Mittel auf {k(tvt[0], 2)} Protonen, {k(100 * tvt[1], 0)} % genau. "
        f"Im Tal liegt die Kippstelle x = 1 bei Z = {k(kt[0], 1)}, A = {k(kt[1], 1)}; auf einem zweiten Weg, je A das ganze Z mit der "
        f"kleinsten Faltungsenergie, bei Z = {kg[0]}, A = {kg[1]} bis {kg[2]}.\n",
        "**Valley of stability.** For fixed A the nucleus with the smallest fold energy Z·E_H + N·E_n − B is the most stable, with "
        "E = 3c⁴·Δs/G for the hydrogen atom and the neutron. Against the lightest measured isobar per A (A = 20 to 238, "
        f"{tvt[2]} mass numbers) the valley of the liquid-drop formula, taking for each A the integer Z of smallest fold energy, is accurate to {k(tvt[0], 2)} protons on average, {k(100 * tvt[1], 0)} % exact. "
        f"In the valley the tipping point x = 1 lies at Z = {k(kt[0], 1)}, A = {k(kt[1], 1)}; by a second route, taking for each A the integer Z "
        f"with the smallest fold energy, at Z = {kg[0]}, A = {kg[1]} to {kg[2]}.\n"))
    w(T("**Schalenglied.** Die Tröpfchenformel kennt keine Schließungen. Mit den Schließzahlen des Modells "
        f"({', '.join(map(str, KE.schliesszahlen_kern()[:-1]))}) als Schalenglied, in der Form von Myers und Swiatecki, mit Wigner- und "
        "Oberflächen-Asymmetrieglied, ebenso einmal angepasst, trifft die Formel die Daten besser: die Bindung mit "
        f"{k(KE.anpassung(SC)[1], 2)} statt {k(KE.anpassung(TR)[1], 2)} MeV Streuung, das leichteste gemessene Isobar je A auf "
        f"{k(tvs[0], 2)} statt {k(tvt[0], 2)} Protonen ({k(100 * tvs[1], 0)} statt {k(100 * tvt[1], 0)} % genau). Ihr Oberflächenglied ist größer, "
        f"a_S = {k(KE.a_S(SC), 2)} MeV, a_C = {k(KE.a_C(SC), 4)} MeV; die Kippstelle im Tal liegt damit bei Z = {k(ks6[0], 1)}, A = {k(ks6[1], 1)} "
        f"(ganzzahlig Z = {ks6g[0]}). Die Kippstelle selbst ist an den Daten nicht prüfbar: Gemessen sind Massen bis Z = {Zd['z_gem_max']}, "
        f"Spaltungen bis x = {k(KE.spaltung_gegen_x(TR)[3], 3)} (x der Tröpfchenformel). Die Tabellen geben darum beide Formeln.\n",
        "**Shell term.** The liquid-drop formula knows no closures. With the closure numbers of the model "
        f"({', '.join(map(str, KE.schliesszahlen_kern()[:-1]))}) as a shell term, in the form of Myers and Swiatecki, with Wigner and "
        "surface-asymmetry terms, likewise fitted once, the formula meets the data better: the binding with "
        f"{k(KE.anpassung(SC)[1], 2)} instead of {k(KE.anpassung(TR)[1], 2)} MeV scatter, the lightest measured isobar per A to "
        f"{k(tvs[0], 2)} instead of {k(tvt[0], 2)} protons ({k(100 * tvs[1], 0)} instead of {k(100 * tvt[1], 0)} % exact). Its surface term is larger, "
        f"a_S = {k(KE.a_S(SC), 2)} MeV, a_C = {k(KE.a_C(SC), 4)} MeV; the tipping point in the valley thus lies at Z = {k(ks6[0], 1)}, A = {k(ks6[1], 1)} "
        f"(integer Z = {ks6g[0]}). The tipping point itself cannot be tested on the data: masses are measured up to Z = {Zd['z_gem_max']}, "
        f"fissions up to x = {k(KE.spaltung_gegen_x(TR)[3], 3)} (x of the liquid-drop formula). The tables therefore give both formulas.\n"))
    w("---\n")
    w(T("## 2 Die Leiter der Perioden\n", "## 2 The Ladder of the Periods\n"))
    w(T("Die Länge der Periode p ist die Zahl der Zustände mit n + l = p + 1, 2⌈(p + 1)/2⌉². Jede Stufe 2k² kommt zweimal, weil n + l = 2k − 1 "
        "und 2k dieselbe größte Windung l = k − 1 haben.\n",
        "The length of period p is the number of states with n + l = p + 1, 2⌈(p + 1)/2⌉². Every step 2k² appears twice, because n + l = 2k − 1 "
        "and 2k have the same largest winding l = k − 1.\n"))
    w(T("| Periode | gemessen schließt (Edelgas) | Länge gemessen | 2⌈(p + 1)/2⌉² | Schalen nach n + l |",
        "| Period | measured closure (noble gas) | measured length | 2⌈(p + 1)/2⌉² | shells by n + l |"))
    w("|---|---|---|---|---|")
    edel = (0,) + Zd["edel"]
    for p in range(1, 9):
        lang = 2 * math.ceil((p + 1) / 2)**2
        schalen = " ".join(f"{n}{L[l]}" for n, l, a, b in K.fuellung(8) if (edel[p - 1] if p <= 7 else 118) < a <= (edel[p] if p <= 7 else 168))
        if p <= 7:
            w(f"| {p} | {edel[p]} ({sym[edel[p]]}) | {edel[p] - edel[p - 1]} | {lang} | {schalen} |")
        else:
            w(f"| 8 | {T('(168, gerechnet)', '(168, computed)')} | – | {lang} | {schalen} |")
    w("")
    import huelle_3_eigenschaften as SD3
    w(T("**Die Leiter als Spiegelform.** Ein Feld, das in log r um einen Radius R gespiegelt ist (r ↔ R²/r), in log r v·f(1 − f) mit "
        "f = 1/(eˣ + 1), bindet an der Schwelle E = 0 genau bei den Kopplungen N(N + 1), N = n + l, jede Stufe mit 2⌈N/2⌉² Zuständen. "
        "Die Zustände des Spiegelfelds an der Schwelle sind Kugelfunktionen, u = √r·P_N^{2l+1}((r − R)/(r + R)) (exakt; numerisch "
        "Überlapp über 0,9999): Grad N = n + l, Ordnung 2l + 1, Pole Kern und Ferne, Äquator R.\n",
        "**The ladder as a mirror form.** A field mirrored in log r about a radius R (r ↔ R²/r), in log r v·f(1 − f) with "
        "f = 1/(eˣ + 1), binds at the threshold E = 0 exactly at the couplings N(N + 1), N = n + l, each step with 2⌈N/2⌉² states. "
        "The states of the mirror field at the threshold are spherical harmonics, u = √r·P_N^{2l+1}((r − R)/(r + R)) (exact; numerically "
        "overlap above 0.9999): degree N = n + l, order 2l + 1, poles nucleus and distance, equator R.\n"))
    w("")
    w("---\n")
    w(T("## 3 Die bekannten Elemente\n", "## 3 The Known Elements\n"))
    w(T("### 3.1 Messdaten: Atom und Stoff (PubChem)\n", "### 3.1 Measured Data: Atom and Substance (PubChem)\n"))
    w(T("Masse in u, Ionisierungsenergie IE in eV, Atomradius r in pm, Elektronegativität EN (Pauling), Schmelz- und Siedepunkt in K, "
        "Dichte ρ in g/cm³. Ein Strich: kein Messwert in PubChem.\n",
        "Mass in u, ionization energy IE in eV, atomic radius r in pm, electronegativity EN (Pauling), melting and boiling point in K, "
        "density ρ in g/cm³. A dash: no measured value in PubChem.\n"))
    w(T("| Z | | Name | Masse | Konfiguration | IE | r | EN | T_Schmelz | T_Sied | ρ |",
        "| Z | | Name | Mass | Configuration | IE | r | EN | T_melt | T_boil | ρ |"))
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        f = lambda s: s if s not in ("", None) else "–"
        z = (f"| {r['AtomicNumber']} | {r['Symbol']} | {r['Name']} | {f(r['AtomicMass'])} | {r['ElectronConfiguration']} | "
             f"{f(r['IonizationEnergy'])} | {f(r['AtomicRadius'])} | {f(r['Electronegativity'])} | {f(r['MeltingPoint'])} | "
             f"{f(r['BoilingPoint'])} | {f(r['Density'])} |")
        w(z.replace(".", ",") if sprache == "de" else z)
    w("")
    w(T("### 3.2 Messdaten: Kern (AME2020)\n", "### 3.2 Measured Data: Nucleus (AME2020)\n"))
    w(T("Zahl der Kerne (gemessen von gesamt), Bereich der gemessenen Massenzahlen, festester gemessener Kern mit seiner Bindung je Nukleon "
        "B/A in keV und seinem N/Z.\n",
        "Number of nuclei (measured of total), range of measured mass numbers, most tightly bound measured nucleus with its binding per nucleon "
        "B/A in keV and its N/Z.\n"))
    w(T("| Z | | Kerne | gemessen A | festester Kern | B/A | N/Z |", "| Z | | nuclei | measured A | most bound nucleus | B/A | N/Z |"))
    w("|---|---|---|---|---|---|---|")
    for Z in range(1, 119):
        e = kz.get(Z)
        if e is None or not e["gem"]:
            w(f"| {Z} | {sym[Z]} | 0 {von} {e['ges'] if e else 0} | – | – | – | – |")
            continue
        A, B = e["fest"]
        w(f"| {Z} | {sym[Z]} | {e['gem']} {von} {e['ges']} | {e['amin']}–{e['amax']} | {sym[Z]}-{A} | {k(B, 1)} | {k((A - Z) / Z, 3)} |")
    w("")
    w(T("### 3.3 Rechnung des Modells\n", "### 3.3 Calculation of the Model\n"))
    w(T("Periode und Platz, letzte Schale nach n + l, Vergleich mit der gemessenen Konfiguration, Kopplung der Hülle Zα, "
        "Bindung je Nukleon des festesten gemessenen Kerns nach der Tröpfchenformel (gemessen in 3.2), Tal der Stabilität (A, N) und Spaltbarkeit x im Tal, "
        "nach der Tröpfchenformel und mit Schalenglied (mit dessen Tal).\n",
        "Period and place, last shell by n + l, comparison with the measured configuration, coupling of the shell Zα, "
        "binding per nucleon of the most bound measured nucleus by the liquid-drop formula (measured in 3.2), valley of stability (A, N) and fissility x "
        "in the valley, by the liquid-drop formula and with the shell term (with its valley).\n"))
    w(T("Magnetismus des Atoms: Grundterm nach Hund aus der gemessenen Konfiguration (größtes S, dann größtes L je offener Unterschale, "
        "summiert; J = |L − S| unter halber Füllung der offenen Schale größter Windung, sonst L + S) und Moment g·J in Bohr-Magnetonen "
        "(Landé-Faktor, g_s = 2). Gemessene Momente enthalten unsere Daten nicht.\n",
        "Magnetism of the atom: ground term by Hund from the measured configuration (largest S, then largest L per open subshell, "
        "summed; J = |L − S| below half filling of the open shell of largest winding, otherwise L + S) and moment g·J in Bohr magnetons "
        "(Landé factor, g_s = 2). Our data contain no measured moments.\n"))
    w(T("| Z | | Periode · Platz | Schale | n + l wie gemessen | Entstehung wie gemessen | zwei Zählweisen wie gemessen | Grundterm · Moment [μ_B] | Zα | B/A gerechnet | A_Tal | N_Tal | x Tropfen | x mit Schale |",
        "| Z | | period · place | shell | n + l as measured | formation as measured | two countings as measured | ground term · moment [μ_B] | Zα | B/A computed | A_valley | N_valley | x drop | x with shell |"))
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    ent = entstehung()
    wie = lambda Z: "–" if Z not in ent else (T("ja", "yes") if ent[Z] == {q: v for q, v in gem[Z].items() if v} else konfig_text(Z, ent[Z]))
    zw = zaehlweisen()
    wie2 = lambda Z: "–" if Z not in zw else (T("ja", "yes") if zw[Z] == {q: v for q, v in gem[Z].items() if v} else konfig_text(Z, zw[Z]))
    xs = x_schale
    for Z in range(1, 119):
        p, pl, lg = periode(Z)
        n = tal_n(Z)
        A = n + Z
        gl = T("ja", "yes") if Z not in abw else art[abw[Z]]
        e = kz.get(Z)
        bg = k(b_je_a(Z, e["fest"][0]), 1) if e and e["fest"] and e["fest"][0] >= 12 else "–"
        w(f"| {Z} | {sym[Z]} | {p} · {pl} {von} {lg} | {letzte_schale(Z)} | {gl} | {wie(Z)} | {wie2(Z)} | {term_text(SD3.grundterm(gem[Z]))} | {k(Z * K.alpha, 3)} | {bg} | {k(A, 1)} | {k(n, 1)} | {k(x_kipp(Z, A), 3)} | {k(xs(Z), 3) if xs(Z) == xs(Z) else '–'} |")
    w("")
    w(T(f"### 3.4 Die {len(abw)} Abweichungen von n + l\n", f"### 3.4 The {len(abw)} Deviations from n + l\n"))
    w(T(f"{Zd['n_el'] - len(abw)} der {Zd['n_el']} Elemente füllen genau nach n + l. Alle Abweichungen liegen in d- und f-Schalen innerhalb "
        "der Perioden; keine ändert eine Schließzahl.\n",
        f"{Zd['n_el'] - len(abw)} of the {Zd['n_el']} elements fill exactly by n + l. All deviations lie in d and f shells within "
        "the periods; none changes a closure number.\n"))
    w(T("| Z | | gemessen (PubChem) | nach n + l | Art | Entstehung | zwei Zählweisen |", "| Z | | measured (PubChem) | by n + l | kind | formation | two countings |"))
    w("|---|---|---|---|---|---|---|")
    gleich = lambda x: x if x != T("ja", "yes") else T("wie gemessen", "as measured")
    for Z, a in Zd["abw"]:
        w(f"| {Z} | {sym[Z]} | {rows[Z - 1]['ElectronConfiguration']} | {konfig_nl(Z)} | {art[a]} | {gleich(wie(Z))} | {gleich(wie2(Z))} |")
    w("")
    zahl = {a: sum(1 for _, b in Zd["abw"] if b == a) for a in art}
    w(T(f"{zahl['schliesst']} schließen eine halbe oder volle Schale (d⁵, d¹⁰, f⁷), {zahl['d_vor_f']} beginnen die f-Schale mit einem "
        f"d-Elektron, der kleineren Windung, {zahl['s_an_d']} geben ein s-Elektron an die d-Schale.\n",
        f"{zahl['schliesst']} close a half or full shell (d⁵, d¹⁰, f⁷), {zahl['d_vor_f']} begin the f shell with a "
        f"d electron, the smaller winding, {zahl['s_an_d']} pass an s electron to the d shell.\n"))
    gz = [Z for Z in sorted(ent) if Z <= 102]
    tr = [Z for Z in gz if wie(Z) == T("ja", "yes")]
    az = [Z for Z, _ in Zd["abw"] if Z in gz]
    ta = [Z for Z in az if Z in tr]
    da = ", ".join(sym[Z] for Z in gz if Z not in tr)
    w(T("**Die Hülle aus ihrer Entstehung.** Gerechnet wird, wie die Hülle entsteht, während δ des Ganzen fällt: Der Kern, ein gehaltener "
        "Trichter, fängt mit k gebundenen Elektronen das nächste dort, wo die Faltung der ganzen Hülle am tiefsten wird. Jede Bahn schließt je "
        "Umlauf in der Neigung des Kerns und der übrigen Windungen, mit dem Takt der Abwicklung, den das Zentrum setzt. Gleichsinnige Windungen "
        "sperren einander ihren Zustand (Schließungsoperator 𝒫), gezählt über die Moden der Kugelfunktionen. Geschützt ist die Zählung, nicht die "
        "einzelne Schleife: Wo mindestens zwei Windungen einrasten und die Zählregel schließt (d⁵, d¹⁰, f⁷, volle Schalen), bleibt gebunden; die "
        "übrigen wandern nach jedem Fang, solange die Faltung des Ganzen sinkt. Kein Wert ist angepasst; die Messung prüft nur am Ende. "
        f"Von {len(gz)} gerechneten Elementen mit gemessener Konfiguration (PubChem) trifft die Entstehung {len(tr)}, von den {len(az)} "
        f"Abweichungen von n + l darunter {len(ta)}. Daneben liegen {da}. Offen sind die Paarung gegensinniger Windungen in derselben Mode und die "
        "Bindung von f gegen d.\n",
        "**The shell from its formation.** What is computed is how the shell forms while δ of the whole falls: the nucleus, a held funnel, "
        "holding k electrons captures the next one where the folding of the whole shell becomes deepest. Every orbit closes per revolution in "
        "the slope of the nucleus and of the other windings, at the clock rate of unwinding set by the centre. Windings of the same sense block "
        "each other's state (closure operator 𝒫), counted over the modes of the spherical harmonics. Protected is the count, not the single loop: "
        "where at least two windings lock and the counting rule closes (d⁵, d¹⁰, f⁷, full shells), it stays bound; the others move after every "
        "capture as long as the folding of the whole falls. No value is fitted; the measurement only checks at the end. "
        f"Of {len(gz)} computed elements with measured configuration (PubChem) the formation meets {len(tr)}, of the {len(az)} deviations from "
        f"n + l among them {len(ta)}. Off are {da}. Open are the pairing of windings of opposite sense in the same mode and the binding of f "
        "against d.\n"))
    zz = [Z for Z in sorted(zw) if Z in gem]
    zt = [Z for Z in zz if wie2(Z) == T("ja", "yes")]
    za = [Z for Z, _ in Zd["abw"] if Z in zt]
    w(T("**Die Hülle aus zwei Zählweisen.** Die Kopplung der Windungen wird auf zwei Arten gezählt: örtlich, aus der Dichte am Ort, und je "
        "Umlaufsinn, jede Windung gegen jede gleichsinnige. Welche Zählung gilt, setzt die Füllung der Spinhälfte der wechselnden Schale größter "
        "Windung, gewichtet mit dem Verhältnis der Tempi der beiden Bahnen, zwischen denen das Elektron wechselt (Fermi-Form). Der Abstand zweier "
        "Konfigurationen ist um die Mitte gespiegelt, das Mittel beider Blickrichtungen, am eingerasteten Niveau J beider Umlaufsinne. Von n + l "
        "aus wechselt ein Elektron nach dem anderen zur tieferen Nachbarin, ohne Kenntnis der Messung. So gerechnet treffen "
        f"{len(zt)} von {len(zz)} Elementen mit Konfiguration in PubChem (1–102 gemessen, 103–118 vorhergesagt) die Konfiguration, darunter "
        f"{len(za)} der {len(Zd['abw'])} Abweichungen von n + l. Gegen das Gitter null gerechnet liegt einzig Np knapp, um 0,01 eV, auf der "
        "Seite von 5f⁵7s².\n",
        "**The shell from two countings.** The coupling of the windings is counted in two ways: locally, from the density at the place, and per "
        "sense of revolution, each winding against each of the same sense. Which counting holds is set by the filling of the spin half of the "
        "changing shell of largest winding, weighted by the ratio of the tempos of the two orbits between which the electron changes (Fermi form). "
        "The distance of two configurations is mirrored about the middle, the mean of both directions of view, at the locked level J of both "
        "senses. From n + l one electron after the other moves to the lower neighbour, without knowledge of the measurement. Computed so, "
        f"{len(zt)} of {len(zz)} elements with a configuration in PubChem (1–102 measured, 103–118 predicted) meet the configuration, among them "
        f"{len(za)} of the {len(Zd['abw'])} deviations from n + l. Computed towards zero grid step, only Np lies narrowly, by 0.01 eV, on the "
        "side of 5f⁵7s².\n"))
    w(T("### 3.5 Was die Rechnung über die bekannten sagt\n", "### 3.5 What the Calculation Says about the Known Elements\n"))
    x294 = x_kipp(118, 294)
    abw_b = [b_je_a(Z, kz[Z]["fest"][0]) - kz[Z]["fest"][1] for Z in kz if kz[Z]["fest"] and kz[Z]["fest"][0] >= 12]
    mb = k(sum(abs(x) for x in abw_b) / len(abw_b), 0)
    zg = Zd["z_gem_max"]
    ng = max(x[1] for x in D.ame2020_markiert() if x[0] == zg and x[4])
    w(T(f"- Die Edelgase (PubChem) schließen nach {', '.join(map(str, Zd['edel']))}; nach n + l gerechnet genau dieselben Zahlen.\n"
        f"- Die Formel trifft die Bindung je Nukleon der festesten gemessenen Kerne (A ≥ 12, {len(abw_b)} Elemente) im Mittel auf "
        f"{mb} keV; die Abweichungen tragen die Schließungen, die die Formel nicht enthält. Unter A = 12 ist die Formel nicht angepasst.\n"
        f"- Die schwersten bekannten Kerne sind neutronenarm: Og-294 (Masse geschätzt) hat N = 176, das Tal verlangt bei Z = 118 "
        f"N = {k(tal_n(118), 0)}. Darum erreicht Og-294 nach der Tröpfchenformel die Kippstelle schon, x = {k(x294, 2)}; im Tal läge Z = 118 bei "
        f"x = {k(x_kipp(118, 118 + tal_n(118)), 2)}, mit Schalenglied bei {k(xs(118), 2)} (Og-294: {k(KE.x(118, 294, SC), 2)}).\n"
        f"- Gemessen sind die Massen bis Z = {zg}. Der neutronenreichste gemessene Kern dort hat N = {ng}, das Tal verlangt N = {k(tal_n(zg), 0)}.\n",
        f"- The noble gases (PubChem) close after {', '.join(map(str, Zd['edel']))}; computed by n + l, exactly the same numbers.\n"
        f"- The formula meets the binding per nucleon of the most bound measured nuclei (A ≥ 12, {len(abw_b)} elements) on average to "
        f"{mb} keV; the deviations carry the closures, which the formula does not contain. Below A = 12 the formula is not fitted.\n"
        f"- The heaviest known nuclei are neutron-poor: Og-294 (mass estimated) has N = 176, the valley requires N = {k(tal_n(118), 0)} "
        f"at Z = 118. This is why Og-294 already reaches the tipping point by the liquid-drop formula, x = {k(x294, 2)}; in the valley, Z = 118 would lie at x = "
        f"{k(x_kipp(118, 118 + tal_n(118)), 2)}, with the shell term at {k(xs(118), 2)} (Og-294: {k(KE.x(118, 294, SC), 2)}).\n"
        f"- Masses are measured up to Z = {zg}. The most neutron-rich measured nucleus there has N = {ng}, the valley requires N = {k(tal_n(zg), 0)}.\n"))
    w("---\n")
    w(T("## 4 Neue Elemente\n", "## 4 New Elements\n"))
    kmin, kmax = k(min(v[0] for v in ks.values()), 0), k(max(v[0] for v in ks.values()), 0)
    ktr = k(gr["tropfen"][0], 0)
    w(T("Über 118 hinaus setzt die Hülle die achte Periode fort: erst 8s mit Ununennium (Uue, 119) und Unbinilium (Ubn, 120), dann "
        "nach n + l die erste g-Schale (l = 4), eine Art von Element, die es noch nicht gibt. "
        f"Die Periode hat nach n + l 50 Plätze und schlösse bei {K.schliesszahlen_huelle(8)[-1]}. Zwei Grenzen begrenzen die Elemente, beide gerechnet:\n",
        "Beyond 118 the shell continues the eighth period: first 8s with ununennium (Uue, 119) and unbinilium (Ubn, 120), then "
        "by n + l the first g shell (l = 4), a kind of element that does not yet exist. "
        f"By n + l the period has 50 places and would close at {K.schliesszahlen_huelle(8)[-1]}. Two limits bound the elements, both computed:\n"))
    w(T(f"- **Hülle:** Die Bindung des innersten Elektrons erreicht die doppelte Faltungsenergie des Elektrons, 6c⁴·Δs_e/G, bei Z = "
        f"{k(min(zkr), 0)} bis {k(max(zkr), 0)} (4.3). Das ist robust: Es hängt kaum an der Formel des Kerns.\n"
        f"- **Kern:** Die Kippstelle x = 1 im Tal liegt nach der Tröpfchenformel bei Z = {ktr}, mit dem Schalenglied der Schließzahlen bei "
        f"Z = {kmin} bis {kmax} (je nach nächster Protonen-Schließung 126, 120 oder 114). Der Unterschied kommt aus dem Oberflächenglied "
        f"(a_S = {k(t['aS'], 1)} gegen {k(KE.a_S(SC), 1)} MeV). Die Schalenformel trifft Bindung und Tal der gemessenen Kerne besser (1), Q_α ab Z = 100 schlechter (4.3); "
        f"die Kippstelle selbst ist nicht gemessen (Massen bis Z = {Zd['z_gem_max']}), und beide x ordnen die gemessenen Spalt-Lebensdauern gleich gut.\n",
        f"- **Shell:** The binding of the innermost electron reaches twice the fold energy of the electron, 6c⁴·Δs_e/G, at Z = "
        f"{k(min(zkr), 0)} to {k(max(zkr), 0)} (4.3). This is robust: it hardly depends on the formula of the nucleus.\n"
        f"- **Nucleus:** The tipping point x = 1 in the valley lies at Z = {ktr} by the liquid-drop formula, and at Z = {kmin} to {kmax} "
        f"with the shell term of the closure numbers (depending on the next proton closure 126, 120 or 114). The difference comes from the surface term "
        f"(a_S = {k(t['aS'], 1)} versus {k(KE.a_S(SC), 1)} MeV). The shell formula meets binding and valley of the measured nuclei better (1), Q_α from Z = 100 worse (4.3); "
        f"the tipping point itself is not measured (masses up to Z = {Zd['z_gem_max']}), and both x order the measured fission lifetimes equally well.\n"))
    w(T(f"Gesichert, nach beiden Formeln im Tal gehalten, sind darum die Elemente Uue (119) bis {systemname(zsicher)[1]} ({zsicher}). Nach der "
        f"Formel, die die Daten besser treffen, mit Schalenglied, reichen sie bis zur Grenze der Hülle, {systemname(zmax)[1]} ({zmax}); dann schlösse die achte Periode bei 168. Die Tabellen "
        f"reichen deshalb bis Z = {zmax} und geben x nach beiden Formeln.\n",
        f"Secured, held in the valley by both formulas, are therefore the elements Uue (119) to {systemname(zsicher)[1]} ({zsicher}). By the "
        f"formula that meets the data better, with the shell term, they reach the limit of the shell, {systemname(zmax)[1]} ({zmax}); then the eighth period would close at 168. The tables "
        f"therefore run to Z = {zmax} and give x by both formulas.\n"))
    w(T("Die Namen sind die systematischen nach IUPAC, aus den Ziffern gebildet, bis ein Element nachgewiesen und benannt ist.\n",
        "The names are the systematic IUPAC names, formed from the digits, until an element is discovered and named.\n"))
    w(T("### 4.1 Hülle und Lage im Tal\n", "### 4.1 Shell and Position in the Valley\n"))
    w(T("| Z | Name | Symbol | Konfiguration nach n + l | aus zwei Zählweisen | Grundterm · Moment [μ_B] | Schale nach n + l | Zα | A_Tal | N_Tal | N/Z | x Tropfen | x mit Schale | Neutronenlücke der Verschmelzung |",
        "| Z | Name | Symbol | configuration by n + l | from two countings | ground term · moment [μ_B] | shell by n + l | Zα | A_valley | N_valley | N/Z | x drop | x with shell | neutron gap of fusion |"))
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    wie3 = lambda Z: "–" if Z not in zw else (T("wie n + l", "as n + l") if zw[Z] == {q: v for q, v in nl(Z).items() if v}
                                             else konfig_text(Z, zw[Z]))
    for Z in neu:
        n = tal_n(Z)
        nm, sy = systemname(Z)
        xs = x_schale(Z)
        w(f"| {Z} | {nm if sprache == 'de' else nm.lower()} | {sy} | {konfig_nl(Z)} | {wie3(Z)} | {term_text(SD3.grundterm(zw.get(Z, nl(Z))))} | {letzte_schale(Z)} | {k(Z * K.alpha, 3)} | {k(n + Z, 1)} | "
          f"{k(n, 1)} | {k(n / Z, 3)} | {k(x_kipp(Z, n + Z), 3)} | {k(xs, 3)} | {k(KE.fusionsluecke(Z, TR), 1)} |")
    w("")
    w(T("Konfiguration nach n + l und aus den zwei Zählweisen (Suche ohne Messung wie in 3.3; – noch nicht gerechnet); Grundterm aus der "
        "gerechneten Konfiguration, sonst nach n + l. "
        "A_Tal und N_Tal nach der Tröpfchenformel; das Tal mit Schalenglied liegt bis auf wenige Nukleonen gleich. Die Neutronenlücke ist N_Tal(Z) "
        "minus die größte Neutronenzahl, die zwei Kerne aus dem Tal (je Z ≥ 20) beim Verschmelzen mitbringen.\n",
        "Configuration by n + l and from the two countings (search without measurement as in 3.3; – not yet computed); ground term from the "
        "computed configuration, otherwise by n + l. "
        "A_valley and N_valley by the liquid-drop formula; the valley with the shell term lies the same to within a few nucleons. The neutron gap is "
        "N_valley(Z) minus the largest neutron number that two nuclei from the valley (each Z ≥ 20) bring along when they fuse.\n"))
    pq, ps = KE.pruefung("qa", 100, TR), KE.pruefung("sn", 90, TR)
    w(T("### 4.2 Eigenschaften des Kerns im Tal\n", "### 4.2 Properties of the Nucleus in the Valley\n"))
    w(T("Für den Kern im Tal (A auf ganze Zahl gerundet) gibt die Formel Bindung je Nukleon B/A, Energie des Alpha-Zerfalls Q_α und "
        "Neutronen-Ablösung S_n. Geprüft an gemessenen Kernen (AME2020): Q_α an "
        f"{pq[0]} Zerfällen ab Z = 100 im Mittel {k(pq[1], 2)} MeV daneben, Streuung {k(pq[2], 2)} MeV; S_n an {ps[0]} Kernen ab Z = 90 "
        f"im Mittel {k(ps[1], 2)} MeV, Streuung {k(ps[2], 2)} MeV. Die Lebensdauer lg T_α (4.3) gilt für gerade-gerade Kerne; jenseits der "
        "Kippstelle der jeweiligen Formel ist der Kern nicht gehalten.\n",
        "For the nucleus in the valley (A rounded to an integer) the formula gives the binding per nucleon B/A, the alpha-decay energy Q_α and "
        "the neutron separation energy S_n. Tested on measured nuclei (AME2020): Q_α on "
        f"{pq[0]} decays from Z = 100 off by {k(pq[1], 2)} MeV on average, scatter {k(pq[2], 2)} MeV; S_n on {ps[0]} nuclei from Z = 90 "
        f"off by {k(ps[1], 2)} MeV on average, scatter {k(ps[2], 2)} MeV. The lifetime lg T_α (4.3) holds for even-even nuclei; beyond the "
        "tipping point of the respective formula the nucleus is not held.\n"))
    w(T("| Z | Symbol | A | B/A [keV] | Q_α [MeV] | Q_α mit Schale [MeV] | S_n [MeV] | lg T_α [s] |",
        "| Z | Symbol | A | B/A [keV] | Q_α [MeV] | Q_α with shell [MeV] | S_n [MeV] | lg T_α [s] |"))
    w("|---|---|---|---|---|---|---|---|")
    for Z in neu:
        A = round(Z + tal_n(Z))
        q1, q2 = KE.q_alpha(Z, A, TR), KE.q_alpha(Z, A, SC)
        t1, t2 = sorted((KE.lg_t_alpha(Z, q1), KE.lg_t_alpha(Z, q2)))
        w(f"| {Z} | {systemname(Z)[1]} | {A} | {k(b_je_a(Z, A), 0)} | {k(q1, 2)} | {k(q2, 2)} | {k(KE.s_n(Z, A, TR), 2)} | {k(t1, 1)} {T('bis', 'to')} {k(t2, 1)} |")
    w("")
    g = KE.b_gemessen()
    qa_ds = g[(108, 266)] + g[(2, 4)] - g[(110, 270)]
    sicher = range(119, zsicher + 1)
    qs = [KE.q_alpha(Z, round(Z + tal_n(Z)), TR) for Z in sicher]
    qs2 = [KE.q_alpha(Z, round(Z + tal_n(Z)), SC) for Z in sicher]
    sn = [KE.s_n(Z, round(Z + tal_n(Z)), TR) for Z in neu]
    w(T(f"Im gesicherten Bereich bis Z = {zsicher} liegt Q_α bei {k(min(qs), 1)} bis {k(max(qs), 1)} MeV, mit Schalenglied bei {k(min(qs2), 1)} bis "
        f"{k(max(qs2), 1)} MeV; der schwerste gemessene Kern Ds-270 hat {k(qa_ds, 2)} MeV. S_n bleibt bis Z = {zmax} mit {k(min(sn), 1)} bis "
        f"{k(max(sn), 1)} MeV positiv: Die Kerne im Tal halten ihre Neutronen. Das Wechseln von S_n folgt der Paarung.\n",
        f"In the secured range up to Z = {zsicher}, Q_α lies at {k(min(qs), 1)} to {k(max(qs), 1)} MeV, with the shell term at {k(min(qs2), 1)} to "
        f"{k(max(qs2), 1)} MeV; the heaviest measured nucleus Ds-270 has {k(qa_ds, 2)} MeV. S_n stays positive up to Z = {zmax}, at {k(min(sn), 1)} to "
        f"{k(max(sn), 1)} MeV: the nuclei in the valley hold their neutrons. The alternation of S_n follows the pairing.\n"))
    w(T("### 4.3 Schalen, Relativität, Lebensdauer\n", "### 4.3 Shells, Relativity, Lifetime\n"))
    o0, o1, nk = KE.anpassung(TR)[1], KE.anpassung(SC)[1], KE.anpassung(SC)[2]
    sq, ss = KE.pruefung("qa", 100, SC), KE.pruefung("sn", 90, SC)
    sq84, pq84 = KE.pruefung("qa", 84, SC), KE.pruefung("qa", 84, TR)
    sch = {p: KE.streuung_ab(("schale", p)) for p in (114, 120, 126)}
    nlo, nhi = tal_n(119), tal_n(zsicher)
    w(T(f"**Schalenkorrektur im Kern.** Die Tröpfchenformel kennt keine Schließungen. Setzt man die Schließzahlen des Modells "
        f"(2, 8, 20, 28, 50, 82, 126, 184, 258) als Schalenglied ein, in der Form von Myers und Swiatecki, mit Wigner- und "
        f"Oberflächen-Asymmetrieglied, sinkt die Streuung der Bindung an {nk} gemessenen Kernen von {k(o0, 2)} auf {k(o1, 2)} MeV. "
        f"S_n wird besser ({k(ps[2], 2)} → {k(ss[2], 2)} MeV ab Z = 90), Q_α ab Z = 84 auch ({k(pq84[2], 2)} → {k(sq84[2], 2)} MeV), "
        f"ab Z = 100 aber nicht ({k(pq[2], 2)} → {k(sq[2], 2)} MeV): Die Verformung der Kerne zwischen den Schalen fehlt der Form. "
        "Darum stehen in 4.2 beide Werte; der Unterschied ist die Unsicherheit. Ob die nächste Protonen-Schließung bei 114, 120 oder "
        f"126 liegt, entscheiden die gemessenen Kerne nicht: Streuung ab Z = 90 {k(sch[114], 2)}, {k(sch[120], 2)} und {k(sch[126], 2)} MeV. "
        f"Die Kerne im Tal von Uue bis {systemname(zsicher)[1]} liegen mit N ≈ {k(nlo, 0)}–{k(nhi, 0)} zwischen den Neutronen-Schließungen 184 und 258.\n",
        f"**Shell correction in the nucleus.** The liquid-drop formula knows no closures. Inserting the closure numbers of the model "
        f"(2, 8, 20, 28, 50, 82, 126, 184, 258) as a shell term, in the form of Myers and Swiatecki, with Wigner and "
        f"surface-asymmetry terms, lowers the scatter of the binding on {nk} measured nuclei from {k(o0, 2)} to {k(o1, 2)} MeV. "
        f"S_n improves ({k(ps[2], 2)} → {k(ss[2], 2)} MeV from Z = 90), Q_α from Z = 84 too ({k(pq84[2], 2)} → {k(sq84[2], 2)} MeV), "
        f"but not from Z = 100 ({k(pq[2], 2)} → {k(sq[2], 2)} MeV): the form lacks the deformation of the nuclei between the shells. "
        "This is why 4.2 gives both values; their difference is the uncertainty. Whether the next proton closure lies at 114, 120 or "
        f"126 is not decided by the measured nuclei: scatter from Z = 90 {k(sch[114], 2)}, {k(sch[120], 2)} and {k(sch[126], 2)} MeV. "
        f"The valley nuclei from Uue to {systemname(zsicher)[1]}, with N ≈ {k(nlo, 0)}–{k(nhi, 0)}, lie between the neutron closures 184 and 258.\n"))
    w(T(f"**Relativität der Hülle.** Nach Dirac ist die Energie des innersten Elektrons um einen punktförmigen Kern (3c⁴·Δs_e/G)·√(1 − (Zα)²), "
        "mit der Faltung Δs_e des Elektrons; bei Zα = 1, Z = 137, endet sie. Um den ausgedehnten Kern mit R = λ_H·A^{1/3} bleibt der Zustand "
        f"gebunden und erreicht die Bindung 6c⁴·Δs_e/G, die doppelte Faltungsenergie des Elektrons, erst bei Z = {k(min(zkr), 1)} bis {k(max(zkr), 1)} "
        "(Dirac-Gleichung für eine gleichmäßig geladene Kugel, auf zwei Wegen gerechnet; die Spanne kommt aus der Halterlänge und dem Tal ohne "
        "und mit Schalenglied).\n",
        f"**Relativity of the shell.** By Dirac, the energy of the innermost electron around a point nucleus is (3c⁴·Δs_e/G)·√(1 − (Zα)²), "
        "with the fold Δs_e of the electron; at Zα = 1, Z = 137, it ends. Around the extended nucleus with R = λ_H·A^{1/3} the state stays bound "
        f"and reaches the binding 6c⁴·Δs_e/G, twice the fold energy of the electron, only at Z = {k(min(zkr), 1)} to {k(max(zkr), 1)} "
        "(Dirac equation for a uniformly charged sphere, computed by two routes; the range comes from the holder length and the valley without "
        "and with the shell term).\n"))
    (la, lb, lc), lst, ln, lu = KE.alpha_lebensdauer()
    sfn, sfr, sfs, sfx = KE.spaltung_gegen_x(TR)
    paar = [(KE.lg_t_alpha(Z, KE.q_alpha(Z, round(Z + tal_n(Z)), TR)), KE.lg_t_alpha(Z, KE.q_alpha(Z, round(Z + tal_n(Z)), SC))) for Z in sicher]
    t_lo, t_hi = min(min(p) for p in paar), max(max(p) for p in paar)
    w(T(f"**Lebensdauer.** Der Alpha-Zerfall ist ein Durchgang durch die Coulomb-Schwelle, der Tunnelfaktor e^{{−√(E_G/E)}} des Modells. "
        f"Daraus folgt lg T_α = a·Z/√Q_α + b·Z + c. An {ln} gerade-geraden Alpha-Strahlern ab Z = 84 (NUBASE2020, AME2020) ist "
        f"a = {k(la, 3)}, b = {k(lb, 4)}, c = {k(lc, 2)} (T in s), Streuung {k(lst, 2)} Dekaden; U-238 trifft sie auf {k(lu, 2)} Dekaden. "
        "Mit den beiden Q_α aus 4.2 gibt sie die Spanne lg T_α der Tabelle, für gerade-gerade Kerne; ungerade leben länger. "
        f"Im gesicherten Bereich, Uue bis {systemname(zsicher)[1]}, hielten die Kerne im Tal gegen den Alpha-Zerfall 10^{{{k(t_lo, 1)}}} bis "
        f"10^{{{k(t_hi, 1)}}} s, {dauer(t_lo, sprache)} bis {dauer(t_hi, sprache)}; je Element steht die Spanne in 4.2. "
        f"Gegen Spaltung gibt es keine solche Form: An {sfn} gemessenen Kernen fällt die Spalt-Lebensdauer zwar mit x (Korrelation "
        f"{k(sfr, 2)}), streut aber um {k(sfs, 1)} Dekaden, und gemessen ist nur bis x = {k(sfx, 3)} (Tröpfchenformel). Wie lange die neuen "
        "Elemente gegen Spaltung halten, sagen Modell und Daten nicht voraus; nahe x = 1 verschwindet die Rückstellung.\n",
        f"**Lifetime.** Alpha decay is a passage through the Coulomb barrier, the tunnelling factor e^{{−√(E_G/E)}} of the model. "
        f"From it follows lg T_α = a·Z/√Q_α + b·Z + c. On {ln} even-even alpha emitters from Z = 84 (NUBASE2020, AME2020), "
        f"a = {k(la, 3)}, b = {k(lb, 4)}, c = {k(lc, 2)} (T in s), scatter {k(lst, 2)} decades; it meets U-238 to {k(lu, 2)} decades. "
        "With the two Q_α of 4.2 it gives the range lg T_α of the table, for even-even nuclei; odd ones live longer. "
        f"In the secured range, Uue to {systemname(zsicher)[1]}, the valley nuclei would hold against alpha decay for 10^{{{k(t_lo, 1)}}} to "
        f"10^{{{k(t_hi, 1)}}} s, {dauer(t_lo, sprache)} to {dauer(t_hi, sprache)}; the range per element is given in 4.2. "
        f"Against fission there is no such form: on {sfn} measured nuclei the fission lifetime does fall with x (correlation "
        f"{k(sfr, 2)}), but scatters by {k(sfs, 1)} decades, and measurements reach only x = {k(sfx, 3)} (liquid drop). How long the new "
        "elements hold against fission is predicted neither by the model nor by the data; near x = 1 the restoring vanishes.\n"))
    w(T("### 4.4 Wie, wo, wie lange\n", "### 4.4 How, Where, How Long\n"))
    wtr, wsc = KE.entstehungsweg(TR), KE.entstehungsweg(SC)
    w(T("**Wie sie entstehen können.** Das Tal krümmt sich: N_Tal wächst schneller als Z. Zwei Kerne aus dem Tal bringen beim Verschmelzen "
        f"darum immer zu wenige Neutronen mit, für Z = 119 bis {zmax} fehlen nach der Tabelle {k(KE.fusionsluecke(119, TR), 0)} bis {k(KE.fusionsluecke(zmax, TR), 0)}; "
        "Kerne, die aus Verschmelzen stammen, liegen neutronenarm und kippen früher, wie Og-294. Ins Tal führt nur das schrittweise Anlagern "
        f"von Neutronen an schwere Kerne, bis N ≈ {k(tal_n(119), 0)} bei Uue (Z = 119), N ≈ {k(tal_n(zsicher), 0)} bei "
        f"{systemname(zsicher)[1]} (Z = {zsicher}) und N ≈ {k(tal_n(zmax), 0)} bei {systemname(zmax)[1]} (Z = {zmax}) erreicht ist, "
        "mit Umwandlungen n → p, die Z erhöhen. "
        f"Neutronenreiche Kerne halten besser, x = a_C·Z²/(2a_S·A) fällt mit A: Zwischen Tal und Abbruchkante (das nächste Neutronenpaar "
        f"nicht mehr gebunden) bleibt der Kern nach der Tröpfchenformel gehalten bis Z = {wtr[0]} (Kante bei N = {wtr[1][wtr[0]][1]}, "
        f"x = {k(wtr[1][wtr[0]][3], 3)}), mit Schalenglied über Z = {max(wsc[1])} hinaus; das Tal selbst nur bis {max(Z for Z, v in wtr[1].items() if v[2] < 1)} "
        f"bzw. {max(Z for Z, v in wsc[1].items() if v[2] < 1)}. Auf neutronenreichem Weg setzt also die Hülle (Z = {k(KE.huellen_grenze(TR), 1)}) das Ende, "
        "nicht der Kern. Die Kanten liegen weit jenseits der gemessenen Massen; sie sind Fortsetzung der angepassten Formeln.\n",
        "**How they can form.** The valley curves: N_valley grows faster than Z. Two nuclei from the valley therefore always bring too few "
        f"neutrons when they fuse; for Z = 119 to {zmax} the table shows {k(KE.fusionsluecke(119, TR), 0)} to {k(KE.fusionsluecke(zmax, TR), 0)} missing. "
        "Nuclei made by fusion lie neutron-poor and tip earlier, like Og-294. Only the stepwise addition of neutrons to heavy nuclei leads into "
        f"the valley, until N ≈ {k(tal_n(119), 0)} is reached for Uue (Z = 119), N ≈ {k(tal_n(zsicher), 0)} for "
        f"{systemname(zsicher)[1]} (Z = {zsicher}) and N ≈ {k(tal_n(zmax), 0)} for {systemname(zmax)[1]} (Z = {zmax}), "
        "with conversions n → p that raise Z. "
        f"Neutron-rich nuclei hold better, x = a_C·Z²/(2a_S·A) falls with A: between the valley and the drip edge (the next neutron pair "
        f"no longer bound) the nucleus stays held by the liquid-drop formula up to Z = {wtr[0]} (edge at N = {wtr[1][wtr[0]][1]}, "
        f"x = {k(wtr[1][wtr[0]][3], 3)}), with the shell term beyond Z = {max(wsc[1])}; the valley itself only up to {max(Z for Z, v in wtr[1].items() if v[2] < 1)} "
        f"and {max(Z for Z, v in wsc[1].items() if v[2] < 1)}. On a neutron-rich path the shell (Z = {k(KE.huellen_grenze(TR), 1)}) thus sets the end, "
        "not the nucleus. The edges lie far beyond the measured masses; they are a continuation of the fitted formulas.\n"))
    zb, ab_ = KE.ende_brennleiter(TR)
    zs_, as_ = KE.ende_brennleiter(SC)
    bk = KE.bekannte()
    kt_co = K.decke_entartet(0.875)[0]
    w(T("**Wo die Elemente entstehen.** Jeder Halter des Trichters macht seinen Teil. *Wärme:* Im Stern brennen die Kerne als Leiter in der Tiefe, "
        f"solange die Bindung je Nukleon wächst, bis x = ¼; im Tal ist das nach der Tröpfchenformel Z = {k(zb, 1)}, A = {k(ab_, 1)}, zwischen "
        f"{sym[bk['m_min'][0]]}-{bk['m_min'][1]} (kleinste Masse je Nukleon) und {sym[bk['b_max'][0]]}-{bk['b_max'][1]} (größte Bindung je Nukleon, AME2020); "
        f"mit Schalenglied Z = {k(zs_, 1)}. Bis zur Eisengruppe entstehen die Elemente durch Brennen. *Zählregel der Elektronen:* Ein Trichter, dessen Tiefe "
        f"sie hält, erreicht höchstens kT = {k(K.decke_entartet(1.0)[0] / (K.m_e * K.c ** 2), 4)}·(μ/μ_e)·m_ec², für Kohlenstoff und Sauerstoff "
        f"{k(kt_co / K.e / 1e3, 0)} keV ({k(kt_co / K.k_B / 1e8, 1)}·10⁸ K), an der Grenzmasse "
        f"{k(K.grenzmasse() / K.M_sun, 3)} M☉. *Zählregel der Nukleonen:* Über die Eisengruppe hinaus führt nur der Einfang von Neutronen (Wie, Wo).\n",
        "**Where the elements form.** Each holder of the funnel does its part. *Heat:* in the star the nuclei burn as a ladder in the depth "
        f"as long as the binding per nucleon grows, up to x = ¼; in the valley this is Z = {k(zb, 1)}, A = {k(ab_, 1)} by the liquid-drop formula, between "
        f"{sym[bk['m_min'][0]]}-{bk['m_min'][1]} (smallest mass per nucleon) and {sym[bk['b_max'][0]]}-{bk['b_max'][1]} (largest binding per nucleon, AME2020); "
        f"with the shell term Z = {k(zs_, 1)}. Up to the iron group the elements form by burning. *Counting rule of the electrons:* a funnel whose depth "
        f"it holds reaches at most kT = {k(K.decke_entartet(1.0)[0] / (K.m_e * K.c ** 2), 4)}·(μ/μ_e)·m_ec², for carbon and oxygen "
        f"{k(kt_co / K.e / 1e3, 0)} keV ({k(kt_co / K.k_B / 1e8, 1)}·10⁸ K), at the limiting mass "
        f"{k(K.grenzmasse() / K.M_sun, 3)} M☉. *Counting rule of the nucleons:* beyond the iron group only neutron capture leads on (How, Where).\n"))
    ms = SD3.makro_schwellen()                                                              # wo die Ionisierung gemessen ist
    koerper = [(T("Erde", "Earth"), K.G * K.M_earth / (6.371e6 * K.c ** 2)), (T("Sonne", "Sun"), K.G * K.M_sun / (6.957e8 * K.c ** 2)),
               (T("Neutronenstern 1,4 M☉, 12 km", "neutron star 1.4 M☉, 12 km"), K.G * 1.4 * K.M_sun / (12e3 * K.c ** 2))]
    dichte = {int(x["AtomicNumber"]): float(x["Density"]) * 1000 for x in rows if x["Density"] and x["StandardState"] != "Gas"}
    mkrit = {Z: (ms[Z]["d_P"] * K.c ** 2 / K.G) ** 1.5 * (3 / (4 * math.pi * dichte[Z])) ** 0.5 for Z in dichte}
    f3 = K.polytrope(3.0)[2]
    hm = lambda d: min(ms, key=lambda Z: abs(math.log(ms[Z]["d_P"] / d)))
    w(T("**Makrozustand: wo die Hülle hält.** Die Atome falten als Masse den Körper, δ = GM/(Rc²) = 3Δs/R; der Körper gibt jedem Atom "
        "Druck, P·V_A = ½·m_A·c²·δ (Zentrum der homogenen Kugel), und als wärmegehaltener Trichter Wärme, kT = f·δ·m̄c² "
        f"(f = {k(f3, 4)} aus der Polytrope n = 3). Die Hülle hält, solange beides unter der Bindung des äußeren Elektrons IE bleibt (PubChem): "
        "bis δ_P = 2·IE/(m_A·c²) und δ_T = IE/(f·m_A·c²). Gegenüber stehen die Körper: "
        + ", ".join(f"{n} δ = {k(d * 10 ** -math.floor(math.log10(d)), 2)}·10^{math.floor(math.log10(d))}" for n, d in koerper)
        + f". Die Schwellen aller Elemente liegen zwischen {k(min(v['d_P'] for v in ms.values()) * 1e11, 1)}·10⁻¹¹ und "
        f"{k(max(v['d_P'] for v in ms.values()) * 1e8, 1)}·10⁻⁸: Wo die Elemente entstehen, im Stern und im Neutronenstern, ist jede Hülle offen; "
        "sie entstehen als Kerne, und ihre Hülle, ihre Konfiguration und ihre Chemie bilden sich erst in flachen Bereichen, im Grenzfall δ → 0 des "
        "Gases, wo die Konfigurationen gemessen sind. Umgekehrt setzt ein Körper aus dem Element mit seiner gemessenen Dichte ρ "
        "in seinem Zentrum δ_P ab der Masse M_krit = (δ_P·c²/G)^{3/2}·(3/(4πρ))^{1/2} (Elemente, die fest oder flüssig sind). Die erste Spalte "
        "nennt den Bereich des Halters (Brennleiter bis x = ¼, darüber Einfang von Neutronen), nicht den Weg jedes Kerns. "
        "Für alle Elemente mit gemessener Ionisierungsenergie.\n",
        "**Macro state: where the shell holds.** The atoms as mass fold the body, δ = GM/(Rc²) = 3Δs/R; the body gives every atom "
        "pressure, P·V_A = ½·m_A·c²·δ (centre of the homogeneous sphere), and as a heat-held funnel heat, kT = f·δ·m̄c² "
        f"(f = {k(f3, 4)} from the polytrope n = 3). The shell holds as long as both stay below the binding of the outer electron IE (PubChem): "
        "up to δ_P = 2·IE/(m_A·c²) and δ_T = IE/(f·m_A·c²). Against them stand the bodies: "
        + ", ".join(f"{n} δ = {k(d * 10 ** -math.floor(math.log10(d)), 2)}·10^{math.floor(math.log10(d))}" for n, d in koerper)
        + f". The thresholds of all elements lie between {k(min(v['d_P'] for v in ms.values()) * 1e11, 1)}·10⁻¹¹ and "
        f"{k(max(v['d_P'] for v in ms.values()) * 1e8, 1)}·10⁻⁸: where the elements form, in the star and the neutron star, every shell is open; "
        "they form as nuclei, and their shell, configuration and chemistry form only in flat regions, in the limit δ → 0 of the gas, "
        "where the configurations are measured. Conversely a body of the element with its measured density ρ reaches δ_P "
        "in its centre from the mass M_krit = (δ_P·c²/G)^{3/2}·(3/(4πρ))^{1/2} (elements that are solid or liquid). The first column "
        "names the range of the holder (burning ladder up to x = ¼, beyond it neutron capture), not the path of each nucleus. "
        "For all elements with a measured ionisation energy.\n"))
    w(T("| Z | | Bereich des Halters | IE [eV] | δ_P | δ_T | M_krit [Erdmassen] | Hülle in der Erde (Zentrum) |",
        "| Z | | range of the holder | IE [eV] | δ_P | δ_T | M_krit [Earth masses] | shell in the Earth (centre) |"))
    w("|---|---|---|---|---|---|---|---|")
    ze = KE.ende_brennleiter(TR)[0]
    for Z in range(1, 171):
        if Z not in ms:
            w(f"| {Z} | {systemname(Z)[1]} | {T('Neutroneneinfang', 'neutron capture')} | – | – | – | – | – |")
            continue
        v = ms[Z]
        ent = T("Ausgang der Leiter", "start of the ladder") if Z == 1 else (T("Brennen", "burning") if Z <= ze else T("Neutroneneinfang", "neutron capture"))
        e10 = lambda d: f"{k(d * 10 ** -math.floor(math.log10(d)), 2)}·10^{math.floor(math.log10(d))}"
        w(f"| {Z} | {sym[Z] if Z in sym else systemname(Z)[1]} | {ent} | {k(v['ie'], 2)} | {e10(v['d_P'])} | {e10(v['d_T'])} | "
          f"{k(mkrit[Z] / K.M_earth, 3) if Z in mkrit else '–'} | {T('hält', 'holds') if v['d_P'] > koerper[0][1] else T('offen', 'open')} |")
    w("")
    w(T("**Wo.** Dafür braucht es freie Neutronen in großer Dichte. Im Modell gibt es sie dort, wo die Zählregel der Nukleonen selbst die Materie "
        f"hält, der gehaltene Trichter mit Spreizung von der Ordnung eins. Ob und wie schnell Kerne dort bis A ≈ {k(at, 0)} und darüber wachsen, "
        "rechnet das Modell nicht.\n",
        "**Where.** This requires free neutrons at high density. In the model they exist where the counting rule of the nucleons itself holds the "
        f"matter, the held funnel with a spread of order one. Whether and how fast nuclei grow there up to A ≈ {k(at, 0)} and beyond the model "
        "does not compute.\n"))
    w(T("**Wie lange.** Die Rückstellung gegen Verformung ist ∝ 1 − x. Im Tal ist sie nach der Tröpfchenformel bei Uue "
        f"{k(1 - x_kipp(119, 119 + tal_n(119)), 3)}, bei {systemname(zsicher)[1]} {k(1 - x_kipp(zsicher, zsicher + tal_n(zsicher)), 3)}, an der Kippstelle null; mit Schalenglied bei Uue "
        f"{k(1 - x_schale(119), 3)}. Zum Vergleich U-238: {k(1 - x_kipp(92, 238), 3)} nach der Tröpfchenformel, "
        f"{k(1 - KE.x(92, 238, SC), 3)} mit Schalenglied.\n",
        "**How long.** The restoring against deformation is ∝ 1 − x. In the valley, by the liquid-drop formula it is "
        f"{k(1 - x_kipp(119, 119 + tal_n(119)), 3)} for Uue and {k(1 - x_kipp(zsicher, zsicher + tal_n(zsicher)), 3)} for {systemname(zsicher)[1]}, zero at the tipping point; with the shell term "
        f"{k(1 - x_schale(119), 3)} for Uue. For comparison, U-238: {k(1 - x_kipp(92, 238), 3)} by the liquid-drop formula, "
        f"{k(1 - KE.x(92, 238, SC), 3)} with the shell term.\n"))
    x126t, x126s = x_kipp(126, 310), KE.x(126, 310, SC)
    w(T(f"**Was offen bleibt und was nicht möglich ist.** Über Z = {zmax} endet die Hülle. Zwischen {zsicher} und {zmax} hängt alles am Schalenglied. "
        f"Die Doppelschließung Unbihexium (Ubh, Z = 126) mit N = 184 ist neutronenarm: x = {k(x126t, 2)} nach der Tröpfchenformel, "
        f"{k(x126s, 2)} mit Schalenglied.\n",
        f"**What remains open and what is not possible.** Above Z = {zmax} the shell ends. Between {zsicher} and {zmax} everything depends on the shell term. "
        f"The double closure unbihexium (Ubh, Z = 126) with N = 184 is neutron-poor: x = {k(x126t, 2)} by the liquid-drop formula, "
        f"{k(x126s, 2)} with the shell term.\n"))
    w(T("**Grenze der Rechnung.** Die Kerne über A = 295 sind Fortsetzung der angepassten Formeln, keine Messung.\n",
        "**Limit of the calculation.** The nuclei above A = 295 are a continuation of the fitted formulas, not a measurement.\n"))
    return "\n".join(o)


def schreibe():
    for sp in ("de", "en"):
        with open(ZIEL[sp], "w", encoding="utf-8") as f:
            f.write(text(sp))


if __name__ == "__main__":
    schreibe()
    print(ZIEL)
