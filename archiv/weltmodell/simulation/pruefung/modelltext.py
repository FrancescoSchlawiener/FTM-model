# -*- coding: utf-8 -*-
"""
Werkzeug für den Modelltext (modell/*.md): Dateien je Baustein, Sätze, Zahlen und ihre Stellen.

Zahl im Modelltext: deutsches Komma, Tausender mit Leerzeichen, Zehnerpotenz mit ·10ⁿ oder allein 10ⁿ (auch mitten im
Satz, etwa „rund 10⁻²⁹“), Vorzeichen −. Kleine ganze Zahlen ≤ 12 ohne Komma und ohne Potenz gelten nicht als Zahl
(Nummern, Exponenten, Zählwörter).
Geltende Stellen: Nullen am Ende einer ganzen Zahl ohne Komma gelten nicht als Stelle (2780 → halbe Stelle 5, 1000 →
500); die halbe letzte Stelle ist die Toleranz der Regel.
Fundstelle: jede Zahl an ihrem Ort (Anfang, Ende im Text); stellen() gibt sie, satz_spannen() die Sätze als Spannen.
Satz: Zeile oder Abschnitt bis zu einem Punkt, dem ein Großbuchstabe folgt; Abkürzungen (z. B., d. h., Mrd. …)
und Tabellenzeilen trennen nicht.
"""
import glob
import os
import re

REGISTERMODUL = True
MODELL = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
                      "modell")
_HOCH = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
# erste Form: Mantisse mit optionalem ·10ⁿ; sie endet nicht vor einer Ziffer oder Hochzahl, sonst ist es ein bloßes 10ⁿ
# (zweite Form). Ohne diese Sperre läse „10⁻⁶“ als „10“ und fiele als kleine Zahl weg.
ZAHL = re.compile(r"(?<![\w,·^.])([−-]?)(\d{1,3}(?:[  ]\d{3})+|\d+(?:,\d+)?)(?:\s?·\s?10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+))?"
                  r"(?![\d⁰¹²³⁴⁵⁶⁷⁸⁹⁻])"
                  r"|(?<![\w·])10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)")
_ABK = ("z. B.", "d. h.", "u. a.", "z. T.", "bzw.", "Mrd.", "Mio.", "ca.", "vgl.", "Nr.", "Abb.", "usw.", "etc.", "Z.",
        "S.", "Kap.", "Gl.")


def datei(baustein):
    """Dateiname in modell/ zum Baustein "01" … "09"."""
    treffer = glob.glob(os.path.join(MODELL, f"{baustein}-*.md"))
    if len(treffer) != 1:
        raise FileNotFoundError(f"Baustein {baustein}: {treffer}")
    return os.path.basename(treffer[0])


def lies(dateiname):
    with open(os.path.join(MODELL, dateiname), encoding="utf-8") as f:
        return f.read()


def saetze(text):
    """Liste der Sätze (siehe Moduldoc)."""
    geschuetzt = text
    for i, a in enumerate(_ABK):
        geschuetzt = geschuetzt.replace(a, a.replace(".", f"\x00{i}\x00"))
    teile = []
    for zeile in geschuetzt.split("\n"):
        teile += re.split(r"(?<=[.!?])\s+(?=[A-ZÄÖÜ„(*])", zeile)
    return [re.sub(r"\x00\d+\x00", ".", t) for t in teile if t.strip()]


def zahl_aus_text(s):
    """(Wert, halbe letzte Stelle) einer Zahl, wie sie im Modelltext steht; halbe Stelle None bei bloßem 10ⁿ."""
    m = ZAHL.fullmatch(s.strip())
    if m is None:
        raise ValueError(f"keine Zahl: {s!r}")
    if m.group(4):
        return 10.0 ** int(m.group(4).translate(_HOCH)), None
    mant = m.group(2).replace(" ", "").replace(" ", "")
    if "," in mant:
        stelle = -len(mant.split(",")[1])
    else:
        stelle = len(mant) - len(mant.rstrip("0")) if mant.strip("0") else 0   # Nullen am Ende gelten nicht
    potenz = 10.0 ** int(m.group(3).translate(_HOCH)) if m.group(3) else 1.0
    vz = -1.0 if m.group(1) else 1.0
    return vz * float(mant.replace(",", ".")) * potenz, 0.5 * 10.0 ** stelle * potenz


def stellen(text):
    """Alle Zahlen im Text als Fundstellen (Anfang, Ende, normierte Zeichenfolge), in Reihenfolge; der Anfang liegt
    hinter einem Vorzeichen (die Stelle von „−0,13“ ist „0,13“)."""
    aus = []
    for m in ZAHL.finditer(text):
        if m.group(4):
            aus.append((m.start(), m.end(), "10" + m.group(4)))
            continue
        mant = m.group(2)
        if not m.group(3) and "," not in mant and float(mant.replace(" ", "").replace(" ", "")) <= 12:
            continue
        aus.append((m.start(2), m.end(), norm(m.group(0))))
    return aus


def zahlen(text):
    """Alle Zahlen im Text als normierte Zeichenfolgen (ohne Vorzeichen, ohne Leerzeichen), in Reihenfolge."""
    return [z for _, _, z in stellen(text)]


def satz_spannen(text):
    """Die Sätze (siehe saetze) als Spannen (Anfang, Ende) im Text."""
    aus, pos = [], 0
    for s in saetze(text):
        i = text.find(s, pos)
        aus.append((i, i + len(s)))
        pos = i + len(s)
    return aus


def norm(s):
    return s.replace(" ", "").replace(" ", "").replace("−", "").replace("-", "").strip()
