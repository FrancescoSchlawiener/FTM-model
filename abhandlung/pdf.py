"""Die Abhandlung und die Beispiele als PDF, aus <name>_<sprache>.md: ohne Semikolon und Gedankenstrich (Komma, bei
Zahlenbereichen „bis“), Formeln zeilenweise gesetzt, Hoch- und Tiefstellungen gesetzt, Seitenzahlen; Bilder (![Text](bilder/…))
mit Bildunterschrift, nur beschriebene Bilder (Zitate) fallen weg. Gesetzt wird mit dem Chromium des Rechners.
Lauf: python abhandlung/pdf.py [de|en] [abhandlung|rechnen]; schreibt <name>_<sprache>.pdf daneben."""
import glob
import html
import os
import pathlib
import re
import subprocess
import sys
import tempfile

HIER = pathlib.Path(__file__).resolve().parent
BIS = {"de": " bis ", "en": " to "}
UND = {"de": " und ", "en": " and "}
BILD = {"de": "Bild", "en": "Figure"}
HOCH = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ⁿⁱᵐ", "0123456789+−nim")
TIEF = str.maketrans("₀₁₂₃₄₅₆₇₈₉₊₋ᵢ", "0123456789+−i")

STIL = """
@page { size: A4; margin: 24mm 22mm 26mm 22mm;
        @bottom-center { content: counter(page); font-family: "FreeSerif"; font-size: 9.5pt; color: #444; } }
@page titel { @bottom-center { content: none; } }
html { font-family: "FreeSerif", "DejaVu Serif", serif; font-size: 11pt; line-height: 1.42; color: #111; }
body { margin: 0; text-align: justify; hyphens: auto; }
.titelseite { page: titel; height: 230mm; display: flex; flex-direction: column; justify-content: center; text-align: center;
              break-after: page; }
.titelseite h1 { font-size: 30pt; font-weight: bold; margin: 0 0 16pt; }
.titelseite p { font-size: 14pt; margin: 0 10mm; line-height: 1.4; }
.titelseite p.autor { font-size: 13pt; font-variant: small-caps; letter-spacing: 0.04em; margin-top: 28pt; }
p.autor { text-align: center; font-variant: small-caps; letter-spacing: 0.04em; margin: -12pt 0 22pt; }
h1.teil { break-after: avoid; font-size: 21pt; font-weight: bold; margin: 34pt 0 14pt; text-align: left;
          border-top: 0.7pt solid #000; padding-top: 10pt; }
h2 { font-size: 13.5pt; font-weight: bold; margin: 18pt 0 6pt; break-after: avoid; text-align: left; }
h1.teil + h2 { margin-top: 6pt; }
h2.neu { break-before: page; margin-top: 0; }
p, li { orphans: 3; widows: 3; }
p { margin: 0 0 7pt; }
p.vor { break-after: avoid; }
p.motto { font-style: italic; text-align: center; margin: 18pt 10mm 0; }
h1.kopf { font-size: 24pt; font-weight: bold; margin: 4pt 0 8pt; text-align: center; }
p.unter { font-size: 13pt; text-align: center; margin: 0 8mm 22pt; line-height: 1.4; }
ul, ol { margin: 0 0 7pt; padding-left: 15pt; }
ol { padding-left: 18pt; }
li { margin: 0 0 3pt; }
figure { margin: 10pt 0 12pt; text-align: center; break-inside: avoid; }
figure img { width: 122mm; }
figcaption { font-size: 9.6pt; line-height: 1.3; text-align: left; margin: 4pt 9mm 0; }
table.k { break-inside: avoid; }
table.f { margin: 5pt auto 9pt; border-collapse: collapse; break-inside: avoid; break-before: avoid; }
table.f td { padding: 1.5pt 9pt; white-space: nowrap; text-align: left; vertical-align: baseline; }
table.f td.r { padding-top: 7pt; }
table.t { width: 100%; border-collapse: collapse; font-size: 9.6pt; line-height: 1.3; margin: 6pt 0 10pt;
          text-align: left; border: 0.8pt solid #000; }
table.t thead { display: table-header-group; }
table.t tr { break-inside: avoid; }
table.t th { font-weight: bold; vertical-align: bottom; padding: 3.5pt 6pt; border: 0.5pt solid #000; background: #efefef; }
table.t td { vertical-align: top; padding: 3pt 6pt; border: 0.5pt solid #000; }
table.t td:first-child, table.t th:first-child { white-space: nowrap; width: 1%; }
sup, sub { font-size: 70%; line-height: 0; }
sup sup, sup sub, sub sup, sub sub { font-size: 85%; }
sup { vertical-align: 0.55em; }
sub { vertical-align: -0.25em; }
"""


def zeichen(t, sprache):
    """Semikolon zu Komma. Gedankenstrich: zwischen Zahlen „bis“, zwischen zwei Wörtern „und“, sonst Komma."""
    t = t.replace("; ", ", ").replace(";", ",")
    t = re.sub(r"(?<=[0-9%])\s?[–—]\s?(?=[0-9])", BIS[sprache], t)
    t = re.sub(r"(?<=\w)[–—](?=\w)", UND[sprache], t)
    t = re.sub(r"\s*[–—]\s*", ", ", t)
    return t


def hochtief(t):
    """^{…}, ^x, _{…}, _x und die hoch- und tiefgestellten Zeichen des Unicode als <sup> und <sub>."""
    t = re.sub(r"\^\{([^{}]*)\}", r"<sup>\1</sup>", t)
    t = re.sub(r"\^([A-Za-z0-9∞α-ω])", r"<sup>\1</sup>", t)
    t = re.sub(r"_\{([^{}]*)\}", r"<sub>\1</sub>", t)
    t = re.sub(r"_([A-Za-z0-9α-ωϑ-ϖΑ-Ω∞*ᵢ]+)", r"<sub>\1</sub>", t)
    t = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ⁿⁱᵐ]+", lambda m: "<sup>" + m.group(0).translate(HOCH) + "</sup>", t)
    t = re.sub(r"[₀₁₂₃₄₅₆₇₈₉₊₋ᵢ]+", lambda m: "<sub>" + m.group(0).translate(TIEF) + "</sub>", t)
    return t


OPERATOREN = "= − + ± → ↔ ≈ ∝ ≤ ≥ ≪ ≫ ⇒ ∗ &lt; &gt;".split()


def _trennung(t, sprache):
    """Weiche Trennstellen in Wörtern des Fließtexts (nur Buchstaben, ab sieben Zeichen); Formelzeichen bleiben unberührt."""
    try:
        import pyphen
    except ImportError:
        return t
    w = _trennung.woerter.setdefault(sprache, pyphen.Pyphen(lang={"de": "de_DE", "en": "en_US"}[sprache]))

    def trenne(m):                                           # je Seite drei Buchstaben, keine Silbe aus einem Buchstaben, kein „c-k“,
        wort, stellen, vorige = m.group(0), [], 0            # kein Bruch zwischen zwei Vokalen (Dimensi-onen), keine Namen, keine
        if re.search(r"[a-zäöü][A-Z]", wort) or wort in NAMEN or wort in _trennung.schutz:   # Wörter, die am Seitenende brächen
            return wort
        for p in w.positions(wort):
            if 3 <= p <= len(wort) - 3 and p - vorige >= 2 and wort[p:p + 2] != "ck" \
                    and not (wort[p - 1] in VOKALE and wort[p] in VOKALE):
                stellen.append(p)
            vorige = p
        return "".join(wort[a:b] for a, b in zip([0] + stellen, stellen + [None])).replace("", "") if not stellen else \
            "­".join(wort[a:b] for a, b in zip([0] + stellen, stellen + [None]))
    return re.sub(r"(?<![\w<>/&; ⁠])[A-Za-zÄÖÜäöüß]{7,}(?![\w; ⁠])", trenne, t)


_trennung.woerter = {}
_trennung.schutz = set()
VOKALE = "aeiouyäöüAEIOUÄÖÜ"
NAMEN = {"Keplerellipse", "Keplerbahn", "Keplerbahnen", "Pantheon", "PubChem", "Planck", "Riemann", "Bertrand", "Eisenstein",
         "Schrödinger", "Coulomb", "Milchstraße", "Smaragdina", "Tabula", "Chandrasekhar", "Eddington", "Carter", "Lyapunov",
         "Preisach", "Onsager", "Casimir", "Thomas", "Fermi", "Dedekind", "Unruh", "Bernoulli", "Rayleigh", "Legendre", "Gamow",
         "Boltzmann", "Newtonian", "Newtonsche", "Newtonschen", "Hubble"}


def _kette(m):
    """Formelkette ohne Umbruch; ist sie länger als 24 Zeichen, darf sie außerhalb von Klammern vor + und − und nach dem
    Gleichheitszeichen umbrechen, sonst sperrt die Zeile davor zu weit."""
    k = m.group(0)
    if len(re.sub(r"<[^>]+>|&[a-z]+;", "", k)) <= 24:
        return k
    aus, tiefe, i = [], 0, 0
    while i < len(k):
        z = k[i]
        tiefe += (z == "(") - (z == ")")
        if tiefe == 0 and z == "\u00a0" and i + 2 < len(k) and k[i + 2] == "\u00a0":
            if k[i + 1] in "+−±↔∝→⇒":
                aus.append(" ")                                    # Umbruch vor dem Zeichen
                i += 1
                continue
            if k[i + 1] in "=≈":
                aus.append("\u00a0" + k[i + 1] + " ")              # Umbruch nach dem Gleichheitszeichen
                i += 3
                continue
        aus.append(z)
        i += 1
    return "".join(aus)


def inline(t, sprache, trennen=True):
    """Fließtext einer Zeile als HTML: Zeichen, fett, kursiv, hoch und tief; Formeln brechen nicht an ihren Zeichen um;
    mit trennen weiche Trennstellen in den Wörtern."""
    t = html.escape(zeichen(t, sprache), quote=False)         # erst die Zeichen, dann &lt; und &gt; (mit Semikolon)
    t = re.sub(r"(?<=\s)-(?=\w)", "‑", t)                     # „-Masse“: kein Umbruch nach dem Bindestrich
    for op in OPERATOREN:
        t = t.replace(f" {op} ", f" {op} ")
    t = re.sub(r"[^ ]+", _kette, t)                                # eine lange Kette darf vor + und − umbrechen
    t = re.sub(r"(?<=\S)/(?=\S)", "/⁠", t)               # kein Umbruch nach dem Bruchstrich (a/b, km/s)
    t = t.replace("ℓ", '<span class="l">ℓ</span>')              # das ℓ steht in der Schrift zu weit ab
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*([A-Za-zÀ-ž][^*\n]*?[\w.:])\*(?=[\s.,:)]|$)", r"<i>\1</i>", t)
    t = hochtief(t)
    return _trennung(t, sprache) if trennen else t


def formel(zeilen, sprache):
    """Formelblock: jede Zeile eine Reihe, Spalten an drei und mehr Leerzeichen getrennt (eine eingerückte Zeile beginnt in der
    zweiten Spalte, ein Gleichheitszeichen am Spaltenanfang gehört zur vorigen Spalte); eine Leerzeile macht Abstand."""
    reihen, abstand = [], False
    for z in zeilen:
        if not z.strip():
            abstand = True
            continue
        zellen = ([""] if z.startswith("   ") else []) + [c.strip() for c in re.split(r"\s{3,}", z.strip())]
        for k in range(len(zellen) - 1, 0, -1):
            if zellen[k].startswith("=") and zellen[k - 1]:
                zellen[k - 1:k + 1] = [zellen[k - 1] + " " + zellen[k]]
        if len(zellen) > 1:
            zellen = [c.rstrip(",") for c in zellen]
        reihen.append((zellen, abstand))
        abstand = False
    breite = max(len(z) for z, _ in reihen)
    aus = ['<table class="f">']
    for zellen, ab in reihen:
        kl = ' class="r"' if ab else ""
        if len(zellen) == 1 and breite > 1:
            aus.append(f'<tr><td{kl} colspan="{breite}">{inline(zellen[0], sprache, False)}</td></tr>')
        else:
            aus.append("<tr>" + "".join(f"<td{kl}>{inline(c, sprache, False)}</td>" for c in zellen) + "</tr>")
    aus.append("</table>")
    return "\n".join(aus)


def tabelle(zeilen, sprache):
    """Tabelle: erste Zeile Kopf, die Trennzeile fällt weg."""
    reihen = [[c.strip() for c in z.strip().strip("|").split("|")] for z in zeilen if not re.match(r"^\|\s*-", z)]
    kopf, rumpf = reihen[0], reihen[1:]
    kl = "t k" if len(rumpf) <= 12 else "t"                     # eine kurze Tabelle bleibt auf einer Seite
    aus = [f'<table class="{kl}"><thead><tr>' + "".join(f"<th>{inline(c, sprache)}</th>" for c in kopf) + "</tr></thead><tbody>"]
    for r in rumpf:
        aus.append("<tr>" + "".join(f"<td>{inline(c, sprache)}</td>" for c in r) + "</tr>")
    aus.append("</tbody></table>")
    return "\n".join(aus)


def seiten(text, sprache, titelseite=True):
    """Markdown als HTML-Rumpf: Titelseite (oder Kopf), Teile, Abschnitte, Formeln, Tabellen, Listen, Bilder, Absätze; die
    Bildbeschreibungen (Zitate) und die Trennlinien fallen weg."""
    aus, absatz, liste, nliste, block, titel, bild = [], [], [], [], None, None, 0
    zeilen = text.split("\n")

    def schliesse(vor_formel=False):
        if absatz:
            t = " ".join(absatz)
            kl = ' class="motto"' if t.startswith("*") and t.endswith(")") else (' class="vor"' if vor_formel else "")
            aus.append(f"<p{kl}>{inline(t, sprache)}</p>")
            absatz.clear()
        if liste:
            aus.append("<ul>" + "".join(f"<li>{inline(e, sprache)}</li>" for e in liste) + "</ul>")
            liste.clear()
        if nliste:
            aus.append("<ol>" + "".join(f"<li>{inline(e, sprache)}</li>" for e in nliste) + "</ol>")
            nliste.clear()

    i = 0
    while i < len(zeilen):
        z = zeilen[i]
        if block is not None:
            if z.startswith("```"):
                aus.append(formel(block, sprache))
                block = None
            else:
                block.append(z)
            i += 1
            continue
        if z.startswith("```"):
            schliesse(vor_formel=True)                       # der Satz vor der Formel bleibt bei ihr
            block = []
        elif z.startswith("> ") or z.strip() == "---" or z.startswith(">"):
            schliesse()
        elif z.startswith("|"):
            schliesse()
            j = i
            while j < len(zeilen) and zeilen[j].startswith("|"):
                j += 1
            aus.append(tabelle(zeilen[i:j], sprache))
            i = j
            continue
        elif z.startswith("# "):
            schliesse()
            if titel is None:
                titel = z[2:].strip()
                j = i + 1
                while j < len(zeilen) and not zeilen[j].strip():
                    j += 1
                unter = zeilen[j].strip().strip("*") if j < len(zeilen) and zeilen[j].startswith("**") else ""
                unter = inline(unter, sprache).replace(": ", ":<br>").replace("-", "‑")   # Zeile und Bindestrich ohne Umbruch
                k = j + 1                                    # Autorzeile nach dem Untertitel, vor der Trennlinie
                while unter and k < len(zeilen) and not zeilen[k].strip():
                    k += 1
                autor = zeilen[k].strip() if unter and k < len(zeilen) and zeilen[k].strip() not in ("", "---") \
                    and not zeilen[k].startswith("#") else ""
                zeile = f'<p class="autor">{html.escape(autor)}</p>' if autor else ""
                if titelseite:
                    aus.append(f'<div class="titelseite"><h1>{html.escape(titel)}</h1><p>{unter}</p>{zeile}</div>')
                else:
                    aus.append(f'<h1 class="kopf">{html.escape(titel)}</h1><p class="unter">{unter}</p>{zeile}')
                i = (k + 1 if autor else j + 1) if unter else i + 1
                continue
            aus.append(f'<h1 class="teil">{inline(z[2:].strip(), sprache)}</h1>')
        elif z.startswith("## "):
            schliesse()
            kl = ' class="neu"' if re.match(r"## (Anhang A|Appendix A)", z) else ""
            aus.append(f"<h2{kl}>{inline(z[3:].strip(), sprache)}</h2>")
        elif z.startswith("!["):
            schliesse()
            m = re.match(r"!\[(.*)\]\((.*)\)", z.strip())
            bild += 1
            aus.append(f'<figure><img src="{(HIER / m.group(2)).resolve().as_uri()}"><figcaption><b>{BILD[sprache]} {bild}.</b> '
                       f"{inline(m.group(1), sprache)}</figcaption></figure>")
        elif z.startswith("- "):
            if absatz:
                schliesse()
            liste.append(z[2:].strip())
        elif re.match(r"\d+\. ", z):
            if absatz:
                schliesse()
            nliste.append(re.sub(r"^\d+\. ", "", z).strip())
        elif not z.strip():
            schliesse()
        else:
            if liste or nliste:
                schliesse()
            absatz.append(z.strip())
        i += 1
    schliesse()
    return titel, "\n".join(aus)


def chromium():
    """Das Chromium des Rechners (Playwright-Ablage oder im Pfad)."""
    wege = glob.glob(os.environ.get("PLAYWRIGHT_BROWSERS_PATH", "/opt/pw-browsers") + "/chromium-*/chrome-linux/chrome")
    for w in sorted(wege, reverse=True) + ["chromium", "chromium-browser", "google-chrome"]:
        if os.path.exists(w) or any(os.path.exists(os.path.join(p, w)) for p in os.environ.get("PATH", "").split(":")):
            return w
    raise RuntimeError("kein Chromium gefunden")


def schreibe(sprache="de", name="abhandlung"):
    """Setzt <name>_<sprache>.md und schreibt <name>_<sprache>.pdf (die Abhandlung mit Titelseite, die Beispiele mit Kopf).
    Gibt den Pfad der PDF-Datei."""
    quelle = HIER / f"{name}_{sprache}.md"
    ziel = HIER / f"{name}_{sprache}.pdf"
    titel, rumpf = seiten(quelle.read_text(encoding="utf-8"), sprache, titelseite=name == "abhandlung")
    seite = (f'<!DOCTYPE html><html lang="{sprache}"><head><meta charset="utf-8"><title>{html.escape(titel)}</title>'
             f"<style>{STIL}</style></head><body>{rumpf}</body></html>")
    _setze(seite, ziel)
    for _ in range(3):                                       # Wörter, die mit Trennstrich über den Seitenwechsel brechen, schützen
        neu = _seitenbrueche(ziel) - _trennung.schutz
        if not neu:
            break
        _trennung.schutz |= neu
        titel, rumpf = seiten(quelle.read_text(encoding="utf-8"), sprache, titelseite=name == "abhandlung")
        _setze(seite.replace(seite[seite.find("<body>"):], f"<body>{rumpf}</body></html>"), ziel)
    return ziel


def _setze(seite, ziel):
    with tempfile.TemporaryDirectory() as d:
        quelle_html = pathlib.Path(d) / "abhandlung.html"
        quelle_html.write_text(seite, encoding="utf-8")
        subprocess.run([chromium(), "--headless", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
                        f"--print-to-pdf={ziel}", quelle_html.as_uri()], check=True, capture_output=True, timeout=600)


def _seitenbrueche(ziel):
    """Wörter, deren Trennstrich am Ende einer Seite steht (Rest auf der nächsten): die ganzen Wörter."""
    import pypdfium2 as pdfium
    pdf = pdfium.PdfDocument(str(ziel))
    seiten_ = [pdf[i].get_textpage().get_text_range().replace("\r", "") for i in range(len(pdf))]
    woerter = set()
    for i in range(len(seiten_) - 1):
        zeilen = [z for z in seiten_[i].strip().split("\n") if z.strip() and not re.fullmatch(r"\d+", z.strip())]
        rest = seiten_[i + 1].strip().split("\n")
        if zeilen and rest and re.search(r"[A-Za-zÄÖÜäöüß]+[‐\u00ad]$", zeilen[-1].strip()):
            anfang = re.search(r"([A-Za-zÄÖÜäöüß]+)[‐\u00ad]$", zeilen[-1].strip()).group(1)
            ende = re.match(r"[A-Za-zÄÖÜäöüß]+", rest[0].strip())
            if ende:
                woerter.add(anfang + ende.group(0))
    return woerter


def _glatt(t, sprache):
    """Text ohne Auszeichnung, Hoch- und Tiefstellung und Leerraum, zum Vergleich von Quelle und PDF."""
    t = zeichen(t, sprache).replace("**", "").replace("^{", "").replace("_{", "").replace("}", "")
    t = t.replace("^", "").replace("_", "").replace("ﬁ", "fi").replace("ﬂ", "fl")
    t = t.translate(HOCH).translate(TIEF)
    return re.sub(r"[\s|*\-‐‑‒­￾ ⁠]", "", t)   # ohne Leerraum, Trenn- und Bindestriche


def pruefe(sprache="de", name="abhandlung"):
    """Jede Zeile der Quelle (Absätze, Listen, Tabellenzeilen, Bildunterschriften, jede Spalte jeder Formelzeile) muss im
    Text der PDF stehen, bis auf Leerraum, Auszeichnung und Hoch- und Tiefstellung. Gibt die Zeilen, die fehlen."""
    import pypdfium2 as pdfium
    pdf = pdfium.PdfDocument(str(HIER / f"{name}_{sprache}.pdf"))
    seiten = [pdf[i].get_textpage().get_text_range() for i in range(len(pdf))]
    text = _glatt("".join(re.sub(r"(\r?\n[0-9 ]+)+\s*$", "", re.sub(rf"^{i + 1}\r?\n", "", s)) + "\n"
                          for i, s in enumerate(seiten)), sprache)   # ohne Seitenzahl und lose Hochzahlen am Seitenende
    fehlt, block = [], False
    for z in (HIER / f"{name}_{sprache}.md").read_text(encoding="utf-8").split("\n"):
        if z.startswith("```"):
            block = not block
            continue
        if not z.strip() or z.startswith(">") or z.strip() == "---" or re.match(r"^\|\s*-", z):
            continue
        if z.startswith("!["):
            z = re.match(r"!\[(.*)\]\(", z).group(1)
        z = re.sub(r"^\d+\. ", "", z)
        teile = [c.strip() for c in re.split(r"\s{3,}", z.strip())] if block else [z]
        if block and len(teile) > 1:
            teile = [c.rstrip(",") for c in teile]
        for t in teile:
            g = _glatt(t.lstrip("#- "), sprache)
            hoch = list(re.finditer(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻]+", t))   # eine Hochzahl kann ans Seitenende rutschen
            ohne = (_glatt((t[:h.start()] + t[h.end():]).lstrip("#- "), sprache) for h in hoch)
            if g and g not in text and not any(o in text for o in ohne):
                fehlt.append(t.strip())
    return fehlt


if __name__ == "__main__":
    sprache = sys.argv[1] if len(sys.argv) > 1 else "de"
    name = sys.argv[2] if len(sys.argv) > 2 else "abhandlung"
    print(schreibe(sprache, name))
    fehlt = pruefe(sprache, name)
    print(f"Prüfung: {len(fehlt)} Zeilen fehlen" + ("" if not fehlt else ":\n  " + "\n  ".join(fehlt)))
