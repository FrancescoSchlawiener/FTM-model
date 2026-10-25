"""The Folded Space (Text des Autors, the_folded_space.md, deutsch the_folded_space_de.md) als PDF neben der Quelle und im
Hauptverzeichnis: Formeln mit KaTeX gesetzt (npm, abgelegt in ~/.cache/ftm_katex), gedruckt mit dem Chromium des Rechners.
Wie in allen PDFs ohne Semikolon und Gedankenstrich (pdf.zeichen: Komma, bei Zahlenbereichen „to“ oder „bis“).
Lauf: python abhandlung/folded_space.py [en|de]; schreibt abhandlung/the_folded_space.pdf (deutsch the_folded_space_de.pdf)
und dieselbe Datei als The_Folded_Space.pdf (The_Folded_Space_de.pdf) ins Hauptverzeichnis."""
import html
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile

from pdf import chromium, zeichen

HIER = pathlib.Path(__file__).resolve().parent
SPRACHE = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ("en", "de") else "en"
ENDE = "" if SPRACHE == "en" else "_de"
QUELLE = HIER / f"the_folded_space{ENDE}.md"
ZIEL = HIER / f"the_folded_space{ENDE}.pdf"
HAUPT = HIER.parent / f"The_Folded_Space{ENDE}.pdf"
KOPF = {"en": "Abstract", "de": "Zusammenfassung"}[SPRACHE]
KATEX = pathlib.Path.home() / ".cache" / "ftm_katex"
VERSION = "0.16.11"

STIL = """
@page { size: A4; margin: 26mm 25mm 28mm 25mm;
        @bottom-center { content: counter(page); font-family: "KaTeX_Main", serif; font-size: 9.5pt; color: #333; } }
html { font-family: "KaTeX_Main", "FreeSerif", "DejaVu Serif", serif; font-size: 11pt; line-height: 1.45; color: #111; }
body { margin: 0; text-align: justify; hyphens: auto; }
h1 { font-size: 19pt; font-weight: bold; text-align: center; line-height: 1.25; margin: 6mm 8mm 5mm; hyphens: none; }
.autor { text-align: center; font-size: 12.5pt; font-variant: small-caps; letter-spacing: 0.04em; margin: 0 0 10mm; }
.abstract { margin: 0 11mm 9mm; font-size: 10pt; line-height: 1.42; }
.abstract h2 { text-align: center; font-size: 10.5pt; font-variant: small-caps; letter-spacing: 0.06em; margin: 0 0 2.5mm; }
.abstract p { margin: 0; text-indent: 0; }
hr { border: 0; border-top: 0.6pt solid #888; width: 30%; margin: 0 auto 8mm; }
h2 { font-size: 12.5pt; font-weight: bold; margin: 7mm 0 2.5mm; break-after: avoid; hyphens: none; }
p { margin: 0 0 2.6mm; orphans: 3; widows: 3; }
.formel { text-align: center; margin: 3.5mm 0 4mm; }
.katex { font-size: 1.04em; }
.inline { white-space: nowrap; }
.katex-display { margin: 0; }
"""


def katex():
    """Pfad des KaTeX-Pakets; fehlt es, wird es mit npm geholt."""
    paket = KATEX / "node_modules" / "katex"
    if not (paket / "dist" / "katex.min.css").exists():
        KATEX.mkdir(parents=True, exist_ok=True)
        subprocess.run(["npm", "install", "--silent", f"katex@{VERSION}", "--prefix", str(KATEX)], check=True)
    return paket


def setze_formeln(formeln):
    """Jede Formel (TeX, als Anzeige?) mit KaTeX als HTML."""
    skript = ("const k=require(process.argv[1]);let s='';process.stdin.on('data',d=>s+=d);process.stdin.on('end',()=>{"
              "const f=JSON.parse(s);process.stdout.write(JSON.stringify(f.map(([t,a])=>k.renderToString(t,"
              "{displayMode:a,throwOnError:true,output:'html'}))));});")
    r = subprocess.run(["node", "-e", skript, str(katex())], input=json.dumps(formeln), capture_output=True, text=True,
                       check=True)
    return json.loads(r.stdout)


def absatz(text, formeln):
    """Text mit $…$-Formeln: Platzhalter für die Formeln (mit anliegender Klammer und Satzzeichen ohne Umbruch), der Rest
    geschützt."""
    aus = []
    for t in re.split(r"(\(?\$[^$]+\$[.,;:)]*)", text):
        m = re.fullmatch(r"(\(?)\$([^$]+)\$([.,;:)]*)", t)
        if m:
            formeln.append((m.group(2), False))
            aus.append(f'<span class="inline">{m.group(1)}\0{len(formeln) - 1}\0{html.escape(zeichen(m.group(3), SPRACHE))}</span>')
        else:
            aus.append(html.escape(zeichen(t, SPRACHE)))
    return "".join(aus)


def seite():
    """HTML der ganzen Arbeit aus der Quelle."""
    bloecke = [b.strip() for b in QUELLE.read_text(encoding="utf-8").split("\n\n") if b.strip()]
    formeln, aus, im_abstract = [], [], False
    for i, b in enumerate(bloecke):
        if b.startswith("# "):
            aus.append(f"<h1>{absatz(b[2:], formeln)}</h1>")
        elif i == 1:
            aus.append(f'<div class="autor">{html.escape(b)}</div>')
        elif b == f"## {KOPF}":
            aus.append(f'<div class="abstract"><h2>{KOPF}</h2>')
            im_abstract = True
        elif b.startswith("## "):
            if im_abstract:
                aus.append("</div><hr>")
                im_abstract = False
            aus.append(f"<h2>{absatz(b[3:], formeln)}</h2>")
        elif b.startswith("$$") and b.endswith("$$"):
            formeln.append((b[2:-2].strip(), True))
            aus.append(f'<div class="formel">\0{len(formeln) - 1}\0</div>')
        else:
            aus.append(f"<p>{absatz(b, formeln)}</p>")
    rumpf = "\n".join(aus)
    for n, h in enumerate(setze_formeln(formeln)):
        rumpf = rumpf.replace(f"\0{n}\0", h)
    css = (katex() / "dist" / "katex.min.css").as_uri()
    titel = "The Folded Space" if SPRACHE == "en" else "Der gefaltete Raum"
    return (f'<!DOCTYPE html><html lang="{SPRACHE}"><head><meta charset="utf-8"><title>{titel}</title>'
            f'<link rel="stylesheet" href="{css}"><style>{STIL}</style></head><body>{rumpf}</body></html>')


def schreibe():
    with tempfile.TemporaryDirectory() as d:
        quelle_html = pathlib.Path(d) / f"the_folded_space{ENDE}.html"
        quelle_html.write_text(seite(), encoding="utf-8")
        subprocess.run([chromium(), "--headless", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
                        "--allow-file-access-from-files", f"--print-to-pdf={ZIEL}", quelle_html.as_uri()],
                       check=True, capture_output=True, timeout=600)
    import shutil
    shutil.copyfile(ZIEL, HAUPT)
    return ZIEL


if __name__ == "__main__":
    os.chdir(HIER)
    print(schreibe())
