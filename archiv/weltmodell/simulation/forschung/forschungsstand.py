# -*- coding: utf-8 -*-
"""
Werkzeug: die Forschungszeilen aus weltmodell/forschung/00-forschungsstand.md.

Forschungszeile = Aufzählungspunkt oder Tabellenzeile, in der "Rechnung" steht. Kurztitel = Text vor der ersten
Klammer (Aufzählung, ohne Fettdruck) bzw. erste Zelle (Tabelle); Skripte = in der Zeile genannte Namen von Archivskripten.
Ein Versuch im Themenregister und ein Katalogeintrag "interessant" nennen eine Zeile über ihren Kurztitel.
"""
import os
import re

REGISTERMODUL = True
DATEI = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "forschung",
                     "00-forschungsstand.md")
ARCHIV = os.path.join(os.path.dirname(DATEI), "..", "archiv", "rechnungen")
_NAME = re.compile(r"[a-z][a-z0-9_]*[a-z0-9]")


def archivskripte(ordner=ARCHIV):
    """Namen (ohne .py) aller Skripte im Archiv der Rechnungen."""
    return sorted(d[:-3] for d in os.listdir(ordner) if d.endswith(".py"))


def zeilen(pfad=DATEI):
    """Liste von dict(nr, kurztitel, skripte, text); skripte = in der Zeile genannte Archivskripte."""
    namen = set(archivskripte())
    aus = []
    with open(pfad, encoding="utf-8") as f:
        for nr, z in enumerate(f, 1):
            z = z.rstrip("\n")
            if "Rechnung" not in z or not (z.startswith("- ") or z.startswith("| ")):
                continue
            if z.startswith("| "):
                zellen = [x.strip() for x in z.split("|")[1:-1]]
                kurz = zellen[1] if zellen[0].isdigit() else zellen[0]
                if kurz in ("Symbol", "Rang", "Vorhersage", "Befund") or set(kurz) <= set("-: "):
                    continue
            else:
                kurz = z[2:].split(" (")[0].split(":")[0]
            kurz = kurz.replace("**", "").strip()
            skripte = sorted({n for n in _NAME.findall(z) if n in namen})
            aus.append(dict(nr=nr, kurztitel=kurz, skripte=skripte, text=z))
    return aus
