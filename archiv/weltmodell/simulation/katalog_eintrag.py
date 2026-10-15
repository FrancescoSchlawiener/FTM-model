# -*- coding: utf-8 -*-
"""
Katalog der Archivdateien (weltmodell/archiv/rechnungen): jedes Skript einmal eingeordnet,
dazu jede verwaiste Datei (Ausgabe .txt/.npy ohne gleichnamiges Skript, übersetztes Modul in __pycache__ ohne Quelle).

KatalogEintrag(skript, klasse, ziel, grund)
  skript  Skript: Name ohne .py; verwaiste Datei: Name mit Endung (relativ zum Archiv).
  klasse  "übernommen": ziel ist die Kennung einer Formel, Prüfzahl, eines Gesetzes oder einer Gegenprobe;
          "interessant": ziel ist der Kurztitel der Forschungszeile (Versuch im Themenregister);
          "nutzlos": ziel ist None, grund sagt warum (ersetzt durch …, trägt nicht und abgeschlossen).
"""
import os

REGISTERMODUL = True
KLASSEN = {"übernommen", "interessant", "nutzlos"}
AUSGABE = ("txt", "npy")   # Endungen der Ausgaben; nur Namen, gelesen wird keine Datei


def verwaiste(ordner=None):
    """Verwaiste Dateien des Archivs: .txt/.npy ohne gleichnamiges Skript, .pyc in __pycache__ ohne Quelle."""
    from forschung.forschungsstand import ARCHIV
    ordner = ARCHIV if ordner is None else ordner
    namen = os.listdir(ordner)
    skripte = {d[:-3] for d in namen if d.endswith(".py")}
    aus = [d for d in namen if d.rsplit(".", 1)[-1] in AUSGABE and d.rsplit(".", 1)[0] not in skripte]
    cache = os.path.join(ordner, "__pycache__")
    if os.path.isdir(cache):
        aus += [f"__pycache__/{d}" for d in os.listdir(cache) if d.endswith(".pyc") and d.split(".")[0] not in skripte]
    return sorted(aus)


def archivdateien(ordner=None):
    """Alles, was der Katalog einordnet: Skripte (ohne .py) und verwaiste Dateien (mit Endung)."""
    from forschung.forschungsstand import ARCHIV, archivskripte
    ordner = ARCHIV if ordner is None else ordner
    return sorted(archivskripte(ordner) + verwaiste(ordner))


class KatalogEintrag:
    __slots__ = ("skript", "klasse", "ziel", "grund")

    def __init__(self, skript, klasse, ziel, grund):
        if klasse not in KLASSEN:
            raise ValueError(f"{skript}: Klasse {klasse!r}")
        if (klasse == "nutzlos") != (ziel is None) or not grund or skript.endswith(".py"):
            raise ValueError(f"{skript}: ziel nur bei übernommen/interessant, Grund immer, Name ohne .py")
        self.skript, self.klasse, self.ziel, self.grund = skript, klasse, ziel, grund


class Katalog(dict):
    """Skript → KatalogEintrag; ein Skript zweimal wirft einen Fehler."""

    def neu(self, *args):
        e = KatalogEintrag(*args)
        if e.skript in self:
            raise ValueError(f"Katalog doppelt: {e.skript}")
        self[e.skript] = e
        return e


KATALOG = Katalog()
