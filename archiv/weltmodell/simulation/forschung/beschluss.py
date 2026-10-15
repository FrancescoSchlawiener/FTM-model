# -*- coding: utf-8 -*-
"""
Beschlüsse gegen Hin und Her: jede Aussage, die ins Modell kommt oder es verlässt, bekommt einen Beschluss.

Beschluss(kennung, aussage, modell_ort, grundlage, datum, richtung, ersetzt=None, neue_grundlage=None)
  modell_ort      (Datei in modell/, Stichwort); beim Herausnehmen das Stichwort vor der Änderung.
  grundlage       Rechnung, Daten oder Setzung, die den Beschluss trägt.
  richtung        "hinein" oder "heraus".
  ersetzt         Kennung eines früheren Beschlusses, der damit wieder geöffnet wird; dann ist neue_grundlage Pflicht
                  (neue Rechnung, neue Daten), sonst schlägt der Wächter an.
lauf.py prüft, dass jede Änderung an modell/ seit dem letzten Lauf einen Beschluss mit Datum ab diesem Lauf hat.
"""
REGISTERMODUL = True
RICHTUNGEN = {"hinein", "heraus"}


class Beschluss:
    __slots__ = ("kennung", "aussage", "modell_ort", "grundlage", "datum", "richtung", "ersetzt", "neue_grundlage")

    def __init__(self, kennung, aussage, modell_ort, grundlage, datum, richtung, ersetzt=None, neue_grundlage=None):
        if richtung not in RICHTUNGEN:
            raise ValueError(f"{kennung}: Richtung {richtung!r}")
        if len(modell_ort) != 2 or not modell_ort[0].endswith(".md"):
            raise ValueError(f"{kennung}: modell_ort = (Datei, Stichwort)")
        if not aussage or not grundlage or len(datum) != 10:
            raise ValueError(f"{kennung}: Aussage, Grundlage oder Datum (JJJJ-MM-TT) fehlt")
        self.kennung, self.aussage, self.modell_ort, self.grundlage = kennung, aussage, tuple(modell_ort), grundlage
        self.datum, self.richtung, self.ersetzt, self.neue_grundlage = datum, richtung, ersetzt, neue_grundlage


class BeschlussRegister(dict):
    """Kennung → Beschluss; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        b = Beschluss(*args, **kw)
        if b.kennung in self:
            raise ValueError(f"Beschluss doppelt: {b.kennung}")
        self[b.kennung] = b
        return b


BESCHLUESSE = BeschlussRegister()
