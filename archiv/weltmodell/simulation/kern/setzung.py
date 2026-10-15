# -*- coding: utf-8 -*-
"""
Setzungs-Register: gesetzte Punkte der tragenden Kette (Standardwerte, Lesarten, Annahmen des Modells).

Setzung(kennung, aussage, modell_ort, wert=None, einheit=None, zeichen=None, weitere_orte=())
  aussage     die Setzung in einem Satz.
  modell_ort  (Datei in modell/, Stichwort, das dort wörtlich steht): wo die Setzung steht.
  weitere_orte  ((Datei, Zeichenfolge), …): weitere Stellen, an denen der Modelltext die Setzung benutzt (Rückrichtung).
  modell_orte gibt alle Stellen.
  wert        Zahl, wenn die Setzung ein Wert ist (dann mit einheit und zeichen aus dem Symbol-Register).
Eine Setzung ist kein Messwert und kein [G]; sie wird im Bericht als Setzung genannt.
"""
from .symbole import SYMBOLE

REGISTERMODUL = True


class Setzung:
    __slots__ = ("kennung", "aussage", "modell_ort", "wert", "einheit", "zeichen", "weitere_orte")

    def __init__(self, kennung, aussage, modell_ort, wert=None, einheit=None, zeichen=None, weitere_orte=()):
        if not aussage:
            raise ValueError(f"{kennung}: Aussage fehlt")
        if len(modell_ort) != 2 or not modell_ort[0].endswith(".md"):
            raise ValueError(f"{kennung}: modell_ort = (Datei, Stichwort)")
        if (wert is None) != (einheit is None) or (wert is None) != (zeichen is None):
            raise ValueError(f"{kennung}: wert, einheit und zeichen zusammen")
        if zeichen is not None and zeichen not in SYMBOLE:
            raise KeyError(f"{kennung}: Symbol {zeichen!r} nicht im Register")
        self.kennung, self.aussage, self.modell_ort = kennung, aussage, tuple(modell_ort)
        self.wert, self.einheit, self.zeichen = wert, einheit, zeichen
        if any(len(x) != 2 or not x[0].endswith(".md") for x in weitere_orte):
            raise ValueError(f"{kennung}: weitere_orte = ((Datei, Zeichenfolge), …)")
        self.weitere_orte = tuple(tuple(x) for x in weitere_orte)

    @property
    def modell_orte(self):
        return (self.modell_ort,) + self.weitere_orte

    def __repr__(self):
        return f"Setzung({self.kennung!r})"


class SetzungRegister(dict):
    """Kennung → Setzung; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        s = Setzung(*args, **kw)
        if s.kennung in self:
            raise ValueError(f"Setzung doppelt: {s.kennung}")
        self[s.kennung] = s
        return s


SETZUNGEN = SetzungRegister()
