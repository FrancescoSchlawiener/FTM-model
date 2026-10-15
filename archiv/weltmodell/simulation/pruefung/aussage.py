# -*- coding: utf-8 -*-
"""
Aussage-Register: die Aussagen des Modelltexts, auf die sich Prüfzahlen und Gesetze beziehen (Rückrichtung vom
Modell zur Prüfung).

Aussage(kennung, modell_ort, art)
  modell_ort  (Datei in modell/, Stichwort, das dort wörtlich steht).
  art         "Zahl" (geprüft durch Prüfzahlen), "Beziehung" (durch Gesetze), "Vorhersage", "Setzung".
geprueft_durch(kennung) sammelt Prüfzahlen, Gesetze und Gegenproben, die die Aussage nennen; es wird nicht
gespeichert, damit es nur eine Quelle gibt.
"""
REGISTERMODUL = True
ARTEN = {"Zahl", "Beziehung", "Vorhersage", "Setzung"}


class Aussage:
    __slots__ = ("kennung", "modell_ort", "art")

    def __init__(self, kennung, modell_ort, art):
        if art not in ARTEN:
            raise ValueError(f"{kennung}: Art {art!r}")
        if len(modell_ort) != 2 or not modell_ort[0].endswith(".md"):
            raise ValueError(f"{kennung}: modell_ort = (Datei, Stichwort)")
        self.kennung, self.modell_ort, self.art = kennung, tuple(modell_ort), art

    def __repr__(self):
        return f"Aussage({self.kennung!r})"


class AussageRegister(dict):
    """Kennung → Aussage; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        a = Aussage(*args, **kw)
        if a.kennung in self:
            raise ValueError(f"Aussage doppelt: {a.kennung}")
        self[a.kennung] = a
        return a


AUSSAGEN = AussageRegister()


def geprueft_durch(kennung):
    from .pruefzahl import PRUEFZAHLEN
    from .gesetz import GESETZE
    return ([p.kennung for p in PRUEFZAHLEN.values() if p.aussage == kennung]
            + [g.kennung for g in GESETZE.values() if g.bezug == kennung])
