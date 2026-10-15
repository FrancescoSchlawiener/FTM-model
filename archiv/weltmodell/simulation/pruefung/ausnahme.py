# -*- coding: utf-8 -*-
"""
Ausnahmen der Rückrichtung: Zahlen im Modelltext, die weder Prüfzahl, Messwert, Eingabe noch Setzung sind
(Namen wie IC 2574, Nummern, Zahlenfaktoren einer Formel, Vertrauensniveaus). Jede mit Grund.

Ausnahme(datei, zahl, grund, stichwort=None)
  datei      Datei in modell/.
  zahl       die Zahl, wie sie im Modelltext steht.
  stichwort  Fundstelle: Die Ausnahme gilt nur in Sätzen, in denen das Stichwort steht. Ohne Stichwort gilt sie für
             jede Stelle der Zahl in der Datei (nur für Namen und Faktoren, die überall dasselbe sind).
Eine behauptete Zahl ohne Kernformel oder nur aus der Forschung ist nie eine Ausnahme, sondern eine offene Prüfzahl
(Pruefzahl mit offen = "fehlt_formel" bzw. "nur_forschung"); der Grund darf darum nicht so beginnen.
"""
import re

from .modelltext import norm

REGISTERMODUL = True
OFFEN_MARKE = re.compile(r"^\s*(fehlt[ _]formel|nur[ _]forschung|offen\b)", re.I)


class Ausnahme:
    __slots__ = ("datei", "zahl", "grund", "stichwort", "text")

    def __init__(self, datei, zahl, grund, stichwort=None):
        if not datei.endswith(".md") or not grund:
            raise ValueError(f"Ausnahme {datei} {zahl}: Datei oder Grund fehlt")
        if OFFEN_MARKE.match(grund):
            raise ValueError(f"Ausnahme {datei} {zahl}: offene Zahl ist eine Prüfzahl mit offen=…, keine Ausnahme")
        self.datei, self.zahl, self.grund, self.stichwort, self.text = datei, norm(zahl), grund, stichwort, zahl


class AusnahmeListe(list):
    def neu(self, *args, **kw):
        a = Ausnahme(*args, **kw)
        if any(b.datei == a.datei and b.zahl == a.zahl and b.stichwort == a.stichwort for b in self):
            raise ValueError(f"Ausnahme doppelt: {a.datei} {a.zahl} {a.stichwort or ''}")
        self.append(a)
        return a


AUSNAHMEN = AusnahmeListe()
