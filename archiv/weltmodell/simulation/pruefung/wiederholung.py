# -*- coding: utf-8 -*-
"""
Wiederholungen: Eine Zahl steht an einer weiteren Fundstelle im Modelltext und meint dort dieselbe Größe wie ein
Eintrag, der sie schon trägt (Prüfzahl, Messwert, Eingabe oder Setzung). Die Rückrichtung zählt die Stelle wie den
Eintrag: gedeckt, oder offen, wenn die Prüfzahl offen ist. So bleibt jede Prüfzahl einmal, und jede Fundstelle ist
benannt.

Wiederholung(datei, zahl, stichwort, von)
  datei      Datei in modell/.
  zahl       die Zahl, wie sie an der Stelle steht.
  stichwort  Wort, das im selben Satz steht (Fundstelle).
  von        Kennung des Eintrags (Prüfzahl, Messwert, Eingabe oder Setzung), den die Stelle wiederholt.
"""
REGISTERMODUL = True


class Wiederholung:
    __slots__ = ("datei", "zahl", "stichwort", "von")

    def __init__(self, datei, zahl, stichwort, von):
        if not datei.endswith(".md") or not zahl or not stichwort or not von:
            raise ValueError(f"Wiederholung {datei} {zahl}: Datei, Zahl, Stichwort oder Eintrag fehlt")
        self.datei, self.zahl, self.stichwort, self.von = datei, zahl, stichwort, von

    def eintrag(self):
        """(Art des Eintrags, Eintrag) oder (None, None), wenn die Kennung in keinem Register steht."""
        from kern.eingabe import EINGABEN
        from kern.messwert import MESSWERTE
        from kern.setzung import SETZUNGEN
        from .pruefzahl import PRUEFZAHLEN
        for art, reg in (("Prüfzahl", PRUEFZAHLEN), ("Messwert", MESSWERTE), ("Eingabe", EINGABEN),
                         ("Setzung", SETZUNGEN)):
            if self.von in reg:
                return art, reg[self.von]
        return None, None


class WiederholungListe(list):
    def neu(self, *args):
        w = Wiederholung(*args)
        if any((v.datei, v.zahl, v.stichwort) == (w.datei, w.zahl, w.stichwort) for v in self):
            raise ValueError(f"Wiederholung doppelt: {w.datei} {w.zahl} {w.stichwort}")
        self.append(w)
        return w


WIEDERHOLUNGEN = WiederholungListe()
