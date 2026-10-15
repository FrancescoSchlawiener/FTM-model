# -*- coding: utf-8 -*-
"""
Themenregister gegen Mehrfachrechnung und blinde Agenten.

Invariante(form, ebene, triade)   VOR jeder Rechnung ausgefüllt: welche zeitfeste Form oder Größe, welche Ebene
                                  L/I/Z/G (Brücke mit mehreren Buchstaben), welche Rolle der Triade (Hintergrund,
                                  Wirbel, Zentrum; "keine" nur mit Grund in form).
Versuch(kennung, datum, weg, ergebnis, haengt_an=(), forschungszeile=None, skripte=())
  ergebnis        "im Modell", "Teil im Modell", "noch nicht im Modell", "trägt nicht", "offen".
  haengt_an       Messwert-Kennungen, Symbol-Zeichen oder "frei:<Größe>" (frei gewählt), an denen das Ergebnis hängt.
  forschungszeile Kurztitel der Zeile in 00-forschungsstand.md; skripte: Archivskripte (Name ohne .py).
  urteil          wie ergebnis, außer: "trägt nicht" auf [G] oder einer frei gewählten Größe wird
                  "offen, hängt an …" (kein Urteil auf ungeprüftem Grund).
Thema(kennung, frage, invariante, versuche, gesperrt, entscheider, stand)
  gesperrt        {Weg: Grund}: nicht erneut rechnen.
  entscheider     was die Frage entscheidet (Messung mit nötiger Genauigkeit, Herleitung …).
auftrag(thema, minuten) erzeugt den Kopf jedes Agentenauftrags; ohne Thema-Kennung kein Agent.
"""
from kern.formeln import TRIADE
from kern.herkunft_g import abhaengigkeit_g
from kern.symbole import EBENEN

REGISTERMODUL = True
ERGEBNISSE = {"im Modell", "Teil im Modell", "noch nicht im Modell", "trägt nicht", "offen"}


class Invariante:
    __slots__ = ("form", "ebene", "triade")

    def __init__(self, form, ebene, triade):
        if not form:
            raise ValueError("Invariante ohne Form")
        if not ebene or not set(ebene) <= EBENEN:
            raise ValueError(f"Ebene {ebene!r}")
        if triade not in TRIADE | {"keine"}:
            raise ValueError(f"Triade {triade!r}")
        self.form, self.ebene, self.triade = form, ebene, triade


class Versuch:
    __slots__ = ("kennung", "datum", "weg", "ergebnis", "haengt_an", "forschungszeile", "skripte")

    def __init__(self, kennung, datum, weg, ergebnis, haengt_an=(), forschungszeile=None, skripte=()):
        if ergebnis not in ERGEBNISSE:
            raise ValueError(f"{kennung}: Ergebnis {ergebnis!r}")
        if not weg or len(datum) != 10:
            raise ValueError(f"{kennung}: Weg fehlt oder Datum nicht JJJJ-MM-TT")
        self.kennung, self.datum, self.weg, self.ergebnis = kennung, datum, weg, ergebnis
        self.haengt_an, self.forschungszeile, self.skripte = tuple(haengt_an), forschungszeile, tuple(skripte)

    @property
    def urteil(self):
        if self.ergebnis != "trägt nicht":
            return self.ergebnis
        grund = [n for n in self.haengt_an if n.startswith("frei:") or abhaengigkeit_g(n)]
        return f"offen, hängt an {', '.join(grund)}" if grund else self.ergebnis


class Thema:
    __slots__ = ("kennung", "frage", "invariante", "versuche", "gesperrt", "entscheider", "stand")

    def __init__(self, kennung, frage, invariante, versuche, gesperrt, entscheider, stand):
        if not isinstance(invariante, Invariante):
            raise TypeError(f"{kennung}: erst die Invariante, dann rechnen")
        if not frage or not entscheider or not stand:
            raise ValueError(f"{kennung}: Frage, Entscheider oder Stand fehlt")
        kennungen = [v.kennung for v in versuche]
        if len(kennungen) != len(set(kennungen)):
            raise ValueError(f"{kennung}: Versuch doppelt")
        self.kennung, self.frage, self.invariante, self.versuche = kennung, frage, invariante, list(versuche)
        self.gesperrt, self.entscheider, self.stand = dict(gesperrt), entscheider, stand


class ThemenRegister(dict):
    """Kennung → Thema; doppelte Kennung oder ein Versuch in zwei Themen wirft einen Fehler."""

    def neu(self, *args, **kw):
        t = Thema(*args, **kw)
        if t.kennung in self:
            raise ValueError(f"Thema doppelt: {t.kennung}")
        alle = {v.kennung for u in self.values() for v in u.versuche}
        doppelt = alle & {v.kennung for v in t.versuche}
        if doppelt:
            raise ValueError(f"{t.kennung}: Versuch schon in einem anderen Thema: {doppelt}")
        self[t.kennung] = t
        return t


THEMEN = ThemenRegister()


def auftrag(thema, minuten=25):
    """Kopf eines Agentenauftrags aus dem Themeneintrag, mit allen früheren Versuchen und gesperrten Wegen."""
    t = THEMEN[thema] if isinstance(thema, str) else thema
    iv = t.invariante
    z = [f"Thema {t.kennung}: {t.frage}",
         f"Invariante: {iv.form} (Ebene {iv.ebene}, Triade {iv.triade})",
         f"Stand: {t.stand}", f"Entscheidet: {t.entscheider}", "Frühere Versuche:"]
    z += [f"  {v.kennung} ({v.datum}): {v.weg} → {v.urteil}" for v in t.versuche] or ["  keine"]
    z += ["Gesperrte Wege (nicht erneut rechnen):"]
    z += [f"  {w}: {g}" for w, g in t.gesperrt.items()] or ["  keine"]
    z += [f"Zeitgrenze: {minuten} Minuten. Formeln und Konstanten nur aus kern/. Ergebnis zuerst als Versuch in dieses "
          "Thema, dann Commit. Keine Zwischenstände fremder Agenten committen."]
    return "\n".join(z)
