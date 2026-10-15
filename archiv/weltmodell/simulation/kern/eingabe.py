# -*- coding: utf-8 -*-
"""
Eingabe-Register: gewählte Rechenpunkte (Massen, Radien, Aufbauten), nie versteckt im Code.

Eingabe(kennung, zeichen, wert, einheit, grund, modell_text=None, wahl=False)
  zeichen     Symbol aus dem Register (M, R, r, a_D …).
  wert        Zahl in `einheit` (SI oder die Einheit, in der sie im Modelltext steht, z. B. "M☉").
  grund       warum dieser Punkt (Beispiel im Modelltext, Annahme eines Aufbaus …).
  modell_text (Datei in modell/, Zeichenfolge), wenn die Eingabe im Modelltext steht (Rückrichtung); steht sie an
              mehreren Stellen, ein Tupel solcher Paare. modell_texte gibt immer das Tupel der Paare.
  wahl        True: die Eingabe steht nicht im Modelltext und wurde von der Rechnung gewählt (Schwelle, Rechenpunkt,
              Aufbau). Sie trägt [W] (Wahl), das wie [G] an jede Prüfzahl und jedes Gesetz vererbt wird, das sie liest.
Eine Eingabe ist kein Messwert: Sie wird gewählt, nicht gemessen, und trägt kein [G].
Leseprotokoll: Wer eine Rechnung ausführt, kann mit `leseprotokoll()` festhalten, welche Eingaben sie über
EINGABEN[kennung] wirklich liest; daraus folgt [W] ohne Handarbeit (herkunft_w.py).
"""
import contextlib
from .symbole import SYMBOLE

REGISTERMODUL = True


class Eingabe:
    __slots__ = ("kennung", "zeichen", "wert", "einheit", "grund", "modell_texte", "wahl")

    def __init__(self, kennung, zeichen, wert, einheit, grund, modell_text=None, wahl=False):
        if zeichen not in SYMBOLE:
            raise KeyError(f"{kennung}: Symbol {zeichen!r} nicht im Register")
        if not grund or not einheit:
            raise ValueError(f"{kennung}: Grund oder Einheit fehlt")
        paare = () if modell_text is None else (
            (tuple(modell_text),) if isinstance(modell_text[0], str) else tuple(tuple(x) for x in modell_text))
        if any(len(x) != 2 or not x[0].endswith(".md") for x in paare):
            raise ValueError(f"{kennung}: modell_text = (Datei, Zeichenfolge) oder Tupel solcher Paare")
        self.kennung, self.zeichen, self.wert, self.einheit = kennung, zeichen, wert, einheit
        self.grund, self.modell_texte, self.wahl = grund, paare, bool(wahl)

    @property
    def modell_text(self):
        """Erste Stelle im Modelltext (oder None)."""
        return self.modell_texte[0] if self.modell_texte else None

    def __repr__(self):
        return f"Eingabe({self.kennung!r} = {self.wert} {self.einheit})"


_PROTOKOLL = []   # Stapel der offenen Leseprotokolle (Mengen von Kennungen)


@contextlib.contextmanager
def leseprotokoll():
    """Menge der Eingabe-Kennungen, die im Block über EINGABEN[kennung] gelesen werden."""
    gelesen = set()
    _PROTOKOLL.append(gelesen)
    try:
        yield gelesen
    finally:
        _PROTOKOLL.remove(gelesen)


class EingabeRegister(dict):
    """Kennung → Eingabe; doppelte Kennung wirft einen Fehler. Jeder Zugriff über [] landet in den offenen
    Leseprotokollen."""

    def __getitem__(self, kennung):
        for g in _PROTOKOLL:
            g.add(kennung)
        return super().__getitem__(kennung)

    def neu(self, *args, **kw):
        e = Eingabe(*args, **kw)
        if e.kennung in self:
            raise ValueError(f"Eingabe doppelt: {e.kennung}")
        self[e.kennung] = e
        return e


EINGABEN = EingabeRegister()
