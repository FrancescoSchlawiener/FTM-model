# -*- coding: utf-8 -*-
"""
Gegenprobe-Register: ein unabhängiger Weg zu einer Formel des Kerns (zweiter Rechenweg, Grenzfall oder Handwert).
Derselbe Code zählt nicht: Der unabhängige Weg darf die geprüfte Formel nicht aufrufen, auch nicht verschachtelt
(Wächter über die Aufrufverfolgung). Dasselbe gilt für jede Formel des Modells, die die geprüfte Formel selbst
transitiv aufruft: ruft der unabhängige Weg sie auch, ist die Probe nicht unabhängig (Feld "teilt", zählt nicht für die
Abdeckung).

Gegenprobe(kennung, prueft, art, eingaben, unabhaengig, toleranz=1e-9, grund=None)
  prueft       Kennung der geprüften Formel (Register FORMELN).
  art          "zweiter_weg", "grenzfall", "handwert" oder "umkehr". "umkehr" kennzeichnet Umkehrproben und Rundläufe
               (die Umkehrung oder Umstellung derselben Beziehung): sie bestehen wie die anderen, zählen aber nicht für
               die Abdeckung (Wächter gegenproben_abdeckung).
  eingaben     Liste von Argument-Tupeln; je Tupel wird prueft(*args) gegen unabhaengig(*args) verglichen. Gibt die
               Formel mehrere Werte (Tupel, Feld), vergleicht die Probe jeden Wert; Zahl und Form müssen gleich sein.
  unabhaengig  Funktion mit denselben Argumenten, ohne Aufruf der geprüften Formel.
  toleranz     relativer Abstand; größer als 1e-9 nur mit grund (Grenzfall, Näherung, Handwert mit Stellen).
"""
import inspect
import math

import numpy as np

from kern.formeln import FORMELN
from kern.verfolgung import verfolge

REGISTERMODUL = True
ARTEN = {"zweiter_weg", "grenzfall", "handwert", "umkehr"}


def _zaehlt_als_formel(eintrag):
    """Formel des Modells mit mindestens einem Pflichtargument; Konstanten (ohne Argument) und Hilfen (Lader, Löser) nicht."""
    if eintrag.status != "modell" or eintrag.funktion is None:
        return False
    return any(p.default is p.empty for p in inspect.signature(eintrag.funktion).parameters.values())


class Gegenprobe:
    __slots__ = ("kennung", "prueft", "art", "eingaben", "unabhaengig", "toleranz", "grund")

    def __init__(self, kennung, prueft, art, eingaben, unabhaengig, toleranz=1e-9, grund=None):
        if art not in ARTEN:
            raise ValueError(f"{kennung}: Art {art!r}")
        if not eingaben or not callable(unabhaengig):
            raise ValueError(f"{kennung}: eingaben oder unabhaengig fehlt")
        if toleranz > 1e-9 and not grund:
            raise ValueError(f"{kennung}: Toleranz über 1e-9 nur mit Grund")
        self.kennung, self.prueft, self.art = kennung, prueft, art
        self.eingaben, self.unabhaengig, self.toleranz, self.grund = [tuple(a) for a in eingaben], unabhaengig, toleranz, grund

    def pruefe(self, register=None):
        """dict(ok, groesster_abstand, ruft_selbst, teilt): ruft_selbst ist wahr, wenn der unabhängige Weg die Formel
        ruft; teilt nennt die Formeln des Modells (mit Argument), die die geprüfte Formel transitiv aufruft und der
        unabhängige Weg ebenfalls (dann ist die Probe nicht unabhängig; Konstanten, Messwerte, Hilfen sind erlaubt)."""
        register = FORMELN if register is None else register
        f = register[self.prueft].funktion
        d_max, ruft, teilt = 0.0, False, set()
        for args in self.eingaben:
            a, innen = verfolge(f, *args, register=register)
            a = np.asarray(a, dtype=float)
            b, gerufen = verfolge(self.unabhaengig, *args, register=register)
            teilt |= {k for k in (innen & gerufen) - {self.prueft} if _zaehlt_als_formel(register[k])}
            b = np.asarray(b, dtype=float)
            ruft = ruft or self.prueft in gerufen
            if a.shape != b.shape:
                d_max = math.inf
                continue
            for x, y in zip(a.ravel(), b.ravel()):
                d_max = max(d_max, abs(y / x - 1) if x != 0 else abs(y))
        return dict(ok=d_max <= self.toleranz and not ruft, groesster_abstand=d_max, ruft_selbst=ruft, teilt=sorted(teilt))

    def __repr__(self):
        return f"Gegenprobe({self.kennung!r})"


class GegenprobeRegister(dict):
    """Kennung → Gegenprobe; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        g = Gegenprobe(*args, **kw)
        if g.kennung in self:
            raise ValueError(f"Gegenprobe doppelt: {g.kennung}")
        self[g.kennung] = g
        return g


GEGENPROBEN = GegenprobeRegister()
