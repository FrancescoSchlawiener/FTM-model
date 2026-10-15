# -*- coding: utf-8 -*-
"""
Vererbung von [W] (Wahl), wie [G] in herkunft_g.py, aber für gewählte Größen statt gemessener.

[W] trägt, was nicht im Modelltext steht und von der Rechnung gewählt wurde:
  - eine Eingabe mit wahl=True (Schwelle, Rechenpunkt, Aufbau),
  - ein Parametersatz mit wahl (eine Lesart, die der Modelltext nicht setzt, etwa Υ_b mitskaliert),
  - ein Gesetz mit wahl (Zuordnung, die der Nutzer noch bestätigen muss, etwa Tafel 5).
Eine Rechnung (Prüfzahl, Gesetz, Gegenprobe) erbt [W] aus den Eingaben, die sie über EINGABEN[…] wirklich liest
(Leseprotokoll, kern/eingabe.py), aus den ausdrücklich genannten Abhängigkeiten (haengt_an) und aus ihrem
Parametersatz. Ergebnis ist immer die Menge der Wurzeln (Eingabe-Kennung, Satz-Kennung oder Text der Wahl), damit
der Bericht sagt, woran eine Zahl hängt. Eingaben, die im Modelltext stehen, tragen kein [W]; ebenso Setzungen
und Messwerte.
"""
from .eingabe import EINGABEN
from .parametersatz import PARAMETERSAETZE

REGISTERMODUL = True


def eingabe_w(kennung):
    """Wurzel [W] einer Eingabe: sie selbst, wenn sie eine Wahl ist."""
    return {kennung} if EINGABEN[kennung].wahl else set()


def abhaengigkeit_w(name):
    """Wurzeln [W] einer ausdrücklich genannten Abhängigkeit; nur Eingaben können [W] tragen, andere Namen
    (Messwerte, Formeln, Symbole) tragen keines."""
    return eingabe_w(name) if name in EINGABEN else set()


def satz_w(satz):
    """Wurzel [W] eines Parametersatzes: seine Kennung, wenn er eine Wahl ist."""
    return {satz.kennung} if satz is not None and satz.wahl else set()


def rechnung_w(gelesen=(), haengt_an=(), satz=None, wahl=None):
    """Wurzeln [W] einer Rechnung: gelesene Eingaben, genannte Abhängigkeiten, Wahl des Parametersatzes und die
    eigene Wahl (Text) der Rechnung."""
    w = set()
    for k in tuple(gelesen) + tuple(haengt_an):
        w |= abhaengigkeit_w(k)
    w |= satz_w(satz)
    if wahl:
        w.add(wahl)
    return w


def eingaben_w():
    """Alle Eingaben, die eine Wahl sind (für den Bericht)."""
    return sorted(k for k, e in EINGABEN.items() if e.wahl)


def sind_wahl(satzkennungen=None):
    """Kennungen der Parametersätze, die eine Wahl sind."""
    return sorted(k for k, p in PARAMETERSAETZE.items() if p.wahl)
