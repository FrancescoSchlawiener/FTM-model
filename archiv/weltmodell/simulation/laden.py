# -*- coding: utf-8 -*-
"""
Lädt alle Module von kern/, forschung/ und pruefung/ sowie katalog.py, ohne feste Liste: Jede Datei dort füllt beim
Import ihr Register. tempsim.py wird nie geladen. Ergebnis: Liste der geladenen Module.
"""
import importlib
import os
import pkgutil
import sys

SIM = os.path.dirname(os.path.abspath(__file__))
if SIM not in sys.path:
    sys.path.insert(0, SIM)
PAKETE = ("kern", "forschung", "pruefung")
EINZELN = ("katalog_eintrag", "katalog")


def module_von(paket):
    p = importlib.import_module(paket)
    return [importlib.import_module(f"{paket}.{m.name}") for m in pkgutil.iter_modules(p.__path__)]


def lade_alles():
    mods = []
    for p in PAKETE:
        mods += module_von(p)
    mods += [importlib.import_module(n) for n in EINZELN]
    return mods
