# -*- coding: utf-8 -*-
"""Pfad: weltmodell/simulation ins Suchverzeichnis; alle Register einmal laden."""
import os
import sys

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SIM not in sys.path:
    sys.path.insert(0, SIM)

import laden  # noqa: E402

laden.lade_alles()


def pytest_configure(config):
    config.addinivalue_line("markers", "schwer: rechnet lange (Gegenproben, Doppelformeln, alle Gesetze); der Hook lässt "
                                       "sie aus (-m 'not schwer'), lauf.py rechnet sie in eigenen Prozessen")
