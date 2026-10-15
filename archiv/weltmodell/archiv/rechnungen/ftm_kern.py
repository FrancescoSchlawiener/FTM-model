# -*- coding: utf-8 -*-
"""Lader: stellt den Kern weltmodell/simulation/kern bereit (alter Inhalt: ftm_kern_alt.py)."""
import os
import sys

_SIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "simulation")
if _SIM not in sys.path:
    sys.path.insert(0, _SIM)
from kern import *  # noqa: F401,F403,E402
