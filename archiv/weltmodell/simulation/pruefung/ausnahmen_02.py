# -*- coding: utf-8 -*-
"""Ausnahmen der Rückrichtung zu Baustein 02 (Klasse in ausnahme.py)."""
from .ausnahme import AUSNAHMEN
from .modelltext import datei

_D = datei("02")
AUSNAHMEN.neu(_D, "45", "Zahlenfaktor in der Formel √𝒞 = (2π³/45)·…, keine Größe")
