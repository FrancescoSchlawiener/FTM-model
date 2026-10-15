# -*- coding: utf-8 -*-
"""
Aufrufverfolgung: welche registrierten Formeln eine Rechnung wirklich aufruft (sys.setprofile, auch verschachtelt).
Grundlage für die automatische [G]-Vererbung einer Prüfzahl und für den Wächter „Gegenprobe ruft die geprüfte
Funktion nicht auf“.
"""
import sys

from .formeln import FORMELN

REGISTERMODUL = True


def formel_codes(register=None):
    """Code-Objekt → Kennung aller Formeln im Register (Standard: FORMELN)."""
    register = FORMELN if register is None else register
    return {f.funktion.__code__: k for k, f in register.items() if f.funktion is not None}


def verfolge(rechnung, *args, register=None, **kw):
    """Führt rechnung(*args, **kw) aus; gibt (Ergebnis, Menge der aufgerufenen Formel-Kennungen)."""
    codes = formel_codes(register)
    gerufen = set()
    alt = sys.getprofile()

    def profil(frame, ereignis, arg):
        if ereignis == "call":
            k = codes.get(frame.f_code)
            if k is not None:
                gerufen.add(k)

    sys.setprofile(profil)
    try:
        ergebnis = rechnung(*args, **kw)
    finally:
        sys.setprofile(alt)
    return ergebnis, gerufen
