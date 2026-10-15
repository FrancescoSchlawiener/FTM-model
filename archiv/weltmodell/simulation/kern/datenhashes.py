# -*- coding: utf-8 -*-
"""
Hashes der Datendateien in weltmodell/daten (SHA-256), gespeichert in daten_hashes.json neben dieser Datei.
Ein Messwert der Art "daten" nennt eine Datei oder einen Ordner daraus. Der Wächter vergleicht die gespeicherten mit
den heutigen Hashes; ändert sich eine Datei gewollt, schreibt `python -m kern.datenhashes` die Tafel neu (mit Commit).
"""
import hashlib
import json
import os

REGISTERMODUL = True
DATEN = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "daten")
TAFEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "daten_hashes.json")


def berechne(wurzel=DATEN):
    """Relativer Pfad (mit /) → SHA-256 jeder Datei unter wurzel."""
    h = {}
    for ordner, _, dateien in os.walk(wurzel):
        for d in sorted(dateien):
            p = os.path.join(ordner, d)
            with open(p, "rb") as f:
                h[os.path.relpath(p, wurzel).replace(os.sep, "/")] = hashlib.sha256(f.read()).hexdigest()
    return dict(sorted(h.items()))


def gespeichert(tafel=TAFEL):
    with open(tafel, encoding="utf-8") as f:
        return json.load(f)


def abweichungen(ist, soll):
    """Liste der Pfade, die fehlen, neu sind oder einen anderen Hash haben."""
    return sorted([f"fehlt: {p}" for p in soll.keys() - ist.keys()] + [f"neu: {p}" for p in ist.keys() - soll.keys()]
                  + [f"geändert: {p}" for p in soll.keys() & ist.keys() if soll[p] != ist[p]])


def enthaelt(datei, soll):
    """Datei oder Ordner (Präfix) steht in der Tafel."""
    return datei in soll or any(p.startswith(datei.rstrip("/") + "/") for p in soll)


if __name__ == "__main__":
    with open(TAFEL, "w", encoding="utf-8") as f:
        json.dump(berechne(), f, indent=0, ensure_ascii=False)
        f.write("\n")
    print(f"{TAFEL}: {len(berechne())} Dateien")
