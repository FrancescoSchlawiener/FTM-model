"""Ablage der Hüllenrechnung in simulation/zwischen: Ergebnisse je Etappe (lesen, schreiben) und Werte je Rechenschritt
für alle Elemente (merke), sicher gegen halb geschriebene Dateien, bei neuer Fassung überschrieben."""
import fcntl
import json
import os
import numpy as np
from huelle_1_loeser import _gitter, _scf, kerne_belegung_nl


ZWISCHEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zwischen")
# Je Etappe eine Ergebnisdatei etappe_<name>.json (die Ergebnisse je Element der ersten Etappe in etappe_konfiguration.json),
# je Rechenschritt eine Datei schritt_<art>.json mit den Werten aller Elemente {Z: {Konfiguration: Wert}}, das Startfeld aller
# Elemente in schritt_startfeld.npz. Jede Datei wird über eine Zwischendatei geschrieben, nie halb.
FASSUNG_GEKOPPELT = "gekoppelt-2: LSDA, ZORA, Gitter h = 0,016 bis 60 a₀, festes K und Startfeld aus n + l, Spin-Bahn je j"


def _pfad(name):
    return os.path.join(ZWISCHEN, name)


def _sicher(pfad, text):
    with open(pfad + ".neu", "w") as f:
        f.write(text)
    os.replace(pfad + ".neu", pfad)


def lesen():
    """Alle Ergebnisse der Hüllenrechnung, aus den Dateien je Etappe."""
    alle = {}
    if not os.path.isdir(ZWISCHEN):
        return alle
    for d in sorted(os.listdir(ZWISCHEN)):
        if d.startswith("etappe_") and d.endswith(".json"):
            name = d[len("etappe_"):-len(".json")]
            inhalt = json.load(open(_pfad(d)))
            if name == "konfiguration":
                alle.update(inhalt)
            else:
                alle["_" + name] = inhalt
    return alle


def schreiben(alle):
    """Ergebnisse ablegen, je Etappe in ihre Datei; nur geänderte Dateien werden geschrieben."""
    os.makedirs(ZWISCHEN, exist_ok=True)
    gruppen = {"konfiguration": {k: v for k, v in alle.items() if not k.startswith("_")}}
    gruppen.update({k[1:]: v for k, v in alle.items() if k.startswith("_")})
    for name, inhalt in gruppen.items():
        pfad = _pfad(f"etappe_{name}.json")
        text = json.dumps(inhalt, indent=0)
        if not os.path.exists(pfad) or open(pfad).read() != text:
            _sicher(pfad, text)
    for d in os.listdir(ZWISCHEN):                            # entfernte Etappen: ihre Datei geht mit
        if d.startswith("etappe_") and d.endswith(".json") and d[len("etappe_"):-len(".json")] not in gruppen:
            os.remove(_pfad(d))


def _startfeld(Z, g=""):
    """Startfeld des Elements (Ve je Spin, festes K) aus schritt_startfeld.npz oder None; g: Kennung des Gitters."""
    p = _pfad("schritt_startfeld.npz")
    if not os.path.exists(p):
        return None
    d = np.load(p)
    if str(d["fassung"]) != FASSUNG_GEKOPPELT or f"{Z}{g}_kf" not in d:
        return None
    return [d[f"{Z}{g}_ve0"], d[f"{Z}{g}_ve1"]], d[f"{Z}{g}_kf"]


def _zwischen(Z, R, h=0.016):
    """Zwischenspeicher eines Elements: Gitter (Schritt h) und Startfeld (Ve, K aus der Füllung nach n + l, einmal gerechnet und
    in schritt_startfeld.npz gespeichert); die Werte jedes Rechenschritts liest und schreibt merke."""
    os.makedirs(ZWISCHEN, exist_ok=True)
    g = "" if h == 0.016 else f"_h{h}"
    r, h = _gitter(h)
    sf = _startfeld(Z, g)
    if sf is None:
        _, _, Ve0, Kf, _, _ = _scf(Z, kerne_belegung_nl(Z), r, h, R)
        p = _pfad("schritt_startfeld.npz")
        with open(p + ".sperre", "w") as sp:
            fcntl.flock(sp, fcntl.LOCK_EX)
            alt = dict(np.load(p)) if os.path.exists(p) else {}
            alt.update({"fassung": FASSUNG_GEKOPPELT, f"{Z}{g}_ve0": Ve0[0], f"{Z}{g}_ve1": Ve0[1], f"{Z}{g}_kf": Kf})
            np.savez(p + ".neu.npz", **alt)
            os.replace(p + ".neu.npz", p)
    else:
        Ve0, Kf = sf
    return dict(Z=Z, R=R, r=r, h=h, Ve0=Ve0, Kf=Kf, g=g)


_SCHRITTE = {}


def _schritt(art, neu_lesen=False):
    """Werte eines Rechenschritts aller Elemente {Z: {Konfiguration: Wert}} (im Prozess gehalten, bei Bedarf neu gelesen)."""
    if neu_lesen or art not in _SCHRITTE:
        p = _pfad(f"schritt_{art}.json")
        d = json.load(open(p)) if os.path.exists(p) else {}
        _SCHRITTE[art] = d["werte"] if d.get("fassung") == FASSUNG_GEKOPPELT else {}
    return _SCHRITTE[art]


def gespeichert(zw, art):
    """Alle gespeicherten Werte des Rechenschritts art für das Element: {Konfiguration: Wert}."""
    return _schritt(art).get(str(zw["Z"]), {})


def _kennung(b):
    return " ".join(f"{n}{l}:{v}" for (n, l), v in sorted(b.items()))


def merke(zw, art, b, rechnung):
    """Wert des Rechenschritts art (etwa "niveau", "exakt", "spreizung2", "hund", "jj") für die Konfiguration b des Elements;
    aus schritt_<art>.json, fehlt er, wird rechnung() einmal ausgeführt und sofort gespeichert."""
    Z, k = str(zw["Z"]), _kennung(b)
    if k in _schritt(art).get(Z, {}) or k in _schritt(art, neu_lesen=True).get(Z, {}):
        return _SCHRITTE[art][Z][k]
    wert = rechnung()
    p = _pfad(f"schritt_{art}.json")
    with open(p + ".sperre", "w") as sp:                     # mehrere Prozesse: unter Sperre neu lesen, ergänzen
        fcntl.flock(sp, fcntl.LOCK_EX)
        werte = _schritt(art, neu_lesen=True)
        werte.setdefault(Z, {})[k] = wert
        _sicher(p, json.dumps(dict(fassung=FASSUNG_GEKOPPELT, werte=werte)))
    return wert


def _energie_gekoppelt(zw, b):
    """Energie eines Zustands mit gekoppeltem Spin und Bahn, je Konfiguration die tiefere Kopplung (hund, jj); aus dem
    Zwischenspeicher, sonst gerechnet und sofort gespeichert. Gibt (E, Kopplung)."""
    return min((merke(zw, kp, b, lambda: _scf(zw["Z"], b, zw["r"], zw["h"], zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"],
                                              kopplung=kp, maxit=300)[0]), kp) for kp in ("hund", "jj"))
