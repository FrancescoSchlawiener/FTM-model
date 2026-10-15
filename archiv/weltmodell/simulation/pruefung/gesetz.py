# -*- coding: utf-8 -*-
"""
Gesetz-Register: ob das Modell die bekannten Gesetze trägt (eine Beziehung, nicht eine einzelne Zahl).

Gesetz(kennung, klasse, aussage, pruefung, quelle, ebene, triade, haengt_an=(), satz="standard", bezug=None,
         wahl=None)
  klasse    A Grenzgesetze der Physik (Lehrbuchform gegen Faltungsform), B Galaxiengesetze der Literatur,
            C Zustandsraum-Gesetze (Simulation oder zweiter Rechenweg), Z Tafeln 1–6 des Autors, M eigene Gesetze aus
            01–09, E Messwerte mit Quelle gegen das Modell (Abstand in σ).
  aussage   das Gesetz in einer Zeile.
  pruefung  Funktion p ↦ Ergebnis; p ist der Parametersatz `satz`. Ergebnis über ident, sim, mess, mess_wert,
            schranke, wahr oder hinweis.
  quelle    Autor, Veröffentlichung oder Baustein (bei Z die Tafel).
  ebene     L, I, Z, G oder Brücke ("LI" …).
  triade    Rolle in der Triade der Ebene: Hintergrund, Wirbel, Zentrum (01); auf der Ebene Z die Rollen des
            Zustandsraums tangential, transversal, parametrisch (05). Pflicht.
  haengt_an wie bei Pruefzahl (für [G]); mess() und mess_wert() tragen ihre Messwerte selbst ein.
  bezug     Kennung einer Aussage (Art "Beziehung"), wenn das Gesetz eine Aussage des Modelltexts prüft.
  wahl      None, oder ein kurzer Text: das Gesetz beruht auf einer Zuordnung, die nicht im Modelltext steht und der
            Nutzer noch bestätigen muss (Tafel 5). Das Gesetz trägt dann [W] mit diesem Text als Wurzel.
Ergebnis ∈ OK, SPANNUNG, WIDERSPRUCH, HINWEIS.

Gegenprobe-Regel: Modell- und Vergleichsseite werden als Funktionen ohne Argument übergeben (lambda: …) und getrennt
verfolgt. Ein Gesetz darf die Funktion, die es prüft, nicht als einzigen Weg nutzen: Ruft die Vergleichsseite jede
Formel der Modellseite auch selbst auf, ist es derselbe Code (regel_gegenprobe). Zahlen sind auf der Vergleichsseite
erlaubt, wenn sie von Hand gerechnet sind (Lehrbuchform); Messwerte nur über das Register.
"""
import math

import numpy as np

from kern.formeln import TRIADE
from kern.eingabe import leseprotokoll
from kern.herkunft_g import rechnung_g
from kern.herkunft_w import rechnung_w
from kern.messwert import MESSWERTE
from kern.parametersatz import PARAMETERSAETZE
from kern.symbole import EBENEN
from kern.verfolgung import verfolge
from kern.zustandsraum import ROLLEN_ZUSTANDSRAUM

from .aussage import AUSSAGEN

REGISTERMODUL = True
KLASSEN = {"A", "B", "C", "Z", "M", "E"}
ERGEBNISSE = ("OK", "SPANNUNG", "WIDERSPRUCH", "HINWEIS")
ROLLEN = TRIADE | set(ROLLEN_ZUSTANDSRAUM)


class Ergebnis:
    __slots__ = ("ergebnis", "modell", "vergleich", "mass", "messwerte", "formeln_modell", "formeln_vergleich", "art")

    def __init__(self, ergebnis, modell, vergleich, mass, messwerte=(), formeln_modell=None, formeln_vergleich=None,
                 art=""):
        if ergebnis not in ERGEBNISSE:
            raise ValueError(f"Ergebnis {ergebnis!r}")
        self.ergebnis, self.modell, self.vergleich, self.mass = ergebnis, modell, vergleich, mass
        self.messwerte, self.formeln_modell, self.formeln_vergleich = tuple(messwerte), formeln_modell, formeln_vergleich
        self.art = art


def _seite(x):
    """Wert und aufgerufene Formeln einer Seite; eine Zahl (von Hand) ruft keine Formel."""
    if callable(x):
        w, f = verfolge(x)
        return w, f
    return x, None


def _abstand(m, v):
    """Größter relativer Abstand |m/v − 1| (bei v = 0 absolut) über alle Stellen; gibt (Abstand, m, v) der Stelle."""
    m, v = np.broadcast_arrays(np.asarray(m, dtype=float), np.asarray(v, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(v != 0, np.abs(m / np.where(v != 0, v, 1.0) - 1), np.abs(m))
    i = int(np.argmax(d))
    return float(d.flat[i]), float(m.flat[i]), float(v.flat[i])


def ident(modell, vergleich, tol=1e-9):
    """Identität: relativer Abstand ≤ tol → OK, sonst WIDERSPRUCH. Bei Vergleichswert 0 absolut; bei Feldern die
    ungünstigste Stelle."""
    m, fm = _seite(modell)
    v, fv = _seite(vergleich)
    d, m, v = _abstand(m, v)
    return Ergebnis("OK" if d <= tol else "WIDERSPRUCH", m, v, f"{d:.1e}", (), fm, fv, "ident")


def sim(modell, simuliert, tol):
    """Simulation oder zweiter Rechenweg gegen die Modellform: relativer Abstand ≤ tol → OK."""
    m, fm = _seite(modell)
    s, fv = _seite(simuliert)
    d, s, m = _abstand(s, m)
    return Ergebnis("OK" if d <= tol else "WIDERSPRUCH", m, s, f"{d:.1e}", (), fm, fv, "sim")


def _urteil(n):
    return "OK" if abs(n) <= 2 else ("SPANNUNG" if abs(n) <= 3 else "WIDERSPRUCH")


def mess(modell, messwert, faktor=1.0):
    """Modellwert gegen Messwert (Kennung im Register), in dessen Einheit; faktor rechnet den Modellwert um.
    Fehler (stat, sys) werden quadratisch addiert. |n| ≤ 2 OK, ≤ 3 SPANNUNG, sonst WIDERSPRUCH."""
    m = MESSWERTE[messwert]
    f = math.hypot(*m.fehler) if isinstance(m.fehler, tuple) else m.fehler
    x, fm = _seite(modell)
    x = float(x) * faktor
    n = (x - m.wert) / f
    return Ergebnis(_urteil(n), x, m.wert, f"{n:+.1f}σ", (messwert,), fm, None, "mess")


def mess_wert(modell, wert, fehler, messwerte):
    """Modellwert gegen einen Wert mit Fehler, die aus Teilen von Messwert-Einträgen stammen (Glied eines Tupels,
    Fehler aus einem eigenen Eintrag); messwerte nennt die Kennungen (für [G])."""
    for k in messwerte:
        if k not in MESSWERTE:
            raise KeyError(f"Messwert {k!r} nicht im Register")
    x, fm = _seite(modell)
    n = (float(x) - wert) / fehler
    return Ergebnis(_urteil(n), float(x), wert, f"{n:+.1f}σ", messwerte, fm, None, "mess")


def schranke(modell, grenze, messwerte=()):
    """Ungleichung modell ≤ grenze: OK, sonst WIDERSPRUCH; Maß ist der Abstand relativ zur Grenze."""
    m, fm = _seite(modell)
    g, fv = _seite(grenze)
    m, g = float(m), float(g)
    d = (g - m) / abs(g)
    return Ergebnis("OK" if m <= g else "WIDERSPRUCH", m, g, f"Abstand {d:+.1e}", messwerte, fm, fv, "schranke")


def wahr(ok, mass="alle"):
    """Strukturaussage (Erhaltung, Symmetrie, Vorzeichen): ok wahr → OK, sonst WIDERSPRUCH."""
    return Ergebnis("OK" if ok else "WIDERSPRUCH", 1.0 if ok else 0.0, 1.0, mass, (), None, None, "wahr")


def hinweis(modell, vergleich, text, messwerte=()):
    m, fm = _seite(modell)
    v, fv = _seite(vergleich)
    return Ergebnis("HINWEIS", float(m), float(v), text, messwerte, fm, fv, "hinweis")


class Gesetz:
    __slots__ = ("kennung", "klasse", "aussage", "pruefung", "quelle", "ebene", "triade", "haengt_an", "satz", "bezug", "wahl")

    def __init__(self, kennung, klasse, aussage, pruefung, quelle, ebene, triade, haengt_an=(), satz="standard",
                 bezug=None, wahl=None):
        if klasse not in KLASSEN:
            raise ValueError(f"{kennung}: Klasse {klasse!r}")
        if not kennung.startswith(klasse + "_"):
            raise ValueError(f"{kennung}: Kennung beginnt nicht mit {klasse}_")
        if not ebene or not set(ebene) <= EBENEN:
            raise ValueError(f"{kennung}: Ebene {ebene!r}")
        if triade not in ROLLEN:
            raise ValueError(f"{kennung}: Triadenrolle {triade!r}")
        if not callable(pruefung) or not aussage or not quelle:
            raise ValueError(f"{kennung}: pruefung, aussage oder quelle fehlt")
        if bezug is not None and bezug not in AUSSAGEN:
            raise KeyError(f"{kennung}: Aussage {bezug!r} nicht angelegt")
        if satz not in PARAMETERSAETZE:
            raise KeyError(f"{kennung}: Parametersatz {satz!r}")
        self.kennung, self.klasse, self.aussage, self.pruefung, self.quelle = kennung, klasse, aussage, pruefung, quelle
        self.ebene, self.triade, self.haengt_an, self.satz, self.bezug = ebene, triade, tuple(haengt_an), satz, bezug
        self.wahl = wahl

    def pruefe(self):
        """dict(ergebnis, modell, vergleich, mass, g, w, formeln, formeln_modell, formeln_vergleich, art)."""
        p = PARAMETERSAETZE[self.satz]
        p.zuruecksetzen()
        with leseprotokoll() as gelesen:
            e, gerufen = verfolge(self.pruefung, p)
        seiten = (e.formeln_modell or set()) | (e.formeln_vergleich or set())
        g = rechnung_g(gerufen | seiten, self.haengt_an + e.messwerte, p, rechnung=self.pruefung)
        w = rechnung_w(gelesen, self.haengt_an, p, self.wahl)
        return dict(ergebnis=e.ergebnis, modell=e.modell, vergleich=e.vergleich, mass=e.mass, g=sorted(g), w=sorted(w),
                    formeln=sorted(gerufen | seiten), art=e.art,
                    formeln_modell=None if e.formeln_modell is None else sorted(e.formeln_modell),
                    formeln_vergleich=None if e.formeln_vergleich is None else sorted(e.formeln_vergleich))

    def __repr__(self):
        return f"Gesetz({self.kennung!r})"


class GesetzRegister(dict):
    """Kennung → Gesetz; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        g = Gesetz(*args, **kw)
        if g.kennung in self:
            raise ValueError(f"Gesetz doppelt: {g.kennung}")
        self[g.kennung] = g
        return g


GESETZE = GesetzRegister()


def regel_gegenprobe(ergebnisse):
    """Funde der Gegenprobe-Regel aus {kennung: pruefe()-dict}: ident, sim und schranke brauchen zwei verfolgte Seiten;
    ruft die Vergleichsseite jede Formel der Modellseite selbst auf, ist es derselbe Code."""
    funde = []
    for k, r in ergebnisse.items():
        if r["art"] not in ("ident", "sim", "schranke"):
            continue
        fm, fv = r["formeln_modell"], r["formeln_vergleich"]
        if fm is None or fv is None:
            funde.append(f"{k}: Seiten nicht als Funktionen übergeben (keine Verfolgung)")
        elif fm and set(fm) <= set(fv):
            funde.append(f"{k}: Vergleichsseite ruft {sorted(fm)} selbst auf (derselbe Code)")
    return funde
