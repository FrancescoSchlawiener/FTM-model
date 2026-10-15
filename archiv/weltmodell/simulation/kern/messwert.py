# -*- coding: utf-8 -*-
"""
Messwert-Register: der einzige Ort für Literatur-, Planck-, CODATA- und Datenwerte. Kein Messwert steht als Literal
in einer Formel, Prüfzahl oder einem Gesetz; dort steht MESSWERTE["kennung"].wert.

Messwert(kennung, zeichen, wert, fehler, einheit, quelle, art, datei=None, modell_texte=())
  kennung       Name im Code (bei Kernwerten der Name in konstanten.py, z. B. "H0_KMS").
  zeichen       Symbol aus dem Register, dessen Wert gemessen ist (σ_v für Crater II usw.).
  wert          Zahl in `einheit`, bei einer Spanne (unten, oben).
  fehler        None, 1σ, oder (stat, sys) bzw. (σ_b, σ_c) wie in der Quelle; Art des Fehlers steht in quelle.
  einheit       Einheit des Werts ("m/s²", "km/s/Mpc", "1" …).
  quelle        Veröffentlichung oder Datensatz, einmal benannt.
  art           "G": aus dem Gedächtnis, nicht nachgeprüft (Literatur im Netz gesperrt) – vererbt [G] an alles, was
                davon abhängt; "daten": aus weltmodell/daten (SPARC, AME2020, PubChem), `datei` nennt Datei oder Ordner,
                deren Hash in kern/datenhashes.py steht; "literatur": im Netz gegen die genannte Quelle geprüft (Autor, Jahr,
                Tabelle bzw. URL in `quelle`), kein [G].
  modell_texte  ((Datei in modell/, Zeichenfolge), …): wo der Messwert im Modelltext steht (Rückrichtung).
"""
from .symbole import SYMBOLE

REGISTERMODUL = True   # Register, keine Formeln
ARTEN = {"G", "daten", "literatur"}


class Messwert:
    __slots__ = ("kennung", "zeichen", "wert", "fehler", "einheit", "quelle", "art", "datei", "modell_texte")

    def __init__(self, kennung, zeichen, wert, fehler, einheit, quelle, art, datei=None, modell_texte=()):
        if zeichen not in SYMBOLE:
            raise KeyError(f"{kennung}: Symbol {zeichen!r} nicht im Register")
        if art not in ARTEN:
            raise ValueError(f"{kennung}: Art {art!r}")
        if (art == "daten") != (datei is not None):
            raise ValueError(f"{kennung}: datei genau bei Art 'daten'")
        if not quelle or not einheit:
            raise ValueError(f"{kennung}: Quelle oder Einheit fehlt")
        self.kennung, self.zeichen, self.wert, self.fehler = kennung, zeichen, wert, fehler
        self.einheit, self.quelle, self.art, self.datei = einheit, quelle, art, datei
        self.modell_texte = tuple(modell_texte)

    @property
    def ist_g(self):
        return self.art == "G"

    def __repr__(self):
        return f"Messwert({self.kennung!r} = {self.wert} {self.einheit}{' [G]' if self.ist_g else ''})"


class MesswertRegister(dict):
    """Kennung → Messwert; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        m = Messwert(*args, **kw)
        if m.kennung in self:
            raise ValueError(f"Messwert doppelt: {m.kennung}")
        self[m.kennung] = m
        return m


MESSWERTE = MesswertRegister()
