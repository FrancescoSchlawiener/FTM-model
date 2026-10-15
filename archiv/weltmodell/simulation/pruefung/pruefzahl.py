# -*- coding: utf-8 -*-
"""
Prüfzahl-Register: jede Zahl des Modelltexts einmal, frisch aus dem Kern gerechnet.

Pruefzahl(kennung, baustein, aussage, modell_text, rechnung, einheit, modell_satz,
          haengt_an=(), satz="standard", toleranz=None, toleranz_grund=None, dex=False, offen=None, beschreibung=None)
  baustein     "01" … "09" (Datei modell/0X-*.md).
  aussage      Kennung einer Aussage (Art "Zahl"), vorher im selben Abschnitt angelegt.
  modell_text  die Zahl genau wie im Modelltext (deutsches Komma, ·10ⁿ, −); bei Spannen zwei Prüfzahlen.
  rechnung     Funktion p ↦ Wert in `einheit`; p ist der Parametersatz `satz`. Nur Kernformeln, Messwerte aus
               MESSWERTE, Eingaben aus EINGABEN, keine Zahlenliterale von Messwerten.
  einheit      Einheit des Werts, wie im Modelltext ("m", "AE", "kg", "1" …); bei offenen Prüfzahlen darf sie fehlen.
  modell_satz  Stichwort, das im selben Satz wie modell_text steht (Wächter Modellsatz-Bindung).
  haengt_an    Messwert-, Eingabe-Kennungen oder Symbol-Zeichen, die die Rechnung als Wert nutzt, ohne sie über eine
               Formel zu rufen (etwa M_sun als Eingang); daraus und aus den gerufenen Formeln folgt [G] automatisch.
  toleranz     None: halbe letzte gedruckte Stelle, absolut (Regel; auch bei „rund“/„etwa“ mit bloßer Mantisse: 4·10¹¹ → ±0,5·10¹¹). Sonst relativ, bei dex=True absolut in dex;
               dann ist toleranz_grund Pflicht.
  offen        None, oder "fehlt_formel" (die Kernformel fehlt) bzw. "nur_forschung" (die Zahl beruht im Modell nur auf
               einer Forschungsrechnung). Dann ist rechnung None und beschreibung nennt die benötigte Formel bzw. die
               Rechnung. Ergebnis OFFEN, weder OK noch ABWEICHUNG; die Rückrichtung zählt die Zahl als erfasst, nicht
               als gedeckt. Nichts im Prüfcode nachbauen.
Ergebnis: "OK", "ABWEICHUNG" oder "OFFEN".
"""
from kern.eingabe import leseprotokoll
from kern.herkunft_g import rechnung_g
from kern.herkunft_w import rechnung_w
from kern.parametersatz import PARAMETERSAETZE
from kern.verfolgung import verfolge

from .aussage import AUSSAGEN
from .modelltext import zahl_aus_text

REGISTERMODUL = True
BAUSTEINE = {f"0{i}" for i in range(1, 10)}
OFFEN = {"fehlt_formel": "fehlt Formel", "nur_forschung": "nur Forschung"}


class Pruefzahl:
    __slots__ = ("kennung", "baustein", "aussage", "modell_text", "rechnung", "einheit", "modell_satz", "haengt_an",
                 "satz", "toleranz", "toleranz_grund", "dex", "offen", "beschreibung")

    def __init__(self, kennung, baustein, aussage, modell_text, rechnung, einheit, modell_satz, haengt_an=(),
                 satz="standard", toleranz=None, toleranz_grund=None, dex=False, offen=None, beschreibung=None):
        if baustein not in BAUSTEINE:
            raise ValueError(f"{kennung}: Baustein {baustein!r}")
        if aussage not in AUSSAGEN:
            raise KeyError(f"{kennung}: Aussage {aussage!r} nicht angelegt")
        if offen is not None:
            if offen not in OFFEN:
                raise ValueError(f"{kennung}: offen {offen!r}")
            if rechnung is not None or not beschreibung:
                raise ValueError(f"{kennung}: offene Prüfzahl ohne Rechnung, mit Beschreibung")
        elif not callable(rechnung) or beschreibung:
            raise TypeError(f"{kennung}: rechnung muss eine Funktion p ↦ Wert sein (Beschreibung nur bei offen)")
        if not modell_satz or (not einheit and offen is None):
            raise ValueError(f"{kennung}: Einheit oder Modellsatz fehlt")
        _, halb = zahl_aus_text(modell_text)
        if (toleranz is None) == bool(toleranz_grund):
            raise ValueError(f"{kennung}: eigene Toleranz nur mit Grund, Grund nur mit eigener Toleranz")
        if toleranz is None and halb is None and offen is None:
            raise ValueError(f"{kennung}: bloßes 10ⁿ braucht eine eigene Toleranz")
        if satz not in PARAMETERSAETZE:
            raise KeyError(f"{kennung}: Parametersatz {satz!r}")
        self.kennung, self.baustein, self.aussage, self.modell_text = kennung, baustein, aussage, modell_text
        self.rechnung, self.einheit, self.modell_satz = rechnung, einheit, modell_satz
        self.haengt_an, self.satz = tuple(haengt_an), satz
        self.toleranz, self.toleranz_grund, self.dex = toleranz, toleranz_grund, dex
        self.offen, self.beschreibung = offen, beschreibung

    def toleranz_absolut(self):
        soll, halb = zahl_aus_text(self.modell_text)
        if self.toleranz is None:
            return halb
        return self.toleranz if self.dex else self.toleranz * abs(soll)

    def pruefe(self):
        """dict(ergebnis, wert, soll, abweichung, toleranz, ok, g, w, formeln); offen: ohne Rechnung."""
        soll, _ = zahl_aus_text(self.modell_text)
        if self.offen:
            return dict(ergebnis="OFFEN", offen=self.offen, beschreibung=self.beschreibung, wert=None, soll=soll,
                        abweichung=None, toleranz=None, ok=None, g=[], w=[], formeln=[])
        p = PARAMETERSAETZE[self.satz]
        p.zuruecksetzen()
        with leseprotokoll() as gelesen:
            wert, gerufen = verfolge(self.rechnung, p)
        wert = float(wert)
        tol = self.toleranz_absolut()
        g = rechnung_g(gerufen, self.haengt_an, p, rechnung=self.rechnung)
        w = rechnung_w(gelesen, self.haengt_an, p)
        ok = abs(wert - soll) <= tol * (1 + 1e-12)
        return dict(ergebnis="OK" if ok else "ABWEICHUNG", wert=wert, soll=soll, abweichung=wert - soll, toleranz=tol,
                    ok=ok, g=sorted(g), w=sorted(w), formeln=sorted(gerufen))

    def __repr__(self):
        return f"Pruefzahl({self.kennung!r})"


class PruefzahlRegister(dict):
    """Kennung → Prüfzahl; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        p = Pruefzahl(*args, **kw)
        if p.kennung in self:
            raise ValueError(f"Prüfzahl doppelt: {p.kennung}")
        self[p.kennung] = p
        return p


PRUEFZAHLEN = PruefzahlRegister()
