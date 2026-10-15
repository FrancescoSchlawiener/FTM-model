# -*- coding: utf-8 -*-
"""
Parametersatz: der benannte Zustand, in dem eine Rechnung läuft. Jede Prüfzahl, jedes Gesetz und jede Gegenprobe
bekommt ihren Satz als Argument p; der Bericht nennt ihn.

Parametersatz(kennung, a0_lesart, H0, Upsilon_d, z, zustand, beobachter, H0_quelle=None, sonnenkreis=None,
                ups_b_skaliert=False, wahl=None)
  a0_lesart   "standard" (Setzung A0), "horizont" (cH₀/2π, Übereinstimmung von heute), "literatur" (Messwert A0_LIT).
  H0          Atemrate heute in km/s/Mpc (Messwert H0_KMS oder gewählt).
  Upsilon_d   Masse-zu-Licht der Scheibe bei 3,6 μm (SPARC).
  z           Epoche als Rotverschiebung.
  zustand     "gebunden" (seit der Übergabe fest, Galaxien, Haufen) oder "mitatmend" (Hintergrund, Zelle).
  beobachter  für wen die Größe gilt (z. B. "hier, heute": eigener Horizont).
  sonnenkreis None, "punktmasse" (Milchstraße als Punktmasse, Sonnenkreis als Rechenpunkt) oder "gemessen" (Sonnenort
              aus den Messwerten v₀, R₀, Oort A und B, ρ_b); Rechnungen am Sonnenort prüfen, dass der Satz ihn nennt.
  ups_b_skaliert  False: Υ_b = 0,7 fest (Setzung SPARC_UPS_B); True: Υ_b wird mit Υ_d mitskaliert (Υ_b/Υ_d wie im
              Standardsatz), eine Lesart des Archivs, die der Modelltext nicht setzt.
  wahl        None, oder Text, warum der Satz eine Wahl ist (nicht im Modelltext gesetzt); jede Rechnung in diesem
              Satz trägt dann [W] (herkunft_w.py).
p.a0 ist die Schwelle nach der Lesart. Jeder Zugriff auf a0 oder H0 wird notiert (p.zugriffe); daraus erbt die
Rechnung [G], wenn der Wert aus einem Messwert [G] stammt.
"""

from .messwerte import MESSWERTE
from .setzungen import SETZUNGEN

REGISTERMODUL = True
LESARTEN = {"standard", "horizont", "literatur"}
ZUSTAENDE = {"gebunden", "mitatmend"}
SONNENKREISE = {None, "punktmasse", "gemessen"}


class Parametersatz:
    _FELDER = ("kennung", "a0_lesart", "H0", "Upsilon_d", "z", "zustand", "beobachter", "H0_quelle", "sonnenkreis",
               "ups_b_skaliert", "wahl")

    def __init__(self, kennung, a0_lesart, H0, Upsilon_d, z, zustand, beobachter, H0_quelle=None, sonnenkreis=None,
                 ups_b_skaliert=False, wahl=None):
        if a0_lesart not in LESARTEN:
            raise ValueError(f"{kennung}: a0_lesart {a0_lesart!r}")
        if zustand not in ZUSTAENDE:
            raise ValueError(f"{kennung}: zustand {zustand!r}")
        if sonnenkreis not in SONNENKREISE:
            raise ValueError(f"{kennung}: sonnenkreis {sonnenkreis!r}")
        if H0_quelle is not None and H0_quelle not in MESSWERTE:
            raise KeyError(f"{kennung}: H0_quelle {H0_quelle!r} kein Messwert")
        d = dict(kennung=kennung, a0_lesart=a0_lesart, H0=H0, Upsilon_d=Upsilon_d, z=z, zustand=zustand,
                 beobachter=beobachter, H0_quelle=H0_quelle, sonnenkreis=sonnenkreis, ups_b_skaliert=bool(ups_b_skaliert),
                 wahl=wahl)
        object.__setattr__(self, "_werte", d)
        object.__setattr__(self, "zugriffe", set())

    def __getattr__(self, name):
        werte = object.__getattribute__(self, "_werte")
        if name == "a0":
            object.__getattribute__(self, "zugriffe").add("a0")
            if self.a0_lesart == "standard":
                return SETZUNGEN["A0"].wert
            if self.a0_lesart == "literatur":
                return MESSWERTE["A0_LIT"].wert
            object.__getattribute__(self, "zugriffe").add("H0")
            from .kosmos import a0_modell          # die eine Formel a₀ = cH₀/2π (spät geladen: kosmos lädt Setzungen)
            return a0_modell(werte["H0"])
        if name in werte:
            if name in ("H0", "Upsilon_d", "z"):
                object.__getattribute__(self, "zugriffe").add(name)
            return werte[name]
        raise AttributeError(name)

    def __setattr__(self, name, wert):
        raise AttributeError("Parametersatz ist fest")

    def g_quellen(self):
        """Messwerte [G], aus denen die seit dem letzten Zurücksetzen gelesenen Felder stammen."""
        q = set()
        if "H0" in self.zugriffe and self.H0_quelle and MESSWERTE[self.H0_quelle].ist_g:
            q.add(self.H0_quelle)
        if "a0" in self.zugriffe and self.a0_lesart == "literatur":
            q.add("A0_LIT")
        return q

    def zuruecksetzen(self):
        self.zugriffe.clear()

    def text(self):
        w = self._werte
        return (f"{w['kennung']}: a₀ {w['a0_lesart']} ({self.a0:.3g} m/s²), H₀ = {w['H0']} km/s/Mpc, "
                f"Υ_d = {w['Upsilon_d']}, z = {w['z']}, {w['zustand']}, Beobachter {w['beobachter']}"
                + (f", Sonnenkreis {w['sonnenkreis']}" if w["sonnenkreis"] else "")
                + (", Υ_b mitskaliert" if w["ups_b_skaliert"] else "") + (" [W]" if w["wahl"] else ""))

    def __repr__(self):
        return f"Parametersatz({self._werte['kennung']!r})"


class ParametersatzRegister(dict):
    """Kennung → Parametersatz; doppelte Kennung wirft einen Fehler."""

    def neu(self, *args, **kw):
        p = Parametersatz(*args, **kw)
        if p.kennung in self:
            raise ValueError(f"Parametersatz doppelt: {p.kennung}")
        self[p.kennung] = p
        return p


PARAMETERSAETZE = ParametersatzRegister()
