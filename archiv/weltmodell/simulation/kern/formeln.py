# -*- coding: utf-8 -*-
"""
Formel-Register: jede Funktion des Kerns trägt genau einen Eintrag Formel, gebunden über den Dekorator @formel (Modul formeln, Dekorator formel: verschiedene Namen, kein Umbinden).
Die Funktion bleibt unverändert (keine Hülle); der Eintrag hängt als Attribut .formel an ihr.

status: modell (steht in einer Datei in modell/), forschung (Rechnung, noch nicht im Modell), hilfe (Rechenhilfe,
keine Formel des Modells: Einheiten, Standardintegrale, Integratoren, Stoffwerte).
modell_ort: (Datei in modell/, Stichwort, das dort wörtlich steht); nur bei status "modell".
aus, ein: Zeichen aus dem Register; ein nennt die Größen, von denen die Ausgabe abhängt, ohne die
Naturkonstanten G, c, ħ, k_B; Steuergrößen ohne Symbol (Lage, Schrittweite) entfallen.
alternative: Kennung der Formel, die dasselbe Symbol auf anderem Weg ausgibt; sonst gibt jedes Symbol genau
eine Formel aus.
gleiche_form: Kennungen von Formeln derselben Form auf einer anderen Ebene (Selbstähnlichkeit, etwa ein Verhältnis
zweier Größen als Schleifenverstärkung am Ort und als Kac-Zeit im Zustandsraum); nur bei verschiedener Ebene gültig.
triade: Rolle in der Triade der Ebene (Hintergrund, Wirbel, Zentrum) oder None, wenn die Formel keine Rolle trägt.
Ort: status "modell" und "hilfe" in kern/, status "forschung" in forschung/ (Wächter).
"""
from .symbole import SYMBOLE, EBENEN

REGISTERMODUL = True
STATUS = {"modell", "forschung", "hilfe"}
TRIADE = {"Hintergrund", "Wirbel", "Zentrum"}
FORMELN = {}


class Formel:
    __slots__ = ("kennung", "gleichung", "aus", "ein", "ebene", "status", "herkunft", "modell_ort", "alternative",
                 "triade", "gleiche_form", "funktion")

    def __init__(self, kennung, gleichung, aus, ein, ebene, status, herkunft, modell_ort=None, alternative=None,
                 triade=None, gleiche_form=()):
        aus, ein = tuple(aus), tuple(ein)
        for z in aus + ein:
            if z not in SYMBOLE:
                raise KeyError(f"{kennung}: Symbol {z!r} nicht im Register")
        if not ebene or not set(ebene) <= EBENEN:
            raise ValueError(f"{kennung}: Ebene {ebene!r}")
        if status not in STATUS:
            raise ValueError(f"{kennung}: Status {status!r}")
        if not herkunft:
            raise ValueError(f"{kennung}: Herkunft fehlt")
        if triade is not None and triade not in TRIADE:
            raise ValueError(f"{kennung}: Triade {triade!r}")
        if (status == "modell") != (modell_ort is not None):
            raise ValueError(f"{kennung}: modell_ort genau bei Status 'modell'")
        if modell_ort is not None and (len(modell_ort) != 2 or not modell_ort[0].endswith(".md")):
            raise ValueError(f"{kennung}: modell_ort = (Datei, Stichwort)")
        self.kennung, self.gleichung, self.aus, self.ein = kennung, gleichung, aus, ein
        self.ebene, self.status, self.herkunft = ebene, status, herkunft
        self.modell_ort, self.alternative, self.triade, self.funktion = modell_ort, alternative, triade, None
        self.gleiche_form = (gleiche_form,) if isinstance(gleiche_form, str) else tuple(gleiche_form)

    def __repr__(self):
        return f"Formel({self.kennung!r}: {self.gleichung})"


def _registriere(f):
    if f.kennung in FORMELN:
        raise ValueError(f"Kennung doppelt: {f.kennung}")
    if f.alternative is None:
        for g in FORMELN.values():
            if g.alternative is None and set(g.aus) & set(f.aus):
                raise ValueError(f"{f.kennung} und {g.kennung} geben {set(g.aus) & set(f.aus)} aus, "
                                 "ohne Alternative")
    FORMELN[f.kennung] = f


def formel(kennung, gleichung, aus, ein, ebene, status, herkunft, modell_ort=None, alternative=None, triade=None,
           gleiche_form=()):
    """Dekorator: bindet einen Formel-Eintrag an die Funktion; Kennung = Funktionsname."""
    f = Formel(kennung, gleichung, aus, ein, ebene, status, herkunft, modell_ort, alternative, triade, gleiche_form)

    def binde(fn):
        if fn.__name__ != kennung:
            raise ValueError(f"Kennung {kennung!r} ≠ Funktionsname {fn.__name__!r}")
        if hasattr(fn, "formel"):
            raise ValueError(f"{kennung}: zweiter Formel-Eintrag")
        f.funktion = fn
        fn.formel = f
        _registriere(f)
        return fn
    return binde


def pruefe_alternativen(register=None):
    """Jede Alternative verweist auf eine Formel, die mindestens ein gleiches Symbol ausgibt; jede gleiche Form auf eine
    vorhandene Formel einer anderen Ebene. Gibt Fehlerliste."""
    register = FORMELN if register is None else register
    fehler = []
    for f in register.values():
        for k in f.gleiche_form:
            g = register.get(k)
            if g is None:
                fehler.append(f"{f.kennung}: gleiche Form {k} fehlt")
            elif g.ebene == f.ebene:
                fehler.append(f"{f.kennung}: gleiche Form {k} auf derselben Ebene {f.ebene} (Doppelformel)")
        if f.alternative is not None:
            g = register.get(f.alternative)
            if g is None:
                fehler.append(f"{f.kennung}: Alternative {f.alternative} fehlt")
            elif not set(g.aus) & set(f.aus):
                fehler.append(f"{f.kennung}: Alternative {f.alternative} gibt kein gleiches Symbol aus")
    return fehler
