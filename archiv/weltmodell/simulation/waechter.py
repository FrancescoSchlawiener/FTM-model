# -*- coding: utf-8 -*-
"""
Wächter der Simulation. Jede Funktion prüft eine Regel des Bauplans und gibt die Liste der Funde zurück (leer = grün).
Jede nimmt ihre Eingaben als Argumente, damit der Test eine eingebaute Fehlerstelle zeigen kann (Negativprobe);
ohne Argument prüft sie den Bestand.

Streng (ein Fund stoppt lauf.py und den Commit): alle außer ZAEHLEND.
Zählend (Lücken werden im Bericht gezählt, bis Schritt 3c/3e sie schließt): ZAEHLEND.
"""
import ast
import io
import math
import os
import re
import tokenize

import laden  # noqa: F401  (Suchpfad)
from kern import konstanten as K
from kern.formeln import FORMELN, TRIADE  # noqa: F401
from kern.symbole import SYMBOLE

SIM = os.path.dirname(os.path.abspath(__file__))
ZAEHLEND = {"rueckrichtung", "symbolkollision", "themen_luecken", "gegenproben_abdeckung", "uhr_muster"}


# ------------------------------------------------------------------ Dateien
def py_dateien(wurzel=SIM):
    """Alle .py-Dateien der Simulation (ohne Caches)."""
    aus = []
    for ordner, unter, dateien in os.walk(wurzel):
        unter[:] = [u for u in unter if not u.startswith((".", "__"))]
        aus += [os.path.join(ordner, d) for d in sorted(dateien) if d.endswith(".py")]
    return aus


def _baum(pfad):
    with open(pfad, encoding="utf-8") as f:
        return ast.parse(f.read(), pfad)


def _rel(pfad):
    return os.path.relpath(pfad, SIM)


# ------------------------------------------------------------------ 1 Symbole eindeutig
def symbole_eindeutig(register=SYMBOLE):
    """Zeichen und Code-Namen eindeutig; ein Eintrag mit Grundzeichen hat die Form Grundzeichen_Index."""
    funde = []
    zeichen = [s.zeichen for s in register.values()]
    codes = [s.code_name for s in register.values()]
    funde += [f"Zeichen doppelt: {z}" for z in sorted({z for z in zeichen if zeichen.count(z) > 1})]
    funde += [f"Code-Name doppelt: {c}" for c in sorted({c for c in codes if codes.count(c) > 1})]
    funde += [f"Schlüssel {k} ≠ Zeichen {s.zeichen}" for k, s in register.items() if k != s.zeichen]
    for s in register.values():
        if s.stamm is not None and (s.stamm not in register or not s.zeichen.startswith(s.stamm + "_")):
            funde.append(f"{s.zeichen}: Grundzeichen {s.stamm} fehlt oder falsche Form")
    return funde


# ------------------------------------------------------------------ 2 Registrierung
def registrierung(module=None, register=FORMELN):
    """Jede Funktion in kern/ und forschung/ (außer Register- und Werkzeugmodulen) trägt genau ihren Formel-Eintrag,
    Status passt zum Ort, kein Funktionsname doppelt, jede Formel im Register gehört zu einer Funktion."""
    if module is None:
        module = [m for m in laden.module_von("kern") + laden.module_von("forschung")]
    funde, gefunden, namen = [], {}, []
    for mod in module:
        if getattr(mod, "REGISTERMODUL", False):
            continue
        ort = mod.__name__.split(".")[0]
        for n, o in vars(mod).items():
            if not (callable(o) and getattr(o, "__module__", None) == mod.__name__ and hasattr(o, "__code__")):
                continue
            namen.append(n)
            f = getattr(o, "formel", None)
            if f is None:
                funde.append(f"{mod.__name__}.{n}: ohne Formel-Eintrag")
                continue
            if f.kennung != n or register.get(n) is not f or f.funktion is not o:
                funde.append(f"{mod.__name__}.{n}: Eintrag gehört nicht zu dieser Funktion")
            erlaubt = {"modell", "hilfe"} if ort == "kern" else {"forschung"}
            if f.status not in erlaubt:
                funde.append(f"{mod.__name__}.{n}: Status {f.status} gehört nicht nach {ort}/")
            gefunden[n] = mod.__name__
    funde += [f"Funktion doppelt definiert: {n}" for n in sorted({n for n in namen if namen.count(n) > 1})]
    if module is not None:
        geprueft = {m.__name__.split(".")[0] for m in module}
        fehlend = [k for k, f in register.items()
                   if k not in gefunden and f.funktion is not None and f.funktion.__module__.split(".")[0] in geprueft]
        funde += [f"Formel {k} ohne Funktion im Modul" for k in fehlend]
    return funde


# ------------------------------------------------------------------ 3 Literale von Naturkonstanten
KONSTANTEN_WERTE = {"G": K.G, "c": K.c, "c in km/s": K.c / 1e3, "ħ": K.hbar, "h": K.hbar * 2 * math.pi, "k_B": K.kB,
                    "e": K.e_el, "m_p": K.m_p, "m_e": K.m_e, "m_u": K.m_u, "ε₀": K.eps0, "N_A": K.N_A, "M☉": K.M_sun, "GM☉": K.GM_sun, "R☉": K.R_sun,
                    "M⊕": K.M_earth, "AE": K.AU, "pc": K.pc, "kpc": K.kpc, "Mpc": K.Mpc, "Lj": K.ly, "Jahr": K.yr,
                    "Jahr (365 d)": 365 * 86400.0, "u in eV": 931.494e6, "k_B in eV/K": K.kB / K.e_el}
LITERAL_ERLAUBT = {os.path.join("kern", "konstanten.py"), os.path.join("kern", "messwerte.py"), "waechter.py"}


def literale(dateien=None, tol=2e-3):
    """Kein Zahlenliteral einer Naturkonstante oder Einheit außerhalb von konstanten.py und messwerte.py."""
    funde = []
    for p in dateien if dateien is not None else py_dateien():
        if _rel(p) in LITERAL_ERLAUBT:
            continue
        with open(p, encoding="utf-8") as f:
            toks = tokenize.generate_tokens(io.StringIO(f.read()).readline)
            for tok in toks:
                if tok.type != tokenize.NUMBER:
                    continue
                try:
                    x = float(tok.string.replace("_", ""))
                except ValueError:
                    continue
                for name, v in KONSTANTEN_WERTE.items():
                    if x != 0 and abs(x - v) <= tol * abs(v):
                        funde.append(f"{_rel(p)}:{tok.start[0]} {tok.string} ≈ {name}")
    return funde


# ------------------------------------------------------------------ 4 keine .txt-Ausgaben, kein Archiv; 8 tempsim
ARCHIV_MODULE = {"ftm_kern", "ftm_pruefstand", "ftm_gesetze", "ftm_gesamtlauf"}


def _archivnamen():
    from forschung.forschungsstand import archivskripte
    return set(archivskripte()) | ARCHIV_MODULE


def _importe(baum):
    for k in ast.walk(baum):
        if isinstance(k, ast.Import):
            for a in k.names:
                yield k.lineno, a.name
        elif isinstance(k, ast.ImportFrom) and k.module and k.level == 0:
            yield k.lineno, k.module


LESEFUNKTIONEN = {"open", "read_text", "read_bytes", "loadtxt", "genfromtxt", "load", "fromfile", "read_csv",
                  "read_table", "read", "readlines", "readline", "fromregex", "memmap"}


def _txt_konstante(knoten):
    """Erste Zeichenkette mit .txt im Ausdruck (jede Tiefe), sonst None."""
    for k in ast.walk(knoten):
        if isinstance(k, ast.Constant) and isinstance(k.value, str) and ".txt" in k.value:
            return k
    return None


def _txt_lesezugriffe(baum):
    """Echte Lesezugriffe auf .txt: Aufruf einer Lesefunktion (open, Path.read_text, np.loadtxt, …), deren Argumente
    oder Empfänger eine .txt-Zeichenkette enthalten oder eine Variable, der eine .txt-Zeichenkette zugewiesen wurde.
    Eine bloße Zeichenkette (Beschreibung, Katalogtext, Wörterbucheintrag) ist kein Zugriff. Gibt (Zeile, Text)."""
    verseucht = {}   # Variable -> Text der Zeichenkette
    for k in ast.walk(baum):
        if isinstance(k, ast.Assign):
            c = _txt_konstante(k.value)
            if c is not None:
                for t in k.targets:
                    for n in ast.walk(t):
                        if isinstance(n, ast.Name):
                            verseucht[n.id] = c.value
    for k in ast.walk(baum):
        if not isinstance(k, ast.Call):
            continue
        f = k.func
        name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
        if name not in LESEFUNKTIONEN:
            continue
        teile = list(k.args) + [w.value for w in k.keywords]
        if isinstance(f, ast.Attribute):
            teile.append(f.value)                      # Path("x.txt").read_text(), open("x.txt").read()
        for t in teile:
            c = _txt_konstante(t)
            if c is not None:
                yield k.lineno, c.value
                break
            v = next((n.id for n in ast.walk(t) if isinstance(n, ast.Name) and n.id in verseucht), None)
            if v is not None:
                yield k.lineno, verseucht[v]
                break


def quellen(dateien=None, archiv=None):
    """Kein Lesen von .txt-Ausgaben (nur echte Lesezugriffe zählen, nicht bloße Zeichenketten), kein Import aus dem
    Archiv (Skriptnamen, archiv/, rechnungen/, auch nicht über sys.path)."""
    archiv = _archivnamen() if archiv is None else archiv
    funde = []
    for p in dateien if dateien is not None else py_dateien():
        if os.path.basename(p) == "waechter.py" or os.sep + "tests" + os.sep in p:
            continue
        baum = _baum(p)
        for z, name in _importe(baum):
            teile = name.split(".")
            if teile[0] in archiv or "archiv" in teile or "rechnungen" in teile:
                funde.append(f"{_rel(p)}:{z} Import aus dem Archiv: {name}")
        for z, text in sorted(set(_txt_lesezugriffe(baum))):    # eine Zeile mit open(…).read() zählt einmal
            funde.append(f"{_rel(p)}:{z} Lesezugriff auf {text!r}")
        for k in ast.walk(baum):
            if isinstance(k, ast.Call) and ast.unparse(k.func) in ("sys.path.insert", "sys.path.append"):
                quelle = ast.unparse(k)
                if "rechnungen" in quelle or "archiv" in quelle:
                    funde.append(f"{_rel(p)}:{k.lineno} Suchpfad ins Archiv")
    return funde


def tempsim(dateien=None):
    """tempsim.py wird von keinem Modul importiert."""
    funde = []
    for p in dateien if dateien is not None else py_dateien():
        if os.path.basename(p) == "tempsim.py":
            continue
        for z, name in _importe(_baum(p)):
            if name.split(".")[-1] == "tempsim" or name.split(".")[0] == "tempsim":
                funde.append(f"{_rel(p)}:{z} importiert tempsim")
    return funde


# ------------------------------------------------------------------ 5 Modellsatz-Bindung
def modellsatz(pruefzahlen=None, aussagen=None):
    """Jede Prüfzahl: modell_text steht im selben Satz wie modell_satz in der Datei des Bausteins. Jede Aussage,
    Setzung, Eingabe mit Modelltext und jeder Messwert mit Modelltext: das Stichwort steht in der genannten Datei."""
    from kern.eingabe import EINGABEN
    from kern.messwert import MESSWERTE
    from kern.setzung import SETZUNGEN
    from pruefung.aussage import AUSSAGEN
    from pruefung.modelltext import datei, lies, saetze
    from pruefung.pruefzahl import PRUEFZAHLEN
    pruefzahlen = PRUEFZAHLEN if pruefzahlen is None else pruefzahlen
    aussagen = AUSSAGEN if aussagen is None else aussagen
    funde = []
    for p in pruefzahlen.values():
        s = saetze(lies(datei(p.baustein)))
        if not any(p.modell_text in x and p.modell_satz in x for x in s):
            funde.append(f"{p.kennung}: „{p.modell_text}“ nicht im selben Satz wie „{p.modell_satz}“ ({p.baustein})")
    orte = [(f"Aussage {a.kennung}", a.modell_ort) for a in aussagen.values()]
    orte += [(f"Setzung {s.kennung}", o) for s in SETZUNGEN.values() for o in s.modell_orte]
    orte += [(f"Eingabe {e.kennung}", t) for e in EINGABEN.values() for t in e.modell_texte]
    orte += [(f"Messwert {m.kennung}", t) for m in MESSWERTE.values() for t in m.modell_texte]
    for wer, (d, stichwort) in orte:
        if stichwort not in lies(d):
            funde.append(f"{wer}: „{stichwort}“ fehlt in {d}")
    from pruefung.ausnahme import AUSNAHMEN
    from pruefung.wiederholung import WIEDERHOLUNGEN
    for w in WIEDERHOLUNGEN:
        if w.eintrag()[0] is None:
            funde.append(f"Wiederholung {w.datei} {w.zahl}: Eintrag {w.von} in keinem Register")
        if not any(w.zahl in x and w.stichwort in x for x in saetze(lies(w.datei))):
            funde.append(f"Wiederholung {w.datei} {w.zahl}: nicht im selben Satz wie „{w.stichwort}“")
    for a in AUSNAHMEN:
        if a.stichwort and not any(a.text in x and a.stichwort in x for x in saetze(lies(a.datei))):
            funde.append(f"Ausnahme {a.datei} {a.text}: nicht im selben Satz wie „{a.stichwort}“")
    return funde


# ------------------------------------------------------------------ Rückrichtung (zählend)
class Deckung:
    """Ein Eintrag der Rückrichtung: deckt (art "gedeckt") oder erfasst (art "offen") die Zahlen in `text`, an jeder
    Stelle, an der `text` in `datei` steht, und nur in Sätzen mit `stichwort`. ueberall=True (Ausnahme ohne Stichwort):
    jede Stelle derselben Zahl in der Datei."""
    __slots__ = ("datei", "text", "stichwort", "art", "ueberall", "wer")

    def __init__(self, datei, text, stichwort=None, art="gedeckt", ueberall=False, wer=""):
        self.datei, self.text, self.stichwort, self.art, self.ueberall, self.wer = \
            datei, text, stichwort, art, ueberall, wer


def deckungen():
    """Alle Einträge der Rückrichtung aus den Registern: Prüfzahlen (offene als "offen") mit ihrem Satz, Messwerte,
    Eingaben und Setzungen mit ihren Orten, Ausnahmen mit ihrer Fundstelle."""
    from kern.eingabe import EINGABEN
    from kern.messwert import MESSWERTE
    from kern.setzung import SETZUNGEN
    from pruefung.ausnahme import AUSNAHMEN
    from pruefung.modelltext import datei
    from pruefung.pruefzahl import PRUEFZAHLEN
    from pruefung.wiederholung import WIEDERHOLUNGEN
    aus = [Deckung(datei(p.baustein), p.modell_text, p.modell_satz, "offen" if p.offen else "gedeckt", wer=p.kennung)
           for p in PRUEFZAHLEN.values()]
    aus += [Deckung(d, t, wer=m.kennung) for m in MESSWERTE.values() for d, t in m.modell_texte]
    aus += [Deckung(d, t, wer=e.kennung) for e in EINGABEN.values() for d, t in e.modell_texte]
    aus += [Deckung(d, t, wer=s.kennung) for s in SETZUNGEN.values() for d, t in s.modell_orte]
    aus += [Deckung(a.datei, a.text, a.stichwort, ueberall=a.stichwort is None, wer=f"Ausnahme {a.zahl}")
            for a in AUSNAHMEN]
    for w in WIEDERHOLUNGEN:
        _, e = w.eintrag()
        aus.append(Deckung(w.datei, w.zahl, w.stichwort, "offen" if getattr(e, "offen", None) else "gedeckt",
                           wer=f"Wiederholung {w.von}"))
    return aus


def fundstellen(texte=None, eintraege=None):
    """Datei → Liste (Anfang, Zahl, Art) aller Zahlen in modell/01–09; Art "gedeckt", "offen" oder None (weder
    Prüfzahl, Messwert, Eingabe, Setzung noch Ausnahme). Gedeckt geht vor offen."""
    from pruefung.modelltext import MODELL, lies, norm, satz_spannen, stellen
    if texte is None:
        texte = {d: lies(d) for d in sorted(os.listdir(MODELL)) if re.match(r"0[1-9]-.*\.md$", d)}
    eintraege = deckungen() if eintraege is None else eintraege
    aus = {}
    for d, t in texte.items():
        st = stellen(t)
        saetze_ = satz_spannen(t)
        art = [None] * len(st)

        def setze(i, a):
            if art[i] != "gedeckt":
                art[i] = a
        for e in (x for x in eintraege if x.datei == d):
            if e.ueberall:
                for i, (_, _, z) in enumerate(st):
                    if z == norm(e.text):
                        setze(i, e.art)
                continue
            pos = t.find(e.text)
            while pos >= 0:
                ende = pos + len(e.text)
                satz = next(((a, b) for a, b in saetze_ if a <= pos < b), (0, len(t)))
                if e.stichwort is None or e.stichwort in t[satz[0]:satz[1]]:
                    for i, (a, b, _) in enumerate(st):
                        if pos <= a and b <= ende:
                            setze(i, e.art)
                pos = t.find(e.text, pos + 1)
        aus[d] = [(a, z, art[i]) for i, (a, _, z) in enumerate(st)]
    return aus


def rueckrichtung(texte=None, eintraege=None):
    """Datei → Liste der Zahlen (je Fundstelle) in modell/01–09, die weder Prüfzahl (auch offene), Messwert, Eingabe,
    Setzung noch Ausnahme erfassen."""
    return {d: [z for _, z, a in f if a is None] for d, f in fundstellen(texte, eintraege).items()}


def rueckrichtung_bilanz(texte=None, eintraege=None):
    """Datei → dict(zahlen, erfasst, gedeckt, offen): erfasst = gedeckt + offen; offen zählt nie als gedeckt."""
    aus = {}
    for d, f in fundstellen(texte, eintraege).items():
        g = sum(a == "gedeckt" for _, _, a in f)
        o = sum(a == "offen" for _, _, a in f)
        aus[d] = dict(zahlen=len(f), erfasst=g + o, gedeckt=g, offen=o)
    return aus


# ------------------------------------------------------------------ Laufzeit
def laufzeit_warnung(dauer_s, budget_s, was="Lauf"):
    """Laufzeitgrenzen sind kein Abbruch (die Dauer hängt an der Last der Maschine): Text der Warnung mit Budget, sonst None."""
    return None if dauer_s <= budget_s else f"{was}: {dauer_s:.1f} s über dem Budget von {budget_s:g} s"


# ------------------------------------------------------------------ Gegenproben
def gegenproben(proben=None, register=FORMELN, ergebnisse=None):
    """Jede Gegenprobe besteht, und ihr unabhängiger Weg ruft die geprüfte Formel nicht auf (Aufrufverfolgung).
    ergebnisse: schon gerechnete pruefe()-Ergebnisse je Kennung (lauf.py rechnet sie einmal); sonst wird gerechnet."""
    from pruefung.gegenprobe import GEGENPROBEN
    proben = GEGENPROBEN if proben is None else proben
    funde = []
    for g in proben.values():
        if g.prueft not in register:
            funde.append(f"{g.kennung}: Formel {g.prueft} fehlt")
            continue
        r = g.pruefe(register) if ergebnisse is None else ergebnisse[g.kennung]
        if r["ruft_selbst"]:
            funde.append(f"{g.kennung}: unabhängiger Weg ruft {g.prueft} auf (derselbe Code)")
        elif not r["ok"]:
            funde.append(f"{g.kennung}: Abstand {r['groesster_abstand']:.1e} > {g.toleranz:g}")
    return funde


def _unabhaengig(g, ergebnisse):
    """Eine Probe zählt für die Abdeckung, wenn sie keine Umkehrprobe ist und keine Formel mit der geprüften teilt."""
    if g.art == "umkehr":
        return False
    r = g.pruefe() if ergebnisse is None else ergebnisse[g.kennung]
    return not r.get("teilt")


def gegenproben_abdeckung(proben=None, register=FORMELN, ergebnisse=None):
    """(zählend) Formeln des Modells ohne unabhängige Probe: ohne Gegenprobe, nur mit Umkehrprobe (art "umkehr") oder nur
    mit Proben, die eine Formel mit der geprüften teilen. ergebnisse: pruefe()-Ergebnisse je Kennung (lauf.py)."""
    from pruefung.gegenprobe import GEGENPROBEN
    proben = GEGENPROBEN if proben is None else proben
    haben = {g.prueft for g in proben.values() if _unabhaengig(g, ergebnisse)}
    return sorted(k for k, f in register.items() if f.status == "modell" and k not in haben)


def gegenproben_abdeckung_zahlen(proben=None, register=FORMELN, ergebnisse=None):
    """dict(formeln, ohne, nur_umkehr, umkehrproben, teilende): Modellformeln, davon ohne unabhängige Probe, davon mit
    Umkehrprobe oder mit teilender Probe (dann nicht unabhängig)."""
    from pruefung.gegenprobe import GEGENPROBEN
    proben = GEGENPROBEN if proben is None else proben
    ohne = gegenproben_abdeckung(proben, register, ergebnisse)
    umk = {g.prueft for g in proben.values() if g.art == "umkehr"}
    teil = {g.prueft for g in proben.values() if g.art != "umkehr" and not _unabhaengig(g, ergebnisse)}
    return dict(formeln=sum(f.status == "modell" for f in register.values()), ohne=len(ohne),
                nur_umkehr=len([k for k in ohne if k in umk]), umkehrproben=sum(g.art == "umkehr" for g in proben.values()),
                teilende=len([k for k in ohne if k in teil]))


# ------------------------------------------------------------------ Inhaltliche Doppelformeln
def _probe_argumente(fn, rng, n=8):
    """n Tupel log-gleichverteilter Zahlen über 10⁻¹⁵ … 10³ (alle Skalen des Kerns), je Pflichtargument eine."""
    import inspect
    pflicht = [p for p in inspect.signature(fn).parameters.values()
               if p.default is inspect.Parameter.empty and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
    return [tuple(10 ** rng.uniform(-15, 3) for _ in pflicht) for _ in range(n)]


def _werte(fn, argumente):
    import numpy as np
    aus = []
    for a in argumente:
        try:
            with np.errstate(all="ignore"):
                w = fn(*a)
            w = np.asarray(w, dtype=float).ravel()
        except Exception:
            return None
        if w.size == 0 or not np.all(np.isfinite(w)):
            return None
        aus.append(w)
    return aus


def doppelformeln(register=FORMELN, seed=20260930):
    """Zwei Formeln mit gleicher Einheit der Ausgabe und gleichen Werten an zufälligen Probe-Eingaben (gleiche Zahl
    der Pflichtargumente), ohne dass eine die andere als Alternative oder als gleiche Form auf anderer Ebene nennt."""
    import numpy as np
    rng = np.random.default_rng(seed)
    proben, funde = {}, []
    for k, f in register.items():
        fn = f.funktion
        if fn is None:
            continue
        args = _probe_argumente(fn, rng)
        proben[k] = (len(args[0]), args)
    werte = {}
    gruppen = {}
    for k, (n, _) in proben.items():
        gruppen.setdefault(n, []).append(k)
    for n, kennungen in gruppen.items():
        args = proben[kennungen[0]][1]
        for k in kennungen:
            werte[k] = _werte(register[k].funktion, args)
        for i, a in enumerate(kennungen):
            for b in kennungen[i + 1:]:
                wa, wb = werte[a], werte[b]
                fa, fb = register[a], register[b]
                if wa is None or wb is None or fa.alternative == b or fb.alternative == a:
                    continue
                if (b in fa.gleiche_form or a in fb.gleiche_form) and fa.ebene != fb.ebene:
                    continue   # dieselbe Form auf anderer Ebene
                if SYMBOLE[fa.aus[0]].einheit_si != SYMBOLE[fb.aus[0]].einheit_si:
                    continue
                verschieden = {round(float(x[0]), 12) for x in wa}
                if len(verschieden) < 3:
                    continue   # (fast) konstant an den Proben: sagt nichts
                if all(x.shape == y.shape and np.allclose(x, y, rtol=1e-9, atol=0) for x, y in zip(wa, wb)):
                    funde.append(f"{a} und {b} geben dieselben Werte ohne Alternative-Eintrag")
    return funde


# ------------------------------------------------------------------ Symbolkollision im Modelltext (zählend)
_GRIECH = "αβγδεζηθικλμνξοπρστυφχψωϖϑϕΓΔΘΛΞΠΣΦΨΩ"
_OPERATOR = set("πΣΔΠ∂")   # Rechenzeichen, keine Größen
_TOKEN = re.compile(rf"(?<![\w])([{_GRIECH}])(_\{{[^}}]*\}}|_[\w⊥☉⊕]+|[₀-₉ₐ-ₜ]+|[²³⁴])?")


def symbolkollision(texte=None, register=SYMBOLE):
    """(zählend) Griechische Zeichen im Modelltext mit Index oder Hochzahl, die nicht im Register stehen, obwohl ihr
    Grundbuchstabe dort eine andere Bedeutung hat (Kollision), und Grundbuchstaben ohne Registereintrag."""
    from pruefung.modelltext import MODELL, lies
    if texte is None:
        texte = {d: lies(d) for d in sorted(os.listdir(MODELL)) if re.match(r"0[0-9]-.*\.md$", d)}
    funde = []
    for d, t in texte.items():
        gesehen = set()
        for m in _TOKEN.finditer(t):
            grund, zusatz = m.group(1), (m.group(2) or "").rstrip("²³⁴")
            token = grund + zusatz
            if token in gesehen or token in register:
                continue
            gesehen.add(token)
            if grund in _OPERATOR and not zusatz:
                continue
            if grund in register:
                funde.append(f"{d}: {token} kollidiert mit {grund} = {register[grund].bedeutung}")
            elif not zusatz:
                funde.append(f"{d}: {grund} nicht im Register")
    return funde


# ------------------------------------------------------------------ Katalog und Themen
def katalog(katalog_=None, skripte=None):
    """Jeder Katalogeintrag ist gültig: Datei im Archiv; Ziel einer Übernahme ist eine Kennung (Formel, Prüfzahl,
    Gesetz, Gegenprobe); Ziel von "interessant" ist der Kurztitel einer Forschungszeile."""
    from forschung.forschungsstand import zeilen
    from katalog_eintrag import KATALOG, archivdateien
    from pruefung.gegenprobe import GEGENPROBEN
    from pruefung.gesetz import GESETZE
    from pruefung.pruefzahl import PRUEFZAHLEN
    katalog_ = KATALOG if katalog_ is None else katalog_
    skripte = set(archivdateien()) if skripte is None else set(skripte)
    titel = {z["kurztitel"] for z in zeilen()}
    kennungen = set(FORMELN) | set(PRUEFZAHLEN) | set(GESETZE) | set(GEGENPROBEN)
    funde = []
    for e in katalog_.values():
        if e.skript not in skripte:
            funde.append(f"Katalog {e.skript}: kein Archivskript")
        if e.klasse == "übernommen" and e.ziel not in kennungen:
            funde.append(f"Katalog {e.skript}: Ziel {e.ziel} ist keine Kennung")
        if e.klasse == "interessant" and e.ziel not in titel:
            funde.append(f"Katalog {e.skript}: Ziel „{e.ziel}“ ist keine Forschungszeile")
    return funde


def katalog_luecken(katalog_=None, skripte=None):
    """Archivdateien (Skripte und verwaiste Dateien) ohne Katalogeintrag."""
    from katalog_eintrag import KATALOG, archivdateien
    katalog_ = KATALOG if katalog_ is None else katalog_
    skripte = archivdateien() if skripte is None else skripte
    return sorted(s for s in skripte if s not in katalog_)


def themen(themen_=None, titel=None):
    """Jeder Versuch nennt eine vorhandene Forschungszeile (oder keine)."""
    from forschung.forschungsstand import zeilen
    from forschung.thema import THEMEN
    themen_ = THEMEN if themen_ is None else themen_
    titel = {z["kurztitel"] for z in zeilen()} if titel is None else titel
    return [f"Versuch {v.kennung}: Forschungszeile „{v.forschungszeile}“ fehlt"
            for t in themen_.values() for v in t.versuche if v.forschungszeile and v.forschungszeile not in titel]


def themen_luecken(themen_=None, zeilen_=None):
    """(zählend) Forschungszeilen mit „Rechnung“, die keinem Versuch zugeordnet sind."""
    from forschung.forschungsstand import zeilen
    from forschung.thema import THEMEN
    themen_ = THEMEN if themen_ is None else themen_
    zeilen_ = zeilen() if zeilen_ is None else zeilen_
    belegt = {v.forschungszeile for t in themen_.values() for v in t.versuche}
    return [f"Z. {z['nr']}: {z['kurztitel']}" for z in zeilen_ if z["kurztitel"] not in belegt]


# ------------------------------------------------------------------ Beobachterfehler: Epoche
EPOCHE_EINGANG = {"z", "H", "H₀"}


def epochen_eingabe(register=FORMELN):
    """Eine Formel, die eine Größe mit Zustand „epoche“ ausgibt, hat die Epoche als Eingang (z, H, H₀ oder eine
    Eingangsgröße, die selbst mit der Epoche atmet)."""
    funde = []
    for f in register.values():
        if not any(SYMBOLE[z].zustand == "epoche" for z in f.aus):
            continue
        if not any(z in EPOCHE_EINGANG or SYMBOLE[z].zustand == "epoche" for z in f.ein):
            funde.append(f"{f.kennung}: gibt {[z for z in f.aus if SYMBOLE[z].zustand == 'epoche']} (epoche) "
                         f"ohne Epoche im Eingang {f.ein}")
    return funde


def uhr_muster(texte=None, register=SYMBOLE):
    """(zählend) Sätze im Modelltext, die eine feste Größe der inneren Richtung oder des Ganzen (Zustand fest oder
    seit_z10) einer Uhrgröße (H, H₀) gleichsetzen, ohne „heute“ zu sagen (etwa a₀ ≈ cH₀/2π, Λ ≈ πc/H₀)."""
    from pruefung.modelltext import MODELL, lies, saetze
    if texte is None:
        texte = {d: lies(d) for d in sorted(os.listdir(MODELL)) if re.match(r"0[0-9]-.*\.md$", d)}
    fest = sorted((s.zeichen for s in register.values() if s.zustand in ("fest", "seit_z10") and set(s.ebene) & set("IG")
                   and len(s.zeichen) <= 3), key=len, reverse=True)
    muster = re.compile(r"(?<![\w_])(" + "|".join(map(re.escape, fest)) + r")\s*(?:≈|=)\s*[^;,.:]{0,30}?"
                        r"(?<![A-Zabd-zα-ω_])(H₀|H)(?![\w_(])")   # c davor erlaubt (cH₀)
    funde = []
    for d, t in texte.items():
        for s in saetze(t):
            for m in muster.finditer(s):
                if "heute" not in s.lower():
                    funde.append(f"{d}: „{m.group(0)}“ ohne „heute“")
    return funde


# ------------------------------------------------------------------ Zufall
def zufall(dateien=None):
    """Zufall nur über np.random.default_rng(seed) innerhalb einer Funktion (je Aufruf ein eigener Generator):
    kein Modul random, kein np.random.<Verteilung>, kein default_rng ohne Seed, kein Generator auf Modulebene."""
    funde = []
    for p in dateien if dateien is not None else py_dateien():
        if os.path.basename(p) == "waechter.py" or os.sep + "tests" + os.sep in p:
            continue
        baum = _baum(p)
        for z, name in _importe(baum):
            if name.split(".")[0] == "random":
                funde.append(f"{_rel(p)}:{z} Modul random")
        in_funktion = set()
        for k in ast.walk(baum):
            if isinstance(k, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                in_funktion |= {id(x) for x in ast.walk(k)}
        for k in ast.walk(baum):
            if not isinstance(k, ast.Call) or not isinstance(k.func, ast.Attribute):
                continue
            kette = ast.unparse(k.func)
            if ".random." not in "." + kette + "." and not kette.startswith("random."):
                continue
            if kette.endswith("default_rng"):
                if not k.args and not k.keywords:
                    funde.append(f"{_rel(p)}:{k.lineno} default_rng ohne Seed")
                if id(k) not in in_funktion:
                    funde.append(f"{_rel(p)}:{k.lineno} Generator auf Modulebene (geteilt)")
            elif not kette.endswith("Generator"):
                funde.append(f"{_rel(p)}:{k.lineno} {kette}: geteilter Zufall ohne eigenen Generator")
    return funde


# ------------------------------------------------------------------ Beschlüsse
def beschluesse(register=None):
    """Wiederöffnen braucht eine neue Grundlage: ein Beschluss mit ersetzt nennt einen vorhandenen früheren Beschluss
    und neue_grundlage."""
    from forschung.beschluss import BESCHLUESSE
    register = BESCHLUESSE if register is None else register
    funde = []
    for b in register.values():
        if b.ersetzt is None:
            continue
        if b.ersetzt not in register:
            funde.append(f"{b.kennung}: ersetzt unbekannten Beschluss {b.ersetzt}")
        if not b.neue_grundlage:
            funde.append(f"{b.kennung}: öffnet {b.ersetzt} wieder ohne neue Grundlage")
    return funde


# ------------------------------------------------------------------ Daten
def daten(ist=None, soll=None):
    """Hashes der Datendateien unverändert; jeder Messwert der Art "daten" nennt eine Datei oder einen Ordner der Tafel."""
    from kern import datenhashes as dh
    from kern.messwert import MESSWERTE
    ist = dh.berechne() if ist is None else ist
    soll = dh.gespeichert() if soll is None else soll
    funde = [f"Daten {x}" for x in dh.abweichungen(ist, soll)]
    funde += [f"Messwert {m.kennung}: Datei {m.datei} nicht in der Tafel" for m in MESSWERTE.values()
              if m.art == "daten" and not dh.enthaelt(m.datei, soll)]
    return funde


# ------------------------------------------------------------------ alle
STRENG = ("symbole_eindeutig", "registrierung", "literale", "quellen", "tempsim", "modellsatz", "gegenproben",
          "doppelformeln", "katalog", "katalog_luecken", "themen", "epochen_eingabe", "zufall", "beschluesse", "daten")


SCHWER = ("gegenproben", "gegenproben_abdeckung", "doppelformeln")   # rechnen viel; lauf.py gibt ihnen eigene Prozesse


def alle(ausser=()):
    """Name → Funde aller Wächter (streng und zählend) für den Bestand; lädt vorher alle Register. `ausser` nennt
    Wächter, die der Aufrufer selbst rechnet (SCHWER)."""
    laden.lade_alles()
    aus = {n: globals()[n]() for n in STRENG if n not in ausser}
    aus.update({n: globals()[n]() for n in sorted(ZAEHLEND) if n not in ausser})
    return aus
