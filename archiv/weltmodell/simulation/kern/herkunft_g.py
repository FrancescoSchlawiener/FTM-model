# -*- coding: utf-8 -*-
"""
Vererbung von [G] (Wert aus dem Gedächtnis, nicht nachgeprüft).

[G] kommt nur aus Messwert-Einträgen: aus einem Messwert der Art "G" oder einem Symbol, dessen Herkunft einen
solchen Messwert nennt ("Messwert … [G]"). Naturkonstanten und festgelegte astronomische Konstanten (Quelle CODATA
oder IAU: G, m_p, m_e, α, GM☉ …) sind Konstanten, nicht [G]; sonst trüge jede Zahl [G], und die Markierung sagte
nichts. Abgeleitete Konstanten (ℓ_P, m_Pl aus ħ, G, c) erben darum nichts.
Formel: [G] über `ein` (die Eingangssymbole mit [G]). Woher ein Eingang wirklich kommt, zeigt erst die Rechnung:
Rechnung (Prüfzahl, Gesetz, Gegenprobe): die aufgerufenen Formeln (Aufrufverfolgung, also auch die Kette über
Zwischengrößen wie H aus H₀), die ausdrücklich genannten Abhängigkeiten (Symbol-Zeichen, Messwert-Kennungen,
Formel-Kennungen), die gelesenen Felder des Parametersatzes und die Messwerte, die ihr Code unmittelbar liest
(K.H0_KMS, K.T_CMB …). Ergebnis ist immer die Menge der Wurzeln, damit der Bericht sagt, woran eine Zahl hängt.
"""
from .eingabe import EINGABEN
from .formeln import FORMELN
from .messwert import MESSWERTE
from .setzung import SETZUNGEN
from .symbole import SYMBOLE

REGISTERMODUL = True
KONSTANTEN_QUELLEN = ("CODATA", "IAU")   # Quellen fester Konstanten: kein [G]


def ist_konstante(quelle):
    """Quelle einer Naturkonstante oder festgelegten astronomischen Konstante (CODATA, IAU)."""
    return any(k in quelle for k in KONSTANTEN_QUELLEN)


def messwert_g(kennung):
    """Wurzeln [G] eines Messwert-Eintrags: er selbst, wenn er Art "G" hat und keine Konstante ist."""
    m = MESSWERTE[kennung]
    return {kennung} if m.ist_g and not ist_konstante(m.quelle) else set()


def messwerte_g(register=None):
    """Kennungen aller Messwerte, die [G] tragen (Zählung des Berichts): Art "G", ohne Konstanten (CODATA, IAU)."""
    register = MESSWERTE if register is None else register
    return sorted(k for k, m in register.items() if m.ist_g and not ist_konstante(m.quelle))


def symbol_wurzeln(zeichen, symbole=None):
    """Wurzeln [G] eines Symbols, dessen Herkunft einen Messwert [G] nennt, der keine Konstante ist: die Messwert-
    Einträge dieses Zeichens (oder von G·Zeichen, etwa GM⊕ für M⊕), die keine Konstanten sind; ohne Eintrag das
    Zeichen selbst. Abgeleitete Symbole ("aus X, Y, …") erben nichts: ihre Bausteine sind Konstanten."""
    symbole = SYMBOLE if symbole is None else symbole
    h = symbole[zeichen].herkunft
    if "[G]" not in h or ist_konstante(h):
        return set()
    eintraege = [m for m in MESSWERTE.values() if m.zeichen in (zeichen, "G" + zeichen)]
    if not eintraege:
        return {zeichen}
    return {m.kennung for m in eintraege if m.ist_g and not ist_konstante(m.quelle)}




def _namen(code):
    """Alle Namen, die ein Code-Objekt und seine inneren Funktionen (Lambdas) lesen."""
    n = set(code.co_names)
    for k in code.co_consts:
        if hasattr(k, "co_names"):
            n |= _namen(k)
    return n


def code_wurzeln(fn):
    """Wurzeln [G] der Messwerte, die eine Funktion im Rumpf unmittelbar liest (K.H0_KMS, K.T_CMB …: Name gleich der
    Messwert-Kennung, oder Code-Name eines Symbols mit [G])."""
    code = getattr(fn, "__code__", None)
    if code is None:
        return set()
    w = set()
    for name in _namen(code):
        if name in MESSWERTE:
            w |= messwert_g(name)
        s = SYMBOLE.nach_code.get(name)
        if s is not None:
            w |= symbol_wurzeln(s.zeichen)
    return w


def formel_g(kennung, register=None):
    """Wurzeln [G] einer Formel: die Wurzeln der Symbole in `ein` und der Messwerte, die ihr Rumpf unmittelbar liest."""
    register = FORMELN if register is None else register
    f = register[kennung]
    w = set()
    for z in f.ein:
        w |= symbol_wurzeln(z)
    return w | code_wurzeln(f.funktion)


def abhaengigkeit_g(name, register=None):
    """Wurzeln [G] einer ausdrücklich genannten Abhängigkeit: Messwert-, Eingabe-, Setzungs- oder Formel-Kennung,
    Symbol-Zeichen oder "frei:<Größe>" (frei gewählt). Eingaben, Setzungen und freie Größen tragen kein [G]."""
    register = FORMELN if register is None else register
    if name.startswith("frei:") or name in EINGABEN or name in SETZUNGEN:
        return set()
    if name in MESSWERTE:
        return messwert_g(name)
    if name in register:
        return formel_g(name, register)
    if name in SYMBOLE:
        return symbol_wurzeln(name)
    raise KeyError(f"Abhängigkeit {name!r} weder Messwert, Formel noch Symbol")


def rechnung_g(formeln_gerufen=(), haengt_an=(), satz=None, register=None, rechnung=None):
    """Wurzeln [G] einer Rechnung: gerufene Formeln, genannte Abhängigkeiten, gelesene Felder des Satzes und die
    Messwerte, die die Rechnung selbst unmittelbar liest."""
    wurzeln = code_wurzeln(rechnung) if rechnung is not None else set()
    for k in formeln_gerufen:
        wurzeln |= formel_g(k, register)
    for n in haengt_an:
        wurzeln |= abhaengigkeit_g(n, register)
    if satz is not None:
        wurzeln |= {q for q in satz.g_quellen() if messwert_g(q)}
    return wurzeln
