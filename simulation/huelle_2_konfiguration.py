"""Etappe 2: die Konfiguration der Hülle als Ganzes, die mit der kleinsten Faltungsenergie. Gesucht wird durch Verschieben
einzelner Elektronen zwischen den äußeren Unterschalen; mit Spin-Bahn-Glied erster Ordnung und gekoppelt je j. Dazu der
zweite Weg (Gitter, Rand, Kern) und die Treffer gegen PubChem."""
import time
from huelle_1_loeser import FM, _aussen, _gitter, _scf, kerne_belegung_nl
from huelle_ablage import _energie_gekoppelt, _zwischen, gespeichert, lesen, schreiben


def huelle(Z, R, schritte=30):
    """Konfiguration kleinster Gesamtenergie ohne und mit Spin-Bahn-Glied, Suche von n + l aus."""
    r, h = _gitter()
    start = {k: v for k, v in sorted(kerne_belegung_nl(Z).items())}
    E0, Es0, Ve0, Kf, _, sch0 = _scf(Z, start, r, h, R)
    gesehen = {tuple(sorted(start.items())): (E0, Es0)}
    aus = {}
    for art in (0, 1):
        jetzt = start
        for _ in range(schritte):
            aussen = _aussen(Z)
            besser = None
            for a in aussen:
                if jetzt.get(a, 0) == 0:
                    continue
                for b in aussen:
                    if b == a or jetzt.get(b, 0) >= 2 * (2 * b[1] + 1):
                        continue
                    neu = dict(jetzt)
                    neu[a] -= 1
                    neu[b] = neu.get(b, 0) + 1
                    neu = {k: v for k, v in neu.items() if v}
                    sch = tuple(sorted(neu.items()))
                    if sch not in gesehen:
                        e, es, *_ = _scf(Z, neu, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0)
                        gesehen[sch] = (e, es)
                    if besser is None or gesehen[sch][art] < gesehen[besser][art]:
                        besser = sch
            if besser is None or gesehen[besser][art] >= gesehen[tuple(sorted(jetzt.items()))][art] - 1e-6:
                break
            jetzt = dict(besser)
        aus[art] = jetzt
    return aus[0], aus[1], len(gesehen)

def _huelle_ein(arg):
    Z, R = arg
    t = time.time()
    try:
        a, b, n = huelle(Z, R)
    except Exception as fehler:                              # vermerken, Lauf geht weiter
        return Z, {"fehler": repr(fehler)[:200], "sekunden": round(time.time() - t, 1)}
    return Z, {"ohne": [[k[0], k[1], v] for k, v in a.items()], "mit": [[k[0], k[1], v] for k, v in b.items()],
               "konfigurationen": n, "sekunden": round(time.time() - t, 1)}

def rechne_huelle(z_bereich, prozesse=4, neu=False):
    """Rechnet die Hülle für alle Z im Bereich, speichert jedes Ergebnis sofort."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    alle = lesen() if not neu else {}
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM) for Z in z_bereich
             if str(Z) not in alle or "fehler" in alle[str(Z)]]
    with Pool(prozesse) as p:
        for Z, erg in p.imap_unordered(_huelle_ein, auftr):
            alle[str(Z)] = erg
            schreiben(alle)
            print(Z, erg["sekunden"], "s", erg.get("konfigurationen", erg.get("fehler")), flush=True)

def treffer_huelle(z_bereich, art):
    """Zahl der Elemente, deren gemessene Konfiguration die Rechnung trifft (art: 'ohne' oder 'mit' Spin-Bahn)."""
    import kerne
    gem = kerne.belegung_gemessen()
    alle = lesen()
    return int(sum({(n, l): v for n, l, v in alle[str(Z)][art]} == gem[Z] for Z in z_bereich
                   if str(Z) in alle and art in alle[str(Z)]))

GITTER_Z = (22, 29, 46, 64, 79, 92, 103, 121, 122, 124, 125, 129, 133, 138, 141, 146, 150, 157, 164, 168)

def _gitter_ein(arg):
    Z, R = arg
    _gitter.__defaults__ = (0.008, 1e-6, 60.0)
    return Z, [[k[0], k[1], v] for k, v in huelle(Z, R)[1].items()]

def gitter_pruefung(prozesse=4):
    """Zweiter Weg: dieselbe Suche mit halber Gitterweite für Stichproben; gespeichert unter '_gitter'."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM) for Z in GITTER_Z]
    with Pool(prozesse) as p:
        erg = dict(p.imap_unordered(_gitter_ein, auftr))
    alle = lesen()
    alle["_gitter"] = {str(Z): erg[Z] for Z in sorted(erg)}
    schreiben(alle)

def huelle_ergebnis():
    """Treffer an den gemessenen Konfigurationen, gerechnete Konfiguration je Z und Bestätigung durch das feinere Gitter."""
    import kerne
    alle = lesen()
    gem = kerne.belegung_gemessen()
    konf = {int(Z): {(n, l): v for n, l, v in e["mit"]} for Z, e in alle.items() if not Z.startswith("_")}
    wege = {}                                                 # zweiter Weg: halbe Gitterweite, an einzelnen Z auch Rand und Kern
    for Z, e in alle.get("_gitter", {}).items():
        wege.setdefault(int(Z), []).append({(n, l): v for n, l, v in e})
    for Z, ve in alle.get("_varianten", {}).items():
        wege.setdefault(int(Z), []).extend({(n, l): v for n, l, v in e} for e in ve.values())
    nl = lambda b: sum(kerne_belegung_nl(Z) == gem[Z] for Z in b)
    skal = {int(Z): {(n, l): v for n, l, v in e["mit"]} for Z, e in alle.get("_skaliert", {}).items()}
    anders = {Z for Z in wege if any(w != konf[Z] for w in wege[Z])}
    io = alle.get("_ionisierung", {})
    ie_neg = sorted(int(Z) for Z, v in io.items() if v[0] <= 0)       # Ion tiefer als Atom: Spin-Bahn-Glied erster Ordnung trägt nicht
    traegt = (min(ie_neg) if ie_neg else 171) - 1
    fest = [Z for Z in range(119, traegt + 1) if Z in konf and skal.get(Z) == konf[Z] and Z not in anders]   # beide Fassungen, zweiter Weg
    return dict(mit=treffer_huelle(range(1, 119), "mit"), ohne=treffer_huelle(range(1, 119), "ohne"),
                nl=nl(range(1, 119)), mit_hinten=treffer_huelle(range(87, 119), "mit"), nl_hinten=nl(range(87, 119)),
                daneben=[Z for Z in range(1, 119) if konf[Z] != gem[Z]], konf=konf,
                erst_7d=[Z for Z in fest if Z < 141 and konf[Z].get((7, 2)) and not konf[Z].get((5, 4)) and not konf[Z].get((6, 3))],
                ab_8p=min(Z for Z in range(119, 169) if all(konf[y].get((8, 1), 0) >= 2 for y in range(Z, 169))),
                fest=fest, ie_negativ=ie_neg, traegt_bis=traegt,
                schluss=sorted(Z for Z in konf if all(konf[Z].get(k) == 2 * (2 * k[1] + 1) for k in ((5, 4), (6, 3), (7, 2), (8, 0), (8, 1)))),
                gitter_gleich=sorted(Z for Z in wege if all(w == konf[Z] for w in wege[Z])),
                gitter_anders=sorted(Z for Z in wege if any(w != konf[Z] for w in wege[Z])))

VARIANTEN = {"gitter_halb": (0.008, 60.0, 1.0), "rand_weit": (0.016, 80.0, 1.0), "kern_gross": (0.016, 60.0, 1.2)}

def _variante_ein(arg):
    Z, R, name = arg
    h, rmax, fr = VARIANTEN[name]
    _gitter.__defaults__ = (h, 1e-6, rmax)
    return Z, name, [[k[0], k[1], v] for k, v in huelle(Z, R * fr)[1].items()]

def varianten(z_liste, prozesse=4):
    """Dritter Weg für die Stellen, an denen das Muster hängt: Gitter, Rand und Kernradius geändert; unter '_varianten'."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, v) for Z in z_liste for v in VARIANTEN]
    alle = lesen()
    erg = alle.get("_varianten", {})
    with Pool(prozesse) as p:
        for Z, name, k in p.imap_unordered(_variante_ein, auftr):
            erg.setdefault(str(Z), {})[name] = k
            print(Z, name, flush=True)
    alle = lesen()
    alle["_varianten"] = erg
    schreiben(alle)

def _huelle_gekoppelt_ein(arg):
    """Konfiguration kleinster Energie mit selbstkonsistent gekoppeltem Spin und Bahn (je Konfiguration die tiefere Kopplung,
    hund oder jj); festes K und Startfeld aus der Füllung nach n + l wie in huelle(). Suche wie dort: von der Füllung nach n + l
    aus je ein Elektron zwischen den äußeren Unterschalen, solange die Energie sinkt. Keine Messung, kein früheres Ergebnis."""
    Z, R = arg
    start = kerne_belegung_nl(Z)
    t = time.time()
    zw = _zwischen(Z, R)
    energie = lambda b: _energie_gekoppelt(zw, b)
    jetzt = dict(start)
    for _ in range(30):
        besser = None
        for a in _aussen(Z):
            if not jetzt.get(a):
                continue
            for b in _aussen(Z):
                if b == a or jetzt.get(b, 0) >= 2 * (2 * b[1] + 1):
                    continue
                neu = {k: v for k, v in {**jetzt, a: jetzt[a] - 1, b: jetzt.get(b, 0) + 1}.items() if v}
                if besser is None or energie(neu)[0] < energie(besser)[0]:
                    besser = neu
        if besser is None or energie(besser)[0] >= energie(jetzt)[0] - 1e-6:
            break
        jetzt = besser
    E, kp = energie(jetzt)
    return Z, dict(konf=[[n, l, v] for (n, l), v in sorted(jetzt.items())], kopplung=kp, E=E, zahl=len(gespeichert(zw, "hund")),
                   sekunden=round(time.time() - t, 1))


def rechne_gekoppelt(z_bereich=range(1, 171), prozesse=4):
    """Hülle als Ganzes mit gekoppeltem Spin und Bahn für alle Z; unter '_gekoppelt' (jedes Ergebnis sofort gespeichert)."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    fertig = lesen().get("_gekoppelt", {})
    auftr = sorted([(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM) for Z in z_bereich
                    if str(Z) not in fertig], key=lambda a: -a[0])
    with Pool(prozesse) as p:
        for Z, e in p.imap_unordered(_huelle_gekoppelt_ein, auftr):
            alle = lesen()
            alle.setdefault("_gekoppelt", {})[str(Z)] = e
            schreiben(alle)
            print(Z, e["kopplung"], e["sekunden"], flush=True)

def gekoppelt_ergebnis():
    """Hülle mit gekoppeltem Spin und Bahn gegen PubChem: Treffer (1 bis 118), Fehltreffer, Konfigurationen, Kopplung."""
    import kerne
    gem = kerne.belegung_gemessen()
    a = lesen().get("_gekoppelt", {})
    konf = {int(Z): {(n, l): v for n, l, v in e["konf"]} for Z, e in a.items()}
    bek = [Z for Z in range(1, 119) if Z in konf]
    return dict(n=len(bek), treffer=sum(konf[Z] == gem[Z] for Z in bek), daneben=[Z for Z in bek if konf[Z] != gem[Z]],
                konf=konf, kopplung={int(Z): e["kopplung"] for Z, e in a.items()})
