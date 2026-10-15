# -*- coding: utf-8 -*-
"""
Einziger Aufruf der Simulation: python lauf.py  (mit --frisch ohne den Zwischenspeicher zu lesen)

Aufgaben: 1. Wächter (pytest in tests/, dazu die schweren Wächter; ein strenger Fund stoppt), 2. Kern-Selbsttest
(Register, Konstanten gegen Register, jede Formel mit Standardargumenten endlich; ein Fehler stoppt), 3. Abschnitte
aus den Registern, ohne feste Liste: Prüfzahlen, Gesetze je Klasse, Gegenproben, Themen, Beschlüsse, Katalog.
Alle Aufgaben laufen in parallelen Prozessen (fork, ein Prozess je Kern; die Prüfzahlen in Scheiben von 8 Stück in
Registerfolge, die längsten Aufgaben des vorigen Laufs zuerst). Zwischenspeicher (.zwischenspeicher/, nicht im Repo):
je Aufgabe ein Ergebnis unter einem Schlüssel aus dem Hash des ganzen Codes der Simulation, der Datenhashes, des
Modelltexts und der Versionen; ändert sich irgendetwas davon, rechnet alles frisch.
Ausgabe: bericht.md (Kopf mit Commit, Parametersätzen, Laufzeit; je Abschnitt eine Tafel mit [G]- und [W]-Spalte; Lage der
Wächter; Änderungen seit dem letzten Lauf mit Beschluss-Prüfung) und register.json (Stand dieses Laufs, Grundlage des
nächsten Vergleichs). Laufzeit über 120 s gibt eine Warnung.
"""
import concurrent.futures as cf
import datetime
import glob
import hashlib
import inspect
import json
import math
import multiprocessing
import os
import pickle
import subprocess
import sys
import time

SIM = os.path.dirname(os.path.abspath(__file__))
if SIM not in sys.path:
    sys.path.insert(0, SIM)

import laden  # noqa: E402

BERICHT = os.path.join(SIM, "bericht.md")
REGISTER_JSON = os.path.join(SIM, "register.json")
BUDGET_S = 120
ZWISCHENSPEICHER = os.path.join(SIM, ".zwischenspeicher")
SCHEIBE = 8            # Prüfzahlen je Aufgabe (Registerfolge, damit Zwischenergebnisse im Prozess wiederverwendet werden)


# ------------------------------------------------------------------ Hilfen
def git(*args):
    return subprocess.run(["git", *args], cwd=SIM, capture_output=True, text=True).stdout.strip()


def zahl(x):
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return f"{x:.6g}"
    return str(x)


def g_text(g):
    return f"[G: {', '.join(g)}]" if g else ""


def w_text(w):
    return f"[W: {', '.join(w)}]" if w else ""


def gw_text(v):
    """Spalte [G]/[W] einer Prüfzahl oder eines Gesetzes."""
    return " ".join(t for t in (g_text(v.get("g")), w_text(v.get("w"))) if t)


def tafel(kopf, zeilen):
    if not zeilen:
        return ["(noch leer)", ""]
    aus = ["| " + " | ".join(kopf) + " |", "|" + "---|" * len(kopf)]
    aus += ["| " + " | ".join(str(z).replace("|", "¦") for z in zeile) + " |" for zeile in zeilen]
    return aus + [""]


# ------------------------------------------------------------------ 2 Kern-Selbsttest
def kern_selbsttest():
    import kern as K
    from kern.eingabe import EINGABEN
    from kern.formeln import FORMELN
    from kern.herkunft_g import formel_g, messwerte_g
    from kern.messwert import MESSWERTE
    from kern.parametersatz import PARAMETERSAETZE
    from kern.setzung import SETZUNGEN
    from kern.symbole import SYMBOLE
    fehler = []
    for name in ("G", "m_p", "m_e", "alpha", "r_e", "a_B", "R_p", "GM_sun", "R_sun", "L_sun", "TEFF_sun", "GM_earth",
                 "R_earth", "H0_KMS", "OMEGA_B_H2", "OMEGA_C_H2", "T_CMB", "N_EFF_STD", "KREISPAAR_GRENZE", "A0_LIT",
                 "A0_SPARC", "OMEGA_R"):
        if getattr(K, name) != MESSWERTE[name].wert:
            fehler.append(f"Konstante {name} ≠ Messwert")
    if K.A0 != SETZUNGEN["A0"].wert or K.g_earth != SETZUNGEN["G_ERDE"].wert:
        fehler.append("Konstante ≠ Setzung")
    ohne_arg = 0
    for k, f in FORMELN.items():
        sig = inspect.signature(f.funktion)
        if all(p.default is not inspect.Parameter.empty for p in sig.parameters.values()):
            ohne_arg += 1
            try:
                w = f.funktion()
                werte = w.values() if isinstance(w, dict) else (w if isinstance(w, tuple) else (w,))
                if not all(callable(x) or math.isfinite(float(x)) for x in werte):
                    fehler.append(f"{k}(): nicht endlich")
            except Exception as e:   # noqa: BLE001
                fehler.append(f"{k}(): {type(e).__name__}: {e}")
    status = {}
    for f in FORMELN.values():
        status[f.status] = status.get(f.status, 0) + 1
    zeilen = [("Symbole", len(SYMBOLE)), ("Formeln", ", ".join(f"{s} {n}" for s, n in sorted(status.items()))),
              ("Formeln mit [G] im Eingang", sum(1 for k in FORMELN if formel_g(k))),
              ("Messwerte", f"{len(MESSWERTE)} (davon [G] {len(messwerte_g())}; Konstanten CODATA und IAU zählen nicht)"),
              ("Setzungen", len(SETZUNGEN)), ("Eingaben", len(EINGABEN)), ("Parametersätze", len(PARAMETERSAETZE)),
              ("Formeln mit Standardargumenten ausgewertet", ohne_arg)]
    return fehler, zeilen


# ------------------------------------------------------------------ 1 und 3: Aufgaben in parallelen Prozessen
def _pytest(args):
    """pytest in eigenem Prozess; gibt (Rückgabecode, letzte Zeilen der Ausgabe)."""
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *args], cwd=SIM,
                       capture_output=True, text=True)
    return r.returncode, "\n".join(r.stdout.strip().splitlines()[-25:])


def _pruefzahl(k):
    from pruefung.pruefzahl import PRUEFZAHLEN
    p = PRUEFZAHLEN[k]
    try:
        r = p.pruefe()
    except Exception as e:   # noqa: BLE001
        r = dict(ergebnis="ABWEICHUNG", wert=None, soll=None, abweichung=None, toleranz=None, ok=False, g=[], w=[],
                 formeln=[], fehler=f"{type(e).__name__}: {e}")
    return dict(r, baustein=p.baustein, einheit=p.einheit or "", modell_text=p.modell_text, satz=p.satz)


def _gesetz(k):
    from pruefung.gesetz import GESETZE
    g = GESETZE[k]
    try:
        r = g.pruefe()
    except Exception as e:   # noqa: BLE001
        r = dict(ergebnis="WIDERSPRUCH", modell=None, vergleich=None, mass=f"{type(e).__name__}: {e}", g=[], w=[],
                 formeln=[])
    return dict(r, klasse=g.klasse, aussage=g.aussage, quelle=g.quelle, ebene=g.ebene,
                modell=_float(r["modell"]), vergleich=_float(r["vergleich"]))


def _gegenprobe(k):
    from pruefung.gegenprobe import GEGENPROBEN
    g = GEGENPROBEN[k]
    return dict(g.pruefe(), prueft=g.prueft, art=g.art)


def _leichte_waechter():
    """Alle Wächter außer den schweren und die Rückrichtung; dazu Themen, Beschlüsse und Katalog für den Bericht."""
    import waechter as W
    from forschung.beschluss import BESCHLUESSE
    from forschung.thema import THEMEN
    from katalog_eintrag import KATALOG
    return dict(funde=W.alle(ausser=W.SCHWER), rueckrichtung=W.rueckrichtung_bilanz(),
                themen={k: dict(frage=t.frage, stand=t.stand, versuche={v.kennung: v.urteil for v in t.versuche})
                        for k, t in THEMEN.items()},
                beschluesse={k: dict(aussage=b.aussage, ort=list(b.modell_ort), datum=b.datum, richtung=b.richtung)
                             for k, b in BESCHLUESSE.items()},
                katalog={k: dict(klasse=e.klasse, ziel=e.ziel) for k, e in KATALOG.items()})


_AUFGABEN = []   # wird vor dem fork gefüllt; die Prozesse erben die Liste und alle geladenen Register


def _fuehre_aus(nr):
    art, name, arg = _AUFGABEN[nr]
    t0 = time.time()
    if art == "pytest":
        r = _pytest(arg)
    elif art == "pruefzahlen":
        r = {k: _pruefzahl(k) for k in arg}
    elif art == "gesetze":
        r = {k: _gesetz(k) for k in arg}
    elif art == "gegenproben":
        r = {k: _gegenprobe(k) for k in arg}
    elif art == "leicht":
        r = _leichte_waechter()
    elif art == "doppelformeln":
        import waechter as W
        r = W.doppelformeln()
    elif art == "selbsttest":
        r = kern_selbsttest()
    else:
        raise ValueError(art)
    return nr, r, time.time() - t0


def schluessel(wurzel=SIM, modell=os.path.join(SIM, "..", "..", "modell")):
    """Schlüssel des Zwischenspeichers: Hash des ganzen Codes der Simulation (ohne Bericht, Register, Caches), der
    Datenhashes, des Modelltexts und der Versionen. Ändert sich eines davon, ist alles neu."""
    import numpy
    import scipy
    h = hashlib.sha256()
    for p in sorted(glob.glob(os.path.join(wurzel, "**", "*"), recursive=True)):
        rel = os.path.relpath(p, wurzel)
        if (not os.path.isfile(p) or rel.startswith((".zwischenspeicher", "__pycache__")) or "__pycache__" in rel
                or rel in ("bericht.md", "register.json") or rel.endswith(".pyc")):
            continue
        h.update(rel.encode() + b"\0" + open(p, "rb").read())
    for p in sorted(glob.glob(os.path.join(modell, "*.md"))):
        h.update(os.path.basename(p).encode() + b"\0" + open(p, "rb").read())
    h.update(f"{sys.version_info[:3]} {numpy.__version__} {scipy.__version__}".encode())
    return h.hexdigest()[:16]


def _aufgaben_plan():
    """Liste (art, name, argument) aller Aufgaben; Prüfzahlen in Scheiben in Registerfolge."""
    from pruefung.gegenprobe import GEGENPROBEN
    from pruefung.gesetz import GESETZE
    from pruefung.pruefzahl import PRUEFZAHLEN
    plan = [("pytest", "pytest schnell", ["-m", "not schwer", os.path.join(SIM, "tests")]),
            ("pytest", "pytest Gesetze", ["-m", "schwer", os.path.join(SIM, "tests", "test_gesetze.py")]),
            ("leicht", "Wächter leicht, Themen, Katalog", None),
            ("doppelformeln", "Wächter Doppelformeln", None),
            ("selbsttest", "Kern-Selbsttest", None)]
    ks = list(PRUEFZAHLEN)
    plan += [("pruefzahlen", f"Prüfzahlen {i // SCHEIBE + 1:02d}", ks[i:i + SCHEIBE]) for i in range(0, len(ks), SCHEIBE)]
    for kl in "ABCZME":
        ge = [k for k, g in GESETZE.items() if g.klasse == kl]
        if kl in "ZM":
            plan += [("gesetze", f"Gesetze {kl} {n + 1}", ge[n::2]) for n in range(2)]
        else:
            plan.append(("gesetze", f"Gesetze {kl}", ge))
    gp = list(GEGENPROBEN)
    plan += [("gegenproben", f"Gegenproben {n + 1}", gp[n::4]) for n in range(4)]
    return plan


def fuehre_aufgaben_aus(frisch=False):
    """Alle Aufgaben, parallel, mit Zwischenspeicher (frisch=True liest ihn nicht). Gibt (Ergebnisse je Aufgabennummer,
    Zeiten, Statistik)."""
    global _AUFGABEN
    _AUFGABEN = _aufgaben_plan()
    key = schluessel()
    ordner = os.path.join(ZWISCHENSPEICHER, key)
    os.makedirs(ordner, exist_ok=True)
    zeitdatei = os.path.join(ZWISCHENSPEICHER, "zeiten.json")
    try:
        with open(zeitdatei, encoding="utf-8") as f:
            alt_zeit = json.load(f)
    except (OSError, ValueError):
        alt_zeit = {}
    erg, zeiten, aus_speicher = {}, {}, 0
    offen = []
    for nr, (art, name, _) in enumerate(_AUFGABEN):
        datei = os.path.join(ordner, f"{nr:03d}.pkl")
        if os.path.exists(datei) and not frisch:
            try:
                with open(datei, "rb") as f:
                    erg[nr], zeiten[nr] = pickle.load(f)
                aus_speicher += 1
                continue
            except Exception:   # noqa: BLE001  (beschädigter Eintrag: frisch rechnen)
                pass
        offen.append(nr)
    offen.sort(key=lambda n: -alt_zeit.get(_AUFGABEN[n][1], 5.0))          # längste zuerst
    kerne = max(1, min(os.cpu_count() or 1, 8))
    if offen:
        with cf.ProcessPoolExecutor(max_workers=kerne, mp_context=multiprocessing.get_context("fork")) as ex:
            for fut in cf.as_completed([ex.submit(_fuehre_aus, n) for n in offen]):
                nr, r, dt = fut.result()
                erg[nr], zeiten[nr] = r, dt
                if not (_AUFGABEN[nr][0] == "pytest" and r[0] != 0):       # rote Wächter werden nie gespeichert
                    with open(os.path.join(ordner, f"{nr:03d}.pkl"), "wb") as f:
                        pickle.dump((r, dt), f)
    alt_zeit.update({_AUFGABEN[n][1]: round(zeiten[n], 1) for n in offen})
    with open(zeitdatei, "w", encoding="utf-8") as f:
        json.dump(alt_zeit, f, ensure_ascii=False, indent=0)
    return erg, zeiten, dict(aufgaben=len(_AUFGABEN), gespeichert=aus_speicher, kerne=kerne, schluessel=key)


def sammle(erg):
    """Ergebnisse der Aufgaben zu Wächter-Lage, Kern-Selbsttest und den Abschnitten des Berichts."""
    import waechter as W
    code, letzte = 0, []
    abschnitt = {"pruefzahlen": {}, "gesetze": {}, "gegenproben": {}}
    leicht = doppel = selbst = None
    for nr, r in erg.items():
        art, name, _ = _AUFGABEN[nr]
        if art == "pytest":
            code = max(code, r[0])
            if r[0]:
                letzte.append(f"{name}:\n{r[1]}")
        elif art in abschnitt:
            abschnitt[art].update(r)
        elif art == "leicht":
            leicht = r
        elif art == "doppelformeln":
            doppel = r
        elif art == "selbsttest":
            selbst = r
    from pruefung.gegenprobe import GEGENPROBEN
    from pruefung.gesetz import GESETZE
    from pruefung.pruefzahl import PRUEFZAHLEN
    ordnung = dict(pruefzahlen=PRUEFZAHLEN, gesetze=GESETZE, gegenproben=GEGENPROBEN)
    erg_ab = {a: {k: abschnitt[a][k] for k in ordnung[a]} for a in abschnitt}
    funde = dict(leicht["funde"])
    funde["gegenproben"] = W.gegenproben(ergebnisse=erg_ab["gegenproben"])
    funde["gegenproben_abdeckung"] = W.gegenproben_abdeckung(ergebnisse=erg_ab["gegenproben"])
    abdeckung = W.gegenproben_abdeckung_zahlen(ergebnisse=erg_ab["gegenproben"])
    funde["doppelformeln"] = doppel
    reihe = [n for n in W.STRENG] + sorted(W.ZAEHLEND)
    funde = {n: funde[n] for n in reihe}
    erg_ab.update(themen=leicht["themen"], beschluesse=leicht["beschluesse"], katalog=leicht["katalog"],
                  rueckrichtung=leicht["rueckrichtung"], abdeckung=abdeckung)
    return code, letzte, funde, selbst, erg_ab


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def stand_register():
    from kern.eingabe import EINGABEN
    from kern.formeln import FORMELN
    from kern.herkunft_g import formel_g
    from kern.messwert import MESSWERTE
    from kern.setzung import SETZUNGEN
    from kern.symbole import SYMBOLE
    return dict(
        symbole={z: dict(code=s.code_name, bedeutung=s.bedeutung, einheit=s.einheit_si, ebene=s.ebene,
                         zustand=s.zustand, herkunft=s.herkunft) for z, s in SYMBOLE.items()},
        formeln={k: dict(gleichung=f.gleichung, aus=list(f.aus), ein=list(f.ein), ebene=f.ebene, status=f.status,
                         modell_ort=list(f.modell_ort) if f.modell_ort else None, g=sorted(formel_g(k)))
                 for k, f in FORMELN.items()},
        messwerte={k: dict(wert=m.wert, fehler=m.fehler, einheit=m.einheit, art=m.art) for k, m in MESSWERTE.items()},
        setzungen={k: dict(wert=s.wert, aussage=s.aussage) for k, s in SETZUNGEN.items()},
        eingaben={k: dict(wert=e.wert, einheit=e.einheit) for k, e in EINGABEN.items()})


# ------------------------------------------------------------------ Änderungen seit dem letzten Lauf
def _flach(d, praefix=""):
    aus = {}
    for k, v in d.items():
        if isinstance(v, dict) and k not in ("fehler",):
            aus.update(_flach(v, f"{praefix}{k}."))
        else:
            aus[f"{praefix}{k}"] = v
    return aus


VERGLEICH = ("pruefzahlen", "gesetze", "gegenproben", "messwerte", "setzungen", "eingaben", "formeln")


def aenderungen(alt, neu):
    zeilen = []
    for teil in VERGLEICH:
        a, n = alt.get(teil, {}), neu.get(teil, {})
        for k in sorted(set(a) | set(n)):
            if k not in a:
                zeilen.append((teil, k, "neu", ""))
            elif k not in n:
                zeilen.append((teil, k, "entfernt", ""))
            else:
                fa, fn = _flach(a[k]) if isinstance(a[k], dict) else {"": a[k]}, \
                    _flach(n[k]) if isinstance(n[k], dict) else {"": n[k]}
                diff = [f"{f}: {zahl(fa.get(f))} → {zahl(fn.get(f))}" for f in sorted(set(fa) | set(fn))
                        if f not in ("formeln",) and fa.get(f) != fn.get(f)]
                if diff:
                    zeilen.append((teil, k, "geändert", "; ".join(diff)[:300]))
    return zeilen


def modell_diff(alter_commit, alte_zeit):
    """(geänderte Dateien in modell/ seit dem Commit des letzten Laufs, Dateien ohne Beschluss, diff --stat)."""
    from forschung.beschluss import BESCHLUESSE
    if not alter_commit:
        return [], [], ""
    dateien = [os.path.basename(p) for p in git("diff", "--name-only", alter_commit, "--", "../../modell").split()]
    stat = git("diff", "--stat", alter_commit, "--", "../../modell")
    datum = (alte_zeit or "")[:10]
    ohne = [d for d in dateien if not any(b.modell_ort[0] == d and b.datum >= datum for b in BESCHLUESSE.values())]
    return dateien, ohne, stat


# ------------------------------------------------------------------ Bericht
def ergebnis_text(v):
    """OK, ABWEICHUNG (mit Fehler) oder OFFEN (fehlt Formel: … / nur Forschung: …)."""
    from pruefung.pruefzahl import OFFEN
    if v["ergebnis"] == "OFFEN":
        return f"OFFEN ({OFFEN[v['offen']]}: {v['beschreibung']})"
    return v["ergebnis"] + (f" ({v['fehler']})" if v.get("fehler") else "")


def herkunft_zeilen(erg):
    """Gesetze je Ergebnis; Zahl der [G]- und [W]-Markierungen (Prüfzahlen, Gesetze) mit den Wurzeln von [W]."""
    from kern.herkunft_g import messwerte_g
    ges = erg["gesetze"]
    art = {e: sum(v["ergebnis"] == e for v in ges.values()) for e in ("OK", "SPANNUNG", "WIDERSPRUCH", "HINWEIS")}
    z = ["Gesetze: " + ", ".join(f"{e} {n}" for e, n in art.items()) + f" (zusammen {len(ges)}).", ""]
    pz, gs = erg["pruefzahlen"], erg["gesetze"]
    wurzeln = sorted({w for v in list(pz.values()) + list(gs.values()) for w in v.get("w", [])})
    zeilen = [("Prüfzahlen", len(pz), sum(bool(v["g"]) for v in pz.values()), sum(bool(v.get("w")) for v in pz.values())),
              ("Gesetze", len(gs), sum(bool(v["g"]) for v in gs.values()), sum(bool(v.get("w")) for v in gs.values()))]
    zeilen.append(("alle", *(sum(r[i] for r in zeilen) for i in (1, 2, 3))))
    z += ["Herkunftsmarken ([G] nur aus Messwert-Einträgen, ohne Konstanten; [W] aus gewählten Eingaben, Lesarten und "
          f"Zuordnungen; {len(messwerte_g())} Messwerte tragen [G]):", ""]
    z += tafel(("Art", "Zahl", "mit [G]", "mit [W]"), zeilen)
    z += ["Wurzeln von [W]: " + (", ".join(wurzeln) if wurzeln else "keine") + ".", ""]
    return z


def zusammenfassung(erg, funde):
    """Zählung der Prüfzahlen je Baustein (OK, ABWEICHUNG, OFFEN je Art), Rückrichtung (erfasst, gedeckt) und
    zählende Wächter."""
    import waechter as W
    z = ["## Zusammenfassung", ""]
    zeilen, summe = [], [0] * 6
    for b in sorted({v["baustein"] for v in erg["pruefzahlen"].values()}):
        vs = [v for v in erg["pruefzahlen"].values() if v["baustein"] == b]
        n = [len(vs), sum(v["ergebnis"] == "OK" for v in vs), sum(v["ergebnis"] == "ABWEICHUNG" for v in vs),
             sum(v["ergebnis"] == "OFFEN" for v in vs),
             sum(v.get("offen") == "fehlt_formel" for v in vs), sum(v.get("offen") == "nur_forschung" for v in vs)]
        summe = [a + x for a, x in zip(summe, n)]
        zeilen.append((b, *n))
    zeilen.append(("alle", *summe))
    z += tafel(("Baustein", "Prüfzahlen", "OK", "ABWEICHUNG", "OFFEN", "davon fehlt Formel", "davon nur Forschung"),
               zeilen)
    rr = erg.get("rueckrichtung", {})
    rz = [(d[:2], v["zahlen"], v["erfasst"], v["gedeckt"], v["offen"], v["zahlen"] - v["erfasst"])
          for d, v in rr.items()]
    s = [sum(r[i] for r in rz) for i in range(1, 6)]
    rz.append(("alle", *s))
    z += ["Rückrichtung je Fundstelle (erfasst = Prüfzahl, auch offen, Messwert, Eingabe, Setzung, Ausnahme oder "
          "Wiederholung; gedeckt = erfasst ohne offene Prüfzahlen):", ""]
    z += tafel(("Datei", "Zahlen", "erfasst", "gedeckt", "offen", "nicht erfasst"), rz)
    z += herkunft_zeilen(erg)
    z += ["Zählende Wächter: " + ", ".join(
        f"{n} {sum(len(v) for v in f.values()) if isinstance(f, dict) else len(f)}"
        for n, f in funde.items() if n in W.ZAEHLEND) + ".", ""]
    return z


def bericht(kopf, waechter, selbst, erg, aend, mdiff):
    from kern.parametersatz import PARAMETERSAETZE
    import waechter as W
    z = ["# Bericht der FTM-Simulation", "",
         f"Commit {kopf['commit']}{' (mit nicht committeten Änderungen)' if kopf['dirty'] else ''} · "
         f"{kopf['zeit']} · Laufzeit {kopf['laufzeit']:.1f} s"
         + (f" · **Warnung: {W.laufzeit_warnung(kopf['laufzeit'], BUDGET_S)}**" if kopf["laufzeit"] > BUDGET_S else ""), "",
         f"Aufgaben {kopf['statistik']['aufgaben']} in {kopf['statistik']['kerne']} parallelen Prozessen, davon "
         f"{kopf['statistik']['gespeichert']} aus dem Zwischenspeicher (Schlüssel {kopf['statistik']['schluessel']}); "
         f"Rechenzeit aller Aufgaben {kopf['statistik']['rechenzeit']:.0f} s; "
         "längste frisch gerechnete Aufgaben: " + ", ".join(f"{n} {t:.0f} s" for t, n in kopf["lang"]) + ".", "",
         "Parametersätze:", ""]
    z += [f"- {p.text()}" for p in PARAMETERSAETZE.values()] + [""]
    z += ["[G] = hängt an einem Messwert aus dem Gedächtnis (nicht nachgeprüft); Naturkonstanten (CODATA, IAU) tragen "
          "kein [G]. [W] = hängt an einer Wahl, die nicht im Modelltext steht (gewählte Eingabe, Lesart, Zuordnung, die "
          "der Nutzer noch bestätigen muss). Beide werden vererbt; die Spalte nennt die Wurzeln.", ""]
    # Wächter
    z += ["## Wächter", "", f"pytest: {'grün' if waechter[0] == 0 else 'ROT (Lauf gestoppt)'}", ""]
    for t in kopf["pytest_rot"]:
        z += ["```", t, "```", ""]
    zeilen = []
    for n, f in waechter[1].items():
        anzahl = sum(len(v) for v in f.values()) if isinstance(f, dict) else len(f)
        art = "zählend" if n in W.ZAEHLEND else "streng"
        beispiel = ""
        if isinstance(f, dict):
            beispiel = ", ".join(f"{d[:2]}: {len(v)}" for d, v in f.items())
        elif f:
            beispiel = "; ".join(f[:3])[:200]
        zeilen.append((n, art, anzahl, beispiel))
    z += tafel(("Wächter", "Art", "Funde", "je Datei / erste Funde"), zeilen)
    # Selbsttest
    z += ["## Kern-Selbsttest", ""] + tafel(("Register", "Stand"), selbst[1])
    if selbst[0]:
        z += ["Fehler:", ""] + [f"- {x}" for x in selbst[0]] + [""]
    if erg is None:
        return z
    z += zusammenfassung(erg, waechter[1])
    # Prüfzahlen je Baustein
    z += ["## Prüfzahlen", "", "OFFEN: die Kernformel fehlt (fehlt Formel) oder die Zahl beruht nur auf einer "
          "Forschungsrechnung (nur Forschung); weder OK noch ABWEICHUNG.", ""]
    for b in sorted({v["baustein"] for v in erg["pruefzahlen"].values()}):
        z += [f"### Baustein {b}", ""]
        rows = [(k, v["modell_text"], zahl(v["wert"]) if v["wert"] is not None else "", v["einheit"],
                 zahl(v["abweichung"]) if v["abweichung"] is not None else "",
                 zahl(v["toleranz"]) if v["toleranz"] is not None else "", v["satz"], ergebnis_text(v),
                 gw_text(v))
                for k, v in erg["pruefzahlen"].items() if v["baustein"] == b]
        z += tafel(("Prüfzahl", "Modell", "Rechnung", "Einheit", "Abweichung", "Toleranz", "Satz", "Ergebnis", "[G] / [W]"),
                   rows)
    # Gesetze je Klasse
    z += ["## Gesetze", ""]
    for kl in "ABCZME":
        rows = [(k, v["aussage"], zahl(v["modell"]), zahl(v["vergleich"]), v["mass"], v["ergebnis"], gw_text(v))
                for k, v in erg["gesetze"].items() if v["klasse"] == kl]
        z += [f"### Klasse {kl}", ""] + tafel(("Gesetz", "Aussage", "Modell", "Vergleich", "Maß", "Ergebnis", "[G] / [W]"),
                                              rows)
    z += ["## Gegenproben", ""]
    a = erg.get("abdeckung")
    if a:
        z += [f"Wächter gegenproben_abdeckung: {a['ohne']} von {a['formeln']} Modellformeln ohne unabhängige Probe "
              f"(davon {a['nur_umkehr']} nur mit Umkehrprobe, {a['teilende']} nur mit Proben, die eine Formel mit der "
              f"geprüften teilen); Umkehrproben (art \"umkehr\") und teilende Proben zählen nicht.", ""]
    z += tafel(("Gegenprobe", "prüft", "Art", "größter Abstand", "ruft sich selbst", "Ergebnis"),
               [(k, v["prueft"], v["art"], f"{v['groesster_abstand']:.1e}", "ja" if v["ruft_selbst"] else "nein",
                 "OK" if v["ok"] else "ABWEICHUNG") for k, v in erg["gegenproben"].items()])
    z += ["## Themen", ""]
    z += tafel(("Thema", "Frage", "Stand", "Versuche (Urteil)"),
               [(k, v["frage"], v["stand"], "; ".join(f"{a}: {u}" for a, u in v["versuche"].items()))
                for k, v in erg["themen"].items()])
    z += ["## Beschlüsse", ""]
    z += tafel(("Beschluss", "Aussage", "Ort", "Datum", "Richtung"),
               [(k, v["aussage"], f"{v['ort'][0]}: {v['ort'][1]}", v["datum"], v["richtung"])
                for k, v in erg["beschluesse"].items()])
    z += ["## Katalog", ""]
    z += tafel(("Skript", "Klasse", "Ziel"), [(k, v["klasse"], v["ziel"] or "") for k, v in erg["katalog"].items()])
    # Änderungen
    z += ["## Änderungen seit dem letzten Lauf", ""]
    if kopf["voriger"] is None:
        z += ["Erster Lauf: kein Vergleich.", ""]
    else:
        z += [f"Voriger Lauf: Commit {kopf['voriger']['commit']}, {kopf['voriger']['zeit']}.", ""]
        z += ["### Werte", ""] + (tafel(("Teil", "Kennung", "Art", "Änderung"), aend) if aend else ["keine", ""])
        dateien, ohne, stat = mdiff
        z += ["### Modelltext (git diff modell/)", ""]
        z += ([f"```\n{stat}\n```", ""] if stat else ["keine Änderung", ""])
        if dateien:
            z += [f"Ohne Beschluss seit dem letzten Lauf: {', '.join(ohne) if ohne else 'keine'}", ""]
    return z


def main():
    t0 = time.time()
    voriger = None
    if os.path.exists(REGISTER_JSON):
        with open(REGISTER_JSON, encoding="utf-8") as f:
            voriger = json.load(f)
    laden.lade_alles()
    erg_a, zeiten, statistik = fuehre_aufgaben_aus(frisch="--frisch" in sys.argv)
    import waechter as W
    code, letzte, funde, selbst, erg = sammle(erg_a)
    streng_rot = code != 0 or any(funde[n] for n in W.STRENG)
    if streng_rot or selbst[0]:
        erg = None
    stand = stand_register()
    kopf = dict(commit=git("rev-parse", "--short", "HEAD"), dirty=any(not z.endswith(("bericht.md", "register.json"))
                              for z in git("status", "--porcelain").splitlines()),
                zeit=datetime.datetime.now().isoformat(timespec="seconds"), laufzeit=0.0, statistik=dict(statistik, rechenzeit=sum(zeiten.values())),
                lang=sorted(((zeiten[n], _AUFGABEN[n][1]) for n in zeiten), reverse=True)[:5], pytest_rot=letzte,
                voriger=None if voriger is None else dict(commit=voriger.get("commit"), zeit=voriger.get("zeit")))
    neu = dict(commit=kopf["commit"], zeit=kopf["zeit"], **stand, **(erg or {}),
               waechter={n: (sum(len(v) for v in f.values()) if isinstance(f, dict) else len(f))
                         for n, f in funde.items()})
    neu = json.loads(json.dumps(neu, ensure_ascii=False, default=str))   # gleiche Form wie der gelesene Vorlauf
    aend = aenderungen(voriger or {}, neu) if voriger else []
    mdiff = modell_diff(voriger.get("commit") if voriger else None, voriger.get("zeit") if voriger else None)
    kopf["laufzeit"] = time.time() - t0
    neu["laufzeit_s"] = round(kopf["laufzeit"], 1)
    with open(BERICHT, "w", encoding="utf-8") as f:
        f.write("\n".join(bericht(kopf, (code, funde), selbst, erg, aend, mdiff)).rstrip() + "\n")
    if erg is not None:
        with open(REGISTER_JSON, "w", encoding="utf-8") as f:
            json.dump(neu, f, ensure_ascii=False, indent=1, default=str)
            f.write("\n")
    ok = erg is not None
    if ok:
        import doku_erzeugen
        doku_erzeugen.main()
    print(f"lauf.py: {'fertig' if ok else 'gestoppt'} in {kopf['laufzeit']:.1f} s → {os.path.relpath(BERICHT)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
