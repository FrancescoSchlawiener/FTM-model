# -*- coding: utf-8 -*-
"""Ein Lauf der Simulation: alle Szenarien, Prüfung der Modellzahlen, ein nummerierter Bericht.

Aufruf: python simulation/lauf.py, mit --voll ohne den Speicher der Prüfwerte
Schreibt berichte/bericht_NNNN.md (nie überschrieben) und berichte/bericht_NNNN/ mit den Bahndaten
(je Szenario eine CSV-Datei und szene.json) für eine spätere Darstellung.
"""
from functools import lru_cache
import datetime
import json
import math
import os
import re
import subprocess
import sys
import time

HIER = os.path.dirname(os.path.abspath(__file__))
WURZEL = os.path.dirname(HIER)
sys.path.insert(0, HIER)

import kern as K            # noqa: E402
import szenarien            # noqa: E402
import daten as D           # noqa: E402
import turm as T            # noqa: E402
import beispiele as B       # noqa: E402
import kerne as KE          # noqa: E402

R_ERDE, R_SONNE = K.R_ERDE, K.R_SONNE

# ---------------------------------------------------------------- Speicher der Prüfwerte
# Jeder Wert der Prüfliste, jede Bahn und die Elemente werden nur neu gerechnet, wenn sich ihr Ausdruck, der Syntaxbaum einer
# Funktion oder Größe, die sie (auch mittelbar) benutzen, oder die gelesenen Daten geändert haben. --voll rechnet alles neu.
import ast            # noqa: E402
import hashlib        # noqa: E402

SPEICHER = os.path.join(HIER, ".pruefspeicher", "werte.json")


@lru_cache(None)
def _module():
    """Je Modul der Simulation: seine Definitionen auf oberster Ebene (Name → Knoten) und seine Kürzel für andere Module
    der Simulation (import kern as K, from kern import a0)."""
    aus = {}
    for f in os.listdir(HIER):
        if not f.endswith(".py"):
            continue
        name = f[:-3]
        baum = ast.parse(open(os.path.join(HIER, f), encoding="utf-8").read())
        defs, kurz = {}, {}
        for k in baum.body:
            if isinstance(k, (ast.FunctionDef, ast.ClassDef)):
                defs[k.name] = k
            elif isinstance(k, (ast.Assign, ast.AnnAssign)):
                for z in (k.targets if isinstance(k, ast.Assign) else [k.target]):
                    for n in ast.walk(z):
                        if isinstance(n, ast.Name):
                            defs[n.id] = k
        for k in ast.walk(baum):
            if isinstance(k, ast.Import):
                for a in k.names:
                    kurz[a.asname or a.name] = (a.name, None)
            elif isinstance(k, ast.ImportFrom) and k.module:
                for a in k.names:
                    kurz[a.asname or a.name] = (k.module, a.name)
        aus[name] = (defs, kurz)
    return aus


@lru_cache(None)
def _datenstand(ordner):
    """Stand eines Ordners, aus dem gelesen wird: Pfad und Inhalt jeder Datei darin. Der Hash des Inhalts wird je Datei
    nur neu gebildet, wenn sich Größe oder Zeit geändert haben (.pruefspeicher/dateien.json)."""
    pfad = os.path.join(os.path.dirname(SPEICHER), "dateien.json")
    try:
        with open(pfad, encoding="utf-8") as f:
            alt = json.load(f)
    except (OSError, ValueError):
        alt = {}
    h = hashlib.sha256()
    for wurzel, _, dateien in sorted(os.walk(ordner)):
        for f in sorted(dateien):
            if f.endswith(".teil"):                              # halb geholte Dateien zählen nicht
                continue
            datei = os.path.join(wurzel, f)
            s = os.stat(datei)
            merk = alt.get(datei)
            if not merk or merk[:2] != [s.st_size, s.st_mtime_ns]:
                inhalt = hashlib.sha1()
                with open(datei, "rb") as d:
                    for stueck in iter(lambda: d.read(1 << 20), b""):
                        inhalt.update(stueck)
                merk = alt[datei] = [s.st_size, s.st_mtime_ns, inhalt.hexdigest()]
            h.update(f"{os.path.relpath(datei, ordner)}:{merk[2]};".encode())
    os.makedirs(os.path.dirname(pfad), exist_ok=True)
    with open(pfad, "w", encoding="utf-8") as f:
        json.dump(alt, f)
    return h.hexdigest()


def _bezuege(knoten, modul):
    """Namen der Simulation, die ein Knoten benutzt: (Modul, Name)."""
    module = _module()
    defs, kurz = module[modul]
    aus = set()
    for n in ast.walk(knoten):
        if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name) and n.value.id in kurz:
            ziel, einzel = kurz[n.value.id]
            if einzel is None and ziel in module:
                aus.add((ziel, n.attr))
        elif isinstance(n, ast.Name):
            if n.id in defs:
                aus.add((modul, n.id))
            elif n.id in kurz and kurz[n.id][1] and kurz[n.id][0] in module:
                aus.add((kurz[n.id][0], kurz[n.id][1]))
    return aus


@lru_cache(None)
def _bezuege_def(m, n):
    """Namen der Simulation, die die Definition n im Modul m benutzt."""
    return frozenset(_bezuege(_module()[m][0][n], m))


@lru_cache(None)
def _def_hash(m, n):
    """Hash der Definition n im Modul m, aus ihrem Syntaxbaum (Kommentare und Form des Textes zählen nicht)."""
    return hashlib.sha256(ast.dump(_module()[m][0][n]).encode()).hexdigest()


@lru_cache(None)
def _datenordner(gesehen):
    """Die Ordner der Messdaten, die die Definitionen gesehen lesen: jeder Unterordner von archiv/weltmodell/daten, dessen
    Name als Text in einer Definition steht, die ORT benutzt. Nennt eine solche Definition keinen, zählt der ganze Ordner."""
    namen = set(os.listdir(D.ORT))
    aus = set()
    for m, n in gesehen:
        knoten = _module()[m][0][n]
        if not any(isinstance(k, ast.Name) and k.id == "ORT" and isinstance(k.ctx, ast.Load)
                   or isinstance(k, ast.Attribute) and k.attr == "ORT" for k in ast.walk(knoten)):
            continue
        hier = {k.value for k in ast.walk(knoten) if isinstance(k, ast.Constant) and k.value in namen}
        aus |= {os.path.join(D.ORT, o) for o in hier} if hier else {D.ORT}
    return tuple(sorted(aus))


def _schluessel(ausdruck=None, start=()):
    """Hash eines Prüfwerts: sein Ausdruck (oder die Definitionen start, (Modul, Name)) und alle Definitionen, die er
    mittelbar benutzt, mit dem Stand der Daten."""
    module = _module()
    offen, gesehen = list(start) + (list(_bezuege(ausdruck, "lauf")) if ausdruck is not None else []), set()
    while offen:
        m, n = offen.pop()
        if (m, n) in gesehen or n not in module[m][0]:
            continue
        gesehen.add((m, n))
        offen += _bezuege_def(m, n)
    h = hashlib.sha256(ast.dump(ausdruck).encode() if ausdruck is not None else b"")
    for m, n in sorted(gesehen):
        h.update(f"\0{m}.{n}\0{_def_hash(m, n)}".encode())
    ordner = _datenordner(frozenset(gesehen))                     # Messdaten: nur die gelesenen Unterordner
    for o in ordner:
        h.update(f"{o}:{_datenstand(o)};".encode())
    if any(m == "huelle_ablage" for m, _ in gesehen):            # Zwischenspeicher der Hülle
        h.update(_datenstand(os.path.join(HIER, "zwischen")).encode())
    return h.hexdigest()


def _gespeichert(was, start, rechne, voll=False):
    """Ergebnis von rechne(), abgelegt in .pruefspeicher/was.pkl unter dem Schlüssel der Definitionen start (mittelbar)
    und dem Stand der Daten. Gerechnet wird nur, wenn sich der Schlüssel geändert hat oder voll."""
    import pickle
    schluessel = _schluessel(start=start)
    pfad = os.path.join(os.path.dirname(SPEICHER), was + ".pkl")
    if not voll and os.path.exists(pfad):
        with open(pfad, "rb") as f:
            alt, erg = pickle.load(f)
        if alt == schluessel:
            return erg
    erg = rechne()
    os.makedirs(os.path.dirname(pfad), exist_ok=True)
    with open(pfad, "wb") as f:
        pickle.dump((schluessel, erg), f)
    return erg


def _pruefliste():
    """Die Prüfliste ohne sie auszurechnen: je Eintrag (Baustein, Stelle, Ausdruck, Schlüssel, Soll)."""
    quelle = open(os.path.abspath(__file__), encoding="utf-8").read()
    fn = next(k for k in ast.parse(quelle).body if isinstance(k, ast.FunctionDef) and k.name == "_modellzahlen")
    liste = next(k for k in fn.body if isinstance(k, ast.Return)).value
    aus = []
    for t in liste.elts:
        nr, stelle, wert, soll = t.elts
        aus.append((nr.value, stelle.value, wert, _schluessel(wert), eval(compile(ast.Expression(soll), "", "eval"))))
    return aus

# Zahlen des Modelltexts: (Datei, Textstelle wie im Modell, Wert aus dem Kern, Stellen im Text)
@lru_cache(None)
def _ght():
    return __import__("huelle_4_formen").gehaltener_trichter()


@lru_cache(None)
def _kug():
    return __import__("huelle_4_formen").kugelform()


def _om_or():
    H0, Ode = K.atmen_heute()
    o_r = K.dichte_heute()[0] / (K.rho_krit100 * (H0 / 100)**2)
    return 1 - Ode - o_r, o_r


def _om_or_a0(a0):
    H0, Ode = K.atmen_heute(a0)
    o_r = K.dichte_heute()[0] / (K.rho_krit100 * (H0 / 100)**2)
    return 1 - Ode - o_r, o_r


@lru_cache(maxsize=None)
def _om_or_q(Q):
    def f(a0):
        rho_r, rho_m = K.dichte_heute()[0], K.omega_b * (1 + Q) * K.rho_krit100
        rho_c = rho_r + rho_m + K.dichte_dunkle_energie(a0, K.B1_ZELLE)
        return rho_m / rho_c, rho_r / rho_c
    return f


def _h0_q(Q, a0):
    rho = K.dichte_heute()[0] + K.omega_b * (1 + Q) * K.rho_krit100 + K.dichte_dunkle_energie(a0, K.B1_ZELLE)
    return math.sqrt(8 * math.pi * K.G / 3 * rho) * K.Mpc / 1e3


def _v_mw(r_kpc):
    """Milchstraße als Punktmasse (6·10¹⁰ M☉) nach der Regel: v [km/s] bei r."""
    return K.bahntempo_regel(6e10 * K.M_sun, r_kpc * K.kpc) / 1e3


def _abst():
    return D.abstaende_pruefung(*_om_or())


@lru_cache(maxsize=None)
def _bao_a0(Q=None):
    """a₀ aus der absoluten Skala von DESI DR2 (r_d des Kerns, H₀ und Ω_m aus dem Atmen): Bestwert, Δχ² = 1, χ²_min.
    Q = None: ω_c gemessen; sonst ω_c = Q·ω_b (dann nur ω_b, T_CMB und a₀)."""
    import numpy as np
    from scipy.optimize import brentq
    gr = np.linspace(1.06, 1.14, 41)
    if Q is None:
        cs = np.array([D.bao_absolut(K.atmen_heute(a * 1e-10)[0], *_om_or_a0(a * 1e-10), K.schallhorizont()) for a in gr])
    else:
        cs = np.array([D.bao_absolut(_h0_q(Q, a * 1e-10), *_om_or_q(Q)(a * 1e-10), K.schallhorizont(Q)) for a in gr])
    i = int(np.argmin(cs))
    p = np.polyfit(gr[i - 5:i + 6], cs[i - 5:i + 6], 2)
    b = -p[1] / (2 * p[0])
    m = np.polyval(p, b)
    f = lambda a: np.polyval(p, a) - m - 1
    return float(b), brentq(f, b - 0.05, b), brentq(f, b, b + 0.05), float(m)


def _ganzes_chi(a0, ob, Q, spitzen, w0=-1.0, wa=0.0):
    """χ² des Ganzen (daten.geometrie_reste) mit ω_b nach Planck, spitzen: dazu ω_c = Q·ω_b gegen die Höhe der Spitzen."""
    r = D.geometrie_reste(ob, a0, Q, w0, wa, "spitzen", spitzen)
    return float(r @ r)


def _s2():
    """Masse in Sgr A* aus der Bahn von S2 (Halbachse am Himmel mal Abstand, Umlaufzeit) und die Halbachse [kg, m]."""
    a = K.S2["a"] * K.S2["R0"] * K.AU
    return K.masse_aus_bahn(a, K.S2["P"] * K.yr), a


@lru_cache(maxsize=None)
def _naht_sparc():
    """SPARC jenseits der Naht (a_N < (−ln q)²·a₀, q = Q/(1 + Q), Q gemessen): Zahl der Punkte, mittlerer Rest gegen die
    Regel und sein Fehler, mittlerer Rest gegen die bei 1 + Q gedeckelte Verstärkung [dex]."""
    import numpy as np
    yd = D.sparc_eichung(K.kpc)
    a, _ = D._bestes(D.sparc_profil("Bose", K.kpc, yd).sum(0))
    x, r = D.sparc_reste(K.kpc, a, yd)
    m = x < math.log10(math.log((1 + K.Q_kosm) / K.Q_kosm)**2)
    nu = 1 / (-np.expm1(-np.sqrt(10**x[m])))
    gedeckelt = r[m] + np.log10(nu) - np.log10(np.minimum(nu, 1 + K.Q_kosm))
    return int(m.sum()), float(r[m].mean()), float(r[m].std() / math.sqrt(m.sum())), float(gedeckelt.mean())


@lru_cache(maxsize=None)
def _sparc_standard():
    """Standardmodell an SPARC: je Galaxie die Faltung der Baryonen plus ein NFW-Halo (zwei freie Zahlen) gegen die eine
    Schwelle der Regel, dieselben Störgrößen; χ², freie Zahlen, BIC und der größte gestapelte Rest der Regel je halber
    Dekade der Faltung der Baryonen [dex]."""
    import numpy as np
    yd = D.sparc_eichung(K.kpc)
    C = D.sparc_profil("Bose", K.kpc, yd)
    a, chi_r = D._bestes(C.sum(0))
    H = D.sparc_halo(K.kpc, K.H0_kms, yd, True)
    chi_h = float(H[:, 0].sum())
    am_rand = int(((H[:, 1] <= 0.51) | (H[:, 1] >= 2.99) | (H[:, 2] <= 0.01) | (H[:, 2] >= 1.99)).sum())
    mb = np.log10([1.33 * p["MHI"] + yd * p["L36"] for _, p, _, _, _ in D.sparc_auswahl(yd)])
    m200 = 3 * H[:, 1]                                       # lg M₂₀₀ bis auf eine Konstante (M₂₀₀ ∝ V₂₀₀³)
    innen = (H[:, 1] > 0.51) & (H[:, 1] < 2.99) & (H[:, 2] > 0.01) & (H[:, 2] < 1.99)
    gerade = np.polyfit(mb[innen], m200[innen], 1)
    masse_streuung = float(np.std(m200[innen] - np.polyval(gerade, mb[innen])))
    n = sum(len(g[2]["R"]) for g in D.sparc_auswahl(yd))
    x, r = D.sparc_reste(K.kpc, a, yd)
    rand = np.arange(-2.0, 1.51, 0.5)
    stapel = max(abs(r[(x >= lo) & (x < hi)].mean()) for lo, hi in zip(rand[:-1], rand[1:]) if ((x >= lo) & (x < hi)).sum() > 10)
    return dict(chi_regel=chi_r, chi_halo=chi_h, zahlen=2 * len(C), bic_regel=chi_r + math.log(n),
                bic_halo=chi_h + 2 * len(C) * math.log(n), stapel=stapel, galaxien=len(C), punkte=n, rand=am_rand,
                masse_streuung=masse_streuung)


@lru_cache(maxsize=None)
def _sparc_bahn_feld():
    """Bahnen gegen Feldform an SPARC (Υ_d aus der Eichung): χ² der Regel ν·g_N und der wirbelfreien Feldform als Bahnregel
    je im besten a₀, Gang der Reste von ν·g_N gegen die Feldform (Steigung, Streuung, Galaxien), dazu Feldform gegen ν·g_N
    bei a₀ (Linse gegen Dynamik) innerhalb 2R_d (Median aller, der flächenschwachen und der flächenhellen Hälfte nach der
    hellsten Flächenhelligkeit der Scheibe) und außerhalb 4R_d (Median) [dex], ebenso die Stufe der Bahn knapp über der
    Schicht gegen ihre Mitte, lg[ν(√(g_N² + (2πGΣ)²)/a₀)/ν(g_N/a₀)] mit Σ = Υ_d·SBdisk + Σ_Gas (stufe_...), und das Mittel
    über die Schicht, ν gleichverteilt über den Anteil m der Masse unter der Höhe, √(g_N² + (2πGΣ·m)²) (schicht_...)."""
    import numpy as np
    yd = D.sparc_eichung(K.kpc)
    a, chi_v = D._bestes(D.sparc_profil("Bose", K.kpc, yd).sum(0))
    chi_w = D._bestes(D.sparc_profil("Bose", K.kpc, yd, None, True).sum(0))[1]
    st, ds, n = D.sparc_wirbelfrei_gang(K.kpc, a, yd)
    G, inn, aus, sb, s_inn, s_aus, m_inn = D.sparc_galaxien(), [], [], [], [], [], []
    nu = lambda y: 1 / -np.expm1(-np.sqrt(y))
    for name, p, k, _, _ in D.sparc_auswahl(yd, None, True):
        w, R, Rd = D._wirbelfrei_an(k["wf"], yd, K.a0), k["R"], G[name]["Rd"]
        gb = np.maximum(k["Vgas"] * abs(k["Vgas"]) + yd * k["Vdisk"]**2 + D.YB * k["Vbul"]**2, 1e-6) * 1e6 / (R * K.kpc)
        gz = 2 * math.pi * K.G * K.M_sun / K.pc**2 * (yd * k["SBdisk"] + np.interp(R, *[D.sparc_flaeche(name)[x] for x in ("R", "Sg")]))
        stufe = np.log10(nu(np.hypot(gb, gz) / K.a0) / nu(gb / K.a0))
        m = (np.arange(200) + 0.5) / 200
        mitte = np.log10(nu(np.hypot(gb[:, None], gz[:, None] * m[None, :]) / K.a0).mean(1) / nu(gb / K.a0))
        if (R <= 2 * Rd).sum() >= 2:
            inn.append(w[R <= 2 * Rd].mean()); sb.append(k["SBdisk"].max()); s_inn.append(stufe[R <= 2 * Rd].mean())
            m_inn.append(mitte[R <= 2 * Rd].mean())
        if (R >= 4 * Rd).sum() >= 2:
            aus.append(w[R >= 4 * Rd].mean()); s_aus.append(stufe[R >= 4 * Rd].mean())
    inn, sb, s_inn, m_inn = np.array(inn), np.array(sb), np.array(s_inn), np.array(m_inn)
    lo = sb < np.median(sb)
    return dict(chi_v=chi_v, chi_w=chi_w, steigung=st, streuung=ds, galaxien=n, innen=float(np.median(inn)),
                lsb=float(np.median(inn[lo])), hsb=float(np.median(inn[~lo])), aussen=float(np.median(aus)),
                stufe_innen=float(np.median(s_inn)), stufe_lsb=float(np.median(s_inn[lo])),
                stufe_hsb=float(np.median(s_inn[~lo])), stufe_aussen=float(np.median(s_aus)),
                schicht_innen=float(np.median(m_inn)), schicht_lsb=float(np.median(m_inn[lo])),
                schicht_hsb=float(np.median(m_inn[~lo])))


@lru_cache(maxsize=None)
def _achsen_zweiter_weg():
    """Feldform in der Ebene auf zwei Wegen: Kuzmin-Scheibe im Löser (welt.achsen_feld) gegen ν am Betrag GM/(R² + a²) an
    der Oberfläche, größte Abweichung [dex] bei R = 0,5, 2, 8 kpc; NGC 2403 (Υ_d = 0,67) im Löser gegen das Gitter in drei
    Richtungen (welt.gitter_feld, 128³ und 256³ Zellen, Kante 2,5·R_max, Körper aus derselben Flächendichte, Fehler
    quadratisch in der Zelle hochgerechnet), Mittel über das mittlere Drittel der Punkte [dex]; Bahn knapp über der Schicht
    einer Scheibe an der Schwelle im Löser (Höhe 5h und 10h, auf 0 gezogen) gegen die Oberfläche nach Freeman mit 2πGΣ,
    größte Abweichung bei 0,5 bis 3·R_d [dex]."""
    import numpy as np
    import welt as W
    a, M = 2 * K.kpc, 3e10 * K.M_sun
    R = np.array([0.5, 2.0, 8.0]) * K.kpc
    S = lambda r: a * M / (2 * math.pi * (r**2 + a**2)**1.5)
    gN, _, w = W.achsen_feld([(S, 0.005 * K.kpc)], None, R, schritt=(0.05 * K.kpc, 0.005 * K.kpc / 8), stufe=1.02)
    nu = lambda y: 1 / -np.expm1(-np.sqrt(y))
    kuzmin = float(np.max(np.abs(np.log10(w[0] / (nu(K.G * M / (R**2 + a**2) / K.a0) * K.G * M * R / (R**2 + a**2)**1.5)))))
    n, yd = "NGC2403", 0.67
    k, f, Rd = D.sparc_kurve(n), D.sparc_flaeche(n), D.sparc_galaxien()[n]["Rd"]
    Rk = k["R"]
    m = slice(len(Rk) // 3, 2 * (len(Rk) // 3))
    loeser = float(D.sparc_wirbelfrei(n, yd, (K.a0,))[0][m].mean())
    Rg = np.linspace(0, 1.5 * Rk[-1], 3000)
    S = (yd * np.where(Rg <= Rk[-1], np.interp(Rg, Rk, k["SBdisk"]), k["SBdisk"][-1] * np.exp(-(Rg - Rk[-1]) / Rd))
         + np.interp(Rg, Rk, f["Sg"], right=0))
    dm = S * 2 * np.pi * Rg * np.gradient(Rg)
    rng, gitter = np.random.default_rng(1), []
    for N, nk in ((128, 200000), (256, 400000)):
        Rs = np.abs(rng.choice(Rg, nk, p=dm / dm.sum()) + rng.uniform(-0.5, 0.5, nk) * (Rg[1] - Rg[0]))
        ph, z = rng.uniform(0, 2 * np.pi, nk), rng.laplace(0, 0.15 * Rd, nk)
        x = np.concatenate([np.stack([Rs * np.cos(ph), Rs * np.sin(ph), z], 1), np.stack([Rk, 0 * Rk, 0 * Rk], 1)]) * K.kpc
        mm = np.concatenate([np.full(nk, dm.sum() * 1e6 * K.M_sun / nk), np.zeros(len(Rk))])
        L = 2.5 * Rk[-1] * K.kpc
        gitter.append(float(np.log10(W.gitter_feld(x, mm, N, L, "wirbelfrei")[nk:, 0]
                                     / W.gitter_feld(x, mm, N, L, "vektor")[nk:, 0])[m].mean()))
    Ms = 6e10 * K.M_sun
    Rd = math.sqrt(K.G * Ms / K.a0 / 2)
    S, h = (lambda r: Ms / (2 * math.pi * Rd**2) * np.exp(-r / Rd)), 0.002 * Rd
    R = np.array([0.5, 1.0, 2.0, 3.0]) * Rd
    gR, gz = W.achsen_feld([(S, h)], None, R, schritt=(0.02 * Rd, h / 8), stufe=1.02, hoehen=np.array([0, 5, 10]) * h)
    b = np.log10(nu(np.hypot(gR, gz) / K.a0) * gR)
    from scipy.special import i0, i1, k0, k1
    y = R / (2 * Rd)
    gF = 4 * math.pi * K.G * Ms / (2 * math.pi * Rd**2) * Rd * y**2 * (i0(y) * k0(y) - i1(y) * k1(y)) / R
    flaeche = np.log10(nu(np.hypot(gF, 2 * math.pi * K.G * S(R)) / K.a0) / nu(gF / K.a0))
    stufe = float(np.max(np.abs((2 * (b[1] - b[0]) - (b[2] - b[0])) - flaeche)))
    return kuzmin, loeser, (4 * gitter[1] - gitter[0]) / 3, stufe


@lru_cache(maxsize=None)
def _ganzes_zustand(Q, z_max=None):
    """Zustandsteil des Ganzen bei festem Q (mit den Spitzen): bester Gang w₀, w_a über a₀ und ω_b, Δχ² gegen den Formteil
    allein; mit z_max ist der Gang so gebunden, dass die Dunkle Energie dort ihr Maximum hat (w = −1)."""
    from scipy.optimize import minimize
    if z_max is None:
        f = lambda p: _ganzes_chi(p[0], p[1], Q, True, p[2], p[3])
        x0 = [1.07, K.omega_b, -0.85, -0.5]
    else:
        u = z_max / (1 + z_max)
        f = lambda p: _ganzes_chi(p[0], p[1], Q, True, p[2], (-1 - p[2]) / u)
        x0 = [1.07, K.omega_b, -0.86]
    r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-6, maxiter=8000, maxfev=8000))
    w0 = float(r.x[2])
    wa = float(r.x[3]) if z_max is None else (-1 - w0) / u
    return _ganzes_bestes(Q, True)[0] - float(r.fun), w0, wa, float(r.fun)


@lru_cache(maxsize=None)
def _ganzes_maximum():
    """Rotverschiebung des Maximums der Dunklen Energie bei Q = 2π − 1 (w = −1 im Gang w₀ + w_a·z/(1 + z)), Bestwert und
    Δχ² = 1-Grenzen aus dem Profil; dazu im Bestwert des Ganzen die Gleichheit von Formteil und Faltung und der Beginn der
    beschleunigten Dehnung (ρ_DE = ρ_m/2, Strahlung vernachlässigt)."""
    from scipy.optimize import brentq
    Q = 2 * math.pi - 1
    _, w0, wa, m = _ganzes_zustand(Q)
    zb = K.maximum_dunkle_energie(w0, wa)
    prof = lambda z: _ganzes_zustand(Q, round(z, 4))[3] - m - 1
    _, a, ob = _ganzes_bestes(Q, True)
    alt = K.omega_b
    K.omega_b = ob
    try:
        rho_m = K.dichte_heute(Q)[2]
    finally:
        K.omega_b = alt
    z_eq = (K.dichte_dunkle_energie(a * 1e-10, K.B1_ZELLE) / rho_m)**(1 / 3) - 1
    z_besch = (2 * K.dichte_dunkle_energie(a * 1e-10, K.B1_ZELLE) / rho_m)**(1 / 3) - 1
    return zb, brentq(prof, 0.15, zb, xtol=1e-3), brentq(prof, zb, 1.0, xtol=1e-3), z_eq, z_besch


@lru_cache(maxsize=None)
def _kids():
    """Linse um isolierte Galaxien (KiDS-1000): χ² der Regel ohne äußeres Feld, mit 10⁻⁴·a₀ im Betrag, jenseits der Naht
    auf 1 + Q begrenzt."""
    return D.kids_chi2()[0], D.kids_chi2(1e-4 * K.a0)[0], D.kids_chi2(deckel=True)[0]


@lru_cache(maxsize=None)
def _zwerge(gruppe, mit_feld):
    return D.zwerge_reste(gruppe, mit_feld)


@lru_cache(maxsize=None)
def _zwerge_m(gruppe, tief=False):
    """Reste der Zwerge mit dem Gas, das die Metallizität ihrer Sterne liest (daten.zwerg_gas, geschichte="metall")."""
    return D.zwerge_reste(gruppe, False, geschichte="metall", tief=tief)


@lru_cache(maxsize=None)
def _zwerge_g(gruppe, mit_feld=False, tief=False):
    """Reste der Zwerge mit dem Gas, das sie einmal trugen (daten.zwerg_gas)."""
    return D.zwerge_reste(gruppe, mit_feld, geschichte=True, tief=tief)


def _xcop(groesse):
    """Median, kleinster und größter Wert einer Größe der X-COP-Haufen (daten.xcop_bilanz)."""
    import numpy as np
    w = [v[groesse] for v in D.xcop_bilanz().values()]
    return float(np.median(w)), float(min(w)), float(max(w))


@lru_cache(maxsize=None)
def _cmb():
    """Hintergrundstrahlung und übrige Präzisionsdaten, ΛCDM gegen das Modell (ω_c = (2π − 1)·ω_b, H₀ aus a₀): je
    Datenstand die besten Punkte aus der Ablage; dazu der Gesamtwert am besten Punkt des Modells an allen Daten
    nachgerechnet; mit Zustandsteil (w₀, w_a) das Modell, das Modell bei a₀ = 1,06 und ΛCDM, das Maximum der Dunklen Energie
    mit Fehler aus der Kovarianz, Gleichheit und Beginn der beschleunigten Dehnung mit dem Formteil."""
    import numpy as np
    Q = K.q_je_radiant()
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    aus = {n: (D.planck_anpassung(None, m), D.planck_anpassung(Q, m))
           for n, m in (("spektren", ("spektren",)), ("linse", ("spektren", "linse")), ("alle", alle))}
    aus["fest"] = D.planck_anpassung(Q, alle, fest=("a0",))
    r = D.planck_reste(aus["alle"][1]["punkt"], Q)
    aus["nach"] = sum(float(r[k] @ r[k]) for k in alle + ("vor",))
    aus["zustand"] = z = D.planck_anpassung(Q, alle, zustand=True)
    aus["zustand_fest"] = D.planck_anpassung(Q, alle, fest=("a0",), zustand=True)
    aus["zustand_lcdm"] = D.planck_anpassung(None, alle, zustand=True)
    def mit_fehler(e, f):                                    # Wert und Fehler aus der Kovarianz von (w₀, w_a)
        p, n = e["punkt"], e["kovarianz"]["namen"]
        C = np.array(e["kovarianz"]["werte"])[np.ix_([n.index("w0"), n.index("wa")], [n.index("w0"), n.index("wa")])]
        x, h = f(p["w0"], p["wa"]), 1e-6
        g = np.array([f(p["w0"] + h, p["wa"]) - x, f(p["w0"], p["wa"] + h) - x]) / h
        return x, math.sqrt(g @ C @ g)
    hoehe = lambda w0, wa: K.gang_dunkle_energie(w0, wa, K.maximum_dunkle_energie(w0, wa)) - 1
    aus["z_max"] = mit_fehler(z, K.maximum_dunkle_energie)
    aus["hoehe"], aus["hoehe_fest"] = mit_fehler(z, hoehe), mit_fehler(aus["zustand_fest"], hoehe)
    p = z["punkt"]
    rho_m = (p["wb"] + p["wc"] + 0.06 / 93.14) * K.rho_krit100      # Formteil gegen Faltung, Strahlung vernachlässigt
    rho_de = K.dichte_dunkle_energie(p["a0"] * 1e-10)
    aus["z_eq"], aus["z_besch"] = (rho_de / rho_m)**(1 / 3) - 1, (2 * rho_de / rho_m)**(1 / 3) - 1
    aus["saett"] = s = D.planck_anpassung(Q, alle, zustand="saettigung", start={"wa": -0.56})
    aus["saett_fest"] = D.planck_anpassung(Q, alle, fest=("a0",), zustand="saettigung", start={"wa": -0.56})
    aus["z_saett"] = K.maximum_dunkle_energie(s["punkt"]["w0"], s["punkt"]["wa"])
    def hoehe_saett(pp):                                     # w₀ folgt aus w_a, a₀ und ω_b (Maximum an der Sättigung)
        a_s = (((1 + Q) * pp["wb"] + 0.06 / 93.14) * K.rho_krit100 / K.dichte_dunkle_energie(pp["a0"] * 1e-10))**(1 / 3)
        w0 = -1 - pp["wa"] * (1 - a_s)
        return hoehe(w0, pp["wa"])
    n, C = s["kovarianz"]["namen"], np.array(s["kovarianz"]["werte"])
    h0 = hoehe_saett(s["punkt"])
    g = np.array([(hoehe_saett(dict(s["punkt"], **{k: s["punkt"][k] + 1e-6 * max(1, abs(s["punkt"][k]))})) - h0)
                  / (1e-6 * max(1, abs(s["punkt"][k]))) for k in n])
    aus["hoehe_saett"] = h0, math.sqrt(g @ C @ g)
    aus["q_lesen"] = D.planck_anpassung("lesen", alle, zustand="saettigung", start={"wa": -0.56})
    a0f = round(K.a0 * 1e10, 4)                              # Formteil mit dem Standardwert, Sättigung bei ihm
    aus["saett_form"] = sf = D.planck_anpassung(Q, alle, zustand="saettigung", start={"wa": -0.56}, form=a0f)
    je = K.zustandsteil_je_gas(Q) / 3                        # Zustandsteil je Formteil und je Anteil F des Gases
    def zustand_form(pp):                                    # Anteil heute, F heute, F am Maximum, z des Maximums, Höhe
        a_s = (((1 + Q) * pp["wb"] + 0.06 / 93.14) * K.rho_krit100 / K.dichte_dunkle_energie(a0f * 1e-10))**(1 / 3)
        w0 = -1 - pp["wa"] * (1 - a_s)
        x, h = (pp["a0"] / a0f)**2 - 1, hoehe(w0, pp["wa"])
        return x, x / je, (x + h * (1 + x)) / je, K.maximum_dunkle_energie(w0, pp["wa"]), h
    n, C = sf["kovarianz"]["namen"], np.array(sf["kovarianz"]["werte"])
    w = zustand_form(sf["punkt"])
    st = {k: 1e-6 * max(1, abs(sf["punkt"][k])) for k in n}
    G = np.array([[(zustand_form(dict(sf["punkt"], **{k: sf["punkt"][k] + st[k]}))[i] - w[i]) / st[k] for k in n]
                  for i in range(5)])
    aus["zustand_form"] = [(w[i], math.sqrt(G[i] @ C @ G[i])) for i in range(5)]
    p = s["punkt"]                                           # Würfelzelle im besten Punkt: nur die Moden k = 2π·n/L
    L = K.zelle_laenge(a0_=p["a0"] * 1e-10, rho_m=(p["wb"] + p["wc"] + 0.06 / 93.14) * K.rho_krit100) / K.Mpc
    _, R, chi = K.daempfung_zelle(p["wb"], p["wc"], p["H0"], p["ns"], p["tau"], L, p["w0"], p["wa"])
    _, Rsw, _ = K.daempfung_zelle(0, 0, 0, 1.0, 0, L, lmax=3, kc_chi=150, sachs_wolfe=chi)   # zweiter Weg: nur Geometrie
    tt, te, ee = K.cmb_spektren(p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["tau"], w0=p["w0"], wa=p["wa"])
    r0 = D.planck_rest(tt, te, ee, p["A"])[:2]
    tt[:len(R)] *= R
    r1 = D.planck_rest(tt, te, ee, p["A"])[:2]
    aus["zelle"] = dict(L_d=L / (2 * chi), R2=float(R[0]), R3=float(R[1]), sw2=float(Rsw[0]), sw3=float(Rsw[1]),
                        dchi2=float(r0 @ r0 - r1 @ r1))
    return aus


@lru_cache(maxsize=None)
def _scherung():
    """KiDS-1000-Scherung im besten Punkt mit Sättigung (Modell) und für ΛCDM: χ² am Punkt, S8, frei gelesenes S8 und χ²."""
    Q = K.q_je_radiant()
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    aus = {}
    for n, e in (("modell", D.planck_anpassung(Q, alle, zustand="saettigung", start={"wa": -0.56})),
                 ("lcdm", D.planck_anpassung(None, alle))):
        p = e["punkt"]
        pt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p.get("w0", -1.0), p.get("wa", 0.0))
        aus[n] = D.kids_scherung_chi2(pt), D.kids_scherung_frei(pt)
    return aus


@lru_cache(maxsize=None)
def _haufen_einfall():
    """Konzentration der Mode aus dem Einfall: Läufe vieler Körper im besten Punkt mit Sättigung (welt.haufen_konzentration,
    R_L = 25 und 18 Mpc, Zufallsfelder 1 und 2). Gibt dict: je Lauf das Ergebnis, Mittel von c (geometrisch), kleinstes und
    größtes c, Q aus HIFLUGCS beim mittleren c mit Fehler (Auslassen je eines Haufens)."""
    import numpy as np
    import welt as W
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    laeufe = [W.haufen_konzentration(punkt, R, saat=s) for R in (25.0, 18.0) for s in (1, 2)]
    cs = np.array([l["c"] for l in laeufe])
    c = float(np.exp(np.log(cs).mean()))
    return dict(laeufe=laeufe, c=c, c_min=float(cs.min()), c_max=float(cs.max()),
                m_min=min(l["M200"] for l in laeufe), m_max=max(l["M200"] for l in laeufe),
                koerper=min(l["koerper"] for l in laeufe), q=D.haufen_q(c=c)[0], q_fehler=D.haufen_q_fehler(c))


@lru_cache(maxsize=None)
def _gas_stabil(weg):
    """Q_T, bei dem im Median sechs Zehntel des Gases der SPARC-Galaxien unter der Gezeitendichte liegen (daten.sparc_gas_gezeit)."""
    from scipy.optimize import brentq
    return brentq(lambda q: D.sparc_gas_gezeit(round(q, 4), weg)[2] - 0.6, 1.0, 2.0, xtol=1e-3)

@lru_cache(maxsize=None)
def _galaxie_drehung():
    """Verwindung auf der Skala der Milchstraße aus den Läufen vieler Körper (R_L = 1,32 Mpc, Zufallsfelder 1 bis 4,
    welt.haufen_drehung) und was die Scheiben an SPARC verlangen (daten.sparc_drehung): dict je Größe (kleinster, größter
    Wert): M200 [10¹¹ M☉], lam (Mode am Ende), lam_lin (bis zur Umkehr), dazu sparc (Zahl, Median, 16 %, 84 %)."""
    import welt as W
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    d = [W.haufen_drehung(punkt, 1.32, saat=s) for s in (1, 2, 3, 4)]
    li = [W.haufen_drehung(punkt, 1.32, saat=s, linear=True) for s in (1, 2, 3, 4)]
    mm = lambda xs: (min(xs), max(xs))
    return dict(M200=mm([x["M200"] / 1e11 for x in d]), lam=mm([x["lam"] for x in d]), lam_lin=mm([x["lam"] for x in li]),
                sparc=D.sparc_drehung(D.sparc_eichung(K.kpc), p["H0"], K.q_je_radiant()))

@lru_cache(maxsize=None)
def _scheibe_austausch():
    """Drehimpuls einer Scheibe an der Schwelle nach 3 Gyr (welt.scheibe_lauf, aus der Ablage): (kleinster und größter
    Verlust unter der Regel über die sieben Läufe, Abweichung nach Newton, Zahl der Läufe)."""
    import welt as W
    regel = [W.scheibe_lauf(saat=s, n=n)["jz"][-1] for n in (20000, 40000) for s in (7, 8, 9)] + [W.scheibe_lauf(saat=7, n=80000)["jz"][-1]]
    newton = W.scheibe_lauf(regel=False, saat=7)["jz"][-1]
    return 1 - max(regel), 1 - min(regel), abs(newton - 1), len(regel)

@lru_cache(maxsize=None)
def _scheibe_formen():
    """Halbmassenradius einer kalten Scheibe an der Schwelle nach 3 Gyr gegen den Anfang (welt.scheibe_lauf, Zufall 10):
    (mit der Regel als Vektor auf dem Gitter, mit der Regel als Vektor im Baum, wirbelfrei Anfang und Ende [r*], größte
    Abweichung des Drehimpulses wirbelfrei)."""
    import welt as W
    v_g, v_b, w = W.scheibe_lauf(feld="vektor", saat=10), W.scheibe_lauf(saat=10), W.scheibe_lauf(feld="wirbelfrei", saat=10)
    return (v_g["rh"][-1] / v_g["rh"][0], v_b["rh"][-1] / v_b["rh"][0], w["rh"][0], w["rh"][-1],
            max(abs(x - 1) for x in w["jz"]))

def _schicht_gleichverteilt():
    """Mittel über eine Schicht nach Masse für das isotherme Profil sech²(z/z₀) und das exponentielle e^{−|z|/h}, direkt
    über die Höhe integriert, gegen das gleichverteilte Mittel über den Anteil m der Masse unter der Höhe: größte
    Abweichung von ν(√(g² + u²)) mit u = 2πGΣ·m, für g = 0,1 bis 3·a₀ und 2πGΣ = 0,1 bis 3·a₀."""
    import numpy as np
    from scipy.integrate import quad
    nu = lambda y: 1 / -math.expm1(-math.sqrt(y))
    fehler = 0.0
    for g in (0.1, 1.0, 3.0):
        for U in (0.1, 1.0, 3.0):
            gleich = quad(lambda m: nu(math.hypot(g, U * m)), 0, 1)[0]
            sech = quad(lambda t: nu(math.hypot(g, U * math.tanh(t))) / math.cosh(t)**2, 0, 40)[0]
            expo = quad(lambda t: nu(math.hypot(g, U * (1 - math.exp(-t)))) * math.exp(-t), 0, 60)[0]
            fehler = max(fehler, abs(sech / gleich - 1), abs(expo / gleich - 1))
    return fehler


@lru_cache(maxsize=None)
def _breiten():
    """Breitenfunktion der Bereiche mit allen Baryonen als Galaxie (kern.breitenfunktion im Punkt der Sättigung, Linie der
    Scheibe an der Schwelle mit der gelesenen Streuung bis lg M_HI = 9,5) gegen
    ALFALFA α.100 (daten.alfalfa_hi_dichte, H₀ des Punkts) je Feld: Verhältnis bei lg W50 = 1,7 bis 1,8 und 2,55 bis 2,65
    (Frühling, Herbst), dazu die größte Abweichung des Anteils der Masse in Bereichen über 10¹⁰, 10¹², 10¹⁴ M☉ von
    erfc(δ_c/√(2S))."""
    import numpy as np
    from scipy.special import erfc
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    ka = np.array([1.7, 1.8, 2.55, 2.65])
    m = K.breitenfunktion(punkt, ka, sigma=_gas_streuung()["leicht"])
    v = [m[i] / D.alfalfa_hi_dichte(punkt[2], feld=f, breiten=ka)[0][i] for i in (0, 2) for f in ("fruehling", "herbst")]
    lg = np.linspace(6.0, 15.5, 381)
    dn = K.bereiche_massenfunktion(punkt, lg)
    rho = (punkt[0] + punkt[1]) * K.rho_krit100 / K.M_sun * K.Mpc**3
    k, P = K.spektrum_linear(punkt)
    sig = lambda M: math.sqrt(np.trapezoid(P * (lambda x: 3 * (np.sin(x) - x * np.cos(x)) / x**3)(k * (3 * M / (4 * math.pi * rho))**(1 / 3))**2 * k**3, np.log(k)) / (2 * math.pi**2))
    fehler = max(abs(np.trapezoid((10**lg * dn)[lg >= m0], lg[lg >= m0]) / rho - erfc(0.15 * (12 * math.pi)**(2 / 3) / sig(10**m0) / math.sqrt(2))) for m0 in (10.0, 12.0, 14.0))
    return tuple(v) + (fehler,)


@lru_cache(maxsize=None)
def _gas_streuung():
    """Streuung des Gases aus den zwei Breiten der Linie (daten.hi_streuung, über die Linie der Scheibe an der Schwelle und
    von Hand über den Rand einer Stufe) und der Anteil des Gases unter der Gezeitendichte in SPARC mit ihr
    (daten.sparc_gas_gezeit mit sigma): Median σ bis lg M_HI = 9,5 und bei 10,5 bis 11, je (Median, nach Masse) aus der
    Flächendichte, die übrigen Wege, die Spanne für ein Drittel und das Dreifache der Streuung, alle Anteile."""
    import numpy as np
    aus, alle, andere, spanne = {}, [], [], []
    for rand, z in ((False, ""), (True, "_rand")):
        lgM, _, sg = D.hi_streuung(rand=rand)
        aus["leicht" + z] = float(np.median(sg[lgM < 9.5]))
        aus["schwer" + z] = float(np.median(sg[(lgM >= 10.5) & (lgM < 11.0)]))
        f = D.hi_streuung_gang(rand=rand)[3]
        aus["flaeche" + z] = D.sparc_gas_gezeit(weg="flaeche", sigma=f)[2:]
        for weg in ("vgas", "hi"):
            andere += list(D.sparc_gas_gezeit(weg=weg, sigma=f)[2:])
        for k in (1 / 3, 3.0):
            spanne += list(D.sparc_gas_gezeit(weg="flaeche", sigma=lambda x, k=k: k * f(x))[2:])
        alle += list(aus["flaeche" + z])
    aus.update(andere=andere, spanne=spanne, alle=alle + andere)
    return aus


@lru_cache(maxsize=None)
def _fornax(weg="summe"):
    """Zeit [Gyr], in der GC1, GC2 und GC5 von Fornax (daten.fornax_haufen) mit der Reibung der Regel von ihrem Abstand,
    dem projizierten und dem 1,5-fachen, bis in den Kern R_e sinken. weg "summe": kern.sinkzeit, "bahn": dieselbe Bahn
    als Anfangswertaufgabe (solve_ivp bis r = R_e). Gibt {Name: (Zeit beim projizierten, beim 1,5-fachen Abstand)}."""
    from scipy.integrate import solve_ivp
    f = D.fornax_haufen()
    M, a, s = f["M"] * K.M_sun, f["R_e"] * K.kpc, f["sigma"] * 1e3
    aN = lambda r: K.G * M * r / (r * r + a * a)**1.5
    acc = lambda r: K.verstaerkung(aN(r)) * aN(r)
    rho = lambda r: 3 * M / (4 * math.pi * a**3) * (1 + r * r / (a * a))**-2.5

    def rate(t, y, m, rk):
        r = y[0]
        h, v = 1e-4 * r, math.sqrt(acc(r) * r)
        dL = (math.sqrt(acc(r + h) * (r + h)) * (r + h) - math.sqrt(acc(r - h) * (r - h)) * (r - h)) / (2 * h)
        ell = (K.G * m / max(-(acc(r + h) - acc(r - h)) / (2 * h), 1e-60))**(1 / 3)
        bmin = max(rk, K.G * m / v**2)
        return [-r * K.reibung(m, v, rho(r), s, bmin, max(r, 1.0001 * bmin), aN=aN(r), ell=ell) / dL]

    def bahn(m, r0, rk):
        ende = lambda t, y, *arg: y[0] - a
        ende.terminal = True
        return solve_ivp(rate, (0, 200 * K.Gyr), [r0], args=(m, rk), events=ende, rtol=1e-8,
                         max_step=0.05 * K.Gyr).t_events[0][0]
    aus = {}
    for h in f["haufen"]:
        if h["r_p"] < f["R_e"]:
            continue
        m, rk = h["m"] * K.M_sun, h["r_halb"] * K.pc
        zeit = (lambda r0: K.sinkzeit(m, r0, M, a, s, rk, r_ende=a)) if weg == "summe" else (lambda r0: bahn(m, r0, rk))
        aus[h["name"].split("-")[-1]] = tuple(zeit(fak * h["r_p"] * K.kpc) / K.Gyr for fak in (1.0, 1.5))
    return aus


@lru_cache(maxsize=None)
def _gezeit_satelliten(geschichte=False, weg="wende"):
    """Die drei Satelliten der Milchstraße mit dem kleinsten r_J/R_e im Perizentrum (daten.satelliten_gezeit): ihre r_J/R_e,
    das nächste der übrigen, ihre Ränge im Rest gegen die Regel und die Wahrscheinlichkeit, dass drei zufällige so weit
    vorn stehen wie der hinterste von ihnen, C(k, 3)/C(N, 3)."""
    a = D.satelliten_gezeit(weg, geschichte)
    rang = {r[2]: i + 1 for i, r in enumerate(a)}
    nach_rj = sorted(a, key=lambda r: r[0])
    raenge = sorted(rang[r[2]] for r in nach_rj[:3])
    return dict(zahl=len(a), eng=[r[0] for r in nach_rj[:3]], naechster=nach_rj[3][0], raenge=raenge,
                zufall=math.comb(raenge[-1], 3) / math.comb(len(a), 3))


@lru_cache(maxsize=None)
def _phasen_bahn(hell=-8.0, faktor=1.0, haufen=False, r_min=0.0):
    """Zweiter Weg zu den Phasen der Satelliten (daten.satelliten_phasen): die Bahn in der Zeit, r̈ = L²/r³ − ν·GM/r²
    mit M = faktor·kern.M_B_MILCHSTRASSE, als Anfangswertaufgabe (solve_ivp) bis zum Perizentrum (rückwärts, wenn
    v_r > 0) und von dort bis zum Apozentrum. haufen, r_min wie dort. Gibt {Name: Phase}."""
    from scipy.integrate import solve_ivp
    GM = faktor * K.G * K.M_B_MILCHSTRASSE * K.M_sun
    aus = {}
    for x in D.satelliten_schar("mw", haufen=haufen):
        if x["wirt"] != "mw" or not x["M_V"] < hell or x["r"] <= r_min:
            continue
        r, vr, L = x["r"] * K.kpc, x["v_r"] * 1e3, x["r"] * K.kpc * x["v_t"] * 1e3
        rate = lambda t, y: [y[1], L**2 / y[0]**3 - K.verstaerkung(GM / y[0]**2) * GM / y[0]**2]
        wende = lambda t, y: y[1]
        wende.terminal = True
        lauf = lambda y0, t: solve_ivp(rate, (0, t), y0, events=wende, rtol=1e-10, atol=1e-6,
                                       max_step=0.01 * K.Gyr).t_events[0]
        tp = abs(lauf([r, vr], -20 * K.Gyr if vr > 0 else 20 * K.Gyr)[0])
        y = solve_ivp(rate, (0, -tp if vr > 0 else tp), [r, vr], rtol=1e-10, atol=1e-6, max_step=0.01 * K.Gyr).y[:, -1]
        aus[x["name"]] = tp / lauf([y[0], abs(y[1]) + 1e-6], 20 * K.Gyr)[0]
    return aus


@lru_cache(maxsize=None)
def _lmc_gruppe_bahn(schritt=2e-4):
    """Zweiter Weg zur Gruppe der LMC (daten.lmc_gruppe_rand): dieselben Körper im Rahmen der Milchstraße, die LMC im
    verstärkten Feld der Milchstraße, die SMC in dem beider, die Zwerge in allen dreien, Milchstraße und LMC erfahren
    den Gegenzug zum Fall ihrer Satelliten, m_Sat/m_Wirt mal so groß (die Milchstraße als Rahmen: abgezogen), dazu die
    Rückstellung H_∞, im Bocksprung (Stoß, Drift, Stoß) mit festem Schritt [Gyr] 3 Gyr rückwärts. Gibt Zeit des letzten
    Apozentrums [Gyr], D und r_J des Paars [kpc] heute und dort, kleinstes und größtes D, wer dort und wer heute
    innerhalb liegt, die Zwerge der LMC dort (Zahl, von) und die übrigen (Zahl, von)."""
    import numpy as np
    s = D.satelliten_schar("mw", wolken=True)
    lmc = next(x for x in s if x["name"] == "lmc")
    rest = sorted((x for x in s if x["name"] != "lmc"), key=lambda x: x["name"] != "smc")
    eigen = np.array([x["wirt"] == "lmc" for x in rest])
    M = K.G * K.M_B_MILCHSTRASSE * K.M_sun
    mL, mS = K.G * D.satellit("lmc")["M_b"] * K.M_sun, K.G * D.satellit("smc")["M_b"] * K.M_sun
    H = K.groesste_umkehr(K.M_B_MILCHSTRASSE * K.M_sun)[1]
    nu = np.vectorize(K.verstaerkung)
    zug = lambda GM, d: nu(GM / np.sum(d * d, axis=-1, keepdims=True)) * GM * d / np.linalg.norm(d, axis=-1, keepdims=True)**3
    gegen = lambda GM, d, GW: GM / GW * zug(GW, d)               # Gegenzug des Wirts GW zum Fall des Satelliten GM
    X = np.array([lmc["x"]] + [x["x"] for x in rest]) * K.kpc
    V = np.array([lmc["v"]] + [x["v"] for x in rest]) * 1e3

    def beschl(X):
        a = zug(M, -X) + H**2 * X
        a[1:] += zug(mL, X[0] - X[1:])                           # SMC und Zwerge im Feld der LMC
        a[2:] += zug(mS, X[1] - X[2:])                           # Zwerge im Feld der SMC
        a[0] += gegen(mS, X[1] - X[0], mL)
        return a - gegen(mL, X[:1], M)[0] - gegen(mS, X[1:2], M)[0]
    dt = -schritt * K.Gyr
    a = beschl(X)
    bahn = [X.copy()]
    for i in range(int(round(3.0 / schritt))):
        V = V + 0.5 * dt * a
        X = X + dt * V
        a = beschl(X)
        V = V + 0.5 * dt * a
        bahn.append(X.copy())
    bahn = np.array(bahn) / K.kpc
    t = np.arange(len(bahn)) * dt / K.Gyr
    Dd = np.linalg.norm(bahn[:, 0], axis=1)
    rJ = Dd * ((mL + mS) / M)**0.25 / math.sqrt(2)
    k = next(j for j in range(1, len(Dd) - 1) if Dd[j - 1] < Dd[j] > Dd[j + 1])
    innen = lambda i: np.linalg.norm(bahn[i, 1:] - bahn[i, 0], axis=1) < rJ[i]
    return dict(t_apo=float(t[k]), D=(Dd[0], Dd[k]), r_J=(rJ[0], rJ[k]), D_grenzen=(float(Dd.min()), float(Dd.max())),
                dort=[rest[i]["name"] for i in np.where(innen(k))[0]], heute=[rest[i]["name"] for i in np.where(innen(0))[0]],
                eigen=(int(np.sum(innen(k) & eigen)), int(eigen.sum())), uebrige=(int(np.sum(innen(k) & ~eigen)), int((~eigen).sum())))


@lru_cache(maxsize=None)
def _satelliten_phasen():
    """Phasen der Bahnen der Satelliten der Milchstraße im Feld der Regel (daten.satelliten_phasen): Zahl n der klassischen
    (M_V < −8), ihr Mittel auf beiden Wegen mit dem Fehler 1/√(12n), die gelesenen Baryonen M mit Grenzen und die Mittel
    des zweiten Wegs dort (±1/√(12n) zurückgerechnet, je ½), wie viele mit (1 + Q)·Newton und mit Newton ungebunden
    wären, das Mittel der lichtschwachen mit Fehler, und dasselbe für die Kugelsternhaufen über 20 kpc und heller als
    M_V = −5 (haufen_…)."""
    k, a = D.satelliten_phasen(-8.0), D.satelliten_phasen(0.0)
    mittel = lambda d: sum(d.values()) / len(d)
    e = 1 / math.sqrt(12 * len(k))
    M = D.satelliten_phasen_masse()
    zweiter = [mittel(_phasen_bahn(faktor=m / K.M_B_MILCHSTRASSE)) + z * e for m, z in zip(M[:3], (0, 1, -1))]
    sats = [x for x in D.satelliten_schar("mw") if x["wirt"] == "mw" and x["M_V"] < -8]
    GM = K.G * K.M_B_MILCHSTRASSE * K.M_sun
    frei = lambda f: sum((x["v_r"]**2 + x["v_t"]**2) * 1e6 > 2 * f * GM / (x["r"] * K.kpc) for x in sats)
    schwach = [v for n, v in a.items() if n not in k and math.isfinite(v)]
    h = D.satelliten_phasen(-5.0, 1.0, True, 20.0)
    eh = 1 / math.sqrt(12 * len(h))
    Mh = D.satelliten_phasen_masse(-5.0, True, 20.0)
    zweiter_h = [mittel(_phasen_bahn(-5.0, m / K.M_B_MILCHSTRASSE, True, 20.0)) + z * eh for m, z in zip(Mh[:3], (0, 1, -1))]
    return dict(n=len(k), mittel=mittel(k), mittel_bahn=mittel(_phasen_bahn()), e=e, M=M, zweiter=zweiter,
                frei_q=frei(1 + K.Q_kosm), frei_newton=frei(1.0), schwach=sum(schwach) / len(schwach),
                e_schwach=1 / math.sqrt(12 * len(schwach)), haufen_n=len(h), haufen_mittel=mittel(h),
                haufen_mittel_bahn=mittel(_phasen_bahn(-5.0, 1.0, True, 20.0)), haufen_e=eh, haufen_M=Mh,
                haufen_zweiter=zweiter_h)


@lru_cache(maxsize=None)
def _formteil_gas():
    """Alle Präzisionsdaten mit F₀ aus dem Gas gehalten (daten.planck_anpassung, Gas im Kreis, fest F₀ beim kleinsten,
    mittleren und größten Wert aus _gas_raum), frei nur k von der Dunklen Energie, gegen ΛCDM: (χ² je F₀, χ² ΛCDM)."""
    import numpy as np
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    F = _gas_raum()[1]
    chi = [D.planck_anpassung(K.q_je_radiant(), alle, fest=("F0",), zustand="gas_bera",
                              start={"F0": float(x) / 100})["chi2"]["gesamt"] for x in (min(F), np.mean(F), max(F))]
    return chi, D.planck_anpassung(None, alle)["chi2"]["gesamt"]


@lru_cache(maxsize=None)
def _formteil_lesen():
    """Formteil a₀, frei aus allen Präzisionsdaten gelesen, mit F₀ aus dem Gas gehalten (daten.planck_anpassung, Lesart
    gezählt mit Startwert F0, beim kleinsten, mittleren und größten Wert aus _gas_raum): (Formteile, ihre Fehler,
    Formteil mit allem gezählten Gas)."""
    import numpy as np
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    F = _gas_raum()[1]
    r = [D.planck_anpassung(K.q_je_radiant(), alle, zustand="gezaehlt", start={"F0": float(x) / 100})
         for x in (min(F), np.mean(F), max(F))]
    return ([x["punkt"]["form"] for x in r], [x["fehler"]["a0"] for x in r],
            D.planck_anpassung(K.q_je_radiant(), alle, zustand="gezaehlt")["punkt"]["form"])


@lru_cache(maxsize=None)
def _warp_drehung():
    """Warp der Milchstraße unter LMC, SMC und Sgr (kern.warp_ringe, Bahnen gemeinsam): je Fassung (gekoppelt, ohne
    Kopplung, rund mit Newton gekoppelt, die LMC mit 25-fachem Antrieb) und Auflösung (N = 70 und 100) die Neigung bei 14 kpc [°], die Drehung der Knotenlinie [km/s/kpc] und ihr
    Abstand zur Messung von Gaia DR2 in systematischen Fehlern (daten.warp_gaia)."""
    s = lambda n, g=False: D.satellit(n, g)
    q = ((s("lmc")["M_b"] * K.M_sun, s("lmc")["x"], s("lmc")["v"]), (s("smc")["M_b"] * K.M_sun, s("smc")["x"], s("smc")["v"]),
         (s("sagittarius_1", True)["M_b"] * K.M_sun, s("sagittarius_1", True)["x"], s("sagittarius_1", True)["v"],
          s("sagittarius_1")["M_b"] * K.M_sun))
    r, _, e, _ = D.warp_gaia()
    aus = {}
    lmc25 = ((q[0][0], q[0][1], q[0][2], 25 * q[0][0]),) + q[1:]
    for name, qq, kw in (("70", q, {}), ("100", q, dict(N=100, n_phi=240)), ("frei70", q, dict(gekoppelt=False)),
                         ("frei100", q, dict(gekoppelt=False, N=100, n_phi=240)), ("rund", q, dict(kopplung_regel=False)),
                         ("lmc25", lmc25, {})):
        w = K.warp_ringe(qq, K.M_B_MILCHSTRASSE * K.M_sun, gemeinsam=True, eltern=[0, 1, 0], **kw)
        d = w[1] * K.kpc / 1e3
        aus[name] = (math.degrees(w[0]), d, abs(r - d) / e)
    return aus


@lru_cache(maxsize=None)
def _gas_frueh():
    """Das Gas als Gang der Dunklen Energie: alle Präzisionsdaten mit dem Gas im Kreis durch heute gelegt, Lage und Breite
    des Maximums frei, die Verhältnisse bei z = 0,34 (Bera) und z = 1,04 (Chowdhury) als Reste (daten.planck_anpassung,
    gas_frueh), F₀ beim Mittel aus _gas_raum gehalten und frei, dazu der Gang mit dem Maximum an der Sättigung (gas_bera) und
    ΛCDM. Gibt dict: je Fassung χ², F₀ mit Fehler, Verhältnis und Rest in Fehlern bei z = 0,34 und 1,04, z des Maximums."""
    import numpy as np
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    F0 = float(np.mean(_gas_raum()[1])) / 100
    z1, R1, e1, _ = D.hi_frueh_eins()
    _, Rb, eb, zb = D.gas_gezaehlt()
    def lesen(e):
        p = e["punkt"]
        r = lambda z: float(K.gas_kreis(1 / (1 + z), p["F0"], p["k"], p["a_s"])[0] / p["F0"])
        return dict(chi2=e["chi2"]["gesamt"], F0=p["F0"], fehler=e["fehler"].get("F0", 0.0), bera=r(zb),
                    bera_rest=(r(zb) - Rb) / eb, eins=r(z1), eins_rest=(r(z1) - R1) / e1, z_max=1 / p["a_s"] - 1)
    aus = {n: lesen(D.planck_anpassung(K.q_je_radiant(), alle, fest=f, zustand="gas_frueh",
                                       start={"F0": F0, "k": 3.0, "a_s": 0.45})) for n, f in (("fest", ("F0",)), ("frei", ()))}
    aus["saett"] = lesen(D.planck_anpassung(K.q_je_radiant(), alle, fest=("F0",), zustand="gas_bera", start={"F0": F0}))
    g = D.planck_anpassung(K.q_je_radiant(), alle, zustand="gezaehlt_frueh", start={"F0": F0, "k": 0.5, "a_s": 0.2})
    aus["form_frei"] = dict(lesen(g), form=g["punkt"]["form"])
    aus["spektren"] = (D.planck_anpassung(K.q_je_radiant(), alle, fest=("F0",), zustand="gas_frueh",
                                          start={"F0": F0, "k": 3.0, "a_s": 0.45})["chi2"]["spektren"]
                       - D.planck_anpassung(K.q_je_radiant(), alle, fest=("F0",), zustand="gas_bera",
                                            start={"F0": F0})["chi2"]["spektren"])
    aus["lcdm"], aus["F0_gezaehlt"] = D.planck_anpassung(None, alle)["chi2"]["gesamt"], F0
    return aus


@lru_cache(maxsize=None)
def _gas_abstaende():
    """Die Abstände allein (daten.abstaende_gas, Formteil beim Standardwert 1,06, F₀ beim Mittel aus _gas_raum): χ² ohne
    die Reste des Gases für den Formteil allein, das Gas mit dem Maximum an der Sättigung, das Gas durch die Messungen bei
    z = 0,34 und 1,04 (frueh), dasselbe aus der Mode gespeist (kopplung)."""
    import numpy as np
    F0 = float(np.mean(_gas_raum()[1])) / 100
    _, Rb, eb, _ = D.gas_gezaehlt()
    _, R1, e1, _ = D.hi_frueh_eins()
    def ohne(r, frueh):
        return r["chi2"] - (((r["R"] - Rb) / eb)**2 if r["R"] is not None else 0.0) - (
            ((r["R_eins"] - R1) / e1)**2 if frueh else 0.0)
    return dict(form=D.abstaende_gas(form=1.06, F0=0)["chi2"],
                saett=ohne(D.abstaende_gas(form=1.06, F0=F0), False),
                frueh=ohne(D.abstaende_gas(form=1.06, F0=F0, frueh=True), True),
                gekoppelt=ohne(D.abstaende_gas(form=1.06, F0=F0, frueh=True, kopplung=True), True))


@lru_cache(maxsize=None)
def _gas_raum():
    """Anteil des Gases unter der Gezeitendichte, wie das Ganze ihn zählt: je Stück in lg M_HI der Median der SPARC-Galaxien
    (daten.sparc_gas_gezeit mit der gelesenen Streuung, je Galaxie), gewichtet mit der Dichte des neutralen Gases je
    Stück aus ALFALFA (daten.alfalfa_hi_dichte mit massen, H₀ des Kerns), Stücke mit weniger als drei Galaxien ausgelassen.
    Gibt (Anteile für beide Wege der Streuung und drei Wege des Gases, dieselben mal dem gezählten Gas F in %)."""
    import numpy as np
    k = np.array([7.0, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])
    rho = D.alfalfa_hi_dichte(K.H0_kms, massen=k)[0]
    aus = []
    for rand in (False, True):
        f = D.hi_streuung_gang(rand=rand)[3]
        for weg in ("flaeche", "vgas", "hi"):
            lg, a = D.sparc_gas_gezeit(weg=weg, sigma=f, je=True)[4]
            stueck = [(lg >= lo) & (lg < hi) for lo, hi in zip(k[:-1], k[1:])]
            med = np.array([np.median(a[s]) if s.sum() >= 3 else np.nan for s in stueck])
            ok = np.isfinite(med)
            aus.append(float(np.sum(rho[ok] * med[ok]) / np.sum(rho[ok])))
    return aus, [x * D.gasanteil_heute() for x in aus]


@lru_cache(maxsize=None)
def _breiten_schnitt():
    """Überschuss der Bereiche, wenn nur Bereiche mit V₂₀₀ > √3·σ ihr Gas halten (kern.breitenfunktion mit v_min), mit der
    Streuung bis lg M_HI = 9,5 auf beiden Wegen (Linie der Scheibe, Rand von Hand), gegen beide Felder von ALFALFA und
    HIPASS: (Verhältnisse bei W50 = 32 bis 100 km/s, bei 100 bis 251 km/s)."""
    import numpy as np
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    ka = np.array([1.5, 2.0, 2.4])
    gemessen = [D.alfalfa_hi_dichte(punkt[2], feld=f, breiten=ka)[0] for f in ("fruehling", "herbst")]
    gemessen.append(D.hipass_hi_dichte(punkt[2], breiten=ka)[0])
    schmal, mittel = [], []
    for sg in (_gas_streuung()["leicht"], _gas_streuung()["leicht_rand"]):
        m = K.breitenfunktion(punkt, ka, sigma=sg, v_min=math.sqrt(3) * sg)
        schmal += [float(m[0] / g[0]) for g in gemessen]
        mittel += [float(m[1] / g[1]) for g in gemessen]
    return schmal, mittel


@lru_cache(maxsize=None)
def _breiten_anteil():
    """Anteil ε der Baryonen, den jeder Bereich bei W50 = 100 bis 251 km/s als Galaxie halten müsste, damit die Breitenfunktion
    (Linie mit der Streuung bis lg M_HI = 9,5, Schnitt V₂₀₀ > √3·σ) beide Felder von ALFALFA und HIPASS trifft, und χ² der
    blauen isolierten Linsen von KiDS-1000 mit dem Rest 1/ε − 1 als nie abgekühltem Gas mit Mode in der Form nach dem
    Einfall (daten.kids_heiss, lg M_b = 10, c = 4): (ε je Messung, χ² beim kleinsten und größten Rest, ohne, mit dem
    Einfachen)."""
    import numpy as np
    from scipy.optimize import brentq
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    ka, sg = np.array([2.0, 2.4]), _gas_streuung()["leicht"]
    gemessen = [D.alfalfa_hi_dichte(punkt[2], feld=f, breiten=ka)[0][0] for f in ("fruehling", "herbst")]
    gemessen.append(D.hipass_hi_dichte(punkt[2], breiten=ka)[0][0])
    eps = [10**brentq(lambda lg: math.log(K.breitenfunktion(punkt, ka, anteil=10**lg, sigma=sg,
                                                            v_min=math.sqrt(3) * sg)[0] / g), -3, 0, xtol=1e-4)
           for g in gemessen]
    rest = [1 / e - 1 for e in eps]
    chi = lambda f: D.kids_heiss("blau", einfall=(10.0, 4.0), fest=f)[1]
    return eps, chi(min(rest)), chi(max(rest)), chi(1e-9), chi(1.0), rest


@lru_cache(maxsize=None)
def _breiten_hipass():
    """Breitenfunktion der Bereiche mit allen Baryonen als Galaxie (Linie mit der gelesenen Streuung) gegen HIPASS
    (daten.hipass_hi_dichte) und gegen die
    beiden Felder von ALFALFA bei lg W50 = 1,8 bis 1,9: (HIPASS, Frühling, Herbst)."""
    import numpy as np
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    ka = np.array([1.8, 1.9])
    m = K.breitenfunktion(punkt, ka, sigma=_gas_streuung()["leicht"])[0]
    return (m / D.hipass_hi_dichte(punkt[2], breiten=ka)[0][0],) + tuple(
        m / D.alfalfa_hi_dichte(punkt[2], feld=f, breiten=ka)[0][0] for f in ("fruehling", "herbst"))


@lru_cache(maxsize=None)
def _polring():
    """Polarring um eine exponentielle Scheibe (6·10¹⁰ M☉, Dicke 0,1·R_d, viermal dichter, an und viermal dünner als die
    Schwelle, welt.polring): (kleinstes und größtes v_pol/v_eq mit der Regel bei 2 und 4·R_d, ebenso nach Newton, mit der
    Regel bei 8·R_d), dazu die größte Abweichung des Lösers von der Kuzmin-Scheibe analytisch (Feld eines Punktes im
    Spiegelbild) bei r = 1, 2, 4·a, nach Newton und mit der Regel."""
    import numpy as np
    import welt as W
    M = 6e10 * K.M_sun
    rs = math.sqrt(K.G * M / K.a0)
    regel, newton, aussen = [], [], []
    for g in (0.5, 1.0, 2.0):
        Rd = g * rs / math.sqrt(2)
        S = lambda r, Rd=Rd: M / (2 * math.pi * Rd**2) * np.exp(-r / Rd)
        for x in (2.0, 4.0, 8.0):
            n_, r_ = W.polring([(S, 0.1 * Rd)], None, x * Rd, schritt=(0.02 * Rd, 0.1 * Rd / 8))
            if x < 8:
                regel.append(r_); newton.append(n_)
            else:
                aussen.append(r_)
    nu = lambda y: 1 / -np.expm1(-np.sqrt(y))
    a, h = 2 * K.kpc, 0.005 * K.kpc
    S = lambda r: a * M / (2 * math.pi * (r**2 + a**2)**1.5)
    fehler = 0.0
    for x in (1.0, 2.0, 4.0):
        r = x * a
        th = (np.arange(4000) + 0.5) / 4000 * math.pi / 2
        d2 = r**2 + a**2 + 2 * a * r * np.cos(th)
        gr, ge = K.G * M * (r + a * np.cos(th)) / d2**1.5, K.G * M * r / (r**2 + a**2)**1.5
        an = (math.sqrt(gr.mean() / ge), math.sqrt((nu(K.G * M / d2 / K.a0) * gr).mean() / (nu(ge / K.a0) * ge)))
        lo = W.polring([(S, h)], None, r, schritt=(0.02 * a, h / 8), stufe=1.02)
        fehler = max(fehler, abs(an[0] - lo[0]), abs(an[1] - lo[1]))
    return min(regel), max(regel), min(newton), max(newton), min(aussen), max(aussen), fehler


@lru_cache(maxsize=None)
def _sigma_neigung():
    """Neigung von lg σ((1 + Q)·M_b) gegen lg M_b über M_b = 10⁸ bis 10^11,5 M☉ im linearen Spektrum heute (Punkt der Sättigung)."""
    import numpy as np
    p = _saett_punkt()
    k, P = K.spektrum_linear(p)
    rho = (p[0] + p[1]) * K.rho_krit100 / K.M_sun * K.Mpc**3
    lm = np.linspace(8.0, 11.5, 36)
    s = [math.sqrt(np.trapezoid(P * (lambda x: 3 * (np.sin(x) - x * np.cos(x)) / x**3)(k * (3 * (1 + K.q_je_radiant()) * 10**m / (4 * math.pi * rho))**(1 / 3))**2 * k**3, np.log(k)) / (2 * math.pi**2)) for m in lm]
    return float(np.polyfit(lm, np.log10(s), 1)[0])


def _saett_punkt():
    """Punkt der Sättigung aller Präzisionsdaten als Tupel (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀, w_a)."""
    p = _cmb()["saett"]["punkt"]
    return (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])


@lru_cache(maxsize=None)
def _lokaler_fluss():
    """Fluss um die Lokale Gruppe (kern.fluss_um_bereich, Baryonen von Milchstraße und M31, 1 + Q = 2π, Weltalter mit drei
    Vakuumteilen) gegen die Zwerge der Local Volume Database (daten.lokaler_fluss): dict mit R0 (mitatmend, mit der
    Regel) [Mpc], Zahl und Median des Rests mit 16 und 84 % über Ziehen."""
    import numpy as np
    M = list(_bereiche())
    t0 = K.weltalter(lambda a: 3 * K.dichte_vakuumteil())
    rr = np.linspace(0.55, 3.05, 26)
    vm = np.array([K.fluss_um_bereich(M, x * K.Mpc, t0, Q=K.q_je_radiant()) / 1e3 for x in rr])
    vr = np.array([K.fluss_um_bereich(M, x * K.Mpc, t0, "regel") / 1e3 for x in rr])
    d, v = D.lokaler_fluss()
    rest = v - np.interp(d, rr, vm)
    rng = np.random.default_rng(4)
    bs = [np.median(rest[rng.integers(0, len(rest), len(rest))]) for _ in range(4000)]
    return dict(R0=float(np.interp(0, vm, rr)), R0_regel=float(np.interp(0, vr, rr)), zahl=len(d), rest=float(np.median(rest)),
                unten=float(np.percentile(bs, 16)), oben=float(np.percentile(bs, 84)))


@lru_cache(maxsize=None)
def _wirte_satelliten(weg="daten"):
    """σ/v_flach der Satelliten um M31, M81, NGC 3627, M83, Cen A und NGC 4594 bis 400 kpc projiziert, σ um den
    Schwerpunkt der Baryonen von Wirt und Begleitern mit Geschwindigkeit, v_flach aus den Baryonen bis dorthin: weg
    "daten" mit daten.wirt_satelliten, "astropy" mit dem Winkel aus SkyCoord.separation und den Massen und
    Geschwindigkeiten unmittelbar aus der Tafel. Gibt {Wirt: (σ/v, Zahl, Mittel gegen den Wirt, Mittel gegen den
    Schwerpunkt, Fehler des Mittels [km/s])}, dazu unter "median" den Median von σ/v und unter "chi2" Σ 2N·ln²(√3·σ/v)."""
    import csv
    import numpy as np
    wirte = ("m_031", "m_081", "ngc_3627", "m_083", "ngc_5128", "ngc_4594")
    aus = {}
    if weg == "daten":
        for w in wirte:
            r = D.wirt_satelliten(w, 400.0)
            aus[w] = (r["sigma"] / (K.v_flach(r["M_b"] * K.M_sun) / 1e3), r["zahl"], r["mittel"], r["mittel_S"], r["fehler"])
    else:
        import astropy.coordinates as ko
        import astropy.units as u
        z = list(csv.DictReader(open(os.path.join(D.ORT, "lvdb", "comb_all.csv"))))
        f = lambda x: float(x) if x not in ("", "nan") else None
        m = lambda r: 10**f(r["mass_stellar"]) + 1.33 * (10**f(r["mass_HI"]) if f(r["mass_HI"]) else 0) if f(r["mass_stellar"]) else 0.0
        for w in wirte:
            W = next(r for r in z if r["key"] == w)
            cw = ko.SkyCoord(f(W["ra"]) * u.deg, f(W["dec"]) * u.deg)
            mit = [r for r in z if r["host"] == w and (r["table"].startswith("dwarf") or r["table"] == "misc")
                   and cw.separation(ko.SkyCoord(f(r["ra"]) * u.deg, f(r["dec"]) * u.deg)).rad * f(W["distance"]) <= 400.0]
            mw = K.M_B_M31 if w == "m_031" else m(W)
            M = mw + sum(m(r) for r in mit)
            if w == "m_031":
                tafel = {r["key"]: r for r in z}
                s = [(x["v"], m(tafel[x["name"]]), x["tafel"]) for x in D.satelliten_schar("m_031", begleiter=True)
                     if x["R"] <= 400.0]
                v0 = 0.0
            else:
                s = [(f(r["vlos_systemic"]), m(r), r["table"]) for r in mit if f(r["vlos_systemic"]) is not None]
                v0 = f(W["vlos_systemic"])
            v_S = np.average([v0] + [v for v, _, _ in s], weights=[mw] + [g for _, g, _ in s])
            dv = np.array([v for v, _, t in s if t.startswith("dwarf")])
            aus[w] = (math.sqrt(np.mean(np.square(dv - v_S))) / ((K.G * M * K.M_sun * K.a0)**0.25 / 1e3), len(dv),
                      dv.mean() - v0, dv.mean() - v_S, dv.std(ddof=1) / math.sqrt(len(dv)))
    aus["median"] = float(np.median([aus[w][0] for w in wirte]))
    aus["chi2"] = float(sum(2 * aus[w][1] * math.log(math.sqrt(3) * aus[w][0])**2 for w in wirte))
    return aus


@lru_cache(maxsize=None)
def _bereiche():
    """Baryonen der Bereiche von Milchstraße und M31 samt ihren Satelliten (daten.bereich_baryonen) [kg]."""
    return D.bereich_baryonen("mw")[0] * K.M_sun, D.bereich_baryonen("m_031")[0] * K.M_sun


@lru_cache(maxsize=None)
def _bereiche_paar():
    """Zweiter Weg zu daten.m31_paar: die Schwerpunkte der Bereiche von Milchstraße und M31 im Rahmen ICRS statt
    galaktozentrisch, das Zentrum der Milchstraße aus dem Ursprung des galaktozentrischen Rahmens zurückgeführt, die
    Glieder mit Eigenbewegung aus der Tafel. Gibt (Abstand [kpc], Annäherung [km/s], quer [km/s])."""
    import csv
    import numpy as np
    import astropy.coordinates as ko
    import astropy.units as u
    s = K.SONNE_GALAKTISCH
    rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                               galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
    tafel = {z["key"]: z for z in csv.DictReader(open(os.path.join(D.ORT, "lvdb", "comb_all.csv")))}
    f = lambda x: float(x) if x not in ("", "nan") else math.nan
    kv = lambda c: (c.cartesian.xyz.to(u.kpc).value, c.velocity.d_xyz.to(u.km / u.s).value)
    zentrum = ko.SkyCoord(x=0 * u.kpc, y=0 * u.kpc, z=0 * u.kpc, v_x=0 * u.km / u.s, v_y=0 * u.km / u.s,
                          v_z=0 * u.km / u.s, frame=rahmen).transform_to(ko.ICRS())
    icrs = lambda z: ko.SkyCoord(ra=f(z["ra"]) * u.deg, dec=f(z["dec"]) * u.deg, distance=f(z["distance"]) * u.kpc,
                                 pm_ra_cosdec=f(z["pmra"]) * u.mas / u.yr, pm_dec=f(z["pmdec"]) * u.mas / u.yr,
                                 radial_velocity=f(z["vlos_systemic"]) * u.km / u.s, frame="icrs")
    masse = lambda z: (10**f(z["mass_stellar"]) if math.isfinite(f(z["mass_stellar"])) else 0.0) + 1.33 * (
        10**f(z["mass_HI"]) if math.isfinite(f(z["mass_HI"])) else 0.0)
    punkte = []
    for wirt, m0, c0 in (("mw", K.M_B_MILCHSTRASSE, zentrum), ("m_031", K.M_B_M31, icrs(tafel["m_031"]))):
        teile = [(m0,) + kv(c0)] + [(masse(tafel[k]),) + kv(icrs(tafel[k])) for k in D.bereich_baryonen(wirt, glieder=True)
                                    if all(math.isfinite(f(tafel[k][n])) for n in ("pmra", "pmdec", "vlos_systemic", "distance"))]
        m = sum(t[0] for t in teile)
        punkte.append((sum(t[0] * t[1] for t in teile) / m, sum(t[0] * t[2] for t in teile) / m))
    x, v = punkte[1][0] - punkte[0][0], punkte[1][1] - punkte[0][1]
    r = float(np.linalg.norm(x))
    v_r = float(x @ v) / r
    return r, -v_r, math.sqrt(float(v @ v) - v_r**2)


@lru_cache(maxsize=None)
def _licht_einfall():
    """Licht um isolierte Galaxien bei 1,9 Mpc (KiDS-1000, vier Massenbins, Punktmasse M_b): Mittel von lg(g_obs/Regel) und
    lg(g_obs/((1 + Q)·g_N)), dazu ((1 + Q)·g_N + Umgebung)/g_obs mit der linearen Faltung um den Bereich (1 + Q)·M_b
    (kern.linse_umgebung, z = 0,2), und die größte Abweichung der beiden Wege der Umgebung."""
    import numpy as np
    Q, p = K.q_je_radiant(), _saett_punkt()
    a, b_, c, w = [], [], [], 0.0
    for b in D.kids_isoliert():
        i = int(np.argmin(abs(b["R"] - 1915.0)))
        gN = K.G * b["M_b"] * K.M_sun / (b["R"][i] * K.kpc)**2
        u = K.linse_umgebung(p, (1 + Q) * b["M_b"], [b["R"][i] / 1e3])[0]
        w = max(w, abs(u / K.linse_umgebung(p, (1 + Q) * b["M_b"], [b["R"][i] / 1e3], weg="xi")[0] - 1))
        a.append(math.log10(b["g"][i] / K.regel_monopol(gN, 0.0)))
        b_.append(math.log10(b["g"][i] / ((1 + Q) * gN)))
        c.append(((1 + Q) * gN + u) / b["g"][i])
    return float(np.mean(a)), float(np.mean(b_)), float(np.mean(c)), w


@lru_cache(maxsize=None)
def _scheibe_lage():
    """Halbmassenradius kalter Scheiben mit der Regel als Erzeuger (Feldform auf dem Gitter, welt.scheibe_lauf, Zufall 10)
    am Anfang und nach 3 Gyr [r*], je viermal dichter, an der Schwelle und viermal dünner (groesse 0,5, 1, 2)."""
    import welt as W
    aus = []
    for g in (0.5, 1.0, 2.0):
        e = W.scheibe_lauf(feld="wirbelfrei", saat=10, groesse=g)
        aus += [e["rh"][0], e["rh"][-1]]
    return tuple(aus)


@lru_cache(maxsize=None)
def _drehung_grenze():
    """Größte Spinzahl eines gehaltenen Körpers, der sich ganz durch Drehung trägt, j = k/√(3Δs/R): (k der Form n = 3,
    j der Sonne damit, Weißer Zwerg an der Grenzmasse M [M☉], R [km], k, j, Neutronenstern an der Grenze M [M☉], R [km], k, j)."""
    k3 = K.polytrope(3.0, traegheit=True)[3]
    aus = [k3, k3 / math.sqrt(K.G * K.M_sun / (K.R_SONNE * K.c**2))]
    for traeger, masse in ((K.m_e, 2 * K.m_u), (K.m_n, K.m_n)):
        M, R, d, _, k = K.tov_sperre(traeger, m_masse=masse, traegheit=True)
        aus += [M / K.M_sun, R / 1e3, k, k / math.sqrt(d)]
    return tuple(aus)


@lru_cache(maxsize=None)
def _haufen_drehung():
    """Verwindung der Haufen aus denselben Läufen vieler Körper (welt.haufen_drehung): dict je Größe (kleinster, größter
    Wert über die vier Läufe): lam, j, M200 [10¹⁵ M☉], phi (Winkel der Horizonte zur Tiefe, cos φ = 1/j, Grad), dazu linear
    lam_lin, cos und verhaeltnis des Kommutators, und das Tempo der Drehung j/R₂₀₀ = √2·λ'·V₂₀₀ [km/s] des Gases (linear)
    und der Mode (am Ende), v_gas und v_mode."""
    import welt as W
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"])
    d = [W.haufen_drehung(punkt, R, saat=s) for R in (25.0, 18.0) for s in (1, 2)]
    li = [W.haufen_drehung(punkt, R, saat=s, linear=True) for R in (25.0, 18.0) for s in (1, 2)]
    mm = lambda xs: (min(xs), max(xs))
    return dict(lam=mm([x["lam"] for x in d]), j=mm([x["j"] for x in d]), M200=mm([x["M200"] / 1e15 for x in d]),
                phi=mm([math.degrees(math.acos(1 / x["j"])) for x in d]), lam_lin=mm([x["lam"] for x in li]),
                cos=mm([x["cos"] for x in li]), verhaeltnis=mm([x["verhaeltnis"] for x in li]),
                v_gas=mm([math.sqrt(2) * y["lam"] * x["V200"] for x, y in zip(d, li)]),
                v_mode=mm([math.sqrt(2) * x["lam"] * x["V200"] for x in d]))


@lru_cache(maxsize=None)
def _kreise_lesung():
    """Kreispaare der Zelle an Planck SMICA: beim Radius α des besten Punkts (gerundet auf 0,1°), Glättung 0,25°, Gitter der
    Richtungen Nside 128: das beste Paar an Planck und an drei Nullproben (zweiter Kreis um eine zufällige Richtung
    gedreht), die Streuung je Paar an Planck, die beste Dreiheit senkrechter Paare an Planck und an vier Nullproben, die
    Erwartung der Zelle aus der Zerlegung in Monopol und Doppler (kern.kreis_erwartung). Zweiter Weg zur Erwartung: der
    gerechnete Himmel der Würfelzelle (kern.himmel_zelle, ℓ ≤ 100) mit Kreisen von 20°, Glättung 2°, vier Ziehungen, an
    den Achsen der Zelle gemessen, gegen die Zerlegung dort. Gibt dict."""
    import camb
    import healpy as hp
    import numpy as np
    p = _cmb()["saett"]["punkt"]
    punkt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["tau"], p["w0"], p["wa"])
    al = round(_zelle_kreis()[2], 1)
    S, _, an = D.kreise_gegenueber([al], fwhm_grad=0.25, richtungen=128, nside=512, schritt_grad=0.1)
    ok = an[:, 0] > 0.5
    paar = float(S[ok, 0].max())
    streu = float(S[ok, 0].std())
    null = []
    for s in range(3):
        Sn, _, a2 = D.kreise_gegenueber([al], fwhm_grad=0.25, richtungen=128, nside=512, schritt_grad=0.1, zufall=s)
        null.append(float(Sn[a2[:, 0] > 0.5, 0].max()))
    drei = D.kreise_dreier([al], fwhm_grad=0.25, richtungen=128, nside=512)[0][0]
    drei_null = [D.kreise_dreier([al], fwhm_grad=0.25, richtungen=128, nside=512, zufall=s)[0][0] for s in range(4)]
    erw = K.kreis_erwartung(punkt, al, 0.25)
    extra = dict(w=p["w0"], wa=p["wa"], dark_energy_model="ppf")
    par = camb.set_params(ombh2=p["wb"], omch2=p["wc"], H0=p["H0"], As=math.exp(p["lnAs"]) * 1e-10, ns=p["ns"],
                          tau=p["tau"], mnu=0.06, **extra)
    bg = camb.get_background(par)
    chi = bg.conformal_time(0) - bg.tau_maxvis
    alms, _ = K.himmel_zelle(punkt, 2 * chi * math.cos(math.radians(20.0)), lmax=100, saat=[1, 2, 3, 4], k_faktor=2.5)
    gem = []
    for a in alms:
        T = hp.alm2map(hp.smoothalm(np.array(a), fwhm=math.radians(2.0)), 128)
        gem += list(D.kreise_gegenueber([20.0], karte=(T, np.ones_like(T, bool), 128), achsen=np.eye(3))[0][:, 0])
    gem = np.array(gem)
    zwanzig = K.kreis_erwartung(punkt, 20.0, 2.0)
    return dict(alpha=al, paar=paar, streu=streu, null=(min(null), max(null)), drei=drei,
                drei_null=(min(drei_null), max(drei_null)), erwartet=erw[:2], doppler=erw[2:],
                himmel=(float(gem.mean()), float(gem.std() / math.sqrt(gem.size))), zerlegung20=zwanzig[:2])


@lru_cache(maxsize=None)
def _zelle_kreis():
    """Kreispaare der Würfelzelle im besten Punkt mit Sättigung: L/(2χ*) mit Fehler aus der Kovarianz (ω_b, a₀, w_a) und
    der Radius α, auch bei ±1 Fehler."""
    import numpy as np
    s = _cmb()["saett"]
    p, n, C = s["punkt"], s["kovarianz"]["namen"], np.array(s["kovarianz"]["werte"])
    f = lambda pp: K.zelle_kreis(pp["wb"], pp["a0"] * 1e-10, pp["wa"])[0]
    x = f(p)
    g = np.array([(f(dict(p, **{k: p[k] + 1e-4 * max(1, abs(p[k]))})) - x) / (1e-4 * max(1, abs(p[k])))
                  if k in ("wb", "a0", "wa") else 0.0 for k in n])
    e = float(np.sqrt(g @ C @ g))
    grad = lambda y: math.degrees(math.acos(min(y, 1.0)))
    return x, e, grad(x), grad(x + e), grad(x - e)


@lru_cache(maxsize=None)
def _kids_teile(R_lin, R_nl):
    """KiDS-1000-Scherung im besten Punkt aller Präzisionsdaten mit Sättigung, linearer Teil und nichtlinearer Rest mit
    den Amplituden R_lin, R_nl (None frei): wie D.kids_scherung_chi2 mit teile."""
    p = D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung',
                           start={'wa': -0.56})['punkt']
    return D.kids_scherung_chi2((p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']), teile=(R_lin, R_nl))


@lru_cache(maxsize=None)
def _uebergabe():
    """Übergebene Mode gegen die Mode als Staub im besten Punkt mit Sättigung (Q = 2π − 1, a₀ = 1,06): ν/(1 + Q) in R₂₀₀
    bei 10¹¹, 10¹² und 10¹³ M☉ und die Masse M₂₀₀, bei der R₂₀₀ an der Naht liegt [M☉]."""
    from scipy.optimize import brentq
    Q, H = K.q_je_radiant(), _cmb()["saett"]["punkt"]["H0"] * 1e3 / K.Mpc
    f = lambda m: K.uebergabe_linse(m * K.M_sun, H, Q)[2]
    return f(1e11), f(1e12), f(1e13), 10**brentq(lambda x: f(10**x) - 1, 12, 16)


@lru_cache(maxsize=None)
def _ring():
    """Größter Abstand D = D_l·D_ls/D_s vor der letzten Streufläche (Formteil, a₀ = 1,06) über die Lage der Linse, aus dem
    Kern, und derselbe Abstand auf zweitem Weg (Trapez über z, H aus den Dichten): (z_l, D, G·M/(R_E²·a₀), D zweiter Weg)."""
    import numpy as np
    from scipy.optimize import minimize_scalar
    zl = minimize_scalar(lambda z: -K.einsteinring(z)[0], bounds=(0.3, 3), method="bounded").x
    Dk, g = K.einsteinring(zl)
    rr, _, rm = K.dichte_heute()
    zz = np.linspace(0, K.z_stern, 2000001)
    a = 1 / (1 + zz)
    dm = np.concatenate(([0], np.cumsum(np.diff(zz) * 0.5 * (lambda h: h[1:] + h[:-1])(K.c / np.sqrt(8 * np.pi * K.G / 3 * (rr / a**4 + rm / a**3 + 3 * K.dichte_vakuumteil()))))))
    Dl = np.interp(zl, zz, dm)
    return zl, Dk, g / K.a0, Dl / (1 + zl) * (1 - Dl / dm[-1])


@lru_cache(maxsize=None)
def _form_atom(form):
    """Zustandsteil heute und gebundener Gasanteil F, wenn der Formteil a₀ = form ist (Sättigung bei ihm, alle
    Präzisionsdaten): (χ², a₀ des Ganzen, Zustandsteil, Fehler, F, Fehler)."""
    import numpy as np
    Q = K.q_je_radiant()
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    sf = D.planck_anpassung(Q, alle, zustand="saettigung", start={"wa": -0.56}, form=form)
    je = K.zustandsteil_je_gas(Q) / 3
    f = lambda pp: np.array([(pp["a0"] / form)**2 - 1, ((pp["a0"] / form)**2 - 1) / je])
    p, n, C = sf["punkt"], sf["kovarianz"]["namen"], np.array(sf["kovarianz"]["werte"])
    w = f(p)
    G = np.array([(f(dict(p, **{k: p[k] + 1e-6 * max(1, abs(p[k]))})) - w) / (1e-6 * max(1, abs(p[k]))) for k in n]).T
    e = np.sqrt(np.einsum("ij,jk,ik->i", G, C, G))
    return sf["chi2"]["gesamt"], p["a0"], w[0], e[0], w[1], e[1]


@lru_cache(maxsize=None)
def _falte_einfall():
    """Einfall radialer Schalen auf die mittlere Spitze des Spektrums im besten Punkt mit Sättigung (R_L = 25 Mpc, 400
    Schalen): Phase ω·(t_Kreuzung − t_Umkehr) je Schale, die nach ihrer Umkehr gekreuzt wird. Gibt (Median, unteres und
    oberes Viertel, Zahl)."""
    import numpy as np
    p = _cmb()["saett"]["punkt"]
    e = K.einfall_haufen((p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"]), 25.0, n=400)
    ok = ~np.isnan(e["t_kr"]) & ~np.isnan(e["t_um"]) & (e["t_kr"] > e["t_um"])
    ph = np.sqrt(K.G * e["M_schale"][ok] / e["r_um"][ok]**3) * (e["t_kr"][ok] - e["t_um"][ok])
    return float(np.median(ph)), float(np.percentile(ph, 25)), float(np.percentile(ph, 75)), int(ok.sum())


@lru_cache(maxsize=None)
def _puls_hubble():
    """Der Puls des Zustandsteils in H: In der eigenen Form des Modells an allen Präzisionsdaten (zustand "gas"),
    ρ_DE = ρ_s·(3 + 16,8·F(a)) gegen den Formteil 3ρ_s bei gleichen Dichten der Faltung und Strahlung. Gibt (größter
    Überschuss von H in %, z dort, Überschuss heute in %, fällt H stets mit der Zeit (1/0), zweiter Weg linear
    ½·Δρ/ρ in %)."""
    import numpy as np
    g = D.planck_anpassung(K.q_je_radiant(), ("spektren", "linse", "bao", "sn", "wachstum", "bbn"), zustand="gas")["punkt"]
    Q, rs = g["Q"], K.dichte_vakuumteil(g["form"] * 1e-10)
    rm, rr = ((1 + Q) * g["wb"] + 0.06 / 93.14) * K.rho_krit100, K.dichte_heute()[0]
    a = np.geomspace(1e-3, 1, 200001)
    zus = K.zustandsteil_je_gas(Q) * rs * K.zustand_gas(a, g["Fs"], g["k"], (rm / (3 * rs))**(1 / 3))
    form = rr / a**4 + rm / a**3 + 3 * rs
    x = np.sqrt((form + zus) / form) - 1
    i = int(np.argmax(x))
    return (100 * float(x[i]), float(1 / a[i] - 1), 100 * float(x[-1]),
            float(np.all(np.diff(np.sqrt(form + zus)) < 0)), 100 * float(np.max(0.5 * zus / form)))


@lru_cache(maxsize=None)
def _gas_form():
    """Dunkle Energie in der Form des Modells (Zustandsteil im kalten Gas, F nie negativ, Maximum an der Sättigung) an
    allen Präzisionsdaten: χ² gegen ΛCDM und den Gang mit Sättigung (Formteil 1,06), F heute und am Maximum mit Fehler aus
    der Kovarianz, z des Maximums, F des Gangs w₀, w_a bei z = 2."""
    import numpy as np
    Q = K.q_je_radiant()
    alle = ("spektren", "linse", "bao", "sn", "wachstum", "bbn")
    g = D.planck_anpassung(Q, alle, zustand="gas")
    je = K.zustandsteil_je_gas(Q) / 3
    def werte(pp):
        a_s = (((1 + Q) * pp["wb"] + 0.06 / 93.14) * K.rho_krit100 / K.dichte_dunkle_energie(pp["form"] * 1e-10))**(1 / 3)
        return np.array([float(K.zustand_gas(1.0, pp["Fs"], pp["k"], a_s)), pp["Fs"], 1 / a_s - 1])
    p, n, C = g["punkt"], g["kovarianz"]["namen"], np.array(g["kovarianz"]["werte"])
    w = werte(p)
    st = {k: 1e-6 * max(1, abs(p[k])) for k in n}
    G = np.array([(werte(dict(p, **{k: p[k] + st[k]})) - w) / st[k] for k in n]).T
    e = np.sqrt(np.einsum("ij,jk,ik->i", G, C, G))
    sf = _cmb()["saett_form"]["punkt"]
    x = (sf["a0"] / 1.06)**2 - 1
    r2 = (1 / 3)**(-3 * (1 + sf["w0"] + sf["wa"])) * math.exp(-3 * sf["wa"] * (1 - 1 / 3))   # ρ_DE(z = 2)/ρ_DE(0)
    gf = D.planck_anpassung(Q, alle, zustand="gas_frei")
    a_s = (((1 + Q) * p["wb"] + 0.06 / 93.14) * K.rho_krit100 / K.dichte_dunkle_energie(p["form"] * 1e-10))**(1 / 3)
    z1 = D.hi_frueh_eins()[0]                                # F bei z ≈ 1 gegen heute in der eigenen Form
    eins = float(K.zustand_gas(1 / (1 + z1), p["Fs"], p["k"], a_s) / K.zustand_gas(1.0, p["Fs"], p["k"], a_s))
    return dict(frei=(gf["chi2"]["gesamt"], 1 / gf["punkt"]["a_s"] - 1, gf["punkt"]["Fs"], gf["fehler"]["Fs"]),
                chi2=g["chi2"]["gesamt"], lcdm=_cmb()["alle"][0]["chi2"]["gesamt"],
                saett=_cmb()["saett_form"]["chi2"]["gesamt"], F0=(w[0], e[0]), Fs=(w[1], e[1]), z_s=w[2],
                F_z2=((1 + x) * r2 - 1) / je, F_eins=eins)


@lru_cache(maxsize=None)
def _haufen_staub():
    """Konzentration der Mode in der Form nach dem Einfall (NFW), die X-COP und HIFLUGCS lesen."""
    return D.haufen_c("xcop"), D.haufen_c("hiflugcs")


@lru_cache(maxsize=None)
def _sparc_ganzes():
    """SPARC mit den Abständen des Hubble-Flusses aus dem H₀ des Ganzen (Q = 2π − 1, mit den Spitzen): a₀, 16–84 %,
    Υ_d nach Stern und Gas, H₀, Zahl der Galaxien mit Hubble-Fluss-Abstand und aller."""
    H0 = _ganzes_a0(True)[3]
    yd = D.sparc_eichung(K.kpc, H0)
    a, (lo, hi) = D.sparc_schwelle(K.kpc, yd=yd, H0=H0)
    aus = D.sparc_auswahl(yd, H0)
    return a / 1e-10, lo / 1e-10, hi / 1e-10, yd, H0, sum(p["fD"] == 1 for _, p, *_ in aus), len(aus)


@lru_cache(maxsize=None)
def _ganzes_bestes(Q, spitzen):
    """Kleinstes χ² des Ganzen bei festem Q über a₀ und ω_b: (χ², a₀, ω_b)."""
    from scipy.optimize import minimize
    r = minimize(lambda p: _ganzes_chi(p[0], p[1], Q, spitzen), [1.085, K.omega_b], method="Nelder-Mead",
                 options=dict(xatol=1e-7, fatol=1e-6))
    return float(r.fun), float(r.x[0]), float(r.x[1])


@lru_cache(maxsize=None)
def _ganzes_a0(spitzen):
    """a₀ des Ganzen bei Q = 2π − 1 (ω_b mitgeschätzt): Bestwert, Δχ² = 1-Grenzen, H₀ und Ω_m im Bestwert."""
    from scipy.optimize import brentq, minimize_scalar
    Q = 2 * math.pi - 1
    m, a, ob = _ganzes_bestes(Q, spitzen)
    g = lambda x: minimize_scalar(lambda o: _ganzes_chi(x, o, Q, spitzen), bounds=(0.0218, 0.0230), method="bounded",
                                  options=dict(xatol=1e-8)).fun - m - 1
    alt = K.omega_b
    K.omega_b = ob
    try:
        H0, om = _h0_q(Q, a * 1e-10), _om_or_q(Q)(a * 1e-10)[0]
    finally:
        K.omega_b = alt
    return a, brentq(g, a - 0.05, a), brentq(g, a, a + 0.05), H0, om


def _sigma(dchi2, frei):
    from scipy.stats import chi2, norm
    return norm.isf(chi2.sf(dchi2, frei) / 2)


@lru_cache(maxsize=None)
def _modellzahlen():
    """Die Prüfliste: (Baustein, Stelle im Text, gerechneter Wert, Wert im Text). Erst beim Prüflauf gerechnet."""
    return [
        ("A:wardenclyffe_de", "Wheelers Näherung 1/(1 + 0,9r/h) weicht davon um 0,2 % ab", 100 * (T.wendel()["L"] / T.wendel()["L_wheeler"] - 1), 0.2),
        ("A:wardenclyffe_de", "Kupferrohr (10,7 km, 5 cm", T.wendel(aus=T.formen()[0])["draht"] / 1e3, 10.7),
        ("A:wardenclyffe_de", "bleibt fest bei 342 m²", T.TURM["hoehe"] * T.AUSLEGUNG["wendel"], 342),
        ("A:wardenclyffe_de", "Beide Wege stimmen auf 12 Stellen", round(1e12 * abs(T.kopplung_fluss(T.reichweite(aus=T.formen()[0]), T.formen()[0]) / T.kopplung(T.reichweite(aus=T.formen()[0]), T.formen()[0])[0] - 1)), 0),
        ("A:wardenclyffe_de", "in 10 km dasselbe auf 0,2 %", 100 * abs(T.kopplung(1e4, T.formen()[0])[1] / T.kopplung(1e4, T.formen()[0])[0] - 1), 0.2),
        ("A:wardenclyffe_de", "Die geschlossene Form stimmt mit der vollen Rechnung auf 0,7 %", 100 * max(abs(T.guete_geschlossen(T.bauform(x)) / (2 * math.pi * T.wendel(aus=T.bauform(x))["f"] * T.wendel(aus=T.bauform(x))["L"] / T.verluste(aus=T.bauform(x))["R_D"]) - 1) for x in (0.3, 0.5, 1.0, 2.0, 3.0)), 0.7),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.bauform(0.3)["wendel"], 10.1),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.bauform(0.3)["hoehe"], 33.8),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.wendel(aus=T.bauform(0.3))["N"], 337),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.wendel(aus=T.bauform(0.3))["f"] / 1e3, 8.46),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.verluste(aus=T.bauform(0.3))["Q"], 9396),
        ("A:wardenclyffe_de", "| 0,3 | 10,1 m | 33,8 m | 337 | 8,46 kHz | 9396 | 91 m |", T.reichweite(aus=T.bauform(0.3)), 91),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.bauform(0.5)["wendel"], 13.1),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.bauform(0.5)["hoehe"], 26.2),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.wendel(aus=T.bauform(0.5))["N"], 261),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.wendel(aus=T.bauform(0.5))["f"] / 1e3, 7.15),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.verluste(aus=T.bauform(0.5))["Q"], 10280),
        ("A:wardenclyffe_de", "| 0,5 | 13,1 m | 26,2 m | 261 | 7,15 kHz | 10280 | 106 m |", T.reichweite(aus=T.bauform(0.5)), 106),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.bauform(1.0)["wendel"], 18.5),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.bauform(1.0)["hoehe"], 18.5),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.wendel(aus=T.bauform(1.0))["N"], 184),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.wendel(aus=T.bauform(1.0))["f"] / 1e3, 5.62),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.verluste(aus=T.bauform(1.0))["Q"], 10737),
        ("A:wardenclyffe_de", "| 1,0 | 18,5 m | 18,5 m | 184 | 5,62 kHz | 10737 | 129 m |", T.reichweite(aus=T.bauform(1.0)), 129),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.bauform(2.0)["wendel"], 26.2),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.bauform(2.0)["hoehe"], 13.1),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.wendel(aus=T.bauform(2.0))["N"], 130),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.wendel(aus=T.bauform(2.0))["f"] / 1e3, 4.43),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.verluste(aus=T.bauform(2.0))["Q"], 10208),
        ("A:wardenclyffe_de", "| 2,0 | 26,2 m | 13,1 m | 130 | 4,43 kHz | 10208 | 157 m |", T.reichweite(aus=T.bauform(2.0)), 157),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.bauform(3.0)["wendel"], 32.0),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.bauform(3.0)["hoehe"], 10.7),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.wendel(aus=T.bauform(3.0))["N"], 106),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.wendel(aus=T.bauform(3.0))["f"] / 1e3, 3.89),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.verluste(aus=T.bauform(3.0))["Q"], 9508),
        ("A:wardenclyffe_de", "| 3,0 | 32,0 m | 10,7 m | 106 | 3,89 kHz | 9508 | 176 m |", T.reichweite(aus=T.bauform(3.0)), 176),
        ("A:wardenclyffe_de", "Das Maximum liegt bei D/h = 1,01", T.form_guete()[0], 1.01),
        ("A:wardenclyffe_de", "| Durchmesser der Wendel | 6 m | 18,6 m | 13,6 m |", T.formen()[0]["wendel"], 18.6),
        ("A:wardenclyffe_de", "| Durchmesser der Wendel | 6 m | 18,6 m | 13,6 m |", T.formen()[1]["wendel"], 13.6),
        ("A:wardenclyffe_de", "| Höhe der Wendel | 57 m | 18,4 m | 25,1 m |", T.formen()[0]["hoehe"], 18.4),
        ("A:wardenclyffe_de", "| Höhe der Wendel | 57 m | 18,4 m | 25,1 m |", T.formen()[1]["hoehe"], 25.1),
        ("A:wardenclyffe_de", "| Windungen | 570 | 184 | 251 |", T.wendel(aus=T.formen()[0])["N"], 184),
        ("A:wardenclyffe_de", "| Windungen | 570 | 184 | 251 |", T.wendel(aus=T.formen()[1])["N"], 251),
        ("A:wardenclyffe_de", "| Kuppel | 20,7 m | 18,6 m | 13,6 m |", T.formen()[0]["kuppel"], 18.6),
        ("A:wardenclyffe_de", "| Kuppel | 20,7 m | 18,6 m | 13,6 m |", T.formen()[1]["kuppel"], 13.6),
        ("A:wardenclyffe_de", "| Erdrohre | 16 | 16 | 64 |", T.SUPRA["staebe"], 64),
        ("A:wardenclyffe_de", "| Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz |", T.wendel()["f"] / 1e3, 8.31),
        ("A:wardenclyffe_de", "| Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz |", T.wendel(aus=T.formen()[0])["f"] / 1e3, 5.57),
        ("A:wardenclyffe_de", "| Frequenz | 8,31 kHz | 5,57 kHz | 6,94 kHz |", T.wendel(aus=T.formen()[1])["f"] / 1e3, 6.94),
        ("A:wardenclyffe_de", "| Güte Q_s | 5969 | 10765 | 7,33·10⁵ |", T.verluste()["Q"], 5969),
        ("A:wardenclyffe_de", "| Güte Q_s | 5969 | 10765 | 7,33·10⁵ |", T.verluste(aus=T.formen()[0])["Q"], 10765),
        ("A:wardenclyffe_de", "| Güte Q_s | 5969 | 10765 | 7,33·10⁵ |", T.verluste(aus=T.formen()[1])["Q"] / 1e5, 7.33),
        ("A:wardenclyffe_de", "| Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ |", T.betrieb()["nu"], 1900),
        ("A:wardenclyffe_de", "| Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ |", T.betrieb(aus=T.formen()[0])["nu"], 3427),
        ("A:wardenclyffe_de", "| Rückkehr ν = Q_s/π | 1900 | 3427 | 2,33·10⁵ |", T.betrieb(aus=T.formen()[1])["nu"] / 1e5, 2.33),
        ("A:wardenclyffe_de", "| Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag |", T.betrieb()["W"] / 1e3, 22.9),
        ("A:wardenclyffe_de", "| Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag |", T.betrieb(aus=T.formen()[0])["W"] / 1e3, 61.5),
        ("A:wardenclyffe_de", "| Energie bei 200 kW | 22,9 kJ | 61,5 kJ | über dem Durchschlag |", float(T.betrieb(aus=T.formen()[1])["E_anteil"] > 1), 1),
        ("A:wardenclyffe_de", "| Kuppel bei 200 kW | 4,91 MV | 8,06 MV | über dem Durchschlag |", T.betrieb()["V"] / 1e6, 4.91),
        ("A:wardenclyffe_de", "| Kuppel bei 200 kW | 4,91 MV | 8,06 MV | über dem Durchschlag |", T.betrieb(aus=T.formen()[0])["V"] / 1e6, 8.06),
        ("A:wardenclyffe_de", "| Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW |", T.grenzleistung()[0] / 1e6, 7.98),
        ("A:wardenclyffe_de", "| Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW |", T.grenzleistung(T.formen()[0])[0] / 1e6, 2.39),
        ("A:wardenclyffe_de", "| Durchschlag bei | 7,98 MW | 2,39 MW | 17,6 kW |", T.grenzleistung(T.formen()[1])[0] / 1e3, 17.6),
        ("A:wardenclyffe_de", "| Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ |", T.grenzleistung()[1] / 1e6, 0.91),
        ("A:wardenclyffe_de", "| Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ |", T.grenzleistung(T.formen()[0])[1] / 1e6, 0.74),
        ("A:wardenclyffe_de", "| Energie am Durchschlag | 0,91 MJ | 0,74 MJ | 0,30 MJ |", T.grenzleistung(T.formen()[1])[1] / 1e6, 0.30),
        ("A:wardenclyffe_de", "| 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m |", T.reichweite(), 54),
        ("A:wardenclyffe_de", "| 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m |", T.reichweite(aus=T.formen()[0]), 130),
        ("A:wardenclyffe_de", "| 50 % an einen gleichen Bau bis | 54 m | 130 m | 454 m |", T.reichweite(aus=T.formen()[1]), 454),
        ("A:wardenclyffe_de", "das 1,80-fache der Güte und das 2,69-fache der gespeicherten Energie", T.verluste(aus=T.formen()[0])["Q"] / T.verluste()["Q"], 1.80),
        ("A:wardenclyffe_de", "das 1,80-fache der Güte und das 2,69-fache der gespeicherten Energie", T.betrieb(aus=T.formen()[0])["W"] / T.betrieb()["W"], 2.69),
        ("A:wardenclyffe_de", "abgestrahlt werden nur noch 26 W", T.betrieb(aus=T.formen()[0])["P_S"], 26),
        ("A:wardenclyffe_de", "64 Erdrohre geben 21 mΩ", T.verluste(aus=T.formen()[1])["R_E"] * 1e3, 21),
        ("A:wardenclyffe_de", "liegt ihr Maximum bei D/h = 0,54", T.form_guete(T.SUPRA)[0], 0.54),
        ("A:wardenclyffe_de", "sie ist 68-mal größer als in der Form aus Kupfer", T.verluste(aus=T.formen()[1])["Q"] / T.verluste(aus=T.formen()[0])["Q"], 68),
        ("A:wardenclyffe_de", "dann laufen 895 A in der Wendel", T.betrieb(P=T.grenzleistung(T.formen()[1])[0], aus=T.formen()[1])["I"], 895),
        ("A:wardenclyffe_de", "Eine doppelt so weite Kuppel trägt 1,82 MJ und erreicht den Durchschlag bei 107 kW", T.grenzleistung({**T.formen()[1], "kuppel": 2 * T.formen()[1]["kuppel"]})[1] / 1e6, 1.82),
        ("A:wardenclyffe_de", "Eine doppelt so weite Kuppel trägt 1,82 MJ und erreicht den Durchschlag bei 107 kW", T.grenzleistung({**T.formen()[1], "kuppel": 2 * T.formen()[1]["kuppel"]})[0] / 1e3, 107),
        ("A:wardenclyffe_de", "Die Amplitude läuft in 34 s auf ihren Endwert zu, die Resonanz ist 9,5 mHz breit", T.anschwingen(T.formen()[1])[0], 34),
        ("A:wardenclyffe_de", "Die Amplitude läuft in 34 s auf ihren Endwert zu, die Resonanz ist 9,5 mHz breit", T.anschwingen(T.formen()[1])[1] * 1e3, 9.5),
        ("A:wardenclyffe_de", "(Kupfer, D = h: 0,61 s und 0,52 Hz)", T.anschwingen(T.formen()[0])[0], 0.61),
        ("A:wardenclyffe_de", "(Kupfer, D = h: 0,61 s und 0,52 Hz)", T.anschwingen(T.formen()[0])[1], 0.52),
        ("A:wardenclyffe_de", "Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m", T.betrieb(aus=T.bauform(T.form_guete()[0], groesse=2))["W"] / 1e3, 254),
        ("A:wardenclyffe_de", "Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m", T.reichweite(aus=T.bauform(T.form_guete()[0], groesse=2)), 263),
        ("A:wardenclyffe_de", "Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m", T.betrieb(aus=T.bauform(T.form_guete()[0], groesse=4))["W"] / 1e6, 1.03),
        ("A:wardenclyffe_de", "Doppelte Größe gibt 254 kJ und 263 m, vierfache 1,03 MJ und 528 m", T.reichweite(aus=T.bauform(T.form_guete()[0], groesse=4)), 528),
        ("A:wardenclyffe_de", "Die Güte des Kupfers Q_D = π·K_N·D/(2gδ) bleibt gleich", T.guete_geschlossen(T.bauform(T.form_guete()[0], groesse=4)) / T.guete_geschlossen(T.formen()[0]), 1.00),
        ("A:wardenclyffe_de", "Die Wendel aus Kupfer schwingt 526-mal, die supraleitende 655-mal schneller", T.einrasten(T.wendel(aus=T.formen()[0])["f"])[0], 526),
        ("A:wardenclyffe_de", "Die Wendel aus Kupfer schwingt 526-mal, die supraleitende 655-mal schneller", T.einrasten(T.wendel(aus=T.formen()[1])["f"])[0], 655),
        ("A:wardenclyffe_de", "gegen 86 mΩ der Erdung, 99,8 %", T.verluste()["R_E"] * 1e3, 86),
        ("A:wardenclyffe_de", "gegen 86 mΩ der Erdung, 99,8 %", 100 * T.viertelwelle_strahlung()[0], 99.8),
        ("A:wardenclyffe_de", "In 1000 km kommen 4,2 mV/m an, das 21-fache der Wendel", T.viertelwelle_strahlung()[1] * 1e3, 4.2),
        ("A:wardenclyffe_de", "In 1000 km kommen 4,2 mV/m an, das 21-fache der Wendel", T.viertelwelle_strahlung()[1] / T.signal(1e6), 21),
        ("A:wardenclyffe_de", "Supraleitend wächst die Güte um das 68-fache", T.verluste(aus=T.formen()[1])["Q"] / T.verluste(aus=T.formen()[0])["Q"], 68),
        ("A:wardenclyffe_de", "reicht von 54 m über 130 m bis 454 m", T.reichweite(aus=T.formen()[1]), 454),
        ("A:wardenclyffe_de", "hat die Kuppel 1,15 nF", T.kuppel()["C"] * 1e9, 1.15),
        ("A:wardenclyffe_de", "trüge sie 31 MV", T.kuppel()["V"] / 1e6, 31),
        ("A:wardenclyffe_de", "die Ladung 0,036 C", T.kuppel()["Q"], 0.036),
        ("A:wardenclyffe_de", "die Energie 0,56 MJ", T.kuppel()["W"] / 1e6, 0.56),
        ("A:wardenclyffe_de", "Turm und Schacht (93,6 m)", T.TURM["hoehe"] + T.TURM["schacht"], 93.6),
        ("A:wardenclyffe_de", "als Viertelwelle bei 801 kHz", T.turm_schwinger()[0] / 1e3, 801),
        ("A:wardenclyffe_de", "mit 570 Windungen", T.wendel()["N"], 570),
        ("A:wardenclyffe_de", "Die Wendel hat 0,194 H", T.wendel()["L"], 0.194),
        ("A:wardenclyffe_de", "ihre Eigenkapazität ist 742 pF", T.wendel()["C_W"] * 1e12, 742),
        ("A:wardenclyffe_de", "mit der Kuppel zusammen 1,89 nF", T.wendel()["C"] * 1e9, 1.89),
        ("A:wardenclyffe_de", "Ihr Draht ist 10,7 km lang", T.wendel()["draht"] / 1e3, 10.7),
        ("A:wardenclyffe_de", "schwingt sie bei 8,31 kHz", T.wendel()["f"] / 1e3, 8.31),
        ("A:wardenclyffe_de", "schlösse die stehende Welle bei 6,98 kHz", T.wendel()["f_viertel"] / 1e3, 6.98),
        ("A:wardenclyffe_de", "k = 1/√(Q_p·Q_s) = 0,0018", T.primaerkreis()["k"], 0.0018),
        ("A:wardenclyffe_de", "1/(2kf) = 33 ms", T.primaerkreis()["t"] * 1e3, 33),
        ("A:wardenclyffe_de", "Im Primärkreis laufen 192 J um", T.primaerkreis()["W_p"], 192),
        ("A:wardenclyffe_de", "153 nF und 2,39 mH", T.primaerkreis()["C_p"] * 1e9, 153),
        ("A:wardenclyffe_de", "153 nF und 2,39 mH", T.primaerkreis()["L_p"] * 1e3, 2.39),
        ("A:wardenclyffe_de", "nur 0,72 mm tief ein", T.verluste()["delta"] * 1e3, 0.72),
        ("A:wardenclyffe_de", "Das gibt 1,60 Ω", T.verluste()["R_D"], 1.60),
        ("A:wardenclyffe_de", "160π²(h/λ)²: 3,9 mΩ", T.verluste()["R_S"] * 1e3, 3.9),
        ("A:wardenclyffe_de", "86 mΩ als untere Grenze", T.verluste()["R_E"] * 1e3, 86),
        ("A:wardenclyffe_de", "Q_s = ωL/R = 5969", T.verluste()["Q"], 5969),
        ("A:wardenclyffe_de", "ν = Q_s/π = 1900", T.betrieb()["nu"], 1900),
        ("A:wardenclyffe_de", "Die Erde hat als Kugel 709 µF", T.erde_gegenplatte()[0] * 1e6, 709),
        ("A:wardenclyffe_de", "trägt die Kuppel 5,7 mC", T.kuppel()["C"] * T.betrieb()["V"] * 1e3, 5.7),
        ("A:wardenclyffe_de", "um 8,0 V", T.erde_gegenplatte(T.kuppel()["C"] * T.betrieb()["V"])[1], 8.0),
        ("A:wardenclyffe_de", "mit 344 A durch die Rohre", T.betrieb()["I"], 344),
        ("A:wardenclyffe_de", "fließen in der Wendel 344 A", T.betrieb()["I"], 344),
        ("A:wardenclyffe_de", "mit 4,91 MV Scheitelspannung", T.betrieb()["V"] / 1e6, 4.91),
        ("A:wardenclyffe_de", "mit 4,91 MV Scheitelspannung", T.betrieb()["V_energie"] / 1e6, 4.91),
        ("A:wardenclyffe_de", "W = Q_s·P/ω = 22,9 kJ", T.betrieb()["W"] / 1e3, 22.9),
        ("A:wardenclyffe_de", "Das Feld an der Kuppel ist 475 kV/m, 16 %", T.betrieb()["E"] / 1e3, 475),
        ("A:wardenclyffe_de", "Das Feld an der Kuppel ist 475 kV/m, 16 %", 100 * T.betrieb()["E_anteil"], 16),
        ("A:wardenclyffe_de", "189 kW im Kupfer warm", T.betrieb()["P_D"] / 1e3, 189),
        ("A:wardenclyffe_de", "10,1 kW in der Erdung", T.betrieb()["P_E"] / 1e3, 10.1),
        ("A:wardenclyffe_de", "465 W strahlt der Turm ab", T.betrieb()["P_S"], 465),
        ("A:wardenclyffe_de", "| 100 | 47,4 kHz | ", T.wendel(100)["f"] / 1e3, 47.4),
        ("A:wardenclyffe_de", "| 100 | 47,4 kHz | ", T.verluste(100)["Q"], 2003),
        ("A:wardenclyffe_de", "| 100 | 47,4 kHz | ", T.betrieb(100)["V"] / 1e6, 1.19),
        ("A:wardenclyffe_de", "| 100 | 47,4 kHz | ", T.betrieb(100)["P_S"] / 1e3, 28.9),
        ("A:wardenclyffe_de", "| 200 | 23,7 kHz | ", T.wendel(200)["f"] / 1e3, 23.7),
        ("A:wardenclyffe_de", "| 200 | 23,7 kHz | ", T.verluste(200)["Q"], 3322),
        ("A:wardenclyffe_de", "| 200 | 23,7 kHz | ", T.betrieb(200)["V"] / 1e6, 2.17),
        ("A:wardenclyffe_de", "| 200 | 23,7 kHz | ", T.betrieb(200)["P_S"] / 1e3, 6.0),
        ("A:wardenclyffe_de", "| 400 | 11,8 kHz | ", T.wendel(400)["f"] / 1e3, 11.8),
        ("A:wardenclyffe_de", "| 400 | 11,8 kHz | ", T.verluste(400)["Q"], 4935),
        ("A:wardenclyffe_de", "| 400 | 11,8 kHz | ", T.betrieb(400)["V"] / 1e6, 3.74),
        ("A:wardenclyffe_de", "| 400 | 11,8 kHz | ", T.betrieb(400)["P_S"] / 1e3, 1.11),
        ("A:wardenclyffe_de", "| 570 | 8,31 kHz | ", T.wendel(570)["f"] / 1e3, 8.31),
        ("A:wardenclyffe_de", "| 570 | 8,31 kHz | ", T.verluste(570)["Q"], 5969),
        ("A:wardenclyffe_de", "| 570 | 8,31 kHz | ", T.betrieb(570)["V"] / 1e6, 4.91),
        ("A:wardenclyffe_de", "| 570 | 8,31 kHz | ", T.betrieb(570)["P_S"] / 1e3, 0.47),
        ("A:wardenclyffe_de", "in 100 km 2,0 mV/m", T.signal(1e5) * 1e3, 2.0),
        ("A:wardenclyffe_de", "in 1000 km 0,20 mV/m", T.signal(1e6) * 1e3, 0.2),
        ("A:wardenclyffe_de", "in 0,1335 s einmal um die Erde", T.erdschale()[0], 0.1335),
        ("A:wardenclyffe_de", "gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz", T.erdschale()[1][0], 7.49),
        ("A:wardenclyffe_de", "gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz", T.erdschale()[1][1], 14.98),
        ("A:wardenclyffe_de", "gibt das 7,49 Hz und die Vielfachen 14,98 und 22,47 Hz", T.erdschale()[1][2], 22.47),
        ("A:wardenclyffe_de", "10,59, 18,34, 25,94 und 33,49 Hz", T.erdschale()[2][0], 10.59),
        ("A:wardenclyffe_de", "10,59, 18,34, 25,94 und 33,49 Hz", T.erdschale()[2][1], 18.34),
        ("A:wardenclyffe_de", "10,59, 18,34, 25,94 und 33,49 Hz", T.erdschale()[2][2], 25.94),
        ("A:wardenclyffe_de", "10,59, 18,34, 25,94 und 33,49 Hz", T.erdschale()[2][3], 33.49),
        ("A:wardenclyffe_de", "Der Turm schwingt 784-mal schneller", T.einrasten()[0], 784),
        ("A:wardenclyffe_de", "die Kopplung müsste halb so schnell sein", 2 * T.einrasten()[1], 1),
        ("A:wardenclyffe_de", "bei 10,59 Hz schwänge, bräuchte sie eine Spule von 1,96·10⁵ H", T.spule_fuer(T.erdschale()[2][0]) / 1e5, 1.96),
        ("A:wardenclyffe_de", "Bei 200 kW kommen auf der Erdoberfläche 3,9·10⁻¹⁰ W/m² an", T.schale_fuer_energie(P=T.AUSLEGUNG["P"])[0] * 1e10, 3.9),
        ("A:wardenclyffe_de", "müsste der Turm 5,1·10¹⁴ W einspeisen", T.schale_fuer_energie()[1] / 1e14, 5.1),
        ("A:wardenclyffe_de", "ist g = 9,82 m/s²", T.regel_am_turm()[0], 9.82),
        ("A:wardenclyffe_de", "s = √(g/a₀) = 3,04·10⁵", T.regel_am_turm()[1] / 1e5, 3.04),
        ("A:wardenclyffe_de", "lg q = −1,3·10⁵", T.regel_am_turm()[2] / 1e5, -1.3),
        ("A:wardenclyffe_de", "schwingt die Kuppel bei 8,31 kHz mit 4,91 MV", T.wendel()["f"] / 1e3, 8.31),
        ("A:wardenclyffe_de", "schwingt die Kuppel bei 8,31 kHz mit 4,91 MV", T.betrieb()["V"] / 1e6, 4.91),
        ("A:wardenclyffe_de", "der Turm schwingt 784-mal zu schnell", T.einrasten()[0], 784),
        ("01", "Erde Δs = 1,478 mm", K.faltung(K.M_earth) * 1e3, 1.478),
        ("01", "Sonne Δs = 492,2 m", K.faltung(K.M_sun), 492.2),
        ("01", "Erde 8,87 mm", K.grenzradius(K.M_earth) * 1e3, 8.87),
        ("01", "Sonne 2,95 km", K.grenzradius(K.M_sun) / 1e3, 2.95),
        ("02", "5,906·10⁻³⁹", K.kopplung_grav(K.m_p) * 1e39, 5.906),
        ("02", "1/137,036", 1 / K.alpha, 137.036),
        ("02", "0,529·10⁻¹⁰ m", K.bohr_radius() * 1e10, 0.529),
        ("02", "6,58·10¹⁵ Hz", K.alpha * K.c / (2 * math.pi * K.bohr_radius()) / 1e15, 6.58),
        ("02", "Erde 29,78 km/s", K.bahntempo(K.M_sun, K.AU) / 1e3, 29.78),
        ("02", "1,853 M☉", K.sternmasse() / K.M_sun, 1.853),
        ("03", "Δs = 197,4 AE", K.faltung(6e10 * K.M_sun) / K.AU, 197.4),
        ("03", "ν = 1,582", K.verstaerkung(K.a0), 1.582),
        ("03", "r* = 8,88 kpc", K.schwellenradius(6e10 * K.M_sun) / K.kpc, 8.88),
        ("03", "170,4 km/s", K.v_flach(6e10 * K.M_sun) / 1e3, 170.4),
        ("03", "Λ = 13,7 Gpc", K.schliesslaenge() / K.Gpc, 13.7),
        ("03", "89,6 Mrd. Jahre", K.rueckkehrzeit() / K.Gyr, 89.6),
        ("03", "L = 27,9 Gpc", None, 27.9),
        ("03", "= 5,36", K.Q_kosm, 5.36),
        ("03", "q = 0,843", K.Q_kosm / (1 + K.Q_kosm), 0.843),
        ("03", "Θ = 2πa₀/(cH₀) = 1,020", K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc), 1.020),
        ("03", "Ω_s = 0,228", K.vakuumteil(H=K.atmen_heute()[0] * 1e3 / K.Mpc), 0.228),
        ("03", "H₀·√Ω_s = 32,1 km/s/Mpc", K.H0_kms * math.sqrt(K.vakuumteil()), 32.1),
        ("03", "k_s = π³/30 = 1,034", math.pi**3 / 30, 1.034),
        ("03", "Ω_DE/Ω_s = 3,03", (1 - K.Omega_b - K.Omega_c) / K.vakuumteil(), 3.03),
        ("03", "a₀ = 1,06·10⁻¹⁰", K.a0 * math.sqrt((1 - K.Omega_b - K.Omega_c) / (3 * K.vakuumteil())) * 1e10, 1.06),
        ("02", "ω₃ = 2,018 der Beiwert der Form n = 3", K.lane_emden(3.0)[1], 2.018),
        ("02", "m_Pl³/(μ_e·m_u)² = 1,456 M☉", K.grenzmasse() / K.M_sun, 1.456),
        ("02", "im Tal nach der Tröpfchenformel Z = 26,7, zwischen ⁵⁶Fe und ⁶²Ni", float(round(__import__("kerne").ende_brennleiter(__import__("kerne").TROPFEN)[0], 1) == 26.7 and __import__("kerne").bekannte()["m_min"] == (26, 56) and __import__("kerne").bekannte()["b_max"] == (28, 62)), 1.0),
        ("02", "höchstens kT = 0,096·(μ/μ_e)·m_ec²", K.decke_entartet(1.0)[0] / (K.m_e * K.c**2), 0.096),
        ("02", "für Kohlenstoff und Sauerstoff 5,0·10⁸ K", K.decke_entartet(0.875)[0] / K.k_B / 1e8, 5.0),
        ("03", "an den blauen, die noch Sterne bilden, χ² = 21 in 15 Bins", D.kids_typ_reste("blau")[3], 21),
        ("F:lesungen", "Sechs gasreiche (Mancera Piña", float(D.udg_reste()[0]), 6),
        ("F:lesungen", "äußerer Punkt bei a_N = 0,018 bis 0,039·a₀", D.udg_reste()[3], 0.018),
        ("F:lesungen", "äußerer Punkt bei a_N = 0,018 bis 0,039·a₀", D.udg_reste()[4], 0.039),
        ("F:lesungen", "drehen wie Newton allein, im Mittel +0,01 dex", D.udg_reste()[2], 0.01),
        ("F:lesungen", "Die Regel verlangt 0,40 dex mehr Geschwindigkeit", -D.udg_reste()[1], 0.40),
        ("F:lesungen", "einen 2,5-mal kleineren Sinus der Neigung", 10**(-D.udg_reste()[1]), 2.5),
        ("03", "die Zelle erreicht die letzte Streufläche bei Θ = 1,020, für jedes a₀", K.theta_zelle(strahlung=False), 1.020),
        ("03", "Mit der Strahlung ist es Θ = 1,024 bis 1,030 (a₀ = 0,9 bis 1,2)", K.theta_zelle(0.9e-10), 1.024),
        ("03", "Mit der Strahlung ist es Θ = 1,024 bis 1,030 (a₀ = 0,9 bis 1,2)", K.theta_zelle(1.2e-10), 1.030),
        ("03", "das Mittel aus Neigung (parametrisch) und Schwelle (transversal), a = √(a_N·a₀)", K.verstaerkung(1e-6 * K.a0) * 1e-6 * K.a0 / math.sqrt(1e-6 * K.a0 * K.a0), 1.0),
        ("03", "Ihr Produkt ist fest, v_flach², gleich sind beide bei r*", float(max(abs(math.sqrt(math.sqrt(K.G * 6e10 * K.M_sun / r) * math.sqrt(K.a0 * r)) / K.v_flach(6e10 * K.M_sun) - 1) for r in (K.kpc, 30 * K.kpc, 300 * K.kpc)) < 1e-12), 1),
        ("03", "Ihr Produkt ist fest, v_flach², gleich sind beide bei r*", (lambda r: math.sqrt(K.G * 6e10 * K.M_sun / r) / math.sqrt(K.a0 * r))(K.schwellenradius(6e10 * K.M_sun)), 1),
        ("03", "1,012 (Zelle und innen, L/ℓ₀ = Θ^{2/3}/(π·Ω_m)^{1/3})", K.schleifen_ganzes()[1], 1.012),
        ("03", "1,012 (Zelle und innen, L/ℓ₀ = Θ^{2/3}/(π·Ω_m)^{1/3})", (lambda th, a: K.zelle_laenge(a0_=a, rho_m=3 * (2 * math.pi * a / (K.c * th))**2 * (1 - math.pi**2 / 15 * th**2) / (8 * math.pi * K.G)) * a / K.c**2)(K.schleifen_ganzes()[1], 0.9e-10), 1),
        ("03", "und 1,023 (Zelle und Horizont, L·H/2πc = (π·Θ·Ω_m)^{−1/3})", K.schleifen_ganzes()[2], 1.023),
        ("03", "und 1,023 (Zelle und Horizont, L·H/2πc = (π·Θ·Ω_m)^{−1/3})", (lambda th, a: K.zelle_laenge(a0_=a, rho_m=3 * (2 * math.pi * a / (K.c * th))**2 * (1 - math.pi**2 / 15 * th**2) / (8 * math.pi * K.G)) * (2 * math.pi * a / (K.c * th)) / (2 * math.pi * K.c))(K.schleifen_ganzes()[2], 1.2e-10), 1),
        ("03", "liegt das Ganze in diesem Fenster von 2 %", 100 * (K.schleifen_ganzes()[2] - K.schleifen_ganzes()[0]), 2),
        ("03", "dauert das Fenster 0,68 Mrd. Jahre", K.zeit_zwischen_theta(1.0, K.schleifen_ganzes()[2]) / 3.15576e16, 0.68),
        ("03", "Es öffnete sich vor 0,58 Mrd. Jahren (z = 0,041)", K.zeit_zwischen_theta(1.0, K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc)) / 3.15576e16, 0.58),
        ("03", "Es öffnete sich vor 0,58 Mrd. Jahren (z = 0,041)", ((15 / math.pi**2 - 1) / (15 / (math.pi**2 * K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc)**2) - 1))**(1 / 3) - 1, 0.041),
        ("03", "und schließt sich in 0,10 Mrd. Jahren", K.zeit_zwischen_theta(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc), K.schleifen_ganzes()[2]) / 3.15576e16, 0.10),
        ("03", "Im besten Punkt ist L/D_S = 0,9984 ± 0,0007", _zelle_kreis()[0], 0.9984),
        ("03", "Im besten Punkt ist L/D_S = 0,9984 ± 0,0007", _zelle_kreis()[1], 0.0007),
        ("03", "also α = 3,2° (2,4° bis 3,9°)", _zelle_kreis()[2], 3.2),
        ("03", "also α = 3,2° (2,4° bis 3,9°)", _zelle_kreis()[3], 2.4),
        ("03", "also α = 3,2° (2,4° bis 3,9°)", _zelle_kreis()[4], 3.9),
        ("03", "(vier Zufallsfelder, M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉", _haufen_einfall()["m_min"] / 1e15, 1.0),
        ("03", "(vier Zufallsfelder, M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉", _haufen_einfall()["m_max"] / 1e15, 3.5),
        ("03", "je über 31 000 Körper in R₂₀₀), liest das Profil", float(_haufen_einfall()["koerper"] // 1000 * 1000), 31000),
        ("03", "liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2)", _haufen_einfall()["c"], 4.0),
        ("03", "liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2)", _haufen_einfall()["c_min"], 3.3),
        ("03", "liest das Profil wie die Haufen c = 4,0 (3,3 bis 5,2)", _haufen_einfall()["c_max"], 5.2),
        ("03", "Mit ihm ist der Rest an X-COP im Kern −0,02 dex", D.xcop_staub(_haufen_einfall()["c"])[2], -0.02),
        ("03", "mit der Konzentration aus dem Einfall Q = 5,6 ± 0,2", _haufen_einfall()["q"], 5.6),
        ("03", "mit der Konzentration aus dem Einfall Q = 5,6 ± 0,2", _haufen_einfall()["q_fehler"], 0.2),
        ("F:lesungen", "(SH0ES 2022, 3492 Gleichungen) liest H₀ = 73,0 ± 1,0", D.sh0es_leiter()[0], 73.0),
        ("F:lesungen", "(SH0ES 2022, 3492 Gleichungen) liest H₀ = 73,0 ± 1,0", D.sh0es_leiter()[1], 1.0),
        ("03", "kreuzt die anderen zuerst in der Mitte, bei ω·t = π/(2√2) = 1,11", K.falte_phase(), 1.11),
        ("03", "kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11", _falte_einfall()[0], 1.10),
        ("03", "kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11", _falte_einfall()[1], 1.08),
        ("03", "kreuzt sie im Median bei 1,10 (Viertel 1,08 bis 1,11", _falte_einfall()[2], 1.11),
        ("F:lesungen", "liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36°", D.udg_neigung()[0], 31),
        ("F:lesungen", "liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36°", D.udg_neigung()[1], 28),
        ("F:lesungen", "liest die Scheibe der Sterne unabhängig vom Gas i = 31° (28° bis 36°", D.udg_neigung()[2], 36),
        ("03", "HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5", D.hiflugcs_gang()[0], 9.8),
        ("03", "HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5", D.hiflugcs_gang()[1], 6.5),
        ("03", "X-COP liest c = 4,2", _haufen_staub()[0], 4.2),
        ("F:lesungen", "An KiDS-1000 (ξ₊ und ξ₋, 225 Werte", float(_scherung()["modell"][0][3]), 225),
        ("F:lesungen", "liegen Modell und ΛCDM gleich, χ² = 307,8", _scherung()["modell"][0][0], 307.8),
        ("F:lesungen", "liegen Modell und ΛCDM gleich, χ² = 307,8", _scherung()["lcdm"][0][0], 307.8),
        ("F:lesungen", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][1][0], 0.73),
        ("F:lesungen", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][0][2], 0.82),
        ("F:lesungen", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][0][0] - _scherung()["modell"][1][1], 38),
        ("F:lesungen", "Bis M₂₀₀ = 2,5·10¹⁴ M☉ liegt R₂₀₀ jenseits der Naht", _uebergabe()[3] / 1e14, 2.5),
        ("03", "liest HIFLUGCS Q = 5,5 ± 0,3", D.haufen_q()[0], 5.5),
        ("03", "liest HIFLUGCS Q = 5,5 ± 0,3", D.haufen_q_fehler(), 0.3),
        ("03", "(Rest im Kern −0,04 dex, in der Form des Gases +0,41)", D.xcop_staub(_haufen_staub()[0])[2], -0.04),
        ("03", "(Rest im Kern −0,04 dex, in der Form des Gases +0,41)", D.xcop_staub()[2], 0.41),
        ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", _haufen_staub()[1], 5.5),
        ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", D.hiflugcs_staub(_haufen_staub()[1])[0], 447),
        ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", D.hiflugcs_staub()[0], 601),
        ("03", "innerhalb 300 kpc −0,04 statt +0,21 dex", D.hiflugcs_staub(_haufen_staub()[1])[2], -0.04),
        ("03", "innerhalb 300 kpc −0,04 statt +0,21 dex", D.hiflugcs_staub()[2], 0.21),
        ("03", "im Kern das 6,8-fache des Gases tragen", D.xcop_sterne()[2][0], 6.8),
        ("F:lesungen", "das Maß des Ganzen 1 + Q, das dort dasselbe gibt, 0,39 dex", 0.5 * math.log10(1 + K.q_je_radiant()) - D.udg_reste()[2], 0.39),
        ("F:lesungen", "(c = 4,0, R₂₀₀ um 45 kpc)", _haufen_einfall()["c"], 4.0),
        ("F:lesungen", "(c = 4,0, R₂₀₀ um 45 kpc)", __import__("statistics").median([K.uebergabe_linse((1 + K.q_je_radiant()) * 10**m * K.M_sun, 68.5e3 / K.Mpc, K.q_je_radiant())[0] / K.kpc for m in (9.21, 9.21, 9.36, 9.05, 9.25, 9.17)]), 45),
        ("F:lesungen", "und sie verlangt dort 0,13 dex mehr", -D.udg_reste(einfall_c=_haufen_einfall()["c"])[1], 0.13),
        ("F:lesungen", "Die gasreichen isolierten Zwerge liegen dagegen 0,37 dex über dieser Form", D.zwerge_reste("gasreich", False, einfall=(_haufen_einfall()["c"], K.q_je_radiant(), 68.5e3 / K.Mpc))[1], 0.37),
        ("03", "χ² = 21 in 15 Bins", float(D.kids_typ_reste("blau")[0]), 15),
        ("03", "mit 0,003·a₀ im Betrag schon 534", D.kids_typ_reste("blau", g_e=0.003 * K.a0)[3], 534),
        ("03", "χ² = 27483 statt 503", _kids()[2], 27483),
        ("03", "Naht (30 bis 66 kpc)", K.naht_radius(D.kids_isoliert()[0]["M_b"] * K.M_sun, K.q_je_radiant()) / K.kpc, 30),
        ("03", "Naht (30 bis 66 kpc)", K.naht_radius(D.kids_isoliert()[3]["M_b"] * K.M_sun, K.q_je_radiant()) / K.kpc, 66),
        ("03", "Die 22 klassischen Satelliten", float(_zwerge("außen", False)[0]), 22),
        ("03", "liegen ohne Wirtsfeld +0,15 ± 0,04 dex über der Regel", _zwerge("außen", False)[1], 0.15),
        ("03", "liegen ohne Wirtsfeld +0,15 ± 0,04 dex über der Regel", _zwerge("außen", False)[2], 0.04),
        ("03", "mit ihm im Betrag +0,46 ± 0,05", _zwerge("außen", True)[1], 0.46),
        ("03", "mit ihm im Betrag +0,46 ± 0,05", _zwerge("außen", True)[2], 0.05),
        ("F:galaxien", "gasreichen isolierten 0,00 ± 0,04 dex", _zwerge("gasreich", False)[1], 0.00),
        ("F:galaxien", "gasreichen isolierten 0,00 ± 0,04 dex", _zwerge("gasreich", False)[2], 0.04),
        ("F:galaxien", "die gasarmen isolierten +0,13 ± 0,07", _zwerge("gasarm", False)[1], 0.13),
        ("F:galaxien", "die gasarmen isolierten +0,13 ± 0,07", _zwerge("gasarm", False)[2], 0.07),
        ("F:galaxien", "die 36 klassischen Satelliten +0,11 ± 0,03", float(_zwerge("klassisch", False)[0]), 36),
        ("F:galaxien", "die 36 klassischen Satelliten +0,11 ± 0,03", _zwerge("klassisch", False)[1], 0.11),
        ("F:galaxien", "die 36 klassischen Satelliten +0,11 ± 0,03", _zwerge("klassisch", False)[2], 0.03),
        ("F:galaxien", "die ultraschwachen +0,37 dex", _zwerge("ultraschwach", False)[1], 0.37),
        ("F:galaxien", "liest Υ_V = 1,96 für die gasreichen", D.zwerge_yv("gasreich"), 1.96),
        ("F:galaxien", "5,2 für die klassischen Satelliten", D.zwerge_yv("klassisch"), 5.2),
        ("F:galaxien", "weder das Verhältnis der Dichten zum Wirt", D.zwerge_gezeit(), 0.0),
        ("F:galaxien", "noch das Perizentrum der Bahn", D.satelliten_perizentrum()[0], 0.0),
        ("F:galaxien", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", D.zwerge_gas_beziehung()[0], 2.09),
        ("F:galaxien", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", -D.zwerge_gas_beziehung()[1], 0.284),
        ("F:galaxien", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", float(D.zwerge_gas_beziehung()[2]), 188),
        ("F:galaxien", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", _zwerge_g("klassisch")[1], 0.011),
        ("F:galaxien", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", _zwerge_g("klassisch")[2], 0.028),
        ("F:galaxien", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", D.zwerge_yv("klassisch", True), 2.49),
        ("F:galaxien", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", float(D.zwerge_ausbeute()[2]), 7),
        ("F:galaxien", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", D.zwerge_ausbeute()[0], -1.03),
        ("F:galaxien", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", D.zwerge_ausbeute()[1], 0.20),
        ("F:galaxien", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", float(_zwerge_m("klassisch")[0]), 27),
        ("F:galaxien", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", _zwerge_m("klassisch")[1], -0.03),
        ("F:galaxien", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", _zwerge_m("klassisch")[2], 0.03),
        ("F:galaxien", "die 27 ultraschwachen +0,07 ± 0,04", float(_zwerge_m("ultraschwach")[0]), 27),
        ("F:galaxien", "die 27 ultraschwachen +0,07 ± 0,04", _zwerge_m("ultraschwach")[1], 0.07),
        ("F:galaxien", "die 27 ultraschwachen +0,07 ± 0,04", _zwerge_m("ultraschwach")[2], 0.04),
        ("F:galaxien", "die gasarmen isolierten +0,05 ± 0,08", _zwerge_m("gasarm")[1], 0.05),
        ("F:galaxien", "die gasarmen isolierten +0,05 ± 0,08", _zwerge_m("gasarm")[2], 0.08),
        ("F:galaxien", "ebenso −0,03 und +0,06", _zwerge_m("klassisch", True)[1], -0.03),
        ("F:galaxien", "ebenso −0,03 und +0,06", _zwerge_m("ultraschwach", True)[1], 0.06),
        ("F:galaxien", "(36 Satelliten der Milchstraße mit Eigenbewegung)", float(D.satelliten_perizentrum()[1]), 36),
        ("03", "X-COP (12 Haufen", float(len(D.xcop())), 12),
        ("03", "a_N = 0,036·a₀ (0,031 bis 0,065)", _xcop("aN_200")[0], 0.036),
        ("03", "a_N = 0,036·a₀ (0,031 bis 0,065)", _xcop("aN_200")[1], 0.031),
        ("03", "a_N = 0,036·a₀ (0,031 bis 0,065)", _xcop("aN_200")[2], 0.065),
        ("03", "gegen 0,029·a₀. Dort", math.log(K.Q_kosm / (1 + K.Q_kosm))**2, 0.029),
        ("03", "M/M_Gas = 7,9 (4,7 bis 10,5)", _xcop("M_Mgas_200")[0], 7.9),
        ("03", "M/M_Gas = 7,9 (4,7 bis 10,5)", _xcop("M_Mgas_200")[1], 4.7),
        ("03", "M/M_Gas = 7,9 (4,7 bis 10,5)", _xcop("M_Mgas_200")[2], 10.5),
        ("03", "gegen 1 + Q = 6,28, die Bilanz", 1 + K.q_je_radiant(), 6.28),
        ("03", "Sternen von einem Viertel des Gases", _xcop("M_Mgas_200")[0] / (1 + K.q_je_radiant()) - 1, 0.3),
        ("03", "in R₂₀₀ das 1,7-fache des Gases", _xcop("verlangt_200")[0], 1.7),
        ("03", "(0,1·R₅₀₀) das 12-fache", _xcop("verlangt_kern")[0], 12),
        ("03", "über alle 108 Punkte +0,03 dex auf den Daten (χ² = 601)", D.hiflugcs_reste("form")[1], 0.03),
        ("03", "über alle 108 Punkte +0,03 dex auf den Daten (χ² = 601)", D.hiflugcs_reste("form")[2], 601),
        ("03", "über alle 108 Punkte +0,03 dex", float(D.hiflugcs_reste("form")[0]), 108),
        ("03", "die örtliche Regel 0,17 dex darunter (χ² = 4173)", D.hiflugcs_reste("regel")[1], 0.17),
        ("03", "die örtliche Regel 0,17 dex darunter (χ² = 4173)", D.hiflugcs_reste("regel")[2], 4173),
        ("03", "das Maß 1 + Q dort nur das 2,9-fache", D.xcop_form(), 2.9),
        ("03", "R_Q = 52,0 kpc", K.naht_radius(6e10 * K.M_sun) / K.kpc, 52.0),
        ("03", "a_N = 0,029·a₀", math.log(K.Q_kosm / (1 + K.Q_kosm))**2, 0.029),
        ("03", "χ = Q·(−ln q) = 0,917", K.Q_kosm * -math.log(K.Q_kosm / (1 + K.Q_kosm)), 0.917),
        ("03", "r*/Λ = 0,971", K.form_zelle(), 0.971),
        ("09", "μ_B = 9,274·10⁻²⁴ J/T", K.magneton(K.m_e) / 1e-24, 9.274),
        ("09", "m_p/m_e = 1836,2", K.m_p / K.m_e, 1836.2),
        ("09", "beim Elektron 2,9·10⁴⁴", K.spinzahl_teilchen(K.m_e) / 1e44, 2.9),
        ("09", "m_Pl/√2 = 1,54·10⁻⁸ kg", K.spinzahl_teilchen(K.m_Pl / math.sqrt(2)) * K.m_Pl / math.sqrt(2) / 1e-8, 1.54),
        ("09", "a_J = 3,27 m", K.verwindung(K.J_ERDE, K.M_earth), 3.27),
        ("09", "3Δs = 4,4 mm", 3 * K.faltung(K.M_earth) * 1e3, 4.4),
        ("09", "also j = 738", K.verwindung(K.J_ERDE, K.M_earth) / (3 * K.faltung(K.M_earth)), 738),
        ("F:galaxien", "Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43)", float(D.zwerge_form()[0]), 65),
        ("F:galaxien", "Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43)", float(D.zwerge_form()[1]), -0.1),
        ("F:galaxien", "Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43)", float(D.zwerge_form()[3]), 0.36),
        ("F:galaxien", "Die Elliptizität der 65 Satelliten hängt nicht am Verhältnis ihrer Dichte zur Gezeitendichte mit der Regel (Rangkorrelation −0,10), und sie sind nicht gestreckter als die isolierten Zwerge (Median 0,36 gegen 0,43)", float(D.zwerge_form()[4]), 0.43),
        ("07", "gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern)", 100 * _scheibe_austausch()[0], 4),
        ("07", "gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern)", 100 * _scheibe_austausch()[1], 29),
        ("07", "gibt 4 bis 29 % ihres Drehimpulses ab (sieben Läufe mit 20 000 bis 80 000 Körpern)", float(_scheibe_austausch()[3]), 7),
        ("07", "in 3 Gyr auf das 2- bis 3-fache ihres Halbmassenradius", _scheibe_formen()[0], 2),
        ("07", "in 3 Gyr auf das 2- bis 3-fache ihres Halbmassenradius", _scheibe_formen()[1], 3),
        ("07", "ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton", _scheibe_formen()[2], 1.17),
        ("07", "ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton", _scheibe_formen()[3], 0.96),
        ("07", "ihr Halbmassenradius geht von 1,17·r* auf 0,96·r*, und ihr Drehimpuls bleibt auf 0,2 % erhalten wie nach Newton", 100 * _scheibe_formen()[4], 0.2),
        ("03", "Über die ganze Schicht gemittelt liest die Bahn das senkrechte Feld für jedes Profil gleichverteilt zwischen 0 und 2πGΣ", round(1e6 * _schicht_gleichverteilt()), 0),
        ("03", "innerhalb 2R_d im Median um 0,03 dex der Beschleunigung, 4 % des Tempos", -_sparc_bahn_feld()["schicht_innen"], 0.03),
        ("03", "innerhalb 2R_d im Median um 0,03 dex der Beschleunigung, 4 % des Tempos", 100 * (1 - 10**(_sparc_bahn_feld()["schicht_innen"] / 2)), 4),
        ("03", "4 % des Tempos (flächenschwach 0,05, flächenhell 0,01)", -_sparc_bahn_feld()["schicht_lsb"], 0.05),
        ("03", "4 % des Tempos (flächenschwach 0,05, flächenhell 0,01)", -_sparc_bahn_feld()["schicht_hsb"], 0.01),
        ("F:galaxien", "zählte die Massenfunktion bei W50 = 50 bis 63 km/s das 218- und 508-fache der Galaxien", _breiten()[0], 218),
        ("F:galaxien", "zählte die Massenfunktion bei W50 = 50 bis 63 km/s das 218- und 508-fache der Galaxien", _breiten()[1], 508),
        ("F:galaxien", "bei 355 bis 447 km/s das 2,0- und 4,4-fache", _breiten()[2], 2.0),
        ("F:galaxien", "bei 355 bis 447 km/s das 2,0- und 4,4-fache", _breiten()[3], 4.4),
        ("F:galaxien", "der Linie einer Scheibe an der Schwelle mit der gelesenen Streuung von 12 km/s", _gas_streuung()["leicht"], 12),
        ("F:galaxien", "unter W50 = 100 km/s auf das 0,4- bis 4,6-fache", min(_breiten_schnitt()[0]), 0.4),
        ("F:galaxien", "unter W50 = 100 km/s auf das 0,4- bis 4,6-fache", max(_breiten_schnitt()[0]), 4.6),
        ("F:galaxien", "bei 100 bis 250 km/s bleibt das 7- bis 21-fache", min(_breiten_schnitt()[1]), 7),
        ("F:galaxien", "bei 100 bis 250 km/s bleibt das 7- bis 21-fache", max(_breiten_schnitt()[1]), 21),
        ("F:galaxien", "verlangte ALFALFA ε = 0,02 bis 0,05", min(_breiten_anteil()[0]), 0.02),
        ("F:galaxien", "verlangte ALFALFA ε = 0,02 bis 0,05", max(_breiten_anteil()[0]), 0.05),
        ("F:galaxien", "das 20- bis 56-fache ihrer Baryonen als nie abgekühltes Gas", min(_breiten_anteil()[5]), 20),
        ("F:galaxien", "das 20- bis 56-fache ihrer Baryonen als nie abgekühltes Gas", max(_breiten_anteil()[5]), 56),
        ("F:galaxien", "hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32", _breiten_anteil()[3], 21),
        ("F:galaxien", "hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32", _breiten_anteil()[1], 340),
        ("F:galaxien", "hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32", _breiten_anteil()[2], 702),
        ("F:galaxien", "hebt es χ² von 21 auf 340 bis 702, schon das Einfache auf 32", _breiten_anteil()[4], 32),
        ("F:galaxien", "liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache", _breiten_hipass()[0], 177),
        ("F:galaxien", "liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache", _breiten_hipass()[1], 101),
        ("F:galaxien", "liest HIPASS (Südhimmel) bei 63 bis 79 km/s das 177-fache, ALFALFA dort das 101- und 304-fache", _breiten_hipass()[2], 304),
        ("03", "In Bereichen schwerer als M liegt der Anteil erfc(δ_c/√(2S))", round(1e3 * _breiten()[4]), 0),
        ("03", "mit 0,90 bis 0,92 des Tempos in der Ebene", _polring()[0], 0.90),
        ("03", "mit 0,90 bis 0,92 des Tempos in der Ebene", _polring()[1], 0.92),
        ("03", "nach Newton mit den Baryonen allein mit 0,86 bis 0,88", _polring()[2], 0.86),
        ("03", "nach Newton mit den Baryonen allein mit 0,86 bis 0,88", _polring()[3], 0.88),
        ("03", "Bei 8·R_d sind es 0,96 bis 0,97", _polring()[4], 0.96),
        ("03", "Bei 8·R_d sind es 0,96 bis 0,97", _polring()[5], 0.97),
        ("03", "Polarringe laufen langsamer als die Scheibe ihrer Galaxie", round(100 * _polring()[6]), 0),
        ("F:galaxien", "fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11}", D.rc100_groesse(_saett_punkt())[0], -0.32),
        ("F:galaxien", "fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11}", D.rc100_groesse(_saett_punkt())[1], 0.11),
        ("F:galaxien", "fällt R_e bei gleicher Masse wie H^{−0,32 ± 0,11}", round(10 * abs(D.rc100_groesse(_saett_punkt(), 2000)[6] - D.rc100_groesse(_saett_punkt())[0])), 0),
        ("F:galaxien", "Ohne Gang läge es 2,8 σ daneben", -D.rc100_groesse(_saett_punkt())[0] / D.rc100_groesse(_saett_punkt())[1], 2.8),
        ("F:galaxien", "R_d ∝ R₂₀₀ ∝ H^{−2/3}, 3,0 σ", (2 / 3 + D.rc100_groesse(_saett_punkt())[0]) / D.rc100_groesse(_saett_punkt())[1], 3.0),
        ("F:lokale_gruppe", "ruht nach 13,8 Gyr die Schale bei 0,89 Mpc", _lokaler_fluss()["R0"], 0.89),
        ("F:lokale_gruppe", "die 53 Zwerge zwischen 0,6 und 3 Mpc", float(_lokaler_fluss()["zahl"]), 53),
        ("F:lokale_gruppe", "liegen im Median 10 km/s darunter (4 bis 17)", -_lokaler_fluss()["rest"], 10),
        ("F:lokale_gruppe", "liegen im Median 10 km/s darunter (4 bis 17)", -_lokaler_fluss()["oben"], 4),
        ("F:lokale_gruppe", "liegen im Median 10 km/s darunter (4 bis 17)", -_lokaler_fluss()["unten"], 17),
        ("F:lokale_gruppe", "ruhte erst die Schale bei 2,7 Mpc", _lokaler_fluss()["R0_regel"], 2.7),
        ("F:lokale_gruppe", "ruht nach 13,8 Gyr die Schale bei 0,89 Mpc", K.weltalter(lambda a: 3 * K.dichte_vakuumteil()) / K.Gyr, 13.8),
        ("F:lokale_gruppe", "um isolierte Galaxien bei 1,9 Mpc 0,03 dex unter ihr", -_licht_einfall()[0], 0.03),
        ("F:lokale_gruppe", "die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc", ((1 + K.q_je_radiant()) * K.G * _bereiche()[0] / K.groesste_umkehr(K.M_sun)[1]**2)**(1 / 3) / K.Mpc, 0.83),
        ("F:lokale_gruppe", "die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc", ((1 + K.q_je_radiant()) * K.G * sum(_bereiche()) / K.groesste_umkehr(K.M_sun)[1]**2)**(1 / 3) / K.Mpc, 1.23),
        ("F:lokale_gruppe", "die Milchstraße mit ihren Satelliten bis 0,83 Mpc, die Lokale Gruppe bis 1,23 Mpc", float(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: (1 + K.q_je_radiant()) * K.G * sum(_bereiche()) / x**2 - K.groesste_umkehr(K.M_sun)[1]**2 * x, 0.1 * K.Mpc, 10 * K.Mpc) / K.Mpc), 1.23),
        ("F:galaxien", "an SPARC λ' ∝ M_b^{−0,10 ± 0,02}", D.sparc_drehung_gang(D.sparc_eichung(K.kpc), _cmb()["saett"]["punkt"]["H0"], K.q_je_radiant())[0], -0.10),
        ("F:galaxien", "an SPARC λ' ∝ M_b^{−0,10 ± 0,02}", D.sparc_drehung_gang(D.sparc_eichung(K.kpc), _cmb()["saett"]["punkt"]["H0"], K.q_je_radiant())[1], 0.02),
        ("F:galaxien", "wie die Schwankung σ(M) des Spektrums über dieselben Bereiche (−0,11)", _sigma_neigung(), -0.11),
        ("F:galaxien", "der Halbmassenradius geht von 0,59 auf 0,46·r*", _scheibe_lage()[0], 0.59),
        ("F:galaxien", "der Halbmassenradius geht von 0,59 auf 0,46·r*", _scheibe_lage()[1], 0.46),
        ("F:galaxien", "an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r*", _scheibe_lage()[2], 1.17),
        ("F:galaxien", "an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r*", _scheibe_lage()[3], 0.96),
        ("F:galaxien", "an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r*", _scheibe_lage()[4], 2.35),
        ("F:galaxien", "an der Schwelle von 1,17 auf 0,96·r*, dünner von 2,35 auf 2,09·r*", _scheibe_lage()[5], 2.09),
        ("07", "für die Kuzmin-Scheibe exakt, über der das Feld das eines Punktes im Spiegelbild unter ihr ist: a = ν(√(g_R² + (2πGΣ)²)/a₀)·g_R", round(100 * _achsen_zweiter_weg()[0]), 0),
        ("07", "In einer dünnen Scheibe liest die Feldform ν am vollen Betrag der Neigung an der Oberfläche der Schicht", round(1000 * abs(_achsen_zweiter_weg()[1] - _achsen_zweiter_weg()[2])), 0),
        ("07", "An SPARC (141 Galaxien, Gas aus V_gas", float(_sparc_bahn_feld()["galaxien"]), 141),
        ("07", "liegt die Feldform als Bahnregel Δχ² = 368 über ν·g_N", _sparc_bahn_feld()["chi_w"] - _sparc_bahn_feld()["chi_v"], 368),
        ("07", "mit der Steigung 0,05 ± 0,12, wo sie 1 verlangt", _sparc_bahn_feld()["steigung"], 0.05),
        ("07", "mit der Steigung 0,05 ± 0,12, wo sie 1 verlangt", _sparc_bahn_feld()["streuung"], 0.12),
        ("03", "innerhalb 2R_d im Median um 0,05 dex", -_sparc_bahn_feld()["innen"], 0.05),
        ("03", "(flächenschwach 0,08, flächenhell 0,02", -_sparc_bahn_feld()["lsb"], 0.08),
        ("03", "(flächenschwach 0,08, flächenhell 0,02", -_sparc_bahn_feld()["hsb"], 0.02),
        ("03", "außerhalb 4R_d und bei Kugeln sind beide gleich", round(100 * abs(_sparc_bahn_feld()["aussen"])), 0),
        ("03", "über einer Kuzmin-Scheibe genau den Wert, den das Licht in der Ebene liest", round(100 * _achsen_zweiter_weg()[3]), 0),
        ("03", "Innerhalb 2R_d der SPARC-Scheiben sind das im Median 0,07 dex der Beschleunigung", -_sparc_bahn_feld()["stufe_innen"], 0.07),
        ("03", "Innerhalb 2R_d der SPARC-Scheiben sind das im Median 0,07 dex der Beschleunigung, 8 % des Tempos", 100 * (1 - 10**(_sparc_bahn_feld()["stufe_innen"] / 2)), 8),
        ("03", "(flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01", -_sparc_bahn_feld()["stufe_lsb"], 0.11),
        ("03", "(flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01", -_sparc_bahn_feld()["stufe_hsb"], 0.02),
        ("03", "(flächenschwach 0,11, flächenhell 0,02), außerhalb 4R_d 0,01", -_sparc_bahn_feld()["stufe_aussen"], 0.01),
        ("F:galaxien", "Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s", _haufen_drehung()["v_gas"][0], 63),
        ("F:galaxien", "Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s", _haufen_drehung()["v_gas"][1], 132),
        ("F:galaxien", "Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s", _haufen_drehung()["v_mode"][0], 29),
        ("F:galaxien", "Es dreht mit j/R₂₀₀ = √2·λ'·V₂₀₀ = 63 bis 132 km/s, die gemischte Mode mit 29 bis 63 km/s", _haufen_drehung()["v_mode"][1], 63),
        ("F:galaxien", "In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124", _galaxie_drehung()["M200"][0], 3.8),
        ("F:galaxien", "In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124", _galaxie_drehung()["M200"][1], 5.4),
        ("F:galaxien", "In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124", _galaxie_drehung()["lam_lin"][0], 0.016),
        ("F:galaxien", "In Läufen vieler Körper auf der Skala der Milchstraße (Gebiet (1 + Q)·M_b, M₂₀₀ = 3,8 bis 5,4·10¹¹ M☉) gibt der Kommutator bis zur Umkehr λ' = 0,016 bis 0,124", _galaxie_drehung()["lam_lin"][1], 0.124),
        ("F:galaxien", "am Ende trägt die gemischte Mode nur 0,013 bis 0,043", _galaxie_drehung()["lam"][0], 0.013),
        ("F:galaxien", "am Ende trägt die gemischte Mode nur 0,013 bis 0,043", _galaxie_drehung()["lam"][1], 0.043),
        ("F:galaxien", "an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135)", _galaxie_drehung()["sparc"][1], 0.073),
        ("F:galaxien", "an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135)", _galaxie_drehung()["sparc"][2], 0.047),
        ("F:galaxien", "an SPARC im Median 0,073 (16 bis 84 %: 0,047 bis 0,135)", _galaxie_drehung()["sparc"][3], 0.135),
        ("F:galaxien", "aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157)", D.sparc_drehung(D.sparc_eichung(K.kpc), _cmb()["saett"]["punkt"]["H0"], K.q_je_radiant(), True)[1], 0.098),
        ("F:galaxien", "aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157)", D.sparc_drehung(D.sparc_eichung(K.kpc), _cmb()["saett"]["punkt"]["H0"], K.q_je_radiant(), True)[2], 0.057),
        ("F:galaxien", "aus den Profilen von Sternen und Gas mit der gemessenen Kurve 0,098 (0,057 bis 0,157)", D.sparc_drehung(D.sparc_eichung(K.kpc), _cmb()["saett"]["punkt"]["H0"], K.q_je_radiant(), True)[3], 0.157),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["M200"][0], 1.0),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["M200"][1], 3.5),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["cos"][0], 0.95),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["cos"][1], 0.99),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["verhaeltnis"][0], 0.78),
        ("F:galaxien", "An den Anfängen der Läufe vieler Körper (M₂₀₀ = 1,0 bis 3,5·10¹⁵ M☉) trifft der Kommutator die Richtung des ganzen linearen Drehimpulses auf cos = 0,95 bis 0,99, den Betrag auf 0,78 bis 1,20", _haufen_drehung()["verhaeltnis"][1], 1.2),
        ("F:galaxien", "bis zur Umkehr gerechnet ist der Anteil der Drehung λ' = J/(√2·M₂₀₀·V₂₀₀·R₂₀₀) = 0,022 bis 0,053", _haufen_drehung()["lam_lin"][0], 0.022),
        ("F:galaxien", "bis zur Umkehr gerechnet ist der Anteil der Drehung λ' = J/(√2·M₂₀₀·V₂₀₀·R₂₀₀) = 0,022 bis 0,053", _haufen_drehung()["lam_lin"][1], 0.053),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["lam"][0], 0.013),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["lam"][1], 0.03),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["j"][0], 2.7),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["j"][1], 9.4),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["phi"][0], 68),
        ("F:galaxien", "Am Ende der Läufe trägt die Drehung in R₂₀₀ λ' = 0,013 bis 0,03, die Spinzahl ist j = √2·λ'·c/V₂₀₀ = 2,7 bis 9,4, und die Horizonte stehen 68° bis 84° zur Tiefe", _haufen_drehung()["phi"][1], 84),
        ("03", "In einer Scheibe an der Schwelle liegt Gas an der Grenze Q_T = 1 innerhalb 4,2·R_d über der Gezeitendichte", D.sparc_gas_gezeit()[0], 4.2),
        ("03", "In SPARC liegen dort vier Zehntel des Gases, aus V_gas im Median 0,37, nach Masse 0,40", D.sparc_gas_gezeit()[2], 0.37),
        ("03", "In SPARC liegen dort vier Zehntel des Gases, aus V_gas im Median 0,37, nach Masse 0,40", D.sparc_gas_gezeit()[3], 0.40),
        ("03", "unabhängig aus M_HI und R_HI 0,36 und 0,40", D.sparc_gas_gezeit(1.0, "hi")[2], 0.36),
        ("03", "unabhängig aus M_HI und R_HI 0,36 und 0,40", D.sparc_gas_gezeit(1.0, "hi")[3], 0.40),
        ("03", "Sechs Zehntel verlangen Gas wenig über der Grenze, Q_T = 1,2 bis 1,4", _gas_stabil("vgas"), 1.2),
        ("03", "Sechs Zehntel verlangen Gas wenig über der Grenze, Q_T = 1,2 bis 1,4", _gas_stabil("hi"), 1.4),
        ("03", "Q_T = 1,2 bis 1,4, aus der Flächendichte 1,28", _gas_stabil("flaeche"), 1.28),
        ("03", "aus der Flächendichte, die V_gas als dünne Scheibe verlangt, ebenso 0,37 und 0,40", D.sparc_gas_gezeit(1.0, "flaeche")[2], 0.37),
        ("03", "aus der Flächendichte, die V_gas als dünne Scheibe verlangt, ebenso 0,37 und 0,40", D.sparc_gas_gezeit(1.0, "flaeche")[3], 0.40),
        ("03", "Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ", (lambda w: (w[1] - w[0]) / 0.1)(K.hi_linie(0.1)), 1.5),
        ("03", "Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ", (lambda w: (w[1] - w[0]) / 0.1)(K.hi_linie(0.1, h_gas=3.0)), 1.5),
        ("03", "Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ", 2 * __import__("scipy.stats", fromlist=["norm"]).norm.ppf(0.8), 1.68),
        ("03", "Eine Scheibe an der Schwelle gibt W20 − W50 = 1,5·σ, der Rand einer Stufe 1,68·σ", 2 * __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: 0.5 * math.erfc(x / math.sqrt(2)) - 0.2, 0.0, 3.0), 1.68),
        ("03", "im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s", min(D.hi_rand_gang("alfalfa")[0]), 19),
        ("03", "im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s", max(D.hi_rand_gang("alfalfa")[0]), 28),
        ("03", "im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s", min(D.hi_rand_gang("hipass")[0]), 26),
        ("03", "im Median bei 19 bis 28 km/s, an HIPASS bei 26 bis 30 km/s", max(D.hi_rand_gang("hipass")[0]), 30),
        ("03", "Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s", _gas_streuung()["leicht_rand"], 11),
        ("03", "Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s", _gas_streuung()["leicht"], 12),
        ("03", "Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s", _gas_streuung()["schwer_rand"], 19),
        ("03", "Bis lg M_HI = 9,5 streut das Gas mit 11 bis 12 km/s, das schwerste mit 19 bis 22 km/s", _gas_streuung()["schwer"], 22),
        ("03", "aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5", _gas_streuung()["flaeche"][0], 0.6),
        ("03", "aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5", _gas_streuung()["flaeche_rand"][0], 0.6),
        ("03", "aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5", _gas_streuung()["flaeche"][1], 0.5),
        ("03", "aus der Flächendichte im Median 0,6 des Gases unter der Gezeitendichte, nach Masse 0,5", _gas_streuung()["flaeche_rand"][1], 0.5),
        ("03", "aus V_gas und aus M_HI und R_HI 0,5 bis 0,6", min(_gas_streuung()["andere"]), 0.5),
        ("03", "aus V_gas und aus M_HI und R_HI 0,5 bis 0,6", max(_gas_streuung()["andere"]), 0.6),
        ("03", "zwischen einem Drittel und dem Dreifachen der gelesenen sind es 0,3 bis 0,8", min(_gas_streuung()["spanne"]), 0.3),
        ("03", "zwischen einem Drittel und dem Dreifachen der gelesenen sind es 0,3 bis 0,8", max(_gas_streuung()["spanne"]), 0.8),
        ("03", "wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen", min(_gas_raum()[0]), 0.5),
        ("03", "wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen", max(_gas_raum()[0]), 0.6),
        ("03", "wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen", min(_gas_raum()[1]), 0.6),
        ("03", "wie das Ganze es zählt, sind es ebenso 0,5 bis 0,6, also F = 0,6 bis 0,8 % der Baryonen", max(_gas_raum()[1]), 0.8),
        ("03", "Mit diesem F₀ gehalten lesen alle Präzisionsdaten χ² = 2025,2 bis 2026,3", min(_formteil_gas()[0]), 2025.2),
        ("03", "Mit diesem F₀ gehalten lesen alle Präzisionsdaten χ² = 2025,2 bis 2026,3", max(_formteil_gas()[0]), 2026.3),
        ("03", "Δχ² = 1,5 bis 2,6 unter ΛCDM, mit einer freien Zahl weniger", _formteil_gas()[1] - max(_formteil_gas()[0]), 1.5),
        ("03", "Δχ² = 1,5 bis 2,6 unter ΛCDM, mit einer freien Zahl weniger", _formteil_gas()[1] - min(_formteil_gas()[0]), 2.6),
        ("03", "Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003)", min(_formteil_lesen()[0]), 1.058),
        ("03", "Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003)", max(_formteil_lesen()[0]), 1.064),
        ("03", "Frei gelassen liest dann auch der Formteil a₀ = 1,058 bis 1,064 (± 0,003)", max(_formteil_lesen()[1]), 0.003),
        ("F:galaxien", "An den 36 Satelliten der Milchstraße mit Eigenbewegung", float(_gezeit_satelliten()["zahl"]), 36),
        ("F:galaxien", "Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9)", min(_gezeit_satelliten()["eng"]), 0.4),
        ("F:galaxien", "Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9)", max(_gezeit_satelliten()["eng"]), 1.3),
        ("F:galaxien", "Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9)", min(_gezeit_satelliten(True)["eng"]), 0.7),
        ("F:galaxien", "Antlia II und Crater II mit 0,4 bis 1,3·R_e (mit dem Gas, das sie einmal trugen, 0,7 bis 1,9)", max(_gezeit_satelliten(True)["eng"]), 1.9),
        ("F:galaxien", "bei allen übrigen liegt er über 2,6·R_e", _gezeit_satelliten()["naechster"], 2.6),
        ("F:galaxien", "im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit", float(_gezeit_satelliten()["raenge"][0]), 1),
        ("F:galaxien", "im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit", float(_gezeit_satelliten()["raenge"][1]), 2),
        ("F:galaxien", "im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit", float(_gezeit_satelliten()["raenge"][2]), 6),
        ("F:galaxien", "im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit", 100 * _gezeit_satelliten()["zufall"], 0.3),
        ("F:galaxien", "im Rest gegen die Regel auf Rang 1, 2 und 6, zufällig mit 0,3 % Wahrscheinlichkeit", 100 * _gezeit_satelliten(weg="bahn")["zufall"], 0.3),
        ("F:galaxien", "GC5 den Kern in 6,7 bis 21 Gyr", _fornax("summe")["gc5"][0], 6.7),
        ("F:galaxien", "GC5 den Kern in 6,7 bis 21 Gyr", _fornax("summe")["gc5"][1], 21),
        ("F:galaxien", "GC1 in 27 bis 84 Gyr", _fornax("summe")["gc1"][0], 27),
        ("F:galaxien", "GC1 in 27 bis 84 Gyr", _fornax("summe")["gc1"][1], 84),
        ("F:galaxien", "in 1,2 bis 4,8 Gyr", _fornax("summe")["gc2"][0], 1.2),
        ("F:galaxien", "in 1,2 bis 4,8 Gyr", _fornax("summe")["gc2"][1], 4.8),
        ("F:galaxien", "GC5 den Kern in 6,7 bis 21 Gyr", _fornax("bahn")["gc5"][0], 6.7),
        ("F:galaxien", "GC5 den Kern in 6,7 bis 21 Gyr", _fornax("bahn")["gc5"][1], 21),
        ("F:galaxien", "GC1 in 27 bis 84 Gyr", _fornax("bahn")["gc1"][0], 27),
        ("F:galaxien", "GC1 in 27 bis 84 Gyr", _fornax("bahn")["gc1"][1], 84),
        ("F:galaxien", "in 1,2 bis 4,8 Gyr", _fornax("bahn")["gc2"][0], 1.2),
        ("F:galaxien", "in 1,2 bis 4,8 Gyr", _fornax("bahn")["gc2"][1], 4.8),
        ("09", "Ein Stern der Form n = 3 (k = 0,075) erreicht höchstens j = 52 (Sonne)", _drehung_grenze()[0], 0.075),
        ("09", "Ein Stern der Form n = 3 (k = 0,075) erreicht höchstens j = 52 (Sonne)", _drehung_grenze()[1], 52),
        ("09", "der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1", _drehung_grenze()[2], 1.417),
        ("09", "der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1", _drehung_grenze()[3] / 1e3, 1.03),
        ("09", "der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1", _drehung_grenze()[4], 0.095),
        ("09", "der Weiße Zwerg an seiner Grenzmasse (1,417 M☉, R = 1,03·10³ km, k = 0,095) 2,1", _drehung_grenze()[5], 2.1),
        ("09", "der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51", _drehung_grenze()[6], 0.710),
        ("09", "der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51", _drehung_grenze()[7], 9.16),
        ("09", "der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51", _drehung_grenze()[8], 0.174),
        ("09", "der Neutronenstern an der Grenze seines Halters (0,710 M☉, R = 9,16 km, k = 0,174) nur 0,51", _drehung_grenze()[9], 0.51),
        ("09", "1 + α/(2π) = 1,001161", K.erste_rueckkehr_windung(), 1.001161),
        ("09", "Gemessen ist 1,001160", K.g_e_CODATA / 2, 1.001160),
        ("F:ganzes", "+(α/2π)² = +1,35·10⁻⁶", (K.rueckkehr_windung_reihe(2) - K.rueckkehr_windung_reihe(1)) / 1e-6, 1.35),
        ("F:ganzes", "gezählt +0,5, −0,25, +0,125, −0,0625, +0,031", sum(abs(K.rueckkehr_windung_wechselnd(n) - w) < 5e-4 for n, w in zip(range(1, 6), (0.5, -0.25, 0.125, -0.0625, 0.031))), 5),
        ("F:ganzes", "+0,5, −0,328, +1,181, −1,91, +9,2", sum(abs(a - w) < 0.06 for a, w in zip(K.QED_BEIWERTE, (0.5, -0.328, 1.181, -1.91, 9.2))), 5),
        ("F:ganzes", "alle fünf Vorzeichen wie die halbe Drehung", sum(math.copysign(1, a) == math.copysign(1, K.rueckkehr_windung_wechselnd(n)) for n, a in zip(range(1, 6), K.QED_BEIWERTE)), 5),
        ("F:ganzes", "Rückkehr −1,76·10⁻⁶", (K.g_e_CODATA / 2 - K.rueckkehr_windung_reihe(1)) / 1e-6, -1.76),
        ("F:ganzes", "also −0,33·(α/π)²", (K.g_e_CODATA / 2 - K.rueckkehr_windung_reihe(1)) / (K.alpha / math.pi)**2, -0.33),
        ("09", "B = M²c²/(eħ) = 4,41·10⁹ T", K.feld_treffpunkt(K.m_e) / 1e9, 4.41),
        ("09", "E = B·c = 1,32·10¹⁸ V/m", K.feld_treffpunkt(K.m_e) * K.c / 1e18, 1.32),
        ("09", "= −2,99·10⁻¹¹ m³/mol", K.diamagnetismus_h() / 1e-11, -2.99),
        ("09", "das sind −2,38·10⁻⁶ cm³/mol", K.diamagnetismus_h() * 1e6 / (4 * math.pi) / 1e-6, -2.38),
        ("09", "ν(s) = ½ + 1/s + ½·L(s/2)", float(max(abs(K.verstaerkung(K.a0 * x * x) - (0.5 + 1 / x + 0.5 * K.langevin(x / 2))) for x in (0.2, 1.0, 3.0, 8.0)) < 1e-9), 1.0),
        ("09", "m = tanh(μB/kT) = 1 − 2f(2μB/kT)", float(abs(K.ausrichtung_windung(0.7) - (1 - 2 / (math.exp(1.4) + 1))) < 1e-12), 1.0),
        ("09", "m ≈ √(3(1 − T/T_C))", K.einrasten_magnet(0.999) / math.sqrt(3 * 0.001), 1.0),
        ("09", "h_c ≈ (2/3)(1 − T/T_C)^{3/2}", K.hysteresefenster(0.999) / (2 / 3 * 0.001**1.5), 1.0),
        ("09", "eingerastet (ohne Spannung) für |I| ≤ I_c", float(K.josephson_einrasten(0.5e-3, 1e-3, 1.0)[2] and not K.josephson_einrasten(2e-3, 1e-3, 1.0)[2]), 1.0),
        ("09", "Φ₀ = h/(2e) = 2,068·10⁻¹⁵ Wb", K.flussquant() / 1e-15, 2.068),
        ("09", "sind das 17 nm", K.london_tiefe(1e29) * 1e9, 17),
        ("05", "M = m_p·(α/α_G)^{3/2} = 2,30·10²⁷ kg", K.masse_windung_faltung()[0] / 1e27, 2.30),
        ("05", "Tiefe dort δ = α²·m_e/m_p = 2,9·10⁻⁸", K.masse_windung_faltung()[1] / 1e-8, 2.9),
        ("05", "(Erde δ = 7,0·10⁻¹⁰)", K.trichter(K.M_earth, K.R_ERDE)["tiefe"] / 1e-10, 7.0),
        ("05", "(Sonne δ = 2,1·10⁻⁶)", K.trichter(K.M_sun, K.R_SONNE)["tiefe"] / 1e-6, 2.1),
        ("05", "Verhältnis α^{3/2} = 6,2·10⁻⁴", K.masse_windung_faltung()[0] / K.sternmasse() / 1e-4, 6.2),
        ("05", "a₀/(πG) = 242 M☉/pc²", K.flaechendichte_schwelle() * K.pc**2 / K.M_sun, 242),
        ("05", "Wasserstoff: (8π/3)·α³ = 3,3·10⁻⁶", K.verlust_windung(K.bohr_radius(), K.alpha * K.c) / 1e-6, 3.3),
        ("05", "Erde um die Sonne: 2,3·10⁻²⁴", K.verlust_faltung(K.M_earth, K.M_sun, K.AU) / 1e-24, 2.3),
        ("05", "in 1,6·10⁻¹¹ s auf den Kern", K.fall_offene_windung() / 1e-11, 1.6),
        ("05", "trägt die Windung 7/8 der Faltung", K.strahlung_je_sinn(2) / K.strahlung_je_sinn(1), 0.875),
        ("F:ganzes", "ρ_DE = ρ_s·(3 + 16,8·F)", K.zustandsteil_je_gas(K.q_je_radiant()), 16.8),
        ("F:ganzes", "das Maß 17,2 inaktiv", K.zustandsteil_je_gas(), 17.2),
        ("F:lesungen", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[1], 0.01),
        ("F:lesungen", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[2], 0.02),
        ("F:lesungen", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[3], 21),
        ("F:lesungen", "Über alle 15 Bins", float(D.kids_typ_reste("blau")[0]), 15),
        ("F:lesungen", "die späten +0,02 ± 0,02", D.kids_typ_reste("spaet")[1], 0.02),
        ("F:lesungen", "die späten +0,02 ± 0,02", D.kids_typ_reste("spaet")[2], 0.02),
        ("F:lesungen", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[1], 0.23),
        ("F:lesungen", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[2], 0.01),
        ("F:lesungen", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[3], 387),
        ("F:lesungen", "die frühen +0,21 ± 0,01", D.kids_typ_reste("frueh")[1], 0.21),
        ("F:lesungen", "die frühen +0,21 ± 0,01", D.kids_typ_reste("frueh")[2], 0.01),
        ("F:lesungen", "liegen die roten +0,13 ± 0,02", D.kids_typ_reste("rot", -12.5)[1], 0.13),
        ("F:lesungen", "liegen die roten +0,13 ± 0,02", D.kids_typ_reste("rot", -12.5)[2], 0.02),
        ("F:lesungen", "die frühen +0,14 ± 0,02 über der Regel", D.kids_typ_reste("frueh", -12.5)[1], 0.14),
        ("F:lesungen", "die frühen +0,14 ± 0,02 über der Regel", D.kids_typ_reste("frueh", -12.5)[2], 0.02),
        ("F:lesungen", "wächst es auf +0,29 und +0,26", D.kids_typ_reste("rot", -14.5, -13.0)[1], 0.29),
        ("F:lesungen", "wächst es auf +0,29 und +0,26", D.kids_typ_reste("frueh", -14.5, -13.0)[1], 0.26),
        ("F:lesungen", "am Gasanteil (Rangkorrelation −0,04)", D.sparc_eigenschaften(K.kpc, K.a0)["Gasanteil"][0], -0.04),
        ("F:lesungen", "noch am Typ (−0,10)", D.sparc_eigenschaften(K.kpc, K.a0)["Typ"][0], -0.10),
        ("F:lesungen", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_lcdm("blau")[1], -0.31),
        ("F:lesungen", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_lcdm("blau")[2], 169),
        ("F:lesungen", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_reste("blau")[3], 21),
        ("F:lesungen", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_lcdm("spaet")[1], -0.26),
        ("F:lesungen", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_lcdm("spaet")[2], 157),
        ("F:lesungen", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_reste("spaet")[3], 14),
        ("F:lesungen", "trifft die roten im Mittel (χ² 146 gegen 387)", D.kids_typ_lcdm("rot")[2], 146),
        ("F:lesungen", "trifft die roten im Mittel (χ² 146 gegen 387)", D.kids_typ_reste("rot")[3], 387),
        ("F:lesungen", "und die frühen (85 gegen 268)", D.kids_typ_lcdm("frueh")[2], 85),
        ("F:lesungen", "und die frühen (85 gegen 268)", D.kids_typ_reste("frueh")[3], 268),
        ("F:lesungen", "mit χ² 34 bei 20 Werten", D.mistele_reste()[0], 34),
        ("03", "bei R₅₀₀ 3,4 bis 7,0 keV", D.xcop_abkuehlung()[2], 3.4),
        ("03", "bei R₅₀₀ 3,4 bis 7,0 keV", D.xcop_abkuehlung()[3], 7.0),
        ("F:ganzes", "im Median nur 0,4 % des Gases innerhalb R₅₀₀", D.xcop_abkuehlung()[0] * 100, 0.4),
        ("03", "höchstens 3,7 % des Gases abkühlbar", D.xcop_abkuehlung()[1] * 100, 3.7),
        ("03", "bei 0,1·R₅₀₀ 4,8 bis 9,6 keV", D.xcop_abkuehlung(ort=0.1)[2], 4.8),
        ("03", "bei 0,1·R₅₀₀ 4,8 bis 9,6 keV", D.xcop_abkuehlung(ort=0.1)[3], 9.6),
        ("03", "Abkühlzeit von 13,6 bis 38 Gyr", D.xcop_abkuehlung(ort=0.1)[4], 13.6),
        ("03", "Abkühlzeit von 13,6 bis 38 Gyr", D.xcop_abkuehlung(ort=0.1)[5], 38),
        ("F:ganzes", "in 13,8 Gyr im Median nur", K.weltalter(lambda a: K.dichte_dunkle_energie()) / K.Gyr, 13.8),
        ("03", "in 13,8 Gyr höchstens 3,7 %", K.weltalter(lambda a: K.dichte_dunkle_energie()) / K.Gyr, 13.8),
        ("F:ganzes", "liegt die Dichte des neutralen Wasserstoffs in sternbildenden Galaxien bei (6,5 ± 1,0 ± 1,2)·10⁷ M☉/Mpc³", D.hi_frueh()[1] / 1e7, 6.5),
        ("F:ganzes", "das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14)", _hi_gang()[0], 1.10),
        ("F:ganzes", "das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14)", _hi_gang()[1], 0.27),
        ("F:ganzes", "das 1,10 ± 0,27-fache von heute (ALFALFA mit demselben H₀, HIPASS 1,14)", D.hi_frueh()[1] / D.hipass_hi_dichte(H0=70.0)[0], 1.14),
        ("F:ganzes", "Der Zustandsteil verlangt am Maximum das 1,74-fache des heutigen gebundenen Gases (Sättigung beim Standardwert) oder das 2,2-fache (eigene Form des Modells)", _hi_gang()[2], 1.74),
        ("F:ganzes", "Der Zustandsteil verlangt am Maximum das 1,74-fache des heutigen gebundenen Gases (Sättigung beim Standardwert) oder das 2,2-fache (eigene Form des Modells)", _hi_gang()[4], 2.2),
        ("F:ganzes", "das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien", D.hi_frueh_eins()[1], 3.5),
        ("F:ganzes", "ihr Verhältnis ist r_e/(3Δs_p) = 2,269·10³⁹", K.atom_und_schliessung()[1] / 1e39, 2.269),
        ("F:ganzes", "ihr Verhältnis ist r_e/(3Δs_p) = 2,269·10³⁹", K.alpha * K.hbar * K.c / (K.G * K.m_e * K.m_p) / 1e39, 2.269),
        ("F:ganzes", "Das Elektron liegt von der Schließung ℓ₀/λ̄_e = 2,196·10³⁹ entfernt, 3 % weniger", K.atom_und_schliessung()[0] / 1e39, 2.196),
        ("F:ganzes", "Das Elektron liegt von der Schließung ℓ₀/λ̄_e = 2,196·10³⁹ entfernt, 3 % weniger", 100 * (1 - K.atom_und_schliessung()[0] / K.atom_und_schliessung()[1]), 3),
        ("F:ganzes", "Gemessen ist r*_p/√(λ̄_e·r_e) = 0,984", K.atom_und_schliessung()[2], 0.984),
        ("F:ganzes", "a₀ = G·m_p/(α³·a_B²) = 1,026·10⁻¹⁰ m/s²", K.atom_und_schliessung()[3] / 1e-10, 1.026),
        ("F:ganzes", "a₀ = G·m_p/(α³·a_B²) = 1,026·10⁻¹⁰ m/s²", K.G * K.m_e**2 * K.m_p * K.c**2 / (K.alpha * K.hbar**2) / 1e-10, 1.026),
        ("F:ganzes", "N_b = (3/(2(1 + Q)))·α·(m_Pl²/(m_e·m_p))² = 1,68·10⁸⁰", K.atom_und_schliessung()[4] / 1e80, 1.68),
        ("F:ganzes", "Ist der Formteil 1,026, lesen alle Präzisionsdaten χ² = 2020,5 (mit 1,06: 2019,7)", _form_atom(1.026)[0], 2020.5),
        ("F:ganzes", "Ist der Formteil 1,026, lesen alle Präzisionsdaten χ² = 2020,5 (mit 1,06: 2019,7)", _form_atom(1.06)[0], 2019.7),
        ("F:ganzes", "der Zustandsteil trägt heute 13,3 ± 0,9 % der Dunklen Energie", 100 * _form_atom(1.026)[2], 13.3),
        ("F:ganzes", "der Zustandsteil trägt heute 13,3 ± 0,9 % der Dunklen Energie", 100 * _form_atom(1.026)[3], 0.9),
        ("F:ganzes", "Das verlangt F = 2,37 ± 0,17 % der Baryonen als dünnes gebundenes Gas", 100 * _form_atom(1.026)[4], 2.37),
        ("F:ganzes", "Das verlangt F = 2,37 ± 0,17 % der Baryonen als dünnes gebundenes Gas", 100 * _form_atom(1.026)[5], 0.17),
        ("F:ganzes", "verlangt diese Grenze vom Formteil a₀ ≥ 1,054, dort ist F = 1,21 ± 0,13 %", 100 * _form_atom(1.054)[4], 1.21),
        ("F:ganzes", "verlangt diese Grenze vom Formteil a₀ ≥ 1,054, dort ist F = 1,21 ± 0,13 %", 100 * _form_atom(1.054)[5], 0.13),
        ("F:ganzes", "das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien", D.hi_frueh_eins()[2], 0.5),
        ("F:ganzes", "das 3,5 ± 0,5-fache des neutralen Gases von heute (GMRT-CATz1, 11 419 Galaxien", float(D.hi_frueh_eins()[3]), 11419),
        ("F:ganzes", "der Zustandsteil verlangt in der eigenen Form des Modells dort das 0,34-fache", _gas_form()["F_eins"], 0.34),
        ("F:ganzes", "bei z = 1,04 ist ihre Flächendichte das 1,48-fache der Schwelle", D.hi_frueh_eins()[0], 1.04),
        ("F:ganzes", "bei z = 1,04 ist ihre Flächendichte das 1,48-fache der Schwelle", D.gas_duenn_frueh(_saett_punkt())[0], 1.48),
        ("F:ganzes", "dann 3 % weniger des Gases unter der Gezeitendichte", round(100 * (1 - D.gas_duenn_frueh(_saett_punkt())[3])), 3),
        ("F:ganzes", "dann 3 % weniger des Gases unter der Gezeitendichte", round(100 * (1 - D.gas_duenn_frueh(_saett_punkt(), weg="vgas")[3])), 3),
        ("F:ganzes", "dann 3 % weniger des Gases unter der Gezeitendichte", round(100 * (1 - D.gas_duenn_frueh(_saett_punkt(), weg="hi")[3])), 3),
        ("F:ganzes", "bei z = 0,34 1 % weniger", round(100 * (1 - D.gas_duenn_frueh(_saett_punkt(), z=0.34)[3])), 1),
        ("F:ganzes", "bei der dreifachen auf das 1,4-fache", D.gas_duenn_frueh(_saett_punkt(), faktor=3.0)[3], 1.4),
        ("F:ganzes", "mit dem gezählten F₀ = 0,67 % beim Standardwert", 100 * _gas_frueh()["F0_gezaehlt"], 0.67),
        ("F:ganzes", "nur auf das 1,58-fache, 3,7 Fehler unter der Messung", _gas_frueh()["fest"]["eins"], 1.58),
        ("F:ganzes", "nur auf das 1,58-fache, 3,7 Fehler unter der Messung", -_gas_frueh()["fest"]["eins_rest"], 3.7),
        ("F:ganzes", "Δχ² = 4,6 über dem Gang mit dem Maximum an der Sättigung, 2,0 über ΛCDM", _gas_frueh()["fest"]["chi2"] - _gas_frueh()["saett"]["chi2"], 4.6),
        ("F:ganzes", "Δχ² = 4,6 über dem Gang mit dem Maximum an der Sättigung, 2,0 über ΛCDM", _gas_frueh()["fest"]["chi2"] - _gas_frueh()["lcdm"], 2.0),
        ("F:ganzes", "gäbe bei z = 1,04 das 0,67-fache, 5,6 Fehler darunter", _gas_frueh()["saett"]["eins"], 0.67),
        ("F:ganzes", "gäbe bei z = 1,04 das 0,67-fache, 5,6 Fehler darunter", -_gas_frueh()["saett"]["eins_rest"], 5.6),
        ("F:ganzes", "(das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber)", _gas_frueh()["frei"]["eins"], 2.67),
        ("F:ganzes", "(das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber)", -_gas_frueh()["frei"]["eins_rest"], 1.6),
        ("F:ganzes", "(das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber)", _gas_frueh()["frei"]["bera"], 1.70),
        ("F:ganzes", "(das 2,67-fache bei z = 1,04, 1,6 Fehler darunter, das 1,70-fache bei z = 0,34, 2,2 Fehler darüber)", _gas_frueh()["frei"]["bera_rest"], 2.2),
        ("F:ganzes", "nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter", 100 * _gas_frueh()["frei"]["F0"], 0.38),
        ("F:ganzes", "nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter", 100 * _gas_frueh()["frei"]["fehler"], 0.07),
        ("F:ganzes", "nur mit F₀ = 0,38 ± 0,07 %, vier Fehler unter", round((_gas_frueh()["F0_gezaehlt"] - _gas_frueh()["frei"]["F0"]) / _gas_frueh()["frei"]["fehler"]), 4),
        ("F:ganzes", "und liegt dann 0,3 über ΛCDM", _gas_frueh()["frei"]["chi2"] - _gas_frueh()["lcdm"], 0.3),
        ("F:ganzes", "= −0,328, ein rationaler Teil 1,368", K.zweite_rueckkehr_windung()[0], -0.328),
        ("F:ganzes", "= −0,328, ein rationaler Teil 1,368", K.zweite_rueckkehr_windung(numerisch=True)[0], -0.328),
        ("F:ganzes", "= −0,328, ein rationaler Teil 1,368", K.zweite_rueckkehr_windung()[1], 1.368),
        ("F:ganzes", "ζ(n) = Σ1/kⁿ, mit +1,724", K.zweite_rueckkehr_windung()[2], 1.724),
        ("F:ganzes", "ζ(n) = Σ1/kⁿ, mit +1,724", K.zweite_rueckkehr_windung(numerisch=True)[2], 1.724),
        ("F:ganzes", "ln 2 = Σ(−1)^{k+1}/k, mit −3,421", K.zweite_rueckkehr_windung()[3], -3.421),
        ("F:ganzes", "ln 2 = Σ(−1)^{k+1}/k, mit −3,421", K.zweite_rueckkehr_windung(numerisch=True)[3], -3.421),
        ("F:ganzes", "Δ₀/kT_c = π·e^{−γ_E} = 1,764", K.paar_einrasten(0.15)[2], 1.764),
        ("F:ganzes", "Δ₀/kT_c = π·e^{−γ_E} = 1,764", K.paar_einrasten(0.15, numerisch=True)[2], 1.764),
        ("F:ganzes", "das Maximum liegt jenseits z = 1", float(min(_gas_frueh()[n]["z_max"] for n in ("fest", "frei", "form_frei")) > 1), 1.0),
        ("F:ganzes", "χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein", _gas_abstaende()["frueh"], 1414.1),
        ("F:ganzes", "χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein", _gas_abstaende()["saett"], 1416.9),
        ("F:ganzes", "χ² = 1414,1 gegen 1416,9 mit dem Maximum an der Sättigung und 1432,6 mit dem Formteil allein", _gas_abstaende()["form"], 1432.6),
        ("F:ganzes", "den Preis zahlen die Spektren der Hintergrundstrahlung, Δχ² = 7,2", _gas_frueh()["spektren"], 7.2),
        ("F:ganzes", "Ist der Formteil frei, liest er a₀ = 1,046", _gas_frueh()["form_frei"]["form"], 1.046),
        ("F:ganzes", "das Gas folgt mit 1,9 Fehlern an beiden Punkten", _gas_frueh()["form_frei"]["bera_rest"], 1.9),
        ("F:ganzes", "das Gas folgt mit 1,9 Fehlern an beiden Punkten", -_gas_frueh()["form_frei"]["eins_rest"], 1.9),
        ("F:ganzes", "und die Präzisionsdaten liegen 1,2 über ΛCDM", _gas_frueh()["form_frei"]["chi2"] - _gas_frueh()["lcdm"], 1.2),
        ("F:ganzes", "Dann liegen die Abstände mit dem Gas 2,5 höher als ungekoppelt", _gas_abstaende()["gekoppelt"] - _gas_abstaende()["frueh"], 2.5),
        ("F:ganzes", "Die Messung liegt 1,3 und 0,9 Fehler darunter", _hi_gang()[3], 1.3),
        ("F:ganzes", "Die Messung liegt 1,3 und 0,9 Fehler darunter", _hi_gang()[5], 0.9),
        ("07", "Mit D = GM/c ist die Tiefe der Rückkehr ϖ_R = √(D·τ_D) = √(GM/a₀) = r*", float(max(abs(math.sqrt(K.G * m / K.c * K.rueckkehrzeit()) / K.schwellenradius(m) - 1) for m in (K.m_p, K.M_sun, 1e11 * K.M_sun)) < 1e-12), 1),
        ("F:ganzes", "bei 30 kAU (1,5 M☉) ν = 3,8", K.verstaerkung(K.G * 1.5 * K.M_sun / (30e3 * K.AU)**2), 3.8),
        ("F:ganzes", "die 1,95-fache Geschwindigkeit", K.verstaerkung(K.G * 1.5 * K.M_sun / (30e3 * K.AU)**2)**0.5, 1.95),
        ("F:ganzes", "(Gaia DR3, bis 130 pc, 1 bis 30 kAU)", float(K.WEITE_DOPPELSTERNE["abstand_pc"]), 130),
        ("F:ganzes", "der Zuwachs von 20 % nach MOND fehlt", K.WEITE_DOPPELSTERNE["mond_zuwachs"] * 100, 20),
        ("03", "die Glättungslänge ℓ = 349 kAU", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3, 349),
        ("03", "eines Paares von 1 M☉ (277 kAU)", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3 / 2**(1 / 3), 277),
        ("F:lesungen", "mit χ² 34 bei 20 Werten", float(D.mistele_reste()[1]), 20),
        ("F:lesungen", "im Mittel je Bin −0,03 bis +0,03 dex", min(D.mistele_reste()[2]), -0.03),
        ("F:lesungen", "im Mittel je Bin −0,03 bis +0,03 dex", max(D.mistele_reste()[2]), 0.03),
        ("F:ganzes", "verlangt F = 1,8 %", 0.1 * 3 / K.zustandsteil_je_gas(K.q_je_radiant()) * 100, 1.8),
        ("F:ganzes", "mit a₀ = 1,06 10,2 ± 1,8 %", 100 * _cmb()["hoehe_fest"][0], 10.2),
        ("F:ganzes", "mit a₀ = 1,06 10,2 ± 1,8 %", 100 * _cmb()["hoehe_fest"][1], 1.8),
        ("F:ganzes", "also F = 1,82 ± 0,32 %", 300 * _cmb()["hoehe_fest"][0] / K.zustandsteil_je_gas(K.q_je_radiant()), 1.82),
        ("F:ganzes", "also F = 1,82 ± 0,32 %", 300 * _cmb()["hoehe_fest"][1] / K.zustandsteil_je_gas(K.q_je_radiant()), 0.32),
        ("F:ganzes", "mit freiem a₀ 7,9 ± 3,2 %", 100 * _cmb()["hoehe"][0], 7.9),
        ("F:ganzes", "mit freiem a₀ 7,9 ± 3,2 %", 100 * _cmb()["hoehe"][1], 3.2),
        ("F:ganzes", "F = 1,4 ± 0,6 %", 300 * _cmb()["hoehe"][0] / K.zustandsteil_je_gas(K.q_je_radiant()), 1.4),
        ("F:ganzes", "F = 1,4 ± 0,6 %", 300 * _cmb()["hoehe"][1] / K.zustandsteil_je_gas(K.q_je_radiant()), 0.6),
        ("F:ganzes", "ist die ganze Dunkle Energie heute Formteil (a₀ = 1,085)", _cmb()["saett"]["punkt"]["a0"], 1.085),
        ("F:ganzes", "ist die Höhe 4,6 ± 1,8 %", 100 * _cmb()["hoehe_saett"][0], 4.6),
        ("F:ganzes", "ist die Höhe 4,6 ± 1,8 %", 100 * _cmb()["hoehe_saett"][1], 1.8),
        ("F:ganzes", "F = 0,81 ± 0,31 %", 300 * _cmb()["hoehe_saett"][0] / K.zustandsteil_je_gas(K.q_je_radiant()), 0.81),
        ("F:ganzes", "F = 0,81 ± 0,31 %", 300 * _cmb()["hoehe_saett"][1] / K.zustandsteil_je_gas(K.q_je_radiant()), 0.31),
        ("F:ganzes", "Baryonen um z = 0,4 als Gas", _cmb()["z_max"][0], 0.4),
        ("03", "gibt er a₀ = 1,08", K.a0_aus_schallwinkel(K.q_je_radiant()) / 1e-10, 1.08),
        ("03", "(H₀ = 68,0)", K.atmen_heute(K.a0_aus_schallwinkel(K.q_je_radiant()), Q=K.q_je_radiant())[0], 68.0),
        ("03", "Δχ² = 5,9 neben ihm", D.bao_absolut(_h0_q(K.q_je_radiant(), K.a0_aus_schallwinkel(K.q_je_radiant())), *_om_or_q(K.q_je_radiant())(K.a0_aus_schallwinkel(K.q_je_radiant())), K.schallhorizont(K.q_je_radiant())) - _bao_a0(K.q_je_radiant())[3], 5.9),
        ("03", "statt 21,4 mit dem gemessenen ω_c", D.bao_absolut(K.atmen_heute(K.a0_aus_schallwinkel())[0], *_om_or_a0(K.a0_aus_schallwinkel()), K.schallhorizont()) - _bao_a0()[3], 21.4),
        ("03", "zusammen lesen a₀ = 1,088", _ganzes_a0(False)[0], 1.088),
        ("03", "(1,084 bis 1,092)", _ganzes_a0(False)[1], 1.084),
        ("03", "1,084 bis 1,092)", _ganzes_a0(False)[2], 1.092),
        ("03", "bei H₀ = 68,2", _ganzes_a0(True)[3], 68.2),
        ("03", "und Ω_m = 0,30", _ganzes_a0(True)[4], 0.30),
        ("03", "Abstände von 79 der 141 Galaxien", _sparc_ganzes()[5], 79),
        ("03", "79 der 141 Galaxien", _sparc_ganzes()[6], 141),
        ("03", "Hubble-Fluss mit H₀ = 73", D.H0_SPARC, 73),
        ("03", "H₀ = 73: a₀ = 1,085", _sparc_ganzes()[0], 1.085),
        ("03", "(0,98 bis 1,19)", _sparc_ganzes()[1], 0.98),
        ("03", "0,98 bis 1,19), bei", _sparc_ganzes()[2], 1.19),
        ("03", "bei Υ_d = 0,66", _sparc_ganzes()[3], 0.66),
        ("03", "w₀ = −0,85", _ganzes_zustand(2 * math.pi - 1)[1], -0.85),
        ("03", "w_a = −0,49", _ganzes_zustand(2 * math.pi - 1)[2], -0.49),
        ("03", "Δχ² = 8,5 gegen den Formteil allein", _ganzes_zustand(2 * math.pi - 1)[0], 8.5),
        ("03", "liegt dann bei z = 0,43", _ganzes_maximum()[0], 0.43),
        ("03", "(0,34 bis 0,55)", _ganzes_maximum()[1], 0.34),
        ("03", "0,34 bis 0,55), zwischen", _ganzes_maximum()[2], 0.55),
        ("03", "Faltung bei z = 0,32", _ganzes_maximum()[3], 0.32),
        ("03", "beschleunigten Dehnung bei z = 0,67", _ganzes_maximum()[4], 0.67),
        ("03", "Der Standardwert 1,06", K.a0 / 1e-10, 1.06),
        ("09", "g_p = 5,586 ist Inhalt", K.g_p_CODATA, 5.586),
        ("09", "Ohne Rückkehr gibt sie 1418,95 MHz", K.hyperfein_wasserstoff(False) / 1e6, 1418.95),
        ("09", "mit der ersten 1420,60 MHz", K.hyperfein_wasserstoff() / 1e6, 1420.60),
        ("09", "gemessen sind 1420,41 MHz", K.NU_21CM / 1e6, 1420.41),
        ("09", "vom g-Faktor zu 0,88·α/2π", (K.NU_21CM / K.hyperfein_wasserstoff(False) - 1) / (K.alpha / (2 * math.pi)), 0.88),
        ("09", "Die restlichen 1,4·10⁻⁴", (K.hyperfein_wasserstoff() / K.NU_21CM - 1) / 1e-4, 1.4),
        ("F:ganzes", "schließt nach 2, 8, 20, 40, 70 (κ = 0,1", float(KE.schliesszahlen_verwindung(0.1, 0.5) == (2, 8, 20, 40, 70)), 1.0),
        ("F:ganzes", "Mit μ = 0,02 statt 0,03 erscheint 126, dann fehlt 28", float(KE.schliesszahlen_verwindung(0.12, 0.35, mu=0.02) == (2, 8, 20, 40, 50, 82, 126)), 1.0),
        ("F:ganzes", "356 Spins (42 %), 716 Paritäten (84 %)", sum(x * y for x, y in zip(KE.grundspin_pruefung(0.12, mu=0.03, paarfenster=0.1)[:2], (1000, 1))), 356716),
        ("F:ganzes", "325 Spins und 558 Paritäten (54 %, 92 %)", sum(x * y for x, y in zip(KE.grundspin_pruefung(0.12, mu=0.03, kugel=True, paarfenster=0.1)[:2], (1000, 1))), 325558),
        ("06", "ν = 1,1567", K.rc_kette_verstaerkung(20, 0.01), 1.1567),
        ("06", "s ≈ 3·10⁵", K.s_(9.81) / 1e5, 3),
        ("08", "Sonne 4,43 und 8,86 km", 9 * K.faltung(K.M_sun) / 1e3, 4.43),
        ("08", "8,86 km", 18 * K.faltung(K.M_sun) / 1e3, 8.86),
        ("02", "Über 93 Elemente", D.pubchem_fenster()[0], 93),
        ("02", "von 0,81 (As", D.pubchem_fenster()[1][1], 0.81),
        ("02", "bis 8,2 (Ga)", D.pubchem_fenster()[2][1], 8.2),
        ("02", "Ionisierungsenergie −0,48", D.pubchem_fenster()[3], -0.48),
        ("03", "im Median 0,71·r*", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[0], 0.71),
        ("03", "(16 bis 84 %: 0,37 bis 1,28)", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[1][0], 0.37),
        ("03", "0,37 bis 1,28) über vier", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[1][1], 1.28),
        ("03", "M_b^{0,31 ± 0,02}", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[2], 0.31),
        ("03", "0,31 ± 0,02}", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[3], 0.02),
        ("01", "für Wasserstoff 2,661·10⁻⁵", K.alpha**2 / 2 * K.m_p / (K.m_e + K.m_p) / 1e-5, 2.661),
        ("03", "Q = 2π − 1 = 5,283", K.q_je_radiant(), 5.283),
        ("03", "gemessen 5,36 ± 0,06", K.q_messung()[1], 0.06),
        ("03", "± 0,06, 1,3 σ", (K.q_messung()[0] - K.q_je_radiant()) / K.q_messung()[1], 1.3),
        ("01", "⁵⁶Fe trägt 0,94 % weniger", 100 * [e for e in D.ame2020() if e[0] == 26 and e[2] == 56][0][3]
            / (K.m_u * K.c**2 / K.e / 1e3), 0.94),
        ("03", "Materie: Υ_d = 0,67", D.sparc_eichung(K.kpc), 0.67),
        ("03", "liest Pantheon+ dagegen H₀ = 73,5", D.pantheon_h0(*_om_or(), 0.023, 0.15)[0], 73.5),
        ("03", "H₀ = 73,5 ± 1,0", D.pantheon_h0(*_om_or(), 0.023, 0.15)[1], 1.0),
        ("03", "(496 Supernovae, 0,023 < z < 0,15)", float(D.pantheon_h0(*_om_or(), 0.023, 0.15)[2]), 496),
        ("03", "zwischen 73,0 und 74,5", min(D.pantheon_h0(*_om_or(), lo, hi)[0] for lo, hi in ((0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0))), 73.0),
        ("03", "73,0 und 74,5. Der", max(D.pantheon_h0(*_om_or(), lo, hi)[0] for lo, hi in ((0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0))), 74.5),
        ("03", "Dipol ist ΔH/H = 0,9 %", D.pantheon_richtung(*_om_or(), 0.023)[0] * 100, 0.9),
        ("03", "0,9 % (0,7 σ)", D.pantheon_richtung(*_om_or(), 0.023)[1], 0.7),
        ("03", "erst bei z > 0,1", float(D.pantheon_richtung(*_om_or(), 0.01, 0.1)[2] < 1 < D.pantheon_richtung(*_om_or(), 0.1)[2]), 1.0),
        ("03", "Teilkatalog auf 1,8 σ", D.pantheon_richtung(*_om_or(), 0.1, 3.0, True)[2], 1.8),
        ("03", "allein, ohne freie Zahl (Ω_m = 0,315)", _om_or()[0], 0.315),
        ("03", "SPARC, 141 Galaxien", float(len(D.sparc_auswahl(D.sparc_eichung(K.kpc)))), 141),
        ("03", "141 Galaxien, 2688 Punkte", float(D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[1]), 2688),
        ("03", "In SPARC (175 Galaxien)", float(len(D.sparc_galaxien())), 175),
        ("03", "v = 249 km/s bei 5 kpc", _v_mw(5), 249),
        ("03", "214 bei r*", _v_mw(K.schwellenradius(6e10 * K.M_sun) / K.kpc), 214),
        ("03", "196 bei 15 kpc", _v_mw(15), 196),
        ("03", "178 bei 50 kpc", _v_mw(50), 178),
        ("03", "kleinste gemessene Rotation 17,8 km/s", min(D.sparc_kurve(n)["V"].max() for n in D.sparc_galaxien()), 17.8),
        ("03", "die kleinste flache 18,6 km/s", min(p["Vflat"] for p in D.sparc_galaxien().values() if p["Vflat"] > 0), 18.6),
        ("03", "Pantheon+, 1590 Supernovae", float(len(D.pantheon()[0])), 1590),
        ("03", "nur Δχ² = 2,0 über", _abst()["modell"] - _abst()["lcdm"], 2.0),
        ("03", "ΛCDM (Ω_m = 0,304)", _abst()["om_lcdm"], 0.304),
        ("03", "Ω_m = 0,304), 1,4 σ", _sigma(_abst()["modell"] - _abst()["lcdm"], 1), 1.4),
        ("03", "Der Formteil trifft es mit χ² = 6,8", B.kosmos()["chi2"], 6.8),
        ("08", "(10⁶ M☉: s² = 1,4·10¹⁷", K.c**4 / (4 * K.G * 1e6 * K.M_sun) / K.a0 / 1e17, 1.4),
        ("08", "am Horizont, r* = 36 pc)", K.schwellenradius(1e6 * K.M_sun) / K.pc, 36),
        ("08", "= 1 − ⅔√2 = 5,7 %", 100 * (1 - 2 / 3 * math.sqrt(2)), 5.7),
        ("08", "sechsmal so viel wie Kernbindung bis zum Eisen", (1 - 2 / 3 * math.sqrt(2)) / ([e for e in D.ame2020() if e[0] == 26 and e[2] == 56][0][3] / (K.m_u * K.c**2 / K.e / 1e3)), 6),
        ("08", "(j = 1) sind es 1 − 1/√3 = 42 %", 100 * (1 - 1 / math.sqrt(3)), 42),
        ("08", "M_irr = M_rot = M/√2: Das ist j = 1", K.falte_gedreht(K.M_sun, K.G * K.M_sun**2 / K.c)["M_irr"] * math.sqrt(2) / K.M_sun, 1.0),
        ("08", "M_irr·M_rot = Jc/(2G)", K.falte_gedreht(K.M_sun, 0.5 * K.G * K.M_sun**2 / K.c)["M_irr"] * K.falte_gedreht(K.M_sun, 0.5 * K.G * K.M_sun**2 / K.c)["M_rot"] / (0.5 * K.M_sun**2 / 2), 1.0),
        ("08", "4π·(r₊² + a_J²) = 4π·r₀(M_irr)²", 4 * math.pi * (K.horizonte(0.5 * K.G * K.M_sun**2 / K.c, K.M_sun)[0].real**2 + K.verwindung(0.5 * K.G * K.M_sun**2 / K.c, K.M_sun)**2) / (4 * math.pi * K.grenzradius(K.falte_gedreht(K.M_sun, 0.5 * K.G * K.M_sun**2 / K.c)["M_irr"])**2), 1.0),
        ("08", "am Treffpunkt j = 1 sind es 1 − 1/√2 = 29 %", 100 * K.falte_gedreht(K.M_sun, K.G * K.M_sun**2 / K.c)["herausholbar"], 29),
        ("08", "die Entropie ist dort halb so groß", K.falte_gedreht(K.M_sun, K.G * K.M_sun**2 / K.c)["entropie"] / K.entropie_falte(K.M_sun), 0.5),
        ("02", "P/(ρc²) = (dτ/dt|_R)/(dτ/dt|_0) − 1", K.mitte_gehalten(K.M_sun, 1.2 * K.grenzradius(K.M_sun), numerisch=True)[1] / K.mitte_gehalten(K.M_sun, 1.2 * K.grenzradius(K.M_sun))[1], 1.0),
        ("02", "schwach ihr Unterschied ½·δ, die halbe Tiefe", K.mitte_gehalten(K.M_sun, 1e6 * K.faltung(K.M_sun), numerisch=True)[1] / (1 - K.takt(K.M_sun, 1e6 * K.faltung(K.M_sun))), 0.5),
        ("02", "(3/2)·√(1 − r₀/R) − ½. Bei gleichem R", K.mitte_gehalten(1.4 * K.M_sun, 12e3, numerisch=True)[0] / (1.5 * K.takt(1.4 * K.M_sun, 12e3) - 0.5), 1.0),
        ("02", "von 1,4 M☉ und 12 km höchstens mit 0,71", K.mitte_gehalten(1.4 * K.M_sun, 12e3, numerisch=True)[0], 0.71),
        ("02", "Hülle jedes solchen Körpers schneller als ⅓", K.takt(K.M_sun, __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda R: K.mitte_gehalten(K.M_sun, R)[0], 1.0001 * K.grenzradius(K.M_sun), 3 * K.grenzradius(K.M_sun))), 1 / 3),
        ("02", "höchstens um z = 2 gerötet an", 1 / K.takt(K.M_sun, 1.125 * K.grenzradius(K.M_sun)) - 1, 2),
        ("02", "außerhalb von (9/8)·r₀ = 27Δs/4", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda R: K.mitte_gehalten(K.M_sun, R)[0], 1.0001 * K.grenzradius(K.M_sun), 3 * K.grenzradius(K.M_sun)) / K.faltung(K.M_sun), 6.75),
        ("02", "also 3Δs/R unter 4/9 gegen ½", 3 * K.faltung(K.M_sun) / (1.125 * K.grenzradius(K.M_sun)), 4 / 9),
        ("02", "bei 1,4 M☉ über 4,65 km (Buchdahl)", 9 / 8 * K.grenzradius(1.4 * K.M_sun) / 1e3, 4.65),
        ("02", "3Δs/r* = ½·√(6/(1 + Q)) = 0,49", 3 * K.faltung(1.6e80 * K.m_p) / math.sqrt(K.G * 1.6e80 * K.m_p * 2 * (1 + K.q_je_radiant()) * 1.6e80 * K.faltung(K.m_p) / K.c**2), 0.49),
        ("08", "die n-te mit 1/n²", 9 * __import__("scipy.integrate", fromlist=["quad"]).quad(lambda u: u * math.exp(-3 * u), 0, 60)[0], 1.0),
        ("08", "wächst die Ablenkung so um (π²/6)·(b/r*)² für b ≪ r*", (K.ablenkung_regel(K.M_sun, 1e-3 * K.schwellenradius(K.M_sun)) / K.ablenkung(K.M_sun, 1e-3 * K.schwellenradius(K.M_sun)) - 1) / 1e-6 / sum(1 / n**2 for n in range(1, 200000)), 1.0),
        ("08", "Bei 4D_E = 5,1 Gpc wächst die Ablenkung dort um 32 %", 100 * (K.ablenkung_regel(1e12 * K.M_sun, math.sqrt(2 * K.grenzradius(1e12 * K.M_sun) * K.einsteinring(1.0)[0])) / K.ablenkung(1e12 * K.M_sun, math.sqrt(2 * K.grenzradius(1e12 * K.M_sun) * K.einsteinring(1.0)[0])) - 1), 32),
        ("08", "Bei 4D_E = 5,1 Gpc wächst die Ablenkung dort um 32 %", 100 * ((lambda M, b: 4 / K.c**2 * (__import__("scipy.integrate", fromlist=["quad"]).quad(lambda l: K.verstaerkung(K.G * M / (b * b + l * l)) * K.G * M * b / (b * b + l * l)**1.5, 0, K.schwellenradius(M), limit=400)[0] + __import__("scipy.integrate", fromlist=["quad"]).quad(lambda l: K.verstaerkung(K.G * M / (b * b + l * l)) * K.G * M * b / (b * b + l * l)**1.5, K.schwellenradius(M), 1e4 * K.schwellenradius(M), limit=400)[0] + K.v_flach(M)**2 * b / (1e4 * K.schwellenradius(M))))(1e12 * K.M_sun, math.sqrt(2 * K.grenzradius(1e12 * K.M_sun) * K.einsteinring(1.0)[0])) / K.ablenkung(1e12 * K.M_sun, math.sqrt(2 * K.grenzradius(1e12 * K.M_sun) * K.einsteinring(1.0)[0])) - 1), 32),
        ("08", "und der Ring liegt bei 1,20·R_E", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda b: b - K.einsteinring(1.0)[0] * K.ablenkung_regel(1e10 * K.M_sun, b), 0.5 * math.sqrt(2 * K.grenzradius(1e10 * K.M_sun) * K.einsteinring(1.0)[0]), 3 * math.sqrt(2 * K.grenzradius(1e10 * K.M_sun) * K.einsteinring(1.0)[0])) / math.sqrt(2 * K.grenzradius(1e10 * K.M_sun) * K.einsteinring(1.0)[0]), 1.20),
        ("08", "Am Trichter vorbei verlängert sie sich um Δt = (r₀/c)·ln(4r₁r₂/b²)", 2 * K.G * K.M_sun / K.c**3 * (math.asinh(K.AU / K.R_SONNE) * 2) / (K.grenzradius(K.M_sun) / K.c * math.log(4 * K.AU**2 / K.R_SONNE**2)), 1.0),
        ("08", "c·|∂Δt/∂b| = 2r₀/b = 12Δs/b", K.grenzradius(K.M_sun) * (math.log(4 * K.AU**2 / (0.999999 * K.R_SONNE)**2) - math.log(4 * K.AU**2 / (1.000001 * K.R_SONNE)**2)) / (2e-6 * K.R_SONNE) / K.ablenkung(K.M_sun, K.R_SONNE), 1.0),
        ("08", "R_E²/(2c·D_E) = r₀/c", K.G * K.M_sun / K.einsteinring(0.5, 2.0)[1] / (2 * K.c * K.einsteinring(0.5, 2.0)[0]) / (K.grenzradius(K.M_sun) / K.c), 1.0),
        ("08", "bei der Sonne 9,85 µs, am Ort der Linse gezählt", 2 * K.G * K.M_sun / K.c**3 * 1e6, 9.85),
        ("08", "Ihre Summe ist die Lage der Quelle, ihr Produkt −1", (lambda w: w[0] * w[1])(__import__("numpy").roots([1, -0.5, -1]).real), -1),
        ("08", "y = √(b₊/|b₋|) − √(|b₋|/b₊)", (lambda w: math.sqrt(max(w) / abs(min(w))) - math.sqrt(abs(min(w)) / max(w)))(__import__("numpy").roots([1, -0.5, -1]).real) / 0.5, 1.0),
        ("08", "(r₀/c)·(y·√(y² + 4) + 2·ln((√(y² + 4) + y)/(√(y² + 4) − y)))", (lambda T, d: (T(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: d(x, 0.5), -5, -1e-3), 0.5) - T(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: d(x, 0.5), 1e-3, 5), 0.5)) / (0.5 * math.sqrt(4.25) + 2 * math.log((math.sqrt(4.25) + 0.5) / (math.sqrt(4.25) - 0.5))))(lambda x, y: (x - y)**2 - 2 * math.log(abs(x)), lambda x, y: 2 * (x - y) - 2 / x), 1.0),
        ("05", "für Masse und Ladung des Protons 1,1·10¹⁸", math.sqrt(K.alpha / K.kopplung_grav(K.m_p)) / 1e18, 1.1),
        ("05", "für Masse und Ladung des Protons 1,1·10¹⁸", (lambda n, T: math.sqrt(K.k_B * T / (4 * math.pi * K.G * n * K.m_p**2)) / math.sqrt(K.eps0 * K.k_B * T / (n * K.e**2)))(1e6, 1e4) / 1e18, 1.1),
        ("05", "ℓ_J/ℓ_D = √(α/α_G), bei jeder Dichte und Temperatur", (lambda n, T: math.sqrt(K.k_B * T / (4 * math.pi * K.G * n * K.m_p**2)) / math.sqrt(K.eps0 * K.k_B * T / (n * K.e**2)))(1e30, 1e7) / (lambda n, T: math.sqrt(K.k_B * T / (4 * math.pi * K.G * n * K.m_p**2)) / math.sqrt(K.eps0 * K.k_B * T / (n * K.e**2)))(1.0, 100.0), 1.0),
        ("05", "N^{2/3}·α_G = α, dieselbe Bedingung", (lambda n, T: (math.sqrt(K.k_B * T / (4 * math.pi * K.G * n * K.m_p**2)) / math.sqrt(K.eps0 * K.k_B * T / (n * K.e**2)))**2 * K.kopplung_grav(K.m_p) / K.alpha)(1.0, 100.0), 1.0),
        ("05", "ρ·ℓ_J³ = M·n·ℓ_D³ mit M = m·(α/α_G)^{3/2}", (lambda n, T: n * K.m_p * math.sqrt(K.k_B * T / (4 * math.pi * K.G * n * K.m_p**2))**3 / (K.masse_windung_faltung()[0] * n * math.sqrt(K.eps0 * K.k_B * T / (n * K.e**2))**3))(1e30, 1e7), 1.0),
        ("04", "also C = −(3/2)·Nk", (lambda si, qd: (lambda s: -1.5 * 3 / 2.5 * qd(lambda x: max(s.sol(x)[0], 0)**2.5 * x**2, 1e-6, s.t_events[0][0], limit=200)[0] / qd(lambda x: -x**3 * s.sol(x)[1] * max(s.sol(x)[0], 0)**1.5, 1e-6, s.t_events[0][0], limit=200)[0])(si.solve_ivp(lambda x, y: [y[1], -max(y[0], 0)**1.5 - 2 * y[1] / x], [1e-6, 50], [1, 0], events=(lambda f: (setattr(f, 'terminal', True), f)[1])(lambda x, y: y[0]), rtol=1e-12, atol=1e-14, dense_output=True)))(__import__('scipy.integrate', fromlist=['solve_ivp']), __import__('scipy.integrate', fromlist=['quad']).quad), -1.5),
        ("04", "also C = −(3/2)·Nk", (lambda P, e, x, h: (e(x + h) - e(x - h)) / ((x + h) / 4 - P(x + h) - (x - h) / 4 + P(x - h)))(K.druck_elektronen, lambda x: 3 / (8 * x**3) * (x * (2 * x**2 + 1) * math.sqrt(1 + x**2) - math.asinh(x)) - 1 + 1.5 * (x / 4 - K.druck_elektronen(x)) - 0.75 * x, 0.003, 1e-6), -1.5),
        ("04", "an der Grenzmasse bei ħk_F = 1,06·m_e·c", K.decke_entartet(1.0)[1], 1.06),
        ("04", "an der Grenzmasse bei ħk_F = 1,06·m_e·c", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: (K.druck_elektronen(x + 1e-6) - K.druck_elektronen(x - 1e-6)) / 2e-6 - 0.25, 0.5, 2, xtol=1e-12), 1.06),
        ("04", "weil M·m̃ = m_Pl²/2 fest ist", K.M_sun * K.spiegelmasse(K.M_sun) / K.m_Pl**2, 0.5),
        ("04", "C = −Mc²/T = −4πk·M/m̃ = −2S", (lambda M: M * 1e-6 * K.c**2 / (K.k_B * (K.falte_temperatur(M * (1 + 1e-6)) - K.falte_temperatur(M))))(K.M_sun) / K.entropie_falte(K.M_sun), -2),
        ("04", "C = −Mc²/T = −4πk·M/m̃ = −2S", (-4 * math.pi * K.M_sun / K.spiegelmasse(K.M_sun)) / K.entropie_falte(K.M_sun), -2),
        ("07", "⟨Δn²⟩ = n_B·ν = n_B + n_B² = 1/(4·sinh²(s/2))", max(abs((lambda s: (lambda q: sum(k * k * (1 - q) * q**k for k in range(20000)) - sum(k * (1 - q) * q**k for k in range(20000))**2)(math.exp(-s)))(s) / (K.rueckkehr_summe(s) * (1 + K.rueckkehr_summe(s))) - 1) for s in (0.2, math.log(2), 1.0, 3.0)) * 1e8 + max(abs(K.rueckkehr_summe(s) * (1 + K.rueckkehr_summe(s)) * 4 * math.sinh(s / 2)**2 - 1) for s in (0.2, 1.0, 3.0)) * 1e8, 0),
        ("07", "⟨Δn²⟩ = −dν/ds = −d(½·coth(s/2))/ds", max(abs((K.verstaerkung((s - 1e-5)**2 * K.a0) - K.verstaerkung((s + 1e-5)**2 * K.a0)) / 2e-5 / (lambda s: (lambda q: sum(k * k * (1 - q) * q**k for k in range(20000)) - sum(k * (1 - q) * q**k for k in range(20000))**2)(math.exp(-s)))(s) - 1) + abs((0.5 / math.tanh((s - 1e-5) / 2) - 0.5 / math.tanh((s + 1e-5) / 2)) / 2e-5 / (lambda s: (lambda q: sum(k * k * (1 - q) * q**k for k in range(20000)) - sum(k * (1 - q) * q**k for k in range(20000))**2)(math.exp(-s)))(s) - 1) for s in (0.2, 1.0, 3.0)) * 1e6, 0),
        ("07", "d ln ν/d ln s = −s·⟨Δn²⟩/ν = −χ", max(abs((math.log(K.verstaerkung(x * (1 + 1e-6)**2 * K.a0)) - math.log(K.verstaerkung(x * (1 - 1e-6)**2 * K.a0))) / (math.log(1 + 1e-6) - math.log(1 - 1e-6)) + K.s_(x * K.a0) * (lambda s: (lambda q: sum(k * k * (1 - q) * q**k for k in range(20000)) - sum(k * (1 - q) * q**k for k in range(20000))**2)(math.exp(-s)))(K.s_(x * K.a0)) / K.verstaerkung(x * K.a0)) + abs(K.s_(x * K.a0) * (lambda s: (lambda q: sum(k * k * (1 - q) * q**k for k in range(20000)) - sum(k * (1 - q) * q**k for k in range(20000))**2)(math.exp(-s)))(K.s_(x * K.a0)) / K.verstaerkung(x * K.a0) / K.schalter(x * K.a0) - 1) for x in (0.1, 0.7, 4.0)) * 1e6, 0),
        ("07", "q = ½, bei a_N = (ln 2)²·a₀ = 0,48·a₀", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: K.verstaerkung(x * K.a0) - 2, 0.1, 3, xtol=1e-14), 0.48),
        ("07", "q = ½, bei a_N = (ln 2)²·a₀ = 0,48·a₀", math.log(2)**2, 0.48),
        ("07", "so viel wie Newton, ν = 2", 1 / (1 - K.rueckkehranteil(math.log(2)**2 * K.a0)), 2),
        ("07", "n_F·(1 − n_F) = −dn_F/ds ≤ ¼", max(K.rueckkehr_summe(x / 1000, 2) * (1 - K.rueckkehr_summe(x / 1000, 2)) for x in range(0, 5001)), 0.25),
        ("07", "n_F·(1 − n_F) = −dn_F/ds ≤ ¼", max(abs((K.rueckkehr_summe(s - 1e-5, 2) - K.rueckkehr_summe(s + 1e-5, 2)) / 2e-5 / (K.rueckkehr_summe(s, 2) * (1 - K.rueckkehr_summe(s, 2))) - 1) for s in (0.2, 1.0, 3.0)) * 1e6, 0),
        ("03", "= π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³", math.pi**3 * K.a0**2 / (5 * K.G * K.c**2) / (K.M_sun / K.Mpc**3) / 1e11, 1.7),
        ("03", "= π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³", (lambda R: 3 * K.v_flach(6e10 * K.M_sun)**2 * R / K.G / (4 * math.pi * R**3))(K.groesste_umkehr(6e10 * K.M_sun)[0]) / (K.M_sun / K.Mpc**3) / 1e11, 1.7),
        ("03", "= π³a₀²/(5Gc²) = 1,7·10¹¹ M☉/Mpc³", 2 * K.dichte_dunkle_energie() / (K.M_sun / K.Mpc**3) / 1e11, 1.7),
        ("03", "2b₁·k_s = π³/5 = 6,20", 2 * K.dichte_dunkle_energie() * K.G * K.c**2 / K.a0**2, 6.20),
        ("03", "2b₁·k_s = π³/5 = 6,20", math.pi**3 / 5, 6.20),
        ("03", "3H_∞²/(8πG) = b₁·ρ_s", 3 * K.groesste_umkehr(K.M_sun)[1]**2 / (8 * math.pi * K.G) / (K.B1_ZELLE * K.dichte_vakuumteil()), 1.0),
        ("03", "2π·max_x (1 − e^{−x}I₀(x))/√x = 3,36 statt π", -2 * math.pi * __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda x: -(1 - __import__("scipy.special", fromlist=["ive"]).ive(0, x)) / math.sqrt(x), bounds=(0.01, 10.0), method="bounded", options={"xatol": 1e-10}).fun, 3.36),
        ("03", "2π·max_x (1 − e^{−x}I₀(x))/√x = 3,36 statt π", 2 * math.pi * max((1 - math.exp(-x / 1000) * __import__("scipy.special", fromlist=["i0"]).i0(x / 1000)) / math.sqrt(x / 1000) for x in range(10, 10000)), 3.36),
        ("03", "für eine Scheibe aus Gas bei R_d 0,259·v", K.toomre_schwelle(1.0)["sigma_v"], 0.259),
        ("03", "für eine Scheibe aus Gas bei R_d 0,259·v", K.toomre_schwelle(1.0, numerisch=True)["sigma_v"], 0.259),
        ("03", "bei 2R_d 0,144·v (dort", K.toomre_schwelle(2.0)["sigma_v"], 0.144),
        ("03", "bei 2R_d 0,144·v (dort", K.toomre_schwelle(2.0, numerisch=True)["sigma_v"], 0.144),
        ("03", "(dort ν = 2,16, längs", K.toomre_schwelle(2.0)["nu"], 2.16),
        ("03", "längs ν(1 − χ/2) = 1,382)", K.toomre_schwelle(2.0, numerisch=True)["laengs"], 1.382),
        ("03", "Streut das Gas mit mehr als 0,224·v_flach", K.toomre_schwelle(1.0)["sigma_vf"], 0.224),
        ("03", "Streut das Gas mit mehr als 0,224·v_flach", max(K.toomre_schwelle(1 + i / 50)["sigma_vf"] for i in range(451)), 0.224),
        ("F:galaxien", "liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD", D.sparc_spiralgrenze()[0], 26),
        ("F:galaxien", "liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD", D.sparc_spiralgrenze()[1], 26),
        ("F:galaxien", "liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD", D.sparc_spiralgrenze(rand=True)[0], 20),
        ("F:galaxien", "liegen in SPARC 26 Galaxien über dieser Grenze, vom Rand von Hand 20, alle vom Typ Sm, Im oder BCD", D.sparc_spiralgrenze(rand=True)[1], 20),
        ("F:galaxien", "am nächsten UGC 1281 (Sdm) mit 0,22·v_flach, von Hand 0,19·v_flach", D.sparc_spiralgrenze()[3], 0.22),
        ("F:galaxien", "am nächsten UGC 1281 (Sdm) mit 0,22·v_flach, von Hand 0,19·v_flach", D.sparc_spiralgrenze(rand=True)[3], 0.19),
        ("03", "bei 2R_d mindestens 26,7 km/s", K.toomre_schwelle(2.0, 6e10 * K.M_sun)["sigma"] / 1e3, 26.7),
        ("03", "bei 2R_d mindestens 26,7 km/s", K.toomre_schwelle(2.0, 6e10 * K.M_sun, numerisch=True)["sigma"] / 1e3, 26.7),
        ("03", "Er rechnet bei 2R_d 0,104·v", K.toomre_schwelle(2.0, numerisch=True)["halo_v"], 0.104),
        ("03", "(19,3 km/s) und liest Q_T", K.toomre_schwelle(2.0, 6e10 * K.M_sun)["halo"] / 1e3, 19.3),
        ("03", "für die Milchstraße ohne Nachbarn 2,2 Mpc", K.letzte_stabile_bahn(6e10 * K.M_sun)[0] / K.Mpc, 2.2),
        ("03", "für die Milchstraße ohne Nachbarn 2,2 Mpc", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: (3 + (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 0.9999))) / 2e-4) * (lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r) / r - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 0.5 * K.Mpc, 3 * K.Mpc) / K.Mpc, 2.2),
        ("03", "Um eine Galaxie, deren Mode übergeben ist (d_eff = 2)", 1 - (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: (3 + (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 0.9999))) / 2e-4) * (lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r) / r - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 0.5 * K.Mpc, 3 * K.Mpc) * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: (3 + (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 0.9999))) / 2e-4) * (lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r) / r - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 0.5 * K.Mpc, 3 * K.Mpc) * 0.9999))) / 2e-4, 2),
        ("03", "Ω = H_∞, einmal in", math.sqrt((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: (3 + (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 0.9999))) / 2e-4) * (lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r) / r - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 0.5 * K.Mpc, 3 * K.Mpc)) / __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: (3 + (math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 1.0001)) - math.log((lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r * 0.9999))) / 2e-4) * (lambda x: K.verstaerkung(K.G * 6e10 * K.M_sun / x**2) * K.G * 6e10 * K.M_sun / x**2)(r) / r - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 0.5 * K.Mpc, 3 * K.Mpc) - K.groesste_umkehr(K.M_sun)[1]**2) / K.groesste_umkehr(K.M_sun)[1], 1),
        ("03", "Ω = H_∞, einmal in", K.letzte_stabile_bahn(6e10 * K.M_sun)[1] / K.groesste_umkehr(K.M_sun)[1], 1),
        ("03", "2π/H_∞ = Θ_∞·τ_D = 110 Mrd. Jahren", 2 * math.pi / K.groesste_umkehr(K.M_sun)[1] / K.Gyr, 110),
        ("03", "2π/H_∞ = Θ_∞·τ_D = 110 Mrd. Jahren", math.sqrt(15) / math.pi * K.rueckkehrzeit() / K.Gyr, 110),
        ("03", "bei M₂₀₀ = 10¹⁵ M☉ 7,0 Mpc, mit", K.letzte_stabile_bahn(1e15 * K.M_sun, regel=False)[0] / K.Mpc, 7.0),
        ("03", "bei M₂₀₀ = 10¹⁵ M☉ 7,0 Mpc, mit", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.G * 1e15 * K.M_sun / r**3 - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 1 * K.Mpc, 20 * K.Mpc) / K.Mpc, 7.0),
        ("03", "mit Ω = √3·H_∞. Zwischen", (lambda r: math.sqrt(K.G * 1e15 * K.M_sun / r**3 - K.groesste_umkehr(K.M_sun)[1]**2) / K.groesste_umkehr(K.M_sun)[1])(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.G * 1e15 * K.M_sun / r**3 - 4 * K.groesste_umkehr(K.M_sun)[1]**2, 1 * K.Mpc, 20 * K.Mpc)), 1.732),
        ("02", "Das Mitfalten senkt die Grenzmasse nur auf 1,417 M☉", K.tov_sperre(K.m_e, None, 2 * K.m_u)[0] / K.M_sun, 1.417),
        ("02", "wenn der Halter in der Mitte schon mit 0,999·c läuft", (lambda x: x / math.sqrt(1 + x * x))(K.tov_sperre(K.m_e, None, 2 * K.m_u)[3]), 0.999),
        ("02", "Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉", K.tov_sperre(K.m_n)[0] / ((K.hbar * K.c / K.G)**1.5 / K.m_n**2), 0.384),
        ("02", "Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉", K.tov_sperre(K.m_n)[0] / K.M_sun, 0.710),
        ("02", "Die Grenzmasse ist 0,384·m_Pl³/m_n² = 0,710 M☉", (lambda m: (lambda A: (lambda eps, P, xh: (lambda st: (lambda xc: st(xc) + (xc,))(__import__('scipy.optimize', fromlist=['o']).minimize_scalar(lambda x: -st(x)[0], bounds=(0.6, 1.2), method='bounded', options=dict(xatol=1e-7)).x))(lambda xc: (lambda s: (s.y[1][-1], s.y[0][-1]))(__import__('scipy.integrate', fromlist=['i']).solve_ivp(lambda h, y: (lambda x, dr: [dr, 4*math.pi*y[0]**2*eps(x)/K.c**2*dr])(xh(h), -y[0]*(y[0] - 2*K.G*y[1]/K.c**2)*K.c**2/(K.G*(y[1] + 4*math.pi*y[0]**3*P(xh(h))/K.c**2))), (math.log(math.sqrt(1 + xc*xc)), 0), [1e-3, 4*math.pi/3*1e-9*eps(xc)/K.c**2], rtol=1e-11, atol=0, method='LSODA'))))(lambda x: A/8*(x*(2*x*x + 1)*math.sqrt(1 + x*x) - math.asinh(x)), lambda x: A/24*(x*(2*x*x - 3)*math.sqrt(1 + x*x) + 3*math.asinh(x)), lambda h: math.sqrt(max(math.exp(2*h) - 1, 0))))(m**4*K.c**5/(math.pi**2*K.hbar**3)))(K.m_n)[0] / K.M_sun, 0.710),
        ("02", "mit R = 9,16 km und 3Δs/R = 0,114", K.tov_sperre(K.m_n)[1] / 1e3, 9.16),
        ("02", "mit R = 9,16 km und 3Δs/R = 0,114", (lambda m: (lambda A: (lambda eps, P, xh: (lambda st: (lambda xc: st(xc) + (xc,))(__import__('scipy.optimize', fromlist=['o']).minimize_scalar(lambda x: -st(x)[0], bounds=(0.6, 1.2), method='bounded', options=dict(xatol=1e-7)).x))(lambda xc: (lambda s: (s.y[1][-1], s.y[0][-1]))(__import__('scipy.integrate', fromlist=['i']).solve_ivp(lambda h, y: (lambda x, dr: [dr, 4*math.pi*y[0]**2*eps(x)/K.c**2*dr])(xh(h), -y[0]*(y[0] - 2*K.G*y[1]/K.c**2)*K.c**2/(K.G*(y[1] + 4*math.pi*y[0]**3*P(xh(h))/K.c**2))), (math.log(math.sqrt(1 + xc*xc)), 0), [1e-3, 4*math.pi/3*1e-9*eps(xc)/K.c**2], rtol=1e-11, atol=0, method='LSODA'))))(lambda x: A/8*(x*(2*x*x + 1)*math.sqrt(1 + x*x) - math.asinh(x)), lambda x: A/24*(x*(2*x*x - 3)*math.sqrt(1 + x*x) + 3*math.asinh(x)), lambda h: math.sqrt(max(math.exp(2*h) - 1, 0))))(m**4*K.c**5/(math.pi**2*K.hbar**3)))(K.m_n)[1] / 1e3, 9.16),
        ("02", "mit R = 9,16 km und 3Δs/R = 0,114", K.tov_sperre(K.m_n)[2], 0.114),
        ("02", "mit R = 9,16 km und 3Δs/R = 0,114", (lambda w: K.G * w[0] / (w[1] * K.c**2))((lambda m: (lambda A: (lambda eps, P, xh: (lambda st: (lambda xc: st(xc) + (xc,))(__import__('scipy.optimize', fromlist=['o']).minimize_scalar(lambda x: -st(x)[0], bounds=(0.6, 1.2), method='bounded', options=dict(xatol=1e-7)).x))(lambda xc: (lambda s: (s.y[1][-1], s.y[0][-1]))(__import__('scipy.integrate', fromlist=['i']).solve_ivp(lambda h, y: (lambda x, dr: [dr, 4*math.pi*y[0]**2*eps(x)/K.c**2*dr])(xh(h), -y[0]*(y[0] - 2*K.G*y[1]/K.c**2)*K.c**2/(K.G*(y[1] + 4*math.pi*y[0]**3*P(xh(h))/K.c**2))), (math.log(math.sqrt(1 + xc*xc)), 0), [1e-3, 4*math.pi/3*1e-9*eps(xc)/K.c**2], rtol=1e-11, atol=0, method='LSODA'))))(lambda x: A/8*(x*(2*x*x + 1)*math.sqrt(1 + x*x) - math.asinh(x)), lambda x: A/24*(x*(2*x*x - 3)*math.sqrt(1 + x*x) + 3*math.asinh(x)), lambda h: math.sqrt(max(math.exp(2*h) - 1, 0))))(m**4*K.c**5/(math.pi**2*K.hbar**3)))(K.m_n)), 0.114),
        ("02", "der Halter läuft in der Mitte mit 0,64·c", (lambda x: x / math.sqrt(1 + x * x))(K.tov_sperre(K.m_n)[3]), 0.64),
        ("02", "der Halter läuft in der Mitte mit 0,64·c", (lambda x: x / math.sqrt(1 + x * x))((lambda m: (lambda A: (lambda eps, P, xh: (lambda st: (lambda xc: st(xc) + (xc,))(__import__('scipy.optimize', fromlist=['o']).minimize_scalar(lambda x: -st(x)[0], bounds=(0.6, 1.2), method='bounded', options=dict(xatol=1e-7)).x))(lambda xc: (lambda s: (s.y[1][-1], s.y[0][-1]))(__import__('scipy.integrate', fromlist=['i']).solve_ivp(lambda h, y: (lambda x, dr: [dr, 4*math.pi*y[0]**2*eps(x)/K.c**2*dr])(xh(h), -y[0]*(y[0] - 2*K.G*y[1]/K.c**2)*K.c**2/(K.G*(y[1] + 4*math.pi*y[0]**3*P(xh(h))/K.c**2))), (math.log(math.sqrt(1 + xc*xc)), 0), [1e-3, 4*math.pi/3*1e-9*eps(xc)/K.c**2], rtol=1e-11, atol=0, method='LSODA'))))(lambda x: A/8*(x*(2*x*x + 1)*math.sqrt(1 + x*x) - math.asinh(x)), lambda x: A/24*(x*(2*x*x - 3)*math.sqrt(1 + x*x) + 3*math.asinh(x)), lambda h: math.sqrt(max(math.exp(2*h) - 1, 0))))(m**4*K.c**5/(math.pi**2*K.hbar**3)))(K.m_n)[2]), 0.64),
        ("02", "ω₃·(√(3π)/2)·m_Pl³/m_n² = 5,73 M☉", K.grenzmasse(mu_e=K.m_n / K.m_u) / K.M_sun, 5.73),
        ("02", "ω₃·(√(3π)/2)·m_Pl³/m_n² = 5,73 M☉", K.tov_sperre(K.m_n, 300.0, faltend=False)[0] / K.M_sun, 5.73),
        ("02", "das Mitfalten schneidet sie um den Faktor 8", K.grenzmasse(mu_e=K.m_n / K.m_u) / K.tov_sperre(K.m_n)[0], 8),
        ("02", "ist 0,384 eine reine Formzahl, gleich für jede Masse des Trägers", K.tov_sperre(K.m_p)[0] / ((K.hbar * K.c / K.G)**1.5 / K.m_p**2), 0.384),
        ("02", "ist 0,384 eine reine Formzahl, gleich für jede Masse des Trägers", K.tov_sperre(K.m_e)[0] / ((K.hbar * K.c / K.G)**1.5 / K.m_e**2), 0.384),
        ("02", "M = 9,726·√(1 − β)/(μ²β²)·m_Pl³/m_H²", 4 * K.lane_emden(3.0)[1] * math.sqrt(45) / math.pi**1.5, 9.726),
        ("02", "M = 9,726·√(1 − β)/(μ²β²)·m_Pl³/m_H²", 4 * 2.01824 * math.sqrt(45) / math.pi**1.5, 9.726),
        ("02", "Die 45 ist die 3 des Drittels mal der 15", 3 * math.pi**4 / (math.gamma(4) * sum(1 / k**4 for k in range(1, 100000))), 45),
        ("02", "Die 45 ist die 3 des Drittels mal der 15", 3 * math.pi**4 / K.strahlung_je_sinn(1, 3), 45),
        ("02", "Licht und Gas halten gleich, β = ½, bei 146 M☉", K.grenzmasse(beta=0.5, mu=(1 / (2 * (1 - K.YP_GEMESSEN[0]) + 0.75 * K.YP_GEMESSEN[0]))) / K.M_sun, 146),
        ("02", "Licht und Gas halten gleich, β = ½, bei 146 M☉", 9.726 * math.sqrt(0.5) / ((1 / (2 * (1 - K.YP_GEMESSEN[0]) + 0.75 * K.YP_GEMESSEN[0]))**2 * 0.25) * (K.hbar * K.c / K.G)**1.5 / K.m_H**2 / K.M_sun, 146),
        ("02", "(Y_p = 0,2449, μ", K.YP_GEMESSEN[0], 0.2449),
        ("02", "μ = 1/(2X + 3Y_p/4) = 0,590", (1 / (2 * (1 - K.YP_GEMESSEN[0]) + 0.75 * K.YP_GEMESSEN[0])), 0.590),
        ("02", "Die Zählregel der Elektronen zählt nach Fermi, b = (3π²)^{1/3}/4", K.grenzmasse() / (K.lane_emden(3.0)[1] * math.sqrt(3 * math.pi) / 2 * (K.hbar * K.c / K.G)**1.5 / (2 * K.m_u)**2), 1.0),
        ("08", "um einen Haufen mit (1 + Q)·M_b = 10¹⁵ M☉ 6,4·10⁻⁶ unter eins", 1e6 * (1 - K.falte_endzustand(1e15 * K.M_sun)["takt_max"]), 6.4),
        ("08", "um einen Haufen mit (1 + Q)·M_b = 10¹⁵ M☉ 6,4·10⁻⁶ unter eins", 1e6 * (lambda M: (lambda R: 1 - math.sqrt(1 - 6 * K.faltung(M) / R - (K.groesste_umkehr(K.M_sun)[1] * R / K.c)**2))((K.G * M / K.groesste_umkehr(K.M_sun)[1]**2)**(1 / 3)))(1e15 * K.M_sun), 6.4),
        ("08", "√(1 − 9Δs/R_max) = √(1 − (M/M_N)^{2/3})", 1e9 * abs((lambda M: 9 * K.faltung(M) / (K.G * M / K.groesste_umkehr(K.M_sun)[1]**2)**(1 / 3) / (M / (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])))**(2 / 3))(1e15 * K.M_sun) - 1), 0),
        ("08", "Die Bahnuhr bleibt √(1 − 9Δs/r), die Photonensphäre", 1e3 * abs((lambda M, r: (lambda f: (1 - (f(r) - r / 2 * (f(r * (1 + 1e-6)) - f(r * (1 - 1e-6))) / (2e-6 * r))) / (9 * K.faltung(M) / r) - 1)(lambda x: 1 - 6 * K.faltung(M) / x - (K.groesste_umkehr(K.M_sun)[1] * x / K.c)**2))(1e15 * K.M_sun, 5e23)), 0),
        ("08", "Ihre Summe ist null und ihre Quadratsumme 2", 1e9 * abs(sum(__import__("numpy").roots([1, 0, -1, 6 * K.faltung(1e22 * K.M_sun) * K.groesste_umkehr(K.M_sun)[1] / K.c]).real)) + sum(x**2 for x in __import__("numpy").roots([1, 0, -1, 6 * K.faltung(1e22 * K.M_sun) * K.groesste_umkehr(K.M_sun)[1] / K.c]).real), 2),
        ("08", "Umkreisradius 2/√3", math.sqrt(2 / 3 * sum(x**2 for x in __import__("numpy").roots([1, 0, -1, 6 * K.faltung(1e22 * K.M_sun) * K.groesste_umkehr(K.M_sun)[1] / K.c]).real)) * math.sqrt(3) / 2, 1),
        ("08", "mit cos 3θ = −M/M_N, die Projektionen", 1e9 * abs(math.cos(3 * math.acos(max(__import__("numpy").roots([1, 0, -1, 6 * K.faltung(1e22 * K.M_sun) * K.groesste_umkehr(K.M_sun)[1] / K.c]).real) * math.sqrt(3) / 2)) + 1e22 * K.M_sun / (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1]))), 0),
        ("08", "von θ = 30° ohne Masse", math.degrees(math.acos(max(__import__("numpy").roots([1, 0, -1, 0]).real) * math.sqrt(3) / 2)), 30),
        ("08", "bis 60° bei M_N dreht", math.degrees(math.acos(K.falte_endzustand(0.999999 * (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])))["r_ganzes"] / K.c * K.groesste_umkehr(K.M_sun)[1] * math.sqrt(3) / 2)), 60),
        ("08", "Ohne Masse sind sie u = 0 und ±1", 1e12 * sum(abs(a - b) for a, b in zip(sorted(__import__("numpy").roots([1, 0, -1, 0]).real), [-1, 0, 1])), 0),
        ("08", "r_F² + r_F·r_H + r_H² = (c/H_∞)²", 1e9 * abs((lambda e: (e["r_falte"]**2 + e["r_falte"] * e["r_ganzes"] + e["r_ganzes"]**2) * K.groesste_umkehr(K.M_sun)[1]**2 / K.c**2)(K.falte_endzustand(1e22 * K.M_sun)) - 1), 0),
        ("08", "r_F² + r_F·r_H + r_H² = (c/H_∞)²", 1e9 * abs((lambda w: (w[1]**2 + w[1] * w[2] + w[2]**2) * K.groesste_umkehr(K.M_sun)[1]**2 / K.c**2)(sorted(__import__("numpy").roots([1, 0, -(K.c / K.groesste_umkehr(K.M_sun)[1])**2, 6 * K.faltung(1e22 * K.M_sun) * (K.c / K.groesste_umkehr(K.M_sun)[1])**2]).real)) - 1), 0),
        ("08", "M_N = c³/(3√3·G·H_∞) = (√5/6π²)·M_U", 1e12 * abs((K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])) / K.schnittmassen()[3] - math.sqrt(5) / (6 * math.pi**2)), 0),
        ("08", "= 2,2·10²² M☉,     r_N", K.falte_endzustand(K.M_sun)["M_N"] / K.M_sun / 1e22, 2.2),
        ("08", "r_N = c/(√3·H_∞) = 3,11 Gpc", K.falte_endzustand(K.M_sun)["r_N"] / K.Gpc, 3.11),
        ("08", "fallen Horizont der Falte, Photonensphäre, R_max und Horizont des Ganzen in r_N zusammen", 1e6 * (lambda e: max(abs(9 * K.faltung(e["M_N"]) / e["r_N"] - 1), abs((K.G * e["M_N"] / K.groesste_umkehr(K.M_sun)[1]**2)**(1 / 3) / e["r_N"] - 1), abs(K.falte_endzustand(e["M_N"] * (1 - 1e-12))["r_falte"] / e["r_N"] - 1) * 1e-2))(K.falte_endzustand(K.M_sun)), 0),
        ("08", "H_∞·τ_D = 2π/Θ_∞ = 2π²/√15", 1e9 * abs(K.groesste_umkehr(K.M_sun)[1] * K.rueckkehrzeit() / (2 * math.pi**2 / math.sqrt(15)) - 1), 0),
        ("08", "bleibt die Wärmelänge 2π·r_N", K.falte_endzustand((K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])) * (1 - 1e-9))["T_falte"] * 2 * math.pi * K.k_B * (K.c / (math.sqrt(3) * K.groesste_umkehr(K.M_sun)[1])) / (K.hbar * K.c), 1),
        ("08", "bleibt die Wärmelänge 2π·r_N", K.falte_endzustand((K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])) * (1 - 1e-9))["T_ganzes"] * 2 * math.pi * K.k_B * (K.c / (math.sqrt(3) * K.groesste_umkehr(K.M_sun)[1])) / (K.hbar * K.c), 1),
        ("08", "Gegen t werden dort beide kalt", (lambda e, u0: (lambda w: abs(u0 / w[1]**2 - 2 * w[1]) / math.sqrt(3))(sorted(__import__("numpy").roots([1, 0, -1, u0]).real)))(1e-9, 6 * K.faltung((1 - 1e-9) * (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1]))) * K.groesste_umkehr(K.M_sun)[1] / K.c), 0),
        ("08", "Der Radius r_N liegt bei 0,58·r*", (K.c / (math.sqrt(3) * K.groesste_umkehr(K.M_sun)[1])) / K.schwellenradius((K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1]))), 0.58),
        ("08", "um ν = 1,22 verstärkte", K.verstaerkung(K.G * (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])) / (K.c / (math.sqrt(3) * K.groesste_umkehr(K.M_sun)[1]))**2), 1.22),
        ("08", "in der Ecke bei M_U zusammen", 1e12 * max(abs(K.leiter_treffen((1, 0), (1 / 3, 2 / 3))[1] - 1), abs(K.leiter_treffen((1 / 3, 2 / 3), (0, 1))[1] - 1)), 0),
        ("08", "M_N liegt davor, bei x = 0,98", K.leiter_ort((K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1]))), 0.98),
        ("08", "M_N liegt davor, bei x = 0,98", 1 + math.log(math.sqrt(5) / (6 * math.pi**2)) / math.log(1 / K.leiterzahl()), 0.98),
        ("08", "liegt 6,6-mal darüber", K.schnittmassen()[3] / 4 / (K.c**3 / (3 * math.sqrt(3) * K.G * K.groesste_umkehr(K.M_sun)[1])), 6.6),
        ("08", "liegt 6,6-mal darüber", 1.5 * math.pi**2 / math.sqrt(5), 6.6),
        ('03', "ε_b = ρ_b·c²/n_γ = 0,574 eV", K.omega_b * K.rho_krit100 * K.c**2 / (K.strahlung_je_sinn(1, 2) / math.pi**2 * (K.k_B * K.T_CMB / (K.hbar * K.c))**3) / K.e, 0.574),
        ('03', "ε_b = ρ_b·c²/n_γ = 0,574 eV", K.omega_b * K.rho_krit100 / K.dichte_photonen() * K.strahlung_je_sinn(1, 3) / K.strahlung_je_sinn(1, 2) * K.k_B * K.T_CMB / K.e, 0.574),
        ('03', "je Photon (π⁴/15)/(2ζ(3))·kT = 2,70·kT", K.strahlung_je_sinn(1, 3) / K.strahlung_je_sinn(1, 2), 2.7),
        ('03', "je Photon (π⁴/15)/(2ζ(3))·kT = 2,70·kT", math.pi**4 / 15 / (2 * __import__('scipy.special', fromlist=['zeta']).zeta(3)), 2.7),
        ('03', "N_eff = 3,046", K.N_eff, 3.046),
        ('03', "je 7/8 bei (4/11)^{1/3} der Temperatur", K.strahlung_je_sinn(2, 3) / K.strahlung_je_sinn(1, 3), 0.875),
        ('03', "je 7/8 bei (4/11)^{1/3} der Temperatur", 1 / K.licht_gegen_neutrinos(), 0.714),
        ('03', "das 1,69-fache", 1 + 7 / 8 * K.N_eff * (4 / 11)**(4 / 3), 1.69),
        ('03', "das 1,69-fache", K.dichte_heute()[0] / K.dichte_photonen(), 1.69),
        ('03', "bei kT_eq = 0,80 eV", (1 + K.Q_kosm) * K.omega_b * K.rho_krit100 * K.c**2 / (K.strahlung_je_sinn(1, 2) / math.pi**2 * (K.k_B * K.T_CMB / (K.hbar * K.c))**3) / (K.strahlung_je_sinn(1, 3) / K.strahlung_je_sinn(1, 2) * (1 + 7 / 8 * K.N_eff * (4 / 11)**(4 / 3))) / K.e, 0.8),
        ('03', "bei kT_eq = 0,80 eV", (1 + K.z_gleichheit()) * K.k_B * K.T_CMB / K.e, 0.8),
        ('03', "ℓ = 2,0·10³·λ_H", (K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]) / (K.omega_b * K.rho_krit100))**(1 / 3) * K.k_B * K.T_CMB / (K.hbar * K.c) / 1e3, 2.0),
        ('03', "ℓ = 2,0·10³·λ_H", (math.pi**2 / K.strahlung_je_sinn(1, 2) * K.dichte_photonen() * K.c**2 / (K.strahlung_je_sinn(1, 3) / K.strahlung_je_sinn(1, 2) * K.k_B * K.T_CMB) / (K.omega_b * K.rho_krit100 * (1 - 0.75 * K.YP_GEMESSEN[0]) / K.m_H))**(1 / 3) / 1e3, 2.0),
        ('03', "μ = 1/(X + Y_p/4) = 1,225", 1 / (1 - 0.75 * K.YP_GEMESSEN[0]), 1.225),
        ('03', "μ = 1/(X + Y_p/4) = 1,225", 1 / ((1 - K.YP_GEMESSEN[0]) + K.YP_GEMESSEN[0] / 4), 1.225),
        ('03', "mit 1 + Q gelesen M_J = 1,0·10³ M☉", (lambda m, l: l**1.5 * K.m_Pl**3 / ((4 * math.pi)**1.5 * math.sqrt(1 + K.Q_kosm) * m**2) / K.M_sun / 1e3)(K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]), (K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]) / (K.omega_b * K.rho_krit100))**(1 / 3) * K.k_B * K.T_CMB / (K.hbar * K.c)), 1.0),
        ('03', "mit 1 + Q gelesen M_J = 1,0·10³ M☉", (lambda a, m: K.dichte_heute()[2] / a**3 * (K.k_B * K.T_CMB / a / (4 * math.pi * K.G * K.dichte_heute()[2] / a**3 * m))**1.5 / K.M_sun / 1e3)(1 / 300, K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0])), 1.0),
        ('03', "und die Jeans-Masse ist zeitfest", (lambda F: F(1 / 1100) / F(1 / 200))(lambda a, m=K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]): K.dichte_heute()[2] / a**3 * (K.k_B * K.T_CMB / a / (4 * math.pi * K.G * K.dichte_heute()[2] / a**3 * m))**1.5), 1.0),
        ('03', "an der letzten Streufläche ist M_J = 2,8·10¹⁶ M☉", (lambda a, R: K.dichte_heute()[2] / a**3 * (K.c**2 / (3 * (1 + R)) / (4 * math.pi * K.G * K.dichte_heute()[2] / a**3))**1.5 / K.M_sun / 1e16)(1 / (1 + K.z_stern), 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern))), 2.8),
        ('03', "an der letzten Streufläche ist M_J = 2,8·10¹⁶ M☉", (lambda m, R: (m * K.c**2 / (3 * (1 + R) * K.k_B * K.T_CMB * (1 + K.z_stern)))**1.5 * ((m / (K.omega_b * K.rho_krit100))**(1 / 3) * K.k_B * K.T_CMB / (K.hbar * K.c))**1.5 * K.m_Pl**3 / ((4 * math.pi)**1.5 * math.sqrt(1 + K.Q_kosm) * m**2) / K.M_sun / 1e16)(K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]), 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern))), 2.8),
        ('03', "ℓ_J = 0,43·c/H", (lambda z, R, r: math.sqrt(K.c**2 / (3 * (1 + R)) / (4 * math.pi * K.G * r[2] * (1 + z)**3)) * math.sqrt(8 * math.pi * K.G / 3 * (r[2] * (1 + z)**3 + r[0] * (1 + z)**4)) / K.c)(K.z_stern, 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern)), K.dichte_heute()), 0.43),
        ('03', "ℓ_J = 0,43·c/H", (lambda z, R, r: math.sqrt(2 * (1 + r[0] * (1 + z) / r[2]) / (9 * (1 + R))))(K.z_stern, 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern)), K.dichte_heute()), 0.43),
        ('03', "= 2,8·10¹³", (lambda m, R: (m * K.c**2 / (3 * (1 + R) * K.k_B * K.T_CMB * (1 + K.z_stern)))**1.5 / 1e13)(K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]), 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern))), 2.8),
        ('03', "= 2,8·10¹³", (lambda m, R, r, a: r / a**3 * (K.c**2 / (3 * (1 + R)) / (4 * math.pi * K.G * r / a**3))**1.5 / (r * 300**3 * (K.k_B * K.T_CMB * 300 / (4 * math.pi * K.G * r * 300**3 * m))**1.5) / 1e13)(K.m_H / (1 - 0.75 * K.YP_GEMESSEN[0]), 3 * K.omega_b * K.rho_krit100 / (4 * K.dichte_photonen() * (1 + K.z_stern)), K.dichte_heute()[2], 1 / (1 + K.z_stern)), 2.8),
        ('08', "398 M☉/pc²", K.flaeche_strahl()*K.pc**2/K.M_sun, 398),
        ('08', "398 M☉/pc²", __import__("scipy.integrate", fromlist=["quad"]).quad(lambda u: (K.verstaerkung(K.a0*math.exp(-2*u)) - 1)*K.a0*math.exp(-2*u), -8, 12, points=[0], limit=400)[0]/(math.pi*K.G)*K.pc**2/K.M_sun, 398),
        ('08', "Linse bei z = 0,5", 0.5, 0.5),
        ('08', "Quelle bei z = 2", 2.0, 2),
        ('08', "ist κ = 0,20", K.flaeche_strahl(K.einsteinring(0.5, 2.0)[0])[1], 0.2),
        ('08', "ist κ = 0,20", math.pi**2/6*2*K.einsteinring(0.5, 2.0)[0]/K.schliesslaenge(), 0.2),
        ('08', "die Bilder lesen 25 % mehr Faltung", 100*(1/(1 - K.flaeche_strahl(K.einsteinring(0.5, 2.0)[0])[1]) - 1), 25),
        ('08', "die Bilder lesen 25 % mehr Faltung", 100*(1/(1 - math.pi**2/6*2*K.einsteinring(0.5, 2.0)[0]/K.schliesslaenge()) - 1), 25),
        ('08', "−½c² und dem Drehimpuls J/M je Masse", (lambda M,J: (0.5*(J/M/K.horizonte(J,M)[0].real)**2-K.G*M/K.horizonte(J,M)[0].real)/K.c**2)(10*K.M_sun,0.8*K.G*(10*K.M_sun)**2/K.c), -0.5),
        ('08', "Perizentrum und Apozentrum sind r₋ und r₊", (lambda M,J: (0.5*(J/M/K.horizonte(J,M)[1].real)**2-K.G*M/K.horizonte(J,M)[1].real)/K.c**2)(10*K.M_sun,0.8*K.G*(10*K.M_sun)**2/K.c), -0.5),
        ('08', "ihre Exzentrizität e = √(1 − j²) = (r₊ − r₋)/r₀", (lambda M,J: (K.horizonte(J,M)[0]-K.horizonte(J,M)[1]).real/K.grenzradius(M)/math.sqrt(1-K.horizonte(J,M)[2]**2))(10*K.M_sun,0.3*K.G*(10*K.M_sun)**2/K.c), 1.0),
        ('08', "wo das Tempo der Flucht c erreicht", math.sqrt(2*K.G*K.M_sun/K.grenzradius(K.M_sun))/K.c, 1.0),
        ('08', "J/M = 3Δs·c", K.horizonte(3*K.faltung(K.M_sun)*K.c*K.M_sun, K.M_sun)[2], 1.0),
        ('08', "mit c durchlaufen", math.sqrt(K.G*K.M_sun/(3*K.faltung(K.M_sun)))/K.c, 1.0),
        ('08', "−½c²/j²", (lambda M,j: min(-K.G*M/r+0.5*(j*K.G*M/K.c)**2/r**2 for r in [3*K.faltung(M)*x/1000 for x in range(1,20000)])/K.c**2*j**2)(10*K.M_sun,1.5), -0.5),
        ('08', "Abbildung der Scheitel auf die Wendepunkte, r₀(M_irr)² = r₀·r₊", (lambda M,J: K.grenzradius(K.falte_gedreht(M,J)['M_irr'])**2/(K.grenzradius(M)*K.horizonte(J,M)[0].real))(10*K.M_sun,0.8*K.G*(10*K.M_sun)**2/K.c), 1.0),
        ('08', "a_J·r₀ = r₀(M_irr)·r₀(M_rot)", (lambda M,J: K.verwindung(J,M)*K.grenzradius(M)/(K.grenzradius(K.falte_gedreht(M,J)['M_irr'])*K.grenzradius(K.falte_gedreht(M,J)['M_rot'])))(10*K.M_sun,0.8*K.G*(10*K.M_sun)**2/K.c), 1.0),
        ('08', "j der Drehimpuls gegen den des Kreises gleicher Energie", (0.6*K.G*(10*K.M_sun)**2/K.c)/(10*K.M_sun)/math.sqrt(K.G*10*K.M_sun*3*K.faltung(10*K.M_sun)), 0.6),
        ('08', "j der Drehimpuls gegen den des Kreises gleicher Energie", (lambda d: 2*d['M_irr']*d['M_rot']/(d['M_irr']**2+d['M_rot']**2))(K.falte_gedreht(10*K.M_sun,0.6*K.G*(10*K.M_sun)**2/K.c)), 0.6),
        ('08', "e = 0,8", (lambda M,J: (K.horizonte(J,M)[0]-K.horizonte(J,M)[1]).real/K.grenzradius(M))(10*K.M_sun,0.6*K.G*(10*K.M_sun)**2/K.c), 0.8),
        ('08', "T = 0,89·T₀", K.falte_gedreht(10*K.M_sun,0.6*K.G*(10*K.M_sun)**2/K.c)['temperatur']/K.falte_temperatur(10*K.M_sun), 0.89),
        ('08', "S = 0,9·S₀", K.falte_gedreht(10*K.M_sun,0.6*K.G*(10*K.M_sun)**2/K.c)['entropie']/K.entropie_falte(10*K.M_sun), 0.9),
        ('08', "−GM/r − ½H_∞²r²", (lambda M: (lambda e: (-K.G*M/e['r_ganzes']-0.5*(K.c/(math.sqrt(3)*e['r_N']))**2*e['r_ganzes']**2)/K.c**2)(K.falte_endzustand(M)))(1e15*K.M_sun), -0.5),
        ('08', "−GM/r − ½H_∞²r²", (lambda M: (lambda e: (-K.G*M/e['r_falte']-0.5*(K.c/(math.sqrt(3)*e['r_N']))**2*e['r_falte']**2)/K.c**2)(K.falte_endzustand(M)))(1e20*K.M_sun), -0.5),
        ('03', "an der Schwelle 2,14-mal", K.verstaerkung(K.a0)**2 * (1 - K.schalter(K.a0) / 4), 2.14),
        ('03', "an der Schwelle 2,14-mal", (lambda x: (lambda q: sum((n + 1) * q**n for n in range(200000)))(math.exp(-math.sqrt(x))) * sum(1 - K.schalter(x * K.a0) / 2 * math.cos(2 * math.pi * k / 1000)**2 for k in range(1000)) / 1000)(1.0), 2.14),
        ('03', "bei a_N = 0,1·a₀ 10,7-mal", K.verstaerkung(0.1 * K.a0)**2 * (1 - K.schalter(0.1 * K.a0) / 4), 10.7),
        ('03', "bei a_N = 0,1·a₀ 10,7-mal", (lambda x: (lambda q: sum((n + 1) * q**n for n in range(200000)))(math.exp(-math.sqrt(x))) * sum(1 - K.schalter(x * K.a0) / 2 * math.cos(2 * math.pi * k / 1000)**2 for k in range(1000)) / 1000)(0.1), 10.7),
        ('03', "bei a_N = 0,1·a₀ 10,7-mal", K.reibung(1e5 * K.M_sun, 3e6, 1e-21, 1e3, 1e15, 1e20, 0.1 * K.a0, 1e15) / K.reibung(1e5 * K.M_sun, 3e6, 1e-21, 1e3, 1e15, 1e20), 10.7),
        ('03', "0 für schnelle Körper und 3/5 für langsame", (lambda X: 3 / (2 * X * X) - math.exp(-X * X) / (2 * X * X * (math.sqrt(math.pi) * math.erf(X) / (4 * X**3) - math.exp(-X * X) / (2 * X * X))))(30.0), 0),
        ('03', "0 für schnelle Körper und 3/5 für langsame", (lambda X: __import__('scipy.integrate', fromlist=['quad']).quad(lambda u: u**4 * math.exp(-X * X * u * u), 0, 1)[0] / __import__('scipy.integrate', fromlist=['quad']).quad(lambda u: u**2 * math.exp(-X * X * u * u), 0, 1)[0])(30.0), 0),
        ('03', "0 für schnelle Körper und 3/5 für langsame", (lambda X: 3 / (2 * X * X) - math.exp(-X * X) / (2 * X * X * (math.sqrt(math.pi) * math.erf(X) / (4 * X**3) - math.exp(-X * X) / (2 * X * X))))(0.01), 0.6),
        ('03', "0 für schnelle Körper und 3/5 für langsame", (lambda X: __import__('scipy.integrate', fromlist=['quad']).quad(lambda u: u**4 * math.exp(-X * X * u * u), 0, 1)[0] / __import__('scipy.integrate', fromlist=['quad']).quad(lambda u: u**2 * math.exp(-X * X * u * u), 0, 1)[0])(0.01), 0.6),
        ('03', "ν² = Σ(n + 1)·qⁿ", (lambda q: sum((n + 1) * q**n for n in range(100000)) / (1 / (1 - q))**2)(math.exp(-0.5)), 1.0),
        ('03', "in einer Kugel ν(1 − χ/2)·ρ_b + χΩ²/(4πG)", (lambda M, a, r: (lambda g: (((r * 1.0001)**2 * g(r * 1.0001) - (r * 0.9999)**2 * g(r * 0.9999)) / (2e-4 * r) / (4 * math.pi * K.G * r * r)) / (K.verstaerkung(K.G * M * r / (r * r + a * a)**1.5) * (1 - K.schalter(K.G * M * r / (r * r + a * a)**1.5) / 2) * 3 * M / (4 * math.pi * a**3) * (1 + r * r / (a * a))**-2.5 + K.schalter(K.G * M * r / (r * r + a * a)**1.5) * g(r) / r / (4 * math.pi * K.G)))(lambda x: K.verstaerkung(K.G * M * x / (x * x + a * a)**1.5) * K.G * M * x / (x * x + a * a)**1.5))(1e9 * K.M_sun, 2 * K.kpc, 3 * K.kpc), 1.0),
        ("08", "weiter die Faltung der Linse selbst", (lambda M, D, y: (lambda r: r[0] / (K.grenzradius(M) / K.c * (r[2] * math.sqrt(r[2]**2 + 4) + 2 * math.log((math.sqrt(r[2]**2 + 4) + r[2]) / (math.sqrt(r[2]**2 + 4) - r[2])))))(K.laufzeit_regel(M, D, y)))(1e11 * K.M_sun, K.einsteinring(0.5, 2.0)[0], 0.3), 1),
        ("03", "ϑ = 3,49", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda x: 8 * (math.cosh(x) - 1)**2 - 9 * math.sinh(x) * (math.sinh(x) - x), 1.0, 10.0), 3.49),
        ("03", "linear δ_v = −2,72", K.leerraum_rand()[0], -2.72),
        ("03", "linear δ_v = −2,72", (lambda si, ev: (lambda s: [-1e-4 * s.t_events[0][0]**(2 / 3), 4.5 * s.t_events[0][0]**2 / s.y_events[0][0][0]**3, s.y_events[0][0][1] / s.y_events[0][0][0] * 1.5 * s.t_events[0][0]])(si(lambda t, y: [y[1], -1 / y[0]**2, y[3], -(1 + 1e-5) / y[2]**2], (1.0, 1e12), [4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3) * (1 + 1e-4 / 3), 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3) * (1 + 1e-4 / (3 * (1 + 1e-5)))], events=ev, rtol=1e-12, atol=1e-14, method='DOP853')))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, type('Rand', (), {'__call__': lambda self, t, y: y[2] - y[0], 'terminal': True, 'direction': -1})())[0], -2.72),
        ("03", "1 + δ = 0,20", K.leerraum_rand()[1], 0.2),
        ("03", "1 + δ = 0,20", (lambda si, ev: (lambda s: [-1e-4 * s.t_events[0][0]**(2 / 3), 4.5 * s.t_events[0][0]**2 / s.y_events[0][0][0]**3, s.y_events[0][0][1] / s.y_events[0][0][0] * 1.5 * s.t_events[0][0]])(si(lambda t, y: [y[1], -1 / y[0]**2, y[3], -(1 + 1e-5) / y[2]**2], (1.0, 1e12), [4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3) * (1 + 1e-4 / 3), 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3) * (1 + 1e-4 / (3 * (1 + 1e-5)))], events=ev, rtol=1e-12, atol=1e-14, method='DOP853')))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, type('Rand', (), {'__call__': lambda self, t, y: y[2] - y[0], 'terminal': True, 'direction': -1})())[1], 0.2),
        ("03", "1,70-fachen", K.leerraum_rand()[2], 1.7),
        ("03", "1,70-fachen", (lambda si, ev: (lambda s: [-1e-4 * s.t_events[0][0]**(2 / 3), 4.5 * s.t_events[0][0]**2 / s.y_events[0][0][0]**3, s.y_events[0][0][1] / s.y_events[0][0][0] * 1.5 * s.t_events[0][0]])(si(lambda t, y: [y[1], -1 / y[0]**2, y[3], -(1 + 1e-5) / y[2]**2], (1.0, 1e12), [4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3) * (1 + 1e-4 / 3), 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3) * (1 + 1e-4 / (3 * (1 + 1e-5)))], events=ev, rtol=1e-12, atol=1e-14, method='DOP853')))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, type('Rand', (), {'__call__': lambda self, t, y: y[2] - y[0], 'terminal': True, 'direction': -1})())[1]**(-1 / 3), 1.7),
        ("03", "H_v = (4/3)·H", 3 * K.leerraum_rand()[3], 4),
        ("03", "H_v = (4/3)·H", 3 * (lambda si, ev: (lambda s: [-1e-4 * s.t_events[0][0]**(2 / 3), 4.5 * s.t_events[0][0]**2 / s.y_events[0][0][0]**3, s.y_events[0][0][1] / s.y_events[0][0][0] * 1.5 * s.t_events[0][0]])(si(lambda t, y: [y[1], -1 / y[0]**2, y[3], -(1 + 1e-5) / y[2]**2], (1.0, 1e12), [4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * (1 / (1 - 1e-4))**(1 / 3) * (1 + 1e-4 / 3), 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3), 2 / 3 * 4.5**(1 / 3) * ((1 + 1e-5) / (1 - 1e-4 / (1 + 1e-5)))**(1 / 3) * (1 + 1e-4 / (3 * (1 + 1e-5)))], events=ev, rtol=1e-12, atol=1e-14, method='DOP853')))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, type('Rand', (), {'__call__': lambda self, t, y: y[2] - y[0], 'terminal': True, 'direction': -1})())[2], 4),
        ("03", "bei denselben δ_v und 1 + δ", K.leerraum_rand(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc))[0], -2.72),
        ("03", "bei denselben δ_v und 1 + δ", (lambda qu, bq, si: (lambda tk, ti: (lambda ai, Hi: (lambda sch, tvr: (lambda kreuz: (lambda di: (lambda k, s: [di * s.y[0][-1] / ai, (math.sinh(1.5 * k[0])**(2 / 3) / k[1])**3, k[2] / k[1] * math.tanh(1.5 * k[0]), s.y[1][-1] / s.y[0][-1] * math.tanh(1.5 * tk)])(kreuz(di), si(lambda t, y: [y[1], -2 / math.tanh(1.5 * t) * y[1] + 1.5 / math.sinh(1.5 * t)**2 * y[0]], (ti, tk), [ai, Hi * ai], rtol=1e-12, atol=1e-20, method='DOP853')))(bq(lambda d: kreuz(d)[0] - tk, -2.5e-3, -1.4e-3, xtol=1e-13)))(lambda d: (lambda A, B: (lambda rc: (tvr(A, rc), rc, math.sqrt(2 * A[1] + 2 * A[2] / rc + rc * rc)))(bq(lambda r: tvr(A, r) - tvr(B, r), B[0] * (1 + 1e-9), 1e3, xtol=1e-13, rtol=1e-14)))(sch(1.0, d), sch(1 + 1e-6, d / (1 + 1e-6)))))(lambda M, d: (lambda ri: (ri, 0.5 * (Hi * ri * (1 - d / 3))**2 - 0.5 * M / ri - 0.5 * ri**2, 0.5 * M))(ai * (M / (1 + d))**(1 / 3)), lambda sh, r: ti + qu(lambda x: 1 / math.sqrt(2 * sh[1] + 2 * sh[2] / x + x * x), sh[0], r, epsabs=0, epsrel=1e-11, limit=800)[0]))(math.sinh(1.5 * ti)**(2 / 3), 1 / math.tanh(1.5 * ti)))(2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)), 1e-5))(__import__('scipy.integrate', fromlist=['quad']).quad, __import__('scipy.optimize', fromlist=['brentq']).brentq, __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp)[0], -2.72),
        ("03", "bei denselben δ_v und 1 + δ", K.leerraum_rand(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc))[1], 0.2),
        ("03", "bei denselben δ_v und 1 + δ", (lambda qu, bq, si: (lambda tk, ti: (lambda ai, Hi: (lambda sch, tvr: (lambda kreuz: (lambda di: (lambda k, s: [di * s.y[0][-1] / ai, (math.sinh(1.5 * k[0])**(2 / 3) / k[1])**3, k[2] / k[1] * math.tanh(1.5 * k[0]), s.y[1][-1] / s.y[0][-1] * math.tanh(1.5 * tk)])(kreuz(di), si(lambda t, y: [y[1], -2 / math.tanh(1.5 * t) * y[1] + 1.5 / math.sinh(1.5 * t)**2 * y[0]], (ti, tk), [ai, Hi * ai], rtol=1e-12, atol=1e-20, method='DOP853')))(bq(lambda d: kreuz(d)[0] - tk, -2.5e-3, -1.4e-3, xtol=1e-13)))(lambda d: (lambda A, B: (lambda rc: (tvr(A, rc), rc, math.sqrt(2 * A[1] + 2 * A[2] / rc + rc * rc)))(bq(lambda r: tvr(A, r) - tvr(B, r), B[0] * (1 + 1e-9), 1e3, xtol=1e-13, rtol=1e-14)))(sch(1.0, d), sch(1 + 1e-6, d / (1 + 1e-6)))))(lambda M, d: (lambda ri: (ri, 0.5 * (Hi * ri * (1 - d / 3))**2 - 0.5 * M / ri - 0.5 * ri**2, 0.5 * M))(ai * (M / (1 + d))**(1 / 3)), lambda sh, r: ti + qu(lambda x: 1 / math.sqrt(2 * sh[1] + 2 * sh[2] / x + x * x), sh[0], r, epsabs=0, epsrel=1e-11, limit=800)[0]))(math.sinh(1.5 * ti)**(2 / 3), 1 / math.tanh(1.5 * ti)))(2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)), 1e-5))(__import__('scipy.integrate', fromlist=['quad']).quad, __import__('scipy.optimize', fromlist=['brentq']).brentq, __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp)[1], 0.2),
        ("03", "auf 0,1 %", (lambda w: 100 * abs(w[3] / (1 + w[4] / 3) - 1))(K.leerraum_rand(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc))), 0.1),
        ("03", "auf 0,1 %", (lambda w: 100 * abs(w[2] / (1 + w[3] / 3) - 1))((lambda qu, bq, si: (lambda tk, ti: (lambda ai, Hi: (lambda sch, tvr: (lambda kreuz: (lambda di: (lambda k, s: [di * s.y[0][-1] / ai, (math.sinh(1.5 * k[0])**(2 / 3) / k[1])**3, k[2] / k[1] * math.tanh(1.5 * k[0]), s.y[1][-1] / s.y[0][-1] * math.tanh(1.5 * tk)])(kreuz(di), si(lambda t, y: [y[1], -2 / math.tanh(1.5 * t) * y[1] + 1.5 / math.sinh(1.5 * t)**2 * y[0]], (ti, tk), [ai, Hi * ai], rtol=1e-12, atol=1e-20, method='DOP853')))(bq(lambda d: kreuz(d)[0] - tk, -2.5e-3, -1.4e-3, xtol=1e-13)))(lambda d: (lambda A, B: (lambda rc: (tvr(A, rc), rc, math.sqrt(2 * A[1] + 2 * A[2] / rc + rc * rc)))(bq(lambda r: tvr(A, r) - tvr(B, r), B[0] * (1 + 1e-9), 1e3, xtol=1e-13, rtol=1e-14)))(sch(1.0, d), sch(1 + 1e-6, d / (1 + 1e-6)))))(lambda M, d: (lambda ri: (ri, 0.5 * (Hi * ri * (1 - d / 3))**2 - 0.5 * M / ri - 0.5 * ri**2, 0.5 * M))(ai * (M / (1 + d))**(1 / 3)), lambda sh, r: ti + qu(lambda x: 1 / math.sqrt(2 * sh[1] + 2 * sh[2] / x + x * x), sh[0], r, epsabs=0, epsrel=1e-11, limit=800)[0]))(math.sinh(1.5 * ti)**(2 / 3), 1 / math.tanh(1.5 * ti)))(2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)), 1e-5))(__import__('scipy.integrate', fromlist=['quad']).quad, __import__('scipy.optimize', fromlist=['brentq']).brentq, __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp)), 0.1),
        ("03", "umkehrt (Θ = 1,020)", K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc), 1.02),
        ("03", "11,8-fache", K.umkehr_heute()[1], 11.8),
        ("03", "11,8-fache", math.sinh(1.5 * 2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)))**2 / (lambda si, bq, t0, ev: (lambda tu: tu(bq(lambda d: tu(d)[0] - t0, 1e-5, 1e-2, xtol=1e-14))[1])(lambda d: (lambda s: (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (5 * t0, None))(si(lambda t, y: [y[1], -0.5 / y[0]**2 + y[0]], (1e-5, 5 * t0), [math.sinh(1.5e-5)**(2 / 3) * (1 / (1 + d))**(1 / 3), math.sinh(1.5e-5)**(2 / 3) * (1 / (1 + d))**(1 / 3) / math.tanh(1.5e-5) * (1 - d / 3)], events=ev, rtol=1e-11, atol=1e-16, method='DOP853'))))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, __import__('scipy.optimize', fromlist=['brentq']).brentq, 2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)), type('Umkehr', (), {'__call__': lambda self, t, y: y[1], 'terminal': True, 'direction': -1})())**3, 11.8),
        ("03", "(3π/4)² = 5,55", (3 * math.pi / 4)**2, 5.55),
        ("03", "(3π/4)² = 5,55", K.kollaps_kugel()[2], 5.55),
        ("03", "0,72·R_max", K.umkehr_heute()[0], 0.72),
        ("03", "0,72·R_max", 2**(1 / 3) * (lambda si, bq, t0, ev: (lambda tu: tu(bq(lambda d: tu(d)[0] - t0, 1e-5, 1e-2, xtol=1e-14))[1])(lambda d: (lambda s: (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (5 * t0, None))(si(lambda t, y: [y[1], -0.5 / y[0]**2 + y[0]], (1e-5, 5 * t0), [math.sinh(1.5e-5)**(2 / 3) * (1 / (1 + d))**(1 / 3), math.sinh(1.5e-5)**(2 / 3) * (1 / (1 + d))**(1 / 3) / math.tanh(1.5e-5) * (1 - d / 3)], events=ev, rtol=1e-11, atol=1e-16, method='DOP853'))))(__import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp, __import__('scipy.optimize', fromlist=['brentq']).brentq, 2 / 3 * math.atanh(K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc) * math.pi / math.sqrt(15)), type('Umkehr', (), {'__call__': lambda self, t, y: y[1], 'terminal': True, 'direction': -1})()), 0.72),
        ("03", "bei 2R_d das 2,08-fache", (lambda gR,u: (lambda F: (F((gR**2+u**2)**0.25)-F(gR**0.5))/(u**2/4))(lambda s: s**4/4+sum(6/k**4*(1-math.exp(-k*s)*(1+k*s+(k*s)**2/2+(k*s)**3/6)) for k in range(1,20000))))((-math.log(1-1/K.toomre_schwelle(2.0)['nu']))**2, 2*math.exp(-2)), 2.08),
        ("03", "0,888·σ²/√(2πGΣa₀)", 1.5*sum(0.5**(1.5*n+1)/(1.5*n+1) for n in range(400)), 0.888),
        ("03", "0,549·σ²/(πGΣ)", math.atanh(0.5), 0.549),
        ("03", "= ρ_s·c² = 1,74·10⁻¹⁰ Pa", K.dichte_vakuumteil()*K.c**2*1e10, 1.74),
        ("03", "höchstens π³a₀²/(30G)", K.a0**2/(2*math.pi*K.G)*K.strahlung_je_sinn(1,3)*1e10, 1.74),
        ("03", "bei R_d 2,16", K.toomre_schwelle(1.0)["arme"], 2.16),
        ("03", "bei R_d 2,16", K.toomre_schwelle(1.0, numerisch=True)["arme"], 2.16),
        ("03", "bei 2R_d 3,45", K.toomre_schwelle(2.0)["arme"], 3.45),
        ("03", "bei 2R_d 3,45", K.toomre_schwelle(2.0, numerisch=True)["arme"], 3.45),
        ("03", "(dort d_eff = 1,58)", K.toomre_schwelle(2.0, numerisch=True)["d_eff"], 1.58),
        ("03", "gäbe dort 2,98", K.toomre_schwelle(2.0)["arme_newton"], 2.98),
        ("03", "bei 2R_d 7,45", K.toomre_schwelle(2.0)["arme_halo"], 7.45),
        ("03", "bei 2R_d 7,45", K.toomre_schwelle(2.0, numerisch=True)["arme_halo"], 7.45),
        ("03", "x₀ = 3,830", K.sunjajew_seldowitsch()[0], 3.83),
        ("03", "bei 217,5 GHz", K.sunjajew_seldowitsch()[1] / 1e9, 217.5),
        ("03", "um den Anteil 4y", K.sunjajew_seldowitsch()[2], 4),
        ("03", "die Stufe d = 2, bleibt", K.sunjajew_seldowitsch()[3], 0),
        ("03", "ΔT/T = −2y", -K.sunjajew_seldowitsch()[4], 2),
        ("03", "∫s^d·⟨Δn²⟩ds = d·Γ(d)ζ(d)", float(max(abs(__import__('scipy.integrate', fromlist=['quad']).quad(lambda s: s**d * math.exp(-s) / math.expm1(-s)**2, 1e-12, 100, limit=400)[0] / (d * math.gamma(d) * __import__('scipy.special', fromlist=['zeta']).zeta(d)) - 1) for d in (2, 3, 4, 5)) < 1e-9), 1),
        ("03", "s·coth(s/2) = s·(2ν − 1) = d", float(max(abs(s / math.tanh(s / 2) - s * (2 * K.verstaerkung(s * s * K.a0) - 1)) for s in (0.3, 1.0, 3.83)) < 1e-9), 1),
        ("03", "= 1,09·√(GM_b·a₀)", (lambda M: K.verstaerkung(K.G * M / K.naht_radius(M)**2) * K.G * M / K.naht_radius(M) / K.v_flach(M)**2)(1e12 * K.M_sun), 1.09),
        ("03", "heute 2,3·10³ Gyr", K.abkuehlung_licht() / K.Gyr / 1e3, 2.3),
        ("03", "noch 1,4·10² Gyr", K.abkuehlung_licht(1.0) / K.Gyr / 1e2, 1.4),
        ("03", "heute kreuzt (Θ = 1,020)", K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc), 1.02),
        ('08', "mit der Determinante 1 und der Spur 2·cos(2π·Re τ_M)", (lambda m: m[0][0]*m[1][1] - m[0][1]*m[1][0])(K.umlauf_abbildung(K.takte_kreisbahn(0.01)['radial'])['matrix']), 1),
        ('08', "mit der Determinante 1 und der Spur 2·cos(2π·Re τ_M)", K.umlauf_abbildung(K.takte_kreisbahn(0.01)['radial'])['spur'] / math.cos(2*math.pi*math.sqrt(K.takte_kreisbahn(0.01)['radial'])), 2),
        ('08', "λ² − Spur·λ + 1 = 0", abs(__import__('numpy').prod(__import__('numpy').linalg.eigvals(__import__('numpy').array(K.umlauf_abbildung(K.takte_kreisbahn(0.3)['radial'])['matrix'])))), 1),
        ('08', "liegen die Eigenwerte auf dem Einheitskreis", abs(K.umlauf_abbildung(K.takte_kreisbahn(0.1)['radial'])['eigen'][0]), 1),
        ('08', "Bei 18Δs fallen sie in eins zusammen", 3 / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.takte_kreisbahn(d)['radial'], 0.1, 0.3), 18),
        ('08', "(δu, δu′) → (δu + 2π·δu′, δu′)", (lambda k: math.sin(2*math.pi*k) / k / math.pi)(math.sqrt(K.takte_kreisbahn(1/6 - 1e-12)['radial'])), 2),
        ('08', "e^{±2π·√(6δ − 1)}", K.umlauf_abbildung(K.takte_kreisbahn(0.25)['radial'])['eigen'][0].real / math.exp(2*math.pi*math.sqrt(6*0.25 - 1)), 1),
        ('08', "an der Photonensphäre mit der Spur 2·cosh 2π = 535,5", K.umlauf_abbildung(K.takte_kreisbahn(1/3)['radial'])['spur'], 535.5),
        ('08', "Bei 24Δs = 8d ist Re τ_M = ½", 3 / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.takte_kreisbahn(d)['radial'] - 0.25, 0.05, 0.16), 24),
        ('08', "Bei 24Δs = 8d ist Re τ_M = ½", 1 / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.takte_kreisbahn(d)['radial'] - 0.25, 0.05, 0.16), 8),
        ('08', "Bei 24Δs = 8d ist Re τ_M = ½", math.sqrt(K.takte_kreisbahn(1/8)['radial']), 0.5),
        ('08', "kehrt nach einem Umlauf mit umgekehrtem Vorzeichen zurück", K.umlauf_abbildung(K.takte_kreisbahn(1/8)['radial'])['matrix'][1][1], -1),
        ('08', "schließt nach zwei Umläufen", 1 / math.sqrt(K.takte_kreisbahn(1/8)['radial']), 2),
        ('08', "Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität", K.umlaufphase(1e4*K.a0), 1),
        ('08', "Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität", (lambda g, r, h: math.sqrt(3 + (math.log(g(r*(1+h))) - math.log(g(r*(1-h)))) / (math.log(1+h) - math.log(1-h))))(lambda r: r, 2.0, 1e-6), 2),
        ('08', "Newton und der gehaltene Körper, Re τ_M = 1 und 2, haben die Identität", round(1e9 * max(abs(K.umlauf_abbildung(k2)['matrix'][i][j] - (i == j)) for k2 in (1.0, 4.0) for i in (0, 1) for j in (0, 1))), 0),
        ('08', "Die Regel allein senkt die Spur von 2 höchstens bis 2·cos(2π√2) = −1,72", K.umlauf_abbildung(K.umlaufphase(1e4*K.a0)**2)['spur'], 2),
        ('08', "Die Regel allein senkt die Spur von 2 höchstens bis 2·cos(2π√2) = −1,72", K.umlauf_abbildung(K.umlaufphase(1e-14*K.a0)**2)['spur'], -1.72),
        ('08', "bringt um eine Punktmasse keine Kreisbahn aus der Drehung", float(all(K.umlauf_abbildung(K.umlaufphase(x*K.a0)**2)['art'] == 'Drehung' and abs(K.umlauf_abbildung(K.umlaufphase(x*K.a0)**2)['spur']) < 2 for x in (1e-12, 1e-4, 1.0, 10.0))), 1),
        ('08', "An der letzten stabilen Kreisbahn um einen gebundenen Bereich liegt die zweite Scherung", round(1e9 * abs((lambda M: (lambda R, H: (2*(K.v_flach(M)/R)**2 - 4*H**2) / ((K.v_flach(M)/R)**2 - H**2))(K.letzte_stabile_bahn(M)[0], K.groesste_umkehr(M)[1]))(6e10*K.M_sun))), 0),
        ('08', "An der letzten stabilen Kreisbahn um einen gebundenen Bereich liegt die zweite Scherung", round(1e9 * abs((lambda M: (lambda R, H: (K.G*M/R**3 - 4*H**2) / (K.G*M/R**3 - H**2))(K.letzte_stabile_bahn(M, regel=False)[0], K.groesste_umkehr(M)[1]))(1e15*K.M_sun))), 0),
        ('08', "und bis R_max dehnt die Abbildung", float((lambda M: (lambda Rm, H: K.umlauf_abbildung((2*(K.v_flach(M)/(0.75*Rm))**2 - 4*H**2) / ((K.v_flach(M)/(0.75*Rm))**2 - H**2))['spur'])(*K.groesste_umkehr(M)))(6e10*K.M_sun) > 2), 1),
        ('08', "mit der Identität bei R_max/√3", (lambda M: (lambda Rm, H: (Rm / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: (2*(K.v_flach(M)/R)**2 - 4*H**2) / ((K.v_flach(M)/R)**2 - H**2) - 1, 0.3*Rm, 0.65*Rm))**2)(*K.groesste_umkehr(M)))(6e10*K.M_sun), 3),
        ('08', "der halben Drehung bei √(7/15)·R_max", (lambda M: (lambda Rm, H: 15 * (__import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: (2*(K.v_flach(M)/R)**2 - 4*H**2) / ((K.v_flach(M)/R)**2 - H**2) - 0.25, 0.6*Rm, 0.7*Rm) / Rm)**2)(*K.groesste_umkehr(M)))(6e10*K.M_sun), 7),
        ('08', "liegt die halbe Drehung bei R_max/∛5", (lambda M: (lambda H: (lambda Rm: (Rm / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: (K.G*M/R**3 - 4*H**2) / (K.G*M/R**3 - H**2) - 0.25, 0.5*Rm, 0.62*Rm))**3)((K.G*M/H**2)**(1/3)))(K.groesste_umkehr(M)[1]))(1e15*K.M_sun), 5),
        ('01', "hat ihre Abbildung die Determinante 1", (lambda ab: ab[0][0]*ab[1][1] - ab[0][1]*ab[1][0])((lambda Kk, s: (lambda f: (f([1, 0]), f([0, 1])))(lambda y0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda t, y: [y[1], -Kk*y[0]], (0, s), y0, rtol=1e-12, atol=1e-14).y[:, -1]))(0.7, 2.3)), 1),
        ('01', "und die Spur 2·cos(√K·λ)", (lambda ab: (ab[0][0] + ab[1][1])/(2*math.cos(math.sqrt(0.7)*2.3)))((lambda Kk, s: (lambda f: (f([1, 0]), f([0, 1])))(lambda y0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda t, y: [y[1], -Kk*y[0]], (0, s), y0, rtol=1e-12, atol=1e-14).y[:, -1]))(0.7, 2.3)), 1),
        ('01', "für K < 0 ist sie 2·cosh(√−K·λ)", (lambda ab: (ab[0][0] + ab[1][1])/(2*math.cosh(math.sqrt(0.7)*2.3)))((lambda Kk, s: (lambda f: (f([1, 0]), f([0, 1])))(lambda y0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda t, y: [y[1], -Kk*y[0]], (0, s), y0, rtol=1e-12, atol=1e-14).y[:, -1]))(-0.7, 2.3)), 1),
        ('01', "für K = 0 schert sie (Spur 2)", (lambda ab: ab[0][0] + ab[1][1])((lambda Kk, s: (lambda f: (f([1, 0]), f([0, 1])))(lambda y0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda t, y: [y[1], -Kk*y[0]], (0, s), y0, rtol=1e-12, atol=1e-14).y[:, -1]))(0.0, 2.3)), 2),
        ('01', "q ist am größten bei |k| = κ/σ, dort gleich 1/Q_T", (lambda ka, si, A: max(2*A*k/(ka**2 + si**2*k**2) for k in [i*1e-4 for i in range(1, 200000)])*ka*si/A)(1.3, 0.7, 0.5), 1),
        ('01', "q ist am größten bei |k| = κ/σ, dort gleich 1/Q_T", (lambda ka, si, A: max(range(1, 200000), key=lambda i: 2*A*i*1e-4/(ka**2 + si**2*(i*1e-4)**2))*1e-4*si/ka)(1.3, 0.7, 0.5), 1),
        ('01', "die Jeans-Länge, Q_T = 1", (lambda ka, si: ka*si/__import__('scipy.optimize', fromlist=['brentq']).brentq(lambda A: min(ka**2 - 2*A*k + si**2*k*k for k in [i*1e-3 for i in range(1, 20000)]), 0.1, 3.0, xtol=1e-12))(1.3, 0.7), 1),
        ('01', "Ihre Antwort Σqⁿ = 1/(1 − q) ist ω_frei²/ω²", (lambda T, n, m, l: ((K.k_B*T/m)/l**2/((K.k_B*T/m)/l**2 - 4*math.pi*K.G*n*m))*(1 - l**2*4*math.pi*K.G*n*m**2/(K.k_B*T)))(10.0, 1e10, K.m_p, 0.6*math.sqrt(K.k_B*10.0/(4*math.pi*K.G*1e10*K.m_p**2))), 1),
        ('01', "Ihre Antwort Σqⁿ = 1/(1 − q) ist ω_frei²/ω²", (lambda T, n, m, l: ((K.k_B*T/m)/l**2/((K.k_B*T/m)/l**2 + n*K.e**2/(K.eps0*m)))*(1 + l**2*n*K.e**2/(K.eps0*K.k_B*T)))(1e4, 1e6, K.m_e, 3.0*math.sqrt(K.eps0*K.k_B*1e4/(1e6*K.e**2))), 1),
        ('01', "(dz/dλ)² − (E/2)·z² = GM/2", (lambda r, E: r*2*(E + 1/r)/4 - E/2*r)(0.37, -0.3), 0.5),
        ('01', "K = −E/2", (lambda y: y[0]/(math.sqrt(1/0.3)*math.sin(math.sqrt(0.15)*1.3))**2)((lambda E, s0, s1: (lambda k, A: (lambda z0, zp0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda s, y: [y[0]*y[1], -1/y[0], y[0]], (s0, s1), [z0*z0, 2*zp0/z0, 0], rtol=1e-13, atol=1e-15, method='DOP853').y[:, -1])(A*math.sin(k*s0), A*k*math.cos(k*s0)))(math.sqrt(-E/2), math.sqrt(1/(-E))))(-0.3, 0.2, 1.3)), 1),
        ('01', "der rollende Kreis mit θ = 2√K·λ", (lambda y, th: y[2]/(math.sqrt((1/0.3)**3/8)*((th(1.3) - math.sin(th(1.3))) - (th(0.2) - math.sin(th(0.2))))))((lambda E, s0, s1: (lambda k, A: (lambda z0, zp0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda s, y: [y[0]*y[1], -1/y[0], y[0]], (s0, s1), [z0*z0, 2*zp0/z0, 0], rtol=1e-13, atol=1e-15, method='DOP853').y[:, -1])(A*math.sin(k*s0), A*k*math.cos(k*s0)))(math.sqrt(-E/2), math.sqrt(1/(-E))))(-0.3, 0.2, 1.3), lambda s: 2*math.sqrt(0.15)*s), 1),
        ('01', "steht z mit umgekehrtem Vorzeichen", (lambda k, s0: __import__('scipy.integrate', fromlist=['solve_ivp']).solve_ivp(lambda s, y: [y[1], -k*k*y[0]], (s0, s0 + math.pi/k), [math.sin(k*s0), k*math.cos(k*s0)], rtol=1e-12, atol=1e-14).y[0, -1]/math.sin(k*s0))(math.sqrt(0.15), 0.4), -1),
        ('01', "r³ = (9/2)·GM·t²", 1/__import__('scipy.integrate', fromlist=['quad']).quad(lambda x: 1/math.sqrt(2/x), 0, 1)[0]**2, 4.5),
        ('01', "r³ = (9/2)·GM·t²", (lambda a2, s: (a2*s*s)**3/(a2*s**3/3)**2)(0.5, 1.3), 4.5),
        ('03', "2b₁·ρ_s ist ihre Gezeitendichte 3𝒯/(4πG)", K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2)/(2*K.dichte_dunkle_energie()), 1),
        ('03', "R_max³ = GM_dyn/H_∞²", (lambda M: (K.G*(K.groesste_umkehr(M)[0]*K.v_flach(M)**2/K.G)/K.groesste_umkehr(M)[1]**2)**(1/3)/K.groesste_umkehr(M)[0])(6e10*K.M_sun), 1),
        ('03', "Θ = Θ_∞/√3", K.gezeit_ganzes(math.sqrt(15)/(math.pi*math.sqrt(3)))[1] - K.uhr_theta(math.sqrt(15)/(math.pi*math.sqrt(3)))[2], 1),
        ('03', "Θ = 1,020", K.theta_(K.a0, K.atmen_heute()[0]*1e3/K.Mpc), 1.02),
        ('03', "1,63·R_max", K.gezeit_ganzes()[2], 1.63),
        ('03', "4,34-fache", K.gezeit_ganzes()[1], 4.34),
        ('02', "A² = 3c²/(8πGρ) = R³/r₀", K.kappe(K.M_earth, K.R_ERDE)["A"] / math.sqrt(K.R_ERDE**3 / K.grenzradius(K.M_earth)), 1),
        ('02', "für die Erde 1,14 AE", K.kappe(K.M_earth, K.R_ERDE)["A"] / K.AU, 1.14),
        ('02', "sin²ψ_R = 2δ", math.sin(K.kappe(K.M_sun, K.R_SONNE)["chi_R"])**2 / (2 * 3 * K.faltung(K.M_sun) / K.R_SONNE), 1),
        ('02', "an der Hülle ist das der Takt", K.kappe(1.4 * K.M_sun, 12e3)["takt"] / K.takt(1.4 * K.M_sun, 12e3), 1),
        ('02', "(3/2)·cos ψ_R − ½·cos ψ", K.kappe(K.M_sun, 1.2 * K.grenzradius(K.M_sun))["uhr"] / K.mitte_gehalten(K.M_sun, 1.2 * K.grenzradius(K.M_sun), numerisch=True)[0], 1),
        ('02', "In der Mitte steht sie bei cos ψ_R = ⅓", K.kappe(K.M_sun, __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda R: K.kappe(K.M_sun, R)["uhr"], 1.0001 * K.grenzradius(K.M_sun), 3 * K.grenzradius(K.M_sun)))["takt"], 0.333333),
        ('02', "70,5° vom Pol", math.degrees(K.kappe(K.M_sun, __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda R: K.kappe(K.M_sun, R)["uhr"], 1.0001 * K.grenzradius(K.M_sun), 3 * K.grenzradius(K.M_sun)))["chi_R"]), 70.5),
        ('02', "bis 36° vom Pol", math.degrees(K.kappe(1.4 * K.M_sun, 12e3)["chi_R"]), 36),
        ('02', "R³/(6A²) = Δs", K.kappe(K.M_earth, K.R_ERDE)["ueberschuss"] / K.faltung(K.M_earth), 1),
        ('02', "Bei 90° läge die Hülle am Äquator, R = A = r₀", math.degrees(K.kappe(K.M_sun, K.grenzradius(K.M_sun))["chi_R"]), 90),
        ('02', "c²/A² = (8πG/3)·ρ + H_∞²", K.kappe(1e15 * K.M_sun, 1e22, H=K.groesste_umkehr(K.M_sun)[1])["takt"] / math.sqrt(1 - K.grenzradius(1e15 * K.M_sun) / 1e22 - (K.groesste_umkehr(K.M_sun)[1] * 1e22 / K.c)**2), 1),
        ('02', "die halbe Kugel vom Radius c/H_∞", K.kappe(0.0, 1e22, H=K.groesste_umkehr(K.M_sun)[1])["A"] / K.falte_endzustand(1e-30 * K.M_sun)["r_ganzes"], 1),
        ('02', "ihr Äquator der Horizont des Ganzen", math.degrees(K.kappe(0.0, K.falte_endzustand(1e-30 * K.M_sun)["r_ganzes"], H=K.groesste_umkehr(K.M_sun)[1])["chi_R"]), 90),
        ('02', "im Augenblick seiner Umkehr", K.kappe(1e15 * K.M_sun, 1e22)["A"]**2 * 2 * K.G * 1e15 * K.M_sun / (K.c**2 * 1e22**3), 1),
        ('08', "= (Θ_∞/2π)²·ℓ₀ = 1,06 Gpc", K.a0/K.groesste_umkehr(K.M_sun)[1]**2/K.Gpc, 1.06),
        ('08', "= 8,5·10²⁰ M☉, x = 0,95", K.a0**3/(K.G*K.groesste_umkehr(K.M_sun)[1]**4)/K.M_sun/1e20, 8.5),
        ('08', "= 8,5·10²⁰ M☉, x = 0,95", K.leiter_ort(K.a0**3/(K.G*K.groesste_umkehr(K.M_sun)[1]**4)), 0.95),
        ('08', "bei v_flach = a₀/H_∞ = 0,196·c", K.a0/K.groesste_umkehr(K.M_sun)[1]/K.c, 0.196),
        ('08', "(6·10¹⁰ M☉) nach der Regel 345", K.groesste_umkehr(6e10*K.M_sun)[0]/K.schwellenradius(6e10*K.M_sun), 345),
        ("08", "·τ_D = 0,45 Mrd. Jahre", K.wachsen_grenze(1 - K.takte_kreisbahn(1/6)['energie'])['t_edd']/K.Gyr, 0.45),
        ("08", "ungedreht 3/(2√2) − 1 = 0,061", 1/K.takte_kreisbahn(1/6)['energie'] - 1, 0.061),
        ("08", "0,061 mit t_e = 27 Mio. Jahre", K.wachsen_grenze(1 - K.takte_kreisbahn(1/6)['energie'])['t_e']/K.Gyr*1e3, 27),
        ("08", "(j = 1) √3 − 1 = 0,73", (1 - 1/math.sqrt(3))/(1/math.sqrt(3)), 0.73),
        ("08", "0,73 mit t_e = 330 Mio. Jahre", K.wachsen_grenze(1 - 1/math.sqrt(3))['t_e']/K.Gyr*1e3, 330),
        ("08", "Θ = 0,1, nach 0,95 Mrd. Jahren", K.uhr_theta(0.1)[0]/K.Gyr, 0.95),
        ("08", "um den Faktor e³⁵", K.wachsen_grenze(1 - K.takte_kreisbahn(1/6)['energie'], 0.1)['stufen'], 35),
        ("08", "bei j = 1 um den Faktor 18", math.exp(K.wachsen_grenze(1 - 1/math.sqrt(3), 0.1)['stufen']), 18),
        ("05", "Bis zum Fall vergeht t = E/((l + 2)·P)", K.zeit_abgabe(1, (lambda r_e: K.m_e*K.c**2*r_e/(2*K.bohr_radius()))(K.coulomb_konstante()*K.e**2/(K.m_e*K.c**2)), (lambda r_e: 2/3*r_e**3*K.m_e*K.c**3/K.bohr_radius()**4)(K.coulomb_konstante()*K.e**2/(K.m_e*K.c**2)))/K.fall_offene_windung(), 1),
        ("05", "Bis zum Fall vergeht t = E/((l + 2)·P)", K.zeit_abgabe(2, K.G*K.M_sun*K.M_earth/(2*K.AU), K.leistung_quadrupol(K.M_earth, K.M_sun, K.AU))/(5/256*K.c**5*K.AU**4/(K.G**3*K.M_earth*K.M_sun*(K.M_earth + K.M_sun))), 1),
        ("05", "Die Erde fiele so in 1,1·10²³ Jahren auf die Sonne", K.zeit_abgabe(2, K.G*K.M_sun*K.M_earth/(2*K.AU), K.leistung_quadrupol(K.M_earth, K.M_sun, K.AU))/K.yr/1e23, 1.1),
        ("05", "Die Erde fiele so in 1,1·10²³ Jahren auf die Sonne", 2*math.pi*math.sqrt(K.AU**3/K.GM_sun)/(4*K.verlust_faltung(K.M_earth, K.M_sun, K.AU))/K.yr/1e23, 1.1),
        ("05", "t = (5/512)·(3Δs/c)·δ⁻⁴", K.verschmelzen_gleich(30*K.M_sun, 0.003)[0]/(lambda a: 3/8*2*math.pi*math.sqrt(a**3/(2*K.G*30*K.M_sun))/abs(K.doppelpulsar(30*K.M_sun, 30*K.M_sun, 2*math.pi*math.sqrt(a**3/(2*K.G*30*K.M_sun)), 0.0)))(3*K.faltung(30*K.M_sun)/0.003), 1),
        ("05", "N = δ^{−5/2}/(64√2·π)", K.verschmelzen_gleich(1.0, 0.01)[1]/(lambda d, n: (lambda a0: sum((lambda x: (K.G/(2*x**2))/K.leistung_quadrupol(1.0, 1.0, x)/(2*math.pi*math.sqrt(x**3/(2*K.G))))(a0*((i + 0.5)/n)**4)*4*((i + 0.5)/n)**3*a0/n for i in range(n)))(3*K.faltung(1.0)/d))(0.01, 200000), 1),
        ("05", "von δ = 0,001 an 1,1·10⁵ Umläufe", K.verschmelzen_gleich(K.M_sun, 0.001)[1]/1e5, 1.1),
        ("05", "von δ = 0,001 an 1,1·10⁵ Umläufe", (lambda d, n: (lambda a0: sum((lambda x: (K.GM_sun*K.M_sun/(2*x**2))/K.leistung_quadrupol(K.M_sun, K.M_sun, x)/(2*math.pi*math.sqrt(x**3/(2*K.GM_sun))))(a0*((i + 0.5)/n)**4)*4*((i + 0.5)/n)**3*a0/n for i in range(n)))(3*K.faltung(K.M_sun)/d))(0.001, 200000)/1e5, 1.1),
        ("05", "bei je einer Sonnenmasse in 13 Stunden", K.verschmelzen_gleich(K.M_sun, 0.001)[0]/3600, 13),
        ("05", "bei je einer Sonnenmasse in 13 Stunden", (lambda a: 3/8*2*math.pi*math.sqrt(a**3/(2*K.GM_sun))/abs(K.doppelpulsar(K.M_sun, K.M_sun, 2*math.pi*math.sqrt(a**3/(2*K.GM_sun)), 0.0)))(3*K.faltung(K.M_sun)/0.001)/3600, 13),
        ("05", "n = 3 über die Windung", (lambda P: math.log(P(2.0)/2.0/P(1.0))/math.log(2.0))(lambda W: K.mu0*1e30**2*W**4/(6*math.pi*K.c**3)), 3),
        ("05", "n = 5 über eine Verformung", (lambda P: math.log(P(2.0)/2.0/P(1.0))/math.log(2.0))(lambda W: 32/5*K.G/K.c**5*(1e38*1e-6)**2*W**6), 5),
        ("05", "seit höchstens t = Ω/((n − 1)·|Ω̇|) = E/(l·P)", (lambda l: sum(math.exp(-2*l*(i + 0.5)*math.log(1e4)/100000) for i in range(100000))*math.log(1e4)/100000/K.zeit_abgabe(l, 0.5, 1.0, drehung=True))(1), 1),
        ("05", "seit höchstens t = Ω/((n − 1)·|Ω̇|) = E/(l·P)", (lambda l: sum(math.exp(-2*l*(i + 0.5)*math.log(1e4)/100000) for i in range(100000))*math.log(1e4)/100000/K.zeit_abgabe(l, 0.5, 1.0, drehung=True))(2), 1),
        ("03", "r* = 0,36 kpc", K.schwellenradius(1e8 * K.M_sun) / K.kpc, 0.36),
        ("03", "bis 6,8 kpc", K.gezeitenrand(1e8 * K.M_sun, 6e10 * K.M_sun, 50 * K.kpc) / K.kpc, 6.8),
        ("03", "mit 36,0 km/s", K.bahntempo_regel(1e8 * K.M_sun, 2 * K.kpc) / 1e3, 36.0),
        ("03", "ν = 6,03", K.verstaerkung(K.G * 1e8 * K.M_sun / (2 * K.kpc)**2), 6.03),
        ("03", "Newton, 14,7 km/s", math.sqrt(K.G * 1e8 * K.M_sun / (2 * K.kpc)) / 1e3, 14.7),
        ("03", "= 120,5 km/s", K.v_flach(6e10 * K.M_sun) / math.sqrt(2) / 1e3, 120.5),
        ("03", "als 20,7 km/s", math.sqrt(2 * K.G * 1e8 * K.M_sun / (2 * K.kpc)) / 1e3, 20.7),
        ("03", "7,2·10⁶ M☉", (2 * K.kpc)**3 * K.schalter(K.G * 6e10 * K.M_sun / (50 * K.kpc)**2) * K.verstaerkung(K.G * 6e10 * K.M_sun / (50 * K.kpc)**2) * 6e10 * K.M_sun / (50 * K.kpc)**3 / 3 / K.M_sun / 1e6, 7.2),
        ("03", "5,0·10⁸ M☉", K.verstaerkung(K.G * 1e8 * K.M_sun / (2 * K.kpc)**2) - 1, 5.0),
        ("F:lokale_gruppe", "laufen die Schwerpunkte mit 42,1 km/s", D.m31_paar(True)[2], 42.1),
        ("F:lokale_gruppe", "laufen die Schwerpunkte mit 42,1 km/s", _bereiche_paar()[2], 42.1),
        ("F:lokale_gruppe", "je Richtung auf ±48 und ±44 km/s", D.m31_paar(True)[3], 48),
        ("F:lokale_gruppe", "je Richtung auf ±48 und ±44 km/s", D.m31_paar(True)[4], 44),
        ("F:lokale_gruppe", "vor 5,3 Gyr bei 54 kpc", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)[1] / K.Gyr, 5.3),
        ("F:lokale_gruppe", "vor 5,3 Gyr bei 54 kpc", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True, numerisch=True)[1] / K.Gyr, 5.3),
        ("F:lokale_gruppe", "vor 5,3 Gyr bei 54 kpc", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)[3] / K.kpc, 54),
        ("F:lokale_gruppe", "vor 5,3 Gyr bei 54 kpc", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True, numerisch=True)[3] / K.kpc, 54),
        ("F:lokale_gruppe", "nur 0,2 % der Baryonen der Milchstraße", (lambda x: 100 * (1 + x) * math.exp(-x))(K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)[3] * math.sqrt(2) / K.schwellenradius(K.M_B_MILCHSTRASSE * K.M_sun)), 0.2),
        ("F:lokale_gruppe", "2,5 % der von M31", (lambda x: 100 * (1 + x) * math.exp(-x))(K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)[3] * math.sqrt(2) / K.schwellenradius(K.M_B_M31 * K.M_sun)), 2.5),
        ("F:lokale_gruppe", "In 3,0 Gyr kommen sich beide wieder so nahe", (lambda z: (z[2] - 2 * z[1]) / K.Gyr)(K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)), 3.0),
        ("F:lokale_gruppe", "In 3,0 Gyr kommen sich beide wieder so nahe", (lambda z: (z[2] - 2 * z[1]) / K.Gyr)(K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True, numerisch=True)), 3.0),
        ("F:lokale_gruppe", "stehen ihre 12 hellsten Satelliten", float(_satelliten_phasen()["n"]), 12),
        ("F:lokale_gruppe", "bei 0,44 ± 0,08", _satelliten_phasen()["mittel"], 0.44),
        ("F:lokale_gruppe", "bei 0,44 ± 0,08", _satelliten_phasen()["mittel_bahn"], 0.44),
        ("F:lokale_gruppe", "bei 0,44 ± 0,08", _satelliten_phasen()["e"], 0.08),
        ("F:lokale_gruppe", "frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰)", _satelliten_phasen()["M"][0] / 1e10, 7.9),
        ("F:lokale_gruppe", "frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰)", _satelliten_phasen()["M"][1] / 1e10, 5.2),
        ("F:lokale_gruppe", "frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰)", _satelliten_phasen()["M"][2] / 1e10, 9.9),
        ("F:lokale_gruppe", "frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰)", min(_satelliten_phasen()["zweiter"]), 0.5),
        ("F:lokale_gruppe", "frei gelesen verlangen sie 7,9·10¹⁰ M☉ (5,2 bis 9,9·10¹⁰)", max(_satelliten_phasen()["zweiter"]), 0.5),
        ("F:lokale_gruppe", "gegen 6,6·10¹⁰ M☉ ihres Bereichs", _bereiche()[0] / K.M_sun / 1e10, 6.6),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("daten")["m_081"][2], 93),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("daten")["m_081"][4], 17),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("daten")["m_081"][3], 10),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", min(_wirte_satelliten("daten")[w][0] for w in ("m_031", "m_081", "ngc_3627", "m_083", "ngc_5128", "ngc_4594")), 0.39),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", max(_wirte_satelliten("daten")[w][0] for w in ("m_031", "m_081", "ngc_3627", "m_083", "ngc_5128", "ngc_4594")), 1.16),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", _wirte_satelliten("daten")["median"], 0.59),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("astropy")["m_081"][2], 93),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("astropy")["m_081"][4], 17),
        ("F:lokale_gruppe", "stehen die Satelliten im Mittel 93 ± 17 km/s über M81 und 10 km/s über dem Schwerpunkt mit M82", _wirte_satelliten("astropy")["m_081"][3], 10),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", min(_wirte_satelliten("astropy")[w][0] for w in ("m_031", "m_081", "ngc_3627", "m_083", "ngc_5128", "ngc_4594")), 0.39),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", max(_wirte_satelliten("astropy")[w][0] for w in ("m_031", "m_081", "ngc_3627", "m_083", "ngc_5128", "ngc_4594")), 1.16),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", _wirte_satelliten("astropy")["median"], 0.59),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", 1 / math.sqrt(3), 0.58),
        ("F:lokale_gruppe", "liegt sie bei σ/v_flach = 0,39 bis 1,16, im Median bei 0,59, 1/√3 = 0,58", math.sqrt((1 - 2 * 0.37 / 3) / (3 - 2 * 0.37)), 0.58),
        ("F:lokale_gruppe", "läuft die LMC zwischen 45 und 182 kpc", min(d for t, d in D.lmc_gruppe_rand()["wende"]), 45),
        ("F:lokale_gruppe", "läuft die LMC zwischen 45 und 182 kpc", _lmc_gruppe_bahn()["D_grenzen"][0], 45),
        ("F:lokale_gruppe", "läuft die LMC zwischen 45 und 182 kpc", max(d for t, d in D.lmc_gruppe_rand()["wende"]), 182),
        ("F:lokale_gruppe", "läuft die LMC zwischen 45 und 182 kpc", _lmc_gruppe_bahn()["D_grenzen"][1], 182),
        ("F:lokale_gruppe", "liegt er bei 19 kpc und nur Pictor II darin", float(D.lmc_gruppe_rand()["r_J"][0]), 19),
        ("F:lokale_gruppe", "liegt er bei 19 kpc und nur Pictor II darin", float(_lmc_gruppe_bahn()["r_J"][0]), 19),
        ("F:lokale_gruppe", "Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc", -D.lmc_gruppe_rand()["t_apo"], 1.47),
        ("F:lokale_gruppe", "Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc", -_lmc_gruppe_bahn()["t_apo"], 1.47),
        ("F:lokale_gruppe", "Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc", float(D.lmc_gruppe_rand()["r_J"][1]), 69),
        ("F:lokale_gruppe", "Im letzten Apozentrum, vor 1,47 Gyr, lag der Rand bei 69 kpc", float(_lmc_gruppe_bahn()["r_J"][1]), 69),
        ("F:lokale_gruppe", "und 7 der 8 Zwerge", float(D.lmc_gruppe_rand()["eigen"][0]), 7),
        ("F:lokale_gruppe", "und 7 der 8 Zwerge", float(_lmc_gruppe_bahn()["eigen"][0]), 7),
        ("F:lokale_gruppe", "und 7 der 8 Zwerge", float(D.lmc_gruppe_rand()["eigen"][1]), 8),
        ("F:lokale_gruppe", "liegt er bei 19 kpc und nur Pictor II darin", float(len(_lmc_gruppe_bahn()["heute"])), 1),
        ("F:lokale_gruppe", "liegt er bei 19 kpc und nur Pictor II darin", float(D.lmc_gruppe_rand()["heute"] == ["pictor_2"]), 1),
        ("F:lesungen", "bei 14 kpc das 2,2-fache von Newton", K.toomre_schwelle(14 * K.kpc / (K.schwellenradius(K.M_B_MILCHSTRASSE * K.M_sun) / math.sqrt(2)))["nu"], 2.2),
        ("F:lesungen", "bei 14 kpc das 2,2-fache von Newton", K.toomre_schwelle(14 * K.kpc / (K.schwellenradius(K.M_B_MILCHSTRASSE * K.M_sun) / math.sqrt(2)), numerisch=True)["nu"], 2.2),
        ("03", "voll ionisiert in der Urmischung 1,93", (2*(1-K.YP_GEMESSEN[0])+0.75*K.YP_GEMESSEN[0])*2/(2-K.YP_GEMESSEN[0]), 1.93),
        ("03", "kT_ab = 1,4 meV", K.kuehlung_licht()['kT_ab']*1e3, 1.4),
        ("03", "kT_ab = 1,4 meV", (135/(64*math.pi**3*K.alpha**2)*(2*(1-K.YP_GEMESSEN[0])+0.75*K.YP_GEMESSEN[0])*2/(2-K.YP_GEMESSEN[0])*math.sqrt(16*1.2020569031595942/(3*math.pi)*(1+K.Q_kosm)*(K.m_e/K.m_Pl)**2*K.faltung_je_photon()['eps_b']*K.e/(K.m_e*K.c**2)))**0.4*K.m_e*K.c**2/K.e*1e3, 1.4),
        ("03", "(16 K)", K.kuehlung_licht()['T_ab'], 16),
        ("03", "z_ab = 5,0", K.kuehlung_licht()['z_ab'], 5.0),
        ("03", "z_ab = 5,0", K.kuehlung_licht()['z_ab_form'], 5.0),
        ("03", "noch 0,68 der Temperatur des Lichts", math.exp(0.4)*0.4**0.4*(math.gamma(0.6)-0.4**0.6*math.exp(-0.4)*sum(0.4**n/math.prod(0.6+k for k in range(n+1)) for n in range(40))), 0.68),
        ("03", "2,8·10⁻¹⁰ der Energie des Lichts", 1.5*K.faltung_je_photon()['eps_b']*K.e/(K.faltung_je_photon()['mu']*K.m_H*K.c**2)/(K.strahlung_je_sinn(1,3)/K.strahlung_je_sinn(1,2))*1e10, 2.8),
        ("03", "0,025/(1 + z)", K.h*K.NU_21CM/(K.k_B*K.T_CMB), 0.025),
        ("08", "wo γ² = 2", 1 / (1 - 2 * 3 * K.faltung(K.M_sun) / (12 * K.faltung(K.M_sun))), 2),
        ("08", "𝒯 = 5GM/r³", K.gezeiten_falte(K.M_sun, K.R_SONNE)['faktor'], 5),
        ("08", "M_V = √5·(¼/δ_m)^{3/2}·m", float(abs(K.gezeiten_falte(K.M_sun, K.R_SONNE, K.gezeiten_falte(K.M_sun, K.R_SONNE)['M_V'])['r_T'] / (12 * K.faltung(K.gezeiten_falte(K.M_sun, K.R_SONNE)['M_V'])) - 1) < 1e-9), 1),
        ("08", "M_V = 9,0·10⁷ M☉", K.gezeiten_falte(K.M_sun, K.R_SONNE)['M_V'] / K.M_sun / 1e6, 90),
        ("08", "bei 0,6 AE", K.gezeiten_falte(K.M_sun, K.R_SONNE, 1e6 * K.M_sun)['r_T'] / K.AU, 0.6),
        ("08", "gleich 4GM/r³", K.gezeiten_falte(K.M_sun, K.R_SONNE, kreisbahn=True)['faktor'], 4),
        ("08", "M_V = 2·(⅙/δ_m)^{3/2}·m", float(abs(K.gezeiten_falte(K.M_sun, K.R_SONNE, K.gezeiten_falte(K.M_sun, K.R_SONNE, kreisbahn=True)['M_V'], kreisbahn=True)['r_T'] / (18 * K.faltung(K.gezeiten_falte(K.M_sun, K.R_SONNE, kreisbahn=True)['M_V'])) - 1) < 1e-9), 1),
        ("08", "für die Sonne 4,4·10⁷ M☉", K.gezeiten_falte(K.M_sun, K.R_SONNE, kreisbahn=True)['M_V'] / K.M_sun / 1e6, 44),
        ("03", "tief ¾", K.relaxation_regel(1e-12*K.a0)[0], 0.75),
        ("03", "zwei Drittel nah", math.log(K.glaettungslaenge(K.M_sun, K.gezeit_flach(1e4, 3e19))/(K.G*K.M_sun/1e4**2))/math.log(3e19/(K.G*K.M_sun/1e4**2)), 0.6667),
        ("03", "tief 2·ln(R/b_min)/N", 8*K.relaxation_regel(1e-12*K.a0)[1], 2),
        ("03", "ein Viertel des Werts nach Newton", K.relaxation_regel(1e-12*K.a0)[1], 0.25),
        ("03", "(ν = 3,3) 3,5-mal", K.verstaerkung(0.13*K.a0), 3.3),
        ("03", "(ν = 3,3) 3,5-mal", K.relaxation_regel(0.13*K.a0)[2], 3.5),
        ("03", "(ν = 10,5) 29-mal", K.verstaerkung(0.01*K.a0), 10.5),
        ("03", "(ν = 10,5) 29-mal", K.relaxation_regel(0.01*K.a0)[2], 29),
        ("03", "dort 1,1 und 2,9", K.relaxation_regel(0.13*K.a0)[3], 1.1),
        ("03", "dort 1,1 und 2,9", K.relaxation_regel(0.01*K.a0)[3], 2.9),
        ("03", "(X² = 6)", 4*3/2, 6),
        ("03", "sind das 0,0074", K.verdampfen_anteil(math.sqrt(6)), 0.0074),
        ("03", "(γ = 9/2)", (K.v_flach(K.M_sun)/K.streuung_kugel(K.M_sun))**2, 4.5),
        ("03", "je Relaxation 1,0·10⁻³", K.verdampfen_anteil(math.sqrt(4.5/2*(2*math.log(10)-1+0.01)))*1e3, 1.0),
        ("03", "von der Schwelle aus 1,2·10⁻¹⁵", K.verdampfen_anteil(K.flucht_regel(1e7*K.M_sun, K.schwellenradius(1e7*K.M_sun))/(math.sqrt(2)*K.streuung_kugel(1e7*K.M_sun)))*1e15, 1.2),
        ("03", "√5/(3π) = 0,080", 1000 * K.netz_achsen()["knoten"], 80),
        ("03", "√5/(3π) = 0,080", 1000 * K.netz_achsen(numerisch=True)["knoten"], 80),
        ("03", "zwei bei 0,420", 1000 * K.netz_achsen(numerisch=True)["faden"], 420),
        ("03", "einer bei 0,420", 1000 * K.netz_achsen(numerisch=True)["wand"], 420),
        ("03", "keiner bei 0,080", 1000 * K.netz_achsen(numerisch=True)["leer"], 80),
        ("03", "mindestens zwei positive Eigenwerte", 2 * (K.netz_achsen(numerisch=True)["knoten"] + K.netz_achsen(numerisch=True)["faden"]), 1),
        ("03", "3σ/√(10π) = 0,535·σ", 1000 * K.netz_achsen()["l1"], 535),
        ("03", "3σ/√(10π) = 0,535·σ", 1000 * K.netz_achsen(numerisch=True)["l1"], 535),
        ("03", "trägt 1/3 von Σλ²", 3 * K.netz_achsen(numerisch=True)["spur"], 1),
        ("03", "die Scherung 2/3", 1.5 * (1 - K.netz_achsen(numerisch=True)["spur"]), 1),
        ("03", "(σ²/15)", 3 / __import__('scipy.integrate', fromlist=['quad']).quad(lambda m: m**4 / 2, -1, 1)[0], 15),
        ("03", "je ein Achtel", 8 * 0.5**3, 1),
        ("03", "fällt sie nur mit λ > 0,71", K.wachstum_ende(1.020)[2], 0.71),
        ("03", "d_eff = 3 − χ, tief 2", 3 - K.schalter(1e-8 * K.a0), 2),
        ("03", "in der Tiefe 1, 2, 2", 3 - K.schalter(1e-8 * K.a0), 2),
        ("03", "v² = 2(1 + Q)·Gμ_b", __import__('scipy.integrate', fromlist=['quad']).quad(lambda z: (1 + z * z)**-1.5, -math.inf, math.inf)[0], 2),
        ("01", "(29 + 6√6)/(29 − 6√6) = 3,055", (29 + 6 * math.sqrt(6)) / (29 - 6 * math.sqrt(6)), 3.055),
        ("01", "(29 + 6√6)/(29 − 6√6) = 3,055", (lambda n: n['faden'] / n['knoten'])(K.netz_achsen('stellen', numerisch=True)), 3.055),
        ("01", "im Mittel 6,11 Fäden", (lambda n: 2 * n['faden'] / n['knoten'])(K.netz_achsen('stellen')), 6.11),
        ("01", "im Mittel 6,11 Fäden", (lambda n: 2 * n['faden'] / n['knoten'])(K.netz_achsen('stellen', numerisch=True)), 6.11),
        ("01", "2b₁ = 6", 2 * K.B1_ZELLE, 6),
        ("03", "nicht 1, 1,5, 2", 1 + (2 - 1) * (1 - K.schalter(1e-8 * K.a0) / 2), 1.5),
        ("10", "Je hundertfach kleineren Verlust wächst die Verweilzeit im Kanal zehnfach", (1 / (1 - K.rueckkehr_dimension(1, 2.0, 0.01 / math.sqrt(100)))) / (1 / (1 - K.rueckkehr_dimension(1, 2.0, 0.01))), 10),
        ("03", "im Feld der Sterne c/v_flach ≤ 0,58·c/σ_Gas", 1 / math.sqrt(3), 0.58),
        ("03", "im Feld der Sterne c/v_flach ≤ 0,58·c/σ_Gas", (lambda M: (K.G * M * K.a0 / 9)**0.25 / K.v_flach(M))(1e10 * K.M_sun), 0.58),
        ("02", "Carter, 5,9·10⁻³⁹", K.kopplung_grav(K.m_p) * 1e39, 5.9),
        ("02", "gegen 2,0·10⁻³⁹", K.alpha**12 * (K.m_e / K.m_p)**4 * 1e40, 20),
        ("01", "auf jede Stelle gleich", 1e9 * max(abs(K.exponent_inhalt(f, x)) for f in (lambda: K.verstaerkung(0.1 * K.a0), K.a0_aus_atmen, K.q_je_radiant, K.theta_, lambda: float(sum(D.sparc_reste(K.kpc, K.a0, D.YD)[1]**2)), lambda: D.kids_chi2()[0], lambda: D.hiflugcs_reste()[2], lambda: D.xcop_staub()[0], lambda: abs(D.udg_reste()[2]), lambda: D.abstaende_pruefung(1 - K.atmen_heute()[1] - K.dichte_heute()[0] / (K.rho_krit100 * (K.atmen_heute()[0] / 100)**2), K.dichte_heute()[0] / (K.rho_krit100 * (K.atmen_heute()[0] / 100)**2), frei=False)["modell"]) for x in ("alpha", "m_e", "m_p")), 0),
        ("01", "kT_ab ∝ α^{−4/5}·m_e^{6/5}", K.exponent_inhalt(lambda: K.kuehlung_licht()["kT_ab_form"], "alpha"), -0.8),
        ("01", "kT_ab ∝ α^{−4/5}·m_e^{6/5}", K.exponent_inhalt(lambda: K.kuehlung_licht()["kT_ab_form"], "m_e"), 1.2),
        ("01", "bis auf den Faktor 3", K.kopplung_grav(K.m_p) / (K.alpha**12 * (K.m_e / K.m_p)**4), 3),
        ("01", "3 % auseinander", 100 * (1 - K.atom_und_schliessung()[0] / K.atom_und_schliessung()[1]), 3),
        ("03", "(D/H und Y_p, ω_b = 0,0219 ± 0,0002)", D.kernbildung_wb()[0], 0.0219),
        ("03", "(D/H und Y_p, ω_b = 0,0219 ± 0,0002)", D.kernbildung_wb()[1], 0.0002),
        ("03", "ist (1 + Q)·ω_b = 0,1384 ± 0,0009", D.geometrie_q('kernbildung')['wm'], 0.1384),
        ("03", "ist (1 + Q)·ω_b = 0,1384 ± 0,0009", D.geometrie_q('kernbildung')['wm_fehler'], 0.0009),
        ("03", "und Q = 5,32 ± 0,06", D.geometrie_q('kernbildung')['Q'], 5.32),
        ("03", "und Q = 5,32 ± 0,06", D.geometrie_q('kernbildung')['Q_fehler'], 0.06),
        ("03", "2π − 1 liegt dann Δχ² = 0,3 daneben, das gemessene Verhältnis Δχ² = 0,6", D.geometrie_q('kernbildung')['d_radiant'], 0.3),
        ("03", "2π − 1 liegt dann Δχ² = 0,3 daneben, das gemessene Verhältnis Δχ² = 0,6", D.geometrie_q('kernbildung')['d_spitzen'], 0.6),
        ("03", "ist Q = 5,40 ± 0,08", 100 * D.geometrie_q('kernbildung', True)['Q'], 540),
        ("03", "ist Q = 5,40 ± 0,08", D.geometrie_q('kernbildung', True)['Q_fehler'], 0.08),
        ("03", "Der Zustandsteil hebt Q um 0,08", D.geometrie_q('kernbildung', True)['Q'] - D.geometrie_q('kernbildung')['Q'], 0.08),
        ("03", "Der Zustandsteil hebt Q um 0,08", K.q_messung()[0] - K.q_je_radiant(), 0.08),
        ("03", "Gegen das Atmen aller Präzisionsdaten (H₀ = 68,25)", D.atmen_praezision()[1], 68.25),
        ("03", "linear H_v = (1 − f·δ_L/3)·H mit f = 0,52 heute", K.leerraum_noetig(D.atmen_praezision()[0], D.sh0es_leiter()[0], 0.15)[0], 0.52),
        ("03", "linear mindestens δ_L = −0,41 tief", K.leerraum_noetig(D.atmen_praezision()[0], D.sh0es_leiter()[0], 0.15)[1], -0.41),
        ("03", "0,59 Gpc im Radius messen", K.leerraum_noetig(D.atmen_praezision()[0], D.sh0es_leiter()[0], 0.15)[2] / 1000, 0.59),
        ("03", "Diese Tiefe ist das 43-fache der Schwankung", K.leerraum_noetig(D.atmen_praezision()[0], D.sh0es_leiter()[0], 0.15)[5], 43),
        ("03", "für jeden Rand bis z = 0,3 mindestens Δχ² = 41 daneben", (lambda a: min(D.pantheon_leerraum(a[3], a[4], r / 100, D.sh0es_leiter()[0] / a[1] - 1, a[0][5], a[0][6])[2] for r in range(2, 31)))(D.atmen_praezision()), 41),
        ("03", "(r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2", K.schallhorizont(K.q_je_radiant(), D.atmen_praezision()[0][0]), 147.4),
        ("03", "(r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2", (lambda a: D.bao_absolut(None, a[3], a[4], K.schallhorizont(K.q_je_radiant(), a[0][0]), a[0][5], a[0][6], zmax=1.0)[0])(D.atmen_praezision()), 68.3),
        ("03", "(r_d = 147,4 Mpc) H₀ = 68,3 ± 0,2", (lambda a: D.bao_absolut(None, a[3], a[4], K.schallhorizont(K.q_je_radiant(), a[0][0]), a[0][5], a[0][6], zmax=1.0)[1])(D.atmen_praezision()), 0.2),
        ("03", "die Cepheiden 73,2 ± 1,1", (lambda a: D.pantheon_h0(a[3], a[4], 0.25, 1.0, a[0][5], a[0][6])[0])(D.atmen_praezision()), 73.2),
        ("03", "die Cepheiden 73,2 ± 1,1", (lambda a: D.pantheon_h0(a[3], a[4], 0.25, 1.0, a[0][5], a[0][6])[1])(D.atmen_praezision()), 1.1),
        ("03", "um 5·lg(73,2/68,3) = 0,15 mag dunkler", (lambda a: 5 * math.log10(D.pantheon_h0(a[3], a[4], 0.25, 1.0, a[0][5], a[0][6])[0] / D.bao_absolut(None, a[3], a[4], K.schallhorizont(K.q_je_radiant(), a[0][0]), a[0][5], a[0][6], zmax=1.0)[0]))(D.atmen_praezision()), 0.15),
        ("F:lesungen", "68,25 ± 0,10, sind es 4,7 σ", D.atmen_praezision()[2], 0.1),
        ("F:lesungen", "68,25 ± 0,10, sind es 4,7 σ", (D.sh0es_leiter()[0] - D.atmen_praezision()[1]) / math.hypot(D.sh0es_leiter()[1], D.atmen_praezision()[2]), 4.7),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["trgb"][0], 70.4),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["trgb"][1], 1.9),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", float(D.cchp_h0()["trgb"][2]), 24),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["trgb_jwst"][0], 68.8),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["trgb_jwst"][1], 2.2),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["jagb_jwst"][0], 67.8),
        ("F:lesungen", "liest H₀ = 70,4 ± 1,9 (TRGB mit 24 Supernovae), mit JWST allein 68,8 ± 2,2 (TRGB) und 67,8 ± 2,7 (JAGB)", D.cchp_h0()["jagb_jwst"][1], 2.7),
        ("F:lesungen", "1,1, 0,3 und 0,2 Fehler neben dem Atmen", abs(D.cchp_h0()["trgb"][0] - D.atmen_praezision()[1]) / math.hypot(D.cchp_h0()["trgb"][1], D.atmen_praezision()[2]), 1.1),
        ("F:lesungen", "1,1, 0,3 und 0,2 Fehler neben dem Atmen", abs(D.cchp_h0()["trgb_jwst"][0] - D.atmen_praezision()[1]) / math.hypot(D.cchp_h0()["trgb_jwst"][1], D.atmen_praezision()[2]), 0.3),
        ("F:lesungen", "1,1, 0,3 und 0,2 Fehler neben dem Atmen", abs(D.cchp_h0()["jagb_jwst"][0] - D.atmen_praezision()[1]) / math.hypot(D.cchp_h0()["jagb_jwst"][1], D.atmen_praezision()[2]), 0.2),
        ("F:lesungen", "das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4", D.erosita_cgm()["q"][1][1] / D.erosita_cgm()["sf"][1][1], 2.7),
        ("F:lesungen", "das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4", D.erosita_cgm()["q"][1][0], 11.1),
        ("F:lesungen", "das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4", D.erosita_cgm()["q"][2][1] / D.erosita_cgm()["sf"][2][1], 7.5),
        ("F:lesungen", "das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4", D.erosita_cgm()["sf"][2][0], 11.3),
        ("F:lesungen", "das 2,7-fache der Leuchtkraft bei lg M* ≈ 11,1 und das 7,5-fache bei 11,3 bis 11,4", D.erosita_cgm()["q"][2][0], 11.4),
        ("03", "ν_G(−ln q) − 1 = 16,8 bei Q = 2π − 1", K.zustandsteil_je_gas(K.q_je_radiant()), 16.8),
        ("03", "heute 1,21 % der Baryonen als neutrales Gas in Galaxien (ALFALFA α.100)", D.gasanteil_heute(), 1.21),
        ("03", "unabhängig davon 1,17 % (HIPASS)", D.gasanteil_hipass(), 1.17),
        ("03", "bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019)", D.gas_gezaehlt()[3], 0.34),
        ("03", "bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019)", D.gas_gezaehlt()[1], 1.12),
        ("03", "bei z = 0,34 das 1,12 ± 0,27-fache von heute (Bera u. a. 2019)", D.gas_gezaehlt()[2], 0.27),
        ("03", "verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 %", K.a0 * 1e10, 1.06),
        ("03", "verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 %", 100 * D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['punkt']['F0'], 0.73),
        ("03", "verlangt der Standardwert 1,06 als Formteil heute F₀ = 0,73 ± 0,11 %", 100 * D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['fehler']['F0'], 0.11),
        ("03", "sechs Zehntel des gezählten Gases", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['punkt']['F0'] / D.gas_gezaehlt()[0], 0.6),
        ("03", "Δχ² = 2,6 unter ΛCDM", D.planck_anpassung(None, ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'))['chi2']['gesamt'] - D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['chi2']['gesamt'], 2.6),
        ("03", "den Formteil a₀ = 1,043 ± 0,003", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse'), zustand='gezaehlt')['punkt']['form'], 1.043),
        ("03", "den Formteil a₀ = 1,043 ± 0,003", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse'), zustand='gezaehlt')['fehler']['a0'], 0.003),
        ("03", "(DESI DR2 mit r_d aus ω_b, Pantheon+, D/H) 1,05 ± 0,01", D.abstaende_gas()['form'], 1.05),
        ("03", "(DESI DR2 mit r_d aus ω_b, Pantheon+, D/H) 1,05 ± 0,01", D.abstaende_gas()['fehler_form'], 0.01),
        ("03", "alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gezaehlt')['punkt']['form'], 1.047),
        ("03", "alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gezaehlt')['fehler']['a0'], 0.003),
        ("03", "alle zusammen 1,047 ± 0,003, Δχ² = 4,1 unter ΛCDM", D.planck_anpassung(None, ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'))['chi2']['gesamt'] - D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gezaehlt')['chi2']['gesamt'], 4.1),
        ("03", "Der Standardwert mit allem gezählten Gas liegt Δχ² = 17 über ΛCDM", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), fest=('a0',), zustand='gezaehlt')['chi2']['gesamt'] - D.planck_anpassung(None, ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'))['chi2']['gesamt'], 17),
        ("03", "je Baryon für v_flach = 200 bis 20 km/s das 27- bis 266-fache", (lambda v: K.inaktives_mass(v**4 / (K.G * K.a0), K.atmen_heute()[0] * 1e3 / K.Mpc, K.a0, K.umkehr_heute()[1] * (1 - math.pi**2 / 15 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2) / (2 * (1 + K.q_je_radiant()))) / K.zustandsteil_je_gas(K.q_je_radiant()))(200e3), 27),
        ("03", "je Baryon für v_flach = 200 bis 20 km/s das 27- bis 266-fache", (lambda v: K.inaktives_mass(v**4 / (K.G * K.a0), K.atmen_heute()[0] * 1e3 / K.Mpc, K.a0, K.umkehr_heute()[1] * (1 - math.pi**2 / 15 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2) / (2 * (1 + K.q_je_radiant()))) / K.zustandsteil_je_gas(K.q_je_radiant()))(20e3), 266),
        ("03", "nur 0,003 bis 0,03 % der Baryonen als dünnes Gas", 100 * D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['punkt']['F0'] / (lambda v: K.inaktives_mass(v**4 / (K.G * K.a0), K.atmen_heute()[0] * 1e3 / K.Mpc, K.a0, K.umkehr_heute()[1] * (1 - math.pi**2 / 15 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2) / (2 * (1 + K.q_je_radiant()))) / K.zustandsteil_je_gas(K.q_je_radiant()))(20e3), 0.003),
        ("03", "nur 0,003 bis 0,03 % der Baryonen als dünnes Gas", 100 * D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['punkt']['F0'] / (lambda v: K.inaktives_mass(v**4 / (K.G * K.a0), K.atmen_heute()[0] * 1e3 / K.Mpc, K.a0, K.umkehr_heute()[1] * (1 - math.pi**2 / 15 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2) / (2 * (1 + K.q_je_radiant()))) / K.zustandsteil_je_gas(K.q_je_radiant()))(200e3), 0.03),
        ("F:lesungen", "bis z = 0,34 einen Anstieg auf das 2,2-fache", D.gas_gezaehlt()[3], 0.34),
        ("F:lesungen", "bis z = 0,34 einen Anstieg auf das 2,2-fache", (lambda p: K.zustand_gas(1 / (1 + D.gas_gezaehlt()[3]), p['Fs'], p['k'], p['a_s']) / K.zustand_gas(1.0, p['Fs'], p['k'], p['a_s']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas')['punkt']), 2.2),
        ("F:lesungen", "mit dem gezählten Gang liegen sie Δχ² = 1,6 höher", D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['chi2']['gesamt'] + D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas_bera')['chi2']['bera'] - D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='gas')['chi2']['gesamt'], 1.6),
        ("F:lesungen", "den linearen Teil bei 0,83 ± 0,07", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(None, None))[4], 0.83),
        ("F:lesungen", "den linearen Teil bei 0,83 ± 0,07", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(None, None))[6], 0.07),
        ("F:lesungen", "den Rest bei 0,68 ± 0,06 (χ² = 269)", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(None, None))[5], 0.68),
        ("F:lesungen", "den Rest bei 0,68 ± 0,06 (χ² = 269)", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(None, None))[7], 0.06),
        ("F:lesungen", "den Rest bei 0,68 ± 0,06 (χ² = 269)", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(None, None))[0], 269),
        ("F:lesungen", "bei k = 1/Mpc 0,99, bei 5/Mpc 0,87", float(D.xcop_gas_faltung(1.0)[0]), 0.99),
        ("F:lesungen", "bei k = 1/Mpc 0,99, bei 5/Mpc 0,87", float(D.xcop_gas_faltung(5.0)[0]), 0.87),
        ("F:lesungen", "10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio", D.warp_gaia()[0], 10.86),
        ("F:lokale_gruppe", "(Local Volume Database, Σm_i/M = 0,08)", D.reflex_scheibe()[3], 0.08),
        ("F:lokale_gruppe", "sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°)", D.reflex_scheibe()[0], 22.5),
        ("F:lokale_gruppe", "sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°)", D.reflex_scheibe(weg="galaktisch")[0], 22.5),
        ("F:lokale_gruppe", "sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°)", D.reflex_scheibe()[1], 82),
        ("F:lokale_gruppe", "sind das 22,5 km/s zum Ziel (l, b) = (82°, −43°)", D.reflex_scheibe(weg="galaktisch")[2], -43),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", D.reflex_gemessen()["v"], 32),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", D.reflex_gemessen()["v_unten"], 4),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", D.reflex_gemessen()["l"], 56),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", D.reflex_gemessen()["l_fehler"], 9),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", D.reflex_gemessen()["b"], -34),
        ("F:lokale_gruppe", "mit 32 ± 4 km/s zum Ziel (56 ± 9°, −34 ± 10°)", round(D.reflex_gemessen()["b_fehler"]), 10),
        ("F:lokale_gruppe", "trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben", 100 * D.reflex_scheibe()[0] / D.reflex_gemessen()["v"], 70),
        ("F:lokale_gruppe", "trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben", (D.reflex_gemessen()["v"] - D.reflex_scheibe()[0]) / D.reflex_gemessen()["v_unten"], 2.4),
        ("F:lokale_gruppe", "trägt der Impuls der Baryonen so 70 % davon, 2,4 Fehler darunter, und die Richtung liegt 22° daneben", (lambda a, m: math.degrees(math.acos(math.sin(math.radians(a[2])) * math.sin(math.radians(m["b"])) + math.cos(math.radians(a[2])) * math.cos(math.radians(m["b"])) * math.cos(math.radians(a[1] - m["l"])))))(D.reflex_scheibe(), D.reflex_gemessen()), 22),
        ("F:lesungen", "10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio", D.warp_gaia()[1], 0.03),
        ("F:lesungen", "10,86 ± 0,03 ± 3,20 km/s/kpc (Poggio", D.warp_gaia()[2], 3.2),
        ("F:lesungen", "an 12 Millionen Riesensternen", D.warp_gaia()[3] / 1e6, 12),
        ("F:lesungen", "6,5 bis 6,7 Fehler neben der Messung", _warp_drehung()["70"][2], 6.5),
        ("F:lesungen", "6,5 bis 6,7 Fehler neben der Messung", _warp_drehung()["100"][2], 6.7),
        ("F:lesungen", "rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung", _warp_drehung()["frei70"][1], 9.4),
        ("F:lesungen", "rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung", _warp_drehung()["frei100"][1], 9.4),
        ("F:lesungen", "rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung", _warp_drehung()["frei100"][2], 0.4),
        ("F:lesungen", "rechtläufig mit 9,4 km/s/kpc, 0,4 bis 0,5 Fehler neben der Messung", _warp_drehung()["frei70"][2], 0.5),
        ("F:lesungen", "kippen die Satelliten den Rand um 0,06°", _warp_drehung()["frei70"][0], 0.06),
        ("F:lesungen", "kippen die Satelliten den Rand um 0,06°", _warp_drehung()["frei100"][0], 0.06),
        ("F:lesungen", "dreht die Knotenlinie rückläufig mit 13,9 km/s/kpc", -_warp_drehung()["rund"][1], 13.9),
        ("F:lesungen", "mit dem 25-fachen Antrieb der LMC rückläufig mit 21,8 km/s/kpc, und der Rand kippt dann um 0,49°", -_warp_drehung()["lmc25"][1], 21.8),
        ("F:lesungen", "mit dem 25-fachen Antrieb der LMC rückläufig mit 21,8 km/s/kpc, und der Rand kippt dann um 0,49°", _warp_drehung()["lmc25"][0], 0.49),
        ("F:lesungen", "nur um 0,04° gegen die innere Scheibe", math.degrees(K.warp_ringe(((D.satellit('lmc')['M_b'] * K.M_sun, D.satellit('lmc')['x'], D.satellit('lmc')['v']), (D.satellit('smc')['M_b'] * K.M_sun, D.satellit('smc')['x'], D.satellit('smc')['v']), (D.satellit('sagittarius_1', True)['M_b'] * K.M_sun, D.satellit('sagittarius_1', True)['x'], D.satellit('sagittarius_1', True)['v'], D.satellit('sagittarius_1')['M_b'] * K.M_sun)), K.M_B_MILCHSTRASSE * K.M_sun, gemeinsam=True, eltern=[0, 1, 0])[0]), 0.04),
        ("F:lesungen", "nur um 0,04° gegen die innere Scheibe", math.degrees(K.warp_ringe(((D.satellit('lmc')['M_b'] * K.M_sun, D.satellit('lmc')['x'], D.satellit('lmc')['v']), (D.satellit('smc')['M_b'] * K.M_sun, D.satellit('smc')['x'], D.satellit('smc')['v']), (D.satellit('sagittarius_1', True)['M_b'] * K.M_sun, D.satellit('sagittarius_1', True)['x'], D.satellit('sagittarius_1', True)['v'], D.satellit('sagittarius_1')['M_b'] * K.M_sun)), K.M_B_MILCHSTRASSE * K.M_sun, gemeinsam=True, eltern=[0, 1, 0], N=100, n_phi=240)[0]), 0.04),
        ("F:lesungen", "rückläufig, mit 10 bis 11 km/s/kpc", -K.warp_ringe(((D.satellit('lmc')['M_b'] * K.M_sun, D.satellit('lmc')['x'], D.satellit('lmc')['v']), (D.satellit('smc')['M_b'] * K.M_sun, D.satellit('smc')['x'], D.satellit('smc')['v']), (D.satellit('sagittarius_1', True)['M_b'] * K.M_sun, D.satellit('sagittarius_1', True)['x'], D.satellit('sagittarius_1', True)['v'], D.satellit('sagittarius_1')['M_b'] * K.M_sun)), K.M_B_MILCHSTRASSE * K.M_sun, gemeinsam=True, eltern=[0, 1, 0])[1] * K.kpc / 1e3, 10),
        ("F:lesungen", "rückläufig, mit 10 bis 11 km/s/kpc", -K.warp_ringe(((D.satellit('lmc')['M_b'] * K.M_sun, D.satellit('lmc')['x'], D.satellit('lmc')['v']), (D.satellit('smc')['M_b'] * K.M_sun, D.satellit('smc')['x'], D.satellit('smc')['v']), (D.satellit('sagittarius_1', True)['M_b'] * K.M_sun, D.satellit('sagittarius_1', True)['x'], D.satellit('sagittarius_1', True)['v'], D.satellit('sagittarius_1')['M_b'] * K.M_sun)), K.M_B_MILCHSTRASSE * K.M_sun, gemeinsam=True, eltern=[0, 1, 0], N=100, n_phi=240)[1] * K.kpc / 1e3, 11),
        ("F:lesungen", "(v/c)² = 4·10⁻⁷ des Umlaufs", (K.bahntempo_regel(K.M_B_MILCHSTRASSE * K.M_sun, 14 * K.kpc) / K.c)**2 * 1e7, 4),
        ("F:lesungen", "27 liegen dort −0,03 ± 0,02 dex auf der Regel", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[0], 27),
        ("F:lesungen", "27 liegen dort −0,03 ± 0,02 dex auf der Regel", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[1], -0.03),
        ("F:lesungen", "27 liegen dort −0,03 ± 0,02 dex auf der Regel", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[2], 0.02),
        ("F:lesungen", "die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[4], 7),
        ("F:lesungen", "die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[5], -0.01),
        ("F:lesungen", "die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[6], 0.03),
        ("F:lesungen", "die 7 unter 45° Neigung −0,01 ± 0,03, keine tiefer als −0,27", D.sparc_gasreich_aussen(D.udg_reste()[3], D.udg_reste()[4])[3], -0.27),
        ("F:lesungen", "die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°)", (D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt'))[5] - D.udg_neigung()[3]) / D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt'))[6], 7),
        ("F:lesungen", "die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°)", D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt'))[5], 32),
        ("F:lesungen", "die Regel verlangt 10°, sieben Fehler unter der Neigung des Gases (32 ± 3°)", D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt'))[6], 3),
        ("F:lesungen", "über alle 15 Bins um einen festen Faktor höher", float(D.kids_typ_reste("rot", bezug="blau")[0]), 15),
        ("F:lesungen", "+0,22 ± 0,03 dex (χ² = 17)", D.kids_typ_reste("rot", bezug="blau")[1], 0.22),
        ("F:lesungen", "+0,22 ± 0,03 dex (χ² = 17)", D.kids_typ_reste("rot", bezug="blau")[2], 0.03),
        ("F:lesungen", "+0,22 ± 0,03 dex (χ² = 17)", (lambda v: v[3] - (v[1] / v[2])**2)(D.kids_typ_reste("rot", bezug="blau")), 17),
        ("F:lesungen", "frei angepasstem f = 3,1 ± 0,2", D.kids_typ_gedaechtnis("rot")[0], 3.1),
        ("F:lesungen", "frei angepasstem f = 3,1 ± 0,2", D.kids_typ_gedaechtnis("rot")[1], 0.2),
        ("F:lesungen", "(χ² = 27, die Regel allein 387)", D.kids_typ_gedaechtnis("rot")[2], 27),
        ("F:lesungen", "(χ² = 27, die Regel allein 387)", D.kids_typ_gedaechtnis("rot")[3], 387),
        ("F:lesungen", "bei den blauen f = 1,05 ± 0,12", D.kids_typ_gedaechtnis("blau")[0], 1.05),
        ("F:lesungen", "bei den blauen f = 1,05 ± 0,12", D.kids_typ_gedaechtnis("blau")[1], 0.12),
        ("03", "beide aus dem Auslassen je eines Haufens auf ± 0,8 und ± 2,5 genau", D.haufen_c_fehler("xcop"), 0.8),
        ("03", "beide aus dem Auslassen je eines Haufens auf ± 0,8 und ± 2,5 genau", D.haufen_c_fehler("hiflugcs"), 2.5),
        ("03", "Um die Spitze ist die mittlere lineare Faltung innen flach", float(max(abs(K.faltung_spitze((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 0.02 * R, R) / K.faltung_spitze((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 0.2 * R, R) - 1) for R in (25.0, 18.0)) < 0.02), 1),
        ("03", "und die inneren Schalen kehren zugleich um", float((lambda e: (lambda t: max(t) / min(t) < 1.05)([x for x, m in zip(e["t_um"], e["M_schale"]) if m < 5e13 * K.M_sun]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=300))), 1),
        ("03", "Fallen sie radial, wird das Profil innen steiler als NFW", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=300)), 68.24915)[0] > 20), 1),
        ("03", "wird das Profil innen flacher als NFW", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=300, mischen=8)), 68.24915)[0] < 1), 1),
        ("07", "An SPARC (2688 Punkte)", sum(len(x[2]['R']) for x in D.sparc_auswahl()), 2688),
        ("07", "sind das höchstens 2·10⁻⁵ dex", 1e5 * max(v[1] for v in D.nome_wahl('baryonen').values()), 2),
        ("07", "Δχ² = 0,01 an SPARC", D.nome_wahl('baryonen')['sparc'][2], 0.01),
        ("07", "0,06 an HIFLUGCS", D.nome_wahl('baryonen')['hiflugcs'][2], 0.06),
        ("07", "als Faltung 1·10⁻⁴ dex", 1e4 * max(v[1] for v in D.nome_wahl('feld').values()), 1),
        ("07", "unter 1/6", K.faltung(K.M_sun) / K.grenzradius(K.M_sun), 0.166667),
        ("07", "Sonne 4,5", math.log(3 * K.c**2 / K.v_flach(K.M_sun)**2) / (2 * math.pi), 4.5),
        ("07", "Milchstraße 2,6", math.log(3 * K.c**2 / K.v_flach(K.M_B_MILCHSTRASSE * K.M_sun)**2) / (2 * math.pi), 2.6),
        ("08", "Steigung 1 hält die Spreizung GM/(Rc²)", (lambda f: f(K.M_sun) / f(1e12 * K.M_sun))(lambda M: K.G * M / (18 * K.faltung(M) * K.c**2)), 1),
        ("08", "½ die Neigung GM/R²", (lambda f: f(K.M_sun) / f(1e12 * K.M_sun))(lambda M: K.G * M / K.schwellenradius(M)**2), 1),
        ("08", "⅓ die Gezeit GM/R³ und mit ihr die mittlere Dichte", (lambda f: f(1e20) / f(1e26))(lambda M: K.G * M / (K.bohr_radius() * (M / K.m_p)**(1 / 3))**3), 1),
        ("08", "im Tiefbereich ¼ den Umlauf, Ω⁴ = GMa₀/R⁴", (lambda f: f(K.M_sun) / f(1e12 * K.M_sun))(lambda M: K.G * M * K.a0 / K.groesste_umkehr(M)[0]**4), 1),
        ("08", "im Tiefbereich ¼ den Umlauf, Ω⁴ = GMa₀/R⁴", (lambda M: (K.v_flach(M) / K.groesste_umkehr(M)[0])**4 / (K.G * M * K.a0 / K.groesste_umkehr(M)[0]**4))(6e10 * K.M_sun), 1),
        ("08", "Steigung 0 hält eine Zeit, ℓ₀ = c·τ_D", K.c * K.rueckkehrzeit() / (K.c**2 / K.a0), 1),
        ("08", "Auf −1 hält λ̄ die Wirkung M·λ̄·c = ħ", (lambda f: f(K.m_e) / f(K.M_sun))(lambda M: M * K.compton(M) * K.c / K.hbar), 1),
        ("08", "auf −½ hält r̃ die Zeit der eigenen Quantenbewegung, M·r̃²/ħ = τ_D", (lambda f: f(K.m_p) / f(K.M_sun))(lambda M: M * K.schwelle_spiegel(M)**2 / K.hbar / K.rueckkehrzeit()), 1),
        ("08", "liegt bei der Masse ihres Treffpunkts mal Z^{1/(k₁ − k₂)}", (lambda M: K.leiter_treffen((1, 0), (0.5, 0.5), Z=3 * K.faltung(M) / K.schwellenradius(M))[0] / M)(6e10 * K.M_sun), 1),
        ("08", "mit dem Exponenten 2 zwischen Spreizung und Neigung", (lambda M: (K.schwellenradius(M) / (3 * K.faltung(M)))**2 / (K.schnittmassen()[3] / M))(6e10 * K.M_sun), 1),
        ("08", "3/2 zwischen Spreizung und Gezeit", (lambda L, lh: K.leiter_treffen((1, math.log(lh / K.compton(K.m_p)) / L), (1 / 3, math.log(lh / K.l_P * (K.m_Pl / K.m_p)**(1 / 3)) / L))[0] / K.sternmasse())(math.log(1 / K.leiterzahl()), K.bohr_radius()), 1),
        ("08", "6 zwischen Neigung und Gezeit", (lambda L, lh: K.leiter_treffen((0.5, 0.5), (1 / 3, math.log(lh / K.l_P * (K.m_Pl / K.m_p)**(1 / 3)) / L))[0] / K.kreuzung_gehalten(lh)[1])(math.log(1 / K.leiterzahl()), K.bohr_radius()), 1),
        ("08", "4 zwischen Neigung und Umlauf im Tiefbereich", (lambda M, H: K.groesste_umkehr(M)[0] / K.schwellenradius(M) / (K.a0**3 / (K.G * H**4) / M)**0.25)(6e10 * K.M_sun, K.groesste_umkehr(K.M_sun)[1]), 1),
        ("08", "in m_p·(α/α_G)^{3/2} bei δ = α²·m_e/m_p", (lambda L, lh, d: K.leiter_treffen((1, math.log(1 / d) / L), (1 / 3, math.log(lh / K.l_P * (K.m_Pl / K.m_p)**(1 / 3)) / L))[0] / (K.m_p * (K.alpha / K.kopplung_grav(K.m_p))**1.5))(math.log(1 / K.leiterzahl()), K.bohr_radius(), K.alpha**2 * K.m_e / K.m_p), 1),
        ("08", "darum liegt ihr Treffpunkt für jeden Halter bei M_*", (lambda L: (lambda t: t(K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3)) / K.sternmasse() * t(1e-6) / t(K.bohr_radius()))(lambda lh: K.leiter_treffen((1, math.log(lh / K.compton(K.m_p)) / L), (1 / 3, math.log(lh / K.l_P * (K.m_Pl / K.m_p)**(1 / 3)) / L))[0]))(math.log(1 / K.leiterzahl())), 1),
        ("08", "M_V = √5·m_p·(λ_H/(12Δs_p))^{3/2}", (lambda lh: K.gezeiten_falte(10 * K.M_sun, lh * (10 * K.M_sun / K.m_p)**(1 / 3))['M_V'] / (math.sqrt(5) * K.m_p * (lh / (12 * K.faltung(K.m_p)))**1.5))(K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3)), 1),
        ("08", "M_V = √5·m_p·(λ_H/(12Δs_p))^{3/2}", (lambda lh, MV: K.gezeiten_falte(10 * K.M_sun, lh * (10 * K.M_sun / K.m_p)**(1 / 3), MV)['r_T'] / (12 * K.faltung(MV)))(K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3), math.sqrt(5) * K.m_p * (K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3) / (12 * K.faltung(K.m_p)))**1.5), 1),
        ("08", "fällt dagegen die Masse des Sterns heraus", (lambda lh: K.gezeiten_falte(0.3 * K.M_sun, lh * (0.3 * K.M_sun / K.m_p)**(1 / 3))['M_V'] / K.gezeiten_falte(10 * K.M_sun, lh * (10 * K.M_sun / K.m_p)**(1 / 3))['M_V'])(K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3)), 1),
        ("08", "auf der Kreisbahn 2·m_p·(λ_H/(18Δs_p))^{3/2}", (lambda lh, MV: K.gezeiten_falte(3 * K.M_sun, lh * (3 * K.M_sun / K.m_p)**(1 / 3), MV, kreisbahn=True)['r_T'] / (18 * K.faltung(MV)))(K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3), 2 * K.m_p * (K.R_SONNE / (K.M_sun / K.m_p)**(1 / 3) / (18 * K.faltung(K.m_p)))**1.5), 1),
        ("08", "parallel zu den Kanten der Reihe nach M_U", K.leiter_treffen((1 / 3, 2 / 3), (0.5, 0.5))[0] / K.schnittmassen()[3], 1),
        ("08", "der Reihe nach M_U, M_min", K.leiter_treffen((-0.5, 0.5), (0, 1))[0] / K.schnittmassen()[0], 1),
        ("08", "M_min und m_Pl", K.leiter_treffen((-3, 0), (3, 0))[0] / K.m_Pl, 1),
        ("08", "parallel zur Seitenhalbierenden von λ̄ nach r* M_c", K.leiter_treffen((-0.5, 0.5), (0.25, 0.75))[0] / K.schnittmassen()[1], 1),
        ("08", "zu der von 3Δs nach r̃ m_Pl²/M_c", K.leiter_treffen((3, 0), (0, 1))[0] / (K.m_Pl**2 / K.schnittmassen()[1]), 1),
        ("08", "M_U = c⁴/(G·a₀) ohne ħ", K.c**4 / (K.G * K.a0) / K.schnittmassen()[3], 1),
        ("08", "M_min = ħa₀/c³ ohne G", K.hbar * K.a0 / K.c**3 / K.schnittmassen()[0], 1),
        ("08", "m_Pl = √(ħc/G) ohne a₀", math.sqrt(K.hbar * K.c / K.G) / K.m_Pl, 1),
        ("08", "GMa₀/c⁴ = M/M_U", K.G * K.M_sun * K.a0 / K.c**4 / (K.M_sun / K.schnittmassen()[3]), 1),
        ("08", "nur um einen Faktor der Form", 6 * K.faltung(K.schnittmassen()[3] / 4) / K.schwellenradius(K.schnittmassen()[3] / 4), 1),
        ("08", "eine Dekade der Masse ist in x 0,016", math.log(10) / math.log(1 / K.leiterzahl()), 0.016),
        ("04", "g ist darin die halbe Steigung von c²·(dτ/dt)² über r", (lambda M, r, h: 0.5 * K.c**2 * (K.takt(M, r + h)**2 - K.takt(M, r - h)**2) / (2 * h) / K.neigung(M, r))(K.M_sun, 3 * K.grenzradius(K.M_sun), 1.0), 1),
        ("04", "Sie spürt die Neigung als g/(dτ/dt)", (lambda M, r, h: K.c**2 * (K.takt(M, r + h) - K.takt(M, r - h)) / (2 * h) / K.waerme_im_takt(M, r, 1.0)[1])(K.M_sun, 3 * K.grenzradius(K.M_sun), 1.0), 1),
        ("04", "Im gleichmäßig gehaltenen Körper ist T/T_R = 1 + P/(ρc²)", (lambda m: (1 + m[1]) * m[0] / K.takt(1.4 * K.M_sun, 12e3))(K.mitte_gehalten(1.4 * K.M_sun, 12e3, numerisch=True)), 1),
        ("04", "mindestens 13 % wärmer als an der Hülle", 100 * (K.takt(1.4 * K.M_sun, 12e3) / K.mitte_gehalten(1.4 * K.M_sun, 12e3)[0] - 1), 13),
        ("04", "mindestens 13 % wärmer als an der Hülle", 100 * K.mitte_gehalten(1.4 * K.M_sun, 12e3, numerisch=True)[1], 13),
        ("01", "Laufen beide parallel, d = n + 2", (lambda n, h: __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: (lambda E: E(1 + h) - 2 * E(1.0) + E(1 - h))(lambda r: r**-n / n - r**(2 - d) / (d - 2)), 2.5, 6.0) - n)(1.5, 1e-4), 2),
        ("01", "Laufen beide parallel, d = n + 2", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.atmen_gehalten([(d, 1.0)], n=1.5), 0.0, 10.0) - 1.5, 2),
        ("01", "ihr Halter ist der Umlauf, L²/r³ mit n = 2", -math.log(2.0 * 2.0**-3) / math.log(2.0), 2),
        ("01", "Die 4 in κ² = (4 − d_eff)·Ω² ist n + 2", (lambda d, h: (lambda F: -(F(1 + h) - F(1 - h)) / (2 * h) + d)(lambda r: r**-3 - r**(1 - d)))(2.5, 1e-6), 4),
        ("01", "Die 4 in κ² = (4 − d_eff)·Ω² ist n + 2", K.atmen_gehalten([(2.5, 1.0)], n=2.0) + 2.5, 4),
        ("01", "darum zählt die Rückstellung, d = 0, vierfach", (lambda H2, h: (1 + (lambda F: (F(1 + h) - F(1 - h)) / (2 * h))(lambda r: (1 - H2) / r**3 - 1 / r**2 + H2 * r)) / H2)(0.09, 1e-6), 4),
        ("01", "darum zählt die Rückstellung, d = 0, vierfach", -K.atmen_gehalten([(0, -1.0)], n=2.0), 4),
        ("01", "Nach Newton, d = 3, ist K = (n − 1)·⟨r·g⟩/⟨r²⟩", (lambda n, h: (lambda E: (E(1 + h) - 2 * E(1.0) + E(1 - h)) / h**2 / 0.6)(lambda r: 0.6 / n * r**-n - 0.6 / r))(1.5, 1e-4), 0.5),
        ("01", "Nach Newton, d = 3, ist K = (n − 1)·⟨r·g⟩/⟨r²⟩", K.atmen_gehalten([(3, 1.0)], n=1.5), 0.5),
        ("01", "die Bindung je Masse ist −((n − 1)/n)·⟨r·g⟩", (lambda n: (0.6 / n - 0.6) / 0.6)(1.5), -0.333),
        ("01", "Gas, n = 2, dreht", K.atmen_gehalten([(3, 1.0)], beta=1.0) + 1, 2),
        ("01", "bindet mit der Hälfte", (lambda n: -(0.6 / n - 0.6) / 0.6)(2.0), 0.5),
        ("01", "wie die Epizyklen an seiner Hülle, K = GM/R³", (lambda M, R, h: (lambda E: (E(R * (1 + h)) - 2 * E(R) + E(R * (1 - h))) / (R * h)**2)(lambda r: 0.3 * K.G * M**2 / R * (R / r)**2 - 0.6 * K.G * M**2 / r) / (0.6 * M) / (K.G * M / R**3))(K.M_sun, K.R_SONNE, 1e-3), 1),
        ("01", "Der Halter bei c, n = 1", K.atmen_gehalten([(3, 1.0)], beta=0.0) + 1, 1),
        ("01", "läuft parallel zu Newton", 1e6 * abs((lambda n, h: (lambda E: (E(1 + h) - 2 * E(1.0) + E(1 - h)) / h**2)(lambda r: 0.6 / n * r**-n - 0.6 / r))(1.0, 1e-4)), 0),
        ("01", "seine Bindung ist null", (lambda n: 0.6 / n - 0.6)(1.0), 0),
        ("01", "5/18 bei β = ½", K.atmen_gehalten([(3, 1.0)], beta=0.5), 0.278),
        ("01", "5/18 bei β = ½", (lambda b, h, br: (lambda r0: (lambda s0: (lambda P: 3 * (math.log(P(r0 * (1 + h))) - math.log(P(r0 * (1 - h)))) / (math.log(1 + h) - math.log(1 - h)) - 4)(lambda r: (lambda T: r * T + T**4)(br(lambda T: 1.5 * math.log(T) - math.log(r) + 4 * T**3 / r - s0, 0.5, 2.0, xtol=1e-15))))(-math.log(r0) + 4 / r0))(b / (1 - b)))(0.5, 1e-5, __import__('scipy.optimize', fromlist=['brentq']).brentq), 0.278),
        ("01", "β/2 bei kleinem β", K.atmen_gehalten([(3, 1.0)], beta=1e-4) / 1e-4, 0.5),
        ("01", "β/2 bei kleinem β", (lambda b, h, br: (lambda r0: (lambda s0: (lambda P: (3 * (math.log(P(r0 * (1 + h))) - math.log(P(r0 * (1 - h)))) / (math.log(1 + h) - math.log(1 - h)) - 4) / b)(lambda r: (lambda T: r * T + T**4)(br(lambda T: 1.5 * math.log(T) - math.log(r) + 4 * T**3 / r - s0, 0.5, 2.0, xtol=1e-15))))(-math.log(r0) + 4 / r0))(b / (1 - b)))(1e-3, 1e-4, __import__('scipy.optimize', fromlist=['brentq']).brentq), 0.5),
        ("01", "n − 1 = β·(4 − 3β)/(8 − 7β)", (lambda b, h, br: (lambda r0: (lambda s0: (lambda P: (3 * (math.log(P(r0 * (1 + h))) - math.log(P(r0 * (1 - h)))) / (math.log(1 + h) - math.log(1 - h)) - 4) / (b * (4 - 3 * b) / (8 - 7 * b)))(lambda r: (lambda T: r * T + T**4)(br(lambda T: 1.5 * math.log(T) - math.log(r) + 4 * T**3 / r - s0, 0.5, 2.0, xtol=1e-15))))(-math.log(r0) + 4 / r0))(b / (1 - b)))(0.3, 1e-5, __import__('scipy.optimize', fromlist=['brentq']).brentq), 1),
        ("01", "d über 3", float((lambda x, h: 1 - (math.log((1 + h)**-2 / (1 - x / (1 + h))) - math.log((1 - h)**-2 / (1 - x / (1 - h)))) / (math.log(1 + h) - math.log(1 - h)))(0.1, 1e-6) > 3), 1),
        ("01", "und dort ist K = 0", (lambda b, c, h, mn: (lambda Mx, E: (lambda x: (lambda M, R: 1e3 * (E(R * (1 + h), M) - 2 * E(R, M) + E(R * (1 - h), M)) / h**2 / (M * M / R))(Mx(x), Mx(x)**(1 / 3) / x))(math.exp(mn(lambda lx: -Mx(math.exp(lx)), bounds=(-5, 8), method='bounded', options={'xatol': 1e-12}).x)))(lambda x: ((-b / x + math.sqrt(b * b / x / x + 8 * c / math.sqrt(1 + x * x))) / (4 * c))**1.5, lambda R, M: M * (math.sqrt(1 + M**(2 / 3) / R**2) - 1) - b * M * M / R - c * M**3 / R**2))(1.0, 0.01, 1e-4, __import__('scipy.optimize', fromlist=['minimize_scalar']).minimize_scalar), 0),
        ("01", "Mit der Regel ist d = 3 − χ̄", (lambda M, R, h, quad: (lambda V, X: 2 - (math.log(V(R * (1 + h))) - math.log(V(R * (1 - h)))) / (math.log(1 + h) - math.log(1 - h)) + X(R) / V(R))(lambda R: quad(lambda x: 3 * x**2 * (x * R) * K.verstaerkung(K.G * M * x / R**2) * K.G * M * x / R**2, 0, 1, limit=200)[0], lambda R: quad(lambda x: 3 * x**2 * (x * R) * K.verstaerkung(K.G * M * x / R**2) * K.G * M * x / R**2 * K.schalter(K.G * M * x / R**2), 0, 1, limit=200)[0]))(1e8 * K.M_sun, K.schwellenradius(1e8 * K.M_sun), 1e-4, __import__('scipy.integrate', fromlist=['quad']).quad), 3),
        ("01", "Im Tiefbereich, d = 2", (lambda M, R, h, quad: (lambda V: 2 - (math.log(V(R * (1 + h))) - math.log(V(R * (1 - h)))) / (math.log(1 + h) - math.log(1 - h)))(lambda R: quad(lambda x: 3 * x**2 * (x * R) * K.verstaerkung(K.G * M * x / R**2) * K.G * M * x / R**2, 0, 1, limit=200)[0]))(1e8 * K.M_sun, 1e4 * K.schwellenradius(1e8 * K.M_sun), 1e-4, __import__('scipy.integrate', fromlist=['quad']).quad), 2),
        ("01", "hängt r·g nicht am Radius", (lambda M, R, quad: (lambda V: V(2 * R) / V(R))(lambda R: quad(lambda x: 3 * x**2 * (x * R) * K.verstaerkung(K.G * M * x / R**2) * K.G * M * x / R**2, 0, 1, limit=200)[0]))(1e8 * K.M_sun, 1e4 * K.schwellenradius(1e8 * K.M_sun), __import__('scipy.integrate', fromlist=['quad']).quad), 1),
        ("01", "bei R_max/∛4", (lambda M, H, h: __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R0: (lambda A: (lambda E: (E(R0 * (1 + h)) - 2 * E(R0) + E(R0 * (1 - h))) / (R0 * h)**2)(lambda r: A / r**2 - 0.6 * K.G * M**2 / r - 0.3 * H**2 * M * r**2))((0.6 * K.G * M**2 / R0**2 - 0.6 * H**2 * M * R0) * R0**3 / 2), 0.3 * (K.G * M / H**2)**(1 / 3), 0.9 * (K.G * M / H**2)**(1 / 3)) / (K.G * M / H**2)**(1 / 3))(1e12 * K.M_sun, K.groesste_umkehr(K.M_sun)[1], 1e-4), 0.63),
        ("01", "wie die letzte stabile Kreisbahn um einen Haufen", (lambda M, H: __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: K.atmen_gehalten([(3, K.G * M / R**3), (0, -H**2)], n=2.0), 1e19, 1e25) / K.letzte_stabile_bahn(M, regel=False)[0])(1e12 * K.M_sun, K.groesste_umkehr(K.M_sun)[1]), 1),
        ("01", "das Vierfache der Gezeitendichte des Ganzen", (lambda M, H: (lambda R: 3 * M / (4 * math.pi * R**3))(__import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: K.atmen_gehalten([(3, K.G * M / R**3), (0, -H**2)], n=2.0), 1e19, 1e25)) / K.gezeitendichte(H**2))(1e12 * K.M_sun, K.groesste_umkehr(K.M_sun)[1]), 4),
        ("08", "Jedes Glied F_i ∝ u^{d_i − 3}", (lambda aN, h: 3 + (math.log(K.verstaerkung(aN * (1 + h)**2) * aN) - math.log(K.verstaerkung(aN * (1 - h)**2) * aN)) / (math.log(1 + h) - math.log(1 - h)) + K.schalter(aN))(0.5 * K.a0, 1e-5), 3),
        ("08", "Σ w_i = 1", (0.1 - 3 * 0.1**2) / 0.1 + 3 * 0.1, 1),
        ("08", "eine Abweichung schwingt mit (Re τ_M)² = 4 − d_eff", (lambda u0, h: 1 - ((u0 - 3 * u0**2 + 3 * (u0 + h)**2) - (u0 - 3 * u0**2 + 3 * (u0 - h)**2)) / (2 * h))(0.1, 1e-6) + K.flussdimension_bahn([(3, 1 - 0.3), (5, 0.3)])[0], 4),
        ("08", "Newton d = 3", 1 - (math.log(K.G * K.M_sun / 2e11**2) - math.log(K.G * K.M_sun / 1e11**2)) / math.log(2.0), 3),
        ("08", "die Regel 3 − χ", (lambda aN, h: 1 - (math.log(K.verstaerkung(aN * (1 + h)**-2) * aN * (1 + h)**-2) - math.log(K.verstaerkung(aN * (1 - h)**-2) * aN * (1 - h)**-2)) / (math.log(1 + h) - math.log(1 - h)) + K.schalter(aN))(0.5 * K.a0, 1e-5), 3),
        ("08", "die Rückstellung d = 0 mit negativem Anteil", 1 - (math.log(K.groesste_umkehr(K.M_sun)[1]**2 * 2.0) - math.log(K.groesste_umkehr(K.M_sun)[1]**2 * 1.0)) / math.log(2.0), 0),
        ("08", "die Rückstellung d = 0 mit negativem Anteil", (lambda M: (lambda R, H: -H**2 / (K.v_flach(M)**2 / R**2 - H**2))(K.letzte_stabile_bahn(M)[0], K.groesste_umkehr(M)[1]))(6e10 * K.M_sun), -1),
        ("08", "Die Faltung fügt das Glied 3GMu²/c² hinzu, mit d = 5", 3 + (math.log(3 * 0.2**2) - math.log(3 * 0.1**2)) / math.log(2.0), 5),
        ("08", "darum trägt Newton nur 1 − 3δ, das Quadrat der umlaufenden Uhr", (0.1 - 3 * 0.1**2) / 0.1 / K.takte_kreisbahn(0.1)['umlaufend']**2, 1),
        ("08", "und die Faltung den Rest 3δ", (1 - K.takte_kreisbahn(0.1)['umlaufend']**2) / 0.1, 3),
        ("08", "So ist d_eff = 3 + 6δ", (K.flussdimension_bahn([(3, K.takte_kreisbahn(0.1)['umlaufend']**2), (5, 1 - K.takte_kreisbahn(0.1)['umlaufend']**2)])[0] - 3) / 0.1, 6),
        ("08", "der Punkt beim Bruchteil 3δ auf der Strecke von 3 nach 5", (K.flussdimension_bahn([(3, 1 - 0.3), (5, 0.3)])[0] - 3) / (5 - 3) / 0.1, 3),
        ("08", "und (Re τ_M)² = 1 − 6δ", K.flussdimension_bahn([(3, 1 - 0.3), (5, 0.3)])[1] / K.takte_kreisbahn(0.1)['radial'], 1),
        ("08", "und (Re τ_M)² = 1 − 6δ", (1 - K.flussdimension_bahn([(3, 1 - 0.21), (5, 0.21)])[1]) / 0.07, 6),
        ("08", "unter 4 die Drehung, bei 4 die Scherung, darüber die Dehnung", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.umlauf_abbildung(K.flussdimension_bahn([(d, 1.0)])[1])['spur'] - 2, 3.9, 4.1), 4),
        ("08", "Die Identität liegt bei 3 und 0", round(sum(abs(K.umlauf_abbildung(K.flussdimension_bahn([(d, 1.0)])[1])['matrix'][i][j] - (i == j)) for d in (0, 3) for i in (0, 1) for j in (0, 1)), 9), 0),
        ("08", "hält Re τ_M·Re τ_M′ = 2 fest", math.sqrt(K.flussdimension_bahn([(3, 1.0)])[1] * K.flussdimension_bahn([(0, 1.0)])[1]), 2),
        ("08", "sein Spiegelpunkt ist der tiefe Bereich bei 2", K.flussdimension_bahn([(3 - K.schalter(1e-12 * K.a0), 1.0)])[0], 2),
        ("08", "sein Spiegelpunkt ist der tiefe Bereich bei 2", K.flussdimension_bahn([(3 - K.schalter(1e-12 * K.a0), 1.0)])[1], 2),
        ("08", "Die halbe Drehung liegt bei 15/4", 4 - K.takte_kreisbahn(1 / 8)['radial'], 3.75),
        ("08", "Die halbe Drehung liegt bei 15/4", K.umlauf_abbildung(K.flussdimension_bahn([(3, 1 - 3 / 8), (5, 3 / 8)])[1])['spur'], -2),
        ("08", "beide über 15/4 hinweg auf 4", (lambda M: (lambda R, H: K.flussdimension_bahn([(2, K.v_flach(M)**2 / R**2), (0, -H**2)])[0])(K.letzte_stabile_bahn(M)[0], K.groesste_umkehr(M)[1]))(6e10 * K.M_sun), 4),
        ("08", "beide über 15/4 hinweg auf 4", (lambda M: (lambda R, H: K.flussdimension_bahn([(3, K.G * M / R**3), (0, -H**2)])[0])(K.letzte_stabile_bahn(M, False)[0], K.groesste_umkehr(M)[1]))(1e15 * K.M_sun), 4),
        ("08", "beide über 15/4 hinweg auf 4", (lambda M: (lambda Rm, H: K.flussdimension_bahn([(2, K.v_flach(M)**2 / (math.sqrt(7 / 15) * Rm)**2), (0, -H**2)])[0])(*K.groesste_umkehr(M)))(6e10 * K.M_sun), 3.75),
        ("08", "Innen liegt die Scherung in der Mitte der Strecke", 3 / __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda d: K.flussdimension_bahn([(3, 1 - 3 * d), (5, 3 * d)])[0] - (3 + 5) / 2, 0.01, 0.3), 18),
        ("08", "bei 18Δs tragen Newton und Faltung je die Hälfte", K.takte_kreisbahn(3 / 18)['umlaufend']**2, 0.5),
        ("08", "bei 18Δs tragen Newton und Faltung je die Hälfte", K.takte_kreisbahn(3 / 18)['radial'], 0),
        ("08", "der Photonensphäre 9Δs", 3 * K.G * K.M_sun / K.c**2 / K.faltung(K.M_sun), 9),
        ("08", "ist d_eff = 5 und (Re τ_M)² = −1", K.flussdimension_bahn([(3, K.takte_kreisbahn(3 / 9)['umlaufend']**2), (5, 1 - K.takte_kreisbahn(3 / 9)['umlaufend']**2)])[0], 5),
        ("08", "ist d_eff = 5 und (Re τ_M)² = −1", (lambda u0, h: 1 - (3 * (u0 + h)**2 - 3 * (u0 - h)**2) / (2 * h))(1 / 3, 1e-6), -1),
        ("08", "Newton trägt nichts", round(K.takte_kreisbahn(3 / 9)['umlaufend']**2, 12), 0),
        ("05", "t = E/(|k − 1|·P) auch für k unter eins", (lambda P, M, n: (1 - math.log(P(2*M)/P(M))/math.log(2))*P(M)*sum(K.c**2/P(M*(i + 0.5)/n)*M/n for i in range(n))/(M*K.c**2))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.M_sun, 20000), 1),
        ("05", "Steigung 1 hält die Zeit E/P fest", (lambda L: (2*K.M_sun/L(2*K.M_sun))/(K.M_sun/L(K.M_sun)))(lambda M: 4*math.pi*K.G*M*K.m_p*K.c/(8*math.pi/3*(K.alpha*K.compton(K.m_e))**2)), 1),
        ("05", "P ∝ M⁻², k = −2", (lambda P: math.log(P(2*K.M_sun)/P(K.M_sun))/math.log(2))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2), -2),
        ("05", "und lebt E/(3P)", (lambda P, M, n: M*K.c**2/(P(M)*sum(K.c**2/P(M*(i + 0.5)/n)*M/n for i in range(n))))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.M_sun, 20000), 3),
        ("05", "An der Grenze wächst sie mit k = 1", (lambda L: math.log(L(2*K.M_sun)/L(K.M_sun))/math.log(2))(lambda M: 4*math.pi*K.G*M*K.m_p*K.c/(8*math.pi/3*(K.alpha*K.compton(K.m_e))**2)), 1),
        ("05", "um den Faktor e je Zeit E/P", (lambda eps, M: K.wachsen_grenze(eps)['t_e']*(1 - eps)/eps*4*math.pi*K.G*M*K.m_p*K.c/(8*math.pi/3*(K.alpha*K.compton(K.m_e))**2)/(M*K.c**2))(1 - 2*math.sqrt(2)/3, 1e8*K.M_sun), 1),
        ("05", "Ẽ^{2 − k}", (lambda P, E: (lambda Et: math.log(Et(1.0)/E(1.0)*P(1.0)/(Et(2.0)/E(2.0)*P(2.0)))/math.log(Et(1.0)/Et(2.0)) + math.log(P(1.0)/P(2.0))/math.log(E(1.0)/E(2.0)))(lambda a: (K.m_Pl*K.c**2)**2/(2*E(a))))(lambda a: K.leistung_quadrupol(K.M_sun, K.M_sun, a*K.AU), lambda a: K.G*K.M_sun**2/(2*a*K.AU)), 2),
        ("05", "kehrt k − 1 um", (lambda P, T, M: (math.log(T(M)/M*P(M)/(T(2*M)/(2*M)*P(2*M)))/math.log(T(M)/T(2*M)) - 1)/(math.log(P(2*M)/P(M))/math.log(2) - 1))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.falte_temperatur, K.M_sun), -1),
        ("05", "und hält die Zeit fest", (lambda P, M: K.spiegelmasse(M)*K.c**2/(3*(K.spiegelmasse(M)/M)*P(M))/K.lebensdauer_falte(M))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.M_sun), 1),
        ("05", "Fest bleibt nur das Wachsen an der Grenze", (lambda P: math.log(K.spiegelmasse(K.M_sun)/K.M_sun*P(K.M_sun)/(K.spiegelmasse(2*K.M_sun)/(2*K.M_sun)*P(2*K.M_sun)))/math.log(K.spiegelmasse(K.M_sun)/K.spiegelmasse(2*K.M_sun)))(lambda M: 4*math.pi*K.G*M*K.m_p*K.c/(8*math.pi/3*(K.alpha*K.compton(K.m_e))**2)), 1),
        ("05", "nimmt ihre Spiegelmasse in derselben Zeit um den Faktor e ab", (lambda eps, M: (lambda P: K.spiegelmasse(M)*K.c**2/(K.spiegelmasse(M)/M*P)/K.wachsen_grenze(eps)['t_e'])((1 - eps)/eps*4*math.pi*K.G*M*K.m_p*K.c/(8*math.pi/3*(K.alpha*K.compton(K.m_e))**2)))(1 - 2*math.sqrt(2)/3, 1e8*K.M_sun), 1),
        ("05", "wird aus k = −2 gespiegelt 2 − k = 4", 2 - (lambda P: math.log(P(2*K.M_sun)/P(K.M_sun))/math.log(2))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2), 4),
        ("05", "steigt mit Ṫ ∝ T⁴", (lambda P, T, M: math.log(T(M)/M*P(M)/(T(2*M)/(2*M)*P(2*M)))/math.log(T(M)/T(2*M)))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.falte_temperatur, K.M_sun), 4),
        ("05", "nach derselben Zeit E/(3P)", (lambda P, M, n: M*K.c**2/(P(M)*sum(K.falte_temperatur(M)/((i + 0.5)/n)**2/(K.falte_temperatur(M*(i + 0.5)/n)/(M*(i + 0.5)/n)*P(M*(i + 0.5)/n)/K.c**2)/n for i in range(n))))(lambda M: math.pi**2/60*(K.k_B*K.falte_temperatur(M))**4/(K.hbar**3*K.c**2)*4*math.pi*(6*K.faltung(M))**2, K.M_sun, 20000), 3),
        ("08", "nur bei der Spur 0 oder ±1", float(sorted({q for q in range(1, 25) for p in range(1, q) if math.gcd(p, q) == 1 and abs(K.umlauf_abbildung((p/q)**2)['spur'] - round(K.umlauf_abbildung((p/q)**2)['spur'])) < 1e-9}) == [2, 3, 4, 6]), 1),
        ("08", "bei Re τ_M = ¼ und ¾, die Spur null", sum(abs(K.umlauf_abbildung(f*f)['spur']) for f in (0.25, 0.75)), 0),
        ("08", "jede Richtung kehrt nach zwei Umläufen wieder", sum((lambda m: (lambda m2: abs(m2[0][1]) + abs(m2[1][0]) + abs(m2[0][0] - m2[1][1]))([[sum(m[i][k]*m[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]))(K.umlauf_abbildung(f*f)['matrix']) for f in (0.25, 0.75)), 0),
        ("08", "die Bahn nahe am Kreis schließt nach vier", (lambda v: float(v[0]) if v[0] == v[1] else 0.0)([min(n for n in range(1, 13) if (lambda P: max(abs(P[i][j] - (i == j)) for i in (0, 1) for j in (0, 1)))(__import__('functools').reduce(lambda A, B: [[sum(A[i][k]*B[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)], [K.umlauf_abbildung(f*f)['matrix']]*n)) < 1e-9) for f in (0.25, 0.75)]), 4),
        ("08", "bei ⅙ und ⅓, ⅔ und ⅚, ist die Spur ±1", sum(abs(abs(K.umlauf_abbildung(f*f)['spur']) - 1) for f in (1/6, 1/3, 2/3, 5/6)), 0),
        ("08", "jede Richtung kehrt nach drei Umläufen wieder", sum((lambda m: (lambda m2: (lambda m3: abs(m3[0][1]) + abs(m3[1][0]) + abs(m3[0][0] - m3[1][1]))([[sum(m2[i][k]*m[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]))([[sum(m[i][k]*m[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]))(K.umlauf_abbildung(f*f)['matrix']) for f in (1/6, 1/3, 2/3, 5/6)), 0),
        ("08", "die Bahn schließt nach drei oder sechs", float([min(n for n in range(1, 13) if (lambda P: max(abs(P[i][j] - (i == j)) for i in (0, 1) for j in (0, 1)))(__import__('functools').reduce(lambda A, B: [[sum(A[i][k]*B[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)], [K.umlauf_abbildung(f*f)['matrix']]*n)) < 1e-9) for f in (1/6, 1/3, 2/3, 5/6)] == [6, 3, 3, 6]), 1),
        ("08", "zwei solche Folgen", (lambda sp: sum(1 for a, b in zip(sp, sp[1:]) if a*b < 0))([K.umlauf_abbildung(K.takte_kreisbahn(1e-5 + i*(1/6 - 2e-5)/20000)['radial'])['spur'] for i in range(20001)]), 2),
        ("08", "zwei solche Folgen", (lambda sp: sum(1 for a, b in zip(sp, sp[1:]) if (abs(a) - 1)*(abs(b) - 1) < 0)/2)([K.umlauf_abbildung(K.takte_kreisbahn(1e-5 + i*(1/6 - 2e-5)/20000)['radial'])['spur'] for i in range(20001)]), 2),
        ("08", "bis zur halben Drehung bei 24Δs", 3/__import__('scipy.optimize', fromlist=['brentq']).brentq(lambda x: K.takte_kreisbahn(x)['radial'] - 0.25, 1e-6, 1/6), 24),
        ("08", "Re τ_M von 1 bis √2", K.umlaufphase(1e3*K.a0), 1),
        ("08", "Re τ_M von 1 bis √2", K.umlaufphase(1e-12*K.a0)**2, 2),
        ("08", "durchläuft eine Folge", (lambda v: sum(1 for a, b in zip(v, v[1:]) if a*b < 0))([K.umlauf_abbildung(K.umlaufphase(10**(-3 + i*6/20000)*K.a0)**2)['spur'] for i in range(20001)]), 1),
        ("08", "erreicht die Identität bei 3/2 nicht", 1.5*float(round(K.umlauf_abbildung(1.5**2)['spur']) == -2 and K.umlaufphase(1e-12*K.a0) < 1.5), 1.5),
        ("08", "die Dreierdrehung bei a_N = 3,17·a₀", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda a: K.umlauf_abbildung(K.umlaufphase(a*K.a0)**2)['spur'] - 1, 1.5, 10.0), 3.17),
        ("08", "den Spiegel bei 1,12·a₀", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda a: K.umlauf_abbildung(K.umlaufphase(a*K.a0)**2)['spur'], 0.5, 3.0), 1.12),
        ("08", "und die Dreierdrehung bei 0,23·a₀", __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda a: K.umlauf_abbildung(K.umlaufphase(a*K.a0)**2)['spur'] + 1, 0.05, 0.5), 0.23),
        ("08", "d²u/dφ² + u = K_el·W/(p_φ·c)² + α_φ²·u", (lambda a: (lambda r: a * (math.sqrt(1 / r**2 + 1) - a / r) * r)(math.sqrt(1 / a**2 - 1)) / (1 - a**2))(0.6), 1),
        ("08", "Auf dem Kreis ist α_φ genau sein Tempo gegen c, bei jedem Tempo", (lambda a: (lambda r: (1 / r) / math.sqrt(1 / r**2 + 1))(math.sqrt(1 / a**2 - 1)) / a)(0.9), 1),
        ("08", "das Quadrat der umlaufenden Uhr, 1 − α_φ²", (lambda a: (lambda r: a * (math.sqrt(1 / r**2 + 1) - a / r) * r / (1 - ((1 / r) / math.sqrt(1 / r**2 + 1))**2))(math.sqrt(1 / a**2 - 1)))(0.7), 1),
        ("08", "das Quadrat der umlaufenden Uhr, 1 − α_φ²", K.bahn_windung(0.6, 0.05, numerisch=True)[1] / (1 - 0.6**2), 1),
        ("08", "den Rest α_φ² mit d = 4", 3 + (K.bahn_windung(0.6, 0.05, numerisch=True)[0] - 3) / 0.36, 4),
        ("08", "(Re τ_M)² = 1 − α_φ²", K.bahn_windung(0.9, 0.5, numerisch=True)[1] + 0.81, 1),
        ("08", "dreht als Rosette um denselben Winkel", K.bahn_windung(0.6, 0.8, numerisch=True)[1] / K.bahn_windung(0.6, 0.05, numerisch=True)[1], 1),
        ("08", "ein Drittel der 3δ der Faltung", (K.bahn_windung(math.sqrt(0.01))[0] - 3) / (1 - K.takte_kreisbahn(0.01)["umlaufend"]**2), 0.333),
        ("08", "die den Takt 2δ mitträgt", (1 - K.takte_kreisbahn(0.01)["ruhend"]**2) / 0.01, 2),
        ("08", "von 3 nach 4", K.bahn_windung(0.0)[0], 3),
        ("08", "von 3 nach 4", K.bahn_windung(1.0)[0], 4),
        ("08", "ist halb so lang wie die von 3 nach 5", (K.bahn_windung(1.0)[0] - 3) / (K.flussdimension_bahn([(5, 1.0)])[0] - 3), 0.5),
        ("08", "schwach um πα_φ² vor", K.bahn_windung(0.01, 0.3, numerisch=True)[2] / (math.pi * 0.01**2), 1),
        ("08", "ein Sechstel der Faltung", K.bahn_windung(1e-3)[2] / (2 * math.pi * (1 / math.sqrt(K.takte_kreisbahn(1e-6)["radial"]) - 1)), 0.167),
        ("08", "Die Scherung liegt bei K_el = p_φ·c", float(K.umlauf_abbildung(K.bahn_windung(1.0)[1])["art"] == "Scherung"), 1),
        ("08", "Die Scherung liegt bei K_el = p_φ·c", K.bahn_windung(0.99, 0.05, numerisch=True)[1] / (1 - 0.99**2), 1),
        ("08", "mit dem Drehimpuls ħ bei |Q₁Q₂| = q_P²", (K.e / math.sqrt(K.alpha))**2 / (4 * math.pi * K.eps0 * K.hbar * K.c), 1),
        ("08", "im Treffpunkt der Windung, Zα = 1", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda Z: K.bahn_windung(Z * K.e**2 / (4 * math.pi * K.eps0 * K.hbar * K.c))[1], 50, 200) * K.alpha, 1),
        ("08", "Darüber dehnt die Abbildung", float(K.umlauf_abbildung(K.bahn_windung(1.05)[1])["art"] == "Dehnung"), 1),
        ("08", "So erreicht die Bahn um eine Windung c an ihrer Scherung", (lambda a: (lambda r: (1 / r) / math.sqrt(1 / r**2 + 1))(math.sqrt(1 / a**2 - 1)))(1 - 1e-9), 1),
        ("08", "die um eine Faltung erst am Ende ihrer Strecke", K.takte_kreisbahn(1 / 3)["tempo_ort"], 1),
        ("08", "die um eine Faltung erst am Ende ihrer Strecke", K.flussdimension_bahn([(3, K.takte_kreisbahn(1 / 3)["umlaufend"]**2), (5, 1 - K.takte_kreisbahn(1 / 3)["umlaufend"]**2)])[0], 5),
        ("08", "(dτ/dt)² = 1 − 6Δs/r + r_Q·3Δs/r²", round(1e9 * abs((lambda M, Q, r: K.takt(M, r, Q)**2 / (1 - 2 * K.G * (M - __import__('scipy.integrate', fromlist=['quad']).quad(lambda x: K.eps0 / 2 * (Q / (4 * math.pi * K.eps0 * x * x))**2 * 4 * math.pi * x * x, r, math.inf)[0] / K.c**2) / (K.c**2 * r)))(10 * K.M_sun, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun, 5 * K.G * 10 * K.M_sun / K.c**2) - 1)), 0),
        ("08", "r_Q·3Δs = (ℓ_P·Q_el/q_P)²", round(1e9 * max(abs(K.e**2 / (4 * math.pi * K.eps0 * M * K.c**2) * 3 * K.faltung(M) / (K.l_P**2 * K.alpha) - 1) for M in (K.m_e, K.m_p, 10 * K.M_sun))), 0),
        ("08", "r² − 6Δs·r + a_J² + r_Q·3Δs = 0", round(1e9 * abs((lambda M, J, Q: K.horizonte(J, M, Q)[0].real / (K.G * M / K.c**2 + math.sqrt((K.G * M / K.c**2)**2 - (J / (M * K.c))**2 - K.G * Q**2 / (4 * math.pi * K.eps0 * K.c**4))))(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "Die Höhe im Thaleskreis wird die Hypotenuse aus a_J und √(r_Q·3Δs)", round(1e9 * abs((lambda M, J, Q: (lambda h: (h[0] * h[1]).real / (K.verwindung(J, M)**2 + Q**2 / (4 * math.pi * K.eps0 * M * K.c**2) * 3 * K.faltung(M)))(K.horizonte(J, M, Q)))(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "beide Horizonte treffen sich bei j² + r_Q/(3Δs) = 1", round(1e9 * abs((lambda M: (lambda d: (lambda J, Q2: ((K.G * M / K.c**2)**2 - (J / (M * K.c))**2 - K.G * Q2 / (4 * math.pi * K.eps0 * K.c**4)) / (K.G * M / K.c**2)**2)(0.6 * d * M * K.c, 0.64 * d * 4 * math.pi * K.eps0 * M * K.c**2))(3 * K.faltung(M)))(10 * K.M_sun))), 0),
        ("08", "beide Horizonte treffen sich bei j² + r_Q/(3Δs) = 1", round(1e6 * abs((lambda M: (lambda h: h[0] - h[1])(K.horizonte(0.6 * K.G * M**2 / K.c, M, 0.8 * math.sqrt(4 * math.pi * K.eps0 * K.G) * M)))(10 * K.M_sun)) / (3 * K.faltung(10 * K.M_sun))), 0),
        ("08", "r_Q/(3Δs) = (Q_el/q_P)²/α_G", round(1e9 * abs((lambda M, Q: Q**2 / (4 * math.pi * K.eps0 * M * K.c**2) / (3 * K.faltung(M)) / ((Q * math.sqrt(K.alpha) / K.e)**2 / K.kopplung_grav(M)))(K.m_p, 3 * K.e) - 1)), 0),
        ("08", "Ohne Drehung liegt der Treffpunkt bei Q_el/q_P = M/m_Pl", round(1e9 * abs((lambda M: (lambda Q: Q**2 / (4 * math.pi * K.eps0 * M * K.c**2) / (3 * K.faltung(M)))(K.e / math.sqrt(K.alpha) * M / K.m_Pl))(10 * K.M_sun) - 1)), 0),
        ("08", "r_e/(3Δs_e) = α/α_G(m_e) = 4,2·10⁴²", K.e**2 / (4 * math.pi * K.eps0 * K.m_e * K.c**2) / (3 * K.faltung(K.m_e)) / 1e42, 4.2),
        ("08", "r_e/(3Δs_e) = α/α_G(m_e) = 4,2·10⁴²", K.alpha / (K.m_e / K.m_Pl)**2 / 1e42, 4.2),
        ("08", "Ohne Drehung ist M_irr die Faltung innerhalb des Horizonts", round(1e9 * abs((lambda M, Q: (lambda rp: (M - Q**2 / (8 * math.pi * K.eps0 * K.c**2 * rp)) / K.falte_gedreht(M, 0.0, Q)['M_irr'])(K.horizonte(0.0, M, Q)[0].real))(10 * K.M_sun, 0.7 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "M = M_irr + Q_el²/(8πε₀c²r₊)", round(1e9 * abs((lambda M, Q: (K.falte_gedreht(M, 0.0, Q)['M_irr'] + Q**2 / (8 * math.pi * K.eps0 * K.c**2 * (K.G * M / K.c**2 + math.sqrt((K.G * M / K.c**2)**2 - K.G * Q**2 / (4 * math.pi * K.eps0 * K.c**4))))) / M)(10 * K.M_sun, 0.7 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "M² = (M_irr + M_el)² + M_rot²", round(1e9 * abs((lambda M, d: ((d['M_irr'] + d['M_el'])**2 + d['M_rot']**2) / M**2)(10 * K.M_sun, K.falte_gedreht(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun)) - 1)), 0),
        ("08", "M_irr·M_el = (m_Pl·Q_el/(2q_P))²", round(1e9 * abs((lambda Q, d: d['M_irr'] * d['M_el'] / (K.m_Pl * Q * math.sqrt(K.alpha) / (2 * K.e))**2)(0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun, K.falte_gedreht(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun)) - 1)), 0),
        ("08", "r₀(M_irr)·r₀(M_el) = r_Q·3Δs", round(1e9 * abs((lambda M, Q, d: K.grenzradius(d['M_irr']) * K.grenzradius(d['M_el']) / (Q**2 / (4 * math.pi * K.eps0 * M * K.c**2) * 3 * K.faltung(M)))(10 * K.M_sun, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun, K.falte_gedreht(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun)) - 1)), 0),
        ("08", "M_irr = M_el = M/2", round(1e6 * (lambda M, d: abs(d['M_irr'] / M - 0.5) + abs(d['M_el'] / M - 0.5))(10 * K.M_sun, K.falte_gedreht(10 * K.M_sun, 0.0, math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun))), 0),
        ("08", "Dort lässt sich die Hälfte der Ruheenergie herausholen", K.falte_gedreht(K.M_sun, 0.0, math.sqrt(4 * math.pi * K.eps0 * K.G) * K.M_sun)['herausholbar'], 0.5),
        ("08", "die Entropie ist ein Viertel der ungeladenen Falte gleicher Masse", K.falte_gedreht(K.M_sun, 0.0, math.sqrt(4 * math.pi * K.eps0 * K.G) * K.M_sun)['entropie'] / K.entropie_falte(K.M_sun), 0.25),
        ("08", "ħc·(r₊ − r₋)/(4πk·r₀(M_irr)²)", round(1e6 * abs((lambda M, Q: (lambda d, rQ, rp: K.hbar * K.c**2 * ((1 - 2 * d / (rp * (1 + 1e-7)) + rQ * d / (rp * (1 + 1e-7))**2) - (1 - 2 * d / (rp * (1 - 1e-7)) + rQ * d / (rp * (1 - 1e-7))**2)) / (2e-7 * rp) / 2 / (2 * math.pi * K.c * K.k_B) / K.falte_gedreht(M, 0.0, Q)['temperatur'])(3 * K.faltung(M), Q**2 / (4 * math.pi * K.eps0 * M * K.c**2), K.horizonte(0.0, M, Q)[0].real))(10 * K.M_sun, 0.7 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "ħc·(r₊ − r₋)/(4πk·r₀(M_irr)²)", round(1e6 * abs((lambda M, J, Q: (lambda m: K.c**2 * (m(1 + 1e-6) - m(1 - 1e-6)) / (2e-6 * 4 * math.pi * (K.falte_gedreht(M, J, Q)['M_irr'] / K.m_Pl)**2) / K.k_B / K.falte_gedreht(M, J, Q)['temperatur'])(lambda f: (lambda Mi: math.sqrt((Mi + Q**2 / (16 * math.pi * K.eps0 * K.G * Mi))**2 + (J * K.c / (2 * K.G * Mi))**2))(K.falte_gedreht(M, J, Q)['M_irr'] * math.sqrt(f))))(10 * K.M_sun, 0.5 * K.G * (10 * K.M_sun)**2 / K.c, 0.6 * math.sqrt(4 * math.pi * K.eps0 * K.G) * 10 * K.M_sun) - 1)), 0),
        ("08", "ohne Windung mit r₀(M_irr)² = r₀·r₊", round(1e9 * abs((lambda M, J: K.grenzradius(K.falte_gedreht(M, J)['M_irr'])**2 / (K.grenzradius(M) * K.horizonte(J, M)[0].real))(10 * K.M_sun, 0.6 * K.G * (10 * K.M_sun)**2 / K.c) - 1)), 0),
        ("09", "das Mittel aus Schwelle und Tiefe, a_J² = r*·3Δs", (K.v_flach(6e10 * K.M_sun) * K.schwellenradius(6e10 * K.M_sun) / K.c)**2 / (K.schwellenradius(6e10 * K.M_sun) * 3 * K.faltung(6e10 * K.M_sun)), 1),
        ("09", "für die Milchstraße (6·10¹⁰ M☉) 1759", K.c / K.v_flach(6e10 * K.M_sun), 1759),
        ("09", "ihre Horizonte stehen 89,97° zur Tiefe", math.degrees(math.acos(K.v_flach(6e10 * K.M_sun) / K.c)), 89.97),
        ("09", "für die Milchstraße (6·10¹⁰ M☉) 1759", K.v_flach(6e10 * K.M_sun) * K.schwellenradius(6e10 * K.M_sun) / K.c / (3 * K.faltung(6e10 * K.M_sun)), 1759),
        ("09", "eine Scheibe an der Schwelle mit J/M = 2R_d·v_flach hat √2·c/v_flach", 2 * K.schwellenradius(6e10 * K.M_sun) / math.sqrt(2) * K.v_flach(6e10 * K.M_sun) / K.c / (3 * K.faltung(6e10 * K.M_sun)) / (math.sqrt(2) * K.c / K.v_flach(6e10 * K.M_sun)), 1),
        ("09", "Ihre Spinzahl fällt mit der Masse wie M^{−1/4}", math.log((K.c / K.v_flach(1e12 * K.M_sun)) / (K.c / K.v_flach(1e10 * K.M_sun))) / math.log(100), -0.25),
        ("09", "erreicht eins erst am Treffpunkt M_U", K.c / K.v_flach(K.c**4 / (K.G * K.a0)), 1),
        ("09", "ist sie 2π/Θ_∞ = 5,1", K.c / K.v_flach(K.schwelle_umkehr()[0]), 5.1),
        ("09", "ist sie 2π/Θ_∞ = 5,1", 2 * math.pi / (math.sqrt(15) / math.pi), 5.1),
        ("03", "im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.verstaerkung(K.G * 6e10 * K.M_sun / r**2) * K.G * 6e10 * K.M_sun / r**2 - 3 * K.G * 1e15 * K.M_sun / K.Mpc**3 * r, K.kpc, 1e4 * K.kpc) / K.kpc, 50),
        ("03", "im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc", (lambda Om: (lambda f: f(f(f(f(f(f(f(f(50 * K.kpc)))))))))(lambda r: K.bahntempo_regel(6e10 * K.M_sun, r) / (math.sqrt(3) * Om)))(math.sqrt(K.G * 1e15 * K.M_sun / K.Mpc**3)) / K.kpc, 50),
        ("03", "im Abstand von 1 Mpc bis 50 kpc statt bis 3,1 Mpc", K.groesste_umkehr(6e10 * K.M_sun)[0] / K.Mpc, 3.1),
        ("03", "verliert ihren Tiefbereich innerhalb von 0,33 Mpc", (3 * K.G * 1e15 * K.M_sun * K.v_flach(6e10 * K.M_sun)**2 / K.a0**2)**(1 / 3) / K.Mpc, 0.33),
        ("03", "verliert ihren Tiefbereich innerhalb von 0,33 Mpc", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda D: K.a0 / (math.sqrt(3) * math.sqrt(K.G * 1e15 * K.M_sun / D**3) * K.v_flach(6e10 * K.M_sun)) - 1, 0.01 * K.Mpc, 10 * K.Mpc) / K.Mpc, 0.33),
        ("05", "t = (4 − d)·E_kin/((l·(d − 2) + 2)·P)", (lambda d, l: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda r: ((0.5*(r*1.000001)**(2 - d) + (math.log(r*1.000001) if d == 2 else (r*1.000001)**(2 - d)/(2 - d))) - (0.5*(r*0.999999)**(2 - d) + (math.log(r*0.999999) if d == 2 else (r*0.999999)**(2 - d)/(2 - d))))/(2e-6*r)/(r**(2*l)*r**(-d*(l + 1))), 0, 1, limit=200)[0]/K.zeit_abgabe(l, 0.5, 1.0, d=d))(2.5, 2), 1),
        ("05", "t = (4 − d)·E_kin/((l·(d − 2) + 2)·P)", (lambda d, l: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda r: ((0.5*(r*1.000001)**(2 - d) + (math.log(r*1.000001) if d == 2 else (r*1.000001)**(2 - d)/(2 - d))) - (0.5*(r*0.999999)**(2 - d) + (math.log(r*0.999999) if d == 2 else (r*0.999999)**(2 - d)/(2 - d))))/(2e-6*r)/(r**(2*l)*r**(-d*(l + 1))), 0, 1, limit=200)[0]/K.zeit_abgabe(l, 0.5, 1.0, d=d))(3.5, 1), 1),
        ("05", "Nach Newton, d = 3, ist das E/((l + 2)·P)", K.zeit_abgabe(2, 0.5, 1.0, d=3)/K.zeit_abgabe(2, 0.5, 1.0), 1),
        ("05", "Nach Newton, d = 3, ist das E/((l + 2)·P)", K.zeit_abgabe(1, 0.5, 1.0, d=3)/K.zeit_abgabe(1, 0.5, 1.0), 1),
        ("05", "je Faktor m·v_flach²", (lambda M: (lambda r: (0.5*K.bahntempo_regel(M, math.e*r)**2 - 0.5*K.bahntempo_regel(M, r)**2 + __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.bahntempo_regel(M, x)**2/x, r, math.e*r, limit=200)[0])/K.v_flach(M)**2)(1e3*K.schwellenradius(M)))(6e10*K.M_sun), 1),
        ("05", "jede Bahn fällt in t = E_kin/P", (lambda m, v, r0: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda r: m*v**2/r/K.leistung_quadrupol(m, v**2*r/K.G, r), 0, r0, limit=200)[0]/K.zeit_abgabe(2, 0.5*m*v**2, K.leistung_quadrupol(m, v**2*r0/K.G, r0), d=2))(1e30, 2e5, 3e19), 1),
        ("05", "jede Bahn fällt in t = E_kin/P", (lambda d, l: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda r: ((0.5*(r*1.000001)**(2 - d) + (math.log(r*1.000001) if d == 2 else (r*1.000001)**(2 - d)/(2 - d))) - (0.5*(r*0.999999)**(2 - d) + (math.log(r*0.999999) if d == 2 else (r*0.999999)**(2 - d)/(2 - d))))/(2e-6*r)/(r**(2*l)*r**(-d*(l + 1))), 0, 1, limit=200)[0]/K.zeit_abgabe(l, 0.5, 1.0, d=d))(2, 1), 1),
        ("F:lesungen", "bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache", _uebergabe()[0], 3.5),
        ("F:lesungen", "bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache", _uebergabe()[1], 2.4),
        ("F:lesungen", "bei 10¹¹ M☉ das 3,5-fache, bei 10¹² M☉ das 2,4-fache, bei 10¹³ M☉ das 1,7-fache", _uebergabe()[2], 1.7),
        ("F:lesungen", "sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']))[0] - D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), gas=True)[0], 6),
        ("F:lesungen", "sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(1.0, None), gas=True)[5], 0.69),
        ("F:lesungen", "sinkt χ² nur um 6, und der Rest liest weiter 0,69 ± 0,06", D.kids_scherung_chi2((lambda p: (p['wb'], p['wc'], p['H0'], p['lnAs'], p['ns'], p['w0'], p['wa']))(D.planck_anpassung(K.q_je_radiant(), ('spektren', 'linse', 'bao', 'sn', 'wachstum', 'bbn'), zustand='saettigung', start={'wa': -0.56})['punkt']), teile=(1.0, None), gas=True)[7], 0.06),
        ("03", "a₀ = 1,10·10⁻¹⁰ m/s² (1,09 bis", _bao_a0()[1], 1.09),
        ("03", "hier (1 + Q)² = 40,5-mal Chandrasekhar", K.reibung_haufen(6e10 * K.M_sun, 50 * K.kpc)[0], 40.5),
        ("03", "hier (1 + Q)² = 40,5-mal Chandrasekhar", sum((n + 1) * (K.Q_kosm / (1 + K.Q_kosm))**n for n in range(2000)), 40.5),
        ("03", "reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.verstaerkung(K.G * 6e10 * K.M_sun / r**2) * K.G * 6e10 * K.M_sun / r**2 - 3 * K.G * 1e15 * K.M_sun / K.Mpc**3 * r, K.kpc, K.Mpc) / K.kpc, 49.6),
        ("03", "reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc", (lambda Om: (lambda f: f(f(f(f(f(f(f(f(f(f(f(f(50 * K.kpc)))))))))))))(lambda r: K.bahntempo_regel(6e10 * K.M_sun, r) / (math.sqrt(3) * Om)))(math.sqrt(K.G * 1e15 * K.M_sun / K.Mpc**3)) / K.kpc, 49.6),
        ("03", "reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc", K.naht_radius(6e10 * K.M_sun) / K.kpc, 52.0),
        ("03", "reicht bis 49,6 kpc, fast bis an ihre Naht bei 52,0 kpc", K.schwellenradius(6e10 * K.M_sun) / math.log((1 + K.Q_kosm) / K.Q_kosm) / K.kpc, 52.0),
        ("03", "und das Verhältnis ist 1,04", K.reibung_haufen(6e10 * K.M_sun, __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.verstaerkung(K.G * 6e10 * K.M_sun / r**2) * K.G * 6e10 * K.M_sun / r**2 - 3 * K.G * 1e15 * K.M_sun / K.Mpc**3 * r, K.kpc, K.Mpc))[1], 1.04),
        ("03", "und das Verhältnis ist 1,04", (1 + K.Q_kosm) * -math.expm1(-K.schwellenradius(6e10 * K.M_sun) / (lambda Om: (lambda f: f(f(f(f(f(f(f(f(f(f(f(f(50 * K.kpc)))))))))))))(lambda r: K.bahntempo_regel(6e10 * K.M_sun, r) / (math.sqrt(3) * Om)))(math.sqrt(K.G * 1e15 * K.M_sun / K.Mpc**3))), 1.04),
        ("03", "ist es (1 + Q)·(1 − 1/e) = 4,0", K.reibung_haufen(6e10 * K.M_sun, K.schwellenradius(6e10 * K.M_sun))[1], 4.0),
        ("03", "ist es (1 + Q)·(1 − 1/e) = 4,0", (1 + K.Q_kosm) * (1 - 1 / math.e), 4.0),
        ("03", "mit γ = 5/3 und Sternen von einem Viertel des Gases, ist q = 0,024", K.welle_haufen(0.8, 1, 2)[0], 0.024),
        ("03", "mit γ = 5/3 und Sternen von einem Viertel des Gases, ist q = 0,024", (lambda vc, r: 4 * math.pi * K.G * (1 + K.Q_kosm) * 0.8 * vc**2 / (4 * math.pi * K.G * (1 + K.Q_kosm) * r**2) / (5 / 3 * vc**2 / 2 * (2 * math.pi / r)**2))(1e6, K.Mpc), 0.024),
        ("03", "und das Gas verstärkt sie um 2,5 %", 100 * (K.welle_haufen(0.8, 1, 2)[1] - 1), 2.5),
        ("03", "und das Gas verstärkt sie um 2,5 %", 100 * (1 / (1 - (lambda vc, r: 4 * math.pi * K.G * (1 + K.Q_kosm) * 0.8 * vc**2 / (4 * math.pi * K.G * (1 + K.Q_kosm) * r**2) / (5 / 3 * vc**2 / 2 * (2 * math.pi / r)**2))(1e6, K.Mpc)) - 1), 2.5),
        ("03", "schaukelt sie erst bei γ_P = 4π²·γ/3 = 22 auf", K.welle_haufen(0.8, 1, 2)[2], 22),
        ("03", "schaukelt sie erst bei γ_P = 4π²·γ/3 = 22 auf", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda gP: K.welle_haufen(1.0, 3, gP)[0] - 1, 1, 100), 22),
        ("02", "die Erde bis 1,50·10⁶ km", K.glaettungslaenge(K.M_earth, 3 * K.GM_sun / K.AU**3) / 1e9, 1.50),
        ("02", "die Erde bis 1,50·10⁶ km", K.AU * (K.M_earth / (3 * K.M_sun))**(1 / 3) / 1e9, 1.50),
        ("02", "bei der mittleren Dichte der Erde 2,4 Stunden", K.zerreiss_umlauf(K.M_earth / (4 / 3 * math.pi * K.R_ERDE**3)) / 3600, 2.4),
        ("02", "bei der mittleren Dichte der Erde 2,4 Stunden", 2 * math.pi / __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda Om: K.gezeitendichte(3 * Om**2) - K.M_earth / (4 / 3 * math.pi * K.R_ERDE**3), 1e-6, 1e-2) / 3600, 2.4),
        ("02", "die kürzer läuft als 4,8 Stunden", K.zerreiss_umlauf(K.M_sun / (4 / 3 * math.pi * K.R_SONNE**3)) / 3600, 4.8),
        ("02", "die kürzer läuft als 4,8 Stunden", 2 * math.pi / __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda Om: K.gezeitendichte(3 * Om**2) - K.M_sun / (4 / 3 * math.pi * K.R_SONNE**3), 1e-6, 1e-2) / 3600, 4.8),
        ("02", "um jeden Planeten 3·√(π/(Gρ))", K.zerreiss_umlauf(1000.0) / (3 * math.sqrt(math.pi / (K.G * 1000.0))), 1),
        ("03", "am Ort der Sonne (8,25 kpc) 1,7 Wasserstoffatome je cm³", (lambda M, D: K.gezeitendichte(K.gezeit_quelle(M, D)[0] - K.verstaerkung(K.G * M / D**2) * K.G * M / D**3) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, 8.25 * K.kpc), 1.7),
        ("03", "am Ort der Sonne (8,25 kpc) 1,7 Wasserstoffatome je cm³", (lambda M, D: K.gezeitendichte(-(K.bahntempo_regel(M, D * 1.0001)**2 / (D * 1.0001) - K.bahntempo_regel(M, D * 0.9999)**2 / (D * 0.9999)) / (0.0002 * D)) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, 8.25 * K.kpc), 1.7),
        ("03", "bei 2R_d 0,57, bei 20 kpc 0,18", (lambda M, D: K.gezeitendichte(K.gezeit_quelle(M, D)[0] - K.verstaerkung(K.G * M / D**2) * K.G * M / D**3) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, math.sqrt(2) * K.schwellenradius(6e10 * K.M_sun)), 0.57),
        ("03", "bei 2R_d 0,57, bei 20 kpc 0,18", (lambda M, D: K.gezeitendichte(-(K.bahntempo_regel(M, D * 1.0001)**2 / (D * 1.0001) - K.bahntempo_regel(M, D * 0.9999)**2 / (D * 0.9999)) / (0.0002 * D)) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, math.sqrt(2) * K.schwellenradius(6e10 * K.M_sun)), 0.57),
        ("03", "bei 2R_d 0,57, bei 20 kpc 0,18", (lambda M, D: K.gezeitendichte(K.gezeit_quelle(M, D)[0] - K.verstaerkung(K.G * M / D**2) * K.G * M / D**3) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, 20 * K.kpc), 0.18),
        ("03", "bei 2R_d 0,57, bei 20 kpc 0,18", (lambda M, D: K.gezeitendichte(-(K.bahntempo_regel(M, D * 1.0001)**2 / (D * 1.0001) - K.bahntempo_regel(M, D * 0.9999)**2 / (D * 0.9999)) / (0.0002 * D)) * (1 - K.YP_GEMESSEN[0]) / K.m_p / 1e6)(6e10 * K.M_sun, 20 * K.kpc), 0.18),
        ("F:lesungen", "Bei 10° ist das Achsenverhältnis der Sterne, 0,86, keine Neigung", D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt'))[1], 0.86),
        ("F:lesungen", "mindestens um 13 % oval", (lambda z, i: 100 * (1 - z[1] / math.sqrt(math.cos(i)**2 + (0.5 * (z[3] + z[4]))**2 * math.sin(i)**2)))(D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt')), math.radians(D.udg_neigung()[3])), 13),
        ("F:lesungen", "mindestens um 13 % oval", (lambda z, i: 100 * (1 - __import__('scipy.optimize', fromlist=['brentq']).brentq(lambda p: (lambda R: (lambda B: (lambda M: math.sqrt(min(D.np.linalg.eigvalsh(M)) / max(D.np.linalg.eigvalsh(M))))(B[:2, :2] - D.np.outer(B[:2, 2], B[2, :2]) / B[2, 2]))(R @ D.np.diag([1, 1 / p**2, 1 / (0.5 * (z[3] + z[4]))**2]) @ R.T))(D.np.array([[1, 0, 0], [0, math.cos(i), -math.sin(i)], [0, math.sin(i), math.cos(i)]])) - z[1], 0.5, 1.0)))(D.np.loadtxt(D.os.path.join(D.ORT, 'linsen_haufen', 'udg_agc114905_sterne2024.txt')), math.radians(D.udg_neigung()[3])), 13),
        ("F:lesungen", "das 2,1-fache ihrer heutigen Baryonen einmal gebunden", D.kids_typ_gedaechtnis("rot")[0] - 1, 2.1),
        ("F:lesungen", "für jede Masse höchstens das 0,3-fache ihrer Baryonen dazu", (lambda x: x / (1 - x))(4 / 3 * math.pi * (1 + K.q_je_radiant()) * K.G * K.dichte_heute()[1] / K.groesste_umkehr(K.M_sun)[1]**2), 0.3),
        ("F:lesungen", "für jede Masse höchstens das 0,3-fache ihrer Baryonen dazu", (lambda M, rho: (lambda R: 4 / 3 * math.pi * R**3 * rho / M)(__import__('scipy.optimize', fromlist=['brentq']).brentq(lambda R: (1 + K.q_je_radiant()) * K.G * (M + 4 / 3 * math.pi * R**3 * rho) / R**2 - K.groesste_umkehr(M)[1]**2 * R, 0.01 * K.Mpc, 30 * K.Mpc)))(10**10.5 * K.M_sun, K.omega_b * 3 * (1e5 / K.Mpc)**2 / (8 * math.pi * K.G)), 0.3),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", math.sqrt(2 / K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 18.0, n=1000, quer="scherung")["quer"] - 2), 2.9),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 18.0, n=1000, quer=0.2)["DL"] / (lambda kP: math.sqrt(__import__("scipy.integrate", fromlist=["quad"]).quad(lambda lk: (lambda k: D.np.interp(math.log(k), D.np.log(kP[0]), kP[1]) * (3 * (math.sin(k * 18.0) - k * 18.0 * math.cos(k * 18.0)) / (k * 18.0)**3)**2 * k**3)(math.exp(lk)), math.log(1e-4), math.log(10.0), limit=400)[0] / (2 * math.pi**2)))(K.spektrum_linear((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748))), 2.9),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", math.sqrt(2 / K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=1000, quer="scherung")["quer"] - 2), 3.8),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=1000, quer=0.2)["DL"] / (lambda kP: math.sqrt(__import__("scipy.integrate", fromlist=["quad"]).quad(lambda lk: (lambda k: D.np.interp(math.log(k), D.np.log(kP[0]), kP[1]) * (3 * (math.sin(k * 25.0) - k * 25.0 * math.cos(k * 25.0)) / (k * 25.0)**3)**2 * k**3)(math.exp(lk)), math.log(1e-4), math.log(10.0), limit=400)[0] / (2 * math.pi**2)))(K.spektrum_linear((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748))), 3.8),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 18.0, quer="scherung")), 68.24915)[0]), 4.4),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 18.0, n=1500, quer=2 / ((K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 18.0, n=1000, quer=0.2)["DL"] / (lambda kP: math.sqrt(__import__("scipy.integrate", fromlist=["quad"]).quad(lambda lk: (lambda k: D.np.interp(math.log(k), D.np.log(kP[0]), kP[1]) * (3 * (math.sin(k * 18.0) - k * 18.0 * math.cos(k * 18.0)) / (k * 18.0)**3)**2 * k**3)(math.exp(lk)), math.log(1e-4), math.log(10.0), limit=400)[0] / (2 * math.pi**2)))(K.spektrum_linear((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748))))**2 + 2))), 68.24915)[0]), 4.4),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, quer="scherung")), 68.24915)[0]), 5.3),
        ("03", "(ν = 2,9 und 3,8) c = 4,4 und 5,3", float(K.konzentration_profil(*(lambda e: (e["r"], e["M"]))(K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=1500, quer=2 / ((K.einfall_haufen((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748), 25.0, n=1000, quer=0.2)["DL"] / (lambda kP: math.sqrt(__import__("scipy.integrate", fromlist=["quad"]).quad(lambda lk: (lambda k: D.np.interp(math.log(k), D.np.log(kP[0]), kP[1]) * (3 * (math.sin(k * 25.0) - k * 25.0 * math.cos(k * 25.0)) / (k * 25.0)**3)**2 * k**3)(math.exp(lk)), math.log(1e-4), math.log(10.0), limit=400)[0] / (2 * math.pi**2)))(K.spektrum_linear((0.02243, 0.1185, 68.24915, 3.05217, 0.96881, -0.88719, -0.46748))))**2 + 2))), 68.24915)[0]), 5.3),
        ("F:lesungen", "(1 − 1/2π)² = 0,71", (K.q_je_radiant() / (1 + K.q_je_radiant()))**2, 0.71),
        ("F:lesungen", "(1 − 1/2π)² = 0,71", (1 - 1 / (2 * math.pi))**2, 0.71),
        ("F:lesungen", "ist χ² = 275 statt 307,8", _kids_teile(1.0, (K.q_je_radiant() / (1 + K.q_je_radiant()))**2)[0], 275),
        ("F:lesungen", "ist χ² = 275 statt 307,8", _kids_teile(1.0, None)[0] + (((K.q_je_radiant() / (1 + K.q_je_radiant()))**2 - _kids_teile(1.0, None)[5]) / _kids_teile(1.0, None)[7])**2, 275),
        ("F:lesungen", "ist χ² = 275 statt 307,8", _kids_teile(1.0, 1.0)[0], 307.8),
        ("F:lesungen", "Frei liest der Rest dann 0,67 ± 0,06", _kids_teile(1.0, None)[5], 0.67),
        ("F:lesungen", "Frei liest der Rest dann 0,67 ± 0,06", (lambda t: t[5] + t[8] * t[7] / t[6] * (1 - t[4]))(_kids_teile(None, None)), 0.67),
        ("F:lesungen", "Frei liest der Rest dann 0,67 ± 0,06", _kids_teile(1.0, None)[7], 0.06),
        ("F:lesungen", "Frei liest der Rest dann 0,67 ± 0,06", (lambda t: t[7] * math.sqrt(1 - t[8]**2))(_kids_teile(None, None)), 0.06),
        ("F:lesungen", "q² liegt Δχ² = 0,4 daneben", (((K.q_je_radiant() / (1 + K.q_je_radiant()))**2 - _kids_teile(1.0, None)[5]) / _kids_teile(1.0, None)[7])**2, 0.4),
        ("F:lesungen", "q² liegt Δχ² = 0,4 daneben", _kids_teile(1.0, (K.q_je_radiant() / (1 + K.q_je_radiant()))**2)[0] - _kids_teile(1.0, None)[0], 0.4),
        ("F:galaxien", "(Rangkorrelation −0,75)", (lambda t: __import__("scipy.stats", fromlist=["spearmanr"]).spearmanr(t[0], t[1])[0])((lambda yv: (lambda H: (lambda sr, so: (sr, D.np.log10(so / sr)))(D.np.array([K.streuung_im_feld(yv * z["L"] * K.M_sun, z["R_e"], 0.0) / 1e3 for z in H]), D.np.array([z["sigma"] for z in H])))([z for z in D.kugelhaufen() if z["lg_g"] < 0]))(D.kugelhaufen_reste()[0])), -0.75),
        ("F:galaxien", "(Rangkorrelation −0,75)", (lambda t: D.np.corrcoef(D.np.argsort(D.np.argsort(t[0])), D.np.argsort(D.np.argsort(t[1])))[0, 1])((lambda yv: (lambda H: (lambda sr, so: (sr, D.np.log10(so / sr)))(D.np.array([K.streuung_im_feld(yv * z["L"] * K.M_sun, z["R_e"], 0.0) / 1e3 for z in H]), D.np.array([z["sigma"] for z in H])))([z for z in D.kugelhaufen() if z["lg_g"] < 0]))(D.kugelhaufen_reste()[0])), -0.75),
        ("F:galaxien", "an den dichten mitgeeicht (Υ_V = 0,8)", D.kugelhaufen_reste(zusatz=1.5)[0], 0.8),
        ("F:galaxien", "an den dichten mitgeeicht (Υ_V = 0,8)", __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda y: sum((1 / ((z["fehler"] / z["sigma"] / math.log(10))**2 + 0.01)) * math.log10(z["sigma"] / math.hypot(K.streuung_im_feld(y * z["L"] * K.M_sun, z["R_e"], 0.0) / 1e3, 1.5))**2 for z in D.kugelhaufen() if z["lg_g"] > math.log10(3)), bounds=(0.1, 5), method="bounded").x, 0.8),
        ("F:galaxien", "die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33)", D.kugelhaufen_reste(zusatz=1.5)[3], 0.00),
        ("F:galaxien", "die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33)", D.kugelhaufen_reste(zusatz=1.5)[4], 0.03),
        ("F:galaxien", "die dünnen liegen dann 0,00 ± 0,03 dex auf der Regel (χ² 33)", D.kugelhaufen_reste(zusatz=1.5)[5], 33),
        ("F:galaxien", "Newton bleibt damit bei +0,12 und verlangt rund 2 km/s", D.kugelhaufen_reste(True, zusatz=1.5)[3], 0.12),
        ("F:galaxien", "Newton bleibt damit bei +0,12 und verlangt rund 2 km/s", __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda s: D.kugelhaufen_reste(True, zusatz=s)[5], bounds=(0, 5), method="bounded").x, 2),
        ("F:galaxien", "Newton bleibt damit bei +0,12 und verlangt rund 2 km/s", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda s: D.kugelhaufen_reste(True, zusatz=s)[3], 0.5, 2.35), 2),
        ("09", "Ω = 2GJ/(c²r³) = 6Δs·c·a_J/r³", K.mitfuehrung(K.J_ERDE, K.R_GPB, bahnebene=True) / (6 * K.faltung(K.M_earth) * K.c * K.verwindung(K.J_ERDE, K.M_earth) / K.R_GPB**3), 1.0),
        ("09", "3Δs·c·a_J/(2r³), einem Viertel davon", K.mitfuehrung(K.J_ERDE, K.R_GPB) / K.mitfuehrung(K.J_ERDE, K.R_GPB, bahnebene=True), 0.25),
        ("08", "Die Temperatur ħc·(r₊ − r₋)/(4πk·r₀·r₊)", K.falte_gedreht(K.M_sun, 0.0)["temperatur"] / K.falte_temperatur(K.M_sun), 1.0),
        ("08", "Erst bei M_U/4 = 1,4·10²³ M☉", K.schnittmassen()[3] / 4 / K.M_sun / 1e23, 1.4),
        ("08", "Am Horizont ist s² = M_U/(4M)", K.c**4 / (4 * K.G * K.M_sun) / K.a0 / (K.schnittmassen()[3] / (4 * K.M_sun)), 1.0),
        ("08", "kT = m̃c²/(4π). Strahlt", K.falte_temperatur(K.M_sun) * 24 * math.pi * K.k_B * K.faltung(K.M_sun) / (K.hbar * K.c), 1.0),
        ("03", "(r_d = 147,1 Mpc aus ω_b und ω_c)", K.schallhorizont(), 147.1),
        ("03", "liest DESI dann a₀ = 1,10", _bao_a0()[0], 1.10),
        ("03", "m/s² (1,09 bis 1,11", _bao_a0()[2], 1.11),
        ("03", "mit H₀ = 68,9 bei χ² = 10,4", K.atmen_heute(_bao_a0()[0] * 1e-10)[0], 68.9),
        ("03", "bei χ² = 10,4 für 13 Werte", _bao_a0()[3], 10.4),
        ("03", "für 13 Werte", float(len(D.desi()[0])), 13),
        ("03", "an sechs Rotverschiebungen", float(len(B.kosmos()["zeilen"])), 6),
        ("01", "(PubChem, 13,598 eV)", B.atom_und_kern()["IE_H"], 13.598),
        ("03", "allein (ω_b, ω_c nach Planck), ist a₀ = 1,09", D.abstaende_a0(_om_or_a0)[0], 1.09),
        ("03", "(1,07 bis 1,11 bei Δχ² = 1)", D.abstaende_a0(_om_or_a0)[1], 1.07),
        ("03", "1,07 bis 1,11 bei Δχ²", D.abstaende_a0(_om_or_a0)[2], 1.11),
        ("03", "mit H₀ = 68,4: Rotationskurven", K.atmen_heute(D.abstaende_a0(_om_or_a0)[0] * 1e-10)[0], 68.4),
        ("02", "nach Z = 2, 10, 18, 36, 54, 86", float(D.schliesszahlen()[2] == (2, 10, 18, 36, 54, 86)), 1.0),
        ("02", "bei N = 16, 20, 28, 50, 82, 126", float(D.schliesszahlen()[0] == (16, 20, 28, 50, 82, 126)), 1.0),
        ("02", "Z = 14, 20, 28, 38, 50, 82 (AME", float(D.schliesszahlen()[1] == (14, 20, 28, 38, 50, 82)), 1.0),
        ("02", "Faktor 4,6 bei He", D.schliesszahlen()[3], 4.6),
        ("02", "Halterlänge λ_H = 1,25 fm", KE.halterlaenge(KE.a_C(KE.TROPFEN)), 1.25),
        ("02", "v/c ≈ λ̄_p/λ_H = 0,17", K.hbar / (K.m_p * K.c) * 1e15 / KE.halterlaenge(KE.a_C(KE.TROPFEN)), 0.17),
        ("02", "Formel zu λ_H = 1,33 fm", KE.halterlaenge(D.spiegelkerne()[1]), 1.33),
        ("02", "unabhängig davon 0,14", math.sqrt(2 * KE.bekannte()["b_je_a"] / (K.energie(K.m_u) / K.e / 1e6)), 0.14),
        ("02", "angepasst an die 2509 gemessenen Kerne", float(KE.anpassung(KE.TROPFEN)[2]), 2509),
        ("02", "AME2020: x = 0,25", KE.bekannte()["x_max"], 0.25),
        ("02", "(Z = 118, A = 294: x = 0,98)", KE.bekannte()["ende"][2], 0.98),
        ("02", "(Z = 118, A = 294", float(KE.bekannte()["ende"][:2] == (118, 294)), 1.0),
        ("02", "Die Hülle hält bis Z = 170 bis 171", min(KE.grenzen()["huelle"]), 170),
        ("02", "Hülle hält bis Z = 170 bis 171:", max(KE.grenzen()["huelle"]), 171),
        ("02", "(36 Paare, AME2020)", float(D.spiegelkerne()[0]), 36),
        ("02", "Zα = 1 (Z = 137) gilt für einen punktförmigen Kern", math.floor(1 / K.alpha), 137),
        ("02", "Mit dem Edelgas Og schließt die siebte Periode bei Z = 118 (PubChem)", float(KE.konfigurationen()[0] == (2, 10, 18, 36, 54, 86, 118) == K.schliesszahlen_huelle(7)), 1.0),
        ("02", "Die Perioden sind 2, 8, 8, 18, 18, 32, 32 lang", float(all(b - a == 2 * math.ceil((p + 1) / 2)**2 for p, (a, b) in enumerate(zip((0,) + KE.konfigurationen()[0], KE.konfigurationen()[0]), 1))), 1.0),
        ("02", "bindet an der Schwelle E = 0 genau bei den Kopplungen N(N + 1)", float(max(abs(v - (n + l) * (n + l + 1)) for (n, l), v in __import__("huelle_4_formen").schliessform().items()) < 0.01), 1.0),
        ("02", "Sein Feld folgt der Spiegelform innen auf 0,6 %", 100 * _ght()["innen"], 0.6),
        ("02", "Seine Höhe v = 3,44·Z^{2/3} liegt von Ne bis Og im Fenster", float(round(_ght()["k"], 2) == 3.44 and all(e["N"] == e["p"] + 1 for Z, e in _ght()["edel"].items() if Z >= 10)), 1.0),
        ("02", "Seine Zustände an der Schwelle sind Kugelfunktionen, u = √r·P_N^{2l+1}((r − R)/(r + R))", float(_kug()["rest"] < 1e-10 and _kug()["ueberlapp"] > 0.9999), 1.0),
        ("01", "Spiegelfeld der Hülle (Kugelfunktionen auf S², Grad n + l, Ordnung 2l + 1)", float(_kug()["rest"] < 1e-10 and _kug()["ueberlapp"] > 0.9999), 1.0),
        ("02", "erst bei Zα = 0,86", 118 * K.alpha, 0.86),
        ("02", "(2n² je Schale: 2, 10, 28, 60)", float(tuple(sum(2 * n * n for n in range(1, k + 1)) for k in range(1, 5)) == (2, 10, 28, 60)), 1.0),
        ("02", "nach 2, 8, 20, 40, 70, 112", float(tuple(sum((n + 1) * (n + 2) for n in range(k)) for k in range(1, 7)) == (2, 8, 20, 40, 70, 112)), 1.0),
        ("02", "also nach 28, 50, 82, 126 (k = 2 bis 5)", float(tuple((k + 1) * (k + 2) * (k + 3) // 3 + 2 * (k + 2) for k in range(2, 6)) == (28, 50, 82, 126) and {28, 50, 82, 126} <= set(D.schliesszahlen()[0])), 1.0),
        ("02", "Die nächste Schließung wäre 184", (6 + 1) * (6 + 2) * (6 + 3) // 3 + 2 * (6 + 2), 184),
        ("02", "schließen nach 2, 10, 18, 36, 54, 86 (PubChem)", float(K.schliesszahlen_huelle()[:6] == D.schliesszahlen()[2]), 1.0),
        ("02", "2,8 bei Rn", D.schliesszahlen()[4], 2.8),
        ("03", "40 km/s 72 % der", D.sparc_gasbruch(D.sparc_eichung(K.kpc), 0, 40) * 100, 72),
        ("03", "150 km/s 12 %", D.sparc_gasbruch(D.sparc_eichung(K.kpc), 150, 1e9) * 100, 12),
        ("03", "Schwelle zu a₀ = 1,11·10⁻¹⁰", D.sparc_schwelle(K.kpc, yd=D.sparc_eichung(K.kpc))[0] / 1e-10, 1.11),
        ("03", "(16 bis 84 %: 1,00 bis 1,21)", D.sparc_schwelle(K.kpc, yd=D.sparc_eichung(K.kpc))[1][0] / 1e-10, 1.00),
        ("03", "1,00 bis 1,21). Die Bose", D.sparc_schwelle(K.kpc, yd=D.sparc_eichung(K.kpc))[1][1] / 1e-10, 1.21),
        ("03", "die einfache (Δχ² = 53)", D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["einfach"]
            - D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["Bose"], 53),
        ("03", "Standard-Form (Δχ² = 41)", D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["standard"]
            - D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["Bose"], 41),
        ("07", "Dimension 0,92 ± 0,05", D.sparc_dimension(K.kpc, yd=D.sparc_eichung(K.kpc))[0], 0.92),
        ("03", "nach Planck a₀ = 1,065·10⁻¹⁰", K.a0_aus_atmen() / 1e-10, 1.065),
        ("03", "m/s² H₀ = 67,2 km/s/Mpc", K.atmen_heute(1.06e-10)[0], 67.2),
        ("03", "und Ω_DE = 0,685, mit", K.atmen_heute(1.06e-10)[1], 0.685),
        ("03", "1,1·10⁻¹⁰ H₀ = 68,9", K.atmen_heute(1.1e-10)[0], 68.9),
        ("03", "68,9 und 0,700", K.atmen_heute(1.1e-10)[1], 0.700),
        ("03", "erst bei z = −0,10 ein", ((K.vakuumteil() / (K.Omega_b + K.Omega_c))**(1 / 3)) - 1, -0.1),
        ("03", "= 0,47·v_flach", K.haltegrenze(1.0), 0.47),
        ("03", "unter v_flach = 17 bis 21", 8e3 / K.haltegrenze(1e3), 17),
        ("03", "17 bis 21 km/s. In SPARC", 10e3 / K.haltegrenze(1e3), 21),
        ("03", "bei z_eq = 3402", K.z_gleichheit(), 3402),
        ("03", "bei 3 bis 4·10¹²", K.grenze_verschmelzen("EH") / 1e12, 3),
        ("03", "3 bis 4·10¹² M☉ (mit", K.grenze_verschmelzen("BBKS") / 1e12, 4),
        ("03", "Bereiche ab 6,4·10¹⁰", K.bindungsmasse(100e3, traegt_selbst=False) / K.M_sun / 1e10, 6.4),
        ("03", "Bereiche ab 6,4·10¹⁰", (lambda s: (lambda M: K.bahntempo_regel(M, 1e4 * K.schwellenradius(M)) / math.sqrt(3) / s)(9 * s**4 / (K.G * K.a0)))(100e3) * 6.4, 6.4),
        ("03", "bis 5,2·10¹² M☉", K.bindungsmasse(300e3, traegt_selbst=False) / K.M_sun / 1e12, 5.2),
        ("03", "bei M > 9·σ⁴/(Ga₀), trägt es sein Feld selbst, bei M > (81/4)·σ⁴/(Ga₀)", K.bindungsmasse(1e5) / K.bindungsmasse(1e5, traegt_selbst=False) / (81 / 36), 1),
        ("03", "M_Gas = 0,6 bis 1,4·10⁷ M☉", K.bindungsmasse(8e3) / K.M_sun / 1e7, 0.6),
        ("03", "1,4·10⁷ M☉. Kaltes", K.bindungsmasse(10e3) / K.M_sun / 1e7, 1.4),
        ("03", "zum Endwert 1,23", 2 * math.pi * K.a0 / (K.c * math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie())), 1.23),
        ("03", "drei Teile tun es bei z = 0,30", (K.dichte_dunkle_energie() / K.dichte_heute()[2])**(1 / 3) - 1, 0.3),
        ("03", "Mit drei Teilen in H_∞ = c·√Λ_s = 55,6 km/s/Mpc", math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) * K.Mpc / 1e3, 55.6),
        ("03", "9,1 km/s bei 10⁷ M☉", K.streuung_kugel(1e7 * K.M_sun) / 1e3, 9.1),
        ("03", "19,3 km/s bei 2·10⁸ M☉", K.streuung_kugel(2e8 * K.M_sun) / 1e3, 19.3),
        ("03", "10¹² M☉/Mpc: 2,6 kpc", K.filament_schwelle(1e12 * K.M_sun / K.Mpc) / K.kpc, 2.6),
        ("03", "(ν = 1,007)", K.verstaerkung(25 * K.a0), 1.007),
        ("03", "(ν = 3,3)", K.verstaerkung(0.13 * K.a0), 3.3),
        ("F:galaxien", "Mit Υ_V = 0,91, geeicht", D.kugelhaufen_reste()[0], 0.91),
        ("F:galaxien", "an den dichten (g_N > 3·a₀, 28)", D.kugelhaufen_reste()[1], 28),
        ("F:galaxien", "dünnen (g_N < a₀, 17)", D.kugelhaufen_reste()[2], 17),
        ("F:galaxien", "mit der Regel +0,10 ± 0,03 dex über", D.kugelhaufen_reste()[3], 0.10),
        ("F:galaxien", "mit der Regel +0,10 ± 0,03 dex über", D.kugelhaufen_reste()[4], 0.03),
        ("F:galaxien", "ohne Halo, +0,41 (", D.kugelhaufen_reste(True)[3], 0.41),
        ("F:galaxien", "(Υ_V = 1,08, χ² 54", D.kugelhaufen_reste(True)[0], 1.08),
        ("F:galaxien", "(Υ_V = 1,08, χ² 54", D.kugelhaufen_reste()[5], 54),
        ("F:galaxien", "χ² 54 gegen 235)", D.kugelhaufen_reste(True)[5], 235),
        ("F:galaxien", "Die gut gemessenen (11)", D.kugelhaufen_reste(False, "gut")[2], 11),
        ("F:galaxien", "liegen +0,03 gegen +0,34", D.kugelhaufen_reste(False, "gut")[3], 0.03),
        ("F:galaxien", "liegen +0,03 gegen +0,34", D.kugelhaufen_reste(True, "gut")[3], 0.34),
        ("F:galaxien", "10 kpc vom Zentrum (9) +0,06", D.kugelhaufen_reste(False, "halo")[2], 9),
        ("F:galaxien", "10 kpc vom Zentrum (9) +0,06", D.kugelhaufen_reste(False, "halo")[3], 0.06),
        ("F:galaxien", "(9) +0,06 gegen +0,43", D.kugelhaufen_reste(True, "halo")[3], 0.43),
        ("F:galaxien", "Database (52 mit gemessener", len(D.kugelhaufen()), 52),
        ("F:galaxien", "v_esc(8,3 kpc) = 528 (+24, −25) km/s", D.gaia_flucht()["v"], 528),
        ("F:galaxien", "v_esc(8,3 kpc) = 528 (+24, −25) km/s", D.gaia_flucht()["oben"], 24),
        ("F:galaxien", "v_esc(8,3 kpc) = 528 (+24, −25) km/s", D.gaia_flucht()["unten"], 25),
        ("F:galaxien", "v_esc(8,3 kpc) = 528 (+24, −25) km/s", D.gaia_flucht()["r0"], 8.3),
        ("F:galaxien", "bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter", D.gaia_flucht()["rand"] / K.Mpc, 0.42),
        ("F:galaxien", "bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter", K.flucht_regel(6e10 * K.M_sun, D.gaia_flucht()["r0"] * K.kpc, rand=D.gaia_flucht()["rand"]) / 1e3, 511),
        ("F:galaxien", "bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter", (lambda M, r, R: math.sqrt(2 * __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.verstaerkung(K.G * M / x**2) * K.G * M / x**2 - K.groesste_umkehr(M)[1]**2 * x, r, R, limit=200)[0]) / 1e3)(6e10 * K.M_sun, D.gaia_flucht()["r0"] * K.kpc, D.gaia_flucht()["rand"]), 511),
        ("F:galaxien", "bis zu derselben Grenze, 2·r₂₀₀ = 0,42 Mpc, 511 km/s, 0,7 Fehler darunter", (D.gaia_flucht()["v"] - K.flucht_regel(6e10 * K.M_sun, D.gaia_flucht()["r0"] * K.kpc, rand=D.gaia_flucht()["rand"]) / 1e3) / D.gaia_flucht()["unten"], 0.7),
        ("F:galaxien", "Bis zur größten Umkehr gezählt sind es 590 km/s", K.flucht_regel(6e10 * K.M_sun, D.gaia_flucht()["r0"] * K.kpc) / 1e3, 590),
        ("F:galaxien", "für ihre Masse von 10¹² M☉", D.gaia_flucht()["M200"] / 1e12, 1),
        ("F:galaxien", "mit γ = 9/2, der Flussdimension der Kugel", K.isotherme_kugel()[0], 4.5),
        ("F:galaxien", "in jede Richtung, (σ/c)² = (2/9)·(v_flach/c)²", K.streuung_kugel(K.M_sun)**2 / math.sqrt(K.G * K.M_sun * K.a0), 0.222),
        ("F:galaxien", "in jede Richtung, (σ/c)² = (2/9)·(v_flach/c)²", (lambda M: K.streuung_im_feld(M, 1e4 * K.schwellenradius(M), 0.0)**2 / math.sqrt(K.G * M * K.a0))(1e4 * K.M_sun), 0.222),
        ("F:galaxien", "Sternhaufen im Median bei R_h = 0,38·r*", D.stufe_zwerge()[0], 0.38),
        ("F:galaxien", "Sternhaufen im Median bei R_h = 0,38·r*", D.stufe_zwerge(weg="feld")[0], 0.38),
        ("F:galaxien", "38 der 52 innerhalb", D.stufe_zwerge()[1], 38),
        ("F:galaxien", "38 der 52 innerhalb", D.stufe_zwerge(weg="feld")[2], 52),
        ("F:galaxien", "die Zwerge bei 16·r*, 76 der 77 außerhalb", D.stufe_zwerge()[3], 16),
        ("F:galaxien", "die Zwerge bei 16·r*, 76 der 77 außerhalb", D.stufe_zwerge(weg="feld")[3], 16),
        ("F:galaxien", "die Zwerge bei 16·r*, 76 der 77 außerhalb", D.stufe_zwerge()[4], 76),
        ("F:galaxien", "die Zwerge bei 16·r*, 76 der 77 außerhalb", D.stufe_zwerge(weg="feld")[5], 77),
        ("F:galaxien", "(Rangkorrelation −0,51 an 14 Haufen", D.kugelhaufen_gezeit()[1], -0.51),
        ("F:galaxien", "(Rangkorrelation −0,51 an 14 Haufen", D.kugelhaufen_gezeit()[0], 14),
        ("F:galaxien", "am heutigen Abstand zu 0,13 bis 0,64 füllt", D.kugelhaufen_gezeit()[3], 0.13),
        ("F:galaxien", "am heutigen Abstand zu 0,13 bis 0,64 füllt", D.kugelhaufen_gezeit()[4], 0.64),
        ("F:galaxien", "zu 0,13 bis 0,64 füllt (+0,31)", D.kugelhaufen_gezeit()[2], 0.31),
        ("03", "a₀/(πG) = 242 M☉/pc²", K.flaechendichte_schwelle() * K.pc**2 / K.M_sun, 242),
        ("03", "a_N = 2,54·a₀", K.maximum_zusatz()[0], 2.54),
        ("03", "beträgt 0,648·a₀", K.maximum_zusatz()[1], 0.648),
        ("07", "Rosetten mit 254,6°", 360 / K.umlaufphase(1e-300), 254.6),
        ("07", "e^{−7480}", K.schwellenradius(K.M_sun) / K.AU, 7480),
        ("08", "bei der Milchstraße 0,42″", math.degrees(K.ablenkung_grenz(6e10 * K.M_sun)) * 3600, 0.42),
        ("08", "Für die eingeschlossene Masse im Bild ist am Ring deshalb G·M/R_E² = c²/(4D_E), für jede Masse und jede runde Verteilung", (lambda D, M: K.G * M / (4 * K.G * M / K.c**2 * D) / (K.c**2 / (4 * D)))(1.0 * K.Gpc, 1e12 * K.M_sun), 1),
        ("08", "Für die eingeschlossene Masse im Bild ist am Ring deshalb G·M/R_E² = c²/(4D_E), für jede Masse und jede runde Verteilung", (lambda D, s: math.pi * s**2 / (4 * math.pi * s**2 * D / K.c**2) / (K.c**2 / (4 * D)))(1.0 * K.Gpc, 250e3), 1),
        ("08", "Für einen Punkt ist im Newton-Bereich R_E² = 2r₀·D_E wie r*² = r₀·Λ, also (R_E/r*)² = 2D_E/Λ", (lambda D, M: 4 * K.G * M * D / K.c**2 / K.schwellenradius(M)**2 / (2 * D / K.schliesslaenge()))(1.0 * K.Gpc, 1e11 * K.M_sun), 1),
        ("08", "In Flächendichte gelesen sind Ring und Schwelle dieselbe Form, g = πG·Σ̄", (lambda M: K.G * M / K.schwellenradius(M)**2 / (math.pi * K.G * M / (math.pi * K.schwellenradius(M)**2)))(1e11 * K.M_sun), 1),
        ("08", "In Flächendichte gelesen sind Ring und Schwelle dieselbe Form, g = πG·Σ̄", (lambda M: M / (math.pi * K.schwellenradius(M)**2) / (K.a0 / (math.pi * K.G)))(1e9 * K.M_sun), 1),
        ("08", "Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06)", 4 * _ring()[1] / K.Gpc, 5.1),
        ("08", "Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06)", 4 * _ring()[3] / K.Gpc, 5.1),
        ("08", "Auch vor der Hintergrundstrahlung bleibt 4D_E unter 5,1 Gpc (Linse bei z = 1, im Formteil mit a₀ = 1,06)", _ring()[0], 1),
        ("08", "Jeder Ring liegt innerhalb der Schwelle der Masse, die er einschließt, bei G·M/R_E² ≥ 5,3·a₀", _ring()[2], 5.3),
        ("08", "Jeder Ring liegt innerhalb der Schwelle der Masse, die er einschließt, bei G·M/R_E² ≥ 5,3·a₀", K.c**2 / K.a0 / (4 * _ring()[3]), 5.3),
        ("08", "ε_L = ℓ_P/ℓ₀ = 1,9·10⁻⁶²", K.leiterzahl() / 1e-62, 1.9),
        ("08", "ε_L^{1/3}·m_Pl = 5,8·10⁻²⁹ kg", K.schnittmassen()[1] / 1e-29, 5.8),
        ("08", "r̃ ist die Schwelle der Spiegelmasse, r̃ = √2·r*(m̃)", K.schwelle_spiegel(K.M_sun) / (math.sqrt(2) * K.schwellenradius(K.spiegelmasse(K.M_sun))), 1),
        ("08", "und r*/r̃ = M/m_Pl", K.schwellenradius(K.M_earth) / K.schwelle_spiegel(K.M_earth) / (K.M_earth / K.m_Pl), 1),
        ("08", "Beide kreuzen sich bei m_Pl in der Länge des Vakuumteils", K.schwelle_spiegel(K.m_Pl) / math.sqrt(K.l_P * K.c**2 / K.a0), 1),
        ("08", "m_Pl²/M_c = ε_L^{−1/3}·m_Pl = 8,1·10¹² kg", K.m_Pl**2 / K.schnittmassen()[1] / 1e12, 8.1),
        ("08", "m_Pl²/M_c = ε_L^{−1/3}·m_Pl = 8,1·10¹² kg", (K.hbar * K.c**5 / (K.G**2 * K.a0))**(1 / 3) / 1e12, 8.1),
        ("08", "mit 3Δs = r̃, in derselben Länge ℓ_P^{2/3}·ℓ₀^{1/3}", (lambda M: 3 * K.faltung(M) / K.schwelle_spiegel(M))(K.m_Pl**2 / K.schnittmassen()[1]), 1),
        ("08", "mit 3Δs = r̃, in derselben Länge ℓ_P^{2/3}·ℓ₀^{1/3}", (lambda M: 3 * K.faltung(M) / (K.l_P**(2 / 3) * (K.c**2 / K.a0)**(1 / 3)))(K.m_Pl**2 / K.schnittmassen()[1]), 1),
        ("08", "liegen Teilchen und Atome unter der Linie r̃ (Wasserstoff R/r̃ = 10⁻¹⁶)", round(math.log10(K.bohr_radius() / K.schwelle_spiegel(K.m_p + K.m_e))), -16),
        ("08", "Erde und Galaxie darüber (10²⁷ und 10⁴⁹)", round(math.log10(K.R_ERDE / K.schwelle_spiegel(K.M_earth))), 27),
        ("08", "Erde und Galaxie darüber (10²⁷ und 10⁴⁹)", round(math.log10(15 * K.kpc / K.schwelle_spiegel(6e10 * K.M_sun))), 49),
        ("08", "kreuzt r̃ bei 2,0·10⁻⁸ kg und r* bei 3,1·10⁻⁸ kg", K.kreuzung_gehalten(K.bohr_radius())[0] / 1e-8, 2.0),
        ("08", "kreuzt r̃ bei 2,0·10⁻⁸ kg und r* bei 3,1·10⁻⁸ kg", K.kreuzung_gehalten(K.bohr_radius())[1] / 1e-8, 3.1),
        ("08", "Leichtere Körner liegen wie die Atome unter r̃ und über r*", float((lambda M: K.bohr_radius() * (M / K.m_p)**(1 / 3) / K.schwelle_spiegel(M))(1e-10) < 1 < (lambda M: K.bohr_radius() * (M / K.m_p)**(1 / 3) / K.schwellenradius(M))(1e-10)), 1),
        ("08", "schwerere wie Erde und Stern über r̃ und unter r*", float(K.R_ERDE / K.schwellenradius(K.M_earth) < 1 < K.R_ERDE / K.schwelle_spiegel(K.M_earth) and K.R_SONNE / K.schwellenradius(K.M_sun) < 1), 1),
        ("08", "ihre Summe ist die Tiefe, ihr Produkt die Verwindung", (lambda h: (h[0] + h[1]).real / (6 * K.faltung(K.M_sun)))(K.horizonte(0.5 * K.G * K.M_sun**2 / K.c, K.M_sun)), 1),
        ("08", "ihre Summe ist die Tiefe, ihr Produkt die Verwindung", (lambda h: (h[0] * h[1]).real / K.verwindung(K.J_ERDE, K.M_earth)**2)(K.horizonte(K.J_ERDE, K.M_earth)), 1),
        ("08", "und liegen auf dem Kreis |r±| = a_J um den Knoten, im Winkel cos φ = 1/j zur Tiefe", abs(K.horizonte(K.J_ERDE, K.M_earth)[0]) / K.verwindung(K.J_ERDE, K.M_earth), 1),
        ("08", "und liegen auf dem Kreis |r±| = a_J um den Knoten, im Winkel cos φ = 1/j zur Tiefe", (lambda h: h[0].real / abs(h[0]) * h[2])(K.horizonte(K.J_ERDE, K.M_earth)), 1),
        ("08", "(Erde φ = 89,92°)", math.degrees(math.acos(1 / K.horizonte(K.J_ERDE, K.M_earth)[2])), 89.92),
        ("09", "die Spinzahl ist diese Zahl durch die Kopplung, j = (J/ħ)/α_G", K.J_ERDE / K.hbar / K.kopplungen(K.M_earth)[0] / K.horizonte(K.J_ERDE, K.M_earth)[2], 1),
        ("09", "Die Drehung erreicht c, wenn die Zahl ihrer Quanten α_G erreicht, J = GM²/c", K.horizonte(K.G * K.M_sun**2 / K.c, K.M_sun)[2], 1),
        ("01", "denn bei λ̄/α_G schließt die Quantenlänge nach einem Umlauf um die eigene Faltung", (lambda m: m * K.kopplungen(m)[0] * K.c * K.compton(m) / K.kopplungen(m)[0] / K.hbar)(K.m_p), 1),
        ("01", "denn bei λ̄/α_G schließt die Quantenlänge nach einem Umlauf um die eigene Faltung", (lambda m: (K.kopplungen(m)[0] * K.c)**2 * K.compton(m) / K.kopplungen(m)[0] / (K.G * m))(K.m_p), 1),
        ("01", "(beim Proton 1,15 Mpc)", K.compton(K.m_p) / K.kopplungen(K.m_p)[0] / K.Mpc, 1.15),
        ("08", "Zugleich ist r̃² = ħ·τ_D/M", K.hbar * K.rueckkehrzeit() / K.m_p / K.schwelle_spiegel(K.m_p)**2, 1),
        ("08", "r̃ = √(2D·τ_D), wie seine Faltung mit D = GM/c bis r* = √(D·τ_D)", math.sqrt(K.hbar / (2 * K.m_e) * K.rueckkehrzeit()) / (K.schwelle_spiegel(K.m_e) / math.sqrt(2)), 1),
        ("07", "Eine Ausdehnung mit der Rate H gibt H·τ_D = 2π·T_H/T_I = 2π/Θ", K.H0 * K.rueckkehrzeit() / (2 * math.pi * (K.hbar * K.H0 / (2 * math.pi)) / (K.hbar * K.a0 / K.c)), 1),
        ("07", "Eine Ausdehnung mit der Rate H gibt H·τ_D = 2π·T_H/T_I = 2π/Θ", K.H0 * K.rueckkehrzeit() / (2 * math.pi / K.theta_(K.a0, K.H0)), 1),
        ("07", "mit der Temperatur des Horizonts kT_H = ħH/(2π)", K.hbar * K.c / (12 * math.pi * K.c / (6 * K.H0)) / (K.hbar * K.H0 / (2 * math.pi)), 1),
        ("08", "Beim Elektron ist λ̄_e = 510·r*", K.compton(K.m_e) / K.schwellenradius(K.m_e), 510),
        ("08", "sein Feld liegt außerhalb der Quantenlänge ganz unter a₀", float(K.G * K.m_e / K.compton(K.m_e)**2 < K.a0), 1),
        ("08", "beim Proton r* = 154·λ̄_p", K.schwellenradius(K.m_p) / K.compton(K.m_p), 154),
        ("08", "Unter m_Pl²/M_c reicht die eigene Quantenbewegung weiter als die Tiefe, r̃ > 3Δs", float(K.schwelle_spiegel(1e12) > 3 * K.faltung(1e12) and K.schwelle_spiegel(1e14) < 3 * K.faltung(1e14)), 1),
        ("08", "Ihr eigenes Spiegelbild, m̃ = M, ist die Falte bei m_Pl/√2, wo r₀ = λ̄ ist", K.spiegelmasse(K.m_Pl / math.sqrt(2)) / (K.m_Pl / math.sqrt(2)), 1),
        ("08", "Ihr eigenes Spiegelbild, m̃ = M, ist die Falte bei m_Pl/√2, wo r₀ = λ̄ ist", K.compton(K.m_Pl / math.sqrt(2)) / K.grenzradius(K.m_Pl / math.sqrt(2)), 1),
        ("08", "bei der ein Teilchen mit Spin ½ die Spinzahl j = 1 erreicht", K.spinzahl_teilchen(K.m_Pl / math.sqrt(2)), 1),
        ("08", "seine Baryonen tragen (3/4π)·M_U, mit r̃ = 2,05·ℓ_P", K.schwelle_spiegel(3 / (4 * math.pi) * K.schnittmassen()[3]) / K.l_P, 2.05),
        ("08", "Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1", K.leiter_laenge(-1, 0, K.M_sun) / K.compton(K.M_sun), 1),
        ("08", "Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1", K.leiter_laenge(1, 0, K.M_sun) / (3 * K.faltung(K.M_sun)), 1),
        ("08", "Jede Länge der Leiter ist λ̄^p·(3Δs)^q·ℓ₀^r mit p + q + r = 1", K.leiter_laenge(0, 1, K.M_sun) / (K.c**2 / K.a0), 1),
        ("08", "r* (½, ½) und r̃ (−½, ½)", K.leiter_laenge(0.5, 0.5, K.m_e) / K.schwellenradius(K.m_e), 1),
        ("08", "r* (½, ½) und r̃ (−½, ½)", K.leiter_laenge(-0.5, 0.5, K.m_e) / K.schwelle_spiegel(K.m_e), 1),
        ("08", "der Schwerpunkt ist ℓ_P^{2/3}·ℓ₀^{1/3} (0, ⅓)", K.leiter_laenge(0, 1 / 3, K.M_sun) / (K.compton(K.M_sun) * 3 * K.faltung(K.M_sun) * K.c**2 / K.a0)**(1 / 3), 1),
        ("08", "darum teilt M_c die Leiter im Drittel", K.leiter_treffen((-1, 0), (0, 1 / 3))[0] / K.schnittmassen()[1], 1),
        ("08", "darum teilt M_c die Leiter im Drittel", K.leiter_treffen((-1, 0), (0.5, 0.5))[0] / K.schnittmassen()[1], 1),
        ("08", "Die Treffmassen sind x = ±1, wo λ̄ und 3Δs die Ecke ℓ₀ erreichen", K.leiter_treffen((1, 0), (0, 1))[0] / K.schnittmassen()[3], 1),
        ("08", "Die Treffmassen sind x = ±1, wo λ̄ und 3Δs die Ecke ℓ₀ erreichen", K.leiter_treffen((-1, 0), (0, 1))[0] / K.schnittmassen()[0], 1),
        ("08", "λ̄·3Δs = ℓ_P² und ebenso r*·r̃ = ℓ_P·ℓ₀", K.schwellenradius(1e40) * K.schwelle_spiegel(1e40) / (K.l_P * K.c**2 / K.a0), 1),
        ("08", "λ̄·3Δs = ℓ_P² und ebenso r*·r̃ = ℓ_P·ℓ₀", K.schwellenradius(K.m_e) * K.schwelle_spiegel(K.m_e) / (K.l_P * K.c**2 / K.a0), 1),
        ("08", "um λ̄ die Tiefe auf λ̄/α_G (−3, 0)", K.leiter_laenge(-3, 0, K.m_p) / (K.compton(K.m_p) / K.kopplungen(K.m_p)[0]), 1),
        ("08", "Tully-Fisher (1, −1) und α_I (−1, −1)", (K.M_sun / K.m_Pl) * K.leiterzahl() / K.kopplungen(K.M_sun)[2], 1),
        ("08", "Tully-Fisher (1, −1) und α_I (−1, −1)", (K.M_sun / K.m_Pl)**-1 * K.leiterzahl() / K.kopplungen(K.M_sun)[1], 1),
        ("08", "r_e, λ̄_e und a_B sind n = 1, 0 und −1", K.bohr_radius() * K.alpha / K.compton(K.m_e), 1),
        ("08", "α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)²", (K.schwellenradius(K.M_sun) / K.schwelle_spiegel(K.M_sun))**2 / K.kopplungen(K.M_sun)[0], 1),
        ("08", "α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)²", (K.l_P / K.schwellenradius(K.M_sun))**2 / K.kopplungen(K.M_sun)[1], 1),
        ("08", "α_G = (r*/r̃)², α_I = (ℓ_P/r*)² und Tully-Fisher = (ℓ_P/r̃)²", (K.l_P / K.schwelle_spiegel(K.M_sun))**2 / K.kopplungen(K.M_sun)[2], 1),
        ("08", "Zwei Spiegel nacheinander strecken jede Länge um das Quadrat des Abstands ihrer Mitten", (lambda a, b, l: b * b / (a * a / l) / l)(K.l_P, K.schwellenradius(K.M_sun), 1.234e5) * K.kopplungen(K.M_sun)[1], 1),
        ("01", "Am Ort gelesen, gegen den Abstand r, wird es zum Nome: s = r*/r", K.schwellenradius(K.M_sun) / (1e4 * K.AU) / K.s_(K.G * K.M_sun / (1e4 * K.AU)**2), 1),
        ("08", "Das Elektron liegt bei x = −0,36, links von M_c, das Proton bei −0,31", K.leiter_ort(K.m_e), -0.36),
        ("08", "Das Elektron liegt bei x = −0,36, links von M_c, das Proton bei −0,31", K.leiter_ort(K.m_p), -0.31),
        ("08", "Erde und Sonne bei 0,53 und 0,62", K.leiter_ort(K.M_earth), 0.53),
        ("08", "Erde und Sonne bei 0,53 und 0,62", K.leiter_ort(K.M_sun), 0.62),
        ("08", "die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99", K.leiter_ort(6e10 * K.M_sun), 0.79),
        ("08", "die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99", K.leiter_ort(1e15 * K.M_sun), 0.86),
        ("08", "die Milchstraße bei 0,79, ein Haufen von 10¹⁵ M☉ bei 0,86 und das Ganze bei 0,99", K.leiter_ort(3 / (4 * math.pi) * K.schnittmassen()[3]), 0.99),
        ("08", "r* und r̃ sind seine Diagonalen und kreuzen sich in der Mitte", math.log(K.schwellenradius(K.schnittmassen()[0]) / K.l_P) + math.log(K.schwelle_spiegel(K.schnittmassen()[3]) / K.l_P) + 1, 1),
        ("01", "Um den Ort gespiegelt wird die Tiefe zur Länge der Neigung, g = c²/(r²/3Δs)", (1e4 * K.AU)**2 / (3 * K.faltung(K.M_sun)) / (K.c**2 / (K.G * K.M_sun / (1e4 * K.AU)**2)), 1),
        ("01", "die Schwelle ist der Ort, an dem dieses Spiegelbild die Schließung ℓ₀ erreicht", K.schwellenradius(K.M_sun)**2 / (3 * K.faltung(K.M_sun)) / (K.c**2 / K.a0), 1),
        ("08", "Aus M, ħ, G, c und a₀ lassen sich nur zwei reine Zahlen bilden, M/m_Pl und ε_L", 5 - __import__("numpy").linalg.matrix_rank(__import__("numpy").array([[1, 1, -1, 0, 0], [0, 2, 3, 1, 1], [0, -1, -2, -1, -2]])), 2),
        ("08", "der Horizont 2·3Δs, die Wärmelänge der Falte 8π·3Δs und c mal ihrer Lebensdauer 5120π·α_G·3Δs", K.hbar * K.c / (K.k_B * K.falte_temperatur(K.M_sun)) / K.leiter_laenge(1, 0, K.M_sun) / (8 * math.pi), 1),
        ("08", "der Horizont 2·3Δs, die Wärmelänge der Falte 8π·3Δs und c mal ihrer Lebensdauer 5120π·α_G·3Δs", K.c * K.lebensdauer_falte(K.M_sun) / K.leiter_laenge(3, 0, K.M_sun) / (5120 * math.pi), 1),
        ("01", "Kreuzung in der Mitte nach dem ganzen, ω·t = π/(2√2) nach der Umkehr", __import__("scipy.integrate", fromlist=["quad"]).quad(lambda r: 1 / math.sqrt(2 * (1 / r - 1)), 0, 1)[0] / (math.pi / (2 * math.sqrt(2))), 1),
        ("01", "der gebundene Bereich als rollender Kreis, r = (r_U/2)·(1 − cos θ)", (lambda th: 0.5 * (1 - math.cos(th)) - (lambda x: x)(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda u: 1 / math.sqrt(2 * (1 / u - 1)), r, 1)[0] - (math.pi - (th - math.sin(th))) / (2 * math.sqrt(2)), 1e-9, 1 - 1e-12)) + 1)(2.0), 1),
        ("03", "Mit g = (1 + Q)·g_bar ist am Rand R_Δ jedes Haufens a_N = Δ·H²·R/(2(1 + Q))", (lambda M, Q: (lambda R: K.G * M / (1 + Q) / R**2 / (200 * K.H0**2 * R / (2 * (1 + Q))))(K.uebergabe_linse(M, K.H0, Q)[0]))(7e14 * K.M_sun, K.q_je_radiant()), 1),
        ("03", "R = 2(1 + Q)·(ln q)²·a₀/(Δ·H²), mit 1 + Q = 2π also R = 2Θ·(ln q)²·(c/H)/Δ, für Δ = 200 heute 1,3 Mpc", (lambda H: K.naht_haufen(200, H, K.q_je_radiant())[0] / K.uebergabe_linse(_uebergabe()[3] * K.M_sun, H, K.q_je_radiant())[0])(_cmb()["saett"]["punkt"]["H0"] * 1e3 / K.Mpc), 1),
        ("03", "R = 2(1 + Q)·(ln q)²·a₀/(Δ·H²), mit 1 + Q = 2π also R = 2Θ·(ln q)²·(c/H)/Δ, für Δ = 200 heute 1,3 Mpc", (lambda H: 2 * K.theta_(K.a0, H) * math.log(1 - 1 / (2 * math.pi))**2 * K.c / H / 200 / K.Mpc)(_cmb()["saett"]["punkt"]["H0"] * 1e3 / K.Mpc), 1.3),
        ("03", "Dort liegt die Masse M = 8·(ln q)⁶·Θ³·M_H/Δ² = 2,5·10¹⁴ M☉", (lambda H: 8 * math.log(1 - 1 / (2 * math.pi))**6 * K.theta_(K.a0, H)**3 * K.c**3 / (2 * K.G * H) / 200**2 / K.M_sun / 1e14)(_cmb()["saett"]["punkt"]["H0"] * 1e3 / K.Mpc), 2.5),
        ("03", "Dort liegt die Masse M = 8·(ln q)⁶·Θ³·M_H/Δ² = 2,5·10¹⁴ M☉", _uebergabe()[3] / 1e14, 2.5),
        ("01", "ohne Rückstellung linear δ_c = (3/20)·(12π)^{2/3} = 1,686", K.kollaps_kugel()[0], 1.686),
        ("01", "ohne Rückstellung linear δ_c = (3/20)·(12π)^{2/3} = 1,686", 3 / 20 * (12 * math.pi)**(2 / 3), 1.686),
        ("01", "und auf das 18π²-fache des Mittels verdichtet", K.kollaps_kugel()[1] / (18 * math.pi**2), 1),
        ("01", "Druck des Lichts p = u/3", __import__("scipy.integrate", fromlist=["quad"]).quad(lambda m: m * m, -1, 1)[0] / 2 * 3, 1),
        ("03", "R_max = v_flach/H_∞, mit H_∞ = c·√Λ_s = 55,6 km/s/Mpc", K.groesste_umkehr(K.M_sun)[1] * K.Mpc / 1e3, 55.6),
        ("03", "R_max = v_flach/H_∞, mit H_∞ = c·√Λ_s = 55,6 km/s/Mpc", K.H0 * math.sqrt(3 * K.vakuumteil(K.a0, K.H0)) * K.Mpc / 1e3, 55.6),
        ("03", "für die Milchstraße (6·10¹⁰ M☉) sind es 3,1 Mpc", K.groesste_umkehr(6e10 * K.M_sun)[0] / K.Mpc, 3.1),
        ("03", "für die Milchstraße (6·10¹⁰ M☉) sind es 3,1 Mpc", (lambda M, H: __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: K.bahntempo_regel(M, r)**2 / r - H**2 * r, 1e3 * K.kpc, 100 * K.Mpc))(6e10 * K.M_sun, K.groesste_umkehr(K.M_sun)[1]) / K.Mpc, 3.1),
        ("03", "R_max³ = (1 + Q)·GM_b/H_∞², bei M₂₀₀ = 10¹⁵ M☉ 11 Mpc", (lambda M, H: (K.G * M / H**2)**(1 / 3))(1e15 * K.M_sun, K.groesste_umkehr(K.M_sun)[1]) / K.Mpc, 11),
        ("03", "R_max³ = (1 + Q)·GM_b/H_∞², bei M₂₀₀ = 10¹⁵ M☉ 11 Mpc", (lambda M, Q, H: __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda x: (1 + Q) * K.G * M / (1 + Q) / (x * K.Mpc) + H**2 * (x * K.Mpc)**2 / 2, bounds=(0.1, 100), method="bounded", options={"xatol": 1e-9}).x)(1e15 * K.M_sun, K.q_je_radiant(), K.groesste_umkehr(K.M_sun)[1]), 11),
        ("03", "denn −d(ln n_B)/ds = ν und damit dΦ/dr = ν·a_N", (lambda s, h: -(math.log(1 / math.expm1(s + h)) - math.log(1 / math.expm1(s - h))) / (2 * h) / K.verstaerkung(s**2 * K.a0))(0.7, 1e-5), 1),
        ("03", "denn −d(ln n_B)/ds = ν und damit dΦ/dr = ν·a_N", (lambda M, r: (lambda f: (f(r * 1.0001) - f(r * 0.9999)) / (0.0002 * r) / (K.verstaerkung(K.G * M / r**2) * K.G * M / r**2))(lambda x: K.v_flach(M)**2 * math.log(1 / math.expm1(K.schwellenradius(M) / x))))(6e10 * K.M_sun, 20 * K.kpc), 1),
        ("03", "Newton allein ist v_flach²·ln q, der Logarithmus des Rückkehranteils", (lambda M, r: -K.G * M / r / (K.v_flach(M)**2 * math.log(K.rueckkehranteil(K.G * M / r**2))))(6e10 * K.M_sun, 20 * K.kpc), 1),
        ("03", "Die Rückkehr fügt den Logarithmus der Summe der Rückläufe hinzu, q·ν = n_B", (lambda M, r: (-K.G * M / r + K.v_flach(M)**2 * math.log(K.verstaerkung(K.G * M / r**2))) / (K.v_flach(M)**2 * math.log(1 / math.expm1(K.schwellenradius(M) / r))))(6e10 * K.M_sun, 20 * K.kpc), 1),
        ("03", "Die Bindung zwischen zwei Orten ist so v_flach² mal dem Logarithmus des Verhältnisses ihrer Rückkehrzahlen", (lambda M, r1, r2: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.verstaerkung(K.G * M / x**2) * K.G * M / x**2, r1, r2, limit=500)[0] / (K.v_flach(M)**2 * math.log(math.expm1(K.schwellenradius(M) / r1) / math.expm1(K.schwellenradius(M) / r2))))(6e10 * K.M_sun, 5 * K.kpc, 500 * K.kpc), 1),
        ("03", "mit n_B(r*/R_max) ≈ R_max/r* = (√15/2π²)·c/v_flach", (lambda M: (1 / math.expm1(K.schwellenradius(M) / K.groesste_umkehr(M)[0])) / (K.groesste_umkehr(M)[0] / K.schwellenradius(M)))(6e10 * K.M_sun), 1),
        ("03", "an der Schwelle 3,43·v_flach = 585 km/s", (lambda M: K.flucht_regel(M, K.schwellenradius(M)) / K.v_flach(M))(6e10 * K.M_sun), 3.43),
        ("03", "an der Schwelle 3,43·v_flach = 585 km/s", (lambda M: math.sqrt(2 * math.log(math.expm1(1.0) * (math.sqrt(15) / (2 * math.pi**2) * K.c / K.v_flach(M))) - 1))(6e10 * K.M_sun), 3.43),
        ("03", "an der Schwelle 3,43·v_flach = 585 km/s", (lambda M: (lambda rs, H: (lambda Re: math.sqrt(2 * __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.verstaerkung(K.G * M / x**2) * K.G * M / x**2 - H**2 * x, rs, Re, limit=500)[0]) / 1e3)(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: K.verstaerkung(K.G * M / x**2) * K.G * M / x**2 - H**2 * x, 10 * K.kpc, 100 * K.Mpc)))(K.schwellenradius(M), K.groesste_umkehr(M)[1]))(6e10 * K.M_sun), 585),
        ("03", "am Ort der Sonne (8,25 kpc) 591 km/s", K.flucht_regel(6e10 * K.M_sun, K.SONNE_GALAKTISCH["R0"] * K.kpc) / 1e3, 591),
        ("03", "am Ort der Sonne (8,25 kpc) 591 km/s", K.SONNE_GALAKTISCH["R0"], 8.25),
        ("03", "√2·v_flach = 241 km/s", math.sqrt(2 * K.G * 6e10 * K.M_sun / K.schwellenradius(6e10 * K.M_sun)) / 1e3, 241),
        ("03", "√2·v_flach = 241 km/s", math.sqrt(2) * K.v_flach(6e10 * K.M_sun) / 1e3, 241),
        ("03", "Sie teilt 4Ω² zwischen Gezeit mit Umlauf und Epizyklen", (lambda g, D, h: (lambda Om2: (-(g(D * (1 + h)) - g(D * (1 - h))) / (2 * h * D) + Om2 + 4 * Om2 + D * (g(D * (1 + h)) / (D * (1 + h)) - g(D * (1 - h)) / (D * (1 - h))) / (2 * h * D)) / Om2)(g(D) / D))((lambda R: K.verstaerkung(K.G * 6e10 * K.M_sun / R**2) * K.G * 6e10 * K.M_sun / R**2), 50 * K.kpc, 1e-4), 4),
        ("03", "Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc", K.gezeitenrand(1e7 * K.M_sun, 6e10 * K.M_sun, 50 * K.kpc) / K.kpc, 3.8),
        ("03", "Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc", (lambda g, ge, D, Om2: (__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: g(D - r) - ge(r) - Om2 * (D - r), 1e-2 * K.kpc, 20 * K.kpc) + __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda r: g(D + r) + ge(r) - Om2 * (D + r), 1e-2 * K.kpc, 20 * K.kpc)) / 2 / K.kpc)((lambda R: K.verstaerkung(K.G * 6e10 * K.M_sun / R**2) * K.G * 6e10 * K.M_sun / R**2), lambda r: K.verstaerkung(K.G * 1e7 * K.M_sun / r**2) * K.G * 1e7 * K.M_sun / r**2, 50 * K.kpc, (lambda R: K.verstaerkung(K.G * 6e10 * K.M_sun / R**2) * K.G * 6e10 * K.M_sun / R**2)(50 * K.kpc) / (50 * K.kpc)), 3.8),
        ("03", "Ein Bereich von 10⁷ M☉ im Abstand von 50 kpc um die Milchstraße (Punktmasse) bindet bis 3,8 kpc, nach Newton mit denselben Baryonen bis 1,9 kpc", 50 * K.kpc * (1e7 * K.M_sun / (3 * 6e10 * K.M_sun))**(1 / 3) / K.kpc, 1.9),
        ("03", "Verliert er neun Zehntel davon als Gas, rückt der Rand auf 3,7 kpc, ohne dieses Gedächtnis auf 2,1 kpc", K.gezeitenrand(1e6 * K.M_sun, 6e10 * K.M_sun, 50 * K.kpc, 1e7 * K.M_sun) / K.kpc, 3.7),
        ("03", "Verliert er neun Zehntel davon als Gas, rückt der Rand auf 3,7 kpc, ohne dieses Gedächtnis auf 2,1 kpc", K.gezeitenrand(1e6 * K.M_sun, 6e10 * K.M_sun, 50 * K.kpc) / K.kpc, 2.1),
        ("03", "also r_J = v_flach/(√d_eff·Ω), linear im Tempo", (lambda aN: K.v_flach(1e7 * K.M_sun) / (math.sqrt(3 - K.schalter(aN)) * K.bahntempo_regel(6e10 * K.M_sun, 50 * K.kpc) / (50 * K.kpc)) / K.gezeitenrand(1e7 * K.M_sun, 6e10 * K.M_sun, 50 * K.kpc))(K.G * 6e10 * K.M_sun / (50 * K.kpc)**2), 1),
        ("03", "gibt die Jeans-Gleichung σ² = v_flach²/γ an jedem Radius", (lambda M: K.bahntempo_regel(M, 1e3 * K.schwellenradius(M))**2 / K.v_flach(M)**2)(1e10 * K.M_sun), 1),
        ("03", "hat nur mit σ < v_flach/√3 = 0,58·v_flach eine endliche Masse", (lambda M: K.bahntempo_regel(M, 1e4 * K.schwellenradius(M)) / (math.sqrt(3) * K.v_flach(M)))(1e10 * K.M_sun), 0.58),
        ("03", "Trägt sich eine isotherme Kugel ganz im Tiefbereich selbst, ist σ⁴ = v_flach⁴/γ² = (4/81)·GMa₀, also γ = 9/2", (lambda M: (K.v_flach(M) / K.streuung_kugel(M))**2)(1e10 * K.M_sun), 4.5),
        ("03", "Sie fällt außen wie r^{−9/2}, Virialsatz und Profil lesen dieselbe Zahl", K.isotherme_kugel()[0], 4.5),
        ("03", "Sie fällt außen wie r^{−9/2}, Virialsatz und Profil lesen dieselbe Zahl", K.isotherme_kugel()[1], 1),
        ("03", "Darum ist 9/2 die Obergrenze, mit Newtonschem Kern fällt die Kugel flacher", float(K.isotherme_kugel(1.0)[0] < 4.5 and K.isotherme_kugel(1e-6)[0] < 4.5), 1),
        ("03", "Sein Spiegelbild um r* ist c_s²/a₀ = ℓ₀·(c_s/c)²", (lambda M, cs: K.G * M / cs**2 * (cs**2 / K.a0) / K.schwellenradius(M)**2)(1e10 * K.M_sun, 10e3), 1),
        ("03", "Das Bahntempo der Regel fällt nie unter v_flach", float(min(K.bahntempo_regel(1e10 * K.M_sun, x * K.schwellenradius(1e10 * K.M_sun)) for x in (0.1, 1, 10, 100, 1e3, 1e4)) > K.v_flach(1e10 * K.M_sun)), 1),
        ("03", "f(S) = δ_c/√(2π)·S^{−3/2}·e^{−δ_c²/(2S)} hat die Form T^{−3/2}·e^{−τ_D/(4T)}, mit S/(2δ_c²) statt T/τ_D", (lambda dc, S: K.rueckkehr_dichte(K.rueckkehrzeit() * S / (2 * dc**2)) * K.rueckkehrzeit() / (2 * dc**2) / (dc / math.sqrt(2 * math.pi) * S**-1.5 * math.exp(-dc**2 / (2 * S))))(K.kollaps_kugel()[0], 1.7), 1),
        ("03", "Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c, der typisch bindenden Masse M_z, und dort sind es 0,317", K.rueckkehr_angekommen(K.rueckkehrzeit() / 2), 0.317),
        ("03", "Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c, der typisch bindenden Masse M_z, und dort sind es 0,317", __import__("scipy.integrate", fromlist=["quad"]).quad((lambda S: K.kollaps_kugel()[0] / math.sqrt(2 * math.pi) * S**-1.5 * math.exp(-K.kollaps_kugel()[0]**2 / (2 * S))), 0, K.kollaps_kugel()[0]**2)[0], 0.317),
        ("03", "Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c", (lambda dc: math.exp(-__import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda lM: -(dc * math.exp(-lM / 6))**2 * dc / math.sqrt(2 * math.pi) * (dc * math.exp(-lM / 6))**-3 * math.exp(-dc**2 / (2 * (dc * math.exp(-lM / 6))**2)), bounds=(-10, 10), method="bounded", options={"xatol": 1e-10}).x / 6))(K.kollaps_kugel()[0]), 1),
        ("03", "Je logarithmischem Schritt der Masse gezählt liegt der Gipfel für jede Neigung des Spektrums bei σ = δ_c", (lambda dc: math.exp(-__import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda lM: -(dc * math.exp(-lM / 2))**2 * dc / math.sqrt(2 * math.pi) * (dc * math.exp(-lM / 2))**-3 * math.exp(-dc**2 / (2 * (dc * math.exp(-lM / 2))**2)), bounds=(-10, 10), method="bounded", options={"xatol": 1e-10}).x / 2))(K.kollaps_kugel()[0]), 1),
        ("03", "Ebenso hat die Rückkehr je logarithmischem Schritt der Zeit ihren Gipfel bei T = τ_D/2", 1 / __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda T: -T * K.rueckkehr_dichte(T * K.rueckkehrzeit()), bounds=(0.01, 5), method="bounded", options={"xatol": 1e-10}).x, 2),
        ("03", "Der Faktor 2 von Press und Schechter ist der Spiegel an der Schwelle", __import__("scipy.integrate", fromlist=["quad"]).quad((lambda S: K.kollaps_kugel()[0] / math.sqrt(2 * math.pi) * S**-1.5 * math.exp(-K.kollaps_kugel()[0]**2 / (2 * S))), 0, K.kollaps_kugel()[0]**2)[0] / (0.5 * math.erfc(1 / math.sqrt(2))), 2),
        ("03", "Ohne ihn fände über alle Skalen nur die Hälfte der Baryonen einen Bereich", __import__("scipy.integrate", fromlist=["quad"]).quad((lambda S: K.kollaps_kugel()[0] / math.sqrt(2 * math.pi) * S**-1.5 * math.exp(-K.kollaps_kugel()[0]**2 / (2 * S))), 0, math.inf)[0] / 2, 0.5),
        ("F:lokale_gruppe", "liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database)", D.m31_paar()[0], 784),
        ("F:lokale_gruppe", "liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database)", D.m31_paar()[1], 90),
        ("F:lokale_gruppe", "liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database)", _bereiche_paar()[0], 784),
        ("F:lokale_gruppe", "liegen 784 kpc auseinander und nähern sich mit 90 km/s (Local Volume Database)", _bereiche_paar()[1], 90),
        ("F:lokale_gruppe", "Newton verlangt für 90 km/s das 19-fache der Baryonen", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda f: K.zeitargument_newton(f * sum(_bereiche()), D.m31_paar()[0] * K.kpc, K.weltalter(lambda a: 3 * K.dichte_vakuumteil()), K.groesste_umkehr(K.M_sun)[1]) - D.m31_paar()[1] * 1e3, 10, 60), 19),
        ("F:lokale_gruppe", "aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉)", D.bereich_baryonen("mw")[0] / 1e10, 6.6),
        ("F:lokale_gruppe", "aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉)", (K.M_B_MILCHSTRASSE + D.satellit("lmc")["M_b"] + D.satellit("smc")["M_b"]) / 1e10, 6.6),
        ("F:lokale_gruppe", "aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉)", D.bereich_baryonen("m_031")[0] / 1e11, 1.49),
        ("F:lokale_gruppe", "aus den Baryonen der Bereiche (6,6·10¹⁰ und 1,49·10¹¹ M☉)", (K.M_B_M31 + 10**9.82 + 10**8.86 + 10**8.84 + 1.33 * 10**5.61 + 10**8.77 + 10**8.47 + 1.33 * 10**5.09 + 10**8.21 + 1.33 * 10**7.67) / 1e11, 1.49),
        ("F:lokale_gruppe", "v₂ = 214 km/s aus den Baryonen der Bereiche", K.bahntempo_regel(_bereiche()[1], 1e4 * K.Mpc) / 1e3, 214),
        ("F:lokale_gruppe", "v₂ = 214 km/s aus den Baryonen der Bereiche", K.v_flach(_bereiche()[1]) / 1e3, 214),
        ("F:lokale_gruppe", "ist das Paar radial seit 13,1 Gyr auf der Bahn", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, gegenzug=True)[2] / K.Gyr, 13.1),
        ("F:lokale_gruppe", "ist das Paar radial seit 13,1 Gyr auf der Bahn", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, gegenzug=True, numerisch=True)[2] / K.Gyr, 13.1),
        ("F:lokale_gruppe", "mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True)[2] / K.Gyr, 13.5),
        ("F:lokale_gruppe", "mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr", K.zeitargument(*_bereiche(), D.m31_paar()[0] * K.kpc, D.m31_paar()[1] * 1e3, v_t=D.m31_paar(True)[2] * 1e3, gegenzug=True, numerisch=True)[2] / K.Gyr, 13.5),
        ("F:lokale_gruppe", "mit dieser Bewegung ist das Paar seit 13,5 Gyr auf der Bahn, gegen das Weltalter von 13,8 Gyr", K.weltalter(lambda a: K.dichte_de_planck()) / K.Gyr, 13.8),
        ("03", "also H = H_∞·coth(3H_∞t/2) und Θ = Θ_∞·tanh(3H_∞t/2) mit Θ_∞ = √15/π", (lambda Hi, rm: (lambda a: (lambda t: Hi / math.tanh(1.5 * Hi * t) / math.sqrt(8 * math.pi * K.G / 3 * (rm / a**3 + 3 * K.dichte_vakuumteil())))(__import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: 1 / (x * math.sqrt(8 * math.pi * K.G / 3 * (rm / x**3 + 3 * K.dichte_vakuumteil()))), 1e-12, a, limit=400)[0]))(0.5))(K.groesste_umkehr(K.M_sun)[1], K.dichte_heute()[2]), 1),
        ("03", "Mit H_∞·τ_D = 2π/Θ_∞ ist die Zeit", K.groesste_umkehr(K.M_sun)[1] * K.rueckkehrzeit() / (2 * math.pi / (math.sqrt(15) / math.pi)), 1),
        ("03", "bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95", K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc), 1.020),
        ("03", "bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95", K.uhr_theta(K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc))[0] / K.Gyr, 13.8),
        ("03", "bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95", K.weltalter(lambda a: K.dichte_dunkle_energie()) / K.Gyr, 13.8),
        ("03", "bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95", K.uhr_theta(K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc))[1], 0.95),
        ("03", "bei Θ = 1,020 heute 13,8 Mrd. Jahre, H₀·t₀ = 0,95", K.weltalter(lambda a: K.dichte_dunkle_energie()) * K.atmen_heute()[0] * 1e3 / K.Mpc, 0.95),
        ("03", "im Formteil ½ − (3/2)·(Θ/Θ_∞)², heute −0,53", K.uhr_theta(K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc))[2], -0.53),
        ("03", "im Formteil ½ − (3/2)·(Θ/Θ_∞)², heute −0,53", (lambda rr, rb, rm: (rr + rm / 2 - K.dichte_dunkle_energie()) / (rr + rm + K.dichte_dunkle_energie()))(*K.dichte_heute()), -0.53),
        ("03", "bei Θ = Θ_∞/√3 = 0,712", math.sqrt(15) / math.pi / math.sqrt(3), 0.712),
        ("03", "bei Θ = Θ_∞/√3 = 0,712", (lambda rr, rb, rm: (lambda a: 2 * math.pi * K.a0 / (K.c * math.sqrt(8 * math.pi * K.G / 3 * (rr / a**4 + rm / a**3 + K.dichte_dunkle_energie()))))(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda a: rr / a**4 + rm / (2 * a**3) - K.dichte_dunkle_energie(), 0.1, 2.0)))(*K.dichte_heute()), 0.712),
        ("03", "Formteil und Faltung sind gleich bei Θ_∞/√2 = 0,872", math.sqrt(15) / math.pi / math.sqrt(2), 0.872),
        ("03", "Formteil und Faltung sind gleich bei Θ_∞/√2 = 0,872", (lambda rr, rb, rm: (lambda a: 2 * math.pi * K.a0 / (K.c * math.sqrt(8 * math.pi * K.G / 3 * (rr / a**4 + rm / a**3 + K.dichte_dunkle_energie()))))(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda a: rm / a**3 - K.dichte_dunkle_energie(), 0.1, 2.0)))(*K.dichte_heute()), 0.872),
        ("03", "Beide Stufen hängen nur an Θ, nicht an a₀", (lambda a0s: max(abs(K.theta_(a0s, math.sqrt(8 * math.pi * K.G / 3 * (K.dichte_heute()[2] / ab**3 + 3 * math.pi**3 * a0s**2 / (30 * K.G * K.c**2))) ) - math.sqrt(15) / math.pi / math.sqrt(3)) for ab in [__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda a: K.dichte_heute()[2] / (2 * a**3) - 3 * math.pi**3 * a0s**2 / (30 * K.G * K.c**2), 0.05, 5.0)]))(1.2e-10) < 1e-9, 1),
        ("03", "dauert das Fenster 0,68 Mrd. Jahre", (lambda th_inf, a, b: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda u: 1 / (3 * u * 2 * math.pi * K.a0 / (K.c * th_inf) * math.sqrt(1 + u)), 15 / (math.pi**2 * b**2) - 1, 15 / (math.pi**2 * a**2) - 1)[0])(math.sqrt(15) / math.pi, 1.0, K.schleifen_ganzes()[2]) / 3.15576e16, 0.68),
        ("03", "im Formteil heute 0,51·L und 0,18·L", K.lichtweg_zelle()[0], 0.51),
        ("03", "im Formteil heute 0,51·L und 0,18·L", (lambda lo, hi, rr, rm: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.c / (math.exp(x) * math.sqrt(8 * math.pi * K.G / 3 * (rr * math.exp(-4 * x) + rm * math.exp(-3 * x) + K.dichte_dunkle_energie()))), lo, hi, limit=800, epsabs=0, epsrel=1e-11)[0] / K.zelle_laenge())(-40, 0, K.dichte_heute()[0], K.dichte_heute()[2]), 0.51),
        ("03", "im Formteil heute 0,51·L und 0,18·L", K.lichtweg_zelle()[1], 0.18),
        ("03", "im Formteil heute 0,51·L und 0,18·L", __import__("scipy.integrate", fromlist=["quad"]).quad(lambda a: K.c / (a * a * math.sqrt(8 * math.pi * K.G / 3 * (K.dichte_heute()[0] / a**4 + K.dichte_heute()[2] / a**3 + K.dichte_dunkle_energie()))), 1, math.inf, limit=800)[0] / K.zelle_laenge(), 0.18),
        ("03", "c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70", (math.sqrt(15) / math.pi / (8 * math.pi**2))**(1 / 3) * __import__("scipy.integrate", fromlist=["quad"]).quad(lambda u: 1 / math.sqrt(u + u**4), 0, math.inf, limit=400)[0], 0.70),
        ("03", "c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70", (math.sqrt(15) / math.pi / (8 * math.pi**2))**(1 / 3) * math.gamma(1 / 6) * math.gamma(1 / 3) / math.gamma(0.5) / 3, 0.70),
        ("03", "c·η_∞/L = (Θ_∞/8π²)^{1/3}·∫₀^∞ du/√(u + u⁴) = 0,70", __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.c / (math.exp(x) * math.sqrt(8 * math.pi * K.G / 3 * (K.dichte_heute()[2] / math.exp(3 * x) + K.dichte_dunkle_energie(1.2e-10)))), -40, 60, limit=800, epsabs=0, epsrel=1e-11)[0] / K.zelle_laenge(a0_=1.2e-10), 0.70),
        ("03", "Mit der Strahlung ist es 0,69", K.lichtweg_zelle()[2], 0.69),
        ("03", "Mit der Strahlung ist es 0,69", (lambda lo, hi, rr, rm: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.c / (math.exp(x) * math.sqrt(8 * math.pi * K.G / 3 * (rr * math.exp(-4 * x) + rm * math.exp(-3 * x) + K.dichte_dunkle_energie()))), lo, hi, limit=800, epsabs=0, epsrel=1e-11)[0] / K.zelle_laenge())(-40, 60, K.dichte_heute()[0], K.dichte_heute()[2]), 0.69),
        ("03", "Der Radius der letzten Streufläche wächst dabei auf 0,68·L", K.lichtweg_zelle()[3], 0.68),
        ("03", "Der Radius der letzten Streufläche wächst dabei auf 0,68·L", (lambda lo, hi, rr, rm: __import__("scipy.integrate", fromlist=["quad"]).quad(lambda x: K.c / (math.exp(x) * math.sqrt(8 * math.pi * K.G / 3 * (rr * math.exp(-4 * x) + rm * math.exp(-3 * x) + K.dichte_dunkle_energie()))), lo, hi, limit=800, epsabs=0, epsrel=1e-11)[0] / K.zelle_laenge())(-math.log(1 + K.z_entkopplung()), 60, K.dichte_heute()[0], K.dichte_heute()[2]), 0.68),
        ("03", "die drei Kreispaare erscheinen und wachsen bis höchstens α = 43°", K.lichtweg_zelle()[4], 43),
        ("03", "die drei Kreispaare erscheinen und wachsen bis höchstens α = 43°", max(K.lichtweg_zelle(a0_=x)[4] for x in (1.0e-10, 1.1e-10, 1.2e-10)), 43),
        ("03", "Ihr Durchmesser, 1,36·L, bleibt unter dem Abstand √2·L der Abbilder über die Diagonalen der Flächen", 2 * K.lichtweg_zelle()[3], 1.36),
        ("03", "Ihr Durchmesser, 1,36·L, bleibt unter dem Abstand √2·L der Abbilder über die Diagonalen der Flächen", float(2 * max(K.lichtweg_zelle(a0_=x)[3] for x in (1.0e-10, 1.2e-10)) < math.sqrt(2)), 1),
        ("03", "= (5/6)·B(5/6, 2/3)·a_× = 1,44·a_×", K.wachstum_ende(1.0)[0], 1.44),
        ("03", "= (5/6)·B(5/6, 2/3)·a_× = 1,44·a_×", 5 / 6 * math.gamma(5 / 6) * math.gamma(2 / 3) / math.gamma(1.5), 1.44),
        ("03", "Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 %", 100 * K.wachstum_ende(1.0)[1], 61),
        ("03", "Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 %", (lambda s, ax: 100 * s.sol(math.log(ax))[0] / s.y[0][-1])((lambda rm, rde, bis: (lambda s: s)(__import__("scipy.integrate", fromlist=["solve_ivp"]).solve_ivp(lambda x, y: [y[1], -(2 - 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde)) * y[1] + 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde) * y[0]], (math.log(1e-4), bis), [1e-4, 1e-4], rtol=1e-10, atol=1e-14, dense_output=True)))(K.dichte_heute()[2], K.dichte_dunkle_energie(), 25.0), (K.dichte_heute()[2] / K.dichte_dunkle_energie())**(1 / 3)), 61),
        ("03", "Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 %", 100 * K.wachstum_ende(K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc))[2], 71),
        ("03", "Bei a_× sind es 61 %, heute, bei Θ = 1,020, 71 %", (lambda s: 100 * s.sol(0.0)[0] / s.y[0][-1])((lambda rm, rde, bis: (lambda s: s)(__import__("scipy.integrate", fromlist=["solve_ivp"]).solve_ivp(lambda x, y: [y[1], -(2 - 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde)) * y[1] + 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde) * y[0]], (math.log(1e-4), bis), [1e-4, 1e-4], rtol=1e-10, atol=1e-14, dense_output=True)))(K.dichte_heute()[2], K.dichte_dunkle_energie(), 25.0)), 71),
        ("03", "Linear wächst Struktur von heute an im Formteil nur noch um 41 %", 100 / K.wachstum_ende(K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc))[2] - 100, 41),
        ("03", "Linear wächst Struktur von heute an im Formteil nur noch um 41 %", (lambda s: 100 * s.y[0][-1] / s.sol(0.0)[0] - 100)((lambda rm, rde, bis: (lambda s: s)(__import__("scipy.integrate", fromlist=["solve_ivp"]).solve_ivp(lambda x, y: [y[1], -(2 - 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde)) * y[1] + 1.5 * rm * math.exp(-3 * x) / (rm * math.exp(-3 * x) + rde) * y[0]], (math.log(1e-4), bis), [1e-4, 1e-4], rtol=1e-10, atol=1e-14, dense_output=True)))(K.dichte_heute()[2], K.dichte_dunkle_energie(), 25.0)), 41),
        ("03", "denn (a_×/a)³ = (Θ_∞/Θ)² − 1 mit Θ_∞ = √15/π", (K.dichte_heute()[2] / K.dichte_dunkle_energie()) / ((math.sqrt(15) / math.pi / K.theta_(K.a0, math.sqrt(8 * math.pi * K.G / 3 * (K.dichte_heute()[2] + K.dichte_dunkle_energie()))))**2 - 1), 1),
        ("08", "Bei 6Δs steht die ruhende Uhr", 3 / __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda d: K.takte_kreisbahn(d)["ruhend"] - 1e-9, 0.1, 0.5), 6),
        ("08", "Bei 9Δs steht die umlaufende", 3 * __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda r: -(1 - 2 / r) / r**2, bounds=(2.1, 10), method="bounded", options={"xatol": 1e-10}).x, 9),
        ("08", "dort läuft die Bahn für die ruhende Uhr am Ort mit c", K.takte_kreisbahn(1 / 3 - 1e-12)["tempo_ort"], 1),
        ("08", "Bei 18Δs steht die radiale Schwingung", 3 * __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda r: r**2 / (r - 3), bounds=(3.1, 20), method="bounded", options={"xatol": 1e-10}).x, 18),
        ("08", "die Bahn läuft dort für die ruhende Uhr mit c/2", 1 / K.takte_kreisbahn(1 / 6)["tempo_ort"], 2),
        ("08", "am kleinsten bei 18Δs und bei 12Δs = 4d genau eins", 3 / __import__("scipy.optimize", fromlist=["minimize_scalar"]).minimize_scalar(lambda d: K.takte_kreisbahn(d)["energie"], bounds=(0.01, 0.33), method="bounded", options={"xatol": 1e-10}).x, 18),
        ("08", "am kleinsten bei 18Δs und bei 12Δs = 4d genau eins", 3 / __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda d: K.takte_kreisbahn(d)["energie"] - 1, 0.17, 0.33), 12),
        ("08", "und läuft gegen die Fernen je Umlauf um 2π·(1 − √(1 − 3δ)) vor, schwach 3πδ", K.takte_kreisbahn(1e-7)["kreisel"] / 1e-7 / math.pi, 3),
        ("08", "und läuft gegen die Fernen je Umlauf um 2π·(1 − √(1 − 3δ)) vor, schwach 3πδ", (lambda M, r: K.kreisel_bahn(M, r) * 2 * math.pi / math.sqrt(K.G * M / r**3) / (3 * math.pi * 3 * K.faltung(M) / r))(K.M_sun, K.AU), 1),
        ("08", "In erster Ordnung ist Re τ_M = 1 − 3δ: Um 3δ verfehlt die Bahn den vollen Umlauf, die Drehung des Perihels", (lambda M, r: 2 * math.pi * (1 / math.sqrt(K.takte_kreisbahn(3 * K.faltung(M) / r)["radial"]) - 1) / K.perihel_je_umlauf(M, r, 0.0))(K.M_sun, K.AU), 1),
        ("08", "An der Photonensphäre ist (Re τ_M)² = −1", K.takte_kreisbahn(1 / 3)["radial"], -1),
        ("08", "mit b_c = 9√3·Δs, dem Rand des Schattens", K.schatten(K.M_sun, 1.0) / 2 / K.faltung(K.M_sun) / math.sqrt(3), 9),
        ("08", "sein Umlauf c/(2π·b_c) bei einer Sonnenmasse 6,2 kHz", K.c / (2 * math.pi * 9 * math.sqrt(3) * K.faltung(K.M_sun)) / 1e3, 6.2),
        ("08", "sein Umlauf c/(2π·b_c) bei einer Sonnenmasse 6,2 kHz", math.sqrt(K.G * K.M_sun / (9 * K.faltung(K.M_sun))**3) / (2 * math.pi) / 1e3, 6.2),
        ("01", "Wo es 2r₀ erreicht, liegt der Horizont der Falte", (lambda r0: K.c**2 / (K.G * 1e6 * K.M_sun / r0**2) / (2 * r0))(6 * K.faltung(1e6 * K.M_sun)), 1),
        ("01", "wo es 4D_E erreicht, der Einsteinring", (lambda D, M: K.c**2 / (K.G * M / (4 * K.G * M * D / K.c**2)) / (4 * D))(1.0 * K.Gpc, 1e11 * K.M_sun), 1),
        ("05", "q = e^{−r*/R} ≈ 0, Sonne r*/R = 1,6·10⁶", K.schwellenradius(K.M_sun) / K.R_SONNE / 1e6, 1.6),
        ("05", "δ(R₂₀₀) = 2,4·10⁻⁵ (X-COP)", __import__("statistics").median([v["tiefe_200"] for v in D.xcop_bilanz().values()]) / 1e-5, 2.4),
        ("08", "bis zur Naht, 5,8·r*, ist sie 0,025 der Höhe breit", K.naht_radius(K.M_sun, K.q_je_radiant()) / K.schwellenradius(K.M_sun), 5.8),
        ("08", "bis zur Naht, 5,8·r*, ist sie 0,025 der Höhe breit", math.log(2 * math.pi * K.naht_radius(K.M_sun, K.q_je_radiant()) / K.schwellenradius(K.M_sun)) / math.log(1 / K.leiterzahl()), 0.025),
        ("F:lesungen", "Mit θ̃ = H ist θ̃·τ_D = 2π/Θ = 6,16", K.atmen_heute()[0] * 1e3 / K.Mpc * K.rueckkehrzeit(), 6.16),
        ("F:lesungen", "Mit θ̃ = H ist θ̃·τ_D = 2π/Θ = 6,16", 2 * math.pi / K.theta_(H=K.atmen_heute()[0] * 1e3 / K.Mpc), 6.16),
        ("F:lesungen", "also ν ≤ 1,09 an jedem Ort", 1 / -math.expm1(-math.sqrt(K.atmen_heute()[0] * 1e3 / K.Mpc * K.rueckkehrzeit())), 1.09),
        ("F:lesungen", "lägen dann −0,01 dex auf der Rechnung", D.udg_reste(mitatmen=K.atmen_heute()[0] * 1e3 / K.Mpc)[1], -0.01),
        ("F:lesungen", "haben fünf der sechs weniger Nachbarn", float(sum(D.alfalfa_umgebung(a)[1] < 0.25 for a in (114905, 122966, 219533, 248945, 334315, 749290))), 5),
        ("F:lesungen", "Ihr Rang liegt im Mittel bei 0,27 statt 0,5, zwei Fehler darunter", sum(D.alfalfa_umgebung(a)[1] for a in (114905, 122966, 219533, 248945, 334315, 749290)) / 6, 0.27),
        ("F:lesungen", "Ihr Rang liegt im Mittel bei 0,27 statt 0,5, zwei Fehler darunter", (0.5 - sum(D.alfalfa_umgebung(a)[1] for a in (114905, 122966, 219533, 248945, 334315, 749290)) / 6) / math.sqrt(1 / 72), 2),
        ("F:lesungen", "in 1 Mpc bei 0,26", sum(D.alfalfa_umgebung(a, R=1.0)[1] for a in (114905, 122966, 219533, 248945, 334315, 749290)) / 6, 0.26),
        ("F:lesungen", "Nur AGC 219533 liegt dicht (0,87)", D.alfalfa_umgebung(219533)[1], 0.87),
        ("F:lesungen", "ihr W50 liegt bei gleicher HI-Masse +0,03 ± 0,02 dex wie das der übrigen", D.alfalfa_breite_umgebung()[1], 0.03),
        ("F:lesungen", "ihr W50 liegt bei gleicher HI-Masse +0,03 ± 0,02 dex wie das der übrigen", D.alfalfa_breite_umgebung()[2], 0.02),
        ("03", "In den Größen des Modells ist es a₀/(√κ·H·v_flach), mit R_ta³ = GM/(κ·H²)", K.inaktives_mass(6e10 * K.M_sun, 68.5e3 / K.Mpc, kappa=2.0) / (K.a0 / (math.sqrt(2.0) * 68.5e3 / K.Mpc * K.v_flach(6e10 * K.M_sun))), 1),
        ("03", "die Umlaufrate an der Schwelle, v_flach/r* = a₀/v_flach", K.v_flach(K.M_sun) / K.schwellenradius(K.M_sun) / (K.a0 / K.v_flach(K.M_sun)), 1),
        ("03", "also (c/v_flach)·Θ/(2π·√κ), tangential gegen parametrisch", K.inaktives_mass(1e8 * K.M_sun, 68.5e3 / K.Mpc) / ((K.c / K.v_flach(1e8 * K.M_sun)) * K.theta_(H=68.5e3 / K.Mpc) / (2 * math.pi)), 1),
        ("03", "also (c/v_flach)·Θ/(2π·√κ), tangential gegen parametrisch", K.inaktives_mass(1e8 * K.M_sun, 68.5e3 / K.Mpc, kappa=9 * math.pi**2 / (64 * math.pi)) / ((K.c / K.v_flach(1e8 * K.M_sun)) * K.theta_(H=68.5e3 / K.Mpc) / (2 * math.pi * math.sqrt(9 * math.pi**2 / (64 * math.pi)))), 1),
        ("03", "Kehrt er heute um, ist mit der Rückstellung κ = 0,3", K.groesste_umkehr(K.M_sun)[1]**2 / (K.umkehr_heute()[0]**3 * (1 + K.q_je_radiant()) * (K.atmen_heute()[0] * 1e3 / K.Mpc)**2), 0.3),
        ("03", "Kehrt er heute um, ist mit der Rückstellung κ = 0,3", K.umkehr_heute()[1] * (1 - math.pi**2 / 15 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2) / (2 * (1 + K.q_je_radiant())), 0.3),
        ("03", "mit Δ = (3π/4)² nach dem halben Umlauf, κ = 0,44", K.kollaps_kugel()[2] / 2 / (1 + K.q_je_radiant()), 0.44),
        ("03", "mit Δ = (3π/4)² nach dem halben Umlauf, κ = 0,44", (3 * math.pi / 4)**2 / (2 * 2 * math.pi), 0.44),
        ("03", "mit Δ = (3π/4)² nach dem halben Umlauf", K.kollaps_kugel()[2] / (3 * math.pi / 4)**2, 1),
        ("08", "bei M_min λ̄ = r̃ = ℓ₀ und r* = ℓ_P", K.schwelle_spiegel(K.schnittmassen()[0]) * K.a0 / K.c**2, 1),
        ("08", "bei M_min λ̄ = r̃ = ℓ₀ und r* = ℓ_P", K.schwellenradius(K.schnittmassen()[0]) / K.l_P, 1),
        ("08", "bei M_U 3Δs = r* = ℓ₀ und r̃ = ℓ_P", K.schwelle_spiegel(K.schnittmassen()[3]) / K.l_P, 1),
        ("08", "Das Dreieck schließt, α_G·α_I = (v_flach/c)⁴", float(max(abs((lambda k: k[0] * k[1])(K.kopplungen(m)) / (K.v_flach(m) / K.c)**4 - 1) for m in (K.m_e, K.m_p, K.M_sun, 1e11 * K.M_sun)) < 1e-12), 1),
        ("08", "α_I = λ̄/ℓ₀ = kT_I/(Mc²) = M_min/M", K.kopplungen(K.m_p)[1] / (K.hbar * K.a0 / K.c / (K.m_p * K.c**2)), 1),
        ("08", "α_I = λ̄/ℓ₀ = kT_I/(Mc²) = M_min/M", K.kopplungen(K.m_p)[1] / (K.schnittmassen()[0] / K.m_p), 1),
        ("08", "ist ihr Mittel √(α_G·α_I)", 3 * K.faltung(K.M_sun) / K.schwellenradius(K.M_sun) / math.sqrt(K.kopplungen(K.M_sun)[0] * K.kopplungen(K.M_sun)[1]), 1),
        ("08", "Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U", K.kopplungen(K.m_Pl)[0], 1),
        ("08", "Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U", K.kopplungen(K.schnittmassen()[0])[1], 1),
        ("08", "Jede Kopplung wird eins an einer Sprosse, α_G bei m_Pl, α_I bei M_min, Tully-Fisher bei M_U", K.kopplungen(K.schnittmassen()[3])[2], 1),
        ("08", "bei M_c sind beide Quantenkopplungen gleich, α_G = α_I = ε_L^{2/3}", K.kopplungen(K.schnittmassen()[1])[0] / K.leiterzahl()**(2 / 3), 1),
        ("08", "bei M_c sind beide Quantenkopplungen gleich, α_G = α_I = ε_L^{2/3}", K.kopplungen(K.schnittmassen()[1])[1] / K.leiterzahl()**(2 / 3), 1),
        ("08", "kehrt α_G um und tauscht die beiden anderen, α_I(m_Pl²/M) = (v_flach/c)⁴", K.kopplungen(K.m_Pl**2 / K.M_sun)[0] * K.kopplungen(K.M_sun)[0], 1),
        ("08", "kehrt α_G um und tauscht die beiden anderen, α_I(m_Pl²/M) = (v_flach/c)⁴", K.kopplungen(K.m_Pl**2 / K.M_sun)[1] / (K.v_flach(K.M_sun) / K.c)**4, 1),
        ("08", "Beim Elektron ist α_I = 4,6·10⁻⁴⁰", K.kopplungen(K.m_e)[1] / 1e-40, 4.6),
        ("04", "kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π)", K.hbar * K.a0 / K.c / (K.M_sun * K.c**2) / K.kopplungen(K.M_sun)[1], 1),
        ("04", "kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π)", K.falte_temperatur(K.M_sun) * K.k_B / (K.M_sun * K.c**2) * 8 * math.pi * K.kopplungen(K.M_sun)[0], 1),
        ("04", "kT_I/(Mc²) = α_I, bei der Falte kT/(Mc²) = 1/(8π·α_G) und T/T_I = (c/v_flach)⁴/(8π)", K.falte_temperatur(K.M_sun) * K.k_B / (K.hbar * K.a0 / K.c) / ((K.c / K.v_flach(K.M_sun))**4 / (8 * math.pi)), 1),
        ("04", "Ihre Ruheenergie erreicht die Falte bei M = m_Pl/√(8π) = 4,3·10⁻⁹ kg", (lambda M: K.falte_temperatur(M) * K.k_B / (M * K.c**2))(K.m_Pl / math.sqrt(8 * math.pi)), 1),
        ("04", "Ihre Ruheenergie erreicht die Falte bei M = m_Pl/√(8π) = 4,3·10⁻⁹ kg", K.m_Pl / math.sqrt(8 * math.pi) / 1e-9, 4.3),
        ("04", "Am Horizont ist das 4π·r₀: Die Wärme der Falte ist die Wärme ihrer Neigung am Horizont", (lambda M: 2 * math.pi * (6 * K.faltung(M))**2 / (3 * K.faltung(M)) * K.k_B * K.falte_temperatur(M) / (K.hbar * K.c))(K.M_sun), 1),
        ("04", "Am Horizont ist das 4π·r₀: Die Wärme der Falte ist die Wärme ihrer Neigung am Horizont", (lambda M: 2 * math.pi * K.c**2 / (K.G * M / (6 * K.faltung(M))**2) * K.k_B * K.falte_temperatur(M) / (K.hbar * K.c))(1e11), 1),
        ("04", "Die Falte erreicht diese Grenze, mit R = r₀ ist sie 4πk·α_G", (lambda M: 2 * math.pi * 6 * K.faltung(M) / K.compton(M) / K.entropie_falte(M))(K.M_sun), 1),
        ("04", "Die Falte erreicht diese Grenze, mit R = r₀ ist sie 4πk·α_G", (lambda M: 4 * math.pi * K.kopplung_grav(M) / K.entropie_falte(M))(1e40), 1),
        ("08", "Zwei gleiche Massen im Abstand a strahlen P = (64/5)·(c⁵/G)·δ⁵ mit δ = 3Δs/a", (lambda M, a: K.leistung_quadrupol(M, M, a) / (64 / 5 * K.c**5 / K.G * (3 * K.faltung(M) / a)**5))(30 * K.M_sun, 1e6), 1),
        ("08", "ihr Treffpunkt c⁵/G, die Ruheenergie je Lichtzeit der Tiefe, hängt an keiner Masse", (lambda M: M * K.c**2 / (3 * K.faltung(M) / K.c) / (K.c**5 / K.G))(1e3 * K.M_sun), 1),
        ("08", "Die Falte selbst strahlt als schwarzer Körper (c⁵/G)/(15360π·α_G)", (lambda M: M * K.c**2 / (3 * K.lebensdauer_falte(M)) / (K.c**5 / K.G / (15360 * math.pi * K.kopplung_grav(M))))(1e11), 1),
        ("08", "lebt sie t = 5120π·α_G·3Δs/c, eine Sonnenmasse 2,1·10⁶⁷ Jahre", (lambda M: K.lebensdauer_falte(M) * K.c / (5120 * math.pi * K.kopplung_grav(M) * 3 * K.faltung(M)))(1e11), 1),
        ("08", "lebt sie t = 5120π·α_G·3Δs/c, eine Sonnenmasse 2,1·10⁶⁷ Jahre", K.lebensdauer_falte(K.M_sun) / K.yr / 1e67, 2.1),
        ("08", "Die Länge c·t ist der Punkt (3, 0) der Leiter, um ℓ_P gespiegelt zu λ̄/α_G (−3, 0)", (lambda M: K.c * K.lebensdauer_falte(M) * K.compton(M) / K.kopplung_grav(M) / (5120 * math.pi * K.l_P**2))(K.M_sun), 1),
        ("08", "Die Länge c·t ist der Punkt (3, 0) der Leiter, um ℓ_P gespiegelt zu λ̄/α_G (−3, 0)", (lambda M: K.c * K.lebensdauer_falte(M) / K.leiter_laenge(3, 0, M) / (5120 * math.pi))(1e20), 1),
        ("03", "S_Λ·Λ_s·ℓ_P² = 3πk wie bei de Sitter", 45 / math.pi**3 * (K.schliesslaenge() / K.l_P)**2 * math.pi**4 / (15 * K.schliesslaenge()**2) * K.l_P**2 / (3 * math.pi), 1),
        ("03", "S_Λ = (45/4π³)·k/ε_L² = 1,0·10¹²³ k", 45 / (4 * math.pi**3) / K.leiterzahl()**2 / 1e123, 1.0),
        ("03", "S_Λ = (45/4π³)·k/ε_L² = 1,0·10¹²³ k", 45 / math.pi**3 * (K.schliesslaenge() / K.l_P)**2 / 1e123, 1.0),
        ("09", "Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c", K.c / K.compton(K.schnittmassen()[0]) * K.rueckkehrzeit(), 1),
        ("09", "Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c", K.c / (3 * K.faltung(K.schnittmassen()[3])) * K.rueckkehrzeit(), 1),
        ("09", "Sie trifft den Compton-Takt bei M_min, den Faltungstakt bei M_U und die Faltungsrate bei M_c", K.kopplungen(K.schnittmassen()[1])[0] * K.c / K.compton(K.schnittmassen()[1]) * K.rueckkehrzeit(), 1),
        ("08", "c³τ_D/G = 1,1·10⁵⁴ kg", K.schnittmassen()[3] / 1e54, 1.1),
        ("08", "Δs(M_U) = ⅔·Λ", K.faltung(K.schnittmassen()[3]) / K.schliesslaenge() * 3, 2),
        ("08", "(4π⁴/15)·ε_L² = 0,94·10⁻¹²²", 4 * math.pi**4 / 15 * K.leiterzahl()**2 / 1e-122, 0.94),
        ("02", "r* = 32 fm", K.schwellenradius(K.m_p) * 1e15, 32),
        ("02", "a_B = 1631·r*", K.bohr_radius() / K.schwellenradius(K.m_p), 1631),
        ("04", "Friedmanns H² = 8πG·ρ_krit/3 sagt, dass diese Kugel ihre eigene Falte ist, 6Δs_H = c/H", (lambda H: 6 * K.faltung(4 * math.pi / 3 * 3 * H**2 / (8 * math.pi * K.G) * (K.c / H)**3) / (K.c / H))(K.H0), 1),
        ("04", "Friedmanns H² = 8πG·ρ_krit/3 sagt, dass diese Kugel ihre eigene Falte ist, 6Δs_H = c/H", (lambda H: K.faltung(4 * math.pi / 3 * 3 * H**2 / (8 * math.pi * K.G) * (K.c / H)**3) / (K.c / (6 * H)))(2 * K.H0), 1),
        ("04", "Sonne 1,05·10⁷⁷ k", K.entropie_falte(K.M_sun) / 1e77, 1.05),
        ("01", "0,022 s", (1 - K.takt(K.M_earth, R_ERDE)) * K.yr, 0.022),
        ("01", "Sonnenoberfläche 67 s", (1 - K.takt(K.M_sun, R_SONNE)) * K.yr, 67),
        ("01", "auf der Erdbahn sind es 0,31 s", 3 * K.faltung(K.M_sun) / K.AU * K.yr, 0.31),
        ("01", "Takt 0,81", K.takt(1.4 * K.M_sun, 12e3), 0.81),
        ("08", "1,7512″", math.degrees(K.ablenkung(K.M_sun, R_SONNE)) * 3600, 1.7512),
        ("F:ganzes", "γ − 1 = (2,1 ± 2,3)·10⁻⁵", K.CASSINI_GAMMA[0] / 1e-5, 2.1),
        ("F:ganzes", "(2,1 ± 2,3)·10⁻⁵", K.CASSINI_GAMMA[1] / 1e-5, 2.3),
        ("F:ganzes", "(gemessen 43,11 ± 0,45″)", K.MERKUR_GEMESSEN[0], 43.11),
        ("F:ganzes", "43,11 ± 0,45″)", K.MERKUR_GEMESSEN[1], 0.45),
        ("F:ganzes", "geht 4,46·10⁻¹⁰ vor", K.uhr_bahn(K.M_earth, K.R_ERDE, K.GPS_BAHN) / 1e-10, 4.46),
        ("F:ganzes", "38,5 µs am Tag", K.uhr_bahn(K.M_earth, K.R_ERDE, K.GPS_BAHN) * 86400e6, 38.5),
        ("F:ganzes", "| 4,465·10⁻¹⁰ (Werkseinstellung", K.GPS_WERK / 1e-10, 4.465),
        ("F:ganzes", "| 6605 mas/a |", K.kreisel_bahn(K.M_earth, K.R_GPB) * K.yr * 180 / math.pi * 3600e3, 6605),
        ("F:ganzes", "6602 ± 18 mas/a", K.GPB_KREISEL[0], 6602),
        ("F:ganzes", "6602 ± 18", K.GPB_KREISEL[1], 18),
        ("F:ganzes", "| 40,8 mas/a |", K.mitfuehrung(K.J_ERDE, K.R_GPB) * K.yr * 180 / math.pi * 3600e3, 40.8),
        ("F:ganzes", "37,2 ± 7,2 mas/a", K.GPB_MITFUEHRUNG[0], 37.2),
        ("F:ganzes", "37,2 ± 7,2", K.GPB_MITFUEHRUNG[1], 7.2),
        ("F:ganzes", "Ṗ_b = −2,40·10⁻¹²", K.doppelpulsar(K.PSR_1913["m1"] * K.M_sun, K.PSR_1913["m2"] * K.M_sun, K.PSR_1913["Pb"], K.PSR_1913["e"]) / 1e-12, -2.40),
        ("F:ganzes", "0,9983 ± 0,0016 davon", K.PSR_1913["verhaeltnis"][0], 0.9983),
        ("F:ganzes", "0,9983 ± 0,0016", K.PSR_1913["verhaeltnis"][1], 0.0016),
        ("F:ganzes", "| 39,7 µas |", K.schatten(K.M87["M"] * K.M_sun, K.M87["D"] * K.Mpc) * 180 / math.pi * 3600e6, 39.7),
        ("F:ganzes", "42 ± 3 µas", K.M87["ring"][0], 42),
        ("F:ganzes", "gleich c auf 10⁻¹⁵", float(max(abs(x) for x in K.GW170817_DC) <= 3e-15), 1.0),
        ("07", "½, ¾ und ⅞ für d = 1, 2, 3", sum(abs(K.strahlung_je_sinn(2, d) / K.strahlung_je_sinn(1, d) - w) < 1e-9 for d, w in ((1, 0.5), (2, 0.75), (3, 0.875))), 3),
        ("07", "∫n_F ds = ln 2 = 0,693", K.strahlung_je_sinn(2, 0), 0.693),
        ("02", "alle 884 gerade-geraden Kerne haben 0⁺", float(KE.paare_grundspin()[1] == KE.paare_grundspin()[0]) * KE.paare_grundspin()[0], 884),
        ("02", "sind 849 von 850 Spins halbzahlig", KE.paare_grundspin()[3], 849),
        ("02", "849 von 850 Spins", KE.paare_grundspin()[2], 850),
        ("02", "bei den 352 ungerade-ungeraden alle ganzzahlig", float(KE.paare_grundspin()[5] == KE.paare_grundspin()[4]) * KE.paare_grundspin()[4], 352),
        ("F:ganzes", "(11/4)^{1/3} = 1,401", K.licht_gegen_neutrinos(), 1.401),
        ("F:ganzes", "T_ν = 1,945 K", K.T_CMB / K.licht_gegen_neutrinos(), 1.945),
        ("07", "∏(1 + qⁿ), tief ln Z = π²/(12s) − ½·ln 2 + s/24", float(all(abs(K.alle_vielfachen(s, z)[0] - K.alle_vielfachen(s, z)[1]) < 1e-9 for s in (0.05, 0.2, 0.5) for z in (1, 2))), 1.0),
        ("02", "positiv, 10,7·A^{−1/2} MeV", KE.anpassung(KE.NIVEAUS)[0][4], 10.7),
        ("02", "Stabil sind nach NUBASE2020 253 Kerne", float(KE.stabile_nach_paaren()[0]), 253),
        ("02", "253 Kerne, 147 gerade-gerade", float(KE.stabile_nach_paaren()[1]), 147),
        ("02", "102 mit ungeradem A", float(KE.stabile_nach_paaren()[2]), 102),
        ("02", "nur 4 ungerade-ungerade", float(KE.stabile_nach_paaren()[3] * KE.stabile_nach_paaren()[4]), 4),
        ("07", "(Euler, bei q = ½ beide 2,384)", K.sperre_und_ungerade(0.5)[0], 2.384),
        ("07", "bei q = ½ beide 2,384)", K.sperre_und_ungerade(0.5)[1], 2.384),
        ("01", "Abgabe der Windung quer zur Schwingung 2/3", K.kugelmittel_quer(), 2 / 3),
        ("F:ganzes", "denselben vollen Umlauf geben, ν = 2π", K.verstaerkung(K.a0 * math.log((1 + K.q_je_radiant()) / K.q_je_radiant())**2) / (2 * math.pi), 1.0),
        ("F:ganzes", "d/π = 0,955", 6 / (1 + K.q_je_radiant()), 0.955),
        ("F:ganzes", "√(d/π) = 0,977", math.sqrt(6 / (1 + K.q_je_radiant())), 0.977),
        ("F:ganzes", "δ = 2·(ln q)⁴·Θ²/Δ = 9,4·10⁻⁶ heute", (lambda H: (lambda R, M: K.G * M / (R * K.c**2))(*K.naht_haufen(200, H, K.q_je_radiant())))(K.atmen_heute()[0] * 1e3 / K.Mpc) * 1e6, 9.4),
        ("F:ganzes", "δ = 2·(ln q)⁴·Θ²/Δ = 9,4·10⁻⁶ heute", 2 * math.log(1 - 1 / (2 * math.pi))**4 * K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)**2 / 200 * 1e6, 9.4),
        ("F:ganzes", "Θ_∞/2π = √15/(2π²) = 0,196", math.sqrt(15) / (2 * math.pi**2), 0.196),
        ("F:ganzes", "Θ_∞/2π = √15/(2π²) = 0,196", K.a0 / (K.groesste_umkehr(K.M_sun)[1] * K.c), 0.196),
        ("F:ganzes", "für die Milchstraße 1759·0,196 = 345", K.c / (K.G * 6e10 * K.M_sun * K.a0)**0.25, 1759),
        ("F:ganzes", "für die Milchstraße 1759·0,196 = 345", K.c / (K.G * 6e10 * K.M_sun * K.a0)**0.25 * math.sqrt(15) / (2 * math.pi**2), 345),
        ("F:ganzes", "für die Milchstraße 1759·0,196 = 345", K.groesste_umkehr(6e10 * K.M_sun)[0] / K.schwellenradius(6e10 * K.M_sun), 345),
        ("F:ganzes", "Spinzahl weit über 2π/Θ_∞ = 5,1", 2 * math.pi**2 / math.sqrt(15), 5.1),
        ("F:ganzes", "Bei 1 + Q = 6·(9/8)² = 7,59 ist sie 4/9", 6 * (9 / 8)**2, 7.59),
        ("F:ganzes", "Bei 1 + Q = 6·(9/8)² = 7,59 ist sie 4/9", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda x: 0.5 * math.sqrt(6 / x) - 4 / 9, 6, 20), 7.59),
        ("F:ganzes", "Im Fenster 5 < Q < 6,59 dazwischen", 6 * (9 / 8)**2 - 1, 6.59),
        ("F:ganzes", "das gemessene Verhältnis, 1 + Q = 6,36", 1 + K.q_messung()[0], 6.36),
        ("F:ganzes", "ganz bleibt nur der Knoten, beim Anteil 0,080 der Masse", 0.5 - math.atan(math.sqrt(5)) / (2 * math.pi) - math.sqrt(5) / (3 * math.pi), 0.080),
        ("F:galaxien", "heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11", K.scheibe_achsen(50e3, K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)), 0.064),
        ("F:galaxien", "heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11", K.scheibe_achsen(50e3, K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc), numerisch=True), 0.064),
        ("F:galaxien", "heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11", K.scheibe_achsen(250e3, K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc)), 0.11),
        ("F:galaxien", "heute bei v_flach = 50 bis 250 km/s 0,064 bis 0,11", K.scheibe_achsen(250e3, K.theta_(K.a0, K.atmen_heute()[0] * 1e3 / K.Mpc), numerisch=True), 0.11),
        ("F:galaxien", "wächst R_d wie M^{5/12 − 0,11} = M^{0,30}", 5 / 12 + _sigma_neigung(), 0.30),
        ("F:galaxien", "wächst R_d wie M^{5/12 − 0,11} = M^{0,30}", _sigma_neigung(), -0.11),
        ("F:galaxien", "an SPARC wachsen die Scheiben wie M^{0,31 ± 0,02}", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[2], 0.31),
        ("F:galaxien", "an SPARC wachsen die Scheiben wie M^{0,31 ± 0,02}", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[3], 0.02),
        ("F:ganzes", "atmet die Sonne in 2π·√(R³/GM) = 2,8 h", math.tau / math.sqrt(K.atmen_gehalten(((3, K.GM_sun / K.R_SONNE**3),))) / 3600, 2.8),
        ("F:ganzes", "atmet die Sonne in 2π·√(R³/GM) = 2,8 h", math.tau * math.sqrt(K.R_SONNE**3 / K.GM_sun) / 3600, 2.8),
        ("F:ganzes", "mit κ = 1,26·Ω in 0,20 Gyr", math.sqrt(K.flussdimension_bahn(((3 - K.schalter(K.a0), 1.0),))[1]), 1.26),
        ("F:ganzes", "mit κ = 1,26·Ω in 0,20 Gyr", math.sqrt(1 + 1 / (math.e - 1)), 1.26),
        ("F:ganzes", "mit κ = 1,26·Ω in 0,20 Gyr", (lambda M: (lambda rs: math.tau / math.sqrt(K.atmen_gehalten(((3 - K.schalter(K.G * M / rs**2), (K.bahntempo_regel(M, rs) / rs)**2),))) / K.Gyr)(K.schwellenradius(M)))(6e10 * K.M_sun), 0.2),
        ("F:ganzes", "mit κ = 1,26·Ω in 0,20 Gyr", (lambda M: (lambda rs: math.tau * rs / (K.bahntempo_regel(M, rs) * math.sqrt(1 + 1 / (math.e - 1))) / K.Gyr)(K.schwellenradius(M)))(6e10 * K.M_sun), 0.2),
        ("F:ganzes", "K = 100·H₀² − 4·H_∞², in 9,3 Gyr", math.tau / math.sqrt(K.atmen_gehalten(((3, 100 * (K.atmen_heute()[0] * 1e3 / K.Mpc)**2), (0, -K.groesste_umkehr(K.M_sun)[1]**2)))) / K.Gyr, 9.3),
        ("F:ganzes", "K = 100·H₀² − 4·H_∞², in 9,3 Gyr", math.tau / math.sqrt(100 * (K.atmen_heute()[0] * 1e3 / K.Mpc)**2 - 4 * math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie())**2) / K.Gyr, 9.3),
        ("F:ganzes", "= −3·H_∞² für jeden Träger", K.atmen_ganzes(1.0), -3),
        ("F:ganzes", "= −3·H_∞² für jeden Träger", K.atmen_ganzes(2.0), -3),
        ("F:ganzes", "= −3·H_∞² für jeden Träger", K.atmen_ganzes(numerisch=True), -3),
        ("F:ganzes", "(1/(√3·H_∞) = 10,2 Gyr)", 1 / (math.sqrt(3) * K.groesste_umkehr(K.M_sun)[1]) / K.Gyr, 10.2),
        ("F:ganzes", "(1/(√3·H_∞) = 10,2 Gyr)", 1 / (math.sqrt(3) * math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie())) / K.Gyr, 10.2),
        ("F:ganzes", "H_∞ = 2π²a₀/(√15·c) = 55,6 km/s/Mpc", 2 * math.pi**2 * K.a0 / (math.sqrt(15) * K.c) * K.Mpc / 1e3, 55.6),
        ("F:ganzes", "H_∞ = 2π²a₀/(√15·c) = 55,6 km/s/Mpc", math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) * K.Mpc / 1e3, 55.6),
        ("F:ganzes", "1/H_∞ = 17,6 Gyr", 1 / K.groesste_umkehr(K.M_sun)[1] / K.Gyr, 17.6),
        ("F:ganzes", "1/H_∞ = 17,6 Gyr", 1 / math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) / K.Gyr, 17.6),
        ("F:ganzes", "Heute atmet der Formteil mit H₀ = 67,2", K.atmen_heute()[0], 67.2),
        ("F:ganzes", "Heute atmet der Formteil mit H₀ = 67,2", math.sqrt(8 * math.pi * K.G / 3 * (K.dichte_heute()[0] + K.dichte_heute()[2] + K.dichte_dunkle_energie())) * K.Mpc / 1e3, 67.2),
        ("F:ganzes", "τ_D = c/a₀ = 89,6 Gyr für jede Masse", K.rueckkehrzeit() / K.Gyr, 89.6),
        ("F:ganzes", "τ_D = c/a₀ = 89,6 Gyr für jede Masse", K.c / K.a0 / K.Gyr, 89.6),
        ("F:ganzes", "H_∞·τ_D = 2π/Θ_∞ = 5,10", K.groesste_umkehr(K.M_sun)[1] * K.rueckkehrzeit(), 5.1),
        ("F:ganzes", "H_∞·τ_D = 2π/Θ_∞ = 5,10", 2 * math.pi / (math.sqrt(15) / math.pi), 5.1),
        ("F:ganzes", "Eine Rückkehr dauert 5,1 Zeiten des Atmens", math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) * K.c / K.a0, 5.1),
        ("F:ganzes", "um den Faktor 163", math.exp(K.groesste_umkehr(K.M_sun)[1] * K.rueckkehrzeit()), 163),
        ("F:ganzes", "um den Faktor 163", math.exp(2 * math.pi**2 / math.sqrt(15)), 163),
        ("F:ganzes", "gleichwertig Θ_∞/2π = 0,196", math.sqrt(15) / (2 * math.pi**2), 0.196),
        ("F:ganzes", "gleichwertig Θ_∞/2π = 0,196", K.a0 / math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) / K.c, 0.196),
        ("F:ganzes", "hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22", _puls_hubble()[0], 1.8),
        ("F:ganzes", "hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22", _puls_hubble()[4], 1.8),
        ("F:ganzes", "hebt er H gegen den Formteil um höchstens 1,8 % bei z = 0,22", _puls_hubble()[1], 0.22),
        ("F:ganzes", "heute um 1,1 %, und H fällt dabei stets", _puls_hubble()[2], 1.1),
        ("F:ganzes", "heute um 1,1 %, und H fällt dabei stets", _puls_hubble()[3], 1),
        ("F:ganzes", "ρ_T = 3H_∞²/(4πG) = 1,2·10⁻²⁶ kg/m³", K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2) / 1e-26, 1.2),
        ("F:ganzes", "ρ_T = 3H_∞²/(4πG) = 1,2·10⁻²⁶ kg/m³", 2 * K.dichte_dunkle_energie() / 1e-26, 1.2),
        ("F:ganzes", "sieben Wasserstoffatome je Kubikmeter", K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2) / K.m_H, 7),
        ("F:ganzes", "sieben Wasserstoffatome je Kubikmeter", 2 * K.dichte_dunkle_energie() / K.m_H, 7),
        ("F:ganzes", "die Sonne mit ρ̄/ρ_T = 1,2·10²⁹", K.GM_sun / K.G / (4 / 3 * math.pi * K.R_SONNE**3) / K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2) / 1e29, 1.2),
        ("F:ganzes", "die Sonne mit ρ̄/ρ_T = 1,2·10²⁹", K.GM_sun / K.R_SONNE**3 / K.groesste_umkehr(K.M_sun)[1]**2 / 1e29, 1.2),
        ("F:ganzes", "ein Haufen in R₂₀₀ heute mit 146", 100 * ((K.atmen_heute()[0] * 1e3 / K.Mpc) / K.groesste_umkehr(K.M_sun)[1])**2, 146),
        ("F:ganzes", "ein Haufen in R₂₀₀ heute mit 146", 200 * 3 * (K.atmen_heute()[0] * 1e3 / K.Mpc)**2 / (8 * math.pi * K.G) / (2 * K.dichte_dunkle_energie()), 146),
        ("F:ganzes", "im Formteil bei z = 0,63", (2 * K.dichte_dunkle_energie() / K.dichte_heute()[2])**(1 / 3) - 1, 0.63),
        ("F:ganzes", "im Formteil bei z = 0,63", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda z: K.dichte_heute()[2] * (1 + z)**3 + 2 * K.dichte_heute()[0] * (1 + z)**4 - 2 * K.dichte_dunkle_energie(), 0, 3), 0.63),
        ("F:ganzes", "kT = ħH/(2π), heute 2,6·10⁻³⁰ K", K.hbar * (K.atmen_heute()[0] * 1e3 / K.Mpc) / (2 * math.pi * K.k_B) / 1e-30, 2.6),
        ("F:ganzes", "kT = ħH/(2π), heute 2,6·10⁻³⁰ K", K.hbar * K.c / (12 * math.pi * K.k_B * K.c / (6 * (K.atmen_heute()[0] * 1e3 / K.Mpc))) / 1e-30, 2.6),
        ("F:ganzes", "im Endzustand 2,2·10⁻³⁰ K", K.hbar * K.groesste_umkehr(K.M_sun)[1] / (2 * math.pi * K.k_B) / 1e-30, 2.2),
        ("F:ganzes", "im Endzustand 2,2·10⁻³⁰ K", K.hbar * K.c / (12 * math.pi * K.k_B * K.c / (6 * K.groesste_umkehr(K.M_sun)[1])) / 1e-30, 2.2),
        ("F:ganzes", "1/(30·S_Λ/k) = 3,3·10⁻¹²⁵", math.pi**2 / 15 * (K.hbar * K.groesste_umkehr(K.M_sun)[1] / (2 * math.pi))**4 / (K.hbar * K.c)**3 / K.c**2 / K.dichte_dunkle_energie() / 1e-125, 3.3),
        ("F:ganzes", "1/(30·S_Λ/k) = 3,3·10⁻¹²⁵", 1 / (30 * 45 / math.pi**3 * (K.schliesslaenge() / K.l_P)**2) / 1e-125, 3.3),
        ("F:ganzes", "kT_I = ħa₀/c = 2,7·10⁻³⁰ K", K.hbar * K.a0 / (K.c * K.k_B) / 1e-30, 2.7),
        ("F:ganzes", "kT_I = ħa₀/c = 2,7·10⁻³⁰ K", K.hbar * K.c / (K.k_B * 2 * K.schliesslaenge()) / 1e-30, 2.7),
        ("F:ganzes", "Heute ist das Atmen 2 % kühler als die Rückkehr", 100 * (1 - (K.hbar * (K.atmen_heute()[0] * 1e3 / K.Mpc) / (2 * math.pi)) / (K.hbar * K.a0 / K.c)), 2),
        ("F:ganzes", "Heute ist das Atmen 2 % kühler als die Rückkehr", 100 * (1 - 1 / K.theta_(K.a0, (K.atmen_heute()[0] * 1e3 / K.Mpc))), 2),
        ("F:ganzes", "im Endzustand um den Faktor Θ_∞ = 1,23", math.sqrt(15) / math.pi, 1.23),
        ("F:ganzes", "im Endzustand um den Faktor Θ_∞ = 1,23", (K.hbar * K.a0 / K.c) / (K.hbar * K.groesste_umkehr(K.M_sun)[1] / (2 * math.pi)), 1.23),
        ("F:ganzes", "Heute liegt das Mittel bei 0,23·ρ_T", K.dichte_heute()[2] / (2 * K.dichte_dunkle_energie()), 0.23),
        ("F:ganzes", "Heute liegt das Mittel bei 0,23·ρ_T", K.dichte_heute()[2] / K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2), 0.23),
        ("F:ganzes", "erst ab dem 17-fachen des Mittels", 8 * K.dichte_dunkle_energie() / K.dichte_heute()[2], 17),
        ("F:ganzes", "erst ab dem 17-fachen des Mittels", 4 * K.gezeitendichte(K.groesste_umkehr(K.M_sun)[1]**2) / K.dichte_heute()[2], 17),
        ("F:lesungen", "im besten Punkt α = 3,2°", _kreise_lesung()["alpha"], 3.2),
        ("F:lesungen", "im besten Punkt α = 3,2°", _zelle_kreis()[2], 3.2),
        ("F:lesungen", "gibt das beste Paar S = 0,78", _kreise_lesung()["paar"], 0.78),
        ("F:lesungen", "die beste Dreiheit senkrechter Paare S₃ = 0,58", _kreise_lesung()["drei"], 0.58),
        ("F:lesungen", "gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63", _kreise_lesung()["null"][0], 0.76),
        ("F:lesungen", "gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63", _kreise_lesung()["null"][1], 0.82),
        ("F:lesungen", "gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63", _kreise_lesung()["drei_null"][0], 0.58),
        ("F:lesungen", "gibt dasselbe, 0,76 bis 0,82 und 0,58 bis 0,63", _kreise_lesung()["drei_null"][1], 0.63),
        ("F:lesungen", "Ein Kreis von 20° Umfang", 360 * math.sin(math.radians(_kreise_lesung()["alpha"])), 20),
        ("F:lesungen", "Streuung σ = 0,23 je Paar nur 1/σ² = 18 unabhängige Flecken", _kreise_lesung()["streu"], 0.23),
        ("F:lesungen", "Streuung σ = 0,23 je Paar nur 1/σ² = 18 unabhängige Flecken", 1 / _kreise_lesung()["streu"]**2, 18),
        ("F:lesungen", "am ersten Gipfel (ℓ = 220) 11 % von C_ℓ, im Tal (ℓ = 400) 47 %", 100 * _kreise_lesung()["doppler"][0], 11),
        ("F:lesungen", "am ersten Gipfel (ℓ = 220) 11 % von C_ℓ, im Tal (ℓ = 400) 47 %", 100 * _kreise_lesung()["doppler"][1], 47),
        ("F:lesungen", "die Zelle erwartet S = 0,25 bis 0,37", _kreise_lesung()["erwartet"][0], 0.25),
        ("F:lesungen", "die Zelle erwartet S = 0,25 bis 0,37", _kreise_lesung()["erwartet"][1], 0.37),
        ("F:lesungen", "liest die Messung an den Achsen S = 0,25 ± 0,06", _kreise_lesung()["himmel"][0], 0.25),
        ("F:lesungen", "liest die Messung an den Achsen S = 0,25 ± 0,06", _kreise_lesung()["himmel"][1], 0.06),
        ("F:lesungen", "die Zerlegung dort 0,28 bis 0,37", _kreise_lesung()["zerlegung20"][0], 0.28),
        ("F:lesungen", "die Zerlegung dort 0,28 bis 0,37", _kreise_lesung()["zerlegung20"][1], 0.37),
        ("F:ganzes", "(CAMB, x_e = 2,9·10⁻⁴)", K.signal_21cm()["x_e"] / 1e-4, 2.9),
        ("F:ganzes", "trägt das Gas die Wärme des Lichts bis z_t = 123", K.signal_21cm()["z_t"], 123),
        ("F:ganzes", "in der geschlossenen Form 122", K.kuehlung_licht(x_e=K.signal_21cm()["x_e"])["z_t_form"], 122),
        ("F:ganzes", "z = 17 auf 7,1 K gegen 49,1 K des Lichts", K.signal_21cm()["T_scharf"], 7.1),
        ("F:ganzes", "z = 17 auf 7,1 K gegen 49,1 K des Lichts", K.signal_21cm()["T_licht"], 49.1),
        ("F:ganzes", "z = 17 auf 7,1 K gegen 49,1 K des Lichts", K.T_CMB * 18, 49.1),
        ("F:ganzes", "mit dem Übergang y 7,3 K, in der vollen Rechnung von CAMB 6,9 K", K.signal_21cm()["T_y"], 7.3),
        ("F:ganzes", "mit dem Übergang y 7,3 K, in der vollen Rechnung von CAMB 6,9 K", K.signal_21cm()["T_camb"], 6.9),
        ("F:ganzes", "A₁₀ = μ₀ω³μ_B²/(3πħc³) = 2,87·10⁻¹⁵ s⁻¹", K.signal_21cm()["A10"] / 1e-15, 2.87),
        ("F:ganzes", "A₁₀ = μ₀ω³μ_B²/(3πħc³) = 2,87·10⁻¹⁵ s⁻¹", 16 * math.pi**3 * K.mu0 * K.NU_21CM**3 * (K.e * K.hbar / (2 * K.m_e))**2 / (3 * K.h * K.c**3) / 1e-15, 2.87),
        ("F:ganzes", "bei z = 17 (79 MHz) τ = 0,09", K.signal_21cm()["MHz"], 79),
        ("F:ganzes", "bei z = 17 (79 MHz) τ = 0,09", K.signal_21cm()["tau"], 0.09),
        ("F:ganzes", "das Signal gegen das Licht −0,20 K (−196 bis −210 mK)", K.signal_21cm()["dT_scharf"] / 1e3, -0.2),
        ("F:ganzes", "das Signal gegen das Licht −0,20 K (−196 bis −210 mK)", K.signal_21cm()["dT_y"], -196),
        ("F:ganzes", "das Signal gegen das Licht −0,20 K (−196 bis −210 mK)", K.signal_21cm()["dT_camb"], -210),
        ("F:ganzes", "EDGES liest bei 78 MHz eine Aufnahme von 0,5 K", D.signal_21cm_messung()["mitte"], 78),
        ("F:ganzes", "EDGES liest bei 78 MHz eine Aufnahme von 0,5 K", D.signal_21cm_messung()["tiefe"], 0.5),
        ("F:ganzes", "das 2,4- bis 2,5-fache dessen, was das Modell dort zulässt", 1e3 * D.signal_21cm_messung()["tiefe"] / abs(K.signal_21cm()["dT_camb"]), 2.4),
        ("F:ganzes", "das 2,4- bis 2,5-fache dessen, was das Modell dort zulässt", 1e3 * D.signal_21cm_messung()["tiefe"] / abs(K.signal_21cm()["dT_y"]), 2.5),
        ("F:ganzes", "SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 %", D.signal_21cm_messung()["band"][0], 55),
        ("F:ganzes", "SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 %", D.signal_21cm_messung()["band"][1], 85),
        ("F:ganzes", "SARAS 3 (55 bis 85 MHz) verwirft dieses Profil im selben Band mit 95,3 %", D.signal_21cm_messung()["verworfen"], 95.3),
        ("F:ganzes", "bei x = 0,25000", __import__("kerne").ende_bindung(False)[2], 0.25),
        ("F:ganzes", "bei x = 0,25000", 1 / (2 * (3 - 1)), 0.25),
        ("F:ganzes", "senkt es um 2Z/(2Z + 1) auf 0,245", __import__("kerne").ende_bindung()[2], 0.245),
        ("F:ganzes", "senkt es um 2Z/(2Z + 1) auf 0,245", (lambda Z: 0.25 * 2 * Z / (2 * Z + 1))(__import__("kerne").ende_bindung()[1]), 0.245),
        ("F:ganzes", "Bei z = 10 fehlen 0,29 %", 100 * (1 - K.schwelle_in_der_zeit(10)), 0.29),
        ("F:ganzes", "Bei z = 10 fehlen 0,29 %", 100 * (1 - K.schwelle_in_der_zeit(10, photon=True)), 0.29),
        ("F:ganzes", "ein Prozent erst bei z = 34", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda z: 0.99 - K.schwelle_in_der_zeit(z), 1, 100), 34),
        ("F:ganzes", "ein Prozent erst bei z = 34", __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda z: 0.99 - K.schwelle_in_der_zeit(z, photon=True), 1, 100), 34),
        ("F:ganzes", "an der letzten Streufläche 24 %", 100 * (1 - K.schwelle_in_der_zeit(K.z_stern)), 24),
        ("F:ganzes", "an der letzten Streufläche 24 %", 100 * (1 - K.schwelle_in_der_zeit(K.z_stern, photon=True)), 24),
        ("F:ganzes", "dort trug der Vakuumteil das 0,57-fache von heute", K.schwelle_in_der_zeit(K.z_stern)**2, 0.57),
        ("F:ganzes", "dort trug der Vakuumteil das 0,57-fache von heute", K.schwelle_in_der_zeit(K.z_stern, photon=True)**2, 0.57),
        ("F:ganzes", "Seit z ≈ 10 ist die Schwelle so auf 0,3 % fest", 100 * (1 - K.schwelle_in_der_zeit(10)), 0.3),
        ("F:ganzes", "noch um a_eq = 0,03 %", 100 * (K.schwelle_in_der_zeit(-1) - 1), 0.03),
        ("F:ganzes", "noch um a_eq = 0,03 %", 100 / (1 + K.z_gleichheit()), 0.03),
        ("F:ganzes", "erst bei v_flach = π²c/√15 = 2,55·c", math.pi**2 / math.sqrt(15), 2.55),
        ("F:ganzes", "erst bei v_flach = π²c/√15 = 2,55·c", K.c**2 / (2 * K.a0) * math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) / K.c, 2.55),
        ("F:ganzes", "trägt die mittlere Achse 19/30 von ⟨J²⟩", 30 * K.netz_drehung()[1], 19),
        ("F:ganzes", "trägt die mittlere Achse 19/30 von ⟨J²⟩", 30 * K.netz_drehung(numerisch=True)[1], 19),
        ("F:ganzes", "jede äußere 11/60", 60 * K.netz_drehung()[0], 11),
        ("F:ganzes", "jede äußere 11/60", 60 * K.netz_drehung(numerisch=True)[0], 11),
        ("F:ganzes", "jede äußere 11/60", 60 * K.netz_drehung(numerisch=True)[2], 11),
        ("F:ganzes", "⟨(λ₁ − λ₃)²⟩ = (19/15)·σ²", 15 * K.netz_drehung()[3], 19),
        ("F:ganzes", "⟨(λ₁ − λ₃)²⟩ = (19/15)·σ²", 15 * K.netz_drehung(numerisch=True)[3], 19),
        ("F:ganzes", "gegen je (11/30)·σ²", 30 * (2 - K.netz_drehung()[3]) / 2, 11),
        ("F:ganzes", "die Schwelle des Bereichs, δ_c = 1,686", K.zwei_schwellen()[2], 1.686),
        ("F:ganzes", "die des Leerraums, δ_v = −2,72", K.zwei_schwellen()[3], -2.72),
        ("F:ganzes", "mit δ_c/(δ_c + |δ_v|) = 0,38 in einem Leerraum", K.zwei_schwellen()[0], 0.38),
        ("F:ganzes", "mit δ_c/(δ_c + |δ_v|) = 0,38 in einem Leerraum", K.zwei_schwellen(numerisch=True)[0], 0.38),
        ("F:ganzes", "und mit 0,62 in einem Bereich", K.zwei_schwellen()[1], 0.62),
        ("F:ganzes", "und mit 0,62 in einem Bereich", K.zwei_schwellen(numerisch=True)[1], 0.62),
        ("F:ganzes", "Mit 1 + δ = 0,20 am Rand", K.leerraum_rand()[1], 0.20),
        ("F:ganzes", "verlangten diese Leerräume das 1,9-fache des Raums", K.zwei_schwellen()[0] / K.leerraum_rand()[1], 1.9),
        ("F:ganzes", "verlangten diese Leerräume das 1,9-fache des Raums", K.zwei_schwellen(numerisch=True)[0] * K.leerraum_rand()[2]**3, 1.9),
        ("F:ganzes", "Knoten und Leerraum tragen je 0,080 der Masse", K.netz_achsen(numerisch=True)["knoten"], 0.080),
        ("F:ganzes", "Knoten und Leerraum tragen je 0,080 der Masse", K.netz_achsen(numerisch=True)["leer"], 0.080),
        ("F:ganzes", "Snider, 22,5 m)", K.POUND["h"], 22.5),
        ("F:ganzes", "| 2,46·10⁻¹⁵ |", 3 * K.faltung(K.M_earth) * K.POUND["h"] / K.R_ERDE**2 / 1e-15, 2.46),
        ("F:ganzes", "gemessen auf 1 % genau", K.POUND["genau"] * 100, 1),
        ("F:ganzes", "Masse aus der Bahn 4,30·10⁶ M☉", _s2()[0] / K.M_sun / 1e6, 4.30),
        ("F:ganzes", "| 12,2′ je Umlauf |", math.degrees(K.perihel_je_umlauf(_s2()[0], _s2()[1], K.S2["e"])) * 60, 12.2),
        ("F:ganzes", "1,10 ± 0,19 davon", K.S2["f"][0], 1.10),
        ("F:ganzes", "1,10 ± 0,19 davon (GRAVITY)", K.S2["f"][1], 0.19),
        ("F:ganzes", "| 53,4 µas |", math.degrees(K.schatten(_s2()[0], K.S2["R0"] * K.pc)) * 3600e6, 53.4),
        ("F:ganzes", "51,8 ± 2,3 µas", K.SGRA_RING[0], 51.8),
        ("F:ganzes", "± 2,3 µas (Event", K.SGRA_RING[1], 2.3),
        ("F:ganzes", "Ṗ_b = −1,248·10⁻¹²", K.doppelpulsar(K.J0737["m1"] * K.M_sun, K.J0737["m2"] * K.M_sun, K.J0737["Pb"], K.J0737["e"]) / 1e-12, -1.248),
        ("F:ganzes", "| 0,75 ± 0,13 ± 0,16 g |", K.ALPHA_G[0], 0.75),
        ("F:ganzes", "0,75 ± 0,13 ± 0,16 g", K.ALPHA_G[1], 0.13),
        ("F:ganzes", "± 0,13 ± 0,16 g |", K.ALPHA_G[2], 0.16),
        ("F:ganzes", "bestätigt auf 1,3·10⁻⁴", K.J0737["genau"] / 1e-4, 1.3),
        ("03", "615 Bins, τ", float(len(D.planck_reste(D.PLANCK_START)["spektren"])), 615),
        ("03", "τ = 0,054 ± 0,007 und Eichung", K.TAU_PLANCK[0], 0.054),
        ("03", "0,054 ± 0,007 und Eichung", K.TAU_PLANCK[1], 0.007),
        ("03", "allein liegt das Modell Δχ² = 1,0 über", _cmb()["spektren"][1]["chi2"]["gesamt"] - _cmb()["spektren"][0]["chi2"]["gesamt"], 1.0),
        ("03", "an allen Daten 0,01.", _cmb()["alle"][1]["chi2"]["gesamt"] - _cmb()["alle"][0]["chi2"]["gesamt"], 0.01),
        ("03", "an allen Daten 0,01.", _cmb()["nach"] - _cmb()["alle"][0]["chi2"]["gesamt"], 0.01),
        ("03", "DESI mit χ² = 18,6", _cmb()["spektren"][1]["chi2"]["bao"], 18.6),
        ("03", "18,6 (13 Werte)", float(len(D.desi()[0])), 13),
        ("03", "ΛCDM mit 31,4", _cmb()["spektren"][0]["chi2"]["bao"], 31.4),
        ("03", "gewinnt das Modell Δχ² = 9,6", _cmb()["alle"][1]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 9.6),
        ("03", "w₀ = −0,84 ± 0,05", _cmb()["zustand"]["punkt"]["w0"], -0.84),
        ("03", "w₀ = −0,84 ± 0,05", _cmb()["zustand"]["fehler"]["w0"], 0.05),
        ("03", "w_a = −0,56 ± 0,18", _cmb()["zustand"]["punkt"]["wa"], -0.56),
        ("03", "w_a = −0,56 ± 0,18", _cmb()["zustand"]["fehler"]["wa"], 0.18),
        ("03", "bei z = 0,40 ± 0,07", _cmb()["z_max"][0], 0.40),
        ("03", "bei z = 0,40 ± 0,07", _cmb()["z_max"][1], 0.07),
        ("03", "die Schwelle liest 1,070 ± 0,013", _cmb()["zustand"]["punkt"]["a0"], 1.070),
        ("03", "die Schwelle liest 1,070 ± 0,013", _cmb()["zustand"]["fehler"]["a0"], 0.013),
        ("03", "der Standardwert nur Δχ² = 0,8 daneben", _cmb()["zustand_fest"]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 0.8),
        ("F:lesungen", "nur auf 1,3 % genau", 100 * _cmb()["zustand"]["fehler"]["a0"] / _cmb()["zustand"]["punkt"]["a0"], 1.3),
        ("03", "die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["punkt"]["a0"], 1.085),
        ("03", "die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["fehler"]["a0"], 0.003),
        ("03", "Das Maximum liegt dann bei z = 0,32", _cmb()["z_saett"], 0.32),
        ("03", "wie ΛCDM und liegt Δχ² = 8,7 darunter", _cmb()["alle"][0]["chi2"]["gesamt"] - _cmb()["saett"]["chi2"]["gesamt"], 8.7),
        ("03", "der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils", 100 * _cmb()["zustand_form"][0][0], 5.4),
        ("03", "der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils", 100 * _cmb()["zustand_form"][0][1], 0.7),
        ("F:ganzes", "F = 0,97 ± 0,12 %", 100 * _cmb()["zustand_form"][1][0], 0.97),
        ("F:ganzes", "F = 0,97 ± 0,12 %", 100 * _cmb()["zustand_form"][1][1], 0.12),
        ("F:ganzes", "am Maximum (z = 0,30) F = 1,68 ± 0,36 %", 100 * _cmb()["zustand_form"][2][0], 1.68),
        ("F:ganzes", "am Maximum (z = 0,30) F = 1,68 ± 0,36 %", 100 * _cmb()["zustand_form"][2][1], 0.36),
        ("F:ganzes", "am Maximum (z = 0,30) F = 1,68", _cmb()["zustand_form"][3][0], 0.30),
        ("F:ganzes", "bei z = 2 läge F bei −3,0 %", 100 * _gas_form()["F_z2"], -3.0),
        ("F:ganzes", "lesen alle Präzisionsdaten χ² = 2023,9", _gas_form()["chi2"], 2023.9),
        ("F:ganzes", "3,9 unter ΛCDM bei gleich vielen Zahlen, 4,2 über dem Gang", _gas_form()["lcdm"] - _gas_form()["chi2"], 3.9),
        ("F:ganzes", "3,9 unter ΛCDM bei gleich vielen Zahlen, 4,2 über dem Gang", _gas_form()["chi2"] - _gas_form()["saett"], 4.2),
        ("10", "die Faltung das 8,1·10⁻³⁷-fache der Windung", K.kopplung_grav(K.m_p) / K.alpha * 1e37, 8.1),
    ("10", "numerisch auf 10⁻⁵ gleich", float(abs(K.rueckkehr_diffusion(1.0) / math.exp(-1.0) - 1) < 1e-5), 1),
    ("10", "auf das 1,582-fache", K.verstaerkung(K.a0), 1.582),
    ("10", "Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache",
     K.rueckkehr_lesart(0.5, 1), 1),
    ("10", "Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache",
     K.rueckkehr_lesart(0.9, 1), 9),
    ("10", "Bei der halben höchsten Rate ist c = K_m, bei 90 % das Neunfache, bei 99 % das 99-fache",
     K.rueckkehr_lesart(0.99, 1), 99),
    ("10", "von 90 auf 99 % der Rate elffach", K.rueckkehr_lesart(0.99, 1) / K.rueckkehr_lesart(0.9, 1), 11),
    ("10", "Aus dem Gleichgewicht numerisch gelöst gleich auf 10⁻¹⁰", float(max(
        abs(K.rueckkehr_lesart(q, 1, umkehr=True) / K.rueckkehr_lesart(q, 1) - 1) for q in (0.5, 0.9, 0.99)) < 1e-10), 1),
    ("10", "bei s = 0,1 um 0,508 und 0,492 Gaben", 1 / (1 - math.exp(-0.1)) - 10, 0.508),
    ("10", "bei s = 0,1 um 0,508 und 0,492 Gaben", 10 - (1 / (1 - math.exp(-0.1)) - 1), 0.492),
    ("10", "Gabe für Gabe summiert gleich auf 10⁻¹²", float(max(abs((1 + K.rueckkehr_lesart(math.exp(-s), 1, True))
        / (0.5 + 0.5 / math.tanh(s / 2)) - 1) for s in (0.1, 1.0)) < 1e-12), 1),
    ("10", "in einer gewürfelten Folge von Binden und Lösen gezählt gleich auf 0,2 %", float(max(
        abs(K.wechsel_rate(a, b, True) / K.wechsel_rate(a, b) - 1) for a, b in ((1.0, 1.0), (3.0, 0.5), (0.2, 5.0))) < 2e-3), 1),
    ("10", "W = S/(S + k_Z) = τ·[L]/(K_d + (1 + τ)·[L])", float(max(abs(
        K.wechsel_rate(L, kab) / (K.wechsel_rate(L, kab) + 1.0) / ((kab / 1.0) * L / (kab + (1 + kab / 1.0) * L)) - 1)
        for L in (0.1, 1.0, 10.0) for kab in (0.1, 1.0, 10.0)) < 1e-12), 1),
    ("10", "in einer gewürfelten Folge mit Grundtakt gezählt gleich auf 0,2 %", float(max(
        abs(K.wechsel_grundtakt(t, s, k, True) / K.wechsel_grundtakt(t, s, k) - 1)
        for t, s, k in ((0.5, 1.0, 3.0), (0.8, 1.0, 0.2), (0.3, 2.0, 2.0))) < 2e-3), 1),
    ("10", "Wirkstoffs mit Δω = +0,3·K auf das 0,39-fache", K.allosterie(0.3, 0.1), 0.39),
    ("10", "Δω = −0,3·K, auf das 2,79-fache", K.allosterie(-0.3, 0.1), 2.79),
    ("10", "Ein genau passender Wirkstoff verliert immer, auf das 0,30-fache", K.allosterie(0.0, 0.1), 0.30),
    ("10", "aus dem abgetasteten Potential gleich auf 10⁻⁹", float(max(
        abs(K.allosterie(d, 0.1, numerisch=True) / K.allosterie(d, 0.1) - 1) for d in (0.3, -0.3, 0.0)) < 1e-9), 1),
    ("10", "ΔV ≈ 2K − π·|Δω|", float(abs((2 - K.einrasten_schwelle(1e-6, 1.0)) / 1e-6 - math.pi) < 1e-4), 1),
    ("10", "bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11", K.hill_zahl(4 / 2.0), 2),
    ("10", "bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11", K.hill_zahl(4 / 1.25), 5),
    ("10", "bei T = 2·T_K ist n_H = 2, bei 1,25·T_K schon 5, bei 1,1·T_K 11", K.hill_zahl(4 / 1.1), 11),
    ("10", "aus der selbstkonsistenten Besetzung abgeleitet gleich auf 10⁻⁹", float(max(abs(
        (math.log(K.besetzung(-2 / t + 1e-6, 4 / t) / (1 - K.besetzung(-2 / t + 1e-6, 4 / t)))
         - math.log(K.besetzung(-2 / t - 1e-6, 4 / t) / (1 - K.besetzung(-2 / t - 1e-6, 4 / t)))) / 2e-6
        / K.hill_zahl(4 / t) - 1) for t in (2.0, 1.25, 1.1)) < 1e-9), 1),
    ("10", "zwischen 0,99·T_K und 0,999·T_K mit dem Exponenten 1,505",
     math.log(K.besetzung_hysterese(4 / 0.99) / K.besetzung_hysterese(4 / 0.999)) / math.log(10), 1.505),
    ("10", "über die Form der eingerasteten Windungen gleich",
     float(abs(4 / 0.99 * K.hysteresefenster(0.99) / K.besetzung_hysterese(4 / 0.99) - 1) < 1e-12), 1),
    ("10", "bei T = 2·T_N zwei Drittel", K.hill_zahl(-4 / 2.0), 0.667),
    ("10", "Gruppe für Gruppe im selben Zustand wie nach −7, −1, zu 96 % besetzt",
     float((K.gewebe_gedaechtnis([-7, -2, -5, -3, -4, -1]) == K.gewebe_gedaechtnis([-7, -1])).all()), 1),
    ("10", "Gruppe für Gruppe im selben Zustand wie nach −7, −1, zu 96 % besetzt", 100 * K.gewebe_gedaechtnis([-7, -1]).mean(), 96),
    ("10", "wieder im Zustand nach −7, −1, −5, zu 13 % besetzt",
     float((K.gewebe_gedaechtnis([-7, -1, -5, -3, -5]) == K.gewebe_gedaechtnis([-7, -1, -5])).all()), 1),
    ("10", "wieder im Zustand nach −7, −1, −5, zu 13 % besetzt", 100 * K.gewebe_gedaechtnis([-7, -1, -5]).mean(), 13),
    ("10", "dazwischen bei −3 waren es 40 %", 100 * K.gewebe_gedaechtnis([-7, -1, -5, -3]).mean(), 40),
    ("10", "sobald der Zustand an jedem Extremwert einschwingt", float(all(
        (K.gewebe_gedaechtnis(f, True) == K.gewebe_gedaechtnis(f)).all()
        for f in ([-7, -2, -5, -3, -4, -1], [-7, -1], [-7, -1, -5], [-7, -1, -5, -3, -5], [-7, -1, -5, -3]))), 1),
    ("10", "aus finiten Differenzen gleich auf 10⁻⁸", float(max(
        abs(K.rueckkehr_membran(x, True) / K.rueckkehr_membran(x) - 1) for x in (0.5, 1.0, 10.0)) < 1e-8), 1),
    ("10", "bei x = 1 halb so groß", 1 - K.rueckkehr_membran(1.0), 0.5),
    ("10", "bei x = 10 elfmal so lang", 1 / (1 - K.rueckkehr_membran(10.0)), 11),
    ("10", "im Kanal zehnfach, 10,5, 100,5, 1000,5", 1 / (1 - K.rueckkehr_dimension(1, 2.0, 0.1)), 10.5),
    ("10", "im Kanal zehnfach, 10,5, 100,5, 1000,5", 1 / (1 - K.rueckkehr_dimension(1, 2.0, 0.01)), 100.5),
    ("10", "im Kanal zehnfach, 10,5, 100,5, 1000,5", 1 / (1 - K.rueckkehr_dimension(1, 2.0, 0.001)), 1000.5),
    ("10", "auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1", 1 / (1 - K.rueckkehr_dimension(2, 2.0, 0.1)), 3.6),
    ("10", "auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1", 1 / (1 - K.rueckkehr_dimension(2, 2.0, 0.01)), 6.8),
    ("10", "auf der Fläche nur um einen festen Schritt, 3,6, 6,8, 10,1", 1 / (1 - K.rueckkehr_dimension(2, 2.0, 0.001)), 10.1),
    ("10", "in freier Lösung kaum, 1,83 und 1,98", 1 / (1 - K.rueckkehr_dimension(3, 2.0, 0.1)), 1.83),
    ("10", "in freier Lösung kaum, 1,83 und 1,98", 1 / (1 - K.rueckkehr_dimension(3, 2.0, 0.01)), 1.98),
    ("10", "Aus finiten Differenzen gleich auf 10⁻⁴", float(max(abs(K.rueckkehr_dimension(d, 2.0, k, True)
        / K.rueckkehr_dimension(d, 2.0, k) - 1) for d in (1, 2, 3) for k in (0.1, 0.01)) < 1e-4), 1),
    ("10", "bei s = 1, also E = 1 − 1/e = 0,63", 1 - math.exp(-1), 0.63),
    ("10", "hält der Körper den Wirkstoff ν = 1,582 Durchgänge", 1 + K.rueckkehr_lesart(math.exp(-1), 1, True), 1.582),
    ("10", "dieselbe Verstärkung wie die Regel an r*", float(abs(K.verstaerkung(K.a0) * (1 - math.exp(-1)) - 1) < 1e-12), 1),
    ("10", "beide mit der Rate (1 − q)/τ, numerisch integriert gleich auf 10⁻¹¹", float(max(
        abs(x / y - 1) for q in (-1.0, -3.0, 0.5, 0.9) for x, y in zip(K.regelkreis(q), K.regelkreis(q, numerisch=True)))
        < 1e-11), 1),
    ("10", "bei |q| = 1 je die Hälfte", K.regelkreis(-1.0)[0], 0.5),
    ("10", "bei |q| = 1 je die Hälfte", -K.regelkreis(-1.0)[1], 0.5),
    ("10", "bei |q| = 3 bleibt ein Viertel", K.regelkreis(-3.0)[0], 0.25),
    ("10", "Bei q = 0,9 wächst die Wirkung auf das Zehnfache", K.regelkreis(0.9)[0], 10),
    ("10", "läuft nach dem Absetzen mit dem Neunfachen nach", K.regelkreis(0.9)[1], 9),
    ("10", "klingt mit einem Zehntel der Rate ab", K.regelkreis(0.9)[2], 0.1),
    ("10", "für m = 1 genau K/π = 0,064 breit", K.resonanz_fenster(0, 1, 0.2), 0.064),
    ("10", "für m = 2 noch 0,0032, für m = 3 noch 0,00028", K.resonanz_fenster(1, 2, 0.2), 0.0032),
    ("10", "für m = 2 noch 0,0032, für m = 3 noch 0,00028", K.resonanz_fenster(1, 3, 0.2), 0.00028),
    ("10", "mit den Exponenten 1,000, 1,996 und 2,994",
     math.log(K.resonanz_fenster(0, 1, 0.2) / K.resonanz_fenster(0, 1, 0.1), 2), 1.000),
    ("10", "mit den Exponenten 1,000, 1,996 und 2,994",
     math.log(K.resonanz_fenster(1, 2, 0.2) / K.resonanz_fenster(1, 2, 0.1), 2), 1.996),
    ("10", "mit den Exponenten 1,000, 1,996 und 2,994",
     math.log(K.resonanz_fenster(1, 3, 0.2) / K.resonanz_fenster(1, 3, 0.1), 2), 2.994),
    ("10", "Aus der Umlaufzahl in 10⁶ Schritten gleich auf 10⁻⁵", float(max(
        abs(K.resonanz_fenster(p, q, k, True) / K.resonanz_fenster(p, q, k) - 1)
        for p, q, k in ((0, 1, 0.2), (1, 2, 0.2), (1, 3, 0.3))) < 1e-5), 1),
    ("10", "Bei E = 0,1 sind es 10 Durchgänge gegen 9,49", 1 + K.rueckkehr_lesart(0.9, 1, True), 10),
    ("10", "Bei E = 0,1 sind es 10 Durchgänge gegen 9,49", -1 / math.log(0.9), 9.49),
    ("10", "bei E = 0,01 ist der Unterschied 0,501", 1 / 0.01 + 1 / math.log(0.99), 0.501),
    ("10", "Gliedweise summiert gleich auf 10⁻¹²", float(max(
        max(abs((1 + K.rueckkehr_lesart(1 - E, 1, True)) * E - 1),
            abs((0.5 + 0.5 / math.tanh(-math.log(1 - E) / 2)) * E - 1)) for E in (0.5, 0.1, 0.01)) < 1e-12), 1),
    ("10", "Gliedweise summiert sind beide gleich auf 10⁻¹⁵", float(max(
        abs(K.rueckkehr_lesart(q, s, True) / K.rueckkehr_lesart(q, s) - 1)
        for q, s in ((0.3, 1), (0.9, 1), (0.3, 2), (0.9, 2), (2.5, 2))) < 1e-15), 1),
    ("10", "bei n_H = 2 mit der halben Rate, numerisch integriert gleich auf 10⁻⁵", K.besetzung_rueckkehr(-1.0, 2.0), 0.5),
    ("10", "bei n_H = 2 mit der halben Rate, numerisch integriert gleich auf 10⁻⁵", float(max(
        abs(K.besetzung_rueckkehr(-J / 2, J, True) / K.besetzung_rueckkehr(-J / 2, J) - 1) for J in (2.0, 3.0)) < 1e-5), 1),
    ("10", "bei J = 4,01 das 0,333-fache von (J − 4)^{3/2}", K.besetzung_hysterese(4.01) / 0.01**1.5, 0.333),
    ("10", "Durch Zählen der Lösungen gleich auf 10⁻⁸", float(max(
        abs(K.besetzung_hysterese(J, True) / K.besetzung_hysterese(J) - 1) for J in (4.01, 5.0, 8.0)) < 1e-8), 1),
    ("10", "Vor jeder Falte fällt die Rate der Rückkehr mit der Wurzel des Abstands", float(abs(math.log(
        K.besetzung_rueckkehr(-2.5 + K.besetzung_hysterese(5.0) / 2 - 1e-4, 5.0)
        / K.besetzung_rueckkehr(-2.5 + K.besetzung_hysterese(5.0) / 2 - 1e-6, 5.0)) / math.log(100) - 0.5) < 0.01), 1),
    ("10", "und ihr Verhältnis K_d bleibt", K.bindung_feld(0.9)[0] / K.bindung_feld(0.1)[0], 1),
    ("10", "geben stationäre Verteilung und mittlere Übergangszeiten dasselbe, numerisch gleich auf 10⁻¹²", float(max(
        abs(x / y - 1) for q in (0.1, 0.5, 0.9) for x, y in zip(K.bindung_feld(q), K.bindung_feld(q, numerisch=True))) < 1e-12), 1),
    ("10", "Bei c_eff = 10·K_d hält der Wirkstoff 11-mal so lange", 1 / (1 - K.besetzung(math.log(10.0))), 11),
    ("10", "(bei s = 0,2 das 1,90-fache)", K.verstaerkung(0.01 * K.a0) / K.verstaerkung(0.04 * K.a0), 1.90),
    ("10", "J = 2·kT gibt n_H = 2", K.hill_zahl(2.0), 2),
    ("10", "numerisch abgeleitet gleich auf 10⁻⁶", float(max(abs((math.log(K.besetzung(-J / 2 + 1e-6, J) / (1 - K.besetzung(-J / 2 + 1e-6, J)))
     - math.log(K.besetzung(-J / 2 - 1e-6, J) / (1 - K.besetzung(-J / 2 - 1e-6, J)))) / 2e-6 / K.hill_zahl(J) - 1) for J in (1.0, 2.0, 3.0)) < 1e-6), 1),
    ("10", "bei W = 10·k sind es 7", math.floor(K.zyklen_im_band(10.0)), 7),
    ("F:andocken", "An 68 Kinasehemmern gegen 442 Kinasen", D.davis().shape[0], 68),
    ("F:andocken", "An 68 Kinasehemmern gegen 442 Kinasen", D.davis().shape[1], 442),
    ("F:andocken", "haben die Rangkorrelation 0,43, mit 1 µM als Schwelle 0,41", D.davis_fenster()[2], 0.43),
    ("F:andocken", "haben die Rangkorrelation 0,43, mit 1 µM als Schwelle 0,41", D.davis_form(schwelle=6.0)[0], 0.41),
    ("F:andocken", "ihr oberes Zehntel beginnt bei 0,69 der stärksten statt bei 0,85", D.davis_form()[0], 0.69),
    ("F:andocken", "im ähnlichsten Viertel 38 %, im fernsten 24 %", 100 * D.davis_sequenz()[0], 38),
    ("F:andocken", "im ähnlichsten Viertel 38 %, im fernsten 24 %", 100 * D.davis_sequenz()[3], 24),
    ("F:andocken", "mit Gruppen aus zwei Aminosäuren 35 % und 25 %", 100 * D.davis_sequenz(k=2)[0], 35),
    ("F:andocken", "mit Gruppen aus zwei Aminosäuren 35 % und 25 %", 100 * D.davis_sequenz(k=2)[3], 25),
    ("F:andocken", "ihr oberes Zehntel beginnt bei 0,69 der stärksten statt bei 0,85", D.davis_form()[1], 0.85),
    ("10", "θ(h) + θ(−h) = 1", K.besetzung(1.3) + K.besetzung(-1.3), 1.0),
    ("10", "Bei Δω = 1,01·K dauert ein Schlupf 44/K, bei Δω = 2·K nur 3,6/K", K.einrasten_zeit(1.01, 1.0), 44),
    ("10", "Bei Δω = 1,01·K dauert ein Schlupf 44/K, bei Δω = 2·K nur 3,6/K", K.einrasten_zeit(2.0, 1.0), 3.6),
    ("10", "bei voller Passung 2K, bei Δω = K/2 nur noch 34 % davon", K.einrasten_schwelle(0.0, 1.0), 2),
    ("10", "bei voller Passung 2K, bei Δω = K/2 nur noch 34 % davon", 100 * K.einrasten_schwelle(0.5, 1.0) / 2, 34),
    ("10", "Bei [L] = K_d braucht er das Doppelte seiner eigenen Bindungskonstante", K.hemmer_haelfte(1.0), 2),
    ("10", "bei [L] = 3·K_d das Vierfache", K.hemmer_haelfte(3.0), 4),
    ("10", "Aus der stationären Verteilung der Übergänge numerisch gleich auf 10⁻¹²",
     float(max(abs(K.hemmer_haelfte(x) / (1 + x) - 1) for x in (0.5, 1.0, 3.0)) < 1e-12), 1),
    ("10", "drei Punkte legt immer eine Drehung auf ihr Spiegelbild, es bleibt kein Abstand",
     float(K.spiegel_deckung([(1, 1, 1), (1, -1, -1), (-1, 1, -1)]) < 1e-12), 1),
    ("10", "ein mittlerer Abstand von 0,71 der Kantenlänge, über alle Drehungen abgetastet ebenso",
     K.spiegel_deckung([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]) / math.sqrt(8), 0.71),
    ("10", "ein mittlerer Abstand von 0,71 der Kantenlänge, über alle Drehungen abgetastet ebenso",
     K.spiegel_deckung([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], gitter=36) / math.sqrt(8), 0.71),
    ("10", "Geschlossene Formen und numerische Integration stimmen auf 10⁻⁵", float(max(
        [abs(K.einrasten_zeit(d, 1.0) / K.einrasten_numerisch(d, 1.0) - 1) for d in (0.5, 0.9, 1.01, 2.0)]
        + [abs(K.einrasten_schwelle(d, 1.0) / K.einrasten_numerisch(d, 1.0, "schwelle") - 1) for d in (0.0, 0.5, 0.9)]) < 1e-5), 1),
    ("F:andocken", "ist am größten bei N_e = (3c/a)^{3/2}", K.bindung_je_atom_groesse(1.0, 2.0, numerisch=True) / 6**1.5, 1),
    ("F:andocken", "gleich bei 68 von 68", float(sum(D.davis_atome() == D.davis_atome("pubchem"))), 68),
    ("F:andocken", "die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37", D.davis_groesse()[0], 0.33),
    ("F:andocken", "die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37", D.davis_groesse()[1], 0.006),
    ("F:andocken", "die Rangkorrelation 0,33 (p = 0,006), unter 1 µM 0,37", D.davis_groesse(6.0)[0], 0.37),
    ("F:andocken", "Größe und stärkste Bindung 0,25", D.davis_groesse()[2], 0.25),
    ("F:andocken", "wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8}", D.davis_groesse(5.0, 4000)[4], 1.9),
    ("F:andocken", "wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8}", D.davis_groesse(5.0, 4000)[5], 0.7),
    ("F:andocken", "wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8}", D.davis_groesse(6.0, 4000)[4], 2.4),
    ("F:andocken", "wächst wie N^{1,9 ± 0,7}, unter 1 µM wie N^{2,4 ± 0,8}", D.davis_groesse(6.0, 4000)[5], 0.8),
    ("F:andocken", "An den K_i und K_d von ChEMBL (5985 Spiegelpaare", D.chembl_spiegel()["paare"], 5985),
    ("F:andocken", "über vertauschtes @ und @@ in den SMILES zu 95 % bestätigt, 10 246 Paare",
     100 * D.chembl_spiegel()["smiles_gleich"], 95),
    ("F:andocken", "über vertauschtes @ und @@ in den SMILES zu 95 % bestätigt, 10 246 Paare", D.chembl_spiegel()["gemessen"], 10246),
    ("F:andocken", "im Median um 0,51 in pK", D.chembl_spiegel()["median"], 0.51),
    ("F:andocken", "die Rangkorrelation mit der Zahl der Atome ist −0,07", D.chembl_spiegel()["rho"], -0.07),
    ("F:andocken", "bis 15 Atome 14 % der Paare, bei 31 bis 40 Atomen 22 %", 100 * D.chembl_spiegel()["bereiche"][0][2], 14),
    ("F:andocken", "bis 15 Atome 14 % der Paare, bei 31 bis 40 Atomen 22 %", 100 * D.chembl_spiegel()["bereiche"][4][2], 22),
    ("F:andocken", "An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare", D.wirkweise_treffer()[0], 643),
    ("F:andocken", "An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare", D.wirkweise_treffer()[1], 30),
    ("F:andocken", "An ChEMBL (643 Messungen von k_off in s⁻¹ an 30 Zielen) treffen nur 15 Paare", D.wirkweise_treffer()[2], 15),
    ("F:andocken", "alle Antagonisten", float(set(D.wirkweise_treffer()[3]) == {"anta"}), 1),
    ("F:andocken", "führt nur 22 k_off mit Ziel und Wirkstoff, keines trifft eine Wirkweise", len(D.bindingdb_kinetik()), 22),
    ("F:andocken", "führt nur 22 k_off mit Ziel und Wirkstoff, keines trifft eine Wirkweise",
     float(sum((u, k) in D.gtopdb_wirkweise() for u, k, _, _, _ in D.bindingdb_kinetik())), 0),
    ("F:andocken", "235 599 Moleküle", len(D.chembl_molekuele()[0]), 235599),
    ("F:andocken", "an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1",
     D.chembl_atome_vergleich()[0], 2436),
    ("F:andocken", "an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1",
     100 * D.chembl_atome_vergleich()[1], 96),
    ("F:andocken", "an den 2436 Molekülen der oberen Grenze gegen die gezählten Atome bei 96 % gleich, bei 99,9 % auf ±1",
     100 * D.chembl_atome_vergleich()[2], 99.9),
    ("F:andocken", "zwischen 12 und 13 Atomen von 9,5 auf 10,4", float(D.chembl_huelle()[1][D.chembl_huelle()[0] == 12][0]), 9.5),
    ("F:andocken", "zwischen 12 und 13 Atomen von 9,5 auf 10,4", float(D.chembl_huelle()[1][D.chembl_huelle()[0] == 13][0]), 10.4),
    ("F:andocken", "bis 66 Atome zwischen 9,9 und 10,9", float(D.chembl_huelle()[0].max()), 66),
    ("F:andocken", "bis 66 Atome zwischen 9,9 und 10,9", float(D.chembl_huelle()[1][D.chembl_huelle()[0] >= 13].min()), 9.9),
    ("F:andocken", "bis 66 Atome zwischen 9,9 und 10,9", float(D.chembl_huelle()[1][D.chembl_huelle()[0] >= 13].max()), 10.9),
    ("F:ganzes", "Dann ist F heute 0,55 ± 0,27 %", 100 * _gas_form()["F0"][0], 0.55),
        ("F:ganzes", "Dann ist F heute 0,55 ± 0,27 %", 100 * _gas_form()["F0"][1], 0.27),
        ("F:ganzes", "schmale Spitze bei z = 0,61, χ² = 2017,9", _gas_form()["frei"][1], 0.61),
        ("F:ganzes", "schmale Spitze bei z = 0,61, χ² = 2017,9", _gas_form()["frei"][0], 2017.9),
        ("F:ganzes", "(F_s = 5 ± 5 %)", 100 * _gas_form()["frei"][2], 5),
        ("F:ganzes", "(F_s = 5 ± 5 %)", 100 * _gas_form()["frei"][3], 5),
        ("F:ganzes", "am Maximum (z = 0,30) 1,21 ± 0,33 %", 100 * _gas_form()["Fs"][0], 1.21),
        ("F:ganzes", "am Maximum (z = 0,30) 1,21 ± 0,33 %", 100 * _gas_form()["Fs"][1], 0.33),
        ("F:ganzes", "am Maximum (z = 0,30) 1,21 ± 0,33 %", _gas_form()["z_s"], 0.30),
        ("03", "lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006", _cmb()["q_lesen"]["punkt"]["Q"], 5.36),
        ("03", "lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006", _cmb()["q_lesen"]["fehler"]["Q"], 0.05),
        ("03", "lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006", _cmb()["q_lesen"]["punkt"]["a0"], 1.076),
        ("03", "lesen alle Präzisionsdaten Q = 5,36 ± 0,05 und a₀ = 1,076 ± 0,006", _cmb()["q_lesen"]["fehler"]["a0"], 0.006),
        ("03", "2π − 1 liegt Δχ² = 2,4 daneben.", _cmb()["saett"]["chi2"]["gesamt"] - _cmb()["q_lesen"]["chi2"]["gesamt"], 2.4),
        ("03", "ist L das 0,998-fache des Durchmessers", _cmb()["zelle"]["L_d"], 0.998),
        ("03", "den Quadrupol der Hintergrundstrahlung um 13 %", 100 * (1 - _cmb()["zelle"]["R2"]), 13),
        ("03", "hebt den Oktupol um 2 %", 100 * (_cmb()["zelle"]["R3"] - 1), 2),
        ("03", "Sachs-Wolfe: 16 % und 2 %", 100 * (1 - _cmb()["zelle"]["sw2"]), 16),
        ("03", "Sachs-Wolfe: 16 % und 2 %", 100 * (_cmb()["zelle"]["sw3"] - 1), 2),
        ("03", "ändert das χ² nur um 0,08", _cmb()["zelle"]["dchi2"], 0.08),
        ("03", "gegen 1434 mit 282", _sparc_standard()["chi_halo"], 1434),
        ("03", "1434 mit 282 freien Zahlen", float(_sparc_standard()["zahlen"]), 282),
        ("03", "χ² = 2918 gegen", _sparc_standard()["chi_regel"], 2918),
        ("03", "liegt die Regel 736 vorn", _sparc_standard()["bic_halo"] - _sparc_standard()["bic_regel"], 736),
        ("03", "vorn, 2926 gegen", _sparc_standard()["bic_regel"], 2926),
        ("03", "2926 gegen 3661", _sparc_standard()["bic_halo"], 3661),
        ("03", "höchstens 0,02 dex", _sparc_standard()["stapel"], 0.02),
        ("03", "An SPARC liegen 131 Punkte jenseits der Naht", float(_naht_sparc()[0]), 131),
        ("03", "ist +0,013 ± 0,006 dex", _naht_sparc()[1], 0.013),
        ("03", "+0,013 ± 0,006 dex", _naht_sparc()[2], 0.006),
        ("03", "lägen sie 0,094 dex darüber", _naht_sparc()[3], 0.094),
        ("03", "bei 31 der 141 Galaxien am Rand", float(_sparc_standard()["rand"]), 31),
        ("03", "um 0,58 dex: Der Halo", _sparc_standard()["masse_streuung"], 0.58),
        ("02", "(κ = 0,12, μ = 0,03: Schließzahlen", KE.NIVEAUS[1], 0.12),
        ("02", "μ = 0,03: Schließzahlen", KE.NIVEAUS[2], 0.03),
        ("02", "Schließzahlen 2, 8, 20, 28, 50, 82), gegen", float(KE.schliesszahlen_verwindung(0.12, 0.35, mu=0.03)[:6] == (2, 8, 20, 28, 50, 82)), 1.0),
        ("02", "A^{−1/3} = 38,0 MeV", KE.schalenabstand(KE.halterlaenge(KE.a_C(KE.TROPFEN))), 38.0),
        ("02", "liest unabhängig davon 39,3", KE.anpassung(KE.NIVEAUS + ("frei",))[0][7], 39.3),
        ("02", "39,3 (3 %)", (KE.anpassung(KE.NIVEAUS + ("frei",))[0][7] / KE.schalenabstand(KE.halterlaenge(KE.a_C(KE.TROPFEN))) - 1) * 100, 3),
        ("02", "also λ_H = 1,23 fm", KE.halterlaenge(KE.a_C(KE.TROPFEN)) * math.sqrt(KE.schalenabstand(KE.halterlaenge(KE.a_C(KE.TROPFEN)))
                                                                                    / KE.anpassung(KE.NIVEAUS + ("frei",))[0][7]), 1.23),
        ("02", "streut die Bindung um 1,83 MeV", KE.anpassung(KE.NIVEAUS)[1], 1.83),
        ("02", "um 1,26 statt 1,93 MeV", KE.streuung_ab(KE.NIVEAUS, 0, True), 1.26),
        ("02", "1,26 statt 1,93 MeV", KE.streuung_ab(KE.SCHALE, 0, True), 1.93),
        ("F:ganzes", "(82 gerade-gerade Strahler", float(KE.alpha_tunnel()[2]), 82),
        ("F:ganzes", "streut um 0,52 Dekaden", KE.alpha_tunnel()[1], 0.52),
        ("F:ganzes", "drei Zahlen (0,51)", KE.alpha_lebensdauer()[1], 0.51),
        ("F:ganzes", "liegt fest 1,24 Dekaden zu schnell", -KE.alpha_tunnel()[0], 1.24),
        ("F:ganzes", "(N = 118 bis 126, 13 Strahler)", float(KE.alpha_tunnel((118, 126))[2]), 13),
        ("F:ganzes", "1,69 ± 0,07 Dekaden zu schnell", -KE.alpha_tunnel((118, 126))[0], 1.69),
        ("F:ganzes", "1,69 ± 0,07 Dekaden", KE.alpha_tunnel((118, 126))[1], 0.07),
        ("F:ganzes", "(N = 128 bis 136, 23 Strahler)", float(KE.alpha_tunnel((128, 136))[2]), 23),
        ("F:ganzes", "nur 0,94 ± 0,03", -KE.alpha_tunnel((128, 136))[0], 0.94),
        ("F:ganzes", "0,94 ± 0,03.", KE.alpha_tunnel((128, 136))[1], 0.03),
        ("F:ganzes", "schließt sich 5,5-mal leichter", 10**(KE.alpha_tunnel((128, 136))[0] - KE.alpha_tunnel((118, 126))[0]), 5.5),
        ("F:ganzes", "in einem von 17 Anläufen", 10**-KE.alpha_tunnel()[0], 17),
        ("08", "0,8756″", math.degrees(K.ablenkung(K.M_sun, R_SONNE)) * 3600 / 2, 0.8756),
        ("08", "42,98″", math.degrees(K.perihel_je_umlauf(K.M_sun, 0.387098 * K.AU, 0.205630)) * 3600 * 100 * K.yr
         / (2 * math.pi * math.sqrt((0.387098 * K.AU)**3 / K.GM_sun)), 42.98),
        ("08", "r* = 7480 AE", K.schwellenradius(K.M_sun) / K.AU, 7480),
        ("F:galaxien", "An RC100 (100 sternbildende Galaxien", float(len(D.rc100())), 100),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle()[0], 1.00),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle()[1], 0.94),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle()[2], 1.08),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle(0.0, 1.2)[0], 0.87),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle(1.2, 2.0)[0], 1.21),
        ("F:galaxien", "liest er a₀ = 1,00 (0,94 bis 1,08)·10⁻¹⁰ m/s², zwischen z = 0,6 und 1,2, 1,2 und 2, 2 und 2,5 getrennt 0,87, 1,21 und 1,01", D.rc100_schwelle(2.0, 3.0)[0], 1.01),
        ("F:galaxien", "Das g_N der Tabelle ist das Feld der Baryonen aus Masse und Radius: V_bar²·R_e/(GM_b) liegt im Median bei 0,60, die dünne Scheibe gibt bei R_e 0,63", _rc100_baryonen()[0], 0.60),
        ("F:galaxien", "Das g_N der Tabelle ist das Feld der Baryonen aus Masse und Radius: V_bar²·R_e/(GM_b) liegt im Median bei 0,60, die dünne Scheibe gibt bei R_e 0,63", _rc100_baryonen()[1], 0.63),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", float(D.rc100_schwelle(fehler="voll")[0]), 1.28),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", float(D.rc100_schwelle(fehler="voll")[1]), 1.17),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", float(D.rc100_schwelle(fehler="voll")[2]), 1.39),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", float(D.rc100_schwelle(fehler="voll")[4]), 110),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", float(D.rc100_schwelle(fehler="voll")[3]), 100),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", _rc100_einzeln()[0], 1.39),
        ("F:galaxien", "Mit den Fehlern von M_b und R_e fortgepflanzt liest er 1,28 (1,17 bis 1,39) bei χ² = 110 für 100 Galaxien, aus jeder Galaxie einzeln, a₀ = g_N/(ln q)², im Median 1,39, 59 von 100 über 1,06", _rc100_einzeln()[1], 59),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang()[0], 0.14),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang()[1], -0.08),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang()[2], 0.36),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang("voll")[0], -0.29),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang("voll")[1], -0.58),
        ("F:galaxien", "Als a₀ ∝ H(z)^m gelesen ist m = 0,14 (−0,08 bis 0,36), mit allen Fehlern −0,29 (−0,58 bis −0,01)", D.rc100_gang("voll")[2], -0.01),
        ("F:galaxien", "Wüchse a₀ mit H(z), läge χ² um 15 und 20 höher", D.rc100_schwelle(m=1.0)[4] - D.rc100_schwelle()[4], 15),
        ("F:galaxien", "Wüchse a₀ mit H(z), läge χ² um 15 und 20 höher", D.rc100_schwelle(m=1.0, fehler="voll")[4] - D.rc100_schwelle(fehler="voll")[4], 20),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[0], 26),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[1] / 1e-10, 0.96),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[2][0] / 1e-10, 0.81),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[2][1] / 1e-10, 1.11),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[3] / 1e-10, 1.12),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[4][0] / 1e-10, 1.01),
        ("F:galaxien", "In SPARC liest die Ursa-Major-Gruppe (26 Galaxien) a₀ = 0,96 (0,81 bis 1,11)·10⁻¹⁰ m/s², die übrigen Galaxien 1,12 (1,01 bis 1,23)", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[4][1] / 1e-10, 1.23),
        ("F:galaxien", "liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[5], -0.03),
        ("F:galaxien", "liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[6], 0.05),
        ("F:galaxien", "liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[7], -0.04),
        ("F:galaxien", "liegen ihre Punkte im Mittel −0,03 ± 0,05 dex neben der Regel, die der übrigen −0,04 ± 0,02", D.sparc_gruppe(K.kpc, D.sparc_eichung(K.kpc))[8], 0.02),
        ("07", "über die Verteilung integriert gleich auf 10⁻¹⁵", float(max(abs(K.erster_durchgang(x) - K.erster_durchgang(x, True)) for x in (0.01, 0.1, 1.0, 4.0, 25.0)) < 1e-15), 1),
        ("08", "r*(m_Pl) = √(ℓ_P·ℓ₀) = 0,117 mm", 1e3 * K.schwellenradius(K.m_Pl), 0.117),
        ("08", "r*(m_Pl) = √(ℓ_P·ℓ₀) = 0,117 mm", K.schwellenradius(K.m_Pl) / math.sqrt(K.l_P * K.c**2 / K.a0), 1),
        ("08", "ρ_s·c² = (π³/30)·ħc/r*(m_Pl)⁴", (math.pi**3 * K.a0**2 / (30 * K.G)) / (math.pi**3 / 30 * K.hbar * K.c / K.schwellenradius(K.m_Pl)**4), 1),
        ("08", "mit 1 + Q = 2π also ε_L·N_b = (3/4π)·m_Pl/m_p", K.leiterzahl() * 3 * K.m_Pl / (2 * (1 + K.q_je_radiant()) * K.leiterzahl() * K.m_p) * K.m_p / K.m_Pl * 4 * math.pi / 3, 1),
        ("08", "Die Zelle trägt N_b = 1,6·10⁸⁰ Baryonen, aus ρ_b·L³ gezählt ebenso", 3 * K.m_Pl / (2 * (1 + K.q_je_radiant()) * K.leiterzahl() * K.m_p) / 1e80, 1.6),
        ("08", "Die Zelle trägt N_b = 1,6·10⁸⁰ Baryonen, aus ρ_b·L³ gezählt ebenso", K.omega_b * K.rho_krit100 * (_zelle() * K.Gpc)**3 / K.m_p / 1e80, 1.6),
        ("03", "Θ ist das Verhältnis zweier Schleifen, der Wärmelänge des Horizonts 2πc/H₀ zur Schleife der inneren Richtung ℓ₀ = 2Λ", (2 * math.pi * K.c / K.H0) / (K.c**2 / K.a0) / K.theta_(), 1),
        ("03", "Gleich lang sind sie, wenn drei Vakuumteile Ω_DE = π²/15 = 0,658 tragen", 3 * K.vakuumteil(K.a0, 2 * math.pi * K.a0 / K.c), 0.658),
        ("03", "Gleich lang sind sie, wenn drei Vakuumteile Ω_DE = π²/15 = 0,658 tragen", math.pi**2 / 15, 0.658),
        ("01", "Beide hängen über das Ganze zusammen: δ·(R/r*)² = R/ℓ₀ mit δ = 3Δs/R", float(max(abs(3 * K.faltung(M) / R * (R / K.schwellenradius(M))**2 / (R * K.a0 / K.c**2) - 1)
                                 for M, R in ((K.M_sun, K.R_SONNE), (K.M_earth, K.R_ERDE), (1e11 * K.M_sun, 10 * K.kpc))) < 1e-12), 1),
        ("01", "Bei der Verwindung ist es die Spinzahl, j = k·c·v/v_K² mit k = I/(MR²)", (lambda om: (K.J_ERDE / (K.M_earth * K.c)) / (3 * K.faltung(K.M_earth)) / (K.J_ERDE / (K.M_earth * K.R_ERDE**2 * om) * K.c * om * K.R_ERDE / (K.G * K.M_earth / K.R_ERDE)))(7.292e-5), 1),
        ("08", "M_c teilt die Leiter im Drittel: Dort ist λ̄ = r* = ℓ_P^{2/3}·ℓ₀^{1/3}", (lambda Mc: K.hbar / (Mc * K.c) / (K.l_P**(2 / 3) * (K.c**2 / K.a0)**(1 / 3)))(K.schnittmassen()[1]), 1),
        ("08", "M_c teilt die Leiter im Drittel: Dort ist λ̄ = r* = ℓ_P^{2/3}·ℓ₀^{1/3}", (lambda Mc: K.schwellenradius(Mc) / (K.hbar / (Mc * K.c)))(K.schnittmassen()[1]), 1),
        ("F:ganzes", "Die Schließung muss das Sechseck umfassen, 1 + Q ≥ 2d, und 1 + Q = 2π ist der Kreis durch seine Ecken, r₀/Λ = d/π", float(all((math.sqrt(6 / (1 + Q)) <= 1) == (1 + Q >= 2 * 3) for Q in (4.0, 4.9, 5.0, 5.28, 6.0))), 1),
        ("F:ganzes", "Die Schließung muss das Sechseck umfassen, 1 + Q ≥ 2d, und 1 + Q = 2π ist der Kreis durch seine Ecken, r₀/Λ = d/π", float(abs(6 * K.faltung(1e50) / ((2 * math.pi) * K.faltung(1e50)) - 3 / math.pi) < 1e-15 and abs(2 * 6 * math.sin(math.pi / 6) - 6) < 1e-12), 1),
        ("F:ganzes", "Die Schwelle des Ganzen, r* = √(r₀·Λ), ist der Spiegel zwischen Sechseck und Kreis", float(abs(math.sqrt(6 * 2 * math.pi) / 6 - math.sqrt(math.pi / 3)) < 1e-15), 1),
        ("02", "An ihrer Schwelle ist δ = 3Δs/r* = r*/ℓ₀ = (v_flach/c)²", float(max(abs(3 * K.faltung(M) / K.schwellenradius(M) / (K.schwellenradius(M) * K.a0 / K.c**2) - 1) + abs((K.v_flach(M) / K.c)**2 / (3 * K.faltung(M) / K.schwellenradius(M)) - 1)
                                 for M in (1e8 * K.M_sun, 6e10 * K.M_sun, 1e14 * K.M_sun)) < 1e-12), 1),
        ("02", "an derselben Stelle, r*/Λ = √(6/(1 + Q)) = 0,977, und die Menge fällt heraus", math.sqrt(6 / (1 + K.q_je_radiant())), 0.977),
        ("02", "an derselben Stelle, r*/Λ = √(6/(1 + Q)) = 0,977, und die Menge fällt heraus", float(max(abs((lambda Nb: math.sqrt(K.G * Nb * K.m_p / (K.c**2 / (2 * (1 + K.q_je_radiant()) * Nb * K.faltung(K.m_p)))) / ((1 + K.q_je_radiant()) * Nb * K.faltung(K.m_p)))(Nb) - math.sqrt(6 / (1 + K.q_je_radiant()))) for Nb in (1e78, 1.6e80, 1e82)) < 1e-12), 1),
        ("03", "In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt", float(sum(g < x for _, g, x in D.haufenstoss_lage())), 4),
        ("03", "In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt", min(x for _, g, x in D.haufenstoss_lage()), 3.8),
        ("03", "In allen vier Teilhaufen liegt der Gipfel der Linse näher bei den Galaxien als beim Gas und 3,8 bis 13,8 Fehler vom Gasgipfel entfernt", max(x for _, g, x in D.haufenstoss_lage()), 13.8),
        ("03", "In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2", D.haufenstoss()["bradac2006"]["haupt_galaxie"][4], 16.3),
        ("03", "In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2", D.haufenstoss()["bradac2006"]["neben_galaxie"][4], 4.3),
        ("03", "In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2", K.q_je_radiant() * D.haufenstoss()["bradac2006"]["haupt_gas"][2], 10.6),
        ("03", "In 1E0657−56 tragen gleich große Bereiche um die Galaxien 16,3 und 4,3·10¹³ M☉, das Q-fache des Gases an den Gasgipfeln ist 10,6 und 2,2", K.q_je_radiant() * D.haufenstoss()["bradac2006"]["neben_gas"][2], 2.2),
        ("03", "Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05", D.haufenstoss()["clowe"]["haupt_gas"][4], 0.05),
        ("03", "Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05", D.haufenstoss()["clowe"]["neben_gas"][4], 0.02),
        ("03", "Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05", D.haufenstoss()["clowe"]["haupt_galaxie"][4], 0.36),
        ("03", "Am Gas bleibt nach Abzug der Galaxien kaum Linse, κ = 0,05 ± 0,06 und 0,02 ± 0,06, bei den Galaxien 0,36 ± 0,06 und 0,20 ± 0,05", D.haufenstoss()["clowe"]["neben_galaxie"][4], 0.20),
        ("03", "Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7", D.haufenstoss()["bradac2006"]["feld"][4] / D.haufenstoss()["bradac2006"]["feld"][2], 7.0),
        ("03", "Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7", (lambda f: f[4] / f[2] * math.hypot(f[5] / f[4], f[3] / f[2]))(D.haufenstoss()["bradac2006"]["feld"]), 0.9),
        ("03", "Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7", (1 + K.q_je_radiant()) / (2 * math.pi), 1),
        ("03", "Über das ganze Feld trägt die Linse das 7,0 ± 0,9-fache des Gases: das Maß 1 + Q = 2π, mit den Sternen 6,7", (1 + K.q_je_radiant()) * (1 + sum(v[2] for v in D.haufenstoss()["clowe"].values())
                                                 / sum(v[0] for v in D.haufenstoss()["clowe"].values())), 6.7),
        ("06", "gilt Σ1/⟨T⟩ = 1", float(abs(K.kac_wiederkehr(numerisch=True)[0] - 1) < 1e-3 and abs(K.kac_wiederkehr(0.2, 0.05, True)[0] - 1) < 1e-3), 1),
        ("07", "in drei kehrte die Irrfahrt nur zu 0,341 zurück (Pólya)", K.irrfahrt_rueckkehr(3), 0.341),
        ("07", "in drei kehrte die Irrfahrt nur zu 0,341 zurück (Pólya)", float(abs(K.irrfahrt_rueckkehr(3) - K.irrfahrt_rueckkehr(3, False)) < 1e-12), 1),
        ("07", "die Summe bliebe bei 1/(1 − 0,341) = 1,516", 1 / (1 - K.irrfahrt_rueckkehr(3, False)), 1.516),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", float(len(D.sparc_auswahl())), 141),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[0], 2.15),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[2][0], 1.56),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[2][1], 2.65),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[1], 0.75),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[3][0], 0.64),
        ("03", "An SPARC (141 Galaxien) liegt der Gipfel bei 2,15·a₀ (16 bis 84 % über Ziehen der Galaxien: 1,56 bis 2,65) mit 0,75·a₀ (0,64 bis 1,05)", _gipfel()[3][1], 1.05),
        ("03", "Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87", math.sqrt(32), 5.7),
        ("03", "Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87", _gipfel()[4], 0.57),
        ("03", "Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87", (K.verstaerkung(math.sqrt(32) * K.a0) - 1) * math.sqrt(32), 0.58),
        ("03", "Bei 5,7·a₀ liest SPARC 0,57·a₀, die Bose-Form gibt 0,58, die einfache Form, die keinen Gipfel hat, 0,87", math.sqrt(32) * (0.5 + math.sqrt(0.25 + 1 / math.sqrt(32)) - 1), 0.87),
        ("03", "Die Standard-Form läge mit ihrem Gipfel bei 0,49·a₀", _gipfel_standard()[0], 0.49),
        ("F:ganzes", "Im Katalog α.100 von ALFALFA (Code 1, bis cz = 15 000 km/s, 23 956 Quellen über der Grenze, an der die Hälfte erfasst ist)", float(D.alfalfa_hi_dichte()[1]), 23956),
        ("F:ganzes", "Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 %", D.gasanteil_heute(), 1.21),
        ("F:ganzes", "Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 %", min(D.gasanteil_heute(z, g) for z in (1.0, 2.0) for g in ('90', '50', '25')), 1.12),
        ("F:ganzes", "Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 %", max(D.gasanteil_heute(z, g) for z in (1.0, 2.0) for g in ('90', '50', '25')), 1.47),
        ("F:ganzes", "Mit Helium (Y_p = 0,2449) liegt heute F = 1,21 % der Baryonen als neutrales Gas in Galaxien, 1,12 bis 1,47 %", K.YP_GEMESSEN[0], 0.2449),
        ("F:ganzes", "(gezählt in Zellen von 2° oder 1°: 6801 oder 6095 Quadratgrad)", D.alfalfa_flaeche(2.0), 6801),
        ("F:ganzes", "(gezählt in Zellen von 2° oder 1°: 6801 oder 6095 Quadratgrad)", D.alfalfa_flaeche(1.0), 6095),
        ("F:ganzes", "Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 %", float(D.hipass_hi_dichte()[1]), 1171),
        ("F:ganzes", "Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 %", D.gasanteil_hipass(), 1.17),
        ("F:ganzes", "Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 %", min(D.gasanteil_hipass(a, b) for a, b in ((0.068, 7.4), (0.084, 9.4), (0.1, 12.0), (0.12, 15.0))), 1.17),
        ("F:ganzes", "Unabhängig davon gibt HIPASS (Parkes, Südhimmel, 1171 Quellen über den Grenzen, an denen 99 % erfasst sind, Abstand aus der Geschwindigkeit im Rahmen der Lokalen Gruppe) F = 1,17 %, mit strengeren oder lockereren Grenzen 1,17 bis 1,22 %", max(D.gasanteil_hipass(a, b) for a, b in ((0.068, 7.4), (0.084, 9.4), (0.1, 12.0), (0.12, 15.0))), 1.22),
    ]


@lru_cache(None)
def _gipfel():
    return D.sparc_zusatz_gipfel(K.kpc)


def _gipfel_standard():
    from scipy.optimize import minimize_scalar
    import numpy as np
    r = minimize_scalar(lambda x: -x * (np.sqrt(0.5 + np.sqrt(0.25 + 1 / x**2)) - 1), bounds=(0.05, 10), method="bounded",
                        options={"xatol": 1e-10})
    return float(r.x), float(-r.fun)


def _rc100_einzeln():
    """Schwelle aus jeder Galaxie von RC100 einzeln, a₀ = g_N/(ln f_DM)² [10⁻¹⁰ m/s²], der zweite Weg: (Median, Zahl
    der Galaxien über 1,06)."""
    import numpy as np
    d = D.rc100()
    gN = (1 - d[:, 5]) * (d[:, 7] * 1e3)**2 / (d[:, 3] * 3.0857e19)
    a = gN / np.log(d[:, 5])**2 / 1e-10
    return float(np.median(a)), float(np.sum(a > 1.06))


def _rc100_baryonen():
    """Das g_N von RC100 als Feld der Baryonen: (Median von V_bar²·R_e/(GM_b) mit V_bar² = (1 − f_DM)·V_c², derselbe Wert
    der dünnen Scheibe am Radius der halben Masse, x = R_e/R_d aus 1 − (1 + x)·e^{−x} = ½, 4y³·[I₀K₀ − I₁K₁](y), y = x/2)."""
    import numpy as np
    from scipy import special as S
    from scipy.optimize import brentq
    d = D.rc100()
    r = (1 - d[:, 5]) * (d[:, 7] * 1e3)**2 * d[:, 3] * K.kpc / (K.G * 10**d[:, 1] * K.M_sun)
    y = brentq(lambda x: 1 - (1 + x) * math.exp(-x) - 0.5, 0.5, 3.0) / 2
    return float(np.median(r)), float(4 * y**3 * (S.i0(y) * S.k0(y) - S.i1(y) * S.k1(y)))


def _hi_gang():
    """Gang des neutralen Gases von heute bis z ≈ 0,34 gegen den, den der Zustandsteil verlangt: (gemessenes Verhältnis,
    Fehler, verlangt bei Sättigung beim Standardwert, Abstand in Fehlern, verlangt in der eigenen Form, Abstand)."""
    z, rho, e1, e2 = D.hi_frueh()
    heute = D.alfalfa_hi_dichte(H0=70.0, flaeche=D.alfalfa_flaeche(2.0))[0]
    r = rho / heute
    er = r * math.hypot(e1, e2) / rho
    aus = [r, er]
    for (f0, e0), (fm, em) in ((_cmb()["zustand_form"][1], _cmb()["zustand_form"][2]), (_gas_form()["F0"], _gas_form()["Fs"])):
        v = fm / f0
        ev = v * math.hypot(em / fm, e0 / f0)
        aus += [v, (v - r) / math.hypot(ev, er)]
    return tuple(aus)


def _zelle():
    """L = (3c⁴/(2G·ρ_m·a₀))^{1/3}, ρ_m = Ω_m·3H₀²/(8πG)."""
    rho_m = (K.Omega_b + K.Omega_c) * 3 * K.H0**2 / (8 * math.pi * K.G)
    return (3 * K.c**4 / (2 * K.G * rho_m * K.a0))**(1 / 3) / K.Gpc


def _modelltext(nr):
    """Text des Bausteins nr aus modell/; „F:name“ ist die Forschungsnotiz forschung/name.md (Rechnung, noch nicht im Modell),
    „A:name“ die Ausarbeitung abhandlung/name.md."""
    if nr.startswith("F:"):
        pfad = os.path.join(WURZEL, "forschung", nr[2:] + ".md")
        return open(pfad, encoding="utf-8").read() if os.path.exists(pfad) else ""
    if nr.startswith("A:"):                                       # Ausarbeitung in abhandlung/ (Wardenclyffe-Turm)
        pfad = os.path.join(WURZEL, "abhandlung", nr[2:] + ".md")
        return open(pfad, encoding="utf-8").read() if os.path.exists(pfad) else ""
    for f in os.listdir(os.path.join(WURZEL, "modell")):
        if f.startswith(nr):
            return open(os.path.join(WURZEL, "modell", f), encoding="utf-8").read()
    return ""


def _stellen(x):
    s = f"{x:g}"
    return len(s.split(".")[1]) if "." in s else 0


def pruefe_modellzahlen(voll=False):
    """Prüft jede Zahl der Prüfliste gegen den Text. Gerechnet wird nur, was sich seit dem letzten Lauf geändert hat
    (Speicher der Prüfwerte), voll: alles neu."""
    alt = {}
    if not voll and os.path.exists(SPEICHER):
        with open(SPEICHER, encoding="utf-8") as f:
            alt = json.load(f)
    neu, gerechnet, zeilen = {}, 0, []
    for nr, stelle, knoten, schluessel, soll in _pruefliste():
        if schluessel in alt:
            wert = alt[schluessel]
        else:
            wert = eval(compile(ast.Expression(knoten), "<prüfliste>", "eval"), globals())
            wert = _zelle() if wert is None else float(wert)
            gerechnet += 1
        neu[schluessel] = wert
        steht = stelle in _modelltext(nr)
        rund = round(wert, _stellen(soll))
        passt = abs(rund - soll) <= 10**(-_stellen(soll)) * 0.5 + 1e-12
        zeilen.append((nr, stelle, wert, steht and passt, steht))
    os.makedirs(os.path.dirname(SPEICHER), exist_ok=True)
    with open(SPEICHER, "w", encoding="utf-8") as f:
        json.dump(neu, f)
    print(f"Prüfliste: {len(zeilen)} Werte, {gerechnet} neu gerechnet")
    zeilen += pruefe_abhandlung()
    zeilen += pruefe_rechnen(voll)
    zeilen += pruefe_elemente(voll)
    return zeilen


def pruefe_abhandlung():
    """Jede Dezimalzahl der Abhandlung steht im geprüften Modell (Zahlen nur aus dem Modell), Namen nach Himmelslage
    (J0025.4) ausgenommen."""
    modell = "".join(open(os.path.join(WURZEL, "modell", f), encoding="utf-8").read()
                     for f in os.listdir(os.path.join(WURZEL, "modell")))
    zeilen = []
    for name, muster, komma in (("abhandlung_de.md", r"(?<![\dJ])\d+,\d+", False), ("abhandlung_en.md", r"(?<![\dJ])\d+\.\d+", True)):
        pfad = os.path.join(WURZEL, "abhandlung", name)
        if not os.path.exists(pfad):
            continue
        zahlen = sorted(set(re.findall(muster, open(pfad, encoding="utf-8").read())))
        fehlt = [z for z in zahlen if (z.replace(".", ",") if komma else z) not in modell]
        zeilen += [("AB", f"{name}: Zahl {z} steht im Modell", 0.0, False, False) for z in fehlt] or \
                  [("AB", f"{name}: alle {len(zahlen)} Zahlen stehen im Modell", float(len(zahlen)), True, True)]
    return zeilen


def pruefe_rechnen(voll=False):
    """Jede gerechnete Zahl der Beispiele steht in beiden Sprachen; jede Dezimalzahl dort ist gerechnet oder im Modell."""
    gerechnet = _gespeichert("beispiele", [("beispiele", "zahlen")], B.zahlen, voll)
    modell = "".join(open(os.path.join(WURZEL, "modell", f), encoding="utf-8").read()
                     for f in os.listdir(os.path.join(WURZEL, "modell")))
    zeilen = []
    for name, muster, komma in (("rechnen_de.md", r"\d+,\d+", True), ("rechnen_en.md", r"\d+\.\d+", False)):
        text = open(os.path.join(WURZEL, "abhandlung", name), encoding="utf-8").read()
        soll = [z.replace(".", ",") if komma else z for z in gerechnet]
        fehlt = [z for z in soll if not re.search(r"(?<![\d.,])" + re.escape(z) + r"(?![\d])", text)]
        fremd = [z for z in set(re.findall(muster, text)) if z not in soll and z.replace(".", ",") not in modell]
        zeilen += [("RB", f"{name}: gerechnete Zahl {z} fehlt", 0.0, False, False) for z in fehlt]
        zeilen += [("RB", f"{name}: Zahl {z} ist weder gerechnet noch im Modell", 0.0, False, False) for z in fremd]
        if not fehlt and not fremd:
            zeilen.append(("RB", f"{name}: alle {len(set(soll))} gerechneten Zahlen stehen", float(len(set(soll))), True, True))
    return zeilen


def pruefe_elemente(voll=False):
    """abhandlung/elemente_de.md ist genau die Rechnung von elemente.py."""
    import elemente
    aus = []
    for sp in ("de", "en"):
        gerechnet = _gespeichert(f"elemente_{sp}", [("elemente", "text")], lambda: elemente.text(sp), voll)
        gleich = open(elemente.ZIEL[sp], encoding="utf-8").read() == gerechnet
        aus.append(("EL", f"elemente_{sp}.md entspricht der Rechnung" if gleich else f"elemente_{sp}.md weicht von der Rechnung ab", float(gleich), gleich, gleich))
    return aus


def naechste_nummer(ordner):
    os.makedirs(ordner, exist_ok=True)
    nummern = [int(m.group(1)) for f in os.listdir(ordner) if (m := re.match(r"bericht_(\d{4})(\.md)?$", f))]
    return max(nummern, default=0) + 1


HOCH = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def zahl(x):
    """Zahl mit fünf gültigen Stellen, Komma, Zehnerpotenz mit Hochzahlen."""
    if x == 0:
        return "0"
    if 1e-3 <= abs(x) < 1e5:
        return f"{x:.5g}".replace(".", ",")
    m, ex = f"{x:.4e}".split("e")
    m = m.rstrip("0").rstrip(".").replace(".", ",")
    return f"{m}·10{str(int(ex)).translate(HOCH)}"


def speichere_daten(ordner, name, erg):
    """CSV: t [s], je Körper x, y, z [m]; szene.json beschreibt Körper und Einheiten."""
    T, X, w = erg["T"], erg["X"], erg["welt"]
    kopf = ["t_s"] + [f"{k.name.replace(',', '.').replace(' ', '_')}_{a}_m" for k in w.koerper for a in "xyz"]
    with open(os.path.join(ordner, f"{name}.csv"), "w", encoding="utf-8") as f:
        f.write(",".join(kopf) + "\n")
        for t, x in zip(T, X):
            f.write(",".join([repr(float(t))] + [repr(float(v)) for v in x.reshape(-1)]) + "\n")
    return dict(datei=f"{name}.csv", titel=erg["titel"], text=erg["text"],
                koerper=[dict(name=k.name, masse_kg=k.masse, ladung_C=k.ladung) for k in w.koerper],
                hintergrund_masse_kg=w.hintergrund, regel=w.regel, takt_glied=w.takt_glied,
                einheiten=dict(zeit="s", ort="m"), proben=len(T))


def main():
    t0 = time.time()
    berichte = os.path.join(WURZEL, "berichte")
    nr = naechste_nummer(berichte)
    name = f"bericht_{nr:04d}"
    daten = os.path.join(berichte, name)
    os.makedirs(daten, exist_ok=False)
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=WURZEL, capture_output=True,
                                text=True).stdout.strip()
    except OSError:
        commit = "?"

    L = [f"# Bericht {nr}", "",
         f"{datetime.datetime.now().isoformat(timespec='seconds')} · Commit {commit} · a₀ = {zahl(K.a0)} m/s², "
         f"H₀ = {zahl(K.H0_kms)} km/s/Mpc", ""]
    szene = []
    fehler = 0
    voll = "--voll" in sys.argv
    for fn in szenarien.ALLE:
        erg = _gespeichert(f"szene_{fn.__name__}", [("szenarien", fn.__name__)], fn, voll)
        szene.append(speichere_daten(daten, fn.__name__, erg))
        L += [f"## {erg['titel']}", "", erg["text"], "",
              "| Aussage | Modell | Simulation | Einheit | Urteil |", "|---|---|---|---|---|"]
        for aussage, soll, ist, einheit, tol in erg["pruefungen"]:
            if tol is None:
                ok = ist <= soll
            else:
                ok = abs(ist - soll) <= tol * abs(soll)
            fehler += not ok
            L.append(f"| {aussage} | {zahl(soll)} | {zahl(ist)} | {einheit} | {'OK' if ok else 'ABWEICHUNG'} |")
        L.append("")
    with open(os.path.join(daten, "szene.json"), "w", encoding="utf-8") as f:
        json.dump(dict(bericht=nr, szenarien=szene), f, ensure_ascii=False, indent=1)

    L += ["## Zahlen im Modelltext", "", "| Baustein | Stelle | Kern | Urteil |", "|---|---|---|---|"]
    for nrb, stelle, wert, ok, steht in pruefe_modellzahlen(voll):
        fehler += not ok
        urteil = "OK" if ok else ("ABWEICHUNG" if steht else "FEHLT IM TEXT")
        L.append(f"| {nrb} | {stelle} | {zahl(wert)} | {urteil} |")
    L += ["", f"**Ergebnis:** {'alles OK' if fehler == 0 else f'{fehler} Abweichung(en)'} · "
              f"Bahndaten in `berichte/{name}/` · Laufzeit {time.time() - t0:.0f} s", ""]
    with open(os.path.join(berichte, f"{name}.md"), "x", encoding="utf-8") as f:
        f.write("\n".join(L))
    print(f"{name}.md: {'alles OK' if fehler == 0 else f'{fehler} Abweichung(en)'} ({time.time() - t0:.0f} s)")
    return fehler


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
