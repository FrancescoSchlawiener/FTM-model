# -*- coding: utf-8 -*-
"""Unsere Daten, nur gelesen: SPARC (Lelli, McGaugh, Schombert 2016) und AME2020 (Huang u. a. 2021)."""
import os

import numpy as np

ORT = os.path.join(os.path.dirname(__file__), "..", "archiv", "weltmodell", "daten")


def sparc_galaxien():
    """Tafel der 175 Galaxien: Name → dict(T, D, Inc, L36, Rd, MHI, Vflat, Q)."""
    aus = {}
    zeilen = open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1").readlines()
    for z in zeilen:
        p = z.split()
        if len(p) < 18:
            continue
        try:
            aus[p[0]] = dict(T=int(p[1]), D=float(p[2]), fD=int(p[4]), Inc=float(p[5]), L36=float(p[7]), Rd=float(p[11]),
                             MHI=float(p[13]), Vflat=float(p[15]), Q=int(p[17]))
        except ValueError:
            continue
    return aus


def sparc_kurve(name):
    """Rotationskurve: Felder R [kpc], Vobs, eV, Vgas, Vdisk, Vbul [km/s], SBdisk, SBbul [L☉/pc²]."""
    d = np.loadtxt(os.path.join(ORT, "sparc", "Rotmod_LTG", name + "_rotmod.dat"), comments="#")
    return dict(R=d[:, 0], V=d[:, 1], eV=d[:, 2], Vgas=d[:, 3], Vdisk=d[:, 4], Vbul=d[:, 5], SBdisk=d[:, 6], SBbul=d[:, 7])


def ame2020():
    """Kerne aus AME2020: Liste (Z, N, A, Bindung je Nukleon [keV], atomare Masse [µu])."""
    aus = []
    for z in open(os.path.join(ORT, "elemente", "mass_1.mas20.txt"), encoding="latin-1"):
        if len(z) < 130 or not z[1:4].strip().lstrip("-").isdigit():
            continue
        try:
            N, Z, A = int(z[4:9]), int(z[9:14]), int(z[14:19])
            B = float(z[54:67].replace("#", "").strip())
            m = float((z[106:109] + z[110:123]).replace("#", "").replace(" ", ""))
        except ValueError:
            continue
        aus.append((Z, N, A, B, m))
    return aus


# ---------------------------------------------------------------- Prüfungen an SPARC
import math
from functools import lru_cache

YD, YB, STREUUNG = 0.5, 0.7, 0.05      # Masse-zu-Licht Scheibe, Bulge; Streuung der Regel selbst [dex]
FORMEN = {"Bose": lambda x: 1 / (-np.expm1(-np.sqrt(x))),
          "einfach": lambda x: 0.5 + np.sqrt(0.25 + 1 / x),
          "standard": lambda x: np.sqrt(0.5 + np.sqrt(0.25 + 1 / x**2))}
GITTER = np.linspace(-10.2, -9.7, 11)  # log₁₀ a₀


H0_SPARC = 73.0                        # km/s/Mpc, mit dem SPARC die Abstände aus dem Hubble-Fluss setzt (f_D = 1)


@lru_cache(maxsize=None)
def sparc_auswahl(yd=YD, H0=None, wirbelfrei=False):
    """Galaxien mit Q ≤ 2 und Neigung ≥ 30°, Punkte mit V > 0, g_bar > 0, Fehler ≤ 10 %; mit Fehlern von Abstand und Neigung.
    H0: Abstände aus dem Hubble-Fluss mit diesem H₀ statt 73 (R ∝ D, V_bar ∝ √D, V bleibt). wirbelfrei: je Punkt die Tafel
    lg(wirbelfrei/(ν·g_N)) über Υ_d und a₀ (sparc_wirbelfrei_tafel) im Feld "wf"."""
    tafel = sparc_wirbelfrei_tafel() if wirbelfrei else None
    eD, eI = {}, {}
    for z in open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
        p = z.split()
        if len(p) >= 18:
            try:
                eD[p[0]], eI[p[0]] = float(p[3]), float(p[6])
            except ValueError:
                pass
    aus = []
    for n, p in sparc_galaxien().items():
        if p["Q"] > 2 or p["Inc"] < 30:
            continue
        k, f = sparc_kurve(n), 1.0
        if H0 is not None and p["fD"] == 1:
            f = H0_SPARC / H0
            p, k = dict(p, D=p["D"] * f), dict(k, R=k["R"] * f, **{x: k[x] * math.sqrt(f) for x in ("Vgas", "Vdisk", "Vbul")})
        vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + YB * k["Vbul"]**2
        m = (k["V"] > 0) & (vb2 > 0) & (k["eV"] / np.maximum(k["V"], 1e-9) <= 0.1)
        if m.sum() >= 3:
            aus.append((n, p, dict({x: k[x][m] for x in k}, **({"wf": tafel[n][:, :, m]} if wirbelfrei else {})),
                        eD[n] * f, eI[n]))
    return tuple(aus)


def _wirbelfrei_an(wf, ud, a0):
    """lg(wirbelfrei/(ν·g_N)) je Punkt bei Υ_d = ud und a₀, linear in lg Υ_d und lg a₀ zwischen den Stellen der Tafel."""
    i = float(np.interp(math.log10(ud), np.log10(YD_TAFEL), np.arange(len(YD_TAFEL))))
    j = float(np.interp(math.log10(a0), GITTER, np.arange(len(GITTER))))
    i0, j0 = min(int(i), len(YD_TAFEL) - 2), min(int(j), len(GITTER) - 2)
    ti, tj = i - i0, j - j0
    return ((1 - ti) * ((1 - tj) * wf[i0, j0] + tj * wf[i0, j0 + 1])
            + ti * ((1 - tj) * wf[i0 + 1, j0] + tj * wf[i0 + 1, j0 + 1]))


def _reste(v, a0, nu, p, k, eD, eI, kpc, yd, halo=None, roh=False):
    """Reste der Regel in log g mit Abstand, Neigung und Υ_d als Störgrößen (Vorgaben als Gauß-Glieder).
    halo = H₀ [km/s/Mpc]: statt der Regel das Standardmodell, Faltung der Baryonen plus NFW-Halo mit v[3] = lg V₂₀₀ [km/s],
    v[4] = lg c (R₂₀₀ = V₂₀₀/(10·H₀)). roh: Reste je Punkt in dex und lg(g_bar/a₀), ohne Störglieder. Trägt k das Feld
    "wf", gilt die wirbelfreie Regel."""
    R = k["R"] * kpc
    inc = math.radians(p["Inc"])
    ld, di, lu = v[:3]
    dl, i2, ud = 10**ld, inc + di, yd * 10**lu
    vb2 = k["Vgas"] * np.abs(k["Vgas"]) + ud * k["Vdisk"]**2 + YB * k["Vbul"]**2
    gb = np.maximum(vb2, 1e-6) * 1e6 / R
    f = math.sin(inc) / math.sin(i2)
    Vo, eV = k["V"] * f, k["eV"] * f
    go = Vo**2 * 1e6 / (R * dl)
    sg = np.sqrt((2 * eV / (Vo * math.log(10)))**2 + STREUUNG**2)
    if halo is None:
        gm = nu(gb / a0) * gb * (10**_wirbelfrei_an(k["wf"], ud, a0) if "wf" in k else 1.0)
    else:
        v2, cc = 10**v[3], 10**v[4]
        x = k["R"] * dl / (v2 / (10 * halo / 1e3))
        m = lambda y: np.log1p(y) - y / (1 + y)
        gm = gb + v2**2 * m(cc * x) / (x * m(cc)) * 1e6 / (R * dl)
    if roh:
        return np.log10(go) - np.log10(gm), (np.log10(gb / a0) if a0 else None)
    return np.concatenate([(np.log10(go) - np.log10(gm)) / sg,
                           [(dl - 1) / (eD / p["D"]), di / math.radians(max(eI, 1)), lu / 0.1]])


@lru_cache(maxsize=None)
def sparc_profil(form, kpc, yd=YD, H0=None, wirbelfrei=False):
    """χ² je Galaxie über dem Gitter von a₀, Störgrößen je Punkt des Gitters ausgeglichen (warmer Start)."""
    from scipy.optimize import least_squares
    nu = FORMEN[form]
    C = np.zeros((len(sparc_auswahl(yd, H0, wirbelfrei)), len(GITTER)))
    for gi, (n, p, k, eD, eI) in enumerate(sparc_auswahl(yd, H0, wirbelfrei)):
        inc = math.radians(p["Inc"])
        lo = np.array([-0.4, max(-0.6, 0.25 - inc), -0.6])
        hi = np.array([0.4, max(1e-6, math.pi / 2 - inc) + 1e-9, 0.6])
        x0 = np.zeros(3)
        for j, la in enumerate(GITTER):
            r = least_squares(_reste, np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi),
                              args=(10**la, nu, p, k, eD, eI, kpc, yd), xtol=1e-10, ftol=1e-10)
            C[gi, j], x0 = 2 * r.cost, r.x
    return C


@lru_cache(maxsize=None)
def sparc_halo(kpc, H0, yd=YD, zahlen=False):
    """Standardmodell an SPARC: je Galaxie die Faltung der Baryonen plus ein NFW-Halo mit zwei freien Zahlen (V₂₀₀, c),
    dieselben Störgrößen wie die Regel; χ² je Galaxie (bester von mehreren Starts). zahlen: je Galaxie (χ², lg V₂₀₀, lg c)."""
    from scipy.optimize import least_squares
    aus = []
    for n, p, k, eD, eI in sparc_auswahl(yd):
        inc = math.radians(p["Inc"])
        lo = np.array([-0.4, max(-0.6, 0.25 - inc), -0.6, 0.5, 0.0])
        hi = np.array([0.4, max(1e-6, math.pi / 2 - inc) + 1e-9, 0.6, 3.0, 2.0])
        vm = math.log10(max(k["V"].max(), 5.0))
        best, x = math.inf, None
        for lv in (vm - 0.2, vm, vm + 0.3):
            for lc in (0.3, 1.0, 1.6):
                x0 = np.clip([0, 0, 0, lv, lc], lo + 1e-9, hi - 1e-9)
                r = least_squares(_reste, x0, bounds=(lo, hi), args=(None, None, p, k, eD, eI, kpc, yd, H0),
                                  xtol=1e-10, ftol=1e-10)
                if 2 * r.cost < best:
                    best, x = 2 * r.cost, r.x
        aus.append((best, x[3], x[4]) if zahlen else best)
    return np.array(aus)


def sparc_reste(kpc, a0, yd=YD, wirbelfrei=False):
    """Reste der Regel je Punkt in dex bei der Schwelle a₀, Störgrößen je Galaxie ausgeglichen; dazu lg(g_bar/a₀)."""
    from scipy.optimize import least_squares
    nu, x, r = FORMEN["Bose"], [], []
    for n, p, k, eD, eI in sparc_auswahl(yd, None, wirbelfrei):
        inc = math.radians(p["Inc"])
        lo = np.array([-0.4, max(-0.6, 0.25 - inc), -0.6])
        hi = np.array([0.4, max(1e-6, math.pi / 2 - inc) + 1e-9, 0.6])
        v = least_squares(_reste, np.zeros(3), bounds=(lo, hi), args=(a0, nu, p, k, eD, eI, kpc, yd)).x
        d, g = _reste(v, a0, nu, p, k, eD, eI, kpc, yd, roh=True)
        r += list(d); x += list(g)
    return np.array(x), np.array(r)


@lru_cache(maxsize=None)
def sparc_wirbelfrei_gang(kpc, a0, yd=YD, ziehungen=500):
    """Gang der Reste der Regel ν(g_N)·g_N gegen die wirbelfreie Form je Punkt, beide innerhalb der Galaxien zentriert
    (Störgrößen je Galaxie mit der Regel als Vektor ausgeglichen): Steigung 1 heißt, die Daten folgen der wirbelfreien
    Form, 0 der Regel als Vektor. Gibt (Steigung, Streuung über Ziehen der Galaxien, Zahl der Galaxien)."""
    from scipy.optimize import least_squares
    nu, x, y = FORMEN["Bose"], [], []
    for n, p, k, eD, eI in sparc_auswahl(yd, None, True):
        kv = {s: k[s] for s in k if s != "wf"}
        inc = math.radians(p["Inc"])
        lo = np.array([-0.4, max(-0.6, 0.25 - inc), -0.6])
        hi = np.array([0.4, max(1e-6, math.pi / 2 - inc) + 1e-9, 0.6])
        v = least_squares(_reste, np.zeros(3), bounds=(lo, hi), args=(a0, nu, p, kv, eD, eI, kpc, yd)).x
        d = _reste(v, a0, nu, p, kv, eD, eI, kpc, yd, roh=True)[0]
        w = _wirbelfrei_an(k["wf"], yd * 10**v[2], a0)
        x.append(w - w.mean()); y.append(d - d.mean())
    st = lambda i: np.polyfit(np.concatenate([x[j] for j in i]), np.concatenate([y[j] for j in i]), 1)[0]
    rng = np.random.default_rng(8)
    return st(range(len(x))), float(np.std([st(rng.integers(0, len(x), len(x))) for _ in range(ziehungen)])), len(x)


def _bestes(v):
    from scipy.interpolate import CubicSpline
    from scipy.optimize import minimize_scalar
    s = CubicSpline(GITTER, v)
    r = minimize_scalar(s, bounds=(GITTER[0], GITTER[-1]), method="bounded")
    return 10**r.x, float(s(r.x))


@lru_cache(maxsize=None)
def sparc_schwelle(kpc, form="Bose", ziehungen=1000, yd=YD, H0=None, wirbelfrei=False):
    """a₀ aus SPARC mit der Regel: bester Wert und 16–84 % über Ziehen der Galaxien."""
    C = sparc_profil(form, kpc, yd, H0, wirbelfrei)
    rng = np.random.default_rng(5)
    bs = [_bestes(C[rng.integers(0, len(C), len(C))].sum(0))[0] for _ in range(ziehungen)]
    return _bestes(C.sum(0))[0], tuple(np.percentile(bs, [16, 84]))


@lru_cache(maxsize=None)
def sparc_formen(kpc, yd=YD):
    """χ² der drei Formen im jeweils besten a₀, je Punkt der Daten."""
    n = sum(len(g[2]["R"]) for g in sparc_auswahl(yd))
    return {f: _bestes(sparc_profil(f, kpc, yd).sum(0))[1] for f in FORMEN}, n


def _einfache_punkte(kpc, yd=YD, H0=None):
    gb, go, gal, fg = [], [], [], []
    for n, p, k, _, _ in sparc_auswahl(yd, H0):
        R = k["R"] * kpc
        vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + YB * k["Vbul"]**2
        gb += list(vb2 * 1e6 / R); go += list(k["V"]**2 * 1e6 / R); gal += [n] * len(R)
        fg += list(k["Vgas"] * np.abs(k["Vgas"]) / vb2)
    return np.array(gb), np.array(go), np.array(gal), np.array(fg)


@lru_cache(maxsize=None)
def sparc_dimension(kpc, ziehungen=200, yd=YD):
    """Dimension der Rückkehr: q = exp(−(a_N/a₀)^p), d_s = 2(1 − p); Mittel und Streuung über Ziehen der Galaxien."""
    from scipy.optimize import minimize
    gb, go, gal, _ = _einfache_punkte(kpc, yd)
    namen = np.unique(gal)
    idx = {n: np.where(gal == n)[0] for n in namen}

    def fit(w):
        f = lambda v: np.sum((np.log10(go[w]) - np.log10(gb[w] / (-np.expm1(-(gb[w] / 10**v[0])**v[1]))))**2)
        return 2 * (1 - minimize(f, [-10, 0.5], method="Nelder-Mead").x[1])
    rng = np.random.default_rng(6)
    bs = [fit(np.concatenate([idx[n] for n in rng.choice(namen, len(namen))])) for _ in range(ziehungen)]
    return fit(np.arange(len(gb))), float(np.std(bs))


@lru_cache(maxsize=None)
def sparc_gasanteil(kpc, a0, ziehungen=500, yd=YD, H0=None):
    """Steigung des Rests gegen den örtlichen Gasanteil innerhalb der Galaxien (je Galaxie zentriert), mit Streuung."""
    gb, go, gal, fg = _einfache_punkte(kpc, yd, H0)
    r = np.log10(go) - np.log10(FORMEN["Bose"](gb / a0) * gb)
    namen = [n for n in np.unique(gal) if (gal == n).sum() >= 5]
    x = {n: fg[gal == n] - fg[gal == n].mean() for n in namen}
    y = {n: r[gal == n] - r[gal == n].mean() for n in namen}
    st = lambda ns: np.polyfit(np.concatenate([x[n] for n in ns]), np.concatenate([y[n] for n in ns]), 1)[0]
    rng = np.random.default_rng(7)
    return st(namen), float(np.std([st(rng.choice(namen, len(namen))) for _ in range(ziehungen)]))


@lru_cache(maxsize=None)
def sparc_gruppe(kpc, yd=YD, fD=4, ziehungen=1000):
    """Galaxien einer Gruppe (fD = 4: Abstand aus der Ursa-Major-Gruppe) gegen die übrigen: a₀ je Teil mit 16–84 % über
    Ziehen der Galaxien, und der mittlere Rest der Punkte mit g_bar < a₀/10 bei der Schwelle aller, je Teil mit Streuung
    über Ziehen der Galaxien: (Zahl der Gruppe, a₀ Gruppe, (16, 84), a₀ übrige, (16, 84), Rest Gruppe, Streuung, Rest
    übrige, Streuung, Zahl der äußeren Punkte der Gruppe)."""
    C = sparc_profil("Bose", kpc, yd)
    auswahl = sparc_auswahl(yd)
    drin = np.array([p["fD"] == fD for _, p, _, _, _ in auswahl])
    rng = np.random.default_rng(8)
    aus = [int(drin.sum())]
    for w in (drin, ~drin):
        Cw = C[w]
        bs = [_bestes(Cw[rng.integers(0, len(Cw), len(Cw))].sum(0))[0] for _ in range(ziehungen)]
        aus += [_bestes(Cw.sum(0))[0], tuple(np.percentile(bs, [16, 84]))]
    a0 = _bestes(C.sum(0))[0]
    gb, go, gal, _ = _einfache_punkte(kpc, yd)
    r = np.log10(go) - np.log10(FORMEN["Bose"](gb / a0) * gb)
    aussen = gb < a0 / 10
    gruppe = {n for (n, p, _, _, _) in auswahl if p["fD"] == fD}
    for w in (True, False):
        namen = [n for n in np.unique(gal[aussen]) if (n in gruppe) == w]
        je = {n: r[aussen & (gal == n)] for n in namen}
        mittel = lambda ns: float(np.concatenate([je[n] for n in ns]).mean())
        aus += [mittel(namen), float(np.std([mittel(rng.choice(namen, len(namen))) for _ in range(ziehungen)]))]
    aus.append(int(sum((aussen & (gal == n)).sum() for n in gruppe)))
    return tuple(aus)


@lru_cache(maxsize=None)
def sparc_durchlaeufe(kpc, a0, yd=YD, alter_gyr=13.8, je_radiant=False, ziehungen=1000):
    """Baut sich die Verstärkung in Durchläufen auf? Nach N Rückläufen ist die Summe erst (1 − q^N)/(1 − q). Gezählt
    als Umläufe N = Alter/(2πR/V) (je_radiant: 2π·N): vorhergesagt ist der Rest δ = lg(1 − q^N) gegen die stationäre
    Regel, q = e^{−√(g_bar/a₀)}. Steigung k des gemessenen Rests gegen δ, beide je Galaxie zentriert (k = 0: stationär,
    k = 1: ein Rücklauf je Umlauf), mit Streuung über Ziehen der Galaxien, dazu die Zahl der Punkte mit δ < −0,05, deren
    mittlerer Rest und mittleres δ."""
    gb, go, gal, R, V = [], [], [], [], []
    for n, p, k, _, _ in sparc_auswahl(yd):
        vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + YB * k["Vbul"]**2
        gb += list(vb2 * 1e6 / (k["R"] * kpc)); go += list(k["V"]**2 * 1e6 / (k["R"] * kpc))
        gal += [n] * len(k["R"]); R += list(k["R"] * kpc); V += list(k["V"] * 1e3)
    gb, go, gal, R, V = map(np.array, (gb, go, gal, R, V))
    N = alter_gyr * 3.156e16 / (2 * math.pi * R / V) * (2 * math.pi if je_radiant else 1)
    q = np.exp(-np.sqrt(gb / a0))
    d = np.log10(-np.expm1(N * np.log(q)))
    r = np.log10(go) - np.log10(FORMEN["Bose"](gb / a0) * gb)
    namen = [n for n in np.unique(gal) if (gal == n).sum() >= 3]
    x = {n: d[gal == n] - d[gal == n].mean() for n in namen}
    y = {n: r[gal == n] - r[gal == n].mean() for n in namen}
    st = lambda ns: float(np.polyfit(np.concatenate([x[n] for n in ns]), np.concatenate([y[n] for n in ns]), 1)[0])
    rng = np.random.default_rng(9)
    w = d < -0.05
    return (st(namen), float(np.std([st(rng.choice(namen, len(namen))) for _ in range(ziehungen)])), int(w.sum()),
            float(r[w].mean()), float(d[w].mean()))


# ---------------------------------------------------------------- Prüfungen an AME2020


@lru_cache(maxsize=None)
def sparc_eigenschaften(kpc, a0, yd=YD):
    """Rangkorrelation des mittleren Rests je Galaxie mit baryonischer Masse, Gasanteil, Hubble-Typ und Flächenhelligkeit der Scheibe: (ρ, p)."""
    from scipy.stats import spearmanr
    gb, go, gal, _ = _einfache_punkte(kpc, yd)
    r = np.log10(go) - np.log10(FORMEN["Bose"](gb / a0) * gb)
    tafel = {n: p for n, p, _, _, _ in sparc_auswahl(yd)}
    namen = list(tafel)
    rest = [r[gal == n].mean() for n in namen]
    mb = [math.log10(1.33 * tafel[n]["MHI"] + yd * tafel[n]["L36"]) for n in namen]
    fg = [1.33 * tafel[n]["MHI"] / (1.33 * tafel[n]["MHI"] + yd * tafel[n]["L36"]) for n in namen]
    sb = {}
    for z in open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
        t = z.split()
        if t[:1] and t[0] in tafel:
            sb[t[0]] = float(t[12])
    typ = [tafel[n]["T"] for n in namen]
    hell = [sb[n] for n in namen]
    return {g: tuple(float(v) for v in spearmanr(rest, w))
            for g, w in (("Masse", mb), ("Gasanteil", fg), ("Typ", typ), ("Flächenhelligkeit", hell))}


@lru_cache(maxsize=None)
def sparc_eichung(kpc, H0=None):
    """Masse-zu-Licht der Scheibe so, dass Stern und Gas derselben Regel folgen (kein Gang des Rests mit dem Gasanteil
    innerhalb der Galaxien): Die Regel wirkt auf die Faltung, nicht auf den Zustand der Materie."""
    from scipy.optimize import brentq
    def gang(yd):
        return sparc_gasanteil(kpc, sparc_schwelle(kpc, ziehungen=200, yd=yd, H0=H0)[0], yd=yd, H0=H0)[0]
    return round(brentq(gang, 0.5, 0.8, xtol=2e-3), 3)


def pubchem():
    """Periodensystem aus PubChem: Symbol → dict(Z, Ionisierungsenergie [eV], Schmelz- und Siedepunkt [K])."""
    import csv
    aus = {}
    with open(os.path.join(ORT, "elemente", "pubchem_periodensystem.csv"), encoding="utf-8") as f:
        for z in csv.DictReader(f):
            zahl = lambda s: float(s) if s not in ("", None) else None
            aus[z["Symbol"]] = dict(Z=int(z["AtomicNumber"]), IE=zahl(z["IonizationEnergy"]),
                                    Ts=zahl(z["MeltingPoint"]), Tb=zahl(z["BoilingPoint"]))
    return aus


@lru_cache(maxsize=None)
def sparc_scheiben(kpc, yd, M_sonne, a0, G, ziehungen=1000):
    """Scheibenlänge gegen Schwellenradius r* = √(GM_b/a₀): Median und 16–84 % von R_d/r*, Steigung von log R_d gegen log M_b
    (Mittel und Streuung über Ziehen der Galaxien), Q ≤ 2."""
    x, y = [], []
    for n, p in sparc_galaxien().items():
        if p["Q"] > 2:
            continue
        z = [s for s in open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1") if s.split()[:1] == [n]][0].split()
        Rd = float(z[11])
        Mb = (yd * p["L36"] + 1.33 * p["MHI"]) * 1e9
        x.append(math.log10(Mb)); y.append(math.log10(Rd))
    x, y = np.array(x), np.array(y)
    rs = np.log10(np.sqrt(G * 10**x * M_sonne / a0) / kpc)
    rng = np.random.default_rng(8)
    st = [np.polyfit(x[i], y[i], 1)[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(ziehungen))]
    q = 10**(y - rs)
    return float(np.median(q)), tuple(np.percentile(q, [16, 84])), float(np.polyfit(x, y, 1)[0]), float(np.std(st))


@lru_cache(maxsize=None)
def pubchem_fenster():
    """Breite des Fensters T_Sied/T_Schmelz über alle Elemente mit Werten, schmalstes und weitestes, Rangkorrelation mit der
    Ionisierungsenergie."""
    from scipy.stats import spearmanr
    w = [(s, v["Tb"] / v["Ts"], v["IE"]) for s, v in pubchem().items() if v["Ts"] and v["Tb"] and v["IE"]]
    lo, hi = min(w, key=lambda x: x[1]), max(w, key=lambda x: x[1])
    rho, p = spearmanr([x[1] for x in w], [x[2] for x in w])
    return len(w), (lo[0], lo[1]), (hi[0], hi[1]), float(rho), float(p)


@lru_cache(maxsize=None)
def sparc_gas_gezeit(QT=1.0, weg="vgas", sigma=None, je=False, dichte=1.0):
    """Anteil des Gases unter der Gezeitendichte seiner Scheibe, gelesen wie in einer Scheibe an der Schwelle bei der
    Stabilität Q_T: Die Mitte der Schicht liegt außerhalb x_c·R_d unter ρ_J (kern.toomre_schwelle, gezeit = Q_T²).
    weg "vgas": Gas innerhalb x_c·R_d aus V_gas²·R/G an diesem Ort (Galaxien, deren Kurve dort Punkte hat), "flaeche":
    aus der Flächendichte des Gases (sparc_flaeche, innen bis zum ersten Punkt gleich dem ersten Wert), "hi":
    exponentielle Scheibe aus M_HI und R_HI (Σ = 1 M☉/pc² bei R_HI, die steilere Lösung). Gas = 1,33·M_HI.
    sigma: Streuung des Gases [km/s] als Zahl oder Funktion von lg M_HI [M☉] statt eines festen Q_T. Dann ist Q_T am Ort
    σ/σ₁(x) mit der Streuung σ₁ an Q_T = 1, und x_c liegt je Galaxie, wo σ/v_flach = (σ₁/v_flach)·√gezeit (v_flach aus
    der Tafel, sonst das größte gemessene Tempo).
    dichte: die Scheibe bei gleicher Masse und gleichem v_flach dichter als die Schwelle, Σ₀ = dichte·a₀/(πG), Gas und Sterne
    um denselben Faktor enger (kern.toomre_schwelle mit dichte), die Lage in Einheiten von R_d bleibt.
    Gibt (x_c, Zahl, Median, nach Gasmasse gewichtet), mit sigma statt x_c den Median der x_c. je: dazu je Galaxie
    (lg M_HI [M☉], Anteil)."""
    import kern as K
    from scipy.optimize import brentq
    if sigma is None:
        xc = brentq(lambda x: K.toomre_schwelle(x, dichte=dichte)["gezeit"] / QT**2 - 1, 1.05, 30.0)
    else:
        xs = np.linspace(1.05, 30.0, 600)
        ys = np.array([(lambda t: t["sigma_vf"] * math.sqrt(t["gezeit"]))(K.toomre_schwelle(x, dichte=dichte)) for x in xs])
        alle_xc = []
    tab = {}
    for z in open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
        p = z.split()
        if len(p) >= 18:
            try:
                tab[p[0]] = dict(Rd=float(p[11]), MHI=float(p[13]), RHI=float(p[14]), Vf=float(p[15]))
            except ValueError:
                pass
    a, w = [], []
    for n, t in tab.items():
        if t["MHI"] <= 0 or t["Rd"] <= 0:
            continue
        if sigma is not None:
            vf = t["Vf"] if t["Vf"] > 0 else float(np.max(sparc_kurve(n)["V"]))
            sg = sigma(math.log10(t["MHI"] * 1e9)) if callable(sigma) else sigma
            xc = float(np.interp(-sg / vf, -ys, xs))
        Rc = xc * t["Rd"]
        if weg in ("vgas", "flaeche"):
            k = sparc_kurve(n)
            if not (k["R"][0] <= Rc or sigma is not None) or Rc >= k["R"][-1]:
                continue
            if weg == "vgas":
                innen = max(float(np.interp(Rc, np.r_[0.0, k["R"]], np.r_[0.0, k["Vgas"]])), 0.0)**2 * 1e6 * Rc * K.kpc / K.G / K.M_sun
            else:
                r = np.linspace(0.0, Rc, 400)
                innen = float(np.trapezoid(2 * math.pi * r * np.interp(r, k["R"], sparc_flaeche(n)["Sg"]), r)) * 1e6
            a.append(max(0.0, 1 - innen / (1.33e9 * t["MHI"])))
        else:
            g = lambda Rg: t["MHI"] * 1e9 / (2 * math.pi * Rg**2) * math.exp(-t["RHI"] / Rg) - 1e6
            if t["RHI"] <= 0 or g(t["RHI"] / 2) <= 0:
                continue
            y = Rc / brentq(g, 0.02 * t["RHI"], t["RHI"] / 2)
            a.append((1 + y) * math.exp(-y))
        w.append(t["MHI"])
        if sigma is not None:
            alle_xc.append(xc)
    a, w = np.array(a), np.array(w)
    aus = (xc if sigma is None else float(np.median(alle_xc))), len(a), float(np.median(a)), float(np.sum(a * w) / w.sum())
    return aus + ((np.log10(w * 1e9), a),) if je else aus


def gas_duenn_frueh(punkt, z=None, weg="flaeche", faktor=1.0):
    """Anteil des Gases unter der Gezeitendichte früher: Bei gleicher Masse ist eine Scheibe enger, R_d ∝ H^{−1/3} (an RC100
    gelesen), ihre Flächendichte das (H(z)/H₀)^{2/3}-fache der Schwelle, H(z) im Atmen des Punkts (kern.atmen_punkt). Die
    Streuung des Gases wie heute (hi_streuung_gang) mal faktor, sparc_gas_gezeit auf dem Weg weg. z = None: z der Stapelung
    bei z ≈ 1 (hi_frueh_eins). Gibt (Dichte gegen die Schwelle, Median des Anteils heute, früher, Verhältnis, Zahl früher)."""
    import kern as K
    z = hi_frueh_eins()[0] if z is None else z
    H = K.atmen_punkt(punkt)[4]
    d = (H(1 / (1 + z)) / H(1.0))**(2 / 3)
    f = hi_streuung_gang()[3]
    heute = sparc_gas_gezeit(weg=weg, sigma=f)[2]
    frueh = sparc_gas_gezeit(weg=weg, sigma=lambda x: faktor * f(x), dichte=d)
    return d, heute, frueh[2], frueh[2] / heute, frueh[1]


def sparc_spiralgrenze(rand=False):
    """Welche SPARC-Galaxien hält ihr Gas außerhalb R_d stabil: Eine Scheibe an der Schwelle ist dort stabil, wenn das Gas
    mit mehr als dem größten σ₁/v_flach außerhalb R_d streut (kern.toomre_schwelle bei x = 1). Die Streuung je Galaxie
    aus ihrem M_HI (hi_streuung_gang, rand: zweiter Weg von Hand), v_flach aus der Tafel, sonst das größte gemessene
    Tempo. Gibt (Zahl über der Grenze, davon vom Typ T ≥ 9, davon mit v_flach aus der Tafel, größtes σ/v_flach einer
    Spirale T ≤ 8, ihr Name, Zahl aller, davon T ≥ 9)."""
    import kern as K
    grenze = K.toomre_schwelle(1.0)["sigma_vf"]
    f = hi_streuung_gang(rand=rand)[3]
    zeilen = []
    for n, t in sparc_galaxien().items():
        if t["MHI"] <= 0:
            continue
        vf = t["Vflat"] if t["Vflat"] > 0 else float(np.max(sparc_kurve(n)["V"]))
        zeilen.append((f(math.log10(t["MHI"] * 1e9)) / vf, t["T"], t["Vflat"] > 0, n))
    ueber = [z for z in zeilen if z[0] > grenze]
    spirale = max(z for z in zeilen if z[1] <= 8)
    return (len(ueber), sum(z[1] >= 9 for z in ueber), sum(z[2] for z in ueber), spirale[0], spirale[3], len(zeilen),
            sum(z[1] >= 9 for z in zeilen))


def sparc_flaeche(name, h_rel=0.15, knoten=16, glatt=1e-3, neu=False):
    """Flächendichte der Baryonen einer SPARC-Galaxie an den Punkten der Kurve [M☉/pc²]: Sterne SBdisk (Υ = 1), Gas aus
    V_gas als nicht negative Summe von Dreiecken in R (knoten Knoten bis 1,5·R_max, Krümmung mit glatt gedämpft), deren
    Feld in der Ebene (welt.achsen_feld, senkrecht e^{−|z|/h}, h = h_rel·R_d) V_gas·|V_gas| trifft. Abgelegt unter
    "_sparc_flaeche". Gibt dict(R, Sd, Sg, Vd2 = V_disk² aus SBdisk mit derselben Dicke, Vg2 = getroffenes V_gas²)."""
    import kern as K
    import welt as W
    import huelle_ablage as HA
    from scipy.optimize import nnls
    schluessel = f"{name} h_rel={h_rel} knoten={knoten} glatt={glatt}"
    alle = HA.lesen()
    if not neu and schluessel in alle.get("_sparc_flaeche", {}):
        return {x: np.array(y) for x, y in alle["_sparc_flaeche"][schluessel].items()}
    k, Rd = sparc_kurve(name), sparc_galaxien()[name]["Rd"]
    R, kpc, s = k["R"], K.kpc, K.M_sun / K.pc**2
    T = np.concatenate([[0.0], np.geomspace(R[0] / 2, 1.5 * R[-1], knoten)])
    dreieck = lambda j: (lambda r: s * np.interp(r / kpc, T, np.eye(len(T))[j], right=0.0))
    stern = lambda r: s * np.where(r / kpc <= R[-1], np.interp(r / kpc, R, k["SBdisk"]),
                                   k["SBdisk"][-1] * np.exp(-(r / kpc - R[-1]) / Rd))
    h = h_rel * Rd * kpc
    g = W.achsen_feld([(stern, h)] + [(dreieck(j), h) for j in range(len(T))], None, R * kpc, je_teil=True,
                      schritt=(min(R[0] / 2, Rd / 8, (T[2] - T[1]) / 2) * kpc, h / 8))
    A = np.array(g[1:]).T * R[:, None] * kpc / 1e6
    D2 = np.diff(np.eye(len(T)), 2, axis=0)
    b = k["Vgas"] * np.abs(k["Vgas"])
    c = nnls(np.vstack([A, math.sqrt(glatt) * np.linalg.norm(A) / np.linalg.norm(D2) * D2]),
             np.concatenate([b, np.zeros(len(D2))]))[0]
    erg = dict(R=R, Sd=k["SBdisk"], Sg=np.interp(R, T, c), Vd2=g[0] * R * kpc / 1e6, Vg2=A @ c)
    alle = HA.lesen()
    alle.setdefault("_sparc_flaeche", {})[schluessel] = {x: [float(v) for v in y] for x, y in erg.items()}
    HA.schreiben(alle)
    return erg


def sparc_wirbelfrei(name, yd=YD, a0s=None, h_rel=0.15, fein=1):
    """Wirbelfreie Regel in der Ebene einer SPARC-Galaxie gegen ν(g_N)·g_N: Sterne yd·SBdisk, Gas aus sparc_flaeche,
    Bulge als Kugel aus V_bul (Υ_b = YB), gelöst mit welt.achsen_feld. Gibt je a₀ (a0s [m/s²]) lg(wirbelfrei/(ν·g_N)) an
    den Punkten der Kurve. fein: Gitter fein-mal feiner (Schritt kleiner von R_min/2 und R_d/8, senkrecht h/8,
    Wachstum 1 + 0,03/fein)."""
    import kern as K
    import welt as W
    a0s = tuple(a0s or (K.a0,))
    k, Rd, f = sparc_kurve(name), sparc_galaxien()[name]["Rd"], sparc_flaeche(name, h_rel)
    R, kpc, s = k["R"], K.kpc, K.M_sun / K.pc**2
    stern = lambda r: yd * s * np.where(r / kpc <= R[-1], np.interp(r / kpc, R, k["SBdisk"]),
                                        k["SBdisk"][-1] * np.exp(-(r / kpc - R[-1]) / Rd))
    gas = lambda r: s * np.interp(r / kpc, R, f["Sg"], right=0.0)
    kugel = None
    if k["Vbul"].max() > 0:
        Mb = YB * k["Vbul"]**2 * 1e6 * R * kpc / K.G
        kugel = lambda r: np.interp(r / kpc, np.concatenate([[0.0], R]), np.concatenate([[0.0], Mb]))
    h = h_rel * Rd * kpc
    gN, v, w = W.achsen_feld([(stern, h), (gas, h)], kugel, R * kpc, a0s=a0s, stufe=1 + 0.03 / fein,
                             schritt=(min(R[0] / 2, Rd / 8) * kpc / fein, h / 8 / fein))
    return [np.log10(w[i] / v[i]) for i in range(len(a0s))]


YD_TAFEL = (0.25, 0.35, 0.5, 0.7, 1.0, 1.4)


def sparc_wirbelfrei_tafel(neu=False, melden=None):
    """sparc_wirbelfrei aller 175 Galaxien je Υ_d aus YD_TAFEL und a₀ aus GITTER, abgelegt unter "_sparc_wirbelfrei"
    (je Galaxie, sobald fertig). Gibt Name → Feld (Υ_d, a₀, Punkt) in dex."""
    import huelle_ablage as HA
    alle = HA.lesen()
    tafel = {} if neu else alle.get("_sparc_wirbelfrei", {})
    for n in sparc_galaxien():
        if n in tafel:
            continue
        tafel[n] = [[[round(float(x), 5) for x in a] for a in sparc_wirbelfrei(n, y, tuple(10**GITTER))] for y in YD_TAFEL]
        alle = HA.lesen()
        alle["_sparc_wirbelfrei"] = tafel
        HA.schreiben(alle)
        if melden:
            melden(n)
    return {n: np.array(v) for n, v in tafel.items()}


@lru_cache(maxsize=None)
def sparc_drehung_gang(yd, H0, Q, ziehungen=2000):
    """Gang des Drehanteils λ' = √2·R_d·v_flach/(V₂₀₀·R₂₀₀) der SPARC-Galaxien (Q ≤ 2, wie sparc_drehung) mit der
    baryonischen Masse: (Steigung von lg λ' gegen lg M_b, Streuung der Steigung über Ziehen der Galaxien, Zahl)."""
    import kern as K
    rc = 3 * (H0 * 1e3 / K.Mpc)**2 / (8 * math.pi * K.G)
    lm, ll = [], []
    for p in sparc_galaxien().values():
        if p["Q"] > 2 or p["Vflat"] <= 0 or p["Rd"] <= 0:
            continue
        Mb = (yd * p["L36"] + 1.33 * p["MHI"]) * 1e9
        M = (1 + Q) * Mb * K.M_sun
        R2 = (3 * M / (4 * math.pi * 200 * rc))**(1 / 3)
        lm.append(math.log10(Mb))
        ll.append(math.log10(math.sqrt(2) * p["Rd"] * K.kpc * p["Vflat"] * 1e3 / (math.sqrt(K.G * M / R2) * R2)))
    lm, ll = np.array(lm), np.array(ll)
    rng = np.random.default_rng(1)
    bs = [np.polyfit(lm[i], ll[i], 1)[0] for i in (rng.integers(0, len(lm), len(lm)) for _ in range(ziehungen))]
    return float(np.polyfit(lm, ll, 1)[0]), float(np.std(bs)), len(lm)


def sparc_drehung(yd, H0, Q, profile=False):
    """Drehanteil, den die Scheiben verlangen: J/M = 2R_d·v_flach gegen das Gebiet (1 + Q)·M_b, gelesen wie in den Läufen
    vieler Körper, λ' = √2·R_d·v_flach/(V₂₀₀·R₂₀₀) mit R₂₀₀ bei 200-facher kritischer Dichte (H0 [km/s/Mpc]),
    M_b = yd·L[3.6] + 1,33·M_HI, Galaxien mit Q ≤ 2. Gibt (Zahl, Median, 16 %, 84 %). profile: J/M aus den Profilen
    von Sternen (yd·SBdisk) und Gas (sparc_flaeche) mit der gemessenen Kurve über die Punkte, ∫Σ·V·R²dR/∫Σ·R dR."""
    import kern as K
    rc = 3 * (H0 * 1e3 / K.Mpc)**2 / (8 * math.pi * K.G)
    lam = []
    for z in open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
        p = z.split()
        if len(p) < 18:
            continue
        try:
            L, Rd, MHI, Vf, q = float(p[7]), float(p[11]), float(p[13]), float(p[15]), int(p[17])
        except ValueError:
            continue
        if Vf <= 0 or Rd <= 0 or q > 2:
            continue
        M = (1 + Q) * (yd * L + 1.33 * MHI) * 1e9 * K.M_sun
        R2 = (3 * M / (4 * math.pi * 200 * rc))**(1 / 3)
        j = 2 * Rd * Vf
        if profile:
            k = sparc_kurve(p[0])
            S = yd * k["SBdisk"] + np.interp(k["R"], *[sparc_flaeche(p[0])[x] for x in ("R", "Sg")])
            j = np.trapezoid(S * k["V"] * k["R"]**2, k["R"]) / np.trapezoid(S * k["R"], k["R"])
        lam.append(j * K.kpc * 1e3 / (math.sqrt(2) * math.sqrt(K.G * M / R2) * R2))
    return len(lam), float(np.median(lam)), float(np.percentile(lam, 16)), float(np.percentile(lam, 84))


def sparc_gasbruch(yd, v_unten, v_oben):
    """Median des Gasanteils 1,33·M_HI/M_b der SPARC-Galaxien mit v_unten ≤ v_flach < v_oben."""
    f = [1.33 * p["MHI"] / (yd * p["L36"] + 1.33 * p["MHI"]) for p in sparc_galaxien().values() if p["Vflat"] > 0 and v_unten <= p["Vflat"] < v_oben]
    return float(np.median(f))


@lru_cache(maxsize=None)
def schliesszahlen(anzahl=6):
    """Schlüsse der Zählregeln: die stärksten Sprünge der Zwei-Nukleonen-Bindung (Median über die Kette, gerade N bzw. Z,
    lokales Maximum; AME2020) und die stärksten Abfälle der Ionisierungsenergie IE(Z)/IE(Z+1) (PubChem)."""
    from collections import defaultdict
    B = {(Z, N): b * A / 1000 for Z, N, A, b, gem in ame2020_markiert() if A > 4 and gem}
    kern = {}
    for achse in (0, 1):
        sprung = defaultdict(list)
        for (Z, N), b in B.items():
            k = (Z, N)[1 - achse]
            if k % 2:
                continue
            d = (0, 2) if achse == 0 else (2, 0)
            m, p = (Z - d[0], N - d[1]), (Z + d[0], N + d[1])
            if m in B and p in B:
                sprung[k].append((b - B[m]) - (B[p] - b))
        ks = sorted(k for k in sprung if len(sprung[k]) >= 3 and k >= 6)
        v = [float(np.median(sprung[k])) for k in ks]
        spitzen = [(v[i], ks[i]) for i in range(1, len(ks) - 1) if v[i] > v[i - 1] and v[i] > v[i + 1]]
        kern["NZ"[achse]] = tuple(sorted(k for _, k in sorted(spitzen, reverse=True)[:anzahl]))
    ie = {v["Z"]: v["IE"] for v in pubchem().values() if v["IE"]}
    abfall = sorted(((ie[z] / ie[z + 1], z) for z in ie if z + 1 in ie), reverse=True)[:anzahl]
    atom = tuple(sorted(z for _, z in abfall))
    return kern["N"], kern["Z"], atom, ie[2] / ie[3], ie[86] / ie[87]


@lru_cache(maxsize=None)
def pantheon():
    """Pantheon+SH0ES: (zHD, zHEL, m_b_corr) und Kovarianz (Statistik + Systematik), nur z > 0,01."""
    import gzip
    p = os.path.join(ORT, "pantheon")
    t = np.genfromtxt(os.path.join(p, "Pantheon+SH0ES.dat"), names=True, dtype=None, encoding=None)
    with gzip.open(os.path.join(p, "Pantheon+SH0ES_STAT+SYS.cov.gz"), "rt") as f:
        n = int(f.readline())
        C = np.loadtxt(f).reshape(n, n)
    w = t["zHD"] > 0.01
    return t["zHD"][w], t["zHEL"][w], t["m_b_corr"][w], C[np.ix_(w, w)]


def desi():
    """DESI DR2 BAO: z, Wert, Größe (DV/DM/DH über r_d) und Kovarianz."""
    p = os.path.join(ORT, "desi")
    z, v, g = [], [], []
    for s in open(os.path.join(p, "desi_dr2_bao_mean.txt")):
        if s.strip() and not s.startswith("#"):
            a, b, c = s.split()
            z.append(float(a)); v.append(float(b)); g.append(c)
    return np.array(z), np.array(v), g, np.loadtxt(os.path.join(p, "desi_dr2_bao_cov.txt"))


@lru_cache(maxsize=None)
def abstaende_pruefung(om_modell, or_, genau=False, frei=True):
    """Form von H(z) an Pantheon+ und DESI DR2, Nullpunkte (M, H₀·r_d) frei, flach.
    Gibt χ² des Modells (Ω_m fest), von ΛCDM und w₀w_a (frei) für SN+BAO, die Bestwerte und H₀·r_d des Modells [km/s]."""
    from scipy.integrate import quad
    from scipy.optimize import minimize
    zs, zh, m, C = pantheon()
    Ci = np.linalg.inv(C)
    zb, vb, gb, Cb = desi()
    Cbi = np.linalg.inv(Cb)
    stuetz = np.concatenate([[0.0], np.sort(np.unique(np.concatenate([zs, zb])))])

    def E(om, w0=-1.0, wa=0.0):
        return lambda z: math.sqrt(or_ * (1 + z)**4 + om * (1 + z)**3
                                   + (1 - om - or_) * (1 + z)**(3 * (1 + w0 + wa)) * math.exp(-3 * wa * z / (1 + z)))

    def chi(e):
        if genau:
            stueck = [quad(lambda x: 1 / e(x), a, b)[0] for a, b in zip(stuetz[:-1], stuetz[1:])]
            dm = dict(zip(stuetz, np.concatenate([[0.0], np.cumsum(stueck)])))
        else:
            g = np.union1d(stuetz, np.linspace(0, stuetz[-1], 20001))
            f = 1 / np.vectorize(e)(g)
            I = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
            dm = dict(zip(stuetz, np.interp(stuetz, g, I)))
        r = m - 5 * np.log10((1 + zh) * np.array([dm[x] for x in zs]))
        eins = np.ones_like(r)
        sn = r @ Ci @ r - (eins @ Ci @ r)**2 / (eins @ Ci @ eins)
        th = np.array([{"DM_over_rs": dm[x], "DH_over_rs": 1 / e(x), "DV_over_rs": (x * dm[x]**2 / e(x))**(1 / 3)}[g]
                       for x, g in zip(zb, gb)])
        A, B = th @ Cbi @ th, th @ Cbi @ vb
        return sn + vb @ Cbi @ vb - B * B / A, 299792.458 * A / B

    modell, h0rd = chi(E(om_modell))
    if not frei:
        return dict(modell=float(modell), h0rd=float(h0rd))
    lcdm = minimize(lambda p: chi(E(p[0]))[0], [0.30], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3))
    w0wa = minimize(lambda p: chi(E(*p))[0], [0.30, -0.9, -0.3], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3, maxiter=3000))
    return dict(modell=float(modell), lcdm=float(lcdm.fun), om_lcdm=float(lcdm.x[0]), w0wa=float(w0wa.fun),
                w=tuple(float(x) for x in w0wa.x), h0rd=float(h0rd))


@lru_cache(maxsize=None)
def abstaende_a0(om_or):
    """a₀ aus den Abständen allein: om_or(a₀) gibt (Ω_m, Ω_r) des Formteils; Bestwert und Δχ² = 1-Grenzen [10⁻¹⁰ m/s²]."""
    from scipy.optimize import brentq, minimize_scalar
    f = lambda a: abstaende_pruefung(*om_or(a * 1e-10), frei=False)["modell"]
    r = minimize_scalar(f, bounds=(0.95, 1.3), method="bounded", options=dict(xatol=1e-4))
    return float(r.x), float(brentq(lambda a: f(a) - r.fun - 1, 0.9, r.x, xtol=1e-4)), float(brentq(lambda a: f(a) - r.fun - 1, r.x, 1.4, xtol=1e-4))


@lru_cache(maxsize=None)
def spiegelkerne():
    """Spiegelkerne (N = Z + 1 gegen Z = N + 1, AME2020): Bindungsunterschied = Coulomb-Unterschied a_C·2Z/A^{1/3};
    Zahl der Paare und a_C [MeV] (kleinste Quadrate durch den Ursprung)."""
    B = {(Z, N): b * A / 1000 for Z, N, A, b, gem in ame2020_markiert() if A > 4 and gem}
    x, y = [], []
    for (Z, N), b in B.items():
        if N == Z + 1 and (N, Z) in B:
            x.append(2 * Z / (Z + N)**(1 / 3)); y.append(b - B[(N, Z)])
    x, y = np.array(x), np.array(y)
    return len(x), float(x @ y / (x @ x))


@lru_cache(maxsize=None)
def pantheon_richtung(om, or_, zmin, zmax=3.0, kataloge=False):
    """Richtung in den Resten von Pantheon+ (volle Kovarianz, Formteil mit Ω_m, Ω_r): Dipol als ΔH/H und seine σ,
    Quadrupol (5 Beiwerte) und seine σ; wahlweise mit eigenem Nullpunkt je Teilkatalog."""
    import gzip
    from scipy.stats import chi2, norm
    p = os.path.join(ORT, "pantheon")
    t = np.genfromtxt(os.path.join(p, "Pantheon+SH0ES.dat"), names=True, dtype=None, encoding=None)
    with gzip.open(os.path.join(p, "Pantheon+SH0ES_STAT+SYS.cov.gz"), "rt") as f:
        n = int(f.readline())
        C = np.loadtxt(f).reshape(n, n)
    E = lambda z: np.sqrt(or_ * (1 + z)**4 + om * (1 + z)**3 + 1 - om - or_)
    g = np.linspace(0, 2.5, 30001)
    I = np.concatenate([[0], np.cumsum((1 / E(g[1:]) + 1 / E(g[:-1])) / 2 * np.diff(g))])
    w = (t["zHD"] > zmin) & (t["zHD"] < zmax)
    r = t["m_b_corr"][w] - 5 * np.log10((1 + t["zHEL"][w]) * np.interp(t["zHD"][w], g, I))
    Ci = np.linalg.inv(C[np.ix_(w, w)])
    ra, de = np.radians(t["RA"][w]), np.radians(t["DEC"][w])
    nx, ny, nz = np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)
    ids = t["IDSURVEY"][w]
    if kataloge:
        B = [(ids == i).astype(float) for i in np.unique(ids) if (ids == i).sum() >= 5]
        rest = ~np.any(np.column_stack(B) > 0, 1)
        B = np.column_stack(B + ([rest.astype(float)] if rest.any() else []))
    else:
        B = np.ones((len(r), 1))

    def gls(X):
        A = X.T @ Ci @ X
        q = np.linalg.lstsq(A, X.T @ Ci @ r, rcond=None)[0]
        e = r - X @ q
        return q, e @ Ci @ e
    sig = lambda d, k: float(norm.isf(chi2.sf(d, k) / 2))
    _, c0 = gls(B)
    q1, c1 = gls(np.column_stack([B, nx, ny, nz]))
    _, c2 = gls(np.column_stack([B, nx, ny, nz, nx * ny, nx * nz, ny * nz, nx * nx - ny * ny, 3 * nz * nz - 1]))
    return float(np.linalg.norm(q1[-3:]) * np.log(10) / 5), sig(c0 - c1, 3), sig(c1 - c2, 5), int(w.sum())


@lru_cache(maxsize=None)
@lru_cache(maxsize=None)
def sh0es_leiter(anker="alle"):
    """Die ganze Leiter von SH0ES 2022 (Riess u. a. 2022, Gl. 6: y = L·q mit Kovarianz C, 3492 Gleichungen, 47 Zahlen),
    nach kleinsten Quadraten gelöst. anker: "alle", "N4258", "LMC" oder "MW", dann bleiben die Bedingungen der anderen Anker
    weg (Zeilen 3207 und 3208 Milchstraße, 3213 NGC 4258, 3214 LMC). Gibt (H₀, Fehler, M_B, Fehler)."""
    from astropy.io import fits
    p = os.path.join(ORT, "sh0es")
    y = fits.getdata(os.path.join(p, "ally_shoes_ceph_topantheonwt6.0_112221.fits")).astype(float)
    L = fits.getdata(os.path.join(p, "alll_shoes_ceph_topantheonwt6.0_112221.fits")).astype(float).T
    C = fits.getdata(os.path.join(p, "allc_shoes_ceph_topantheonwt6.0_112221.fits")).astype(float)
    weg = {"alle": (), "N4258": (3207, 3208, 3214), "LMC": (3207, 3208, 3213), "MW": (3213, 3214)}[anker]
    k = np.setdiff1d(np.arange(len(y)), np.array(weg, int))
    Ci = np.linalg.inv(C[np.ix_(k, k)])
    F = L[k].T @ Ci @ L[k]
    q = np.linalg.solve(F, L[k].T @ Ci @ y[k])
    e = np.sqrt(np.diag(np.linalg.inv(F)))
    H0 = 10**(q[46] / 5)
    return float(H0), float(H0 * math.log(10) / 5 * e[46]), float(q[42]), float(e[42])


def pantheon_h0(om, or_, zlo, zhi, w0=-1.0, wa=0.0):
    """H₀ aus Pantheon+ mit den Cepheiden-Eichungen (volle Kovarianz, Form des Atmens mit Ω_m, Ω_r, Dunkle Energie mit
    w(a) = w₀ + w_a(1 − a)): (H₀ [km/s/Mpc], Fehler, Zahl der Supernovae im Fluss, Helligkeit M)."""
    import gzip
    p = os.path.join(ORT, "pantheon")
    t = np.genfromtxt(os.path.join(p, "Pantheon+SH0ES.dat"), names=True, dtype=None, encoding=None)
    with gzip.open(os.path.join(p, "Pantheon+SH0ES_STAT+SYS.cov.gz"), "rt") as f:
        n = int(f.readline())
        C = np.loadtxt(f).reshape(n, n)
    E = lambda z: np.sqrt(or_ * (1 + z)**4 + om * (1 + z)**3 + (1 - om - or_) * _de_gang(z, w0, wa))
    g = np.linspace(0, 2.5, 30001)
    I = np.concatenate([[0], np.cumsum((1 / E(g[1:]) + 1 / E(g[:-1])) / 2 * np.diff(g))])
    cal = t["IS_CALIBRATOR"] == 1
    fl = (~cal) & (t["zHD"] > zlo) & (t["zHD"] < zhi)
    w = cal | fl
    Ci = np.linalg.inv(C[np.ix_(w, w)])
    k = cal[w]
    dl = (1 + t["zHEL"][w]) * np.interp(t["zHD"][w], g, I) * 299792.458
    y = t["m_b_corr"][w] - np.where(k, t["CEPH_DIST"][w], 5 * np.log10(dl) + 25)
    X = np.column_stack([np.ones_like(y), np.where(k, 0.0, -5.0)])
    A = X.T @ Ci @ X
    q = np.linalg.solve(A, X.T @ Ci @ y)
    H = 10**q[1]
    return float(H), float(H * np.log(10) * np.sqrt(np.linalg.inv(A)[1, 1])), int(fl.sum()), float(q[0])


def atmen_praezision():
    """Das Atmen aller Präzisionsdaten im Modell (ω_c = (2π − 1)·ω_b, H₀ aus a₀, Maximum der Dunklen Energie an der
    Sättigung des Kreises), aus der Ablage: (Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀, w_a), H₀ [km/s/Mpc], Fehler von H₀
    aus der Kovarianz von ω_b und a₀, Ω_m, Ω_r)."""
    import kern as K
    Q = K.q_je_radiant()
    e = planck_anpassung(Q, ("spektren", "linse", "bao", "sn", "wachstum", "bbn"), zustand="saettigung", start={"wa": -0.56})
    p, n, C = e["punkt"], e["kovarianz"]["namen"], np.array(e["kovarianz"]["werte"])
    H = lambda a0, wb: K.atmen_heute(a0 * 1e-10, Q=Q, wb=wb, wnu=0.06 / 93.14)[0]
    g = np.zeros(len(n))
    g[n.index("a0")] = (H(p["a0"] + 1e-5, p["wb"]) - H(p["a0"] - 1e-5, p["wb"])) / 2e-5
    g[n.index("wb")] = (H(p["a0"], p["wb"] + 1e-7) - H(p["a0"], p["wb"] - 1e-7)) / 2e-7
    h2 = (p["H0"] / 100)**2
    return ((p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p["w0"], p["wa"]), p["H0"], float(math.sqrt(g @ C @ g)),
            (p["wb"] + p["wc"] + 0.06 / 93.14) / h2, K.dichte_heute()[0] / K.rho_krit100 / h2)


def pantheon_leerraum(om, or_, z_rand, eps, w0=-1.0, wa=0.0):
    """Sprung ε = H_v/H − 1 eines Leerraums um den Beobachter mit Rand bei z_rand, ohne Eichung aus der Form von
    Pantheon+ gelesen (pantheon_form mit leerraum): (bester Sprung, Fehler bei Δχ² = 1, Δχ² des Sprungs eps gegen den
    besten)."""
    from scipy.optimize import brentq, minimize_scalar
    chi = lambda e: pantheon_form(om, or_, w0, wa, leerraum=(z_rand, e))
    b = minimize_scalar(chi, bounds=(-0.3, 0.3), method="bounded", options=dict(xatol=1e-6))
    lo = brentq(lambda e: chi(e) - b.fun - 1, -0.3, b.x, xtol=1e-6)
    hi = brentq(lambda e: chi(e) - b.fun - 1, b.x, 0.3, xtol=1e-6)
    return float(b.x), (hi - lo) / 2, chi(eps) - b.fun


@lru_cache(maxsize=None)
@lru_cache(maxsize=None)
def planck_lite():
    """Planck 2018, plik_lite v22: gebinnte Spektren C_b [µK²] für TT (215 Bins, ℓ = 30–2508), TE und EE (je 199 Bins,
    ℓ = 30–1996) mit ihrer Kovarianz (Cholesky-Faktor), Grenzen und Gewichte der Bins; dazu TT unter ℓ = 30 als zwei
    Bins (ℓ = 2–15, 16–29, nahezu gaußisch, Prince und Dunkley 2019)."""
    import scipy.linalg
    from scipy.io import FortranFile
    d = os.path.join(ORT, "planck")
    _, X, _ = np.genfromtxt(os.path.join(d, "cl_cmb_plik_v22.dat"), unpack=True)
    lese = lambda n, t=float: np.loadtxt(os.path.join(d, n)).astype(t)
    cov = FortranFile(os.path.join(d, "c_matrix_plik_v22.dat"), "r").read_reals(dtype=float).reshape((613, 613))
    L = scipy.linalg.cholesky(np.tril(cov) + np.tril(cov, -1).T, lower=True)
    _, Xt, st = np.genfromtxt(os.path.join(d, "CTT_bin_low_ell_2018.dat"), unpack=True)
    tief = (Xt, st, lese("blmin_low_ell.dat", int), lese("blmax_low_ell.dat", int), lese("bweight_low_ell.dat"))
    return X, L, lese("blmin.dat", int), lese("blmax.dat", int), lese("bweight.dat"), tief


def planck_rest(Dtt, Dte, Dee, eichung=1.0, lmin=2):
    """Geweißter Rest (Daten − Modell)/σ der Spektren D_ℓ (ab ℓ = lmin) gegen plik_lite, gebinnt wie bei Planck;
    die Summe der Quadrate ist χ²."""
    import scipy.linalg
    X, L, blmin, blmax, bw, (Xt, st, tmin, tmax, tw) = planck_lite()
    ls = np.arange(len(Dtt)) + lmin
    fak = ls * (ls + 1) / (2 * math.pi)
    binne = lambda C, ab, b0, b1, w, n: [np.sum(C[b0[i] + ab - lmin:b1[i] + ab + 1 - lmin] * w[b0[i]:b1[i] + 1]) for i in range(n)]
    m = np.array(binne(Dtt / fak, 30, blmin, blmax, bw, 215) + binne(Dte / fak, 30, blmin, blmax, bw, 199)
                 + binne(Dee / fak, 30, blmin, blmax, bw, 199))
    mt = np.array(binne(Dtt / fak, 2, tmin, tmax, tw, 2))
    return np.concatenate([(Xt - mt / eichung**2) / st, scipy.linalg.solve_triangular(L, X - m / eichung**2, lower=True)])


def planck_chi2(Dtt, Dte, Dee, eichung=1.0, lmin=2):
    """χ² der Spektren gegen plik_lite mit den zwei Bins unter ℓ = 30 (ohne Vorgaben)."""
    r = planck_rest(Dtt, Dte, Dee, eichung, lmin)
    return float(r @ r)


@lru_cache(maxsize=None)
@lru_cache(maxsize=None)
def planck_karte(fwhm_grad=1.0, nside=256):
    """Karte der Hintergrundstrahlung (Planck SMICA, planck_karte/smica_t_512.npz): Temperatur [µK] und Maske, RING,
    galaktisch. Außerhalb der Maske null gesetzt, mit einem Gaußstrahl fwhm_grad geglättet, durch die ebenso geglättete
    Maske geteilt, frei bleibt, wo die geglättete Maske über 0,9 liegt; dann auf nside gemittelt. Gibt (T, frei, nside)."""
    import healpy as hp
    z = np.load(os.path.join(ORT, "planck_karte", "smica_t_512.npz"))
    T, M = z["T_uK"].astype(float), z["maske"].astype(float)
    T = (T - np.sum(T * M) / M.sum()) * M
    if fwhm_grad:
        f = math.radians(fwhm_grad)
        T, M = hp.smoothing(T, fwhm=f, iter=1), hp.smoothing(M, fwhm=f, iter=1)
    frei = M > 0.9
    T = np.where(frei, T / np.maximum(M, 1e-3), 0.0)
    if nside != 512:
        T, frei = hp.ud_grade(T, nside), hp.ud_grade(frei.astype(float), nside) > 0.999
    return T, frei, nside


@lru_cache(maxsize=None)
def planck_alm(lmax=30):
    """Große Formen der Karte (Planck SMICA, ganzer Himmel ohne Maske, Monopol und Dipol entfernt): a_ℓm bis lmax in der
    Ordnung von healpy, µK, galaktisch."""
    import healpy as hp
    z = np.load(os.path.join(ORT, "planck_karte", "smica_t_512.npz"))
    T = hp.remove_dipole(z["T_uK"].astype(float), fitval=False)
    return hp.map2alm(T, lmax=lmax, iter=3)


def kreise_gegenueber(alphas_grad, fwhm_grad=1.0, richtungen=16, nside=256, zufall=None, schritt_grad=0.25, karte=None, achsen=None):
    """Kreispaare der geschlossenen Zelle: Liegt das Abbild des Beobachters im Abstand L < 2χ* in Richtung n, sieht er
    dieselben Punkte der letzten Streufläche auf dem Kreis um n und auf dem Kreis um −n, beide mit dem Radius α,
    cos α = L/(2χ*), Punkt für Punkt d₁ = cos α·n + sin α·e(φ) ↔ d₂ = −cos α·n + sin α·e(φ). Für jede Richtung n (Mitten
    des Gitters nside = richtungen auf der halben Kugel) und jedes α: S = 2Σ T₁T₂/Σ(T₁² + T₂²) über die Punkte, die beide
    frei liegen (planck_karte). zufall: Saat der Nullprobe, der zweite Kreis liegt dann um eine zufällige Richtung mit
    zufälliger Drehung. karte: (T, frei, nside) statt Planck, etwa ein gerechneter Himmel. achsen: vorgegebene Richtungen
    (Einheitsvektoren) statt des Gitters. Gibt (S[Richtung, α], Richtungen
    als Einheitsvektoren, Anteil freier Punkte)."""
    import healpy as hp
    T, frei, ns = planck_karte(fwhm_grad, nside) if karte is None else karte
    if achsen is None:
        n = np.array(hp.pix2vec(richtungen, np.arange(hp.nside2npix(richtungen)))).T
        n = n[(n[:, 2] > 1e-9) | ((abs(n[:, 2]) <= 1e-9) & (n[:, 1] > 0))]
    else:
        n = np.array(achsen, dtype=float)
    rng = np.random.default_rng(zufall) if zufall is not None else None
    m = -n if rng is None else rng.normal(size=n.shape)
    m /= np.linalg.norm(m, axis=1)[:, None]

    def basis(v):
        a = np.where(abs(v[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
        e1 = np.cross(v, a)
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        return e1, np.cross(v, e1)
    e1, e2 = basis(n)
    f1, f2 = (e1, e2) if rng is None else basis(m)
    dreh = np.zeros(len(n)) if rng is None else rng.uniform(0, 2 * math.pi, len(n))
    S = np.zeros((len(n), len(alphas_grad)))
    anteil = np.zeros_like(S)
    for j, ag in enumerate(alphas_grad):
        a = math.radians(ag)
        k = max(32, int(2 * math.pi * math.sin(a) / math.radians(schritt_grad)))
        phi = np.linspace(0, 2 * math.pi, k, endpoint=False)
        c, s_ = math.cos(a), math.sin(a)
        d1 = c * n[:, None, :] + s_ * (np.cos(phi)[None, :, None] * e1[:, None, :] + np.sin(phi)[None, :, None] * e2[:, None, :])
        p2 = phi[None, :] + dreh[:, None]
        d2 = c * m[:, None, :] + s_ * (np.cos(p2)[:, :, None] * f1[:, None, :] + np.sin(p2)[:, :, None] * f2[:, None, :])
        i1 = hp.vec2pix(ns, d1[..., 0], d1[..., 1], d1[..., 2])
        i2 = hp.vec2pix(ns, d2[..., 0], d2[..., 1], d2[..., 2])
        w = frei[i1] & frei[i2]
        t1, t2 = np.where(w, T[i1], 0.0), np.where(w, T[i2], 0.0)
        S[:, j] = 2 * np.sum(t1 * t2, axis=1) / np.maximum(np.sum(t1 * t1 + t2 * t2, axis=1), 1e-30)
        anteil[:, j] = w.mean(axis=1)
    return S, n, anteil


def kreise_dreier(alphas_grad, fwhm_grad=0.5, richtungen=64, nside=512, zufall=None, schritt_grad=1.0, karte=None):
    """Die Würfelzelle hat drei senkrechte Achsen, jede trägt ein Kreispaar vom selben Radius: S je Achse aus
    kreise_gegenueber (Gitter richtungen), dann für jede Achse n₁ und jede Drehung ψ um sie (Schritt schritt_grad) die
    Achsen n₂ ⊥ n₁ und n₃ = n₁ × n₂, S₃ = (S(n₁) + S(n₂) + S(n₃))/3 am nächsten Gitterpunkt. Gibt je α das größte S₃
    und die drei Achsen dazu (galaktisch, Einheitsvektoren). karte wie bei kreise_gegenueber."""
    import healpy as hp
    S, n, anteil = kreise_gegenueber(alphas_grad, fwhm_grad, richtungen, nside, zufall, 0.25, karte)
    S = np.where(anteil > 0.5, S, -1.0)
    voll = np.full((hp.nside2npix(richtungen), len(alphas_grad)), -1.0)
    for v, sv in zip(n, S):
        for w in (v, -v):
            voll[hp.vec2pix(richtungen, *w)] = sv
    aus = []
    psi = np.radians(np.arange(0.0, 90.0, schritt_grad))
    for j in range(len(alphas_grad)):
        bestes, achsen = -1.0, None
        for v in n:
            a = np.array([0.0, 0.0, 1.0]) if abs(v[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            e1 = np.cross(v, a)
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(v, e1)
            n2 = np.cos(psi)[:, None] * e1 + np.sin(psi)[:, None] * e2
            n3 = np.cross(v, n2)
            s3 = (voll[hp.vec2pix(richtungen, *v), j] + voll[hp.vec2pix(richtungen, n2[:, 0], n2[:, 1], n2[:, 2]), j]
                  + voll[hp.vec2pix(richtungen, n3[:, 0], n3[:, 1], n3[:, 2]), j]) / 3
            k = int(np.argmax(s3))
            if s3[k] > bestes:
                bestes, achsen = float(s3[k]), (v, n2[k], n3[k])
        aus.append((bestes, achsen))
    return aus


def planck_linse():
    """Planck PR4 (NPIPE), Linsenpotential über die Spektren gemittelt (Carron, Mirmelstein, Lewis 2022): 9 Bandleistungen
    [L(L + 1)]²C_L^φφ/2π für L = 8–400, Kovarianz (Cholesky-Faktor), Fenster je Bin und Fenster der linearen Korrektur
    (L = 2–2500) mit ihrem Bezugswert."""
    import scipy.linalg
    d = os.path.join(ORT, "planck_linse")
    stamm = os.path.join(d, "pp_consext8_npipe_smicaed_TiPi_jTP_pre30T_kfilt_rdn0cov_PS1")
    X = np.loadtxt(stamm + "_bandpowers.dat")[:, 4]
    L = scipy.linalg.cholesky(np.loadtxt(stamm + "_CMBmarged_cov.dat"), lower=True)
    fenster = np.zeros((2, 9, 2499))
    for b in range(9):
        for k, n in enumerate(("win", "kor")):
            w = np.loadtxt(os.path.join(d, f"{n}{b + 1}.dat"))
            fenster[k, b, w[:, 0].astype(int) - 2] = w[:, 1]
    bezug = np.loadtxt(stamm + "_CMBmarged_lensing_fiducial_correction.dat")[:, 1]
    return X, L, fenster, bezug


def planck_linse_rest(Cpp):
    """Geweißter Rest der Linsenbandleistungen gegen das Potential Cpp = [L(L + 1)]²C_L^φφ/2π (L = 2–2500)."""
    import scipy.linalg
    X, L, fenster, bezug = planck_linse()
    m = fenster[0] @ Cpp + fenster[1] @ Cpp - bezug
    return scipy.linalg.solve_triangular(L, X - m, lower=True)


@lru_cache(maxsize=None)
def wachstum():
    """Wachstum f·σ₈ aus der vollen Form der Galaxienverteilung: SDSS DR12 (z = 0,38; 0,51), eBOSS DR16 LRG (0,698) und
    Quasare (1,48), je über die Abstände derselben Messung gemittelt (Teilblock der Kovarianz). Gibt (z fallend, Wert,
    Cholesky-Faktor)."""
    d = os.path.join(ORT, "wachstum")
    z, v, bl = [], [], []
    for n in ("LRG", "QSO"):
        stamm = os.path.join(d, f"sdss_DR16_BAOplus_{n}_FSBAO_DMDHfs8")
        zeilen = [s.split() for s in open(stamm + ".dat") if s.strip()]
        i = [k for k, s in enumerate(zeilen) if s[2] == "f_sigma8"]
        z += [float(zeilen[k][0]) for k in i]
        v += [float(zeilen[k][1]) for k in i]
        bl.append(np.loadtxt(stamm + "_covtot.txt")[np.ix_(i, i)])
    import scipy.linalg
    C = scipy.linalg.block_diag(*bl)
    o = np.argsort(z)[::-1]
    return tuple(np.array(z)[o]), np.array(v)[o], np.linalg.cholesky(C[np.ix_(o, o)])


def bbn_rest(wb):
    """Geweißter Rest von D/H und Y_p gegen die Urknall-Nukleosynthese bei ω_b (PRIMAT-Tafel aus CAMB, N_eff = 3,044;
    Fehler der Messung und der Kernraten zusammen)."""
    import kern as K
    from camb import bbn
    p = bbn.get_predictor()
    dh = (K.DH_GEMESSEN[0] - p.DH(wb, 0.0)) / math.hypot(K.DH_GEMESSEN[1], p.get("sig(D/H)", wb, 0.0))
    yp = (K.YP_GEMESSEN[0] - p.Y_p(wb, 0.0)) / math.hypot(K.YP_GEMESSEN[1], p.get("sig(Yp^BBN)", wb, 0.0))
    return np.array([dh, yp])


@lru_cache(maxsize=None)
def zwerge():
    """Druckgestützte Zwerge mit gemessener Sichtlinienstreuung (Local Volume Database): Satelliten der Milchstraße und von
    M31 ab 20 kpc Abstand zum Wirt und isolierte Zwerge; ohne die drehenden, gasreichen (HI über 10^7,5 M☉). Je Zwerg dict
    mit Name, Wirt, L_V [L☉], M_HI [M☉], R_e [kpc] (kreisgemittelt), Abstand zum Wirt [kpc] (isoliert None), σ und Fehler
    [km/s], [Fe/H] der Sterne (spektroskopisch, sonst wie in der Tafel; None, wenn nicht gemessen)."""
    import csv
    import kern as K
    aus = []
    for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))):
        if z["confirmed_galaxy"] != "1" or not z["vlos_sigma"] or z["vlos_sigma_ul"] or z["host"] not in ("mw", "m_031", ""):
            continue
        hi = 10**float(z["mass_HI"]) if z["mass_HI"] else 0.0
        dh = float(z["distance_host"]) if z["host"] else None
        if hi > 10**7.5 or (dh is not None and dh < 20):
            continue
        e = float(z["ellipticity"]) if z["ellipticity"] else 0.0
        zahl = lambda k: float(z[k]) if z[k] else None
        aus.append(dict(name=z["key"], wirt=z["host"], L=10**(-0.4 * (float(z["M_V"]) - K.M_V_SONNE)), M_HI=hi,
                        R_e=float(z["rhalf"]) / 60 * math.pi / 180 * float(z["distance"]) * math.sqrt(1 - e),
                        d_wirt=dh, sigma=float(z["vlos_sigma"]),
                        fehler=(float(z["vlos_sigma_em"]) + float(z["vlos_sigma_ep"])) / 2,
                        ra=float(z["ra"]), dec=float(z["dec"]), abstand=float(z["distance"]), pmra=zahl("pmra"),
                        pmdec=zahl("pmdec"), vlos=zahl("vlos_systemic"),
                        feh=zahl("metallicity_spectroscopic") if z["metallicity_spectroscopic"] else zahl("metallicity")))
    return tuple(aus)


def fornax_haufen(yv=2.0, yv_haufen=2.0):
    """Fornax und seine Kugelsternhaufen (Local Volume Database): Fornax mit Masse Υ_V·L_V [M☉] (Υ_V = yv, ohne Gas),
    kreisgemitteltem R_e [kpc], Abstand [kpc] und Streuung σ [km/s]; je Haufen (Wirt fornax_1) Name, Masse yv_haufen·L_V
    [M☉], projizierter Abstand zur Mitte von Fornax [kpc] und Radius r_½ [pc]."""
    import csv
    import kern as K
    zeilen = list(csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))))
    f = next(z for z in zeilen if z["key"] == "fornax_1")
    D = float(f["distance"])
    ra0, de0 = math.radians(float(f["ra"])), math.radians(float(f["dec"]))
    licht = lambda z: 10**(-0.4 * (float(z["M_V"]) - K.M_V_SONNE))
    haufen = []
    for z in zeilen:
        if z["host"] != "fornax_1" or not z["M_V"]:
            continue
        ra, de = math.radians(float(z["ra"])), math.radians(float(z["dec"]))
        w = math.acos(min(1.0, math.sin(de0) * math.sin(de) + math.cos(de0) * math.cos(de) * math.cos(ra - ra0)))
        haufen.append(dict(name=z["key"], m=yv_haufen * licht(z), r_p=D * w,
                           r_halb=float(z["rhalf"]) / 60 * math.pi / 180 * float(z["distance"]) * 1e3))
    e = float(f["ellipticity"]) if f["ellipticity"] else 0.0
    return dict(M=yv * licht(f), R_e=float(f["rhalf"]) / 60 * math.pi / 180 * D * math.sqrt(1 - e), D=D,
                sigma=float(f["vlos_sigma"]), haufen=sorted(haufen, key=lambda h: h["name"]))


def zwerge_form(yv=2.0):
    """Form der Zwerge in der Gezeit (Local Volume Database): Elliptizität der Satelliten gegen lg(ρ/ρ_J), ρ die mittlere
    Dichte in 4/3·R_e mit der Regel (Masse aus Licht mit Υ_V = yv und Gas), ρ_J = 3·d_eff·ν·g_N/(4πG·D) die Gezeitendichte
    des Wirts samt Umlauf mit der Regel, d_eff = 3 − χ. Gibt (Zahl der Satelliten, Rangkorrelation, p, Median der
    Elliptizität der Satelliten, der isolierten Zwerge, Zahl der isolierten)."""
    import csv
    import kern as K
    from scipy.stats import spearmanr
    ell = {z["key"]: float(z["ellipticity"]) for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv")))
           if z["ellipticity"]}
    e, x, iso = [], [], []
    for z in zwerge():
        if z["name"] not in ell:
            continue
        if not z["wirt"]:
            iso.append(ell[z["name"]])
            continue
        M = (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun
        rs = 4 / 3 * z["R_e"] * K.kpc
        rho = K.verstaerkung(K.G * M / 2 / rs**2) * M / 2 / (4 * math.pi / 3 * rs**3)
        mw = {"mw": K.M_B_MILCHSTRASSE, "m_031": K.M_B_M31}[z["wirt"]] * K.M_sun
        D = z["d_wirt"] * K.kpc
        gN = K.G * mw / D**2
        e.append(ell[z["name"]])
        x.append(math.log10(rho / (3 * (3 - K.schalter(gN)) * K.verstaerkung(gN) * gN / D / (4 * math.pi * K.G))))
    r = spearmanr(e, x)
    return len(e), float(r[0]), float(r[1]), float(np.median(e)), float(np.median(iso)), len(iso)


@lru_cache(maxsize=None)
def m31_paar(quer=False, bereich=True):
    """Die Bereiche von Milchstraße und M31 aus der Local Volume Database (v1.1.1): Abstand der Schwerpunkte ihrer
    Baryonen (bereich_schwerpunkt) [kpc] und Annäherung [km/s], quer: dazu die Bewegung quer zur Verbindung [km/s] und die
    Fehler der Eigenbewegung von M31 (HST, van der Marel u. a. 2012) je Richtung als Geschwindigkeit [km/s].
    bereich=False: Milchstraße und M31 selbst, Abstand zum Zentrum der Milchstraße und Annäherung im Ruhesystem der
    Milchstraße (−velocity_gsr der Zeile m_031), quer aus Abstand, Sichtliniengeschwindigkeit und Eigenbewegung (Sonne
    nach kern.SONNE_GALAKTISCH)."""
    import csv
    import kern as K
    z = next(z for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))) if z["key"] == "m_031")
    if not (quer or bereich):
        return float(z["distance_gc"]), -float(z["velocity_gsr"])
    if bereich:
        (x0, v0, _), (x1, v1, _) = bereich_schwerpunkt("mw"), bereich_schwerpunkt("m_031")
        x, v = x1 - x0, v1 - v0
    else:
        x, v = _galaktisch(z)
    r = float(np.linalg.norm(x))
    v_r = x @ v / r
    if not quer:
        return r, -float(v_r)
    je = 4.740_47 * float(z["distance"])                          # km/s je mas/a im Abstand [kpc]
    return (r if bereich else float(z["distance_gc"]), -float(v_r) if bereich else -float(z["velocity_gsr"]),
            float(math.sqrt(v @ v - v_r**2)), je * float(z["pmra_em"]), je * float(z["pmdec_em"]))


def _galaktisch(z):
    """Ort [kpc] und Geschwindigkeit [km/s] einer Zeile der Local Volume Database im Rahmen der Milchstraße aus Abstand,
    Sichtliniengeschwindigkeit und Eigenbewegung (Sonne nach kern.SONNE_GALAKTISCH)."""
    import astropy.coordinates as ko
    import astropy.units as u
    import kern as K
    s = K.SONNE_GALAKTISCH
    rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                               galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
    c = ko.SkyCoord(ra=float(z["ra"]) * u.deg, dec=float(z["dec"]) * u.deg, distance=float(z["distance"]) * u.kpc,
                    pm_ra_cosdec=float(z["pmra"]) * u.mas / u.yr, pm_dec=float(z["pmdec"]) * u.mas / u.yr,
                    radial_velocity=float(z["vlos_systemic"]) * u.km / u.s).transform_to(rahmen)
    return (np.array([c.x.value, c.y.value, c.z.value]),
            np.array([c.v_x.to(u.km / u.s).value, c.v_y.to(u.km / u.s).value, c.v_z.to(u.km / u.s).value]))


@lru_cache(maxsize=None)
def bereich_schwerpunkt(wirt="mw"):
    """Schwerpunkt der Baryonen eines Bereichs (mw oder m_031) im Rahmen der Milchstraße: Ort [kpc] und Geschwindigkeit
    [km/s] aus dem Wirt (die Milchstraße ruhend im Zentrum, M31 mit Eigenbewegung) und allen Gliedern seines Bereichs
    (bereich_baryonen) mit Abstand, Sichtliniengeschwindigkeit und Eigenbewegung, je Sterne und 1,33·M_HI, dazu der
    Anteil der Baryonen des Bereichs darin."""
    import csv
    import kern as K
    zeilen = {z["key"]: z for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv")))}
    f = lambda x: float(x) if x not in ("", "nan") else math.nan
    masse = lambda z: sum(c * 10**f(z[k]) for c, k in ((1.0, "mass_stellar"), (1.33, "mass_HI")) if math.isfinite(f(z[k])))
    if wirt == "mw":
        teile = [(K.M_B_MILCHSTRASSE, np.zeros(3), np.zeros(3))]
    else:
        teile = [(K.M_B_M31,) + _galaktisch(zeilen[wirt])]
    teile += [(masse(zeilen[k]),) + _galaktisch(zeilen[k]) for k in bereich_baryonen(wirt, glieder=True)
              if all(math.isfinite(f(zeilen[k][s])) for s in ("pmra", "pmdec", "vlos_systemic", "distance"))]
    m = np.array([t[0] for t in teile])
    return (sum(t[0] * t[1] for t in teile) / m.sum(), sum(t[0] * t[2] for t in teile) / m.sum(),
            float(m.sum() / bereich_baryonen(wirt)[0]))


def lokaler_fluss(von=0.6, bis=3.0):
    """Zwerge um die Lokale Gruppe aus der Local Volume Database (v1.1.1, Tafeln dwarf_local_field und
    dwarf_local_field_distant) mit Abstand von bis bis [Mpc] vom Schwerpunkt der Gruppe (distance_lg) und der
    Geschwindigkeit im Rahmen der Gruppe (velocity_lg): (Abstand [Mpc], radiale Geschwindigkeit [km/s])."""
    import csv
    zeilen = list(csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))))
    f = lambda x: float(x) if x not in ("", "nan") else math.nan
    a = np.array([(f(z["distance_lg"]) / 1e3, f(z["velocity_lg"])) for z in zeilen
                  if z["table"] in ("dwarf_local_field", "dwarf_local_field_distant")])
    a = a[np.isfinite(a).all(axis=1) & (a[:, 0] >= von) & (a[:, 0] <= bis)]
    return a[:, 0], a[:, 1]


def bereich_baryonen(wirt="mw", glieder=False):
    """Baryonen eines Bereichs [M☉]: der Wirt (kern.M_B_MILCHSTRASSE, M_B_M31, sonst aus der Tafel misc, Sterne und
    1,33·M_HI) und alle Zwerge und Galaxien der Local
    Volume Database (Tafeln dwarf_* und misc), deren Kette der Wirte zu ihm führt, je Sterne und 1,33·M_HI. Gibt (Summe,
    Wirt, Satelliten, Zahl), mit glieder die Namen der Glieder."""
    import csv
    import kern as K
    zeilen = [z for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv")))
              if z["table"].startswith("dwarf") or z["table"] == "misc"]
    wirte = {z["key"]: z["host"] for z in zeilen}

    def kette(k):
        while k:
            k = wirte.get(k, "")
            if k == wirt:
                return True
        return False
    if glieder:
        return [z["key"] for z in zeilen if kette(z["key"])]
    f = lambda x: 10**float(x) if x not in ("", "nan") else 0.0
    sat = [f(z["mass_stellar"]) + 1.33 * f(z["mass_HI"]) for z in zeilen if kette(z["key"])]
    z = next((z for z in zeilen if z["key"] == wirt), None)
    eigen = {"mw": K.M_B_MILCHSTRASSE, "m_031": K.M_B_M31}.get(wirt, f(z["mass_stellar"]) + 1.33 * f(z["mass_HI"]) if z else 0.0)
    return eigen + sum(sat), eigen, sum(sat), len(sat)


def wirt_satelliten(wirt, R_max=400.0):
    """Satelliten um einen Wirt der Local Volume Database (Zwerge mit Sichtliniengeschwindigkeit, deren Wirt er ist) im
    projizierten Abstand bis R_max [kpc] (Winkel zum Wirt mal dessen Abstand): Streuung σ [km/s] als quadratisches
    Mittel von v − v_S um den Schwerpunkt der Baryonen v_S von Wirt und Begleitern mit Geschwindigkeit bis R_max (bei
    m_031 mit Eigenbewegung, satelliten_schar), Zahl, Mittel der Zwerge gegen den Wirt und gegen v_S [km/s], Median des
    projizierten Abstands [kpc] und die Baryonen [M☉] des Wirts und aller seiner Begleiter (Tafeln dwarf_* und misc)
    bis R_max."""
    import csv
    zeilen = list(csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))))
    w = next(z for z in zeilen if z["key"] == wirt)
    f = lambda x: float(x) if x not in ("", "nan") else None
    ra0, de0, D = math.radians(f(w["ra"])), math.radians(f(w["dec"])), f(w["distance"])
    R = lambda z: D * math.acos(min(1.0, math.sin(de0) * math.sin(math.radians(f(z["dec"])))
                                    + math.cos(de0) * math.cos(math.radians(f(z["dec"])))
                                    * math.cos(math.radians(f(z["ra"])) - ra0)))
    mit = [z for z in zeilen if z["host"] == wirt and (z["table"].startswith("dwarf") or z["table"] == "misc")
           and R(z) <= R_max]
    masse = lambda z: (10**f(z["mass_stellar"]) if f(z["mass_stellar"]) else 0.0) + 1.33 * (
        10**f(z["mass_HI"]) if f(z["mass_HI"]) else 0.0)
    M_b = bereich_baryonen(wirt)[1] + sum(masse(z) for z in mit)
    if wirt == "m_031":
        zeile = {z["key"]: z for z in zeilen}
        s = [(x["v"], x["R"], masse(zeile[x["name"]]), x["tafel"]) for x in satelliten_schar("m_031", begleiter=True)
             if x["R"] <= R_max]
    else:
        s = [(f(z["vlos_systemic"]) - f(w["vlos_systemic"]), R(z), masse(z), z["table"]) for z in mit
             if f(z["vlos_systemic"]) is not None]
    v_S = sum(m * v for v, _, m, _ in s) / (bereich_baryonen(wirt)[1] + sum(m for _, _, m, _ in s))
    dv = np.array([v for v, _, _, t in s if t.startswith("dwarf")])
    Rp = [r for _, r, _, t in s if t.startswith("dwarf")]
    return dict(sigma=float(np.sqrt(np.mean(np.square(dv - v_S)))), zahl=len(dv), mittel=float(dv.mean()),
                mittel_S=float(dv.mean() - v_S), fehler=float(dv.std(ddof=1) / math.sqrt(len(dv))),
                R_median=float(np.median(Rp)), M_b=M_b)


@lru_cache(maxsize=None)
def satellit(name="lmc", geschichte=False):
    """Ein Satellit der Milchstraße aus der Local Volume Database (v1.1.1), auch die gasreichen (LMC, SMC): Ort [m] und
    Geschwindigkeit [m/s] aus Abstand, Sichtliniengeschwindigkeit und Eigenbewegung im Rahmen der Milchstraße (Sonne nach
    kern.SONNE_GALAKTISCH), gedreht so, dass z längs des Drehimpulses der Scheibe zeigt (x, −y, −z des galaktozentrischen
    Rahmens), und seine Baryonen [M☉], Sterne und 1,33·M_HI. geschichte=True: das Gas, das er einmal trug, wenn es mehr ist
    als heute, wie isolierte Zwerge gleicher Helligkeit (zwerge_gas_beziehung, Gedächtnis der Bindung). Gibt dict(x, v,
    M_b)."""
    import csv
    import astropy.coordinates as ko
    import astropy.units as u
    import kern as K
    z = next(z for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))) if z["key"] == name)
    s = K.SONNE_GALAKTISCH
    rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                               galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
    c = ko.SkyCoord(ra=float(z["ra"]) * u.deg, dec=float(z["dec"]) * u.deg, distance=float(z["distance"]) * u.kpc,
                    pm_ra_cosdec=float(z["pmra"]) * u.mas / u.yr, pm_dec=float(z["pmdec"]) * u.mas / u.yr,
                    radial_velocity=float(z["vlos_systemic"]) * u.km / u.s).transform_to(rahmen)
    x = (c.x.to(u.m).value, -c.y.to(u.m).value, -c.z.to(u.m).value)
    v = (c.v_x.to(u.m / u.s).value, -c.v_y.to(u.m / u.s).value, -c.v_z.to(u.m / u.s).value)
    hi = 10**float(z["mass_HI"]) if z["mass_HI"] else 0.0
    if geschichte:
        a, b = zwerge_gas_beziehung()[:2]
        L = 10**(-0.4 * (float(z["M_V"]) - K.M_V_SONNE))
        hi = max(hi, L * 10**(a + b * math.log10(L)))
    M_b = 10**float(z["mass_stellar"]) + 1.33 * hi
    return dict(x=tuple(float(a) for a in x), v=tuple(float(a) for a in v), M_b=M_b)


def zwerge_gas_beziehung():
    """Gas der isolierten Zwerge mit gemessenem HI (Local Volume Database, ohne obere Grenzen): Gerade
    lg(M_HI/L_V) = a + b·lg L_V. Gibt (a, b, Zahl, Streuung [dex], kleinstes und größtes lg L_V)."""
    import csv
    import kern as K
    x, y = [], []
    for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))):
        if z["confirmed_galaxy"] != "1" or not z["mass_HI"] or not z["M_V"] or z["mass_HI_ul"] or z["host"]:
            continue
        lg_l = -0.4 * (float(z["M_V"]) - K.M_V_SONNE)
        x.append(lg_l)
        y.append(float(z["mass_HI"]) - lg_l)
    b, a = np.polyfit(x, y, 1)
    return float(a), float(b), len(x), float(np.std(np.array(y) - a - b * np.array(x))), min(x), max(x)


def _sterne_metall(mu):
    """Mittlere Metallizität der Sterne im geschlossenen Kasten, in Einheiten der Ausbeute p, beim Gasanteil μ."""
    return 1 + mu * math.log(mu) / (1 - mu)


@lru_cache(maxsize=None)
def zwerge_ausbeute(yv=2.0):
    """Ausbeute p (Einheiten der Sonne) des geschlossenen Kastens, geeicht an den gasreichen isolierten Zwergen: [Fe/H] der
    Sterne gegen den heutigen Gasanteil μ = 1,33·M_HI/(1,33·M_HI + Υ_V·L_V). Gibt (lg p als Median, Streuung, Zahl)."""
    w = [z["feh"] - math.log10(_sterne_metall(1.33 * z["M_HI"] / (1.33 * z["M_HI"] + yv * z["L"])))
         for z in zwerge_gruppe("gasreich") if z["feh"] is not None]
    return float(np.median(w)), float(np.std(w)), len(w)


def zwerg_gas(z, geschichte=False, yv=2.0):
    """M_HI [M☉] des Zwergs z: heute oder das Gas, das er einmal trug, wenn es mehr ist als heute. geschichte=True: das Gas
    isolierter Zwerge gleicher Helligkeit (zwerge_gas_beziehung); geschichte="metall": das Gas beim Gasanteil μ, den die
    Metallizität seiner Sterne im geschlossenen Kasten liest (Ausbeute aus zwerge_ausbeute), M_Gas = Υ_V·L_V·μ/(1 − μ)."""
    from scipy.optimize import brentq
    if not geschichte:
        return z["M_HI"]
    if geschichte == "metall":
        x = 10**(z["feh"] - zwerge_ausbeute(yv)[0])
        mu = 0.0 if x >= 1 else brentq(lambda m: _sterne_metall(m) - x, 1e-12, 1 - 1e-12)
        return max(z["M_HI"], yv * z["L"] * mu / (1 - mu) / 1.33)
    a, b = zwerge_gas_beziehung()[:2]
    return max(z["M_HI"], z["L"] * 10**(a + b * math.log10(z["L"])))


def zwerg_streuung(z, yv, mit_feld=True, geschichte=False, tief=False, einfall=None):
    """σ [km/s] des Zwergs z nach der Regel auf das Gesamtfeld (Plummer, Masse Υ_V·L_V + 1,33·M_HI), im Newtonfeld des Wirts
    als Punktmasse; mit_feld=False: allein; geschichte: das Gas nach zwerg_gas liest nur die Rückkehr (Gedächtnis der
    Bindung), die Neigung am Ort liest das heutige; tief: allein ganz im Tiefbereich, σ⁴ = (4/81)·GMa₀ mit der Masse der
    Rückkehr, unabhängig von Profil und Radius; einfall = (c, Q, H₀): die nicht übergebene Mode in der Form nach dem Einfall
    statt der Regel (kern.streuung_im_feld)."""
    import kern as K
    M = (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun
    if einfall is not None:
        return K.streuung_im_feld(M, z["R_e"] * K.kpc, einfall=einfall) / 1e3
    M_B = (yv * z["L"] + 1.33 * zwerg_gas(z, geschichte, yv)) * K.M_sun
    if tief:
        return K.streuung_kugel(M_B) / 1e3
    mw = {"mw": K.M_B_MILCHSTRASSE, "m_031": K.M_B_M31}.get(z["wirt"])
    g_e = K.G * mw * K.M_sun / (z["d_wirt"] * K.kpc)**2 if (mw and mit_feld) else 0.0
    return K.streuung_im_feld(M, z["R_e"] * K.kpc, g_e, M_B=M_B) / 1e3


@lru_cache(maxsize=None)
def xcop():
    """X-COP (Eckert u. a. 2022): je Haufen R [kpc] und eingeschlossene Masse [M☉], gesamt (hydrostatisch) und Gas, mit
    unterer und oberer Grenze: {Name: dict(R, M, M_lo, M_hi, Mgas, Mgas_lo, Mgas_hi)}."""
    roh = {}
    for z in open(os.path.join(ORT, "linsen_haufen", "xcop_profiles.txt")):
        if z.startswith("#") or not z.strip():
            continue
        s = z.split()
        roh.setdefault(s[0], []).append([float(x) for x in s[3:]])
    aus = {}
    for n, w in roh.items():
        w = np.array(sorted({tuple(x) for x in w}))                 # doppelte Zeilen der Datei nur einmal
        R, dm, dlo, dhi, g, glo, ghi = w.T
        aus[n] = dict(R=R, M=dm + g, M_lo=dlo + glo, M_hi=dhi + ghi, Mgas=g, Mgas_lo=glo, Mgas_hi=ghi)
    return aus


@lru_cache(maxsize=None)
def kids_isoliert():
    """KiDS-1000, isolierte Linsen (Brouwer u. a. 2021): je Bin der Sternmasse M_* [M☉], R [kpc], g_obs mit 16/84-%-Band
    [m/s²]; dazu die baryonische Masse wie in der Messung, M_*·(1 + f_kalt), f_kalt = 10^(−0,69·lg M_* + 6,63)."""
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "kids1000_rar_isoliert.txt"))
    aus = []
    for b in range(1, 5):
        w = d[d[:, 0] == b]
        m = w[0, 1]
        aus.append(dict(M_stern=m, M_b=m * (1 + 10**(-0.69 * math.log10(m) + 6.63)), R=w[:, 2] * 1e3,
                        g=w[:, 3], g_lo=w[:, 4], g_hi=w[:, 5]))
    return aus


@lru_cache(maxsize=None)
def kids_typ(typ):
    """KiDS-1000, isolierte Linsen eines Typs ("blau", "rot" nach Farbe, "spaet", "frueh" nach Sérsic-Index), Brouwer u. a. 2021,
    Abb. 8: g_bar [m/s²] (Sterne und kaltes Gas nach Boselli 2014, für alle Typen gleich), g_obs = 4G·ΔΣ/(1 + K) und Fehler."""
    import kern as K
    datei = {"blau": "farbe_blau", "rot": "farbe_rot", "spaet": "sersic_spaet", "frueh": "sersic_frueh"}[typ]
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "kids1000_rar_" + datei + ".txt"))
    f = 4 * K.G * K.M_sun / K.pc**2 / d[:, 4]
    ok = d[:, 1] > 0
    return d[ok, 0], f[ok] * d[ok, 1], f[ok] * d[ok, 3]


def kids_typ_reste(typ, lg_von=-20.0, lg_bis=0.0, g_e=0.0, bezug=None):
    """Rest lg(g_obs/Regel(g_bar)) der Linsen eines Typs in den Bins mit lg g_bar zwischen lg_von und lg_bis, im gleichförmigen
    äußeren Feld g_e [m/s²]. bezug: statt gegen die Regel gegen die Linsen dieses Typs in denselben Bins, ohne die Form der
    Regel, mit den Fehlern beider. Gibt (Zahl, gewichtetes Mittel, Fehler, χ² gegen null)."""
    import kern as K
    gb, go, ge = kids_typ(typ)
    w = (np.log10(gb) >= lg_von) & (np.log10(gb) <= lg_bis)
    if bezug is None:
        r = np.log10(go[w] / np.array([K.regel_monopol(x, g_e) for x in gb[w]]))
        e = ge[w] / go[w] / math.log(10)
    else:
        gb_, gob, geb = kids_typ(bezug)
        assert np.allclose(gb_, gb)
        r = np.log10(go[w] / gob[w])
        e = np.hypot(ge[w] / go[w], geb[w] / gob[w]) / math.log(10)
    g = 1 / e**2
    return int(w.sum()), float(np.sum(g * r) / g.sum()), float(1 / math.sqrt(g.sum())), float(np.sum((r / e)**2))


def kids_typ_gedaechtnis(typ, lg_von=-20.0, lg_bis=0.0):
    """Gedächtnis der Bindung an den Linsen eines Typs: a = g_bar + (ν(f·g_bar) − 1)·f·g_bar, f = einmal gebundene durch
    heutige Baryonen (gleich verteilt). Gibt (f, Fehler, χ² mit f, χ² mit f = 1)."""
    from scipy.optimize import minimize_scalar
    import kern as K
    gb, go, ge = kids_typ(typ)
    w = (np.log10(gb) >= lg_von) & (np.log10(gb) <= lg_bis)
    gb, go, e = gb[w], go[w], ge[w] / go[w] / math.log(10)
    chi = lambda f: float(np.sum((np.log10(go / np.array([x + (K.verstaerkung(f * x) - 1) * f * x for x in gb])) / e)**2))
    r = minimize_scalar(lambda lf: chi(10**lf), bounds=(-1, 1.5), method="bounded", options=dict(xatol=1e-6))
    h = 1e-3
    krumm = (chi(10**(r.x + h)) - 2 * chi(10**r.x) + chi(10**(r.x - h))) / h**2
    return float(10**r.x), float(10**r.x * math.log(10) * math.sqrt(2 / krumm)), float(r.fun), chi(1.0)


def kids_typ_lcdm(typ):
    """ΛCDM an denselben Punkten: MICE2, gestapelt wie bei KiDS (linsen_haufen/mice_rar_*.txt), je Massenbin eine Kurve in
    lg g_bar; am Punkt das Mittel der Kurven, die ihn überdecken. Gibt (Zahl, mittlerer Rest lg(g_obs/ΛCDM), χ²)."""
    datei = {"blau": "farbe_blau", "rot": "farbe_rot", "spaet": "sersic_spaet", "frueh": "sersic_frueh"}[typ]
    m = np.loadtxt(os.path.join(ORT, "linsen_haufen", "mice_rar_" + datei + ".txt"))
    kurven = [k[np.argsort(k[:, 0])] for k in np.split(m, len(m) // 15)]
    gb, go, ge = kids_typ(typ)
    r, e = [], []
    for x, y, f in zip(np.log10(gb), go, ge):
        v = [np.interp(x, k[:, 0], k[:, 1]) for k in kurven if k[0, 0] <= x <= k[-1, 0]]
        if v:
            r.append(math.log10(y) - np.mean(v))
            e.append(f / y / math.log(10))
    r, e = np.array(r), np.array(e)
    return len(r), float(r.mean()), float(np.sum((r / e)**2))


def mistele_reste():
    """Kreisgeschwindigkeiten aus der entprojizierten Linse um isolierte Galaxien (KiDS-1000, Mistele u. a. 2024, 4 Bins × 5
    Radien, 34 bis 115 kpc) gegen die Regel auf die Baryonen der Messung als Punktmasse: (χ² mit Fehler aus Statistik und
    System, Zahl, mittlerer Rest lg(v/v_Regel) je Bin)."""
    import kern as K
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "mistele_gg_wl.txt"))
    chi, mittel = 0.0, []
    for b in range(1, 5):
        w = d[d[:, 0] == b]
        R = 10**w[:, 3] * K.kpc
        gN = K.G * 10**w[0, 2] * K.M_sun / R**2
        vp = np.sqrt(np.array([K.verstaerkung(x) * x for x in gN]) * R) / 1e3
        chi += float(np.sum(((w[:, 4] - vp) / np.hypot(w[:, 5], w[:, 6]))**2))
        mittel.append(float(np.mean(np.log10(w[:, 4] / vp))))
    return chi, len(d), mittel


def _kids_vorhersage(b, g_e=0.0, deckel=False, a0_=None, umgerechnet=False):
    """Regel um den Bin b der isolierten Linsen (Punktmasse M_b) an seinen Radien: direkt die Beschleunigung, oder mit
    umgerechnet wie die Messung über ΔΣ (kern.linse_umgerechnet); deckel: jenseits der Naht R_Q auf (1 + Q)·g_N begrenzt."""
    import kern as K
    a0_ = K.a0 if a0_ is None else a0_
    Q = K.q_je_radiant()
    R = np.geomspace(0.01 * K.kpc, 1e6 * K.kpc, 6000) if umgerechnet else b["R"] * K.kpc
    gN = K.G * b["M_b"] * K.M_sun / R**2
    p = np.array([K.regel_monopol(x, g_e, a0_) for x in gN])
    if deckel:
        p = np.where(R > K.naht_radius(b["M_b"] * K.M_sun, Q, a0_), np.minimum(p, (1 + Q) * gN), p)
    return K.linse_umgerechnet(R, p, b["R"] * K.kpc) if umgerechnet else p


def kids_gang(umgerechnet=False):
    """Rest lg(g_obs/Regel) der isolierten Linsen (KiDS-1000, vier Massenbins, Punktmasse M_b), je Radius über die Bins gemittelt:
    (Radien [kpc], mittlerer Rest [dex]); umgerechnet: die Regel über ΔΣ wie die Messung."""
    reste = [np.log10(b["g"] / _kids_vorhersage(b, umgerechnet=umgerechnet)) for b in kids_isoliert()]
    return kids_isoliert()[0]["R"], np.mean(reste, axis=0)


def udg_neigung():
    """AGC 114905: Neigung aus der Scheibe der Sterne (b/a, Eigendicke q0, cos²i = ((b/a)² − q0²)/(1 − q0²)), kleinste und
    größte über b/a ± Fehler und q0 im Bereich, dazu die Neigung, die die Regel verlangt: sin i = sin i_Gas·v_obs/v_Regel am
    äußeren Punkt (v_obs mit i_Gas gelesen). Gibt (i bei b/a und q0 = 0,2, kleinste, größte, verlangte) in Grad."""
    import kern as K
    z = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_agc114905_sterne2024.txt"))
    ba, e, q_lo, q_hi, i_gas = z[1], z[2], z[3], z[4], z[5]
    i = lambda b, q: math.degrees(math.acos(math.sqrt((b * b - q * q) / (1 - q * q))))
    werte = [i(b, q) for b in (ba - e, ba + e) for q in (q_lo, q_hi)]
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_mancera_pina2022.txt"))
    z0 = d[d[:, 0] == z[0]][0]
    R = z0[2] * K.kpc
    gN = (z0[6] * 1e3)**2 / R
    v_regel = math.sqrt(K.verstaerkung(gN) * gN * R) / 1e3
    soll = math.degrees(K.neigung_verlangt(math.sin(math.radians(i_gas)) * z0[3], v_regel))
    return i(ba, 0.5 * (q_lo + q_hi)), min(werte), max(werte), soll


def udg_reste(einfall_c=None, H0=68.5, mitatmen=None):
    """Gasreiche ultradiffuse Galaxien (Mancera Piña u. a. 2022, 6, äußerer Punkt): lg(v_obs/v_Regel) und lg(v_obs/v_Newton) je
    Galaxie, die Regel auf das Newtonsche Feld der Baryonen am selben Radius. einfall_c: statt der Regel die nicht übergebene
    Mode Q·M_b in der Form nach dem Einfall (NFW mit Konzentration c, M₂₀₀ = (1 + Q)·M_b, 200·ρ_krit bei H0, Q = 2π − 1) zum
    Newtonschen Feld. mitatmen = H [1/s]: der Bereich atmet noch mit, die Summe der Verluste s² = a_N/a₀ + H·τ_D. Gibt (Zahl,
    Mittel gegen die Regel oder die Mode, Mittel gegen Newton, kleinstes und größtes a_N/a₀)."""
    import kern as K
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_mancera_pina2022.txt"))
    R = d[:, 2] * K.kpc
    gN = (d[:, 6] * 1e3)**2 / R
    vr = np.sqrt(np.array([K.verstaerkung(x) * x for x in gN]) * R) / 1e3
    if mitatmen is not None:
        s = np.sqrt(gN / K.a0 + mitatmen * K.rueckkehrzeit())
        vr = np.sqrt(gN * R / -np.expm1(-s)) / 1e3
    if einfall_c is not None:
        Q, rc = K.q_je_radiant(), 3 * (H0 * 1e3 / K.Mpc)**2 / (8 * math.pi * K.G)
        m = lambda x: math.log(1 + x) - x / (1 + x)
        Mb = 10**d[:, 1] * K.M_sun
        R200 = (3 * (1 + Q) * Mb / (4 * math.pi * 200 * rc))**(1 / 3)
        Mm = np.array([Q * mb * m(r * einfall_c / r2) / m(einfall_c) for mb, r, r2 in zip(Mb, R, R200)])
        vr = np.sqrt(gN * R + K.G * Mm / R) / 1e3
    return len(d), float(np.mean(np.log10(d[:, 3] / vr))), float(np.mean(np.log10(d[:, 3] / d[:, 6]))), \
        float(gN.min() / K.a0), float(gN.max() / K.a0)


def udg_breite(lg_von=8.9, lg_bis=9.3, D_von=60.0, D_bis=110.0):
    """Gasreiche ultradiffuse Galaxien (Mancera Piña u. a. 2022) in ALFALFA α.100 (Code 1 und 2), unabhängig vom
    Interferometer: je Galaxie die größte Neigung, unter der sie die Regel trifft, sin i ≤ W50/(2·v_Regel) (W50 enthält
    auch die Streuung des Gases), v_Regel am äußeren Punkt (kern.neigung_verlangt), und der Anteil der Quellen gleicher
    HI-Masse (lg M_HI von lg_von bis lg_bis) und gleichen Abstands (D_von bis D_bis [Mpc]), die schmaler sind, W50 je
    v_flach aus dem Gas allein, v_flach⁴ = G·1,33·M_HI·a₀. Gibt (größte Neigungen [°], Zahl der Quellen, Anteile)."""
    import kern as K
    a = alfalfa(code=None, cz_max=20000.0)
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_mancera_pina2022.txt"))
    x = lambda i: a["W50"][i] / K.v_flach(1.33 * 10**a["lgM"][i] * K.M_sun)
    w = (a["D"] > D_von) & (a["D"] < D_bis) & (a["lgM"] > lg_von) & (a["lgM"] < lg_bis)
    neigung, anteil = [], []
    for z in d:
        i = int(np.where(a["AGC"] == int(z[0]))[0][0])
        R = z[2] * K.kpc
        v_regel = K.bahntempo_regel((z[6] * 1e3)**2 * R / K.G, R)
        neigung.append(math.degrees(K.neigung_verlangt(a["W50"][i] * 1e3 / 2, v_regel)))
        anteil.append(float(np.mean(x(w) < x(i))))
    return neigung, int(w.sum()), anteil


def udg_verdichtung():
    """Gasreiche ultradiffuse Galaxien (Mancera Piña u. a. 2022): Verdichtung ihrer Baryonen innerhalb des äußeren Punkts
    gegen ihr Mittel heute, eingeschlossen nach Newton v_N²·R/G (kern.verdichtung_baryonen). Gibt (kleinste, größte)."""
    import kern as K
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_mancera_pina2022.txt"))
    v = [K.verdichtung_baryonen((z[6] * 1e3)**2 * z[2] * K.kpc / K.G, z[2] * K.kpc) for z in d]
    return float(min(v)), float(max(v))


def sparc_gasreich_aussen(a_von, a_bis, neigung=45.0):
    """Gasreiche Galaxien aus SPARC (Qualität 1 und 2, jede Neigung, 1,33·M_HI ≥ Υ_d·L[3,6] mit Υ_d aus sparc_eichung),
    deren äußerer Punkt bei a_N zwischen a_von·a₀ und a_bis·a₀ liegt: Rest lg(V/V_Regel) dort, mit den Baryonen der
    Massenmodelle, wie bei den ultradiffusen Galaxien. Gibt (Zahl, Mittel, Fehler des Mittels, kleinster Rest, Zahl,
    Mittel und Fehler des Mittels der Galaxien unter der Neigung [°])."""
    import kern as K
    yd = sparc_eichung(K.kpc)
    r, inc = [], []
    for n, p in sparc_galaxien().items():
        if p["Q"] > 2 or 1.33 * p["MHI"] < yd * p["L36"]:
            continue
        k = sparc_kurve(n)
        vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + YB * k["Vbul"]**2
        m = np.where((k["V"] > 0) & (vb2 > 0))[0]
        if len(m) < 3:
            continue
        j, R = m[-1], k["R"][m[-1]] * K.kpc
        M = vb2[j] * 1e6 * R / K.G
        if not a_von <= K.G * M / R**2 / K.a0 <= a_bis:
            continue
        r.append(math.log10(k["V"][j] / (K.bahntempo_regel(M, R) / 1e3)))
        inc.append(p["Inc"])
    r, unten = np.array(r), np.array(inc) < neigung
    f = lambda x: float(x.std() / math.sqrt(len(x)))
    return len(r), float(r.mean()), f(r), float(r.min()), int(unten.sum()), float(r[unten].mean()), f(r[unten])


def kids_chi2(g_e=0.0, deckel=False, a0_=None, umgerechnet=False):
    """χ² der Linse um isolierte Galaxien (KiDS-1000, 60 Werte) gegen die Regel auf die Faltung der Baryonen (Punktmasse),
    im gleichförmigen äußeren Feld g_e [m/s²] oder, mit deckel, jenseits der Naht R_Q auf 1 + Q begrenzt (Q = 2π − 1);
    umgerechnet: die Regel über ΔΣ wie die Messung. Fehler je Seite aus dem 16/84-%-Band in lg g. Gibt (χ², mittlerer Rest
    [dex] je Bin innen und außen)."""
    chi, rest = 0.0, []
    for b in kids_isoliert():
        r = np.log10(b["g"] / _kids_vorhersage(b, g_e, deckel, a0_, umgerechnet))
        e = np.where(r > 0, np.log10(b["g"] / b["g_lo"]), np.log10(b["g_hi"] / b["g"]))
        chi += float(np.sum((r / e)**2))
        rest.append((float(r[:5].mean()), float(r[-5:].mean())))
    return chi, rest


def zwerg_feld(z, yv=2.0):
    """Eigenes Newtonsches Feld eines Zwergs am entprojizierten Halblichtradius (4/3)·R_e [m/s²], Sterne Υ_V·L_V und Gas
    1,33·M_HI, die Hälfte der Masse innen."""
    import kern as K
    return K.G * (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun / 2 / (4 / 3 * z["R_e"] * K.kpc)**2


def zwerge_gruppe(gruppe, yv=2.0):
    """Zwerge einer Gruppe: "isoliert", davon "gasreich" (M_HI ≥ 0,05·L_V) und "gasarm", "innen" (Satelliten, eigenes Feld am
    Halblichtradius überwiegt das des Wirts), "außen" (klassische Satelliten über 10⁵ L☉, Wirtsfeld überwiegt), "klassisch"
    (alle Satelliten über 10⁵ L☉), "ultraschwach" (Satelliten bis 10⁵ L☉)."""
    import kern as K
    g_i = lambda z: zwerg_feld(z, yv)
    def g_e(z):
        return K.G * {"mw": K.M_B_MILCHSTRASSE, "m_031": K.M_B_M31}[z["wirt"]] * K.M_sun / (z["d_wirt"] * K.kpc)**2
    wahl = {"isoliert": lambda z: not z["wirt"],
            "gasreich": lambda z: not z["wirt"] and z["M_HI"] >= 0.05 * z["L"],
            "gasarm": lambda z: not z["wirt"] and z["M_HI"] < 0.05 * z["L"],
            "innen": lambda z: z["wirt"] and g_e(z) < g_i(z),
            "außen": lambda z: z["wirt"] and z["L"] > 1e5 and g_e(z) >= g_i(z),
            "klassisch": lambda z: z["wirt"] and z["L"] > 1e5,
            "ultraschwach": lambda z: z["wirt"] and z["L"] <= 1e5}[gruppe]
    return [z for z in zwerge() if wahl(z)]


def zwerge_gas_faktor(yv=2.0):
    """Baryonen je Sternmasse der gasreichen isolierten Zwerge, (Υ_V·L_V + 1,33·M_HI)/(Υ_V·L_V)."""
    return sorted((yv * z["L"] + 1.33 * z["M_HI"]) / (yv * z["L"]) for z in zwerge_gruppe("gasreich"))


def zwerge_yv(gruppe, geschichte=False):
    """Υ_V, bei dem der mittlere Rest der Gruppe (Regel ohne Wirtsfeld) verschwindet; geschichte wie in zwerg_streuung."""
    from scipy.optimize import brentq
    G = zwerge_gruppe(gruppe)
    mittel = lambda lg: np.mean([math.log10(z["sigma"] / zwerg_streuung(z, 10**lg, False, geschichte)) for z in G])
    return 10**brentq(mittel, -1, 3)


def zwerge_reste(gruppe, mit_feld, yv=2.0, streuung=0.1, geschichte=False, tief=False, einfall=None):
    """Rest lg(σ_gemessen/σ_Regel) der Zwerge einer Gruppe: "isoliert", "innen" (Satelliten, eigenes Feld am Halblichtradius
    überwiegt das des Wirts), "außen" (klassische Satelliten über 10⁵ L☉, Wirtsfeld überwiegt), "ultraschwach" (Satelliten
    bis 10⁵ L☉), auch "gasreich", "gasarm", "klassisch"; geschichte und tief wie in zwerg_streuung. Gibt (Zahl, Mittel,
    Fehler des Mittels, χ² mit eigener Streuung)."""
    import kern as K
    G = [z for z in zwerge_gruppe(gruppe, yv) if geschichte != "metall" or z["feh"] is not None]
    r = np.array([math.log10(z["sigma"] / zwerg_streuung(z, yv, mit_feld, geschichte, tief, einfall)) for z in G])
    e = np.array([z["fehler"] / z["sigma"] / math.log(10) for z in G])
    return len(G), float(r.mean()), float(r.std() / math.sqrt(len(G))), float(np.sum(r**2 / (e**2 + streuung**2)))


@lru_cache(maxsize=None)
def kugelhaufen():
    """Sternhaufen mit gemessener Sichtlinienstreuung (Local Volume Database, ohne obere Grenzen): je Haufen dict mit Name,
    Wirt, Abstand zum Wirt [kpc], L_V [L☉], R_e [m] (kreisgemittelt), σ und Fehler [km/s], lg(g_N/a₀) am Halblichtradius
    der Kugel (4/3·R_e) für Υ_V = 1."""
    import csv
    import kern as K
    aus = []
    for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))):
        if (z["confirmed_star_cluster"] != "1" or not z["vlos_sigma"] or z["vlos_sigma_ul"] or not z["vlos_sigma_em"]
                or not z["M_V"] or not z["rhalf"]):
            continue
        e = float(z["ellipticity"]) if z["ellipticity"] else 0.0
        L = 10**(-0.4 * (float(z["M_V"]) - K.M_V_SONNE))
        R = float(z["rhalf"]) / 60 * math.pi / 180 * float(z["distance"]) * math.sqrt(1 - e) * K.kpc
        aus.append(dict(name=z["key"], wirt=z["host"], d_wirt=float(z["distance_host"] or 0), L=L, R_e=R,
                        sigma=float(z["vlos_sigma"]), fehler=(float(z["vlos_sigma_em"]) + float(z["vlos_sigma_ep"])) / 2,
                        lg_g=math.log10(K.G * L * K.M_sun / (4 / 3 * R)**2 / K.a0)))
    return tuple(aus)


@lru_cache(maxsize=None)
def kugelhaufen_reste(newton=False, auswahl="alle", streuung=0.1, zusatz=0.0):
    """Dünne Sternhaufen gegen dichte: Υ_V geeicht an den dichten (g_N > 3·a₀, dort ν < 1,22), dann Rest lg(σ_gemessen/σ)
    der dünnen (g_N < a₀) mit der Regel oder Newton (kern.streuung_im_feld, ohne äußeres Feld). auswahl: "alle", "gut"
    (Fehler der Streuung unter 30 %), "halo" (Milchstraße, über 10 kpc vom Zentrum). zusatz [km/s]: fester Teil der Streuung,
    quadratisch addiert, wie ihn die Bahnen der Doppelsterne geben, die Einzeltrichter sind. Gibt (Υ_V, Zahl dicht, Zahl
    dünn, gewichtetes Mittel, Fehler, χ² mit eigener Streuung)."""
    import kern as K
    from scipy.optimize import brentq
    a = 1e-40 if newton else K.a0
    H = kugelhaufen()
    rest = lambda z, yv: math.log10(z["sigma"] / math.hypot(K.streuung_im_feld(yv * z["L"] * K.M_sun, z["R_e"], 0.0, a) / 1e3,
                                                             zusatz))
    gewicht = lambda z: 1 / ((z["fehler"] / z["sigma"] / math.log(10))**2 + streuung**2)
    dicht = [z for z in H if z["lg_g"] > math.log10(3)]
    yv = brentq(lambda y: sum(gewicht(z) * rest(z, y) for z in dicht), 0.05, 20.0)
    duenn = [z for z in H if z["lg_g"] < 0 and (auswahl == "alle" or (auswahl == "gut" and z["fehler"] < 0.3 * z["sigma"])
                                                or (auswahl == "halo" and z["wirt"] == "mw" and z["d_wirt"] > 10))]
    r = np.array([rest(z, yv) for z in duenn])
    w = np.array([gewicht(z) for z in duenn])
    return yv, len(dicht), len(duenn), float(np.sum(w * r) / np.sum(w)), float(1 / math.sqrt(np.sum(w))), float(np.sum(w * r**2))


@lru_cache(maxsize=None)
def kugelhaufen_gezeit():
    """Dünne Sternhaufen der Milchstraße (g_N < a₀): Rangkorrelation des Rests nach Newton mit lg g_N und mit der Füllung des
    Jacobi-Radius r_h/r_J (heutiger Abstand, flache Kurve, r_J = ℓ/2^{1/3}, Masse mit Υ_V der Regel). Gibt (Zahl,
    Korrelation mit lg g_N, mit lg r_h/r_J, kleinste und größte Füllung)."""
    import kern as K
    from scipy.stats import spearmanr
    yr, yn = kugelhaufen_reste()[0], kugelhaufen_reste(True)[0]
    v = K.SONNE_GALAKTISCH["v"][1] * 1e3
    lg_g, lg_t, r = [], [], []
    for z in kugelhaufen():
        if z["lg_g"] >= 0 or z["wirt"] != "mw":
            continue
        r_j = K.glaettungslaenge(yr * z["L"] * K.M_sun, K.gezeit_flach(v, z["d_wirt"] * K.kpc)) / 2**(1 / 3)
        lg_g.append(z["lg_g"])
        lg_t.append(math.log10(4 / 3 * z["R_e"] / r_j))
        r.append(math.log10(z["sigma"] / (K.streuung_im_feld(yn * z["L"] * K.M_sun, z["R_e"], 0.0, 1e-40) / 1e3)))
    return len(r), float(spearmanr(lg_g, r)[0]), float(spearmanr(lg_t, r)[0]), 10**min(lg_t), 10**max(lg_t)


@lru_cache(maxsize=None)
def stufe_zwerge(yv=2.0, weg="radius"):
    """Sternhaufen und Zwerge gegen ihre Schwelle (Local Volume Database, Masse Υ_V·L_V, bei den Zwergen mit 1,33·M_HI):
    R_h/r* mit R_h = 4/3·R_e, dem Halblichtradius der Kugel. weg "radius": aus R_h und r* = √(GM/a₀), "feld": aus g_N am
    Halblichtradius, R_h/r* = √(a₀/g_N), bei den Sternhaufen aus dem lg g_N der Tafel. Gibt (Median der Sternhaufen, Zahl
    innerhalb r*, Zahl, Median der Zwerge, Zahl außerhalb r*, Zahl)."""
    import kern as K
    if weg == "radius":
        h = [4 / 3 * z["R_e"] / K.schwellenradius(yv * z["L"] * K.M_sun) for z in kugelhaufen()]
        w = [4 / 3 * z["R_e"] * K.kpc / K.schwellenradius((yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun) for z in zwerge()]
    else:
        h = [10**(-z["lg_g"] / 2) / math.sqrt(yv) for z in kugelhaufen()]
        w = [math.sqrt(K.a0 / (K.G * (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun / (4 / 3 * z["R_e"] * K.kpc)**2))
             for z in zwerge()]
    h, w = np.array(h), np.array(w)
    return float(np.median(h)), int((h < 1).sum()), len(h), float(np.median(w)), int((w > 1).sum()), len(w)


@lru_cache(maxsize=None)
def satelliten_gezeit(weg="wende", geschichte=False, yv=2.0):
    """Gezeit am Perizentrum der Satelliten der Milchstraße mit Eigenbewegung: r_J im Perizentrum, beide tief,
    r_J = r_peri·(m/M)^{1/4}/√2 mit M = kern.M_B_MILCHSTRASSE und m = Υ_V·L_V + 1,33·M_HI (geschichte: mit dem Gas, das er
    einmal trug, zwerg_gas), gegen den Halblichtradius R_e, dazu der Rest lg(σ/σ_Regel) (zwerg_streuung). weg "wende":
    Perizentrum aus den Wendepunkten (satelliten_bahnen), "bahn": aus der Bahn, 3 Gyr rückwärts gerechnet (solve_ivp)
    im Feld der Punktmasse nach der Regel. Gibt [(r_J/R_e, Rest, Name, L_V)] nach dem Rest geordnet."""
    import kern as K
    if weg == "wende":
        peri = {n: p for n, (r, p, a) in satelliten_bahnen().items()}
    else:
        import astropy.coordinates as ko
        import astropy.units as u
        from scipy.integrate import solve_ivp
        s = K.SONNE_GALAKTISCH
        rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                                   galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
        GM = K.G * K.M_B_MILCHSTRASSE * K.M_sun

        def rate(t, y):
            r = np.linalg.norm(y[:3])
            return np.concatenate([y[3:], -K.verstaerkung(GM / r**2) * GM / r**3 * y[:3]])
        peri = {}
        for z in zwerge():
            if z["wirt"] != "mw" or z["pmra"] is None or z["vlos"] is None:
                continue
            c = ko.SkyCoord(ra=z["ra"] * u.deg, dec=z["dec"] * u.deg, distance=z["abstand"] * u.kpc,
                            pm_ra_cosdec=z["pmra"] * u.mas / u.yr, pm_dec=z["pmdec"] * u.mas / u.yr,
                            radial_velocity=z["vlos"] * u.km / u.s).transform_to(rahmen)
            y0 = [c.x.to(u.m).value, c.y.to(u.m).value, c.z.to(u.m).value,
                  c.v_x.to(u.m / u.s).value, c.v_y.to(u.m / u.s).value, c.v_z.to(u.m / u.s).value]
            sol = solve_ivp(rate, (0, -3 * K.Gyr), y0, rtol=1e-9, atol=1.0, max_step=0.002 * K.Gyr)
            peri[z["name"]] = float(np.min(np.linalg.norm(sol.y[:3], axis=0))) / K.kpc
    aus = []
    for z in zwerge():
        if z["name"] not in peri:
            continue
        m = yv * z["L"] + 1.33 * zwerg_gas(z, geschichte, yv)
        rJ = peri[z["name"]] * (m / K.M_B_MILCHSTRASSE)**0.25 / math.sqrt(2)
        aus.append((rJ / z["R_e"], math.log10(z["sigma"] / zwerg_streuung(z, yv, False, False, False, None)), z["name"],
                    z["L"]))
    return sorted(aus, key=lambda a: a[1])


def satelliten_bahnen():
    """Satelliten der Milchstraße mit Eigenbewegung: heutiger Abstand, Peri- und Apozentrum [kpc] der Bahn im Feld der
    Milchstraße als Punktmasse nach der Regel (Sonne nach kern.SONNE_GALAKTISCH). Gibt {Name: (r, r_peri, r_apo)}."""
    import astropy.coordinates as ko
    import astropy.units as u
    import kern as K
    s = K.SONNE_GALAKTISCH
    rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                               galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
    aus = {}
    for z in zwerge():
        if z["wirt"] != "mw" or z["pmra"] is None or z["vlos"] is None:
            continue
        c = ko.SkyCoord(ra=z["ra"] * u.deg, dec=z["dec"] * u.deg, distance=z["abstand"] * u.kpc,
                        pm_ra_cosdec=z["pmra"] * u.mas / u.yr, pm_dec=z["pmdec"] * u.mas / u.yr,
                        radial_velocity=z["vlos"] * u.km / u.s).transform_to(rahmen)
        x = np.array([c.x.to(u.m).value, c.y.to(u.m).value, c.z.to(u.m).value])
        v = np.array([c.v_x.to(u.m / u.s).value, c.v_y.to(u.m / u.s).value, c.v_z.to(u.m / u.s).value])
        r = np.linalg.norm(x)
        v_r = x @ v / r
        v_t = math.sqrt(max(v @ v - v_r**2, 0.0))
        peri, apo = K.bahn_wendepunkte(K.M_B_MILCHSTRASSE * K.M_sun, r, v_r, v_t)
        aus[z["name"]] = (r / K.kpc, peri / K.kpc, apo / K.kpc)
    return aus


@lru_cache(maxsize=None)
def satelliten_schar(wirt="mw", ziehungen=400, haufen=False, wolken=False, begleiter=False):
    """Satelliten als Schar von Probekörpern im Feld ihres Wirts (Local Volume Database, bestätigte Zwerge mit Wirt mw,
    lmc oder m_031, ohne LMC und SMC selbst, wolken: mit ihnen; haufen: statt der Zwerge die Kugelsternhaufen der
    Milchstraße, gc_harris und gc_mw_new, mit ihrem Wirt). mw (und lmc): galaktozentrisch r [kpc], v_r und v_t [km/s]
    (Sonne nach kern.SONNE_GALAKTISCH) und die Varianz der Messfehler von v_r und v_t² aus ziehungen Ziehungen von
    Entfernungsmodul, Eigenbewegung und Sichtliniengeschwindigkeit, dazu Ort x [kpc] und Geschwindigkeit v [km/s] und ihre
    Ziehungen x_z, v_z. m_031: projizierter Abstand R [kpc] und Sichtliniengeschwindigkeit gegen
    M31 v_gsr − v_M31·n [km/s] (M31 mit Eigenbewegung) mit Fehler, begleiter: dazu die Galaxien der Tafel misc mit Wirt
    m_031 (M33). Gibt Liste von dict mit name, wirt, M_V, tafel."""
    import csv
    import astropy.coordinates as ko
    import astropy.units as u
    import kern as K
    s = K.SONNE_GALAKTISCH
    rahmen = ko.Galactocentric(galcen_distance=s["R0"] * u.kpc, z_sun=s["z"] * u.kpc,
                               galcen_v_sun=ko.CartesianDifferential(list(s["v"]) * u.km / u.s))
    zeilen = list(csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))))
    f = lambda x: float(x) if x not in ("", "nan") else math.nan
    rng = np.random.default_rng(1)

    def galaktisch(z, n):                                        # Index 0: die Bestwerte, danach die Ziehungen
        rausch = lambda: np.concatenate([[0.0], rng.standard_normal(n)[1:]])
        dm = f(z["distance_modulus"]) + f(z["distance_modulus_em"] or 0) * rausch()
        zufall = lambda k: f(z[k]) + 0.5 * (f(z[k + "_em"] or 0) + f(z[k + "_ep"] or 0)) * rausch()
        c = ko.SkyCoord(ra=np.full(n, f(z["ra"])) * u.deg, dec=np.full(n, f(z["dec"])) * u.deg,
                        distance=10**(dm / 5 - 2) * u.kpc, pm_ra_cosdec=zufall("pmra") * u.mas / u.yr,
                        pm_dec=zufall("pmdec") * u.mas / u.yr,
                        radial_velocity=zufall("vlos_systemic") * u.km / u.s).transform_to(rahmen)
        return (np.array([c.x.value, c.y.value, c.z.value]).T,
                np.array([c.v_x.to(u.km / u.s).value, c.v_y.to(u.km / u.s).value, c.v_z.to(u.km / u.s).value]).T)
    aus = []
    if wirt == "mw":
        for z in zeilen:
            if ((z["table"] not in ("gc_harris", "gc_mw_new") if haufen else
                 (z["table"] != "dwarf_mw" or z["confirmed_galaxy"] != "1" or (z["key"] in ("lmc", "smc") and not wolken)))
                    or not all(math.isfinite(f(z[k])) for k in ("pmra", "pmdec", "vlos_systemic", "distance_modulus"))):
                continue
            x, v = galaktisch(z, ziehungen + 1)
            r = np.linalg.norm(x, axis=1)
            v_r = np.einsum("ij,ij->i", x, v) / r
            v_t2 = np.einsum("ij,ij->i", v, v) - v_r**2
            aus.append(dict(name=z["key"], wirt=z["host"], M_V=f(z["M_V"]), r=float(r[0]), v_r=float(v_r[0]),
                            v_t=float(math.sqrt(max(v_t2[0], 0.0))), e_r2=float(np.var(v_r[1:])),
                            e_t2=float(np.trace(np.cov(v[1:].T)) - np.var(v_r[1:])), x=tuple(x[0]), v=tuple(v[0]),
                            x_z=x[1:], v_z=v[1:]))
    else:
        m = next(z for z in zeilen if z["key"] == "m_031")
        xm, vm = galaktisch(m, 1)
        sonne = np.array([-s["R0"], 0.0, s["z"]])
        for z in zeilen:
            if (z["host"] != "m_031" or not ((z["table"].startswith("dwarf") and z["confirmed_galaxy"] == "1")
                                              or (begleiter and z["table"] == "misc"))
                    or not math.isfinite(f(z["velocity_gsr"]))):
                continue
            o = ko.SkyCoord(ra=f(z["ra"]) * u.deg, dec=f(z["dec"]) * u.deg, distance=f(m["distance"]) * u.kpc)
            o = o.transform_to(rahmen)
            n = np.array([o.x.value, o.y.value, o.z.value]) - sonne
            n /= np.linalg.norm(n)
            R = f(m["distance"]) * math.acos(min(1.0, float(n @ (xm[0] - sonne)) / np.linalg.norm(xm[0] - sonne)))
            aus.append(dict(name=z["key"], wirt="m_031", M_V=f(z["M_V"]), R=R, v=float(f(z["velocity_gsr"]) - vm[0] @ n),
                            e2=(0.5 * (f(z["vlos_systemic_em"] or 0) + f(z["vlos_systemic_ep"] or 0)))**2,
                            tafel=z["table"]))
    return aus


def lmc_gruppe(dauer=3.0, geschichte=False, schritte=3001, alle=False):
    """Milchstraße, LMC, SMC und die übrigen Zwerge mit Wirt lmc (satelliten_schar) gemeinsam rückwärts in der Zeit,
    dauer [Gyr] (kern.bahnen_gemeinsam): Milchstraße (kern.M_B_MILCHSTRASSE), LMC und SMC (satellit, mit geschichte das
    Gas, das sie einmal trugen) ziehen einander, die LMC im Bereich der Milchstraße, die SMC in dem der LMC, die Zwerge
    als Probekörper in allen Feldern. alle: dazu alle übrigen Satelliten der Milchstraße mit Eigenbewegung (zur
    Gegenprobe). Gibt Zeiten [Gyr], Abstand LMC zur Milchstraße D [kpc], Jacobi-Radius des Paars r_J = D·(m/M)^{1/4}/√2
    mit m = m_LMC + m_SMC [kpc], m [M☉] und je Zwerg Abstand [kpc] und Tempo [km/s] gegen die LMC."""
    import kern as K
    s = satelliten_schar("mw", wolken=True)
    lmc = next(x for x in s if x["name"] == "lmc")
    gruppe = sorted((x for x in s if x["wirt"] == "lmc" or (alle and x["name"] != "lmc")), key=lambda x: x["name"] != "smc")
    wahl = lambda x, k: np.array(x[k]) * (K.kpc if k == "x" else 1e3)
    masse = lambda n, g: satellit(n, g)["M_b"] * K.M_sun
    M = K.M_B_MILCHSTRASSE * K.M_sun
    koerper = [(masse("lmc", geschichte), masse("lmc", False), wahl(lmc, "x"), wahl(lmc, "v"))] + [
        (masse("smc", geschichte), masse("smc", False), wahl(g, "x"), wahl(g, "v")) if g["name"] == "smc"
        else (0.0, 0.0, wahl(g, "x"), wahl(g, "v")) for g in gruppe]
    t, x, v, _ = K.bahnen_gemeinsam(M, koerper, dauer * K.Gyr, schritte, eltern=[0] + [1 if g["name"] == "smc" else 0
                                                                                     for g in gruppe])
    x, v = x / K.kpc, v / 1e3
    D = np.linalg.norm(x[0], axis=0)
    m = masse("lmc", geschichte) + masse("smc", geschichte)
    return dict(t=t / K.Gyr, D=D, r_J=D * (m / M)**0.25 / math.sqrt(2), m=m / K.M_sun,
                abstand={g["name"]: np.linalg.norm(x[1 + i] - x[0], axis=0) for i, g in enumerate(gruppe)},
                tempo={g["name"]: np.linalg.norm(v[1 + i] - v[0], axis=0) for i, g in enumerate(gruppe)})


@lru_cache(maxsize=None)
def lmc_gruppe_rand():
    """Wer beim letzten Apozentrum der LMC (rückwärts das erste Maximum von D, lmc_gruppe mit alle) und wer heute innerhalb
    des Rands des Paars liegt: Zeit [Gyr], D [kpc] und r_J [kpc] heute und im Apozentrum, die Wendepunkte von D der
    letzten 3 Gyr (Zeit [Gyr], D [kpc]), die Zwerge der LMC innerhalb (Zahl, von), die übrigen (Zahl, von) und die Namen
    derer, die heute innerhalb liegen."""
    wirt = {x["name"]: x["wirt"] for x in satelliten_schar("mw", wolken=True)}
    wende = lambda f: [k for k in range(1, len(f) - 1) if (f[k] - f[k - 1]) * (f[k + 1] - f[k]) < 0]
    g = lmc_gruppe(3.0, alle=True)
    k = next(k for k in wende(g["D"]) if g["D"][k] > g["D"][k - 1])
    a = [n for n, d in g["abstand"].items() if d[k] < g["r_J"][k]]
    eigen = [n for n in g["abstand"] if wirt[n] == "lmc"]
    return dict(t_apo=float(g["t"][k]), D=(g["D"][0], g["D"][k]), r_J=(g["r_J"][0], g["r_J"][k]),
                wende=[(float(g["t"][j]), float(g["D"][j])) for j in wende(g["D"])],
                eigen=(sum(n in a for n in eigen), len(eigen)), uebrige=(sum(n not in eigen for n in a), len(g["abstand"]) - len(eigen)),
                heute=sorted(n for n, d in g["abstand"].items() if d[0] < g["r_J"][0]))


@lru_cache(maxsize=None)
def satelliten_phasen(hell=-8.0, faktor=1.0, haufen=False, r_min=0.0):
    """Phasen der Bahnen der Satelliten der Milchstraße heller als M_V = hell (satelliten_schar, ohne die der LMC; haufen:
    die Kugelsternhaufen ohne die von Sagittarius) ab r_min [kpc] im Feld der Regel von faktor·kern.M_B_MILCHSTRASSE
    (kern.bahn_phase). Gibt {Name: Phase}, inf für ungebunden."""
    import kern as K
    M = faktor * K.M_B_MILCHSTRASSE * K.M_sun
    return {x["name"]: K.bahn_phase(M, x["r"] * K.kpc, x["v_r"] * 1e3, x["v_t"] * 1e3)
            for x in satelliten_schar("mw", haufen=haufen) if x["wirt"] == "mw" and x["M_V"] < hell and x["r"] > r_min}


def satelliten_phasen_masse(hell=-8.0, haufen=False, r_min=0.0):
    """Baryonen der Milchstraße [M☉], die die Phasen der Satelliten (satelliten_phasen) lesen: Mittel der Phasen ½,
    dazu die Grenzen bei ½ ∓ 1/√(12·n). Gibt (M, M_unten, M_oben, n)."""
    import kern as K
    from scipy.optimize import brentq
    n = len(satelliten_phasen(hell, 1.0, haufen, r_min))
    mittel = lambda l: float(np.mean(list(satelliten_phasen(hell, round(math.exp(l), 4), haufen, r_min).values())))
    e = 1 / math.sqrt(12 * n)
    f = [math.exp(brentq(lambda l: mittel(l) - z, math.log(0.2), math.log(3.0), xtol=1e-3)) for z in (0.5, 0.5 - e, 0.5 + e)]
    return f[0] * K.M_B_MILCHSTRASSE, f[1] * K.M_B_MILCHSTRASSE, f[2] * K.M_B_MILCHSTRASSE, n


def satelliten_perizentrum(yv=2.0):
    """Korrelation des Rests (Regel ohne Wirtsfeld) der Satelliten der Milchstraße mit lg r_peri ihrer Bahn."""
    Z = {z["name"]: z for z in zwerge()}
    b = satelliten_bahnen()
    x = [math.log10(b[n][1]) for n in b]
    r = [math.log10(Z[n]["sigma"] / zwerg_streuung(Z[n], yv, False)) for n in b]
    return float(np.corrcoef(x, r)[0, 1]), len(b)


def zwerge_gezeit(yv=2.0):
    """Korrelation des Rests lg(σ_gemessen/σ_Regel) der Satelliten (Regel ohne Wirtsfeld) mit dem Verhältnis ihrer Dichte
    in 4/3·R_e zu der des Wirts in ihrem Abstand, beide mit der Regel verstärkt (Gezeit)."""
    import kern as K
    x, r = [], []
    for z in zwerge():
        if not z["wirt"]:
            continue
        M = (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun
        rs, D = 4 / 3 * z["R_e"] * K.kpc, z["d_wirt"] * K.kpc
        mw = {"mw": K.M_B_MILCHSTRASSE, "m_031": K.M_B_M31}[z["wirt"]] * K.M_sun
        rho_s = K.verstaerkung(K.G * M / 2 / rs**2) * M / 2 / rs**3
        rho_w = K.verstaerkung(K.G * mw / D**2) * mw / D**3
        x.append(math.log10(rho_s / rho_w))
        r.append(math.log10(z["sigma"] / zwerg_streuung(z, yv, False)))
    return float(np.corrcoef(x, r)[0, 1])


def hi_frueh():
    """Neutraler Wasserstoff in Galaxien bei z ≈ 0,34 (Bera u. a. 2019, H0 = 70): (z, ρ_HI [M☉/Mpc³], Fehler der
    Statistik, Fehler der Feldschwankung)."""
    z = [float(v) for v in open(os.path.join(ORT, "hi_frueh", "bera2019.txt"), encoding="utf-8").readlines()[-1].split()]
    return z[0], z[1], z[2], z[3]


def hi_frueh_eins():
    """Neutrales Gas sternbildender Galaxien bei z ≈ 1 gegen heute bei gleicher Sternmasse (Chowdhury u. a. 2022):
    (z, Verhältnis, Fehler, Zahl der Galaxien)."""
    z, m, e, m0, e0, n = np.loadtxt(os.path.join(ORT, "hi_frueh", "chowdhury2022.txt"))
    return float(z), float(m / m0), float(m / m0 * math.hypot(e / m, e0 / m0)), int(n)


def signal_21cm_messung():
    """Das mittlere 21-cm-Signal des Himmels (Ordner hi21): EDGES (Bowman u. a. 2018) mit Mitte und Breite [MHz] und Tiefe
    [K], SARAS 3 (Singh u. a. 2022) mit Band [MHz] und Vertrauen, mit dem es das Profil von EDGES verwirft [%]."""
    m, b, t = np.loadtxt(os.path.join(ORT, "hi21", "edges2018.txt"))
    lo, hi, v = np.loadtxt(os.path.join(ORT, "hi21", "saras3.txt"))
    return dict(mitte=float(m), breite=float(b), tiefe=float(t), band=(float(lo), float(hi)), verworfen=float(v))


def reflex_gemessen():
    """Bewegung der Scheibe gegen den äußeren Halo (Petersen und Peñarrubia 2021, Ordner reflex): dict(v, v_oben, v_unten
    [km/s], l, l_fehler, b, b_fehler [°], Fehler je Mittel aus oben und unten)."""
    v, vo, vu, l, lo, lu, b, bo, bu = np.loadtxt(os.path.join(ORT, "reflex", "petersen2021.txt"))
    return dict(v=float(v), v_oben=float(vo), v_unten=float(vu), l=float(l), l_fehler=float((lo + lu) / 2), b=float(b),
                b_fehler=float((bo + bu) / 2))


def reflex_scheibe(namen=("lmc", "smc", "sagittarius_1"), weg="rahmen"):
    """Rückstoß der Scheibe: Läuft der Bereich mit dem Schwerpunkt seiner Baryonen, bewegt sich die Scheibe (Baryonen
    kern.M_B_MILCHSTRASSE) gegen ihn mit v = −Σm_i·v_i/(M + Σm_i), m_i und v_i der Satelliten (satellit). weg "rahmen": im
    Rahmen von satellit (z längs des Drehimpulses, zurück nach galaktisch gespiegelt), "galaktisch": zweiter Weg, Ort und
    Geschwindigkeit direkt im galaktozentrischen Rahmen von astropy. Gibt (|v| [km/s], l, b [°] des Ziels, Σm_i/M)."""
    import kern as K
    p, m = np.zeros(3), 0.0
    for n in namen:
        s = satellit(n)
        if weg == "rahmen":
            v = np.array(s["v"]) * np.array([1.0, -1.0, -1.0])
        else:
            import csv
            import astropy.coordinates as ko
            import astropy.units as u
            z = next(z for z in csv.DictReader(open(os.path.join(ORT, "lvdb", "comb_all.csv"))) if z["key"] == n)
            g = K.SONNE_GALAKTISCH
            c = ko.SkyCoord(ra=float(z["ra"]) * u.deg, dec=float(z["dec"]) * u.deg, distance=float(z["distance"]) * u.kpc,
                            pm_ra_cosdec=float(z["pmra"]) * u.mas / u.yr, pm_dec=float(z["pmdec"]) * u.mas / u.yr,
                            radial_velocity=float(z["vlos_systemic"]) * u.km / u.s).transform_to(
                ko.Galactocentric(galcen_distance=g["R0"] * u.kpc, z_sun=g["z"] * u.kpc,
                                  galcen_v_sun=ko.CartesianDifferential(list(g["v"]) * u.km / u.s)))
            v = np.array([c.v_x.to(u.m / u.s).value, c.v_y.to(u.m / u.s).value, c.v_z.to(u.m / u.s).value])
        p += s["M_b"] * v
        m += s["M_b"]
    w = -p / (K.M_B_MILCHSTRASSE + m) / 1e3
    r = float(np.linalg.norm(w))
    return r, math.degrees(math.atan2(w[1], w[0])) % 360, math.degrees(math.asin(w[2] / r)), m / K.M_B_MILCHSTRASSE


def gaia_flucht():
    """Fluchtgeschwindigkeit am Ort der Sonne (Deason u. a. 2019, Ordner gaia_flucht): dict(r0 [kpc], v [km/s], oben,
    unten, grenze in r₂₀₀, M200 [M☉], rand [m] = grenze·r₂₀₀ mit r₂₀₀ aus M₂₀₀ bei 200-facher kritischer Dichte heute,
    H₀ des Kerns)."""
    import kern as K
    r0, v, hi, lo, g, m = np.loadtxt(os.path.join(ORT, "gaia_flucht", "deason2019.txt"))
    rho_c = 3 * (K.H0_kms * 1e3 / K.Mpc)**2 / (8 * math.pi * K.G)
    r200 = (3 * m * 1e12 * K.M_sun / (4 * math.pi * 200 * rho_c))**(1 / 3)
    return dict(r0=float(r0), v=float(v), oben=float(hi), unten=float(lo), grenze=float(g), M200=float(m) * 1e12,
                rand=float(g) * r200)


def erosita_cgm():
    """Röntgenleuchtkraft des heißen Gases um Zentralgalaxien bis R_500c (eROSITA, Zhang u. a. 2024, Ordner erosita_cgm):
    {"sf": [(Median lg M* [M☉], L_X [10⁴⁰ erg/s]), …], "q": […]} für sternbildende und erloschene Galaxien, je drei Stücke."""
    aus = {"sf": [], "q": []}
    for z in open(os.path.join(ORT, "erosita_cgm", "zhang2024.txt")):
        if z.startswith("#") or not z.strip():
            continue
        t = z.split()
        aus[t[0]].append((float(t[1]), float(t[2])))
    return aus


def cchp_h0():
    """H₀ der Leiter über die Spitze der Riesenäste und die Kohlenstoffsterne (Freedman u. a. 2024, Ordner cchp): je Weg
    ("trgb" mit HST und JWST, "trgb_jwst", "jagb_jwst") (H₀ [km/s/Mpc], Fehler aus Statistik, System und Supernovae
    quadratisch zusammen, Zahl der Supernovae)."""
    z = np.loadtxt(os.path.join(ORT, "cchp", "freedman2024.txt"))
    return {n: (float(r[1]), float(math.sqrt(r[2]**2 + r[3]**2 + r[4]**2)), int(r[5]))
            for n, r in zip(("trgb", "trgb_jwst", "jagb_jwst"), z)}


def warp_gaia():
    """Präzession des Warps der Milchstraße aus Gaia DR2 (Poggio u. a. 2020): (Rate [km/s/kpc], rechtläufig positiv,
    statistischer Fehler, systematischer Fehler, Zahl der Sterne)."""
    r, e1, e2, n = np.loadtxt(os.path.join(ORT, "warp", "poggio2020.txt"))
    return float(r), float(e1), float(e2), int(n)


@lru_cache(maxsize=None)
def rc100():
    """RC100 (Nestor Shachar u. a. 2023, Tabelle 3): Felder z, lg M_baryon, Fehler, R_e [kpc], Fehler, f_DM(<R_e),
    Fehler, V_c(R_e) [km/s], Fehler als Spalten eines Arrays (100 Zeilen)."""
    return np.loadtxt(os.path.join(ORT, "rc100", "rc100_tab3.txt"), usecols=range(1, 10))


def rc100_groesse(punkt, ziehungen=0):
    """Größe der Scheiben in RC100 über die Zeit: lg R_e = c + b·(lg M_baryon − 10,8) + β·lg(H(z)/H₀), kleinste Quadrate
    über alle 100, H(z) im Atmen des Punkts (kern.atmen_punkt). Gibt (β, Fehler, b, Fehler, Streuung, Zahl), mit
    ziehungen > 0 dazu das Mittel und die Streuung von β über Ziehen der Galaxien."""
    import kern as K
    H = K.atmen_punkt(punkt)[4]
    d = rc100()
    z, lm, y = d[:, 0], d[:, 1], np.log10(d[:, 3])
    A = np.vstack([np.ones(len(z)), lm - 10.8, np.log10([H(1 / (1 + x)) / H(1.0) for x in z])]).T
    c, res = np.linalg.lstsq(A, y, rcond=None)[:2]
    s = math.sqrt(res[0] / (len(z) - 3))
    f = s * np.sqrt(np.diag(np.linalg.inv(A.T @ A)))
    aus = (float(c[2]), float(f[2]), float(c[1]), float(f[1]), s, len(z))
    if ziehungen:
        rng = np.random.default_rng(3)
        bs = [np.linalg.lstsq(A[i], y[i], rcond=None)[0][2] for i in (rng.integers(0, len(z), len(z)) for _ in range(ziehungen))]
        aus += (float(np.mean(bs)), float(np.std(bs)))
    return aus


@lru_cache(maxsize=None)
def rc100_schwelle(von=0.0, bis=3.0, m=0.0, fehler="f", Om=0.315):
    """Schwelle a₀ aus dem Anteil der Dunklen Materie innerhalb R_e: Im Modell ist er der Rückkehranteil,
    f_DM = V_DM²/V_c² = (a − a_N)/a = q = e^{−√(g_N/a₀)}, mit g_N = (1 − f_DM)·V_c²/R_e. Angepasst über χ² an f_DM für
    die Galaxien mit von ≤ z < bis. m: a₀ wüchse wie (H(z)/H₀)^m (Ω_m), angepasst wird dann a₀ heute. fehler: f nur
    der Fehler von f_DM, voll dazu die Fehler von M_baryon und R_e in g_N fortgepflanzt, denn g_N ist das Feld der
    Baryonen aus Masse und Radius, V_c steckt nur in f_DM. (a₀ [10⁻¹⁰ m/s²], unten, oben aus Δχ² = 1, Zahl, χ²,
    χ² bei a₀ = 1,06)."""
    from scipy.optimize import minimize_scalar
    d = rc100()
    z, eM, R, eR, f, ef, V = d[:, 0], d[:, 2], d[:, 3], d[:, 4], d[:, 5], d[:, 6], d[:, 7]
    w = (z >= von) & (z < bis)
    z, eM, R, eR, f, ef, V = z[w], eM[w], R[w], eR[w], f[w], ef[w], V[w]
    gN = (1 - f) * (V * 1e3)**2 / (R * 3.0857e19)
    E = np.sqrt(Om * (1 + z)**3 + 1 - Om)**m
    slg = np.hypot(math.log(10) * eM, 2 * eR / R)

    def chi2(a0):
        y = gN / (a0 * 1e-10 * E)
        fp = np.exp(-np.sqrt(y))
        if fehler == "f":
            sig = ef
        else:
            sig = np.hypot(ef, np.sqrt(y) / 2 * fp * slg)
        return float(np.sum(((f - fp) / sig)**2))

    r = minimize_scalar(lambda x: chi2(math.exp(x)), bounds=(math.log(1e-4), math.log(100.0)), method="bounded",
                        options={"xatol": 1e-8})
    a, c = math.exp(r.x), r.fun
    from scipy.optimize import brentq
    unten = brentq(lambda x: chi2(x) - c - 1, a / 1e3, a)
    oben = brentq(lambda x: chi2(x) - c - 1, a, a * 1e3)
    return a, unten, oben, int(w.sum()), c, chi2(1.06)


@lru_cache(maxsize=None)
def rc100_gang(fehler="f"):
    """Exponent m in a₀ ∝ (H(z)/H₀)^m aus RC100, a₀ heute je m angepasst: (m, unten, oben aus Δχ² = 1)."""
    from scipy.optimize import minimize_scalar, brentq
    c = lambda m: rc100_schwelle(m=m, fehler=fehler)[4]
    r = minimize_scalar(c, bounds=(-2.0, 3.0), method="bounded", options={"xatol": 1e-6})
    return r.x, brentq(lambda m: c(m) - r.fun - 1, -2.0, r.x), brentq(lambda m: c(m) - r.fun - 1, r.x, 3.0)


@lru_cache(maxsize=None)
def haufenstoss():
    """Zwei Haufenstöße, aus den Tabellen abgeschrieben (QUELLE in haufenstoss): dict Datei → dict Teil → Zahlen.
    clowe: 1E0657−56 in 100 kpc (M_X, Fehler, M_stern, Fehler, κ, Fehler), bradac2006: 1E0657−56 in Ellipsen
    (Ra, Rb, M_gas, Fehler, M_total, Fehler) und Lagen (dx, Fehler, dy, Fehler), bradac2008: MACS J0025.4−1222, Lagen und
    Massen in 500 kpc (M, unten, oben)."""
    aus = {}
    for k, n in (("clowe", "clowe2006_tab2.txt"), ("bradac2006", "bradac2006_tab.txt"), ("bradac2008", "bradac2008_tab.txt")):
        aus[k] = {z.split()[0]: tuple(float(v) for v in z.split()[1:])
                  for z in open(os.path.join(ORT, "haufenstoss", n), encoding="utf-8") if z.strip() and not z.startswith("#")}
    return aus


def haufenstoss_lage():
    """Abstand des Gipfels der Linse von den Galaxien und vom Gasgipfel in Einheiten seines Fehlers, √χ² über beide
    Richtungen, für die vier Teilhaufen (1E0657−56 Haupt- und Nebenhaufen gegen die Zentralgalaxie, MACS J0025.4−1222
    Südost und Nordwest gegen den Schwerpunkt der Galaxien): Liste (Name, zu den Galaxien, zum Gas)."""
    d = haufenstoss()
    b, m = d["bradac2006"], d["bradac2008"]
    abst = lambda p, q: math.sqrt(((p[0] - q[0])**2 / (p[1]**2 + q[1]**2)) + ((p[2] - q[2])**2 / (p[3]**2 + q[3]**2)))
    null = (0.0, 0.0, 0.0, 0.0)
    aus = []
    for name, l, g, x in (("1E0657 Haupt", b["lage_haupt_linse"], null, b["lage_haupt_gas"]),
                          ("1E0657 Neben", b["lage_neben_linse"], null, b["lage_neben_gas"]),
                          ("MACS J0025 SO", m["lage_so_linse"], m["lage_so_galaxien"], m["lage_so_gas"]),
                          ("MACS J0025 NW", m["lage_nw_linse"], m["lage_nw_galaxien"], m["lage_nw_gas"])):
        aus.append((name, abst(l, g), abst(l, x)))
    return aus


@lru_cache(maxsize=None)
def hiflugcs():
    """HIFLUGCS, 10 Haufen (Li u. a. 2023, Massen aus der Bewegung der Haufengalaxien; über Monjo und Banik 2025): Liste
    (Name, R [kpc], lg g_bar, lg g_tot, Fehler unten, oben [dex])."""
    aus = []
    for z in open(os.path.join(ORT, "linsen_haufen", "clusterRAR.dat")):
        if z.startswith("#") or not z.strip():
            continue
        t = z.split()
        aus.append((t[0],) + tuple(float(v) for v in t[2:7]))
    return tuple(aus)


def hiflugcs_reste(lesart="form"):
    """Rest lg(g_tot/Vorhersage) aller HIFLUGCS-Punkte: lesart "form" (1 + Q)·g_bar, das Maß des Ganzen, "regel" die örtliche
    Regel auf g_bar. Gibt (Zahl, gewichtetes Mittel, χ²) mit dem mittleren Fehler je Punkt."""
    import kern as K
    Q = K.q_je_radiant()
    r, e = [], []
    for _, _, lgb, lgt, lo, hi in hiflugcs():
        gb = 10**lgb
        p = (1 + Q) * gb if lesart == "form" else K.verstaerkung(gb) * gb
        r.append(lgt - math.log10(p))
        e.append((lo + hi) / 2)
    r, e = np.array(r), np.array(e)
    return len(r), float(np.sum(r / e**2) / np.sum(1 / e**2)), float(np.sum((r / e)**2))


@lru_cache(maxsize=None)
def nome_wahl(tiefe="baryonen"):
    """Regel mit dem Rückkehranteil e^{−s}, a_N·(1 + n_B), gegen die mit dem Kern-Nome q_Δs = Δs/r, a_N·(1 + n_B·∏(1 − q_Δsⁿ)²⁴),
    an SPARC (Υ_d = 0,5, ohne Störgrößen), KiDS-1000 (isolierte Linsen, Punktmasse M_b), HIFLUGCS und X-COP (am Nome des Ganzen,
    n_B = Q). tiefe "baryonen": q_Δs = g_bar·R/(3c²), "feld": aus der ganzen Masse des Feldes, g_obs·R/(3c²) (X-COP hydrostatisch).
    Gibt {Daten: (Zahl, größter Unterschied |lg(a_Δs/a_e)| [dex], Δχ² = χ²(Δs/r) − χ²(e^{−s}))}, Fehler wie kids_chi2 und
    hiflugcs_reste, SPARC mit der Streuung der Regel, X-COP ohne χ²."""
    import kern as K
    Q, feld, aus = K.q_je_radiant(), tiefe == "feld", {}

    def auswerten(name, punkte):
        d = dchi = 0.0
        for gb, go, R, fehler, ganzes in punkte:
            n_b = Q if ganzes else K.rueckkehr_summe(K.s_(gb))
            q = (go if feld else gb) * R / (3 * K.c**2)
            a_e, a_d = (1 + n_b) * gb, (1 + n_b * K.spitzenfaktor(q)) * gb
            d = max(d, abs(math.log10(a_d / a_e)))
            if fehler:
                dchi += (math.log10(go / a_d) / fehler(go / a_d))**2 - (math.log10(go / a_e) / fehler(go / a_e))**2
        aus[name] = (len(punkte), d, float(dchi))

    p = []
    for _, _, k, _, _ in sparc_auswahl():
        R = k["R"] * K.kpc
        gb = (k["Vgas"] * np.abs(k["Vgas"]) + YD * k["Vdisk"]**2 + YB * k["Vbul"]**2) * 1e6 / R
        go = k["V"]**2 * 1e6 / R
        sg = np.sqrt((2 * k["eV"] / (k["V"] * math.log(10)))**2 + STREUUNG**2)
        p += [(gb[i], go[i], R[i], lambda x, e=sg[i]: e, False) for i in range(len(R))]
    auswerten("sparc", p)
    p = []
    for b in kids_isoliert():
        R = b["R"] * K.kpc
        gN = K.G * b["M_b"] * K.M_sun / R**2
        for i in range(len(R)):
            lo, hi = math.log10(b["g"][i] / b["g_lo"][i]), math.log10(b["g_hi"][i] / b["g"][i])
            p.append((gN[i], b["g"][i], R[i], lambda x, lo=lo, hi=hi: lo if x > 1 else hi, False))
    auswerten("kids", p)
    auswerten("hiflugcs", [(10**lgb, 10**lgt, Rk * K.kpc, lambda x, e=(lo + hi) / 2: e, True)
                           for _, Rk, lgb, lgt, lo, hi in hiflugcs()])
    p = []
    for h in xcop().values():
        R = h["R"] * K.kpc
        p += [(K.G * h["Mgas"][i] * K.M_sun / R[i]**2, K.G * h["M"][i] * K.M_sun / R[i]**2, R[i], None, True)
              for i in range(len(R))]
    auswerten("xcop", p)
    return aus


@lru_cache(maxsize=None)
def davis():
    """Davis u. a. 2011: 68 Kinasehemmer gegen 442 Kinasen, K_d [nM], 10000 heißt keine Bindung gemessen. Das Feld liegt
    als NumPy-Pickle vor und wird nur mit einem gesperrten Leser geladen, der allein die Bausteine von Feldern zulässt.
    Gibt das Feld der pK_d = 9 − lg K_d (NaN ohne Bindung)."""
    import pickle

    class Gesperrt(pickle.Unpickler):
        def find_class(self, mod, name):
            if (mod.replace("numpy._core", "numpy.core"), name) in {("numpy.core.multiarray", "_reconstruct"),
                                                                     ("numpy", "ndarray"), ("numpy", "dtype")}:
                return getattr(np.core.multiarray if name == "_reconstruct" else np, name)
            raise pickle.UnpicklingError(f"nicht zugelassen: {mod}.{name}")
    with open(os.path.join(ORT, "davis", "Y"), "rb") as f:
        Y = np.asarray(Gesperrt(f, encoding="latin1").load(), dtype=float)
    return np.where(Y < 10000, 9 - np.log10(Y), np.nan)


@lru_cache(maxsize=None)
def gtopdb_wirkweise(nach="InChIKey"):
    """Guide to Pharmacology (GtoPdb 2026.3.1, interactions.csv und ligands.csv): je Paar aus Ziel (UniProt) und Wirkstoff
    (nach InChIKey oder ChEMBL-Nummer, nach="ChEMBL ID") die Wirkweise, "voll" (Full agonist), "teil" (Partial agonist),
    "anta" (Antagonist) oder "invers" (Inverse agonist). Paare mit mehr als einer Wirkweise bleiben fort. Gibt
    {(UniProt, Wirkstoff): Wirkweise}."""
    import csv
    art = {"Full agonist": "voll", "Partial agonist": "teil", "Antagonist": "anta", "Inverse agonist": "invers"}
    with open(os.path.join(ORT, "gtopdb", "ligands.csv"), encoding="utf-8") as f:
        next(f)
        schluessel = {z["Ligand ID"]: z[nach] for z in csv.DictReader(f) if z[nach]}
    paare = {}
    with open(os.path.join(ORT, "gtopdb", "interactions.csv"), encoding="utf-8") as f:
        next(f)
        for z in csv.DictReader(f):
            a, k = art.get(z["Action"]), schluessel.get(z["Ligand ID"])
            if a and k and z["Target UniProt ID"]:
                for u in z["Target UniProt ID"].split("|"):
                    paare.setdefault((u, k), set()).add(a)
    return {p: next(iter(a)) for p, a in paare.items() if len(a) == 1}


def chembl_kinetik(neu=False):
    """ChEMBL: alle Messungen von k_off in s⁻¹ (standard_type koff oder k_off), je Zeile Molekül, Ziel als UniProt-Nummer
    (über die Querverweise von UniProt auf ChEMBL), k_off und Dokument. Abgelegt in archiv/weltmodell/daten/chembl/kinetik.csv, nur
    geholt, wenn die Datei fehlt oder neu. Gibt eine Liste von (UniProt, ChEMBL-Nummer, k_off, None, Dokument)."""
    import csv
    import json
    import time
    import urllib.parse
    import urllib.request

    def hole(url):
        for versuch in range(8):
            try:
                return json.load(urllib.request.urlopen(url, timeout=300))
            except Exception:
                time.sleep(2 ** versuch)
        raise RuntimeError(f"nicht geholt: {url}")
    pfad = os.path.join(ORT, "chembl", "kinetik.csv")
    if neu or not os.path.exists(pfad):
        basis = "https://www.ebi.ac.uk/chembl/api/data/"
        zeilen, start = [], 0
        while True:
            d = hole(basis + "activity.json?standard_type__in=koff,k_off&standard_value__isnull=false&standard_units=s-1"
                     f"&only=molecule_chembl_id,target_chembl_id,standard_value,document_chembl_id&limit=500&offset={start}")
            zeilen += d["activities"]
            if not d["page_meta"]["next"]:
                break
            start += 500
        uni = {}                                                  # Ziel → UniProt, über die Querverweise von UniProt
        ziele = sorted({a["target_chembl_id"] for a in zeilen})
        for i in range(0, len(ziele), 40):
            frage = urllib.parse.quote(" OR ".join(f"xref:chembl-{t}" for t in ziele[i:i + 40]))
            for versuch in range(8):
                try:
                    text = urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/search?query=" + frage
                                                  + "&fields=accession,xref_chembl&format=tsv&size=500", timeout=300)
                    text = text.read().decode()
                    break
                except Exception:
                    time.sleep(2 ** versuch)
            else:
                raise RuntimeError("UniProt nicht erreicht")
            for zeile in text.splitlines()[1:]:
                acc, verweise = (zeile.split("\t") + [""])[:2]
                for t in verweise.strip(";").split(";"):
                    if t.strip() in ziele:
                        uni.setdefault(t.strip(), []).append(acc)
        uni = {t: "|".join(sorted(v)) for t, v in uni.items()}
        with open(pfad, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["molekuel", "ziel", "uniprot", "koff", "dokument"])
            for a in zeilen:
                w.writerow([a["molecule_chembl_id"], a["target_chembl_id"], uni.get(a["target_chembl_id"], ""), a["standard_value"],
                            a["document_chembl_id"]])
    with open(pfad, encoding="utf-8") as f:
        return [(u, z["molekuel"], float(z["koff"]), None, z["dokument"])
                for z in csv.DictReader(f) if float(z["koff"]) > 0 for u in z["uniprot"].split("|") if u]


@lru_cache(maxsize=None)
def bindingdb_kinetik():
    """BindingDB (aus Artikeln, Stand 10.2026): alle Messungen mit k_off. Die Tafel wird aus dem Zip gelesen, ohne
    auszupacken. Gibt eine Liste von (UniProt, InChIKey, k_off [1/s], k_on [1/(M·s)] oder None, DOI oder PMID)."""
    import csv
    import io
    import zipfile
    pfad = os.path.join(ORT, "bindingdb", "BindingDB_BindingDB_Articles_202610_tsv.zip")
    zahl = lambda x: float(x) if x and x.strip()[0] not in "<>~" else None
    aus = []
    with zipfile.ZipFile(pfad) as z:
        name = next(n for n in z.namelist() if n.endswith(".tsv") and "/" not in n and ".." not in n)
        with z.open(name) as f:
            r = csv.reader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t")
            kopf = next(r)
            spalte = lambda anfang: next(i for i, k in enumerate(kopf) if k.startswith(anfang))
            i_ab, i_an, i_key = spalte("koff"), spalte("kon"), spalte("Ligand InChI Key")
            i_uni, i_doi, i_pmid = spalte("UniProt (SwissProt) Primary ID of Target Chain"), spalte("Article DOI"), spalte("PMID")
            for z_ in r:
                if len(z_) <= max(i_ab, i_uni) or not z_[i_ab].strip():
                    continue
                try:
                    ab = zahl(z_[i_ab])
                    an = zahl(z_[i_an]) if z_[i_an].strip() else None
                except ValueError:
                    continue
                if ab and ab > 0 and z_[i_uni] and z_[i_key]:
                    aus.append((z_[i_uni].strip(), z_[i_key].strip(), ab, an, z_[i_doi].strip() or z_[i_pmid].strip()))
    return aus


def wirkweise_treffer():
    """Was die Kinetik von ChEMBL (chembl_kinetik) für die Leiter der Takte hergibt: Zahl der Messungen, der Ziele und der
    Paare aus Ziel und Wirkstoff mit einer Wirkweise in GtoPdb, dazu die Wirkweisen dieser Paare. Gibt (Messungen, Ziele,
    Paare, {Wirkweise: Zahl})."""
    import csv
    from collections import Counter
    chembl_kinetik()
    with open(os.path.join(ORT, "chembl", "kinetik.csv"), encoding="utf-8") as f:
        zeilen = list(csv.DictReader(f))
    w = gtopdb_wirkweise("ChEMBL ID")
    paare = {(u, z["molekuel"]) for z in zeilen for u in z["uniprot"].split("|") if (u, z["molekuel"]) in w}
    return len(zeilen), len({z["ziel"] for z in zeilen}), len(paare), dict(Counter(w[p] for p in paare))


def wirkweise_takt(je_arbeit=False, quelle="chembl"):
    """Prüfung der Leiter der Takte: k_off (quelle "chembl" oder "bindingdb") gegen die Wirkweise (GtoPdb) an demselben
    Ziel. Je Ziel (und mit je_arbeit innerhalb derselben Arbeit) alle Paare zweier Wirkweisen verglichen, das Mittel von
    lg k_off je Wirkstoff. Gibt je Paar der Wirkweisen (schneller, langsamer laut Modell) die Zahl der Ziele, in denen das
    Modell öfter recht hat, in denen es öfter unrecht hat, und die Zahl der verglichenen Paare."""
    w = gtopdb_wirkweise("ChEMBL ID" if quelle == "chembl" else "InChIKey")
    lg = {}
    for u, k, ab, _, quelle in (chembl_kinetik() if quelle == "chembl" else bindingdb_kinetik()):
        if (u, k) in w:
            lg.setdefault((u, quelle if je_arbeit else "", k), []).append(math.log10(ab))
    gruppe = {}
    for (u, q, k), v in lg.items():
        gruppe.setdefault((u, q), {}).setdefault(w[(u, k)], []).append(float(np.mean(v)))
    aus = {}
    for schnell, langsam in (("voll", "teil"), ("teil", "anta"), ("voll", "anta"), ("anta", "invers"), ("voll", "invers")):
        recht = unrecht = paare = 0
        for g in gruppe.values():
            if schnell in g and langsam in g:
                d = [a - b for a in g[schnell] for b in g[langsam]]
                paare += len(d)
                gut, schlecht = sum(x > 0 for x in d), sum(x < 0 for x in d)
                recht += gut > schlecht
                unrecht += schlecht > gut
        aus[(schnell, langsam)] = (recht, unrecht, paare)
    return aus


def davis_atome(weg="smiles"):
    """Zahl der schweren Atome der 68 Hemmer von Davis u. a. 2011, in der Reihenfolge der Zeilen von davis(). weg "smiles":
    aus den SMILES gezählt (Atome in Klammern ohne Wasserstoff, Cl, Br, B, C, N, O, P, S, F, I und aromatisch b, c, n, o,
    p, s), "pubchem": HeavyAtomCount je PubChem-Nummer von PubChem, abgelegt in davis/atome_pubchem.csv."""
    import csv
    import json
    import re
    import urllib.request
    with open(os.path.join(ORT, "davis", "ligands_can.txt"), encoding="utf-8") as f:
        lig = json.load(f)                                       # Reihenfolge der Datei ist die der Zeilen
    if weg == "smiles":
        z = lambda s: sum(1 for t in re.findall(r"\[[^\]]+\]|Br|Cl|[BCNOPSFI]|[bcnops]", s)
                          if not re.fullmatch(r"\[\d*H[^a-z]*\]", t))
        return np.array([z(s) for s in lig.values()])
    pfad = os.path.join(ORT, "davis", "atome_pubchem.csv")
    if not os.path.exists(pfad):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/HeavyAtomCount/CSV"
        text = urllib.request.urlopen(url.format(",".join(lig)), timeout=300).read().decode()
        with open(pfad, "w", encoding="utf-8") as f:
            f.write(text)
    with open(pfad, encoding="utf-8") as f:
        da = {z["CID"]: int(z["HeavyAtomCount"]) for z in csv.DictReader(f)}
    return np.array([da[c] for c in lig])


def davis_groesse(schwelle=5.0, ziehungen=0):
    """Größe gegen Fenster an Davis u. a. 2011: je Hemmer die Zahl N der schweren Atome (davis_atome), die Zahl seiner
    Partner (pK_d ≥ schwelle) und die stärkste Bindung. Gibt (Rangkorrelation N mit Partnern, p, N mit pK_max, p,
    Exponent von Partnern über N in log), mit ziehungen > 0 dazu die halbe Breite des mittleren 68 % des Exponenten über
    so viele Ziehungen der Hemmer mit Zurücklegen."""
    from scipy.stats import spearmanr
    P, N = davis(), davis_atome()
    n, pmax = np.sum(P >= schwelle, axis=1), np.nanmax(P, axis=1)
    r1, r2 = spearmanr(N, n), spearmanr(N, pmax)
    gut = n > 0
    x, y = np.log(N[gut]), np.log(n[gut])
    aus = (float(r1.correlation), float(r1.pvalue), float(r2.correlation), float(r2.pvalue), float(np.polyfit(x, y, 1)[0]))
    if ziehungen:
        z = np.random.default_rng(1).integers(0, len(x), (ziehungen, len(x)))
        e = np.percentile([np.polyfit(x[i], y[i], 1)[0] for i in z], [16, 84])
        aus += (float(e[1] - e[0]) / 2,)
    return aus


def davis_fenster():
    """Je Hemmer die stärkste Bindung pK_d,max und die Zahl der gebundenen Kinasen (K_d < 10 µM), dazu die
    Rangkorrelation beider über alle Hemmer. Gibt (pK_max, Zahl, Rangkorrelation, p-Wert)."""
    from scipy.stats import spearmanr
    P = davis()
    pmax, n = np.nanmax(P, axis=1), np.sum(np.isfinite(P), axis=1)
    r = spearmanr(pmax, n)
    return pmax, n, float(r.correlation), float(r.pvalue)


def davis_form(boden=5.0, anteil=0.9, schwelle=None):
    """Form der Bindungen je Hemmer gegen das Einrasten bei gleichmäßig verteilter Verstimmung: y = (pK_d − Boden)/(pK_max −
    Boden) aller Nebenbindungen, gesammelt über die Hemmer mit mindestens 5 Bindungen, gegen y = ΔV/2K mit x = Δω/K
    gleichverteilt (Quantil zum Anteil a: y bei x = 1 − a). schwelle: statt der Form die Rangkorrelation der stärksten
    Bindung mit der Zahl der Bindungen über dieser pK_d. Gibt (Quantil der Daten, Quantil des Einrastens) oder (ρ, p)."""
    import kern as K
    from scipy.stats import spearmanr
    P = davis()
    if schwelle is not None:
        r = spearmanr(np.nanmax(P, axis=1), np.sum(P >= schwelle, axis=1))
        return float(r.correlation), float(r.pvalue)
    y = []
    for z in P:
        v = z[np.isfinite(z)]
        if len(v) >= 5:
            y += list(((v - boden) / (v.max() - boden))[v < v.max()])
    return float(np.quantile(y, anteil)), K.einrasten_schwelle(1 - anteil, 1.0) / 2


def davis_sequenz(k=3, grenze=0.9):
    """Verstimmung aus der Sequenz: Ähnlichkeit jeder Kinase zur stärksten Partnerin eines Hemmers (Kosinus der Häufigkeiten
    der Gruppen aus k Aminosäuren), nahezu gleiche Sequenzen (Ähnlichkeit ≥ grenze, Mutanten) ausgelassen. Gibt den Anteil
    gebundener Kinasen je Viertel der Ähnlichkeit, das ähnlichste zuerst, gemittelt über die Hemmer."""
    import json
    from collections import Counter
    seqs = list(json.load(open(os.path.join(ORT, "davis", "proteins.txt"), encoding="utf-8")).values())
    zahl = [Counter(q[i:i + k] for i in range(len(q) - k + 1)) for q in seqs]
    nr = {g: i for i, g in enumerate(sorted(set().union(*zahl)))}
    M = np.zeros((len(seqs), len(nr)))
    for i, c in enumerate(zahl):
        for g, v in c.items():
            M[i, nr[g]] = v
    M /= np.linalg.norm(M, axis=1, keepdims=True)
    S = M @ M.T
    teile = []
    for z in davis():
        b = int(np.nanargmax(z))
        o = np.argsort(-S[b])
        o = o[(o != b) & (S[b][o] < grenze)]
        teile.append([float(np.isfinite(z[g]).mean()) for g in np.array_split(o, 4)])
    return [float(x) for x in np.mean(teile, axis=0)]


def chembl_bindung(seiten=700, neu=False):
    """ChEMBL: alle Messungen von K_i und K_d mit pChEMBL (seiten·1000 Zeilen), je Zeile Molekül, Ziel, Art, pChEMBL und
    die Ligandeneffizienz LE = 1,37·pChEMBL/N (N schwere Atome). Geholt über die Schnittstelle von ChEMBL und abgelegt in
    archiv/weltmodell/daten/chembl/bindung.csv, nur geholt, wenn die Datei fehlt oder neu. Gibt die Zeilen."""
    import csv
    import json
    import time
    import urllib.request
    pfad = os.path.join(ORT, "chembl", "bindung.csv")
    if neu or not os.path.exists(pfad):
        os.makedirs(os.path.dirname(pfad), exist_ok=True)
        url = ("https://www.ebi.ac.uk/chembl/api/data/activity.json?pchembl_value__isnull=false&standard_type__in=Ki,Kd"
               "&only=molecule_chembl_id,target_chembl_id,standard_type,pchembl_value,ligand_efficiency&limit=1000&offset={}")
        from concurrent.futures import ThreadPoolExecutor

        def seite(i):
            for versuch in range(6):
                try:
                    return json.load(urllib.request.urlopen(url.format(1000 * i), timeout=300))["activities"]
                except Exception:
                    time.sleep(2 ** versuch)
            raise RuntimeError(f"Seite {i} nicht geholt")
        with open(pfad + ".teil", "w", newline="", encoding="utf-8") as f, ThreadPoolExecutor(8) as pool:
            w = csv.writer(f)
            w.writerow(["molekuel", "ziel", "art", "pchembl", "le", "bei"])
            for zeilen in pool.map(seite, range(seiten)):           # in der Reihenfolge der Seiten
                for a in zeilen:
                    le = a.get("ligand_efficiency") or {}
                    w.writerow([a["molecule_chembl_id"], a["target_chembl_id"], a["standard_type"], a["pchembl_value"],
                                le.get("le"), le.get("bei")])
                f.flush()
        os.replace(pfad + ".teil", pfad)
        with open(os.path.join(ORT, "chembl", "QUELLE.txt"), "w", encoding="utf-8") as f:
            f.write("ChEMBL über www.ebi.ac.uk/chembl/api/data/activity: alle K_i und K_d mit pChEMBL, je Zeile Molekül, Ziel, "
                    "Art, pChEMBL, Ligandeneffizienz LE = 1,37·pChEMBL/N und BEI. Geholt mit daten.chembl_bindung am "
                    + time.strftime("%d.%m.%Y") + ".\natome.csv: Zahl der schweren Atome (molecule_properties.heavy_atoms) "
                    "je Molekül, geholt mit daten.chembl_atome.\nstrukturen.csv: InChI und SMILES der Moleküle mit festgelegter "
                    "räumlicher Anordnung, geholt mit daten.chembl_strukturen.\n")
    with open(pfad, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def chembl_atome(molekuele):
    """ChEMBL: Zahl der schweren Atome je Molekül aus der Molekülschnittstelle (molecule_properties.heavy_atoms), abgelegt
    in archiv/weltmodell/daten/chembl/atome.csv (leer, wo ChEMBL keine Zahl führt), nur Fehlendes wird geholt. Gibt
    {Molekül: N}."""
    import csv
    import json
    import time
    import urllib.request
    pfad = os.path.join(ORT, "chembl", "atome.csv")
    da = {}
    if os.path.exists(pfad):
        with open(pfad, encoding="utf-8") as f:
            da = {z["molekuel"]: int(z["atome"]) if z["atome"] else None for z in csv.DictReader(f)}
    fehlt = sorted(set(molekuele) - set(da))
    da.update({m: None for m in fehlt})
    url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in={}"
           "&only=molecule_chembl_id,molecule_properties&limit=100")
    for i in range(0, len(fehlt), 100):
        for versuch in range(6):
            try:
                d = json.load(urllib.request.urlopen(url.format(",".join(fehlt[i:i + 100])), timeout=300))
                break
            except Exception:
                time.sleep(2 ** versuch)
        else:
            raise RuntimeError(f"Moleküle ab {i} nicht geholt")
        for m in d["molecules"]:
            p = m.get("molecule_properties") or {}
            if p.get("heavy_atoms") is not None:
                da[m["molecule_chembl_id"]] = int(p["heavy_atoms"])
    if fehlt:                                                    # nur schreiben, wenn etwas geholt wurde
        with open(pfad, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["molekuel", "atome"])
            w.writerows(sorted((m, "" if n is None else n) for m, n in da.items()))
    return {m: n for m, n in da.items() if n is not None}


def chembl_strukturen(neu=False):
    """ChEMBL: Strukturen der Moleküle aus chembl_bindung über die Molekülschnittstelle (standard_inchi, canonical_smiles),
    abgelegt in archiv/weltmodell/daten/chembl/strukturen.csv nur für Moleküle mit festgelegter räumlicher Anordnung
    (InChI mit Schicht /m). Gibt {Molekül: (InChI, SMILES)}."""
    import csv
    import json
    import time
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor
    pfad = os.path.join(ORT, "chembl", "strukturen.csv")
    if neu or not os.path.exists(pfad):
        namen = sorted({z["molekuel"] for z in chembl_bindung()})
        url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in={}"
               "&only=molecule_chembl_id,molecule_structures&limit=100")

        def stueck(i):
            for versuch in range(6):
                try:
                    d = json.load(urllib.request.urlopen(url.format(",".join(namen[i:i + 100])), timeout=300))
                    break
                except Exception:
                    time.sleep(2 ** versuch)
            else:
                raise RuntimeError(f"Moleküle ab {i} nicht geholt")
            aus = []
            for m in d["molecules"]:
                s = m.get("molecule_structures") or {}
                inchi = s.get("standard_inchi") or ""
                if "/m" in inchi:
                    aus.append((m["molecule_chembl_id"], inchi, s.get("canonical_smiles") or ""))
            return aus
        with open(pfad + ".teil", "w", newline="", encoding="utf-8") as f, ThreadPoolExecutor(8) as pool:
            w = csv.writer(f)
            w.writerow(["molekuel", "inchi", "smiles"])
            for zeilen in pool.map(stueck, range(0, len(namen), 100)):
                w.writerows(zeilen)
        os.replace(pfad + ".teil", pfad)
    with open(pfad, encoding="utf-8") as f:
        return {z["molekuel"]: (z["inchi"], z["smiles"]) for z in csv.DictReader(f)}


@lru_cache(maxsize=None)
def chembl_molekuele(art=("Ki", "Kd")):
    """Je Molekül in chembl_bindung die Zahl N der schweren Atome aus LE = 1,37·pChEMBL/N (LE auf 0,01 gerundet: die
    Bereiche aller Zeilen eines Moleküls geschnitten, die ganze Zahl nächst dem Mittel) und die stärkste Bindung über alle
    Ziele. Gibt (Moleküle, N, pK_max)."""
    unten, oben, mitte, best = {}, {}, {}, {}
    for z in chembl_bindung():
        if z["art"] not in art or not z["le"] or float(z["le"]) <= 0.005:
            continue
        m, p, le = z["molekuel"], float(z["pchembl"]), float(z["le"])
        unten[m] = max(unten.get(m, 0.0), 1.37 * p / (le + 0.005))
        oben[m] = min(oben.get(m, 1e9), 1.37 * p / (le - 0.005))
        mitte.setdefault(m, []).append(1.37 * p / le)
        best[m] = max(best.get(m, -1e9), p)
    namen, N = list(best), []
    for m in namen:
        ganz = np.arange(math.ceil(unten[m]), math.floor(oben[m]) + 1)
        mi = float(np.mean(mitte[m]))
        N.append(int(ganz[np.argmin(np.abs(ganz - mi))]) if len(ganz) else int(round(mi)))
    return namen, np.array(N), np.array([best[m] for m in namen])


@lru_cache(maxsize=None)
def chembl_spiegelpaare(art=("Ki", "Kd")):
    """Spiegelbilder in ChEMBL: Paare von Molekülen, deren InChI bis auf die Schicht /m0 gegen /m1 gleich ist (dieselbe
    Anordnung, gespiegelt). Zweiter Weg: in der SMILES des einen @ und @@ vertauscht gibt die SMILES des anderen (Anteil
    der Paare). Je Paar und Ziel, an dem beide gemessen sind, der Unterschied |ΔpK| der Mittel aller Messungen. Gibt dict
    mit N, dpK (je Paar und Ziel), paare, smiles_gleich (Anteil), zentren (Zahl der Zentren je Paar und Ziel)."""
    import re
    namen, Nm, _ = chembl_molekuele(tuple(art))
    groesse = dict(zip(namen, Nm))
    st = chembl_strukturen()
    gruppe = {}
    for m, (inchi, smi) in st.items():
        t = re.search(r"/m([01])", inchi)
        if t and m in groesse:
            gruppe.setdefault(re.sub(r"/m[01]", "", inchi), {})[t.group(1)] = m
    paare = [(g["0"], g["1"]) for g in gruppe.values() if len(g) == 2]
    tausch = lambda s: s.replace("@@", "\0").replace("@", "@@").replace("\0", "@")
    gleich = float(np.mean([tausch(st[a][1]) == st[b][1] for a, b in paare]))
    werte = {}
    for z in chembl_bindung():
        if z["art"] in art:
            werte.setdefault((z["molekuel"], z["ziel"]), []).append(float(z["pchembl"]))
    ziele = {}
    for (m, t) in werte:
        ziele.setdefault(m, set()).add(t)
    N, dpK, zentren = [], [], []
    for a, b in paare:
        for t in ziele.get(a, set()) & ziele.get(b, set()):
            N.append(groesse[a])
            dpK.append(abs(np.mean(werte[a, t]) - np.mean(werte[b, t])))
            zentren.append(len(re.search(r"/t([^/]*)", st[a][0]).group(1).split(",")))
    return dict(N=np.array(N), dpK=np.array(dpK), paare=len(paare), smiles_gleich=gleich, zentren=np.array(zentren))


def chembl_spiegel(grenzen=(15, 20, 25, 30, 40)):
    """Spiegelbilder nach Größe: je Bereich der Zahl schwerer Atome (bis grenzen[0], dazwischen, über grenzen[-1]) der
    Median von |ΔpK|, die Zahl der Paare mit Ziel und der Anteil ohne Unterschied (|ΔpK| < 0,1), dazu der Median über
    alle und die Rangkorrelation von |ΔpK| mit N. Gibt dict."""
    from scipy.stats import spearmanr
    s = chembl_spiegelpaare()
    k = np.digitize(s["N"], grenzen, right=True)
    r = spearmanr(s["N"], s["dpK"])
    return dict(bereiche=[(float(np.median(s["dpK"][k == i])), int((k == i).sum()), float(np.mean(s["dpK"][k == i] < 0.1)))
                          for i in range(len(grenzen) + 1)], median=float(np.median(s["dpK"])),
                rho=float(r.correlation), p=float(r.pvalue),
                paare=s["paare"], gemessen=len(s["N"]), smiles_gleich=s["smiles_gleich"])


def chembl_huelle(quantil=0.99, art=("Ki", "Kd"), mindest=200):
    """Obere Hülle der Bindung über der Größe (ChEMBL): je N (chembl_molekuele) das quantil der stärksten Bindung der
    Moleküle, nur N mit mindestens mindest Molekülen. Gibt (N, pK_q, Zahl der Moleküle)."""
    _, N, P = chembl_molekuele(tuple(art))
    stufen = [n for n in np.unique(N) if (N == n).sum() >= mindest]
    return (np.array(stufen), np.array([np.quantile(P[N == n], quantil) for n in stufen]),
            np.array([(N == n).sum() for n in stufen]))


def chembl_atome_vergleich(quantil=0.99, art=("Ki", "Kd"), mindest=200):
    """Zweiter Weg zu N: die Moleküle auf und über der Hülle (chembl_huelle) gegen die Zahl der schweren Atome aus der
    Molekülschnittstelle (chembl_atome). Gibt (Zahl der Moleküle, Anteil gleich, Anteil auf ±1 gleich)."""
    namen, N, P = chembl_molekuele(tuple(art))
    st, q, _ = chembl_huelle(quantil, art, mindest)
    grenze = dict(zip(st, q))
    wahl = [i for i in range(len(namen)) if N[i] in grenze and P[i] >= grenze[N[i]]]
    da = chembl_atome([namen[i] for i in wahl])
    paar = np.array([(N[i], da[namen[i]]) for i in wahl if namen[i] in da])
    d = np.abs(paar[:, 0] - paar[:, 1])
    return len(paar), float(np.mean(d == 0)), float(np.mean(d <= 1))


def chembl_huelle_form(quantil=0.99, b=None, art=("Ki", "Kd"), mindest=200):
    """Die obere Hülle (chembl_huelle) in der Form des gehaltenen Trichters, pK = min(a·N^b − c, P) (kern.bindung_huelle),
    kleinste Quadrate über die Stufen. b = None: Exponent frei. Gibt dict mit a, b, c, P, Treffpunkt N_T, Größe der
    stärksten Bindung je Atom N_e, Streuung."""
    import numpy as np
    from scipy.optimize import least_squares
    import kern as K
    N, q, _ = chembl_huelle(quantil, art, mindest)
    fest = b is not None

    def rest(x):
        a, c, P = x[:3]
        return K.bindung_huelle(N, a, c, P, b if fest else x[3]) - q
    x0 = [3.0, 4.0, float(q.max())] + ([] if fest else [2 / 3])
    r = least_squares(rest, x0, bounds=([0, -50, 0] + ([] if fest else [0.05]), [100, 100, 20] + ([] if fest else [3])))
    a, c, P = r.x[:3]
    be = b if fest else r.x[3]
    return dict(a=a, b=be, c=c, P=P, N_T=((P + c) / a)**(1 / be), N_e=K.bindung_je_atom_groesse(a, c, be),
                streuung=float(np.sqrt(np.mean(r.fun**2))), stufen=len(N))


def kids_scherung(winkel=(0.0, 1e9), art="+-"):
    """KiDS-1000 Scherung (Asgari u. a. 2021): ξ₊ für 0,5 bis 300′, ξ₋ für 4 bis 300′ (Schnitte wie bei KiDS) je Paar der 5
    Bins, mit der Kovarianz und den n(z); winkel: nur θ in diesem Bereich [′], art: "+", "-" oder beide. Gibt (Zeilen
    [(+/−, i, j, θ, Wert)], Kovarianz, (z, [n_1 … n_5]))."""
    from astropy.io import fits
    f = fits.open(os.path.join(ORT, "kids_scherung", "xipm_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_"
                               "svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"))
    zeilen, idx = [], []
    for name, vz, ab in (("xiP", "+", 0.5), ("xiM", "-", 4.0)):
        for k, r in enumerate(f[name].data):
            if r["ANG"] >= ab and winkel[0] <= r["ANG"] <= winkel[1] and vz in art:
                zeilen.append((vz, int(r["BIN1"]) - 1, int(r["BIN2"]) - 1, float(r["ANG"]), float(r["VALUE"])))
                idx.append(k + (0 if vz == "+" else len(f["xiP"].data)))
    C = np.array(f["COVMAT"].data)[np.ix_(idx, idx)]
    nz = f["NZ_SOURCE"].data
    return zeilen, C, (np.array(nz["Z_MID"]), [np.array(nz[f"BIN{i}"]) for i in range(1, 6)])


@lru_cache(maxsize=None)
def kids_scherung_xi(punkt, gas=False, winkel=(0.0, 1e9), art="+-"):
    """ξ± der Faltung im Punkt je Paar der Bins an den Winkeln von kids_scherung, je als (GG, GI, II) für die ganze Faltung,
    den linearen Teil und den nichtlinearen Rest; gas: der Rest mal der Faltung eines Haufens, dessen Gas nie abkühlte
    (xcop_gas_faltung). Gibt ({(i, j): {"ganz", "lin", "rest"}}, Winkel, S8)."""
    import kern as K
    zeilen, _, nz = kids_scherung(winkel, art)
    kg = np.geomspace(0.01, 100.0, 120)
    ell, cl, S8 = K.scherung_cl(*punkt[:5], nz, *punkt[5:], teile=True,
                                daempfung=(kg, xcop_gas_faltung(kg)) if gas else None)
    th = sorted({z[3] for z in zeilen})
    xi = lambda s: [K.scherung_xi(ell, c, th) for c in s]
    return {ij: dict(ganz=xi([a + b for a, b in zip(*t)]), lin=xi(t[0]), rest=xi(t[1])) for ij, t in cl.items()}, th, S8


def kids_scherung_chi2(punkt, verhaeltnis=False, teile=None, gas=False, **auswahl):
    """χ² der KiDS-1000-Scherung (kids_scherung) für den Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀, w_a): Faltung der Mode wie
    druckfreier Staub, Ausrichtung der Bilder mit freier Amplitude A; auswahl: winkel, art wie bei kids_scherung. Gibt (χ², A,
    S8, Zahl der Werte); verhaeltnis: dazu Messung gegen Faltung bei diesem A, linear geschätzt mit der vollen Kovarianz,
    R = gᵀC⁻¹(d − m_A)/gᵀC⁻¹g mit dem Fehler 1/√(gᵀC⁻¹g). teile = (R_lin, R_nl): der lineare Teil der Faltung und der
    nichtlineare Rest, den vor allem die gebundenen Bereiche tragen, je mit dieser Amplitude, None frei (bei festem A linear
    gelöst), gibt dann dazu R_lin, R_nl, ihre Fehler und ihre Korrelation. gas: der Rest mal der Faltung eines Haufens, dessen
    Gas nie abkühlte (kids_scherung_xi)."""
    from scipy.optimize import minimize_scalar
    zeilen, C, _ = kids_scherung(**auswahl)
    xi, th, S8 = kids_scherung_xi(tuple(punkt), gas, **auswahl)
    stuecke, R = (("lin", "rest"), teile) if teile is not None else (("ganz",), (1.0,))
    d = np.array([z[4] for z in zeilen])
    Ci = np.linalg.inv(C)
    def stueck(A, s, nur_faltung=False):                           # Stück s bei der Ausrichtung A
        m = []
        for vz, i, j, t, _ in zeilen:
            k, e = th.index(t), 0 if vz == "+" else 1
            g, gi, ii = xi[(min(i, j), max(i, j))][s]
            m.append(g[e][k] + (0 if nur_faltung else A * gi[e][k] + A * A * ii[e][k]))
        return np.array(m)
    def loesung(A):                                                # freie Amplituden bei festem A linear gelöst
        b = [stueck(A, s) for s in stuecke]
        frei = [n for n in range(len(R)) if R[n] is None]
        rest = d - sum(R[n] * b[n] for n in range(len(R)) if R[n] is not None)
        F = np.array([[b[n] @ Ci @ b[m] for m in frei] for n in frei])
        x = np.linalg.solve(F, [b[n] @ Ci @ rest for n in frei]) if frei else np.zeros(0)
        r = rest - sum(x[i] * b[n] for i, n in enumerate(frei))
        return float(r @ Ci @ r), frei, x, np.linalg.inv(F) if frei else F
    o = minimize_scalar(lambda A: loesung(A)[0], bounds=(-3, 6), method="bounded")
    chi, frei, x, Cx = loesung(o.x)
    aus = (chi, float(o.x), float(S8), len(d))
    if teile is not None:
        Rw = [float(x[frei.index(n)]) if n in frei else R[n] for n in (0, 1)]
        e = [math.sqrt(Cx[frei.index(n), frei.index(n)]) if n in frei else 0.0 for n in (0, 1)]
        kor = float(Cx[0, 1] / math.sqrt(Cx[0, 0] * Cx[1, 1])) if len(frei) == 2 else 0.0
        return aus + (Rw[0], Rw[1], e[0], e[1], kor)
    if verhaeltnis:
        g = stueck(o.x, "ganz", nur_faltung=True)
        F = float(g @ Ci @ g)
        return aus + (float(g @ Ci @ (d - stueck(o.x, "ganz") + g)) / F, 1 / math.sqrt(F))
    return aus


def kids_scherung_frei(punkt, **auswahl):
    """KiDS-1000-Scherung mit freier Amplitude ln10¹⁰A_s um den Punkt (sonst wie kids_scherung_chi2): (S8, χ², A)."""
    from scipy.optimize import minimize_scalar
    f = lambda d: kids_scherung_chi2((punkt[0], punkt[1], punkt[2], punkt[3] + d) + tuple(punkt[4:]), **auswahl)
    r = minimize_scalar(lambda d: f(d)[0], bounds=(-0.5, 0.2), method="bounded", options=dict(xatol=1e-3))
    c, A, S8, _ = f(r.x)
    return S8, c, A


def hiflugcs_gang(innen=300.0, aussen=1000.0):
    """HIFLUGCS: Verhältnis g_tot/g_bar innen (r < innen [kpc]) und außen (r > aussen) je Haufen. Gibt (Median innen, Median
    außen, Zahl der Haufen mit Punkten innen und außen, davon innen höher)."""
    h = hiflugcs()
    ri, ra, n, hoch = [], [], 0, 0
    for name in sorted({x[0] for x in h}):
        P = [x for x in h if x[0] == name]
        i = [10**(x[3] - x[2]) for x in P if x[1] < innen]
        o = [10**(x[3] - x[2]) for x in P if x[1] > aussen]
        ri += i
        ra += o
        if i and o:
            n += 1
            hoch += int(np.median(i) > np.median(o))
    return float(np.median(ri)), float(np.median(ra)), n, hoch


def xcop_staub(c=None, von=0.05, ohne=()):
    """X-COP von von·R₅₀₀ bis R₂₀₀: die hydrostatische Masse gegen Baryonen plus Mode, beide an R₂₀₀ mit dem Maß des Ganzen
    geeicht, M(R₂₀₀) = (1 + Q)·M_b(R₂₀₀). Die Baryonen liegen wie das Gas; die Mode mit dem Gewicht Q·M_b(R₂₀₀) entweder
    ebenso wie das Gas (c = None) oder wie druckfreier Staub nach dem Einfall, NFW-Form mit der Konzentration c = R₂₀₀/r_s.
    ohne: ausgelassene Haufen. Gibt (χ² in lg M mit dem halben Band als Fehler, Zahl der Punkte, mittlerer Rest innerhalb 0,15·R₅₀₀)."""
    import kern as K
    Q = K.q_je_radiant()
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    nfw = lambda x: math.log(1 + x) - x / (1 + x)
    chi, n, kern_rest = 0.0, 0, []
    for name, v in xcop().items():
        if name in ohne:
            continue
        R = v["R"] * K.kpc
        rho = 3 * v["M"] * K.M_sun / (4 * math.pi * R**3)
        i5, i2 = np.argmin(abs(np.log(rho / (500 * rc)))), np.argmin(abs(np.log(rho / (200 * rc))))
        Mb2 = v["M"][i2] / (1 + Q)                          # Baryonen an R₂₀₀ aus der Bilanz
        b = v["Mgas"] / v["Mgas"][i2] * Mb2                  # Baryonen wie das Gas
        if c is None:
            mode = Q * b
        else:
            mode = Q * Mb2 * np.array([nfw(c * r / R[i2]) for r in R]) / nfw(c)
        w = (R >= von * R[i5]) & (R <= R[i2])
        r = np.log10(v["M"][w] / (b[w] + mode[w]))
        e = np.log10(v["M_hi"][w] / v["M_lo"][w]) / 2
        chi += float(np.sum((r / e)**2))
        n += int(w.sum())
        kern_rest += list(r[R[w] <= 0.15 * R[i5]])
    return chi, n, float(np.mean(kern_rest))


def xcop_gas_faltung(k, c=4.0):
    """X-COP: Faltung eines Haufens, dessen Gas nie abkühlte, gegen die Mode als Staub, je k [1/Mpc]. Die Mode trägt
    Q·M_b(R₂₀₀) wie druckfreier Staub nach dem Einfall (NFW, Konzentration c), die Baryonen M_b(R₂₀₀) = M(R₂₀₀)/(1 + Q)
    liegen wie das Gas (innerhalb des ersten Radius mit fester Dichte), verglichen mit (1 + Q)·M_b(R₂₀₀) ganz wie Staub:
    (u/u_Staub)² der Fourier-Bilder bis R₂₀₀, Median der 12 Haufen."""
    import kern as K
    Q = K.q_je_radiant()
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    nfw = lambda x: np.log(1 + x) - x / (1 + x)
    k = np.atleast_1d(np.asarray(k, float)) / K.Mpc
    aus = []
    for v in xcop().values():
        R = v["R"] * K.kpc
        rho = 3 * v["M"] * K.M_sun / (4 * math.pi * R**3)
        i2 = np.argmin(abs(np.log(rho / (200 * rc))))
        r = np.geomspace(1e-3 * R[i2], R[i2], 6000)
        g = np.where(r < R[0], v["Mgas"][0] * (r / R[0])**3, np.interp(np.log(r), np.log(R), v["Mgas"]))
        staub = nfw(c * r / R[i2]) / nfw(c)
        aus.append((K.profil_fourier(r, Q * staub + g / g[-1], k) / K.profil_fourier(r, staub, k))**2)
    return np.median(aus, axis=0)


def hiflugcs_staub(c=None, Q=None, ohne=()):
    """HIFLUGCS: g_tot gegen Baryonen (g_bar wie gemessen) plus Mode mit dem Gewicht Q·M_b(R₂₀₀). Die Mode liegt wie die
    Baryonen (c = None, dann (1 + Q)·g_bar) oder wie druckfreier Staub nach dem Einfall, NFW-Form mit c = R₂₀₀/r_s; R₂₀₀ aus
    g_tot (mittlere Dichte 200-mal die kritische, jenseits des letzten Punkts mit der Neigung der letzten zwei Punkte
    fortgesetzt); Q: das Gewicht der Mode (sonst 2π − 1); ohne: Namen ausgelassener Haufen. Gibt (χ² mit dem mittleren Fehler je Punkt, Zahl, mittlerer Rest
    innerhalb 300 kpc, außerhalb 1 Mpc)."""
    import kern as K
    Q = K.q_je_radiant() if Q is None else Q
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    nfw = lambda x: math.log(1 + x) - x / (1 + x)
    h = hiflugcs()
    chi, n, ri, ra = 0.0, 0, [], []
    for name in sorted({x[0] for x in h} - set(ohne)):
        P = sorted([x for x in h if x[0] == name], key=lambda x: x[1])
        r = np.array([x[1] for x in P]) * K.kpc
        Mt = 10**np.array([x[3] for x in P]) * r**2 / K.G
        Mb = 10**np.array([x[2] for x in P]) * r**2 / K.G
        lr, lt, lb = np.log(r), np.log(Mt), np.log(Mb)
        if len(r) > 1:
            st, sb = (lt[-1] - lt[-2]) / (lr[-1] - lr[-2]), (lb[-1] - lb[-2]) / (lr[-1] - lr[-2])
        else:
            st, sb = 1.0, 1.0
        x = np.linspace(lr[0], lr[-1] + 3, 3000)
        Mtx = np.exp(np.where(x <= lr[-1], np.interp(x, lr, lt), lt[-1] + st * (x - lr[-1])))
        Mbx = np.exp(np.where(x <= lr[-1], np.interp(x, lr, lb), lb[-1] + sb * (x - lr[-1])))
        i2 = np.argmin(abs(np.log(3 * Mtx / (4 * math.pi * np.exp(3 * x)) / (200 * rc))))
        R2, Mb2 = math.exp(x[i2]), Mbx[i2]
        mode = Q * Mb if c is None else Q * Mb2 * np.array([nfw(c * ri_ / R2) for ri_ in r]) / nfw(c)
        res = np.log10(Mt / (Mb + mode))
        e = np.array([(x_[4] + x_[5]) / 2 for x_ in P])
        chi += float(np.sum((res / e)**2))
        n += len(r)
        ri += list(res[r < 300 * K.kpc])
        ra += list(res[r > 1000 * K.kpc])
    return chi, n, float(np.mean(ri)), float(np.mean(ra))


def haufen_q(ohne=(), c=None):
    """Q und Konzentration c, die HIFLUGCS mit der Mode in der Form nach dem Einfall liest (beide frei, hiflugcs_staub);
    ohne: ausgelassene Haufen; c: Konzentration fest, nur Q frei. Gibt (Q, c, χ²)."""
    from scipy.optimize import minimize, minimize_scalar
    if c is not None:
        r = minimize_scalar(lambda q: hiflugcs_staub(c, Q=q, ohne=ohne)[0], bounds=(3.0, 9.0), method="bounded",
                            options=dict(xatol=1e-5))
        return float(r.x), float(c), float(r.fun)
    r = minimize(lambda p: hiflugcs_staub(p[0], Q=p[1], ohne=ohne)[0], [5.0, 5.5], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-4))
    return float(r.x[1]), float(r.x[0]), float(r.fun)


def haufen_q_fehler(c=None):
    """Fehler von Q aus HIFLUGCS durch Auslassen je eines Haufens (Jackknife über die 10 Haufen); c wie bei haufen_q."""
    namen = sorted({x[0] for x in hiflugcs()})
    q = np.array([haufen_q((n,), c)[0] for n in namen])
    return float(math.sqrt((len(q) - 1) / len(q) * np.sum((q - q.mean())**2)))


def haufen_c(welche, ohne=()):
    """Konzentration c der Mode in der Form nach dem Einfall, die X-COP (welche = "xcop", xcop_staub) oder HIFLUGCS
    ("hiflugcs", hiflugcs_staub) liest, kleinstes χ² über c, ohne: ausgelassene Haufen."""
    from scipy.optimize import minimize_scalar
    f = {"xcop": xcop_staub, "hiflugcs": hiflugcs_staub}[welche]
    return float(minimize_scalar(lambda c: f(c, ohne=ohne)[0], bounds=(0.5, 30), method="bounded").x)


def haufen_c_fehler(welche):
    """Fehler der Konzentration aus dem Auslassen je eines Haufens (Jackknife), welche wie bei haufen_c."""
    namen = sorted(xcop()) if welche == "xcop" else sorted({x[0] for x in hiflugcs()})
    c = np.array([haufen_c(welche, (n,)) for n in namen])
    return float(math.sqrt((len(c) - 1) / len(c) * np.sum((c - c.mean())**2)))


def xcop_form(von=0.05, bis=0.15):
    """X-COP: Median von M/M_Gas über alle Radien zwischen von·R₅₀₀ und bis·R₅₀₀ aller Haufen, geteilt durch 1 + Q."""
    import kern as K
    Q = K.q_je_radiant()
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    w = []
    for v in xcop().values():
        R = v["R"] * K.kpc
        i5 = np.argmin(abs(np.log(3 * v["M"] * K.M_sun / (4 * math.pi * R**3) / (500 * rc))))
        x = v["R"] / v["R"][i5]
        w += list((v["M"] / v["Mgas"] / (1 + Q))[(x >= von) & (x < bis)])
    return float(np.median(w))


def xcop_sterne(ort=0.1):
    """X-COP bei ort·R₅₀₀: das Gas wirkt mit dem Maß des Ganzen (1 + Q)·g_Gas, den Rest der hydrostatischen Beschleunigung
    tragen Sterne mit der örtlichen Regel, ν(g_*)·g_* = g − (1 + Q)·g_Gas. Gibt (Radius [kpc], Sternmasse [M☉] und
    Sterne je Gas, je als Median, kleinster und größter Wert)."""
    import kern as K
    from scipy.optimize import brentq
    Q = K.q_je_radiant()
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    r_, m_, f_ = [], [], []
    for v in xcop().values():
        R = v["R"] * K.kpc
        i5 = np.argmin(abs(np.log(3 * v["M"] * K.M_sun / (4 * math.pi * R**3) / (500 * rc))))
        r = ort * R[i5]
        M, Mg = np.interp(r, R, v["M"]) * K.M_sun, np.interp(r, R, v["Mgas"]) * K.M_sun
        rest = K.G * (M - (1 + Q) * Mg) / r**2
        gs = brentq(lambda x: K.verstaerkung(x) * x - rest, 1e-20, rest) if rest > 0 else 0.0
        r_.append(r / K.kpc)
        m_.append(gs * r**2 / K.G / K.M_sun)
        f_.append(m_[-1] * K.M_sun / Mg)
    drei = lambda x: (float(np.median(x)), float(min(x)), float(max(x)))
    return drei(r_), drei(m_), drei(f_)


def kids_heiss(typ, lg_von=-20.0, lg_bis=0.0, einfall=None, fest=None):
    """Linsen eines Typs (kids_typ), der Überschuss über die Regel auf Sterne und kaltes Gas gelesen als nie abgekühltes Gas.
    Ohne einfall mit dem Maß des Ganzen am Ort: M_heiß/M_b = (g_obs − ν(g_bar)·g_bar)/((1 + Q)·g_bar) je Bin zwischen lg_von
    und lg_bis, gibt (Zahl, Median, kleinster, größter Wert). einfall = (lg M_b, c): das Gas mit seiner nicht übergebenen
    Mode, zusammen M₂₀₀ = (1 + Q)·M_heiß, liegt als druckfreier Staub in der Form nach dem Einfall (NFW mit Konzentration c,
    R₂₀₀ bei 200·ρ_krit heute) um die Regel der Baryonen M_b der Messung als Punktmasse, umgerechnet über ΔΣ wie die Messung
    (kern.linse_umgerechnet), je Bin am Radius R = √(G·M_b/g_bar), M_heiß über die Bins angepasst. Gibt dann (M_heiß/M_b, χ²).
    fest: M_heiß/M_b vorgegeben statt angepasst."""
    import kern as K
    Q = K.q_je_radiant()
    gb, go, ge = kids_typ(typ)
    w = (np.log10(gb) >= lg_von) & (np.log10(gb) <= lg_bis)
    if einfall is None:
        f = (go[w] - np.array([K.verstaerkung(x) * x for x in gb[w]])) / ((1 + Q) * gb[w])
        return int(w.sum()), float(np.median(f)), float(f.min()), float(f.max())
    from scipy.optimize import minimize_scalar
    lg_Mb, c = einfall
    gb, go, e = gb[w], go[w], ge[w] / go[w] / math.log(10)
    M = 10**lg_Mb * K.M_sun
    r = np.geomspace(0.01 * K.kpc, 1e6 * K.kpc, 6000)
    regel = np.array([K.verstaerkung(x) * x for x in K.G * M / r**2])
    R = np.sqrt(K.G * M / gb)
    rc = 3 * K.H0**2 / (8 * math.pi * K.G)
    nfw = lambda y: np.log(1 + y) - y / (1 + y)

    def chi(lf):
        M200 = (1 + Q) * 10**lf * M
        R200 = (3 * M200 / (4 * math.pi * 200 * rc))**(1 / 3)
        p = K.linse_umgerechnet(r, regel + K.G * M200 * nfw(c * r / R200) / nfw(c) / r**2, R)
        return float(np.sum((np.log10(go / p) / e)**2))
    if fest is not None:
        return float(fest), chi(math.log10(fest))
    b = minimize_scalar(chi, bounds=(-1.5, 2.5), method="bounded", options=dict(xatol=1e-4))
    return float(10**b.x), float(b.fun)


def xcop_abkuehlung(alter_gyr=None, X=0.75, ort=1.0):
    """X-COP: Temperatur des Gases aus dem hydrostatischen Gleichgewicht (Druck von außen nach innen integriert, null am letzten
    Punkt), Abkühlzeit durch Bremsstrahlung t = (3/2)(n_e + n_i)kT/(n_e·n_i·Λ) bei Wasserstoffanteil X. Je Haufen der Massenanteil
    des Gases innerhalb R₅₀₀, das sich in alter_gyr [Gyr] hätte abkühlen können (Vorgabe das Weltalter im Atmen des
    Formteils), Temperatur [keV] und Abkühlzeit [Gyr] bei ort·R₅₀₀. Gibt (Median und größter Anteil, kleinste und größte
    Temperatur, kürzeste und längste Abkühlzeit)."""
    import kern as K
    alter_gyr = K.weltalter(lambda a: K.dichte_dunkle_energie()) / K.Gyr if alter_gyr is None else alter_gyr
    mu, mue = 1 / (2 * X + 0.75 * (1 - X)), 2 / (1 + X)
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    anteil, temp, zeit = [], [], []
    for v in xcop().values():
        R, Mg, M = v["R"] * K.kpc, v["Mgas"] * K.M_sun, v["M"] * K.M_sun
        i5 = np.argmin(abs(np.log(3 * M / (4 * math.pi * R**3) / (500 * rc))))
        rho = np.gradient(Mg, R) / (4 * math.pi * R**2)
        g = K.G * M / R**2
        P = np.zeros_like(R)
        for i in range(len(R) - 2, -1, -1):
            P[i] = P[i + 1] + 0.5 * (rho[i] * g[i] + rho[i + 1] * g[i + 1]) * (R[i + 1] - R[i])
        T = np.maximum(P * mu * K.m_p / (rho * K.k_B), 1.0)
        ne, n = rho / (mue * K.m_p), rho / (mu * K.m_p)
        tc = 1.5 * n * K.k_B * T / (ne * (n - ne) * K.BREMSSTRAHLUNG * np.sqrt(T)) / K.Gyr
        dM = np.gradient(Mg, R) * np.gradient(R)
        innen = np.arange(len(R)) <= i5
        anteil.append(float(np.sum(dM[innen & (tc < alter_gyr)]) / np.sum(dM[innen])))
        io = np.argmin(abs(R - ort * R[i5]))
        temp.append(float(T[io] * K.k_B / (1e3 * K.e)))
        zeit.append(float(tc[io]))
    return float(np.median(anteil)), max(anteil), min(temp), max(temp), min(zeit), max(zeit)


def xcop_bilanz():
    """X-COP je Haufen bei R₅₀₀ und R₂₀₀ (mittlere Dichte 500 und 200 mal die kritische heute): M/M_Gas und das Verhältnis
    des hydrostatischen Felds zur Regel auf das Gas allein; dazu dasselbe Verhältnis bei 0,1·R₅₀₀ und die Tiefe
    δ = GM/(R·c²) bei R₂₀₀."""
    import kern as K
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    aus = {}
    for n, v in xcop().items():
        R = v["R"] * K.kpc
        rho = 3 * v["M"] * K.M_sun / (4 * math.pi * R**3)
        gN, g = K.G * v["Mgas"] * K.M_sun / R**2, K.G * v["M"] * K.M_sun / R**2
        verh = g / np.array([K.verstaerkung(x) * x for x in gN])
        i5, i2 = np.argmin(abs(np.log(rho / (500 * rc)))), np.argmin(abs(np.log(rho / (200 * rc))))
        ik = np.argmin(abs(v["R"] - 0.1 * v["R"][i5]))
        def verlangt(i):                                      # Baryonen, die die Regel für das hydrostatische Feld verlangt
            from scipy.optimize import brentq
            x = brentq(lambda gn: K.verstaerkung(gn) * gn - g[i], 1e-16, g[i])
            return x * R[i]**2 / K.G / K.M_sun / v["Mgas"][i]
        aus[n] = dict(M_Mgas_500=v["M"][i5] / v["Mgas"][i5], M_Mgas_200=v["M"][i2] / v["Mgas"][i2],
                      regel_500=verh[i5], regel_200=verh[i2], regel_kern=verh[ik], aN_200=gN[i2] / K.a0,
                      verlangt_200=verlangt(i2), verlangt_kern=verlangt(ik), tiefe_200=g[i2] * R[i2] / K.c**2)
    return aus


PLANCK_START = dict(wb=0.02237, wc=0.1200, H0=67.36, a0=1.06, lnAs=3.044, ns=0.9649, tau=0.0544, A=1.0, w0=-1.0, wa=0.0,
                    Q=2 * math.pi - 1, Fs=0.017, k=10.0, a_s=0.77, F0=0.0097)
PLANCK_SCHRITT = dict(wb=0.00015, wc=0.0012, H0=0.5, a0=0.01, lnAs=0.014, ns=0.004, tau=0.007, A=0.0025, w0=0.05, wa=0.2,
                      Q=0.05, Fs=0.003, k=2.0, a_s=0.03, F0=0.003)


def _planck_punkt(p, Q):
    """Die sechs Zahlen für CAMB aus dem Punkt p: ΛCDM (Q = None) mit ω_c und H₀, sonst ω_c = Q·ω_b, H₀ aus a₀;
    Q = "lesen": Q aus dem Punkt."""
    import kern as K
    if Q is None:
        return p["wb"], p["wc"], p["H0"]
    q = p["Q"] if Q == "lesen" else Q
    return p["wb"], q * p["wb"], K.atmen_heute(p["a0"] * 1e-10, Q=q, wb=p["wb"], wnu=0.06 / 93.14)[0]


_CAMB = {}


def _de_gas(p):
    """Dunkle Energie mit dem Zustandsteil im kalten Gas (kern.zustand_gas) für den Punkt p (Fs, k, form, a_s): de(z) =
    ρ_DE(z)/ρ_DE(0) und die Tafel (a, w) für CAMB."""
    import kern as K
    je = K.zustandsteil_je_gas(p.get("Q", K.q_je_radiant())) / 3
    r = lambda a: 1 + je * K.zustand_gas(a, p["Fs"], p["k"], p["a_s"])
    a = np.concatenate([np.geomspace(1e-5, 0.05, 60), np.linspace(0.05, 1.0, 400)[1:]])
    d = 1e-4
    w = -1 - (np.log(r(a * (1 + d))) - np.log(r(a * (1 - d)))) / (3 * math.log((1 + d) / (1 - d)))
    return (lambda z: r(1 / (1 + z)) / r(1.0)), (a, w)


def planck_reste(p, Q=None):
    """Geweißte Reste des Punkts p je Datensatz: Spektren (plik_lite mit TT unter ℓ = 30), Linse (PR4), DESI DR2,
    Pantheon+ (Nullpunkt frei), Wachstum f·σ₈, Urknall-Nukleosynthese und die Vorgaben für τ und die Eichung; Abstände
    mit Ω_m aus ω_b, ω_c und der Neutrinomasse, r_d aus CAMB; dazu σ₈ (kein Rest). Die Dunkle Energie heute setzt a₀
    (oder H₀), w₀ und w_a geben ihren Gang."""
    import kern as K
    wb, wc, H = _planck_punkt(p, Q)
    w0, wa = p.get("w0", -1.0), p.get("wa", 0.0)
    de, tafel = None, None
    if "form" in p:                                                  # Zustandsteil im kalten Gas, Form des Modells
        de, tafel = _de_gas(p)
    kenn = (wb, wc, H, p["lnAs"], p["ns"], p["tau"], w0, wa) + ((p["Fs"], p["k"], p["form"]) if de else ())
    if kenn not in _CAMB:
        _CAMB.clear()
        _CAMB[kenn] = K.cmb_spektren(*kenn[:6], alles=True, wachstum=wachstum()[0], w0=w0, wa=wa, w_tafel=tafel)
    tt, te, ee, pp, rd, fs8, s8 = _CAMB[kenn]
    _, v8, L8 = wachstum()
    h2 = (H / 100)**2
    om, or_ = (wb + wc + 0.06 / 93.14) / h2, K.dichte_heute()[0] / K.rho_krit100 / h2
    return dict(spektren=planck_rest(tt, te, ee, p["A"]), linse=planck_linse_rest(pp),
                bao=bao_absolut(H, om, or_, rd, w0, wa, rest=True, de=de), sn=pantheon_form(om, or_, w0, wa, rest=True, de=de),
                wachstum=np.linalg.solve(L8, v8 - fs8), bbn=bbn_rest(wb), sigma8=s8,
                vor=np.array([(p["tau"] - K.TAU_PLANCK[0]) / K.TAU_PLANCK[1], (p["A"] - K.A_PLANCK[0]) / K.A_PLANCK[1]]))


def planck_anpassung(Q=None, mit=("spektren",), fest=(), zustand=False, neu=False, start=None, form=None):
    """Bester Punkt an den Datensätzen mit (aus planck_reste; die Vorgaben für τ und Eichung immer). Q = None: ΛCDM, frei
    ω_b, ω_c, H₀, A_s, n_s, τ; sonst ω_c = Q·ω_b und H₀ aus a₀ (Formteil, flach), frei ω_b, a₀, A_s, n_s, τ, mit
    Q = "lesen" dazu Q; fest: Namen,
    die auf dem Startwert bleiben; zustand: dazu w₀ und w_a frei; zustand="saettigung": nur w_a frei, das Maximum der
    Dunklen Energie liegt dort, wo drei Teile die Faltung einholen (w₀ = −1 − w_a·(1 − a_s)). Abgelegt in
    zwischen/etappe_cmb.json, gerechnet nur, wenn er fehlt oder neu. Gibt dict(chi2 je Datensatz und gesamt, punkt,
    fehler, kovarianz der freien Zahlen); start: Startwerte einzelner freier Zahlen (sonst PLANCK_START), mit im Namen;
    form: a₀ des Formteils [10⁻¹⁰ m/s²] für die Sättigung (drei Teile des Formteils holen die Faltung ein), a₀ des Punkts
    setzt dann die ganze Dunkle Energie heute, Formteil und Zustandsteil; zustand="gas": der Zustandsteil im kalten Gas
    (kern.zustand_gas, frei F_s und k, Maximum an der Sättigung des Formteils, form sonst a₀ des Kerns), a₀ des Punkts folgt
    aus Formteil und Gas heute; zustand="gas_frei": dazu die Lage a_s des Maximums frei; zustand="gas_bera": das Gas im Kreis
    durch seinen heutigen Wert F₀ gelegt (kern.gas_kreis), frei F₀ und k ≥ 0, das Verhältnis bei z = 0,34 zu heute (Bera)
    als Rest, Formteil form; zustand="gezaehlt": F₀ das gezählte Gas (gas_gezaehlt) oder der Startwert F0, frei k ≥ 0 mit
    demselben Rest und a₀ als Formteil (mit fest=("a0",) der Standardwert); zustand="gas_frueh": wie "gas_bera", dazu die
    Lage a_s des Maximums frei und das Verhältnis bei z = 1,04 zu heute (Chowdhury, hi_frueh_eins) als zweiter Rest;
    zustand="gezaehlt_frueh": wie "gezaehlt" mit dem Startwert F0 und a₀ als Formteil, dazu wie "gas_frueh" die Lage a_s frei
    und beide Reste des Gases."""
    import huelle_ablage as HA
    from scipy.optimize import least_squares
    name = f"Q={'frei' if Q is None else Q if isinstance(Q, str) else round(Q, 6)} mit={'+'.join(mit)} fest={','.join(fest)}" + (
        f" {zustand}" if zustand in ("saettigung", "gas", "gas_frei", "gas_bera", "gezaehlt", "gas_frueh", "gezaehlt_frueh")
        else " zustand" if zustand else "") + (
        f" start={sorted(start.items())}" if start else "") + (f" form={form}" if form else "")
    alle = HA.lesen()
    if not neu and name in alle.get("_cmb", {}):
        return alle["_cmb"][name]
    kreis = zustand in ("gas_bera", "gezaehlt", "gas_frueh", "gezaehlt_frueh")
    frueh = zustand in ("gas_frueh", "gezaehlt_frueh")
    st = dict(dict(PLANCK_START, k=3.0) if kreis else PLANCK_START, **{k: v for k, v in (start or {}).items() if k in fest})
    frei = [k for k in (("wb", "wc", "H0") if Q is None else ("wb",) + (("a0",) if zustand not in ("gas", "gas_bera", "gas_frueh") else ())
                        + (("Q",) if Q == "lesen" else ()))
            + ("lnAs", "ns", "tau", "A")
            + (("wa",) if zustand == "saettigung" else ("Fs", "k") if zustand == "gas" else ("Fs", "k", "a_s")
               if zustand == "gas_frei" else ("F0", "k") if zustand == "gas_bera" else ("F0", "k", "a_s")
               if zustand == "gas_frueh" else ("k", "a_s") if zustand == "gezaehlt_frueh" else ("k",) if zustand == "gezaehlt"
               else ("w0", "wa") if zustand else ())
            if k not in fest]
    import kern as K
    gas = zustand in ("gas", "gas_frei", "gas_bera", "gas_frueh")
    je = K.zustandsteil_je_gas(Q) / K.B1_ZELLE if kreis else None
    form = (form or round(K.a0 * 1e10, 4)) if gas else form

    def punkt(y):
        p = dict(st, **{k: st[k] + PLANCK_SCHRITT[k] * v for k, v in zip(frei, y)})
        if zustand == "saettigung":                                  # drei Teile holen die Faltung ein, im Punkt selbst
            q = p["Q"] if Q == "lesen" else Q
            a_s = (((1 + q) * p["wb"] + 0.06 / 93.14) * K.rho_krit100
                   / K.dichte_dunkle_energie((form or p["a0"]) * 1e-10))**(1 / 3)
            p["w0"] = -1 - p["wa"] * (1 - a_s)
        if kreis:                                                    # Gas im Kreis, durch heute gelegt
            gez = zustand in ("gezaehlt", "gezaehlt_frueh")
            p["form"] = p["a0"] if gez else form
            p["F0"] = (start or {}).get("F0", gas_gezaehlt()[0]) if gez else p["F0"]
            if not frueh:
                p["a_s"] = (((1 + Q) * p["wb"] + 0.06 / 93.14) * K.rho_krit100
                            / K.dichte_dunkle_energie(p["form"] * 1e-10))**(1 / 3)
            u = 1 / p["a_s"]
            p["Fs"] = p["F0"] / (u**p["k"] * math.exp(p["k"] * (1 - u)))
            p["a0"] = p["form"] * math.sqrt(1 + je * p["F0"])
        elif gas:                                                    # Maximum des Gases an der Sättigung des Formteils
            q = p["Q"] if Q == "lesen" else Q
            if zustand == "gas":
                p["a_s"] = (((1 + q) * p["wb"] + 0.06 / 93.14) * K.rho_krit100
                            / K.dichte_dunkle_energie(form * 1e-10))**(1 / 3)
            p["form"] = form
            p["a0"] = form * math.sqrt(1 + K.zustandsteil_je_gas(q) / 3 * K.zustand_gas(1.0, p["Fs"], p["k"], p["a_s"]))
        return p

    def rest(y):
        p = punkt(y)
        r = planck_reste(p, Q)
        aus = np.concatenate([r[k] for k in mit + ("vor",)])
        return np.append(aus, [_bera_rest(p)] + ([_frueh_rest(p)] if frueh else [])) if kreis else aus

    y0 = np.array([((start or {}).get(k, st[k]) - st[k]) / PLANCK_SCHRITT[k] for k in frei])
    unten = np.array([-st[k] / PLANCK_SCHRITT[k] if kreis and k == "k" else (0.05 - st[k]) / PLANCK_SCHRITT[k]
                      if k == "a_s" else -np.inf for k in frei])
    erg = least_squares(rest, y0, bounds=(unten, np.full(len(frei), np.inf)), diff_step=0.05, x_scale=1.0, xtol=1e-6,
                        ftol=1e-8)
    s = np.array([PLANCK_SCHRITT[k] for k in frei])
    kov = np.linalg.pinv(erg.jac.T @ erg.jac) * np.outer(s, s)          # Pseudoinverse: gleich der Inversen, wenn nicht entartet
    p = punkt(erg.x)
    p["wb"], p["wc"], p["H0"] = _planck_punkt(p, Q)
    r = planck_reste(p, Q)
    sigma8 = float(r.pop("sigma8"))
    chi2 = {k: float(v @ v) for k, v in r.items()}
    chi2["gesamt"] = sum(chi2[k] for k in mit + ("vor",))
    if kreis:
        chi2["bera"] = float(_bera_rest(p)**2)
    if frueh:
        chi2["frueh"] = float(_frueh_rest(p)**2)
    aus = dict(chi2=chi2, sigma8=sigma8, punkt=p, fehler={k: float(math.sqrt(kov[i, i])) for i, k in enumerate(frei)},
               kovarianz=dict(namen=frei, werte=kov.tolist()))
    alle = HA.lesen()
    alle.setdefault("_cmb", {})[name] = aus
    HA.schreiben(alle)
    return aus


def bao_absolut(H0, om, or_, rd, w0=-1.0, wa=0.0, rest=False, de=None, zmax=None):
    """χ² von DESI DR2 mit fester Skala: H₀ [km/s/Mpc], Ω_m, Ω_r, flach, Schallhorizont r_d [Mpc] fest;
    Dunkle Energie mit w(a) = w₀ + w_a(1 − a) oder de(z) = ρ_DE(z)/ρ_DE(0). rest: der geweißte Rest statt χ².
    zmax: H₀ frei, aus den Werten bis zmax nach kleinsten Quadraten gelesen, gibt (H₀, Fehler), H0 bleibt unbenutzt."""
    z, v, g, C = desi()
    x = np.linspace(0, 2.4, 24001)
    E = np.sqrt(or_ * (1 + x)**4 + om * (1 + x)**3 + (1 - om - or_) * _de_gang(x, w0, wa, de))
    I = np.concatenate([[0.0], np.cumsum((1 / E[1:] + 1 / E[:-1]) / 2 * np.diff(x))])
    dh0 = 299792.458 / (1.0 if zmax is not None else H0)
    dm, dh = dh0 * np.interp(z, x, I), dh0 / np.interp(z, x, E)
    th = np.array([{"DM_over_rs": a, "DH_over_rs": b, "DV_over_rs": (zz * a * a * b)**(1 / 3)}[gg]
                   for zz, a, b, gg in zip(z, dm, dh, g)]) / rd
    if zmax is not None:                                          # Wert = th/H₀, gelesen u = 1/H₀
        s = z <= zmax
        Ci = np.linalg.inv(C[np.ix_(s, s)])
        A = th[s] @ Ci @ th[s]
        u = th[s] @ Ci @ v[s] / A
        return float(1 / u), float(1 / (u * u * math.sqrt(A)))
    r = v - th
    if rest:
        import scipy.linalg
        return scipy.linalg.solve_triangular(np.linalg.cholesky(C), r, lower=True)
    return float(r @ np.linalg.inv(C) @ r)


def _de_gang(z, w0=-1.0, wa=0.0, de=None):
    """ρ_DE(z)/ρ_DE(0): mit w(a) = w₀ + w_a(1 − a) oder aus de(z)."""
    if de is not None:
        return de(z)
    return (1 + z)**(3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))


def pantheon_form(om, or_, w0=-1.0, wa=0.0, rest=False, de=None, leerraum=None):
    """χ² der Pantheon+-Form (Nullpunkt frei, flach) bei Ω_m, Ω_r; Dunkle Energie mit w(a) = w₀ + w_a(1 − a) oder
    de(z) = ρ_DE(z)/ρ_DE(0). rest: der geweißte Rest statt χ² (Nullpunkt herausprojiziert). leerraum: (z_rand, ε), ein
    Leerraum um den Beobachter, innen mit H·(1 + ε) bis zum Rand bei z_rand, außen das Atmen. Jenseits des Randes ist
    1 + z = (1 + z_rand)·(1 + z')/(1 + z_H), z_H die Rotverschiebung, die das Atmen dem Rand gäbe, z' die der Quelle."""
    zs, zh, m, C = pantheon()
    Ci = _pantheon_inv()
    g = np.linspace(0, 3.0, 30001) if leerraum else np.linspace(0, 2.4, 24001)
    f = 1 / np.sqrt(or_ * (1 + g)**4 + om * (1 + g)**3 + (1 - om - or_) * _de_gang(g, w0, wa, de))
    I = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
    d = np.interp(zs, g, I)
    if leerraum:
        zr, e = leerraum
        zH = np.interp(np.interp(zr, g, I) / (1 + e), I, g)
        d = np.where(zs <= zr, d / (1 + e), np.interp((1 + zs) * (1 + zH) / (1 + zr) - 1, g, I))
    r = m - 5 * np.log10((1 + zh) * d)
    if rest:
        return _pantheon_wurzel() @ r
    eins = np.ones_like(r)
    return float(r @ Ci @ r - (eins @ Ci @ r)**2 / (eins @ Ci @ eins))


@lru_cache(maxsize=None)
def _pantheon_inv():
    return np.linalg.inv(pantheon()[3])


@lru_cache(maxsize=None)
def _pantheon_wurzel():
    """M mit MᵀM = C⁻¹ − C⁻¹11ᵀC⁻¹/(1ᵀC⁻¹1): weißt den Rest bei freiem Nullpunkt."""
    Ci = _pantheon_inv()
    u = Ci @ np.ones(len(Ci))
    lam, V = np.linalg.eigh(Ci - np.outer(u, u) / u.sum())
    return np.sqrt(np.clip(lam, 0, None))[:, None] * V.T


@lru_cache(maxsize=None)
def kernbildung_wb():
    """ω_b, den die Kernbildung allein liest (D/H und Y_p, bbn_rest): Bestwert und 1σ aus Δχ² = 1."""
    from scipy.optimize import brentq, minimize_scalar
    f = lambda w: float(bbn_rest(w) @ bbn_rest(w))
    r = minimize_scalar(f, bounds=(0.020, 0.024), method="bounded", options=dict(xatol=1e-8))
    lo = brentq(lambda w: f(w) - r.fun - 1, r.x - 0.002, r.x)
    hi = brentq(lambda w: f(w) - r.fun - 1, r.x, r.x + 0.002)
    return float(r.x), float(hi - lo) / 2


@lru_cache(maxsize=None)
def _schallwinkel_planck():
    """Schallwinkel θ_* in Plancks bestem Modell, der Bezug für die Lage der Spitzen."""
    import kern as K
    return K.schallwinkel(lambda a: K.dichte_de_planck())


def geometrie_reste(wb, a0, Q, w0=-1.0, wa=0.0, baryonen="kernbildung", spitzen=False):
    """Geweißte Reste der Geometrie des Ganzen bei ω_b, a₀ [10⁻¹⁰ m/s²] und Q: Lage der Spitzen (Schallwinkel θ_* gegen
    Plancks bestes Modell, 1σ nach Planck), Schallwelle (DESI DR2, feste Skala r_d aus ω_b und Q·ω_b), Form der Abstände
    (Pantheon+, Nullpunkt frei). a₀ setzt die Dunkle Energie heute, H₀ folgt aus ihr und (1 + Q)·ω_b, w₀ und w_a geben
    ihren Gang w = w₀ + w_a(1 − a), w = −1 ist der Formteil allein. Der offene Teil ω_b wird eigens gelesen:
    baryonen="kernbildung" aus D/H und Y_p, "spitzen" nach Planck, None gar nicht; spitzen: dazu ω_c = Q·ω_b gegen die
    Höhe der Spitzen (Planck)."""
    import kern as K
    rho0 = K.dichte_dunkle_energie(a0 * 1e-10)
    H0 = K.atmen_heute(a0 * 1e-10, Q=Q, wb=wb)[0]
    h2 = (H0 / 100)**2
    om, or_ = (1 + Q) * wb / h2, K.dichte_heute()[0] / K.rho_krit100 / h2
    th_ref = _schallwinkel_planck()
    alt = K.omega_b
    K.omega_b = wb
    try:
        th = K.schallwinkel(lambda x: rho0 * x**(-3 * (1 + w0 + wa)) * math.exp(-3 * wa * (1 - x)), Q)
        rd = K.schallhorizont(Q)
    finally:
        K.omega_b = alt
    r = [[(th / th_ref - 1) / K.theta_stern_rel], bao_absolut(H0, om, or_, rd, w0, wa, rest=True),
         pantheon_form(om, or_, w0, wa, rest=True)]
    if baryonen == "kernbildung":
        r.append(bbn_rest(wb))
    elif baryonen == "spitzen":
        r.append([(wb - K.omega_b) / K.omega_b_fehler])
    if spitzen:
        r.append([(Q * wb - K.omega_c) / K.omega_c_fehler])
    return np.concatenate(r)


@lru_cache(maxsize=None)
def geometrie_q(baryonen="kernbildung", saettigung=False, Q=None):
    """Q aus der Geometrie des Ganzen (geometrie_reste), frei ω_b, a₀ und Q, mit saettigung dazu w_a, das Maximum der
    Dunklen Energie an der Sättigung des Kreises (kern.w0_saettigung), sonst der Formteil. Q gegeben: Q fest, dann liest
    die Geometrie den offenen Teil ω_b. Gibt dict(Q, wb, wb_fehler, chi2), bei freiem Q dazu Q_fehler, wm = (1 + Q)·ω_b,
    wm_fehler, d_radiant (Δχ² von Q = 2π − 1) und d_spitzen (Δχ² des aus der Höhe der Spitzen gemessenen ω_c/ω_b)."""
    import kern as K
    from scipy.optimize import least_squares

    def anpassen(q):
        x0 = [0.0220, 1.085] + ([5.3] if q is None else []) + ([-0.4] if saettigung else [])
        s = [2e-4, 0.01] + ([0.05] if q is None else []) + ([0.2] if saettigung else [])

        def f(y):
            qq = y[2] if q is None else q
            wa = y[-1] if saettigung else 0.0
            w0 = K.w0_saettigung(wa, (1 + qq) * y[0], y[1] * 1e-10) if saettigung else -1.0
            return geometrie_reste(y[0], y[1], qq, w0, wa, baryonen)
        e = least_squares(f, x0, x_scale=s, diff_step=1e-3, xtol=1e-10, ftol=1e-12)
        return e, float(e.fun @ e.fun), np.linalg.inv(e.jac.T @ e.jac)

    e, chi, C = anpassen(Q)
    q = float(e.x[2]) if Q is None else Q
    aus = dict(Q=q, wb=float(e.x[0]), wb_fehler=float(math.sqrt(C[0, 0])), chi2=chi)
    if Q is None:
        g = np.zeros(len(e.x))
        g[0], g[2] = 1 + q, e.x[0]
        aus.update(Q_fehler=float(math.sqrt(C[2, 2])), wm=(1 + q) * float(e.x[0]), wm_fehler=float(math.sqrt(g @ C @ g)),
                   d_radiant=anpassen(K.q_je_radiant())[1] - chi, d_spitzen=anpassen(K.omega_c / K.omega_b)[1] - chi)
    return aus


@lru_cache(maxsize=None)
def ame2020_markiert():
    """Kerne aus AME2020 mit Kennung: Liste (Z, N, A, Bindung je Nukleon [keV], gemessen?); '#' = geschätzt."""
    aus = []
    for z in open(os.path.join(ORT, "elemente", "mass_1.mas20.txt"), encoding="latin-1"):
        if len(z) < 130 or not z[1:4].strip().lstrip("-").isdigit():
            continue
        try:
            N, Z, A = int(z[4:9]), int(z[9:14]), int(z[14:19])
            B = float(z[54:67].replace("#", "").strip())
        except ValueError:
            continue
        aus.append((Z, N, A, B, "#" not in z[28:67]))
    return aus


_EINHEIT = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15, "ps": 1e-12, "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0,
            "m": 60.0, "h": 3600.0, "d": 86400.0, "y": 3.15576e7, "ky": 3.15576e10, "My": 3.15576e13, "Gy": 3.15576e16,
            "Ty": 3.15576e19, "Py": 3.15576e22, "Ey": 3.15576e25, "Zy": 3.15576e28, "Yy": 3.15576e31}


@lru_cache(maxsize=None)
def nubase_spins():
    """NUBASE2020, Grundzustände mit eindeutig zugewiesenem Spin und Parität (keine Klammer, kein #, keine Wahl):
    {(Z, A): (2j, Parität ±1)}."""
    import re
    aus = {}
    for z in open(os.path.join(ORT, "elemente", "nubase_4.mas20.txt"), encoding="latin-1"):
        if z.startswith("#") or len(z) < 100 or z[7] != "0":
            continue
        jp = z[88:102].split()
        if not jp or "#" in jp[0] or "(" in jp[0]:
            continue
        t = jp[0].rstrip("*")                                      # * = unmittelbar gemessen
        m = re.fullmatch(r"(\d+)/2([+-])", t) or re.fullmatch(r"(\d+)([+-])", t)
        if not m:
            continue
        j2 = int(m.group(1)) if "/" in t else 2 * int(m.group(1))
        aus[(int(z[4:7]), int(z[0:3]))] = (j2, 1 if m.group(2) == "+" else -1)
    return aus


def nubase(stabil=False):
    """NUBASE2020, Grundzustände mit gemessener Halbwertszeit: {(Z, A): (T [s], {Zerfallsart: Anteil [%]})};
    stabil: die Menge der stabilen Grundzustände (Z, A)."""
    import re
    aus = {}
    for z in open(os.path.join(ORT, "elemente", "nubase_4.mas20.txt"), encoding="latin-1"):
        if z.startswith("#") or len(z) < 100 or z[7] != "0":
            continue
        A, Z = int(z[0:3]), int(z[4:7])
        t, u = z[69:78].strip(), z[78:80].strip()
        if stabil:
            if t == "stbl":
                aus[(Z, A)] = None
            continue
        if not t or t in ("stbl", "p-unst") or "#" in t or t[0] in "<>~" or u not in _EINHEIT:
            continue
        arten = {}
        for m in re.finditer(r"(?:^|;)([A-Za-z0-9+-]+)=([\d.eE+-]+)", z[119:].strip()):
            try:
                arten[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
        aus[(Z, A)] = (float(t) * _EINHEIT[u], arten)
    return set(aus) if stabil else aus


def sparc_zusatz_gipfel(kpc, a0=1.11e-10, ziehungen=1000, yd=YD):
    """Gipfel der Zusatzbeschleunigung g_obs − g_bar an SPARC: Mediane in Fenstern von 0,2 dex in g_bar/a₀, eine Parabel in
    lg(g_bar/a₀) durch die Fenster zwischen 0,5·a₀ und 8·a₀, Lage und Höhe ihres Gipfels in Einheiten von a₀; Streuung
    durch Ziehen der Galaxien. Gibt (Lage, Höhe, 16–84 % der Lage, 16–84 % der Höhe, Median in 4 bis 8·a₀)."""
    gb, go, gal, _ = _einfache_punkte(kpc, yd)
    x, ex = np.log10(gb / a0), (go - gb) / a0
    namen = np.unique(gal)
    kanten = np.arange(np.log10(0.5), np.log10(8.0) + 1e-9, 0.2)

    def gipfel(w):
        m = [(a + b) / 2 for a, b in zip(kanten[:-1], kanten[1:])]
        y = [np.median(ex[w & (x >= a) & (x < b)]) for a, b in zip(kanten[:-1], kanten[1:])]
        c = np.polyfit(m, y, 2)
        lx = -c[1] / (2 * c[0])
        return 10**lx, float(np.polyval(c, lx))
    lage, hoehe = gipfel(np.ones(len(x), bool))
    rng = np.random.default_rng(3)
    bs = []
    for _ in range(ziehungen):
        wahl = rng.choice(namen, len(namen))
        zahl = dict(zip(*np.unique(wahl, return_counts=True)))
        idx = np.concatenate([np.repeat(np.where(gal == n)[0], c) for n, c in zahl.items()])
        xs, es = x[idx], ex[idx]
        m = [(a + b) / 2 for a, b in zip(kanten[:-1], kanten[1:])]
        y = [np.median(es[(xs >= a) & (xs < b)]) for a, b in zip(kanten[:-1], kanten[1:])]
        c = np.polyfit(m, y, 2)
        lx = -c[1] / (2 * c[0])
        bs.append((10**lx, np.polyval(c, lx)))
    bs = np.array(bs)
    hoch = np.median(ex[(x >= np.log10(4)) & (x < np.log10(8))])
    return float(lage), float(hoehe), tuple(np.percentile(bs[:, 0], [16, 84])), tuple(np.percentile(bs[:, 1], [16, 84])), float(hoch)


def alfalfa(code=1, cz_max=15000.0):
    """ALFALFA α.100 (Haynes u. a. 2018): Quellen des Codes code (None: 1 und 2) bis cz_max: dict mit RA, Dec [°], cz [km/s],
    W50 [km/s], S [Jy km/s], D [Mpc, H0 = 70], lgM [lg M☉], AGC (Nummer), W20 [km/s], SNR, sigW [km/s]."""
    d = np.genfromtxt(os.path.join(ORT, "alfalfa", "a100.code12.table2.190808.csv"), delimiter=",", skip_header=1,
                      usecols=(2, 3, 6, 7, 10, 14, 16, 18, 0, 9, 12, 8))
    w = ((d[:, 7] == code) if code is not None else True) & (d[:, 2] <= cz_max) & (d[:, 2] > 0) & (d[:, 5] > 0)
    d = d[w]
    return dict(RA=d[:, 0], Dec=d[:, 1], cz=d[:, 2], W50=d[:, 3], S=d[:, 4], D=d[:, 5], lgM=d[:, 6], AGC=d[:, 8].astype(int),
                W20=d[:, 9], SNR=d[:, 10], sigW=d[:, 11])


def hi_streuung(snr_min=20.0, h_gas=2.0, rand=False):
    """Streuung des neutralen Gases je Quelle von ALFALFA α.100 (Code 1, SNR über snr_min) aus ihren zwei Breiten: W20/W50
    gibt über die Linie einer Scheibe an der Schwelle s = σ/(v_flach·sin i), daraus V = v_flach·sin i = W50/w50(s) und
    σ = s·V. rand: zweiter Weg von Hand, der Rand einer Stufe, mit einer Gaußglocke verbreitert, fällt von der Hälfte auf
    ein Fünftel über Φ⁻¹(0,8)·σ, also σ = (W20 − W50)/(2Φ⁻¹(0,8)) = (W20 − W50)/1,68, V = W50/2. Die Auflösung ist nicht
    abgezogen, sie verbreitert nur, ebenso Verwerfung und Bewegung quer zur Bahn: σ ist eine obere Grenze.
    Gibt (lg M_HI [M☉] bei H0 = 70, V [km/s], σ [km/s])."""
    import kern as K
    from scipy.stats import norm
    a = alfalfa()
    w = (a["SNR"] > snr_min) & (a["W20"] > a["W50"])
    W50, W20, lgM = a["W50"][w], a["W20"][w], a["lgM"][w]
    if rand:
        return lgM, W50 / 2, (W20 - W50) / (2 * norm.ppf(0.8))
    s, w5, w2 = K.hi_linie_tafel(h_gas)
    q = w2 / w5
    m = np.r_[True, q[1:] > np.maximum.accumulate(q)[:-1]]
    s, w5, q = s[m], w5[m], q[m]
    r = np.clip(W20 / W50, q[0], q[-1])
    sq = np.interp(r, q, s)
    V = W50 / np.interp(sq, s, w5)
    return lgM, V, sq * V


def hi_rand_gang(quelle="alfalfa", kanten=(30.0, 80.0, 150.0, 250.0, 400.0, 600.0), snr_min=20.0):
    """Median von W20 − W50 [km/s] je Stück in W50: ALFALFA α.100 (Code 1, SNR über snr_min) oder HIPASS (HICAT, Breiten
    an der Spitze, W20max − W50max, Spitzenfluss über snr_min-mal dem Rauschen). Gibt (Median je Stück, Zahl je Stück)."""
    if quelle == "alfalfa":
        a = alfalfa()
        W50, d, w = a["W50"], a["W20"] - a["W50"], a["SNR"] > snr_min
    else:
        h = hipass(("W50max", "W20max", "Speak", "RMS"))
        W50, d = h["W50max"], h["W20max"] - h["W50max"]
        w = (h["Speak"] > snr_min * h["RMS"]) & np.isfinite(d)
    k = np.asarray(kanten)
    sel = [w & (W50 >= lo) & (W50 < hi) for lo, hi in zip(k[:-1], k[1:])]
    return [float(np.median(d[s])) for s in sel], [int(s.sum()) for s in sel]


def hi_streuung_gang(kanten=(7.0, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0), snr_min=20.0, rand=False):
    """Median der Streuung des Gases (hi_streuung) je Stück in lg M_HI. Gibt (Mitten, Median σ [km/s], Zahl) und die
    Funktion σ(lg M_HI), dazwischen linear, außen der Rand."""
    lgM, _, sg = hi_streuung(snr_min, rand=rand)
    k = np.asarray(kanten)
    mi, me, n = [], [], []
    for lo, hi in zip(k[:-1], k[1:]):
        s = (lgM >= lo) & (lgM < hi)
        if s.sum() >= 10:
            mi.append(0.5 * (lo + hi))
            me.append(float(np.median(sg[s])))
            n.append(int(s.sum()))
    mi, me = np.array(mi), np.array(me)
    return mi, me, n, (lambda x: float(np.interp(x, mi, me)))


def alfalfa_umgebung(agc, R=2.0, dv=500.0, n=300, saat=1):
    """Umgebung einer Quelle in α.100 (Code 1 und 2): Zahl der Nachbarn im Abstand R [Mpc] auf dem Himmel (bei ihrem Abstand)
    und mit |Δcz| < dv, dazu ihr Rang unter n zufälligen Vergleichsquellen mit Abstand ±10 Mpc und lg M_HI ±0,5 (Anteil mit
    weniger Nachbarn, Gleichstand halb). Gibt (Nachbarn, Rang, Median der Vergleichsquellen)."""
    a = alfalfa(code=None, cz_max=20000.0)
    ra, de = np.radians(a["RA"]), np.radians(a["Dec"])
    def zahl(i):
        c = np.sin(de[i]) * np.sin(de) + np.cos(de[i]) * np.cos(de) * np.cos(ra - ra[i])
        return int(np.sum((a["D"][i] * np.arccos(np.clip(c, -1, 1)) < R) & (abs(a["cz"] - a["cz"][i]) < dv))) - 1
    i = int(np.where(a["AGC"] == agc)[0][0])
    wahl = np.where((abs(a["D"] - a["D"][i]) < 10) & (abs(a["lgM"] - a["lgM"][i]) < 0.5) & (a["AGC"] != agc))[0]
    wahl = np.random.default_rng(saat).choice(wahl, size=min(n, len(wahl)), replace=False)
    nc, k = np.array([zahl(j) for j in wahl]), zahl(i)
    return k, float(np.mean(nc < k) + 0.5 * np.mean(nc == k)), float(np.median(nc))


def alfalfa_breite_umgebung(D_von=60.0, D_bis=110.0, lg_von=8.8, lg_bis=9.4, n=1500, R=2.0, dv=500.0, saat=2):
    """α.100 (Code 1 und 2): Quellen mit Abstand D_von bis D_bis [Mpc] und lg M_HI zwischen lg_von und lg_bis, zufällig n davon,
    Rest lg W50 gegen eine Gerade in lg M_HI, getrennt nach der Zahl der Nachbarn (in R [Mpc] auf dem Himmel, |Δcz| < dv). Gibt
    (Zahl ohne Nachbarn, ihr mittlerer Rest, sein Fehler, mittlerer Rest der übrigen)."""
    a = alfalfa(code=None, cz_max=20000.0)
    ra, de = np.radians(a["RA"]), np.radians(a["Dec"])
    idx = np.where((a["D"] > D_von) & (a["D"] < D_bis) & (a["lgM"] > lg_von) & (a["lgM"] < lg_bis))[0]
    idx = np.random.default_rng(saat).choice(idx, size=min(n, len(idx)), replace=False)
    def zahl(i):
        c = np.sin(de[i]) * np.sin(de) + np.cos(de[i]) * np.cos(de) * np.cos(ra - ra[i])
        return int(np.sum((a["D"][i] * np.arccos(np.clip(c, -1, 1)) < R) & (abs(a["cz"] - a["cz"][i]) < dv))) - 1
    k, m = np.array([zahl(i) for i in idx]), a["lgM"][idx]
    A = np.vstack([np.ones_like(m), m]).T
    r = np.log10(a["W50"][idx]) - A @ np.linalg.lstsq(A, np.log10(a["W50"][idx]), rcond=None)[0]
    leer = k == 0
    return int(leer.sum()), float(r[leer].mean()), float(r[leer].std() / math.sqrt(leer.sum())), float(r[~leer].mean())


def alfalfa_flaeche(zelle=1.0):
    """Fläche des Feldes von α.100 in Quadratgrad: die Zellen von zelle Grad im Quadrat innerhalb der Grenzen des Artikels
    (07h30m bis 16h30m und 22h bis 03h, 0° bis 36°), in denen mindestens eine Quelle liegt, mit ihrer wahren Fläche."""
    a = alfalfa(code=1, cz_max=20000.0)
    ra, dec = a["RA"], a["Dec"]
    belegt = set(zip(np.floor(ra / zelle).astype(int), np.floor(dec / zelle).astype(int)))
    flaeche = 0.0
    for i, j in belegt:
        r0, d0 = i * zelle, j * zelle
        drin = (112.5 <= r0 < 247.5 or r0 >= 330.0 or r0 < 45.0) and 0.0 <= d0 < 36.0
        if drin:
            flaeche += zelle * np.degrees(np.sin(np.radians(d0 + zelle)) - np.sin(np.radians(d0)))
    return float(flaeche)


def alfalfa_hi_dichte(H0=70.0, grenze="50", feld=None, flaeche=None, breiten=None, massen=None):
    """Mittlere Dichte des neutralen Wasserstoffs aus ALFALFA α.100 mit 1/V_max: Code 1, cz bis 15 000 km/s, nur Quellen
    über der Grenze der Vollständigkeit (Haynes u. a. 2011, Gl. 4 und 5) für ihre Breite W50, V_max bis zum Abstand, an dem
    die Quelle an die Grenze fiele, höchstens bis zu dem von cz = 15 000 km/s. Abstände von H0 = 70 auf H0 umgerechnet.
    feld: "fruehling" oder "herbst" für eines der beiden Felder. Gibt (ρ_HI in M☉/Mpc³, Zahl der Quellen). breiten: Kanten
    in lg W50 [km/s], dann die Breitenfunktion, Zahl der Galaxien je Dekade der Breite [1/Mpc³] und Zahl der Quellen je Stück.
    massen: Kanten in lg M_HI [M☉], dann ρ_HI je Stück [M☉/Mpc³] und Zahl der Quellen je Stück."""
    a = alfalfa()
    lw = np.log10(a["W50"])
    s90 = np.where(lw < 2.5, 0.5 * lw - 1.14, lw - 2.39)
    sg = 10**(s90 - {"90": 0.0, "50": 0.067, "25": 0.102}[grenze])
    w = a["S"] > sg
    if feld == "fruehling":
        w &= (a["RA"] > 100) & (a["RA"] < 260)
    elif feld == "herbst":
        w &= (a["RA"] < 100) | (a["RA"] > 260)
    f = 70.0 / H0
    D = a["D"][w] * f
    M = 10**a["lgM"][w] * f**2
    dmax = np.minimum(D * np.sqrt(a["S"][w] / sg[w]), 15000.0 / H0)
    om = (alfalfa_flaeche() if flaeche is None else flaeche) * (np.pi / 180)**2
    if feld is not None:
        om *= {"fruehling": 1.0, "herbst": 1.0}[feld] * _alfalfa_anteil(feld)
    V = om / 3 * dmax**3
    if breiten is not None:
        lw, kanten = np.log10(a["W50"][w]), np.asarray(breiten)
        return np.histogram(lw, kanten, weights=1 / V)[0] / np.diff(kanten), np.histogram(lw, kanten)[0]
    if massen is not None:
        lm, kanten = np.log10(M), np.asarray(massen)
        return np.histogram(lm, kanten, weights=M / V)[0], np.histogram(lm, kanten)[0]
    return float(np.sum(M / V)), int(w.sum())


def _alfalfa_anteil(feld):
    """Anteil eines Feldes an der belegten Fläche (Zellen von 1°)."""
    a = alfalfa(code=1, cz_max=20000.0)
    zellen = set(zip(np.floor(a["RA"]).astype(int), np.floor(a["Dec"]).astype(int)))
    fl = lambda z: np.degrees(np.sin(np.radians(z[1] + 1)) - np.sin(np.radians(z[1])))
    drin = [z for z in zellen if (112 <= z[0] < 248 or z[0] >= 330 or z[0] < 45) and 0 <= z[1] < 36]
    teil = [z for z in drin if (100 < z[0] < 260) == (feld == "fruehling")]
    return sum(map(fl, teil)) / sum(map(fl, drin))


def gasanteil_heute(zelle=2.0, grenze="50"):
    """Anteil F der Baryonen, der heute als neutrales Gas (Wasserstoff mit Helium) in Galaxien gebunden ist, aus ALFALFA:
    ρ_HI mit 1/V_max bei H0 des Modells, geteilt durch X = 1 − Y_p und durch ρ_b = ω_b·ρ_krit,100. zelle: Seitenlänge der
    Zellen, aus denen die Fläche des Feldes gezählt wird. Gibt F in Prozent."""
    import kern as K
    h = K.H0_kms / 100
    rho, _ = alfalfa_hi_dichte(H0=K.H0_kms, grenze=grenze, flaeche=alfalfa_flaeche(zelle))
    rho_b = K.omega_b * K.rho_krit100 / K.M_sun * K.Mpc**3
    return float(100 * rho / (1 - K.YP_GEMESSEN[0]) / rho_b)


def hipass(spalten=None):
    """HICAT (Meyer u. a. 2004): (Geschwindigkeit im Rahmen der Lokalen Gruppe [km/s], im galaktischen Ruhesystem [km/s],
    Spitzenfluss [Jy], integrierter Fluss [Jy km/s], Breite W50 [km/s] als Mittel aus W50max und W50min). spalten: Namen
    der Spalten, dann ein dict mit ihnen (leere Felder NaN)."""
    zeilen = [z.rstrip("\n").split("\t") for z in open(os.path.join(ORT, "hipass", "hicat.tsv"), encoding="latin-1")
              if z[:1] not in ("#", "\n") and z.split("\t")[0].strip().isdigit()]
    kopf = [z.rstrip("\n").split("\t") for z in open(os.path.join(ORT, "hipass", "hicat.tsv"), encoding="latin-1")
            if z.startswith("recno")][0]
    namen = ("RVLG", "RVGSR", "Speak", "Sint", "W50max", "W50min") if spalten is None else tuple(spalten)
    i = {k: kopf.index(k) for k in namen}
    wert = lambda z, k: float(z[i[k]]) if z[i[k]].strip() else np.nan
    d = np.array([[wert(z, k) for k in namen] for z in zeilen])
    if spalten is not None:
        return {k: d[:, j] for j, k in enumerate(namen)}
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3], (d[:, 4] + d[:, 5]) / 2


def hipass_hi_dichte(H0=70.0, sp=0.084, si=9.4, breiten=None):
    """Mittlere Dichte des neutralen Wasserstoffs aus HICAT mit 1/V_max: nur Quellen über den Grenzen, an denen 99 % erfasst
    sind (Spitzenfluss sp in Jy, integrierter Fluss si in Jy km/s), Abstand v_LG/H0, M_HI = 2,356·10⁵·D²·S_int,
    V_max bis zum Abstand, an dem die Quelle an die erste der beiden Grenzen fiele, höchstens bis 12 700 km/s, ohne die
    Kugel unter 300 km/s, Fläche südlich von +2°. Gibt (ρ_HI in M☉/Mpc³, Zahl der Quellen). breiten: Kanten in lg W50
    [km/s], dann die Breitenfunktion, Zahl der Galaxien je Dekade der Breite [1/Mpc³] und Zahl der Quellen je Stück."""
    v, vg, p, s, W = hipass()
    w = (p >= sp) & (s >= si) & (vg >= 300) & (v > 300) & (v <= 12700)
    D = v[w] / H0
    M = 2.356e5 * D**2 * s[w]
    dmax = np.minimum(D * np.minimum(np.sqrt(s[w] / si), np.sqrt(p[w] / sp)), 12700.0 / H0)
    om = 2 * np.pi * (1 + np.sin(np.radians(2.0)))
    V = om / 3 * (dmax**3 - (300.0 / H0)**3)
    if breiten is not None:
        lw, kanten = np.log10(W[w]), np.asarray(breiten)
        return np.histogram(lw, kanten, weights=1 / V)[0] / np.diff(kanten), np.histogram(lw, kanten)[0]
    return float(np.sum(M / V)), int(w.sum())


def gasanteil_hipass(sp=0.084, si=9.4):
    """Anteil F der Baryonen als neutrales Gas mit Helium aus HIPASS, wie gasanteil_heute aus ALFALFA. Gibt F in Prozent."""
    import kern as K
    rho, _ = hipass_hi_dichte(H0=K.H0_kms, sp=sp, si=si)
    rho_b = K.omega_b * K.rho_krit100 / K.M_sun * K.Mpc**3
    return float(100 * rho / (1 - K.YP_GEMESSEN[0]) / rho_b)


@lru_cache(maxsize=None)
def gas_gezaehlt():
    """Das gezählte gebundene Gas: F₀, der Anteil der Baryonen, der heute als neutrales Gas mit Helium in Galaxien liegt,
    im Mittel aus ALFALFA (Zellen von 2°, Grenze 50 %) und HIPASS, dazu das Verhältnis der Dichte des neutralen Wasserstoffs
    bei z = 0,34 (Bera u. a. 2019) zu heute (Mittel beider Zählungen bei H₀ = 70) mit seinem Fehler (Statistik und Feld).
    Für das neutrale Gas die obere Grenze, gebundenes ionisiertes Gas ist nicht gezählt. Gibt (F₀, R, σ_R, z)."""
    z, rho, e1, e2 = hi_frueh()
    heute = (alfalfa_hi_dichte(H0=70.0, flaeche=alfalfa_flaeche(2.0))[0] + hipass_hi_dichte(H0=70.0)[0]) / 2
    return (gasanteil_heute() + gasanteil_hipass()) / 200, rho / heute, math.hypot(e1, e2) / heute, z


def _bera_rest(p):
    """Rest des Gases im Kreis bei z = 0,34 gegen das gezählte Verhältnis zu heute (Bera), p mit F0, k und a_s."""
    import kern as K
    _, R_m, R_e, z_b = gas_gezaehlt()
    return (K.gas_kreis(1 / (1 + z_b), p["F0"], p["k"], p["a_s"])[0] / p["F0"] - R_m) / R_e


def _frueh_rest(p):
    """Rest des Gases im Kreis bei z = 1,04 gegen das Verhältnis des neutralen Gases sternbildender Galaxien bei gleicher
    Verteilung der Sternmasse zu heute (Chowdhury, hi_frueh_eins), p mit F0, k und a_s."""
    import kern as K
    z, R, e, _ = hi_frueh_eins()
    return (K.gas_kreis(1 / (1 + z), p["F0"], p["k"], p["a_s"])[0] / p["F0"] - R) / e


@lru_cache(maxsize=None)
def abstaende_gas(form=None, F0="gezaehlt", bera=True, kopplung=False, frueh=False):
    """Die Abstände allein lesen den Kreis aus Rückstellung, Dehnung und Bindung: DESI DR2 mit fester Skala (r_d aus ω_b
    und ω_c = Q·ω_b, kern.schallhorizont), Pantheon+ in der Form, D/H und Y_p für ω_b, Q = 2π − 1. Die Dunkle Energie ist
    der Formteil mal (1 + j·F(a)) mit dem gebundenen Gas im Kreis (kern.gas_kreis, Maximum an der Sättigung des Formteils,
    k ≥ 0 frei), H₀ aus Formteil und Gas heute; bera: das Verhältnis bei z = 0,34 zu heute als Rest. form: a₀ des Formteils
    [10⁻¹⁰ m/s²], None frei; F0: Gas heute, "gezaehlt" aus gas_gezaehlt, None frei, 0 der Formteil allein. frueh: dazu
    die Lage a_s des Maximums frei und das Verhältnis bei z = 1,04 (hi_frueh_eins) als zweiter Rest. kopplung: der
    Zustandsteil wird aus der Mode gespeist, die Mode verliert, was er gewinnt, so wirkt er auf die Dehnung mit seiner
    über das Volumen gemittelten Geschichte, ρ_Z,eff(a) = 3a⁻³·∫₀^a ρ_Z(a₁)·a₁²·da₁ (kern.gas_gekoppelt).
    Gibt dict(chi2 mit den Resten des Gases, form, fehler_form, F0, fehler_F0, k, a_s, R, R_eins, wb)."""
    import kern as K
    from scipy.optimize import least_squares
    q = K.q_je_radiant()
    j = K.zustandsteil_je_gas(q) / K.B1_ZELLE
    wnu = 0.06 / 93.14
    F_z, R_m, R_e, z_b = gas_gezaehlt()
    if F0 == "gezaehlt":
        F0 = F_z
    frei = ["wb"] + (["form"] if form is None else []) + (["F0"] if F0 is None else []) + (["k"] if F0 != 0 else []) + (
        ["a_s"] if frueh and F0 != 0 else [])
    st = dict(wb=0.02237, form=1.06, F0=0.0097, k=3.0, a_s=0.45)
    sch = dict(wb=1e-4, form=0.01, F0=0.003, k=2.0, a_s=0.03)
    z_e, R_1, e_1, _ = hi_frueh_eins()

    def lage(v):
        return v["a_s"] if frueh else (((1 + q) * v["wb"] + wnu) * K.rho_krit100
                                        / K.dichte_dunkle_energie(v["form"] * 1e-10))**(1 / 3)

    def werte(y):
        v = dict(st, **{n: st[n] + sch[n] * x for n, x in zip(frei, y)})
        v["form"] = v["form"] if form is None else form
        v["F0"] = v["F0"] if F0 is None else F0
        if F0 == 0:
            v["k"] = 0.0
        return v

    def rest(y):
        v = werte(y)
        wb, f, F, k = v["wb"], v["form"] * 1e-10, v["F0"], v["k"]
        a_s = lage(v)
        if kopplung:
            de, F_heute = K.gas_gekoppelt(F, k, a_s)
        else:
            de, F_heute = (lambda zz: K.gas_kreis(1 / (1 + np.asarray(zz)), F, k, a_s)[1]), F
        H0 = K.atmen_heute(f * math.sqrt(1 + j * F_heute), Q=q, wb=wb, wnu=wnu)[0]
        h2 = (H0 / 100)**2
        om, or_ = ((1 + q) * wb + wnu) / h2, K.dichte_heute()[0] / K.rho_krit100 / h2
        t = [bao_absolut(H0, om, or_, K.schallhorizont(q, wb), rest=True, de=de), pantheon_form(om, or_, rest=True, de=de),
             bbn_rest(wb)]
        if bera and F0 != 0:
            t.append([(K.gas_kreis(1 / (1 + z_b), F, k, a_s)[0] / F - R_m) / R_e])
        if frueh and F0 != 0:
            t.append([(K.gas_kreis(1 / (1 + z_e), F, k, a_s)[0] / F - R_1) / e_1])
        return np.concatenate(t)

    lo = [0.0 if n == "k" else 0.05 if n == "a_s" else -np.inf for n in frei]
    erg = least_squares(rest, np.zeros(len(frei)), bounds=([(l - st[n]) / sch[n] for l, n in zip(lo, frei)],
                                                           [np.inf] * len(frei)), diff_step=1e-3)
    kov = np.linalg.pinv(erg.jac.T @ erg.jac)
    v = werte(erg.x)
    feh = {n: float(sch[n] * math.sqrt(kov[i, i])) for i, n in enumerate(frei)}
    a_s = lage(v)
    R = float(K.gas_kreis(1 / (1 + z_b), v["F0"], v["k"], a_s)[0] / v["F0"]) if v["F0"] else None
    R1 = float(K.gas_kreis(1 / (1 + z_e), v["F0"], v["k"], a_s)[0] / v["F0"]) if v["F0"] else None
    return dict(chi2=float(erg.fun @ erg.fun), form=float(v["form"]), fehler_form=feh.get("form", 0.0), F0=float(v["F0"]),
                fehler_F0=feh.get("F0", 0.0), k=float(v["k"]), a_s=float(a_s), R=R, R_eins=R1, wb=float(v["wb"]))
