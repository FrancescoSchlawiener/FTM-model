# -*- coding: utf-8 -*-
"""Unsere Daten, nur gelesen: SPARC (Lelli, McGaugh, Schombert 2016) und AME2020 (Huang u. a. 2021)."""
import os

import numpy as np

ORT = os.path.join(os.path.dirname(__file__), "..", "archiv", "weltmodell", "daten")


def sparc_galaxien():
    """Tafel der 175 Galaxien: Name → dict(T, D, Inc, L36, MHI, Vflat, Q)."""
    aus = {}
    zeilen = open(os.path.join(ORT, "sparc", "SPARC_Lelli2016c.mrt"), encoding="latin-1").readlines()
    for z in zeilen:
        p = z.split()
        if len(p) < 18:
            continue
        try:
            aus[p[0]] = dict(T=int(p[1]), D=float(p[2]), fD=int(p[4]), Inc=float(p[5]), L36=float(p[7]), MHI=float(p[13]),
                             Vflat=float(p[15]), Q=int(p[17]))
        except ValueError:
            continue
    return aus


def sparc_kurve(name):
    """Rotationskurve: Felder R [kpc], Vobs, eV, Vgas, Vdisk, Vbul [km/s]."""
    d = np.loadtxt(os.path.join(ORT, "sparc", "Rotmod_LTG", name + "_rotmod.dat"), comments="#")
    return dict(R=d[:, 0], V=d[:, 1], eV=d[:, 2], Vgas=d[:, 3], Vdisk=d[:, 4], Vbul=d[:, 5])


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
def sparc_auswahl(yd=YD, H0=None):
    """Galaxien mit Q ≤ 2 und Neigung ≥ 30°, Punkte mit V > 0, g_bar > 0, Fehler ≤ 10 %; mit Fehlern von Abstand und Neigung.
    H0: Abstände aus dem Hubble-Fluss mit diesem H₀ statt 73 (R ∝ D, V_bar ∝ √D, V bleibt)."""
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
            aus.append((n, p, {x: k[x][m] for x in k}, eD[n] * f, eI[n]))
    return tuple(aus)


def _reste(v, a0, nu, p, k, eD, eI, kpc, yd, halo=None, roh=False):
    """Reste der Regel in log g mit Abstand, Neigung und Υ_d als Störgrößen (Vorgaben als Gauß-Glieder).
    halo = H₀ [km/s/Mpc]: statt der Regel das Standardmodell, Faltung der Baryonen plus NFW-Halo mit v[3] = lg V₂₀₀ [km/s],
    v[4] = lg c (R₂₀₀ = V₂₀₀/(10·H₀)). roh: Reste je Punkt in dex und lg(g_bar/a₀), ohne Störglieder."""
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
        gm = nu(gb / a0) * gb
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
def sparc_profil(form, kpc, yd=YD, H0=None):
    """χ² je Galaxie über dem Gitter von a₀, Störgrößen je Punkt des Gitters ausgeglichen (warmer Start)."""
    from scipy.optimize import least_squares
    nu = FORMEN[form]
    C = np.zeros((len(sparc_auswahl(yd, H0)), len(GITTER)))
    for gi, (n, p, k, eD, eI) in enumerate(sparc_auswahl(yd, H0)):
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


def sparc_reste(kpc, a0, yd=YD):
    """Reste der Regel je Punkt in dex bei der Schwelle a₀, Störgrößen je Galaxie ausgeglichen; dazu lg(g_bar/a₀)."""
    from scipy.optimize import least_squares
    nu, x, r = FORMEN["Bose"], [], []
    for n, p, k, eD, eI in sparc_auswahl(yd):
        inc = math.radians(p["Inc"])
        lo = np.array([-0.4, max(-0.6, 0.25 - inc), -0.6])
        hi = np.array([0.4, max(1e-6, math.pi / 2 - inc) + 1e-9, 0.6])
        v = least_squares(_reste, np.zeros(3), bounds=(lo, hi), args=(a0, nu, p, k, eD, eI, kpc, yd)).x
        d, g = _reste(v, a0, nu, p, k, eD, eI, kpc, yd, roh=True)
        r += list(d); x += list(g)
    return np.array(x), np.array(r)


def _bestes(v):
    from scipy.interpolate import CubicSpline
    from scipy.optimize import minimize_scalar
    s = CubicSpline(GITTER, v)
    r = minimize_scalar(s, bounds=(GITTER[0], GITTER[-1]), method="bounded")
    return 10**r.x, float(s(r.x))


@lru_cache(maxsize=None)
def sparc_schwelle(kpc, form="Bose", ziehungen=1000, yd=YD, H0=None):
    """a₀ aus SPARC mit der Regel: bester Wert und 16–84 % über Ziehen der Galaxien."""
    C = sparc_profil(form, kpc, yd, H0)
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
def pantheon_h0(om, or_, zlo, zhi):
    """H₀ aus Pantheon+ mit den Cepheiden-Eichungen (volle Kovarianz, Form des Atmens mit Ω_m, Ω_r):
    (H₀ [km/s/Mpc], Fehler, Zahl der Supernovae im Fluss, Helligkeit M)."""
    import gzip
    p = os.path.join(ORT, "pantheon")
    t = np.genfromtxt(os.path.join(p, "Pantheon+SH0ES.dat"), names=True, dtype=None, encoding=None)
    with gzip.open(os.path.join(p, "Pantheon+SH0ES_STAT+SYS.cov.gz"), "rt") as f:
        n = int(f.readline())
        C = np.loadtxt(f).reshape(n, n)
    E = lambda z: np.sqrt(or_ * (1 + z)**4 + om * (1 + z)**3 + 1 - om - or_)
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


@lru_cache(maxsize=None)
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


def zwerg_streuung(z, yv, mit_feld=True, geschichte=False, tief=False):
    """σ [km/s] des Zwergs z nach der Regel auf das Gesamtfeld (Plummer, Masse Υ_V·L_V + 1,33·M_HI), im Newtonfeld des Wirts
    als Punktmasse; mit_feld=False: allein; geschichte: das Gas nach zwerg_gas liest nur die Rückkehr (Gedächtnis der
    Bindung), die Neigung am Ort liest das heutige; tief: allein ganz im Tiefbereich, σ⁴ = (4/81)·GMa₀ mit der Masse der
    Rückkehr, unabhängig von Profil und Radius."""
    import kern as K
    M = (yv * z["L"] + 1.33 * z["M_HI"]) * K.M_sun
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


def kids_typ_reste(typ, lg_von=-20.0, lg_bis=0.0, g_e=0.0):
    """Rest lg(g_obs/Regel(g_bar)) der Linsen eines Typs in den Bins mit lg g_bar zwischen lg_von und lg_bis, im gleichförmigen
    äußeren Feld g_e [m/s²]: (Zahl, gewichtetes Mittel, Fehler, χ² gegen null)."""
    import kern as K
    gb, go, ge = kids_typ(typ)
    w = (np.log10(gb) >= lg_von) & (np.log10(gb) <= lg_bis)
    r = np.log10(go[w] / np.array([K.regel_monopol(x, g_e) for x in gb[w]]))
    e = ge[w] / go[w] / math.log(10)
    g = 1 / e**2
    return int(w.sum()), float(np.sum(g * r) / g.sum()), float(1 / math.sqrt(g.sum())), float(np.sum((r / e)**2))


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


def udg_reste():
    """Gasreiche ultradiffuse Galaxien (Mancera Piña u. a. 2022, 6, äußerer Punkt): lg(v_obs/v_Regel) und lg(v_obs/v_Newton) je
    Galaxie, die Regel auf das Newtonsche Feld der Baryonen am selben Radius. Gibt (Zahl, Mittel gegen die Regel, Mittel gegen
    Newton, kleinstes und größtes a_N/a₀)."""
    import kern as K
    d = np.loadtxt(os.path.join(ORT, "linsen_haufen", "udg_mancera_pina2022.txt"))
    R = d[:, 2] * K.kpc
    gN = (d[:, 6] * 1e3)**2 / R
    vr = np.sqrt(np.array([K.verstaerkung(x) * x for x in gN]) * R) / 1e3
    return len(d), float(np.mean(np.log10(d[:, 3] / vr))), float(np.mean(np.log10(d[:, 3] / d[:, 6]))), \
        float(gN.min() / K.a0), float(gN.max() / K.a0)


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


def zwerge_reste(gruppe, mit_feld, yv=2.0, streuung=0.1, geschichte=False, tief=False):
    """Rest lg(σ_gemessen/σ_Regel) der Zwerge einer Gruppe: "isoliert", "innen" (Satelliten, eigenes Feld am Halblichtradius
    überwiegt das des Wirts), "außen" (klassische Satelliten über 10⁵ L☉, Wirtsfeld überwiegt), "ultraschwach" (Satelliten
    bis 10⁵ L☉), auch "gasreich", "gasarm", "klassisch"; geschichte und tief wie in zwerg_streuung. Gibt (Zahl, Mittel,
    Fehler des Mittels, χ² mit eigener Streuung)."""
    import kern as K
    G = [z for z in zwerge_gruppe(gruppe, yv) if geschichte != "metall" or z["feh"] is not None]
    r = np.array([math.log10(z["sigma"] / zwerg_streuung(z, yv, mit_feld, geschichte, tief)) for z in G])
    e = np.array([z["fehler"] / z["sigma"] / math.log(10) for z in G])
    return len(G), float(r.mean()), float(r.std() / math.sqrt(len(G))), float(np.sum(r**2 / (e**2 + streuung**2)))


@lru_cache(maxsize=None)
def kugelhaufen():
    """Sternhaufen mit gemessener Sichtlinienstreuung (Local Volume Database, ohne obere Grenzen): je Haufen dict mit Name,
    Wirt, Abstand zum Wirt [kpc], L_V [L☉], R_e [m] (kreisgemittelt), σ und Fehler [km/s], lg(g_N/a₀) am Halblichtradius
    der Kugel (4/3·R_e) für Υ_V = 2."""
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
def kugelhaufen_reste(newton=False, auswahl="alle", streuung=0.1):
    """Dünne Sternhaufen gegen dichte: Υ_V geeicht an den dichten (g_N > 3·a₀, dort ν < 1,22), dann Rest lg(σ_gemessen/σ)
    der dünnen (g_N < a₀) mit der Regel oder Newton (kern.streuung_im_feld, ohne äußeres Feld). auswahl: "alle", "gut"
    (Fehler der Streuung unter 30 %), "halo" (Milchstraße, über 10 kpc vom Zentrum). Gibt (Υ_V, Zahl dicht, Zahl dünn,
    gewichtetes Mittel, Fehler, χ² mit eigener Streuung)."""
    import kern as K
    a = 1e-40 if newton else K.a0
    H = kugelhaufen()
    rest = lambda z, yv: math.log10(z["sigma"] / (K.streuung_im_feld(yv * z["L"] * K.M_sun, z["R_e"], 0.0, a) / 1e3))
    gewicht = lambda z: 1 / ((z["fehler"] / z["sigma"] / math.log(10))**2 + streuung**2)
    dicht = [z for z in H if z["lg_g"] > math.log10(3)]
    yv = 10**(2 * sum(gewicht(z) * rest(z, 1.0) for z in dicht) / sum(gewicht(z) for z in dicht))
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


def kids_scherung_chi2(punkt, verhaeltnis=False, **auswahl):
    """χ² der KiDS-1000-Scherung (kids_scherung) für den Punkt (ω_b, ω_c, H₀, ln10¹⁰A_s, n_s, w₀, w_a): Faltung der Mode wie
    druckfreier Staub, Ausrichtung der Bilder mit freier Amplitude A; auswahl: winkel, art wie bei kids_scherung. Gibt (χ², A,
    S8, Zahl der Werte); verhaeltnis: dazu Messung gegen Faltung bei diesem A, linear geschätzt mit der vollen Kovarianz,
    R = gᵀC⁻¹(d − m_A)/gᵀC⁻¹g mit dem Fehler 1/√(gᵀC⁻¹g)."""
    import kern as K
    from scipy.optimize import minimize_scalar
    zeilen, C, nz = kids_scherung(**auswahl)
    ell, cl, S8 = K.scherung_cl(*punkt[:5], nz, *punkt[5:])
    th = sorted({z[3] for z in zeilen})
    xi = {}
    for (i, j), (gg, gi, ii) in cl.items():
        xi[(i, j)] = [K.scherung_xi(ell, c, th) for c in (gg, gi, ii)]
    d = np.array([z[4] for z in zeilen])
    Ci = np.linalg.inv(C)
    def modell(A):
        m = []
        for vz, i, j, t, _ in zeilen:
            k, s = th.index(t), 0 if vz == "+" else 1
            g, gi, ii = xi[(min(i, j), max(i, j))]
            m.append(g[s][k] + A * gi[s][k] + A * A * ii[s][k])
        return np.array(m)
    chi = lambda A: float((d - modell(A)) @ Ci @ (d - modell(A)))
    r = minimize_scalar(chi, bounds=(-3, 6), method="bounded")
    if verhaeltnis:
        g = np.array([xi[(min(i, j), max(i, j))][0][0 if vz == "+" else 1][th.index(t)] for vz, i, j, t, _ in zeilen])
        F = float(g @ Ci @ g)
        return r.fun, float(r.x), float(S8), len(d), float(g @ Ci @ (d - modell(r.x) + g)) / F, 1 / math.sqrt(F)
    return r.fun, float(r.x), float(S8), len(d)


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


def xcop_staub(c=None, von=0.05):
    """X-COP von von·R₅₀₀ bis R₂₀₀: die hydrostatische Masse gegen Baryonen plus Mode, beide an R₂₀₀ mit dem Maß des Ganzen
    geeicht, M(R₂₀₀) = (1 + Q)·M_b(R₂₀₀). Die Baryonen liegen wie das Gas; die Mode mit dem Gewicht Q·M_b(R₂₀₀) entweder
    ebenso wie das Gas (c = None) oder wie druckfreier Staub nach dem Einfall, NFW-Form mit der Konzentration c = R₂₀₀/r_s.
    Gibt (χ² in lg M mit dem halben Band als Fehler, Zahl der Punkte, mittlerer Rest innerhalb 0,15·R₅₀₀)."""
    import kern as K
    Q = K.q_je_radiant()
    rc = K.rho_krit100 * (K.H0_kms / 100)**2
    nfw = lambda x: math.log(1 + x) - x / (1 + x)
    chi, n, kern_rest = 0.0, 0, []
    for v in xcop().values():
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


def haufen_q(ohne=()):
    """Q und Konzentration c, die HIFLUGCS mit der Mode in der Form nach dem Einfall liest (beide frei, hiflugcs_staub);
    ohne: ausgelassene Haufen. Gibt (Q, c, χ²)."""
    from scipy.optimize import minimize
    r = minimize(lambda p: hiflugcs_staub(p[0], Q=p[1], ohne=ohne)[0], [5.0, 5.5], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-4))
    return float(r.x[1]), float(r.x[0]), float(r.fun)


def haufen_q_fehler():
    """Fehler von Q aus HIFLUGCS durch Auslassen je eines Haufens (Jackknife über die 10 Haufen)."""
    namen = sorted({x[0] for x in hiflugcs()})
    q = np.array([haufen_q((n,))[0] for n in namen])
    return float(math.sqrt((len(q) - 1) / len(q) * np.sum((q - q.mean())**2)))


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


def kids_heiss(typ, lg_von=-20.0, lg_bis=0.0):
    """Linsen eines Typs (kids_typ), der Überschuss über die Regel auf Sterne und kaltes Gas gelesen wie im Haufen, als nie
    abgekühltes Gas mit dem Maß des Ganzen: M_heiß/M_b = (g_obs − ν(g_bar)·g_bar)/((1 + Q)·g_bar) je Bin zwischen lg_von und
    lg_bis. Gibt (Zahl, Median, kleinster, größter Wert)."""
    import kern as K
    Q = K.q_je_radiant()
    gb, go, _ = kids_typ(typ)
    w = (np.log10(gb) >= lg_von) & (np.log10(gb) <= lg_bis)
    f = (go[w] - np.array([K.verstaerkung(x) * x for x in gb[w]])) / ((1 + Q) * gb[w])
    return int(w.sum()), float(np.median(f)), float(f.min()), float(f.max())


def xcop_abkuehlung(alter_gyr=13.8, X=0.75, ort=1.0):
    """X-COP: Temperatur des Gases aus dem hydrostatischen Gleichgewicht (Druck von außen nach innen integriert, null am letzten
    Punkt), Abkühlzeit durch Bremsstrahlung t = (3/2)(n_e + n_i)kT/(n_e·n_i·Λ) bei Wasserstoffanteil X. Je Haufen der Massenanteil
    des Gases innerhalb R₅₀₀, das sich in alter_gyr [Gyr] hätte abkühlen können, Temperatur [keV] und Abkühlzeit [Gyr] bei
    ort·R₅₀₀. Gibt (Median und größter Anteil, kleinste und größte Temperatur, kürzeste und längste Abkühlzeit)."""
    import kern as K
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
    des hydrostatischen Felds zur Regel auf das Gas allein; dazu dasselbe Verhältnis bei 0,1·R₅₀₀."""
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
                      verlangt_200=verlangt(i2), verlangt_kern=verlangt(ik))
    return aus


PLANCK_START = dict(wb=0.02237, wc=0.1200, H0=67.36, a0=1.06, lnAs=3.044, ns=0.9649, tau=0.0544, A=1.0, w0=-1.0, wa=0.0,
                    Q=2 * math.pi - 1)
PLANCK_SCHRITT = dict(wb=0.00015, wc=0.0012, H0=0.5, a0=0.01, lnAs=0.014, ns=0.004, tau=0.007, A=0.0025, w0=0.05, wa=0.2,
                      Q=0.05)


def _planck_punkt(p, Q):
    """Die sechs Zahlen für CAMB aus dem Punkt p: ΛCDM (Q = None) mit ω_c und H₀, sonst ω_c = Q·ω_b, H₀ aus a₀;
    Q = "lesen": Q aus dem Punkt."""
    import kern as K
    if Q is None:
        return p["wb"], p["wc"], p["H0"]
    q = p["Q"] if Q == "lesen" else Q
    return p["wb"], q * p["wb"], K.atmen_heute(p["a0"] * 1e-10, Q=q, wb=p["wb"], wnu=0.06 / 93.14)[0]


_CAMB = {}


def planck_reste(p, Q=None):
    """Geweißte Reste des Punkts p je Datensatz: Spektren (plik_lite mit TT unter ℓ = 30), Linse (PR4), DESI DR2,
    Pantheon+ (Nullpunkt frei), Wachstum f·σ₈, Urknall-Nukleosynthese und die Vorgaben für τ und die Eichung; Abstände
    mit Ω_m aus ω_b, ω_c und der Neutrinomasse, r_d aus CAMB; dazu σ₈ (kein Rest). Die Dunkle Energie heute setzt a₀
    (oder H₀), w₀ und w_a geben ihren Gang."""
    import kern as K
    wb, wc, H = _planck_punkt(p, Q)
    w0, wa = p.get("w0", -1.0), p.get("wa", 0.0)
    kenn = (wb, wc, H, p["lnAs"], p["ns"], p["tau"], w0, wa)
    if kenn not in _CAMB:
        _CAMB.clear()
        _CAMB[kenn] = K.cmb_spektren(*kenn[:6], alles=True, wachstum=wachstum()[0], w0=w0, wa=wa)
    tt, te, ee, pp, rd, fs8, s8 = _CAMB[kenn]
    _, v8, L8 = wachstum()
    h2 = (H / 100)**2
    om, or_ = (wb + wc + 0.06 / 93.14) / h2, K.dichte_heute()[0] / K.rho_krit100 / h2
    return dict(spektren=planck_rest(tt, te, ee, p["A"]), linse=planck_linse_rest(pp),
                bao=bao_absolut(H, om, or_, rd, w0, wa, rest=True), sn=pantheon_form(om, or_, w0, wa, rest=True),
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
    setzt dann die ganze Dunkle Energie heute, Formteil und Zustandsteil."""
    import huelle_ablage as HA
    from scipy.optimize import least_squares
    name = f"Q={'frei' if Q is None else Q if isinstance(Q, str) else round(Q, 6)} mit={'+'.join(mit)} fest={','.join(fest)}" + (
        " saettigung" if zustand == "saettigung" else " zustand" if zustand else "") + (
        f" start={sorted(start.items())}" if start else "") + (f" form={form}" if form else "")
    alle = HA.lesen()
    if not neu and name in alle.get("_cmb", {}):
        return alle["_cmb"][name]
    frei = [k for k in (("wb", "wc", "H0") if Q is None else ("wb", "a0") + (("Q",) if Q == "lesen" else ()))
            + ("lnAs", "ns", "tau", "A")
            + (("wa",) if zustand == "saettigung" else ("w0", "wa") if zustand else ()) if k not in fest]
    import kern as K

    def punkt(y):
        p = dict(PLANCK_START, **{k: PLANCK_START[k] + PLANCK_SCHRITT[k] * v for k, v in zip(frei, y)})
        if zustand == "saettigung":                                  # drei Teile holen die Faltung ein, im Punkt selbst
            q = p["Q"] if Q == "lesen" else Q
            a_s = (((1 + q) * p["wb"] + 0.06 / 93.14) * K.rho_krit100
                   / K.dichte_dunkle_energie((form or p["a0"]) * 1e-10))**(1 / 3)
            p["w0"] = -1 - p["wa"] * (1 - a_s)
        return p

    def rest(y):
        r = planck_reste(punkt(y), Q)
        return np.concatenate([r[k] for k in mit + ("vor",)])

    y0 = np.array([((start or {}).get(k, PLANCK_START[k]) - PLANCK_START[k]) / PLANCK_SCHRITT[k] for k in frei])
    erg = least_squares(rest, y0, diff_step=0.05, x_scale=1.0, xtol=1e-6, ftol=1e-8)
    s = np.array([PLANCK_SCHRITT[k] for k in frei])
    kov = np.linalg.inv(erg.jac.T @ erg.jac) * np.outer(s, s)
    p = punkt(erg.x)
    p["wb"], p["wc"], p["H0"] = _planck_punkt(p, Q)
    r = planck_reste(p, Q)
    sigma8 = float(r.pop("sigma8"))
    chi2 = {k: float(v @ v) for k, v in r.items()}
    chi2["gesamt"] = sum(chi2[k] for k in mit + ("vor",))
    aus = dict(chi2=chi2, sigma8=sigma8, punkt=p, fehler={k: float(math.sqrt(kov[i, i])) for i, k in enumerate(frei)},
               kovarianz=dict(namen=frei, werte=kov.tolist()))
    alle = HA.lesen()
    alle.setdefault("_cmb", {})[name] = aus
    HA.schreiben(alle)
    return aus


def bao_absolut(H0, om, or_, rd, w0=-1.0, wa=0.0, rest=False):
    """χ² von DESI DR2 mit fester Skala: H₀ [km/s/Mpc], Ω_m, Ω_r, flach, Schallhorizont r_d [Mpc] fest;
    Dunkle Energie mit w(a) = w₀ + w_a(1 − a). rest: der geweißte Rest statt χ²."""
    z, v, g, C = desi()
    x = np.linspace(0, 2.4, 24001)
    E = np.sqrt(or_ * (1 + x)**4 + om * (1 + x)**3
                + (1 - om - or_) * (1 + x)**(3 * (1 + w0 + wa)) * np.exp(-3 * wa * x / (1 + x)))
    I = np.concatenate([[0.0], np.cumsum((1 / E[1:] + 1 / E[:-1]) / 2 * np.diff(x))])
    dh0 = 299792.458 / H0
    dm, dh = dh0 * np.interp(z, x, I), dh0 / np.interp(z, x, E)
    th = np.array([{"DM_over_rs": a, "DH_over_rs": b, "DV_over_rs": (zz * a * a * b)**(1 / 3)}[gg]
                   for zz, a, b, gg in zip(z, dm, dh, g)]) / rd
    r = v - th
    if rest:
        import scipy.linalg
        return scipy.linalg.solve_triangular(np.linalg.cholesky(C), r, lower=True)
    return float(r @ np.linalg.inv(C) @ r)


def pantheon_form(om, or_, w0=-1.0, wa=0.0, rest=False):
    """χ² der Pantheon+-Form (Nullpunkt frei, flach) bei Ω_m, Ω_r; Dunkle Energie mit w(a) = w₀ + w_a(1 − a).
    rest: der geweißte Rest statt χ² (Nullpunkt herausprojiziert)."""
    zs, zh, m, C = pantheon()
    Ci = _pantheon_inv()
    g = np.linspace(0, 2.4, 24001)
    f = 1 / np.sqrt(or_ * (1 + g)**4 + om * (1 + g)**3
                    + (1 - om - or_) * (1 + g)**(3 * (1 + w0 + wa)) * np.exp(-3 * wa * g / (1 + g)))
    I = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(g))])
    r = m - 5 * np.log10((1 + zh) * np.interp(zs, g, I))
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


