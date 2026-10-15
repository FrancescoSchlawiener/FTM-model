"""Simulation der Kerne und Elemente.

Eine Rechnung für alles, was Kerne und Elemente betrifft: Hülle (Füllung nach n + l, Kopplung Zα, Grenze nach Dirac),
Kern (Tröpfchenformel und Schalenglied, beide einmal an den gemessenen Kernen von AME2020 angepasst; Tal, Spaltbarkeit x,
Kippstelle, Halterlänge, Q_α, S_n), Lebensdauer (Tunnelform an NUBASE2020) und die möglichen neuen Elemente.
Alle anderen Teile (Text, Lauf, Tests) rufen nur diese Funktionen auf.

Energien in MeV, als Faltungsenergie E = 3c⁴·Δs/G.
"""
import csv
import math
import os
import re
from functools import lru_cache

import numpy as np

import kern as K
import daten as D

MEV = K.e * 1e6
E_N = K.energie(K.m_n) / MEV                 # Faltungsenergie des Neutrons [MeV]
E_H = K.energie(K.m_H) / MEV                 # des Wasserstoffatoms
LAM_E = K.hbar / (K.m_e * K.c) * 1e15        # Quantenlänge des Elektrons [fm]

TROPFEN = ("tropfen",)
SCHALE = ("schale", 126)                     # Schalenglied mit den Schließzahlen des Modells


# ---------------------------------------------------------------- Schließzahlen
def schliesszahlen_kern(k_max=8):
    """Schließzahlen des gehaltenen Trichters: O(k) = (k+1)(k+2)(k+3)/3 bis 20, dann O(k) + 2(k + 2)."""
    O = lambda k: (k + 1) * (k + 2) * (k + 3) // 3
    return tuple(O(k) for k in range(3)) + tuple(O(k) + 2 * (k + 2) for k in range(2, k_max + 1))


def niveaus_verwindung(kappa, n_max=7, form="fest", mu=0.0):
    """Gehaltener Trichter: Schale N trägt l = N, N − 2, …; die Verwindung koppelt an die Bahnwindung, E = N − κ·⟨l·s⟩ in
    Einheiten der Schale, ⟨l·s⟩ = l/2 für j = l + ½ und −(l + 1)/2 für j = l − ½; form „je_schale“: Beiwert κ/(N + 1).
    mu: Abflachung des gehaltenen Trichters, die große Windungen senkt, −μ·(l(l + 1) − N(N + 3)/2) (Mittel der Schale null).
    Gibt (E, N, l, 2j), nach E geordnet."""
    aus = []
    for N in range(n_max + 1):
        for l in range(N, -1, -2):
            for j2 in ((2 * l + 1,) if l == 0 else (2 * l + 1, 2 * l - 1)):
                ls = l / 2 if j2 == 2 * l + 1 else -(l + 1) / 2
                k = kappa / (N + 1) if form == "je_schale" else kappa
                aus.append((N - k * ls - mu * (l * (l + 1) - N * (N + 3) / 2), N, l, j2))
    return sorted(aus)


def schliesszahlen_verwindung(kappa, luecke=0.5, form="fest", mu=0.0):
    """Schließzahlen aus der Füllung der Niveaus: Besetzungen, vor denen die Lücke zum nächsten Niveau mindestens
    luecke Schalen beträgt."""
    n, aus, niv = 0, [], niveaus_verwindung(kappa, form=form, mu=mu)
    for k, (E, N, l, j2) in enumerate(niv):
        n += j2 + 1
        if k + 1 < len(niv) and niv[k + 1][0] - E >= luecke:
            aus.append(n)
    return tuple(aus)


def grundspin(Z, A, kappa, form="fest", mu=0.0, niv=None, paarfenster=0.0):
    """Grundspin (2j) und Parität eines ungeraden Kerns aus der Füllung der Niveaus des gehaltenen Trichters mit
    Verwindung: Paare koppeln zu null, das ungepaarte Nukleon trägt j und (−1)^l; sonst None.
    paarfenster > 0: liegt innerhalb dieses Fensters (Schalen) ein leeres oder gerade besetztes Niveau kleinerer Windung,
    geht das Paar in die größere Windung und das einzelne Nukleon in die kleinere."""
    niv = niveaus_verwindung(kappa, form=form, mu=mu) if niv is None else niv
    def ungepaart(n):
        if n % 2 == 0:
            return None
        bes = []
        for E, N, l, j2 in niv:
            bes.append(min(n, j2 + 1))
            n -= bes[-1]
            if n == 0:
                break
        i = len(bes) - 1
        E0, j2 = niv[i][0], niv[i][3]
        kand = [k for k in range(len(niv)) if abs(niv[k][0] - E0) <= paarfenster and niv[k][3] < j2
                and (k >= len(bes) or bes[k] % 2 == 0 and bes[k] <= niv[k][3])]
        if paarfenster > 0 and kand:
            i = min(kand, key=lambda k: niv[k][3])
        return (niv[i][3], (-1)**niv[i][2])
    if Z % 2 == 1 and (A - Z) % 2 == 0:
        return ungepaart(Z)
    if (A - Z) % 2 == 1 and Z % 2 == 0:
        return ungepaart(A - Z)
    return None


def grundspin_pruefung(kappa, form="fest", mu=0.0, kugel=False, paarfenster=0.0):
    """Ungerade Kerne mit eindeutigem Spin (NUBASE2020) gegen die Füllung: (Treffer Spin und Parität, Treffer Parität, Zahl).
    kugel: nur außerhalb der verformten Gebiete 150 ≤ A ≤ 190 und A ≥ 220."""
    beide = par = n = 0
    niv = niveaus_verwindung(kappa, form=form, mu=mu)
    for (Z, A), (j2, p) in D.nubase_spins().items():
        if kugel and (150 <= A <= 190 or A >= 220):
            continue
        v = grundspin(Z, A, kappa, form, mu, niv, paarfenster)
        if v is None:
            continue
        n += 1
        par += v[1] == p
        beide += v == (j2, p)
    return beide, par, n


def paare_grundspin():
    """Grundspins nach der Zahl halber Drehungen (NUBASE2020, eindeutig zugewiesen): gerade-gerade (Zahl, davon 0⁺),
    ungerades A (Zahl, davon halbzahlig), ungerade-ungerade (Zahl, davon ganzzahlig)."""
    s = D.nubase_spins()
    gg = [v for (z, a), v in s.items() if z % 2 == 0 and (a - z) % 2 == 0]
    ua = [v for (z, a), v in s.items() if a % 2 == 1]
    uu = [v for (z, a), v in s.items() if z % 2 == 1 and (a - z) % 2 == 1]
    return (len(gg), sum(v == (0, 1) for v in gg), len(ua), sum(v[0] % 2 == 1 for v in ua), len(uu),
            sum(v[0] % 2 == 0 for v in uu))


def stabile_nach_paaren():
    """Stabile Grundzustände (NUBASE2020): Zahl gesamt, gerade-gerade, ungerades A, ungerade-ungerade, und ob alle
    ungerade-ungeraden N = Z haben."""
    s = D.nubase(stabil=True)
    uu = [(z, a) for z, a in s if z % 2 == 1 and (a - z) % 2 == 1]
    return (len(s), sum(1 for z, a in s if z % 2 == 0 and (a - z) % 2 == 0), sum(1 for z, a in s if a % 2 == 1), len(uu),
            all(a == 2 * z for z, a in uu))


def magie(f):
    """Schließzahlen (Neutronen, Protonen) einer Formel; bei SCHALE mit wählbarer nächster Protonen-Schließung."""
    m = (0,) + schliesszahlen_kern()
    if f[0] != "schale" or f[1] == 126:
        return m, m
    return m, tuple(f[1] if z == 126 else z for z in m)


def _schale_f(n, mag):
    """Abweichung von der glatten Füllung zwischen zwei Schließungen (Form nach Myers und Swiatecki)."""
    n = np.atleast_1d(np.asarray(n, float))
    aus = np.zeros_like(n)
    for a, b in zip(mag[:-1], mag[1:]):
        q = 0.6 * (b**(5 / 3) - a**(5 / 3)) / (b - a)
        m = (n >= a) & (n < b)
        aus[m] = q * (n[m] - a) - 0.6 * (n[m]**(5 / 3) - a**(5 / 3))
    return aus


@lru_cache(maxsize=None)
def _fuellung(kappa, mu, n_max=12, breite=1.2):
    """Summe der tiefsten n Niveaus des gehaltenen Trichters mit Verwindung [Schalen], n = 0 … 400, das letzte Niveau
    teilweise besetzt; dazu die glatte Füllung aus der über breite Schalen gemittelten Niveaudichte (Gauß mit
    Krümmungskorrektur vierter Ordnung, 15/8 − 5u²/2 + u⁴/2), deren Fermi-Kante n Plätze fasst."""
    from scipy.special import erf
    niv = niveaus_verwindung(kappa, n_max=n_max, mu=mu)
    S, s = [0.0], 0.0
    for E, N, l, j2 in niv:
        for _ in range(j2 + 1):
            s += E; S.append(s)
    S = np.array(S[:401])
    e = np.array([v[0] for v in niv]); d = np.array([v[3] + 1 for v in niv], float)
    def zahl_energie(lam):
        u = (lam - e[:, None]) / breite
        g = np.exp(-u**2) / math.sqrt(math.pi)
        n_k = 0.5 * (1 + erf(u)) + g * (7 / 8 * u - u**3 / 4)                     # Plätze bis lam
        e_k = e[:, None] * n_k - breite * g / 2 * (3 / 8 - 3 / 2 * u**2 + u**4 / 2)  # Energie bis lam
        return (d[:, None] * n_k).sum(0), (d[:, None] * e_k).sum(0)
    lam = np.linspace(e.min() - 4 * breite, n_max - 2.0, 20001)
    nz, ez = zahl_energie(lam)
    glatt = np.interp(np.arange(len(S)), nz, ez)
    return S, glatt


NIVEAUS = ("niveaus", 0.12, 0.03)           # Schalenenergie aus den Niveaus mit Verwindung, Abstand aus der Halterlänge


def schalen_niveaus(n, kappa=0.12, mu=0.03):
    """Schalenenergie einer Sorte [Schalen]: Füllung der Niveaus weniger glatte Füllung. In der Bindung zählt nur der
    bindende Teil der Summe beider Sorten: Wo sie nicht bindet, verformt sich der Kern und nimmt ihn weg."""
    S, glatt = _fuellung(kappa, mu)
    n = np.atleast_1d(np.asarray(n, int))
    return S[n] - glatt[n]


def schalenabstand(lam):
    """Abstand der Schalen ħω·A^{1/3} [MeV] des gehaltenen Trichters R = λ_H·A^{1/3}: gleicher mittlerer Radius wie die
    gefüllten Schalen, ħω = (5/6)(3/2)^{4/3}·ħ²/(m_N λ_H²)·A^{−1/3}."""
    hc = K.hbar * K.c / MEV * 1e15
    mn = K.m_p * K.c**2 / MEV
    return 5 / 6 * 1.5**(4 / 3) * hc**2 / (mn * lam**2)


# ---------------------------------------------------------------- gemessene Kerne
@lru_cache(maxsize=None)
def b_gemessen():
    """Gemessene Bindungen (AME2020, ohne geschätzte) [MeV] je (Z, A)."""
    return {(z, a): b * a / 1e3 for z, n, a, b, g in D.ame2020_markiert() if g and z >= 1}


@lru_cache(maxsize=None)
def _feld():
    g = [(z, a, b) for (z, a), b in b_gemessen().items() if a >= 12]
    return tuple(np.array(x, float) for x in zip(*g))


# ---------------------------------------------------------------- Bindung
def bindung(Z, A, f, p=None):
    """Bindung [MeV] nach Formel f; p: Parameter (sonst die angepassten)."""
    p = anpassung(f)[0] if p is None else p
    Z = np.atleast_1d(np.asarray(Z, float)); A = np.atleast_1d(np.asarray(A, float)); N = A - Z
    par = np.where((Z % 2 == 0) & (N % 2 == 0), 1, np.where((Z % 2 == 1) & (N % 2 == 1), -1, 0))
    B = p[0] * A - p[1] * A**(2 / 3) - p[2] * Z * (Z - 1) / A**(1 / 3) - p[3] * (N - Z)**2 / A + p[4] * par / np.sqrt(A)
    if f[0] == "schale":
        mn, mz = magie(f)
        sh = (_schale_f(N, mn) + _schale_f(Z, mz)) / (A / 2)**(2 / 3)
        B = B - p[5] * np.abs(N - Z) + p[6] * (N - Z)**2 / A**(4 / 3) - p[7] * (sh - p[8] * A**(1 / 3))
    if f[0] == "niveaus":
        sh = np.minimum(schalen_niveaus(N, f[1], f[2]) + schalen_niveaus(Z, f[1], f[2]), 0) / A**(1 / 3)
        hw = p[7] if len(f) > 3 else schalenabstand(halterlaenge(a_C(TROPFEN)))
        B = B - p[5] * np.abs(N - Z) + p[6] * (N - Z)**2 / A**(4 / 3) - hw * sh
    return B


@lru_cache(maxsize=None)
def anpassung(f):
    """Einmalige Anpassung an die gemessenen Kerne (A ≥ 12): Parameter, Streuung [MeV], Zahl der Kerne."""
    from scipy.optimize import least_squares
    Z, A, B = _feld()
    start = [15.8, 18, 0.7, 23, 12] + ([0.5, 30, 3, 0.3] if f[0] == "schale" else [0.5, 30] if f[0] == "niveaus" else [])
    start += [40] if f[0] == "niveaus" and len(f) > 3 else []
    p = least_squares(lambda q: bindung(Z, A, f, q) - B, start).x
    r = bindung(Z, A, f, p) - B
    return tuple(p), float(np.sqrt(np.mean(r**2))), len(B)


def a_S(f):
    return anpassung(f)[0][1]


def a_C(f):
    return anpassung(f)[0][2]


def halterlaenge(aC):
    """λ_H = (3/5)·αħc/a_C [fm] aus dem Coulomb-Glied a_C [MeV]."""
    return 0.6 * K.alpha * K.hbar * K.c / MEV * 1e15 / aC


def masse(Z, A, f):
    """Faltungsenergie des Atoms [MeV]: Z·E_H + N·E_n − B."""
    Z = np.asarray(Z, float)
    return Z * E_H + (A - Z) * E_N - bindung(Z, A, f)


def x(Z, A, f):
    """Spaltbarkeit x = E_C/(2E_S) = a_C·Z²/(2a_S·A)."""
    return a_C(f) * Z**2 / (2 * a_S(f) * A)


def q_alpha(Z, A, f):
    """Energie des Alpha-Zerfalls [MeV]: B(Z−2, A−4) + B(⁴He, gemessen) − B(Z, A)."""
    return float(bindung(Z - 2, A - 4, f)[0] + b_gemessen()[(2, 4)] - bindung(Z, A, f)[0])


def s_n(Z, A, f):
    """Neutronen-Ablösung [MeV]."""
    return float(bindung(Z, A, f)[0] - bindung(Z, A - 1, f)[0])


def pruefung(groesse, zmin, f):
    """Formel gegen Messung: Zahl, mittlere Abweichung, Streuung [MeV] für Q_α ('qa') oder S_n ('sn') ab Z = zmin."""
    g = b_gemessen()
    if groesse == "qa":
        d = [q_alpha(z, a, f) - (g[(z - 2, a - 4)] + g[(2, 4)] - b) for (z, a), b in g.items() if z >= zmin and (z - 2, a - 4) in g]
    else:
        d = [s_n(z, a, f) - (b - g[(z, a - 1)]) for (z, a), b in g.items() if z >= zmin and (z, a - 1) in g]
    return len(d), sum(d) / len(d), math.sqrt(sum(v * v for v in d) / len(d))


# ---------------------------------------------------------------- Tal der Stabilität und Kippstelle
def z_tal(A, f):
    """Tal: Kernladung kleinster Masse bei festem A, kontinuierlich."""
    from scipy.optimize import minimize_scalar
    return minimize_scalar(lambda z: float(masse(z, A, f)[0]), bounds=(0.25 * A, 0.55 * A), method="bounded",
                           options=dict(xatol=1e-6)).x


@lru_cache(maxsize=None)
def a_tal(Z, f):
    """Massenzahl im Tal bei Z (Umkehrung von z_tal)."""
    from scipy.optimize import brentq
    return brentq(lambda a: z_tal(a, f) - Z, max(2.0 * Z - 5, 1.0), 4.5 * Z + 10, xtol=1e-6)


def n_tal(Z, f):
    return a_tal(Z, f) - Z


@lru_cache(maxsize=None)
def tal_ganz(f, a_max=800):
    """Tal, je A das ganze Z mit kleinster Masse: {A: Z}."""
    aus = {}
    for A in range(20, a_max):
        zs = np.arange(1, min(A, 300))
        aus[A] = int(zs[int(np.argmin(masse(zs, A, f)))])
    return aus


@lru_cache(maxsize=None)
def kippstelle(f):
    """Kippstelle x = 1 im Tal, kontinuierlich: (Z, A)."""
    from scipy.optimize import brentq
    a = brentq(lambda a: x(z_tal(a, f), a, f) - 1, 150, 800, xtol=1e-4)
    return z_tal(a, f), a


@lru_cache(maxsize=None)
def kippstelle_ganz(f):
    """Zweiter Weg: erstes ganzes Z im Tal (je A das ganze Z kleinster Masse) mit x ≥ 1, und dessen A-Bereich."""
    t = tal_ganz(f)
    erste = next(t[A] for A in sorted(t) if x(t[A], A, f) >= 1)
    As = [A for A in sorted(t) if t[A] == erste and x(erste, A, f) >= 1]
    return erste, min(As), max(As)


@lru_cache(maxsize=None)
def tal_vergleich(f):
    """Leichtestes gemessenes Isobar je A (20 ≤ A ≤ 238) gegen das ganzzahlige Tal: mittlere |ΔZ|, Anteil genau, Zahl."""
    leicht = {}
    for (z, a), b in b_gemessen().items():
        m = z * E_H + (a - z) * E_N - b
        if 20 <= a <= 238 and (a not in leicht or m < leicht[a][1]):
            leicht[a] = (z, m)
    t = tal_ganz(f)
    d = np.abs(np.array([leicht[a][0] - t[a] for a in leicht]))
    return float(d.mean()), float(np.mean(d == 0)), len(d)


# ---------------------------------------------------------------- Hülle
def huellen_grenze(f):
    """Kernladung, bei der die 1s-Bindung (Dirac) die doppelte Faltungsenergie des Elektrons 6c⁴·Δs_e/G erreicht,
    um den Kern mit R = λ_H·A^{1/3} im Tal der Formel f."""
    lam = halterlaenge(a_C(f))
    return K.z_kritisch(lambda Z: lam * a_tal(Z, f)**(1 / 3) / LAM_E)


@lru_cache(maxsize=None)
def konfigurationen():
    """Gemessene Konfigurationen (PubChem) gegen die Füllung nach n + l: Edelgase, Abweichungen (Z, Art), Zahl der Elemente.
    Art: 'schliesst' (d⁵, d¹⁰, f⁷, f¹⁴ erreicht), 'd_vor_f', 's_an_d'."""
    gem = belegung_gemessen()
    zeilen = _pubchem_zeilen()
    edel = tuple(int(z["AtomicNumber"]) for z in zeilen if z["GroupBlock"] == "Noble gas")
    abw = []
    for Z in sorted(gem):
        ist, soll = gem[Z], belegung_nl(Z)
        if ist == soll:
            continue
        mehr = {k: ist.get(k, 0) - soll.get(k, 0) for k in set(ist) | set(soll) if ist.get(k, 0) != soll.get(k, 0)}
        if any((l == 2 and ist.get((n, l)) in (5, 10)) or (l == 3 and ist.get((n, l)) in (7, 14)) for n, l in mehr):
            art = "schliesst"
        elif any(l == 3 for n, l in mehr):
            art = "d_vor_f"
        else:
            art = "s_an_d"
        abw.append((Z, art))
    return edel, tuple(abw), len(zeilen)


def belegung_nl(Z):
    """Belegung (n, l) → Zahl nach n + l."""
    return {(n, l): min(Z, b) - a + 1 for n, l, a, b in K.fuellung(8) if a <= Z}


@lru_cache(maxsize=None)
def _pubchem_zeilen():
    with open(os.path.join(D.ORT, "elemente", "pubchem_periodensystem.csv"), encoding="utf-8") as f:
        return tuple(csv.DictReader(f))


@lru_cache(maxsize=None)
def belegung_gemessen():
    """Gemessene Belegung (n, l) → Zahl je Element (PubChem)."""
    cfg = {z["Symbol"]: z["ElectronConfiguration"] for z in _pubchem_zeilen()}

    def belegung(s):
        s = re.sub(r"\(.*?\)", "", s).strip()
        m = re.match(r"\[(\w+)\]\s*(.*)", s)
        occ = {}
        if m:
            occ, s = dict(belegung(cfg[m.group(1)])), m.group(2)
        for n, l, k in re.findall(r"(\d+)([spdfg])(\d+)", s):
            occ[(int(n), "spdfg".index(l))] = occ.get((int(n), "spdfg".index(l)), 0) + int(k)
        return {k: v for k, v in occ.items() if v}
    return {int(z["AtomicNumber"]): belegung(z["ElectronConfiguration"]) for z in _pubchem_zeilen()}


# ---------------------------------------------------------------- Lebensdauer
@lru_cache(maxsize=None)
@lru_cache(maxsize=None)
def _alpha_strahler():
    """Alle gerade-geraden Alpha-Strahler ab Z = 84 mit gemessener Halbwertszeit (NUBASE2020) und gemessenen Massen
    (AME2020): (Z, A, Q_α [MeV], lg T_α [s])."""
    g = b_gemessen()
    X = []
    for (Z, A), (T, arten) in D.nubase().items():
        if Z >= 84 and Z % 2 == 0 and (A - Z) % 2 == 0 and arten.get("A", 0) > 0 and (Z, A) in g and (Z - 2, A - 4) in g:
            Q = g[(Z - 2, A - 4)] + g[(2, 4)] - g[(Z, A)]
            if Q > 0:
                X.append((Z, A, Q, math.log10(T * 100 / arten["A"])))
    return tuple(X)


def alpha_tunnel(neutronen=None):
    """Alpha-Zerfall als Tunnel durch den Rand des gehaltenen Trichters, ohne angepasste Zahl: Rand R = λ_H·(A_d^{1/3} + 4^{1/3})
    mit λ_H aus dem Coulomb-Glied, Wall 2Z_d·αħc/r bis zum Austritt bei Q_α, Anlauf v/(2R).
    Gibt mittleren Versatz gerechnet − gemessen [Dekaden], Streuung um ihn bzw. Fehler des Mittels bei neutronen [Dekaden],
    Zahl; neutronen = (N_min, N_max): nur die Strahler in diesem Bereich der Neutronenzahl."""
    from scipy.integrate import quad
    lam = halterlaenge(a_C(TROPFEN))
    hc = K.hbar * K.c / MEV * 1e15                                   # MeV·fm
    mu = K.m_u * K.c**2 / MEV
    d = []
    for Z, A, Q, lgT in _alpha_strahler():
        if neutronen and not neutronen[0] <= A - Z <= neutronen[1]:
            continue
        m = mu * 4 * (A - 4) / A
        R = lam * ((A - 4)**(1 / 3) + 4**(1 / 3))
        b = 2 * (Z - 2) * K.alpha * hc / Q
        I = quad(lambda r: math.sqrt(max(2 * m * (2 * (Z - 2) * K.alpha * hc / r - Q), 0.0)), R, b, limit=200)[0] / hc
        anlauf = K.c * 1e15 * math.sqrt(2 * Q / m) / (2 * R)
        d.append(math.log10(math.log(2) / anlauf) + 2 * I / math.log(10) - lgT)
    d = np.array(d)
    return float(d.mean()), float(d.std() / (math.sqrt(len(d)) if neutronen else 1)), len(d)


def alpha_lebensdauer():
    """Tunnelform des Alpha-Zerfalls, lg T_α = a·Z/√Q_α + b·Z + c (T in s), angepasst an alle gerade-geraden Alpha-Strahler
    ab Z = 84 mit gemessener Halbwertszeit (NUBASE2020) und gemessenen Massen (AME2020).
    Gibt (a, b, c), Streuung [Dekaden], Zahl, Abweichung bei U-238 [Dekaden]."""
    g = b_gemessen()
    X = np.array([(Z, Q, lgT) for Z, A, Q, lgT in _alpha_strahler()])
    M = np.c_[X[:, 0] / np.sqrt(X[:, 1]), X[:, 0], np.ones(len(X))]
    p, *_ = np.linalg.lstsq(M, X[:, 2], rcond=None)
    rest = M @ p - X[:, 2]
    qu = g[(90, 234)] + g[(2, 4)] - g[(92, 238)]
    tu = D.nubase()[(92, 238)]
    du = lg_t_alpha(92, qu, tuple(p)) - math.log10(tu[0] * 100 / tu[1]["A"])
    return tuple(p), float(np.sqrt(np.mean(rest**2))), len(X), float(du)


def lg_t_alpha(Z, Q, p=None):
    a, b, c = alpha_lebensdauer()[0] if p is None else p
    return a * Z / math.sqrt(Q) + b * Z + c


@lru_cache(maxsize=None)
def spaltung_gegen_x(f):
    """Gerade-gerade Kerne ab Z = 90 mit gemessener Spaltung (NUBASE2020): Zahl, Korrelation lg T_Spalt mit x,
    Streuung um die beste Gerade [Dekaden], größtes gemessenes x."""
    X = np.array([(x(Z, A, f), math.log10(T * 100 / arten["SF"])) for (Z, A), (T, arten) in D.nubase().items()
                  if Z >= 90 and Z % 2 == 0 and (A - Z) % 2 == 0 and arten.get("SF", 0) > 0])
    p = np.polyfit(X[:, 0], X[:, 1], 1)
    return len(X), float(np.corrcoef(X[:, 0], X[:, 1])[0, 1]), float(np.sqrt(np.mean((np.polyval(p, X[:, 0]) - X[:, 1])**2))), float(X[:, 0].max())


# ---------------------------------------------------------------- bekannte Kerne
@lru_cache(maxsize=None)
def bekannte():
    """Gemessene Kerne: kleinste Masse je Nukleon, größte Bindung je Nukleon (A ≥ 12), schwerstes gemessenes Z;
    schwerstes Element in AME2020 mit seinem festesten Eintrag und dessen x; x am Maximum der Bindung je Nukleon (Tröpfchenformel)."""
    g = {k: v for k, v in b_gemessen().items() if k[1] >= 12}
    m_min = min(g, key=lambda k: (k[0] * E_H + (k[1] - k[0]) * E_N - g[k]) / k[1])
    b_max = max(g, key=lambda k: g[k] / k[1])
    alle = [(z, a, b) for z, n, a, b, gem in D.ame2020_markiert() if z >= 1]
    z_end = max(z for z, a, b in alle)
    a_end = max((v for v in alle if v[0] == z_end), key=lambda v: v[2])[1]
    bb = max(((float(bindung(z, a, TROPFEN)[0]) / a, z, a) for a in range(21, 121, 2) for z in range(1, a)))   # ungerades A: ohne Paarglied
    return dict(m_min=m_min, b_max=b_max, b_je_a=g[b_max] / b_max[1], z_gem_max=max(z for z, a in b_gemessen()),
                ende=(z_end, a_end, float(x(z_end, a_end, TROPFEN))), x_max=float(x(bb[1], bb[2], TROPFEN)))


# ---------------------------------------------------------------- neue Elemente
def fusionsluecke(Z, f, zmin=20):
    """Fehlende Neutronen, wenn zwei Kerne aus dem Tal zu Z verschmelzen: N_Tal(Z) − max(N_Tal(z₁) + N_Tal(Z − z₁))."""
    return n_tal(Z, f) - max(n_tal(z, f) + n_tal(Z - z, f) for z in range(zmin, Z - zmin + 1))


@lru_cache(maxsize=None)
def grenzen():
    """Grenzen der Elemente: Kippstelle beider Formeln (Schalenglied mit nächster Protonen-Schließung 114, 120, 126),
    Grenze der Hülle beider Formeln, gesichertes und letztes mögliches Z."""
    kt = kippstelle(TROPFEN)
    ks = {p: kippstelle(("schale", p)) for p in (114, 120, 126)}
    hg = (huellen_grenze(TROPFEN), huellen_grenze(SCHALE))
    z_sicher = int(min(kt[0], min(v[0] for v in ks.values())))
    z_max = int(min(min(hg), max(v[0] for v in ks.values())))
    return dict(tropfen=kt, schale=ks, huelle=hg, z_sicher=z_sicher, z_max=z_max)


@lru_cache(maxsize=None)
def streuung_ab(f, zmin=90, verformt=None):
    """Streuung der Bindung [MeV] der gemessenen Kerne ab Z = zmin; verformt True/False: nur die Kerne in den verformten
    Gebieten 150 ≤ A < 190 und A ≥ 220 bzw. nur die übrigen."""
    Z, A, B = _feld()
    m = Z >= zmin
    if verformt is not None:
        m &= (((A >= 150) & (A < 190)) | (A >= 220)) == verformt
    return float(np.sqrt(np.mean((bindung(Z[m], A[m], f) - B[m])**2)))


def abbruchkante(Z, f):
    """Neutronen-Abbruchkante bei Z: größtes A, bis zu dem das nächste Neutronenpaar gebunden ist (S_n zweier Schritte > 0)."""
    A = int(a_tal(Z, f))
    while s_n(Z, A + 1, f) + s_n(Z, A + 2, f) > 0:
        A += 1
    return A


def entstehungsweg(f, z_bereich=range(92, 200)):
    """Weg der Entstehung durch Anlagern von Neutronen: zwischen Tal und Abbruchkante; je Z die Spaltbarkeit x im Tal und an der
    Kante. Gibt das letzte Z, bei dem die Kante noch gehalten ist (x < 1), und je Z (N_Tal, N_Kante, x_Tal, x_Kante)."""
    weg = {}
    for Z in z_bereich:
        At, Ak = a_tal(Z, f), abbruchkante(Z, f)
        weg[Z] = (At - Z, Ak - Z, float(x(Z, At, f)), float(x(Z, Ak, f)))
    return max(Z for Z, w in weg.items() if w[3] < 1), weg


def ende_brennleiter(f):
    """Im Tal die Kernladung, bei der x = ¼ (größte Bindung je Nukleon): (Z, A)."""
    from scipy.optimize import brentq
    z = brentq(lambda z: x(z, a_tal(z, f), f) - 0.25, 5, 80, xtol=1e-6)
    return z, a_tal(z, f)
