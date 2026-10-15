"""Etappe 5: andere Fassungen des Austauschs und der Korrelation im Vergleich, an den Fehltreffern: Ausrichtung und ohne
Selbstwirkung, exakter Austausch (Hartree-Fock, kugelgemittelt), Paarkorrelation zweiter Ordnung, PBE."""
import math
import time
import numpy as np
import kern as K
from huelle_1_loeser import FM, _aussen, _ck, _eigen, _gitter, _hartree, _matrix, _scf, _skalierung, _xc, c, kerne_belegung_nl, _xc_korrelation
from huelle_2_konfiguration import huelle_ergebnis
from huelle_ablage import _zwischen, gespeichert, lesen, merke, schreiben


FASSUNGEN = {"grund": (False, False), "ausrichtung": (True, False), "ohne_selbst": (False, True), "beides": (True, True)}
SKALIERT = False      # Takt aus dem Tempo: jede äußere Bahn ε/(1 + S) (für ein Elektron genau Dirac); der Rumpf bleibt im
                      # festen Takt: bei festem K rückt er mit dem Feld genau wie bei Dirac, ε/(1 + S) verfälschte das
_SK = lambda: "_sk" if SKALIERT else ""
_SKZ = lambda Z: set(_aussen(Z)) if SKALIERT else False
KORRELATION = True    # angepasste Korrelation (Perdew-Zunger) in beiden Zählweisen; False: nur die Zählung der Sperre, örtlich und je Sinn
_KO = lambda: "" if KORRELATION else "_ok"
PAARUNG = False       # Einrasten gegensinniger Windungen (06) in beiden Zählweisen, aus dem Kreis der Entstehung
EINRASTEN = "band"    # "band": alle gegensinnigen Paare im Kopplungsfenster; "mode": je Mode ein Paar (jeder Zyklus rastet einmal ein)

def _fassung_ein(arg):
    Z, R, gem, ger = arg
    r, h = _gitter()
    _, _, Ve0, Kf, _, sch0 = _scf(Z, kerne_belegung_nl(Z), r, h, R)      # festes K aus der Füllung nach n + l, wie in der Suche
    aus = {}
    for name, (term, selbst) in FASSUNGEN.items():
        Eg = _scf(Z, gem, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0, term=term, selbst=selbst)[1]
        Er = _scf(Z, ger, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0, term=term, selbst=selbst)[1]
        aus[name] = round((Eg - Er) * K.alpha ** 2 * K.energie(K.m_e) / K.e, 4)   # eV; < 0: gemessene tiefer
    return Z, aus

def fassungen_vergleich(z_liste, prozesse=4):
    """Für Elemente, deren gemessene Konfiguration die Rechnung nicht trifft: E(gemessen) − E(gerechnet) in jeder Fassung
    (Ausrichtung der Windungen, ohne Selbstwirkung, beides); gespeichert unter '_fassungen', jedes Ergebnis sofort."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem, konf = kerne.belegung_gemessen(), huelle_ergebnis()["konf"]
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, gem[Z], konf[Z]) for Z in z_liste]
    with Pool(prozesse) as p:
        for Z, erg in p.imap_unordered(_fassung_ein, auftr):
            alle = lesen()
            alle.setdefault("_fassungen", {})[str(Z)] = erg
            schreiben(alle)
            print(Z, erg, flush=True)

# ---------------------------------------------------------------------------------------------------------------
# Fassung „Austausch exakt“: der Tausch der Windungen exakt (Hartree-Fock, kugelgemittelt je Spin), so dass das eigene Feld
# genau null auf sich selbst wirkt; nur die Korrelation örtlich (Perdew-Zunger), ihr Anteil jedes Elektrons mit sich selbst
# abgezogen. Ohne angepassten Wert.
# ---------------------------------------------------------------------------------------------------------------
_KERN = {}

def _kern_k(r, k):
    """Matrix r<^k / r>^(k+1) auf dem Gitter."""
    key = (len(r), float(r[0]), float(r[-1]), k)
    if key not in _KERN:
        a, b = np.minimum.outer(r, r), np.maximum.outer(r, r)
        _KERN[key] = a ** k / b ** (k + 1)
    return _KERN[key]

def _scf_hf(Z, bel, R, h=0.024, rmax=60.0, maxit=200, Kfest=None, korrelation=True, anteil=1.0, skaliert=False):
    """Selbstkonsistente Hülle mit exaktem Austausch und örtlicher Korrelation; Gesamtenergie [Hartree]."""
    from scipy.linalg import eigh
    from sympy.physics.wigner import wigner_3j
    r, h = _gitter(h, 1e-6, rmax)
    B = r * r * h
    Vn = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r)
    bes = {k: (min(m, 2 * k[1] + 1), m - min(m, 2 * k[1] + 1)) for k, m in bel.items()}
    ls = sorted({l for _, l in bel})
    A3 = {(la, k, lb): float(wigner_3j(la, k, lb, 0, 0, 0)) ** 2 for la in range(5) for lb in range(5) for k in range(9)}
    phi = {}
    Vloc = [np.zeros_like(r), np.zeros_like(r)]
    Kx = {(l, s): np.zeros((len(r), len(r))) for l in ls for s in (0, 1)}
    Kt, E_alt, eps = {}, 0.0, {}
    for it in range(maxit):
        Vm = Vn + 0.5 * (Vloc[0] + Vloc[1])
        Kf = Kfest if Kfest is not None else 1 / (1 - Vm / (2 * c * c))
        summe, neu = 0.0, {}
        for s in (0, 1):
            for l in ls:
                ns = sorted(n for (n, ll) in bel if ll == l)
                d, e, _ = _matrix(Vn + Vloc[s], Kf, r, h, l, -1)
                A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1) - anteil * Kx[(l, s)]
                sigma = -1.6 * Z * Z - 5

                def tiefste(M, k):
                    mu, v = eigh(np.diag(B), M - sigma * np.diag(B), subset_by_index=(len(r) - k, len(r) - 1))
                    o = np.argsort(-mu)
                    v = v[:, o]
                    return sigma + 1 / mu[o], v / np.sqrt(np.sum(v * v * B[:, None], axis=0))
                w, v = tiefste(A, max(ns) - l)
                w, v = list(w), [v[:, j] for j in range(v.shape[1])]
                if (l, s) in Kt:                                  # offene Bahn mit Korrektur, senkrecht zum Rumpf gleicher l
                    na = [n for n in ns if 0 < bes[(n, l)][s] < 2 * l + 1][0]
                    M = A + anteil * Kt[(l, s)]
                    for j in range(na - l - 1):
                        bv = B * v[j]
                        M = M + (10 * Z * Z + 100) * np.outer(bv, bv)
                    wa, va = tiefste(M, 1)
                    w[na - l - 1], v[na - l - 1] = wa[0], va[:, 0]
                for n in ns:
                    neu[(n, l, s)] = v[n - l - 1]
                    summe += bes[(n, l)][s] * w[n - l - 1]
                    eps[(n, l, s)] = w[n - l - 1]
        phi = neu if not phi else {k: phi.get(k, neu[k]) * 0.0 + neu[k] for k in neu}
        q = [sum(bes[(n, l)][s] * phi[(n, l, s)] ** 2 * B for (n, l) in bel) for s in (0, 1)]
        qs = q[0] + q[1]
        VH = _hartree(qs, r)
        dichte = [qq / (4 * math.pi * r ** 3 * h) for qq in q]
        ec, vcu, vcd = _xc_korrelation(*dichte) if korrelation else (np.zeros_like(r),) * 3
        # Austausch: Kern je (l, Spin) aus allen besetzten Bahnen desselben Spins
        Kneu = {}
        for (l, s) in Kx:
            M = np.zeros((len(r), len(r)))
            for (nb, lb) in bel:
                fb = bes[(nb, lb)][s]
                if not fb:
                    continue
                g = phi[(nb, lb, s)] * B
                for k in range(abs(l - lb), l + lb + 1, 2):
                    M += fb * A3[(l, k, lb)] * np.outer(g, g) * _kern_k(r, k)
            Kneu[(l, s)] = M
        # teilweise gefüllte Unterschale: Paare f(f − 1)/2 statt f²/2 (Mittel über alle Determinanten, das eigene Feld
        # hebt sich ganz auf); als Operator c·(J_a − K_a) mit c = (f − 1)·(2l + 1)/(2l) − f
        Ktneu, Eteil = {}, 0.0
        for (n, l) in bel:
            for s in (0, 1):
                f = bes[(n, l)][s]
                if l == 0 or not 0 < f < 2 * l + 1:
                    continue
                g = phi[(n, l, s)] * B
                Ja = np.diag(_hartree(phi[(n, l, s)] ** 2 * B, r) * B)
                Ka = sum(A3[(l, k, l)] * np.outer(g, g) * _kern_k(r, k) for k in range(0, 2 * l + 1, 2))
                cc = (f - 1) * (2 * l + 1) / (2 * l) - f
                Ktneu[(l, s)] = Ktneu.get((l, s), 0.0) + cc * (Ja - Ka)
                Eteil += 0.5 * f * cc * float(phi[(n, l, s)] @ (Ja - Ka) @ phi[(n, l, s)])
        Ex = -0.5 * sum(bes[(na, la)][s] * float(phi[(na, la, s)] @ Kneu[(la, s)] @ phi[(na, la, s)])
                        for (na, la) in bel for s in (0, 1) if bes[(na, la)][s])
        Ec = float(np.sum(ec * qs))
        Esic = 0.0                                           # Korrelation eines Elektrons mit sich selbst abgezogen
        for (n, l) in bel:
            for s in (0, 1):
                f = bes[(n, l)][s]
                if f:
                    q1 = phi[(n, l, s)] ** 2 * B
                    if korrelation:
                        e1, _, _ = _xc_korrelation(q1 / (4 * math.pi * r ** 3 * h), np.zeros_like(r))
                        Esic -= f * float(np.sum(e1 * q1))
        exl = -0.75 * (6 / math.pi) ** (1 / 3) * (np.maximum(dichte[0], 0) ** (4 / 3) + np.maximum(dichte[1], 0) ** (4 / 3)) / np.maximum(dichte[0] + dichte[1], 1e-30)
        vxl = [-(6 / math.pi * np.maximum(dd, 0)) ** (1 / 3) for dd in dichte]
        Exl = float(np.sum(exl * qs))                             # örtlicher Austausch, Anteil 1 − anteil
        Eloc = [VH + vcu + (1 - anteil) * vxl[0], VH + vcd + (1 - anteil) * vxl[1]]
        Eop = anteil * sum(bes[(n, l)][s] * float(phi[(n, l, s)] @ Kt[(l, s)] @ phi[(n, l, s)])
                  for (n, l) in bel for s in (0, 1) if (l, s) in Kt and bes[(n, l)][s] and 0 < bes[(n, l)][s] < 2 * l + 1)
        E = summe - float(np.sum(Vloc[0] * q[0] + Vloc[1] * q[1])) + 0.5 * float(np.sum(VH * qs)) - anteil * Ex + Ec \
            - Eop + anteil * Eteil + (1 - anteil) * Exl
        Vloc = [0.5 * Vloc[s] + 0.5 * Eloc[s] for s in (0, 1)] if it else Eloc
        Kx = {k: 0.5 * Kx[k] + 0.5 * Kneu[k] for k in Kx} if it else Kneu
        Kt = {k: 0.5 * Kt.get(k, 0.0) + 0.5 * Ktneu[k] for k in Ktneu} if it else Ktneu
        if abs(E - E_alt) < 1e-7 and it > 5:
            break
        E_alt = E
    if skaliert:                                             # Takt aus dem Tempo: jede Bahn ε/(1 + S); als Menge nur diese
        E += sum(bes[(n, l)][s] * (eps[(n, l, s)] / (1 + _skalierung(Kf, r, h, l, -1, phi[(n, l, s)])) - eps[(n, l, s)])
                 for (n, l) in bel for s in (0, 1) if bes[(n, l)][s] and (skaliert is True or (n, l) in skaliert))
    _scf_hf.korrelation_selbst = Esic                         # Abzug der Korrelation je Elektron mit sich selbst

    def kx(l, s):                                             # Austauschkern für beliebiges l aus den besetzten Bahnen
        M = np.zeros((len(r), len(r)))
        for (nb, lb) in bel:
            fb = bes[(nb, lb)][s]
            if fb:
                g = phi[(nb, lb, s)] * B
                for k in range(abs(l - lb), l + lb + 1, 2):
                    M += fb * float(wigner_3j(l, k, lb, 0, 0, 0)) ** 2 * np.outer(g, g) * _kern_k(r, k)
        return M
    _scf_hf.zustand = dict(r=r, h=h, Vn=Vn, Vloc=Vloc, Kf=Kf, Kx=kx, phi=phi, bes=bes, E=E + Esic)
    return E + Esic, it

def _hf_paar_ein(arg):
    """E(gemessen) − E(Gegenspielerin) [eV] mit exaktem Austausch; festes K aus der Füllung nach n + l."""
    Z, R, gem = arg
    p = max(n for (n, l) in kerne_belegung_nl(Z) if l == 0)
    geg = dict(gem)
    if gem.get((p, 0), 0) == 2:                              # s² gemessen: Gegenspielerin d^(n+1) s¹
        geg[(p, 0)], geg[(p - 1, 2)] = 1, geg.get((p - 1, 2), 0) + 1
    else:                                                     # s¹ oder s⁰ gemessen: Gegenspielerin mit einem s mehr
        geg[(p, 0)], geg[(p - 1, 2)] = geg.get((p, 0), 0) + 1, geg[(p - 1, 2)] - 1
    geg = {k: v for k, v in geg.items() if v}
    r, h = _gitter(0.024, 1e-6, 60.0)
    Kf = _scf(Z, kerne_belegung_nl(Z), r, h, R)[3]
    Eg = _scf_hf(Z, gem, R, Kfest=Kf)[0]
    Eb = _scf_hf(Z, geg, R, Kfest=Kf)[0]
    return Z, round((Eg - Eb) * K.alpha ** 2 * K.energie(K.m_e) / K.e, 4)

def _hf_gegen_ein(arg):
    """E(gemessen) − E(bisher gerechnet) [eV] mit exaktem Austausch; Spin-Bahn-Glied erster Ordnung aus der Grundrechnung
    derselben Konfigurationen; festes K aus der Füllung nach n + l."""
    Z, R, gem, ger = arg
    r, h = _gitter(0.024, 1e-6, 60.0)
    _, _, Ve0, Kf, _, sch0 = _scf(Z, kerne_belegung_nl(Z), r, h, R)
    so = {}
    for name, b in (("g", gem), ("r", ger)):
        e, es = _scf(Z, b, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0)[:2]
        so[name] = es - e
    Eg = _scf_hf(Z, gem, R, Kfest=Kf)[0] + so["g"]
    Er = _scf_hf(Z, ger, R, Kfest=Kf)[0] + so["r"]
    return Z, round((Eg - Er) * K.alpha ** 2 * K.energie(K.m_e) / K.e, 4)

def hf_gegen(z_liste, prozesse=4):
    """Exakter Austausch mit Spin-Bahn: gemessene gegen bisher gerechnete Konfiguration; < 0 heißt richtig. '_hf_gegen'."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem, konf = kerne.belegung_gemessen(), huelle_ergebnis()["konf"]
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, gem[Z], konf[Z]) for Z in z_liste]
    with Pool(prozesse) as p:
        for Z, d in p.imap_unordered(_hf_gegen_ein, auftr):
            alle = lesen()
            alle.setdefault("_hf_gegen", {})[str(Z)] = d
            schreiben(alle)
            print(Z, d, flush=True)

def hf_paare(z_liste, prozesse=4):
    """Exakter Austausch: gemessene Konfiguration gegen ihre s/d-Gegenspielerin; < 0 heißt richtig. Unter '_hf_paare'."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem = kerne.belegung_gemessen()
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, gem[Z]) for Z in z_liste]
    with Pool(prozesse) as p:
        for Z, d in p.imap_unordered(_hf_paar_ein, auftr):
            alle = lesen()
            alle.setdefault("_hf_paare", {})[str(Z)] = d
            schreiben(alle)
            print(Z, d, flush=True)

def _hf_pruef_ein(arg):
    """Ist die gemessene Konfiguration mit exaktem Austausch die tiefste unter ihren Nachbarn (ein Elektron zwischen den
    äußeren Unterschalen verschoben)? Vorauswahl der Nachbarn mit der Grundrechnung (bis 3 eV über der gemessenen)."""
    Z, R, gem = arg
    t = time.time()
    r, h = _gitter(0.024, 1e-6, 60.0)
    _, _, Ve0, Kf, _, sch0 = _scf(Z, kerne_belegung_nl(Z), r, h, R)
    aussen = _aussen(Z)
    nachbarn = []
    for a in aussen:
        if not gem.get(a):
            continue
        for b in aussen:
            if b != a and gem.get(b, 0) < 2 * (2 * b[1] + 1):
                x = dict(gem)
                x[a] -= 1
                x[b] = x.get(b, 0) + 1
                nachbarn.append({k: v for k, v in x.items() if v})
    ev = K.alpha ** 2 * K.energie(K.m_e) / K.e
    lsda, term = {}, {}
    for x in [gem] + nachbarn:
        e, es = _scf(Z, x, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0)[:2]
        lsda[tuple(sorted(x.items()))] = (e, es)
    g0 = tuple(sorted(gem.items()))
    kand = [x for x in nachbarn if (lsda[tuple(sorted(x.items()))][1] - lsda[g0][1]) * ev < 3.0]
    def hf(x):                                               # exakter Austausch + Spin-Bahn + Ausrichtung der Windungen
        k = tuple(sorted(x.items()))
        dt = _scf(Z, x, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0, term=True)[0] - lsda[k][0]
        return _scf_hf(Z, x, R, Kfest=Kf)[0] + (lsda[k][1] - lsda[k][0]) + dt
    Eg = hf(gem)
    abst = {"".join(f"{n}{'spdfg'[l]}{v}" for (n, l), v in sorted(x.items()) if (n, l) in aussen):
            round((hf(x) - Eg) * ev, 4) for x in kand}
    return Z, {"nachbarn": len(nachbarn), "gerechnet": len(kand), "abstand": abst,
               "tiefste": all(d > 0 for d in abst.values()), "sekunden": round(time.time() - t)}

def hf_pruefung(z_liste, prozesse=4, schluessel="_hf_pruefung"):
    """Prüfstufe für den exakten Austausch an allen gemessenen Elementen; jedes Ergebnis sofort gespeichert."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem = kerne.belegung_gemessen()
    fertig = lesen().get(schluessel, {})
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, gem[Z]) for Z in sorted(z_liste, reverse=True)
             if str(Z) not in fertig]
    with Pool(prozesse) as p:
        for Z, erg in p.imap_unordered(_hf_pruef_ein, auftr):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = erg
            schreiben(alle)
            print(Z, erg["tiefste"], erg["gerechnet"], erg["sekunden"], "s", min(erg["abstand"].values(), default=None), flush=True)

# ---------------------------------------------------------------------------------------------------------------
# Korrelation der Paare, zweite Ordnung: aus jeder besetzten Determinante (ganze Besetzung je m und Spin) springen zwei
# Elektronen in unbesetzte Bahnen; E2 = Σ |⟨ij||ab⟩|² / (ε_i + ε_j − ε_a − ε_b). Bahnen und Energien aus der Rechnung mit
# exaktem Austausch (kugelgemittelt), ohne angepassten Wert.
# ---------------------------------------------------------------------------------------------------------------
def _fock_bahnen(Z, bel, R, l_max, n_virt, h=0.024, rmax=60.0, Kfest=None):
    """Selbstkonsistente Hülle mit exaktem Austausch, dann je l und Spin besetzte und n_virt unbesetzte Bahnen des
    Fock-Operators (Energie, φ). Gibt r, h, Bahnen {(l, s): (w, v)}."""
    from scipy.linalg import eigh
    _scf_hf(Z, bel, R, h=h, rmax=rmax, Kfest=Kfest, korrelation=False)
    st = _scf_hf.zustand
    r, h, B = st["r"], st["h"], st["r"] ** 2 * st["h"]
    aus = {}
    for l in range(l_max + 1):
        for s in (0, 1):
            d, e, _ = _matrix(st["Vn"] + st["Vloc"][s], st["Kf"], r, h, l, -1)
            A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1) - st["Kx"](l, s)
            sigma = -1.6 * Z * Z - 5
            k = n_virt + 6
            mu, v = eigh(np.diag(B), A - sigma * np.diag(B), subset_by_index=(len(r) - k, len(r) - 1))
            o = np.argsort(-mu)
            aus[(l, s)] = (sigma + 1 / mu[o], v[:, o] / np.sqrt(np.sum(v[:, o] ** 2 * B[:, None], axis=0)))
    return r, h, aus

def mp2(Z, bel, R, l_max=4, n_virt=20, Kfest=None, rmax=60.0):
    """Korrelation der Paare in zweiter Ordnung für die Determinante mit ganzer Besetzung (je Unterschale die ersten m,
    Spin nach Hund). Spinbahnen aus den Fock-Bahnen; brute-force über m mit Gaunt-Beiwerten. Gibt E2 [Hartree]."""
    r, h, bahnen = _fock_bahnen(Z, bel, R, l_max, n_virt, Kfest=Kfest, rmax=rmax)
    B = r * r * h
    bes = {k: (min(m, 2 * k[1] + 1), m - min(m, 2 * k[1] + 1)) for k, m in bel.items()}
    bes_m = {}                                                # besetzte Spinbahnen: (n, l, m, s)
    for (n, l), (fu, fd) in bes.items():
        ms = list(range(l, -l - 1, -1))
        for s, f in ((0, fu), (1, fd)):
            for m in ms[:f]:
                bes_m[(n, l, m, s)] = True
    occ, virt = [], []
    for (l, s), (w, v) in bahnen.items():
        for j in range(len(w)):
            n = l + 1 + j
            for m in range(-l, l + 1):
                eintrag = (w[j], v[:, j], l, m, s)
                if (n, l, m, s) in bes_m:
                    occ.append(eintrag)
                elif (n, l) not in bel or j < len(w):        # alle übrigen Zustände als unbesetzt
                    if not (n, l) in bel or (n, l, m, s) not in bes_m:
                        if j < n_virt + (max([nn for (nn, ll) in bel if ll == l], default=l) - l):
                            virt.append(eintrag)
    RK = {}

    def rk(k, a, b, c, d):
        q1, q2 = a[1] * c[1] * B, b[1] * d[1] * B
        y = (np.cumsum(q1 * r ** k) - 0.5 * q1 * r ** k) / r ** (k + 1) + \
            r ** k * (np.cumsum((q1 / r ** (k + 1))[::-1])[::-1] - 0.5 * q1 / r ** (k + 1))
        return float(np.sum(y * q2))

    def zwei(a, b, c, d):                                     # ⟨ab|cd⟩ über Spinbahnen
        if a[4] != c[4] or b[4] != d[4] or a[3] + b[3] != c[3] + d[3]:
            return 0.0
        return sum(_ck(k, a[2], a[3], c[2], c[3]) * _ck(k, d[2], d[3], b[2], b[3]) * rk(k, a, b, c, d)
                   for k in range(abs(a[2] - c[2]), a[2] + c[2] + 1, 2) if abs(d[2] - b[2]) <= k <= d[2] + b[2])
    E2 = 0.0
    for i in range(len(occ)):
        for j in range(i + 1, len(occ)):
            a_, b_ = occ[i], occ[j]
            for x in range(len(virt)):
                for y in range(x + 1, len(virt)):
                    c_, d_ = virt[x], virt[y]
                    if a_[3] + b_[3] != c_[3] + d_[3] or a_[4] + b_[4] != c_[4] + d_[4]:
                        continue
                    if (a_[2] + b_[2] + c_[2] + d_[2]) % 2:
                        continue
                    t = zwei(a_, b_, c_, d_) - zwei(a_, b_, d_, c_)
                    if t:
                        E2 += t * t / (a_[0] + b_[0] - c_[0] - d_[0])
    return E2

def mp2_schnell(Z, bel, R, l_max=4, n_virt=20, Kfest=None, rmax=60.0, aktiv=None):
    """Wie mp2, die Summe über m als Tensorrechnung je Unterschalen-Quartett. aktiv: korrelierte besetzte Unterschalen
    (Standard alle). Anregungen, bei denen beide Elektronen in offene Unterschalen der Determinante gehen, fehlen; sie sind
    die Ausrichtung der Windungen (_term)."""
    r, h, bahnen = _fock_bahnen(Z, bel, R, l_max, n_virt, Kfest=Kfest, rmax=rmax)
    B = r * r * h
    bes = {k: (min(m, 2 * k[1] + 1), m - min(m, 2 * k[1] + 1)) for k, m in bel.items()}
    offen = {k for k, (u, d) in bes.items() if 0 < u + d < 2 * (2 * k[1] + 1)}
    aktiv = set(aktiv or bel)
    nmax = {l: max([n for (n, ll) in bel if ll == l], default=l) for l in range(l_max + 1)}
    # besetzte Spinbahnen-Gruppen: (n, l, s) -> Maske der besetzten m (Hund: zuerst m = l, l−1, …)
    def maske(l, f):
        x = np.zeros(2 * l + 1, bool)
        x[:f] = True                                            # Index 0 ↔ m = l
        return x
    occ = {(n, l, s): maske(l, bes[(n, l)][s]) for (n, l) in bel for s in (0, 1) if bes[(n, l)][s] and (n, l) in aktiv}
    # freie Radialbahnen je (l, s): alle bis n_virt über der höchsten besetzten, dazu offene Unterschalen mit ihren freien m
    frei = {}
    for (l, s), (w, v) in bahnen.items():
        for j in range(min(len(w), nmax[l] - l + n_virt)):
            n = l + 1 + j
            f = bes.get((n, l), (0, 0))[s]
            if f < 2 * l + 1:
                frei[(n, l, s)] = ~maske(l, f)
    def bahn(n, l, s):
        w, v = bahnen[(l, s)]
        return w[n - l - 1], v[:, n - l - 1]
    def ck_mat(k, la, lb):                                      # C[ma, mb] = c^k(la ma, lb mb), Index 0 ↔ m = l
        return np.array([[_ck(k, la, ma, lb, mb) for mb in range(lb, -lb - 1, -1)] for ma in range(la, -la - 1, -1)])
    def yk(q, k):
        return (np.cumsum(q * r ** k) - 0.5 * q * r ** k) / r ** (k + 1) + \
            r ** k * (np.cumsum((q / r ** (k + 1))[::-1])[::-1] - 0.5 * q / r ** (k + 1))
    E2 = 0.0
    okeys = sorted(occ)
    fkeys = sorted(frei)
    for ii, I in enumerate(okeys):
        for J in okeys[ii:]:
            (ni, li, si), (nj, lj, sj) = I, J
            ei, pi = bahn(*I)
            ej, pj = bahn(*J)
            for A_ in fkeys:
                for B_ in fkeys:
                    (na, la, sa), (nb, lb, sb) = A_, B_
                    if (sa, sb) != (si, sj) or (li + lj + la + lb) % 2:
                        continue
                    if ((na, la) in offen and (nb, lb) in offen):
                        continue
                    ea, pa = bahn(*A_)
                    eb, pb = bahn(*B_)
                    den = ei + ej - ea - eb
                    # direkt ⟨ij|ab⟩: Σ_k c^k(i,a) c^k(b,j) R^k(ia, jb)
                    def direkt(la_, lb_, pa_, pb_, li_, lj_, pi_, pj_):
                        T = 0.0
                        for k in range(abs(li_ - la_), li_ + la_ + 1, 2):
                            if not abs(lb_ - lj_) <= k <= lb_ + lj_:
                                continue
                            R = float(np.sum(yk(pi_ * pa_ * B, k) * pj_ * pb_ * B))
                            T = T + R * np.einsum("ia,bj->ijab", ck_mat(k, li_, la_), ck_mat(k, lb_, lj_))
                        return T
                    D = direkt(la, lb, pa, pb, li, lj, pi, pj)
                    if isinstance(D, float):
                        continue
                    mi = np.arange(li, -li - 1, -1); mj = np.arange(lj, -lj - 1, -1)
                    ma = np.arange(la, -la - 1, -1); mb = np.arange(lb, -lb - 1, -1)
                    erh = (mi[:, None, None, None] + mj[None, :, None, None] == ma[None, None, :, None] + mb[None, None, None, :])
                    D = D * erh
                    if si == sj:                                # Austausch: ⟨ij|ba⟩
                        X = 0.0
                        for k in range(abs(li - lb), li + lb + 1, 2):
                            if not abs(la - lj) <= k <= la + lj:
                                continue
                            R = float(np.sum(yk(pi * pb * B, k) * pj * pa * B))
                            X = X + R * np.einsum("ib,aj->ijab", ck_mat(k, li, lb), ck_mat(k, la, lj))
                        D = D - (X * erh if not isinstance(X, float) else 0.0)
                    m = occ[I][:, None, None, None] & occ[J][None, :, None, None] & frei[A_][None, None, :, None] & frei[B_][None, None, None, :]
                    if I == J:                                  # gleiche Gruppe: Paare i < j, sonst doppelt
                        m = m & np.triu(np.ones((2 * li + 1, 2 * li + 1), bool), 1)[:, :, None, None]
                    w2 = (D * D * m).sum()
                    if si == sj:
                        w2 *= 0.5                               # a, b geordnet
                    E2 += w2 / den
    return E2

def fehltreffer_befund():
    """Befund zu den Fehltreffern, gespeichert unter '_befund': E(gemessen) − E(Gegenspielerin) bei V und Tb für die Anteile
    0, ½, 1 an exaktem Austausch, die Grenzanteile (Vorzeichenwechsel, linear zwischen den Stützstellen), und für Tb der
    Abstand mit exaktem Austausch ohne Korrelation sowie der Beitrag der Paarkorrelation (l ≤ 5, 15 freie Bahnen)."""
    import kerne
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem = kerne.belegung_gemessen()
    ev = K.alpha ** 2 * K.energie(K.m_e) / K.e
    fall = {23: {(4, 0): 1, (3, 2): 4}, 65: {(4, 3): 8, (5, 2): 1}}
    aus = {}
    for Z, wechsel in fall.items():
        R = lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM
        r, h = _gitter(0.024, 1e-6, 60.0)
        Kf = _scf(Z, kerne_belegung_nl(Z), r, h, R)[3]
        g = gem[Z]
        b = dict(g)
        b.update(wechsel)
        b = {k: v for k, v in b.items() if v}
        d = {a: (_scf_hf(Z, g, R, Kfest=Kf, anteil=a)[0] - _scf_hf(Z, b, R, Kfest=Kf, anteil=a)[0]) * ev for a in (0.0, 0.5, 1.0)}
        st = sorted(d)
        grenze = next((a0 + (a1 - a0) * d[a0] / (d[a0] - d[a1]) for a0, a1 in zip(st, st[1:]) if d[a0] * d[a1] < 0), None)
        aus[str(Z)] = dict(abstand={str(a): round(x, 4) for a, x in d.items()}, grenzanteil=round(grenze, 3))
        if Z == 65:
            hf = (_scf_hf(Z, g, R, Kfest=Kf, korrelation=False)[0] - _scf_hf(Z, b, R, Kfest=Kf, korrelation=False)[0]) * ev
            ak = lambda x: [k for k in x if k[0] >= 5 or k == (4, 3)]
            e2 = (mp2_schnell(Z, g, R, l_max=5, n_virt=15, rmax=15, aktiv=ak(g))
                  - mp2_schnell(Z, b, R, l_max=5, n_virt=15, rmax=15, aktiv=ak(b))) * ev
            aus[str(Z)].update(austausch=round(hf, 4), paarkorrelation=round(e2, 4))
    alle = lesen()
    alle["_befund"] = aus
    schreiben(alle)
    return aus

def _fassung_gradient_ein(arg):
    Z, R, gem, ger = arg
    r, h = _gitter()
    _, _, Ve0, Kf, _, sch0 = _scf(Z, kerne_belegung_nl(Z), r, h, R, gradient=True, maxit=300)
    aus = {}
    for name, term in (("pbe", False), ("pbe_ausrichtung", True)):
        Eg = _scf(Z, gem, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0, term=term, gradient=True, maxit=300)[1]
        Er = _scf(Z, ger, r, h, R, Ve=Ve0, Kfest=Kf, schaetzung=sch0, term=term, gradient=True, maxit=300)[1]
        aus[name] = round((Eg - Er) * K.alpha ** 2 * K.energie(K.m_e) / K.e, 4)
    return Z, aus

def fassungen_gradient(z_liste=(22, 23, 27, 28, 58, 92, 93), prozesse=4):
    """E(gemessen) − E(gerechnet) [eV] mit Austausch-Korrelation nach PBE (Dichte und ihre Steigung), ohne und mit Ausrichtung
    der Windungen; ergänzt in '_fassungen'."""
    import kerne
    from multiprocessing import Pool
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    gem, konf = kerne.belegung_gemessen(), huelle_ergebnis()["konf"]
    auftr = [(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, gem[Z], konf[Z]) for Z in z_liste]
    with Pool(prozesse) as p:
        for Z, erg in p.imap_unordered(_fassung_gradient_ein, auftr):
            alle = lesen()
            alle.setdefault("_fassungen", {}).setdefault(str(Z), {}).update(erg)
            schreiben(alle)


PROBEN = (22, 23, 27, 28, 43, 58, 64, 74, 93, 96,                      # Fehltreffer (gekoppelt oder erste Ordnung)
          92, 78, 41, 79, 44, 24, 29, 45, 57, 40, 47,                     # Grenzfälle der d-Reihen
          59, 65, 66, 67, 68, 97, 98)                                     # Grenzfälle der f-Reihen (kippen mit exaktem Austausch)


def _zw(Z, h=0.016):
    """Zwischenspeicher des Elements auf dem Gitter h (Startfeld aus n + l einmal gerechnet)."""
    import kerne
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    return _zwischen(Z, lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM, h)


def _bahnen(zw, bel):
    """Feld und Bahnen der Konfiguration (Spin-Bahn erster Ordnung), vom gespeicherten Startfeld aus."""
    r, h, R, Z = zw["r"], zw["h"], zw["R"], zw["Z"]
    E0, Es0, Ve, Kf, it, sch = _scf(Z, bel, r, h, R, Ve=zw["Ve0"], Kfest=zw["Kf"])
    return E0, Ve, Kf, sch


def selbstwirkung(Z, bel):
    """Selbstwirkung je Elektron der äußeren Unterschalen im Feld der Konfiguration (Spin-Bahn erster Ordnung):
    U[ρ_i] + E_xc[ρ_i, 0] [eV], ρ_i die Dichte einer Bahn; dazu ⟨r⟩ [a₀]. Was das örtliche Feld je Elektron zu viel rechnet."""
    zw = _zw(Z)
    r, h, R = zw["r"], zw["h"], zw["R"]
    E0, Ve, Kf, sch = _bahnen(zw, bel)
    Vm = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
    aus = {}
    for (n, l) in _aussen(Z):
        if not bel.get((n, l)):
            continue
        w, v = _eigen(Vm, Kf, r, h, l, -1, n - l, -1.6 * Z * Z - 5,
                      [0.5 * (sch[(m, l, 0)] + sch[(m, l, 1)]) for m in range(l + 1, n + 1)])
        q1 = v[:, n - l - 1] ** 2 * r * r * h
        e1, _, _ = _xc(q1 / (4 * math.pi * r ** 3 * h), np.zeros_like(r))
        si = 0.5 * np.sum(_hartree(q1, r) * q1) + np.sum(e1 * q1)
        aus[f"{n}{'spdfg'[l]}"] = dict(si=float(si) * 27.211386, r=float(np.sum(q1 * r)))
    return aus


def selbstwirkung_summe(Z, bel):
    """Σ f·(U[ρ_i] + E_xc[ρ_i, 0]) über die Elektronen der äußeren Unterschalen der Konfiguration [eV]; gespeichert."""
    return merke(_zw(Z), "si", bel, lambda: sum(bel[(int(k[:-1]), "spdfg".index(k[-1]))] * v["si"]
                                                  for k, v in selbstwirkung(Z, bel).items()))


def ausrichtung(Z, bel):
    """Ausrichtung der Windungen (tiefster Term gegen das Mittel der Determinanten) der Konfiguration [eV], aus den Bahnen der
    Rechnung mit Spin-Bahn-Glied erster Ordnung, vom gespeicherten Startfeld; gespeichert."""
    zw = _zw(Z)

    def rechnung():
        e0 = _scf(Z, bel, zw["r"], zw["h"], zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"])[0]
        e1 = _scf(Z, bel, zw["r"], zw["h"], zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"], term=True)[0]
        return (e1 - e0) * 27.211386
    return merke(zw, "term", bel, rechnung)


def konkurrentin(Z):
    """Gemessene Konfiguration und ihre tiefste Konkurrentin unter den gespeicherten gekoppelten Zuständen; Energieabstand [eV]
    (> 0: die gemessene liegt höher)."""
    import kerne
    zw = _zw(Z)
    gem = kerne.belegung_gemessen()[Z]
    E = {}
    for art in ("hund", "jj"):
        for c, v in gespeichert(zw, art).items():
            E[c] = min(E.get(c, 1e99), v)
    kg = " ".join(f"{n}{l}:{v}" for (n, l), v in sorted(gem.items()))
    c, e = min(((c, e) for c, e in E.items() if c != kg), key=lambda x: x[1])
    b = {(int(t.split(":")[0][:-1]), int(t.split(":")[0][-1])): int(t.split(":")[1]) for t in c.split()}
    return gem, b, (E[kg] - e) * 27.211386


def _rk(Pa, Pb, Pc, Pd, r, h, k):
    """Slater-Integral R^k(ab, cd) = ∫∫ P_a P_c (r<^k / r>^(k+1)) P_b P_d, P = u auf dem Gitter (∫P² dr = 1)."""
    w = r * h
    q1, q2 = Pa * Pc * w, Pb * Pd * w
    y = (np.cumsum(q1 * r ** k) - 0.5 * q1 * r ** k) / r ** (k + 1) + r ** k * (np.cumsum((q1 / r ** (k + 1))[::-1])[::-1]
                                                                             - 0.5 * q1 / r ** (k + 1))
    return float(np.sum(y * q2))


def paarmischung(Z, bel):
    """Mischung des geschlossenen Paars ns² mit dem Paar (n − 1)d² (beide Spins gepaart, Gesamtbahn 0): zwei Zustände
    desselben Spins, Kopplung V = R²(ns ns; d d)/√5 aus den Bahnen der Konfiguration, Abstand aus der Rechnung erster Ordnung
    mit dem d²-Paar gepaart. Gibt die Absenkung
    [eV] (≤ 0; 0 ohne Paar oder ohne Platz in der d-Schale); gespeichert."""
    p = max(n for (n, l) in bel if l == 0)
    s, d = (p, 0), (p - 1, 2)
    if bel.get(s) != 2 or bel.get(d, 0) > 8 or p < 4:
        return 0.0
    zw = _zw(Z)

    def rechnung():
        r, h = zw["r"], zw["h"]
        E0, Ve, Kf, sch = _bahnen(zw, bel)
        Vm = np.where(r < zw["R"], -Z / (2 * zw["R"]) * (3 - (r / zw["R"]) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
        P = {}
        for (n, l) in (s, d):
            st = [0.5 * (sch[(m, l, 0)] + sch[(m, l, 1)]) for m in range(l + 1, n + 1)] if (n, l, 1) in sch else None
            w, v = _eigen(Vm, Kf, r, h, l, -1, n - l, -1.6 * Z * Z - 5, st)
            P[(n, l)] = v[:, n - l - 1] * np.sqrt(r)        # u = φ·√r, ∫u²dr = Σφ²r²h = 1
        V = _rk(P[s], P[s], P[d], P[d], r, h, 2) / math.sqrt(5)
        partner = {k: v for k, v in {**bel, s: 0, d: bel.get(d, 0) + 2}.items() if v}
        m = bel.get(d, 0)                                    # Partner mit demselben Spin: das d²-Paar gepaart
        sb = {d: (min(m, 5) + 1, m - min(m, 5) + 1)}
        e0 = _scf(Z, bel, r, h, zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"])[1]
        e1 = _scf(Z, partner, r, h, zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"], spinbes=sb)[1]
        dE = e1 - e0
        return (dE / 2 - math.sqrt((dE / 2) ** 2 + V * V)) * 27.211386
    return merke(zw, "paar2", bel, rechnung)


def austausch_schritt(Z, bel, schale):
    """Fehler des örtlichen Austauschs gegen den exakten für n gleichsinnige Elektronen einer Unterschale, aus deren Bahn in
    der Konfiguration: exakt −½·[n·F⁰ + n(n − 1)·K̄] (K̄ mittlerer Austausch zweier paralleler Elektronen, aus F^k), örtlich
    n^{4/3}·X₁ (X₁ örtlicher Austausch der Dichte einer Bahn). Gibt je n = 1 … 2l + 1 den Zuwachs des Fehlers [eV]."""
    from huelle_1_loeser import _ck
    zw = _zw(Z)
    r, h, R = zw["r"], zw["h"], zw["R"]
    E0, Ve, Kf, sch = _bahnen(zw, bel)
    n_, l = schale
    Vm = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
    w, v = _eigen(Vm, Kf, r, h, l, -1, n_ - l, -1.6 * Z * Z - 5,
                  [0.5 * (sch[(m, l, 0)] + sch[(m, l, 1)]) for m in range(l + 1, n_ + 1)])
    P = v[:, n_ - l - 1] * np.sqrt(r)
    F = {k: _rk(P, P, P, P, r, h, k) for k in range(0, 2 * l + 1, 2)}
    ms = range(-l, l + 1)
    K = sum(_ck(k, l, m, l, m2) ** 2 * F[k] for m in ms for m2 in ms if m != m2 for k in F) / ((2 * l + 1) * 2 * l)
    rho1 = P * P / (4 * math.pi * r * r)
    X1 = -0.75 * (6 / math.pi) ** (1 / 3) * np.sum(rho1 ** (4 / 3) * 4 * math.pi * r ** 3 * h)
    fehler = lambda n: n ** (4 / 3) * X1 + 0.5 * (n * F[0] + n * (n - 1) * K)
    return [(fehler(n + 1) - fehler(n)) * 27.211386 for n in range(1, 2 * l + 1)], K * 27.211386


def _takt_bahn(Vm, Kf, r, h, l, kappa, n, st, sigma):
    """Energie der Bahn (n, l, κ) im Feld Vm; mit SKALIERT im Takt aus dem Tempo, ε/(1 + S)."""
    w, v = _eigen(Vm, Kf, r, h, l, kappa, n - l, sigma, st)
    e = w[n - l - 1]
    return e / (1 + _skalierung(Kf, r, h, l, kappa, v[:, n - l - 1])) if SKALIERT else e


def niveau(Z, bel, h=0.016):
    """Tiefstes Niveau J der Konfiguration statt ihres Mittels: Energie ohne Spin-Bahn-Glied (Spin-Bahn erster Ordnung
    abgezogen) und je offener Unterschale die Spin-Bahn-Lage ihres Niveaus nach Hund (J = |L − S| unter halber Füllung, L + S
    darüber, E = (λ/2)[J(J+1) − L(L+1) − S(S+1)], λ = ±ζ/(2S), ζ aus der Aufspaltung der eigenen Bahn). Gibt
    (E_ohne [eV], Spin-Bahn erster Ordnung [eV], Spin-Bahn des Niveaus [eV]); gespeichert."""
    zw = _zw(Z, h)

    def rechnung():
        r, h, R = zw["r"], zw["h"], zw["R"]
        E, Es, Ve, Kf, it, sch = _scf(Z, bel, r, h, R, Ve=zw["Ve0"], Kfest=zw["Kf"], skaliert=_SKZ(Z), korrelation=KORRELATION)
        Vm = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
        so = 0.0
        for (n, l), m in bel.items():
            if l == 0 or m in (0, 2 * (2 * l + 1)) or m == 2 * l + 1:
                continue
            st = [0.5 * (sch[(k, l, 0)] + sch[(k, l, 1)]) for k in range(l + 1, n + 1)]
            wm, wp = (_takt_bahn(Vm, Kf, r, h, l, ka, n, st, -1.6 * Z * Z - 5) for ka in (l, -l - 1))
            zeta = (wp - wm) / (l + 0.5)
            auf, ab = min(m, 2 * l + 1), m - min(m, 2 * l + 1)
            S = 0.5 * (auf - ab)
            L = abs(sum(range(l, l - auf, -1)) + sum(range(l, l - ab, -1)))
            J = abs(L - S) if m < 2 * l + 1 else L + S
            lam = zeta / (2 * S) if m < 2 * l + 1 else -zeta / (2 * S)
            so += lam / 2 * (J * (J + 1) - L * (L + 1) - S * (S + 1))
        return [E * 27.211386, (Es - E) * 27.211386, so * 27.211386]
    return merke(zw, "niveau" + zw["g"] + _SK() + _KO(), bel, rechnung)


def niveau_gekoppelt(Z, bel, h=0.016):
    """Tiefstes Niveau mit beiden Umlaufsinnen gekoppelt: Abstoßung der Windungen und ζ·l·s aller offenen Unterschalen in einer
    Matrix (ζ aus der Aufspaltung der eigenen Bahn), gegen den tiefsten Term ohne Spin-Bahn [eV]; ersetzt die Lage nach Hund
    in LS. Gespeichert."""
    zw = _zw(Z, h)

    def rechnung():
        from huelle_1_loeser import _term
        r, h, R = zw["r"], zw["h"], zw["R"]
        E, Es, Ve, Kf, it, sch = _scf(Z, bel, r, h, R, Ve=zw["Ve0"], Kfest=zw["Kf"], korrelation=KORRELATION)
        Vm = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
        zeta = {}
        for (n, l), m in bel.items():
            if 0 < m < 2 * (2 * l + 1):
                if l == 0:
                    zeta[(n, l)] = 0.0
                    continue
                st = [0.5 * (sch[(k, l, 0)] + sch[(k, l, 1)]) for k in range(l + 1, n + 1)]
                wm, wp = (_takt_bahn(Vm, Kf, r, h, l, ka, n, st, -1.6 * Z * Z - 5) for ka in (l, -l - 1))
                zeta[(n, l)] = (wp - wm) / (l + 0.5)
        bes = {k: (min(m, 2 * k[1] + 1), m - min(m, 2 * k[1] + 1)) for k, m in bel.items()}
        if not any(zeta.values()):
            return 0.0
        return _term(Z, bel, bes, r, h, R, Ve, Kf, sch, zeta=zeta) * 27.211386
    return merke(zw, "niveau_gekoppelt" + zw["g"] + _SK() + _KO(), bel, rechnung)


def _lese_kurz(s, Z):
    """Konfiguration aus der Kurzschreibweise ab dem Edelgasrumpf (etwa "3d34s1") als {(n, l): Zahl}."""
    import re
    rumpf = max([e for e in (2, 10, 18, 36, 54, 86, 118) if e < Z], default=0)
    b = dict(kerne_belegung_nl(rumpf)) if rumpf else {}
    for n, l, v in re.findall(r"(\d)([spdfg])(\d+?)(?=\d[spdfg]|$)", s):
        b[(int(n), "spdfg".index(l))] = int(v)
    return {k: v for k, v in b.items() if v}


def _haelfte(n, l):
    """Füllung der Spinhälfte (ein Umlaufsinn, 2l + 1 Plätze); eine geschlossene Hälfte zählt ganz."""
    h = 2 * l + 1
    return 1.0 if n in (h, 2 * h) else (n % h) / h


def _gewicht(Z, a, c, gitter=0.016, tempo=None):
    """Gewicht der Zählung je Umlaufsinn beim Wechsel a → c, als Chance zweier Spiegel: Füllung h der Spinhälfte der Schale
    größter Windung, die wechselt (eine geschlossene Hälfte zählt ganz), mal Verhältnis der Tempi v/c = √δ der Bahnen,
    zwischen denen das Elektron wechselt, im Feld von a auf dem Gitter (die Windung ist Phase, linear im Tempo).
    w/(1 − w) = h/(1 − h) · √(δ_hin/δ_weg) (Fermi-Form). tempo: δ je Unterschale {"nl": δ} im Feld von a aus einer anderen
    Rechnung (etwa dem Kreis der Entstehung), sonst spreizung."""
    d = {k for k in set(a) | set(c) if a.get(k, 0) != c.get(k, 0)}
    l = max(k[1] for k in d)
    k = max((x for x in d if x[1] == l), key=lambda x: x[0])
    h = max(_haelfte(a.get(k, 0), l), _haelfte(c.get(k, 0), l))
    if h in (0.0, 1.0):
        return h
    weg = max(d, key=lambda x: a.get(x, 0) - c.get(x, 0))
    hin = max(d, key=lambda x: c.get(x, 0) - a.get(x, 0))
    s = tempo if tempo is not None else spreizung(Z, a, h=gitter)
    vh, vw = (s.get(f"{k[0]}{k[1]}", 0.0) ** 0.5 for k in (hin, weg))   # ungebunden: kein Tempo
    if vw == 0:
        return h if vh == 0 else 1.0
    q = h / (1 - h) * vh / vw
    return q / (1 + q)


def zwei_zaehlweisen(Z, gekoppelt=False, h=0.024):
    """Gemessene Konfiguration gegen jede Nachbarin der Prüfstufe hf_pruefung, nach den zwei Zählweisen (_abstand_zwei).
    Gibt {Konkurrentin: Abstand [eV]}, > 0: die gemessene liegt tiefer."""
    import kerne
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    return {s: _abstand_paar(Z, g, _lese_kurz(s, Z), gekoppelt, h) for s in lesen()["_hf_pruefung"][str(Z)]["abstand"]}


def _energie_exakt(Z, bel, h=0.024):
    """Energie mit exaktem Austausch [eV] (festes K wie in der Prüfstufe hf_pruefung, Gitterschritt h); gespeichert."""
    import kerne
    zw = _zw(Z)

    def rechnung():
        lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
        R = lam * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM
        r, hg = _gitter(h, 1e-6, 60.0)
        Kf = _scf(Z, kerne_belegung_nl(Z), r, hg, R)[3]
        return float(_scf_hf(Z, bel, R, Kfest=Kf, h=h, skaliert=_SKZ(Z), korrelation=KORRELATION)[0]) * 27.211386
    return merke(zw, ("exakt" if h == 0.024 else f"exakt_h{h}") + _SK() + _KO(), bel, rechnung)


def rechne_zwei_zaehlweisen(z_bereich=range(1, 119), prozesse=4, gekoppelt=False, h=0.024, schluessel="_zwei_zaehlweisen"):
    """zwei_zaehlweisen (Lage des Niveaus gekoppelt oder nach Hund, exakter Teil auf dem Gitter h) für alle gemessenen Elemente;
    unter schluessel, jedes Ergebnis sofort gespeichert, schon gerechnete übersprungen. Gibt die Paare mit Abstand ≤ 0."""
    from multiprocessing import Pool
    with Pool(prozesse) as p:
        fertig = lesen().get(schluessel, {})
        for Z, d in p.imap_unordered(_zz_ein, [(Z, gekoppelt, h) for Z in z_bereich if str(Z) not in fertig]):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = d
            schreiben(alle)
    return sorted((int(Z), s, round(v, 3)) for Z, d in lesen()[schluessel].items() for s, v in d.items() if v <= 0)


def _zz_ein(arg):
    return arg[0], zwei_zaehlweisen(*arg)


def _exakt(Z, bel, h):
    """Exakter Teil [eV] auf dem Gitter h; h = 0: gegen das Gitter null extrapoliert aus h = 0,016 und 0,012 (Fehler ∝ h²)."""
    if h:
        return _energie_exakt(Z, bel, h)
    e1, e2 = _energie_exakt(Z, bel, 0.016), _energie_exakt(Z, bel, 0.012)
    return e2 + (e2 - e1) * 0.012 ** 2 / (0.016 ** 2 - 0.012 ** 2)


def einrasten(Z, bel):
    """Bindung der eingerasteten Paare der Konfiguration [eV] (Einrasten gegensinniger Windungen, Zählung EINRASTEN), aus dem
    Kreis der Entstehung; in beiden Zählweisen gleich; gespeichert."""
    import huelle_2_entstehung as K2
    return merke(_zw(Z), "einrasten_" + EINRASTEN, bel,
                 lambda: K2.paar_probe(Z, bel, (K2.N_FREI,), (None,), EINRASTEN)[(K2.N_FREI, None)])


def _abstand_zwei(Z, a, c, gekoppelt=False, h=0.024):
    """E(c) − E(a) nach den zwei Zählweisen (Gewicht _gewicht), am eingerasteten Niveau J [eV]; mit gekoppelt die Lage des
    Niveaus aus niveau_gekoppelt statt nach Hund in LS; mit PAARUNG dazu das Einrasten der Paare, in beiden Zählweisen gleich."""
    Ea, ja, la = niveau(Z, a)
    Ec, jc, lc = niveau(Z, c)
    if gekoppelt:
        la, lc = niveau_gekoppelt(Z, a), niveau_gekoppelt(Z, c)
    if PAARUNG:
        la, lc = la + einrasten(Z, a), lc + einrasten(Z, c)
    w = _gewicht(Z, a, c)
    oertlich = (Ec + lc) - (Ea + la)
    if w == 0:
        return oertlich
    exakt = _exakt(Z, c, h) - _exakt(Z, a, h) + lc - la      # beide Zählweisen kugelgemittelt, ohne Ausrichtung
    return (1 - w) * oertlich + w * exakt


def zaehlweisen_streit():
    """Paare (gemessene gegen Nachbarin der Prüfstufe), in denen örtliche und exakte Zählung verschiedene Vorzeichen geben: je
    Paar die wechselnden Unterschalen mit Füllung vorher → nachher, örtlich, exakt (ohne Ausrichtung), Ausrichtung, Lage des
    Niveaus, Gewicht w und das Gewicht w*, bei dem der Abstand null wäre (w > w*: die exakte Zählung muss gewinnen)."""
    import kerne
    gem = kerne.belegung_gemessen()
    aus = []
    for Z in range(1, 119):
        g = {k: v for k, v in gem[Z].items() if v}
        Eg, jg, lg = niveau(Z, g)
        for s in lesen()["_hf_pruefung"][str(Z)]["abstand"]:
            c = _lese_kurz(s, Z)
            Ec, jc, lc = niveau(Z, c)
            ort = (Ec + lc) - (Eg + lg)
            ex = _energie_exakt(Z, c) - _energie_exakt(Z, g) + lc - lg
            au = ausrichtung(Z, c) - ausrichtung(Z, g)
            if (ort > 0) == (ex + au > 0):
                continue
            w = _gewicht(Z, g, c)
            wechsel = " ".join(f"{n}{'spdfg'[l]}{g.get((n, l), 0)}→{c.get((n, l), 0)}" for (n, l) in sorted(set(g) | set(c))
                               if g.get((n, l), 0) != c.get((n, l), 0))
            aus.append(dict(Z=Z, wechsel=wechsel, oertlich=round(ort, 3), exakt=round(ex, 3), ausrichtung=round(au, 3),
                            niveau=round(lc - lg, 3), w=round(w, 3), w_null=round(ort / (ort - ex - au), 3)))
    return aus


def abstand_gitter(Z, kurz, hs=(0.016, 0.012, 0.008)):
    """Abstand der gemessenen Konfiguration zur Nachbarin kurz, beide Zählweisen und die Lage des Niveaus auf demselben Gitter
    h, gespiegelt um die Mitte; je h und gegen das Gitter null (Fehler ∝ h², aus den beiden feinsten), dazu die Probe, ob der
    Fehler wie h² geht (Verhältnis der Schritte). Gibt {h: Abstand [eV]}."""
    import kerne
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    c = _lese_kurz(kurz, Z)

    def ein(a, b, h):
        Ea, ja, la = niveau(Z, a, h)
        Eb, jb, lb = niveau(Z, b, h)
        la, lb = niveau_gekoppelt(Z, a, h), niveau_gekoppelt(Z, b, h)
        w = _gewicht(Z, a, b, h)
        oertlich = (Eb + lb) - (Ea + la)
        return oertlich if w == 0 else (1 - w) * oertlich + w * (_energie_exakt(Z, b, h) - _energie_exakt(Z, a, h) + lb - la)
    aus = {h: 0.5 * (ein(g, c, h) - ein(c, g, h)) for h in hs}
    h1, h2 = sorted(hs)[:2]
    aus[0] = aus[h1] + (aus[h1] - aus[h2]) * h1 ** 2 / (h2 ** 2 - h1 ** 2)
    return aus


def _abstand_paar(Z, a, c, gekoppelt=False, h=0.024):
    """Abstand E(c) − E(a) eines Paares, gespiegelt um die Mitte: Mittel der beiden Blickrichtungen,
    ½·[_abstand_zwei(a → c) − _abstand_zwei(c → a)]; antisymmetrisch, unabhängig davon, von welcher Seite gesucht wird [eV]."""
    return 0.5 * (_abstand_zwei(Z, a, c, gekoppelt, h) - _abstand_zwei(Z, c, a, gekoppelt, h))


def suche_zwei_zaehlweisen(Z, fenster=3.0, schritte=30, gekoppelt=False, h=0.024):
    """Konfiguration aus den zwei Zählweisen, ohne Kenntnis der Messung: von n + l aus je ein Elektron zwischen den äußeren
    Unterschalen; Nachbarinnen bis fenster [eV] über der jetzigen in der örtlichen Zählung werden nach beiden Zählweisen
    verglichen (Lage des Niveaus gekoppelt oder nach Hund, exakter Teil auf dem Gitter h); weiter zur tiefsten, solange eine
    tiefer liegt."""
    jetzt = kerne_belegung_nl(Z)
    for _ in range(schritte):
        Ej = sum(niveau(Z, jetzt)[::2])
        besser, bd = None, 0.0
        for a in _aussen(Z):
            if not jetzt.get(a):
                continue
            for b in _aussen(Z):
                if b == a or jetzt.get(b, 0) >= 2 * (2 * b[1] + 1):
                    continue
                c = {k: v for k, v in {**jetzt, a: jetzt[a] - 1, b: jetzt.get(b, 0) + 1}.items() if v}
                if sum(niveau(Z, c)[::2]) - Ej > fenster:
                    continue
                d = _abstand_paar(Z, jetzt, c, gekoppelt, h)
                if d < bd - 1e-6:
                    besser, bd = c, d
        if besser is None:
            return jetzt
        jetzt = besser
    return jetzt


def pruefe_gemessene(z_bereich=range(1, 119), fenster=3.0, h=0.024, gekoppelt=False):
    """Prüfung an der Messung: Abstand (_abstand_zwei) der gemessenen Konfiguration zu jeder Nachbarin (ein Elektron zwischen
    den äußeren Unterschalen, bis fenster [eV] in der örtlichen Zählung). Gibt je Z (kleinster Abstand [eV], Nachbarin)."""
    import kerne
    gem = kerne.belegung_gemessen()
    aus = {}
    for Z in z_bereich:
        g = {k: v for k, v in gem[Z].items() if v}
        Eg = sum(niveau(Z, g)[::2])
        best = (None, None)
        for a in _aussen(Z):
            if not g.get(a):
                continue
            for b in _aussen(Z):
                if b == a or g.get(b, 0) >= 2 * (2 * b[1] + 1):
                    continue
                c = {k: v for k, v in {**g, a: g[a] - 1, b: g.get(b, 0) + 1}.items() if v}
                if sum(niveau(Z, c)[::2]) - Eg > fenster:
                    continue
                d = _abstand_paar(Z, g, c, gekoppelt, h)
                if best[0] is None or d < best[0]:
                    best = (d, " ".join(f"{n}{'spdfg'[l]}{v}" for (n, l), v in sorted(c.items()) if (n, l) in _aussen(Z)))
        aus[Z] = best
    return aus


def rechne_pruefung(h=0.016, z_bereich=range(1, 119), prozesse=4, gekoppelt=False, schluessel=None, skaliert=False,
                    korrelation=True, paarung=False, einrasten="band"):
    """pruefe_gemessene (exakter Teil auf dem Gitter h, Lage des Niveaus gekoppelt oder nach Hund) für alle gemessenen
    Elemente; je Element sofort unter schluessel gespeichert, schon geprüfte übersprungen. Gibt die Elemente mit Abstand ≤ 0."""
    from multiprocessing import Pool
    global SKALIERT, KORRELATION, PAARUNG, EINRASTEN
    SKALIERT, KORRELATION, PAARUNG, EINRASTEN = skaliert, korrelation, paarung, einrasten
    schl = schluessel or f"_pruefung_h{h}"
    fertig = lesen().get(schl, {})
    with Pool(prozesse) as p:
        for Z, (d, c) in p.imap_unordered(_pr_ein, [(Z, h, gekoppelt) for Z in z_bereich if str(Z) not in fertig]):
            alle = lesen()
            alle.setdefault(schl, {})[str(Z)] = [d, c]
            schreiben(alle)
            print(Z, d, c, flush=True)
    return sorted(int(Z) for Z, (d, c) in lesen()[schl].items() if d is not None and d <= 0)


def _pr_ein(arg):
    Z, h, gekoppelt = arg
    return Z, pruefe_gemessene([Z], h=h, gekoppelt=gekoppelt)[Z]


def rechne_suche_zwei(z_bereich=range(1, 119), prozesse=4, gekoppelt=False, h=0.024, schluessel="_suche_zwei"):
    """Suche nach den zwei Zählweisen; unter schluessel (jedes Ergebnis sofort gespeichert; schon gespeicherte Elemente
    werden übersprungen, bei geänderter Rechnung ist der Eintrag vorher zu löschen)."""
    from multiprocessing import Pool
    fertig = lesen().get(schluessel, {})
    with Pool(prozesse) as p:
        for Z, b in p.imap_unordered(_sz_ein, [(Z, gekoppelt, h) for Z in z_bereich if str(Z) not in fertig]):
            if b is None:
                continue
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = [[n, l, v] for (n, l), v in sorted(b.items())]
            schreiben(alle)
            print(Z, flush=True)


def rechne_ueber_118(z_start=119, prozesse=4, gekoppelt=False, h=0.024, schluessel="_suche_zwei"):
    """Suche ohne Messung über Z = 118 hinaus, ohne vorgegebenes Ende: je prozesse Elemente zugleich, gespeichert unter
    schluessel; das Ende setzt die Rechnung selbst: Das erste Element, dessen gefundene Konfiguration eine besetzte äußere
    Bahn nicht mehr bindet (keine Spreizung) oder dessen Hülle sich nicht rechnen lässt, gibt es nicht. Gibt dieses Z."""
    Z = z_start
    while True:
        block = list(range(Z, Z + prozesse))
        rechne_suche_zwei(block, prozesse, gekoppelt, h, schluessel)
        fertig = lesen().get(schluessel, {})
        for z in block:
            if str(z) not in fertig:
                return z
            b = {(n, l): v for n, l, v in fertig[str(z)]}
            s = spreizung(z, b)
            if any(f"{n}{l}" not in s for (n, l) in _aussen(z) if b.get((n, l))):
                return z
        Z += prozesse


def _sz_ein(arg):
    Z, gekoppelt, h = arg
    try:
        return Z, suche_zwei_zaehlweisen(Z, gekoppelt=gekoppelt, h=h)
    except Exception as fehler:                              # Hülle nicht rechenbar: kein Eintrag
        print(Z, "nicht rechenbar:", fehler, flush=True)
        return Z, None


def spreizung(Z, bel, spin=False, h=0.016):
    """Spreizung δ = (v/c)² = 2|ε|/(m_ec²) jeder äußeren Unterschale, besetzt oder leer, im Feld der Konfiguration (Spin-Bahn
    erster Ordnung, Gitter h): über die Spins gemittelt, oder mit spin je Spinhälfte (0 die volle Hälfte); gespeichert.
    Gibt {"nl": δ} bzw. {"nl": [δ₀, δ₁]}."""
    zw = _zw(Z, h)

    def rechnung():
        r, h, R = zw["r"], zw["h"], zw["R"]
        E0, Ve, Kf, sch = _bahnen(zw, bel)
        Vk = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r)
        felder = (Vk + Ve[0], Vk + Ve[1]) if spin else (Vk + 0.5 * (Ve[0] + Ve[1]),)
        out = {}
        for (n, l) in _aussen(Z):                            # auch leere Bahnen, im Feld der Konfiguration
            d = [2 * abs(e) * K.alpha ** 2 if e < 0 else None
                 for e in (_eigen(V, Kf, r, h, l, -1, n - l, -1.6 * Z * Z - 5)[0][n - l - 1] for V in felder)]
            if None not in d:
                out[f"{n}{l}"] = d if spin else d[0]
        return out
    return merke(zw, ("spreizung_spin" if spin else "spreizung2") + zw["g"], bel, rechnung)


def takt_abstand(Z, kurz):
    """Wirkung des Takts aus dem Tempo auf den Abstand der gemessenen Konfiguration zur Nachbarin kurz (örtliche Zählung, h =
    0,016, festes K): Energie mit jeder Bahn ε/(1 + S) (für ein Elektron genau Dirac) gegen ε, und mit ε/(1 + S) nur für die
    äußeren Unterschalen (der Rumpf im festen Takt). Gibt (Abstand mit festem Takt, mit Takt aus dem Tempo, nur außen) [eV]."""
    import kerne
    zw = _zw(Z)
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    c = _lese_kurz(kurz, Z)
    E = {}
    for name, b in (("g", g), ("c", c)):
        for sk in (False, True, set(_aussen(Z))):
            E[(name, str(sk))] = _scf(Z, b, zw["r"], zw["h"], zw["R"], Ve=zw["Ve0"], Kfest=zw["Kf"], skaliert=sk)[0] * 27.211386
    return tuple(E[("c", str(sk))] - E[("g", str(sk))] for sk in (False, True, set(_aussen(Z))))
