# -*- coding: utf-8 -*-
"""Ein Lauf der Simulation: alle Szenarien, Prüfung der Modellzahlen, ein nummerierter Bericht.

Aufruf: python simulation/lauf.py
Schreibt berichte/bericht_NNNN.md (nie überschrieben) und berichte/bericht_NNNN/ mit den Bahndaten
(je Szenario eine CSV-Datei und szene.json) für eine spätere Darstellung.
"""
from functools import lru_cache
import datetime
import json
import math
import os
import re
import statistics
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


@lru_cache(maxsize=None)
def _abst_q():
    """Q aus den Abständen allein: DESI DR2 absolut (r_d aus ω_b und Q·ω_b) + Pantheon+-Form, a₀ je Q bestmöglich.
    Bestwert, Δχ² = 1, und Δχ² von 2π − 1."""
    import numpy as np
    from scipy.optimize import brentq, minimize_scalar

    def chi(Q, a):
        om, or_ = _om_or_q(Q)(a * 1e-10)
        return D.bao_absolut(_h0_q(Q, a * 1e-10), om, or_, K.schallhorizont(Q)) + D.pantheon_form(om, or_)

    prof = lambda Q: minimize_scalar(lambda a: chi(Q, a), bounds=(1.0, 1.2), method="bounded", options=dict(xatol=5e-4)).fun
    gr = np.linspace(4.9, 6.1, 7)
    cs = np.array([prof(q) for q in gr])
    p = np.polyfit(gr, cs, 2)
    b = -p[1] / (2 * p[0])
    m = np.polyval(p, b)
    f = lambda q: np.polyval(p, q) - m - 1
    return float(b), brentq(f, b - 1, b), brentq(f, b, b + 1), float(np.polyval(p, 2 * math.pi - 1) - m)


def _ganzes_chi(a0, ob, Q, spitzen, w0=-1.0, wa=0.0):
    """χ² des Ganzen bei a₀ [10⁻¹⁰ m/s²], ω_b und Q: Schallwinkel (gegen Plancks bestes Modell, 1σ nach Planck),
    DESI DR2 mit fester Skala, Pantheon+-Form, ω_b nach Planck; spitzen: ω_c = Q·ω_b gegen die Höhe der Spitzen (Planck).
    a₀ setzt die Dunkle Energie heute; w₀, w_a geben ihren Gang, w = −1 ist der Formteil allein."""
    alt = K.omega_b
    th_ref = K.schallwinkel(lambda a: K.dichte_de_planck())
    K.omega_b = ob
    try:
        rho0 = K.dichte_dunkle_energie(a0 * 1e-10, K.B1_ZELLE)
        th = K.schallwinkel(lambda x: rho0 * x**(-3 * (1 + w0 + wa)) * math.exp(-3 * wa * (1 - x)), Q)
        om, or_ = _om_or_q(Q)(a0 * 1e-10)
        return (((th - th_ref) / (K.theta_stern_rel * th_ref))**2
                + D.bao_absolut(_h0_q(Q, a0 * 1e-10), om, or_, K.schallhorizont(Q), w0, wa) + D.pantheon_form(om, or_, w0, wa)
                + ((ob - alt) / K.omega_b_fehler)**2 + (((Q * ob - K.omega_c) / K.omega_c_fehler)**2 if spitzen else 0.0))
    finally:
        K.omega_b = alt


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
    from scipy.optimize import brentq, minimize_scalar
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
def _kaltes_gas():
    """Anteil der Baryonen im kalten Gas heute (ALFALFA, xCOLD GASS) im Punkt mit Formteil 1,06 und Sättigung."""
    p = _cmb()["saett_form"]["punkt"]
    return K.anteil_kaltes_gas(p["wb"], p["H0"])


@lru_cache(maxsize=None)
def _scherung_winkel():
    """KiDS-1000 nach Winkeln im besten Punkt des Modells (mit Sättigung): je Teil (ξ₊ bis 10′, 10′ bis 60′, ab 60′, ξ₋ bis
    30′, ab 30′) χ², A, S8, Zahl der Werte und Messung gegen Faltung R ± Fehler (linear geschätzt)."""
    Q = K.q_je_radiant()
    p = D.planck_anpassung(Q, ("spektren", "linse", "bao", "sn", "wachstum", "bbn"), zustand="saettigung",
                           start={"wa": -0.56})["punkt"]
    pt = (p["wb"], p["wc"], p["H0"], p["lnAs"], p["ns"], p.get("w0", -1.0), p.get("wa", 0.0))
    return {n: D.kids_scherung_chi2(pt, verhaeltnis=True, winkel=w, art=a)
            for n, w, a in (("p_klein", (0, 10), "+"), ("p_mitte", (10, 60), "+"), ("p_gross", (60, 1e9), "+"),
                            ("m_klein", (0, 30), "-"), ("m_gross", (30, 1e9), "-"))}


@lru_cache(maxsize=None)
def _haufen_staub():
    """Konzentration der Mode in der Form nach dem Einfall (NFW), die X-COP und HIFLUGCS lesen."""
    from scipy.optimize import minimize_scalar
    return tuple(minimize_scalar(lambda c: f(c)[0], bounds=(1, 20), method="bounded").x
                 for f in (D.xcop_staub, D.hiflugcs_staub))


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
def _ganzes_q(spitzen):
    """Q aus dem Ganzen (Schallwinkel, Schallwelle, Abstände; mit spitzen auch ω_c): Bestwert, 1σ aus der Krümmung,
    Δχ² von 2π − 1 und vom gemessenen ω_c/ω_b."""
    from scipy.optimize import minimize_scalar
    f = lambda q: _ganzes_bestes(round(q, 6), spitzen)[0]
    b = minimize_scalar(f, bounds=(5.0, 5.6), method="bounded", options=dict(xatol=1e-4)).x
    d = 0.05
    k = (f(b + d) - 2 * f(b) + f(b - d)) / d**2
    return float(b), float(math.sqrt(2 / k)), f(2 * math.pi - 1) - f(b), f(K.omega_c / K.omega_b) - f(b)


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


@lru_cache(maxsize=None)
def _lebensdauer_spanne():
    """lg T_α [s] im gesicherten Bereich, kleinster und größter Wert beider Formeln."""
    lg = [KE.lg_t_alpha(Z, KE.q_alpha(Z, round(KE.a_tal(Z, KE.TROPFEN)), f)) for Z in range(119, KE.grenzen()["z_sicher"] + 1)
          for f in (KE.TROPFEN, KE.SCHALE)]
    return min(lg), max(lg)


def _sigma(dchi2, frei):
    from scipy.stats import chi2, norm
    return norm.isf(chi2.sf(dchi2, frei) / 2)


MODELLZAHLEN = [
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
    ("03", "(KiDS-1000, 60 Werte bis 2,6 Mpc)", float(sum(len(b["R"]) for b in D.kids_isoliert())), 60),
    ("03", "(KiDS-1000, 60 Werte bis 2,6 Mpc)", max(b["R"][-1] for b in D.kids_isoliert()) / 1e3, 2.6),
    ("03", "ohne äußeres Feld, χ² = 503", _kids()[0], 503),
    ("03", "umgerechnet wie die Messung: 495)", D.kids_chi2(umgerechnet=True)[0], 495),
    ("03", "10⁻⁴·a₀ im Betrag kostet Δχ² = 783", _kids()[1] - _kids()[0], 783),
    ("03", "bei 35 kpc 0,25 dex darunter", -D.kids_gang()[1][0], 0.25),
    ("03", "bei 35 kpc 0,25 dex darunter", float(D.kids_gang()[0][0]), 35),
    ("03", "bei 561 kpc 0,15 darüber", float(max(D.kids_gang()[1])), 0.15),
    ("03", "bei 561 kpc 0,15 darüber", float(D.kids_gang()[0][list(D.kids_gang()[1]).index(max(D.kids_gang()[1]))]), 561),
    ("03", "bei 2,6 Mpc 0,11 darunter", -D.kids_gang()[1][-1], 0.11),
    ("03", "(15 Bins in g_bar, Brouwer 2021): ohne Feld χ² = 21", D.kids_typ_reste("blau")[3], 21),
    ("03", "die blauen Linsen +0,01 ± 0,02 dex auf der Regel", D.kids_typ_reste("blau")[1], 0.01),
    ("03", "Sechs gasreiche (Mancera Piña", float(D.udg_reste()[0]), 6),
    ("03", "äußerer Punkt bei a_N = 0,018 bis 0,039·a₀", D.udg_reste()[3], 0.018),
    ("03", "äußerer Punkt bei a_N = 0,018 bis 0,039·a₀", D.udg_reste()[4], 0.039),
    ("03", "drehen wie Newton allein, im Mittel +0,01 dex", D.udg_reste()[2], 0.01),
    ("03", "Die Regel verlangt 0,40 dex mehr Geschwindigkeit", -D.udg_reste()[1], 0.40),
    ("03", "einen 2,5-mal kleineren Sinus der Neigung", 10**(-D.udg_reste()[1]), 2.5),
    ("03", "HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5", D.hiflugcs_gang()[0], 9.8),
    ("03", "HIFLUGCS liest g/g_bar innerhalb 300 kpc 9,8, außerhalb 1 Mpc 6,5", D.hiflugcs_gang()[1], 6.5),
    ("03", "X-COP liest c = 4,2", _haufen_staub()[0], 4.2),
    ("03", "An KiDS-1000 (ξ₊ und ξ₋, 225 Werte", float(_scherung()["modell"][0][3]), 225),
    ("03", "liegen Modell und ΛCDM gleich, χ² = 307,8", _scherung()["modell"][0][0], 307.8),
    ("03", "liegen Modell und ΛCDM gleich, χ² = 307,8", _scherung()["lcdm"][0][0], 307.8),
    ("03", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][1][0], 0.73),
    ("03", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][0][2], 0.82),
    ("03", "frei S8 = 0,73 statt 0,82 (Δχ² = 38)", _scherung()["modell"][0][0] - _scherung()["modell"][1][1], 38),
    ("03", "ξ₊ unter 10′ (60 Werte) bei 0,76 ± 0,05", _scherung_winkel()["p_klein"][4], 0.76),
    ("03", "ξ₊ unter 10′ (60 Werte) bei 0,76 ± 0,05", _scherung_winkel()["p_klein"][5], 0.05),
    ("03", "von 10′ bis 60′ bei 0,92 ± 0,07, ab 60′ bei 0,94 ± 0,21", _scherung_winkel()["p_mitte"][4], 0.92),
    ("03", "von 10′ bis 60′ bei 0,92 ± 0,07, ab 60′ bei 0,94 ± 0,21", _scherung_winkel()["p_mitte"][5], 0.07),
    ("03", "von 10′ bis 60′ bei 0,92 ± 0,07, ab 60′ bei 0,94 ± 0,21", _scherung_winkel()["p_gross"][4], 0.94),
    ("03", "von 10′ bis 60′ bei 0,92 ± 0,07, ab 60′ bei 0,94 ± 0,21", _scherung_winkel()["p_gross"][5], 0.21),
    ("03", "liegt unter 30′ bei 0,79 ± 0,06", _scherung_winkel()["m_klein"][4], 0.79),
    ("03", "liegt unter 30′ bei 0,79 ± 0,06", _scherung_winkel()["m_klein"][5], 0.06),
    ("03", "bei 0,81 ± 0,06. Einen Gang", _scherung_winkel()["m_gross"][4], 0.81),
    ("03", "bei 0,81 ± 0,06. Einen Gang", _scherung_winkel()["m_gross"][5], 0.06),
    ("03", "ξ₊ unter 10′ (60 Werte) bei 0,76 ± 0,05", float(_scherung_winkel()["p_klein"][3]), 60),
    ("03", "liest HIFLUGCS Q = 5,5 ± 0,3", D.haufen_q()[0], 5.5),
    ("03", "liest HIFLUGCS Q = 5,5 ± 0,3", D.haufen_q_fehler(), 0.3),
    ("00", "die Haufen unabhängig 5,5 ± 0,3", D.haufen_q()[0], 5.5),
    ("03", "(Rest im Kern −0,04 dex, in der Form des Gases +0,41)", D.xcop_staub(_haufen_staub()[0])[2], -0.04),
    ("03", "(Rest im Kern −0,04 dex, in der Form des Gases +0,41)", D.xcop_staub()[2], 0.41),
    ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", _haufen_staub()[1], 5.5),
    ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", D.hiflugcs_staub(_haufen_staub()[1])[0], 447),
    ("03", "HIFLUGCS c = 5,5 (χ² = 447 statt 601", D.hiflugcs_staub()[0], 601),
    ("03", "innerhalb 300 kpc −0,04 statt +0,21 dex", D.hiflugcs_staub(_haufen_staub()[1])[2], -0.04),
    ("03", "innerhalb 300 kpc −0,04 statt +0,21 dex", D.hiflugcs_staub()[2], 0.21),
    ("03", "im Kern das 6,8-fache des Gases tragen", D.xcop_sterne()[2][0], 6.8),
    ("F:ganzes", "bei den roten Linsen innen das 0,6-fache, außen das 12-fache", D.kids_heiss("rot", -12.5)[1], 0.6),
    ("F:ganzes", "bei den roten Linsen innen das 0,6-fache, außen das 12-fache", D.kids_heiss("rot", -14.5, -13)[1], 12),
    ("F:ganzes", "bei den frühen 0,3 und 11", D.kids_heiss("frueh", -12.5)[1], 0.3),
    ("F:ganzes", "bei den frühen 0,3 und 11", D.kids_heiss("frueh", -14.5, -13)[1], 11),
    ("F:ganzes", "ist χ² = 167 gegen 315", D.kids_typ_reste("blau")[3] + D.kids_typ_lcdm("rot")[2], 167),
    ("F:ganzes", "ist χ² = 167 gegen 315", D.kids_typ_lcdm("blau")[2] + D.kids_typ_lcdm("rot")[2], 315),
    ("F:ganzes", "und 408 für die Regel bei beiden", D.kids_typ_reste("blau")[3] + D.kids_typ_reste("rot")[3], 408),
    ("F:ganzes", "nach dem Sérsic-Index 100 gegen 243 und 282", D.kids_typ_reste("spaet")[3] + D.kids_typ_lcdm("frueh")[2], 100),
    ("F:ganzes", "nach dem Sérsic-Index 100 gegen 243 und 282", D.kids_typ_lcdm("spaet")[2] + D.kids_typ_lcdm("frueh")[2], 243),
    ("F:ganzes", "nach dem Sérsic-Index 100 gegen 243 und 282", D.kids_typ_reste("spaet")[3] + D.kids_typ_reste("frueh")[3], 282),
    ("03", "das Maß des Ganzen 1 + Q, das dort dasselbe gibt, 0,39 dex", 0.5 * math.log10(1 + K.q_je_radiant()) - D.udg_reste()[2], 0.39),
    ("03", "Die sieben gasreichen isolierten Zwerge an derselben Stelle", float(len(D.zwerge_gruppe("gasreich"))), 7),
    ("03", "(a_N im Median 0,026·a₀)", statistics.median(D.zwerg_feld(z) / K.a0 for z in D.zwerge_gruppe("gasreich")), 0.026),
    ("03", "die blauen Linsen +0,01 ± 0,02 dex auf der Regel", D.kids_typ_reste("blau")[2], 0.02),
    ("03", "die roten innen +0,13 ± 0,02 darüber", D.kids_typ_reste("rot", -12.5)[1], 0.13),
    ("03", "die roten innen +0,13 ± 0,02 darüber", D.kids_typ_reste("rot", -12.5)[2], 0.02),
    ("03", "(15 Bins in g_bar, Brouwer 2021)", float(D.kids_typ_reste("blau")[0]), 15),
    ("03", "mit 0,003·a₀ 534", D.kids_typ_reste("blau", g_e=0.003 * K.a0)[3], 534),
    ("03", "mit 0,01·a₀ 996", D.kids_typ_reste("blau", g_e=0.01 * K.a0)[3], 996),
    ("03", "χ² = 27483 statt 503", _kids()[2], 27483),
    ("03", "Naht (30–66 kpc)", K.naht_radius(D.kids_isoliert()[0]["M_b"] * K.M_sun, K.q_je_radiant()) / K.kpc, 30),
    ("03", "Naht (30–66 kpc)", K.naht_radius(D.kids_isoliert()[3]["M_b"] * K.M_sun, K.q_je_radiant()) / K.kpc, 66),
    ("03", "isolierten +0,04 ± 0,04 dex", _zwerge("isoliert", False)[1], 0.04),
    ("03", "isolierten +0,04 ± 0,04 dex", _zwerge("isoliert", False)[2], 0.04),
    ("03", "Die 22 klassischen Satelliten", float(_zwerge("außen", False)[0]), 22),
    ("03", "+0,15 ± 0,04 dex ohne Wirtsfeld", _zwerge("außen", False)[1], 0.15),
    ("03", "+0,15 ± 0,04 dex ohne Wirtsfeld", _zwerge("außen", False)[2], 0.04),
    ("03", "+0,46 ± 0,05 mit ihm", _zwerge("außen", True)[1], 0.46),
    ("03", "+0,46 ± 0,05 mit ihm", _zwerge("außen", True)[2], 0.05),
    ("03", "(χ² = 76 gegen 341)", _zwerge("außen", False)[3], 76),
    ("03", "(χ² = 76 gegen 341)", _zwerge("außen", True)[3], 341),
    ("03", "gasreichen isolierten 0,00 ± 0,04 dex", _zwerge("gasreich", False)[1], 0.00),
    ("03", "gasreichen isolierten 0,00 ± 0,04 dex", _zwerge("gasreich", False)[2], 0.04),
    ("03", "die gasarmen isolierten +0,13 ± 0,07", _zwerge("gasarm", False)[1], 0.13),
    ("03", "die gasarmen isolierten +0,13 ± 0,07", _zwerge("gasarm", False)[2], 0.07),
    ("03", "die 36 klassischen Satelliten +0,11 ± 0,03", float(_zwerge("klassisch", False)[0]), 36),
    ("03", "die 36 klassischen Satelliten +0,11 ± 0,03", _zwerge("klassisch", False)[1], 0.11),
    ("03", "die 36 klassischen Satelliten +0,11 ± 0,03", _zwerge("klassisch", False)[2], 0.03),
    ("03", "die ultraschwachen +0,37 dex", _zwerge("ultraschwach", False)[1], 0.37),
    ("03", "liest Υ_V = 1,96 für die gasreichen", D.zwerge_yv("gasreich"), 1.96),
    ("03", "5,2 für die klassischen Satelliten", D.zwerge_yv("klassisch"), 5.2),
    ("03", "weder das Verhältnis der Dichten zum Wirt", D.zwerge_gezeit(), 0.0),
    ("03", "noch das Perizentrum der Bahn", D.satelliten_perizentrum()[0], 0.0),
    ("03", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", D.zwerge_gas_beziehung()[0], 2.09),
    ("03", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", -D.zwerge_gas_beziehung()[1], 0.284),
    ("03", "lg(M_HI/L_V) = 2,09 − 0,284·lg L_V (188 Zwerge)", float(D.zwerge_gas_beziehung()[2]), 188),
    ("03", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", _zwerge_g("klassisch")[1], 0.011),
    ("03", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", _zwerge_g("klassisch")[2], 0.028),
    ("03", "liegen dann +0,011 ± 0,028 dex auf der Regel und lesen Υ_V = 2,49", D.zwerge_yv("klassisch", True), 2.49),
    ("03", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", float(D.zwerge_ausbeute()[2]), 7),
    ("03", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", D.zwerge_ausbeute()[0], -1.03),
    ("03", "die Ausbeute an den 7 gasreichen geeicht (lg p = −1,03, Streuung 0,20)", D.zwerge_ausbeute()[1], 0.20),
    ("03", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", float(_zwerge_m("klassisch")[0]), 27),
    ("03", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", _zwerge_m("klassisch")[1], -0.03),
    ("03", "Die 27 klassischen Satelliten mit gemessener Metallizität liegen dann −0,03 ± 0,03", _zwerge_m("klassisch")[2], 0.03),
    ("03", "die 27 ultraschwachen +0,07 ± 0,04", float(_zwerge_m("ultraschwach")[0]), 27),
    ("03", "die 27 ultraschwachen +0,07 ± 0,04", _zwerge_m("ultraschwach")[1], 0.07),
    ("03", "die 27 ultraschwachen +0,07 ± 0,04", _zwerge_m("ultraschwach")[2], 0.04),
    ("03", "die gasarmen isolierten +0,05 ± 0,08", _zwerge_m("gasarm")[1], 0.05),
    ("03", "die gasarmen isolierten +0,05 ± 0,08", _zwerge_m("gasarm")[2], 0.08),
    ("03", "ebenso −0,03 und +0,06", _zwerge_m("klassisch", True)[1], -0.03),
    ("03", "ebenso −0,03 und +0,06", _zwerge_m("ultraschwach", True)[1], 0.06),
    ("03", "(36 Satelliten der Milchstraße mit Eigenbewegung)", float(D.satelliten_perizentrum()[1]), 36),
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
    ("F:ganzes", "ν_G(−ln q) − 1 = 16,8 bei Q = 2π − 1", K.zustandsteil_je_gas(K.q_je_radiant()), 16.8),
    ("F:ganzes", "17,2 mit dem gemessenen Q", K.zustandsteil_je_gas(), 17.2),
    ("F:ganzes", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[1], 0.01),
    ("F:ganzes", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[2], 0.02),
    ("F:ganzes", "die blauen +0,01 ± 0,02 dex auf der Regel (χ² 21)", D.kids_typ_reste("blau")[3], 21),
    ("F:ganzes", "Über alle 15 Bins", float(D.kids_typ_reste("blau")[0]), 15),
    ("F:ganzes", "die späten +0,02 ± 0,02", D.kids_typ_reste("spaet")[1], 0.02),
    ("F:ganzes", "die späten +0,02 ± 0,02", D.kids_typ_reste("spaet")[2], 0.02),
    ("F:ganzes", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[1], 0.23),
    ("F:ganzes", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[2], 0.01),
    ("F:ganzes", "die roten +0,23 ± 0,01 darüber (χ² 387)", D.kids_typ_reste("rot")[3], 387),
    ("F:ganzes", "die frühen +0,21 ± 0,01", D.kids_typ_reste("frueh")[1], 0.21),
    ("F:ganzes", "die frühen +0,21 ± 0,01", D.kids_typ_reste("frueh")[2], 0.01),
    ("F:ganzes", "liegen die roten +0,13 ± 0,02", D.kids_typ_reste("rot", -12.5)[1], 0.13),
    ("F:ganzes", "liegen die roten +0,13 ± 0,02", D.kids_typ_reste("rot", -12.5)[2], 0.02),
    ("F:ganzes", "die frühen +0,14 ± 0,02 über der Regel", D.kids_typ_reste("frueh", -12.5)[1], 0.14),
    ("F:ganzes", "die frühen +0,14 ± 0,02 über der Regel", D.kids_typ_reste("frueh", -12.5)[2], 0.02),
    ("F:ganzes", "wächst es auf +0,29 und +0,26", D.kids_typ_reste("rot", -14.5, -13.0)[1], 0.29),
    ("F:ganzes", "wächst es auf +0,29 und +0,26", D.kids_typ_reste("frueh", -14.5, -13.0)[1], 0.26),
    ("F:ganzes", "am Gasanteil (Rangkorrelation −0,04)", D.sparc_eigenschaften(K.kpc, K.a0)["Gasanteil"][0], -0.04),
    ("F:ganzes", "noch am Typ (−0,10)", D.sparc_eigenschaften(K.kpc, K.a0)["Typ"][0], -0.10),
    ("F:ganzes", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_lcdm("blau")[1], -0.31),
    ("F:ganzes", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_lcdm("blau")[2], 169),
    ("F:ganzes", "verfehlt umgekehrt die blauen um −0,31 dex (χ² 169, die Regel 21)", D.kids_typ_reste("blau")[3], 21),
    ("F:ganzes", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_lcdm("spaet")[1], -0.26),
    ("F:ganzes", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_lcdm("spaet")[2], 157),
    ("F:ganzes", "und die späten um −0,26 (157 gegen 14)", D.kids_typ_reste("spaet")[3], 14),
    ("F:ganzes", "trifft die roten im Mittel (χ² 146 gegen 387)", D.kids_typ_lcdm("rot")[2], 146),
    ("F:ganzes", "trifft die roten im Mittel (χ² 146 gegen 387)", D.kids_typ_reste("rot")[3], 387),
    ("F:ganzes", "und die frühen (85 gegen 268)", D.kids_typ_lcdm("frueh")[2], 85),
    ("F:ganzes", "und die frühen (85 gegen 268)", D.kids_typ_reste("frueh")[3], 268),
    ("F:ganzes", "mit χ² 34 bei 20 Werten", D.mistele_reste()[0], 34),
    ("03", "bei R₅₀₀ 3,4 bis 7,0 keV", D.xcop_abkuehlung()[2], 3.4),
    ("03", "bei R₅₀₀ 3,4 bis 7,0 keV", D.xcop_abkuehlung()[3], 7.0),
    ("F:ganzes", "im Median nur 0,4 % des Gases innerhalb R₅₀₀", D.xcop_abkuehlung()[0] * 100, 0.4),
    ("03", "höchstens 3,7 % des Gases abkühlbar", D.xcop_abkuehlung()[1] * 100, 3.7),
    ("03", "bei 0,1·R₅₀₀ 4,8 bis 9,6 keV", D.xcop_abkuehlung(ort=0.1)[2], 4.8),
    ("03", "bei 0,1·R₅₀₀ 4,8 bis 9,6 keV", D.xcop_abkuehlung(ort=0.1)[3], 9.6),
    ("03", "Abkühlzeit von 13,6 bis 38 Gyr", D.xcop_abkuehlung(ort=0.1)[4], 13.6),
    ("03", "Abkühlzeit von 13,6 bis 38 Gyr", D.xcop_abkuehlung(ort=0.1)[5], 38),
    ("F:ganzes", "in 13,8 Gyr im Median nur", 13.8, 13.8),
    ("03", "in 13,8 Gyr höchstens 3,7 %", 13.8, 13.8),
    ("F:ganzes", "die Glättungslänge ℓ = 349 kAU", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3, 349),
    ("F:ganzes", "(ℓ/2^{1/3} = 277 kAU)", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3 / 2**(1 / 3), 277),
    ("F:ganzes", "bei 30 kAU (1,5 M☉) ν = 3,8", K.verstaerkung(K.G * 1.5 * K.M_sun / (30e3 * K.AU)**2), 3.8),
    ("F:ganzes", "die 1,95-fache Geschwindigkeit", K.verstaerkung(K.G * 1.5 * K.M_sun / (30e3 * K.AU)**2)**0.5, 1.95),
    ("F:ganzes", "(Gaia DR3, bis 130 pc, 1 bis 30 kAU)", float(K.WEITE_DOPPELSTERNE["abstand_pc"]), 130),
    ("F:ganzes", "der Zuwachs von 20 % nach MOND fehlt", K.WEITE_DOPPELSTERNE["mond_zuwachs"] * 100, 20),
    ("03", "die Glättungslänge ℓ = 349 kAU", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3, 349),
    ("03", "eines Paares von 1 M☉ (277 kAU)", K.glaettungslaenge(K.M_sun, K.gezeit_flach(K.SONNE_GALAKTISCH["v"][1] * 1e3, K.SONNE_GALAKTISCH["R0"] * K.kpc)) / K.AU / 1e3 / 2**(1 / 3), 277),
    ("F:ganzes", "mit χ² 34 bei 20 Werten", float(D.mistele_reste()[1]), 20),
    ("F:ganzes", "im Mittel je Bin −0,03 bis +0,03 dex", min(D.mistele_reste()[2]), -0.03),
    ("F:ganzes", "im Mittel je Bin −0,03 bis +0,03 dex", max(D.mistele_reste()[2]), 0.03),
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
    ("03", "liest Q = 5,22", _ganzes_q(False)[0], 5.22),
    ("03", "Q = 5,22 ± 0,05", _ganzes_q(False)[1], 0.05),
    ("03", "2π − 1 liegt Δχ² = 1,7 daneben", _ganzes_q(False)[2], 1.7),
    ("03", "Verhältnis Δχ² = 8,7", _ganzes_q(False)[3], 8.7),
    ("03", "ist Q = 5,26", _ganzes_q(True)[0], 5.26),
    ("03", "Q = 5,26 ± 0,04", _ganzes_q(True)[1], 0.04),
    ("03", "2π − 1 liegt Δχ² = 0,2 daneben", _ganzes_q(True)[2], 0.2),
    ("03", "das gemessene Verhältnis 5,0", _ganzes_q(True)[3], 5.0),
    ("03", "zusammen lesen a₀ = 1,088", _ganzes_a0(False)[0], 1.088),
    ("03", "(1,084–1,092)", _ganzes_a0(False)[1], 1.084),
    ("03", "1,084–1,092)", _ganzes_a0(False)[2], 1.092),
    ("03", "mit den Spitzen 1,085", _ganzes_a0(True)[0], 1.085),
    ("03", "(1,082–1,089)", _ganzes_a0(True)[1], 1.082),
    ("03", "1,082–1,089)", _ganzes_a0(True)[2], 1.089),
    ("03", "bei H₀ = 68,2", _ganzes_a0(True)[3], 68.2),
    ("03", "und Ω_m = 0,30", _ganzes_a0(True)[4], 0.30),
    ("03", "Abstände von 79 der 141 Galaxien", _sparc_ganzes()[5], 79),
    ("03", "79 der 141 Galaxien", _sparc_ganzes()[6], 141),
    ("03", "Hubble-Fluss mit H₀ = 73", D.H0_SPARC, 73),
    ("03", "H₀ = 73: a₀ = 1,085", _sparc_ganzes()[0], 1.085),
    ("03", "(0,98–1,19)", _sparc_ganzes()[1], 0.98),
    ("03", "0,98–1,19), bei", _sparc_ganzes()[2], 1.19),
    ("03", "bei Υ_d = 0,66", _sparc_ganzes()[3], 0.66),
    ("03", "w₀ = −0,85", _ganzes_zustand(2 * math.pi - 1)[1], -0.85),
    ("03", "w_a = −0,49", _ganzes_zustand(2 * math.pi - 1)[2], -0.49),
    ("03", "Δχ² = 8,5 gegen den Formteil allein", _ganzes_zustand(2 * math.pi - 1)[0], 8.5),
    ("03", "gemessenen Q Δχ² = 12,7", _ganzes_zustand(K.omega_c / K.omega_b)[0], 12.7),
    ("03", "liegt dann bei z = 0,43", _ganzes_maximum()[0], 0.43),
    ("03", "(0,34–0,55)", _ganzes_maximum()[1], 0.34),
    ("03", "0,34–0,55), zwischen", _ganzes_maximum()[2], 0.55),
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
    ("02", "850 ungerade Kerne mit eindeutig gemessenem Spin", KE.grundspin_pruefung(0.12, mu=0.03)[2], 850),
    ("02", "321 Spins mit Parität (38 %), 641 Paritäten (75 %)", KE.grundspin_pruefung(0.12, mu=0.03)[0] * 1000 + KE.grundspin_pruefung(0.12, mu=0.03)[1], 321641),
    ("02", "606 Kerne) 299 Spins und 520 Paritäten (49 %, 86 %)", sum(x * y for x, y in zip(KE.grundspin_pruefung(0.12, mu=0.03, kugel=True), (10**6, 1000, 1))), 299520606),
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
    ("03", "(16–84 %: 0,37–1,28)", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[1][0], 0.37),
    ("03", "0,37–1,28) über vier", D.sparc_scheiben(K.kpc, D.sparc_eichung(K.kpc), K.M_sun, K.a0, K.G)[1][1], 1.28),
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
    ("03", "(H₀ = 68,9) sind es 0,14 mag", 5 * math.log10(D.pantheon_h0(*_om_or(), 0.023, 0.15)[0] / K.atmen_heute(_bao_a0()[0] * 1e-10)[0]), 0.14),
    ("03", "= 0,20 mag heller", 5 * math.log10(D.pantheon_h0(*_om_or(), 0.023, 0.15)[0] / K.atmen_heute()[0]), 0.20),
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
    ("00", "1,3 σ neben Plancks Wert", abs(K.q_je_radiant() - K.q_messung()[0]) / K.q_messung()[1], 1.3),
    ("03", "Pantheon+, 1590 Supernovae", float(len(D.pantheon()[0])), 1590),
    ("03", "nur Δχ² = 2,0 über", _abst()["modell"] - _abst()["lcdm"], 2.0),
    ("03", "ΛCDM (Ω_m = 0,304)", _abst()["om_lcdm"], 0.304),
    ("03", "Ω_m = 0,304), 1,4 σ", _sigma(_abst()["modell"] - _abst()["lcdm"], 1), 1.4),
    ("03", "(w₀ = −0,89,", _abst()["w"][1], -0.89),
    ("03", "w_a = −0,21)", _abst()["w"][2], -0.21),
    ("03", "Δχ² = 6,9 unter", _abst()["modell"] - _abst()["w0wa"], 6.9),
    ("03", "dem Formteil, 1,8 σ", _sigma(_abst()["modell"] - _abst()["w0wa"], 3), 1.8),
    ("03", "r_d = 149,0 Mpc", _abst()["h0rd"] / K.atmen_heute()[0], 149.0),
    ("03", "Der Formteil trifft es mit χ² = 6,8", B.kosmos()["chi2"], 6.8),
    ("08", "(10⁶ M☉: s² = 1,4·10¹⁷", K.c**4 / (4 * K.G * 1e6 * K.M_sun) / K.a0 / 1e17, 1.4),
    ("08", "am Horizont, r* = 36 pc)", K.schwellenradius(1e6 * K.M_sun) / K.pc, 36),
    ("08", "= 1 − ⅔√2 = 5,7 %", 100 * (1 - 2 / 3 * math.sqrt(2)), 5.7),
    ("08", "sechsmal so viel wie Kernbindung bis zum Eisen", (1 - 2 / 3 * math.sqrt(2)) / ([e for e in D.ame2020() if e[0] == 26 and e[2] == 56][0][3] / (K.m_u * K.c**2 / K.e / 1e3)), 6),
    ("08", "(j = 1) sind es 1 − 1/√3 = 42 %", 100 * (1 - 1 / math.sqrt(3)), 42),
    ("08", "ist T/T_I = M_U/(8πM)", K.falte_temperatur(K.M_sun) * K.c * K.k_B / (K.hbar * K.a0) / (K.schnittmassen()[3] / (8 * math.pi * K.M_sun)), 1.0),
    ("08", "r₀ = Λ/2π, dieselbe Zählung", 2 * K.G * K.schnittmassen()[3] / (8 * math.pi) / K.c**2 / (K.schliesslaenge() / (2 * math.pi)), 1.0),
    ("08", "eine 2,3·10²² mal kürzere Periode", K.falte_temperatur(K.M_sun) * K.c * K.k_B / (K.hbar * K.a0) / 1e22, 2.3),
    ("08", "Erst bei M_U/4 = 1,4·10²³ M☉", K.schnittmassen()[3] / 4 / K.M_sun / 1e23, 1.4),
    ("08", "Am Horizont ist s² = M_U/(4M)", K.c**4 / (4 * K.G * K.M_sun) / K.a0 / (K.schnittmassen()[3] / (4 * K.M_sun)), 1.0),
    ("08", "kT = m̃c²/(4π). Gegen", K.falte_temperatur(K.M_sun) * 24 * math.pi * K.k_B * K.faltung(K.M_sun) / (K.hbar * K.c), 1.0),
    ("03", "(r_d = 147,1 Mpc aus ω_b und ω_c)", K.schallhorizont(), 147.1),
    ("03", "liest DESI dann a₀ = 1,10", _bao_a0()[0], 1.10),
    ("03", "a₀ = 1,10·10⁻¹⁰ m/s² (1,09–", _bao_a0()[1], 1.09),
    ("03", "m/s² (1,09–1,11", _bao_a0()[2], 1.11),
    ("03", "mit H₀ = 68,9 bei χ² = 10,4", K.atmen_heute(_bao_a0()[0] * 1e-10)[0], 68.9),
    ("03", "bei χ² = 10,4 für 13 Werte", _bao_a0()[3], 10.4),
    ("03", "für 13 Werte", float(len(D.desi()[0])), 13),
    ("03", "Der Standardwert 1,06 liegt Δχ² = 27", D.bao_absolut(K.atmen_heute()[0], *_om_or(), K.schallhorizont()) - _bao_a0()[3], 27),
    ("03", "Mit Q = 2π − 1 statt des gemessenen ebenso 1,10", _bao_a0(2 * math.pi - 1)[0], 1.10),
    ("03", "(r_d = 147,6 Mpc, H₀ = 68,9", K.schallhorizont(2 * math.pi - 1), 147.6),
    ("03", "r_d = 147,6 Mpc, H₀ = 68,9", _h0_q(2 * math.pi - 1, _bao_a0(2 * math.pi - 1)[0] * 1e-10), 68.9),
    ("03", "H₀ = 68,9, χ² = 10,3)", _bao_a0(2 * math.pi - 1)[3], 10.3),
    ("03", "ist Q = 5,5 (5,3–5,7)", _abst_q()[0], 5.5),
    ("03", "Q = 5,5 (5,3–", _abst_q()[1], 5.3),
    ("03", "(5,3–5,7). 2π − 1", _abst_q()[2], 5.7),
    ("03", "2π − 1 liegt Δχ² = 0,9 daneben", _abst_q()[3], 0.9),
    ("03", "an sechs Rotverschiebungen", float(len(B.kosmos()["zeilen"])), 6),
    ("01", "(PubChem, 13,598 eV)", B.atom_und_kern()["IE_H"], 13.598),
    ("03", "allein (ω_b, ω_c nach Planck), ist a₀ = 1,09", D.abstaende_a0(_om_or_a0)[0], 1.09),
    ("03", "(1,07–1,11 bei Δχ² = 1)", D.abstaende_a0(_om_or_a0)[1], 1.07),
    ("03", "1,07–1,11 bei Δχ²", D.abstaende_a0(_om_or_a0)[2], 1.11),
    ("03", "mit H₀ = 68,4: Rotationskurven", K.atmen_heute(D.abstaende_a0(_om_or_a0)[0] * 1e-10)[0], 68.4),
    ("03", "lesen die Abstände a₀ = 1,08", D.abstaende_a0(_om_or_q(2 * math.pi - 1))[0], 1.08),
    ("03", "1,08 (1,06–1,10)", D.abstaende_a0(_om_or_q(2 * math.pi - 1))[1], 1.06),
    ("03", "(1,06–1,10), H₀", D.abstaende_a0(_om_or_q(2 * math.pi - 1))[2], 1.10),
    ("03", "H₀ = 68,0, bei gleichem", _h0_q(2 * math.pi - 1, D.abstaende_a0(_om_or_q(2 * math.pi - 1))[0] * 1e-10), 68.0),
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
    ("02", "auf 0,32 Protonen genau", KE.tal_vergleich(KE.TROPFEN)[0], 0.32),
    ("02", "219 Massenzahlen bis A = 238", float(KE.tal_vergleich(KE.TROPFEN)[2]), 219),
    ("02", "erst bei Z = 134 (A = 370)", KE.kippstelle(KE.TROPFEN)[0], 134),
    ("02", "Z = 134 (A = 370)", KE.kippstelle(KE.TROPFEN)[1], 370),
    ("02", "erst bei Z = 134", float(KE.kippstelle_ganz(KE.TROPFEN)[0]), 134),
    ("02", "Schließzahlen des Modells bei Z = 181 bis 191", KE.kippstelle(KE.SCHALE)[0], 181),
    ("02", "bei Z = 181 bis 191", KE.kippstelle(("schale", 114))[0], 191),
    ("02", "sinkt von 3,39 auf 1,72 MeV", KE.anpassung(KE.TROPFEN)[1], 3.39),
    ("02", "auf 1,72 MeV, das leichteste", KE.anpassung(KE.SCHALE)[1], 1.72),
    ("02", "trifft sie auf 0,26 statt 0,32", KE.tal_vergleich(KE.SCHALE)[0], 0.26),
    ("02", "Massen nur bis Z = 110", float(KE.bekannte()["z_gem_max"]), 110),
    ("02", "Die Hülle hält bis Z = 170 bis 171", min(KE.grenzen()["huelle"]), 170),
    ("02", "Hülle hält bis Z = 170 bis 171:", max(KE.grenzen()["huelle"]), 171),
    ("02", "sind die Elemente bis Z = 133", float(KE.grenzen()["z_sicher"]), 133),
    ("02", "mit N ≈ 202–235", float(round(KE.n_tal(119, KE.TROPFEN)) == 202 and round(KE.n_tal(KE.grenzen()["z_sicher"], KE.TROPFEN)) == 235), 1.0),
    ("02", "reichen sie bis zur Grenze der Hülle, Z = 170", float(KE.grenzen()["z_max"]), 170),
    ("02", "(x = 1,06 ohne", KE.x(126, 310, KE.TROPFEN), 1.06),
    ("02", "0,84 mit Schalenglied)", KE.x(126, 310, KE.SCHALE), 0.84),
    ("02", "gegen den Alpha-Zerfall Sekunden bis Stunden", float(_lebensdauer_spanne()[0] >= 0 and 3600 <= 10**_lebensdauer_spanne()[1] < 86400), 1.0),
    ("02", "(36 Paare, AME2020)", float(D.spiegelkerne()[0]), 36),
    ("02", "Zα = 1 (Z = 137) gilt für einen punktförmigen Kern", math.floor(1 / K.alpha), 137),
    ("02", "An 82 gerade-geraden Alpha-Strahlern", float(KE.alpha_lebensdauer()[2]), 82),
    ("02", "trifft das auf 0,51 Dekaden", KE.alpha_lebensdauer()[1], 0.51),
    ("02", "Mit dem Edelgas Og schließt die siebte Periode bei Z = 118 (PubChem)", float(KE.konfigurationen()[0] == (2, 10, 18, 36, 54, 86, 118) == K.schliesszahlen_huelle(7)), 1.0),
    ("02", "Die Perioden sind 2, 8, 8, 18, 18, 32, 32 lang", float(all(b - a == 2 * math.ceil((p + 1) / 2)**2 for p, (a, b) in enumerate(zip((0,) + KE.konfigurationen()[0], KE.konfigurationen()[0]), 1))), 1.0),
    ("02", "99 der 118 Elemente füllen genau nach n + l", float(KE.konfigurationen()[2] - len(KE.konfigurationen()[1])), 99),
    ("02", "Die 19 Abweichungen", float(len(KE.konfigurationen()[1])), 19),
    ("02", "Acht schließen dabei eine halbe oder volle Schale", float(sum(x[1] == "schliesst" for x in KE.konfigurationen()[1])), 8),
    ("02", "sieben beginnen die f-Schale mit einem d-Elektron", float(sum(x[1] == "d_vor_f" for x in KE.konfigurationen()[1])), 7),
    ("02", "vier geben ein s-Elektron an die d-Schale", float(sum(x[1] == "s_an_d" for x in KE.konfigurationen()[1])), 4),
    ("02", "mit 8s (Z = 119, 120) und der ersten g-Schale", float([(x[2], x[3]) for x in K.fuellung(8) if x[:2] in ((8, 0), (5, 4))] == [(119, 120), (121, 138)]), 1.0),
    ("02", "und schlösse bei 168", K.schliesszahlen_huelle(8)[-1], 168),
    ("02", "bindet an der Schwelle E = 0 genau bei den Kopplungen N(N + 1)", float(max(abs(v - (n + l) * (n + l + 1)) for (n, l), v in __import__("huelle_4_formen").schliessform().items()) < 0.01), 1.0),
    ("02", "Sein Feld folgt der Spiegelform innen auf 0,6 %", 100 * _ght()["innen"], 0.6),
    ("02", "Seine Höhe v = 3,44·Z^{2/3} liegt von Ne bis Og im Fenster", float(round(_ght()["k"], 2) == 3.44 and all(e["N"] == e["p"] + 1 for Z, e in _ght()["edel"].items() if Z >= 10)), 1.0),
    ("02", "(104,8 zwischen 90 und 110)", float(round(_ght()["v168"], 1) == 104.8 and 90 < _ght()["v168"] < 110), 1.0),
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
    ("03", "(16–84 %: 1,00–1,21)", D.sparc_schwelle(K.kpc, yd=D.sparc_eichung(K.kpc))[1][0] / 1e-10, 1.00),
    ("03", "1,00–1,21). Die Bose", D.sparc_schwelle(K.kpc, yd=D.sparc_eichung(K.kpc))[1][1] / 1e-10, 1.21),
    ("03", "die einfache (Δχ² = 53)", D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["einfach"]
        - D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["Bose"], 53),
    ("03", "Standard-Form (Δχ² = 41)", D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["standard"]
        - D.sparc_formen(K.kpc, D.sparc_eichung(K.kpc))[0]["Bose"], 41),
    ("03", "Rückkehr ist 0,92 ± 0,05", D.sparc_dimension(K.kpc, yd=D.sparc_eichung(K.kpc))[0], 0.92),
    ("03", "0,92 ± 0,05. Der Rest", D.sparc_dimension(K.kpc, yd=D.sparc_eichung(K.kpc))[1], 0.05),
    ("07", "Dimension 0,92 ± 0,05", D.sparc_dimension(K.kpc, yd=D.sparc_eichung(K.kpc))[0], 0.92),
    ("03", "nach Planck a₀ = 1,065·10⁻¹⁰", K.a0_aus_atmen() / 1e-10, 1.065),
    ("03", "m/s² H₀ = 67,2 km/s/Mpc", K.atmen_heute(1.06e-10)[0], 67.2),
    ("03", "und Ω_DE = 0,685, mit", K.atmen_heute(1.06e-10)[1], 0.685),
    ("03", "1,1·10⁻¹⁰ H₀ = 68,9", K.atmen_heute(1.1e-10)[0], 68.9),
    ("03", "68,9 und 0,700", K.atmen_heute(1.1e-10)[1], 0.700),
    ("03", "erst bei z = −0,10 ein", ((K.vakuumteil() / (K.Omega_b + K.Omega_c))**(1 / 3)) - 1, -0.1),
    ("03", "= 0,47·v_flach", K.haltegrenze(1.0), 0.47),
    ("03", "unter v_flach = 17–21", 8e3 / K.haltegrenze(1e3), 17),
    ("03", "17–21 km/s. In SPARC", 10e3 / K.haltegrenze(1e3), 21),
    ("03", "bei z_eq = 3402", K.z_gleichheit(), 3402),
    ("03", "bei 3–4·10¹²", K.grenze_verschmelzen("EH") / 1e12, 3),
    ("03", "3–4·10¹² M☉ (mit", K.grenze_verschmelzen("BBKS") / 1e12, 4),
    ("03", "Bereiche ab 1,4·10¹¹", K.bindungsmasse(100e3) / K.M_sun / 1e11, 1.4),
    ("03", "bis 1,2·10¹³ M☉", K.bindungsmasse(300e3) / K.M_sun / 1e13, 1.2),
    ("03", "M_Gas = 0,6–1,4·10⁷ M☉", K.bindungsmasse(8e3) / K.M_sun / 1e7, 0.6),
    ("03", "1,4·10⁷ M☉. Kaltes", K.bindungsmasse(10e3) / K.M_sun / 1e7, 1.4),
    ("03", "zum Endwert 1,23", 2 * math.pi * K.a0 / (K.c * math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie())), 1.23),
    ("03", "drei Teile tun es bei z = 0,30", (K.dichte_dunkle_energie() / K.dichte_heute()[2])**(1 / 3) - 1, 0.3),
    ("03", "Mit drei Teilen in 55,6 km/s/Mpc", math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) * K.Mpc / 1e3, 55.6),
    ("03", "9,1 km/s bei 10⁷ M☉", K.streuung_kugel(1e7 * K.M_sun) / 1e3, 9.1),
    ("03", "19,3 km/s bei 2·10⁸ M☉", K.streuung_kugel(2e8 * K.M_sun) / 1e3, 19.3),
    ("03", "10¹² M☉/Mpc: 2,6 kpc", K.filament_schwelle(1e12 * K.M_sun / K.Mpc) / K.kpc, 2.6),
    ("03", "(ν = 1,007)", K.verstaerkung(25 * K.a0), 1.007),
    ("03", "(ν = 3,3)", K.verstaerkung(0.13 * K.a0), 3.3),
    ("03", "Mit Υ_V = 0,93, geeicht", D.kugelhaufen_reste()[0], 0.93),
    ("03", "an den dichten (g_N > 3·a₀, 28)", D.kugelhaufen_reste()[1], 28),
    ("03", "dünnen (g_N < a₀, 17)", D.kugelhaufen_reste()[2], 17),
    ("03", "mit der Regel +0,10 ± 0,03 dex über", D.kugelhaufen_reste()[3], 0.10),
    ("03", "mit der Regel +0,10 ± 0,03 dex über", D.kugelhaufen_reste()[4], 0.03),
    ("03", "ohne Halo, +0,41 (", D.kugelhaufen_reste(True)[3], 0.41),
    ("03", "(Υ_V = 1,08, χ² 53", D.kugelhaufen_reste(True)[0], 1.08),
    ("03", "(Υ_V = 1,08, χ² 53", D.kugelhaufen_reste()[5], 53),
    ("03", "χ² 53 gegen 235)", D.kugelhaufen_reste(True)[5], 235),
    ("03", "Die gut gemessenen (11)", D.kugelhaufen_reste(False, "gut")[2], 11),
    ("03", "liegen +0,03 gegen +0,34", D.kugelhaufen_reste(False, "gut")[3], 0.03),
    ("03", "liegen +0,03 gegen +0,34", D.kugelhaufen_reste(True, "gut")[3], 0.34),
    ("03", "10 kpc vom Zentrum (9) +0,06", D.kugelhaufen_reste(False, "halo")[2], 9),
    ("03", "10 kpc vom Zentrum (9) +0,06", D.kugelhaufen_reste(False, "halo")[3], 0.06),
    ("03", "(9) +0,06 gegen +0,43", D.kugelhaufen_reste(True, "halo")[3], 0.43),
    ("03", "Database (52 mit gemessener", len(D.kugelhaufen()), 52),
    ("03", "(Rangkorrelation −0,51 an 14 Haufen", D.kugelhaufen_gezeit()[1], -0.51),
    ("03", "(Rangkorrelation −0,51 an 14 Haufen", D.kugelhaufen_gezeit()[0], 14),
    ("03", "am heutigen Abstand zu 0,13 bis 0,63 füllt", D.kugelhaufen_gezeit()[3], 0.13),
    ("03", "am heutigen Abstand zu 0,13 bis 0,63 füllt", D.kugelhaufen_gezeit()[4], 0.63),
    ("03", "zu 0,13 bis 0,63 füllt (+0,31)", D.kugelhaufen_gezeit()[2], 0.31),
    ("03", "a₀/(πG) = 242 M☉/pc²", K.flaechendichte_schwelle() * K.pc**2 / K.M_sun, 242),
    ("03", "a_N = 2,54·a₀", K.maximum_zusatz()[0], 2.54),
    ("03", "beträgt 0,648·a₀", K.maximum_zusatz()[1], 0.648),
    ("07", "Rosetten mit 254,6°", 360 / K.umlaufphase(1e-300), 254.6),
    ("07", "e^{−7480}", K.schwellenradius(K.M_sun) / K.AU, 7480),
    ("08", "bei der Milchstraße 0,42″", math.degrees(K.ablenkung_grenz(6e10 * K.M_sun)) * 3600, 0.42),
    ("08", "ε_L = ℓ_P/ℓ₀ = 1,9·10⁻⁶²", K.leiterzahl() / 1e-62, 1.9),
    ("08", "ε_L^{1/3}·m_Pl = 5,8·10⁻²⁹ kg", K.schnittmassen()[1] / 1e-29, 5.8),
    ("08", "c³τ_D/G = 1,1·10⁵⁴ kg", K.schnittmassen()[3] / 1e54, 1.1),
    ("08", "Δs(M_U) = ⅔·Λ", K.faltung(K.schnittmassen()[3]) / K.schliesslaenge() * 3, 2),
    ("08", "(4π⁴/15)·ε_L² = 0,94·10⁻¹²²", 4 * math.pi**4 / 15 * K.leiterzahl()**2 / 1e-122, 0.94),
    ("02", "r* = 32 fm", K.schwellenradius(K.m_p) * 1e15, 32),
    ("02", "a_B = 1631·r*", K.bohr_radius() / K.schwellenradius(K.m_p), 1631),
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
    ("03", "(ω_b = 0,02238", _cmb()["linse"][0]["punkt"]["wb"], 0.02238),
    ("03", "ω_c = 0,1199, H₀", _cmb()["linse"][0]["punkt"]["wc"], 0.1199),
    ("03", "H₀ = 67,37). An", _cmb()["linse"][0]["punkt"]["H0"], 67.37),
    ("03", "trifft Plancks Werte", abs(_cmb()["linse"][0]["punkt"]["wb"] - K.omega_b) / K.omega_b_fehler, 0),
    ("03", "trifft Plancks Werte", abs(_cmb()["linse"][0]["punkt"]["wc"] - K.omega_c) / K.omega_c_fehler, 0),
    ("03", "allein liegt das Modell Δχ² = 1,0 über", _cmb()["spektren"][1]["chi2"]["gesamt"] - _cmb()["spektren"][0]["chi2"]["gesamt"], 1.0),
    ("03", "an allen Daten 0,01.", _cmb()["alle"][1]["chi2"]["gesamt"] - _cmb()["alle"][0]["chi2"]["gesamt"], 0.01),
    ("03", "an allen Daten 0,01.", _cmb()["nach"] - _cmb()["alle"][0]["chi2"]["gesamt"], 0.01),
    ("03", "DESI mit χ² = 18,6", _cmb()["spektren"][1]["chi2"]["bao"], 18.6),
    ("03", "18,6 (13 Werte)", float(len(D.desi()[0])), 13),
    ("03", "ΛCDM mit 31,4", _cmb()["spektren"][0]["chi2"]["bao"], 31.4),
    ("03", "liest die Schwelle 1,081 ± 0,003", _cmb()["alle"][1]["punkt"]["a0"], 1.081),
    ("03", "liest die Schwelle 1,081 ± 0,003", _cmb()["alle"][1]["fehler"]["a0"], 0.003),
    ("03", "der Standardwert liegt Δχ² = 52,6", _cmb()["fest"]["chi2"]["gesamt"] - _cmb()["alle"][1]["chi2"]["gesamt"], 52.6),
    ("03", "gewinnt das Modell Δχ² = 9,6", _cmb()["alle"][1]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 9.6),
    ("00", "w_a < 0, Δχ² = 9,6)", _cmb()["alle"][1]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 9.6),
    ("03", "w₀ = −0,84 ± 0,05", _cmb()["zustand"]["punkt"]["w0"], -0.84),
    ("03", "w₀ = −0,84 ± 0,05", _cmb()["zustand"]["fehler"]["w0"], 0.05),
    ("03", "w_a = −0,56 ± 0,18", _cmb()["zustand"]["punkt"]["wa"], -0.56),
    ("03", "w_a = −0,56 ± 0,18", _cmb()["zustand"]["fehler"]["wa"], 0.18),
    ("03", "bei z = 0,40 ± 0,07", _cmb()["z_max"][0], 0.40),
    ("03", "bei z = 0,40 ± 0,07", _cmb()["z_max"][1], 0.07),
    ("03", "beschleunigten Dehnung (z = 0,65)", _cmb()["z_besch"], 0.65),
    ("03", "Faltung (z = 0,31)", _cmb()["z_eq"], 0.31),
    ("03", "DESI liegt dann bei χ² = 9,3", _cmb()["zustand"]["chi2"]["bao"], 9.3),
    ("03", "die Schwelle liest 1,070 ± 0,013", _cmb()["zustand"]["punkt"]["a0"], 1.070),
    ("03", "die Schwelle liest 1,070 ± 0,013", _cmb()["zustand"]["fehler"]["a0"], 0.013),
    ("03", "der Standardwert nur Δχ² = 0,8 daneben", _cmb()["zustand_fest"]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 0.8),
    ("03", "liegt Δχ² = 8,8 darunter", _cmb()["alle"][0]["chi2"]["gesamt"] - _cmb()["zustand_fest"]["chi2"]["gesamt"], 8.8),
    ("03", "liegt noch 2,5 tiefer", _cmb()["zustand_fest"]["chi2"]["gesamt"] - _cmb()["zustand_lcdm"]["chi2"]["gesamt"], 2.5),
    ("03", "nur auf 1,3 % genau", 100 * _cmb()["zustand"]["fehler"]["a0"] / _cmb()["zustand"]["punkt"]["a0"], 1.3),
    ("03", "treffen die Daten das ebenso (Δχ² = 0,9)", _cmb()["saett"]["chi2"]["gesamt"] - _cmb()["zustand"]["chi2"]["gesamt"], 0.9),
    ("03", "die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["punkt"]["a0"], 1.085),
    ("03", "die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["fehler"]["a0"], 0.003),
    ("03", "Das Maximum liegt dann bei z = 0,32", _cmb()["z_saett"], 0.32),
    ("03", "geht der Standardwert 1,06 nicht (Δχ² = 48)", _cmb()["saett_fest"]["chi2"]["gesamt"] - _cmb()["zustand_fest"]["chi2"]["gesamt"], 48),
    ("03", "wie ΛCDM und liegt Δχ² = 8,7 darunter", _cmb()["alle"][0]["chi2"]["gesamt"] - _cmb()["saett"]["chi2"]["gesamt"], 8.7),
    ("03", "treffen sie ebenso (χ² = 2019,7 gegen 2019,1)", _cmb()["saett_form"]["chi2"]["gesamt"], 2019.7),
    ("03", "treffen sie ebenso (χ² = 2019,7 gegen 2019,1)", _cmb()["saett"]["chi2"]["gesamt"], 2019.1),
    ("03", "der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils", 100 * _cmb()["zustand_form"][0][0], 5.4),
    ("03", "der Zustandsteil ist heute 5,4 ± 0,7 % des Formteils", 100 * _cmb()["zustand_form"][0][1], 0.7),
    ("00", "Mit dem Standardwert als Formteil ist der Zustandsteil heute 5,4 %", 100 * _cmb()["zustand_form"][0][0], 5.4),
    ("F:ganzes", "F = 0,97 ± 0,12 %", 100 * _cmb()["zustand_form"][1][0], 0.97),
    ("F:ganzes", "F = 0,97 ± 0,12 %", 100 * _cmb()["zustand_form"][1][1], 0.12),
    ("F:ganzes", "am Maximum (z = 0,30) F = 1,68 ± 0,36 %", 100 * _cmb()["zustand_form"][2][0], 1.68),
    ("F:ganzes", "am Maximum (z = 0,30) F = 1,68 ± 0,36 %", 100 * _cmb()["zustand_form"][2][1], 0.36),
    ("F:ganzes", "am Maximum (z = 0,30) F = 1,68", _cmb()["zustand_form"][3][0], 0.30),
    ("F:ganzes", "Ω_HI·h₇₀ = (3,9 ± 0,6)·10⁻⁴", K.KALTES_GAS["hi"][0] * 1e4, 3.9),
    ("F:ganzes", "Ω_HI·h₇₀ = (3,9 ± 0,6)·10⁻⁴", math.hypot(*K.KALTES_GAS["hi"][1:]) * 1e4, 0.6),
    ("F:ganzes", "Ω_H₂·h = 5,3·10⁻⁵", K.KALTES_GAS["h2"][0] * 1e5, 5.3),
    ("F:ganzes", "F = 1,32 ± 0,17 % gegen 0,97 ± 0,12 %", 100 * _kaltes_gas()[0], 1.32),
    ("F:ganzes", "F = 1,32 ± 0,17 % gegen 0,97 ± 0,12 %", 100 * _kaltes_gas()[1], 0.17),
    ("F:ganzes", "dieselbe Größe, 1,7 σ auseinander", (_kaltes_gas()[0] - _cmb()["zustand_form"][1][0])
     / math.hypot(_kaltes_gas()[1], _cmb()["zustand_form"][1][1]), 1.7),
    ("F:ganzes", "1,7-mal so viel kaltes Gas wie heute", _cmb()["zustand_form"][2][0] / _cmb()["zustand_form"][1][0], 1.7),
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
    ("00", "lesen alle Präzisionsdaten Q = 5,36 ± 0,05, 2π − 1", _cmb()["q_lesen"]["punkt"]["Q"], 5.36),
    ("00", "2π − 1 liegt Δχ² = 2,4 daneben, die Haufen", _cmb()["saett"]["chi2"]["gesamt"] - _cmb()["q_lesen"]["chi2"]["gesamt"], 2.4),
    ("00", "lesen alle Präzisionsdaten die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["punkt"]["a0"], 1.085),
    ("00", "lesen alle Präzisionsdaten die Dunkle Energie heute als a₀ = 1,085 ± 0,003", _cmb()["saett"]["fehler"]["a0"], 0.003),
    ("00", "ΛCDM und Δχ² = 8,7 darunter", _cmb()["alle"][0]["chi2"]["gesamt"] - _cmb()["saett"]["chi2"]["gesamt"], 8.7),
    ("03", "liest ΛCDM selbst Q = 5,278", _cmb()["alle"][0]["punkt"]["wc"] / _cmb()["alle"][0]["punkt"]["wb"], 5.278),
    ("03", "an den Spektren allein 5,367", _cmb()["spektren"][0]["punkt"]["wc"] / _cmb()["spektren"][0]["punkt"]["wb"], 5.367),
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
    ("02", "bei 126 (2,03 statt", KE.streuung_ab(KE.NIVEAUS, 0, False), 2.03),
    ("02", "2,03 statt 1,62 MeV", KE.streuung_ab(KE.SCHALE, 0, False), 1.62),
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
]


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


def pruefe_modellzahlen():
    zeilen = []
    for nr, stelle, wert, soll in MODELLZAHLEN:
        if wert is None:
            wert = _zelle()
        steht = stelle in _modelltext(nr)
        rund = round(wert, _stellen(soll))
        passt = abs(rund - soll) <= 10**(-_stellen(soll)) * 0.5 + 1e-12
        zeilen.append((nr, stelle, wert, steht and passt, steht))
    zeilen += pruefe_abhandlung()
    zeilen += pruefe_rechnen()
    zeilen += pruefe_elemente()
    return zeilen


def pruefe_abhandlung():
    """Jede Dezimalzahl der Abhandlung steht im geprüften Modell (Zahlen nur aus dem Modell)."""
    modell = "".join(open(os.path.join(WURZEL, "modell", f), encoding="utf-8").read()
                     for f in os.listdir(os.path.join(WURZEL, "modell")))
    zeilen = []
    for name, muster, komma in (("abhandlung_de.md", r"\d+,\d+", False), ("abhandlung_en.md", r"\d+\.\d+", True)):
        pfad = os.path.join(WURZEL, "abhandlung", name)
        if not os.path.exists(pfad):
            continue
        zahlen = sorted(set(re.findall(muster, open(pfad, encoding="utf-8").read())))
        fehlt = [z for z in zahlen if (z.replace(".", ",") if komma else z) not in modell]
        zeilen += [("AB", f"{name}: Zahl {z} steht im Modell", 0.0, False, False) for z in fehlt] or \
                  [("AB", f"{name}: alle {len(zahlen)} Zahlen stehen im Modell", float(len(zahlen)), True, True)]
    return zeilen


def pruefe_rechnen():
    """Jede gerechnete Zahl der Beispiele steht in beiden Sprachen; jede Dezimalzahl dort ist gerechnet oder im Modell."""
    gerechnet = B.zahlen()
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


def pruefe_elemente():
    """abhandlung/elemente_de.md ist genau die Rechnung von elemente.py."""
    import elemente
    aus = []
    for sp in ("de", "en"):
        gleich = open(elemente.ZIEL[sp], encoding="utf-8").read() == elemente.text(sp)
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
    for fn in szenarien.ALLE:
        erg = fn()
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
    for nrb, stelle, wert, ok, steht in pruefe_modellzahlen():
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
