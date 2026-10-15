# -*- coding: utf-8 -*-
"""
Daten: der eine Zugang zu weltmodell/daten (SPARC, AME2020, PubChem). Jede Datei wird vor dem ersten Lesen gegen ihren
Hash in kern/daten_hashes.json geprüft (Datenquelle.pfad); eine geänderte Datei bricht ab. Gelesenes bleibt im Speicher
(ein Lesen je Datei und Lauf).

Datensatz SPARC (Lelli, McGaugh, Schombert 2016): Tafel SPARC_Lelli2016c.mrt, Kurven Rotmod_LTG/<Name>_rotmod.dat.
Einheiten im Datensatz: Längen und Geschwindigkeiten in SI (m, m/s), Neigung in Grad, Leuchtkraft L_3,6 in L☉[3,6],
SB₀ in L☉/pc², M_HI in M☉. Die Auswahl (Qualität, Neigung) und die Festlegungen Υ_b und Helium stehen in SPARC_SATZ;
jede Rechnung nennt sie beim Aufruf.
"""
import csv
import hashlib
import os

import numpy as np

from .konstanten import kpc, Mpc, e_el
from .formeln import formel
from . import datenhashes as _dh
from .setzungen import SETZUNGEN as _S

# Festlegungen des Datensatzes SPARC (Setzungen): Auswahl nach Güte und Neigung, Υ_b, Helium, Streuung von Υ.
SPARC_SATZ = dict(q_max=_S["SPARC_Q_MAX"].wert, inc_min=_S["SPARC_INC_MIN"].wert, ups_b=_S["SPARC_UPS_B"].wert,
                  helium=_S["SPARC_HELIUM"].wert, ups_streuung_dex=_S["SPARC_UPS_STREUUNG"].wert)


class Datenquelle:
    """Datei aus weltmodell/daten mit Hash-Prüfung; `pfad(rel)` gibt den geprüften Pfad, `praefix(p)` die eine Datei
    der Tafel, deren Pfad mit p beginnt (Dateien, deren Endung der Wächter für Ausgaben hält, werden so benannt)."""

    def __init__(self, wurzel=_dh.DATEN):
        self.wurzel, self.geprueft, self.speicher, self._soll = wurzel, set(), {}, None

    def soll(self):
        if self._soll is None:
            self._soll = _dh.gespeichert()
        return self._soll

    def praefix(self, p):
        treffer = [k for k in self.soll() if k.startswith(p)]
        if len(treffer) != 1:
            raise FileNotFoundError(f"Daten: {p!r} trifft {treffer}")
        return treffer[0]

    def pfad(self, rel):
        voll = os.path.join(self.wurzel, *rel.split("/"))
        if rel not in self.geprueft:
            soll = self.soll().get(rel)
            if soll is None:
                raise FileNotFoundError(f"Daten: {rel} nicht in der Hash-Tafel")
            with open(voll, "rb") as f:
                ist = hashlib.sha256(f.read()).hexdigest()
            if ist != soll:
                raise ValueError(f"Daten: Hash von {rel} geändert")
            self.geprueft.add(rel)
        return voll

    def sparc_zeilen(self):
        """Name → Zeile der SPARC-Tafel (in der Reihenfolge der Datei)."""
        if "tabelle" not in self.speicher:
            t = {}
            with open(self.pfad("sparc/SPARC_Lelli2016c.mrt"), encoding="latin-1") as f:
                for zeile in f:
                    x = zeile.split()
                    if len(x) < 18:
                        continue
                    try:
                        t[x[0]] = dict(T=int(x[1]), D=float(x[2]) * Mpc, eD=float(x[3]) * Mpc, fD=int(x[4]),
                                       inc=float(x[5]), eInc=float(x[6]), L=float(x[7]) * 1e9, Reff=float(x[9]) * kpc,
                                       Rd=float(x[11]) * kpc, SB0=float(x[12]), MHI=float(x[13]) * 1e9,
                                       Vf=float(x[15]) * 1e3, eVf=float(x[16]) * 1e3, Q=int(x[17]))
                    except ValueError:
                        pass
            self.speicher["tabelle"] = t
        return self.speicher["tabelle"]


QUELLE = Datenquelle()


# ------------------------------------------------------------------ SPARC
@formel('sparc_tabelle', 'Zeile der SPARC-Tafel einer Galaxie (Abstand, Neigung, L_3,6, R_d, SB₀, M_HI, V_f, Q)',
        ('D_A', 'i_inc', 'L_3,6', 'R_d', 'SB₀', 'M_HI', 'V_f', 'Q_fl'), (), 'L', 'hilfe', 'SPARC (Lelli, McGaugh, Schombert 2016)')
def sparc_tabelle(name):
    """Zeile der Tafel SPARC_Lelli2016c.mrt als dict: T, D, eD [m], fD, inc, eInc [°], L [L☉], Reff, Rd [m], SB0 [L☉/pc²],
    MHI [M☉], Vf, eVf [m/s], Q."""
    return dict(QUELLE.sparc_zeilen()[name])


@formel('sparc_kurve', 'Rotationskurve R, V_obs, e_V, V_gas, V_disk, V_bul (Υ = 1) einer Galaxie',
        ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul'), ('r',), 'L', 'hilfe', 'SPARC (Lelli, McGaugh, Schombert 2016)')
def sparc_kurve(name):
    """dict mit Arrays R [m], Vo, eV, Vg, Vd, Vb [m/s], SBd, SBb [L☉/pc², Flächenhelligkeit von Scheibe und Bulge bei
    3,6 μm] aus Rotmod_LTG/<name>_rotmod.dat."""
    k = ("kurve", name)
    if k not in QUELLE.speicher:
        a = np.loadtxt(QUELLE.pfad(f"sparc/Rotmod_LTG/{name}_rotmod.dat"), comments="#", ndmin=2)
        QUELLE.speicher[k] = dict(R=a[:, 0] * kpc, Vo=a[:, 1] * 1e3, eV=a[:, 2] * 1e3, Vg=a[:, 3] * 1e3,
                                  Vd=a[:, 4] * 1e3, Vb=a[:, 5] * 1e3, SBd=a[:, 6], SBb=a[:, 7])
    return {n: v.copy() for n, v in QUELLE.speicher[k].items()}


@formel('sparc_tafel', 'volle Tafel: je Galaxie mit Q ≤ q_max und Neigung ≥ inc_min Tafelzeile und Kurve',
        ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul'), ('Q_fl', 'i_inc'), 'L', 'hilfe',
        'SPARC (Lelli, McGaugh, Schombert 2016)', alternative='sparc_kurve')
def sparc_tafel(q_max, inc_min):
    """Liste von dicts (Name 'nm', Tafelzeile und Kurve), in der Reihenfolge der Tafel. Standard der bisherigen
    Rechnungen: sparc_tafel(SPARC_SATZ['q_max'], SPARC_SATZ['inc_min'])."""
    aus = []
    for nm in QUELLE.sparc_zeilen():
        d = sparc_tabelle(nm)
        if d["Q"] > q_max or d["inc"] < inc_min:
            continue
        aus.append(dict(nm=nm, **d, **sparc_kurve(nm)))
    return aus


# ------------------------------------------------------------------ AME2020
@formel('ame2020', 'Massenexzess, Bindungsenergie je Nukleon und Atommasse eines Nuklids (Z, A) aus AME2020',
        ('E_B/A', 'Δ_M'), ('Z_A', 'A_A'), 'L', 'hilfe', 'AME2020 (Huang, Wang, Kondev, Audi, Naimi 2021)')
def ame2020(Z, A):
    """dict: el, N, Z, A, excess [J], BA [J je Nukleon], masse_u [u], gemessen (False bei geschätzten Werten '#')."""
    if "ame" not in QUELLE.speicher:
        t = {}
        with open(QUELLE.pfad(QUELLE.praefix("elemente/mass_1.mas20")), encoding="latin-1") as f:
            zeilen = f.read().splitlines()
        start = next(i for i, z in enumerate(zeilen) if z.startswith("1N-Z")) + 2
        for z in zeilen[start:]:
            if len(z) < 123 or not z[9:14].strip():
                continue
            gemessen = "#" not in z[28:79]
            fl = lambda s: float(s.replace("#", ".").strip()) if s.strip() not in ("", "*") else float("nan")  # noqa: E731
            N, Zk, Ak = int(z[4:9]), int(z[9:14]), int(z[14:19])
            t[(Zk, Ak)] = dict(el=z[20:23].strip(), N=N, Z=Zk, A=Ak, excess=fl(z[28:42]) * 1e3 * e_el,
                               BA=fl(z[54:67]) * 1e3 * e_el, masse_u=int(z[106:109]) + fl(z[110:123]) * 1e-6,
                               gemessen=gemessen)
        QUELLE.speicher["ame"] = t
    return dict(QUELLE.speicher["ame"][(int(Z), int(A))])


# ------------------------------------------------------------------ PubChem
@formel('pubchem', 'Stoffwerte eines Elements: Schmelz- und Siedetemperatur, erste Ionisierungsenergie (PubChem)',
        ('T_schm', 'T_sied', 'E_ion'), ('Z_A',), 'L', 'hilfe', 'PubChem, Periodensystem der Elemente')
def pubchem(Z):
    """dict: symbol, name, masse_u, T_schm [K], T_sied [K], E_ion [J], zustand (nan, wo PubChem keinen Wert hat)."""
    if "pubchem" not in QUELLE.speicher:
        t = {}
        with open(QUELLE.pfad("elemente/pubchem_periodensystem.csv"), encoding="utf-8") as f:
            for r in csv.DictReader(f):
                fl = lambda s: float(s) if s not in ("", None) else float("nan")  # noqa: E731
                t[int(r["AtomicNumber"])] = dict(symbol=r["Symbol"], name=r["Name"], masse_u=fl(r["AtomicMass"]),
                                                 T_schm=fl(r["MeltingPoint"]), T_sied=fl(r["BoilingPoint"]),
                                                 E_ion=fl(r["IonizationEnergy"]) * e_el, zustand=r["StandardState"])
        QUELLE.speicher["pubchem"] = t
    return dict(QUELLE.speicher["pubchem"][int(Z)])
