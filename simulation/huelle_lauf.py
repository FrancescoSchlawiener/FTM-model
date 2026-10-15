"""Lauf der Hüllenrechnung in Etappen. Jede Etappe liest, was die vorigen abgelegt haben (huelle_ablage), und überspringt,
was schon gerechnet ist. Ohne Angabe: Stand jeder Etappe. Mit Namen: diese Etappen rechnen, in der Reihenfolge der Liste.

    python huelle_lauf.py                         Stand
    python huelle_lauf.py entstehung kreis
    python huelle_lauf.py entstehung 24 29        nur diese Elemente, neu gerechnet"""
import sys

import huelle_2_entstehung as K2
import kerne as K
import huelle_5_fassungen as K5
import huelle_3_eigenschaften as K3
from huelle_ablage import lesen

# (Name, Schlüssel der Ablage, Rechnung)
ETAPPEN = [
    ("entstehung", "_entstehung", lambda z=None: K2.rechne_entstehung(z, neu=True) if z else K2.rechne_entstehung()),                          # Kreis: Hülle beim Fallen von δ
    ("paar_he", "_paar_he", lambda z=None: K2.paar_probe(2, {(1, 0): 2})),                                   # Einrasten 1s²
    ("paar_pr", "_paar_pr", lambda z=None: K2.paar_probe(59, K.belegung_gemessen()[59], (200,), (6,), "band")),  # Einrasten Pr
    ("pruefung_kreis", "_pruefung_kreis", lambda z=None: K2.rechne_pruefe_zwei(                            # gemessen gegen
        z or (59, 65, 66, 67, 68, 69, 70, 71, 74, 91, 97))),                                                 # Kreis, zwei Zählweisen
    ("pruefung_kreis_j", "_pruefung_kreis_j", lambda z=None: K2.rechne_pruefe_zwei(                       # dazu die Lage
        z or (59, 65, 66, 67, 68, 69, 70, 71, 74, 91, 97), schluessel="_pruefung_kreis_j", verwindung=True)),  # des Niveaus J
    ("entstehung_zwei", "_entstehung_zwei", lambda z=None: K2.rechne_entstehung(     # Kreis, zwei Zählweisen, Niveau J
        z or range(1, 103), prozesse=2, schluessel="_entstehung_zwei", zaehlweisen=True, verwindung=True)),
    ("abstand", "_abstand", lambda z=None: K2.rechne_abstand(z, prozesse=2)),                  # Prüfung: gemessen gegen gerechnet
    ("rand", "_rand", lambda z: K2.rechne_rand(z)),
    ("rand_paar", "_rand_paar", lambda z: K2.rechne_rand(z, schluessel="_rand_paar", paarung=True)),       # mit eingerasteten Paaren
    ("rand_band", "_rand_band", lambda z: K2.rechne_rand(z, schluessel="_rand_band", paarung=True, einrasten="band")),
    ("rand_ganz", "_rand_ganz", lambda z: K2.rechne_rand(z, schluessel="_rand_ganz", paarung=True, einrasten="band", grundterm=True,
                                                         verwindung=True)),          # Einrasten, Grundterm, Verwindung                                         # Prüfung: Nachbarn der gemessenen
    ("abstand_windung", "_abstand_windung", lambda z: K2.rechne_abstand(z, prozesse=2, schluessel="_abstand_windung", je_windung=True)),
    ("stern", "_stern", lambda: K2.stern()),                                                # Kette im Stern: Tiefe, Druck, Zelle
    ("kreis", "_kreis", lambda: K3.kreis()),
    # Stand 05.10.: Konfiguration aus zwei Zählweisen (örtlich und je Umlaufsinn, Gewicht Füllung der Spinhälfte mal Tempo,
    # gespiegelt um die Mitte), am eingerasteten Niveau J (gekoppelte Lage), exakter Teil auf h = 0,016 bzw. gegen das Gitter null
    ("suche", "_suche", lambda z=None: K5.rechne_suche_zwei(gekoppelt=True, h=0.016, schluessel="_suche")),
    ("pruefung", "_pruefung", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(gekoppelt=True, h=0.016, schluessel="_pruefung"))),
    ("pruefung_grenz", "_pruefung_grenz", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(gekoppelt=True, h=0, schluessel="_pruefung_grenz"))),
    ("gitter_np", "_gitter_np", lambda z=None: print("Np gegen 5f5 7s2 je Gitter:", K5.abstand_gitter(93, "5f57s2"))),
    ("takt", "_takt", lambda z=None: [print(Z, k, K5.takt_abstand(Z, k), flush=True) for Z, k in               # fester Takt,
        ((74, "4f145d56s1"), (29, "3d94s2"), (23, "3d44s1"), (93, "5f57s2")) if not z or Z in z]),                # Tempo, nur außen
    ("pruefung_takt", "_pruefung_takt", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(
        gekoppelt=True, h=0.016, z_bereich=z or (93, 29, 23),
        schluessel="_pruefung_takt", skaliert=True))),        # Takt aus dem Tempo, die drei knappsten
    ("pruefung_ok", "_pruefung_ok", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(     # ohne angepasste Korrelation:
        gekoppelt=True, h=0.016, z_bereich=z or (59, 66, 46, 74, 93, 29, 23, 57, 64, 65, 70, 90, 91, 58, 78, 77),  # nur die
        prozesse=2, schluessel="_pruefung_ok", korrelation=False))),                                        # Zählung der Sperre
    ("pruefung_ok_paar", "_pruefung_ok_paar", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(   # ohne angepasste
        gekoppelt=True, h=0.016, z_bereich=z or range(1, 119), prozesse=2,                                # Zahl, mit dem
        schluessel="_pruefung_ok_paar", korrelation=False, paarung=True))),                               # Einrasten der Paare
    ("pruefung_ok_mode", "_pruefung_ok_mode", lambda z=None: print("Abstand ≤ 0:", K5.rechne_pruefung(   # ohne angepasste
        gekoppelt=True, h=0.016, prozesse=2, schluessel="_pruefung_ok_mode", korrelation=False, paarung=True,   # Zahl, je Mode
        einrasten="mode", z_bereich=z or (*(59, 28, 45, 27, 43, 74, 46, 29, 23, 57, 58, 65, 66, 70, 90, 93, 64, 77, 78, 91),
                                          *range(1, 119))))),                                           # ein Paar
    ("ueber_118", "_suche", lambda z=None: print("Ende bei Z =", K5.rechne_ueber_118(gekoppelt=True, h=0.016, schluessel="_suche"))),                                                # Entstehung, Kern, Hülle, Rückwirkung
]


def stand():
    alle = lesen()
    for name, schluessel, _ in ETAPPEN:
        d = {k: v for k, v in alle.items() if not k.startswith("_")} if schluessel is None else alle.get(schluessel, {})
        print(f"{name:24s} {len(d):4d} Einträge")


if __name__ == "__main__":
    namen = dict((n, f) for n, _, f in ETAPPEN)
    if len(sys.argv) == 1:
        stand()
    z = [int(a) for a in sys.argv[1:] if a.isdigit()]
    for n in [a for a in sys.argv[1:] if not a.isdigit()]:
        print("Etappe", n, flush=True)
        namen[n](z) if z else namen[n]()
