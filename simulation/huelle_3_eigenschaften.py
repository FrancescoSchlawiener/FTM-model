"""Etappe 3: Eigenschaften aus der Konfiguration: Grundterm und Moment nach Hund, Schwellen der Kopplung der Hülle an das
Ganze (Makrozustand) aus den gemessenen Daten."""
import math
import os
import kern as K


def grundterm(bel):
    """Grundterm der Konfiguration nach Hund über alle offenen Unterschalen: S und L je Schale (größtes S, dann größtes L),
    summiert; J = |L − S| unter halber Füllung der offenen Schale größter Windung, sonst L + S (ist sie genau halb, entscheidet
    die nächste). Landé-Faktor g (g_s = 2), Moment g·J und g·√(J(J + 1)) in Bohr-Magnetonen."""
    S = L = 0.0
    offen = []
    for (n, l), m in bel.items():
        cap = 2 * (2 * l + 1)
        if 0 < m < cap:
            auf, ab = min(m, 2 * l + 1), m - min(m, 2 * l + 1)
            S += 0.5 * (auf - ab)
            L += abs(sum(range(l, l - auf, -1)) + sum(range(l, l - ab, -1)))
            offen.append((l, m / cap))
    if not offen:
        return dict(S=0.0, L=0.0, J=0.0, g=0.0, mu=0.0, mu_eff=0.0)
    unter = None
    for l, a in sorted(offen, reverse=True):
        if a != 0.5:
            unter = a < 0.5
            break
    J = abs(L - S) if unter in (True, None) else L + S
    g = 1 + (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1)) if J > 0 else 0.0
    return dict(S=S, L=L, J=J, g=g, mu=g * J, mu_eff=g * math.sqrt(J * (J + 1)))

def makro_schwellen():
    """Kopplung der Hülle an das Ganze, je Element mit gemessener Ionisierungsenergie (PubChem): Die Atome falten als Masse den
    Körper (δ = GM/(Rc²) = 3Δs/R); der Körper gibt jedem Atom Druck, P·V_A = ½·m_A·c²·δ (homogene Kugel, Zentrum), und als
    wärmegehaltener Trichter Wärme, kT = f·δ·m̄c² (f aus der Polytrope n = 3, m̄ = m_A). Die Hülle hält, solange beides unter
    ihrer Bindung IE bleibt: δ_P = 2·IE/(m_A·c²), δ_T = IE/(f·m_A·c²). m_A und IE gemessen (PubChem)."""
    import csv
    import daten
    zeilen = csv.DictReader(open(os.path.join(daten.ORT, "elemente", "pubchem_periodensystem.csv"), encoding="utf-8"))
    f3 = K.polytrope(3.0)[2]
    aus = {}
    for x in zeilen:
        if not x["IonizationEnergy"]:
            continue
        Z, A, ie = int(x["AtomicNumber"]), float(x["AtomicMass"]), float(x["IonizationEnergy"])
        mc2 = A * K.m_u * K.c ** 2
        aus[Z] = dict(A=A, ie=ie, d_P=2 * ie * K.e / mc2, d_T=ie * K.e / (f3 * mc2))
    return aus


def kreis(z_bereich=range(1, 103)):
    """Der ganze Kreis je Element, mit eigenen Zahlen: Ort der Entstehung des Kerns (Brennleiter bis x = ¼, darüber Einfang
    von Neutronen), Kern im Tal (A, Kippstelle x, Radius R = λ_H·A^{1/3}), Hülle aus der Etappe entstehung mit der Bindung
    des zuletzt gebundenen Elektrons im Spiegelfeld der übrigen, E_b, und ihre Rückwirkung auf das Ganze: δ_P = 2·E_b/(m_A·c²),
    δ_T = E_b/(f·m_A·c²). Gespeichert unter '_kreis'."""
    import kerne
    from huelle_ablage import lesen, schreiben
    alle = lesen()
    ent = alle.get("_entstehung", {})
    f3 = K.polytrope(3.0)[2]
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN))
    z_leiter = kerne.ende_brennleiter(kerne.TROPFEN)[0]
    aus = {}
    for Z in z_bereich:
        if str(Z) not in ent:
            continue
        A = float(kerne.a_tal(Z, kerne.TROPFEN))
        eb = ent[str(Z)]["bindung_ev"]
        mc2 = A * K.m_u * K.c ** 2
        aus[str(Z)] = dict(ort="Brennleiter" if Z <= z_leiter else "Neutroneneinfang", A=A, x=float(kerne.x(Z, A, kerne.TROPFEN)),
                           R_fm=lam * A ** (1 / 3), konf=ent[str(Z)]["konf"], bindung_ev=eb,
                           d_P=2 * eb * K.e / mc2, d_T=eb * K.e / (f3 * mc2))
    alle["_kreis"] = aus
    schreiben(alle)
    return aus
