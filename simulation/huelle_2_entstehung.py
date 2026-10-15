"""Etappe 2: Entstehung der Hülle als Kreis, während δ des Ganzen fällt. Der Kern Z ist entstanden (gehaltener Trichter,
R_K = λ_H·A^{1/3}); er hält k Elektronen und fängt das nächste. Jede Bahn schließt je Umlauf, tief im Trichter mit dem
langsameren Takt der Abwicklung, gesetzt vom Zentrum (A3: K = 1/(1 − V_K/2c²)), in der Neigung des Kerns und der
Windungen der übrigen Gebundenen (1/r² je Windung, ohne die eigene); gleichsinnige Windungen am selben Ort sperren einander
ihren Zustand (𝒫), das gesperrte Maß wird nicht gespürt; Kreis: Bahnen → Windungen → Neigung → Bahnen, bis beides
steht (Rückkehr Σqⁿ). Bindung sperrt Moden (Schließungsoperator 𝒫): Ein einmal gebundenes Elektron bleibt gebunden; geschützt
ist die Zählung, nicht die einzelne Schleife. Wo die Zählregel schließt (ein Umlaufsinn einer Unterschale ganz oder beide), ist
fest gebunden; die übrigen wandern nach jedem Fang, solange die Faltung der ganzen Hülle sinkt. Fällt δ, fängt der Kern mit k
Gebundenen das nächste Elektron dort, wo die Faltung der ganzen Hülle am tiefsten wird (das Ganze steuert das Einzelne), wenn kT = f·δ·m_A·c² unter dessen Bindung fällt.
Die Konfiguration ist die Spur dieser Fänge. Höhe v = 8·max F und Ort R des Spiegelfelds (F = −r²V) aus dem fertigen Feld.
Ohne angepassten Wert."""
import math
from functools import lru_cache

import numpy as np
from huelle_1_loeser import FM, _eigen, _gitter, _matrix, c
from huelle_ablage import lesen, schreiben

N_MAX, L_MAX = 8, 4
N_FREI = 200                      # freie Moden je Windung l für das Einrasten
PAARUNG = False                   # Bindung eingerasteter Paare in der Faltung der Hülle
L_PAAR = None                     # höchste Windung der freien Moden beim Einrasten (ohne Angabe: 2·l der äußeren Schalen)
GRUNDTERM = False                 # Windungen offener Schalen in den Moden des tiefsten Zustands (Hund)
VERWINDUNG = False                # Verwindung des tiefsten Niveaus
EINRASTEN = "mode"                # eingerastet: "mode" (beide Umlaufsinne einer Mode) oder "band" (alle gegensinnigen im Fenster)
HARTREE_EV = 27.211386
_START = {}
KREUZ_JE_WINDUNG = False          # Sperre gegen andere Schalen: je Windung (True) oder Mittel des Sinns am Ort (False)
ZAEHLWEISEN = False               # Wahl nach zwei Zählweisen: je Umlaufsinn (der Kreis) und örtlich, gewichtet und gespiegelt


def _kern(Z, r):
    """Neigung des Kerns als gehaltener Trichter: innen der gleichmäßig gehaltene Körper, außen −Z/r; R = λ_H·A^{1/3} im Tal."""
    import kerne
    R = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN)) * kerne.a_tal(Z, kerne.TROPFEN) ** (1 / 3) * FM
    return np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r)


def _takt(V, e):
    """Takt einer Windung aus ihrem eigenen Tempo am Ort (A3, auf jeder Bahn δ = v²/c²): v² = 2(ε − V), die Schließung je Umlauf
    läuft mit K = 1/(1 + (ε − V)/(2c²)). Die Phase ändert den Takt nicht: Die Tiefe der Neigung allein setzt ihn nicht, nur das
    Tempo, das sie der Windung gibt; eine gleichmäßige Verschiebung von V und ε lässt ihn stehen."""
    return 1 / (1 + (e - V) / (2 * c * c))


def _neigung(q, r):
    """Neigung einer kugelförmig verteilten Windung (q je Gitterpunkt, Σq = 1): Q(r)/r + Σ_außen q/r'."""
    return (np.cumsum(q) - 0.5 * q) / r + np.cumsum((q / r)[::-1])[::-1] - 0.5 * q / r


def _bahn(Z, V, r, h, s, e=None):
    """Bahn (n, l) in der Neigung V: (ε, φ) mit Σφ²r²h = 1, im Takt ihres Tempos bei der Faltung e; ohne e, bis ε steht."""
    for _ in range(1 if e is not None else 20):
        f = e if e is not None else (_START.get((Z, s)) or [-Z * Z / (2.0 * s[0] ** 2)])[-1]
        w, v = _eigen(V, _takt(V, f), r, h, s[1], -1, s[0] - s[1], -2.0 * Z * Z - 5, _START.get((Z, s)))
        _START[(Z, s)] = list(w)
        if e is None and abs(w[-1] - f) < 1e-10 * max(1.0, abs(f)):
            break
    return float(w[-1]), v[:, -1]


def _sinne(besetzt):
    """Elektronen je Umlaufsinn jeder Unterschale: zuerst ein Sinn ganz (2l + 1 Zustände), dann der andere."""
    return {s: (min(m, 2 * s[1] + 1), m - min(m, 2 * s[1] + 1)) for s, m in besetzt.items()}


def _yk(q, r, k):
    """Neigung der Stufe k (Kugelfunktion vom Grad k) einer Verteilung q je Gitterpunkt: r^{−k−1}∫₀ʳ q r'^k + r^k ∫ᵣ^∞ q r'^{−k−1}."""
    a, b = q * r ** k, q / r ** (k + 1)
    return (np.cumsum(a) - 0.5 * a) / r ** (k + 1) + r ** k * (np.cumsum(b[::-1])[::-1] - 0.5 * b)


@lru_cache(maxsize=None)
def _kopplung(l1, k, l2):
    """Gewicht, mit dem zwei Kugelfunktionen der Windungen l1, l2 die Stufe k bilden: (l1 k l2; 0 0 0)²."""
    from sympy.physics.wigner import wigner_3j
    return float(wigner_3j(l1, k, l2, 0, 0, 0)) ** 2


@lru_cache(maxsize=None)
def _gaunt(k, l1, m1, l2, m2):
    """Anteil der Stufe k an der Überlagerung der Mode m1 einer Kugelfunktion vom Grad l1 mit der Mode m2 vom Grad l2
    (Gaunt-Zahl c^k(l1 m1, l2 m2))."""
    from sympy.physics.wigner import wigner_3j
    return ((-1) ** m1 * math.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * float(wigner_3j(l1, k, l2, 0, 0, 0))
            * float(wigner_3j(l1, k, l2, -m1, m1 - m2, m2)))


@lru_cache(maxsize=None)
def _moden(l):
    """Gezählt über die Moden m der Kugelfunktion vom Grad l: b_k der Sperre je gleichsinnigem Paar im Mittel über alle Moden
    (⟨J − K⟩ − F⁰ für m ≠ m')."""
    ks = range(2, 2 * l + 1, 2)
    J = lambda m, mp, k: _gaunt(k, l, m, l, m) * _gaunt(k, l, mp, l, mp)
    K = lambda m, mp, k: _gaunt(k, l, m, l, mp) ** 2
    ms = range(-l, l + 1)
    paare = [(m, mp) for m in ms for mp in ms if m != mp]
    b = {k: sum(J(m, mp, k) - K(m, mp, k) for m, mp in paare) / len(paare) for k in ks} if l else {}
    return b


def _sperre(besetzt, bahn, r, h):
    """Bindung sperrt Moden (𝒫): Wo eine Windung eines Umlaufsinns läuft, ist ihr Zustand für die übrigen desselben Sinns gesperrt;
    die gesperrte Überlagerung zweier gleichsinniger Windungen a, b am selben Ort, gezählt über die Moden ihrer Kugelfunktionen,
    wird nicht gespürt. Gibt (Beitrag zur Faltung [Hartree, ≤ 0], Neigung je Zustand {s: V}). Jede Windung spürt die Moden, die sie
    selbst sperrt: mit den gleichsinnigen ihrer eigenen Unterschale an ihrem Ort, −(n_σ − 1)·d_aa/q_a; mit denen anderer Schalen
    das Mittel ihres Sinns am Ort, 2·e_σ/ρ_σ."""
    sinn = _sinne(besetzt)
    schalen = sorted(besetzt)
    quer = [np.zeros_like(r), np.zeros_like(r)]
    kreuz = {a: [np.zeros_like(r), np.zeros_like(r)] for a in schalen}
    E = 0.0
    eigen = {}
    for i, a in enumerate(schalen):
        for c in schalen[i:]:
            q = bahn[a][1] * bahn[c][1] * r * r * h
            if a == c:
                d = -sum(bk * q * _yk(q, r, k) for k, bk in _moden(a[1]).items()) if a[1] else 0.0 * r
                paare = [n * (n - 1) / 2 for n in sinn[a]]
                eigen[a] = d
            else:
                d = sum(_kopplung(a[1], k, c[1]) * q * _yk(q, r, k) for k in range(abs(a[1] - c[1]), a[1] + c[1] + 1, 2))
                paare = [sinn[a][0] * sinn[c][0], sinn[a][1] * sinn[c][1]]
                for t in (0, 1):
                    quer[t] = quer[t] - paare[t] * d
                    kreuz[a][t] = kreuz[a][t] - sinn[c][t] * d
                    kreuz[c][t] = kreuz[c][t] - sinn[a][t] * d
            E -= (paare[0] + paare[1]) * float(np.sum(d))
    rho = [sum(sinn[s][t] * bahn[s][1] ** 2 * r * r * h for s in schalen) for t in (0, 1)]
    Vq = [np.where(rho[t] > 1e-300, 2 * quer[t] / np.maximum(rho[t], 1e-300), 0.0) for t in (0, 1)]
    V = {}
    for a in schalen:
        qa = np.maximum(bahn[a][1] ** 2 * r * r * h, 1e-300)
        if KREUZ_JE_WINDUNG:                                 # die Moden, die sie selbst in den anderen Schalen sperrt
            eta = 1e-3 * float(np.max(qa))
            Vk = [kreuz[a][t] * qa / (qa * qa + eta * eta) for t in (0, 1)]
        else:
            Vk = Vq
        V[a] = sum(sinn[a][t] * (Vk[t] - max(sinn[a][t] - 1, 0) * eigen[a] / qa) for t in (0, 1)) / besetzt[a]
    return E, V


def _gesperrt(besetzt, bahn, V, e, r, h):
    """Bindung sperrt Moden (𝒫, L_eff = 𝒫L): Die gebundenen inneren Bahnen gleicher Windung sind für jede äußere gesperrt; sie
    schließt als tiefste Schwingung in den übrigen Moden (Σφ_aφ_b r²h = 0), in ihrer Neigung V[s] und im Takt ihres Tempos bei
    e[s]. Von innen nach außen; inverse Iteration im Raum ohne die gesperrten Moden. Gibt die Bahnen {s: (ε, φ)}."""
    from scipy.linalg import solve_banded
    aus = {}
    B = r * r * h
    for s in sorted(besetzt, key=lambda x: (x[1], x[0])):
        innen = [aus[t][1] for t in sorted(aus) if t[1] == s[1] and t[0] < s[0]]
        d, n, _ = _matrix(V[s], _takt(V[s], e[s]), r, h, s[1], -1)
        mal = lambda x: d * x + np.append(n * x[1:], 0) + np.append(0, n * x[:-1])
        x = bahn[s][1].copy()
        for t in innen:
            x = x - float(np.sum(x * t * B)) * t
        x = x / np.sqrt(float(np.sum(x * x * B)))
        E = float(x @ mal(x))
        if innen:
            U = np.array([t * B for t in innen]).T
            for _ in range(60):
                ab = np.vstack([np.append(0, n), d - (E - 1e-9 * max(1.0, abs(E))) * B, np.append(n, 0)])
                Y = solve_banded((1, 1), ab, np.column_stack([B * x, U]), check_finite=False)
                y, Yu = Y[:, 0], Y[:, 1:]
                y = y - Yu @ np.linalg.solve(U.T @ Yu, U.T @ y)
                x = y / np.sqrt(float(np.sum(y * y * B)))
                En = float(x @ mal(x))
                if abs(En - E) < 1e-12 * max(1.0, abs(En)):
                    E = En
                    break
                E = En
        aus[s] = (E, x)
    return aus


def _kreis(Z, besetzt, r, h, start=None, runden=300):
    """Kreis der gebundenen Hülle: jedes Elektron in der Neigung des Kerns, der Windungen der übrigen und der Sperre gleichsinniger
    Windungen (_sperre). Gibt (Bahnen, VW, Faltung der ganzen Hülle E [Hartree])."""
    Vk = _kern(Z, r)
    q = lambda b: b[1] ** 2 * r * r * h
    bahn = {s: (start or {}).get(s) or _bahn(Z, Vk, r, h, s) for s in besetzt}
    eins = {s: _neigung(q(bahn[s]), r) for s in besetzt}
    VW = sum(m * eins[s] for s, m in besetzt.items())
    Vs = _sperre(besetzt, bahn, r, h)[1]
    alt = None
    for _ in range(runden):
        Vb = {s: Vk + VW - eins[s] + Vs[s] for s in besetzt}
        e = {s: bahn[s][0] for s in besetzt}
        bahn = _gesperrt(besetzt, {s: _bahn(Z, Vb[s], r, h, s, e[s]) for s in besetzt}, Vb, e, r, h)
        eins = {s: _neigung(q(bahn[s]), r) for s in besetzt}
        VW_neu = sum(m * eins[s] for s, m in besetzt.items())
        steht = float(np.max(np.abs(VW_neu - VW) * r))           # der Kreis steht erst, wenn auch die Neigung der Windungen steht
        VW = 0.5 * VW + 0.5 * VW_neu
        neu = _sperre(besetzt, bahn, r, h)[1]
        Vs = {s: 0.5 * Vs[s] + 0.5 * neu[s] for s in besetzt}
        jetzt = np.array([bahn[s][0] for s in sorted(besetzt)])
        if alt is not None and np.max(np.abs(jetzt - alt)) < 1e-9 and steht < 1e-7:
            break
        alt = jetzt
    Es, Vs = _sperre(besetzt, bahn, r, h)
    E = Es + sum(m * (bahn[s][0] - float(np.sum((0.5 * (VW - eins[s]) + Vs[s]) * q(bahn[s])))) for s, m in besetzt.items())
    if PAARUNG:
        E += _paarung(Z, besetzt, bahn, {s: Vk + VW - eins[s] + Vs[s] for s in besetzt}, r, h)
    if GRUNDTERM:
        E += _grundterm(besetzt, bahn, r, h)
    if VERWINDUNG:                                           # ζ je Windung aus dem Gefälle ihres Takts
        z = {}
        for s in besetzt:
            K = _takt(Vk + VW - eins[s] + Vs[s], bahn[s][0])
            z[s] = (bahn[s][0], bahn[s][1], float(np.sum(bahn[s][1] ** 2 * np.gradient(K, np.log(r)) * h)))
        E += _verwindung(besetzt, z)
    return bahn, VW, E


def _oertlich(besetzt, bahn, r, h):
    """Örtliche Zählung derselben Bahnen [Hartree]. Gibt (örtlich gegen die Sperre je Umlaufsinn, Kopplung gegensinniger
    Windungen): Örtlich wird jede Windung in der ganzen Neigung gerechnet und Sperre und Kopplung aus der Dichte jedes Umlaufsinns am
    Ort gelesen, ½Σ m·⟨eigene Neigung⟩ − Sperre + E_ort[ρ↑, ρ↓]; die Kopplung, die die Sperre gleichsinniger Windungen nicht
    zählt (vor allem gegensinniger), wird in beiden Zählweisen am Ort gelesen, E_ort − E_Sperre-örtlich."""
    from huelle_1_loeser import _xc
    from huelle_5_fassungen import _xc_korrelation
    sinn = _sinne(besetzt)
    q = {s: bahn[s][1] ** 2 * r * r * h for s in besetzt}
    selbst = 0.5 * sum(m * float(np.sum(_neigung(q[s], r) * q[s])) for s, m in besetzt.items())
    rho = [sum(sinn[s][t] * q[s] for s in besetzt) / (4 * math.pi * r ** 3 * h) for t in (0, 1)]
    dV = (rho[0] + rho[1]) * 4 * math.pi * r ** 3 * h
    exc, ec = _xc(rho[0], rho[1])[0], _xc_korrelation(rho[0], rho[1])[0]
    return selbst - _sperre(besetzt, bahn, r, h)[0] + float(np.sum(exc * dV)), float(np.sum(ec * dV))


def _zustand(Z, besetzt, r, h, start=None):
    """Kreis der Konfiguration mit beiden Zählweisen: dict(bes, bahn, VW, je = Faltung je Umlaufsinn mit der örtlich gelesenen
    Kopplung, ort = örtlich gezählt, tempo = δ = 2|ε|/c² jeder gebundenen Unterschale im Feld der Konfiguration, besetzt oder
    frei) [Hartree]."""
    b, VW, E = _kreis(Z, besetzt, r, h, start)
    tempo = {f"{n}{l}": 2 * abs(e) / (c * c) for (n, l), e in _frei(Z, besetzt, VW, r, h).items()}
    tempo.update({f"{n}{l}": 2 * abs(b[(n, l)][0]) / (c * c) for (n, l) in besetzt if b[(n, l)][0] < 0})
    ort, kopplung = _oertlich(besetzt, b, r, h)
    return dict(bes=besetzt, bahn=b, VW=VW, je=E + kopplung, ort=E + ort, tempo=tempo)


def _abstand_zwei(Z, A, C):
    """E(C) − E(A) nach den zwei Zählweisen, um die Mitte gespiegelt: (1 − w̄)·Δ_örtlich + w̄·Δ_je Umlaufsinn, w̄ das Mittel der
    Gewichte beider Blickrichtungen (Füllung der Spinhälfte mal Verhältnis der Tempi, je im Feld der Ausgangskonfiguration)."""
    from huelle_5_fassungen import _gewicht
    w = 0.5 * (_gewicht(Z, A["bes"], C["bes"], tempo=A["tempo"]) + _gewicht(Z, C["bes"], A["bes"], tempo=C["tempo"]))
    return (1 - w) * (C["ort"] - A["ort"]) + w * (C["je"] - A["je"])


def _hund(besetzt):
    """Windungen offener Unterschalen im tiefsten Zustand (Spiegel der Schale um ihre Mitte, f(x) + f(−x) = 1): alle im selben
    Umlaufsinn zuerst, dann die Moden m = l, l − 1, …; S am größten in der Mitte der Schale, L(lⁿ) = L(l^{2(2l+1)−n}).
    Gibt [(Unterschale, m, Sinn)]."""
    sinn = _sinne(besetzt)
    offen = [s for s, m in besetzt.items() if 0 < m < 2 * (2 * s[1] + 1)]
    return [(s, m, t) for s in offen for t in (0, 1) for m in range(s[1], s[1] - sinn[s][t], -1)]


def _grundterm(besetzt, bahn, r, h):
    """Überlagerung der Windungen offener Unterschalen in den Moden des tiefsten Zustands (_hund) gegen das Mittel über die Moden,
    das Kreis und Sperre schon zählen [Hartree]; geschlossene Schalen sind kugelrund und tragen nichts bei."""
    el = _hund(besetzt)
    offen = sorted({s for s, _, _ in el})
    q = {s: bahn[s][1] ** 2 * r * r * h for s in offen}
    F, G = {}, {}
    for i, a in enumerate(offen):
        for c in offen[i:]:
            for k in range(0, 2 * min(a[1], c[1]) + 1, 2):
                F[(a, c, k)] = F[(c, a, k)] = float(np.sum(q[a] * _yk(q[c], r, k)))
            if a != c:
                qa = bahn[a][1] * bahn[c][1] * r * r * h
                for k in range(abs(a[1] - c[1]), a[1] + c[1] + 1, 2):
                    G[(a, c, k)] = G[(c, a, k)] = float(np.sum(qa * _yk(qa, r, k)))
    E = 0.0
    for i, (a, ma, ta) in enumerate(el):
        for c, mc, tc in el[i + 1:]:
            la, lc = a[1], c[1]
            E += sum(_gaunt(k, la, ma, la, ma) * _gaunt(k, lc, mc, lc, mc) * F[(a, c, k)]
                     for k in range(2, 2 * min(la, lc) + 1, 2))
            if ta != tc:
                continue
            if a == c:
                E -= sum(_gaunt(k, la, ma, la, mc) ** 2 * F[(a, a, k)] + _moden(la)[k] * F[(a, a, k)]
                         for k in range(2, 2 * la + 1, 2))
            else:
                E -= sum((_gaunt(k, la, ma, lc, mc) ** 2 - _kopplung(la, k, lc)) * G[(a, c, k)]
                         for k in range(abs(la - lc), la + lc + 1, 2))
    return E


def _verwindung(besetzt, bahn):
    """Verwindung des tiefsten Niveaus [Hartree]: Jede Windung koppelt ihren Sinn an ihren Umlauf, ζ·m·m_s (ζ aus dem Gefälle ihres
    Takts); im tiefsten Zustand (_hund) mit L = Σm, S = Σm_s ist A = Σζ·m·m_s/(L·S), und das tiefste J ist |L − S| für A > 0,
    L + S für A < 0."""
    el = _hund(besetzt)
    L = sum(m for _, m, _ in el)
    S = sum(0.5 - t for _, _, t in el)
    X = sum(bahn[s][2] * m * (0.5 - t) for s, m, t in el)
    if L == 0 or S == 0:
        return 0.0
    A = X / (L * S)
    if A < 0:
        return A * L * S
    J = abs(L - S)
    return A / 2 * (J * (J + 1) - L * (L + 1) - S * (S + 1))


def _moden_frei(V, r, h, l, i0, anzahl):
    """Moden der Windung l in der Neigung V von der Stufe i0 an (Knotenzahl), im geschlossenen Raum des Gitters: (ε, φ) mit
    Σφ²r²h = 1, Takt aus dem Tempo an der Schwelle."""
    w, v = _eigen(V, _takt(V, 0.0), r, h, l, -1, i0 + anzahl, -1e3)
    return np.asarray(w[i0:]), v[:, i0:].T


def _paarung(Z, besetzt, bahn, Vb, r, h):
    """Einrasten (06): Eine Mode in beiden Umlaufsinnen hat gleiche Frequenz, sie rastet ein und wird ein Zyklus (T² → T¹), das
    Paar; gegensinnige Windungen verschiedener Moden sind kein Paar. Gezählt werden die Moden in beiden Umlaufsinnen je
    Unterschale. Die Kopplung ist der Schwankungsteil der gegenseitigen Neigung beider Windungen der Mode m, ⟨cd|V|mm⟩ über die
    freien Moden c, d (Kugelfunktionen mit m_c + m_d = 2m, gemittelt über m); Rauschen und Dissipation stammen aus derselben
    Dämpfung. Je Paar freier Moden ist die Bindung die des Zweierkreises,
    (√(Δ² + 4|V|²) − Δ)/2 mit Δ = ε_c + ε_d − ε_a − ε_b: weit außerhalb des Fensters |V|²/Δ, eingerastet |V|; Moden unter dem
    Paar sind nicht frei. Gezählt in den äußeren Unterschalen (ns, np, (n − 1)d, (n − 2)f),
    die zwischen den Konfigurationen wechseln. Freie Moden sind die Schwingungen im Feld aller Gebundenen, gesperrt gegen die
    gebundenen Bahnen gleicher Windung (Bindung sperrt Moden), dazu die freien Plätze offener Unterschalen. Gibt die Bindung aller eingerasteten Paare [Hartree, ≤ 0]."""
    sinn = _sinne(besetzt)
    n = max(k[0] + max(k[1] - 1, 0) for k in besetzt if k[1] <= 3)
    aussen = [y for y in [(n, 0), (n, 1), (n - 1, 2), (n - 2, 3)] if besetzt.get(y)]
    B = r * r * h
    voll = lambda y: 2 * (2 * y[1] + 1)
    L_PAAR = globals()["L_PAAR"] or max(4, 2 * max(y[1] for y in aussen))   # zwei Windungen erreichen den Grad 2l
    Vn = _kern(Z, r) + sum(m * _neigung(bahn[y][1] ** 2 * B, r) for y, m in besetzt.items())
    moden = {}
    for l in range(L_PAAR + 1):                               # freie Moden: Feld aller Gebundenen, gesperrt gegen die gebundenen
        geb = [bahn[y][1] for y in sorted(besetzt) if y[1] == l]
        w, P = _moden_frei(Vn, r, h, l, 0, N_FREI + len(geb))
        for g in geb:
            P = P - np.outer(P @ (g * B), g)
        S = (P * B) @ P.T                                    # im Raum ohne die gebundenen Moden neu schwingen lassen
        d, e, _ = _matrix(Vn, _takt(Vn, 0.0), r, h, l, -1)
        HP = P * d + np.hstack([P[:, 1:] * e, np.zeros((len(P), 1))]) + np.hstack([np.zeros((len(P), 1)), P[:, :-1] * e])
        Hm = HP @ P.T
        from scipy.linalg import eigh
        s_w, s_v = np.linalg.eigh(S)
        halten = s_w > 1e-2                                  # nur Richtungen, die nach dem Sperren merklich bleiben
        X = s_v[:, halten] / np.sqrt(s_w[halten])
        ew, ev = eigh(X.T @ Hm @ X)
        moden[l] = (ew[:N_FREI], (ev[:, :N_FREI].T @ X.T) @ P)
    frei = {}
    for a in aussen:
        frei[a] = {}
        for l in range(L_PAAR + 1):
            offen = [y for y in sorted(besetzt) if y[1] == l and besetzt[y] < voll(y) and y != a]
            w = np.concatenate([[bahn[y][0] for y in offen], moden[l][0]])
            P = np.vstack([bahn[y][1] for y in offen] + [moden[l][1]]) if offen else moden[l][1]
            platz = np.array([1 - besetzt[y] / voll(y) for y in offen] + [1.0] * len(moden[l][0]))
            frei[a][l] = (w, P, platz)
    E = 0.0
    for i, a in enumerate(aussen):
        for b in aussen[i:]:
            if EINRASTEN == "mode":
                zahl = sinn[a][1] if a == b else 0          # nur Moden in beiden Umlaufsinnen
            else:                                            # Kopplungsfenster: alle gegensinnigen Windungen im Band
                zahl = sinn[a][0] * sinn[a][1] if a == b else sinn[a][0] * sinn[b][1] + sinn[a][1] * sinn[b][0]
            if zahl == 0:
                continue
            la, lb, ea, eb = a[1], b[1], bahn[a][0], bahn[b][0]
            e2 = k2 = 0.0
            for lc in range(L_PAAR + 1):
                Ec, Pc, fc = frei[a][lc]
                Qc = Pc * bahn[a][1] * B
                for ld in range(L_PAAR + 1):
                    Ed, Pd, fd = frei[b][ld]
                    Qd = Pd * bahn[b][1] * B
                    ks = [k for k in range(max(abs(la - lc), abs(lb - ld)), min(la + lc, lb + ld) + 1)
                          if not (la + k + lc) % 2 and not (lb + k + ld) % 2]
                    if not ks:
                        continue
                    R = {k: Qc @ np.array([_yk(x, r, k) for x in Qd]).T for k in ks}
                    D = Ec[:, None] + Ed[None, :] - ea - eb
                    F = fc[:, None] * fd[None, :] * (D > 0)             # Moden unter dem Paar sind nicht frei
                    paare = [(m, m) for m in range(-la, la + 1)] if (EINRASTEN == "mode" and a == b) else \
                        [(ma, mb) for ma in range(-la, la + 1) for mb in range(-lb, lb + 1)]
                    for ma, mb in paare:                     # die beiden Windungen in den Moden ma, mb
                        for mc in range(-lc, lc + 1):
                            md = ma + mb - mc
                            if abs(md) > ld:
                                continue
                            M = sum(_gaunt(k, lc, mc, la, ma) * _gaunt(k, ld, md, lb, mb) * R[k] for k in ks)
                            W = M * M * F
                            e2 += float(np.sum(np.where(D > 0, (np.sqrt(D * D + 4 * W) - D) / 2, 0.0))) / len(paare)
                            k2 += float(np.sum(W)) / len(paare)
            if a != b:
                _paarung.fenster.append((a, b, abs(ea - eb) * HARTREE_EV, 2 * np.sqrt(k2) * HARTREE_EV, zahl * e2 * HARTREE_EV))
            if a != b and abs(ea - eb) > 2 * np.sqrt(k2):
                continue                                     # außerhalb des Fensters: frei, kein Einrasten
            E -= zahl * e2
    return E


_paarung.fenster = []                                        # Protokoll der Fenster: (a, b, |Δε|, 2k, Bindung) [eV]


def paar_probe(Z, besetzt, n_frei=(100, 200, 400), l_paar=(4, 6, 8), einrasten="mode"):
    """Probe des Einrastens an einer Konfiguration: Bindung aller eingerasteten Paare [eV] je Zahl der freien Moden und höchster
    Windung der freien Moden (None: 2·l der äußeren Schalen); zeigt, ob die Summe über die freien Moden steht.
    Gibt {(n_frei, l_paar): eV}."""
    global N_FREI, L_PAAR, EINRASTEN
    r, h = _gitter()
    EINRASTEN = einrasten
    bahn, VW, E = _kreis(Z, besetzt, r, h)
    eins = {s: _neigung(bahn[s][1] ** 2 * r * r * h, r) for s in besetzt}
    Vs = _sperre(besetzt, bahn, r, h)[1]
    Vb = {s: _kern(Z, r) + VW - eins[s] + Vs[s] for s in besetzt}
    aus = {}
    for nf in n_frei:
        for lp in l_paar:
            N_FREI, L_PAAR = nf, lp
            aus[(nf, lp)] = _paarung(Z, besetzt, bahn, Vb, r, h) * HARTREE_EV
            print(Z, nf, lp, aus[(nf, lp)], flush=True)
    return aus


def _geschlossen(s, m):
    """Die Zählregel schließt: ein Umlaufsinn der Unterschale ganz (m = 2l + 1) oder beide (m = 2(2l + 1)), mit mindestens zwei
    Windungen; einrasten können nur gekoppelte Zyklen (Tᵐ → T¹, m ≥ 2), eine einzelne Windung ist schon ein Zyklus."""
    return m >= 2 and m in (2 * s[1] + 1, 2 * (2 * s[1] + 1))


def _ziele(besetzt, frei):
    """Offene Plätze: jede teilweise besetzte Unterschale und die vier tiefsten leeren Zustände; welche bindet, entscheidet die
    Faltung des Ganzen (im Feld ohne Sperre liegt eine teilweise gefüllte Schale zu hoch, um nach ihrer Lage gewählt zu werden)."""
    leer = [y for y in sorted(frei, key=frei.get) if y not in besetzt][:4]
    return [y for y, m in besetzt.items() if m < 2 * (2 * y[1] + 1)] + leer


def _frei(Z, besetzt, VW, r, h):
    """Freie Zustände im Feld des Kerns und aller Gebundenen: {s: ε} (nur gebundene, ε < 0)."""
    aus = {}
    V = _kern(Z, r) + VW
    for l in range(L_MAX + 1):
        w, _ = _eigen(V, _takt(V, 0.0), r, h, l, -1, N_MAX - l, -2.0 * Z * Z - 5)
        w = [float(_eigen(V, _takt(V, x), r, h, l, -1, i + 1, -2.0 * Z * Z - 5, list(w[:i + 1]))[0][-1])
             for i, x in enumerate(w)]                       # jeder Zustand im Takt seines Tempos
        for i, e in enumerate(w):
            s = (l + 1 + i, l)
            if e < 0 and besetzt.get(s, 0) < 2 * (2 * l + 1):
                aus[s] = float(e)
    return aus


def entstehung(Z):
    """Hülle von Z, Fang für Fang. Gibt (Konfiguration {(n, l): Zahl}, Fänge [(n, l, ε [Hartree], δ beim Fang)], Bindung des
    am schwächsten gebundenen ε [Hartree], Höhe v und Ort R des Spiegelfelds)."""
    import kerne
    import kern as K
    r, h = _gitter()
    f3 = K.polytrope(3.0)[2]
    mc2 = float(kerne.a_tal(Z, kerne.TROPFEN)) * K.m_u * K.c ** 2 / K.e        # eV
    besetzt, faenge, bahn, VW = {}, [], {}, np.zeros_like(r)
    for k in range(Z):
        frei = _frei(Z, besetzt, VW, r, h)
        if not frei:
            break                                            # kein gebundener Platz mehr: die Hülle ist hier zu Ende
        wahl = []
        for s in _ziele(besetzt, frei):             # das Ganze steuert: dort, wo die Faltung der ganzen Hülle am tiefsten wird
            neu = {**besetzt, s: besetzt.get(s, 0) + 1}
            if ZAEHLWEISEN:
                wahl.append((s, _zustand(Z, neu, r, h, bahn)))
                continue
            b, V, E = _kreis(Z, neu, r, h, bahn)
            wahl.append((b[s][0], s, neu, b, V, E))
        if ZAEHLWEISEN:
            s, jetzt = wahl[0]
            for t, X in wahl[1:]:
                if _abstand_zwei(Z, jetzt, X) < 0:
                    s, jetzt = t, X
            e, besetzt, bahn, VW, E = jetzt["bahn"][s][0], jetzt["bes"], jetzt["bahn"], jetzt["VW"], jetzt
        else:
            e, s, besetzt, bahn, VW, E = min(wahl, key=lambda x: x[5])
        faenge.append((s[0], s[1], e, -e * HARTREE_EV / (f3 * mc2)))
        while True:                                          # Rückkehr: geschützt ist die Zählung, nicht die einzelne Schleife
            oben = sorted([x for x in besetzt if not _geschlossen(x, besetzt[x])], key=lambda x: bahn[x][0])[-3:]
            ziele = [(y, 0) for y in _ziele(besetzt, _frei(Z, besetzt, VW, r, h))]
            bester = None
            for x in oben:
                for y, _ in ziele:
                    if y == x:
                        continue
                    neu = {p: m for p, m in {**besetzt, x: besetzt[x] - 1, y: besetzt.get(y, 0) + 1}.items() if m}
                    if ZAEHLWEISEN:
                        X = _zustand(Z, neu, r, h, bahn)
                        if _abstand_zwei(Z, E, X) < -1e-7 and (bester is None or _abstand_zwei(Z, bester[0], X) < 0):
                            bester = (X, neu, X["bahn"], X["VW"])
                        continue
                    b, V, En = _kreis(Z, neu, r, h, bahn)
                    if En < E - 1e-7 and (bester is None or En < bester[0]):
                        bester = (En, neu, b, V)
            if bester is None:
                break
            E, besetzt, bahn, VW = bester
    F = -r * r * (_kern(Z, r) + VW)
    i = int(np.argmax(F))
    return besetzt, faenge, max(bahn[s][0] for s in besetzt), 8 * float(F[i]), float(r[i])


def _ein(Z):
    return Z, entstehung(Z)


def rechne_entstehung(z_bereich=range(1, 103), prozesse=4, schluessel="_entstehung", neu=False, zaehlweisen=False,
                      verwindung=False):
    """entstehung für alle Z; unter schluessel je Element sofort gespeichert, schon gerechnete übersprungen (neu: auch diese,
    jedes Element ersetzt, sobald es fertig ist); zaehlweisen: Wahl nach den zwei Zählweisen; verwindung: mit der Lage des
    tiefsten Niveaus."""
    from multiprocessing import Pool
    global ZAEHLWEISEN, VERWINDUNG
    ZAEHLWEISEN, VERWINDUNG = zaehlweisen, verwindung
    fertig = {} if neu else lesen().get(schluessel, {})
    with Pool(prozesse) as p:
        for Z, (bes, faenge, oben, v, R) in p.imap_unordered(_ein, [Z for Z in z_bereich if str(Z) not in fertig]):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = dict(
                konf=[[n, l, m] for (n, l), m in sorted(bes.items())], faenge=[list(x) for x in faenge],
                bindung_ev=-oben * HARTREE_EV, v=v, R=R)
            schreiben(alle)


def entstehung_ergebnis(schluessel="_entstehung"):
    """Prüfung am Ende gegen die gemessenen Konfigurationen (PubChem 1–102): Treffer und Abweichungen."""
    import kerne
    gem = kerne.belegung_gemessen()
    e = lesen().get(schluessel, {})
    aus = {}
    for Z, d in e.items():
        b = {(n, l): m for n, l, m in d["konf"]}
        g = {k: v for k, v in gem[int(Z)].items() if v}
        aus[int(Z)] = b == g
    return dict(treffer=sum(aus.values()), zahl=len(aus), daneben=sorted(Z for Z, t in aus.items() if not t))


def pruefe_zwei(Z, schluessel="_entstehung"):
    """Prüfung am Ende: gemessene Konfiguration (PubChem) gegen die der Entstehung unter schluessel, beide als Kreis, nach den zwei
    Zählweisen (_abstand_zwei) und je Umlaufsinn allein [eV]; > 0: die gemessene liegt tiefer."""
    import kerne
    r, h = _gitter()
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    b = {(n, l): m for n, l, m in lesen()[schluessel][str(Z)]["konf"]}
    if b == g:
        return None
    G, B = _zustand(Z, g, r, h), _zustand(Z, b, r, h)
    return _abstand_zwei(Z, G, B) * HARTREE_EV, (B["je"] - G["je"]) * HARTREE_EV, (B["ort"] - G["ort"]) * HARTREE_EV


def _pz_ein(Z):
    return Z, pruefe_zwei(Z)


def rechne_pruefe_zwei(z_bereich, prozesse=4, schluessel="_pruefung_kreis", verwindung=False):
    """pruefe_zwei für die Elemente z_bereich, je Element sofort unter schluessel gespeichert; verwindung: mit der Lage des
    tiefsten Niveaus."""
    from multiprocessing import Pool
    global VERWINDUNG
    VERWINDUNG = verwindung
    with Pool(prozesse) as p:
        for Z, d in p.imap_unordered(_pz_ein, z_bereich):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = d
            schreiben(alle)
            print(Z, d, flush=True)


def _paare(besetzt):
    """Gegensinnige Paare in derselben Mode: je Unterschale die Windungen über einem ganzen Umlaufsinn."""
    return sum(max(0, m - (2 * s[1] + 1)) for s, m in besetzt.items())


def _abstand(Z):
    """Faltung der gemessenen Konfiguration über der gerechneten (eV) und der Unterschied ihrer Paare."""
    import kerne
    r, h = _gitter()
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    b = {(n, l): m for n, l, m in lesen()["_entstehung"][str(Z)]["konf"]}
    hoehe = lambda VW: 8 * float(np.max(-r * r * (_kern(Z, r) + VW)))   # Höhe v des Spiegelfelds der Konfiguration
    _, Vb, Eb = _kreis(Z, b, r, h)
    _, Vg, Eg = _kreis(Z, g, r, h)
    return Z, dict(ev=(Eg - Eb) * HARTREE_EV, paare=_paare(g) - _paare(b), v_gemessen=hoehe(Vg), v_gerechnet=hoehe(Vb))


def _rand(Z):
    """Rand der gemessenen Konfiguration: jede Nachbarkonfiguration, in der ein Elektron zwischen den äußeren Schalen (ns, np,
    (n − 1)d, (n − 2)f) umzieht, mit ihrer Faltung über der gemessenen (eV) und dem Unterschied ihrer Paare."""
    import kerne
    r, h = _gitter()
    g = {k: v for k, v in kerne.belegung_gemessen()[Z].items() if v}
    n = max(k[0] + max(k[1] - 1, 0) for k in g if k[1] <= 3)
    aussen = [y for y in [(n, 0), (n, 1), (n - 1, 2), (n - 2, 3)] if y[0] > y[1]]
    bahn, _, Eg = _kreis(Z, g, r, h)
    aus = {}
    for x in [x for x in aussen if g.get(x)]:
        for y in [y for y in aussen if y != x and g.get(y, 0) < 2 * (2 * y[1] + 1)]:
            neu = {p: m for p, m in {**g, x: g[x] - 1, y: g.get(y, 0) + 1}.items() if m}
            En = _kreis(Z, neu, r, h, bahn)[2]
            aus[f"{x[0]}{'spdfg'[x[1]]}→{y[0]}{'spdfg'[y[1]]}"] = dict(ev=(En - Eg) * HARTREE_EV, paare=_paare(neu) - _paare(g))
    return Z, aus


def rechne_rand(z_bereich, prozesse=4, schluessel="_rand", paarung=False, einrasten="mode", grundterm=False, verwindung=False):
    """_rand für die Elemente z_bereich; je Element sofort gespeichert."""
    from multiprocessing import Pool
    global PAARUNG, EINRASTEN, GRUNDTERM, VERWINDUNG
    PAARUNG, EINRASTEN, GRUNDTERM, VERWINDUNG = paarung, einrasten, grundterm, verwindung
    with Pool(prozesse) as p:
        for Z, d in p.imap_unordered(_rand, z_bereich):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = d
            schreiben(alle)


def rechne_abstand(z_bereich=None, prozesse=4, schluessel="_abstand", je_windung=False):
    """Für jedes Element, dessen Hülle anders entsteht als gemessen: wie weit die gemessene Konfiguration in der Rechnung über der
    gerechneten liegt; je Element sofort gespeichert."""
    from multiprocessing import Pool
    global KREUZ_JE_WINDUNG
    KREUZ_JE_WINDUNG = je_windung
    z = z_bereich or entstehung_ergebnis()["daneben"]
    with Pool(prozesse) as p:
        for Z, d in p.imap_unordered(_abstand, z):
            alle = lesen()
            alle.setdefault(schluessel, {})[str(Z)] = d
            schreiben(alle)


def stern(M=None, z_max=102, schluessel="_stern"):
    """Kette der Elemente im Stern: Ein Stern der Masse M (ohne Angabe die Sonnenmasse) besteht aus Element Z und verschmilzt es mit
    dem Proton zu Z + 1, wenn seine Tiefe die Zündung erreicht: kT = f·δ·m̄c² (Polytrope n = 3, m̄ = m_A/(Z + 1)), gezündet bei
    gleichem Tunnelfaktor e^{−√(E_G/E)} wie der Wasserstoff in der Tiefe der Sonne (E_G = 2μc²(παZ)²). Mit jedem schwereren
    Element wächst die Tiefe δ_Z, der Stern zieht sich zusammen (R = GM/(δc²)), Dichte und Druck steigen (P = ½·ρ·c²·δ im Zentrum).
    Hält die Zählregel der Elektronen den Trichter, endet die Kette an der Decke kT_max, außer über der Grenzmasse. Über x = ¼ führt
    nur der Einfang von Neutronen, im von den Nukleonen gehaltenen Trichter: δ = λ̄_p/λ_H, Zelle λ_H·A^{1/3}.
    Danach fällt δ (Ausdehnung des Sterns, ähnlich, R ∝ 1/δ, Längen ∝ 1/δ), und die Hülle fängt ihre Elektronen, jedes bei der Tiefe,
    bei der die Wärme unter seine Bindung fällt, δ_k = ε_k/(f·m̄_k·c²) (Fänge der Etappe entstehung); dort hat das Atom die Zelle
    R_A(δ_k) = R_A(δ_Z)·δ_Z/δ_k. Gibt je Z die Tiefe und den Zustand der Entstehung und das kleinste Verhältnis der Zelle zur Bahn
    beim Fang, R_A/⟨r⟩ (⟨r⟩ = (3n² − l(l + 1))/(2n·√(2ε)) [a_B], Bahn mit der Ladung, die ihre Bindung gibt)."""
    import kerne
    import kern as K
    x1, dt1, f = K.polytrope(3.0)
    dichte_c = x1 / (3 * dt1)                                         # Dichte im Zentrum über dem Mittel, n = 3
    M = M or K.M_sun
    aB = K.bohr_radius()
    zeilen = {int(z["AtomicNumber"]): float(z["AtomicMass"]) for z in kerne._pubchem_zeilen()}
    mc2 = lambda Z: zeilen[Z] * K.m_u * K.c ** 2
    E_G = lambda Z: 2 * (K.m_p * zeilen[Z] * K.m_u / (K.m_p + zeilen[Z] * K.m_u)) * K.c ** 2 * (math.pi * K.alpha * Z) ** 2
    d_sonne = K.GM_sun / (6.957e8 * K.c ** 2)
    kT_1 = f * d_sonne * mc2(1) / 2                                   # Wasserstoff zündet in der Tiefe der Sonne
    lam = kerne.halterlaenge(kerne.a_C(kerne.TROPFEN)) * 1e-15
    z_leiter = kerne.ende_brennleiter(kerne.TROPFEN)[0]
    ent = lesen().get("_entstehung", {})
    aus = {}
    for Z in range(1, z_max + 1):
        m_A = zeilen[Z] * K.m_u
        if Z <= z_leiter:
            kT = kT_1 * E_G(Z) / E_G(1)
            d = kT / (f * mc2(Z) / (Z + 1))
            R = K.G * M / (d * K.c ** 2)
            rho = dichte_c * 3 * M / (4 * math.pi * R ** 3)
            zelle = (3 * m_A / (4 * math.pi * rho)) ** (1 / 3)
            mu_e = zeilen[Z] / Z
            decke = K.decke_entartet((zeilen[Z] / (Z + 1)) / mu_e)[0]
            ort = "Brennleiter" if (kT <= decke or M > K.grenzmasse(mu_e)) else "nicht erreicht (Decke der Zählregel)"
        else:
            d = K.hbar / (K.m_p * K.c) / lam
            zelle = lam * zeilen[Z] ** (1 / 3)
            rho = m_A / (4 * math.pi * zelle ** 3 / 3)
            kT = f * d * mc2(Z) / (Z + 1)
            ort = "Neutroneneinfang"
        P = 0.5 * rho * K.c ** 2 * d
        verh = None
        if str(Z) in ent:
            k = 0
            for n, l, e, _ in ent[str(Z)]["faenge"]:
                eps = -e                                              # Hartree
                dk = eps * HARTREE_EV * K.e / (f * mc2(Z) / (1 + Z - k))
                r = (3 * n * n - l * (l + 1)) / (2 * n * math.sqrt(2 * eps))
                v = zelle * d / dk / aB / r
                verh = v if verh is None else min(verh, v)
                k += 1
        aus[str(Z)] = dict(ort=ort, delta=d, kT_keV=kT / K.e / 1e3, rho_g_cm3=rho / 1e3, P_Pa=P, zelle_aB=zelle / aB, zelle_zu_bahn=verh)
    alle = lesen()
    alle[schluessel] = aus
    schreiben(alle)
    return aus
