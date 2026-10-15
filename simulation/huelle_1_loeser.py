"""Etappe 1: der Löser der Hülle. Gitter in log r, ZORA-Radialgleichung in φ = u/√r, Austausch und Korrelation
(LSDA nach Perdew-Zunger, PBE), Hartree-Feld, selbstkonsistente Hülle mit Spin nach Hund und Spin-Bahn-Glied (erster
Ordnung oder gekoppelt je j), Ausrichtung der Windungen. Kern mit R = λ_H·A^{1/3} im Tal; atomare Einheiten."""
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import kern as K


c = 1 / K.alpha                                          # Lichtgeschwindigkeit in atomaren Einheiten

FM = 1e-15 / (K.hbar / (K.m_e * K.c * K.alpha))        # 1 fm in Bohr-Radien

def _gitter(h=0.016, rmin=1e-6, rmax=60.0):
    x = np.arange(math.log(rmin), math.log(rmax) + h / 2, h)
    return np.exp(x), h

def _matrix(V, Kf, r, h, l, kappa):
    """ZORA-Radialgleichung in φ = u/√r auf logarithmischem Gitter: Tridiagonale A (Diagonale, Nebendiagonale) und B."""
    Kx = np.gradient(Kf, np.log(r))
    cm = 0.5 * (Kf[:-1] + Kf[1:]) / h
    W = Kf / 8 - Kx / 4 + Kf * l * (l + 1) / 2 - kappa * Kx / 2 + V * r * r
    d = np.zeros_like(r)
    d[:-1] += cm / 2
    d[1:] += cm / 2
    return d + W * h, -cm / 2, r * r * h

def _sturm(d, e, B, E):
    """Zahl der Eigenwerte von (A, B) unter jedem E (Trägheit von A − E·B, Sturm-Folge)."""
    p = d[0] - E * B[0]
    zahl = (p < 0).astype(int)
    e2 = e * e
    for i in range(1, len(d)):
        p = d[i] - E * B[i] - e2[i - 1] / np.where(p == 0, 1e-300, p)
        zahl += p < 0
    return zahl

def _eigen(V, Kf, r, h, l, kappa, anzahl, sigma, schaetzung=None, rueckfall=False):
    """Niedrigste Eigenwerte und φ (Σφ²r²h = 1). Mit Schätzwerten inverse Iteration je Niveau (Knotenzahl geprüft),
    sonst Lanczos auf dem ganzen Spektrum."""
    from scipy.linalg import solve_banded
    d, e, B = _matrix(V, Kf, r, h, l, kappa)
    if schaetzung is not None and schaetzung != "lanczos":
        w, vs = [], []
        for i, E in enumerate(schaetzung[:anzahl]):
            x = np.ones_like(r)
            for _ in range(40):
                ab = np.vstack([np.append(0, e), d - (E - 1e-9 * max(1.0, abs(E))) * B, np.append(e, 0)])
                y = solve_banded((1, 1), ab, B * x, check_finite=False)
                x = y / math.sqrt(np.sum(y * y * B))
                Ax = d * x
                Ax[:-1] += e * x[1:]
                Ax[1:] += e * x[:-1]
                En = float(x @ Ax)
                if abs(En - E) < 1e-11 * max(1.0, abs(En)):
                    E = En
                    break
                E = En
            w.append(E)
            vs.append(x)
        w = np.array(w)
        dE = 1e-7 * np.maximum(1.0, np.abs(w))
        if not np.array_equal(_sturm(d, e, B, np.concatenate([w - dE, w + dE])),
                              np.concatenate([np.arange(anzahl), np.arange(1, anzahl + 1)])):
            return _eigen(V, Kf, r, h, l, kappa, anzahl, sigma, None if not rueckfall else "lanczos",
                          rueckfall=rueckfall)                # Startwerte neu, dann Lanczos
        return w, np.array(vs).T
    if schaetzung is None:                                   # Startwerte: Tridiagonalform B^(-1/2)·A·B^(-1/2)
        from scipy.linalg import eigh_tridiagonal
        start = eigh_tridiagonal(d / B, e / np.sqrt(B[:-1] * B[1:]), eigvals_only=True, select="i",
                                 select_range=(0, anzahl - 1), lapack_driver="stemr")
        return _eigen(V, Kf, r, h, l, kappa, anzahl, sigma, list(start), rueckfall=True)
    try:
        M = sp.diags([e, d - sigma * B, e], [-1, 0, 1], format="csc")
        mu, v = eigsh(sp.diags(B, format="csc"), k=anzahl, M=M, which="LA", tol=1e-12, maxiter=20000)
        o = np.argsort(-mu)
        w, v = sigma + 1 / mu[o], v[:, o]
    except Exception:
        from scipy.linalg import eigh
        w, v = eigh(np.diag(d) + np.diag(e, 1) + np.diag(e, -1), np.diag(B), subset_by_index=(0, anzahl - 1))
    return w, v / np.sqrt(np.sum(v * v * B[:, None], axis=0))

def _xc(nu, nd):
    """LSDA: Slater-Austausch, Korrelation nach Perdew-Zunger; (ε_xc je Elektron, v↑, v↓)."""
    nu, nd = np.maximum(nu, 0), np.maximum(nd, 0)
    n = np.maximum(nu + nd, 1e-30)
    z = np.clip((nu - nd) / n, -1, 1)
    ex = -0.75 * (6 / math.pi) ** (1 / 3) * (nu ** (4 / 3) + nd ** (4 / 3)) / n
    rs = (3 / (4 * math.pi * n)) ** (1 / 3)

    def pz(g, b1, b2, A, B, Cc, D):
        sr, lr = np.sqrt(rs), np.log(rs)
        e = np.where(rs >= 1, g / (1 + b1 * sr + b2 * rs), A * lr + B + Cc * rs * lr + D * rs)
        de = np.where(rs >= 1, -g * (b1 / (2 * sr) + b2) / (1 + b1 * sr + b2 * rs) ** 2, A / rs + Cc * (lr + 1) + D)
        return e, de
    eu, deu = pz(-0.1423, 1.0529, 0.3334, 0.0311, -0.048, 0.0020, -0.0116)
    ep, dep = pz(-0.0843, 1.3981, 0.2611, 0.01555, -0.0269, 0.0007, -0.0048)
    f = ((1 + z) ** (4 / 3) + (1 - z) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    df = (4 / 3) * ((1 + z) ** (1 / 3) - (1 - z) ** (1 / 3)) / (2 ** (4 / 3) - 2)
    ec, dr, dz = eu + f * (ep - eu), deu + f * (dep - deu), df * (ep - eu)
    vc = ec - rs / 3 * dr
    return ex + ec, -(6 / math.pi * nu) ** (1 / 3) + vc + dz * (1 - z), -(6 / math.pi * nd) ** (1 / 3) + vc - dz * (1 + z)

def _xc_korrelation(nu, nd):
    """Nur die Korrelation der LSDA (Perdew-Zunger): (ε_c je Elektron, v_c↑, v_c↓)."""
    nu, nd = np.maximum(nu, 0), np.maximum(nd, 0)
    n = np.maximum(nu + nd, 1e-30)
    exc, vu, vd = _xc(nu, nd)
    ex = -0.75 * (6 / math.pi) ** (1 / 3) * (nu ** (4 / 3) + nd ** (4 / 3)) / n
    return exc - ex, vu + (6 / math.pi * nu) ** (1 / 3), vd + (6 / math.pi * nd) ** (1 / 3)

def _pw92(rs, z):
    """Korrelation des gleichförmigen Mediums nach Perdew-Wang (ε_c je Elektron)."""
    def G(A, a1, b1, b2, b3, b4):
        s = np.sqrt(rs)
        return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / (2 * A * (b1 * s + b2 * rs + b3 * rs * s + b4 * rs * rs)))
    e0 = G(0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
    e1 = G(0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    ac = -G(0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
    f = ((1 + z) ** (4 / 3) + (1 - z) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    return e0 + ac * f / 1.709921 * (1 - z ** 4) + (e1 - e0) * f * z ** 4

def _pbe(nu, nd, gu, gd):
    """Austausch-Korrelation je Volumen nach Perdew-Burke-Ernzerhof aus den Spindichten und ihren Steigungen dρ/dr."""
    e = 0.0
    for n1, g1 in ((nu, gu), (nd, gd)):                       # Austausch: Spinskalierung
        n2 = np.maximum(2 * n1, 1e-30)
        kf = (3 * math.pi ** 2 * n2) ** (1 / 3)
        s2 = (2 * g1 / (2 * kf * n2)) ** 2
        fx = 1 + 0.804 - 0.804 / (1 + 0.2195149727645171 * s2 / 0.804)
        e = e - 0.5 * 0.75 * (3 / math.pi) ** (1 / 3) * n2 ** (4 / 3) * fx
    n = np.maximum(nu + nd, 1e-30)
    z = np.clip((nu - nd) / n, -1 + 1e-12, 1 - 1e-12)
    rs = (3 / (4 * math.pi * n)) ** (1 / 3)
    ec = _pw92(rs, z)
    phi = 0.5 * ((1 + z) ** (2 / 3) + (1 - z) ** (2 / 3))
    kf = (3 * math.pi ** 2 * n) ** (1 / 3)
    ks = np.sqrt(4 * kf / math.pi)
    t2 = ((gu + gd) / (2 * phi * ks * n)) ** 2
    beta, gam = 0.06672455060314922, (1 - math.log(2)) / math.pi ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        A = beta / gam / np.expm1(np.minimum(-ec / (gam * phi ** 3), 700))
        H = gam * phi ** 3 * np.log(1 + beta / gam * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))
    return e + n * (ec + np.nan_to_num(H))

def _xc_gga(nu, nd, r, h):
    """PBE: (ε_xc je Elektron, v↑, v↓); v = ∂f/∂ρ − (1/r²)·d/dr(r²·∂f/∂ρ'), Ableitungen nach ρ und ρ' als Differenzen."""
    nu, nd = np.maximum(nu, 0), np.maximum(nd, 0)
    gu, gd = np.gradient(nu, h) / r, np.gradient(nd, h) / r
    f = _pbe(nu, nd, gu, gd)
    arg = [nu, nd, gu, gd]
    d = []
    for i in range(4):
        dl = 1e-4 * np.abs(arg[i]) + 1e-12
        p, m = list(arg), list(arg)
        p[i], m[i] = arg[i] + dl, arg[i] - dl
        if i < 2:
            m[i] = np.maximum(m[i], 0)
            d.append((_pbe(*p) - _pbe(*m)) / (p[i] - m[i]))
        else:
            d.append((_pbe(*p) - _pbe(*m)) / (2 * dl))
    leer = (nu + nd) < 1e-10
    v = []
    for s in (0, 1):
        q = r * r * np.where(leer, 0, d[2 + s])
        v.append(np.where(leer, 0, d[s] - np.gradient(q, h) / r / (r * r)))
    n = np.maximum(nu + nd, 1e-30)
    return np.where(leer, 0, f / n), v[0], v[1]

def _hartree(q, r):
    """Hartree-Potential aus der Ladung je Gitterpunkt."""
    return (np.cumsum(q) - 0.5 * q) / r + np.cumsum((q / r)[::-1])[::-1] - 0.5 * q / r

def _besetzung_j(bel, kopplung):
    """Besetzung je Bahn (n, l, Spin s, j = l ∓ ½ als −1/+1). "hund": je Spin zuerst j = l − ½ (l Plätze je Spin), dann
    j = l + ½; "jj": zuerst j = l − ½ mit beiden Spins (2l Plätze), dann j = l + ½, je j die Spins nach Hund."""
    occ = {}
    for (n, l), m in bel.items():
        if l == 0:
            occ[(n, 0, 0, 1)], occ[(n, 0, 1, 1)] = min(m, 1), m - min(m, 1)
        elif kopplung == "hund":
            fu = min(m, 2 * l + 1)
            for s, fs in ((0, fu), (1, m - fu)):
                occ[(n, l, s, -1)], occ[(n, l, s, 1)] = min(fs, l), fs - min(fs, l)
        else:
            nm = min(m, 2 * l)
            occ[(n, l, 0, -1)], occ[(n, l, 1, -1)] = min(nm, l), nm - min(nm, l)
            occ[(n, l, 0, 1)], occ[(n, l, 1, 1)] = min(m - nm, l + 1), m - nm - min(m - nm, l + 1)
    return occ

def _scf(Z, bel, r, h, R, Ve=None, Kfest=None, maxit=150, schaetzung=None, selbst=False, term=False, skaliert=False, gradient=False,
         kopplung=None, spinbes=None, korrelation=True):
    """Selbstkonsistente Hülle; Gesamtenergie ohne und mit Spin-Bahn-Glied erster Ordnung (jj, j = l − ½ zuerst).
    selbst=True: das eigene Feld wirkt nicht auf sich selbst (je Elektron Hartree und Austausch-Korrelation der eigenen
    Dichte abgezogen, Perdew-Zunger).
    kopplung="hund" oder "jj": Spin und Bahn selbstkonsistent gekoppelt, je Bahn eigene Radialfunktion für j = l ∓ ½ (κ = l bzw.
    −l − 1), besetzt nach _besetzung_j; dann ist das Spin-Bahn-Glied enthalten (zweiter Rückgabewert = erster)."""
    occ = _besetzung_j(bel, kopplung) if kopplung else None
    # spinbes: Spinbesetzung (↑, ↓) einzelner Unterschalen vorgegeben statt nach Hund
    Vn = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r)
    B = r * r * h
    bes = {k: (min(m, 2 * k[1] + 1), m - min(m, 2 * k[1] + 1)) for k, m in bel.items()}
    bes.update(spinbes or {})
    ls = sorted({l for _, l in bel})
    Ve = [np.zeros_like(r), np.zeros_like(r)] if Ve is None else [v.copy() for v in Ve]
    sigma = -1.6 * Z * Z - 5
    hist_in, hist_res, E_alt = [], [], 0.0
    schaetzung = dict(schaetzung or {})
    vsic, esic_alt = {}, 0.0
    for it in range(maxit):
        Vs = [Vn + Ve[0], Vn + Ve[1]]
        Kf = Kfest if Kfest is not None else 1 / (1 - 0.5 * (Vs[0] + Vs[1]) / (2 * c * c))
        q, summe, eigen_q, vsic_erw = [np.zeros_like(r), np.zeros_like(r)], 0.0, {}, 0.0
        for s in (0, 1):
            for l, kappa, j in [(l, ka, j) for l in ls for ka, j in
                                (((-1, None),) if not occ or l == 0 else ((l, -1), (-l - 1, 1)))]:
                ns = sorted(n for (n, ll) in bel if ll == l)
                sk = (lambda n: (n, l, s)) if j is None else (lambda n: (n, l, s, j))
                alt = [schaetzung.get(sk(n)) for n in range(l + 1, max(ns) + 1)]
                w, v = _eigen(Vs[s], Kf, r, h, l, kappa, max(ns) - l, sigma, None if None in alt else alt)
                for n in range(l + 1, max(ns) + 1):
                    schaetzung[sk(n)] = w[n - l - 1]
                for n in ns:
                    f = bes[(n, l)][s] if not occ else occ[(n, l, s, 1 if j is None else j)]
                    if not f:
                        continue
                    wi, vi = w[n - l - 1], v[:, n - l - 1]
                    if selbst and sk(n) in vsic:                     # eigenes Feld abgezogen: Bahn im eigenen Potential
                        ws, vs_ = _eigen(Vs[s] + vsic[sk(n)], Kf, r, h, l, kappa, n - l, sigma,
                                         [schaetzung[sk(m)] for m in range(l + 1, n + 1)])
                        wi, vi = ws[n - l - 1], vs_[:, n - l - 1]
                        vsic_erw += f * float(np.sum(vsic[sk(n)] * vi ** 2 * B))
                    q[s] += f * vi ** 2 * B
                    summe += f * wi
                    eigen_q[sk(n)] = (f, vi ** 2 * B)
        dichte = [qq / (4 * math.pi * r ** 3 * h) for qq in q]
        qs = q[0] + q[1]
        VH = _hartree(qs, r)
        exc, vu, vd = _xc_gga(*dichte, r, h) if gradient else _xc(*dichte)
        if not korrelation:                                  # nur der örtlich gezählte Austausch (Slater), ohne angepasste Zahl
            ec, vcu, vcd = _xc_korrelation(*dichte)
            exc, vu, vd = exc - ec, vu - vcu, vd - vcd
        esic = 0.0
        if selbst:                                               # Selbstwirkung je Elektron: U[ρ_i] + E_xc[ρ_i, 0]
            for k, (f, q1) in eigen_q.items():
                VH1 = _hartree(q1, r)
                e1, v1, _ = _xc(q1 / (4 * math.pi * r ** 3 * h), np.zeros_like(r))
                esic -= f * (0.5 * np.sum(VH1 * q1) + np.sum(e1 * q1))
                neu = -(VH1 + v1)
                vsic[k] = neu if k not in vsic else 0.5 * (vsic[k] + neu)
        E = summe - np.sum(Ve[0] * q[0] + Ve[1] * q[1]) - vsic_erw + 0.5 * np.sum(VH * qs) + np.sum(exc * qs) + esic
        res = np.concatenate([VH + vu - Ve[0], VH + vd - Ve[1]])
        x_in = np.concatenate(Ve)
        if np.max(np.abs(res) * np.tile(r, 2)) < 1e-3 and abs(E - E_alt) < 1e-8 and abs(esic - esic_alt) < 1e-8:
            break
        E_alt, esic_alt = E, esic
        hist_in, hist_res = (hist_in + [x_in])[-6:], (hist_res + [res])[-6:]
        if len(hist_in) > 1:                              # Pulay-Mischung
            Rm = np.array(hist_res) * np.tile(r, 2)
            A = Rm @ Rm.T
            A = np.vstack([np.hstack([A, np.ones((len(A), 1))]), np.append(np.ones(len(A)), 0)])
            try:
                co = np.linalg.solve(A + 1e-12 * np.eye(len(A)), np.append(np.zeros(len(A) - 1), 1))[:-1]
            except np.linalg.LinAlgError:
                co = np.eye(len(hist_in))[-1]
            x_neu = sum(cc * (xi + 0.3 * ri) for cc, xi, ri in zip(co, hist_in, hist_res))
        else:
            x_neu = x_in + 0.3 * res
        Ve = [x_neu[:len(r)], x_neu[len(r):]]
    # skalierte ZORA: jede besetzte Bahn mit ε/(1 + S) statt ε (für ein Elektron exakt Dirac); skaliert als Menge von (n, l):
    # nur diese Unterschalen
    if skaliert:
        Vs = [Vn + Ve[0], Vn + Ve[1]]
        for s in (0, 1):
            for l in ls:
                ns = sorted(n for (n, ll) in bel if ll == l)
                w, v = _eigen(Vs[s], Kf, r, h, l, -1, max(ns) - l, sigma, [schaetzung[(n, l, s)] for n in range(l + 1, max(ns) + 1)])
                for n in ns:
                    f = bes[(n, l)][s]
                    if f and (skaliert is True or (n, l) in skaliert):
                        S_ = _skalierung(Kf, r, h, l, -1, v[:, n - l - 1])
                        E += f * (w[n - l - 1] / (1 + S_) - w[n - l - 1])
    if occ:                                                  # Spin-Bahn schon in den Bahnen
        return E, E, Ve, Kf, it, schaetzung
    # Spin-Bahn erster Ordnung in der spingemittelten Lage, offene Unterschalen l ≥ 1
    so = 0.0
    Vm = Vn + 0.5 * (Ve[0] + Ve[1])
    for (n, l), m in bel.items():
        if l == 0 or m == 2 * (2 * l + 1):
            continue
        sch = [0.5 * (schaetzung[(k, l, 0)] + schaetzung[(k, l, 1)]) for k in range(l + 1, n + 1)]
        wm, vm_ = _eigen(Vm, Kf, r, h, l, l, n - l, sigma, sch)                    # j = l − ½ (κ = l)
        wp, vp_ = _eigen(Vm, Kf, r, h, l, -l - 1, n - l, sigma, sch)               # j = l + ½
        if skaliert:
            wm = wm[n - l - 1] / (1 + _skalierung(Kf, r, h, l, l, vm_[:, n - l - 1]))
            wp = wp[n - l - 1] / (1 + _skalierung(Kf, r, h, l, -l - 1, vp_[:, n - l - 1]))
        else:
            wm, wp = wm[n - l - 1], wp[n - l - 1]
        mittel = (2 * l * wm + (2 * l + 2) * wp) / (4 * l + 2)
        unten = min(m, 2 * l)
        so += unten * (wm - mittel) + (m - unten) * (wp - mittel)
    if term:                                                 # Ausrichtung der Windungen (größtes L) der offenen Schalen
        dt = _term(Z, bel, bes, r, h, R, Ve, Kf, schaetzung)
        E += dt
    return E, E + so, Ve, Kf, it, schaetzung

_CK = {}

def _ck(k, l, m, l2, m2):
    """Gaunt-Beiwert c^k(l m, l2 m2) der Kugelflächenfunktionen."""
    key = (k, l, m, l2, m2)
    if key not in _CK:
        from sympy.physics.wigner import wigner_3j
        _CK[key] = float((-1) ** m * math.sqrt((2 * l + 1) * (2 * l2 + 1)) * wigner_3j(l, k, l2, 0, 0, 0)
                         * wigner_3j(l, k, l2, -m, m - m2, m2))
    return _CK[key]

def _term(Z, bel, bes, r, h, R, Ve, Kf, sch, zeta=None):
    """Ausrichtung der Windungen in den offenen Unterschalen (Hundsche Regel des größten L): tiefster Term gegen das Mittel
    aller Determinanten mit derselben Spinbesetzung, aus den Slater-Integralen der eigenen Bahnen [Hartree].
    Mit zeta ({(n, l): ζ}) beide Umlaufsinne gekoppelt: Abstoßung und ζ·l·s in einer Matrix über alle Spinbesetzungen (Block
    kleinstes |M_J|); zurück das tiefste Niveau gegen den tiefsten Term der Hund-Besetzung ohne Spin-Bahn [Hartree]."""
    from itertools import combinations, product
    offen = [k for k, m in bel.items() if 0 < m < 2 * (2 * k[1] + 1)]
    if not offen:
        return 0.0
    Vm = np.where(r < R, -Z / (2 * R) * (3 - (r / R) ** 2), -Z / r) + 0.5 * (Ve[0] + Ve[1])
    phi = {}
    for (n, l) in offen:
        w, v = _eigen(Vm, Kf, r, h, l, -1, n - l, -1.6 * Z * Z - 5,
                      [0.5 * (sch[(m, l, 0)] + sch[(m, l, 1)]) for m in range(l + 1, n + 1)])
        phi[(n, l)] = v[:, n - l - 1]
    w2 = r * r * h
    RK = {}

    def rk(k, a, b, c, d):                                    # R^k(ab, cd) = ∫∫ P_a P_c r<^k/r>^(k+1) P_b P_d
        key = (k, a, b, c, d)
        if key not in RK:
            q1, q2 = phi[a] * phi[c] * w2, phi[b] * phi[d] * w2
            y = (np.cumsum(q1 * r ** k) - 0.5 * q1 * r ** k) / r ** (k + 1) + \
                r ** k * (np.cumsum((q1 / r ** (k + 1))[::-1])[::-1] - 0.5 * q1 / r ** (k + 1))
            RK[key] = float(np.sum(y * q2))
        return RK[key]
    so = [(sh, m, s) for sh in offen for m in range(-sh[1], sh[1] + 1) for s in (0, 1)]

    def zwei(p, q, t, u):                                     # ⟨pq|tu⟩ über Spinbahnen
        (a, ma, sa), (b, mb, sb), (c, mc, sc), (d, md, sd) = so[p], so[q], so[t], so[u]
        if sa != sc or sb != sd or ma + mb != mc + md:
            return 0.0
        la, lb, lc, ld = a[1], b[1], c[1], d[1]
        return sum(_ck(k, la, ma, lc, mc) * _ck(k, ld, md, lb, mb) * rk(k, a, b, c, d)
                   for k in range(abs(la - lc), la + lc + 1, 2) if abs(ld - lb) <= k <= ld + lb)
    teile, hund = [], []
    for sh in offen:
        idx = {s_: [i for i, x in enumerate(so) if x[0] == sh and x[2] == s_] for s_ in (0, 1)}
        hund.append([tuple(sorted(cu + cd)) for cu in combinations(idx[0], bes[sh][0])
                     for cd in combinations(idx[1], bes[sh][1])])
        teile.append(hund[-1] if zeta is None else list(combinations(idx[0] + idx[1], sum(bes[sh]))))
    dets = [tuple(sorted(sum(t, ()))) for t in product(*teile)]
    if zeta is not None:                                       # nur der Block kleinsten |M_J| (2·M_J ganzzahlig)
        mj2 = {D: sum(2 * so[i][1] + (1 - 2 * so[i][2]) for i in D) for D in dets}
        dets = [D for D in dets if mj2[D] == len(dets[0]) % 2]
    if len(dets) == 1 and zeta is None or len(dets) > 3000 and zeta is None:
        return 0.0

    def eins(p, t):                                           # ⟨p| ζ l·s |t⟩, s = 0: +½, s = 1: −½
        (a, ma, sa), (c, mc, sc) = so[p], so[t]
        if a != c:
            return 0.0
        z, l = zeta[a], a[1]
        if sa == sc and ma == mc:
            return z * ma * (0.5 - sa)
        if sa == 1 and sc == 0 and ma == mc + 1:              # l₊s₋ : (m, ↑) → (m + 1, ↓)
            return 0.5 * z * math.sqrt(l * (l + 1) - mc * (mc + 1))
        if sa == 0 and sc == 1 and ma == mc - 1:
            return 0.5 * z * math.sqrt(l * (l + 1) - mc * (mc - 1))
        return 0.0
    ml = {D: sum(so[i][1] for i in D) for D in dets}

    def element(D1, D2):
        if D1 == D2:
            return sum(zwei(p, q, p, q) - zwei(p, q, q, p) for i, p in enumerate(D1) for q in D1[i + 1:]) + \
                (sum(eins(p, p) for p in D1) if zeta is not None else 0.0)
        nur1 = [x for x in D1 if x not in D2]
        nur2 = [x for x in D2 if x not in D1]
        if len(nur1) > 2:
            return 0.0
        vz = (-1) ** (sum(D1.index(x) for x in nur1) + sum(D2.index(x) for x in nur2))
        if len(nur1) == 1:
            p, t = nur1[0], nur2[0]
            return vz * (sum(zwei(p, q, t, q) - zwei(p, q, q, t) for q in D1 if q != p) +
                         (eins(p, t) if zeta is not None else 0.0))
        (p, q), (t, u) = nur1, nur2
        return vz * (zwei(p, q, t, u) - zwei(p, q, u, t))
    if zeta is not None:                                       # gegen den tiefsten Term der Hund-Besetzung (größtes S)
        mit = _tiefstes_niveau(dets, len(so), zwei, eins)
        zeta = None
        dh = [tuple(sorted(sum(t, ()))) for t in product(*hund)]
        mlh = {D: sum(so[i][1] for i in D) for D in dh}
        tief_h = min(np.linalg.eigvalsh(np.array([[element(D1, D2) for D2 in blk] for D1 in blk]))[0]
                     for L in set(mlh.values()) for blk in [[D for D in dh if mlh[D] == L]])
        return float(mit - tief_h)
    tief, spur = 0.0, 0.0
    for L in set(ml.values()):
        blk = [D for D in dets if ml[D] == L]
        H = np.array([[element(D1, D2) for D2 in blk] for D1 in blk])
        spur += np.trace(H)
        tief = min(tief, np.linalg.eigvalsh(H)[0]) if tief else np.linalg.eigvalsh(H)[0]
    return float(tief - spur / len(dets))

def _tiefstes_niveau(dets, P, zwei, eins):
    """Tiefster Eigenwert von Abstoßung und Spin-Bahn im Raum der Determinanten dets (Tupel der besetzten Spinbahnen 0 … P−1):
    Matrix dünn nach den Slater-Condon-Regeln (nur Paare, die sich in höchstens zwei Bahnen unterscheiden), tiefster
    Eigenwert nach Lanczos; dieselben Elemente wie die volle Matrix."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import eigsh
    A = np.zeros((P, P, P, P))
    for p in range(P):
        for q in range(P):
            for s in range(P):
                for u in range(P):
                    A[p, q, s, u] = zwei(p, q, s, u) - zwei(p, q, u, s)
    E1 = np.array([[eins(p, s) for s in range(P)] for p in range(P)])
    B = np.einsum("pqrq->prq", A)                              # B[p, r, q] = A[p, q, r, q]
    m = np.array([sum(1 << i for i in D) for D in dets], dtype=np.int64)
    O = ((m[:, None] >> np.arange(P)) & 1).astype(float)
    diag = 0.5 * np.einsum("ip,pq,iq->i", O, np.einsum("pqpq->pq", A), O) + O @ np.diag(E1)
    pos = lambda D, b: np.bitwise_count(D & ((np.int64(1) << b) - 1))
    bit = lambda v: np.log2(v.astype(float)).round().astype(np.int64)
    zeilen, spalten, werte = [np.arange(len(m))], [np.arange(len(m))], [diag]
    for i0 in range(0, len(m), 256):
        Di = m[i0:i0 + 256, None]
        x = Di ^ m[None, :]
        c = np.bitwise_count(x)
        ii, jj = np.nonzero(((c == 2) | (c == 4)) & (np.arange(len(m))[None, :] > np.arange(i0, i0 + len(Di))[:, None]))
        a_, b_ = m[ii + i0], m[jj]
        weg, hin = a_ & ~b_, b_ & ~a_
        einf = np.bitwise_count(weg) == 1
        v = np.zeros(len(ii))
        if einf.any():                                         # ein Elektron versetzt
            p, s = bit(weg[einf]), bit(hin[einf])
            vz = (-1.0) ** (pos(a_[einf], p) + pos(b_[einf], s))
            summe = np.einsum("kq,kq->k", O[ii[einf] + i0], B[p, s, :]) - B[p, s, p]
            v[einf] = vz * (E1[p, s] + summe)
        zw = ~einf                                             # zwei Elektronen versetzt
        if zw.any():
            lo_w = weg[zw] & -weg[zw]; lo_h = hin[zw] & -hin[zw]
            p, q, s, u = bit(lo_w), bit(weg[zw] ^ lo_w), bit(lo_h), bit(hin[zw] ^ lo_h)
            vz = (-1.0) ** (pos(a_[zw], p) + pos(a_[zw], q) + pos(b_[zw], s) + pos(b_[zw], u))
            v[zw] = vz * A[p, q, s, u]
        k = v != 0
        zeilen += [ii[k] + i0, jj[k]]
        spalten += [jj[k], ii[k] + i0]
        werte += [v[k], v[k]]
    H = coo_matrix((np.concatenate(werte), (np.concatenate(zeilen), np.concatenate(spalten))), shape=(len(m), len(m))).tocsr()
    if len(m) < 60:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    return float(eigsh(H, k=1, which="SA", tol=1e-10)[0][0])


def _skalierung(Kf, r, h, l, kappa, phi):
    """Skalierte ZORA: S = ⟨φ| σ·p (c²/(2c² − V)²) σ·p |φ⟩ = ⟨φ| ½ σ·p (K²/2c²) σ·p |φ⟩; die Bahnenergie ε/(1 + S) ist für ein
    Elektron exakt die Dirac-Energie."""
    d, e, B = _matrix(np.zeros_like(r), Kf * Kf / (2 * c * c), r, h, l, kappa)
    return float(phi @ (d * phi) + 2 * np.sum(e * phi[:-1] * phi[1:]))

def kerne_belegung_nl(Z):
    return {(n, l): min(Z, b) - a + 1 for n, l, a, b in K.fuellung(9) if a <= Z}

def _aussen(Z):
    """Äußere Unterschalen der Periode p von Z (höchstes s nach n + l): ps, pp, (p−1)d, (p−2)f, (p−3)g."""
    p = max(n for (n, l) in kerne_belegung_nl(Z) if l == 0)
    return [(n, l) for n, l in ((p, 0), (p, 1), (p - 1, 2), (p - 2, 3), (p - 3, 4)) if n > l]
