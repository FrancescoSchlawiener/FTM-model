"""Etappe 4: die Leiter der Perioden als Spiegelform. Feld in log r um R gespiegelt (r ↔ R²/r), Schwellen N(N + 1),
Zustände an der Schwelle als Kugelfunktionen."""
import math
import numpy as np


# ---------------------------------------------------------------------------------------------------------------
# Die Leiter der Perioden als schließende Form: ein Feld, das in log r um R gespiegelt ist (r ↔ R²/r), an der Schwelle
# E = 0. −u'' + l(l+1)/r²·u = v·W(r)·u mit W = 1/(r(1 + r)²) (R = 1): die Kopplungen v hängen nur von N = n + l ab,
# v = N(N + 1) (in log r ein Pöschl-Teller-Topf). Zu jedem N gehören 2⌈N/2⌉² Zustände.
# ---------------------------------------------------------------------------------------------------------------
def schliessform(l_max=3, anzahl=6, punkte=3000):
    """Kopplungen v(n, l) an der Schwelle im gespiegelten Feld; gegen N(N + 1) mit N = n + l."""
    from scipy.linalg import eigh
    x = np.linspace(math.log(1e-6), math.log(1e6), punkte)
    r, h = np.exp(x), x[1] - x[0]
    aus = {}
    for l in range(l_max + 1):
        i = np.arange(len(r))
        A = np.zeros((len(r), len(r)))
        A[i, i] = 2 / h ** 2 + (l + 0.5) ** 2
        A[i[:-1], i[:-1] + 1] = A[i[:-1] + 1, i[:-1]] = -1 / h ** 2
        v = eigh(A, np.diag(r / (1 + r) ** 2), subset_by_index=(0, anzahl - 1), eigvals_only=True)
        for k, vk in enumerate(v):
            aus[(l + 1 + k, l)] = float(vk)
    return aus


def kugelform(l_max=3, anzahl=5, punkte=3000):
    """Die Zustände des Spiegelfelds an der Schwelle sind Kugelfunktionen: y = u/√r = P_N^{2l+1}(t), t = tanh(x/2) = (r − R)/(r + R),
    N = n + l. Zwei Wege: exakt (sympy, Rest der Gleichung −y'' + (l + ½)²y − N(N + 1)/(4cosh²(x/2))·y) und numerisch (Überlapp
    mit dem Eigenvektor der Schwellenrechnung). Gespiegelt r ↔ R²/r ist t ↔ −t, der Äquator der Kugel; Parität (−1)^(n − l − 1)."""
    import sympy as sp
    from scipy.linalg import eigh
    from scipy.special import lpmv
    X = sp.symbols("x", real=True)
    rest = 0
    for l in range(l_max + 1):
        for k in range(anzahl):
            N = 2 * l + 1 + k
            y = sp.assoc_legendre(N, 2 * l + 1, sp.tanh(X / 2))
            g = -sp.diff(y, X, 2) + sp.Rational(2 * l + 1, 2) ** 2 * y - N * (N + 1) / (4 * sp.cosh(X / 2) ** 2) * y
            rest = max(rest, max(abs(float(g.subs(X, xv))) / max(1e-300, abs(float(y.subs(X, xv))) + 1e-12)
                                 for xv in (-3.1, -0.7, 0.4, 2.3)))
    x = np.linspace(math.log(1e-6), math.log(1e6), punkte)
    r, h = np.exp(x), x[1] - x[0]
    t = np.tanh(x / 2)
    ueber = 1.0
    for l in range(l_max + 1):
        i = np.arange(len(r))
        A = np.zeros((len(r), len(r)))
        A[i, i] = 2 / h ** 2 + (l + 0.5) ** 2
        A[i[:-1], i[:-1] + 1] = A[i[:-1] + 1, i[:-1]] = -1 / h ** 2
        v, vec = eigh(A, np.diag(r / (1 + r) ** 2), subset_by_index=(0, anzahl - 1))
        for k in range(anzahl):
            N = 2 * l + 1 + k
            p = lpmv(2 * l + 1, N, t)
            a = vec[:, k]
            ueber = min(ueber, abs(np.dot(a, p)) / (np.linalg.norm(a) * np.linalg.norm(p)))
            assert np.allclose(p[::-1], (-1) ** k * p, atol=1e-9 * np.max(abs(p)))   # Parität am Äquator
    return dict(rest=rest, ueberlapp=float(ueber))


def gehaltener_trichter():
    """Die Hülle als von der Zählregel der Elektronen gehaltener Trichter (wie der Weiße Zwerg): Lane-Emden n = 3/2 um den Kern,
    in der Form von Thomas-Fermi, x·φ'' = φ^{3/2}·√x, R ∝ Z^{−1/3}. Höhe des Felds v = 8·max(−r²V) = k·Z^{2/3}; Abweichung des Felds
    von der Spiegelform v·f(1 − f) innen und außen vom Spiegelradius (feste Breite, ohne freien Wert); Stufe N_v (größtes N mit
    N(N + 1) < v) der Edelgase gegen N = p + 1; Höhe bei Z = 168."""
    from scipy.integrate import solve_ivp
    s = solve_ivp(lambda t, y: [y[1], max(y[0], 0) ** 1.5 / math.sqrt(t)], [1e-10, 300], [1.0, -1.588071022611375],
                  rtol=1e-11, atol=1e-13, dense_output=True)
    x = np.linspace(math.log(1e-4), math.log(200), 20000)
    F = np.exp(x) * s.sol(np.exp(x))[0]
    i = int(np.argmax(F))
    d = (F - F[i] / np.cosh((x - x[i]) / 2) ** 2) / F[i]
    w = F > 0.05 * F[i]
    k = 8 * 0.5 * (3 * math.pi / 4) ** (2 / 3) * float(F[i])
    nv = lambda v: max(n for n in range(14) if n * (n + 1) < v)
    edel = {2: 1, 10: 2, 18: 3, 36: 4, 54: 5, 86: 6, 118: 7}
    return dict(k=k, innen=float(np.sqrt(np.mean(d[w & (x < x[i])] ** 2))), aussen=float(np.sqrt(np.mean(d[w & (x > x[i])] ** 2))),
                edel={Z: dict(p=p, v=k * Z ** (2 / 3), N=nv(k * Z ** (2 / 3))) for Z, p in edel.items()}, v168=k * 168 ** (2 / 3))
