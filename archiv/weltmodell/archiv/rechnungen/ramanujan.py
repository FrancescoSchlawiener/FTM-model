# -*- coding: utf-8 -*-
"""
Forschung 04 - Ramanujan-Werkzeuge fuer das FTM.
Numerische Pruefung aller Formeln aus forschung/04-ramanujan-werkzeuge.md.
Ausgabe: rechnungen/ramanujan.txt   (python ramanujan.py > ramanujan.txt)

Konventionen:
  tau in der oberen Halbebene, q = exp(2 pi i tau).
  FTM-Kopplung (Forschung 01): q = exp(-sqrt(y)), y = a_N/a0  <=>  tau = i*t, t = sqrt(y)/(2 pi).
  Kopplung aus Punkt 10:        tau = i*y^(1/4)/(2 pi).
Nur Standardbibliothek + numpy.
"""
import sys, math, cmath, itertools
from decimal import Decimal as Dm, getcontext
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PI = math.pi
PHI = (1 + math.sqrt(5)) / 2
EULER_GAMMA = 0.5772156649015329
LOG10 = math.log10


def hdr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# --------------------------------------------------------------------------
# Physikalische Eingaben (wie forschung_uebergang.py)
# --------------------------------------------------------------------------
G, C = 6.674e-11, 2.998e8
MSUN, AU, KPC = 1.989e30, 1.496e11, 3.0857e19
A0 = 1.2e-10
RSTAR_MW = math.sqrt(G * 6e10 * MSUN / A0)          # Milchstrasse, 6e10 Msun
Y_EARTH = G * MSUN / AU**2 / A0
Y_SATURN = G * MSUN / (9.58 * AU)**2 / A0
Y_NEPTUN = G * MSUN / (30.1 * AU)**2 / A0

# --------------------------------------------------------------------------
# Modulare Grundfunktionen
# --------------------------------------------------------------------------
def nome(tau):
    return cmath.exp(2j * PI * tau)


def lam(q, e, sq=False, tol=1e-24):
    """sum_{n>=1} n^e q^n/(1-q^n)   (sq=True: .../(1-q^n)^2)"""
    s = 0j
    qn = 1 + 0j
    n = 0
    while True:
        n += 1
        qn *= q
        den = (1 - qn) ** 2 if sq else (1 - qn)
        s += n**e * qn / den
        if n > 8 and abs(n**e * qn) < tol * (1 + abs(s)):
            return s
        if n > 500000:
            raise RuntimeError("Lambert-Reihe konvergiert nicht")


def E2q(q): return 1 - 24 * lam(q, 1)
def E4q(q): return 1 + 240 * lam(q, 3)
def E6q(q): return 1 - 504 * lam(q, 5)
def DE2q(q): return -24 * lam(q, 2, True)      # q d/dq
def DE4q(q): return 240 * lam(q, 4, True)
def DE6q(q): return -504 * lam(q, 6, True)
def E2(t): return E2q(nome(t))
def E4(t): return E4q(nome(t))
def E6(t): return E6q(nome(t))


def eta(tau):
    q = nome(tau)
    p = 1 + 0j
    qn = 1 + 0j
    while True:
        qn *= q
        p *= (1 - qn)
        if abs(qn) < 1e-22:
            break
    return cmath.exp(2j * PI * tau / 24) * p


def Delta(tau):
    return eta(tau) ** 24


def E2_it(t):
    """E2(i t) reell, fuer t<1 ueber E2(i/t) = -t^2 E2(it) + 6t/pi"""
    if t >= 1:
        return E2q(math.exp(-2 * PI * t)).real
    return (6 * t / PI - E2q(math.exp(-2 * PI / t)).real) / t**2


def Delta_it(t):
    if t >= 1:
        return Delta(1j * t).real
    return Delta(1j / t).real / t**12


def sigma_list(N, k=1):
    s = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            s[m] += d**k
    return s


def divcount_list(N):
    return sigma_list(N, 0)


SIG1 = sigma_list(20000, 1)


def nu_rar(y):
    return 1.0 / (1.0 - math.exp(-math.sqrt(y)))


YGRID = [1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1e3]
YDENSE = [10 ** (-2 + 4 * i / 400) for i in range(401)]       # 0.01 ... 100

# ==========================================================================
hdr("1  Serre-Ramanujan-Ableitung D_k = q d/dq - (k/12) E2  (E2 als Zusammenhang)")
# ==========================================================================
tau0 = 0.3 + 1.1j
print(f"Probepunkt tau = {tau0}")
print("Transformationsgesetz E2(g tau) = (c tau+d)^2 E2(tau) + 12 c (c tau+d)/(2 pi i):")
for (a, b, c, d) in [(0, -1, 1, 0), (2, 1, 5, 3), (1, 0, 3, 1), (5, 2, 7, 3)]:
    gt = (a * tau0 + b) / (c * tau0 + d)
    lhs = E2(gt)
    rhs = (c * tau0 + d) ** 2 * E2(tau0) + 12 * c * (c * tau0 + d) / (2j * PI)
    print(f"   g=({a},{b};{c},{d})  Im(g tau)={gt.imag:.4f}   rel. Fehler = {rel(lhs, rhs):.2e}")
print("   Inhomogener Term = (12/2 pi i) (c tau+d)^2 * d/dtau log(c tau+d)  ->  Eichtransformation A -> g A g^-1 + g d g^-1")

a, b, c, d = 2, 1, 5, 3
gt = (a * tau0 + b) / (c * tau0 + d)
j = c * tau0 + d
q0, qg = nome(tau0), nome(gt)
DE4_0, DE4_g = DE4q(q0), DE4q(qg)
E4_0, E2_0, E6_0 = E4q(q0), E2q(q0), E6q(q0)
E4_g, E2_g = E4q(qg), E2q(qg)
S4_0 = DE4_0 - E2_0 * E4_0 / 3
S4_g = DE4_g - E2_g * E4_g / 3
print("\nKovarianz fuer f = E4 (Gewicht 4) unter g=(2,1;5,3):")
print(f"   gewoehnliche Ableitung: DE4(g tau)/((c tau+d)^6 DE4(tau)) = {DE4_g / (j**6 * DE4_0):.6f}  (nicht 1)")
stoer = 4 * c * j**5 * E4_0 / (2j * PI)
print(f"   Stoerterm DE4(g tau)-(c tau+d)^6 DE4(tau) gegen 4c(c tau+d)^5 E4/(2 pi i): rel. Fehler {rel(DE4_g - j**6 * DE4_0, stoer):.2e}")
print(f"   Serre D4E4(g tau)/((c tau+d)^6 D4E4(tau)) = {S4_g / (j**6 * S4_0):.12f}")
print(f"   D4 E4 = -E6/3 :  rel. Fehler {rel(S4_0, -E6_0 / 3):.2e}")
h = 1e-5
dDelta = (Delta(tau0 + h) - Delta(tau0 - h)) / (2 * h) / (2j * PI)
print(f"   D12 Delta = q dDelta/dq - E2 Delta  (Differenzenquotient): |D12 Delta|/|Delta| = {abs(dDelta - E2_0 * Delta(tau0)) / abs(Delta(tau0)):.2e}")

# Nullstelle von E2 auf der imaginaeren Achse -> Maximum des Delta-Kerns entlang r
lo, hi = 0.3, 0.9
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if E2_it(mid) > 0:
        hi = mid
    else:
        lo = mid
t0 = 0.5 * (lo + hi)
y0 = (2 * PI * t0) ** 2
print(f"\nNullstelle E2(i t0) = 0:  t0 = {t0:.10f}   (E2(i) = {E2_it(1.0):.10f} = 3/pi = {3 / PI:.10f})")
print("Mit tau = i r*/(2 pi r):  d ln Delta / d ln r = sqrt(y) * E2(tau)   (aus D12 Delta = 0)")
for y in [1.0, 10.0, 50.0]:
    t = math.sqrt(y) / (2 * PI)
    # r/r* = 1/sqrt(y); ln r = -0.5 ln y ; numerische Ableitung
    hh = 1e-5
    f = lambda yy: math.log(Delta_it(math.sqrt(yy) / (2 * PI)))
    dnum = (f(y * math.exp(hh)) - f(y * math.exp(-hh))) / (2 * hh)   # d ln Delta / d ln y
    slope_num = -2 * dnum
    print(f"   y = {y:5.1f}:  numerisch {slope_num:+.8f}   sqrt(y) E2 = {math.sqrt(y) * E2_it(t):+.8f}")
print(f"Delta(tau(r)) ist maximal, wo E2 = 0:  sqrt(y0) = 2 pi t0 = {2 * PI * t0:.6f},  y0 = a_N/a0 = {y0:.4f}")
print(f"   r_max = r*/(2 pi t0) = {1 / (2 * PI * t0):.4f} r*   (Milchstrasse r* = {RSTAR_MW / KPC:.2f} kpc -> {RSTAR_MW / KPC / (2 * PI * t0):.2f} kpc),  a_N = {y0 * A0:.3e} m/s^2")
ts = [0.5 + 0.001 * i for i in range(1001)]
tmax = max(ts, key=lambda t: t**6 * Delta_it(t))
print(f"   invariantes |Delta| (Im tau)^6: Maximum bei t = {tmax:.3f}  (selbstdual tau = i, y = 4 pi^2 = {4 * PI**2:.3f}, r = r*/(2 pi))")

# ==========================================================================
hdr("2  Ramanujans Differentialsystem (P,Q,R) und Chazy-Gleichung als parametrischer Fluss")
# ==========================================================================
for q in [0.05, 0.3, 0.7, 0.05 + 0.2j]:
    P, Q, R = E2q(q), E4q(q), E6q(q)
    e1 = rel(DE2q(q), (P * P - Q) / 12)
    e2 = rel(DE4q(q), (P * Q - R) / 3)
    e3 = rel(DE6q(q), (P * R - Q * Q) / 2)
    print(f"   q = {q!s:12s}  rel. Fehler  qP' {e1:.1e}   qQ' {e2:.1e}   qR' {e3:.1e}")


def DjP(q, jj, N=4000):
    s = 0j
    qm = 1 + 0j
    for m in range(1, N):
        qm *= q
        s += m**jj * SIG1[m] * qm
        if abs(qm) * m**(jj + 2) < 1e-22:
            break
    return (1 if jj == 0 else 0) - 24 * s if jj == 0 else -24 * s


print("Chazy  y''' = 2 y y'' - 3 y'^2  fuer y = pi i E2(tau), Ableitungen gliedweise aus der q-Reihe:")
for tau in [0.1 + 0.8j, 0.25 + 0.35j, -0.4 + 0.2j]:
    q = nome(tau)
    P = E2q(q)
    yv = PI * 1j * P
    y1 = PI * 1j * (2j * PI) * DjP(q, 1)
    y2 = PI * 1j * (2j * PI) ** 2 * DjP(q, 2)
    y3 = PI * 1j * (2j * PI) ** 3 * DjP(q, 3)
    res = y3 - 2 * yv * y2 + 3 * y1**2
    print(f"   tau = {tau}:  |Residuum|/|y'''| = {abs(res) / abs(y3):.2e}")

J = np.array([[1 / 6, -1 / 12, 0], [1 / 3, 1 / 3, -1 / 3], [1 / 2, -1, 1 / 2]])
ev, V = np.linalg.eig(J)
print("\nLinearisierung am Newton-Fixpunkt (P,Q,R) = (1,1,1):  J =", J.tolist())
print("   Eigenwerte:", np.round(ev.real, 12).tolist(), "  Rang(J) =", np.linalg.matrix_rank(J), " Rang(J^2) =", np.linalg.matrix_rank(J @ J))
i1 = int(np.argmin(abs(ev - 1)))
v1 = V[:, i1].real
print("   Eigenvektor zu 1 (normiert auf -24):", np.round(v1 / v1[0] * -24, 9).tolist(), " = Eisenstein-Koeffizienten (-24, 240, -504)")
print("   Eigenwert 0 doppelt mit Rang(J)=2 -> Jordan-Block: Abweichungen ~ ln q = -sqrt(y) (Potenzgesetz, an der Spitze nicht holomorph)")


def F_ram(X):
    P, Q, R = X
    return np.array([(P * P - Q) / 12, (P * Q - R) / 3, (P * R - Q * Q) / 2])


def rk4_flow(X, s0, s1, n):
    hstep = (s1 - s0) / n
    for _ in range(n):
        k1 = F_ram(X)
        k2 = F_ram(X + 0.5 * hstep * k1)
        k3 = F_ram(X + 0.5 * hstep * k2)
        k4 = F_ram(X + hstep * k3)
        X = X + hstep * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return X


print("\nFluss in s = ln q (= -sqrt(y)) vom Newton-Fixpunkt aus, Start X = (1,1,1) + c*q0*(-24,240,-504), q0 = 1e-9:")
for cc in [1.0, 0.5]:
    q0s = 1e-9
    X0 = np.array([1.0, 1.0, 1.0]) + cc * q0s * np.array([-24.0, 240.0, -504.0])
    for qend in [0.3, 0.6]:
        X = rk4_flow(X0, math.log(q0s), math.log(qend), 20000)
        ex = np.array([E2q(cc * qend).real, E4q(cc * qend).real, E6q(cc * qend).real])
        print(f"   c = {cc}:  q = {qend}:  Fluss P,Q,R = {X[0]:+.10f} {X[1]:+.8f} {X[2]:+.8f}   Reihe bei c*q: max rel. Fehler {max(abs(X - ex) / abs(ex)):.1e}")
print("   -> Die einzige an der Spitze regulaere Loesungsschar ist E(c q): c = Skalenwahl von q (Wahl von a0).")
print("\nProfil entlang der Galaxie (q = exp(-sqrt y)):")
print(f"   {'y':>8s} {'q':>9s} {'P=E2':>10s} {'Q=E4':>11s} {'R=E6':>11s} {'Delta':>11s} {'dlnDelta/dlnr':>13s}")
for y in [100, 4 * PI**2, y0, 10, 1, 0.1]:
    t = math.sqrt(y) / (2 * PI)
    q = math.exp(-math.sqrt(y))
    if t >= 1:
        P, Q, R = E2q(q).real, E4q(q).real, E6q(q).real
    else:
        qq = math.exp(-2 * PI / t)
        P = E2_it(t)
        Q = E4q(qq).real / t**4
        R = -E6q(qq).real / t**6
    print(f"   {y:8.3f} {q:9.3e} {P:+10.5f} {Q:11.4e} {R:+11.4e} {Delta_it(t):11.4e} {math.sqrt(y) * P:+13.5f}")

# ==========================================================================
hdr("3  Hecke-Operatoren, Multiplikativitaet, Ramanujan-Petersson")
# ==========================================================================
NT = 1000


def mul_trunc(a, b, N):
    r = [0] * (N + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        lim = N - i
        for jj in range(0, min(lim, len(b) - 1) + 1):
            bj = b[jj]
            if bj:
                r[i + jj] += ai * bj
    return r


Eu = [0] * (NT + 1)
k = 0
while True:
    g1 = k * (3 * k - 1) // 2
    if g1 > NT:
        break
    Eu[g1] += (-1) ** k
    if k > 0:
        g2 = k * (3 * k + 1) // 2
        if g2 <= NT:
            Eu[g2] += (-1) ** k
    k += 1
E2p = mul_trunc(Eu, Eu, NT)
E4p = mul_trunc(E2p, E2p, NT)
E8p = mul_trunc(E4p, E4p, NT)
E16p = mul_trunc(E8p, E8p, NT)
E24p = mul_trunc(E16p, E8p, NT)
TAU = [0] + E24p[:NT]            # TAU[n] = tau(n), n = 1..NT
print("tau(1..12) =", TAU[1:13])
primes = [p for p in range(2, NT + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
viol = 0
cnt = 0
for m in range(2, NT + 1):
    for n in range(2, NT // m + 1):
        if math.gcd(m, n) == 1:
            cnt += 1
            if TAU[m * n] != TAU[m] * TAU[n]:
                viol += 1
print(f"Multiplikativitaet tau(mn)=tau(m)tau(n), ggT=1, mn<=1000: {cnt} Paare, Verletzungen: {viol}")
viol = 0
cnt = 0
for p in primes:
    pk = [1, p]
    while pk[-1] * p <= NT:
        pk.append(pk[-1] * p)
    for kk in range(1, len(pk) - 1):
        cnt += 1
        if TAU[pk[kk + 1]] != TAU[p] * TAU[pk[kk]] - p**11 * TAU[pk[kk - 1]]:
            viol += 1
print(f"Rekursion tau(p^(k+1)) = tau(p)tau(p^k) - p^11 tau(p^(k-1)): {cnt} Faelle, Verletzungen: {viol}")
xs = [TAU[p] / p**5.5 for p in primes]
print(f"Ramanujan-Petersson |tau(p)| <= 2 p^(11/2), p <= {NT} ({len(primes)} Primzahlen): max |tau(p)|/(2 p^5.5) = {max(abs(x) for x in xs) / 2:.4f}")
dl = divcount_list(NT)
print(f"Deligne |tau(n)| <= d(n) n^(11/2), n <= {NT}: max = {max(abs(TAU[n]) / (dl[n] * n**5.5) for n in range(1, NT + 1)):.4f}")
m2 = sum(x * x for x in xs) / len(xs)
m4 = sum(x**4 for x in xs) / len(xs)
print(f"Sato-Tate-Momente von x_p = tau(p)/p^5.5 = 2cos(theta_p):  <x^2> = {m2:.3f} (Erwartung 1),  <x^4> = {m4:.3f} (Erwartung 2)")
x2 = TAU[2] / 2**5.5


def chebU(kk, x):
    u0, u1 = 1.0, 2 * x
    if kk == 0:
        return u0
    for _ in range(kk - 1):
        u0, u1 = u1, 2 * x * u1 - u0
    return u1


print("Skalenleiter p = 2:  tau(2^k)/2^(5.5k) = U_k(cos theta_2), |.| <= k+1:")
for kk in range(1, 10):
    print(f"   k = {kk}:  {TAU[2**kk] / 2**(5.5 * kk):+.6f}   U_k = {chebU(kk, x2 / 2):+.6f}   Schranke {kk + 1}")


def T_n_lattice(f, tau, n, k):
    s = 0j
    for aa in range(1, n + 1):
        if n % aa:
            continue
        dd = n // aa
        for bb in range(dd):
            s += dd ** (-k) * f((aa * tau + bb) / dd)
    return n ** (k - 1) * s


taut = 0.1 + 1.3j
print(f"Hecke als Gitter-Mittel (Untergitter vom Index n) bei tau = {taut}: T_n Delta / Delta")
for n in [2, 3, 4, 5, 6, 7]:
    val = T_n_lattice(Delta, taut, n, 12) / Delta(taut)
    print(f"   n = {n}:  {val.real:+.6f}{val.imag:+.1e}i   tau(n) = {TAU[n]}")


def hnf_count(n, box):
    seen = set()
    rng = range(-box, box + 1)
    for v1 in itertools.product(rng, rng):
        for v2 in itertools.product(rng, rng):
            det = v1[0] * v2[1] - v1[1] * v2[0]
            if abs(det) != n:
                continue
            # HNF: Basis (a,b),(0,d); a = ggT der ersten Komponenten
            a_ = math.gcd(v1[0], v2[0])
            if a_ == 0:
                continue
            # erweiterter ggT
            def eg(x, y):
                if y == 0:
                    return (abs(x), 1 if x >= 0 else -1, 0)
                g_, s_, t_ = eg(y, x % y)
                return (g_, t_, s_ - (x // y) * t_)
            g_, s_, t_ = eg(v1[0], v2[0])
            b0 = s_ * v1[1] + t_ * v2[1]
            dd = n // a_
            seen.add((a_, b0 % dd, dd))
    return len(seen)


print("Zahl der Untergitter vom Index n in Z^2 (Brute Force) gegen sigma_1(n):")
print("   ", [(n, hnf_count(n, n), SIG1[n]) for n in range(1, 9)])

# ==========================================================================
hdr("4  Rogers-Ramanujan-Kettenbruch R(q)")
# ==========================================================================


def K_cf(q, sign=+1, start=1, depth=None):
    """1 + s q^start/(1 + s q^(start+1)/(1 + ...))  (bottom-up)"""
    if depth is None:
        depth = int(80 / max(1e-12, -math.log(abs(q)))) + 60
    x = 1.0
    for n in range(depth, start - 1, -1):
        x = 1 + sign * q**n / x
    return x


def R_cf(q):
    q5 = math.copysign(abs(q) ** 0.2, q)
    return q5 / K_cf(q)


def R_prod_tau(tau):
    q = nome(tau)
    p = 1 + 0j
    n = 1
    while True:
        p *= (1 - q**(5 * n - 1)) * (1 - q**(5 * n - 4)) / ((1 - q**(5 * n - 2)) * (1 - q**(5 * n - 3)))
        if abs(q**(5 * n - 4)) < 1e-22:
            break
        n += 1
    return cmath.exp(2j * PI * tau / 5) * p


for t in [1.0, 0.5, 0.2]:
    q = math.exp(-2 * PI * t)
    print(f"   q = e^(-2 pi {t}):  Kettenbruch {R_cf(q):.15f}   Produkt {R_prod_tau(1j * t).real:.15f}")
v1_ = math.sqrt((5 + math.sqrt(5)) / 2) - PHI
print(f"R(e^-2pi) = {R_cf(math.exp(-2 * PI)):.15f}   sqrt((5+sqrt5)/2) - phi = {v1_:.15f}")
v2_ = math.sqrt((5 - math.sqrt(5)) / 2) - (math.sqrt(5) - 1) / 2
print(f"R(-e^-pi) = {R_cf(-math.exp(-PI)):.15f}   sqrt((5-sqrt5)/2) - 1/phi = {v2_:.15f}   (Betrag gleich: {abs(abs(R_cf(-math.exp(-PI))) - v2_):.1e})")
taui = 0.1 + 0.8j
Rv = R_prod_tau(taui)
lhs = Rv**-5 - 11 - Rv**5
rhs = (eta(taui) / eta(5 * taui)) ** 6
print(f"Ikosaeder-Identitaet R^-5 - 11 - R^5 = (eta(tau)/eta(5tau))^6 bei tau={taui}: rel. Fehler {rel(lhs, rhs):.1e}")
for t in [0.5, 0.2, 0.1]:
    Rt = R_cf(math.exp(-2 * PI * t))
    Rs = R_cf(math.exp(-2 * PI / t))
    print(f"   S-Transformation t={t}: R(e^-2pi t) = {Rt:.14f}   (1-phi R')/(phi+R') = {(1 - PHI * Rs) / (PHI + Rs):.14f}")
print("Grenzwert q -> 1:")
for eps in [1.0, 0.1, 0.01, 0.001]:
    print(f"   q = e^-{eps}:  R = {R_cf(math.exp(-eps)):.12f}   1/phi - R = {1 / PHI - R_cf(math.exp(-eps)):.3e}   Vorhersage (1+phi^-2) e^(-4pi^2/(5 eps)) = {(1 + PHI**-2) * math.exp(-4 * PI**2 / (5 * eps)):.3e}")
print(f"   1/phi = {1 / PHI:.12f};  (1/phi)^-5 - 11 - (1/phi)^5 = {PHI**5 - 11 - PHI**-5:.1e};  phi^5 = {PHI**5:.10f} = (11+5 sqrt5)/2 = {(11 + 5 * math.sqrt(5)) / 2:.10f}")
print("Gewichtete Fibonacci-Rekursion x_n = x_(n+1) + q^(n+1) x_(n+2): Ebenenverhaeltnis x_0/x_1 = K(q) = q^(1/5)/R(q):")
for q in [0.5, 0.9, 0.99, 0.999]:
    print(f"   q = {q}:  K = {K_cf(q):.10f}   (phi = {PHI:.10f})")
print("Diophantische Guete n^2 |theta - m/n| der besten Naeherungen (Hurwitz: Grenze 1/sqrt5 = %.6f fuer phi):" % (1 / math.sqrt(5)))
F = [1, 1]
while len(F) < 40:
    F.append(F[-1] + F[-2])
print("   phi:", [f"{F[n]}^2|phi-{F[n+1]}/{F[n]}| = {F[n]**2 * abs(PHI - F[n + 1] / F[n]):.6f}" for n in (5, 10, 20)])
for (m, n) in [(22, 7), (333, 106), (355, 113), (103993, 33102)]:
    print(f"   pi: n = {n:6d}, m/n = {m}/{n}:  n^2|pi - m/n| = {n * n * abs(PI - m / n):.6f}")

# --- Uebergangsfunktionen aus verschachtelter Rueckkopplung ----------------


def G_RR(a, q, tol=1e-18):
    s = 0.0
    poch = 1.0
    n = 0
    while True:
        if n > 0:
            poch *= (1 - q**n)
        t = a**n * q**(n * n) / poch
        s += t
        if n > 5 and abs(t) < tol * max(1.0, abs(s)):
            return s
        n += 1


qt = 0.3
cf_minus = 1 / K_cf(qt, sign=-1)
print(f"\nRamanujan: 1/(1+aq/(1+aq^2/...)) = G(aq)/G(a), G(a)=sum a^n q^(n^2)/(q;q)_n.  Test a=-1, q={qt}: KB {cf_minus:.14f}  Quotient {G_RR(-qt, qt) / G_RR(-1, qt):.14f}")
lo, hi = 0.3, 0.9
print(f"   G(-1;q): q=0.3 -> {G_RR(-1, 0.3):+.4f},  q=0.9 -> {G_RR(-1, 0.9):+.4f}")
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if G_RR(-1, mid) > 0:
        lo = mid
    else:
        hi = mid
QC = 0.5 * (lo + hi)
YC = math.log(QC) ** 2
print(f"   Pol der voll positiven Verschachtelung: G(-1;q_c)=0 bei q_c = {QC:.10f}  ->  y_c = (ln q_c)^2 = {YC:.6f}")


def nu_kb_minus(y):
    q = math.exp(-math.sqrt(y))
    return G_RR(-q, q) / G_RR(-1, q)


def nu_kb_plus(y):
    q = math.exp(-math.sqrt(y))
    Cq = G_RR(q * q, q) / G_RR(q, q)
    return 1 / (1 - q * Cq)


qq = 0.4
print(f"   Test KB+: 1/(1-q/(1+q^2/(1+q^3/...))) bei q={qq}: KB {1 / (1 - qq / K_cf(qq, +1, start=2)):.14f}   Quotientenform {1 / (1 - qq * G_RR(qq * qq, qq) / G_RR(qq, qq)):.14f}")
hd = 1e-7
Gp = (G_RR(-1, QC + hd) - G_RR(-1, QC - hd)) / (2 * hd)
Ares = G_RR(-QC, QC) / (-Gp)
SC = QC / Ares
print(f"   Residuum am Pol: nu ~ A/(q_c - q), A = {Ares:.8f};  BTFR-Normierung verlangt s = q_c/A = {SC:.8f}")


def nu_kb_star(y):
    q = QC * math.exp(-math.sqrt(y) / SC)
    return G_RR(-q, q) / G_RR(-1, q)


print(f"\n   {'y':>7s} {'nu_RAR':>9s} {'KB-':>9s} {'dex':>7s} {'KB+':>8s} {'dex':>7s} {'KB*':>9s} {'dex':>7s}")
for y in YGRID:
    r = nu_rar(y)
    km = nu_kb_minus(y) if y > YC else float('nan')
    kp = nu_kb_plus(y)
    ks = nu_kb_star(y)
    dm = LOG10(km / r) if y > YC and km > 0 else float('nan')
    print(f"   {y:7.3f} {r:9.4f} {km:9.4f} {dm:+7.3f} {kp:8.4f} {LOG10(kp / r):+7.3f} {ks:9.4f} {LOG10(ks / r):+7.3f}")
mx_p = max(abs(LOG10(nu_kb_plus(y) / nu_rar(y))) for y in YDENSE)
mx_s = max(abs(LOG10(nu_kb_star(y) / nu_rar(y))) for y in YDENSE)
mx_m = max(abs(LOG10(nu_kb_minus(y) / nu_rar(y))) for y in YDENSE if y > 1.5 * YC)
print(f"   max |dex| auf y in [0.01,100]:  KB+ {mx_p:.3f}   KB* {mx_s:.3f}   KB- (nur y > 1.5 y_c) {mx_m:.3f}")
print(f"   KB+ im tiefen Regime -> 1/(1-1/phi) = phi^2 = {PHI**2:.6f};  nu_KB+(1e-4) = {nu_kb_plus(1e-4):.6f}")
print(f"   KB* tief: nu*sqrt(y) bei y=1e-4: {nu_kb_star(1e-4) * 1e-2:.5f}, y=1e-6: {nu_kb_star(1e-6) * 1e-3:.5f} (BTFR verlangt 1)")

# ==========================================================================
hdr("5  Rogers-Ramanujan-Identitaeten: Ausschlussregeln")
# ==========================================================================


def count_parts_allowed(N, allowed):
    a = np.zeros(N + 1)
    a[0] = 1.0
    for m in range(1, N + 1):
        if allowed(m):
            for r in range(m):
                a[r::m] = np.cumsum(a[r::m])
    return a


def count_parts_int(N, allowed):
    a = [1] + [0] * N
    for m in range(1, N + 1):
        if allowed(m):
            for i in range(m, N + 1):
                a[i] += a[i - m]
    return a


NN = 60
prodG = count_parts_int(NN, lambda m: m % 5 in (1, 4))
prodH = count_parts_int(NN, lambda m: m % 5 in (2, 3))


def sum_side(NN, shift):
    tot = [0] * (NN + 1)
    n = 0
    while n * n + shift * n <= NN:
        inv = [1] + [0] * NN
        for kk in range(1, n + 1):
            for i in range(kk, NN + 1):
                inv[i] += inv[i - kk]
        e = n * n + shift * n
        for i in range(0, NN + 1 - e):
            tot[i + e] += inv[i]
        n += 1
    return tot


sG, sH = sum_side(NN, 0), sum_side(NN, 1)
print(f"G: sum q^(n^2)/(q;q)_n = prod 1/((1-q^(5n+1))(1-q^(5n+4))): Koeffizienten bis q^{NN} gleich: {sG == prodG}")
print(f"H: sum q^(n^2+n)/(q;q)_n = prod 1/((1-q^(5n+2))(1-q^(5n+3))): Koeffizienten bis q^{NN} gleich: {sH == prodH}")
from functools import lru_cache


@lru_cache(maxsize=None)
def gap2(n, m):
    if n == 0:
        return 1
    return sum(gap2(n - p, p + 2) for p in range(m, n + 1))


print("Partitionen mit Abstand >= 2 (direkt gezaehlt) gegen Produktseite, n = 0..25:",
      all(gap2(n, 1) == prodG[n] for n in range(26)), " | mit Teilen >= 2:", all(gap2(n, 2) == prodH[n] for n in range(26)))
print("   n = 10..15:", [gap2(n, 1) for n in range(10, 16)])
for q in [0.5, 0.9]:
    s = G_RR(1, q)
    p = 1.0
    n = 0
    while True:
        f1, f4 = q**(5 * n + 1), q**(5 * n + 4)
        p /= (1 - f1) * (1 - f4)
        if f1 < 1e-18:
            break
        n += 1
    print(f"   q = {q}: Summe {s:.14f}  Produkt {p:.14f}")
print("Harte Kerne (keine zwei benachbarten Plaetze besetzt) auf N Plaetzen = F_(N+2):")
fib_ok = True
for Nn in range(1, 17):
    cntb = sum(1 for m in range(1 << Nn) if not (m & (m >> 1)))
    fib_ok &= (cntb == F[Nn + 1])
print(f"   N = 1..16 geprueft: {fib_ok};  Entropie/Platz ln phi = {math.log(PHI):.6f} statt ln 2 = {math.log(2):.6f};  Verlust {math.log(2) - math.log(PHI):.6f} nats ({(1 - math.log(PHI) / math.log(2)) * 100:.1f} %)")
print(f"Hard-Hexagon kritische Aktivitaet z_c = (11+5 sqrt5)/2 = {(11 + 5 * math.sqrt(5)) / 2:.8f} = phi^5 = {PHI**5:.8f}")

# Andrews-Gordon k=3, i=3 (Brute Force)


def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest


ok = True
AG3 = count_parts_int(30, lambda m: m % 7 not in (0, 3, 4))
for n in range(0, 26):
    cntb = 0
    for pt in partitions(n):
        fr = [0] * (n + 3)
        for p in pt:
            fr[p] += 1
        if fr[1] <= 2 and all(fr[jj] + fr[jj + 1] <= 2 for jj in range(1, n + 1)):
            cntb += 1
    ok &= (cntb == AG3[n])
print(f"Andrews-Gordon k=3: (f_1<=2, f_j+f_(j+1)<=2) = Teile nicht 0,+-3 mod 7, n<=25: {ok}")
print("Asymptotik eps*ln G(e^-eps) -> pi^2/15 = %.8f  (frei: pi^2/6 = %.8f):" % (PI**2 / 15, PI**2 / 6))
for eps in [0.1, 0.01, 0.001]:
    n = np.arange(0, int(60 / eps / 5) + 10)
    lg = -np.sum(np.log1p(-np.exp(-(5 * n + 1) * eps))) - np.sum(np.log1p(-np.exp(-(5 * n + 4) * eps)))
    print(f"   eps = {eps}:  eps ln G = {eps * lg:.8f}")

# ==========================================================================
hdr("6  Mock-Theta-Funktion f(q), Schatten und Vervollstaendigung")
# ==========================================================================
NF = 400
ser = [0] * (NF + 1)
n = 0
while n * n <= NF:
    inv = [1] + [0] * NF
    for kk in range(1, n + 1):
        for _ in range(2):                         # 1/(1+q^k)^2
            for i in range(kk, NF + 1):
                inv[i] -= inv[i - kk]
    for i in range(0, NF + 1 - n * n):
        ser[i + n * n] += inv[i]
    n += 1
print("f(q) = sum q^(n^2)/((-q;q)_n)^2 = ", " ".join(f"{c:+d}q^{i}" for i, c in enumerate(ser[:13])))
print("Koeffizienten-Asymptotik (Dragonette/Andrews) a(n) ~ (-1)^(n-1) e^(pi sqrt(n/6-1/144)) / (2 sqrt(n-1/24)):")
for n in [50, 100, 200, 400]:
    pred = (-1) ** (n - 1) * math.exp(PI * math.sqrt(n / 6 - 1 / 144)) / (2 * math.sqrt(n - 1 / 24))
    print(f"   n = {n}:  a(n) = {ser[n]}   Verhaeltnis {ser[n] / pred:.6f}")


def f_b_dec(eps, prec=90):
    getcontext().prec = prec
    q = -(Dm(-eps).exp())
    tol = Dm(10) ** (-(prec - 20))
    s = Dm(0)
    den = Dm(1)
    n = 0
    while True:
        if n > 0:
            den *= (1 + q**n)
        t = q**(n * n) / (den * den)
        s += t
        if n > 5 and abs(t) < tol:
            break
        n += 1
    b = Dm(1)
    n = 1
    while True:
        qn = q**n
        b *= (1 - qn) / ((1 + qn) * (1 + qn))
        if abs(qn) < tol:
            break
        n += 1
    return s, b


print("Ramanujans Grenzwert q -> -1 radial (q = -e^-eps):  f(q) und b(q) = (q;q)_inf/(-q;q)_inf^2 explodieren, f + b bleibt endlich:")
fb_prev = None
for eps in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
    fv, bv = f_b_dec(eps)
    fb = float(fv + bv)
    extra = ""
    if fb_prev is not None:
        e0, v0 = fb_prev
        extra = f"   lineare Extrapolation eps->0: {v0 + (fb - v0) * (0 - e0) / (eps - e0):.5f}"
    fb_prev = (eps, fb)
    print(f"   eps = {eps:4}:  f = {float(fv):+.6e}   b = {float(bv):+.6e}   f + b = {fb:+.10f}{extra}")

# E2* (einfachstes Beispiel)
t_ = 0.3 + 1.1j
E2s = lambda z: E2(z) - 3 / (PI * z.imag)
print(f"\nE2* = E2 - 3/(pi Im tau): E2*(-1/tau)/(tau^2 E2*(tau)) = {E2s(-1 / t_) / (t_**2 * E2s(t_)):.12f}  (holomorpher Teil E2 allein: {E2(-1 / t_) / (t_**2 * E2(t_)):.6f})")

# Zwegers-Vervollstaendigung von f


def f_mock(q):
    s = 0j
    den = 1 + 0j
    n = 0
    while True:
        if n > 0:
            den *= (1 + q**n)
        t = q**(n * n) / den**2
        s += t
        if n > 4 and abs(q**(n * n)) < 1e-24:
            return s
        n += 1


def omega_mock(x):
    s = 0j
    den = 1 + 0j
    n = 0
    while True:
        den *= (1 - x**(2 * n + 1))
        t = x**(2 * n * (n + 1)) / den**2
        s += t
        if n > 4 and abs(x**(2 * n * (n + 1))) < 1e-24:
            return s
        n += 1


def Fvec(tau):
    q = nome(tau)
    qh = cmath.exp(1j * PI * tau)
    return np.array([cmath.exp(-2j * PI * tau / 24) * f_mock(q),
                     2 * cmath.exp(2j * PI * tau / 3) * omega_mock(qh),
                     2 * cmath.exp(2j * PI * tau / 3) * omega_mock(-qh)])


def eichler(tau, coef, shift, alt):
    """int_{-conj tau}^{i inf} g(z)/sqrt(-i(z+tau)) dz,
       g(z) = sum_n c_n exp(3 pi i (n+shift)^2 z), c_n = coef*(+-1)^n (n+shift)"""
    v = tau.imag
    s = 0j
    for n in range(-60, 61):
        a = math.sqrt(3) * (n + shift)
        x = abs(a) * math.sqrt(2 * PI * v)
        if x > 26:
            continue
        cn = coef * ((-1) ** n if alt else 1) * (n + shift)
        s += cn * (1j / abs(a)) * math.erfc(x) * cmath.exp(-1j * PI * a * a * tau)
    return s


def Rvec(tau):
    I0 = eichler(tau, 1.0, 1 / 3, True)       # g0
    I1 = eichler(tau, -1.0, 1 / 6, False)     # g1
    I2 = eichler(tau, 1.0, 1 / 3, False)      # g2
    return 2j * math.sqrt(3) * np.array([I1, I0, -I2])


Mz = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
z24, z3 = cmath.exp(2j * PI / 24), cmath.exp(2j * PI / 3)
print("Zwegers: H = F - R,  F = (q^(-1/24) f(q), 2q^(1/3) w(q^1/2), 2q^(1/3) w(-q^1/2)),  R = Eichler-Integral des Schattens (g1,g0,-g2)")
for tau in [0.2 + 0.9j, -0.35 + 0.8j, 0.05 + 1.3j]:
    Hv = Fvec(tau) - Rvec(tau)
    HS = Fvec(-1 / tau) - Rvec(-1 / tau)
    rhsS = cmath.sqrt(-1j * tau) * (Mz @ Hv)
    HT = Fvec(tau + 1) - Rvec(tau + 1)
    rhsT = np.array([[z24**-1, 0, 0], [0, 0, z3], [0, z3, 0]]) @ Hv
    FS = Fvec(-1 / tau)
    rhsF = cmath.sqrt(-1j * tau) * (Mz @ Fvec(tau))
    print(f"   tau = {tau}:  S: max|H(-1/tau) - sqrt(-i tau) M H| / |H| = {np.max(abs(HS - rhsS)) / np.max(abs(HS)):.1e}"
          f"   T: {np.max(abs(HT - rhsT)) / np.max(abs(HT)):.1e}   nur holomorph F: {np.max(abs(FS - rhsF)) / np.max(abs(FS)):.1e}")
print("Anteil des Schattens |R_0|/|F_0| entlang tau = i t  (t = sqrt(y)/(2 pi) bei FTM-Kopplung):")
for t in [0.15, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]:
    Fv = Fvec(1j * t)
    Rv_ = Rvec(1j * t)
    print(f"   t = {t:5.2f} (y = {(2 * PI * t)**2:8.2f}):  |R0/F0| = {abs(Rv_[0] / Fv[0]):.3e}   exp(-pi t/6)/sqrt(t) = {math.exp(-PI * t / 6) / math.sqrt(t):.3e}")

# ==========================================================================
hdr("7  Dedekind-eta, Partitionen, 1+2+3+... -> -1/12 und q^(1/24)")
# ==========================================================================
for tau in [0.1 + 0.7j, 0.37 + 0.5j, -0.2 + 1.5j]:
    l = eta(-1 / tau)
    r = cmath.sqrt(-1j * tau) * eta(tau)
    qp = lambda z: eta(z) / cmath.exp(2j * PI * z / 24)
    print(f"   tau = {tau}:  eta(-1/tau)/(sqrt(-i tau) eta(tau)) = {l / r:.12f}   ohne q^(1/24): {qp(-1 / tau) / (cmath.sqrt(-1j * tau) * qp(tau)):.6f}")
tz = 0.1 + 0.7j
print(f"   Multiplikator: eta(tau+1)/eta(tau) = {eta(tz + 1) / eta(tz):.12f} = e^(2 pi i/24) = {cmath.exp(2j * PI / 24):.12f};  (eta^24)(tau+1)/(eta^24)(tau) = {(eta(tz + 1) / eta(tz))**24:.12f}")
PN = 10000
pp = [1] + [0] * PN
for n in range(1, PN + 1):
    s = 0
    kk = 1
    while True:
        g1 = kk * (3 * kk - 1) // 2
        if g1 > n:
            break
        sg = 1 if kk % 2 else -1
        s += sg * pp[n - g1]
        g2 = kk * (3 * kk + 1) // 2
        if g2 <= n:
            s += sg * pp[n - g2]
        kk += 1
    pp[n] = s
print("p(0..15) =", pp[:16])
print("Abgeschnittene Summe sum n e^(-n eps) = 1/eps^2 - 1/12 + eps^2/240 - ...:")
for eps in [0.5, 0.1, 0.01]:
    sm = math.exp(-eps) / (1 - math.exp(-eps)) ** 2
    print(f"   eps = {eps}:  sum - 1/eps^2 = {sm - 1 / eps**2:.10f}   (-1/12 = {-1 / 12:.10f})")
print(f"Grundverschiebung: (1/2) sum n -> (1/2)(-1/12) = {-1 / 24:.6f} = Exponent von q in 1/eta = q^(-1/24)/prod(1-q^n)")
HBAR = 1.0546e-34
for L in [1e-6, 1.0]:
    print(f"   Beispiel geschlossene Schleife L = {L} m (masseloses Feld, beide Umlaufrichtungen): E0 = -pi hbar c/(6L) = {-PI * HBAR * C / (6 * L):.3e} J")
print(f"FTM-Kopplung: 1/eta(i t) bei t = 1 (y = 4pi^2): {1 / eta(1j).real:.8f};   eta(i) = Gamma(1/4)/(2 pi^(3/4)) = {math.gamma(0.25) / (2 * PI**0.75):.12f}  vs  {eta(1j).real:.12f}")

# ==========================================================================
hdr("8  Hardy-Ramanujan: Zustandszahl und Entropie")
# ==========================================================================
for n in [10, 100, 1000, 10000]:
    hr = math.exp(PI * math.sqrt(2 * n / 3)) / (4 * n * math.sqrt(3))
    S = math.log(pp[n]) if n < 1000 else (len(str(pp[n])) - 1) * math.log(10) + math.log(int(str(pp[n])[:15]) / 1e14)
    print(f"   n = {n:5d}:  p(n)/HR = {pp[n] / hr if n < 1000 else math.exp(S - math.log(hr)):.6f}   S = ln p = {S:10.4f}   pi sqrt(2n/3) = {PI * math.sqrt(2 * n / 3):10.4f}   pi sqrt(2n/3) - ln(4 sqrt3 n) = {PI * math.sqrt(2 * n / 3) - math.log(4 * math.sqrt(3) * n):10.4f}")


def dedekind_s(h, k):
    s = 0.0
    for r in range(1, k):
        x = h * r / k
        s += (r / k) * (x - math.floor(x) - 0.5)
    return s


def rademacher(n, K):
    lam_ = math.sqrt(n - 1 / 24)
    tot = 0.0
    for k in range(1, K + 1):
        Ak = 0j
        for hh in range(k):
            if math.gcd(hh, k) == 1:
                Ak += cmath.exp(1j * PI * dedekind_s(hh, k) - 2j * PI * n * hh / k)
        Cc = PI * math.sqrt(2 / 3) / k
        deriv = (Cc * math.cosh(Cc * lam_) / lam_**2 - math.sinh(Cc * lam_) / lam_**3) / 2
        tot += Ak.real * math.sqrt(k) * deriv
    return tot / (PI * math.sqrt(2))


for n in [100, 200]:
    print(f"   Rademacher n = {n}: K=1: {rademacher(n, 1):.4f}  K=12: {rademacher(n, 12):.4f}   exakt p(n) = {pp[n]}")
print(f"Cardy: S = 2 pi sqrt(c n/6); c = 1 (freie Moden): 2 pi sqrt(1/6) = {2 * PI * math.sqrt(1 / 6):.10f} = pi sqrt(2/3) = {PI * math.sqrt(2 / 3):.10f}")
print("   Die Grundverschiebung c/24 bestimmt den Entropie-Koeffizienten: c = 24 * (Exponent der Grundverschiebung)")

# ==========================================================================
hdr("9  Ramanujans Master-Theorem und das Gedaechtnisintegral")
# ==========================================================================


def expsinh(f, h=0.004, umin=-4.5, umax=4.0):
    u = np.arange(umin, umax + h / 2, h)
    x = np.exp(PI / 2 * np.sinh(u))
    w = PI / 2 * np.cosh(u) * x * h
    return float(np.sum(w * f(x)))


print("(a) Potenzgesetz-Kern K(s) = (1+s)^-alpha,  phi(k) = (alpha)_k  ->  int s^(sig-1) K = Gamma(sig) Gamma(alpha-sig)/Gamma(alpha)")
al = 2.5
for sg in [0.5, 1.0, 1.7]:
    num = expsinh(lambda x: x**(sg - 1) * (1 + x) ** (-al))
    ex = math.gamma(sg) * math.gamma(al - sg) / math.gamma(al)
    print(f"   alpha = {al}, sigma = {sg}:  Quadratur {num:.12f}   Master-Theorem {ex:.12f}")
print("(b) K(s) = (1-e^-s)/s,  phi(k) = 1/(k+1)  ->  Gamma(sig)/(1-sig)")
for sg in [0.3, 0.6]:
    num = expsinh(lambda x: x**(sg - 2) * (-np.expm1(-x)))
    print(f"   sigma = {sg}:  Quadratur {num:.12f}   Master-Theorem {math.gamma(sg) / (1 - sg):.12f}")


def Li(s, q, N=None):
    if N is None:
        N = int(80 / -math.log(q)) + 10
    lq = math.log(q)
    return sum(math.exp(n * lq - s * math.log(n)) for n in range(1, N))


print("(c) Rueckkopplungskern K(s) = sum_n q^n e^(-n lambda s),  phi(k) = Li_(-k)(q)  ->  Gamma(sig) lambda^-sig Li_sig(q)")
qv = math.exp(-1.0)
for sg in [0.5, 1.0, 2.0]:
    num = expsinh(lambda x: x**(sg - 1) * qv * np.exp(-x) / (1 - qv * np.exp(-x)))
    print(f"   q = e^-1, sigma = {sg}:  Quadratur {num:.12f}   Gamma(sig) Li_sig(q) = {math.gamma(sg) * Li(sg, qv):.12f}")
print(f"   K(0) = Li_0(q) = q/(1-q) = {qv / (1 - qv):.10f} = nu_RAR(1) - 1 = {nu_rar(1.0) - 1:.10f};  int K ds = Li_1(q)/lambda = -ln(1-q)/lambda = {-math.log(1 - qv):.10f}")

gl_x, gl_w = np.polynomial.legendre.leggauss(24)


def integrate_panels(f, L, npan):
    edges = np.linspace(0, L, npan + 1)
    tot = 0.0
    for aa, bb in zip(edges[:-1], edges[1:]):
        xm = 0.5 * (bb - aa) * gl_x + 0.5 * (aa + bb)
        tot += 0.5 * (bb - aa) * np.sum(gl_w * f(xm))
    return tot


print("Kreisbahn-Gedaechtnisintegral I(Omega) = int_0^inf K(s) sin(Omega s) ds:")
print("   Master-Theorem -> I = (1/lambda) sum_m (-1)^m (Omega/lambda)^(2m+1) phi(-2m-2)")
lamb = 1.0
for (Om, qv) in [(0.5, math.exp(-1.0)), (0.3, math.exp(-0.2))]:
    Kf = lambda s: qv * np.exp(-lamb * s) / (1 - qv * np.exp(-lamb * s))
    quad = integrate_panels(lambda s: Kf(s) * np.sin(Om * s), 120.0 / lamb, 600)
    ser_ = sum((-1) ** m * (Om / lamb) ** (2 * m + 1) * Li(2 * m + 2, qv) for m in range(60)) / lamb
    closed = sum(qv**n * Om / (n * n * lamb * lamb + Om * Om) for n in range(1, 4000))
    print(f"   q = {qv:.4f}, Omega/lambda = {Om}:  Quadratur {quad:.12f}   Reihe {ser_:.12f}   sum_n q^n Omega/(n^2 lambda^2+Omega^2) = {closed:.12f}")
qv = math.exp(-1.0)
for ratio in [1e-1, 1e-2, 1e-3]:
    Om = 1.0
    lb = ratio * Om
    val = sum(qv**n * Om / (n * n * lb * lb + Om * Om) for n in range(1, 4000))
    print(f"   langes Gedaechtnis lambda/Omega = {ratio}:  Omega*I = {Om * val:.8f}  ->  q/(1-q) = nu_RAR - 1 = {qv / (1 - qv):.8f}")
print("   Potenzgesetz-Kern (1+s)^-3: nur Momente m < 2 existieren (phi(-sig) = Gamma(3-sig)/2 hat Pol bei sig = 3):")
for Om in [0.1, 0.01]:
    Lp = 4e4 if Om == 0.1 else 2e5
    quad = integrate_panels(lambda s: (1 + s) ** -3.0 * np.sin(Om * s), Lp, int(Lp / 4))
    print(f"   Omega = {Om}:  I/Omega = {quad / Om:.8f}   (erstes Moment phi(-2) = Gamma(1)/Gamma(3) = 0.5)")
print("Mellin-Momente des Delta-Kerns: Lambda(sig) = int_0^inf s^(sig-1) Delta(is) ds = (2 pi)^-sig Gamma(sig) L(Delta,sig)")


def Lam_int(sg, umin=-3.5, umax=3.5, N=1401):
    u = np.linspace(umin, umax, N)
    vals = np.array([math.exp(sg * uu) * Delta_it(math.exp(uu)) for uu in u])
    return float(np.trapezoid(vals, u)) if hasattr(np, "trapezoid") else float(np.trapz(vals, u))


def Lam_dir(sg):
    return (2 * PI) ** -sg * math.gamma(sg) * sum(TAU[n] / n**sg for n in range(1, NT + 1))


for sg in [8, 10, 12]:
    print(f"   sigma = {sg:2d}:  Integral {Lam_int(sg):.10e}   Dirichlet-Reihe (n<=1000) {Lam_dir(sg):.10e}")
print("   Funktionalgleichung Lambda(sig) = Lambda(12-sig), Integral konvergiert fuer alle sig (Delta = 0 an beiden Spitzen):")
for sg in [0, 2, 4, 6]:
    print(f"   Lambda({sg}) = {Lam_int(sg):.12e}   Lambda({12 - sg}) = {Lam_int(12 - sg):.12e}")

# ==========================================================================
hdr("10  Lambert-Reihe: Rueckkopplung ueber alle Harmonischen")
# ==========================================================================
q = 0.3
print(f"q = {q}:  sum sigma(n) q^n = {sum(SIG1[n] * q**n for n in range(1, 200)):.14f}")
print(f"         sum n q^n/(1-q^n) = {lam(q, 1).real:.14f}")
print(f"         (1-E2)/24         = {((1 - E2q(q)) / 24).real:.14f}")
print(f"         sum q^n/(1-q^n)   = {lam(q, 0).real:.14f} = sum d(n) q^n = {sum(divcount_list(200)[n] * q**n for n in range(1, 200)):.14f}  (NICHT (1-E2)/24)")
print("Einzelschleife: q/(1-q) = nu_RAR - 1 ist genau der Term n = 1 von sum_n q^n/(1-q^n); Harmonische n: q^n = exp(-sqrt(n^2 y))")


def nu_d(y):
    return 1 + lam(math.exp(-math.sqrt(y)), 0).real


def nu_sig(y):
    return 1 + lam(math.exp(-math.sqrt(y)), 1).real


def nu_P(y):
    t = y**0.25 / (2 * PI)
    return 1 + (1 - E2_it(t)) / (4 * PI**2)


def nu_Pstar(y):
    t = y**0.25 / (2 * PI)
    return nu_P(y) + 3 / (PI * t) / (4 * PI**2)


print(f"\n   {'y':>7s} {'nu_RAR':>9s} | {'nu_d':>9s} {'dex':>7s} | {'nu_sig':>10s} {'dex':>7s} | {'nu_P':>9s} {'dex':>7s} | {'nu_P*':>9s} {'dex':>7s}")
for y in YGRID:
    r = nu_rar(y)
    vals = [nu_d(y), nu_sig(y), nu_P(y), nu_Pstar(y)]
    print(f"   {y:7.3f} {r:9.4f} | {vals[0]:9.4f} {LOG10(vals[0] / r):+7.3f} | {vals[1]:10.4f} {LOG10(vals[1] / r):+7.3f} | {vals[2]:9.4f} {LOG10(vals[2] / r):+7.3f} | {vals[3]:9.4f} {LOG10(vals[3] / r):+7.3f}")
for nm, fn in [("nu_d", nu_d), ("nu_sig", nu_sig), ("nu_P", nu_P), ("nu_P*", nu_Pstar)]:
    devs = [LOG10(fn(y) / nu_rar(y)) for y in YDENSE]
    i = int(np.argmax(np.abs(devs)))
    print(f"   {nm:7s}: max |dex| auf [0.01,100] = {abs(devs[i]):.4f} bei y = {YDENSE[i]:.3f};  rms = {math.sqrt(sum(d * d for d in devs) / len(devs)):.4f}")
alg = lambda y: 1 + 1 / (4 * PI**2) + y**-0.5 - 3 / PI**2 * y**-0.25
print(f"   nu_P gegen algebraische Form 1 + 1/(4pi^2) + y^-1/2 - (3/pi^2) y^-1/4 auf [0.01,100]: max rel. Abweichung {max(rel(nu_P(y), alg(y)) for y in YDENSE):.1e}")
print("Tiefes Regime (Asymptotik):")
for y in [1e-4, 1e-6, 1e-8]:
    e = math.sqrt(y)
    print(f"   y = {y:.0e}:  nu_P*sqrt(y) = {nu_P(y) * e:.6f}   nu_P gegen 1+1/(4pi^2)+y^-1/2-(3/pi^2)y^-1/4: rel {rel(nu_P(y), 1 + 1 / (4 * PI**2) + 1 / e - 3 / PI**2 * y**-0.25):.1e}"
          f"   nu_sig*y = {nu_sig(y) * y:.6f} (pi^2/6 = {PI**2 / 6:.6f})   nu_d*sqrt(y)/(gamma+ln(1/sqrt y)) = {nu_d(y) * e / (EULER_GAMMA - math.log(e)):.6f}" if y >= 1e-6 else
          f"   y = {y:.0e}:  nu_P*sqrt(y) = {nu_P(y) * e:.6f}")
print("Newton-Regime (Sonnensystem), relative Zusatzbeschleunigung nu - 1:")
for nm, y in [("Erde", Y_EARTH), ("Saturn", Y_SATURN), ("Neptun", Y_NEPTUN)]:
    t = y**0.25 / (2 * PI)
    print(f"   {nm:7s} y = {y:.3e}:  nu_P - 1 = {nu_P(y) - 1:.3e}   (6/pi^2) e^(-y^1/4) = {6 / PI**2 * math.exp(-y**0.25):.3e}   nu_P* - 1 = {nu_Pstar(y) - 1:.3e}   nu_RAR - 1 = {nu_rar(y) - 1:.3e}")

# ==========================================================================
hdr("11  Ramanujans 1/pi-Reihe, singulaere Moduln, E2(i) = 3/pi")
# ==========================================================================
getcontext().prec = 110


def arctan_inv(x):
    x = Dm(x)
    x2 = x * x
    term = 1 / x
    s = term
    n = 1
    while True:
        term /= -x2
        t = term / (2 * n + 1)
        if abs(t) < Dm(10) ** -108:
            break
        s += t
        n += 1
    return s


PI_D = 4 * (4 * arctan_inv(5) - arctan_inv(239))
fact = math.factorial
sacc = Dm(0)
for kk in range(0, 8):
    sacc += Dm(fact(4 * kk) * (1103 + 26390 * kk)) / Dm(fact(kk) ** 4 * 396 ** (4 * kk))
    inv_pi = Dm(2) * Dm(2).sqrt() / Dm(9801) * sacc
    err = abs(1 / inv_pi - PI_D)
    print(f"   K = {kk}:  1/pi aus {kk + 1} Termen, Fehler in pi = {float(err):.3e}  ({-math.log10(float(err)):.1f} Stellen)")
e58 = (PI_D * Dm(58).sqrt()).exp()
print(f"   e^(pi sqrt58) = {e58:.12f}")
print(f"   396^4 - 104   = {Dm(396) ** 4 - 104}   Differenz = {float(Dm(396) ** 4 - 104 - e58):.6e}")
print(f"   E2(i) = {E2_it(1.0):.15f},  3/pi = {3 / PI:.15f};   E2*(i) = E2(i) - 3/pi = {E2_it(1.0) - 3 / PI:.1e}")
x = PI_D
cf = []
for _ in range(14):
    a_ = int(x)
    cf.append(a_)
    x = 1 / (x - a_)
print("   Kettenbruch pi =", cf, "   phi = [1;1,1,1,...]")

# ==========================================================================
hdr("12  Eigene Formeln")
# ==========================================================================
print("12a  Hecke-Skalenleiter: T_p E2 = (1+p) E2  (Gewicht 2, Gitter-Mittel)")
for p in [2, 3, 5, 7]:
    val = T_n_lattice(E2, taut, p, 2) / E2(taut)
    print(f"   p = {p}:  T_p E2 / E2 = {val.real:.12f}{val.imag:+.1e}i   sigma_1(p) = {1 + p}")
print("   FTM-Form (q = e^-sqrt y): l(y) = sum sigma(n) e^(-n sqrt y),  l^(p)(y) = sum_m sigma(pm) e^(-m sqrt y)")
print("   Leiter: l^(p)(y) = (1+p) l(y) - p l(p^2 y)")


def ell(y, p=1, N=20000):
    e = math.sqrt(y)
    s = 0.0
    for m in range(1, N):
        if p * m > 20000:
            break
        t = SIG1[p * m] * math.exp(-m * e)
        s += t
        if m > 10 and t < 1e-18 * s:
            break
    return s


for y in [0.3, 1.0, 5.0]:
    for p in [2, 3]:
        lhs = ell(y, p)
        rhs = (1 + p) * ell(y) - p * ell(p * p * y)
        print(f"   y = {y}, p = {p}:  l^(p) = {lhs:.12f}   (1+p) l(y) - p l(p^2 y) = {rhs:.12f}")
neg = sum(1 for p in primes if TAU[p] < 0)
print(f"   Delta-Leiter: Multiplikatoren tau(p), p <= 1000: {neg} von {len(primes)} negativ;  Eisenstein: sigma_(k-1)(p) > 0 fuer alle p")
print("   Delta-Leiter FTM-Form: sum_m tau(pm) q^m = tau(p) Delta(q) - p^11 Delta(q^p):")
for y in [1.0, 4.0]:
    qv = math.exp(-math.sqrt(y))
    for p in [2, 3]:
        lhs = sum(TAU[p * m] * qv**m for m in range(1, NT // p + 1))
        Dq = sum(TAU[m] * qv**m for m in range(1, NT + 1))
        Dqp = sum(TAU[m] * qv**(p * m) for m in range(1, NT + 1))
        print(f"   y = {y}, p = {p}:  links {lhs:.12e}   rechts {TAU[p] * Dq - p**11 * Dqp:.12e}")

print("\n12b  Vervollstaendigter Kern: Rest im Newton-Regime (Spitze), holomorph + Schatten")
print("   Eisenstein-Schatten (E2*): Zusatz 3/(pi Im tau) faellt nur wie Potenz;  Spitzenform-Schatten: erfc-Faktor, exponentiell")
for nm, y in [("Neptun", Y_NEPTUN), ("Saturn", Y_SATURN), ("Erde", Y_EARTH)]:
    tR = math.sqrt(y) / (2 * PI)
    tP = y**0.25 / (2 * PI)
    print(f"   {nm:7s}: Im tau (RAR) = {tR:9.2f}  -> E2*-Zusatz/24 = {3 / (PI * tR) / 24:.2e}, Schatten ~ exp(-pi t/6) = {math.exp(-PI * tR / 6):.1e};"
          f"   Im tau (y^1/4) = {tP:6.3f} -> E2*-Zusatz/(4pi^2) = {3 / (PI * tP) / (4 * PI**2):.2e}, exp(-pi t/6) = {math.exp(-PI * tP / 6):.2e}")
SS_BOUND = 1e-10
kap = [math.log(1 / SS_BOUND) / (2 * PI * (Y_SATURN**0.25 / (2 * PI))), math.log(1 / SS_BOUND) / (2 * PI * (math.sqrt(Y_SATURN) / (2 * PI)))]
print(f"   Annahme Sonnensystem-Grenze |nu-1| < {SS_BOUND:.0e} bei Saturn -> kleinster Exponent kappa des Schattens (Rest ~ e^(-2 pi kappa Im tau)):")
print(f"      Kopplung y^1/4: kappa > {kap[0]:.3f};   Kopplung sqrt(y): kappa > {kap[1]:.4f};  f-Schatten hat kappa = 1/12 = {1 / 12:.4f}")

print("\n12c  Ausschluss-Entropie: Andrews-Gordon-Klasse k  ->  c_eff = 1 - 3/(2k+1),  S(n) ~ 2 pi sqrt(c_eff n/6)")
NA = 4000
classes = [("frei (Bosonen)", 1.0, None, None),
           ("verschieden (Pauli)", 0.5, 2, [1]),
           ("RR k=2", 2 / 5, 5, [1, 4]),
           ("AG k=3", 4 / 7, 7, [1, 2, 5, 6]),
           ("AG k=4", 6 / 9, 9, [1, 2, 3, 6, 7, 8])]
for nm, ceff, m, res in classes:
    if m is None:
        arr = np.array([float(pp[n]) for n in range(NA + 1)])
        Cpref = None
    else:
        arr = count_parts_allowed(NA, lambda x, m=m, res=res: x % m in res)
        Cpref = 1.0
        for a_ in res:
            Cpref *= math.gamma(a_ / m) / math.sqrt(2 * PI)
    Acoef = PI**2 * ceff / 6
    outs = []
    for n in [500, 2000, 4000]:
        if Cpref is None:
            pred = math.exp(PI * math.sqrt(2 * n / 3)) / (4 * n * math.sqrt(3))
        else:
            pred = Cpref * Acoef**0.25 / (2 * math.sqrt(PI) * n**0.75) * math.exp(2 * math.sqrt(Acoef * n))
        outs.append(f"n={n}: Zahl/Vorhersage {arr[n] / pred:.4f}")
    print(f"   {nm:20s} c_eff = {ceff:.4f}  2pi sqrt(c/6) = {2 * PI * math.sqrt(ceff / 6):.5f}   " + "   ".join(outs))
n = 1000
print(f"   Entropieverlust durch RR-Ausschluss bei n = {n}: pi sqrt(2n/3)(1 - sqrt(2/5)) = {PI * math.sqrt(2 * n / 3) * (1 - math.sqrt(0.4)):.3f} nats (fuehrende Ordnung)")
print("\nEnde.")
