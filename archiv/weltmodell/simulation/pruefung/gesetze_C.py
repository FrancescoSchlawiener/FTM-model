# -*- coding: utf-8 -*-
"""
Gesetze Klasse C: Gesetze des Zustandsraums (Laplace/Lévy, Kac, Pólya/Watson, Bose/Wien, Kippstelle,
Fluktuation-Dissipation, Onsager, Ramanujan, Thetareihen, Barkhausen), durch Simulation oder zweiten Rechenweg.
Aufbau wie gesetze_A.py; Kennungen beginnen mit "C_". Die Modellseite ist die Formel des Kerns, wo es sie gibt;
Zufall nur über einen eigenen Generator mit festem Seed je Gesetz. Größen der Simulationen (Schrittzahl, Länge,
Verlust) sind Rechenpunkte der Gesetze, keine Aussagen des Modells.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, ive, lambertw

import kern as K
from kern.formeln import FORMELN

from .aussage import AUSSAGEN
from .gesetz import GESETZE, ident, sim, wahr

_g = GESETZE.neu
SEED = 29092026


# ------------------------------------------------------------------ C1 Rückkehr mit Verlust
RW_L, RW_P, RW_N = 12, 0.004, 40000          # Abstand, Verlust je Schritt, Irrfahrer (D = ½, Γ = p)


def _rueckkehr_irrfahrt():
    rng = np.random.default_rng(SEED)
    x = np.full(RW_N, RW_L)
    lebt = np.ones(RW_N, bool)
    treffer = np.zeros(RW_N, bool)
    while lebt.any():
        i = np.where(lebt)[0]
        x[i] += rng.choice((-1, 1), size=i.size)
        weg = rng.uniform(0.0, 1.0, i.size) < RW_P
        da = x[i] == 0
        treffer[i[da & ~weg]] = True
        lebt[i[da | weg]] = False
    return treffer.mean()


def _rueckkehr_kern():
    """q = e^{−√(Γ·τ_D)} des Kerns mit τ_D = L²/D (D = ½): a₀ = c/τ_D."""
    return float(K.rueckkehr(0.0, RW_P, 0.0, K.c / (RW_L**2 / 0.5)))


_g("C_rueckkehr_verlust", "C", "Rückkehr mit Verlust: q = e^{−L√(Γ/D)} (Irrfahrt, D = ½, Γ = p, L = 12)",
   lambda p: sim(_rueckkehr_kern, _rueckkehr_irrfahrt, 0.05), "Laplace-Transformierte der Erstrückkehrzeit (Lévy)",
   "I", "Wirbel")
_g("C_rueckkehr_diskret", "C", "diskret exakt: [(1 − √(1 − (1−p)²))/(1 − p)]^L gegen e^{−L√(2p)}",
   lambda p: ident(_rueckkehr_kern, lambda: ((1 - math.sqrt(1 - (1 - RW_P)**2)) / (1 - RW_P))**RW_L, 0.01),
   "Feller, Kap. XIV", "I", "Wirbel")

# ------------------------------------------------------------------ C2 Kac
KAC_N, KAC_LAEUFER = 15, 20000


def _kac_irrfahrt():
    rng = np.random.default_rng(SEED + 1)
    x = np.ones(KAC_LAEUFER, int)
    t = np.ones(KAC_LAEUFER)
    fertig = np.zeros(KAC_LAEUFER, bool)
    while not fertig.all():
        i = np.where(~fertig)[0]
        xn = x[i] + rng.choice((-1, 1), size=i.size)
        xn[xn > KAC_N] = KAC_N - 1
        x[i] = xn
        t[i] += 1
        fertig[i[xn == 0]] = True
    return t.mean()


_g("C_kac_halbgerade", "C", "Kac: mittlere Rückkehrzeit auf der bei N geschlossenen Halbgeraden = μ(X)/μ(0) = 2N (N = 15)",
   lambda p: sim(lambda: K.kac_rueckkehrzeit(KAC_N, 0.5), _kac_irrfahrt, 0.03), "Kac 1947", "Z", "transversal")


def _watson_integral():
    T = 1e5
    return quad(lambda t: ive(0, t / 3)**3, 0, T, limit=2000)[0] + 2 * (3 / (2 * math.pi))**1.5 / math.sqrt(T)


_g("C_watson", "C", "Pólya/Watson: Rückkehrsumme W₃ = ∫e^{−t}I₀(t/3)³dt = (√6/32π³)Γ(1/24)Γ(5/24)Γ(7/24)Γ(11/24)",
   lambda p: ident(_watson_integral, lambda: math.sqrt(6) / (32 * math.pi**3) * gamma(1 / 24) * gamma(5 / 24)
                   * gamma(7 / 24) * gamma(11 / 24), 2e-4), "Watson 1939; Pólya 1921", "Z", "transversal")

for _d in (1, 2, 3):
    _g(f"C_bose_leiter_{_d}", "C", f"Bose-Leiter ∫s^{_d}n_B ds = Γ({_d}+1)ζ({_d}+1)",
       lambda p, d=_d: ident(lambda: K.bose_leiter_integral(d), lambda: quad(lambda x: x**d / math.expm1(x), 0, 200)[0],
                             1e-7),
       "Planck 1900; Stefan-Boltzmann (d = 3)", "Z", "tangential")
_g("C_wien_verschiebung", "C", "Wien: Maximum von s³n_B bei s = 3 + W₀(−3e⁻³)",
   lambda p: ident(lambda: K.bose_stufe_max(3)[0], lambda: 3 + lambertw(-3 * math.exp(-3)).real, 1e-9),
   "Wien 1893 (Verschiebungsgesetz)", "Z", "tangential")

# ------------------------------------------------------------------ C5 Kippstelle, C6 FDT, C7 Onsager
FENSTER_B = (0.4, 0.1)


def _fenster_num(b):
    """Breite des Bereichs von a mit drei reellen Wurzeln von −x³ + bx + a (Eigenwerte der Begleitmatrizen)."""
    aa = np.linspace(-0.2, 0.2, 40001)          # festes Raster, unabhängig von b
    C = np.zeros((aa.size, 3, 3))
    C[:, 0, 1] = b
    C[:, 0, 2] = aa
    C[:, 1, 0] = C[:, 2, 1] = 1.0   # Begleitmatrix von x³ − bx − a
    ev = np.linalg.eigvals(C)
    drei = np.sum(np.abs(ev.imag) < 1e-7, axis=1) == 3
    w = aa[drei]
    return w.max() - w.min()


_g("C_kippstelle_fenster", "C", "Kippstelle: Hysteresefenster ∝ b^{3/2} (Verhältnis b = 0,4 zu 0,1 = 8)",
   lambda p: sim(lambda: (FENSTER_B[0] / FENSTER_B[1])**1.5, lambda: _fenster_num(FENSTER_B[0]) / _fenster_num(FENSTER_B[1]),
                 2e-3), "Thom 1972 (Spitzenkatastrophe)", "Z", "parametrisch")

FDT_KT, FDT_K, FDT_N, FDT_SCHRITTE = 1.0, 2.0, 2000, 20000


def _fdt(gam, seed):
    rng = np.random.default_rng(seed)
    dt = 0.02 * gam / FDT_K
    x = np.zeros(FDT_N)
    acc, n = 0.0, 0
    for i in range(FDT_SCHRITTE):
        x += -FDT_K / gam * x * dt + math.sqrt(2 * FDT_KT / gam * dt) * rng.standard_normal(FDT_N)
        if i > FDT_SCHRITTE // 2 and i % 10 == 0:
            acc += np.mean(x * x)
            n += 1
    return acc / n


for _gam, _k in ((1.0, "C_fdt_daempfung_1"), (10.0, "C_fdt_daempfung_10")):
    _g(_k, "C", f"Fluktuation-Dissipation ⟨x²⟩ = kT/k, Dämpfung γ = {_gam:g} (Langevin)",
       lambda p, gam=_gam: sim(lambda: FDT_KT / FDT_K, lambda: _fdt(gam, SEED + int(gam)), 0.02),
       "Einstein 1905; Kubo 1966", "Z", "transversal")


def _onsager():
    rng = np.random.default_rng(SEED + 2)
    Ls = rng.normal(size=(4, 4))
    Ls = Ls @ Ls.T
    La = rng.normal(size=(4, 4))
    La = La - La.T
    X = rng.normal(size=4)
    return X @ (Ls + La) @ X, X @ Ls @ X


_g("C_onsager_casimir", "C", "Onsager-Casimir: Entropieerzeugung X·(L_s + L_a)·X = X·L_s·X (nur der symmetrische Teil)",
   lambda p: ident(lambda: _onsager()[0], lambda: _onsager()[1], 1e-12), "Onsager 1931; Casimir 1945", "Z",
   "transversal")


# ------------------------------------------------------------------ C8 Ramanujan, C9 Thetareihen
def _ramanujan():
    N = 14
    poly = [1] + [0] * N
    for n in range(1, N + 1):
        for _ in range(24):
            neu = poly[:]
            for k in range(n, N + 1):
                neu[k] -= poly[k - n]
            poly = neu
    tau = {n: poly[n - 1] for n in range(1, N + 1)}

    def sig11(n):
        return sum(d**11 for d in range(1, n + 1) if n % d == 0)
    ok = (tau[2] == -24 and tau[3] == 252 and tau[4] == -1472 and tau[5] == 4830 and tau[6] == tau[2] * tau[3]
          and tau[4] == tau[2]**2 - 2**11 and all((tau[n] - sig11(n)) % 691 == 0 for n in range(1, N + 1))
          and all(abs(tau[q]) <= 2 * q**5.5 for q in (2, 3, 5, 7, 11, 13)))
    return wahr(ok, f"τ(2) = {tau[2]}")


_g("C_ramanujan_tau", "C", "Ramanujan τ(n): Werte, Multiplikativität, Hecke, Kongruenz mod 691, Deligne-Schranke",
   lambda p: _ramanujan(), "Ramanujan 1916; Mordell 1917; Deligne 1974", "Z", "tangential")


def _bilder():
    def zaehl(n, sechs):
        return sum(1 for m in range(-6, 7) for k in range(-6, 7) if (m * m + m * k + k * k if sechs else m * m + k * k) == n)

    def teiler(n, mod, r):
        return sum(1 for d in range(1, n + 1) if n % d == 0 and d % mod == r)
    ok = all(zaehl(n, False) == 4 * (teiler(n, 4, 1) - teiler(n, 4, 3))
             and zaehl(n, True) == 6 * (teiler(n, 3, 1) - teiler(n, 3, 2)) for n in range(1, 26))
    return wahr(ok, "n = 1…25; quadratisch 4, 4, 0, 4, 8, sechseckig 6, 0, 6, 6, 0 für n = 1…5")


_g("C_thetareihen_bilder", "C",
   "Bilder im Kaleidoskop (Zelle als Torus): r₂(n) = 4(d₁ − d₃), sechseckig 6(d₁ − d₂) mod 3, n ≤ 25",
   lambda p: _bilder(), "Jacobi 1829; Eisenstein", "G", "Hintergrund")


# ------------------------------------------------------------------ C10 Barkhausen, C11 Kac für eine Zerlegung
VDP_MU, VDP_T, VDP_DT = 0.05, 600.0, 0.005


def _vdp():
    x, v, amp = 0.1, 0.0, 0.0
    n = int(VDP_T / VDP_DT)
    ab = int((VDP_T - 20) / VDP_DT)
    for i in range(n):
        v += (VDP_MU * (1 - x * x) * v - x) * VDP_DT
        x += v * VDP_DT
        if i > ab and abs(x) > amp:
            amp = abs(x)
    return amp


_g("C_barkhausen", "C", "Barkhausen: Amplitude des Grenzzyklus A* = 2 (van der Pol, μ = 0,05; Schleifenverstärkung 1 im Mittel)",
   lambda p: sim(lambda: 2.0, _vdp, 0.01), "van der Pol 1926; Barkhausen 1935", "Z", "transversal")

KAC_KETTEN, KAC_SCHRITTE = (3, 5, 8), 100000


def _kac_ketten():
    """Σ 1/⟨T_i⟩ über die Zustände dreier zufälliger geschlossener Markov-Ketten (3, 5, 8 Zustände), gezählt."""
    aus = []
    for j, n in enumerate(KAC_KETTEN):
        rng = np.random.default_rng(SEED + 10 + j)
        P = rng.uniform(0.0, 1.0, (n, n))
        P /= P.sum(1, keepdims=True)
        kum = np.cumsum(P, axis=1).tolist()
        u = rng.uniform(0.0, 1.0, KAC_SCHRITTE).tolist()
        zust, zuletzt, summe, zahl = 0, [-1] * n, [0] * n, [0] * n
        for t in range(KAC_SCHRITTE):
            if zuletzt[zust] >= 0:
                summe[zust] += t - zuletzt[zust]
                zahl[zust] += 1
            zuletzt[zust] = t
            reihe, x, k = kum[zust], u[t], 0
            while k < n - 1 and x > reihe[k]:
                k += 1
            zust = k
        aus.append(sum(zahl[i] / summe[i] for i in range(n)))
    return np.array(aus)


AUSSAGEN.neu("C_kac_budget", FORMELN["kac_budget"].modell_ort, "Beziehung")
_g("C_kac_budget_ketten", "C", "Kac für jede Zerlegung: Σ 1/⟨T_i⟩ = 1 (drei geschlossene Ketten, 3, 5, 8 Zustände)",
   lambda p: sim(lambda: K.kac_budget(1.0, [0.5, 0.5]), _kac_ketten, 0.01), "Kac 1947; Baustein 06", "Z", "transversal",
   bezug="C_kac_budget")
