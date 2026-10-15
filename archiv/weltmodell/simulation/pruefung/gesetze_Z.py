# -*- coding: utf-8 -*-
"""
Gesetze Klasse Z: Zustandsraum des Autors, Tafeln 1–6, durch Simulation oder zweiten Rechenweg. Aufbau wie
gesetze_A.py; Kennungen beginnen mit "Z_". Die Quelle nennt die Tafel:
  Tafel 1  Übersicht und Wirkfelder: die Schleifen Gravitation, Wärme, Strömung, Akustik/Abstrahlung, Chemie/Alterung,
           Ermüdung, Kriechen, Rauschen (05, 06)
  Tafel 2  Generatoren von Flussrichtungen: drei Rollen, GENERIC (05)
  Tafel 3  geschlossene Wirklogik: Maxwell, Stoffgesetze, Lorentz, Poynting, Wärme (05)
  Tafel 4  Topologie: Schnitt, Rand, Schließung, Faserbündel (05)
  Tafel 5  Bindungsoperator: L_eff = 𝒫L, Gedächtniskern, Quelle gegen Operator (05)
  Tafel 6  Trennbeobachtungen, Tests 1–6 (05)
dazu die Kopplungsgesetze aus 06 (Einrasten, Zungen, Kuramoto, Koexistenz, Amplitudentod, Kippstelle, Dreieck).
Die Triadenrolle ist auf dieser Ebene die Rolle im Zustandsraum (tangential, transversal, parametrisch).
Zufall nur über einen eigenen Generator mit festem Seed je Gesetz; Größen der Simulationen sind Rechenpunkte.
"""
import functools
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

import kern as K
from kern.eingabe import EINGABEN
from kern.messwert import MESSWERTE

from kern.zustandsraum import TRIADE_GLIEDER

from .aussage import AUSSAGEN
from .gesetz import GESETZE, ident, sim, wahr
from .sonnenort import SONNENORT, gezeit

_g = GESETZE.neu
SEED = 29092026
T1, T2, T3, T4, T5, T6 = (f"Autor, Tafel {n}" for n in range(1, 7))
# Die Quelle des Autors liegt nicht im Repo: Die Zuordnung der Tafel 5 zum Bindungsoperator ist gewählt, bis der Nutzer sie
# bestätigt. Alle Gesetze der Tafel 5 tragen [W] mit diesem Text.
W_TAFEL5 = "Tafel 5 = Bindungsoperator (Zuordnung für Gedächtnis, Maxwell-Element und Operator, Nutzer bestätigt noch; Poisson-Inversion und Additivität bestätigt)"


def _E(k):
    return EINGABEN[k].wert


def _rk4(f, y, h, n, beob=None):
    for i in range(n):
        k1 = f(y)
        k2 = f(y + h / 2 * k1)
        k3 = f(y + h / 2 * k2)
        k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if beob is not None:
            beob(i, y)
    return y


# ================================================================== Tafel 3: Maxwell, Lorentz, Poynting, Wärme
POY_N = 400


def _poynting():
    """Feldverlust und Joule-Wärme (∂u/∂t + ∇·S = −J·E) auf einer Linie mit leitendem Stück, RK4."""
    x = np.arange(POY_N)
    sg = np.where((x > 250) & (x < 330), 0.02, 0.0)

    def f(y):
        E, H = y[:POY_N], y[POY_N:2 * POY_N - 1]
        dE = np.zeros(POY_N)
        dE[1:-1] = H[1:] - H[:-1] - sg[1:-1] * E[1:-1]
        return np.concatenate([dE, E[1:] - E[:-1], [np.sum(sg * E * E)]])
    E0 = np.exp(-((x - 100) / 15.0)**2)
    E0[0] = E0[-1] = 0
    y = np.concatenate([E0, -np.exp(-((x[:-1] + 0.5 - 100) / 15.0)**2), [0.0]])
    W0 = 0.5 * np.sum(y[:-1]**2)
    y = _rk4(f, y, 0.1, 4000)
    return W0 - 0.5 * np.sum(y[:-1]**2), y[-1]


_g("Z_poynting_joule", "Z", "Poynting ∂u/∂t + ∇·S = −J·E: Feldverlust = Joule-Wärme",
   lambda p: ident(lambda: _poynting()[0], lambda: _poynting()[1], 1e-6), T3 + " (II/III); 05", "Z", "transversal")


def _maxwell():
    rng = np.random.default_rng(SEED)
    n = 12

    def fw(f, a):
        return np.roll(f, -1, a) - f

    def bw(f, a):
        return f - np.roll(f, 1, a)
    B = np.zeros((3, n, n, n))
    D = np.zeros((3, n, n, n))
    rho = np.zeros((n, n, n))
    for _ in range(40):
        Ef, Hf, J = rng.normal(size=(3, 3, n, n, n))
        B -= 0.1 * np.array([fw(Ef[2], 1) - fw(Ef[1], 2), fw(Ef[0], 2) - fw(Ef[2], 0), fw(Ef[1], 0) - fw(Ef[0], 1)])
        D += 0.1 * (np.array([bw(Hf[2], 1) - bw(Hf[1], 2), bw(Hf[0], 2) - bw(Hf[2], 0), bw(Hf[1], 0) - bw(Hf[0], 1)]) - J)
        rho -= 0.1 * sum(bw(J[i], i) for i in range(3))
    divB = abs(sum(fw(B[i], i) for i in range(3))).max()
    gauss = abs(sum(bw(D[i], i) for i in range(3)) - rho).max()
    return wahr(divB < 1e-10 and gauss < 1e-10, f"{max(divB, gauss):.0e}")


_g("Z_maxwell_struktur", "Z", "Maxwell als Strukturbedingung: ∇·B = 0 und ∇·D = ρ bleiben erhalten (Gitter 12³, Zufallsfelder)",
   lambda p: _maxwell(), T3 + " (I); 05", "Z", "tangential")


def _kreuz(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _lorentz():
    v = (1.0, 0.3, 0.2)
    t = (0.2 * 0.025, 0.5 * 0.025, 1.0 * 0.025)
    f = 2 / (1 + sum(x * x for x in t))
    s = tuple(f * x for x in t)
    v0 = math.sqrt(sum(x * x for x in v))
    for _ in range(20000):
        c = _kreuz(v, t)
        vp = (v[0] + c[0], v[1] + c[1], v[2] + c[2])
        c = _kreuz(vp, s)
        v = (v[0] + c[0], v[1] + c[1], v[2] + c[2])
    return math.sqrt(sum(x * x for x in v)), v0


_g("Z_lorentz_ohne_arbeit", "Z", "Lorentzkraft f = ρE + J×B tangential: |v| im Magnetfeld erhalten (Boris-Schritt)",
   lambda p: ident(lambda: _lorentz()[0], lambda: _lorentz()[1], 1e-12), T3 + " (Lorentz); 09", "Z", "tangential")

# Tafel 3, Wärme: Durchgehen bei q = 1


def _durchgeht(I2):
    """Läuft T über 10 (Durchgehen) oder kommt es zur Ruhe (Rate unter 10⁻¹¹, fester Punkt)?"""
    def ueber10(t, T):
        return T[0] - 10
    ueber10.terminal = True

    def ruhe(t, T):
        return I2 * math.exp(T[0]) - T[0] - 1e-11
    ruhe.terminal = True
    s = solve_ivp(lambda t, T: I2 * np.exp(T) - T, (0, 1e7), [0.0], rtol=1e-10, atol=1e-12, events=(ueber10, ruhe))
    return len(s.t_events[0]) > 0


def _durchgehen():
    lo, hi = 0.30, 0.40
    for _ in range(25):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if _durchgeht(mid) else (mid, hi)
    return lo


_g("Z_durchgehen", "Z", "Thermisches Durchgehen bei Schleifenverstärkung q = 1: I²R₀α/G = 1/e (R = R₀e^{αT})",
   lambda p: sim(lambda: 1 / math.e, _durchgehen, 1e-6), T3 + " (Wärme); 09", "Z", "transversal")
Q_WAERME = 0.2


def _verlangsamung():
    Ts = brentq(lambda T: Q_WAERME * math.exp(T) - T, 0, 1)
    s = solve_ivp(lambda t, T: Q_WAERME * np.exp(T) - T, (0, 3), [Ts + 1e-6], rtol=1e-12, atol=1e-16)
    return -3 / math.log((s.y[0, -1] - Ts) / 1e-6), Ts


_g("Z_verstaerkung_verlangsamung", "Z", "Verstärkung ist Verlangsamung: τ = τ₀/(1 − q) = τ₀·Σqⁿ (Wärmeschleife)",
   lambda p: sim(lambda: K.schleifensumme(_verlangsamung()[1]), lambda: _verlangsamung()[0], 1e-4), T3 + "; 05", "Z",
   "transversal")


# ================================================================== Tafel 4: Schnitt, Rand, Schließung, Faserbündel
RING_N = 400


@functools.lru_cache(maxsize=None)
def _ring_strecke():
    """Kleinste Eigenwerte von −D∂² auf dem Ring (periodisch, dünn besetzt, Shift-Invert) und auf der Strecke mit
    festen Rändern (tridiagonal); D = L = 1, Gitter mit RING_N Punkten."""
    from scipy.linalg import eigvalsh_tridiagonal
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh
    h = 1.0 / RING_N
    ring = diags([-1.0, 2.0, -1.0, -1.0, -1.0], [-1, 0, 1, RING_N - 1, 1 - RING_N], shape=(RING_N, RING_N)) / h**2
    r = np.sort(eigsh(ring.tocsc(), k=3, sigma=-1.0, which="LM", return_eigenvectors=False))
    hs = 1.0 / (RING_N + 1)
    st = np.sort(eigvalsh_tridiagonal(np.full(RING_N, 2 / hs**2), np.full(RING_N - 1, -1 / hs**2),
                                      select="i", select_range=(0, 2)))
    return r, st


def _schliessung():
    r, s = _ring_strecke()
    return wahr(abs(r[0]) < 1e-6 and s[0] > 1.0, f"{r[0]:.0e} / {s[0]:.2f}")


_g("Z_schliessung_nullmode", "Z", "Schließung: Ring hat einen Nullmode (Gesamtmenge bleibt), Strecke mit Rand keinen",
   lambda p: _schliessung(), T4 + "; 05", "Z", "transversal")
_g("Z_spektralluecke_kreis", "Z", "Spektrallücke des geschlossenen Kreises 4π²D/L² (D = L = 1)",
   lambda p: sim(lambda: 4 * math.pi**2, lambda: _ring_strecke()[0][1], 1e-4), T4 + "; 05", "Z", "transversal")
_g("Z_fluchtrate_strecke", "Z", "Fluchtrate der aufgeschnittenen Strecke π²D/L²",
   lambda p: sim(lambda: math.pi**2, lambda: _ring_strecke()[1][0], 1e-4), T4 + "; 05", "Z", "transversal")


def _parametrisch():
    def g(t):
        return 1 + 3 * (3 * (t / 3000)**2 - 2 * (t / 3000)**3) if t < 3000 else 4.0
    s = solve_ivp(lambda t, z: [z[1], -g(t) * z[0]], (0, 3000), [1.0, 0.0], method="DOP853", rtol=1e-9, atol=1e-11)
    x, v = s.y[:, -1]
    return (0.5 * v**2 + 2 * x**2) / 2.0


_g("Z_parametrisch_wirkung", "Z", "Parametrisch: Wirkung E/ω bleibt bei langsamer Änderung von g (g × 4 über ≈ 700 Perioden)",
   lambda p: sim(lambda: 0.5, _parametrisch, 1e-3), T4 + " (Faserbündel); 05", "Z", "parametrisch")


# ================================================================== Tafel 1: Schleifen der Wirkfelder
BALKEN_N = 300


def _vorspannung():
    """Gelenkig gelagerter Balken: kleinster Eigenwert von ∂⁴ + P∂² (Differenzen, Bandmatrix)."""
    from scipy.linalg import eig_banded
    n, f = BALKEN_N, (BALKEN_N + 1)**2
    D2 = (np.eye(n, k=1) + np.eye(n, k=-1) - 2 * np.eye(n)) * f

    def w2(P):
        A = D2 @ D2 + P * D2
        band = np.array([np.concatenate([np.zeros(k), np.diag(A, k)]) for k in (2, 1, 0)])
        return eig_banded(band, eigvals_only=True, select="i", select_range=(0, 0))[0]
    return w2(0.5 * math.pi**2) / w2(0.0)


_g("Z_vorspannung_resonanz", "Z",
   "Gravitation: Vorspannung verschiebt Resonanz ω²(P)/ω²(0) = 1 − P/P_E (P = P_E/2; Knicken bei P_E = π²EI/L²)",
   lambda p: sim(lambda: 0.5, _vorspannung, 1e-3), T1 + " (Gravitation); 05", "Z", "parametrisch")

SELBST_A, SELBST_R, SELBST_B = 0.06, 0.01, 1.0     # aktives Element a, Verlust R, Sättigung b (ω = m = 1)


def _selbsterregung():
    """Rayleigh-Schwinger ẍ + x + (R − a)ẋ + bẋ³ = 0: Geschwindigkeitsamplitude des Grenzzyklus."""
    tt = np.linspace(500, 600, 5001)
    s = solve_ivp(lambda t, z: [z[1], -z[0] - (SELBST_R - SELBST_A) * z[1] - SELBST_B * z[1]**3], (0, 600),
                  [0.01, 0.0], rtol=1e-9, atol=1e-12, t_eval=tt)
    return np.max(np.abs(s.y[1]))


_g("Z_selbsterregung", "Z",
   "Strömung/Selbsterregung: Grenzzyklus genau bei q_self = a/(R + ¾bV²) = 1, V*² = 4(a − R)/(3b)",
   lambda p: sim(lambda: math.sqrt(4 * (SELBST_A - SELBST_R) / (3 * SELBST_B)), _selbsterregung, 1e-2),
   T1 + " (Strömung); 05", "Z", "transversal")
ABSTR_R, ABSTR_L = 0.08, 1.0                        # Abstrahlungswiderstand, Induktivität (ω = 1)


def _abstrahlung():
    Td = 2 * math.pi / math.sqrt(1 - (ABSTR_R / (2 * ABSTR_L))**2)       # Periode der gedämpften Schwingung
    s = solve_ivp(lambda t, z: [z[1], -z[0] - ABSTR_R / ABSTR_L * z[1]], (0, 11 * Td), [1.0, 0.0], rtol=1e-11,
                  atol=1e-13, dense_output=True)
    e = lambda t: 0.5 * (s.sol(t)[0]**2 + s.sol(t)[1]**2)                  # noqa: E731  Energie, gemittelt über Perioden
    t1, t2 = 5 * Td, 9 * Td
    m1 = quad(e, t1, t1 + Td)[0]
    m2 = quad(e, t2, t2 + Td)[0]
    return math.log(m1 / m2) / (t2 - t1) / 2


_g("Z_abstrahlung", "Z", "Akustik/Abstrahlung: reiner Verlust (q < 0), Beitrag −R_rad/2L zur Querrate λ_⊥",
   lambda p: sim(lambda: ABSTR_R / (2 * ABSTR_L), _abstrahlung, 1e-6), T1 + " (Akustik, Abstrahlung); 05", "Z",
   "transversal")
ERM_M, ERM_D0 = 2.0, 1.0                            # Exponent, Anfangsrate Ḋ₀ (D₀ = 0)


def _ermuedung():
    """Zeit bis D = 1 aus dem Schädigungsgesetz: t = ∫₀¹ dD/Ḋ(D) (numerisch)."""
    return quad(lambda D: (1 - D)**ERM_M / ERM_D0, 0, 1)[0]


_g("Z_ermuedung", "Z", "Ermüdung (Energie → Zustand → Struktur): Ḋ ∝ (1 − D)^{−m} bricht nach t_f = 1/((m + 1)Ḋ₀) (m = 2)",
   lambda p: sim(lambda: 1 / ((ERM_M + 1) * ERM_D0), _ermuedung, 1e-5), T1 + " (Ermüdung); 05", "Z", "parametrisch")
ALT_EA, ALT_T = 30.0, 2.0                           # Aktivierung und Temperatur in derselben Einheit


def _alterung():
    rate = lambda T: math.exp(-ALT_EA / T)                                   # noqa: E731
    h = 1e-5
    return (math.log(rate(ALT_T + h)) - math.log(rate(ALT_T - h))) / (2 * h)


_g("Z_alterung_arrhenius", "Z", "Chemie/Alterung: Empfindlichkeit ∂ln Ȧ/∂T = E_a/T² (Arrhenius), einseitige Drift",
   lambda p: ident(_alterung, lambda: ALT_EA / ALT_T**2, 1e-8), T1 + " (Chemie); 05", "Z", "parametrisch")
RAUSCH_LAM, RAUSCH_D, RAUSCH_N, RAUSCH_DT, RAUSCH_SCHRITTE = (1.0, 0.25), 1.0, 2000, 0.01, 12000


def _rauschen():
    """Banddicke² ⟨x²⟩ eines Ornstein-Uhlenbeck-Bands für zwei Querraten (Euler-Maruyama); Verhältnis."""
    var = []
    for j, lam in enumerate(RAUSCH_LAM):
        rng = np.random.default_rng(SEED + 20 + j)
        x = np.zeros(RAUSCH_N)
        acc, n = 0.0, 0
        for i in range(RAUSCH_SCHRITTE):
            x += -lam * x * RAUSCH_DT + math.sqrt(2 * RAUSCH_D * RAUSCH_DT) * rng.standard_normal(RAUSCH_N)
            if i > RAUSCH_SCHRITTE // 2 and i % 20 == 0:
                acc += np.mean(x * x)
                n += 1
        var.append(acc / n)
    return var[1] / var[0]


_g("Z_rauschen_banddicke", "Z", "Rauschen: Banddicke ∝ √(D_⊥/|λ_⊥|), wächst für λ_⊥ → 0 (⟨x²⟩-Verhältnis bei λ = 1 und ¼)",
   lambda p: sim(lambda: RAUSCH_LAM[0] / RAUSCH_LAM[1], _rauschen, 0.04), T1 + " (Rauschen); 05", "Z", "transversal")
KRIECH_E, KRIECH_ETA, KRIECH_T = 2.0, 5.0, 4.0


def _kriechen():
    """Spannung bei fester Dehnung (Feder und Dämpfer in Reihe): σ̇ = −(E/η)σ, explizit in kleinen Schritten."""
    n = 400000
    dt = KRIECH_T / n
    s = 1.0
    for _ in range(n // 1000):
        s *= (1 - KRIECH_E / KRIECH_ETA * dt)**1000
    return s


_g("Z_kriechen", "Z", "Kriechen: Spannung relaxiert σ = σ₀e^{−Et/η} (q < 0, stabilisierend, irreversibel)",
   lambda p: sim(lambda: math.exp(-KRIECH_E / KRIECH_ETA * KRIECH_T), _kriechen, 1e-5), T1 + " (Kriechen); 05", "Z",
   "transversal")


# ================================================================== Tafel 2: drei Rollen (GENERIC)
def _generic():
    zet, Cw = 0.3, 1.0
    s = solve_ivp(lambda t, z: [z[1], -z[0] - zet * z[1], zet * z[1]**2], (0, 60), [1.0, 0.0, 1.0], rtol=1e-12,
                  atol=1e-14, dense_output=True)
    zz = s.sol(np.linspace(0, 60, 3000))
    E = 0.5 * zz[1]**2 + 0.5 * zz[0]**2 + zz[2]
    S = Cw * np.log(zz[2])
    rng = np.random.default_rng(SEED + 3)
    ok = True
    for _ in range(100):
        x, pp, e = rng.normal(), rng.normal(), rng.uniform() + 0.1
        Lm = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
        vv = np.array([0, 1, -pp])
        Mm = zet * e / Cw * np.outer(vv, vv)
        gE, gS = np.array([x, pp, 1]), np.array([0, 0, Cw / e])
        ok &= np.allclose(Lm @ gS, 0) and np.allclose(Mm @ gE, 0) and np.allclose(Lm, -Lm.T)
        ok &= min(np.linalg.eigvalsh(Mm)) > -1e-12 and np.allclose(Lm @ gE + Mm @ gS, [pp, -x - zet * pp, zet * pp**2])
    return wahr(ok and np.ptp(E) < 1e-9 and np.all(np.diff(S) > -1e-12), f"ΔE {np.ptp(E):.0e}, ΔS {S[-1] - S[0]:+.3f}")


_g("Z_generic_rollen", "Z", "GENERIC: 𝓛 = −𝓛ᵀ, 𝓛∇S = 0, 𝓜 = 𝓜ᵀ ≥ 0, 𝓜∇E = 0; E fest, S steigt (Oszillator mit Wärme)",
   lambda p: _generic(), T2 + " (drei Rollen); 05; Öttinger 2005", "Z", "tangential")


# ================================================================== 06: Einrasten, Zungen, Kuramoto, Koexistenz, Tod
def _einrasten():
    dw = np.linspace(0.05, 1.95, 39)
    ph = np.array([0, math.pi / 3, math.pi / 2, 2 * math.pi / 3])
    DW, PH = np.meshgrid(dw, ph)
    k = 0.5
    zw = {}

    def f(psi):
        return DW - 2 * k * np.cos(PH) * np.sin(psi)

    def beob(i, y):
        if i == 20000:
            zw["m"] = y.copy()
    psi = _rk4(f, np.zeros_like(DW), 0.02, 40000, beob)
    rast = np.abs(psi - zw["m"]) / 400 < 1e-3
    grenze = np.vectorize(K.rastbreite)(k, PH)
    erw = np.abs(DW) <= grenze
    rand = np.abs(np.abs(DW) - grenze) < 0.03
    return wahr(np.all((rast == erw) | rand), f"{int(rast.sum())} von 156 gerastet")


_g("Z_einrasten", "Z", "Einrasten genau für |Δω| ≤ 2k|cos φ|; bei φ = π/2 nie (156 Paare Δω, φ)", lambda p: _einrasten(),
   "06; Adler 1946", "Z", "transversal")


def _rastrate():
    psi = math.asin(0.6)
    s = solve_ivp(lambda t, x: 0.6 - np.sin(x), (0, 5), [psi + 1e-6], rtol=1e-12, atol=1e-18)
    return math.log((s.y[0, -1] - psi) / 1e-6) / 5


_g("Z_rastrate", "Z", "Rate am Rastpunkt −√((2k cos φ)² − Δω²) (Δω = 0,6, k = 0,5, φ = 0; Abklingen einer Störung)",
   lambda p: sim(lambda: -math.sqrt(K.rastbreite(0.5, 0.0)**2 - 0.36), _rastrate, 1e-5), "06", "Z", "transversal")


def _zunge(k, m=3, n=2):
    det = np.linspace(0, 2 * (m + n) * k, 201)
    th1 = np.zeros_like(det)
    th2 = np.zeros_like(det)
    w1 = 2.0
    w2 = (m * w1 + det) / n
    h = 0.01
    pm = None
    for i in range(60000):
        ps = n * th2 - m * th1
        sp = k * np.sin(ps)
        th1 = th1 + h * (w1 + sp)
        th2 = th2 + h * (w2 - sp)
        if i == 30000:
            pm = n * th2 - m * th1
    frei = np.abs(n * th2 - m * th1 - pm) > 1.0
    return det[np.argmax(frei)] / ((m + n) * k)


for _k in (0.02, 0.04):
    _g(f"Z_einrastzunge_{int(_k * 100):02d}", "Z",
       f"Einrastzunge 3:2: Rand bei |nω₂ − mω₁| = (m + n)k (k = {_k:g}{', Breite linear' if _k > 0.03 else ''})",
       lambda p, k=_k: sim(lambda: 1.0, lambda: _zunge(k), 0.015), "06; Arnold 1965", "Z", "tangential")


def _kuramoto_selbstkonsistent(N):
    from scipy.optimize import minimize_scalar
    om = np.linspace(-0.5, 0.5, N)
    r = minimize_scalar(lambda u: u / np.mean(np.sqrt(np.clip(1 - (om / u)**2, 0, None))), bounds=(0.5, 2),
                        method="bounded")
    return (N - 1) / N * r.fun


KURA_N = 40


def _kuramoto_sim():
    rng = np.random.default_rng(SEED + 4)
    om = np.linspace(-0.5, 0.5, KURA_N)
    kv = np.linspace(0.55, 0.70, 61) / (KURA_N - 1)
    th = np.tile(rng.uniform(0.0, 1.0, KURA_N) * 2 * math.pi, (kv.size, 1))
    acc = np.zeros_like(th)
    kk = kv[:, None] * KURA_N
    for i in range(20000):
        c, s = np.cos(th), np.sin(th)
        d = om + kk * (s.mean(1, keepdims=True) * c - c.mean(1, keepdims=True) * s)
        th = th + 0.1 * d
        if i >= 15000:
            acc += d
    gerastet = np.ptp(acc / 5000, 1) < 1e-3
    return (KURA_N - 1) * kv[np.argmax(gerastet)]


_g("Z_kuramoto", "Z", "Kuramoto N = 40: volles Einrasten (N − 1)k/W, Simulation gegen Selbstkonsistenz",
   lambda p: sim(lambda: _kuramoto_selbstkonsistent(KURA_N), _kuramoto_sim, 5e-3), "06; Kuramoto 1975", "Z", "transversal")
_g("Z_torusdimension_grenze", "Z", "Grenze N_max ≈ 1 + 2W/(πk): (N − 1)k_c/W → 2/π (N = 4000)",
   lambda p: ident(lambda: _kuramoto_selbstkonsistent(4000), lambda: 2 / math.pi, 2e-3), "06; Kuramoto: K_c = 2/(πg(0))",
   "Z", "transversal")


def _koexistenz():
    def koex(a1, a2, b1=1.2, b2=0.9, t12=0.3, t21=0.8):
        s = solve_ivp(lambda t, I: [2 * I[0] * (a1 - b1 * I[0] - t12 * I[1]), 2 * I[1] * (a2 - b2 * I[1] - t21 * I[0])],
                      (0, 400), [0.1, 0.1], rtol=1e-9)
        return bool(min(s.y[:, -1]) > 1e-3)
    return wahr(all(koex(r, 1.0) == (0.3 / 0.9 < r < 1.2 / 0.8) for r in (0.2, 0.5, 1.0, 1.4, 2.0)), "5 Verhältnisse")


_g("Z_koexistenz", "Z", "Koexistenz zweier Zyklen genau für θ₁₂/β₂ < α₁/α₂ < β₁/θ₂₁", lambda p: _koexistenz(), "06", "Z",
   "transversal")


def _amplitudentod():
    def tod(Kd, dw, mu=1.0):
        w = np.array([0.0, dw])
        s = solve_ivp(lambda t, z: (mu + 1j * w - abs(z)**2) * z + Kd * (z[::-1] - z), (0, 150),
                      np.array([1.0 + 0j, 0.5j]), rtol=1e-7, atol=1e-14)
        return abs(s.y[:, -1]).max() < 1e-6
    return wahr(tod(2.0, 6.0) and not tod(2.0, 3.0) and not tod(0.8, 6.0), "3 Punkte")


_g("Z_amplitudentod", "Z", "Amplitudentod genau für K_d > μ und |Δω| > 2√(μ(2K_d − μ))", lambda p: _amplitudentod(),
   "06; Aronson, Ermentrout, Kopell 1990", "Z", "transversal")


def _preisach():
    rng = np.random.default_rng(SEED + 5)
    al, be = rng.uniform(-1, 1, 3000), rng.uniform(-1, 1, 3000)
    al, be = np.maximum(al, be), np.minimum(al, be)
    w = rng.uniform(0.0, 1.0, 3000)

    def y(folge):
        return K.preisach_ausgang(folge, al, be, w)
    ok = (y([-1, 0.8, -0.3, 0.6, 0.2]) == y([-1, 0.8, -0.1, 0.4, -0.3, 0.6, 0.2])
          and y([-1, 0.8, -0.3, 0.5]) == y([-1, 0.8, -0.3, 0.5, 0.1, 0.5]))
    return wahr(ok, "Folge der Extremwerte")


_g("Z_preisach", "Z", "Preisach: gleiche Folge der Extremwerte ⇒ gleicher Zustand; Rückkehrpunkt-Gedächtnis",
   lambda p: _preisach(), "06; Preisach 1935", "Z", "parametrisch")
KIPP_A, KIPP_D = 1.0, 1e-8


def _kipp_numerisch():
    """1 − q = f'(x*) am stabilen Punkt vor der Falte (Wurzel numerisch)."""
    x = min(np.real(np.roots([1, 0, -KIPP_A, -(2 * (KIPP_A / 3)**1.5 - KIPP_D)])))
    return 3 * x**2 - KIPP_A


_g("Z_kippstelle_wurzel", "Z", "Kippstelle: 1 − q = 2(3a)^{1/4}δ^{1/2} vor der Falte (a = 1, δ = 10⁻⁸)",
   lambda p: sim(lambda: 1 - K.kipp_schleife(KIPP_A, KIPP_D), _kipp_numerisch, 1e-3), "06", "Z", "parametrisch")
_g("Z_kippstelle_sprung", "Z", "Kippstelle: Sprung an der Falte √(3a)",
   lambda p: ident(lambda: max(np.real(np.roots([1, 0, -KIPP_A, -2 * (KIPP_A / 3)**1.5]))) + math.sqrt(KIPP_A / 3),
                   lambda: math.sqrt(3 * KIPP_A), 1e-6), "06", "Z", "parametrisch")


def _q3_krit(t3):
    def re_max(Q):
        co = np.poly([-1 / t for t in t3]) * np.prod(t3)
        co[-1] += Q
        return max(np.real(np.roots(co)))
    return brentq(re_max, 1e-6, 1e6, xtol=1e-14)


for _t3, _k in (((1.0, 1.0, 1.0), "Z_dreieck_gleich"), ((0.1, 1.0, 30.0), "Z_dreieck_getrennt")):
    _g(_k, "Z", f"Oszillatorischer Einsatz im Dreieck Energie–Struktur–Zustand, Zeiten {_t3}: |q₃| = (Σt)(Σ1/t) − 1",
       lambda p, t3=_t3: ident(lambda: _q3_krit(t3), lambda: sum(t3) * sum(1 / t for t in t3) - 1, 1e-8), "06", "Z",
       "transversal")


def _netz():
    rng = np.random.default_rng(SEED + 6)
    A = rng.uniform(0.0, 1.0, (6, 6))
    A *= 0.9 / max(abs(np.linalg.eigvals(A)))
    S, P = np.zeros((6, 6)), np.eye(6)
    for _ in range(600):
        S += P
        P = P @ A
    return np.linalg.inv(np.eye(6) - A).sum(), S.sum()


_g("Z_ruecklaeufe_netz", "Z", "Rückläufe im Netz: (I − A)⁻¹ = ΣAⁿ (Spektralradius 0,9)",
   lambda p: ident(lambda: _netz()[0], lambda: _netz()[1], 1e-9), "05; Leontief, Hawkins-Simon", "Z", "transversal")


# ================================================================== Tafel 5: Bindungsoperator
def _spiel_umkehr():
    v, gam = _E("SPIEL_V"), _E("SPIEL_GAMMA")
    ereignis = lambda t, z: z[1]                                              # noqa: E731
    s = solve_ivp(lambda t, z: [z[1], -1 / z[0]**2 - gam * z[2], z[1]], (0, 500), [1.0, v, 0.0], rtol=1e-12,
                  atol=1e-12, events=ereignis)
    return s.y_events[0][0][0]


def _spiel_erhaltung():
    """Umkehrpunkt aus H = v²/2 − 1/x + γy²/2 bei t_K → ∞ (y = x − 1): v₀²/2 − 1 + 1/x − γ(x − 1)²/2 = 0."""
    v, gam = _E("SPIEL_V"), _E("SPIEL_GAMMA")
    return brentq(lambda x: v * v / 2 - 1 + 1 / x - gam / 2 * (x - 1)**2, 2, 50)


_g("Z_gedaechtnis_bindet", "Z",
   "Gedächtniskern bindet ohne Energieverlust: Umkehrpunkt aus der Erhaltung H = v²/2 − 1/x + γy²/2 (v > √2)",
   lambda p: ident(_spiel_erhaltung, _spiel_umkehr, 1e-8), T5 + "; 05", "Z", "parametrisch",
   haengt_an=("SPIEL_V", "SPIEL_GAMMA"), wahl=W_TAFEL5)


def _maxwell_element(wt):
    re = quad(lambda u: math.exp(-u) * math.cos(wt * u), 0, 60, limit=400)[0]
    im = quad(lambda u: math.exp(-u) * math.sin(wt * u), 0, 60, limit=400)[0]
    return 1j * wt * (re - 1j * im)


for _wt in (0.5, 2.0):
    _g(f"Z_maxwell_speicher_{str(_wt).replace('.', '_')}", "Z", f"Maxwell-Element ωt_K = {_wt:g}: Speicheranteil (ωt_K)²/(1 + (ωt_K)²)",
       lambda p, wt=_wt: ident(lambda: _maxwell_element(wt).real, lambda: wt**2 / (1 + wt**2), 1e-9), T5 + "; 05", "Z",
       "tangential", wahl=W_TAFEL5)
    _g(f"Z_maxwell_verlust_{str(_wt).replace('.', '_')}", "Z", f"Maxwell-Element ωt_K = {_wt:g}: Verlustanteil ωt_K/(1 + (ωt_K)²)",
       lambda p, wt=_wt: ident(lambda: _maxwell_element(wt).imag, lambda: wt / (1 + wt**2), 1e-9), T5 + "; 05", "Z",
       "transversal", wahl=W_TAFEL5)
PQ_Q, PQ_N = 0.9, 64


def _operator():
    """𝒫_q = (1 − q)(I − qT)⁻¹ vertauscht mit der Translation T (kein Term im Generator) und hat vollen Rang (löscht
    keine Freiheitsgrade); für q → 1 wird er der Projektor 𝒫 = (1/N)ΣT^m auf den Nullmode."""
    T = np.roll(np.eye(PQ_N), 1, 0)
    P = (1 - PQ_Q) * np.linalg.inv(np.eye(PQ_N) - PQ_Q * T)
    P1 = (1 - 0.999999) * np.linalg.inv(np.eye(PQ_N) - 0.999999 * T)
    proj = np.ones((PQ_N, PQ_N)) / PQ_N
    ok = np.allclose(P @ T, T @ P) and np.linalg.matrix_rank(P) == PQ_N and np.allclose(P1, proj, atol=1e-4)
    return wahr(ok, f"Rang {np.linalg.matrix_rank(P)} von {PQ_N}")


_g("Z_bindungsoperator", "Z", "L_eff = 𝒫L: 𝒫_q vertauscht mit dem Generator, löscht keine Freiheitsgrade (voller Rang), "
   "q → 1 gibt den Projektor 𝒫", lambda p: _operator(), T5 + "; 05", "Z", "transversal", wahl=W_TAFEL5)
POISSON_R = 30.0                                    # in r*: tief im Bereich der flachen Kurve


def _poisson_numerisch(p):
    """ρ_X = ∇·(a − a_N)/(4πG) der Punktmasse, numerisch über (1/r²)d(r²(a − a_N))/dr."""
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    r = POISSON_R * K.r_stern(M, p.a0)
    zus = lambda x: x * x * (K.a_regel(K.G * M / x**2, p.a0) - K.G * M / x**2)          # noqa: E731
    h = 1e-4 * r
    return (zus(r + h) - zus(r - h)) / (2 * h) / r**2 / (4 * math.pi * K.G)


def _poisson_kern(p):
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    r = POISSON_R * K.r_stern(M, p.a0)
    return K.dichte_gelesen(K.G * M / r**2, r, 0.0, p.a0)


_g("Z_poisson_inversion", "Z",
   "Poisson-Inversion: die gelesene Dichte ρ_X = ∇·(a − a_N)/(4πG) ist Ausgabe der Regel (außen isothermer Halo)",
   lambda p: ident(lambda: _poisson_kern(p), lambda: _poisson_numerisch(p), 1e-6), T5 + "; 05", "Z", "parametrisch",
   haengt_an=("M_MILCHSTRASSE",))


def _additiv(p):
    m1, m2 = 1e9 * K.M_sun, 3e9 * K.M_sun
    fal = abs(K.ds(m1 + m2) / (K.ds(m1) + K.ds(m2)) - 1)
    r = 10 * K.r_stern(m1 + m2, p.a0)
    a1, a2 = K.G * m1 / r**2, K.G * m2 / r**2
    reg = K.a_regel(a1 + a2, p.a0) / (K.a_regel(a1, p.a0) + K.a_regel(a2, p.a0)) - 1
    return wahr(fal < 1e-14 and abs(reg) > 0.1, f"Faltung {fal:.0e}, Regel {reg:+.2f}")


_g("Z_faltung_additiv", "Z", "Faltung additiv (ΣΔs), die Regel nicht additiv je Teilchen (Operator statt Quelle)",
   _additiv, T5 + "; 01, 05", "Z", "parametrisch")


# ================================================================== Tafel 6: Trennbeobachtungen
def _linse_dyn(bb):
    """Test 1: Linsenmasse der Kernformel α(b) gegen die dynamische Masse im Zylinder b (M_dyn(r) = r²a(r)/G)."""
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    rs = K.r_stern(M)
    b = bb * rs
    a = lambda r: float(K.a_regel(K.G * M / r**2))                         # noqa: E731
    Md = lambda r: r * r * a(r) / K.G                                        # noqa: E731

    def dM(r):
        return (Md(r * 1.0001) - Md(r / 1.0001)) / (r * 1.0001 - r / 1.0001)
    Mz = Md(b) + quad(lambda u: dM(b * math.exp(u)) * (1 - math.sqrt(1 - math.exp(-2 * u))) * b * math.exp(u), 0, 40,
                      limit=400)[0]
    return 4 * K.G * Mz / (K.c**2 * b)


for _bb in (0.3, 1.0, 5.0):
    _g(f"Z_linsentest_kugel_{str(_bb).replace('.', '_')}", "Z",
       f"Test 1: Kugel, Linsenmasse = dynamische Masse, η = 1 (b = {_bb:g} r*)",
       lambda p, bb=_bb: ident(lambda: K.ablenkung_linse(EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun,
                                                          bb * K.r_stern(EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun)),
                               lambda: _linse_dyn(bb), 1e-5), T6 + " (Test 1); 05, 08", "Z", "Wirbel",
       haengt_an=("M_MILCHSTRASSE",))


WEIT_B = 200.0                                      # in r*


def _ablenkung_weit(p):
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    return K.ablenkung_linse(M, WEIT_B * K.r_stern(M, p.a0), p.a0)


def _ablenkung_weit_hand(p):
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    b = WEIT_B * math.sqrt(K.G * M / p.a0)
    return 2 * math.pi * math.sqrt(K.G * M * p.a0) / K.c**2 + 2 * K.G * M / (K.c**2 * b)


_g("Z_ablenkung_weit_aussen", "Z", "Ablenkung weit außen α → 2π(v_f/c)² + 6Δs/b (b = 200 r*)",
   lambda p: ident(lambda: _ablenkung_weit(p), lambda: _ablenkung_weit_hand(p), 1e-5), T6 + " (Test 1); 08", "Z",
   "Wirbel", haengt_an=("M_MILCHSTRASSE",))


def _wachstum(p):
    """Test 2: die Verstärkung einer Störung hängt an Skala (a_N) und Epoche (a₀(z))."""
    x = [0.1 * p.a0, 10 * p.a0]
    skala = K.nu(x[0], p.a0) / K.nu(x[1], p.a0)
    epoche = K.nu(x[0], p.a0 * K.a0_frueh_relativ(EINGABEN["Z_REKOMBINATION"].wert)) / K.nu(x[0], p.a0)
    return wahr(abs(skala - 1) > 0.1 and abs(epoche - 1) > 0.01, f"Skala {skala:.2f}, Epoche {epoche:.3f}")


_g("Z_trenntest_wachstum", "Z", "Test 2: Wachstum mit Skalen- und Epochenabhängigkeit (Quelle: nur über die Dichte)",
   _wachstum, T6 + " (Test 2); 05", "Z", "parametrisch", haengt_an=("Z_REKOMBINATION",))


def _umfeld(p):
    """Test 3: die gelesene Dichte hängt am äußeren Feld (Umfeld), nicht nur am eigenen Profil."""
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    ge = _E("GE_003") * p.a0
    r = _E("R_50KPC") * K.kpc

    def rho(aussen):
        def zus(x):
            aN = K.G * M / x**2
            n = K.nu_aussenfeld(aN, ge, "senkrecht", p.a0) if aussen else float(K.nu(aN, p.a0))
            return x * x * (n - 1) * aN
        h = 1e-4 * r
        return (zus(r + h) - zus(r - h)) / (2 * h) / r**2
    v = rho(True) / rho(False)
    return wahr(v < 0.99, f"ρ_X mit/ohne Außenfeld {v:.2f}")


_g("Z_trenntest_umfeld", "Z", "Test 3: Form der gelesenen Verteilung an das Umfeld gebunden (äußeres Feld 0,03·a₀, 50 kpc)",
   _umfeld, T6 + " (Test 3); 05", "Z", "parametrisch", haengt_an=("M_MILCHSTRASSE", "GE_003", "R_50KPC"))


def _gedaechtnis(p):
    """Test 4: gleicher Zustand heute (θ = 3H₀), verschiedene Geschichte → verschiedenes Gewicht w (Preisach)."""
    H = K.H0_si(p.H0)
    w1 = K.gewicht_w_preisach([6 * H, 3 * H], p.a0)
    w2 = K.gewicht_w_preisach([6 * H, 0.5 * H, 3 * H], p.a0)
    return wahr(abs(w1 - w2) > 1e-3, f"w {w1:.4f} gegen {w2:.4f}")


_g("Z_trenntest_gedaechtnis", "Z", "Test 4: Pfadabhängigkeit über w (Gedächtnis statt Markow)", _gedaechtnis,
   T6 + " (Test 4); 05", "Z", "parametrisch")


def _aequivalenz(p):
    """Test 5: Erde und Mond haben verschiedene Glättungslängen (ℓ³ = GM/𝒯), fallen am selben Ort aber gleich: die
    Regel liest nur das Newton-Feld am Ort, keine Eigenschaft des Körpers."""
    T = gezeit(p)
    lE, lM = K.glaettungslaenge(K.M_earth, T), K.glaettungslaenge(MESSWERTE["M_MOND"].wert, T)
    a = [K.a_regel(K.GM_sun / K.AU**2, p.a0) for _ in (K.M_earth, MESSWERTE["M_MOND"].wert)]
    return wahr(abs(lE / lM - 1) > 0.5 and a[0] == a[1], f"ℓ Erde/Mond {lE / lM:.2f}, a gleich")


_g("Z_trenntest_aequivalenz", "Z",
   "Test 5: Äquivalenzprinzip: alle Materie folgt der Regel gleich, trotz verschiedener Glättungslänge (Erde, Mond)",
   _aequivalenz, T6 + " (Test 5); 05, 07", "Z", "parametrisch", haengt_an=SONNENORT + ("M_MOND",),
   satz="Sonnenkreis_gemessen")


def _kanaele(p):
    """Test 6: in der Ebene einer Scheibe sieht Licht (Feldform, Wirbelanteil fehlt) innen weniger als die Bahn."""
    M = EINGABEN["M_MILCHSTRASSE"].wert * K.M_sun
    Rd = _E("RD_MILCHSTRASSE") * K.kpc
    S0 = M / (2 * math.pi * Rd**2)

    def v(R):
        gR = float(K.g_expscheibe(M, Rd, R))
        S = S0 * math.exp(-R / Rd)
        return K.a_scheibe_ebene(gR, S, 1.0, p.a0) / K.a_scheibe_ebene(gR, S, 0.0, p.a0)
    innen, aussen = v(Rd / 4), v(5 * Rd)
    return wahr(innen < 1 and aussen > 0.99, f"Licht/Bahn innen {innen:.2f}, bei 5R_d {aussen:.3f}")


_g("Z_trenntest_kanaele", "Z", "Test 6: kanalspezifische Reste: in Scheiben innen η = M_dyn/M_lens > 1, außen 1",
   _kanaele, T6 + " (Test 6); 05, 07", "Z", "Wirbel", haengt_an=("M_MILCHSTRASSE", "RD_MILCHSTRASSE"))


def _hodge():
    """Triade als Hodge-Zerlegung auf dem Torus (periodisches Gitter): Feld = harmonisch + Gradient + Rotation,
    eindeutig und orthogonal; harmonisch ist auf dem geschlossenen Raum das gleichförmige Feld."""
    n = 63                                          # ungerade: kein Nyquist-Mode ohne Vorzeichen
    rng = np.random.default_rng(SEED + 7)
    F = rng.normal(size=(2, n, n))
    k = np.fft.fftfreq(n) * 2 * math.pi
    kx, ky = np.meshgrid(k, k, indexing="ij")
    k2 = kx**2 + ky**2
    Fh = np.fft.fft2(F)
    harm = np.zeros_like(Fh)
    harm[:, 0, 0] = Fh[:, 0, 0]
    k2s = np.where(k2 == 0, 1.0, k2)
    proj = (kx * Fh[0] + ky * Fh[1]) / k2s
    grad = np.array([kx * proj, ky * proj])
    grad[:, 0, 0] = 0
    rot = Fh - harm - grad
    teile = [np.real(np.fft.ifft2(x)) for x in (harm, grad, rot)]
    summe = np.abs(sum(teile) - F).max()
    ortho = max(abs(np.sum(teile[i] * teile[j])) for i, j in ((0, 1), (0, 2), (1, 2))) / np.sum(F * F)
    div_rot = np.abs(kx * rot[0] + ky * rot[1]).max()
    rot_grad = np.abs(kx * grad[1] - ky * grad[0]).max()
    ok = summe < 1e-12 and ortho < 1e-12 and div_rot < 1e-9 and rot_grad < 1e-9 and max(np.ptp(x) for x in teile[0]) < 1e-12
    return wahr(ok, f"Rest {summe:.0e}, Überlapp {ortho:.0e}, div {div_rot:.0e}, rot {rot_grad:.0e}, "
                f"harmonisch {max(np.ptp(x) for x in teile[0]):.0e}")


AUSSAGEN.neu("Z_triade_form", TRIADE_GLIEDER["Hintergrund"][1], "Beziehung")
_g("Z_triade_hodge", "Z", "Triade als Form: jedes Feld auf dem geschlossenen Raum = Hintergrund (harmonisch, gleichförmig) "
   "+ Zentrum (Gradient) + Wirbel (Rotation), eindeutig und orthogonal", lambda p: _hodge(), "Baustein 01", "L",
   "Hintergrund", bezug="Z_triade_form")
