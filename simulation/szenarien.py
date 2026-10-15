# -*- coding: utf-8 -*-
"""Die Szenarien der Simulation. Jedes baut Körper auf, lässt sie laufen und prüft das Ergebnis gegen das Modell.

Rückgabe je Szenario: dict mit
  titel, text          – Überschrift und ein Satz, was gerechnet wird
  welt, T, X           – die Welt, Zeiten [s], Orte [m] (Proben × Körper × 3) für die Darstellung
  pruefungen           – Liste (Aussage, Modell, Simulation, Einheit, Toleranz relativ)
"""
import math
import numpy as np

import kern as K
from welt import Koerper, Welt


def _winkel(T, X, i, zentrum=None):
    """Mittlere Winkelgeschwindigkeit von Körper i um das Zentrum (oder den Ursprung)."""
    rel = X[:, i, :2] - (X[:, zentrum, :2] if zentrum is not None else 0.0)
    phi = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
    return (phi[-1] - phi[0]) / (T[-1] - T[0])


def atom():
    """Wasserstoff als klassische Kreisbahn: Proton und Elektron, Coulomb (Windung)."""
    r = K.bohr_radius()
    mu = K.m_e * K.m_p / (K.m_e + K.m_p)
    v_rel = math.sqrt(K.coulomb_konstante() * K.e**2 / (mu * r))
    M = K.m_e + K.m_p
    p = Koerper("Proton", K.m_p, [-r * K.m_e / M, 0, 0], [0, -v_rel * K.m_e / M, 0], ladung=K.e)
    el = Koerper("Elektron", K.m_e, [r * K.m_p / M, 0, 0], [0, v_rel * K.m_p / M, 0], ladung=-K.e)
    w = Welt([p, el])
    umlauf = 2 * math.pi * r / v_rel
    T, X, V, _ = w.lauf(3 * umlauf, 3000, proben=600)
    f = _winkel(T, X, 1, 0) / (2 * math.pi)
    v = np.mean(np.linalg.norm(V[:, 1] - V[:, 0], axis=1))
    f_modell = K.alpha * K.c / (2 * math.pi * K.bohr_radius())
    return dict(titel="Atom: Wasserstoff", welt=w, T=T, X=X,
                text="Elektron und Proton auf der Kreisbahn um den gemeinsamen Schwerpunkt, drei Umläufe, "
                     "nur Coulomb: die Neigung der Faltung ist dort 2,3·10³⁹-mal schwächer als die Windung.",
                pruefungen=[("Umlauffrequenz der Hülle", f_modell, f, "Hz", 1e-3),
                            ("Tempo der Hülle v/c = α", K.alpha, v / K.c, "", 1e-3),
                            ("Coulomb gegen Neigung", 2.27e39, K.coulomb_konstante() * K.e**2 / (K.G * K.m_e * K.m_p), "", 1e-3)])


PLANETEN = [("Merkur", 0.387098), ("Venus", 0.723332), ("Erde", 1.0), ("Mars", 1.523679), ("Jupiter", 5.2026)]


def sonnensystem():
    """Planeten auf Kreisbahnen im Trichter der Sonne; die Eigenzeit läuft mit."""
    kk = [Koerper(n, 1.0, [a * K.AU, 0, 0], [0, K.bahntempo(K.M_sun, a * K.AU), 0]) for n, a in PLANETEN]
    w = Welt(kk, hintergrund=K.M_sun)
    dauer = 12 * K.yr
    T, X, V, tau = w.lauf(dauer, 60000, proben=2000)
    iE = 2
    om = _winkel(T, X, iE)
    jahr = 2 * math.pi / om
    v = np.mean(np.linalg.norm(V[:, iE], axis=1))
    uhr = (dauer - tau[iE]) / (dauer / K.yr)
    pr = [("Erde: Umlaufzeit", 2 * math.pi * math.sqrt(K.AU**3 / K.GM_sun) / 86400, jahr / 86400, "Tage", 1e-5),
          ("Erde: Tempo v = c·√(3Δs☉/r)", K.bahntempo(K.M_sun, K.AU) / 1e3, v / 1e3, "km/s", 1e-5),
          ("Erde: Uhr im Sonnentrichter je Jahr", 3 * K.faltung(K.M_sun) / K.AU * K.yr, uhr, "s", 1e-3)]
    for i, (n, a) in enumerate(PLANETEN):
        Tsim = 2 * math.pi / _winkel(T, X, i) / K.yr
        pr.append((f"{n}: Umlaufzeit aus T² = 4π²r³/(3c²Δs☉)",
                   2 * math.pi * math.sqrt((a * K.AU)**3 / (3 * K.c**2 * K.faltung(K.M_sun))) / K.yr, Tsim, "Jahre", 1e-4))
    return dict(titel="Sonnensystem", welt=w, T=T, X=X,
                text="Fünf Planeten auf Kreisbahnen im Trichter der Sonne, zwölf Jahre; die Sonne ist ein Einzeltrichter "
                     "(Newton-Bereich), die Eigenzeit der Planeten läuft mit dem Takt √(1 − 6Δs☉/r) mit.",
                pruefungen=pr)


def merkur():
    """Merkur auf der Ellipse, mit und ohne Takt- und Faltungsglied; gemessen wird die Drehung des Perihels."""
    a, ecc = 0.387098 * K.AU, 0.205630
    rp = a * (1 - ecc)
    vp = math.sqrt(K.GM_sun * (1 + ecc) / rp)
    umlauf = 2 * math.pi * math.sqrt(a**3 / K.GM_sun)
    n_um = 10
    erg = {}
    for glied in (False, True):
        w = Welt([Koerper("Merkur", 1.0, [rp, 0, 0], [0, vp, 0])], takt_glied=glied, hintergrund=K.M_sun)
        T, X, V, _ = w.lauf(n_um * umlauf, n_um * 20000, proben=2000)
        x, v = X[-1, 0], V[-1, 0]
        # Runge-Lenz-Vektor zeigt auf das Perihel
        L = np.cross(x, v)
        A = np.cross(v, L) - K.GM_sun * x / np.linalg.norm(x)
        erg[glied] = (math.atan2(A[1], A[0]), w, T, X)
    je_jh = 100 * K.yr / umlauf / n_um
    bogen = math.degrees(1) * 3600
    mit, ohne = erg[True][0] * je_jh * bogen, erg[False][0] * je_jh * bogen
    modell = math.degrees(K.perihel_je_umlauf(K.M_sun, a, ecc)) * 3600 * 100 * K.yr / umlauf
    _, w, T, X = erg[True]
    return dict(titel="Merkur: Drehung des Perihels", welt=w, T=T, X=X,
                text="Merkur auf seiner Ellipse, zehn Umläufe; Takt und Faltungsstrecke als Glied erster Ordnung. "
                     "Gegenprobe ohne das Glied: das Perihel steht.",
                pruefungen=[("Perihel mit Takt und Faltung je Jahrhundert", modell, mit, "″", 2e-3),
                            ("Gegenprobe ohne Glied, Betrag unter", 0.01, abs(ohne), "″", None)])


GALAXIE_RADIEN = [5.0, 8.88, 15.0, 50.0]      # kpc (8,88 = r*)
M_GALAXIE = 6e10                              # M☉, Milchstraße als Punktmasse


def galaxie():
    """Teststerne um eine Galaxie als Punktmasse, einmal mit der Rückkehrregel, einmal nach Newton."""
    M = M_GALAXIE * K.M_sun
    def start(regel):
        kk = []
        for r in GALAXIE_RADIEN:
            aN = K.G * M / (r * K.kpc)**2
            a = K.verstaerkung(aN) * aN if regel else aN
            kk.append(Koerper(f"Stern {r:g} kpc".replace(".", ","), 1.0, [r * K.kpc, 0, 0], [0, math.sqrt(r * K.kpc * a), 0]))
        return kk
    dauer = 2 * K.Gyr
    w = Welt(start(True), regel=True, hintergrund=M)
    w_n = Welt(start(False), regel=False, hintergrund=M)
    T, X, V, _ = w.lauf(dauer, 20000, proben=2000)
    Tn, Xn, Vn, _ = w_n.lauf(dauer, 20000, proben=200)
    pr = []
    modell = {5.0: 249.25, 8.88: 214.41, 15.0: 196.21, 50.0: 178.07}
    for i, r in enumerate(GALAXIE_RADIEN):
        v = np.mean(np.linalg.norm(V[:, i], axis=1)) / 1e3
        rm = np.mean(np.linalg.norm(X[:, i], axis=1)) / K.kpc
        pr.append((f"Rotation bei {r:g} kpc (Regel)".replace(".", ","), modell[r], v, "km/s", 3e-3))
        pr.append((f"Bahnradius bleibt {r:g} kpc".replace(".", ","), r, rm, "kpc", 1e-3))
    v50n = np.mean(np.linalg.norm(Vn[:, 3], axis=1)) / 1e3
    pr.append(("Newton bei 50 kpc (ohne Regel, fällt wie Kepler)", math.sqrt(K.G * M / (50 * K.kpc)) / 1e3, v50n, "km/s", 3e-3))
    pr.append(("tiefer Wert v⁴ = GMa₀", 170.45, K.v_flach(M) / 1e3, "km/s", 1e-3))
    return dict(titel="Galaxie: Rotationskurve", welt=w, T=T, X=X,
                text="Vier Teststerne auf Kreisbahnen um eine Milchstraße als Punktmasse (6·10¹⁰ M☉), zwei Milliarden Jahre, "
                     "mit der Rückkehrregel a = a_N/(1 − q); Gegenprobe nach Newton.",
                pruefungen=pr)


def labor_kette():
    """Leiterkette aus Widerständen R (längs) und Leitwerten G (quer) mit Rückführung: V(0) = V_ein + V(n_R).
    Das ist die Erneuerung der Rückkehr u(0) = g_N + u(ϖ_R) als Schaltkreis; gelöst wird das Knotennetz."""
    R, Gq, n_R, N = 1e3, 1e-5, 20, 4000          # 1 kΩ, R·G = 0,01, Abgriff 20, Kette lang gegen n_R
    A = np.zeros((N + 1, N + 1))
    b = np.zeros(N + 1)
    A[0, 0], A[0, n_R], b[0] = 1.0, -1.0, 1.0    # Rückführung, V_ein = 1 V
    for k in range(1, N):
        A[k, k - 1] = A[k, k + 1] = 1 / R
        A[k, k] = -2 / R - Gq
    A[N, N - 1], A[N, N] = 1 / R, -1 / R - Gq
    V = np.linalg.solve(A, b)
    T = np.arange(N + 1, dtype=float)             # Knotennummer statt Zeit
    X = np.zeros((N + 1, 1, 3)); X[:, 0, 0] = V
    w = Welt([Koerper("Spannung", 1.0, [0, 0, 0], [0, 0, 0])])
    nu = V[0]
    return dict(titel="Labor: RC-Kette mit Rückführung", welt=w, T=T[:400], X=X[:400],
                text="Leiterkette (R = 1 kΩ, R·G = 0,01) mit Rückführung des Abgriffs 20 auf den Eingang, als Knotennetz gelöst; "
                     "die Spannung am Eingang ist die Verstärkung ν. Die Daten tragen statt der Zeit die Knotennummer.",
                pruefungen=[("Verstärkung ν = 1/(1 − e^{−s}), s = n·arcosh(1 + RG/2)", K.rc_kette_verstaerkung(n_R, R * Gq), nu, "", 1e-6),
                            ("Abfall je Glied e^{−arcosh(1 + RG/2)}", math.exp(-math.acosh(1 + R * Gq / 2)), V[2] / V[1], "", 1e-6)])


def _mc_fluchtmass(s, n=400_000, saat=1):
    """Zweiter Weg zu ν_G: Monte-Carlo im Phasenraum der homogenen Kugel; je Probe wird geprüft, ob die Bahn in R
    bleibt (Apozentrum ≤ R ⇔ E ≤ Φ(R) + L²/2R²), Regel gegen Newton."""
    from scipy.integrate import cumulative_trapezoid
    rng = np.random.default_rng(saat)
    a0_ = 1 / s**2
    r = np.linspace(1e-6, 1.0, 20001)
    def mass(regel):
        nu = 1 / (-np.expm1(-np.sqrt(r / a0_))) if regel else 1.0
        phi = cumulative_trapezoid(nu * r, r, initial=0.0)
        x = rng.uniform(-1, 1, (n, 3)); rr = np.linalg.norm(x, axis=1); m = rr < 1
        x, rr = x[m], rr[m]
        vmax = math.sqrt(2 * phi[-1])
        v = rng.uniform(-vmax, vmax, x.shape)
        vr = np.sum(x * v, 1) / rr; v2 = np.sum(v * v, 1)
        E = 0.5 * v2 + np.interp(rr, r, phi)
        return np.count_nonzero(E <= phi[-1] + rr**2 * (v2 - vr**2) / 2) * (2 * vmax)**3 / n
    return mass(True) / mass(False)


def kosmos():
    """Atmen der Zelle mit beiden Kandidaten der Dunklen Energie, verbunden wie bei der Dunklen Materie:
    Grundglied ist die Rückstellung (Vakuumteil ρ_s), die Rückkehrsumme das Fluchtmaß des Ganzen nach Kac, ν_G = μ_R/μ_N.
    Lesart L: am Raum, die Hülle der Zelle wächst bei fester Faltung (s ∝ 1/a). Lesart I: zeitfest in der inneren Richtung
    (s = r*/Λ). Integriert wird das Atmen a(t) aus den Dichten heute (ω_b, ω_c, T_CMB), H₀ folgt."""
    from scipy.integrate import solve_ivp, quad
    rho_r, rho_b, rho_m = K.dichte_heute()
    rho_s = K.dichte_vakuumteil()
    s0, sI = K.schwelle_zelle(), K.form_zelle()
    nuG = K.fluchtmass_ganzes_glatt
    les = {"L": lambda a: nuG(s0 / a), "I": lambda a: nuG(sI), "T": lambda a: float(K.B1_ZELLE)}
    H = {k: (lambda a, f=f: math.sqrt(8 * math.pi * K.G / 3 * (rho_r / a**4 + rho_m / a**3 + f(a) * rho_s)))
         for k, f in les.items()}
    lnA = np.linspace(math.log(1e-4), math.log(3.0), 1200)
    T, A, aus = {}, {}, {}
    for k in les:
        t_start = quad(lambda x: 1 / H[k](math.exp(x)), math.log(1e-12), lnA[0])[0]
        s = solve_ivp(lambda x, y: [1 / H[k](math.exp(x))], (lnA[0], lnA[-1]), [t_start], t_eval=lnA, rtol=1e-10, atol=1)
        T[k], A[k] = s.y[0], np.exp(lnA)
        t0 = float(np.interp(0.0, lnA, s.y[0]))
        t0_quad = quad(lambda x: 1 / H[k](math.exp(x)), math.log(1e-12), 0.0, limit=400, epsabs=0, epsrel=1e-9)[0]
        H0k = H[k](1.0)
        w = lambda a, f=les[k]: -1 - (math.log(f(a * 1.001)) - math.log(f(a / 1.001))) / (3 * 2 * math.log(1.001))
        aus[k] = dict(nu=les[k](1.0), H0=H0k * K.Mpc / 1e3, t0=t0 / K.Gyr, t0_quad=t0_quad / K.Gyr,
                      Om_DE=les[k](1.0) * rho_s / (3 * H0k**2 / (8 * math.pi * K.G)), w0=w(1.0), w1=w(0.5))
    # Gesehen statt jetzt: der Abstand zur letzten Streufläche nimmt jede Schale in ihrer Epoche; gleichwertiger fester Anteil
    rho_L_planck = 3 * K.H0**2 / (8 * math.pi * K.G) - rho_m - rho_r
    D_planck = K.abstand_mitbewegt(lambda a: rho_L_planck)
    gesehen = {k: (K.abstand_mitbewegt(lambda a, f=f: f(a) * rho_s) / D_planck, K.anteil_gesehen(lambda a, f=f: f(a) * rho_s))
               for k, f in les.items()}
    # Darstellung: zwei „Körper“, deren x-Koordinate a(t) auf gemeinsamer Zeitachse ist
    Tg = np.linspace(max(T["L"][0], T["I"][0]), min(T["L"][-1], T["I"][-1]), 1500)
    Tg = np.linspace(max(T[k][0] for k in T), min(T[k][-1] for k in T), 1500)
    X = np.zeros((len(Tg), 3, 3))
    for i, k in enumerate(("L", "I", "T")):
        X[:, i, 0] = np.interp(Tg, T[k], A[k])
    w_ = Welt([Koerper("a(t) Lesart " + k, 1.0, [0, 0, 0], [0, 0, 0]) for k in ("L", "I", "T")])
    soll = (1 - K.Omega_b - K.Omega_c - rho_r / (3 * K.H0**2 / (8 * math.pi * K.G))) / K.vakuumteil()
    L_, I_ = aus["L"], aus["I"]
    z = lambda x, n=3: f"{x:.{n}f}".replace(".", ",").replace("-", "−")
    text = (f"Dunkle Energie als Rückstellung mal Rückkehrzahl des Ganzen, ρ_DE = ν_G·ρ_s. Lesart L (r*/R = {z(s0)} heute): "
            f"ν_G = {z(L_['nu'])}, verlangt sind nach Planck {z(soll)}; H₀ = {z(L_['H0'], 2)} km/s/Mpc, Weltalter {z(L_['t0'], 2)} Mrd. J., "
            f"Ω_DE = {z(L_['Om_DE'])}, w = {z(L_['w0'], 2)} heute und {z(L_['w1'], 2)} bei z = 1 (die Hülle wächst bei fester Faltung, "
            f"ν_G steigt). Lesart I (r*/Λ = {z(sI)}, zeitfest): ν_G = {z(I_['nu'])}, H₀ = {z(I_['H0'], 2)} km/s/Mpc, Weltalter "
            f"{z(I_['t0'], 2)} Mrd. J., Ω_DE = {z(I_['Om_DE'])}, w = −1. Gesehen, wie die Messung liest (Abstand zur letzten Streufläche, "
            f"jede Schale in ihrer früheren Epoche): Lesart L wirkt wie ein fester Anteil {z(gesehen['L'][1])} statt {z(soll)}, "
            f"der Abstand ist {z(100 * (gesehen['L'][0] - 1), 1)} % länger als nach Planck; Lesart I {z(gesehen['I'][1])}, "
            f"{z(100 * (gesehen['I'][0] - 1), 1)} %. Triade des Ganzen (b₁ = 3 Windungsklassen, fester Anteil 3): H₀ = {z(aus['T']['H0'], 2)} km/s/Mpc, "
            f"Weltalter {z(aus['T']['t0'], 2)} Mrd. J., Ω_DE = {z(aus['T']['Om_DE'])}, Abstand zur letzten Streufläche {z(100 * (gesehen['T'][0] - 1), 1)} % "
            f"gegen Planck. Die Daten tragen a(t) der Lesarten L, I, T als x-Koordinate.")
    return dict(titel="Kosmos: Atmen mit Rückstellung und Fluchtmaß", welt=w_, T=Tg, X=X, text=text,
                pruefungen=[("ν_G Lesart L: Quadratur gegen Monte-Carlo", _mc_fluchtmass(s0), L_["nu"], "", 0.02),
                            ("ν_G Lesart I: Quadratur gegen Monte-Carlo", _mc_fluchtmass(sI, saat=2), I_["nu"], "", 0.02),
                            ("Weltalter L: Integral gegen Lauf", L_["t0_quad"], L_["t0"], "Mrd. J.", 1e-4),
                            ("Weltalter I: Integral gegen Lauf", I_["t0_quad"], I_["t0"], "Mrd. J.", 1e-4),
                            ("Gesehen: fester Anteil von Lesart I gibt sich selbst zurück", I_["nu"], gesehen["I"][1], "", 1e-6),
                            ("H₀ Triade aus dem Kern gegen den Lauf", K.atmen_heute()[0], aus["T"]["H0"], "km/s/Mpc", 1e-9),
                            ("H₀ L von Hand 100·√(ω_m + ω_r + ν_G·ω_s)",
                             100 * math.sqrt((rho_m + rho_r + L_["nu"] * rho_s) / K.rho_krit100), L_["H0"], "km/s/Mpc", 1e-6)])


def daten():
    """Prüfung an unseren Daten: SPARC (Regel, Form, Dimension, Unabhängigkeit) und AME2020 (Kerne)."""
    import daten as D
    yd = D.sparc_eichung(K.kpc)
    a0s, band = D.sparc_schwelle(K.kpc, yd=yd)
    formen, n = D.sparc_formen(K.kpc, yd)
    ds, eds = D.sparc_dimension(K.kpc, yd=yd)
    gs, egs = D.sparc_gasanteil(K.kpc, a0s, yd=yd)
    eig = D.sparc_eigenschaften(K.kpc, a0s, yd=yd)
    import kerne as KE
    ame = KE.bekannte()
    gb, go, _, _ = D._einfache_punkte(K.kpc, yd)
    X = np.zeros((len(gb), 1, 3)); X[:, 0, 0], X[:, 0, 1] = np.log10(gb), np.log10(go)
    s = (band[1] - band[0]) / 2
    z = lambda x, k=2: f"{x:.{k}f}".replace(".", ",").replace("-", "−")
    text = (f"SPARC, {len(D.sparc_auswahl(yd))} Galaxien, {n} Punkte, Abstand, Neigung und Masse-zu-Licht je Galaxie mitgeschätzt, "
            f"Masse-zu-Licht der Scheiben so, dass Stern und Gas derselben Regel folgen (Υ_d = {z(yd)}): "
            f"a₀ = {z(a0s / 1e-10)}·10⁻¹⁰ m/s² (16–84 %: {z(band[0] / 1e-10)}–{z(band[1] / 1e-10)}); χ² Bose {z(formen['Bose'], 1)}, "
            f"einfach {z(formen['einfach'], 1)}, Standard {z(formen['standard'], 1)}; Dimension der Rückkehr {z(ds)} ± {z(eds)}; "
            f"Rest gegen Gasanteil innerhalb der Galaxien {z(gs)} ± {z(egs)} dex; je Galaxie gegen Masse ρ = {z(eig['Masse'][0])}, "
            f"gegen Gasanteil ρ = {z(eig['Gasanteil'][0])}, gegen Hubble-Typ ρ = {z(eig['Typ'][0])}, gegen Flächenhelligkeit ρ = {z(eig['Flächenhelligkeit'][0])}. AME2020: kleinste Masse je Nukleon A = {ame['m_min'][1]}, größte Bindung "
            f"A = {ame['b_max'][1]}, Spaltbarkeit am Maximum der Tröpfchenformel {z(ame['x_max'], 3)}. Die Daten tragen log g_bar, log g_obs je Punkt.")
    return dict(titel="Daten: SPARC und AME2020", welt=Welt([Koerper("Punkt der Beschleunigungsrelation", 1.0, [0, 0, 0], [0, 0, 0])]),
                T=np.arange(len(gb), dtype=float), X=X, text=text,
                pruefungen=[("a₀ aus SPARC gegen a₀ aus dem Atmen (Toleranz 1,5 σ)", K.a0_aus_atmen() / 1e-10, a0s / 1e-10, "10⁻¹⁰ m/s²",
                             1.5 * s / K.a0_aus_atmen()),
                            ("Dimension der Rückkehr (Toleranz 2 σ)", 1.0, ds, "", 2 * eds),
                            ("Bose-Form vor Standard-Form", 1.0, float(formen["standard"] > formen["Bose"]), "", 0),
                            ("kleinste Masse je Nukleon bei A", 56, ame["m_min"][1], "", 0),
                            ("größte Bindung je Nukleon bei A", 62, ame["b_max"][1], "", 0),
                            ("Spaltbarkeit am Maximum der Bindung", 0.25, ame["x_max"], "", 0.05),
                            ("Bindung des Wasserstoffs ist die halbe Spreizung α²/2 (PubChem)", K.alpha**2 / 2 * K.m_p / (K.m_e + K.m_p),
                             D.pubchem()["H"]["IE"] * K.e / (K.m_e * K.c**2), "", 1e-3)]
                + [(f"Rest je Galaxie unabhängig von {g} (p > 0,05)", 1.0, float(eig[g][1] > 0.05), "", 0)
                   for g in ("Masse", "Gasanteil", "Typ", "Flächenhelligkeit")])


ALLE = [atom, sonnensystem, merkur, galaxie, labor_kette, kosmos, daten]
