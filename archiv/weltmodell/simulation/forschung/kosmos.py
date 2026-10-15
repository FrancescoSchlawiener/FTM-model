# -*- coding: utf-8 -*-
"""
Kosmos (Ebene G, Brücke Z–G), Rechnungen, noch nicht im Modell: Lauf des eigenen Atmens Ẋ_an(t).

Kette (03, 05): Bindung macht Fluchtmoden inaktiv; je Bereich zählt das Maß der ganzen Tori in der Naht R_ta gegen
Newton am selben Rand (kern.fluchtmass_bereich), je Baryon im größten Bereich (geschachtelt), ab der Übergabe bei der
ersten Faltung (t_F = 1,8·t_ta). 𝒫 erhält das Maß: Die inaktiven Richtungen laufen durch die innere Richtung (Dichte
erster Rückkehrzeiten kern.kern_zeit, τ_D = 2Λ/c) und kommen als Atmen der Raumschleifen an, H_eff = H_b + ⅓·Ẋ_an
(kern.atemrate_eff). Das Atmen wirkt über H auf Zeit, Wachstum und Umkehr zurück (Varianten).
Blick: X ist ein Maß des Zustandsraums je Bereich und je Baryon (Z, G des Bereichs), keine Dichte am Ort; der Übergang
läuft in I (Rückkehrzeit, keine Strecke); angekommen ist er Atmen des Ganzen (G). Epoche überall als Eingang (a, t).
Grundlage: Forschungsrechnungen de_fluchtmoden_mass (noch nicht im Modell), de_fluchtmoden_verzoegert (trägt nicht in
der Menge), de_rest_anstieg, zelle_form. σ(M) aus CAMB mit Messwerten [G] (A_s, n_s, τ_re, Σm_ν); camb nur hier, erst
beim Aufruf geladen.
Laufzeit: CAMB ≈ 15–20 s einmal je Prozess, V0 und S1 wenige Sekunden, S2 mit Umkehr im eigenen Atmen je Runde einige
Sekunden. Jede Formel prüft zuerst die Art ihrer Eingänge (Proben mit Zahlen enden sofort).
"""
import math

import numpy as np

from kern.konstanten import *
from kern.formeln import formel
from kern.eingabe import EINGABEN
from kern.kosmos import (E_hintergrund, hintergrund_modell, rho_m_mitbewegt, d_LSS, zellgrenze, weltzeit,
                         atemrate_eff)
from kern.naht_haufen import fluchtmass_bereich
from kern.rueckkehr import kern_zeit

VARIANTEN = ("V0", "S1", "S2")
_SPEKTRUM = {}


@formel('streuung_masse', 'σ(M) = [∫d ln k k³P(k)/(2π²)·W²(kR)]^{1/2}, W = 3(sin x − x cos x)/x³, R = (3M/(4πρ_m))^{1/3}',
        ('σ_M',), ('M', 'ω_b', 'ω_c', 'H₀', 'A_s', 'n_s', 'τ_re', 'Σm_ν'), 'G', 'forschung',
        'Rechnung de_fluchtmoden (σ(M) aus CAMB)', None)
def streuung_masse(lg_M):
    """Lineare Streuung der Dichte heute auf der Masse M_ges = 10^lg_M M☉ (Feld), Kugelfilter, Leistungsspektrum aus
    CAMB mit ω_b, ω_c, H₀ und den Messwerten A_s, n_s, τ_re, Σm_ν [G]; das Spektrum wird je Prozess einmal gerechnet."""
    if not (isinstance(lg_M, np.ndarray)):
        raise TypeError("streuung_masse: Feld der log₁₀ M_ges erwartet")
    if "P0" not in _SPEKTRUM:
        import camb
        from kern.messwert import MESSWERTE as M
        p = camb.CAMBparams()
        p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=M["M_NU"].wert, tau=M["TAU_RE"].wert)
        p.InitPower.set_params(As=M["A_S"].wert, ns=M["N_S"].wert)
        p.set_matter_power(redshifts=[0.0], kmax=400)
        pk = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False, kmax=400,
                                                var1='delta_tot', var2='delta_tot', zmax=1)
        _SPEKTRUM["k"] = np.exp(np.linspace(math.log(1e-4), math.log(400.0), 5000))
        _SPEKTRUM["P0"] = pk.P(0.0, _SPEKTRUM["k"])
    k, P0 = _SPEKTRUM["k"], _SPEKTRUM["P0"]
    aus = []
    for lg in lg_M:
        R = (3 * 10**lg * M_sun / (4 * math.pi * rho_m_mitbewegt()))**(1 / 3) / Mpc
        x = k * R
        W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
        aus.append(math.sqrt(np.trapezoid(k**3 * P0 / (2 * math.pi**2) * W**2, np.log(k))))
    return np.array(aus)


@formel('wachstum_linear', "D″ + (2 + d ln E/d ln a)D′ = ³⁄₂Ω_m a⁻³/E²·D, D = D′ = a bei a = 10⁻³ (′ = d/d ln a)",
        ('D_+',), ('E(z)', 'Ω_m'), 'G', 'forschung', 'Rechnung de_rest_anstieg', None)
def wachstum_linear(E, a_ende=3.5):
    """Lineares Wachstum im Hintergrund E(a) (Funktion von a); gibt die Funktion a ↦ (D, D′)."""
    if not (callable(E)):
        raise TypeError("wachstum_linear: Hintergrund E(a) als Funktion erwartet")
    from scipy.integrate import solve_ivp

    def rechts(x, y):
        a, h = math.exp(x), 1e-5
        dl = (math.log(E(a * math.exp(h))) - math.log(E(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / E(a)**2 * y[0]]
    s = solve_ivp(rechts, (math.log(1e-3), math.log(a_ende)), [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-13)
    return lambda a: s.sol(math.log(a))


@formel('kugel_umkehr', 'R″ = −GM/R² + A(a)·R (Kugel im Hintergrund E, Beschleunigung A der Schleifen); gibt je a_ta '
        '(1 + Δ_ta, δ_ta) mit δ_ta = δ_i·D(a_ta)/a_i',
        ('Δ_ta', 'δ_ta'), ('ln a', 'E(z)', 'Ω_m'), 'LG', 'forschung', 'Rechnung de_rest_anstieg (Kugelmodell)', None,
        alternative='dichte_umkehr_heute')
def kugel_umkehr(a_ta, E, beschleunigung, n_start=90):
    """Umkehr einer Kugel, die bis zur Umkehr mit den Schleifen des Ganzen mitatmet (θ̃ ≠ 0): R̈ = −GM/R² + A(a)·R in
    Einheiten H₀ = 1, Start a_i = 10⁻³ im Hubble-Fluss (1 − δ_i/3). Je Anfangskontrast δ_i ein Lauf in ln a bis Ṙ = 0;
    a_ta(δ_i) wird umgekehrt (Spline in ln). Gibt Feld der Zeilen (Dichte der Kugel gegen ρ̄_m bei der Umkehr, lineare
    Schwelle δ_ta = δ_i·D(a_ta)/a_i mit dem Wachstum im selben Hintergrund)."""
    if not (callable(E) and callable(beschleunigung)):
        raise TypeError("kugel_umkehr: E(a) und A(a) als Funktionen erwartet")
    from scipy.integrate import solve_ivp
    from scipy.interpolate import CubicSpline
    ai = 1e-3
    D = wachstum_linear(E)

    def lauf(di):
        GM = 0.5 * OMEGA_M * (1 + di)

        def rechts(x, y):
            a = math.exp(x)
            H = E(a)
            return [y[1] / H, (-GM / y[0]**2 + beschleunigung(a) * y[0]) / H]
        ev = lambda x, y: y[1]
        ev.terminal, ev.direction = True, -1
        s = solve_ivp(rechts, (math.log(ai), math.log(20.0)), [ai, E(ai) * ai * (1 - di / 3)], events=ev, rtol=1e-10,
                      atol=1e-14)
        return (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (None, None)
    zeilen = []
    for ld in np.linspace(-4.3, 0.2, n_start):
        x, R = lauf(10**ld)
        if x is not None:
            zeilen.append((x, ld, math.log(R)))
    x, ld, lR = map(np.array, zip(*sorted(zeilen)))
    sd, sR = CubicSpline(x, ld), CubicSpline(x, lR)
    aus = []
    for a in np.atleast_1d(a_ta):
        di, Rt = 10**float(sd(math.log(a))), math.exp(float(sR(math.log(a))))
        aus.append(((1 + di) * (a / Rt)**3, di * D(a)[0] / ai))
    return np.array(aus)


@formel('fluchtmass_je_baryon', 'X(a) = ∫(−dF/d ln M)·ln B(M, a) d ln M, F = erfc(δ_ta(a)/(√2σ(M)D(a))), M ≥ M_min '
        '(jedes Baryon im größten gebundenen Bereich)',
        ('X',), ('ln a', 'σ_M', 'δ_ta', 'X_b', 'D_+'), 'ZG', 'forschung', 'Rechnung de_fluchtmoden_mass', None)
def fluchtmass_je_baryon(a, lg_M, sigma, delta_ta, D, lnB, lg_M_min):
    """Inaktives Fluchtmaß je Baryon zur Epoche a: Massenanteil in Bereichen, die bis a umgekehrt sind (erfc mit der
    linearen Schwelle δ_ta und σ(M)·D(a)), gewichtet mit dem Maß je Bereich ln B (Regel gegen Newton in der Naht)."""
    if not (isinstance(lg_M, np.ndarray) and isinstance(lnB, np.ndarray)):
        raise TypeError("fluchtmass_je_baryon: Felder erwartet")
    from scipy.special import erfc
    F = erfc(delta_ta / (math.sqrt(2) * sigma * D))
    dF = -np.gradient(F, lg_M * math.log(10))
    m = lg_M >= lg_M_min
    return float(np.trapezoid((dF * lnB)[m], lg_M[m] * math.log(10)))


@formel('ankunft_rueckkehr', 'Ẋ_an(t) = ∫S(t′)·K(t − t′)dt′, S = dX/dt zur Zeit t/(t_F/t_ta), K Dichte erster Rückkehrzeiten; '
        'X_an = ∫Ẋ_an dt',
        ('Ẋ_an', 'X_an'), ('X', 't_F/t_ta', 'τ_D'), 'IG', 'forschung', 'Rechnung de_fluchtmoden_verzoegert', None)
def ankunft_rueckkehr(t_ta, X, t_F_zu_t_ta, n=3000):
    """Verzögerte Ankunft: Das bei der Umkehr (Zeit t_ta, Feld in Gyr) inaktiv gewordene Maß X geht ab der Übergabe
    t_F = (t_F/t_ta)·t_ta in die innere Richtung und kommt mit der Dichte erster Rückkehrzeiten an (Γ = 0, θ̃ = 0 im
    gebundenen Bereich); der Bestand vor dem ersten Gitterpunkt als Puls. Gibt (Zeitgitter [Gyr], Ẋ_an [1/Gyr], X_an)."""
    if not (isinstance(t_ta, np.ndarray) and isinstance(X, np.ndarray)):
        raise TypeError("ankunft_rueckkehr: Felder erwartet")
    from scipy.integrate import cumulative_trapezoid
    from scipy.interpolate import CubicSpline
    K_Gyr = lambda T: np.array([kern_zeit(x * Gyr, 0.0) * Gyr if x > 0 else 0.0 for x in np.atleast_1d(T)])
    TG = np.linspace(t_ta[0], t_ta[-1], n)
    dt = TG[1] - TG[0]
    Ki = K_Gyr(np.arange(n) * dt)
    sp = CubicSpline(t_ta, X)
    tq = TG / t_F_zu_t_ta
    S = sp(np.clip(tq, t_ta[0], None), 1) * (tq >= t_ta[0]) / t_F_zu_t_ta
    t_puls = t_ta[0] * t_F_zu_t_ta
    R = np.array([np.trapezoid(Ki[j::-1] * S[:j + 1], dx=dt) for j in range(n)])
    R += X[0] * K_Gyr(TG - t_puls) * (TG > t_puls)
    return TG, R, cumulative_trapezoid(R, TG, initial=0.0)


@formel('lauf_eigenes_atmen', 'H_eff = H_b + ⅓Ẋ_an, Ẋ_an aus X(a) über die Rückkehr; V0: Bildung im Vergleichshintergrund, '
        'S1: Zeit und Wachstum im eigenen Atmen (iteriert), S2: dazu Umkehr im eigenen Atmen',
        ('H_eff',), ('ln a', 'X', 'Ẋ_an', 'H_b', 'τ_D', 't_F/t_ta', 'L_min'), 'G', 'forschung',
        'Rechnungen de_rest_anstieg, zelle_form', None, alternative='atemrate_eff')
def lauf_eigenes_atmen(variante, runden=8, grenze=2e-5):
    """Lauf des eigenen Atmens. Gibt dict: E_eff (a ↦ H_eff/H₀), t0 [Gyr], X0 (inaktiv heute), X_an0 (angekommen heute),
    anteil (X_an/X heute), Omega_X0 = E_eff(1)² − E_b(1)² (als Energie gelesen), dLSS_faktor und dLSS (Durchmesser der
    letzten Streufläche für den Beobachter heute, ρ_b fest), grenze_wuerfel ((1 + Q)·a₀-Grenze der Zelle mit L_min =
    0,985·d_LSS), runden, aenderung (letzte max|ΔE/E|)."""
    if not (variante in VARIANTEN):
        raise TypeError(f"lauf_eigenes_atmen: Variante aus {VARIANTEN}")
    from scipy.interpolate import CubicSpline
    H0g = H0_si() * Gyr
    Om, Or, NK = OMEGA_M, OMEGA_R, 1 + Q_KOSM
    tF = EINGABEN["T_F_ZU_T_TA"].wert
    lg_min = math.log10(EINGABEN["FLUCHT_M_MIN"].wert)
    E_L = lambda a: E_hintergrund(1 / a - 1)
    E_b = lambda a: hintergrund_modell(1 / a - 1) / H0_si()
    # Gitter: Epochen, Massen (M_ges in M☉) und Umkehr im Vergleichshintergrund
    AA = np.exp(np.linspace(math.log(1 / 31), math.log(3.0), 56))
    LM = np.linspace(6.0, 16.5, 106)
    LMg = np.linspace(6.0, 16.5, 22)
    AK = np.array([1 / 31, 1 / 11, 1 / 7, 1 / 4, 1 / 3, 1 / 2.5, 1 / 2, 1 / 1.5, 0.8, 1.0, 1.4, 2.0, 3.0])
    OL = 1 - Om
    KL = kugel_umkehr(AK, lambda a: math.sqrt(Om * a**-3 + OL), lambda a: OL)         # Naht: ohne Strahlung
    KR = kugel_umkehr(AK, E_L, lambda a: (1 - Om - Or) - Or * a**-4)                    # Schwelle: mit Strahlung
    s_Delta, s_delta = CubicSpline(np.log(AK), np.log(KL[:, 0])), CubicSpline(np.log(AK), KR[:, 1])
    SIG = streuung_masse(LM)
    TB, TS = {}, {}
    for a in AA:
        b1, b2 = [], []
        for lg in LMg:
            Mb = 10**lg * M_sun / NK
            R = (3 * NK * Mb / (4 * math.pi * rho_m_mitbewegt() * a**-3 * math.exp(float(s_Delta(math.log(a))))))**(1 / 3)
            x1 = fluchtmass_bereich(Mb, n=4001, R_ta=R)
            b1.append(x1)
            b2.append((fluchtmass_bereich(Mb, n=4001, R_ta=1.05 * R) - x1) / math.log(1.05))
        TB[a], TS[a] = np.array(b1), np.array(b2)
    D_L1 = wachstum_linear(E_L)(1.0)[0]

    def D_M(Ef):
        zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, 1090.0, 4001)[1:]])
        return np.trapezoid(1 / np.array([Ef(1 / (1 + x)) for x in zs]), zs)

    def schritt(E_str, umkehr_eigen):
        TA = np.array([weltzeit(math.log(a), E=lambda u: E_str(math.exp(u))) for a in AA]) / Gyr
        g = wachstum_linear(E_str)
        dfak, rfak = np.ones(len(AA)), np.ones(len(AA))
        if umkehr_eigen:
            lg = np.linspace(math.log(1e-3), math.log(3.6), 2000)
            E2 = np.array([E_str(math.exp(x))**2 for x in lg])
            OX = E2 - Om * np.exp(-3 * lg) - Or * np.exp(-4 * lg) - (E_b(1.0)**2 - Om - Or)
            A_v = (E_b(1.0)**2 - Om - Or) + OX + 0.5 * np.gradient(OX, lg)
            E_v = np.sqrt(np.clip(E2 - Or * np.exp(-4 * lg), 1e-30, None))
            sA, sE = CubicSpline(lg, A_v), CubicSpline(lg, E_v)
            KK = kugel_umkehr(AK, lambda a: float(sE(math.log(a))), lambda a: float(sA(math.log(a))))
            dfak = np.interp(np.log(AA), np.log(AK), KK[:, 1] / KL[:, 1])
            rfak = np.interp(np.log(AA), np.log(AK), (KL[:, 0] / KK[:, 0])**(1 / 3))
        X = np.array([fluchtmass_je_baryon(a, LM, SIG, float(s_delta(math.log(a))) * dfak[j], g(a)[0] / D_L1,
                                           np.interp(LM, LMg, TB[a] + TS[a] * math.log(rfak[j])), lg_min)
                      for j, a in enumerate(AA)])
        TG, R, Xan = ankunft_rueckkehr(TA, X, tF)
        AG = np.exp(np.interp(TG, TA, np.log(AA)))
        E_neu = atemrate_eff(np.array([E_b(a) for a in AG]), R / H0g)            # H_eff/H₀ = E_b + ⅓Ẋ_an/H₀
        t0 = weltzeit(0.0, E=lambda u: E_str(math.exp(u))) / Gyr
        return dict(X=X, TG=TG, AG=AG, Xan=Xan, E=E_neu, t0=t0)

    def E_aus(d):
        AG, E = d["AG"], d["E"]
        return lambda a: float(np.interp(a, AG, E)) if a >= AG[0] else E_b(a)

    d = schritt(E_L, False)
    n_runden, aend = 0, 0.0
    if variante != "V0":
        probe = np.exp(np.linspace(math.log(AA[0]), math.log(2.0), 200))
        for n_runden in range(1, runden + 1):
            E_alt = E_aus(d)
            d = schritt(E_alt, variante == "S2")
            E_neu = E_aus(d)
            aend = max(abs(E_neu(a) / E_alt(a) - 1) for a in probe)
            if aend < grenze:
                break
    Ef = E_aus(d)
    faktor = D_M(Ef) / D_M(E_L)
    X0 = float(np.interp(0.0, np.log(AA), d["X"]))
    Xan0 = float(np.interp(d["t0"], d["TG"], d["Xan"]))
    dL = d_LSS() * faktor
    return dict(E_eff=Ef, t0=d["t0"], X0=X0, X_an0=Xan0, anteil=Xan0 / X0, Omega_X0=Ef(1.0)**2 - E_b(1.0)**2,
                dLSS_faktor=faktor, dLSS=dL, grenze_wuerfel=zellgrenze(1.0, KREISPAAR_GRENZE * dL), runden=n_runden,
                aenderung=aend)


@formel('winkelmoden_zelle', 'C_ℓ/C_ℓ,∞ = 2ℓ(ℓ + 1)·[Σ_k 2π²j_ℓ²(kχ)/(V k³) + ∫_{kmax χ}^∞ j_ℓ²(x)dx/x], k ∈ 2π·Γ* (Sachs-Wolfe, '
        'skaleninvariant); Kovarianz ⟨a_ℓm a_ℓm′⟩ im Achsensystem der Zelle',
        ('C_ℓ/C_ℓ∞',), ('Gitter', 'L', 'd_LSS', 'ℓ_W'), 'G', 'forschung', 'Rechnung zelle_form', None)
def winkelmoden_zelle(gitter, V, chi, lmax=20, x_max=150.0):
    """Winkelmoden der Hintergrundstrahlung in der geschlossenen Zelle (Gitter aus kern.GITTER_ZELLE, Volumen V,
    Abstand χ = d_LSS/2 des Beobachters heute): nur Geometrie (Sachs-Wolfe, skaleninvariant, ohne spätes Atmen). Die
    fehlenden Moden k < k₁ senken die großen Winkel (Schließung, keine Dichte). Gibt dict: C (ℓ ↦ C_ℓ/C_ℓ∞),
    S_halb (∫_{−1}^{1/2}C(θ)²dcosθ gegen den unendlichen Raum), M2, M3 (Kovarianz ℓ = 2, 3 im Achsensystem, reelle Y,
    m = −ℓ … ℓ, Spur 2ℓ + 1)."""
    from kern.kosmos import GITTER_ZELLE
    if gitter not in GITTER_ZELLE:
        raise TypeError(f"winkelmoden_zelle: Gitter aus {tuple(GITTER_ZELLE)}")
    from numpy.polynomial import legendre as Lg
    from scipy.integrate import quad
    from scipy.special import spherical_jn, sph_harm_y
    A = np.array(GITTER_ZELLE[gitter], dtype=float)
    A = A / abs(np.linalg.det(A))**(1 / 3)
    B = 2 * math.pi * np.linalg.inv(A).T / V**(1 / 3)
    kmax = x_max / chi
    N = int(kmax / np.min(np.linalg.norm(B, axis=1)) * 1.8) + 2
    r = np.arange(-N, N + 1)
    n = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    kv = n @ B
    k = np.linalg.norm(kv, axis=1)
    m = (k > 0) & (k < kmax)
    kv, k = kv[m], k[m]
    x, w = k * chi, 2 * math.pi**2 / (V * k**3)
    rest = lambda l: quad(lambda y: spherical_jn(l, y)**2 / y, kmax * chi, np.inf, limit=400)[0]
    C = {l: 2 * l * (l + 1) * (np.sum(w * spherical_jn(l, x)**2) + rest(l)) for l in range(2, lmax + 1)}
    th, ph = np.arccos(np.clip(kv[:, 2] / k, -1, 1)), np.arctan2(kv[:, 1], kv[:, 0])
    kov = {}
    for l in (2, 3):
        Y = []
        for mm in range(-l, l + 1):
            y = sph_harm_y(l, abs(mm), th, ph)
            Y.append(math.sqrt(2) * (-1)**mm * (y.imag if mm < 0 else y.real) if mm != 0 else y.real)
        Y = np.array(Y)
        Mk = (Y * (w * spherical_jn(l, x)**2)) @ Y.T * 4 * math.pi / (2 * l + 1) + np.eye(2 * l + 1) * rest(l) / (2 * l + 1)
        kov[l] = Mk / np.trace(Mk) * (2 * l + 1)
    mu = np.linspace(-1, 0.5, 3001)

    def korrelation(Rf):
        ko = np.zeros(lmax + 1)
        for l in range(2, lmax + 1):
            ko[l] = (2 * l + 1) / (4 * math.pi) * Rf(l) / (l * (l + 1))
        return Lg.legval(mu, ko)
    S_halb = np.trapezoid(korrelation(lambda l: C[l])**2, mu) / np.trapezoid(korrelation(lambda l: 1.0)**2, mu)
    return dict(C=C, S_halb=float(S_halb), M2=kov[2], M3=kov[3])
