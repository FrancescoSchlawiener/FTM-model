# -*- coding: utf-8 -*-
"""
Naht und Haufen (03; Brücke L–G): R_Q, wo das örtliche q das kosmische erreicht, M_b,Q, Umkehr heute (Top-Hat mit
Cache _UMKEHR), Umkehrradius, Dauer des Übergangs von w vor der Umkehr, Fluchtmaß je Bereich, Fenster von Q.
Blickrichtung: R_Q, R_ta sind Radien im Raum (L), gesetzt vom Zustand des Ganzen (G: q_kosm, Hintergrund heute);
das Fluchtmaß zählt Tori im Zustandsraum (Z) innerhalb der Naht.
"""
from .konstanten import *
from .formeln import formel
from .faltung import r_stern
from .rueckkehr import nu, gewicht_w
from .kosmos import rho_m_mitbewegt, E_hintergrund, zellgrenze, gewicht_hintergrund, theta_


@formel('r_Q', 'R_Q = r*/(−ln q_kosm)',
        ('R_Q',), ('M', 'a₀', 'q_kosm'), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'R_Q = r*/(−ln q)'))
def r_Q(M, a0=A0, qk=q_KOSM):
    """R_Q = r*/(−ln q_kosm): dort erreicht das örtliche q das kosmische, a_N/a₀ = (ln q)²."""
    return r_stern(M, a0) / (-math.log(qk))


@formel('masse_bQ', 'M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴)',
        ('M_b,Q',), ('a₀', 'H', 'Q'), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴)'))
def masse_bQ(a0=A0, H=None, Qk=Q_KOSM):
    """M_b,Q = (1 + Q)²(ln q)⁶·a₀³/(10⁴·G·H⁴): Baryonenmasse, bei der R_Q die Grenze R₂₀₀ erreicht."""
    H = H0_si() if H is None else H
    lq = math.log(Qk / (1 + Qk))
    return (1 + Qk)**2 * lq**6 * a0**3 / (1e4 * G * H**4)


_UMKEHR = {}


@formel('dichte_umkehr_heute', '1 + Δ_ta aus Top-Hat im flachen Hintergrund mit Vakuumteil 1 − Ω_m',
        ('Δ_ta',), ('Ω_m',), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Bis zum Umkehrradius (θ̃ = 0) ist w = 0'))
def dichte_umkehr_heute(Om=OMEGA_M):
    """1 + Δ_ta: mittlere Dichte (gegen ρ_m) einer Kugelschale, die heute umkehrt (Top-Hat im flachen Hintergrund
    E(z) mit Ω_r = 0 und Vakuumteil 1 − Ω_m, Start a = 10⁻³ mit Hubble-Fluss (1 − δ/3)). Unabhängig von der Masse."""
    if Om in _UMKEHR:
        return _UMKEHR[Om]
    from scipy.integrate import solve_ivp, quad
    from scipy.optimize import brentq
    H0, OL = H0_si(), 1 - Om
    Hb = lambda a: H0 * E_hintergrund(1 / a - 1, Om, 0.0)
    t_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
    t0, ai = t_a(1.0), 1e-3
    ti = t_a(ai)

    def lauf(di):   # Einheit: Radius ri = 1, GM = (H_i²/2)·Ω_m(a_i)·(1 + δ)
        Hi = Hb(ai); GM = 0.5 * Hi**2 * (Om / ai**3) / (Om / ai**3 + OL) * (1 + di)
        f = lambda t, y: [y[1], -GM / y[0]**2 + OL * H0**2 * y[0]]
        ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
        sol = solve_ivp(f, (ti, 2 * t0), [1.0, Hi * (1 - di / 3)], events=ev, rtol=1e-9, atol=1e-12)
        return (sol.t_events[0][0], sol.y_events[0][0][0]) if len(sol.t_events[0]) else (2 * t0, None)
    ld = brentq(lambda x: lauf(10**x)[0] - t0, -4, 0, xtol=1e-10)
    di = 10**ld; R = lauf(di)[1]
    _UMKEHR[Om] = (1 + di) / (R**3 * ai**3)     # Masse/Volumen gegen ρ_m heute: (1 + δ)·ρ_m(a_i)·r_i³/R³
    return _UMKEHR[Om]


@formel('umkehrradius_heute', 'R_ta = (3M_ges/(4πρ_m(1 + Δ_ta)))^{1/3}',
        ('R_ta',), ('M_ges', 'Δ_ta', 'ρ_m'), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', '319 / 580 / 1482 kpc'))
def umkehrradius_heute(M_ges):
    """R_ta = (3M/(4πρ_m(1 + Δ_ta)))^{1/3}, M_ges = (1 + Q)·M_b (Maß des Ganzen)."""
    return (3 * M_ges / (4 * math.pi * rho_m_mitbewegt() * dichte_umkehr_heute()))**(1 / 3)


@formel('schalen_dauer', 'Dauer von w = 0,9 bis w = 0,1 vor der Umkehr: Ṙ² = 2μ(1/R − 1/R_ta) − Ω_Λ(R_ta² − R²), '
        'θ = 3Ṙ/R, w = 1 − e^{−√(θ/θ_c)}, Kollaps (R = 0,01·R_ta) bei z_c',
        ('t_w',), ('z_c', 'a₀', 'Ω_m', 'H₀'), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'fällt w vor dem Umkehrpunkt von 0,9 auf 0,1'))
def schalen_dauer(z_c, a0=A0, w_oben=0.9, w_unten=0.1):
    """Zeit in s, in der das Gewicht einer Schale vor ihrer Umkehr von w_oben auf w_unten fällt (Preisach: θ̃ = θ bis
    zur Umkehr). Energieintegral der Schale im flachen Hintergrund mit Vakuumteil 1 − Ω_m ohne Strahlung, Einheit
    1/H₀; die Schale kollabiert (R = 0,01·R_ta) bei z_c."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    H0, Om = H0_si(), OMEGA_M
    OL, MU, thc = 1 - Om, Om / 2, a0 / c / H0_si()
    x_oben, x_unten = math.log(1 - w_oben)**2, math.log(1 - w_unten)**2      # θ/θ_c aus w = 1 − e^{−√x}

    def t_von_a(a):
        return 2 / (3 * math.sqrt(OL)) * math.asinh(math.sqrt(OL / Om) * a**1.5)

    def Rdot(R, Rta):
        return math.sqrt(max(2 * MU * (1 / R - 1 / Rta) - OL * (Rta**2 - R**2), 0.0))

    def t_bis(s1, s2, Rta):                  # Zeit von s1·R_ta bis s2·R_ta, Substitution s = 1 − u²
        f = lambda u: 2 * u * Rta / max(Rdot(Rta * (1 - u * u), Rta), 1e-300)
        return quad(f, math.sqrt(1 - s2), math.sqrt(1 - s1), limit=400, epsrel=1e-11)[0]

    def theta_s(s, Rta):
        return 3 * Rdot(s * Rta, Rta) / (s * Rta)

    tc = t_von_a(1 / (1 + z_c))
    Rta = brentq(lambda R: 2 * t_bis(0.0, 1.0, R) - t_bis(0.0, 0.01, R) - tc, 1e-3, 1.5, xtol=1e-13)
    s_o, s_u = (brentq(lambda s: theta_s(s, Rta) / thc - x, 1e-6, 1 - 1e-15, xtol=1e-15) for x in (x_oben, x_unten))
    return t_bis(s_o, s_u, Rta) / H0


@formel('gewicht_schale_aussen', 'w = 1 − e^{−√(θ/θ_c)}, θ = 3Ṙ/R heute aus Ṙ² = 2μ/R + Ω_Λ R² − K mit Alter t₀, '
        '2μ/Ω_m = (1 + Δ_ta)(R_ta/R)³ (Punktmasse, Einheiten H₀ = 1, R = 1 heute)',
        ('w',), ('R/R_ta', 'Δ_ta', 'a₀', 'Ω_m', 'H₀'), 'LG', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Außerhalb steigt w rasch'), alternative='gewicht_w')
def gewicht_schale_aussen(R_zu_Rta, a0=A0):
    """Gewicht einer Schale knapp außerhalb der Naht, die heute noch ausdehnt (θ̃ = θ > 0, Preisach): Ihre Ausdehnungsrate
    folgt aus dem Energieintegral im flachen Hintergrund mit Vakuumteil 1 − Ω_m ohne Strahlung (Alter t₀ wie der
    Hintergrund), die mittlere Dichte innen aus der Masse, die bei R_ta heute umkehrt. Nach außen nähert w den Wert des
    mitatmenden Hintergrunds (θ̃ = 3H)."""
    if not R_zu_Rta > 1:
        raise ValueError("gewicht_schale_aussen: nur außerhalb der Naht (R/R_ta > 1)")
    from scipy.integrate import quad
    from scipy.optimize import brentq
    Om = OMEGA_M
    OL = 1 - Om
    T0 = 2 / (3 * math.sqrt(OL)) * math.asinh(math.sqrt(OL / Om))
    mu = dichte_umkehr_heute(Om) * R_zu_Rta**-3 * Om / 2

    def alter(K):
        return quad(lambda R: 1 / math.sqrt(max(2 * mu / R + OL * R * R - K, 1e-300)), 0, 1, limit=400, epsrel=1e-12)[0]
    K = brentq(lambda K: alter(K) - T0, -50.0, 2 * mu + OL - 1e-12, xtol=1e-14)
    return gewicht_w(3 * math.sqrt(2 * mu + OL - K) * H0_si(), a0)


@formel('rand_schale_aussen', 'R/R_ta, an dem gewicht_schale_aussen den Wert w erreicht (Umkehrung)',
        ('R/R_ta',), ('w', 'Δ_ta', 'a₀', 'Ω_m', 'H₀'), 'LG', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'w wechselt in einer Schale bis 1,16·R_ta'))
def rand_schale_aussen(w_wert, a0=A0):
    """Radius (in R_ta), ab dem eine heute ausdehnende Schale außerhalb der Naht das Gewicht w_wert erreicht."""
    from scipy.optimize import brentq
    if not 0 < w_wert < gewicht_hintergrund(theta_(a0, H0_si())):
        raise ValueError("rand_schale_aussen: w zwischen 0 und dem Wert des mitatmenden Hintergrunds")
    return brentq(lambda x: gewicht_schale_aussen(x, a0) - w_wert, 1 + 1e-9, 20.0, xtol=1e-12)


@formel('uebergabe_vor_umkehr', 'I_vor = ∫₀^{R_ta}(1 − w)·ω dR/Ṙ, ω² = GM/R³, θ = 3Ṙ/R (Schale, die heute umkehrt)',
        ('I_vor',), ('Δ_ta', 'a₀', 'Ω_m', 'H₀'), 'LG', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Die Übergabe füllt die Schale bis heute nicht'))
def uebergabe_vor_umkehr(a0=A0):
    """Übergabeintegral einer Schale, die heute umkehrt, von Anfang bis zur Umkehr: Solange sie mitatmet (w ≈ 1), geht
    nichts über; erst kurz vor der Umkehr fällt w. Unabhängig von der Masse (Top-Hat); Einheiten H₀ = 1, R_ta = 1."""
    from scipy.integrate import quad
    Om = OMEGA_M
    OL = 1 - Om
    mu = dichte_umkehr_heute(Om) * Om / 2
    thc = a0 / c / H0_si()

    def f(u):                                   # R = 1 − u², dR = 2u du
        R = 1 - u * u
        Rd = math.sqrt(max(2 * mu * (1 / R - 1) - OL * (1 - R * R), 1e-300))
        return 2 * u * math.exp(-math.sqrt(3 * Rd / R / thc)) * math.sqrt(mu / R**3) / Rd
    return quad(f, 0, 1, limit=400, epsrel=1e-10)[0]


@formel('fluchtmass_bereich', 'X_b = ln(μ_Regel/μ_Newton), μ = ∫4πr²·(4π/3)v_R³/(1 − r²/R_ta²)dr, v_R² = 2∫_r^{R_ta} a dr′',
        ('X_b',), ('M', 'a₀', 'R_ta'), 'LZ', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Milchstraße ln 284 = 5,65'))
def fluchtmass_bereich(M_b, a0=A0, n=200001, R_ta=None):
    """Maß der ganzen Tori innerhalb der Naht R_ta (Umkehrradius der Masse (1 + Q)·M_b) nach der Regel gegen Newton am
    selben Rand, als ln (Punktmasse). Ein Torus liegt in R_ta, wenn v_r² + v_t²(1 − r²/R²) ≤ v_R²(r). Newton analytisch:
    μ_N = (16π²/3)(2GM_bR)^{3/2}·π(3/2 − √2). Ohne R_ta die Naht eines heute umkehrenden Bereichs; für einen Bereich,
    der zu einer anderen Epoche umkehrt, dessen Naht R_ta."""
    R = umkehrradius_heute((1 + Q_KOSM) * M_b) if R_ta is None else R_ta
    GM = G * M_b
    # Gitter in x = r/R: dicht bei 0 und bei 1 (Wurzelverhalten (1 − x)^{1/2} des Integranden)
    x = np.unique(np.concatenate([np.geomspace(1e-9, 0.5, n // 2), 1 - np.geomspace(1e-12, 0.5, n // 2)]))
    aN = GM / (x * R)**2
    a = nu(aN, a0) * aN
    from scipy.integrate import cumulative_trapezoid
    I = cumulative_trapezoid(a[::-1], -x[::-1] * R, initial=0.0)[::-1]      # ∫_r^R a dr′
    v2 = 2 * I
    v2[-1] = 0.0
    f = x**2 * v2**1.5 / ((1 - x) * (1 + x))
    mu_R = (16 * math.pi**2 / 3) * R**3 * np.trapezoid(f, x)
    mu_N = (16 * math.pi**2 / 3) * (2 * GM * R)**1.5 * math.pi * (1.5 - math.sqrt(2))
    return math.log(mu_R / mu_N)


@formel('q_fenster', 'Q ≥ 5 aus r*(M_b,Zelle)/Λ = √(6/ν_kosm) ≤ 1; Q ≤ (1 + Q)a₀-Grenze/a₀ − 1',
        ('Q_min', 'Q_max'), ('a₀', 'L_min'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'r*(M_b) = Λ·√(6/ν_kosm) ≤ Λ, also Q ≥ 5'))
def q_fenster(a0=A0, L_min=None):
    """Fenster der kosmischen Rückkehrzahl: unten r*(M_b,Zelle)/Λ = √(6/ν) ≤ 1 ⇒ Q ≥ 5; oben die Zellgrenze
    (1 + Q)a₀ ≤ 3c⁴/(2Gρ_bL_min³) (Würfel) mit L ≥ 0,985·d_LSS (keine Kreispaare)."""
    return 5.0, zellgrenze(1.0, L_min) / a0 - 1


@formel('haufenstoss_trennzeit', 't_T = I_F/ω_orb (Übergabe beim Gas erst nach I_F Radiant seit dem Stoß)',
        ('t_T',), ('I_F', 'ω_orb'), 'L', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'ω_orb·t ≈ I_F'))
def haufenstoss_trennzeit(I_F, omega_orb):
    """Zeit, bis die Bahnen seit dem Stoß I_F Radiant durchlaufen haben und die Mode beim Gas übergeht."""
    return I_F / omega_orb


_UEBERGABE = {}


@formel('uebergabe_faltung', 'I_F = ∫₁^{τ_F} ω dτ, ω² = (π²/8)/Λ³ (Masse innen fest bis zur Faltung); Λ″ = −(π²/8)τ^p𝓜(Λτ^{−s})/Λ² '
        '+ (π²/8)ε_j/Λ³, s = 2/(9e) + 2/3, p = 2/(3e), 𝓜 aus der Bahn selbst; τ_F: Λτ^{−s} = λ_c (äußerster Rücklauf)',
        ('I_F',), ('e_M', 'ε_j'), 'L', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'I_F = ∫ω dt ≈ 1,1–1,25 Radiant nach der Umkehr'))
def uebergabe_faltung(e, eps_j, n_t=1500, n_l=600, rtol=1e-7, grenze=1e-3, n_max=60, misch=0.7):
    """Übergabe der Mode ξ bei der ersten Faltung des Phasenraumblatts am Ort der Schale (selbstähnlicher Einfall im
    mitatmenden Hintergrund, eine Bahn für alle Schalen): Anfangsprofil δM/M ∝ M^{−e} (Keim e = 1), Drehimpuls ab der
    Umkehr j² = ε_j·GM·r_ta. Zeit in Einheiten der Umkehrzeit, Radius in r_ta(t) ∝ t^s. Das Massenprofil 𝓜(λ) folgt
    aus der Bahn selbst (Fixpunkt, gemischt mit dem Faktor misch); gefaltet ist die Schale, wenn ihr λ den äußersten
    Rücklauf λ_c nach dem ersten Perizentrum erreicht. Gibt (I_F in Radiant, t_F/t_ta, r_F/R_ta). Ebene: Bahnen im
    gebundenen Bereich (w = 0), die Übergabe ändert die Form der Mode, nicht ihr Maß."""
    if not (0.1 <= e <= 1.5 and 1e-3 <= eps_j <= 0.5):
        raise ValueError("uebergabe_faltung: Profil 0,1 ≤ e ≤ 1,5 und Drehimpuls 10⁻³ ≤ ε_j ≤ 0,5 (gebundene Bahnen)")
    schluessel = (e, eps_j, n_t, n_l, rtol, grenze, n_max, misch)
    if schluessel in _UEBERGABE:
        return _UEBERGABE[schluessel]
    from scipy.integrate import solve_ivp, quad
    from scipy.optimize import brentq, minimize_scalar
    K8 = math.pi**2 / 8
    p, s = 2 / (3 * e), 2 / (9 * e) + 2 / 3
    tmax = max(8.0, 0.15 ** (-1 / p))              # Schalen jenseits (Masse ≤ 15 % von M_ta) als Zentrum
    tg = np.geomspace(1, tmax, n_t)
    lg = np.geomspace(1e-6, 1.0, n_l)
    l0, dl = math.log(lg[0]), math.log(lg[1] / lg[0])
    Mtab = np.minimum(2.0 * lg, 1.0)
    kanten = np.concatenate([[1.0], np.sqrt(tg[1:] * tg[:-1]), [tmax]])
    gew, rest = kanten[:-1]**(-p) - kanten[1:]**(-p), tmax**(-p)

    def bahn(M_tafel, dicht=False):
        lM = np.log(np.maximum(M_tafel, 1e-300)).tolist()
        n = len(lM) - 1

        def rechts(t, y):
            L, V = y
            lam = L * t**(-s)
            if lam >= 1:
                Mv = math.exp(lM[-1])
            elif lam <= lg[0]:
                Mv = math.exp(lM[0]) * lam / lg[0]
            else:
                x = (math.log(lam) - l0) / dl
                i = min(int(x), n - 1)
                Mv = math.exp(lM[i] + (x - i) * (lM[i + 1] - lM[i]))
            return [V, -K8 * t**p * Mv / L**2 + K8 * eps_j / L**3]
        return solve_ivp(rechts, (1, tmax), [1.0, 0.0], t_eval=None if dicht else tg, dense_output=dicht, rtol=rtol,
                         atol=rtol * 1e-3, method="DOP853")
    for _ in range(n_max):
        lamp = bahn(Mtab).y[0] * tg**(-s)
        o = np.argsort(lamp)
        cm = rest + np.cumsum(gew[o])
        Mneu = np.interp(lg, lamp[o], cm, left=rest, right=cm[-1])
        m = lg > 0.05
        d = np.max(np.abs(Mneu - Mtab)[m] / Mneu[m])
        Mtab = misch * Mneu + (1 - misch) * Mtab
        if d < grenze:
            break
    sol = bahn(Mtab, dicht=True)
    L = lambda t: sol.sol(t)[0]
    tt = np.geomspace(1, tmax, 20000)
    Y = sol.sol(tt)
    ip = int(np.argmax((Y[1][:-1] < 0) & (Y[1][1:] >= 0)))
    tp = brentq(lambda t: sol.sol(t)[1], tt[ip], tt[ip + 1], xtol=1e-12)
    j = int(np.argmax(np.where(tt > tp, Y[0] * tt**(-s), -1.0)))
    r = minimize_scalar(lambda t: -L(t) * t**(-s), bounds=(tt[max(j - 1, 0)], tt[min(j + 1, len(tt) - 1)]),
                        method="bounded", options={"xatol": 1e-12})
    lam_c = -r.fun
    tF = brentq(lambda t: L(t) * t**(-s) - lam_c, 1.0 + 1e-12, tp, xtol=1e-13)
    I_F = quad(lambda t: math.sqrt(K8 / L(t)**3), 1.0, tF, limit=200, epsrel=1e-10)[0]
    _UEBERGABE[schluessel] = (I_F, tF, float(L(tF)))
    return _UEBERGABE[schluessel]


# ------------------------------------------------------------------ Linse jenseits der Naht, Rest im Haufen
@formel('linse_punktmasse', 'ΔΣ(R) = M/(πR²) (Linsensignal einer Punktmasse im Abstand R)',
        ('ΔΣ',), ('M_ges', 'R'), 'L', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Milchstraße bei 1 Mpc ΔΣ = 0,12 M☉/pc²'))
def linse_punktmasse(M, R):
    """Σ̄(<R) − Σ(R) einer Punktmasse: jenseits 1,05·R_ta liest die Linse das Maß des Ganzen M_ges = (1 + w·Q)·M_b mit
    w ≈ 1, unabhängig von N*."""
    return M / (math.pi * R**2)


@formel('haufen_baryonen', 'M_b(<r) = M_g·(r′ − r_c·atan(r′/r_c))/(r_n − r_c·atan(r_n/r_c)), r′ = min(r, r_max) '
        '(β-Gas, β = ⅔) + Hernquist (Zentralgalaxie) + Hernquist (Galaxien)',
        ('M_b(<r)',), ('r', 'M', 'b_H'), 'L', 'hilfe', 'Baustein 03', None)
def haufen_baryonen(r, M_gas, r_c, r_norm, r_max, M_bcg, b_bcg, M_gal, b_gal):
    """Baryonen eines Haufens: Gas mit ρ ∝ 1/(1 + (r/r_c)²) bis r_max, M_gas innerhalb r_norm; dazu Zentralgalaxie und
    Galaxien als Hernquist-Kugeln."""
    from .galaxie import hernquist_masse
    prof = lambda x: x - r_c * math.atan(x / r_c)
    return (M_gas * prof(min(r, r_max)) / prof(r_norm) + hernquist_masse(M_bcg, b_bcg, r)
            + hernquist_masse(M_gal, b_gal, r))


@formel('masse_linse_nfw', 'M_L(<r) = M₂₀₀·μ(c·r/r₂₀₀)/μ(c), μ(x) = ln(1 + x) − x/(1 + x), r₂₀₀ aus 200·ρ_krit(z)',
        ('M_L(<r)',), ('r', 'M', 'c_NFW', 'z'), 'L', 'hilfe', 'Baustein 03', None)
def masse_linse_nfw(r, M200, c_nfw, z):
    """Gelesene Linsenmasse nach dem Referenzprofil (NFW) einer Linse bei z; ρ_krit(z) = 3H(z)²/(8πG) des
    Vergleichshintergrunds. Ein Maß der gelesenen Faltung, kein Stoff."""
    rho_k = 3 * (H0_si() * E_hintergrund(z))**2 / (8 * math.pi * G)
    r200 = (3 * M200 / (4 * math.pi * 200 * rho_k))**(1 / 3)
    mu = lambda x: math.log1p(x) - x / (1 + x)
    return M200 * mu(c_nfw * r / r200) / mu(c_nfw)


@formel('haufen_rest', 'F_H = M_L(<r)/(ν(GM_b/r²)·M_b(<r)) (gelesene Linse gegen die örtliche Regel allein)',
        ('F_H',), ('M_L(<r)', 'M_b(<r)', 'r', 'a₀'), 'L', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Mit der Regel allein bleibt in Haufen ein Rest'))
def haufen_rest(r, M_b, M_L, a0=A0):
    """Faktor, um den die gelesene Linsenmasse im Haufen die örtliche Regel auf den eigenen Baryonen übersteigt; die
    übergebene Mode trägt dort keinen Rest (gebunden w = 0)."""
    return M_L / (float(nu(G * M_b / r**2, a0)) * M_b)

