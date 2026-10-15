# -*- coding: utf-8 -*-
"""
Rückkehr (Ebene I, innere Richtung orthogonal zum Raum): die eine Rückkehrgleichung
∂_t u + (Γ + θ)u = D∂_ϖ²u, am Ort projiziert û = e^{−√((p + Γ + θ)·τ_D)}.
Daraus Rückkehranteil q, s, ν, χ, n_B, Q̃, Rückkehrzeit τ_D, Schließlänge Λ, ℓ₀, Gewicht w, Verlust Γ.
"""
from .konstanten import *
from .formeln import formel


@formel('verlust_gamma', 'Γ = a_N/c',
        ('Γ',), ('a_N',), 'LI', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'Γ = c|∇σ| = a_N/c'))
def verlust_gamma(aN):
    return aN / c


@formel('tau_D', 'τ_D = c/a₀',
        ('τ_D',), ('a₀',), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'τ_D = ϖ_R²/D = c/a₀'))
def tau_D(a0=A0):
    return c / a0


@formel('rueckkehr_exponent', 's_p = √((p + Γ + θ)·τ_D)',
        ('s_p',), ('p', 'Γ', 'θ', 'a₀'), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'û(ϖ_R) = e^{−√((p + Γ + θ)·τ_D)}'))
def rueckkehr_exponent(p, Gamma, theta=0.0, a0=A0):
    """Die eine Rückkehrgleichung ∂_t u + (Γ + θ)u = D∂_ϖ²u (07), am Ort projiziert:
    û(ϖ_R) = e^{−√((p + Γ + θ)·τ_D)}. Gibt den Exponenten √((p + Γ + θ)·τ_D) (p komplex erlaubt).
    p = θ = 0: Regel (q = e^{−s}); Γ = 0: Rest 1 − w des Gewichts; p ≠ 0: Kern K̂(p)."""
    x = (p + Gamma + theta) * tau_D(a0)
    return np.sqrt(x)


@formel('rueckkehr', 'û = e^{−√((p + Γ + θ)·τ_D)}',
        ('û',), ('p', 'Γ', 'θ', 'a₀'), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'û(ϖ_R) = e^{−√((p + Γ + θ)·τ_D)}'))
def rueckkehr(p, Gamma, theta=0.0, a0=A0):
    """û = e^{−√((p + Γ + θ)·τ_D)}: Rückkehranteil q, Rest 1 − w und Kern K̂ in einem Ausdruck."""
    return np.exp(-rueckkehr_exponent(p, Gamma, theta, a0))


@formel('s_', 's = √(Γ·τ_D) = √(a_N/a₀)',
        ('s',), ('a_N', 'a₀'), 'LI', 'modell', 'Baustein 07', ('07-feldgleichung.md', 's² = Γ·τ_D = a_N/a₀'))
def s_(aN, a0=A0):
    """s = √(Γ·τ_D) mit Γ = a_N/c: Exponent der einen Rückkehr ohne Ausdehnung (rueckkehr_exponent)."""
    return rueckkehr_exponent(0.0, np.asarray(aN) / c, 0.0, a0)


@formel('q_', 'q = e^{−√(a_N/a₀)}',
        ('q',), ('a_N', 'a₀'), 'LI', 'modell', 'Autor (FTM-PDF)', ('03-von-der-galaxie-bis-zum-kosmos.md', 'q = e^{−√(a_N/a₀)} = e^{−r*/r}'))
def q_(aN, a0=A0):
    return np.exp(-s_(aN, a0))


@formel('nu', 'ν = 1/(1 − e^{−s}) = 1 + n_B(s)',
        ('ν',), ('a_N', 'a₀'), 'LI', 'modell', 'Autor (FTM-PDF)', ('07-feldgleichung.md', 'also ν = 1 + n_B(s)'))
def nu(aN, a0=A0):
    """Verstärkung ν = 1/(1 − e^{−s}) = 1 + n_B(s)."""
    return 1 / (-np.expm1(-s_(aN, a0)))


@formel('chi', 'χ = s/(e^s − 1)',
        ('χ',), ('a_N', 'a₀'), 'LI', 'modell', 'Autor (FTM-PDF)', ('07-feldgleichung.md', 'χ = s/(e^s − 1),     s = √(a_N/a₀)'))
def chi(aN, a0=A0):
    s = s_(aN, a0)
    return s / np.expm1(s)


@formel('n_B', 'n_B = 1/(e^s − 1)',
        ('n_B',), ('s',), 'LI', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'n_B = 1/(e^s − 1)'))
def n_B(s):
    return 1 / np.expm1(s)


@formel('Q_tilde', 'Q̃(s) = −4∫_s^∞ x³n_B(x)dx',
        ('Q̃',), ('s',), 'LG', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'Q̃ = −4∫_s^∞ x³n_B(x)dx'))
def Q_tilde(s):
    """Strukturwirkung Q̃(s) = −4∫_s^∞ x³ n_B(x) dx = −4[s³Li₁ + 3s²Li₂ + 6s·Li₃ + 6Li₄](e^{−s}) (07); Q̃(0) = −4π⁴/15."""
    from scipy.integrate import quad
    return -4 * quad(lambda x: x**3 * math.exp(-x) / -math.expm1(-x) if x > 0 else 0.0, s, np.inf)[0]


@formel('bose_stufe_max', 'd(1 − e^{−s}) = s: Maximum von s^d·n_B; gibt (s_d, a_N/a₀, Wert)',
        ('s_d',), ('d',), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'd(1 − e^{−s}) = s'))
def bose_stufe_max(d):
    """Maximum von s^d·n_B(s): d(e^s − 1) = s·e^s; gibt (s, a_N/a₀, Wert)."""
    from scipy.optimize import brentq
    s = brentq(lambda x: d * np.expm1(x) - x * math.exp(x), 0.1, 20)
    return s, s * s, s**d * n_B(s)


@formel('schliesslaenge', 'Λ = c²/(2a₀)',
        ('Λ',), ('a₀',), 'I', 'modell', 'Autor (FTM-PDF)', ('07-feldgleichung.md', 'ist a₀ = c²/(2Λ)'))
def schliesslaenge(a0=A0):
    """Λ = c²/(2a₀)."""
    return c**2 / (2 * a0)


@formel('ell0', 'ℓ₀ = 2Λ = c²/a₀',
        ('ℓ₀',), ('a₀',), 'I', 'modell', 'Baustein 08', ('08-licht-und-trichter.md', 'ℓ₀ = c²/a₀'))
def ell0(a0=A0):
    """ℓ₀ = 2Λ = c²/a₀; Leiter r*² = 3Δs·ℓ₀."""
    return 2 * schliesslaenge(a0)


@formel('a0_aus_schliesslaenge', 'a₀ = c²/(2Λ)',
        ('a₀',), ('Λ',), 'I', 'modell', 'Autor (FTM-PDF)', ('07-feldgleichung.md', 'ist a₀ = c²/(2Λ)'))
def a0_aus_schliesslaenge(Lam):
    return c**2 / (2 * Lam)


@formel('gewicht_w', 'w = 1 − e^{−√(θ̃/θ_c)}, θ_c = a₀/c',
        ('w',), ('θ̃', 'a₀'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'w = 1 − e^{−√(θ̃/θ_c)}'))
def gewicht_w(theta, a0=A0):
    """Schließungsregel aus derselben Rückkehr wie die Regel: 1 − w = E[e^{−θ₊T}] über die Dichte
    erster Rückkehrzeiten, also w = 1 − e^{−√(θ₊/θ_c)}, θ_c = a₀/c = 1/τ_D."""
    th = max(theta, 0.0)
    return -math.expm1(-float(rueckkehr_exponent(0.0, 0.0, th, a0)))


@formel('gewicht_w_preisach', 'w(θ̃), θ̃ = max(min θ der Geschichte, 0)',
        ('w',), ('θ', 'a₀'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'θ̃ = max(min_{t′≤t} θ(t′), 0)'), alternative='gewicht_w')
def gewicht_w_preisach(theta_geschichte, a0=A0):
    """Schließungsregel mit Gedächtnis: w(θ̃), θ̃ = max(min θ der Geschichte, 0)."""
    return gewicht_w(min(theta_geschichte), a0)


# ------------------------------------------------------------------ Kern K in Zeit und Frequenz (07)
@formel('kern_zeit', 'K(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/4T − ΓT}',
        ('K(T)',), ('T_R', 'Γ', 'a₀'), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'K(T) = √(τ_D/4π)·T^{−3/2}·e^{−τ_D/4T − ΓT}'))
def kern_zeit(T, Gamma, a0=A0):
    """Dichte erster Rückkehrzeiten aus der Tiefe ϖ_R (Irrfahrt mit einem Takt), gewichtet mit dem Verlust e^{−ΓT};
    T ist die Rückkehrzeit in der inneren Richtung. Ihr Laplace-Bild ist rueckkehr(p, Γ)."""
    tD = tau_D(a0)
    return math.sqrt(tD / (4 * math.pi)) * T**-1.5 * math.exp(-tD / (4 * T) - Gamma * T)


@formel('kern_bahnfrequenz', '|K̂(iΩ)|/K̂(0) = e^{−Re√((Γ + iΩ)τ_D) + √(Γτ_D)}',
        ('|K̂(iΩ)|/K̂(0)',), ('Ω', 'Γ', 'a₀'), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', '|K̂(iΩ)|/K̂(0) ≈ e^{−√(Ωτ_D/2) + s}'))
def kern_bahnfrequenz(Omega, Gamma, a0=A0):
    """Was die Bahn vom Kern sieht: K̂ bei der Bahnfrequenz p = iΩ gegen den stationären Wert K̂(0) = q.
    Für Ωτ_D ≫ s² gilt e^{−√(Ωτ_D/2) + s}; die Bahn sieht nur die stationäre Regel."""
    return abs(complex(rueckkehr(1j * Omega, Gamma, 0.0, a0))) / float(rueckkehr(0.0, Gamma, 0.0, a0))


# ------------------------------------------------------------------ Halbgerade, Glättung, Erneuerung (07, 01)
@formel('glaettung_innen', 'P_ϖ = e^{−|k|ϖ}, P_a·P_b = P_{a+b}',
        ('P_ϖ',), ('k_ϖ', 'ϖ'), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'setzt das Feld harmonisch fort, e^{−|k|ϖ}'), triade='Wirbel')
def glaettung_innen(k, varpi):
    """Harmonische Fortsetzung einer Mode k in die innere Richtung, (∇² + ∂_ϖ²)δσ = 0: Faktor e^{−|k|ϖ}. Beim Aufsteigen
    klingt nur die eigene aufgelöste Struktur (Wirbel) ab; der Fluss GM (k = 0) bleibt. ϖ ≥ 0 ist eine Halbgerade."""
    if varpi < 0:
        raise ValueError("ϖ ist eine Halbgerade, ϖ ≥ 0")
    return math.exp(-abs(k) * varpi)


@formel('erneuerung', 'u(ϖ) = u(0)·e^{−ϖ√(Γ/D)}, u(0) = g_N + u(ϖ_R) = g_N/(1 − q)',
        ('u',), ('g_N', 'ϖ', 'M', 'a₀'), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'u(x, 0) = g_N + u(x, ϖ_R)'))
def erneuerung(gN, varpi, M, a0=A0):
    """Stationäre Rückkehr am Ort einer Punktmasse: ∂_ϖ²u = (Γ/D)·u, beschränkt, mit der Erneuerung
    u(0) = g_N + u(ϖ_R); Γ/D = a_N/GM (a_N = |g_N|), ϖ_R = r*. Am Ort (ϖ = 0) ist u = a = ν·g_N."""
    k = math.sqrt(abs(gN) / (G * M))
    return gN / (1 - float(q_(abs(gN), a0))) * math.exp(-k * varpi)


# ------------------------------------------------------------------ Schluss bei Λ (07)
@formel('eigenwert_schluss', 'λ_Λ = Dπ²/(4Λ²)',
        ('λ_Λ',), ('D', 'Λ'), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'kleinsten Eigenwert Dπ²/(4Λ²)'))
def eigenwert_schluss(D, Lam):
    """Kleinster Eigenwert der bei Λ geschlossenen Rückkehr (Ort ϖ = 0 fest, bei Λ reflektierend)."""
    return D * math.pi**2 / (4 * Lam**2)


@formel('rueckkehr_geschlossen', 'q_Λ = cosh(s(Λ/ϖ_R − 1))/cosh(s·Λ/ϖ_R)',
        ('q_Λ',), ('s', 'Λ/ϖ_R'), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'Die bei Λ geschlossene Rückkehr ist exakt der thermische Propagator der Periode 2Λ/c'))
def rueckkehr_geschlossen(s, Lam_zu_varpiR):
    """Rückkehranteil der bei Λ geschlossenen Halbgeraden: G'' = k²G, G(0) = 1, G'(Λ) = 0, gelesen bei ϖ_R (kϖ_R = s).
    Gleich dem thermischen Propagator der Periode 2Λ; Λ/ϖ_R → ∞ gibt q = e^{−s}."""
    L = Lam_zu_varpiR
    if s * L > 700:
        return math.exp(-s)
    return math.cosh(s * (L - 1)) / math.cosh(s * L)


@formel('abweichung_schluss', 'ε_Λ = 1/tanh(kΛ) − 1',
        ('ε_Λ',), ('kΛ',), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'dort ist die Abweichung 1/tanh(π/2) − 1 = 9,0 %'))
def abweichung_schluss(kLam):
    """Abweichung der Regel von der offenen Rückkehr im tiefen Bereich (s → 0 bei festem kΛ): ν_Λ/ν − 1."""
    return 1 / math.tanh(kLam) - 1


@formel('radius_schluss', 'r_Λ = (2/π)·Λ·r*/ϖ_R',
        ('r_Λ',), ('Λ', 'r*', 'ϖ_R'), 'LI', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'r = (2/π)·Λ·r*/ϖ_R'))
def radius_schluss(Lam, rstern, varpiR):
    """Ort, an dem die Verlustrate Γ auf den kleinsten Eigenwert Dπ²/(4Λ²) fällt (Punktmasse Γ/D = 1/r², ϖ_R-Skala)."""
    return 2 / math.pi * Lam * rstern / varpiR


@formel('temperatur_innen', 'kT_I = ħa₀/c = ħ/(Periode 2Λ/c)',
        ('T_I',), ('a₀',), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'kT_I = ħa₀/c'))
def temperatur_innen(a0=A0):
    """Die Periode 2Λ/c der geschlossenen inneren Richtung, als Temperatur geschrieben; eine Periode, keine Wärme."""
    return hbar * a0 / (c * kB)


@formel('kac_zeit_schluss', '⟨T⟩ = 2Λ/c = τ_D',
        ('τ_D',), ('Λ',), 'I', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', '⟨T⟩ = 2Λ/c, unabhängig von der Schrittweite'), alternative='tau_D')
def kac_zeit_schluss(Lam):
    """Kac: mittlere Rückkehrzeit zum Ort auf der bei Λ geschlossenen Halbgeraden mit c-Schritten, für jede Masse."""
    return 2 * Lam / c


# ------------------------------------------------------------------ Strukturteil des Energiefunktionals (07)
@formel('strukturteil', 'Q_s(x) = ∫₀^x (ν(√x′) − 1) dx′, x = s⁴',
        ('Q_s',), ('x_V',), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'Q_s(x) = ∫₀^x (ν(√x′) − 1) dx′'))
def strukturteil(x):
    """Strukturteil der Hierarchie-Wirkung; ν(√x′) − 1 = n_B(x′^{1/4}). Q_s(0) = 0 je Zone (kein örtlicher Vakuumterm)."""
    from scipy.integrate import quad
    if x <= 0:
        return 0.0
    return quad(lambda y: float(n_B(y**0.25)) if y > 0 else 0.0, 0.0, x, limit=200)[0]


# ------------------------------------------------------------------ Form der Rückkehr: Gitter, Drehung (07)
@formel('watson_rueckkehr', 'G_d(1) = ∫₀^∞ e^{−t}·I₀(t/d)^d dt = 1/(1 − p_d)',
        ('G_d(1)',), ('d',), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'in drei ist ν ≤ 1,516 (Watsons Integral)'))
def watson_rueckkehr(d):
    """Summe aller Rückläufe der Irrfahrt auf dem Gitter ℤ^d ohne Verlust (q → 1): endlich erst ab d = 3 (Pólya);
    zugleich die größte Verstärkung, die eine Rückkehr in d Dimensionen geben kann."""
    from scipy.integrate import quad
    from scipy.special import ive
    d = int(d)
    if not 3 <= d <= 50:
        raise ValueError("3 ≤ d ≤ 50 (d ≤ 2 divergent)")
    f = lambda t: float(ive(0, t / d))**d
    return quad(f, 0, 50, limit=400, epsabs=0, epsrel=1e-12)[0] + quad(f, 50, np.inf, limit=400, epsabs=0,
                                                                         epsrel=1e-12)[0]


@formel('rueckkehrzeiten_drehung', 'T_μ = {n ≥ 1: erste Rückkehr von x + nα_D (mod 1) in [0, μ)}',
        ('T_μ',), ('α_D', 'μ(A)'), 'I', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'golden, Zelle 1 %: 55, 89, 144'))
def rueckkehrzeiten_drehung(alpha, mu):
    """Alle Werte der ersten Rückkehrzeit einer Drehung um α_D in den Bogen [0, μ), aus dem Abstand aufeinanderfolgender
    Besuche einer Bahn (höchstens drei Werte, der größte die Summe der beiden anderen)."""
    if not (0 < mu < 0.5 and 0 < alpha < 1):
        raise ValueError("0 < μ < ½, 0 < α < 1")
    n = np.arange(int(400 / mu))
    besuche = np.nonzero(np.mod(n * alpha, 1.0) < mu)[0]
    return np.unique(np.diff(besuche)).astype(float)


# ------------------------------------------------------------------ Windungskorrelation am Knoten (04)
@formel('chsh_oberschwingungen', 'S = max_θ |3E(θ) − E(3θ)|, E(Δ) = −Σ_{m<N_H} w_m cos((2m+1)Δ), w_m ∝ (2m+1)⁻²',
        ('S_CHSH',), ('N_H',), 'I', 'modell', 'Baustein 04', ('04-zeit-verschraenkung-syntropie-gedaechtnis.md',
                                                             'zwei Harmonische des Dreiecks: 2,36'))
def chsh_oberschwingungen(N_H):
    """Korrelation der Windung mit den ersten N_H Harmonischen des Dreiecks (auf E(0) = −1 normiert), Winkel
    a = 0, b = θ, a′ = 2θ, b′ = 3θ der drehinvarianten Korrelation. N_H = 1: Grundschwingung, 2√2; N_H → ∞: Dreieck
    (Eigenschaft, die die Quelle mitgibt), 2."""
    from scipy.optimize import minimize_scalar
    N_H = int(N_H)
    if not 1 <= N_H <= 10000:
        raise ValueError("1 ≤ N_H ≤ 10⁴")
    k = 2 * np.arange(N_H) + 1.0
    w = 1 / k**2
    w /= w.sum()
    S = lambda t: -abs(3 * np.dot(w, np.cos(k * t)) - np.dot(w, np.cos(3 * k * t)))
    th = np.linspace(0, math.pi / 2, 4001)
    t0 = th[int(np.argmin([S(t) for t in th]))]
    r = minimize_scalar(S, bounds=(max(t0 - 1e-3, 0), min(t0 + 1e-3, math.pi / 2)), method='bounded',
                        options=dict(xatol=1e-12))
    return -float(r.fun)
