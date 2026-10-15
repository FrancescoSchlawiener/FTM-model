# -*- coding: utf-8 -*-
"""
Galaxie (Ebene L mit der Regel aus I): örtliche Regel a = ν(a_N)·a_N, Kreisbahn, flache Rotation,
äußeres Feld, Grenzflächendichte, Bahnintegrator im achsensymmetrischen Newton-Feld, Hernquist-Masse.
Lesung an Rotationskurven (SPARC): Rest je Galaxie mit dem Druck des Gases und äußerem Feld, Schwelle je Galaxie,
Streuung aus den Messfehlern (Monte Carlo, eigener Generator je Aufruf), Toomre-Integral, Bose-Form an Bin-Medianen,
Spektraldimension. Scheiben: Exponential-, Kuzmin- und Miyamoto-Nagai-Feld, Regel in der Ebene mit dem senkrechten
Scheibenfeld, Vertikalfrequenzen, Präzession geneigter Ringe. System als Ganzes im äußeren Feld: Massenmittel A,
Einseitigkeit m = 1, Plummer-Kugel mit Jeans-Gleichung. Gezeit und Hierarchie der Trichter: 𝒯 aus der spurfreien
Hesse-Form am Sonnenort, 𝒯_H des Hill-Radius (eine andere Größe); 𝒯, ℓ, ρ_J, M×, R× und die Flussdimension
stehen in faltung.
"""
from .konstanten import *
from .formeln import formel
from .rueckkehr import nu, chi
from .faltung import ds, gezeit_betrag
from .daten import SPARC_SATZ

_F03 = "03-von-der-galaxie-bis-zum-kosmos.md"
_F07 = "07-feldgleichung.md"
_UPS_B = SPARC_SATZ["ups_b"]


@formel('a_regel', 'a = ν(a_N)·a_N = a_N/(1 − q)',
        ('a',), ('a_N', 'a₀'), 'LI', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'a = a_N · Σ qⁿ = a_N/(1 − q)'))
def a_regel(aN, a0=A0):
    """Örtliche Regel, Betrag: a = ν(a_N)·a_N."""
    return nu(aN, a0) * aN


@formel('v_kreis_punkt', 'v = √(ν(a_N)·a_N·r), a_N = GM/r²',
        ('v',), ('M', 'r', 'a₀'), 'L', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Milchstraße als Punktmasse: v = 250 km/s bei 5 kpc'))
def v_kreis_punkt(M, r, a0=A0):
    aN = G * M / r**2
    return math.sqrt(a_regel(aN, a0) * r)


@formel('v_flach', 'v_flach⁴ = GMa₀',
        ('v_flach',), ('M', 'a₀'), 'L', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'v⁴ = GMa₀  (flache Rotationskurve, Tully-Fisher)'))
def v_flach(M, a0=A0):
    return (G * M * a0) ** 0.25


@formel('nu_aussenfeld', 'ν(|g|), senkrecht s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e²',
        ('ν',), ('a_i', 'a_e', 'a₀'), 'L', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Steht das äußere Feld senkrecht, gilt s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e²'), alternative='nu')
def nu_aussenfeld(a_i, a_e, lage, a0=A0):
    """ν bei äußerem Feld: senkrecht |g| = √(a_i² + a_e²) (s⁴ = s_i⁴ + s_e⁴), gleichgerichtet |g| = a_i + a_e."""
    g = math.hypot(a_i, a_e) if lage == "senkrecht" else a_i + a_e
    return float(nu(g, a0))


@formel('sigma_kollektiv', 'Σ_a₀ = a₀/(πG), Rückgabe in M☉/pc²',
        ('Σ_a₀',), ('a₀',), 'L', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Σ < a₀/(πG) ≈ 251 M☉/pc²'))
def sigma_kollektiv(a0=A0):
    """Grenzflächendichte a₀/(πG) in M☉/pc²."""
    return a0 / (math.pi * G) * pc**2 / M_sun


@formel('hernquist_masse', 'M(<r) = M·r²/(r + b_H)²',
        ('M(<r)',), ('M', 'b_H', 'r'), 'L', 'hilfe', 'Baustein 08', None)
def hernquist_masse(M, a, r):
    return M * r * r / (r + a) ** 2




# ================================================================== Lesung an Rotationskurven (SPARC)
@formel('v_bar_quadrat', 'V_bar² = V_gas|V_gas| + Υ_d·V_disk² + Υ_b·V_bul²',
        ('V_bar',), ('V_gas', 'V_disk', 'V_bul', 'Υ_d', 'Υ_b'), 'L', 'modell', 'Baustein 03',
        (_F03, 'mit den eigenen Baryonen'))
def v_bar_quadrat(Vg, Vd, Vb, ups_d, ups_b=_UPS_B):
    """Quadrat der baryonischen Rotation mit Vorzeichen (V_gas < 0: Gas zieht nach außen); g_N = V_bar²/R."""
    return Vg * np.abs(Vg) + ups_d * Vd**2 + ups_b * Vb**2


@formel('g_gemessen', 'g_obs = (V_obs² + σ_v²·R/R_d)/R',
        ('g_obs',), ('V_obs', 'σ_v', 'r', 'R_d'), 'L', 'modell', 'Baustein 03',
        (_F03, 'bis man den Druck des Gases einrechnet'))
def g_gemessen(R, Vo, Rd, sig):
    """Gemessene Zentripetalbeschleunigung mit dem Druck des Gases (Streuung σ_v, Gas ∝ e^{−R/R_d}):
    V_c² = V_obs² + σ_v²·R/R_d."""
    return (Vo**2 + sig**2 * R / Rd) / R


@formel('rest_galaxie', 'ε_R = Median log₁₀(g_obs/(ν(|g_N + g_e|)·g_N)), g_N = V_bar²/R, g_e senkrecht',
        ('ε_R',), ('V_obs', 'V_gas', 'V_disk', 'V_bul', 'r', 'R_d', 'Υ_d', 'Υ_b', 'σ_v', 'a_e', 'a₀'), 'L', 'modell',
        'Baustein 03', (_F03, 'hängt der Rest je Galaxie an keiner der geprüften Größen'))
def rest_galaxie(R, Vo, Vg, Vd, Vb, Rd, ups_d, sig=0.0, ge=0.0, a0=A0, ups_b=_UPS_B, teil="alle"):
    """Rest einer Galaxie gegen die örtliche Regel mit ihren eigenen Baryonen (Punkte mit V_bar² > 0 und V_obs > 0,
    mindestens 5; sonst nan). ge: gemeinsames äußeres Feld senkrecht zur Scheibe (m/s²). teil: "alle", "innen"
    (R < 2R_d) oder "aussen" (R > 2R_d)."""
    Vb2 = v_bar_quadrat(Vg, Vd, Vb, ups_d, ups_b)
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        return float("nan")
    gN = Vb2[m] / R[m]
    n_ = nu(gN, a0) if ge == 0 else np.array([nu_aussenfeld(x, ge, "senkrecht", a0) for x in gN])
    r = np.log10(g_gemessen(R[m], Vo[m], Rd, sig) / (n_ * gN))
    if teil != "alle":
        r = r[R[m] < 2 * Rd] if teil == "innen" else r[R[m] > 2 * Rd]
    return float(np.median(r)) if r.size else float("nan")


@formel('a0_je_galaxie', 'a₀_gal: ε_R(a₀) = 0 je Galaxie',
        ('a₀_gal',), ('V_obs', 'V_gas', 'V_disk', 'V_bul', 'r', 'R_d', 'Υ_d', 'Υ_b', 'σ_v'), 'LI', 'modell',
        'Baustein 03', (_F03, 'Je Galaxie gelesen (mit dem Druck des Gases)'))
def a0_je_galaxie(R, Vo, Vg, Vd, Vb, Rd, ups_d, sig=0.0, ups_b=_UPS_B):
    """Schwelle, bei der der Rest der Galaxie null ist (Nullstelle in log a₀ zwischen 10⁻¹³ und 10⁻⁷ m/s²; sonst nan)."""
    from scipy.optimize import brentq
    f = lambda la: rest_galaxie(R, Vo, Vg, Vd, Vb, Rd, ups_d, sig, 0.0, 10**la, ups_b)  # noqa: E731
    try:
        return 10 ** brentq(f, -13, -7)
    except ValueError:
        return float("nan")


@formel('streuung_messfehler', 'σ_ε = Streuung einer Lesung über Ziehungen von Abstand, Neigung und Υ',
        ('σ_ε',), ('D_A', 'i_inc', 'Υ_d'), 'L', 'modell', 'Baustein 03',
        (_F03, 'aus den Messfehlern von Abstand, Neigung und Masse-zu-Licht'))
def streuung_messfehler(gal, lesung, n_mc, seed, ups_streuung=SPARC_SATZ["ups_streuung_dex"]):
    """Monte Carlo über die Fehler einer Galaxie (Tafelzeile mit D, eD, inc, eInc und Kurve R, Vo, Vg, Vd, Vb, Rd):
    Abstand D(1 + δ), δ ~ N(0, e_D/D), Faktor ≥ 0,2; Neigung i + N(0, e_i) in [20°, 90°]; Υ·10^{N(0, ups_streuung)}.
    Skalierung: R, R_d ∝ D; V_obs ∝ sin i/sin i′; V_gas² ∝ D; V_disk², V_bul² ∝ Υ·D. lesung(R, Vo, Vg, Vd, Vb, Rd)
    gibt eine Zahl je Ziehung; Ergebnis ist deren Streuung (ohne nan). Eigener Generator je Aufruf mit seed."""
    rng = np.random.default_rng(seed)
    werte = []
    for _ in range(int(n_mc)):
        fD = max(1 + rng.normal() * gal["eD"] / gal["D"], 0.2)
        inc = min(max(gal["inc"] + rng.normal() * gal["eInc"], 20.0), 90.0)
        fi = math.sin(math.radians(gal["inc"])) / math.sin(math.radians(inc))
        fU = 10 ** (rng.normal() * ups_streuung)
        werte.append(lesung(gal["R"] * fD, gal["Vo"] * fi, gal["Vg"] * math.sqrt(fD), gal["Vd"] * math.sqrt(fU * fD),
                            gal["Vb"] * math.sqrt(fU * fD), gal["Rd"] * fD))
    return float(np.nanstd(werte))


@formel('toomre_sigma_dyn', 'Σ_dyn(0) = (1/2πG)∫₀^∞ V²/R² dR, innen V ∝ R, außen V = V_letzter',
        ('Σ_dyn',), ('v', 'r'), 'L', 'modell', 'Baustein 03',
        (_F03, 'Zentrale Flächendichte aus der Rotationskurve (Toomre-Integral)'))
def toomre_sigma_dyn(R, V):
    """Toomre-Integral einer Kurve (Arrays, R aufsteigend): innen bis zum ersten Punkt V ∝ R (Beitrag V₀²/R₀), zwischen
    den Punkten Trapez, außen V konstant (Beitrag V_n²/R_n). In kg/m²."""
    return (V[0]**2 / R[0] + np.trapezoid((V / R)**2, R) + V[-1]**2 / R[-1]) / (2 * math.pi * G)


# ================================================================== Scheiben: Felder und Regel in der Ebene
@formel('g_expscheibe', 'g_R = GM·R/(2R_d³)·[I₀K₀ − I₁K₁](R/2R_d)',
        ('g_R',), ('M', 'R_d', 'r'), 'L', 'modell', 'Baustein 03', (_F03, 'Milchstraße als Exponentialscheibe'))
def g_expscheibe(M, Rd, R):
    """Newton-Feld einer dünnen Exponentialscheibe in ihrer Ebene (Freeman), Betrag; skalierte Bessel-Funktionen."""
    from scipy.special import i0e, i1e, k0e, k1e
    y = np.asarray(R, float) / (2 * Rd)
    return G * M * np.asarray(R, float) / (2 * Rd**3) * (i0e(y) * k0e(y) - i1e(y) * k1e(y))


@formel('sigma_dyn_expscheibe', 'Σ_dyn(0) = Σ₀·∫₀^∞ ν(2πGΣ₀·f(y))·f(y) dy/y, f(y) = y[I₀K₀ − I₁K₁](y), y = R/2R_d',
        ('Σ_dyn',), ('Σ₀', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'eine ideale Exponentialscheibe läge'),
        alternative='toomre_sigma_dyn')
def sigma_dyn_expscheibe(Sigma0, a0=A0):
    """Toomre-Integral der idealen dünnen Exponentialscheibe unter der örtlichen Regel in der Ebene (radiale Form);
    hängt nur an Σ₀ und a₀. Jenseits y = 10⁴ trägt f = 1/(4y²) (Punktmasse)."""
    from scipy.integrate import quad
    gs = 2 * math.pi * G * Sigma0
    f = lambda y: float(g_expscheibe(1.0, 1.0, 2 * y)) / G       # noqa: E731  f(y) = y[I₀K₀ − I₁K₁](y)
    innen = quad(lambda u: float(nu(gs * f(math.exp(u)), a0)) * f(math.exp(u)), math.log(1e-8), math.log(1e4),
                 limit=400)[0]
    aussen = quad(lambda y: float(nu(gs / (4 * y * y), a0)) / (4 * y**3), 1e4, np.inf, limit=200)[0]
    return Sigma0 * (innen + aussen)


@formel('kuzmin_feld', 'Kuzmin an der Scheibe (z = 0⁺): g_R = GMR/(R² + b²)^{3/2}, g_z = GMb/(R² + b²)^{3/2} = 2πGΣ',
        ('g_R', 'g_z'), ('M', 'b_K', 'r'), 'L', 'modell', 'Baustein 07', (_F07, 'Für die Kuzmin-Scheibe gilt exakt'),
        alternative='g_expscheibe')
def kuzmin_feld(M, b, R, komponente):
    """Newton-Feld der Kuzmin-Scheibe knapp über der Ebene, Betrag der Komponente "R" oder "z". Oberhalb der Ebene ist es
    das Feld einer Punktmasse bei z = −b; darum gilt dort in der Feldform exakt g = ν(|g_N|)·g_N."""
    n = G * M / (R * R + b * b) ** 1.5
    return n * (R if komponente == "R" else b)


@formel('a_scheibe_ebene', 'a_R = ν(√(g_R² + f·(2πGΣ)²))·g_R',
        ('a_R',), ('g_R', 'Σ', 'f_Σ', 'a₀'), 'L', 'modell', 'Baustein 07', (_F07, 'a = ν(√(g_R² + f·(2πGΣ)²))·g_R'))
def a_scheibe_ebene(gR, Sigma, f, a0=A0):
    """Radiale Beschleunigung in der Ebene einer Scheibe: f = 0 örtliche Bahnregel in der Mittelebene (radiale Form),
    f = 1 Feldform (das Scheibenfeld wirkt wie ein senkrechtes Außenfeld)."""
    return nu(np.sqrt(gR**2 + f * (2 * math.pi * G * Sigma)**2), a0) * gR


@formel('g_miyamoto_nagai', 'g_R = GMR/D, g_z = GMz(a + ζ)/(ζD), ζ = √(z² + b²), D = (R² + (a + ζ)²)^{3/2}',
        ('g_R', 'g_z'), ('M', 'a_MN', 'b_MN', 'r', 'z_h'), 'L', 'modell', 'Baustein 07', (_F07, 'Miyamoto-Nagai-Scheibe'),
        alternative='g_expscheibe')
def g_miyamoto_nagai(M, a, b, R, z, komponente):
    """Newton-Feld der Miyamoto-Nagai-Scheibe Φ = −GM/√(R² + (a + √(z² + b²))²), Komponente "R" oder "z" (zur Mitte
    bzw. zur Ebene hin; g_z mit dem Vorzeichen von z)."""
    zb = math.sqrt(z * z + b * b)
    D = (R * R + (a + zb) ** 2) ** 1.5
    return G * M * R / D if komponente == "R" else G * M * z * (a + zb) / (zb * D)


@formel('vertikalfrequenz', 'Kernform ω_z² = ω_z,N² + (ν − 1)Ω_N², örtliche Regel ω_z² = ν·ω_z,N²',
        ('ω_z',), ('ω_z,N', 'Ω_N', 'ν'), 'L', 'modell', 'Baustein 07',
        (_F07, 'Verschieden sind auch die Vertikalfrequenzen'))
def vertikalfrequenz(omega_zN, Omega_N, nu_, form):
    """Vertikalfrequenz in der Mittelebene: form "kern" (radialer Kern e_r⊗e_r) oder "ort" (a = ν(|g_N|)·g_N)."""
    if form == "kern":
        return math.sqrt(omega_zN**2 + (nu_ - 1) * Omega_N**2)
    return math.sqrt(nu_) * omega_zN


@formel('knotenpraezession_klein', 'ω_p = Ω − ω_z (Neigung → 0), Frequenzen aus dem Feld am Ort',
        ('ω_p',), ('M', 'a_MN', 'b_MN', 'r', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'die √ν-fache der Newtonschen'), alternative='praezession_ring')
def knotenpraezession_klein(M, a, b, R, regel=True, a0=A0, h=None):
    """Knotenpräzession eines schwach geneigten Rings in der Miyamoto-Nagai-Scheibe (negativ = rückläufig):
    Ω² = a_R/R, ω_z² = a_z(R, h)/h mit a = ν(|g_N|)·g_N (regel) oder g_N."""
    h = 1e-4 * kpc if h is None else h

    def a_(Rr, z, k):
        gR, gz = g_miyamoto_nagai(M, a, b, Rr, z, "R"), g_miyamoto_nagai(M, a, b, Rr, z, "z")
        f = float(nu(math.hypot(gR, gz), a0)) if regel else 1.0
        return f * (gR if k == "R" else gz)
    return math.sqrt(a_(R, 0.0, "R") / R) - math.sqrt(a_(R, h, "z") / h)


@formel('praezession_ring', 'ω_p = mittlere Drehrate des Drehimpulses eines geneigten Kreisrings (3D-Bahn, a = ν(|g_N|)·g_N)',
        ('ω_p',), ('M', 'a_MN', 'b_MN', 'r', 'i_inc', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'Präzessionsraten geneigter Ringe'))
def praezession_ring(M, a, b, R, inc, regel=True, a0=A0, n_orb=12):
    """Knotenpräzession [s⁻¹] eines anfangs kreisförmigen Rings mit Neigung inc (rad) in der Miyamoto-Nagai-Scheibe:
    Bahn über n_orb Umläufe (DOP853, rtol 10⁻¹⁰), Steigung des Azimuts des Drehimpulses."""
    from scipy.integrate import solve_ivp

    def acc(x):
        Rr = math.hypot(x[0], x[1])
        gR, gz = g_miyamoto_nagai(M, a, b, Rr, x[2], "R"), g_miyamoto_nagai(M, a, b, Rr, x[2], "z")
        f = float(nu(math.hypot(gR, gz), a0)) if regel else 1.0
        return np.array([-f * gR * x[0] / Rr, -f * gR * x[1] / Rr, -f * gz])
    Om = math.sqrt(-acc(np.array([R, 0.0, 0.0]))[0] / R)
    v = Om * R
    y0 = np.array([R, 0.0, 0.0, 0.0, v * math.cos(inc), v * math.sin(inc)])
    T = 2 * math.pi / Om
    sol = solve_ivp(lambda t, y: np.concatenate([y[3:], acc(y[:3])]), (0, n_orb * T), y0, method="DOP853",
                    rtol=1e-10, atol=1e-6, dense_output=True)
    ts = np.linspace(0, n_orb * T, 40 * n_orb + 1)
    Y = sol.sol(ts)
    L = np.cross(Y[:3].T, Y[3:].T)
    return float(np.polyfit(ts, np.unwrap(np.arctan2(L[:, 1], L[:, 0])), 1)[0])


# ================================================================== System als Ganzes im äußeren Feld
@formel('einseitigkeit_m1', 'A_m1 = 2|F₁|/F₀ von v_c(φ), v_c² = −R·(ν(|g|)g − B)·e_R, g = −g_i·e_R + g_e in der Ebene',
        ('A_m1',), ('a_i', 'r', 'a_e', 'A', 'a₀'), 'L', 'modell', 'Baustein 03',
        (_F03, 'Die Rotation schwankt mit dem Azimut (m = 1)'))
def einseitigkeit_m1(g_i, R, g_e, bezug, a0=A0, n_phi=720):
    """m = 1-Anteil der Rotation auf dem Kreis R, wenn das äußere Feld g_e in der Scheibenebene liegt; g_i Betrag des
    inneren Felds bei R (zur Mitte), bezug: Beschleunigung des Systems als Ganzes längs g_e (g_e, ν(g_e)·g_e oder A)."""
    ph = np.arange(n_phi) * 2 * math.pi / n_phi
    eR = np.stack([np.cos(ph), np.sin(ph)])
    gt = -g_i * eR + np.array([g_e, 0.0])[:, None]
    acc = nu(np.linalg.norm(gt, axis=0), a0) * gt - np.array([bezug, 0.0])[:, None]
    vc = np.sqrt(np.clip(-R * np.sum(acc * eR, axis=0), 0, None))
    F = np.fft.rfft(vc) / n_phi
    return float(2 * abs(F[1]) / F[0].real)


@formel('massenmittel_A_plummer', 'A = ⟨ν(|g_i + g_e|)(g_i + g_e)⟩_M, Plummer-Kugel',
        ('A',), ('M', 'R_e', 'a_e', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'A = ⟨ν(|g_i + g_e|)(g_i + g_e)⟩_M'))
def massenmittel_A_plummer(M, b, g_e, a0=A0, n_mu=200):
    """Massenmittel der Regel über eine Plummer-Kugel (Skala b) im homogenen Feld g_e, Betrag längs g_e (m/s²):
    je Schale Gauß-Legendre über μ = cos, Schalen mit dM/M = 3x²(1 + x²)^{−5/2}dx, x = r/b."""
    from scipy.integrate import quad
    mu, w = np.polynomial.legendre.leggauss(n_mu)
    X = G * M / (b * b)

    def schale(x):
        gi = X * x / (1 + x * x) ** 1.5
        gz = g_e - gi * mu
        return 3 * x * x / (1 + x * x) ** 2.5 * 0.5 * np.sum(w * nu(np.hypot(gz, gi * np.sqrt(1 - mu**2)), a0) * gz)
    return quad(schale, 0, np.inf, limit=400)[0]


@formel('massenmittel_A_scheibe', 'A = ⟨ν(|g_i + g_e|)(g_i + g_e)⟩_M, dünne Exponentialscheibe, g_e in der Ebene',
        ('A',), ('M', 'R_d', 'a_e', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'Eine Milchstraße als Exponentialscheibe'),
        alternative='massenmittel_A_plummer')
def massenmittel_A_scheibe(M, Rd, g_e, a0=A0, n_je=40, n_phi=720):
    """Massenmittel der Regel über eine dünne Exponentialscheibe mit g_e in ihrer Ebene (Betrag längs g_e, m/s²):
    dM/M = u·e^{−u}du·dφ/2π, u = R/R_d, Gauß-Legendre je Teilstück, φ gleichmäßig."""
    kanten = np.array([0, 0.25, 0.5, 1, 2, 3, 5, 8, 12, 20, 35, 60])
    x, w = np.polynomial.legendre.leggauss(n_je)
    u = np.concatenate([0.5 * (b_ - a_) * x + 0.5 * (b_ + a_) for a_, b_ in zip(kanten[:-1], kanten[1:])])
    wu = np.concatenate([0.5 * (b_ - a_) * w for a_, b_ in zip(kanten[:-1], kanten[1:])])
    ph = np.arange(n_phi) * 2 * math.pi / n_phi
    g = g_expscheibe(M, Rd, u * Rd)
    gx = g_e - g[:, None] * np.cos(ph)[None, :]
    gy = -g[:, None] * np.sin(ph)[None, :]
    return float(np.einsum("i,ij->", wu * u * np.exp(-u), nu(np.hypot(gx, gy), a0) * gx) / n_phi)


@formel('sigma_plummer_jeans', 'σ_v² = ∫r²p·w dr/∫r²ρ·w dr, p = ∫_r^∞ρ·g dr, g = ⟨ν(|g_i + g_e|)(g_i + g_e)⟩ über die Kugel, Fenster R < R_e',
        ('σ_v',), ('M', 'R_e', 'a_e', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'Plummer-Kugel, Jeans-Gleichung, Außenfeld'))
def sigma_plummer_jeans(M, Re, g_e, a0=A0, regel=True, n_r=1500, n_mu=96, lage="kugel"):
    """Sichtlinienstreuung einer Plummer-Kugel (b = R_e), isotrope Jeans-Gleichung, gemittelt im Fenster R < R_e
    (Anteil w jeder Schale im Zylinder). Feld: Monopol der örtlichen Regel im Gesamtfeld mit homogenem g_e (regel,
    lage "kugel"), oder die Regel mit dem Betrag des Gesamtfelds bei fester Lage des äußeren Felds zum inneren
    (lage "gleich": |g| = g_i + g_e, "senkrecht": |g| = √(g_i² + g_e²)), oder Newton g_i (regel False). In m/s."""
    r = np.logspace(math.log10(Re) - 4, math.log10(Re) + 4.5, n_r)
    gi = G * M * r / (r * r + Re * Re) ** 1.5
    rho = (1 + (r / Re) ** 2) ** -2.5
    if regel and lage != "kugel":
        g = np.array([nu_aussenfeld(x, g_e, lage, a0) for x in gi]) * gi
    elif regel:
        mu, wm = np.polynomial.legendre.leggauss(n_mu)
        x = np.sqrt(gi[:, None]**2 + g_e**2 + 2 * gi[:, None] * g_e * mu[None, :])
        g = 0.5 * np.sum(wm * nu(np.maximum(x, 1e-300), a0) * (gi[:, None] + g_e * mu[None, :]), axis=1)
    else:
        g = gi
    f = rho * g
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r) ** 2, 0, 1)))
    return math.sqrt(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r))


@formel('sigma_v_tief', 'σ_v⁴ = (4/81)·GMa₀',
        ('σ_v',), ('M', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'σ_v⁴ = (4/81)·GMa₀'), alternative='sigma_plummer_jeans')
def sigma_v_tief(M, a0=A0):
    """Massengewichtete isotrope Streuung einer druckgestützten Kugel ganz im Tiefbereich, für jedes Profil."""
    return (4 / 81 * G * M * a0) ** 0.25


# ================================================================== Gezeit und Hierarchie der Trichter
@formel('gezeit_sonnenort', '𝒯 aus diag(Φ_RR, Φ_φφ, Φ_zz): a(g_N) = v₀²/R₀, d ln g_N/d ln R = (2 d ln v/d ln R − 1)/(1 − χ/2), Φ_zz = 4πGρ − Φ_RR − Φ_φφ',
        ('𝒯',), ('v', 'r', 'A_O', 'B_O', 'ρ', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'im Sonnensystem ist das die galaktische'), alternative='gezeit_betrag')
def gezeit_sonnenort(v0, R0, A_O, B_O, rho, a0=A0):
    """Gezeit des geglätteten Newton-Felds der Galaxie am Ort der Sonne aus der Kreisbahn (v₀, R₀), den Oort-Konstanten
    (d ln v/d ln R = −(A + B)/(A − B)) und der Dichte der Baryonen am Ort; das Newton-Feld folgt aus der Regel rückwärts,
    d ln ν/d ln a_N = −χ/2."""
    from scipy.optimize import brentq
    g_obs = v0**2 / R0
    gN = brentq(lambda x: a_regel(x, a0) - g_obs, g_obs * 1e-6, g_obs)
    dlnv = -(A_O + B_O) / (A_O - B_O)
    dlngN = (2 * dlnv - 1) / (1 - float(chi(gN, a0)) / 2)
    Prr, Ppp = gN / R0 * dlngN, gN / R0
    return gezeit_betrag(np.diag([Prr, Ppp, 4 * math.pi * G * rho - Prr - Ppp]))


@formel('gezeit_hill', '𝒯_H = 4Ω² − κ² = 4A_O(A_O − B_O)',
        ('𝒯_H',), ('A_O', 'B_O'), 'L', 'modell', 'Baustein 02', ('02-vom-teilchen-bis-zum-stern.md', 'Hill-Radius 1,37 pc'))
def gezeit_hill(A_O, B_O):
    """Gezeit des Hill-Radius r_H = (GM/𝒯_H)^{1/3} (ungeglättetes Feld in der drehenden Scheibe); eine andere Größe als 𝒯."""
    return 4 * A_O * (A_O - B_O)


# ================================================================== Bose-Form und Spektraldimension an den Punkten


@formel('fq_form', 'F_q = Σ((Ȳ_bin − log₁₀ν_Form(g_N; a₀))/σ)², Ȳ_bin Median von log₁₀(g_obs/g_N)',
        ('F_q',), ('Ȳ_bin', 'a₀'), 'L', 'modell', 'Baustein 07', (_F07, 'Fehlerquadratsumme 3- bis 8-fach'))
def fq_form(B, a0, form=None, a0_bins=A0):
    """Fehlerquadratsumme einer Form ν(g_N, a₀) (Standard: die Bose-Form ν des Kerns) an Bin-Medianen B (bin_mediane
    von log₁₀(g_obs/g_N), Mitten bezogen auf a0_bins)."""
    form = nu if form is None else form
    return float(np.sum(((B[:, 1] - np.log10(form(10 ** B[:, 0] * a0_bins, a0))) / B[:, 3]) ** 2))


@formel('a0_an_punkten', 'a₀_pkt = argmin F_q(a₀), a₀ in 10^{±0,5}·10⁻¹⁰',
        ('a₀_pkt',), ('Ȳ_bin',), 'LI', 'modell', 'Baustein 07', (_F07, 'An den SPARC-Punkten (g_N ≤ 3a₀, a₀ frei)'))
def a0_an_punkten(B, form=None, a0_bins=A0):
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda la: fq_form(B, 10**la * 1e-10, form, a0_bins), bounds=(-0.5, 0.5), method="bounded")
    return 10 ** r.x * 1e-10


@formel('nu_ds', 'ν_{d_s} = 1/(1 − e^{−(a_N/a₀)^{1−d_s/2}}); d_s = 1 Bose-Form, tief a ∝ a_N^{d_s/2}',
        ('ν',), ('a_N', 'd_s', 'a₀'), 'LI', 'modell', 'Baustein 07', (_F07, '1 − q ∝ Γ^{1−d_s/2}'), alternative='nu')
def nu_ds(aN, d_s, a0=A0):
    return 1 / (-np.expm1(-(np.asarray(aN) / a0) ** (1 - d_s / 2)))


@formel('ds_anpassung', 'd_s aus ln g_obs = ln(ν_{d_s}(g_N; a₀)·g_N), robust (soft_l1, Skala 0,1), a₀ frei',
        ('d_s',), ('g_obs', 'g_N'), 'LI', 'modell', 'Baustein 07',
        (_F07, 'An den SPARC-Kurven gibt die Form der Regel d_s'))
def ds_anpassung(gN, g, a0=A0):
    """Spektraldimension (d_s in [0,2; 1,9]) mit freiem Faktor der Schwelle (10^{±1})."""
    from scipy.optimize import least_squares
    y = np.log(g)
    r = least_squares(lambda p: np.log(nu_ds(gN, p[0], a0 * 10 ** p[1]) * gN) - y, [1.0, 0.0],
                      bounds=([0.2, -1], [1.9, 1]), loss="soft_l1", f_scale=0.1)
    return float(r.x[0])


@formel('spektraldimension_regel', 'd_s = 2·d ln a/d ln a_N = 2 − χ',
        ('d_s',), ('a_N', 'a₀'), 'LI', 'modell', 'Baustein 03', (_F03, 'd ln a/d ln a_N = 1 − χ/2'),
        alternative='ds_anpassung')
def spektraldimension_regel(aN, a0=A0):
    """Örtliche Steigung der Regel als Spektraldimension; 1 erst für s → 0, 2 in der Newton-Grenze."""
    return 2 - chi(aN, a0)


# ================================================================== Flussdimension, gelesene Dichte, Faltungsmesser
@formel('dichte_gelesen', '4πGρ_X = χΩ² + 4πGρ·(ν(1 − χ/2) − 1), Ω² = νa_N/r',
        ('ρ_X',), ('a_N', 'r', 'ρ', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, '4πGρ_X = χΩ² + 4πGρ_b·∂(a − a_N)/∂a_N'))
def dichte_gelesen(aN, r, rho, a0=A0):
    """Dichte, die ein Beobachter liest, der Newton erwartet (kugelsymmetrisch); eine Ausgabe seiner Rechnung, kein
    Stoff. ∂(a − a_N)/∂a_N = ν(1 − χ/2) − 1 wechselt das Vorzeichen bei a_N = 2,54·a₀."""
    n_, x = nu(aN, a0), chi(aN, a0)
    return (x * n_ * aN / r + 4 * math.pi * G * rho * (n_ * (1 - x / 2) - 1)) / (4 * math.pi * G)


@formel('r_stern_filament', 'r*_fil = 2Gμ/a₀',
        ('r*_fil',), ('μ_fil', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'r*_fil = 2Gμ/a₀'))
def r_stern_filament(mu, a0=A0):
    return 2 * G * mu / a0


@formel('faltungsmesser', 'Δs = v_flach⁴/(3c²a₀)',
        ('Δs',), ('v_flach', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'Δs = v_flach⁴/(3c²a₀)'), alternative='ds')
def faltungsmesser(v, a0=A0):
    """Faltung (und damit baryonische Masse) aus der flachen Rotationsgeschwindigkeit."""
    return v**4 / (3 * c**2 * a0)


@formel('ds_eff', 'Δs_eff = ν·Δs = a·r²/(3c²)',
        ('Δs_eff',), ('M', 'r', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'Δs_eff = a·r²/(3c²) = Δs·ν'))
def ds_eff(M, r, a0=A0):
    """Gezählte Faltung: Faltung mal Rückkehrzahl; die Regel faltet nicht."""
    return nu(G * M / r**2, a0) * ds(M)


@formel('masse_aus_temperatur', 'M_b = (S_T·k_BT/(μ m_p))²/(G·a₀)',
        ('M_b',), ('T', 'S_T', 'μ_m', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'M_b = (S·kT/μm_p)²/(G·a₀)'))
def masse_aus_temperatur(T, S_T, mu, a0=A0):
    """Baryonische Masse eines isothermen Gashalos im Tiefbereich aus seiner Temperatur (∝ T², Newton ∝ T^{3/2})."""
    return (S_T * kB * T / (mu * m_p)) ** 2 / (G * a0)


# ------------------------------------------------------------------ Einsetzen, das der Modelltext als Zahl nennt (03)
@formel('pulsar_rate', 'Ṗ/P = a_r/c (Doppler der Beschleunigung längs der Sichtlinie)',
        ('Ṗ/P',), ('a',), 'L', 'modell', 'Baustein 03', (_F03, 'radial Ṗ/P'), gleiche_form='verlust_gamma')
def pulsar_rate(a):
    """Relative Änderung der Pulsperiode eines Binärpulsars durch die Beschleunigung a längs der Sichtlinie; dieselbe
    Form a/c wie der Verlust Γ = a_N/c der Rückkehr (Ebene I)."""
    return a / c


@formel('warp_viertelzeit', 't_W = (π/2)/ω_p',
        ('t_W',), ('ω_p',), 'L', 'modell', 'Baustein 03', (_F03, 'bei der gemessenen Rate nach 0,14 Mrd. Jahren'))
def warp_viertelzeit(omega_p):
    """Zeit, bis die Knotenlinie bei der Präzessionsrate ω_p um 90° gedreht ist."""
    return (math.pi / 2) / omega_p


@formel('verstaerkung_druck', 'V_σ = (σ_v/σ_v,N)²',
        ('V_σ',), ('σ_v', 'σ_v,N'), 'L', 'modell', 'Baustein 03', (_F03, 'σ_v²/σ_v,N²'))
def verstaerkung_druck(sigma, sigma_N):
    """Verstärkung, die ein druckgestütztes System gegen Newton allein verlangt (σ ∝ √M bei gleichem Aufbau)."""
    return (sigma / sigma_N) ** 2


# ================================================================== Feldform ∇²Φ = ∇·[ν(|∇Φ_N|)∇Φ_N] achsensymmetrisch
# Einheiten der Löser: G = a₀ = 1; Längen in einer gewählten Einheit (Scheiben: r* der Gesamtmasse, dann M = 1).
_FELDFORM_SPEICHER = {}


@formel('feldform_loeser', 'Finite Volumen in (R, z ≥ 0): ∇²Φ_N = 4πGρ (+ Schicht 2πGΣ bei z = 0), dann '
        '∇²Φ = ∇·[ν(|∇Φ_N|)∇Φ_N]; Rand: Punktmasse Φ_N = −M/r, Φ = −M/r − √M·ln(1 − e^{−√M/r})',
        ('g_R', 'g_z', 'a_R'), ('Σ', 'ρ', 'a₀'), 'L', 'hilfe', 'Baustein 07', alternative='g_expscheibe')
def feldform_loeser(quelle, n=200, L=80.0, h0=None):
    """Löser der Feldform für eine achsensymmetrische Quelle (Einheiten G = a₀ = 1). quelle: dict mit "m_kum" (M(<R)
    einer dünnen Schicht in der Ebene, dazu "sigma" Σ(R)) und/oder "rho_punkt" (ρ(R, z) im Volumen), "name" für den
    Speicher.
    Gitter: n × n Zellen, sinh-gedehnt bis L, erste Zelle h0 (Standard min(0,004; L/2·10⁴)). Gibt dict: R, z
    (Zellmitten), gR, gz (Newton, zur Mitte bzw. zur Ebene hin positiv), aR, az (Feldform), M (Gittermasse), sigma.
    Ergebnisse bleiben je (name, n, L, h0) im Speicher."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    h0 = min(0.004, L / 2e4) if h0 is None else h0
    schluessel = (quelle.get("name"), n, L, h0)
    if quelle.get("name") is not None and schluessel in _FELDFORM_SPEICHER:
        return _FELDFORM_SPEICHER[schluessel]
    lo, hi = 1e-3, 60.0
    for _ in range(200):                                  # Dehnung so, dass die erste Zelle h0 breit ist
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if L * math.sinh(mid / n) / math.sinh(mid) > h0 else (lo, mid)
    Rf = L * np.sinh(mid * np.arange(n + 1) / n) / math.sinh(mid)
    zf = Rf.copy()
    Rc, zc = 0.5 * (Rf[1:] + Rf[:-1]), 0.5 * (zf[1:] + zf[:-1])
    dz, dR2 = np.diff(zf), np.diff(Rf**2)
    I, J = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    idx = I + n * J
    cR = 2 * math.pi * Rf[1:-1, None] * dz[None, :] / (Rc[1:] - Rc[:-1])[:, None]
    cZ = math.pi * dR2[:, None] / (zc[1:] - zc[:-1])[None, :]
    cRb = 2 * math.pi * Rf[-1] * dz / (Rf[-1] - Rc[-1])
    cZb = math.pi * dR2 / (zf[-1] - zc[-1])
    zeilen, spalten, werte = [], [], []
    for a_, b_, c_ in ((idx[:-1].ravel(), idx[1:].ravel(), cR.ravel()),
                       (idx[:, :-1].ravel(), idx[:, 1:].ravel(), cZ.ravel())):
        zeilen += [a_, a_, b_, b_]
        spalten += [a_, b_, b_, a_]
        werte += [-c_, c_, -c_, c_]
    zeilen += [idx[-1, :], idx[:, -1]]
    spalten += [idx[-1, :], idx[:, -1]]
    werte += [-cRb, -cZb]
    A = sp.csc_matrix((np.concatenate(werte), (np.concatenate(zeilen), np.concatenate(spalten))), shape=(n * n, n * n))
    lu = splu(A)
    rhs = np.zeros((n, n))
    sig = np.zeros(n)
    M = 0.0
    if "m_kum" in quelle:
        ring = quelle["m_kum"](Rf[1:]) - quelle["m_kum"](Rf[:-1])
        rhs[:, 0] += 2 * math.pi * ring                   # Fluss durch die Ebene: ∂Φ_N/∂z(0⁺) = 2πGΣ
        sig = quelle["sigma"](Rc)
        M += ring.sum()
    if "rho_punkt" in quelle:                           # Zellmittel in z über 4 Gauß-Punkte, in R die Zellmitte
        xg, wg = np.polynomial.legendre.leggauss(4)
        rq = sum(wi / 2 * quelle["rho_punkt"](Rc[:, None], (zf[:-1] + (xi + 1) / 2 * dz)[None, :]) for xi, wi in zip(xg, wg))
        dm = rq * (math.pi * dR2[:, None] * dz[None, :])
        rhs += 4 * math.pi * dm
        M += 2 * dm.sum()                                 # Halbraum, gespiegelt
    Rb, zt, wM = Rf[-1], zf[-1], math.sqrt(M)
    phiN_rand = lambda R, z: -M / np.hypot(R, z)                                            # noqa: E731
    phi_rand = lambda R, z: -M / np.hypot(R, z) - wM * np.log(-np.expm1(-wM / np.hypot(R, z)))  # noqa: E731
    nu_rand = lambda R, z: nu(M / (R * R + z * z), 1.0)                                    # noqa: E731
    b = rhs.copy()
    b[-1, :] -= cRb * phiN_rand(Rb, zc)
    b[:, -1] -= cZb * phiN_rand(Rc, zt)
    F = lu.solve(b.T.ravel()).reshape(n, n).T
    gR, gz = np.empty_like(F), np.empty_like(F)
    gR[1:-1] = (F[2:] - F[:-2]) / (Rc[2:] - Rc[:-2])[:, None]
    gR[0], gR[-1] = (F[1] - F[0]) / (Rc[1] - Rc[0]), gR[-2]
    gz[:, 1:-1] = (F[:, 2:] - F[:, :-2]) / (zc[2:] - zc[:-2])[None, :]
    gz[:, -1] = gz[:, -2]
    gz[:, 0] = 0.5 * (2 * math.pi * sig + (F[:, 1] - F[:, 0]) / (zc[1] - zc[0]))
    nR = nu(np.hypot((F[1:] - F[:-1]) / (Rc[1:] - Rc[:-1])[:, None], 0.5 * (gz[1:] + gz[:-1])), 1.0)
    nZ = nu(np.hypot((F[:, 1:] - F[:, :-1]) / (zc[1:] - zc[:-1])[None, :], 0.5 * (gR[:, 1:] + gR[:, :-1])), 1.0)
    r2 = np.zeros((n, n))
    fl = nR * cR * (F[1:] - F[:-1])
    r2[:-1] += fl
    r2[1:] -= fl
    fl = nZ * cZ * (F[:, 1:] - F[:, :-1])
    r2[:, :-1] += fl
    r2[:, 1:] -= fl
    r2[-1, :] += nu_rand(Rb, zc) * cRb * (phiN_rand(Rb, zc) - F[-1, :]) - cRb * phi_rand(Rb, zc)
    r2[:, -1] += nu_rand(Rc, zt) * cZb * (phiN_rand(Rc, zt) - F[:, -1]) - cZb * phi_rand(Rc, zt)
    P = lu.solve(r2.T.ravel()).reshape(n, n).T
    aR, az = np.empty_like(P), np.empty_like(P)
    aR[1:-1] = (P[2:] - P[:-2]) / (Rc[2:] - Rc[:-2])[:, None]
    aR[0], aR[-1] = (P[1] - P[0]) / (Rc[1] - Rc[0]), aR[-2]
    az[:, 1:-1] = (P[:, 2:] - P[:, :-2]) / (zc[2:] - zc[:-2])[None, :]
    az[:, -1] = az[:, -2]
    az[:, 0] = 0.5 * (2 * math.pi * sig * nu(np.hypot(gR[:, 0], 2 * math.pi * sig), 1.0)
                      + (P[:, 1] - P[:, 0]) / (zc[1] - zc[0]))
    aus = dict(R=Rc, z=zc, gR=gR, gz=gz, aR=aR, az=az, M=M, sigma=sig)
    if quelle.get("name") is not None:
        _FELDFORM_SPEICHER[schluessel] = aus
    return aus


@formel('scheibe_quelle', 'dünne Scheibe der Masse 1 mit Skala h: exponentiell Σ = e^{−R/h}/(2πh²), '
        'Kuzmin Σ = h/(2π(R² + h²)^{3/2})', ('Σ',), ('M', 'R_d', 'r'), 'L', 'hilfe', 'Baustein 07')
def scheibe_quelle(profil, h):
    """Quelle für feldform_loeser (Masse 1, Längen in r*): dict mit m_kum, sigma, name; profil "exp" oder "kuzmin"."""
    if profil == "exp":
        return dict(name=("exp", h), m_kum=lambda R: 1 - (1 + R / h) * np.exp(-R / h),
                    sigma=lambda R: np.exp(-R / h) / (2 * math.pi * h * h))
    if profil == "kuzmin":
        return dict(name=("kuzmin", h), m_kum=lambda R: 1 - h / np.sqrt(R * R + h * h),
                    sigma=lambda R: h / (2 * math.pi * (R * R + h * h) ** 1.5))
    raise ValueError(f"Profil {profil!r}")


@formel('eta_scheibe', 'η = ν(g_R)·g_R / a_R,Feldform in der Ebene einer dünnen Scheibe bei R = x·h',
        ('η',), ('R_d', 'r*', 'r', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'Exponentialscheibe, Feldform gegen radiale Regel'))
def eta_scheibe(profil, h, x, n=200):
    """M_dyn/M_lens in der Ebene einer dünnen Scheibe (Masse 1, h in r*): Bahn nach der radialen Regel ν(g_R)·g_R
    gegen die Feldform, die die Metrik und damit das Licht trägt. Licht sieht das 1/η-Fache der Bahnbeschleunigung."""
    f = feldform_loeser(scheibe_quelle(profil, h), n)
    gR, aF = f["gR"][:, 0], f["aR"][:, 0]
    return float(np.interp(x * h, f["R"], nu(gR, 1.0) * gR / aF))


@formel('gleichstand_scheibe', 'x_gl: η(x_gl·R_d) = 1 zwischen R_d und 8R_d (Feldform gleich radialer Regel)',
        ('x_gl',), ('R_d', 'r*', 'a₀'), 'L', 'modell', 'Baustein 07', (_F07, 'Gleichstand bei etwa'))
def gleichstand_scheibe(h, n=200):
    """Lage des Gleichstands von Feldform und radialer Regel in der Ebene einer Exponentialscheibe (h = R_d in r*), in
    Einheiten R_d; lineare Interpolation zwischen den Zellmitten."""
    f = feldform_loeser(scheibe_quelle("exp", h), n)
    x = f["R"] / h
    eta = nu(f["gR"][:, 0], 1.0) * f["gR"][:, 0] / f["aR"][:, 0]
    k = np.where((x[:-1] > 1) & (x[:-1] < 8) & (eta[:-1] >= 1) & (eta[1:] < 1))[0]
    if not k.size:
        return float("nan")
    i = k[0]
    return float(x[i] + (eta[i] - 1) / (eta[i] - eta[i + 1]) * (x[i + 1] - x[i]))


@formel('feldform_naeherung', 'ε_F = a_R,Feldform / (ν(√(g_R² + (2πGΣ)²))·g_R) in der Ebene bei R = x·h',
        ('ε_F',), ('R_d', 'r*', 'r', 'a₀'), 'L', 'modell', 'Baustein 07', (_F07, 'die Feldform auf 0,02 dex'))
def feldform_naeherung(profil, h, x, n=200):
    """Gelöste Feldform gegen die örtliche Näherung mit dem senkrechten Scheibenfeld (f = 1) in der Ebene."""
    f = feldform_loeser(scheibe_quelle(profil, h), n)
    gR = f["gR"][:, 0]
    naeh = nu(np.hypot(gR, 2 * math.pi * f["sigma"]), 1.0) * gR
    return float(np.interp(x * h, f["R"], f["aR"][:, 0] / naeh))


# ================================================================== Sonnenort: Newton-Feld, Dynamikmasse, Sattelpunkt
@formel('a_N_aus_kreisbahn', 'a_N: ν(a_N)·a_N = v²/r (die Regel rückwärts)',
        ('a_N',), ('v', 'r', 'a₀'), 'L', 'modell', 'Baustein 03', (_F03, 'das geglättete Feld der Galaxie ist dort'),
        alternative='neigung_stoerung')
def a_N_aus_kreisbahn(v, r, a0=A0):
    """Newton-Feld der Baryonen am Ort einer Kreisbahn (v, r), aus der gemessenen Beschleunigung v²/r über die Regel."""
    from scipy.optimize import brentq
    g = v * v / r
    return brentq(lambda x: a_regel(x, a0) - g, g * 1e-8, g, xtol=g * 1e-15, rtol=1e-14)


@formel('sattelpunkt', 'd_S = D·√M₂/(√M₁ + √M₂): GM₂/d_S² = GM₁/(D − d_S)²',
        ('r',), ('M', 'r'), 'L', 'modell', 'Baustein 03', (_F03, 'Am Sattelpunkt Sonne–Erde'))
def sattelpunkt(GM1, GM2, D):
    """Abstand des Nullpunkts der Newton-Felder zweier Körper vom leichteren (GM₂) auf der Verbindungslinie."""
    return D * math.sqrt(GM2) / (math.sqrt(GM1) + math.sqrt(GM2))


@formel('dynamikmasse_ort', 'Σ_z = a_z(z₁)/(2πG): Kern ν_z² = ν_z,N² + (ν − 1)Ω_N² → Σ_b + (ν(g_R) − 1)g_R·z₁/(2πGR); '
        'örtliche Regel a_z = ν(|g_N|)·g_z → ν(√(g_R² + (2πGΣ_b)²))·Σ_b (Schicht) bis ν(g_R)·Σ_b (Ebene)',
        ('Σ_z',), ('Σ', 'g_R', 'z_h', 'r', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'Die Dynamikmasse am Sonnenort'))
def dynamikmasse_ort(Sigma_b, gR, z1, R, form, a0=A0):
    """Gelesene Flächendichte der Säule |z| < z₁ am Ort R einer Scheibe (kg/m²): Baryonen Σ_b, Newton-Feld g_R der
    Scheibe in der Ebene. form "kern" (radialer Kern e_r⊗e_r), "ort_schicht" (örtliche Regel mit dem Feld der Säule),
    "ort_ebene" (örtliche Regel mit g_R allein)."""
    gz = 2 * math.pi * G * Sigma_b
    if form == "kern":
        return Sigma_b + (float(nu(gR, a0)) - 1) * gR * z1 / R / (2 * math.pi * G)
    if form == "ort_schicht":
        return float(nu(math.hypot(gR, gz), a0)) * Sigma_b
    if form == "ort_ebene":
        return float(nu(gR, a0)) * Sigma_b
    raise ValueError(f"Form {form!r}")


@formel('knotenpraezession_scheibe', 'ω_p = Ω − ω_z in der Ebene einer Exponentialscheibe mit Schicht der Höhe h: '
        'ω_z,N² = 4πGΣ/(2h) − (dv_N²/dR)/R, Ω² = ν(g_R)·Ω_N², ω_z nach Kern- oder Feldform',
        ('ω_p',), ('M', 'R_d', 'z_h', 'r', 'a₀'), 'L', 'modell', 'Baustein 07',
        (_F07, 'Präzessionsraten geneigter Ringe unterscheiden Kern- und Feldform'), alternative='praezession_ring')
def knotenpraezession_scheibe(M, Rd, h, R, form, a0=A0):
    """Knotenpräzession schwach geneigter Ringe (s⁻¹, negativ = rückläufig) in der Ebene einer dünnen Exponentialscheibe
    (Masse M, Skala R_d, Mittelebenendichte Σ/(2h)); form "kern" (Kernform) oder "ort" (örtliche Regel = Feldform in
    der Ebene, ω_z² = ν·ω_z,N²)."""
    dR = 1e-4 * R
    v2 = lambda r: float(g_expscheibe(M, Rd, r)) * r                       # noqa: E731
    dv2 = (v2(R + dR) - v2(R - dR)) / (2 * dR)
    rho0 = M / (2 * math.pi * Rd * Rd) * math.exp(-R / Rd) / (2 * h)
    Om_N = math.sqrt(v2(R)) / R
    om_zN = math.sqrt(4 * math.pi * G * rho0 - dv2 / R)
    return knotenpraezession_frequenzen(Om_N, om_zN, float(nu(v2(R) / R, a0)), form)


@formel('knotenpraezession_frequenzen', 'ω_p = √ν·Ω_N − ω_z(ω_z,N, Ω_N, ν; Form)',
        ('ω_p',), ('Ω_N', 'ω_z,N', 'ν'), 'L', 'modell', 'Baustein 07',
        (_F07, 'die √ν-fache der Newtonschen'), alternative='praezession_ring')
def knotenpraezession_frequenzen(Om_N, om_zN, nu_, form):
    """Knotenpräzession schwach geneigter Ringe in der Mittelebene aus den Newton-Frequenzen am Ort; Ω² = ν·Ω_N² in
    beiden Formen, ω_z nach vertikalfrequenz ("kern" oder "ort"; "ort" gibt genau √ν·(Ω_N − ω_z,N))."""
    return math.sqrt(nu_) * Om_N - vertikalfrequenz(om_zN, Om_N, nu_, form)


@formel('milchstrasse_quelle', 'ρ = Bulge ρ₀(1 + r′/r₀)^{−α}e^{−(r′/r_c)²} + Σ Sternscheiben Σ₀/(2z_d)e^{−R/R_d−|z|/z_d} '
        '+ Σ Gasscheiben Σ₀/(4z_d)e^{−R_m/R−R/R_d}sech²(z/2z_d)', ('ρ',), ('Σ₀', 'R_d', 'z_h'), 'L', 'hilfe',
        'McMillan 2017 (baryonisch)')
def milchstrasse_quelle(teile, name="milchstrasse", faktor_scheibe=1.0):
    """Quelle für feldform_loeser in den Einheiten G = a₀ = 1, Länge 1 kpc (Masse in a₀·kpc²/G). teile: dict mit
    "bulge" (ρ₀ [M☉/kpc³], r₀, r_cut, α, q), "duenn", "dick" (Σ₀ [M☉/kpc²], R_d, z_d), "HI", "H2" (Σ₀, R_d, R_m, z_d);
    Längen in kpc. faktor_scheibe skaliert die beiden Sternscheiben. In z je Zelle über 4 Gauß-Punkte gemittelt."""
    mu = M_sun / (A0 * kpc**2 / G)                                           # M☉ in Masseneinheiten des Lösers
    b = teile["bulge"]

    def rho_punkt(R, z):
        rp = np.sqrt(R * R + (z / b[4])**2)
        r = b[0] / (1 + rp / b[1])**b[3] * np.exp(-(rp / b[2])**2)
        for k in ("duenn", "dick"):
            S0, Rd, zd = teile[k]
            r = r + faktor_scheibe * S0 / (2 * zd) * np.exp(-R / Rd - np.abs(z) / zd)
        for k in ("HI", "H2"):
            S0, Rd, Rm, zd = teile[k]
            e = np.exp(-np.abs(z) / zd)
            r = r + S0 / (4 * zd) * np.exp(-Rm / np.maximum(R, 1e-6) - R / Rd) * 4 * e / (1 + e)**2
        return r * mu
    return dict(name=(name, faktor_scheibe), rho_punkt=rho_punkt)


@formel('frequenzen_loesung', 'Ω_N² = g_R/R, ω_z,N² = 4πGρ(R, 0) − (1/R)·d(R g_R)/dR in der Mittelebene einer Lösung',
        ('Ω_N', 'ω_z,N'), ('g_R', 'ρ', 'r'), 'L', 'hilfe', 'Baustein 07', alternative='umlauf_kepler')
def frequenzen_loesung(loesung, rho0, R):
    """(Ω_N, ω_z,N, g_R) am Radius R in der Mittelebene aus einer Lösung von feldform_loeser (Einheiten des Lösers);
    rho0 Dichte in der Mittelebene bei R. d(R g_R)/dR über ±5 % in R."""
    Rg, gR = loesung["R"], loesung["gR"][:, 0]
    Rg_ = lambda x: x * np.interp(x, Rg, gR)                                 # noqa: E731
    d = 0.05 * R
    ab = (Rg_(R + d) - Rg_(R - d)) / (2 * d)
    g = float(np.interp(R, Rg, gR))
    return math.sqrt(g / R), math.sqrt(4 * math.pi * rho0 - ab / R), g


@formel('saeule_loesung', 'Σ_z = a_z(R, z₁)/(2πG) aus einer Lösung (Feldform oder Newton)',
        ('Σ_z',), ('g_z', 'z_h', 'r'), 'L', 'modell', 'Baustein 07', (_F07, 'mit der Feldform bei 69–92'),
        alternative='dynamikmasse_ort')
def saeule_loesung(loesung, R, z1, feld="feldform"):
    """Gelesene Flächendichte der Säule |z| < z₁ bei R aus einer Lösung von feldform_loeser (Einheiten des Lösers:
    Masse je Längeneinheit², G = 1); feld "feldform" (a_z) oder "newton" (g_z). Bilineare Interpolation."""
    from scipy.interpolate import RegularGridInterpolator
    A = loesung["az"] if feld == "feldform" else loesung["gz"]
    return float(RegularGridInterpolator((loesung["R"], loesung["z"]), A)([[R, z1]])[0]) / (2 * math.pi)


@formel('warp_aufwicklung', 't_Wick = Δφ/|ω_p(R₁) − ω_p(R₂)|', ('t_Wick',), ('ω_p',), 'L', 'modell', 'Baustein 03',
        (_F03, 'Die freie Präzession wickelt den Warp zudem auf'))
def warp_aufwicklung(omega1, omega2, dphi):
    """Zeit, bis die Knotenlinien zweier Ringe mit den Präzessionsraten ω₁, ω₂ um Δφ auseinanderlaufen."""
    return dphi / abs(omega1 - omega2)


@formel('quadrupol_glaettung', 'Q₂ = −(3/2)(g_e/r_e)·J(β, λ), β = g_e/a₀, λ = ℓ/r_e, r_e = √(GM/g_e); '
        'J = ½∫∫(ν(β|h e_r + e_z|) − 1)·k(ρ, θ)·ρ sin θ dρ dθ + (4/15)(−χ/2)ν/λ³, h = ρ/(ρ² + λ²)^{3/2}',
        ('𝒯_Q',), ('M', 'a_e', 'ℓ', 'a₀'), 'L', 'modell', 'Baustein 07', (_F07, 'Im Planetensystem bleibt ein Restfeld'))
def quadrupol_glaettung(GM, g_e, ell, a0=A0, n=240):
    """Quadrupol (s⁻²) des Zusatzfelds der Regel um einen Körper, dessen Feld über seine Glättungslänge ℓ (Plummer)
    geglättet ist, im gleichförmigen Feld g_e der Elternebene; Kern k = −3(3cos²θ − 1)h/ρ² + (9cosθ − 15cos³θ)/ρ².
    Gauß-Legendre in ln ρ (−10 … 12, Stützstellen um ln λ) und θ (vier Teilstücke)."""
    re = math.sqrt(GM / g_e)
    beta, lam = g_e / a0, ell / re
    xg, wg = np.polynomial.legendre.leggauss(n)
    L = math.log(lam)
    bps = sorted(b for b in {-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 3.5, 6.0, 9.0, 12.0,
                             L - 1.5, L - 0.7, L - 0.3, L, L + 0.3, L + 0.7, L + 1.5} if -10.0 <= b <= 12.0)
    tseg = [(0.0, 0.5 * math.pi), (0.5 * math.pi, 0.9 * math.pi), (0.9 * math.pi, 0.99 * math.pi), (0.99 * math.pi, math.pi)]
    J = 0.0
    for u0, u1 in zip(bps[:-1], bps[1:]):
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0)
        wu = 0.5 * (u1 - u0) * wg
        rho_ = np.exp(u)
        for t0, t1 in tseg:
            th = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0)
            wt = 0.5 * (t1 - t0) * wg
            R_, TH = np.meshgrid(rho_, th, indexing="ij")
            cc = np.cos(TH)
            h = R_ / (R_**2 + lam**2)**1.5
            gm = np.sqrt(np.maximum(h**2 + 2 * cc * h + 1, 0.0))
            k = -3 * (3 * cc**2 - 1) * h / R_**2 + (9 * cc - 15 * cc**3) / R_**2
            J += np.einsum("i,j,ij->", wu, wt, 0.5 * (nu(beta * gm, 1.0) - 1.0) * k * R_ * np.sin(TH))
    J += (4 / 15) * (-float(chi(beta, 1.0)) / 2) * float(nu(beta, 1.0)) / lam**3
    return -1.5 * (g_e / re) * J


# ================================================================== SPARC: profilierte Lesung der Schwelle
_PROFIL_SPEICHER = {}


@formel('sparc_profil_satz', 'Galaxien der Tafel mit V_obs > 0, e_V > 0 an mindestens 5 Punkten; Fehler e_D/D ≥ 0,01, e_i ≥ 1°',
        ('V_obs', 'e_V'), ('Q_fl', 'i_inc'), 'L', 'hilfe', 'SPARC (Lelli, McGaugh, Schombert 2016)', alternative='sparc_kurve')
def sparc_profil_satz(tafel):
    """Aufbereitete Galaxien für die profilierte Lesung: gepolsterte Felder (Galaxien × Punkte) mit Maske. tafel: Liste
    aus sparc_tafel. Gibt dict: name, fD, lgo = ln(V_obs²/R), sv = 2e_V/V_obs, gas, disk, bul (V²/R in m/s², Gas mit
    Vorzeichen), eD, inc, ei (rad), maske, SB0."""
    gal = []
    for g in tafel:
        m = (g["Vo"] > 0) & (g["R"] > 0) & (g["eV"] > 0)
        if m.sum() >= 5:
            gal.append((g, m))
    n, p = len(gal), max(int(m.sum()) for _, m in gal)
    z = lambda: np.zeros((n, p))                                              # noqa: E731
    aus = dict(name=[g["nm"] for g, _ in gal], fD=np.array([g["fD"] for g, _ in gal]), lgo=z(), sv=np.ones((n, p)),
               gas=z(), disk=z(), bul=z(), maske=np.zeros((n, p), bool), SB0=np.array([g["SB0"] for g, _ in gal]),
               eD=np.array([max(g["eD"] / g["D"], 0.01) for g, _ in gal]),
               inc=np.radians([g["inc"] for g, _ in gal]), ei=np.radians([max(g["eInc"], 1.0) for g, _ in gal]))
    for k, (g, m) in enumerate(gal):
        j = int(m.sum())
        R = g["R"][m]
        aus["lgo"][k, :j] = np.log(g["Vo"][m]**2 / R)
        aus["sv"][k, :j] = 2 * g["eV"][m] / g["Vo"][m]
        aus["gas"][k, :j] = g["Vg"][m] * np.abs(g["Vg"][m]) / R
        aus["disk"][k, :j] = g["Vd"][m]**2 / R
        aus["bul"][k, :j] = g["Vb"][m]**2 / R
        aus["maske"][k, :j] = True
    return aus


@formel('profil_chi2', 'χ²_prof = min über (u_D, δ_i, u_Υ) von Σ[(ln g_obs − u_D + 2 ln(sin i/sin(i + δ_i)) '
        '− ln ν(g_N)g_N)/σ]² + (u_D/e_D)² + (δ_i/e_i)² + (u_Υ/σ_Υ)², σ² = (2e_V/V)² + σ_int², Υ = Υ_d·e^{u_Υ}',
        ('F_q',), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'D_A', 'i_inc', 'Υ_d', 'a₀'), 'LI', 'modell', 'Baustein 03',
        (_F03, 'mit Abstand, Neigung und Υ_d als Störgrößen'), alternative='fq_form')
def profil_chi2(satz, a0, ups_d, sigma_int, verlust="linear", sigma_ups=None, start=None, n_it=60, ups_b=_UPS_B,
                mit_stoer=False):
    """Profilierte Fehlerquadratsumme je Galaxie bei der Schwelle a0 (Levenberg-Marquardt, je Galaxie drei Störgrößen,
    alle Galaxien zugleich). verlust "linear" (Summe r²) oder "soft_l1" (Summe 2(√(1 + r²) − 1), auch für die Vorgaben).
    Abstand verschiebt ln g_obs (g_N hängt nicht am Abstand), Neigung skaliert V_obs, Υ skaliert Scheibe und Bulge
    (Υ_b Standard aus SPARC_SATZ). sigma_ups: Streuung von ln Υ (Standard 0,1 dex); start: Störgrößen (Galaxien × 3).
    Gibt χ²_prof je Galaxie, mit mit_stoer das Paar (χ²_prof, Störgrößen)."""
    su = math.log(10) * SPARC_SATZ["ups_streuung_dex"] if sigma_ups is None else sigma_ups
    mk, n = satz["maske"], len(satz["name"])
    sig = np.sqrt(satz["sv"]**2 + sigma_int**2)
    x = np.zeros((n, 3)) if start is None else start.copy()
    ab = np.stack([satz["eD"], satz["ei"], np.full(n, su)], 1)

    def rechne(x):
        uD, di, uY = x[:, 0:1], x[:, 1:2], x[:, 2:3]
        inc = np.clip(satz["inc"][:, None] + di, 0.2, math.pi / 2)
        st = np.exp(uY) * (ups_d * satz["disk"] + ups_b * satz["bul"])
        gb_roh = satz["gas"] + st
        gb = np.maximum(gb_roh, 1e-14)
        lgo = satz["lgo"] - uD + 2 * np.log(np.sin(satz["inc"])[:, None] / np.sin(inc))
        r = np.where(mk, (lgo - np.log(a_regel(gb, a0))) / sig, 0.0)
        dla = 1 - chi(gb, a0) / 2
        frei = (satz["inc"][:, None] + di > 0.2) & (satz["inc"][:, None] + di < math.pi / 2)
        Jp = np.stack([-np.ones_like(r), np.broadcast_to(np.where(frei, -2 / np.tan(inc), 0.0), r.shape),
                       -dla * np.where(gb_roh > 1e-14, st / gb, 0.0)], -1) / sig[..., None]
        Jp = np.where(mk[..., None], Jp, 0.0)
        rv = x / ab
        return r, Jp, rv

    def kosten(r, rv):
        if verlust == "linear":
            return np.sum(r * r, 1) + np.sum(rv * rv, 1)
        return np.sum(2 * (np.sqrt(1 + r * r) - 1) * mk, 1) + np.sum(2 * (np.sqrt(1 + rv * rv) - 1), 1)

    r, Jp, rv = rechne(x)
    c_alt = kosten(r, rv)
    lam = np.full(n, 1e-3)
    for _ in range(n_it):
        w = np.ones_like(r) if verlust == "linear" else 1 / np.sqrt(1 + r * r)
        wv = np.ones_like(rv) if verlust == "linear" else 1 / np.sqrt(1 + rv * rv)
        H = np.einsum("np,npi,npj->nij", w, Jp, Jp) + np.einsum("ni,ij->nij", wv / ab**2, np.eye(3))
        gr = np.einsum("np,np,npi->ni", w, r, Jp) + wv * rv / ab
        Hd = H + lam[:, None, None] * np.einsum("nii->ni", H)[:, :, None] * np.eye(3)
        dx = -np.linalg.solve(Hd, gr[..., None])[..., 0]
        xn = x + dx
        rn, Jn, rvn = rechne(xn)
        cn = kosten(rn, rvn)
        besser = cn <= c_alt
        gewinn = np.max(np.where(besser, c_alt - cn, 0.0))
        x = np.where(besser[:, None], xn, x)
        r, Jp, rv = (np.where(besser[:, None], rn, r), np.where(besser[:, None, None], Jn, Jp),
                     np.where(besser[:, None], rvn, rv))
        c_alt = np.where(besser, cn, c_alt)
        lam = np.where(besser, lam / 3, lam * 10)
        if gewinn < 1e-9:
            break
    return (c_alt, x) if mit_stoer else c_alt


@formel('a0_profiliert', 'a₀_prof = argmin Σ_Galaxien χ²_prof(a₀) (Parabel im Gitter ln a₀), Fehler: Bootstrap über '
        'Galaxien', ('a₀_prof',), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'D_A', 'i_inc', 'Υ_d'), 'LI', 'modell',
        'Baustein 03', (_F03, 'mit Abstand, Neigung und Υ_d als Störgrößen'))
def a0_profiliert(satz, ups_d, sigma_int, verlust="linear", auswahl=None, n_boot=1000, seed=7, sigma_ups=None,
                  ups_b=_UPS_B):
    """Schwelle aus allen Kurven zugleich mit Abstand, Neigung und Υ je Galaxie als Störgrößen (profil_chi2), Gitter
    a₀ = 1,1·10⁻¹⁰·10^{−0,25 … 0,20} (46 Punkte). auswahl: Indizes der Galaxien (Standard alle). Gibt (a₀, Streuung
    über n_boot Ziehungen ganzer Galaxien, eigener Generator mit seed; n_boot = 0 gibt nur a₀). Tafeln bleiben
    je Aufruf im Speicher."""
    gitter = A0 * 10 ** np.linspace(-0.25, 0.20, 46)
    schluessel = (tuple(satz["name"]), ups_d, sigma_int, verlust, sigma_ups, ups_b)
    if schluessel not in _PROFIL_SPEICHER:
        T, x = [], None
        for a0 in gitter:
            c2, x = profil_chi2(satz, a0, ups_d, sigma_int, verlust, sigma_ups, start=x, ups_b=ups_b, mit_stoer=True)
            T.append(c2)
        _PROFIL_SPEICHER[schluessel] = np.array(T).T
    T = _PROFIL_SPEICHER[schluessel]
    idx = np.arange(T.shape[0]) if auswahl is None else np.asarray(auswahl)
    lx, h = np.log(gitter), math.log(gitter[1] / gitter[0])

    def best(c2):
        k = min(max(int(np.argmin(c2)), 1), len(c2) - 2)
        y0, y1, y2 = c2[k - 1:k + 2]
        return math.exp(lx[k] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2)))
    b = best(T[idx].sum(0))
    if not n_boot:
        return b
    rng = np.random.default_rng(seed)
    bs = [best(T[idx[rng.integers(0, idx.size, idx.size)]].sum(0)) for _ in range(int(n_boot))]
    return b, float(np.std(bs))


@formel('a0_je_drittel', 'a₀ je Drittel von Σ₀ = Υ_d·SB₀: argmin der Summe der χ²_prof (Abstand, Neigung als Störgrößen, '
        'Υ_d fest, σ = 2e_V/V); Anstieg d log a₀/d log Σ₀ über die drei Drittel, Fehler Bootstrap über Galaxien',
        ('a₀_gem', 'k_a'), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'D_A', 'i_inc', 'SB₀', 'Υ_d'), 'LI', 'modell',
        'Baustein 03', (_F03, 'Gemeinsam je Drittel der Flächendichte angepasst'))
def a0_je_drittel(satz, ups_d, n_boot=300, seed=29092026):
    """Gemeinsame Schwelle je Drittel der zentralen Flächendichte (Galaxien mit SB₀ > 0), Gitter 0,3–4·10⁻¹⁰ (49
    Punkte, Parabel in ln a₀). Gibt (a₀ unteres, mittleres, oberes Drittel, Anstieg je Dekade, Streuung des Anstiegs
    über n_boot Ziehungen ganzer Galaxien mit der ursprünglichen Einteilung; n_boot = 0 ohne diese)."""
    gitter = np.exp(np.linspace(math.log(0.3e-10), math.log(4e-10), 49))
    schluessel = ("drittel", tuple(satz["name"]), ups_d)
    if schluessel not in _PROFIL_SPEICHER:
        T, x = [], None
        for a0 in gitter:
            c2, x = profil_chi2(satz, a0, ups_d, 0.0, "linear", sigma_ups=1e-6, start=x, mit_stoer=True)
            T.append(c2)
        _PROFIL_SPEICHER[schluessel] = np.array(T).T
    T = _PROFIL_SPEICHER[schluessel]
    sel = np.where(satz["SB0"] > 0)[0]
    xs = np.log10(ups_d * satz["SB0"][sel])
    grp = np.digitize(xs, np.percentile(xs, [100 / 3, 200 / 3]))
    lx, h = np.log(gitter), math.log(gitter[1] / gitter[0])

    def lg_best(c2):
        k = min(max(int(np.argmin(c2)), 1), len(c2) - 2)
        y0, y1, y2 = c2[k - 1:k + 2]
        return (lx[k] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2))) / math.log(10)

    def drittel(j):
        a = np.array([lg_best(T[sel[j[grp[j] == k]]].sum(0)) for k in range(3)])
        xm = np.array([xs[j[grp[j] == k]].mean() for k in range(3)])
        return a, np.polyfit(xm, a, 1)[0]
    a, k = drittel(np.arange(sel.size))
    if not n_boot:
        return 10 ** a[0], 10 ** a[1], 10 ** a[2], float(k)
    rng = np.random.default_rng(seed)
    bs = [drittel(rng.integers(0, sel.size, sel.size))[1] for _ in range(int(n_boot))]
    return 10 ** a[0], 10 ** a[1], 10 ** a[2], float(k), float(np.std(bs))


@formel('a0_senkung_aussenfeld', 'a₀_e = argmin Σ w·[ln(ν(|g_N + g_e|)·g_N) − ln(ν(g_N; a₀_e)·g_N)]² an den SPARC-Punkten',
        ('a₀_e',), ('V_gas', 'V_disk', 'V_bul', 'e_V', 'Υ_d', 'a_e', 'a₀'), 'LI', 'modell', 'Baustein 03',
        (_F03, 'Das äußere Feld senkt den gelesenen Wert'))
def a0_senkung_aussenfeld(satz, ups_d, g_e, lage, gewichtet, a0=A0):
    """Schwelle, die man an den Punkten der SPARC-Kurven (g_N aus den eigenen Baryonen) ohne äußeres Feld liest, wenn
    die Regel dort mit dem äußeren Feld g_e (lage "senkrecht" oder "gleich") wirkt; Gewicht (2e_V/V)⁻² oder gleich."""
    from scipy.optimize import minimize_scalar
    mk = satz["maske"]
    gb = np.maximum((satz["gas"] + ups_d * satz["disk"] + _UPS_B * satz["bul"])[mk], 1e-15)
    w = satz["sv"][mk] ** -2 if gewichtet else np.ones(gb.size)
    obs = np.log(np.array([nu_aussenfeld(x, g_e, lage, a0) for x in gb]) * gb)
    r = minimize_scalar(lambda la: np.sum(w * (obs - np.log(nu(gb, a0 * math.exp(la)) * gb))**2), bounds=(-3, 1),
                        method="bounded", options=dict(xatol=1e-10))
    return a0 * math.exp(r.x)


@formel('sigma_dyn_rest', 'ε_Σ = ⟨log₁₀(Σ_dyn,gem/Σ_dyn,Regel)⟩ je Galaxie (Toomre-Integral), Anstieg k_Σ gegen log₁₀Σ₀',
        ('ε_Σ', 'k_Σ'), ('V_obs', 'V_gas', 'V_disk', 'V_bul', 'SB₀', 'Υ_d', 'a₀'), 'L', 'modell', 'Baustein 03',
        (_F03, 'Zentrale Flächendichte aus der Rotationskurve'))
def sigma_dyn_rest(tafel, ups_d, a0=A0):
    """Galaxie für Galaxie: Toomre-Integral der gemessenen Kurve gegen das der Regel mit den eigenen Baryonen
    (radiale Form, V_bar² ≥ 10⁻⁶ m²/s²). Gibt (Mittel, Fehler des Mittels, Anstieg je Dekade Σ₀ = Υ_d·SB₀, Fehler des
    Anstiegs aus der Streuung um die Gerade)."""
    S0, res = [], []
    for g in tafel:
        Vb2 = v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], ups_d)
        gN = np.clip(Vb2, 1e-6, None) / g["R"]
        So = toomre_sigma_dyn(g["R"], g["Vo"])
        Sm = toomre_sigma_dyn(g["R"], np.sqrt(nu(gN, a0) * gN * g["R"]))
        if ups_d * g["SB0"] > 0 and So > 0 and Sm > 0:
            S0.append(ups_d * g["SB0"])
            res.append(math.log10(So / Sm))
    x, y = np.log10(S0), np.array(res)
    A = np.vstack([x - x.mean(), np.ones_like(x)]).T
    cf, rs, *_ = np.linalg.lstsq(A, y, rcond=None)
    return (float(y.mean()), float(y.std(ddof=1) / math.sqrt(y.size)), float(cf[0]),
            float(math.sqrt(rs[0] / (y.size - 2) / np.sum((x - x.mean())**2))))


@formel('rest_gasanteil', 'k_gas = d(Rest − Median der Galaxie)/d f_gas über alle Punkte, f_gas = V_gas|V_gas|/V_bar²',
        ('k_gas',), ('V_obs', 'V_gas', 'V_disk', 'V_bul', 'Υ_d', 'σ_v', 'a₀'), 'L', 'modell', 'Baustein 03',
        (_F03, 'Die Regel hängt nicht am Zustand der Materie'))
def rest_gasanteil(tafel, ups_d, sig, a0=A0, n_boot=300, seed=1):
    """Änderung des Rests log₁₀(g_obs/(ν(g_N)g_N)) mit dem örtlichen Gasanteil innerhalb der Galaxien (Nullpunkt je
    Galaxie abgezogen); Punkte mit V_bar² > 0, V_obs > 0, mindestens 5 je Galaxie; f_gas auf [−0,5; 1,5] begrenzt.
    Gibt (Anstieg, Streuung über n_boot Ziehungen ganzer Galaxien); n_boot = 0 nur den Anstieg."""
    F, Rr, G_ = [], [], []
    for j, g in enumerate(tafel):
        gas = g["Vg"] * np.abs(g["Vg"])
        Vb2 = v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], ups_d)
        m = (Vb2 > 0) & (g["Vo"] > 0)
        if m.sum() < 5:
            continue
        gN = Vb2[m] / g["R"][m]
        r = np.log10(g_gemessen(g["R"][m], g["Vo"][m], g["Rd"], sig) / (nu(gN, a0) * gN))
        F.append(np.clip(gas[m] / Vb2[m], -0.5, 1.5))
        Rr.append(r - np.median(r))
    steig = lambda fs, rs: np.polyfit(np.concatenate(fs) - np.concatenate(fs).mean(), np.concatenate(rs), 1)[0]  # noqa: E731
    if not n_boot:
        return float(steig(F, Rr))
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(int(n_boot)):
        k = rng.integers(0, len(F), len(F))
        bs.append(steig([F[i] for i in k], [Rr[i] for i in k]))
    return float(steig(F, Rr)), float(np.std(bs))


@formel('rest_jenseits', 'ε_Q = Median über Galaxien des Medians von log₁₀(g_obs/(ν(g_N)g_N)) an Punkten a_N < a_Q; '
        'Fehler: Bootstrap über Galaxien ⊕ Monte Carlo (Abstand, Neigung, Υ, e_V)',
        ('ε_Q',), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'D_A', 'i_inc', 'Υ_d', 'σ_v', 'a_N', 'a₀'), 'L', 'modell',
        'Baustein 03', (_F03, 'jenseits folgen die SPARC-Kurven der Regel'))
def rest_jenseits(tafel, ups_d, sig, a_Q, a0=A0, n_boot=1000, n_mc=200, seed=29092026):
    """Rest der Regel an den Punkten mit a_N < a_Q (etwa jenseits R_Q, a_Q = (ln q_kosm)²a₀), je Galaxie der Median,
    dann der Median über die Galaxien mit solchen Punkten; Galaxien mit V_bar² > 0, V_obs > 0 an mindestens 5 Punkten,
    e_V ≥ 1 km/s. Gibt (Rest, Fehler); n_boot = 0 nur den Rest."""
    gal = []
    for g in tafel:
        Vb2 = v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], ups_d)
        m = (Vb2 > 0) & (g["Vo"] > 0)
        if m.sum() >= 5 and (Vb2[m] / g["R"][m] < a_Q).any():
            gal.append(dict(g, R=g["R"][m], Vo=g["Vo"][m], eV=np.maximum(g["eV"][m], 1e3), Vg=g["Vg"][m], Vd=g["Vd"][m],
                            Vb=g["Vb"][m]))

    def rest(g, fD=1.0, fi=1.0, fU=1.0, dV=0.0):
        R = g["R"] * fD
        Vb2 = fD * (g["Vg"] * np.abs(g["Vg"]) + fU * (ups_d * g["Vd"]**2 + _UPS_B * g["Vb"]**2))
        gN = np.clip(Vb2, 1e-30, None) / R
        r = np.log10(g_gemessen(R, (g["Vo"] + dV) * fi, g["Rd"] * fD, sig) / (nu(gN, a0) * gN))
        k = gN < a_Q
        return np.median(r[k]) if k.any() else np.nan
    rng = np.random.default_rng(seed)
    r0 = np.array([rest(g) for g in gal])
    rv = r0[np.isfinite(r0)]
    if not n_boot:
        return float(np.median(rv))
    boot = np.std([np.median(rv[rng.integers(0, rv.size, rv.size)]) for _ in range(int(n_boot))])
    mc = []
    for _ in range(int(n_mc)):
        z = []
        for g in gal:
            fD = max(1 + rng.normal() * g["eD"] / g["D"], 0.2)
            inc = min(max(g["inc"] + rng.normal() * g["eInc"], 20.0), 90.0)
            z.append(rest(g, fD, math.sin(math.radians(g["inc"])) / math.sin(math.radians(inc)),
                          10 ** (rng.normal() * SPARC_SATZ["ups_streuung_dex"]), rng.normal(size=g["Vo"].size) * g["eV"]))
        mc.append(np.nanmedian(z))
    return float(np.median(rv)), float(math.hypot(boot, np.std(mc)))


@formel('reststeigung_aussen', 'k_aus = Σ XY/Σ X², X, Y: log₁₀(a_N/a₀) und Rest je Galaxie zentriert, Punkte a_N < 0,1·a₀ '
        'der Galaxien mit äußerstem a_N ≤ a_rand; Regel mit äußerem Feld ν(|g_N + g_e|)',
        ('k_aus',), ('V_obs', 'V_gas', 'V_disk', 'V_bul', 'Υ_d', 'σ_v', 'a_e', 'a₀'), 'L', 'modell', 'Baustein 03',
        (_F03, 'Reststeigung'))
def reststeigung_aussen(tafel, ups_d, sig, g_e, lage, a_rand, a0=A0, n_boot=500, seed=29092026):
    """Steigung des Rests log₁₀(g_obs/(ν(|g_N + g_e|)g_N)) gegen log₁₀ a_N innerhalb der Galaxien an ihren Punkten
    unter 0,1·a₀ (mindestens 3), in den Galaxien, deren äußerster Punkt a_N ≤ a_rand hat. Gibt (Steigung, Streuung
    über n_boot Ziehungen ganzer Galaxien); n_boot = 0 nur die Steigung."""
    X, Y = [], []
    for g in tafel:
        Vb2 = v_bar_quadrat(g["Vg"], g["Vd"], g["Vb"], ups_d)
        m = (Vb2 > 0) & (g["Vo"] > 0)
        if m.sum() < 5:
            continue
        R = g["R"][m]
        gN = np.clip(Vb2[m], 1e-30, None) / R
        if gN[-1] > a_rand:
            continue
        n_ = np.array([nu_aussenfeld(x, g_e, lage, a0) for x in gN]) if g_e else nu(gN, a0)
        r = np.log10(g_gemessen(R, g["Vo"][m], g["Rd"], sig) / (n_ * gN))
        lx = np.log10(gN / a0)
        k = lx < -1
        if k.sum() >= 3:
            X.append(lx[k] - lx[k].mean())
            Y.append(r[k] - r[k].mean())
    steig = lambda j: sum(np.sum(X[i] * Y[i]) for i in j) / sum(np.sum(X[i]**2) for i in j)  # noqa: E731
    if not n_boot:
        return float(steig(range(len(X))))
    rng = np.random.default_rng(seed)
    bs = [steig(rng.integers(0, len(X), len(X))) for _ in range(int(n_boot))]
    return float(steig(range(len(X)))), float(np.std(bs))


@formel('sparc_ziehungen', 'Ziehungen der SPARC-Punkte: Galaxien mit Zurücklegen, Abstand, Neigung und Υ je Galaxie aus '
        'ihren Fehlern; g_N = V_bar²/R, g_obs mit Gasdruck', ('g_obs',), ('V_obs', 'D_A', 'i_inc', 'Υ_d', 'σ_v'), 'L',
        'hilfe', 'SPARC (Lelli, McGaugh, Schombert 2016)', alternative='g_gemessen')
def sparc_ziehungen(tafel, n, seed, ups_d, sig, boot=True, mc=True):
    """Liste von n Ziehungen (g_N, g_obs) aller Punkte mit V_bar² > 0 und V_obs > 0; n = 0 gibt die Nennwerte.
    Abstand D(1 + δ), δ ~ N(0, e_D/D), Faktor ≥ 0,2; Neigung i + N(0, e_i) in [20°, 90°]; Υ·10^{N(0; 0,1)}; eigener
    Generator mit seed."""
    rng = np.random.default_rng(seed)
    aus = []
    for k in range(max(int(n), 1)):
        idx = rng.integers(0, len(tafel), len(tafel)) if (boot and n) else np.arange(len(tafel))
        N, C = [], []
        for j in idx:
            g = tafel[j]
            fD = fi = fU = 1.0
            if mc and n:
                fD = max(1 + rng.normal() * g["eD"] / g["D"], 0.2)
                inc = min(max(g["inc"] + rng.normal() * g["eInc"], 20.0), 90.0)
                fi = math.sin(math.radians(g["inc"])) / math.sin(math.radians(inc))
                fU = 10 ** (rng.normal() * SPARC_SATZ["ups_streuung_dex"])
            R, Vo = g["R"] * fD, g["Vo"] * fi
            Vb2 = fD * (g["Vg"] * np.abs(g["Vg"]) + fU * (ups_d * g["Vd"]**2 + _UPS_B * g["Vb"]**2))
            m = (Vb2 > 0) & (Vo > 0)
            N.append(Vb2[m] / R[m])
            C.append(g_gemessen(R[m], Vo[m], g["Rd"] * fD, sig))
        aus.append((np.concatenate(N), np.concatenate(C)))
    return aus


@formel('ds_streuung', 'σ(d_s) = Streuung von d_s (ds_anpassung) über Ziehungen von Galaxien, Abstand, Neigung und Υ',
        ('σ_ds',), ('V_obs', 'D_A', 'i_inc', 'Υ_d', 'σ_v', 'a₀'), 'LI', 'modell', 'Baustein 07',
        (_F07, 'An den SPARC-Kurven gibt die Form der Regel d_s'))
def ds_streuung(tafel, ups_d, sig, n=400, seed=29092026, a0=A0):
    """Fehler der Spektraldimension aus n Ziehungen (sparc_ziehungen mit Bootstrap und Monte Carlo)."""
    return float(np.std([ds_anpassung(gN, gc, a0) for gN, gc in sparc_ziehungen(tafel, n, seed, ups_d, sig)]))


@formel('boden_anpassung', 'B_a: ln g_obs = ln(ν(g_N)·g_N + B_a·a₀) an den Punkten s < s_max, robust (soft_l1, 0,1)',
        ('B_a',), ('g_obs', 'g_N', 'a₀'), 'LI', 'modell', 'Baustein 07', (_F07, 'an SPARC liegt ein Boden unter'))
def boden_anpassung(gN, g, s_max, a0=A0):
    """Boden der Regel a = ν(a_N)·a_N + B (B in a₀, auch negativ zugelassen) an den tiefen Punkten."""
    from scipy.optimize import least_squares
    m = np.sqrt(gN / a0) < s_max
    x, y = gN[m], np.log(g[m])
    r = least_squares(lambda p: np.log(np.clip(a_regel(x, a0) + p[0] * a0, 1e-40, None)) - y, [0.0], loss="soft_l1",
                      f_scale=0.1)
    return float(r.x[0])


@formel('boden_grenze', 'B_95 = 95-%-Quantil von B_a über Ziehungen (mindestens 0)', ('B_a',),
        ('V_obs', 'D_A', 'i_inc', 'Υ_d', 'σ_v', 'a₀'), 'LI', 'modell', 'Baustein 07',
        (_F07, 'an SPARC liegt ein Boden unter'), alternative='boden_anpassung')
def boden_grenze(tafel, ups_d, sig, s_max, n=400, seed=29092026, a0=A0):
    """Obere Grenze (95 %) des Bodens aus n Ziehungen (sparc_ziehungen)."""
    B = [boden_anpassung(gN, gc, s_max, a0) for gN, gc in sparc_ziehungen(tafel, n, seed, ups_d, sig)]
    return max(float(np.percentile(B, 95)), 0.0)


@formel('boden_drift', 'B = 2Pe·a₀ (Drift zum Torus hin)', ('B_a',), ('Pe', 'a₀'), 'I', 'modell', 'Baustein 07',
        (_F07, 'Eine Drift zum Torus hin gäbe einen Boden a → 2Pe·a₀'), alternative='boden_anpassung')
def boden_drift(Pe, a0=A0):
    """Boden, den eine Drift der Rückkehr mit der Péclet-Zahl Pe gäbe (in m/s²)."""
    return 2 * Pe * a0


@formel('boden_quer', 'B = √(2/π)·(ϖ_R/b)·a₀ (kontrahierende Querrichtung der Bandbreite b)', ('B_a',),
        ('ϖ_R', 'b_q', 'a₀'), 'I', 'modell', 'Baustein 07', (_F07, 'eine kontrahierende Querrichtung'),
        alternative='boden_anpassung')
def boden_quer(varpi_R, b, a0=A0):
    """Boden, den eine kontrahierende Querrichtung der Bandbreite b gäbe (in m/s²)."""
    return math.sqrt(2 / math.pi) * varpi_R / b * a0


@formel('sigma_aus_rotation', 'Σ(R) = (1/π²G)[(1/R)∫₀^R (dV²/dr)K(r/R)dr + ∫_R^∞ (dV²/dr)K(R/r)dr/r] (Toomre 1963)',
        ('Σ',), ('v', 'r'), 'L', 'hilfe', 'Toomre 1963 (dünne Scheibe)', alternative='scheibe_quelle')
def sigma_aus_rotation(R, V2, n=4000):
    """Flächendichte einer dünnen Scheibe aus ihrem Rotationsbeitrag V² (mit Vorzeichen) an den Radien R (m, m²/s²);
    innen linear bis 0, außen keplersch bis 40·R_letzter. In kg/m², nicht negativ."""
    from scipy.special import ellipk
    r = np.linspace(1e-3 * R[-1] / 1e3, R[-1], n)
    yi = np.interp(r, np.r_[0, R], np.r_[0, V2])
    rr = np.r_[r, np.linspace(R[-1] * 1.001, 40 * R[-1], n)]
    yy = np.r_[yi, V2[-1] * R[-1] / rr[n:]]
    dy = np.gradient(yy, rr)
    aus = []
    for Rk in R:
        i_, o_ = rr < Rk * (1 - 1e-4), rr > Rk * (1 + 1e-4)
        v = (np.trapezoid(dy[i_] * ellipk(np.minimum((rr[i_] / Rk)**2, 1 - 1e-12)), rr[i_]) / Rk
             + np.trapezoid(dy[o_] * ellipk(np.minimum((Rk / rr[o_])**2, 1 - 1e-12)) / rr[o_], rr[o_]))
        aus.append(v / (math.pi**2 * G))
    return np.maximum(np.array(aus), 0.0)


_F_GITTER = np.round(np.arange(-0.5, 1.51, 0.02), 3)


@formel('gewicht_f_tafel', 'χ²(f, a₀) je Galaxie: ln g_obs gegen ln[ν(√(g_R² + |f|(2πGΣ)²))·g_R] (f < 0 gespiegelt), '
        'σ² = (2e_V/V)² + (e_D/D)² + (2e_i cot i)² + (∂ln a/∂ln Υ·σ_Υ)², Σ = Υ_d·SB_disk + Σ_gas',
        ('F_q',), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'SB₀', 'Υ_d', 'f_Σ', 'a₀'), 'L', 'hilfe', 'Baustein 07',
        alternative='fq_form')
def gewicht_f_tafel(tafel, ups_d):
    """Tafel (Galaxien × f × a₀) der Fehlerquadratsummen; f von −0,5 bis 1,5 (Schritt 0,02), a₀ von 0,7 bis 1,9·10⁻¹⁰
    (61 Punkte, logarithmisch). Punkte mit V_obs > 0, R > 0, e_V > 0; Gas aus V_gas (sigma_aus_rotation), Sterne aus der
    Flächenhelligkeit der Scheibe. Gibt (Tafel, Namen, a₀-Gitter); bleibt im Speicher."""
    schluessel = ("f", tuple(g["nm"] for g in tafel), ups_d)
    if schluessel in _PROFIL_SPEICHER:
        return _PROFIL_SPEICHER[schluessel]
    AG = np.exp(np.linspace(math.log(0.7e-10), math.log(1.9e-10), 61))
    su = math.log(10) * SPARC_SATZ["ups_streuung_dex"]
    msp = M_sun / pc**2
    T, namen = [], []
    for g in tafel:
        m = (g["Vo"] > 0) & (g["R"] > 0) & (g["eV"] > 0)
        if not m.any():
            continue
        R, Vo = g["R"][m], g["Vo"][m]
        Sg = sigma_aus_rotation(g["R"], g["Vg"] * np.abs(g["Vg"]))[m]

        def lna(Y, f, a0):
            gR = np.maximum((g["Vg"][m] * np.abs(g["Vg"][m]) + Y * g["Vd"][m] * np.abs(g["Vd"][m])
                             + _UPS_B * g["Vb"][m]**2) / R, 1e-14)
            gz = 2 * math.pi * G * (Y * g["SBd"][m] * msp + Sg)
            l0 = np.log(nu(gR, a0))
            l1 = np.log(nu(np.sqrt(gR[None, :]**2 + np.abs(f)[:, None] * gz[None, :]**2), a0))
            return np.log(gR)[None, :] + np.where(f[:, None] >= 0, l1, 2 * l0[None, :] - l1)
        d = np.log(Vo**2 / R)
        null = np.zeros(1)
        dY = (lna(ups_d * math.exp(su), null, A0) - lna(ups_d * math.exp(-su), null, A0))[0] / 2
        inc = math.radians(g["inc"])
        s2 = (2 * g["eV"][m] / Vo)**2 + (g["eD"] / g["D"])**2 + (2 * math.radians(g["eInc"]) / math.tan(inc))**2 + dY**2
        T.append(np.array([np.sum((d[None, :] - lna(ups_d, _F_GITTER, a0))**2 / s2[None, :], 1) for a0 in AG]).T)
        namen.append(g["nm"])
    _PROFIL_SPEICHER[schluessel] = (np.array(T), namen, AG)
    return _PROFIL_SPEICHER[schluessel]


@formel('gewicht_f', 'f̂ = argmin_f min_{a₀} Σ_Galaxien χ²(f, a₀) (Parabel im Gitter), Fehler: Jackknife über Galaxien',
        ('f_Σ',), ('V_obs', 'e_V', 'V_gas', 'V_disk', 'V_bul', 'SB₀', 'Υ_d'), 'L', 'modell', 'Baustein 07',
        (_F07, 'geht das senkrechte Scheibenfeld mit dem Gewicht'))
def gewicht_f(tafel, ups_d, ohne=(), mit_fehler=True):
    """Gewicht des senkrechten Scheibenfelds in a = ν(√(g_R² + f(2πGΣ)²))·g_R an allen SPARC-Kurven der Tafel, ohne die
    genannten Galaxien. Gibt (f̂, Jackknife-Fehler), ohne mit_fehler nur f̂."""
    T, namen, _ = gewicht_f_tafel(tafel, ups_d)
    sel = np.array([n not in ohne for n in namen])
    h = _F_GITTER[1] - _F_GITTER[0]

    def best(C):
        Cf = C.min(1)
        j = int(np.argmin(Cf))
        if 0 < j < len(Cf) - 1:
            y0, y1, y2 = Cf[j - 1:j + 2]
            return _F_GITTER[j] + h * (y0 - y2) / (2 * (y0 - 2 * y1 + y2))
        return float(_F_GITTER[j])
    S = T[sel].sum(0)
    if not mit_fehler:
        return float(best(S))
    jk = np.array([best(S - T[i]) for i in np.where(sel)[0]])
    n = jk.size
    return float(best(S)), float(math.sqrt((n - 1) / n * np.sum((jk - jk.mean())**2)))

