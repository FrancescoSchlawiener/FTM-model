# -*- coding: utf-8 -*-
"""
Kosmos (Ebene G, das Ganze, der Torus von außen bzw. die Zelle): Hintergrund E(z) (einziger Ort), eigener Hintergrund
H_b und Fluchtmaß, Abstände und Horizonte, Rückstellbilanz, Vakuumteil Λ_s, Ω_s, H_s, w_s, Θ = 2πa₀/(cH), Λ als ganze
innere Faltung der Zelle, Zellgrenze, a₀ früh, N_eff, Übergabezeit, Entropie S_H und S_Λ, Kreislauf.
Blickrichtung: Λ, Λ_s, H_s, S_Λ gehören zur inneren Richtung (seit z ≈ 10 fest); H, Θ, Ω_s, E(z), Horizonte und d_LSS
sind Zustand der Epoche und des Beobachters und haben die Epoche (z oder H) als Eingang.
"""
from .konstanten import *
from .formeln import formel
from .rueckkehr import gewicht_w, schliesslaenge
from .setzungen import SETZUNGEN


@formel('a0_modell', 'a_GH = cH/2π (Horizonttemperatur als Beschleunigung; heute nahe a₀)',
        ('a_GH',), ('H₀',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'a₀ ≈ c·kT_GH/ħ = cH₀/2π'))
def a0_modell(H0=H0_KMS):
    """cH/2π: Temperatur des Hubble-Horizonts als Beschleunigung. Fällt mit H; gleich a₀ nur heute (Θ ≈ 1)."""
    return c * H0_si(H0) / (2 * math.pi)


@formel('E_hintergrund', 'E = √(Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_X), a = 1/(1 + z) = e^{ln a}, Ω_X fest (ohne Angabe 1 − Ω_m − Ω_r)',
        ('E(z)',), ('z', 'ln a', 'Ω_m', 'Ω_r', 'Ω_X'), 'G', 'hilfe', 'Baustein 03', None)
def E_hintergrund(z, Om=OMEGA_M, Or=OMEGA_R, Ox=None, ln_a=None):
    """Einziger Ort von E = H/H₀ (Vergleichshintergrund mit festem Rest Ω_X, ohne Angabe flach 1 − Ω_m − Ω_r); z Zahl
    oder Feld. Die Epoche steht als z oder, für die Zukunft (z → −1 ist in Gleitkomma nicht darstellbar), als ln a
    (dann z = None)."""
    Ox = 1 - Om - Or if Ox is None else Ox
    la = -np.log1p(np.asarray(z, dtype=float)) if ln_a is None else np.asarray(ln_a, dtype=float)
    e = np.sqrt(Om * np.exp(-3 * la) + Or * np.exp(-4 * la) + Ox)
    return float(e) if e.ndim == 0 else e


@formel('weltzeit', 't_A(a) = ∫₀^a da′/(a′H) = (1/H₀)∫_{−∞}^{ln a} d ln a′/E(a′) (Anfang bis a = 10⁻⁸ analytisch)',
        ('t_A',), ('ln a', 'H₀', 'Ω_m', 'Ω_r', 'Ω_X', 'E(z)'), 'G', 'modell', 'Baustein 04',
        ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'nach rund 4·10¹¹ bis 10¹² Jahren'))
def weltzeit(ln_a, Om=OMEGA_M, Or=OMEGA_R, Ox=None, H0=H0_KMS, E=None):
    """Weltzeit bis zur Epoche ln a (Atemmaß des Ganzen) längs der Atemrate des Vergleichshintergrunds E (fester Rest
    Ω_X) oder einer gegebenen Atemrate E(ln a) = H/H₀ (etwa des eigenen Atmens, früh gleich dem Vergleich). Bis a = 10⁻⁸
    trägt allein die Strahlung, t = a²/(2H₀√Ω_r) (ohne Strahlung die Materie, t = ⅔a^{3/2}/(H₀√Ω_m)); danach
    Quadratur in ln a in Stücken von höchstens 4."""
    from scipy.integrate import quad
    u0 = math.log(1e-8)
    t = math.exp(2 * u0) / (2 * math.sqrt(Or)) if Or > 0 else 2 / 3 * math.exp(1.5 * u0) / math.sqrt(Om)
    f = (lambda u: 1 / E_hintergrund(None, Om, Or, Ox, ln_a=u)) if E is None else (lambda u: 1 / E(u))
    grenzen = np.append(np.arange(u0, ln_a, 4.0), ln_a)
    rel = 1e-11 if E is None else 1e-8          # eine gegebene Atemrate ist meist stückweise (Gitter): weicher
    t += sum(quad(f, a, b, limit=400, epsabs=0, epsrel=rel)[0] for a, b in zip(grenzen, grenzen[1:]))
    return t / H0_si(H0)


@formel('horizontwachstum', 'Ṡ_H = −2S_H·Ḣ/H, −Ḣ = H₀²(³⁄₂Ω_m a⁻³ + 2Ω_r a⁻⁴) (Verdünnung der Materie und Strahlung)',
        ('Ṡ_H',), ('ln a', 'H₀', 'Ω_m', 'Ω_r', 'Ω_X'), 'G', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Solange die Verdünnung der Materie den Horizont wachsen lässt'))
def horizontwachstum(ln_a=0.0, Om=OMEGA_M, Or=OMEGA_R, Ox=None, H0=H0_KMS):
    """Wachstum der Horizontentropie S_H = πc²/(ℓ_P²H²) in k_B/s zur Epoche ln a: Nur die Verdünnung (Materie,
    Strahlung) senkt H; der feste Rest Ω_X trägt zu Ḣ nichts bei."""
    H = H0_si(H0) * E_hintergrund(None, Om, Or, Ox, ln_a=ln_a)
    mHdot = H0_si(H0)**2 * (1.5 * Om * math.exp(-3 * ln_a) + 2 * Or * math.exp(-4 * ln_a))
    return 2 * horizontentropie(H) * mHdot / H


@formel('kreislauf_versiegen', 'ln a mit Ṡ_H(a) = Ṡ_i, d. h. −Ḣ = ℓ_P²H³Ṡ_i/(2πc²) (ab dort muss neue Faltung den Horizont tragen)',
        ('ln a',), ('Ṡ_i', 'H₀', 'Ω_m', 'Ω_r', 'Ω_X'), 'G', 'modell', 'Baustein 04',
        ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'ab etwa 2,6·10¹¹ Jahren verletzt'))
def kreislauf_versiegen(S_i, Om=OMEGA_M, Or=OMEGA_R, Ox=None, H0=H0_KMS):
    """Epoche (ln a, nach heute), ab der die Verdünnung den Horizont langsamer wachsen lässt, als Entropie entsteht
    (Ṡ_i in k_B/s je Horizont, fest gehalten). Mit festem Rest Ω_X fällt Ṡ_H ∝ a⁻³ und schneidet Ṡ_i genau einmal."""
    from scipy.optimize import brentq
    return brentq(lambda u: math.log(horizontwachstum(u, Om, Or, Ox, H0) / S_i), 0.0, 200.0, xtol=1e-12)


@formel('hubble_rate', 'H = H₀·E(z) im Parametersatz p (p.H0, p.z)',
        ('H',), ('H₀', 'z'), 'G', 'hilfe', 'Baustein 03', None, alternative='H0_si')
def hubble_rate(p):
    """Atemrate der Epoche des Parametersatzes: H₀ und z aus p, E(z) aus E_hintergrund."""
    return H0_si(p.H0) * E_hintergrund(p.z)


@formel('hintergrund_modell', 'H_b² = (8πG/3)(ρ_m(1+z)³ + ρ_r(1+z)⁴) + Λ_s c²/3, ρ = ω·3(100 km/s/Mpc)²/(8πG)',
        ('H_b',), ('z', 'ω_b', 'ω_c', 'Ω_r', 'a₀'), 'G', 'modell', 'Baustein 04',
        ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'Rate von Bremsglied der Materie und Vakuumteil π⁴/(15Λ²)'))
def hintergrund_modell(z, a0=A0):
    """Eigener Hintergrund ohne Fluchtmaß: fest sind die physikalischen Dichten (ω_b, ω_c, ω_r) und Λ_s = π⁴/(15Λ²);
    kein Anteil Ω wird gesetzt (Anteile sind Folge von H). a₀ seit z ≈ 10 fest; früher trägt der Vakuumteil nichts."""
    H100 = H0_si(100.0)
    om, orad = OMEGA_B_H2 + OMEGA_C_H2, OMEGA_R * h_**2
    zz = 1 + np.asarray(z, dtype=float)
    H = np.sqrt(H100**2 * (om * zz**3 + orad * zz**4) + lambda_vakuum(schliesslaenge(a0)) * c**2 / 3)
    return float(H) if H.ndim == 0 else H


@formel('atemrate_eff', 'H_eff = H_b + ⅓·dX_an/dt',
        ('H_eff',), ('H_b', 'Ẋ_an'), 'G', 'hilfe', 'Rechnung lauf_eigenes_atmen (Forschung, verzögerte Ankunft widerlegt)', None)
def atemrate_eff(H_b, dX_dt):
    """Atemrate mit dem über die Rückkehr angekommenen Fluchtmaß X_an (Atmen der Raumschleifen, nicht Λ)."""
    return H_b + dX_dt / 3


@formel('abstand_mitbewegt', 'd_C = (c/H₀)∫₀^z dz′/E(z′)',
        ('d_C',), ('z', 'H₀'), 'G', 'hilfe', 'Baustein 03', None)
def abstand_mitbewegt(z, H0=H0_KMS, n=20000):
    zs = np.linspace(0, z, n)
    return c / H0_si(H0) * np.trapezoid(1 / E_hintergrund(zs), zs)


@formel('d_LSS', 'd_LSS = 2(c/H₀)∫₀^{z} dz′/E(z′), z = 1090',
        ('d_LSS',), ('z',), 'G', 'hilfe', 'Baustein 03', None)
def d_LSS(z=1090.0):
    """Durchmesser der letzten Streufläche (mitbewegt) für den Beobachter heute."""
    zs = np.concatenate([np.linspace(0, 10, 20001), np.geomspace(10, z, 20001)[1:]])
    return 2 * c / H0_si() * np.trapezoid(1 / E_hintergrund(zs), zs)


@formel('abstand_zwischen', 'd_C = (c/H₀)∫_{z₁}^{z₂} dz′/E(z′) (mitbewegt, über u = ln(1 + z))',
        ('d_C',), ('z', 'Ω_m', 'Ω_r'), 'G', 'hilfe', 'Baustein 03', None, alternative='abstand_mitbewegt')
def abstand_zwischen(z1, z2, Om=OMEGA_M, Or=OMEGA_R):
    """Mitbewegter Abstand zwischen z1 und z2 (heutige Einheiten), quad in u = ln(1 + z), bei z = 10 und 10⁴ geteilt."""
    from scipy.integrate import quad
    f = lambda u: math.exp(u) / E_hintergrund(math.expm1(u), Om, Or)
    u1, u2 = math.log1p(z1), math.log1p(z2)
    teile = sorted({u1, u2} | {u for u in (math.log1p(10.0), math.log1p(1e4)) if u1 < u < u2})
    return c / H0_si() * sum(quad(f, a, b, limit=400, epsabs=0, epsrel=1e-12)[0] for a, b in zip(teile, teile[1:]))


@formel('d_LSS_beobachter', 'd_LSS(z_B) = 2(c/H₀)∫_{z_B}^{1090} dz′/E(z′) (mitbewegt, heutige Einheiten)',
        ('d_LSS',), ('z',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'bei z = 1 wäre L/d_LSS = 1,32'), alternative='d_LSS')
def d_LSS_beobachter(z_beob, z=1090.0):
    """Durchmesser der letzten Streufläche für einen Beobachter bei z_beob (mitbewegt): wächst mit dem Alter des
    Beobachters; die Zelle L ist mitbewegt fest."""
    return 2 * abstand_zwischen(z_beob, z)


@formel('teilchenhorizont', 'd_H = (c/H₀)∫_z^∞ dz′/E(z′) (mitbewegt)',
        ('d_H',), ('z', 'Ω_m', 'Ω_r'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Kugel des eigenen Teilchenhorizonts heute'))
def teilchenhorizont(z, Om=OMEGA_M, Or=OMEGA_R):
    """Mitbewegter Teilchenhorizont eines Beobachters bei z (heutige Einheiten); physikalisch d_H/(1 + z)."""
    return abstand_zwischen(z, math.expm1(40.0), Om, Or)


@formel('rueckstellbilanz', 'ΔL/R = ΣΔs_b/R = 4πGρ_b(z)R²/(9c²), R = d_H/(1 + z)',
        ('ΔL/R',), ('z', 'ρ_b'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'ΔL_Raum/R = ΣΔs/R'))
def rueckstellbilanz(z=0.0):
    """σ = ΣΔs_Zentren − ΔL_Raum = 0: Dehnung der Raumschleifen je Radius = Faltung der Baryonen je Radius in der
    Kugel des eigenen Teilchenhorizonts (Beobachter bei z). Mit ihren Rückläufen mal ν_kosm = 1 + Q."""
    R = teilchenhorizont(z) / (1 + z)
    return 4 * math.pi * G * rho_b_mitbewegt() * (1 + z)**3 * R**2 / (9 * c**2)


@formel('faltung_horizont', 'Δs_H = r_H/6 = c/(6H)',
        ('Δs_H',), ('H',), 'G', 'modell', 'Baustein 04', ('04-zeit-verschraenkung-syntropie-gedaechtnis.md', 'Δs_H = r_H/6'))
def faltung_horizont(H):
    """Globale Engstelle: Horizont r_H = c/H jedes Beobachters zur jeweiligen Epoche als Faltung gelesen."""
    return c / (6 * H)


@formel('theta_', 'Θ = 2πa₀/(cH)',
        ('Θ',), ('a₀', 'H'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Θ = 2πa₀/(cH)'))
def theta_(a0, H):
    return 2 * math.pi * a0 / (c * H)


@formel('omega_vakuum', 'Ω_s = (π²/45)·Θ²',
        ('Ω_s',), ('a₀', 'H₀'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Ω_s = (π²/45)·Θ²'))
def omega_vakuum(a0=A0, H0=H0_KMS):
    """Vakuumteil als Anteil, Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH); H0 ist die Atemrate der Epoche in km/s/Mpc."""
    return K_S * theta_(a0, H0_si(H0))**2


@formel('lambda_vakuum', 'Λ_s = π⁴/(15Λ²)',
        ('Λ_s',), ('Λ',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Λ_s = π⁴/(15·Λ²)'))
def lambda_vakuum(Lam):
    """Λ_s = π⁴/(15Λ²) = 4π⁴/(15ℓ₀²), ℓ₀ = 2Λ."""
    return math.pi**4 / (15 * Lam**2)


@formel('H_vakuum', 'H_s² = Λ_s·c²/3',
        ('H_s',), ('Λ',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'H_s² = Λ_s c²/3'))
def H_vakuum(Lam):
    """Atemrate des Vakuumteils allein: H_s² = Λ_s·c²/3 = π⁴c²/(45Λ²) (Endzustand bei festem Λ)."""
    return math.sqrt(lambda_vakuum(Lam) * c**2 / 3)


@formel('zustand_vakuum', '1 + w_s = −⅔·d ln a₀/d ln a',
        ('Δw_Z',), ('d ln a₀/d ln a',), 'IG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'w_s = −1 − (2/3)·d ln a₀/d ln a'))
def zustand_vakuum(dlna0_dlna):
    """Zustandsgleichung des Vakuumteils als Abstand von −1: Λ_s ∝ a₀², also 1 + w_s = −⅔·d ln a₀/d ln a (a₀ fest: 0).
    Direkt ausgegeben, weil w_s selbst in Gleitkomma bei −1 bleibt (1 + w_s ≈ 10⁻¹⁹ heute)."""
    return -2 / 3 * dlna0_dlna


@formel('gewicht_hintergrund', '1 − w = e^{−√(6π/Θ)} (θ̃ = 3H)',
        ('w',), ('Θ',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', '1 − w = e^{−√(6π/Θ)}'), alternative='gewicht_w')
def gewicht_hintergrund(Theta):
    """Gewicht der Schließungsregel im mitatmenden Hintergrund (θ̃ = 3H): 3Hτ_D = 6π/Θ."""
    return -math.expm1(-math.sqrt(6 * math.pi / Theta))


@formel('rho_m_mitbewegt', 'ρ_m = Ω_m·3H₀²/(8πG)',
        ('ρ_m',), ('Ω_m', 'H₀'), 'G', 'hilfe', 'Baustein 03', None)
def rho_m_mitbewegt():
    return OMEGA_M * 3 * H0_si()**2 / (8 * math.pi * G)


@formel('rho_b_mitbewegt', 'ρ_b = Ω_b·3H₀²/(8πG)',
        ('ρ_b',), ('Ω_b', 'H₀'), 'G', 'hilfe', 'Baustein 03', None)
def rho_b_mitbewegt():
    return OMEGA_B * 3 * H0_si()**2 / (8 * math.pi * G)


@formel('a0_aus_zelle', 'a₀ = 3c⁴/(2Gρ_mL³), ρ_m = ν_kosm·ρ_b',
        ('a₀',), ('L', 'ρ_m'), 'IG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'a₀ = c²/(2Λ) = 3c⁴/(2G·ν_kosm·ρ_b·L³)'), alternative='a0_aus_schliesslaenge')
def a0_aus_zelle(L):
    """Λ als ganze innere Faltung der Zelle: Λ = Gρ_mL³/(3c²), a₀ = 3c⁴/(2Gρ_mL³)."""
    return 3 * c**4 / (2 * G * rho_m_mitbewegt() * L**3)


@formel('zelle_aus_a0', 'L = (3c⁴/(2Gρ_m a₀))^{1/3}',
        ('L',), ('a₀', 'ρ_m'), 'IG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'verlangt L = 28,1–26,9 Gpc'))
def zelle_aus_a0(a0):
    return (3 * c**4 / (2 * G * rho_m_mitbewegt() * a0)) ** (1 / 3)


@formel('zelle_hubble', 'L·H/c = (8π²/(ν_kosm·Ω_b·Θ))^{1/3}',
        ('L·H/c',), ('Θ', 'Ω_m'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'L·H/c = (8π²/(ν_kosm·Ω_b·Θ))^{1/3}'))
def zelle_hubble(Theta, Om=OMEGA_M):
    """Zelle in Hubble-Längen der Epoche; Om = ν_kosm·Ω_b der Epoche."""
    return (8 * math.pi**2 / (Om * Theta)) ** (1 / 3)


@formel('zellgrenze', '(1 + Q)·a₀ ≤ f³·3c⁴/(2Gρ_b L_min³), L_min = 0,985·d_LSS',
        ('(1+Q)a₀',), ('f³', 'L_min', 'ρ_b'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', '(1 + Q)·a₀ ≤ 7,17·10⁻¹⁰ m/s²·f³'))
def zellgrenze(f3=1.0, L_min=None):
    """Obere Grenze von (1 + Q)·a₀ aus der Zelle ohne Kreispaare (L ≥ L_min), Zellform f³ (Würfel 1). L_min hängt am
    Alter des Beobachters (d_LSS wächst); ohne Angabe 0,985·d_LSS für den Beobachter heute."""
    L_min = KREISPAAR_GRENZE * d_LSS() if L_min is None else L_min
    return f3 * 3 * c**4 / (2 * G * rho_b_mitbewegt() * L_min**3)


@formel('a0_frueh_relativ', 'a₀(z)/a₀ = 1/(1 + Ω_r/(Ω_m a)), a = 1/(1 + z)',
        ('a₀(z)/a₀',), ('z', 'Ω_r', 'Ω_m'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'bei der Rekombination 0,76·a₀'))
def a0_frueh_relativ(z):
    """Strahlung erhöht die innere Faltung: Λ(a) = Λ_m(1 + Ω_r/(Ω_m a))."""
    a = 1 / (1 + z)
    return 1 / (1 + OMEGA_R / (OMEGA_M * a))


@formel('delta_N_eff', 'ΔN_eff = f_N/(1 − f_N)·(1 + 0,2271·N_eff)/0,2271',
        ('ΔN_eff',), ('f_N', 'N_eff'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'im Strahlungszeitalter wie ΔN_eff = 2,09'))
def delta_N_eff(f):
    """Fester Anteil f am Gesamtbudget im Strahlungszeitalter als ΔN_eff."""
    return f / (1 - f) * (1 + 0.2271 * N_EFF_STD) / 0.2271


@formel('uebergabezeit', 'ω_orb² = Δ_B·H²/2, 1/λ_m = N*/ω_orb; gibt (1/ω_orb, 1/λ_m min, max)',
        ('ω_orb', 'λ_m'), ('Δ_B', 'z', 'N*'), 'LG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'λ(w) = (1 − w)·ω_orb/N*'))
def uebergabezeit(dichte_faktor, z=0.0, N_stern=None):
    """1/λ_m = N*/ω_orb, ω_orb² = (4πG/3)·Δ·ρ_krit(z); gibt (1/ω_orb, 1/λ_m min, max) in s. N* (untere, obere Grenze)
    ist die Setzung NSTERN_MIN/NSTERN_MAX, wenn nicht übergeben."""
    if N_stern is None:
        N_stern = (SETZUNGEN["NSTERN_MIN"].wert, SETZUNGEN["NSTERN_MAX"].wert)
    Hz = H0_si() * E_hintergrund(z)
    w = math.sqrt(dichte_faktor / 2) * Hz          # (4πG/3)·Δ·3H²/(8πG) = ΔH²/2
    return 1 / w, N_stern[0] / w, N_stern[1] / w


@formel('horizontentropie', 'S_H = πc²/(ℓ_P²H²) (in k_B)',
        ('S_H',), ('H',), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'S_H = πkc²/(ℓ_P²H²)'))
def horizontentropie(H):
    """S_max = πkc²/(ℓ_P²H²) in Einheiten von k."""
    return math.pi * c**2 / (l_P**2 * H**2)


@formel('entropie_vakuum', 'S_Λ = S_H(H_s) = (45/π³)(Λ/ℓ_P)² (in k_B)',
        ('S_Λ',), ('Λ',), 'IG', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'S_Λ = (45/π³)·k·(Λ/ℓ_P)²'))
def entropie_vakuum(Lam):
    """Horizontentropie zum Vakuumteil allein: S_Λ = S_H bei H = H_s; hängt nur am Quadrat der Faltungslänge Λ."""
    return horizontentropie(H_vakuum(Lam))


@formel('kreislauf_mindestrate', 'Λ̇/Λ ≥ Ṡ_i/(2S_Λ); gibt (Λ̇/Λ, S_Λ, H_s)',
        ('Λ̇/Λ',), ('Λ', 'Ṡ_i'), 'G', 'modell', 'Baustein 03', ('03-von-der-galaxie-bis-zum-kosmos.md', 'Λ̇/Λ ≥ Ṡ_i/(2S_Λ)'))
def kreislauf_mindestrate(Lam, sigma_k=None):
    """Λ̇/Λ ≥ Ṡ_i/(2S_Λ) im Endzustand; Entropieproduktion Ṡ_i (sigma_k) in k/s je Horizont, ohne Angabe der Messwert
    ENTROPIEPRODUKTION."""
    if sigma_k is None:
        from .messwert import MESSWERTE
        sigma_k = MESSWERTE["ENTROPIEPRODUKTION"].wert
    S = entropie_vakuum(Lam)
    return sigma_k / (2 * S), S, H_vakuum(Lam)


@formel('omega_de', 'Ω_DE = 1 − Ω_m − Ω_r (flacher Raum, im Standardmodell gelesen)',
        ('Ω_DE',), ('Ω_m', 'Ω_r'), 'G', 'hilfe', 'Baustein 03', None)
def omega_de(Om=OMEGA_M, Or=OMEGA_R):
    """Anteil der Dunklen Energie, wenn das Ausdehnungsgesetz als Standardmodell gelesen wird (flach)."""
    return 1 - Om - Or


# ------------------------------------------------------------------ Zelle: Form, Drehungen, Winkelmoden
_S3 = math.sqrt(3)
GITTER_ZELLE = {"Würfel": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                "sechseckig": ((1, 0, 0), (0.5, _S3 / 2, 0), (0, 0, 1)),
                "raumzentriert": ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5)),
                "flächenzentriert": ((0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0))}
# Reine Verschiebungsgitter der Zelle (Zeilen = Gittervektoren), je Familie das beste Glied; sechseckig mit Höhe = Seite.


@formel('formfaktor_zelle', 'f³ = ℓ_min³/V (kürzeste geschlossene Schleife hoch drei durch das Volumen der Zelle)',
        ('f³',), ('Gitter',), 'Z', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Würfel f³ = 1, sechseckig 1,15'))
def formfaktor_zelle(gitter):
    """Formfaktor der Zelle ohne Kreispaare: Bei gleichem Volumen V ist die kürzeste Schleife ℓ_min = f·V^{1/3}; die
    Zellgrenze wächst mit f³. Eine Formgröße (Z), unabhängig von Epoche und Beobachter."""
    A = np.array(GITTER_ZELLE[gitter], dtype=float)
    n = np.array([(i, j, k) for i in range(-2, 3) for j in range(-2, 3) for k in range(-2, 3) if (i, j, k) != (0, 0, 0)])
    return float(np.min(np.linalg.norm(n @ A, axis=1)))**3 / abs(np.linalg.det(A))


@formel('drehwinkel_zelle', 'φ_D = 360°/n mit e^{2πi/n}·(ℤ + τℤ) = ℤ + τℤ (Drehsymmetrien des flachen Torus ℂ/(ℤ + τℤ))',
        ('φ_D',), ('τ_Z',), 'Z', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Das sind genau die Drehsymmetrien flacher Tori'))
def drehwinkel_zelle(tau, n_max=12):
    """Drehwinkel (Grad, absteigend), um die der flache Torus mit Modul τ in sich übergeht: e^{2πi/n} bildet 1 und τ
    wieder auf Gitterpunkte x + yτ (x, y ganz) ab. Allgemein nur 180°; τ = i dazu 90°, τ = e^{2πi/3} dazu 120° und 60°."""
    tau = complex(tau)

    def im_gitter(w):
        y = w.imag / tau.imag
        x = w.real - y * tau.real
        return abs(x - round(x)) < 1e-9 and abs(y - round(y)) < 1e-9
    aus = []
    for n in range(2, n_max + 1):
        zeta = complex(math.cos(2 * math.pi / n), math.sin(2 * math.pi / n))
        if im_gitter(zeta) and im_gitter(zeta * tau):
            aus.append(360.0 / n)
    return tuple(aus)


@formel('streuung_kosmisch', 'σ_ℓ = √(2/(2ℓ + 1)) (ein Himmel: 2ℓ + 1 Moden je Winkelordnung)',
        ('σ_ℓ',), ('ℓ_W',), 'G', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'weniger als seine kosmische Streuung'))
def streuung_kosmisch(ell):
    """Relative Streuung der Leistung einer Winkelordnung ℓ, die ein einzelner Beobachter an seinem Himmel liest."""
    return math.sqrt(2 / (2 * ell + 1))


@formel('anteil_isotrop', 'Anteil = n/(2ℓ + 1) (n von 2ℓ + 1 Formen einer Winkelordnung im unendlichen, isotropen Raum)',
        ('Anteil_ℓ',), ('ℓ_W',), 'G', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', '(unendlicher Raum 40 %)'))
def anteil_isotrop(n, ell):
    """Ohne Zelle hat keine Richtung Vorrang: Jede der 2ℓ + 1 Formen einer Winkelordnung trägt gleich viel Leistung;
    Vergleich für die Ausrichtung in der geschlossenen Zelle."""
    return n / (2 * ell + 1)


# ------------------------------------------------------------------ Faltung der Strahlung, Störungen bei der Rekombination
@formel('gang_schwelle_strahlung', 'd ln a₀/d ln a = x/(1 + x), x = Ω_r(1 + z)/Ω_m (Strahlung verliert Faltung durch Rotverschiebung)',
        ('d ln a₀/d ln a',), ('z', 'Ω_r', 'Ω_m'), 'IG', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Zerstrahlung senkt Λ'))
def gang_schwelle_strahlung(z=0.0):
    """Gang der Schwelle aus der Faltung der Strahlung: Λ(a) = Λ_m(1 + Ω_r/(Ω_m a)) fällt mit a, a₀ = c²/2Λ steigt;
    heute x ≈ Ω_r/Ω_m = 2,9·10⁻⁴, der Anteil verschwindet mit 1/a."""
    x = OMEGA_R * (1 + z) / OMEGA_M
    return x / (1 + x)


@formel('wachstum_uebergabe', 'D″ + (2 + d ln E/d ln a)D′ = ³⁄₂Ω_m(a)(1 − ε)D, dε/d ln a = (1 − w)√(Ω_m(a)/2)(1 − ε)/N*, '
        '1 − w = e^{−√(3H/θ_c(a))}; sofort: Quelle ³⁄₂Ω_m(a)(1 − (1 − w)Q/(1 + Q))D; gibt Δ(fσ₈)/fσ₈ = ΔD′/D′ bei z',
        ('Δ(fσ₈)/fσ₈',), ('z', 'N*', 'a₀', 'Ω_m', 'Q'), 'G', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Peaks der Hintergrundstrahlung und Strukturwachstum folgen'))
def wachstum_uebergabe(z=0.0, N_stern=1.0, a0=A0, Om=OMEGA_M, sofort=False):
    """Lineares Wachstum mit der Schließungsregel im mitatmenden Hintergrund (θ̃ = 3H): Der Anteil 1 − w der Mode atmet
    nicht mit und geht mit λ = (1 − w)ω_orb/N* (ω_orb = H√(Ω_m(a)/2)) in die Operatorform über, die im homogenen
    Hintergrund nicht klumpt (ε); sofort (obere Schranke) klumpt der Anteil (1 − w)·Q/(1 + Q) gar nicht. a₀(z) mit der
    Faltung der Strahlung (a0_frueh_relativ). Hintergrund flach mit Materie und festem Rest, ohne Strahlung; D = a ab
    a = 10⁻³ (ε ab 10⁻⁴). Gibt die relative Änderung von fσ₈ ∝ dD/d ln a gegen dasselbe Wachstum ohne Übergabe."""
    from scipy.integrate import solve_ivp
    Q = OMEGA_C_H2 / OMEGA_B_H2

    def rechte_seite(n, y, mit):
        a = math.exp(n)
        E2 = E_hintergrund(None, Om, 0.0, ln_a=n)**2
        Oma = Om * a**-3 / E2
        D, Dp, eps = y
        omw = 1 - gewicht_w(3 * H0_si() * math.sqrt(E2), a0 * a0_frueh_relativ(1 / a - 1)) if mit else 0.0
        if sofort:
            return [Dp, -(2 - 1.5 * Oma) * Dp + 1.5 * Oma * (1 - omw * Q / (1 + Q)) * D, 0.0]
        return [Dp, -(2 - 1.5 * Oma) * Dp + 1.5 * Oma * (1 - eps) * D, omw * math.sqrt(Oma / 2) * (1 - eps) / N_stern]

    def lauf(mit):
        e3 = 0.0
        if mit and not sofort:
            s0 = solve_ivp(lambda n, y: [0.0, 0.0, rechte_seite(n, [1.0, 1.0, y[2]], True)[2]],
                           (math.log(1e-4), math.log(1e-3)), [0.0, 0.0, 0.0], rtol=1e-10, atol=1e-14)
            e3 = s0.y[2][-1]
        s = solve_ivp(lambda n, y: rechte_seite(n, y, mit), (math.log(1e-3), -math.log1p(z)), [1e-3, 1e-3, e3],
                      rtol=1e-11, atol=1e-16)
        return s.y[1][-1]
    return lauf(True) / lauf(False) - 1


@formel('neigung_stoerung', 'a_N = (Φ/c²)·c²·(1 + z)/d (Neigung einer linearen Störung der mitbewegten Skala d zur Epoche z)',
        ('a_N',), ('Φ/c²', 'd_C', 'z'), 'G', 'modell', 'Baustein 03',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'Bei der Rekombination liegen lineare Störungen über der Schwelle'))
def neigung_stoerung(phi_c2, d, z):
    """Newton-Neigung einer linearen Störung mit Potentialtiefe Φ/c² über ihre physikalische Länge d/(1 + z); ob die
    Regel auf sie wirkt, entscheidet die Summe der Verluste (w), nicht dieser Wert allein."""
    return phi_c2 * c**2 * (1 + z) / d

