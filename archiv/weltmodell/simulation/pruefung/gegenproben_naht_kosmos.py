# -*- coding: utf-8 -*-
"""
Gegenproben der Module kosmos und naht_haufen (Felder in gegenprobe.py). Jeder unabhängige Weg ruft die geprüfte
Formel nicht auf (Wächter). Kennungen "<Formel>_<Weg>". Handwerte sind Zahlen, die auf anderem Weg gerechnet wurden
(von Hand oder mit anderem Verfahren), mit der Toleranz ihrer Stellen.
"""
import math

import kern as K

from .gegenprobe import GEGENPROBEN

_g = GEGENPROBEN.neu
A0 = K.A0
H0 = K.H0_si()
KMS = 1e3 / K.Mpc                                  # 1 km/s/Mpc in s⁻¹
LAM = K.schliesslaenge(A0)
MW = 6e10 * K.M_sun                                # Milchstraße als Punktmasse (Eingabe M_MILCHSTRASSE)


def _a0_horizont(H0_kms):
    """a₀ mit Θ = 1 bei H₀ (von Hand, ohne Kern-Formel)."""
    return K.c * H0_kms * 1e3 / K.Mpc / (2 * math.pi)


# ------------------------------------------------------------------ Vakuumteil, Entropie, Zustandsgleichung
_g("omega_vakuum_theta_eins", "omega_vakuum", "grenzfall",
   [(_a0_horizont(67.4), 67.4), (_a0_horizont(73.0), 73.0)], lambda a0, H0_: math.pi**2 / 45)
_g("lambda_vakuum_dichte_hand", "lambda_vakuum", "handwert", [(LAM,)],
   lambda L: 8 * math.pi * K.G * 2.085e-27 / K.c**2, toleranz=3e-4,
   grund="ρ_s = π³a₀²/(30Gc²) = 2,085·10⁻²⁷ kg/m³ bei a₀ = 1,1·10⁻¹⁰, von Hand auf vier Stellen")
_g("lambda_vakuum_dichte", "lambda_vakuum", "zweiter_weg",
   [(K.schliesslaenge(a),) for a in (1.0e-10, 1.1e-10, 1.2e-10)],
   lambda L: 8 * math.pi * K.G / K.c**2 * math.pi**3 * (K.c**2 / (2 * L))**2 / (30 * K.G * K.c**2))
_g("H_vakuum_H0_wurzel_omega", "H_vakuum", "zweiter_weg",
   [(K.schliesslaenge(a),) for a in (1.0e-10, 1.1e-10, 1.2e-10)],
   lambda L: H0 * math.sqrt(K.omega_vakuum(K.c**2 / (2 * L), K.H0_KMS)))
_g("H_vakuum_33", "H_vakuum", "handwert", [(LAM,)], lambda L: 33.3 * KMS, toleranz=2e-3,
   grund="H_s = π²c/(√45·Λ) = 33,3 km/s/Mpc bei Λ = 13,24 Gpc, von Hand auf drei Stellen")
_g("entropie_vakuum_hand", "entropie_vakuum", "zweiter_weg",
   [(K.schliesslaenge(a),) for a in (1.0e-10, 1.1e-10, 1.2e-10)], lambda L: 45 / math.pi**3 * (L / K.l_P)**2)
_g("horizontentropie_faltung", "horizontentropie", "zweiter_weg", [(H0,), (3 * H0,), (1e-10,)],
   lambda H: 36 * math.pi * (K.c / (6 * H) / K.l_P)**2)


def _w_aus_dichte(x):
    """1 + w = −⅓·d ln ρ_s/d ln a mit ρ_s ∝ Λ_s(Λ(a₀(a))), a₀ ∝ a^x (Erhaltung der Energie im Ausdehnungsgesetz)."""
    rho = lambda ln_a: K.lambda_vakuum(K.schliesslaenge(A0 * math.exp(x * ln_a)))
    h = 1e-3
    return -(math.log(rho(h)) - math.log(rho(-h))) / (2 * h) / 3


_g("zustand_vakuum_erhaltung", "zustand_vakuum", "zweiter_weg", [(0.0,), (0.1,), (-0.3,), (1.0,)], _w_aus_dichte,
   toleranz=1e-9)


# ------------------------------------------------------------------ Hintergrund, Atmen, Horizonte
_g("hintergrund_modell_de_sitter", "hintergrund_modell", "grenzfall",
   [(-1 + 1e-6, 1.1e-10), (-1 + 1e-6, 1.06e-10)], lambda z, a0: K.H_vakuum(K.schliesslaenge(a0)))
_g("hintergrund_modell_hand", "hintergrund_modell", "handwert", [(0.0, 1.1e-10), (0.0, 1.06e-10)],
   lambda z, a0: {1.1e-10: 50.34, 1.06e-10: 49.55}[a0] * KMS, toleranz=2e-4,
   grund="H_b = 100·√(ω_m + ω_r + ω_s) km/s/Mpc, im Lauf des eigenen Atmens auf vier Stellen (50,34 bzw. 49,55)")
_g("atemrate_eff_ohne_flucht", "atemrate_eff", "grenzfall", [(H0, 0.0), (0.5 * H0, 0.0)], lambda Hb, X: Hb)
_g("atemrate_eff_hand", "atemrate_eff", "handwert", [(H0, 3 * H0)], lambda Hb, X: 2 * H0)
_g("d_LSS_beobachter_gitter", "d_LSS_beobachter", "zweiter_weg", [(0.0,)], lambda z: K.d_LSS(), toleranz=1e-7,
   grund="Trapez auf festem Gitter (d_LSS) gegen quad in ln(1 + z)")
_g("teilchenhorizont_eds", "teilchenhorizont", "grenzfall", [(0.0, 1.0, 0.0)], lambda z, Om, Or: 2 * K.c / H0,
   toleranz=5e-9, grund="Einstein-de Sitter 2c/H₀; das Integral endet bei 1 + z = e⁴⁰")


def _bilanz_faltung(z):
    """ΣΔs_b/R = Δs(M_b)/R mit M_b = (4π/3)ρ_b(z)R³ in der Kugel des eigenen Teilchenhorizonts."""
    R = K.teilchenhorizont(z) / (1 + z)
    return K.ds(4 * math.pi / 3 * K.rho_b_mitbewegt() * (1 + z)**3 * R**3) / R


_g("rueckstellbilanz_faltung", "rueckstellbilanz", "zweiter_weg", [(0.0,), (1.0,), (10.0,)], _bilanz_faltung)
_g("faltung_horizont_hubble_masse", "faltung_horizont", "zweiter_weg", [(H0,), (20 * H0,)],
   lambda H: K.ds(4 * math.pi / 3 * (K.c / H)**3 * 3 * H**2 / (8 * math.pi * K.G)))


# ------------------------------------------------------------------ Θ, Gewicht, Zelle
_g("theta_temperaturen", "theta_", "zweiter_weg", [(A0, H0), (1.2e-10, 2 * H0)],
   lambda a0, H: (K.hbar * a0 / (K.c * K.kB)) / (K.hbar * H / (2 * math.pi * K.kB)))
_g("a0_modell_theta", "a0_modell", "zweiter_weg", [(67.4,), (73.0,)],
   lambda H0_: A0 / K.theta_(A0, K.H0_si(H0_)))
_g("gewicht_hintergrund_rueckkehr", "gewicht_hintergrund", "zweiter_weg", [(1.0,), (K.theta_(A0, H0),), (0.3,)],
   lambda Th: K.gewicht_w(3 * 2 * math.pi * A0 / (K.c * Th), A0), toleranz=1e-12)
_g("zelle_hubble_zelle", "zelle_hubble", "zweiter_weg",
   [(K.theta_(a, H0), K.OMEGA_M) for a in (1.04e-10, 1.1e-10, 1.19e-10)],
   lambda Th, Om: K.zelle_aus_a0(Th * K.c * H0 / (2 * math.pi)) * H0 / K.c)
_g("a0_aus_zelle_faltung", "a0_aus_zelle", "zweiter_weg", [(27.6 * K.Gpc,), (20 * K.Gpc,)],
   lambda L: K.a0_aus_schliesslaenge((1 + K.Q_KOSM) * K.ds(K.rho_b_mitbewegt() * L**3)))
_g("zelle_aus_a0_umkehr", "zelle_aus_a0", "umkehr", [(1.04e-10,), (1.1e-10,), (1.19e-10,)],
   lambda a0: K.Gpc * (K.a0_aus_zelle(K.Gpc) / a0)**(1 / 3))
_L_MIN = K.KREISPAAR_GRENZE * 27.786 * K.Gpc
_g("zellgrenze_a0_zelle", "zellgrenze", "zweiter_weg", [(1.0, _L_MIN), (1.15, _L_MIN), (1.41, 30 * K.Gpc)],
   lambda f3, L: f3 * (1 + K.Q_KOSM) * K.a0_aus_zelle(L))
_g("zellgrenze_hand", "zellgrenze", "handwert", [(1.0, _L_MIN)], lambda f3, L: 7.17e-10, toleranz=7e-4,
   grund="(1 + Q)a₀ ≤ 7,17·10⁻¹⁰ m/s² (Würfel), drei Stellen")
_g("a0_frueh_relativ_hand", "a0_frueh_relativ", "handwert", [(1090.0,), (10.0,)],
   lambda z: {1090.0: 0.7586, 10.0: 0.99686}[z], toleranz=2e-4,
   grund="1/(1 + Ω_r(1 + z)/Ω_m) von Hand: Ω_r/Ω_m = 2,916·10⁻⁴")
_g("delta_N_eff_hand", "delta_N_eff", "handwert", [(math.pi**2 / 45,)], lambda f: 2.09, toleranz=3e-3,
   grund="ΔN_eff = 2,09 bei Ω_s = π²/45, drei Stellen")


# ------------------------------------------------------------------ Naht und Haufen
_g("r_Q_milchstrasse", "r_Q", "handwert", [(MW,)], lambda M: 51.0 * K.kpc, toleranz=1e-3,
   grund="R_Q(Milchstraße) = 8,727 kpc/0,17107 = 51,0 kpc, von Hand auf drei Stellen")


def _masse_bQ_naht(a0, H, Qk):
    """M_b,Q aus R_Q(M_b) = R₂₀₀ mit M_ges = (1 + Q)M_b (Nullstelle)."""
    from scipy.optimize import brentq
    rho_k = 3 * H**2 / (8 * math.pi * K.G)
    r200 = lambda Mb: (3 * (1 + Qk) * Mb / (4 * math.pi * 200 * rho_k))**(1 / 3)
    lm = brentq(lambda x: math.log(K.r_Q(10**x, a0, Qk / (1 + Qk)) / r200(10**x)), 38, 50, xtol=1e-14, rtol=1e-15)
    return 10**lm


_g("masse_bQ_naht", "masse_bQ", "zweiter_weg", [(1.1e-10, H0, K.Q_KOSM), (1.2e-10, 2 * H0, 5.0)], _masse_bQ_naht,
   toleranz=1e-11, grund="Nullstelle in log₁₀ M auf 10⁻¹⁴")
_g("dichte_umkehr_heute_eds", "dichte_umkehr_heute", "grenzfall", [(1.0,)], lambda Om: 9 * math.pi**2 / 16,
   toleranz=1e-8, grund="Einstein-de Sitter 9π²/16; Start a = 10⁻³ im Wachstumsmodus, rtol 10⁻⁹")
_g("umkehrradius_heute_trichterkette", "umkehrradius_heute", "zweiter_weg",
   [((1 + K.Q_KOSM) * m * K.M_sun,) for m in (1e10, 6e10, 1e12)],
   lambda M: (3 * M / (4 * math.pi * K.rho_m_mitbewegt() * 11.816))**(1 / 3), toleranz=3e-5,
   grund="1 + Δ_ta = 11,816 aus eigenem Top-Hat mit Strahlung (Start a = 10⁻², Zeit in 1/H₀), fünf Stellen")


def _zykloide_dauer(z_c, a0=A0, w_oben=0.9, w_unten=0.1):
    """Dauer w_oben → w_unten vor der Umkehr in der Zykloide (Einstein-de Sitter): R ∝ 1 − cos η, t = B(η − sin η),
    Kollaps t_c = 2πB = Alter bei z_c (mit Strahlung), θ = 3 sin η/(B(1 − cos η)²)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    t_c = quad(lambda z: 1 / ((1 + z) * H0 * K.E_hintergrund(z)), z_c, math.inf, limit=400)[0]
    B = t_c / (2 * math.pi)
    theta = lambda eta: 3 * math.sin(eta) / (B * (1 - math.cos(eta))**2)
    t = [B * (e - math.sin(e)) for e in (brentq(lambda eta: theta(eta) - math.log(1 - w)**2 * a0 / K.c, 1e-6, math.pi)
                                          for w in (w_oben, w_unten))]
    return t[1] - t[0]


_g("schalen_dauer_zykloide", "schalen_dauer", "grenzfall", [(3.4,), (2.0,)], _zykloide_dauer, toleranz=7e-3,
   grund="Zykloide ohne Vakuumteil gegen Energieintegral mit Ω_Λ: 0,5 % Unterschied bei z_c ≥ 2")
_g("fluchtmass_bereich_wirkung", "fluchtmass_bereich", "handwert", [(MW,)], lambda M: math.log(284.10),
   toleranz=1e-4, grund="Maß über das Wirkungsintegral J_r (anderer Weg): 284,10 für die Milchstraße")


# ------------------------------------------------------------------ Formeln mit mehreren Werten
def _mindestrate(Lam, sigma_k):
    """S_Λ = (45/π³)(Λ/ℓ_P)², H_s = π²c/(√45·Λ) aus H_s² = Λ_s c²/3, Λ_s = π⁴/(15Λ²); Rate Ṡ_i/(2S_Λ)."""
    S = 45 / math.pi**3 * (Lam / K.l_P)**2
    return sigma_k / (2 * S), S, math.pi**2 * K.c / (math.sqrt(45) * Lam)


_g("kreislauf_mindestrate_hand", "kreislauf_mindestrate", "zweiter_weg", [(LAM, 7.1e86), (2 * LAM, 1e80)], _mindestrate,
   toleranz=1e-12, grund="Rundung der geschlossenen Formen")


def _fenster(a0, L_min):
    """Q_min aus √(6/(1 + Q)) = 1 (Nullstelle); Q_max aus (1 + Q)a₀ = 3c⁴/(2Gρ_bL_min³), ρ_b = Ω_b·3H₀²/(8πG)."""
    from scipy.optimize import brentq
    rho_b = K.OMEGA_B * 3 * (K.H0_KMS * KMS)**2 / (8 * math.pi * K.G)
    return (brentq(lambda Q: math.sqrt(6 / (1 + Q)) - 1, 0, 100, xtol=1e-14),
            3 * K.c**4 / (2 * K.G * rho_b * L_min**3 * a0) - 1)


_g("q_fenster_hand", "q_fenster", "zweiter_weg", [(A0, 0.985 * 28e9 * K.pc), (1.2e-10, 30 * K.Gpc)], _fenster,
   toleranz=1e-12, grund="Nullstelle mit brentq (xtol 10⁻¹⁴)")


def _uebergabe(Delta, z, N_stern):
    """ω_orb² = (4πG/3)·Δ·ρ_krit(z), ρ_krit = 3H(z)²/(8πG); 1/λ_m = N*/ω_orb."""
    H = K.H0_KMS * KMS * K.E_hintergrund(z)
    w = math.sqrt(4 * math.pi * K.G / 3 * Delta * 3 * H**2 / (8 * math.pi * K.G))
    return 1 / w, N_stern[0] / w, N_stern[1] / w


_g("uebergabezeit_dichte", "uebergabezeit", "zweiter_weg", [(200.0, 0.0, (0.3, 2.4)), (18 * math.pi**2, 1.0, (1.0, 3.0))],
   _uebergabe, toleranz=1e-12, grund="Rundung (4πG/3)·3/(8πG) = 1/2")


# ------------------------------------------------------------------ Einsetzen als Kernformel (03)
_g("haufenstoss_trennzeit_umlauf", "haufenstoss_trennzeit", "grenzfall", [(2 * math.pi, 1e-16), (2 * math.pi, 3e-15)],
   lambda I, w: 2 * math.pi / w, grund=None)
_g("omega_de_hand", "omega_de", "handwert", [(0.3, 0.0), (0.3134, 9.14e-5)],
   lambda Om, Or: {(0.3, 0.0): 0.7, (0.3134, 9.14e-5): 0.6865086}[(Om, Or)], toleranz=1e-7,
   grund="1 − 0,3 = 0,7 und 1 − 0,3134 − 0,0000914 = 0,6865086 von Hand")


# ------------------------------------------------------------------ Weltzeit, Horizontwachstum, Kreislauf
def _zeit_materie_lambda(ln_a, Om, Or, Ox, H0_kms):
    """Materie und fester Rest ohne Strahlung: t = 2/(3H₀√Ω_X)·arsinh(√(Ω_X/Ω_m)·a^{3/2})."""
    return 2 / (3 * H0_kms * KMS * math.sqrt(Ox)) * math.asinh(math.sqrt(Ox / Om) * math.exp(1.5 * ln_a))


def _zeit_materie_strahlung(ln_a, Om, Or, Ox, H0_kms):
    """Materie und Strahlung ohne Rest: t = 2/(3H₀Ω_m²)·[(Ω_m a − 2Ω_r)√(Ω_m a + Ω_r) + 2Ω_r^{3/2}]."""
    a = math.exp(ln_a)
    return 2 / (3 * H0_kms * KMS * Om**2) * ((Om * a - 2 * Or) * math.sqrt(Om * a + Or) + 2 * Or**1.5)


_g("weltzeit_materie_lambda", "weltzeit", "grenzfall",
   [(0.0, 0.3134, 0.0, 0.6866, 67.4), (math.log(1e6), 0.3134, 0.0, 0.6866, 67.4), (math.log(0.1), 0.3, 0.0, 0.7, 70.0)],
   _zeit_materie_lambda, toleranz=1e-9, grund="Quadratur in ln a (epsrel 10⁻¹¹) gegen die geschlossene Form")
_g("weltzeit_materie_strahlung", "weltzeit", "grenzfall", [(0.0, 0.3, 1e-4, 0.0, 67.4), (math.log(0.01), 0.3, 1e-4, 0.0, 67.4)],
   _zeit_materie_strahlung, toleranz=1e-8, grund="Anfang bis a = 10⁻⁸ als reine Strahlung (Fehler Ω_m·10⁻⁸/Ω_r)")


def _E_materie_lambda(u):
    """Gegebene Atemrate E(ln a) für Materie und festen Rest (0,3134; 0,6866), ohne Strahlung."""
    return math.sqrt(0.3134 * math.exp(-3 * u) + 0.6866)


_g("weltzeit_gegebene_atemrate", "weltzeit", "grenzfall",
   [(0.0, 0.3134, 0.0, 0.6866, 67.4, _E_materie_lambda), (math.log(3.0), 0.3134, 0.0, 0.6866, 67.4, _E_materie_lambda)],
   lambda ln_a, Om, Or, Ox, H0_kms, E: _zeit_materie_lambda(ln_a, Om, Or, Ox, H0_kms), toleranz=1e-8,
   grund="gegebene Atemrate mit Quadratur epsrel 10⁻⁸ gegen die geschlossene Form")


def _horizont_ableitung(ln_a, Om, Or, Ox, H0_kms):
    """Ṡ_H = dS_H/d ln a · H mit S_H = πc²/(ℓ_P²H²) und H(ln a) von Hand; Ableitung als Differenz (h = 10⁻⁴)."""
    H = lambda u: H0_kms * KMS * math.sqrt(Om * math.exp(-3 * u) + Or * math.exp(-4 * u) + Ox)
    S = lambda u: math.pi * K.c**2 / (K.l_P**2 * H(u)**2)
    h = 1e-4
    return (S(ln_a + h) - S(ln_a - h)) / (2 * h) * H(ln_a)


_g("horizontwachstum_ableitung", "horizontwachstum", "zweiter_weg",
   [(0.0, K.OMEGA_M, K.OMEGA_R, 1 - K.OMEGA_M - K.OMEGA_R, 67.4), (2.0, 0.3, 1e-4, 0.24, 67.4), (-3.0, 0.3, 1e-4, 0.7, 70.0)],
   _horizont_ableitung, toleranz=1e-7, grund="zentrale Differenz, Fehler von der Ordnung h² = 10⁻⁸")


def _versiegen_fern(S_i, Om, Or, Ox, H0_kms):
    """Fern von heute trägt allein die Materie −Ḣ = ³⁄₂Ω_m a⁻³H₀² und H = H₀√Ω_X: a³ = 3πc²Ω_m/(ℓ_P²H₀Ω_X^{3/2}Ṡ_i)."""
    H0_ = H0_kms * KMS
    return math.log(3 * math.pi * K.c**2 * Om / (K.l_P**2 * H0_ * Ox**1.5 * S_i)) / 3


_g("kreislauf_versiegen_fern", "kreislauf_versiegen", "grenzfall",
   [(7.1e86, K.OMEGA_M, 0.0, 1 - K.OMEGA_M, 67.4), (1e85, 0.3, 0.0, 0.24, 67.4)], _versiegen_fern, toleranz=1e-10,
   grund="H = H₀√Ω_X genähert (Fehler Ω_m a⁻³/Ω_X ≈ 10⁻¹⁸), Nullstelle auf 10⁻¹² in ln a")


# ------------------------------------------------------------------ Zelle: Form, Drehungen, Winkelmoden
_g("formfaktor_zelle_hand", "formfaktor_zelle", "handwert",
   [("Würfel",), ("sechseckig",), ("raumzentriert",), ("flächenzentriert",)],
   lambda g: {"Würfel": 1.0, "sechseckig": 2 / math.sqrt(3), "raumzentriert": 3 * math.sqrt(3) / 4,
              "flächenzentriert": math.sqrt(2)}[g], toleranz=1e-12,
   grund="geschlossene Formen von Hand: Seite³/Volumen, 2/√3, (√3/2)³·2, (√2/2)³·4")


def _drehwinkel_bilder(tau, n_max=12):
    """Zahl c der nächsten Bilder (kürzeste Gittervektoren m + nτ) ist die größte Drehordnung; Drehungen 360°/n für n | c."""
    tau = complex(tau)
    laengen = [abs(m + n * tau) for m in range(-3, 4) for n in range(-3, 4) if (m, n) != (0, 0)]
    kurz = min(laengen)
    zahl = sum(1 for x in laengen if abs(x - kurz) < 1e-9 * kurz)
    return tuple(360.0 / n for n in range(2, n_max + 1) if zahl % n == 0)


_g("drehwinkel_zelle_bilder", "drehwinkel_zelle", "zweiter_weg",
   [(1j,), (complex(-0.5, math.sqrt(3) / 2),), (complex(0.3, 1.7),), (2j,)], _drehwinkel_bilder, toleranz=1e-12,
   grund="Zählung der nächsten Bilder statt Abbildung des Gitters (Rundung der Winkel)")


def _streuung_chi2(ell):
    """Leistung einer Winkelordnung an einem Himmel ist χ²-verteilt mit 2ℓ + 1 Freiheitsgraden: σ/Mittel."""
    from scipy.stats import chi2
    k = 2 * ell + 1
    return chi2(k).std() / k


_g("streuung_kosmisch_chi2", "streuung_kosmisch", "zweiter_weg", [(2,), (3,), (10,)], _streuung_chi2, toleranz=1e-12,
   grund="Rundung der Verteilungsfunktion")
_g("anteil_isotrop_hand", "anteil_isotrop", "handwert", [(2, 2), (1, 2), (1, 3)],
   lambda n, ell: {(2, 2): 0.4, (1, 2): 0.2, (1, 3): 1 / 7}[(n, ell)], toleranz=1e-12,
   grund="2/5, 1/5 und 1/7 von Hand")


# ------------------------------------------------------------------ Schale außerhalb der Naht, Übergabe vor der Umkehr
_D_TA = 11.8157                                     # 1 + Δ_ta aus der Umkehr heute (Rechnung gewicht_rueckkehr, 5 Stellen)
_g("gewicht_schale_aussen_hintergrund", "gewicht_schale_aussen", "grenzfall", [(_D_TA**(1 / 3),)],
   lambda x: K.gewicht_hintergrund(K.theta_(A0, H0)), toleranz=2e-5,
   grund="Schale mit der mittleren Dichte des Hintergrunds (1 + Δ = 1) dehnt mit θ = 3H₀; 1 + Δ_ta auf fünf Stellen")
_g("gewicht_schale_aussen_tafel", "gewicht_schale_aussen", "handwert",
   [((_D_TA / 5)**(1 / 3),), ((_D_TA / 8)**(1 / 3),)], lambda x: {5: 0.95240, 8: 0.88758}[round(_D_TA / x**3)],
   toleranz=2e-5, grund="Tafel 1 + Δ = 5 bzw. 8 der Rechnung gewicht_rueckkehr (w_R, fünf Stellen)")
_g("rand_schale_aussen_tafel", "rand_schale_aussen", "handwert", [(0.9,)], lambda w: 1.157, toleranz=5e-4,
   grund="Rechnung gewicht_rueckkehr: Schale w 0,9 → 0,1 in R/R_ta 1,157 → 1,000 (drei Stellen)")
_g("uebergabe_vor_umkehr_bahn", "uebergabe_vor_umkehr", "handwert", [(A0,)], lambda a0: 0.1098, toleranz=5e-4,
   grund="Bahnrechnung naht_rq_umkehr (Top-Hat in der Zeit, Milchstraße): I_vor = 0,1098")


# ------------------------------------------------------------------ Übergabe bei der ersten Faltung
_g("uebergabe_faltung_archiv", "uebergabe_faltung", "handwert", [(1.0, 0.2)], lambda e, ej: (1.2553, 1.870, 0.572),
   toleranz=2e-3, grund="anderes Verfahren (Rechnung uebergabe_faltung: 6000 Zeitpunkte, Faltung am Gitter, Mischung "
   "½), drei bis vier Stellen")


# ------------------------------------------------------------------ Wachstum mit der Schließungsregel
def _wachstum_linear_antwort(z, N_stern):
    """Lineare Antwort statt voller Kopplung: D₀ ohne Übergabe, ε′ = (1 − w)√(Ω_m(a)/2)/N* (ohne den Faktor 1 − ε),
    δD″ + (2 − ³⁄₂Ω_m(a))δD′ = ³⁄₂Ω_m(a)(δD − εD₀); Δ(fσ₈)/fσ₈ = δD′/D₀′. Gilt für ε ≪ 1 (große N*)."""
    from scipy.integrate import solve_ivp
    Om, H0_ = K.OMEGA_M, K.H0_KMS * KMS

    def rechts(n, y):
        a = math.exp(n)
        E2 = Om * a**-3 + 1 - Om
        Oma = Om * a**-3 / E2
        omw = math.exp(-math.sqrt(3 * H0_ * math.sqrt(E2) * K.c / (A0 * K.a0_frueh_relativ(1 / a - 1))))
        D0, D0p, dD, dDp, eps = y
        return [D0p, -(2 - 1.5 * Oma) * D0p + 1.5 * Oma * D0, dDp, -(2 - 1.5 * Oma) * dDp + 1.5 * Oma * (dD - eps * D0),
                omw * math.sqrt(Oma / 2) / N_stern]
    s = solve_ivp(rechts, (math.log(1e-3), -math.log1p(z)), [1e-3, 1e-3, 0.0, 0.0, 0.0], rtol=1e-11, atol=1e-18)
    return s.y[3][-1] / s.y[1][-1]


_g("wachstum_uebergabe_lineare_antwort", "wachstum_uebergabe", "grenzfall", [(0.0, 20.0), (0.5, 20.0)],
   _wachstum_linear_antwort, toleranz=5e-4, grund="lineare Antwort, Fehler von der Ordnung ε ≈ 2·10⁻⁴ (N* = 20); "
   "ohne den Anteil vor a = 10⁻³")
_g("wachstum_uebergabe_archiv", "wachstum_uebergabe", "handwert", [(0.0, 1.0), (0.0, 0.3), (0.0, 2.4)],
   lambda z, N: {1.0: -1.18e-3, 0.3: -3.94e-3, 2.4: -4.94e-4}[N], toleranz=5e-3,
   grund="anderes Verfahren (Rechnung gewicht_rueckkehr, w_R, Lesart Übergabe), drei Stellen")


# ------------------------------------------------------------------ Faltung der Strahlung, Störungen
def _gang_differenz(z):
    """d ln a₀/d ln a als Differenz von ln a0_frueh_relativ in ln a (h = 10⁻⁴)."""
    f = lambda u: math.log(K.a0_frueh_relativ(math.exp(-u) - 1))
    u, h = -math.log1p(z), 1e-4
    return (f(u + h) - f(u - h)) / (2 * h)


_g("gang_schwelle_strahlung_differenz", "gang_schwelle_strahlung", "zweiter_weg", [(0.0,), (10.0,), (1090.0,)],
   _gang_differenz, toleranz=1e-7, grund="zentrale Differenz, Fehler von der Ordnung h² = 10⁻⁸")
_g("neigung_stoerung_hand", "neigung_stoerung", "handwert", [(1e-5, 147 * K.Mpc, 1090.0)], lambda P, d, z: 2.1617e-10,
   toleranz=5e-5, grund="10⁻⁵·8,98755·10¹⁶·1091/(147·3,08568·10²² m) = 2,1617·10⁻¹⁰ m/s² von Hand")


# ------------------------------------------------------------------ Linse jenseits der Naht, Rest im Haufen
def _linse_plummer(M, R, b_rel=1e-6):
    """Plummer-Linse ΔΣ = Σ̄(<R) − Σ(R) = M R²/(π(R² + b²)²) im Grenzfall b → 0."""
    b = b_rel * R
    return M * R**2 / (math.pi * (R**2 + b**2)**2)


_g("linse_punktmasse_plummer", "linse_punktmasse", "grenzfall",
   [(MW * (1 + K.Q_KOSM), K.Mpc), (1e15 * K.M_sun, 0.3 * K.Mpc)], _linse_plummer, toleranz=1e-9,
   grund="Plummer mit b = 10⁻⁶R: Abstand 2·10⁻¹²")


def _rest_bose(r, M_b, M_L, a0):
    """ν = 1/(1 − e^{−s}), s = √(a_N/a₀), von Hand; F = M_L/(ν·M_b)."""
    s = math.sqrt(K.G * M_b / r**2 / a0)
    return M_L * (1 - math.exp(-s)) / M_b


_g("haufen_rest_bose", "haufen_rest", "zweiter_weg",
   [(100 * K.kpc, 2.1e13 * K.M_sun, 1.7e14 * K.M_sun, A0), (K.Mpc, 1.1e14 * K.M_sun, 5.9e14 * K.M_sun, A0)], _rest_bose,
   toleranz=1e-12, grund="Rundung von e^{−s}")


def _baryonen_quad(r, M_gas, r_c, r_norm, r_max, M_bcg, b_bcg, M_gal, b_gal):
    """β-Gas als Integral der Dichte ∝ x²/(1 + (x/r_c)²) (quad), Hernquist M·r²/(r + b)² von Hand."""
    from scipy.integrate import quad
    f = lambda x: x * x / (1 + (x / r_c)**2)
    gas = M_gas * quad(f, 0, min(r, r_max), epsrel=1e-12)[0] / quad(f, 0, r_norm, epsrel=1e-12)[0]
    return gas + M_bcg * r**2 / (r + b_bcg)**2 + M_gal * r**2 / (r + b_gal)**2


_AUFBAU = (1e14 * K.M_sun, 150 * K.kpc, K.Mpc, 3 * K.Mpc, 1.5e12 * K.M_sun, 20 * K.kpc, 1e13 * K.M_sun, 300 * K.kpc)
_g("haufen_baryonen_quad", "haufen_baryonen", "zweiter_weg", [(r,) + _AUFBAU for r in (100 * K.kpc, K.Mpc, 4 * K.Mpc)],
   _baryonen_quad, toleranz=1e-10, grund="Quadratur (epsrel 10⁻¹²) gegen geschlossene Form")


def _nfw_quad(r, M200, c_nfw, z):
    """NFW als Integral der Dichte ∝ 1/(x(1 + x)²) (quad), r₂₀₀ aus H(z) von Hand."""
    from scipy.integrate import quad
    H = K.H0_KMS * KMS * math.sqrt(K.OMEGA_M * (1 + z)**3 + K.OMEGA_R * (1 + z)**4 + 1 - K.OMEGA_M - K.OMEGA_R)
    r200 = (3 * M200 / (4 * math.pi * 200 * 3 * H**2 / (8 * math.pi * K.G)))**(1 / 3)
    f = lambda x: x / (1 + x)**2
    return M200 * quad(f, 0, c_nfw * r / r200, epsrel=1e-12)[0] / quad(f, 0, c_nfw, epsrel=1e-12)[0]


_g("masse_linse_nfw_quad", "masse_linse_nfw", "zweiter_weg",
   [(100 * K.kpc, 1e15 * K.M_sun, 4.0, 0.3), (K.Mpc, 1e15 * K.M_sun, 4.0, 0.3), (K.Mpc, 1e14 * K.M_sun, 7.0, 0.0)],
   _nfw_quad, toleranz=1e-10, grund="Quadratur (epsrel 10⁻¹²) gegen geschlossene Form")

