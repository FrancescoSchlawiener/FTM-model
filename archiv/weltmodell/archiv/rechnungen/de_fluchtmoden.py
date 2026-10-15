# -*- coding: utf-8 -*-
"""
Dunkle Energie als Kehrseite der Bindung, im Zustandsraum (Idee des Autors): Die Fluchtmoden, die die Regel (𝒫)
gegenüber Newton inaktiv macht, verschwinden im geschlossenen Ganzen nicht (Liouville/Kac); gehen sie in die
parametrische Richtung (Atmen)?
Blick: Das Maß μ ist das zugängliche Einteilchen-Phasenraummaß eines gebundenen Bereichs (je Bereich, G des Bereichs),
keine Dichte am Ort. Das Kac-Maß 1 + Q des Modells zählt in der inneren Richtung (I). Atmen = parametrisch (G, Zelle).
1 Je Bereich (Punktmasse, Hernquist): μ_N (Newton: gebunden bis ∞), μ_R (Regel: gebunden bis zur Naht R_ta),
  μ_inakt = ∫(v_R³ − v_N³)₊, je Baryon am Ort ⟨1 − (v_N/v_R)³⟩_M; gegen Masse und Epoche; Liouville mitbewegt.
2 Erhaltung: Zuwachs des parametrischen Maßes = inaktives Fluchtmaß der Bereiche (Massenfunktion), zwei Zählungen
  (A additiv, Norm = Maß des Bereichs im Hubble-Fluss; B multiplikativ, Liouville-Produkt je Baryon, normfrei);
  Übersetzung ins Atmen δ ln V = δ ln μ_par bei festem mitbewegtem Impuls ⇒ H_eff = H_b·(1 + ⅓ dX/d ln a).
3 Verriegelung T² → T¹: Frequenzverhältnis Ω_r/Ω_φ unter Regel und Newton, gesperrte Umläufe je Baryon, Rate.
4 Beobachterfehler und Doppelzählung.  5 Von Hand.
Aus ftm_kern: nu, r_stern, v_flach, hernquist_masse, umkehrradius_heute, dichte_umkehr_heute, omega_vakuum,
E_hintergrund, rho_m_mitbewegt, Konstanten.
Nicht im Kern (lokal, gemeldet): Umkehrdichte 1 + Δ_ta zur Epoche a (wie de_kehrseite_bindung), lineares Wachstum D(a),
σ(M) aus CAMB [G], Massenanteil erfc(δ_ta/√2σD), Tabelle δ_ta(a) (de_uebergabe_anteil), Phasenraummaß und
Epizykelfrequenz (Integrale über die Kern-Regel).
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.special import erfc, gamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, kpc, Mpc, OMEGA_M, OMEGA_B, OMEGA_R, OMEGA_DE, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM,
                      q_KOSM, H0_KMS, H0_si, E_hintergrund, rho_m_mitbewegt, omega_vakuum, nu, r_stern, v_flach,
                      hernquist_masse, umkehrradius_heute, dichte_umkehr_heute)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

H0 = H0_si()
NK = 1 + Q_KOSM
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
OM_C = OMEGA_M - OMEGA_B
RHO0 = rho_m_mitbewegt()
K6 = 16 * math.pi**2 / 3          # 4πr²dr · (4π/3)v³

# ------------------------------------------------ Umkehrdichte zur Epoche a (Top-Hat wie im Kern)
_UK = {}
def umkehr_dichte(a_ta):
    key = round(a_ta, 6)
    if key in _UK:
        return _UK[key]
    Om, OL = OMEGA_M, 1 - OMEGA_M
    Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
    t_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
    tz, ai = t_a(a_ta), 1e-3
    ti = t_a(ai)
    def lauf(di):
        Hi = Hb(ai); GM = 0.5 * Hi**2 * (Om / ai**3) / (Om / ai**3 + OL) * (1 + di)
        f = lambda t, y: [y[1], -GM / y[0]**2 + OL * H0**2 * y[0]]
        ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
        s = solve_ivp(f, (ti, 2 * tz), [1.0, Hi * (1 - di / 3)], events=ev, rtol=1e-9, atol=1e-12)
        return (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (2 * tz, None)
    ld = brentq(lambda x: lauf(10**x)[0] - tz, -4, 0.3, xtol=1e-10)
    di = 10**ld; R = lauf(di)[1]
    _UK[key] = (1 + di) / (R**3 * ai**3) * a_ta**3
    return _UK[key]

def R_ta(M_b, a):
    return (3 * NK * M_b / (4 * math.pi * RHO0 * a**-3 * umkehr_dichte(a)))**(1 / 3)

def Hz(a):
    return H0 * E_hintergrund(1 / a - 1)

# ------------------------------------------------ Phasenraummaß eines Bereichs
def masse(M_b, a, b=None, n=4000):
    """Zugängliches Einteilchenmaß bis R_ta. b = None: Punktmasse, sonst Hernquist-Skala b (M(<∞) = M_b).
    Gibt dict mit μ_N, μ_R, μ_inakt, μ_akt (Newton gebunden, Regel über die Naht), f_ort (massengewichtet), μ_fluss."""
    R = R_ta(M_b, a)
    r = np.geomspace(R * 1e-7, R, n)
    Mr = M_b * np.ones_like(r) if b is None else np.array([hernquist_masse(M_b, b, x) for x in r])
    aN = G * Mr / r**2
    acc = nu(aN) * aN
    I = cumulative_trapezoid(acc[::-1], -r[::-1], initial=0.0)[::-1]    # ∫_r^R a dr'
    vR3 = (2 * I)**1.5
    phiN = G * M_b / r if b is None else G * M_b / (r + b)
    vN3 = (2 * phiN)**1.5
    vNR3 = np.clip(2 * (phiN - (G * M_b / R if b is None else G * M_b / (R + b))), 0, None)**1.5   # Newton bis R_ta
    w = K6 * r**2
    muR = np.trapezoid(w * vR3, r); muN = np.trapezoid(w * vN3, r)
    muI = np.trapezoid(w * np.clip(vR3 - vN3, 0, None), r)
    muA = np.trapezoid(w * np.clip(vN3 - vR3, 0, None), r)
    muNR = np.trapezoid(w * vNR3, r)
    RL = R * umkehr_dichte(a)**(1 / 3)                                   # Lagrange-Radius (Maß im Fluss)
    muF = (4 * math.pi / 3)**2 * RL**6 * Hz(a)**3
    out = dict(R=R, muR=muR, muN=muN, muI=muI, muA=muA, muNR=muNR, muF=muF)
    if b is not None:
        dM = np.gradient(Mr, r)
        vRs = np.where(vR3 > 0, vR3, np.inf)
        out['f_ort'] = float(np.trapezoid(np.clip(1 - vN3 / vRs, 0, None) * dM, r) / Mr[-1])
        out['V_ort'] = float(np.trapezoid(vR3 / vN3 * dM, r) / Mr[-1])
    return out

# ============================================================== 1 Je Bereich
kopf("1  Inaktives Fluchtmaß je gebundenem Bereich (Einteilchen-Phasenraum bis zur Naht R_ta, Kern-Regel ν)")
Mw = 6e10 * M_sun
P("   Probe: R_ta(Milchstraße, a = 1) = %.3f Mpc (Kern %.3f), 1 + Δ_ta = %.3f (Kern %.3f)"
  % (R_ta(Mw, 1.0) / Mpc, umkehrradius_heute(NK * Mw) / Mpc, umkehr_dichte(1.0), dichte_umkehr_heute()))
P("   μ = ∫4πr²·(4π/3)v_esc³ dr; Newton v_N² = 2GM/r (E > 0 entkommt), Regel v_R² = 2∫_r^{R_ta}ν(a_N)a_N dr' (Naht = Rand).")
P("   Tief (R_ta ≫ r*): μ_R/μ_N → (3/2)·Γ(5/2)/3^{5/2}·(R_ta/r*)^{3/2} = %.4f·(R_ta/r*)^{3/2}" % (1.5 * gamma(2.5) / 3**2.5))
P()
P("   Punktmasse")
P("   M_b [M☉]   a     R_ta[kpc]  R_ta/r*   μ_R/μ_N   Tiefformel  μ_inakt/μ_R   μ_akt/μ_R   μ_R/μ_N(bis R_ta)  ln(μ_R/μ_N)  μ_inakt/μ_Fluss")
TAB1 = {}
for lm in (8, 9, 10, 10.778, 12, 13, 14):
    Mb = 10**lm * M_sun
    for a in (1.0, 0.5, 0.25):
        m = masse(Mb, a); x = m['R'] / r_stern(Mb)
        TAB1[(lm, a)] = m
        P("   10^%-6.3f  %.2f  %8.0f  %7.1f   %8.2f   %8.2f     %.4f       %.2e      %8.2f          %6.3f       %.3f"
          % (lm, a, m['R'] / kpc, x, m['muR'] / m['muN'], 1.5 * gamma(2.5) / 3**2.5 * x**1.5, m['muI'] / m['muR'],
             m['muA'] / m['muR'], m['muR'] / m['muNR'], math.log(m['muR'] / m['muN']), m['muI'] / m['muF']))
P("   μ_R/μ_N ∝ (R_ta/r*)^{3/2} ∝ M_b^{−1/4}: kleine und späte Bereiche machen je Baryon mehr Flucht inaktiv.")
P("   Kac: μ_R/μ_N ist die mittlere Rückkehrzeit (in Einheiten des Umlaufs) in den Newton-gebundenen Teil; 1 + Q = %.3f" % NK)
P("   liegt nur bei R_ta/r* = %.1f (M_b ≈ 10^%.1f M☉ heute) – ein Einzelwert, keine Invariante."
  % ((NK / (1.5 * gamma(2.5) / 3**2.5))**(2 / 3),
     brentq(lambda l: masse(10**l * M_sun, 1.0)['muR'] / masse(10**l * M_sun, 1.0)['muN'] - NK, 10, 20)))
P()
P("   Profil (Hernquist, Skala b): Maß des Bereichs und je Baryon am Ort (Geschwindigkeitsraum an der Stelle der Masse)")
P("   M_b [M☉]  b[kpc]  a     μ_R/μ_N   μ_inakt/μ_R   ⟨1 − (v_N/v_R)³⟩_M   ⟨(v_R/v_N)³⟩_M")
TAB1P = {}
for lm, bk in ((9, 1.0), (10.778, 3.0), (12, 10.0)):
    Mb = 10**lm * M_sun
    for a in (1.0, 0.5, 0.25):
        m = masse(Mb, a, b=bk * kpc); TAB1P[(lm, a)] = m
        P("   10^%-6.3f %5.1f  %.2f   %8.2f     %.4f           %.3f               %.3f"
          % (lm, bk, a, m['muR'] / m['muN'], m['muI'] / m['muR'], m['f_ort'], m['V_ort']))
P("   Am Ort der Baryonen (ν klein) ist die inaktive Flucht klein; das Bereichsmaß sitzt fast ganz jenseits von r*, wo keine Baryonen sind.")
P("   Liouville: nach der Umkehr sind R_ta und v_R physikalisch fest; mitbewegt (R/a)³(a·v)³ = R³v³ – das Maß je Bereich ist")
P("   in beiden Koordinaten gleich und zeitlich fest. Es wächst nur durch neue Bindung (dF/dt) und größere Bereiche.")

# ============================================================== 2 Massenfunktion und Atmen
kopf("2  Erhaltung im Ganzen: Zuwachs des parametrischen Maßes = inaktives Fluchtmaß (Messwerte nur Test [G])")
def wachstum():
    def f(x, y):
        a = math.exp(x); h = 1e-5
        E = lambda aa: E_hintergrund(1 / aa - 1)
        dlnE = (math.log(E(a * math.exp(h))) - math.log(E(a * math.exp(-h)))) / (2 * h)
        Oma = OMEGA_M * a**-3 / E(a)**2
        return [y[1], -(2 + dlnE) * y[1] + 1.5 * Oma * y[0]]
    x0 = math.log(1e-3)
    s = solve_ivp(f, (x0, 0.0), [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-12)
    D1 = s.sol(0.0)[0]
    return lambda a: s.sol(math.log(a))[0] / D1
Dg = wachstum()

import camb
p = camb.CAMBparams()
p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=0.06, tau=0.0544)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=400)
PD = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False, kmax=400,
                                        var1='delta_tot', var2='delta_tot', zmax=1)
KS = np.exp(np.linspace(math.log(1e-4), math.log(400.0), 5000)); P0 = PD.P(0.0, KS)
def sig0(Mg):
    R = (3 * Mg / (4 * math.pi * RHO0))**(1 / 3) / Mpc
    x = KS * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(np.trapezoid(KS**3 * P0 / (2 * math.pi**2) * W**2, np.log(KS)))

LM = np.linspace(6.0, 16.5, 106)                     # log10 M_ges
SIG = np.array([sig0(10**l * M_sun) for l in LM])
P("   σ(M_ges, 0) [G, CAMB]: 10⁸ %.3f, 10¹² %.3f, 10¹⁴ %.3f" % tuple(float(np.interp(l, LM, SIG)) for l in (8, 12, 14)))
ZT = np.array([10, 6, 3, 2, 1, 0.5, 0.4, 0.0]); DT = np.array([1.0741, 1.0709, 1.0734, 1.0815, 1.1132, 1.1624, 1.1785, 1.2766])
def delta_ta(a): return float(np.interp(-math.log(a), np.log(1 + ZT[::-1]), DT[::-1]))

AA = np.exp(np.linspace(math.log(1 / 11), math.log(1.1), 23))      # bis a = 1,1 für die Ableitung bei a = 1
LMg = np.linspace(6.0, 16.5, 22)
XT = {}
for a in AA:
    rows = []
    for l in LMg:
        m = masse(10**l * M_sun / NK, a, n=1500)
        rows.append((m['muI'] / m['muF'], math.log(m['muR'] / m['muN']), m['muI'] / m['muR']))
    XT[a] = np.array(rows)

def X_von_a(a, lm_min):
    Fg = erfc(delta_ta(a) / (math.sqrt(2) * SIG * Dg(a)))
    dF = -np.gradient(Fg, LM * math.log(10))
    t = XT[a]
    xa = np.interp(LM, LMg, t[:, 0]); xb = np.interp(LM, LMg, t[:, 1]); xf = np.interp(LM, LMg, t[:, 2])
    msk = LM >= lm_min
    I = lambda y: float(np.trapezoid((dF * y)[msk], LM[msk] * math.log(10)))
    return I(xa), I(xb), I(xf), float(np.interp(lm_min, LM, Fg))

LA = np.log(AA)
EB = lambda a: math.sqrt(OMEGA_M * a**-3 + OMEGA_R * a**-4 + OM_S)            # Hintergrund ohne die fehlenden 0,442
ELCDM = lambda a: E_hintergrund(1 / a - 1)
i1 = int(np.argmin(abs(AA - 1.0)))
P("   F = erfc(δ_ta(a)/(√2σD)) (Bindung ab der Umkehr, Regel setzt dort ein); je Bereich μ bei R_ta(M, a).")
P("   (A) additiv: X_A = Σ dF/dlnM · μ_inakt/μ_Fluss, μ_Fluss = (4π/3)R_L³·(4π/3)(H R_L)³ (Maß des Bereichs im Hubble-Fluss).")
P("   (B) Liouville-Produkt je Baryon: X_B = Σ dF/dlnM · ln(μ_R/μ_N) (normfrei).")
P("   Atmen: δ ln μ_par = X bei festem mitbewegtem Impuls ⇒ δ ln V = X, a_eff = a·e^{X/3}, H_eff = H_b(1 + ⅓ dX/d ln a).")
P("   H_b: Ω_m + Ω_s = %.4f (E_b(1) = %.3f); nötig ⅓dX/dln a (a=1) = %.3f, damit E_eff(1) = 1."
  % (OMEGA_M + OM_S, EB(1.0), 1 / EB(1.0) - 1))
P()
P("   M_min   F(0)   ⟨inakt/μ_R⟩  X_A(z=3)  X_A(1)  X_A(0)  dX_A/dlna(0)  Ω_X,A(0)  /0,442   |  X_B(3)  X_B(1)  X_B(0)  dX_B/dlna(0)  Ω_X,B(0)  /0,442")
RES = {}
for lmin in (6, 8, 10, 12, 13):
    v = np.array([X_von_a(a, lmin) for a in AA])
    XA, XB, XF, FF = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    dA = np.gradient(XA, LA); dB = np.gradient(XB, LA)
    EA = np.array([EB(a) for a in AA]) * (1 + dA / 3); EBB = np.array([EB(a) for a in AA]) * (1 + dB / 3)
    OXA = EA**2 - np.array([EB(a)**2 for a in AA]); OXB = EBB**2 - np.array([EB(a)**2 for a in AA])
    RES[lmin] = (XA, XB, dA, dB, EA, EBB, OXA, OXB, FF)
    zi = lambda z: int(np.argmin(abs(AA - 1 / (1 + z))))
    P("   10^%-3d  %.3f   %.4f      %.4f   %.4f  %.4f   %.4f       %.4f    %.3f    |  %.3f   %.3f   %.3f   %.3f        %.3f     %.2f"
      % (lmin, FF[i1], XF[i1] / FF[i1], XA[zi(3)], XA[zi(1)], XA[i1], dA[i1], OXA[i1], OXA[i1] / OM_X,
         XB[zi(3)], XB[zi(1)], XB[i1], dB[i1], OXB[i1], OXB[i1] / OM_X))
P()
P("   Verlauf mit a (B, M_min = 10⁸ bzw. 10¹²): E_eff/E_ΛCDM [G] und Ω_X,B(a) = E_eff² − E_b²")
for lmin in (8, 12):
    XA, XB, dA, dB, EA, EBB, OXA, OXB, FF = RES[lmin]
    s = "   10^%-3d" % lmin
    for z in (6, 3, 2, 1, 0.5, 0):
        i = int(np.argmin(abs(AA - 1 / (1 + z))))
        s += "  z=%-3g: %.3f (Ω_X %.3f)" % (z, EBB[i] / ELCDM(AA[i]), OXB[i])
    P(s)
P("   Vorzeichen: X > 0, wächst mit der Bindung ⇒ zusätzliche Ausdehnung. Verlauf: dX/dln a folgt dF/dln a, Gipfel früh")
P("   (kleine Bereiche binden bei z ≈ 2–6), heute fallend – also mitlaufend mit der Strukturbildung, nicht fest wie w ≈ −1.")
P("   (B) gibt heute 0,18–0,33 (41–75 % von 0,442), hebt aber H(z) bei z = 1–3 um 22–44 % über ΛCDM [G] (Ω_X(z=2) ≈ 10):")
P("   ein Glied, das in die Vergangenheit steil wächst, statt heute zu überwiegen. Mit größerer Norm wird z = 1–3 nur schlimmer.")
P("   Ω_X/Ω_c heute (Ω_DE/Ω_Mode gemessen %.3f/%.3f = %.2f):" % (OMEGA_DE, OM_C, OMEGA_DE / OM_C))
for lmin in (8, 12):
    XA, XB, dA, dB, EA, EBB, OXA, OXB, FF = RES[lmin]
    P("   10^%-3d  (A) %.4f   (B) (Ω_s + Ω_X,B)/Ω_c = %.2f" % (lmin, (OM_S + OXA[i1]) / OM_C, (OM_S + OXB[i1]) / OM_C))
P("   Norm: (A) hängt an der Wahl des Flussmaßes (Kugel R_L, Geschwindigkeit H·R_L, Faktoren (4π/3)²) – aus dem Modell nicht fest;")
P("   μ_Fluss ∝ (a·ȧ)³ ist kein Liouville-Invariant (der Fluss ist kalt, keine freie Bewegung): X_A driftet mit der Epoche,")
P("   dX_A/dln a < 0 heute – das Vorzeichen setzt die Norm, nicht die Bindung. (A) ist keine tragfähige Zählung.")
P("   (B) ist normfrei, verlangt aber, dass das Maß je Baryon als Faktor (Produktraum) ganz ins Zellvolumen geht.")

# ============================================================== 3 Verriegelung
kopf("3  Verriegelung T² → T¹: Frequenzverhältnis der Bahnen und gesperrte Umläufe je Baryon")
def kappa_zu_omega(M, r):
    h = 1e-4
    Om2 = lambda x: float(nu(G * M / x**2) * G * M / x**2) / x
    dOm2 = (Om2(r * (1 + h)) - Om2(r * (1 - h))) / (2 * h * r)
    return math.sqrt(r * dOm2 + 4 * Om2(r)) / math.sqrt(Om2(r))
P("   Fast kreisförmig, Punktmasse: Ω_r/Ω_φ = κ/Ω (Kepler exakt 1: Radial- und Umlauf verriegelt, Bahn geschlossen, T¹)")
s = "   r/r* = "
for x in (0.03, 0.3, 1, 3, 10, 30, 100):
    s += " %g: %.4f " % (x, kappa_zu_omega(Mw, x * r_stern(Mw)))
P(s)
P("   Tief → √2 (irrational): die Regel löst die Kepler-Verriegelung, die Bahn füllt einen T² (Rosette).")
P("   Gesperrte Umläufe je gebundenem Baryon (kugelsymmetrisch, 3 Freiheitsgrade): Newton 2 (Ebene + Kepler), Regel 1 (Ebene).")
P("   Gegenüber Newton gibt die Regel je Baryon einen Umlauf frei (−1), gegenüber dem ungebundenen Fluss (keine Umläufe)")
P("   schafft Bindung Umläufe, sie sperrt keine. Anzahl je Baryon heute: F(>M_min, 0) · 1:")
for lmin in (8, 12):
    FF = RES[lmin][8]
    dF = np.gradient(FF, LA)
    P("   M_min 10^%-3d  F(0) = %.3f, dF/dln a (0) = %.3f;  als Rate dN/dt/N_b = %.3f·H₀" % (lmin, FF[i1], dF[i1], dF[i1]))
vf = v_flach(Mw)
P("   Übersetzung in Ausdehnung braucht Wirkung oder Energie je Umlauf; das Modell setzt keine. Größenordnungen:")
P("   Wirkung je Masse: Bahn v_f·r* = %.2e m²/s gegen Atmen der Zelle H₀·Λ² = %.2e m²/s (Verhältnis %.1e)"
  % (vf * r_stern(Mw), H0 * (c**2 / (2 * A0))**2, vf * r_stern(Mw) / (H0 * (c**2 / (2 * A0))**2)))
P("   Energie: Ω_b·F·(v_f/c)² = %.4f·%.2f·%.2e = %.1e gegen 0,442" % (OMEGA_B, RES[12][8][i1], (vf / c)**2, OMEGA_B * RES[12][8][i1] * (vf / c)**2))

# ============================================================== 4 Blick
kopf("4  Beobachterfehler und Doppelzählung")
P("   μ ist das Maß des Zustandsraums eines Bereichs (G des Bereichs); über die Zelle summiert G. Nirgends eine Dichte am Ort.")
P("   Am Ort der Baryonen (Profil) ist die inaktive Flucht klein; groß wird sie nur im leeren Außenteil bis R_ta: das")
P("   Bereichsmaß ist eine Eigenschaft der Regel (Naht), kein Stoff.")
P("   Kac: 1 + Q = %.3f zählt die Rückkehr in der inneren Richtung (I) und ist als Mode ξ schon als Staub (Ω_c) gebucht;" % NK)
P("   das Bahnmaß μ_R/μ_N (L/G des Bereichs) ist davon verschieden (Faktor 1–10³, massenabhängig). Wer μ_inakt zusätzlich ins")
P("   Atmen bucht, zählt die Bindung zweimal, falls sie dasselbe Maß ist wie q = %.3f; der Vakuumteil hängt nur an Λ (I)." % q_KOSM)
P("   𝒫 löscht keine Freiheitsgrade: die inaktiven Richtungen bleiben Richtungen desselben Bereichs (Rückkehr), sie")
P("   verlassen ihn nicht. Ein Übergang in die Rolle der Zelle ist zusätzlich gesetzt, nicht aus Liouville gefolgert.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand (Milchstraße, M_b = 6·10¹⁰ M☉, heute)")
rs = r_stern(Mw); R = R_ta(Mw, 1.0)
P("   r* = √(GM/a₀) = √(%.4e·6e10/1,1e-10) = %.3e m = %.2f kpc; R_ta = %.1f kpc; R_ta/r* = %.2f" % (G * M_sun, rs, rs / kpc, R / kpc, R / rs))
P("   ∫₀¹u²(ln 1/u)^{3/2}du = Γ(5/2)/3^{5/2} = 1,32934/15,5885 = %.5f" % (gamma(2.5) / 3**2.5))
P("   μ_R/μ_N ≈ (3/2)·0,08528·(R_ta/r*)^{3/2} = 0,12792·%.1f = %.1f  (Code, voller Kern-ν: %.1f)"
  % ((R / rs)**1.5, 0.12792 * (R / rs)**1.5, TAB1[(10.778, 1.0)]['muR'] / TAB1[(10.778, 1.0)]['muN']))
rq = TAB1[(10.778, 1.0)]['muR'] / TAB1[(10.778, 1.0)]['muN']
P("   Rest zur Tiefformel aus ν ≈ 1/s + ½ + …: v_R² = 2v_f²ln(R/r) + GM(1/r − 1/R) + …; relativ %.1f %%." % (100 * (rq / (0.12792 * (R / rs)**1.5) - 1)))
P("   ln(μ_R/μ_N) = ln %.1f = %.3f." % (rq, math.log(rq)))
P("   Rotation: κ/Ω tief = √(r·d(v_f²/r²)/dr/(v_f²/r²) + 4) = √(−2 + 4) = √2 = 1,4142.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_fluchtmoden.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
