# -*- coding: utf-8 -*-
"""
Dunkle Energie als Kehrseite der Bindung (Idee des Autors): Erscheint der Überschuss der gezählten Faltung
Δs_eff = ν·Δs gebundener Bereiche über das Maß des Ganzen (1 + Q)·Δs_b als Dehnung der Raumschleifen?
Blick: Δs_eff ist gelesen am Ort (L); 1 + Q ist Kac-Maß (I, zählt im Ganzen G); ΔL_Raum ist G.
1 Überschuss je Bereich: e = ν(a_N(R_ta)) − (1 + Q) (an R_ta gelesen, gleich für jedes Profil innerhalb R_ta);
  massengewichtet ∫(ν(r) − (1 + Q))dM/M_b (Profil Hernquist); Lesung von außen (jenseits 1,05·R_ta): (1 + wQ) − (1 + Q).
2 Über die Massenfunktion: x(a) = Σ dF/dlnM · e(M, a)/(1 + Q), F aus σ(M) (CAMB [G]), Schwelle δ_F bzw. δ_ta.
3 Ansatz ΔL_Raum = X, zwei Übersetzungen ins Ausdehnungsgesetz (Bremsglied bzw. Λ), Lesart L und G.
4 Atmen der Zentren (Takt, Engstelle): Größenordnung.
5 Von Hand.
Aus ftm_kern: nu, ds, takt, r_stern, r_Q, umkehrradius_heute, dichte_umkehr_heute, hernquist_masse, gewicht_w,
omega_vakuum, schliesslaenge, E_hintergrund, rho_m_mitbewegt, Konstanten.
Nicht im Kern (hier lokal, gemeldet): Umkehrdichte 1 + Δ_ta zur Epoche a (Verallgemeinerung von dichte_umkehr_heute),
lineares Wachstum D(a), σ(M) (CAMB wie de_uebergabe_anteil), Massenanteil erfc(δ/√2σD) und seine Ableitung.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import erfc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, kpc, Mpc, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM, H0_KMS,
                      H0_si, E_hintergrund, rho_m_mitbewegt, omega_vakuum, nu, ds, takt, r_stern, r_Q,
                      umkehrradius_heute, dichte_umkehr_heute, hernquist_masse, gewicht_w, schliesslaenge)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 110); P(t); P("=" * 110)

H0 = H0_si()
NK = 1 + Q_KOSM
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
RHO0 = rho_m_mitbewegt()

# ------------------------------------------------ Umkehrdichte zur Epoche a (Top-Hat wie im Kern, Umkehr bei a_ta)
_UK = {}
def umkehr_dichte(a_ta):
    """1 + Δ_ta gegen ρ_m(a_ta) für eine Schale, die bei a_ta umkehrt (gleiches Verfahren wie dichte_umkehr_heute)."""
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
    """Umkehrradius (physikalisch) des Bereichs mit Maß (1 + Q)·M_b, der bei a umkehrt."""
    return (3 * NK * M_b / (4 * math.pi * RHO0 * a**-3 * umkehr_dichte(a)))**(1 / 3)

def e_Rta(M_b, a):
    """Überschuss der gezählten Faltung an R_ta je Baryonenfaltung: ν(a_N(R_ta)) − (1 + Q)."""
    return float(nu(G * M_b / R_ta(M_b, a)**2)) - NK

# ============================================================== 1 Je Bereich
kopf("1  Überschuss der gezählten Faltung je gebundenem Bereich (Einheit Δs_b = GM_b/3c²)")
Mw = 6e10 * M_sun
P("   Probe Kern: 1 + Δ_ta(a=1) = %.3f (Kern %.3f), R_ta(Milchstraße) = %.3f Mpc (Kern %.3f)"
  % (umkehr_dichte(1.0), dichte_umkehr_heute(), R_ta(Mw, 1.0) / Mpc, umkehrradius_heute(NK * Mw) / Mpc))
w_aus = gewicht_w(3 * H0)
P("   Lesung von außen (jenseits 1,05·R_ta, Gauß): (1 + wQ)·Δs_b, w = %.4f (mitbewegt, θ̃ = 3H) ⇒ e_G = −(1 − w)Q = %+.4f" % (w_aus, -(1 - w_aus) * Q_KOSM))
P()
P("   M_b [M☉]   a     1+Δ_ta   R_ta [kpc]  a_N/a₀ an R_ta   ν(R_ta)   e_L = ν − (1+Q)   e_L/(1+Q)   R_Q [kpc]  R_Q/R_ta")
for lm in (8, 9, 10, 10.778, 11, 12, 13, 14):
    Mb = 10**lm * M_sun
    for a in (1.0, 0.5, 0.25):
        R = R_ta(Mb, a); aN = G * Mb / R**2; e = e_Rta(Mb, a)
        P("   10^%-6.3f  %.2f  %6.2f   %8.0f    %.2e        %7.2f   %8.2f          %6.2f      %6.1f     %.3f"
          % (lm, a, umkehr_dichte(a), R / kpc, aN / A0, e + NK, e, e / NK, r_Q(Mb) / kpc, r_Q(Mb) / R))
P("   Tief (a_N ≪ a₀): ν ≈ 1/s = R_ta/r* ∝ M_b^{−1/6}; kleinere Bereiche zählen je Baryon mehr.")

P()
P("   Massengewichtet (ν am Ort jeder Baryonenschale, Hernquist, Skala b): ⟨ν⟩_M − (1 + Q) = ∫(ν(r) − (1+Q)) dM/M_b bis R_ta")
for lm, b in ((9, 1.0), (10.778, 3.0), (12, 10.0)):
    Mb = 10**lm * M_sun; b *= kpc; R = R_ta(Mb, 1.0)
    Mtot = Mb * (R + b)**2 / R**2           # so, dass M(<R_ta) = M_b
    rr = np.geomspace(1e-3 * b, R, 4000)
    Mr = np.array([hernquist_masse(Mtot, b, r) for r in rr])
    dM = np.gradient(Mr, rr)
    nv = nu(G * Mr / rr**2)
    mv = np.trapezoid(nv * dM, rr) / Mb
    P("   M_b = 10^%.2f, b = %4.1f kpc: ⟨ν⟩_M = %.3f ⇒ Überschuss %+.3f (%+.3f·(1+Q))" % (lm, b / kpc, mv, mv - NK, (mv - NK) / NK))
P("   Nach der Masse gewichtet liegt ν bei den Baryonen unter 1 + Q: Überschuss nur jenseits R_Q, dort ohne Baryonen.")

# ============================================================== 2 Massenfunktion
kopf("2  Über die Massenfunktion: x(a) = Σ_M (dF/dlnM)·e(M, a)/(1 + Q)   (Anteil der gebundenen Masse je lnM)")
def E(a): return E_hintergrund(1 / a - 1)
def wachstum():
    # D'' + (2 + dlnH/dlna) D' = 1,5 Ω_m(a) D  in ln a
    def f(x, y):
        a = math.exp(x); h = 1e-5
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

LM = np.linspace(5.0, 16.5, 231)                          # log10 M_ges [M☉]
SIG = np.array([sig0(10**l * M_sun) for l in LM])
P("   σ(M_ges, 0) [G, CAMB]: 10⁸ %.3f, 10¹² %.3f, 10¹⁴ %.3f (de_uebergabe_anteil: 5,898 / 2,245 / 1,014)"
  % tuple(float(np.interp(l, LM, SIG)) for l in (8, 12, 14)))
# Schwellen aus de_uebergabe_anteil (Kern-Hintergrund): δ_F(I_F = 1,10) ≈ 1,560–1,575; δ_ta 1,07–1,28
DELTA_F = 1.565
ZT = np.array([10, 6, 3, 2, 1, 0.5, 0.4, 0.0]); DT = np.array([1.0741, 1.0709, 1.0734, 1.0815, 1.1132, 1.1624, 1.1785, 1.2766])
def delta_ta(a): return float(np.interp(-math.log(a), np.log(1 + ZT[::-1]), DT[::-1]))

AA = np.array([1 / (1 + z) for z in (10, 6, 4, 3, 2, 1.5, 1, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0)])
EM = {}
def e_tab(a):
    if a not in EM:
        EM[a] = np.array([e_Rta(10**l * M_sun / NK, a) for l in LM])
    return EM[a]

def x_von_a(a, lm_min, schwelle, e_fn=None):
    """x(a) = ∫_{M_min} (−dF(>M)/dlnM)·e/(1+Q) dlnM, F(>M) = erfc(δ/(√2 σ D)); dazu F(>M_min)."""
    d = DELTA_F if schwelle == "F" else delta_ta(a)
    Fg = erfc(d / (math.sqrt(2) * SIG * Dg(a)))
    dF = -np.gradient(Fg, LM * math.log(10))
    e = e_tab(a) if e_fn is None else e_fn(a)
    m = LM >= lm_min
    return float(np.trapezoid((dF * e / NK)[m], LM[m] * math.log(10))), float(np.interp(lm_min, LM, Fg))

P()
P("   Lesart L (an R_ta gelesen), Schwelle δ_F = %.3f (Übergabe I_F):" % DELTA_F)
P("   M_min(ges)   x(z=6)    x(3)     x(2)     x(1)     x(0,5)   x(0)     F(>M_min, 0)   ⟨e/(1+Q)⟩ der Gebundenen heute")
XL = {}
for lmin in (6, 8, 10, 12, 13):
    xs = [x_von_a(a, lmin, "F")[0] for a in AA]
    XL[lmin] = np.array(xs)
    F0 = x_von_a(1.0, lmin, "F")[1]
    P("   10^%-3d     %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f     %.3f          %.2f"
      % (lmin, xs[1], xs[3], xs[4], xs[6], xs[8], xs[-1], F0, xs[-1] / F0))
P("   Schwelle δ_ta (Bindung ab der Umkehr, w = 0):")
XT = {}
for lmin in (8, 12, 13):
    xs = [x_von_a(a, lmin, "T")[0] for a in AA]; XT[lmin] = np.array(xs)
    P("   10^%-3d     x(1) = %.3f  x(0,5) = %.3f  x(0) = %.3f" % (lmin, xs[6], xs[8], xs[-1]))
P("   Lesart G (von außen jenseits 1,05·R_ta): x_G(a) = −(1 − w)Q/(1 + Q)·F(a), |x_G| ≤ %.4f" % ((1 - w_aus) * Q_KOSM / NK))

# ============================================================== 3 Ausdehnung
kopf("3  Ansatz ΔL_Raum = X: Beitrag zur Ausdehnung (Messwerte nur Test [G])")
LA = np.log(AA)
def w0wa(rho):
    lr = np.log(np.abs(rho)); w = -1 - np.gradient(lr, LA) / 3
    wa = -(w[-1] - w[-2]) / (AA[-1] - AA[-2])
    return w[-1], wa
P("   (A) Bremsglied: σ = ΣΔs − ΔL = 0 wird im Ausdehnungsgesetz als Staub gelesen (01); Überschuss als Zusatz Ω_X(a) = Ω_m·x(a)·a⁻³.")
P("       Nötig heute Ω_X = %.3f ⇒ x(0) = %.3f." % (OM_X, OM_X / OMEGA_M))
P("   M_min     Ω_X(0)    Ω_X(0)/0,442   ρ_DE = Ω_s + Ω_X: w₀, w_a     ρ_X(z=1)/ρ_X(0)")
for lmin in (6, 8, 10, 12, 13):
    ox = OMEGA_M * XL[lmin] * AA**-3
    w0, wa = w0wa(OM_S + ox)
    P("   10^%-3d   %7.3f    %7.2f        %+.3f, %+.3f              %.2f" % (lmin, ox[-1], ox[-1] / OM_X, w0, wa, ox[6] / ox[-1]))
P("   (B) Λ-Übersetzung: Λ = ν_kosm·ΣΔs_b ≙ ΔL_Raum (03); mit Überschuss Λ(a) = Λ·(1 + x(a)) ⇒ a₀ ∝ 1/(1 + x), Ω_s ∝ (1 + x)⁻².")
for lmin in (6, 12, 13):
    x = XL[lmin]
    P("   10^%-3d   a₀(0)/a₀(z=6) = %.3f   Ω_s(0) = %.4f statt %.4f   (Modell: a₀ seit z ≈ 10 fest auf 0,3 %%)"
      % (lmin, (1 + x[1]) / (1 + x[-1]), OM_S / (1 + x[-1])**2 * (1 + x[1])**2, OM_S))
P("   Vorzeichen (B): Λ wächst ⇒ Vakuumteil fällt ⇒ w > −1, aber Betrag verschiebt a₀ um Faktoren; ausgeschlossen.")
P("   Lesart G: Ω_X = Ω_m·x_G·a⁻³ ≥ %.4f (negativ, Staubverlust) bzw. Λ-Änderung ≤ %.1e relativ." % (-OMEGA_M * (1 - w_aus) * Q_KOSM / NK, (1 - w_aus) * Q_KOSM / NK))

# ============================================================== 3b Blick
kopf("3b Welche Lesart verlangt das Modell")
P("   Δs_eff = a·r²/(3c²) ist die an r gelesene Faltung (L); G₀₀ = (6/r²)·dΔs_eff/dr: der Überschuss innerhalb R_ta ist der")
P("   Strukturterm der Regel, keine Faltung (die Regel faltet nicht). Integriert über den ganzen Bereich einschließlich der")
P("   Naht gibt Gauß nur Δs_eff(außen) − 0 = (1 + wQ)·Δs_b: der Überschuss hebt sich in der Naht auf (dΔs_eff/dr < 0 dort).")
P("   Triade: Beim Aufsteigen in der inneren Richtung bleiben Fluss GM und Hintergrund der Eltern, die eigene aufgelöste")
P("   Struktur klingt ab. Die Ebene darüber (G) sieht (1 + Q)·M_b. Kac: ν_kosm = 1/(1 − q) mittelt über den ganzen Zustandsraum;")
P("   Übergabe ändert die Form, nicht das Maß. Λ zählt in I, wo nichts atmet. ⇒ Das Modell verlangt Lesart G: X_G ≈ 0.")
xm = e_tab(1.0)
P("   Probe (L): Summe Σ_ges e/(1+Q) über alle Ebenen der Hierarchie zählte jede Stufe erneut: mit M_min → klein x → %.1f (10⁶)." % XL[6][-1])

# ============================================================== 4 Zentren
kopf("4  Atmen der innersten Ebene (Zentren, Engstelle r₀ = 6Δs)")
P("   Faltung eines Zentrums ist Δs = GM/3c², additiv; Einfall in ein Zentrum verschiebt Faltung, ändert ΣΔs nicht.")
P("   ΣΔs ändert sich nur, wo Ruhe-Energie als Strahlung weggeht (Strahlung atmet mit, a⁻⁴): Faltung nimmt ab.")
# Größenordnungen [G]: Sternanteil der Baryonen f*, verarbeiteter Anteil in Fusion, Wirkungsgrad; Anteil Zentralmassen
f_st, f_fus, eps_fus = 0.07, 0.1, 0.007       # [G]
f_bh, eps_bh = 1e-4, 0.1                      # [G] Zentralmassen/Baryonen, Abstrahlung beim Einfall
dSig = f_st * f_fus * eps_fus + f_bh * eps_bh / (1 - eps_bh)
P("   Abgestrahlt: Sterne %.0e·%.1f·%.3f + Zentren %.0e·%.2f ⇒ ΔΣΔs/ΣΔs ≈ −%.1e [G]" % (f_st, f_fus, eps_fus, f_bh, eps_bh / (1 - eps_bh), dSig))
P("   Als Λ-Änderung: Λ fällt ⇒ Ω_s steigt um ≈ %.1e·Ω_s = %.1e (w < −1, Betrag 10⁻⁴ von 0,442)." % (2 * dSig, 2 * dSig * OM_S))
# Takt: massengewichtete Tiefe δ = 1 − dτ/dt
Msun = M_sun
d_sun = 1 - takt(Msun, 6.957e8); d_wd = 1 - takt(0.6 * Msun, 0.012 * 6.957e8); d_ns = 1 - takt(1.4 * Msun, 1.2e4)
d_gal = 1 - takt(6e10 * Msun, 8.2 * kpc)
P("   Takt δ = 1 − dτ/dt: Sonne %.1e, Weißer Zwerg %.1e, Neutronenstern %.2f, Galaxie am Ort der Sonne (Kern) %.1e" % (d_sun, d_wd, d_ns, d_gal))
dmean = f_st * (0.9 * d_sun + 0.08 * d_wd + 0.01 * d_ns + 0.01 * 1.0) + f_bh * 1.0
P("   Massenmittel über Baryonen [G: 7 %% in Sternen, davon 8 %% WD, 1 %% NS, 1 %% Falten; 10⁻⁴ Zentralmassen]: ⟨δ⟩ ≈ %.1e" % dmean)
P("   Engstelle r₀: Faltung je Zentrum fest (Δs), das Atmen der Zentren ist ihr Wachstum (Umlagerung, ΣΔs fest bis auf Abstrahlung).")
P("   Größenordnung jedes Glieds der innersten Ebene: ≤ 10⁻³ von ΣΔs, gegen x(0) = 1,41 nötig.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand")
Rm = R_ta(Mw, 1.0)
P("   Milchstraße: M_ges = 6,364·6·10¹⁰ = 3,82·10¹¹ M☉; ρ_m = %.3e kg/m³, 1 + Δ_ta = %.2f" % (RHO0, umkehr_dichte(1.0)))
P("   R_ta = (3·3,82e11·1,989e30/(4π·%.3e·%.2f))^{1/3} = %.3f Mpc" % (RHO0, umkehr_dichte(1.0), Rm / Mpc))
aN = G * Mw / Rm**2
P("   a_N = GM_b/R² = 7,96e30/(%.3e)² = %.3e m/s² = %.3e·a₀; s = √ = %.4f; ν = 1/(1 − e^{−s}) ≈ 1/s + ½ + s/12 = %.2f"
  % (Rm, aN, aN / A0, math.sqrt(aN / A0), 1 / math.sqrt(aN / A0) + 0.5 + math.sqrt(aN / A0) / 12))
P("   Code ν = %.2f; e_L/(1+Q) = (ν − 6,364)/6,364 = %.2f; von außen −(1 − 0,9854)·5,364/6,364 = %.4f"
  % (float(nu(aN)), (float(nu(aN)) - NK) / NK, -(1 - w_aus) * Q_KOSM / NK))

open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
