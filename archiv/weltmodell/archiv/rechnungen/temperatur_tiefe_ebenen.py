# -*- coding: utf-8 -*-
"""
Temperatur und Tiefe auf den Ebenen über dem Stern: Gilt kT = f·δ·m̄c² auch für heißes Gas in Galaxienhalos,
Gruppen und Haufen, wenn δ die Tiefe nach der Regel ist (auf Bahnen δ = a·r/c² = v²/c², a = ν(a_N)·a_N)?
Rechnet nur mit ftm_kern.  [G] = Messwert aus dem Gedächtnis, nicht nachgeprüft;  [A] = Annahme.
Blickrichtung: T lokal gemessen (L, Spektrum am Ort, eigene Uhr); a(r) ist die örtliche Regel im Gesamtfeld (L);
δ = a·r/c² ist die wirksame Tiefe der Bahn am Ort (Δs_eff = ν·Δs, δ_eff = 3Δs_eff/r). Keine zweite Dichte:
die Regel ändert die Bahn, nicht den Stoff.
Teile: 1 Form   2 M_b–T (Gruppen, Haufen)   3 Elliptische Galaxien   4 Tiefenleiter   5 Gedächtnis (Kühlzeit)   6 Handproben
"""
import os, sys, math
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

out = []
def P(s=""):
    out.append(s)

m_u = 1.66053906660e-27
eV = e_el; keV = 1e3 * eV; MeV = 1e6 * eV
MU = 0.60            # [A] mittleres Teilchengewicht des voll ionisierten Plasmas (X = 0,72, Z = 0,3 Z☉)
MU_E = 1.14          # [A] Masse je Elektron in m_p
mupc2 = MU * m_p * c**2
kpc_ = 1e3 * pc
h70 = 0.70 / h_      # Umrechnung h70⁻¹-Massen auf H0 des Kerns

def rho_c(z=0.0):
    H = H0_si() * E_hintergrund(z)
    return 3 * H * H / (8 * math.pi * G)

def aN_aus_a(a):
    """Umkehr der Regel: a = ν(a_N)·a_N  →  a_N."""
    return brentq(lambda x: float(a_regel(x)) - a, a * 1e-8, a)

P("TEMPERATUR UND TIEFE AUF DEN EBENEN")
P("kT = f·δ·m̄c²,  δ = a·r/c² (Tiefe der Bahn nach der Regel),  m̄ = μ·m_p,  μ = %.2f [A]" % MU)
P()

# ============================================================== 1 Form
P("1) FORM FÜR GAS IN a(r) = ν(a_N)·a_N")
P("   Hydrostatik am Ort: d(nkT)/dr = −ρ·a  ⇒  kT(r) = μm_p·a·r/S,  S = −d ln ρ/d ln r − d ln T/d ln r")
P("   ⇒  kT = f·δ·m̄c²,  δ = a·r/c² = v_c²/c²,  f = 1/S  (örtlich).  Integral: kT(r) = (μm_p/ρ(r))·∫_r^∞ ρ·a dr'.")
P("   Im Stern dieselbe Bilanz, dort über das Potential integriert: f = θ_n = Δδ_innen/((n+1)δ).")
P("   Die Regel geht nur über δ ein: δ_eff = ν·δ_N; im Tiefbereich a·r = √(GM_b a₀) = v_flach² fest ⇒ δ fest ⇒ Gas isotherm.")
# Lane-Emden zur Kontrolle der Sternwerte
def lane_emden(n):
    def f(x, y):
        return [y[1], -max(y[0], 0) ** n - 2 * y[1] / x]
    x0 = 1e-6
    ev = lambda x, y: y[0]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, [x0, 50], [1 - x0**2 / 6, -x0 / 3], events=ev, rtol=1e-11, atol=1e-13)
    x1 = s.t_events[0][0]; th1 = s.y_events[0][0][1]
    w = -x1**2 * th1
    # Δδ_innen/δ = (Φ(R)−Φ(0))/(GM/R) = (n+1)·x1/w ... f = θ_n = x1/((n+1)·w)·(n+1) → f = x1·θ'-Form
    f_ = 1.0 / ((n + 1) * (-x1 * th1))
    return x1, w, f_
for n in (1.5, 3.0):
    x1, w, f_ = lane_emden(n)
    P("   Stern-Kontrolle Polytrope n = %.1f: ξ₁ = %.4f, f = 1/((n+1)·ξ₁|θ'(ξ₁)|) = %.4f" % (n, x1, f_))
P("   f für Gas (örtlich, f = 1/S):")
P("     singuläre isotherme Kugel ρ ∝ r⁻²: S = 2, f = 0,500")
for beta in (0.5, 2 / 3, 0.8):
    fs = ["x=%g: %.3f" % (x, (1 + x * x) / (3 * beta * x * x)) for x in (1, 3, 10)]
    P("     Beta-Profil β = %.3f (isotherm): f(x = r/r_c) = (1+x²)/(3βx²):  %s;  x→∞: %.3f" % (beta, ", ".join(fs), 1 / (3 * beta)))
# numerische Kontrolle: Beta-Profil im reinen Tiefbereich a = v²/r, Integralform
beta = 2 / 3
rho_b = lambda x: (1 + x * x) ** (-1.5 * beta)
for x in (1.0, 3.0, 10.0):
    I = quad(lambda t: rho_b(t) / t, x, np.inf, limit=200)[0] / rho_b(x)
    P("     Kontrolle Integralform, β = 2/3, Tiefbereich: kT/(μm_p v_flach²) bei x = %4.0f: %.4f  (Hydrostatik ohne Randdruck)" % (x, I))
# globale Virialform im Tiefbereich: Plummer-Kugel, Jeans isotrop
def virial_tief(b=1.0, M=1.0, Ga0=1.0):
    rho = lambda r: 3 * M / (4 * math.pi * b**3) * (1 + r * r / b / b) ** -2.5
    Mr = lambda r: M * r**3 / (r * r + b * b) ** 1.5
    a = lambda r: math.sqrt(Ga0 * Mr(r)) / r
    sig2 = lambda r: quad(lambda t: rho(t) * a(t), r, np.inf, limit=200)[0] / rho(r)
    num = quad(lambda r: 4 * math.pi * r * r * rho(r) * sig2(r), 0, 200 * b, limit=400, points=[1, 10])[0]
    return num / M / math.sqrt(Ga0 * M)
vt = virial_tief()
P("   Global (isoliert, ganz im Tiefbereich, jedes Profil): ⟨σ²⟩_M = ⟨kT⟩_M/μm_p = (2/9)·√(GM_b a₀)  ⇔  σ⁴ = (4/81)GM_b a₀")
P("     Plummer numerisch: ⟨σ_r²⟩_M/√(GMa₀) = %.4f  (2/9 = %.4f);  f_global = 2/9 bezogen auf δ_flach = v_flach²/c²" % (vt, 2 / 9))
P("   Blickrichtung: T lokal (L); δ = a·r/c² am Ort der Bahn (L), dieselbe Zahl, die der ferne Beobachter (G) als Taktunterschied")
P("   3Δs/r nur im Newton-Bereich liest; im Tiefbereich ist δ = 3Δs_eff/r mit Δs_eff = ν·Δs. Kein Stoff kommt hinzu.")
P()

# ============================================================== 2 M_b–T
P("2) GRUPPEN UND HAUFEN: M_b gegen T")
P("   Tiefbereich (Form): v_flach⁴ = GM_b a₀, kT = μm_p v_flach²/S  ⇒  M_b = (S·kT/μm_p)²/(G a₀)   ∝ T²")
P("   Newton: M = S·kT·r/(Gμm_p), r_Δ ∝ M^{1/3}  ⇒  M ∝ T^{3/2}")
Tl = [0.5, 1, 2, 5, 10]
P("   T[keV]   M_b,tief(S=2)[M☉]   M_b,tief(S=2,4)   M_b,tief global 2/9 (⟨kT⟩_M)")
for T in Tl:
    u = T * keV / (MU * m_p)
    Mb2 = (2 * u) ** 2 / (G * A0) / M_sun
    Mb24 = (2.4 * u) ** 2 / (G * A0) / M_sun
    Mbg = (u / (2 / 9)) ** 2 / (G * A0) / M_sun
    P("   %5.1f     %.3e          %.3e        %.3e" % (T, Mb2, Mb24, Mbg))
P("   (S = 2: β = 2/3 außen; S = 2,4: gemessenes Dichte- plus Temperaturgefälle bei r₅₀₀ [G, typisch 2–3])")
P()
# Messwert-Skalen [G]
P("   [G] M₅₀₀–T: Arnaud, Pointecouteau, Pratt 2005: M₅₀₀ = 3,84·10¹⁴ h70⁻¹ M☉ (T/5 keV)^1,71 (alle), 4,10·10¹⁴ (T/5)^1,49 (T > 3,5 keV);")
P("       Vikhlinin 2006: 3,32·10¹⁴ h⁻¹ (T/5)^1,47.  [G] f_gas,500 (Vikhlinin 2006, Sun 2009, Lovisari 2015): 0,04 / 0,06 / 0,08 / 0,11 / 0,125;")
P("       [G] f_*,500 inkl. BCG+ICL (Lin 2003, Gonzalez 2013): 0,05 / 0,035 / 0,025 / 0,015 / 0,010  (je ±30–50 %)")
P("       [G] gemessene Steigungen: M_gas ∝ T^1,8–2,1 (Arnaud 2007, Sun 2009, Lovisari 2015), M_tot ∝ T^1,5–1,7")
fgas = {0.5: 0.04, 1: 0.06, 2: 0.08, 5: 0.11, 10: 0.125}
fst = {0.5: 0.05, 1: 0.035, 2: 0.025, 5: 0.015, 10: 0.010}
def M500_arnaud(T):
    return (3.84e14 * (T / 5) ** 1.71 if T < 3.5 else 4.10e14 * (T / 5) ** 1.49) * h70
def M500_vikh(T):
    return 3.32e14 / h_ * (T / 5) ** 1.47
def r_delta(M, D=500, z=0.0):
    return (3 * M * M_sun / (4 * math.pi * D * rho_c(z))) ** (1 / 3)
P()
P("   Lesung am Rand r₅₀₀ (a_tot = GM₅₀₀/r₅₀₀² aus der Hydrostatik [G]), Regel umgekehrt → nötige Baryonen M_b,Regel:")
P("   T    Quelle    M₅₀₀[M☉]   r₅₀₀[kpc]  a_tot/a₀  a_N/a₀    ν     S_impl  M_b,Regel   M_b,gem[G]  Regel/gem  M_b,tief(S_impl)")
rows = []
for T in Tl:
    for name, Mf in (("Arnaud", M500_arnaud), ("Vikhl.", M500_vikh)):
        M5 = Mf(T); r5 = r_delta(M5)
        at = G * M5 * M_sun / r5**2
        aN = aN_aus_a(at)
        Mreg = aN * r5**2 / G / M_sun
        Mgem = (fgas[T] + fst[T]) * M5
        S_impl = G * M5 * M_sun * MU * m_p / (T * keV * r5)
        Mtief = (S_impl * T * keV / (MU * m_p)) ** 2 / (G * A0) / M_sun
        P("   %4.1f %-8s %.3e  %7.0f    %.3f    %.4f  %5.2f  %5.2f   %.3e   %.3e   %5.2f     %.3e" %
          (T, name, M5, r5 / kpc_, at / A0, aN / A0, float(nu(aN)), S_impl, Mreg, Mgem, Mreg / Mgem, Mtief))
        if name == "Arnaud":
            rows.append((T, Mreg, Mgem, M5))
lT = np.log([r[0] for r in rows])
p_reg = np.polyfit(lT, np.log([r[1] for r in rows]), 1)[0]
p_gem = np.polyfit(lT, np.log([r[2] for r in rows]), 1)[0]
p_tot = np.polyfit(lT, np.log([r[3] for r in rows]), 1)[0]
P("   Steigungen 0,5–10 keV (Arnaud): M_b,Regel ∝ T^%.2f,  M_b,gem ∝ T^%.2f [G],  M₅₀₀ ∝ T^%.2f [G]" % (p_reg, p_gem, p_tot))
P("   Die Regel verlangt die Steigung 2 (Tiefbereich, a_N/a₀ ≈ 0,03–0,1 bei r₅₀₀); gemessen ist ~2 für die Baryonen.")
# Grenzmasse M_Q in Temperatur
c_nfw = 4.0
mnfw = lambda x: math.log(1 + x) - x / (1 + x)
# M500/M200 für NFW c = 4
def m500_zu_m200(cc):
    # finde r500/r_s: mittlere Dichte 500 ρc
    f = lambda y: mnfw(y) / y**3 - (500 / 200) * mnfw(cc) / cc**3
    y5 = brentq(f, 0.01, cc)
    return mnfw(y5) / mnfw(cc)
r52 = m500_zu_m200(c_nfw)
MQ200 = 2.8e14
MQ500 = r52 * MQ200
TQ = brentq(lambda T: M500_arnaud(T) - MQ500, 0.1, 20)
P("   Grenzmasse der Mode (Modell): M_Q ≈ 2,8·10¹⁴ M☉ (M₂₀₀) → M₅₀₀ = %.3f·M₂₀₀ (NFW c = 4) = %.2e M☉ → T_Q ≈ %.1f keV (Arnaud [G])" % (r52, MQ500, TQ))
P("   Nach dem Modell wirkt unter T_Q nur der Operator: Dort müsste die Regel allein (Spalte Regel/gem) ≈ 1 geben.")
P("   Gerechnet: 0,5–2 keV (unter T_Q) 1,3–2,0, 5–10 keV (über T_Q) 2,3–2,7. Die Lücke wächst mit T (wie der Rest der Mode),")
P("   verschwindet unter T_Q aber nicht. Mit hydrostatischem Defizit [G, 10–30 %] wüchse sie noch; mit Baryonen außerhalb")
P("   der Zählung (warmes Gas, f_* und f_gas der Gruppen je ±30–50 % [G]) schrumpfte sie. Bei 0,5–1 keV im Fehler, bei 2 keV nicht.")
P()

# Haufenprofil: T(r) aus Hydrostatik mit Newton, Regel und Gesamtwirkung (Regel + Mode), Baryonen wie in der Rechnung zur Schließungsregel
P("   Haufenprofil (Baryonen [A]: Gas β = 2/3, r_c = 150 kpc, 10¹⁴ M☉ in 1 Mpc, bis 3 Mpc; BCG 1,5·10¹² (Hernquist 20 kpc);")
P("   Galaxien 10¹³ (Hernquist 300 kpc); Gesamtwirkung Regel + Mode vertreten durch NFW M₂₀₀ = 10¹⁵ M☉, c = 4, z = 0,3,")
P("   die das Modell mit dem Rest der Mode auf 1–5 % trifft)")
rc_g = 150 * kpc_
rho_shape = lambda r: (1 + (r / rc_g) ** 2) ** -1.0
norm = 1e14 * M_sun / quad(lambda r: 4 * math.pi * r * r * rho_shape(r), 0, 1000 * kpc_)[0]
rg = np.logspace(math.log10(1 * kpc_), math.log10(3000 * kpc_), 3000)
dM = 4 * math.pi * rg**2 * norm * rho_shape(rg)
Mg = np.concatenate([[0], np.cumsum(0.5 * (dM[1:] + dM[:-1]) * np.diff(rg))]) + norm * 4 / 3 * math.pi * rg[0] ** 3
Mgas = lambda r: np.interp(r, rg, Mg)
Mb_h = lambda r: Mgas(r) + hernquist_masse(1.5e12 * M_sun, 20 * kpc_, r) + hernquist_masse(1e13 * M_sun, 300 * kpc_, r)
z_h = 0.3
r200 = r_delta(1e15, 200, z_h); rs = r200 / c_nfw
Mnfw = lambda r: 1e15 * M_sun * mnfw(r / rs) / mnfw(c_nfw)
aNh = lambda r: G * Mb_h(r) / r**2
acc = {"Newton": lambda r: aNh(r), "Regel": lambda r: float(a_regel(aNh(r))), "Regel+Mode": lambda r: G * Mnfw(r) / r**2}
def T_profil(afun, r, Rout=3000 * kpc_):
    I = quad(lambda t: rho_shape(t) * afun(t), r, Rout, limit=300)[0]
    return MU * m_p * I / rho_shape(r) / keV
P("   r[kpc]   a_N/a₀    ν     kT_Newton  kT_Regel  kT_Regel+Mode [keV]   δ_Regel=a r/c²  f = kT/(δ m̄c²) (Regel)")
Tp = {}
for rk in (50, 100, 250, 500, 1000, 1500):
    r = rk * kpc_
    Ts = {k: T_profil(f, r) for k, f in acc.items()}
    dR = acc["Regel"](r) * r / c**2
    Tp[rk] = Ts
    P("   %5d   %.4f  %5.2f    %6.2f     %6.2f      %6.2f            %.3e        %.3f" %
      (rk, aNh(r) / A0, float(nu(aNh(r))), Ts["Newton"], Ts["Regel"], Ts["Regel+Mode"], dR, Ts["Regel"] * keV / (dR * mupc2)))
# emissionsgewichtet 0,15–1 r500
M5h = brentq(lambda M: M * M_sun - Mnfw(r_delta(M, 500, z_h)), 1e13, 1e16)
r5h = r_delta(M5h, 500, z_h)
def T_ew(afun):
    rr = np.linspace(0.15 * r5h, r5h, 40)
    w = np.array([rho_shape(x) ** 2 * x * x for x in rr])
    t = np.array([T_profil(afun, x) for x in rr])
    return float(np.sum(w * t) / np.sum(w))
Tew = {k: T_ew(f) for k, f in acc.items()}
P("   emissionsgewichtet 0,15–1 r₅₀₀ (r₅₀₀ = %.0f kpc, M₅₀₀ = %.2e): Newton %.2f, Regel %.2f, Regel+Mode %.2f keV" %
  (r5h / kpc_, M5h, Tew["Newton"], Tew["Regel"], Tew["Regel+Mode"]))
P("   [G] ein Haufen mit M₅₀₀ ≈ 7·10¹⁴ M☉ hat T ≈ 7–8 keV (Arnaud 2005: %.1f keV).  Regel allein: T um Faktor %.2f zu kalt,"
  % (brentq(lambda T: M500_arnaud(T) - M5h, 0.5, 30), Tew["Regel+Mode"] / Tew["Regel"]))
P("   grob gleichbedeutend M_b um den Faktor %.2f (T² im Tiefbereich) zu klein: die bekannte Lücke; der Rest der Mode schließt sie."
  % ((Tew["Regel+Mode"] / Tew["Regel"]) ** 2))
P()

# ============================================================== 3 Elliptische Galaxien
P("3) ELLIPTISCHE GALAXIEN: T des heißen Halos gegen σ der Sterne")
P("   Gas und Sterne stehen im selben a(r): σ_r² = (1/ρ_*)∫ρ_* a dr, kT/μm_p = (1/ρ_g)∫ρ_g a dr.")
P("   Örtlich: σ² ≈ a·r/S_*, kT/μm_p = a·r/S_g  ⇒  β_spec = μm_pσ²/kT ≈ S_g/S_* (Verhältnis der Profile), wenn a·r fest ist.")
rcg, bg = 1.0 * kpc_, 0.5     # [A] Gas β = 0,5, r_c = 1 kpc; Gasmasse vernachlässigt (≲ 5 % der Sterne in 50 kpc)
rhoG = lambda r: (1 + (r / rcg) ** 2) ** (-1.5 * bg)
Rmax = 500 * kpc_
res3 = {}
for Ms, Re in ((1.5e11, 6.0), (4e11, 10.0)):      # [A] Sternmasse, effektiver Radius
    aH = Re / 1.8153 * kpc_
    Mstar = lambda r, Ms=Ms, aH=aH: hernquist_masse(Ms * M_sun, aH, r)
    rhoS = lambda r, Ms=Ms, aH=aH: Ms * M_sun * aH / (2 * math.pi * r * (r + aH) ** 3)
    aE = {"Newton": lambda r, Mstar=Mstar: G * Mstar(r) / r**2,
          "Regel": lambda r, Mstar=Mstar: float(a_regel(G * Mstar(r) / r**2))}
    def sig2(afun, r):
        return quad(lambda t: rhoS(t) * afun(t), r, Rmax, limit=300, points=[aH, 10 * aH])[0] / rhoS(r)
    def kTg(afun, r):
        return MU * m_p * quad(lambda t: rhoG(t) * afun(t), r, Rmax, limit=300, points=[aH, 10 * aH])[0] / rhoG(r) / keV
    def sig_ap(afun):
        # 3D-Mittel im Radius R_e/2, gewichtet mit ρ_* (Näherung der zentralen Blende)
        rr = np.logspace(math.log10(0.01 * aH), math.log10(0.5 * Re * kpc_), 60)
        w = np.array([rhoS(x) * x**3 for x in rr]); s_ = np.array([sig2(afun, x) for x in rr])
        return math.sqrt(np.sum(w * s_) / np.sum(w))
    def kT_ew(afun, r1, r2):
        rr = np.logspace(math.log10(r1), math.log10(r2), 40)
        w = np.array([rhoG(x) ** 2 * x**3 for x in rr]); t = np.array([kTg(afun, x) for x in rr])
        return float(np.sum(w * t) / np.sum(w))
    P("   Modell [A]: Hernquist M_* = %.1e M☉, R_e = %.0f kpc (v_flach = %.0f km/s); Gas β = 0,5, r_c = 1 kpc; T emissionsgewichtet 1–5 R_e"
      % (Ms, Re, v_flach(Ms * M_sun) / 1e3))
    P("   Lesart   σ_*(<R_e/2)[km/s]   kT(1–5R_e)[keV]  kT(R_e)  kT(5R_e)  kT(20R_e)   β_spec = μm_pσ²/kT")
    for k, f in aE.items():
        sg = sig_ap(f); kt = kT_ew(f, Re * kpc_, 5 * Re * kpc_)
        res3[(Ms, k)] = (sg, kt)
        P("   %-7s      %6.1f              %.3f         %.3f    %.3f     %.3f        %.3f" %
          (k, sg / 1e3, kt, kTg(f, Re * kpc_), kTg(f, 5 * Re * kpc_), kTg(f, 20 * Re * kpc_), MU * m_p * sg**2 / (kt * keV)))
    for rk in (1, 5, 20):
        r = rk * Re * kpc_; aN = G * Mstar(r) / r**2
        P("     bei %2d R_e: a_N/a₀ = %.3f, ν = %.2f" % (rk, aN / A0, float(nu(aN))))
    P("   Regel/Newton: kT(1–5 R_e) × %.2f, σ_* × %.3f" % (res3[(Ms, "Regel")][1] / res3[(Ms, "Newton")][1], res3[(Ms, "Regel")][0] / res3[(Ms, "Newton")][0]))
    if Ms == 1.5e11:
        aE1, kTg1, Re1 = aE, kTg, Re
P("   Örtlich kürzt sich a in β_spec nur, wo a·r fest ist (Tiefbereich, isotherm): β_spec → S_g/S_*. Newton: a·r fällt außen,")
P("   T fällt nach außen steil und β_spec > 1. Die Regel gibt flaches T und β_spec ≈ 0,6.")
P("   [G] gemessen (O'Sullivan, Ponman, Collins 2003; Davis & White 1996): σ ≈ 250 km/s ↔ kT_X ≈ 0,6–0,9 keV, β_spec ≈ 0,5–1;")
P("   T flach oder nach außen steigend (NGC 720, NGC 1521 [G]).")
bs = res3[(4e11, "Regel")]; bsp = MU * m_p * bs[0]**2 / (bs[1] * keV)
P("   Mit β_spec(Regel) = %.3f: σ = 250 km/s → kT = μm_pσ²/β_spec = %.2f keV (μm_pσ² = %.3f keV)." %
  (bsp, MU * m_p * 250e3**2 / keV / bsp, MU * m_p * 250e3**2 / keV))
P("   Faber-Jackson [G]: σ ≈ 250 km/s gehört zu M_* ≈ 2–5·10¹¹ M☉.")
P()

# ============================================================== 4 Tiefenleiter
P("4) TIEFENLEITER ÜBER DIE EBENEN: δ und kT/(m̄c²)")
m_n = 1.67492749804e-27
Rsun = R_sun
lad = []
dsun = delta(M_sun, Rsun) if False else 3 * ds(M_sun) / Rsun
lad.append(("Sonnenkern (Hauptreihe)", dsun, 1.57e7, 0.85 * m_u, "δ = 3Δs/R; T_c [G]; nicht entartet"))
lad.append(("Weißer Zwerg 0,6 M☉ (1e7 K)", 3 * ds(0.6 * M_sun) / (0.0125 * Rsun), 1e7, 12 / 7 * m_u, "entartet: Zählregel hält die Tiefe"))
lad.append(("Neutronenstern alt (1e8 K)", 3 * ds(1.4 * M_sun) / 12e3, 1e8, m_n, "entartet (Neutronen)"))
lad.append(("Proto-Neutronenstern (3e11 K)", 3 * ds(1.4 * M_sun) / 20e3, 3e11, m_n, "[G] ~30 MeV, R ≈ 20 km"))
aRe = aE1["Regel"](3 * Re1 * kpc_)
lad.append(("Elliptische 1,5e11 (Halo, 3R_e)", aRe * 3 * Re1 * kpc_ / c**2, kTg1(aE1["Regel"], 3 * Re1 * kpc_) * keV / kB, MU * m_p, "δ = a r/c² nach der Regel"))
for T in (1.0, 5.0):
    M5 = M500_arnaud(T); r5 = r_delta(M5)
    at = G * M5 * M_sun / r5**2
    lad.append(("%s %.0f keV (r₅₀₀)" % ("Gruppe" if T < 2 else "Haufen", T), at * r5 / c**2, T * keV / kB, MU * m_p, "δ = a_tot r/c² [G-Skala]"))
rH = Tew["Regel+Mode"]
lad.append(("Haufen 1e15 (Profil, 0,15–1 r₅₀₀)", G * Mnfw(0.5 * r5h) / (0.5 * r5h) / c**2, rH * keV / kB, MU * m_p, "δ bei r₅₀₀/2, Regel+Mode"))
P("   Ebene                                δ             T [K]       kT/(m̄c²)     f = kT/(δ m̄c²)   Bemerkung")
for n_, d, T, mb, bem in lad:
    x = kB * T / (mb * c**2)
    P("   %-35s  %.3e    %.3e    %.3e    %.3g            %s" % (n_, d, T, x, x / d, bem))
P("   Nicht entartete Ebenen (Stern, Halo, Gruppe, Haufen): f = 0,3–1,4; entartete: f ≪ 1 (Zählregel statt Wärme).")
P("   δ der Gas-Ebenen (10⁻⁶ bis 3·10⁻⁵) liegt im selben Band wie die Brennleiter der Sterne (7·10⁻⁷ pp … 1,7·10⁻⁴ Si):")
P("   Gruppen in der Tiefe des Wasserstoffbrennens, reiche Haufen in der des Heliumbrennens (Tatsache der Tabelle, nicht hergeleitet).")
P()
P("   Kosmische Ebene (Blickrichtung G): Im homogenen Hintergrund gibt es kein Gefälle, a_N = 0, keine Rückkehr (Γ = 0), keine Bahn;")
P("   δ = a·r/c² ist dort null. Die Temperatur des Hintergrunds folgt dem Atmen, T ∝ 1/a (Rotverschiebung), nicht einer Tiefe.")
Trec = 3000.0
P("   Rekombination: kT/(m_p c²) = %.2e bei T = 3000 K; Tiefe der Störungen |Φ|/c² ~ 10⁻⁵ [G]: gibt nur ΔT/T = Φ/(3c²) (Sachs-Wolfe)," % (kB * Trec / (m_p * c**2)))
P("   also f nicht von der Ordnung 1: Die globale Ebene hat keine Tiefe in diesem Sinn; die Tiefe tritt dort nur als Störung auf.")
P("   Horizonttemperatur (04, global): T = ħc/(12πkΔs_H) enthält ħ, nicht m̄c²·δ; sie ist eine andere Form (Quantenlänge statt Masse).")
P()

# ============================================================== 5 Gedächtnis: Kühlzeit
P("5) GEDÄCHTNIS: Kühlzeit gegen Alter")
P("   t_cool = (3/2)(n_e+n_i)kT/(n_e n_H Λ(T)) = 2,875·kT/(n_H Λ),  n_e = 1,2 n_H, n_ges = 2,3 n_H [A]")
P("   Λ(T) [G] Tozzi & Norman 2001 (Z = 0,3 Z☉): Λ = (8,6·10⁻³ T^−1,7 + 5,8·10⁻² T^0,5 + 6,3·10⁻²)·10⁻²² erg cm³/s, T in keV")
def Lam(T):
    return (8.6e-3 * T ** -1.7 + 5.8e-2 * T**0.5 + 6.3e-2) * 1e-22 * 1e-7 * 1e-6   # → W m³
def tcool(ne_cm3, T):
    nH = ne_cm3 / 1.2 * 1e6
    return 2.875 * T * keV / (nH * Lam(T)) / Gyr
def ne_r500(T, fg, Sg=2.7):
    # örtliche Dichte am Rand: ρ(r₅₀₀) = (S_g/3)·f_gas·500ρ_c  (Potenzprofil, S_g [G] 2,5–3)
    return Sg / 3 * fg * 500 * rho_c() / (MU_E * m_p) / 1e6
P("   Ort                                   n_e[cm⁻³]      kT[keV]   t_cool[Gyr]   K = kT n_e^−2/3 [keV cm²]")
cases = [("Kühlkern-Zentrum (10 kpc) [G]", 0.05, 3.0), ("Haufenkern ohne Kühlkern [G]", 5e-3, 7.0),
         ("Gruppenzentrum [G]", 0.03, 0.8), ("Elliptischer Halo, 5 R_e [G]", 1e-3, 0.7),
         ("Gruppe 1 keV, r₅₀₀ (aus f_gas)", ne_r500(1, fgas[1]), 1.0), ("Haufen 5 keV, r₅₀₀ (aus f_gas)", ne_r500(5, fgas[5]), 5.0)]
for n_, ne, T in cases:
    P("   %-36s  %.2e      %5.1f      %8.3g        %.3g" % (n_, ne, T, tcool(ne, T), T * ne ** (-2 / 3)))
z_c = 1.0
def alter(z):
    return quad(lambda zz: 1 / ((1 + zz) * E_hintergrund(zz)), 0, z)[0] / H0_si() / Gyr
cs = math.sqrt(5 / 3 * 5 * keV / (MU * m_p))
r5_5 = r_delta(M500_arnaud(5))
P("   Alter seit Bildung bei z = 1 bzw. 0,5: %.1f bzw. %.1f Gyr.  Schallzeit bis r₅₀₀ (5 keV): %.2f Gyr (c_s = %.0f km/s)" %
  (alter(1.0), alter(0.5), r5_5 / cs / Gyr, cs / 1e3))
tb = 0.05; ne_eq = brentq(lambda n: tcool(n, 3.0) - alter(1.0), 1e-4, 1)
P("   Kühlzeit = Alter (%.1f Gyr) bei 3 keV für n_e = %.1e cm⁻³: innen kühlt das Gas, außen nicht." % (alter(1.0), ne_eq))
# Kühlkern-Temperatur gegen örtliche Tiefe: v_c des BCG-Trichters bei 10–30 kpc (Regel), S ≈ 1 (flacher Kern)
P("   Kühlkern: Temperatur gegen örtliche Tiefe (BCG 1,5·10¹² M☉ + Haufen wie oben, Regel+Mode):")
for rk in (10, 30):
    r = rk * kpc_
    aT = G * (Mnfw(r) + hernquist_masse(1.5e12 * M_sun, 20 * kpc_, r)) / r**2
    vc2 = aT * r
    P("     r = %2d kpc: v_c = %.0f km/s, δ = %.2e, μm_p v_c² = %.2f keV;  kT = f·δ·m̄c² mit f = 1/S, S = 1–2: %.2f–%.2f keV"
      % (rk, math.sqrt(vc2) / 1e3, vc2 / c**2, MU * m_p * vc2 / keV, MU * m_p * vc2 / keV / 2, MU * m_p * vc2 / keV))
P("   [G] Kühlkerne: T fällt innen auf 0,3–0,5 T_mittel (2–4 keV bei 7 keV), nicht tiefer; K₀ ≈ 10–30 keV cm² (Cavagnolo 2009).")
P()

# ============================================================== 6 Handproben
P("6) UNABHÄNGIGE PROBEN (von Hand, anderer Weg)")
u1 = 1 * keV / (MU * m_p)
P("   (a) M_b,tief(1 keV, S = 2) von Hand: kT/μm_p = 1,6022·10⁻¹⁶/(0,6·1,6726·10⁻²⁷) = 1,5965·10¹¹ m²/s²; ×2 → 3,193·10¹¹, ² = 1,0195·10²³;")
P("       G a₀ = 6,6743·10⁻¹¹·1,1·10⁻¹⁰ = 7,3417·10⁻²¹ → 1,3886·10⁴³ kg = 6,98·10¹² M☉.  Skript: %.3e M☉" % ((2 * u1) ** 2 / (G * A0) / M_sun))
P("   (b) T im Tiefbereich aus M_b: 10¹⁴ M☉: v_f⁴ = G·1,989·10⁴⁴·1,1·10⁻¹⁰ = 1,4602·10²⁴ → v_f² = 1,2084·10¹² → v_f = 1099 km/s;")
vf2 = math.sqrt(G * 1e14 * M_sun * A0)
P("       kT = μm_p v_f²/2 = 0,6·1,6726·10⁻²⁷·1,2084·10¹²/2 = 6,063·10⁻¹⁶ J = 3,78 keV.  Skript: %.3f keV" % (MU * m_p * vf2 / 2 / keV))
P("   (c) t_cool(n_e = 0,05, 3 keV) von Hand: Λ = (8,6·10⁻³·0,1546 + 5,8·10⁻²·1,732 + 0,063)·10⁻²² = 1,648·10⁻²³ erg cm³/s;")
P("       n_H = 0,04167; 2,875·4,806·10⁻⁹ erg/(0,04167·1,648·10⁻²³) = 2,012·10¹⁶ s = 0,638 Gyr.  Skript: %.3f Gyr" % tcool(0.05, 3.0))
P("   (d) δ eines 5-keV-Haufens mit f = ½: 5·10³ eV/(0,5·0,6·938,27·10⁶ eV) = 1,776·10⁻⁵.  Skript: %.4e" % (5 * keV / (0.5 * mupc2)))
P("   (e) Elliptische, Tiefbereich-Grenze: v_f = (G·1,5·10¹¹·1,989·10³⁰·1,1·10⁻¹⁰)^¼ = (2,190·10²¹)^¼ = 216,3 km/s;")
P("       kT(S = 1,5) = 0,6·1,6726·10⁻²⁷·(2,163·10⁵)²/1,5 = 3,130·10⁻¹⁷ J = 0,195 keV.  Skript v_flach: %.1f km/s" % (v_flach(Ms * M_sun) / 1e3))

txt = "\n".join(out)
print(txt)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "temperatur_tiefe_ebenen.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
