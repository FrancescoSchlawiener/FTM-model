# -*- coding: utf-8 -*-
"""
Trichterkette mit Gewicht w: Welche Regel gilt auf welcher Ebene zwischen Lokaler Gruppe, Virgo, Großem Attraktor, Shapley?
Frühere Kette (trichterkette_sonne): örtliche Regel ν(a_N) auf Mpc-Skalen → 1 500–5 400 km/s.
Hier: Die örtliche Regel wirkt nur in umgekehrten Bereichen (θ̃ → 0, w → 0); mitatmend (w ≈ 1) wirkt die Mode linear wie
druckfreier Staub mit 1 + Q (Standardmodell mit Ω_m = (1 + Q)Ω_b).
Abschnitte
  1  Kugelschale (Top-Hat) im Kern-Hintergrund: θ, θ̃, w gegen den heutigen mittleren Kontrast 1 + Δ; Umkehr heute Δ_ta
  2  Bereiche um uns, Virgo, Großen Attraktor, Shapley (1–200 Mpc): 1 + Δ, θ/3H₀, θ̃/θ_c, w, umgekehrt?
  3  Pekuliargeschwindigkeit der Lokalen Gruppe linear (f H₀/4π ∫δ r̂/r² dV) und kugelig nichtlinear für Virgo; Richtung
  4  Anteile: linear 1/r² (Überschussmasse) gegen frühere 1/r-Lesart; frühere Zahlen
  5  Wechsel örtlich → mitatmend: Umkehrradius gegen R_Q; Zeitargument (Lynden-Bell) mit Mode (1 + Q)M_b gegen ν(a_N)
  6  Unabhängige Prüfung von Hand
Blickrichtung: θ, w [G am Ort L], a_N, ν [L], θ_c = a₀/c [I].  [G] Wert aus dem Gedächtnis, ungeprüft; [A] Annahme.
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []
def P(s=""):
    print(s); OUT.append(s)
def H(t):
    P(); P("=" * 110); P(t); P("=" * 110)

G, c, Mpc, kpc, Msun = K.G, K.c, K.Mpc, K.kpc, K.M_sun
H0 = K.H0_si()
Om, Or = K.OMEGA_M, K.OMEGA_R
OL = 1 - Om - Or
Q, qk = K.Q_KOSM, K.q_KOSM
thc = K.A0 / c                       # θ_c [I]
rho_m = K.rho_m_mitbewegt()          # kg/m³
rho_m_Msun_Mpc3 = rho_m * Mpc**3 / Msun
fb = 1 / (1 + Q)

def einheit(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)])
def lb(v):
    r = np.linalg.norm(v)
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2] / r))
def winkel(u, v):
    return math.degrees(math.acos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))

CMB = einheit(276.0, 30.0); V_CMB = 620.0            # [G] Dipol der LG

# ------------------------------------------------------------------ 1 Top-Hat
def tophat(eps, Om_=Om, Or_=Or, OL_=None, a_i=1e-2, dicht=4000):
    """Schale mit Massenüberschuss eps, Anfang im Wachstumsmodus (Materiezeit). Zeit in 1/H₀.
    ä = −(Ω_m/2)a⁻² − Ω_r a⁻³ + Ω_Λ a;  ÿ = −(Ω_m/2)(1+eps)y⁻² − Ω_r a⁻⁴ y + Ω_Λ y.
    Rückgabe heute (a = 1): 1+Δ, θ/3H₀, θ̃/3H₀ (min der Geschichte, ≥ 0), umgekehrt?"""
    if OL_ is None:
        OL_ = 1 - Om_ - Or_
    E_i = math.sqrt(Om_ * a_i**-3 + Or_ * a_i**-4 + OL_)
    adot = a_i * E_i
    y0, ydot = a_i, adot * (1 - eps / 3)
    def rhs(t, s):
        a, ad, y, yd = s
        return [ad, -0.5 * Om_ / a**2 - Or_ / a**3 + OL_ * a,
                yd, -0.5 * Om_ * (1 + eps) / y**2 - Or_ / a**4 * y + OL_ * y]
    def heute(t, s): return s[0] - 1.0
    heute.terminal = True
    def kollaps(t, s): return s[2] - 1e-4 * a_i
    kollaps.terminal = True
    sol = solve_ivp(rhs, (0, 5), [a_i, adot, y0, ydot], events=[heute, kollaps], rtol=1e-10, atol=1e-13, dense_output=True)
    if sol.t_events[1].size:                                      # vor heute kollabiert
        return dict(D=np.inf, th=np.nan, thmin=0.0, um=True, t0=np.nan)
    t0 = sol.t_events[0][0]
    ts = np.linspace(sol.t[0], t0, dicht)
    a, ad, y, yd = sol.sol(ts)
    th = yd / y                                                 # θ/3H₀
    return dict(D=(1 + eps) * (1 / y[-1])**3, th=th[-1], thmin=max(th.min(), 0.0), um=bool(yd.min() <= 0), t0=t0)

def w_von(th3):
    """w aus θ̃/3H₀ (θ̃ = 3H₀·th3), Kern gewicht_w."""
    return K.gewicht_w(3 * H0 * th3)

H("1   Kugelschale im Kern-Hintergrund: θ̃ und w gegen den heutigen mittleren Kontrast 1 + Δ (Blick G am Ort L)")
P(f"   Ω_m = (1+Q)Ω_b = {Om:.4f}, Ω_r = {Or:.2e}, Ω_Λ = {OL:.4f};  H₀ = {H0:.4e} s⁻¹;  θ_c = a₀/c = {thc:.4e} s⁻¹ [I]")
P(f"   θ_c/3H₀ = {thc/(3*H0):.5f}  (3H₀τ_D = {3*H0*c/K.A0:.2f});  w = ½ bei θ̃ = θ_c, w = 0,9 bei θ̃ = 3θ_c, w = 0,1 bei θ̃ = θ_c/3")
P(f"   ρ̄_m = {rho_m:.4e} kg/m³ = {rho_m_Msun_Mpc3:.4e} M☉/Mpc³")
eps_grid = np.concatenate([np.linspace(-0.009, 0.0, 10), np.geomspace(1e-4, 0.05, 60)])
tab = []
for e in eps_grid:
    r = tophat(e)
    if np.isfinite(r["D"]) and not r["um"]:
        tab.append((r["D"], r["th"], r["thmin"], e))
tab = np.array(sorted(tab))
t0 = tophat(0.0)["t0"]
P(f"   Alter t₀ = {t0:.5f}/H₀ = {t0/H0/K.Gyr:.3f} Gyr;  Probe ε = 0: 1+Δ = {tophat(0.0)['D']:.6f}, θ/3H₀ = {tophat(0.0)['th']:.6f} (soll 1, E(0) = {K.E_hintergrund(0):.6f})")
e_ta = brentq(lambda e: tophat(e)["th"] if not tophat(e)["um"] else -1.0, 1e-3, 0.05, xtol=1e-9)
r_ta = tophat(e_ta)
D_ta = r_ta["D"]
P(f"   Umkehr heute (θ = 0): 1 + Δ_ta = {D_ta:.3f}")
eds = tophat(e_ta, Om_=1.0, Or_=0.0, OL_=0.0)
e_eds = brentq(lambda e: tophat(e, 1.0, 0.0, 0.0)["th"] if not tophat(e, 1.0, 0.0, 0.0)["um"] else -1.0, 1e-3, 0.05, xtol=1e-9)
P(f"   Probe Einstein-de Sitter: 1 + Δ_ta = {tophat(e_eds, 1.0, 0.0, 0.0)['D']:.4f}  (analytisch 9π²/16 = {9*math.pi**2/16:.4f})")
def D_bei_w(wz):
    thz = (thc / (3 * H0)) * math.sqrt(wz / (1 - wz))
    e = brentq(lambda e: tophat(e)["th"] - thz if not tophat(e)["um"] else -1.0, 0.0, e_ta, xtol=1e-10)
    return tophat(e)["D"], thz
P("   w heute     θ̃/3H₀     1 + Δ       (1+Δ)/(1+Δ_ta)   v_pec/(H₀R) = θ/3H₀ − 1")
for wz in (0.99, 0.9, 0.5, 0.1):
    D, thz = D_bei_w(wz)
    P(f"   {wz:6.3f}   {thz:9.5f}   {D:9.4f}    {D/D_ta:9.4f}        {thz-1:+.4f}")
P("   θ fällt in jeder noch expandierenden Schale monoton (θ̃ = θ heute); Tabelle 1 + Δ → θ/3H₀, w:")
for D in (1.0, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0):
    if D < tab[-1, 0]:
        th = np.interp(D, tab[:, 0], tab[:, 1]); thm = np.interp(D, tab[:, 0], tab[:, 2])
        f_lin = Om**0.55
        P(f"      1+Δ = {D:5.2f}: θ/3H₀ = {th:.4f} (linear 1 − fΔ/3 = {1 - f_lin*(D-1)/3:+.4f}), θ̃/θ_c = {3*H0*thm/thc:7.2f}, 1 − w = {1-w_von(thm):.2e}")
def lies(D):
    if D >= D_ta:
        return 0.0, 0.0, 0.0, True
    th = np.interp(D, tab[:, 0], tab[:, 1]); thm = np.interp(D, tab[:, 0], tab[:, 2])
    return th, thm, w_von(thm), False

# ------------------------------------------------------------------ Massenmodell [G]/[A]
# Überschussmassen (Plummer, Radius a [A]) — Name, M [M☉], d [Mpc], l, b, a [Mpc]
QUELLEN = [
    ("Milchstraße",        0.8e12,   0.0,   0.0,   0.0, 0.1),
    ("M31",                1.2e12,   0.78, 121.2, -21.6, 0.1),
    ("M81-Gruppe",         1.0e12,   3.6, 142.1,  40.9, 0.3),
    ("Cen-A-Gruppe",       3.0e12,   3.7, 309.5,  19.4, 0.3),
    ("Virgo",              1.2e15,  16.5, 283.8,  74.5, 2.0),
    ("Großer Attraktor",   5.0e15,  65.0, 307.0,   9.0, 10.0),
    ("Perseus-Pisces",     5.0e15,  73.0, 150.6, -13.3, 10.0),
    ("Coma",               1.5e15, 100.0,  58.1,  88.0, 3.0),
    ("Shapley",            1.0e16, 200.0, 312.0,  30.7, 10.0),
]
def orte(quellen):
    return {n: d * einheit(l, b) for n, M, d, l, b, a in quellen}

rng = np.random.default_rng(1)
NMC = 200000
U = rng.random((NMC, 3))
r_pl = 1 / np.sqrt(U[:, 0]**(-2 / 3) - 1)                        # Plummer, a = 1
cth = 2 * U[:, 1] - 1; ph = 2 * math.pi * U[:, 2]; sth = np.sqrt(1 - cth**2)
EINH = np.stack([sth * np.cos(ph), sth * np.sin(ph), cth], 1) * r_pl[:, None]

def M_in(quellen, zentrum, R):
    """Überschussmasse der Quellen in der Kugel (Zentrum [Mpc], R [Mpc]), Plummer ausgewertet mit Monte Carlo."""
    s = 0.0
    for n, M, d, l, b, a in quellen:
        x = d * einheit(l, b) + a * EINH
        s += M * np.mean(np.linalg.norm(x - zentrum, axis=1) < R)
    return s

H("2   Bereiche: mittlerer Kontrast 1 + Δ(<R) [G, Massenmodell unten], θ/3H₀, θ̃/θ_c, w (Top-Hat aus 1)")
P("   Überschussmassen [G], Plummer-Radius [A]:")
for n, M, d, l, b, a in QUELLEN:
    P(f"      {n:18s} {M:8.2e} M☉  d = {d:6.1f} Mpc  (l, b) = ({l:5.1f}°, {b:5.1f}°)  a = {a:4.1f} Mpc")
P("   Nicht enthalten: Leeren (Lokale Leere u. a., negativer Überschuss); 1 + Δ großer Kugeln ist daher eine obere Schätzung.")
O = orte(QUELLEN)
LG_bary = (0.8e12 * O["Milchstraße"] + 1.2e12 * O["M31"]) / 2.0e12
ZENTREN = [("uns (Milchstraße)", np.zeros(3)), ("Virgo", O["Virgo"]), ("Großer Attraktor", O["Großer Attraktor"]), ("Shapley", O["Shapley"])]
RADIEN = (1, 3, 10, 30, 70, 200)
P("   Zentrum               R[Mpc]   M_über[M☉]   1+Δ       θ/3H₀    θ̃/θ_c      w         umgekehrt")
ERG2 = {}
for zn, zp in ZENTREN:
    for R in RADIEN:
        Mx = M_in(QUELLEN, zp, R)
        D = 1 + Mx / (rho_m_Msun_Mpc3 * 4 / 3 * math.pi * R**3)
        th, thm, w, um = lies(D)
        ERG2[(zn, R)] = (D, th, w, um)
        P(f"   {zn:20s} {R:6.0f}   {Mx:10.3e}   {D:8.3f}   {th:7.4f}   {3*H0*thm/thc:8.2f}   {w:.5f}   {'ja (w = 0)' if um else 'nein'}")
def R_um(M):
    return (3 * M / (4 * math.pi * rho_m_Msun_Mpc3 * D_ta))**(1 / 3)
P(f"   Umkehrradius (Punktmasse, 1+Δ_ta = {D_ta:.2f}):  R_ta = (3M/(4πρ̄_m(1+Δ_ta)))^(1/3)")
for n, M in (("Lokale Gruppe 2·10¹²", 2e12), ("Virgo 1,2·10¹⁵", 1.2e15), ("Großer Attraktor 5·10¹⁵", 5e15), ("Shapley 10¹⁶", 1e16)):
    Rt = R_um(M)
    th1, thm1, w1, _ = lies(1 + M / (rho_m_Msun_Mpc3 * 4 / 3 * math.pi * (1.5 * Rt)**3))
    P(f"      {n:26s} R_ta = {Rt:6.2f} Mpc;  bei 1,5·R_ta: w = {w1:.4f};  Schale w = ½ bei R = {R_um(M)*(D_ta/ (D_bei_w(0.5)[0]))**(1/3):6.2f} Mpc, w = 0,9 bei {R_um(M)*(D_ta/(D_bei_w(0.9)[0]))**(1/3):6.2f} Mpc")
P(f"   Gemessen [G]: Nullgeschwindigkeitsfläche Lokale Gruppe ≈ 0,96 Mpc, Virgo ≈ 7,2 Mpc (Karachentsev).  Abstand LG–Virgo 16,5 Mpc.")

# ------------------------------------------------------------------ 3 Geschwindigkeit
H("3   Pekuliargeschwindigkeit der Lokalen Gruppe: linear mitatmend (Schließungsregel = Standardmodell mit Ω_m = (1+Q)Ω_b)")
def wachstum_f():
    """f = d ln D/d ln a heute aus der Wachstumsgleichung im Kern-Hintergrund (Staub + Vakuumteil, Strahlung)."""
    def E2(a): return Om * a**-3 + Or * a**-4 + OL
    def dlnE_dlna(a): return (-3 * Om * a**-3 - 4 * Or * a**-4) / (2 * E2(a))
    def rhs(x, s):
        a = math.exp(x); D, Dp = s
        return [Dp, -(2 + dlnE_dlna(a)) * Dp + 1.5 * Om * a**-3 / E2(a) * D]
    a_i = 1e-3
    sol = solve_ivp(rhs, (math.log(a_i), 0.0), [a_i, a_i], rtol=1e-10, atol=1e-14)
    return sol.y[1, -1] / sol.y[0, -1]
f_ode = wachstum_f(); f_55 = Om**0.55
P(f"   f = d ln D/d ln a heute: Wachstumsgleichung {f_ode:.4f}, Ω_m^0,55 = {f_55:.4f}")
P("   v = (f H₀/4π)∫δ r̂/r² dV = (2f/(3H₀Ω_m))·Σ G·ΔM_i·d_i/(d_i² + a_i²)^{3/2} (Plummer)")
def g_ueber(quellen, ort, ohne=("Milchstraße", "M31")):
    g = np.zeros(3); je = {}
    for n, M, d, l, b, a in quellen:
        if n in ohne:
            continue
        x = (d * einheit(l, b) - ort) * Mpc
        r = np.linalg.norm(x); gi = G * M * Msun * x / (r**2 + (a * Mpc)**2)**1.5
        g += gi; je[n] = gi
    return g, je
VF = 2 * f_ode / (3 * H0 * Om)
def zeige(quellen, titel):
    g, je = g_ueber(quellen, LG_bary)
    v = VF * g / 1e3
    P(f"   {titel}")
    for n, gi in je.items():
        vi = VF * gi / 1e3
        P(f"      {n:18s} |v| = {np.linalg.norm(vi):7.1f} km/s, Anteil längs Σ = {np.dot(vi, v)/np.dot(v, v)*100:6.1f} %")
    l_, b_ = lb(v)
    P(f"      Σ: |v| = {np.linalg.norm(v):6.1f} km/s nach (l, b) = ({l_:5.1f}°, {b_:5.1f}°), Winkel zum CMB-Dipol {winkel(v, CMB):5.1f}°;  Komponente längs Dipol {np.dot(v, CMB):6.1f} km/s (gemessen {V_CMB:.0f} [G])")
    return v, je
v_basis, je_basis = zeige(QUELLEN, "Grundmassen [G]:")
VAR = [list(q) for q in QUELLEN]
for q in VAR:
    if q[0] == "Großer Attraktor": q[1], q[5] = 2e16, 20.0
    if q[0] == "Shapley": q[1], q[5] = 5e16, 20.0
v_var, _ = zeige([tuple(q) for q in VAR], "Variante: Großer Attraktor 2·10¹⁶ (a = 20 Mpc), Shapley 5·10¹⁶ (a = 20 Mpc) [G, Obergrenze]:")
LV = (-0.8 * rho_m_Msun_Mpc3 * 4 / 3 * math.pi * 15**3)
v_lv, _ = zeige(QUELLEN + [("Lokale Leere [A]", LV, 25.0, 30.0, 15.0, 5.0)],
                f"Grundmassen + Lokale Leere als Fehlmasse {LV:.2e} M☉ (δ = −0,8, R = 15 Mpc, 25 Mpc bei (30°, 15°)) [A]:")
# Virgo kugelig nichtlinear
D_V = ERG2[("Virgo", 30)][0]
MV = M_in(QUELLEN, O["Virgo"], 16.5)
D_V16 = 1 + MV / (rho_m_Msun_Mpc3 * 4 / 3 * math.pi * 16.5**3)
thV, _, wV, _ = lies(D_V16)
vV_nl = (1 - thV) * H0 * 16.5 * Mpc / 1e3
vV_lin = f_ode * (D_V16 - 1) / 3 * H0 * 16.5 * Mpc / 1e3
P(f"   Virgo kugelig: 1 + Δ(<16,5 Mpc um Virgo) = {D_V16:.3f}; Einfall nichtlinear (Top-Hat) {vV_nl:.0f} km/s, linear fΔ/3·H₀R {vV_lin:.0f} km/s; w der Schale an der LG = {wV:.4f}")
P(f"   [G] gemessener Einfall der LG auf Virgo ≈ 185–250 km/s.")

# ------------------------------------------------------------------ 4 Anteile, frühere Lesart
H("4   Anteile Virgo / Großer Attraktor / Perseus-Pisces / Shapley: linear 1/r² (Überschussmasse) gegen frühere 1/r-Lesart")
drei = [q for q in QUELLEN if q[0] in ("Virgo", "Großer Attraktor", "Perseus-Pisces", "Shapley")]
def anteile(gew):
    vs = {n: gew(M, d) * einheit(l, b) for n, M, d, l, b, a in drei}
    S = sum(vs.values())
    return {n: np.dot(v, S) / np.dot(S, S) * 100 for n, v in vs.items()}, S
lin, Slin = anteile(lambda M, d: M / d**2)
alt, Salt = anteile(lambda M, d: math.sqrt(M) / d)
P("   Quelle              linear 1/r²    frühere 1/r (√M/r)")
for n in lin:
    P(f"   {n:18s} {lin[n]:9.1f} %    {alt[n]:9.1f} %")
for nm, S in (("linear", Slin), ("1/r", Salt)):
    l_, b_ = lb(S); P(f"   Richtung {nm:6s}: ({l_:5.1f}°, {b_:5.1f}°), Winkel zum CMB-Dipol {winkel(S, CMB):5.1f}°")
P("   Frühere Kette: ν(a_N) auf das Gesamtfeld, Skala a/H₀ = 1 529 km/s (M_b), 3 945 (M_beob), 5 403 (Σ 1/r, M_beob).")
g_Mb = sum(G * fb * M * Msun / (d * Mpc)**2 * einheit(l, b) for n, M, d, l, b, a in drei if n != "Perseus-Pisces")
gN = np.linalg.norm(g_Mb); a_nu = float(K.a_regel(gN))
P(f"   Dieselben Quellen (Virgo, GA, Shapley; M_b = f_b·M): ν-Regel in die lineare Formel statt G·ΔM/r²: v = 2f/(3H₀Ω_m)·ν g_N = {VF*a_nu/1e3:.0f} km/s (ν = {a_nu/gN:.1f})")
g_M = sum(G * M * Msun / (d * Mpc)**2 * einheit(l, b) for n, M, d, l, b, a in drei if n != "Perseus-Pisces")
P(f"   mitatmend mit (1+Q)M_b = M: v = {VF*np.linalg.norm(g_M)/1e3:.0f} km/s.  Verhältnis ν/(1+Q) = {a_nu/gN/(1+Q):.2f}: Mode und Regel zugleich zählt den Zusatz doppelt.")

# ------------------------------------------------------------------ 5 Wechselradius
H("5   Wechsel örtlich (L) → mitatmend (G): Umkehrradius gegen R_Q; Zeitargument mit Mode gegen örtliche Regel")
lnq = -math.log(qk)
P(f"   R_Q = r*/(−ln q_kosm), −ln q_kosm = {lnq:.4f}")
for n, Mb in (("Milchstraße 6·10¹⁰", 6e10), ("Lokale Gruppe 1,8·10¹¹ (Sterne) [G]", 1.8e11), ("Lokale Gruppe 3·10¹¹ (mit Gas) [A]", 3e11),
              ("Virgo f_b·1,2·10¹⁵", fb * 1.2e15)):
    rs = K.r_stern(Mb * Msun)
    P(f"      {n:36s} r* = {rs/kpc:7.1f} kpc, R_Q = {rs/lnq/kpc:8.1f} kpc;  (1+Q)M_b = {(1+Q)*Mb:.2e} M☉, R_ta dazu = {R_um((1+Q)*Mb):.2f} Mpc")
t0_s = t0 / H0
P(f"   Zeitargument: radiale Bahn von R = 0 bei t = 0 bis zur Umkehr R_ta genau heute (t₀ = {t0_s/K.Gyr:.3f} Gyr), Kraft a(R) − Ω_Λ H₀² R")
def R_ta_zeit(a_fn, OL=OL):
    """a_fn(R) nach innen gerichteter Betrag [m/s²]; Φ(R_a) − Φ(R) = ∫_R^{R_a} a dR − ½Ω_ΛH₀²(R_a² − R²). Löse t(R_a) = t₀."""
    def dPhi(R, Ra):
        return quad(a_fn, R, Ra, limit=200)[0] - 0.5 * OL * H0**2 * (Ra**2 - R**2)
    def t_bis(Ra):
        f = lambda u: 2 * Ra * u / math.sqrt(max(2 * dPhi(Ra * (1 - u * u), Ra), 1e-300))   # R = Ra(1 − u²)
        return quad(f, 1e-9, 1.0, limit=200)[0]
    R_eq = brentq(lambda R: a_fn(R) - OL * H0**2 * R, 1e-3 * Mpc, 1e3 * Mpc) if OL > 0 else 30 * Mpc   # Kraft null (Λ gleich Zug)
    return brentq(lambda Ra: t_bis(Ra) - t0_s, 1e-3 * R_eq, 0.999 * R_eq, xtol=1e-5 * Mpc)
def newton(M):  return lambda R: G * M * Msun / R**2
def regel(Mb):  return lambda R: float(K.a_regel(G * Mb * Msun / R**2))
M_LG_R0 = brentq(lambda M: R_ta_zeit(newton(M)) - 0.96 * Mpc, 1e11, 1e14, xtol=1e9)
P(f"   Lokale Gruppe, R₀ = 0,96 Mpc [G] ⇒ Newton-Masse M = {M_LG_R0:.3e} M☉ (Top-Hat-Formel: {4/3*math.pi*rho_m_Msun_Mpc3*D_ta*0.96**3:.3e})")
for n, Mb in (("1,8·10¹¹ [G]", 1.8e11), ("3·10¹¹ [A]", 3e11)):
    Ra_mode = R_ta_zeit(newton((1 + Q) * Mb)); Ra_nu = R_ta_zeit(regel(Mb))
    P(f"      M_b = {n:13s}: mitatmend/Mode (1+Q)M_b → R_ta = {Ra_mode/Mpc:.3f} Mpc;  örtliche Regel ν(a_N) bis R_ta → R_ta = {Ra_nu/Mpc:.3f} Mpc")
Mb_nu = brentq(lambda Mb: R_ta_zeit(regel(Mb)) - 0.96 * Mpc, 1e8, 1e12, xtol=1e6)
Mb_mode = M_LG_R0 / (1 + Q)
P(f"      R₀ = 0,96 Mpc verlangt: mit Mode M_b = {Mb_mode:.2e} M☉;  mit ν(a_N) bis R₀ M_b = {Mb_nu:.2e} M☉  (Sterne MW + M31 ≈ 1,8·10¹¹ [G])")
def regel_efe(Mb, ge):
    """Regel im äußeren Feld, gleichgerichtet s² = s_i² + s_e² (kleinstes ν, günstigster Fall für die Regel)."""
    def a(R):
        gi = G * Mb * Msun / R**2
        s = math.sqrt((gi + ge) / K.A0)
        return gi / (-math.expm1(-s))
    return a
ge_dip = V_CMB * 1e3 / VF
P(f"   Äußeres Feld der LG aus dem Dipol, linear: g_e = v/(2f/(3H₀Ω_m)) = {ge_dip:.3e} m/s² = {ge_dip/K.A0:.4f} a₀ (ν(g_e) = {float(K.nu(ge_dip)):.2f})")
for n, ge in (("Virgo+GA+Shapley M_b-Newton", 9.843e-14), ("Virgo 1,2·10¹⁵ linear", 6.1e-13), ("aus dem Dipol 620 km/s", ge_dip)):
    Ra = R_ta_zeit(regel_efe(1.8e11, ge)); Mb_e = brentq(lambda Mb: R_ta_zeit(regel_efe(Mb, ge)) - 0.96 * Mpc, 1e8, 1e13, xtol=1e6)
    P(f"      mit äußerem Feld g_e = {ge:.2e} m/s² ({n}) [s² = s_i² + s_e²]: M_b = 1,8·10¹¹ → R_ta = {Ra/Mpc:.3f} Mpc; R₀ = 0,96 verlangt M_b = {Mb_e:.2e} M☉")
def regel_efe_senk(Mb, ge):
    def a(R):
        gi = G * Mb * Msun / R**2
        s = ((gi / K.A0)**2 + (ge / K.A0)**2)**0.25
        return gi / (-math.expm1(-s))
    return a
Ra = R_ta_zeit(regel_efe_senk(1.8e11, ge_dip)); Mb_e = brentq(lambda Mb: R_ta_zeit(regel_efe_senk(Mb, ge_dip)) - 0.96 * Mpc, 1e8, 1e13, xtol=1e6)
P(f"      dasselbe senkrecht [s⁴ = s_i⁴ + s_e⁴]: M_b = 1,8·10¹¹ → R_ta = {Ra/Mpc:.3f} Mpc; R₀ = 0,96 verlangt M_b = {Mb_e:.2e} M☉")
ge_Q = (math.log(1 - 1 / (1 + Q)))**2 * K.A0
P(f"      ν(g_e) = 1 + Q verlangte g_e = {ge_Q:.2e} m/s² = {ge_Q/K.A0:.4f} a₀ (die Schwelle von R_Q)")
MV_R0 = brentq(lambda M: R_ta_zeit(newton(M)) - 7.2 * Mpc, 1e13, 1e17, xtol=1e11)
MbV_nu = brentq(lambda Mb: R_ta_zeit(regel(Mb)) - 7.2 * Mpc, 1e9, 1e16, xtol=1e8)
P(f"   Virgo, R₀ = 7,2 Mpc [G] ⇒ Newton-Masse {MV_R0:.3e} M☉; Baryonen dazu mit Mode {MV_R0/(1+Q):.2e}; mit ν(a_N) bis R₀ verlangt M_b = {MbV_nu:.2e} M☉")
aLG = float(K.a_regel(G * 1.8e11 * Msun / (0.96 * Mpc)**2)); aN = G * 1.8e11 * Msun / (0.96 * Mpc)**2
P(f"   An R₀ der LG (M_b = 1,8·10¹¹): a_N/a₀ = {aN/K.A0:.2e}, ν = {aLG/aN:.1f} gegen 1 + Q = {1+Q:.2f};  Wirkmasse νM_b = {aLG/aN*1.8e11:.2e} M☉")
P("   Einfall der Zwerge auf die LG außerhalb R₀: w der Schale (Masse 2·10¹²):")
for R in (1.0, 1.2, 1.5, 2.0, 3.0):
    D = 1 + 2e12 / (rho_m_Msun_Mpc3 * 4 / 3 * math.pi * R**3)
    th, thm, w, um = lies(D)
    P(f"      R = {R:3.1f} Mpc: 1+Δ = {D:6.2f}, θ/3H₀ = {th:.3f}, w = {w:.4f}{'  (umgekehrt)' if um else ''}")

P("   Kette (Regel je Ebene):")
P("      Sonne in der Milchstraße (8,2 kpc < R_Q = 51 kpc): umgekehrt, w = 0, örtliche Regel ν(|g_N|) auf das Gesamtfeld")
P("      Milchstraße in der LG (0,47 Mpc vom Schwerpunkt < R_ta ≈ 1 Mpc): umgekehrt, w = 0, örtliche Regel auf das Gesamtfeld mit g_e")
P("      LG gegen Virgo (16,5 Mpc > R_ta,Virgo ≈ 7–8,5 Mpc), GA, Perseus-Pisces, Shapley: w = 0,995–0,997, linear, Staub mit 1 + Q; A = g_e")

# ------------------------------------------------------------------ 6 von Hand
H("6   Unabhängige Prüfung von Hand (Zahlen eingetippt, keine Kernfunktionen)")
H0h = 67.4e3 / 3.0857e22
P(f"   θ_c/3H₀ = (1,1e-10/2,998e8)/(3·{H0h:.4e}) = {(1.1e-10/2.998e8)/(3*H0h):.5f};  w(θ = 3H₀) = 1/(1 + 0,05600²) → 1 − w = {0.056**2/(1+0.056**2):.3e}")
g_h = 6.674e-11 * 1.2e15 * 1.989e30 / (16.5 * 3.0857e22)**2
v_h = 2 * 0.5275 / (3 * H0h * 0.3134) * g_h / 1e3
P(f"   Virgo als Punkt 1,2·10¹⁵ M☉ bei 16,5 Mpc: g = {g_h:.4e} m/s², v = 2f/(3H₀Ω_m)·g = {v_h:.0f} km/s (f = 0,5275)")
rs_h = math.sqrt(6.674e-11 * 1.8e11 * 1.989e30 / 1.1e-10) / 3.0857e19
P(f"   R_Q(LG, 1,8·10¹¹) = r*/(−ln(5,364/6,364)) = {rs_h:.2f} kpc/{-math.log(5.364/6.364):.4f} = {rs_h/(-math.log(5.364/6.364)):.1f} kpc")
t0h = t0_s
M_lb = math.pi**2 / 8 * (0.96 * 3.0857e22)**3 / (6.674e-11 * t0h**2) / 1.989e30
P(f"   Lynden-Bell ohne Λ: M = π²R₀³/(8Gt₀²) = {M_lb:.3e} M☉ (Skript mit Λ {M_LG_R0:.3e}; Λ erhöht M)")
vf = (6.674e-11 * 1.8e11 * 1.989e30 * 1.1e-10)**0.25
P(f"   Tiefe Regel ohne Λ: Φ = v_f² ln R, t = R_a√π/(v_f√2) ⇒ R_a = t₀v_f√(2/π) = {t0_s*vf*math.sqrt(2/math.pi)/3.0857e22:.3f} Mpc (v_f = {vf/1e3:.1f} km/s); "
  f"Skript mit ν, Ω_Λ = 0: {R_ta_zeit(regel(1.8e11), OL=0.0)/Mpc:.3f} Mpc")
P(f"   Newton ohne Λ, M = 1,734·10¹²: R₀ = (8GMt₀²/π²)^(1/3) = {(8*6.674e-11*1.734e12*1.989e30*t0_s**2/math.pi**2)**(1/3)/3.0857e22:.3f} Mpc; Skript Ω_Λ = 0: {R_ta_zeit(newton(1.734e12), OL=0.0)/Mpc:.3f} Mpc")
P(f"   Einstein-de Sitter 1+Δ_ta = 9π²/16 = 5,5517 (Skript {tophat(e_eds, 1.0, 0.0, 0.0)['D']:.4f}); tiefe Wirkmasse an 0,96 Mpc: r·√(M_b a₀/G) = "
  f"{0.96*3.0857e22*math.sqrt(1.8e11*1.989e30*1.1e-10/6.674e-11)/1.989e30:.2e} M☉")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "trichterkette_gewicht.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
