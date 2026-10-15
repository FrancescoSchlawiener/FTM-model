# -*- coding: utf-8 -*-
"""
Weiße Zwerge als Gedächtnis der Tiefe: Folgen die Grenzen He | C/O | O/Ne aus der Leiter der Zündtiefen?
Ein-Zonen-Trichter (wie fusion_trichter): kT_c/(μ m_u c²) = f·δ − (m_e/(μ_e m_u))·Π(x),  δ = 3Δs/R, f = θ_n (Polytrope).
Rechnet nur mit ftm_kern (Δs, Konstanten). [G] = Wert aus dem Gedächtnis, nicht nachgeprüft; [A] = Annahme.
Blickrichtung: T lokal gemessen (eigene Uhr am Ort); δ = Taktunterschied gegen den flachen Außenraum. Für die Wärme im Kern zählt
nur der Unterschied innerhalb des Kerns (Δδ_innen); die Hülle verschiebt den Takt des ganzen Kerns gleich und heizt ihn nicht.
Teile: 1 Kurve T_max(M/M_Grenze)  2 Grenzmassen  3 Blitz oder ruhig  4 Welche Tiefe hält die Asche?  5 Gedächtnis  6 Gegenproben
"""
import os, sys, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s=""):
    out.append(s)

m_u = 1.66053906660e-27          # CODATA 2018, fehlt im Kern
keV = 1e3 * e_el
muc2 = m_u * c**2
mec2 = m_e * c**2
lam_e = compton_red(m_e)

# ------------------------------------------------------------------ Polytropen
def lane_emden(n):
    def rhs(x, y):
        return [y[1], -max(y[0], 0)**n - 2 * y[1] / x]
    x0 = 1e-6
    ev = lambda x, y: y[0]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(rhs, [x0, 20], [1 - x0**2 / 6, -x0 / 3], events=ev, rtol=1e-11, atol=1e-13)
    x1 = s.t_events[0][0]; mw = -x1**2 * s.y_events[0][0][1]
    return dict(theta=x1 / ((n + 1) * mw), D=x1**3 / (3 * mw))
LE = {1.5: lane_emden(1.5), 3.0: lane_emden(3.0)}

def Pi_(x):
    """Π = P_e/(n_e m_e c²) des kalten Elektronengases (exakt)."""
    return (x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x)) / (8 * x**3)

def zustand(M, R, mu, mue, n):
    """kT_c (lokal), x = λ̄_e k_F, ρ_c, δ, Wärmetiefe δ_W = kT/(f μ m_u c²)."""
    t = LE[n]
    rho = t["D"] * 3 * M / (4 * math.pi * R**3)
    x = lam_e * (3 * math.pi**2 * rho / (mue * m_u)) ** (1 / 3)
    d = 3 * ds(M) / R
    kT = mu * muc2 * (t["theta"] * d - m_e / (mue * m_u) * Pi_(x))
    return kT, x, rho, d, kT / (t["theta"] * mu * muc2)

def Tmax(M, mu, mue, n):
    r = minimize_scalar(lambda lr: -zustand(M, 10**lr, mu, mue, n)[0], bounds=(3, 11), method="bounded", options=dict(xatol=1e-8))
    return zustand(M, 10**r.x, mu, mue, n), 10**r.x

def M_grenze(n, mue=2.0):
    Kur = (3 * math.pi**2) ** (1 / 3) * hbar * c / 4 / (mue * m_u) ** (4 / 3)
    return (Kur * (3 * LE[n]["D"] / (4 * math.pi)) ** (1 / 3) / (LE[n]["theta"] * G)) ** 1.5

def mu_mix(X):   # X: {(Z, A): Massenanteil}, voll ionisiert
    return 1 / sum(x * (Z + 1) / A for (Z, A), x in X.items())

MG3 = M_grenze(3.0); MG15 = M_grenze(1.5)
MS = M_sun
mu_He = mu_mix({(2, 4): 1.0})
mu_CO = mu_mix({(6, 12): 0.5, (8, 16): 0.5})                       # [A] C:O = 1:1 nach Masse
mu_ONe = mu_mix({(8, 16): 0.6, (10, 20): 0.3, (12, 24): 0.1})      # [A] O/Ne/Mg 0,6/0,3/0,1
stufen = [("He-Kern → 3α", mu_He, 1.0e8), ("C/O-Kern → C", mu_CO, 6.0e8), ("O/Ne-Kern → Ne", mu_ONe, 1.2e9)]

P("WEISSE ZWERGE ALS GEDÄCHTNIS DER TIEFE")
P("kT_c/(μ m_u c²) = f·δ − (m_e/(μ_e m_u))·Π(x),  δ = 3Δs/R,  f = θ_n;  kT lokal, δ gegen den flachen Außenraum")
P(f"μ: He {mu_He:.3f}, C/O {mu_CO:.3f} [A 1:1], O/Ne {mu_ONe:.3f} [A 0,6/0,3/0,1]; μ_e = 2.  M_Grenze: n = 3 → {MG3/MS:.3f} M☉, n = 1,5 → {MG15/MS:.3f} M☉")
P()

# ================================================================== 1 universelle Kurve
P("1) T_max HÄNGT NUR AN m = M/M_Grenze:  kT_max = (μ/μ_e)·m_ec²·g(m)")
P("   (Tiefe ∝ M/R und Zählregel ∝ M^{2/3}/R² laufen mit derselben Strukturzahl; g ist für jede Polytrope dieselbe Funktion.)")
# M_Grenze skaliert mit μ_e⁻²; darum g mit μ_e = 2 rechnen und durch (μ/μ_e) teilen
def g2(m, n=3.0):
    (kT, x, rho, d, dW), R = Tmax(m * M_grenze(n, 2.0), 1.0, 2.0, n)
    return kT / (0.5 * mec2), x, rho, R
P(f"   {'m':>6s} {'g(n=3)':>9s} {'g(n=1,5)':>9s} {'x bei T_max':>11s}")
for m in (0.1, 0.2, 0.3, 0.323, 0.4, 0.5, 0.6, 0.7, 0.72, 0.8, 0.9, 0.95, 0.99, 0.999):
    a, xa, _, _ = g2(m, 3.0); b, _, _, _ = g2(m, 1.5)
    P(f"   {m:6.3f} {a:9.5f} {b:9.5f} {xa:11.3f}")
gmax = -minimize_scalar(lambda x: -(x / 4 - Pi_(x)), bounds=(0.05, 20), method="bounded").fun
P(f"   Decke g(1) = max_x[x/4 − Π(x)] = {gmax:.4f}")
P()

# ================================================================== 2 Grenzmassen
P("2) GRENZMASSEN: Bei welcher Kernmasse erreicht der von der Zählregel gehaltene Kern die Zündtiefe der nächsten Stufe?")
P("   Kern unter der Grenzmasse: T_max(M_K) < T_zünd → die Stufe wird nie erreicht, der Kern bleibt als Zwerg dieser Asche.")
P(f"   {'Stufe':16s} {'T_zünd':>8s} {'δ_zünd':>9s} {'g nötig':>8s} {'m':>7s} {'M_K(n=3)':>9s} {'M_K(n=1,5)':>10s} {'ρ_c[kg/m³]':>11s} {'x':>5s}")
grenz = {}
for nm, mu, T in stufen:
    gneed = kB * T / ((mu / 2) * mec2)
    dz = kB * T / (LE[3.0]["theta"] * mu * muc2)
    if gneed >= gmax:
        P(f"   {nm:16s} {T:8.1e} {dz:9.2e} {gneed:8.4f}   nie (g nötig > Decke {gmax:.4f}): nur ein Kern, der nicht entartet, zündet")
        grenz[nm] = None
        continue
    m = brentq(lambda m: g2(m)[0] - gneed, 0.01, 0.9999)
    _, x, rho, R = g2(m)
    grenz[nm] = m
    P(f"   {nm:16s} {T:8.1e} {dz:9.2e} {gneed:8.4f} {m:7.4f} {m*MG3/MS:9.3f} {m*MG15/MS:10.3f} {rho:11.2e} {x:5.2f}")
P("   [G] Vergleich: He-Blitz bei Kern 0,47 M☉ (m = 0,32); Grenze C/O|O/Ne ≈ 1,05 M☉ (m = 0,72); O/Ne-Kerne über M_Grenze kollabieren (Elektroneneinfang).")
P()
P("   Kohlenstoff: Die Zündtemperatur hängt von der Dichte ab (Abschirmung) [G: 6·10⁸ K bei 10⁸–10⁹ kg/m³, tiefer bei höherer Dichte].")
P(f"   {'T_C':>8s} {'g nötig':>8s} {'m':>7s} {'M_K(n=3)':>9s} {'ρ_c':>10s}")
for T in (3.5e8, 4.0e8, 4.5e8, 4.8e8, 4.9e8):
    gneed = kB * T / ((mu_CO / 2) * mec2)
    m = brentq(lambda m: g2(m)[0] - gneed, 0.01, 0.99999)
    _, x, rho, R = g2(m)
    P(f"   {T:8.1e} {gneed:8.4f} {m:7.4f} {m*MG3/MS:9.3f} {rho:10.2e}")
P()
# Eichung an der He-Grenze: ein Faktor κ auf T (Ein-Zonen-Fehler), dann C-Grenze vorhersagen
m_He_G = 0.47 * MS / MG3
kappa = (kB * 1e8 / ((mu_He / 2) * mec2)) / g2(m_He_G)[0]
P(f"   Eichung am He-Blitz [G 0,47 M☉ → m = {m_He_G:.3f}]: Die Ein-Zonen-Rechnung gibt dort g = {g2(m_He_G)[0]:.4f}; wirklich nötig {kB*1e8/((mu_He/2)*mec2):.4f}")
P(f"   → Faktor κ = T_wirklich/T_Ein-Zone = {kappa:.3f} (die Rechnung unterschätzt T bei gegebener Masse).")
P(f"   Mit demselben κ: Decke C/O = κ·(μ/μ_e)·g(1)·m_ec² = {kappa*(mu_CO/2)*gmax*mec2/kB:.2e} K.")
for T in (5.0e8, 6.0e8, 7.0e8):
    gneed = kB * T / ((mu_CO / 2) * mec2) / kappa
    if gneed < gmax:
        m = brentq(lambda m: g2(m)[0] - gneed, 0.01, 0.99999)
        P(f"   C zündet bei {T:.0e} K → m = {m:.3f}, M_K = {m*MG3/MS:.3f} M☉  [G 1,05]")
for T in (1.2e9,):
    gneed = kB * T / ((mu_ONe / 2) * mec2) / kappa
    P(f"   Ne zündet bei {T:.1e} K: nötig g = {gneed:.4f} gegen Decke {gmax:.4f} → " + ("nie" if gneed >= gmax else "erreichbar"))
m105 = 1.05 * MS / MG3
T105 = kappa * g2(m105)[0] * (mu_CO / 2) * mec2 / kB
P(f"   Umgekehrt: Grenze 1,05 M☉ (m = {m105:.3f}) verlangte mit κ eine C-Zündung bei {T105:.2e} K, ohne κ bei {T105/kappa:.2e} K.")
P("   Eine Eichung an He trifft C nicht: C/O|O/Ne landet bei 1,19–1,45 statt 1,05. Die Grenze folgt nur als 'nahe unter M_Grenze'.")
P()

# ================================================================== 3 Blitz oder ruhig
P("3) BLITZ ODER RUHIG: Entartungsgrad ψ = kT/E_F im Augenblick der Zündung (erster Schnitt T(R) = T_zünd beim Zusammenziehen)")
def zuend_psi(M, mu, T, n=3.0):
    (kTm, xm, _, _, _), Rm = Tmax(M, mu, 2.0, n)
    if kTm < kB * T:
        return None
    lR = brentq(lambda lr: zustand(M, 10**lr, mu, 2.0, n)[0] - kB * T, math.log10(Rm), 12)
    kT, x, rho, d, dW = zustand(M, 10**lR, mu, 2.0, n)
    EF = mec2 * (math.sqrt(1 + x * x) - 1)
    return kT / EF, rho, d
for nm, mu, T in stufen[:2]:
    row = []
    for Mk in (0.35, 0.47, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0, 4.0):
        r = zuend_psi(Mk * MS, mu, T)
        row.append(f"{Mk}:{'—' if r is None else f'{r[0]:.2f}'}")
    P(f"   {nm:16s} ψ bei Zündung (M_K in M☉): " + "  ".join(row))
for nm, mu, T in stufen[:2]:
    try:
        Mb = brentq(lambda m: (zuend_psi(m * MS, mu, T) or (0,))[0] - 1.0, 0.3, 20)
        P(f"   {nm:16s} ψ = 1 (Grenze entartet/ruhig) bei M_K = {Mb:.2f} M☉")
    except ValueError:
        P(f"   {nm:16s} ψ = 1 nicht im Bereich")
P("   [G] He zündet im Blitz für Anfangsmassen ≲ 2 M☉ (Kern 0,45–0,5), ruhig darüber; C zündet entartet (Flamme außermittig) in Kernen 1,05–1,4 M☉.")
P()

# ================================================================== 4 welche Tiefe hält die Asche?
P("4) WELCHE TIEFE HÄLT DIE ASCHE?  Gesamttiefe δ gegen Wärmetiefe δ_W = kT/(f μ m_u c²) = δ − (m_e/(f μ_e m_u))·Π(x)")
def R_kalt(M, mu, n=3.0):
    (kTm, xm, _, _, _), Rm = Tmax(M, mu, 2.0, n)
    return 10**brentq(lambda lr: zustand(M, 10**lr, mu, 2.0, n)[0], 3, math.log10(Rm))
dz = {nm: kB * T / (LE[3.0]["theta"] * mu * muc2) for nm, mu, T in stufen}
P(f"   Zündtiefen (n = 3): 3α {dz['He-Kern → 3α']:.2e}, C {dz['C/O-Kern → C']:.2e}, Ne {dz['O/Ne-Kern → Ne']:.2e}")
R_G = {0.4: 0.0165, 0.6: 0.0125, 1.0: 0.0081, 1.2: 0.0059, 1.35: 0.0040}   # [G ±15 %] kalte Radien (Hamada-Salpeter-Kurve)
for Mk, mu, nm in ((0.4, mu_He, "He-Zwerg"), (0.6, mu_CO, "C/O-Zwerg"), (1.0, mu_CO, "C/O-Zwerg"), (1.2, mu_ONe, "O/Ne-Zwerg"), (1.35, mu_ONe, "O/Ne-Zwerg")):
    Rk = R_kalt(Mk * MS, mu)
    d = 3 * ds(Mk * MS) / Rk
    dG = 3 * ds(Mk * MS) / (R_G[Mk] * R_sun)
    (kTm, xm, _, dm, dWm), Rm = Tmax(Mk * MS, mu, 2.0, 3.0)
    P(f"   {nm:11s} {Mk:4.2f} M☉: δ_kalt Ein-Zone {d:.2e} (R {Rk/1e3:5.0f} km), mit R [G] {R_G[Mk]*R_sun/1e3:5.0f} km: δ = {dG:.2e}"
      f" = {dG/dz['O/Ne-Kern → Ne']:4.1f}·δ_zünd(Ne);  δ_W,max = {dWm:.2e}")
P("   Mit den gemessenen Radien liegt die Gesamttiefe ab 0,6 M☉ über der Ne-Zündtiefe, die des He-Zwergs 0,4 über der C-Zündtiefe (1,4-fach);")
P("   die Ein-Zonen-Radien sind 2,2–2,7-fach zu groß (Ein-Zonen-Fehler), ändern das Urteil ab 1 M☉ nicht: δ allein sortiert die Asche nicht.")
P("   δ_W,max bleibt bei jedem entarteten Kern unter δ_zünd(C) (Decke); O/Ne-Asche entsteht nur, wo der Kern bei der C-Zündung noch nicht voll entartet war.")
P("   Die Asche hält das Maximum der Wärmetiefe δ_W (der Teil der Tiefe, den nicht die Zählregel trägt): Stufe = max_t δ_W(t) ≥ δ_zünd.")
P()

# ================================================================== 5 Gedächtnis
P("5) GEDÄCHTNIS: Endzusammensetzung = f(max_t δ_W), nicht f(Endzustand)")
P("   Preisach-Muster: Die Asche ist der Extremwert der Geschichte der Wärmetiefe; spätere Abnahme (Kühlung, Massenverlust) löscht nichts.")
P("   [G] Beispiele:")
P("   – Asymptotischer Riesenast: Der C/O-Kern wächst, die Hülle wird abgeblasen; Endmasse ≈ Kernmasse am Ende, Asche = größte Wärmetiefe des Kerns.")
P("   – Unterzwerge sdB: Hülle an der Spitze des Riesenasts verloren, Kern 0,47 hatte He gezündet → C/O-Zwerge mit 0,45–0,5 M☉ (unter dem Gipfel 0,6).")
P("   – Doppelsterne: vor dem He-Blitz abgestreift → He-Zwerge 0,15–0,45 M☉ (in der Einzelsternentwicklung zu jung für das Weltalter).")
P("   – Neon-Novae (z. B. V1974 Cyg, QU Vul): Auswurf reich an Ne, der Zwerg verliert Masse, die Asche verrät die O/Ne-Stufe.")
P("   – Kalte Verschmelzung zweier C/O-Zwerge über 1,05 M☉ ohne Zündung bleibt C/O: große δ, kleine δ_W (bestätigt Teil 4).")
P("   Die Grenzen in der Masse der Zwerge sind daher nur für Einzelsterne scharf; die Zuordnung Asche ↔ max δ_W gilt immer.")
P("   Anfangs-End-Massen-Beziehung: nicht aus der Tiefe allein. Sie braucht den Kernzuwachs gegen den Massenverlust der Hülle (Wind).")
P("   Aus der Leiter folgt nur die Ordnung: M_WD ≥ Kernmasse bei der letzten Zündung (He-Blitz 0,47 → Gipfel der Verteilung bei 0,6 [G]);")
P("   O/Ne-Zwerge nur, wo der Kern m ≈ 0,72 erreicht hat, und nie über M_Grenze.")
P("   [G] Gaia (Kilic 2020, Tremblay 2016; arXiv nicht erreichbar): Gipfel 0,59–0,62 M☉, Nebengipfel/Überschuss ≈ 0,8, He-Zwerge wenige Prozent;")
P("   [G] Anfangsmassen: He-Blitz ≲ 2 M☉, O/Ne ab ≈ 7–8 (Super-AGB) bis ≈ 9–10 M☉.")
P()

# ================================================================== 6 Gegenproben
P("6) GEGENPROBEN (anderer Weg)")
# (a) Decke He von Hand
P(f"   (a) Decke He: (μ/μ_e)·g(1)·m_ec² = (4/3)/2·0,0964·511,0 keV = 32,84 keV → /8,617·10⁻⁵ eV/K = 3,81·10⁸ K;"
  f"  Skript {(mu_He/2)*gmax*mec2/kB:.3e} K")
# (b) He-Grenze über die geschlossene nichtrelativistische Form (kT_max = μ m_u (θGM)²/(4B)), n = 1,5 → n = 3 gleich über m
th = LE[3.0]["theta"]; D = LE[3.0]["D"]
Knr = (3 * math.pi**2) ** (2 / 3) * hbar**2 / (5 * m_e)
def Tnr(M, mu=mu_He, mue=2.0):
    B = Knr * (3 * D * M / (4 * math.pi * mue * m_u)) ** (2 / 3) / (mue * m_u)
    return mu * m_u * (th * G * M) ** 2 / (4 * B) / kB
Mnr = brentq(lambda m: Tnr(m * MS) - 1e8, 0.01, 2)
P(f"   (b) He-Grenze nichtrelativistisch geschlossen (T_max ∝ M^{{4/3}}): {Mnr:.3f} M☉ (n = 3); exakter Druck {grenz['He-Kern → 3α']*MG3/MS:.3f} M☉"
  " – nahe beieinander, da x ≈ 0,4 noch fast nichtrelativistisch.")
P(f"       Von Hand: T_max ∝ M^{{4/3}}; aus T_nr(1 M☉) = {Tnr(MS):.3e} K → M = (10⁸/{Tnr(MS):.3e})^{{3/4}} = {(1e8/Tnr(MS))**0.75:.3f} M☉")
P()

txt = "\n".join(out)
print(txt)
with open(os.path.join(HERE, "zwerge_tiefengedaechtnis.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
