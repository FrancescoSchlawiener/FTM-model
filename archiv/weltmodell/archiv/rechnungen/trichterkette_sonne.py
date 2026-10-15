# -*- coding: utf-8 -*-
"""
Trichterkette um die Sonne: Liegt das Sonnensystem im Trichter von „Sagittarius“, und ist dieser unser Großer Attraktor?
Lesarten von „Sagittarius“: (L1) Sgr A*, (L2) Zentrum der Milchstraße (Richtung Schütze), (L3) Sagittarius-Zwerggalaxie,
(L4) „Sagittarius“ als Großer Attraktor der Bewegung der Milchstraße / Lokalen Gruppe.
Abschnitte
  1  Kette der Trichter: Δs, r*, g_N am Ort der Sonne, g_N/a₀, ν, Regel, Richtung (l, b); Gesamtfeld am Ort der Sonne
  2  Sgr A*: kompakter Trichter? r*/r₀, Einflussradius GM/σ², Glättungslänge ℓ³ = GM/T; Anteil an der Beschleunigung der Sonne
  3  Bahn der Sonne (v²/R) gegen Regel; Bewegung der Milchstraße als Ganzes (Massenmittel A); Lokale Gruppe in Virgo, GA, Shapley
  4  Anteile im tiefen Bereich (1/r) gegen Newton (1/r²); Richtungen gegen CMB-Dipol und Sagittarius-Richtung
  5  Netzwerk: Sattel zwischen Milchstraße (bzw. Lokaler Gruppe) und Virgo / Großem Attraktor
  6  Unabhängige Prüfung von Hand (eigene Zahlen, ohne Kernfunktionen)
Regel: a = ν(|g_N|)·g_N auf das Gesamtfeld (örtliche Regel), ν = 1/(1 − e^{−s}), s = √(g_N/a₀). Tiefer Bereich je Quelle: √(GMa₀)/r.
[A] Annahme, [G] Wert aus dem Gedächtnis, ungeprüft.  a₀ = A0 aus dem Kern.
"""
import io, os, sys, math
import numpy as np
from scipy.special import i0e, i1e, k0e, k1e
from ftm_kern import G, c, M_sun, GM_sun, pc, kpc, Mpc, A0, H0_si, OMEGA_B, OMEGA_M, ds, r_stern, r0, nu, v_flach

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []
def P(s=""): out.append(s)
def H(t): P(); P("=" * 104); P(t); P("=" * 104)

def einheit(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)])
def lb(v):
    v = np.asarray(v, float); r = np.linalg.norm(v)
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2] / r))
def winkel(u, v):
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))))
def nuf(g): return float(nu(g, A0))
FB = OMEGA_B / OMEGA_M                                   # kosmischer Baryonanteil (Kern), für Haufen-Faltungsmasse [A]

# ------------------------------------------------------------------ Eingaben
R0 = 8.2 * kpc                                           # Abstand Sonne–Zentrum [G]
V_SONNE = 230e3                                          # [G]
# Name, M_beob [M☉] (Beobachter-Masse, meist dynamisch) [G], M_b (Faltung, baryonisch) [G]/[A], Abstand von der Sonne [m], l, b [G]
QUELLEN = [
    ("Sgr A*",            4.3e6,  4.3e6,        8.2 * kpc,    0.0,     0.0),
    ("Milchstraße (Pkt)", 6e10,   6e10,         8.2 * kpc,    0.0,     0.0),
    ("Sgr-Zwerg 18 kpc",  1e9,    1e8,          18 * kpc,     5.6,   -14.2),
    ("Sgr-Zwerg 26,5 kpc",1e9,    1e8,          26.5 * kpc,   5.6,   -14.2),
    ("Große Mag. Wolke",  1.5e11, 3e9,          50 * kpc,   280.5,   -32.9),
    ("M31",               1.5e12, 1.2e11,       780 * kpc,  121.2,   -21.6),
    ("Virgo",             1e15,   FB * 1e15,    16.5 * Mpc, 283.8,    74.5),
    ("Großer Attraktor",  5e15,   FB * 5e15,    65 * Mpc,   307.0,     9.0),
    ("Shapley",           1e16,   FB * 1e16,    200 * Mpc,  312.0,    30.7),
]
CMB_LG = einheit(276.0, 30.0)                            # Bewegung der Lokalen Gruppe gegen den Hintergrund, ~620 km/s [G]
SGR_RICHT = [("Sgr A* / Zentrum (0°, 0°)", einheit(0.0, 0.0)), ("Sgr-Zwerg (5,6°, −14,2°)", einheit(5.6, -14.2))]
NORMA = ("Norma-Haufen A3627 (325,3°, −7,3°) [G]", einheit(325.3, -7.3))

P("Trichterkette um die Sonne – Lesarten von „Sagittarius“")
P(f"a₀ = {A0:.2e} m/s²;  Baryonanteil der Haufen f_b = Ω_b/Ω_m = {FB:.3f} [A];  R₀ = 8,2 kpc, v☉ = 230 km/s [G]")
P("Massen: M_beob = Wert der Punkteliste (Beobachter-Masse; bei LMC, Sgr-Zwerg, M31/LG, Haufen dynamisch gelesen, enthält also die")
P("Zusatzbeschleunigung schon) [G];  M_b = Faltungsmasse: Sgr-Zwerg 10⁸ (Sterne) [G], LMC 3·10⁹ [G], M31 1,2·10¹¹ [G], Haufen f_b·M_beob [A].")
P("Die Regel mit M_beob zählt den Zusatz doppelt; sie steht nur als obere Grenze.")

# ============================================================================ 1
H("1   Kette der Trichter: Werte am Ort der Sonne, je Quelle für sich (die Regel wirkt auf das Gesamtfeld, Zeile Σ)")
P(f"   {'Quelle':19s} {'M_b [M☉]':>9s} {'Δs [m]':>9s} {'r*':>11s} {'d/r*':>8s} {'g_N [m/s²]':>10s} {'g_N/a₀':>9s} {'ν':>7s} "
  f"{'ν·g_N':>9s} {'√(GMa₀)/d':>10s} {'l, b':>13s} | mit M_beob: {'g_N/a₀':>9s} {'√(GMa₀)/d':>10s}")
def rtxt(r):
    return f"{r/pc:8.1f} pc" if r < kpc else (f"{r/kpc:7.2f} kpc" if r < Mpc else f"{r/Mpc:7.2f} Mpc")
VEK = {}
for nm, Mo, Mb, d, l, b in QUELLEN:
    gN = G * Mb * M_sun / d**2; gNo = G * Mo * M_sun / d**2
    VEK[nm] = (gN * einheit(l, b), gNo * einheit(l, b))
    P(f"   {nm:19s} {Mb:9.2e} {ds(Mb*M_sun):9.2e} {rtxt(r_stern(Mb*M_sun)):>11s} {d/r_stern(Mb*M_sun):8.2f} {gN:10.3e} {gN/A0:9.2e} "
      f"{nuf(gN):7.3f} {nuf(gN)*gN:9.2e} {math.sqrt(G*Mb*M_sun*A0)/d:10.3e} {l:6.1f},{b:6.1f} | {'':11s}{gNo/A0:9.2e} {math.sqrt(G*Mo*M_sun*A0)/d:10.3e}")
# Milchstraße als Exponentialscheibe (R_d = 3 kpc [A]), Feld in der Ebene
def expd_g(M, Rd, R):
    y = R / (2 * Rd)
    return G * M * R / (2 * Rd**3) * (i0e(y) * k0e(y) - i1e(y) * k1e(y))
gdisk = float(expd_g(6e10 * M_sun, 3 * kpc, R0))
P(f"   Milchstraße als Exponentialscheibe 6·10¹⁰ M☉, R_d = 3 kpc [A]: g_N(R₀) = {gdisk:.3e} m/s² = {gdisk/A0:.3f} a₀, ν = {nuf(gdisk):.3f}, "
  f"a = {nuf(gdisk)*gdisk:.3e}")
# Gesamtfeld am Ort der Sonne (Milchstraße als Punkt; Sgr-Zwerg mit 26,5 kpc)
TEIL = ["Sgr A*", "Milchstraße (Pkt)", "Sgr-Zwerg 26,5 kpc", "Große Mag. Wolke", "M31", "Virgo", "Großer Attraktor", "Shapley"]
gtot = sum(VEK[n][0] for n in TEIL); gt = np.linalg.norm(gtot)
P(f"   Σ Gesamtfeld (M_b): |g_N| = {gt:.4e} m/s² = {gt/A0:.4f} a₀, ν = {nuf(gt):.4f}, a = {nuf(gt)*gt:.4e} m/s², Richtung "
  f"(l, b) = ({lb(gtot)[0]:.3f}°, {lb(gtot)[1]:.3f}°)")
P("   Anteile am Newton-Gesamtfeld (Komponente längs Σ) – dieselben Anteile trägt a = ν(|g_N|)g_N, weil ν ein gemeinsamer Faktor ist:")
for n in TEIL:
    P(f"      {n:19s} {np.dot(VEK[n][0], gtot)/gt**2*100:10.4f} %   (Betrag/|Σ| {np.linalg.norm(VEK[n][0])/gt:.2e})")
P("   Gezeitenanteil (Feld an der Sonne minus Feld am Zentrum der Milchstraße) – nur er ändert die Bahn der Sonne um das Zentrum:")
GC = R0 * einheit(0, 0)
for nm, Mo, Mb, d, l, b in QUELLEN[2:]:
    x = d * einheit(l, b); xg = x - GC
    gS = G * Mb * M_sun * x / d**3; gG = G * Mb * M_sun * xg / np.linalg.norm(xg)**3
    P(f"      {nm:19s} |g(☉) − g(GC)| = {np.linalg.norm(gS-gG):.3e} m/s² = {np.linalg.norm(gS-gG)/gdisk*100:.2e} % des Feldes der Scheibe")

# ============================================================================ 2
H("2   Sgr A*: kompakter Trichter, Reichweite des Einzeltrichters, Anteil an der Beschleunigung der Sonne")
M_bh = 4.3e6 * M_sun; GMb = G * M_bh
rs_bh = r_stern(M_bh); r0_bh = r0(M_bh)
P(f"   Δs = {ds(M_bh):.4e} m,  r₀ = 6Δs = {r0_bh:.4e} m,  r* = {rs_bh:.4e} m = {rs_bh/pc:.2f} pc,  r*/r₀ = {rs_bh/r0_bh:.3e} (kompakt: r* ≫ Hülle)")
P(f"   Sonne: d/r* = {R0/rs_bh:.1f}  → Sonne liegt weit außerhalb r* von Sgr A*")
P("   Einflussradius r_h = GM/σ² und Glättungslänge ℓ³ = GM/T mit der Elterngezeit T = √(3/2)|E^TF| des Sternfelds um Sgr A*.")
P("   Sternfeld als isotherme Kugel (σ [G]): Φ = 2σ²ln r, E^TF = diag(−8/3, 4/3, 4/3)σ²/r², T(r) = 4σ²/r²; ℓ aus GM/ℓ³ = T(ℓ) → ℓ = GM/(4σ²).")
for sig in (70e3, 100e3, 120e3):
    rh = GMb / sig**2
    # ℓ numerisch aus Fixpunkt mit numerischer Hesse-Form (anderer Weg als die Formel)
    def T_num(r, h=1e-4):
        Phi = lambda x: 2 * sig**2 * math.log(np.linalg.norm(x))
        x0 = np.array([r, 0.0, 0.0]); Hm = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                ei = np.eye(3)[i] * h * r; ej = np.eye(3)[j] * h * r
                Hm[i, j] = (Phi(x0+ei+ej) - Phi(x0+ei-ej) - Phi(x0-ei+ej) + Phi(x0-ei-ej)) / (4 * h * h * r * r)
        E = Hm - np.trace(Hm) / 3 * np.eye(3)
        return math.sqrt(1.5) * np.linalg.norm(E)
    l_ = rh / 4
    for _ in range(60): l_ = (GMb / T_num(l_))**(1 / 3)
    P(f"      σ = {sig/1e3:5.0f} km/s: r_h = {rh/pc:5.2f} pc,  ℓ = GM/4σ² = {rh/4/pc:5.3f} pc (Fixpunkt numerisch {l_/pc:5.3f} pc),  "
      f"ℓ/r* = {l_/rs_bh:.4f},  g_N(ℓ)/a₀ = {GMb/l_**2/A0:.2e}")
P("   Damit endet der Einzeltrichter von Sgr A* bei ~0,5 pc, tief im Newton-Bereich (g ≫ a₀); schon dort geht er ins kollektive Feld")
P("   des Zentrums über. Sein r* = 74 pc hat für ihn keine eigene Bedeutung mehr (ℓ < r*: die Rückkehr setzt erst im Gesamtfeld ein).")
g_bh = GMb / R0**2
P(f"   Am Ort der Sonne: g_N(Sgr A*) = {g_bh:.3e} m/s² = {g_bh/A0:.2e} a₀;  Anteil am Gesamtfeld {g_bh/gt*100:.4f} % (Punkt), "
  f"{g_bh/gdisk*100:.4f} % (Scheibe)")
P(f"   Las man Sgr A* fälschlich als freien kollektiven Trichter (√(GMa₀)/d = {math.sqrt(GMb*A0)/R0:.3e}): {math.sqrt(GMb*A0)/R0/(nuf(gt)*gt)*100:.2f} % — "
  "gilt nicht, denn die Regel sieht nur das Gesamtfeld.")

# ============================================================================ 3
H("3   Welcher Trichter bestimmt die Bahn der Sonne, welcher die Bewegung der Milchstraße / Lokalen Gruppe")
a_obs = V_SONNE**2 / R0
MMW = 6e10 * M_sun
g_pt = G * MMW / R0**2
P(f"   Bahn der Sonne: v²/R₀ = {a_obs:.4e} m/s² (v = 230 km/s [G]) = {a_obs/A0:.3f} a₀")
for lab, g in (("Milchstraße Punkt", g_pt), ("Milchstraße Exp.-Scheibe", gdisk), ("Scheibe + Sgr A*", gdisk + g_bh)):
    a = nuf(g) * g
    P(f"      {lab:26s} g_N = {g:.4e}, ν = {nuf(g):.4f}, a = {a:.4e} → v = {math.sqrt(a*R0)/1e3:6.1f} km/s, a/(v²/R₀) = {a/a_obs:.3f}")
P("   → Die Bahn der Sonne bestimmt der kollektive Trichter der Milchstraße (Zentrum Richtung Schütze, l = 0°).")

# Massenmittel A einer Exponentialscheibe im homogenen äußeren Feld (Lage β zur Ebene); Regel auf |g_i + g_e|
xL, wL = np.polynomial.laguerre.laggauss(160)            # ∫ e^{−u} u f(u) du, u = R/R_d
PH = np.arange(720) * 2 * math.pi / 720
def A_scheibe(M, Rd, ge_vec):
    """A = ⟨ν(|g_i + g_e|)(g_i + g_e)⟩_M, dünne Exponentialscheibe in der Ebene z = 0, g_i in der Ebene zum Zentrum"""
    R = xL * Rd; gi = expd_g(M, Rd, R)
    ex, ey = np.cos(PH), np.sin(PH)
    gx = -gi[:, None] * ex[None, :] + ge_vec[0]; gy = -gi[:, None] * ey[None, :] + ge_vec[1]; gz = np.zeros_like(gx) + ge_vec[2]
    gb = np.sqrt(gx**2 + gy**2 + gz**2); n_ = np.asarray(nu(gb, A0))
    w = (wL * xL)[:, None] / len(PH)                        # Massengewicht u·e^{−u} du dφ/2π, Summe 1
    return np.array([np.sum(w * n_ * gx), np.sum(w * n_ * gy), np.sum(w * n_ * gz)])
chk = A_scheibe(MMW, 3 * kpc, np.array([0.03 * A0, 0, 0]))
P(f"   Massenmittel A der Scheibe (6·10¹⁰ M☉, R_d = 3 kpc): Probe g_e = 0,03 a₀ in der Ebene: A/g_e = {chk[0]/(0.03*A0):.3f} "
  f"(Satz des Modells: 1,27), ν(g_e) = {nuf(0.03*A0):.3f}")

def ext_at_GC(namen, massen="b"):
    k = 1 if massen == "b" else 2
    g = np.zeros(3)
    for nm, Mo, Mb, d, l, b in QUELLEN:
        if nm in namen:
            M = (Mb if massen == "b" else Mo) * M_sun
            xg = d * einheit(l, b) - GC; g += G * M * xg / np.linalg.norm(xg)**3
    return g
P("   Milchstraße als Ganzes im äußeren Feld (am Zentrum; Scheibe liegt in der galaktischen Ebene):")
P(f"   {'äußere Quellen':34s} {'Massen':6s} {'|g_e|/a₀':>9s} {'(l, b) von g_e':>16s} {'A/|g_e|':>8s} {'ν(g_e)':>7s} {'Winkel A,g_e':>12s}")
for namen, lab in ((["Große Mag. Wolke"], "LMC"), (["M31"], "M31"), (["Virgo", "Großer Attraktor", "Shapley"], "Virgo + GA + Shapley"),
                   (["Große Mag. Wolke", "M31", "Virgo", "Großer Attraktor", "Shapley"], "alle")):
    for ms in ("b", "o"):
        ge = ext_at_GC(namen, ms); A = A_scheibe(MMW, 3 * kpc, ge)
        P(f"   {lab:34s} {'M_b' if ms=='b' else 'M_beob':6s} {np.linalg.norm(ge)/A0:9.2e} ({lb(ge)[0]:6.1f}°,{lb(ge)[1]:5.1f}°) "
          f"{np.linalg.norm(A)/np.linalg.norm(ge):8.3f} {nuf(np.linalg.norm(ge)):7.2f} {winkel(A, ge):10.2f}°")
P("   A liegt nach dem Modell zwischen g_e (dichter Kern) und ν(g_e)·g_e; die Milchstraße sitzt fast ganz im dichten Kern (A ≈ g_e).")

# Lokale Gruppe in Virgo, GA, Shapley
H("3b  Beschleunigung der Lokalen Gruppe durch Virgo, Großen Attraktor, Shapley (Ort der Sonne ≈ Ort der LG auf dieser Skala)")
FERN = ["Virgo", "Großer Attraktor", "Shapley"]
def fern_tab(ms):
    rows = []
    for nm, Mo, Mb, d, l, b in QUELLEN:
        if nm in FERN:
            M = (Mb if ms == "b" else Mo) * M_sun; u = einheit(l, b)
            rows.append((nm, M, d, u, G * M / d**2, math.sqrt(G * M * A0) / d, nuf(G * M / d**2) * G * M / d**2))
    return rows
t0 = 1 / H0_si()
for ms, lab in (("b", "M_b = f_b·M_beob (Faltungsmasse) [A]"), ("o", "M_beob [G] (obere Grenze, Zusatz doppelt)")):
    rows = fern_tab(ms)
    P(f"   {lab}")
    P(f"   {'Quelle':18s} {'M [M☉]':>9s} {'d':>9s} {'g_N':>10s} {'g_N/a₀':>9s} {'√(GMa₀)/d':>10s} {'ν(g_N)g_N':>10s} {'Regel/Newton':>12s}")
    for nm, M, d, u, gN, gd, gr in rows:
        P(f"   {nm:18s} {M/M_sun:9.2e} {d/Mpc:6.1f} Mpc {gN:10.3e} {gN/A0:9.2e} {gd:10.3e} {gr:10.3e} {gr/gN:12.2f}")
    gN_v = sum(r[4] * r[3] for r in rows); gd_v = sum(r[5] * r[3] for r in rows)
    gnt = np.linalg.norm(gN_v); ar = nuf(gnt) * gnt
    P(f"      Newton Σ:              |g| = {gnt:.3e} m/s², Richtung ({lb(gN_v)[0]:6.1f}°, {lb(gN_v)[1]:5.1f}°), Winkel zum CMB-Dipol {winkel(gN_v, CMB_LG):5.1f}°, "
      f"Skala |g|/H₀ = {gnt*t0/1e3:6.0f} km/s")
    P(f"      Regel auf Gesamtfeld:  |a| = {ar:.3e} m/s² (ν = {nuf(gnt):.2f}), Richtung = Newton-Richtung, Skala a/H₀ = {ar*t0/1e3:6.0f} km/s")
    P(f"      Σ je Quelle 1/r:       |a| = {np.linalg.norm(gd_v):.3e} m/s², Richtung ({lb(gd_v)[0]:6.1f}°, {lb(gd_v)[1]:5.1f}°), Winkel zum CMB-Dipol "
      f"{winkel(gd_v, CMB_LG):5.1f}°, Skala a/H₀ = {np.linalg.norm(gd_v)*t0/1e3:6.0f} km/s")
    Am = A_scheibe(MMW, 3 * kpc, ext_at_GC(FERN, ms))
    P(f"      Massenmittel A (Milchstraße als Scheibe): |A| = {np.linalg.norm(Am):.3e} m/s² = {np.linalg.norm(Am)/np.linalg.norm(ext_at_GC(FERN, ms)):.3f}·|g_e|, "
      f"Richtung ({lb(Am)[0]:6.1f}°, {lb(Am)[1]:5.1f}°), Winkel zum CMB-Dipol {winkel(Am, CMB_LG):5.1f}°, Skala A/H₀ = {np.linalg.norm(Am)*t0/1e3:6.0f} km/s")
    for nm2, u2 in SGR_RICHT + [NORMA]:
        P(f"         Winkel {nm2:36s}: zu Newton Σ {winkel(gN_v, u2):5.1f}°, zu 1/r-Σ {winkel(gd_v, u2):5.1f}°")
P("   Skala a/H₀ [A]: nur Größenordnung (Beschleunigung × Hubble-Zeit), keine Wachstumsrechnung. CMB-Dipol der LG: 620 km/s, (276°, 30°) [G].")
P("   Für die LG als Ganzes gilt wie oben: A zwischen g_e und ν(g_e)g_e; die Galaxien der LG liegen in ihren dichten Kernen (A ≈ g_e).")

# ============================================================================ 4
H("4   Ordnung der Anteile: Regel im tiefen Bereich (∝ √M/r) gegen Newton (∝ M/r²); Richtungen")
for ms, lab in (("b", "M_b"), ("o", "M_beob")):
    rows = fern_tab(ms); sN = sum(r[4] for r in rows); sD = sum(r[5] for r in rows)
    P(f"   {lab}:  " + " | ".join(f"{r[0]}: Newton {r[4]/sN*100:5.1f} %, 1/r {r[5]/sD*100:5.1f} %" for r in rows))
P("   Winkel zwischen den Richtungen:")
U = {n: einheit(l, b) for n, Mo, Mb, d, l, b in QUELLEN}
for n in FERN:
    P(f"      {n:18s}: zum CMB-Dipol {winkel(U[n], CMB_LG):5.1f}°, zu l = 0° (Sgr A*) {winkel(U[n], SGR_RICHT[0][1]):5.1f}°, "
      f"zum Sgr-Zwerg {winkel(U[n], SGR_RICHT[1][1]):5.1f}°")
P(f"      Norma A3627       : zum CMB-Dipol {winkel(NORMA[1], CMB_LG):5.1f}°, zu l = 0° {winkel(NORMA[1], SGR_RICHT[0][1]):5.1f}°")
P(f"      CMB-Dipol         : zu l = 0° {winkel(CMB_LG, SGR_RICHT[0][1]):5.1f}°, zum Sgr-Zwerg {winkel(CMB_LG, SGR_RICHT[1][1]):5.1f}°")

# ============================================================================ 5
H("5   Netzwerk: Sattel zwischen zwei kollektiven Trichtern (Felder gleich, längs der Verbindung)")
P("   Newton (1/r²):  x/D = √M₁/(√M₁ + √M₂);   tiefe Form je Quelle (1/r, Φ_X = v_flach²·ln r):  x/D = M₁^{1/4}/(M₁^{1/4} + M₂^{1/4}).")
P("   Örtliche Regel auf das Gesamtfeld: a = ν(|g_N|)g_N verschwindet genau dort, wo g_N verschwindet → Sattel am Newton-Punkt.")
PAARE = [("Milchstraße 6·10¹⁰", 6e10, "Virgo", 16.5 * Mpc), ("Lokale Gruppe 1,8·10¹¹ (MW + M31) [A]", 1.8e11, "Virgo", 16.5 * Mpc),
         ("Milchstraße 6·10¹⁰", 6e10, "Großer Attraktor", 65 * Mpc), ("Lokale Gruppe 1,8·10¹¹ (MW + M31) [A]", 1.8e11, "Großer Attraktor", 65 * Mpc)]
MB = {n: Mb for n, Mo, Mb, d, l, b in QUELLEN}; MO = {n: Mo for n, Mo, Mb, d, l, b in QUELLEN}
for n1, M1, n2, D in PAARE:
    for lab, M2 in (("M_b", MB[n2]), ("M_beob", MO[n2])):
        xN = D * math.sqrt(M1) / (math.sqrt(M1) + math.sqrt(M2)); xT = D * M1**0.25 / (M1**0.25 + M2**0.25)
        gxN = G * M1 * M_sun / xN**2; gxT = math.sqrt(G * M1 * M_sun * A0) / xT
        P(f"   {n1:38s} – {n2:16s} ({lab:6s}): Newton-Sattel {xN/kpc:7.0f} kpc (g_N = {gxN/A0:.1e} a₀), "
          f"tiefe Form {xT/kpc:7.0f} kpc (a = {gxT/A0:.1e} a₀);  r*₁ = {r_stern(M1*M_sun)/kpc:5.1f} kpc")
P("   v_flach (Φ_X = v_flach²·ln ν): " + ", ".join(f"{n} {v_flach(MB[n]*M_sun)/1e3:.0f} km/s" for n in ["Milchstraße (Pkt)", "M31", "Virgo", "Großer Attraktor", "Shapley"]))

# ============================================================================ 6
H("6   Unabhängige Prüfung von Hand (Konstanten eingetippt, keine Kernfunktionen)")
GMs, a0h, pch = 1.32712e20, 1.1e-10, 3.0857e16
rsh = math.sqrt(GMs * 4.3e6 / a0h)
P(f"   r*(Sgr A*) = √(1,32712e20·4,3e6/1,1e-10) = {rsh:.4e} m = {rsh/pch:.2f} pc   (Skript: {rs_bh/pc:.2f} pc)")
Dv = 16.5e6 * pch; GMv = GMs * 1e15
gNv = GMv / Dv**2; sv = math.sqrt(gNv / a0h); nuv = 1 / (1 - math.exp(-sv))
P(f"   Virgo (M_beob 10¹⁵) am Ort der LG: g_N = {gNv:.4e}, s = {sv:.4f}, ν = {nuv:.3f}, a = {nuv*gNv:.4e};  tief √(GMa₀)/d = {math.sqrt(GMv*a0h)/Dv:.4e}")
P(f"      Reihe ν ≈ 1/s + 1/2 + s/12 = {1/sv + 0.5 + sv/12:.3f} (Skript ν = {nuf(G*1e15*M_sun/(16.5*Mpc)**2):.3f})")
gmw = GMs * 6e10 / (8.2e3 * pch)**2; smw = math.sqrt(gmw / a0h)
P(f"   Milchstraße (Punkt) an der Sonne: g_N = {gmw:.4e} = {gmw/a0h:.3f} a₀, ν = {1/(1-math.exp(-smw)):.4f}, v = {math.sqrt(gmw/(1-math.exp(-smw))*8.2e3*pch)/1e3:.1f} km/s")
P(f"   Winkel CMB-Dipol – l = 0°: cos = cos 30°·cos 276° = {math.cos(math.radians(30))*math.cos(math.radians(276)):.4f} → "
  f"{math.degrees(math.acos(math.cos(math.radians(30))*math.cos(math.radians(276)))):.1f}°")
P(f"   ℓ(Sgr A*, σ = 100 km/s) = GM/(4σ²) = {GMs*4.3e6/(4*1e10)/pch:.3f} pc")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "trichterkette_sonne.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
