# -*- coding: utf-8 -*-
"""
Gesamtlauf: das Modell als Ganzes aus dem einen Kern (ftm_kern.py), in der Reihenfolge der tragenden Kette.
Aufruf: python ftm_gesamtlauf.py   → schreibt ftm_gesamtlauf.txt
  a Setzungen und Konstanten    b Rückkehr (eine Gleichung)    c Ebenen/Leiter 3Δs : r* : ℓ₀
  d Galaxien (SPARC)            e Kosmos                       f Naht und Haufen
  g Vorhersagen (03)            h Prüfstand und Gesetze        i Skripte ohne Kern
Alle Formeln aus dem Kern; hier nur Eingaben (Massen, Radien), Datenlesen und Tabellen.
"""
import io, os, sys, math, time, subprocess, re
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
import ftm_kern as K
from ftm_kern import (G, c, M_sun, m_p, M_earth, AU, pc, kpc, Mpc, Gpc, Gyr, yr, A0, A0_LIT, Q_KOSM, q_KOSM,
                      OMEGA_B, OMEGA_M, OMEGA_R, OMEGA_DE, H0_KMS, H0_si)
import sparc_flaechendichte as sf, sparc_universal as su      # nur Datenlesen

t_start = time.time()
OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 100); P(t); P("=" * 100)
def tab(zeilen, kopfzeile):
    b = [max(len(str(r[i])) for r in [kopfzeile] + zeilen) for i in range(len(kopfzeile))]
    f = lambda r: " | ".join(str(x).ljust(b[i]) if i == 0 else str(x).rjust(b[i]) for i, x in enumerate(r))
    P(f(kopfzeile)); P("-+-".join("-" * x for x in b))
    for r in zeilen:
        P(f(r))
g3 = lambda x: f"{x:.3g}"
e3 = lambda x: f"{x:.3e}"

H0 = H0_si()
LAM = K.schliesslaenge(A0)
MW = 6e10 * M_sun
P("FTM-Gesamtlauf – das Modell aus dem Kern ftm_kern.py (a₀ = %.2e m/s², H₀ = %.1f km/s/Mpc)" % (A0, H0_KMS))

# ======================================================================== a Setzungen und Konstanten
kopf("a  Setzungen und Konstanten")
TH = K.theta_(A0, H0)
tab([["a₀ (Standardwert)", e3(A0), "m/s²"],
     ["a₀ Literatur (Vergleich)", e3(A0_LIT), "m/s²"],
     ["cH₀/2π", e3(K.a0_modell()), "m/s²"],
     ["Θ = 2πa₀/(cH₀)", f"{TH:.4f}", ""],
     ["Λ = c²/2a₀", f"{LAM / Gpc:.3f}", "Gpc"],
     ["Λ/(πc/H₀)", f"{LAM / (math.pi * c / H0):.4f}", ""],
     ["ℓ₀ = 2Λ", f"{K.ell0() / Gpc:.3f}", "Gpc"],
     ["τ_D = c/a₀", f"{K.tau_D() / Gyr:.1f}", "Mrd. J."],
     ["θ_c = a₀/c = 1/τ_D", e3(1 / K.tau_D()), "1/s"],
     ["Q = Ω_c/Ω_b (gemessen eingesetzt)", f"{Q_KOSM:.4f}", ""],
     ["q_kosm = Q/(1 + Q)", f"{q_KOSM:.4f}", ""],
     ["2π − 1 (Zählung je Radiant)", f"{2 * math.pi - 1:.4f}", ""],
     ["Ω_b, Ω_m, Ω_r, Ω_DE", f"{OMEGA_B:.4f}, {OMEGA_M:.4f}, {OMEGA_R:.2e}, {OMEGA_DE:.4f}", ""]],
    ["Größe", "Wert", "Einheit"])

# ======================================================================== b Rückkehr
kopf("b  Rückkehr: û = e^{−√((p + Γ + θ)τ_D)} – Regel (p = θ = 0), Gewicht (Γ = 0), Kern (p = iω)")
P("Regel: Γ = a_N/c; q = rueckkehr(0, Γ), ν = 1/(1 − q), χ = s·q/(1 − q), Q_lok = q/(1 − q), d ln a/d ln a_N = 1 − χ/2")
zeilen = []
for x in (100, 25, 10, 1, 0.13, (math.log(q_KOSM))**2, 0.01, 0.001):
    Gam = x * A0 / c
    q = float(K.rueckkehr(0.0, Gam)); s = -math.log(q)
    nu_ = 1 / (1 - q); chi_ = s * q / (1 - q)
    assert abs(nu_ - float(K.nu(x * A0))) < 1e-9 * nu_ and abs(chi_ - float(K.chi(x * A0))) < 1e-9
    zeilen.append([f"{x:.4g}", f"{s:.4f}", f"{q:.4f}", f"{nu_:.4f}", f"{chi_:.4f}", f"{q / (1 - q):.3f}", f"{1 - chi_ / 2:.4f}"])
tab(zeilen, ["a_N/a₀", "s", "q", "ν", "χ", "Q_lok", "Steigung"])
P()
P("Gewicht: 1 − w = rueckkehr(0, 0, θ̃), w = ½ bei θ̃ = θ_c·(ln 2)²")
zeilen = []
H1090 = H0 * K.E_hintergrund(1090)
a0_1090 = A0 * K.a0_frueh_relativ(1090)
for lab, th, a0x in (("θ̃ = θ_c", 1 / K.tau_D(), A0), ("θ̃ = θ_c(ln 2)²", math.log(2)**2 / K.tau_D(), A0),
                     ("heute 3H₀", 3 * H0, A0), ("heute 3H₀, a₀ = cH₀/2π", 3 * H0, K.a0_modell()),
                     ("Rekombination 3H(1090), a₀(1090)", 3 * H1090, a0_1090), ("gebunden (θ̃ = 0)", 0.0, A0)):
    rest = float(K.rueckkehr(0.0, 0.0, th, a0x)); ex = float(K.rueckkehr_exponent(0.0, 0.0, th, a0x))
    zeilen.append([lab, f"{th * K.tau_D(a0x):.4g}", f"{ex:.4g}", f"{rest:.4g}", f"{1 - rest:.4g}"])
tab(zeilen, ["Zustand", "θ̃τ_D", "Exponent", "1 − w", "w"])
P()
P("Kern bei Bahnfrequenz: K̂(iω) = rueckkehr(iω, Γ); ωτ_D ≫ Γτ_D ⇒ die Regel ist stationär")
zeilen = []
for lab, M, r in (("Sonne, 1 AE", M_sun, AU), ("Sonne, r*", M_sun, K.r_stern(M_sun)), ("Milchstraße, r*", MW, K.r_stern(MW)),
                  ("Milchstraße, R_Q", MW, K.r_Q(MW)), ("Milchstraße, 580 kpc", MW, 580 * kpc)):
    aN = G * M / r**2; om = math.sqrt(float(K.a_regel(aN)) / r)
    ex = complex(K.rueckkehr_exponent(1j * om, aN / c))
    zeilen.append([lab, f"{aN * K.tau_D() / c:.3g}", f"{om * K.tau_D():.3g}", f"{-ex.real / math.log(10):.4g}", f"{-ex.imag:.4g}",
                   f"{float(K.q_(aN)):.4f}"])
tab(zeilen, ["Ort", "Γτ_D", "ωτ_D", "log₁₀|K̂|", "arg K̂ [rad]", "q = K̂(0)"])

# ======================================================================== c Ebenen/Leiter
kopf("c  Ebenen/Leiter: 3Δs : r* : ℓ₀ (r*² = 3Δs·ℓ₀, ℓ₀ = c²/a₀ für alle Ebenen gleich)")
Lz = K.zelle_aus_a0(A0)
rho_b = OMEGA_B * 3 * H0**2 / (8 * math.pi * G)
M_zelle = rho_b * Lz**3
zeilen = []
for lab, M in (("Proton", m_p), ("Erde", M_earth), ("Sonne", M_sun), ("Milchstraße 6·10¹⁰ M☉", MW),
               ("Haufen 1,7·10¹⁴ M☉", 1.7e14 * M_sun), ("Zelle (Baryonen ρ_b·L³)", M_zelle)):
    d3 = 3 * K.ds(M); rs = K.r_stern(M); l0 = K.ell0()
    zeilen.append([lab, e3(M / M_sun), e3(d3), e3(rs), f"{rs / d3:.3e}", f"{rs * rs / (d3 * l0):.6f}", f"{math.log10(rs / d3):.2f}"])
tab(zeilen, ["Ebene", "M [M☉]", "3Δs [m]", "r* [m]", "r*/3Δs = ℓ₀/r*", "r*²/(3Δs·ℓ₀)", "log₁₀"])
P(f"Zelle: L = {Lz / Gpc:.2f} Gpc, M_b = {M_zelle / M_sun:.3e} M☉, Λ = ν_kosm·Δs_b = {(1 + Q_KOSM) * K.ds(M_zelle) / Gpc:.3f} Gpc "
  f"(gegen c²/2a₀ = {LAM / Gpc:.3f} Gpc), r*/Λ = {K.r_stern(M_zelle) / LAM:.4f}")

# ======================================================================== d Galaxien (SPARC)
kopf("d  Galaxien: SPARC (175 Kurven), Regel a = ν(a_N)·a_N, Υ_d = 0,5, Υ_b = 0,7 fest")
tv = su.tab_voll()
alle, ausw = [], []
per_gal = []
for nm, d in tv.items():
    try:
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
    except OSError:
        continue
    Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() == 0:
        continue
    gN, go = Vb2[m] / R[m], Vo[m]**2 / R[m]
    alle.append((gN, go))
    if d['Q'] <= 2 and d['inc'] >= 30 and m.sum() >= 5:
        ausw.append((gN, go)); per_gal.append(np.median(np.log10(go / K.a_regel(gN))))
def statistik(paare, a0=A0):
    gN = np.concatenate([p[0] for p in paare]); go = np.concatenate([p[1] for p in paare])
    r = np.log10(go / K.a_regel(gN, a0)); rn = np.log10(go / gN)
    return len(paare), gN.size, np.median(r), r.mean(), r.std(), np.median(rn), rn.std(), gN, go
def a0_best(paare):
    gN = np.concatenate([p[0] for p in paare]); go = np.concatenate([p[1] for p in paare])
    x = np.linspace(-0.4, 0.3, 141); f = [np.sum(np.log10(go / K.a_regel(gN, A0 * 10**y))**2) for y in x]
    return A0 * 10**x[int(np.argmin(f))]
zeilen = []
for lab, pp in (("alle Kurven, alle Punkte", alle), ("Q ≤ 2, Neigung ≥ 30°, ≥ 5 Punkte", ausw)):
    ng, npkt, med, mit, sd, medN, sdN, _, _ = statistik(pp)
    zeilen.append([lab, ng, npkt, f"{med:+.3f}", f"{mit:+.3f}", f"{sd:.3f}", f"{medN:+.3f}", f"{sdN:.3f}", f"{a0_best(pp) * 1e10:.3f}"])
tab(zeilen, ["Auswahl", "Gal.", "Punkte", "Rest Median", "Mittel", "Streuung", "Newton Med.", "Streuung", "a₀ gelesen [1e-10]"])
pg = np.array(per_gal)
P(f"Rest je Galaxie (Auswahl, Median der Punkte): Median {np.median(pg):+.3f} dex, Streuung {pg.std():.3f} dex, n = {pg.size}")
# Rest nach Beschleunigungsbereich
_, _, _, _, _, _, _, gN, go = statistik(ausw)
r = np.log10(go / K.a_regel(gN)); x = np.log10(gN / A0)
zeilen = []
for lo, hi in ((-3, -1.5), (-1.5, -1), (-1, -0.5), (-0.5, 0), (0, 0.5), (0.5, 2)):
    k = (x >= lo) & (x < hi)
    if k.sum():
        zeilen.append([f"[{lo:+.1f}, {hi:+.1f})", int(k.sum()), f"{np.median(r[k]):+.3f}", f"{r[k].std():.3f}"])
tab(zeilen, ["log a_N/a₀", "Punkte", "Rest Median", "Streuung"])
# Tully-Fisher
V, Mb, Vr = [], [], []
for nm, d in tv.items():
    if d['Vf'] > 0 and d['Q'] <= 2 and d['inc'] >= 30:
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
        k = (Vb2 > 0)[-3:]
        if not k.any():
            continue
        V.append(d['Vf'] * 1e3); Mb.append((0.5 * d['L'] + 1.33 * d['MHI']) * 1e9 * M_sun)
        Vr.append(np.mean(np.sqrt(K.a_regel(Vb2[-3:][k] / R[-3:][k]) * R[-3:][k])))   # Regel an den drei äußersten Punkten
V, Mb, Vr = np.array(V), np.array(Mb), np.array(Vr)
lV, lM = np.log10(V / 1e3), np.log10(Mb / M_sun)
Ntf = lM - 4 * lV
steig, achse = np.polyfit(lV, lM, 1)
a0_tf = np.median(V**4 / (G * Mb))
P(f"Tully-Fisher (V_flach, M_b = 0,5·L + 1,33·M_HI, n = {V.size}): Nullpunkt log M_b − 4 log V = {np.median(Ntf):.3f} "
  f"(Modell log(1/(Ga₀)) = {np.log10(1e12 / (G * A0 * M_sun)):.3f}), Streuung {Ntf.std():.3f} dex; freie Steigung {steig:.2f}; "
  f"a₀ aus V⁴/(GM_b) Median {a0_tf * 1e10:.3f}·10⁻¹⁰")
P(f"   dieselbe Größe mit der Regel an den drei äußersten Punkten statt V_flach: Nullpunkt {np.median(lM - 4 * np.log10(Vr / 1e3)):.3f}, "
  f"Streuung {(lM - 4 * np.log10(Vr / 1e3)).std():.3f} dex; log V_flach − log V_Regel: Median {np.median(lV - np.log10(Vr / 1e3)):+.3f}, "
  f"Streuung {(lV - np.log10(Vr / 1e3)).std():.3f} dex")

# ======================================================================== e Kosmos
kopf("e  Kosmos: Gewicht, Vakuumteil, Zelle, Q-Fenster")
Ls = K.lambda_vakuum(LAM); Hs = K.H_vakuum(LAM)
dL = K.d_LSS()
Qmin, Qmax = K.q_fenster()
nu_k = 1 + Q_KOSM
tab([["1 − w heute (θ̃ = 3H₀)", f"{1 - K.gewicht_w(3 * H0):.4e}"],
     ["1 − w Rekombination: Exponent √(θ̃τ_D)", f"{float(K.rueckkehr_exponent(0, 0, 3 * H1090, a0_1090)):.1f}"],
     ["Ω_s = (π²/45)Θ²", f"{K.omega_vakuum():.4f}"],
     ["Ω_s zweiter Weg Λ_s c²/(3H₀²)", f"{Ls * c**2 / (3 * H0**2):.4f}"],
     ["Ω_s/Ω_DE", f"{K.omega_vakuum() / OMEGA_DE:.3f}"],
     ["Ω_s/(1 − w) heute", f"{K.omega_vakuum() / (1 - K.gewicht_w(3 * H0)):.2f}"],
     ["Λ_s = π⁴/(15Λ²) [1/m²]", e3(Ls)],
     ["H_s = √(Λ_s c²/3) [km/s/Mpc]", f"{Hs * Mpc / 1e3:.2f}"],
     ["r*(M_b,Zelle)/Λ = √(6/ν_kosm)", f"{math.sqrt(6 / nu_k):.4f}"],
     ["Zelle L = (3c⁴/(2Gρ_m a₀))^{1/3} [Gpc]", f"{Lz / Gpc:.2f}"],
     ["d_LSS [Gpc]", f"{dL / Gpc:.2f}"],
     ["L/d_LSS", f"{Lz / dL:.4f}"],
     ["a₀ höchstens (L ≥ 0,985·d_LSS) [1e-10]", f"{K.a0_aus_zelle(K.KREISPAAR_GRENZE * dL) * 1e10:.4f}"],
     ["a₀(z)/a₀: z = 1 / 10 / 1090", f"{K.a0_frueh_relativ(1):.5f} / {K.a0_frueh_relativ(10):.4f} / {K.a0_frueh_relativ(1090):.3f}"],
     ["Ω_s(1090) (Θ mit a₀(1090), H(1090))", e3(K.K_S * K.theta_(a0_1090, H1090)**2)],
     ["Q-Fenster (unten Schluss, oben Zelle Würfel)", f"{Qmin:.2f} ≤ Q ≤ {Qmax:.3f}"],
     ["Q gemessen / 2π − 1", f"{Q_KOSM:.3f} / {2 * math.pi - 1:.3f}"]],
    ["Größe", "Wert"])

# ======================================================================== f Naht und Haufen
kopf("f  Naht und Haufen: R_Q = r*/(−ln q_kosm), R_ta (Umkehr heute, Maß (1 + Q)·M_b), M_b,Q")
lq = -math.log(q_KOSM)
D_ta = K.dichte_umkehr_heute()
P(f"a_N/a₀ an R_Q = (ln q)² = {lq**2:.4f}; R_Q/r* = {1 / lq:.3f} (≈ Q + ½ = {Q_KOSM + 0.5:.3f}); χ(R_Q) = {float(K.chi(lq**2 * A0)):.4f}; "
  f"ν(R_Q) − 1 = {float(K.nu(lq**2 * A0)) - 1:.4f} = Q; 1 + Δ_ta = {D_ta:.3f}")
zeilen = []
for Mb_ in (1e10, 6e10, 1e12, 4.5e13):
    M = Mb_ * M_sun; rta = K.umkehrradius_heute(nu_k * M); aN = G * M / rta**2
    zeilen.append([f"{Mb_:.2g}", f"{K.r_stern(M) / kpc:.2f}", f"{K.r_Q(M) / kpc:.1f}", f"{rta / kpc:.0f}", f"{rta / K.r_Q(M):.2f}",
                   f"{aN / A0:.2e}", f"{float(K.nu(aN)):.1f}"])
tab(zeilen, ["M_b [M☉]", "r* [kpc]", "R_Q [kpc]", "R_ta [kpc]", "R_ta/R_Q", "a_N/a₀(R_ta)", "ν(R_ta)"])
MbQ = K.masse_bQ()
from scipy.optimize import brentq
M2 = brentq(lambda lm: math.log(K.r_Q(10**lm * M_sun) / (G * nu_k * 10**lm * M_sun / (100 * H0**2))**(1 / 3)), 10, 17)
P(f"M_b,Q = (1 + Q)²(ln q)⁶a₀³/(10⁴GH₀⁴) = {MbQ / M_sun:.3e} M☉; zweiter Weg R_Q = R₂₀₀ numerisch: {10**M2:.3e} M☉; "
  f"M_b,Q ∝ Θ⁴: bei a₀ = cH₀/2π {K.masse_bQ(K.a0_modell()) / M_sun:.3e} M☉")

# ======================================================================== g Vorhersagen (03)
kopf("g  Vorhersagen (Tabelle in 03) als Zahlen")
aNs = G * MW / (8.2 * kpc)**2
def hq(r, D, M_cl=1e14 * M_sun, a=3 * kpc):
    gi = G * K.hernquist_masse(MW, a, r) / r**2; ge = G * M_cl * min(1.0, D / Mpc) / D**2
    return math.sqrt(float(K.a_regel(math.hypot(gi, ge))) / math.hypot(gi, ge) * gi * r) / math.sqrt(float(K.a_regel(gi)) * r)
a50 = G * MW / (50 * kpc)**2
zeilen = [
    ["Weite Doppelsterne 5 / 30 kAE: a_N/a₀ (ν, wäre die Regel wirksam)",
     f"{G * M_sun / (5e3 * AU)**2 / A0:.2f} ({float(K.nu(G * M_sun / (5e3 * AU)**2)):.2f}) / {G * M_sun / (3e4 * AU)**2 / A0:.3f} ({float(K.nu(G * M_sun / (3e4 * AU)**2)):.2f})",
     "Newton (Einzeltrichter)"],
    ["Pulsare am Sonnenort Ṗ/P radial: Regel / Newton [s⁻¹]", f"{float(K.a_regel(aNs)) / c:.2e} / {aNs / c:.2e}", "6,3 / 4,1·10⁻¹⁹"],
    ["Galaxie im Haufen, 80 kpc: v/v_frei bei 0,5 / 2 Mpc", f"{hq(80 * kpc, 0.5 * Mpc):.2f} / {hq(80 * kpc, 2 * Mpc):.2f}", "0,50 / 0,77"],
    ["Massenänderung: Verzögerung ν·r/c bei 50 kpc [Mio. J.]", f"{float(K.nu(a50)) * 50 * kpc / c / yr / 1e6:.2f}", "≈ 1"],
    ["Endantwort d ln a/d ln M = 1 − χ/2 bei 15 / 30 / 50 kpc",
     " / ".join(f"{1 - float(K.chi(G * MW / (r * kpc)**2)) / 2:.3f}" for r in (15, 30, 50)), "0,54–0,59"],
    ["Schwelle wandert nicht: a₀(z = 1)/a₀ fest gegen Λ ∝ a", f"{K.a0_frueh_relativ(1):.4f} gegen 2", "keine Verdopplung"],
    ["Zelle = Durchmesser der letzten Streufläche: L/d_LSS", f"{Lz / dL:.3f}", "≈ 1"],
    ["Kosmische Streuung des Quadrupols √(2/5)", f"{math.sqrt(2 / 5):.2f}", "0,63"],
    ["Galaxien bei R_Q: a/a_RAR, χ", f"1 / {float(K.chi(lq**2 * A0)):.3f}", "1 / 0,917"],
    ["Tiefgrenze σ⁴ = (4/81)GMa₀: 10⁷ / 2·10⁸ M☉ [km/s]",
     f"{(4 / 81 * K.v_flach(1e7 * M_sun)**4)**0.25 / 1e3:.1f} / {(4 / 81 * K.v_flach(2e8 * M_sun)**4)**0.25 / 1e3:.1f}", "9,2 / 19,5"],
]
tab(zeilen, ["Vorhersage", "Lauf", "Modell"])
P("Nicht im Gesamtlauf gerechnet (eigene Rechnungen): Quadrupol im Sonnensystem, vertikale Flächendichte am Sonnenort, Feldform innerer Punkte,")
P("frisch kollabierte Galaxien, Kollaps nach Umkehr, M_dyn/M_lens, Isobaren, Einseitigkeit (m = 1), gegenläufige Scheiben, Linsen um isolierte")
P("Galaxien (N*), Strukturwachstum, Haufenstöße, Toomre-Integral, Sattelpunkt Sonne–Erde.")

# ======================================================================== h Prüfstand und Gesetze
kopf("h  Prüfstand und Gesetze")
for skript in ("ftm_pruefstand.py", "ftm_gesetze.py"):
    t0 = time.time()
    r_ = subprocess.run([sys.executable, os.path.join(HIER, skript)], cwd=HIER, capture_output=True, text=True, encoding="utf-8")
    zeilen_ = [z for z in r_.stdout.strip().split("\n") if z.strip()]
    summe = zeilen_[-1] if zeilen_ else f"(keine Ausgabe, Fehler: {r_.stderr.strip()[-200:]})"
    P(f"{skript:20s}: {summe}   [{time.time() - t0:.0f} s]")
    for z in zeilen_:
        if re.search(r"\|\s*(ABWEICHUNG|WIDERSPRUCH|SPANNUNG|HINWEIS)\s*(\||$)", z):
            P("   " + re.sub(r"\s{2,}", " ", z.strip())[:160])

# ======================================================================== i Skripte ohne Kern
kopf("i  Hinweis: Skripte in weltmodell/rechnungen, die ftm_kern nicht importieren")
ohne = sorted(f for f in os.listdir(HIER) if f.endswith(".py") and "ftm_kern" not in open(os.path.join(HIER, f), encoding="utf-8", errors="ignore").read())
P(f"{len(ohne)} von {sum(1 for f in os.listdir(HIER) if f.endswith('.py'))}:")
for i in range(0, len(ohne), 5):
    P("  " + ", ".join(ohne[i:i + 5]))

P(); P(f"Laufzeit {time.time() - t_start:.0f} s")
open(os.path.join(HIER, "ftm_gesamtlauf.txt"), "w", encoding="utf-8").write("\n".join(OUT) + "\n")
