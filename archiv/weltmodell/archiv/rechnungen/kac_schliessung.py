# -*- coding: utf-8 -*-
"""
Kac-Schließung der inneren Richtung: folgen (a) τ_D = 2Λ/c, (b) θ_c = 1/τ_D und die Form von w,
(c) 1 − q_kosm = Ω_b/Ω_m (je Radiant 1/2π), (d) Λ = ν_kosm·ΣΔs_b aus einem Prinzip?

Prinzip: Die innere Richtung ist je Falte ein Kreis 2Λ ohne Umlauf (spiegelsymmetrisch: Halbgerade [0, Λ],
bei Λ reflektierend); das Ganze hat kein Außen; die Kac-Zählung (mittlere Rückkehrzeit, Maß des offenen
Bereichs) gilt für die Faltung wie für die Zelle.

Abschnitte
  0 Grundgrößen
  1 (a) Kac-Zeit: unabhängig von der Dynamik; Massenunabhängigkeit; was zusätzlich nötig ist
  2 (b) Gewicht w als Laplace-Transformierte der Rückkehr bei der Rate θ̃; Formen; Folgen (1 − w, Breiten, Umkehrkante)
  3 (c) Maß des offenen Bereichs: Kac-Budget, Zählung je Radiant, Alternativen
  4 (d) Λ aus Bilanz und Kac-Lemma (Turm); Zellform und a₀(L, Form)
  5 Γ = 0 gegen ν_kosm: Kac-Maß ohne Gefälle (drei Dynamiken)
  6 Atmen als Verlustrate: Zelle als eine Falte, Kandidaten Γ = 0 / 3H / H, Irrfahrt, w(z), Q̃(s_bg), Fixpunkt
  7 Ergebnis
  8 Unabhängige Prüfungen (von Hand, zweiter Weg)
Blickrichtung: θ, H [G]; Γ, a_N [L]; Λ, D, ϖ_R, τ_D [I]. Q_KOSM ist eingesetzte Zahl (Kern), keine Vorhersage.
Rechnet nur mit ftm_kern. Ausgabe: kac_schliessung.txt
"""
import io, os, sys, math
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []
def P(s=""):
    print(s); OUT.append(s)
def Hd(t):
    P(); P("=" * 110); P(t); P("=" * 110)

c, G = K.c, K.G
A0 = K.A0
H0 = K.H0_si()
Lam = K.schliesslaenge(A0)
tauD = K.tau_D(A0)
thc = A0 / c
Q, qk = K.Q_KOSM, K.q_KOSM
nuk = 1 + Q
rng = np.random.default_rng(20260929)

# ------------------------------------------------------------------ 0
Hd("0  Grundgrößen")
P(f"   a₀ = {A0:.3e} m/s²,  Λ = c²/(2a₀) = {Lam/K.Gpc:.4f} Gpc,  τ_D = c/a₀ = {tauD/K.Gyr:.2f} Gyr,  θ_c = a₀/c = {thc:.4e} s⁻¹")
P(f"   H₀ = {H0:.4e} s⁻¹,  3H₀τ_D = {3*H0*tauD:.3f},  H₀τ_D = {H0*tauD:.3f}")
P(f"   eingesetzt (Kern): Q_KOSM = {Q:.4f}, q_KOSM = {qk:.4f}, 1 − q = {1-qk:.5f}, ν_kosm = {nuk:.4f}")

# ------------------------------------------------------------------ 1 (a)
Hd("1  (a) τ_D = 2Λ/c: Kac-Zeit des Kreises 2Λ")
P("   Kac: ⟨T_A⟩ = μ(X)/μ(A) für jede ergodische, maßerhaltende Dynamik. Kreis 2Λ in Zellen der Länge ℓ, Schrittzeit ℓ/c:")
P("   ⟨T_Ort⟩ = (2Λ/ℓ)·(ℓ/c) = 2Λ/c, unabhängig von ℓ (also von der Masse) und von der Dynamik. Prüfung (Kreis mit M Zellen, Ort = 1 Zelle):")
def kac_mittel(M, dyn, n_start=4000, tmax=10**6):
    """mittlere Rückkehrzeit in Zelle 0 (Start in 0) auf dem Kreis mit M Zellen."""
    if dyn == "Bernoulli":
        T = rng.geometric(1 / M, size=n_start)
        return T.mean(), T.std() / math.sqrt(n_start)
    if dyn == "Drehung":
        al = (math.sqrt(5) - 1) / 2
        x0 = rng.random(n_start) / M              # gleichverteilt in Zelle [0, 1/M)
        T = np.zeros(n_start); x = x0.copy(); aktiv = np.ones(n_start, bool); t = 0
        while aktiv.any() and t < tmax:
            t += 1; x = (x + al) % 1.0
            hit = aktiv & (x < 1 / M); T[hit] = t; aktiv &= ~hit
        return T.mean(), T.std() / math.sqrt(n_start)
    if dyn == "Irrfahrt":
        # exakt: Stationärmaß gleichverteilt (symmetrischer Lauf auf dem Kreis) ⇒ ⟨T⟩ = M; Monte-Carlo zur Probe
        pos = np.zeros(n_start, int); T = np.zeros(n_start); aktiv = np.ones(n_start, bool); t = 0
        while aktiv.any() and t < tmax:
            t += 1; pos = (pos + rng.choice((-1, 1), size=n_start)) % M
            hit = aktiv & (pos == 0); T[hit] = t; aktiv &= ~hit
        return T[~aktiv].mean(), T[~aktiv].std() / math.sqrt(n_start)
for M in (20, 60):
    r = [kac_mittel(M, d) for d in ("Bernoulli", "Drehung", "Irrfahrt")]
    P(f"   M = {M:3d}: Bernoulli {r[0][0]:6.2f} ± {r[0][1]:.2f}   goldene Drehung {r[1][0]:6.2f} ± {r[1][1]:.2f}   Irrfahrt {r[2][0]:6.2f} ± {r[2][1]:.2f}   (Kac M = {M})")
P("   Massenunabhängigkeit: Schritt r₀ = 6Δs, D = r₀c/2 = GM/c; die Laplace-Form verlangt τ_D = ϖ_R²/D. τ_D = 2Λ/c ⇔ ϖ_R² = r₀Λ ⇔ a₀ = c²/2Λ:")
for nm, M in (("Sonne", K.M_sun), ("Milchstraße", 6e10 * K.M_sun), ("Haufen", 1.7e14 * K.M_sun)):
    D = G * M / c; rs = K.r_stern(M); r0 = K.r0(M)
    P(f"   {nm:12s} ϖ_R²/D = {rs**2/D/K.Gyr:7.3f} Gyr,  2Λ/c = {2*Lam/c/K.Gyr:7.3f} Gyr,  r*²/(r₀Λ) = {rs**2/(r0*Lam):.9f}")
P("   Massenfreie Zeiten der geschlossenen Halbgeraden mit c-Schritten: Λ/c (Überquerung), 2Λ/c (Kac, jede Dynamik), 4Λ²/(π²D) (unterste Mode, ∝ 1/M).")
P("   Nur 2Λ/c ist massen- und dynamikfrei ⇒ τ_D ∝ Λ/c folgt aus der universellen Schwelle; der Faktor 2 folgt aus Kac.")
P("   Was Kac nicht liefert: dass die Tiefe der Laplace-Form die Tiefe ist, die die Diffusion in einer Kac-Zeit erreicht,")
P("   ϖ_R² = D·⟨T⟩ (Zusatzannahme X_a). Die Rückkehr an den Ort selbst (Start einen Schritt tief) gibt s² = Γ·2r₀/c, nicht Γτ_D:")
M = 6e10 * K.M_sun; r0 = K.r0(M); aN = A0
P(f"   Milchstraße bei a_N = a₀: Start ϖ = r₀: s = √(Γ·2r₀/c) = {math.sqrt(aN/c*2*r0/c):.3e};  Start ϖ = r*: s = √(Γτ_D) = {math.sqrt(aN/c*tauD):.4f}")
P("   Ohne unabhängige Festlegung von Λ ist (a) die Definition Λ := c²/2a₀; Gehalt bekommt sie erst mit (d).")

# ------------------------------------------------------------------ 2 (b)
Hd("2  (b) Gewicht w: Anteil der Rückkehr, die innerhalb der Ausdehnungszeit 1/θ̃ NICHT schließt")
P("   1 − w = E[e^{−θ̃T}] (Rückkehr vor dem Verlust durch Ausdehnung mit Rate θ̃), T erste Rückkehrzeit aus der Tiefe ϖ_R")
P("   auf [0, Λ], reflektierend bei Λ:  1 − w = cosh(k(Λ − ϖ_R))/cosh(kΛ), k = √(θ̃/D), kϖ_R = √x, x = θ̃τ_D.")
P("   Für Λ ≫ ϖ_R (jede Falte unterhalb der Zelle): 1 − w = e^{−√x}, also w = 1 − e^{−√x} = 1/ν(a_N → cθ̃): dieselbe Funktion wie die Regel,")
P("   mit θ̃ statt Γ = a_N/c. Die Rate setzt τ_D: θ_c = 1/τ_D ohne neue Zahl.")
levy = lambda t, t0: math.sqrt(t0 / math.pi) * t**-1.5 * math.exp(-t0 / t)      # Dichte erster Rückkehrzeiten, t0 = τ_D/4
def lap_levy(x):
    t0 = 0.25
    f = lambda t: levy(t, t0) * math.exp(-x * t)
    return quad(f, 0, 1, limit=200)[0] + quad(f, 1, np.inf, limit=200)[0]
P("   Prüfung Laplace-Integral der Dichte (τ_D = 1) gegen e^{−√x}:  " + ",  ".join(f"x={x}: {lap_levy(x):.8f}/{math.exp(-math.sqrt(x)):.8f}" for x in (0.01, 1.0, 17.86)))
forms = {
    "Kac/Lévy, reelle Rate  1 − e^{−√x}": lambda x: -math.expm1(-math.sqrt(x)),
    "Maxwell (Modell)  x²/(1+x²)": lambda x: x * x / (1 + x * x),
    "Exponentialkern, reelle Rate  x/(1+x)": lambda x: x / (1 + x),
    "Lévy bei Frequenz iθ̃  1 − e^{−y}cos y, y = √(x/2)": lambda x: 1 - math.exp(-math.sqrt(x / 2)) * math.cos(math.sqrt(x / 2)),
    "Kreis 2Λ, eine Falte = Zelle  tanh-/cosh-Form (Abschn. 6)": None,
}
P("   Form                                                   w(0,1)   w(1/3)    w(1)    w(3)    w(17,86)   x bei w = ½")
for nm, f in forms.items():
    if f is None:
        continue
    xh = brentq(lambda x: f(x) - 0.5, 1e-6, 50)
    P(f"   {nm:52s} {f(0.1):7.4f}  {f(1/3):7.4f}  {f(1):7.4f}  {f(3):7.4f}  {f(3*H0*tauD):8.5f}   {xh:.4f}")
P("   Maxwell x²/(1+x²) ist Re[1 − Ĝ(iω)] eines Kerns mit EINER Zeit bei FREQUENZ ω = θ̃. Die Kac-Rückkehr hat keine eine Zeit")
P("   (Dichte ∝ t^{−3/2}), und θ̃ ist eine Rate, keine Frequenz: Die Maxwell-Form folgt aus dem Prinzip nicht.")
P("   Die Lévy-Form bei Frequenz überschießt: max w = " + f"{max(1 - math.exp(-y)*math.cos(y) for y in np.linspace(0, 10, 20001)):.4f} (> 1, kein Gewicht).")

wK = lambda x: -math.expm1(-math.sqrt(max(x, 0.0)))
wM = lambda x: x * x / (1 + x * x)
x0 = 3 * H0 * tauD
x0m = 3 * H0 * K.tau_D(K.a0_modell())
z_r = 1090.0
a0r = A0 * K.a0_frueh_relativ(z_r)
xr = 3 * H0 * K.E_hintergrund(z_r) * c / a0r
P()
P("   Folgen im Hintergrund (θ̃ = 3H):")
P(f"   heute, a₀ = 1,1:     x = {x0:.3f}   1 − w: Maxwell {1-wM(x0):.4e}   Kac e^(−√x) {math.exp(-math.sqrt(x0)):.4e}   Verhältnis {math.exp(-math.sqrt(x0))/(1-wM(x0)):.2f}")
P(f"   heute, a₀ = cH₀/2π:  x = {x0m:.3f} (= 6π)  1 − w: Maxwell {1-wM(x0m):.4e}   Kac {math.exp(-math.sqrt(x0m)):.4e}")
P(f"   z = 1090 (a₀(z) = {K.a0_frueh_relativ(z_r):.4f}·a₀): x = {xr:.4e}   1 − w: Maxwell {1-wM(xr):.2e}   Kac e^(−{math.sqrt(xr):.1f}) = 10^{-math.sqrt(xr)/math.log(10):.0f}")
P(f"   w heute im Hintergrund: Maxwell {wM(x0):.5f}, Kac {wK(x0):.5f};  „außerhalb w ≥ 0,99“ hält mit Kac nicht (w = 0,99 bei x = {math.log(100)**2:.2f} > {x0:.2f})")
P(f"   Beziehung Ω_s/(1 − w): Maxwell (1 − w ≈ (Θ/6π)²) → 36π⁴/45 = {36*math.pi**4/45:.2f} in jeder Epoche;  Kac: 1 − w = e^(−√(6π/Θ)), Ω_s/(1 − w) heute = {K.omega_vakuum()/math.exp(-math.sqrt(x0)):.2f}, nicht fest")

# Top-Hat (wie gewicht_preisach: Kern-Hintergrund, Zeit in 1/H₀)
OM, ORr, OL = K.OMEGA_M, K.OMEGA_R, K.OMEGA_DE
E = lambda a: K.E_hintergrund(1 / a - 1)
def t_of_a(a):
    return quad(lambda y: 1 / (y * E(y)), 1e-9, a, limit=200)[0]
ai = 1e-3
def rhs(t, y):
    a, x, v = y
    return [a * E(a), v, -OM / (2 * x * x) - ORr * a**-4 * x + OL * x]
def tophat(di, tmax):
    x = ai * (1 - di / 3); v = x * E(ai) * (1 - di / (3 * (1 + di)))
    def ev_ta(t, y): return y[2]
    ev_ta.direction = -1
    s1 = solve_ivp(rhs, (t_of_a(ai), tmax), [ai, x, v], events=ev_ta, rtol=1e-10, atol=1e-13, max_step=2e-3)
    if len(s1.t_events[0]) == 0:
        return None
    t_ta, x_ta = s1.t_events[0][0], s1.y_events[0][0][1]
    def ev_col(t, y): return y[1] - 0.01 * x_ta
    ev_col.terminal = True; ev_col.direction = -1
    s2 = solve_ivp(rhs, (t_of_a(ai), tmax), [ai, x, v], events=ev_col, rtol=1e-10, atol=1e-14, max_step=2e-3, dense_output=True)
    return dict(t_ta=t_ta, x_ta=x_ta, t_c=s2.t_events[0][0] if len(s2.t_events[0]) else None, sol=s2)
def collapse_at(zc):
    tc = t_of_a(1 / (1 + zc)); lo, hi = 1e-4, 0.05
    for _ in range(40):
        mid = math.sqrt(lo * hi); r = tophat(mid, 1.6 * tc + 0.05)
        if r is None or r["t_c"] is None or r["t_c"] > tc:
            lo = mid
        else:
            hi = mid
    return tophat(hi, 1.6 * tc + 0.05)
thcH = thc / H0
def t_bei_x(r, xw):
    """Zeit vor der Umkehr, zu der θ̃ = x·θ_c (θ = 3ẋ/x fällt monoton bis zur Umkehr)."""
    th = lambda t: 3 * r["sol"].sol(t)[2] / r["sol"].sol(t)[1] - xw * thcH
    return brentq(th, r["t_ta"] * 0.3, r["t_ta"], xtol=1e-12)
tH = 1 / H0 / K.Gyr
x9M, x1M = 3.0, 1 / 3
x9K, x1K = math.log(10)**2, math.log(1 / 0.9)**2
P()
P(f"   Übergang w 0,9 → 0,1 vor der Umkehr (Top-Hat im Kern-Hintergrund). Maxwell: x = {x9M:.3f} → {x1M:.4f};  Kac: x = {x9K:.4f} → {x1K:.5f}")
P("   z_c    Maxwell [Mio. J.]  (gewicht_preisach: 7,24 / 22,48 / 72,00 / 156,40)   Kac [Mio. J.]   Verhältnis   Kac: w = 0,9 vor der Umkehr [Mio. J.]")
breiten = {}
for zc in (3.4, 2.0, 1.0, 0.5):
    r = collapse_at(zc)
    dM = (t_bei_x(r, x1M) - t_bei_x(r, x9M)) * tH * 1e3
    dK = (t_bei_x(r, x1K) - t_bei_x(r, x9K)) * tH * 1e3
    vor = (r["t_ta"] - t_bei_x(r, x9K)) * tH * 1e3
    breiten[zc] = (dM, dK)
    P(f"   {zc:4.1f}   {dM:8.2f}                                                            {dK:8.2f}        {dK/dM:.3f}       {vor:8.2f}")

# Umkehrkante: Schale mit Überschuss eps (wie trichterkette_gewicht), heute 1 + Δ gegen θ̃
def schale(eps, a_i=1e-2):
    E_i = math.sqrt(OM * a_i**-3 + ORr * a_i**-4 + OL)
    def f(t, s):
        a, ad, y, yd = s
        return [ad, -0.5 * OM / a**2 - ORr / a**3 + OL * a, yd, -0.5 * OM * (1 + eps) / y**2 - ORr / a**4 * y + OL * y]
    def heute(t, s): return s[0] - 1.0
    heute.terminal = True
    def kollaps(t, s): return s[2] - 1e-4 * a_i
    kollaps.terminal = True
    sol = solve_ivp(f, (0, 5), [a_i, a_i * E_i, a_i, a_i * E_i * (1 - eps / 3)], events=[heute, kollaps], rtol=1e-10, atol=1e-13)
    if sol.t_events[1].size or not sol.t_events[0].size:
        return np.inf, -1.0
    y, yd = sol.y_events[0][0][2], sol.y_events[0][0][3]
    return (1 + eps) / y**3, yd / y          # 1 + Δ, θ/3H₀
e_ta = brentq(lambda e: schale(e)[1], 1e-3, 0.03, xtol=1e-11)
D_ta = schale(e_ta)[0]
def D_bei(xw):
    thz = xw * thc / (3 * H0)
    return schale(brentq(lambda e: schale(e)[1] - thz, 0.0, e_ta, xtol=1e-12))[0]
P()
P(f"   Umkehrkante heute (Schale, Kern-Hintergrund): Umkehr 1 + Δ_ta = {D_ta:.3f}")
P("   w      Maxwell: x, 1 + Δ       Kac: x, 1 + Δ")
for w in (0.99, 0.9, 0.5, 0.1):
    xm = math.sqrt(w / (1 - w)); xk = math.log(1 / (1 - w))**2
    dm = D_bei(xm) if xm < x0 else float("nan")
    dk = D_bei(xk) if xk < x0 else float("nan")
    P(f"   {w:4.2f}   {xm:7.4f}  {dm:7.3f}        {xk:8.5f}  " + (f"{dk:7.3f}" if xk < x0 else "Hintergrund liegt darunter (w_bg = %.4f)" % wK(x0)))

# ------------------------------------------------------------------ 3 (c)
Hd("3  (c) Maß des offenen Bereichs")
P("   Kac-Budget: Σ_i μ(A_i)/μ(X) = 1 für jede Zerlegung; es legt die Größe eines Teils nicht fest (jede Zerlegung erfüllt es).")
P("   Beispiel: Kreis in 2π Radiant, offener Teil A beliebig: 1/⟨T_A⟩ + 1/⟨T_Rest⟩ = μ(A)/2π + (2π − μ(A))/2π = 1 für jedes μ(A).")
for muA in (0.5, 1.0, 2.0):
    P(f"     μ(A) = {muA:.1f} rad: ⟨T_A⟩ = {2*math.pi/muA:.4f}, ⟨T_Rest⟩ = {2*math.pi/(2*math.pi-muA):.4f}, Summe der Kehrwerte = {muA/(2*math.pi)+(2*math.pi-muA)/(2*math.pi):.12f}")
P("   Die Phasenidentifikation θ ~ θ + 2π gibt das Gesamtmaß 2π, nicht das Maß von A. Mit A = ein Radiant (Zählung je Radiant, Zusatzannahme X_c):")
P(f"     Ω_b/Ω_m = 1/2π = {1/(2*math.pi):.6f},  ν = 2π = {2*math.pi:.5f},  Q = 2π − 1 = {2*math.pi-1:.5f},  q = 1 − 1/2π = {1-1/(2*math.pi):.6f}")
P(f"     (eingesetzt im Kern: Q_KOSM = {Q:.4f}, 1 − q = {1-qk:.5f}; Abstand der Modellzahl {100*((2*math.pi-1)/Q-1):+.2f} % in Q)")
P("   Andere ebenso „natürliche“ Wahlen desselben Kreises (keine wird von Kac ausgezeichnet):")
for nm, m in (("ein Radiant auf dem Kreis 2Λ (Halbgerade: Λ/2π)", 1 / (2 * math.pi)), ("ein Radiant auf der Halbgeraden Λ (Bogen Λ/π)", 1 / math.pi),
              ("Hälfte (Spiegel)", 0.5), ("Ort = Schritt r₀ = 6Δs der Zelle (bei Λ = ν·ΣΔs_b)", 6 / nuk)):
    P(f"     {nm:52s} 1 − q = {m:.5f}  Q = {1/m-1:7.4f}")
P("   Spiegelbild: Halbgerade [0, Λ] = Kreis 2Λ; Ortsstück [0, L] ↔ Bogen 2L ⇒ 1 − q = L/Λ (Kac = gambler's ruin).")
P("   Je Radiant heißt L = Λ/2π: die doppelte Faltung ist der Radius des inneren Kreises (2L = 2Λ/2π).")

# ------------------------------------------------------------------ 4 (d)
Hd("4  (d) Λ = ν_kosm·ΣΔs_b: Kac-Lemma (Turm) und kein Außen")
P("   Kac-Lemma: μ(A)·⟨T_A⟩ = μ(⋃_n T^n A) ≤ μ(X); Gleichheit genau dann, wenn der Turm über A den ganzen Raum füllt (kein Teil ohne Rückkehr).")
P("   Mit μ(A) = ΣΔs_b (Faltung der Baryonen der Zelle), μ(X) = Λ, ⟨T_A⟩ = ν_kosm: ν_kosm·ΣΔs_b ≤ Λ; „kein Außen“ macht daraus Λ = ν_kosm·ΣΔs_b.")
def turm(M1, M2, m, zus):
    """Irrfahrt auf Kreis M1 (enthält A = m Zellen); zus: zweiter Kreis M2 (nicht erreichbar) oder verbunden."""
    if not zus:
        Mges = M1 + M2    # zweiter Kreis ohne Verbindung: Turm füllt nur M1
        return m * (M1 / m), Mges
    Mges = M1 + M2        # verbunden zu einem Kreis M1 + M2
    return m * (Mges / m), Mges
for zus in (False, True):
    t, X = turm(40, 60, 5, zus)
    P(f"   Kreis 40 + Kreis 60 Zellen, A = 5 Zellen, {'verbunden' if zus else 'getrennt (Außen)'}: μ(A)·⟨T_A⟩ = {t:.0f}, μ(X) = {X}")
Mc = 100; mA = 5; nst = 20000
pos = rng.integers(0, mA, nst); T = np.zeros(nst); akt = np.ones(nst, bool); t = 0
while akt.any():
    t += 1; pos = (pos + rng.choice((-1, 1), size=nst)) % Mc
    hit = akt & (pos < mA); T[hit] = t; akt &= ~hit
P(f"   Probe Irrfahrt verbunden (Kreis 100, A = 5 Zellen, 2·10⁴ Läufe): μ(A)·⟨T_A⟩ = {mA*T.mean():.2f} (μ(X) = 100)")
P("   Gehalt: Das Lemma macht aus „ν = Rückkehrzeit des Bereichs ΣΔs_b“ genau Λ = ν·ΣΔs_b; den Wert von ν legt es nicht fest.")
P("   Zusatzannahme X_d: Die Rückkehrzahl der Mode ξ (Ω_m/Ω_b) ist die Kac-Rückkehrzeit des Bereichs ΣΔs_b in der inneren Richtung der Zelle.")
P("   Die Bilanz σ = 0 setzt die gezählte Faltung gleich der Dehnung der Raumschleifen; zu Λ trägt sie keine eigene Gleichung bei.")
P()
P("   Zellform: kürzeste geschlossene Schleife L (Kreispaar-Grenze L ≥ 0,985·d_LSS); Volumen V = f·L³")
dL = K.d_LSS()
formen = (("würfelförmig (T³, kubisch)", 1.0), ("sechseckig (Prisma, Höhe = L)", math.sqrt(3) / 2),
          ("raumzentriert (kürzester Vektor L)", 4 / (3 * math.sqrt(3))), ("flächenzentriert (dichteste, V = L³/√2)", 1 / math.sqrt(2)))
rho_b = K.rho_m_mitbewegt() / nuk
P(f"   d_LSS = {dL/K.Gpc:.3f} Gpc, L_min = {K.KREISPAAR_GRENZE*dL/K.Gpc:.3f} Gpc;  a₀ = 3c⁴/(2G·ν·ρ_b·V)")
P("   Form                                         f        a₀(L_min) ν = 1+Q   a₀(L_min) ν = 2π   L(a₀ = 1,1) ν = 1+Q   L(1,1)/d_LSS")
for nm, f in formen:
    a0f = lambda L, nu_: 3 * c**4 / (2 * G * nu_ * rho_b * f * L**3)
    Lmin = K.KREISPAAR_GRENZE * dL
    L11 = (3 * c**4 / (2 * G * nuk * rho_b * f * A0))**(1 / 3)
    P(f"   {nm:44s} {f:.4f}   {a0f(Lmin, nuk)*1e10:6.3f}·10⁻¹⁰        {a0f(Lmin, 2*math.pi)*1e10:6.3f}·10⁻¹⁰        {L11/K.Gpc:7.3f} Gpc          {L11/dL:.4f}")
P(f"   Probe Kern (a0_aus_zelle, würfelförmig) bei L_min: {K.a0_aus_zelle(K.KREISPAAR_GRENZE*dL)*1e10:.4f}·10⁻¹⁰")
P("   Flächenzentriert ist die dichteste Gitterpackung (Hermite γ₃ = 2^{1/3}): kleinstes V bei gegebenem L ⇒ größtes a₀ jeder flachen Zelle T³.")
P(f"   a₀(z) aus festem ΣΔs_b (mitbewegt, Baryonen erhalten), früh mit Strahlung (Kern): z = 10: {K.a0_frueh_relativ(10):.4f}, z = 1090: {K.a0_frueh_relativ(1090):.4f}")

# ------------------------------------------------------------------ 5
Hd("5  Γ = 0 im homogenen Hintergrund gegen ν_kosm als Rückkehrzahl")
P("   Ohne Verlust (Γ = 0) kehrt auf dem geschlossenen Raum jede Bahn zurück; das Mittel der Rückkehrzeit ist das Kac-Maß, für jede Dynamik.")
P("   Kreis 200 Zellen, A = 32 Zellen (μ(A) = 0,16), Start gleichverteilt in A, 2·10⁴ Läufe:")
Mc, mA, nst = 200, 32, 20000
def rueck(dyn):
    if dyn == "Bernoulli":
        return rng.geometric(mA / Mc, size=nst).astype(float)
    x = rng.random(nst) * mA / Mc if dyn == "Drehung" else None
    pos = rng.integers(0, mA, nst) if dyn == "Irrfahrt" else None
    T = np.zeros(nst); akt = np.ones(nst, bool); t = 0; al = (math.sqrt(5) - 1) / 2
    while akt.any():
        t += 1
        if dyn == "Drehung":
            x = (x + al) % 1.0; inA = x < mA / Mc
        else:
            pos = (pos + rng.choice((-1, 1), size=nst)) % Mc; inA = pos < mA
        hit = akt & inA; T[hit] = t; akt &= ~hit
    return T
for dyn in ("Bernoulli", "Drehung", "Irrfahrt"):
    T = rueck(dyn)
    qg = 1 - mA / Mc
    P(f"   {dyn:9s}: ⟨T⟩ = {T.mean():6.3f} ± {T.std()/math.sqrt(nst):.3f} (Kac 1/μ = {Mc/mA:.3f});  P(T=1) = {np.mean(T==1):.3f} (geometrisch 1 − q = {1-qg:.3f});  P(T>5) = {np.mean(T>5):.3f} (geometrisch q⁵ = {qg**5:.3f})")
P("   Das Mittel ν = 1/μ(A) = 1/(1 − q) ist Topologie und Maß; die geometrische Verteilung Σqⁿ (Bernoulli) verlangt zusätzlich Mischung.")
P("   Kein Widerspruch: Γ = 0 heißt kein Verlust; die Rückkehrgleichung D·G'' = Γ·G gibt dann q = 1 (jede Bahn kehrt zurück),")
P("   das Kac-Maß zählt, wie oft der offene Bereich besucht wird. Das örtliche q(Γ) ist Dynamik (Laplace), das kosmische q Maß.")

# ------------------------------------------------------------------ 6 Atmen als Verlustrate
Hd("6  Atmen als Verlustrate: die Zelle als eine Falte (Kandidaten vorab: (a) Γ = 0, (b) Γ = θ = 3H, (c) Γ = H)")
L_zelle = dL                                     # Zelle = d_LSS (Annahme, unabhängig von Q)
SDb = G * rho_b * L_zelle**3 / (3 * c**2)        # ΣΔs_b (Baryonen, Ω_b h² aus dem Kern)
P(f"   Zelle L = d_LSS = {L_zelle/K.Gpc:.3f} Gpc (würfelförmig, Annahme), ΣΔs_b = {SDb/K.Gpc:.4f} Gpc, Λ = {Lam/K.Gpc:.4f} Gpc (a₀ = 1,1), ΣΔs_b/Λ = {SDb/Lam:.4f}")
def q_zelle(Gam, Ds, Lam_, cell=True):
    """Rückkehr aus der Tiefe ϖ_R = √(6·Ds·Λ) auf [0, Λ] reflektierend, D = 3c·Ds (c-Schritte der Länge 6Ds)."""
    D = 3 * c * Ds; w_ = math.sqrt(6 * Ds * Lam_)
    w_ = w_ % (2 * Lam_); w_ = 2 * Lam_ - w_ if w_ > Lam_ else w_      # Tiefe auf dem Kreis 2Λ
    if Gam == 0:
        return 1.0, w_
    k = math.sqrt(Gam / D)
    return math.cosh(k * (Lam_ - w_)) / math.cosh(k * Lam_), w_
def q_irrfahrt(Gam, Ds, Lam_, N=4000):
    """zweiter Weg: Kette 0 … N (Schritt ℓ = Λ/N, Dauer h = ℓ²/2D), Überleben e^{−Γh} je Schritt, Spiegel bei N; exakt gelöst."""
    D = 3 * c * Ds; w_ = math.sqrt(6 * Ds * Lam_)
    w_ = w_ % (2 * Lam_); w_ = 2 * Lam_ - w_ if w_ > Lam_ else w_
    l = Lam_ / N; h = l * l / (2 * D); p = math.exp(-Gam * h)
    # u_j = p/2 (u_{j−1} + u_{j+1}), u_0 = 1, u_N = p·u_{N−1}
    from scipy.linalg import solve_banded
    n = N
    ab = np.zeros((3, n)); b = np.zeros(n)
    ab[1, :] = 1.0; ab[0, 1:] = -p / 2; ab[2, :-1] = -p / 2
    b[0] = p / 2
    ab[2, n - 2] = -p      # letzte Zeile: u_N − p·u_{N−1} = 0
    u = solve_banded((1, 1), ab, b)
    m = int(round(w_ / l))
    return u[m - 1] if m >= 1 else 1.0
P("   Zelle als eine Falte: Δs = ΣΔs_b (baryonisch) oder ν·ΣΔs_b (gezählt); ϖ_R = r*(Zelle) = √(6Δs·Λ), Tiefe auf dem Kreis 2Λ gefaltet.")
P("   Kandidat          Δs            ϖ_R/Λ (gefaltet)   kΛ        q_bg (cosh)   q_bg (Irrfahrt)   1 − q_bg    gegen 1 − q = 0,1571")
kand = (("(a) Γ = 0", 0.0), ("(b) Γ = 3H₀", 3 * H0), ("(c) Γ = H₀", H0))
for nm, Gam in kand:
    for lab, Ds in (("baryonisch", SDb), ("gezählt", nuk * SDb)):
        qb, wr = q_zelle(Gam, Ds, Lam)
        qi = q_irrfahrt(Gam, Ds, Lam) if Gam > 0 else 1.0
        kL = math.sqrt(Gam / (3 * c * Ds)) * Lam
        P(f"   {nm:13s}  {lab:10s} {math.sqrt(6*Ds*Lam)/Lam:5.3f} ({wr/Lam:5.3f})      {kL:7.3f}    {qb:10.6f}    {qi:10.6f}      {1-qb:8.5f}    {(1-qb)/(1-qk):7.3f}×")
P("   Monte-Carlo-Probe (b) baryonisch, Kette N = 200, 2·10⁴ Läufe:")
N = 200; Ds = SDb; D = 3 * c * Ds; w_ = math.sqrt(6 * Ds * Lam); l = Lam / N; h = l * l / (2 * D); p = math.exp(-3 * H0 * h)
m = int(round(w_ / l)); nst = 20000
pos = np.full(nst, m); alive = np.ones(nst, bool); ret = np.zeros(nst, bool)
while alive.any():
    idx = np.where(alive)[0]
    kill = rng.random(idx.size) > p
    alive[idx[kill]] = False
    idx = idx[~kill]
    st = rng.choice((-1, 1), size=idx.size)
    newp = pos[idx] + st
    newp = np.where(newp > N, 2 * N - newp, newp)
    pos[idx] = newp
    arr = idx[newp == 0]; ret[arr] = True; alive[arr] = False
P(f"     q_MC = {ret.mean():.4f} ± {math.sqrt(ret.mean()*(1-ret.mean())/nst):.4f};  Kette exakt (N = 200) {q_irrfahrt(3*H0, SDb, Lam, N=200):.4f};  cosh {q_zelle(3*H0, SDb, Lam)[0]:.4f}")
P("   Keiner der drei Kandidaten gibt ohne Anpassung 1 − q = 0,157: (a) q = 1 (ν = ∞); (b) 1 − q = 0,97 (baryonisch) bzw. 0,50 (gezählt);")
P("   (c) 0,83 bzw. 0,26. Nächster Wert (c, gezählt) liegt 1,6-fach daneben.")
P()
P("   Zeitabhängig (Γ = θ = 3H(z), x = θτ_D(z)): 1 − w = q_bg; Bedingung 1 − w(1090) ≲ 10⁻¹⁰")
P("   z        x = 3Hτ_D     1 − w: Falte ≪ Λ e^(−√x)    Zelle als Falte (baryonisch)    Maxwell")
for z in (0.0, 1.0, 10.0, 1090.0):
    a0z = A0 * K.a0_frueh_relativ(z); Lz = c**2 / (2 * a0z); tz = c / a0z
    Hz = H0 * K.E_hintergrund(z); xz = 3 * Hz * tz
    # Zelle als Falte: kϖ = √x, kΛ = √x·Λ/ϖ_R, ϖ_R = √(6ΣΔs_b Λ)
    wr = math.sqrt(6 * SDb * Lz); rat = Lz / wr
    lq = math.sqrt(xz) * (rat - 1) + math.log1p(math.exp(-2 * math.sqrt(xz) * (rat - 1))) - (math.sqrt(xz) * rat + math.log1p(math.exp(-2 * math.sqrt(xz) * rat)))
    P(f"   {z:6.0f}   {xz:10.4e}    10^{-math.sqrt(xz)/math.log(10):9.3f}             10^{lq/math.log(10):9.3f}                    {1-wM(xz):.3e}")
P("   Beide Rückkehrformen erfüllen die Bedingung bei z = 1090 um Hunderte Dekaden; heute ist 1 − w 1,5 % (Falte) bzw. 2,6 % (Zelle).")
P()
P("   Vakuumteil mit Q̃(s_bg) statt Q̃(0): |Q̃(s)| = 4F(s), F(s) = ∫_s^∞ x³n_B dx, F(0) = π⁴/15; Ω_s → Ω_s·F(s_bg)/F(0)")
F = lambda s: quad(lambda x: x**3 * math.exp(-x) / -math.expm1(-x) if x > 0 else 0.0, s, np.inf, limit=200)[0]
F0 = math.pi**4 / 15
P(f"   Probe F(0) = {F(1e-12):.6f} (π⁴/15 = {F0:.6f})")
Oms0 = K.omega_vakuum()
for nm, fak in (("(b) Γ = 3H", 3.0), ("(c) Γ = H", 1.0)):
    sb = lambda z: math.sqrt(fak * H0 * K.E_hintergrund(z) * tauD)
    rows = []
    for z in (0.0, 0.5, 1.0, 2.0):
        dz = 1e-4
        la = lambda zz: math.log(F(sb(zz)))
        # w_s = −1 − (1/3) d ln ρ_s/d ln a,  d ln a = −dz/(1+z)
        dlnF_dlna = -(la(z + dz) - la(z - dz)) / (2 * dz) * (1 + z)
        rows.append((z, sb(z), F(sb(z)) / F0, -1 - dlnF_dlna / 3))
    P(f"   {nm}: " + ";  ".join(f"z={r[0]:.1f}: s_bg={r[1]:.3f}, F/F(0)={r[2]:.4f}, w_s={r[3]:+.3f}" for r in rows))
    P(f"      Ω_s heute {Oms0:.4f} → {Oms0*rows[0][2]:.4f}")
P("   w_s < −1: F wächst, wenn H fällt; der Vakuumteil nähme zu (Phantom-Seite), gemessen ist w₀ > −1.")
P()
P("   Fixpunkt Λ = ν_kosm(Λ)·ΣΔs_b mit ν_kosm = 1/(1 − q_bg(Λ)); a₀ = c²/2Λ (τ_D = 2Λ/c); ΣΔs_b fest (Zelle = d_LSS)")
P(f"   Hand: r*(Zelle)/Λ = √(6ΣΔs_b/Λ) = √(6/ν_kosm): ν = 1+Q: {math.sqrt(6/nuk):.5f};  ν = 2π: {math.sqrt(6/(2*math.pi)):.5f};  gezählt (Δs = Λ): √6 = {math.sqrt(6):.4f}")
for nm, Gam in kand:
    g = lambda L_: (1 - q_zelle(Gam, SDb, L_)[0]) * L_ / SDb - 1
    Ls = np.geomspace(1.0001 * SDb, 200 * SDb, 4000)
    gv = np.array([g(L_) for L_ in Ls])
    roots = [brentq(g, Ls[i], Ls[i + 1]) for i in range(len(Ls) - 1) if np.sign(gv[i]) != np.sign(gv[i + 1])]
    if Gam == 0:
        P(f"   {nm}: q_bg = 1 ⇒ (1 − q)Λ = 0 < ΣΔs_b: kein Fixpunkt. Reine Kac-Lesung 1 − q = ΣΔs_b/Λ: Gleichung gilt für jedes Λ (Identität, Λ offen).")
    else:
        P(f"   {nm}: Fixpunkte Λ/ΣΔs_b = " + (", ".join(f"{r/SDb:.4f} (Λ = {r/K.Gpc:.3f} Gpc, a₀ = {c**2/(2*r)*1e10:.3f}·10⁻¹⁰, ν = {r/SDb:.3f})" for r in roots) if roots else "keiner"))
P("   (mit Γ = θ ist τ_D selbst an Λ gebunden; hier Γ fest in H₀ gerechnet, k = √(Γ/D) hängt nicht an Λ)")

# ------------------------------------------------------------------ 7
Hd("7  Ergebnis")
P("   Setzung                         folgt?                                      Zusatzannahme")
P("   (a) τ_D = 2Λ/c                  Faktor und Massenfreiheit: folgt (Kac)       X_a: Tiefe der Laplace-Form = Diffusionstiefe einer Kac-Zeit (ϖ_R² = D⟨T⟩)")
P("   (b) θ_c = 1/τ_D                 folgt (dieselbe Rückkehr, keine neue Zahl)   –")
P("       Form w = x²/(1+x²)          folgt nicht; es folgt w = 1 − e^{−√x}         (Maxwell verlangt eine Zeit und Frequenz statt Rate)")
P("   (c) 1 − q = 1/2π                folgt nicht aus Kac/Budget                   X_c: offener Bereich = ein Radiant des Kreises 2Λ")
P("   (d) Λ = ν·ΣΔs_b                 folgt als Kac-Lemma mit „kein Außen“         X_d: ν_kosm (Ω_m/Ω_b) ist die Kac-Rückkehrzeit des Bereichs ΣΔs_b")
P("   Atmen als Verlust (Γ = 3H, H)   folgt nicht: 1 − q_bg = 0,97/0,50 (3H) bzw. 0,83/0,26 (H), baryonisch/gezählt, statt 0,157; Γ = 0 gibt q = 1")
P("   Fixpunkt Λ = ν(Λ)·ΣΔs_b          Γ = 0: Identität (Λ offen); Γ = 3H₀ bzw. H₀: a₀ = 2,9 bzw. 2,4·10⁻¹⁰ (ν = 2,3 bzw. 2,9)")

# ------------------------------------------------------------------ 8
Hd("8  Unabhängige Prüfungen (von Hand, zweiter Weg)")
P(f"   1) √17,86 = 4,2261; e^(−4,2261) = 0,01461 (Tabelle {math.exp(-math.sqrt(x0)):.5f}); Maxwell 1/(1 + 17,86²) = 1/320,0 = 0,003125 (Kern {1-K.gewicht_w(3*H0):.6f})")
P(f"   2) w = 0,9: √x = ln 10 = 2,302585 → x = 5,30190; w = 0,1: √x = ln(1/0,9) = 0,105361 → x = 0,011101;")
P(f"      linear um die Umkehr Breite ∝ Δx: (5,30190 − 0,01110)/(3 − 1/3) = 5,29080/2,66667 = 1,9841 (Top-Hat: " + ", ".join(f"{breiten[z][1]/breiten[z][0]:.3f}" for z in breiten) + ")")
P(f"   3) Zellform: a₀ ∝ 1/V; 1,13/(√3/2) = 1,305; 1,13/(1/√2) = 1,598 (Tabelle Abschn. 4 aus dem Volumen)")
P(f"   4) √(6/6,3643) = √0,94276 = 0,97096 (Rechnung {math.sqrt(6/nuk):.5f})")
qc = math.cosh(math.sqrt(H0/(3*c*SDb))*(Lam-math.sqrt(6*SDb*Lam)))/math.cosh(math.sqrt(H0/(3*c*SDb))*Lam)
P(f"   5) (c) Γ = H₀ von Hand: k = √(H₀/(3cΣΔs_b)) = √(2,1843e-18/(3·2,9979e8·{SDb:.4e})) = {math.sqrt(H0/(3*c*SDb)):.4e} m⁻¹;")
P(f"      kΛ = {math.sqrt(H0/(3*c*SDb))*Lam:.4f}, kϖ_R = √(H₀τ_D) = √{H0*tauD:.4f} = {math.sqrt(H0*tauD):.4f}; q = cosh(kΛ − kϖ_R)/cosh(kΛ) = {qc:.5f}")
P(f"   6) Q = 2π − 1 = 6,283185 − 1 = 5,283185; 1/2π = 0,159155")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kac_schliessung.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
