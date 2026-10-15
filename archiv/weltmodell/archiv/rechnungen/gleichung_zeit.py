# -*- coding: utf-8 -*-
"""
Die eine Gleichung mit Zeit: Projektionen für einen Bereich, der sich mit der Rate θ im Raum dehnt.
(Rechnung, noch nicht im Modell)

Statisch (07): (∇² + ∂_ϖ²)δσ = 0, ϖ ≥ 0 innere Richtung, Ort ϖ = 0 (G(0) = 1), geschlossen bei Λ, gelesen bei ϖ_R = r*.
Die statische Gleichung legt die Zeitform nicht fest. Geprüft werden drei Zeitformen mit derselben statischen Grenze:
  [H] parabolisch (Irrfahrt der Ebenen, Schritt r₀ mit c, D = r₀c/2):  ∂_t u + θu = D(∇² + ∂_ϖ²)u
  [W] hyperbolisch (Welle in t, ϖ raumartig):                           (∂_t + θ)²u/c² = (∇² + ∂_ϖ²)u
  [T] Telegraph (beides, Irrfahrt mit Nachwirkung):                     (∂_t+θ)²u/c² + (∂_t+θ)u/D = (∇² + ∂_ϖ²)u
θ = ∇·v ist die Ausdehnungsrate des Bereichs im Raum (L); sie tritt über die Kontinuität ∂_t u + ∇·(u v) ein.
Raummode k (Kreisbahnmode k = 1/r), Laplace-Variable p der reellen Zeit:
  [H] û(ϖ_R) = exp(−ϖ_R √(k² + (p+θ)/D)) = exp(−√((p + Γ + θ)τ_D)),  Γ = Dk², τ_D = ϖ_R²/D
  [W] û(ϖ_R) = exp(−ϖ_R √(k² + (p+θ)²/c²))
1  Statik gleich, Zeitformen verschieden: q, 1 − w, Kern K̂ je Form; Masse im Gewicht; FD-Lösung der Randwertaufgabe
2  Summe der Verluste unter einer Wurzel: s² = Γτ_D + θτ_D (FD)
3  Bereich mit Kollaps: θ(t) aus der Kugel (Umkehr, Kollaps, Virial bei R_ta/2), Lösung in (t, ϖ) von [H] in reeller Zeit
   (Feynman-Kac-Integral und Crank-Nicolson), Projektion auf den Ort, gegen w(θ(t)) und Preisach w(θ̃)
4  Übergabe: Relaxationsraten der Gleichung gegen λ_m = ω_orb/N*
5  Mode ξ: Eikonal einer in 5D lichtartigen Phase, Nullmode und erste geschlossene Mode
6  Handprüfungen
Nur Kern-Konstanten (ftm_kern).
"""
import math, os, sys, io
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import erfc
from scipy.linalg import solve_banded
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import ftm_kern as K
from ftm_kern import G, c, Gyr, kpc, Gpc, M_sun, A0, OMEGA_M, gewicht_w

OUT = []
def P(s=""):
    OUT.append(s); print(s)
def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)

TD = K.tau_D()                 # s
LAM = K.schliesslaenge()       # m
H0 = K.H0_si()
TDg = TD / Gyr
H0g = H0 * Gyr

def levy(T, td=TDg):
    T = np.asarray(T, dtype=float)
    out = np.zeros_like(T)
    m = T > 0
    out[m] = np.sqrt(td / (4 * np.pi)) * T[m]**-1.5 * np.exp(-td / (4 * T[m]))
    return out

def fd_rand(alpha, X=60.0, N=60000, xr=1.0):
    """u'' = alpha·u auf [0, X] (x = ϖ/ϖ_R), u(0) = 1, u'(X) = 0 (geschlossen), Wert bei x = xr."""
    h = X / N
    n = N                                  # Unbekannte u_1 … u_N
    ab = np.zeros((3, n))
    ab[1, :] = -(2 + alpha * h * h)
    ab[0, 1:] = 1.0
    ab[2, :-1] = 1.0
    ab[0, n - 1] = 1.0                     # obere Nebendiagonale am letzten Punkt ungenutzt
    ab[2, n - 2] = 2.0                     # Neumann: u_{N+1} = u_{N−1}
    rhs = np.zeros(n); rhs[0] = -1.0
    u = solve_banded((1, 1), ab, rhs)
    i = int(round(xr / h)) - 1
    return u[i]

# ============================================================================================ 1
hdr("1  Statik gleich, Zeitformen verschieden (Milchstraße 6·10¹⁰ M☉ bei 8 kpc; Sonne; Haufen 10¹⁴ M☉ bei 300 kpc)")
cases = [("Sonne, 1000 AE", M_sun, 1000 * K.AU), ("Milchstraße, 8 kpc", 6e10 * M_sun, 8 * kpc),
         ("Haufen, 300 kpc", 1e14 * M_sun, 300 * kpc)]
P("   D = r₀c/2 (Schritt r₀ = 6Δs mit c) gegen GM/c;  τ_D = ϖ_R²/D mit ϖ_R = r* gegen c/a₀;  Γ = Dk² (k = 1/r) gegen a_N/c")
for name, M, r in cases:
    D = K.r0(M) * c / 2
    wr = K.r_stern(M)
    k = 1 / r
    aN = G * M / r**2
    P(f"   {name:20s} D = {D:.4e} (GM/c {G*M/c:.4e}) m²/s   τ_D = {wr**2/D/Gyr:.3f} Gyr (c/a₀ {TDg:.3f})"
      f"   Γ = {D*k*k:.4e} (a_N/c {aN/c:.4e}) s⁻¹")
P()
P("   Stationär (p = θ = 0): alle Formen exp(−kϖ_R) = e^{−s};  FD-Randwert u'' = s²u, u(0) = 1, geschlossen bei X = 60 ϖ_R:")
for s in (0.3, 1.0, 3.0):
    P(f"      s = {s}:  FD {fd_rand(s*s):.6f}   e^(−s) {math.exp(-s):.6f}")
P()
P("   Gewicht (k = 0, θ = 3H₀, p = 0): 1 − w je Zeitform")
th = 3 * H0
for name, M, r in cases:
    D = K.r0(M) * c / 2
    wr = K.r_stern(M)
    oneH = math.exp(-wr * math.sqrt(th / D))
    oneW = math.exp(-wr * th / c)
    oneT = math.exp(-wr * math.sqrt(th / D + th**2 / c**2))
    P(f"      {name:20s} [H] {oneH:.6e} (w = {1-oneH:.5f})   [W] w = {1-oneW:.3e} (∝ ϖ_R ∝ √M)   [T] w = {1-oneT:.5f}")
P(f"      Modellformel gewicht_w(3H₀) = {gewicht_w(th):.5f}; nur [H]/[T] geben θ_c = D/ϖ_R² = a₀/c für jede Masse.")
P()
P("   Kern bei Bahnfrequenz, |û(iΩ)|/û(0) (Milchstraße, v = 220 km/s):")
M, r = 6e10 * M_sun, 8 * kpc
D = K.r0(M) * c / 2; wr = K.r_stern(M); k = 1 / r; Om = 220e3 / r
def uH(p): return np.exp(-wr * np.sqrt(k * k + p / D + 0j))
def uW(p): return np.exp(-wr * np.sqrt(k * k + (p / c)**2 + 0j))
def uT(p): return np.exp(-wr * np.sqrt(k * k + p / D + (p / c)**2 + 0j))
for nm, fn in (("[H]", uH), ("[W]", uW), ("[T]", uT)):
    P(f"      {nm}: {abs(fn(1j*Om))/abs(fn(0)):.4e}")
P(f"      [W]: Verzögerung ϖ_R/c = {wr/c/K.yr:.0f} Jahre (Laufzeit in ϖ mit c); [H]: Lévy-Kern, Modus τ_D/6 = {TDg/6:.1f} Gyr")
P(f"      Übergang Welle/Irrfahrt bei |p| ~ c²/D = 2c/r₀ = {2*c/K.r0(M):.3e} s⁻¹ (Milchstraße), {2*c/K.r0(M_sun):.3e} s⁻¹ (Sonne)")

# ============================================================================================ 2
hdr("2  Summe der Verluste unter einer Wurzel (FD, [H]): u(ϖ_R) = exp(−√(Γτ_D + θτ_D))")
for sG2, th2 in ((1.0, 0.0), (0.0, 17.86), (1.0, 17.86), (0.1, 1.0)):
    fd = fd_rand(sG2 + th2)
    P(f"   Γτ_D = {sG2:5.2f}, θτ_D = {th2:6.2f}:  FD {fd:.6e}   exp(−√Σ) {math.exp(-math.sqrt(sG2+th2)):.6e}"
      f"   Produkt q·(1 − w) {math.exp(-math.sqrt(sG2))*math.exp(-math.sqrt(th2)):.6e}")
P("   Die Raumrichtung (∇² → Γ = Dk²) und die Zeitrichtung (∂_t + ∇·v → θ) addieren sich unter einer Wurzel, nicht als Produkt.")

# ============================================================================================ 3
hdr("3  Bereich mit Kollaps: Lösung in (t, ϖ) von [H] in reeller Zeit gegen w(θ) und Preisach w(θ̃)")
OM, OL = OMEGA_M, 1 - OMEGA_M
def Hbg(a): return math.sqrt(OM * a**-3 + OL)            # H/H₀, ohne Strahlung (wie die Kugel)
def t_of_a(a):
    return quad(lambda x: 1 / (x * Hbg(x)), 0, a, limit=200)[0] / H0g   # Gyr
T0 = t_of_a(1.0)
a_i = 1e-3
t_i = t_of_a(a_i)

def kugel(delta_i, t_end):
    """R̈ = −(Ω_m/2)/R² + Ω_L R (Einheiten H₀ = 1, R = a für δ = 0)."""
    R0 = a_i * (1 - delta_i / 3)
    V0 = Hbg(a_i) * a_i * (1 - 2 * delta_i / 3)
    def rhs(t, y): return [y[1], -0.5 * OM / y[0]**2 + OL * y[0]]
    def hit(t, y): return y[0] - 1e-4 * R0 / a_i * 1e-2
    hit.terminal = True; hit.direction = -1
    sol = solve_ivp(rhs, [t_i * H0g, t_end * H0g], [R0, V0], rtol=1e-10, atol=1e-13,
                    events=hit, dense_output=True, max_step=0.002)
    return sol

def kollapszeit(delta_i):
    sol = kugel(delta_i, 60.0)
    return sol.t_events[0][0] / H0g if len(sol.t_events[0]) else np.inf

def bereich(z_c):
    t_c = t_of_a(1 / (1 + z_c))
    d = brentq(lambda x: kollapszeit(x) - t_c, 1e-4, 0.5, xtol=1e-12)
    sol = kugel(d, t_c * 0.99999)
    tt = np.linspace(t_i * H0g, sol.t[-1], 400001)
    R, Rd = sol.sol(tt)
    i_ta = int(np.argmax(R))
    R_ta = R[i_ta]
    j = i_ta + int(np.argmax(R[i_ta:] <= R_ta / 2))
    t_g = tt / H0g
    Rv = R.copy(); th = 3 * Rd / R * H0g
    Rv[j:] = R_ta / 2; th[j:] = 0.0
    return dict(z_c=z_c, t_c=t_c, t_ta=t_g[i_ta], t_vir=t_g[j], d=d, t=t_g[:j + 1], R=Rv[:j + 1], th=th[:j + 1], R_vir=R_ta / 2)

def lnV_fn(B, extra=None):
    """ln V(t) des Bereichs (Gyr): vor t_i mitbewegt (a ∝ t^{2/3}), nach dem Virial konstant, optional Zusatzgeschichte."""
    t, R = B["t"], B["R"]
    lnR = np.log(R)
    def f(tq):
        tq = np.asarray(tq, dtype=float)
        out = np.empty_like(tq)
        m0 = tq < t[0]
        out[m0] = 3 * (lnR[0] + (2 / 3) * np.log(np.maximum(tq[m0], 1e-30) / t[0]))
        m1 = (~m0) & (tq <= t[-1])
        out[m1] = 3 * np.interp(tq[m1], t, lnR)
        m2 = tq > t[-1]
        out[m2] = 3 * math.log(B["R_vir"]) if extra is None else extra(tq[m2])
        return out
    return f

def eins_minus_w_H(lnV, t, n=400000):
    """[H] in reeller Zeit, Feynman-Kac: 1 − w(t) = ∫₀^t f(T)·V(t−T)/V(t) dT (vor dem Anfang V → 0)."""
    x = np.linspace(0, 1, n + 1)[1:]
    T = t * x**2                          # dichter bei kleinen T
    wgt = 2 * t * x
    val = levy(T) * np.exp(lnV(t - T) - lnV(np.array([t]))[0]) * wgt
    return np.trapezoid(np.concatenate([[0.0], val]), np.concatenate([[0.0], x]))

def w_quasi(th_per_gyr):
    return gewicht_w(max(th_per_gyr, 0.0) / Gyr)

# Crank-Nicolson als zweiter Weg: u_t = u_xx − θτ_D u (x = ϖ/ϖ_R, Zeit in τ_D), u(0) = 1, u_x(X) = 0
def cn_loesung(theta_fn, t_start, t_list, X=40.0, N=2000, dt_g=0.002):
    h = X / N
    x = np.linspace(h, X, N)
    th0 = max(theta_fn(t_start), 0.0) * TDg
    u = np.exp(-x * math.sqrt(th0))
    res = {}
    t = t_start
    targets = sorted(t_list)
    k = 0
    while k < len(targets):
        dt = min(dt_g, targets[k] - t)
        tm = t + dt / 2
        a = theta_fn(tm) * TDg
        r = (dt / TDg) / h**2
        # Matrix A = I − ½dt·L, L u = u_xx − a u
        main = np.full(N, 1 + r + 0.5 * dt / TDg * a)
        off = np.full(N - 1, -0.5 * r)
        ab = np.zeros((3, N)); ab[1] = main; ab[0, 1:] = off; ab[2, :-1] = off
        ab[2, N - 2] = -r                # Neumann am Ende
        Lu = np.empty(N)
        Lu[0] = (1 - 2 * u[0] + u[1]) / h**2
        Lu[1:-1] = (u[:-2] - 2 * u[1:-1] + u[2:]) / h**2
        Lu[-1] = (2 * u[-2] - 2 * u[-1]) / h**2
        rhs = u + 0.5 * dt / TDg * (Lu - a * u)
        rhs[0] += 0.5 * r * 1.0          # Randwert u(0) = 1 (beide Hälften)
        u = solve_banded((1, 1), ab, rhs)
        t += dt
        if abs(t - targets[k]) < 1e-12:
            res[targets[k]] = u[int(round(1.0 / h)) - 1]
            k += 1
    return res

# Hintergrund
_ag = np.geomspace(1e-8, 1.0, 4000)
_tg = np.array([t_of_a(a) for a in _ag])
def a_of_t(t): return float(np.exp(np.interp(math.log(max(t, 1e-30)), np.log(_tg), np.log(_ag))))
def lnV_bg(tq):
    tq = np.atleast_1d(np.asarray(tq, dtype=float))
    return 3 * np.interp(np.log(np.maximum(tq, 1e-30)), np.log(_tg), np.log(_ag))
P(f"   Weltalter (ohne Strahlung) t₀ = {T0:.3f} Gyr, τ_D = {TDg:.3f} Gyr")
wbg = 1 - eins_minus_w_H(lnV_bg, T0)
P(f"   Hintergrund heute: [H] in reeller Zeit (rückwärts, wirkliches a(t)) w = {wbg:.5f};  quasistatisch w(3H₀) = {w_quasi(3*H0g):.5f};"
  f"  Obergrenze bei wachsendem V: 1 − w ≤ P(T < t₀) = erfc(√(τ_D/4t₀)) = {erfc(math.sqrt(TDg/(4*T0))):.5f}")
P()
P("   Kugel mit Kollaps bei z_c; Virial bei R_ta/2 (danach θ = 0).  w_Q = w(θ(t)), w_P = w(θ̃) (Preisach), w_H = [H] reell (Feynman-Kac).")
P("   z_c   t_ta   t_vir  t_c    | t(w_P=0,9) t(w_P=0,1) Δt [Mio. J.] | w_H bei: t_ta/2   t_ta    t_vir   t₀      | max(1−w_H) bis t₀")
rows = {}
for z_c in (3.4, 2.0, 1.0, 0.5):
    B = bereich(z_c)
    rows[z_c] = B
    lnV = lnV_fn(B)
    tt, th = B["t"], B["th"]
    wP = np.array([w_quasi(x) for x in np.minimum.accumulate(th)])
    i9 = np.argmax(wP <= 0.9); i1 = np.argmax(wP <= 0.1)
    pts = [B["t_ta"] / 2, B["t_ta"], B["t_vir"], T0]
    wH = [1 - eins_minus_w_H(lnV, t) for t in pts]
    grid = np.linspace(B["t_ta"] * 0.3, T0, 60)
    mx = max(eins_minus_w_H(lnV, t) for t in grid)
    B["wH"] = dict(zip(("ta2", "ta", "vir", "t0"), wH))
    P(f"   {z_c:3.1f}  {B['t_ta']:5.3f}  {B['t_vir']:5.3f}  {B['t_c']:5.3f} | {tt[i9]:7.4f}   {tt[i1]:7.4f}   {1e3*(tt[i1]-tt[i9]):7.1f}      |"
      f"        {wH[0]:.4f}  {wH[1]:.4f}  {wH[2]:.4f}  {wH[3]:.4f} | {mx:.4f}")
P("   (Modell 03: w_P von 0,9 auf 0,1 in 14 / 45 / 143 / 309 Mio. J. bei z_c = 3,4 / 2 / 1 / 0,5)")

# zweiter Weg: Crank-Nicolson für z_c = 2
B = rows[2.0]
def th_fn(t):
    if t < B["t"][0]:
        return 2 / t        # 3·(2/3)/t im mitbewegten Anfang
    if t > B["t"][-1]:
        return 0.0
    return float(np.interp(t, B["t"], B["th"]))
pts = [B["t_ta"], B["t_vir"], T0]
cn = cn_loesung(th_fn, 0.05, pts)
lnV = lnV_fn(B)
P()
P("   Zweiter Weg (Crank-Nicolson auf (t, ϖ), Start t = 0,05 Gyr im stationären Zustand), z_c = 2:")
for t in pts:
    P(f"      t = {t:6.3f} Gyr:  CN 1 − w = {cn[t]:.5f}   Feynman-Kac {eins_minus_w_H(lnV, t):.5f}")

# Nachgeschichte eines gebundenen Bereichs (z_c = 2): (i) Schwingung, (ii) Wiederausdehnung mit 3H ab 6 Gyr
P()
P("   Nachgeschichte (z_c = 2, virialisiert): Preisach hält w = 0 in beiden Fällen.")
Rv3 = 3 * math.log(B["R_vir"])
om = math.sqrt(200 / 2) * H0g * 1.0
def ex_osc(tq): return Rv3 + np.log1p(0.3 * np.sin(om * (tq - B["t_vir"])))
def ex_re(tq):
    t1 = 6.0
    return np.where(tq < t1, Rv3, Rv3 + (lnV_bg(tq) - lnV_bg(np.full_like(tq, t1))))
for nm, ex in (("Schwingung ΔV/V = 0,3 mit ω_orb(Δ=200)", ex_osc), ("Wiederausdehnung mit 3H ab 6 Gyr", ex_re)):
    f = lnV_fn(B, ex)
    vals = [(t, 1 - eins_minus_w_H(f, t)) for t in (5.0, 6.0, 8.0, 10.0, T0)]
    P("      " + nm + ":  " + "   ".join(f"t={t:.1f}: w_H={w:.4f}" for t, w in vals))
P(f"      quasistatisch bei Wiederausdehnung heute: w(3H₀) = {w_quasi(3*H0g):.4f}")
thmax = 0.3 * om / math.sqrt(1 - 0.09)
P(f"      quasistatisch ohne Gedächtnis bei der Schwingung: θ zwischen ±{thmax:.3f}/Gyr ⇒ w(θ(t)) springt zwischen 0 und {w_quasi(thmax):.4f}")

# ============================================================================================ 4
hdr("4  Übergabe: Relaxationsraten der Gleichung [H] gegen λ_m = ω_orb/N*")
M, r = 6e10 * M_sun, 8 * kpc
D = K.r0(M) * c / 2; aN = G * M / r**2; Om = 220e3 / r
lam_min = D * math.pi**2 / (4 * LAM**2)
P(f"   Milchstraße 8 kpc: Γ = {aN/c:.3e}, 1/τ_D = {1/TD:.3e}, Dπ²/4Λ² = {lam_min:.3e} s⁻¹;  λ_m = (0,01–0,13)Ω = {0.01*Om:.2e}–{0.13*Om:.2e} s⁻¹")
P(f"   langsamste Rate der Gleichung am Ort (θ = 0): Γ + Dπ²/4Λ² = {aN/c+lam_min:.3e} s⁻¹ = λ_m/{0.01*Om/(aN/c+lam_min):.0f} … λ_m/{0.13*Om/(aN/c+lam_min):.0f}")
for Dl in (200, 500):
    w_orb, l1, l2 = K.uebergabezeit(Dl)
    P(f"   Bereich Δ = {Dl}: 1/ω_orb = {w_orb/Gyr:.3f} Gyr, 1/λ_m = {l1/Gyr:.2f}–{l2/Gyr:.2f} Gyr;  Zeitskalen der Gleichung: τ_D = {TDg:.1f}, Modus τ_D/6 = {TDg/6:.1f} Gyr")
P("   Die Gleichung hat keine Rate ω_orb/N*: ihre Raten sind Γ, θ und D/Länge² in ϖ (alle ≤ 1/τ_D-Skala am Ort ϖ_R).")

# ============================================================================================ 5
hdr("5  Mode ξ: Eikonal einer in (t, x, ϖ) lichtartigen Phase S = ξ(t, x) + p_ϖ·ϖ (nur Zeitform [W]/[T])")
P("   (∂S)²₅ = (∂ξ)²₄ + p_ϖ² = 0  ⇒  (∂ξ)² = −p_ϖ²;  Flusslinien = Geodäten, druckfrei (Kaluza-Klein-Bild).")
P(f"   Neumann-Moden auf [0, Λ]: p_n = nπ/Λ;  n = 0 (Nullmode): (∂ξ)² = 0 (lichtartig), nicht −1;  n = 1: p₁ = π/Λ = {math.pi/LAM:.3e} m⁻¹, c·p₁ = {c*math.pi/LAM:.3e} s⁻¹ = {c*math.pi/LAM/H0:.3f} H₀")
P("   Mit G(0) = 1 und beschränkter Lösung (Halbgerade) gibt es stationär keine reelle ϖ-Welle (nur e^{−kϖ}): dort ist (∂ξ)² = +k² (raumartig).")
P("   In [H] gibt es kein Eikonal; in [T] gilt es erst oberhalb |p| ~ 2c/r₀ (Abschnitt 1).")

# ============================================================================================ 6
hdr("6  Handprüfungen")
P(f"   τ_D = c/a₀ = 2,99792458e8/1,1e-10 = 2,72539e18 s / 3,15576e16 = 86,362 Gyr (Kern {TDg:.3f})")
x = math.sqrt(TDg / (4 * T0))
P(f"   Obergrenze 1 − w in reeller Zeit: erfc(√(τ_D/4t₀)), √(86,362/(4·{T0:.3f})) = {x:.4f}")
P("      erf(x) = 2/√π·Σ (−1)ⁿ x^{2n+1}/(n!(2n+1))  (Maclaurin, Glied für Glied)")
terms = [(-1)**n * x**(2*n+1) / (math.factorial(n) * (2*n+1)) for n in range(14)]
P("      Glieder: " + ", ".join(f"{t_:+.5f}" for t_ in terms[:8]) + ", …")
P(f"      erfc = 1 − 2/√π·{sum(terms):.6f} = {1 - 2/math.sqrt(math.pi)*sum(terms):.5f}   scipy {erfc(x):.5f}")
P("   Milchstraße, [W]-Gewicht: w ≈ ϖ_R·3H₀/c = r*·3H₀/c;  r* = 8,35 kpc·√(1,2/1,1) ≈ 8,72 kpc = 2,69e20 m;")
P(f"      3H₀/c = 3·2,1843e-18/2,9979e8 = 2,186e-26 m⁻¹ ⇒ w ≈ 5,9e-6   (Code {1-math.exp(-K.r_stern(6e10*M_sun)*3*H0/c):.2e})")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gleichung_zeit.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
