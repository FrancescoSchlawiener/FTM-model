# -*- coding: utf-8 -*-
"""
Ebenengrenze der Halbgeraden bei druckgestützten Systemen: Ordnet sie Pal 14, NGC 2419, Crater II, Draco, Fornax, DF2, DF4?

1  Formel (Notiz 53, Abschnitt 3b; halbgerade_hierarchie.py): D·G'' = Γ·G auf ϖ ≥ 0, G(0) = 1, G beschränkt.
   Unter der Tiefe ℓ (Glättungslänge im Elternfeld, ℓ = d·(M/3M_p)^{1/3}) ist D = GM/c, darüber D·μ, μ = M_p/M; Γ je Zeit gleich;
   Fluss D·G′ und G stetig ⇒ R = (1 − √μ)/(1 + √μ);  q = e^{−kL}(1 + R e^{−2k(ℓ−L)})/(1 + R e^{−2kℓ}),  k = s/L,  L = ϖ_R = r*.
   Zweiter Weg: Schießen von oben (abklingende Lösung über ℓ, Integration nach unten mit solve_ivp), Grenzfälle von Hand
   (μ = 1: R = 0, e^{−s};  μ → ∞: G(ℓ) = 0, q = sinh(k(ℓ − L))/sinh(kℓ)).
2  Sieben Systeme (Eingaben wie glaettung_stetig.py SYS, [G]); Plummer b = R_e, isotrope Jeans-Gleichung, σ im Fenster R < R_e,
   Monopol ⟨ν(|g_i + g_e|)(g_i + g_e cos θ)⟩ wie gedaechtnis_umgebung.py; mit Ebenengrenze ν → ν_b = 1/(1 − q_b).
   Varianten: Tiefe ϖ_J = 0,75·ℓ (Poisson-Kern, Notiz 70), R = (1 − μ)/(1 + μ) (Γ = D/r² je Ebene), ℓ aus der wirksamen Gezeit.
3  Schranken: ν_b ≥ 1 (G > 0) und ν_b ≤ ν für μ ≥ 1 (R ≤ 0).
4  Verlangte Tiefe für Fornax und DF2 und der zugehörige Abstand.
[A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.  a₀ = A0 aus ftm_kern (1,1·10⁻¹⁰ m/s²).
"""
import io, os, sys, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ftm_kern import *

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(all="ignore")
out = []
def P(s=""): out.append(s)
def H(t): P(); P("=" * 110); P(t); P("=" * 110)

# ------------------------------------------------------------------ Einheiten: kpc, km/s; Felder in a₀
A0K = A0 * kpc / 1e6                     # a₀ in (km/s)²/kpc
GK = G * M_sun / kpc / 1e6               # G in kpc (km/s)²/M☉
GMu = lambda M: GK * M / A0K             # GM in kpc²·a₀  (r* = √GMu in kpc)
kms = lambda s2: np.sqrt(np.maximum(s2, 0) * A0K)

# ------------------------------------------------------------------ 1 Formel
def q_b(s, x, mu, Rform="sqrt"):
    """Rückkehranteil mit Ebenengrenze. s = √(|g|/a₀), x = ℓ/L (Tiefe der Grenze in Einheiten von ϖ_R = r*), μ = M_p/M."""
    s = np.asarray(s, float)
    R = (1 - math.sqrt(mu)) / (1 + math.sqrt(mu)) if Rform == "sqrt" else (1 - mu) / (1 + mu)
    if x >= 1:
        return np.exp(-s) * (1 + R * np.exp(-2 * s * (x - 1))) / (1 + R * np.exp(-2 * s * x))
    kp = 1 / math.sqrt(mu) if Rform == "sqrt" else 1.0      # k' / k über der Grenze
    Gl = np.exp(-s * x) * (1 + R) / (1 + R * np.exp(-2 * s * x))
    return Gl * np.exp(-kp * s * (1 - x))
def nu_b(s, x, mu, Rform="sqrt"): return 1 / (1 - q_b(s, x, mu, Rform))

H("1   Formel und zweiter Weg")
P("   D·G'' = Γ·G, G(0) = 1, G beschränkt; unter ϖ = ℓ: D, darüber μD (Faltungen je Ebene additiv, D_Ebene = GM_Ebene/c), Γ gleich.")
P("   Stetig: G und Fluss D·G′ (Schritte mit c). Unten G = A e^{−kϖ} + B e^{kϖ}, oben C e^{−kϖ/√μ}:")
P("   B/A = R e^{−2kℓ}, R = (1 − √μ)/(1 + √μ);  q = G(L) = e^{−kL}(1 + R e^{−2k(ℓ−L)})/(1 + R e^{−2kℓ}),  k = s/L,  L = ϖ_R = r*")
P("   Tiefe der Grenze: ℓ = d·(M/3M_p)^{1/3} (ℓ³ = GM/T, T = 3GM_p/d³, Elternmasse als Punkt). Mit s = r*/r (Punktmasse) ist kℓ = ℓ/r.")
P("   Tief (s ≪ 1): ν_b/ν = (1 + R e^{−2ℓ/r})/(1 − R e^{−2ℓ/r});  R → −1: ν_b/ν → tanh(ℓ/r)")
def shoot(s, x, mu):
    """zweiter Weg: oben abklingend, von ϖ = x nach unten integrieren (L = 1, D = 1 unten), normieren auf G(0) = 1"""
    k = s; kp = k / math.sqrt(mu)
    G0, F0 = 1.0, -mu * kp * 1.0                              # G(ℓ⁺) = 1, Fluss μD·G′ = −μ k′
    sol = solve_ivp(lambda w, y: [y[1], k * k * y[0]], [x, 0.0], [G0, F0 / 1.0], rtol=1e-12, atol=1e-14, dense_output=True)
    return sol.sol(1.0)[0] / sol.sol(0.0)[0]
P("   Prüfung Formel gegen Schießen (solve_ivp, unabhängig vom Finite-Volumen-Code in halbgerade_hierarchie.py):")
P("       s      ℓ/L     μ        q Formel       q Schießen     rel. Abw.")
mx = 0
for s, x, mu in ((0.7, 5.0, 100.0), (0.1, 3.0, 600.0), (0.05, 30.0, 1e4), (1.5, 2.0, 4.0), (0.3, 1.2, 1e6)):
    qa = float(q_b(s, x, mu)); qn = shoot(s, x, mu); mx = max(mx, abs(qn / qa - 1))
    P(f"   {s:6.2f} {x:7.1f} {mu:8.0f}   {qa:.10f}   {qn:.10f}   {qn/qa-1:+.1e}")
P(f"   max. Abweichung {mx:.1e}")
P("   Grenzfälle von Hand: μ = 1 ⇒ R = 0 ⇒ q = e^{−s} (keine Grenze);  μ → ∞ ⇒ R = −1, G(ℓ) = 0 (Spiegelladung mit Vorzeichen −):")
s, x = 0.4, 3.0
P(f"      q(μ = 1) = {float(q_b(s, x, 1.0)):.8f}, e^{{−s}} = {math.exp(-s):.8f};  q(μ = 10¹²) = {float(q_b(s, x, 1e12)):.8f},"
  f" sinh(k(ℓ − L))/sinh(kℓ) = {math.sinh(s*(x-1))/math.sinh(s*x):.8f}  (s = 0,4, ℓ = 3L)")

H("2   Sieben Systeme: Abstand, ℓ, Tiefe ℓ/ϖ_R, R, Verstärkung, σ")
SYS = [  # Name, M [M☉], R_e [kpc], Wirt [M☉], Abstand [kpc], σ gemessen, Fehler  (wie glaettung_stetig.py SYS, [G])
    ("Pal 14",    2e4,   0.0247, 6e10, 71.0,  0.38, 0.12),
    ("NGC 2419",  9e5,   0.019,  6e10, 90.0,  4.1,  0.5),
    ("Crater II", 3.2e5, 1.066,  6e10, 117.5, 2.7,  0.3),
    ("Draco",     5.8e5, 0.221,  6e10, 76.0,  9.1,  1.2),
    ("Fornax",    4.0e7, 0.71,   6e10, 147.0, 11.7, 0.9),
    ("DF2",       2.0e8, 2.2,    1e11, 80.0,  8.5,  2.3),
    ("DF4",       1.5e8, 1.6,    1e11, 165.0, 4.2,  2.2)]
MU, WMU = np.polynomial.legendre.leggauss(200)
def plummer(M, Re, n=3000):
    b = Re; GM = GMu(M); r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, n)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5; Grho = 3 * GM / (4 * np.pi * b**3) * (1 + (r / b)**2)**-2.5
    return r, GMenc / r**2, Grho
def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return float(kms(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r)))
def mono(gi, ge, nufn):
    """Monopol ⟨ν(|g|)(g_i + g_e cos θ)⟩ über die Kugel (Satz von Gauß, wie gedaechtnis_umgebung.py)"""
    x = np.sqrt(np.maximum(gi[:, None]**2 + ge**2 + 2 * gi[:, None] * ge * MU[None, :], 1e-300))
    return 0.5 * np.sum(WMU * nufn(np.sqrt(x)) * (gi[:, None] + ge * MU[None, :]), axis=1)
nu0 = lambda s: 1 / (-np.expm1(-s))
def sigma(M, Re, ge, x=None, mu=None, Rform="sqrt"):
    r, gi, rho = plummer(M, Re)
    if x is None: return jeans(r, rho, mono(gi, ge, nu0), Re)
    return jeans(r, rho, mono(gi, ge, lambda s: nu_b(s, x, mu, Rform)), Re)
def ell_of(M, Mh, D, T=None):
    if T is None: T = 3 * GMu(Mh) / D**3
    return (GMu(M) / T)**(1 / 3)
def T_eff(Mh, D):
    """Gezeit des verstärkten Wirtsfelds √(3/2)|E^TF| = |g′ − g/D| (wie glaettung_stetig.py tides)"""
    GMh = GMu(Mh); gD = lambda y: nu0(math.sqrt(GMh / y**2)) * GMh / y**2; h = 1e-4 * D
    return abs((gD(D + h) - gD(D - h)) / (2 * h) - gD(D) / D)

P("   Eingaben [G]: M (Υ = 2), R_e, Wirt Milchstraße 6·10¹⁰ bzw. NGC 1052 10¹¹ M☉ [A], Abstand; r½ = 1,305·R_e (Plummer)")
P(f"   {'System':10s} {'d[kpc]':>7s} {'μ':>9s} {'R':>8s} {'r*[pc]':>7s} {'ℓ[kpc]':>7s} {'ℓ/ϖ_R':>7s} {'ℓ/r½':>6s} {'g_e/a₀':>7s} "
  f"{'ν_b/ν(r½)':>9s}   {'σ_N':>6s} {'isol.':>6s} {'bisher':>6s} {'Grenze':>6s} {'Gr.+iso':>7s}  {'gemessen':>11s}")
ROW = {}
for nm, M, Re, Mh, D, so, er in SYS:
    GM = GMu(M); rs = math.sqrt(GM); mu = Mh / M; R = (1 - math.sqrt(mu)) / (1 + math.sqrt(mu))
    ell = ell_of(M, Mh, D); x = ell / rs; rh = 1.305 * Re; ge = GMu(Mh) / D**2
    gih = GM * rh / (rh**2 + Re**2)**1.5; sh = math.sqrt(math.hypot(gih, ge))
    ratio = float(nu_b(sh, x, mu) / nu0(sh))
    r, gi, rho = plummer(M, Re); sN = jeans(r, rho, gi, Re)
    siso = sigma(M, Re, 0.0); sF = sigma(M, Re, ge); sB = sigma(M, Re, ge, x, mu); sBi = sigma(M, Re, 0.0, x, mu)
    ROW[nm] = dict(M=M, Re=Re, Mh=Mh, D=D, so=so, er=er, sN=sN, siso=siso, sF=sF, sB=sB, x=x, mu=mu, ell=ell, rs=rs, ge=ge, rh=rh)
    P(f"   {nm:10s} {D:7.1f} {mu:9.3g} {R:8.4f} {rs*1e3:7.1f} {ell:7.3f} {x:7.1f} {ell/rh:6.2f} {ge:7.4f} {ratio:9.5f}   "
      f"{sN:6.2f} {siso:6.2f} {sF:6.2f} {sB:6.2f} {sBi:7.2f}  {so:5.2f} ± {er:<4.2f}")
P("   σ in km/s. bisher = Regel mit äußerem Feld (Monopol); Grenze = dieselbe mit Ebenengrenze; Gr.+iso = Ebenengrenze ohne äußeres Feld.")
P("   ν_b/ν(r½): Verhältnis der Verstärkung am Halbmassenradius (Gesamtfeld im Argument).")
P()
P("   Varianten (σ mit äußerem Feld und Grenze, km/s): Tiefe 0,75·ℓ (Poisson) / R = (1 − μ)/(1 + μ) / ℓ aus wirksamer Gezeit / alle drei")
for nm, d in ROW.items():
    x = d["x"]; mu = d["mu"]; xe = ell_of(d["M"], d["Mh"], d["D"], T_eff(d["Mh"], d["D"])) / d["rs"]
    v1 = sigma(d["M"], d["Re"], d["ge"], 0.75 * x, mu); v2 = sigma(d["M"], d["Re"], d["ge"], x, mu, "lin")
    v3 = sigma(d["M"], d["Re"], d["ge"], xe, mu); v4 = sigma(d["M"], d["Re"], d["ge"], 0.75 * xe, mu, "lin")
    d["xe"] = xe
    P(f"   {nm:10s} bisher {d['sF']:6.2f}:  {v1:6.2f} / {v2:6.2f} / {v3:6.2f} (ℓ_w = {xe*d['rs']:.3f} kpc) / {v4:6.2f}   gemessen {d['so']}")
P()
P("   Verlangte Richtung (σ_gem/σ_bisher)² gegen die Wirkung der Grenze (σ_Grenze/σ_bisher)², nach ℓ/r½ geordnet:")
for nm in sorted(ROW, key=lambda n: ROW[n]["ell"] / ROW[n]["rh"]):
    d = ROW[nm]
    P(f"      {nm:10s} ℓ/r½ = {d['ell']/d['rh']:6.2f}   verlangt {(d['so']/d['sF'])**2:6.3f}   Grenze {(d['sB']/d['sF'])**2:7.4f}")

H("3   Schranken: unter Newton, über isoliert")
P("   G ist auf der ganzen Halbgeraden positiv (Γ > 0, Maximumprinzip), also 0 < q_b < 1 und ν_b > 1: nie unter Newton.")
P("   Für μ ≥ 1 ist R ≤ 0 und q_b ≤ e^{−s}: ν_b ≤ ν, die Grenze hebt nie über den Wert ohne Grenze. Heben (R > 0) verlangte μ < 1.")
ss = np.logspace(-3, 1.5, 400); mn = 1e9; mxr = 0
for x in (1.0, 1.5, 3, 10, 100):
    for mu in (1.0, 10.0, 1e3, 1e6, 1e12):
        nb = nu_b(ss, x, mu); mn = min(mn, float(nb.min())); mxr = max(mxr, float((nb / nu0(ss)).max()))
P(f"   Gitter s = 10⁻³…30, ℓ/L = 1…100, μ = 1…10¹²: min ν_b = {mn:.6f} (> 1), max ν_b/ν = {mxr:.6f} (≤ 1)")
for nm in ("Pal 14", "NGC 2419", "DF4", "Draco"):
    d = ROW[nm]
    P(f"   {nm:9s}: σ_N = {d['sN']:.2f}, mit Grenze {d['sB']:.2f} (≥ σ_N), gemessen {d['so']};  isoliert {d['siso']:.2f}, Grenze allein {sigma(d['M'], d['Re'], 0.0, d['x'], d['mu']):.2f}")
P("   Selbst R = −1 (Grenze absorbierend) und Tiefe ℓ → L (Grenze an ϖ_R): σ bleibt ≥ σ_N:")
for nm in ("Pal 14", "NGC 2419", "DF4"):
    d = ROW[nm]
    P(f"      {nm:9s} ℓ = 1,0·ϖ_R, μ = 10¹²: σ = {sigma(d['M'], d['Re'], d['ge'], 1.0, 1e12):.3f} km/s  (σ_N = {d['sN']:.3f}, gemessen {d['so']})")

H("4   Verlangte Tiefe (Fornax, DF2: gemessen über σ_N, unter der Regel)")
for nm in ("Fornax", "DF2", "NGC 2419"):
    d = ROW[nm]
    if d["so"] <= d["sN"]:
        P(f"   {nm}: gemessen unter σ_N, keine Tiefe reicht."); continue
    f = lambda lx: sigma(d["M"], d["Re"], d["ge"], math.exp(lx), d["mu"]) - d["so"]
    xr = math.exp(brentq(f, math.log(1.0001), math.log(d["x"] * 10), xtol=1e-6))
    Dn = d["D"] * xr / d["x"]
    P(f"   {nm}: verlangt ℓ = {xr*d['rs']:.3f} kpc = {xr:.1f}·ϖ_R = {xr*d['rs']/d['rh']:.2f}·r½ (Formel: {d['ell']:.2f} kpc);"
      f" bei festem Wirt Abstand {Dn:.1f} kpc statt {d['D']} (g_e dann {GMu(d['Mh'])/Dn**2:.3f} a₀)")

H("5   Von Hand (zweiter Weg)")
d = ROW["Crater II"]
P(f"   Crater II: μ = 6·10¹⁰/3,2·10⁵ = {6e10/3.2e5:.4g}, √μ = {math.sqrt(6e10/3.2e5):.2f}, R = (1 − 433,0)/(1 + 433,0) = {(1-433.0)/(1+433.0):.5f}")
P(f"      ℓ = 117,5·(3,2·10⁵/1,8·10¹¹)^{{1/3}} = 117,5·{(3.2e5/1.8e11)**(1/3):.5f} = {117.5*(3.2e5/1.8e11)**(1/3):.3f} kpc (Skript {d['ell']:.3f})")
Ph = 6e10 / 3.2e5
P("   Draco, tief, Punktmasse bei r½ = 0,288 kpc, ℓ = 1,12 kpc: e^{−2ℓ/r} = e^{−7,78} = 4,2·10⁻⁴, ν_b/ν ≈ 1 − 2·|R|·4,2·10⁻⁴ = 0,99917")
xD = ROW["Draco"]["ell"] / ROW["Draco"]["rh"]; RD = (1 - math.sqrt(6e10/5.8e5)) / (1 + math.sqrt(6e10/5.8e5)); e = math.exp(-2 * xD)
P(f"      Skript: ℓ/r½ = {xD:.3f}, (1 + R e)/(1 − R e) = {(1+RD*e)/(1-RD*e):.6f}")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ebenengrenze_zwerge.txt"), "w", encoding="utf-8").write(txt)
print(txt)
