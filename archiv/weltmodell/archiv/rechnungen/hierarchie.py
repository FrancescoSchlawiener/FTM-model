# -*- coding: utf-8 -*-
"""
FTM - Einzeltrichter-/Hierarchie-Variante.  Rechnungen zu forschung/13-einzeltrichter-hierarchie.md
Aufruf:  python hierarchie.py > hierarchie.txt        (benoetigt nur numpy)

Varianten (Bezeichnungen wie in der Notiz)
  G   Gesamtfeld (QUMOND-Grenzfall von FTM-AS, Notiz 06): nu vom Newtonschen Gesamtfeld
  E1  Einzeltrichter mit eigenem nu bis zur Gezeitengrenze r_J (volle Verstaerkung ab r*)
  E2  Einzeltrichter mit eigenem nu bis r_e (Gleichstand eigene/aeussere Neigung, ~ r*), dann Elternfeld
  E0  Einzeltrichter Newtonsch; nur das kollektive, ueber r_J geglaettete Feld wird verstaerkt
  L   lokale Fassung von E1 ueber einen Gezeitenschaetzer (Abschnitt 2)

Abschnitte
  0  Grundgroessen, Sonnenumgebung, Hierarchie-Tabelle
  1  Konsistenz: actio=reactio (naiv vs. Energiefunktional), l-Kraefte, WEP-Abschaetzung
  2  Lokale Variante L (Gezeitenschaetzer): Identitaeten, Testteilchen-Gesetz
  3  Quadrupol Q2 (Cassini) je Variante, Galaktische Gezeiten, Periheldrehungen, Zusatzbeschleunigung
  4  Weite Doppelsterne
  5  Aeussere Oortsche Wolke: scheinbare Bahnenergie, Gezeiten-Einschuss
  6  Zwerggalaxien im aeusseren Feld
  7  Verträglichkeits-Uebersicht
Eintraege mit [G] sind Literatur-/Beobachtungswerte aus dem Gedaechtnis (nicht verifiziert).
"""
import math, sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
np.seterr(all="ignore")

# ------------------------------------------------------------------ Konstanten
G    = 6.67430e-11
c    = 299792458.0
AU   = 1.495978707e11
PC   = 3.0856775814913673e16
KPC  = 1e3 * PC
MPC  = 1e6 * PC
YR   = 3.15576e7
MAS  = 180 / math.pi * 3600e3
H0   = 70e3 / MPC
a0H  = c * H0 / (2 * math.pi)       # FTM: a0 = c H0 / 2pi
a0E  = 1.2e-10                      # empirischer RAR-Wert
MSUN = 1.98847e30
GMs  = 1.32712440018e20
KK   = (1e3 / KPC)                  # km/s/kpc -> 1/s

def H(t): print("\n" + "=" * 96 + "\n" + t + "\n" + "=" * 96)

def nu(x):
    x = np.maximum(np.asarray(x, dtype=float), 1e-300)
    return 1.0 + 1.0 / np.expm1(np.sqrt(x))

def xdnudx(x):
    """x dnu/dx = dnu/dln g = -nu(nu-1) s/2"""
    s = math.sqrt(x); n = float(nu(x)); return -n * (n - 1) * s / 2

# Q_s(x) = int_0^{x^2} (nu-1) dy = 4 int_0^s t^3/(e^t-1) dt,   s = sqrt(x),   Q_s(0)=0
_tg = np.linspace(0.0, 150.0, 600001)
_f = np.zeros_like(_tg); _f[1:] = 4 * _tg[1:] ** 3 / np.expm1(_tg[1:])
_cum = np.concatenate([[0.0], np.cumsum(0.5 * (_f[1:] + _f[:-1]) * (_tg[1] - _tg[0]))])
def Qs(x): return np.interp(np.sqrt(np.maximum(x, 0)), _tg, _cum)

def inv_rar(gob, a0):
    lo, hi = 1e-8, 1e4
    for _ in range(200):
        m = math.sqrt(lo * hi)
        if m * float(nu(m)) < gob / a0: lo = m
        else: hi = m
    return math.sqrt(lo * hi) * a0

MU, WMU = np.polynomial.legendre.leggauss(2000)

def g_avg(r, GM, gNe, a0, own=None, sm=None):
    """kugelgemittelte Einwaertsbeschleunigung (exakt ueber Flussidentitaet):
       Newton + eigene Struktur W(r)(nu(g_N)-1)g_N + Fluss der Elternstruktur (nu(|gNe z + g_sm r|)-1)(...)
       own: None (keine eigene Struktur), np.inf (Fenster 1), Zahl (Plummer-Uebergabe bei l=own)
       sm : None (Elternfeld sieht Punktmasse = QUMOND), Zahl (Plummer-geglaettet ueber l=sm)
       gNe: None -> kein Elternfeld"""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    gN = GM / r ** 2
    tot = gN.copy()
    if own is not None:
        W = np.ones_like(r) if np.isinf(own) else 1 - r ** 3 / (r ** 2 + own ** 2) ** 1.5
        tot = tot + W * (nu(gN / a0) - 1) * gN
    if gNe is not None:
        gsm = gN if sm is None else GM * r / (r ** 2 + sm ** 2) ** 1.5
        gs = gsm[:, None]; m = MU[None, :]
        gm = np.sqrt(np.maximum(gs ** 2 + gNe ** 2 + 2 * gs * gNe * m, 0.0))
        tot = tot + 0.5 * np.sum(WMU[None, :] * (nu(gm / a0) - 1) * (gs + gNe * m), axis=1)
    return tot

# ============================================================================ 0
H("0   Grundgroessen, Sonnenumgebung, Hierarchie")
v0, R0 = 233e3, 8.178 * KPC                      # [G] v0 ~ 229-233 km/s, R0 = 8,178 kpc (GRAVITY 2019)
gobs = v0 ** 2 / R0
print(f"   Sonnenumgebung [G]: v0 = {v0/1e3:.0f} km/s, R0 = {R0/KPC:.3f} kpc  ->  g_obs = v0²/R0 = {gobs:.3e} m/s²")
print(f"   (mit v0 = 229 km/s: {229e3**2/R0:.3e} m/s²)")
GNE = {}
for nm, a0 in [("a0=cH0/2π", a0H), ("a0=1,2e-10", a0E)]:
    gNe = inv_rar(gobs, a0); GNE[a0] = gNe
    print(f"   {nm:11s} a0 = {a0:.4e}: g_obs/a0 = {gobs/a0:.3f};  Newtonsches Feld aus RAR-Umkehr g_Ne = {gNe:.3e} "
          f"= {gNe/a0:.3f} a0;  nu_e = {float(nu(gNe/a0)):.4f};  K = dln nu/dln g = {xdnudx(gNe/a0)/float(nu(gNe/a0)):.4f}")
gNeE = GNE[a0E]; gNeH = GNE[a0H]
A_, B_ = 15.3 * KK, -11.9 * KK                   # [G] Oort-Konstanten (Bovy 2017)
rho0 = 0.10 * MSUN / PC ** 3                     # [G] dynamische lokale Dichte ~0,1 M☉/pc³
PhiRR = (A_ - B_) * (-3 * A_ - B_); Phiyy = (A_ - B_) ** 2; Phizz = 4 * math.pi * G * rho0 + 2 * (A_ ** 2 - B_ ** 2)
Tgal = np.diag([PhiRR, Phiyy, Phizz]); trT = np.trace(Tgal); TgalTF = Tgal - trT / 3 * np.eye(3)
TgalF = math.sqrt(np.sum(Tgal ** 2))
tideJ = 4 * A_ * (A_ - B_)
print(f"   Galaktische Gezeiten [G: A=15,3, B=-11,9 km/s/kpc, ρ0=0,10 M☉/pc³]: Φ_RR = {PhiRR:.3e}, Φ_φφ = {Phiyy:.3e}, "
      f"Φ_zz = {Phizz:.3e} s⁻²;  4A(A-B) = {tideJ:.3e} s⁻²;  |T|_F = {TgalF:.3e}")
def rJ_sun(GM): return (GM / tideJ) ** (1 / 3)
rs_E = math.sqrt(GMs / a0E); re_E = math.sqrt(GMs / gNeE); rJs = rJ_sun(GMs)
print(f"   Sonne: r* = {rs_E/AU:.0f} AE,  r_e = √(GM/g_Ne) = {re_E/AU:.0f} AE (= {re_E/rs_E:.3f} r*),  "
      f"r_J = (GM/4A(A-B))^(1/3) = {rJs/AU:.3e} AE = {rJs/PC:.3f} pc")

print("\n   Hierarchie-Tabelle (a0 = 1,2e-10):  r* = √(GM/a0),  r_e: eigenes Newtonfeld = Elternfeld,  r_J: Gezeitengrenze")
print(f"   {'Trichter (in Eltern)':34s} {'GM [m³/s²]':>11s} {'r*':>14s} {'r_e':>14s} {'r_J':>14s} {'r_J/r*':>9s} {'g_Eltern/a0':>11s}")
GMe, GMj = 3.986004418e14, 1.26686534e17
gsun1, gsun5 = GMs / AU ** 2, GMs / (5.2026 * AU) ** 2
rows = [("Erde (Sonne, 1 AE)", GMe, gsun1, AU * (GMe / (3 * GMs)) ** (1 / 3), AU),
        ("Jupiter (Sonne, 5,2 AE)", GMj, gsun5, 5.2026 * AU * (GMj / (3 * GMs)) ** (1 / 3), AU),
        ("Sonne (Milchstrasse)", GMs, gNeE, rJs, AU),
        ("Doppelstern 1,5 M☉ (Milchstrasse)", 1.5 * GMs, gNeE, rJ_sun(1.5 * GMs), AU)]
GMmw, GMm31, Dlg = G * 6e10 * MSUN, G * 1.5e11 * MSUN, 780 * KPC   # [G] M31 ~1,5e11 M☉ baryonisch, D = 780 kpc
rows.append(("Milchstrasse (Lokale Gruppe) [G]", GMmw, GMm31 / Dlg ** 2, Dlg * (GMmw / (3 * GMm31)) ** (1 / 3), KPC))
for nm, GM, gp, rj, un in rows:
    rs = math.sqrt(GM / a0E); re = math.sqrt(GM / gp); u = "AE" if un == AU else "kpc"
    print(f"   {nm:34s} {GM:11.3e} {rs/un:10.4g} {u:>3s} {re/un:10.4g} {u:>3s} {rj/un:10.4g} {u:>3s} {rj/rs:9.3g} {gp/a0E:11.3g}")

# ============================================================================ 1
H("1a  actio = reactio: naive Einzeltrichter-Regel gegen Energiefunktional")
m1, m2, s_ = 1.0, 0.5, 1e4 * AU
g1, g2 = G * m1 * MSUN / s_ ** 2, G * m2 * MSUN / s_ ** 2
F21 = m2 * MSUN * float(nu(g1 / a0E)) * g1      # Kraft auf 2 im verstaerkten Feld von 1 (nu der QUELLE)
F12 = m1 * MSUN * float(nu(g2 / a0E)) * g2
acm = (F21 - F12) / ((m1 + m2) * MSUN)
print(f"   Sonne (1 M☉) + Begleiter (0,5 M☉), s = 10 kAE, naive Regel 'Feld der Quelle mit nu(eigenes Feld der Quelle)':")
print(f"   nu(g1) = {float(nu(g1/a0E)):.4f}, nu(g2) = {float(nu(g2/a0E)):.4f};  F21/F12 = {F21/F12:.4f};  "
      f"Selbstbeschleunigung des Schwerpunkts = {acm:.3e} m/s² = {acm/gobs:.3f} g_obs")
print(f"   ueber 1 Gyr: Δv = {acm*1e9*YR/1e3:.1f} km/s  (Energiefunktional: beide Glieder in EINEM Trichter, F21 = F12 exakt)")

H("1b  l-Kraefte: Abhaengigkeit der Strukturenergie vom Uebergaberadius l (Parameter p)")
def V_own(GM, ell, a0, N=6000):
    r = np.logspace(math.log10(math.sqrt(GM / a0)) - 4, math.log10(ell) + 5, N)
    W = 1 - r ** 3 / (r ** 2 + ell ** 2) ** 1.5
    f = W * Qs(GM / (r ** 2 * a0)) * 4 * math.pi * r ** 2
    return -(a0 ** 2 / (8 * math.pi * G)) * np.trapezoid(f, r)
MU4, WMU4 = np.polynomial.legendre.leggauss(400)
XQ, WQ = np.polynomial.legendre.leggauss(12)
def V_coll(GM, ell, gNe, a0, Nr=1500):
    """-(a0²/8πG)∫[Q_s(y)-Q_s(y_e)] d³x; linearer Term (nu_e-1)(y-y_e) analytisch gemittelt (sonst Ausloeschung),
       Rest ΔQ = ∫_{y_e}^{y}(nu(√y')-nu_e)dy' per Gauss-Legendre"""
    r = np.logspace(math.log10(ell) - 5, math.log10(ell) + 5, Nr)
    gs = (GM * r / (r ** 2 + ell ** 2) ** 1.5)[:, None]; m = MU4[None, :]
    ye = (gNe / a0) ** 2; nue = float(nu(gNe / a0))
    y = (gs ** 2 + gNe ** 2 + 2 * gs * gNe * m) / a0 ** 2
    yp = ye + (y - ye)[..., None] * (1 + XQ) / 2
    dQ = (y - ye) / 2 * np.sum(WQ * (nu(np.sqrt(np.maximum(yp, 0))) - nue), axis=-1)
    dens = (0.5 * np.sum(WMU4[None, :] * dQ, axis=1) + (nue - 1) * gs[:, 0] ** 2 / a0 ** 2) * 4 * math.pi * r ** 2
    return -(a0 ** 2 / (8 * math.pi * G)) * np.trapezoid(dens, r)
Ms = GMs / G
for lab, ell in [("E1 (eigene Struktur bis r_J)", rJs)]:
    d = 1e-3
    dV = (V_own(GMs, ell * (1 + d), a0E) - V_own(GMs, ell * (1 - d), a0E)) / (2 * d * ell)
    Fm = -dV / Ms * (2 / 3) * ell / R0
    ana = (4 / 9) * math.sqrt(GMs * a0E) / R0
    print(f"   {lab}: dV/dl = {dV:.4e} J/m;  l ∝ R^(2/3)  ->  Kraft/Masse = {Fm:.3e} m/s² (tiefe MOND-Naeherung (4/9)√(GMa0)/R0 = {ana:.3e})"
          f";  relativ zu g_obs: {Fm/gobs:.2e}")
d = 1e-2
dVc = (V_coll(GMs, rJs * (1 + d), gNeE, a0E) - V_coll(GMs, rJs * (1 - d), gNeE, a0E)) / (2 * d * rJs)
Fc = -dVc / Ms * (2 / 3) * rJs / R0
Wpl = (3 * math.pi / 32) * G * Ms ** 2 / rJs
print(f"   E0 (kollektives Feld, Sonne ueber r_J geglaettet): dV/dl = {dVc:.4e} J/m -> Kraft/Masse = {Fc:.3e} m/s², relativ {Fc/gobs:.2e}")
print(f"      Vergleich (nu_e-1)·|W_Plummer|/l: {(float(nu(gNeE/a0E))-1)*Wpl/rJs:.3e} J/m")
print(f"   Normierung: Q_s(0) = 0 (kein Vakuumterm je Trichter). Mit der Normierung Q~ aus 06 (Q~(0) = -4π⁴/15 = {-_cum[-1]:.4f})")
print(f"      waere die Vakuumenergie je Zone ∝ l³ und divergiert fuer das Plummer-Fenster (W ~ 3l²/2r²).")

H("1c  WEP-Abschaetzung in E0 (verschiedene Glaettungsradien im galaktischen Strukturpotential)")
rho_ph = 0.01 * MSUN / PC ** 3                    # [G] lokale 'DM'-Dichte ~0,01 M☉/pc³ = Strukturdichte im FTM
grad = 4 * math.pi * G * rho_ph / R0
for nm, GMx in [("Jupiter", GMj), ("Komet 1e13 kg", G * 1e13), ("Stern 0,1 M☉", 0.1 * GMs)]:
    lx = (GMx / tideJ) ** (1 / 3)
    da = 0.5 * abs(rJs ** 2 - lx ** 2) * grad
    print(f"   Sonne (l = {rJs/PC:.2f} pc) gegen {nm:14s} (l = {lx/AU:10.4g} AE): Δa ≈ ½Δ(l²)·∇(4πGρ_ph) = {da:.2e} m/s²,  Δa/g_obs = {da/gobs:.1e}")

# ============================================================================ 2
H("2   Lokale Variante L: g~ = (√10/2)|∂∂Φ_N|²/|∂∂∂Φ_N|")
def tau(n): return 3 * np.outer(n, n) - np.eye(3)
def Aten(n):
    d = np.eye(3)
    return 5 * np.einsum("i,j,k->ijk", n, n, n) - (np.einsum("ij,k->ijk", d, n) + np.einsum("ik,j->ijk", d, n) + np.einsum("jk,i->ijk", d, n))
x0 = np.array([0.37, -0.52, 0.81]); r0n = np.linalg.norm(x0)
def PhiN(x): return -1.0 / np.linalg.norm(x)
def d2(f, x, i, j, h=1e-3):
    e = np.eye(3)
    return (f(x + h * e[i] + h * e[j]) - f(x + h * e[i] - h * e[j]) - f(x - h * e[i] + h * e[j]) + f(x - h * e[i] - h * e[j])) / (4 * h * h)
T_num = np.array([[d2(PhiN, x0, i, j) for j in range(3)] for i in range(3)])
def d3(f, x, i, j, k, h=1e-3):
    e = np.eye(3)
    return (d2(f, x + h * e[k], i, j) - d2(f, x - h * e[k], i, j)) / (2 * h)
C_num = np.array([[[d3(PhiN, x0, i, j, k) for k in range(3)] for j in range(3)] for i in range(3)])
nT, nC = np.sqrt(np.sum(T_num ** 2)), np.sqrt(np.sum(C_num ** 2))
print(f"   Punktmasse GM=1, r = {r0n:.4f}: |T|·r³ = {nT*r0n**3:.5f} (√6 = {math.sqrt(6):.5f}), |C|·r⁴ = {nC*r0n**4:.4f} (√90 = {math.sqrt(90):.4f})")
print(f"   g~ = (√10/2)|T|²/|C| = {math.sqrt(10)/2*nT**2/nC:.6f},  GM/r² = {1/r0n**2:.6f}  (homogenes aeusseres Feld aendert T, C nicht)")
Ff = lambda r: np.exp(-r) / r ** 2
def FT(x): r = np.linalg.norm(x); return Ff(r) * tau(x / r)
def FA(x): r = np.linalg.norm(x); return Ff(r) * Aten(x / r)
s2 = sum(d2(lambda y, i=i, j=j: FT(y)[i, j], x0, i, j) for i in range(3) for j in range(3))
rr = r0n; hh = 1e-4
Fp = (Ff(rr + hh) - Ff(rr - hh)) / (2 * hh); Fpp = (Ff(rr + hh) - 2 * Ff(rr) + Ff(rr - hh)) / hh ** 2
print(f"   Identitaet ∂i∂j[F τ_ij] = 2F''+10F'/r+6F/r²:  numerisch {s2:.6f},  Formel {2*Fpp+10*Fp/rr+6*Ff(rr)/rr**2:.6f}")
s3 = sum(d3(lambda y, i=i, j=j, k=k: FA(y)[i, j, k], x0, i, j, k) for i in range(3) for j in range(3) for k in range(3))
Fh = lambda r: (np.exp(-r) / r ** 2) * (-1 - 2 / r) + 4 * Ff(r) / r     # F' + 4F/r
Fhp = (Fh(rr + hh) - Fh(rr - hh)) / (2 * hh); Fhpp = (Fh(rr + hh) - 2 * Fh(rr) + Fh(rr - hh)) / hh ** 2
print(f"   Identitaet ∂i∂j∂k[F A_ijk] = L[F'+4F/r]:      numerisch {s3:.5f},  Formel {2*Fhpp+10*Fhp/rr+6*Fh(rr)/rr**2:.5f}")
print("   Testteilchen-Gesetz fuer Punktmasse (h = (nu-1)GM):  Φ_s' = -h''/15 + h'/(5r) + 4h/(5r²)   [QUMOND: Φ_s' = h/r²]")
rg = np.logspace(-3, 4, 200001)                  # Einheiten GM = 1, a0 = 1 -> r* = 1, x = 1/r²
hg = nu(1 / rg ** 2) - 1
hp = np.gradient(hg, rg); hpp = np.gradient(hp, rg)
nuhat = 1 + rg ** 2 * (-hpp / 15 + hp / (5 * rg) + 4 * hg / (5 * rg ** 2))
print(f"   {'x=g_N/a0':>9s} {'nu (FTM)':>10s} {'nu^ (L)':>10s} {'nu^/nu-1':>10s}")
for xv in [100, 30, 10, 5, 3, 2, 1, 0.3, 0.1, 0.01, 1e-3, 1e-4]:
    rv = 1 / math.sqrt(xv); i = np.searchsorted(rg, rv)
    print(f"   {xv:9.0e} {float(nu(xv)):10.5f} {nuhat[i]:10.5f} {nuhat[i]/float(nu(xv))-1:+10.4f}")
rL = (math.sqrt(6) * GMs / TgalF) ** (1 / 3)
print(f"   Zone von L fuer die Sonne: |T_eigen|_F = |T_gal|_F bei r_L = (√6 GM/|T_gal|)^(1/3) = {rL/AU:.3e} AE = {rL/PC:.3f} pc (r_J = {rJs/PC:.3f} pc)")
def gtilde_ratio(r, dirv):
    n = dirv / np.linalg.norm(dirv)
    Tn = -GMs * tau(n) / r ** 3 + Tgal; Cn = 3 * GMs * Aten(n) / r ** 4
    return (math.sqrt(10) / 2 * np.sum(Tn ** 2) / math.sqrt(np.sum(Cn ** 2))) / (GMs / r ** 2)
print("   g~/g_eigen fuer die Sonne im galaktischen Gezeitenfeld:  " + ";  ".join(
    f"r = {f_:.2f} r_L: x-Richtung {gtilde_ratio(f_*rL, np.array([1.,0,0])):.4f}, z-Richtung {gtilde_ratio(f_*rL, np.array([0,0,1.])):.4f}"
    for f_ in [0.01, 0.1, 0.3, 1.0]))

# ============================================================================ 3
H("3a  Quadrupol Q2 am Ort der Sonne: Elternstruktur sieht die Sonne ueber l geglaettet (Plummer), λ = l/r_e")
print("   J(β,λ) = ½∫∫(nu(β|ĝ|)-nu(β))[-3(3c²-1)h/ρ² + (9c-15c³)/ρ²]dρdc + (4/15)·β nu'(β)/λ³ ,  h = ρ/(ρ²+λ²)^(3/2) (λ=0: h=1/ρ²)")
print("   Q2 = -(3/2)(g_e/r_e)·J,  r_e = √(GM/g_e)   (λ=0 ist exakt die Rechnung aus 06)")
def J_quad(beta, lam, Nn=240):
    xg, wg = np.polynomial.legendre.leggauss(Nn)
    bps = {-10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 3.5, 6.0, 9.0, 12.0}
    if lam > 0:
        L = math.log(lam); bps |= {L - 1.5, L - 0.7, L - 0.3, L, L + 0.3, L + 0.7, L + 1.5}
    bps = sorted(b for b in bps if -10.0 <= b <= 12.0)
    tseg = [(0.0, 0.5 * math.pi), (0.5 * math.pi, 0.9 * math.pi), (0.9 * math.pi, 0.99 * math.pi), (0.99 * math.pi, math.pi)]
    nue = float(nu(beta)); tot = 0.0
    for u0, u1 in zip(bps[:-1], bps[1:]):
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0); wu = 0.5 * (u1 - u0) * wg; rho = np.exp(u)
        for t0, t1 in tseg:
            th = 0.5 * (t1 - t0) * xg + 0.5 * (t1 + t0); wt = 0.5 * (t1 - t0) * wg
            R_, TH = np.meshgrid(rho, th, indexing="ij"); cc = np.cos(TH)
            h = 1 / R_ ** 2 if lam == 0 else R_ / (R_ ** 2 + lam ** 2) ** 1.5
            gm = np.sqrt(np.maximum(h ** 2 + 2 * cc * h + 1, 0.0))
            kern = -3 * (3 * cc ** 2 - 1) * h / R_ ** 2 + (9 * cc - 15 * cc ** 3) / R_ ** 2
            # (nu-1) statt (nu-nu_e): der Anteil (nu_e-1)ĝ traegt weder zum Volumen- noch zum Oberflaechenterm bei,
            # vermeidet aber Ausloeschungsrauschen bei ρ -> 0 (Punktmasse)
            F = 0.5 * (nu(beta * gm) - 1.0) * kern * R_ * np.sin(TH)
            tot += np.einsum("i,j,ij->", wu, wt, F)
    surf = 0.0 if lam == 0 else (4 / 15) * xdnudx(beta) / lam ** 3
    return tot, surf
def J_grid(beta, lam, Nu=3000, Nt=3000):
    """unabhaengig: Strukturdichte S = ∇·[(nu-nu_e)ĝ] auf Gitter, dann Hauptwert -(1/4π)∫S ∂z²(1/r)"""
    u = np.linspace(-8.0, 9.0, Nu); th = np.linspace(1e-6, math.pi - 1e-6, Nt)
    R_, TH = np.meshgrid(np.exp(u), th, indexing="ij"); cc, sn = np.cos(TH), np.sin(TH)
    h = R_ / (R_ ** 2 + lam ** 2) ** 1.5
    gr, gt = h + cc, -sn
    gm = np.sqrt(gr ** 2 + gt ** 2)
    dn = nu(beta * gm) - float(nu(beta))
    Vr, Vt = dn * gr, dn * gt
    S = np.gradient(R_ ** 2 * Vr, u, axis=0) / R_ / R_ ** 2 + np.gradient(sn * Vt, th, axis=1) / (R_ * sn)
    K = (3 * cc ** 2 - 1) / R_ ** 3
    return -(1 / (4 * math.pi)) * 2 * math.pi * np.trapezoid(np.trapezoid(S * K * R_ ** 3 * sn, th, axis=1), u)
v1, s1 = J_quad(1.0, 0.0)
print(f"   Probe gegen 06 (β = 1, Punktmasse): J = {v1+s1:.6f}   (06: -0,169873)")
for lam in [0.5, 1.0, 3.0]:
    v, s = J_quad(1.0, lam); v2, s2_ = J_quad(1.0, lam, Nn=400)
    print(f"   Probe λ = {lam}: Quadratur J = {v+s:.6f} (Volumen {v:.6f}, Oberflaeche {s:.6f}); N=400: {v2+s2_:.6f};  Gitter-Hauptwert: {J_grid(1.0, lam):.6f}")
lamJ = rJs / re_E
print(f"\n   Sonne mit g_e = g_Ne aus RAR-Umkehr; λ_J = r_J/r_e = {lamJ:.1f};  Cassini: Q2 = (3 ± 3)·10⁻²⁷ s⁻² (Hees et al. 2014)")
print(f"   {'λ = l/r_e':>10s} {'l [AE]':>9s} {'J':>10s} {'Q2 [s⁻²] a0=1,2e-10':>21s} {'(Q2-3e-27)/3e-27':>17s} {'Q2 [s⁻²] a0=cH0/2π':>20s}")
Q2tab = {}
reH = math.sqrt(GMs / gNeH)
for lam in [0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, lamJ]:
    v, s = J_quad(gNeE / a0E, lam); J = v + s; Q2 = -1.5 * (gNeE / re_E) * J
    vH, sH = J_quad(gNeH / a0H, lam); Q2H = -1.5 * (gNeH / reH) * (vH + sH)
    Q2tab[lam] = Q2
    print(f"   {lam:10.2f} {lam*re_E/AU:9.0f} {J:10.5f} {Q2:21.3e} {(Q2-3e-27)/3e-27:17.1f} {Q2H:20.3e}")
lams = np.array(sorted(Q2tab)); qs = np.array([Q2tab[l] for l in lams])
def lam_for(Qlim):
    for i in range(len(lams) - 1):
        if qs[i] >= Qlim > qs[i + 1]:
            return math.exp(np.interp(math.log(Qlim), [math.log(qs[i + 1]), math.log(qs[i])], [math.log(max(lams[i + 1], 1e-9)), math.log(max(lams[i], 1e-9))]))
    return float("nan")
for nm, Ql in [("1σ (6e-27)", 6e-27), ("2σ (9e-27)", 9e-27), ("3σ (1,2e-26)", 1.2e-26)]:
    lm = lam_for(Ql); print(f"   Cassini-Grenze {nm}: λ ≥ {lm:.2f}  ->  l ≥ {lm*re_E/AU:.0f} AE = {lm*re_E/rs_E:.2f} r*")
print("   Zum Vergleich 06: g_e = 1,2e-10 (β=1): Q2 = " + f"{-1.5*(1.2e-10/math.sqrt(GMs/1.2e-10))*sum(J_quad(1.0,0.0)):.3e} s⁻²")
Q2tideH = -1.5 * TgalTF[0, 0]
print(f"   Galaktische Gezeit (effektiv, in jeder Variante vorhanden): spurfreier Tensor diag = ({TgalTF[0,0]:.2e}, {TgalTF[1,1]:.2e}, {TgalTF[2,2]:.2e}) s⁻²;")
print(f"      Projektion auf Hees-Form (Richtung GZ): Q2_Gezeit = -(3/2)T^TF_xx = {Q2tideH:.2e} s⁻²;  g_obs/R0 = {gobs/R0:.2e};  GM/r_J³ = {GMs/rJs**3:.2e}")

H("3b  Periheldrehung der Planeten durch Q2 (Gauss-Gleichungen, Bahnmittel)")
eps = math.radians(23.4393); ra, de = math.radians(266.405), math.radians(-28.936)
xq = np.array([math.cos(de) * math.cos(ra), math.cos(de) * math.sin(ra), math.sin(de)])
egc = np.array([xq[0], xq[1] * math.cos(eps) + xq[2] * math.sin(eps), -xq[1] * math.sin(eps) + xq[2] * math.cos(eps)])
lam_e, bet_e = math.degrees(math.atan2(egc[1], egc[0])) % 360, math.degrees(math.asin(egc[2]))
print(f"   Richtung Galaktisches Zentrum in Ekliptik-Koordinaten: λ = {lam_e:.2f}°, β = {bet_e:.2f}°")
def dvarpi(aAU, e, inc, Om, om, acc, GM=GMs, N=40000):
    a = aAU * AU; n = math.sqrt(GM / a ** 3); p = a * (1 - e * e)
    M = np.linspace(0, 2 * math.pi, N, endpoint=False); E = M.copy()
    for _ in range(50): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    f = 2 * np.arctan2(math.sqrt(1 + e) * np.sin(E / 2), math.sqrt(1 - e) * np.cos(E / 2)); r = a * (1 - e * np.cos(E))
    i, O, w = map(math.radians, (inc, Om, om)); u = w + f
    rh = np.stack([math.cos(O) * np.cos(u) - math.sin(O) * np.sin(u) * math.cos(i),
                   math.sin(O) * np.cos(u) + math.cos(O) * np.sin(u) * math.cos(i), np.sin(u) * math.sin(i)], 1)
    wh = np.array([math.sin(O) * math.sin(i), -math.cos(O) * math.sin(i), math.cos(i)])
    sh = np.cross(wh[None, :], rh)
    av = acc(rh * r[:, None])
    Rc, Sc, Wc = np.sum(av * rh, 1), np.sum(av * sh, 1), av @ wh
    d = math.sqrt(1 - e * e) / (n * a * e) * (-Rc * np.cos(f) + Sc * (1 + r / p) * np.sin(f)) \
        + math.tan(i / 2) * r * np.sin(u) * Wc / (n * a * a * math.sqrt(1 - e * e))
    return float(np.mean(d))
Atest = 1e-12
t_num = dvarpi(1.52368, 0.0934, 1.85, 49.56, 286.5, lambda X: Atest * X / np.linalg.norm(X, axis=1)[:, None])
a_ = 1.52368 * AU; t_ana = Atest * math.sqrt(1 - 0.0934 ** 2) / (math.sqrt(GMs / a_ ** 3) * a_)
print(f"   Probe (konstante Radialbeschleunigung 1e-12 an Mars): numerisch {t_num:.6e}, analytisch A√(1-e²)/(na) = {t_ana:.6e} rad/s")
planets = [("Merkur", 0.387098, 0.205630, 7.0050, 48.331, 29.124, -0.020, 0.030),
           ("Venus", 0.723332, 0.006773, 3.3946, 76.680, 54.884, 0.026, 0.016),
           ("Erde", 1.000000, 0.016709, 0.0, 0.0, 102.937, 0.0019, 0.0019),
           ("Mars", 1.523679, 0.093405, 1.8497, 49.558, 286.502, -0.00020, 0.00037),
           ("Jupiter", 5.2026, 0.048498, 1.3033, 100.464, 273.867, 0.587, 0.283),
           ("Saturn", 9.5549, 0.055546, 2.4845, 113.665, 339.392, -0.0032, 0.0047),
           ("Uranus", 19.2184, 0.046381, 0.7732, 74.006, 96.999, -3.89, 3.90),
           ("Neptun", 30.1104, 0.009456, 1.7700, 131.784, 276.336, -4.44, 5.40)]   # Bahnelemente J2000, EPM2011 [G]
def accQ(Q2): return lambda X: Q2 * (np.outer(X @ egc, egc) - X / 3)
Q2G, Q2E2, Q2E0 = Q2tab[0.0], Q2tab[1.0], Q2tab[lamJ]
print(f"   Q2-Werte: G = {Q2G:.2e}, E2 (λ=1) = {Q2E2:.2e}, E1/E0 (λ_J) = {Q2E0:.2e} (+ Gezeit {Q2tideH:.1e}) s⁻²")
print(f"   {'Planet':8s} {'Δϖ̇ je 1e-26 s⁻²':>17s} {'G':>11s} {'E2':>11s} {'E1/E0+Gez.':>12s} {'EPM2011 [G]':>17s} {'G/σ':>7s} {'E2/σ':>7s} {'Q2-Grenze 1σ':>13s}  [mas/Jahr]")
for nm, a, e, i, O, w, obs, err in planets:
    unit = dvarpi(a, e, i, O, w, accQ(1e-26)) * YR * MAS
    vG, vE2, vE0 = unit * Q2G / 1e-26, unit * Q2E2 / 1e-26, unit * (Q2E0 + Q2tideH) / 1e-26
    print(f"   {nm:8s} {unit:17.4e} {vG:11.3e} {vE2:11.3e} {vE0:12.2e} {obs:8.4g} ± {err:<7.2g} {abs(vG)/err:7.1f} {abs(vE2)/err:7.1f} {err/abs(unit)*1e-26:13.1e}")
print("\n   Zusatzbeschleunigungen: Quadrupol |Q2|·a und eigener Strukturterm (nu-1)g_N (E1: nu vom eigenen Feld)")
print(f"   {'Planet':8s} {'|Q2_G|·a':>11s} {'|Q2_E2|·a':>11s} {'|Q2_E0|·a':>11s} {'log10[(nu-1)g_N]':>17s}  [m/s²]")
for nm, a, *_ in planets:
    aa = a * AU; gN = GMs / aa ** 2; s = math.sqrt(gN / a0E)
    print(f"   {nm:8s} {abs(Q2G)*aa:11.2e} {abs(Q2E2)*aa:11.2e} {abs(Q2E0+Q2tideH)*aa:11.2e} {(-s/math.log(10)+math.log10(gN)):17.1f}")

# ============================================================================ 4
H("4   Weite Doppelsterne: γ = <g>/g_N (kugelgemittelt), v/v_N = √γ")
Mb = 1.5; GMb = Mb * GMs
reb, rJb = math.sqrt(GMb / gNeE), rJ_sun(GMb)
print(f"   M_tot = {Mb} M☉, a0 = 1,2e-10: r* = {math.sqrt(GMb/a0E)/AU:.0f} AE, r_e = {reb/AU:.0f} AE, r_J = {rJb/AU:.3e} AE, g_Ne = {gNeE:.3e}")
VARS = [("Newton", None, None, None),
        ("G (g_Ne RAR)", gNeE, None, None),
        ("G (1,5 a0)", 1.5 * a0E, None, None),
        ("G (2 a0)", 2.0 * a0E, None, None),
        ("E2 (l=r_e)", gNeE, reb, reb),
        ("E1 (l=r_J)", gNeE, rJb, rJb),
        ("E0 (l=r_J)", gNeE, None, rJb)]
svals = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30]) * 1e3 * AU
print(f"   {'s [kAE]':>8s} {'g_N [m/s²]':>11s} " + " ".join(f"{v[0]:>13s}" for v in VARS))
WB = {}
for nm, gne, own, sm in VARS:
    WB[nm] = g_avg(svals, GMb, gne, a0E, own=own, sm=sm) / (GMb / svals ** 2) if gne is not None else np.ones_like(svals)
for k, s in enumerate(svals):
    print(f"   {s/AU/1e3:8.0f} {GMb/s**2:11.3e} " + " ".join(f"{WB[v[0]][k]:13.4f}" for v in VARS))
print("   v/v_N = √γ:")
for k, s in enumerate(svals):
    print(f"   {s/AU/1e3:8.0f} {'':11s} " + " ".join(f"{math.sqrt(WB[v[0]][k]):13.4f}" for v in VARS))
print("   a0 = cH0/2π (g_Ne entsprechend):")
rebH = math.sqrt(GMb / gNeH)
for nm, gne, own, sm in [("G", gNeH, None, None), ("E2", gNeH, rebH, rebH), ("E1", gNeH, rJb, rJb), ("E0", gNeH, None, rJb)]:
    gg = g_avg(svals, GMb, gne, a0H, own=own, sm=sm) / (GMb / svals ** 2)
    print(f"   {nm:3s} γ: " + " ".join(f"{x:.3f}" for x in gg))
be = gNeE / a0E; ne = float(nu(be)); Ke = xdnudx(be) / ne
print(f"   Lineares EFE-Regime (G, s ≫ r_e): γ(θ) = nu_e(1 + (K/2)sin²θ), nu_e = {ne:.4f}, K = {Ke:.4f}: "
      f"γ(θ=0) = {ne:.3f}, γ(90°) = {ne*(1+Ke/2):.3f}, Mittel nu_e(1+K/3) = {ne*(1+Ke/3):.3f}")
print("   Tief-MOND-Zweikoerperfaktor a_rel/√(G M a0)/s = (2/3)[1-p^1,5-(1-p)^1,5]/(p(1-p)) [Milgrom, G]:  " +
      ", ".join(f"p={p}: {(2/3)*(1-p**1.5-(1-p)**1.5)/(p*(1-p)):.3f}" for p in [0.5, 0.3, 0.1, 0.01]))
print(f"   Gezeitenbeitrag |T|_F·s/g_N bei 30 kAE: {TgalF*30e3*AU/(GMb/(30e3*AU)**2):.2e}")
print("   Daten [G]: Chae 2023 γ ≈ 1,43 ± 0,06 fuer g_N ≲ 1e-10; Hernandez 2023 aehnlich; Banik et al. 2024 γ ≈ 1 (MOND mit EFE stark verworfen).")
print("\n   Familie E(λ): Uebergabe bei l = λ r_e (eigene Struktur innen, Elternstruktur mit geglaetteter Masse aussen)")
print(f"   {'λ':>6s} {'Q2 Sonne':>10s} {'Q2/σ':>6s} {'γ(5 kAE)':>9s} {'γ(10 kAE)':>10s} {'γ(20 kAE)':>10s} {'γ(30 kAE)':>10s}")
for lam in [0.0, 1.0, 1.5, 2.0, 3.0, 5.0, lamJ]:
    sv = np.array([5e3, 1e4, 2e4, 3e4]) * AU
    if lam == 0: gg = g_avg(sv, GMb, gNeE, a0E)
    else: gg = g_avg(sv, GMb, gNeE, a0E, own=lam * reb, sm=lam * reb)
    gg = gg / (GMb / sv ** 2)
    print(f"   {lam:6.2f} {Q2tab[lam]:10.2e} {(Q2tab[lam]-3e-27)/3e-27:6.1f} {gg[0]:9.3f} {gg[1]:10.3f} {gg[2]:10.3f} {gg[3]:10.3f}")

# ============================================================================ 5
H("5   Aeussere Oortsche Wolke (Sonne, a0 = 1,2e-10)")
rg5 = np.logspace(math.log10(100 * AU), math.log10(5e6 * AU), 5000)
OVARS = [("Newton", None, None, None), ("G", gNeE, None, None), ("E2", gNeE, re_E, re_E),
         ("E1", gNeE, rJs, rJs), ("E0", gNeE, None, rJs)]
PHI = {}
for nm, gne, own, sm in OVARS:
    ga = GMs / rg5 ** 2 if gne is None else g_avg(rg5, GMs, gne, a0E, own=own, sm=sm)
    extra = ga - GMs / rg5 ** 2
    dPhi = np.concatenate([[0.0], np.cumsum(0.5 * (extra[1:] + extra[:-1]) * np.diff(rg5))])
    PHI[nm] = (ga, dPhi)
print("   γ = <g>/g_N (kugelgemittelt); Verhaeltnis galaktische Gezeit |T|·r zu Eigenbeschleunigung (Newton bzw. E1)")
print(f"   {'r [AE]':>8s} " + " ".join(f"{'γ '+v[0]:>10s}" for v in OVARS) + f" {'T r/g (N)':>10s} {'T r/g (E1)':>11s}")
for rA in [1e4, 2e4, 5e4, 1e5, 2e5]:
    r = rA * AU; gam = {v[0]: float(np.interp(r, rg5, PHI[v[0]][0])) / (GMs / r ** 2) for v in OVARS}
    print(f"   {rA:8.0f} " + " ".join(f"{gam[v[0]]:10.3f}" for v in OVARS) + f" {TgalF*r/(GMs/r**2):10.2e} {TgalF*r/(gam['E1']*GMs/r**2):11.2e}")
print("\n   Scheinbare 'urspruengliche' Bahnenergie: Bahnbestimmung bei ~250 AE ist Newtonsch (Strukturpotential dort konstant).")
print("   E_N = -GM/Q + ΔΦ_s(Q),  1/a_app = -2E_N/GM  fuer einen radialen Kometen mit wahrem Aphel Q")
print(f"   {'Q [AE]':>8s} " + " ".join(f"{'1/a '+v[0]:>13s}" for v in OVARS) + "   [1e-5 AE⁻¹]")
for QA in [5e3, 1e4, 2e4, 3e4, 5e4, 1e5, 2e5]:
    Q = QA * AU
    vals = [(-2 * (-GMs / Q + float(np.interp(Q, rg5, PHI[v[0]][1]))) / GMs) * AU / 1e-5 for v in OVARS]
    print(f"   {QA:8.0f} " + " ".join(f"{x:13.2f}" for x in vals))
print("   Beobachtete Oort-Spitze [G]: 1/a_orig ≈ 2…5·10⁻⁵ AE⁻¹ (gebunden); stark hyperbolische Ursprungsbahnen fehlen.")
for nm in ["G", "E2", "E1"]:
    ga, dP = PHI[nm]
    f0 = -GMs / rg5 + dP
    i0 = np.argmax(f0 > 0)
    print(f"   {nm}: Aphel, ab dem ein Oortkomet hyperbolisch erscheint (E_N = 0): Q0 = {rg5[i0]/AU:.0f} AE;  "
          f"ΔΦ_s(∞) - ΔΦ_s(0) ≈ {dP[-1]:.3e} m²/s²  (vgl. E_N der Spitze bei 1/a = 3e-5: {-GMs*3e-5/AU/2:.3e})")

nuz2 = Phizz
xt, wt = np.polynomial.legendre.leggauss(3000); xt = 0.5 * (xt + 1); wt = 0.5 * wt
def radial_orbit(nm, Q):
    ga, dP = PHI[nm]
    Phi = lambda r: -GMs / r + np.interp(r, rg5, dP, left=0.0)
    r = Q * (1 - xt ** 2); dr = 2 * Q * xt
    v = np.sqrt(2 * np.maximum(Phi(Q) - Phi(r), 1e-300))
    P = 2 * np.sum(wt * dr / v)
    r2 = 2 * np.sum(wt * r ** 2 * dr / v) / P
    inv_r_out = 2 * np.sum(wt * np.where(r > re_E, 1 / r, 0.0) * dr / v)      # ∫ dt/r fuer r>r_e (EFE-Drehmoment)
    return P, r2, inv_r_out
print("\n   Drehmoment je Masse: EFE-Anisotropie (G, E2; lineares Regime) nu_e|K|GM/(2r) gegen galaktische Gezeit ½|T|_F r²")
for rA in [6442, 1e4, 3e4, 1e5]:
    r = rA * AU
    print(f"   r = {rA:7.0f} AE: EFE {ne*abs(Ke)*GMs/(2*r):.3e} m²/s², Gezeit {0.5*TgalF*r**2:.3e} m²/s², Verhaeltnis {ne*abs(Ke)*GMs/(2*r)/(0.5*TgalF*r**2):.2e}")
sq2GM = math.sqrt(2 * GMs); sAU = math.sqrt(AU)
def dsq(nm, Q):
    P, r2, ir = radial_orbit(nm, Q)
    tide = 0.5 * nuz2 * r2 * P / sq2GM / sAU
    efe = (ne * abs(Ke) * GMs / 2) * ir / sq2GM / sAU if nm in ("G", "E2") else 0.0
    return P, tide, efe
print(f"\n   Gezeiten-Einschuss: Δ√q je Umlauf ≈ ½ν_z²<r²>P/√(2GM) (+ EFE-Drehmoment nu_e|K|GM/2r fuer r>r_e bei G, E2);"
      f" Sprung ueber die Jupiter-Saturn-Barriere (15 → 5 AE) verlangt Δ√q ≥ √15-√5 = {math.sqrt(15)-math.sqrt(5):.2f} √AE")
print(f"   {'Q [AE]':>8s} " + " ".join(f"{'P '+v[0]+' [Myr]':>14s}" for v in OVARS))
for QA in [1e4, 3e4, 6e4, 1e5]:
    print(f"   {QA:8.0f} " + " ".join(f"{radial_orbit(v[0], QA*AU)[0]/YR/1e6:14.3f}" for v in OVARS))
print(f"   {'Variante':8s} {'Q_krit [AE]':>12s} {'Δ√q Gezeit':>11s} {'Δ√q EFE':>9s} {'1/a_app(Q_krit) [1e-5 AE⁻¹]':>28s}")
need = math.sqrt(15) - math.sqrt(5)
for nm, *_ in OVARS:
    lo, hi = 2e3 * AU, 3e6 * AU
    for _ in range(60):
        mid = math.sqrt(lo * hi); _, t_, e_ = dsq(nm, mid)
        if t_ + e_ < need: lo = mid
        else: hi = mid
    Qc = math.sqrt(lo * hi); _, t_, e_ = dsq(nm, Qc)
    ia = (-2 * (-GMs / Qc + float(np.interp(Qc, rg5, PHI[nm][1]))) / GMs) * AU / 1e-5
    print(f"   {nm:8s} {Qc/AU:12.0f} {t_:11.3f} {e_:9.3f} {ia:28.2f}")
print(f"   Probe Newton: Q_krit/2 = a_krit (Literatur [G]: ~2–3·10⁴ AE)")

# ============================================================================ 6
H("6   Zwerggalaxien im aeusseren Feld: σ_los = √(r_½ g(r_½)/3) (Wolf-Schaetzer), M(<r_½) = M/2, a0 = 1,2e-10")
dwarfs = [("Crater II (MW)", 3.2e5, 1.066, 6e10, 117.5, 2.7, 0.3, 0.3),
          ("And XIX (M31) [unsicher]", 8.0e5, 1.6, 1.5e11, 115.0, 4.7, 1.6, 1.4),
          ("NGC1052-DF2", 2.0e8, 2.2, 1.0e11, 80.0, 8.5, 2.3, 3.1),
          ("NGC1052-DF4", 1.5e8, 1.6, 1.0e11, 165.0, 4.2, 4.4, 2.2)]
print("   Eingaben [G]: M = M*(Υ=2), R_e, Wirt-Baryonenmasse, Abstand (3D ≈ projiziert), σ_obs.")
print(f"   {'System':26s} {'r_½[kpc]':>8s} {'g_Ne/a0':>8s} {'r_e[kpc]':>8s} {'l_J[kpc]':>8s} {'Newton':>7s} {'isol.':>7s} {'G=E0':>7s} {'E2':>7s} {'E1':>7s} {'σ_obs':>13s}  [km/s]")
for nm, M, Re, Mh, D, so, ep, em in dwarfs:
    GM = G * M * MSUN; rh = 4 / 3 * Re * KPC; GMh = GM / 2
    gNe = G * Mh * MSUN / (D * KPC) ** 2; ge = float(nu(gNe / a0E)) * gNe
    lJ = (GM * D * KPC / (2 * ge)) ** (1 / 3); re_d = math.sqrt(GM / gNe)
    def sig(g): return math.sqrt(rh * g / 3) / 1e3
    gN = GMh / rh ** 2
    sN = sig(gN); sI = sig(float(nu(gN / a0E)) * gN)
    sG = sig(float(g_avg(rh, GMh, gNe, a0E)[0]))
    sE2 = sig(float(g_avg(rh, GMh, gNe, a0E, own=re_d, sm=re_d)[0]))
    sE1 = sig(float(g_avg(rh, GMh, gNe, a0E, own=lJ, sm=lJ)[0]))
    print(f"   {nm:26s} {rh/KPC:8.2f} {gNe/a0E:8.4f} {re_d/KPC:8.2f} {lJ/KPC:8.2f} {sN:7.2f} {sI:7.2f} {sG:7.2f} {sE2:7.2f} {sE1:7.2f} {so:6.1f} +{ep}/-{em}")
print("   E1 mit l_J aus der eigenen EFFEKTIVEN Neigung (√(GM a0)/r = 2 g_e r/D):")
for nm, M, Re, Mh, D, so, ep, em in dwarfs:
    GM = G * M * MSUN; rh = 4 / 3 * Re * KPC; GMh = GM / 2
    gNe = G * Mh * MSUN / (D * KPC) ** 2; ge = float(nu(gNe / a0E)) * gNe
    lJe = math.sqrt(math.sqrt(GM * a0E) * D * KPC / (2 * ge))
    sE1e = math.sqrt(rh * float(g_avg(rh, GMh, gNe, a0E, own=lJe, sm=lJe)[0]) / 3) / 1e3
    print(f"   {nm:26s} l_J,eff = {lJe/KPC:7.2f} kpc -> σ(E1) = {sE1e:.2f} km/s")
print("   Literatur [G]: Crater II σ = 2,7 ± 0,3 (Caldwell 2017), MOND-EFE-Vorhersage 2,1 (McGaugh 2016);")
print("   DF2 σ = 8,5 +2,3/-3,1 (Danieli 2019), MOND isoliert ~20, mit EFE ~13–14 (Famaey, McGaugh & Milgrom 2018); DF4 σ = 4,2 +4,4/-2,2 (van Dokkum 2019).")

# ============================================================================ 7
H("7   Uebersicht (a0 = 1,2e-10; Sonne mit g_Ne aus RAR-Umkehr; Doppelstern 1,5 M☉)")
def ia_at(nm, QA): return (-2 * (-GMs / (QA * AU) + float(np.interp(QA * AU, rg5, PHI[nm][1]))) / GMs) * AU / 1e-5
print(f"   {'Variante':8s} {'Q2 [s⁻²]':>10s} {'Q2/σ':>6s} {'γ(10kAE)':>9s} {'γ(20kAE)':>9s} {'1/a_app(Q=5e4 AE)':>18s}")
for nm, key, lam in [("G", "G (g_Ne RAR)", 0.0), ("E2", "E2 (l=r_e)", 1.0), ("E1", "E1 (l=r_J)", lamJ), ("E0", "E0 (l=r_J)", lamJ)]:
    q = Q2tab[lam] + (Q2tideH if lam == lamJ else 0.0)
    print(f"   {nm:8s} {q:10.2e} {(q-3e-27)/3e-27:6.1f} {WB[key][5]:9.3f} {WB[key][7]:9.3f} {ia_at(nm, 5e4):18.2f}")
print(f"   Newton   {'-':>10s} {'-':>6s} {1.0:9.3f} {1.0:9.3f} {ia_at('Newton', 5e4):18.2f}")
