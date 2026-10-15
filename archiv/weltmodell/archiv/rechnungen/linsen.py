# -*- coding: utf-8 -*-
"""Gravitationslinsen im FTM-Operator (Punkt II.2 des Autors).

Operator (Notiz 01):  ∇²Φ = ∇·[ν(|∇Φ_N|/a0) ∇Φ_N],  ν = 1/(1-q),  q = exp(-sqrt(a_N/a0))
Faltung: Δs = GM/(3c²),  Schwelle r* = sqrt(GM/a0).

Teile:
 1  Linsengleichung, Optionen A / B_halb / B_strikt, Punktmasse in Δs und r*
 2  Galaxie-Galaxie-Linsen: ΔΣ(R) 10 kpc - 1 Mpc, Vergleich mit der Linsen-RAR (Brouwer et al. 2021)
 3  Massereicher Haufen: Einstein-Radius, M_2D(<500 kpc), Rest gegenüber NFW-Referenz
 4  Bullet Cluster: 3D-QUMOND auf einem Gitter (FFT-Poisson, isolierte Randbedingung), κ(x,y)
 5  Gedächtnis-Variante: Operatorfeld läuft der Zustandsgeschichte nach (Nachlaufzeit 1/λ)
Nur numpy. Ausgabe: linsen.txt (python linsen.py > linsen.txt)
"""
import math, sys, time
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------------------------------------------------------- Konstanten
G = 6.6743e-11; c = 299792458.0
MS = 1.98847e30; PC = 3.085677581e16; KPC = 1e3 * PC; MPC = 1e6 * PC
YR = 3.15576e7; MYR = 1e6 * YR; GYR = 1e9 * YR
ARCSEC = math.pi / (180 * 3600)
A0 = 1.2e-10
H0 = 70e3 / MPC; OM, OL = 0.3, 0.7
A0_H0 = c * H0 / (2 * math.pi)
SIG_PC2 = MS / PC**2          # 1 M☉/pc² in kg/m²

def E(z): return math.sqrt(OM * (1 + z)**3 + OL)

def nu_ftm(y):
    """ν = 1/(1-q), q = exp(-sqrt(y)) (Rückkopplungssumme, Notiz 01)."""
    s = np.sqrt(np.maximum(np.asarray(y, float), 1e-300))
    return 1.0 / (-np.expm1(-np.maximum(s, 1e-12)))

def nu_simple(y):
    y = np.maximum(np.asarray(y, float), 1e-300); return 0.5 + np.sqrt(0.25 + 1 / y)

def nu_standard(y):
    y = np.maximum(np.asarray(y, float), 1e-300); return np.sqrt(0.5 + np.sqrt(0.25 + 1 / y**2))

XG, WG = np.polynomial.legendre.leggauss(3000)

def Mproj(R, gfun, zmax=None):
    """Projizierte Linsenmasse M_2D(<R) = (R²/G) ∫_0^{θmax} g(R/cosθ)/cosθ dθ  (z = R tanθ).
    Folgt aus α(R) = (2/c²)∫ g_⊥ dz = 4 G M_2D(<R)/(c² R) für axialsymmetrische Linsen."""
    R = np.atleast_1d(np.asarray(R, float))[:, None]
    tmax = math.pi / 2 if zmax is None else np.arctan(zmax / R)
    th = 0.5 * tmax * (XG[None, :] + 1); w = 0.5 * tmax * WG[None, :]
    ct = np.cos(th)
    return (R[:, 0]**2 / G) * np.sum(w * gfun(R / ct) / ct, axis=1)

def Sigma_from_Mproj(R, Mfun):
    h = 1e-3 * R
    dM = (Mfun(R + h) - Mfun(R - h)) / (2 * h)
    return dM / (2 * math.pi * R)

def DeltaSigma(R, Mfun):
    return Mfun(R) / (math.pi * R**2) - Sigma_from_Mproj(R, Mfun)

def DC(z, n=4000):
    zz = np.linspace(0, z, n); return c / H0 * np.trapezoid(1 / np.sqrt(OM * (1 + zz)**3 + OL), zz)

def Sigma_crit(zl, zs):
    Dl = DC(zl) / (1 + zl); Ds = DC(zs) / (1 + zs); Dls = (DC(zs) - DC(zl)) / (1 + zs)
    return c**2 / (4 * math.pi * G) * Ds / (Dl * Dls), Dl, Ds, Dls

def fmt(x, p=3): return f"{x:.{p}e}"

print("=" * 100)
print("FTM – Gravitationslinsen quantitativ (linsen.py)")
print(f"G = {G}, c = {c}, M☉ = {MS}, kpc = {KPC:.6e} m, a0 = {A0}, c·H0/2π (H0=70) = {A0_H0:.4e} m/s²")
print("Kosmologie nur für Winkeldurchmesser-Abstände: flach, H0 = 70 km/s/Mpc, Ωm = 0.3 (Annahme)")
print("=" * 100)

# =============================================================================== 1
print("\n1) LINSENGLEICHUNG")
print("   Metrik schwaches Feld: ds² = -(1+2Φ/c²)c²dt² + (1-2Ψ/c²)dx²")
print("   Ablenkung: α(b) = (1/c²)∫∇⊥(Φ+Ψ)dl   [mit Faktor 1/c², nicht 2/c²; Φ=Ψ ergibt 4GM/(c²b)]")
Msun_ds = G * MS / (3 * c**2); Rsun = 6.957e8
print(f"   Kontrolle Sonne: Δs☉ = {Msun_ds:.4f} m;  12Δs/R☉ = {12*Msun_ds/Rsun/ARCSEC:.4f}\"  (GR 1.751\")")
print(f"   nur Takt-Anteil (c_eff aus 08, Brechzahl n = 1/sqrt(1-6Δs/r)): 6Δs/R☉ = {6*Msun_ds/Rsun/ARCSEC:.4f}\"  (halber Wert)")

def I_A(x, nu=nu_ftm):
    """α_A(b) = (12Δs/b)·I(b/r*),  I(x) = ∫_0^{π/2} cosθ·ν(cos²θ/x²) dθ"""
    th = 0.25 * math.pi * (XG + 1); w = 0.25 * math.pi * WG
    return float(np.sum(w * np.cos(th) * nu(np.cos(th)**2 / x**2)))

M = 1e11 * MS; ds = G * M / (3 * c**2); rs = math.sqrt(G * M / A0)
print(f"\n   Punktmasse M = 1e11 M☉: Δs = {ds:.4e} m = {ds/1.495978707e11:.1f} AE, r* = {rs/KPC:.3f} kpc")
print(f"   Tiefgrenze Option A: α_∞ = 6πΔs/r* = 2π·sqrt(GM a0)/c² = {6*math.pi*ds/rs/ARCSEC:.4f}\"")
print(f"   {'b/r*':>7s} {'b[kpc]':>8s} {'I(b/r*)':>9s} {'α_N=α_Bs':>10s} {'α_B_halb':>10s} {'α_A':>10s} {'Näherung A':>11s}   [arcsec]")
for x in [0.1, 0.3, 1, 3, 10, 30, 100]:
    b = x * rs; I = I_A(x); aN = 12 * ds / b
    approx = 6 * math.pi * ds / rs + 6 * ds / b
    print(f"   {x:7.1f} {b/KPC:8.2f} {I:9.4f} {aN/ARCSEC:10.5f} {0.5*aN*(1+I)/ARCSEC:10.5f} {aN*I/ARCSEC:10.5f} {approx/ARCSEC:11.5f}")
print("   Näherung A für b >> r*: α_A ≈ 6πΔs/r* + 6Δs/b")
print("\n   η(r) = M_dyn/M_lens (3D, Punktmasse), y = (r*/r)²:")
print(f"   {'r/r*':>7s} {'ν':>9s} {'η_A':>6s} {'η_B_halb=2ν/(1+ν)':>19s} {'η_B_strikt=ν':>13s} {'r/r*+1/2':>9s}")
for x in [0.1, 0.5, 1, 2, 5, 10, 30, 100]:
    n = float(nu_ftm(1 / x**2))
    print(f"   {x:7.1f} {n:9.4f} {1.0:6.2f} {2*n/(1+n):19.4f} {n:13.4f} {x+0.5:9.3f}")

# =============================================================================== 2
print("\n" + "=" * 100)
print("2) GALAXIE-GALAXIE-LINSEN: ΔΣ(R) für isolierte Galaxien (Baryonen als Hernquist-Kugel)")
print("   Annahme: a_H = 1.5 kpc·(M/1e10)^0.3 (Sterne + kaltes Gas, ohne heißen Halo). Punktmasse wäre ab 20 kpc fast gleich.")
Rk = np.array([10, 20, 30, 50, 100, 200, 300, 500, 700, 1000.])
R_si = Rk * KPC
for Msol in [1e10, 3e10, 1e11]:
    Mb = Msol * MS; a = 1.5 * KPC * (Msol / 1e10)**0.3
    gN = lambda r, Mb=Mb, a=a: G * Mb / (r + a)**2
    gA = lambda r, gN=gN: nu_ftm(gN(r) / A0) * gN(r)
    MN = lambda R, gN=gN: Mproj(R, gN)
    MA = lambda R, gA=gA: Mproj(R, gA)
    MBh = lambda R: 0.5 * (MN(R) + MA(R))
    dsN = DeltaSigma(R_si, MN) / SIG_PC2; dsA = DeltaSigma(R_si, MA) / SIG_PC2; dsB = DeltaSigma(R_si, MBh) / SIG_PC2
    rstar = math.sqrt(G * Mb / A0); vf = (G * Mb * A0)**0.25
    print(f"\n   M_bar = {Msol:.0e} M☉, a_H = {a/KPC:.2f} kpc, r* = {rstar/KPC:.2f} kpc, v_flach = (GMa0)^1/4 = {vf/1e3:.1f} km/s")
    print(f"   {'R[kpc]':>7s} {'ΔΣ_Newton':>10s} {'ΔΣ_B_halb':>10s} {'ΔΣ_A':>10s}  [M☉/pc²] {'g_bar/a0':>9s} {'4GΔΣ_A/g_RAR':>13s} {'4GΔΣ_Bh/g_RAR':>14s} {'4GΔΣ_N/g_RAR':>13s} {'A/N':>7s}")
    for i, R in enumerate(R_si):
        Mb3 = Mb * R**2 / (R + a)**2
        gbar = G * Mb3 / R**2; grar = float(nu_ftm(gbar / A0)) * gbar
        f = 4 * G * SIG_PC2 / grar
        print(f"   {Rk[i]:7.0f} {dsN[i]:10.4f} {dsB[i]:10.4f} {dsA[i]:10.4f}           {gbar/A0:9.2e} {f*dsA[i]:13.3f} {f*dsB[i]:14.3f} {f*dsN[i]:13.4f} {dsA[i]/dsN[i]:7.2f}")
print("\n   g_RAR = ν(g_bar/a0)·g_bar mit g_bar = G·M_bar(<R)/R² (3D, wie in der Linsen-RAR).")
print("   4GΔΣ ist die SIS-Umrechnung ΔΣ → g_obs (exakt für ρ ∝ r⁻²). Werte ≈1: Linsen folgen der RAR.")

# =============================================================================== 3
print("\n" + "=" * 100)
print("3) MASSEREICHER GALAXIENHAUFEN: Einstein-Radius und Linsenmasse")
zl = 0.3
Hz = H0 * E(zl); rhoc = 3 * Hz**2 / (8 * math.pi * G)
print(f"   z_l = {zl}; ρ_c(z_l) = {rhoc:.4e} kg/m³")
# baryonisches Modell
rc_g, rt_g, Mg1 = 150 * KPC, 3000 * KPC, 1e14 * MS
fbeta = lambda x: x - np.arctan(x)
Kg = Mg1 / fbeta(1000 / 150)
def Mgas(r): return Kg * fbeta(np.minimum(r, rt_g) / rc_g)
M_bcg, a_bcg = 1.5e12 * MS, 20 * KPC
M_gal, a_gal = 1.0e13 * MS, 300 * KPC
def Mbar(r): return Mgas(r) + M_bcg * r**2 / (r + a_bcg)**2 + M_gal * r**2 / (r + a_gal)**2
print("   Baryonen: Gas β-Modell β=2/3, r_c = 150 kpc, M_gas(<1 Mpc) = 1e14 M☉, abgeschnitten bei 3 Mpc;")
print("             BCG (mit ICL) Hernquist 1.5e12 M☉, a = 20 kpc; Mitgliedsgalaxien Hernquist 1e13 M☉, a = 300 kpc (Annahmen)")

def gN_cl(r): return G * Mbar(r) / r**2
def make_gA(nu, a0): return lambda r: nu(gN_cl(r) / a0) * gN_cl(r)
gA_cl = make_gA(nu_ftm, A0)

# NFW-Referenz (typische Linsenmasse eines massereichen Haufens)
def nfw(M200, cc):
    r200 = (3 * M200 / (4 * math.pi * 200 * rhoc))**(1 / 3); rs_ = r200 / cc
    dc = 200 / 3 * cc**3 / (math.log(1 + cc) - cc / (1 + cc)); rhos = dc * rhoc
    def M3(r):
        x = r / rs_; return 4 * math.pi * rhos * rs_**3 * (np.log(1 + x) - x / (1 + x))
    def M2(R):
        x = np.atleast_1d(R / rs_).astype(float); h = np.empty_like(x)
        lo = x < 1; hi = x > 1; eq = ~(lo | hi)
        h[lo] = np.log(x[lo] / 2) + 2 / np.sqrt(1 - x[lo]**2) * np.arctanh(np.sqrt((1 - x[lo]) / (1 + x[lo])))
        h[hi] = np.log(x[hi] / 2) + 2 / np.sqrt(x[hi]**2 - 1) * np.arctan(np.sqrt((x[hi] - 1) / (x[hi] + 1)))
        h[eq] = 1 + math.log(0.5)
        return 4 * math.pi * rhos * rs_**3 * h
    return M3, M2, r200

def einstein_R(M2fun, Sc, Rgrid=None):
    if Rgrid is None: Rgrid = np.geomspace(0.3, 800, 1500) * KPC
    mbar = M2fun(Rgrid) / (math.pi * Rgrid**2) - Sc
    idx = np.where(mbar > 0)[0]
    if len(idx) == 0: return 0.0
    i = idx[-1]
    if i == len(Rgrid) - 1: return Rgrid[-1]
    lo, hi = Rgrid[i], Rgrid[i + 1]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if M2fun(np.array([mid]))[0] / (math.pi * mid**2) - Sc > 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

M2_N = lambda R: Mproj(R, gN_cl)
M2_A = lambda R: Mproj(R, gA_cl)
M2_Bh = lambda R: 0.5 * (M2_N(R) + M2_A(R))
refs = [("NFW 1e15, c=4", *nfw(1e15 * MS, 4.0)), ("NFW 1e15, c=6", *nfw(1e15 * MS, 6.0)), ("NFW 2e15, c=5", *nfw(2e15 * MS, 5.0))]

Rt = np.array([50, 100, 250, 500, 1000.]) * KPC
print(f"\n   3D-Profile: {'r[kpc]':>7s} {'M_bar':>10s} {'g_N/a0':>8s} {'ν':>6s} {'M_dyn,A=νM_bar':>15s} {'M_NFW(1e15,c4)':>15s} {'Rest NFW/A':>10s}")
for r in Rt:
    y = gN_cl(r) / A0; n = float(nu_ftm(y))
    print(f"               {r/KPC:7.0f} {Mbar(r)/MS:10.3e} {y:8.3f} {n:6.3f} {n*Mbar(r)/MS:15.3e} {refs[0][1](r)/MS:15.3e} {refs[0][1](r)/(n*Mbar(r)):10.2f}")

print(f"\n   Projizierte Linsenmasse M_2D(<R) [M☉] (Sichtlinie bis ∞):")
print(f"   {'R[kpc]':>7s} {'Baryonen(N)':>12s} {'B_halb':>11s} {'Option A':>11s} {'A/bar':>6s} " + " ".join(f"{n[0]:>14s}" for n in refs) + f" {'Rest(1e15,c4)':>13s}")
M2N = M2_N(Rt); M2A = M2_A(Rt)
for i, R in enumerate(Rt):
    rr = [f[2](np.array([R]))[0] for f in refs]
    print(f"   {R/KPC:7.0f} {M2N[i]/MS:12.3e} {0.5*(M2N[i]+M2A[i])/MS:11.3e} {M2A[i]/MS:11.3e} {M2A[i]/M2N[i]:6.2f} " + " ".join(f"{x/MS:14.3e}" for x in rr) + f" {rr[0]/M2A[i]:13.2f}")

print("\n   Einstein-Radien (Radius, in dem die mittlere Konvergenz 1 ist):")
for zs in [1.0, 2.0]:
    Sc, Dl, Ds, Dls = Sigma_crit(zl, zs)
    print(f"   z_s = {zs}: Σ_crit = {Sc/SIG_PC2:.1f} M☉/pc², D_l = {Dl/MPC:.1f} Mpc, 1\" = {Dl*ARCSEC/KPC:.3f} kpc")
    for name, f in [("Baryonen Newton (=B_strikt)", M2_N), ("Option B_halb", M2_Bh), ("Option A (ν_FTM, a0=1.2e-10)", M2_A)] + [(n[0], n[2]) for n in refs]:
        RE = einstein_R(f, Sc)
        s = f"θ_E = {RE/Dl/ARCSEC:6.2f}\"  R_E = {RE/KPC:6.1f} kpc" if RE > 0 else "kein Einstein-Ring (mittlere Konvergenz < 1 überall)"
        print(f"      {name:34s} {s}")

Sc2, Dl2, _, _ = Sigma_crit(zl, 2.0)
print("\n   Nötige Linsenmasse für einen beobachteten Einstein-Radius (z_s = 2) gegenüber Option A:")
for th in [10, 20, 30, 45]:
    RE = th * ARCSEC * Dl2; Mreq = math.pi * RE**2 * Sc2; mA = M2_A(np.array([RE]))[0]; mN = M2_N(np.array([RE]))[0]
    print(f"   θ_E = {th:2d}\" (R_E = {RE/KPC:5.1f} kpc): M_2D nötig = {Mreq/MS:.3e} M☉, Option A = {mA/MS:.3e}, Baryonen = {mN/MS:.3e}"
          f"  ⇒  Rest = {Mreq/mA:4.2f}, nötig/Baryonen = {Mreq/mN:5.2f}")

print("\n   Vergleich der Übergangsfunktionen und von a0 (Option A), M_2D(<500 kpc) und Rest gegenüber NFW(1e15, c=4):")
M500ref = refs[0][2](np.array([500 * KPC]))[0]
a0z = A0_H0 * E(zl)
for lab, nu, a0v in [("ν_FTM = 1/(1-e^-√y), a0 = 1.2e-10", nu_ftm, A0), ("ν_simple, a0 = 1.2e-10", nu_simple, A0),
                     ("ν_standard, a0 = 1.2e-10", nu_standard, A0), ("ν_FTM, a0 = cH0/2π", nu_ftm, A0_H0),
                     (f"ν_FTM, a0(z=0.3) = cH(z)/2π = {a0z:.3e}", nu_ftm, a0z)]:
    g = make_gA(nu, a0v)
    m5 = Mproj(np.array([500 * KPC]), g)[0]; m1 = Mproj(np.array([100 * KPC]), g)[0]
    m1r = refs[0][2](np.array([100 * KPC]))[0]
    print(f"   {lab:46s} M_2D(500) = {m5/MS:.3e}  Rest(500) = {M500ref/m5:5.2f}   M_2D(100) = {m1/MS:.3e}  Rest(100) = {m1r/m1:5.2f}")
print("   Rest := M_2D,NFW-Referenz / M_2D,Option A. Die Referenz steht für die beobachtete Linsenmasse (Annahme).")

# =============================================================================== 4
print("\n" + "=" * 100)
print("4) BULLET CLUSTER 1E 0657-56: 3D-Operator (QUMOND-artig) auf dem Gitter, projiziert zu κ(x,y)")
Gk = G * MS / KPC / 1e6                 # kpc (km/s)² / M☉
a0k = A0 * KPC / 1e6                    # (km/s)²/kpc
Nx, Ny, Nz = 192, 192, 96
hx = hy = 20.0; hz = 40.0               # kpc
xc = 360.0
xs = xc + (np.arange(Nx) - Nx // 2) * hx
ys = (np.arange(Ny) - Ny // 2) * hy
zs_ = (np.arange(Nz) - Nz // 2) * hz
dV = hx * hy * hz
print(f"   Gitter {Nx}x{Ny}x{Nz}, Zellen {hx}x{hy}x{hz} kpc, Box {Nx*hx/1e3:.2f} x {Ny*hy/1e3:.2f} x {Nz*hz/1e3:.2f} Mpc,")
print("   isolierte Randbedingung durch Nullauffüllung (Hockney-Eastwood), 1/r-Kern im doppelten Gitter.")
X, Y, Z = np.meshgrid(xs, ys, zs_, indexing="ij")

def comp_density(kind, x0, y0, M, p, rt):
    r = np.sqrt((X - x0)**2 + (Y - y0)**2 + Z**2)
    if kind == "beta":            # β = 2/3
        rho = 1.0 / (1 + (r / p)**2)
    else:                         # Hernquist, weich im Zentrum (halbe Zelle)
        rs2 = np.sqrt(r**2 + (0.5 * hx)**2)
        rho = 1.0 / (rs2 * (rs2 + p)**3)
    rho[r > rt] = 0.0
    return rho * (M / (rho.sum() * dV))

# Green-Funktion (einmal)
t0 = time.time()
px = np.arange(2 * Nx); px = np.where(px < Nx, px, px - 2 * Nx) * hx
py = np.arange(2 * Ny); py = np.where(py < Ny, py, py - 2 * Ny) * hy
pz = np.arange(2 * Nz); pz = np.where(pz < Nz, pz, pz - 2 * Nz) * hz
PX, PY, PZ = np.meshgrid(px, py, pz, indexing="ij")
Rk_ = np.sqrt(PX**2 + PY**2 + PZ**2)
Rk_[0, 0, 0] = 1.0
Kern = 1.0 / Rk_
Kern[0, 0, 0] = 2.38 / (hx * hy * hz)**(1 / 3)     # Zellmittel von 1/r (Würfelnäherung)
del PX, PY, PZ, Rk_
Kf = np.fft.rfftn(Kern); del Kern

def potential(rho):
    pad = np.zeros((2 * Nx, 2 * Ny, 2 * Nz)); pad[:Nx, :Ny, :Nz] = rho * dV
    phi = np.fft.irfftn(np.fft.rfftn(pad) * Kf, s=pad.shape, axes=(0, 1, 2))[:Nx, :Ny, :Nz]
    return -Gk * phi

def phantom_sigma(rho, nu=nu_ftm, a0=a0k):
    """Σ_ph(x,y) = ∫ρ_ph dz, ρ_ph = -∇·[(ν-1) g_N]/(4πG). Zusatz aus Struktur, nicht aus Quelle."""
    phi = potential(rho)
    gx, gy, gz = np.gradient(-phi, hx, hy, hz)
    gm = np.sqrt(gx**2 + gy**2 + gz**2)
    f = nu(gm / a0) - 1.0
    div = np.gradient(f * gx, hx, axis=0) + np.gradient(f * gy, hy, axis=1) + np.gradient(f * gz, hz, axis=2)
    rph = -div / (4 * math.pi * Gk)
    return rph.sum(axis=2) * hz, rph, gm

print(f"   FFT-Kern vorbereitet in {time.time()-t0:.1f} s")

# --- Validierung an einer Kugel (β-Gas + BCG) gegen die 1D-Rechnung
print("\n   Validierung: kugelsymmetrischer Haufen (β-Gas r_c = 150 kpc, 1e14 M☉ bis 1.5 Mpc; BCG 1e12, a=25) im Gitterzentrum")
rho_v = comp_density("beta", xc, 0, 1e14, 150.0, 1500.0) + comp_density("hern", xc, 0, 1e12, 25.0, 1500.0)
Sph_v, _, _ = phantom_sigma(rho_v)
SN_v = rho_v.sum(axis=2) * hz
XX, YY = np.meshgrid(xs, ys, indexing="ij")
def aperture(S, x0, y0, R): return float(S[(XX - x0)**2 + (YY - y0)**2 <= R**2].sum() * hx * hy)
# 1D-Referenz: Masse im Zylinder (Radius R, |z| < Zmax) über den Gauß-Satz, g aus der 1D-Masse
Zmax = Nz * hz / 2
def M1d(r):
    x = np.minimum(r, 1500.0) / 150.0
    Kv = 1e14 / fbeta(1500.0 / 150.0)
    return Kv * fbeta(x) + 1e12 * r**2 / (r + 25.0)**2
def cyl_mass(R, gfun):
    zz = 0.5 * Zmax * (XG + 1); wz = 0.5 * Zmax * WG
    rr = np.sqrt(R**2 + zz**2); side = 2 * math.pi * R * 2 * np.sum(wz * gfun(rr) * R / rr)
    Rp = 0.5 * R * (XG + 1); wR = 0.5 * R * WG
    rr2 = np.sqrt(Rp**2 + Zmax**2); cap = 2 * np.sum(wR * 2 * math.pi * Rp * gfun(rr2) * Zmax / rr2)
    return (side + cap) / (4 * math.pi * Gk)
g1N = lambda r: Gk * M1d(r) / r**2
g1A = lambda r: nu_ftm(g1N(r) / a0k) * g1N(r)
print(f"   {'R[kpc]':>7s} {'M_N Gitter':>11s} {'M_N 1D':>11s} {'M_A Gitter':>11s} {'M_A 1D':>11s} {'Abw. A':>7s}")
for R in [100, 250, 500]:
    mng = aperture(SN_v, xc, 0, R); mag = mng + aperture(Sph_v, xc, 0, R)
    mn1 = cyl_mass(R, g1N); ma1 = cyl_mass(R, g1A)
    print(f"   {R:7.0f} {mng:11.3e} {mn1:11.3e} {mag:11.3e} {ma1:11.3e} {100*(mag/ma1-1):6.1f}%")
del rho_v, Sph_v

# --- Baryonisches Modell des Bullet Clusters
print("\n   Baryonisches Modell (Achse Ost→West = +x; Geometrie gerundet nach Clowe et al. 2006, Abb. 1; ±50 kpc unsicher):")
comps = [
    # name,               art,   x,     y, Masse[M☉], Parameter[kpc], r_t[kpc], Gruppe
    ("Haupt-BCG",          "hern",   0.0, 0.0, 1.0e12,  25.0, 1500.0, "gal_main"),
    ("Haupt-Galaxien",     "hern",   0.0, 0.0, 1.0e13, 250.0, 1500.0, "gal_main"),
    ("Sub-BCG",            "hern", 720.0, 0.0, 6.0e11,  20.0, 1500.0, "gal_sub"),
    ("Sub-Galaxien",       "hern", 720.0, 0.0, 5.0e12, 150.0, 1500.0, "gal_sub"),
    ("Hauptgas",           "beta", 200.0, 0.0, 1.5e14, 200.0, 1500.0, "gas_main"),
    ("Bullet-Gas",         "beta", 560.0, 0.0, 2.0e13,  60.0,  500.0, "gas_sub"),
]
for n, k, x0, y0, Mc, p, rt, gr in comps:
    print(f"     {n:15s} {'β=2/3, r_c' if k=='beta' else 'Hernquist a'} = {p:5.0f} kpc  x = {x0:5.0f} kpc  M = {Mc:.2e} M☉ (r < {rt:.0f} kpc)")
Mgas_tot = sum(cc[4] for cc in comps if cc[7].startswith("gas")); Mst_tot = sum(cc[4] for cc in comps if cc[7].startswith("gal"))
print(f"     Gas gesamt {Mgas_tot:.2e} M☉, Sterne gesamt {Mst_tot:.2e} M☉, Gas/Sterne = {Mgas_tot/Mst_tot:.1f}")
print("     Massen sind Modellannahmen in der Größenordnung von Clowe et al. 2006 / Bradač et al. 2006 (nicht exakt übernommen).")

def build(shift_gas_phi=None, star_f=1.0, gas_only=None):
    """Dichte; mit shift_gas_phi = φ liegen die Gaskomponenten an X_gal + φ(X_gas - X_gal) (Gedächtnislage).
    star_f skaliert die Sternmassen; gas_only=True/False liefert nur Gas bzw. nur Sterne."""
    rho = np.zeros((Nx, Ny, Nz))
    gal_x = {"gas_main": 0.0, "gas_sub": 720.0}
    for n, k, x0, y0, Mc, p, rt, gr in comps:
        isgas = gr.startswith("gas")
        if gas_only is True and not isgas: continue
        if gas_only is False and isgas: continue
        xx = x0
        if shift_gas_phi is not None and isgas:
            xx = gal_x[gr] + shift_gas_phi * (x0 - gal_x[gr])
        rho += comp_density(k, xx, y0, Mc * (1.0 if isgas else star_f), p, rt)
    return rho

t0 = time.time()
rho_b = build()
Sig_N = rho_b.sum(axis=2) * hz
Sig_ph, rph, gm = phantom_sigma(rho_b)
print(f"   QUMOND-Lauf in {time.time()-t0:.1f} s; ρ_ph min/max = {rph.min():.3e} / {rph.max():.3e} M☉/kpc³ "
      f"(Anteil negativer Phantommasse: {abs(rph[rph<0].sum())/rph[rph>0].sum()*100:.2f} %)")
j0 = Ny // 2; k0 = Nz // 2
for lab, x0 in [("Haupt-BCG", 0), ("Hauptgas", 200), ("Mitte", 360), ("Bullet-Gas", 560), ("Sub-BCG", 720)]:
    i = int(round((x0 - xs[0]) / hx))
    print(f"     g_N/a0 in der Mittelebene bei {lab:10s} (x={x0:4d}): {gm[i,j0,k0]/a0k:6.3f}   ν = {float(nu_ftm(gm[i,j0,k0]/a0k)):5.2f}")
del rph

zl_b = 0.296
Sc_b, Dl_b, _, _ = Sigma_crit(zl_b, 1.0)
Sc_b_k = Sc_b / (MS / KPC**2)
print(f"   z_l = {zl_b}, z_s = 1: Σ_crit = {Sc_b_k:.3e} M☉/kpc², 1\" = {Dl_b*ARCSEC/KPC:.3f} kpc")

KX = np.fft.fftfreq(Nx, d=hx) * 2 * math.pi; KY = np.fft.fftfreq(Ny, d=hy) * 2 * math.pi
KX2, KY2 = np.meshgrid(KX**2, KY**2, indexing="ij")
def smooth(S, sig):
    """Gauß-Glättung (σ in kpc) – Auflösung einer Linsen-Rekonstruktion."""
    if sig <= 0: return S
    return np.real(np.fft.ifft2(np.fft.fft2(S) * np.exp(-0.5 * sig**2 * (KX2 + KY2))))

SIG_S = 50.0      # Hauptwert der Glättung [kpc]
def peaks_row(S, sig=SIG_S):
    row = smooth(S, sig)[:, j0]
    out = {}
    for lab, lo, hi in [("main", -300, 360), ("sub", 360, 1000)]:
        m = (xs >= lo) & (xs < hi); idx = np.where(m)[0]
        i = idx[np.argmax(row[idx])]
        y1, y2, y3 = row[i - 1], row[i], row[i + 1]
        d = 0.5 * (y1 - y3) / (y1 - 2 * y2 + y3) if (y1 - 2 * y2 + y3) != 0 else 0.0
        out[lab] = (xs[i] + d * hx, y2 / Sc_b_k)
    return out

maps = {"Newton (Baryonen) = B_strikt": Sig_N, "Option B_halb": Sig_N + 0.5 * Sig_ph, "Option A (Operator)": Sig_N + Sig_ph,
        "nur Phantomanteil (A - Newton)": Sig_ph}
print("\n   Konvergenzmaxima auf der Stoßachse (y = 0); Karte mit Gauß σ geglättet (Auflösung der Rekonstruktion):")
print("   (Galaxien bei x = 0 / 720 kpc, Gas bei x = 200 / 560 kpc)")
for sg in [0.0, 50.0, 100.0]:
    print(f"   σ = {sg:5.0f} kpc  {'Karte':32s} {'x_peak Haupt':>12s} {'κ':>6s} {'x_peak Sub':>11s} {'κ':>6s}")
    for k_, S in maps.items():
        p = peaks_row(S, sg)
        print(f"                {k_:32s} {p['main'][0]:12.1f} {p['main'][1]:6.3f} {p['sub'][0]:11.1f} {p['sub'][1]:6.3f}")
print("   σ = 0: Pixelspitze der BCG (20 kpc) – kein Linsenobservabel; maßgeblich sind σ = 50–100 kpc.")

# Empfindlichkeit: Sternmassen skaliert
print("   Empfindlichkeit gegenüber der Sternmasse (Option A; σ = 50 kpc, in Klammern σ = 100 kpc):")
for sf in [0.25, 0.5, 1.0, 2.0]:
    rho_s = build(star_f=sf); SN_s = rho_s.sum(axis=2) * hz; Sph_s, _, _ = phantom_sigma(rho_s); del rho_s
    p = peaks_row(SN_s + Sph_s, 50.0); p2 = peaks_row(SN_s + Sph_s, 100.0)
    Sst = SN_s - build(gas_only=True).sum(axis=2) * hz
    print(f"     Sterne x{sf:<4} (Sterne in 100 kpc um Haupt-/Sub-BCG: {aperture(Sst, 0, 0, 100):.1e} / {aperture(Sst, 720, 0, 100):.1e} M☉)"
          f"  x_peak Haupt = {p['main'][0]:6.1f} ({p2['main'][0]:6.1f})   x_peak Sub = {p['sub'][0]:6.1f} ({p2['sub'][0]:6.1f})")
    del SN_s, Sph_s, Sst
print("   Clowe et al. 2006 geben in 100-kpc-Aperturen um die BCGs Sternmassen von nur ~5e11 M☉ an (UNSICHER) ≈ Faktor x0.25 hier.")

# Rest: 100-kpc-Aperturen
print("\n   Aperturmassen (R = 100 kpc) [M☉]:")
clowe = {"Haupt-BCG": (0, 2.8e13), "Hauptgas": (200, 2.3e13), "Bullet-Gas": (560, 1.8e13), "Sub-BCG": (720, 2.0e13)}
print(f"   {'Ort':11s} {'Gas':>10s} {'Sterne':>10s} {'Newton':>10s} {'Option A':>10s} {'A/Bar':>6s} {'M_κ Clowe*':>11s} {'Rest':>6s}")
Sig_gas = sum(comp_density(k, x0, y0, Mc, p, rt).sum(axis=2) * hz for n, k, x0, y0, Mc, p, rt, gr in comps if gr.startswith("gas"))
Sig_st = Sig_N - Sig_gas
for lab, (x0, mk) in clowe.items():
    mg = aperture(Sig_gas, x0, 0, 100); mst = aperture(Sig_st, x0, 0, 100)
    mn = aperture(Sig_N, x0, 0, 100); ma = mn + aperture(Sig_ph, x0, 0, 100)
    print(f"   {lab:11s} {mg:10.2e} {mst:10.2e} {mn:10.2e} {ma:10.2e} {ma/mn:6.2f} {mk:11.2e} {mk/ma:6.2f}")
print("   *M_κ: Richtwerte der Größenordnung nach Clowe et al. 2006, Tab. 2 (aus dem Gedächtnis, UNSICHER, am Original prüfen).")
print("    Beobachtet: M_κ(100 kpc) ~ 2–3e13 M☉ an beiden Galaxienpeaks, Plasma dort nur ~(3–6)e12, Sterne ~5e11.")

print("\n   Profil κ(x, y=0) (z_s = 1), geglättet mit σ = 50 kpc:")
print(f"   {'x[kpc]':>7s} {'κ_Gas':>7s} {'κ_Sterne':>8s} {'κ_bar':>7s} {'κ_A':>7s}")
sGas = smooth(Sig_gas, 50.0); sSt = smooth(Sig_st, 50.0); sA = smooth(Sig_N + Sig_ph, 50.0)
for x0 in range(-300, 1041, 60):
    i = int(round((x0 - xs[0]) / hx))
    print(f"   {x0:7d} {sGas[i,j0]/Sc_b_k:7.4f} {sSt[i,j0]/Sc_b_k:8.4f} {(sGas[i,j0]+sSt[i,j0])/Sc_b_k:7.4f} {sA[i,j0]/Sc_b_k:7.4f}")

# =============================================================================== 5
print("\n" + "=" * 100)
print("5) GEDÄCHTNIS-VARIANTE: Operatorfeld läuft der Zustandsgeschichte nach")
print("   Gedächtnisgeschwindigkeit u(t) = λ∫e^{-λ(t-t')} v(t') dt' (Kern aus 07), Gedächtnislage dX_mem/dt = u.")
print("   Struktur wird aus ρ_mem gebildet:  ∇²Φ = 4πGρ_bar + ∇·[(ν(|∇Φ_N[ρ_mem]|/a0) - 1) ∇Φ_N[ρ_mem]]")
print("   Stationär gilt ρ_mem = ρ_bar (Galaxien, Sonnensystem unverändert). Stoßfreie Galaxien: v = const ⇒ X_mem = X_gal.")
print("   Gas, zur Zeit 0 abgebremst:  (X_mem,gas - X_gal)/(X_gas - X_gal) = φ(λt)")
print("     sofortige Abbremsung:    φ_s(x) = 1 - (1 - e^{-x})/x          (x = λt)")
print("     konstante Abbremsung:    φ_k(x) = 1 - 2/x + 2(1 - e^{-x})/x²")

def phi_s(x): return 1 - (-np.expm1(-x)) / x if x > 1e-8 else x / 2
def phi_k(x): return 1 - 2 / x + 2 * (-np.expm1(-x)) / x**2 if x > 1e-4 else x / 3
def inv(fun, target):
    lo, hi = 1e-6, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if fun(mid) < target: lo = mid
        else: hi = mid
    return math.sqrt(lo * hi)

phis = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
t_obs = 150 * MYR; v_coll = 4500e3
print(f"   Durchgangsweg 4500 km/s · 150 Myr = {v_coll*t_obs/KPC:.0f} kpc; Kreuzungszeit 720 kpc/4500 km/s = {720*KPC/v_coll/MYR:.0f} Myr")

def crossing(vals, crit):
    for i in range(len(phis) - 1):
        if vals[i] <= crit < vals[i + 1]:
            return phis[i] + (crit - vals[i]) / (vals[i + 1] - vals[i]) * (phis[i + 1] - phis[i])
    return None

results = {}
for sf in [1.0, 0.25]:
    res = []
    SN_true = build(star_f=sf).sum(axis=2) * hz
    print(f"\n   Scan über φ, Sternmassen x{sf}; Karten mit σ = 50 kpc geglättet (in Klammern β bei σ = 100 kpc)")
    print(f"   {'φ':>5s} {'x_peak Haupt':>12s} {'β_Haupt':>8s} {'x_peak Sub':>11s} {'β_Sub':>7s} {'κ_Haupt':>8s} {'κ_Sub':>7s} {'(β_H, β_S) σ=100':>18s}   β = (x_peak - x_gal)/(x_gas - x_gal)")
    for ph in phis:
        rho_m = build(shift_gas_phi=ph, star_f=sf)
        Sph_m, _, _ = phantom_sigma(rho_m)
        S = SN_true + Sph_m
        p = peaks_row(S); p2 = peaks_row(S, 100.0)
        bm = p["main"][0] / 200.0; bs = (p["sub"][0] - 720) / (560 - 720)
        bm2 = p2["main"][0] / 200.0; bs2 = (p2["sub"][0] - 720) / (560 - 720)
        res.append((ph, bm, bs, bm2, bs2))
        print(f"   {ph:5.2f} {p['main'][0]:12.1f} {bm:8.3f} {p['sub'][0]:11.1f} {bs:7.3f} {p['main'][1]:8.3f} {p['sub'][1]:7.3f}   ({bm2:5.2f}, {bs2:5.2f})")
        if ph == 0.0:
            print("         Aperturmassen R=100 kpc bei φ=0: " + ", ".join(
                f"{lab} {aperture(S, x0, 0, 100):.2e} (Rest {mk/aperture(S, x0, 0, 100):.2f})" for lab, (x0, mk) in clowe.items() if "BCG" in lab))
        del rho_m, Sph_m
    res = np.array(res); results[sf] = res
    for lab, col in [("σ = 50 kpc", (1, 2)), ("σ = 100 kpc", (3, 4))]:
        bmax = np.maximum(res[:, col[0]], res[:, col[1]])
        for bc in [0.25, 0.4]:
            pc_ = crossing(bmax, bc)
            if pc_ is None:
                print(f"   {lab}, β ≤ {bc}: " + ("schon bei φ = 0 verletzt" if bmax[0] > bc else "für alle φ erfüllt"))
            else:
                xs_ = inv(phi_s, pc_); xk_ = inv(phi_k, pc_)
                print(f"   {lab}, β ≤ {bc}: φ ≤ {pc_:.3f} ⇒ λt ≤ {xs_:.3f} (sofort) / {xk_:.3f} (konstant) ⇒ "
                      f"1/λ ≥ {150/xs_:5.0f} Myr (sofort) / {150/xk_:5.0f} Myr (konstant) bei t = 150 Myr")
    del SN_true

print("\n   Zusatz: Rest (fehlender Faktor) ebenfalls als Gedächtnisstruktur an den stoßfreien Komponenten, φ = 0:")
print("   κ = κ_N[ρ_bar] + s·κ_ph[ρ_mem];  s = 1 ist der reine Operator, s > 1 enthält den Rest als Struktur")
for sf in [1.0, 0.25]:
    SN_true = build(star_f=sf).sum(axis=2) * hz
    Sph0, _, _ = phantom_sigma(build(shift_gas_phi=0.0, star_f=sf))
    for s in [1.0, 1.5, 2.0, 3.0]:
        S = SN_true + s * Sph0
        p = peaks_row(S); p2 = peaks_row(S, 100.0)
        bm = p["main"][0] / 200.0; bs = (p["sub"][0] - 720) / (560 - 720)
        bm2 = p2["main"][0] / 200.0; bs2 = (p2["sub"][0] - 720) / (560 - 720)
        print(f"   Sterne x{sf:<4} s = {s:3.1f}: β_Haupt = {bm:5.2f} ({bm2:5.2f}), β_Sub = {bs:5.2f} ({bs2:5.2f})   "
              f"M(100 kpc) Haupt-BCG = {aperture(S, 0, 0, 100):.2e}, Sub-BCG = {aperture(S, 720, 0, 100):.2e}")
    del SN_true, Sph0

res = results[1.0]; bmax = np.maximum(res[:, 1], res[:, 2])
phi_crit = crossing(bmax, 0.25)
def beta_of_phi(ph): return float(np.interp(ph, res[:, 0], bmax))
print("\n   Vorhersage β(t) = Versatz(κ-Galaxien)/Versatz(Gas-Galaxien) als Funktion der Zeit seit dem Stoß")
print("   (Bullet-Geometrie, Sterne x1, σ = 50 kpc, sofortige Abbremsung; ΛCDM-stoßfrei: β = 0; ohne Gedächtnis: β(φ=1))")
tl = [0.05, 0.1, 0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
lams = [300, 600, 1000, 3000]
print("   t[Gyr] " + " ".join(f"{'1/λ='+str(L)+' Myr':>15s}" for L in lams) + f" {'ohne Gedächtnis':>16s} {'ΛCDM(stoßfrei)':>15s}")
for t in tl:
    print(f"   {t:6.2f} " + " ".join(f"{beta_of_phi(phi_s(t*1e3/L)):15.3f}" for L in lams) + f" {beta_of_phi(1.0):16.3f} {0.0:15.3f}")
print("   Versatz κ–Gas = (1 - β)·Versatz Gas–Galaxien.")

print("\n   Zeitmarken (Literatur, UNSICHER, prüfen): Bullet ~0.15 Gyr; El Gordo ~0.5 Gyr; Musket Ball ~0.7 Gyr (Dawson 2013).")
for name, t in [("Bullet", 0.15), ("El Gordo", 0.5), ("Musket Ball", 0.7), ("1 Gyr", 1.0)]:
    print(f"   {name:12s} t = {t:4.2f} Gyr: " + "  ".join(f"1/λ={L:4d} Myr: φ={phi_s(t*1e3/L):.2f}, β={beta_of_phi(phi_s(t*1e3/L)):.2f}" for L in lams))
print("\n   Nötige Nachlaufzeit, damit β ≤ 0.25 (Sterne x1, σ = 50 kpc) auch t Gyr nach dem Stoß gilt:")
if phi_crit is not None:
    for t in [0.15, 0.5, 0.7, 1.0, 2.0]:
        print(f"   t = {t:4.2f} Gyr  ⇒  1/λ ≥ {t*1e3/inv(phi_s, phi_crit):6.0f} Myr (sofort) / {t*1e3/inv(phi_k, phi_crit):6.0f} Myr (konstant)")
print("\nFertig.")
