# -*- coding: utf-8 -*-
"""
Kovariante Lücken der Wirkung (00, Punkt 7): Wirbelanteil der örtlichen Regel, Gewicht w und Übergabe λ_m.
Rechnung, noch nicht im Modell. Alle Formeln und Konstanten aus ftm_kern.

A  Wirbelanteil: W = ν(|g_N|)·g_N − g_F (g_F = −∇Φ_F Feldform, der wirbelfreie Teil).
   Kovariant: J_μ = h_μ^ν(ν∇_νσ − a_ν) ist schon ein Feld der Wirkung; die Faltungsgleichung sagt D·J = 0,
   statisch c²J = −W. Der Wirbelanteil ist also genau der Strom der Faltungsgleichung.
   Milchstraße: Exponentialscheibe 5·10¹⁰ M☉, R_d = 3 kpc, h = 0,3 kpc, Hernquist-Bulge 10¹⁰ M☉ (0,6 kpc).
   Löser: Finite Volumen in (R, z ≥ 0), gleiche Matrix für Newton und Feldform; Prüfung an der Plummer-Kugel (W = 0).
   Energie: E = v²/2 + Φ_F, dE/dt = W·v ist genau die mit dem Medium getauschte Arbeit.
B  Gewicht w in der Wirkung: Faktor (1 − w) vor Q_s bzw. ν − 1; Größen bei Rekombination und heute; Preisach als
   lineares Programm (Hilfsfeld mit einseitiger Nebenbedingung) gegen das laufende Minimum.
C  Übergabe λ_m: Energierate der Mode gegen den expliziten w-Term; gelesene Masse von Mode und Regel um R_Q.
"""
import sys, io, math, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import linprog
import mpmath as mp
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ftm_kern import (G, c, pc, M_sun, A0, Gyr, nu, gewicht_w, tau_D, H0_si, E_hintergrund, a0_frueh_relativ,
                      Q_KOSM, q_KOSM, OMEGA_B, OMEGA_M, omega_vakuum, r_stern, uebergabezeit)

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []
def P(s=""):
    out.append(s); print(s, flush=True)

kpc = 1e3 * pc
GK = G * M_sun / (kpc * 1e6)                  # G in kpc·(km/s)²/M☉
A0K = A0 * kpc / 1e6                          # a₀ in (km/s)²/kpc
TU = kpc / 1e3 / (1e6 * 3.15576e7)            # Zeiteinheit kpc/(km/s) in Myr
nuK = lambda g: nu(np.maximum(g, 1e-300), A0K)

# ============================================================================ A  Löser
def stretched(N, L, h0):
    lo, hi = 1e-3, 30.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if L * math.sinh(m / N) / math.sinh(m) > h0: lo = m
        else: hi = m
    return L * np.sinh(m * np.arange(N + 1) / N) / math.sinh(m)

class Gitter:
    def __init__(self, Rf, zf):
        self.Rf, self.zf = Rf, zf
        self.Rc, self.zc = 0.5 * (Rf[1:] + Rf[:-1]), 0.5 * (zf[1:] + zf[:-1])
        self.NR, self.NZ = len(self.Rc), len(self.zc)
        self.dz, self.dR2 = np.diff(zf), np.diff(Rf**2)
        self.vol = np.pi * self.dR2[:, None] * self.dz[None, :]
        # innere Flächen: radial (i|i+1, j) und axial (i, j|j+1); Kopplung A/d
        self.cR = 2 * np.pi * Rf[1:-1, None] * self.dz[None, :] / np.diff(self.Rc)[:, None]      # (NR−1, NZ)
        self.cZ = np.pi * self.dR2[:, None] / np.diff(self.zc)[None, :]                          # (NR, NZ−1)
        idx = np.arange(self.NR * self.NZ).reshape(self.NR, self.NZ)
        self.idx = idx
        rows, cols, vals = [], [], []
        for (a, b, cc) in ((idx[:-1, :], idx[1:, :], self.cR), (idx[:, :-1], idx[:, 1:], self.cZ)):
            a, b, cc = a.ravel(), b.ravel(), cc.ravel()
            rows += [a, a, b, b]; cols += [a, b, b, a]; vals += [-cc, cc, -cc, cc]
        L = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.NR * self.NZ,) * 2).tolil()
        L[0, :] = 0; L[0, 0] = 1.0                    # Konstante festlegen (reiner Neumann-Rand)
        self.lu = splu(L.tocsc())

    def solve(self, rhs):
        r = rhs.ravel().copy(); r[0] = 0.0
        return self.lu.solve(r).reshape(self.NR, self.NZ)

    def newton(self, rho, Mtot):
        """Σ_innen (A/d)(Φ_nb − Φ_k) = 4πG·M_k − Randfluss; Randfluss aus dem Monopol, auf Gauß normiert."""
        M = rho * self.vol
        rhs = 4 * np.pi * GK * M
        bR = np.zeros_like(M); bZ = np.zeros_like(M)
        Rm, zm = self.Rf[-1], self.zf[-1]
        r = np.hypot(Rm, self.zc); bR[-1, :] = 2 * np.pi * Rm * self.dz * GK * Mtot / r**2 * (Rm / r)
        r = np.hypot(self.Rc, zm); bZ[:, -1] = np.pi * self.dR2 * GK * Mtot / r**2 * (zm / r)
        b = bR + bZ; b *= 4 * np.pi * GK * M.sum() / b.sum()
        return self.solve(rhs - b)

    def face_g(self, phi):
        """Normale Feldkomponenten auf den inneren Flächen, g = −∇Φ, und |g| dort (Tangente gemittelt)."""
        gRf = -(phi[1:, :] - phi[:-1, :]) / np.diff(self.Rc)[:, None]
        gZf = -(phi[:, 1:] - phi[:, :-1]) / np.diff(self.zc)[None, :]
        gRc, gZc = self.center_g(phi)
        tZ_R = 0.5 * (gZc[1:, :] + gZc[:-1, :]); tR_Z = 0.5 * (gRc[:, 1:] + gRc[:, :-1])
        return gRf, gZf, np.hypot(gRf, tZ_R), np.hypot(gZf, tR_Z)

    def center_g(self, phi):
        gR = -np.gradient(phi, self.Rc, axis=0)
        gZ = -np.gradient(phi, self.zc, axis=1)
        gZ[:, 0] = -(phi[:, 1] - phi[:, 0]) / (self.zc[1] + self.zc[0])      # Spiegelzelle bei −z₀
        return gR, gZ

    def feldform(self, phiN):
        """Σ_innen (A/d)(Φ_F,nb − Φ_F,k) = Σ_innen (A/d)·ν_f·(Φ_N,nb − Φ_N,k)."""
        gRf, gZf, aR, aZ = self.face_g(phiN)
        fR = self.cR * nuK(aR) * (phiN[1:, :] - phiN[:-1, :])
        fZ = self.cZ * nuK(aZ) * (phiN[:, 1:] - phiN[:, :-1])
        rhs = np.zeros_like(phiN)
        rhs[:-1, :] += fR; rhs[1:, :] -= fR
        rhs[:, :-1] += fZ; rhs[:, 1:] -= fZ
        return self.solve(rhs)

def felder(gi, rho, Mtot):
    phiN = gi.newton(rho, Mtot)
    phiF = gi.feldform(phiN)
    gR, gZ = gi.center_g(phiN)
    n = nuK(np.hypot(gR, gZ))
    fR, fZ = gi.center_g(phiF)
    return dict(phiN=phiN, phiF=phiF, gR=gR, gZ=gZ, aR=n * gR, aZ=n * gZ, fR=fR, fZ=fZ, WR=n * gR - fR, WZ=n * gZ - fZ)

P("=" * 100)
P("A  Wirbelanteil der örtlichen Regel und Arbeitsaustausch mit dem Medium")
P("=" * 100)
P(f"   a₀ = {A0:.2e} m/s² = {A0K:.4f} (km/s)²/kpc,  G = {GK:.5e} kpc(km/s)²/M☉,  Zeiteinheit {TU:.1f} Myr")

# ---- Prüfung an der Plummer-Kugel: Feldform = örtliche Regel, W = 0
Rf = stretched(220, 400.0, 0.03); zf = stretched(220, 400.0, 0.03)
gi = Gitter(Rf, zf)
Mp, bp = 6e10, 3.0
Rc, zc = np.meshgrid(gi.Rc, gi.zc, indexing="ij")
rho_p = 3 * Mp / (4 * np.pi * bp**3) * (1 + (Rc**2 + zc**2) / bp**2) ** -2.5
rho_p *= (Mp / 2) / (rho_p * gi.vol).sum()
F = felder(gi, rho_p, Mp)
P("\nA0 Prüfung des Lösers: Plummer-Kugel 6·10¹⁰ M☉, b = 3 kpc (Feldform = Regel, W = 0)")
P("     R [kpc]   g_N Gitter / exakt     a Regel     a Feldform    |W|/a")
for R in (1.0, 3.0, 8.0, 20.0, 60.0):
    i = np.searchsorted(gi.Rc, R); j = 0
    r = math.hypot(gi.Rc[i], gi.zc[j]); gex = GK * Mp * r / (r**2 + bp**2) ** 1.5
    a = math.hypot(F["aR"][i, j], F["aZ"][i, j]); f = math.hypot(F["fR"][i, j], F["fZ"][i, j])
    P(f"     {gi.Rc[i]:7.2f}   {abs(F['gR'][i, j]) / gex:.5f}          {a:9.3f}   {f:9.3f}     {math.hypot(F['WR'][i, j], F['WZ'][i, j]) / a:.2e}")

# ---- Prüfung an der dünnen Exponentialscheibe ohne Bulge (Vergleichswerte des Modells, 07: R_d/4 0,82 / 0,53, R_d 0,92 / 0,79)
from scipy.special import i0, i1, k0, k1
P("\nA0b Dünne Exponentialscheibe 5·10¹⁰ M☉ ohne Bulge (h = 0,01 R_d), Feldform/Regel in der Mittelebene")
rs5 = r_stern(5e10 * M_sun) / kpc
for Rdr in (0.3, 1.0, 3.0):
    Rd_ = Rdr * rs5
    g2 = Gitter(stretched(260, 60 * Rd_, 0.004), stretched(260, 60 * Rd_, 0.004))
    R2, Z2 = np.meshgrid(g2.Rc, g2.zc, indexing="ij")
    r2 = np.exp(-R2 / Rd_) * np.exp(-Z2 / (0.01 * Rd_)); r2 *= 2.5e10 / (r2 * g2.vol).sum()
    F2 = felder(g2, r2, 5e10)
    row = []
    for x in (0.25, 1.0, 3.5):
        i = np.searchsorted(g2.Rc, x * Rd_); y = g2.Rc[i] / (2 * Rd_)
        gfr = GK * 5e10 / Rd_**2 * 2 * y * y * (i0(y) * k0(y) - i1(y) * k1(y)) / (g2.Rc[i] / Rd_)
        row.append(f"{x} R_d: g_N/Freeman {-F2['gR'][i, 0] / gfr:.3f}, Feldform/Regel {F2['fR'][i, 0] / F2['aR'][i, 0]:.3f}")
    P(f"     R_d = {Rdr} r*:  " + ";  ".join(row))

# ---- Milchstraße
Md, Rd, hz, Mb, ab = 5e10, 3.0, 0.3, 1e10, 0.6
Mtot = Md + Mb
rd = Md / (2 * np.pi * Rd**2) * np.exp(-Rc / Rd) * np.exp(-zc / hz) / (2 * hz)
rd *= (Md / 2) / (rd * gi.vol).sum()
rr = np.maximum(np.hypot(Rc, zc), 1e-4)
rb = Mb * ab / (2 * np.pi * rr * (rr + ab) ** 3)
rb *= (Mb / 2) / (rb * gi.vol).sum()
rho = rd + rb
t0 = time.time()
F = felder(gi, rho, Mtot)
P(f"\nA1 Milchstraße (Scheibe {Md:.0e}, R_d = {Rd} kpc, h = {hz} kpc; Bulge {Mb:.0e}, a = {ab} kpc), r* = "
  f"{r_stern(Mtot * M_sun) / kpc:.2f} kpc, R_d/r* = {Rd / (r_stern(Mtot * M_sun) / kpc):.3f}   ({time.time() - t0:.1f} s)")
P("   Mittelebene (unterste Zelle z = %.3f kpc): Feldform gegen örtliche Regel, Wirbelanteil W_R = a_R − g_F,R" % gi.zc[0])
P("     R/R_d   R [kpc]   a_N        a Regel    a Feldform   Feldform/Regel   W_R/a_Regel   v Regel / Feldform [km/s]")
mid = {}
for x in (0.25, 0.5, 1.0, 2.0, 2.67, 3.5, 5.0, 8.0):
    i = np.searchsorted(gi.Rc, x * Rd)
    aN, aR, fR = -F["gR"][i, 0], -F["aR"][i, 0], -F["fR"][i, 0]
    mid[x] = fR / aR
    P(f"     {x:5.2f}   {gi.Rc[i]:6.2f}   {aN:9.2f}   {aR:9.2f}   {fR:9.2f}     {fR / aR:.3f}           {1 - fR / aR:+.3f}"
      f"        {math.sqrt(aR * gi.Rc[i]):6.1f} / {math.sqrt(max(fR, 0) * gi.Rc[i]):6.1f}")
# Divergenzfreiheit des Wirbelanteils (Probe, Flächenfluss) und Größe in der Schicht
gRf, gZf, aRf, aZf = gi.face_g(F["phiN"])
fRf = -(F["phiF"][1:, :] - F["phiF"][:-1, :]) / np.diff(gi.Rc)[:, None]
fZf = -(F["phiF"][:, 1:] - F["phiF"][:, :-1]) / np.diff(gi.zc)[None, :]
WRf = nuK(aRf) * gRf - fRf; WZf = nuK(aZf) * gZf - fZf
div = np.zeros_like(F["phiN"])
flR = 2 * np.pi * gi.Rf[1:-1, None] * gi.dz[None, :] * WRf; flZ = np.pi * gi.dR2[:, None] * WZf
div[:-1, :] += flR; div[1:, :] -= flR; div[:, :-1] += flZ; div[:, 1:] -= flZ
ref = np.abs(flR).max()
P(f"   Probe: Flächenfluss von W je Zelle, größter Betrag / größter Einzelfluss = {np.abs(div).max() / ref:.1e} (divergenzfrei)")
P("   Wirbelanteil über der Ebene (R = R_d): |W|/|a_Regel| bei z = 0 / h / 3h / 1 kpc / 3 kpc")
i = np.searchsorted(gi.Rc, Rd)
vals = []
for z in (0.0, hz, 3 * hz, 1.0, 3.0):
    j = min(np.searchsorted(gi.zc, z), gi.NZ - 1)
    vals.append(math.hypot(F["WR"][i, j], F["WZ"][i, j]) / math.hypot(F["aR"][i, j], F["aZ"][i, j]))
P("     " + " / ".join(f"{v:.3f}" for v in vals))

# Spannung des Austauschterms: senkrechte Bilanz ∂_zΣ_zz = ρ·W_z, Σ_zz(z) = −∫_z^∞ ρW_z dz'
P("\nA2 Größe des Austauschterms T_Austausch (senkrechte Bilanz ∂_zΣ_zz = ρW_z, radial Σ_RR ~ ∫ρW_R dz·h)")
P("     R/R_d   Σ_zz(0)/P_z(0)   Σ_R/P_z(0)   Σ/(ρ₀c²)      (P_z(0) = ∫₀^∞ ρ|a_z| dz Gewicht der Säule nach der Regel)")
ckm = c / 1e3
for x in (0.25, 1.0, 2.0, 3.5):
    i = np.searchsorted(gi.Rc, x * Rd)
    dz = gi.dz
    Szz = np.sum(rho[i, :] * F["WZ"][i, :] * dz)
    Pz = np.sum(rho[i, :] * np.abs(F["aZ"][i, :]) * dz)
    SR = np.sum(rho[i, :] * np.abs(F["WR"][i, :]) * dz) * hz
    P(f"     {x:5.2f}   {abs(Szz) / Pz:10.3f}       {SR / Pz:8.3f}     {max(abs(Szz), SR) / (rho[i, 0] * ckm**2):.1e}")

# ---- Bahnen: getauschte Arbeit auf 10 Mrd. Jahre
P("\nA3 Arbeit mit dem Medium auf Sternbahnen, 10 Mrd. Jahre: A(t) = ∫W·v dt (W = a_Regel − g_F; mit E = v²/2 + Φ_F gilt dE/dt = W·v)")
P("   säkular = Mittel von A über die letzte minus über die erste Mrd. Jahre (die Schwingung längs der Bahn mittelt sich heraus);")
P("   Kontrolle: dieselben Bahnen im Feld g_F allein, E ebenso über die Fenster gemittelt (Fehler der Interpolation)")
# Felder einmal auf ein gleichmäßiges Feingitter (R 0,5–30 kpc, z 0–6 kpc), dann schnelle lineare Interpolation
from scipy.ndimage import map_coordinates
Ru = np.arange(0.5, 30.0001, 0.02); zu = np.arange(0.0, 6.0001, 0.005)
PU = np.stack(np.meshgrid(Ru, zu, indexing="ij"), -1).reshape(-1, 2)
def fein(A):
    f = RegularGridInterpolator((gi.Rc, gi.zc), A, bounds_error=False, fill_value=None)
    return f(PU).reshape(len(Ru), len(zu))
UF = {k: fein(F[k]) for k in ("aR", "aZ", "fR", "fZ", "phiF", "WR", "WZ")}
def mk(A):
    def f(p):
        return map_coordinates(A, [(p[:, 0] - Ru[0]) / 0.02, np.minimum(p[:, 1], zu[-1]) / 0.005], order=1, mode="nearest")
    return f
IaR, IaZ, IfR, IfZ, IpF, IWR, IWZ = (mk(UF[k]) for k in ("aR", "aZ", "fR", "fZ", "phiF", "WR", "WZ"))
vcR = lambda R: np.sqrt(np.maximum(-R * IaR(np.c_[R, 0 * R]), 0))

def bahnen(R0, vR0, vz0, Lz, regel, T=1e4 / TU, dt=0.5 / TU):
    R, z, vR, vz = R0.copy(), np.zeros_like(R0), vR0.copy(), vz0.copy()
    IR, IZ = (IaR, IaZ) if regel else (IfR, IfZ)
    def acc(R, z):
        p = np.c_[R, np.abs(z)]
        return IR(p) + Lz**2 / R**3, np.sign(z) * IZ(p)
    E = lambda R, z, vR, vz: 0.5 * (vR**2 + vz**2) + 0.5 * Lz**2 / R**2 + IpF(np.c_[R, np.abs(z)])
    aR, aZ = acc(R, z); n = int(T / dt); nw = int(1e3 / TU / dt)
    A = np.zeros_like(R); sA = [np.zeros_like(R), np.zeros_like(R)]; sE = [np.zeros_like(R), np.zeros_like(R)]; zs = [0.0, 0.0]
    for k in range(n):
        vR += 0.5 * dt * aR; vz += 0.5 * dt * aZ
        p = np.c_[R, np.abs(z)]
        w0 = IWR(p) * vR + np.sign(z) * IWZ(p) * vz
        R += dt * vR; z += dt * vz
        p = np.c_[R, np.abs(z)]
        A += 0.5 * dt * (w0 + IWR(p) * vR + np.sign(z) * IWZ(p) * vz)      # Trapezregel über den Driftschritt
        aR, aZ = acc(R, z)
        vR += 0.5 * dt * aR; vz += 0.5 * dt * aZ
        m = 0 if k < nw else (1 if k >= n - nw else -1)
        if m >= 0:
            sA[m] += A / nw; sE[m] += E(R, z, vR, vz) / nw; zs[m] += np.sqrt(np.mean(z**2)) / nw
    return sA[1] - sA[0], sE[1] - sE[0], zs[0], zs[1], A

rng = np.random.default_rng(7)
for (Rlo, Rhi, sz) in ((6.0, 10.0, 18.0), (6.0, 10.0, 30.0), (3.0, 6.0, 30.0)):
    N = 400
    R0 = rng.uniform(Rlo, Rhi, N); Lz = R0 * vcR(R0); vR0 = rng.normal(0, 30.0, N); vz0 = rng.normal(0, sz, N)
    t0 = time.time()
    dA, dE, z0, z1, Aend = bahnen(R0, vR0, vz0, Lz, True)
    _, dEk, zk0, zk1, _ = bahnen(R0, vR0, vz0, Lz, False)
    P(f"   R {Rlo:.0f}–{Rhi:.0f} kpc, σ_z0 = {sz:.0f} km/s: Arbeit säkular Mittel {dA.mean():+6.2f} ± {dA.std() / math.sqrt(N):.2f}, "
      f"Streuung {dA.std():5.2f} (km/s)² = {dA.mean() / sz**2 * 100:+.2f} % / {dA.std() / sz**2 * 100:.2f} % von σ_z0²; "
      f"ΔE {dE.mean():+.2f}; Schwingung längs der Bahn (A am Ende) Streuung {Aend.std() / sz**2 * 100:.1f} %; rms z {z0:.3f} → {z1:.3f} kpc ({(z1 / z0 - 1) * 100:+.1f} %)")
    P(f"        Kontrolle g_F: ΔE Mittel {dEk.mean():+.2f}, Streuung {dEk.std():.2f} (km/s)², rms z {(zk1 / zk0 - 1) * 100:+.1f} %   ({time.time() - t0:.0f} s)")

# ============================================================================ B  Gewicht w in der Wirkung
P("\n" + "=" * 100)
P("B  Gewicht w in der Wirkung: ℒ_𝒮 = 2hBB + (2/ℓ₀²)[(1 − w)·Q_s(q) − 4π⁴/15],  ν_eff = 1 + (1 − w)(ν − 1)")
P("=" * 100)
def Qs(s):
    q = mp.e**(-s)
    Qt = -4 * (s**3 * mp.polylog(1, q) + 3 * s**2 * mp.polylog(2, q) + 6 * s * mp.polylog(3, q) + 6 * mp.polylog(4, q))
    return float(Qt + 4 * mp.pi**4 / 15)
def one_minus_w(z):
    th = 3 * H0_si() * E_hintergrund(z)
    x = th * tau_D(A0 * a0_frueh_relativ(z))
    return mp.e**(-mp.sqrt(x)), x
for z in (1090.0, 0.0):
    omw, x = one_minus_w(z)
    P(f"   z = {z:6.0f}: θτ_D = {x:.4g} (a₀(z) = {a0_frueh_relativ(z):.3f}·a₀),  1 − w = {mp.nstr(omw, 3)}")
P(f"   Probe der Kernfunktion: 1 − gewicht_w(3H₀) = {1 - gewicht_w(3 * H0_si()):.5f}")
P("   Lineare Störungen, ν_eff − 1 = (1 − w)(ν − 1) bei |δg| = 10⁻³ / 10⁻² / 0,1 / 1 a₀:")
for z in (1090.0, 0.0):
    omw, _ = one_minus_w(z)
    P(f"     z = {z:6.0f}: " + " / ".join(mp.nstr(omw * float(nu(f * A0) - 1), 3) for f in (1e-3, 1e-2, 0.1, 1.0)))
Om_s = omega_vakuum()
omw0 = float(one_minus_w(0.0)[0])
P(f"   Faktor vor ganz Q̃ statt vor Q_s: Vakuumteil heute Ω_s = {Om_s:.4f} → (1 − w)·Ω_s = {omw0 * Om_s:.2e} im mitatmenden Kosmos,"
  f" in jeder gebundenen Zone dagegen voll (Vakuum je Zone)")
P("   Q_s(s) = Q̃(s) + 4π⁴/15 (Probe klein: 4s³/3):  " + ",  ".join(f"s = {s}: {Qs(s):.4g} ({4 * s**3 / 3:.4g})" for s in (0.05, 0.171, 1.0)))
dw = lambda x: math.exp(-math.sqrt(x)) / (2 * math.sqrt(x))
P("   ∂w/∂(θτ_D) = e^{−√x}/(2√x) bei x = 10⁻⁴ / 10⁻² / 1 / 17,9: " + " / ".join(f"{dw(x):.3g}" for x in (1e-4, 1e-2, 1.0, 17.9))
  + "  (am Umkehrpunkt θ → 0 unbeschränkt, integrierbar)")

P("\n   Preisach als Hilfsfeld mit einseitiger Nebenbedingung: max Σθ̃_n unter θ̃_n ≤ θ_n, θ̃_{n+1} ≤ θ̃_n (lineares Programm)")
rng = np.random.default_rng(3)
t = np.linspace(0, 6, 120)
th = 1.0 - 0.25 * t + 0.6 * np.sin(2.3 * t) + 0.15 * rng.normal(size=t.size)
n = t.size
A = np.vstack([np.eye(n), (np.eye(n, k=1) - np.eye(n))[:-1]]); b = np.concatenate([th, np.zeros(n - 1)])
res = linprog(-np.ones(n), A_ub=A, b_ub=b, bounds=[(None, None)] * n, method="highs")
rm = np.minimum.accumulate(th)
kap = res.ineqlin.marginals[:n]
P(f"     Lösung gegen laufendes Minimum: größte Abweichung {np.abs(res.x - rm).max():.1e};  "
  f"Komplementarität max|κ_n(θ_n − θ̃_n)| = {np.abs(kap * (th - res.x)).max():.1e};  aktiv an {np.sum(np.abs(kap) > 1e-9)} von {n} Schritten")
P("     Die Lösung ist kausal (θ̃_n hängt nur an θ_{m≤n}): das Extremum mit einseitigen Bedingungen ergibt das Gedächtnis.")

# ============================================================================ C  Übergabe λ_m
P("\n" + "=" * 100)
P("C  Übergabe λ_m: Energie der Mode ξ (Multiplikator-Staub) gegen den Strukturteil")
P("=" * 100)
Mmw = 6e10 * M_sun
rs = r_stern(Mmw); RQ = rs / (-math.log(q_KOSM))
P(f"   Milchstraße 6·10¹⁰ M☉: r* = {rs / kpc:.2f} kpc, R_Q = r*/(−ln q) = {RQ / kpc:.1f} kpc (q = {q_KOSM:.4f})")
P("   gelesene Masse der Regel (Punktmasse) (ν(r) − 1)·M gegen die der Mode Q·M:  r = R_Q / 100 / 200 / 300 kpc:")
P("     " + " / ".join(f"{(float(nu(G * Mmw / r**2)) - 1) / Q_KOSM:.3f}" for r in (RQ, 100 * kpc, 200 * kpc, 300 * kpc)))
rho_c = 3 * H0_si()**2 / (8 * math.pi * G)
for D, s in ((200.0, 0.171), (500.0, 1.0)):
    rb = D * rho_c * OMEGA_B / OMEGA_M
    w_orb, lm_lo, lm_hi = uebergabezeit(D)
    e_xi = (1 / lm_hi) * Q_KOSM * rb * c**2          # kleinste Rate (N* = 2,4)
    t_w = 14e6 * 3.15576e7                            # schnellster Übergang von w (14 Mio. Jahre)
    e_w = A0**2 * Qs(s) / (8 * math.pi * G) / t_w
    P(f"   {D:.0f}-fache kritische Dichte, s = {s}: ρ_b = {rb:.2e} kg/m³, 1/λ_m = {lm_lo / 3.156e16:.2f}–{lm_hi / 3.156e16:.2f} Mrd. J.;"
      f" Rate der Mode λ_m·Qρ_bc² ≥ {e_xi:.2e} W/m³, expliziter w-Term a₀²Q_s/(8πG)·ẇ ≤ {e_w:.2e} W/m³, Verhältnis {e_w / e_xi:.1e}")
P(f"   Aufbau der gelesenen Masse der Regel ohne Nachlauf: ν·r/c bei 50 kpc = {float(nu(G * Mmw / (50 * kpc)**2)) * 50 * kpc / c / 3.156e13:.2f} Mio. J.")

# ============================================================================ Prüfung von Hand
P("\n" + "=" * 100)
P("Unabhängige Prüfung (von Hand, nicht mit diesem Code)")
P("=" * 100)
P("   1) R_Q Milchstraße: GM = 6,674e-11·1,989e30·6e10 = 7,965e30; /1,1e-10 = 7,241e40; √ = 2,691e20 m = 8,72 kpc;")
P("      q = 5,364/6,364 = 0,8429, −ln q = 0,1709; R_Q = 51,0 kpc  (Code: %.1f kpc)" % (RQ / kpc))
P("   2) 1 − w heute: θ = 3·67,4/3,0857e19 = 6,553e-18 s⁻¹; τ_D = c/a₀ = 2,7254e18 s; x = 17,86; √x = 4,226; e^{−4,226} = 0,0146"
  "  (Code: %.4f)" % omw0)
P("   3) Größe des Austauschterms: W ≈ 0,1·v²/R bei R_d: 0,1·(220e3)²/9,26e19 = 5,2e-11 m/s²; mal h = 9,26e18 m: 4,8e8 m²/s²"
  " ≈ (22 km/s)², also Σ ~ ρσ_z², Σ/(ρc²) ≈ 5e-9")

open(__file__.replace(".py", ".txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
