"""Warp-Präzession der Milchstraßenscheibe im FTM (offene Kernfrage B7).

Einheiten: kpc, km/s, M☉.  G = 4,30091e-6 kpc (km/s)²/M☉.
Schwelle a0 = 1,2e-10 m/s² (gemessen); Vergleichszeile mit cH0/2π (H0 = 70).

Teile
 1  Baryonisches Massenmodell (McMillan 2017 [G], ohne Halo)
 2  Axialsymmetrischer Finite-Volumen-Poisson-Löser (R,z), Newton + QUMOND-Feldform
    ∇²Φ = ∇²Φ_N + ∇·[(ν−1)∇Φ_N],  ν = 1/(1−q),  q = exp(−√(|∇Φ_N|/a0))
 3  Frequenzen Ω, κ, ν_z; Kontrollen (Rotationskurve, Σ_dyn am Sonnenort)
 4  Präzession eines geneigten Rings / des geneigten Außenteils aus der
    zweiten Ableitung des Energiefunktionals V[ρ] (statischer + mitkippender Teil)
 5  Gedächtnis-Nachlauf 1/λ_m: Frequenzabhängigkeit des mitkippenden Teils
 6  Äußeres Feld (EFE) der Lokalen Gruppe: l=2-Anteil des Feldes, Drehmoment
 7  Große Magellansche Wolke: lineare Antwort der Feldgleichung, l=2-Gezeit
 8  Verwindung a_J, Lense-Thirring, Knotenzahl ζ und Re τ
"""
import sys, io, math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import RectBivariateSpline

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

G = 4.30091e-6                 # kpc (km/s)^2 / Msun
UACC = 1e6 / 3.0857e19         # 1 (km/s)^2/kpc in m/s^2
A0_SI = float(__import__("os").environ.get("FTM_A0", "1.2e-10"))
C_KMS = 299792.458
KMSKPC_PER_GYR = 1.0227        # 1 km/s/kpc = 1,0227 /Gyr
NCHK = [0, 0]


def chk(name, x, y, tol):
    rel = abs(x - y) / max(abs(y), 1e-300)
    ok = rel < tol
    NCHK[0] += 1
    NCHK[1] += (not ok)
    print(f"   [{'OK' if ok else 'ABW'}] {name}: {x:.6g} vs {y:.6g} (rel. {rel:.1e})")


def nu_of(y):
    y = np.maximum(y, 1e-30)
    q = np.exp(-np.sqrt(y))
    return 1.0 / (1.0 - q)


def dnu_dy(y):
    y = np.maximum(y, 1e-30)
    s = np.sqrt(y)
    q = np.exp(-s)
    return -q / (2.0 * s * (1.0 - q) ** 2)


# ---------------------------------------------------------------------------
# 1  Massenmodell (McMillan 2017, MNRAS 465, 76) [G]
# ---------------------------------------------------------------------------
MOD = dict(
    bulge=dict(rho0=98.4e9, r0=0.075, rcut=2.1, alpha=1.8, q=0.5, M=9.23e9),
    thin=dict(S0=896e6, Rd=2.50, zd=0.300),
    thick=dict(S0=183e6, Rd=3.02, zd=0.900),
    HI=dict(S0=53.1e6, Rd=7.0, Rm=4.0, zd=0.085),
    H2=dict(S0=2180e6, Rd=1.5, Rm=12.0, zd=0.045),
)


def rho_components(R, z):
    b = MOD["bulge"]
    rp = np.sqrt(R ** 2 + (z / b["q"]) ** 2)
    rb = b["rho0"] / (1 + rp / b["r0"]) ** b["alpha"] * np.exp(-(rp / b["rcut"]) ** 2)
    out = {"bulge": rb}
    for k in ("thin", "thick"):
        d = MOD[k]
        out[k] = d["S0"] / (2 * d["zd"]) * np.exp(-R / d["Rd"] - np.abs(z) / d["zd"])
    for k in ("HI", "H2"):
        d = MOD[k]
        Rs = np.maximum(R, 1e-6)
        e = np.exp(-np.abs(z) / d["zd"])                      # sech²(z/2z_d) = 4e/(1+e)²
        out[k] = d["S0"] / (4 * d["zd"]) * np.exp(-d["Rm"] / Rs - Rs / d["Rd"]) * 4 * e / (1 + e) ** 2
    return out


# ---------------------------------------------------------------------------
# 2  Gitter und Löser (Halbraum z ≥ 0, gerade Parität)
# ---------------------------------------------------------------------------
def make_faces(d0, xu, growth, xmax):
    f = list(np.arange(0.0, xu + 1e-9, d0))
    d = d0
    while f[-1] < xmax:
        d *= growth
        f.append(f[-1] + d)
    return np.array(f)


Rf = make_faces(0.1, 30.0, 1.05, 1500.0)
zf = np.concatenate([np.arange(0.0, 0.40 + 1e-9, 0.01)])
d = 0.01
while zf[-1] < 1500.0:
    d *= 1.045
    zf = np.append(zf, zf[-1] + d)
Rc = 0.5 * (Rf[1:] + Rf[:-1])
zc = 0.5 * (zf[1:] + zf[:-1])
NR, NZ = len(Rc), len(zc)
dRc = Rf[1:] - Rf[:-1]
dZc = zf[1:] - zf[:-1]
RR, ZZ = np.meshgrid(Rc, zc, indexing="ij")
VOL = np.pi * (Rf[1:] ** 2 - Rf[:-1] ** 2)[:, None] * dZc[None, :]   # Halbraum
idx = np.arange(NR * NZ).reshape(NR, NZ)

# Flächen und Abstände
AR = 2 * np.pi * Rf[1:-1][:, None] * dZc[None, :]              # radiale Innenflächen (NR-1, NZ)
DR = (Rc[1:] - Rc[:-1])[:, None]
AZ = np.pi * (Rf[1:] ** 2 - Rf[:-1] ** 2)[:, None] * np.ones((1, NZ - 1))
DZ = (zc[1:] - zc[:-1])[None, :]


def build_interior_laplacian():
    rows, cols, vals = [], [], []
    wR = (AR / DR).ravel()
    i0 = idx[:-1, :].ravel()
    i1 = idx[1:, :].ravel()
    wZ = (AZ / DZ).ravel()
    j0 = idx[:, :-1].ravel()
    j1 = idx[:, 1:].ravel()
    for a, b, w in ((i0, i1, wR), (j0, j1, wZ)):
        rows += [a, b, a, b]
        cols += [b, a, a, b]
        vals += [w, w, -w, -w]
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    return sp.csr_matrix((vals, (rows, cols)), shape=(NR * NZ, NR * NZ))


L_INT = build_interior_laplacian()

# Dirichlet-Außenränder (Newton): Randwert am Außenrand aus Monopol + Quadrupol
Rout, zout = Rf[-1], zf[-1]
wBR = 2 * np.pi * Rout * dZc / (Rf[-1] - Rc[-1])                   # Rand R = Rout
wBZ = np.pi * (Rf[1:] ** 2 - Rf[:-1] ** 2) / (zf[-1] - zc[-1])     # Rand z = zout
diagB = np.zeros((NR, NZ))
diagB[-1, :] += wBR
diagB[:, -1] += wBZ
L_DIR = (L_INT - sp.diags(diagB.ravel())).tocsc()
LU_DIR = spla.splu(L_DIR)

# Neumann (Phantom): Innenflächen-Laplace, eine Zelle festgehalten
L_NEU = L_INT.tolil()
pin = idx[-1, -1]
L_NEU[pin, :] = 0
L_NEU[pin, pin] = 1.0
LU_NEU = spla.splu(L_NEU.tocsc())


def multipole_boundary(rho):
    """Monopol und Quadrupol Q = ∫ρ r² P2 dV (voller Raum)."""
    m = rho * VOL * 2
    M = m.sum()
    r2 = RR ** 2 + ZZ ** 2
    P2 = 0.5 * (3 * ZZ ** 2 / np.maximum(r2, 1e-30) - 1)
    Q = (m * r2 * P2).sum()
    return M, Q


def phi_multipole(M, Q, R, z):
    r = np.sqrt(R ** 2 + z ** 2)
    P2 = 0.5 * (3 * z ** 2 / r ** 2 - 1)
    return -G * M / r - G * Q * P2 / r ** 3


def solve_newton(rho):
    M, Q = multipole_boundary(rho)
    rhs = 4 * np.pi * G * rho * VOL
    b = np.zeros((NR, NZ))
    b[-1, :] -= wBR * phi_multipole(M, Q, Rout, zc)
    b[:, -1] -= wBZ * phi_multipole(M, Q, Rc, zout)
    phi = LU_DIR.solve((rhs - b).ravel()).reshape(NR, NZ)
    return phi, M


def face_fields(phi):
    """Normalgradienten an Innenflächen und Tangentialkomponente (gemittelt)."""
    gR_c, gz_c = grad_c(phi, +1)
    nR = (phi[1:, :] - phi[:-1, :]) / DR                     # dΦ/dR an R-Flächen
    tR = 0.5 * (gz_c[1:, :] + gz_c[:-1, :])
    nZ = (phi[:, 1:] - phi[:, :-1]) / DZ                     # dΦ/dz an z-Flächen
    tZ = 0.5 * (gR_c[:, 1:] + gR_c[:, :-1])
    return nR, tR, nZ, tZ


def solve_phantom(phiN, a0):
    nR, tR, nZ, tZ = face_fields(phiN)
    yR = np.sqrt(nR ** 2 + tR ** 2) / a0
    yZ = np.sqrt(nZ ** 2 + tZ ** 2) / a0
    FR = (nu_of(yR) - 1) * nR * AR
    FZ = (nu_of(yZ) - 1) * nZ * AZ
    div = np.zeros((NR, NZ))
    div[:-1, :] += FR
    div[1:, :] -= FR
    div[:, :-1] += FZ
    div[:, 1:] -= FZ
    rhs = div.ravel().copy()
    rhs[pin] = 0.0
    return LU_NEU.solve(rhs).reshape(NR, NZ)


def mirror(f, parity):
    return np.concatenate([parity * f[:, ::-1], f], axis=1)


ZE = np.concatenate([-zc[::-1], zc])


def grad_c(f, parity):
    fe = mirror(f, parity)
    gR = np.gradient(fe, Rc, axis=0)[:, NZ:]
    gz = np.gradient(fe, ZE, axis=1)[:, NZ:]
    return gR, gz


# ---------------------------------------------------------------------------
print("=" * 78)
print("1  Massenmodell und Gitter")
print("=" * 78)
comp = rho_components(RR, ZZ)
Mgrid = {k: (v * VOL * 2).sum() for k, v in comp.items()}
comp["bulge"] *= MOD["bulge"]["M"] / Mgrid["bulge"]            # Bulge auf Nennmasse
Mgrid["bulge"] = MOD["bulge"]["M"]
Mana = {"thin": 2 * np.pi * MOD["thin"]["S0"] * MOD["thin"]["Rd"] ** 2,
        "thick": 2 * np.pi * MOD["thick"]["S0"] * MOD["thick"]["Rd"] ** 2}
print(f"   Gitter: {NR} × {NZ} Zellen, R bis {Rf[-1]:.0f} kpc, z bis {zf[-1]:.0f} kpc")
print(f"   kleinste Zellen: ΔR = {dRc[0]:.3f} kpc, Δz = {dZc[0]:.3f} kpc")
for k in ("bulge", "thin", "thick", "HI", "H2"):
    print(f"   {k:6s}: M = {Mgrid[k]:.4e} M☉")
chk("dünne Scheibe, Gittermasse gegen 2πΣ0Rd²", Mgrid["thin"], Mana["thin"], 5e-3)
chk("dicke Scheibe, Gittermasse gegen 2πΣ0Rd²", Mgrid["thick"], Mana["thick"], 5e-3)
rho_b = sum(comp.values())
Mb = sum(Mgrid.values())
print(f"   baryonische Gesamtmasse M_b = {Mb:.4e} M☉")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("2  Newton- und Feldform-Potential (axialsymmetrisch, numerisch)")
print("=" * 78)
# Test: homogene Kugel
Rt = 5.0
rho_t = np.where(RR ** 2 + ZZ ** 2 < Rt ** 2, 1e8, 0.0)
phi_t, Mt = solve_newton(rho_t)
sp_t = RectBivariateSpline(Rc[:300], ZE[NZ - 60:NZ + 60], mirror(phi_t, 1)[:300, NZ - 60:NZ + 60])
for rr in (2.0, 10.0, 20.0):
    gnum = sp_t(rr, 0.0, dx=1)[0, 0]
    Min = Mt * min(1.0, (rr / Rt) ** 3)
    chk(f"Testkugel dΦ/dR bei R = {rr} kpc", gnum, G * Min / rr ** 2, 2e-2)


def field_set(a0_si):
    a0 = a0_si / UACC
    phiN, M = solve_newton(rho_b)
    phiP = solve_phantom(phiN, a0)
    return a0, phiN, phiP


a0, phiN, phiP = field_set(A0_SI)
phiF = phiN + phiP
print(f"   a0 = {A0_SI:.3g} m/s² = {a0:.1f} (km/s)²/kpc")

NRs = np.searchsorted(Rc, 60.0)
NZs = np.searchsorted(zc, 12.0)
Rsub = Rc[:NRs]
Zsub = ZE[NZ - NZs:NZ + NZs]


def spline(phi):
    return RectBivariateSpline(Rsub, Zsub, mirror(phi, 1)[:NRs, NZ - NZs:NZ + NZs], kx=3, ky=3)


spN, spF = spline(phiN), spline(phiF)


def freqs(spl, R):
    PR = spl(R, 0.0, dx=1)[0, 0]
    PRR = spl(R, 0.0, dx=2)[0, 0]
    Pzz = spl(R, 0.0, dy=2)[0, 0]
    Om2 = PR / R
    kap2 = PRR + 3 * Om2
    return math.sqrt(Om2), math.sqrt(kap2), math.sqrt(Pzz), PR


# Feldform-Probe: im Außenbereich der Kugel muss ∇Φ = ν ∇Φ_N gelten (Kugelsymmetrie, curl-frei)
phiPt = solve_phantom(phi_t, a0)
spt2 = RectBivariateSpline(Rc[:300], ZE[NZ - 60:NZ + 60], mirror(phi_t + phiPt, 1)[:300, NZ - 60:NZ + 60])
for rr in (10.0, 20.0):
    gN = G * Mt / rr ** 2
    chk(f"Testkugel Feldform dΦ/dR = ν·g_N bei {rr} kpc", spt2(rr, 0.0, dx=1)[0, 0], nu_of(gN / a0) * gN, 2e-2)

print()
print("   Rotationskurve (Mittelebene)")
print("   R[kpc]  v_N[km/s]  v_FTM[km/s]  g_N/a0   ν=v²/v_N²   Beob. [G]")
obs = {8.2: "229 ± 1 (Eilers 2019) [G]", 12.0: "≈ 222 [G]", 14.0: "≈ 218 [G]", 16.0: "≈ 214 [G]", 20.0: "≈ 209 [G]"}
for R in (4.0, 8.2, 12.0, 14.0, 16.0, 20.0, 25.0):
    vN = math.sqrt(R * freqs(spN, R)[3])
    vF = math.sqrt(R * freqs(spF, R)[3])
    print(f"   {R:5.1f}   {vN:8.1f}   {vF:9.1f}   {freqs(spN, R)[3] / a0:6.3f}   {(vF / vN) ** 2:6.3f}   {obs.get(R, '')}")

# Dynamikmasse am Sonnenort: Σ_dyn(|z|<1,1 kpc) = K_z/(2πG)
Kz_F = spF(8.2, 1.1, dy=1)[0, 0]
Kz_N = spN(8.2, 1.1, dy=1)[0, 0]
Sig_bar = sum(((comp[k][np.searchsorted(Rc, 8.2) - 1, :] * 2 * dZc)[zc < 1.1]).sum() for k in comp)
print(f"\n   Σ(|z|<1,1 kpc) bei R = 8,2 kpc: baryonisch {Sig_bar / 1e6:.1f}, "
      f"Newton {Kz_N / (2 * np.pi * G) / 1e6:.1f}, FTM {Kz_F / (2 * np.pi * G) / 1e6:.1f} M☉/pc²  "
      f"(Beob. 71 ± 6 [G])")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("3  Frequenzen und freie Ringpräzession (Ring im vollen Potential)")
print("=" * 78)
print("   Ω, κ, ν_z in km/s/kpc; ω_p = Ω − ν_z; ω_lin = −(ν_z² − Ω²)/(2Ω)")
print("   Form     R   Ω      κ/Ω    ν_z/Ω   ω_p      ω_lin")
ring_full = {}
for name, spl in (("Newton", spN), ("FTM", spF)):
    for R in (12.0, 14.0, 16.0):
        Om, kap, nz, _ = freqs(spl, R)
        ring_full[(name, R)] = (Om, kap, nz)
        print(f"   {name:7s} {R:4.0f} {Om:6.2f} {kap / Om:6.3f} {nz / Om:6.3f} {Om - nz:8.2f} {-(nz ** 2 - Om ** 2) / (2 * Om):8.2f}")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("4  Präzession des geneigten Außenteils aus dem Energiefunktional")
print("=" * 78)
print("""   V[ρ] = −(1/8πG)∫[|∇Φ_N|² + a0² Q_s(|∇Φ_N|²/a0²)] d³x.  Kippung des Außenteils ρ_o
   um θ:  V = V0 + ½Kθ²,  ω_p = −K/L_o.
   K = K_stat + K_mit
   K_stat = ∫ρ_o ⟨D²(Φ − Φ_N,o)⟩ dV        (Außenteil im ruhenden Gesamtfeld, ohne eigenes Newton-Feld)
   K_mit  = −(1/4πG)∫[(ν−1)⟨|∇ψ|²⟩ + yν'(y)⟨(n̂·∇ψ)²⟩] dV,  ψ = D Φ_N,o  (Feld kippt mit)
   D = −z∂_y + y∂_z (Drehung um x),  ⟨⟩ = Mittel über φ.""")

y_c = None


def cell_fields(phi):
    gR, gz = grad_c(phi, +1)
    return gR, gz


gNR, gNz = cell_fields(phiN)
absgN = np.sqrt(gNR ** 2 + gNz ** 2)
y_c = absgN / a0
nu_c = nu_of(y_c)
ynup = y_c * dnu_dy(y_c)
nR_c = gNR / np.maximum(absgN, 1e-30)
nz_c = gNz / np.maximum(absgN, 1e-30)


def D2avg(psi):
    """⟨D²Ψ⟩_φ für axialsymmetrisches Ψ(R,z)."""
    PR, Pz = grad_c(psi, +1)
    h = RR * Pz - ZZ * PR                               # ungerade in z
    hR = np.gradient(mirror(h / RR, -1), Rc, axis=0)[:, NZ:]
    hz = np.gradient(mirror(h, -1), ZE, axis=1)[:, NZ:]
    return -ZZ * h / RR - 0.5 * ZZ * RR * hR + 0.5 * RR * hz


def h_fields(phi_o):
    PR, Pz = grad_c(phi_o, +1)
    h = RR * Pz - ZZ * PR
    he = mirror(h, -1)
    hR = np.gradient(he, Rc, axis=0)[:, NZ:]
    hz = np.gradient(he, ZE, axis=1)[:, NZ:]
    return h, hR, hz


# Rotationsgeschwindigkeit für den Drehimpuls
vF_c = np.sqrt(np.maximum(RR * grad_c(phiF, +1)[0], 0))
vN_c = np.sqrt(np.maximum(RR * grad_c(phiN, +1)[0], 0))
disc_rho = comp["thin"] + comp["thick"] + comp["HI"] + comp["H2"]


def tilt_K(mask_R, phiF_, phiN_, nu_c_, ynup_, a0_):
    rho_o = disc_rho * mask_R[:, None]
    phi_o, Mo = solve_newton(rho_o)
    Kstat = 2 * (rho_o * D2avg(phiF_ - phi_o) * VOL).sum()
    KN = 2 * (rho_o * D2avg(phiN_ - phi_o) * VOL).sum()
    h, hR, hz = h_fields(phi_o)
    grad2 = 0.5 * (hR ** 2 + hz ** 2 + (h / RR) ** 2)
    ndot2 = 0.5 * (nR_c * hR + nz_c * hz) ** 2
    Kmit = -(2 / (4 * np.pi * G)) * (((nu_c_ - 1) * grad2 + ynup_ * ndot2) * VOL).sum()
    # Kontrolle: (1/4πG)∫⟨|∇ψ|²⟩ = ∫ρ_o⟨D²Φ_o⟩ (Newton-Identität)
    lhs = (2 / (4 * np.pi * G)) * (grad2 * VOL).sum()
    rhs = 2 * (rho_o * D2avg(phi_o) * VOL).sum()
    LoF = 2 * (rho_o * RR * vF_c * VOL).sum()
    LoN = 2 * (rho_o * RR * vN_c * VOL).sum()
    return dict(Mo=Mo, Kstat=Kstat, Kmit=Kmit, KN=KN, LoF=LoF, LoN=LoN, idl=lhs, idr=rhs)


cases = [("Außenteil R > 10", (Rc > 10.0)), ("Außenteil R > 12", (Rc > 12.0)),
         ("Außenteil R > 14", (Rc > 14.0)), ("Außenteil R > 16", (Rc > 16.0)),
         ("Ring 11–13", (Rc > 11) & (Rc < 13)), ("Ring 13–15", (Rc > 13) & (Rc < 15)),
         ("Ring 15–17", (Rc > 15) & (Rc < 17))]
TILT = {}
print()
print("   Teil               M_o[M☉]    ω_N      ω_stat   ω_mit    ω_FTM    (km/s/kpc; <0 = rückläufig)")
for name, m in cases:
    r = tilt_K(m.astype(float), phiF, phiN, nu_c, ynup, a0)
    TILT[name] = r
    wN = -r["KN"] / r["LoN"]
    wS = -r["Kstat"] / r["LoF"]
    wM = -r["Kmit"] / r["LoF"]
    print(f"   {name:18s} {r['Mo']:.3e} {wN:8.2f} {wS:8.2f} {wM:8.2f} {wS + wM:8.2f}")
print()
print("   Kontrolle der Newton-Identität (1/4πG)∫⟨|∇ψ|²⟩ = ∫ρ_o⟨D²Φ_o⟩:")
for name in ("Außenteil R > 12", "Ring 13–15"):
    r = TILT[name]
    chk(f"   {name}", r["idl"], r["idr"], 3e-2)

# Test des Formalismus: dünner Ring bei R_t in Φ_N einer Punktmasse -> ω = 0 (Kugelsymmetrie)
rho_pm = np.where(RR ** 2 + ZZ ** 2 < 1.0, 1e10, 0.0)
phi_pm, _ = solve_newton(rho_pm)
test = 2 * (disc_rho * ((Rc > 13) & (Rc < 15))[:, None] * D2avg(phi_pm) * VOL).sum()
print(f"   Test: Außenring im Feld einer Kugel, K/K_typ = {test / TILT['Ring 13–15']['KN']:.2e} (Soll 0)")

# Vergleich mit a0 = cH0/2π
a0m_si = 299792458.0 * 70e3 / 3.0857e22 / (2 * np.pi)
a0m, phiNm, phiPm = field_set(a0m_si)
gNRm, gNzm = grad_c(phiNm, +1)
ym = np.sqrt(gNRm ** 2 + gNzm ** 2) / a0m
print(f"\n   Vergleich a0 = cH0/2π = {a0m_si:.4g} m/s² (H0 = 70):")
for name in ("Außenteil R > 12", "Ring 13–15"):
    m = [c for c in cases if c[0] == name][0][1].astype(float)
    r = tilt_K(m, phiNm + phiPm, phiNm, nu_of(ym), ym * dnu_dy(ym), a0m)
    print(f"   {name:18s} ω_stat = {-r['Kstat'] / r['LoF']:7.2f}, ω_mit = {-r['Kmit'] / r['LoF']:7.2f}, "
          f"ω_FTM = {-(r['Kstat'] + r['Kmit']) / r['LoF']:7.2f} km/s/kpc")

# Variante leichtere Sternscheiben (Faktor 0,8), Empfindlichkeit auf das Massenmodell
print()
print("   Variante: dünne und dicke Scheibe × 0,8 (Massenmodell-Empfindlichkeit)")
rho_b_save, disc_save = rho_b, disc_rho
rho_b = rho_b_save - 0.2 * (comp["thin"] + comp["thick"])
disc_rho = disc_save - 0.2 * (comp["thin"] + comp["thick"])
phiNv, _ = solve_newton(rho_b)
phiPv = solve_phantom(phiNv, a0)
spFv = spline(phiNv + phiPv)
gRv, gzv = grad_c(phiNv, +1)
yv = np.sqrt(gRv ** 2 + gzv ** 2) / a0
vF_save = vF_c
vF_c = np.sqrt(np.maximum(RR * grad_c(phiNv + phiPv, +1)[0], 0))
nR_save, nz_save = nR_c, nz_c
nR_c = gRv / np.maximum(np.sqrt(gRv ** 2 + gzv ** 2), 1e-30)
nz_c = gzv / np.maximum(np.sqrt(gRv ** 2 + gzv ** 2), 1e-30)
Sig_v = Sig_bar - 0.2 * sum(((comp[k][np.searchsorted(Rc, 8.2) - 1, :] * 2 * dZc)[zc < 1.1]).sum() for k in ("thin", "thick"))
print(f"   M_b = {2 * (rho_b * VOL).sum():.3e} M☉, v_FTM(8,2) = {math.sqrt(8.2 * freqs(spFv, 8.2)[3]):.1f}, "
      f"v_FTM(14) = {math.sqrt(14 * freqs(spFv, 14.0)[3]):.1f} km/s, Σ_bar(1,1) = {Sig_v / 1e6:.1f}, "
      f"Σ_dyn,FTM(1,1) = {spFv(8.2, 1.1, dy=1)[0, 0] / (2 * np.pi * G) / 1e6:.1f} M☉/pc²")
for name in ("Außenteil R > 12", "Ring 13–15"):
    m = [c for c in cases if c[0] == name][0][1].astype(float)
    r = tilt_K(m, phiNv + phiPv, phiNv, nu_of(yv), yv * dnu_dy(yv), a0)
    Om = freqs(spFv, 14.0)[0]
    print(f"   {name:18s} ω_stat = {-r['Kstat'] / r['LoF']:7.2f}, ω_mit = {-r['Kmit'] / r['LoF']:7.2f}, "
          f"ω_FTM = {-(r['Kstat'] + r['Kmit']) / r['LoF']:7.2f} km/s/kpc;  Ring bei 14: Ω − ν_z = {Om - freqs(spFv, 14.0)[2]:.2f}")
rho_b, disc_rho, vF_c, nR_c, nz_c = rho_b_save, disc_save, vF_save, nR_save, nz_save

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("5  Gedächtnis-Nachlauf: der mitkippende Teil folgt mit Kern λ e^{−λΔt}")
print("=" * 78)
print("""   Kippung Θ ∝ e^{st}; der mitkippende Teil folgt mit dem Faktor λ/(λ+s):
   L s = −i[K_stat + K_mit λ/(λ+s)]   ⇒   s² + (λ + i a) s + iλ(a + b) = 0,
   a = K_stat/L, b = K_mit/L.  Präzession = Im s (> 0 rechtläufig), Anwachsen = Re s.""")
print("   Teil               1/λ_m[Myr]   Mode 1: Im s   Re s[1/Gyr]  |  Mode 2: Im s   Re s[1/Gyr]")
MEM = {}
for name in ("Außenteil R > 12", "Ring 13–15"):
    r = TILT[name]
    a_ = r["Kstat"] / r["LoF"]
    b_ = r["Kmit"] / r["LoF"]
    for tl in (0.0, 50.0, 100.0, 170.0, 200.0, 260.0, 500.0, 1e12):
        if tl == 0.0:
            roots = [complex(0, -(a_ + b_)), complex(-np.inf, 0)]
        else:
            lam = 1.0 / (tl / 1000.0) / KMSKPC_PER_GYR           # km/s/kpc
            B = lam + 1j * a_
            Cc = 1j * lam * (a_ + b_)
            disc = np.sqrt(B * B - 4 * Cc)
            roots = sorted([(-B + disc) / 2, (-B - disc) / 2], key=lambda z: -z.real)
        MEM[(name, tl)] = roots
        tls = "∞" if tl > 1e8 else f"{tl:.0f}"
        r1, r2 = roots
        m2 = "   (Nachlaufmode fehlt)" if not np.isfinite(r2.real) else f"{r2.imag:12.2f}  {r2.real * KMSKPC_PER_GYR:11.3f}"
        print(f"   {name:18s} {tls:>8s}   {r1.imag:12.2f}  {r1.real * KMSKPC_PER_GYR:11.3f}  |  {m2}")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("6  Äußeres Feld der Lokalen Gruppe (EFE): l=2-Anteil und Drehmoment")
print("=" * 78)
# sphärisches Hintergrundfeld der Milchstraße aus der Gittermasse
rcell = np.sqrt(RR ** 2 + ZZ ** 2).ravel()
mcell = (rho_b * VOL * 2).ravel()
order = np.argsort(rcell)
r_sorted = rcell[order]
M_cum = np.cumsum(mcell[order])


def M_enc(r):
    return np.interp(r, r_sorted, M_cum, left=0.0, right=M_cum[-1])


def rho_sph(r):
    dr = 1e-3 * r + 1e-4
    return (M_enc(r + dr) - M_enc(r - dr)) / (2 * dr) / (4 * np.pi * r ** 2)


def gN_mw(r):
    return G * M_enc(r) / r ** 2


def yN_from_g(gobs_over_a0):
    lo, hi = 1e-12, 1e3
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if nu_of(mid) * mid < gobs_over_a0:
            lo = mid
        else:
            hi = mid
    return mid


rg = np.logspace(-1, 5, 6001)
mu_n, mu_w = np.polynomial.legendre.leggauss(400)
P2 = 0.5 * (3 * mu_n ** 2 - 1)


def l2_potential_coeff(rho2, rr, Rq):
    """C(R) = Φ_2(R)/R²,  Φ_2(R) = −(4πG/5)[R^{-3}∫_0^R ρ2 r^4 dr + R²∫_R^∞ ρ2/r dr]  (Koeffizient von P2)."""
    out = []
    for R in Rq:
        m1 = rr <= R
        I1 = np.trapezoid(rho2[m1] * rr[m1] ** 4, rr[m1])
        I2 = np.trapezoid(rho2[~m1] / rr[~m1], rr[~m1])
        out.append(-(4 * np.pi * G / 5) * (I1 / R ** 5 + I2))    # C = Φ_2(R)/R²
    return np.array(out)


def efe_l2(ge_obs_over_a0):
    yNe = yN_from_g(ge_obs_over_a0)
    gNe = yNe * a0
    r = rg[:, None]
    mu = mu_n[None, :]
    gr = gN_mw(rg)[:, None]
    dgr = np.gradient(gN_mw(rg), rg)[:, None]
    # g = −g_r r̂ + g_e ê,  |g|² = g_r² + g_e² − 2 g_r g_e μ
    s2 = gr ** 2 + gNe ** 2 - 2 * gr * gNe * mu
    ab = np.sqrt(s2)
    ds2_dr = 2 * gr * dgr - 2 * dgr * gNe * mu
    ds2_dmu = -2 * gr * gNe
    sin2 = 1 - mu ** 2
    g_dot_grads2 = (-gr + gNe * mu) * ds2_dr + (gNe * sin2 / r) * ds2_dmu
    g_dot_gradabs = g_dot_grads2 / (2 * ab)
    y = ab / a0
    # ∇·[(ν−1)g] = (ν−1)∇·g + ν'(y)/a0 (g·∇|g|),  ∇·g = −4πGρ_b ;  4πGρ_ph = −∇·[(ν−1)g]
    rb = rho_sph(rg)[:, None]
    div = (nu_of(y) - 1) * (-4 * np.pi * G * rb) + dnu_dy(y) / a0 * g_dot_gradabs
    rho_ph = -div / (4 * np.pi * G)
    rho2 = 2.5 * (rho_ph * P2[None, :] * mu_w[None, :]).sum(axis=1)
    return gNe, rho2


print("   g_e (beob.)  g_N,e/a0    r_e[kpc]   C(14)[(km/s/kpc)²]  Ω_EFE(12/14/16) max [km/s/kpc]")
EFE = {}
for ge in (0.01, 0.02, 0.03):
    gNe, rho2 = efe_l2(ge)
    Cq = l2_potential_coeff(rho2, rg, [12.0, 14.0, 16.0])
    re = rg[np.argmin(np.abs(gN_mw(rg) - gNe))]
    om = [1.5 * abs(Cq[i]) / ring_full[("FTM", R)][0] for i, R in enumerate((12.0, 14.0, 16.0))]
    EFE[ge] = (Cq, om)
    print(f"   {ge:5.2f} a0    {gNe / a0:.3e}  {re:8.0f}   {Cq[1]:+.3e}        "
          f"{om[0]:.2e} / {om[1]:.2e} / {om[2]:.2e}")
# Kontrolle: ohne äußeres Feld verschwindet der l=2-Anteil
_, rho2_0 = efe_l2(1e-9)
print(f"   Kontrolle g_e → 0: max|ρ2| = {np.abs(rho2_0).max():.2e} M☉/kpc³ (Soll ≈ 0)")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("7  Große Magellansche Wolke: lineare Antwort der Feldgleichung (l = 2)")
print("=" * 78)
M_LMC = 3.2e9          # Sterne 2,7e9 + Gas ~0,5e9 [G]
d_LMC = 50.0           # kpc [G]
b_LMC = -33.0          # galaktische Breite [G]


def linear_response_l2(Mp, dp, ge_obs=0.02, alpha_const=None, l=2):
    """(r²u')' − l(l+1)u = (r²αf')' − l(l+1)(ν−1)f  auf log-Gitter (Finite Volumen)."""
    r = np.logspace(-1.5, 5, 8001)
    lr = np.log(r)
    rf = np.exp(0.5 * (lr[1:] + lr[:-1]))
    rfe = np.concatenate([[r[0] * r[0] / rf[0]], rf, [r[-1] * r[-1] / rf[-1]]])
    dr = rfe[1:] - rfe[:-1]
    yNe = yN_from_g(ge_obs)

    def f_of(x):
        return np.where(x < dp, -G * Mp * x ** l / dp ** (l + 1), -G * Mp * dp ** l / x ** (l + 1))

    def fp_of(x):
        return np.where(x < dp, -l * G * Mp * x ** (l - 1) / dp ** (l + 1), (l + 1) * G * Mp * dp ** l / x ** (l + 2))

    def bg(x):
        yy = np.sqrt(gN_mw(x) ** 2 + (yNe * a0) ** 2) / a0
        return nu_of(yy) - 1, (nu_of(yy) - 1) + yy * dnu_dy(yy)

    num1_c, _ = bg(r)
    num1_f, al_f = bg(rf)
    if alpha_const is not None:
        num1_c = np.full_like(r, alpha_const); al_f = np.full_like(rf, alpha_const)
    wf = rf ** 2 / (r[1:] - r[:-1])
    n = len(r)
    main = -l * (l + 1) * dr
    A = sp.diags([main], [0], shape=(n, n)).tolil()
    for k in range(n - 1):
        A[k, k] -= wf[k]; A[k, k + 1] += wf[k]
        A[k + 1, k + 1] -= wf[k]; A[k + 1, k] += wf[k]
    flux = rf ** 2 * al_f * fp_of(rf)
    rhs = -l * (l + 1) * num1_c * f_of(r) * dr
    rhs[:-1] += flux
    rhs[1:] -= flux
    A[0, :] = 0; A[0, 0] = 1; rhs[0] = 0
    A[-1, :] = 0; A[-1, -1] = 1; rhs[-1] = 0
    u = spla.spsolve(A.tocsr(), rhs)
    return r, u, f_of(r)


# Kontrolle: konstantes α = ν−1 = c  ⇒  u = c·f
r_t, u_t, f_t = linear_response_l2(M_LMC, d_LMC, alpha_const=1.5)
k14 = np.argmin(np.abs(r_t - 14.0))
chk("lineare Antwort mit konstantem ν−1 = 1,5: u/f bei 14 kpc", u_t[k14] / f_t[k14], 1.5, 3e-3)

r_l, u_l, f_l = linear_response_l2(M_LMC, d_LMC, ge_obs=0.02)
cb = abs(math.cos(math.radians(90 - b_LMC)))
sb = abs(math.sin(math.radians(90 - b_LMC)))
print(f"   M_LMC = {M_LMC:.2e} M☉ (baryonisch [G]), d = {d_LMC} kpc, b = {b_LMC}°, |cos β| = {cb:.3f}")
print("   R[kpc]  C_N[(km/s/kpc)²]  C_FTM     Verstärkung  Ω_prec,N   Ω_prec,FTM  Knotenrate θ_w=0,1  [km/s/kpc]")
LMC = {}
for R in (12.0, 14.0, 16.0):
    k = np.argmin(np.abs(r_l - R))
    CN = f_l[k] / r_l[k] ** 2
    CF = (f_l[k] + u_l[k]) / r_l[k] ** 2
    Om = ring_full[("FTM", R)][0]
    OmN = ring_full[("Newton", R)][0]
    wN = 1.5 * abs(CN) * cb / OmN
    wF = 1.5 * abs(CF) * cb / Om
    wnode = 1.5 * abs(CF) * cb * sb / Om / 0.1
    LMC[R] = (CN, CF, wF, wnode)
    print(f"   {R:5.0f}   {CN:+.4e}      {CF:+.4e}  {CF / CN:6.2f}     {wN:.2e}   {wF:.2e}    {wnode:.2e}")
# Welche baryonische LMC-Masse bräuchte die Knotenrate 10,9 km/s/kpc (lineare Skalierung)?
Mreq = M_LMC * 10.9 / LMC[14.0][3]
print(f"   Für eine Knotenrate von 10,9 km/s/kpc (θ_w = 0,1) bei 14 kpc wäre eine baryonische LMC-Masse von "
      f"{Mreq:.2e} M☉ nötig ({Mreq / Mb:.2f} M_b der Milchstraße).")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("8  Verwindung a_J, Lense-Thirring, Knotenzahl ζ, Re τ")
print("=" * 78)
J = 2 * (disc_rho * RR * vF_c * VOL).sum()
Mtot = Mb
aJ = J / (Mtot * C_KMS)                       # kpc
ds = G * Mtot / (3 * C_KMS ** 2)              # kpc
print(f"   J = {J:.3e} M☉ kpc km/s,  a_J = J/(Mc) = {aJ * 1000:.2f} pc,  Δs = {ds * 3.0857e19:.3e} m,  "
      f"χ = a_J/(3Δs) = {aJ / (3 * ds):.0f}")
rstar = math.sqrt(G * Mtot / a0)
print(f"   r* = √(GM/a0) = {rstar:.2f} kpc")
print("   R    ω_LT=2GJ/(c²R³)  ×ν(R)     Im τ_M  Re τ_M=κ/Ω  Ω−κ     ν_z/Ω   ζ_starr  ζ_Ring  ζ_Gaia")
for R, key, rk in ((12.0, "Außenteil R > 10", "Ring 11–13"), (14.0, "Außenteil R > 12", "Ring 13–15"),
                   (16.0, "Außenteil R > 14", "Ring 15–17")):
    wLT = 2 * G * J / (C_KMS ** 2 * R ** 3)
    Om, kap, nz = ring_full[("FTM", R)]
    yR = freqs(spN, R)[3] / a0
    r = TILT[key]; r2 = TILT[rk]
    zeta = -(r["Kstat"] + r["Kmit"]) / r["LoF"] / Om
    zeta2 = -(r2["Kstat"] + r2["Kmit"]) / r2["LoF"] / Om
    print(f"   {R:4.0f}  {wLT:.2e}         {wLT * nu_of(yR):.2e}  {rstar / (2 * np.pi * R):.4f}  {kap / Om:.4f}    "
          f"{Om - kap:6.2f}  {nz / Om:.3f}   {zeta:+.3f}   {zeta2:+.3f}  {10.9 / Om:+.3f}")

print()
print(f"Prüfungen: {NCHK[0]}, Abweichungen über Toleranz: {NCHK[1]}")
