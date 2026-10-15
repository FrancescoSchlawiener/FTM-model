# Radiale Beschränkung der Rückkehr: Welche Lesart des Verlusts Γ trägt Rotationskurven (SPARC) und Vertikaldynamik zugleich?
# ν = 1/(1 − e^{−s}), s = √(Y/a₀); Y ist das Gefälle, das den Verlust Γ = Y/c setzt. Jede Lesart legt zwei Dinge fest:
#   Argument (welches Gefälle setzt Γ) und Anwendung (auf welche Richtung der Bewegung wirkt die Verstärkung).
#   F  Feldform:   ∇²Φ = ∇·[ν(|∇Φ_N|)∇Φ_N]   (Wirkung mit Y = |∇σ|², wirbelfreie Lösung)
#   a  Kernform:   Y = g_R, Anwendung längs e_r (Kugelradial vom Zentrum)
#   b  Tensor:     je Torusrichtung eine Rückkehr, a_i = ν(|g_i|)·g_i
#   c  Torusmittel: eine Rückkehr je Bahntorus, Γ über die Bahn gemittelt (τ_D ≫ Bahnperioden), Anwendung längs g_N
#   d  Bahnregel:  n⊗n mit n = Dσ/|Dσ| am Ort des Sterns, a = ν(|g_N(x)|)·g_N(x), ohne wirbelfreie Ergänzung
# Daten: SPARC (Lelli, McGaugh, Schombert 2016): daten/sparc/Rotmod_LTG/*.dat und SPARC_Lelli2016c.mrt.
# [A] Annahmen: Υ_d = 0,5, Υ_b = 0,7 fest; Gas als dünne Scheibe aus V_gas (Masse an 1,33·M_HI gebunden); Auswahl Q ≤ 2, Neigung ≥ 30°.
# Sonnenort: Σ_bar(|z| < 1,1 kpc) = 47,1 M☉/pc², g_N,R = 1,374·10⁻¹⁰ (Scheibenmodell) bzw. 1,244·10⁻¹⁰ m/s² (Punktmasse), R₀ = 8,2 kpc.
import numpy as np, sys, io, os, time
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.optimize import lsq_linear, minimize_scalar
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(invalid="ignore", divide="ignore")
HIER = os.path.dirname(os.path.abspath(__file__)); DAT = os.path.join(HIER, "..", "daten", "sparc")
out = []
def P(s=""):
    out.append(s); print(s, flush=True)
t00 = time.time()

G = 4.30091e-6                  # kpc (km/s)²/M☉
A10 = 1e-10*3.0857e13           # 1·10⁻¹⁰ m/s² in (km/s)²/kpc
ZPG = 2*np.pi*G*1e6             # 2πG·(1 M☉/pc²) in (km/s)²/kpc
A0 = float(__import__("os").environ.get("FTM_A0", "1.2e-10"))/1e-10*A10
def nu(y): return 1.0/(-np.expm1(-np.sqrt(np.maximum(y, 1e-300))))       # y = |g|/a₀
UG, WG = np.polynomial.legendre.leggauss(48); UG = 0.5*(UG+1); WG = 0.5*WG   # u ∈ [0, 1]
def mbar(x):                                                              # ∫₀¹ √(1 + x²u²) du
    x = np.asarray(x, float); return 0.5*(np.sqrt(1+x*x) + np.where(x > 1e-8, np.arcsinh(x)/np.maximum(x, 1e-300), 1.0))
# Verhältnis zur radialen Form ν(g_R)·g_R in der Schichtnäherung, g_z(u) = 2πGΣ·u, u gleichverteilt (Tracer = Massenschicht)
def f_F(gR, gz, a0): return nu(np.hypot(gR, gz)/a0)/nu(gR/a0)
def f_d(gR, gz, a0):
    gR, gz = np.asarray(gR, float), np.asarray(gz, float)
    return (nu(np.hypot(gR[..., None], gz[..., None]*UG)/a0) @ WG)/nu(gR/a0)
def f_c(gR, gz, a0):
    gR, gz = np.asarray(gR, float), np.asarray(gz, float)
    return nu(gR*mbar(gz/gR)/a0)/nu(gR/a0)

# =====================================================================================
P("0) Einzelwert von Hand: a_N = 0,873 a₀ → s = 0,93434, e^{−s} = 0,39284, ν = 1,6470")
P(f"   im Skript: ν(0,873) = {nu(0.873):.5f};  ∫₀¹√(1 + x²u²)du bei x = 5: von Hand ½(√26 + asinh 5/5) = 2,7807, im Skript {mbar(5.0):.4f}")
P()

# =====================================================================================
P("1) Mittel des Vertikalfelds über die Dichte einer selbstgravitierenden Schicht: g_z = 2πGΣ·(2M(z)/Σ − 1), M(z)/Σ ist unter der")
P("   Dichte gleichverteilt, also |g_z|/2πGΣ = u gleichverteilt auf [0, 1]: ⟨|g_z|⟩ = ½·2πGΣ, ⟨g_z²⟩ = ⅓·(2πGΣ)², für jedes Profil.")
z = np.linspace(0, 40, 400001)
for lab, rho, gzf in (("sech²", 1/np.cosh(z)**2, np.tanh(z)), ("exponentiell", np.exp(-z), 1-np.exp(-z))):
    w = rho/np.trapezoid(rho, z)
    P(f"   {lab:12s}: ⟨|g_z|⟩/2πGΣ = {np.trapezoid(w*gzf, z):.5f},  ⟨g_z²⟩/(2πGΣ)² = {np.trapezoid(w*gzf**2, z):.5f}")
P()

# =====================================================================================
def gitter(N, Lmax, h0):
    lo, hi = 0.1, 40
    for _ in range(200):
        mid = 0.5*(lo+hi); w = Lmax*np.sinh(mid/N)/np.sinh(mid)
        if w > h0: lo = mid
        else: hi = mid
    return Lmax*np.sinh(mid*np.arange(N+1)/N)/np.sinh(mid)

class Loeser:
    """Finite Volumen in (R, z ≥ 0); Φ_N und Feldform ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N]; Rand: Punktmasse."""
    def __init__(s, N, Lmax, h0):
        s.Rf = gitter(N, Lmax, h0); s.zf = s.Rf.copy(); s.NR = s.NZ = N
        s.Rc = 0.5*(s.Rf[1:]+s.Rf[:-1]); s.zc = 0.5*(s.zf[1:]+s.zf[:-1])
        dz = np.diff(s.zf); dR2 = np.diff(s.Rf**2); s.Aring = np.pi*dR2
        s.cR = 2*np.pi*s.Rf[1:-1][:, None]*dz[None, :]/(s.Rc[1:]-s.Rc[:-1])[:, None]
        s.cz = np.pi*dR2[:, None]/(s.zc[1:]-s.zc[:-1])[None, :]
        s.cRb = 2*np.pi*s.Rf[-1]*dz/(s.Rf[-1]-s.Rc[-1]); s.czb = np.pi*dR2/(s.zf[-1]-s.zc[-1])
        NR, NZ = N, N; k = lambda i, j: i + NR*j
        I, J = np.meshgrid(np.arange(NR-1), np.arange(NZ), indexing="ij"); a1, b1, c1 = k(I, J).ravel(), k(I+1, J).ravel(), s.cR.ravel()
        I, J = np.meshgrid(np.arange(NR), np.arange(NZ-1), indexing="ij"); a2, b2, c2 = k(I, J).ravel(), k(I, J+1).ravel(), s.cz.ravel()
        A = np.concatenate([a1, a2]); B = np.concatenate([b1, b2]); C = np.concatenate([c1, c2])
        diag = np.zeros(NR*NZ); np.add.at(diag, A, -C); np.add.at(diag, B, -C)
        diag[k(NR-1, np.arange(NZ))] -= s.cRb; diag[k(np.arange(NR), NZ-1)] -= s.czb
        rows = np.concatenate([A, B, np.arange(NR*NZ)]); cols = np.concatenate([B, A, np.arange(NR*NZ)]); vals = np.concatenate([C, C, diag])
        s.LU = splu(sp.csc_matrix((vals, (rows, cols)), shape=(NR*NZ, NR*NZ)))
        s.rb_out = np.hypot(s.Rf[-1], s.zc); s.rb_top = np.hypot(s.Rc, s.zf[-1])
    def flat(s, F): return F.ravel(order="F")
    def feld(s, x): return x.reshape((s.NR, s.NZ), order="F")
    def poisson(s, Mcell, Mring0, Mtot):
        rhs = 4*np.pi*G*np.asarray(Mcell, float)*np.ones((s.NR, s.NZ)); rhs[:, 0] += 2*np.pi*G*Mring0
        rhs[-1, :] -= s.cRb*(-G*Mtot/s.rb_out); rhs[:, -1] -= s.czb*(-G*Mtot/s.rb_top)
        return s.feld(s.LU.solve(s.flat(rhs)))
    def grads(s, F, dphidz0):
        Rc, zc = s.Rc, s.zc; gR = np.zeros_like(F); gz = np.zeros_like(F)
        gR[1:-1, :] = -(F[2:, :]-F[:-2, :])/(Rc[2:]-Rc[:-2])[:, None]; gR[0, :] = -(F[1, :]-F[0, :])/(Rc[1]-Rc[0]); gR[-1, :] = gR[-2, :]
        gz[:, 1:-1] = -(F[:, 2:]-F[:, :-2])/(zc[2:]-zc[:-2])[None, :]; gz[:, -1] = gz[:, -2]
        gz[:, 0] = -0.5*(dphidz0 + (F[:, 1]-F[:, 0])/(zc[1]-zc[0]))
        return gR, gz
    def feldform(s, FN, dphidz0, Mtot, a0):
        gR, gz = s.grads(FN, dphidz0)
        gn = -(FN[1:, :]-FN[:-1, :])/(s.Rc[1:]-s.Rc[:-1])[:, None]; gt = 0.5*(gz[1:, :]+gz[:-1, :]); nR = nu(np.hypot(gn, gt)/a0)
        gn = -(FN[:, 1:]-FN[:, :-1])/(s.zc[1:]-s.zc[:-1])[None, :]; gt = 0.5*(gR[:, 1:]+gR[:, :-1]); nz = nu(np.hypot(gn, gt)/a0)
        rhs = np.zeros_like(FN)
        fl = nR*s.cR*(FN[1:, :]-FN[:-1, :]); rhs[:-1, :] += fl; rhs[1:, :] -= fl
        fl = nz*s.cz*(FN[:, 1:]-FN[:, :-1]); rhs[:, :-1] += fl; rhs[:, 1:] -= fl
        GM = G*Mtot; rs = np.sqrt(GM/a0); vf2 = np.sqrt(GM*a0)
        fN = lambda r: -GM/r; ft = lambda r: -GM/r - vf2*np.log(-np.expm1(-rs/r)); nb = lambda r: nu(GM/r**2/a0)
        rhs[-1, :] += nb(s.rb_out)*s.cRb*(fN(s.rb_out)-FN[-1, :]) - s.cRb*ft(s.rb_out)
        rhs[:, -1] += nb(s.rb_top)*s.czb*(fN(s.rb_top)-FN[:, -1]) - s.czb*ft(s.rb_top)
        return s.feld(s.LU.solve(s.flat(rhs)))

L = Loeser(440, 3000.0, 0.004)
def spalte(F2, R): return np.array([np.interp(np.log(R), np.log(L.Rc), F2[:, j]) for j in range(L.NZ)])

P("2) Dicke Exponentialscheibe (sech², z₀ = 0,15 R_d), 10¹⁰ M☉, a₀ = 1,2: Verhältnis zur radialen Form ν(g_R)·g_R in der Mittelebene")
P("   In der Mittelebene ist g_z = 0 (Spiegelung), also n = Dσ/|Dσ| = e_R: die örtliche Regel ist dort exakt radial.")
P("   F = Feldform (Löser), ört. = ν(|g_N|)·g_R in der untersten Zelle, d/c = Bahnregel bzw. Torusmittel mit Tracer = Massenschicht")
P("   (2D-Feld des Lösers, Dichtegewicht über z); Schicht: dieselben aus g_z = 2πGΣ·u (u gleichverteilt), F-Näherung ν(√(g_R² + (2πGΣ)²)).")
M = 1e10; rs = np.sqrt(G*M/A0)
P(f"   r* = {rs:.3f} kpc.   {'R_d/r*':>6s} {'R/R_d':>5s} {'x':>6s} | {'F':>6s} {'ört.':>6s} {'d':>6s} {'c':>6s} | Schicht: {'F':>6s} {'d':>6s} {'c':>6s} | lg-Anteil d/F, c/F")
for fac in (0.3, 1.0, 3.0):
    Rd = fac*rs; z0 = 0.15*Rd
    ring = M*((1+L.Rf[:-1]/Rd)*np.exp(-L.Rf[:-1]/Rd)-(1+L.Rf[1:]/Rd)*np.exp(-L.Rf[1:]/Rd))
    wz = 0.5*np.diff(np.tanh(L.zf/z0))
    FN = L.poisson(ring[:, None]*wz[None, :], 0.0, M)
    F = L.feldform(FN, np.zeros(L.NR), M, A0)
    gR2, gz2 = L.grads(FN, np.zeros(L.NR)); gR2, gz2 = np.abs(gR2), np.abs(gz2)
    aF2 = np.abs(L.grads(F, np.zeros(L.NR))[0])
    for xr in (0.25, 0.5, 1.0, 2.0):
        R = xr*Rd; gRz, gzz, aFz = spalte(gR2, R), spalte(gz2, R), spalte(aF2, R)
        alg0 = nu(gRz[0]/A0)*gRz[0]
        rF = aFz[0]/alg0
        rloc = nu(np.hypot(gRz[0], gzz[0])/A0)/nu(gRz[0]/A0)
        g = np.hypot(gRz, gzz); algm = np.sum(wz*nu(gRz/A0)*gRz)
        rd = np.sum(wz*nu(g/A0)*gRz)/algm
        rc = nu(np.sum(wz*g)/wz.sum()/A0)*np.sum(wz*gRz)/wz.sum()/(algm/wz.sum())
        Sig = M*np.exp(-R/Rd)/(2*np.pi*Rd*Rd); gzs = 2*np.pi*G*Sig; x = gzs/gRz[0]
        sF, sd, sc = f_F(gRz[0], gzs, A0), f_d(gRz[0], gzs, A0), f_c(gRz[0], gzs, A0)
        P(f"   {'':>9s}{fac:6.1f} {xr:5.2f} {x:6.2f} | {rF:6.3f} {rloc:6.3f} {rd:6.3f} {rc:6.3f} | {'':8s} {sF:6.3f} {float(sd):6.3f} {float(sc):6.3f} |"
          f" {np.log(rd)/np.log(rF):5.2f} {np.log(rc)/np.log(rF):5.2f}")
P("   Die Feldform senkt auch bei dicker Scheibe, obwohl das Gefälle in der Mittelebene radial ist: Die Tiefe kommt aus der wirbelfreien")
P("   Ergänzung von ν∇Φ_N (die abgesenkte ν über der Schicht wird in die Ebene getragen), nicht aus der Richtung der Rückkehr am Stern.")
P()

# =====================================================================================
tab = {}
for line in open(os.path.join(DAT, "SPARC_Lelli2016c.mrt"), encoding="utf-8"):
    f = line.split()
    if len(f) < 18: continue
    try:
        tab[f[0]] = dict(Inc=float(f[5]), Rd=float(f[11]), MHI=float(f[13]), RHI=float(f[14]), Q=int(f[17]))
    except ValueError:
        continue
rot = {n: np.loadtxt(os.path.join(DAT, "Rotmod_LTG", f"{n}_rotmod.dat"), comments="#", ndmin=2) for n in tab}
sel = [n for n in sorted(tab) if tab[n]["Q"] <= 2 and tab[n]["Inc"] >= 30]
DL = 0.15; KN = np.exp(np.arange(np.log(0.03), np.log(250), DL))
def dach(k, R): return np.maximum(0, 1-np.abs(np.log(R)-np.log(KN[k]))/DL)
GK = []; MK = []
for k in range(len(KN)):
    S = dach(k, L.Rc)*1e6; ring = S*L.Aring
    GK.append(-L.grads(L.poisson(0.0, ring, ring.sum()), 2*np.pi*G*S)[0][:, 0]); MK.append(ring.sum())
GK = np.array(GK); MK = np.array(MK)
def gas(a, t):
    R, Vg = a[:, 0], a[:, 3]
    if np.all(Vg == 0): return np.zeros(L.NR)
    y = Vg*np.abs(Vg); Rmax = 1.3*max(R[-1], t["RHI"]); Mt = 1.33e9*t["MHI"]
    use = np.where((KN > 0.5*R[0]) & (KN < Rmax))[0]; n = len(use)
    A = np.array([R*np.interp(R, L.Rc, GK[k]) for k in use]).T
    sg = 0.02*np.max(np.abs(y)) + 0.02*np.abs(y)
    D2 = np.zeros((n-2, n))
    for i in range(n-2): D2[i, i:i+3] = [1, -2, 1]
    sc = np.median(np.abs(A).max(axis=0)/np.max(sg))
    rA = [A/sg[:, None], 0.1*sc*D2]; rb = [y/sg, np.zeros(n-2)]
    if Mt > 0: rA.append((MK[use]/(0.1*Mt))[None, :]); rb.append([10.0])
    w = lsq_linear(np.vstack(rA), np.concatenate(rb), bounds=(0, np.inf)).x
    return sum(w[j]*dach(k, L.Rc) for j, k in enumerate(use))
PT = {k: [] for k in ("gal", "V", "eV", "gob", "gg", "gd", "gbu", "SB", "Sg")}
for ig, n in enumerate(sel):
    a = rot[n]; R = a[:, 0]; Sgp = np.interp(R, L.Rc, gas(a, tab[n]))
    for key, v in (("gal", np.full(len(R), ig)), ("V", a[:, 1]), ("eV", a[:, 2]), ("gob", a[:, 1]**2/R), ("gg", a[:, 3]*np.abs(a[:, 3])/R),
                   ("gd", a[:, 4]*np.abs(a[:, 4])/R), ("gbu", a[:, 5]**2/R), ("SB", a[:, 6]), ("Sg", Sgp)):
        PT[key].append(v)
for k in PT: PT[k] = np.concatenate(PT[k])
gi = PT["gal"].astype(int)
gb = PT["gg"] + 0.5*PT["gd"] + 0.7*PT["gbu"]
ok = (gb > 0) & (PT["V"] > 0) & (PT["eV"] > 0)
x = ZPG*(0.5*PT["SB"] + PT["Sg"])/gb
lnob = np.log(PT["gob"]); sig = 2*PT["eV"]/PT["V"]; W = 1/sig**2
ugal = np.unique(gi[ok]); IDX = {g: np.where(ok & (gi == g))[0] for g in ugal}
first = np.zeros(len(gb), bool)
for g in ugal: first[IDX[g][:3]] = True

P("3) SPARC: 153 Galaxien, Υ_d = 0,5, Υ_b = 0,7, x = 2πGΣ/g_bar (Σ = Υ_d·SBdisk + Σ_gas). Lesarten in der Schichtnäherung mit g_R = g_bar;")
P("   c und d mit Tracer = Massenschicht (obere Grenze); ein Tracer in der Mittelebene gibt bei c und d genau die radiale Form.")
P(f"   Punkte: {ok.sum()} in {len(ugal)} Galaxien; x ≥ 1: {np.sum(ok & (x >= 1))}, x ≥ 3: {np.sum(ok & (x >= 3))}")
l = np.log10(gb/A0)
dD = (lnob - np.log(nu(gb/A0)*gb))/np.log(10)
gzs = x*gb
dF = np.log10(f_F(gb, gzs, A0)); dd = np.log10(f_d(gb, gzs, A0)); dc = np.log10(f_c(gb, gzs, A0))
def beta(idx, y, dFv):
    X = np.column_stack([np.ones(len(idx)), l[idx], l[idx]**2, dFv[idx]])
    return np.linalg.lstsq(X, y[idx], rcond=None)[0][3]
rng = np.random.default_rng(1); BOOT = [rng.choice(ugal, len(ugal)) for _ in range(1000)]
def bpts(bs, msk): return np.concatenate([IDX[g][msk[IDX[g]]] for g in bs])
P("   β = Anteil der Feldform-Tiefe (Regression Δ = c₀ + c₁l + c₂l² + β·Δ_F; F-Näherung als Δ_F); Daten mit Fehler durch Ziehen von Galaxien:")
for lab, msk in (("alle Punkte", ok), ("ohne die drei innersten je Galaxie", ok & ~first)):
    ii = np.where(msk)[0]
    bb = [beta(bpts(bs, msk), dD, dF) for bs in BOOT]
    P(f"   {lab:36s}: Daten β = {beta(ii, dD, dF):+.2f} ± {np.std(bb):.2f} | F 1,00 | a, b 0 | c {beta(ii, dc, dF):.2f} | d {beta(ii, dd, dF):.2f}"
      f"   (Punkte mit Δ_F < −0,05: {np.sum(msk & (dF < -0.05))})")
ref = ok & (x < 0.3)
def detrend(d):
    c = np.polyfit(l[ref], d[ref], 2); return d - np.polyval(c, l)
eD, eF, ed, ec = detrend(dD), detrend(dF), detrend(dd), detrend(dc)
P("   δ = Abweichung bei gleichem g_bar (Trend aus x < 0,3 abgezogen):   x-Bereich | Punkte | δ Daten (±) | F | c | d")
for lo, hi in ((0, 0.3), (0.3, 1), (1, 3), (3, 1e9)):
    mm = ok & (x >= lo) & (x < hi)
    be = [np.mean(eD[bpts(bs, mm)]) for bs in BOOT[:300]]
    P(f"   {'':>62s}{lo:>4}–{('∞' if hi > 1e8 else hi):<4}| {mm.sum():5d} | {eD[mm].mean():+.3f} ({np.std(be):.3f}) | {eF[mm].mean():+.3f} | {ec[mm].mean():+.3f} | {ed[mm].mean():+.3f}")
def chi(a, fac):
    return np.sum(W[ok]*(lnob[ok] - np.log(nu(gb[ok]/a)*gb[ok]*fac(gb[ok], gzs[ok], a)))**2)
def fit(fac):
    r = minimize_scalar(lambda la: chi(np.exp(la)*A10, fac), bounds=(np.log(0.3), np.log(5)), method="bounded", options={"xatol": 1e-6})
    return np.exp(r.x), r.fun
one = lambda gR, gz, a: np.ones_like(gR)
res = {lab: fit(fc) for lab, fc in (("radial (a, b)", one), ("Feldform-Näherung", f_F), ("c (Tracer = Schicht)", f_c), ("d (Tracer = Schicht)", f_d))}
c0 = res["radial (a, b)"][1]
P("   Ein a₀ für alle Punkte (Fehler allein aus errV):")
for lab, (av, cv) in res.items():
    P(f"      {lab:22s}: bestes a₀ = {av:.3f}·10⁻¹⁰, χ²/N = {cv/ok.sum():.2f}, χ² − χ²(radial) = {cv-c0:+.0f}")
P()

# =====================================================================================
P("4) Sonnenort, Dynamikmasse Σ_dyn(|z| < 1,1 kpc) = a_z/(2πG) mit Σ_bar = 47,1 M☉/pc², z = 1,1 kpc, R₀ = 8,2 kpc (Literaturwert 71 ± 6)")
S1 = 47.1; gz1 = S1*ZPG; z1, R0 = 1.1, 8.2
P("   Vertikalfrequenz in der Ebene: a, Kernform (e_r): ν_z² = ν_z,N² + (ν−1)ω_N²;  F, c, d: ν_z² = ν·ν_z,N² (a_z = ν·g_z, g_z(0) = 0);  b: a_z = ν(|g_z|)·g_z, ∝ √z für z → 0")
for lab, gNR_SI in (("Scheibenmodell", 1.374), ("Punktmasse", 1.244)):
    gNR = gNR_SI*A10; nR = nu(gNR/A0); nT = nu(np.hypot(gNR, gz1)/A0)
    K = S1 + (nR-1)*gNR*z1/R0/ZPG
    P(f"   {lab:14s} (g_N,R = {gNR_SI}·10⁻¹⁰, 2πGΣ = {gz1/A10:.3f}·10⁻¹⁰): a Kern e_r {K:5.1f} | a mit Anwendung längs g_N {nR*S1:5.1f} |"
      f" b {nu(gz1/A0)*S1:5.1f} | c {nT*S1:5.1f}–{nR*S1:5.1f} | d {nT*S1:5.1f} | F (örtlich genähert) {nT*S1:5.1f}")
P()

# =====================================================================================
P("5) Kugeln: Auf der Kreisbahn ist g_N radial und gleich groß; F, a, c, d geben a = ν(g_N)·g_N, b ebenso (g_z = 0 in der Bahnebene).")
for y in (1.0, 0.01):
    g = y*A0; a_ok = nu(y)*g; s_ = f"   b, g_N = {y} a₀: "
    for thd in (45, 30):
        th = np.radians(thd); gRs, gzs_ = g*np.sin(th), g*np.cos(th); aR, az = nu(gRs/A0)*gRs, nu(gzs_/A0)*gzs_
        s_ += f"bei {thd}° zur Achse |a_b|/|a| = {np.hypot(aR, az)/a_ok:.3f}, Richtung {np.degrees(np.arctan2(aR, az)-th):+.1f}° gegen radial;  "
    P(s_ + (f"tief bei 45°: 2^¼ = {2**0.25:.3f}" if y < 0.1 else ""))
P("   Die Tensorregel hängt an den Achsen: Kugeln bleiben nicht kugelförmig.")
P("   c: für exzentrische Bahnen setzt das Bahnmittel statt des Ortswerts das Argument; druckgestützte Kugeln ändern sich (nicht gerechnet).")
P()

# =====================================================================================
P("6) Äußeres Feld, Milchstraße (6·10¹⁰ M☉) bei 50 kpc, g_e = 0,03 a₀: ν frei / senkrecht zur Scheibe / in der Scheibenebene")
gi_ = G*6e10/50**2; ge = 0.03*A0
phi = np.linspace(0, 2*np.pi, 200001)[:-1]
gav = np.mean(np.sqrt(gi_**2 + ge**2 + 2*gi_*ge*np.cos(phi)))
P(f"   g_i = {gi_/A0:.4f} a₀, frei ν = {nu(gi_/A0):.2f}")
P(f"   F, d (örtlicher Betrag):          senkrecht {nu(np.hypot(gi_, ge)/A0):.2f}, gleichgerichtet {nu((gi_+ge)/A0):.2f}")
P(f"   a, b (nur das radiale Gefälle):  senkrecht {nu(gi_/A0):.2f} (ohne Wirkung), gleichgerichtet {nu((gi_+ge)/A0):.2f}")
P(f"   c (Γ über die Kreisbahn gemittelt): senkrecht {nu(np.hypot(gi_, ge)/A0):.2f}, in der Bahnebene ⟨|g|⟩ = {gav/A0:.4f} a₀ → {nu(gav/A0):.2f}")
P()
P(f"(Rechenzeit {time.time()-t00:.0f} s)")
open(os.path.join(HIER, "radiale_beschraenkung.txt"), "w", encoding="utf-8").write("\n".join(out)+"\n")
