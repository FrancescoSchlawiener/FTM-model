# Gewicht des Scheibenfelds in der Regel:  a = ν(√(g_R² + f·(2πGΣ)²)/a₀)·g_R  an den SPARC-Rotationskurven.
# f = 0: radiale Form ν(g_R)·g_R;  f = 1: Feldform der dünnen Scheibe (in der Ebene auf 0,02 dex).  ν = 1/(1 − e^{−s}), s = √(|g|/a₀).
# Daten (gemessen): SPARC (Lelli, McGaugh, Schombert 2016), daten/sparc/Rotmod_LTG/*.dat und SPARC_Lelli2016c.mrt.
# Σ = Υ_d·SBdisk + Σ_gas; Σ_gas aus Vgas zurückgerechnet (dünne Scheibe, Masse an 1,33·M_HI gebunden), Verfahren aus sparc_test.py kopiert.
# [A] Annahmen: Υ_b = 0,7; Υ_d = 0,4 / 0,5 / 0,6; Υ streut je Galaxie um 0,1 dex; die Neigung wirkt nur auf Vobs;
#     Abstand und Neigung streuen je Galaxie mit e_D und e_Inc der Tabelle (Gauß); Bulge kugelförmig, nicht in Σ.
import numpy as np, sys, io, os, time
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.optimize import lsq_linear, least_squares, brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(invalid="ignore", divide="ignore", over="ignore")
HIER = os.path.dirname(os.path.abspath(__file__)); DAT = os.path.join(HIER, "..", "daten", "sparc")
out = []
def P(s=""):
    out.append(s); print(s, flush=True)
def log(s): print(s, file=sys.stderr, flush=True)
t00 = time.time()

G = 4.30091e-6                  # kpc (km/s)²/M☉
A10 = 1e-10*3.0857e13           # 1·10⁻¹⁰ m/s² in (km/s)²/kpc
ZPG = 2*np.pi*G*1e6             # 2πG·(1 M☉/pc²) in (km/s)²/kpc
def nu(y): return 1.0/(-np.expm1(-np.sqrt(np.maximum(y, 1e-300))))   # y = |g|/a₀
YB = 0.7

# =====================================================================================
P("0) Einzelwert von Hand: a_N = 0,873 a₀ → s = 0,93434, e^{−s} = 0,39284, ν = 1,6470, a = 1,4378 a₀")
P(f"   im Skript: ν(0,873) = {nu(0.873):.5f}, a = {0.873*nu(0.873):.5f} a₀")
P()

# ---------- Gitter und Poisson-Löser (aus sparc_test.py, nur für die Gasrückrechnung) ----------
def gitter(N, Lmax, h0):
    lo, hi = 0.1, 40
    for _ in range(200):
        mid = 0.5*(lo+hi); w = Lmax*np.sinh(mid/N)/np.sinh(mid)
        if w > h0: lo = mid
        else: hi = mid
    return Lmax*np.sinh(mid*np.arange(N+1)/N)/np.sinh(mid)

class Loeser:
    def __init__(s, N, Lmax, h0):
        s.Rf = gitter(N, Lmax, h0); s.zf = s.Rf.copy(); s.NR = s.NZ = N
        s.Rc = 0.5*(s.Rf[1:]+s.Rf[:-1]); s.zc = 0.5*(s.zf[1:]+s.zf[:-1])
        dz = np.diff(s.zf); dR2 = np.diff(s.Rf**2)
        s.V = np.pi*dR2[:, None]*dz[None, :]; s.Aring = np.pi*dR2
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

L = Loeser(440, 3000.0, 0.004)

# ---------- Daten ----------
tab = {}
for line in open(os.path.join(DAT, "SPARC_Lelli2016c.mrt"), encoding="utf-8"):
    f = line.split()
    if len(f) < 18: continue
    try:
        tab[f[0]] = dict(T=int(f[1]), D=float(f[2]), eD=float(f[3]), fD=int(f[4]), Inc=float(f[5]), eInc=float(f[6]), L=float(f[7]),
                         Rd=float(f[11]), SB0=float(f[12]), MHI=float(f[13]), RHI=float(f[14]), Q=int(f[17]))
    except ValueError:
        continue
rot = {n: np.loadtxt(os.path.join(DAT, "Rotmod_LTG", f"{n}_rotmod.dat"), comments="#", ndmin=2) for n in tab}
sel = [n for n in sorted(tab) if tab[n]["Q"] <= 2 and tab[n]["Inc"] >= 30]
NG = len(sel)

# ---------- Gas als dünne Scheibe aus Vgas (aus sparc_test.py) ----------
DL = 0.15; KN = np.exp(np.arange(np.log(0.03), np.log(250), DL))
def dach(k, R): return np.maximum(0, 1-np.abs(np.log(R)-np.log(KN[k]))/DL)
GK = []; MK = []
for k in range(len(KN)):
    S = dach(k, L.Rc)*1e6; ring = S*L.Aring
    GK.append(-L.grads(L.poisson(0.0, ring, ring.sum()), 2*np.pi*G*S)[0][:, 0]); MK.append(ring.sum())
GK = np.array(GK); MK = np.array(MK)
def gas(a, t, mit_masse=True):
    R, Vg = a[:, 0], a[:, 3]
    if np.all(Vg == 0): return np.zeros(L.NR), 0.0, np.nan
    y = Vg*np.abs(Vg); Rmax = 1.3*max(R[-1], t["RHI"]); Mt = 1.33e9*t["MHI"]
    use = np.where((KN > 0.5*R[0]) & (KN < Rmax))[0]; n = len(use)
    A = np.array([R*np.interp(R, L.Rc, GK[k]) for k in use]).T
    sg = 0.02*np.max(np.abs(y)) + 0.02*np.abs(y)
    D2 = np.zeros((n-2, n))
    for i in range(n-2): D2[i, i:i+3] = [1, -2, 1]
    sc = np.median(np.abs(A).max(axis=0)/np.max(sg))
    rA = [A/sg[:, None], 0.1*sc*D2]; rb = [y/sg, np.zeros(n-2)]
    if mit_masse and Mt > 0: rA.append((MK[use]/(0.1*Mt))[None, :]); rb.append([10.0])
    w = lsq_linear(np.vstack(rA), np.concatenate(rb), bounds=(0, np.inf)).x
    ym = A@w; vm = np.sign(ym)*np.sqrt(np.abs(ym))
    S = sum(w[j]*dach(k, L.Rc) for j, k in enumerate(use))
    return S, np.sqrt(np.mean((vm-Vg)**2)), ((w*MK[use]).sum()/Mt if Mt > 0 else np.nan)

PT = {k: [] for k in ("gal", "R", "V", "eV", "gg", "gd", "gbu", "SB", "Sg")}
RMS = []
for ig, n in enumerate(sel):
    a = rot[n]; t = tab[n]; R = a[:, 0]
    Sg, rmsg, _ = gas(a, t); RMS.append(rmsg)
    for key, v in (("gal", np.full(len(R), ig)), ("R", R), ("V", a[:, 1]), ("eV", a[:, 2]), ("gg", a[:, 3]*np.abs(a[:, 3])/R),
                   ("gd", a[:, 4]*np.abs(a[:, 4])/R), ("gbu", a[:, 5]**2/R), ("SB", a[:, 6]), ("Sg", np.interp(R, L.Rc, Sg))):
        PT[key].append(v)
for k in PT: PT[k] = np.concatenate(PT[k])
gi = PT["gal"].astype(int)
gb05 = PT["gg"] + 0.5*PT["gd"] + YB*PT["gbu"]
ok = (gb05 > 0) & (PT["V"] > 0) & (PT["eV"] > 0)
x05 = ZPG*(0.5*PT["SB"] + PT["Sg"])/gb05
fDg = np.array([tab[n]["fD"] for n in sel]); Dg = np.array([tab[n]["D"] for n in sel]); eDg = np.array([tab[n]["eD"] for n in sel])
INCg = np.deg2rad([tab[n]["Inc"] for n in sel]); eINCg = np.deg2rad([tab[n]["eInc"] for n in sel])
lnob = np.log(PT["V"]**2/PT["R"]); sig = 2*PT["eV"]/PT["V"]
first = np.zeros(len(gi), bool)
for g in range(NG): first[np.where(ok & (gi == g))[0][:3]] = True
P(f"1) Daten: Auswahl Q ≤ 2, Neigung ≥ 30°: {NG} Galaxien, {ok.sum()} verwendbare Punkte;"
  f" Primärabstand (f_D = 2–5) {np.sum(fDg >= 2)} Galaxien, Hubble-Fluss (f_D = 1) {np.sum(fDg == 1)}")
P(f"   Gas zurückgerechnet: Vgas-Abweichung je Galaxie Median {np.nanmedian(RMS):.2f} km/s; Punkte mit x = 2πGΣ/g_R ≥ 1 (Υ_d = 0,5): {np.sum(ok & (x05 >= 1))}"
  f", ≥ 3: {np.sum(ok & (x05 >= 3))}")
P(f"   e_D/D Median: Primärabstand {np.median((eDg/Dg)[fDg >= 2]):.3f}, Hubble-Fluss {np.median((eDg/Dg)[fDg == 1]):.3f};"
  f"  e_Inc Median {np.rad2deg(np.median(eINCg)):.1f}°")
P()

# ---------- Physikalisch begründete Gewichte ohne Anpassung ----------
P("2) Begründete Werte von f: Vertikalfeld einer symmetrischen Schicht g_z(z) = 2πGΣ·F(z), F = Anteil der Säule innerhalb |z|")
P("   Über die Dichte gemittelt ist F gleichverteilt auf [0, 1]; daher ⟨|g_z|⟩ = πGΣ (f = 1/4) und ⟨g_z²⟩ = (2πGΣ)²/3 (f = 1/3) für jedes Profil.")
zz = np.linspace(-40, 40, 400001); dz_ = zz[1]-zz[0]
for lab, rho in (("gleichförmig", (np.abs(zz) < 1).astype(float)), ("sech²", 1/np.cosh(zz)**2), ("exponentiell", np.exp(-np.abs(zz)))):
    rho = rho/(rho.sum()*dz_); F = np.abs(np.cumsum(rho)*dz_*2 - 1)
    P(f"   {lab:13s}: ⟨F⟩² = {(np.sum(rho*F)*dz_)**2:.4f},  ⟨F²⟩ = {np.sum(rho*F*F)*dz_:.4f}")
for kap in (0.5, 0.25):
    rs_ = 1/np.cosh(zz)**2; rs_ = rs_/(rs_.sum()*dz_); F = np.abs(np.cumsum(rs_)*dz_*2 - 1)
    rg = 1/np.cosh(zz/kap)**2; rg = rg/(rg.sum()*dz_)
    P(f"   Gas in eigener sech²-Schicht mit h_gas = {kap} z₀ im Feld der Sterne: ⟨F⟩² = {(np.sum(rg*F)*dz_)**2:.4f}, ⟨F²⟩ = {np.sum(rg*F*F)*dz_:.4f} [A]")
GLx, GLw = np.polynomial.legendre.leggauss(16); FX = 0.5*(GLx+1); FW = 0.5*GLw
def nu_mittel(gR, Z, a0k): return sum(w*nu(np.sqrt(gR**2 + (xx*Z)**2)/a0k) for xx, w in zip(FX, FW))
def f_eff(gR, Z, a0k):
    ziel = nu_mittel(gR, Z, a0k); r = np.empty(len(gR))
    for i in range(len(gR)):
        r[i] = brentq(lambda f: nu(np.sqrt(gR[i]**2 + f*Z[i]**2)/a0k) - ziel[i], 0, 1)
    return r
s_ = "   Tiefe ln[ν(g_R)/ν(√(g_R² + f(2πGΣ)²))] relativ zu f = 1, f = ¼ / ⅓ bei (g_R/a₀, x): "
for gr_ in (0.3, 1.0):
    for xx_ in (1, 3, 5):
        t1 = np.log(nu(gr_)/nu(gr_*np.sqrt(1+xx_**2)))
        s_ += f"({gr_}, {xx_}) {np.log(nu(gr_)/nu(gr_*np.sqrt(1+xx_**2/4)))/t1:.2f} / {np.log(nu(gr_)/nu(gr_*np.sqrt(1+xx_**2/3)))/t1:.2f}  "
P(s_)
m = ok & (x05 >= 1)
fe = f_eff(gb05[m], ZPG*(0.5*PT["SB"][m] + PT["Sg"][m]), 1.2*A10)
P(f"   Regel über die Dichte gemittelt, ⟨ν(√(g_R² + (2πGΣF)²))⟩: wirksames f an den Punkten mit x ≥ 1 (Υ_d = 0,5, a₀ = 1,2) Median {np.median(fe):.3f}"
  f", 5–95 % {np.percentile(fe, 5):.3f}–{np.percentile(fe, 95):.3f}")
P()

# ---------- Modell und Anpassung ----------
class Satz:
    def __init__(s, Yd, mk):
        s.k = np.where(mk)[0]; k = s.k; s.g = gi[k]
        s.gg = PT["gg"][k]; s.gs = Yd*PT["gd"][k] + YB*PT["gbu"][k]; s.SY = Yd*PT["SB"][k]; s.Sg = PT["Sg"][k]
        s.d = lnob[k]; s.sig = sig[k]; s.inc = INCg[s.g]
def lnM(S, dU, a0, f):
    eU = np.exp(dU); gb = np.maximum(S.gg + eU*S.gs, 1e-3); Z = ZPG*(eU*S.SY + S.Sg); a0k = a0*A10
    if isinstance(f, str): v = nu_mittel(gb, Z, a0k)
    elif f >= 0: v = nu(np.sqrt(gb**2 + f*Z**2)/a0k)
    else: v = nu(gb/a0k)**2/nu(np.sqrt(gb**2 - f*Z**2)/a0k)       # f < 0: Spiegelbild in ln ν, nur als Maß der Richtung
    return np.log(v*gb)
TAU = np.column_stack([eDg/Dg, eINCg, np.full(NG, 0.1*np.log(10))])          # Vorgaben je Galaxie: ln D, Neigung [rad], ln Υ
PLO = np.column_stack([np.full(NG, -3.0), np.deg2rad(5) - INCg, np.full(NG, -0.7)])
PHI = np.column_stack([np.full(NG, 3.0), np.pi/2 - INCg, np.full(NG, 0.7)])
def werte(S, p, a0, f):
    pk = p[S.g]; it_ = S.inc + pk[:, 1]
    m = lnM(S, pk[:, 2], a0, f) + pk[:, 0] + 2*(np.log(np.sin(it_)) - np.log(np.sin(S.inc)))
    r = (S.d - m)/S.sig
    return r, np.bincount(S.g, r*r, minlength=NG) + ((p/TAU)**2).sum(1), pk, it_
def nuis(S, a0, f, it=40):
    """Abstand, Neigung, Υ je Galaxie mit Gauß-Vorgabe mitgeschätzt; drei Startwerte in ln Υ, je Galaxie der beste."""
    best = None
    for u0 in (0.0, -0.35, 0.35):
        c, p = nuis1(S, a0, f, u0, it)
        if best is None: best = (c, p)
        else:
            b_ = c < best[0]; best = (np.where(b_, c, best[0]), np.where(b_[:, None], p, best[1]))
    return best
def nuis1(S, a0, f, u0, it):
    p = np.zeros((NG, 3)); p[:, 2] = u0; lam = np.full(NG, 1e-3)
    r, c, pk, it_ = werte(S, p, a0, f)
    for _ in range(it):
        h = 1e-4; dm = (lnM(S, pk[:, 2]+h, a0, f) - lnM(S, pk[:, 2]-h, a0, f))/(2*h)
        J = np.column_stack([np.ones_like(dm), 2/np.tan(it_), dm])/S.sig[:, None]
        A = np.empty((NG, 3, 3)); b = np.empty((NG, 3))
        for i in range(3):
            b[:, i] = np.bincount(S.g, J[:, i]*r, minlength=NG) - p[:, i]/TAU[:, i]**2
            for j in range(i, 3):
                A[:, i, j] = A[:, j, i] = np.bincount(S.g, J[:, i]*J[:, j], minlength=NG)
        A[:, [0, 1, 2], [0, 1, 2]] += 1/TAU**2
        Ad = A.copy(); Ad[:, [0, 1, 2], [0, 1, 2]] *= (1 + lam[:, None])
        dp = np.linalg.solve(Ad, b[..., None])[..., 0]
        fest = ((p <= PLO + 1e-9) & (dp < 0)) | ((p >= PHI - 1e-9) & (dp > 0))       # am Rand und nach außen: festhalten
        if fest.any():
            Ad[fest[:, :, None] | fest[:, None, :]] = 0; Ad[:, [0, 1, 2], [0, 1, 2]] += fest; b[fest] = 0
            dp = np.linalg.solve(Ad, b[..., None])[..., 0]
        pn = np.clip(p + dp, PLO, PHI)
        _, cn, _, _ = werte(S, pn, a0, f)
        acc = cn <= c
        p = np.where(acc[:, None], pn, p); lam = np.where(acc, lam/3, lam*5)
        r, c, pk, it_ = werte(S, p, a0, f)
    return c, p
def diag_chi(S, a0, f, sI):
    r = (S.d - lnM(S, np.zeros(len(S.d)), a0, f))/sI
    return np.bincount(S.g, r*r, minlength=NG)

# Prüfung des Mitschätzens je Galaxie gegen scipy (Υ_d = 0,4, a₀ = 0,8, f = 0 und Υ_d = 0,5, a₀ = 1,2, f = 0,5)
s_ = "   Prüfung gegen scipy least_squares je Galaxie: "
for Yd, a0t, ft in ((0.4, 0.8, 0.0), (0.5, 1.2, 0.5)):
    S = Satz(Yd, ok); c_, p_ = nuis(S, a0t, ft); dd = []
    for g in range(NG):
        Sg1 = Satz(Yd, ok & (gi == g))
        def res(q):
            pp = np.zeros((NG, 3)); pp[g] = q
            r, _, _, _ = werte(Sg1, pp, a0t, ft); return np.concatenate([r, q/TAU[g]])
        dd.append(c_[g] - 2*least_squares(res, np.zeros(3), bounds=(PLO[g], PHI[g]), x_scale=TAU[g]).cost)
    s_ += (f"(Υ_d {Yd}, a₀ {a0t}, f {ft}) Σχ² {c_.sum():.1f}, größte Abweichung {max(dd):+.3f}, kleinste {min(dd):+.3f};"
           f" ln Υ am Rand ±0,3 dex: {np.sum(np.abs(p_[:, 2]) > 0.699)} Galaxien.  ")
P("3) Anpassung. Modell: ln g_obs = ln[ν(√(g_R² + f·(2πGΣ)²)/a₀)·g_R] + δ_D + 2·ln[sin(i + δ_i)/sin i],  g_R und Σ mit Υ·e^{δ_Υ}")
P("   Zwei Fehlerbehandlungen für Abstand, Neigung und Υ (Vorgabe σ(ln D) = e_D/D und σ(i) = e_Inc aus der Tabelle, σ(ln Υ) = 0,1 dex [A]):")
P("   P  Fehler je Punkt: σ² = (2errV/V)² + (e_D/D)² + (2e_i·cot i)² + (∂ln a/∂ln Υ·0,23)², δ_D = δ_i = δ_Υ = 0;")
P("   N  je Galaxie mitgeschätzt: σ = 2errV/V, δ_D, δ_i, δ_Υ je Galaxie frei mit diesen Vorgaben (|δ_Υ| ≤ 0,3 dex).")
P("   Grenzen in beiden Fällen durch Ziehen ganzer Galaxien, weil die drei Anteile je Galaxie gemeinsam wirken.")
P(s_)

AGR = np.geomspace(0.6, 2.8, 41); LAG = np.log(AGR); DLA = LAG[1]-LAG[0]
FGR = np.array([-1, -0.7, -0.5, -0.35, -0.25, -0.15, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1/3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0])
FLIST = list(FGR) + ["nu"]
def tabelle(Yd, mk, art):
    S = Satz(Yd, mk); T = np.zeros((NG, len(AGR), len(FLIST)))
    if art == "P":
        eU = 0.1*np.log(10); h = 1e-4
        dU = (lnM(S, np.full(len(S.d), h), 1.2, 0.0) - lnM(S, np.full(len(S.d), -h), 1.2, 0.0))/(2*h)
        sI = np.sqrt(S.sig**2 + (eDg/Dg)[S.g]**2 + (2*eINCg[S.g]/np.tan(S.inc))**2 + (dU*eU)**2)
    for ia, a0 in enumerate(AGR):
        for jf, f in enumerate(FLIST):
            if art == "N": T[:, ia, jf] = nuis(S, a0, f)[0]
            elif art == "P": T[:, ia, jf] = diag_chi(S, a0, f, sI)
            else: T[:, ia, jf] = diag_chi(S, a0, f, S.sig)
    return T, len(S.k)

def prof_a(T3):
    """T3[..., ia, jf] → Minimum über a₀ (Parabel in ln a₀) je f: (χ², a₀)."""
    j = np.clip(np.argmin(T3, axis=-2), 1, len(AGR)-2)
    y0 = np.take_along_axis(T3, (j-1)[..., None, :], -2)[..., 0, :]; y1 = np.take_along_axis(T3, j[..., None, :], -2)[..., 0, :]
    y2 = np.take_along_axis(T3, (j+1)[..., None, :], -2)[..., 0, :]
    d = np.clip(0.5*(y0-y2)/(y0-2*y1+y2), -1, 1)
    return y1 - 0.25*(y0-y2)*d, np.exp(LAG[j] + d*DLA)
def best_f(C):
    """C[..., jf] über FGR → bestes f (Parabel durch drei Nachbarn, am Rand der Randwert)."""
    j = np.argmin(C, axis=-1); jj = np.clip(j, 1, len(FGR)-2)
    x0, x1, x2 = FGR[jj-1], FGR[jj], FGR[jj+1]
    y0 = np.take_along_axis(C, (jj-1)[..., None], -1)[..., 0]; y1 = np.take_along_axis(C, jj[..., None], -1)[..., 0]
    y2 = np.take_along_axis(C, (jj+1)[..., None], -1)[..., 0]
    num = (x1-x0)**2*(y1-y2) - (x1-x2)**2*(y1-y0); den = (x1-x0)*(y1-y2) - (x1-x2)*(y1-y0)
    xv = np.clip(x1 - 0.5*num/np.where(den == 0, 1e-30, den), x0, x2)
    return np.where(j == 0, FGR[0], np.where(j == len(FGR)-1, FGR[-1], xv))

GRP = {"Primär": fDg >= 2, "Hubble": fDg == 1, "alle": np.ones(NG, bool)}
rng = np.random.default_rng(7); NB = 1000
BOOT = {k: [np.bincount(rng.choice(np.where(v)[0], v.sum()), minlength=NG) for _ in range(NB)] for k, v in GRP.items()}
NF = len(FGR); JW = [int(np.argmin(abs(FGR - v))) for v in (0, 0.25, 1/3, 1)]
def auswerten(T, npt_mask):
    res = {}
    for k, gm in GRP.items():
        Tt = T[gm].sum(0); C, A0 = prof_a(Tt); Cf, A0f = C[:NF], A0[:NF]
        fb = float(best_f(Cf)); Cmin = Cf.min(); npt = np.sum(npt_mask & gm[gi]); s2 = Cmin/(npt-2)
        W = np.array(BOOT[k]).astype(float); Tb = np.tensordot(W, T, axes=(1, 0)); Cb, Ab = prof_a(Tb)
        fbb = best_f(Cb[:, :NF])
        ab = np.array([np.exp(np.interp(fbb[i], FGR, np.log(Ab[i, :NF]))) for i in range(NB)])
        a_best = float(np.exp(np.interp(fb, FGR, np.log(A0f))))
        # je Galaxie: bestes f beim gemeinsamen a₀ (nächster Gitterpunkt)
        ia = int(np.argmin(abs(LAG - np.log(a_best)))); gidx = np.where(gm)[0]
        fg = FGR[np.argmin(T[gidx, ia, :NF], axis=1)]
        dg = T[gidx, ia, JW[0]] - T[gidx, ia, JW[1]]                 # χ²(f = 0) − χ²(f = ¼) je Galaxie
        top = np.argsort(-np.abs(dg))[:3]
        res[k] = dict(f=fb, q=np.percentile(fbb, [2.5, 16, 50, 84, 97.5]), s2=s2, npt=npt, a=a_best, aq=np.percentile(ab, [16, 84]),
                      aw=[A0f[j] for j in JW], anu=A0[NF], d=np.array([Cf[j]-Cmin for j in JW] + [C[NF]-Cmin])/s2,
                      p0=np.mean(fbb <= 1e-9), da=np.percentile(np.log(Ab[:, JW[3]]/Ab[:, JW[0]]), [16, 50, 84]),
                      fg=np.percentile(fg, [16, 50, 84]), nfg=(np.sum(fg < 0), np.sum(fg == 0), np.sum(fg > 0)),
                      top=", ".join(f"{sel[gidx[i]]} {dg[i]:+.0f}" for i in top), dsum=dg.sum())
    return res
def lab_(key): return f"Υ_d = {key[0]}, {key[1]}, {'je Galaxie mitgeschätzt' if key[2] == 'N' else 'Fehler je Punkt'}"
def zeile(lab, res):
    for k in ("Primär", "Hubble", "alle"):
        r = res[k]; q = r["q"]
        P(f"   {lab:48s} {k:7s} {r['npt']:5d} | f = {r['f']:+.2f}  (16–84 % {q[1]:+.2f} … {q[3]:+.2f}; 2,5–97,5 % {q[0]:+.2f} … {q[4]:+.2f}; f ≤ 0 in {100*r['p0']:.0f} %)"
          f" | χ²/FG {r['s2']:.2f} | Δχ²/(χ²/FG) bei f = 0, ¼, ⅓, 1, ⟨ν⟩: " + " ".join(f"{v:.1f}" for v in r["d"]))
def zeile_g(lab, res):
    for k in ("Primär", "Hubble", "alle"):
        r = res[k]
        P(f"   {lab:48s} {k:7s} | je Galaxie bestes f: 16/50/84 % {r['fg'][0]:+.2f} / {r['fg'][1]:+.2f} / {r['fg'][2]:+.2f}; Galaxien mit f < 0 / = 0 / > 0:"
          f" {r['nfg'][0]} / {r['nfg'][1]} / {r['nfg'][2]} | χ²(0) − χ²(¼) Summe {r['dsum']:+.0f}, größte: {r['top']}")
def zeile_a(lab, res):
    for k in ("Primär", "Hubble", "alle"):
        r = res[k]; a = r["aw"]
        P(f"   {lab:48s} {k:7s} | a₀ bei f = 0: {a[0]:.3f}, ¼: {a[1]:.3f}, ⅓: {a[2]:.3f}, 1: {a[3]:.3f}, ⟨ν⟩: {r['anu']:.3f}"
          f" | beim besten f {r['a']:.3f} (16–84 % {r['aq'][0]:.3f}–{r['aq'][1]:.3f}) | ln[a₀(1)/a₀(0)] = {r['da'][1]:+.3f} ({r['da'][0]:+.3f} … {r['da'][2]:+.3f})")

# Gegenprobe ohne Mitschätzen, Fehler nur errV, Υ_d = 0,5 (f = 0 radial, f = 1 Näherung der Feldform)
T, _ = tabelle(0.5, ok, "E"); RE = auswerten(T, ok)
P(f"   Gegenprobe nur errV, Υ fest 0,5: bestes a₀ bei f = 0 / 1: alle {RE['alle']['aw'][0]:.3f} / {RE['alle']['aw'][3]:.3f};"
  f" Primär {RE['Primär']['aw'][0]:.3f} / {RE['Primär']['aw'][3]:.3f}; Hubble {RE['Hubble']['aw'][0]:.3f} / {RE['Hubble']['aw'][3]:.3f}")
P(f"   nur errV: bestes f alle {RE['alle']['f']:+.2f}, Primär {RE['Primär']['f']:+.2f}, Hubble {RE['Hubble']['f']:+.2f} (χ²/FG {RE['alle']['s2']:.1f})")
P("   f < 0 ist das Spiegelbild in ln ν (die eigene Flächendichte höbe ν) und dient nur als Maß der Richtung.")
P()

ERG = {}; TT = {}
for Yd in (0.4, 0.5, 0.6):
    for mk_lab, mk in (("alle Punkte", ok), ("ohne innerste drei", ok & ~first)):
        for art in ("P", "N"):
            T, _ = tabelle(Yd, mk, art); ERG[(Yd, mk_lab, art)] = auswerten(T, mk); TT[(Yd, mk_lab, art)] = T
            log(f"  {Yd} {mk_lab} {art}  {time.time()-t00:.0f} s")
for art, titel in (("P", "4a) Verlangtes f, Fehler je Punkt: bestes f (± durch Ziehen ganzer Galaxien, 1000-mal) | χ²/FG | Abstand der Werte f = 0, ¼, ⅓, 1, ⟨ν⟩ in χ²/(χ²/FG)"),
                   ("N", "4b) Verlangtes f, Abstand, Neigung, Υ je Galaxie mitgeschätzt (gleiche Spalten)")):
    P(titel)
    for key, r in ERG.items():
        if key[2] == art: zeile(lab_(key), r)
    P()
P("4c) Vorlieben je Galaxie (beim gemeinsamen a₀) und größte Beiträge zu χ²(f = 0) − χ²(f = ¼)")
for key, r in ERG.items():
    zeile_g(lab_(key), r)
P()
P("4d) Ohne die Galaxie mit dem größten Beitrag zu χ²(f = 0) − χ²(f = ¼) (Υ_d = 0,5, alle Punkte)")
for art in ("P", "N"):
    T = TT[(0.5, "alle Punkte", art)]
    for k in ("Primär", "Hubble", "alle"):
        gm = GRP[k].copy(); C, A0 = prof_a(T[gm].sum(0)); ia = int(np.argmin(abs(LAG - np.log(A0[JW[0]]))))
        dg = np.where(gm, T[:, ia, JW[0]] - T[:, ia, JW[1]], 0); g0 = int(np.argmax(np.abs(dg))); gm[g0] = False
        C, A0 = prof_a(T[gm].sum(0)); fb = float(best_f(C[:NF]))
        W = np.array([np.bincount(rng.choice(np.where(gm)[0], gm.sum()), minlength=NG) for _ in range(NB)]).astype(float)
        fbb = best_f(prof_a(np.tensordot(W, T, axes=(1, 0)))[0][:, :NF]); q = np.percentile(fbb, [2.5, 16, 84, 97.5])
        P(f"   {'Fehler je Punkt' if art == 'P' else 'je Galaxie mitgeschätzt':24s} {k:7s} ohne {sel[g0]} ({dg[g0]:+.0f}): f = {fb:+.2f}"
          f" (16–84 % {q[1]:+.2f} … {q[2]:+.2f}; 2,5–97,5 % {q[0]:+.2f} … {q[3]:+.2f}; f ≤ 0 in {100*np.mean(fbb <= 1e-9):.0f} %), a₀ bei f = 0 / ¼ / 1:"
          f" {A0[JW[0]]:.3f} / {A0[JW[1]]:.3f} / {A0[JW[3]]:.3f}")
P()
P("5) a₀ [10⁻¹⁰ m/s²] als Funktion von f, getrennt nach Abstandsmethode")
for key, r in ERG.items():
    zeile_a(lab_(key), r)
P()
P(f"(Rechenzeit {time.time()-t00:.0f} s)")
open(os.path.join(HIER, "gewicht_scheibenfeld.txt"), "w", encoding="utf-8").write("\n".join(out)+"\n")
