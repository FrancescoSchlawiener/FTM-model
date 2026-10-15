# SPARC-Scheibentest: gemessene Rotationskurven gegen die algebraische Relation g = ν(g_bar/a₀)·g_bar
# und gegen die Feldform ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N];  ν = 1/(1 − e^{−s}), s = √(|g_N|/a₀).
# Daten (gemessen): SPARC, astroweb.cwru.edu/SPARC (Lelli, McGaugh, Schombert 2016):
#   daten/sparc/Rotmod_LTG/*.dat (Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul; Vdisk, Vbul für Υ = 1) und SPARC_Lelli2016c.mrt.
# Löser: Finite Volumen in (R, z ≥ 0), Aufbau wie in schwelle_anpassung.py (dort unverändert), hier in kpc, km/s, M☉;
#   ein Gitter (440 × 440 Zellen bis 3000 kpc, erste Zelle 4 pc) für alle Galaxien.
# [A] Annahmen: Υ_d = 0,5 und Υ_b = 0,7 (3,6 µm), fest; Bulge kugelförmig (aus Vbul); Gas als dünne Scheibe
#   (aus Vgas und M_HI rückgerechnet); Sternscheibe dünn oder sech² mit z₀ = 0,15 R_d.
import numpy as np, sys, io, os, time
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.optimize import lsq_linear, minimize_scalar
from scipy.stats import spearmanr
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(invalid="ignore", divide="ignore")
HIER = os.path.dirname(os.path.abspath(__file__)); DAT = os.path.join(HIER, "..", "daten", "sparc")
out = []
def P(s=""):
    out.append(s); print(s, flush=True)
def log(s): print(s, file=sys.stderr, flush=True)
t00 = time.time()

G = 4.30091e-6                  # kpc (km/s)²/M☉
A10 = 1e-10*3.0857e13           # 1·10⁻¹⁰ m/s² in (km/s)²/kpc
C_KMS = 299792.458
def a_H(H0): return C_KMS*H0/1000/(2*np.pi)/A10          # cH₀/2π in 10⁻¹⁰ m/s²
ZPG = 2*np.pi*G*1e6             # 2πG·(1 M☉/pc²) in (km/s)²/kpc
def nu(y): return 1.0/(-np.expm1(-np.sqrt(np.maximum(y, 1e-300))))   # y = |g|/a₀
A_MESS, A_674, A_73 = 1.2, a_H(67.4), a_H(73.0)

# =====================================================================================
P("0) Einzelwert von Hand: a_N = 0,873 a₀ → s = √0,873 = 0,93434, e^{−s} = 0,39284, ν = 1/(1 − 0,39284) = 1,6470, a = 1,4378 a₀")
P(f"   im Skript: ν(0,873) = {nu(0.873):.5f}, a = {0.873*nu(0.873):.5f} a₀.   cH₀/2π = {A_674:.4f} (H₀ = 67,4) bzw. {A_73:.4f} (H₀ = 73) ·10⁻¹⁰ m/s²;"
  f"  a₀/(2πG) = {A_MESS*A10/ZPG:.1f} M☉/pc² für 1,2·10⁻¹⁰")
P()

# ---------- Gitter und Löser ----------
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
        """∇²Φ_N = 4πGρ; Mcell Masse je Zelle (obere Hälfte), Mring0 Masse der dünnen Scheibe je Ring."""
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
        """∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N]; a0 in (km/s)²/kpc; Rand: exakte Lösung der Punktmasse."""
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
def ebene(F): return -L.grads(F, np.zeros(L.NR))[0][:, 0]          # Radialbeschleunigung in der Ebene (Betrag)
def an(R, g): return np.interp(np.log(R), np.log(L.Rc), g)

P("1) Prüfung des Lösers in kpc, km/s, M☉ (a₀ = 1,2·10⁻¹⁰)")
a0 = A_MESS*A10; M = 1e10; rs = np.sqrt(G*M/a0); b = 0.3*rs
rho = (1+(L.Rc[:, None]**2+L.zc[None, :]**2)/b**2)**-2.5*L.V; Mc = rho*(0.5*M/rho.sum())
FN = L.poisson(Mc, np.zeros(L.NR), M); F = L.feldform(FN, np.zeros(L.NR), M, a0)
gN, gF = ebene(FN), ebene(F)
i = np.argmin(abs(L.Rc-rs)); R = L.Rc[i]; aex = G*M*R/(R*R+b*b)**1.5
P(f"   Plummer-Kugel 10¹⁰ M☉, b = 0,3 r*, bei R = {R/rs:.4f} r*: a_N/a₀ exakt {aex/a0:.4f}, Gitter {gN[i]/a0:.4f};"
  f" Feldform {gF[i]/a0:.4f}, von Hand ν(a_N)·a_N = {nu(aex/a0)*aex/a0:.4f}")
dev = []
for Rt in np.geomspace(0.1, 25, 12):
    i = np.argmin(abs(L.Rc-Rt*rs)); dev.append(abs(gF[i]/(nu(gN[i]/a0)*gN[i])-1))
P(f"   Plummer, Feldform gegen algebraisch 0,1–25 r*: größte Abweichung {max(dev):.1e}")
h = 1.0; M = 1e9
ring = M*(h/np.sqrt(L.Rf[:-1]**2+h*h)-h/np.sqrt(L.Rf[1:]**2+h*h)); Sig = M*h/(2*np.pi*(L.Rc**2+h*h)**1.5)
FN = L.poisson(0.0, ring, M); F = L.feldform(FN, 2*np.pi*G*Sig, M, a0); gF = ebene(F)
s = "   Kuzmin-Scheibe 10⁹ M☉, h = 1 kpc (exakt a = ν(|g_N|)·g_R): Gitter/exakt bei R/h = "
for x in (0.25, 0.5, 1, 2, 5, 8):
    R = x*h; ex = nu(G*M/(R*R+h*h)/a0)*G*M*R/(R*R+h*h)**1.5
    s += f"{x}: {an(np.array([R]), gF)[0]/ex:.4f}  "
P(s)
M = 1e10; rs = np.sqrt(G*M/a0)
for fac in (1.0, 3.0):
    Rd = fac*rs; Sig = M*np.exp(-L.Rc/Rd)/(2*np.pi*Rd*Rd)
    ring = M*((1+L.Rf[:-1]/Rd)*np.exp(-L.Rf[:-1]/Rd)-(1+L.Rf[1:]/Rd)*np.exp(-L.Rf[1:]/Rd))
    FN = L.poisson(0.0, ring, M); F = L.feldform(FN, 2*np.pi*G*Sig, M, a0)
    xx = np.array([0.25, 0.5, 1, 2, 3, 5, 8]); gN = an(xx*Rd, ebene(FN)); gF = an(xx*Rd, ebene(F))
    P(f"   Exponentialscheibe R_d = {fac:.0f} r*, Feldform/algebraisch bei R/R_d = 0,25 … 8: " + " ".join(f"{v:.3f}" for v in gF/(nu(gN/a0)*gN)))
P()

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
P(f"2) Daten: {len(tab)} Galaxien, Auswahl Q ≤ 2 und Neigung ≥ 30°: {len(sel)} Galaxien, {sum(len(rot[n]) for n in sel)} Punkte")

# ---------- Massenmodelle ----------
def sigma_sterne(a, Rd, Rg):
    """Υ = 1: SBdisk an den Messradien, log-linear dazwischen, innen und außen exponentiell mit R_d der Tabelle [L☉/pc²]."""
    R, SB = a[:, 0], a[:, 6]; m = SB > 0; R, lS = R[m], np.log(SB[m])
    o = np.interp(Rg, R, lS); o = np.where(Rg < R[0], lS[0]+(R[0]-Rg)/Rd, o); o = np.where(Rg > R[-1], lS[-1]-(Rg-R[-1])/Rd, o)
    return np.exp(o)
DL = 0.15; KN = np.exp(np.arange(np.log(0.03), np.log(250), DL))          # Knoten der Gasringe (Dach in ln R)
def dach(k, R): return np.maximum(0, 1-np.abs(np.log(R)-np.log(KN[k]))/DL)
GK = []; MK = []
for k in range(len(KN)):
    S = dach(k, L.Rc)*1e6; ring = S*L.Aring
    GK.append(-L.grads(L.poisson(0.0, ring, ring.sum()), 2*np.pi*G*S)[0][:, 0]); MK.append(ring.sum())
GK = np.array(GK); MK = np.array(MK)
def gas(a, t, mit_masse=True):
    """Σ_gas ≥ 0 aus Vgas|Vgas| (dünne Scheibe), glatt, Gesamtmasse weich an 1,33·M_HI (10 %) gebunden."""
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
def bulge(a):
    """Kugel aus Vbul (Υ = 1): M(<r) = Vbul²·r/G; Zellmassen der oberen Hälfte, Gesamtmasse."""
    R, Vb = a[:, 0], a[:, 5]
    if np.all(Vb == 0): return 0.0, 0.0
    Mr = np.maximum.accumulate(Vb**2*R/G)
    r = np.geomspace(1e-3, 3000, 4000)
    Mf = np.where(r < R[0], Mr[0]*(r/R[0])**3, np.interp(r, R, Mr))
    rho = np.maximum(np.gradient(Mf, r), 0)/(4*np.pi*r*r)
    Mc = np.interp(np.hypot(L.Rc[:, None], L.zc[None, :]), r, rho)*L.V
    return Mc*(0.5*Mr[-1]/Mc.sum()), Mr[-1]

AG = np.array([0.25, 0.4, 0.6, 0.8, 0.95, A_674, A_73, 1.2, 1.35, 1.6, 2.0, 2.8, 4.0])      # a₀-Gitter [10⁻¹⁰ m/s²]
AGV = np.array([0.8, 0.95, A_674, 1.2, 1.35, 1.6, 2.0])                                          # kleineres Gitter für die Υ-Varianten
YD, YB, Z0 = 0.5, 0.7, 0.15
PT = {k: [] for k in ("gal", "R", "V", "eV", "gob", "gg", "gd", "gbu", "Ss", "Sg", "gNd", "gNk", "Rd", "Rk", "R4", "R6", "loc12", "loc10")}
GAL = []
for ig, n in enumerate(sel):
    a = rot[n]; t = tab[n]; R = a[:, 0]
    Ss = sigma_sterne(a, t["Rd"], L.Rc)*1e6; rs_ = Ss*L.Aring; Ls = rs_.sum()
    FSd = L.poisson(0.0, rs_, Ls)                                                         # Υ = 1, dünn
    FSk = L.poisson(rs_[:, None]*0.5*np.diff(np.tanh(L.zf/(Z0*t["Rd"])))[None, :], 0.0, Ls)  # Υ = 1, sech²
    Sg, rmsg, mrg = gas(a, t); Sg6 = Sg*1e6; rg = Sg6*L.Aring; Mg = rg.sum()
    FG = L.poisson(0.0, rg, Mg) if Mg > 0 else np.zeros((L.NR, L.NZ))
    Mcb, Mb = bulge(a)
    FB = L.poisson(Mcb, 0.0, Mb) if Mb > 0 else np.zeros((L.NR, L.NZ))
    gd, gg, gbu = a[:, 4]*np.abs(a[:, 4])/R, a[:, 3]*np.abs(a[:, 3])/R, a[:, 5]**2/R
    def modell(yd, dick):
        FNm = yd*(FSk if dick else FSd) + FG + YB*FB
        d0 = 2*np.pi*G*(Sg6 + (0 if dick else yd*Ss)); Mt = yd*Ls + Mg + YB*Mb
        return FNm, d0, Mt
    rat = {}
    for key, yd, dick, grid_ in (("Rd", YD, False, AG), ("Rk", YD, True, AG), ("R4", 0.4, False, AGV), ("R6", 0.6, False, AGV)):
        FNm, d0, Mt = modell(yd, dick); gN = an(R, ebene(FNm))
        if key == "Rd": PT["gNd"].append(gN)
        if key == "Rk": PT["gNk"].append(gN)
        cols = []
        for av in grid_:
            gF = an(R, ebene(L.feldform(FNm, d0, Mt, av*A10)))
            cols.append(gF/(nu(np.abs(gN)/(av*A10))*np.abs(gN)))
        rat[key] = np.array(cols).T
    Sgp = np.interp(R, L.Rc, Sg)
    gNd = PT["gNd"][-1]
    for key, av in (("loc12", A_MESS), ("loc10", A_674)):
        gz = ZPG*(YD*a[:, 6] + Sgp); gr = np.abs(gNd)
        PT[key].append(nu(np.hypot(gr, gz)/(av*A10))/nu(gr/(av*A10)))
    for key, v in (("gal", np.full(len(R), ig)), ("R", R), ("V", a[:, 1]), ("eV", a[:, 2]), ("gob", a[:, 1]**2/R), ("gg", gg), ("gd", gd),
                   ("gbu", gbu), ("Ss", a[:, 6]), ("Sg", Sgp), ("Rd", rat["Rd"]), ("Rk", rat["Rk"]), ("R4", rat["R4"]), ("R6", rat["R6"])):
        PT[key].append(v)
    GAL.append(dict(name=n, rmsg=rmsg, mrg=mrg, Mb=Mb, Ls=Ls, Mg=Mg))
    if ig % 10 == 0: log(f"  {ig}/{len(sel)} {n}  {time.time()-t00:.0f} s")
for k in PT: PT[k] = np.concatenate(PT[k]) if PT[k][0].ndim == 1 else np.vstack(PT[k])
gi = PT["gal"].astype(int)

# ---------- Punkte, Auswahl ----------
def gbar(yd=YD, yb=YB): return PT["gg"] + yd*PT["gd"] + yb*PT["gbu"]
gb = gbar()
ok = (gb > 0) & (PT["gNd"] > 0) & (PT["gNk"] > 0) & (PT["V"] > 0) & (PT["eV"] > 0) & np.all(PT["Rd"] > 0, axis=1) & np.all(PT["Rk"] > 0, axis=1)
P(f"   verwendbar (g_bar > 0, Vobs > 0, errV > 0): {ok.sum()} Punkte in {len(np.unique(gi[ok]))} Galaxien; ohne: {(~ok).sum()} Punkte (meist innen, Gasloch)")
fD = np.array([tab[sel[i]]["fD"] for i in gi]); Inc = np.array([tab[sel[i]]["Inc"] for i in gi]); eInc = np.array([tab[sel[i]]["eInc"] for i in gi])
Qp = np.array([tab[sel[i]]["Q"] for i in gi])
x = ZPG*(YD*PT["Ss"] + PT["Sg"])/gb                      # 2πGΣ/g_bar, Σ = Υ_d·SBdisk + Σ_gas
xs = ZPG*YD*PT["Ss"]/gb                                  # nur Sterne
lnob = np.log(PT["gob"]); sig = 2*PT["eV"]/PT["V"]; W = 1/sig**2

rms = np.array([g["rmsg"] for g in GAL]); mrg = np.array([g["mrg"] for g in GAL])
P("   Massenmodell, geprüft gegen die Daten:")
P(f"   Gas: Vgas aus der rückgerechneten Scheibe, Abweichung (rms je Galaxie) Median {np.nanmedian(rms):.2f}, 90 % unter {np.nanpercentile(rms, 90):.2f},"
  f" größte {np.nanmax(rms):.1f} km/s; Masse/1,33·M_HI Median {np.nanmedian(mrg):.3f}, außerhalb 0,8–1,25: {np.sum((mrg < 0.8) | (mrg > 1.25))} Galaxien")
mr0 = []
for n in sel[::1]:
    _, _, m0 = gas(rot[n], tab[n], mit_masse=False); mr0.append(m0)
P(f"   Gas ohne Bindung an M_HI: Masse/1,33·M_HI Median {np.nanmedian(mr0):.2f} (16–84 %: {np.nanpercentile(mr0, 16):.2f}–{np.nanpercentile(mr0, 84):.2f});"
  " die Gesamtmasse ist aus Vgas allein schlecht bestimmt (Gas jenseits des letzten Punkts, unstetige Vgas-Kurven)")
m = ok & (PT["gd"] > 0)
for key, lab in (("gNd", "dünn"), ("gNk", f"sech², z₀ = {Z0} R_d")):
    r_ = PT[key][ok]/gb[ok]
    P(f"   Newton-Feld des Modells / g_bar der Daten, Sternscheibe {lab}: Median {np.median(r_):.3f}, 16–84 % {np.percentile(r_, 16):.3f}–{np.percentile(r_, 84):.3f}")
P(f"   Anteil des Gases an Σ bei x ≥ 1: Median {np.median((PT['Sg']/(YD*PT['Ss']+PT['Sg']))[ok & (x >= 1)]):.2f}")
lr = np.log10(PT["Rd"][:, list(AG).index(A_MESS)]); ll = np.log10(PT["loc12"])
P(f"   Näherung ν(√(g_R² + (2πGΣ)²))/ν(g_R) gegen den Löser (dünn, 1,2): mittlere Abweichung {np.mean(ll[ok]-lr[ok]):+.3f} dex, 90 % innerhalb {np.percentile(abs(ll[ok]-lr[ok]), 90):.3f} dex")
P()

# ---------- Hilfen ----------
LAG = np.log(AG)
def lnrat(RT, a, grid_=AG):
    la = np.log(np.clip(a, grid_[0], grid_[-1])); lg = np.log(grid_)
    j = int(np.clip(np.searchsorted(lg, la)-1, 0, len(grid_)-2)); t = (la-lg[j])/(lg[j+1]-lg[j])
    return (1-t)*np.log(RT[:, j]) + t*np.log(RT[:, j+1])
def m_alg(gbv, a): return np.log(nu(gbv/(a*A10))*gbv)
def fit(f, lo=0.05, hi=20):
    r = minimize_scalar(lambda la: f(np.exp(la)), bounds=(np.log(lo), np.log(hi)), method="bounded", options={"xatol": 1e-6})
    return np.exp(r.x), r.fun
def chi_alg(msk, a, lo_=None, gbv=None):
    gbv = gb if gbv is None else gbv; lo_ = lnob if lo_ is None else lo_
    return np.sum(W[msk]*(lo_[msk]-m_alg(gbv[msk], a))**2)
def chi_F(msk, a, RT, lo_=None, gbv=None, grid_=AG):
    gbv = gb if gbv is None else gbv; lo_ = lnob if lo_ is None else lo_
    return np.sum(W[msk]*(lo_[msk]-m_alg(gbv[msk], a)-lnrat(RT[msk], a, grid_))**2)
ugal = np.unique(gi[ok])
rng = np.random.default_rng(1)
BOOT = [rng.choice(ugal, len(ugal)) for _ in range(1000)]
def boot_idx(bs, msk):
    return np.concatenate([np.where(msk & (gi == g))[0] for g in bs])
IDX = {g: np.where(ok & (gi == g))[0] for g in ugal}
def boot_pts(bs): return np.concatenate([IDX[g] for g in bs])
def beta(idx, dD, dF, l):
    X = np.column_stack([np.ones(len(idx)), l[idx], l[idx]**2, dF[idx]])
    return np.linalg.lstsq(X, dD[idx], rcond=None)[0][3]
def bet_str(mk, dD, dF, l):
    bb = []
    for bs in BOOT[:500]:
        ii = boot_pts(bs); bb.append(beta(ii[mk[ii]], dD, dF, l))
    return f"β = {beta(np.where(mk)[0], dD, dF, l):+.2f} ± {np.std(bb):.2f} ({np.sum(mk & (dF < -0.05))} Punkte mit Δ_F < −0,05)"

# =====================================================================================
P("3) Vorhersage A: Residuen gegen x = 2πGΣ/g_bar (Σ = Υ_d·SBdisk + Σ_gas; Bulge nicht in Σ)")
P("   Δ = log g_obs − log[ν(g_bar/a₀)·g_bar]; Feldform: Δ_F = log(Feldform/algebraisch) aus dem Löser je Galaxie.")
P("   δ = Δ minus Trend in l = log(g_bar/a₀) (quadratisch, angepasst an die Punkte mit x < 0,3): Abweichung bei gleichem g_bar.")
P("   β: Regression Δ = c₀ + c₁l + c₂l² + β·Δ_F über alle Punkte (Feldform β = 1, algebraisch β = 0); Fehler durch Ziehen von Galaxien.")
BIN = [(0, 0.3), (0.3, 1), (1, 3), (3, 1e9)]
okI = np.where(ok)[0]
for av in (A_MESS, A_674):
    j = list(AG).index(av)
    l = np.log10(gb/(av*A10)); dD = (lnob - m_alg(gb, av))/np.log(10)
    dFd = np.log10(PT["Rd"][:, j]); dFk = np.log10(PT["Rk"][:, j])
    ref = ok & (x < 0.3)
    def detrend(d):
        c = np.polyfit(l[ref], d[ref], 2); return d - np.polyval(c, l)
    eD, eFd, eFk = detrend(dD), detrend(dFd), detrend(dFk)
    P(f"   a₀ = {av:.3f}·10⁻¹⁰:   x-Bereich | Punkte | Δ Daten (±) | Δ_F dünn / dick | δ Daten (±) | δ_F dünn / dick")
    for lo, hi in BIN:
        mm = ok & (x >= lo) & (x < hi)
        bI = [boot_pts(bs) for bs in BOOT[:300]]; bI = [ii[mm[ii]] for ii in bI]
        bD = [np.mean(dD[ii]) for ii in bI]; bE = [np.mean(eD[ii]) for ii in bI]
        P(f"      {lo:>4}–{('∞' if hi > 1e8 else hi):<4} | {mm.sum():5d} | {dD[mm].mean():+.3f} ({np.std(bD):.3f}) | {dFd[mm].mean():+.3f} / {dFk[mm].mean():+.3f}"
          f" | {eD[mm].mean():+.3f} ({np.std(bE):.3f}) | {eFd[mm].mean():+.3f} / {eFk[mm].mean():+.3f}")
    for dF, lab in ((dFd, "dünn"), (dFk, "dick")):
        b0 = beta(okI, dD, dF, l); bb = [beta(boot_pts(bs), dD, dF, l) for bs in BOOT]
        bsyn = beta(okI, dF + 0*dD, dF, l)
        P(f"      β ({lab}) = {b0:+.2f} ± {np.std(bb):.2f}   (Probe: Feldform-Daten ohne Rauschen geben β = {bsyn:.2f})")
    # Nur Punkte, an denen die Feldform deutlich abweicht
    mm = ok & (dFd < -0.05)
    P(f"      Punkte mit Δ_F(dünn) < −0,05 dex: {mm.sum()} in {len(np.unique(gi[mm]))} Galaxien; δ Daten {eD[mm].mean():+.3f}, δ_F dünn {eFd[mm].mean():+.3f}, dick {eFk[mm].mean():+.3f}")
    P("      nach der Vorhersage Δ_F (dünn) geordnet:  Δ_F-Bereich | Punkte | Galaxien | Median x | δ Daten (±) | δ_F dünn / dick")
    for lo, hi in ((-0.02, 1), (-0.05, -0.02), (-0.1, -0.05), (-0.2, -0.1), (-9, -0.2)):
        mm = ok & (dFd > lo) & (dFd <= hi)
        bI = [boot_pts(bs) for bs in BOOT[:300]]; bE = [np.mean(eD[ii[mm[ii]]]) for ii in bI]
        P(f"      {lo:+.2f} … {hi:+.2f} | {mm.sum():5d} | {len(np.unique(gi[mm])):4d} | {np.median(x[mm]):5.2f} | {eD[mm].mean():+.3f} ({np.std(bE):.3f}) | {eFd[mm].mean():+.3f} / {eFk[mm].mean():+.3f}")
P()

# =====================================================================================
P("4) Vorhersage B: a₀ je Galaxie bei festem Υ (algebraische Relation), gegen Σ₀ = Υ_d·SBdisk(zentral, Tabelle)")
jM = list(AG).index(A_MESS)
rows = []
for g in ugal:
    idx = IDX[g]
    if len(idx) < 5: continue
    mk = np.zeros(len(gb), bool); mk[idx] = True
    aA, _ = fit(lambda a: chi_alg(mk, a))
    aFd, _ = fit(lambda a: chi_F(mk, a, PT["Rd"])); aFk, _ = fit(lambda a: chi_F(mk, a, PT["Rk"]))
    pr = []
    for RT in (PT["Rd"], PT["Rk"]):
        for av in (A_MESS, A_674):
            syn = m_alg(gb, av) + lnrat(RT, av)
            pr.append(fit(lambda a: chi_alg(mk, a, lo_=syn))[0]/av)
    rows.append((g, 0.5*tab[sel[g]]["SB0"], aA, aFd, aFk, *pr))
rows = np.array(rows); S0 = rows[:, 1]; lS = np.log10(S0)
def steig(y, ii=None):
    ii = np.arange(len(y)) if ii is None else ii
    return np.polyfit(lS[ii], np.log10(y[ii]), 1)[0]
bi = [rng.choice(len(S0), len(S0)) for _ in range(1000)]
P(f"   {len(rows)} Galaxien mit ≥ 5 Punkten. Steigung d log a₀/d log Σ₀ (± durch Ziehen), Spearman ρ, Median a₀ je Σ₀-Bereich [M☉/pc²]:")
SB = [(0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e9)]
P("                                     Steigung        ρ     |  <30   30–100  100–300  300–1000  >1000")
P("   Anzahl Galaxien                                          |  " + "  ".join(f"{np.sum((S0 >= lo) & (S0 < hi)):5d}  " for lo, hi in SB))
for col, lab in ((2, "Daten, algebraisch angepasst"), (3, "Daten, Feldform dünn angepasst"), (4, "Daten, Feldform dick angepasst"),
                 (5, "Vorhersage Feldform dünn, a₀ 1,2 (÷a₀)"), (6, "Vorhersage Feldform dünn, cH₀/2π (÷a₀)"),
                 (7, "Vorhersage Feldform dick, a₀ 1,2 (÷a₀)"), (8, "Vorhersage Feldform dick, cH₀/2π (÷a₀)")):
    y = rows[:, col]; st = steig(y); sb = np.std([steig(y, ii) for ii in bi]); rho_ = spearmanr(S0, y)[0]
    P(f"   {lab:40s} {st:+.3f} ± {sb:.3f}  {rho_:+.2f}  | " + "  ".join(f"{np.median(y[(S0 >= lo) & (S0 < hi)]):6.3f}" for lo, hi in SB))
P("   algebraische Regel als Vorhersage: Steigung 0, a₀ je Galaxie = wahres a₀.")
dsl = [steig(rows[:, 2], ii)-steig(rows[:, 5], ii) for ii in bi]
P(f"   Steigung Daten (algebraisch angepasst) minus Vorhersage Feldform dünn (1,2): {steig(rows[:, 2])-steig(rows[:, 5]):+.3f} ± {np.std(dsl):.3f}")
fDg = np.array([tab[sel[int(g)]]["fD"] for g in rows[:, 0]]); i2 = np.where(fDg >= 2)[0]; i1 = np.where(fDg == 1)[0]
P(f"   nach Abstandsmethode: f_D ≥ 2 ({len(i2)} Galaxien) Steigung Daten {steig(rows[:, 2], i2):+.3f}, Vorhersage {steig(rows[:, 5], i2):+.3f};"
  f"  f_D = 1 ({len(i1)}) Daten {steig(rows[:, 2], i1):+.3f}, Vorhersage {steig(rows[:, 5], i1):+.3f};  Median Σ₀ f_D ≥ 2 / = 1: {np.median(S0[i2]):.0f} / {np.median(S0[i1]):.0f} M☉/pc²")
P(f"   Anpassungen am Rand (0,05 oder 20): {np.sum((rows[:, 2:5] < 0.051) | (rows[:, 2:5] > 19.9))}")
P(f"   Streuung der Daten: log a₀ (algebraisch) je Galaxie σ = {np.std(np.log10(rows[:, 2])):.2f} dex")
P()

# =====================================================================================
P("5) Gesamt: ein a₀ für alle Punkte, Υ fest; χ² mit Fehlern allein aus errV")
N = ok.sum()
aA, cA = fit(lambda a: chi_alg(ok, a)); aFd, cFd = fit(lambda a: chi_F(ok, a, PT["Rd"])); aFk, cFk = fit(lambda a: chi_F(ok, a, PT["Rk"]))
P(f"   {N} Punkte, {len(ugal)} Galaxien.          bestes a₀ [10⁻¹⁰] | χ²/N beim besten | χ²/N bei 1,2 | bei cH₀/2π (67,4) | bei cH₀/2π (73)")
for lab, av, cv, fn in (("algebraisch", aA, cA, lambda a: chi_alg(ok, a)), ("Feldform dünn", aFd, cFd, lambda a: chi_F(ok, a, PT["Rd"])),
                        ("Feldform dick", aFk, cFk, lambda a: chi_F(ok, a, PT["Rk"]))):
    P(f"   {lab:15s}                    {av:.3f}        |     {cv/N:.2f}       |    {fn(A_MESS)/N:.2f}    |     {fn(A_674)/N:.2f}         |    {fn(A_73)/N:.2f}")
# wo trennen sie sich
rA = W*(lnob-m_alg(gb, aA))**2
for RT, aF, lab in ((PT["Rd"], aFd, "dünn"), (PT["Rk"], aFk, "dick")):
    rF = W*(lnob-m_alg(gb, aF)-lnrat(RT, aF))**2
    s = f"   χ²(alg) − χ²(Feldform {lab}) je x-Bereich, jeweils beim besten a₀:"
    for lo, hi in BIN:
        mm = ok & (x >= lo) & (x < hi); s += f"  {lo}–{('∞' if hi > 1e8 else hi)}: {rA[mm].sum()-rF[mm].sum():+.0f}"
    P(s)
    dg = np.array([rA[IDX[g]].sum()-rF[IDX[g]].sum() for g in ugal])
    top = np.argsort(-np.abs(dg))[:4]
    P(f"      je Galaxie: {np.sum(dg > 0)} bevorzugen die Feldform, {np.sum(dg < 0)} die algebraische Relation; größte |Δχ²|: "
      + ", ".join(f"{sel[ugal[i]]} {dg[i]:+.0f}" for i in top))
P()

# =====================================================================================
P("6) Grenzen, gerechnet")
def glob3(msk=ok, lo_=None, gbv=None, RT4=None, grid_=AG):
    r = [fit(lambda a: chi_alg(msk, a, lo_, gbv))[0]]
    for RT in ((PT["Rd"], PT["Rk"]) if RT4 is None else (RT4,)):
        r.append(fit(lambda a: chi_F(msk, a, RT, lo_, gbv, grid_))[0])
    return r
for yd, key in ((0.4, "R4"), (0.6, "R6")):
    gbv = gbar(yd); mk = ok & (gbv > 0) & np.all(PT[key] > 0, axis=1)
    r = glob3(mk, gbv=gbv, RT4=PT[key], grid_=AGV)
    j = list(AGV).index(1.2); l = np.log10(gbv/(1.2*A10)); dD = (lnob-m_alg(gbv, 1.2))/np.log(10); dF = np.log10(PT[key][:, j])
    P(f"   Υ_d = {yd} (Υ_b skaliert nicht): bestes a₀ algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f};  {bet_str(mk, dD, dF, l)}")
fh = np.where(fD == 1, 67.4/73.0, 1.0); lo2 = lnob + np.log(fh)
r = glob3(ok, lo_=lo2)
P(f"   Abstände aus dem Hubble-Fluss (f_D = 1, H₀ = 73 in den Daten) auf H₀ = 67,4 umgerechnet (g_obs ∝ 1/D, g_bar und x fest):"
  f" bestes a₀ algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f}")
mk = ok & (fD >= 2); r = glob3(mk)
P(f"   nur Abstände aus TRGB, Cepheiden, UMa, SN (f_D ≥ 2; {len(np.unique(gi[mk]))} Galaxien): algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f}")
mk = ok & (fD == 1); r = glob3(mk)
P(f"   nur Hubble-Fluss (f_D = 1; {len(np.unique(gi[mk]))} Galaxien): algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f}")
sI = np.hypot(sig, 2*np.deg2rad(eInc)/np.tan(np.deg2rad(Inc))); Wsave = W.copy(); W = 1/sI**2
aA2, cA2 = fit(lambda a: chi_alg(ok, a)); aF2, cF2 = fit(lambda a: chi_F(ok, a, PT["Rd"]))
P(f"   mit Neigungsfehler (2·δi/tan i) im Fehler: bestes a₀ algebraisch {aA2:.3f} (χ²/N {cA2/N:.2f}), Feldform dünn {aF2:.3f} (χ²/N {cF2/N:.2f})")
W = Wsave
mk = ok & (Qp == 1); r = glob3(mk)
l = np.log10(gb/(1.2*A10)); dD = (lnob-m_alg(gb, 1.2))/np.log(10); dF = np.log10(PT["Rd"][:, jM])
P(f"   nur Q = 1 ({len(np.unique(gi[mk]))} Galaxien): algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f};  {bet_str(mk, dD, dF, l)}")
first = np.zeros(len(gb), bool)
for g in ugal: first[IDX[g][:3]] = True
mk = ok & ~first; r = glob3(mk)
P(f"   ohne die drei innersten Punkte je Galaxie: algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f};  {bet_str(mk, dD, dF, l)}")
mk = ok & (sig < 0.1); r = glob3(mk)
P(f"   nur Punkte mit 2·errV/Vobs < 0,1: algebraisch {r[0]:.3f}, Feldform dünn {r[1]:.3f}, dick {r[2]:.3f};  {bet_str(mk, dD, dF, l)}")
xs_ = ZPG*YD*PT["Ss"]/gb; ref = ok & (xs_ < 0.3)
P(f"   Punkte mit x ≥ 1 ohne Gas in Σ: {np.sum(ok & (xs_ >= 1))} statt {np.sum(ok & (x >= 1))}")
P()
P(f"(Rechenzeit {time.time()-t00:.0f} s)")
open(os.path.join(HIER, "sparc_test.txt"), "w", encoding="utf-8").write("\n".join(out)+"\n")
