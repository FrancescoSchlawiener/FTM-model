# Welchen a₀ liest die algebraische Relation a = ν(a_N/a₀)·a_N ab, wenn die Scheiben der Feldform
# ∇²Φ = ∇·[ν(|∇Φ_N|/a₀)∇Φ_N] mit dem wahren a₀ = cH₀/2π folgen?
# Finite-Volumen-Löser in (R, z ≥ 0) übernommen aus scheibe_feldform.py (dort unverändert), hier mit
#   ν = 1/(1 − e^{−s}),  s = √(|g_N|/a₀)   (dort wird ν mit |g_N|² aufgerufen, also s = |g_N|/a₀),
#   Quelle der Kugel mit richtigem Vorzeichen, Rand mit demselben ν; Aufbau der Matrix vektorisiert.
# Einheiten G = M = a₀(wahr) = 1  ⇒  r* = 1, v_f = 1.  Äußeres Feld: senkrecht s⁴ = s_i⁴ + s_e⁴, gleichgerichtet s² = s_i² + s_e².
import numpy as np, sys, io, os, time
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.optimize import minimize_scalar, minimize, brentq
from scipy.integrate import quad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []
def P(s=""):
    out.append(s); print(s, flush=True)

def nu_s(s): return 1.0/(-np.expm1(-np.maximum(s, 1e-12)))
def s_of(g, ge=0.0, mode=None):
    g = np.abs(g)
    if mode is None or ge == 0: return np.sqrt(g)
    if mode == "senk": return (g*g + ge*ge)**0.25
    return np.sqrt(g + ge)                       # gleichgerichtet
def nu_g(g, ge=0.0, mode=None): return nu_s(s_of(g, ge, mode))

# ---------- Gitter und Löser ----------
def grid(N, Lmax, h0):
    lo, hi = 0.1, 20
    for _ in range(100):
        mid = 0.5*(lo+hi); w = Lmax*np.sinh(mid/N)/np.sinh(mid)
        if w > h0: lo = mid
        else: hi = mid
    return Lmax*np.sinh(mid*np.arange(N+1)/N)/np.sinh(mid)

class Loeser:
    def __init__(s, N, Lmax, h0):
        s.Rf = grid(N, Lmax, h0); s.zf = s.Rf.copy()
        s.NR = s.NZ = N
        s.Rc = 0.5*(s.Rf[1:]+s.Rf[:-1]); s.zc = 0.5*(s.zf[1:]+s.zf[:-1])
        dz = np.diff(s.zf); dR2 = np.diff(s.Rf**2)
        s.V = np.pi*dR2[:, None]*dz[None, :]
        s.cR = 2*np.pi*s.Rf[1:-1][:, None]*dz[None, :]/(s.Rc[1:]-s.Rc[:-1])[:, None]      # (NR−1, NZ)
        s.cz = np.pi*dR2[:, None]/(s.zc[1:]-s.zc[:-1])[None, :]                           # (NR, NZ−1)
        s.cRb = 2*np.pi*s.Rf[-1]*dz/(s.Rf[-1]-s.Rc[-1])                                  # äußerer Rand, je j
        s.czb = np.pi*dR2/(s.zf[-1]-s.zc[-1])                                            # oberer Rand, je i
        NR, NZ = N, N
        k = lambda i, j: i + NR*j
        I, J = np.meshgrid(np.arange(NR-1), np.arange(NZ), indexing="ij")
        a1, b1, c1 = k(I, J).ravel(), k(I+1, J).ravel(), s.cR.ravel()
        I, J = np.meshgrid(np.arange(NR), np.arange(NZ-1), indexing="ij")
        a2, b2, c2 = k(I, J).ravel(), k(I, J+1).ravel(), s.cz.ravel()
        A = np.concatenate([a1, a2]); B = np.concatenate([b1, b2]); C = np.concatenate([c1, c2])
        diag = np.zeros(NR*NZ)
        np.add.at(diag, A, -C); np.add.at(diag, B, -C)
        diag[k(NR-1, np.arange(NZ))] -= s.cRb
        diag[k(np.arange(NR), NZ-1)] -= s.czb
        rows = np.concatenate([A, B, np.arange(NR*NZ)]); cols = np.concatenate([B, A, np.arange(NR*NZ)])
        vals = np.concatenate([C, C, diag])
        s.LU = splu(sp.csc_matrix((vals, (rows, cols)), shape=(NR*NZ, NR*NZ)))
        s.rb_out = np.hypot(s.Rf[-1], s.zc); s.rb_top = np.hypot(s.Rc, s.zf[-1])

    def flat(s, F): return F.ravel(order="F")
    def feld(s, x): return x.reshape((s.NR, s.NZ), order="F")

    def poisson(s, Mcell, Mring0, fN):
        """∇²Φ_N = 4πρ; Mcell Masse je Zelle (Halbraum), Mring0 Masse der dünnen Scheibe je Ring (Fluss 2πM am Boden)."""
        rhs = 4*np.pi*Mcell.copy(); rhs[:, 0] += 2*np.pi*Mring0
        rhs[-1, :] -= s.cRb*fN(s.rb_out); rhs[:, -1] -= s.czb*fN(s.rb_top)
        return s.feld(s.LU.solve(s.flat(rhs)))

    def grads(s, F, dphidz0):
        Rc, zc = s.Rc, s.zc
        gR = np.zeros_like(F); gz = np.zeros_like(F)
        gR[1:-1, :] = -(F[2:, :]-F[:-2, :])/(Rc[2:]-Rc[:-2])[:, None]; gR[0, :] = -(F[1, :]-F[0, :])/(Rc[1]-Rc[0]); gR[-1, :] = gR[-2, :]
        gz[:, 1:-1] = -(F[:, 2:]-F[:, :-2])/(zc[2:]-zc[:-2])[None, :]; gz[:, -1] = gz[:, -2]
        gz[:, 0] = -0.5*(dphidz0 + (F[:, 1]-F[:, 0])/(zc[1]-zc[0]))
        return gR, gz

    def feldform(s, FN, dphidz0, ftot, fN, ge=0.0, mode=None):
        """∇²Φ = ∇·[ν ∇Φ_N]; Bodenflüsse heben sich heraus."""
        gR, gz = s.grads(FN, dphidz0)
        gn = -(FN[1:, :]-FN[:-1, :])/(s.Rc[1:]-s.Rc[:-1])[:, None]; gt = 0.5*(gz[1:, :]+gz[:-1, :])
        nR = nu_g(np.hypot(gn, gt), ge, mode)
        gn = -(FN[:, 1:]-FN[:, :-1])/(s.zc[1:]-s.zc[:-1])[None, :]; gt = 0.5*(gR[:, 1:]+gR[:, :-1])
        nz = nu_g(np.hypot(gn, gt), ge, mode)
        rhs = np.zeros_like(FN)
        fl = nR*s.cR*(FN[1:, :]-FN[:-1, :]); rhs[:-1, :] += fl; rhs[1:, :] -= fl
        fl = nz*s.cz*(FN[:, 1:]-FN[:, :-1]); rhs[:, :-1] += fl; rhs[:, 1:] -= fl
        nb = lambda r: nu_g(1/r**2, ge, mode)
        rhs[-1, :] += nb(s.rb_out)*s.cRb*(fN(s.rb_out)-FN[-1, :]) - s.cRb*ftot(s.rb_out)
        rhs[:, -1] += nb(s.rb_top)*s.czb*(fN(s.rb_top)-FN[:, -1]) - s.czb*ftot(s.rb_top)
        F = s.feld(s.LU.solve(s.flat(rhs)))
        return F, gR

fN = lambda r: -1/r
def ftot_fuer(ge=0.0, mode=None):
    if ge == 0: return lambda r: -1/r - np.log(-np.expm1(-1/r))
    a = lambda r: nu_g(1/r**2, ge, mode)/r**2           # Punktmasse am Rand, Potential = −∫_r^∞ a dr
    rr = np.linspace(70, 130, 241)
    tail = nu_g(1/1e6**2, ge, mode)/1e6
    Phi = np.array([-(quad(a, r, 1e6, limit=400)[0] + tail) for r in rr])
    return lambda r: np.interp(r, rr, Phi)

# ---------- Massenmodelle (Gesamtmasse 1) ----------
def scheibe_ring(Rf, Rd, M):
    return M*((1+Rf[:-1]/Rd)*np.exp(-Rf[:-1]/Rd)-(1+Rf[1:]/Rd)*np.exp(-Rf[1:]/Rd))
def sterne(L, Rd, Ms, BT, z0):
    """Exponentialscheibe (dünn oder sech²-dick mit z0) der Masse Ms, dazu Plummer-Bulge b = 0,15 R_d der Masse BT."""
    ring = scheibe_ring(L.Rf, Rd, Ms)
    Mcell = np.zeros((L.NR, L.NZ)); Mring0 = np.zeros(L.NR); d0 = np.zeros(L.NR)
    if z0 > 0: Mcell += ring[:, None]*0.5*np.diff(np.tanh(L.zf/z0))[None, :]
    else: Mring0 = ring; d0 = 2*np.pi*Ms*np.exp(-L.Rc/Rd)/(2*np.pi*Rd*Rd)
    if BT > 0:
        b = 0.15*Rd
        rho = (1+(L.Rc[:, None]**2+L.zc[None, :]**2)/b**2)**-2.5*L.V
        Mcell += rho*(0.5*BT/rho.sum())
    return Mcell, Mring0, d0

# ---------- Anpassung: g_bar = k·g_Sterne + g_Gas ----------
SA = 0.1                 # Fehler je Punkt in ln a (5 % in v)
SK = 0.1*np.log(10)      # Vorgabe für Υ: 0,1 dex um den wahren Wert
def g_mod(d, a0, k):
    gb = k*d[0]+d[1]; return nu_g(gb/a0)*gb
def chi2(d, a0, k): return np.sum((np.log(d[2])-np.log(g_mod(d, a0, k)))**2)/SA**2
B = (np.log(1e-3), np.log(1e3))
def fit_fest(d):
    r = minimize_scalar(lambda la: chi2(d, np.exp(la), 1.0), bounds=B, method="bounded", options={"xatol": 1e-8})
    return np.exp(r.x)
def fit_k(d, a0, vorgabe):
    f = lambda lk: chi2(d, a0, np.exp(lk)) + (lk/SK)**2*vorgabe
    r = minimize_scalar(f, bounds=B, method="bounded", options={"xatol": 1e-8})
    return np.exp(r.x), r.fun
def fit_frei(d, vorgabe):
    f = lambda p: chi2(d, np.exp(p[1]), np.exp(p[0])) + (p[0]/SK)**2*vorgabe
    best = None
    for lk0 in (-1.0, 0.0, 0.5):
        r = minimize(f, [lk0, 0.0], method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 6000})
        if best is None or r.fun < best.fun: best = r
    return np.exp(best.x[0]), np.exp(best.x[1])
def fit_global(D, modus):
    """modus 0: Υ fest; 1: Υ frei je Galaxie mit Vorgabe 0,1 dex; 2: Υ frei ohne Vorgabe."""
    def F(la):
        a0 = np.exp(la)
        if modus == 0: return sum(chi2(d, a0, 1.0) for d in D)
        return sum(fit_k(d, a0, 1.0 if modus == 1 else 0.0)[1] for d in D)
    r = minimize_scalar(F, bounds=(np.log(0.1), np.log(10)), method="bounded", options={"xatol": 1e-6})
    return np.exp(r.x)
def g3(D): return " / ".join(f"{fit_global(D, m):.3f}" for m in (0, 1, 2))

t0 = time.time()
# =====================================================================
P("Einheiten G = M = a₀ = 1 (r* = 1). ν = 1/(1 − e^{−s}), s = √(|g_N|/a₀).")
P()
# ---------- 1) Prüfung des Lösers ----------
P("1) Prüfung des Lösers (Gitter N × N bis 80 r*, erste Zelle 0,002 r*)")
for N in (300, 420):
    L = Loeser(N, 80.0, 0.002)
    b = 0.3
    rho = (3/(4*np.pi*b**3))*(1+(L.Rc[:, None]**2+L.zc[None, :]**2)/b**2)**-2.5*L.V
    Mc = rho*(0.5/rho.sum())
    FN = L.poisson(Mc, np.zeros(L.NR), fN)
    F, gR = L.feldform(FN, np.zeros(L.NR), ftot_fuer(), fN)
    GR, _ = L.grads(F, np.zeros(L.NR))
    P(f"   N = {N}: Plummer-Kugel b = 0,3 (Feldform = algebraisch exakt), Mittelebene:")
    P("       R/r*    a_N exakt   a_N Gitter   a Feldform   a algebraisch   Verhältnis")
    for Rt in (0.1, 0.3, 1.0, 3.0, 10.0, 25.0):
        i = np.argmin(abs(L.Rc-Rt)); R = L.Rc[i]; aNex = R/(R*R+b*b)**1.5
        alg = nu_g(aNex)*aNex
        P(f"     {R:7.3f}   {aNex:9.5f}   {-gR[i,0]:9.5f}   {-GR[i,0]:10.5f}   {alg:11.5f}   {-GR[i,0]/alg:9.5f}")
    for h in (0.3, 3.0):
        ring = (h/np.sqrt(L.Rf[:-1]**2+h*h)-h/np.sqrt(L.Rf[1:]**2+h*h))
        Sig = h/(2*np.pi*(L.Rc**2+h*h)**1.5)
        FN = L.poisson(np.zeros((L.NR, L.NZ)), ring, fN)
        F, gR = L.feldform(FN, 2*np.pi*Sig, ftot_fuer(), fN)
        GR, _ = L.grads(F, 0)
        P(f"   N = {N}: Kuzmin-Scheibe h = {h} r*: exakt a = ν(|g_N|)·g_N,R, |g_N| = 1/(R²+h²)")
        P("       R/h    a exakt    a Gitter   Gitter/exakt   exakt/algebraisch   √cos θ")
        for x in (0.25, 0.5, 1, 2, 5, 8):
            R = x*h
            if R > 30: continue
            i = np.argmin(abs(L.Rc-R)); R = L.Rc[i]
            gabs = 1/(R*R+h*h); gRex = R/(R*R+h*h)**1.5
            ex = nu_g(gabs)*gRex; alg = nu_g(gRex)*gRex
            P(f"     {R/h:6.2f}   {ex:9.5f}   {-GR[i,0]:9.5f}   {-GR[i,0]/ex:9.4f}   {ex/alg:14.4f}   {np.sqrt(gRex/gabs):9.4f}")
    P()
P(f"   (Rechenzeit bis hier {time.time()-t0:.0f} s)")
P()

# ---------- 2) Scheiben ----------
L = Loeser(420, 80.0, 0.002)
Rds = [0.2, 0.3, 0.45, 0.7, 1.0, 1.5, 2.0, 3.0]
FG = dict(zip(Rds, [0.1, 0.1, 0.15, 0.25, 0.35, 0.5, 0.6, 0.7]))      # Gasanteil [A] für die Variante mit Gas
BTs = [(0.2, 0.2), (0.2, 0.4), (0.3, 0.2), (0.3, 0.4), (0.45, 0.2)]
xs = np.arange(0.5, 8.01, 0.5)                                          # Radien der Rotationskurve in R_d
def interp(R, g): return np.exp(np.interp(np.log(R), np.log(L.Rc), np.log(np.maximum(-g, 1e-300))))
def kurve(Rd, BT=0.0, z0=0.0, fg=0.0, ges=((0.0, None),), x=xs):
    Ms = 1-BT-fg
    Mc, M0, d0 = sterne(L, Rd, Ms, BT, z0)
    FNs = L.poisson(Mc, M0, lambda r: -(Ms+BT)/r)
    FN = FNs.copy(); dg = np.zeros(L.NR)
    if fg > 0:
        Rg = 2.5*Rd
        FNg = L.poisson(np.zeros((L.NR, L.NZ)), scheibe_ring(L.Rf, Rg, fg), lambda r: -fg/r)
        FN += FNg; dg = 2*np.pi*fg*np.exp(-L.Rc/Rg)/(2*np.pi*Rg*Rg)
    R = x*Rd
    gs = interp(R, L.grads(FNs, d0)[0][:, 0])
    gg = interp(R, L.grads(FNg, dg)[0][:, 0]) if fg > 0 else np.zeros_like(R)
    res = {}
    for ge, mode in ges:
        F, gR = L.feldform(FN, d0+dg, ftot_fuer(ge, mode), fN, ge, mode)
        res[(ge, mode)] = (gs, gg, interp(R, L.grads(F, 0)[0][:, 0]))
    return res

P("2) Exponentialscheiben (dünn), wahres a₀ = 1: Verhältnis Feldform / algebraische Relation ν(g_R)·g_R in der Ebene")
P("   R_d/r*  2πGΣ₀/a₀ | R/R_d = 0,25   0,5    1     2     3     4     5     6     8")
x2 = np.array([0.25, 0.5, 1, 2, 3, 4, 5, 6, 8])
for Rd in Rds:
    d = kurve(Rd, x=x2)[(0.0, None)]
    r = d[2]/g_mod(d, 1.0, 1.0)
    P(f"   {Rd:5.2f}  {1/Rd**2:8.2f}  |      " + "  ".join(f"{v:.3f}" for v in r))
P()
P("   Umrechnung [A]: a₀ = cH₀/2π (H₀ = 70) ⇒ a₀/(2πG) = 118 M☉/pc², Σ₀ = 118·(r*/R_d)² M☉/pc².")
P("   Beispiele [A]: M = 6·10¹⁰ M☉, R_d = 2,6 kpc → Σ₀ = 1400 M☉/pc², R_d/r* = 0,29;  M = 10⁹ M☉, R_d = 1,5 kpc → Σ₀ = 71, R_d/r* = 1,3;")
P("                  M = 10⁸ M☉, R_d = 1 kpc → Σ₀ = 16, R_d/r* = 2,7.")
P()

# ---------- 3) Anpassung ohne äußeres Feld ----------
P("3) Anpassung mit der algebraischen Relation, Radien 0,5–8 R_d (16 Punkte, Fehler 0,1 in ln a je Punkt)")
P("   fest: Υ wahr (k = 1), nur a₀;  Vorgabe: k und a₀ je Galaxie, k mit Vorgabe 0,1 dex;  ohne Vorgabe: k frei (nur k·a₀ sicher bestimmt)")
samm = []
for Rd in Rds:
    samm.append((f"R_d = {Rd:4.2f}", kurve(Rd)[(0.0, None)]))
for Rd, BT in BTs:
    samm.append((f"R_d = {Rd:4.2f}, B/T {BT:.1f}", kurve(Rd, BT)[(0.0, None)]))
P("   Galaxie                  a₀ fest | Vorgabe: k      a₀   | ohne Vorgabe: k      a₀      k·a₀")
af = []
for lab, d in samm:
    a1 = fit_fest(d); k2, a2 = fit_frei(d, 1.0); k3, a3 = fit_frei(d, 0.0); af.append(a1)
    P(f"   {lab:22s}  {a1:7.4f} |  {k2:7.4f}  {a2:7.4f}  |   {k3:9.3g} {a3:9.3g}  {k3*a3:7.4f}")
af = np.array(af)
D = [d for _, d in samm]; Dd = D[:len(Rds)]
P(f"   a₀ fest je Galaxie: Median {np.median(af):.3f}, von {af.min():.3f} bis {af.max():.3f}")
P(f"   global, alle Punkte gemeinsam (fest / Vorgabe / ohne Vorgabe): alle 13: {g3(D)};  nur die 8 Scheiben ohne Bulge: {g3(Dd)}")
AG = [fit_global(D, m) for m in (0, 1, 2)]
P("   k je Galaxie beim globalen a₀ mit Vorgabe:  " + " ".join(f"{fit_k(d, AG[1], 1.0)[0]:.3f}" for d in D))
P("   k je Galaxie beim globalen a₀ ohne Vorgabe: " + " ".join(f"{fit_k(d, AG[2], 0.0)[0]:.3f}" for d in D))
P()
def a0_bei_kann(D, kann):
    r = minimize_scalar(lambda la: sum(chi2(d, np.exp(la), kann) for d in D), bounds=(np.log(0.1), np.log(10)), method="bounded", options={"xatol": 1e-7})
    return np.exp(r.x)
def kann_fuer(D, ziel):
    return brentq(lambda k: a0_bei_kann(D, k)-ziel, 0.5, 1.5, xtol=1e-5)
P("   Υ fest, aber Υ_angenommen = k_ann·Υ_wahr: angepasstes a₀ (global, 13 Galaxien)")
P("   " + ";  ".join(f"k_ann = {k:.2f} → {a0_bei_kann(D, k):.3f}" for k in (0.8, 0.85, 0.9, 0.95, 1.0)))
P("   verlangtes Υ_wahr/Υ_angenommen für a₀ = 1,20·10⁻¹⁰ bei H₀ = 67,4 / 70 / 73 (Faktor 1,151 / 1,109 / 1,063): Feldform "
  + " / ".join(f"{1/kann_fuer(D, z):.3f}" for z in (1.151, 1.109, 1.063))
  + ";  algebraisch " + " / ".join(f"{1/kann_fuer([(d[0], d[1], g_mod(d, 1.0, 1.0)) for d in D], z):.3f}" for z in (1.151, 1.109, 1.063)))
P()
P("   Andere Radienbereiche, global (fest / Vorgabe / ohne Vorgabe) | a₀ fest je Scheibe R_d = 0,2 … 3:")
for lo, hi in ((1.0, 8.0), (0.5, 5.0), (2.0, 8.0)):
    sel = (xs >= lo-1e-9) & (xs <= hi+1e-9)
    D2 = [(d[0][sel], d[1][sel], d[2][sel]) for d in D]
    P(f"   {lo:.1f}–{hi:.0f} R_d:  {g3(D2)}   | " + " ".join(f"{fit_fest(d):.3f}" for d in D2[:len(Rds)]))
P()

# ---------- 4) Dicke Scheiben ----------
P("4) Dicke Scheibe (sech², z₀ = 0,1 R_d), Rotationskurve in der Mittelebene")
P("   R_d/r* | Feldform/alg. bei R/R_d = 0,5  1  2  3  5  8 | a₀ fest")
sel6 = [np.argmin(abs(xs-v)) for v in (0.5, 1, 2, 3, 5, 8)]
Dt = []
for Rd in Rds:
    d = kurve(Rd, 0.0, 0.1*Rd)[(0.0, None)]; Dt.append(d)
    r = d[2]/g_mod(d, 1.0, 1.0)
    P(f"   {Rd:5.2f}  | " + "  ".join(f"{r[i]:.3f}" for i in sel6) + f" | {fit_fest(d):.4f}")
P(f"   global (fest / Vorgabe / ohne Vorgabe): dick {g3(Dt)};  dünn {g3(Dd)}")
P()

# ---------- 5) Mit Gasscheibe ----------
P("5) Sterne + dünne Gasscheibe (R_g = 2,5 R_d, Gasanteil [A] 0,1 … 0,7 von hoher zu niedriger Flächendichte), Υ skaliert nur die Sterne")
P("   R_d/r*  Gas | a₀ fest | Vorgabe: k      a₀")
Dg = []
for Rd in Rds:
    d = kurve(Rd, 0.0, 0.0, FG[Rd])[(0.0, None)]; Dg.append(d)
    k2, a2 = fit_frei(d, 1.0)
    P(f"   {Rd:5.2f}  {FG[Rd]:.2f} | {fit_fest(d):.4f}  |  {k2:.4f}  {a2:.4f}")
P(f"   global (fest / Vorgabe / ohne Vorgabe): {g3(Dg)}")
P("   verlangtes Υ_wahr/Υ_angenommen (nur Sterne) für a₀ = 1,20·10⁻¹⁰ bei H₀ = 67,4 / 70 / 73: Feldform "
  + " / ".join(f"{1/kann_fuer(Dg, z):.3f}" for z in (1.151, 1.109, 1.063))
  + ";  algebraisch " + " / ".join(f"{1/kann_fuer([(d[0], d[1], g_mod(d, 1.0, 1.0)) for d in Dg], z):.3f}" for z in (1.151, 1.109, 1.063)))
P()

# ---------- 6) Äußeres Feld ----------
P("6) Äußeres Feld g_e, angepasst ohne äußeres Feld; Feldform dünn (13 Galaxien), daneben algebraisch mit g_e")
P("   g_e/a₀  Lage    | Feldform global: fest / Vorgabe / ohne | algebraisch+g_e: fest / Vorgabe | a₀ fest je Scheibe R_d = 0,2 … 3")
for ge in (0.01, 0.03, 0.05):
    for mode in ("senk", "gleich"):
        DF = []; DA = []
        for Rd in Rds:
            DF.append(kurve(Rd, 0.0, 0.0, 0.0, ((ge, mode),))[(ge, mode)])
        for Rd, BT in BTs:
            DF.append(kurve(Rd, BT, 0.0, 0.0, ((ge, mode),))[(ge, mode)])
        for d in DF:
            gb = d[0]+d[1]; DA.append((d[0], d[1], nu_g(gb, ge, mode)*gb))
        P(f"   {ge:5.2f}   {mode:6s}  |   {g3(DF)}   |   {fit_global(DA, 0):.3f} / {fit_global(DA, 1):.3f}   | "
          + " ".join(f"{fit_fest(d):.3f}" for d in DF[:len(Rds)]))
P()

# ---------- 7) Übertrag in m/s² ----------
c = 2.99792458e8; Mpc = 3.0857e22
P("7) Angepasster Wert in m/s² (Faktor × cH₀/2π) gegen 1,20·10⁻¹⁰, dünne Scheiben ohne äußeres Feld (13 Galaxien)")
for H in (67.4, 70.0, 73.0):
    a = c*H*1e3/Mpc/(2*np.pi)
    P(f"   H₀ = {H}: cH₀/2π = {a*1e10:.3f};  verlangt Faktor {1.2e-10/a:.3f};  Υ fest → {AG[0]*a*1e10:.3f}, Vorgabe → {AG[1]*a*1e10:.3f}, ohne Vorgabe → {AG[2]*a*1e10:.3f}  (·10⁻¹⁰ m/s²)")
P()

# ---------- 8) Abhängigkeit von der Flächendichte ----------
P("8) Prüfbare Folge. a) a₀ fest je Scheibe gegen Σ₀ (H₀ = 70):")
P("   " + ";  ".join(f"{118/Rd**2:.0f} → {fit_fest(d):.3f}" for Rd, d in zip(Rds, Dd)) + "   (Σ₀ in M☉/pc² → a₀)")
P("   b) Abweichung von der algebraischen Relation mit wahrem a₀ und wahrem Υ gegen 2πGΣ/g_R (8 dünne Scheiben, alle Punkte):")
X = []; Y = []
for Rd, d in zip(Rds, Dd):
    Sig = np.exp(-xs)/(2*np.pi*Rd*Rd)
    X.append(2*np.pi*Sig/d[0]); Y.append(np.log10(d[2]/g_mod(d, 1.0, 1.0)))
X = np.concatenate(X); Y = np.concatenate(Y)
for lo, hi in ((0, 0.3), (0.3, 1), (1, 2), (2, 4)):
    m = (X >= lo) & (X < hi)
    if m.sum(): P(f"   2πGΣ/g_R in [{lo}, {hi}):  {m.sum():3d} Punkte,  Δlog g_obs = {Y[m].mean():+.3f} (von {Y[m].min():+.3f} bis {Y[m].max():+.3f})")
P()
P(f"(Rechenzeit {time.time()-t0:.0f} s)")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "schwelle_anpassung.txt"), "w", encoding="utf-8").write("\n".join(out)+"\n")
