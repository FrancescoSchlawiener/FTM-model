# Unabhängige Nachrechnung Notiz 74: Gewicht f des Scheibenfelds an SPARC.
# Eigener Code: Gas-Flächendichte über Toomre-Inversion mit elliptischem Integral,
# Unsicherheit über Jackknife der Galaxien (statt Ziehen), a0 auf Gitter profiliert.
import numpy as np, glob, os
from scipy.special import ellipk
D0 = r"C:\DEV\FTM\weltmodell\daten\sparc"
KPC = 3.0857e19; G = 6.674e-11; MSUN = 1.989e30; PC = 3.0857e16
GK = 1e6/KPC                       # (km/s)^2/kpc -> m/s^2
ZS = 2*np.pi*G*MSUN/PC**2          # 2πGΣ pro (Msun/pc^2) in m/s^2
GKPC = G*MSUN/(PC*1e3)/1e6         # G in (km/s)^2 kpc / Msun

tab = {}
for line in open(os.path.join(D0, "SPARC_Lelli2016c.mrt"), encoding="latin-1"):
    p = line.split()
    if len(p) > 17 and p[1].isdigit() and len(p[0]) > 1:
        try:
            tab[p[0]] = dict(D=float(p[2]), eD=float(p[3]), fD=int(p[4]), inc=float(p[5]), einc=float(p[6]), Q=int(p[17]))
        except ValueError:
            pass

def sigma_gas(R, Vg):
    # Toomre (1963): Σ(R) = 1/(π²G) [ (1/R)∫0^R dV²/dr K(r/R) dr + ∫R^∞ dV²/dr K(R/r) dr/r ]
    y = Vg*np.abs(Vg)
    r = np.linspace(1e-3, R[-1], 4000)
    yi = np.interp(r, np.r_[0, R], np.r_[0, y])
    rr = np.r_[r, np.linspace(R[-1]*1.001, 40*R[-1], 4000)]
    yy = np.r_[yi, y[-1]*R[-1]/rr[len(r):]]             # außen Kepler
    dy = np.gradient(yy, rr)
    S = []
    for Rk in R:
        ins = rr < Rk*(1-1e-4); out = rr > Rk*(1+1e-4)
        k1 = ellipk(np.minimum((rr[ins]/Rk)**2, 1-1e-12))
        k2 = ellipk(np.minimum((Rk/rr[out])**2, 1-1e-12))
        v = np.trapezoid(dy[ins]*k1, rr[ins])/Rk + np.trapezoid(dy[out]*k2/rr[out], rr[out])
        S.append(v/(np.pi**2*GKPC)/1e6)                  # Msun/kpc^2 -> Msun/pc^2
    return np.maximum(np.array(S), 0.0)

gal = []
for fn in sorted(glob.glob(os.path.join(D0, "Rotmod_LTG", "*_rotmod.dat"))):
    name = os.path.basename(fn).replace("_rotmod.dat", "")
    t = tab.get(name)
    if t is None or t["Q"] > 2 or t["inc"] < 30: continue
    d = np.loadtxt(fn)
    R, Vo, eV, Vg, Vd, Vb, SBd = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5], d[:, 6]
    m = (Vo > 0) & (R > 0) & (eV > 0)
    Sg = sigma_gas(R, Vg)
    gal.append(dict(name=name, prim=t["fD"] >= 2, R=R[m], Vo=Vo[m], eV=eV[m], Vg=Vg[m], Vd=Vd[m], Vb=Vb[m],
                    SBd=SBd[m], Sg=Sg[m], eD=t["eD"]/t["D"], ei=np.radians(t["einc"]), inc=np.radians(t["inc"])))
print("Galaxien", len(gal), "Primär", sum(g["prim"] for g in gal), "Punkte", sum(len(g["R"]) for g in gal))

def nu(y): return 1.0/(-np.expm1(-np.sqrt(np.maximum(y, 1e-300))))
def lna(g, Yd, f, a0):
    gR = np.maximum((g["Vg"]*np.abs(g["Vg"]) + Yd*g["Vd"]*np.abs(g["Vd"]) + 0.7*g["Vb"]**2)/g["R"]*GK, 1e-14)
    gz = ZS*(Yd*g["SBd"] + g["Sg"])
    l0 = np.log(nu(gR/a0))
    l1 = np.log(nu(np.sqrt(gR**2 + abs(f)*gz**2)/a0))
    return np.log(gR) + (l1 if f >= 0 else 2*l0 - l1)

AG = np.exp(np.linspace(np.log(0.7e-10), np.log(1.9e-10), 61))
FG = np.round(np.arange(-0.5, 1.51, 0.02), 3)

def chi_table(Yd):
    T = np.zeros((len(gal), len(FG), len(AG)))
    for i, g in enumerate(gal):
        d = np.log(g["Vo"]**2/g["R"]*GK)
        base = [lna(g, Yd, 0.0, 1.2e-10)]
        # Fehler je Punkt (Rechnung P): errV, Abstand, Neigung, Υ (0,1 dex)
        dY = (lna(g, Yd*np.exp(0.23), 0.0, 1.2e-10) - lna(g, Yd*np.exp(-0.23), 0.0, 1.2e-10))/2
        s2 = (2*g["eV"]/g["Vo"])**2 + g["eD"]**2 + (2*g["ei"]/np.tan(g["inc"]))**2 + dY**2
        for j, f in enumerate(FG):
            for k, a0 in enumerate(AG):
                T[i, j, k] = np.sum((d - lna(g, Yd, f, a0))**2/s2)
    return T

def best(T, sel):
    C = T[sel].sum(0)                      # (f, a0)
    Cf = C.min(1); j = int(np.argmin(Cf))
    # Parabel um das Minimum
    if 0 < j < len(FG)-1:
        y0, y1, y2 = Cf[j-1], Cf[j], Cf[j+1]; h = FG[1]-FG[0]
        fb = FG[j] + h*(y0-y2)/(2*(y0-2*y1+y2))
    else: fb = FG[j]
    return fb, C

def a0_at(C, f):
    j = int(np.argmin(abs(FG-f))); return AG[int(np.argmin(C[j]))]

names = np.array([g["name"] for g in gal]); prim = np.array([g["prim"] for g in gal])
npt = sum(len(g["R"]) for g in gal)
for Yd in (0.4, 0.5, 0.6):
    T = chi_table(Yd)
    out = []
    for lab, sel in (("Primär", prim), ("Hubble", ~prim), ("alle", np.ones(len(gal), bool)),
                     ("alle ohne IC2574", names != "IC2574"), ("Primär ohne IC2574", prim & (names != "IC2574"))):
        fb, C = best(T, sel)
        idx = np.where(sel)[0]
        jk = np.array([best(T, np.isin(np.arange(len(gal)), idx[idx != i]))[0] for i in idx])
        n = len(idx); se = np.sqrt((n-1)/n*np.sum((jk-jk.mean())**2))
        s2 = C.min()/(sum(len(gal[i]["R"]) for i in idx) - 2)
        a = [a0_at(C, f)*1e10 for f in (0, 0.25, 1.0)]
        out.append(f"  {lab:20s} f = {fb:+.2f} ± {se:.2f} (Jackknife)  χ²/FG {s2:.2f}   a0(f=0/¼/1) = {a[0]:.2f}/{a[1]:.2f}/{a[2]:.2f}")
    print(f"Υ_d = {Yd}"); print("\n".join(out))
    if Yd == 0.5:
        fg = np.array([best(T, names == nm)[0] for nm in names])
        print("  je Galaxie bestes f: unter 0:", int((fg < 0).sum()), " über 0:", int((fg > 0).sum()), " genau 0:", int((fg == 0).sum()))
        # Einfluss einzelner Galaxien auf f (alle)
        idx = np.arange(len(gal)); fall = best(T, np.ones(len(gal), bool))[0]
        infl = sorted(((best(T, idx != i)[0] - fall, names[i]) for i in idx))
        print("  größter Einfluss (f ohne diese Galaxie − f):", [(n_, round(v, 3)) for v, n_ in infl[:3]], [(n_, round(v, 3)) for v, n_ in infl[-2:]])

# Schicht-Mittel von Hand geprüft: F gleichverteilt → ⟨F⟩ = 1/2 → f = 1/4; ⟨F²⟩ = 1/3 → f = 1/3
F = np.linspace(0, 1, 100001); print("⟨F⟩² =", np.mean(F)**2, " ⟨F²⟩ =", np.mean(F**2))
