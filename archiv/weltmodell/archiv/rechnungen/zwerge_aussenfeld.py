# -*- coding: utf-8 -*-
# Zwerge DF2 und DF4 im homogenen äußeren Feld des Nachbarn.
# Feldform ∇²Φ = ∇·[ν(|∇Φ_N,innen + g_e|)(∇Φ_N,innen + g_e)] − Außenfeld, achsensymmetrisch in (R, z) über die volle z-Achse,
# Finite Volumen (Methode wie scheibe_feldform.py), dazu algebraisch gleichgerichtet (s² = s_i² + s_e²) und senkrecht (s⁴ = s_i⁴ + s_e⁴).
# Streuung aus der Jeans-Gleichung (kugelförmig, isotrop, Monopol des Feldes), Blickrichtung längs/quer zum Feld aus dem Tensor-Virialsatz.
# Einheiten: Länge kpc, Beschleunigung a₀ = 1,2·10⁻¹⁰ m/s² = 1, G·M in kpc²·a₀; σ² [kpc·a₀] → (km/s)² mit a₀ = 3702,8 (km/s)²/kpc.
# [A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.
import numpy as np, sys, io, os
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []; P = lambda s="": out.append(s)

A0 = float(__import__("os").environ.get("FTM_A0", "1.2e-10")) * 3.0856775814913673e19 / 1e6       # a₀ in (km/s)²/kpc
GK = 4.30091e-6                                   # G in kpc (km/s)²/M☉
GMu = lambda M: GK * M / A0                       # G·M in kpc²·a₀
kms = lambda s2: np.sqrt(np.maximum(s2, 0) * A0)  # σ² [kpc·a₀] → km/s
def nu(y):
    y = np.maximum(np.asarray(y, float), 1e-300); return 1.0 / (-np.expm1(-np.sqrt(y)))
def Lnu(y):                                        # d ln ν / d ln g_N
    return -(nu(y) - 1) * np.sqrt(y) / 2

# ---------------------------------------------------------------- Eingaben
# M = Sternmasse mit Υ = 2 [G], R_e [G], Nachbar NGC 1052 baryonisch [A/G], projizierter Abstand [G], σ gemessen [G]
DW = {"DF2": dict(M=2.0e8, Re=2.2, Mh=1.0e11, D=80.0,  so=8.5, ep=2.3, em=3.1),
      "DF4": dict(M=1.5e8, Re=1.6, Mh=1.0e11, D=165.0, so=4.2, ep=4.4, em=2.2)}
gNe = lambda d, D: GMu(d["Mh"]) / D**2             # Newton-Feld des Nachbarn (Punktmasse) in a₀
P("0) Eingaben [G/A]: M (Υ = 2), R_e, Nachbar M_h, projizierter Abstand D, σ gemessen")
for n, d in DW.items():
    P(f"   {n}: M = {d['M']:.1e} M☉, R_e = {d['Re']} kpc, M_h = {d['Mh']:.0e} M☉, D = {d['D']} kpc  ->  g_Ne = {gNe(d, d['D']):.4f} a₀, "
      f"ν(g_Ne) = {nu(gNe(d, d['D'])):.2f};  r* = {np.sqrt(GMu(d['M'])):.3f} kpc;  σ_obs = {d['so']} +{d['ep']}/−{d['em']} km/s")
P()

MU, WMU = np.polynomial.legendre.leggauss(400)
def phantom(g, ge, kind):
    """Monopol der Strukturbeschleunigung (einwärts) bei innerem Newton-Feld g (Betrag, radial) und äußerem g_e.
       feld: Flussmittel über die Kugel, exakt für die Feldform (Gauß); gleich: s² addieren; senk: s⁴ addieren."""
    g = np.atleast_1d(np.asarray(g, float))
    if kind == "newton": return 0 * g
    if kind == "iso":    return (nu(g) - 1) * g
    if kind == "gleich": return (nu(g + ge) - 1) * g
    if kind == "senk":   return (nu(np.hypot(g, ge)) - 1) * g
    gs = g[:, None]; m = MU[None, :]
    gm = np.sqrt(np.maximum(gs**2 + ge**2 + 2 * gs * ge * m, 0))
    return 0.5 * np.sum(WMU * (nu(gm) - 1) * (gs + ge * m), axis=1)

# ---------------------------------------------------------------- 1) Rekonstruktion der 18–20 km/s
P("1) Rekonstruktion der bisherigen Zahl: Wolf-Schätzer σ² = r_½·g(r_½)/3, r_½ = 4/3·R_e, M(<r_½) = M/2 als Punktmasse,")
P("   äußeres Feld = Newton-Feld des Nachbarn (Punktmasse 10¹¹ M☉ im projizierten Abstand), keine Glättung des Zwergs (Sterne ℓ ≪ Zwerg)")
P("        Newton  isoliert  gleichger.  senkrecht  Flussmittel(=Feldform)   [km/s]")
for n, d in DW.items():
    rh = 4 / 3 * d["Re"]; g = GMu(d["M"]) / 2 / rh**2; ge = gNe(d, d["D"])
    row = [kms(rh * (g + phantom(g, ge, k)[0]) / 3) for k in ("newton", "iso", "gleich", "senk", "feld")]
    P(f"   {n}  " + "  ".join(f"{v:7.2f}" for v in row))
P("   Das Flussmittel ist der Monopol der Feldform (Satz von Gauß): die 18,0/19,7 km/s sind schon die Feldform, im Wolf-Schätzer.")
P()

# ---------------------------------------------------------------- Plummer und Jeans
rj = np.logspace(-3, np.log10(3000), 2000)          # radiales Gitter [kpc]
def plummer(d):
    b = d["Re"]; GM = GMu(d["M"])
    GMenc = GM * rj**3 / (rj**2 + b**2)**1.5
    rho = (1 + rj**2 / b**2)**-2.5                  # Form, Normierung kürzt sich
    return b, GM, GMenc, rho
def Fshell(r, rp, l):
    """Anteil einer mit dem Plummer-Kern ℓ geglätteten Kugelschale (Radius r') innerhalb des Radius r (geschlossen)."""
    t1, t2 = r - rp, r + rp
    return ((rp * t1 - l * l) / np.sqrt(1 + t1**2 / l**2) + (l * l + rp * t2) / np.sqrt(1 + t2**2 / l**2)) / (2 * l * rp)
_Fcache = {}
def smooth(GMenc, l):
    if l <= 0: return GMenc
    if l not in _Fcache:
        rm = np.sqrt(rj * np.concatenate([[rj[0] / 1.0075], rj[:-1]]))      # Schalenmitte
        _Fcache[l] = Fshell(rj[:, None], rm[None, :], l)
    dM = np.diff(GMenc, prepend=0.0)
    return _Fcache[l] @ dM
def g_monopol(d, ge, kind, l=0.0):
    """Einwärtsbeschleunigung (Monopol) auf rj: Newton + über ℓ geglättete Struktur aus dem über ℓ geglätteten Feld (Hierarchie E0)."""
    b, GM, GMenc, rho = plummer(d)
    gbar = smooth(GMenc, l) / rj**2
    GMs = rj**2 * phantom(gbar, ge, kind)
    return GMenc / rj**2 + smooth(GMs, l) / rj**2
def jeans(d, gr):
    """isotrope Jeans-Gleichung; σ_los im Fenster R < R_e (leuchtkraftgewichtet), Virialwert, Wolf-Wert"""
    b, GM, GMenc, rho = plummer(d)
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(rj))[::-1])[::-1], [0.0]])
    Rap = d["Re"]; w = np.where(rj <= Rap, 1.0, 1 - np.sqrt(np.clip(1 - (Rap / rj)**2, 0, 1)))
    ap = np.trapezoid(rj**2 * p * w, rj) / np.trapezoid(rj**2 * rho * w, rj)
    vir = np.trapezoid(rj**3 * rho * gr, rj) / (3 * np.trapezoid(rj**2 * rho, rj))
    rh = b / np.sqrt(2**(2 / 3) - 1)
    wolf = rh * np.interp(rh, rj, gr) / 3
    return kms(ap), kms(vir), kms(wolf)
P("   Probe Glättung: Plummer-Kern auf Punktmasse, Anteil in r = ℓ: geschlossen "
  f"{Fshell(1.0, 1e-4, 1.0):.5f}, erwartet 2^(-3/2) = {2**-1.5:.5f}")
P()

# ---------------------------------------------------------------- 2) Feldform auf dem Gitter
def sgrid(N, L, h0):
    lo, hi = 1e-3, 60.0
    for _ in range(200):
        a = 0.5 * (lo + hi); w = L * np.sinh(a / N) / np.sinh(a)
        if w > h0: lo = a
        else: hi = a
    return L * np.sinh(a * np.arange(N + 1) / N) / np.sinh(a)
Lb = 500.0
Rf = sgrid(220, Lb, 0.03); zp = sgrid(220, Lb, 0.03); zf = np.concatenate([-zp[::-1], zp[1:]])
Rc = 0.5 * (Rf[1:] + Rf[:-1]); zc = 0.5 * (zf[1:] + zf[:-1]); NR, NZ = len(Rc), len(zc)
idx = (np.arange(NR)[:, None] + NR * np.arange(NZ)[None, :])            # idx[i, j]
RC, ZC = np.meshgrid(Rc, zc, indexing="ij")
dz = np.diff(zf); dR2 = np.diff(Rf**2); VOL = np.pi * dR2[:, None] * dz[None, :]
faces = []   # (k, m, c, Rface, zface) innen
cR = 2 * np.pi * Rf[1:-1, None] * dz[None, :] / (Rc[1:] - Rc[:-1])[:, None]
faces.append((idx[:-1, :].ravel(), idx[1:, :].ravel(), cR.ravel(),
              np.broadcast_to(Rf[1:-1, None], cR.shape).ravel(), np.broadcast_to(zc[None, :], cR.shape).ravel()))
cz = np.pi * dR2[:, None] / (zc[1:] - zc[:-1])[None, :]
faces.append((idx[:, :-1].ravel(), idx[:, 1:].ravel(), cz.ravel(),
              np.broadcast_to(Rc[:, None], cz.shape).ravel(), np.broadcast_to(zf[1:-1][None, :], cz.shape).ravel()))
bnd = [(idx[-1, :], 2 * np.pi * Rf[-1] * dz / (Rf[-1] - Rc[-1]), np.full(NZ, Rf[-1]), zc),
       (idx[:, 0], np.pi * dR2 / (zc[0] - zf[0]), Rc, np.full(NR, zf[0])),
       (idx[:, -1], np.pi * dR2 / (zf[-1] - zc[-1]), Rc, np.full(NR, zf[-1]))]
rows, cols, vals = [], [], []
for k, m, c, _, _ in faces:
    rows += [k, k, m, m]; cols += [k, m, m, k]; vals += [-c, c, -c, c]
for k, c, _, _ in bnd:
    rows.append(k); cols.append(k); vals.append(-c)
LAP = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(NR * NZ, NR * NZ))
LU = splu(LAP)
kc = idx.ravel(order="C"); Rcell = RC.ravel(); Zcell = ZC.ravel()
order = np.argsort(kc); Rk, Zk = Rcell[order], Zcell[order]          # Zellmitten in Index-Reihenfolge

def solve2d(d, ge, phib):
    """φ = Φ − Φ_außen,  ∇²φ = ∇·[ν(|∇Φ_N|)∇Φ_N − ν_e∇Φ_N,e],  Φ_N,e = g_e·z (Feld zeigt nach −z)."""
    b, GM = d["Re"], GMu(d["M"])
    PhiI = lambda R, z: -GM / np.sqrt(R * R + z * z + b * b)
    def gabs(R, z):
        q = GM / (R * R + z * z + b * b)**1.5
        return np.hypot(q * R, q * z + ge)
    ne = nu(ge) if ge > 0 else 0.0
    Pc = PhiI(Rk, Zk); Pe = ge * Zk
    rhs = np.zeros(NR * NZ)
    for k, m, c, Rf_, zf_ in faces:
        nf = nu(gabs(Rf_, zf_))
        D = c * (nf * (Pc[m] - Pc[k]) + (nf - ne) * (Pe[m] - Pe[k]))
        np.add.at(rhs, k, D); np.add.at(rhs, m, -D)
    for k, c, Rf_, zf_ in bnd:
        nf = nu(gabs(Rf_, zf_))
        D = c * (nf * (PhiI(Rf_, zf_) - Pc[k]) + (nf - ne) * (ge * zf_ - Pe[k]))
        np.add.at(rhs, k, D - c * phib(Rf_, zf_))
    phi = LU.solve(rhs)
    F = phi[idx]                                                      # (NR, NZ)
    dR, dZ = np.gradient(F, Rc, zc)
    return F, dR, dZ
def monopol2d(dR, dZ, r):
    iR = RegularGridInterpolator((Rc, zc), dR, bounds_error=False, fill_value=None)
    iZ = RegularGridInterpolator((Rc, zc), dZ, bounds_error=False, fill_value=None)
    mu = MU[None, :]; rr = r[:, None]; st = np.sqrt(1 - mu**2)
    pts = np.stack([np.broadcast_to(rr * st, (len(r), len(MU))), np.broadcast_to(rr * mu, (len(r), len(MU)))], -1)
    return 0.5 * np.sum(WMU * (st * iR(pts) + mu * iZ(pts)), axis=1)          # ∂φ/∂r gemittelt = Einwärtsbeschleunigung
def far_efe(d, ge):
    GM = GMu(d["M"]); ne = nu(ge); L = Lnu(ge)
    return lambda R, z: -ne * GM / np.sqrt(R * R + z * z) * (1 + 0.5 * L * R * R / (R * R + z * z))
def g2d_on_rj(d, ge, dR, dZ):
    rin = rj[rj < 150]; g = monopol2d(dR, dZ, rin)
    GM = GMu(d["M"]); gout = nu(ge) * (1 + Lnu(ge) / 3) * GM / rj[rj >= 150]**2
    return np.concatenate([g, gout])
def tensor_virial(d, dR, dZ):
    b = d["Re"]; rho = (1 + (RC**2 + ZC**2) / b**2)**-2.5
    Mt = np.sum(rho * VOL)
    szz = np.sum(rho * ZC * dZ * VOL) / Mt; sxx = 0.5 * np.sum(rho * RC * dR * VOL) / Mt
    return szz, sxx

P(f"2) Löser: Gitter {NR} × {NZ} Zellen, R ≤ {Lb:.0f} kpc, |z| ≤ {Lb:.0f} kpc, erste Zelle 0,03 kpc; Rand: Fernfeld")
P("   φ_Rand = −ν_e·GM/r·[1 + (L/2)·sin²θ],  L = d ln ν/d ln g_N bei g_Ne,  θ gegen die Feldrichtung (linearisiert, exakt für r → ∞)")
# Prüfung a) isolierte Plummer-Kugel: Feldform = algebraische Regel exakt
d = DW["DF2"]; b, GM, GMenc, rho = plummer(d)
galg = GMenc / rj**2 * nu(GMenc / rj**2)
Phialg = np.concatenate([[0.0], np.cumsum(0.5 * (galg[1:] + galg[:-1]) * np.diff(rj))])
F, dR, dZ = solve2d(d, 0.0, lambda R, z: np.interp(np.sqrt(R * R + z * z), rj, Phialg))
rt = np.array([0.5, 1.0, 2.2, 5.0, 10.0, 30.0]); gm = monopol2d(dR, dZ, rt); gex = np.interp(rt, rj, galg)
P("   a) DF2 isoliert, Feldform auf dem Gitter / exakt ν(g_N)·g_N:  " + ",  ".join(f"r = {x:g}: {a/e:.4f}" for x, a, e in zip(rt, gm, gex)))
# Prüfung b) mit äußerem Feld: Monopol auf dem Gitter / Flussidentität
ge = gNe(d, d["D"])
F, dR, dZ = solve2d(d, ge, far_efe(d, ge)); gm = monopol2d(dR, dZ, rt)
gfl = np.interp(rt, rj, g_monopol(d, ge, "feld"))
P(f"   b) DF2 bei g_Ne = {ge:.4f} a₀, Monopol Gitter / Flussidentität:  " + ",  ".join(f"r = {x:g}: {a/e:.4f}" for x, a, e in zip(rt, gm, gfl)))
P()

# ---------------------------------------------------------------- 3) Vergleich im projizierten Abstand
P("3) Streuung im projizierten Abstand (Plummer, b = R_e [A]; isotrope Jeans-Gleichung, σ_los im Fenster R < R_e) [km/s]")
P("        Newton  isoliert  gleichger.  senkrecht  Feldform(Gitter)  Feldform(Fluss)   längs Feld  quer Feld   Virial  Wolf   gemessen [G]")
RES = {}
for n, d in DW.items():
    ge = gNe(d, d["D"])
    s = {k: jeans(d, g_monopol(d, ge, k)) for k in ("newton", "iso", "gleich", "senk", "feld")}
    F, dR, dZ = solve2d(d, ge, far_efe(d, ge)); sg = jeans(d, g2d_on_rj(d, ge, dR, dZ))
    szz, sxx = tensor_virial(d, dR, dZ); tot = (szz + 2 * sxx) / 3
    lang, quer = sg[0] * np.sqrt(szz / tot), sg[0] * np.sqrt(sxx / tot)
    RES[n] = dict(s=s, sg=sg, szz=szz, sxx=sxx)
    P(f"   {n}  {s['newton'][0]:6.2f}  {s['iso'][0]:7.2f}  {s['gleich'][0]:9.2f}  {s['senk'][0]:9.2f}  {sg[0]:13.2f}  {s['feld'][0]:14.2f}"
      f"   {lang:9.2f}  {quer:9.2f}   {sg[1]:6.2f}  {sg[2]:5.2f}   {d['so']} +{d['ep']}/−{d['em']}")
P("   längs/quer: Blickrichtung parallel bzw. senkrecht zum äußeren Feld, Verhältnis aus dem Tensor-Virialsatz (σ_zz²/σ_xx² = "
  + ", ".join(f"{n} {RES[n]['szz']/RES[n]['sxx']:.3f}" for n in DW) + ")")
P()

# ---------------------------------------------------------------- 4) Abstand zum Nachbarn
P("4) Abhängigkeit vom 3D-Abstand (≥ projiziert), Nachbar 10¹¹ M☉ [A]; σ_los(R < R_e) [km/s]")
P("        D/D_proj   D[kpc]   g_Ne/a₀   gleichger.  senkrecht  Feldform(Gitter)")
for n, d in DW.items():
    for f in (1.0, 1.5, 2.0, 3.0):
        D = f * d["D"]; ge = gNe(d, D)
        F, dR, dZ = solve2d(d, ge, far_efe(d, ge)); sg = jeans(d, g2d_on_rj(d, ge, dR, dZ))
        P(f"   {n}   {f:5.1f}   {D:7.0f}   {ge:8.4f}   {jeans(d, g_monopol(d, ge, 'gleich'))[0]:9.2f}  {jeans(d, g_monopol(d, ge, 'senk'))[0]:9.2f}  {sg[0]:12.2f}")
d = DW["DF2"]; ge0 = gNe(d, d["D"]); geV = nu(ge0) * ge0
F, dR, dZ = solve2d(d, geV, far_efe(d, geV)); sV = jeans(d, g2d_on_rj(d, geV, dR, dZ))[0]
P(f"   Vergleich DF2: ein Feld vom Betrag des verstärkten Nachbarfelds ν·g_Ne = {geV:.3f} a₀ im Argument gäbe {sV:.2f} km/s")
P()

# ---------------------------------------------------------------- 5) Glättungsradius der Hierarchie
def ellJ(d, D):
    """ℓ³ = GM/(√(3/2)·|E^TF|), E = Hesse-Form des verstärkten Nachbarfelds g(D) = ν(g_N)g_N"""
    GMh = GMu(d["Mh"]); gD = lambda x: nu(GMh / x**2) * GMh / x**2
    h = 1e-4 * D; gp = (gD(D + h) - gD(D - h)) / (2 * h); g = gD(D)
    E = np.array([gp, g / D, g / D]); ETF = E - E.sum() / 3
    return (GMu(d["M"]) / (np.sqrt(1.5) * np.linalg.norm(ETF)))**(1 / 3)
def ell_star(d):
    """ℓ eines Sterns von 1 M☉ bei r_½ im Feld des Zwergs (Feldform-Monopol, projizierter Abstand)"""
    rh = d["Re"] / np.sqrt(2**(2 / 3) - 1); gr = g_monopol(d, gNe(d, d["D"]), "feld")
    g = np.interp(rh, rj, gr); gp = np.interp(rh, rj, np.gradient(gr, rj))
    E = np.array([gp, g / rh, g / rh]); ETF = E - E.sum() / 3
    return (GMu(1.0) / (np.sqrt(1.5) * np.linalg.norm(ETF)))**(1 / 3)
P("5) Glättung ℓ des Zwergs im kollektiven Feld (Hierarchie E0: Argument, Quelle und gespürte Struktur über ℓ geglättet), projizierter Abstand")
P("   Kriterium des Modells: Der Zwerg ist kollektiv (g_N(r_½)/a₀ = "
  + ", ".join(f"{n} {GMu(d['M'])/2/(d['Re']/np.sqrt(2**(2/3)-1))**2:.3f}" for n, d in DW.items())
  + "); geglättet werden nur seine Sterne (ℓ eines Sterns von 1 M☉ bei r_½: "
  + ", ".join(f"{n} {1e3*ell_star(d):.1f} pc" for n, d in DW.items()) + "), also gilt ℓ = 0.")
for n, d in DW.items():
    ge = gNe(d, d["D"]); rh = d["Re"] / np.sqrt(2**(2 / 3) - 1); lJ = ellJ(d, d["D"])
    ls = [0.0, 0.5 * rh, rh, lJ, 2 * lJ, 5 * lJ]
    txt = ",  ".join(f"ℓ = {l:.2f} kpc: {jeans(d, g_monopol(d, ge, 'feld', l))[0]:.2f}" for l in ls)
    P(f"   {n} (r_½ = {rh:.2f} kpc, ℓ_J = {lJ:.2f} kpc):  {txt};  Newton {RES[n]['s']['newton'][0]:.2f} km/s")
P()

# ---------------------------------------------------------------- 6) Was fehlt
P("6) Was für die gemessenen Werte nötig wäre (Feldform, Jeans, Fenster R < R_e)")
for n, d in DW.items():
    ge0 = gNe(d, d["D"]); sN = RES[n]["s"]["newton"][0]
    for tgt, lab in ((d["so"] + d["ep"], "oberer Rand"), (d["so"], "Mittelwert")):
        if tgt <= sN * 1.001:
            M_need = brentq(lambda M: jeans(dict(d, M=M), g_monopol(dict(d, M=M), ge0, "newton"))[0] - tgt, 1e5, 1e10)
            P(f"   {n} {lab} {tgt:.1f} km/s: liegt unter Newton ({sN:.2f}); schon Newton verlangt M ≤ {M_need:.2e} M☉ (Υ ≤ {2*M_need/d['M']:.2f})")
            continue
        gneed = brentq(lambda g: jeans(d, g_monopol(d, g, "feld"))[0] - tgt, 1e-4, 1e3)
        Dneed = np.sqrt(GMu(d["Mh"]) / gneed); Mhneed = gneed * d["D"]**2 * A0 / GK
        Mneed = brentq(lambda M: jeans(dict(d, M=M), g_monopol(dict(d, M=M), ge0, "feld"))[0] - tgt, 1e4, 1e10)
        lneed = brentq(lambda l: jeans(d, g_monopol(d, ge0, "feld", l))[0] - tgt, 0.05, 500)
        P(f"   {n} {lab} {tgt:.1f} km/s: g_Ne ≥ {gneed:.3f} a₀ (= {gneed/ge0:.1f}·heute; bei 10¹¹ M☉: D ≤ {Dneed:.1f} kpc; bei D = {d['D']:.0f} kpc: "
          f"M_h ≥ {Mhneed:.1e} M☉);  oder M ≤ {Mneed:.2e} M☉ (Υ ≤ {2*Mneed/d['M']:.2f});  oder ℓ ≥ {lneed:.1f} kpc")
cr = dict(M=3.2e5, Re=1.066, Mh=6e10, D=117.5, so=2.7, ep=0.3, em=0.3)       # Crater II an der Milchstraße [G/A]
gcr = gNe(cr, cr["D"]); lcr = ellJ(cr, cr["D"])
P(f"   Gegenprobe Crater II [G] (3,2·10⁵ M☉, R_e = 1,07 kpc, Milchstraße 6·10¹⁰ M☉ bei 117,5 kpc, σ = 2,7 ± 0,3): "
  f"Newton {jeans(cr, g_monopol(cr, gcr, 'newton'))[0]:.2f}, Feldform ℓ = 0: {jeans(cr, g_monopol(cr, gcr, 'feld'))[0]:.2f}, "
  f"ℓ = ℓ_J = {lcr:.2f} kpc: {jeans(cr, g_monopol(cr, gcr, 'feld', lcr))[0]:.2f} km/s")
F, dR, dZ = solve2d(DW["DF2"], 0.30, far_efe(DW["DF2"], 0.30))
P(f"   Probe Gitter bei g_Ne = 0,30 a₀ (DF2): {jeans(DW['DF2'], g2d_on_rj(DW['DF2'], 0.30, dR, dZ))[0]:.2f} km/s, "
  f"Flussidentität {jeans(DW['DF2'], g_monopol(DW['DF2'], 0.30, 'feld'))[0]:.2f} km/s")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zwerge_aussenfeld.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
