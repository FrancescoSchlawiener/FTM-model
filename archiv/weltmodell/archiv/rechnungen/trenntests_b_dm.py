# -*- coding: utf-8 -*-
"""
Trenntests Bindung durch Struktur (B: L_eff = 𝒫L, keine Quelle) gegen Dunkle Materie (DM: zusätzliche Quelle).
Rechnung, noch nicht im Modell. Formeln nur aus ftm_kern; eigene Daten nur SPARC.

1  η = M_dyn/M_lens: Kugel (exakt 1), Kuzmin-Scheibe (geschlossen), Exponentialscheibe (Feldform = wirbelfreier Teil von ν·g_N,
   hier über die effektive Quelldichte −∇·(ν g_N)/4πG und Hankel-Transformation, unabhängig vom Finite-Volumen-Löser), Haufen.
   η(R) = a_Bahn/a_Licht in der Mittelebene: Bahn folgt der örtlichen Regel ν(|g_N|)·g_N, Licht der Feldform (07).
2  Wachstum f(k, z)/f_ΛCDM: linear (k < k_ta) mit Gewicht w und Übergabe im Hintergrund; jenseits der Umkehrskala
   Kraftverhältnis gebundener Bereiche μ_b = ν(a_N,b(R_ta/2))/(1 + Q) und Massenanteil jenseits der Umkehr.
3  Gedächtnis: Kern bei Bahnfrequenz, Übergabezeit 1/λ_m, Haufenstoß-Anteil nach der Faltung I_F.
4  Äquivalenzprinzip an SPARC: Rest der Regel gegen den örtlichen Gasanteil; Kopplung κ des Gases in g_N.
5  Umgebung/Geschichte an SPARC: Ursa-Major-Mitglieder (f_D = 4) gegen übrige; Form (außen − innen) gegen Typ, Gasanteil, Σ.
6  Handwerte.
"""
import math, os, sys, io
import numpy as np
from scipy.special import j0, j1, erfc
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (G, c, A0, M_sun, pc, kpc, Mpc, Gyr, nu, chi, s_, Q_KOSM, q_KOSM, OMEGA_M, OMEGA_B, H0_si, h_,
                      E_hintergrund, gewicht_w, tau_D, uebergabezeit, dichte_umkehr_heute, rho_m_mitbewegt, r_stern,
                      nu_aussenfeld)
import sparc_flaechendichte as sf
import sparc_universal as su

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
OUT = []


def P(s=""):
    OUT.append(s); print(s)


def H(t):
    P(); P("=" * 110); P(t); P("=" * 110)


# ======================================================================================= 1  η = M_dyn/M_lens
H("1  Dynamische gegen Linsenmasse, η(R) = a_Bahn/a_Licht in der Mittelebene (Einheiten G = M = a₀ = 1, r* = 1)")
P("   Kugel: ν(|g_N|)·g_N ist radial und hängt nur von r ab ⇒ wirbelfrei, Feldform = örtliche Regel, η = 1 exakt (jedes Profil, jede Tiefe).")

nu1 = lambda y: nu(np.maximum(y, 1e-300), 1.0)
NR, NZ, NK, RMAX = 1600, 400, 6000, 3e3
Rg = np.logspace(-3.5, math.log10(RMAX), NR)
zg = np.logspace(-4.5, math.log10(RMAX), NZ)
kg = np.logspace(-5, 4, NK)
dk = np.gradient(np.log(kg)) * kg
dR = np.gradient(Rg); dz = np.gradient(zg)
J0m = j0(np.outer(Rg, kg)); J1m = j1(np.outer(Rg, kg))
J0f = J0m * np.exp(-np.outer(dR, kg)**2 / 8)      # Vorwärts-Hankel: Quelle je Zelle gaußisch geglättet (kein Aliasing)


def feldform(St, Sig):
    """St(k) = ∫Σ J0 R dR, Sig(R) = Σ(R). Feldform = Newton-Feld der effektiven Quelle ρ_eff = −∇·(ν g_N)/4π:
    Fläche ν₀Σ (ν₀ an |g_N(R, 0±)|) plus Volumen −g_N·∇ν/4π. Gibt |g_N,R|, a_Bahn = ν(|g_N|)|g_N|, a_Licht (Ebene)."""
    Sk = St(kg)
    gR0 = -2 * math.pi * J1m @ (kg * Sk * dk)
    gz0 = -2 * math.pi * Sig(Rg)
    gR = np.empty((NR, NZ)); gz = np.empty((NR, NZ))
    for j, z in enumerate(zg):
        e = np.exp(-kg * z) * kg * Sk * dk
        gR[:, j] = -2 * math.pi * J1m @ e
        gz[:, j] = -2 * math.pi * J0m @ e
    nuv = nu1(np.hypot(gR, gz))
    rho_v = -(gR * np.gradient(nuv, Rg, axis=0) + gz * np.gradient(nuv, zg, axis=1)) / (4 * math.pi)
    nu0 = nu1(np.hypot(gR0, gz0))
    S_eff = J0f.T @ (nu0 * Sig(Rg) * Rg * dR)
    rt = J0f.T @ (rho_v * (Rg * dR)[:, None])
    S_v = 2 * (rt * np.exp(-np.outer(kg, zg))) @ dz
    aF = -2 * math.pi * J1m @ (kg * (S_eff + S_v) * dk)
    return np.abs(gR0), nu1(np.abs(gR0)) * np.abs(gR0), np.abs(aF)


P()
P("   Prüfung an der Kuzmin-Scheibe (geschlossen: Feldform in der Ebene = ν(1/(R²+h²))·R/(R²+h²)^{3/2}, Bahn = ν(g_R)·g_R):")
P("      h/r*    R/h    η exakt   η Hankel   a_Feld Hankel/exakt")
for h in (0.3, 3.0):
    gN, aal, aF = feldform(lambda k: np.exp(-k * h) / (2 * math.pi), lambda R: h / (2 * math.pi * (R * R + h * h)**1.5))
    for x in (0.25, 1.0, 3.0, 10.0):
        R = x * h
        gRx = R / (R * R + h * h)**1.5
        aFx = nu1(1 / (R * R + h * h)) * gRx
        ex = nu1(gRx) * gRx / aFx
        aFh = np.interp(R, Rg, aF); aalh = np.interp(R, Rg, aal)
        P(f"      {h:4.1f}   {x:5.2f}    {ex:.4f}    {aalh/aFh:.4f}      {aFh/aFx:.4f}")

P()
P("   Dünne Exponentialscheibe, η(R) = a_Bahn/a_Licht (Mittelebene):")
Rd_list = (0.1, 0.34, 1.0, 3.0)
xs = (0.25, 0.5, 1.0, 2.0, 3.0, 3.5, 5.0, 8.0, 15.0)
P("      R_d/r*  " + "".join(f"  R={x:>4}R_d" for x in xs))
ETA = {}
for Rd in Rd_list:
    gN, aal, aF = feldform(lambda k, Rd=Rd: 1 / (2 * math.pi) / (1 + (k * Rd)**2)**1.5,
                           lambda R, Rd=Rd: np.exp(-R / Rd) / (2 * math.pi * Rd * Rd))
    et = [np.interp(x * Rd, Rg, aal) / np.interp(x * Rd, Rg, aF) for x in xs]
    ETA[Rd] = et
    P(f"      {Rd:5.2f}   " + "".join(f"   {e:8.3f}" for e in et))
P("      (R_d/r* = 0,34: Milchstraße 6·10¹⁰ M☉, R_d = 3 kpc, r* = %.1f kpc; Feldform-Löser (Finite Volumen) bei R_d = 1 r*: η = 1/0,605, 1/0,829, 1/1,021 bei 0,25/1/5 R_d)"
  % (r_stern(6e10 * M_sun) / kpc))
P("      Genauigkeit (Kuzmin, feineres Gitter 2400×600×8000: R_d = 0,34 / 1 r*: η = 1,268/1,644 bei R_d/4, 1,097/1,201 bei R_d, 0,988/0,989 bei 5R_d): ±1–2 %")
P("   Kreuzungspunkt η = 1 und Unterschwinger außen (Unterschwinger ≲ 2 % liegt in der Rechengenauigkeit):")
for Rd in Rd_list:
    et = np.array(ETA[Rd]); xx = np.array(xs)
    i = np.where(np.diff(np.sign(et - 1)) != 0)[0]
    xc = float(np.interp(1, et[i[0]:i[0] + 2][::-1], xx[i[0]:i[0] + 2][::-1])) if len(i) else float('nan')
    P(f"      R_d = {Rd:4.2f} r*: η = 1 bei R ≈ {xc:.2f} R_d, kleinstes η außen {et[xx > xc].min() if np.any(xx > xc) else float('nan'):.3f}")

H("1c Haufen: Kugel im Gleichgewicht η = 1 (hydrostatische Masse = Linsenmasse). Stoß: Nachlauf der Mode ξ")
for D in (200, 500):
    w_, _, _ = uebergabezeit(D, z=0.3)
    for IF in (1.0, 1.25):
        P(f"   Δ = {D}, z = 0,3: 1/ω_orb = {w_/Gyr:.2f} Mrd. J.; Anteil der Linse beim Gas null bis t_F = I_F/ω = {IF*w_/Gyr:.2f} Mrd. J. (I_F = {IF}),"
          f" danach 1 − I_F/(ωt): bei t = 2 t_F {1-0.5:.2f}, 4 t_F {1-0.25:.2f}")
P("   Bei t ≲ 0,2 Mrd. J. nach dem Stoß (Bullet-artig) ist der Anteil beim Gas 0: Linse folgt den Galaxien wie bei DM; η(Gas-Ort) ≠ 1 nur später.")

# ======================================================================================= 2  Wachstum
H("2  Wachstum: f(k, z)/f_ΛCDM")
Om, Ob, Q = OMEGA_M, OMEGA_B, Q_KOSM
H0 = H0_si()
Ea = lambda a: E_hintergrund(1 / a - 1)


def dlnE(a, d=1e-5):
    return (math.log(Ea(a * (1 + d))) - math.log(Ea(a * (1 - d)))) / (math.log(1 + d) - math.log(1 - d))


def Om_a(a):
    return Om * a**-3 / Ea(a)**2


def wachstum(mu_fn, a_end=1.0, a_i=0.01):
    def rhs(x, y):
        a = math.exp(x); D, Dp, eps = y
        muv, deps = mu_fn(a, eps)
        return [Dp, -(2 + dlnE(a)) * Dp + 1.5 * Om_a(a) * muv * D, deps]
    return solve_ivp(rhs, (math.log(a_i), 0.0), [a_i, a_i, 0.0], dense_output=True, rtol=1e-10, atol=1e-14)


def f_of(sol, z):
    x = -math.log(1 + z); D, Dp, _ = sol.sol(x); return Dp / D, D


ref = wachstum(lambda a, e: (1.0, 0.0))


def mu_w(a, e):
    """(1 − w)(ν − 1)/(1 + Q), Argument der Regel = Summe der Verluste s² = θ̃τ_D (θ̃ = 3H), a_N/a₀ vernachlässigt."""
    th = 3 * H0 * Ea(a)
    w = gewicht_w(th)
    s = math.sqrt(th * tau_D())
    return 1 + (1 - w) / math.expm1(s) / (1 + Q), 0.0


def mu_ueb(Nst):
    def f(a, e):
        th = 3 * H0 * Ea(a); w = gewicht_w(th)
        om = H0 * Ea(a) * math.sqrt(Om_a(a) / 2)
        deps = (1 - w) * (1 - e) * om / Nst / (H0 * Ea(a))   # dε/dln a
        return 1 - e * Q / (1 + Q), deps                     # übergebener Anteil klumpt nicht (konservativ)
    return f


P("   Linear (k unterhalb k_ta, skalenfrei): δf/f = f/f_ΛCDM − 1")
P("      z      f_ΛCDM    w-Glied (1−w)(ν−1)    Übergabe N*=0,3     N*=1        N*=2,4      [fσ₈: N*=0,3 / 1 / 2,4]")
sols = {"w": wachstum(mu_w)}
for N in (0.3, 1.0, 2.4):
    sols[N] = wachstum(mu_ueb(N))
for z in (0.0, 0.5, 1.0, 2.0):
    f0, _ = f_of(ref, z)
    row = [f_of(sols[k], z)[0] / f0 - 1 for k in ("w", 0.3, 1.0, 2.4)]
    fs = [f_of(sols[k], z)[0] * f_of(sols[k], z)[1] / (f0 * f_of(ref, z)[1]) - 1 for k in (0.3, 1.0, 2.4)]
    P(f"     {z:3.1f}    {f0:.4f}      {row[0]:+.2e}          {row[1]:+.2e}      {row[2]:+.2e}    {row[3]:+.2e}     [{fs[0]:+.1e} / {fs[1]:+.1e} / {fs[2]:+.1e}]")
P("      (fσ₈ bei gleicher Anfangsamplitude, also σ₈ selbst mitgeändert)")

# Umkehrskala: linearer Kontrast bei der Umkehr (Top-Hat im Hintergrund des Kerns) und σ(R) (BBKS [G], σ₈ = 0,811 [G])


def delta_ta_lin(z_ta):
    """Linearer Kontrast (heute extrapoliert auf a_ta) einer Kugel, die bei z_ta umkehrt."""
    a_ta = 1 / (1 + z_ta); ai = 1e-3
    OL = 1 - Om
    Hb = lambda a: math.sqrt(Om / a**3 + OL)                  # H₀ = 1
    # lineares Wachstum im selben Hintergrund (ohne Strahlung), D'' + (2 + dlnE)D' = 1,5Ω_m(a)D in ln a
    def rhsD(x, y):
        a = math.exp(x); E2 = Om / a**3 + OL
        return [y[1], -(2 - 1.5 * Om / a**3 / E2) * y[1] + 1.5 * Om / a**3 / E2 * y[0]]
    Dl = solve_ivp(rhsD, (math.log(ai), 0.0), [ai, ai], dense_output=True, rtol=1e-10, atol=1e-14)

    def lauf(di):
        Hi = Hb(ai); GM = 0.5 * Hi**2 * (Om / ai**3) / (Om / ai**3 + OL) * (1 + di)
        f = lambda t, y: [y[1], -GM / y[0]**2 + OL * y[0]]
        ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
        s = solve_ivp(f, (0, 5), [1.0, Hi * (1 - di / 3)], events=ev, rtol=1e-10, atol=1e-13)
        return (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (1e3, None)
    from scipy.integrate import quad
    t_of = lambda a: quad(lambda x: 1 / (x * Hb(x)), ai, a)[0]
    tt = t_of(a_ta)
    di = 10**brentq(lambda x: lauf(10**x)[0] - tt, -4, -0.5, xtol=1e-11)
    xa = math.log(a_ta); x0 = math.log(ai)
    R = lauf(di)[1]
    # linear: Start rein wachsend (Geschwindigkeit 1 − δ/3); nichtlinear: 1 + Δ_ta gegen ρ_m(z_ta)
    return di * Dl.sol(xa)[0] / Dl.sol(x0)[0], (1 + di) / (R * ai)**3 * a_ta**3


Gam = Om * h_ * math.exp(-Ob - math.sqrt(2 * h_) * Ob / Om)
kk = np.logspace(-4, 3, 6000)
qq = kk / Gam
T = np.log(1 + 2.34 * qq) / (2.34 * qq) * (1 + 3.89 * qq + (16.1 * qq)**2 + (5.46 * qq)**3 + (6.71 * qq)**4)**-0.25
Pk = kk**0.9649 * T**2


def sig(R):
    x = kk * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(float(np.trapezoid(kk**2 * Pk * W**2, kk)) / (2 * math.pi**2))


nrm = 0.811 / sig(8.0)
Dz = lambda z: f_of(ref, z)[1] / f_of(ref, 0)[1]
DTA = dichte_umkehr_heute()
P()
P(f"   Umkehrskala: 2σ-Bereiche erreichen die Umkehr, 2σ(R_L, z) = δ_ta,lin(z); 1 + Δ_ta (Kern) = {DTA:.3f}")
P("      z    δ_ta,lin   1 + Δ_ta   R_L [Mpc/h]   k_ta = 1/R_L [h/Mpc]   M [M☉/h]")
DTAz = {z: delta_ta_lin(z) for z in (0.0, 1.0, 2.0)}
for z in (0.0, 1.0, 2.0):
    dl, dnl = DTAz[z]
    R1 = brentq(lambda R: 2 * nrm * sig(R) * Dz(z) - dl, 0.05, 60)
    P(f"     {z:3.1f}   {dl:.4f}    {dnl:6.2f}     {R1:7.3f}        {1/R1:.3f}            {4/3*math.pi*R1**3*Om*2.775e11:.2e}")

P()
P("   Jenseits der Umkehrskala (gebunden, Mode übergeben): Kraftverhältnis gegen ΛCDM-Halo μ_b = ν(a_N,b)/(1 + Q), a_N,b = G f_b M/(R_ta/2)²")
P("   und Massenanteil jenseits der Umkehr F_ta = erfc(δ_ta/(√2σ)). μ_b > 1: schnellerer Kollaps; μ_b < 1 ab a_N,b > (ln q)²·a₀ = %.4f a₀" % math.log(q_KOSM)**2)
P("      k [h/Mpc]    z    M [M☉]      a_N,b/a₀     μ_b      F_ta")
fb = 1 / (1 + Q)
for z in (0.0, 1.0, 2.0):
    dl, DTA = DTAz[z]
    for k in (0.1, 0.3, 1.0, 3.0, 10.0):
        RL = 1 / k
        M = 4 / 3 * math.pi * (RL * Mpc / h_)**3 * rho_m_mitbewegt()
        Rta = (3 * M / (4 * math.pi * rho_m_mitbewegt() * (1 + z)**3 * DTA))**(1 / 3)
        aN = G * fb * M / (Rta / 2)**2
        mub = float(nu(aN)) / (1 + Q)
        F = erfc(dl / (math.sqrt(2) * nrm * sig(RL) * Dz(z)))
        P(f"      {k:5.1f}      {z:3.1f}  {M/M_sun:.2e}   {aN/A0:.4f}    {mub:6.2f}    {F:.3f}")
Mx = {}
for z in (0.0, 1.0, 2.0):
    aNc = math.log(q_KOSM)**2 * A0
    # a_N,b = 4 G f_b M^{1/3} (4πρ_ta/3)^{2/3}
    rta = rho_m_mitbewegt() * (1 + z)**3 * DTAz[z][1]
    Mx[z] = (aNc / (4 * G * fb * (4 * math.pi * rta / 3)**(2 / 3)))**3 / M_sun
P("   Masse, ab der μ_b < 1 (a_N,b(R_ta/2) = 0,029 a₀): z = 0: %.1e, z = 1: %.1e, z = 2: %.1e M☉" % (Mx[0.0], Mx[1.0], Mx[2.0]))

# ======================================================================================= 3  Gedächtnis
H("3  Gedächtnis: wo zeigt sich Geschichte?")
v, R = 230e3, 8.2 * kpc
Om_ = v / R
# a_N am Sonnenkreis aus a = ν(a_N)a_N = v²/R
a_obs = v * v / R
aN_s = brentq(lambda x: float(nu(x)) * x - a_obs, 1e-12, 1e-8)
ss = float(s_(aN_s))
K = math.exp(-math.sqrt(Om_ * tau_D() / 2) + ss)
P(f"   Sonnenkreis (v = 230 km/s, R = 8,2 kpc [G]): Ωτ_D = {Om_*tau_D():.0f}, s = {ss:.3f}, |K̂(iΩ)|/K̂(0) ≈ e^(−√(Ωτ_D/2)+s) = {K:.1e}")
for D in (200, 500):
    for z in (0.0, 2.0):
        w_, lo, hi = uebergabezeit(D, z)
        P(f"   Übergabezeit 1/λ_m = N*/ω_orb (Δ = {D}, z = {z}): {lo/Gyr:.2f}–{hi/Gyr:.2f} Mrd. J.")
P("   Preisach: w hängt an θ̃ = max(min θ, 0) der Geschichte: w(θ̃ = 0) = %.3f bleibt nach Wiederausdehnung (θ = 3H heute: w = %.4f)"
  % (gewicht_w(0.0), gewicht_w(3 * H0)))
aNm = 0.1 * A0
P("   Endantwort der Außenbeschleunigung auf δM/M bei festem r: d ln a/d ln M = 1 − χ/2 = %.3f (a_N = 0,1 a₀) … %.3f (0,01 a₀); ohne Nachlauf"
  % (1 - float(chi(aNm)) / 2, 1 - float(chi(0.01 * A0)) / 2))

# ======================================================================================= 4, 5  SPARC
T_ = su.tab_voll()


def punkte(ups_d=0.5, kappa=1.0, sig_gas=0.0):
    rows = []
    for nm, d in T_.items():
        if d['Q'] > 2 or d['inc'] < 30:
            continue
        R, Vo, Vg, Vd, Vb = sf.kurve(nm)
        gas = Vg * np.abs(Vg); st = ups_d * Vd**2 + 0.7 * Vb**2
        Vb2 = kappa * gas + st
        m = (Vb2 > 0) & (Vo > 0)
        if m.sum() < 5:
            continue
        gN = Vb2[m] / R[m]
        go = (Vo[m]**2 + (sig_gas * 1e3)**2 * (R[m] / kpc) / d['Rd']) / R[m]
        r = np.log10(go / (nu(gN) * gN))
        fl = np.clip(gas[m] / Vb2[m], -0.5, 1.5)
        for i in range(m.sum()):
            rows.append((nm, gN[i] / A0, r[i], fl[i], R[m][i] / kpc / d['Rd'], d['fD'], d['T']))
    return rows


H("4  Äquivalenzprinzip an SPARC: hängt der Rest der Regel am Zustand der Quelle (Gas gegen Sterne)?")
for ud in (0.4, 0.5, 0.6):
    for sg in (0.0, 8.0):
        rows = punkte(ud, 1.0, sg)
        nm = np.array([r[0] for r in rows]); x = np.log10(np.array([r[1] for r in rows])); res = np.array([r[2] for r in rows]); fl = np.array([r[3] for r in rows])
        # Galaxie-Nullpunkt abziehen (Abstand, Neigung, Υ): nur der Verlauf innerhalb der Galaxie zählt
        rd = res.copy()
        for g in np.unique(nm):
            k = nm == g; rd[k] -= np.median(res[k])
        rho, p = spearmanr(fl, rd)
        # gleiches g_N: Differenz gasdominiert (f > 0,6) − sterndominiert (f < 0,2) je Bin
        diffs = []
        for lo in (-2.5, -1.5, -1.0, -0.5, 0.0):
            b = (x >= lo) & (x < lo + (1.0 if lo == -2.5 else 0.5))
            gA, sA = b & (fl > 0.6), b & (fl < 0.2)
            if gA.sum() > 15 and sA.sum() > 15:
                diffs.append((lo, np.median(res[gA]) - np.median(res[sA]), gA.sum(), sA.sum()))
        # Steigung innerhalb der Galaxien (Bootstrap über Galaxien)
        gl = np.unique(nm); rng = np.random.default_rng(1); sl = []
        idx = {g: np.where(nm == g)[0] for g in gl}
        for _ in range(300):
            pick = np.concatenate([idx[g] for g in rng.choice(gl, len(gl))])
            A = fl[pick]; B = rd[pick]; sl.append(np.polyfit(A - A.mean(), B, 1)[0])
        s0 = np.polyfit(fl - fl.mean(), rd, 1)[0]
        P(f"   Υ_d = {ud}, σ_Gas = {sg:.0f} km/s: {len(gl)} Galaxien, {len(rd)} Punkte; innerhalb der Galaxien d(Rest)/d(f_Gas,örtl) = {s0:+.3f} ± {np.std(sl):.3f} dex, ρ = {rho:+.3f} (p = {p:.1e})")
        P("        gleiches g_N (Bin log g_N/a₀ ab): Rest(Gas > 0,6) − Rest(Sterne < 0,2) = " +
          ", ".join(f"[{lo:+.1f}] {d:+.3f} ({a}/{b})" for lo, d, a, b in diffs))


def chi2_kappa(kap, ud=0.5, sg=0.0):
    rows = punkte(ud, kap, sg)
    nm = np.array([r[0] for r in rows]); res = np.array([r[2] for r in rows])
    tot = 0.0
    for g in np.unique(nm):
        k = nm == g; tot += np.sum((res[k] - np.median(res[k]))**2)
    return tot


P("   Kopplung des Gases in g_N (Regel universell ⇒ κ = 1), Form der Kurven je Galaxie (Nullpunkt frei):")
for ud in (0.4, 0.5, 0.6):
    for sg in (0.0, 8.0):
        r_ = minimize_scalar(lambda k: chi2_kappa(k, ud, sg), bounds=(0.3, 2.5), method='bounded', options={'xatol': 1e-3})
        c1 = chi2_kappa(1.0, ud, sg)
        P(f"      Υ_d = {ud}, σ = {sg:.0f}: κ̂ = {r_.x:.2f}  (Σr² {r_.fun:.2f} gegen {c1:.2f} bei κ = 1)")

H("5  Umgebung und Geschichte an SPARC")
o = su.reste()
names = list(o)
uma = np.array([T_[n]['fD'] == 4 for n in names])
rest = np.array([o[n]['rest'] for n in names]); dA = np.array([o[n]['aussen'] - o[n]['innen'] for n in names])
P(f"   Ursa-Major-Haufen (f_D = 4, gleiche Entfernung): n = {uma.sum()}, übrige n = {(~uma).sum()}")
P(f"      Rest je Galaxie Median: UMa {np.median(rest[uma]):+.3f}, übrige {np.median(rest[~uma]):+.3f} dex (Mann-Whitney p = {mannwhitneyu(rest[uma], rest[~uma]).pvalue:.2f})")
m = np.isfinite(dA)
P(f"      Form außen − innen: UMa {np.nanmedian(dA[uma]):+.3f}, übrige {np.nanmedian(dA[~uma]):+.3f} dex (p = {mannwhitneyu(dA[uma & m], dA[~uma & m]).pvalue:.2f})")
rows = punkte(0.5)
nmA = np.array([r[0] for r in rows]); xA = np.array([r[1] for r in rows]); rA = np.array([r[2] for r in rows]); fDA = np.array([r[5] for r in rows])
tief = xA < 0.1
ru, ro = rA[tief & (fDA == 4)], rA[tief & (fDA != 4)]
P(f"      Punkte mit g_N < 0,1 a₀: Rest UMa {np.median(ru):+.3f} ({ru.size}), übrige {np.median(ro):+.3f} ({ro.size}) dex")
# welches senkrechte äußere Feld gäbe diese Senkung an denselben UMa-Punkten?
uv = xA[tief & (fDA == 4)] * A0
for ge in (0.01, 0.03, 0.1):
    dd = np.median([math.log10(nu_aussenfeld(g, ge * A0, "senkrecht") / float(nu(g))) for g in uv])
    P(f"      äußeres Feld senkrecht g_e = {ge} a₀ senkte den Rest dieser UMa-Punkte um {dd:+.3f} dex")
P("   Form (außen − innen) und Rest gegen Geschichte (Spearman):")
for k, lab in (("T", "Hubble-Typ"), ("fgas", "Gasanteil"), ("lSB", "log Σ₀"), ("lMb", "log M_b")):
    xx = np.array([o[n][k] for n in names])
    r1, p1 = spearmanr(xx[m], dA[m]); r2, p2 = spearmanr(xx, rest)
    P(f"      {lab:12s}: außen−innen ρ = {r1:+.3f} (p {p1:.1e});  Rest ρ = {r2:+.3f} (p {p2:.1e})")

# ======================================================================================= 6  Handwerte
H("6  Handwerte (von Hand, nicht mit diesem Code)")
P("   Kuzmin h = 3r*, R = h: g_R = 1/(9·2^{3/2}) = 0,03928, |g| = 1/18 = 0,05556;")
P("   ν(0,03928): s = 0,19819, 1 − e^{−s} = 0,17979, ν = 5,562; ν(0,05556): s = 0,23570, 1 − e^{−s} = 0,20995, ν = 4,763; η = 5,562/4,763 = 1,168")
P("   Kern am Sonnenkreis: Ω = 2,3e5/(8,2·3,0857e19) = 9,09e-16 s⁻¹, τ_D = c/a₀ = 2,725e18 s, Ωτ_D/2 = 1239, √ = 35,2, s ≈ 1,1 ⇒ e^{−34,1} = 1,5e-15")
P(f"   Code: η(Kuzmin, h = 3, R = h) = {float(nu1(1/(9*2**1.5)))/float(nu1(1/18)):.4f}")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "trenntests_b_dm.txt"), "w", encoding="utf-8").write("\n".join(OUT) + "\n")
