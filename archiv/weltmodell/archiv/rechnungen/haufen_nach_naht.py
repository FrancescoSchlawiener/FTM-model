# -*- coding: utf-8 -*-
"""
Haufen nach der Naht: ohne Budget Q·M_b je System.

Gebundener Haufen (w = 0 nach der Umkehr): örtliche Regel ν(a_N)·M_b plus Mode ξ nur, soweit sie je Schale noch
nicht übergeben ist, ξ_R = e^{−I/N*}, I = ∫(1 − w)ω dt (Schale mit w aus ftm_kern.gewicht_w, ω² = GM_L/R³, nach
R_ta/2 fest; Verfahren wie naht_rq_umkehr). Keine Subtraktion Q − (ν − 1), keine Kappe.
  (B) Übergang wie naht_rq_umkehr (ii):   M_gel = M_b + (1 − ξ_R)(ν − 1)M_b + ξ_R·M_ξ
  (A) additiv (obere Grenze):              M_gel = ν·M_b + ξ_R·M_ξ
  Mode vor der Übergabe wie druckfreier Staub mit der Bahngeschichte der Baryonen (Form NFW, stoßfrei), zwei Normierungen:
    (i) Q·M_b(<r200)·f_NFW(r)   (ii) q·M_ref(r)
Haufen wie schliessungsregel / Prüfstand: β-Gas 10¹⁴ M☉ in 1 Mpc (r_c 150 kpc, bis 3 Mpc), BCG 1,5·10¹² (20 kpc),
Galaxien 10¹³ (300 kpc); Linsenmasse NFW 10¹⁵, c = 4 bei z = 0,3 (gelesen, [G]).
Dazu: Kac-Mittel ⟨ν⟩_b je System (G), M_Q in M_b mit R₂₀₀ aus der gelesenen Masse, Prüfung von Hand.
Nur ftm_kern (Formeln, Konstanten). Aufruf: python haufen_nach_naht.py > haufen_nach_naht.txt
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

G, Ms, kpc, Mpc, Gyr = K.G, K.M_sun, K.kpc, K.Mpc, K.Gyr
a0, Q, qk = K.A0, K.Q_KOSM, K.q_KOSM
H0, Om = K.H0_si(), K.OMEGA_M
OL = 1 - Om
rho_c0 = 3 * H0**2 / (8 * math.pi * G)
rho_m0 = Om * rho_c0
NS = (0.3, 1.0, 2.4)
RK = (100, 250, 500, 1000, 1500)


def hdr(t):
    print("\n" + "=" * 110 + "\n" + t + "\n" + "=" * 110)


# ------------------------------------------------------------------ Haufen
zl = 0.3
rho_cz = rho_c0 * K.E_hintergrund(zl)**2
f_g = lambda x: x * x / (1 + (x / (150 * kpc))**2)
N_g = quad(f_g, 0, Mpc)[0]


def Mb_cl(r):
    Mg = 1e14 * Ms * quad(f_g, 0, min(r, 3 * Mpc))[0] / N_g
    return Mg + K.hernquist_masse(1.5e12 * Ms, 20 * kpc, r) + K.hernquist_masse(1e13 * Ms, 300 * kpc, r)


M200 = 1e15 * Ms
r200 = (3 * M200 / (4 * math.pi * 200 * rho_cz)) ** (1 / 3)
mu = lambda x: math.log(1 + x) - x / (1 + x)
Mref = lambda r: M200 * mu(4 * r / r200) / mu(4)
Mb200 = Mb_cl(r200)
nu_r = lambda r: float(K.nu(G * Mb_cl(r) / r**2, a0))
Mxi = {"i": lambda r: Q * Mb200 * Mref(r) / M200, "ii": lambda r: qk * Mref(r)}

# ------------------------------------------------------------------ Schale (wie naht_rq_umkehr.lauf)
Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_von_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
AI = 1e-3; TI = t_von_a(AI); RHO_I = rho_m0 / AI**3
t_obs = t_von_a(1 / (1 + zl))


def lauf(M, di, t_end):
    ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
    vi = Hb(AI) * ri * (1 - di / 3)

    def f(t, y):
        R, V, I = y
        w = K.gewicht_w(3 * V / R, a0)
        return [V, -G * M / R**2 + OL * H0**2 * R, (1 - w) * math.sqrt(G * M / R**3)]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (TI, t_end), [ri, vi, 0.0], events=ev, rtol=1e-8, atol=1e-6)
    if not len(s.t_events[0]):
        return s.y[0][-1], s.y[2][-1], None, None
    t1, (R1, V1, I1) = s.t_events[0][0], s.y_events[0][0]
    ev2 = lambda t, y: y[0] - R1 / 2; ev2.terminal = True; ev2.direction = -1
    s2 = solve_ivp(f, (t1, t_end), [R1, -1e-9, I1], events=ev2, rtol=1e-8, atol=1e-6)
    if not len(s2.t_events[0]):
        return s2.y[0][-1], s2.y[2][-1], t1, None
    t2, I2 = s2.t_events[0][0], s2.y_events[0][0][2]
    return R1 / 2, I2 + (t_end - t2) * math.sqrt(G * M / (R1 / 2)**3), t1, t2


def schale(r):
    """Schale mit M_L = M_ref(<r), die bei z = 0,3 den Radius r hat (virialisiert: R_ta/2 = r)."""
    M = Mref(r)
    fn = lambda ld: lauf(M, 10**ld, t_obs)[0] - r
    ld = brentq(fn, -3, 1.5, xtol=1e-10)
    R, I, t1, t2 = lauf(M, 10**ld, t_obs)
    z_ = lambda t: 1 / brentq(lambda a: t_von_a(a) - t, 1e-4, 1.0) - 1 if t else float("nan")
    return I, z_(t1), z_(t2), math.sqrt(G * M / r**3)


hdr("0  Haufen, Schalen und Übergabe (z = 0,3, a₀ = %.2e)" % a0)
print(f"   r200 = {r200/kpc:.0f} kpc, M_b(<r200) = {Mb200/Ms:.3e} M☉ (Anteil {Mb200/M200:.3f}, kosmisch {1/(1+Q):.3f})")
print("   r[kpc]   M_b/M☉     a_N/a₀    ν      M_ref/(ν M_b)   z_Umkehr  z_R_ta/2   ω[1/Gyr]   I      ξ_R N* = 0,3 / 1 / 2,4")
S = {}
for rk in RK:
    r = rk * kpc
    mb = Mb_cl(r)
    I, zta, zc, om = schale(r)
    xi = [math.exp(-I / n) for n in NS]
    S[rk] = (mb, nu_r(r), I, xi)
    print(f"   {rk:5d}  {mb/Ms:9.3e}  {G*mb/r**2/a0:7.4f}  {nu_r(r):6.3f}   {Mref(r)/(nu_r(r)*mb):6.3f}         {zta:5.2f}     {zc:5.2f}      {om*Gyr:6.3f}   {I:6.3f}   "
          + " / ".join(f"{x:.3f}" for x in xi))

hdr("1  Fehlender Faktor M_ref/M_gel und Anteil der Regel am gelesenen Zusatz")
for form in ("B", "A"):
    for nrm in ("i", "ii"):
        print(f"   ({form}) Normierung ({nrm}):")
        print("     r[kpc]   ξ_R nötig |  M_ref/M_gel  N* = 0,3 / 1 / 2,4   |  Anteil Regel am Zusatz (M_gel − M_b)  N* = 0,3 / 1 / 2,4")
        for rk in RK:
            r = rk * kpc
            mb, n_, I, xi = S[rk]
            mx = Mxi[nrm](r)
            if form == "B":
                gel = lambda x: mb + (1 - x) * (n_ - 1) * mb + x * mx
                reg = lambda x: (1 - x) * (n_ - 1) * mb
                need = (Mref(r) - n_ * mb) / (mx - (n_ - 1) * mb)
            else:
                gel = lambda x: n_ * mb + x * mx
                reg = lambda x: (n_ - 1) * mb
                need = (Mref(r) - n_ * mb) / mx
            fak = [Mref(r) / gel(x) for x in xi]
            ant = [reg(x) / (gel(x) - mb) for x in xi]
            print(f"     {rk:5d}    {need:6.3f}   |  " + " / ".join(f"{v:5.2f}" for v in fak)
                  + "              |  " + " / ".join(f"{v:5.2f}" for v in ant))
print("   ξ_R nötig: Rest der Mode, der die Linsenmasse genau trifft (> 1: auch ohne jede Übergabe zu wenig).")
print("   Alt (Budget Q − (ν − 1)): Operatoranteil 1 − W_rest bei 100 / 1500 kpc = 0,21 / 0,71 ((ii), a₀ = 1,2·10⁻¹⁰).")

# Bullet: Anteil, der den Galaxien folgt (unübergebene Mode) gegen Regel (folgt dem Gas mit Nachlauf)
mb, n_, I, xi = S[100]; mx = Mxi["ii"](100 * kpc)
fr = [x * mx / ((1 - x) * (n_ - 1) * mb + x * mx) for x in xi]
print("   Bullet (100 kpc, (B)(ii)): Anteil der unübergebenen Mode am Zusatz (folgt den Galaxien) = "
      + " / ".join(f"{abs(v):.3f}" for v in fr) + "  (N* = 0,3 / 1 / 2,4; alt 0,79)")

# Robustheit: kleinste mögliche Übergabe bei 100 kpc
r = 100 * kpc
om_nu = math.sqrt(G * S[100][1] * S[100][0] / r**3)
for zc in (0.5, 0.4):
    dt = t_obs - t_von_a(1 / (1 + zc))
    Imin = om_nu * dt
    print(f"   Robust, 100 kpc: ω aus ν·M_b ({om_nu*Gyr:.2f}/Gyr), Kern erst bei z = {zc} gebunden (Δt = {dt/Gyr:.2f} Gyr): I ≥ {Imin:.2f},"
          f" ξ_R ≤ {math.exp(-Imin/2.4):.3f} (N* = 2,4); ξ_R = 0,67 verlangte N* ≥ {Imin/-math.log(0.665):.0f}")

hdr("2  Kac-Mittel ⟨ν⟩_b = ∫ν(a_N(r)) dM_b / M_b je System (G), gegen 1 + Q = %.3f" % (1 + Q))
Mw = 6e10 * Ms


def nu_mittel(Mb_fn, rmax, rmin=1e-3 * kpc, n=4000):
    rr = np.geomspace(rmin, rmax, n)
    mm = np.array([Mb_fn(x) for x in rr])
    nn = np.array([float(K.nu(G * m / x**2, a0)) for m, x in zip(mm, rr)])
    dm = np.diff(mm)
    return float(np.sum(0.5 * (nn[1:] + nn[:-1]) * dm) / (mm[-1] - mm[0]))


for a_h in (3, 10):
    print(f"   Milchstraße Hernquist 6·10¹⁰ M☉, a = {a_h} kpc, bis 100 kpc: ⟨ν⟩_b = {nu_mittel(lambda r: K.hernquist_masse(Mw, a_h * kpc, r), 100 * kpc):.3f}")
print(f"   Haufen bis r200: ⟨ν⟩_b = {nu_mittel(Mb_cl, r200):.3f};  bis 3 Mpc: {nu_mittel(Mb_cl, 3 * Mpc):.3f}")
print(f"   Haufen gelesen in r200 / M_b(<r200) = {M200/Mb200:.3f};  ν(r200)·M_b(<r200)/M_b = ν(r200) = {nu_r(r200):.3f}")

hdr("3  M_Q: R_Q = R₂₀₀ mit der gelesenen Masse ν(R₂₀₀)·M_b (Punktmasse)")
lq = -math.log(qk)
print(f"   ν an R_Q = 1/(1 − q) = {float(K.nu(lq**2 * a0, a0)):.4f} = 1 + Q: gelesen an R_Q ist (1 + Q)·M_b nach der Regel.")
for zz in (0.0, 0.3, 1.0):
    Hz = H0 * K.E_hintergrund(zz)
    MbQ = (1 + Q)**2 * lq**6 * a0**3 / (1e4 * G * Hz**4)
    RQ = K.r_stern(MbQ, a0) / lq
    R200 = (G * (1 + Q) * MbQ / (100 * Hz**2)) ** (1 / 3)
    print(f"   z = {zz}: M_b,Q = (1 + Q)²(ln q)⁶a₀³/(10⁴GH⁴) = {MbQ/Ms:.3e} M☉, gelesen (1 + Q)M_b,Q = {(1+Q)*MbQ/Ms:.3e} M☉;"
          f" R_Q = {RQ/kpc:.0f} kpc, R₂₀₀ = {R200/kpc:.0f} kpc")

hdr("4  Von Hand (Zahlen eingetippt, ohne ftm_kern)")
# 100 kpc: Gas β: ∫x²/(1+x²/rc²) = rc³(u − atan u), u = r/rc
rc = 150.0
gas = lambda rk: rc**3 * (min(rk, 3000) / rc - math.atan(min(rk, 3000) / rc))
Mg = 1e14 * gas(100) / gas(1000)
Mbh = Mg + 1.5e12 * 100**2 / 120**2 + 1e13 * 100**2 / 400**2
GMsun = 1.32712440018e20
aN = GMsun * Mbh / (100 * 3.0856775814913673e19)**2
s = math.sqrt(aN / 1.1e-10)
nuh = 1 / (1 - math.exp(-s))
print(f"   100 kpc: M_gas = {Mg:.4e}, BCG {1.5e12*100**2/120**2:.4e}, Galaxien {1e13*100**2/400**2:.4e}, M_b = {Mbh:.4e} M☉ (Skript {S[100][0]/Ms:.4e})")
print(f"   a_N/a₀ = {aN/1.1e-10:.4f}, s = {s:.4f}, ν = {nuh:.4f} (Skript {S[100][1]:.4f})")
# NFW: ρ_c(z = 0,3) = 3H²/(8πG), H = 67,4·E, E² = 0,3134·1,3³ + 0,6865(+Ω_r)
E2 = 0.3134 * 1.3**3 + 0.6865 + 9.14e-5 * 1.3**4
Hh = 67.4e3 / 3.0856775814913673e22 * math.sqrt(E2)
rhoc = 3 * Hh**2 / (8 * math.pi * 6.6743e-11)
r2 = (3 * 1e15 * 1.98841e30 / (4 * math.pi * 200 * rhoc)) ** (1 / 3) / 3.0856775814913673e19
x = 4 * 100 / r2
Mr = 1e15 * (math.log(1 + x) - x / (1 + x)) / (math.log(5) - 0.8)
print(f"   r200 = {r2:.0f} kpc, M_ref(100 kpc) = {Mr:.4e} M☉;  Faktor M_ref/(ν M_b) = {Mr/(nuh*Mbh):.3f} (Skript {Mref(100*kpc)/(S[100][1]*S[100][0]):.3f})")
GMr = 1.32712440018e20 * Mr
om_h = math.sqrt(GMr / (100 * 3.0856775814913673e19)**3) * 3.15576e16
print(f"   ω(100 kpc) = √(G·M_ref/r³) = {om_h:.3f}/Gyr (Skript {S[100][2] and math.sqrt(G*Mref(100*kpc)/(100*kpc)**3)*Gyr:.3f})")
print(f"   (B)(ii) ohne Übergabe (ξ_R = 1): M_gel = M_b + q·M_ref = {Mbh + 0.84287*Mr:.4e}, M_ref/M_gel = {Mr/(Mbh+0.84287*Mr):.3f}")
