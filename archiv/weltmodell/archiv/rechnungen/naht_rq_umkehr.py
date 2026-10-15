# -*- coding: utf-8 -*-
"""
Naht zwischen R_Q und dem Umkehrradius R_ta: Wie geht die örtliche Zählung ν(a_N) (L) in die Mitnahme mit 1 + Q (G) über?

Milchstraße als Punktmasse M_b = 6·10¹⁰ M☉ (Kernwerte), gelesene Masse M_gel(<r) = r²a(r)/G, 10 kpc – 2 Mpc.
  (i)   örtliche Regel bis R_ta, außen Mitnahme (1 + Q)·M_b
  (ii)  Rest der Mode ξ_R je Schale nach der Übergabe λ = (1 − w)ω/N* (w aus ftm_kern.gewicht_w, ω² = GM_L/R³,
        nach r_ta/2 fest; Verfahren wie linsen_isoliert), N* = 0,3 / 1 / 2,4, in zwei Formen:
          je Schale (L):  M_gel = M_b[1 + (1 − ξ_R)(ν − 1) + ξ_R·Q·m(r)]
          Gaußfluss:      dM_gel/dr = (1 − ξ_R)·M_b dν/dr + ξ_R·Q·M_b·dm/dr
        m(r) = Anteil der Mode innerhalb r. Lagrange-Masse der Schalen zwei Wege:
          Kac:   M_L = M_b(1 + Q·min(r/R_ta, 1)) + Hintergrund jenseits R_ta   (Mode Q·M_b bis R_ta, gleichmäßig in r)
          flach: M_L = v_f² r/G  (Lesung der flachen Kurve, wie linsen_isoliert)
  (iii) Kac-Kappe ν ≤ 1 + Q jenseits R_Q.
Bilanz in R_ta, Kappe gegen SPARC (nur weltmodell/daten/sparc; Auswahl wie sparc_aussenfeld_rand), Prüfung von Hand.
"""
import io, sys, os, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K
import sparc_flaechendichte as sf, sparc_universal as su

G, c, Ms, kpc, Mpc, Gyr = K.G, K.c, K.M_sun, K.kpc, K.Mpc, K.Gyr
a0, Q, qk = K.A0, K.Q_KOSM, K.q_KOSM
H0, Om = K.H0_si(), K.OMEGA_M
OL = 1 - Om
Mb = 6e10 * Ms
rho_c = 3 * H0**2 / (8 * math.pi * G)
rho_m = Om * rho_c
rng = np.random.default_rng(2909)

def hdr(t):
    print("\n" + "=" * 110 + "\n" + t + "\n" + "=" * 110)

# ------------------------------------------------------------------ 0 Grundgrößen
Hb = lambda a: H0 * math.sqrt(Om / a**3 + OL)
t_von_a = lambda a: quad(lambda x: 1 / (x * Hb(x)), 1e-8, a, limit=200)[0]
t0 = t_von_a(1.0)
AI = 1e-3; TI = t_von_a(AI); RHO_I = rho_m / AI**3

def lauf(M, di, t_obs, voll=True):
    """Lagrange-Schale der Masse M (Anfangsüberdichte di bei a = 10⁻³); gibt R(t_obs), I = ∫(1 − w)ω dt, Umkehr (t, R)."""
    ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
    vi = Hb(AI) * ri * (1 - di / 3)
    def f(t, y):
        R, V, I = y
        w = K.gewicht_w(3 * V / R, a0)
        return [V, -G * M / R**2 + OL * H0**2 * R, (1 - w) * math.sqrt(G * M / R**3)]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (TI, t_obs), [ri, vi, 0.0], events=ev, rtol=1e-8, atol=1e-6)
    if not len(s.t_events[0]):
        return s.y[0][-1], s.y[2][-1], None
    t1, (R1, V1, I1) = s.t_events[0][0], s.y_events[0][0]
    ev2 = lambda t, y: y[0] - R1 / 2; ev2.terminal = True; ev2.direction = -1
    s2 = solve_ivp(f, (t1, t_obs), [R1, -1e-9, I1], events=ev2, rtol=1e-8, atol=1e-6)
    if not len(s2.t_events[0]):
        return s2.y[0][-1], s2.y[2][-1], (t1, R1)
    t2, I2 = s2.t_events[0][0], s2.y_events[0][0][2]
    return R1 / 2, I2 + (t_obs - t2) * math.sqrt(G * M / (R1 / 2)**3), (t1, R1)

# Umkehrdichte heute: Schale, die genau heute umkehrt
def umkehr_heute(M):
    fn = lambda ld: (lambda r: (r[2][0] if r[2] else 2 * t0) - t0)(lauf(M, 10**ld, 2 * t0))
    ld = brentq(fn, -4, 0, xtol=1e-9)
    R, I, (tt, Rt) = lauf(M, 10**ld, 2 * t0)
    _, I_ta, _ = lauf(M, 10**ld, tt * (1 - 1e-9))
    return Rt, 3 * M / (4 * math.pi * Rt**3) / rho_m, I_ta

M_mit = (1 + Q) * Mb
R_ta, D_ta, I_ta = umkehr_heute(M_mit)
rs = K.r_stern(Mb, a0)
R_Q = rs / (-math.log(qk))
sQ = math.sqrt(G * Mb / R_Q**2 / a0)
vf = K.v_flach(Mb, a0)

hdr("0  Grundgrößen")
print(f"   M_b = 6·10¹⁰ M☉, a₀ = {a0:.2e}, Q = {Q:.3f}, q_kosm = {qk:.4f}, 1 + Q = {1+Q:.3f};  t₀ = {t0/Gyr:.3f} Mrd. J.")
print(f"   r* = {rs/kpc:.2f} kpc, R_Q = r*/(−ln q) = {R_Q/kpc:.2f} kpc [L]; a_N/a₀ dort = {sQ**2:.4f}, ν = {K.nu(sQ**2*a0):.4f}, χ = {K.chi(sQ**2*a0):.4f}")
print(f"   Umkehr heute, Schale (1 + Q)M_b = {M_mit/Ms:.3e} M☉ (Top-Hat mit Λ): 1 + Δ_ta = {D_ta:.3f}, R_ta = {R_ta/Mpc:.4f} Mpc [L, gesetzt durch w: G-Zustand]")
print(f"   I_vor = ∫(1 − w)ω dt bis zur Umkehr = {I_ta:.4f}  →  ξ_R(Umkehr) = e^(−I/N*) = " +
      ", ".join(f"{math.exp(-I_ta/n):.3f} (N* = {n})" for n in (0.3, 1.0, 2.4)))
x_ta = G * Mb / R_ta**2 / a0
print(f"   an R_ta: a_N/a₀ = {x_ta:.4e}, ν = {K.nu(x_ta*a0):.2f}; v_f = {vf/1e3:.1f} km/s")

# ------------------------------------------------------------------ Rest der Mode je Schale
def M_L_kac(r):
    if r <= R_ta:
        return Mb * (1 + Q * r / R_ta)
    return M_mit + 4 / 3 * math.pi * rho_m * (r**3 - R_ta**3)
M_L_flach = lambda r: vf**2 * r / G
PROF = {"Kac": M_L_kac, "flach": M_L_flach}

def I_bei(r, Mf):
    M = Mf(r)
    fn = lambda ld: lauf(M, 10**ld, t0)[0] - r
    if fn(-6.0) < 0:
        return lauf(M, 1e-6, t0)[1]
    ld = brentq(fn, -6.0, 0.3, xtol=1e-5)
    return lauf(M, 10**ld, t0)[1]

R_XI = np.geomspace(20, 2500, 24) * kpc
I_TAB = {p: np.array([I_bei(r, PROF[p]) for r in R_XI]) for p in PROF}

def xi_R(r, prof, Ns):
    lI = np.interp(np.log(r), np.log(R_XI), np.log(np.maximum(I_TAB[prof], 1e-300)))
    return np.exp(-np.exp(lI) / Ns)

# ------------------------------------------------------------------ 1 Tabelle M_gel
RR = np.geomspace(10 * kpc, 2000 * kpc, 4000)
xN = G * Mb / RR**2 / a0
NU = K.nu(xN * a0)
m_mode = np.minimum(RR / R_ta, 1.0)

def M_i(r, nu_):
    return np.where(r <= R_ta, nu_, 1 + Q)

def M_ii_lokal(r, nu_, xi):
    return 1 + (1 - xi) * (nu_ - 1) + xi * Q * np.minimum(r / R_ta, 1.0)

def M_ii_gauss(r, nu_, xi):
    dnu = np.gradient(nu_, r)
    dm = np.where(r < R_ta, 1 / R_ta, 0.0)
    d = (1 - xi) * dnu + xi * Q * dm
    start = 1 + (1 - xi[0]) * (nu_[0] - 1) + xi[0] * Q * r[0] / R_ta
    return start + np.concatenate([[0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(r))])

def M_iii(r, nu_):
    return np.where(r <= R_Q, nu_, np.minimum(nu_, 1 + Q))

hdr("1  Gelesene Masse M_gel(<r)/M_b = r²a/(G M_b)   (Milchstraße, Punktmasse)")
print(f"   Radien: R_Q = {R_Q/kpc:.1f} kpc, R_ta = {R_ta/kpc:.0f} kpc.  ξ_R je Schale (Lagrange-Masse Kac | flach):")
print("   r[kpc]  " + "  ".join(f"{int(r_):>5d}" for r_ in (30, 51, 100, 150, 200, 300, 400, 500, 580, 700, 1000, 1500, 2000)))
for p in PROF:
    for Ns in (0.3, 1.0, 2.4):
        print(f"   ξ_R {p:5s} N* = {Ns:3.1f}: " + " ".join(f"{float(xi_R(r_*kpc, p, Ns)):6.3f}" for r_ in (30, 51, 100, 150, 200, 300, 400, 500, 580, 700, 1000, 1500, 2000)))

tab_r = np.array([10, 20, 30, R_Q / kpc, 75, 100, 150, 200, 300, 400, 500, 0.999 * R_ta / kpc, 1.001 * R_ta / kpc, 700, 1000, 1500, 2000])
def at(arr, r_):
    return np.interp(np.log(r_ * kpc), np.log(RR), arr)

VAR = [("(i) Regel/Mitnahme", M_i(RR, NU)), ("(iii) Kappe", M_iii(RR, NU))]
for p in PROF:
    for Ns in (0.3, 1.0, 2.4):
        xi = xi_R(RR, p, Ns)
        VAR.append((f"(ii) {p} {Ns} Schale", M_ii_lokal(RR, NU, xi)))
        VAR.append((f"(ii) {p} {Ns} Gauß", M_ii_gauss(RR, NU, xi)))
print("\n   " + f"{'r [kpc]':24s}" + "".join(f"{r_:7.0f}" for r_ in tab_r))
print("   " + f"{'a_N/a₀':24s}" + "".join(f"{G*Mb/(r_*kpc)**2/a0:7.4f}" if r_ < 200 else f" {G*Mb/(r_*kpc)**2/a0:6.0e}" for r_ in tab_r))
print("   " + f"{'ν(a_N)':24s}" + "".join(f"{float(K.nu(G*Mb/(r_*kpc)**2)):7.2f}" for r_ in tab_r))
EXAKT = {"(i) Regel/Mitnahme": lambda r_: float(M_i(r_ * kpc, K.nu(G * Mb / (r_ * kpc)**2))),
         "(iii) Kappe": lambda r_: float(M_iii(r_ * kpc, K.nu(G * Mb / (r_ * kpc)**2)))}
for lab, arr in VAR:
    f_ = EXAKT.get(lab, lambda r_: at(arr, r_))
    print("   " + f"{lab:24s}" + "".join(f"{f_(r_):7.2f}" for r_ in tab_r))
print("   (Spalten 12/13: unmittelbar innen/außen an R_ta)")

print("\n   Stetigkeit: Sprung von M_gel an R_ta (außen/innen), Maximum und größter Abfall d ln M_gel/d ln r:")
for lab, arr in VAR:
    k = np.searchsorted(RR, R_ta)
    j = arr[k + 1] / arr[k - 1]
    dl = np.gradient(np.log(arr), np.log(RR))
    kmax = np.argmax(arr)
    neg = RR[dl < -1e-3]
    zone = f"{neg.min()/kpc:.0f}–{neg.max()/kpc:.0f} kpc" if neg.size else "keiner"
    print(f"   {lab:24s}: außen/innen an R_ta = {j:6.3f};  Maximum {arr[kmax]:6.2f} bei {RR[kmax]/kpc:6.0f} kpc;  "
          f"min d ln M/d ln r = {dl.min():+8.2f};  M_gel fällt (Gauß: negative gelesene Dichte) in {zone};  M_gel(2 Mpc) = {arr[-1]:.2f}")

# ------------------------------------------------------------------ 2 Bilanz
hdr("2  Bilanz in R_ta: gelesene Masse nach der Regel gegen (1 + Q)·M_b (Kac, G)")
nu_ta = float(K.nu(x_ta * a0)); nu_Q = float(K.nu(sQ**2 * a0))
I_X = quad(lambda r: float(np.gradient([K.nu(G*Mb/(r*(1-1e-6))**2), K.nu(G*Mb/(r*(1+1e-6))**2)], [r*(1-1e-6), r*(1+1e-6)])[0]), R_Q, R_ta, limit=200)[0]
print(f"   M_gel,Regel(<R_ta) = ν(R_ta)·M_b = {nu_ta:.2f}·M_b = {nu_ta*Mb/Ms:.3e} M☉;  (1 + Q)·M_b = {1+Q:.3f}·M_b = {M_mit/Ms:.3e} M☉;  Verhältnis {nu_ta/(1+Q):.2f}")
print(f"   Zusatz der Regel (ν − 1)M_b: an R_Q {nu_Q-1:.3f}·M_b = Q (Bilanz exakt), an R_ta {nu_ta-1:.2f}·M_b = {(nu_ta-1)/Q:.2f}·Q")
print(f"   Integral ∫_{{R_Q}}^{{R_ta}} M_b dν/dr dr = {I_X:.2f}·M_b (Kontrolle ν(R_ta) − ν(R_Q) = {nu_ta-nu_Q:.2f}); Mode jenseits R_Q (Kac-Profil) = Q(1 − R_Q/R_ta)·M_b = {Q*(1-R_Q/R_ta):.2f}·M_b")
print(f"   mittlere gelesene Dichte in R_ta nach der Regel: {nu_ta/(1+Q)*D_ta:.1f}·ρ̄_m (Kac/Umkehr: {D_ta:.2f}·ρ̄_m)")
# Umkehr, wenn die Schale an R_ta die Regel der inneren Masse läse (flache Lagrange-Masse)
R_ta_fl = math.sqrt(3 * vf**2 / (4 * math.pi * G * D_ta * rho_m))
print(f"   Läse die Schale an der Umkehr die innere Regel (M_L = v_f²r/G), läge die Umkehr heute bei {R_ta_fl/Mpc:.2f} Mpc statt {R_ta/Mpc:.2f} Mpc "
      f"(ν dort {float(K.nu(G*Mb/R_ta_fl**2)):.0f})")
# Kac-Mittel über die Lagrange-Masse der Mode: mittlere Zählung der Schalen
m_sh = lambda r: Q * Mb / R_ta                                                # dM_ξ/dr
nu_mittel = quad(lambda r: float(K.nu(G*Mb/r**2)) * m_sh(r), 1e-3 * kpc, R_ta, limit=400)[0] / (Q * Mb)
nu_mittel_Q = quad(lambda r: float(K.nu(G*Mb/r**2)) * m_sh(r), 1e-3 * kpc, R_Q, limit=400)[0] / (Q * Mb * R_Q / R_ta)
print(f"   ν gemittelt über die Mode (Kac-Profil) in R_ta: {nu_mittel:.2f}; nur bis R_Q: {nu_mittel_Q:.2f} (gegen 1 + Q = {1+Q:.2f})")
# Punktmasse: Mode-Halo mit v_f (isotherm) endet dort, wo er Q·M_b trägt
print(f"   Isothermer Halo mit der Tully-Fisher-Geschwindigkeit trägt Q·M_b bis r = Q·r* = {Q*rs/kpc:.1f} kpc; Regel: (ν − 1) = Q bei R_Q = {R_Q/kpc:.1f} kpc "
      f"(tief: ν − 1 ≈ r/r* − ½ ⇒ R_Q ≈ (Q + ½)r* = {(Q+0.5)*rs/kpc:.1f} kpc)")

# ------------------------------------------------------------------ 3 SPARC: Kappe jenseits R_Q
hdr("3  SPARC: Kappe ν ≤ 1 + Q jenseits R_Q (a_N < s_kosm²·a₀ = 0,0292·a₀) gegen die Regel ohne Kappe")
S_K2 = math.log(qk)**2
tab = su.tab_voll()
ALLE = []
for nm, d in tab.items():
    if d['Q'] > 2 or d['inc'] < 30:
        continue
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', nm + '_rotmod.dat'), comments='#')
    R, Vo, eV, Vg, Vd, Vb = a[:, 0] * kpc, a[:, 1] * 1e3, a[:, 2] * 1e3, a[:, 3] * 1e3, a[:, 4] * 1e3, a[:, 5] * 1e3
    Vb2 = Vg * np.abs(Vg) + 0.5 * Vd**2 + 0.7 * Vb**2
    m = (Vb2 > 0) & (Vo > 0)
    if m.sum() < 5:
        continue
    ALLE.append(dict(name=nm, R=R[m], Vo=Vo[m], eV=np.maximum(eV[m], 1e3), Vg=Vg[m], Vd=Vd[m], Vb=Vb[m], Rd=d['Rd'],
                     eD=d['eD'] / d['D'], inc=d['inc'], eInc=d['eInc']))

def rechne(g, sig=8.0, fD=1.0, fi=1.0, fU=1.0, dV=None):
    R = g['R'] * fD
    Vo = (g['Vo'] if dV is None else g['Vo'] + dV) * fi
    Vb2 = g['Vg'] * np.abs(g['Vg']) * fD + fU * fD * (0.5 * g['Vd']**2 + 0.7 * g['Vb']**2)
    gN = np.clip(Vb2, 1e-30, None) / R
    go = (Vo**2 + (sig * 1e3)**2 * (R / kpc) / (g['Rd'] * fD)) / R
    return R, gN, go

def modell(gN, kappe):
    n = K.nu(gN, a0)
    return (np.minimum(n, 1 + Q) if kappe else n) * gN

def mc(g):
    fD = max(1 + rng.normal() * g['eD'], 0.2)
    inc = min(max(g['inc'] + rng.normal() * g['eInc'], 20.0), 90.0)
    return dict(fD=fD, fi=math.sin(math.radians(g['inc'])) / math.sin(math.radians(inc)), fU=10**(rng.normal() * 0.1),
                dV=rng.normal(size=g['Vo'].size) * g['eV'])

def kennzahlen(gals, sig, kappe, zz=None):
    """je Galaxie: Median-Rest jenseits R_Q, Sprung (jenseits − diesseits bis 0,1 a₀), Steigungsdifferenz jenseits R_Q"""
    rest, sprung, steig = [], [], []
    for i, g in enumerate(gals):
        kw = zz[i] if zz is not None else {}
        R, gN, go = rechne(g, sig, **kw)
        x = gN / a0
        r = np.log10(go / modell(gN, kappe))
        a_ = x < S_K2; b_ = (x >= S_K2) & (x < 0.1)
        if a_.sum() >= 1:
            rest.append(np.median(r[a_]))
        if a_.sum() >= 2 and b_.sum() >= 2:
            sprung.append(np.median(r[a_]) - np.median(r[b_]))
        if a_.sum() >= 3:
            so = np.polyfit(np.log(R[a_]), 0.5 * np.log(go[a_] * R[a_]), 1)[0]
            sm = np.polyfit(np.log(R[a_]), 0.5 * np.log(modell(gN[a_], kappe) * R[a_]), 1)[0]
            steig.append(so - sm)
    return np.array(rest), np.array(sprung), np.array(steig)

jen = [g for g in ALLE if (rechne(g)[1] / a0 < S_K2).any()]
print(f"   Auswahl wie sparc_aussenfeld_rand (Q ≤ 2, Neigung ≥ 30°, Υ_d = 0,5, Υ_b = 0,7): {len(ALLE)} Galaxien; mit Punkten jenseits R_Q: {len(jen)}; "
      f"Punkte jenseits R_Q: {sum(int((rechne(g)[1]/a0 < S_K2).sum()) for g in ALLE)}")
nmax = max(float(K.nu(rechne(g)[1].min())) for g in jen)
print(f"   ν an diesen Punkten: bis {nmax:.2f}; die Kappe senkt dort um log(ν/(1+Q)) bis {math.log10(nmax/(1+Q)):.3f} dex")
NMC = 200
ZZ = [[mc(g) for g in jen] for _ in range(NMC)]
ERG = {}
for sig in (8.0, 0.0):
    for kappe in (False, True):
        rest, sp, st = kennzahlen(jen, sig, kappe)
        B = lambda v: np.std([np.median(v[rng.integers(0, v.size, v.size)]) for _ in range(1000)])
        MCv = np.array([[np.median(v) for v in kennzahlen(jen, sig, kappe, Z)] for Z in ZZ[:(NMC if sig == 8.0 else 60)]])
        e = [math.hypot(B(v), np.std(MCv[:, k])) for k, v in enumerate((rest, sp, st))]
        ERG[(sig, kappe)] = [(np.median(v), ee, v.size) for v, ee in zip((rest, sp, st), e)]
        lab = ("Kappe" if kappe else "ohne Kappe") + f", σ = {sig:.0f} km/s"
        print(f"   {lab:24s}: Rest jenseits R_Q {np.median(rest):+.3f} ± {e[0]:.3f} dex (n = {rest.size});  Sprung an R_Q {np.median(sp):+.3f} ± {e[1]:.3f} (n = {sp.size});  "
              f"Steigung d ln v/d ln r gemessen − Modell jenseits R_Q {np.median(st):+.3f} ± {e[2]:.3f} (n = {st.size})")
# Bruchteil der Galaxien, die die Kappe besser beschreibt (χ² jenseits R_Q mit Geschwindigkeitsfehlern)
besser = 0; n = 0; dchi = []
for g in jen:
    R, gN, go = rechne(g)
    a_ = gN / a0 < S_K2
    if a_.sum() < 2:
        continue
    ev = g['eV'][a_]; V = np.sqrt(go[a_] * R[a_])
    c0 = np.sum(((V - np.sqrt(modell(gN[a_], False) * R[a_])) / ev)**2)
    c1 = np.sum(((V - np.sqrt(modell(gN[a_], True) * R[a_])) / ev)**2)
    besser += c1 < c0; n += 1; dchi.append(c1 - c0)
print(f"   je Galaxie (≥ 2 Punkte jenseits R_Q, Geschwindigkeitsfehler): Kappe besser in {besser} von {n}; Σ Δχ² (Kappe − ohne) = {np.sum(dchi):.0f}, Median {np.median(dchi):.1f}")
# Υ_d-Spielraum: welches Υ_d bräuchte die Kappe, damit der Rest jenseits R_Q null wird?
def rest_fU(fU, kappe):
    return np.median([np.median(np.log10(go[gN/a0 < S_K2] / modell(gN[gN/a0 < S_K2], kappe)))
                      for g in jen for (R, gN, go) in [rechne(g, 8.0, fU=fU)] if (gN/a0 < S_K2).any()])
for kappe in (False, True):
    print(f"   Rest jenseits R_Q bei Υ_d = 0,25 / 0,5 / 1,0 (× 0,5 / 1 / 2): " + " / ".join(f"{rest_fU(f, kappe):+.3f}" for f in (0.5, 1.0, 2.0))
          + ("  (Kappe)" if kappe else "  (ohne Kappe)"))

# 3b: nur Kurven, die R_Q von innen nach außen kreuzen; nur Punkte außerhalb des Kreuzungsradius, nach Tiefe
def kreuz(g, sig=8.0, **kw):
    R, gN, go = rechne(g, sig, **kw)
    x = gN / a0
    innen = np.where(x >= S_K2)[0]
    if innen.size < 2:
        return None
    k0 = innen.max() + 1                                    # alle Punkte danach liegen jenseits R_Q
    if R.size - k0 < 1:
        return None
    return R, gN, go, k0

QUER = [g for g in ALLE if kreuz(g) is not None]
BIN = ((0.015, S_K2), (0.0075, 0.015), (0.0, 0.0075))
def bins(gals, sig, kappe, zz=None):
    out = [[] for _ in BIN]
    for i, g in enumerate(gals):
        kw = zz[i] if zz is not None else {}
        k = kreuz(g, sig, **kw)
        if k is None:
            continue
        R, gN, go, k0 = k
        x = gN[k0:] / a0; r = np.log10(go[k0:] / modell(gN[k0:], kappe))
        for j, (lo, hi) in enumerate(BIN):
            m = (x >= lo) & (x < hi)
            if m.any():
                out[j].append(np.median(r[m]))
    return [np.array(v) for v in out]

print(f"\n   Nur Kurven, die R_Q kreuzen (≥ 2 Punkte mit a_N ≥ 0,029·a₀ innen), nur Punkte außerhalb der letzten Kreuzung: {len(QUER)} Galaxien, "
      f"{sum(kreuz(g)[0].size - kreuz(g)[3] for g in QUER)} Punkte")
print("   Median-Rest je Galaxie, Bereiche a_N/a₀ = 0,015–0,029 | 0,0075–0,015 | < 0,0075; Fehler Bootstrap ⊕ Monte-Carlo (200); in Klammern log(ν/(1+Q)) im Bereich")
ZQ = [[mc(g) for g in QUER] for _ in range(NMC)]
for sig in (8.0, 0.0):
    for kappe in (False, True):
        b0 = bins(QUER, sig, kappe)
        MCb = [bins(QUER, sig, kappe, Z) for Z in ZQ[:(NMC if sig == 8.0 else 60)]]
        z = []
        for j, (lo, hi) in enumerate(BIN):
            v = b0[j]
            if v.size == 0:
                z.append("—"); continue
            eb = np.std([np.median(v[rng.integers(0, v.size, v.size)]) for _ in range(1000)])
            em = np.std([np.median(m_[j]) for m_ in MCb if m_[j].size])
            sh = math.log10(float(K.nu(math.sqrt(lo * max(hi, 1e-4)) * a0 if lo > 0 else 0.005 * a0)) / (1 + Q))
            z.append(f"{np.median(v):+.3f} ± {math.hypot(eb, em):.3f} (n = {v.size}, Δ ≈ {sh:.2f})")
        print(f"   {('Kappe' if kappe else 'ohne Kappe'):10s} σ = {sig:.0f}: " + " | ".join(z))

# ------------------------------------------------------------------ 4 Lesart (iv): Verlustraten addiert
hdr("4  Lesart (iv): Verlustraten auf derselben Rückkehr addiert, s_tot² = a_N/a₀ + θ̃τ_D;  a/a_N = 1 + w·Q + (1 − w)(ν(s_tot) − 1)")
tauD = c / a0
def theta_schale(r):
    """expandierende Schale (Kac-Lagrange-Masse), die heute bei r liegt: θ(t₀) = 3V/R, θ̃ = min über die Geschichte, 1 + Δ(<r)"""
    M = M_L_kac(r)
    def run(di):
        ri = (3 * M / (4 * math.pi * RHO_I * (1 + di))) ** (1 / 3)
        vi = Hb(AI) * ri * (1 - di / 3)
        f = lambda t, y: [y[1], -G * M / y[0]**2 + OL * H0**2 * y[0]]
        return solve_ivp(f, (TI, t0), [ri, vi], rtol=1e-9, atol=1e-6, dense_output=True)
    fn = lambda ld: run(10**ld).y[0][-1] - r
    ld = brentq(fn, -6.0, math.log10(0.5), xtol=1e-9)
    sol = run(10**ld)
    tt = np.geomspace(TI * 10, t0, 400); Y = sol.sol(tt)
    th = 3 * Y[1] / Y[0]
    if (Y[1] < 0).any():
        return 0.0, 0.0, M / (4 / 3 * math.pi * r**3 * rho_m)
    return 3 * sol.y[1][-1] / sol.y[0][-1], max(th.min(), 0.0), M / (4 / 3 * math.pi * r**3 * rho_m)

def lesart_iv(r, thtil):
    xN_ = G * Mb / r**2 / a0
    w = K.gewicht_w(thtil, a0)
    s_tot = math.sqrt(xN_ + thtil * tauD)
    nu_tot = 1 / (-math.expm1(-s_tot))
    return 1 + w * Q + (1 - w) * (nu_tot - 1), w, s_tot, nu_tot

R4 = np.concatenate([np.geomspace(10, 0.9999 * R_ta / kpc, 12), [1.0005 * R_ta / kpc, 1.002 * R_ta / kpc, 1.005 * R_ta / kpc, 1.01 * R_ta / kpc,
                     1.02 * R_ta / kpc, 1.04 * R_ta / kpc, 1.07 * R_ta / kpc, 1.1 * R_ta / kpc, 1.15 * R_ta / kpc, 1.25 * R_ta / kpc,
                     1.4 * R_ta / kpc, 1.7 * R_ta / kpc, 2000, 3000, 5000]])
print("   r[kpc]   1+Δ(<r)  θ/3H₀    θ̃/θ_c     w        s_N       s_tot    ν(s_tot)  a/a_N (iv)   Regel ν(a_N)")
Z4 = []
for r_ in R4:
    r = r_ * kpc
    if r <= R_ta:
        th, tht, Dl = 0.0, 0.0, M_L_kac(r) / (4 / 3 * math.pi * r**3 * rho_m)
    else:
        th, tht, Dl = theta_schale(r)
    m4, w, st, nt = lesart_iv(r, tht)
    Z4.append((r_, Dl, th, tht, w, m4))
    print(f"   {r_:7.1f}  {Dl:8.3f}  {th/(3*H0):7.4f}  {tht*tauD:8.3f}  {w:7.4f}  {math.sqrt(G*Mb/r**2/a0):8.5f}  {st:8.4f}  {nt:8.3f}  {m4:9.3f}   {float(K.nu(G*Mb/r**2)):8.2f}")
Z4 = np.array(Z4)
au = Z4[:, 0] > R_ta / kpc
k = np.argmin(np.where(au, Z4[:, 5], np.inf))
print(f"   außerhalb R_ta: Minimum von a/a_N (iv) {Z4[k,5]:.3f} bei r = {Z4[k,0]:.1f} kpc (w = {Z4[k,4]:.3f}); weit außen 1 + w·Q → {1 + K.gewicht_w(3*H0, a0)*Q:.3f} "
      f"(w im Hintergrund {K.gewicht_w(3*H0, a0):.4f}); außerhalb R_ta überall < 1 + Q")
# w = ½-Schale und der Wert dort analytisch: θ̃τ_D = (ln 2)² ⇒ s_tot ≥ ln 2
sh = math.log(2)
print(f"   bei w = ½: θ̃τ_D = (ln 2)² = {sh**2:.4f}, s_tot ≈ ln 2 (a_N ≪) ⇒ a/a_N = 1 + Q/2 + ½(1/(1 − e^(−ln 2)) − 1) = {1 + Q/2 + 0.5*(1/(1-0.5)-1):.3f}")
wv = np.linspace(0.001, 0.999, 999)
fv = 1 + wv * Q + (1 - wv) * (1 / wv - 1)                       # a_N → 0: ν(s_tot) = 1/w
print(f"   a_N → 0 (ν(s_tot) = 1/w): a/a_N = 1 + wQ + (1 − w)²/w; Minimum {fv.min():.3f} bei w = {wv[np.argmin(fv)]:.3f} (analytisch w = 1/√(1+Q) = {1/math.sqrt(1+Q):.3f}, Wert 2√(1+Q) − 1 = {2*math.sqrt(1+Q)-1:.3f})")
# Stetigkeit an R_ta
kin = np.where(Z4[:, 0] < R_ta / kpc)[0].max()
print(f"   an R_ta: innen {Z4[kin,5]:.2f}, außen 1,0005·R_ta: {Z4[kin+1,5]:.2f}, 1,002·R_ta: {Z4[kin+2,5]:.2f}, 1,005·R_ta: {Z4[kin+3,5]:.2f}; "
      f"stetig, weil θ̃ → 0 an der Umkehr, aber der Abfall 67 → 4,1 liegt in {(Z4[k,0]-R_ta/kpc):.1f} kpc ({(Z4[k,0]*kpc/R_ta-1)*100:.1f} % von R_ta)")
dl4 = np.diff(np.log(Z4[:, 5])) / np.diff(np.log(Z4[:, 0]))
print(f"   größter Abfall d ln M_gel/d ln r = {dl4.min():+.1f} (Gauß: negative gelesene Dichte in der w-Schale)")
# Gegenprobe θ/3H₀ gegen das Energieintegral der Tabelle in trichterkette_gewicht (andere Rechnung): 1+Δ = 5 → 0,5211, 8 → 0,2717, 3 → 0,7245, 2 → 0,8482
tabE = {1.0: 1.0, 1.5: 0.9192, 2.0: 0.8482, 3.0: 0.7245, 5.0: 0.5211, 8.0: 0.2717}
xs = sorted(tabE); ys = [tabE[x] for x in xs]
print("   Gegenprobe θ/3H₀ gegen die Energieintegral-Tabelle (trichterkette_gewicht, 1+Δ → θ/3H₀, linear interpoliert):")
for r_, Dl, th, tht, w, m4 in Z4:
    if r_ > R_ta / kpc and 1.5 <= Dl <= 8.0:
        print(f"      r = {r_:6.0f} kpc: 1+Δ = {Dl:.3f}, Schale θ/3H₀ = {th/(3*H0):.4f}, Tabelle ≈ {np.interp(Dl, xs, ys):.4f}")
for Dz, soll in ((5.0, 0.5211), (8.0, 0.2717)):
    rz = brentq(lambda r: M_L_kac(r) / (4 / 3 * math.pi * r**3 * rho_m) - Dz, 1.0001 * R_ta, 5 * Mpc)
    print(f"      1+Δ = {Dz:.0f} genau (r = {rz/kpc:.0f} kpc): Schale θ/3H₀ = {theta_schale(rz)[0]/(3*H0):.4f}, Tabelle {soll}")
print("   SPARC: Alle Punkte liegen in umgekehrten Bereichen (θ̃ = 0, w = 0): (iv) ist dort identisch mit der Regel ohne Kappe (Rest jenseits R_Q oben).")

# ------------------------------------------------------------------ 5 Prüfung von Hand
hdr("5  Unabhängige Prüfung von Hand (Zahlen eingetippt, ohne ftm_kern)")
GMb = 6e10 * 1.32712440018e20; kp = 3.0856775814913673e19
rst = math.sqrt(GMb / 1.1e-10); Qh = 0.1200 / 0.02237; qh = Qh / (1 + Qh)
RQh = rst / (-math.log(qh))
print(f"   r* = √(GM/a₀) = {rst/kp:.3f} kpc, q = 5,3643/6,3643 = {qh:.5f}, R_Q = r*/{-math.log(qh):.5f} = {RQh/kp:.2f} kpc (Skript {R_Q/kpc:.2f})")
r580 = 0.58e3 * kp; xh = GMb / r580**2 / 1.1e-10; nuh = 1 / (1 - math.exp(-math.sqrt(xh)))
print(f"   an 0,58 Mpc: a_N/a₀ = {xh:.4e}, s = {math.sqrt(xh):.5f}, ν = 1/(1 − e^(−s)) = {nuh:.2f} ≈ 1/s + ½ = {1/math.sqrt(xh)+0.5:.2f};  ν/(1 + Q) = {nuh/(1+Qh):.2f}")
# Umkehrdichte: Einstein-de Sitter 9π²/16 = 5,55 bei t_ta = t₀; mit Λ oben. R_ta aus 1 + Δ_ta von Hand:
rhm = 3 * (67.4e3 / 3.0856775814913673e22)**2 / (8 * math.pi * 6.6743e-11) * 0.3134
Rh = (3 * (1 + Qh) * 6e10 * 1.98841e30 / (4 * math.pi * rhm * D_ta))**(1 / 3)
print(f"   R_ta = (3(1+Q)M_b/(4π ρ̄_m (1+Δ_ta)))^(1/3) mit 1+Δ_ta = {D_ta:.3f}: {Rh/kp:.1f} kpc (Skript {R_ta/kpc:.1f})")
# SPARC: eine Galaxie mit Punkten jenseits R_Q von Hand
gh = min(QUER, key=lambda g: (lambda k: k[1][k[3]:].min())(kreuz(g)))
rows = [[float(v) for v in l.split()] for l in open(os.path.join(sf.D, 'Rotmod_LTG', gh['name'] + '_rotmod.dat')) if l.strip() and not l.startswith('#')]
rows = [r for r in rows if (r[3] * abs(r[3]) + 0.5 * r[4]**2 + 0.7 * r[5]**2) > 0 and r[1] > 0]
kmin = min(range(len(rows)), key=lambda k: (rows[k][3] * abs(rows[k][3]) + 0.5 * rows[k][4]**2 + 0.7 * rows[k][5]**2) / rows[k][0])
print(f"   {gh['name']} (R_d = {gh['Rd']} kpc), tiefster Punkt außerhalb der Kreuzung von R_Q:")
Rk, V, eVh, Vg_, Vd_, Vb_ = rows[kmin][:6]
vb2 = Vg_ * abs(Vg_) + 0.5 * Vd_**2 + 0.7 * Vb_**2
gn = vb2 * 1e6 / (Rk * kp); x = gn / 1.1e-10; nuh = 1 / (1 - math.exp(-math.sqrt(x)))
vc2 = V**2 + 64.0 * Rk / gh['Rd']
print(f"   R = {Rk} kpc, V_obs = {V}, V_bar² = {vb2:.1f} (km/s)², a_N/a₀ = {x:.5f}, ν = {nuh:.3f}; V_c = {math.sqrt(vc2):.1f}, "
      f"V_Regel = {math.sqrt(nuh*vb2):.1f}, V_Kappe = {math.sqrt((1+Qh)*vb2):.1f} km/s;  Rest {math.log10(vc2/(nuh*vb2)):+.3f} / {math.log10(vc2/((1+Qh)*vb2)):+.3f} dex")
_R, _gN, _go = rechne(gh)
_k = int(np.argmin(np.abs(_R / kpc - Rk)))
print(f"   Skript am selben Punkt: a_N/a₀ = {_gN[_k]/a0:.5f}, Rest {np.log10(_go[_k]/modell(_gN[_k], False)):+.3f} / {np.log10(_go[_k]/modell(_gN[_k], True)):+.3f} dex")
