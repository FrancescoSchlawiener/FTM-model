# Wahl des Imaginärteils von τ_M in der Feldgleichung des Autors
# (A) |q_A| = Δs/r, Im τ = ln(r/Δs)/2π, Kern Δ(τ)·χ,  χ = s/(e^s−1)
# (B) |q_B| = e^{−r*/r}, Im τ = r*/(2πr), Bindung ν = 1/(1−|q_B|)
# Wo unterscheiden sich beide beobachtbar?
import os, sys, io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G = 6.674e-11; c = 2.99792458e8; a0 = 1.2e-10; Msun = 1.989e30; kpc = 3.0857e19; AU = 1.496e11
out = []; P = lambda s="": out.append(s)

def chi(s): return s/np.expm1(s)
def nu(s): return 1/(-np.expm1(-s))

def tau_R(N):                       # Koeffizienten von q∏(1−qⁿ)²⁴
    p = [0]*(N+1); p[0] = 1
    for n in range(1, N+1):
        for _ in range(24):
            for k in range(N, n-1, -1): p[k] -= p[k-n]
    return [p[k-1] for k in range(1, N+1)]
TR = tau_R(10)

def logDelta_it(t):                  # ln Δ(i t), t > 0, mit Δ(i/t) = t¹²Δ(i t)
    if t < 1: return -12*np.log(t) + logDelta_it(1/t)
    q = np.exp(-2*np.pi*t); n = np.arange(1, 200)
    return np.log(q) + 24*np.sum(np.log1p(-q**n))
def Delta_c(tau, N=8):               # Δ(τ) komplex, für kleines |q|
    q = np.exp(2j*np.pi*tau); n = np.arange(1, N+1)
    return q*np.prod((1-q**n)**24)
def R_A(qa):                         # Δ/q auf der imaginären Achse, q = q_A (klein)
    n = np.arange(1, 8); return np.prod((1-qa**n)**24)

def vf(M): return (G*M*Msun*a0)**0.25
def eps(M): return vf(M)**2/(3*c**2)          # Δs/r* = v_f²/(3c²)

P("0) Fourier-Koeffizienten τ_R(n) von Δ: " + ", ".join(str(x) for x in TR))
P("   Bindungssumme von (B): Σ|q|ⁿ, alle Koeffizienten 1")
P()

# 1) Identität in erster Ordnung
P("1) χ(s)·(Δs/r)·(3c²/r*) gegen a_N·|q_B|/(1−|q_B|)   (Punktmasse, s = r*/r)")
worst = 0
for M in [1e8, 6e10, 1e14]:
    GM = G*M*Msun; rs = np.sqrt(GM/a0); ds = GM/(3*c**2)
    for x in [0.01, 0.1, 0.3, 1, 3, 10, 100]:
        r = x*rs; aN = GM/r**2; s = np.sqrt(aN/a0); qB = np.exp(-s)
        A = chi(s)*(ds/r)*(3*c**2/rs); B = aN*qB/(1-qB)
        worst = max(worst, abs(A/B-1))
P(f"   größte relative Abweichung über M = 1e8…1e14 M☉, r/r* = 0,01…100: {worst:.1e}")
P("   also χ(s) = s·|q_B|/(1−|q_B|): der Schalter von (A) enthält den Nome von (B)")
P()

# 2) Höhere Fourier-Glieder
P("2a) (A): Δ(τ_A)/q_A − 1 und a_A/a_B − 1 am Übergang, q_A = (v_f²/3c²)·s")
P("      M [M☉]    v_f [km/s]   q_A(r*)     Im τ_A(r*)   Δ/q_A−1 bei r*   max|a_A/a_B−1|   bei r/r*")
for M in [1e8, 1e10, 6e10, 1e12, 1e14, 1e15]:
    e = eps(M); qs = e; imA = np.log(1/e)/(2*np.pi)
    xs = np.logspace(-2, 2, 4001); sv = 1/xs
    dev = np.array([(R_A(e*s)-1)*(nu(s)-1)/nu(s) for s in sv]); k = np.argmax(abs(dev))
    P(f"   {M:9.0e}   {vf(M)/1e3:9.1f}   {qs:.3e}   {imA:9.3f}    {R_A(qs)-1:+.3e}      {dev[k]:+.3e}      {xs[k]:.2f}")
P("   Näherung: max|a_A/a_B − 1| = 24·(v_f²/3c²)/e = 2,94·(v_f/c)²  bei r = r*")
P()
P("2b) (B): Bindung mit Koeffizienten 1 gegen Δ(τ_B) als Kern,  a/a_N = 1 + Σ|q|ⁿ  bzw. 1 + Δ(i s/2π)")
P("     r/r*     s       |q_B|      Im τ_B    Σ|q|ⁿ = ν−1     Δ(τ_B)         Verhältnis")
for x in [0.1, 0.2, 0.304, 0.5, 1, 3, 10]:
    s = 1/x; q = np.exp(-s); D = np.exp(logDelta_it(s/(2*np.pi)))
    P(f"   {x:6.3f}   {s:6.3f}   {q:.3e}   {s/(2*np.pi):7.4f}   {q/(1-q):.4e}     {D:.4e}     {D/(q/(1-q)):.3e}")
P("   Δ(τ_B) hat sein Maximum bei r = 0,304 r* und verschwindet für r → ∞ (Δ(i/t) = t¹²Δ(it)).")
P("   In (B) trägt die Bindung daher die Reihe mit Koeffizienten 1; Δ ist dort nur das erste Glied.")
P("   Unterschied bei r ≪ r*: Σ|q|ⁿ − Δ = 25|q|² + …")
P()

# 3) Phase des Nomes
P("3) Phase: Re τ = √(1+χ) (Punktmasse), Milchstraße 6e10 M☉")
e = eps(6e10)
P("     r/r*    Re τ    (A) |Δ(τ)|/Δ(i Im τ) − 1  |  (B) ν = 1/(1−|q|)   Re 1/(1−q)   |1/(1−q)|")
for x in [0.3, 1, 3, 10]:
    s = 1/x; ret = np.sqrt(1+chi(s)); qa = e*s
    tau = ret + 1j*np.log(1/qa)/(2*np.pi)
    D = Delta_c(tau); D0 = qa*R_A(qa)
    d1 = abs(D)/D0-1; d2 = (D*np.exp(-2j*np.pi*ret)).real/D0-1
    q = np.exp(-s)*np.exp(2j*np.pi*ret)
    P(f"   {x:6.1f}   {ret:.4f}        {d1:+.2e}             |     {nu(s):7.4f}          {(1/(1-q)).real:7.4f}      {abs(1/(1-q)):7.4f}")
P("   Die Phase in den höheren Gliedern ändert (A) um ≤ 24|q_A|, (B) um Größen der Ordnung 1.")
P()

# 4) Resonanzradien
P("4) Resonanzradien Re τ = m/n (Punktmasse 6e10 M☉): gleich in (A) und (B) bis auf die Δ-Glieder")
M = 6e10; GM = G*M*Msun; rs = np.sqrt(GM/a0); e = eps(M)
def lna(lr, model):
    r = np.exp(lr); s = 1/r; base = -2*lr + np.log(nu(s))
    if model == 'B': return base
    return -2*lr + np.log(1 + (nu(s)-1)*R_A(e*s))
def ret2(lr, model, h=1e-4):
    d = (lna(lr+h, model) - lna(lr-h, model))/(2*h); return 3 + d
P("    m/n    r [kpc]   s      |q_A|       |q_B|    Glied n: |τ_R(n)|·|q_A|ⁿ   |q_B|ⁿ     Verschiebung δr/r (A−B)")
for m, n in [(6, 5), (5, 4), (4, 3), (7, 5)]:
    t2 = (m/n)**2
    lB = brentq(lambda l: ret2(l, 'B')-t2, np.log(0.05), np.log(500), xtol=1e-15)
    lA = brentq(lambda l: ret2(l, 'A')-t2, np.log(0.05), np.log(500), xtol=1e-15)
    r = np.exp(lB); s = 1/r; qa = e*s; qb = np.exp(-s)
    P(f"   {m}/{n}   {r*rs/kpc:7.2f}   {s:5.3f}   {qa:.3e}   {qb:.3f}      {abs(TR[n-1])*qa**n:.1e}          {qb**n:.2e}      {np.exp(lA-lB)-1:+.1e}")
P()

# 5) Sonnensystem
P("5) Sonnensystem (M = 1 M☉, r* = %.0f AE)" % (np.sqrt(G*Msun/a0)/AU))
GM = G*Msun; rs = np.sqrt(GM/a0); ds = GM/(3*c**2)
P("   Ort          q_A = Δs/r   Im τ_A    s = r*/r    Im τ_B     lg|q_B|     lg a_K [m/s²] (A = B)   q_A·c²r/GM")
for name, r in [("Merkur", 0.387*AU), ("Erde", AU), ("Neptun", 30.07*AU), ("r* der Sonne", rs)]:
    s = rs/r; aN = GM/r**2
    lgaK = np.log10(aN) - s/np.log(10) - np.log10(-np.expm1(-s))
    P(f"   {name:12s} {ds/r:.3e}    {np.log(r/ds)/(2*np.pi):6.3f}   {s:9.2f}   {s/(2*np.pi):9.3f}   {-s/np.log(10):10.1f}      {lgaK:10.1f}          {ds/r*c**2*r/GM:.4f}")
P("   (A): q_A = GM/(3c²r) ist ein Glied erster post-Newtonscher Ordnung; unterdrückt wird der Kern nur durch χ.")
P("   (B): der Nome selbst ist e^{−s}; jede Funktion von τ_B ist an der Spitze.")
P()

# 6) Gedächtnis: Schalter vor dem Integral (A) oder im Nome (B)
P("6) Gedächtnisintegral, Einheiten GM = r* = a0 = 1 (v_f = 1, Zeit r*/v_f), Kreisfrequenz Ω")
P("   a = −e_r/r² − C·κ(r)·(e_r·m)·e_r,   ṁ = K·v − λm,   κ(r) = (λ² + ω_c²)/(r ω_c²)  (Kreisbahn am Ort geeicht)")
P("   (A): C = χ(s(t)), K = s   (Schalter zur Zeit t, Nome längs der Bahn)")
P("   (B): C = 1,      K = s·χ(s) = s²|q_B|/(1−|q_B|)   (ganze Bindung längs der Bahn)")
P("   Auf jeder Kreisbahn geben beide exakt a = a_N/(1−e^{−s}).")

def rhs_rot(y, L, lam, model, kfix=None):   # mitrotierend: m = m_r e_r + m_φ e_φ, L erhalten (Kraft radial)
    r, p, mr, mp = y; s = 1/r; W = L/r**2
    w2 = nu(s)*s**3; kap = (lam**2 + w2)/(r*w2) if kfix is None else kfix
    C, K = (chi(s), s) if model == 'A' else (1.0, s*chi(s))
    return np.array([p, L**2/r**3 - s*s - C*kap*mr, K*p - lam*mr + W*mp, K*L/r - lam*mp - W*mr])
def fixpunkt(r0, lf, model):
    s = 1/r0; W = np.sqrt(nu(s)*s**3); L = W*r0**2; lam = lf*W
    K = s if model == 'A' else s*chi(s)
    mr = K*r0*W**2/(lam**2 + W**2)
    return W, L, lam, np.array([r0, 0.0, mr, lam*mr/W])
def moden(r0, lf, model, fest=False):
    W, L, lam, y0 = fixpunkt(r0, lf, model); J = np.zeros((4, 4))
    kf = (lam**2 + W**2)/(r0*W**2) if fest else None
    for j in range(4):
        h = 1e-6*max(1e-3, abs(y0[j])); e_ = np.zeros(4); e_[j] = h
        J[:, j] = (rhs_rot(y0+e_, L, lam, model, kf) - rhs_rot(y0-e_, L, lam, model, kf))/(2*h)
    ev = np.linalg.eigvals(J)/W
    return ev[np.argsort(-ev.imag)][:2]

P()
P("6a) Kreisbahn, lineare Moden (Frequenz/Ω, Anwachsrate/Ω); ohne Gedächtnis: eine Mode √(1+χ)")
P("     r/r*   λ/Ω    √(1+χ)  |  (A) Frequenzen     größte Rate  |  (B) Frequenzen     größte Rate")
for r0 in [0.1, 0.35, 1.0, 3.0]:
    for lf in [0.01, 1.0]:
        eA = moden(r0, lf, 'A'); eB = moden(r0, lf, 'B')
        P(f"   {r0:5.2f}   {lf:4.2f}   {np.sqrt(1+chi(1/r0)):.4f}  |  {eA[0].imag:.4f}  {eA[1].imag:.4f}     {eA.real.max():+.4f}   |  {eB[0].imag:.4f}  {eB[1].imag:.4f}     {eB.real.max():+.4f}")
P()
P("   Anwachsen in (B) über r/r* = 0,01…2 (A: alle Moden gedämpft für λ/Ω = 0,001…100, r/r* = 0,01…100):")
rr = np.logspace(-2, 0.3, 400)
gAmax = max(moden(r, lf, 'A').real.max() for r in np.logspace(-2, 2, 200) for lf in [0.001, 0.01, 0.1, 1, 10, 100])
for lf in [0.001, 0.01, 0.1, 0.3, 1, 2]:
    g = np.array([moden(r, lf, 'B').real.max() for r in rr]); i = np.where(g > 1e-3)[0]
    band = f"{rr[i[0]]:.2f}–{rr[i[-1]]:.2f} r*" if len(i) else "keines"
    P(f"   λ/Ω = {lf:5.3f}:  größte Rate {g.max():+.4f} Ω bei r = {rr[g.argmax()]:.2f} r*,  Rate > 10⁻³Ω auf {band}")
P(f"   (A): größte Rate über alle Fälle {gAmax:+.0e} Ω (Rundungsgrenze)")
g = np.array([moden(r, 0.01, 'B', fest=True).real.max() for r in rr]); i = np.where(g > 1e-3)[0]
P(f"   κ fest an der Ausgangsbahn statt am Ort, λ/Ω = 0,01: (B) größte Rate {g.max():+.4f} Ω bei {rr[g.argmax()]:.2f} r*, Band {rr[i[0]]:.2f}–{rr[i[-1]]:.2f} r*;"
  f" (A) {max(moden(r, 0.01, 'A', fest=True).real.max() for r in rr):+.0e} Ω")
Myr = 3.15576e13; e6 = vf(6e10); rs6 = np.sqrt(G*6e10*Msun/a0)
W35 = np.sqrt(nu(1/0.35)/0.35**3)*e6/rs6
P(f"   Milchstraße 6e10 M☉ (r* = {rs6/kpc:.2f} kpc): Band 0,07–0,55 r* = {0.07*rs6/kpc:.1f}–{0.55*rs6/kpc:.1f} kpc, 1/Ω(0,35 r*) = {1/W35/Myr:.1f} Myr,"
  f" λ/Ω für 1/λ_m = 0,1 Gyr: {1/(100*Myr)/W35:.3f}")

opt = dict(method='DOP853', rtol=1e-11, atol=1e-14)
P()
P("   Nichtlineare Rechnung, r0 = 0,35 r*, λ/Ω = 0,01, Anfangsauslenkung 10⁻⁴: max|r/r0 − 1| nach 0…5, 5…10, 10…15 Umläufen")
for model in 'AB':
    W, L, lam, y0 = fixpunkt(0.35, 0.01, model); y0 = y0.copy(); y0[0] *= 1 + 1e-4
    T = 2*np.pi/W; sol = solve_ivp(lambda t, y: rhs_rot(y, L, lam, model), (0, 15*T), y0, dense_output=True, **opt)
    amp = [np.abs(sol.sol(np.linspace(k*5*T, (k+1)*5*T, 4000))[0]/0.35 - 1).max() for k in range(3)]
    P(f"   ({model}):  " + "   ".join(f"{x:.1e}" for x in amp))

# exzentrische Bahnen: Gedächtnis eingeschwungen auf der Bahn im lokalen Feld, dann volle Dynamik
def Phi(r): return -np.log(np.expm1(1/r))
def Lz(rp, ra): return np.sqrt(2*(Phi(ra)-Phi(rp))/(1/rp**2-1/ra**2))
def rhs_ref(t, y, L, lam):
    r, pr, ph, ax, ay, bx, by = y; s = 1/r
    co, si = np.cos(ph), np.sin(ph); phd = L/r**2
    vx = pr*co - r*phd*si; vy = pr*si + r*phd*co; KA = s; KB = s*chi(s)
    return [pr, L**2/r**3 - nu(s)*s*s, phd, KA*vx-lam*ax, KA*vy-lam*ay, KB*vx-lam*bx, KB*vy-lam*by]
def rhs_mem(t, y, L, lam, model):
    r, pr, ph, mx, my = y; s = 1/r
    co, si = np.cos(ph), np.sin(ph); phd = L/r**2
    vx = pr*co - r*phd*si; vy = pr*si + r*phd*co
    w2 = nu(s)*s**3; kap = (lam**2 + w2)/(r*w2)
    C, K = (chi(s), s) if model == 'A' else (1.0, s*chi(s))
    return [pr, L**2/r**3 - s*s - C*kap*(co*mx + si*my), phd, K*vx-lam*mx, K*vy-lam*my]
def apo(t, y, *a): return y[1]
apo.direction = -1
def peri(t, y, *a): return y[1]
peri.direction = 1
opt = dict(method='DOP853', rtol=1e-10, atol=1e-12)
def exzentrisch(rp, ra, lf, N=20):
    L = Lz(rp, ra); rg = brentq(lambda r: L**2 - r*nu(1/r), 1e-3, 1e4)
    Og = L/rg**2; lam = lf*Og; Tr = 2*np.pi/np.sqrt(1+chi(1/rg))/Og
    y = solve_ivp(rhs_ref, (0, max(12/lam, 40*Tr)), [ra, 0, 0, 0, 0, 0, 0], args=(L, lam), **opt).y[:, -1]
    s2 = solve_ivp(rhs_ref, (0, 3*Tr), y, args=(L, lam), events=apo, **opt)
    y = s2.y_events[0][np.argmax(s2.t_events[0] > 1e-6)]
    res = {}
    for model in 'AB':
        m = y[3:5] if model == 'A' else y[5:7]
        s_ = solve_ivp(rhs_mem, (0, 4*N*Tr), list(y[:3]) + list(m), args=(L, lam, model), events=[apo, peri], **opt)
        ras = np.concatenate([[y[0]], s_.y_events[0][s_.t_events[0] > 1e-6][:, 0]]); rps = s_.y_events[1][:, 0]
        n = min(len(ras), len(rps), N+1); amp = (ras[:n] - rps[:n])/(ras[0] - rps[0])
        res[model] = (amp, ras[n-1]/ras[0], n-1)
    return res
P()
P("6b) Exzentrische Bahnen, volle Dynamik (Gedächtnis zuvor auf der Bahn im lokalen Feld eingeschwungen):")
P("    Hub (r_a − r_p) relativ zum Anfang nach 5 und 20 Radialperioden, r_a am Ende relativ zum Anfang")
P("    r_p–r_a [r*]   λ/Ω   |  (A) nach 5    nach 20   r_a   |  (B) nach 5    nach 20   r_a")
for rp, ra in [(0.2, 1.0), (0.3, 1.5), (1.0, 5.0), (3.0, 15.0)]:
    for lf in [0.01, 0.1]:
        R = exzentrisch(rp, ra, lf)
        def f(model):
            amp, rr_, n = R[model]
            return f"{amp[min(5, n)]:6.3f}   {amp[n]:6.3f}{'' if n == 20 else '(' + str(n) + ')'}  {rr_:5.2f}"
        P(f"   {rp:4.1f}–{ra:4.1f}      {lf:4.2f}  |     {f('A')}   |     {f('B')}")
P()

# 7) Sprunghafte Änderung des Baryonenfeldes
P("7) Plötzliche Änderung a_N → f·a_N an einem Ort (Stoß, Gasverlust): Zusatzbeschleunigung unmittelbar danach")
P("   (A): Schalter springt sofort, χ(s₂)/χ(s₁);  (B): 1 (alles im Gedächtnis, Nachlauf 1/λ_m)")
P("    s₁      f = 0,25    0,5       2         4        | Endwert s₂χ(s₂)/(s₁χ(s₁)) bei f = 0,25 … 4")
for s1 in [0.3, 1.0, 3.0]:
    row = [chi(s1*np.sqrt(f))/chi(s1) for f in [0.25, 0.5, 2, 4]]
    fin = [s1*np.sqrt(f)*chi(s1*np.sqrt(f))/(s1*chi(s1)) for f in [0.25, 0.5, 2, 4]]
    P(f"   {s1:4.1f}    " + "  ".join(f"{x:7.4f}" for x in row) + "   |  " + "  ".join(f"{x:6.3f}" for x in fin))

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "wahl_im_tau.txt"), "w", encoding="utf-8").write(txt)
print(txt)
