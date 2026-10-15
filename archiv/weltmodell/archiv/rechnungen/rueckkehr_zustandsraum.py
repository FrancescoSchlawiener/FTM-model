# Rückkehr im Zustandsraum: trägt q = e^{−√(a_N/a0)} als Rückkehr quer zum Torus oder als Gedächtniskern?
# 1) Fokker-Planck (Rückwärtsgleichung) mit Drift, Kontraktion λ_⊥, Wirkungsdiffusion; Modellgrößen; Echtzeit
# 2) Gedächtniskern: Lévy-Dichte, Zwanzig-Mori-Kern einer eliminierten Kette, Maxwell-Element, druckfreier Einfall
# 3) Kac-Budget: geschlossener Raum mit offener Unterrichtung; Schließung 𝒫 als Bildsumme
import numpy as np, sys, io, os
import mpmath as mp
from scipy import integrate, special, linalg, optimize
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
G = 6.674e-11; c = 2.99792458e8; a0 = 1.2e-10; Msun = 1.989e30; kpc = 3.0857e19; Gyr = 3.156e16
out = []; P = lambda s="": out.append(s)
rng = np.random.default_rng(20260929)
nu = lambda s: 1/(-np.expm1(-s))

def q_fd(s, b, X, h):
    """Rückwärtsgleichung q'' + b(x)q' − s²q = 0 (L = D = 1), q(0) = 1, reflektierend bei X; Wert bei x = L = 1"""
    n = int(round(X/h)); x = h*np.arange(1, n+1); bx = b(x)
    lo = 1/h**2 - bx/(2*h); up = 1/h**2 + bx/(2*h); di = -2/h**2 - s*s*np.ones(n)
    lo[-1] += up[-1]; up[-1] = 0.0
    rhs = np.zeros(n); rhs[0] = -lo[0]
    ab = np.zeros((3, n)); ab[0, 1:] = up[:-1]; ab[1] = di; ab[2, :-1] = lo[1:]
    return linalg.solve_banded((1, 1), ab, rhs)[int(round(1/h))-1]

P("1) Rückwärtsgleichung für q(ρ) = E[e^{−ΓT}], T = erste Rückkehr zum Torus (ρ = 0), Start im Abstand L")
P("   D q'' + b(ρ) q' − Γ q = 0,  q(0) = 1,  s² = ΓL²/D,  ν = 1/(1−q);  Tully-Fisher verlangt s·ν → 1 für s → 0")
P()
P("1a) konstante Drift zum Torus hin, Pe = vL/2D (Pe < 0: vom Torus weg):  −ln q = √(Pe² + s²) − Pe")
P("      Pe       s·ν bei s = 1      0.1         0.01        0.001     |  Tiefgrenze")
for Pe in [-0.1, 0.0, 1e-3, 1e-2, 0.1, 1.0]:
    row = [s/(1-np.exp(Pe-np.hypot(Pe, s))) for s in [1, 0.1, 0.01, 0.001]]
    lim = (f"a → 2Pe·a0 = {2*Pe:.3g}·a0 (Boden)" if Pe > 0 else
           f"a → {1/(1-np.exp(2*Pe)):.2f}·a_N (Newton mal Zahl)" if Pe < 0 else "a → √(a_N a0)")
    P(f"   {Pe:7.3f}   " + "  ".join(f"{x:10.4f}" for x in row) + f"   |  {lim}")
P("   Kontrolle Differenzenverfahren (Gitter 0.002 L):")
for Pe, s in [(0.5, 1.0), (0.0, 0.1), (-0.3, 0.3)]:
    qa = np.exp(Pe-np.hypot(Pe, s)); qn = q_fd(s, lambda x: -2*Pe+0*x, 40+20/s, 0.002)
    P(f"     Pe = {Pe:5.2f}, s = {s}:  geschlossen {qa:.6f}   numerisch {qn:.6f}")
P()

P("1b) Kontraktion quer zum Torus (Normalform ρ̇ = −κρ + Rauschen, κ = |λ_⊥|), Bandbreite σ² = D/κ, y = L/σ:")
P("    q = e^{y²/4}·D_{−Γ/κ}(y)/D_{−Γ/κ}(0)  (parabolische Zylinderfunktion),  Γ/κ = s²/y²")
def q_ou(s, y):
    n = mp.mpf(s)**2/mp.mpf(y)**2
    return float(mp.e**(mp.mpf(y)**2/4)*mp.pcfd(-n, y)/mp.pcfd(-n, 0))
def T_ou(y):   # mittlere erste Rückkehrzeit in Einheiten L²/D
    I = integrate.quad(lambda z: np.exp(z*z/2)*np.sqrt(np.pi/2)*special.erfc(z/np.sqrt(2)), 0, y)[0]
    return I/y**2
P("      y = L/σ   s·ν bei s = 1      0.1         0.01        0.001     |  Boden a/a0 = τ_D/⟨T⟩   √(2/π)·y")
for y in [0.01, 0.1, 1.0, 3.0]:
    row = [s/(1-q_ou(s, y)) for s in [1, 0.1, 0.01, 0.001]]
    P(f"   {y:8.2f}  " + "  ".join(f"{x:10.4f}" for x in row) + f"   |  {1/T_ou(y):.4e}          {np.sqrt(2/np.pi)*y:.4e}")
P("   Kontrolle Differenzenverfahren:")
for s, y in [(1.0, 0.1), (0.1, 1.0), (1.0, 1.0)]:
    qn = q_fd(s, lambda x: -y*y*x, 40+20/s, 0.002)
    P(f"     s = {s}, y = {y}:  Zylinderfunktion {q_ou(s, y):.6f}   numerisch {qn:.6f}")
P("   Die Wurzel gilt nur für s ≫ y = L/σ; darunter 1 − q ≈ Γ⟨T⟩ (linear), a läuft gegen den Boden a0·τ_D/⟨T⟩.")
P()

P("1c) Wirkungsdiffusion zerfallener Tori (Standardabbildung, K ≫ 1): D_p gegen K²/4")
P("      K     D_p/(K²/4)   1 − 2J₂(K) + 2J₂(K)²")
for K in [10.0, 20.0, 40.0, 80.0]:
    N = 40000; n = 400
    x = rng.uniform(0, 2*np.pi, N); p = rng.uniform(0, 2*np.pi, N); p0 = p.copy()
    for _ in range(n):
        p = p + K*np.sin(x); x = (x + p) % (2*np.pi)
    D = np.var(p - p0)/(2*n); J2 = special.jv(2, K)
    P(f"   {K:6.1f}    {D/(K*K/4):.4f}      {1-2*J2+2*J2*J2:.4f}")
P("   D ∝ K² ∝ (a0/a_N)² hängt selbst von a_N ab. Mit D ∝ a_N^α, Γ ∝ a_N^β gilt s ∝ a_N^{(β−α)/2};")
P("   die Wurzel verlangt β − α = 1. Wirkungsdiffusion mit Γ ∝ a_N gibt (β−α)/2 = 3/2.")
P()

P("1d) Modellgrößen (Punktmasse 6e10 M☉, r = 8 kpc; a0 = c²/(2πR_H), R_H = c/H0)")
M = 6e10*Msun; GM = G*M; ds = GM/(3*c**2); rs = np.sqrt(GM/a0); RH = c**2/(2*np.pi*a0); vf = (GM*a0)**0.25
r = 8*kpc; aN = GM/r**2; s = np.sqrt(aN/a0); Gam = aN/c; tauD = c/a0
P(f"   einzige lokale, universelle Rate ∝ a_N: Γ = a_N/c = 3(Δs/r)(c/r) = {Gam*Gyr:.4f} /Gyr  (3Δs c/r²: {3*ds*c/r**2*Gyr:.4f} /Gyr)")
P(f"   dann τ_D = L²/D = c/a0 = 2πR_H/c = 2π/H0 = {tauD/Gyr:.2f} Gyr; L und D einzeln bleiben frei")
P(f"   s² = Γ·τ_D = {Gam*tauD:.6f};  a_N/a0 = {aN/a0:.6f};  6π(Δs/r)(R_H/r) = {6*np.pi*ds/r*RH/r:.6f}")
P(f"   mittlere Dauer überlebender Rückkehren τ_D/(2s) = {tauD/(2*s)/Gyr:.1f} Gyr bei s = {s:.3f}; bei s = 1: {tauD/2/Gyr:.1f} Gyr")
P()

P("1e) Rückkehr in Echtzeit, seit der Zeit t gelaufen: q_t = E[e^{−ΓT}; T < t], τ_D = 2π/H0")
def q_t(s, tH):
    A = np.sqrt(np.pi/(2*tH)); B = s*np.sqrt(tH/(2*np.pi))
    return 0.5*(np.exp(-s)*special.erfc(A - B) + np.exp(s)*special.erfc(A + B))
fL = lambda u, tau=1.0: np.sqrt(tau/(4*np.pi))*u**-1.5*np.exp(-tau/(4*u))
chk = integrate.quad(lambda u: fL(u)*np.exp(-0.01*u), 0, 1/(2*np.pi))[0]
P(f"   Kontrolle (s = 0.1, tH0 = 1): Formel {q_t(0.1, 1):.6f}, Integral {chk:.6f}")
P("     t·H0     s·ν_t bei s = 1     0.1       0.01   |  q_t(s→0) = erfc(√(π/2tH0))   (verlangt: s·ν = 1.582, 1.051, 1.005)")
for tH in [0.5, 1.0, 2.0, 10.0]:
    row = [s/(1-q_t(s, tH)) for s in [1, 0.1, 0.01]]
    P(f"   {tH:6.1f}    " + "   ".join(f"{x:8.4f}" for x in row) + f"   |  {special.erfc(np.sqrt(np.pi/(2*tH))):.4f}")
need = optimize.brentq(lambda tH: 0.1/(1-q_t(0.1, tH)) - 0.9*0.1*nu(0.1), 1, 1e7)
P(f"   s·ν bis auf 10 % bei s = 0.1 verlangt t·H0 = {need:.0f}")
P()

P("1f) Bahnkohärenz: Kern auf der Geschwindigkeit, K̂(p) = e^{−√((Γ+p)τ_D)}; Kreisbahn sieht K̂(−iω), Ruhe K̂(0) = e^{−s}")
P(f"    Punktmasse 6e10 M☉, r* = {rs/kpc:.2f} kpc, v_f = {vf/1e3:.1f} km/s, c/v_f = {c/vf:.0f}")
P("     r/r*     s       v [km/s]   ω·τ_D = (c/v)s²ν     lg(|K̂(−iω)|/K̂(0))")
for x in [0.3, 1, 3, 10, 100]:
    rr = x*rs; s = 1/x; a = GM/rr**2*nu(s); v = np.sqrt(a*rr); wt = v/rr*tauD
    lr = (s - np.sqrt(complex(s*s, -wt)).real)/np.log(10)
    P(f"   {x:6.1f}   {s:6.3f}   {v/1e3:7.1f}     {wt:12.1f}          {lr:9.2f}")
P("   Der Kern hat bei Δt → 0 kein Gewicht (f(0) = 0) und liegt bei Δt ~ τ_D ≫ 1/ω: die Bahnphase mittelt ihn weg.")
P()

P("2a) Lévy-Dichte f(u) = √(τ/4π)u^{−3/2}e^{−τ/4u}: ∫f(u)e^{−Γu}du gegen e^{−√(Γτ)} (τ = 1)")
for s in [0.1, 1.0, 3.0]:
    I = integrate.quad(lambda u: fL(u)*np.exp(-s*s*u), 0, 1)[0] + integrate.quad(lambda u: fL(u)*np.exp(-s*s*u), 1, np.inf)[0]
    P(f"   s = {s}:  Integral {I:.10f}   e^{{−s}} {np.exp(-s):.10f}")
P()

P("2b) Zwanzig-Mori: Kette y_1…y_N (Rate w je Richtung, keine Rückstellung), Schleife = Platz 0, Einspeisung bei n0 = 50")
P("    Kern K(u) = w·[e^{Au}]_{1,n0} = erste Rückkehrdichte; K̂(Γ) = w·[(Γ − A)⁻¹]_{1,n0}")
def q_chain(Gam, n0, N, w=1.0):
    di = -(2*w + Gam)*np.ones(N); di[-1] = -(w + Gam)
    ab = np.zeros((3, N)); ab[0, 1:] = w; ab[1] = di; ab[2, :-1] = w
    rhs = np.zeros(N); rhs[0] = -w
    return linalg.solve_banded((1, 1), ab, rhs)[n0-1]
n0 = 50; Nk = 3000; tau = n0**2
d = -2*np.ones(Nk); d[-1] = -1; lam, V = linalg.eigh_tridiagonal(d, np.ones(Nk-1))
fch = lambda u: (V[0, :]*np.exp(lam*u)) @ V[n0-1, :]
fc = lambda u: n0/np.sqrt(4*np.pi*u**3)*np.exp(-n0**2/(4*u))
P("     u/τ     Kette K(u)·τ    Lévy f(u)·τ    Verhältnis")
for uu in [0.05, 0.25, 1, 4, 16]:
    u = uu*tau; P(f"   {uu:6.2f}    {fch(u)*tau:.6e}   {fc(u)*tau:.6e}   {fch(u)/fc(u):.5f}")
Gt = 1/tau
P(f"   K̂ aus Eigenmoden {np.sum(V[0,:]*V[n0-1,:]/(Gt-lam)):.8f}, aus Rückwärtsgleichung {q_chain(Gt, n0, Nk):.8f}, e^{{−1}} {np.exp(-1):.8f}")
P("    s·ν bei             s = 1       0.1        0.01")
rows = {"Maxwell-Element (eine Stufe, t_K = τ)": lambda s: 1/(1+s*s),
        f"Kette in Reihe, gerichtet ({n0} Stufen)": lambda s: (1+s*s/n0)**(-n0),
        "Diffusionskette, offen (N = 2·10⁶)": lambda s: q_chain(s*s/tau, n0, 2_000_000)}
for k, fq in rows.items():
    P(f"   {k:40s} " + "  ".join(f"{s/(1-fq(s)):9.4f}" for s in [1, 0.1, 0.01]))
P(f"   {'Regel e^{−s}':40s} " + "  ".join(f"{s*nu(s):9.4f}" for s in [1, 0.1, 0.01]))
P("   Endlich viele Stufen: K̂ rational, 1 − q ∝ Γ. Die Wurzel verlangt einen Verzweigungspunkt von K̂ bei Γ = 0,")
P("   also ein lückenloses Kontinuum eliminierter Moden: eine offene Kette ohne Drift und ohne Rückstellung.")
P()

P("2c) Druckfreier Einfall (Einstein-de Sitter): δ̈ + (4/3t)δ̇ − (2/3t²)δ = 0, G(t,t') = (3/5)t'[(t/t')^{2/3} − t'/t]")
sol = integrate.solve_ivp(lambda t, y: [y[1], -4/(3*t)*y[1] + 2/(3*t*t)*y[0]], (1, 1000), [0, 1], rtol=1e-11, atol=1e-13, dense_output=True)
for t in [10, 100, 1000]:
    Gf = 0.6*(t**(2/3) - 1/t); P(f"   t/t' = {t:5d}:  numerisch {sol.sol(t)[0]:.6f}   Formel {Gf:.6f}   G/t^(2/3) = {Gf/t**(2/3):.4f}")
P("   Der Kern wächst, ∫G dt divergiert: er liefert keine Schleifenverstärkung < 1.")
P()

P("3a) Geschlossener Raum mit offener Unterrichtung: Schleife bei 0, Start L, reflektierend bei Λ")
P("    q = cosh((Λ−L)k)/cosh(Λk), k = s/L;  ⟨T⟩ = L(2Λ−L)/(2D) endlich (Kac);  Λ = L: 1/cosh s")
def q_box(s, lr):
    return np.exp(-s)*(1 + np.exp(-2*(lr-1)*s))/(1 + np.exp(-2*lr*s))
for lr in [10, 100]:
    N = n0*lr; s = 0.1; qc = q_chain(s*s/tau, n0, N)
    di = -2*np.ones(N); di[-1] = -1; ab = np.zeros((3, N)); ab[0, 1:] = 1; ab[1] = di; ab[2, :-1] = 1
    Tn = linalg.solve_banded((1, 1), ab, -np.ones(N))[n0-1]; Lam = N + 0.5   # mittlere erste Rückkehrzeit
    P(f"   Kette Λ/L = {lr}: q(s=0.1) = {qc:.6f}, Formel {q_box(s, Lam/n0):.6f};  ⟨T⟩ = {Tn:.1f}, Formel {n0*(2*Lam-n0)/2:.1f}")
P("     Λ/L     s·ν bei s = 1     0.3       0.1       0.03      0.01      0.001  |  Boden a/a0 = 2/(2Λ/L−1)   s bei +10 %")
for lr in [1, 10, 100, 1000, np.inf]:
    if np.isinf(lr):
        row = [s*nu(s) for s in [1, 0.3, 0.1, 0.03, 0.01, 0.001]]; fl = "0 (Tully-Fisher)"; s10 = "–"
    else:
        row = [s/(1-q_box(s, lr)) for s in [1, 0.3, 0.1, 0.03, 0.01, 0.001]]; fl = f"{2/(2*lr-1):.4e}"
        s10 = f"{optimize.brentq(lambda s: (1-np.exp(-s))/(1-q_box(s, lr)) - 1.1, 1e-7, 3):.4f}"
    P(f"   {lr:6}   " + "  ".join(f"{x:8.4f}" for x in row) + f"  |  {fl:>18s}        {s10}")
P("   Wurzel für s ≫ L/Λ; darunter a → Boden, v² ∝ r jenseits r ≈ r*·Λ/L.")
P()

P("3b) Schließung 𝒫 = (1/M)Σ T_{mL} auf die abgeschirmte Ausbreitung (−D∂² + Γ): Ring der Länge L gegen Gerade")
P("    G_Ring(0)/G_Gerade(0) = Σ_{n∈ℤ} e^{−|n|s} = coth(s/2) = 1 + 2Σ_{n≥1}qⁿ = 2ν − 1,  s = L√(Γ/D)")
Mr = 400
for s in [0.1, 1.0, 3.0]:
    Gm = 2*(np.cosh(s/Mr) - 1); k = np.arange(Mr)
    Gr = np.mean(1/(Gm + 2*(1 - np.cos(2*np.pi*k/Mr)))); Gl = 1/np.sqrt(Gm*Gm + 4*Gm)
    P(f"   s = {s}:  Ring/Gerade {Gr/Gl:.8f}   coth(s/2) {1/np.tanh(s/2):.8f}   2ν−1 {2*nu(s)-1:.8f}   Γ·G_Ring·M = {Gm*Gr*Mr:.5f} (Nullmode, Kac)")
P("   Gleicher Tiefgrenzfall (Windungslänge L/2): coth s gegen die Regel ν(s)")
P("      s      coth s     ν(s)      Verhältnis")
for s in [0.1, 0.3, 0.5, 1.0, 2.0, 3.0]:
    P(f"   {s:5.2f}   {1/np.tanh(s):8.4f}   {nu(s):8.4f}   {1/np.tanh(s)/nu(s):.4f}")

txt = "\n".join(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rueckkehr_zustandsraum.txt"), "w", encoding="utf-8").write(txt)
print(txt)
