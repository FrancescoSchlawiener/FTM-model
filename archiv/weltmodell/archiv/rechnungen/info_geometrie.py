"""Informationsgeometrie der Rücklaufzahl, Jacobson/Verlinde im Modell, Entropie der Engstelle.
Rechnung, noch nicht im Modell. Konstanten und Formeln aus ftm_kern."""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from ftm_kern import G, c, hbar, kB, l_P, M_sun, A0, H0_si, ds, nu, chi, n_B, schliesslaenge, tau_D, alpha_G, pc

def lnp(n, s):                      # geometrische Verteilung der Rücklaufzahl, q = e^{-s}
    return math.log(-math.expm1(-s)) - n * s

print("1) Rücklaufzahl n ~ (1-q)q^n, q = e^{-s}: Exponentialfamilie, natürlicher Parameter -s, log-Partition psi(s) = ln nu")
print("   Fisher g(s) = Var(n) = n_B(1+n_B) = 1/(4 sinh²(s/2)) = -dnu/ds;  Mittel <n> = n_B = nu-1")
print("   s      g(Formel)   g(Summe über n, Score²)   chi      -dln nu/dln s   H=ln nu+chi  H(Summe)   L(s->inf)=-ln tanh(s/4)  L(Quadratur)")
for s in [0.1, 0.5, 1.0, 2.0, 5.0]:
    g = n_B(s) * (1 + n_B(s))
    h = 1e-5
    N = np.arange(0, 20000)
    p = np.exp([lnp(n, s) for n in N])
    score = np.array([(lnp(n, s + h) - lnp(n, s - h)) / (2 * h) for n in N])
    g_sum = float(np.sum(p * score**2))
    m = p > 0
    H_sum = float(-np.sum(p[m] * np.log(p[m])))
    ela = -(math.log(float(nu((s * (1 + h))**2 * A0))) - math.log(float(nu((s * (1 - h))**2 * A0)))) / (math.log(1 + h) - math.log(1 - h))
    L = -math.log(math.tanh(s / 4))
    Lq = quad(lambda x: math.exp(-x / 2) / (1 - math.exp(-x)), s, 200)[0]
    print(f"   {s:4.1f}  {g:10.6f}  {g_sum:10.6f}                {float(chi(s*s*A0)):.6f}   {ela:.6f}       {math.log(float(nu(s*s*A0)))+float(chi(s*s*A0)):.6f}   {H_sum:.6f}   {L:.6f}                 {Lq:.6f}")

print("\n   Ist chi eine affine Koordinate einer alpha-Verbindung (dann chi' ∝ g^{(1-alpha)/2}, Exponent konstant)?")
print("   s      dln|chi'|/dln g   (alpha=1: 0, Levi-Civita: 0,5, alpha=-1: 1)")
def dchi(s): return n_B(s) * (1 - s * (1 + n_B(s)))
for s in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    e = 1e-4
    gg = lambda x: n_B(x) * (1 + n_B(x))
    k = (math.log(abs(dchi(s + e))) - math.log(abs(dchi(s - e)))) / (math.log(gg(s + e)) - math.log(gg(s - e)))
    print(f"   {s:4.1f}   {k:8.4f}")
print("   -> nicht konstant: chi ist keine Geodäte (kein affiner Parameter) und nicht die Metrik.")
print("   chi = s·<n> = Produkt der dualen Koordinaten (natürlich × Mittel) = -dln nu/dln s; H = psi + chi (Legendre); dH/ds = -s·g")
print("   Identitäten: g = nu(nu-1) = chi·nu/s;  ln nu = KL(Newton-Zustand n=0 || Verteilung bei s)")
s = 1.0
print(f"   Kontrolle s=1: KL(delta_0||p) = -ln p(0) = {-lnp(0, s):.6f}, ln nu = {math.log(float(nu(A0))):.6f}")

print("\n   Zweite Familie: erste Rückkehrzeit t, Dichte f(t) ∝ t^{-3/2}e^{-tau/4t}, gekippt mit e^{-Gamma t}; log-Partition = -s")
tau = 1.0
def f(t): return math.sqrt(tau / (4 * math.pi * t**3)) * math.exp(-tau / (4 * t))
for s in [0.3, 1.0, 3.0]:
    Gam = s * s / tau
    Z = quad(lambda u: f(math.exp(u)) * math.exp(-Gam * math.exp(u)) * math.exp(u), -30, 30, limit=400)[0]
    m1 = quad(lambda u: f(math.exp(u)) * math.exp(-Gam * math.exp(u)) * math.exp(2 * u), -30, 30, limit=400)[0] / Z
    m2 = quad(lambda u: f(math.exp(u)) * math.exp(-Gam * math.exp(u)) * math.exp(3 * u), -30, 30, limit=400)[0] / Z
    gs = (m2 - m1**2) * (2 * s / tau)**2
    print(f"   s={s}: Laplace = {Z:.6f} (e^-s = {math.exp(-s):.6f}), Fisher in s = {gs:.6f} (1/s = {1/s:.6f}); Abstand 0->s = 2√s")
print("   -> auch hier ist chi weder Metrik (1/s) noch Abstand (2√s).")

print("\n2) Jacobson/Verlinde mit der Periode der inneren Richtung")
Lam = schliesslaenge()
TI = hbar * A0 / (c * kB)
H0 = H0_si()
TGH = hbar * H0 / (2 * math.pi * kB)
print(f"   Kac-Periode 2Λ/c = {2*Lam/c:.4e} s = tau_D = c/a0 = {tau_D():.4e} s = hbar/kT_I = {hbar/(kB*TI):.4e} s")
print(f"   T_I = hbar a0/(c k) = {TI:.4e} K;  Gibbons-Hawking hbar H0/(2πk) = {TGH:.4e} K;  Verhältnis {TI/TGH:.4f}")
print("   Verlust Gamma = a_N/c = 2π k T_U(a_N)/hbar (Unruh des ruhenden Beobachters, Fläche des lokalen Rindler-Horizonts)")
for x in [0.01, 1.0, 100.0]:
    aN = x * A0
    TU = hbar * aN / (2 * math.pi * c * kB)
    print(f"   a_N/a0={x:6g}: s² = Gamma·tau_D = {aN/c*tau_D():.6g},  2π T_U/T_I = {2*math.pi*TU/TI:.6g}")
print("   Sattel s = min_t[Gamma t + tau_D/(4t)] (Verlust gegen Pfadentropie der Rückkehr):")
for x in [0.01, 1.0]:
    Gam = x * A0 / c
    r = minimize_scalar(lambda lt: Gam * math.exp(lt) + tau_D() / (4 * math.exp(lt)), bounds=(20, 60), method="bounded")
    print(f"   a_N/a0={x}: min = {r.fun:.6f}, √(a_N/a0) = {math.sqrt(x):.6f}")
M = 1e11 * M_sun
print("   Tiefbereich (Punktmasse 10^11 M☉): v⁴/(G M a0)")
for rk in [100, 300, 1000]:
    r = rk * 1e3 * pc
    aN = G * M / r**2
    v4 = (float(nu(aN)) * aN * r)**2
    print(f"   r={rk:5d} kpc: Regel {v4/(G*M*A0):.4f};  Deser-Levin mit T_I (a = √(2·2πa0·a_N)): {4*math.pi*aN*A0*r**2/(G*M*A0):.4f};"
          f"  Verlinde F=T ΔS/Δx mit √(T_U T_I): {2*math.pi*aN*A0*r**2/(G*M*A0):.4f};  Boltzmann e^(-2πT_U/T_I) ohne Wurzel: v²/(a0 r) = {float(1+n_B(aN/A0))*aN*r/(A0*r):.4f}")
print(f"   Verlinde 2016 [G]: a_M = cH0/6 = {c*H0/6:.3e};  Modell cH0/2π = {c*H0/(2*math.pi):.3e};  Verhältnis {2*math.pi/6:.4f}")

print("\n3) Entropie der Engstelle")
for name, Mx in [("Sonne", M_sun), ("10 M☉", 10 * M_sun), ("4·10^6 M☉", 4e6 * M_sun)]:
    D = ds(Mx)
    rs = 2 * G * Mx / c**2
    S1 = 36 * math.pi * (D / l_P)**2
    S2 = 4 * math.pi * rs**2 / (4 * l_P**2)
    S3 = 4 * math.pi * alpha_G(Mx)
    T1 = hbar * c / (24 * math.pi * kB * D)
    T2 = hbar * c**3 / (8 * math.pi * G * Mx * kB)
    print(f"   {name}: r_s/Δs = {rs/D:.6f}; 36π(Δs/ℓ_P)² = {S1:.6e}; A/4ℓ_P² = {S2:.6e}; 4πα_G = {S3:.6e}; T = {T1:.4e} K, Hawking {T2:.4e} K")
