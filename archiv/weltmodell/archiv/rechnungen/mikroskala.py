"""Mikroskala im FTM: (C) Ladung als Windungszahl, (A) m_p/m_e und die Dimensionstransmutation als Nome,
(B) Feinstrukturkonstante und Modulformen.  Jede Zahl der Notiz 34 stammt aus diesem Skript."""
import math, itertools
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

mp.mp.dps = 30
rng = np.random.default_rng(1)

# --- Konstanten (CODATA 2022 / PDG, aus dem Gedächtnis eingetragen) ---
c, hbar, G = 299792458.0, 1.054571817e-34, 6.67430e-11
h = 2*math.pi*hbar
e, eps0 = 1.602176634e-19, 8.8541878188e-12
ainv0 = 137.035999177          # α(0)^-1, rel. Unsicherheit 1,5e-10
alpha = 1/ainv0
mu_pe = 1836.152673426         # m_p/m_e
me_MeV, mp_MeV = 0.51099895000, 938.27208816
lP = math.sqrt(hbar*G/c**3)
mPl_GeV = math.sqrt(hbar*c/G)*c**2/e/1e9
hbarc_GeVm = hbar*c/e/1e9      # GeV·m
aG = (mp_MeV/1e3/mPl_GeV)**2
aGe = (me_MeV/1e3/mPl_GeV)**2
line = "-"*78

print(line); print("TEIL C  Ladung als Windungszahl"); print(line)

# C1: Windung auf dem Kreis ist ganzzahlig: diskrete Prüfung
print("C1  Windungszahl w = (1/2π)∮dθ eines stetigen Feldes ψ = |ψ|e^{iθ} ≠ 0 auf dem Kreis")
for nu in [0, 1, -2, 0.5, 1/3]:
    y = np.linspace(0, 2*np.pi, 2001)
    psi = np.exp(1j*nu*y)                         # Versuch: Phase läuft um 2πν
    dth = np.angle(psi[1:]/psi[:-1]).sum()/(2*np.pi)
    mism = abs(psi[-1]-psi[0])                    # Randauslenkung |ψ(2π)−ψ(0)|
    print(f"    ν = {nu:+.3f}: Phasenlauf/2π = {dth:+.4f}, Randauslenkung |ψ(L)−ψ(0)| = {mism:.4f}"
          f"  (2|sin πν| = {2*abs(math.sin(math.pi*nu)):.4f})")

# C2: Nicht-ganzzahlige Windung auf dem Torus: Energie der Schnittlinie ∝ Länge L2 (Gitter, Ginzburg-Landau)
print("\nC2  Torus Nx×Ny mit Verdrehung ψ(x+Nx) = e^{2πiν}ψ(x): minimale Energie (Gitter-GL, Gradientenfluss)")
def gl_min(nu, Nx, Ny, K=1.0, lam=1.0, steps=6000, dt=0.05):
    th0 = rng.uniform(-0.3, 0.3, (Nx, Ny))
    psi = np.exp(1j*th0)*(1+0.05*rng.standard_normal((Nx, Ny)))
    tw = np.exp(2j*np.pi*nu)
    def nb(p):                                 # Nachbarn mit verdrehter Randbedingung in x
        xp = np.roll(p, -1, 0); xp[-1] = p[0]*tw
        xm = np.roll(p, 1, 0);  xm[0] = p[-1]/tw
        return xp, xm, np.roll(p, -1, 1), np.roll(p, 1, 1)
    def energy(p):
        xp, _, yp, _ = nb(p)
        return K*(np.abs(xp-p)**2).sum() + K*(np.abs(yp-p)**2).sum() + lam/4*((np.abs(p)**2-1)**2).sum()
    for _ in range(steps):
        xp, xm, yp, ym = nb(psi)
        grad = K*(4*psi - xp - xm - yp - ym) + lam/2*(np.abs(psi)**2-1)*psi
        psi = psi - dt*grad
    return energy(psi)
print(f"    {'ν':>6s} {'Ny=8':>10s} {'Ny=16':>10s} {'Ny=32':>10s} {'E/Ny (32)':>10s} {'analytisch':>11s}")
Nx = 24
for nu in [0.0, 0.25, 1/3, 0.5, 1.0, 2.0]:
    Es = [gl_min(nu, Nx, Ny) for Ny in (8, 16, 32)]
    d = abs(nu-round(nu))
    ana = 2*Nx*(1-math.cos(2*math.pi*d/Nx))  # gleichmäßige Verteilung der Verdrehung, |ψ|≈1, pro Zeile
    print(f"    {nu:6.3f} {Es[0]:10.4f} {Es[1]:10.4f} {Es[2]:10.4f} {Es[2]/32:10.4f} {ana:11.4f}")
print("    -> E = 0 genau für ν ∈ ℤ; sonst E ∝ Ny (Länge der Schnittlinie).")
w_list = [1.0, 0.5, 0.1]
print("    Scharfe Naht der Breite w (keine Verteilung): E = K|δ|²·L2/w, |δ| = 2|sin πν|:")
for w in w_list:
    print(f"      w = {w:4.2f}: σ(ν=1/2) = {4*1.0/w:.1f} K,  σ(ν=1/3) = {3*1.0/w:.1f} K  (divergiert für w→0)")

# C3: Kaluza-Klein
print("\nC3  Kaluza-Klein: q_n = n·√(16πGε₀)·ħ/(c R)   (SI);   Gauß: q_n = 2n√G·ħ/(cR)")
kappa = math.sqrt(16*math.pi*G*eps0)/c
R_KK = math.sqrt(16*math.pi*G*eps0)*hbar/(c*e)
print(f"    R_KK (q_1 = e)        = {R_KK:.4e} m = {R_KK/lP:.3f} ℓ_P   (Formel 2ℓ_P/√α = {2*lP/math.sqrt(alpha):.4e} m)")
print(f"    2πR_KK (Umfang)       = {2*math.pi*R_KK:.4e} m")
m1_GeV = hbarc_GeVm/R_KK
print(f"    Masse der n=1-Mode    = ħ/(cR) = {m1_GeV:.4e} GeV = (√α/2)·m_Pl = {math.sqrt(alpha)/2*mPl_GeV:.4e} GeV")
print(f"    Ladung²/Gravitation der Mode: e²/(4πε₀Gm_1²) = {alpha/(math.sqrt(alpha)/2)**2:.4f}")
ratio_e = alpha/aGe
print(f"    Elektron: e²/(4πε₀Gm_e²) = α/α_G,e = {ratio_e:.4e};  (e/m_e)/(q_1/m_1) = {math.sqrt(ratio_e/4):.4e}")
print(f"    Faltungstiefe der Mode Δs_1 = ℓ_P²/(3R) = {lP**2/(3*R_KK):.4e} m")
print(f"    Umformung: α = (2ℓ_P/R_KK)² = 12·Δs_1/λ̄_1 (λ̄_1 = R_KK) = {12*lP**2/(3*R_KK)/R_KK:.10f}  (1/α = {R_KK**2/(4*lP**2):.6f})")
print("    (Vergleich: α_G = 3·Δs_p/λ̄_p)")

# C4: Dirac / 't Hooft auf dem Torus
print("\nC4  Dirac-Bedingung  q·Φ_m = 2πnħ  (Φ_m = magnetischer Fluss in Wb)")
Phi0 = h/e
print(f"    Flussquant h/e        = {Phi0:.6e} Wb")
print(f"    Gauß: e·g = nħc/2  ->  g_D/e = 1/(2α) = {1/(2*alpha):.3f}")
print(f"    Torus (Fläche A, konstantes B): Kozykel e^{{ieBA/ħ}} = 1  ->  eBA = 2πnħ; c₁ = eΦ/h ∈ ℤ")
for n in [1, 2]:
    phase = np.exp(1j*e*(n*Phi0)/hbar)
    print(f"      Φ = {n}·h/e: Kozykel-Phase = {phase.real:+.6f}{phase.imag:+.6f}i")
ph = np.exp(1j*e*(0.5*Phi0)/hbar)
print(f"      Φ = 0,5·h/e: Kozykel-Phase = {ph.real:+.6f}{ph.imag:+.6f}i  (verletzt)")

# C5: Drittelladungen und Z3
print("\nC5  Z₃-Zentrum: e^{2πiQ}·ω^t = 1 ⇔ Q + t/3 ∈ ℤ   (t = Trialität: Quarks − Antiquarks mod 3)")
states = [("u", Fr(2,3), 1), ("d", Fr(-1,3), 1), ("ū", Fr(-2,3), -1), ("d̄", Fr(1,3), -1),
          ("e⁻", Fr(-1), 0), ("ν", Fr(0), 0), ("W⁺", Fr(1), 0), ("Gluon", Fr(0), 0),
          ("p = uud", Fr(1), 3), ("n = udd", Fr(0), 3), ("π⁺ = ud̄", Fr(1), 0), ("Δ⁺⁺ = uuu", Fr(2), 3),
          ("Diquark uu", Fr(4,3), 2), ("Ω⁻ = sss", Fr(-1), 3)]
for nm, Q, t in states:
    ok = (Q + Fr(t, 3)).denominator == 1
    free = (t % 3 == 0)
    print(f"    {nm:12s} Q = {str(Q):>5s}  t mod 3 = {t%3}  Q+t/3 ganz: {ok}   farbneutral: {free}   Q ganz: {Q.denominator==1}")
sigma_GeV2 = 0.18
print(f"    Saitenspannung QCD σ ≈ {sigma_GeV2} GeV² (Gitter, aus dem Gedächtnis) = {sigma_GeV2/(hbarc_GeVm*1e15):.3f} GeV/fm")

print(line); print("TEIL A  Massenverhältnis und Dimensionstransmutation"); print(line)
mu_ = 2*2.16 + 4.70
print(f"A0  m_p = {mp_MeV:.3f} MeV;  2m_u+m_d (MSbar, 2 GeV) = {mu_:.2f} MeV = {100*mu_/mp_MeV:.2f} % von m_p")
v_H = 246.22
print(f"    m_e = y_e v/√2 mit v = {v_H} GeV  ->  y_e = {math.sqrt(2)*me_MeV/1e3/v_H:.4e}")
print(f"    α_G = (m_p/m_Pl)² = {aG:.4e},  m_Pl = {mPl_GeV:.5e} GeV,  ln(m_Pl/m_p) = {math.log(mPl_GeV/(mp_MeV/1e3)):.3f}")

MZ, asMZ = 91.1876, 0.1180
mc, mb, mt = 1.27, 4.18, 172.7
b0 = lambda nf: 11 - 2*nf/3
b1 = lambda nf: 102 - 38*nf/3
def nf_of(mu): return 3 if mu < mc else 4 if mu < mb else 5 if mu < mt else 6
def run(a, mu0, mu1, loops=2, n=4000):
    """α_s von mu0 nach mu1, Schwellen stetig (MSbar, führende Ordnung der Anpassung)."""
    L0, L1 = math.log(mu0), math.log(mu1)
    edges = sorted({L0, L1} | {math.log(m) for m in (mc, mb, mt) if min(mu0, mu1) < m < max(mu0, mu1)})
    if L1 < L0: edges = edges[::-1]
    for A, B in zip(edges[:-1], edges[1:]):
        nf = nf_of(math.exp((A+B)/2)); dL = (B-A)/n
        f = lambda x: -(b0(nf)/(2*math.pi))*x**2 - (loops == 2)*(b1(nf)/(8*math.pi**2))*x**3
        for _ in range(n):
            k1 = f(a); k2 = f(a+dL*k1/2); k3 = f(a+dL*k2/2); k4 = f(a+dL*k3)
            a += dL*(k1+2*k2+2*k3+k4)/6
            if not (0 < a < 5): return float('inf')
    return a
MGUT = 2e16
print("\nA1  Lauf von α_s(M_Z) = 0,1180 nach oben (Schwellen m_c, m_b, m_t)")
res = {}
for loops in (1, 2):
    aPl = run(asMZ, MZ, mPl_GeV, loops); aGU = run(asMZ, MZ, MGUT, loops)
    res[loops] = aPl
    print(f"    {loops}-Schleifen: α_s(M_GUT=2e16) = {aGU:.5f} (1/α = {1/aGU:.2f});  α_s(M_Pl) = {aPl:.5f} (1/α = {1/aPl:.2f})")

# Nome: Λ = M·exp(−2π/(β₀α_s(M)))  (1-Schleife, feste n_f)
print("\nA2  Nome der Dimensionstransmutation: Λ/M = |q̃| = e^{−2π Im τ̃},  Im τ̃ = 1/(β₀ α_s(M))")
a_mc = run(asMZ, MZ, mc, 2)
Lam3_1l = mc*math.exp(-2*math.pi/(b0(3)*a_mc))
print(f"    α_s(m_c) = {a_mc:.4f};  Λ₃ (1-Schleifen-Form an m_c) = {Lam3_1l*1e3:.1f} MeV;  C = m_p/Λ₃ = {mp_MeV/1e3/Lam3_1l:.3f}")
a_mc1 = run(asMZ, MZ, mc, 1)
Lam3_pure = mc*math.exp(-2*math.pi/(b0(3)*a_mc1))
# rein einschleifig mit Schwellen: RG-invariant, gleiches Λ₃ aus M_Pl
aPl1 = res[1]
b = lambda nf: b0(nf)/(2*math.pi)
lnL3 = math.log(mc) - (1/aPl1 + b(6)*math.log(mt/mPl_GeV) + b(5)*math.log(mb/mt) + b(4)*math.log(mc/mb))/b(3)
print(f"    rein 1-Schleife: Λ₃ aus α_s(m_c) = {Lam3_pure*1e3:.2f} MeV, aus α_s(M_Pl) mit Schwellenformel = {math.exp(lnL3)*1e3:.2f} MeV")
C1 = mp_MeV/1e3/Lam3_pure
print(f"    C₁ = m_p/Λ₃(1-Schleife) = {C1:.3f}")
print(f"    dlnΛ₃/d(1/α_s(M_Pl)) = −2π/β₀(3) = {-2*math.pi/9:.4f}  (exakt, 1 Schleife)")
qt = math.exp(lnL3)/mPl_GeV
print(f"    |q̃|(M_Pl) = Λ₃/M_Pl = {qt:.4e},  Im τ̃_eff = ln(M_Pl/Λ₃)/2π = {math.log(1/qt)/(2*math.pi):.4f}")
print(f"    Instanton-Nome τ = θ/2π + i/α_s:  |q| = e^{{−2π/α_s(M_Pl)}} = 10^{-2*math.pi/aPl1/math.log(10):.1f}")
print(f"    α_G = C²|q̃|² = {C1**2*qt**2:.4e}   (gemessen {aG:.4e})")
for nf in (3, 6):
    print(f"    Einfache Form exp(−4π/(β₀α_s(M_Pl))) mit n_f={nf} (β₀={b0(nf):.0f}): {math.exp(-4*math.pi/(b0(nf)*aPl1)):.3e}")
beff = 4*math.pi/(aPl1*math.log(1/aG))
beff2 = 4*math.pi/(res[2]*math.log(1/aG))
print(f"    effektives β₀ aus α_G = exp(−4π/(β_eff α_s(M_Pl))): {beff:.3f} (1-Schl.), {beff2:.3f} (2-Schl.)")

print("\nA3  Empfindlichkeit: gleiche Dynamik (C fest), andere Kopplung am Planck-Punkt (2 Schleifen)")
def mu_half(aPl):
    lo, hi = math.log(0.05), math.log(mPl_GeV)     # Skala mit α_s = 0,5 als Λ-Stellvertreter
    a_of = lambda L: run(aPl, mPl_GeV, math.exp(L), 2, 800)
    for _ in range(50):
        mid = (lo+hi)/2
        if a_of(mid) > 0.5: lo = mid
        else: hi = mid
    return math.exp((lo+hi)/2)
aPl2 = res[2]; ref = mu_half(aPl2)
print(f"    α_s(M_Pl) = {aPl2:.5f} -> μ(α_s=0,5) = {ref*1e3:.1f} MeV, C' = m_p/μ = {mp_MeV/1e3/ref:.3f}")
for f in (0.95, 0.98, 1.02, 1.05):
    m = mp_MeV*mu_half(aPl2*f)/ref
    print(f"    α_s(M_Pl)·{f:.2f} = {aPl2*f:.5f}: m_p = {m:10.2f} MeV,  α_G = {(m/1e3/mPl_GeV)**2:.3e}")
# gesuchtes α_s(M_Pl) für m_p: per Konstruktion das gelaufene; dazu Umkehrprobe
print(f"    -> das gemessene m_p verlangt α_s(M_Pl) = {aPl2:.5f} (1/α_s = {1/aPl2:.2f}); ±1 % in α_s(M_Pl) ändert m_p um"
      f" {100*(mu_half(aPl2*1.01)/ref-1):+.1f} % / {100*(mu_half(aPl2*0.99)/ref-1):+.1f} %")

print("\nA4  Zufallsprüfung für m_p/m_e = 1836,15267 (Toleranz 2e-5 relativ)")
lenz = 6*math.pi**5
print(f"    6π⁵ = {lenz:.4f}, Abweichung {(lenz/mu_pe-1):+.3e} relativ = {abs(lenz/mu_pe-1)/6.0e-11:.1e} σ (m_p/m_e rel. Unsicherheit 6e-11)")
rats = sorted({Fr(p, q) for p in range(1, 13) for q in range(1, 13)})
A_exp = [Fr(k, 2) for k in range(-16, 17)]
B_exp = list(range(-4, 5))
vals, labels = [], []
for r in rats:
    for a in A_exp:
        for bb in B_exp:
            vals.append(math.log(r.numerator/r.denominator) + float(a)*math.log(math.pi) + bb)
            labels.append((r, a, bb))
vals = np.array(vals); order = np.argsort(vals); sv = vals[order]
uniq = np.unique(np.round(sv, 12))
def hits(target, tol):
    lt = math.log(target); lo = np.searchsorted(uniq, lt - tol); hi = np.searchsorted(uniq, lt + tol)
    return hi - lo
print(f"    Ausdrücke r·π^a·e^b: r = p/q (p,q ≤ 12), a ∈ {{−8, −7,5, …, 8}}, b ∈ {{−4..4}}: {len(vals)} Formen, {len(uniq)} verschiedene Werte")
lt = math.log(mu_pe)
hl = [labels[i] for i in range(len(vals)) if abs(vals[i]-lt) < 2e-5]
print(f"    Treffer bei 1836,15267: {len(hl)} Formen; {hits(mu_pe, 2e-5)} verschiedene Werte")
for (r, a, bb) in hl:
    v = float(r)*math.pi**float(a)*math.e**bb
    print(f"      {r}·π^{float(a)}·e^{bb} = {v:.4f}  ({v/mu_pe-1:+.2e})")
targets = np.exp(rng.uniform(math.log(1000), math.log(3000), 20000))
H = np.array([hits(t, 2e-5) for t in targets])
print(f"    Zufallsziele in [1000, 3000] (20000): mittlere Trefferzahl {H.mean():.3f}; Anteil mit ≥1 Treffer {100*(H>0).mean():.1f} %")

print("\nA5  Eigenwerte auf dem Quadrat-Torus: λ ∝ n₁²+n₂²,  m ∝ √λ")
N0 = mu_pe**2
print(f"    (m_p/m_e)² = {N0:.2f};  √{round(N0)} = {math.sqrt(round(N0)):.6f} ({math.sqrt(round(N0))/mu_pe-1:+.2e})")
def is_sum2(n):
    k = 2
    while k*k <= n:
        if n % k == 0:
            cnt = 0
            while n % k == 0: n //= k; cnt += 1
            if k % 4 == 3 and cnt % 2: return False
        k += 1
    return not (n % 4 == 3)
lo, hi = int(N0*(1-4e-5)), int(N0*(1+4e-5))
S = [n for n in range(lo, hi+1) if is_sum2(n)]
print(f"    Summen zweier Quadrate N mit |√N/1836,15267 − 1| < 2e-5: {len(S)} (von {hi-lo+1} ganzen Zahlen im Fenster)")
print(f"    nächstes: N = {min(S, key=lambda n: abs(n-N0))}, √N = {math.sqrt(min(S, key=lambda n: abs(n-N0))):.6f}")

print("\nA6  Geometrische Lesarten von 6π⁵: k·V₁ bzw. k·V₁·V₂ bzw. k·V₁/V₂ = 6π⁵ exakt, k ganz ≤ 1000")
def vol_sphere(n):   # Vol(S^n) = 2π^{(n+1)/2}/Γ((n+1)/2) = rat·π^{m}
    m2 = n+1   # Exponent 2·(n+1)/2
    # Γ(k/2): für gerades k rational, für ungerades k rational·√π
    k = m2
    if k % 2 == 0:
        g = Fr(math.factorial(k//2 - 1)); pe = Fr(k, 2)
    else:
        g = Fr(1); x = Fr(1, 2)
        while x < Fr(k, 2): g *= x; x += 1
        pe = Fr(k, 2) - Fr(1, 2)
    return (Fr(2)/g, pe)
def vol_ball(n):
    r, pe = vol_sphere(n-1); return (r/n, pe)
objs = {}
for n in range(1, 13):
    objs[f"S^{n}"] = vol_sphere(n); objs[f"B^{n}"] = vol_ball(n); objs[f"T^{n}"] = (Fr(2**n), Fr(n))
target = (Fr(6), Fr(5))
reps = []
names = list(objs)
def chk(r, pe, desc):
    if pe == 5:
        k = Fr(6)/r
        if k.denominator == 1 and 1 <= k <= 1000: reps.append((int(k), desc))
for a in names:
    chk(*objs[a], a)
for a, b2 in itertools.combinations_with_replacement(names, 2):
    ra, pa = objs[a]; rb, pb = objs[b2]; chk(ra*rb, pa+pb, f"{a}·{b2}")
for a in names:
    for b2 in names:
        if a != b2:
            ra, pa = objs[a]; rb, pb = objs[b2]; chk(ra/rb, pa-pb, f"{a}/{b2}")
print(f"    {len(reps)} Darstellungen; Beispiele (kleinstes k zuerst):")
for k, d in sorted(reps)[:10]:
    print(f"      6π⁵ = {k}·Vol({d})")
chk_ball = vol_ball(10); print(f"    Kontrolle: Vol(B¹⁰) = {chk_ball[0]}·π^{chk_ball[1]}, Vol(S⁹) = {vol_sphere(9)[0]}·π^{vol_sphere(9)[1]}")

print(line); print("TEIL B  Feinstrukturkonstante"); print(line)
wyler = (9/(8*math.pi**4))*(math.pi**5/(2**4*math.factorial(5)))**0.25
print(f"B1  Wyler: α = (9/8π⁴)(π⁵/(2⁴·5!))^{{1/4}} = 1/{1/wyler:.6f};  1/α gemessen {ainv0};"
      f"  Abweichung {(1/wyler)/ainv0-1:+.3e} = {((1/wyler)/ainv0-1)/1.5e-10:.0f} σ")
print(f"    α(M_Z) (MSbar): 1/127,95;  Unterschied zu α(0): Δ(1/α) = {ainv0-127.95:.2f}")
mlep = [me_MeV/1e3, 0.1056583755, 1.77686]
dlep = sum(alpha/(3*math.pi)*(math.log(MZ**2/m**2) - 5/3) for m in mlep)
print(f"    Leptonanteil Δα_lep(M_Z) (1 Schleife) = {dlep:.5f};  1/α nach Leptonen allein = {1/(alpha/(1-dlep)):.2f}")

# Modulformen an speziellen Punkten
def qexp(tau): return mp.exp(2j*mp.pi*tau)
def eta(tau): q = qexp(tau); return mp.exp(2j*mp.pi*tau/24)*mp.qp(q)
def Ek(tau, k, N=200):
    q = qexp(tau); c0 = {2: -24, 4: 240, 6: -504}[k]
    s = mp.mpf(0)
    for n in range(1, N):
        s += n**(k-1)*q**n/(1-q**n)
    return 1 + c0*s
pts = {"i": 1j, "ρ": mp.exp(2j*mp.pi/3), "i√2": 1j*mp.sqrt(2), "i√3": 1j*mp.sqrt(3), "2i": 2j,
       "(1+i√7)/2": (1+1j*mp.sqrt(7))/2, "(1+i√11)/2": (1+1j*mp.sqrt(11))/2, "(1+i√19)/2": (1+1j*mp.sqrt(19))/2,
       "(1+i√43)/2": (1+1j*mp.sqrt(43))/2, "(1+i√67)/2": (1+1j*mp.sqrt(67))/2, "(1+i√163)/2": (1+1j*mp.sqrt(163))/2,
       "i√5": 1j*mp.sqrt(5), "3i": 3j}
base = {}
for nm, t in pts.items():
    et = eta(t); e2, e4, e6 = Ek(t, 2), Ek(t, 4), Ek(t, 6)
    D = et**24; j = e4**3/D
    for lab, v in [("|η|", abs(et)), ("|E₂|", abs(e2)), ("|E₄|", abs(e4)), ("|E₆|", abs(e6)), ("|j|", abs(j)), ("|q|", abs(qexp(t)))]:
        if v > 1e-12 and abs(v-1) > 1e-12: base[f"{lab}({nm})"] = float(v)
base["Γ(1/4)"] = float(mp.gamma(0.25)); base["Γ(1/3)"] = float(mp.gamma(mp.mpf(1)/3))
base["ϖ (Lemniskate)"] = float(mp.gamma(0.25)**2/(2*mp.sqrt(2*mp.pi)))
print(f"\nB2  Kontrolle: η(i) = {float(abs(eta(1j))):.10f},  Γ(1/4)/(2π^{{3/4}}) = {float(mp.gamma(0.25)/(2*mp.pi**0.75)):.10f};"
      f"  E₂(i) = {float(Ek(1j,2).real):.10f} = 3/π = {3/math.pi:.10f};  j(i) = {float((Ek(1j,4)**3/eta(1j)**24).real):.4f}")
print(f"    E₂(ρ) = {complex(Ek(pts['ρ'],2))}  (2√3/π = {2*math.sqrt(3)/math.pi:.10f})")
print(f"    j((1+i√163)/2) = {float(mp.re(Ek(pts['(1+i√163)/2'],4)**3/eta(pts['(1+i√163)/2'])**24)):.6e} (−640320³ = {-640320**3:.6e})")
bn = list(base); bv = np.array([base[k] for k in bn])
pows = [1, 2, 3, 4, 6, 8, 12, 24]; pows = pows + [-p for p in pows] + [1/p for p in pows[1:8]] + [-1/p for p in pows[1:8]]
pis = range(-4, 5)
rq = sorted({Fr(p, q) for p in range(1, 13) for q in range(1, 13)})
lr = np.array([math.log(r.numerator/r.denominator) for r in rq])
logs = []
for i, v in enumerate(bv):
    for p in pows:
        for k in pis:
            logs.append(p*math.log(v) + k*math.log(math.pi) + lr)
logs = np.concatenate(logs)
Ntot = logs.size
uq = np.unique(np.round(np.sort(logs), 11))
print(f"\nB3  Zufallsprüfung Modulformen: {len(bn)} Grundwerte × {len(pows)} Potenzen × {len(pis)} π-Potenzen × {len(rq)} Brüche"
      f" = {Ntot} Formen, {uq.size} verschiedene Werte")
def cnt(t, tol):
    L = math.log(t); return np.searchsorted(uq, L+tol) - np.searchsorted(uq, L-tol)
for tol in (1e-4, 1.5e-10):
    print(f"    Toleranz {tol:.1e}: Treffer bei 137,035999 = {cnt(ainv0, tol)}")
T = np.exp(rng.uniform(math.log(100), math.log(200), 20000))
Hm = np.array([cnt(t, 1e-4) for t in T])
print(f"    Zufallsziele in [100, 200]: mittlere Trefferzahl bei 1e-4 = {Hm.mean():.2f};  erwartet bei 1,5e-10: {Hm.mean()*1.5e-6:.2e}")
# Trefferliste 1e-4 (einige)
hitsl = []
L = math.log(ainv0)
for i, v in enumerate(bv):
    for p in pows:
        for k in pis:
            arr = p*math.log(v) + k*math.log(math.pi) + lr
            for j in np.where(np.abs(arr - L) < 1e-4)[0]:
                hitsl.append((abs(arr[j]-L), bn[i], p, k, rq[j]))
hitsl.sort()
print(f"    beste Treffer (von {len(hitsl)} Formen):")
for dv, nm, p, k, r in hitsl[:5]:
    print(f"      {r}·π^{k}·[{nm}]^{p:.4g}: rel. Abw. {dv:.2e}")

print("\nB4  Modul-Fixpunkte und E₂-Nullstelle, τ = θ/2π + i/α (Montonen-Olive-Normierung, 4π/e² = 1/α)")
print(f"    τ = i  (S-Fixpunkt): α = 1;   τ = ρ (ST-Fixpunkt): Im ρ = √3/2 -> α = 2/√3 = {2/math.sqrt(3):.4f}, θ = −π")
tz = mp.findroot(lambda y: mp.re(Ek(1j*y, 2)), 0.52)
print(f"    E₂(iy) = 0 bei y = {float(tz):.10f}:  α = 1/y = {1/float(tz):.4f};  |q| = e^{{−2πy}} = {math.exp(-2*math.pi*float(tz)):.5f} = 1/{math.exp(2*math.pi*float(tz)):.3f}")
yq = math.log(ainv0)/(2*math.pi)
print(f"    |q| = α verlangt Im τ = ln(137,036)/2π = {yq:.5f}; j(i·{yq:.4f}) = {float(mp.re(Ek(1j*yq,4)**3/eta(1j*yq)**24)):.3f}, E₂ = {float(mp.re(Ek(1j*yq,2))):.5f}")
b_qed = 4/3*(3*1 + 3*3*(4/9) + 3*3*(1/9))
Lp = 2*math.pi*(127.95-1)/b_qed
print(f"    Skala von α = 1 (QED-artiger Lauf, alle SM-Fermionen, b = {b_qed:.3f}): ln(μ/M_Z) = {Lp:.1f}, μ ≈ 10^{math.log10(MZ)+Lp/math.log(10):.1f} GeV (M_Pl = 10^{math.log10(mPl_GeV):.1f})")

print("\nB5  Carter nach Erklärung von α_G:  α ≈ α_G^{1/12}·(m_p/m_e)^{1/3}")
aC = aG**(1/12)*mu_pe**(1/3)
print(f"    = 1/{1/aC:.2f}  (gemessen 1/{ainv0:.3f});  Faktor {alpha/aC:.4f};  α_G/(α¹²μ⁴) = {aG/(alpha**12*mu_pe**-4):.3f}")

print("\nB6  Vereinheitlichung (1 Schleife): 1/α_em = 1/α_Y + 1/α_2, α_1 = (5/3)α_Y; bei M_X: 1/α_em = (8/3)/α_X")
s2 = 0.23122
a1i, a2i, a3i = (3/5)*(1-s2)*127.95, s2*127.95, 1/asMZ
print(f"    M_Z: 1/α_1 = {a1i:.2f}, 1/α_2 = {a2i:.2f}, 1/α_3 = {a3i:.2f}")
bSM = (41/10, -19/6, -7); bMS = (33/5, 1, -3)
for nm, bb in (("SM", bSM), ("MSSM (ab M_Z, hypothetisch)", bMS)):
    t12 = 2*math.pi*(a1i-a2i)/(bb[0]-bb[1]); t23 = 2*math.pi*(a2i-a3i)/(bb[1]-bb[2])
    print(f"    {nm}: α_1=α_2 bei 10^{math.log10(MZ)+t12/math.log(10):.2f} GeV, α_2=α_3 bei 10^{math.log10(MZ)+t23/math.log(10):.2f} GeV;"
          f"  1/α_X ≈ {a2i - bb[1]/(2*math.pi)*t23:.2f}")
# Vorhersage α_em(M_Z) aus α_s, sin²θ_W und MSSM-Vereinigung (3 Gleichungen, Unbekannte α_X, M_X, α_em)
# 1/α_1 = (3/5)(1-s2)/a, 1/α_2 = s2/a  mit a = α_em;  Bedingung: Schnittpunkte gleich
b1_, b2_, b3_ = bMS
# (1/α_1 - 1/α_2)/(b1-b2) = (1/α_2 - 1/α_3)/(b2-b3)
# X/a = ((3/5)(1-s2) - s2)/(b1-b2)  und (s2/a - a3i)/(b2-b3)
A_ = ((3/5)*(1-s2) - s2)/(b1_-b2_); B_ = s2/(b2_-b3_); C_ = a3i/(b2_-b3_)
ainv_pred = C_/(B_-A_)
print(f"    MSSM-Vorhersage aus α_s(M_Z) und sin²θ_W: 1/α_em(M_Z) = {ainv_pred:.2f}  (MSbar gemessen 127,95)")
print("\nFertig.")
