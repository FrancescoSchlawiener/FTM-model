"""Mikroskala mit der einen Struktur (Betrag, Phase, Verwindung am Knoten) und die zwei Schreibweisen des Nomes
in einem dynamischen Vorgang.  Konstanten aus ftm_kern.py, Kerne aus AME2020 (daten/elemente/mass_1.mas20.txt).
[G] = aus dem Gedächtnis, nicht in den Daten."""
import os, math, itertools
import numpy as np
from scipy.integrate import solve_ivp
import ftm_kern as K

HERE = os.path.dirname(os.path.abspath(__file__))
DAT = os.path.join(HERE, "..", "daten", "elemente")
out = []
def P(s=""):
    out.append(s); print(s)

c, G, hbar, A0 = K.c, K.G, K.hbar, K.A0
m_p, m_e, e = K.m_p, K.m_e, K.e_el
lP, mPl = K.l_P, K.m_Pl
MU_P_G, MU_N_G = 2.792847, -1.913043       # μ/μ_N [G], nicht in den Daten

P("=" * 100)
P("1  TRIADE AUF DER TEILCHENEBENE")
P("=" * 100)
# 1a Hodge auf dem Knoten: Phasenfeld dθ mit Windung n auf dem Torus = harmonischer Teil n/L + exakter Teil
L_, Nx = 1.0, 256
x = np.arange(Nx) / Nx * L_
for n in (1, -1, 2):
    theta = 2 * math.pi * n * x / L_ + 0.3 * np.sin(2 * math.pi * 3 * x / L_)   # Windung n + glatte Störung
    dth = np.diff(np.concatenate([theta, [theta[0] + 2 * math.pi * n]])) * Nx / L_
    harm = dth.mean() * L_ / (2 * math.pi)
    P(f"1a Phase mit Windung {n:+d}: harmonischer Teil (Mittel von dθ)·L/2π = {harm:+.12f}; Rest ist exakt (Gradient einer periodischen Funktion)")
P("   ⇒ Auf dem Knoten ist die Windung (Ladung) genau der harmonische Teil der Phase: Hintergrund = Ladung (topologisch, ohne Länge),")
P("     Zentrum = Gradientfluss GM (Faltung 3Δs), Wirbel = Rotationsteil (Verwindung a_J = J/Mc).")
P()

def leiter(m, J_hbar=0.5):
    ds3 = 3 * K.ds(m); lam = hbar / (m * c); aJ = J_hbar * hbar / (m * c)
    return ds3, lam, aJ
P("1b Leiter je Knoten (Längen in m); Gesetze: 3Δs·λ̄ = ℓ_P², 3Δs·a_J = (J/ħ)·ℓ_P²")
P(f"   {'Knoten':10s} {'3Δs':>11s} {'ℓ_P':>11s} {'a_J':>11s} {'λ̄':>11s}   {'λ̄/ℓ_P = ℓ_P/3Δs':>18s}  {'3Δs·λ̄/ℓ_P²':>12s}  {'r*':>10s}")
for name, m in (("Elektron", m_e), ("Proton", m_p), ("m_p/3", m_p / 3), ("Planck", mPl)):
    ds3, lam, aJ = leiter(m)
    P(f"   {name:10s} {ds3:11.4e} {lP:11.4e} {aJ:11.4e} {lam:11.4e}   {lam/lP:18.6e}  {ds3*lam/lP**2:12.9f}  {K.r_stern(m):10.3e}")
P("   Die Mitte ℓ_P ist für jede Masse dieselbe; jeder Knoten hat genau eine Zahl, (m/m_Pl)² = α_G(m).")
P("   Skalierung: 3Δs ∝ m¹, λ̄ ∝ m⁻¹, a_J ∝ m⁻¹, ℓ_P ∝ m⁰, r* ∝ m^{1/2}·a₀^{-1/2}.")
P("   Jede Gleichung zwischen Sprossen zweier Knoten hat die Form (m_p/m_e)^k = Zahl aus ℓ_P, λ̄, a₀ … ;")
P("   ohne a₀ bleibt nur k·ln(m_p/m_e) = 0 (keine Aussage); mit a₀ = c²/2Λ bände sie m_p/m_e an die Schließlänge Λ (Zustand, Epoche).")
# 1c: Prüfung: alle Paare aus {3Δs, a_J, λ̄, r*} von p und e: Exponent k in m_p/m_e und ob a0 eingeht
spr = {"3Δs": (1, 0), "a_J": (-1, 0), "λ̄": (-1, 0), "r*": (0.5, 1)}
P("1c Paare (Sprosse Proton = Sprosse Elektron): Exponent der Masse und ob a₀ eingeht")
for (a, (ka, xa)), (b, (kb, xb)) in itertools.product(spr.items(), repeat=2):
    if a == b and ka != 0:
        continue
    # m_p^ka · f(a0)^xa = m_e^kb · f(a0)^xb  → nur Fälle mit gleicher Dimension sinnvoll; hier Längen
    P(f"   {a:>4s}(p) = {b:<4s}(e): m_p^{ka:+.1f} ~ m_e^{kb:+.1f}{'  (hängt an a₀)' if xa or xb else ''}")
P("   ⇒ Keine Sprosse legt m_p/m_e fest, ohne eine Größe der inneren Richtung (a₀, Λ) als Teilchengröße zu lesen.")
P()

# 1d Magnetisches Moment: elementarer Knoten g = 2, zusammengesetzter Knoten
P("1d Magnetisches Moment. Elementar: Ladung und Faltung am selben imaginären Versatz i·a_J ⇒ g = 2, μ = Q_el·c·a_J.")
aJp = hbar / (2 * m_p * c)
muN = e * hbar / (2 * m_p)
P(f"   Proton als ein Knoten: μ = e·c·a_J = {e*c*aJp:.6e} J/T = {e*c*aJp/muN:.6f} μ_N")
P(f"   gemessen [G] μ_p = {MU_P_G} μ_N ⇒ Versatz der Ladung a_Q = μ_p/(e c) = {MU_P_G:.4f}·a_J: Betrag und Phase sitzen im Proton nicht am selben Ort.")
# zusammengesetzter Knoten: drei Windungen mit g = 2 je Knoten; Spin-Flavor-Zustand explizit gebaut
up, dn = np.array([1, 0.]), np.array([0, 1.])
def kron(*v):
    r = np.array([1.])
    for a in v: r = np.kron(r, a)
    return r
# Spinzustand ½, m = +½ aus zwei gleichen Knoten (1,2) im Spin 1 und Knoten 3: |χ_S⟩
chiS = (2 * kron(up, up, dn) - kron(up, dn, up) - kron(dn, up, up)) / math.sqrt(6)
chiA = (kron(up, dn, up) - kron(dn, up, up)) / math.sqrt(2)
flav = {"u": np.array([1, 0.]), "d": np.array([0, 1.])}
Qf = np.diag([2 / 3, -1 / 3])
def zustand(f1, f2, f3):
    # total symmetrisch in Spin ⊗ Flavor (Farbe antisymmetrisch, Z₃-neutral): Summe über Permutationen
    phiS = None
    tot = np.zeros(64)
    for perm in itertools.permutations(range(3)):
        fl = [(f1, f2, f3)[i] for i in perm]
        for spin_state in (chiS,):
            # Permutation des Spins entsprechend
            s = spin_state.reshape(2, 2, 2).transpose(perm).reshape(8)
            tot += np.kron(s, kron(flav[fl[0]], flav[fl[1]], flav[fl[2]]))
    return tot / np.linalg.norm(tot)
sz = np.diag([1, -1.])
I2 = np.eye(2)
def op_i(i, A, B):
    ops_s = [I2, I2, I2]; ops_f = [I2, I2, I2]
    ops_s[i] = A; ops_f[i] = B
    S = np.kron(np.kron(ops_s[0], ops_s[1]), ops_s[2]); F = np.kron(np.kron(ops_f[0], ops_f[1]), ops_f[2])
    # Reihenfolge des Zustands: Spin(8) ⊗ Flavor(8)
    return np.kron(S, F)
def summe_Q_sz(psi):
    return sum(psi @ op_i(i, sz, Qf) @ psi for i in range(3))
for name, fl, mess in (("Proton uud", ("u", "u", "d"), MU_P_G), ("Neutron udd", ("u", "d", "d"), MU_N_G)):
    psi = zustand(*fl)
    S3 = sum(psi @ op_i(i, sz, I2) @ psi for i in range(3)) / 2
    QQ = sum(psi @ op_i(i, I2, Qf) @ psi for i in range(3))
    m = summe_Q_sz(psi)          # ⟨Σ Q_i σ_zi⟩
    P(f"   {name}: ⟨ΣS_z⟩ = {S3:+.6f}, ΣQ_el = {QQ:+.6f}, ⟨Σ Q_i σ_zi⟩ = {m:+.6f}  ⇒ μ = ⟨ΣQσ⟩·(m_p/m_q) μ_N;  mit m_q = m_p/3: {3*m:+.4f} μ_N  (gemessen [G] {mess:+.4f})")
psiP, psiN = zustand("u", "u", "d"), zustand("u", "d", "d")
rat = summe_Q_sz(psiN) / summe_Q_sz(psiP)
P(f"   μ_n/μ_p = {rat:+.6f} (hängt nicht an m_q, nur an m_u = m_d); gemessen [G] {MU_N_G/MU_P_G:+.5f}, Abweichung {(rat/(MU_N_G/MU_P_G)-1)*100:+.2f} %")
P(f"   m_q, das μ_p trifft: m_p/{MU_P_G:.4f} = {K.m_p*c**2/e/1e6/MU_P_G:.1f} MeV (= m_p/3 · {3/MU_P_G:.4f})")
P("   Erzwungen im Modell: g = 2 je elementarem Knoten, Summe der Momente (Multipolreihen addieren sich), ΣQ_el ganz (Z₃).")
P("   Nicht erzwungen: Drittelladungen (Eingabe), Farb-Antisymmetrie (volles SU(3), das Modell hat nur das Zentrum Z₃),")
P("   gleiche Teilung der Faltung m_q = m_p/3 (die gebundene Faltung liegt großteils nicht auf den Windungen). Kein Weg zu m_e.")
P()

P("=" * 100)
P("2  KERNE (AME2020): Klasse kompakt/kollektiv und Regel im Kern")
P("=" * 100)
L = open(os.path.join(DAT, "mass_1.mas20.txt"), encoding="latin-1").read().split("\n")
start = next(i for i, l in enumerate(L) if "micro-u" in l) + 1
BA = {}
for l in L[start:]:
    if len(l) < 110: continue
    try:
        N_ = int(l[4:9]); Z_ = int(l[9:14]); A_ = int(l[14:19])
    except ValueError: continue
    bas = l[54:67].strip(); mes = l[28:42].strip()
    if "#" in bas or "#" in mes or "*" in bas or A_ < 2: continue
    b = float(bas)
    BA[A_] = max(BA.get(A_, -1e9), b)
As = np.array(sorted(BA)); Bm = np.array([BA[a] for a in As]) / 1000
imax = int(np.argmax(Bm))
P(f"2a gemessene Isobaren: {len(As)}; größtes B/A je A: Maximum {Bm[imax]:.4f} MeV bei A = {As[imax]}")
for frac in (0.8, 0.9, 0.95):
    a_k = As[np.argmax(Bm >= frac * Bm[imax])]
    P(f"   erstes A mit B/A ≥ {frac:.2f}·max: A = {a_k}")
m_u = 1.66053906660e-27
r0n = 1.2e-15   # [A] R = r₀A^{1/3}
P("2b Kern als Trichter: r*/R, g_N(R)/a₀, Rückkehranteil q = e^{-s}, Kern-Nome Δs/R")
P(f"   {'A':>4s} {'R [fm]':>8s} {'r* [fm]':>9s} {'r*/R':>8s} {'g_N(R)/a₀':>11s} {'s':>8s} {'lg q':>9s} {'Δs/R':>10s}")
for A_ in (2, 4, 12, 16, 56, 62, 208, 238):
    M = A_ * m_u; R = r0n * A_ ** (1 / 3); rs = K.r_stern(M); gN = G * M / R**2; s = math.sqrt(gN / A0)
    P(f"   {A_:4d} {R*1e15:8.3f} {rs*1e15:9.2f} {rs/R:8.2f} {gN/A0:11.3e} {s:8.2f} {-s/math.log(10):9.2f} {K.ds(M)/R:10.3e}")
A_x = (r0n / math.sqrt(G * m_u / A0)) ** 6
P(f"   Grenze kompakt/kollektiv r*/R = 1 (r*/R ∝ A^{{1/6}}): A = {A_x:.2e} < 1 ⇒ jeder Kern ist kompakt, die Regel ruht (dazu Erdfeld: s = {math.sqrt(9.81/A0):.0f}).")
P("   Die Triade der Kerne (Zentrum Faltung, Wirbel Nukleonenbewegung, Hintergrund Vakuum) sagt keinen Übergang in B/A voraus;")
P("   Sättigung, Maximum (x = ¼) und Schalen folgen aus Reichweite der Kernbindung, Coulomb und Zählregeln, nicht aus einer Länge des Modells.")
P()

P("=" * 100)
P("3  NOME: |q_Δs| = Δs/r (instantan) gegen |q| = e^{-r*/r} im dynamischen Vorgang")
P("=" * 100)
def Delta(q, nmax=80):
    return q * np.prod([(1 - q**n) ** 24 for n in range(1, nmax)])
def a_A(M, r):   # Punkt 3: a_N + χ(s)·Δ(q_Δs)·3c²/r*
    aN = G * M / r**2; s = math.sqrt(aN / A0); rs = K.r_stern(M)
    chi = s / math.expm1(s)
    return aN + chi * Delta(K.ds(M) / r) * 3 * c**2 / rs
def a_B(M, r):   # Punkt 1: a_N/(1 − e^{-r*/r})
    aN = G * M / r**2; s = K.r_stern(M) / r
    return aN / (-math.expm1(-s))
def a_N_(M, r): return G * M / r**2
rng = np.random.default_rng(3)
dev = []
for _ in range(2000):
    M = 10 ** rng.uniform(30, 46); rr = K.r_stern(M) * 10 ** rng.uniform(-1.5, 2)
    dev.append(abs(a_A(M, rr) / a_B(M, rr) - 1) / (math.sqrt(G * M * A0) / c**2))
P(f"3a Punktweise, 2000 zufällige (M, r), M = 10³⁰–10⁴⁶ kg, r/r* = 0,03–100: max |a_A/a_B − 1| / (v_f/c)² = {max(dev):.4f} (Schranke 24/(3e) = {24/(3*math.e):.4f})")
P("   Beide Formen sind Funktionen des augenblicklichen (M, r); die Identität χ(s)·(Δs/r)·3c²/r* = a_N·q/(1−q) gilt für jedes M(t), r(t).")

# 3b Innere Richtung: Massenänderung ändert das Profil nicht
P("3b Die eine Gleichung ∂_t u + Γu = D ∂_ϖ²u bei Massensprung M → f·M am festen Ort r:")
for f in (0.7, 2.0):
    M1 = 1e9 * K.M_sun; r = K.r_stern(M1)
    D1, Ga1 = G * M1 / c, G * M1 / (c * r**2)
    D2, Ga2 = f * D1, f * Ga1
    w = np.linspace(0, 8 * r, 4001); h = w[1] - w[0]
    u = np.exp(-w / r)                       # stationäres Profil vor dem Sprung
    lap = (u[2:] - 2 * u[1:-1] + u[:-2]) / h**2
    res_rel = np.max(np.abs(D2 * lap - Ga2 * u[1:-1])) / (Ga2 * 1.0)
    fd = (1 - np.exp(-h / r)) / (h / r)  # Diskretisierung bekannt: Rest ~ h²
    P(f"   f = {f}: Γ/D vorher {Ga1/D1*r**2:.6f}/r², nachher {Ga2/D2*r**2:.6f}/r²; ∂_t u nach dem Sprung / Γ = {res_rel:.2e} (Gitterfehler, h/r = {h/r:.0e})")
    P(f"          Lesetiefe ϖ_R = √(D·2Λ/c): {math.sqrt(D1*2*K.schliesslaenge()/c)/r:.6f} r* → {math.sqrt(D2*2*K.schliesslaenge()/c)/r:.6f} r*  (= √f); q: {math.exp(-1):.6f} → {math.exp(-math.sqrt(f)):.6f}")
P("   ⇒ Γ/D = a_N/GM = 1/r² hängt nicht an M: Ein Massensprung lässt das Profil in der inneren Richtung stationär,")
P("     nur die Lesetiefe ϖ_R ∝ √M rückt, und zwar zugleich mit D. Keine Übergangszeit, auch wenn man die Gleichung zeitlich liest.")

# 3c Konkrete Fälle: Bahnen
def bahn(Mfun, r0_, T, afun):
    rs_ = r0_
    def rhs(t, y):
        x_, y_, vx, vy = y; r = math.hypot(x_, y_); a = afun(Mfun(t), r)
        return [vx, vy, -a * x_ / r, -a * y_ / r]
    v0 = math.sqrt(a_B(Mfun(-1.0), r0_) * r0_)
    sol = solve_ivp(rhs, (0, T), [r0_, 0, 0, v0], method="DOP853", rtol=1e-12, atol=1e-6 * r0_, dense_output=True)
    ts = np.linspace(0, T, 2001); Y = sol.sol(ts)
    return ts, np.hypot(Y[0], Y[1]), np.unwrap(np.arctan2(Y[1], Y[0]))
faelle = [
    ("Zwerg 10⁹ M☉, Stern bei r*, Supernova-Wind nimmt 30 % sofort", 1e9, 1.0, lambda M0: (lambda t: M0 if t < 0 else 0.7 * M0), 3.0),
    ("Zwerg 10⁹ M☉, Stern bei 0,3 r*, Gezeitenabriss −50 % über 0,2 Mrd. J.", 1e9, 0.3, lambda M0: (lambda t: M0 * (1 - 0.5 * min(max(t / (0.2 * K.Gyr), 0), 1))), 3.0),
    ("Haufen 10¹⁴ M☉, Galaxie bei r*, Verschmelzung ×2 über 1 Mrd. J.", 1e14, 1.0, lambda M0: (lambda t: M0 * (1 + min(max(t / K.Gyr, 0), 1))), 10.0),
]
P("3c Bahnen nach der Massenänderung (Punktmasse; Anfang Kreisbahn): (A) Δs/r mit voller Spitzenform, (B) e^{-r*/r}, (N) Newton")
for name, Msun, xr, mk, TG in faelle:
    M0 = Msun * K.M_sun; Mf = mk(M0); r0_ = xr * K.r_stern(M0); T = TG * K.Gyr
    tA, rA, pA = bahn(Mf, r0_, T, a_A)
    tB, rB, pB = bahn(Mf, r0_, T, a_B)
    tN, rN, pN = bahn(Mf, r0_, T, a_N_)
    vf = (G * M0 * A0) ** 0.25
    P(f"   {name}: v_f = {vf/1e3:.1f} km/s, (v_f/c)² = {(vf/c)**2:.2e}")
    P(f"      nach {TG:.0f} Mrd. J.: r_max/r₀ (A) {rA.max()/r0_:.9f}  (B) {rB.max()/r0_:.9f}  (N) {rN.max()/r0_:.6f};  r_min/r₀ (A) {rA.min()/r0_:.9f}  (B) {rB.min()/r0_:.9f}")
    P(f"      max |r_A − r_B|/r₀ = {np.max(np.abs(rA-rB))/r0_:.2e},  |φ_A − φ_B| am Ende = {abs(pA[-1]-pB[-1]):.2e} rad;  |r_B − r_N|/r₀ bis {np.max(np.abs(rB-rN))/r0_:.3f}, |φ_B − φ_N| = {abs(pB[-1]-pN[-1]):.2f} rad")
    # Kontrolle: dieselbe Bahn ohne Massenänderung (Unterschied allein aus der Spitzenform, stationär)
    Mk = lambda t, M0=M0: M0
    tA0, rA0, pA0 = bahn(Mk, r0_, T, a_A); tB0, rB0, pB0 = bahn(Mk, r0_, T, a_B)
    P(f"      Kontrolle ohne Massenänderung: |φ_A − φ_B| am Ende = {abs(pA0[-1]-pB0[-1]):.2e} rad (Umläufe {pB0[-1]/(2*math.pi):.2f}); mit Vorgang {abs(pA[-1]-pB[-1]):.2e} rad (Umläufe {pB[-1]/(2*math.pi):.2f})")
P("   ⇒ Der Vorgang selbst ist groß (gegen Newton Ordnung 1); die beiden Schreibweisen trennen sich nur um die Spitzenform-Glieder ≲ 3(v_f/c)².")
P()

P("=" * 100)
P("Von Hand")
P("=" * 100)
P("   s = 1: χ = 1/(e − 1) = 1/1,718282 = 0,581977;  s·q/(1 − q) = 0,367879/0,632121 = 0,581977")
P("   Proton aus drei Knoten: μ_p = (4μ_u − μ_d)/3, μ_u = +2 μ_N, μ_d = −1 μ_N (m_q = m_p/3)  ⇒ (8 + 1)/3 = 3 μ_N;  μ_n = (4·(−1) − 2)/3 = −2 μ_N;  μ_n/μ_p = −2/3")
P(f"   Zwerg 10⁹ M☉: v_f⁴ = GMa₀ = 1,3271·10²⁹ · 1,1·10⁻¹⁰ = 1,4598·10¹⁹ ⇒ v_f = 61,81 km/s; (v_f/c)² = 4,251·10⁻⁸; 24/(3e)·(v_f/c)² = 1,25·10⁻⁷")
P(f"   Kern ⁵⁶Fe: r* = √(G·56·m_u/a₀) = √(6,6743·10⁻¹¹·9,2990·10⁻²⁶/1,1·10⁻¹⁰) = √(5,642·10⁻²⁶) = 2,375·10⁻¹³ m; R = 1,2·56^{{1/3}} fm = 4,59 fm; r*/R = 51,7")

open(os.path.join(HERE, "mikro_nome.txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
