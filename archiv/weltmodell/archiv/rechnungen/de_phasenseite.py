"""Phasenseite des Ganzen und Dunkle Energie.

Frage: Traegt die Wicklung W der Zelle (Phase des Nomes des Ganzen) zum Vakuumteil bei,
Omega_s(W) = (pi^2/45) Theta^2 * R(W),  R(W) = Sum cos(2 pi n W)/n^4 / zeta(4) = -30 B_4(W)?
Kann das die uebrigen zwei Drittel geben, und welches w(z) folgt? Wo geht q -> 1 im Ganzen?
Alle Konstanten aus ftm_kern.py. Messwerte w0, w_a nur als Test [G].
"""
import math
import numpy as np
import ftm_kern as K

zeta4 = math.pi**4 / 90
PHI = (math.sqrt(5) - 1) / 2          # goldene Wicklung (mod 1: 0,618 bzw. 0,382)


def B4(x):
    return x**4 - 2 * x**3 + x**2 - 1 / 30


def R_reihe(W, n=200000):
    k = np.arange(1, n + 1, dtype=float)
    return float(np.sum(np.cos(2 * np.pi * k * W) / k**4) / zeta4)


def R_bern(W):
    return -30 * B4(W % 1.0)


def Li4_betrag(W, n=200000):
    k = np.arange(1, n + 1, dtype=float)
    z = np.sum(np.exp(2j * np.pi * k * W) / k**4)
    return abs(z) / zeta4, z.imag / zeta4


print("=" * 100)
print("1  Phase je Ruecklauf im Vakuumteil: R(W) = Sum cos(2 pi n W)/n^4 / zeta(4) = -30 B4(W) = 1 - 30 W^2 (1-W)^2")
print("=" * 100)
Om_s = K.omega_vakuum()
Om_DE = K.OMEGA_DE
print(f"   Omega_s (W = 0) = {Om_s:.4f},  Omega_DE = {Om_DE:.4f},  noetig R = Omega_DE/Omega_s = {Om_DE / Om_s:.3f}")
print("   W        R Reihe     R Bernoulli   |Li4|/zeta4   Im/zeta4   Omega_s*R")
for W in [0.0, 0.05, 0.1, 0.2, 1 - PHI, 0.25, 1 / 3, 0.4, 0.5, PHI, 1 / math.sqrt(2) % 1, math.sqrt(2) % 1]:
    r1, r2 = R_reihe(W), R_bern(W)
    ab, im = Li4_betrag(W)
    print(f"   {W:.4f}   {r1:+.6f}   {r2:+.6f}     {ab:.6f}      {im:+.4f}   {Om_s * r1:+.4f}")
Ws = np.linspace(0, 1, 100001)
Rs = R_bern(Ws)
print(f"   max R = {Rs.max():.6f} bei W = {Ws[Rs.argmax()]:.3f};  min R = {Rs.min():.6f} bei W = {Ws[Rs.argmin()]:.3f}")
W0 = 0.5 - math.sqrt(0.25 - 1 / math.sqrt(30))
print(f"   R = 0 bei W(1-W) = 1/sqrt(30): W = {W0:.4f} bzw. {1 - W0:.4f}; R < 0 fuer {W0:.4f} < W < {1 - W0:.4f} (golden liegt darin)")
print("   Dreiecksungleichung: |Sum e^{2 pi i n W}/n^4| <= zeta(4) fuer jedes W; die Phase kann den Betrag nur senken.")
print(f"   Golden: Omega_s*R = {Om_s * R_bern(PHI):+.4f}; Rest der Dunklen Energie dann {Om_DE - Om_s * R_bern(PHI):.4f} statt {Om_DE - Om_s:.4f}")

print()
print("   Beide Umlaufsinne (+W und -W) geben den Realteil; ein Sinn allein gibt Li4 komplex (Im != 0 ist keine Energie).")
print("   Bose-Form mit Phase: Int s^3/(e^{s - 2 pi i W} - 1) ds = Gamma(4) Li4(e^{2 pi i W}), Pruefung numerisch:")
from scipy.integrate import quad
for W in [0.0, 1 - PHI, 0.5]:
    f_re = lambda s: (s**3 * (np.exp(s) * math.cos(2 * math.pi * W) - 1) / (np.exp(2 * s) - 2 * np.exp(s) * math.cos(2 * math.pi * W) + 1))
    I, _ = quad(f_re, 1e-12, 80, limit=400)
    print(f"      W = {W:.4f}: Re Int = {I:.6f},  Gamma(4)*zeta4*R = {6 * zeta4 * R_bern(W):.6f}")

print()
print("=" * 100)
print("1b Wicklung als Verdrehung eines Gitters (Re tau), Epstein-Summe s = 2: Sum' (Im tau)^2/|m + n tau|^4")
print("   Chowla-Selberg: Re tau steht nur in Gliedern mit e^{-2 pi |mn| Im tau}; der fuehrende Teil 2 zeta(4) (Im tau)^2 haengt nicht an Re tau.")
print("=" * 100)


def epstein(tau, N=400):
    m = np.arange(-N, N + 1)
    M, Nn = np.meshgrid(m, m)
    z = M + Nn * tau
    a = np.abs(z)**4
    a[N, N] = np.inf
    return float(np.sum(tau.imag**2 / a))


L = K.zelle_aus_a0(K.A0)
Lam = K.schliesslaenge()
for lab, y in [("Im tau = 1 (quadratisch)", 1.0), ("Im tau = L/2Lambda", L / (2 * Lam)), ("Im tau = 2Lambda/L", 2 * Lam / L), ("Im tau = 0,5", 0.5)]:
    E0 = epstein(complex(0.0, y))
    Eg = epstein(complex(1 - PHI, y))
    print(f"   {lab:26s} y = {y:.4f}: E(Re=0) = {E0:.6f}, E(Re=golden) = {Eg:.6f}, rel. Aenderung {Eg / E0 - 1:+.2e}, e^(-2 pi y) = {math.exp(-2 * math.pi * y):.2e}")
print(f"   (L = {L / K.Mpc / 1e3:.2f} Gpc, 2 Lambda = {2 * Lam / K.Mpc / 1e3:.2f} Gpc)")

print()
print("=" * 100)
print("2  Zeitverlauf: W = omega1/omega2 bei gleichem Atmen und bei neuer Faltung")
print("=" * 100)
aa = np.geomspace(1e-3, 10, 7)
L1, L2 = 1.0, PHI
W_t = (L2 * aa) / (L1 * aa)    # Takte omega_i ~ 1/L_i; beide atmen mit a
print(f"   gleiches Atmen, W(a) = {np.array2string(W_t, precision=12)}  -> fest (Streuung {np.ptp(W_t):.1e})")
eps, S, H_end = K.kreislauf_mindestrate(Lam)
H0 = K.H0_si()
print(f"   Mindestrate eps = {eps:.3e} 1/s (Kern), eps/H0 = {eps / H0:.2e}")
Rg = R_bern(PHI)
dR = -30 * (4 * PHI**3 - 6 * PHI**2 + 2 * PHI)
print(f"   golden: R = {Rg:.4f}, dR/dW = {dR:.4f}, |dlnR/dW| = {abs(dR / Rg):.3f}")
print(f"   selbst wenn neue Faltung W mit voller Rate eps verschiebt (dW/dln a <= W*eps/H):")
bound = abs(dR / Rg) * PHI * eps / H0 / 3
print(f"   |1 + w_s| <= (1/3)|dlnR/dW| W eps/H0 = {bound:.1e}  -> w0 = -1, w_a = 0 (Test [G]: w0 ~ -0,75...-0,85, w_a ~ -0,6...-0,9)")
print("   Ein Drittel mit festem R(W) aendert nur die Hoehe des Vakuumteils, nicht die Form: w = -1 zu jeder Zeit nach z ~ 10.")

print()
print("=" * 100)
print("3  Wo geht q -> 1 im Ganzen? q = e^{-s}, s = sqrt((Gamma + theta~) tau_D), Gamma = a_N/c")
print("=" * 100)
tD = K.tau_D()
th_c = 1 / tD
print(f"   theta_c = a0/c = {th_c:.3e} 1/s,  tau_D = {tD / K.Gyr:.1f} Gyr")
print("   s -> 0 verlangt a_N -> 0 (Sattel, Zentrum) UND theta~ -> 0 (gebunden, nach der ersten Umkehr).")
print("   Regel am Grenzwert: (1 - q) Sum q^n = 1 (Projektor P), a = a_N/(1 - q):")
for s in [1.0, 0.1, 0.01, 1e-3, 1e-4]:
    aN = s**2 * K.A0
    q = math.exp(-s)
    print(f"      s = {s:7.0e}: Sum q^n = {1 / (1 - q):10.3e},  a = {aN / (1 - q):.3e} m/s^2,  sqrt(a_N a0) = {math.sqrt(aN * K.A0):.3e},  a/a0 = {aN / (1 - q) / K.A0:.3e}")
print("   Die Summe divergiert, die Beschleunigung geht gegen null: q -> 1 ist der Grenzwert der Rueckkehr (P), keine Falte.")
print("   Kippstelle (06) verlangt Schleifenverstaerkung q >= 1 (Durchgehen); |q| = e^{-s} <= 1 fuer s >= 0: nie ueberschritten.")
# Neue Faltung aus Bindung? Lambda = nu_kosm * Sum Delta s_b, Delta s = GM/3c^2 additiv
Mb_Zelle = K.rho_m_mitbewegt() * K.OMEGA_B / K.OMEGA_M * L**3
Lam_check = (1 + K.Q_KOSM) * K.G * Mb_Zelle / (3 * K.c**2)
print(f"   Lambda = (1+Q) G M_b(Zelle)/(3c^2) = {Lam_check / K.Mpc / 1e3:.2f} Gpc (Kern c^2/2a0 = {Lam / K.Mpc / 1e3:.2f} Gpc)")
print("   Delta s ist additiv (Sum Delta s = G Sum M/3c^2): Bindung, Saettel und q -> 1 aendern Sum Delta s_b nicht; Rate aus q -> 1: null.")
print(f"   Vergleich der Raten: naechste Rate der Rueckkehr theta_c = {th_c:.2e} 1/s, Mindestrate eps = {eps:.2e} 1/s, eps/theta_c = {eps / th_c:.1e}")
print(f"   Kritische Verlangsamung am Sattel: Rate theta_c*(1-q) ~ theta_c*s -> 0; bei s = 1e-3: {th_c * 1e-3:.1e} 1/s")

print()
print("=" * 100)
print("4  Von Hand")
print("=" * 100)
Wg = 0.381966
x = Wg * (1 - Wg)
print(f"   W = 0,381966: W(1-W) = 0,381966*0,618034 = {x:.6f} (= sqrt5 - 2 = 0,236068); (W(1-W))^2 = {x * x:.6f}")
print(f"   R = 1 - 30*{x * x:.6f} = {1 - 30 * x * x:.4f};  Omega_s*R = {Om_s:.4f}*{1 - 30 * x * x:.4f} = {Om_s * (1 - 30 * x * x):+.4f}")
print(f"   Reihe bis n = 3 von Hand: cos(2pi W) + cos(4pi W)/16 + cos(6pi W)/81 mit cos(2,4000) = {math.cos(2 * math.pi * Wg):.4f}, "
      f"cos(4,8000) = {math.cos(4 * math.pi * Wg):.4f}, cos(7,2000) = {math.cos(6 * math.pi * Wg):.4f}")
s3 = math.cos(2 * math.pi * Wg) + math.cos(4 * math.pi * Wg) / 16 + math.cos(6 * math.pi * Wg) / 81
print(f"   Sum_{{n<=3}} cos/n^4 = {s3:.4f}; /zeta4 (1,0823) = {s3 / zeta4:.4f}")
