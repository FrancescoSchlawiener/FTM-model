# -*- coding: utf-8 -*-
"""
Ratenschwelle der Rückstellung: Anteil F(a) des Ganzen, dessen Region-Verlust θ̃ unter der Schwelle
θ_½ = (ln 2)²/τ_D liegt (w < ½), und der Gipfel von dF/d ln a.
  w = 1 − e^{−√(θ̃τ_D)} < ½   ⇔   θ̃ < (ln 2)²·a₀(z)/c        (Kern: gewicht_w, a0_frueh_relativ)
Blick: θ̃ ist die Rate eines Bereichs (Region als Ganzes), τ_D die Rückkehrzeit der inneren Richtung (I),
F ist ein Massen- bzw. Volumenanteil des Ganzen (G), dimensionslos, epochenabhängig; keine Dichte am Ort.

1 Schwelle je Epoche: Top-Hat (Kern-Hintergrund) → lineare Dichte δ_½(a), mit der der Bereich bei a die Schwelle kreuzt;
  Vergleich mit Umkehr δ_ta(a) und Kollaps δ_c(a). θ̃ = kleinstes θ der Geschichte (vor der Umkehr fällt θ monoton).
2 Anteil: F(a; M) = erfc(δ_½(a)/(√2·σ(M,a)))  (Massenanteil, Ausflug-Zählung, Setzung) bzw. ½erfc (ohne Zählung von
  Wolke in Wolke); σ(M) aus CAMB mit A_s = 2,1·10⁻⁹, n_s = 0,9649 [G]; M = gezählte Masse (1 + Q)M_b.
3 Gipfel von dF/d ln a gegen z, F(heute) je Glättungsmasse M.
4 Volumenanteil (grob, Kugel-Abbildung).
5 Zweiter Weg: analytisch (Einstein-de-Sitter-Grenzfall, Wachstum als Integral), und von Hand.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
from scipy.special import erfc, erfcinv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, Mpc, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM, H0_KMS,
                      H0_si, E_hintergrund, a0_frueh_relativ, gewicht_w, rho_m_mitbewegt, omega_vakuum, tau_D)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)

H0 = H0_si()
TH_C = (A0 / c) / H0                      # θ_c in H₀
X5 = math.log(2) ** 2
AI = 1e-3


def E(a):
    return E_hintergrund(1 / a - 1)


def th_half(a):                           # θ_½ in H₀
    return X5 * TH_C * a0_frueh_relativ(1 / a - 1)


def rhs(t, y):
    a, x, v = y
    return [a * E(a), v, -OMEGA_M / (2 * x * x) - OMEGA_R * a**-4 * x + OMEGA_DE * x]


def start(di):
    x = AI * (1 - di / 3)
    v = x * E(AI) * (1 - di / (3 * (1 + di)))
    return [AI, x, v]


T_END = 60.0                               # in 1/H₀, weit in der Zukunft

# ------------------------------------------------------------------ lineares Wachstum aus derselben Gleichung
lin = solve_ivp(rhs, (0.0, T_END), start(1e-7), rtol=1e-11, atol=1e-16, dense_output=True, max_step=0.05)
_tl = np.linspace(0, T_END, 200001); _yl = lin.sol(_tl)
LNA = np.log(_yl[0]); DL = ((_yl[0] / _yl[1])**3 - 1) / 1e-7       # δ_lin(a) je Einheit δ_i
def Dlin(a): return float(np.interp(math.log(a), LNA, DL))
D1 = Dlin(1.0)

# ------------------------------------------------------------------ 1 Schwelle je Epoche
def kreuzung(di):
    """a bei θ = θ_½ (w = ½), bei Umkehr (θ = 0) und bei Kollaps."""
    def ev_h(t, y): return 3 * y[2] / y[1] - th_half(y[0])
    ev_h.direction = -1
    def ev_ta(t, y): return y[2]
    ev_ta.direction = -1
    def ev_c(t, y): return y[1] - 1e-3 * AI
    ev_c.direction = -1; ev_c.terminal = True
    s = solve_ivp(rhs, (0.0, T_END), start(di), events=[ev_h, ev_ta, ev_c], rtol=1e-10, atol=1e-15, max_step=0.05)
    g = lambda k: s.y_events[k][0][0] if len(s.y_events[k]) else None
    # Probe Preisach: θ vor der Umkehr monoton fallend → θ̃ = θ bis zur Kreuzung
    return g(0), g(1), g(2)


DIS = np.geomspace(1.5e-4, 0.3, 240)
rows = []
for di in DIS:
    ah, ata, ac = kreuzung(di)
    rows.append((di, ah, ata, ac))

def schwelle(idx):
    A, Dd = [], []
    for r in rows:
        if r[idx] is not None and 0.02 < r[idx] < 20:
            A.append(r[idx]); Dd.append(r[0] * Dlin(r[idx]))    # lineare Dichte bei der Kreuzung
    A = np.array(A); Dd = np.array(Dd); o = np.argsort(A)
    return A[o], Dd[o]

AH, DH = schwelle(1); ATA, DTA = schwelle(2); AC, DC = schwelle(3)
SH, STA, SC = CubicSpline(np.log(AH), DH), CubicSpline(np.log(ATA), DTA), CubicSpline(np.log(AC), DC)
def d_half(a): return float(SH(math.log(a)))
def d_ta(a): return float(STA(math.log(a)))
def d_c(a): return float(SC(math.log(a)))

P("0  a₀ = %.2e, τ_D = %.2f Gyr, θ_c = a₀/c = %.4f H₀; w = ½ bei θ̃ = (ln 2)²θ_c = %.4f H₀ = %.4f·3H₀ (Hintergrund heute 3H₀: 1 − w = %.4f)"
  % (A0, tau_D() / (1e9 * 3.15576e7), TH_C, X5 * TH_C, X5 * TH_C / 3, 1 - gewicht_w(3 * H0)))
P()
P("1  Lineare Dichte (zur Epoche a) des Bereichs, der bei a kreuzt:  w = ½ | Umkehr θ = 0 | Kollaps;  θ_½/3H(a)")
for z in (6, 4, 3, 2, 1.5, 1, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.3):
    a = 1 / (1 + z)
    P("   z = %5.2f:  δ_½ = %.4f   δ_ta = %.4f   δ_c = %.4f    θ_½/3H = %.4f" % (z, d_half(a), d_ta(a), d_c(a), th_half(a) / (3 * E(a))))
P("   (Einstein-de-Sitter: δ_ta = (3/20)(6π)^{2/3} = %.4f, δ_c = (3/20)(12π)^{2/3} = %.4f)" % (0.15 * (6 * math.pi)**(2 / 3), 0.15 * (12 * math.pi)**(2 / 3)))

# ------------------------------------------------------------------ 2 Spektrum
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=0.06, tau=0.0544)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=200)
PD = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False, kmax=200,
                                        var1='delta_tot', var2='delta_tot', zmax=10)
KS = np.exp(np.linspace(math.log(1e-4), math.log(200.0), 4000)); P0 = PD.P(0.0, KS)
rho_m0 = rho_m_mitbewegt()


def sig0(M):
    R = (3 * M * M_sun / (4 * math.pi * rho_m0)) ** (1 / 3) / Mpc
    x = KS * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(np.trapezoid(KS**3 * P0 / (2 * math.pi**2) * W**2, np.log(KS)))


LA = np.linspace(math.log(1 / 11), math.log(1.0), 2001)
AA = np.exp(LA); ZZ = 1 / AA - 1
DHA = np.array([d_half(a) for a in AA]); GA = np.array([Dlin(a) / D1 for a in AA])


def kurve(M, faktor=1.0, dth=None):
    s = sig0(M)
    nu_ = (DHA if dth is None else dth) / (s * GA)
    F = faktor * erfc(nu_ / math.sqrt(2))
    dF = np.gradient(F, LA)
    i = int(np.argmax(dF))
    return s, F, dF, i


def zgipfel(dF, i):
    if i == 0 or i == len(dF) - 1:
        return None
    y0, y1, y2 = dF[i - 1], dF[i], dF[i + 1]
    lg = LA[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (LA[1] - LA[0])
    return math.exp(-lg) - 1


P()
P("2+3  Massenanteil F(a) = erfc(δ_½/√2σ) (Ausflug-Zählung);  σ(M, 0) aus CAMB [G];  Gipfel von dF/d ln a")
P("   M [M☉]     σ(M,0)   F(z=0)   F(z=1)   z_Gipfel   dF/dln a im Gipfel   (½erfc: F(0), z_Gipfel)")
RES = {}
for lm in np.arange(6, 16.01, 0.5):
    M = 10**lm
    s, F, dF, i = kurve(M)
    _, F2, dF2, i2 = kurve(M, faktor=0.5)
    RES[lm] = (s, F[-1], zgipfel(dF, i))
    zg = ("%6.2f" % zgipfel(dF, i)) if 0 < i < len(ZZ) - 1 else (" ≥10 " if i == 0 else " ≤ 0 ")
    zg2 = ("%5.2f" % zgipfel(dF2, i2)) if 0 < i2 < len(ZZ) - 1 else (" ≥10" if i2 == 0 else " ≤ 0")
    P("   10^%-5.1f   %6.3f   %6.3f   %6.3f    %s      %6.3f              (%5.3f, %s)"
      % (lm, s, F[-1], float(np.interp(-math.log(2), LA, F)), zg, dF[i], F2[-1], zg2))

# Glättungsmasse für Gipfel bei z = 0,4 und für F(0) = 2/3
def zpk(lm):
    s, F, dF, i = kurve(10**lm); return zgipfel(dF, i)
lm04 = brentq(lambda lm: zpk(lm) - 0.4, 11.0, 15.5, xtol=1e-3)
s04, F04, dF04, i04 = kurve(10**lm04)
nu23 = math.sqrt(2) * erfcinv(2 / 3)
lm23 = brentq(lambda lm: kurve(10**lm)[1][-1] - 2 / 3, 6, 15.9, xtol=1e-3)
s23, F23, dF23, i23 = kurve(10**lm23)
P()
P("   Gipfel bei z = 0,4 verlangt M = %.2e M☉ (σ₀ = %.3f): dann F(0) = %.3f, F(z=1) = %.3f"
  % (10**lm04, s04, F04[-1], float(np.interp(-math.log(2), LA, F04))))
P("   F(0) = 2/3 verlangt ν₀ = δ_½/σ = %.3f, M = %.2e M☉ (σ₀ = %.3f): dann Gipfel bei z = %.2f" % (nu23, 10**lm23, s23, zgipfel(dF23, i23)))
P("   ½erfc (ohne Wolke-in-Wolke) bleibt unter ½ für jedes M: F = 2/3 unerreichbar.")
# Gipfel mit Kollaps statt Rate-Schwelle (Vergleich)
DCA = np.array([d_c(a) for a in AA]); DTAA = np.array([d_ta(a) for a in AA])
P("   Vergleich bei M = %.2e M☉: Gipfel mit Umkehr-Schwelle z = %.2f, mit Kollaps-Schwelle z = %.2f, mit w = ½ z = %.2f"
  % (10**lm04, zgipfel(*kurve(10**lm04, dth=DTAA)[2:]), (lambda k: zgipfel(k[2], k[3]) if 0 < k[3] < len(ZZ) - 1 else -0.0)(kurve(10**lm04, dth=DCA)), zgipfel(dF04, i04)))
_, FQ, dFQ, iQ = kurve(2.8e14)
P("   Grenzmasse M_Q = 2,8·10¹⁴ M☉ (heute, 03) als Glättung: Gipfel z = %.2f, F(0) = %.3f" % (zgipfel(dFQ, iQ), FQ[-1]))
P("   Form der Rate: dF/dln a bei z = 3 / 2 / 1 / 0,4 / 0 relativ zum Gipfel: %s"
  % " / ".join("%.2f" % (float(np.interp(-math.log(1 + z), LA, dF04)) / dF04[i04]) for z in (3, 2, 1, 0.4, 0)))
P("   Ω_s(heute) = %.3f; Rest der Dunklen Energie 1 − Ω_m − Ω_s = %.3f = %.2f des Ganzen der Dunklen Energie"
  % (omega_vakuum(), OMEGA_DE - omega_vakuum(), (OMEGA_DE - omega_vakuum()) / OMEGA_DE))

# ------------------------------------------------------------------ 4 Volumenanteil (grob)
def volanteil(M, a):
    """Lagrange-Gauß in δ_L, Volumen je Masse 1/(1+δ_NL): Kugel-Abbildung (1 − δ_L/1,5)^{−1,5} unter δ_½,
       darüber zwischen Umkehr (5,55) und Kollaps (178) linear in δ_L (Setzung)."""
    s = sig0(M) * Dlin(a) / D1
    dL = np.linspace(-8 * s, 8 * s, 40001); pg = np.exp(-0.5 * (dL / s)**2)
    dh, dc = d_half(a), d_c(a)
    rho = np.where(dL < dh, (1 - np.minimum(dL, 1.4999) / 1.5)**-1.5,
                   np.where(dL < dc, 5.55 + (178 - 5.55) * (dL - dh) / (dc - dh), 178.0))
    vol = pg / rho
    return np.sum(vol[dL >= dh]) / np.sum(vol), np.sum(pg[dL >= dh]) / np.sum(pg)
P()
P("4  Volumenanteil (grob, Kugel-Abbildung, Gauß ohne Wolke-in-Wolke) bei M = %.2e M☉ und bei 10¹⁰ M☉:" % 10**lm04)
for M in (10**lm04, 1e10):
    P("   M = %.1e:  z = 1 / 0,4 / 0:  Volumen %s   (Masse ½erfc %s)" % (M,
      " / ".join("%.4f" % volanteil(M, 1 / (1 + z))[0] for z in (1, 0.4, 0)),
      " / ".join("%.3f" % volanteil(M, 1 / (1 + z))[1] for z in (1, 0.4, 0))))

# ------------------------------------------------------------------ 5 Zweiter Weg
P()
P("5  Zweiter Weg")
# (a) Grenzfall hoher Rotverschiebung: θ_½/3H → 0, δ_½ → δ_ta (EdS analytisch)
P("   (a) z → groß: θ_½/3H → 0, also δ_½ → δ_ta = (3/20)(6π)^{2/3} = 1,0624 (EdS). Code z = 6: δ_½ = %.4f, δ_ta = %.4f" % (d_half(1 / 7), d_ta(1 / 7)))
# (b) EdS-Zykloide, Schwelle knapp vor der Umkehr: θ/3H = 1 − (3/2)·(η − sinη)·sinη/(1 − cosη)² ...
def eds_delta_lin_bei_theta(r):
    """EdS-Top-Hat: R ∝ 1 − cos η, t ∝ η − sin η; θ/3H = (Ṙ/R)/(ȧ/a) = (3/2)·sinη(η − sinη)/(1 − cosη)²;
       δ_lin = (3/20)(6(η − sinη))^{2/3}. Gibt δ_lin bei θ/3H = r."""
    f = lambda eta: 1.5 * math.sin(eta) * (eta - math.sin(eta)) / (1 - math.cos(eta))**2 - r
    eta = brentq(f, 1e-3, math.pi)
    return 0.15 * (6 * (eta - math.sin(eta)))**(2 / 3)
r0 = th_half(1.0) / (3 * E(1.0))
def eds_code(r):
    rh = lambda t, y: [y[0] * y[0]**-1.5, y[2], -1 / (2 * y[1]**2)]
    def run(di):
        x = AI * (1 - di / 3); v = x * AI**-1.5 * (1 - di / (3 * (1 + di)))
        ev = lambda t, y: 3 * y[2] / y[1] - r * 3 * y[0]**-1.5; ev.terminal = True; ev.direction = -1
        return solve_ivp(rh, (0, 50), [AI, x, v], events=[ev], rtol=1e-11, atol=1e-16)
    di = brentq(lambda d: (run(d).y_events[0][0][0] - 1) if len(run(d).y_events[0]) else 9, 3e-4, 0.05, xtol=1e-13)
    return di / AI
P("   (b0) Ereignis-Methode im EdS-Hintergrund (ODE, Wachstum D = a), θ/3H = %.4f bei a = 1: δ_lin = %.4f" % (r0, eds_code(r0)))
P("   (b) EdS-Zykloide mit θ/3H = %.4f (Verhältnis heute): δ_lin = %.4f (δ_ta 1,0624); Code heute (ΛCDM-Hintergrund) δ_½ = %.4f, δ_ta = %.4f"
  % (r0, eds_delta_lin_bei_theta(r0), d_half(1.0), d_ta(1.0)))
P("       Verhältnis δ_½/δ_ta: EdS %.4f, Code heute %.4f" % (eds_delta_lin_bei_theta(r0) / 1.0624, d_half(1.0) / d_ta(1.0)))
# (c) Wachstum als Integral (Heath), Schwelle fest 1,0624·(δ_½/δ_ta aus EdS), Gipfel analytisch dF/dlna ∝ ν e^{−ν²/2}·f
def D_int(a):
    Ez = lambda x: math.sqrt(OMEGA_M * x**-3 + OMEGA_DE)
    return 2.5 * OMEGA_M * Ez(a) * quad(lambda x: 1 / (x * Ez(x))**3, 1e-8, a)[0]
Di1 = D_int(1.0)
def pk_int(sig, dthr=None):
    la = np.linspace(math.log(1 / 11), 0, 4001)
    nu_ = np.array([(d_half(math.exp(l)) if dthr is None else dthr) / (sig * D_int(math.exp(l)) / Di1) for l in la])
    F = erfc(nu_ / math.sqrt(2)); dF = np.gradient(F, la)
    return 1 / math.exp(la[int(np.argmax(dF))]) - 1, F[-1]
zi, Fi = pk_int(s04); zf, Ff = pk_int(s04, 1.0624)
P("   (c) Wachstum als Integral D ∝ E∫da/(aE)³ (ohne Strahlung) statt ODE, σ₀ = %.3f, Schwelle δ_½(a): Gipfel z = %.2f, F(0) = %.3f (Code 0,40, %.3f);"
  % (s04, zi, Fi, F04[-1]))
P("       mit fester Schwelle 1,0624 (ohne Anstieg durch den Vakuumteil): Gipfel z = %.2f, F(0) = %.3f" % (zf, Ff))
P("   (d) Von Hand, EdS mit fester Schwelle: F = erfc(ν/√2), ν = δ/(σ₀a) ⇒ dF/dln a = √(2/π)·ν·e^{−ν²/2}, Gipfel bei ν = 1.")
P("       F(0) = 2/3 ⇔ ν₀ = √2·erfc⁻¹(2/3) = %.4f; Gipfel bei ν = 1 ⇔ a_G = ν₀ = 0,431, z_G = %.2f (EdS). Gipfel z = 0,4 ⇔ ν₀ = 1/1,4 = 0,714,"
  % (nu23, 1 / nu23 - 1))
P("       F(0) = erfc(0,714/√2) = %.3f (EdS).  Beides zugleich geht mit einer Schwelle und einer Glättungsmasse nicht." % erfc(1 / 1.4 / math.sqrt(2)))
P("       Code (ΛCDM-Wachstum, Rate-Schwelle): F(0) = 2/3 → z_G = %.2f; z_G = 0,4 → F(0) = %.3f." % (zgipfel(dF23, i23), F04[-1]))
P("   (e) Von Hand θ_½ heute: (ln 2)² = 0,48045; a₀/(cH₀) = 1,1·10⁻¹⁰/(2,9979·10⁸·2,1843·10⁻¹⁸) = 0,16798; θ_½ = 0,08071 H₀ = 0,0269·3H₀")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
