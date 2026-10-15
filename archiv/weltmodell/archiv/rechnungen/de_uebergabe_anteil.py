# -*- coding: utf-8 -*-
"""
Übergebene Mode und Staubanteil des Ganzen.
Frage: Zählt die Mode ξ nach der Übergabe (erste Faltung, I_F = ∫ω dt ≈ 1,1–1,25 rad nach der Umkehr) im Atmen des
Ganzen (G) weiter als Staub (1 + Q)? Wenn nicht: übergebener Massenanteil F(a), H(z), und was ein Beobachter mit
festem Ω_m als Dunkle Energie liest.
Blick: F ist ein Massenanteil des Ganzen (G, epochenabhängig); ν am Ort ist L; ν_kosm = 1 + Q ist Kac-Zählung (I/G).
1 Zählung aus dem Modell: Kac, Bilanz an R_Q, gleiche Zählung in Λ und im Staub; Probe über a₀(z).
2 F(a) = erfc(δ_F(a)/√2σ(M_min, a)), δ_F aus dem Top-Hat im Kern-Hintergrund (Ereignis I = I_F nach der Umkehr);
  σ(M) aus CAMB [G]; Abhängigkeit von der kleinsten gebundenen Masse M_min.
3 Gesetzt, sie zählte nicht: Lücke der Dunklen Energie, gelesenes w₀, w_a, Ω_DE(z) bei festem Ω_m.
4 Von Hand.
"""
import math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import erfc
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, G, c, M_sun, Mpc, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B, OMEGA_B_H2, OMEGA_C_H2, Q_KOSM,
                      H0_KMS, E_hintergrund, rho_m_mitbewegt, omega_vakuum, nu, r_Q, a0_frueh_relativ)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)

F_Q = Q_KOSM / (1 + Q_KOSM)                  # Anteil der Mode am Staub, 0,843
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
IF_WERTE = (1.10, 1.25)
AI = 1e-3
T_END = 40.0

# ============================================================== 1 Zählung
P("=" * 110)
P("1  Zählt die übergebene Mode im Atmen des Ganzen als Staub?")
P("=" * 110)
q_k = Q_KOSM / (1 + Q_KOSM)
P("   Kac: ν_kosm = 1/(1 − q) = 1/μ(offener Bereich) = %.4f = 1 + Q; das Maß hängt nicht daran, wie die Rückläufe" % (1 / (1 - q_k)))
P("   örtlich verteilt sind (Kac mittelt über den ganzen Zustandsraum). Übergabe ändert die Form (Mode → Operator), nicht das Maß.")
Mb = 6e10 * M_sun
RQ = r_Q(Mb)
aN = G * Mb / RQ**2
P("   Bilanz im gebundenen Bereich (Punkt 6·10¹⁰ M☉): an R_Q = %.1f kpc ist ν = %.4f = 1 + Q (%.4f); innen ν > 1 + Q, außen < 1 + Q,"
  % (RQ / (Mpc / 1e3), nu(aN), 1 + Q_KOSM))
P("   jenseits 1,05·R_ta liest man (1 + w·Q)·M_b → (1 + Q)·M_b (03): von außen (Ebene darüber) zählt der Bereich wie vorher.")
P("   Gleiche Zählung in beiden Gliedern der Ausdehnung: Staub (1 + Q)ρ_b und Λ = (1 + Q)·ΣΔs_b (03, Kac in I).")
P("   Probe: fiele der übergebene Anteil aus der Zählung, fiele auch ν_kosm im Anteil F: Λ ∝ 1 − f·F, a₀ ∝ 1/(1 − f·F), f = Q/(1+Q) = %.3f." % F_Q)
P("   (Werte F aus Teil 2; Modell: a₀ seit z ≈ 10 fest auf 0,3 %)")

# ============================================================== 2 Schwelle δ_F(a)
def E(a): return E_hintergrund(1 / a - 1)

def rhs(t, y):
    a, x, v, I = y
    om = math.sqrt(OMEGA_M / (2 * x**3)) if x > 0 else 0.0
    return [a * E(a), v, -OMEGA_M / (2 * x * x) - OMEGA_R * a**-4 * x + OMEGA_DE * x, om if v < 0 else 0.0]

def start(di):
    x = AI * (1 - di / 3)
    v = x * E(AI) * (1 - di / (3 * (1 + di)))
    return [AI, x, v, 0.0]

lin = solve_ivp(lambda t, y: rhs(t, y)[:3] + [0.0], (0, T_END), start(1e-7), rtol=1e-11, atol=1e-16, dense_output=True, max_step=0.05)
_t = np.linspace(0, T_END, 100001); _y = lin.sol(_t)
LNA = np.log(_y[0]); DL = ((_y[0] / _y[1])**3 - 1) / 1e-7
def Dlin(a): return float(np.interp(math.log(a), LNA, DL))
D1 = Dlin(1.0)

def kreuzung(di):
    evs = []
    for IF in IF_WERTE:
        def ev(t, y, IF=IF): return y[3] - IF
        ev.direction = 1
        evs.append(ev)
    def ev_ta(t, y): return y[2]
    ev_ta.direction = -1
    def ev_end(t, y): return y[3] - 1.6
    ev_end.terminal = True
    s = solve_ivp(rhs, (0, T_END), start(di), events=evs + [ev_ta, ev_end], rtol=1e-10, atol=1e-15, max_step=0.02)
    g = lambda k: s.y_events[k][0][0] if len(s.y_events[k]) else None
    return [g(0), g(1), g(2)]

DIS = np.geomspace(2e-4, 0.3, 160)
ROWS = [(di, *kreuzung(di)) for di in DIS]

def schwelle(k):
    A, Dd = [], []
    for r in ROWS:
        if r[k] is not None and 0.02 < r[k] < 12:
            A.append(r[k]); Dd.append(r[0] * Dlin(r[k]))
    A = np.array(A); Dd = np.array(Dd); o = np.argsort(A)
    return CubicSpline(np.log(A[o]), Dd[o])

SF = [schwelle(1), schwelle(2)]; STA = schwelle(3)
def d_F(a, k=0): return float(SF[k](math.log(a)))

P()
P("=" * 110)
P("2  Übergebener Massenanteil F(a) = erfc(δ_F/√2σ(M_min, a))  (Ausflug-Zählung, Masse in Bereichen mit I > I_F)")
P("=" * 110)
P("   Schwelle (lineare Dichte zur Epoche a des Bereichs, der bei a die erste Faltung erreicht), Kern-Hintergrund:")
for z in (10, 6, 3, 2, 1, 0.5, 0.4, 0):
    a = 1 / (1 + z)
    P("   z = %5.2f:  δ_ta = %.4f   δ_F(I_F = 1,10) = %.4f   δ_F(I_F = 1,25) = %.4f" % (z, float(STA(math.log(a))), d_F(a, 0), d_F(a, 1)))
th = lambda IF: 4 * math.atan(math.exp(IF / math.sqrt(2)))
dEdS = lambda IF: 0.15 * (6 * (th(IF) - math.sin(th(IF))))**(2 / 3)
P("   EdS (Zykloide, I = √2·ln tan(θ/4)): δ_F = %.4f / %.4f, t_F/t_ta = %.3f / %.3f"
  % (dEdS(1.10), dEdS(1.25), (th(1.10) - math.sin(th(1.10))) / math.pi, (th(1.25) - math.sin(th(1.25))) / math.pi))

import camb
p = camb.CAMBparams()
p.set_cosmology(H0=H0_KMS, ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, mnu=0.06, tau=0.0544)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=400)
PD = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False, kmax=400,
                                        var1='delta_tot', var2='delta_tot', zmax=10)
KS = np.exp(np.linspace(math.log(1e-4), math.log(400.0), 5000)); P0 = PD.P(0.0, KS)
rho0 = rho_m_mitbewegt()
def sig0(M):
    R = (3 * M * M_sun / (4 * math.pi * rho0))**(1 / 3) / Mpc
    x = KS * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    return math.sqrt(np.trapezoid(KS**3 * P0 / (2 * math.pi**2) * W**2, np.log(KS)))

LA = np.linspace(math.log(1 / 21), 0.0, 1501); AA = np.exp(LA)
GA = np.array([Dlin(a) / D1 for a in AA])
DFA = [np.array([d_F(a, k) for a in AA]) for k in (0, 1)]

def Fkurve(M, k=0):
    s = sig0(M)
    return s, erfc(DFA[k] / (math.sqrt(2) * s * GA))

def Fz(F, z): return float(np.interp(-math.log(1 + z), LA, F))

P()
P("   M_min [M☉]  σ(M,0)   F(z=20)  F(10)   F(3)    F(1)    F(0,4)  F(0)    a₀(0)/a₀(10) falls nicht gezählt   (I_F = 1,25: F(0))")
FMAP = {}
for lm in (6, 7, 8, 9, 10, 11, 12, 13, 14):
    s, F = Fkurve(10.0**lm)
    _, F2 = Fkurve(10.0**lm, 1)
    FMAP[lm] = F
    P("   10^%-3d      %6.3f   %6.3f   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f   %6.3f                             (%6.3f)"
      % (lm, s, Fz(F, 20), Fz(F, 10), Fz(F, 3), Fz(F, 1), Fz(F, 0.4), F[-1], (1 - F_Q * Fz(F, 10)) / (1 - F_Q * F[-1]), F2[-1]))
P("   Hierarchie: jede Masse ab der kleinsten gebundenen zählt; F = Anteil über M_min, mit M_min → klein F → 1 in jeder Epoche.")
P("   Das Modell setzt M_min nicht (Leiter endet unten erst an der Quantenlänge); die Zählung ist nicht eindeutig, F(0) = 0,1–0,9.")

# ============================================================== 3 Gesetzt, sie zählte nicht
P()
P("=" * 110)
P("3  Gesetzt, die übergebene Mode zählte nicht mehr als Staub (Λ fest): Bilanz heute und gelesene Dunkle Energie")
P("=" * 110)
P("   Staub wirksam: Ω_m(1 − f·F(a))·a⁻³. Ω_s = %.4f (Λ), gemessen H₀: nötiger Rest Ω_x = 1 − Ω_r − Ω_m(1 − fF₀) − Ω_s" % OM_S)
P("   Beobachter mit festem Ω_m = %.4f (Hintergrundstrahlung, F = 0) liest ρ_DE(a) = Ω_s + Ω_x − Ω_m·f·F(a)·a⁻³" % OMEGA_M)

def w0_wa(rho):
    """w(a) = −1 − (1/3)·d ln ρ/d ln a, w₀ = w(1), w_a = −dw/da bei a = 1; z, wo ρ_DE gelesen null wird."""
    lr = np.log(np.abs(rho)); w = -1 - np.gradient(lr, LA) / 3
    wa = -(w[-1] - w[-6]) / (AA[-1] - AA[-6])
    neg = np.where(rho <= 0)[0]
    z0 = (1 / AA[neg[-1]] - 1) if len(neg) else None
    return w[-1], wa, z0

P("   M_min     F(0)    Ω_x nötig (statt %.3f)   H₀ ohne Ω_x / gemessen   ρ_DE gelesen z=0/0,5/1/2 (rel. heute)   w₀, w_a   ρ_DE gelesen < 0 ab z" % OM_X)
for lm in (6, 8, 10, 12, 13, 14):
    F = FMAP[lm]
    ox = 1 - OMEGA_R - OMEGA_M * (1 - F_Q * F[-1]) - OM_S
    rho = OM_S + ox - OMEGA_M * F_Q * F * AA**-3
    w0, wa, z0 = w0_wa(rho)
    h_ohne = math.sqrt(OMEGA_M * (1 - F_Q * F[-1]) + OM_S + OMEGA_R)
    rr = [float(np.interp(-math.log(1 + z), LA, rho)) / rho[-1] for z in (0, 0.5, 1, 2)]
    P("   10^%-3d   %6.3f   %6.3f                   %6.3f                  %s        %+.3f, %+.3f   %s"
      % (lm, F[-1], ox, h_ohne, " / ".join("%.3f" % x for x in rr), w0, wa, ("%.2f" % z0) if z0 is not None else "nie (z ≤ 20)"))
P("   Ohne Übergabe (F = 0): Ω_x = %.3f, H₀ ohne Ω_x = %.3f·gemessen." % (OM_X, math.sqrt(OMEGA_M + OM_S + OMEGA_R)))
P("   Vorzeichen: Staubverlust senkt H; gelesen als Dunkle Energie ist er ein negativer, mit der Bindung wachsender Beitrag.")
P("   Die Lücke wächst um Ω_m·f·F₀; gelesen w < −1 (ρ_DE steigt mit der Zeit), gemessen w₀ > −1, w_a < 0 [G].")
P("   Mit Λ in derselben Zählung (Λ ∝ 1 − fF): Ω_s heute durch a₀ heute fest, früher kleiner (w_s < −1); Lücke heute gleich.")

# ============================================================== 4 Von Hand
P()
P("=" * 110)
P("4  Von Hand")
P("=" * 110)
P("   θ_F = 4·atan(e^{1,10/√2}) = 4·atan(2,1772) = 4·1,1406 = 4,5624; θ − sin θ = 4,5624 + 0,9888 = 5,5512")
P("   δ_F(EdS) = 0,15·(6·5,5512)^{2/3} = 0,15·33,307^{2/3} = 0,15·10,352 = 1,553   (Code EdS %.4f)" % dEdS(1.10))
s12 = sig0(1e12)
P("   M_min = 10¹²: σ₀ = 2,245 (de_ratenschwelle), δ_F(0) = 1,571 (Code %.3f): F₀ = erfc(1,571/(1,4142·2,245)) = erfc(0,4948) = 1 − (0,5205 − 0,0052·0,8788) = 0,484 (Code %.3f, σ₀ hier %.3f)"
  % (d_F(1.0), FMAP[12][-1], s12))
P("   Lücke: 0,442 + 0,3134·0,843·0,484 = 0,442 + 0,128 = 0,570;  a₀ wüchse seit z = 10 um 1/(1 − 0,843·0,484) = 1/0,592 = 1,69 (F(10) ≈ 0)")

open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
