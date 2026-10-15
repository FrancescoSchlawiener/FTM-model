# -*- coding: utf-8 -*-
"""
Rest (18 %) und später Anstieg des inaktiven Fluchtmaßes: fehlt ein zwingender Anteil?
1 Selbstkonsistenz: Das Atmen H_eff ändert die Strukturbildung (Zeit t(a), Wachstum D(a), Umkehrdichte Δ_ta, Schwelle δ_ta,
  Massenanteil F) und damit X(a). Bisher: Bildung im Hintergrund des Standardmodells (E_Λ). Iteration H → Bildung → X → R → H.
2 Rückwirkung auf die Übergabe: R_ta und Naht im Atmen, θ̃ = 0 im gebundenen Bereich.
3 Weitere Quellen inaktiven Maßes: Resonanzzungen (07), Haufen jenseits R₂₀₀, Mode ξ bei der Übergabe.
4 Tabelle je Variante [Messwerte nur Test, G].  5 Beobachterfehler.  6 Von Hand.
Aus ftm_kern: OMEGA_*, Q_KOSM, H0_si, Gyr, M_sun, omega_vakuum, E_hintergrund, q_fenster, d_LSS, KREISPAAR_GRENZE.
Aus de_fluchtmoden_mass / de_fluchtmoden_verzoegert / de_fluchtmoden (Hilfsfunktionen): mass_bereich, R_ta, umkehr_dichte,
delta_ta, σ(M), f(T) (Kern-Rückkehrdichte), T_F.
Lokal: Kugelmodell der Umkehr in beliebigem Hintergrund, Faltung R(t) auf eigenem Zeitgitter, Iteration.
"""
import math, os, sys, io, contextlib, warnings
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import erfc
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
from ftm_kern import (M_sun, Gyr, OMEGA_M, OMEGA_R, OMEGA_DE, Q_KOSM, omega_vakuum, E_hintergrund, q_fenster, d_LSS,
                      KREISPAAR_GRENZE)

# ---- Hilfsfunktionen aus de_fluchtmoden_mass (Kopf; enthält den Kopf von de_fluchtmoden_verzoegert)
_src = open(os.path.join(HIER, "de_fluchtmoden_mass.py"), encoding="utf-8").read().split("\n")
_stop = next(i for i, l in enumerate(_src) if l.startswith("# ============================================================== 1 Herleitung"))
M_ = {"__file__": os.path.join(HIER, "de_fluchtmoden_mass.py"), "__name__": "dfmass"}
with contextlib.redirect_stdout(io.StringIO()):
    exec("\n".join(_src[:_stop]), M_)
V = M_['V']; dfm = M_['dfm']; AA = M_['AA']; mass_bereich = M_['mass_bereich']
E_L = M_['E_L']; EB = M_['EB']; H0 = V['H0']; fdich = V['fdich']; T_F = V['T_F']

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

NK = 1 + Q_KOSM
OM_S = omega_vakuum()
LN10 = math.log(10)
HG0 = H0 * Gyr                                              # H₀ in 1/Gyr

# ---- ln B je (Epoche, Masse) bei R_ta und Steigung d ln B / d ln R (für verschobene Naht)
TB = {}; TS = {}
with contextlib.redirect_stdout(io.StringIO()):
    for a in AA:
        r1, r2 = [], []
        for l in dfm.LMg:
            Mb = 10**l * M_sun / NK; R = dfm.R_ta(Mb, a)
            b1 = mass_bereich(Mb, R)["B"]; b2 = mass_bereich(Mb, 1.05 * R)["B"]
            r1.append(b1); r2.append((b2 - b1) / math.log(1.05))
        TB[a] = np.array(r1); TS[a] = np.array(r2)

# ============================================================== Werkzeuge
LNA = np.linspace(math.log(1e-3), math.log(3.2), 4000)
def zeit(Efun):
    """t(a) in Gyr (Materieanfang t ∝ a^{3/2} unter a = 10⁻³)."""
    Ev = np.array([Efun(math.exp(x)) for x in LNA])
    t = cumulative_trapezoid(1 / Ev, LNA, initial=0.0) + (2 / 3) / math.sqrt(OMEGA_M) * 1e-3**1.5 / math.sqrt(1 + OMEGA_R / (OMEGA_M * 1e-3))
    return CubicSpline(LNA, t / HG0)

def wachstum(Efun, a_end=3.2):
    """D(a) mit D = a bei a = 10⁻³ (gleiche frühe Amplitude); gibt (D, dD/dln a)."""
    def f(x, y):
        a = math.exp(x); h = 1e-5
        dl = (math.log(Efun(a * math.exp(h))) - math.log(Efun(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Efun(a)**2 * y[0]]
    s = solve_ivp(f, (math.log(1e-3), math.log(a_end)), [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-13)
    return lambda a: s.sol(math.log(a))

def kugel(Eth, acc, a_list):
    """Umkehr einer Kugel (Einheiten H₀ = 1): R̈ = −GM/R² + acc(a)·R, Hintergrund Eth(a) (ohne Strahlung).
    Gibt je a_ta: Δ_ta = ρ_Kugel/ρ̄ bei der Umkehr, δ_lin(a_ta)."""
    lna = np.linspace(math.log(1e-3), math.log(3.5), 3000)
    Ev = np.array([Eth(math.exp(x)) for x in lna])
    t = cumulative_trapezoid(1 / Ev, lna, initial=0.0) + (2 / 3) / math.sqrt(OMEGA_M) * 1e-3**1.5
    lna_t = CubicSpline(t, lna)
    ai = 1e-3; ti = t[0]; Hi = Eth(ai)
    def f(x, y):
        a = math.exp(x); h = 1e-5
        dl = (math.log(Eth(a * math.exp(h))) - math.log(Eth(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Eth(a)**2 * y[0]]
    gs = solve_ivp(f, (math.log(ai), math.log(3.5)), [ai, ai], dense_output=True, rtol=1e-10, atol=1e-14)
    def lauf(di):
        GM = 0.5 * OMEGA_M * (1 + di)
        rhs = lambda tt, y: [y[1], -GM / y[0]**2 + acc(math.exp(float(lna_t(tt)))) * y[0]]
        ev = lambda tt, y: y[1]; ev.terminal = True; ev.direction = -1
        s = solve_ivp(rhs, (ti, t[-1]), [ai, Hi * ai * (1 - di / 3)], events=ev, rtol=1e-10, atol=1e-14)
        return (s.t_events[0][0], s.y_events[0][0][0]) if len(s.t_events[0]) else (t[-1] * 2, None)
    out = []
    for a_ta in a_list:
        tz = float(np.interp(math.log(a_ta), lna, t))
        ld = brentq(lambda x: lauf(10**x)[0] - tz, -4.5, 0.3, xtol=1e-11)
        di = 10**ld; Rt = lauf(di)[1]
        out.append(((1 + di) * (a_ta / Rt)**3, di * gs.sol(math.log(a_ta))[0] / ai))
    return np.array(out)

A_K = np.array([1 / 31, 1 / 11, 1 / 7, 1 / 4, 1 / 3, 1 / 2.5, 1 / 2, 1 / 1.5, 0.8, 1.0, 1.4, 2.0, 3.0])
KL = kugel(lambda a: math.sqrt(OMEGA_M * a**-3 + 1 - OMEGA_M), lambda a: 1 - OMEGA_M, A_K)   # Standardmodell-Kugel

# ---- Faltung auf eigenem Zeitgitter
NT = 3000
def uebergang(TA, X):
    """R(t) aus der Quelle S = dX/dt, Übergabe ab I_F (t_F = 1,8·t_ta), θ̃ = 0; TA, X auf AA."""
    TG = np.linspace(TA[0], TA[-1], NT); dt = TG[1] - TG[0]
    Fi = fdich(np.arange(NT) * dt)
    sp = CubicSpline(TA, X)
    tq = TG / T_F
    S = sp(np.clip(tq, TA[0], None), 1) * (tq >= TA[0]) / T_F
    t_puls = TA[0] * T_F
    R = np.array([np.trapezoid(Fi[n::-1] * S[:n + 1], dx=dt) for n in range(NT)])
    R += X[0] * fdich(TG - t_puls) * (TG > t_puls)
    return TG, R, cumulative_trapezoid(R, TG, initial=0.0)

def X_von(Dfun, dfak, rfak):
    """X(a) auf AA (geschachtelt, M ≥ 10⁶): F = erfc(δ_ta·dfak/(√2σD)), ln B an der verschobenen Naht R_ta·rfak."""
    X = []
    for k, a in enumerate(AA):
        Fg = erfc(dfm.delta_ta(a) * dfak[k] / (math.sqrt(2) * dfm.SIG * Dfun(a)))
        Lv = np.interp(dfm.LM, dfm.LMg, TB[a] + TS[a] * math.log(rfak[k]))
        dF = -np.gradient(Fg, dfm.LM * LN10)
        msk = dfm.LM >= 6.0
        X.append(float(np.trapezoid((dF * Lv)[msk], dfm.LM[msk] * LN10)))
    return np.array(X)

DL1 = wachstum(E_L)(1.0)[0]

def schritt(E_str, kugelart):
    """Eine Runde: Bildung im Hintergrund E_str → X(a) → R(t) → neues E_eff(a)."""
    tf = zeit(E_str)
    TA = np.array([float(tf(math.log(a))) for a in AA])
    g = wachstum(E_str)
    Dfun = lambda a: g(a)[0] / DL1
    if kugelart == "Λ":
        dfak = np.ones(len(AA)); rfak = np.ones(len(AA)); KK = None
    else:
        lg = np.linspace(math.log(1e-3), math.log(3.6), 2000)
        E2 = np.array([E_str(math.exp(x))**2 for x in lg])
        OXv = E2 - OMEGA_M * np.exp(-3 * lg) - OMEGA_R * np.exp(-4 * lg) - OM_S
        accv = OM_S + OXv + 0.5 * np.gradient(OXv, lg)
        Ethv = np.sqrt(np.clip(E2 - OMEGA_R * np.exp(-4 * lg), 1e-30, None))
        KK = kugel(lambda a: float(np.interp(math.log(a), lg, Ethv)), lambda a: float(np.interp(math.log(a), lg, accv)), A_K)
        dfak = np.interp(np.log(AA), np.log(A_K), KK[:, 1] / KL[:, 1])
        rfak = np.interp(np.log(AA), np.log(A_K), (KL[:, 0] / KK[:, 0])**(1 / 3))
    X = X_von(Dfun, dfak, rfak)
    TG, R, Xarr = uebergang(TA, X)
    AGn = np.exp(np.interp(TG, TA, np.log(AA)))
    E_neu = np.array([EB(a) for a in AGn]) + R / (3 * HG0)
    return dict(TA=TA, X=X, TG=TG, AG=AGn, R=R, Xarr=Xarr, E=E_neu, t0=float(tf(0.0)), D1=Dfun(1.0),
                dfak=dfak, rfak=rfak, KK=KK)

def E_aus(d):
    AG, E = d['AG'], d['E']
    return lambda a: float(np.interp(a, AG, E)) if a >= AG[0] else EB(a)

def iteriere(kugelart, start=None, n_max=25, tol=2e-5):
    d = start if start is not None else schritt(E_L, kugelart)
    verlauf = [d]
    for it in range(n_max):
        Ealt = E_aus(d)
        dn = schritt(Ealt, kugelart)
        En = E_aus(dn)
        dmax = max(abs(En(a) / Ealt(a) - 1) for a in np.exp(np.linspace(math.log(AA[0]), math.log(2.0), 200)))
        d = dn; verlauf.append(d)
        if dmax < tol:
            break
    return d, verlauf, dmax

def ausw(d):
    Ef = E_aus(d)
    a = np.exp(np.linspace(math.log(1 / 4.2), math.log(1.5), 800))
    E = np.array([Ef(x) for x in a])
    OX = E**2 - OMEGA_M * a**-3 - OMEGA_R * a**-4                 # Vakuumteil + Übergang
    lna = np.log(a)
    wz = -1 - np.gradient(np.log(OX), lna) / 3
    m2 = (a >= 1 / 3) & (a <= 1.0)
    A_ = np.vstack([np.ones(m2.sum()), 1 - a[m2]]).T
    w0, wa = np.linalg.lstsq(A_, wz[m2], rcond=None)[0]
    at = lambda arr, z: float(np.interp(1 / (1 + z), a, arr))
    return dict(OX0=at(OX, 0) - OM_S, OD={z: at(OX, z) / at(E, z)**2 for z in (0, 0.5, 1, 2)}, w=at(wz, 0), w0=w0, wa=wa,
                hz={z: Ef(1 / (1 + z)) / E_L(1 / (1 + z)) for z in (0, 0.5, 1, 2, 3)}, E=Ef)

def D_M(Efun, z=1090.0):
    zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, z, 4001)[1:]])
    return np.trapezoid(1 / np.array([Efun(1 / (1 + x)) for x in zs]), zs)
dm_L = D_M(E_L)
gL = wachstum(E_L)
def gesamt(d):
    e = ausw(d); Ef = e['E']
    dm = D_M(Ef)
    gE = wachstum(Ef)
    e['DM'] = dm / dm_L
    e['Q'] = q_fenster(L_min=KREISPAAR_GRENZE * d_LSS() * dm / dm_L)[1]
    e['fs'] = {z: gE(1 / (1 + z))[1] / gL(1 / (1 + z))[1] for z in (0, 0.3, 0.6, 1.0, 1.5)}
    e['X0'] = float(np.interp(1.0, AA, d['X'])); e['Xan0'] = float(np.interp(d['t0'], d['TG'], d['Xarr']))
    e['t0'] = d['t0']; e['D1'] = d['D1']
    return e

# ============================================================== 1 Selbstkonsistenz
kopf("1  Selbstkonsistenz: Bildung im eigenen Atmen statt im Standardhintergrund")
KL_ch = [dfm.umkehr_dichte(a) for a in A_K[[1, 4, 6, 9]]]
P("   Kugelmodell geprüft (Standardhintergrund): Δ_ta lokal %s gegen de_fluchtmoden %s"
  % (", ".join("%.3f" % x for x in KL[[1, 4, 6, 9], 0]), ", ".join("%.3f" % x for x in KL_ch)))
P("   δ_ta lokal (ohne Strahlung) z = %s: %s; Tafel de_uebergabe_anteil: %s (verwendet wird nur das Verhältnis)"
  % (", ".join("%.2f" % (1 / a - 1) for a in A_K[[1, 4, 6, 9]]), ", ".join("%.4f" % x for x in KL[[1, 4, 6, 9], 1]),
     ", ".join("%.4f" % dfm.delta_ta(a) for a in A_K[[1, 4, 6, 9]])))
P("   Bisher (V0): t(a), D(a), σ·D, δ_ta, Δ_ta aus E_Λ; das Atmen E_eff = E_b + R/3H₀ nur als Ergebnis.")
P("   S1: t(a), D(a), σ·D aus E_eff (iteriert); δ_ta, Δ_ta wie bisher.")
P("   S2: zusätzlich Umkehr der Kugel im eigenen Hintergrund: Die Kugel atmet bis zur Umkehr mit (θ̃ ≠ 0, noch nicht gebunden),")
P("       R̈/R = −GM/R³ + H₀²(Ω_s + Ω_X + ½ dΩ_X/dln a) = dieselbe Beschleunigung wie die Schleifen des Ganzen; danach θ̃ = 0.")
P("       Diese Form ist die, aus der auch die Wachstumsgleichung (Atmen nur über H) folgt: erzwungen.")
RES = {}
d0 = schritt(E_L, "Λ"); RES["V0"] = (d0, 0, 0.0)
for nm, art in (("S1", "Λ"), ("S2", "eigen")):
    d, verl, dm_ = iteriere(art, RES["S1"][0] if nm == "S2" else d0)
    RES[nm] = (d, len(verl) - 1, dm_)
    P("   %s: %d Runden, letzte Änderung max|ΔE/E| = %.1e;  Ω_X(0) je Runde: %s"
      % (nm, len(verl) - 1, dm_, " → ".join("%.4f" % ausw(x)['OX0'] for x in verl[:8]) + (" …" if len(verl) > 8 else "")))
ERG = {k: gesamt(v[0]) for k, v in RES.items()}
P()
P("   Variante  t₀[Gyr]  D(1)/D_Λ(1)  X(0)   angek.(0)  δ_ta·(1)  Δ_ta·(1)   Ω_X(0)  Ω_s+Ω_X  /0,687")
for k in ("V0", "S1", "S2"):
    d = RES[k][0]; e = ERG[k]
    P("   %-8s  %6.2f   %7.4f     %6.3f  %7.4f    %7.4f   %7.4f    %7.4f  %7.4f   %5.3f"
      % (k, e['t0'], e['D1'], e['X0'], e['Xan0'], float(np.interp(0.0, np.log(AA), d['dfak'])),
         float(np.interp(0.0, np.log(AA), d['rfak']))**-3, e['OX0'], OM_S + e['OX0'], (OM_S + e['OX0']) / OMEGA_DE))
P("   δ_ta·, Δ_ta·: Faktor gegen die Werte im Standardhintergrund (heute).")
P("   X(a) je Variante (z = 6, 3, 2, 1, 0,5, 0):")
for k in ("V0", "S1", "S2"):
    d = RES[k][0]
    P("   %-4s " % k + "  ".join("%.3f" % float(np.interp(math.log(1 / (1 + z)), np.log(AA), d['X'])) for z in (6, 3, 2, 1, 0.5, 0)))
P("   Richtung: E_eff < E_Λ (z ≲ 3) ⇒ mehr Zeit und mehr Wachstum ⇒ mehr und frühere Bindung; ältere Welt ⇒ mehr Transit")
P("   angekommen (f(T) hängt nur an τ_D). Umkehr im eigenen Atmen (S2): δ_ta −1,5 %, Δ_ta −9 %, X(0) +0,03 – klein.")

# ============================================================== 2 Übergabe
kopf("2  Rückwirkung auf die Übergabe (Naht, Umkehr)")
d = RES["S2"][0]
P("   Faktor der Naht R_ta (S2 gegen Standard) bei z = 3, 1, 0,5, 0: %s"
  % ", ".join("%.4f" % float(np.interp(math.log(1 / (1 + z)), np.log(AA), d['rfak'])) for z in (3, 1, 0.5, 0)))
P("   Faktor der Schwelle δ_ta: %s" % ", ".join("%.4f" % float(np.interp(math.log(1 / (1 + z)), np.log(AA), d['dfak'])) for z in (3, 1, 0.5, 0)))
P("   Im eigenen Atmen fehlt der Kugel ein Teil der Beschleunigung des Standardmodells (Ω_s + Ω_X < 0,687, später Anstieg):")
P("   sie kehrt bei kleinerer linearer Dichte um (größere Bereiche binden früher, nicht später) und etwas weniger dicht.")
P("   ln B verschiebt sich um (d ln B/d ln R)·ln(R-Faktor); d ln B/d ln R = %.2f (10¹² M☉, heute, tief 3/2)."
  % float(np.interp(12.0, dfm.LMg, TS[AA[int(np.argmin(abs(AA - 1)))]])))
P("   θ̃ = 0 nach der Umkehr: R_ta physisch fest, das Atmen dehnt nur noch den Fluss; kein weiterer Anteil.")

# ============================================================== 3 Weitere Quellen
kopf("3  Weitere Quellen inaktiven Maßes (nur was aus dem Modell zwingend folgt)")
P("   (a) Resonanzzungen (07): Einrasten ändert die Gestalt der Bahn (T² → T¹, Insel um die periodische Bahn), nicht die Menge")
P("       der Zustände, deren ganzer Torus in der Naht liegt; das Maß B zählt Inseln schon mit (Liouville: Inseln haben das")
P("       Maß ihrer Zustände). Einrasten zu echter Verriegelung verlangt Dissipation, die nicht 𝒫 ist. Gegenüber Newton hebt die")
P("       Regel die Kepler-Verriegelung aller Bahnen auf (Vorzeichen: frei, nicht inaktiv). Beitrag 0; eigener Ansatz wäre Doppelzählung.")
LMm = dfm.LM; a1 = AA[int(np.argmin(abs(AA - 1)))]
Fg = erfc(dfm.delta_ta(1.0) * float(np.interp(0.0, np.log(AA), d['dfak'])) / (math.sqrt(2) * dfm.SIG * d['D1']))
Lv = np.interp(LMm, dfm.LMg, TB[a1]); dF = -np.gradient(Fg, LMm * LN10)
tot = np.trapezoid(dF * Lv, LMm * LN10)
for lmc in (13.0, 14.0):
    m = LMm >= lmc
    P("   (b) Gruppen/Haufen: Anteil M ≥ 10^%.0f M☉ an X(0) = %.3f (%.1f %%); gezählt bis zur Naht R_ta ≈ 3–4·R₂₀₀ (Δ_ta ≈ 5,6 gegen 200)."
      % (lmc, np.trapezoid((dF * Lv)[m], LMm[m] * LN10), 100 * np.trapezoid((dF * Lv)[m], LMm[m] * LN10) / tot))
P("       Jenseits R_ta ist Fluss (θ̃ ≠ 0), dort ist nichts gebunden; die Schale bis 1,16·R_ta (w > 0) gab +0,016 in Ω_X")
P("       (de_fluchtmoden_mass C), nur Spannweite. Kein fehlender Anteil.")
P("   (c) Mode ξ bei der Übergabe: Kac 1 + Q zählt die Rückkehr in I und ist als Ω_c gebucht; bei I_F wechselt die Mode die")
P("       Form (w: 1 → 0, Regel ν), ihr Maß verschwindet nicht und wird nicht frei. Sie geht nur über M_ges = (1 + Q)M_b in R_ta")
P("       ein (schon gezählt). Als zusätzliches inaktives Maß gebucht: Doppelzählung (Faktor 6,36 auf X, Mode als Stoff). Beitrag 0.")

# ============================================================== 4 Tabelle
kopf("4  Tabelle je Variante [Messwerte nur Test, G: Ω_DE = 0,687; Standardmodell Ω_DE(z) = 0,394 / 0,215 / 0,068 bei z = 0,5 / 1 / 2]")
P("   Var  Ω_X(0)  Ω_DE-Anteil z=0,5  z=1    z=2     w(0)    w₀      w_a     H/H_Λ z=0    0,5    1      2      3     D_M(1090)  Q oben  fσ₈/Λ z=0   0,3    0,6    1,0    1,5")
for k in ("V0", "S1", "S2"):
    e = ERG[k]
    P("   %-3s  %6.3f  %6.3f          %6.3f %6.3f   %6.3f  %6.3f  %6.3f    %6.3f  %6.3f %6.3f %6.3f %6.3f   %6.3f     %5.2f   %6.3f  %6.3f %6.3f %6.3f %6.3f"
      % (k, e['OX0'], e['OD'][0.5], e['OD'][1], e['OD'][2], e['w'], e['w0'], e['wa'], e['hz'][0], e['hz'][0.5], e['hz'][1],
         e['hz'][2], e['hz'][3], e['DM'], e['Q'], e['fs'][0], e['fs'][0.3], e['fs'][0.6], e['fs'][1.0], e['fs'][1.5]))
P("   Q-Fenster unten 5 (ΛCDM oben %.2f). Ω_DE-Anteil = (Ω_s + Ω_X)(z)/E²(z); fσ₈ mit gleicher früher Amplitude." % q_fenster()[1])
eS = ERG["S2"]; e0 = ERG["V0"]
P("   Selbstkonsistent (S2) gegen bisher: Ω_X(0) %+.3f (%+.1f %%), Ω_DE-Anteil z = 1 %+.3f, fσ₈(0) %+.3f."
  % (eS['OX0'] - e0['OX0'], 100 * (eS['OX0'] / e0['OX0'] - 1), eS['OD'][1] - e0['OD'][1], eS['fs'][0] - e0['fs'][0]))
P("   Lücke: Ω_s + Ω_X = %.3f von 0,687 (%.0f %%); fehlend %.3f." % (OM_S + eS['OX0'], 100 * (OM_S + eS['OX0']) / OMEGA_DE, OMEGA_DE - OM_S - eS['OX0']))

# ============================================================== 5 Blick
kopf("5  Beobachterfehler")
P("   X: Maß der Tori je Bereich (G des Bereichs), je Baryon über die Zelle (G). Weg: Rückkehr in I (τ_D = 2Λ/c fest).")
P("   Atmen: Rate des Ganzen (G); es wirkt auf Bildung nur über H und über die Beschleunigung der noch mitatmenden Schleifen")
P("   bis zur Umkehr (θ̃ ≠ 0); im gebundenen Bereich θ̃ = 0. Nirgends Dichte am Ort, kein Stoff, der ballt oder drückt.")
P("   Zustandsabhängig: X, t₀, D, δ_ta, Δ_ta, Ω_X(z) – selbst Folge des Atmens; fest: Λ, a₀, τ_D, Ω_s.")

# ============================================================== 6 Von Hand
kopf("6  Von Hand")
Eeff1 = eS['E'](1.0); EL1 = 1.0
Om_e = OMEGA_M / Eeff1**2
P("   Wachstumsrate heute nach f ≈ Ω_m(a)^0,55 (anderer Weg als die Wachstumsgleichung):")
P("   E_eff(1) = %.4f ⇒ Ω_m,eff = 0,3134/%.4f² = %.4f; f_eff/f_Λ = (%.4f/0,3134)^0,55 = %.4f" % (Eeff1, Eeff1, Om_e, Om_e, (Om_e / OMEGA_M)**0.55))
gE = wachstum(eS['E'])
lnD = lambda Ef: np.trapezoid([(OMEGA_M * x**-3 / Ef(x)**2)**0.55 for x in np.exp(np.linspace(math.log(1e-3), 0, 3000))], np.linspace(math.log(1e-3), 0, 3000))
rD = math.exp(lnD(eS['E']) - lnD(E_L))
P("   D_eff/D_Λ(1) = exp∫(f_eff − f_Λ)dln a = %.4f (Wachstumsgleichung %.4f); fσ₈(0)-Verhältnis ≈ %.4f·%.4f = %.3f (Code %.3f)"
  % (rD, gE(1.0)[0] / gL(1.0)[0], (Om_e / OMEGA_M)**0.55, rD, (Om_e / OMEGA_M)**0.55 * rD, eS['fs'][0]))
P("   Ω_X(0) = E_eff(1)² − E_b(1)² = %.4f² − 0,7467² = %.4f (Code %.4f); E_b(1) = √(0,3134 + 0,2440 + 0,0001)."
  % (Eeff1, Eeff1**2 - 0.7467**2, eS['OX0']))

with open(os.path.join(HIER, "de_rest_anstieg.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
