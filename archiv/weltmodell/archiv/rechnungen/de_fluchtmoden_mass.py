# -*- coding: utf-8 -*-
"""
Welches Maß X je Baryon verlangt das Modell für das inaktive Fluchtmaß? (Fortsetzung von de_fluchtmoden,
de_fluchtmoden_hierarchie, de_fluchtmoden_verzoegert, de_fluchtmoden_erneuerung)
Kandidaten (vor dem Ergebnis aus den Axiomen geprüft, Abschnitt 1):
  A   Ort × Geschwindigkeit bis R_ta, Newton E > 0 entkommt (bisher)
  A'  Ort × Geschwindigkeit bis R_ta, Newton bis zur selben Naht
  B   Wirkungs-Winkel: ganze Tori innerhalb der Naht R_ta (Regel und Newton)
  C   wie B bis zur w-Schale 1,16·R_ta (Regel dort fortgesetzt)
  D   wie B, Bezug Mitnahme (1 + Q)·Newton statt Newton
  E   linear statt ln (mit B)
2 Übergang erzwungen (M_min = 10⁶, geschachtelt, θ̃ = 0, erste Rückkehr f) ab I_F und ab Umkehr: Ω_X, w, H/H_ΛCDM [G].
3 Gesamtbild: Ω_s + Ω_X gegen Ω_DE [G], D_M(1090), Q-Fenster, fσ₈.  4 Beobachterfehler.  5 Von Hand (Wirkungsweg).
Aus ftm_kern: nu, r_stern, v_flach, G, M_sun, kpc, Q_KOSM, OMEGA_DE, omega_vakuum, E_hintergrund, q_fenster, d_LSS.
Aus de_fluchtmoden / de_fluchtmoden_verzoegert (Hilfsfunktionen): R_ta, Massenfunktion, Rückkehrdichte, lauf, auswertung.
Neu (lokal): Torus-Maß ∫4πr²·(4π/3)v_esc³/(1 − r²/R²) dr (Ellipsoid im Geschwindigkeitsraum), Wirkungsintegral J_r.
"""
import math, os, sys, io, contextlib, warnings
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid, solve_ivp
from scipy.optimize import brentq
from scipy.special import gamma
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
from ftm_kern import (G, M_sun, kpc, Mpc, Gyr, OMEGA_M, OMEGA_DE, Q_KOSM, H0_si, E_hintergrund, omega_vakuum, nu,
                      r_stern, v_flach, q_fenster, d_LSS, KREISPAAR_GRENZE)

# ---- Hilfsfunktionen aus de_fluchtmoden_verzoegert (Kopf bis auswertung, ohne dessen Tabellen)
_src = open(os.path.join(HIER, "de_fluchtmoden_verzoegert.py"), encoding="utf-8").read().split("\n")
V = {"__file__": os.path.join(HIER, "de_fluchtmoden_verzoegert.py"), "__name__": "dfv"}
with contextlib.redirect_stdout(io.StringIO()):
    exec("\n".join(_src[:161]), V)
dfm = V['dfm']; AA = V['AA']; AG = V['AG']; i0 = V['i0']; EB = V['EB']; E_L = V['E_L']; H0 = V['H0']; EBG = V['EBG']

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

NK = 1 + Q_KOSM
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
K6 = 16 * math.pi**2 / 3
DEFS = ("A", "A'", "B", "C", "D", "E")

# ------------------------------------------------ Maße eines Bereichs
def v2_regel(M_b, Rb, n=1500):
    r = np.concatenate([np.geomspace(Rb * 1e-7, 0.9 * Rb, n), Rb * (1 - np.geomspace(0.1, 1e-8, n // 2))[1:], [Rb]])
    aN = G * M_b / r**2
    acc = nu(aN) * aN
    I = cumulative_trapezoid(acc[::-1], -r[::-1], initial=0.0)[::-1]      # ∫_r^Rb a dr'
    return r, np.clip(2 * I, 0, None)

def mass_bereich(M_b, R, n=1500):
    """ln-Verhältnisse je Definition (Punktmasse)."""
    out = {}
    for tag, Rb in (("R", R), ("S", 1.16 * R)):
        r, vR2 = v2_regel(M_b, Rb, n)
        x = r / Rb
        w = K6 * r**2
        wt = np.where(x < 1, w / np.clip(1 - x**2, 1e-300, None), 0.0)
        vN2i = 2 * G * M_b / r
        vN2R = np.clip(2 * G * M_b * (1 / r - 1 / Rb), 0, None)
        I = lambda wt_, v2: float(np.trapezoid(wt_ * v2**1.5, r))
        out[tag] = dict(Ra=I(w, vR2), Ni=I(w, vN2i), NRa=I(w, vN2R), Rb=I(wt, vR2), Nb=I(wt, vN2R))
    o, s = out["R"], out["S"]
    rB = o['Rb'] / o['Nb']
    return {"A": math.log(o['Ra'] / o['Ni']), "A'": math.log(o['Ra'] / o['NRa']), "B": math.log(rB),
            "C": math.log(s['Rb'] / s['Nb']), "D": math.log(rB) - 1.5 * math.log(NK), "E": rB - 1.0,
            "_rat": dict(A=o['Ra'] / o['Ni'], Ap=o['Ra'] / o['NRa'], B=rB, Rb_Ra=o['Rb'] / o['Ra'], Nb_NRa=o['Nb'] / o['NRa'])}

# ============================================================== 1 Herleitung
kopf("1  Welches Maß verlangt das Modell? (Axiome: Liouville, kein Außen, 𝒫 maßerhaltend; vor dem Ergebnis)")
P("   (b) = (a) als Maß: Wirkungs-Winkel ist kanonisch, d³x d³p = d³J d³θ (Jacobi 1); Liouville gibt kein anderes Maß.")
P("       Verschieden ist nur die MENGE. 05: Grundraum aus Tori, Freiheitsgrad = Umlaufrichtung; gebunden an den Bereich ist")
P("       ein Zustand, dessen ganzer Torus in der Naht liegt (w = 0 bis R_ta). (a) zählt am Ort, ob v < v_esc(r → R_ta): das")
P("       schneidet ganze Tori ab (fast kreisförmige Bahnen nahe R_ta mit E > Φ(R_ta), die R_ta nie erreichen).")
P("       Torus in R_ta ⇔ v_r² + v_t²(1 − r²/R²) ≤ v_esc,R²(r): Ellipsoid, Volumen (4π/3)v_esc³/(1 − r²/R²).")
P("   Kein Außen: Newton 'E > 0 entkommt' zählt bis ∞; im Modell gibt es kein ∞, der Bereich endet an der Naht. Beide")
P("       Dynamiken sind am selben Rand zu messen: Newton-Tori innerhalb R_ta (A' und B). A zählt Newton-Bahnen, die den")
P("       Bereich verlassen, als gebunden – gemischter Rand.")
P("   (c) Schale bis 1,16·R_ta: dort ist w > 0 (Mode mitatmend, Übergang zum Fluss), die Form ist gesetzt; gebunden heißt")
P("       θ̃ = 0, also bis R_ta. Nur Spannweite.")
P("   (d) Mitnahme (1 + wQ): innerhalb der Naht ist w = 0 ⇒ 1 + wQ = 1, die Mitnahme des Bereichs IST Newton; (1 + Q)")
P("       gilt nur mitatmend (w = 1, außen). (1 + Q)·M_b als bindende Masse innen läse die Mode als Stoff (Beobachterfehler).")
P("       𝒫_q (05): q → 0 gibt I = Newton (offen); inaktiv heißt gegenüber I. Nur Spannweite.")
P("   (e) ln: N Baryonen, Maß des Ganzen = Produkt der Einteilchenmaße (Liouville auf dem Produktraum); das Atmen δ ln V je")
P("       Baryon ist additiv in ln μ. Linear (μ_R/μ_N − 1) addierte Maße verschiedener Bereiche, als wären sie disjunkte")
P("       Teile eines Einteilchenraums – dann fehlt die Norm (de_fluchtmoden A). Nur Spannweite.")
P("   ⇒ Erzwungen: B (ganze Tori in der Naht, Regel gegen Newton, ln).")

# ---- Tabellen
TABM = {}
with contextlib.redirect_stdout(io.StringIO()):
    for a in AA:
        TABM[a] = [mass_bereich(10**l * M_sun / dfm.NK, dfm.R_ta(10**l * M_sun / dfm.NK, a)) for l in dfm.LMg]

P()
P("   Je Bereich (Punktmasse), Verhältnisse μ_R/μ_N:")
P("   M_b [M☉]   a     R_ta/r*   A        A'       B        ln A    ln A'   ln B    ln C    ln D    B: Regel μ_b/μ_a  Newton μ_b/μ_a")
for lm in (8, 10, 10.778, 12, 14):
    for a in (1.0, 0.5, 0.25):
        Mb = 10**lm * M_sun; R = dfm.R_ta(Mb, a)
        m = mass_bereich(Mb, R); rr = m['_rat']
        P("   10^%-6.3f  %.2f  %7.1f   %7.2f  %7.2f  %7.2f   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f      %.3f             %.3f"
          % (lm, a, R / r_stern(Mb), rr['A'], rr['Ap'], rr['B'], m['A'], m["A'"], m['B'], m['C'], m['D'], rr['Rb_Ra'], rr['Nb_NRa']))
P("   Probe A gegen de_fluchtmoden (Milchstraße a = 1): ln A = %.3f (dort 4.259)." % mass_bereich(6e10 * M_sun, dfm.R_ta(6e10 * M_sun, 1.0))['A'])

# ============================================================== 2 Übergang je Definition
kopf("2  X je Definition mit dem erzwungenen Übergang (M_min = 10⁶, geschachtelt, θ̃ = 0, erste Rückkehr f) [Messwerte G]")
from scipy.special import erfc
LN10 = math.log(10)
SEL = {"d": "A"}
def X_def(a, lm_min, art):
    Fg = erfc(dfm.delta_ta(a) / (math.sqrt(2) * dfm.SIG * dfm.Dg(a)))
    Lv = np.interp(dfm.LM, dfm.LMg, np.array([t[SEL['d']] for t in TABM[a]]))
    msk = dfm.LM >= lm_min
    dF = -np.gradient(Fg, dfm.LM * LN10)
    return float(np.trapezoid((dF * Lv)[msk], dfm.LM[msk] * LN10))
V['X_zaehlung'] = X_def

def ox_bei(d, z):
    return float(np.interp(1 / (1 + z), AG, d['OX']))
ERG = {}
P("   Def  ab      X(0)    Ω_X(0)  /%.3f  Ω_X z=0,5  z=1     z=2      w(0)    w₀      w_a     H/H_Λ z=0   0,5    1      2      3     Ω_s+Ω_X"
  % OM_X)
for dfn in DEFS:
    for ab in ("I_F", "Umkehr"):
        SEL['d'] = dfn
        d = V['lauf'](6, "oben", "0", ab); e = V['auswertung'](d)
        hz = {z: float(np.interp(1 / (1 + z), AG, d['Eeff'])) / E_L(1 / (1 + z)) for z in (0, 0.5, 1, 2, 3)}
        ERG[(dfn, ab)] = (d, e, hz)
        P("   %-3s  %-6s  %6.3f  %7.4f  %6.3f   %7.4f   %7.4f  %7.4f   %6.3f  %6.3f  %7.3f   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f   %6.3f"
          % (dfn, ab, e['X0'], e['OX0'], e['OX0'] / OM_X, ox_bei(d, 0.5), ox_bei(d, 1), ox_bei(d, 2), e['w_heute'], e['w0'], e['wa'],
             hz[0], hz[0.5], hz[1], hz[2], hz[3], OM_S + e['OX0']))
P("   Ω_DE gemessen [G] %.3f; Ω_s = %.3f fest. w₀, w_a: Anpassung z ≤ 2 wie de_fluchtmoden_verzoegert." % (OMEGA_DE, OM_S))
P("   E (linear) treibt das Atmen um Größenordnungen; H/H_Λ(0) ≫ 1 ist keine Lösung, nur Grenze.")
# Bedarf
d_, e_, _ = ERG[("B", "I_F")]
ratio_need = None
try:
    base = e_['X0']
    def ox_mit_faktor(k):
        SEL['d'] = "B"
        old = V['X_zaehlung']
        V['X_zaehlung'] = lambda a, l, art: k * X_def(a, l, art)
        dd = V['lauf'](6, "oben", "0", "I_F")
        V['X_zaehlung'] = old
        return V['auswertung'](dd)['OX0']
    ratio_need = brentq(lambda k: ox_mit_faktor(k) - OM_X, 0.5, 4.0, xtol=1e-3)
    P("   Für Ω_X(0) = %.3f bräuchte B (I_F) den Faktor %.2f auf X (X(0) = %.2f statt %.2f je Baryon)."
      % (OM_X, ratio_need, ratio_need * base, base))
except Exception as ex:
    P("   Bedarf nicht bestimmbar: %s" % ex)

# ============================================================== 3 Gesamtbild
kopf("3  Gesamtbild (erzwungene Definition B, und A zum Vergleich)")
def D_M(Efun, z=1090.0):
    zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, z, 4001)[1:]])
    return np.trapezoid(1 / np.array([Efun(1 / (1 + x)) for x in zs]), zs)
dm_L = D_M(E_L)
def wachs(Efun):
    def f(x, y):
        a = math.exp(x); h = 1e-4
        dl = (math.log(Efun(a * math.exp(h))) - math.log(Efun(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Efun(a)**2 * y[0]]
    s = solve_ivp(f, (math.log(AG[0]), 0.0), [AG[0], AG[0]], dense_output=True, rtol=1e-8, atol=1e-12)
    return lambda a: s.sol(math.log(a))
gL = wachs(E_L)
for key in (("B", "I_F"), ("B", "Umkehr"), ("A", "I_F"), ("C", "I_F")):
    d, e, hz = ERG[key]
    Ef = lambda a, d=d: float(np.interp(a, AG, d['Eeff'])) if a >= AG[0] else EB(a)
    dm = D_M(Ef)
    L_min = KREISPAAR_GRENZE * d_LSS() * dm / dm_L
    gE = wachs(Ef)
    fs = ", ".join("z = %.1f: %.3f" % (z, gE(1 / (1 + z))[1] / gL(1 / (1 + z))[1]) for z in (0.0, 0.3, 0.6, 1.0, 1.5))
    P("   %s ab %s: Ω_s + Ω_X = %.3f (Ω_DE %.3f [G]); D_M(1090)/D_M,ΛCDM = %.3f; Q-Fenster oben %.2f (ΛCDM %.2f), unten 5"
      % (key[0], key[1], OM_S + e['OX0'], OMEGA_DE, dm / dm_L, q_fenster(L_min=L_min)[1], q_fenster()[1]))
    P("      fσ₈/fσ₈,ΛCDM (gleiche frühe Amplitude): " + fs)

# ============================================================== 4 Blick
kopf("4  Beobachterfehler")
P("   μ ist das Maß des Zustandsraums eines Bereichs (G des Bereichs), gezählt in ganzen Tori; über die Zelle je Baryon (G).")
P("   Der Faktor 1/(1 − r²/R²) ist die Streckung des Geschwindigkeits-Ellipsoids (Tori, die die Naht nicht erreichen), keine")
P("   Dichte am Ort und kein Raumvolumen; das Mehr sitzt in fast kreisförmigen Bahnen nahe der Naht, wo kaum Baryonen sind.")
P("   (1 + Q)·M_b als Bezug (D) läse die Mode ξ innen als bindenden Stoff; innen ist w = 0, die Mode ist Regel (ν).")
P("   X, Ω_X zustandsabhängig (Epoche, Massenfunktion); τ_D und Λ fest.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand: Newton-Tori analytisch, Regel über den Wirkungsweg (Milchstraße, M_b = 6·10¹⁰ M☉, a = 1)")
In = quad(lambda x: math.sqrt(x * (1 - x)) / (1 + x), 0, 1)[0]
P("   Newton, Punktmasse: μ_b = (16π²/3)(2GM R)^{3/2}·∫₀¹√(x(1−x))/(1+x)dx, Integral = %.5f  [= π(3/2 − √2) = %.5f]"
  % (In, math.pi * (1.5 - math.sqrt(2))))
P("   gegen bis R (A'): ∫₀¹x^{1/2}(1−x)^{3/2}dx = B(3/2, 5/2) = π/16 = %.5f ⇒ μ_b/μ_a (Newton) = 16(3/2 − √2) = %.4f"
  % (math.pi / 16, 16 * (1.5 - math.sqrt(2))))
Mw = 6e10 * M_sun; R = dfm.R_ta(Mw, 1.0); rs = r_stern(Mw)
mw = mass_bereich(Mw, R)
P("   Code Newton μ_b/μ_a = %.4f." % mw['_rat']['Nb_NRa'])
IR = quad(lambda u: u**2 * math.log(1 / u)**1.5 / (1 - u**2), 0, 1, limit=200)[0]
P("   Regel tief (v_R² = 2v_f² ln(R/r)): μ_R,b/μ_N,b = (R/r*)^{3/2}·∫u²(ln 1/u)^{3/2}/(1−u²)du / ∫√(x(1−x))/(1+x)dx")
P("     = %.2f^{1,5}·%.5f/%.5f = %.1f (Code, voller Kern-ν: %.1f; Rest aus ν ≈ 1/s + ½)"
  % (R / rs, IR, In, (R / rs)**1.5 * IR / In, mw['_rat']['B']))
# Wirkungsweg: μ = (2π)³∫₀^{L_c} 2L·J_r(Φ(R) + L²/2R², L) dL
r_, vR2 = v2_regel(Mw, R, 4000)
Phi = lambda x: -0.5 * float(np.interp(x, r_, vR2))          # Φ(r) − Φ(R)
def Jr(E, L, Ph, acc):
    f = lambda x: 2 * (E - Ph(x)) - L**2 / x**2
    rc = brentq(lambda x: L**2 / x**3 - acc(x), R * 1e-9, R) if L > 0 else R * 1e-9
    if f(rc) <= 0:
        return 0.0
    rp = brentq(f, R * 1e-12, rc) if L > 0 else R * 1e-12
    ra = brentq(f, rc, R * (1 - 1e-12)) if f(R * (1 - 1e-12)) < 0 else R
    return quad(lambda x: math.sqrt(max(f(x), 0)), rp, ra, limit=200)[0] / math.pi
def mu_wirkung(Ph, acc):
    Lc = math.sqrt(R**3 * acc(R))
    Ls = np.linspace(0, Lc, 161)[1:-1]
    vals = [2 * L * Jr(Ph(R) + L**2 / (2 * R**2), L, Ph, acc) for L in Ls]
    return (2 * math.pi)**3 * np.trapezoid([0] + vals + [0], np.concatenate([[0], Ls, [Lc]]))
accR = lambda x: float(nu(G * Mw / x**2) * G * Mw / x**2)
muR_w = mu_wirkung(Phi, accR)
PhiN = lambda x: -G * Mw * (1 / x - 1 / R)
muN_w = mu_wirkung(PhiN, lambda x: G * Mw / x**2)
r2, _ = v2_regel(Mw, R, 1500)
x2 = r2 / R; wt2 = np.where(x2 < 1, K6 * r2**2 / np.clip(1 - x2**2, 1e-300, None), 0)
muN_code = float(np.trapezoid(wt2 * np.clip(2 * G * Mw * (1 / r2 - 1 / R), 0, None)**1.5, r2))
muN_ana = K6 * (2 * G * Mw * R)**1.5 * In
P("   Wirkungsweg (2π)³∫2L·J_r dL (anderer Weg, Quadratur in E und L):")
P("     Newton μ_b: Wirkung %.4e, Ellipsoid-Code %.4e, analytisch %.4e" % (muN_w, muN_code, muN_ana))
P("     Regel/Newton μ_b: Wirkung %.2f, Code %.2f; ln %.3f gegen %.3f" % (muR_w / muN_w, mw['_rat']['B'], math.log(muR_w / muN_w), mw['B']))
Fb = float(erfc(dfm.delta_ta(1.0) / (math.sqrt(2) * np.interp(6.0, dfm.LM, dfm.SIG))))
sh = math.log(32 / (3 * math.pi)) + math.log((IR / (gamma(2.5) / 3**2.5)) / (16 * (1.5 - math.sqrt(2))))
P("   Verschiebung je Bereich, tief: ln B − ln A = ln(32/3π) + ln[(%.5f/%.5f)/1,3726] = 1,2225 + %.4f = %.3f (Tabelle 1: 1,39)"
  % (IR, gamma(2.5) / 3**2.5, sh - math.log(32 / (3 * math.pi)), sh))
P("   ⇒ X_B(0) − X_A(0) ≈ %.3f·F(>10⁶, 0) = %.3f·%.3f = %.3f (Code %.3f)." % (sh, sh, Fb, sh * Fb, ERG[("B", "I_F")][1]['X0'] - ERG[("A", "I_F")][1]['X0']))
d_, e_, _ = ERG[("B", "I_F")]
P("   Grob: Ω_X ≈ (E_b + ΔE)² − E_b², ΔE ∝ X; B/A in X(0) = %.3f/%.3f = %.2f ⇒ ΔE_B ≈ %.2f·0,137 = %.3f,"
  % (e_['X0'], ERG[("A", "I_F")][1]['X0'], e_['X0'] / ERG[("A", "I_F")][1]['X0'], e_['X0'] / ERG[("A", "I_F")][1]['X0'],
     0.137 * e_['X0'] / ERG[("A", "I_F")][1]['X0']))
dE = 0.137 * e_['X0'] / ERG[("A", "I_F")][1]['X0']
P("   Ω_X ≈ (0,747 + %.3f)² − 0,747² = %.3f (Code %.3f; Zeitverlauf der Quelle ungleich, daher nur grob)."
  % (dE, (0.747 + dE)**2 - 0.747**2, e_['OX0']))

with open(os.path.join(HIER, "de_fluchtmoden_mass.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
