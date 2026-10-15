# -*- coding: utf-8 -*-
"""
Inaktives Fluchtmaß im hierarchischen Aufbau: Wird ein kleiner Bereich Teil eines größeren, ist sein Maß im Maß
des größeren enthalten (geschachtelt) oder zählt es zusätzlich (unabhängig)?
1 Herleitung der Zählung aus dem Modell (Triade, Projektion in ϖ, Naht, Übergabe, Liouville/Kac).
2 X(a) je Baryon (Liouville-Produkt, normfrei: ln μ_R/μ_N des Bereichs), zwei Zählungen:
  (H) geschachtelt: jedes Baryon zählt im größten Bereich, der es zur Epoche a enthält (Erstüberschreitung im
      Ausflug-Bild [EPS]; der Bereich erreicht seine Umkehr gerade bei a, Naht R_ta(M, a)).
  (U) unabhängig: jeder je entstandene Bereich behält sein Maß (Liouville je Bereich fest), Entstehungsrate
      analytisch (Sasaki: Entstehung = dn/dt + Zerstörung, Zerstörung ∝ n·Ḋ/D) und im Verschmelzungsbaum
      (Monte-Carlo-Ausflüge, neue Generation bei Massenverdopplung).
  Atmen wie de_fluchtmoden: H_eff = H_b(1 + ⅓dX/dln a), Ω_X = E_eff² − E_b², w(a), CPL w₀, w_a.
3 Zusammenhänge: Kac 1 + Q, a₀, (Ω_s + Ω_X)/Ω_Mode gegen 2,60 [G].  4 Beobachterfehler.  5 Von Hand.
Formeln aus ftm_kern (über die Hilfsfunktionen aus de_fluchtmoden: masse, R_ta, umkehr_dichte, wachstum, sig0,
delta_ta – dort lokal, nicht im Kern). σ(M) aus CAMB [G].
"""
import math, os, sys
import numpy as np
from scipy.special import erfc, gamma
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)

# ---- Hilfsfunktionen aus de_fluchtmoden (nur die Definitionen, ohne dessen Tabellen und Ausgabe)
_src = open(os.path.join(HIER, "de_fluchtmoden.py"), encoding="utf-8").read().split("\n")
_ns = {"__file__": os.path.join(HIER, "de_fluchtmoden.py")}
exec("\n".join(_src[:101]), _ns)                      # Kopf, Umkehrdichte, R_ta, Hz, masse
exec("\n".join(_src[141:172]), _ns)                   # wachstum, CAMB σ(M), δ_ta(a)
masse, Dg, sig0, delta_ta = _ns['masse'], _ns['Dg'], _ns['sig0'], _ns['delta_ta']
LM, SIG = _ns['LM'], _ns['SIG']
from ftm_kern import (M_sun, OMEGA_M, OMEGA_B, OMEGA_R, OMEGA_DE, Q_KOSM, q_KOSM, E_hintergrund, omega_vakuum,
                      r_stern, A0, c)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

NK = 1 + Q_KOSM
OM_S = omega_vakuum(); OM_X = OMEGA_DE - OM_S; OM_C = OMEGA_M - OMEGA_B
LN10 = math.log(10)

# ---- Tabelle ln(μ_R/μ_N)(M_ges, a) je Bereich (Punktmasse, Naht R_ta(M, a))
AA = np.exp(np.linspace(math.log(1 / 21), math.log(1.1), 34))
LMg = np.linspace(6.0, 16.5, 22)
LR = np.array([[math.log((lambda m: m['muR'] / m['muN'])(masse(10**l * M_sun / NK, a, n=1500))) for l in LMg] for a in AA])
LA = np.log(AA)
def lnr(lm, a):
    """bilinear in (log M, ln a)"""
    la = math.log(a); j = np.clip(np.searchsorted(LA, la) - 1, 0, len(LA) - 2); t = (la - LA[j]) / (LA[j + 1] - LA[j])
    return (1 - t) * np.interp(lm, LMg, LR[j]) + t * np.interp(lm, LMg, LR[j + 1])

def dlnD(a, h=1e-4):
    return (math.log(Dg(a * math.exp(h))) - math.log(Dg(a * math.exp(-h)))) / (2 * h)

def mf(a):
    """F(>M) und dF/dlnM (PS/EPS mit Umkehrschwelle δ_ta), ν = δ_ta/(σD)"""
    nu_ = delta_ta(a) / (SIG * Dg(a))
    Fg = erfc(nu_ / math.sqrt(2))
    return Fg, -np.gradient(Fg, LM * LN10), nu_

AF = np.exp(np.linspace(math.log(1 / 21), math.log(1.1), 240))     # feines a-Gitter für Integrale
def X_H(a, lmin):
    Fg, dF, _ = mf(a); m = LM >= lmin
    return float(np.trapezoid((dF * lnr(LM, a))[m], LM[m] * LN10))
def rate_U(a, lmin):
    Fg, dF, nu_ = mf(a); m = LM >= lmin
    return float(np.trapezoid((dF * nu_**2 * lnr(LM, a))[m], LM[m] * LN10)) * dlnD(a)
def rate_F(a, lmin):                                                # Entstehungsrate der Masse (Sasaki), ohne Gewicht
    Fg, dF, nu_ = mf(a); m = LM >= lmin
    return float(np.trapezoid((dF * nu_**2)[m], LM[m] * LN10)) * dlnD(a)

EB = lambda a: math.sqrt(OMEGA_M * a**-3 + OMEGA_R * a**-4 + OM_S)
EL = lambda a: E_hintergrund(1 / a - 1)
LAF = np.log(AF); i1 = int(np.argmin(abs(AF - 1.0)))
def atmen(X):
    d = np.gradient(X, LAF)
    Eb = np.array([EB(a) for a in AF]); Ee = Eb * (1 + d / 3); OX = Ee**2 - Eb**2
    return d, Ee, OX
def w_von(OX):
    with np.errstate(divide='ignore', invalid='ignore'):
        return -1 - np.gradient(np.log(np.where(OX > 0, OX, np.nan)), LAF) / 3
def cpl(w):
    m = (AF >= 0.5) & (AF <= 1.0) & np.isfinite(w)
    A = np.vstack([np.ones(m.sum()), 1 - AF[m]]).T
    return np.linalg.lstsq(A, w[m], rcond=None)[0]
zi = lambda z: int(np.argmin(abs(AF - 1 / (1 + z))))

# ============================================================== 1
kopf("1  Welche Zählung verlangt das Modell? (Herleitung, nicht nach Ergebnis)")
for s in (
 "   (a) Ein Zustandsraum (ein Medium, kein Außen): das Einteilchenmaß ist ein Maß auf EINEM Raum; additiv nur über",
 "       disjunkte Mengen. Nach dem Einfall liegt jeder Zustand, der im Kind gebunden war, auch im Eltern gebunden",
 "       (tieferes Potential, Kind innerhalb R_ta der Eltern): A_Kind ⊂ A_Eltern, μ(A_Kind ∪ A_Eltern) = μ(A_Eltern).",
 "   (b) Naht/Übergabe: μ_inakt ist relativ zur Naht definiert (Regel bis R_ta, draußen Fluss). Innerhalb R_ta der",
 "       Eltern ist w = 0, dort gilt die volle Regel der Eltern; die Naht des Kindes (Schale bis 1,16·R_ta) hat keinen",
 "       Fluss mehr außen – sie löst sich auf, das Kind hat kein eigenes Außenmaß mehr, gegen das es inaktiv wäre.",
 "   (c) Triade/Projektion in ϖ: aufsteigend bleiben Fluss GM (Zentrum) und Elternhintergrund, die eigene aufgelöste",
 "       Struktur klingt ab (Gezeit ∝ ϖ⁻³). Das inaktive Maß sitzt fast ganz zwischen r* und R_ta (Struktur der Regel,",
 "       'Wirbel'); jenseits der Übergangstiefe ℓ ≈ ϖ_J gehört dieser Teil zum Hintergrund der Eltern. Von oben liest",
 "       man das Kind nur als Fluss (1 + wQ)·M_b, der in M_b der Eltern eingeht (Δs_ges = ΣΔs_j).",
 "   (d) 𝒫 löscht keine Freiheitsgrade: die inaktiven Richtungen des Kindes verschwinden nicht, sie werden Richtungen",
 "       des Elternbereichs (dieselben Zustände), sie werden nicht ein zweites Mal inaktiv.",
 "   (e) Kac-Budget Σ1/⟨T_A⟩ = 1 gilt für eine Zerlegung; geschachtelte Mengen sind keine Zerlegung.",
 "   ⇒ Das Modell verlangt die GESCHACHTELTE Zählung (H): je Baryon das Maß des größten Bereichs, der es enthält.",
 "     Kompakte Körper im Raum geschachtelt: innen behält das Kind bis ℓ seine Regel, das ist eine Verfeinerung",
 "     innerhalb der Elternmenge, kein zusätzliches Maß. Die unabhängige Zählung (U) zählte dieselben Zustände je",
 "     Generation einmal mehr; sie wird nur zum Vergleich gerechnet."):
    P(s)

# ============================================================== 2
kopf("2  X(a) je Baryon: geschachtelt (H) gegen unabhängig (U); Atmen H_eff = H_b(1 + ⅓dX/dln a)  [Messwerte nur Test, G]")
P("   (H) X_H(a) = ∫_{M ≥ M_min} dF/dlnM(a)·ln(μ_R/μ_N)(M, a) dlnM   (größter Bereich erreicht die Umkehr gerade bei a)")
P("   (U) dX_U/dln a = ∫ dF/dlnM·ν²·dlnD/dlna·ln(μ_R/μ_N) dlnM   (Sasaki-Entstehungsrate, jedes Maß bleibt)")
P("   E_b(1) = %.3f; für E_eff(1) = 1 nötig dX/dln a (a = 1) = %.3f" % (EB(1.0), 3 * (1 / EB(1.0) - 1)))
P()
P("   M_min  Zählung  X(z=3)  X(1)   X(0)   dX/dlna(0)  Gipfel dX/dlna bei z (Wert)  Ω_X(0)  /0,442  E_eff(1)  w₀      w_a     E/E_ΛCDM z=1  z=2   z=3    (Ω_s+Ω_X)/Ω_c")
RES = {}
for lmin in (6, 8, 10, 12):
    XH = np.array([X_H(a, lmin) for a in AF])
    rU = np.array([rate_U(a, lmin) for a in AF])
    XU = np.concatenate([[XH[0]], XH[0] + np.cumsum(0.5 * (rU[1:] + rU[:-1]) * np.diff(LAF))])   # Start: gleich H bei z = 20
    rF = np.array([rate_F(a, lmin) for a in AF])
    RES[lmin] = dict(XH=XH, XU=XU, rF=rF)
    for nm, X in (("H", XH), ("U", XU)):
        d, Ee, OX = atmen(X); w = w_von(OX); w0, wa = cpl(w)
        m = AF <= 1.0; ip = int(np.argmax(np.where(m, d, -1e9)))
        RES[lmin][nm] = (d, Ee, OX, w)
        P("   10^%-3d   %s     %6.3f  %6.3f  %6.3f   %6.3f       z = %5.2f  (%.3f)          %6.3f  %5.2f   %.3f    %6.3f  %6.3f      %.3f    %.3f  %.3f     %.2f"
          % (lmin, nm, X[zi(3)], X[zi(1)], X[i1], d[i1], 1 / AF[ip] - 1, d[ip], OX[i1], OX[i1] / OM_X, Ee[i1], w0, wa,
             Ee[zi(1)] / EL(AF[zi(1)]), Ee[zi(2)] / EL(AF[zi(2)]), Ee[zi(3)] / EL(AF[zi(3)]), (OM_S + OX[i1]) / OM_C))
P()
P("   Verlauf dX/dln a (H | U), M_min = 10⁸ und 10¹²:")
for lmin in (8, 12):
    s = "   10^%-3d" % lmin
    for z in (10, 6, 4, 3, 2, 1, 0.5, 0):
        s += "  z=%-3g %.2f|%.2f" % (z, RES[lmin]['H'][0][zi(z)], RES[lmin]['U'][0][zi(z)])
    P(s)
P("   Durchlauf je Baryon (U, Zahl der Generationen = ∫ Entstehungsrate der Masse dln a bis heute):")
for lmin in (6, 8, 10, 12):
    rF = RES[lmin]['rF']; n = np.concatenate([[0], np.cumsum(0.5 * (rF[1:] + rF[:-1]) * np.diff(LAF))])
    P("   10^%-3d  N_gen(0) = %.2f   F(>M_min, 0) = %.3f" % (lmin, n[i1], float(np.interp(lmin, LM, mf(1.0)[0]))))

# ---- Verschmelzungsbaum: Monte-Carlo-Ausflüge (scharfer k-Filter in S = σ²), anderer Weg für (H) und (U)
P()
P("   Verschmelzungsbaum (Monte Carlo, 20000 Ausflüge in S = σ²(M), Schranke δ_ta(a)/D(a), 240 Epochen):")
rng = np.random.default_rng(1)
LMf = np.linspace(16.5, 6.0, 1400); Sf = np.interp(LMf, LM, SIG)**2       # von groß nach klein, S wachsend
NW = 20000
steps = (rng.standard_normal((NW, len(Sf))) * np.sqrt(np.diff(np.concatenate([[0.0], Sf])))).astype(np.float32)
runmax = np.maximum.accumulate(np.cumsum(steps, axis=1), axis=1); del steps
AMC = AF[AF <= 1.0 + 1e-12]
host = np.full((NW, len(AMC)), -np.inf)
for k, a in enumerate(AMC):
    B = delta_ta(a) / Dg(a)
    hit = runmax >= B
    first = np.where(hit.any(axis=1), hit.argmax(axis=1), -1)
    host[:, k] = np.where(first >= 0, LMf[np.clip(first, 0, None)], -np.inf)
for lmin in (8, 12):
    XHm = np.zeros(len(AMC)); XUm = np.zeros(len(AMC))
    acc = np.zeros(NW); gen_m = np.full(NW, -np.inf); ngen = np.zeros(NW)
    for k, a in enumerate(AMC):
        h = host[:, k]; ok = h >= lmin
        lr = np.where(ok, lnr(np.where(ok, h, lmin), a), 0.0)
        XHm[k] = lr.mean()
        new = ok & (h >= gen_m + math.log10(2))                          # neue Generation: Masse verdoppelt
        acc += np.where(new, lr, 0.0); gen_m = np.where(new, h, gen_m); ngen += new
        XUm[k] = acc.mean()
    kk = [int(np.argmin(abs(AMC - 1 / (1 + z)))) for z in (3, 1, 0)]
    P("   10^%-3d  MC (H): X(3) %.3f  X(1) %.3f  X(0) %.3f  [analytisch %.3f %.3f %.3f]   MC (U, Verdopplung): X(3) %.3f  X(1) %.3f  X(0) %.3f  (Generationen je Baryon %.2f)"
          % (lmin, *XHm[kk], RES[lmin]['XH'][zi(3)], RES[lmin]['XH'][zi(1)], RES[lmin]['XH'][i1], *XUm[kk], ngen.mean()))
    dH = np.gradient(XHm, np.log(AMC)); dU = np.gradient(XUm, np.log(AMC))
    P("           MC dX/dlna(0): H %.3f, U %.3f; Gipfel H bei z = %.2f, U bei z = %.2f"
      % (dH[-1], dU[-1], 1 / AMC[int(np.argmax(dH))] - 1, 1 / AMC[int(np.argmax(dU))] - 1))

P()
P("   Vorzeichen in (H): X_H steigt, weil R_ta(M, a)/r* mit a wächst (Naht des neuen, größeren Bereichs liegt weiter");
P("   draußen), obwohl ln(μ_R/μ_N) bei fester Epoche mit M fällt (∝ −¼ ln M). Beim Einfall sinkt das Maß je Baryon nicht")
P("   zwingend: Kind 10⁸ bei a = 0,25 und Milchstraße bei a = 1 (Abschnitt 5).")

# ============================================================== 3
kopf("3  Zusammenhänge")
P("   Kac: 1 + Q = %.3f zählt die Rückkehr in I (Mode ξ als Staub, Ω_c); X ist das Bahnmaß des Bereichs (L/G). Q geht nur" % NK)
P("   über M_ges = (1 + Q)M_b in R_ta ein; X schreibt nicht auf Q zurück: 1 + Q und q = %.3f unberührt." % q_KOSM)
P("   a₀ = c²/2Λ, Λ = (1 + Q)ΣΔs_b: X ändert keine Faltung, Λ und a₀ = %.2e m/s² fest. Ω_s = π²/45·Θ² als Dichte fest;" % A0)
P("   verschiebt X aber H₀ (E_eff(1) ≠ 1), so verschiebt es Θ = 2πa₀/(cH₀) – Rückwirkung nur über H, nicht über a₀.")
P("   (Ω_s + Ω_X)/Ω_Mode gegen %.2f [G] (Ω_s = %.3f, Ω_c = %.4f): Tabelle Abschnitt 2, letzte Spalte." % (OMEGA_DE / OM_C, OM_S, OM_C))

# ============================================================== 4
kopf("4  Beobachterfehler")
for s in (
 "   μ ist das Maß des Zustandsraums eines Bereichs (G des Bereichs), über die Zelle als Mittel je Baryon summiert (G).",
 "   Der Ebenenwechsel ist Projektion in ϖ: das Kind wird von oben als Fluss gelesen, sein Maß geht als Teilmenge in das",
 "   des Elternbereichs auf – nichts wird am Ort erzeugt oder vernichtet, keine Dichte, kein Stoff.",
 "   X ist zustandsabhängig (Epoche, Massenfunktion, M_min), keine feste Zahl."):
    P(s)

# ============================================================== 5
kopf("5  Von Hand")
rs = r_stern(1e8 * M_sun); R = _ns['R_ta'](1e8 * M_sun, 1.0)
P("   Tiefformel 10⁸ M☉, a = 1: R_ta/r* = %.1f, 0,1279·%.1f^1,5 = %.1f, ln = %.3f (Code, voller Kern-ν: %.3f)"
  % (R / rs, R / rs, 0.12792 * (R / rs)**1.5, math.log(0.12792 * (R / rs)**1.5), float(lnr(8 + math.log10(NK), 1.0))))
l1 = float(lnr(8 + math.log10(NK), 0.25)); l2 = float(lnr(10.778 + math.log10(NK), 1.0))
P("   Einfall 10⁸ (gebunden a = 0,25, ln = %.3f) in Milchstraße (a = 1, ln = %.3f): (H) %.3f, (U) %.3f je Baryon" % (l1, l2, l2, l1 + l2))
nu8 = delta_ta(1.0) / (float(np.interp(12, LM, SIG)) * Dg(1.0))
P("   Sasaki-Faktor ν² bei 10¹² M☉ heute: δ_ta/σ = 1,2766/%.3f = %.3f, ν² = %.3f; dlnD/dlna(1) = %.3f"
  % (float(np.interp(12, LM, SIG)), nu8, nu8**2, dlnD(1.0)))

with open(os.path.join(HIER, "de_fluchtmoden_hierarchie.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
