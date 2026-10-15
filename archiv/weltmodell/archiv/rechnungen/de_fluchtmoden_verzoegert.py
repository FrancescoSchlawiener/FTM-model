# -*- coding: utf-8 -*-
"""
Verzögerter Übergang des inaktiven Fluchtmaßes ins Atmen des Ganzen (Folge des Axioms: geschlossen, kein Außen,
𝒫 löscht keine Freiheitsgrade ⇒ die inaktiven Richtungen gehen ins Ganze über; der Weg ist die Rückkehr in der
inneren Richtung I mit der Dichte erster Rückkehrzeiten, τ_D = c/a₀ = 2Λ/c).
Blick: X je Baryon ist ein Maß des Zustandsraums je Bereich (G des Bereichs); der Übergang läuft in I (Rückkehrzeit T,
keine Strecke im Raum); angekommen wirkt er aufs Atmen (G, Zelle): H_eff = H_b + ⅓·R(t). Nirgends eine Dichte am Ort.
1 Rückkehrdichte f(T) aus dem Kern (Laplace-Inversion von rueckkehr(p, 0, θ)), Faltung R(t) = ∫ S(t') f(t − t') Ü(t', t) dt',
  S = dX/dt (Liouville-Produkt je Baryon, de_fluchtmoden); Gipfel, Ω_X(0), w₀, w_a, H(z).
2 Abhängigkeiten: M_min, Zählung (oberste Ebene / je Ebene), Verlust (θ̃ = 0 / 3H), Übergabe ab Umkehr / ab I_F.
3 Zusammenhänge: Λ, a₀, CMB, Q-Fenster, Wachstum fσ₈.  4 Blick.  5 Von Hand.
Aus ftm_kern: rueckkehr, tau_D, schliesslaenge, E_hintergrund, omega_vakuum, q_fenster, d_LSS, Konstanten.
Aus de_fluchtmoden: masse, X_von_a (mit XT), Dg, delta_ta, LM, SIG, EB.
Nicht im Kern (lokal, gemeldet): Talbot-Inversion, Faltung, Wachstum mit H_eff, t_F/t_ta = 1,8 (uebergabe_faltung, 1,71–1,89).
"""
import math, os, sys, io, contextlib
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import erfc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (A0, c, M_sun, Mpc, Gyr, OMEGA_M, OMEGA_R, OMEGA_DE, OMEGA_B, H0_si, E_hintergrund,
                      omega_vakuum, rueckkehr, tau_D, schliesslaenge, q_fenster, d_LSS, KREISPAAR_GRENZE)
with contextlib.redirect_stdout(io.StringIO()):
    import de_fluchtmoden as dfm

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

H0 = H0_si()
TD = tau_D() / Gyr                           # Gyr
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
E_L = lambda a: E_hintergrund(1 / a - 1)
EB = dfm.EB
T_F = 1.8                                    # t_F/t_ta (uebergabe_faltung: 1,71–1,89)

# ------------------------------------------------ Zeit längs des Hintergrunds (Strukturbildung wie in de_fluchtmoden)
def t_von_a(a):
    return quad(lambda x: 1 / (x * E_L(x)), 1e-9, a, limit=200)[0] / (H0 * Gyr)
T0 = t_von_a(1.0)

# ------------------------------------------------ 1 Rückkehrdichte aus dem Kern
def talbot(Fhat, t, M=48):
    """Inverse Laplace (Abate–Valkó, fester Talbot-Weg); Fhat(p) mit p in 1/Gyr."""
    r = 2 * M / (5 * t)
    k = np.arange(1, M)
    th = k * math.pi / M
    cot = 1 / np.tan(th)
    d = r * th * (cot + 1j)
    g = 1 + 1j * th * (1 + cot**2) - 1j * cot
    s = 0.5 * math.exp(r * t) * Fhat(r) + np.sum((np.exp(t * d) * Fhat(d) * g).real)
    return float((r / M * s).real)

def Khat(p, theta=0.0):      # p, θ in 1/Gyr → Kern rueckkehr(p, Γ = 0, θ) in SI
    return rueckkehr(np.asarray(p, dtype=complex) / Gyr, 0.0, theta / Gyr)

TT = np.concatenate([np.geomspace(0.05, 1, 40)[:-1], np.linspace(1, 60, 600)])
FK = np.array([talbot(lambda p: Khat(p), t) for t in TT])
fdich = lambda T: np.interp(T, TT, FK, left=0.0, right=0.0)
f_form = lambda T: math.sqrt(TD / (4 * math.pi)) * T**-1.5 * math.exp(-TD / (4 * T))
CDF = lambda t: float(talbot(lambda p: Khat(p) / p, t))

kopf("1  Rückkehrdichte f(T) aus dem Kern und verzögerter Übergang")
P("   τ_D = c/a₀ = 2Λ/c = %.2f Gyr (Λ = %.3f Gpc).  f(T) = L⁻¹[e^{−√(pτ_D)}] (Talbot) gegen √(τ_D/4π)T^{−3/2}e^{−τ_D/4T} (07):"
  % (TD, schliesslaenge() / (1e3 * Mpc)))
P("   T [Gyr]:      " + "  ".join("%7.2f" % t for t in (2, 5, 10, TD / 6, 13.8, 30, 60)))
P("   Kern (Talbot):" + "  ".join("%7.4f" % float(talbot(lambda p: Khat(p), t)) for t in (2, 5, 10, TD / 6, 13.8, 30, 60)))
P("   Formel 07:    " + "  ".join("%7.4f" % f_form(t) for t in (2, 5, 10, TD / 6, 13.8, 30, 60)))
Tm = TT[np.argmax(FK)]
P("   Maximum bei T = %.2f Gyr (τ_D/6 = %.2f). Angekommen nach t: ∫₀ᵗf = erfc(√(τ_D/4t)):" % (Tm, TD / 6))
for t in (5.0, 10.0, T0, 30.0, 100.0):
    P("     t = %6.2f Gyr: Kern %.4f, erfc %.4f" % (t, CDF(t), erfc(math.sqrt(TD / (4 * t)))))
th0 = 3 * H0 * Gyr
P("   Mit Verlust θ = 3H₀ (fest): gesamt ∫f·e^{−θT} = K̂(0; θ) = e^{−√(3H₀τ_D)} = %.4f (3H₀τ_D = %.1f); Talbot-Summe %.4f."
  % (float(Khat(0, th0).real), 3 * H0 * tau_D(),
     np.trapezoid([talbot(lambda p: Khat(p, th0), t) for t in TT], TT)))
P("   Zeitabhängiges θ = 3H(t) im Ganzen: räumlich gleich ⇒ Faktor e^{−∫θdt} = (a(t')/a(t))³ auf f (exakt).")

# ------------------------------------------------ X(a) auf einem Gitter bis z = 30 und in die Zukunft
def _wachstum_bis(a_end=2.0):      # wie dfm.wachstum, aber über a = 1 hinaus integriert (keine Extrapolation)
    def f(x, y):
        a = math.exp(x); h = 1e-5
        dlnE = (math.log(E_L(a * math.exp(h))) - math.log(E_L(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dlnE) * y[1] + 1.5 * OMEGA_M * a**-3 / E_L(a)**2 * y[0]]
    s = solve_ivp(f, (math.log(1e-3), math.log(a_end)), [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-12)
    D1 = s.sol(0.0)[0]
    return lambda a: s.sol(math.log(a))[0] / D1
_Dg_alt = dfm.Dg
dfm.Dg = _wachstum_bis(3.5)
assert abs(dfm.Dg(0.5) / _Dg_alt(0.5) - 1) < 1e-6
AA =np.exp(np.linspace(math.log(1 / 31), math.log(3.0), 56))
for a in AA:
    if a not in dfm.XT:
        rows = []
        for l in dfm.LMg:
            m = dfm.masse(10**l * M_sun / dfm.NK, a, n=1500)
            rows.append((m['muI'] / m['muF'], math.log(m['muR'] / m['muN']), m['muI'] / m['muR']))
        dfm.XT[a] = np.array(rows)
LN10 = math.log(10)

def X_zaehlung(a, lm_min, art):
    """oben: X_B = ∫(−dF/dlnM)·L dlnM (jedes Baryon in seinem größten Bereich);
    ebene: ∫F(>M)·L dlnM/ln10 (jedes Baryon in jedem enthaltenden Bereich, eine Ebene je Dekade)."""
    Fg = erfc(dfm.delta_ta(a) / (math.sqrt(2) * dfm.SIG * dfm.Dg(a)))
    L = np.interp(dfm.LM, dfm.LMg, dfm.XT[a][:, 1])
    msk = dfm.LM >= lm_min
    if art == "oben":
        dF = -np.gradient(Fg, dfm.LM * LN10)
        return float(np.trapezoid((dF * L)[msk], dfm.LM[msk] * LN10))
    return float(np.trapezoid((Fg * L)[msk], dfm.LM[msk] * LN10) / LN10)

TA = np.array([t_von_a(a) for a in AA])
NT = 3000
TG = np.linspace(TA[0], TA[-1], NT); dt = TG[1] - TG[0]
AG = np.exp(np.interp(TG, TA, np.log(AA)))
HG = H0 * Gyr * np.array([E_L(a) for a in AG])             # 1/Gyr (Hintergrund der Bildung)
EBG = np.array([EB(a) for a in AG])
i0 = int(np.argmin(abs(TG - T0)))
Fi = fdich(np.arange(NT) * dt)                              # f auf dem Zeitgitter

def lauf(lm_min, art, verlust, ab):
    X = np.array([X_zaehlung(a, lm_min, art) for a in AA])
    sp = CubicSpline(TA, X)
    tq = TG / T_F if ab == "I_F" else TG                     # Bindung bei t_ta ⇒ Übergabe bei 1,8·t_ta
    S = sp(np.clip(tq, TA[0], None), 1) * (tq >= TA[0]) / (T_F if ab == "I_F" else 1.0)
    X0 = X[0]; t_puls = TA[0] * (T_F if ab == "I_F" else 1.0)
    R = np.zeros(NT)
    lna3 = 3 * np.log(AG)
    for n in range(NT):
        w = Fi[n::-1] * S[:n + 1]
        if verlust == "3H":
            w = w * np.exp(lna3[:n + 1] - lna3[n])
        R[n] = np.trapezoid(w, dx=dt)
        tp = TG[n] - t_puls                                  # Anteil vor z = 30 als Puls
        if tp > 0:
            R[n] += X0 * fdich(tp) * ((AG[n] / AG[n]) if verlust == "0" else
                                      (np.interp(t_puls, TG, AG) / AG[n])**3)
    Xarr = cumulative_trapezoid(R, TG, initial=0.0)
    Eeff = EBG + R / (3 * H0 * Gyr)
    OX = Eeff**2 - EBG**2
    rate = R / HG                                            # dX_ank/d ln a
    return dict(X=X, S=S, R=R, Xarr=Xarr, Eeff=Eeff, OX=OX, rate=rate)

def auswertung(d):
    a = AG; OX = d['OX']
    g = int(np.argmax(d['rate']))
    zg = 1 / a[g] - 1
    msk = (a > 1 / 3.2) & (a < 1.05) & (OX > 0)
    lnr = np.log(np.clip(OX, 1e-300, None)); lna = np.log(a)
    wz = -1 - np.gradient(lnr, lna) / 3
    m2 = (a >= 1 / 3) & (a <= 1.0)
    A_ = np.vstack([np.ones(m2.sum()), 1 - a[m2]]).T
    w0, wa = np.linalg.lstsq(A_, wz[m2], rcond=None)[0]
    hz = {z: float(np.interp(1 / (1 + z), a, d['Eeff'])) / E_L(1 / (1 + z)) for z in (0, 1, 2, 3)}
    return dict(OX0=float(OX[i0]), zg=zg, rand=g >= len(a) - 2, OX1=float(np.interp(0.5, a, OX)), w0=w0, wa=wa, w_heute=float(wz[i0]), hz=hz, Xarr0=float(d['Xarr'][i0]),
                X0=float(np.interp(1.0, AA, d['X'])))

# ------------------------------------------------ Grundfall
kopf("1b Grundfall (M_min = 10⁸, oberste Ebene, θ̃ = 0, Übergabe ab Umkehr) gegen sofortige Übergabe (de_fluchtmoden)")
d = lauf(8, "oben", "0", "Umkehr"); e = auswertung(d)
S_ln = d['S'] / HG
gS = int(np.argmax(S_ln))
P("   Quelle dX/dln a: Gipfel bei z = %.2f (sofortige Übergabe).  Verzögert: Gipfel dX_ank/dln a bei z = %.2f%s."
  % (1 / AG[gS] - 1, e['zg'], "  (Zukunft, a = %.2f)" % AG[np.argmax(d['rate'])] if e['zg'] < 0 else ""))
P("   X(0) = %.3f inaktiv, davon bis heute angekommen %.3f (%.1f %%)." % (e['X0'], e['Xarr0'], 100 * e['Xarr0'] / e['X0']))
P("   z:           " + "  ".join("%7g" % z for z in (10, 6, 3, 2, 1, 0.5, 0, -0.3, -0.5, -0.6)))
for nm, arr in (("dX/dln a", S_ln), ("dX_ank/dln a", d['rate']), ("Ω_X(a)", d['OX']), ("E_eff/E_ΛCDM", d['Eeff'] / np.array([E_L(x) for x in AG]))):
    P("   %-13s" % nm + "  ".join("%7.4f" % float(np.interp(1 / (1 + z), AG, arr)) for z in (10, 6, 3, 2, 1, 0.5, 0, -0.3, -0.5, -0.6)))
P("   Ω_X(0) = %.4f (%.0f %% von %.3f), w heute %.3f, CPL-Anpassung z ≤ 2: w₀ = %.3f, w_a = %.3f; E_eff(0) = %.3f."
  % (e['OX0'], 100 * e['OX0'] / OM_X, OM_X, e['w_heute'], e['w0'], e['wa'], float(d['Eeff'][i0])))

# ------------------------------------------------ 2 Tabelle
kopf("2  Abhängigkeiten: Ω_X(0), Gipfel z, w₀, w_a, H/H_ΛCDM bei z = 1, 2, 3 [G], angekommen bis heute")
P("   M_min  Zählung  Verlust  ab      Ω_X(0)  /0,442  Ω_X(z=1)  z_Gipfel  w(0)    w₀      w_a     H/H_Λ z=0   z=1    z=2    z=3    X(0)    angek.")
TAB = {}
for lm in (6, 8, 10, 12, 13):
    for art in ("oben", "ebene"):
        for verl in ("0", "3H"):
            for ab in ("Umkehr", "I_F"):
                d_ = lauf(lm, art, verl, ab); e_ = auswertung(d_); TAB[(lm, art, verl, ab)] = (d_, e_)
                P("   10^%-3d %-7s  %-6s   %-6s  %7.4f  %6.3f   %6.4f   %6.2f%s  %6.3f  %6.3f  %7.3f   %6.3f  %6.3f  %6.3f  %6.3f   %6.3f  %6.4f"
                  % (lm, art, "θ=0" if verl == "0" else "θ=3H", ab, e_['OX0'], e_['OX0'] / OM_X, e_['OX1'], e_['zg'],
                     "*" if e_['rand'] else " ", e_['w_heute'], e_['w0'], e_['wa'], e_['hz'][0], e_['hz'][1], e_['hz'][2],
                     e_['hz'][3], e_['X0'], e_['Xarr0']))
P("   z_Gipfel < 0: Gipfel in der Zukunft; * = am Rand der Rechnung (a = 3), liegt noch später. Bildung nach heute mit")
P("   δ_ta(z = 0) fortgesetzt. w₀, w_a: Anpassung w(a) = w₀ + w_a(1 − a) über z ≤ 2 mit ρ_X ∝ Ω_X(a) = E_eff² − E_b².")
P("   'ebene': eine Ebene je Dekade der Masse (Wahl der Stufung; je halbe Dekade verdoppelt X).")

# ------------------------------------------------ 3 Zusammenhänge
kopf("3  Zusammenhänge (erzwungene Kombination und Spannweite)")
def D_M(Efun, z=1090.0):
    zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, z, 4001)[1:]])
    return np.trapezoid(1 / np.array([Efun(1 / (1 + x)) for x in zs]), zs)
dm_L = D_M(E_L)
for key in ((6, "ebene", "0", "I_F"), (8, "oben", "0", "I_F"), (6, "oben", "0", "I_F"), (12, "oben", "3H", "I_F")):
    d_, e_ = TAB[key]
    Ef = lambda a, d_=d_: float(np.interp(a, AG, d_['Eeff'])) if a >= AG[0] else EB(a) + 0.0
    dm = D_M(Ef)
    P("   %s: angekommen bei z = 10: ΔE/E = %.1e; z = 1090: 0 (vor T ≈ τ_D/40 ist f < e^{−10}); D_M(1090)/D_M,ΛCDM = %.3f"
      % (str(key), float(np.interp(1 / 11, AG, d_['R'] / (3 * H0 * Gyr) / EBG)), dm / dm_L))
    L_min = KREISPAAR_GRENZE * d_LSS() * dm / dm_L
    P("      Q-Fenster oben (Kreispaare, L ≥ 0,985·d_LSS mit H_eff): %.2f (ΛCDM-Abstand %.2f); unten 5." % (q_fenster(L_min=L_min)[1], q_fenster()[1]))
    # Wachstum mit H_eff (das Übergegangene ballt nicht: es ist Atmen, G)
    def wachs(Efun):
        def f(x, y):
            a = math.exp(x); h = 1e-4
            dl = (math.log(Efun(a * math.exp(h))) - math.log(Efun(a * math.exp(-h)))) / (2 * h)
            return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Efun(a)**2 * y[0]]
        s = solve_ivp(f, (math.log(AG[0]), 0.0), [AG[0], AG[0]], dense_output=True, rtol=1e-8, atol=1e-12)
        return lambda a: s.sol(math.log(a))
    gL, gE = wachs(E_L), wachs(Ef)
    fs = []
    for z in (0.0, 0.3, 0.6, 1.0, 1.5):
        yl, ye = gL(1 / (1 + z)), gE(1 / (1 + z))
        fs.append((z, ye[1] / yl[1]))                           # fσ₈ ∝ dD/dln a bei gleicher Anfangsamplitude
    P("      fσ₈/fσ₈,ΛCDM (gleiche frühe Amplitude): " + ", ".join("z = %.1f: %.3f" % x for x in fs))
P("   Λ = c²/2a₀ bleibt unberührt: der Übergang geht ins Atmen (G), nicht in die Schließlänge (I); a₀ = c²/2Λ ist seit z ≈ 10")
P("   fest, auch mit Übergang (τ_D und f(T) hängen nur an Λ).")
P()
P("   Aus dem Modell erzwungen (vor dem Ergebnis festgelegt):")
P("   (a) M_min: das Axiom gilt für jeden gebundenen Bereich, das Modell setzt keine untere Schranke ⇒ so klein wie Bindung")
P("       reicht (hier 10⁶, Rand des Gitters); in der geschachtelten Zählung hängt Ω_X(0) darunter kaum noch an M_min.")
P("   (b) geschachtelt ('oben'): beim Einfall löst sich die Naht des Kindes auf, sein Maß ist Teilmenge des Elternbereichs")
P("       (Projektion in ϖ, ein Zustandsraum) – je Baryon zählt der größte Bereich (de_fluchtmoden_hierarchie).")
P("   (c) θ̃ = 0: θ in der Rückkehrgleichung ist die Ausdehnung des Bereichs, dessen Maß zurückkehrt; er ist gebunden,")
P("       θ̃ = max(min θ, 0) = 0 (Gedächtnis). θ = 3H gehört dem Ganzen, das empfängt, nicht dem Weg; es ließe 98,5 % liegen.")
P("   (d) ab I_F: die Übergabe in die Form des Operators ist dort erzwungen (03); erst dann ist die Flucht inaktiv.")
d_, e_ = TAB[(6, "oben", "0", "I_F")]
P("   ⇒ (10⁶, oben, θ̃ = 0, I_F): Ω_X(0) = %.3f (%.0f %% von %.3f), Ω_X(z = 1) = %.3f, Gipfel z = %.2f%s, w(0) = %.3f,"
  % (e_['OX0'], 100 * e_['OX0'] / OM_X, OM_X, e_['OX1'], e_['zg'], " (Rand)" if e_['rand'] else "", e_['w_heute']))
P("     w₀ = %.2f, w_a = %.2f; H/H_ΛCDM [G] z = 0: %.3f, 1: %.3f, 2: %.3f, 3: %.3f." % (e_['w0'], e_['wa'], e_['hz'][0], e_['hz'][1], e_['hz'][2], e_['hz'][3]))

# ------------------------------------------------ 4 Blick
kopf("4  Beobachterfehler")
P("   X ist das Maß des Zustandsraums je Bereich (G des Bereichs), je Baryon gemittelt über die Zelle (G). T ist Rückkehrzeit in")
P("   der inneren Richtung (I), keine Laufzeit im Raum und keine Strecke. Angekommen ändert es die Atemrate der Zelle (G);")
P("   H_eff = H_b + ⅓R ist eine Rate des Ganzen, keine Dichte am Ort und kein Stoff, der ballt (Wachstum nur über H).")
P("   Zustandsabhängig: X, R und Ω_X hängen an Epoche, Massenfunktion und Übergabe; τ_D = c/a₀ ist fest (Λ).")

# ------------------------------------------------ 5 Von Hand
kopf("5  Von Hand")
P("   τ_D = c/a₀ = 2,998e8/1,1e-10 = 2,725e18 s = %.2f Gyr; Maximum τ_D/6 = %.2f Gyr." % (TD, TD / 6))
P("   angekommen nach 13,8 Gyr (Puls bei t = 0): erfc(√(86,36/55,2)) = erfc(1,2508) = %.4f (Kern-Inversion %.4f)."
  % (erfc(math.sqrt(TD / 55.2)), CDF(13.8)))
P("   f(13,8 Gyr) = √(86,36/4π)·13,8^{−1,5}·e^{−86,36/55,2} = 2,622·0,01951·0,2092 = %.5f /Gyr (Kern %.5f)."
  % (f_form(13.8), float(talbot(lambda p: Khat(p), 13.8))))

d_, e_ = TAB[(6, "oben", "0", "I_F")]
P("   Erzwungene Kombination grob: X(0) = %.2f, Quelle bei t' ≈ 1,8·(1–4) Gyr, Laufzeit T ≈ 7–12 Gyr, f ≈ 0,006–0,010/Gyr ⇒"
  % e_['X0'])
P("   R ≈ 3,35·0,008 = 0,027/Gyr, R/H₀ = 0,027/0,0689 = 0,39 (Code %.3f); ΔE = R/3H₀ = 0,13, Ω_X = (0,747 + 0,13)² − 0,747² = 0,21 (Code %.3f)."
  % (float(d_['R'][i0]) / (H0 * Gyr), e_['OX0']))
P("   Gipfel nahe heute nicht zufällig: τ_D/6 = 2π/(6H₀) = π/(3H₀) = %.1f Gyr bei a₀ = cH₀/2π (t₀ = %.1f Gyr)."
  % (math.pi / (3 * H0 * Gyr), T0))

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_fluchtmoden_verzoegert.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
