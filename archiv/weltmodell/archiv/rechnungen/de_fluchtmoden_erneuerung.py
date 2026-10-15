# -*- coding: utf-8 -*-
"""
Übergang des inaktiven Fluchtmaßes ins Atmen des Ganzen mit der ERNEUERUNG der Rückkehr (07: u(x,0) = g_N + u(x,ϖ_R),
ν = Σqⁿ) statt nur der ersten Rückkehr (de_fluchtmoden_verzoegert).
Blick: X je Baryon ist Maß des Zustandsraums je Bereich (G des Bereichs); der Weg läuft in I (Rückkehrzeit T, keine
Strecke); angekommen wirkt er aufs Atmen (G): H_eff = H_b + ⅓R(t). Keine Dichte am Ort.
1 Kerne aus rueckkehr(): erste Rückkehr f = L⁻¹[K̂], Erneuerungsdichte u = L⁻¹[K̂/(1 − K̂)] = Σ f^{*n}; Grenzverhalten,
  geschlossene Rückkehr (Wand bei Λ/ϖ_R) gegen Kac.
2 Grundfall mit u; Ω_X(z), Gipfel, w, w₀, w_a, H/H_ΛCDM, D_M(1090), Q-Fenster, fσ₈.
3 Tabelle: M_min, Beginn (Umkehr / I_F), Maß je Rücklauf (f einmal / u voll / u anteilig 1 − q), Verlust θ̃ = 0 / 3H.
4 Blick.  5 Von Hand.
Aus ftm_kern: rueckkehr, tau_D, r_stern, E_hintergrund, omega_vakuum, q_fenster, d_LSS, Konstanten.
Aus de_fluchtmoden: masse, XT, Dg, delta_ta, LM, LMg, SIG, EB, R_ta, NK.
Nicht im Kern (lokal, wie de_fluchtmoden_verzoegert): Talbot-Inversion, Faltung, Wachstum mit H_eff, t_F/t_ta = 1,8.
'anteilig': je Rücklauf geht 1 − q ins Ganze, q = e^{−r*/R_ta} (Regel an der Naht des Bereichs), der Rest erneuert;
  gerechnet als Quelle gewichtet mit 1 − q(M, a) und Kern u (Fehler O(1 − q) < 1 %).
"""
import math, os, sys, io, contextlib
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import erfc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import (c, M_sun, Mpc, Gyr, OMEGA_M, OMEGA_DE, H0_si, E_hintergrund, omega_vakuum, rueckkehr, tau_D,
                      schliesslaenge, q_fenster, d_LSS, r_stern, KREISPAAR_GRENZE)
with contextlib.redirect_stdout(io.StringIO()):
    import de_fluchtmoden as dfm

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

H0 = H0_si()
TD = tau_D() / Gyr
OM_X = OMEGA_DE - omega_vakuum()
E_L = lambda a: E_hintergrund(1 / a - 1)
EB = dfm.EB
T_F = 1.8

def t_von_a(a):
    return quad(lambda x: 1 / (x * E_L(x)), 1e-9, a, limit=200)[0] / (H0 * Gyr)
T0 = t_von_a(1.0)

def talbot(Fhat, t, M=48):
    r = 2 * M / (5 * t)
    k = np.arange(1, M)
    th = k * math.pi / M
    cot = 1 / np.tan(th)
    d = r * th * (cot + 1j)
    g = 1 + 1j * th * (1 + cot**2) - 1j * cot
    s = 0.5 * math.exp(r * t) * Fhat(r) + np.sum((np.exp(t * d) * Fhat(d) * g).real)
    return float((r / M * s).real)

def Khat(p):                                   # p in 1/Gyr, Kern rueckkehr(p, Γ = 0, θ = 0)
    return rueckkehr(np.asarray(p, dtype=complex) / Gyr, 0.0, 0.0)
Uhat = lambda p: Khat(p) / (1 - Khat(p))
f_form = lambda T, tau=TD: math.sqrt(tau / (4 * math.pi)) * T**-1.5 * math.exp(-tau / (4 * T)) if T > 0 else 0.0
u_reihe = lambda T: sum(f_form(T, n * n * TD) for n in range(1, 60))      # Σ f^{*n}: K̂ⁿ = e^{−√(p n²τ_D)}

# ============================================================== 1 Kerne
kopf("1  Welcher Kern? Erste Rückkehr f gegen Erneuerungsdichte u = Σf^{*n} (beide aus rueckkehr(), θ̃ = 0)")
P("   τ_D = %.2f Gyr.  K̂ⁿ = e^{−n√(pτ_D)} ⇒ f^{*n} = f mit τ_D → n²τ_D (n-ter Rücklauf braucht n²-fache Zeit)." % TD)
P("   T [Gyr]:        " + "  ".join("%8.2f" % t for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   f (Talbot)     " + "  ".join("%8.5f" % talbot(Khat, t) for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   u (Talbot)     " + "  ".join("%8.5f" % talbot(Uhat, t) for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   u = Σ f_{n²τ}  " + "  ".join("%8.5f" % u_reihe(t) for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   1/√(πτ_D T)    " + "  ".join("%8.5f" % (1 / math.sqrt(math.pi * TD * t)) for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   1/τ_D (Kac)    " + "  ".join("%8.5f" % (1 / TD) for t in (5, 10, 13.8, 30, 86.36, 300, 1000, 1e4)))
P("   ∫₀ᵗu = Σ erfc(n√(τ_D/4t)):  t = 13,8: %.5f (f allein %.5f);  t = 86: %.4f;  t = 1000: %.3f (2√(t/πτ_D) − ½ = %.3f)"
  % (sum(erfc(n * math.sqrt(TD / 55.2)) for n in range(1, 50)), erfc(math.sqrt(TD / 55.2)),
     sum(erfc(n * math.sqrt(TD / (4 * 86.36))) for n in range(1, 80)),
     sum(erfc(n * math.sqrt(TD / 4000)) for n in range(1, 200)), 2 * math.sqrt(1000 / (math.pi * TD)) - 0.5))
P("   Offene Halbgerade (Kern des Modells): ⟨T⟩ = ∞, u → 1/√(πτ_D t) → 0, nicht 1/τ_D; ∫u wächst ohne Grenze (Σqⁿ, q → 1).")
P("   Geschlossen (Wand bei L = Λ/ϖ_R, Start ϖ_R): K̂_L = cosh((L−1)y)/cosh(Ly), y = √(pτ_D); u → 1/⟨T⟩, ⟨T⟩ = (2L − 1)τ_D/2:")
for L in (1.5, 3.0, 10.0):
    KL = lambda p, L=L: np.cosh((L - 1) * np.sqrt(np.asarray(p, dtype=complex) * TD)) / np.cosh(L * np.sqrt(np.asarray(p, dtype=complex) * TD))
    UL = lambda p, KL=KL: KL(p) / (1 - KL(p))
    P("     L = %4.1f: u(13,8) = %.5f, u(86) = %.5f, u(10·L²τ_D) = %.5f, 1/⟨T⟩ = %.5f /Gyr"
      % (L, talbot(UL, 13.8), talbot(UL, 86.36), talbot(UL, 10 * L * L * TD), 2 / ((2 * L - 1) * TD)))
Mg = 1e12 * M_sun / dfm.NK
P("   ⟨T⟩ = τ_D (Kac, 07) nur bei L = 3/2; real L = Λ/r*: 10¹² M☉ %.2e, dort wirkt die Wand erst ab T ~ L²τ_D. Bis t₀ gilt u ≈ f."
  % (schliesslaenge() / r_stern(Mg)))

# ============================================================== X(a), Quellen
def _wachstum_bis(a_end=3.5):
    def fn(x, y):
        a = math.exp(x); h = 1e-5
        dlnE = (math.log(E_L(a * math.exp(h))) - math.log(E_L(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dlnE) * y[1] + 1.5 * OMEGA_M * a**-3 / E_L(a)**2 * y[0]]
    s = solve_ivp(fn, (math.log(1e-3), math.log(a_end)), [1e-3, 1e-3], dense_output=True, rtol=1e-9, atol=1e-12)
    D1 = s.sol(0.0)[0]
    return lambda a: s.sol(math.log(a))[0] / D1
_Dg_alt = dfm.Dg
dfm.Dg = _wachstum_bis()
assert abs(dfm.Dg(0.5) / _Dg_alt(0.5) - 1) < 1e-6
AA = np.exp(np.linspace(math.log(1 / 31), math.log(3.0), 56))
PQ = {}
for a in AA:
    if a not in dfm.XT:
        rows = []
        for l in dfm.LMg:
            m = dfm.masse(10**l * M_sun / dfm.NK, a, n=1500)
            rows.append((m['muI'] / m['muF'], math.log(m['muR'] / m['muN']), m['muI'] / m['muR']))
        dfm.XT[a] = np.array(rows)
    PQ[a] = np.array([1 - math.exp(-r_stern(10**l * M_sun / dfm.NK) / dfm.R_ta(10**l * M_sun / dfm.NK, a)) for l in dfm.LMg])
LN10 = math.log(10)

def X_oben(a, lm_min, anteilig=False):
    Fg = erfc(dfm.delta_ta(a) / (math.sqrt(2) * dfm.SIG * dfm.Dg(a)))
    L = np.interp(dfm.LM, dfm.LMg, dfm.XT[a][:, 1])
    if anteilig:
        L = L * np.interp(dfm.LM, dfm.LMg, PQ[a])
    msk = dfm.LM >= lm_min
    dF = -np.gradient(Fg, dfm.LM * LN10)
    return float(np.trapezoid((dF * L)[msk], dfm.LM[msk] * LN10))

TA = np.array([t_von_a(a) for a in AA])
NT = 3000
TG = np.linspace(TA[0], TA[-1], NT); dt = TG[1] - TG[0]
AG = np.exp(np.interp(TG, TA, np.log(AA)))
HG = H0 * Gyr * np.array([E_L(a) for a in AG])
EBG = np.array([EB(a) for a in AG])
ELG = np.array([E_L(a) for a in AG])
i0 = int(np.argmin(abs(TG - T0)))
TK = np.arange(NT) * dt
KERN = {"f": np.array([f_form(t) for t in TK]), "u": np.array([u_reihe(t) for t in TK])}
KERN["u"][0] = KERN["f"][0] = 0.0

def lauf(lm_min, kern, verlust, ab):
    anteilig = kern == "u·(1−q)"
    Ki = KERN["f"] if kern == "f" else KERN["u"]
    X = np.array([X_oben(a, lm_min, anteilig) for a in AA])
    sp = CubicSpline(TA, X)
    sf = T_F if ab == "I_F" else 1.0
    tq = TG / sf
    S = sp(np.clip(tq, TA[0], None), 1) * (tq >= TA[0]) / sf
    X0 = X[0]; t_puls = TA[0] * sf
    Kf = lambda T: np.interp(T, TK, Ki, left=0.0, right=0.0)
    a_puls = np.interp(t_puls, TG, AG)
    lna3 = 3 * np.log(AG)
    R = np.zeros(NT)
    for n in range(NT):
        w = Ki[n::-1] * S[:n + 1]
        if verlust == "3H":
            w = w * np.exp(lna3[:n + 1] - lna3[n])
        R[n] = np.trapezoid(w, dx=dt)
        tp = TG[n] - t_puls
        if tp > 0:
            R[n] += X0 * Kf(tp) * (1.0 if verlust == "0" else (a_puls / AG[n])**3)
    Eeff = EBG + R / (3 * H0 * Gyr)
    return dict(X=X, S=S, R=R, Xarr=cumulative_trapezoid(R, TG, initial=0.0), Eeff=Eeff, OX=Eeff**2 - EBG**2, rate=R / HG)

def auswertung(d):
    a = AG; OX = d['OX']
    g = int(np.argmax(d['rate']))
    lnr = np.log(np.clip(OX, 1e-300, None)); lna = np.log(a)
    wz = -1 - np.gradient(lnr, lna) / 3
    m2 = (a >= 1 / 3) & (a <= 1.0)
    A_ = np.vstack([np.ones(m2.sum()), 1 - a[m2]]).T
    w0, wa = np.linalg.lstsq(A_, wz[m2], rcond=None)[0]
    oz = {z: float(np.interp(1 / (1 + z), a, OX)) for z in (0.5, 1, 2)}
    hz = {z: float(np.interp(1 / (1 + z), a, d['Eeff'])) / E_L(1 / (1 + z)) for z in (0, 0.5, 1, 2, 3)}
    return dict(OX0=float(OX[i0]), oz=oz, zg=1 / a[g] - 1, rand=g >= len(a) - 2, w0=w0, wa=wa, wh=float(wz[i0]), hz=hz,
                Xarr0=float(d['Xarr'][i0]), X0=float(np.interp(1.0, AA, d['X'])))

# ============================================================== 2 Grundfall
kopf("2  Erneuerung u (M_min = 10⁶, oberste Ebene, θ̃ = 0, ab I_F) gegen erste Rückkehr f")
ERG = {}
for kern in ("f", "u"):
    d = lauf(6, kern, "0", "I_F"); e = auswertung(d); ERG[kern] = (d, e)
    P("   %s: Ω_X(0) = %.4f (%.0f %% von %.3f); Ω_X(z = 0,5/1/2) = %.4f/%.4f/%.4f; Gipfel z = %.2f%s; w(0) = %.3f; w₀ = %.3f, w_a = %.3f"
      % (kern, e['OX0'], 100 * e['OX0'] / OM_X, OM_X, e['oz'][0.5], e['oz'][1], e['oz'][2], e['zg'], "*" if e['rand'] else "",
         e['wh'], e['w0'], e['wa']))
    P("      H/H_ΛCDM [G] z = 0/0,5/1/2/3: %s; X(0) = %.3f, bis heute angekommen %.4f (%.1f %%)"
      % ("/".join("%.3f" % e['hz'][z] for z in (0, 0.5, 1, 2, 3)), e['X0'], e['Xarr0'], 100 * e['Xarr0'] / e['X0']))
du, df_ = ERG["u"][0], ERG["f"][0]
P("   R_u/R_f − 1: z = 2: %.4f, 1: %.4f, 0: %.4f, a = 2: %.4f, a = 3: %.4f"
  % tuple(float(np.interp(a, AG, du['R'])) / float(np.interp(a, AG, df_['R'])) - 1 for a in (1 / 3, 0.5, 1.0, 2.0, 2.99)))

def D_M(Efun, z=1090.0):
    zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, z, 4001)[1:]])
    return np.trapezoid(1 / np.array([Efun(1 / (1 + x)) for x in zs]), zs)
dm_L = D_M(E_L)
def wachs(Efun):
    def fn(x, y):
        a = math.exp(x); h = 1e-4
        dl = (math.log(Efun(a * math.exp(h))) - math.log(Efun(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Efun(a)**2 * y[0]]
    s = solve_ivp(fn, (math.log(AG[0]), 0.0), [AG[0], AG[0]], dense_output=True, rtol=1e-8, atol=1e-12)
    return lambda a: s.sol(math.log(a))
gL = wachs(E_L)
def zusammenhang(d, name):
    Ef = lambda a: float(np.interp(a, AG, d['Eeff'])) if a >= AG[0] else EB(a)
    dm = D_M(Ef)
    L_min = KREISPAAR_GRENZE * d_LSS() * dm / dm_L
    gE = wachs(Ef)
    fs = ", ".join("z = %.1f: %.3f" % (z, gE(1 / (1 + z))[1] / gL(1 / (1 + z))[1]) for z in (0.0, 0.3, 0.6, 1.0, 1.5))
    P("   %s: D_M(1090)/D_M,ΛCDM = %.3f; Q-Fenster oben %.2f (ΛCDM-Abstand %.2f), unten 5; ΔE/E(z = 10) = %.1e"
      % (name, dm / dm_L, q_fenster(L_min=L_min)[1], q_fenster()[1], float(np.interp(1 / 11, AG, d['R'] / (3 * H0 * Gyr) / EBG))))
    P("      fσ₈/fσ₈,ΛCDM (gleiche frühe Amplitude): " + fs)
for kern in ("f", "u"):
    zusammenhang(ERG[kern][0], "Kern " + kern)
# Kac-Grenze als Probe der Annahme u → 1/τ_D (nicht aus dem Kern): R = X_übergeben/τ_D ab I_F
d = ERG["u"][0]
sp = CubicSpline(TA, d['X']); Xg = sp(np.clip(TG / T_F, TA[0], None)) * (TG / T_F >= TA[0])
Ek = EBG + Xg / TD / (3 * H0 * Gyr)
P("   Probe (nicht aus dem Kern): wäre u = 1/τ_D ab Übergabe, R = X/τ_D: Ω_X(0) = %.3f, Ω_X(z = 1) = %.3f, E/E_ΛCDM(0) = %.3f"
  % (float(Ek[i0]**2 - EBG[i0]**2), float(np.interp(0.5, AG, Ek**2 - EBG**2)), float(Ek[i0])))

# ============================================================== 3 Tabelle
kopf("3  Abhängigkeiten (oberste Ebene): M_min, Beginn, Maß je Rücklauf, Verlust")
lm_, a_ = dfm.LMg, 1.0
P("   1 − q an der Naht, q = e^{−r*/R_ta}: " + ", ".join("10^%d: a=0,25 %.4f, a=1 %.4f" % (l, float(np.interp(l, lm_, PQ[AA[np.argmin(abs(AA - 0.25))]])),
                                                              float(np.interp(l, lm_, PQ[AA[np.argmin(abs(AA - 1.0))]]))) for l in (6, 8, 10, 12)))
P("   M_min  ab      Kern     θ̃    Ω_X(0)  /0,442  Ω_X z=0,5   z=1      z=2    z_Gipfel  w(0)    w₀      w_a     H/H_Λ z=0    1      2      3    angek.")
TAB = {}
for lm in (6, 8, 10, 12):
    for ab in ("Umkehr", "I_F"):
        for kern in ("f", "u", "u·(1−q)"):
            for verl in ("0", "3H"):
                d_ = lauf(lm, kern, verl, ab); e_ = auswertung(d_); TAB[(lm, ab, kern, verl)] = e_
                P("   10^%-3d %-6s  %-8s %-3s  %7.4f  %6.3f   %7.4f  %7.4f  %7.4f   %6.2f%s  %6.3f  %6.3f  %7.3f   %6.3f  %6.3f  %6.3f  %6.3f  %.4f"
                  % (lm, ab, kern, "0" if verl == "0" else "3H", e_['OX0'], e_['OX0'] / OM_X, e_['oz'][0.5], e_['oz'][1], e_['oz'][2],
                     e_['zg'], "*" if e_['rand'] else " ", e_['wh'], e_['w0'], e_['wa'], e_['hz'][0], e_['hz'][1], e_['hz'][2], e_['hz'][3], e_['Xarr0']))
P("   f = Maß einmal bei erster Rückkehr; u = jeder Rücklauf gibt X voll ab (Maß nicht erhalten); u·(1−q) = je Rücklauf Anteil 1 − q,")
P("   Rest erneuert (Maß erhalten, ∫₀ᵗ ≤ ∫₀ᵗf). * Gipfel am Rand (a = 3). angek. = bis heute übergeben je Baryon.")

# ============================================================== 4 Blick
kopf("4  Beobachterfehler")
P("   Rückkehr und Erneuerung laufen in I (Rückkehrzeit T, keine Strecke, kein Weg im Raum); X ist Maß des Zustandsraums je")
P("   Bereich (G des Bereichs); angekommen ändert es die Atemrate der Zelle (G), H_eff = H_b + ⅓R, keine Dichte am Ort, kein")
P("   Stoff, der ballt. Die Erneuerung zählt Rückläufe als Ereignisse (Regel), nicht als Menge: 'voll je Rücklauf' läse die")
P("   Zählung der Regel als Maßübertrag und erzeugt Maß (∫u → ∞). τ_D und Λ fest; X, R, Ω_X zustandsabhängig.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand")
P("   u(13,8) = f_τ + f_4τ: 2,622·0,01951·e^{−1,5645} = 0,01070; √(345,4/4π) = 5,243, e^{−6,258} = 1,914e-3 ⇒ 5,243·0,01951·1,914e-3")
P("   = 1,96e-4; Summe 0,01090 /Gyr (Code u(13,8) = %.5f). ∫₀^{13,8}u = erfc(1,2508) + erfc(2,5016) = 0,07694 + 0,00041 = 0,0774" % u_reihe(13.8))
P("   (Code %.4f). Große t: û = 1/(e^{√(pτ)} − 1) ≈ 1/√(pτ) − ½ ⇒ ∫₀ᵗu ≈ 2√(t/πτ_D) − ½, bei t = 1000 Gyr: 2·1,919 − 0,5 = 3,34;"
  % sum(erfc(n * math.sqrt(TD / 55.2)) for n in range(1, 50)))
P("   nicht t/τ_D = 11,6 (das gälte nur geschlossen, ⟨T⟩ endlich). Heute: u/f − 1 = 1,8 %, Ω_X ändert sich um < 1 %.")
eu = ERG["u"][1]
P("   Ω_X(0) mit u grob: R ≈ X·ū, ΔE = R/3H₀ ⇒ Code R(0)/H₀ = %.3f, ΔE = %.3f, Ω_X = (0,747 + ΔE)² − 0,747² = %.3f (Code %.3f)."
  % (float(du['R'][i0]) / (H0 * Gyr), float(du['R'][i0]) / (3 * H0 * Gyr), (EBG[i0] + du['R'][i0] / (3 * H0 * Gyr))**2 - EBG[i0]**2, eu['OX0']))

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_fluchtmoden_erneuerung.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
