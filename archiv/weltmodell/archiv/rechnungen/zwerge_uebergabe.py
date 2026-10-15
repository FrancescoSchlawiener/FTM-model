# -*- coding: utf-8 -*-
"""
Zwerge vor der Übergabe: Liegen frisch geschlossene Systeme näher an Newton, weil die Mode ξ noch nicht übergeben ist?

03 (Gebunden, Naht): gebunden w = 0, Übergaberate λ = (1 − w)ω/N*, erzwungen bei der ersten Faltung I_F = ∫ω dt ≈ 1,1–1,25.
   Rest der Mode ξ_R(I) = e^{−I/N*} für I < I_F, 0 ab I_F;  übergebener Anteil f = 1 − ξ_R.
   a = f·a_L + ξ_R·a_N·(1 + κQ)
     κ = 0: Rest mit dem Gewicht w gezählt (03 Naht: a/a_N → 1 + wQ, gebunden w = 0) ⇒ Rest = Newton
     κ = 1: Rest wie Staub mit den Baryonen (kollaps_nach_umkehr g_hom) ⇒ Rest = (1 + Q)·a_N
   a_L = örtliche Regel mit Außenfeld (Monopol, wie ort_oder_system L).
1  σ_los(f) je System, benötigtes f, I = −N*·ln(1 − f), Zeit seit Schließung t = I/ω; ω² = g(r½)/r½
2  SPARC: I = ω·t an den äußersten Punkten (ω = V_obs/R)
3  Von Hand
Eingaben und Messwerte [G] wie ort_oder_system (aus glaettung_stetig / zwerge_aussenfeld).
"""
import io, sys, os, math
import numpy as np
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import A0, nu, kpc, Q_KOSM, Gyr
import sparc_flaechendichte as sf, sparc_universal as su

np.seterr(all="ignore")
GK = 4.30091e-6                                  # G in kpc (km/s)²/M☉
A0K = A0 * kpc / 1e6                             # a₀ in (km/s)²/kpc
KMSKPC = 1e3 / kpc * Gyr                         # 1 km/s/kpc in 1/Gyr
MU, WMU = np.polynomial.legendre.leggauss(400)
I_F = (1.10, 1.25)
NST = (0.3, 2.4)


def nux(x):
    return nu(np.maximum(np.asarray(x, float), 1e-300) * A0, A0)


def mono(gi, ge):
    g = np.atleast_1d(gi)[:, None]; m = MU[None, :]
    gm = np.sqrt(np.maximum(g**2 + ge**2 + 2 * g * ge * m, 0))
    return 0.5 * np.sum(WMU * nux(gm) * (g + ge * m), axis=1)


def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return math.sqrt(max(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r), 0) * A0K)


SYS = [  # Name, M [M☉], R_e [kpc], Wirt [M☉], Abstand [kpc], σ gemessen, (+, −), Alter der Sterne [G, Gedächtnis]
    ("Pal 14", 2e4, 0.0247, 6e10, 71.0, 0.38, (0.12, 0.12), "≈ 10–11 Gyr"),
    ("Crater II", 3.2e5, 1.066, 6e10, 117.5, 2.7, (0.3, 0.3), "≈ 10–12 Gyr"),
    ("Draco", 5.8e5, 0.221, 6e10, 76.0, 9.1, (1.2, 1.2), "≈ 10–12 Gyr"),
    ("Fornax", 4.0e7, 0.71, 6e10, 147.0, 11.7, (0.9, 0.9), "gemischt, bis 12 Gyr"),
    ("DF2", 2.0e8, 2.2, 1e11, 80.0, 8.5, (2.3, 3.1), "≈ 9 Gyr"),
    ("DF4", 1.5e8, 1.6, 1e11, 165.0, 4.2, (4.4, 2.2), "≈ 9 Gyr")]

out = []
P = out.append
P("0  Form aus 03: ξ_R = e^{−I/N*} (I < I_F), 0 ab I_F;  a = (1 − ξ_R)·a_L + ξ_R·a_N·(1 + κQ)")
P(f"   κ = 0: Rest mit w gezählt, gebunden w = 0 ⇒ Newton;  κ = 1: Rest wie Staub ⇒ (1 + Q)a_N, 1 + Q = {1 + Q_KOSM:.3f}")
P(f"   größter Rest vor der Faltung ξ_R(I_F⁻) = e^{{−I_F/N*}}: N* = 0,3: {math.exp(-I_F[0]/0.3):.3f}–{math.exp(-I_F[1]/0.3):.3f};"
  f"  N* = 2,4: {math.exp(-I_F[0]/2.4):.3f}–{math.exp(-I_F[1]/2.4):.3f}")
P("")
P("1  σ_los(R < R_e) [km/s], Plummer, isotrop; Eingaben [G]; ω = √(g(r½)/r½) in 1/Gyr (Newton | Regel); t_F = I_F/ω")
P(f"   {'System':10s} {'Newton':>7s} {'Staub':>7s} {'Regel':>7s}  {'gemessen':>14s}  {'f nötig κ=0':>11s} {'κ=1':>6s}   {'ω_N':>6s} {'ω_L':>6s}   {'t_F [Gyr] (ω_N, I_F 1,10–1,25)':>31s}  {'t für f (κ=0), N* 0,3 | 2,4 [Gyr]':>34s}  Sterne [G]")
RES = {}
for nm, M, Re, Mh, D, so, er, alt in SYS:
    b = Re; GM = GK * M / A0K
    r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, 6000)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5
    rho = (1 + (r / b)**2)**-2.5
    gN = GMenc / r**2
    ge = GK * Mh / D**2 / A0K
    gL = mono(gN, ge)

    def sig(f, kap):
        return jeans(r, rho, f * gL + (1 - f) * gN * (1 + kap * Q_KOSM), Re)
    sN, sS, sL = sig(0, 0), sig(0, 1), sig(1, 0)
    fn = {}
    for kap in (0, 1):
        lo, hi = sig(0, kap), sig(1, kap)
        if (lo - so) * (hi - so) < 0:
            fn[kap] = brentq(lambda f: sig(f, kap) - so, 0, 1, xtol=1e-6)
        else:
            fn[kap] = None
    rh = b / math.sqrt(2**(2 / 3) - 1)
    i = np.searchsorted(r, rh)
    wN = math.sqrt(gN[i] * A0K / rh) * KMSKPC
    wL = math.sqrt(gL[i] * A0K / rh) * KMSKPC
    tF = (I_F[0] / wN, I_F[1] / wN)

    def tfor(f, Ns):
        if f is None:
            return "–"
        I = -Ns * math.log(1 - f) if f < 1 else float("inf")
        if I >= I_F[0]:
            return f"I={I:.2f}≥I_F"
        wm = (1 - f) * wN + f * wL
        return f"{I / wm:.3f}"
    fs = lambda f: "  –   " if f is None else f"{f:6.3f}"
    RES[nm] = dict(sN=sN, sS=sS, sL=sL, f0=fn[0], f1=fn[1], wN=wN, wL=wL, tF=tF)
    P(f"   {nm:10s} {sN:7.2f} {sS:7.2f} {sL:7.2f}  {so:5.2f} +{er[0]}/−{er[1]}  {fs(fn[0]):>11s} {fs(fn[1]):>6s}   {wN:6.2f} {wL:6.2f}   "
      f"{tF[0]:13.3f}–{tF[1]:.3f}            {tfor(fn[0], 0.3):>14s} | {tfor(fn[0], 2.4):>14s}   {alt}")
P("   (f nötig: übergebener Anteil, der σ gemessen trifft; – = kein f ∈ [0, 1] trifft; Staub = f = 0 mit κ = 1)")
P("   In 2σ bei f = 0 (κ = 0, Newton): " + ", ".join(
    f"{k} {'ja' if abs(v['sN'] - SYS[j][5]) <= 2 * (SYS[j][6][0] if v['sN'] > SYS[j][5] else SYS[j][6][1]) else 'nein'}"
    for j, (k, v) in enumerate(RES.items())))
P("   In 2σ bei f = 1 (Regel):          " + ", ".join(
    f"{k} {'ja' if abs(v['sL'] - SYS[j][5]) <= 2 * (SYS[j][6][0] if v['sL'] > SYS[j][5] else SYS[j][6][1]) else 'nein'}"
    for j, (k, v) in enumerate(RES.items())))

# ------------------------------------------------------------------ 2 SPARC
tab = su.tab_voll()
Imin = []
for n, d in tab.items():
    if d['Q'] > 2 or d['inc'] < 30:
        continue
    a = np.loadtxt(os.path.join(sf.D, 'Rotmod_LTG', n + '_rotmod.dat'), comments='#')
    R, Vo = a[:, 0], a[:, 1]
    m = Vo > 0
    if m.sum() < 5:
        continue
    R, Vo = R[m], Vo[m]
    w = Vo[-1] / R[-1] * KMSKPC
    Imin.append((w, n, R[-1], Vo[-1]))
Imin.sort()
ws = np.array([x[0] for x in Imin])
P("")
P(f"2  SPARC (Q ≤ 2, i ≥ 30°): {len(Imin)} Galaxien, äußerster Punkt ω = V/R [1/Gyr]; I = ω·t")
P(f"   ω am äußersten Punkt: Minimum {ws[0]:.2f} ({Imin[0][1]}, R = {Imin[0][2]:.1f} kpc, V = {Imin[0][3]:.0f} km/s), 5 %-Quantil {np.quantile(ws, 0.05):.2f}, Median {np.median(ws):.2f}")
for t in (0.25, 1.0, 3.0, 10.0):
    P(f"   seit Schließung t = {t:5.2f} Gyr: I < I_F = 1,10 bei {np.sum(ws * t < 1.10)} von {len(ws)} Galaxien; Median I = {np.median(ws) * t:.1f} rad")
P(f"   größte Zeit, bei der noch ein äußerster SPARC-Punkt I < 1,25 hätte: {1.25 / ws[0]:.2f} Gyr")

# ------------------------------------------------------------------ 3 von Hand
P("")
P("3  Von Hand")
rh = 2.2 / math.sqrt(2**(2 / 3) - 1)
w2 = GK * 1e8 / rh**3
P(f"   DF2: r½ = 2,2/√(2^(2/3) − 1) = {rh:.3f} kpc;  ω² = G·(M/2)/r½³ = 4,30091·10⁻⁶·10⁸/{rh**3:.2f} = {w2:.2f} (km/s/kpc)²;"
  f"  ω = {math.sqrt(w2):.3f} km/s/kpc = {math.sqrt(w2) * KMSKPC:.3f}/Gyr (Plummer-Profil innen, Code: {RES['DF2']['wN']:.3f})")
P(f"   t_F = 1,10/{math.sqrt(w2) * KMSKPC:.3f} = {1.10 / (math.sqrt(w2) * KMSKPC):.3f} Gyr")
P(f"   DF2 Staub (κ = 1, f = 0): σ = √(1 + Q)·σ_N = {math.sqrt(1 + Q_KOSM):.3f}·{RES['DF2']['sN']:.2f} = {math.sqrt(1 + Q_KOSM) * RES['DF2']['sN']:.2f} km/s (Code {RES['DF2']['sS']:.2f})")

txt = "\n".join(out)
print(txt)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zwerge_uebergabe.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
