# -*- coding: utf-8 -*-
"""
Bistabilität, Hysterese, Torus-Geburt  (Forschungsnotiz 31)
===========================================================

  1  Torus-Geburt im Labor: VO2-Element auf einer Mikrobrücke im Magnetfeld
     (Spannung, Temperatur, Auslenkung, Geschwindigkeit; Joule, Kühlung,
     Metall-Isolator-Übergang, Lorentzkraft, Bewegungs-EMK, Übergangsdehnung)
     1a Parameter, Fixpunkte, Hopf-Punkt (Eigenwerte)
     1b Lyapunov-Spektren entlang der Vorspannung, Magnetfeld-Pfad, steife Brücke
     1c Floquet-Rate der mechanischen Mode auf dem Grenzzyklus (Neimark-Sacker)
  2  Kippstelle = Divergenz der Rückkopplungssumme
     2a allgemeine Identität dT*/dV = (dP/dV/G)·Σ q^n, Falte bei q = 1
     2b Spitze: Breite ∝ ε^{3/2}, Sprung ∝ ε^{1/2} (VO2 mit Lastwiderstand)
     2c kritische Verlangsamung λ = G(1-q)/C ∝ (V_f - V)^{1/2}
     2d Normalform x^3 - a x - b
     2e Galaxien: q = exp(-sqrt(y)); Relaxationszeit, Empfindlichkeit, Milchstraße
  3  Preisach gegen Gedächtniskern
     3a Rückkehrpunkt-Gedächtnis, Auslöschung, Extremwertstapel
     3b Ratenunabhängigkeit gegen Kern
     3c thermisch aktivierte Hysteronen: Übergang Preisach -> Prony-Kern
     3d w(θ) als Preisach-Operator: Kugelkollaps, Größe des Effekts
  4  Thermoelektrischer Transport im Magnetfeld (Bismut, Zweibandmodell)
     Onsager-Casimir, Zerlegung symmetrisch/antisymmetrisch, Entropieproduktion
  5  Anwendungen
     5a Curie-Punkt als stabile Randlage, Temperaturschwankung
     5b Kramers-Verweilzeiten und Kac-Budget im bistabilen VO2-Element

Nur numpy. Ausgabe zusätzlich nach bistabilitaet.txt (gleicher Ordner).
"""
import io
import math
import os
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

_buf = io.StringIO()


def out(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    _buf.write(s + "\n")


def head(t):
    out("")
    out("=" * 78)
    out(t)
    out("=" * 78)


def sub(t):
    out("")
    out("--- " + t)


T_START = time.time()

# =============================================================================
# VO2-Element (gemeinsam für Teil 1, 2, 5b)
# Einheiten: µs, K, V, kΩ, mA, mW, nF, nJ/K, µm
# =============================================================================
T0 = 300.0      # Substrat [K]
TC = 340.0      # Übergangstemperatur VO2 [K] (Literatur ~341 K)
WTR = 1.5       # Breite des Übergangs [K] (Annahme, glatter Übergang ohne eigene Hysterese)
RI = 20.0       # Widerstand isolierend [kΩ] (Annahme, typisches Mikroelement)
RM = 0.4        # Widerstand metallisch [kΩ] (Annahme)
CTH = 5e-3      # Wärmekapazität [nJ/K] = 5 pJ/K (Annahme)
GTH = 5e-3      # Wärmeleitwert [mW/K] = 5 µW/K (Annahme) -> τ_th = 1 µs
LNRR = math.log(RI / RM)


def sgm(T, sh=0.0):
    return 1.0 / (1.0 + np.exp((T - TC - sh) / WTR))


def Rvo2(T, sh=0.0):
    return np.exp(math.log(RM) + LNRR * sgm(T, sh))


# =============================================================================
head("TEIL 1  Torus-Geburt im Labor: VO2-Element auf einer Mikrobrücke im Magnetfeld")
# =============================================================================
#  C dV/dt   = (Vs - V)/RL - I,          I = (V - Γ v)/R(T, x)
#  Cth dT/dt = I^2 R - G (T - T0)
#  dx/dt = v
#  dv/dt = -(ω/Q) v - ω^2 y (1 + (y/xD)^2),  y = x - x_eq,
#  x_eq  = D_tr·(1 - s(T)) + λ_B·B·I     (Übergangsdehnung + Lorentzkraft)
#  R(T,x): Übergangstemperatur durch Dehnung verschoben, TC -> TC + κ x
RL1 = 10.0     # Lastwiderstand [kΩ]
CEL = 0.02     # Parallelkapazität [nF] = 20 pF


def f1(Z, p):
    V, T, x, v = Z
    Vs, kap, fm, Q, Dtr, lamB, Gam, xD = p
    Rv = Rvo2(T, kap * x)
    I = (V - Gam * v) / Rv
    om = 2 * np.pi * fm
    y = x - (Dtr * (1 - sgm(T)) + lamB * I)
    return np.array([((Vs - V) / RL1 - I) / CEL,
                     (I * I * Rv - GTH * (T - T0)) / CTH,
                     v,
                     -om / Q * v - om * om * y * (1 + (y / xD) ** 2)])


def fp1(p):
    Vs = p[0]
    Ts = np.linspace(300.01, 370, 70001)
    x = 0.0
    for _ in range(100):
        Rv = Rvo2(Ts, p[1] * x)
        V = Vs * Rv / (Rv + RL1)
        g = V * V / Rv - GTH * (Ts - T0)
        idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        if len(idx) != 1:
            return None
        k = idx[0]
        T = Ts[k] - g[k] * (Ts[k + 1] - Ts[k]) / (g[k + 1] - g[k])
        Rv = Rvo2(T, p[1] * x)
        V = Vs * Rv / (Rv + RL1)
        I = V / Rv
        xn = p[4] * (1 - sgm(T)) + p[5] * I
        if abs(xn - x) < 1e-13:
            break
        x = xn
    return np.array([V, T, x, 0.0])


def jac1(z, p, h=1e-6):
    J = np.zeros((4, 4))
    for i in range(4):
        dz = np.zeros(4)
        dz[i] = h * max(1, abs(z[i]))
        J[:, i] = (f1(z + dz, p) - f1(z - dz, p)) / (2 * dz[i])
    return J


# Grundparameter (alle Annahmen, siehe Notiz)
KAP = -0.1      # Verschiebung der Übergangstemperatur pro µm Auslenkung [K/µm]
FM = 0.5        # Eigenfrequenz der Brücke [MHz]
QM = 1000.0     # mechanische Güte (Vakuum)
DTR = 0.2       # Auslenkung durch die Übergangsdehnung über den ganzen Übergang [µm]
LAMB1 = 0.063   # Lorentz-Auslenkung pro mA und Tesla [µm/(mA T)]: Γ/k, ℓ = 100 µm, k = 1,58 N/m
GAM1 = 1e-4     # Γ = Bℓ pro Tesla [V/(m/s)]
XD = 0.3        # Duffing-Länge (geometrische Nichtlinearität) [µm]


def par1(Vs, B=1.0, kap=KAP, fm=FM, xD=XD):
    return (Vs, kap, fm, QM, DTR, LAMB1 * B, GAM1 * B, xD)


sub("1a Parameter und lineare Stabilität des Fixpunkts")
out("VO2: T_C = %.0f K, Breite %.1f K, R_isol = %.0f kΩ, R_met = %.1f kΩ, C_th = %.0f pJ/K, G = %.0f µW/K, τ_th = %.1f µs"
    % (TC, WTR, RI, RM, CTH * 1e3, GTH * 1e3, CTH / GTH))
out("Kreis: R_L = %.0f kΩ, C = %.0f pF;  Brücke: f_m = %.2f MHz, Q = %.0f, D_tr = %.2f µm, κ = %.2f K/µm, x_D = %.2f µm"
    % (RL1, CEL * 1e3, FM, QM, DTR, KAP, XD))
out("Lorentz: λ_B = %.3f µm/(mA·T), Γ = %.0e V·s/m pro T" % (LAMB1, GAM1))
out("Vs [V] |  T* [K] | x* [µm] | elektrothermisches Paar        | mechanisches Paar")
for Vs in [2.6, 3.0, 3.4, 3.5, 3.54, 3.55, 3.6, 4.0, 5.0]:
    p = par1(Vs)
    z = fp1(p)
    if z is None:
        out("%5.2f  | mehrere Fixpunkte" % Vs)
        continue
    ev = np.linalg.eigvals(jac1(z, p))
    ev = ev[ev.imag > 0]
    ev = ev[np.argsort(-np.abs(ev.imag))]
    if len(ev) < 2:
        out("%5.2f  | %7.2f | %.4f  | (reelle Eigenwerte)" % (Vs, z[1], z[2]))
        continue
    out("%5.2f  | %7.2f | %.4f  | Re %+8.4f /µs, f = %.3f MHz | Re %+8.5f /µs, f = %.4f MHz"
        % (Vs, z[1], z[2], ev[0].real, ev[0].imag / 2 / np.pi, ev[1].real, ev[1].imag / 2 / np.pi))


def lead_re(Vs):
    p = par1(Vs)
    z = fp1(p)
    ev = np.linalg.eigvals(jac1(z, p))
    ev = ev[ev.imag > 0]
    ev = ev[np.argsort(-np.abs(ev.imag))]
    return ev[0]


a_, b_ = 3.50, 3.60
for _ in range(40):
    m = 0.5 * (a_ + b_)
    if lead_re(m).real > 0:
        a_ = m
    else:
        b_ = m
VH = 0.5 * (a_ + b_)
evH = lead_re(VH)
out("Hopf-Punkt (elektrothermisches Paar, Re λ = 0): V_H = %.4f V, f_H = %.3f MHz" % (VH, evH.imag / 2 / np.pi))
zH = fp1(par1(VH))
out("   Arbeitspunkt am Hopf-Punkt: T* = %.2f K, R = %.3f kΩ, I = %.3f mA" %
    (zH[1], Rvo2(zH[1], KAP * zH[2]), zH[0] / Rvo2(zH[1], KAP * zH[2])))
# bistabiles Fenster des Kreises ohne Kapazität bei R_L = 10 kΩ (für die Einordnung)
Tg1 = np.linspace(300.5, 375, 200001)
Vb1 = np.sqrt(GTH * (Tg1 - T0)) * (Rvo2(Tg1) + RL1) / np.sqrt(Rvo2(Tg1))
dV1 = np.diff(Vb1)
ix = np.where(np.sign(dV1[:-1]) != np.sign(dV1[1:]))[0]
if len(ix) == 2:
    out("stationär bistabiles Fenster bei R_L = 10 kΩ (ohne Brücke): %.4f V < Vs < %.4f V"
        % (min(Vb1[ix + 1]), max(Vb1[ix + 1])))


def rk4_step(F, Z, p, dt):
    k1 = F(Z, p)
    k2 = F(Z + 0.5 * dt * k1, p)
    k3 = F(Z + 0.5 * dt * k2, p)
    k4 = F(Z + dt * k3, p)
    return Z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def lyap_batch(plist, dt=0.01, Ttr=2500.0, Tm=3000.0, ren=10, d0=1e-7):
    """Benettin: je Parametersatz eine Basisbahn + 4 Nachbarn, QR alle ren Schritte."""
    N = len(plist)
    M = 5
    P = tuple(np.repeat(np.array([pp[k] for pp in plist], dtype=float), M) for k in range(8))
    Z = np.zeros((4, N * M))
    for j, pp in enumerate(plist):
        z = fp1(pp)
        if z is None:
            z = np.array([1.5, 339.0, 0.1, 0.0])
        Z[:, j * M:(j + 1) * M] = (z + np.array([0.01, 0.05, 0.05, 0.0]))[:, None]
    for _ in range(int(Ttr / dt)):
        Z = rk4_step(f1, Z, P, dt)
    Zr = Z.reshape(4, N, M)
    base = Zr[:, :, 0].copy()
    for j in range(4):
        Zr[:, :, j + 1] = base
        Zr[j, :, j + 1] += d0
    Z = Zr.reshape(4, N * M)
    S = np.zeros((N, 4))
    rec = []
    nst = int(Tm / dt)
    for i in range(nst):
        Z = rk4_step(f1, Z, P, dt)
        if i % 5 == 0:
            rec.append(Z.reshape(4, N, M)[:, :, 0].copy())
        if (i + 1) % ren == 0:
            Zr = Z.reshape(4, N, M)
            base = Zr[:, :, 0]
            D = np.stack([(Zr[:, :, j + 1] - base).T for j in range(4)], axis=2) / d0
            Qm, Rr = np.linalg.qr(D)
            S += np.log(np.abs(np.diagonal(Rr, axis1=1, axis2=2)))
            for j in range(4):
                Zr[:, :, j + 1] = base + d0 * Qm[:, :, j].T
            Z = Zr.reshape(4, N * M)
    return S / Tm, np.array(rec), 5 * dt


def peak(s, dts):
    s = s - s.mean()
    S = np.abs(np.fft.rfft(s * np.hanning(len(s))))
    fr = np.fft.rfftfreq(len(s), dts)
    k = np.argmax(S[1:]) + 1
    # Parabel-Interpolation
    if 1 <= k < len(S) - 1:
        a, b, c = np.log(S[k - 1] + 1e-300), np.log(S[k] + 1e-300), np.log(S[k + 1] + 1e-300)
        d = 0.5 * (a - c) / (a - 2 * b + c)
    else:
        d = 0.0
    return fr[k] + d * (fr[1] - fr[0])


def classify(L, ampV, ampX, trend, fx, fV):
    if L[0] > 0.05:
        return "Chaos"
    if ampV < 1e-3:
        return "Fixpunkt" if (L[0] < 0 and trend < 0.97) else "Grenzzyklus (mech.)"
    r = fx / fV
    for qd in (1, 2, 3, 4, 5):
        pn = round(r * qd)
        if pn >= 1 and abs(r - pn / qd) < 2e-3:
            return "eingerastet %d:%d" % (pn, qd)
    if L[1] > -3e-3 and trend > 0.97:
        return "Torus T²"
    return "Grenzzyklus"


SKIP1 = bool(os.environ.get("SKIP1"))
if not SKIP1:
  sub("1b Lyapunov-Spektren (Benettin, RK4, dt = 0,01 µs, Einschwingen 2500 µs, Mittelung 3000 µs)")
  path = [3.60, 3.47, 3.46, 3.44, 3.40, 3.20, 2.95, 2.87, 2.80, 2.70]
  Bpath = [0.0, 2.0, 4.0]
  stiff = [2.90, 2.80]
  plist = [par1(v) for v in path] + [par1(3.47, B=b) for b in Bpath] + \
          [par1(v, kap=-0.05, fm=2.0, xD=0.1) for v in stiff]
  labels = ["Vs=%.2f B=1" % v for v in path] + ["Vs=3.47 B=%g" % b for b in Bpath] + \
           ["steif Vs=%.2f" % v for v in stiff]
  t0 = time.time()
  L, rec, dts = lyap_batch(plist)
  out("(Rechenzeit %.0f s)" % (time.time() - t0))
  out("Lauf            | λ1       λ2       λ3       λ4   [1/µs] | Δx [µm] ΔV [V] | f_x [MHz] f_V [MHz] f_x/f_V | Zustand")
  res1 = {}
  for j, lab in enumerate(labels):
      x = rec[:, 2, j]
      V = rec[:, 0, j]
      h = len(x) // 2
      trend = np.ptp(x[h:]) / max(np.ptp(x[:h]), 1e-12)
      fx = peak(x, dts)
      fV = peak(V, dts) if np.ptp(V) > 1e-3 else float("nan")
      cl = classify(L[j], np.ptp(V), np.ptp(x), trend, fx, fV)
      res1[lab] = (L[j], cl)
      out("%-15s | %+.4f  %+.4f  %+.4f  %+.4f | %7.4f %6.4f | %8.4f  %8.4f  %7.4f | %s"
          % (lab, L[j][0], L[j][1], L[j][2], L[j][3], np.ptp(x), np.ptp(V), fx, fV, fx / fV if fV == fV else float("nan"), cl))
  out("Hinweis: Ein Lyapunov-Exponent gilt hier als null, wenn |λ| < 0,01/µs (endliche Mittelungszeit).")

  sub("1c Mechanische Mode auf dem elektrothermischen Grenzzyklus: Floquet-Rate (Neimark-Sacker)")
  VNS = np.array([3.52, 3.50, 3.48])
  N = len(VNS)
  P = tuple(np.array([par1(v)[k] for v in VNS], dtype=float) for k in range(8))
  Z = np.array([fp1(par1(v)) for v in VNS]).T + np.array([0.01, 0.05, 0.05, 0.0])[:, None]
  dt = 0.01
  xs = []
  for i in range(int(6000 / dt)):
      Z = rk4_step(f1, Z, P, dt)
      if i % 10 == 0:
          xs.append(Z[2].copy())
  xs = np.array(xs).reshape(6, -1, N)
  rates = []
  for k in range(N):
      amp = np.array([np.ptp(xs[s, :, k]) for s in range(6)])
      tt = 500 + 1000 * np.arange(6)
      rr = np.polyfit(tt[1:], np.log(amp[1:]), 1)[0]
      rates.append(rr)
      out("Vs = %.2f V: Amplitude der Brücke je 1000 µs: %s  -> Floquet-Rate %+.5f /µs"
          % (VNS[k], " ".join("%.4f" % a for a in amp), rr))
  cNS = np.polyfit(VNS, rates, 1)
  VNS0 = -cNS[1] / cNS[0]
  out("lineare Extrapolation der Rate auf null: V_NS ≈ %.3f V (Neimark-Sacker, Geburt des T²)" % VNS0)
  out("mechanische Dämpfung am Fixpunkt: -ω_m/(2Q) = %.5f /µs" % (-np.pi * FM / QM))

# =============================================================================
head("TEIL 2  Kippstelle = Divergenz der Rückkopplungssumme")
# =============================================================================
# stationär: G (T - T0) = P(T),  P = Vs^2 R/(R+RL)^2  (VO2 mit Lastwiderstand, ohne Kapazität)


def Pst(T, Vs, RL):
    r = Rvo2(T)
    return Vs ** 2 * r / (r + RL) ** 2


def dlnR(T):
    e = np.exp((T - TC) / WTR)
    return -LNRR * e / (1 + e) ** 2 / WTR


def qst(T, Vs, RL):
    r = Rvo2(T)
    return Vs ** 2 * (RL - r) / (r + RL) ** 3 * r * dlnR(T) / GTH


TG = np.linspace(300.5, 375, 300001)


def Vbranch(RL):
    r = Rvo2(TG)
    return np.sqrt(GTH * (TG - T0)) * (r + RL) / np.sqrt(r)


def folds(RL):
    V = Vbranch(RL)
    d = np.diff(V)
    idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    return [(TG[i + 1], V[i + 1]) for i in idx]


def Tstat(Vs, RL, which):
    V = Vbranch(RL)
    idx = np.where(np.sign(V[:-1] - Vs) != np.sign(V[1:] - Vs))[0]
    roots = [TG[i] + (Vs - V[i]) * (TG[i + 1] - TG[i]) / (V[i + 1] - V[i]) for i in idx]
    return roots[0] if which == "low" else roots[-1]


sub("2a Identität: statische Antwort = Quelle · Σ q^n (isolierender Ast, R_L = 0,8 R_c)")
lo, hi = 1.0, 20.0
for _ in range(60):
    m = 0.5 * (lo + hi)
    if len(folds(m)) >= 2:
        lo = m
    else:
        hi = m
RCUSP = lo
RLb = 0.8 * RCUSP
VFOLD = max(v for _, v in folds(RLb))
for dd in [30e-3, 3e-3, 0.3e-3]:
    Vs = VFOLD - dd
    T = Tstat(Vs, RLb, "low")
    q = qst(T, Vs, RLb)
    dV = 1e-7
    dT = (Tstat(Vs + dV, RLb, "low") - Tstat(Vs - dV, RLb, "low")) / (2 * dV)
    src = 2 * Pst(T, Vs, RLb) / Vs / GTH
    parts = [src * sum(q ** n for n in range(Nn + 1)) for Nn in (0, 1, 5, 50, 500)]
    out("Vs = Vs_f - %.1f mV: T* = %.3f K, q = %.4f, dT*/dVs numerisch = %.3f K/V, (∂P/∂Vs)/G/(1-q) = %.3f K/V"
        % (1e3 * dd, T, q, dT, src / (1 - q)))
    out("      Teilsummen (∂P/∂Vs)/G·Σ_{n≤N} q^n, N = 0, 1, 5, 50, 500: %s" % ", ".join("%.3f" % s_ for s_ in parts))

sub("2b Spitze (Cusp) und Hysterese-Breite")
Pp = GTH * (TG - T0)
Vv = np.sqrt(Pp * Rvo2(TG))
Ii = np.sqrt(Pp / Rvo2(TG))
Rd = np.gradient(Vv, Ii)
out("Spitze bei R_c = %.5f kΩ;  größter negativer differentieller Widerstand des Elements max(-dV/dI) = %.5f kΩ"
    % (RCUSP, (-Rd).max()))
out("=> bistabil genau dann, wenn R_L < max|R_d| (Lastgerade flacher als der steilste NDR-Abschnitt)")
fo = folds(RLb)
out("R_L = 0,8 R_c = %.3f kΩ: Falten bei" % RLb)
for (T, V) in fo:
    out("    T = %.3f K, Vs = %.5f V, Schleifenverstärkung q = %.5f" % (T, V, qst(T, V, RLb)))
eps = np.logspace(-3, -1.5, 7)
W, J, Tf = [], [], []
for e in eps:
    (T1, V1), (T2, V2) = folds(RCUSP * (1 - e))
    W.append(abs(V1 - V2))
    Vf, Tfo = (V1, T1) if V1 > V2 else (V2, T2)
    V = Vbranch(RCUSP * (1 - e))
    cand = TG[np.where(np.diff(np.sign(V - Vf)))[0]]
    J.append(np.max(np.abs(cand - Tfo)))
W = np.array(W)
J = np.array(J)
sW = np.polyfit(np.log(eps), np.log(W), 1)[0]
sJ = np.polyfit(np.log(eps), np.log(J), 1)[0]
out("ε = 1 - R_L/R_c | Breite ΔVs [mV] | Temperatursprung [K]")
for e, w_, j_ in zip(eps, W, J):
    out("   %.2e     |   %10.5f   |  %.4f" % (e, 1e3 * w_, j_))
out("Exponenten (Anpassung log-log): Breite %.3f (Normalform 3/2), Sprung %.3f (Normalform 1/2)" % (sW, sJ))

sub("2c Kritische Verlangsamung vor der Falte")
(T1, V1), (T2, V2) = fo
Vup = max(V1, V2)   # obere Falte: isolierender Ast endet
out("isolierender Ast endet bei Vs_f = %.5f V" % Vup)
out("Vs_f - Vs [mV] | q       | λ = G(1-q)/C_th [1/µs] | λ aus Simulation [1/µs]")
dl, lam_a, lam_s = [], [], []
for d in [100e-3, 30e-3, 10e-3, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5]:
    Vs = Vup - d
    T = Tstat(Vs, RLb, "low")
    q = qst(T, Vs, RLb)
    la = GTH * (1 - q) / CTH
    # Simulation: kleine Störung, Abklingrate
    Tt = T + 1e-4
    dts2 = 0.01 / max(la, 0.01)
    n = int(3.0 / la / dts2)
    for _ in range(n):
        Tt = Tt + dts2 * (Pst(Tt, Vs, RLb) - GTH * (Tt - T0)) / CTH
    ls = -math.log((Tt - T) / 1e-4) / (n * dts2)
    dl.append(d)
    lam_a.append(la)
    lam_s.append(ls)
    out("  %8.3f     | %.5f | %.5f                | %.5f" % (1e3 * d, q, la, ls))
sl = np.polyfit(np.log(dl[-4:]), np.log(lam_s[-4:]), 1)[0]
out("Exponent (Abstände ≤ 1 mV) λ ∝ (Vs_f - Vs)^%.3f (Normalform 1/2);  Relaxationszeit bei 0,3 mV Abstand: %.1f µs, bei 0,03 mV: %.1f µs, gegen τ_th = %.1f µs"
    % (sl, 1 / lam_s[5], 1 / lam_s[-1], CTH / GTH))

sub("2d Normalform x^3 - a x - b = 0")
for a in [0.1, 1.0, 4.0]:
    xf = math.sqrt(a / 3)
    bf = 2 * (a / 3) ** 1.5
    out("a = %.1f: Falten bei b = ±%.5f, Breite 4(a/3)^{3/2} = %.5f, Sprung 3·sqrt(a/3) = %.5f" % (a, bf, 2 * bf, 3 * xf))
a = 1.0
xf = math.sqrt(a / 3)
bf = -2 * (a / 3) ** 1.5    # obere Falte des unteren Astes? prüfe beide numerisch
out("Relaxationsrate am stabilen Ast vor der Falte (a = 1): 1-q = 3x²-a gegen 2(3a)^{1/4}(b_f-b)^{1/2}")
for db in [1e-2, 1e-3, 1e-4]:
    b = 2 * (a / 3) ** 1.5 - db
    rts = np.roots([1, 0, -a, -b])
    rts = np.sort(rts[np.abs(rts.imag) < 1e-9].real)
    xs_ = rts[0]  # unterer Ast (endet bei b_f = +2(a/3)^{3/2})
    out("   b_f - b = %.0e: 1-q = %.6f, Formel %.6f" % (db, 3 * xs_ ** 2 - a, 2 * (3 * a) ** 0.25 * math.sqrt(db)))

sub("2e Galaxien: q = exp(-sqrt(y)), y = g_N/a0")
A0 = 1.2e-10
GN = 6.674e-11
MSUN = 1.989e30
KPC = 3.0857e19
MYR = 3.156e13
yy = np.logspace(-4, 2, 2001)
qq = np.exp(-np.sqrt(yy))
nu = 1 / (1 - qq)
dlog = 1 - qq * np.sqrt(yy) / (2 * (1 - qq))
out("q < 1 für alle y > 0: max q im Gitter = %.6f (bei y = 1e-4);  d ln g/d ln g_N liegt in [%.4f, %.4f]"
    % (qq.max(), dlog.min(), dlog.max()))
out("=> g(g_N) ist streng monoton, keine Falte, keine Bistabilität, kein Hysterese-Ast in der RAR")
out("kleines y: 1 - q = %.5f bei y = 1e-4 gegen sqrt(y) = %.5f" % (1 - math.exp(-0.01), 0.01))
Mb = 6e10 * MSUN
ah = 3.0 * KPC
out("Milchstraße: M_b = 6·10^10 M☉, Hernquist a = 3 kpc (Annahmen)")
out("r [kpc] | y       | q      | ν=1/(1-q) | dlng/dlng_N | dg/dg_N | v_c [km/s] | t_dyn=r/v_c [Myr] | τ_rel=ν·t_dyn [Myr] | ν·r/c [Myr]")
MW = []
for r_k in [10, 20, 30, 40, 50]:
    r = r_k * KPC
    Mr = Mb * r ** 2 / (r + ah) ** 2
    gN = GN * Mr / r ** 2
    y = gN / A0
    q = math.exp(-math.sqrt(y))
    n_ = 1 / (1 - q)
    dl_ = 1 - q * math.sqrt(y) / (2 * (1 - q))
    vc = math.sqrt(n_ * gN * r)
    tdyn = r / vc / MYR
    MW.append((r_k, y, q, n_, dl_, n_ * dl_, vc / 1e3, tdyn, n_ * tdyn, n_ * r / 2.998e8 / MYR))
    out("  %3d   | %.4f  | %.4f | %6.3f    | %.4f      | %6.3f  | %7.1f    | %7.1f           | %8.1f            | %.4f"
        % MW[-1])
out("Antwort der äußeren Beschleunigung auf eine Massenänderung δM/M = 1 % (nach der Relaxation): δg/g = dlng/dlng_N · 1 %:")
out("   " + ", ".join("r = %d kpc: %.3f %%" % (m[0], m[4]) for m in MW))
out("Nachlauf nach einer sprunghaften Änderung: δg(t)/δg(∞) = 1 - exp(-t/τ_rel);  nach t = t_dyn noch Rest exp(-1/ν) = "
    + ", ".join("%.3f" % math.exp(-1 / m[3]) for m in MW))

# =============================================================================
head("TEIL 3  Hysterese als Gedächtnis ohne Speicher: Preisach gegen Gedächtniskern")
# =============================================================================
rng = np.random.default_rng(31)
nG = 120
ag = np.linspace(-1, 1, nG)
AA, BB = np.meshgrid(ag, ag, indexing="ij")
mask = AA >= BB
alpha = AA[mask]
beta = BB[mask]
# Gewichtsverteilung: Mittelpunkt c=(α+β)/2 ~ N(0, 0.3), halbe Breite h=(α-β)/2 ~ N(0.25, 0.15)
cc = 0.5 * (alpha + beta)
hh = 0.5 * (alpha - beta)
mu = np.exp(-cc ** 2 / (2 * 0.3 ** 2) - (hh - 0.25) ** 2 / (2 * 0.15 ** 2))
mu /= mu.sum()
out("Hysteronen: %d (Gitter %d×%d, α ≥ β), μ ∝ exp(-c²/2·0,3² - (h-0,25)²/2·0,15²)" % (len(mu), nG, nG))


def preisach_run(u_seq, state=None):
    s = -np.ones(len(mu)) if state is None else state.copy()
    outp = []
    for u in u_seq:
        s[u >= alpha] = 1.0
        s[u <= beta] = -1.0
        outp.append(float(mu @ s))
    return np.array(outp), s


def reversal_points(u):
    d = np.sign(np.diff(u))
    d[d == 0] = 1
    idx = np.where(d[1:] != d[:-1])[0] + 1
    return np.concatenate([u[idx], u[-1:]])


def alternate(seq):
    res = [seq[0]]
    for r in seq[1:]:
        if r == res[-1]:
            continue
        if len(res) >= 2 and (r - res[-1]) * (res[-1] - res[-2]) > 0:
            res[-1] = r
        else:
            res.append(r)
    return res


def extremum_stack(revs, track=None):
    st = []
    for r in alternate(list(revs)):
        if len(st) >= 2 and (r - st[-1]) * (st[-1] - st[-2]) > 0:
            st[-1] = r
        else:
            st.append(r)
        while len(st) >= 3 and (st[-1] - st[-2]) * (st[-3] - st[-2]) > 0 and abs(st[-1] - st[-2]) >= abs(st[-3] - st[-2]):
            last = st.pop()
            st.pop()
            st.pop()
            if len(st) >= 2 and (last - st[-1]) * (st[-1] - st[-2]) > 0:
                st[-1] = last
            else:
                st.append(last)
        if track is not None:
            track.append(len(st))
    return st


def ramp(points, n=40):
    seq = []
    for a_, b_ in zip(points[:-1], points[1:]):
        seq.extend(np.linspace(a_, b_, n, endpoint=False))
    seq.append(points[-1])
    return np.array(seq)


sub("3a Rückkehrpunkt-Gedächtnis und Auslöschung")
hA = [-1.0, 0.8, -0.6, 0.5, -0.3, 0.2]
hB = [-1.0, 0.8, -0.6, 0.4, -0.2, 0.35, -0.1, 0.5, -0.3, 0.2]
uA = ramp(hA)
uB = ramp(hB)
yA, sA = preisach_run(uA)
yB, sB = preisach_run(uB)
out("Geschichte A (Extremwerte): %s" % hA)
out("Geschichte B (Extremwerte): %s" % hB)
out("Extremwertstapel A: %s" % [round(x, 3) for x in extremum_stack(hA)])
out("Extremwertstapel B: %s" % [round(x, 3) for x in extremum_stack(hB)])
out("Zustand aller Hysteronen gleich: %s;  Ausgang A = %.10f, B = %.10f" % (np.array_equal(sA, sB), yA[-1], yB[-1]))
fut = ramp([0.2] + list(rng.uniform(-0.9, 0.9, 30)), n=20)
fA, _ = preisach_run(fut, sA)
fB, _ = preisach_run(fut, sB)
out("gleiche Zukunft (30 zufällige Umkehrpunkte): max |y_A - y_B| = %.2e" % np.max(np.abs(fA - fB)))
# Rückkehrpunkt
s0 = preisach_run(ramp([-1.0, 0.5]))[1]
y0 = float(mu @ s0)
yl, sl_ = preisach_run(ramp([0.5, -0.2, 0.5]), s0)
out("Rückkehrpunkt: y(0,5) vor Nebenschleife = %.10f, nach 0,5 → -0,2 → 0,5 = %.10f, Zustand gleich: %s"
    % (y0, yl[-1], np.array_equal(s0, sl_)))
# Rekonstruktion aus dem Stapel
T_hist = 200000
u = np.zeros(T_hist)
for i in range(1, T_hist):
    u[i] = np.clip(u[i - 1] + rng.normal(0, 0.03), -0.95, 0.95)
revs = reversal_points(u)
st_len = []
stk = extremum_stack(np.concatenate([[-1.0], revs]), st_len)
_, s_full = preisach_run(np.concatenate([[-1.0], u]))
_, s_red = preisach_run(np.array(stk))
out("Zufallsweg: %d Schritte, %d Umkehrpunkte; Stapel am Ende %d Werte (Maximum unterwegs %d)"
    % (T_hist, len(revs), len(stk), max(st_len)))
out("Zustand aus nur %d Stapelwerten = Zustand aus %d Eingangswerten: %s (%d Hysteronen verglichen)"
    % (len(stk), T_hist, np.array_equal(s_full, s_red), len(mu)))

sub("3b Ratenunabhängigkeit: Preisach gegen linearen Kern K(τ) = e^{-τ/τ0}/τ0")


def kernel_out(useq, dtk, tau0):
    y = 0.0
    ys = []
    a = math.exp(-dtk / tau0)
    for uu in useq:
        y = a * y + (1 - a) * uu
        ys.append(y)
    return np.array(ys)


for n_ in (40, 400):
    uu = ramp(hA, n=n_)
    yp, _ = preisach_run(uu)
    yk = kernel_out(uu, 1.0, 50.0)
    out("Geschichte A mit %3d Schritten je Ast: Preisach-Endwert %.6f, Kern-Endwert (τ0 = 50) %.6f" % (n_, yp[-1], yk[-1]))
ykA = kernel_out(ramp(hA), 1.0, 50.0)
ykB = kernel_out(ramp(hB), 1.0, 50.0)
out("Kern-Endwert A = %.6f, B = %.6f (Preisach: gleich, s. o.)" % (ykA[-1], ykB[-1]))

sub("3c Thermisch aktivierte Hysteronen: wann wird Preisach zum Kern?")
# Zweiniveau-Hysteronen: Barriere E_up = E0 (1-(u-c)/h)^2, E_down = E0 (1+(u-c)/h)^2 für |u-c|<h
# Raten k = ν0 exp(-E/kT), Mastergleichung für p = P(oben). E0 = e0·h.
nc, nh = 24, 24
cgrid = np.linspace(-0.8, 0.8, nc)
hgrid = np.linspace(0.02, 0.6, nh)
Cc, Hh = np.meshgrid(cgrid, hgrid, indexing="ij")
Cc = Cc.ravel()
Hh = Hh.ravel()
wgt = np.exp(-Cc ** 2 / (2 * 0.3 ** 2) - (Hh - 0.25) ** 2 / (2 * 0.15 ** 2))
wgt /= wgt.sum()
NU0 = 1e6
E0 = 40.0   # Barriere in Einheiten der Energieskala pro Einheit h


def rates(uu, kT):
    z = (uu - Cc) / Hh
    Eu = np.where(z < 1, E0 * Hh * np.clip(1 - z, 0, None) ** 2, 0.0)
    Ed = np.where(z > -1, E0 * Hh * np.clip(1 + z, 0, None) ** 2, 0.0)
    return NU0 * np.exp(-Eu / kT), NU0 * np.exp(-Ed / kT)


def ens_run(useq, dtu, kT, p):
    ys = []
    for uu in useq:
        ku, kd = rates(uu, kT)
        kk = ku + kd
        pe = ku / kk
        p = pe + (p - pe) * np.exp(-kk * dtu)
        ys.append(float(wgt @ (2 * p - 1)))
    return np.array(ys), p


out("Ensemble: %d Hysteronen, ν0 = 1e6 pro Zeiteinheit, Barriere E_b = %.0f·h·(1∓(u-c)/h)², Rampen mit 200 Schritten der Dauer 1" % (len(Cc), E0))
out("Grenze aktiv/eingefroren: E_b(u=c) = kT·ln(ν0·t_obs), t_obs = 1e3 -> kT·%.1f" % math.log(NU0 * 1e3))
out("kT   | eingefroren | Remanenz y(0|von +1) - y(0|von -1) | Rückkehrpunkt-Fehler | Ratenabhängigkeit (Rampe 10× langsamer)")
for kT in [0.2, 1.0, 3.0]:
    frozen = np.mean(E0 * Hh > kT * math.log(NU0 * 1e3))
    yp_, _ = ens_run(ramp([-1.0, 1.0, 0.0], n=200), 1.0, kT, np.zeros(len(Cc)))
    ym_, _ = ens_run(ramp([1.0, -1.0, 0.0], n=200), 1.0, kT, np.ones(len(Cc)))
    _, pA = ens_run(ramp([-1.0, 0.5], n=200), 1.0, kT, np.zeros(len(Cc)))
    y0_ = float(wgt @ (2 * pA - 1))
    yl_, _ = ens_run(ramp([0.5, -0.2, 0.5], n=200), 1.0, kT, pA.copy())
    ys_, _ = ens_run(ramp([-1.0, 0.5], n=2000), 1.0, kT, np.zeros(len(Cc)))
    out("%.1f  | %.3f       | %+.4f                             | %.2e             | %.2e"
        % (kT, frozen, yp_[-1] - ym_[-1], abs(yl_[-1] - y0_), abs(ys_[-1] - y0_)))


def lin_test(kT, p_start, amps):
    ku, kd = rates(0.0, kT)
    Kk = ku + kd
    du = 1e-6
    ku2, kd2 = rates(du, kT)
    ku1, kd1 = rates(-du, kT)
    dku = (ku2 - ku1) / (2 * du)
    dkd = (kd2 - kd1) / (2 * du)
    # Zustand bei u = 0 nach der Vorgeschichte relaxieren lassen (Dauer 1e3)
    _, p0 = ens_run(np.zeros(100), 10.0, kT, p_start.copy())
    sgain = dku * (1 - p0) - dkd * p0
    act = Kk > 0
    tt = np.arange(0, 400, 0.1)
    dtu = 0.1
    useq0 = np.sin(2 * np.pi * tt / 37.0) + 0.6 * np.sin(2 * np.pi * tt / 11.3 + 1.0)
    res = []
    for amp in amps:
        useq = amp * useq0
        yN, _ = ens_run(useq, dtu, kT, p0.copy())
        yF, _ = ens_run(np.zeros(len(useq)), dtu, kT, p0.copy())
        yN = yN - yF
        dp = np.zeros(len(Kk))
        yL = []
        ek = np.exp(-Kk * dtu)
        for uu in useq:
            dp = dp * ek + sgain / Kk * (1 - ek) * uu
            yL.append(float(wgt @ (2 * dp)))
        yL = np.array(yL)
        res.append(np.sqrt(np.mean((yN - yL) ** 2)) / np.sqrt(np.mean(yN ** 2)))
    return act.mean(), Kk[act].min() if act.any() else float("nan"), Kk.max(), res


for kT, lab in [(3.0, "Gleichgewicht"), (1.0, "nach Vorgeschichte -1 -> +0,5 -> 0")]:
    if kT == 3.0:
        ps = rates(0.0, kT)[0] / (rates(0.0, kT)[0] + rates(0.0, kT)[1])
    else:
        _, ps = ens_run(ramp([-1.0, 0.5, 0.0], n=200), 1.0, kT, np.zeros(len(Cc)))
    fa, kmin, kmax, errs = lin_test(kT, ps, [0.003, 0.03, 0.3])
    out("Prony-Kern K(τ) = Σ_i 2 w_i s_i e^{-K_i τ}, kT = %.1f (%s): Raten K_i %.2e … %.2e"
        % (kT, lab, kmin, kmax))
    out("   relative Abweichung nichtlinear gegen Kern bei Amplitude 0,003 / 0,03 / 0,3: %s" % " / ".join("%.4f" % e for e in errs))
    out("   Glättungsbreite der Schaltkante kT/|∂E_b/∂u| = kT/(2 E0) = %.4f" % (kT / (2 * E0)))

sub("3d Schließungsregel w(θ) als Preisach-Operator: Kugelkollaps")
out("Kugelkollaps R = A(1-cos η), t = B(η - sin η); θ/(3H_bg) mit H_bg = 2/(3t) (Einstein-de Sitter)")
out("η/π   | R/R_ta | θ/(3H_bg)")
th_list = []
for eta in [0.5, 0.9, 1.0, 1.1, 1.25, 1.4, 1.5]:
    e_ = eta * math.pi
    RR = (1 - math.cos(e_)) / 2
    hr = math.sin(e_) / (1 - math.cos(e_)) ** 2          # Ṙ/R in Einheiten 1/B
    Hb = 2.0 / (3 * (e_ - math.sin(e_)))
    th_list.append((eta, RR, hr / Hb))
    out("%.2f  | %.4f | %+.4f" % (eta, RR, hr / Hb))
THV = th_list[-1][2]
out("Virialisierung (η = 3π/2, R = R_ta/2): θ_min/(3H_bg) = %.4f, unabhängig von der Masse" % THV)
QW = 5.364
out("Preisach-Gewicht: Hysteronen mit Schwellen β (Entkopplung) gleichverteilt in [θ_vir, 0], α = 1 (Wiedereinkopplung nur im Hubble-Fluss)")
out("=> w = 1 - θ_min/θ_vir für θ_vir ≤ θ_min ≤ 0;  Beitrag des projizierten Flusses (lineare Interpolation, Annahme): δg/g = Q·w/ν, Q = %.3f" % QW)
for eta, RR, thv in th_list[2:]:
    wv = 1 - thv / THV
    out("   Schale bei η = %.2fπ (R/R_ta = %.3f): θ_min = %+.3f -> w = %.3f;  δg/g bei ν = 3: %.3f, bei ν = 10: %.3f"
        % (eta, RR, thv, wv, QW * wv / 3, QW * wv / 10))
sig_int = 0.057
lim = 10 ** sig_int - 1
out("RAR-Streuung (intrinsisch ≲ %.3f dex, aus dem Gedächtnis) -> δg/g < %.3f -> w < %.3f·ν/Q:" % (sig_int, lim, lim))
for n_ in [1.0, 2.0, 3.0, 6.0]:
    out("   ν = %.0f: w < %.4f" % (n_, lim * n_ / QW))
out("Auslöschung: ein späterer Verschmelzungsstoß mit θ > θ_min ändert w nicht; nur ein tieferer Kollaps senkt w weiter.")
hist_list = {"früher, ungestörter Kollaps": [1.0, THV, 0.0, -0.5, 0.2, 0.0],
             "derselbe + späte Verschmelzung": [1.0, THV, 0.0, -4.0, 0.3, 0.0],
             "abgebrochener Kollaps (Gezeiten)": [1.0, -3.0, 0.0, -1.0, 0.0]}
for k, hseq in hist_list.items():
    wv = max(1 - min(hseq) / THV, 0.0)
    out("   %-34s θ-Extremwerte %s -> Stapel-Minimum %.3f -> w = %.3f;  δg/g bei ν = 3: %.3f"
        % (k, [round(x, 2) for x in hseq], min(hseq), wv, QW * wv / 3))
out("Rückkopplungssumme im Einfallsbereich eines Systems mit M_b = 6·10^10 M☉ (Punktmasse) und δg/g = Q·w/ν:")
for r_k in [100, 300, 1000]:
    gN = GN * Mb / (r_k * KPC) ** 2
    y = gN / A0
    n_ = 1 / (1 - math.exp(-math.sqrt(y)))
    out("   r = %4d kpc: y = %.2e, ν = %7.2f;  w = 0,803 -> δg/g = %.4f;  w = 0,480 -> δg/g = %.4f"
        % (r_k, y, n_, QW * 0.803 / n_, QW * 0.480 / n_))

# =============================================================================
head("TEIL 4  Thermoelektrischer Transport im Magnetfeld: Onsager-Casimir")
# =============================================================================
E_CH = 1.602e-19
KB = 1.381e-23
L0 = 2.44e-8
TB = 300.0
nb = 2.7e24          # Elektronen = Löcher [m^-3]   (Bi, 300 K, Größenordnung aus dem Gedächtnis)
mue, muh = 1.2, 0.6  # Beweglichkeiten [m²/Vs]     (Bi, 300 K, Größenordnung aus dem Gedächtnis)
Se, Sh = -130e-6, 60e-6   # Band-Seebeck [V/K]  (Annahme, so gewählt, dass S_ges ≈ -67 µV/K)
kL = 2.0             # Gitterwärmeleitung [W/mK] (Annahme)


def band_sigma(n, mu_, s, B):
    b = s * mu_ * B
    return n * E_CH * mu_ / (1 + b * b) * np.array([[1, b], [-b, 1]])


def transport(B):
    sg_e = band_sigma(nb, mue, -1, B)
    sg_h = band_sigma(nb, muh, +1, B)
    sig = sg_e + sg_h
    alp = Se * sg_e + Sh * sg_h
    kbar = kL * np.eye(2) + (L0 * TB + TB * Se ** 2) * sg_e + (L0 * TB + TB * Sh ** 2) * sg_h
    L = np.zeros((4, 4))
    L[:2, :2] = TB * sig
    L[:2, 2:] = TB ** 2 * alp
    L[2:, :2] = TB ** 2 * alp        # α̃(B) = α^T(-B) = α(B) für bandweise α_i = S_i σ_i
    L[2:, 2:] = TB ** 2 * kbar
    return L, sig, alp, kbar


for B in [0.0, 0.1, 1.0, 3.0]:
    L, sig, alp, kbar = transport(B)
    Lm, *_ = transport(-B)
    Ls = 0.5 * (L + L.T)
    La = 0.5 * (L - L.T)
    cas = np.max(np.abs(L.T - Lm)) / np.max(np.abs(L))
    even = np.max(np.abs(Ls - 0.5 * (Lm + Lm.T))) / np.max(np.abs(L))
    odd = np.max(np.abs(La + 0.5 * (Lm - Lm.T))) / np.max(np.abs(L))
    X = rng.normal(size=(4, 20000)) * np.array([1e-3, 1e-3, 1e-6, 1e-6])[:, None]
    sA_ = np.einsum("in,ij,jn->n", X, La, X)
    sS_ = np.einsum("in,ij,jn->n", X, Ls, X)
    rho = np.linalg.inv(sig)
    S_tens = rho @ alp                 # E = ρ α ∇T bei J = 0
    kap_open = kbar - TB * alp @ rho @ alp
    out("")
    out("B = %.1f T:" % B)
    out("   Onsager-Casimir |L^T(B) - L(-B)|/|L| = %.1e;  L_s gerade: %.1e, L_a ungerade: %.1e" % (cas, even, odd))
    out("   Eigenwerte L_s > 0: %s" % (np.all(np.linalg.eigvalsh(Ls) > 0)))
    out("   Entropieproduktion (20000 Zufallskräfte): max|X·L_a X|/max(X·L_s X) = %.1e, min X·L_s X > 0: %s"
        % (np.max(np.abs(sA_)) / np.max(sS_), np.all(sS_ > 0)))
    out("   ρ_xx = %.4e Ωm, Hall ρ_yx = %.3e Ωm, tan θ_H = σ_xy/σ_xx = %+.4f" % (rho[0, 0], rho[1, 0], sig[0, 1] / sig[0, 0]))
    out("   Seebeck S_xx = %+.2f µV/K,  Nernst E_y/∂_xT = %+.3f µV/K  (Koeffizient %s)"
        % (S_tens[0, 0] * 1e6, S_tens[1, 0] * 1e6, ("%.2f µV/(K·T)" % (S_tens[1, 0] * 1e6 / B)) if B > 0 else "-"))
    out("   κ_xx(J=0) = %.3f W/mK, Righi-Leduc κ_yx/κ_xx = %+.4f" % (kap_open[0, 0], kap_open[1, 0] / kap_open[0, 0]))
    if B > 0:
        L0_, sig0, alp0, kbar0 = transport(0.0)
        rho0 = np.linalg.inv(sig0)
        out("   Magnetowiderstand Δρ_xx/ρ_xx = %.3f" % (rho[0, 0] / rho0[0, 0] - 1))
    # Anteil der antisymmetrischen Flüsse
    Xt = np.array([1e-3, 0, 1e-6, 0])
    Ja = La @ Xt
    Js = Ls @ Xt
    out("   Kraft entlang x: |J_a|/|J_s| = %.4f, X·J_a = %.1e, X·J_s = %.3e" % (np.linalg.norm(Ja) / np.linalg.norm(Js), Xt @ Ja, Xt @ Js))
# Bridgman-Beziehung: Ettingshausen P = T N / κ (lineare Ordnung in B)
B = 0.01
L, sig, alp, kbar = transport(B)
rho = np.linalg.inv(sig)
Nn = (rho @ alp)[1, 0] / B
# Ettingshausen: J = (Jx,0), ∇T = (0, gy), J_Q,y = 0 -> gy; E aus J = σE - α∇T
# unbekannte: E (2), gy; Gleichungen: J_y = 0, J_x = 1, J_Q,y = 0
Jx = 1.0


def resid(v):
    Ex, Ey, gy = v
    E = np.array([Ex, Ey])
    gT = np.array([0.0, gy])
    J = sig @ E - alp @ gT
    JQ = TB * alp @ E - kbar @ gT
    return np.array([J[0] - Jx, J[1], JQ[1]])


Mx = np.zeros((3, 3))
r0 = resid(np.zeros(3))
for i in range(3):
    e = np.zeros(3)
    e[i] = 1.0
    Mx[:, i] = resid(e) - r0
sol = np.linalg.solve(Mx, -r0)
Pett = sol[2] / (Jx * B)
kap_open = kbar - TB * alp @ rho @ alp
out("")
out("Bridgman (B = 0,01 T): Ettingshausen P = %.4e K·m/(A·T), T·N/κ = %.4e, Verhältnis %.5f"
    % (Pett, TB * Nn / kap_open[0, 0], Pett / (TB * Nn / kap_open[0, 0])))

# =============================================================================
head("TEIL 5  Anwendungen")
# =============================================================================
sub("5a Curie-Punkt als stabile Randlage (P = P0 (1 - T/T_C)^{2β}, 2β = 0,72)")
TCu = 493.0
T0f = 293.0
bexp = 0.72
Cf = 0.75     # Wärmekapazität [J/K]: 1 g Ferrit, c = 750 J/(kg K) (Annahme)
Gf = 0.01     # Wärmeleitwert [W/K] (Annahme)
out("T_C = %.0f K, T0 = %.0f K, C = %.2f J/K, G = %.3f W/K (τ0 = C/G = %.0f s)" % (TCu, T0f, Cf, Gf, Cf / Gf))
out("P0/G [K] | T* [K]   | T_C-T* [K] | q        | 1/(1-q) | dT*/dlnP0 [K] | τ = τ0/(1-q) [s] | σ_T thermisch [K] | δT bei 1 % Leistungsrauschen [K]")
for PG in [1e2, 1e3, 1e4, 1e5]:
    a_, b_ = T0f, TCu - 1e-12
    for _ in range(200):
        m = 0.5 * (a_ + b_)
        if PG * (1 - m / TCu) ** bexp - (m - T0f) > 0:
            a_ = m
        else:
            b_ = m
    Ts_ = 0.5 * (a_ + b_)
    q = -PG * bexp / TCu * (1 - Ts_ / TCu) ** (bexp - 1)
    dTd = PG * (1 - Ts_ / TCu) ** bexp / (1 - q)
    sT = math.sqrt(KB * Ts_ ** 2 / (Cf * (1 - q)))
    out("%8.0f | %8.3f | %10.4f | %8.2f | %7.4f | %13.4f | %16.3f | %.2e          | %.4f"
        % (PG, Ts_, TCu - Ts_, q, 1 / (1 - q), dTd, Cf / Gf / (1 - q), sT, 0.01 * dTd))

sub("5b Kramers-Verweilzeiten und Kac-Budget im bistabilen VO2-Element (R_L = 0,8 R_c)")
(T1, V1), (T2, V2) = folds(RLb)
Vlo, Vhi = min(V1, V2), max(V1, V2)
out("bistabiles Fenster: %.5f V < Vs < %.5f V (Breite %.2f mV)" % (Vlo, Vhi, 1e3 * (Vhi - Vlo)))
TT = np.linspace(320, 360, 400001)
dTT = TT[1] - TT[0]


def potential(Vs):
    fT = (Pst(TT, Vs, RLb) - GTH * (TT - T0)) / CTH
    U = -np.cumsum(fT) * dTT
    return fT, U


def wells(Vs):
    fT, U = potential(Vs)
    idx = np.where(np.sign(fT[:-1]) != np.sign(fT[1:]))[0]
    return [TT[i] for i in idx], fT, U


# Spannung mit gleich hohen Barrieren suchen
def barrier_diff(Vs):
    r, fT, U = wells(Vs)
    i = [np.argmin(np.abs(TT - x)) for x in r]
    return (U[i[1]] - U[i[0]]) - (U[i[1]] - U[i[2]])


a_, b_ = Vlo + 1e-6, Vhi - 1e-6
for _ in range(60):
    m = 0.5 * (a_ + b_)
    if barrier_diff(a_) * barrier_diff(m) <= 0:
        b_ = m
    else:
        a_ = m
Veq = 0.5 * (a_ + b_)
for Vs, dlist in [(Veq, [3.0, 4.0]), (Veq + 0.25 * (Vhi - Veq), [3.0])]:
    r, fT, U = wells(Vs)
    i = [np.argmin(np.abs(TT - x)) for x in r]
    dU1 = U[i[1]] - U[i[0]]
    dU2 = U[i[1]] - U[i[2]]
    c1 = -np.gradient(fT, dTT)[i[0]]
    c2 = -np.gradient(fT, dTT)[i[2]]
    cb = np.gradient(fT, dTT)[i[1]]
    out("")
    out("Vs = %.5f V: Mulden T1 = %.3f K, T2 = %.3f K, Sattel %.3f K; Barrieren ΔU1 = %.4f, ΔU2 = %.4f K²/µs"
        % (Vs, r[0], r[2], r[1], dU1, dU2))
    out("   Krümmungen U''(T1) = %.4f, U''(T2) = %.4f, |U''(Sattel)| = %.4f /µs" % (c1, c2, cb))
    # thermisches Eigenrauschen
    # Einheiten: k_B T^2 G / C^2 in K^2/s -> K^2/µs.  G [mW/K]=1e-3 W/K, C [nJ/K]=1e-9 J/K
    Dth = KB * r[0] ** 2 * (GTH * 1e-3) / (CTH * 1e-9) ** 2 * 1e-6
    out("   Eigenrauschen D_th = k_B T² G/C² = %.3e K²/µs -> ΔU1/D_th = %.2e (spontane Sprünge praktisch nie)" % (Dth, dU1 / Dth))
    for Dn in [dU1 / dd_ for dd_ in dlist]:
        kr1 = math.sqrt(c1 * cb) / (2 * math.pi) * math.exp(-dU1 / Dn)
        kr2 = math.sqrt(c2 * cb) / (2 * math.pi) * math.exp(-dU2 / Dn)
        # Simulation Euler-Maruyama
        nw = 200
        dts_ = 0.02
        nsteps = 500000
        Tw = np.where(np.arange(nw) % 2 == 0, r[0], r[2]).astype(float)
        cell = (Tw > r[1]).astype(int)
        # Übergang erst an der jeweils anderen Mulde gezählt
        last = cell.copy()
        tlast = np.zeros(nw)
        dw = [[], []]
        samp_every = int(round(1.0 / dts_))
        lab = []
        sq = math.sqrt(2 * Dn * dts_)
        for k in range(nsteps):
            fTw = (Pst(Tw, Vs, RLb) - GTH * (Tw - T0)) / CTH
            Tw = Tw + fTw * dts_ + sq * rng.normal(size=nw)
            now = np.where(Tw < r[0], 0, np.where(Tw > r[2], 1, last))
            ch = np.where(now != last)[0]
            t = (k + 1) * dts_
            for j in ch:
                if tlast[j] > 0:
                    dw[last[j]].append(t - tlast[j])
                tlast[j] = t
            last = now
            if k % samp_every == 0 and t > 1000.0:
                lab.append((Tw > r[1]).astype(np.int8))
        lab = np.array(lab)
        # Kac: mittlere Rückkehrzeit (in Abtastschritten) in jede Zelle
        kac = []
        pis = []
        for c in (0, 1):
            rt = []
            for j in range(nw):
                idxs = np.where(lab[:, j] == c)[0]
                if len(idxs) > 1:
                    rt.extend(np.diff(idxs))
            kac.append(np.mean(rt))
            pis.append(np.mean(lab == c))
        out("   D = ΔU1/%.1f = %.4f K²/µs: Verweilzeit Mulde 1 sim %.1f µs / Kramers %.1f µs; Mulde 2 sim %.1f µs / Kramers %.1f µs (%d bzw. %d Sprünge)"
            % (dU1 / Dn, Dn, np.mean(dw[0]), 1 / kr1, np.mean(dw[1]), 1 / kr2, len(dw[0]), len(dw[1])))
        out("      Kac: <τ_ret,1> = %.4f, <τ_ret,2> = %.4f Abtastschritte (1 µs); 1/<τ1> + 1/<τ2> = %.5f; Anteile p1 = %.4f, p2 = %.4f"
            % (kac[0], kac[1], 1 / kac[0] + 1 / kac[1], pis[0], pis[1]))
        out("      Verhältnis der Verweilzeiten sim %.3f, Kramers %.3f" % (np.mean(dw[0]) / np.mean(dw[1]), kr2 / kr1))
    S_P = 2 * (dU1 / 4.0) * (CTH * 1e-9) ** 2 * 1e6
    out("   nötiges eingespeistes Leistungsrauschen für D = ΔU1/4: S_P = 2 D C² = %.2e W²/Hz (√S_P = %.2e W/√Hz)"
        % (S_P, math.sqrt(S_P)))

out("")
out("Gesamtrechenzeit %.0f s" % (time.time() - T_START))

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bistabilitaet.txt"), "w", encoding="utf-8") as fh:
    fh.write(_buf.getvalue())
