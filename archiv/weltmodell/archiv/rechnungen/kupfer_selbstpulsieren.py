# -*- coding: utf-8 -*-
"""
Selbstpulsieren des eingespannten Kupferdrahts im Magnetfeld (06, Laborprüfungen).
Spannungsquelle V cos(ω_d t), Kreis R(T)·I = V − Γẋ, Mode m ẍ + c ẋ + k(T) x = Γ I, Wärme C Ṫ = R I² + c ẋ² − G ΔT.

Weg A: Hüllkurve X(t) (x = Re X e^{iωt}), 3 Gleichungen (Re X, Im X, ΔT); stationäre Punkte, Eigenwerte (Hopf = Selbstpulsieren),
       Zeitintegration der Hüllkurve bis auf den Grenzzyklus.
Weg B: die volle Schwingung (x, ẋ, ΔT) ohne Mittelung, direkt integriert (unabhängige Prüfung von A).
Stationäre Punkte zusätzlich über die Impedanz Z_m = Γ²/(c + i(mω − k/ω)) (Weg der quasistatischen Rechnung).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ftm_kern import draht

p = draht()
R0, G, C, m, Gam, alR, k, w0, cm = (p[x] for x in ('R0', 'G_th', 'C_th', 'm', 'Gam', 'alR', 'k', 'w0', 'c_m'))
tau_th = C / G


def out(s=""):
    print(s)


def R(dT):
    return R0 * (1 + alR * dT)


# ---------------------------------------------------------------- Weg A: Hüllkurve
def rhs_env(t, y, wd, V):
    X = y[0] + 1j * y[1]
    dT = y[2]
    Rt = R(dT)
    Ic = (V - 1j * wd * Gam * X) / Rt
    dX = (Gam * Ic - (k(dT) - m * wd**2 + 1j * wd * cm) * X) / (2j * wd * m)
    heat = 0.5 * Rt * abs(Ic)**2 + 0.5 * cm * wd**2 * abs(X)**2
    return [dX.real, dX.imag, (heat - G * dT) / C]


def X_stat(dT, wd, V):
    Rt = R(dT)
    return Gam * V / Rt / (k(dT) - m * wd**2 + 1j * wd * cm + 1j * wd * Gam**2 / Rt)


def f_stat(dT, wd, V):
    X = X_stat(dT, wd, V)
    Rt = R(dT)
    Ic = (V - 1j * wd * Gam * X) / Rt
    return 0.5 * Rt * np.abs(Ic)**2 + 0.5 * cm * wd**2 * np.abs(X)**2 - G * dT


def stat_punkte(wd, V, Tmax=3.0, n=30001):
    xs = np.linspace(0, Tmax, n)
    fv = f_stat(xs, wd, V)
    idx = np.nonzero(np.sign(fv[:-1]) != np.sign(fv[1:]))[0]
    return [brentq(f_stat, xs[i], xs[i + 1], args=(wd, V), xtol=1e-14) for i in idx]


def f_impedanz(dT, wd, V):
    Rt = R(dT)
    Zm = Gam**2 / (cm + 1j * (m * wd - k(dT) / wd))
    I = V / (Rt + Zm)
    return 0.5 * abs(I)**2 * (Rt + Zm.real) - G * dT


def jac(y, wd, V, h=1e-9):
    J = np.zeros((3, 3))
    sc = [abs(y[0]) + abs(y[1]) + 1e-12] * 2 + [1e-3]
    for j in range(3):
        e = np.zeros(3); e[j] = h * sc[j] * 1e3
        J[:, j] = (np.array(rhs_env(0, y + e, wd, V)) - np.array(rhs_env(0, y - e, wd, V))) / (2 * e[j])
    return J


out("Kupferdraht: f_m0 = %.2f Hz, τ_th = C/G = %.4f s, R_m/R0 = %.3f" % (w0 / 2 / math.pi, tau_th, Gam**2 / cm / R0))
Ql = m * w0 / cm / (1 + Gam**2 / cm / R0)
out("belastete Güte Q_l = %.1f, Ringzeit der Amplitude τ_m = 2Q_l/ω = %.4f s, Verhältnis τ_th/τ_m = %.3f" % (Ql, 2 * Ql / w0, tau_th / (2 * Ql / w0)))
out("")
out("1) Stationäre Punkte, Weg Hüllkurve gegen Weg Impedanz (V = 10 mV):")
for dfd in (-40.0, -33.0, -20.0):
    wd = w0 + 2 * math.pi * dfd
    ra = stat_punkte(wd, 10e-3)
    rb = [brentq(f_impedanz, x - 1e-3, x + 1e-3, args=(wd, 10e-3), xtol=1e-14) for x in ra]
    out("   f_d − f_m0 = %+6.1f Hz: Hüllkurve %s | Impedanz %s" % (dfd, ", ".join("%.6f" % x for x in ra), ", ".join("%.6f" % x for x in rb)))

out("")
out("2) Eigenwerte der Hüllkurve längs des heißen Asts (3. Lösung) bis zu seinem Faltpunkt, Raster 0,01·(V/V₁)² Hz;")
out("   Hopf = komplexes Paar mit Re > 0 (der Mittelast ist der Sattel). V₁ = 10 mV·√(851/Q_l): gleiche Schleifenverstärkung wie bei Q_m = 1000, 10 mV")
out("   Q_m    | τ_th/τ_m | V [mV] | Faltpunkt heißer Ast [Hz] | größtes Re λ (Paar) [1/s] | Hopf-Grenze [Hz] | Ringfrequenz dort [Hz]")
from ftm_kern import DRAHT, CU
hopf_tab = {}


def setze(Qm):
    global R0, G, C, m, Gam, alR, k, w0, cm, tau_th
    d = dict(DRAHT); d['Qm'] = Qm
    pp = draht(CU, d)
    R0, G, C, m, Gam, alR, k, w0, cm = (pp[x] for x in ('R0', 'G_th', 'C_th', 'm', 'Gam', 'alR', 'k', 'w0', 'c_m'))
    tau_th = C / G


def ew_heiss(dfd, V):
    wd = w0 + 2 * math.pi * dfd
    st = stat_punkte(wd, V)
    X = X_stat(st[-1], wd, V)
    return st, np.linalg.eigvals(jac(np.array([X.real, X.imag, st[-1]]), wd, V))


for Qm in (1e3, 2e3, 3e3, 1e4, 3e4, 1e5):
    setze(Qm)
    Qlx = m * w0 / cm / (1 + Gam**2 / cm / R0)
    for Vf in (1.0, 2.0):
        V = 10e-3 * Vf * math.sqrt(851.0 / Qlx)
        best, drei, grenze, ring, fold = -1e9, False, None, None, None
        for dfd in np.arange(0.0, -400.0, -0.01 * Vf**2):
            st, ev = ew_heiss(dfd, V)
            if len(st) >= 3:
                drei = True
            if drei and len(st) < 3:
                fold = dfd + 0.01 * Vf**2
                break
            if drei:
                paar = [e for e in ev if abs(e.imag) > 1e-9]
                if paar:
                    rp = max(e.real for e in paar)
                    best = max(best, rp)
                    if rp > 0 and grenze is None:
                        grenze = dfd; ring = max(abs(e.imag) for e in paar) / 2 / math.pi
        hopf_tab[(Qm, V)] = grenze
        out("   %6.0f | %8.3f | %6.3f | %+10.2f                | %+9.3f                 | %-16s | %s" % (
            Qm, tau_th / (2 * Qlx / w0), V * 1e3, fold if fold is not None else float('nan'), best,
            "%+.2f" % grenze if grenze is not None else "keine", "%.3f" % ring if ring is not None else "–"))


def heiss(wd, V):
    st = stat_punkte(wd, V)
    return st, X_stat(st[-1], wd, V)


def re_heiss(dfd, V):
    wd = w0 + 2 * math.pi * dfd
    st, X = heiss(wd, V)
    ev = np.linalg.eigvals(jac(np.array([X.real, X.imag, st[-1]]), wd, V))
    return max(ev.real), max(ev.imag), len(st)


def lauf_B(dfd, V, pert, Tend, nper=48):
    """volle Schwingung, RK4 mit fester Schrittweite (nper Schritte je Periode der Anregung)"""
    wd = w0 + 2 * math.pi * dfd
    st, X0 = heiss(wd, V)
    X0 = X0 * pert
    x, v, T = X0.real, -wd * X0.imag, st[-1] * pert
    h = 2 * math.pi / wd / nper

    def f(t, x, v, T):
        Rt = R(T); I = (V * math.cos(wd * t) - Gam * v) / Rt
        return v, (Gam * I - cm * v - k(T) * x) / m, (Rt * I * I + cm * v * v - G * T) / C
    Ts = []; t = 0.0
    for i in range(int(Tend / h)):
        a1 = f(t, x, v, T); a2 = f(t + h / 2, x + h / 2 * a1[0], v + h / 2 * a1[1], T + h / 2 * a1[2])
        a3 = f(t + h / 2, x + h / 2 * a2[0], v + h / 2 * a2[1], T + h / 2 * a2[2]); a4 = f(t + h, x + h * a3[0], v + h * a3[1], T + h * a3[2])
        x += h / 6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0]); v += h / 6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
        T += h / 6 * (a1[2] + 2 * a2[2] + 2 * a3[2] + a4[2]); t += h
        if i % nper == 0:
            Ts.append((t, T))
    return st, np.array(Ts)

out("")
out("3) Art der Hopf-Grenze (Hüllkurve integriert): 0,0005 Hz dahinter mit Störung 0,002 % über 12/Re λ;")
out("   0,002 Hz und 0,03 Hz davor mit Stößen 0,1 / 0,3 / 1 / 3 % über 8/|Re λ| (höchstens 3000 s). Endzustand heiß/kalt.")
out("   Q_m    | V [mV] | Hopf [Hz]  | Faltpunkt [Hz] | Ringfrequenz [Hz] | 0,0005 Hz dahinter | 0,002 Hz davor          | 0,03 Hz davor")


def endzustand(dd, V, pert, Tend):
    wd = w0 + 2 * math.pi * dd
    st, X = heiss(wd, V)
    y0 = [X.real * pert, X.imag * pert, st[-1] * pert]
    sa = solve_ivp(rhs_env, (0, Tend), y0, args=(wd, V), method='LSODA', rtol=1e-10, atol=1e-16, max_step=0.02)
    Te = sa.y[2, -1]
    return "heiß" if abs(Te - st[-1]) < abs(Te - st[0]) else "kalt"


for (Qm, V), gz in hopf_tab.items():
    if gz is None:
        continue
    setze(Qm)
    fh = brentq(lambda d: re_heiss(d, V)[0], gz, gz + 0.05 + 0.04 * (V / 5e-3)**2, xtol=1e-6)
    d = fh
    while re_heiss(d, V)[2] >= 3:
        d -= 0.002
    ring = re_heiss(fh, V)[1] / 2 / math.pi
    dd = fh - 0.0005
    hinter = endzustand(dd, V, 1.00002, 12 / re_heiss(dd, V)[0])
    davor = []
    for dx in (0.002, 0.03):
        dd = fh + dx
        Tl = min(8 / abs(re_heiss(dd, V)[0]), 3000.0)
        davor.append("/".join(endzustand(dd, V, p_, Tl) for p_ in (1.001, 1.003, 1.01, 1.03)))
    out("   %6.0f | %6.3f | %+9.4f  | %+9.3f      | %8.3f          | %-18s | %-23s | %s" % (Qm, V * 1e3, fh, d + 0.002, ring, hinter, davor[0], davor[1]))

out("")
out("4) Gegenprobe mit der vollen Schwingung (ohne Mittelung, RK4), Q_m = 10⁴, Störung 1 %")
setze(1e4)
V = 2 * 10e-3 * math.sqrt(851.0 / (m * w0 / cm / (1 + Gam**2 / cm / R0)))
for dfd in (-33.30, -33.42):
    ra, ia, _ = re_heiss(dfd, V)
    for nper in (48, 96):
        st, Ts = lauf_B(dfd, V, 1.01, 12.0, nper)
        tt, TT = Ts[:, 0], Ts[:, 1]
        sel = (tt > 0.5) & (tt < 4)
        y = TT[sel] - TT[sel][-1]
        zc = tt[sel][1:][(y[:-1] < 0) & (y[1:] >= 0)]
        fr = 1 / np.mean(np.diff(zc)) if len(zc) > 2 else 0.0
        out("   V = %.3f mV, %+.2f Hz, %d Schritte/Periode: Äste %s K, ΔT nach 12 s = %.4f K, Ringfrequenz %.3f Hz (Hüllkurve %.3f Hz, Re λ = %+.3f/s)"
            % (V * 1e3, dfd, nper, ", ".join("%.4f" % x for x in st), TT[-1], fr, ia / 2 / math.pi, ra))
