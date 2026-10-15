# -*- coding: utf-8 -*-
"""
Übergabe der Mode ξ bei der ersten Faltung des Phasenraumblatts: zweiter Weg (selbstähnlich) und Folgen als Stufe.

A  Selbstähnlicher Einfall (Fillmore–Goldreich) auf eine Überdichte im EdS-Hintergrund, eine Bahn statt N Schalen:
   mittlere Anfangsüberdichte δ̄(<r_i) ∝ r_i^{−ε}  ⇔  δM/M ∝ M^{−e}, e = ε/3 (Keim: ε = 3).
   Jede Schale: r = r_ta(M)·Λ(τ), τ = t/t_ta(M), r_ta ∝ t^s, s = 2/(9e) + 2/3, M_ta ∝ t^p, p = 2/(3e).
   Λ'' = −(π²/8)·τ^p·𝓜(Λτ^{−s})/Λ² + (π²/8)·ε_j/Λ³  (ab der Umkehr Drehimpuls j² = ε_j·GM·r_ta wie gedaechtnis_topologie C),
   𝓜(λ) = M(<λ r_ta(t))/M_ta(t) aus der Bahn selbst (alle Schalen sind dieselbe Bahn, verschoben); Iteration bis fest.
   Erste Faltung am Ort der Schale: ihr λ = Λτ^{−s} erreicht den äußersten Rücklauf λ_c = max λ der Schalen nach dem
   ersten Perizentrum (äußerste Kaustik). Bis dahin ist M(<r) = M fest: I_F = ∫ω dt, ω² = GM/r³.
   Radial (ε_j → 0) ist die Bahn bis zur Faltung die Zykloide: I_F = √2·ln tan(η_F/4), (1 − cos η_F)/2 = λ_c·τ_F^s.
B  Vergleich mit dem N-Schalen-Weg (gedaechtnis_topologie C: 0,9–1,2 rad, t_F/t_ta ≈ 1,8, r_F/R_ta ≈ 0,6–0,7).
C  Stufe statt e^{−I/N*} in der Linse um isolierte Galaxien (Verfahren linsen_isoliert: I je Lagrange-Schale).
D  Nicht übergebene Mode heute: Galaxien und Haufen (Stufe gegen Exponent).
E  Handwert.
Nur ftm_kern (über linsen_isoliert). Aufruf: python uebergabe_faltung.py (schreibt uebergabe_faltung.txt).
"""
import io, os, sys, math, contextlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
K8 = math.pi**2 / 8
OUT = []


def P(s=""):
    OUT.append(s); print(s, flush=True)


def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)


# ============================================================================================== A  selbstähnlich
def loese(e, ej, ngrid=6000, nit=40):
    p = 2 / (3 * e); s = 2 / (9 * e) + 2 / 3
    tmax = max(8.0, 0.15 ** (-1 / p))              # Schalen mit τ' > tmax (Masse ≤ 15 % von M_ta) als Zentrum
    tg = np.geomspace(1, tmax, ngrid)
    lg = np.geomspace(1e-6, 1.0, 3000); llg = np.log(lg)
    Mtab = np.minimum(2.0 * lg, 1.0)
    edges = np.concatenate([[1.0], np.sqrt(tg[1:] * tg[:-1]), [tmax]])
    wts = edges[:-1]**(-p) - edges[1:]**(-p); tail = tmax**(-p)
    d = 1.0
    for it in range(nit):
        lM = np.log(np.maximum(Mtab, 1e-300))

        def f(t, y):
            L, V = y
            lam = L * t**(-s)
            if lam >= 1:
                Mv = math.exp(lM[-1])
            elif lam <= lg[0]:
                Mv = math.exp(lM[0]) * lam / lg[0]
            else:
                Mv = math.exp(np.interp(math.log(lam), llg, lM))
            return [V, -K8 * t**p * Mv / L**2 + K8 * ej / L**3]
        sol = solve_ivp(f, (1, tmax), [1.0, 0.0], t_eval=tg, rtol=1e-9, atol=1e-12, method="DOP853")
        Lam = sol.y[0]; lamp = Lam * tg**(-s)
        o = np.argsort(lamp); cm = tail + np.cumsum(wts[o])
        Mnew = np.interp(lg, lamp[o], cm, left=tail, right=cm[-1])
        m = lg > 0.05
        d = np.max(np.abs(Mnew - Mtab)[m] / Mnew[m])
        Mtab = 0.5 * Mtab + 0.5 * Mnew
        if d < 3e-4 and it > 4:
            break
    V = np.gradient(Lam, tg)
    ip = int(np.argmax((V[:-1] < 0) & (V[1:] >= 0)))
    lam_c = lamp[ip:].max()
    iF = int(np.argmax(lamp <= lam_c))
    om = np.sqrt(K8 / Lam[:iF + 1]**3)
    I_F = np.trapezoid(om, tg[:iF + 1])
    I_p = np.trapezoid(np.sqrt(K8 * tg[:ip + 1]**p * np.exp(np.interp(np.log(np.maximum(lamp[:ip + 1], 1e-6)), llg, np.log(Mtab)))
                               / Lam[:ip + 1]**3), tg[:ip + 1])
    # Probe: Masse innerhalb der Schale weicht erst ab der Faltung von M ab
    Min = np.exp(np.interp(np.log(lamp), llg, np.log(Mtab))) * tg**p
    iM = int(np.argmax(np.abs(Min - 1) > 0.01))
    return dict(e=e, ej=ej, lam_c=lam_c, tF=tg[iF], LF=Lam[iF], I_F=I_F, tp=tg[ip], Lp=Lam[ip], I_p=I_p, d=d, it=it,
                tM=tg[iM], s=s)


def I_radial(lam_c, s):
    """Zykloide bis zur Faltung: (1 − cos η)/2 = λ_c·((η − sin η)/π)^s, I = √2 ln tan(η/4)."""
    g = lambda eta: (1 - math.cos(eta)) / 2 - lam_c * ((eta - math.sin(eta)) / math.pi)**s
    eta = brentq(g, math.pi + 1e-9, 2 * math.pi - 1e-9)
    return eta, math.sqrt(2) * math.log(math.tan(eta / 4)), (eta - math.sin(eta)) / math.pi


def _job(a):
    return loese(*a)


if __name__ == "__main__":
    EPS = (0.6, 1.0, 1.5, 2.0, 3.0)                # ε des Anfangsprofils δ̄ ∝ r^{−ε}
    EJ = (0.01, 0.05, 0.2)
    jobs = [(eu / 3, ej) for eu in EPS for ej in EJ] + [(1.0, 0.1)]
    if os.environ.get("SCHNELL"):
        jobs = [(1.0, 0.1)]; EPS = (3.0,); EJ = (0.1,)
    with Pool(4) as pool:
        RES = pool.map(_job, jobs)
    R = {(round(r["e"] * 3, 3), r["ej"]): r for r in RES}

    hdr("A  Selbstähnlicher Einfall, eine Bahn (EdS): erste Faltung am Ort der Schale, I_F = ∫ω dt ab der Umkehr")
    P("   ε (δ̄ ∝ r^−ε)   e = ε/3   ε_j    λ_c = r_c/r_ta(t)   t_F/t_ta   r_F/R_ta   I_F [rad]   I_F radial (Zykloide, λ_c)   t_Peri/t_ta   I bis Peri   Konv.")
    for eu in EPS:
        for ej in EJ:
            r = R[(round(eu, 3), ej)]
            eta, Irad, tr = I_radial(r["lam_c"], r["s"])
            P("   %4.1f          %.3f     %.2f    %.4f             %.3f      %.3f      %.3f       %.3f (η_F = %.3f)          %.3f         %.3f      %.0e"
              % (eu, eu / 3, ej, r["lam_c"], r["tF"], r["LF"], r["I_F"], Irad, eta, r["tp"], r["I_p"], r["d"]))
    P("   Probe: Masse innerhalb der Schale weicht erst ab t/t_ta = " +
      ", ".join("%.3f" % R[(round(eu, 3), EJ[len(EJ) // 2])]["tM"] for eu in EPS) + " um > 1 %% ab (ε_j = %.2f; vgl. t_F oben)." % EJ[len(EJ) // 2])
    Ivals = [R[k]["I_F"] for k in R]
    P("   ⇒ I_F = %.2f–%.2f rad über ε = 0,6–3 und ε_j = 0,01–0,2; die Faltung liegt vor dem eigenen Perizentrum." % (min(Ivals), max(Ivals)))

    hdr("B  Vergleich mit dem N-Schalen-Weg (gedaechtnis_topologie C, Keim: ε = 3)")
    P("   ε_j    selbstähnlich: t_F/t_ta   r_F/R_ta   I_F      N-Schalen (M/M_Keim = 12 / 18 / 26): I_F")
    NS = {0.05: "Anlauf / 0,383 / 0,904", 0.1: "0,948 / 0,971 / 1,045", 0.2: "1,136 / 1,129 / 1,113"}
    for ej in [x for x in (0.05, 0.1, 0.2) if (3.0, x) in R]:
        r = R[(3.0, ej)]
        P("   %.2f          %.3f        %.3f     %.3f    %s" % (ej, r["tF"], r["LF"], r["I_F"], NS[ej]))
    P("   N-Schalen: t_F/t_ta = 1,71–1,84, r_F/R_ta = 0,61–0,69 (ε_j = 0,1–0,2).")

    # ========================================================================================== C  Linse
    import linsen_isoliert as L
    L.MB = [m * L.Ms for m in L.MB]                # wie in linsen_isoliert.__main__ (Massen in kg)
    import ftm_kern as Kk
    G, Ms, kpc, Mpc, Gyr = L.G, L.Ms, L.kpc, L.Mpc, L.Gyr
    _orig = L.xi_R

    def xi_stufe(r, Mb, prof, z, Ns):
        if isinstance(Ns, tuple):                 # ("Stufe", I_F): Rest 1 bis zur Faltung, danach 0
            I = L.I_profil(Mb, prof, z)
            Ir = np.exp(np.interp(np.log(r), np.log(L.R_XI), np.log(np.maximum(I, 1e-300)), left=np.log(I[0]), right=np.log(I[-1])))
            return (Ir < Ns[1]).astype(float)
        return _orig(r, Mb, prof, z, Ns)
    L.xi_R = xi_stufe
    IF_lo, IF_mid, IF_hi = min(Ivals), float(np.median(Ivals)), max(Ivals)
    STUFEN = (("Stufe", IF_lo), ("Stufe", IF_mid), ("Stufe", IF_hi))

    hdr("C  Linse um isolierte Galaxien (Milchstraße 6·10¹⁰ M☉, z = 0): δg/g = g/g_Regel − 1, Stufe bei I_F gegen e^{−I/N*}")
    Mmw = 6e10 * Ms
    RK = (100, 200, 300, 500, 700, 1000)
    for prof in ("flach", "NFW"):
        I = L.I_profil(Mmw, prof, 0.0)
        Ir = [math.exp(np.interp(math.log(rk * kpc), np.log(L.R_XI), np.log(I))) for rk in RK]
        rF = [math.exp(np.interp(-x, -I, np.log(L.R_XI))) / kpc for x in (IF_lo, IF_mid, IF_hi)]
        P("   %-5s I(r) = ∫(1 − w)ω dt bei %s kpc: %s;  Faltungsradius heute r_F (I = %.2f / %.2f / %.2f): %s kpc"
          % (prof, "/".join(map(str, RK)), " / ".join("%.2f" % x for x in Ir), IF_lo, IF_mid, IF_hi, " / ".join("%.0f" % x for x in rF)))
    for form in ("lokal", "kumulativ"):
        P("   Form %s (3D):  r [kpc] = %s" % (form, " / ".join(map(str, RK))))
        for prof in ("flach", "NFW"):
            for Ns in (0.3, 2.4) + STUFEN:
                r = np.array(RK) * kpc
                x = L.g_3D(Mmw, r, Ns=Ns, prof=prof, z=0.0, form=form) / L.g_3D(Mmw, r) - 1
                nm = ("N* = %.1f" % Ns) if not isinstance(Ns, tuple) else ("Stufe I_F = %.2f" % Ns[1])
                P("     %-5s %-17s " % (prof, nm) + " ".join("%+.3f" % v for v in x))
    P("   χ² gegen KiDS (4 Bins, Fenster R ≤ 0,30 / 1,04 / 2,60 Mpc; linsen_isoliert: Regel 104,9 / 429,6 / 461,8)")
    for nm, kw in [("N* = 0,3 flach", dict(Ns=0.3)), ("N* = 2,4 flach", dict(Ns=2.4))] + \
                  [("Stufe %.2f %s" % (st[1], pr), dict(Ns=st, prof=pr)) for st in STUFEN for pr in ("flach", "NFW")]:
        res = [L.chi2_kids(kw, rm) for rm in (0.3035, 1.0365, 2.6035)]
        P("     %-22s " % nm + "  ".join("%7.1f/%d" % (c_, n) for c_, n in res))

    P("   Formtest (je Bin eigene Massenskala f_M wie linsen_isoliert 2c; dort Regel 44,8 / 81,5 / 138,7, N* = 0,3: 44,7 / 80,0 / 101,4)")
    for nm, kw in [("Stufe %.2f %s" % (st[1], pr), dict(Ns=st, prof=pr)) for st in STUFEN[:2] for pr in ("flach", "NFW")]:
        out = []
        for rm in (0.3035, 1.0365, 2.6035):
            ch, n, fs = L.chi2_frei(kw, rm)
            out.append("%6.1f/%d (" % (ch, n) + " ".join("%.2f" % f for f in fs) + ")")
        P("     %-22s " % nm + "  ".join(out))
    P("   (KiDS/Mistele sind Literaturwerte, hier nicht nachgeprüft.)")

    # ========================================================================================== D  Anteil
    hdr("D  Nicht übergebene Mode heute: Faltungsradius r_F gegen Umkehrradius R_ta (heute) und R₂₀₀")

    def R_ta_heute(Mf, t_obs):
        """Schale, die genau bei t_obs umkehrt, mit M = M_prof(R_ta)."""
        def ta(M, di):
            ri = (3 * M / (4 * math.pi * L.RHO_I * (1 + di))) ** (1 / 3); vi = L.Hb(L.AI) * ri * (1 - di / 3)
            f = lambda t, y: [y[1], -G * M / y[0]**2 + L.OL * L.H0**2 * y[0]]
            ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
            s_ = solve_ivp(f, (L.TI, 3 * t_obs), [ri, vi], events=ev, rtol=1e-9, atol=1e-6)
            return (s_.t_events[0][0], s_.y_events[0][0][0]) if len(s_.t_events[0]) else (np.inf, np.nan)
        def gap(r):
            M = Mf(r)
            ld = brentq(lambda ld: ta(M, 10**ld)[0] - t_obs, -4, 0)
            return ta(M, 10**ld)[1] - r
        return brentq(gap, 50 * kpc, 20 * Mpc, xtol=1e-3 * kpc)

    t0 = L.t_von_a(1.0)
    P("   Galaxien (Profile wie linsen_isoliert, z = 0); Anteil = Masse (Profil) zwischen r_F und R_ta an M(<R_ta):")
    P("   M_b/M☉    Profil   R_ta[kpc]   r_F [kpc] (I_F = %.2f/%.2f/%.2f)   Anteil Stufe            Anteil e^{−I/N*} N* = 0,3 / 2,4" % (IF_lo, IF_mid, IF_hi))
    for Mb_ in (1e10, 6e10, 3e11):
        Mb = Mb_ * Ms
        Mf_, _ = L.massen(Mb)
        for prof in ("flach", "NFW"):
            Rt = R_ta_heute(Mf_[prof], t0)
            I = L.I_profil(Mb, prof, 0.0)
            rF = [math.exp(np.interp(-x, -I, np.log(L.R_XI))) for x in (IF_lo, IF_mid, IF_hi)]
            rr = np.geomspace(1 * kpc, Rt, 400)
            dM = np.gradient(np.array([Mf_[prof](x) for x in rr]), rr)
            Ir = np.exp(np.interp(np.log(rr), np.log(L.R_XI), np.log(I)))
            Mt = np.trapezoid(dM, rr)
            st = [np.trapezoid(dM * (Ir < x), rr) / Mt for x in (IF_lo, IF_mid, IF_hi)]
            ex = [np.trapezoid(dM * np.exp(-Ir / n), rr) / Mt for n in (0.3, 2.4)]
            P("   %.0e   %-5s    %6.0f      %s                  %s          %s"
              % (Mb_, prof, Rt / kpc, " / ".join("%4.0f" % (x / kpc) for x in rF), " / ".join("%.2f" % x for x in st),
                 " / ".join("%.3f" % x for x in ex)))
        # Masse (1 + Q)·M_b fest (Naht in 03: Milchstraße R_ta = 0,58 Mpc)
        Mq = (1 + L.Q) * Mb
        Rt = R_ta_heute(lambda r: Mq, t0)
        rr = np.geomspace(30 * kpc, Rt, 40)
        Ir = np.array([L.I_schale(x, Mq, t0)[0] for x in rr])
        rF = [math.exp(np.interp(-x, -Ir, np.log(rr))) if Ir[-1] < x < Ir[0] else float("nan") for x in (IF_lo, IF_mid, IF_hi)]
        P("   %.0e   (1+Q)M_b %6.0f      %s                  (Masse innen fest: kein Anteil je Masse)   I(R_ta/2) = %.2f"
          % (Mb_, Rt / kpc, " / ".join("%4.0f" % (x / kpc) for x in rF), float(np.interp(Rt / 2, rr, Ir))))

    P("   Haufen (Gesamtmasse NFW M₂₀₀ = 10¹⁵ M☉, c = 4; Schalen wie linsen_isoliert):")
    rho_c0 = L.rho_c
    for zc in (0.0, 0.3):
        rho = rho_c0 * Kk.E_hintergrund(zc)**2
        r200 = (3 * 1e15 * Ms / (4 * math.pi * 200 * rho)) ** (1 / 3)
        mu = lambda x: math.log(1 + x) - x / (1 + x)
        Mf = lambda r: 1e15 * Ms * mu(4 * r / r200) / mu(4)
        tob = L.t_von_a(1 / (1 + zc))
        rows = []
        for fr in (0.5, 1.0, 1.5, 2.0, 3.0):
            I_, _ = L.I_schale(fr * r200, Mf(fr * r200), tob)
            rows.append((fr, I_))
        Rt = R_ta_heute(Mf, tob)
        P("     z = %.1f: R₂₀₀ = %.0f kpc, R_ta = %.0f kpc (%.2f R₂₀₀);  I bei r/R₂₀₀ = %s"
          % (zc, r200 / kpc, Rt / kpc, Rt / r200, ", ".join("%.1f: %.2f" % x for x in rows)))
    P("   Stoßnachlauf (03: x = λ_m t) als Stufe: ganz übergeben ab t_F = I_F/ω_orb")
    for D_ in (200, 500):
        iw = Kk.uebergabezeit(D_, 0.3, (1, 1))[0] / Gyr
        P("     Δ = %d, z = 0,3: 1/ω_orb = %.3f Gyr ⇒ t_F = %.2f / %.2f / %.2f Gyr (I_F = %.2f / %.2f / %.2f)"
          % (D_, iw, IF_lo * iw, IF_mid * iw, IF_hi * iw, IF_lo, IF_mid, IF_hi))

    # ========================================================================================== E  Hand
    hdr("E  Handwert (radial, Keim ε = 3, λ_c = 0,35)")
    P("   η = 4,55: cos η = −0,1627 ⇒ Λ = (1 + 0,1627)/2 = 0,5814;  sin η = −0,9867 ⇒ τ = (4,55 + 0,9867)/π = 1,7624")
    P("   τ^{8/9} = e^{0,889·0,5666} = 1,6548;  λ = 0,5814/1,6548 = 0,3513 ≈ λ_c")
    P("   I_F = √2·ln tan(η/4) = 1,4142·ln tan(1,1375) = 1,4142·ln(0,9075/0,4200) = 1,4142·0,7705 = 1,0896")
    eta, Irad, tr = I_radial(0.35, 2 / 9 + 2 / 3)
    P("   Code (I_radial, λ_c = 0,35): η_F = %.4f, τ_F = %.4f, I_F = %.4f" % (eta, tr, Irad))

    with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
