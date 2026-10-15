# -*- coding: utf-8 -*-
"""
Übergabe der Mode ξ als geglättete Stufe aus der Verteilung der ersten Faltung über die Bahnen einer Schale; Linse um
isolierte Galaxien und Nachlauf nach Haufenstößen ohne freie Übergaberate.

A  I_F(ε, ε_j) aus dem selbstähnlichen Einfall (uebergabe_faltung.loese): erste Faltung am Ort der Bahn, I_F = ∫ω dt ab der Umkehr.
   ε: Steigung des Anfangsprofils; im selbstähnlichen Einfall M(<r) ∝ r^{3/(1+ε)}, also ε = 3/α − 1 mit α = d ln M/d ln r
      des Profils am Ort der Schale (flach: α = 1 ⇒ ε = 2; NFW: ε(r) wächst nach außen).
   ε_j: Drehimpuls der Bahn ab der Umkehr, j² = ε_j·GM·r_ta; er streut über die Bahnen einer Schale.
B  Form: Rest S(I) = Anteil der Bahnen der Schale mit I_F > I (noch nicht gefaltet) = geglättete Stufe.
   Verteilungen von ε_j: log-gleich 0,01–0,2 (Bereich der Faltungsrechnung), isotrop tangential (j² exponentiell, 2 Freiheitsgrade)
   mit Mittel 0,05 bzw. 0,1 (abgeschnitten bei 0,4).  Effektives N*: gleiche mittlere Phase ∫S dI = N* (wie ∫e^{−I/N*} dI);
   daneben Median I_½/ln 2.
C  Linse um isolierte Galaxien (Verfahren linsen_isoliert, kumulativ): ξ_R = S(I(r)) statt e^{−I/N*}; χ² gegen KiDS [G].
D  Haufenstoß: Operatorform an der Gedächtnislage m(t) = ∫K(τ)d(t − τ)dτ, K = −dS/dτ, Gas bei d(t) ∝ t:
   φ = m/d = 1 − (1/x)∫_0^x S(I) dI, x = ω_orb·t  (Exponent: 1 − (1 − e^{−x'})/x', x' = x/N*, wie in 03).
E  Handwerte.
Nur ftm_kern (über linsen_isoliert, uebergabe_faltung). Aufruf: python linse_ohne_freiheit.py (schreibt linse_ohne_freiheit.txt).
"""
import io, os, sys, math
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uebergabe_faltung as U

OUT = []


def P(s=""):
    OUT.append(s); print(s, flush=True)


def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)


def _job(a):
    return U.loese(*a)


EPS = (2.0, 3.0, 5.0)
EJ = (0.01, 0.03, 0.1, 0.2, 0.4)
IG = np.linspace(0, 6, 6001)                       # Gitter der Phase I [rad]


def verteilungen():
    """Gewichte über ein feines ln ε_j-Gitter 0,01–0,4."""
    lj = np.linspace(math.log(0.01), math.log(0.4), 801); j = np.exp(lj)
    d = {}
    d["log-gleich 0,01–0,2"] = np.where(j <= 0.2 + 1e-12, 1.0, 0.0)
    for m in (0.05, 0.1):
        d["isotrop, Mittel %.2f" % m] = j * np.exp(-j / m)          # p(ε_j) dε_j = p·ε_j d ln ε_j
    return lj, {k: v / v.sum() for k, v in d.items()}


def form(IFj, lj, wts):
    """S(I) = Σ w·H(I_F(ε_j) − I); IFj: I_F auf dem ln-ε_j-Gitter."""
    return np.array([np.sum(wts * (IFj > I)) for I in IG])


def kennzahlen(S):
    Nm = np.trapezoid(S, IG)
    Ih = float(np.interp(0.5, S[::-1], IG[::-1]))
    I10 = float(np.interp(0.9, S[::-1], IG[::-1])); I90 = float(np.interp(0.1, S[::-1], IG[::-1]))
    return Nm, Ih, I10, I90


if __name__ == "__main__":
    jobs = [(eu / 3, ej) for eu in EPS for ej in EJ]
    import json
    cache = os.environ.get("IF_CACHE", "")           # optional: Zwischenspeicher der Faltungsrechnung
    if cache and os.path.exists(cache):
        RES = json.load(open(cache))
    else:
        with Pool(int(os.environ.get("KERNE", "4"))) as pool:
            RES = pool.map(_job, jobs)
        RES = [{k: float(v) for k, v in r.items()} for r in RES]
        if cache:
            json.dump(RES, open(cache, "w"))
    TAB = {(round(r["e"] * 3, 3), r["ej"]): r for r in RES}

    hdr("A  Erste Faltung I_F(ε, ε_j) [rad] (selbstähnlich, EdS; Verfahren uebergabe_faltung A)")
    P("   ε \\ ε_j   " + "  ".join("%6.2f" % j for j in EJ) + "    (Konvergenz max)")
    for eu in EPS:
        P("   %4.1f      " % eu + "  ".join("%6.3f" % TAB[(eu, j)]["I_F"] for j in EJ)
          + "    %.0e" % max(TAB[(eu, j)]["d"] for j in EJ))

    hdr("B  Rest S(I) über die Bahnen einer Schale: geglättete Stufe; effektives N*")
    lj, VER = verteilungen()
    lEJ = np.log(EJ)
    SFORM = {}
    P("   ε    Verteilung ε_j           N*_eff = ∫S dI   I_½     N*(Median) = I_½/ln2   Breite I(S=0,9)–I(S=0,1)")
    for eu in EPS:
        IFj = np.interp(lj, lEJ, [TAB[(eu, j)]["I_F"] for j in EJ])
        for k, w in VER.items():
            S = form(IFj, lj, w)
            SFORM[(eu, k)] = S
            Nm, Ih, I10, I90 = kennzahlen(S)
            P("   %.1f  %-24s %6.3f          %6.3f   %6.3f                 %.3f–%.3f" % (eu, k, Nm, Ih, Ih / math.log(2), I10, I90))
    GRUND = "log-gleich 0,01–0,2"
    P("   S(I) Grundform (log-gleich) gegen e^{−I/N*_eff}:")
    P("   I [rad]          " + " ".join("%6.2f" % x for x in (0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0)))
    for eu in (2.0, 3.0, 5.0):
        S = SFORM[(eu, GRUND)]; Nm = kennzahlen(S)[0]
        P("   ε = %.0f  Form   " % eu + " ".join("%6.3f" % np.interp(x, IG, S) for x in (0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0)))
        P("          e^{−I/%.2f}" % Nm + " ".join("%6.3f" % math.exp(-x / Nm) for x in (0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0)))

    # ========================================================================================== C  Linse
    import linsen_isoliert as L
    L.MB = [m * L.Ms for m in L.MB]
    G, Ms, kpc, Mpc, Gyr = L.G, L.Ms, L.kpc, L.Mpc, L.Gyr
    _orig = L.xi_R
    ESTUETZ = np.array(EPS)

    def S_von(I, eps, vert=GRUND):
        """S(I) bei ε (linear zwischen den Stützstellen, geklemmt auf 1–5)."""
        eps = np.clip(eps, ESTUETZ[0], ESTUETZ[-1])
        Sk = np.array([np.interp(I, IG, SFORM[(e_, vert)], right=0.0) for e_ in ESTUETZ])   # (nε, nI)
        out = np.empty_like(I)
        for i in range(len(I)):
            out[i] = np.interp(eps[i], ESTUETZ, Sk[:, i])
        return out

    def eps_profil(r, Mb, prof):
        if prof == "flach":
            return np.full_like(r, 2.0)
        Mf, _ = L.massen(Mb)
        a = np.array([(math.log(Mf[prof](x * 1.01)) - math.log(Mf[prof](x / 1.01))) / (2 * math.log(1.01)) for x in r])
        return 3 / a - 1

    def xi_form(r, Mb, prof, z, Ns):
        if isinstance(Ns, tuple):                   # ("Form", Verteilung)
            I = L.I_profil(Mb, prof, z)
            Ir = np.exp(np.interp(np.log(r), np.log(L.R_XI), np.log(np.maximum(I, 1e-300)), left=np.log(I[0]), right=np.log(I[-1])))
            return S_von(Ir, eps_profil(r, Mb, prof), Ns[1])
        return _orig(r, Mb, prof, z, Ns)
    L.xi_R = xi_form
    FORM = ("Form", GRUND)
    NEFF = {e_: kennzahlen(SFORM[(e_, GRUND)])[0] for e_ in EPS}

    hdr("C  Linse um isolierte Galaxien: ξ_R = S(I(r)), ε aus dem Profil (flach: 2; NFW: ε(r) = 3/α − 1)")
    Mmw = 6e10 * Ms
    RK = (100, 200, 300, 400, 500, 700, 1000)
    rr = np.array(RK) * kpc
    P("   ε(r) NFW (Milchstraße) bei %s kpc: %s" % ("/".join(map(str, RK)), " / ".join("%.1f" % x for x in eps_profil(rr, Mmw, "NFW"))))
    P("   Milchstraße 6·10¹⁰ M☉, z = 0, kumulativ, 3D: δg/g = g/g_Regel − 1 bei r [kpc] = %s" % " / ".join(map(str, RK)))
    for prof in ("flach", "NFW"):
        I = L.I_profil(Mmw, prof, 0.0)
        Ir = [math.exp(np.interp(math.log(x), np.log(L.R_XI), np.log(I))) for x in rr]
        P("     %-5s I(r)                " % prof + " ".join("%6.2f" % x for x in Ir))
        P("     %-5s ξ_R Form            " % prof + " ".join("%6.3f" % x for x in xi_form(rr, Mmw, prof, 0.0, FORM)))
        for nm, Ns in (("Form", FORM), ("e^{−I/%.2f}" % NEFF[2.0], NEFF[2.0]), ("N* = 0,3", 0.3), ("N* = 2,4", 2.4)):
            x = L.g_3D(Mmw, rr, Ns=Ns, prof=prof, z=0.0) / L.g_3D(Mmw, rr) - 1
            P("     %-5s %-18s " % (prof, nm) + " ".join("%+6.3f" % v for v in x))

    P()
    P("   Vorhersage 4GΔΣ/g_bar (Bin 3, M_b = %.2e M☉, z = 0,25, flach, kumulativ), 30 kpc–2 Mpc; KiDS [G] zum Vergleich" % (L.MB[2] / Ms))
    b = 2; Mb = L.MB[b]
    gbar = G * Mb / (L.R_MPC * Mpc)**2
    mods = {"Regel": L.g_4GdS(Mb, L.R_MPC * Mpc), "Form": L.g_4GdS(Mb, L.R_MPC * Mpc, Ns=FORM),
            "Form NFW": L.g_4GdS(Mb, L.R_MPC * Mpc, Ns=FORM, prof="NFW"),
            "exp N*eff": L.g_4GdS(Mb, L.R_MPC * Mpc, Ns=NEFF[2.0])}
    P("   R [kpc]   Regel   Form   Form NFW   e^{−I/N*eff}   KiDS [G] (16–84 %)")
    for i, R in enumerate(L.R_MPC):
        if R > 2.0:
            continue
        d = L.KIDS[b][i]
        P("   %6.0f   %6.2f  %6.2f   %6.2f     %6.2f        %6.2f (%.2f–%.2f)" % (R * 1e3, mods["Regel"][i] / gbar[i], mods["Form"][i] / gbar[i],
          mods["Form NFW"][i] / gbar[i], mods["exp N*eff"][i] / gbar[i], d[0] / gbar[i], d[1] / gbar[i], d[2] / gbar[i]))

    P()
    P("   χ² gegen KiDS [G] (4 Bins, Fenster R ≤ 0,30 / 1,04 / 2,60 Mpc; σ = halbe 16–84-Breite, Punkte unabhängig)")
    VAR = [("Regel", {}), ("Form flach", dict(Ns=FORM)), ("Form NFW", dict(Ns=FORM, prof="NFW"))]
    for k in VER:
        if k != GRUND:
            VAR.append(("Form flach, " + k, dict(Ns=("Form", k))))
    VAR += [("e^{−I/%.2f} flach" % NEFF[2.0], dict(Ns=NEFF[2.0])), ("N* = 0,3 flach", dict(Ns=0.3)), ("N* = 2,4 flach", dict(Ns=2.4))]
    for nm, kw in VAR:
        res = [L.chi2_kids(kw, rm) for rm in (0.3035, 1.0365, 2.6035)]
        P("     %-36s " % nm + "  ".join("%7.1f/%d" % (c_, n) for c_, n in res))
    P("   Formtest (je Bin eigene Massenskala f_M wie linsen_isoliert 2c):")
    for nm, kw in VAR[:3]:
        out = []
        for rm in (0.3035, 1.0365, 2.6035):
            ch, n, fs = L.chi2_frei(kw, rm)
            out.append("%6.1f/%d (" % (ch, n) + " ".join("%.2f" % f for f in fs) + ")")
        P("     %-20s " % nm + "  ".join(out))
    P("   (KiDS: Brouwer 2021 aus dem Repositorium in linsen_isoliert, Literaturwerte, nicht nachgeprüft.)")

    # ========================================================================================== D  Stoß
    hdr("D  Haufenstoß: Anteil der Operatorform beim Gas φ(t) = 1 − (1/x)∫_0^x S dI, x = ω_orb·t")
    import ftm_kern as Kk
    TT = (0.15, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)
    B = 0.169                                        # Bullet-Wert harvey_nachlauf
    for D_ in (200, 500):
        iw = Kk.uebergabezeit(D_, 0.3, (1, 1))[0] / Gyr
        P("   Δ = %d, z = 0,3: 1/ω_orb = %.3f Gyr;  t [Gyr] = %s" % (D_, iw, " / ".join("%.2f" % t for t in TT)))
        for eu in (2.0, 3.0):
            S = SFORM[(eu, GRUND)]; cS = np.concatenate([[0], np.cumsum(0.5 * (S[1:] + S[:-1]) * np.diff(IG))])
            ph = [1 - np.interp(t / iw, IG, cS) / (t / iw) for t in TT]
            tF = IG[np.argmax(S < 0.5)] * iw
            P("     Form ε = %.0f        φ: " % eu + " ".join("%.3f" % x for x in ph) + "   (φ > 0 ab %.2f Gyr, Hälfte der Bahnen gefaltet bei %.2f Gyr)"
              % (IG[np.argmax(S < 0.999)] * iw, tF))
        for Ns in (NEFF[3.0], 0.3, 2.4):
            ph = [1 - (-math.expm1(-t / iw / Ns)) / (t / iw / Ns) for t in TT]
            P("     e^{−I/%.2f}        φ: " % Ns + " ".join("%.3f" % x for x in ph))
    P("   Spät (x ≫ I_F): φ → 1 − N*_eff/x, gleich wie der Exponent mit N* = N*_eff; früh (x < I_F,min): φ = 0 statt x/(2N*).")
    P("   Überschuss gegen ΛCDM Δβ = B·φ (B = %.3f, Bullet-Wert harvey_nachlauf; Schranke Δβ ≤ 0,10 [G] ⇔ φ ≤ %.3f)." % (B, 0.10 / B))

    # ========================================================================================== E  Hand
    hdr("E  Handwerte")
    v = (1.193, 1.215, 1.349)                        # uebergabe_faltung.txt A, ε = 2, ε_j = 0,01 / 0,05 / 0,2
    a1, a2 = math.log(5), math.log(4)
    P("   ⟨I_F⟩ (ε = 2, log-gleich 0,01–0,2) aus uebergabe_faltung A (1,193 / 1,215 / 1,349 bei ε_j = 0,01 / 0,05 / 0,2), Trapez über ln ε_j:")
    P("     = [1,6094·(1,193 + 1,215)/2 + 1,3863·(1,215 + 1,349)/2]/2,9957 = %.4f   (Code, eigene Stützstellen, fein: %.4f)"
      % ((a1 * (v[0] + v[1]) / 2 + a2 * (v[1] + v[2]) / 2) / (a1 + a2), NEFF[2.0]))
    P("   Stufe bei I_F: φ(x) = 1 − I_F/x; x = 2·I_F ⇒ φ = 0,5; Code mit S = Stufe bei 1,2: φ(2,4) = %.3f"
      % (1 - np.trapezoid((IG < 1.2).astype(float) * (IG <= 2.4), IG) / 2.4))

    with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
