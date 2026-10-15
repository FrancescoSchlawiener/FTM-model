# -*- coding: utf-8 -*-
"""
Gedächtnis als Topologie: Bindung als Schließen einer Bahn im Phasenraum eines Bereichs.
Fragen: (1) ersetzt die Umlauf-/Windungszahl die kleinste Rate der Geschichte (Preisach, θ̃)?
        (2) ist N* die Zahl der Radiant bis zur Schließung (Fenster 0,3–2,4)?
        (3) ist ξ mit (∂ξ)² = −1 die Phase der Windung?
A  ein Bereich (Zykloide, Top-Hat): Phase η, Gewicht w(η), Umkehr als Kreuzung eines Strahls im Phasenraum
B  Monotonie der Phase, Wiederausdehnung, Aufbrechen der Schleife (Bahn mit Drehimpuls, feste Masse GM = 1)
C  Schalen ohne Stoß (Einfall um einen Keim, Einstein-de Sitter): Umkehr, erste Faltung des Phasenraumblatts, Radiant
D  Fenster: Stufe bei I_c statt e^{−I/N*} gegen die RAR-Enge (gewicht_preisach) und den Nachlauf nach Stößen
E  ξ: eindeutig bis zur ersten Faltung, danach mehrblättrig (Zahl der Ströme)
F  Handprüfungen
Konstanten und w aus ftm_kern.
"""
import math, os, sys, io, contextlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

OUT = []
def P(s=""):
    OUT.append(s); print(s)
def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)

TAU_D = c / A0                          # s
thc = 1 / TAU_D
X9, X1 = math.log(10)**2, math.log(10 / 9)**2     # θ/θ_c bei w = 0,9 und w = 0,1

# ============================================================================================== A
hdr("A  Ein Bereich als Zykloide: r = (R_ta/2)(1 − cos η), t = B(η − sin η), B = t_ta/π; Umkehr bei η = π")
P("   θ = 3ṙ/r = 3 sin η/(B(1 − cos η)²);  w = gewicht_w(θ) (Kern); Umkehr = Kreuzung des Strahls {v = 0, r > r_Gl} im Phasenraum")
P("   t_ta aus gleichung_zeit.txt (Kern-Hintergrund).")
P("   z_c   t_ta[Gyr]  η(w=0,9)   η(w=0,1)   Δη [rad]   Hand Δη      Δt [Mio. J.]   w bei η = π")
for zc, tta in ((3.4, 0.932), (2.0, 1.641), (1.0, 2.929), (0.5, 4.302)):
    B = tta * Gyr / math.pi
    eta = np.linspace(math.pi - 0.3, math.pi + 0.3, 600001)
    th = 3 * np.sin(eta) / (B * (1 - np.cos(eta))**2)
    w = -np.expm1(-np.sqrt(np.maximum(th, 0) / thc))
    e9 = eta[np.argmax(w < 0.9)]; e1 = eta[np.argmax(w < 0.1)]
    hand = (4 / 3) * B * (X9 - X1) * thc
    dt = B * (e1 - e9 - (math.sin(e1) - math.sin(e9)))
    P("   %4.1f   %6.3f     %.5f    %.5f    %.5f    %.5f      %6.1f        %.3f"
      % (zc, tta, e9, e1, e1 - e9, hand, dt / Gyr * 1e3, gewicht_w(0.0)))
P("   ⇒ w ist eine Stufe in der Phase des Bereichs: 1 vor, 0 nach der halben Windung η = π; Breite 0,01–0,1 rad.")
P()
P("   Radiant der Bahn ab der Umkehr, I = ∫ω dt, ω² = GM/r³ (feste Masse):")
f = lambda e: 1 / (math.sqrt(2) * math.sin(e / 2))          # dI/dη = ω·dt/dη = 1/√(1 − cos η)
def simpson(fn, a, b, n=200000):
    h = (b - a) / n; s = fn(a) + fn(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * fn(a + i * h)
    return s * h / 3
I_vir = simpson(f, math.pi, 1.5 * math.pi)
P("   Umkehr → r = R_ta/2 (η = 3π/2):  numerisch %.6f,  √2·ln(1 + √2) = %.6f" % (I_vir, math.sqrt(2) * math.log(1 + math.sqrt(2))))
P("   Umkehr → Kollaps mit ω geklemmt bei R_ta/2 (Konvention gewicht_preisach): %.4f + (π/2 − 1) = %.4f"
  % (I_vir, I_vir + math.pi / 2 - 1))
P("   Umkehr → R_ta/2 in ω_vir gezählt: ω_vir·Δt = π/2 + 1 = %.4f;  bis zum Kollaps ohne Klemme: ∫ divergiert (log)" % (math.pi / 2 + 1))
P("   Eine volle Windung (Umkehr → nächste Umkehr) im virialisierten Bereich: 2π = %.3f rad" % (2 * math.pi))

# ============================================================================================== B
hdr("B  Monotonie der Phase: Bahn mit Drehimpuls um r_Gl, Phase φ = arg(r − r_Gl, −v/ω_Gl); Umkehr = φ kreuzt 0 mod 2π")
P("   Analytisch: dφ/dt = [x·F(x)/ω + ω·y²]/(x² + y²), x = r − r_Gl, y = v/ω, F Rückstellkraft (x·F > 0) ⇒ dφ/dt > 0:")
P("   Die Phase läuft nur vorwärts, jede Umkehr kreuzt den Strahl im selben Sinn; rückwärts nur, wo die Kraft bei v = 0 nach außen zeigt.")


def bahn(szen, T=60.0, dt=2e-4):
    # GM = 1 anfangs; Start innen mit v > 0 (Ausdehnung vor der Umkehr); j² so, dass r_Gl = 1
    GM0, j2 = 1.0, 1.0
    r, v = 0.6, 0.5                     # innen, v > 0, E = −0,153 (gebunden)
    ts, rs, vs, Ms = [], [], [], []
    t = 0.0
    def GM(t):
        if szen == "fest":
            return GM0
        if szen == "langsam":            # 30 % Masseverlust über t = 20…35 (Wiederausdehnung, gebunden)
            return GM0 * (1 - 0.3 * min(max((t - 20) / 15, 0), 1))
        if szen == "plötzlich":          # 90 % Masseverlust bei t = 25: Schleife bricht auf
            return GM0 * (0.1 if t > 25 else 1.0)
    def acc(r, t):
        return -GM(t) / r**2 + j2 / r**3
    a = acc(r, t)
    while t < T:
        v += 0.5 * dt * a; r += dt * v; t += dt; a = acc(r, t); v += 0.5 * dt * a
        ts.append(t); rs.append(r); vs.append(v); Ms.append(GM(t))
    return np.array(ts), np.array(rs), np.array(vs), np.array(Ms), j2


P("   Szenario      Umkehrungen n   min dφ/dt   Energie am Ende   Schleife am Ende   w Preisach   w(θ jetzt) am Ende   max w(θ) nach 1. Umkehr")
for szen in ("fest", "langsam", "plötzlich"):
    ts, rs, vs, Ms, j2 = bahn(szen)
    rgl = j2 / Ms                        # Gleichgewicht je Zeit
    om = np.sqrt(Ms / rgl**3)
    x = rs - rgl; y = vs / om
    phi = np.unwrap(np.arctan2(-y, x))
    dphi = np.diff(phi) / np.diff(ts)
    dphi = dphi[np.diff(Ms) == 0]         # Sprung von r_Gl beim plötzlichen Verlust ausgenommen
    # Umkehrungen = Kreuzungen v: + → − (Strahl y = 0, x > 0), gezählt mit Vorzeichen
    sg = np.sign(vs)
    n_plus = int(np.sum((sg[:-1] > 0) & (sg[1:] < 0)))
    n_minus = int(np.sum((sg[:-1] < 0) & (sg[1:] > 0) & (x[1:] > 0)))
    E = 0.5 * vs[-1]**2 + j2 / (2 * rs[-1]**2) - Ms[-1] / rs[-1]
    th = 3 * vs / rs                     # Rate in Einheiten der Bahn; θ_c ≪ ω: Schwelle 10⁻³
    thcb = 1e-3
    wq = -np.expm1(-np.sqrt(np.maximum(th, 0) / thcb))
    wp = -np.expm1(-np.sqrt(max(np.min(th), 0) / thcb))
    i1 = np.argmax(vs < 0)
    P("   %-12s  %3d (−%d)        %.3f       %+.3f            %-16s   %.3f        %.3f                %.3f"
      % (szen, n_plus, n_minus, dphi.min(), E, "geschlossen" if E < 0 else "offen", wp, wq[-1], wq[i1:].max()))
P("   ⇒ Die gezählten Umkehrungen (Windungen der Geschichte) fallen nie; Preisach w = 0 ⇔ n ≥ 1, auch wenn die Schleife selbst aufbricht.")
P("     Geschützt ist die Zählung (ein Umlaufsinn, dφ/dt > 0), nicht die momentane Schleife: eine geschlossene Bahn kann über die Separatrix öffnen.")

# ============================================================================================== C
hdr("C  Schalen ohne Stoß um einen Keim (EdS, G = 1, t_i = 1): Umkehr, erste Faltung des Phasenraumblatts am Ort der Schale")
P("   Radiale Bahnen; ab der eigenen Umkehr Drehimpuls j² = ε·G·M(<r)·r (Perizentrum ≈ ε/2 … ε·R_ta); innerste Schalen (t_ta < t/60) als Zentralmasse.")
P("   Faltung = erste Änderung des Rangs der Schale (ein zweiter Strom erreicht ihren Ort, Zahl der Ströme 1 → 3).")
P("   I = ∫ω dt ab der Umkehr, ω² = G·M(<r)/r³ am Ort der Schale;  I_½ mit ω geklemmt bei R_ta/2 und M = M(<R_ta) (Konvention A).")


def schalen(eps, N=3000, t_end=1300.0, ziele=(8.0, 12.0, 18.0, 26.0)):
    rho_i = 1 / (6 * math.pi)
    Mtot = 0.25
    dM = Mtot / 150
    m = (Mtot - 0) / N
    Mlag = dM + m * (np.arange(N) + 0.5)                   # eingeschlossene Masse bis zur Schale (Lagrange)
    r = (3 * (Mlag - dM) / (4 * math.pi * rho_i)) ** (1 / 3)
    r = np.maximum(r, 1e-6)
    v = (2 / 3) * r                                        # Hubble-Fluss bei t = 1
    j2 = np.zeros(N)
    aktiv = np.ones(N, bool)
    Mc = dM
    t = 1.0
    tta = np.full(N, np.nan); rta = np.full(N, np.nan); gedreht = np.zeros(N, bool)
    zi = [int(np.argmin(abs(Mlag / dM - z))) for z in ziele]
    rec = {k: dict(I=0.0, I_half=0.0, tfold=None, Ifold=None, Ihalf_fold=None, rfold=None, peri=False, streams=[],
                   Iperi=None, Iapo2=None, napo=0, vprev=1.0) for k in zi}

    def beschl(r, t):
        idx = np.argsort(r[aktiv])
        ra = r[aktiv]
        Min = np.empty(ra.size)
        Min[idx] = Mc + m * np.arange(ra.size)             # Masse innen (ohne sich selbst)
        a = np.zeros(N)
        a[aktiv] = -Min / ra**2 + j2[aktiv] / ra**3
        M_all = np.zeros(N); M_all[aktiv] = Min
        rank = np.full(N, -1); rank[np.where(aktiv)[0][idx]] = np.arange(ra.size)
        return a, M_all, rank

    a, Mall, rank = beschl(r, t)
    rank0 = {k: None for k in zi}
    while t < t_end:
        ra = r[aktiv]; va = v[aktiv]
        om = np.sqrt((Mall[aktiv] + 1e-30) / ra**3)
        dt = min(0.02 * np.min(1 / np.maximum(om, 1e-12)), 0.02 * np.min(ra / np.maximum(abs(va), 1e-12)), 0.002 * t)
        v[aktiv] += 0.5 * dt * a[aktiv]
        r[aktiv] += dt * v[aktiv]
        neg = aktiv & (r < 0)
        r[neg] = -r[neg]; v[neg] = -v[neg]
        t += dt
        a, Mall, rank = beschl(r, t)
        v[aktiv] += 0.5 * dt * a[aktiv]
        # Umkehr: erster Vorzeichenwechsel v: + → −, dann Drehimpuls
        neu = aktiv & (~gedreht) & (v < 0)
        if np.any(neu):
            gedreht[neu] = True; tta[neu] = t; rta[neu] = r[neu]
            j2[neu] = eps * Mall[neu] * r[neu]
        # innerste Schalen einfrieren (Masse ins Zentrum), wenn sie nicht zu Zielen gehören
        alt = aktiv & gedreht & (tta < t / 60)
        for k in zi:
            alt[k] = False
        if np.any(alt):
            Mc += m * np.sum(alt); aktiv[alt] = False
            a, Mall, rank = beschl(r, t)
        for k in zi:
            d = rec[k]
            if rank0[k] is None:
                rank0[k] = rank[k]
            # Rang relativ: Zahl innerer aktiver Schalen mit kleinerem Lagrange-Index
            if gedreht[k]:
                om_k = math.sqrt(Mall[k] / r[k]**3)
                Mta = d.setdefault("Mta", Mall[k])
                om_h = math.sqrt(Mta / max(r[k], 0.5 * rta[k])**3)
                d["I"] += om_k * dt; d["I_half"] += om_h * dt
                if d["vprev"] < 0 and v[k] >= 0 and d["Iperi"] is None:
                    d["Iperi"] = d["I"]
                if d["vprev"] > 0 and v[k] <= 0 and d["Iperi"] is not None and d["Iapo2"] is None:
                    d["Iapo2"] = d["I"]
                d["vprev"] = v[k]
            if d["tfold"] is None:
                innen_lag = np.sum(aktiv[:k])
                if rank[k] != innen_lag:
                    d["tfold"] = t; d["Ifold"] = d["I"]; d["Ihalf_fold"] = d["I_half"]
                    d["rfold"] = r[k] / rta[k] if gedreht[k] else float("nan")
                    d["peri"] = d["Iperi"] is not None
                    d["tfold_rel"] = t / tta[k] if gedreht[k] else float("nan")
    return rec, zi, Mlag / dM, tta, rta, r.copy(), aktiv.copy()


for eps in (0.05, 0.1, 0.2):
    rec, zi, Mrel, tta, rta, r_end, akt_end = schalen(eps)
    LETZT = (zi, rta, r_end, akt_end)
    P("   ε = %.2f" % eps)
    P("     M/M_Keim   t_ta     t_Falt/t_ta   r_Falt/R_ta   Perizentrum vorher   I bis Faltung   I_½ bis Faltung   I bis 1. Perizentrum   I bis 2. Umkehr")
    for k in zi:
        d = rec[k]
        if d["tfold"] is None:
            P("     %6.1f   %7.2f   keine Faltung bis t_end" % (Mrel[k], tta[k])); continue
        if not (d["tfold_rel"] >= 1):
            P("     %6.1f   %7.2f   Faltung schon vor der eigenen Umkehr (Anlauf aus dem Hubble-Start), t = %.1f" % (Mrel[k], tta[k], d["tfold"])); continue
        P("     %6.1f   %7.2f      %.3f         %.3f          %-5s              %.3f            %.3f              %s                %s"
          % (Mrel[k], tta[k], d["tfold_rel"], d["rfold"], "ja" if d["peri"] else "nein", d["Ifold"], d["Ihalf_fold"],
             "%.3f" % d["Iperi"] if d["Iperi"] is not None else "–", "%.3f" % d["Iapo2"] if d["Iapo2"] is not None else "–"))

P("   Auflösung (ε = 0,2, t_end = 700): I bis Faltung je Schale M/M_Keim = 8 / 12 / 18 / 26")
for Nn in (1500, 6000):
    rr = schalen(0.2, N=Nn, t_end=700.0)
    P("     N = %d: %s" % (Nn, " / ".join("%.3f" % rr[0][k]["Ifold"] for k in rr[1])))
P("   ⇒ erste Faltung am Ort der Schale ≈ 0,9–1,2 rad nach der Umkehr (vor dem eigenen Perizentrum, bei r ≈ 0,6–0,7 R_ta);")
P("     halbe Windung (Perizentrum) ≈ 2,9–3,3 rad ≈ π, volle Windung (2. Umkehr) ≈ 5,6–6,4 rad ≈ 2π.")

# ============================================================================================== D
hdr("D  Fenster 0,3–2,4 bei Stufe statt Exponent: Übergabe vollständig, sobald I ≥ I_c (Schließung), vorher nicht")
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    import gewicht_preisach as GP
P("   gewicht_preisach (Konvention A, Preisach): I(t) = ∫(1 − w)ω dt, heute bzw. bei z_obs;  N*(RAR-Enge, Exponent) = %.3f" % GP.NRARp)
cf = max(abs(GP.ratio_rar(rk, 1.0, GP.res[0.5]["rv"]) - 1) for rk in (10, 20, 30, 50))
P("   volle Mode (R = 1) bei z_c = 0,5 heute: |g/g_RAR − 1| bis %.3f; Exponent verlangt e^{−I/N*} ≤ 0,05/%.3f" % (cf, cf))
P("   z_c     z_obs   I(t_obs) [rad]   I_c ≤ I? (1 | 1,25 | π/2 | 2,5 | π | 2π)")
cands = (1.0, math.sqrt(2) * math.log(1 + math.sqrt(2)), math.pi / 2, 2.5, math.pi, 2 * math.pi)
Imin = 1e9
for key, zo in (("Q", 2.0), (3.4, 2.0), (2.0, 1.0), (1.0, 0.5), (0.5, 0.0)):
    tobs = GP.t_of_a(1 / (1 + zo))
    I = GP.I_at(key, tobs, "p")
    Imin = min(Imin, I)
    zcs = GP.zcQ if key == "Q" else key
    P("   %5.2f   %4.1f    %7.3f          %s" % (zcs, zo, I, "  ".join("ja " if I >= ic else "NEIN" for ic in cands)))
P("   Stufe: RAR-Enge bei z_c = 0,5 heute verlangt I_c ≤ %.2f rad (Exponent: N* ≤ %.2f)" % (GP.I_at(0.5, GP.t0, "p"), GP.NRARp))
P("   (frisch kollabierte Bereiche bei z_obs: kleinstes I = %.2f; dort ist die Regel ohnehin mit R_ξ ≈ 0,06–0,14 gerechnet)" % Imin)
P()
P("   Nachlauf nach Stößen (harvey_nachlauf, B = 0,169, Δβ ≤ 0,10): Stufe gibt φ = 1 ab ω·t ≥ I_c ⇒ Δβ = 0,169 > 0,10;")
P("   nötig I_c > ω_orb·t für die Stoßalter t:")
for D_ in (200, 500):
    iw = uebergabezeit(D_, 0.3, (1, 1))[0] / Gyr
    P("     Δ = %d, z = 0,3: 1/ω_orb = %.3f Gyr;  t = 0,2 / 0,5 / 1,0 Gyr ⇒ I_c > %.2f / %.2f / %.2f rad" % (D_, iw, 0.2 / iw, 0.5 / iw, 1.0 / iw))
P("   Exponent zum Vergleich: N* ≥ 0,19–0,29 (t = 0,5 Gyr).")

# ============================================================================================== E
hdr("E  ξ mit (∂ξ)² = −1: Einheitstakt längs der Strömung, eindeutig nur bis zur ersten Faltung")
P("   (∂ξ)² = −1 ⇔ ξ ist Eigenzeit einer druckfreien Kongruenz (Flusslinien Geodäten). Eine Funktion längs der Flusslinien")
P("   ist auf dem Raum genau dann eindeutig, wenn die Lagrange-Abbildung injektiv ist; nach der Faltung hat ξ an jedem Ort so viele Werte")
P("   wie Ströme (1, 3, 5, …). Die Phase der Windung ist Φ = ∫ω dξ, (∂Φ)² = −ω² (ortsabhängig), also nicht ξ selbst.")
# Zahl der Ströme am Ort in C (letzter Lauf ε = 0,2): Lagrange-Abbildung r(M) am Ende
P("   Zahl der Ströme (Faltungen der Lagrange-Abbildung r(M)) im Lauf ε = 0,2 am Ende, je Radius in Einheiten von R_ta der Zielschale 26:")
zi, rta, r_end, akt_end = LETZT
ra = r_end[akt_end]
R0 = rta[zi[-1]]
for fr in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0):
    Rr = fr * R0
    ns = int(np.sum(np.diff(np.sign(ra - Rr)) != 0))
    P("     r = %.2f R_ta: %d Ströme" % (fr, ns))
P("   Die Zahl der Ströme ist ungerade (1 + 2·Faltungen) und wächst nach innen: sie zählt die Windungen des Phasenraumblatts.")

# ============================================================================================== F
hdr("F  Handprüfungen")
P("   √2·ln(1 + √2): ln(2,41421) = 0,881374; ·1,414214 = %.6f" % (1.414214 * 0.881374))
P("   Zykloide: dI/dη = ω·dt/dη; ω = √(GM/r³), r = (R_ta/2)(1 − cos η), dt/dη = B(1 − cos η), 1/B² = 8GM/R_ta³")
P("            ⇒ ω·dt/dη = (1 − cos η)^{−1/2} = 1/(√2 sin(η/2));  ∫_π^{3π/2} = √2·[ln tan(η/4)] = √2·ln tan(3π/8) = √2·ln(1 + √2)")
B2 = 2.0 / math.pi * Gyr / 1.0
P("   Breite z_c = 2 (t_ta = 1,641 Gyr): B = 0,5224 Gyr; Δη = (4/3)·B·(ln²10 − ln²(10/9))/τ_D = 1,3333·0,5224·5,2908/86,36 = %.4f rad"
  % (1.3333 * 0.5224 * 5.2908 / 86.36))
P("            Δt = 2B·Δη = %.1f Mio. J. (03: 45)" % (2 * 0.5224 * 1.3333 * 0.5224 * 5.2908 / 86.36 * 1e3))

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
