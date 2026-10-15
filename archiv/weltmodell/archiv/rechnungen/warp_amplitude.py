# -*- coding: utf-8 -*-
"""
Warp der Milchstraße: Präzession ω_p und Amplitudenänderung Ȧ/A unter der örtlichen Regel a = ν(|g_N|)·g_N.

1  Kinematik der Messgröße: v_z = ∂Z/∂t + Ω ∂Z/∂φ bei Z = A(R,t) sin(φ − φ_w(t)).
2  Welches Ȧ/A bräuchte die freie (rückläufige) Präzession, um die gemessene Rate vorzutäuschen?
3  Kippsteifigkeit der örtlichen Regel: gekoppelte Ringe (Selbstgravitation der Scheibe), Eigenmoden
   und Zeitentwicklung eines anfangs ruhenden Warps.
4  Blickrichtung: Ruhesystem der Galaxie, Eigenbewegung der Sonne.
5  Proben.

Linear um die ebene Scheibe: In der Mittelebene ist g_N,z = 0, |g_N| ändert sich erst in 2. Ordnung.
Die senkrechte Störbeschleunigung ist daher ν(g_R)·δg_N,z (örtliche Regel, am Ort des Sterns).
Ring k (Radius R_k, Kreisfrequenz Ω_k), Auslenkung z = Re[Z_k e^{iφ}] (m = 1), längs der Bahn d/dt = ∂_t + iΩ_k:
    (∂_t + iΩ_k)² Z_k = −ν_k Σ_l K_kl Z_l,   K_kl = δ_kl (Σ_j A_kj + GM_B(R_k)/R_k³) − B_kl,
    A_kl = (G m_l/2π)∮dψ/d³,  B_kl = (G m_l/2π)∮cos ψ dψ/d³,  d² = R_k² + R_l² − 2R_kR_l cos ψ + ε².
Starre Kippung aller Ringe (Z_k = R_k) gibt ν_k(KZ)_k = ν_k g_R,k = Ω_k² R_k: Frequenz 0 (Probe der Kovarianz der Regel).
Vorzeichen: Z ∝ e^{−iωt}, ω > 0 rechtläufig (im Drehsinn).  Einheiten kpc, km/s, M☉; 1 km/s/kpc = 1,0227/Gyr.
"""
import math
import numpy as np
from scipy.linalg import eig, expm
from ftm_kern import G, M_sun, pc, A0, A0_LIT, nu

kpc = 1e3 * pc
GK = G * M_sun / (kpc * 1e6)            # G in kpc (km/s)²/M☉
A0K = A0 * kpc / 1e6                    # a₀ in (km/s)²/kpc
GYR_PER = 1e9 * 3.15576e7 * 1e3 / kpc   # 1 km/s/kpc in 1/Gyr
MYR = lambda w: 1e3 / (abs(w) * GYR_PER)   # Zeitskala 1/w in Myr

# ------------------------------------------------------------------ Massenmodell (McMillan 2017 [G], nur Baryonen)
def sigma_disk(R):
    thin = 896.0e6 * np.exp(-R / 2.50)
    thick = 183.0e6 * np.exp(-R / 3.02)
    HI = 53.1e6 * np.exp(-4.0 / R - R / 7.0)
    H2 = 2180.0e6 * np.exp(-12.0 / R - R / 1.5)
    return thin + thick + HI + H2          # M☉/kpc²


def bulge_M_of_r():
    """Bulge kugelgemittelt (Abplattung 0,5 nicht berücksichtigt), Gesamtmasse 9,23·10⁹ M☉."""
    r = np.linspace(1e-4, 60, 60001)
    rho = (1 + r / 0.075) ** -1.8 * np.exp(-(r / 2.1) ** 2)
    m = np.concatenate([[0], np.cumsum(0.5 * (rho[1:] * r[1:]**2 + rho[:-1] * r[:-1]**2) * np.diff(r))])
    m *= 9.23e9 / m[-1]
    return lambda x: np.interp(x, r, m)


MB = bulge_M_of_r()


class Scheibe:
    def __init__(self, eps=0.3, dR=0.25, Rmax=40.0, a0=A0K, regel=True, npsi=4096):
        self.R = np.arange(dR / 2, Rmax, dR)
        R = self.R
        self.m = 2 * math.pi * R * dR * sigma_disk(R)
        psi = (np.arange(npsi) + 0.5) * 2 * math.pi / npsi
        dpsi = 2 * math.pi / npsi
        A = np.zeros((len(R), len(R)))
        B = np.zeros_like(A)
        cp = np.cos(psi)
        for k, Rk in enumerate(R):
            d3 = (Rk**2 + R[:, None]**2 - 2 * Rk * R[:, None] * cp[None, :] + eps**2) ** -1.5
            A[k] = GK * self.m / (2 * math.pi) * d3.sum(1) * dpsi
            B[k] = GK * self.m / (2 * math.pi) * (d3 * cp[None, :]).sum(1) * dpsi
        self.Ab = GK * MB(R) / R**3
        self.K = np.diag(A.sum(1) + self.Ab) - B
        self.gR = A.sum(1) * R - B @ R + self.Ab * R       # Newton in der Ebene (gleicher Kern)
        self.nuk = nu(self.gR * 1e6 / kpc, a0 * 1e6 / kpc) if regel else np.ones_like(R)
        self.Om = np.sqrt(self.nuk * self.gR / R)
        self.nz_star = np.sqrt(self.nuk * (A.sum(1) + self.Ab))     # Stern ohne Mitgehen der Nachbarn
        self.nz_ring = np.sqrt(self.nuk * (A.sum(1) + self.Ab - np.diag(B)))  # Ring als Ganzes, Nachbarn eben
        N = len(R)
        Om = np.diag(self.Om)
        self.L = np.block([[np.zeros((N, N)), np.eye(N)],
                           [Om @ Om - np.diag(self.nuk) @ self.K, -2j * Om]])

    def idx(self, r):
        return int(np.argmin(abs(self.R - r)))

    def moden(self):
        lam, V = eig(self.L)
        om = 1j * lam                      # Z ∝ e^{λt} = e^{−iωt}
        N = len(self.R)
        return om, V[:N]

    def entwickle(self, Z0, ts):
        N = len(self.R)
        y0 = np.concatenate([Z0.astype(complex), np.zeros(N, complex)])   # anfangs ruhende Warpfläche
        dt = ts[1] - ts[0]
        P = expm(self.L * dt)
        out = [y0[:N]]
        y = y0
        for _ in ts[1:]:
            y = P @ y
            out.append(y[:N])
        return np.array(out)


def linie(t="="):
    print(t * 78)


if __name__ == "__main__":
    linie()
    print("1  Quellen (Suchergebnisse; Volltexte arXiv/Nature/A&A/MNRAS vom Netz gesperrt)")
    linie()
    print("""   Poggio et al. 2020 (NatAs 4, 590): 12·10⁶ Riesensterne Gaia DR2; ω = +10,86 ± 0,03 (stat) ± 3,20 (sys)
     km/s/kpc; Mittel über vier geometrische Warp-Modelle aus der Literatur; mittlere v_z-Karte, v_R = 0,
     eine Rate für alle R; Amplitude zeitlich fest (∂h/∂t = 0) [G, Methoden nicht eingesehen].
   Cheng et al. 2020 (ApJ 905, 49): Gaia DR2 + APOGEE, ω = 13,57 (+0,20/−0,18) km/s/kpc; jüngere Sterne
     (3–6 Gyr) stärkere v_z-Störung als 6–9 Gyr.
   Chrobáková & López-Corredoira 2021 (ApJ 912, 130): Amplitude je Sternalter (γ, γ̇ = −dγ/d Alter);
     ω = 4 (+6/−4) km/s/kpc, verträglich mit 0: mit altersabhängiger Amplitude ist keine Präzession nötig.
   Chrobáková et al. 2022 (A&A, Überriesen EDR3): Warp asymmetrisch, z_w = +0,66 / −0,72 kpc bei R = 19,5–20 kpc.
   Dehnen, Semczuk & Schönrich 2023 (MNRAS 523, 1556), Cepheiden: Neigung ≈ 3° ab 14 kpc, rechtläufig
     ≈ 12 km/s/kpc bei 12 kpc, ≈ 6 ab 14 kpc; Knotenlinie verdreht sich um ≈ 14°/kpc im Drehsinn
     (−15° bei 11–12 kpc bis +45° bei 16,5 kpc).
   Cabrera-Gadea et al. 2024 (MNRAS 528, 4409), Cepheiden Z und v_z gemeinsam: Z_w = 0,00019 R^3,08
     sin(φ − (3,87R − 41,79 + 4,86 t)) kpc, ω = 4,86 ± 0,88 (stat) ± 2,14 (sys) km/s/kpc; Unsymmetrie (m = 2)
     wichtig, Amplitudenänderung nicht.""")

    # ------------------------------------------------------------------ 3 Ringe
    linie()
    print("3  Kippsteifigkeit der örtlichen Regel: gekoppelte Ringe (a₀ = %.2e m/s²)" % A0)
    linie()
    S = Scheibe()
    SN = Scheibe(regel=False)
    print("   Ringe 0,125–39,9 kpc, ΔR = 0,25 kpc, Glättung ε = 0,3 kpc; Bulge kugelgemittelt")
    print("   R    g_N/a0   ν      v[km/s]  Ω      ν_z,Stern  Ω−ν_z(Stern)  Ω−ν_z(Ring)  Newton: Ω  Ω−ν_z(Stern)")
    for r in (8.125, 12.125, 14.125, 16.125, 20.125):
        k = S.idx(r)
        print("  %5.1f  %6.3f  %5.3f  %6.1f  %6.2f  %7.2f    %8.2f      %8.2f     %6.2f   %8.2f" % (
            r, S.gR[k] / A0K, S.nuk[k], S.Om[k] * r, S.Om[k], S.nz_star[k], S.Om[k] - S.nz_star[k],
            S.Om[k] - S.nz_ring[k], SN.Om[k], SN.Om[k] - SN.nz_star[k]))

    # Probe starre Kippung
    Zr = S.R.astype(complex)
    res = (np.diag(S.nuk) @ S.K @ Zr) - S.Om**2 * S.R
    print("\n   Probe starre Kippung aller Ringe: max|ν(KZ) − Ω²R|/Ω²R = %.1e (Soll 0, Frequenz 0)" %
          np.max(abs(res) / (S.Om**2 * S.R)))

    # Eigenmoden
    for lab, Sx in (("örtliche Regel", S), ("Newton (nur Baryonen)", SN)):
        om, V = Sx.moden()
        pk = np.argmax(abs(V), axis=0)
        Ompk = Sx.Om[pk]
        slow = om.real < 0.5 * Ompk
        o = om[slow]; p = pk[slow]
        nz = abs(o) > 1e-3
        print("\n   Eigenmoden, %s: %d langsame Zweigmoden (Re ω < Ω/2 am Maximum)" % (lab, slow.sum()))
        print("     davon Re ω > +0,05 km/s/kpc (rechtläufig, ohne starre Kippung): %d; größtes Re ω = %+.3f" % (
            np.sum(o[nz].real > 0.05), o[nz].real.max()))
        print("     größtes |Im ω| (Wachsen/Abklingen) im langsamen Zweig: %.2e km/s/kpc; Mode ω≈0: %+.2e" % (
            abs(o.imag).max(), o[np.argmin(abs(o))].real))
        sel = (Sx.R[p] > 11.5) & (Sx.R[p] < 17)
        ss = np.sort(o[sel].real)[::-1]
        print("     Moden mit Maximum bei 11,5–17 kpc: %d, Re ω von %+.2f bis %+.2f km/s/kpc" % (len(ss), ss.min(), ss.max()))
        # glatteste Moden: wenige Nulldurchgänge
        zc = np.array([np.sum(np.diff(np.sign(np.real(V[:, j] * np.exp(-1j * np.angle(V[np.argmax(abs(V[:, j])), j]))))) != 0)
                       for j in range(V.shape[1])])[slow]
        print("     glatteste langsame Moden (Knoten in R, Re ω, Maximum):",
              ", ".join("(%d, %+.2f, %.1f kpc)" % (zc[j], o[j].real, Sx.R[p[j]]) for j in np.argsort(zc)[:6]))

    # Zeitentwicklung eines anfangs ruhenden Warps
    linie("-")
    print("   Zeitentwicklung, anfangs ruhende Warpfläche (∂_tZ = 0), Knotenrate ω = −d arg Z/dt, Ȧ/A = d ln|Z|/dt")
    ts = np.linspace(0, 1.0 * GYR_PER, 201)   # 0–1000 Myr, Schritt 5 Myr; t[kpc/(km/s)] = t[Gyr]·1,0227
    formen = {"Cabrera-Gadea 0,00019 R^3,08": lambda R: 0.00019 * R**3.08,
              "(R − 10)² ab 10 kpc": lambda R: np.clip(R - 10, 0, None)**2,
              "(R − 12)^1,5 ab 12 kpc": lambda R: np.clip(R - 12, 0, None)**1.5}
    S2, S5 = Scheibe(eps=0.2), Scheibe(eps=0.5)
    laeufe = [(n, f, "örtlich", S) for n, f in formen.items()] + [(n, f, "Newton", SN) for n, f in formen.items()] + \
             [("Cabrera-Gadea 0,00019 R^3,08", formen["Cabrera-Gadea 0,00019 R^3,08"], "örtl. ε0,2", S2),
              ("Cabrera-Gadea 0,00019 R^3,08", formen["Cabrera-Gadea 0,00019 R^3,08"], "örtl. ε0,5", S5)]
    print("   mittlere Knotenrate ω über 0–25 / 0–200 / 0–1000 Myr; Ȧ/A über 0–25 / 0–200 Myr [km/s/kpc]; momentan: Schritt 5 Myr")
    for name, f, lab, Sx in laeufe:
        Z = Sx.entwickle(f(Sx.R), ts)
        print("   %-28s %-10s" % (name, lab), end="")
        for r in (12.125, 14.125, 16.125):
            k = Sx.idx(r)
            ph = np.unwrap(np.angle(Z[:, k])); am = np.log(abs(Z[:, k]))
            w = [-(ph[i] - ph[0]) / (ts[i] - ts[0]) for i in (5, 40, 200)]
            a = [(am[i] - am[0]) / (ts[i] - ts[0]) for i in (5, 40)]
            print(" | R%2.0f: %+5.1f %+5.1f %+5.1f ; %+5.1f %+5.1f" % (r, w[0], w[1], w[2], a[0], a[1]), end="")
        k12, k16 = Sx.idx(12.125), Sx.idx(16.125)
        dph = np.degrees(np.unwrap(np.angle(Z[:, k12])) - np.unwrap(np.angle(Z[:, k16])))
        wmax = wgut = -1e9; amin = 1.0; anteil = []
        for r in (13.125, 14.125, 15.125, 16.125):
            kk = Sx.idx(r)
            w = -np.diff(np.unwrap(np.angle(Z[:, kk]))) / np.diff(ts)
            Am = 0.5 * (abs(Z[1:, kk]) + abs(Z[:-1, kk])) / abs(Z[0, kk])
            i = np.argmax(w)
            if w[i] > wmax:
                wmax, amax = w[i], Am[i]
            wgut = max(wgut, w[Am > 0.5].max()); amin = min(amin, Am.min())
            anteil.append(np.mean((w > 5) & (Am > 0.5)))
        print(" | Knoten 12/16 nach 1 Gyr %.0f° auseinander | 13–16 kpc in 1 Gyr: größte momentane Rate %+.1f bei A/A₀ = %.2f,"
              " mit A > A₀/2 höchstens %+.1f, Zeitanteil ω > 5 und A > A₀/2 %.0f %%, kleinstes A/A₀ %.2f" % (
            abs(dph[-1] - dph[0]), wmax, amax, wgut, 100 * max(anteil), amin))
    print("   (bei R = 12 kpc und Form (R−12)^1,5 ist A(0) ≈ 0: Phase dort nicht bestimmt)")

    # Empfindlichkeit
    linie("-")
    print("   Empfindlichkeit (Stern- und Ringwerte bei 14 kpc, Moden 11,5–17 kpc):")
    for lab, kw in (("ε = 0,2 kpc", dict(eps=0.2)), ("ε = 0,5 kpc", dict(eps=0.5)),
                    ("a₀ = 1,2e-10", dict(a0=A0_LIT * kpc / 1e6))):
        Sx = Scheibe(**kw)
        k = Sx.idx(14.125)
        om, V = Sx.moden()
        pk = np.argmax(abs(V), axis=0)
        slow = om.real < 0.5 * Sx.Om[pk]
        o = om[slow]; nzm = abs(o) > 1e-3
        print("   %-13s v(14) = %5.1f km/s, Ω = %5.2f, Ω−ν_z Stern %+6.2f, Ring %+6.2f; langsame Moden mit Re ω > 0,05: %d, max Re ω %+.3f, max|Im ω| %.1e" % (
            lab, Sx.Om[k] * Sx.R[k], Sx.Om[k], Sx.Om[k] - Sx.nz_star[k], Sx.Om[k] - Sx.nz_ring[k],
            np.sum(o[nzm].real > 0.05), o[nzm].real.max(), abs(o.imag).max()))

    # ------------------------------------------------------------------ 2 Kinematik
    linie()
    print("2  Kinematik: v_z = Ȧ sin(φ−φ_w) + A(Ω − ω_p) cos(φ−φ_w)")
    linie()
    k = S.idx(14.125)
    Om14 = S.Om[k]
    print("""   Am Knoten (φ = φ_w) fällt der Ȧ-Term heraus: v_z,Knoten = A(Ω − ω_p), unabhängig von Ȧ.
   Eine um den Knoten symmetrische Stichprobe mittelt den Ȧ-Term weg (ungerade in φ − φ_w).
   Liegt die Mitte der Stichprobe um δ neben dem Knoten, gibt eine Anpassung mit Ȧ = 0:
       ω_fit = ω_p − (Ȧ/A)·tan δ    →    Ȧ/A = (ω_p − ω_fit)/tan δ.""")
    print("   Ω(14 kpc) örtliche Regel = %.2f km/s/kpc; Umlaufzeit 2π/Ω = %.0f Myr" % (Om14, 2 * math.pi * MYR(Om14)))
    w_fit = 10.86
    cases = [("einzelner Stern Ω−ν_z", S.Om[k] - S.nz_star[k]), ("einzelner Ring Ω−ν_z", S.Om[k] - S.nz_ring[k]),
             ("Ring Feldform (38)", -13.98), ("starrer Außenteil Feldform (38)", -1.19), ("statischer Warp", 0.0)]
    print("   benötigtes Ȧ/A [km/s/kpc] (e-Zeit in Myr) für ω_fit = +10,86:")
    print("   %-32s %8s |" % ("freie Präzession", "ω_p") + "".join("  δ = %2d°        " % d for d in (5, 15, 30, 45)))
    for lab, wp in cases:
        print("   %-32s %+8.2f |" % (lab, wp) + "".join(
            "  %+7.1f (%4.0f)" % ((wp - w_fit) / math.tan(math.radians(d)), MYR((wp - w_fit) / math.tan(math.radians(d))))
            for d in (5, 15, 30, 45)))
    print("\n   Alternative Lesart (Amplitude der Riesensterne ≠ geometrische Amplitude): v_z,Knoten gleich ⇒")
    print("   A_kin/A_geo = (Ω − ω_fit)/(Ω − ω_p) bei Ω = %.2f:" % Om14)
    for lab, wp in cases:
        print("     %-32s %.3f" % (lab, (Om14 - w_fit) / (Om14 - wp)))

    # ------------------------------------------------------------------ 4 Blickrichtung
    linie()
    print("4  Blickrichtung (Ruhesystem der Galaxie)")
    linie()
    print("""   v_z = ∂Z/∂t + Ω∂Z/∂φ ist Geometrie im Ruhesystem der Galaxie; die Regel wirkt nur über ω_p und Ȧ.
   Der Beobachter geht über die Eigenbewegung der Sonne W☉ ein (≈ 7,25 km/s [G], Schönrich et al. 2010)
   und über Ω aus gemessenem v_φ.  Δω_fit = ΔW☉/A:""")
    for A in (0.5, 1.0):
        for dW in (0.5, 1.0):
            print("     A = %.1f kpc, ΔW☉ = %.1f km/s → Δω = %.1f km/s/kpc" % (A, dW, dW / A))
    for lab, wp in cases[:4]:
        print("   für ω_p = %+6.2f statt +10,86 müsste v_z,Knoten um %.1f·A km/s (A in kpc) verschoben sein" % (wp, w_fit - wp))

    # ------------------------------------------------------------------ 5 Proben
    linie()
    print("5  Proben")
    linie()
    S12 = Scheibe(a0=A0_LIT * kpc / 1e6)
    SN12 = S12
    k = S12.idx(14.125)
    print("   Rotationskurve bei a₀ = 1,2e-10 (Ringkern) gegen Gitterlösung warp_praezession (Feldform):")
    for r, vref, vnref in ((8.125, 220.0, 177.4), (12.125, 216.2, 157.1), (14.125, 213.0, 147.1), (16.125, 210.1, 138.2)):
        kk = S12.idx(r)
        print("     R = %.2f: v_N = %.1f (Gitter %.1f), v örtlich = %.1f (Gitter Feldform %.1f) km/s" % (
            r, math.sqrt(S12.gR[kk] * r), vnref, S12.Om[kk] * r, vref))
    gN14 = S.gR[S.idx(14.125)] / A0K
    s = math.sqrt(gN14); qq = math.exp(-s)
    print("   Hand: g_N/a₀ = %.4f → s = %.4f, q = e^{−s} = %.4f, ν = 1/(1−q) = %.4f (Kern %.4f)" % (gN14, s, qq, 1 / (1 - qq), S.nuk[S.idx(14.125)]))
    wp, d = S.Om[S.idx(14.125)] - S.nz_star[S.idx(14.125)], 15
    print("   Hand: Ȧ/A = (%.2f − 10,86)/tan 15° = %.1f/0,2679 = %.1f km/s/kpc; 1/(Ȧ/A) = %.1f Myr" % (
        wp, wp - 10.86, (wp - 10.86) / math.tan(math.radians(d)), MYR((wp - 10.86) / math.tan(math.radians(d)))))
