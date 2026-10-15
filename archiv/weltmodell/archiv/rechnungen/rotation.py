"""Forschung I.2: Kopplung an rotierende Massen.
Drehfaltung (Kerr in FTM-Größen), präzedierende Eigenfrequenzen (T³),
Siegel-Periodenmatrix Ω und Igusa-Spitzenform χ10, Folgerungen für Scheiben.
Nur numpy. Jede Identität wird mit [OK]/[FEHLER] geprüft."""
import sys, math, cmath
import numpy as np
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

G, c = 6.674e-11, 299792458.0
MS = 1.989e30
AU, PC, KPC = 1.496e11, 3.0857e16, 3.0857e19
YR = 3.15576e7
MAS = math.pi / (180 * 3600 * 1000)          # rad pro Millibogensekunde
KMSKPC = 1e3 / KPC                           # 1 km/s/kpc in 1/s
a0 = 1.2e-10

NCHK = [0, 0]
def chk(name, x, y, tol):
    rel = abs(x - y) / max(abs(y), 1e-300)
    ok = rel < tol
    NCHK[0] += 1
    NCHK[1] += (not ok)
    xs = f"{x:.10g}" if not isinstance(x, complex) else f"{x:.8g}"
    ys = f"{y:.10g}" if not isinstance(y, complex) else f"{y:.8g}"
    print(f"   [{'OK' if ok else 'FEHLER'}] {name}: {xs} vs {ys}  (rel. Abw. {rel:.1e})")

def ds_of(M): return G * M / (3 * c * c)     # Faltungstiefe Δs

# ---------------------------------------------------------------------------
# Kerr in geometrischen Einheiten: Länge in GM/c² = 3Δs, Spin χ = a/(3Δs)
# ---------------------------------------------------------------------------
def kerr(r, th, chi):
    s2, c2 = math.sin(th) ** 2, math.cos(th) ** 2
    Sig = r * r + chi * chi * c2
    Dk = r * r - 2 * r + chi * chi
    A = (r * r + chi * chi) ** 2 - chi * chi * Dk * s2
    gtt = -(1 - 2 * r / Sig)
    gtp = -2 * chi * r * s2 / Sig
    gpp = A * s2 / Sig
    return gtt, gtp, gpp, Sig / Dk, Sig, Sig, Dk, A

def horizons(chi):
    w = math.sqrt(1 - chi * chi)
    return 1 + w, 1 - w

def r_ergo(th, chi):
    return 1 + math.sqrt(1 - chi * chi * math.cos(th) ** 2)

def isco(chi, prograde=True):
    Z1 = 1 + (1 - chi * chi) ** (1 / 3) * ((1 + chi) ** (1 / 3) + (1 - chi) ** (1 / 3))
    Z2 = math.sqrt(3 * chi * chi + Z1 * Z1)
    s = math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    return 3 + Z2 - s if prograde else 3 + Z2 + s

def freqs(r, chi):
    """Okazaki-Kato-Fukue / Kato / Wagoner, prograd, Einheiten c³/(GM)."""
    Op = 1 / (r ** 1.5 + chi)
    Or2 = Op * Op * (1 - 6 / r + 8 * chi * r ** -1.5 - 3 * chi * chi / r ** 2)
    Oz2 = Op * Op * (1 - 4 * chi * r ** -1.5 + 3 * chi * chi / r ** 2)
    return Op, (math.sqrt(Or2) if Or2 > 0 else float("nan")), math.sqrt(Oz2)

print("=" * 100)
print("1) DREHFALTUNG — Kerr-Metrik in FTM-Größen (GM/c² = 3Δs, a = J/(Mc), χ = a/(3Δs))")
print("=" * 100)
print("   Horizonte r± = 3Δs ± sqrt(9Δs² − a²);  Ergofläche r_E(θ) = 3Δs + sqrt(9Δs² − a² cos²θ)")
for chi in (0.3, 0.9, 0.998):
    rp, rm = horizons(chi)
    chk(f"χ={chi}: Δ_K(r+) = 0  (r+ = {rp:.6f}·3Δs)", rp * rp - 2 * rp + chi * chi + 1, 1.0, 1e-13)
    chk(f"χ={chi}: Δ_K(r−) = 0  (r− = {rm:.6f}·3Δs)", rm * rm - 2 * rm + chi * chi + 1, 1.0, 1e-13)
    for thd in (30, 60, 90):
        th = math.radians(thd)
        chk(f"χ={chi}: g_tt(r_E(θ={thd}°)) = 0", kerr(r_ergo(th, chi), th, chi)[0] + 1, 1.0, 1e-13)
    chk(f"χ={chi}: r_E(Äquator) = 2·(3Δs) = 6Δs = r0", r_ergo(math.pi / 2, chi), 2.0, 1e-15)
    chk(f"χ={chi}: r_E(Pol) = r+", r_ergo(0.0, chi), rp, 1e-15)

print("\n   Takt ruhender Uhren (dτ/dt = sqrt(−g_tt)) am Äquator = FTM-Formel sqrt(1 − 6Δs/r), für jedes a:")
for chi in (0.5, 0.9):
    for r in (2.5, 4.0, 10.0):
        chk(f"χ={chi}, r={r}·3Δs: sqrt(−g_tt(θ=90°)) = sqrt(1−2/r)",
            math.sqrt(-kerr(r, math.pi / 2, chi)[0]), math.sqrt(1 - 2 / r), 1e-14)

print("\n   Mitnahme ω = −g_tφ/g_φφ (exakt) gegen schwaches Feld 2GJ/(c²r³) = 6cΔs·a/r³ (Einheiten: 2χ/r³):")
for r in (10.0, 100.0, 1000.0):
    gtt, gtp, gpp, *_ = kerr(r, math.pi / 2, 0.9)
    print(f"   r = {r:7.0f}·3Δs   ω_exakt/ω_schwach = {(-gtp / gpp) / (2 * 0.9 / r ** 3):.8f}")
chk("ω(r=1000)/ω_schwach → 1", (-kerr(1000.0, math.pi / 2, 0.9)[1] / kerr(1000.0, math.pi / 2, 0.9)[2]) / (1.8e-9), 1.0, 1e-5)

print("\n   Verwindung als Uhrzeit: T_φ = 2π(r^{3/2} ± χ)·GM/c³  ⇒  T_pro − T_retro = 4πa/c (unabhängig von r)")
for r in (6.0, 20.0, 1e4):
    chi = 0.9
    Tp = 2 * math.pi / freqs(r, chi)[0]
    Tr = 2 * math.pi * (r ** 1.5 - chi)          # retrograd: |Ω_φ| = 1/(r^{3/2} − χ)
    chk(f"r={r:g}: (T_pro − T_retro)/(GM/c³) = 4πχ", Tp - Tr, 4 * math.pi * chi, 1e-9)

bodies = [
    ("Erde", 5.972e24, 5.86e33, 6.371e6),
    ("Sonne", MS, 1.92e41, 6.957e8),
    ("Sgr A*", 4.30e6 * MS, None, None),
    ("M87*", 6.5e9 * MS, None, None),
    ("Milchstraße", 6e10 * MS, 2 * 6e10 * MS * 2.6 * KPC * 220e3, None),
]
print("\n   Körper: Faltungstiefe Δs, Verwindung a = J/(Mc), χ = a/(3Δs), GJ/c³ = 3Δs·a [m²], Uhreffekt 4πa/c")
print(f"   {'Körper':12s} {'Δs [m]':>11s} {'3Δs=GM/c²':>11s} {'J [kg m²/s]':>12s} {'a [m]':>11s} {'χ':>9s} {'GJ/c³ [m²]':>11s} {'4πa/c [s]':>11s}")
BODY = {}
for name, M, J, R in bodies:
    ds = ds_of(M)
    if J is None:
        chi = 0.9
        J = chi * G * M * M / c
    a = J / (M * c)
    chi = a / (3 * ds)
    BODY[name] = (M, J, R, ds, a, chi)
    print(f"   {name:12s} {ds:11.4e} {3*ds:11.4e} {J:12.4e} {a:11.4e} {chi:9.4g} {G*J/c**3:11.4e} {4*math.pi*a/c:11.4e}")
    chk(f"{name}: GJ/c³ = 3Δs·a", G * J / c ** 3, 3 * ds * a, 1e-12)
print("   Hinweis: GJ/c³ ist eine Fläche. Die Verwindungs-LÄNGE ist a = J/(Mc); a/Δs = 3χ.")

M, J, R, ds, a, chi = BODY["Erde"]
print(f"\n   Erde: Kerr-Quadrupol J2_Kerr = (a/R)² = {(a/R)**2:.3e}  gegen gemessen J2 = 1.0826e-3")
print("   -> Erde/Sonne sind keine Kerr-Körper; universell ist nur der Dipolanteil (Lense-Thirring), Ordnung a¹.")
M, J, R, ds, a, chi = BODY["Sonne"]
print(f"   Sonne: J2_Kerr = (a/R)² = {(a/R)**2:.3e}  gegen helioseismisch J2 ≈ 2.2e-7")

print("\n   Schwarze Löcher (χ = 0.9 angenommen):")
for name in ("Sgr A*", "M87*"):
    M, J, R, ds, a, chi = BODY[name]
    rg = 3 * ds
    rp, rm = horizons(chi)
    OmH = c * a / ((rp * rg) ** 2 + a * a)
    ri, ro = isco(chi, True), isco(chi, False)
    Op_isco = freqs(ri, chi)[0] * c ** 3 / (G * M)
    print(f"   {name}: 3Δs = {rg:.4e} m,  r+ = {rp*rg:.4e} m ({rp:.4f}·3Δs), r− = {rm*rg:.4e} m,"
          f" r_E(Äq) = 6Δs = {2*rg:.4e} m, r_E(Pol) = r+")
    print(f"          Ω_H = {OmH:.4e} rad/s (Periode {2*math.pi/OmH/60:.2f} min),"
          f" ISCO pro/retro = {ri:.4f} / {ro:.4f} ·3Δs,  f_φ(ISCO) = {Op_isco/(2*math.pi):.4e} Hz"
          f" (Periode {2*math.pi/Op_isco/60:.2f} min)")
    chk(f"{name}: κ²(ISCO) = 0 (BPT-ISCO ist Nullstelle der Radialfrequenz)",
        1 - 6 / ri + 8 * chi * ri ** -1.5 - 3 * chi * chi / ri ** 2 + 1, 1.0, 1e-12)

print("\n   Lense-Thirring und geodätische Präzession (schwaches Feld):")
M, J, R, ds, a, chi = BODY["Erde"]
GM = G * M
r = 7027.4e3
geo = 1.5 * GM ** 1.5 / (c * c * r ** 2.5)
fd = G * J / (2 * c * c * r ** 3)
chk("GP-B geodätisch: (3/2)(GM)^{3/2}/(c²r^{5/2}) = (9/2)(Δs/r)·Ω_Kepler",
    geo, 4.5 * (ds / r) * math.sqrt(GM / r ** 3), 1e-12)
print(f"   GP-B (r = 7027.4 km, polar): geodätisch       {geo*YR/MAS:8.1f} mas/Jahr  (gemessen 6601.8 ± 18.3)")
print(f"                                Mitnahme GJ/(2c²r³) {fd*YR/MAS:8.1f} mas/Jahr (Bahnmittel, Kreiselrichtung ⊥ Erdachse)")
print(f"                                × cos(16.84°) (Leitstern IM Peg) = {fd*math.cos(math.radians(16.84))*YR/MAS:6.1f} mas/Jahr"
      f"  (Vorhersage 39.2; gemessen 37.2 ± 7.2)")
chk("Mitnahme-Rate = 3cΔs·a/(2r³)", fd, 3 * c * ds * a / (2 * r ** 3), 1e-12)
for nm, aa, e in (("LAGEOS", 12270e3, 0.0045), ("LAGEOS 2", 12163e3, 0.0135)):
    nod = 2 * G * J / (c * c * aa ** 3 * (1 - e * e) ** 1.5)
    print(f"   {nm:9s} Knoten: 2GJ/(c²a³(1−e²)^{{3/2}}) = {nod*YR/MAS:6.2f} mas/Jahr")
print(f"   Gravitomagnetischer Uhreffekt Erde: T_pro − T_retro = 4πJ/(Mc²) = {4*math.pi*a/c*1e9:.1f} ns pro Umlauf")

M, J, R, ds, a, chi = BODY["Sonne"]
GM = G * M
am, e, Tm = 5.7909e10, 0.2056, 87.969 * 86400
norb = 100 * YR / Tm
per = 6 * math.pi * GM / (c * c * am * (1 - e * e))
chk("Merkur: 6πGM/(c²a(1−e²)) = 18πΔs/(a(1−e²))", per, 18 * math.pi * ds / (am * (1 - e * e)), 1e-12)
print(f"   Merkur Perihel (Schwarzschild): {per*norb/MAS/1000:.2f} ″/Jhdt")
nodM = 2 * G * J / (c * c * am ** 3 * (1 - e * e) ** 1.5)
print(f"   Merkur Lense-Thirring: Knoten {nodM*100*YR/MAS:.3f} mas/Jhdt, Perihellänge −2× = {-2*nodM*100*YR/MAS:.3f} mas/Jhdt")
print(f"   Sonnenoberfläche: ω = 2GJ/(c²R³) = {2*G*J/(c*c*R**3):.3e} rad/s = {2*G*J/(c*c*R**3)*YR/MAS:.1f} mas/Jahr")

M, J, R, ds, a, chi = BODY["Milchstraße"]
for Rk in (8.2, 15.0, 25.0):
    w = 2 * G * J / (c * c * (Rk * KPC) ** 3)
    print(f"   Milchstraße, r = {Rk:4.1f} kpc: 2GJ/(c²r³) = {w/KMSKPC:.3e} km/s/kpc  (Gaia-Warp-Präzession ≈ 10.9 km/s/kpc)")

print("\n   Takt rotierender Trichter, χ = 0.9 (r in 3Δs; dimensionslos, gilt für Sgr A* und M87*):")
print("   statisch: dτ/dt = sqrt(1 − 6Δs·r/(r²+a²cos²θ));  ZAMO: α = sqrt(ΣΔ_K/A)")
chi = 0.9
print(f"   {'r/3Δs':>7s} {'stat. Äquator':>14s} {'stat. Pol':>10s} {'ZAMO Äquator':>13s} {'ZAMO Pol':>10s} {'Pol/Äq (stat.)':>15s}")
for r in (1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 20.0):
    ge = kerr(r, math.pi / 2, chi); gp = kerr(r, 1e-9, chi)
    se = math.sqrt(-ge[0]) if ge[0] < 0 else float("nan")
    sp = math.sqrt(-gp[0]) if gp[0] < 0 else float("nan")
    Sig, Dk, A = ge[4], ge[6], ge[7]
    ze = math.sqrt(Sig * Dk / A) if Dk > 0 else float("nan")
    zp = math.sqrt(gp[4] * gp[6] / gp[7]) if gp[6] > 0 else float("nan")
    print(f"   {r:7.2f} {se:14.6f} {sp:10.6f} {ze:13.6f} {zp:10.6f} {sp/se if se==se else float('nan'):15.6f}")
rp = horizons(chi)[0]
chk("Pol: ZAMO-Lapse = statischer Takt (g_tφ = 0 auf der Achse)",
    math.sqrt(kerr(3.0, 1e-9, chi)[4] * kerr(3.0, 1e-9, chi)[6] / kerr(3.0, 1e-9, chi)[7]),
    math.sqrt(1 - 2 * 3.0 / (9 + chi * chi)), 1e-9)

print("\n   Takt an der Erdoberfläche (Äquator gegen Pol), mitrotierende Uhren:")
GMe, Re, Rp, J2, Om = 3.986004418e14, 6378137.0, 6356752.314, 1.08263e-3, 7.292115e-5
Weq = GMe / Re * (1 + J2 / 2) + 0.5 * Om * Om * Re * Re
Wpo = GMe / Rp * (1 - J2 * (Re / Rp) ** 2)
print(f"   W_Äq = {Weq:.1f} m²/s²   W_Pol = {Wpo:.1f} m²/s²   (Potential inkl. J2 und Fliehkraft)")
print(f"   1 − dτ/dt: Äquator {Weq/c**2:.6e}  Pol {Wpo/c**2:.6e}  Differenz {(Weq-Wpo)/c**2:.2e}")
print(f"   Einzelbeiträge: Fliehkraft Ω²R²/2c² = {0.5*Om*Om*Re*Re/c**2:.3e},  Abplattung GM(1/Rp−1/Re)/c² = {GMe*(1/Rp-1/Re)/c**2:.3e}")
Mz, Jz, Rz, dsz, az, chiz = BODY["Erde"]
print(f"   Kerr-eigene Beiträge: Mitnahme 6Δs·a·Ω/(cR) = {6*dsz*az*Om/(c*Re):.2e},  a²-Term 3Δs·a²/R³ = {3*dsz*az*az/Re**3:.2e}")
print("   -> Das Geoid ist Äquipotentialfläche: die Uhren gehen gleich (Rest = vernachlässigtes J4/Geoidhöhe).")

# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("2) PRÄZEDIERENDE EIGENFREQUENZEN — Ω_φ, Ω_r, Ω_θ (Kerr, äquatoriale Kreisbahn)")
print("=" * 100)
print("   Ω_φ = c³/(GM(r̃^{3/2}+χ)),  Ω_r² = Ω_φ²(1 − 6/r̃ + 8χ r̃^{-3/2} − 3χ²/r̃²),  Ω_θ² = Ω_φ²(1 − 4χ r̃^{-3/2} + 3χ²/r̃²)")

def dr_num(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)

def circ_Omega(r, chi):
    th = math.pi / 2; h = 1e-5 * r
    dgt = dr_num(lambda x: kerr(x, th, chi)[0], r, h)
    dgp = dr_num(lambda x: kerr(x, th, chi)[1], r, h)
    dgf = dr_num(lambda x: kerr(x, th, chi)[2], r, h)
    disc = dgp * dgp - dgf * dgt
    roots = [(-dgp + s * math.sqrt(disc)) / dgf for s in (1, -1)]
    return max(roots)

def EL(r, chi, Op):
    gtt, gtp, gpp, *_ = kerr(r, math.pi / 2, chi)
    ut = 1 / math.sqrt(-(gtt + 2 * gtp * Op + gpp * Op * Op))
    return -(gtt + gtp * Op) * ut, (gtp + gpp * Op) * ut, ut

def ginv(r, th, chi):
    gtt, gtp, gpp, grr, gthth, *_ = kerr(r, th, chi)
    det = gtt * gpp - gtp * gtp
    return gpp / det, -gtp / det, gtt / det, 1 / grr, 1 / gthth

def Vpot(r, th, chi, E, L):
    itt, itp, ipp, _, _ = ginv(r, th, chi)
    return 0.5 * (itt * E * E - 2 * itp * E * L + ipp * L * L)

def lin_freqs(r, chi):
    Op = circ_Omega(r, chi)
    E, L, ut = EL(r, chi, Op)
    hr, ht = 1e-4 * r, 1e-4
    th0 = math.pi / 2
    V0 = Vpot(r, th0, chi, E, L)
    Vrr = (Vpot(r + hr, th0, chi, E, L) - 2 * V0 + Vpot(r - hr, th0, chi, E, L)) / hr ** 2
    Vtt = (Vpot(r, th0 + ht, chi, E, L) - 2 * V0 + Vpot(r, th0 - ht, chi, E, L)) / ht ** 2
    _, _, _, irr, ithth = ginv(r, th0, chi)
    return Op, math.sqrt(irr * Vrr) / ut, math.sqrt(ithth * Vtt) / ut, V0

print("   Prüfung A: Ω_φ aus der Kreisbahnbedingung ∂_r g_tt + 2Ω ∂_r g_tφ + Ω² ∂_r g_φφ = 0 (numerisch)")
print("   Prüfung B: Ω_r, Ω_θ aus Linearisierung der Geodätengleichung (Hamilton, V'' numerisch, /u^t)")
for chi, r in ((0.0, 8.0), (0.5, 7.0), (0.9, 5.0), (0.9, 12.0), (-0.7, 12.0)):
    Opf, Orf, Ozf = freqs(r, chi) if chi >= 0 else (None, None, None)
    if chi < 0:   # retrograd = prograd mit χ → −χ in den Formeln
        Opf = 1 / (r ** 1.5 + chi)
        Orf = Opf * math.sqrt(1 - 6 / r + 8 * chi * r ** -1.5 - 3 * chi * chi / r ** 2)
        Ozf = Opf * math.sqrt(1 - 4 * chi * r ** -1.5 + 3 * chi * chi / r ** 2)
    Op, Or, Oz, V0 = lin_freqs(r, chi)
    chk(f"χ={chi:+.1f}, r̃={r}: Ω_φ", Op, Opf, 1e-8)
    chk(f"χ={chi:+.1f}, r̃={r}: Ω_r", Or, Orf, 1e-6)
    chk(f"χ={chi:+.1f}, r̃={r}: Ω_θ", Oz, Ozf, 1e-6)
    chk(f"χ={chi:+.1f}, r̃={r}: Normierung g^{{μν}}p_μp_ν = −1", 2 * V0, -1.0, 1e-10)

print("\n   Prüfung C: volle nichtlineare Geodäte (RK4), χ = 0.9, r̃ = 8, Auslenkung 1e-4, Frequenzen aus Nulldurchgängen")
def geod_rk4(r0, chi, dr, dth, lam_end, hst):
    Op = circ_Omega(r0, chi)
    E, L, _ = EL(r0, chi, Op)
    def H(r, th, pr, pth):
        itt, itp, ipp, irr, ith = ginv(r, th, chi)
        return 0.5 * (irr * pr * pr + ith * pth * pth + itt * E * E - 2 * itp * E * L + ipp * L * L)
    def rhs(s):
        t, r, th, ph, pr, pth = s
        itt, itp, ipp, irr, ith = ginv(r, th, chi)
        hr, ht = 1e-6 * r, 1e-6
        dHr = (H(r + hr, th, pr, pth) - H(r - hr, th, pr, pth)) / (2 * hr)
        dHt = (H(r, th + ht, pr, pth) - H(r, th - ht, pr, pth)) / (2 * ht)
        return np.array([-itt * E + itp * L, irr * pr, ith * pth, -itp * E + ipp * L, -dHr, -dHt])
    s = np.array([0.0, r0 + dr, math.pi / 2 + dth, 0.0, 0.0, 0.0])
    n = int(lam_end / hst)
    T = np.empty(n + 1); Rr = np.empty(n + 1); Th = np.empty(n + 1); Ph = np.empty(n + 1)
    T[0], Rr[0], Th[0], Ph[0] = s[0], s[1], s[2], s[3]
    for i in range(n):
        k1 = rhs(s); k2 = rhs(s + 0.5 * hst * k1); k3 = rhs(s + 0.5 * hst * k2); k4 = rhs(s + hst * k3)
        s = s + hst / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        T[i + 1], Rr[i + 1], Th[i + 1], Ph[i + 1] = s[0], s[1], s[2], s[3]
    return T, Rr, Th, Ph

def freq_from_crossings(T, X):
    x = X - 0.5 * (X.max() + X.min())
    idx = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = T[idx] - x[idx] * (T[idx + 1] - T[idx]) / (x[idx + 1] - x[idx])
    return 2 * math.pi * (len(tc) - 1) / (tc[-1] - tc[0])

r0, chi = 8.0, 0.9
T, Rr, Th, Ph = geod_rk4(r0, chi, 1e-4 * r0, 1e-4, 2600.0, 0.2)
Opf, Orf, Ozf = freqs(r0, chi)
Or_n = freq_from_crossings(T, Rr); Oz_n = freq_from_crossings(T, Th)
Op_n = (Ph[-1] - Ph[0]) / (T[-1] - T[0])
chk("RK4: Ω_φ (mittlere Winkelrate)", Op_n, Opf, 1e-4)
chk("RK4: Ω_r (Radialschwingung)", Or_n, Orf, 1e-4)
chk("RK4: Ω_θ (Vertikalschwingung)", Oz_n, Ozf, 1e-4)
print(f"   -> Bahn auf T³ mit (Ω_φ, Ω_r, Ω_θ) = ({Opf:.6f}, {Orf:.6f}, {Ozf:.6f}) c³/GM;"
      f" Knoten {Opf-Ozf:.6f}, Periapsis {Opf-Orf:.6f}")

print("\n   Grenzfälle:")
for r in (6.0, 10.0, 50.0):
    Op, Or, Oz = freqs(r, 0.0)
    chk(f"χ=0, r̃={r}: Ω_θ = Ω_φ (Entartung: Vertikalzyklus = Azimutzyklus)", Oz, Op, 1e-15)
for r in (1e3, 1e4, 1e5):
    Op, Or, Oz = freqs(r, 0.9)
    print(f"   r̃ = {r:.0e}: (Ω_φ−Ω_θ)/(2χ/r̃³) = {(Op-Oz)/(1.8/r**3):.6f}   (Ω_φ−Ω_r)/(3/r̃^{{5/2}}) = {(Op-Or)/(3/r**2.5):.6f}")
print("   Die Abweichung von 1 fällt wie r̃^{-1/2}: nächste Ordnung 1 − (3/4)χ r̃^{-1/2} (Knoten), 1 − (4/3)χ r̃^{-1/2} (Periapsis)")
for r in (1e4, 1e5, 1e6):
    Op, Or, Oz = freqs(r, 0.9)
    chk(f"r̃={r:.0e}: Knoten/(2χ/r̃³) = 1 − (3/4)χ/√r̃", (Op - Oz) / (1.8 / r ** 3), 1 - 0.75 * 0.9 / math.sqrt(r), 3 / r)
    chk(f"r̃={r:.0e}: Periapsis/(3/r̃^{{5/2}}) = 1 − (4/3)χ/√r̃", (Op - Or) / (3 / r ** 2.5), 1 - 4 / 3 * 0.9 / math.sqrt(r), 3 / r)

print("\n   Frequenzverhältnisse und Knotenzahl ζ = (Ω_φ−Ω_θ)/Ω_φ, χ = 0.9 (prograd) und χ = 0:")
print(f"   {'r̃':>7s} {'Ω_r/Ω_φ(0.9)':>13s} {'Ω_θ/Ω_φ(0.9)':>13s} {'ζ(0.9)':>11s} {'2χr̃^-3/2':>11s} {'Ω_r/Ω_φ(0)':>12s} {'Ω_θ/Ω_φ(0)':>12s}")
ri = isco(0.9)
for r in (ri, 3.0, 4.0, 6.0, 10.0, 20.0, 100.0):
    Op, Or, Oz = freqs(r, 0.9)
    Op0, Or0, Oz0 = freqs(r, 0.0)
    print(f"   {r:7.3f} {Or/Op:13.6f} {Oz/Op:13.6f} {1-Oz/Op:11.6f} {1.8*r**-1.5:11.6f} {Or0/Op0 if Or0==Or0 else float('nan'):12.6f} {Oz0/Op0:12.6f}")

def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return 0.5 * (lo + hi)
r32 = bisect(lambda r: freqs(r, 0.9)[2] / freqs(r, 0.9)[1] - 1.5, isco(0.9) + 1e-6, 50.0)
Op, Or, Oz = freqs(r32, 0.9)
chk("3:2-Resonanz Ω_θ:Ω_r = 3:2 bei χ=0.9", Oz / Or, 1.5, 1e-12)
print(f"   Radius der 3:2-Resonanz (χ=0.9): r̃ = {r32:.5f}  (Ω_r/Ω_φ = {Or/Op:.5f}, Ω_θ/Ω_φ = {Oz/Op:.5f})")
for name in ("Sgr A*", "M87*"):
    M = BODY[name][0]; u = c ** 3 / (G * M) / (2 * math.pi)
    print(f"   {name}: bei 3:2: f_θ = {Oz*u:.4e} Hz, f_r = {Or*u:.4e} Hz;  f_Knoten = {(Op-Oz)*u:.4e} Hz;"
          f"  f_φ = {Op*u:.4e} Hz")

print("\n   Knotenzahl an der Schwelle r* (schwaches Feld): ζ* = 2χ(3Δs/r*)^{3/2} = 2χ(v_flach/c)³")
for name in ("Erde", "Sonne", "Sgr A*", "M87*", "Milchstraße"):
    M, J, R, ds, a, chi = BODY[name]
    rs = math.sqrt(G * M / a0); v = (G * M * a0) ** 0.25
    z1 = 2 * chi * (3 * ds / rs) ** 1.5
    z2 = 2 * chi * (v / c) ** 3
    print(f"   {name:12s} r* = {rs:.3e} m   v_flach = {v/1e3:8.2f} km/s   χ = {chi:9.4g}   ζ* = {z1:.3e}")
    chk(f"{name}: 2χ(3Δs/r*)^{{3/2}} = 2χ(v_flach/c)³", z1, z2, 1e-12)

# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("3) SIEGEL-PERIODENMATRIX Ω ∈ ℍ₂ UND IGUSA-SPITZENFORM χ10 = −2^{-14} ∏_{10 gerade m} θ_m(Ω)²")
print("=" * 100)

def eta(tau):
    q = cmath.exp(2j * math.pi * tau); p = 1 + 0j; qn = q
    for n in range(1, 20000):
        p *= (1 - qn); qn *= q
        if abs(qn) < 1e-19: break
    return cmath.exp(2j * math.pi * tau / 24) * p

def Delta(tau): return eta(tau) ** 24

def E2(tau):
    q = cmath.exp(2j * math.pi * tau); s = 0j; qn = q
    for n in range(1, 20000):
        s += n * qn / (1 - qn); qn *= q
        if abs(qn) * n < 1e-19: break
    return 1 - 24 * s

def th1g(a, b, w, tau, N=60):
    n = np.arange(-N, N + 1); v = n + a / 2
    return np.exp(1j * math.pi * tau * v * v + 2j * math.pi * v * (w + b / 2)).sum()

def jac_theta1(tau, z):          # Jacobi θ1, θ1'(0) = 2π η³
    return -th1g(1, 1, z, tau)

EVEN = [(a1, a2, b1, b2) for a1 in (0, 1) for a2 in (0, 1) for b1 in (0, 1) for b2 in (0, 1)
        if (a1 * b1 + a2 * b2) % 2 == 0]
_GRID = {}
def theta2(Om, ch):
    lam = np.linalg.eigvalsh(Om.imag).min()
    N = int(math.sqrt(46 / (math.pi * lam))) + 3
    if N not in _GRID:
        n = np.arange(-N, N + 1); _GRID[N] = np.meshgrid(n, n, indexing="ij")
    n1, n2 = _GRID[N]
    v1 = n1 + ch[0] / 2; v2 = n2 + ch[1] / 2
    Q = Om[0, 0] * v1 * v1 + 2 * Om[0, 1] * v1 * v2 + Om[1, 1] * v2 * v2
    return np.exp(1j * math.pi * (Q + v1 * ch[2] + v2 * ch[3])).sum()

def chi10(Om):
    p = 1 + 0j
    for ch in EVEN:
        p *= theta2(Om, ch) ** 2
    return -p / 2 ** 14

def Omat(t1, z, t2): return np.array([[t1, z], [z, t2]], dtype=complex)

print(f"   Anzahl gerader Charakteristiken: {len(EVEN)}")
print("   Geschlecht 1 (Grundlage):")
for tau in (0.1 + 1.1j, -0.3 + 0.7j):
    t2, t3, t4 = th1g(1, 0, 0, tau), th1g(0, 0, 0, tau), th1g(0, 1, 0, tau)
    chk(f"τ={tau}: θ2θ3θ4 = 2η³", t2 * t3 * t4, 2 * eta(tau) ** 3, 1e-12)
    h = 1e-6
    d = (jac_theta1(tau, h) - jac_theta1(tau, -h)) / (2 * h)
    chk(f"τ={tau}: θ1'(0) = 2πη³", d, 2 * math.pi * eta(tau) ** 3, 1e-8)
    chk(f"τ={tau}: Δ(−1/τ) = τ¹²Δ(τ)", Delta(-1 / tau), tau ** 12 * Delta(tau), 1e-10)
    h = 1e-5
    dl = (cmath.log(Delta(tau + h)) - cmath.log(Delta(tau - h))) / (2 * h)
    chk(f"τ={tau}: d log Δ/dτ = 2πi E2", dl, 2j * math.pi * E2(tau), 1e-8)

print("\n   Modularität von χ10 (Gewicht 10) unter Sp(4,ℤ):")
Om0 = Omat(0.3 + 1.1j, 0.2 + 0.3j, -0.1 + 0.9j)
x0 = chi10(Om0)
print(f"   χ10(Ω0) = {x0:.10e}")
chk("J: χ10(−Ω⁻¹) = det(Ω)^10 χ10(Ω)", chi10(-np.linalg.inv(Om0)), np.linalg.det(Om0) ** 10 * x0, 1e-9)
chk("Translation Ω + [[1,0],[0,0]]", chi10(Om0 + np.array([[1, 0], [0, 0]])), x0, 1e-10)
chk("Translation Ω + [[0,1],[1,0]] (Nebendiagonale + 1)", chi10(Om0 + np.array([[0, 1], [1, 0]])), x0, 1e-10)
U = np.array([[1, 0], [1, 1]])
chk("GL(2,ℤ): χ10(UΩUᵀ) = χ10(Ω), U = [[1,0],[1,1]]", chi10(U @ Om0 @ U.T), x0, 1e-10)
U = np.diag([1, -1])
chk("Spiegelung U = diag(1,−1): Nebendiagonale z → −z, χ10 unverändert", chi10(U @ Om0 @ U.T), x0, 1e-10)

print("\n   Diagonaler Rand Ω = diag(τ1, τ2):")
t1, t2 = 0.1 + 0.8j, -0.2 + 1.0j
Od = Omat(t1, 0, t2)
p9 = 1 + 0j; vanish = None
for ch in EVEN:
    th = theta2(Od, ch)
    if ch == (1, 1, 1, 1): vanish = th
    else: p9 *= th ** 2
print(f"   θ[11;11](diag) = {abs(vanish):.2e}  (verschwindet exakt: ungerade × ungerade)")
chk("Produkt der 9 übrigen θ² = 2¹² η(τ1)¹⁸ η(τ2)¹⁸", p9, 2 ** 12 * eta(t1) ** 18 * eta(t2) ** 18, 1e-10)

print("\n   Nähe des Randes: Ω = diag(τ1,τ2) + ε·[[0,1],[1,0]] ;  Q(ε) := χ10/(π² Δ(τ1) Δ(τ2))")
D1, D2 = Delta(t1), Delta(t2)
cs = {}
for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4):
    Q = chi10(Omat(t1, eps, t2)) / (math.pi ** 2 * D1 * D2)
    cs[eps] = (1 - Q / eps ** 2) / eps ** 2
    print(f"   ε = {eps:7.0e}   Q/ε² = {Q/eps**2:.12f}   (1 − Q/ε²)/ε² = {cs[eps]:.8f}")
chk("führende Ordnung: χ10 = π² ε² Δ(τ1)Δ(τ2) (ε = 1e-4)",
    chi10(Omat(t1, 1e-4, t2)) / (math.pi ** 2 * 1e-8 * D1 * D2), 1.0, 1e-7)
e1, e2 = 1e-2, 2e-2
c1 = (1 - chi10(Omat(t1, e1, t2)) / (math.pi ** 2 * D1 * D2 * e1 ** 2)) / e1 ** 2
c2 = (1 - chi10(Omat(t1, e2, t2)) / (math.pi ** 2 * D1 * D2 * e2 ** 2)) / e2 ** 2
c_rich = (4 * c1 - c2) / 3
chk("nächste Ordnung: χ10 = π²ε²Δ1Δ2 [1 − (π²/3) E2(τ1) E2(τ2) ε² + O(ε⁴)]",
    c_rich, math.pi ** 2 / 3 * E2(t1) * E2(t2), 1e-5)
chk("gerade in ε: χ10(+ε) = χ10(−ε)", chi10(Omat(t1, -0.13, t2)), chi10(Omat(t1, 0.13, t2)), 1e-10)
chk("periodisch in ε: χ10(ε+1) = χ10(ε)", chi10(Omat(t1, 1.13, t2)), chi10(Omat(t1, 0.13, t2)), 1e-9)
print(f"   ε = 1 (ganzzahlig): |χ10| = {abs(chi10(Omat(t1, 1.0, t2))):.2e}   (Rand, zerfällt wieder)")

print("\n   Spitze (q → 0): χ10 = q1 q2 sin²(πε) + O(q³);  hier Im τ = 3:")
t1s, t2s = 0.2 + 3j, -0.1 + 3j
q1, q2 = cmath.exp(2j * math.pi * t1s), cmath.exp(2j * math.pi * t2s)
for eps in (0.1, 0.3, 0.5):
    chk(f"ε={eps}: χ10/(q1 q2) = sin²(πε)", chi10(Omat(t1s, eps, t2s)) / (q1 * q2), math.sin(math.pi * eps) ** 2, 1e-6)

print("\n   Fourier-Jacobi-Rand (ein Zyklus fällt weg, q2 → 0): χ10 → ¼ q2 η(τ1)¹⁸ θ1(τ1,z)²")
t1f = 0.15 + 0.9j
for z in (0.05, 0.2 + 0.1j):
    t2f = 0.1 + 4.0j; q2 = cmath.exp(2j * math.pi * t2f)
    chk(f"z={z}: χ10/(q2 η¹⁸ θ1²)", chi10(Omat(t1f, z, t2f)) / (q2 * eta(t1f) ** 18 * jac_theta1(t1f, z) ** 2), 0.25, 1e-8)

print("\n   Genauigkeit im Fernfeld (kleines Im Ω): Q(ε)/ε² bei τ1 = 1.41+iη, τ2 = 1+iη, Erwartung 1 − (π²/3)E2E2ε²")
for et in (0.3, 0.1, 0.05):
    ta, tb = 1.41 + 1j * et, 1.0 + 1j * et
    for eps in (1e-3, 1e-4):
        Q = chi10(Omat(ta, eps, tb)) / (math.pi ** 2 * Delta(ta) * Delta(tb))
        pred = 1 - math.pi ** 2 / 3 * E2(ta) * E2(tb) * eps ** 2
        print(f"   η = {et:4.2f}, ε = {eps:.0e}: Q/ε² = {Q/eps**2:.8f}   Vorhersage {pred:.8f}"
              f"   (π²/3)|E2E2| = {math.pi**2/3*abs(E2(ta)*E2(tb)):.3e}")
    chk(f"η={et}, ε=1e-4: Fernfeld-Numerik", Q / eps ** 2, pred, 1e-6)
print("   -> Die Entwicklung gilt nur für ε² ≪ 3/(π²|E2E2|); im Fernfeld wird |E2| ~ 1/η² groß.")

print("\n   Δ(τ) im Fernfeld mit Realteil: τ = iη (Notiz 01) gegen τ = √2 + iη (flache Kurve, κ/Ω = √2):")
print(f"   {'η':>6s} {'|Δ(iη)|':>11s} {'|Δ(√2+iη)|':>12s} {'η⁶|Δ(√2+iη)|':>14s} {'η⁶|Δ(1.3+iη)|':>14s}")
for et in (1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03):
    d0, d1, d2 = abs(Delta(1j * et)), abs(Delta(math.sqrt(2) + 1j * et)), abs(Delta(1.3 + 1j * et))
    print(f"   {et:6.2f} {d0:11.3e} {d1:12.3e} {et**6*d1:14.3e} {et**6*d2:14.3e}")
tt = 0.2 + 0.3j
chk("Invariante Im(τ)⁶|Δ(τ)| unter τ → −1/τ", 0.3 ** 6 * abs(Delta(tt)),
    ((-1 / tt).imag) ** 6 * abs(Delta(-1 / tt)), 1e-10)

print("\n   Kerr-Periodenmatrix (Setzung): τ_r = Ω_r/Ω_φ + iη, τ_θ = Ω_θ/Ω_φ + iη, Ω12 = ζ = 1 − Ω_θ/Ω_φ")
print("   Für Schwarze Löcher ist η = r*/(2πr) astronomisch groß (q ≈ 0) ⇒ Q = sin²(πζ)/π² exakt bis e^{−2πη}:")
for name in ("Sgr A*",):
    M = BODY[name][0]
    rs = math.sqrt(G * M / a0); rg = G * M / c ** 2
    for r in (isco(0.9), 6.0, 20.0, 1e3):
        Op, Or, Oz = freqs(r, 0.9)
        et = rs / (2 * math.pi * r * rg)
        zt = 1 - Oz / Op
        print(f"   {name} r̃ = {r:8.3f}: η = {et:.3e}, Re τ_r = {Or/Op:.5f}, Re τ_θ = {Oz/Op:.5f}, ζ = {zt:.3e},"
              f" Q = sin²(πζ)/π² = {math.sin(math.pi*zt)**2/math.pi**2:.3e}")
print("   Numerischer Beleg bei moderatem η = 2 (sonst Unterlauf):")
for r in (isco(0.9), 6.0):
    Op, Or, Oz = freqs(r, 0.9); zt = 1 - Oz / Op
    ta, tb = Or / Op + 2j, Oz / Op + 2j
    Q = chi10(Omat(ta, zt, tb)) / (math.pi ** 2 * Delta(ta) * Delta(tb))
    print(f"   r̃ = {r:6.3f}: ζ = {zt:.5f}, Q = {Q.real:.6e}{Q.imag:+.1e}i, sin²(πζ)/π² = {math.sin(math.pi*zt)**2/math.pi**2:.6e}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("4) FOLGERUNGEN FÜR SCHEIBEN — exponentielle Scheibe, M = 6e10 M☉, R_d = 2.6 kpc, Höhe h")
print("=" * 100)

def I_n(n, x, m=4000):
    t = np.linspace(0, math.pi, m + 1)
    f = np.exp(x * np.cos(t)) * np.cos(n * t)
    return (f.sum() - 0.5 * (f[0] + f[-1])) * (math.pi / m) / math.pi

def K_n(n, x, m=8000):
    T = math.acosh(80 / x + 1)
    t = np.linspace(0, T, m + 1)
    f = np.exp(-x * np.cosh(t)) * np.cosh(n * t)
    return (f.sum() - 0.5 * (f[0] + f[-1])) * (T / m)

chk("I0(1)", I_n(0, 1.0), 1.2660658777520082, 1e-12)
chk("I1(1)", I_n(1, 1.0), 0.5651591039924851, 1e-12)
chk("K0(1)", K_n(0, 1.0), 0.42102443824070834, 1e-10)
chk("K1(1)", K_n(1, 1.0), 0.6019072301972346, 1e-10)
chk("I0(2.5) gegen numpy.i0", I_n(0, 2.5), float(np.i0(2.5)), 1e-12)

Md, Rd = 6e10 * MS, 2.6 * KPC
Sig0 = Md / (2 * math.pi * Rd * Rd)
def vN2(R):
    y = R / (2 * Rd)
    return 4 * math.pi * G * Sig0 * Rd * y * y * (I_n(0, y) * K_n(0, y) - I_n(1, y) * K_n(1, y))
def nu_f(y): return 1 / (1 - math.exp(-math.sqrt(y)))

def disk_freqs(R, h):
    dR = 1e-3 * R
    v2 = vN2(R); dv2 = (vN2(R + dR) - vN2(R - dR)) / (2 * dR)
    rho0 = Sig0 * math.exp(-R / Rd) / (2 * h)
    OmN2 = v2 / R ** 2
    kapN2 = dv2 / R + 2 * v2 / R ** 2          # κ² = (1/R³) d(R²v²)/dR ... = dv²/dR/R + 2v²/R²
    nuzN2 = 4 * math.pi * G * rho0 - dv2 / R
    gN = v2 / R; y = gN / a0; nu = nu_f(y)
    dnu = (nu_f(vN2(R + dR) / (R + dR) / a0) - nu_f(vN2(R - dR) / (R - dR) / a0)) / (2 * dR)
    Om2 = nu * OmN2
    kap2 = (nu * dv2 + dnu * v2) / R + 2 * nu * v2 / R ** 2
    nuzK2 = nuzN2 + (nu - 1) * OmN2
    nuzF2 = nu * nuzN2
    return dict(v2=v2, dv2=dv2, rho0=rho0, OmN2=OmN2, kapN2=kapN2, nuzN2=nuzN2, y=y, nu=nu, dnu=dnu,
                Om2=Om2, kap2=kap2, nuzK2=nuzK2, nuzF2=nuzF2)

print("\n   Feldform (QUMOND-artig): ∇²Φ = ∇·[ν ∇Φ_N]. In der Mittelebene mit g_R = ν g_N,R folgt exakt ν_z² = ν·ν_z,N²:")
for R in (6 * KPC, 12 * KPC, 20 * KPC):
    d = disk_freqs(R, 0.3 * KPC)
    lap = d["nu"] * 4 * math.pi * G * d["rho0"] + d["dnu"] * d["v2"] / R
    dR = 1e-3 * R
    def RgR(x):
        return nu_f(vN2(x) / x / a0) * vN2(x)
    radial = (RgR(R + dR) - RgR(R - dR)) / (2 * dR) / R
    chk(f"R = {R/KPC:.0f} kpc: ∇²Φ − (1/R)∂_R(R g_R) = ν·ν_z,N²", lap - radial, d["nu"] * d["nuzN2"], 1e-5)

print("\n   Kernform (radialer Projektor e_r⊗e_r, Kugel-Radialrichtung): ν_z² = ν_z,N² + (ν−1)Ω_N²")
print("   Beide Formen haben dieselbe Rotationskurve (Ω², κ²) und unterscheiden sich NUR vertikal.")
for h in (0.15, 0.3, 0.6):
    print(f"\n   h = {h} kpc   (Ω, ν_z, Ω_p = Ω − ν_z in km/s/kpc; ν = Rückkopplungsfaktor 1/(1−q))")
    print(f"   {'R':>5s} {'y':>7s} {'ν':>6s} {'v':>6s} {'Ω':>7s} {'κ/Ω':>6s} {'ν_z/Ω N':>8s} {'ν_z/Ω K':>8s} {'ν_z/Ω F':>8s}"
          f" {'Ω_p N':>8s} {'Ω_p K':>8s} {'Ω_p F':>8s}")
    for Rk in (4, 8.2, 12, 16, 20, 25, 30):
        d = disk_freqs(Rk * KPC, h * KPC)
        Om = math.sqrt(d["Om2"]); OmN = math.sqrt(d["OmN2"])
        rN, rK, rF = math.sqrt(d["nuzN2"]) / OmN, math.sqrt(d["nuzK2"]) / Om, math.sqrt(d["nuzF2"]) / Om
        print(f"   {Rk:5.1f} {d['y']:7.3f} {d['nu']:6.3f} {math.sqrt(d['v2']*d['nu'])/1e3:6.1f} {Om/KMSKPC:7.2f}"
              f" {math.sqrt(d['kap2'])/Om:6.3f} {rN:8.3f} {rK:8.3f} {rF:8.3f}"
              f" {(OmN-rN*OmN)/KMSKPC:8.2f} {(Om-rK*Om)/KMSKPC:8.2f} {(Om-rF*Om)/KMSKPC:8.2f}")
    d = disk_freqs(20 * KPC, h * KPC)
    chk(f"h={h}: F erhält ν_z/Ω aus Newton", math.sqrt(d["nuzF2"] / d["Om2"]), math.sqrt(d["nuzN2"] / d["OmN2"]), 1e-12)
    chk(f"h={h}: K erhält ν_z² − Ω² aus Newton", d["nuzK2"] - d["Om2"], d["nuzN2"] - d["OmN2"], 1e-9)

print("\n   Aufwickelzeit einer Verbiegung (Warp) zwischen 12 und 20 kpc: t = 2π/|Ω_p(12) − Ω_p(20)|  (h = 0.3 kpc)")
res = {}
for lab in ("N", "K", "F"):
    ops = []
    for Rk in (12, 20):
        d = disk_freqs(Rk * KPC, 0.3 * KPC)
        if lab == "N": Om, nz = math.sqrt(d["OmN2"]), math.sqrt(d["nuzN2"])
        elif lab == "K": Om, nz = math.sqrt(d["Om2"]), math.sqrt(d["nuzK2"])
        else: Om, nz = math.sqrt(d["Om2"]), math.sqrt(d["nuzF2"])
        ops.append(Om - nz)
    t = 2 * math.pi / abs(ops[0] - ops[1]) / (1e9 * YR)
    res[lab] = t
    print(f"   {lab}: Ω_p(12) = {ops[0]/KMSKPC:7.3f}, Ω_p(20) = {ops[1]/KMSKPC:7.3f} km/s/kpc  ->  t_wind = {t:6.2f} Gyr")
print(f"   Verhältnis t_wind(K)/t_wind(F) = {res['K']/res['F']:.2f}")

print("\n   Vertikale Dynamikmasse am Sonnenort (|z| < 1.1 kpc), Σ_bar = 47.1 M☉/pc² (McKee+2015):")
Sbar = 47.1 * MS / PC ** 2
gz_bar = 2 * math.pi * G * Sbar
zc, R0 = 1.1 * KPC, 8.2 * KPC
for lab, gNR in (("Scheibenmodell", vN2(R0) / R0), ("Punktmasse (07)", 1.244e-10)):
    nuR = nu_f(gNR / a0)
    gK = gz_bar + (nuR - 1) * gNR * zc / R0
    gtot = math.hypot(gNR, gz_bar)
    gF = nu_f(gtot / a0) * gz_bar
    print(f"   {lab:16s}: g_N,R = {gNR:.3e} m/s²  ->  Σ_dyn(K) = {gK/(2*math.pi*G)/(MS/PC**2):5.1f},"
          f"  Σ_dyn(F) = {gF/(2*math.pi*G)/(MS/PC**2):5.1f} M☉/pc²   (gemessen 71 ± 6, Kuijken & Gilmore 1991)")

print("\n   Siegel-Größen der Scheibe (h = 0.3 kpc): τ_r = κ/Ω + iη, τ_θ = ν_z/Ω + iη, ζ = 1 − ν_z/Ω, η = r*/(2πR)")
rstar = math.sqrt(G * Md / a0)
print(f"   r* = {rstar/KPC:.3f} kpc")
print(f"   {'R':>5s} {'η':>6s} {'ζ_K':>8s} {'Q_K':>11s} {'sin²/π² K':>11s} {'ζ_F':>8s} {'Q_F':>11s} {'sin²/π² F':>11s}")
for Rk in (4, 8.2, 12, 16, 20):
    d = disk_freqs(Rk * KPC, 0.3 * KPC)
    et = rstar / (2 * math.pi * Rk * KPC)
    Om = math.sqrt(d["Om2"]); kr = math.sqrt(d["kap2"]) / Om
    out = []
    for key in ("nuzK2", "nuzF2"):
        rz = math.sqrt(d[key]) / Om
        zt = 1 - rz
        ta, tb = kr + 1j * et, rz + 1j * et
        Q = chi10(Omat(ta, zt, tb)) / (math.pi ** 2 * Delta(ta) * Delta(tb))
        out.append((zt, Q, math.sin(math.pi * zt) ** 2 / math.pi ** 2))
    (zK, QK, sK), (zF, QF, sF) = out
    print(f"   {Rk:5.1f} {et:6.3f} {zK:8.4f} {abs(QK):11.3e} {sK:11.3e} {zF:8.4f} {abs(QF):11.3e} {sF:11.3e}")

print("\n" + "=" * 100)
print(f"PRÜFUNGEN: {NCHK[0]} durchgeführt, {NCHK[1]} fehlgeschlagen")
print("=" * 100)
