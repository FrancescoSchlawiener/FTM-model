# -*- coding: utf-8 -*-
"""
FTM / Zustandsraum - explizite Feldgleichung (Punkte 1-4 des Autors)
Rechnungen zu 07-feldgleichung.md. Nur Python-Standardbibliothek.
Aufruf:  python feldgleichung.py > feldgleichung.txt   (Variante B = Hauptfassung)
"""
import math, cmath, sys, io

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ---------------------------------------------------------------- Konstanten
# Einheitliche Konstanten des Weltmodells (01-08)
G     = 6.674e-11            # m^3 kg^-1 s^-2
c     = 299792458.0          # m/s
Msun  = 1.989e30             # kg
GMsun = G * Msun             # m^3/s^2
AU    = 1.495978707e11       # m
pc    = 3.0856775814913673e16
kpc   = 1e3 * pc
a0    = 1.2e-10              # m/s^2
v0_MW = 220e3                # m/s  (v0 aus dem Toy Paper des Autors)
M_MW  = 6e10 * Msun          # baryonische Milchstrasse
GM_MW = G * M_MW

def H(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ================================================================ PUNKT 2
H("PUNKT 2: Koeffizienten tau_R(n) von Delta = q * prod_{n>=1} (1-q^n)^24")
NMAX = 30
# Polynom-Multiplikation abgeschnitten bei q^NMAX
poly = [0] * (NMAX + 1); poly[0] = 1           # prod (1-q^n)^24
for n in range(1, NMAX + 1):
    for _ in range(24):
        new = poly[:]
        for i in range(NMAX + 1 - n):
            new[i + n] -= poly[i]
        poly = new
tauR = {m: poly[m - 1] for m in range(1, NMAX + 1)}   # Delta = sum tauR(m) q^m
for m in range(1, NMAX + 1):
    print(f"  tau_R({m:2d}) = {tauR[m]:>22d}")

# Kontrolle: Produkt und Reihe an Probepunkten
print("\n  Kontrolle Produkt vs. Reihe:")
def delta_prod(q, N=400):
    p = q
    for n in range(1, N):
        p *= (1 - q ** n) ** 24
    return p
def delta_series(q):
    return sum(tauR[m] * q ** m for m in range(1, NMAX + 1))
for q in [0.01, 0.05, 0.1, 0.02 + 0.03j]:
    print(f"   q={q}:  Produkt={delta_prod(q):.15g}  Reihe(n<=30)={delta_series(q):.15g}")

# Faktor 24: q^2-Koeffizient aus den 24 Faktoren (1-q)
print("\n  q^2-Koeffizient: aus (1-q)^24 -> -C(24,1) = ", -math.comb(24, 1))
print("  q^3-Koeffizient: C(24,2) - C(24,1) [aus (1-q)^24 und (1-q^2)^24] =",
      math.comb(24, 2) - math.comb(24, 1))

# Gewicht 12: Delta(-1/tau) = tau^12 Delta(tau)
print("\n  Gewicht 12, Kontrolle Delta(-1/tau) / (tau^12 Delta(tau)):")
def delta_tau(tau):
    return delta_prod(cmath.exp(2j * math.pi * tau), N=300)
for tau in [0.3 + 1.1j, -0.2 + 0.9j, 0.1 + 1.5j]:
    r = delta_tau(-1 / tau) / (tau ** 12 * delta_tau(tau))
    print(f"   tau={tau}:  Verhaeltnis = {r.real:.12f} {r.imag:+.2e} i")
    r2 = delta_tau(tau + 1) / delta_tau(tau)
    print(f"               Delta(tau+1)/Delta(tau) = {r2.real:.12f} {r2.imag:+.2e} i")

# ================================================================ PUNKT 1
H("PUNKT 1: Kreisbahn-Auswertung des Gedaechtnisintegrals")
print("""  r(t') = r (cos W t', sin W t'),  rdot(t') = r W e_phi(t'),  s = t - t'
  (e_r(t) (x) e_r(t)) rdot(t') = e_r(t) * r W * sin(W s)
  I = int_0^inf kappa0 M e^{-lam s} r W sin(W s) ds = kappa0 M r W^2/(lam^2+W^2)""")
def I_num(kM, lam, r, W, N=200000):
    smax = 60.0 / lam
    h = smax / N
    tot = 0.0
    for i in range(N + 1):
        s = i * h
        w = 1 if i in (0, N) else (4 if i % 2 else 2)
        tot += w * kM * math.exp(-lam * s) * r * W * math.sin(W * s)
    return tot * h / 3
for kM, lam, r, W in [(1.0, 1.0, 1.0, 1.0), (2.0, 0.3, 1.5, 1.2), (0.5, 3.0, 2.0, 0.7)]:
    num = I_num(kM, lam, r, W)
    ana = kM * r * W**2 / (lam**2 + W**2)
    print(f"   kM={kM}, lam={lam}, r={r}, W={W}:  numerisch={num:.10f}  analytisch={ana:.10f}")

print("\n  Exakte Kreisbahnloesung u = v^2 aus")
print("  u^2 + u (lam^2 r^2 - kappa0 M r^2 - GM/r) - GM lam^2 r = 0")
def u_exact(GM, r, kM, lam):
    B = lam**2 * r**2 - kM * r**2 - GM / r
    return 0.5 * (-B + math.sqrt(B * B + 4 * GM * lam**2 * r))
# Kontrolle: Einsetzen in v^2/r = GM/r^2 + kM (v^2/r)/(lam^2+v^2/r^2)
for GM, r, kM, lam in [(1.0, 1.0, 0.3, 0.8), (1.0, 5.0, 0.01, 0.05), (1.0, 2.0, 0.0, 1.0)]:
    u = u_exact(GM, r, kM, lam)
    lhs = u / r
    rhs = GM / r**2 + kM * (u / r) / (lam**2 + u / r**2)
    print(f"   GM={GM}, r={r}, kM={kM}, lam={lam}: u={u:.10f}  lhs={lhs:.10f} rhs={rhs:.10f}")

# ================================================================ PUNKT 3/4 Zahlen
def ds(GM): return GM / (3 * c**2)

H("PUNKT 3: Faltungstiefe, q(r), eta(r)")
for name, GM in [("Sonne", GMsun), ("Milchstrasse (6e10 Msun)", GM_MW)]:
    d = ds(GM)
    print(f"  {name}: GM = {GM:.6e} m^3/s^2   Delta_s = GM/(3c^2) = {d:.6e} m"
          f" = {d/AU:.6e} AU = {d/pc:.6e} pc")

# Referenzlaenge L in c^2/L.  Variante B (Hauptfassung): L = r*/3, r* = sqrt(GM/a0)
#                            Variante A (Vergleich):   L = GM/(3 v0^2), v0 = 220 km/s
# (Das Symbol r0 bleibt im Modell dem Grenzradius 6 Ds vorbehalten.)
rstar_MW = math.sqrt(GM_MW / a0)
L_B = rstar_MW / 3                            # Variante B, Milchstrasse
L_A = GM_MW / (3 * v0_MW**2)                  # Variante A, Milchstrasse
rstar_S = math.sqrt(GMsun / a0)
L_B_S = rstar_S / 3                           # Variante B, Sonne
print("\n  Referenzlaenge in c^2/L (Milchstrasse):")
print(f"   B: a_eff(r*) = a0  -> L = r*/3 = sqrt(GM/a0)/3 = {L_B:.6e} m = {L_B/kpc:.6f} kpc  (r* = {rstar_MW/kpc:.6f} kpc)")
print(f"      a_eff = sqrt(GM a0)/r,  v_rot = (GM a0)^(1/4) = {(GM_MW*a0)**0.25/1e3:.6f} km/s,  v_rot^4 = G M a0")
print(f"   A: v_rot = v0 = 220 km/s -> L = GM/(3 v0^2) = {L_A:.6e} m = {L_A/kpc:.6f} kpc")
print(f"  Variante B (Sonne): r*/3 = {L_B_S:.6e} m = {L_B_S/AU:.4f} AU (r* = {rstar_S/AU:.4f} AU),"
      f"  v_rot = (GM a0)^(1/4) = {(GMsun*a0)**0.25:.6f} m/s")
print(f"  c^2/L = 3 c^2/r*:  B(MW): {c**2/L_B:.6e} m/s^2   B(Sonne): {c**2/L_B_S:.6e} m/s^2   A(MW): {c**2/L_A:.6e} m/s^2")
print(f"  Kontrolle B (MW): (Ds/r)(3c^2/r*) bei r = r* : {ds(GM_MW)/rstar_MW*3*c**2/rstar_MW:.6e} m/s^2 (= a0)")

def row(GM, r, L):
    d = ds(GM); q = d / r
    eta = math.log(r / d) / (2 * math.pi)
    aN = GM / r**2
    aeff = GM / (3 * r * L)
    V = math.sqrt(GM / (3 * L))
    vN = math.sqrt(GM / r)
    vtot = math.sqrt(GM / r + V**2)
    dq = delta_series(q)
    rel = dq / q - 1
    ratio = math.sqrt((GM / r + 2 * V**2) / (GM / r + V**2))   # omega_rad/omega_azim
    return dict(q=q, eta=eta, aN=aN, aeff=aeff, V=V, vN=vN, vtot=vtot, rel=rel,
                ratio=ratio, reg=("a=a0 Schwelle" if abs(aN / a0 - 1) < 1e-6 else
                     ("a>a0 gekoppelt" if aN > a0 else "a<a0 kritisch")))

hdr = (f"  {'r[kpc]':>7} {'q=Ds/r':>11} {'eta':>7} {'Delta/q-1':>11} {'a_N':>11} {'a_eff':>11}"
       f" {'v_N':>8} {'v_rot':>8} {'v_tot':>8} {'w_r/w_phi':>9}  Regime")
def table(GM, L, radii, unit, uname, vunit=1e3):
    print(hdr.replace("r[kpc]", f"r[{uname}]"))
    for R in radii:
        d = row(GM, R * unit, L)
        print(f"  {R:7.4g} {d['q']:11.4e} {d['eta']:7.4f} {d['rel']:11.3e} {d['aN']:11.4e} {d['aeff']:11.4e}"
              f" {d['vN']/vunit:8.4g} {d['V']/vunit:8.4g} {d['vtot']/vunit:8.4g} {d['ratio']:9.5f}  {d['reg']}")
radii = [5, rstar_MW / kpc, 10, 15, 20, 30, 40, 50]
H("TABELLE Milchstrasse, Variante B (Hauptfassung, a_eff = sqrt(GM a0)/r)")
table(GM_MW, L_B, radii, kpc, "kpc")
H("TABELLE Milchstrasse, Variante A (Vergleich, v_rot = 220 km/s)")
table(GM_MW, L_A, radii, kpc, "kpc")
H("TABELLE Sonne, Variante B; Geschwindigkeiten in km/s")
table(GMsun, L_B_S, [1, 30.07, 100, 1000, rstar_S / AU, 1e4, 1e5], AU, "AU")
print("\n  Sonne/Erdbahn: a_eff/a_N = r/r* =", f"{AU/rstar_S:.4e}")

H("PUNKT 1/3: Kernnorm, Gedaechtnisrate, Kreisbahnbedingung")
print("  Zeitintegral des Kerns: int_0^inf ||K|| dDt = kappa0 M_k / lam   [1/s]")
print("  Kreisbahn: a_K = kappa0 M_k v W/(lam^2+W^2) = kappa0 M_k v / lam_eff, lam_eff = W + lam^2/W")
print("  lam_eff minimal bei lam = W: lam_eff,min = 2W")
print("  Aus v_rot = kappa0 Ds/lam_eff:  kappa0/lam_eff = v_rot/Ds = 3c^2/(r* v_rot)  (Variante B)")
VB = (GM_MW*a0)**0.25
print(f"   Milchstrasse B: kappa0/lam_eff = {VB/ds(GM_MW):.6e} 1/s   (Kontrolle 3c^2/(r* v_rot) = {3*c**2/(rstar_MW*VB):.6e})")
VBS = (GMsun*a0)**0.25
print(f"   Sonne B:        kappa0/lam_eff = {VBS/ds(GMsun):.6e} 1/s")
print(f"   Milchstrasse A: kappa0/lam_eff = {v0_MW/ds(GM_MW):.6e} 1/s")
for R in [5, 10, 20, 50]:
    r = R * kpc
    d = row(GM_MW, r, L_B)
    W = d['vtot'] / r
    print(f"   Milchstrasse B, r={R} kpc: W = v_tot/r = {W:.4e} 1/s  (Umlaufzeit {2*math.pi/W/3.15576e7/1e6:.1f} Myr)")
W_E = math.sqrt(GMsun / AU**3)
d_E = ds(GMsun) / AU
print(f"   Erdbahn: W = {W_E:.6e} 1/s (Umlaufzeit {2*math.pi/W_E/86400:.3f} d);  q = {d_E:.4e};"
      f"  eta = {math.log(1/d_E)/(2*math.pi):.4f};  Ds_Sonne = {ds(GMsun):.4f} m")

print("\n  Grenzfaelle der exakten Kreisbahnloesung mit M_k = q = Ds/r:")
print("   lam -> 0     :  v^2 = GM/r + kappa0 Ds r = GM/r + kappa0 GM r/(3c^2)")
print("   lam^2 r^2 >> :  v^2 = (GM/r) / (1 - kappa0 q/lam^2)")
print("  Bedingung fuer das Profil von Punkt 3 (a_K = V^2/r) auf der Kreisbahn:")
print("   kappa0 Ds W^2/(lam^2+W^2) = V^2/r   <=>   lam^2 = W^2 (kappa0 Ds r / V^2 - 1)")

H("PUNKT 4: Schwellenradius r* = sqrt(GM/a0)")
for name, GM in [("Sonne", GMsun), ("Milchstrasse", GM_MW)]:
    rs = math.sqrt(GM / a0)
    qs = math.sqrt(GM * a0) / (3 * c**2)
    etas = math.log(1 / qs) / (2 * math.pi)
    print(f"  {name}: r* = {rs:.6e} m = {rs/AU:.6e} AU = {rs/pc:.6e} pc = {rs/kpc:.6e} kpc")
    print(f"          q* = sqrt(GM a0)/(3c^2) = {qs:.6e}   (Kontrolle Ds/r* = {ds(GM)/rs:.6e})"
          f"   eta* = {etas:.6f}")
    print(f"          (GM a0)^(1/4) = {(GM*a0)**0.25:.6e} m/s ;  3 c^2 q* = {3*c**2*qs:.6e} m^2/s^2 = sqrt(GM a0)")
print("\n  Lage relativ zur Schwelle (a_N/a0, r/r*):")
for name, GM, r in [("Erdbahn 1 AU", GMsun, AU), ("Neptunbahn 30.07 AU", GMsun, 30.07 * AU),
                    ("MW 5 kpc", GM_MW, 5 * kpc), ("MW 10 kpc", GM_MW, 10 * kpc),
                    ("MW 20 kpc", GM_MW, 20 * kpc), ("MW 50 kpc", GM_MW, 50 * kpc)]:
    aN = GM / r**2; rs = math.sqrt(GM / a0)
    print(f"   {name:22s}: a_N = {aN:.4e} m/s^2  a_N/a0 = {aN/a0:.4e}   r/r* = {r/rs:.4e}")
print(f"   Sonne im Feld der Milchstrasse (R = 8.2 kpc): a_N = {GM_MW/(8.2*kpc)**2:.4e} m/s^2,"
      f" a_N/a0 = {GM_MW/(8.2*kpc)**2/a0:.4f}")

print("\n  Frequenzverhaeltnis w_rad/w_azim fuer v^2 = GM/r + V^2:")
print("   (w_rad/w_azim)^2 = (GM/r + 2V^2)/(GM/r + V^2);  Kepler -> 1, V dominant -> sqrt(2)")
print("   Resonanzradien w_rad/w_azim = m/n (n<=5) im Intervall (1, sqrt2):  r = (GM/V^2)(p^2-1)/(2-p^2)")
fr = sorted({(m, n) for n in range(1, 6) for m in range(n, 2 * n)
             if 1 < m / n < math.sqrt(2) and math.gcd(m, n) == 1}, key=lambda x: x[0] / x[1])
for name, L in [("B", L_B), ("A", L_A)]:
    V2 = GM_MW / (3 * L)
    print(f"   Milchstrasse, Variante {name} (GM/V^2 = {GM_MW/V2/kpc:.4f} kpc):")
    for m, n in fr:
        p2 = (m / n) ** 2
        x = (2 - p2) / (p2 - 1)
        r = GM_MW / (V2 * x)
        print(f"     {m}/{n} = {m/n:.5f}:  r = {r/kpc:9.4f} kpc,  a_N/a0 = {GM_MW/r**2/a0:.4f}")
    rs = math.sqrt(GM_MW / a0)
    p = math.sqrt((GM_MW / rs + 2 * V2) / (GM_MW / rs + V2))
    print(f"     an r* = {rs/kpc:.4f} kpc: w_rad/w_azim = {p:.6f}")
print(f"   Erdbahn (B): w_rad/w_azim = {row(GMsun, AU, L_B_S)['ratio']:.6f}")

H("DIMENSIONSANALYSE (SI-Exponenten m, s)")
dims = {
    "r_ddot": (1, -2), "rdot": (1, -1), "dt'": (0, 1), "grad Phi": (1, -2),
    "kappa0": (0, -2), "M_k": (0, 0), "lam": (0, -1), "K": (0, -2),
    "int K rdot dt'": (1, -2), "Ds=GM/3c^2": (1, 0), "q=Ds/r": (0, 0),
    "||K_F|| = kappa0 q/lam": (0, -1), "a_eff = ||K_F|| v": (1, -2),
    "3c^2/r*": (1, -2), "GM/(r r*)": (1, -2), "sqrt(GM a0)/r": (1, -2), "v_rot = sqrt(a_eff r)": (1, -1),
    "r*=sqrt(GM/a0)": (1, 0), "gamma": (1, -2),
}
for k, (m, s) in dims.items():
    print(f"   [{k:26s}] = m^{m} s^{s}")
# Kontrolle Produktregeln
def mul(*xs): return (sum(x[0] for x in xs), sum(x[1] for x in xs))
print("   Kontrollen:")
print("    K * rdot * dt'         =", mul(dims["K"], dims["rdot"], dims["dt'"]), " (soll (1,-2))")
print("    kappa0 * M_k           =", mul(dims["kappa0"], dims["M_k"]), " (soll (0,-2))")
print("    (kappa0 q/lam) * v     =", mul(dims["kappa0"], (0, 0), (0, 1), dims["rdot"]), " (soll (1,-2))")
print("    GM/(3c^2): m^3 s^-2 / (m^2 s^-2) =", (3 - 2, -2 + 2), " (soll (1,0))")
print("    (Ds/r)(3c^2/r*)        =", mul((0, 0), (2, -2), (-1, 0)), " (soll (1,-2))")
print("    (GM a0)^(1/4): (m^3 s^-2 * m s^-2)^(1/4) = (m^4 s^-4)^(1/4) = m s^-1")

H("STRUKTUR STATT MASSE: Beobachterschema (Reifikation) formal")
print("  (1) global: chi 0->1, P aktiv;  (2) lokal gemessen: delta_a = a_obs - GM/r^2 = sqrt(GM a0)/r  (B)")
print("  (3) reifizierte Lesart: M_schein(r) = r^2 a_obs / G = M (1 + r/r*);  M_schein - M = M r/r* = v_rot^2 r / G")
print("  Feldgleichung: delta_a = q(r) * 3c^2/r*  (Strukturzahl q mal Regelkonstante 3c^2/r*), keine Groesse in kg")
print("  Ruhetest: rdot == 0  =>  int K rdot dt' = 0  (Gedaechtnisterm wirkt nur auf Bewegung)")
for name, L in [("B", L_B), ("A", L_A)]:
    print(f"   Milchstrasse {name}: L = {L/kpc:.4f} kpc, GM/v_rot^2 = 3L = {3*L/kpc:.4f} kpc")
    for R in radii:
        r = R * kpc
        da = GM_MW / (3 * r * L)
        q = ds(GM_MW) / r
        print(f"     r={R:7.4g} kpc: delta_a = {da:.4e} m/s^2 = q*c^2/L = {q*c**2/L:.4e};"
              f"  (M_schein-M)/M = r/(3L) = {r/(3*L):.4f};  M_schein-M = {GM_MW*r/(3*L)/G/Msun:.4e} Msun")
