# -*- coding: utf-8 -*-
"""
FTM - Licht und Trichter (Punkte 1-3 des Autors)
Rechnungen zu 08-licht-und-trichter.md. Nur Python-Standardbibliothek.
Aufruf:  python licht_trichter.py > licht_trichter.txt

Formeln des Autors:
  GM        = 3 c^2 Ds
  c_eff(r)  = c * dtau/dt = c * sqrt(1 - 6 Ds / r)
  r0        = 6 Ds                       (Grenzradius, c_eff -> 0)
  g(r)      = 3 c^2 Ds / r^2 = g(R) (R/r)^2
  delta(r)  = 1 - sqrt(1 - 6 Ds / r)  ~  3 Ds / r = g(r) r / c^2
  r*        = sqrt(GM / a0) = sqrt(3 c^2 Ds / a0)
"""
import math, sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ---------------------------------------------------------------- Konstanten
# Einheitliche Konstanten des Weltmodells (01-08)
G     = 6.674e-11              # m^3 kg^-1 s^-2
c     = 299792458.0            # m/s
M_E   = 5.972e24               # kg  Erde
M_S   = 1.989e30               # kg  Sonne
GM_E  = G * M_E                # m^3/s^2  Erde
GM_S  = G * M_S                # m^3/s^2  Sonne
GM_MW = 6e10 * GM_S            # m^3/s^2  Milchstrasse, baryonisch 6e10 M_S (wie 03, 07)
GM_NS = 1.4 * GM_S             # m^3/s^2  Neutronenstern 1,4 M_S (wie 01)
R_E   = 6.3710e6               # m  mittlerer Erdradius (wie 01)
R_S   = 6.957e8                # m  Sonnenradius 695 700 km (wie 01)
R_NS  = 1.2e4                  # m  Neutronenstern-Radius (wie 01)
AU    = 1.495978707e11         # m
LY    = 9.4607304725808e15     # m
pc    = 3.0856775814913673e16  # m
kpc   = 1e3 * pc
R_MW  = 50000 * LY             # m  Aussenradius Milchstrasse (wie 03)
a0    = 1.2e-10                # m/s^2  Schwellenbeschleunigung
JAHR  = 365.25 * 86400.0       # s
r_ISS = 6.771e6                # m
r_Mond= 3.844e8                # m
r_Nep = 30.07 * AU             # m  (wie 07)
v0_A  = 220e3                  # m/s  Variante A (Vergleich) aus 07

def H(t): print("\n" + "=" * 90 + "\n" + t + "\n" + "=" * 90)
def Ds_of(GM): return GM / (3 * c**2)
def g_of(GM, r): return GM / r**2
def rate(GM, r): return math.sqrt(1 - 6 * Ds_of(GM) / r)
def delta(GM, r):
    # numerisch stabile Form von 1 - sqrt(1-u):  u / (1 + sqrt(1-u))
    u = 6 * Ds_of(GM) / r
    return u / (1 + math.sqrt(1 - u))
def delta1(GM, r): return 3 * Ds_of(GM) / r
def ceff(GM, r): return c * rate(GM, r)
def rstar(GM): return math.sqrt(GM / a0)
def r0(GM): return 6 * Ds_of(GM)

# ================================================================ PUNKT 1
H("PUNKT 1a: GM = 3 c^2 Ds  -  c als Skalierungsfaktor")
print(f"  c                 = {c:.0f} m/s")
print(f"  c^2               = {c**2:.12e} m^2/s^2")
print(f"  3 c^2             = {3*c**2:.12e} m^2/s^2")
print(f"  Ds/M = G/(3c^2)   = {G/(3*c**2):.6e} m/kg")
for n, GM in [("Erde", GM_E), ("Sonne", GM_S), ("Milchstrasse", GM_MW), ("Neutronenstern", GM_NS)]:
    Ds = Ds_of(GM)
    print(f"  {n:16s} GM = {GM:.6e} m^3/s^2   Ds = {Ds:.6e} m   GM/(3 Ds) = {GM/(3*Ds):.12e} = c^2")
gE = g_of(GM_E, R_E)
print(f"\n  Probe Erde: g(R_E) R_E^2 / (3 Ds_E) = {gE*R_E**2/(3*Ds_of(GM_E)):.12e} m^2/s^2  (c^2 = {c**2:.12e})")

H("PUNKT 1b: Lokale Messung  c_lokal = c_eff / (dtau/dt)")
for n, GM, r in [("Erdoberflaeche", GM_E, R_E), ("Sonnenoberflaeche", GM_S, R_S), ("Neutronenstern", GM_NS, R_NS)]:
    print(f"  {n:18s} c_eff = {ceff(GM, r):.6f} m/s   dtau/dt = {rate(GM, r):.15f}   c_eff/(dtau/dt) = {ceff(GM, r)/rate(GM, r):.6f} m/s")

H("PUNKT 1c: Grenzradius r0 = 6 Ds (Nullstelle von c_eff)")
for n, GM, R in [("Erde", GM_E, R_E), ("Sonne", GM_S, R_S), ("Milchstrasse", GM_MW, R_MW), ("Neutronenstern", GM_NS, R_NS)]:
    x0 = r0(GM)
    print(f"  {n:16s} Ds = {Ds_of(GM):.6e} m   r0 = {x0:.6e} m = {x0/AU:.6e} AU   R = {R:.6e} m   R/r0 = {R/x0:.6e}")
    print(f"  {'':16s} g(r0) = c^2/(12 Ds) = {c**2/(12*Ds_of(GM)):.6e} m/s^2   Probe GM/r0^2 = {GM/x0**2:.6e}")

print("\n  Verhalten bei r -> r0,  r = r0 (1 + xi):  dtau/dt = sqrt(xi/(1+xi)),  c_eff ~ c sqrt(xi)")
print(f"  {'xi':>10s} {'dtau/dt':>16s} {'sqrt(xi)':>14s} {'c_eff [m/s]':>16s} {'delta':>14s} {'3Ds/r':>12s}")
GM = GM_S
for xi in [10.0, 1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-6, 0.0]:
    r = r0(GM) * (1 + xi)
    q = 1 - 6 * Ds_of(GM) / r
    rt = math.sqrt(max(q, 0.0))
    print(f"  {xi:>10.0e} {rt:>16.10f} {math.sqrt(xi):>14.8f} {c*rt:>16.4f} {1-rt:>14.10f} {3*Ds_of(GM)/r:>12.8f}")
print("  Ableitung: d c_eff/dr = c (3 Ds / r^2) / sqrt(1 - 6 Ds/r)  ->  unendlich fuer r -> r0")
print("  An r0: 6 Ds/r0 = 1, dtau/dt = 0, delta = 1, 3 Ds/r0 = 1/2, g(r0) r0 / c^2 = 1/2")

H("PUNKT 1d: Topologische Abkuerzung  t_aussen = d_aussen / c,  d_innen = 0")
for n, d in [("Erde - Mond", r_Mond), ("Sonne - Erde (1 AU)", AU), ("Sonne - r* der Sonne", rstar(GM_S)),
             ("100 Lichtjahre", 100*LY), ("Durchmesser Milchstrasse (100 000 Lj)", 1e5*LY)]:
    t = d / c
    print(f"  {n:40s} d_aussen = {d:.6e} m   t_aussen = {t:.6e} s = {t/JAHR:.6e} a   d_innen = 0   t_innen = 0")

# ================================================================ PUNKT 2
H("PUNKT 2: Reihenentwicklung von delta(r) = 1 - sqrt(1 - 6 Ds/r)")
print("  1 - sqrt(1-u) = sum_{n>=1} C(2n,n) u^n / ((2n-1) 4^n),  u = 6 Ds/r")
print("  delta = sum a_n (Ds/r)^n,   a_n = C(2n,n) (3/2)^n / (2n-1)")
print("  delta = sum b_n eps^n,      eps = g r / c^2 = 3 Ds/r,  b_n = C(2n,n) / ((2n-1) 2^n)")
from fractions import Fraction
for n in range(1, 7):
    an = Fraction(math.comb(2*n, n)) * Fraction(3, 2)**n / (2*n - 1)
    bn = Fraction(math.comb(2*n, n), (2*n - 1) * 2**n)
    un = Fraction(math.comb(2*n, n), (2*n - 1) * 4**n)
    print(f"   n={n}:  Koeff. von u^n = {str(un):>8s}   a_n = {str(an):>10s} = {float(an):10.4f}   b_n = {str(bn):>6s}")
s = sum(Fraction(math.comb(2*n, n), (2*n - 1) * 4**n) for n in range(1, 4001))
print(f"  Summe der u^n-Koeffizienten bis n=4000 (u=1, also r=r0): {float(s):.6f}  (Grenzwert 1; Rest ~ 1/sqrt(pi n))")
print(f"  1/sqrt(pi*4000) = {1/math.sqrt(math.pi*4000):.6f}")

print("\n  Pruefung der Reihe an Neutronenstern (6Ds/R gross) und Sonnenoberflaeche:")
for n, GM, r in [("Neutronenstern", GM_NS, R_NS), ("Sonnenoberflaeche", GM_S, R_S)]:
    k = Ds_of(GM) / r
    ex = delta(GM, r)
    partial = 0.0
    for N in range(1, 7):
        an = math.comb(2*N, N) * 1.5**N / (2*N - 1)
        partial += an * k**N
        print(f"   {n:18s} bis n={N}: {partial:.12e}   exakt {ex:.12e}   rel. Rest {(ex-partial)/ex:+.3e}")

# ================================================================ PUNKT 2 Tabelle
H("PUNKT 2C: Tabelle g(r), delta(r), Zeitverzug pro Jahr, c_eff(r)")
rows = [
    ("Erdoberflaeche (R_E)",            GM_E,  R_E,          "6 371 km"),
    ("Erd-Orbit ISS",                   GM_E,  r_ISS,        "6 771 km (+400 km)"),
    ("Mondbahn (Erdtrichter)",          GM_E,  r_Mond,       f"384 400 km = {r_Mond/R_E:.2f} R_E"),
    ("Sonnenoberflaeche (R_S=695700km)",GM_S,  R_S,          "695 700 km"),
    ("Erdbahn im Sonnentrichter (1 AU)",GM_S,  AU,           f"1 AU = {AU/R_S:.2f} R_S"),
    ("Neptunbahn (30,07 AU)",           GM_S,  r_Nep,        f"30,07 AU = {r_Nep/R_S:.1f} R_S"),
    ("Sonne bei r*",                    GM_S,  rstar(GM_S),  f"{rstar(GM_S)/AU:.1f} AU"),
    ("Neutronenstern (1,4 M_S, 12 km)", GM_NS, R_NS,         "12 km"),
    ("Milchstrasse Oberflaeche",        GM_MW, R_MW,         f"50 000 Lj = {R_MW/kpc:.3f} kpc"),
    ("Milchstrasse bei r*",             GM_MW, rstar(GM_MW), f"{rstar(GM_MW)/kpc:.3f} kpc"),
    ("Milchstrasse bei 50 kpc",         GM_MW, 50*kpc,       "50 kpc"),
]
print(f"  {'Ort':34s} {'r [m]':>12s} {'g [m/s^2]':>13s} {'delta exakt':>13s} {'3Ds/r':>13s} {'g r/c^2':>13s} "
      f"{'Verzug/Jahr [s]':>15s} {'dtau/dt':>20s} {'c-c_eff=c*delta':>14s}")
for n, GM, r, lab in rows:
    g = g_of(GM, r); d = delta(GM, r); d1 = delta1(GM, r)
    print(f"  {n:34s} {r:12.5e} {g:13.5e} {d:13.5e} {d1:13.5e} {g*r/c**2:13.5e} {d*JAHR:15.5e} "
          f"{rate(GM, r):20.15f} {c*delta(GM, r):14.6e}")
    print(f"  {'':34s} ({lab};  c_eff = {ceff(GM, r):.6f} m/s)")

print("\n  Probe g(r) = g(R) (R/r)^2 fuer die Erde:")
for r in [r_ISS, r_Mond]:
    print(f"   r = {r:.4e} m:  g(R)(R/r)^2 = {gE*(R_E/r)**2:.6e}   3c^2 Ds/r^2 = {3*c**2*Ds_of(GM_E)/r**2:.6e}")

print("\n  Zeitverzug in anschaulichen Einheiten:")
for n, GM, r, lab in rows:
    d = delta(GM, r) * JAHR
    print(f"   {n:34s} {d:.6e} s = {d*1e3:.6f} ms = {d/86400:.6e} d")

# ================================================================ PUNKT 3
H("PUNKT 3: Schwellenradius r* = sqrt(GM/a0) = sqrt(3 c^2 Ds / a0) = c sqrt(3 Ds / a0)")
print(f"  a0 = {a0:.3e} m/s^2")
print(f"  {'System':16s} {'Ds [m]':>12s} {'r* [m]':>13s} {'r* [AU]':>11s} {'r* [pc]':>11s} {'r* [Lj]':>11s} "
      f"{'r0=6Ds [m]':>12s} {'r*/r0':>12s} {'sqrt(g(r0)/a0)':>14s} {'delta(r*)':>11s} {'a0 r*/c^2':>11s}")
for n, GM in [("Erde", GM_E), ("Sonne", GM_S), ("Milchstrasse", GM_MW), ("Neutronenstern", GM_NS)]:
    Ds = Ds_of(GM); rs = rstar(GM); x0 = r0(GM)
    print(f"  {n:16s} {Ds:12.5e} {rs:13.6e} {rs/AU:11.4e} {rs/pc:11.5f} {rs/LY:11.5f} {x0:12.5e} {rs/x0:12.5e} "
          f"{math.sqrt(c**2/(12*Ds)/a0):14.5e} {delta(GM, rs):11.5e} {a0*rs/c**2:11.5e}")
    print(f"  {'':16s} Probe c*sqrt(3Ds/a0) = {c*math.sqrt(3*Ds/a0):.6e} m;  r*/r0 = c/sqrt(12 a0 Ds) = {c/math.sqrt(12*a0*Ds):.6e};"
          f"  q* = Ds/r* = {Ds/rs:.5e};  c - c_eff(r*) = {c*delta(GM, rs):.6e} m/s;  a0 r*/c = {a0*rs/c:.6e} m/s")
Ds_gl = c**2 / (12 * a0)
print(f"\n  r*/r0 = 1  bei  Ds = c^2/(12 a0) = {Ds_gl:.5e} m,  GM = 3c^2 Ds = {3*c**2*Ds_gl:.5e} m^3/s^2 = {3*c**2*Ds_gl/GM_S:.5e} M_S")
print(f"  Milchstrasse r* = {rstar(GM_MW)/kpc:.4f} kpc = {rstar(GM_MW)/LY:.1f} Lj;  Sonne r* = {rstar(GM_S)/AU:.2f} AU = {rstar(GM_S)/LY:.5f} Lj = {rstar(GM_S)/pc:.5f} pc")

H("PUNKT 3: Uebergang Einzeltrichter -> Netzwerk-Faltungsfeld (Rueckkopplungssumme a = g/(1-q), q = exp(-r*/r))")
def nu_of(GM, r):
    q = math.exp(-rstar(GM) / r)
    return q, 1.0 / (1.0 - q)
print(f"  An der Schwelle: q(r*) = e^-1 = {math.exp(-1):.5f},  nu(r*) = 1/(1-e^-1) = {1/(1-math.exp(-1)):.5f}")

# 3.4 Sonne: kompakter Trichter im Feld der Milchstrasse -> Newtonsch bis zur Gezeitengrenze r_J
A_oort, B_oort = 15.3e3 / kpc, -11.9e3 / kpc          # 1/s
r_J = (GM_S / (4 * A_oort * (A_oort - B_oort)))**(1/3)
rs_S = rstar(GM_S)
print(f"\n  Sonne: g(R_Sonne) = {GM_S/(6.957e8)**2:.3f} m/s^2;  r* = {rs_S/AU:.2f} AE;  "
      f"Gezeitengrenze r_J = {r_J/AU:.3e} AE = {r_J/pc:.3f} pc = {r_J/rs_S:.1f} r*")
print(f"  {'Ort':>22s} {'r/r*':>11s} {'g [m/s^2]':>13s} {'g/a0':>11s} {'a = g':>13s} {'v = v_N [m/s]':>14s}")
for n, r in [("Erdbahn (1 AE)", AU), ("Neptunbahn", r_Nep), ("1 000 AE", 1000*AU), ("0,5 r*", 0.5*rs_S),
             ("r*", rs_S), ("2 r*", 2*rs_S), ("10 r*", 10*rs_S), ("100 000 AE", 1e5*AU)]:
    g = g_of(GM_S, r)
    print(f"  {n:>22s} {r/rs_S:11.4e} {g:13.5e} {g/a0:11.4e} {g:13.5e} {math.sqrt(GM_S/r):14.2f}")

# 3.5 Milchstrasse: ausgedehnt (g(R) < a0) -> Uebergang am eigenen r*
rs_MW = rstar(GM_MW)
print(f"\n  Milchstrasse: r* = {rs_MW/kpc:.4f} kpc,  v_rot = (GM a0)^(1/4) = {(GM_MW*a0)**0.25/1e3:.3f} km/s")
print(f"  {'r [kpc]':>9s} {'r/r*':>8s} {'g':>12s} {'g/a0':>9s} {'q':>8s} {'nu':>9s} {'a = nu g':>12s} {'v_N':>9s} {'v':>9s}")
for r in [2*kpc, 5*kpc, rs_MW, 10*kpc, R_MW, 20*kpc, 50*kpc]:
    g = g_of(GM_MW, r); q, nu = nu_of(GM_MW, r); vN = math.sqrt(GM_MW / r)
    print(f"  {r/kpc:9.4f} {r/rs_MW:8.4f} {g:12.5e} {g/a0:9.5f} {q:8.5f} {nu:9.5f} {nu*g:12.5e} {vN/1e3:9.3f} {math.sqrt(nu)*vN/1e3:9.3f}")
print(f"  Netzwerkterm bei 50 kpc: (nu-1) = {nu_of(GM_MW, 50*kpc)[1]-1:.2f}-fach der Einzeltrichter-Neigung")

H("Einheiten-Kontrolle")
print("  [3c^2 Ds]   = m^2 s^-2 * m        = m^3 s^-2 = [GM]")
print("  [c_eff]     = m s^-1 * 1          = m s^-1")
print("  [g]         = m^2 s^-2 * m / m^2  = m s^-2")
print("  [delta]     = 1 ;  [g r / c^2] = m s^-2 * m / (m^2 s^-2) = 1")
print("  [r*]        = sqrt(m^3 s^-2 / (m s^-2)) = m ;  [r*/r0] = 1")
