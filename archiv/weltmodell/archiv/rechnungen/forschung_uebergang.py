"""Forschung: Variante B der Feldgleichung, Sprung am Schwellenradius r*, glatte Übergänge, a0 aus dem Atmen."""
import math
G, c = 6.674e-11, 2.998e8
MS, AU, KPC, Mpc = 1.989e30, 1.496e11, 3.0857e19, 3.0857e22
a0 = 1.2e-10

print("1) Variante B: r0 = r*/3  ->  a_eff = GM/(3 r r0) = sqrt(GM a0)/r,  v_rot^4 = G M a0")
for M_sol in [1e9, 1e10, 6e10, 2e11]:
    M = M_sol*MS; rs = math.sqrt(G*M/a0); r0 = rs/3
    v = math.sqrt(G*M/(3*r0))
    print(f"   M = {M_sol:8.1e} M☉   r* = {rs/KPC:6.2f} kpc   v_rot = {v/1e3:6.1f} km/s   v^4/(G M a0) = {v**4/(G*M*a0):.6f}")

M = 6e10*MS; rs = math.sqrt(G*M/a0)
print(f"\n2) Harter Schalter χ an r* = {rs/KPC:.2f} kpc (Milchstraße, 6e10 M☉):")
for f in [0.999, 1.001]:
    r = f*rs; aN = G*M/r**2; chi = 1 if aN < a0 else 0
    a = aN + chi*math.sqrt(G*M*a0)/r
    print(f"   r = {f:.3f} r*   a = {a:.3e} m/s²   v = {math.sqrt(a*r)/1e3:6.1f} km/s")
print("   -> Sprung um Faktor 2 in a, um √2 in v.")

print("\n3) Glatte Übergänge, alle mit Newton für a_N >> a0 und sqrt(a_N a0) für a_N << a0:")
def rar(y):   return 1/(1-math.exp(-math.sqrt(y)))            # empirische Kurve (Datenbezug)
def feedback(y):                                               # Rückkopplung: a = a_N + q a, q = exp(-sqrt(y))
    q = math.exp(-math.sqrt(y)); return 1/(1-q)
def add_smooth(y):                                             # a = a_N + sqrt(a_N a0)·(1 - exp(-1/sqrt(y)))
    return 1 + (1-math.exp(-1/math.sqrt(y)))/math.sqrt(y)*1.0
print(f"   {'y=a_N/a0':>9s} {'empirisch':>10s} {'Rückkopplung':>13s} {'additiv glatt':>14s} {'Schalter':>9s}")
for y in [100,10,3,1,0.3,0.1,0.01]:
    sw = 1 + (1/math.sqrt(y) if y < 1 else 0)
    print(f"   {y:9.2f} {rar(y):10.3f} {feedback(y):13.3f} {add_smooth(y):14.3f} {sw:9.3f}")
print("   Rückkopplung Σqⁿ ist mit der empirischen Kurve identisch (gleiche Funktion).")

print("\n4) Schwelle a0 und das Atmen des Torus: a0 = c·H0/(2π)")
for H in [67.4, 70, 73]:
    H0 = H*1e3/Mpc
    print(f"   H0 = {H:5.1f} km/s/Mpc  ->  c·H0/2π = {c*H0/(2*math.pi):.3e} m/s²   (a0 = 1.2e-10)")

print("\n5) Sonnensystem mit Rückkopplungs-Übergang (q = exp(-sqrt(a_N/a0))):")
for name, r in [("Erdbahn",AU),("Saturn",9.58*AU),("Neptun",30.1*AU),("Sonne r*", math.sqrt(G*MS/a0))]:
    aN = G*MS/r**2; y = aN/a0; q = math.exp(-math.sqrt(y))
    print(f"   {name:9s} a_N/a0 = {y:10.3e}  q = {q:.2e}  Zusatz = {aN*q/(1-q):.2e} m/s²")

print("\n6) Ramanujan-Δ auf der imaginären Achse τ = i·t (q = e^{-2πt}), Δ(i/t) = t^12 Δ(it):")
def Delta(t):
    if t < 1: return Delta(1/t)/t**12
    q = math.exp(-2*math.pi*t); p = q
    for n in range(1, 400): p *= (1-q**n)**24
    return p
for t in [3, 2, 1, 0.5, 0.3, 0.2, 0.1]:
    print(f"   t = {t:4.1f}   q = {math.exp(-2*math.pi*t):.3e}   Δ = {Delta(t):.3e}")
print("   -> Δ ist an beiden Enden null: bei q→0 (Nahfeld) und bei q→1 (tiefes Fernfeld).")
print("   Geometrische Reihe Σqⁿ = q/(1-q) wächst dagegen für q→1 unbegrenzt.")

print("\n7) Kepler-Frequenzverhältnis: Punktmasse ω_rad/ω_azim = 1 (rational, geschlossene Bahnen),")
print(f"   flache Rotationskurve κ/Ω = √2 = {math.sqrt(2):.4f} (irrational).")
