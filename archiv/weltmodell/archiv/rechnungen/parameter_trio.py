"""Parameter-Trio (κ0, λ, r0) und Festlegung A vs. B."""
import math
G, c = 6.674e-11, 299792458.0
MS, KPC, MPC = 1.989e30, 3.0857e19, 3.0857e22
H0 = 70e3/MPC; RH = c/H0; a0_obs = 1.2e-10; a0_th = c*H0/(2*math.pi)

print("1) A gegen B für die Milchstraße")
M = 6e10*MS; v0 = 220e3
rA = G*M/(3*v0**2); rB = math.sqrt(G*M/a0_obs)/3
print(f"   A: r0 = GM/(3 v0²)       = {rA/KPC:.3f} kpc   (v0 = 220 km/s, M = 6e10 M☉)")
print(f"   B: r0 = (1/3)·sqrt(GM/a0) = {rB/KPC:.3f} kpc")
print("   Identität: setzt man v0⁴ = G·M·a0 ein, wird A zu B:")
v_B = (G*M*a0_obs)**0.25
print(f"      v0 aus v⁴ = GMa0 bei M = 6e10 M☉: {v_B/1e3:.1f} km/s  ->  A = {G*M/(3*v_B**2)/KPC:.3f} kpc = B")
M_A = v0**4/(G*a0_obs)
print(f"      M aus v⁴ = GMa0 bei v0 = 220 km/s: {M_A/MS:.3e} M☉  ->  B = {math.sqrt(G*M_A/a0_obs)/3/KPC:.3f} kpc = A (mit dieser Masse)")
print("   -> A und B sind dieselbe Formel r0 = r*/3 = GM/(3 v_flach²). Der Unterschied liegt nur im Datenpaar (M, v0).")

print("\n2) r0 aus Naturkonstanten und Torus-Atmen: r0 = r*/3 = (1/3)·sqrt(6π·Δs·R_H),  a0 = c²/(2π R_H) = c·H0/2π")
print(f"   R_H = c/H0 = {RH/MPC:.1f} Mpc,  c²/(2π R_H) = {c*c/(2*math.pi*RH):.4e} m/s²")
for m in [1e9, 1e10, 6e10, 1.47e11, 1e12]:
    Mx = m*MS; ds = G*Mx/(3*c*c)
    r0 = math.sqrt(6*math.pi*ds*RH)/3
    print(f"   M = {m:8.2e} M☉   Δs = {ds:.3e} m   r0 = {r0/KPC:6.3f} kpc   v_flach = {(G*Mx*c*H0/(2*math.pi))**0.25/1e3:6.1f} km/s")

print("\n3) κ0 und λ: stationäres Feldgedächtnis g = g_N + ∫κ0 e^{-λs} g ds  =>  g = g_N/(1 - κ0/λ)")
print("   Vergleich mit g = g_N/(1-q):  κ0/λ = q(r) = exp(-r*/r)  (nur das Verhältnis ist festgelegt)")
for rk in [2,5,8.35,20,50]:
    rs = math.sqrt(G*M/a0_obs); q = math.exp(-rs/(rk*KPC))
    print(f"   Milchstraße r = {rk:5.2f} kpc   κ0/λ = q = {q:.4f}")
print("   λ allein bestimmt nur die Nachlaufzeit 1/λ nach einer Änderung der Massenverteilung:")
for name, lam in [("λ = H0 (Atmen)", H0), ("λ = 1/(100 Mio a)", 1/(1e8*3.156e7)), ("λ = 1/(1 Mrd a)", 1/(1e9*3.156e7))]:
    print(f"   {name:22s} Nachlaufzeit 1/λ = {1/lam/3.156e7/1e9:8.3f} Mrd Jahre")
