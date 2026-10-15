"""Die Kopplungskonstanten als Faltungsverhältnisse; Carter als Längenbeziehung; die Kreuzung als Ratengleichheit."""
import math
G, c, hbar, kB, eV = 6.67430e-11, 299792458.0, 1.054571817e-34, 1.380649e-23, 1.602176634e-19
mp, me = 1.67262192e-27, 9.1093837e-31
alpha = 7.2973525693e-3
MPC = 3.0857e22; H0 = 67.36e3/MPC; Om, Or, OL = 0.3153, 9.14e-5, 0.6846
H = lambda z: H0*math.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)

lp, le = hbar/(mp*c), hbar/(me*c)              # reduzierte Compton-Längen
ds_p = G*mp/(3*c**2)                          # Faltungstiefe des Protons
r_e = alpha*le                                 # klassischer Elektronenradius = α·λ̄_e
aG = G*mp**2/(hbar*c)

print("1) Die beiden Kopplungskonstanten als Längenverhältnisse")
print(f"   Gravitation:  3·Δs_p/λ̄_p = {3*ds_p/lp:.6e}   α_G = G m_p²/(ħc) = {aG:.6e}")
print(f"   Elektromagn.: r_e/λ̄_e    = {r_e/le:.10f}   α = {alpha:.10f}")
print(f"   Δs_p = {ds_p:.4e} m,  λ̄_p = {lp:.4e} m,  r_e = {r_e:.4e} m,  λ̄_e = {le:.4e} m")
print("   -> α_G ist die Faltungstiefe des Protons in Einheiten seiner Compton-Länge (×3);")
print("      α ist der elektromagnetische Radius des Elektrons in Einheiten seiner Compton-Länge.")

print("\n2) Carters Koinzidenz als reine Längenbeziehung")
lhs = 3*ds_p*le**16; rhs = r_e**12*lp**5
print(f"   3·Δs_p·λ̄_e¹⁶ = {lhs:.4e} m¹⁷    r_e¹²·λ̄_p⁵ = {rhs:.4e} m¹⁷    Verhältnis = {lhs/rhs:.3f}")
print("   (gleichwertig zu α_G = α¹²(m_e/m_p)⁴, Verhältnis 2,94)")

print("\n3) Die Kreuzung als Gleichheit zweier Raten")
wp = mp*c**2/hbar
GammaG = aG*wp
print(f"   Compton-Takt des Protons ω_p = m_p c²/ħ = {wp:.4e} 1/s")
print(f"   Faltungsrate Γ_G = α_G·ω_p = g_p(λ̄_p)/c = {GammaG:.4e} 1/s   (Probe G m_p³ c/ħ² = {G*mp**3*c/hbar**2:.4e})")
lo, hi = 0.0, 1e6
for _ in range(200):
    mid = (lo+hi)/2
    if H(mid)/(2*math.pi) < GammaG: lo = mid
    else: hi = mid
print(f"   Atemrate pro Radiant H/2π = Γ_G bei z = {lo:.1f}  (kT = {kB*2.7255*(1+lo)/eV:.4f} eV)")
print("   -> Die Kreuzung ist: Das Universum atmet pro Radiant so schnell, wie das Proton pro Compton-Takt faltet.")

print("\n4) Dasselbe für das Elektron und die Elektromagnetik (zum Vergleich)")
we = me*c**2/hbar
print(f"   Faltungsrate des Elektrons α_G,e·ω_e = G m_e³ c/ħ² = {G*me**3*c/hbar**2:.3e} 1/s  <  H_Λ/2π = {H0*math.sqrt(OL)/(2*math.pi):.3e} 1/s  -> nie")
print(f"   EM-Rate des Elektrons α·ω_e = {alpha*we:.3e} 1/s  (Größenordnung Atomphysik)")

print("\n5) Was die Kreuzung mit der Rekombination verbindet (Kette aus Notiz 19, in Raten):")
print(f"   ħ·Γ_G-Temperatur bei Kreuzung      : kT = {kB*2.7255*(1+lo)/eV:.4f} eV")
print(f"   m_p c²·α_G^(1/4)                   : {mp*c**2*aG**0.25/eV:.4f} eV")
print(f"   α³ m_e c²                          : {alpha**3*me*c**2/eV:.4f} eV")
print(f"   kT(z = 1090, letzte Streuung)      : {kB*2.7255*1091/eV:.4f} eV")
