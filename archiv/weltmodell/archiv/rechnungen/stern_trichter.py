"""Der Stern als Faltungstrichter: Stern-Faltung, Photosphäre und Carter im FTM-Längenbild."""
import math
G, c, hbar, kB, eV = 6.67430e-11, 299792458.0, 1.054571817e-34, 1.380649e-23, 1.602176634e-19
mp, me, alpha = 1.67262192e-27, 9.1093837e-31, 7.2973525693e-3
MS, RS, LS, TS = 1.989e30, 6.957e8, 3.828e26, 5772.0
lp, le = hbar/(mp*c), hbar/(me*c); ds_p = G*mp/(3*c**2); r_e = alpha*le
aG = G*mp**2/(hbar*c); mPl = math.sqrt(hbar*c/G)

print("1) Welche Faltungstiefe hat ein Stern? Die Sternmassen-Skala M_* = m_Pl³/m_p² (Landau/Chandrasekhar)")
M_star = mPl**3/mp**2
ds_star = G*M_star/(3*c**2)
print(f"   M_* = {M_star:.3e} kg = {M_star/MS:.2f} M☉,   Δs_* = {ds_star:.1f} m")
print(f"   FTM-Längenform:  Δs_*² · Δs_p = (λ̄_p/3)³   ->  {ds_star**2*ds_p:.4e} = {(lp/3)**3:.4e}")
print(f"   Gleichwertig: 3Δs_*/λ̄_p = α_G^(-1/2) = {3*ds_star/lp:.4e}  (Proton: 3Δs_p/λ̄_p = α_G^(+1) = {3*ds_p/lp:.3e})")
for n, M in [("Sonne", MS), ("Chandrasekhar-Grenze 1,44 M☉", 1.44*MS)]:
    print(f"   {n:30s} Δs = {G*M/(3*c**2):7.1f} m")
print("   -> Sterne sind die Trichter, deren Faltung die geometrische Mitte zwischen Proton-Faltung und Compton-Länge trifft.")

print("\n2) Photosphäre eines Sterns in Faltungsgrößen (Größenordnung, Thomson-Opazität, T_c ~ α² m_p c²):")
print("   Hydrostatik (parametrisch): kT_c ~ G M m_p/R = 3c²Δs·m_p/R")
print("   Strahlungstransport (tangential/transversal): σ_T = (8π/3)·r_e²  (EM-Radius des Elektrons!)")
print("   L ~ T_c⁴ R⁴ m_p/(σ_T M),  T_s⁴ = L/(4πR²σ_SB)  =>  T_s⁴ ~ α²·m_p³·m_e²/m_Pl   (mit M ~ M_*, T_c ~ α² m_p)")
T_x = mp*c**2*aG**0.25                      # Kreuzungstemperatur Proton (Notiz 19/20)
T_s_scale = math.sqrt(alpha*me*c**2 * T_x)
print(f"   => kT_s ~ √(α·m_e c² · kT_×)  mit kT_× = m_p c² α_G^(1/4) = {T_x/eV:.4f} eV  (Proton-Kreuzung, Rekombination)")
print(f"      √(α m_e c² · kT_×) = {T_s_scale/eV:.1f} eV   (reine Skala ohne Vorfaktoren)")

print("\n3) Die Photosphäre liegt an der Ionisationsgrenze (Übergang Strahlung/Konvektion) genau dann, wenn:")
print("   kT_s ~ α² m_e c²  ⟺  kT_× ~ α³ m_e c²  ⟺  α_G ~ α¹²(m_e/m_p)⁴  (Carter)")
print(f"   α² m_e c² = {alpha**2*me*c**2/eV:.2f} eV;  α³ m_e c² = {alpha**3*me*c**2/eV:.4f} eV;  kT_× = {T_x/eV:.4f} eV")
print("   -> Dieselbe Temperatur kT_× regelt (a) wann das Proton kosmisch zum Einzeltrichter wird und")
print("      (b) über √(α m_e c² · kT_×) die Temperaturskala der Stern-Oberflächen.")

print("\n4) Ehrliche Vorfaktor-Prüfung an der Sonne (Größenordnung):")
Tc_vir = G*MS*mp/(2*RS)/kB
sigma_T = 8*math.pi/3*r_e**2
a_rad = 4*5.670374e-8/c
rho_mean = MS/(4/3*math.pi*RS**3)
L_T = (16*math.pi*a_rad*c/3)*(Tc_vir**4)*RS/(3*(sigma_T/mp)*rho_mean) * 0.1   # grobe Mittelung (Faktor ~0,1 für Gradient)
Ts_T = (L_T/(4*math.pi*RS**2*5.670374e-8))**0.25
print(f"   Virialtemperatur T_c ≈ {Tc_vir:.3e} K  (Modell der Sonne: 1,57e7 K)")
print(f"   L mit Thomson-Opazität ≈ {L_T:.2e} W  (gemessen {LS:.2e} W)  ->  T_s ≈ {Ts_T:.0f} K  (gemessen {TS:.0f} K)")
print("   -> Die Skalen stimmen in der Größenordnung; die Vorfaktoren (reale Opazität, Gradienten, Kernfusion)")
print("      sind mit dieser Grobrechnung nicht genau zu treffen.")
