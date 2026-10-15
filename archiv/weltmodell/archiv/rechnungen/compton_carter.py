"""Was steckt hinter der Proton-Kreuzung? Zerlegung in bekannte Kopplungen: α_G, α, m_e/m_p, η."""
import math
G, c, hbar, kB, eV = 6.67430e-11, 299792458.0, 1.054571817e-34, 1.380649e-23, 1.602176634e-19
mp, me = 1.67262192e-27, 9.1093837e-31
alpha = 7.2973525693e-3
aG = G*mp**2/(hbar*c)
mu = me/mp
eta = 6.13e-10

print("1) Die drei Energieskalen bei z ≈ 1100")
E_cross = mp*c**2*aG**0.25          # Temperatur, bei der H(T) ~ G m_p³/ħ² (Strahlung/Gleichheit)
E_atom3 = alpha**3*me*c**2          # α³ m_e c²
kT_rec  = kB*2.7255*1101            # Rekombination (z* ≈ 1090–1100)
print(f"   m_p c²·α_G^(1/4)  = {E_cross/eV:.4f} eV   (Gravitation + Nukleonmasse)")
print(f"   α³·m_e c²         = {E_atom3/eV:.4f} eV   (Elektromagnetismus + Elektronmasse)")
print(f"   kT bei z = 1100   = {kT_rec/eV:.4f} eV   (gemessene Rekombination)")

print("\n2) Carters Koinzidenz  α_G ≈ α¹²·(m_e/m_p)⁴  (Carter 1974; Bedingung für konvektive Sternhüllen)")
carter = alpha**12*mu**4
print(f"   α_G = {aG:.4e}   α¹²(m_e/m_p)⁴ = {carter:.4e}   Verhältnis = {aG/carter:.3f}")
print(f"   Vierte Wurzel des Verhältnisses = {(aG/carter)**0.25:.4f}  =  (m_p α_G^(1/4)) / (α³ m_e)  = {E_cross/E_atom3:.4f}")
print("   -> 'Proton-Neigung an der Compton-Grenze = Schwelle zur Rekombination' ist gleichbedeutend mit")
print("      m_p α_G^(1/4) ≈ α³ m_e, also mit Carters Koinzidenz (bis Faktor 1,3 in der Energie).")

print("\n3) Warum liegt die Rekombination bei ≈ α³ m_e c²?  Saha: kT_rec ≈ B/L, B = α² m_e c²/2 = 13,6 eV")
B = alpha**2*me*c**2/2
L = B/kT_rec
print(f"   B = {B/eV:.3f} eV,  L = B/kT_rec = {L:.1f};   1/(2α) = {1/(2*alpha):.1f}")
print(f"   ln(1/η) = {math.log(1/eta):.1f}  (+ 1,5·ln(m_e c²/kT) = {1.5*math.log(me*c**2/kT_rec):.1f}) -> L_Saha ≈ {math.log(1/eta)+1.5*math.log(me*c**2/kT_rec):.1f}")
print("   -> kT_rec ≈ α³ m_e c² gilt, weil der Saha-Logarithmus (≈ 50) zufällig nahe 1/(2α) (≈ 69) liegt.")

print("\n4) Warum ist die Neigungs-Schwelle bei der Rekombination m_p α_G^(1/4)?")
print("   Bedingung H(T) = 2π G m_p³ c/ħ²; nahe der Gleichheit Materie/Strahlung gilt H ~ T²/(m_Pl c²)·(c/ħ)·Zahl")
print("   => T² ~ m_p³/m_Pl  =>  T ~ m_p·(m_p/m_Pl)^(1/2) = m_p·α_G^(1/4)")
zeq = 3387; zrec = 1090
print(f"   Dazu gehört: Gleichheit (z_eq = {zeq}) liegt nahe der Rekombination (z = {zrec}); Verhältnis {(1+zeq)/(1+zrec):.2f}")

print("\n5) Zusammenfassung der Kette (jede Stufe mit ihrem Zahlenfaktor):")
print(f"   a) Kreuzung ⇔ T_cross ≈ m_p c² α_G^(1/4) · (Faktoren aus η, R, 2π)      T_cross = {kB*2.7255*1120/eV:.3f} eV (z=1119)")
print(f"   b) m_p α_G^(1/4) = α³ m_e · (α_G/α¹²μ⁴)^(1/4)                             Faktor {(aG/carter)**0.25:.3f}")
print(f"   c) α³ m_e c² ≈ kT_rec, weil L_Saha ≈ 1/(2α)                                Faktor {kT_rec/E_atom3:.3f}")
