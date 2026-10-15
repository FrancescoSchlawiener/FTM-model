"""Logische Folgerungen aus den FTM-Formeln: drei Radien, Einzeltrichter vs. kollektiv, Schwelle und das Ganze."""
import math
G, c = 6.674e-11, 299792458.0
ME, RE, gE = 5.972e24, 6.371e6, 9.81
MS, RS = 1.989e30, 6.957e8
AU, LY, PC, KPC, MPC = 1.496e11, 9.461e15, 3.0857e16, 3.0857e19, 3.0857e22
a0 = 1.2e-10; H0 = 70e3/MPC; RH = c/H0

def dist(m):
    if m < 1e3: return f"{m:.3g} m"
    if m < 1e9: return f"{m/1e3:.3g} km"
    if m < 1e15: return f"{m/AU:.3g} AE"
    if m < 1e19: return f"{m/LY:.3g} Lj"
    if m < 1e22: return f"{m/KPC:.3g} kpc"
    return f"{m/MPC:.3g} Mpc"

K = [("Erde",ME,RE),("Jupiter",1.898e27,6.9911e7),("Sonne",MS,RS),("Beteigeuze",16.5*MS,764*RS),
     ("Weißer Zwerg",1.02*MS,0.0084*RS),("Neutronenstern",1.4*MS,1.2e4),("Kugelsternhaufen",5e5*MS,5*PC),
     ("Zwerggalaxie",1e8*MS,1*KPC),("Milchstraße",6e10*MS,15.33*KPC),("Galaxienhaufen",1.7e14*MS,2*MPC),
     ("Universum",1.5e53,4.4e26)]
print("1) Drei Radien je Körper: r0 = 6Δs (Takt = 0), R (Oberfläche), r* = sqrt(3c²Δs/a0) (Schwelle)")
print(f"   {'Körper':17s} {'r0=6Δs':>11s} {'R':>11s} {'r*':>11s} {'g(R)/a0':>10s} {'Regime':>12s}")
for n,M,R in K:
    ds = G*M/(3*c*c); r0 = 6*ds; rs = math.sqrt(3*c*c*ds/a0); gR = G*M/R**2
    print(f"   {n:17s} {dist(r0):>11s} {dist(R):>11s} {dist(rs):>11s} {gR/a0:10.2e} {'kollektiv' if gR<a0 else 'Einzeltrichter':>12s}")
print("   Kriterium mit g(x,y) des Autors: kollektiv  <=>  x/y²·g⊕ < a0  <=>  y > sqrt(x·g⊕/a0)")
print(f"   Kritische Flächendichte: g = πGΣ < a0  <=>  Σ < a0/(πG) = {a0/(math.pi*G):.3f} kg/m² = {a0/(math.pi*G)/(MS/PC**2):.0f} M☉/pc²")

print("\n2) r* als geometrisches Mittel aus lokalem Horizont und kosmischer Länge:")
print("   r*² = GM/a0 = r0 · c²/(2 a0);  mit a0 = cH0/2π:  r* = sqrt(r0 · π·R_H) = sqrt(6π·Δs·R_H)")
for n,M,R in [K[2],K[8]]:
    ds=G*M/(3*c*c); aH=c*H0/(2*math.pi)
    print(f"   {n:12s} sqrt(6π Δs R_H) = {dist(math.sqrt(6*math.pi*ds*RH))}   sqrt(GM/(cH0/2π)) = {dist(math.sqrt(G*M/aH))}")

print("\n3) Takt-Abweichung an der Schwelle: δ(r*) = 3Δs/r* = sqrt(GM a0)/c² = (v_flach/c)²")
for n,M,R in [K[2],K[7],K[8],K[9]]:
    rs=math.sqrt(G*M/a0); d=G*M/(rs*c*c); vf=(G*M*a0)**0.25
    print(f"   {n:15s} δ(r*) = {d:.4e}   (v_flach/c)² = {(vf/c)**2:.4e}   v_flach = {vf/1e3:7.1f} km/s")

print("\n4) Faltungsmesser: Δs = v_flach⁴/(3c²·a0)  (aus GM = 3c²Δs und v⁴ = GM a0)")
for v in [60,110,176,220,300]:
    ds=(v*1e3)**4/(3*c*c*a0)
    print(f"   v_flach = {v:3d} km/s  ->  Δs = {ds:.3e} m = {ds/AU:8.1f} AE   (M = {3*c*c*ds/G/MS:.2e} M☉)")

print("\n5) Scheinbare Zusatzfaltung jenseits r*: Δs_eff(r) = a(r)·r²/(3c²) = Δs·ν,  ν = 1/(1-q), q = exp(-sqrt(a_N/a0))")
M=6e10*MS; ds=G*M/(3*c*c)
for rk in [5,8.35,15.33,30,50,100]:
    r=rk*KPC; aN=G*M/r**2; q=math.exp(-math.sqrt(aN/a0)); nu=1/(1-q)
    print(f"   r = {rk:6.2f} kpc   ν = {nu:6.3f}   Δs_eff = {ds*nu/AU:8.1f} AE   Zusatzstrecke Δs_eff−Δs = {ds*(nu-1)/AU:8.1f} AE")
print("   Tief jenseits r*: Δs_eff ≈ sqrt(GM a0)·r/(3c²) wächst linear mit r (Zusatzstrecke ohne Zusatzmasse).")

print("\n6) Das Ganze auf seiner eigenen Schwelle:")
MU,RU=1.5e53,4.4e26
print(f"   Hüll-Neigung des beobachtbaren Universums g_U = GM/R² = {G*MU/RU**2:.3e} m/s²  (a0 = 1.2e-10, Verhältnis {G*MU/RU**2/a0:.2f})")
gH=c*H0/2
print(f"   Hubble-Kugel bei kritischer Dichte: g_H = G·M_H/R_H² = c·H0/2 = {gH:.3e} m/s²;  g_H/π = c·H0/2π = {gH/math.pi:.3e} m/s²")
print(f"   r* des Universums = sqrt(GM_U/a0) = {math.sqrt(G*MU/a0)/RU:.2f} · R_U")

print("\n7) Jenseits r* (tiefes Regime) ist das Potential logarithmisch: Φ = v_flach²·ln(r):")
vf=(G*6e10*MS*a0)**0.25
for f in [2,10,100,1000]:
    print(f"   Faltungs-Differential zwischen r* und {f:4d}·r*: Δδ = (v/c)²·ln({f}) = {(vf/c)**2*math.log(f):.3e}")
print("   -> Das Differential wächst ohne Grenze: ein kollektiver Trichter hat keinen flachen Außenraum,")
print("      er muss in die Trichter der Nachbarn übergehen (Netzwerk).")
