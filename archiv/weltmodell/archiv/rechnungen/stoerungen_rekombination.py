"""Wirkt die örtliche Regel a = ν(|g_N|)·g_N auf lineare Störungen bei der Rekombination?

1 a_N/a₀(z) für Φ/c² = 10⁻⁵ und für das gemessene Spektrum (CAMB, Weyl-Potential) bei z = 1090, 100, 10, 0
2 Schalter: Gewicht w(θ) mit θ = 3H bzw. θ_lok = 3H − δ̇, Γ/(3H), ρ̄/ρ_J einer linearen Störung
3 Gegenprobe gedacht: wäre ν(a_N) bei z = 1090 wirksam, Größe von ν − 1 auf den Skalen der Peaks
4 Einsatz der Regel: Umkehrpunkt je Skala (σ(R,z)·ν_pk = 1,062), Übergangsbreite w 0,9 → 0,1
Blickrichtung: L lokal (am Ort), G global (atmendes Ganzes), I innere Richtung (ϖ, Λ, τ_D).
Alle Konstanten aus ftm_kern; Spektrum [G]: A_s = 2,1·10⁻⁹, n_s = 0,9649 (Planck 2018), τ = 0,0544.
"""
import math, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K

OUT = []
def P(s=""):
    print(s); OUT.append(s)

c, G, Mpc = K.c, K.G, K.Mpc
H0 = K.H0_si()
def H(z): return H0 * K.E_hintergrund(z)
def a0z(z, a0=K.A0): return a0 * K.a0_frueh_relativ(z)     # Strahlung erhöht Λ früh; seit z ≈ 10 fest
def Om_z(z):
    return K.OMEGA_M * (1 + z)**3 / K.E_hintergrund(z)**2

ZS = (1090.0, 100.0, 10.0, 0.0)
LAM = (1, 3, 10, 30, 100, 147, 300, 500)                 # mitbewegte Skala in Mpc

P("=" * 100)
P("STÖRUNGEN BEI DER REKOMBINATION: wirkt die örtliche Regel auf lineare Störungen?")
P("=" * 100)
P(f"a₀ = {K.A0:.2e} (Modell), Literatur {K.A0_LIT:.2e};  a₀(z)/a₀ = 1/(1 + Ω_r(1+z)/Ω_m):")
for z in ZS:
    P(f"   z = {z:6.0f}:  a₀(z)/a₀ = {K.a0_frueh_relativ(z):.4f},  H = {H(z):.3e} s⁻¹,  3Hc = {3*H(z)*c:.3e} m/s²  (Blick G)")

# ---------------------------------------------------------------- 1a feste Amplitude
P("\n1a  Φ/c² = 10⁻⁵ fest, a_N = k_phys·Φ   (Blick L; a₀(z) Blick I)")
P("    Spalten: a_N/a₀(z) mit k = 1/λ (Lesart a_N ≈ Φ/λ_phys) | mit k = 2π/λ | ν(k=1/λ)")
Phi = 1e-5 * c**2
hdr = "   λ [Mpc] " + "".join(f"| z = {z:<6.0f}                    " for z in ZS)
P(hdr)
for lam in LAM:
    row = f"   {lam:7.0f}  "
    for z in ZS:
        lp = lam * Mpc / (1 + z)
        y1 = Phi / lp / a0z(z); y2 = 2 * math.pi * y1
        row += f"| {y1:9.3g} {y2:9.3g} ν={float(K.nu(y1, 1.0)):6.3f} "
    P(row)
lp = 147 * Mpc / 1091
P(f"   Nachrechnung 147 Mpc, z = 1090: a_N = Φ/λ_phys = {Phi/lp:.4e} m/s²;  /(0,76·1,2·10⁻¹⁰) = {Phi/lp/a0z(1090, K.A0_LIT):.3f};"
  f"  /(0,76·1,1·10⁻¹⁰) = {Phi/lp/a0z(1090):.3f};  mit 2π/λ: {2*math.pi*Phi/lp/a0z(1090):.2f}")

# ---------------------------------------------------------------- 1b gemessenes Spektrum
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=K.H0_KMS, ombh2=K.OMEGA_B_H2, omch2=K.OMEGA_C_H2, mnu=0.06, tau=0.0544)
p.InitPower.set_params(As=2.1e-9, ns=0.9649)
p.set_matter_power(redshifts=[0.0], kmax=200)
PW = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False,
                                        kmax=200, var1='Weyl', var2='Weyl', zmax=1100)
PD = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=False, k_hunit=False,
                                        kmax=200, var1='delta_tot', var2='delta_tot', zmax=1100)
def DPhi(z, k):   # Δ_Φ pro ln k (Weyl-Potential (Φ+Ψ)/2 in Einheiten c²), k in 1/Mpc
    return math.sqrt(k**3 * PW.P(z, k) / (2 * math.pi**2)) / k**2
def Ddel(z, k):
    return math.sqrt(k**3 * PD.P(z, k) / (2 * math.pi**2))

P("\n1b  Gemessenes Spektrum (CAMB): Δ_Φ(k) pro ln k, a_N = k_phys·Δ_Φ·c²   (Blick L, gelesen aus dem Feld des Ganzen G)")
P("   λ = 2π/k [Mpc] | " + " | ".join(f"z={z:<5.0f} Δ_Φ   a_N/a₀(z)" for z in ZS))
for lam in LAM:
    k = 2 * math.pi / lam
    row = f"   {lam:12.0f}   "
    for z in ZS:
        d = DPhi(z, k); aN = k * (1 + z) / Mpc * d * c**2
        row += f"| {d:8.2e} {aN/a0z(z):9.3g}   "
    P(row)
lk = np.linspace(math.log(1e-4), math.log(100.0), 3000)
P("   Effektivwert des Gradienten über alle k (bis 100/Mpc), √∫k_phys²Δ_Φ² dln k:")
AN_RMS = {}
for z in ZS:
    ks = np.exp(lk)
    integ = np.array([(k * (1 + z) / Mpc * DPhi(z, k) * c**2)**2 for k in ks])
    AN_RMS[z] = math.sqrt(np.trapezoid(integ, lk))
    kpk = ks[np.argmax(integ)]
    P(f"      z = {z:6.0f}:  a_N,rms = {AN_RMS[z]:.3e} m/s² = {AN_RMS[z]/a0z(z):.3g}·a₀(z);  Maximum pro ln k bei λ = {2*math.pi/kpk:.0f} Mpc")
k147 = 2 * math.pi / 147
P(f"   Nachrechnung: Δ_Φ(k = 2π/147 Mpc, z = 1090) = {DPhi(1090, k147):.2e}; a_N = k·(1+z)·Δ_Φ·c²/Mpc = "
  f"{k147*1091/Mpc*DPhi(1090,k147)*c**2:.3e} m/s²")

# ---------------------------------------------------------------- 2 Schalter
P("\n2  Schalter der Regel. w = θ̃²/(θ̃² + θ_c²), θ_c = a₀(z)/c (I), θ = Ausdehnungsrate des Bereichs (G am Ort L)")
P("   Lineare Störung: θ_lok = 3H − δ̇ = 3H(1 − f·δ/3), f = d ln D/d ln a; Γ = a_N/c (L); ρ̄/ρ_J = 1/(3δ) (ebene Welle, ρ_J = 3√(3/2)|E^TF|/4πG)")
def growth_f(z, dz=1e-3):
    k = 0.05
    a1, a2 = 1 / (1 + z + dz), 1 / (1 + max(z - dz, 0))
    d1, d2 = Ddel(z + dz, k), Ddel(max(z - dz, 0), k)
    return math.log(d2 / d1) / math.log(a2 / a1)
for z in ZS:
    th = 3 * H(z); thc = a0z(z) / c
    w = K.gewicht_w(th, a0z(z))
    f = growth_f(z)
    P(f"   z = {z:6.0f}:  3H·τ_D = {th/thc:9.3e},  1 − w(3H) = {(thc/th)**2/(1+(thc/th)**2):.2e},  f = {f:.3f}")
    for lam in (10, 147, 500):
        k = 2 * math.pi / lam
        d = Ddel(z, k); aN = k * (1 + z) / Mpc * DPhi(z, k) * c**2
        thl = th * (1 - f * d / 3)
        P(f"        λ = {lam:4d} Mpc: δ = {d:.2e}, θ_lok/3H = {thl/th:.6f}, 1 − w(θ_lok) = {1-K.gewicht_w(thl, a0z(z)):.2e},"
          f"  Γ/3H = {aN/c/th:.2e},  ρ̄/ρ_J = {1/(3*d):.2e},  ℓ/R = {(1/(3*d))**(1/3):.2f}")
    aN = AN_RMS[z]
    P(f"        Effektivwert: Γ/3H = {aN/c/th:.2e}")
g_mw = K.A0
P(f"   Zum Vergleich gebundene Galaxie (Regel wirkt), a_N = a₀ bei r*: Γ/3H₀ = {g_mw/c/(3*H0):.3f} < 1, θ_lok ≤ 0 ⇒ w = 0.")
P("   ⇒ Γ/3H ist nicht der Schalter (in Galaxien ebenfalls < 1); a_N/a₀ auch nicht (heute ≪ 1 und trotzdem Staub);")
P("     entscheidend ist θ_lok gegen θ_c: Die Regel greift erst, wenn der Bereich nicht mehr mitatmet (θ̃ ≲ θ_c).")
P("   Umrechnung: Γ/θ = (a_N/a₀)·(θ_c/θ) = s²·√((1 − w)/w).")
P("   ρ̄/ρ_J = 1/(3δ) ≫ 1 und Glättungslänge ℓ/R = (ρ̄/ρ_J)^{1/3} > 1: formal läge die Störung ganz in ihrer eigenen")
P("   Newton-Zone (Einzeltrichter); der kosmologische Grenzfall der Körperdefinition ist aber offen, entscheidend bleibt w.")
P("   (z = 0, λ = 10 Mpc: δ > 1, nichtlinear, die linearen Werte gelten dort nicht.)")

# ---------------------------------------------------------------- 3 Gegenprobe gedacht
P("\n3  Gedachte Gegenprobe: ν(a_N) bei z = 1090 wirksam (Spektrum-Niveau, k = 2π/λ), Skalen der Peaks")
dA = K.d_LSS(1090.0) / 2
P(f"   Abstand zur letzten Streufläche (mitbewegt) {dA/Mpc:.0f} Mpc; ℓ ≈ k·d")
for l in (30, 220, 540, 810, 1130, 1450):
    k = l / (dA / Mpc)
    aN = k * 1091 / Mpc * DPhi(1090, k) * c**2
    n = float(K.nu(aN, a0z(1090)))
    P(f"   ℓ = {l:5d}: k = {k:.4f}/Mpc, λ = {2*math.pi/k:6.0f} Mpc, a_N/a₀(z) = {aN/a0z(1090):6.3f}, ν = {n:.3f}, ΔΦ/Φ = {n-1:+.1%}")
P("   Planck [G]: Peakhöhen auf ~1 %, 100θ* = 1,04109 ± 0,00030, r_d = 147,09 ± 0,26 Mpc.")
P("   Abschätzung: +3–4 % im Potential auf den Peaks (Antrieb, C_ℓ etwa doppelt: 6–8 %) und +26 % bei ℓ ≈ 30 (Sachs-Wolfe,")
P("   C_ℓ ~ +50 %) lägen weit über der Planck-Genauigkeit; r_d (Schallgeschwindigkeit, Hintergrund) bliebe, die BAO-Amplitude nicht.")

# ---------------------------------------------------------------- 4 Einsatz der Regel
P("\n4  Einsatz: Ein Bereich verlässt das Atmen am Umkehrpunkt (δ_lin = 1,062); dort fällt w von 0,9 auf 0,1.")
def sigR(z, R):   # Top-Hat, Lagrange-Radius R in Mpc
    ks = np.exp(np.linspace(math.log(1e-4), math.log(200.0), 4000))
    x = ks * R; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    Pk = PD.P(z, ks)
    return math.sqrt(np.trapezoid(ks**3 * Pk / (2 * math.pi**2) * W**2, np.log(ks)))
rho_m0_Msun_Mpc3 = K.rho_m_mitbewegt() * Mpc**3 / K.M_sun
s8 = sigR(0.0, 8 / K.h_)
P(f"   σ₈ (CAMB, A_s [G]) = {s8:.3f}")
zgrid = np.concatenate([np.linspace(0, 10, 201), np.geomspace(10.5, 200, 200)])
def D_rel(z): return Ddel(z, 0.05) / Ddel(0.0, 0.05)
Drel = np.array([D_rel(z) for z in zgrid])
P("   R_L [Mpc]   M [M☉]     σ(R,0)   z des Umkehrpunkts für ν_pk = 1 | 2 | 3")
for R in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 73.5, 147.0):
    s0 = sigR(0.0, R); M = 4 * math.pi / 3 * rho_m0_Msun_Mpc3 * R**3
    zt = []
    for nupk in (1, 2, 3):
        need = 1.062 / (nupk * s0)
        if need > 1: zt.append("  nie  ")
        else: zt.append(f"{np.interp(need, Drel[::-1], zgrid[::-1]):6.2f} ")
    P(f"   {R:8.1f}   {M:9.2e}   {s0:7.4f}   " + " | ".join(zt))

P("   Übergangsbreite (w 0,9 → 0,1) für eine Kugel, die bei z_c kollabiert (Zykloide, θ = 3Ṙ/R, t_c = Alter bei z_c):")
zz = np.geomspace(1e-4, 3000, 40001)
def alter(z):
    zs = zz[zz >= z]; Ez = np.sqrt(K.OMEGA_M*(1+zs)**3 + K.OMEGA_R*(1+zs)**4 + K.OMEGA_DE)
    return np.trapezoid(1 / ((1 + zs) * Ez), zs) / H0
def breite(zc):
    tc = alter(zc); B = tc / (2 * math.pi)
    eta = np.linspace(1e-3, 2 * math.pi - 1e-3, 400001)
    th = 3 * np.sin(eta) / (B * (1 - np.cos(eta))**2)
    t = B * (eta - np.sin(eta))
    ta = np.interp(math.pi * B, t, t)
    zta = np.interp(-ta, -np.array([alter(x) for x in np.geomspace(zc, zc*3+5, 400)]), np.geomspace(zc, zc*3+5, 400))
    thc = a0z(zta) / c
    i9 = np.argmax(th < 3 * thc); i1 = np.argmax(th < thc / 3)
    return (t[i1] - t[i9]) / K.yr / 1e6, zta
for zc in (30, 20, 10, 5, 3.4, 2, 1, 0.5):
    b, zta = breite(zc)
    P(f"      z_c = {zc:5.1f} (Umkehr z = {zta:5.2f}):  Breite {b:8.2f} Mio. Jahre")
P("   (Zykloide ohne Rückstellung im Inneren; bei z_c ≤ 1 überschätzt sie die Breite leicht.)")

# ---------------------------------------------------------------- 5 Handprüfungen
P("\n5  Von Hand (unabhängig vom Code oben):")
P("   a_N(147 Mpc, z = 1090) = 10⁻⁵·(2,998·10⁸)²/(147·3,0857·10²² m/1091) = 8,988·10¹¹/4,158·10²¹ = 2,162·10⁻¹⁰ m/s²;")
P("     /(0,759·1,2·10⁻¹⁰ = 9,11·10⁻¹¹) = 2,37;  /(0,759·1,1·10⁻¹⁰ = 8,35·10⁻¹¹) = 2,59")
P("   H(1090) = 2,184·10⁻¹⁸·√(0,3134·1,2986·10⁹ + 9,14·10⁻⁵·1,4168·10¹²) = 2,184·10⁻¹⁸·2,316·10⁴ = 5,06·10⁻¹⁴ s⁻¹;")
P("     θ_c/3H = (8,35·10⁻¹¹/2,998·10⁸)/(1,518·10⁻¹³) = 2,79·10⁻¹⁹/1,518·10⁻¹³ = 1,84·10⁻⁶ ⇒ 1 − w = 3,4·10⁻¹²")
P("   Γ/3H (147 Mpc, Φ/c² = 10⁻⁵, k = 1/λ) = 2,162·10⁻¹⁰/(2,998·10⁸·1,518·10⁻¹³ = 4,55·10⁻⁵) = 4,8·10⁻⁶")

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
