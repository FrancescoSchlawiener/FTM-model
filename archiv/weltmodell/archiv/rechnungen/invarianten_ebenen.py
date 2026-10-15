# -*- coding: utf-8 -*-
"""
Invarianten nach Ebenen: echte Invarianten des Modells gegen Übereinstimmungen einer Epoche oder eines Beobachters.

Ebenen: L = lokal im Raum, G = global (Torus, Zelle, Horizont von außen), I = innere Richtung (orthogonal zum Raum).
Für jede Zahl: Wert heute, Wert bei z = 0, 1, 10, 1090 (a₀(z) aus der inneren Faltung der Zelle, H(z) ΛCDM-Hintergrund
mit Ω_m aus dem Kern), Wert für andere Beobachter (äußeres Feld 0,03·a₀ und 0,1·a₀ senkrecht/gleichgerichtet,
anderer Ort in der Milchstraße, bewegter Beobachter), Wert für andere Massen.
Klassen: (i) exakt (alle Ebenen, Epochen, Beobachter), (ii) invariant in einer Ebene, (iii) Übereinstimmung einer Epoche,
(iv) beobachterabhängig.
Rechnet nur mit ftm_kern.
"""
import math
import numpy as np
from scipy.integrate import quad
from ftm_kern import *

ZS = [0.0, 1.0, 10.0, 1090.0]
MASSEN = [("Proton", m_p), ("Sonne", M_sun), ("Milchstraße", 6e10 * M_sun), ("Haufen", 1.7e14 * M_sun)]
EPS = 1e-12


def a0_z(z):
    """a₀ bei Rotverschiebung z: Λ als innere Faltung der Zelle, Strahlung früh (Kern: a0_frueh_relativ), normiert auf heute."""
    return A0 * a0_frueh_relativ(z) / a0_frueh_relativ(0.0)


def H_z(z):
    return H0_si() * E_hintergrund(z)


def chi_mitbewegt(z):
    """mitbewegter Abstand 0 → z (heutige Einheiten), in Stücken integriert."""
    f = lambda x: c / H0_si() / E_hintergrund(x)
    if z <= 10:
        return quad(f, 0, z, limit=400)[0]
    return quad(f, 0, 10, limit=400)[0] + quad(f, 10, z, limit=400)[0]


def teilchenhorizont(z):
    """mitbewegter Teilchenhorizont bei z: ∫_z^∞ c dz/H, über ln(1+z)."""
    g = lambda u: c / H0_si() * math.exp(u) / E_hintergrund(math.expm1(u))
    return quad(g, math.log1p(z), 40.0, limit=400)[0]


def bose_integral(d, ge=0.0, lage="senkrecht"):
    """∫(a − a_N)(r*/r)^{d−1} d ln r / a₀ für die Punktmasse, a = ν(|g|)·a_N, |g| mit äußerem Feld ge (in a₀)."""
    def f(s):
        aN = s * s
        if ge == 0:
            g = aN
        elif lage == "senkrecht":
            g = math.hypot(aN, ge)
        else:
            g = aN + ge
        n = 1 / (-math.expm1(-math.sqrt(g)))
        return (n - 1) * aN * s**(d - 2)
    return quad(f, 0, 60, limit=400)[0]


def nu_rstern(ge=0.0, lage="senkrecht"):
    aN = 1.0
    g = aN if ge == 0 else (math.hypot(aN, ge) if lage == "senkrecht" else aN + ge)
    return float(nu(g * A0))


def zeile(name, werte, fmt="{:.4g}"):
    return f"   {name:44s} " + "  ".join(fmt.format(v) for v in werte)


out = []
P = out.append
P("Invarianten nach Ebenen (L lokal / G global / I innere Richtung) – Epoche, Beobachter, Masse")
P(f"Kern: a₀ = {A0:.2e} m/s², H₀ = {H0_KMS} km/s/Mpc, Ω_m = {OMEGA_M:.4f}, Ω_r = {OMEGA_R:.2e}, Q = {Q_KOSM:.4f}, q = {q_KOSM:.4f}")

# ------------------------------------------------------------------ 0 Epoche: Grundgrößen
P("\n0  Hintergrund je Epoche (a₀(z) aus der inneren Faltung der Zelle, Kern a0_frueh_relativ)")
P("   z:                                            " + "  ".join(f"{z:>9g}" for z in ZS))
a0s = [a0_z(z) for z in ZS]
Hs = [H_z(z) for z in ZS]
P(zeile("a₀(z)/a₀(0)", [a / a0s[0] for a in a0s], "{:>9.4f}"))
P(zeile("H(z)/H₀", [h / Hs[0] for h in Hs], "{:>9.4g}"))
Lam = [schliesslaenge(a) for a in a0s]
P(zeile("Λ(z)/Λ(0)", [l / Lam[0] for l in Lam], "{:>9.4f}"))

# ------------------------------------------------------------------ 1 Epochenabhängige kosmische Zahlen
P("\n1  Kosmische Zahlen über die Epoche")
Th = [theta_(a, h) for a, h in zip(a0s, Hs)]
P(zeile("Θ = 2πa₀/(cH)", Th, "{:>9.4g}"))
Oms = [K_S * t**2 for t in Th]
P(zeile("Ω_s = (π²/45)Θ²  [heute π²/45 = 0,2193]", Oms, "{:>9.3g}"))
P(zeile("S_H/S_Λ (= Ω_s)", Oms, "{:>9.3g}"))
ew = [1 - gewicht_w(3 * h, a) for a, h in zip(a0s, Hs)]
P(zeile("1 − w (θ = 3H)", ew, "{:>9.3g}"))
P(zeile("Ω_s/(1 − w)  [36π⁴/45 = 77,93]", [o / e for o, e in zip(Oms, ew)], "{:>9.4f}"))
P(zeile("3H·τ_D  [heute 17,9]", [3 * h * tau_D(a) for a, h in zip(a0s, Hs)], "{:>9.4g}"))
P(zeile("H·τ_D  (Kac-Zeit der Schwelle, Schritt 1/H)", [h * tau_D(a) for a, h in zip(a0s, Hs)], "{:>9.4g}"))
P(zeile("ν_kosm = 1 + Q (Ω_m/Ω_b, alle Epochen)", [1 + Q_KOSM] * 4, "{:>9.4f}"))
P(zeile("(cH/2)/(π·a₀) = 1/Θ  (Randneigung Hubble-Kugel)", [1 / t for t in Th], "{:>9.4g}"))
L0 = zelle_aus_a0(A0)
P(f"   Zelle mitbewegt L = {L0/Mpc/1e3:.2f} Gpc (a₀ = {A0:.2e}); Λ = {Lam[0]/Mpc/1e3:.2f} Gpc; M_U = c⁴/(Ga₀) = {c**4/(G*A0):.3e} kg")
LHc = [L0 / (1 + z) * h / c for z, h in zip(ZS, Hs)]
P(zeile("L·H/c (physikalisch)", LHc, "{:>9.4g}"))
Omz = [OMEGA_M * (1 + z)**3 / E_hintergrund(z)**2 for z in ZS]
P(zeile("(8π²/(Ω_m(z)·Θ))^{1/3} mit a₀(0)", [(8 * math.pi**2 / (om * theta_(A0, h)))**(1 / 3) for om, h in zip(Omz, Hs)], "{:>9.4g}"))
MQ = [(1 + Q_KOSM)**3 * math.log(q_KOSM)**6 * a**3 / (1e4 * G * h**4) / M_sun for a, h in zip(a0s, Hs)]
P(zeile("M_Q [M☉]  (∝ Θ⁴)", MQ, "{:>9.3g}"))
P(zeile("Λ_s·Λ² = π⁴/15 (Vakuumteil der Wirkung)", [lambda_vakuum(l) * l**2 for l in Lam], "{:>9.5f}"))
P(zeile("a₀/c = θ_c [s⁻¹]", [a / c for a in a0s], "{:>9.3g}"))

# Rückstellbilanz im eigenen Horizont und Zelle gegen letzte Streufläche
P("\n   Rückstellbilanz im eigenen Teilchenhorizont und Zelle gegen letzte Streufläche, je Beobachter-Epoche")
rho_b0 = OMEGA_B * 3 * H0_si()**2 / (8 * math.pi * G)
chi1090 = chi_mitbewegt(1090.0)
for zo in [0.0, 0.5, 1.0, 3.0, 10.0, 1090.0]:
    Rp = teilchenhorizont(zo) / (1 + zo)
    rb = rho_b0 * (1 + zo)**3
    dsR = 4 * math.pi * G / (9 * c**2) * rb * Rp**2
    dL = 2 * (chi1090 - chi_mitbewegt(zo))          # mitbewegt, heutige Einheiten
    P(f"   z_Beob = {zo:5.1f}: ΣΔs_b/R (eigener Horizont) = {dsR:.4f}, mit ν_kosm {dsR*(1+Q_KOSM):.3f};  "
      + (f"d_LSS(mitbewegt) = {dL/Mpc/1e3:6.2f} Gpc, L/d_LSS = {L0/dL:.3f}" if dL > 0 else "(Beobachter an der Streufläche)"))
# Zukunft: d_LSS wächst weiter bis zur Grenze des Ereignishorizonts
chi_inf = quad(lambda u: c / H0_si() * math.exp(u) / E_hintergrund(math.expm1(u)), 0, 40)[0]  # = Teilchenhorizont heute
ev = quad(lambda a: c / H0_si() / (a * a * E_hintergrund(1 / a - 1)), 1, 1e4, limit=400)[0]
P(f"   Zukunft: d_LSS wächst bis 2·(χ_1090 + Ereignishorizont) = {2*(chi1090+ev)/Mpc/1e3:.2f} Gpc, L/d_LSS → {L0/(2*(chi1090+ev)):.3f}")

# ------------------------------------------------------------------ 2 Örtliche Zahlen (L, I) und Beobachter
P("\n2  Örtliche Zahlen der Regel: isoliert und mit äußerem Feld g_e (Beobachter in einer Umgebung)")
P("   Größe                                        isoliert   ⊥0,03    ∥0,03    ⊥0,1     ∥0,1  (g_e in a₀)")
cases = [(0.0, "senkrecht"), (0.03, "senkrecht"), (0.03, "gleich"), (0.1, "senkrecht"), (0.1, "gleich")]
P(zeile("ν(r*) = 1/(1 − e⁻¹) = 1,5820", [nu_rstern(g, l) for g, l in cases], "{:>8.4f}"))
for d, ref, nm in [(1, math.pi**2 / 6, "π²/6"), (2, 2 * 1.2020569031595942, "2ζ(3)"), (3, math.pi**4 / 15, "π⁴/15")]:
    P(zeile(f"Leiter d={d}: ∫(a−a_N)(r*/r)^{{d−1}}dlnr/a₀ [{nm}={ref:.4f}]", [bose_integral(d, g, l) for g, l in cases], "{:>8.4f}"))
RQ = 1 / (-math.log(q_KOSM))
P(f"   R_Q/r* = 1/(−ln q) = {RQ:.4f}, a_N(R_Q)/a₀ = (ln q)² = {math.log(q_KOSM)**2:.4f}  (für jede Masse; epochenfest, da Ω_c/Ω_b fest)")
qRQ = [math.exp(-math.sqrt(math.hypot(math.log(q_KOSM)**2, g) if l == "senkrecht" else math.log(q_KOSM)**2 + g)) if g else q_KOSM for g, l in cases]
P(zeile("q_lokal an R_Q (gleich q_kosm nur isoliert)", qRQ, "{:>8.4f}"))

# Glättung: ℓ/r* = (M×/M)^{1/6}; Ort in der Milchstraße (flache Kurve: T ∝ R⁻²)
T_sun = 5.5519e-30          # aus trichter_invarianten.txt (Sonnenort)
P("\n   Glättung ℓ/r* = (M×/M)^{1/6}, M× = a₀³/(GT²): Sonne an anderen Orten der Milchstraße (flache Kurve, T ∝ R⁻², R☉ = 8,2 kpc)")
for Rk in [4.0, 8.2, 16.0, 50.0]:
    T = T_sun * (8.2 / Rk)**2
    Mx = A0**3 / (G * T**2)
    ell = (G * M_sun / T)**(1 / 3)
    P(f"   R = {Rk:5.1f} kpc: M× = {Mx/M_sun:.3e} M☉, ℓ/r* = {(Mx/M_sun)**(1/6):.2f} (direkt {ell/r_stern(M_sun):.2f})")

# Beobachter im Potential und bewegt
MW = 6e10 * M_sun
r8 = 8.2e3 * pc
s8 = r_stern(MW) / r8
dN = G * MW / (c**2 * r8)
dX = (v_flach(MW) / c)**2 * math.log(float(nu((G * MW / r8**2))))
P(f"\n   Beobachter am Sonnenort: Potentialtiefe (Newton + Regel) Φ/c² ≈ {dN:.2e} + {dX:.2e}: jede Taktlesung ändert sich um ~1e-6")
beta = 369.8e3 / c   # [G] CMB-Dipol
P(f"   bewegter Beobachter β = {beta:.2e} [G]: Regel ist Galilei-invariant (a_N, ν, q unverändert); H aus Rotverschiebungen "
  f"richtungsabhängig um β/z (z = 0,05: {beta/0.05*100:.1f} %), Monopol um β² = {beta**2:.1e}")

# ------------------------------------------------------------------ 3 Massen
P("\n3  Massenabhängigkeit (a₀ fest)")
P("   Masse          Δs[m]       r*[m]       r*²λ̄/(ℓ_P²ℓ₀)  3Δsλ̄/ℓ_P²  δ*=(v_f/c)²  r*/r₀·2δ*  α∞/(2πδ*)  ν(r*)")
l0 = c**2 / A0
for nm, M in MASSEN:
    lb = compton_red(M)
    rs = r_stern(M)
    dstar = (v_flach(M) / c)**2
    P(f"   {nm:12s} {ds(M):10.3e} {rs:10.3e}  {rs**2*lb/(l_P**2*l0):14.10f}  {3*ds(M)*lb/l_P**2:10.8f}  {dstar:10.3e}  "
      f"{rs/r0(M)*2*dstar:9.6f}  {ablenkung_kollektiv(M)/(2*math.pi*dstar):9.6f}  {nu_rstern():.5f}")
MU = c**4 / (G * A0)
P(f"   Δs(M_U)/Λ = {ds(MU)/schliesslaenge():.6f} (⅔);  M_Zelle = 3c²Λ/G = {3*c**2*schliesslaenge()/G:.3e} kg = {3*c**2*schliesslaenge()/G/MU:.4f}·M_U")
epsL = A0 / (c**2 / l_P)
Mmin, Mc = epsL * m_Pl, epsL**(1 / 3) * m_Pl
P(f"   ε_L = a₀/a_P = {epsL:.3e}; M_min = {Mmin:.3e} kg, M_c = {Mc:.3e} kg, M_U = {MU:.3e} kg")
P(f"   M_min·M_U/m_Pl² = {Mmin*MU/m_Pl**2:.8f};  M_c³·M_U/m_Pl⁴ = {Mc**3*MU/m_Pl**4:.8f};  Λ_s·ℓ_P² = {lambda_vakuum(schliesslaenge())*l_P**2:.3e} "
  f"= 4(π⁴/15)ε_L² = {4*math.pi**4/15*epsL**2:.3e}")
P("   ε_L über die Epoche: " + ", ".join(f"z={z:g}: {a/(c**2/l_P):.3e}" for z, a in zip(ZS, a0s)))
P(f"   Kac ⟨T⟩ = 2Λ/c = τ_D = {2*schliesslaenge()/c/Gyr:.1f} Mrd. Jahre, für jede Masse (Schrittweite r₀ = 6Δs kürzt sich: r*² = r₀·Λ): "
  + ", ".join(f"{nm} {r_stern(M)**2/(r0(M)*schliesslaenge()):.6f}" for nm, M in MASSEN))
P(f"   Spin ½: a_J·Δs/ℓ_P² (a_J = ħ/(2Mc)) = " + ", ".join(f"{nm} {hbar/(2*M*c)*ds(M)/l_P**2:.6f}" for nm, M in MASSEN[:2]) + "  (1/6 = 0,166667)")
P("   drei Raten ω_C²/(Γ_G·c/3Δs) = " + ", ".join(f"{nm} {(c/compton_red(M))**2/(alpha_G(M)*c/compton_red(M)*c/(3*ds(M))):.6f}" for nm, M in MASSEN))

# ------------------------------------------------------------------ 4 Konstanten-Identitäten
P("\n4  Identitäten der Konstanten (L), sternunabhängig")
Dstar = ds(sternmasse())
C_konst = alpha_G() * (m_p / m_e)**4 / alpha**12
f_s = kB * TEFF_sun / (alpha**2 * m_e * c**2)
dsun = 3 * ds(M_sun) / R_sun
LE = eddington(M_sun)
sqrtC_stern = (2 * math.pi**3 / 45) * f_s**4 * (ds(M_sun) / Dstar) * (alpha**2 / dsun)**2 * LE / L_sun
P(f"   𝒞 = α_G(m_p/m_e)⁴/α¹² = {C_konst:.4f};  aus den Sternzahlen der Sonne (√𝒞)² = {sqrtC_stern**2:.4f}")
sig_SB = math.pi**2 * kB**4 / (60 * hbar**3 * c**2)
for nm, Mf, Rf, Tf in [("Stern 0,5 M☉/0,5 R☉/3800 K", 0.5, 0.5, 3800.0), ("Stern 10 M☉/5 R☉/25000 K", 10.0, 5.0, 25000.0)]:
    M, R = Mf * M_sun, Rf * R_sun
    L = 4 * math.pi * R**2 * sig_SB * Tf**4
    fs = kB * Tf / (alpha**2 * m_e * c**2)
    v = (2 * math.pi**3 / 45) * fs**4 * (ds(M) / Dstar) * (alpha**2 / (3 * ds(M) / R))**2 * eddington(M) / L
    P(f"   {nm}: 𝒞 = {v**2:.4f}  (M, R, T kürzen sich: √𝒞 = ħm_p/(3α⁶m_e²cΔs_*))")
P(f"   Sonne: L_Stefan/L☉ = {4*math.pi*R_sun**2*sig_SB*TEFF_sun**4/L_sun:.4f}")
# Spaltbarkeit: B/A = a_v − a_s A^{-1/3} − a_c k²A^{2/3}, Maximum in A
for a_s, a_c in [(17.8, 0.711), (18.3, 0.714), (16.0, 0.8)]:
    k = 0.46
    A = np.linspace(10, 300, 290001)
    BA = -a_s * A**(-1 / 3) - a_c * k * k * A**(2 / 3)
    Am = A[np.argmax(BA)]
    ES, EC = a_s * Am**(2 / 3), a_c * k * k * Am**(5 / 3)
    P(f"   Tröpfchen a_s = {a_s}, a_c = {a_c}, Z/A = {k}: Maximum B/A bei A = {Am:.1f}, x = E_C/(2E_S) = {EC/(2*ES):.5f} (¼ für alle Koeffizienten)")

# ------------------------------------------------------------------ 5 Formgleichheit über die Ebenen
P("\n5  Formgleichheit über die Ebenen (Wert heute, z = 0)")
Hn = Hs[0]
rH = c / Hn
DsH = rH / 6
k_L = c**2 / (2 * r0(M_sun))
k_G = c**2 / rH
P(f"   Oberflächenschwere κ = c²/(2ℓ): L Sonne ℓ = r₀ = 6Δs: {k_L:.3e};  G ℓ = r_H/2 = 3Δs_H: {k_G:.3e} m/s²;  I: a₀ = c²/(2Λ) = {A0:.3e} m/s²")
T_L = hbar * k_L / (2 * math.pi * c * kB)
T_G = hbar * Hn / (2 * math.pi * kB)
T_I = hbar * c / (2 * Lam[0] * kB)
P(f"   Temperatur T = ħ/(k·Periode): L Periode 2πc/κ_L: {T_L:.3e} K;  G Periode 2π/H: {T_G:.3e} K;  I Periode 2Λ/c: {T_I:.3e} K")
P(f"   a₀ = c·kT_I/ħ = {c*kB*T_I/hbar:.4e} (exakt, jede Epoche);  c·kT_GH/ħ = cH/2π = {c*kB*T_G/hbar:.4e};  T_I/T_GH = Θ = {T_I/T_G:.4f}")
P("   T_I/T_GH über z: " + ", ".join(f"z={z:g}: {math.pi*c/(l*h):.4g}" for z, h, l in zip(ZS, Hs, Lam)))
S_L = 36 * math.pi * (ds(M_sun) / l_P)**2
S_G = 36 * math.pi * (DsH / l_P)**2
S_I = 45 / math.pi**3 * (Lam[0] / l_P)**2
P(f"   Entropie Länge²/ℓ_P²: L 36π(Δs/ℓ_P)² Sonne {S_L:.3e};  G 36π(Δs_H/ℓ_P)² = πr_H²/ℓ_P² {S_G:.3e} (Kern {horizontentropie(Hn):.3e});  "
  f"I (45/π³)(Λ/ℓ_P)² {S_I:.3e};  S_G/S_I = {S_G/S_I:.4f} = Ω_s {Oms[0]:.4f}")
P("   S_G/S_I über z: " + ", ".join(f"z={z:g}: {36*math.pi*(c/h/6/l_P)**2/(45/math.pi**3*(l/l_P)**2):.3g}" for z, h, l in zip(ZS, Hs, Lam)))
# Kac-Zählung
P(f"   Kac-Zählung: L ν = 1/(1−q) an r* {nu_rstern():.4f}, an R_Q {1/(1-q_KOSM):.4f};  G ν_kosm = Ω_m/Ω_b = {OMEGA_M/OMEGA_B:.4f};  "
  f"I ⟨T⟩·H = 2ΛH/c = H·τ_D = {Hn*tau_D():.4f} (heute; z=1: {Hs[1]*tau_D(a0s[1]):.3f})")
# geometrisches Mittel mit Entartung
P(f"   geometrisches Mittel, Entartung aller drei Längen: L bei m_Pl: λ̄ = {compton_red(m_Pl):.4e}, 3Δs = {3*ds(m_Pl):.4e}, ℓ_P = {l_P:.4e} m;  "
  f"G/I bei M_U: r* = {r_stern(MU):.4e}, 3Δs = {3*ds(MU):.4e}, ℓ₀ = 2Λ = {2*schliesslaenge():.4e} m")
P(f"   Hubble-Masse: 6Δs(M_H) = {6*ds(4*math.pi/3*rH**3*3*Hn**2/(8*math.pi*G))/rH:.6f}·r_H (Grenzradius = Horizont, jede Epoche);  "
  f"r*(M_H)/r_H = √(π/Θ) = {r_stern(4*math.pi/3*rH**3*3*Hn**2/(8*math.pi*G))/rH:.4f}")
# Bose-Form auf G
sk = -math.log(q_KOSM)
P(f"   Bose-Form auf G: q_kosm = e^(−s_k), s_k = {sk:.4f}, ν_kosm = 1 + n_B(s_k) = {1+float(n_B(sk)):.4f};  χ(s_k) = {sk*float(n_B(sk)):.4f}")
P(f"   Zählung je Radiant: 2π − 1 = {2*math.pi-1:.4f} gegen Q = {Q_KOSM:.4f};  heute H₀τ_D − 1 = {Hn*tau_D()-1:.4f}")

text = "\n".join(out)
print(text)
with open(__file__.replace(".py", ".txt"), "w", encoding="utf-8") as fh:
    fh.write(text + "\n")
