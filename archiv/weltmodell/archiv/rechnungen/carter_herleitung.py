"""Carters Beziehung α_G ≈ α¹²(m_e/m_p)⁴ aus dem Sterntrichter.

Teil 1  Größen des Trichters in Faltungs- und EM-Längen
Teil 2  Die Kette Hydrostatik -> Zündung -> Transport -> Photosphäre als Monom-Rechnung
        (Exponenten 12 und 4, Herkunft jedes Exponenten, Thomson gegen Kramers)
Teil 3  Einfachste Form: optische Tiefe des Trichters = (T_c/T_ion)⁴
Teil 4  Polytrope n = 3 (Eddington-Standardmodell) für die Strukturzahlen
Teil 5  Vorfaktor 2,94 mit Sonnenzahlen (exakte Identität + Zerlegung)
Teil 6  Grenzmasse, Hayashi-Grenze, Massenfenster, Breite der Bedingung
Teil 7  Wasserstoff-Ionisation und adiabatischer Gradient (Konvektionsschwelle)
Teil 8  Torus-Dynamik: Faltungsrate gegen EM-Raten des Wasserstoffs
"""
import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

# ---------------------------------------------------------------- Konstanten (SI, CODATA 2018)
G, c, hbar, kB, eV = 6.67430e-11, 299792458.0, 1.054571817e-34, 1.380649e-23, 1.602176634e-19
mp, me, alpha = 1.67262192369e-27, 9.1093837015e-31, 7.2973525693e-3
sigSB = 5.670374419e-8
a_rad = 4 * sigSB / c
muB, mu0, g_p = 9.2740100783e-24, 1.25663706212e-6, 5.5856946893
# Sonne (IAU-Nennwerte) und Standard-Sonnenmodell (Bahcall, Serenelli, Basu 2005)
MS, RS, LS, TeffS = 1.98847e30, 6.957e8, 3.828e26, 5772.0
TcS, rhocS = 1.57e7, 1.50e5          # K, kg/m³

lp, le = hbar / (mp * c), hbar / (me * c)
r_e = alpha * le
ds_p = G * mp / (3 * c**2)
aG = G * mp**2 / (hbar * c)
mPl = math.sqrt(hbar * c / G)
Mstar = mPl**3 / mp**2
sigT = 8 * math.pi / 3 * r_e**2
CARTER = aG / (alpha**12 * (me / mp)**4)
Ry = 0.5 * alpha**2 * me * c**2
wp, we = mp * c**2 / hbar, me * c**2 / hbar
wa = alpha**2 * we                    # Umlauffrequenz des Elektrons im Wasserstoff, αc/a_B
GammaG = aG * wp

def h(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

# =============================================================================================
h("Teil 1  Größen des Trichters in Faltungs- und EM-Längen")
ds = lambda M: G * M / (3 * c**2)
print(f"α_G = 3Δs_p/λ̄_p = {3*ds_p/lp:.6e}     α = r_e/λ̄_e = {r_e/le:.10f}")
print(f"Carter-Verhältnis 𝒞 = α_G/(α¹²(m_e/m_p)⁴) = {CARTER:.4f}")
print(f"Längenform 3Δs_p·λ̄_e¹⁶/(r_e¹²·λ̄_p⁵) = {3*ds_p*le**16/(r_e**12*lp**5):.4f}")
print(f"Vierte Wurzel: (3Δs_p/λ̄_p)^(1/4) = {(3*ds_p/lp)**0.25:.5e}   r_e³·λ̄_p/λ̄_e⁴ = {r_e**3*lp/le**4:.5e}"
      f"   Verhältnis {(3*ds_p/lp)**0.25/(r_e**3*lp/le**4):.4f} = 𝒞^(1/4)")
print(f"\nSternmasse M_* = m_Pl³/m_p² = {Mstar:.4e} kg = {Mstar/MS:.4f} M☉;  Δs_* = {ds(Mstar):.1f} m;"
      f"  3Δs_*/λ̄_p = {3*ds(Mstar)/lp:.4e} = α_G^(-1/2) = {aG**-0.5:.4e}")
print(f"Teilchenzahl N_* = M_*/m_p = Δs_*/Δs_p = {Mstar/mp:.4e} = α_G^(-3/2) = {aG**-1.5:.4e}")
LEdd = lambda M: 4 * math.pi * c * G * M * mp / sigT
print(f"Eddington-Leuchtkraft L_Edd = 4πcGMm_p/σ_T = (9/2)·m_p c³·Δs/r_e²:"
      f"  Sonne {LEdd(MS):.4e} W  |  {4.5*mp*c**3*ds(MS)/r_e**2:.4e} W   (L☉/L_Edd = {LS/LEdd(MS):.3e})")
print(f"Ionisationsskala α²m_ec² = ħc·r_e²/λ̄_e³: {alpha**2*me*c**2/eV:.4f} eV | {hbar*c*r_e**2/le**3/eV:.4f} eV;"
      f"  Bohr-Radius a_B = λ̄_e²/r_e = {le**2/r_e:.4e} m")
EG = math.pi**2 * alpha**2 * mp * c**2
print(f"Gamow-Energie p+p: E_G = π²α²m_pc² = π²·ħc·r_e²/(λ̄_e²λ̄_p) = {EG/eV/1e3:.2f} keV"
      f" | {math.pi**2*hbar*c*r_e**2/(le**2*lp)/eV/1e3:.2f} keV")
eps_S = 3 * ds(MS) / RS
print(f"Faltungsdifferential an der Sonnenoberfläche ε = 3Δs/R = {eps_S:.4e}  (Rotverschiebung cε = {c*eps_S:.0f} m/s)")

# =============================================================================================
h("Teil 2  Die Kette als Monom-Rechnung (natürliche Einheiten ħ = c = k = 1)")
BASIS = ['a', 'me', 'mp', 'mPl', 'mu', 'f', 'fs', 'th', 'k', 'mumol', 'CE', 'CK']
NAME = {'a': 'α', 'me': 'm_e', 'mp': 'm_p', 'mPl': 'm_Pl', 'mu': 'μ_M', 'f': 'f_nuc', 'fs': 'f_s',
        'th': 'θ', 'k': 'κ̄/κ_T', 'mumol': 'μ_mol', 'CE': 'C_E', 'CK': 'C_K·D'}
IDX = {b: i for i, b in enumerate(BASIS)}


class Mono:
    def __init__(self, num=1.0, **e):
        self.num = float(num)
        self.e = np.zeros(len(BASIS))
        for k_, v in e.items():
            self.e[IDX[k_]] = v

    def __mul__(self, o):
        o = o if isinstance(o, Mono) else Mono(o)
        r = Mono(self.num * o.num); r.e = self.e + o.e; return r
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, Mono) else Mono(o)
        r = Mono(self.num / o.num); r.e = self.e - o.e; return r

    def __rtruediv__(self, o):
        return Mono(o) / self

    def __pow__(self, p):
        r = Mono(self.num**p); r.e = self.e * p; return r

    def show(self):
        s = f"{self.num:.5g}"
        for b in BASIS:
            v = self.e[IDX[b]]
            if abs(v) > 1e-12:
                s += f"·{NAME[b]}^{v:g}"
        return s


def als_alphaG(cond):
    """cond = 1  ->  α_G = m_p²/m_Pl² = Vorfaktor · Π x^e."""
    v = cond.e
    vP = v[IDX['mPl']]
    s = 2.0 / vP
    ex = v * s
    ex[IDX['mPl']] = 0.0
    ex[IDX['mp']] += 2.0
    r = Mono(cond.num**s); r.e = ex
    return r, s


pi = math.pi
Gm = Mono(mPl=-2)                                   # G = 1/m_Pl²
M = Mono(mu=1, mPl=3, mp=-2)                        # M = μ_M · M_*
Tc = Mono(f=1, a=2, mp=1)                           # Zündung: kT_c = f_nuc·α²m_pc²
GMm = Gm * M * Mono(mp=1)                           # G M m_p  (= μ_M · m_Pl/m_p)
R = Mono(th=1) * GMm / Tc                           # Hydrostatik: kT_c = θ·G M m_p/R
sig = Mono(8 * pi / 3, a=2, me=-2)                  # σ_T = (8π/3) r_e²,  r_e = α/m_e
L_Edd = Mono(4 * pi) * GMm / sig
one_m_beta = Mono(mu=2, mumol=4, CE=-2)             # Eddington: 1-β = (μ_M μ_mol²/C_E)²  (1-β ≪ 1)
lam_T = one_m_beta / Mono(k=1)                      # λ = L/L_Edd = (κ_T/κ̄)(1-β)
sSB = Mono(pi**2 / 60)
Teff4 = lam_T * L_Edd / (Mono(4 * pi) * R**2 * sSB)
target = Mono(fs=1, a=2, me=1)**4                   # Photosphäre an der H-Schwelle: kT_eff = f_s·α²m_ec²
cond = Teff4 / target
print("T_eff⁴ =", Teff4.show())
res, s = als_alphaG(cond)
print("Bedingung T_eff = f_s α² m_e  =>  α_G =", res.show())
print(f"   Vorfaktor = (2π³/45)² = {(2*pi**3/45)**2:.5f}  (Rechnung: {res.num:.5f})")
print(f"   -> Exponent von α: {res.e[IDX['a']]:g},  von m_e: {res.e[IDX['me']]:g},  von m_p: {res.e[IDX['mp']]:g}")

print("\nHerkunft der Exponenten (Beitrag jedes Kettenglieds zum Exponenten in α_G):")
factors = [
    ("Photosphäre  (kT_eff = f_s α² m_e)⁻⁴", 1 / target),
    ("Zündung      (kT_c = f α² m_p)² aus R⁻²", Tc**2),
    ("Thomson      σ_T⁻¹ = (3/8π) m_e²/α²", 1 / sig),
    ("Gravitation  (G M m_p)⁻¹ θ⁻², M = μ_M M_*", Mono(th=-2) / GMm),
    ("Strahlungsdruck 1-β", one_m_beta),
    ("Opazität     κ_T/κ̄", 1 / Mono(k=1)),
    ("Geometrie    4π/(4π σ_SB)", 1 / sSB),
]
prod = Mono(1.0)
for _, fm in factors:
    prod = prod * fm
assert np.allclose(prod.e, cond.e) and abs(prod.num / cond.num - 1) < 1e-12
print(f"   {'Glied':45s} {'α':>6s} {'m_e':>6s} {'m_p':>6s} {'μ_M':>6s}")
for name, fm in factors:
    e = fm.e * s + 0.0
    print(f"   {name:45s} {e[IDX['a']]:6.0f} {e[IDX['me']]:6.0f} {e[IDX['mp']]:6.0f} {e[IDX['mu']]:6.0f}")
print(f"   {'α_G = m_p²/m_Pl² (m_p²)':45s} {0:6.0f} {0:6.0f} {2:6.0f} {0:6.0f}")
print(f"   {'Summe':45s} {res.e[IDX['a']]:6.0f} {res.e[IDX['me']]:6.0f} {res.e[IDX['mp']]:6.0f} {res.e[IDX['mu']]:6.0f}")
print(f"   Faktor 2 = 2/|v_mPl|: α_G tritt als m_p/m_Pl = α_G^(1/2) auf (G M_* m_p = m_Pl/m_p).")
print(f"   Massenabhängigkeit: T_eff⁴ ∝ μ_M^{Teff4.e[IDX['mu']]:g}  ->  M_div/M_* ∝ 𝒞^(-1/2)")

print("\nVariante Kramers-Opazität (frei-frei, Rosseland): κ̄/κ_T = C_K·α·n·λ̄_e³·(m_e/kT_c)^(7/2), n = D·M/(m_p R³)")
n_ch = M / (Mono(mp=1) * R**3)
k_Kr = Mono(CK=1, a=1) * n_ch * Mono(me=-3) * (Mono(me=1) / Tc)**3.5
print("   κ̄/κ_T =", k_Kr.show(), "  (α kürzt sich heraus)")
lam_K = one_m_beta / k_Kr
Teff4_K = lam_K * L_Edd / (Mono(4 * pi) * R**2 * sSB)
res_K, sK = als_alphaG(Teff4_K / target)
print("   α_G =", res_K.show())
print(f"   -> α-Exponent {res_K.e[IDX['a']]:g}, m_e-Exponent {res_K.e[IDX['me']]:g}, m_p-Exponent {res_K.e[IDX['mp']]:g};"
      f"  T_eff⁴ ∝ μ_M^{Teff4_K.e[IDX['mu']]:g}  ->  M_div ∝ 𝒞'^(-1/6)")
print(f"   α¹²(m_e/m_p)⁵ = {alpha**12*(me/mp)**5:.4e};  α_G/(α¹²(m_e/m_p)⁵) = {aG/(alpha**12*(me/mp)**5):.1f}")

# =============================================================================================
h("Teil 3  Einfachste Form: optische Tiefe des Trichters = (T_c/T_ion)⁴")
print("Diffusionsrate der Strahlung Γ_diff = c/(τ R), τ = N σ_T/R²  (transversal);")
print("Rate, mit der eine Oberfläche bei T_ion denselben Vorrat aT_c⁴R³ leert: Γ_surf = (c/R)(T_ion/T_c)⁴.")
tau = Mono(mu=1) * Mono(mPl=3, mp=-3) * sig / R**2      # N = μ_M N_*, N_* = m_Pl³/m_p³
ratio4 = (Tc / Mono(fs=1, a=2, me=1))**4
cond3 = tau / ratio4
res3, _ = als_alphaG(cond3)
print("   τ =", tau.show())
print("   Γ_diff = Γ_surf  ⟺  τ = (T_c/T_ion)⁴  ⟺  α_G =", res3.show())
N_S = MS / mp
tau_S = 3 * N_S * sigT / (4 * pi * RS**2)
print(f"   Sonne: τ_T = 3Nσ_T/(4πR²) = {tau_S:.3e};  (T_c/T_eff)⁴ = {(TcS/TeffS)**4:.3e};"
      f"  Verhältnis {(TcS/TeffS)**4/tau_S:.0f} (Strukturfaktor: Gradient, Opazität)")
print(f"   In Faltungslängen: τ = (Δs/Δs_p)·(8π/3)r_e²/R²;  T_c/T_ion = (f_nuc/f_s)·(λ̄_e/λ̄_p)")

# =============================================================================================
h("Teil 4  Polytrope n = 3 (Eddington-Standardmodell)")


def lane_emden(n):
    def rhs(x, y):
        return [y[1], -y[0]**n - 2 * y[1] / x if y[0] > 0 else -2 * y[1] / x]
    x0 = 1e-6
    y0 = [1 - x0**2 / 6 + n * x0**4 / 120, -x0 / 3 + n * x0**3 / 30]
    ev = lambda x, y: y[0]
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, [x0, 20], y0, events=ev, rtol=1e-12, atol=1e-14)
    x1 = sol.t_events[0][0]
    dth = sol.y_events[0][0][1]
    return x1, dth


x1, dth = lane_emden(3.0)
m3 = -x1**2 * dth
rho_ratio = x1 / (3 * -dth)
Wn = 1 / (4 * pi * 4 * dth**2)
Tc_coef = Wn * (4 * pi / 3) / rho_ratio
C_E = 4 * m3 * math.sqrt(45) / pi**1.5
print(f"ξ₁ = {x1:.5f},  ξ₁²|θ'| = {m3:.5f},  ρ_c/ρ̄ = {rho_ratio:.3f},  kT_c = {Tc_coef:.4f}·μ_mol m_p GM/R")
print(f"Eddington-Quartik: M/M_* = C_E·μ_mol⁻²·√(1-β)/β²,  C_E = 4·ξ₁²|θ'|·√45/π^(3/2) = {C_E:.4f}"
      f"  ({C_E*Mstar/MS:.2f} M☉/μ_mol²)")
x15, dth15 = lane_emden(1.5)
print(f"(n = 1,5: ξ₁ = {x15:.4f}, ρ_c/ρ̄ = {x15/(3*-dth15):.3f})")


def beta_of(Mrel, mumol):
    return brentq(lambda b: C_E * mumol**-2 * math.sqrt(1 - b) / b**2 - Mrel, 1e-9, 1 - 1e-15)


for mu_m in (0.60, 0.62):
    for Mm in (1.0, 1.3, 10.0):
        b = beta_of(Mm * MS / Mstar, mu_m)
        print(f"   μ_mol = {mu_m:.2f}, M = {Mm:5.1f} M☉:  1-β = {1-b:.4e}"
              f"   (Näherung (μ_M μ_mol²/C_E)² = {(Mm*MS/Mstar*mu_m**2/C_E)**2:.4e})")
Mmax = C_E * Mstar / 0.62**2 * math.sqrt(0.5) / 0.25
print(f"   Strahlungsdruck = Gasdruck (β = 1/2) bei M = {Mmax/MS:.0f} M☉ = {Mmax/Mstar:.0f} M_* (μ_mol = 0,62)")

# =============================================================================================
h("Teil 5  Vorfaktor 2,94 mit Sonnenzahlen")
f_s = kB * TeffS / (alpha**2 * me * c**2)
f_nuc = kB * TcS / (alpha**2 * mp * c**2)
theta = kB * TcS * RS / (G * MS * mp)
mu_M = MS / Mstar
lam = LS / LEdd(MS)
print(f"f_s = kT_eff/(α²m_ec²) = {f_s:.5f}   f_nuc = kT_c/(α²m_pc²) = {f_nuc:.5f}   θ = kT_cR/(GMm_p) = {theta:.4f}")
print(f"μ_M = M/M_* = Δs☉/Δs_* = {mu_M:.4f}   λ = L/L_Edd = {lam:.4e}   ε☉ = {eps_S:.4e}"
      f"   (f_nuc/θ)α² = {f_nuc/theta*alpha**2:.4e}")
C_id = ((2 * pi**3 / 45) * f_s**4 * mu_M * (alpha**2 / eps_S)**2 / lam)**2
print(f"\nExakte Identität (nur L = 4πR²σT_eff⁴ umgeschrieben):")
print(f"   √𝒞 = (2π³/45)·f_s⁴·(Δs/Δs_*)·[(r_e/λ̄_e)²/(3Δs/R)]²/(L/L_Edd)  ->  𝒞 = {C_id:.4f}   (direkt: {CARTER:.4f})")
print("   (T_c kürzt sich: θ/f_nuc = α²/ε. Die Identität gilt für jeden Stern mit seinen eigenen Zahlen.)")

print("\nZerlegung mit dem Eddington-Standardmodell (λ = (κ_T/κ̄)(1-β)):")
Gam_log = 3 * (EG / (4 * kB * TcS))**(1 / 3)
print(f"   Gamow-Logarithmus Λ_G = 3(E_G/4kT_c)^(1/3) = {Gam_log:.2f};  f_nuc = π²/(4(Λ_G/3)³) = {pi**2/(4*(Gam_log/3)**3):.5f}")
print(f"   f_s als Logarithmus: kT_eff = (Ry)/Λ_s mit Λ_s = 1/(2 f_s) = {1/(2*f_s):.1f}")
for mu_m in (0.60, 0.62):
    b = beta_of(mu_M, mu_m)
    kbar = (1 - b) / lam
    A = 2 * pi**3 / 45
    B = C_E**2
    Dg = theta**2 / mu_m**4
    E = f_s**4 / f_nuc**2
    F = 1 / mu_M
    # exakte Version mit (1-β) statt Näherung:
    rootC = A * f_s**4 * theta**2 * mu_M / (lam * f_nuc**2)
    print(f"   μ_mol = {mu_m:.2f}: 1-β = {1-b:.4e},  erforderliches κ̄/κ_T = (1-β)/λ = {kbar:.1f}"
          f"  ->  κ̄ = {kbar*sigT/mp*10:.2f} cm²/g  (κ_T = σ_T/m_p = {sigT/mp*10:.3f} cm²/g)")
    print(f"      Gruppen: A=2π³/45={A:.3f}  B=C_E²={B:.1f}  C=κ̄/κ_T={kbar:.1f}  D=θ²/μ_mol⁴={Dg:.3f}"
          f"  E=f_s⁴/f_nuc²={E:.4e}  F=1/μ_M={F:.3f}")
    corr = (mu_M**2 * mu_m**4 / C_E**2) / (1 - b)
    print(f"      A·B·C·D·E·F = {A*B*kbar*Dg*E*F:.4f};  mal (Näherung/exaktes 1-β) = {corr:.4f}:"
          f" {A*B*kbar*Dg*E*F*corr:.4f} = √𝒞 = {math.sqrt(CARTER):.4f}  (Kontrolle {rootC:.4f})")
kappa_ff_c = 6.4e22 * (rhocS / 1e3) * TcS**-3.5
print(f"   Vergleich: Kramers frei-frei (reiner Wasserstoff, Rosseland) im Sonnenzentrum {kappa_ff_c:.2f} cm²/g;")
print("   OPAL-Werte der Sonne: Zentrum rund 1,2 cm²/g, Unterkante der Konvektionszone rund 20 cm²/g.")
print(f"   Kompensation: A·B·D·F = {2*pi**3/45*C_E**2*theta**2/0.61**4/mu_M:.3e} (groß) gegen E = {f_s**4/f_nuc**2:.3e} (Logarithmen)")

print("\nEmpfindlichkeit:")
print(f"   𝒞 ∝ f_s⁸: Faktor 2,94 in 𝒞 = Faktor {CARTER**(1/8):.4f} in T_crit")
print(f"   M_div ∝ 𝒞^(-1/2) (Thomson): Faktor {CARTER**0.5:.3f};   M_div ∝ 𝒞'^(-1/6) (Kramers): Faktor {CARTER**(1/6):.3f}")
print(f"   Reine Skala (alle Strukturzahlen = 1): T_eff(M_*) = [(45/2π³)·α²m_p³m_e²/m_Pl]^(1/4) ="
      f" {((45/(2*pi**3))*alpha**2*(mp*c**2)**3*(me*c**2)**2/(mPl*c**2))**0.25/eV:.1f} eV"
      f" = {((45/(2*pi**3))*alpha**2*(mp*c**2)**3*(me*c**2)**2/(mPl*c**2))**0.25/(alpha**2*me*c**2):.3f} α²m_ec²")

# =============================================================================================
h("Teil 6  Grenzmasse, Hayashi-Grenze, Massenfenster")
print("Grenzmasse aus Homologie ab der Sonne (Strukturzahlen fest): T_eff ∝ M^(1/4) (Thomson), M^(3/4) (Kramers)")
for Tcrit in (6200.0, 6500.0, 7000.0):
    mT = (Tcrit / TeffS)**4
    mK = (Tcrit / TeffS)**(4 / 3)
    print(f"   T_crit = {Tcrit:.0f} K (f_s = {kB*Tcrit/(alpha**2*me*c**2):.4f}):  M_div = {mT:.2f} M☉ (Thomson) = {mT*MS/Mstar:.2f} M_*"
          f"   |  {mK:.2f} M☉ (Kramers) = {mK*MS/Mstar:.2f} M_*")
print("   Beobachtet: Hüllenkonvektion verschwindet oberhalb rund 1,2–1,4 M☉ (T_eff ≈ 6200–7000 K) -> M_div ≈ 0,7 M_*")

print("\nHayashi-Grenze (voll konvektiv, n = 1,5; Photosphäre P = (2/3)g/κ, κ ∝ ρ^a T^b):")
for a_, b_ in ((0.5, 9.0),):
    # Unbekannte [lnP, lnρ, lnT, lnR]; rechte Seite [lnM, lnL]
    A_ = np.array([[1 + 0, 0, 0, 2], [1, -1, -1, 0], [1, 0, -2.5, 1.5], [0, 0, 4, 2]], float)
    A_[0] = [1, a_, b_, 2]
    Bm = np.array([[1, 0], [0, 0], [-0.5, 0], [0, 1]], float)
    X = np.linalg.solve(A_, Bm)
    print(f"   a = {a_}, b = {b_}:  T_eff ∝ M^{X[2,0]:.5f} · L^{X[2,1]:.5f}"
          + ("   (7/51 = 0,13725; 1/102 = 0,00980)" if b_ == 9.0 else ""))
EHm = 0.754 * eV
print(f"   H⁻-Bindung 0,754 eV = {EHm/(alpha**2*me*c**2):.4f}·α²m_ec²;  Hayashi-T_eff 3000–4500 K ->"
      f" f = {kB*3000/(alpha**2*me*c**2):.4f}–{kB*4500/(alpha**2*me*c**2):.4f}")

print("\nMassenfenster der Sterne in Einheiten M_*:")
Mmin_scale = (f_nuc * alpha**2 * mp / me)**0.75
print(f"   untere Grenze (Entartung vor Zündung), Skala (f_nuc α² m_p/m_e)^(3/4) M_* = {Mmin_scale:.4f} M_*"
      f" = {Mmin_scale*Mstar/MS:.3f} M☉;  beobachtet 0,075 M☉ = {0.075*MS/Mstar:.3f} M_*")
print(f"   obere Grenze (β = 1/2): {Mmax/Mstar:.0f} M_*;  beobachtet rund 150 M☉ = {150*MS/Mstar:.0f} M_*")
Mdiv = 1.3
lo = CARTER * (Mdiv / 150)**2
hi = CARTER * (Mdiv / 0.075)**2
print(f"   Schwache Form (M_div liegt im Fenster), Thomson-Skalierung: {lo:.1e} < 𝒞 < {hi:.1e}"
      f"  (Breite {math.log10(hi/lo):.1f} Dekaden in α_G)")

# =============================================================================================
h("Teil 7  Wasserstoff-Ionisation und adiabatischer Gradient")


def nabla_ad(x, y):          # y = χ/kT, reiner Wasserstoff (Kippenhahn-Weigert)
    q = x * (1 - x)
    return (2 + q * (2.5 + y)) / (5 + q * (2.5 + y)**2)


def saha_x(T, nH):
    S = (me * kB * T / (2 * pi * hbar**2))**1.5 * math.exp(-Ry / (kB * T)) / nH
    return (-S + math.sqrt(S * S + 4 * S)) / 2


nH = 1e23
print(f"n_H = {nH:.0e} m⁻³ (ρ ≈ {nH*mp*1e-3:.1e} g/cm³, Photosphärendichte):")
for T in (5000, 6500, 8000, 10000, 12000, 14000, 20000):
    x = saha_x(T, nH)
    print(f"   T = {T:6d} K  (f = {kB*T/(alpha**2*me*c**2):.4f}):  x = {x:.3e}   ∇_ad = {nabla_ad(x, Ry/(kB*T)):.3f}")
Th = brentq(lambda T: saha_x(T, nH) - 0.5, 5e3, 3e4)
print(f"   Halbe Ionisation bei T = {Th:.0f} K, kT = Ry/{Ry/(kB*Th):.1f}")
print("   ∇_rad = (λ/4)(κ/κ_T)(P_gas/P_rad): Faltung (λ, P_gas/P_rad ∝ (Δs_*/Δs)²) und EM-Opazität;")
print("   Konvektion setzt ein bei ∇_rad > ∇_ad, d. h. die Auftriebsfrequenz N² ∝ (∇_ad − ∇) geht durch null.")

# =============================================================================================
h("Teil 8  Torus-Dynamik: Faltungsrate gegen EM-Raten des Wasserstoffs")
print(f"Γ_G = α_G ω_p = {GammaG:.4e} 1/s;  ω_p = {wp:.4e};  ω_e = {we:.4e};  ω_a = α²ω_e = {wa:.4e} 1/s")
print(f"Ratenform der Carter-Beziehung: Γ_G·ω_p³ = 𝒞·(α³ω_e)⁴:  {GammaG*wp**3:.4e} = {CARTER:.4f} × {(alpha**3*we)**4:.4e}")
print(f"   (α_G^(1/4) ω_p = {aG**0.25*wp:.4e} 1/s,  α³ω_e = α·ω_a = {alpha**3*we:.4e} 1/s)")
print(f"\nEinrasten dΔθ/dt = Δω − K sin Δθ für ein einzelnes Atom, K = Γ_G, Δω ~ ω_a:  K/Δω = {GammaG/wa:.3e}")

A_Lya = (2 / 3)**8 * alpha**5 * we
nu21 = 1420.405751768e6
A21_SI = (2 * pi * nu21)**3 * mu0 * muB**2 / (3 * pi * hbar * c**3)
A21_mono = (64 / 81) * g_p**3 * alpha**13 * (me / mp)**3 * we
w_hf = 4 / 3 * g_p * alpha**4 * (me / mp) * we / (1 + me / mp)**3
rates = [
    ("Umlauf (Rydberg)  α²ω_e/2", 2, 0, 0.5 * alpha**2 * we),
    ("Feinstruktur 2p   α⁴ω_e/32", 4, 0, alpha**4 * we / 32),
    ("Lyman-α-Zerfall   (2/3)⁸α⁵ω_e", 5, 0, A_Lya),
    ("2s→1s, zwei Photonen (8,229 1/s)", 8, 0, 8.229),
    ("Hyperfein 21 cm   (4/3)g_p α⁴(m_e/m_p)ω_e", 4, 1, w_hf),
    ("21-cm-Zerfall     (64/81)g_p³α¹³(m_e/m_p)³ω_e", 13, 3, A21_SI),
]
print("\nRaten des Wasserstoffs als α^n (m_e/m_p)^k ω_e  (Carter als Γ_G = Γ_X verlangt n = 12, k = 3):")
for nm, n_, k_, val in rates:
    print(f"   {nm:48s} n = {n_:2d}, k = {k_}:  {val:.4e} 1/s   Γ_G/Γ = {GammaG/val:.3e}")
print(f"   benötigt: α¹²(m_e/m_p)³ω_e = {alpha**12*(me/mp)**3*we:.4e} 1/s (mit 𝒞: {CARTER*alpha**12*(me/mp)**3*we:.4e})")
print(f"\n21-cm-Linie: A₁₀ (SI, gemessenes ν) = {A21_SI:.4e} 1/s;  Monomform = {A21_mono:.4e} 1/s")
print(f"   Γ_G/A₁₀ = {GammaG/A21_SI:.4f};  Zerlegung Γ_G/A₁₀ = 𝒞·81/(64 g_p³ α) = {CARTER*81/(64*g_p**3*alpha):.4f}")
print(f"   (64/81)·g_p³·α = {(64/81)*g_p**3*alpha:.5f}: g_p³ (QCD) liegt zufällig bei 81/(64α); Exponent 13 statt 12.")

print("\nKollektiv: α_G tritt nur über G M m_p = μ_M·α_G^(-1/2)·ħc auf (N_* = α_G^(-3/2)).")
print(f"   N_*·Γ_G = α_G^(-1/2)ω_p = {aG**-0.5*wp:.3e} 1/s  (kein Zyklus, reine Buchhaltung)")
Rdyn = math.sqrt(G * MS / RS**3)
print(f"   Sonne: dynamische Frequenz √(GM/R³) = {Rdyn:.3e} 1/s")
tdiff = 3 * N_S * sigT / (4 * pi * RS**2) * RS / c
print(f"   Thomson-Diffusionszeit τR/c = {tdiff/3.156e7:.2e} a  (Rate {1/tdiff:.3e} 1/s)")
