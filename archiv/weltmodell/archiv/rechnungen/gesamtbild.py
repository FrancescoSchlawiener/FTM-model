"""Gesamtbild des FTM: Variablen, Längenleiter, 2π-Zählung, Querverbindungen, Topologie-Zelle.
Nur numpy. Alle Zahlen der Notiz 29-gesamtbild stammen aus dieser Ausgabe."""
import math
import numpy as np

P = print
def H1(t): P("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

# ---------------------------------------------------------------- Konstanten
G, c, hbar, kB, eV = 6.67430e-11, 299792458.0, 1.054571817e-34, 1.380649e-23, 1.602176634e-19
mp, me, alpha = 1.67262192e-27, 9.1093837e-31, 7.2973525693e-3
MPC, KPC, PC, AU, LY, YR = 3.0857e22, 3.0857e19, 3.0857e16, 1.495978707e11, 9.4607e15, 3.15576e7
Msun = 1.989e30
h = 0.6736; H0 = 100 * h * 1e3 / MPC
Om, Or = 0.3153, 9.14e-5; OL = 1 - Om - Or
zstar = 1089.9
lP = math.sqrt(hbar * G / c**3); mPl = math.sqrt(hbar * c / G)
RH = c / H0
a0 = c * H0 / (2 * math.pi)          # Definition des Modells (A3)
a0_fit = 1.2e-10                     # Anpassung RAR (Daten)
aG = G * mp**2 / (hbar * c)
TCMB = 2.7255
def E(z, Om=Om, Or=Or):
    OL = 1 - Om - Or
    return np.sqrt(Om * (1 + z)**3 + Or * (1 + z)**4 + OL)
def Hz(z): return H0 * E(z)

H1("0) Grundgrößen")
P(f"H0 = {100*h:.2f} km/s/Mpc, R_H = c/H0 = {RH:.4e} m = {RH/MPC/1e3:.4f} Gpc")
P(f"a0 = cH0/2π = c²/(2πR_H) = {a0:.4e} m/s²   (a0_fit = {a0_fit:.2e}; Verhältnis a0_fit/a0 = {a0_fit/a0:.4f})")
P("a0_fit = 1,20 ± 0,02 (stat) ± 0,24 (sys) ·1e-10 m/s² (McGaugh, Lelli, Schombert 2016, Eingabe)")
P(f"Ω_s/Ω_Λ(0,6847) = {math.pi**2/45/0.6847:.3f};  kosmische Varianz σ(C2)/C2 = √(2/5) = {math.sqrt(0.4):.3f};  Bedarf Nachlauf (Bullet Cluster, Eingabe): 1/λ ≥ 170–260 Mio a")
Mmw = 1.193e41
for lab, aa in [("a0_fit = 1,2e-10", 1.2e-10), ("cH0/2π, H0 = 70", c*70e3/MPC/(2*math.pi)), ("cH0/2π, H0 = 67,36", a0)]:
    P(f"  r* mit {lab:20s}: Milchstraße {math.sqrt(G*Mmw/aa)/KPC:.2f} kpc, Sonne {math.sqrt(G*Msun/aa)/AU:.0f} AE = {math.sqrt(G*Msun/aa)/LY:.3f} Lj")
P(f"a0 mit H0 = 70: {c*70e3/MPC/(2*math.pi):.4e};  mit H0 = 73: {c*73e3/MPC/(2*math.pi):.4e}")
P(f"ℓ_P = {lP:.5e} m, m_Pl = {mPl:.5e} kg, α_G = G m_p²/(ħc) = {aG:.5e}, N = R_H/ℓ_P = {RH/lP:.4e}")

# ---------------------------------------------------------------- 1) Längenleiter
H1("1) Die Längenleiter eines Trichters (a0 = cH0/2π)")
# Name, M [kg], R [m] (Hülle), r_J [m] oder None, J/ħ (Teilchen) oder None
Rp = 0.8409e-15
re = alpha * hbar / (me * c)
# Gezeitengrenzen
rJ_earth = AU * (5.972e24 / (3 * Msun))**(1/3)
A_oort, B_oort = 15.3e3 / KPC, -11.9e3 / KPC
rJ_sun = (G * Msun / (4 * A_oort * (A_oort - B_oort)))**(1/3)
M_M31, D_M31 = 1.5e11 * Msun, 765 * KPC
rJ_mw = D_M31 * (1.193e41 / (3 * M_M31))**(1/3)
bodies = [
    ("Elektron", me, re, None, 0.5),
    ("Proton", mp, Rp, None, 0.5),
    ("Erde", 5.972e24, 6.371e6, rJ_earth, None),
    ("Sonne", Msun, 6.957e8, rJ_sun, None),
    ("Milchstraße", 1.193e41, 4.730e20, rJ_mw, None),
    ("Haufen", 3.381e44, 6.171e22, None, None),
    ("Universum", 1.5e53, 4.399e26, None, None),
]
P(f"Gezeitengrenzen: Erde (Hill, 1 AE) {rJ_earth/AU:.4f} AE; Sonne (Oort A={A_oort*KPC/1e3}, B={B_oort*KPC/1e3}) "
  f"{rJ_sun/PC:.3f} pc = {rJ_sun/AU:.4e} AE; Milchstraße (Hill gegen M31, M_bar = 1,5e11 M☉, D = 765 kpc) {rJ_mw/KPC:.0f} kpc")
P(f"\n{'Körper':12s} {'λ̄ [m]':>11s} {'Δs [m]':>11s} {'r0=6Δs':>11s} {'R [m]':>11s} {'r* [m]':>11s} {'r_J [m]':>11s} {'R_H [m]':>11s}")
rows = {}
for n, M, R, rJ, s in bodies:
    lam = hbar / (M * c); ds = G * M / (3 * c**2); r0 = 6 * ds; rs = math.sqrt(G * M / a0)
    rows[n] = (M, lam, ds, r0, R, rs, rJ, s)
    P(f"{n:12s} {lam:11.3e} {ds:11.3e} {r0:11.3e} {R:11.3e} {rs:11.3e} {(rJ if rJ else float('nan')):11.3e} {RH:11.3e}")

P("\nAnschaulich:")
for n in rows:
    M, lam, ds, r0, R, rs, rJ, s = rows[n]
    def f(x):
        if x is None: return "–"
        if x < 1e-3: return f"{x:.3e} m"
        if x < 1e5: return f"{x:.4g} m"
        if x < 1e13: return f"{x/1e3:.4g} km"
        if x < 3e16: return f"{x/AU:.4g} AE"
        if x < 3e19: return f"{x/PC:.4g} pc"
        if x < 3e22: return f"{x/KPC:.4g} kpc"
        return f"{x/MPC:.4g} Mpc"
    P(f"  {n:12s} Δs = {f(ds):>14s} r0 = {f(r0):>14s} R = {f(R):>14s} r* = {f(rs):>14s} r_J = {f(rJ):>14s}")

P("\nUniverselle Beziehungen (für jeden Körper, Quotient links/rechts):")
P(f"{'Körper':12s} {'Δs·λ̄/(ℓP²/3)':>14s} {'r*²/(r0πR_H)':>13s} {'3Δs/λ̄/α_G(M)':>14s} {'δ(r*)/(v/c)²':>13s} {'GM/v²/r*':>9s} {'q(r*)':>8s} {'ν(r*)':>7s}")
for n in rows:
    M, lam, ds, r0, R, rs, rJ, s = rows[n]
    aGM = G * M**2 / (hbar * c)
    v = (G * M * a0)**0.25
    xx = r0 / rs
    delta_ex = xx / (1 + math.sqrt(1 - xx)) if xx < 1 else float('nan')
    P(f"{n:12s} {ds*lam/(lP**2/3):14.10f} {rs**2/(r0*math.pi*RH):13.10f} {3*ds/lam/aGM:14.10f} "
      f"{(3*ds/rs)/(v/c)**2:13.10f} {G*M/v**2/rs:9.6f} {math.exp(-1):8.5f} {1/(1-math.exp(-1)):7.4f}   exakt δ(r*)/(v/c)² = {delta_ex/(v/c)**2:.6f} (r0/r* = {xx:.3e})")
P("  (Spalte δ(r*) in erster Ordnung 3Δs/r*; exakt δ = 1-√(1-6Δs/r*), rechts)")

P("\nMassenabhängige Verhältnisse (μ = M/m_Pl):")
P(f"{'Körper':12s} {'μ':>11s} {'α_G(M)=μ²':>11s} {'r*/λ̄':>11s} {'R/r0':>11s} {'r*/R':>11s} {'r_J/r*':>9s} {'R_H/r*':>11s} {'g(R)/a0':>11s} {'q(R)':>10s}")
for n in rows:
    M, lam, ds, r0, R, rs, rJ, s = rows[n]
    mu = M / mPl; gR = G * M / R**2
    qR = math.exp(-math.sqrt(gR / a0))
    P(f"{n:12s} {mu:11.3e} {mu**2:11.3e} {rs/lam:11.3e} {R/r0:11.3e} {rs/R:11.3e} "
      f"{(rJ/rs if rJ else float('nan')):9.3g} {RH/rs:11.3e} {gR/a0:11.3e} {qR:10.3e}")

P("\nLeiter in Planck-Einheiten: λ̄ = ℓP/μ, Δs = ℓP·μ/3, r0 = 2ℓP·μ, r* = ℓP·√(2π·μ·N), R_H = N·ℓP")
N = RH / lP
for n in rows:
    M, lam, ds, r0, R, rs, rJ, s = rows[n]
    mu = M / mPl
    P(f"  {n:12s} log10(λ̄/ℓP) = {math.log10(lam/lP):8.3f}  log10(Δs/ℓP) = {math.log10(ds/lP):8.3f}  "
      f"log10(r*/ℓP) = {math.log10(rs/lP):8.3f}  Probe r*/(ℓP√(2πμN)) = {rs/(lP*math.sqrt(2*math.pi*mu*N)):.10f}")
P(f"  log10(R_H/ℓP) = {math.log10(N):.3f}")

P("\nBesondere Massen (Gleichheit zweier Längen):")
M_pl_eq = mPl / math.sqrt(2)
P(f"  λ̄ = r0          : μ = 1/√2,  M = {M_pl_eq:.4e} kg")
Mc0 = (hbar**2 * a0 / (G * c**2))**(1/3)
HL = H0 * math.sqrt(OL)
Mc_L = (hbar**2 * (c*HL/(2*math.pi)) / (G * c**2))**(1/3)
P(f"  λ̄ = r*          : M_c = (ħ²a0/(G c²))^(1/3) = {Mc0:.4e} kg = {Mc0*c**2/eV/1e6:.2f} MeV/c² = {Mc0/me:.1f} m_e = {Mc0/mp:.4f} m_p (heute);"
      f"  im de-Sitter-Endzustand {Mc_L*c**2/eV/1e6:.2f} MeV/c²")
P(f"                     Weinberg-Form (ħ²H0/(Gc))^(1/3) = {(hbar**2*H0/(G*c))**(1/3)*c**2/eV/1e6:.2f} MeV/c²  (Pion 139,57 MeV)")
rs_p = math.sqrt(G * mp / a0); lam_p = hbar / (mp * c)
P(f"  Proton heute: r*_p/λ̄_p = {rs_p/lam_p:.3f};  r*_p = λ̄_p ⟺ a0(z) = G m_p³c²/ħ² ⟺ H(z)/2π = α_G ω_p  (Probe: {2*math.pi*G*mp**3*c**2/hbar**2/c:.4e} = 2π·α_G·ω_p = {2*math.pi*aG*mp*c**2/hbar:.4e})")
MH = c**3 / (2 * G * H0)
P(f"  r0 = R_H/π (r* = R_H): M = {c**2*RH/(2*math.pi*G):.4e} kg;  r0 = R_H (Hubble-Masse M_H = c³/2GH0): {MH:.4e} kg, r*(M_H)/R_H = {math.sqrt(G*MH/a0)/RH:.6f} = √π = {math.sqrt(math.pi):.6f}")
P(f"  r* = r0 ⟺ r0 = π·R_H: M = {math.pi*c**2*RH/(2*G):.4e} kg")
MU = 1.5e53; r0U = 2 * G * MU / c**2
P(f"  Beobachtbares Universum: r0/R_H = {r0U/RH:.4f}, r*/r0 = {math.sqrt(G*MU/a0)/r0U:.4f}, r*/R_H = {math.sqrt(G*MU/a0)/RH:.4f}")

P("\nSkalenleiter mit dem Lineal λ̄_p (3Δs/λ̄_p = (M/m_p)·α_G):")
Mstar = mPl**3 / mp**2
for n, M in [("Proton", mp), ("Planck-Masse", mPl), ("Sternskala M_*", Mstar), ("Sonne", Msun),
             ("Milchstraße", 1.193e41), ("Haufen", 3.381e44), ("Hubble-Masse M_H", MH), ("Universum (beob.)", 1.5e53)]:
    x = 3 * (G * M / (3 * c**2)) / lam_p
    P(f"  {n:18s} 3Δs/λ̄_p = {x:10.3e}   Exponent zur Basis α_G: {math.log(x)/math.log(aG):8.4f}")
ds_star = G * Mstar / (3 * c**2); ds_p = G * mp / (3 * c**2)
P(f"  Probe Δs_*²·Δs_p/(λ̄_p/3)³ = {ds_star**2*ds_p/(lam_p/3)**3:.10f};  M_* = {Mstar/Msun:.3f} M☉, Δs_* = {ds_star:.1f} m")
P(f"  R_H/λ̄_p = {RH/lam_p:.4e},  α_G^-1 = {1/aG:.4e},  Quotient (Dirac) = {RH/lam_p*aG:.1f}")

P("\nEntropie, Temperatur, Verwindung in Faltungsgrößen:")
for n in rows:
    M, lam, ds, r0, R, rs, rJ, s = rows[n]
    aGM = G * M**2 / (hbar * c)
    S = 36 * math.pi * (ds / lP)**2
    kTH = hbar * c / (24 * math.pi * ds)
    extra = ""
    if s:
        a_tw = s * lam
        extra = f"  Verwindung a = (J/ħ)λ̄ = {a_tw:.3e} m, a·Δs/(ℓP²/6) = {a_tw*ds/(lP**2/6):.6f}, Spinzahl χ = a/3Δs = {a_tw/(3*ds):.3e} = (J/ħ)/α_G(M) → {(s/aGM):.3e}"
    P(f"  {n:12s} S/k = {S:.4e} = 4π·α_G(M)·[{S/(4*math.pi*aGM):.10f}]  kT_H/(Mc²) = {kTH/(M*c**2):.3e} = 1/(8πα_G)·[{kTH/(M*c**2)*8*math.pi*aGM:.8f}]"
      f"  T_H/T_GH = {(kTH/hbar)/(H0/(2*math.pi)):.3e} = R_H/(2r0)·[{(kTH/hbar)/(H0/(2*math.pi))/(RH/(2*r0)):.8f}]{extra}")
P(f"  S = 1·k bei M = m_Pl/√(4π) = {mPl/math.sqrt(4*math.pi):.4e} kg;  1 bit (S = k ln2) bei {mPl*math.sqrt(math.log(2)/(4*math.pi)):.4e} kg")
M1, M2 = Msun, Msun
S = lambda M: 4 * math.pi * G * M**2 / (hbar * c)
P(f"  Vereinigung zweier Falten (Δs additiv): S(M1+M2) - S(M1) - S(M2) = 8πG M1 M2/(ħc) = {S(M1+M2)-S(M1)-S(M2):.4e} k (2×1 M☉) = {100*(S(M1+M2)-S(M1)-S(M2))/S(M1+M2):.1f} % von S(M1+M2)")

# ---------------------------------------------------------------- 2) 2π-Zählung
H1("2) Die Zahl 2π: Zählung und Schließung")
z4 = sum(1.0 / n**4 for n in range(1, 200000))
P(f"ζ(4) = {z4:.12f}, π⁴/90 = {math.pi**4/90:.12f}, (2π)⁴/1440 = {(2*math.pi)**4/1440:.12f}")
Q0 = 4 * 6 * z4
P(f"∫0^∞ 4t³ Σ_n e^(-nt) dt = 24ζ(4) = {Q0:.10f} = 4π⁴/15 = {4*math.pi**4/15:.10f}")
Ls = Q0 * (a0 / c**2)**2
Os = Ls * c**2 / (3 * H0**2)
P(f"Λ_s = 24ζ(4)(a0/c²)² = {Ls:.4e} m⁻²,  Ω_s = Λ_s c²/(3H0²) = {Os:.10f} = π²/45 = {math.pi**2/45:.10f} = (2π)²/180 = {(2*math.pi)**2/180:.10f}")
P(f"a0 · (2π/H0) = {a0*2*math.pi/H0:.6e} m/s = c ({c:.6e}):  Die Schwelle baut in einer Periode 2π/H0 die Geschwindigkeit c auf.")
TGH = hbar * H0 / (2 * math.pi * kB)
P(f"Gibbons-Hawking T_GH(H0) = ħH0/(2πk) = {TGH:.4e} K;  a0 = c·kT_GH/ħ: {c*kB*TGH/hbar:.4e} = a0 {a0:.4e}")
P(f"Oberflächenschwere des Hubble-Horizonts κ_dS = c²/R_H = {c**2/RH:.4e};  κ_dS/2π = {c**2/RH/(2*math.pi):.4e} = a0")
P(f"Unruh-Temperatur der Schwelle kT_U(a0) = ħa0/(2πc) → T_U/T_GH = {(hbar*a0/(2*math.pi*c*kB))/TGH:.6f} = 1/2π = {1/(2*math.pi):.6f}")
P(f"Lichtumlauf auf Umfang L = 2πR_H: T = L/c = {2*math.pi*RH/c:.4e} s = 2π/H0 = {2*math.pi/H0:.4e} s; Grundfrequenz f1 = c/L = {c/(2*math.pi*RH):.4e} Hz = H0/2π; a0 = c·f1 = {c*c/(2*math.pi*RH):.4e}")
P(f"Spektrallücke Kreis L: D(2π/L)²; mit L = 2πR:  D/R² (kein 2π, wenn der Radius R benutzt wird)")
P(f"Diffusions-Deutung a0 = c·κ_gap mit L = 2πR_H verlangt D = c·R_H/2π = {c*RH/(2*math.pi):.4e} m²/s (aus denselben Größen gebaut)")
P(f"Wellenoperator -c²∂² auf L = 2πR_H: ω1 = 2πc/L = {2*math.pi*c/(2*math.pi*RH):.4e} = H0 = {H0:.4e}")

# Kreuzung: H/2π = α_G ω_p (gemischte Zählung) gegen H = α_G ω_p (einheitlich)
wp = mp * c**2 / hbar; GamG = aG * wp
def zsolve(target):
    lo, hi = 0.0, 1e7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if Hz(mid) < target: lo = mid
        else: hi = mid
    return lo
z_mix = zsolve(2 * math.pi * GamG)
z_uni = zsolve(GamG)
z_loop = zsolve(GamG / (2 * math.pi))
kT = lambda z: kB * TCMB * (1 + z) / eV
z_eq = Om / Or - 1
P(f"\nΓ_G = α_G ω_p = {GamG:.4e} 1/s")
P(f"gemischte Zählung H/2π = α_G ω_p :  z = {z_mix:.1f}, kT = {kT(z_mix):.4f} eV")
P(f"einheitliche Zählung H = α_G ω_p  (≡ H/2π = α_G f_p):  z = {z_uni:.1f}, kT = {kT(z_uni):.4f} eV")
P(f"Umlauf-Zählung (a0 pro Umlauf, volle Compton-Wellenlänge λ = 2πλ̄): H = α_G ω_p/2π:  z = {z_loop:.1f}, kT = {kT(z_loop):.4f} eV")
P(f"Geometrische Form: r*_p = λ̄_p ⟺ gemischte Zählung (z = {z_mix:.1f});  G m_p R_H/c² = λ̄_p² ⟺ H = α_G ω_p (z = {z_uni:.1f})")
P(f"Gleichheit Materie/Strahlung z_eq = Ω_m/Ω_r - 1 = {z_eq:.1f} (kT = {kT(z_eq):.4f} eV); letzte Streuung z* = {zstar} (kT = {kT(zstar):.4f} eV)")
P(f"m_p c² α_G^(1/4) = {mp*c**2*aG**0.25/eV:.4f} eV, α³ m_e c² = {alpha**3*me*c**2/eV:.4f} eV")
P(f"Nome: q = exp(-r*/r) enthält kein 2π; τ = i r*/(2πr) setzt es per Konvention. Selbstdualer Punkt τ = i: r = r*/2π, y = a_N/a0 = 4π² = {4*math.pi**2:.3f}")
yy = np.array([4*math.pi**2*0.5, 4*math.pi**2, 4*math.pi**2*2])
nu = 1 / (1 - np.exp(-np.sqrt(yy)))
P(f"  ν(y) um den selbstdualen Punkt (y = 19,7 / 39,5 / 79,0): {nu[0]:.6f} / {nu[1]:.6f} / {nu[2]:.6f}  (glatt, keine Signatur)")

# ---------------------------------------------------------------- 3) Querverbindungen
H1("3) Querverbindungen")
P("3a) Kac-Budget gegen Rückkopplungssumme")
Ob_h2, Oc_h2, sOb, sOc = 0.02237, 0.1200, 0.00015, 0.0012
Q = Oc_h2 / Ob_h2
sQ = Q * math.sqrt((sOb/Ob_h2)**2 + (sOc/Oc_h2)**2)
q = Q / (1 + Q); sq = sQ / (1 + Q)**2
P(f"  Q = Ω_c/Ω_b = {Q:.4f} ± {sQ:.4f};  q = Q/(1+Q) = Ω_c/Ω_m = {q:.5f} ± {sq:.5f};  1-q = f_b = Ω_b/Ω_m = {1-q:.5f};  ν = 1/(1-q) = Ω_m/Ω_b = {1/(1-q):.4f}")
P(f"  Kac mit zwei Zellen (Maß 1-q und q): ⟨τ_b⟩ = {1/(1-q):.4f}, ⟨τ_S⟩ = {1/q:.4f}, Σ1/⟨τ⟩ = {(1-q)+q:.6f}")
P(f"  Geometrische Flucht mit Rückkehrwahrscheinlichkeit q: mittlere Zahl der Rückläufe vor der Flucht q/(1-q) = {q/(1-q):.4f} = Q")
# Pólya/Watson: Rückkehrwahrscheinlichkeiten
P("  Rückkehrwahrscheinlichkeiten offener Gitter (Pólya): sc 3D 0,340537; bcc 0,282229; fcc 0,256318 (Literatur) — keine nahe q")
# Monte-Carlo-freie Kontrolle: 1D-Irrfahrt auf Ring mit Leckrate ε: Rückkehr-Summe
eps = 1 - q
P(f"  Leck pro Umlauf 1-q = {eps:.4f} ⇒ Summe Σ qⁿ = {sum(q**n for n in range(5000)):.4f}")

P("\n3b) Schleifenverstärkung Dreieck gegen q")
for taus in [(1, 1, 1), (1, 10, 100)]:
    s1, s2 = sum(taus), sum(1 / t for t in taus)
    P(f"  τ = {taus}: |q3|_krit = {s1*s2-1:.4f}")
P(f"  statische Schwelle Σq ≥ 1; kosmisches q = {q:.4f} < 1 (Reserve 1-q = {1-q:.4f}); |q3|_krit/q ≥ {8/q:.2f}")

P("\n3c) Hecke-Leiter gegen Skalenleiter")
P("  T_n: q → qⁿ ⟺ √y → n√y ⟺ y → n²y ⟺ r*→n·r* bei festem r ⟺ M → n²M")
for nm, ratio in [("Proton→Planck", mPl/mp), ("Planck→Sternskala", Mstar/mPl), ("Sonne→Milchstraße", 1.193e41/Msun), ("Milchstraße→Haufen", 3.381e44/1.193e41)]:
    P(f"  {nm:22s} Massenfaktor {ratio:.3e} ⟹ Hecke-Index n = √Faktor = {math.sqrt(ratio):.3e}")

P("\n3d) eV-Ebene gegen Frequenzleiter")
lad = [("m_e c² α⁰", me*c**2), ("m_e c² α¹", me*c**2*alpha), ("m_e c² α² (ω2)", me*c**2*alpha**2),
       ("m_e c² α³", me*c**2*alpha**3), ("m_e c² α⁴", me*c**2*alpha**4),
       ("m_p c² α_G⁰ (ω1)", mp*c**2), ("m_p c² α_G^1/4", mp*c**2*aG**0.25), ("m_p c² α_G^1/2", mp*c**2*aG**0.5), ("m_p c² α_G¹ (Γ_G)", mp*c**2*aG),
       ("ħ·10¹³ rad/s (ω3)", hbar*1e13), ("kT_rek (z*)", kB*TCMB*(1+zstar)), ("ħH0", hbar*H0)]
for nm, En in lad:
    P(f"  {nm:20s} E = {En/eV:12.4e} eV   ω = E/ħ = {En/hbar:11.4e} rad/s   log10(E/eV) = {math.log10(En/eV):8.3f}")
P(f"  Born-Oppenheimer: α² m_e c²·(m_e/m_p)^(1/2) = {alpha**2*me*c**2*(me/mp)**0.5/eV:.4f} eV; ·(m_e/m_p) = {alpha**2*me*c**2*(me/mp)/eV*1e3:.3f} meV;  ħω3 = {hbar*1e13/eV*1e3:.3f} meV")
# Zufallsprüfung Carter
la, lr, lg = math.log(alpha), math.log(me/mp), math.log(aG)
hits = [(a, b) for a in range(0, 25) for b in range(-12, 13) if abs(a*la + b*lr - lg) < math.log(3)]
P(f"  Zufallsprüfung Carter: Paare (a,b), a∈[0,24], b∈[-12,12] mit α^a(m_e/m_p)^b innerhalb Faktor 3 von α_G: {len(hits)} → {hits}")

P("\n3e) Zufallsprüfung q = 0,8429")
consts = {"π": math.pi, "e": math.e, "φ": (1+5**0.5)/2, "√2": 2**0.5, "√3": 3**0.5, "ln2": math.log(2), "γ": 0.5772156649, "ζ3": 1.2020569, "2π": 2*math.pi}
cand = {}
for a in range(1, 21):
    for b in range(a+1, 21):
        cand[f"{a}/{b}"] = a / b
for k, v in consts.items():
    for m in range(1, 13):
        for val, nm in [(v**(-1/m) if v > 1 else v**(1/m), f"{k}^(∓1/{m})"), (math.exp(-v/m), f"exp(-{k}/{m})"),
                        (1 - 1/(m*v), f"1-1/({m}{k})"), (v/(v+m), f"{k}/({k}+{m})"), (m/(m+v), f"{m}/({m}+{k})")]:
            cand[nm] = val
    for m in range(1, 13):
        cand[f"exp(-1/{m})"] = math.exp(-1/m)
cand["erf(1)"] = math.erf(1)
win = [(k, v) for k, v in cand.items() if abs(v - q) < sq]
win2 = [(k, v) for k, v in cand.items() if abs(v - q) < 2*sq]
dens = sum(1 for v in cand.values() if 0.7 < v < 0.95) / 0.25
P(f"  {len(cand)} Kandidaten; im 1σ-Fenster [{q-sq:.4f}, {q+sq:.4f}]: {len(win)}; im 2σ-Fenster: {len(win2)}; erwartet bei gleicher Dichte (1σ): {dens*2*sq:.1f}")
P("  1σ-Treffer: " + ", ".join(f"{k}={v:.5f}" for k, v in sorted(win, key=lambda t: abs(t[1]-q))[:12]))
P(f"  e^(-1/6) = {math.exp(-1/6):.5f} ({(math.exp(-1/6)-q)/sq:+.1f}σ), 1-1/2π = {1-1/(2*math.pi):.5f} ({(1-1/(2*math.pi)-q)/sq:+.1f}σ), erf(1) = {math.erf(1):.5f} ({(math.erf(1)-q)/sq:+.1f}σ)")
P(f"  zugehöriges y = (ln q)² = {math.log(q)**2:.5f}, r/r* = {-1/math.log(q):.3f}")

P("\n3f) Nachlaufzeit als Dynamikzeit an der Schwelle (Hypothese 1/λ = r*/v = v/a0)")
for nm, M in [("Zwerggalaxie 1e8 M☉", 1e8*Msun), ("Milchstraße 6e10", 1.193e41), ("Gruppe 1e13", 1e13*Msun), ("Haufen 1,7e14", 3.381e44), ("Haufen 1e15", 1e15*Msun)]:
    v = (G*M*a0)**0.25; rs = math.sqrt(G*M/a0)
    P(f"  {nm:22s} v_flach = {v/1e3:7.1f} km/s, r* = {rs/KPC:9.2f} kpc, r*/v = v/a0 = {v/a0/YR/1e6:7.1f} Mio a, 2π r*/v = {2*math.pi*v/a0/YR/1e9:6.2f} Mrd a")

# ---------------------------------------------------------------- 4) Topologie
H1("4) a0 als Grundfrequenz einer Torus-Schleife: L = c²/a0 und die letzte Streufläche")
def chi_over_RH(Om_, Or_=Or, z=zstar, n=200000):
    x = np.linspace(0, math.log(1+z), n)
    zz = np.exp(x) - 1
    f = np.exp(x) / E(zz, Om_, Or_)
    return np.trapezoid(f, x) if hasattr(np, "trapezoid") else np.trapz(f, x)
I0 = chi_over_RH(Om)
chi = I0 * RH; dLSS = 2 * chi
P(f"χ*/R_H = I(Ω_m={Om}) = {I0:.5f};  χ* = {chi/MPC/1e3:.4f} Gpc, d_LSS = {dLSS/MPC/1e3:.4f} Gpc;  I/π = {I0/math.pi:.5f}")
L_2pi = 2 * math.pi * RH
P(f"L = 2πR_H = {L_2pi/MPC/1e3:.4f} Gpc;  L/d_LSS = π/I = {math.pi/I0:.5f}  (H0 kürzt sich; nur Ω_m und Ω_r h² gehen ein)")
P("\nAbhängigkeit von Ω_m bei festem h = 0,6736 (Ω_r = 9,14e-5):")
for om in [0.20, 0.25, 0.28, 0.30, 0.3153, 0.33, 0.35, 0.40, 1.0]:
    I = chi_over_RH(om)
    P(f"  Ω_m = {om:6.4f}: I = χ*/R_H = {I:.4f}, L/d_LSS = π/I = {math.pi/I:.4f}")
P(f"  Einstein-de-Sitter analytisch 2(1-1/√(1+z*)) = {2*(1-1/math.sqrt(1+zstar)):.4f}")
# Ω_m bei I = π und bei L/d = 0,985
def solve_om(target, lo=0.1, hi=0.9):
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if chi_over_RH(mid, n=40000) > target: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)
om_pi = solve_om(math.pi); om_lim = solve_om(math.pi/0.985)
P(f"  I = π bei Ω_m = {om_pi:.4f};  L/d_LSS = 0,985 (COMPACT-Grenze) bei Ω_m = {om_lim:.4f}  → Deutung verlangt Ω_m ≥ {om_lim:.4f}")
P(f"  Planck Ω_m = 0,3153 ± 0,0073: Abstand zu {om_lim:.4f} = {(0.3153-om_lim)/0.0073:.2f} σ")
P("\nFeste Kombination Ω_m h² = 0,1430 (CMB), H0 variabel:")
for H0v in [64, 67.36, 70, 73]:
    hv = H0v/100; omv = 0.1430/hv**2; orv = 4.147e-5/hv**2
    I = chi_over_RH(omv, orv)
    P(f"  H0 = {H0v:5.2f}: Ω_m = {omv:.4f}, I = {I:.4f}, L/d_LSS = {math.pi/I:.4f}")
lo, hi = 60.0, 80.0
for _ in range(50):
    mid = 0.5*(lo+hi); hv = mid/100
    if math.pi/chi_over_RH(0.1430/hv**2, 4.147e-5/hv**2, n=40000) > 0.985: lo = mid
    else: hi = mid
P(f"  L/d_LSS = 0,985 bei H0 = {lo:.2f} km/s/Mpc  → Deutung verlangt H0 ≤ {lo:.2f}")
P("\nEpochenabhängigkeit (Beobachter beim Skalenfaktor a_o, gleiche Entkopplung a* = 1/(1+z*)):")
astar = 1/(1+zstar)
def chi_com(a_o, n=200000):
    x = np.linspace(math.log(astar), math.log(a_o), n)
    a = np.exp(x)
    Ea = np.sqrt(Om*a**-3 + Or*a**-4 + OL)
    f = 1/(a*Ea)          # c/H0 · ∫ da/(a² E) = ∫ dln a /(a E)
    return (np.trapezoid(f, x) if hasattr(np, "trapezoid") else np.trapz(f, x))   # in Einheiten R_H(heute)
for a_o in [0.25, 0.5, 1.0, 2.0, 4.0]:
    chic = chi_com(a_o)
    E_o = math.sqrt(Om*a_o**-3 + Or*a_o**-4 + OL)
    ratio_track = math.pi / (a_o*E_o*chic)     # 2πR_H(t)/d_LSS(t), beides physikalisch
    ratio_fixed = math.pi / chic               # L komobil fest = 2πR_H(heute)
    P(f"  a_o = {a_o:4.2f}: 2πR_H(t)/d_LSS(t) = {ratio_track:.4f};  feste Zelle 2πR_H,0/d_LSS(t) = {ratio_fixed:.4f}")

P("\nDiskrete Sachs-Wolfe-Modensumme im kubischen 3-Torus (P ∝ k^-3, ohne ISW):")
def jl(l, x):
    s, co = np.sin(x), np.cos(x)
    if l == 2: return (3/x**3 - 1/x)*s - 3*co/x**2
    if l == 3: return (15/x**4 - 6/x**2)*s - (15/x**3 - 1/x)*co
    if l == 4: return (105/x**5 - 45/x**3 + 1/x)*s - (105/x**4 - 10/x**2)*co
Nmax = 48
g = np.arange(-Nmax, Nmax+1)
n2 = (g[:, None, None]**2 + g[None, :, None]**2 + g[None, None, :]**2).ravel()
n2 = n2[(n2 > 0) & (n2 <= Nmax**2)]
vals, counts = np.unique(n2, return_counts=True)
nn = np.sqrt(vals)
def resid(L, l):
    k = 2*math.pi*nn/L
    S = (2*math.pi/L)**3/(4*math.pi) * np.sum(counts * k**-3 * jl(l, k*chi)**2)
    K = 2*math.pi*Nmax/L
    kk = np.linspace(K, K*400, 400000)
    tail = np.trapezoid(jl(l, kk*chi)**2/kk, kk) if hasattr(np, "trapezoid") else np.trapz(jl(l, kk*chi)**2/kk, kk)
    tail += 1/(4*(K*400*chi)**2)
    return (S + tail) / (1/(2*l*(l+1)))
def crude(L, l):
    x0 = 2*math.pi*chi/L
    xx = np.linspace(x0, 4000, 2000000)
    val = np.trapezoid(jl(l, xx)**2/xx, xx) if hasattr(np, "trapezoid") else np.trapz(jl(l, xx)**2/xx, xx)
    return (val + 1/(4*4000**2)) / (1/(2*l*(l+1)))
P(f"  Kontrolle L = 20·d_LSS: C2-Rest = {resid(20*dLSS, 2):.4f}")
cases = [("0,25·d", 0.25*dLSS), ("0,5·d", 0.5*dLSS), ("0,875·d (a0 = 1,20e-10)", c**2/1.2e-10), ("a0 = 1,10e-10", c**2/1.1e-10), ("0,985·d (COMPACT)", 0.985*dLSS),
         ("L = 2πR_H (a0 = cH0/2π)", L_2pi), ("a0 = 1,00e-10", c**2/1.0e-10), ("1,2·d", 1.2*dLSS), ("1,5·d", 1.5*dLSS), ("2,0·d", 2.0*dLSS)]
for nm, L in cases:
    x = L/dLSS
    alpha_c = math.degrees(math.acos(x)) if x < 1 else float('nan')
    P(f"  {nm:26s} L = {L/MPC/1e3:6.2f} Gpc, L/d = {x:.4f}, Kreisradius α = {alpha_c:5.1f}°, "
      f"C2-Rest = {resid(L,2):.3f}, C3-Rest = {resid(L,3):.3f}, C4-Rest = {resid(L,4):.3f}  | grob (Abschneiden bei k_min) C2 = {crude(L,2):.3f}")
P(f"  a0-Obergrenze aus L ≥ 0,985 d_LSS: a0 ≤ c²/(0,985 d_LSS) = {c**2/(0.985*dLSS):.4e} m/s²")

P("\nZwei Lesarten von a0(z):")
for z in [0.5, 1, 2, 3]:
    Ez = float(E(z))
    P(f"  z = {z}: Atmen a0 ∝ H(z): Faktor {Ez:.3f}, tiefer Bereich g_obs/g_bar ∝ √a0: +{100*(math.sqrt(Ez)-1):.1f} %;  "
      f"feste Zelle a0 = c²/L_phys ∝ (1+z): Faktor {1+z:.3f}, +{100*(math.sqrt(1+z)-1):.1f} %")

P("\nZufallsprüfung zu I ≈ π:")
Is = [chi_over_RH(om, n=40000) for om in (0.2, 0.4)]
P(f"  Ω_m ∈ [0,2; 0,4] überstreicht I ∈ [{Is[1]:.3f}; {Is[0]:.3f}]")
simple = {"3": 3.0, "π": math.pi, "√10": 10**0.5, "e+1/2": math.e+0.5, "2φ": 1+5**0.5, "√3·φ... (φ²+φ/2)": None}
simple = {k: v for k, v in simple.items() if v}
inside = {k: v for k, v in simple.items() if Is[1] <= v <= Is[0]}
P("  einfache Zahlen im Intervall: " + ", ".join(f"{k}={v:.4f}" for k, v in inside.items()))
dI = (chi_over_RH(Om+0.0073, n=40000) - chi_over_RH(Om-0.0073, n=40000))/2
P(f"  Unsicherheit von I aus σ(Ω_m) = 0,0073: ±{abs(dI):.4f};  |I-π| = {abs(I0-math.pi):.4f} = {abs(I0-math.pi)/abs(dI):.2f} σ")

# ---------------------------------------------------------------- 5) Parameter
H1("5) Parameterzählung")
P("Naturkonstanten: G, c, ħ, k (+ m_p, m_e, α für Teilchen); Atmen: H0.")
P("Pro Körper: M (dazu R, J, Umgebung für r_J).")
P("Galaxien (stationär): 0 freie Parameter (a0 = cH0/2π, q = exp(-√(a_N/a0))).")
P("Kosmos: Q = 5,36 (bzw. q), Form von w(θ) (mind. 1 Übergangsgröße), λ (Nachlauf), Rest-Dunkle-Energie (≥1), Zellgröße L (+Orientierung).")
