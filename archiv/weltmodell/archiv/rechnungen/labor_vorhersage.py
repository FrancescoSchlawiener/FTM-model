"""Laborvorhersage: gibt es eine Anordnung, in der nur dieses Modell etwas anderes sagt?

A  Regel im Labor: s = √(|g|/a₀) je Feldebene (L Erde, Sonne; G galaktisch, geglättet)
B  Sattelpunkt Sonne–Erde: Newton-Feld ungeglättet (MOND-artige Lesart) gegen geglättetes Feld (07, Hierarchie)
C  Ort, an dem das Gesamtfeld unter a₀ fällt (Sonne gegen Galaxie)
D  Takt: Uhren im Labor, Zusatz aus der Regel
E  Kippstelle: 1 − q = 2(3a)^{1/4}δ^{1/2} ist die Sattel-Knoten-Normalform (Schleife x = u + 2x − a x³);
   Wurzel der Rückkehr 1 − q ≈ √(Γτ) = Kabelgleichung
Konstanten aus ftm_kern; Werte der Milchstraße (g_N am Sonnenort, Gezeit T, Glättungslängen) wie in wep_glaettung.
"""
import math
import numpy as np
from scipy.optimize import brentq
from ftm_kern import G, c, A0, GM_sun, GM_earth, AU, R_earth, g_earth, nu, chi, s_, bose_stufe_max

out = []
def P(t=""): out.append(t); print(t)
def H(t): P(); P("=" * 96); P(t); P("=" * 96)

GN_GAL = 1.474e-10          # g_N der Galaxie am Sonnenort (Regel umgekehrt aus v₀ = 233 km/s, wep_glaettung)
T_GAL = 5.549e-30           # galaktische Gezeit am Sonnenort [s⁻²] (wep_glaettung)
ELL = {"Sonne": (GM_sun / T_GAL) ** (1 / 3), "Erde": (GM_earth / T_GAL) ** (1 / 3)}

def log10_q(g):
    return -float(s_(g)) / math.log(10)

# ---------------------------------------------------------------- A
H("A  Regel je Feldebene im Labor (a₀ = 1,1e-10)")
P(f"   {'Feld':38s} {'|g| [m/s²]':>11s} {'s':>10s} {'log10 q':>12s}  ν − 1")
for lab, g in [("Erdoberfläche (L, eigenes Feld)", g_earth),
               ("Sonne an der Erdbahn (L)", GM_sun / AU**2),
               ("Galaxie am Sonnenort, geglättet (G)", GN_GAL)]:
    s = float(s_(g)); nm1 = float(nu(g)) - 1
    P(f"   {lab:38s} {g:11.3e} {s:10.4g} {log10_q(g):12.5g}  {nm1:.4g}" if nm1 > 1e-300 else
      f"   {lab:38s} {g:11.3e} {s:10.4g} {log10_q(g):12.5g}  10^{log10_q(g):.5g}")
zus = (float(nu(GN_GAL)) - 1) * GN_GAL
P(f"   Zusatz der Galaxie (ν−1)·g_N = {zus:.3e} m/s², gleichförmig über das Sonnensystem;")
P(f"   Änderung über 1 AE ≤ T·AE·(ν + dν/dln g) ≈ {T_GAL*AU*2:.1e} m/s², über 1 km ≈ {T_GAL*1e3*2:.1e} m/s²")
P(f"   Glättungslängen ℓ = (GM/T)^(1/3): Sonne {ELL['Sonne']/AU:.3g} AE, Erde {ELL['Erde']/AU:.3g} AE")

# ---------------------------------------------------------------- B
H("B  Sattelpunkt Sonne–Erde")
R = AU
f = lambda x: GM_earth / x**2 - GM_sun / (R - x)**2     # x = Abstand von der Erde zur Sonne hin
xs = brentq(f, 1e7, 1e10)
gS = GM_earth / xs**2
Txx = 2 * GM_earth / xs**3 + 2 * GM_sun / (R - xs)**3    # Gradient längs der Achse
Tyy = GM_earth / xs**3 + GM_sun / (R - xs)**3            # quer (Betrag)
P(f"   Sattel {xs:.4e} m von der Erde; dort je Körper |g| = {gS:.3e} m/s²")
P(f"   Gradient längs {Txx:.3e} s⁻², quer {Tyy:.3e} s⁻²")
P(f"   Newton-Nullpunkt mit Galaxie: verschoben um g_gal/T = {GN_GAL/Txx:.2f} m (längs) bis {GN_GAL/Tyy:.2f} m (quer)")
P(f"   Bereich |g_N| < a₀: Halbachsen {A0/Txx:.2f} m (längs), {A0/Tyy:.2f} m (quer)")
smax, xmax, vmax = bose_stufe_max(2)
P(f"   (i) ungeglättete Lesart (Regel auf dem wirklichen Newton-Feld): größter Zusatz (ν−1)a_N = s²n_B·a₀ = {vmax*A0:.3e} m/s² bei a_N = {xmax:.2f} a₀")
for d in [10, 100, 300, 1000]:
    g = Txx * d
    P(f"       d = {d:5d} m längs: |g| = {g:.2e}, Zusatz (ν−1)|g| = {(float(nu(g))-1)*g:.2e} m/s²")
P(f"   (ii) Modell (07 Hierarchie): Erde bis ℓ_E = {ELL['Erde']/AU:.0f} AE, Sonne bis ℓ☉ = {ELL['Sonne']/AU:.3g} AE Newtonsch;")
gbar = GN_GAL + T_GAL * xs
P(f"        geglättetes Feld am Sattel = g_gal + T·d = {gbar:.4e} m/s² (s = {float(s_(gbar)):.4f}); Zusatz gleichförmig {zus:.3e},")
P(f"        sein Gradient ~ T = {T_GAL:.1e} s⁻² → über 1 km {T_GAL*1e3:.1e} m/s²: kein Sattelsignal. Newton-Dynamik am Sattel.")

# ---------------------------------------------------------------- C
H("C  Wo fiele das Gesamtfeld unter a₀?")
rb = math.sqrt(GM_sun / GN_GAL)
P(f"   Newton ungeglättet: Sonne hebt die Galaxie auf bei r = √(GM☉/g_gal) = {rb/AU:.0f} AE (Richtung Zentrum)")
P(f"   geglättet: Feld der Sonnenwolke ≤ GM☉/ℓ☉² = {GM_sun/ELL['Sonne']**2:.2e} ≪ g_gal: |ḡ| ≥ {GN_GAL/A0:.2f} a₀ überall im Sonnensystem")
P(f"   Kleinstes s erreichbar: {float(s_(GN_GAL)):.3f} (q = {math.exp(-float(s_(GN_GAL))):.3f}), und nur als gleichförmiger Zusatz")

# ---------------------------------------------------------------- D
H("D  Takt: Uhren")
Gam = g_earth / c
P(f"   Taktgefälle an der Erdoberfläche Γ = g/c = {Gam:.3e} s⁻¹ = Rotverschiebung {g_earth/c**2:.3e} je m (ART)")
P(f"   Zusatz der Regel am Erdfeld: ν − 1 = 10^{log10_q(g_earth):.6g}")
for h in [1.0, 1e3, 2e7]:
    P(f"   Uhrenpaar Abstand {h:8.0e} m, frei fallend im Galaxiefeld: Gezeitenanteil ≈ T h²/(2c²) = {T_GAL*h*h/2/c**2:.1e}")
P("   Gleichförmiger Zusatz der Galaxie: im frei fallenden Labor ohne Taktunterschied (Uhren am selben Ort fallen gleich).")

# ---------------------------------------------------------------- E
H("E  Kippstelle und Wurzelgesetz")
a = 0.7
xf = math.sqrt((2 - 1) / (3 * a)); uf = (1 - 2) * xf + a * xf**3      # u(x) = (1−k)x + a x³, k = 2; Falte bei du/dx = 0
for dl in [1e-2, 1e-4, 1e-6]:
    x = brentq(lambda x: (1 - 2) * x + a * x**3 - (uf + dl), xf, 5)
    q = 2 - 3 * a * x**2
    P(f"   δ = {dl:.0e}: 1 − q = {1-q:.6e},  2(3a)^(1/4)δ^(1/2) = {2*(3*a)**0.25*dl**0.5:.6e}")
P("   → Sattel-Knoten-Normalform (Standard der Verzweigungstheorie), Fenster ∝ ε^(3/2), Erholung ∝ √Abstand.")
for g in [1e-3, 1e-2, 1e-1]:
    s = math.sqrt(g)
    P(f"   Rückkehr mit Verlust Γτ = {g:.0e}: 1 − e^(−√(Γτ)) = {1-math.exp(-s):.4e}, √(Γτ) = {s:.4e} (Kabelgleichung, Warburg)")

open(__file__.replace(".py", ".txt"), "w").write("\n".join(out) + "\n")
