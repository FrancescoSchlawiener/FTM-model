# -*- coding: utf-8 -*-
"""Q als Formverhältnis: aus welchen Formbeziehungen des Modells folgt ν = 1 + Q = 6,364 ± 0,065 [G]?
Nur mit ftm_kern. Ausgabe q_form.txt."""
import math, itertools, sys
from ftm_kern import *

out = open(__file__.replace('.py', '.txt'), 'w', encoding='utf-8')
def P(*a):
    print(*a); print(*a, file=out)

NU = 1 + Q_KOSM; SIG = Q_fehler()
LO, HI = NU - SIG, NU + SIG
def zs(x): return (x - NU) / SIG
P(f"ν = 1 + Q = {NU:.4f} ± {SIG:.4f} [G], q = {q_KOSM:.5f}; 1σ-Band ν ∈ [{LO:.4f}, {HI:.4f}]")
P(f"Band in q: [{1-1/LO:.5f}, {1-1/HI:.5f}]; Fenster ν 6–6,5 ⇔ q ∈ [{1-1/6:.5f}, {1-1/6.5:.5f}]")

# ---------------------------------------------------------------- 1 flacher Torus
P("\n1 Flacher Torus ℂ/(ℤ+τℤ), Nome q = e^{2πiτ}  [G, Form]")
it = im_tau_aus_nome(q_KOSM)
P(f"  Im τ(q_kosm) = {it:.5f}  (Zyklen 1 : {1/it:.2f}); Fenster ν 6–6,5: Im τ = {im_tau_aus_nome(1-1/6.5):.5f}–{im_tau_aus_nome(1-1/6):.5f}")
P(f"  Nome e^{{iπτ}} statt e^{{2πiτ}}: Im τ = {2*it:.5f}")
formen = {"τ = i (Quadrat)": 1j, "τ = e^{iπ/3} (Sechseck)": complex(0.5, math.sqrt(3)/2),
          "τ = i·√2 (Rechteck 1:√2)": 1j*math.sqrt(2), "τ = (1+i)/2": complex(0.5, 0.5),
          "S: −1/τ(q_kosm)": -1/(1j*it)}
for f3, nm in [(1.0, "Würfel"), (1.15, "sechseckig"), (1.30, "raumzentriert"), (1.41, "flächenzentriert")]:
    formen[f"Zelle {nm}, Im τ = f = f³^(1/3) = {f3**(1/3):.3f}"] = 1j * f3**(1/3)
    formen[f"Zelle {nm}, Im τ = 1/f³ = {1/f3:.3f}"] = 1j / f3
P(f"  {'Form':48s} {'|q| (2πiτ)':>12s} {'ν':>9s} {'|q| (iπτ)':>11s} {'ν':>9s}")
for nm, t in formen.items():
    q1 = math.exp(-2*math.pi*t.imag); q2 = math.exp(-math.pi*t.imag)
    P(f"  {nm:48s} {q1:12.3e} {1/(1-q1):9.4f} {q2:11.3e} {1/(1-q2):9.4f}")
P("  Drehordnungen / nächste Bilder der Zellen: 2, 4, 6 (2D), 6 Würfel, 8 raumz., 12 flächenz.; ν = 6 liegt "
  f"{-zs(6):.1f}σ neben der Messung")
P("  ⇒ keine ausgezeichnete Form gibt |q| nahe 0,84: sie geben ν = 1,000–1,07 (dünner Nome) oder ganze Zahlen.")

# ---------------------------------------------------------------- 2 eingebetteter Torus
P("\n2 Eingebetteter Torus (R, r): offener Teil = Streifen von einem Radiant  [G, Form]")
P("  Streifen im langen Winkel φ (ein Radiant von 2π): Anteil 1/2π für jedes R/r ⇒ ν = 2π = "
  f"{2*math.pi:.4f} ({zs(2*math.pi):+.2f}σ)")
P(f"  {'R/r':>8s} {'ν außen':>10s} {'ν innen':>10s}")
for Rr in [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]:
    P(f"  {Rr:8.2f} {torus_streifen_nu(Rr, 1, 'aussen'):10.4f} {torus_streifen_nu(Rr, 1, 'innen'):10.4f}")
from scipy.optimize import brentq
def Rr_innen(n): return brentq(lambda x: torus_streifen_nu(x, 1, 'innen') - n, 1.0001, 1e6)
P(f"  Horntorus R = r: außen {torus_streifen_nu(1,1,'aussen'):.4f}, innen {torus_streifen_nu(1,1,'innen'):.2f}")
P(f"  ν > 2π nur mit Streifen innen: ν = {NU:.3f} ⇔ R/r = {Rr_innen(NU):.1f}; 1σ-Band R/r = "
  f"{Rr_innen(HI):.1f}–{Rr_innen(LO):.0f} (außen nie > 2π)")
# Rohr-Radiant als Bogen r auf dem inneren/äußeren Äquator (Länge 2π(R∓r))
P("  Bogen r (ein Radiant des Rohres) auf dem Äquator: ν = 2π(R ∓ r)/r ⇒ ν = 6,364 bei R/r = "
  f"{1+NU/(2*math.pi):.4f} (innen) bzw. {NU/(2*math.pi)-1:.4f} (außen, Spindel); R/r = 2 gibt wieder 2π")

# ---------------------------------------------------------------- 3 zwei Beziehungen zusammen
P("\n3 Fenster aus Formbeziehungen  [G, Form]")
lam = schliesslaenge(); L = zelle_aus_a0(A0)
rho_b = OMEGA_B * 3 * H0_si()**2 / (8*math.pi*G)
Mb = rho_b * L**3
rs = r_stern(Mb)
P(f"  Schluss: r*(M_b,Zelle)/Λ = {rs/lam:.4f}, √(6/ν) = {math.sqrt(6/NU):.4f}, bei ν = 2π {math.sqrt(6/(2*math.pi)):.4f}  ⇒ ν ≥ 6")
P("  Zählung je Radiant: ν = 2π (Gleichung, keine Grenze)")
Kf = 3*c**4/(2*G*rho_b*(KREISPAAR_GRENZE*d_LSS())**3)
P(f"  Dritte Beziehung (Zelle ohne Kreispaare): (1+Q)·a₀ ≤ K·f³, K = {Kf:.4e}")
P(f"  {'a₀':>6s} {'Würfel':>8s} {'sechs.':>8s} {'raumz.':>8s} {'flächenz.':>9s}   (obere Grenze für ν)")
for a0 in [1.04e-10, 1.06e-10, 1.10e-10, 1.13e-10, 1.19e-10]:
    P(f"  {a0*1e10:6.2f} " + " ".join(f"{Kf*f3/a0:8.3f}" for f3 in [1, 1.15, 1.30, 1.41]))
P(f"  Kern q_fenster (a₀ = 1,1): Q ∈ [{q_fenster()[0]:.2f}, {q_fenster()[1]:.3f}] ⇒ ν ∈ [6, {q_fenster()[1]+1:.3f}]")
tf = 16*(1.5-math.sqrt(2))
P(f"  Torusfaktor 16(3/2−√2) = {tf:.4f}: Maß der Fluchtmoden (Dunkle Energie); ν/Faktor = {NU/tf:.3f}, "
  f"6·Faktor = {6*tf:.3f}, 2π/Faktor = {2*math.pi/tf:.3f}, e^Faktor = {math.exp(tf):.3f}: kein Glied für Q (Zufallsprobe 5)")
P(f"  Messung {-zs(6):.1f}σ über 6, {-zs(2*math.pi):.2f}σ über 2π, {zs(Kf/1.1e-10):.1f}σ unter der Würfelgrenze (a₀ = 1,1)")

# ---------------------------------------------------------------- 4 Selbstähnlichkeit
P("\n4 Selbstähnlichkeit (nur zeitfeste Verhältnisse)")
for nm, n in [("Messung", NU), ("ν = 2π", 2*math.pi), ("ν = 6", 6.0)]:
    qq = 1 - 1/n; lq = -math.log(qq)
    P(f"  {nm:8s}: R_Q/r* = 1/(−ln q) = {1/lq:.4f}, χ(R_Q) = (ν−1)(−ln q) = {(n-1)*lq:.4f}, f_b = 1/ν = {1/n:.4f}")
P("  Teilchen: ħ = h/2π, λ̄ = λ/2π zählen je Radiant (dieselbe 2π, keine Zahl 6,36); r₀ = 6Δs = 2GM/c² gibt die 6 in √(6/ν).")
P("  Galaxie: R_Q/r* und χ(R_Q) sind Funktionen von q, keine zweite Formbeziehung; Haufen: 1+Q = Maß des Ganzen in R₂₀₀, dieselbe Zahl.")
P("  Die Grenzen 6 (Schwarzschild-Zählung r₀/Δs) und 2π (Radiant) tauchen auf Teilchen- und Sternebene auf; 6,364 selbst nicht.")

# ---------------------------------------------------------------- 5 Zufallsprobe
P("\n5 Zufallsprobe: x∘y aus einfachen Formzahlen, ∘ ∈ +,−,×,÷")
basis = {"1": 1, "2": 2, "3": 3, "4": 4, "6": 6, "8": 8, "12": 12, "π": math.pi, "2π": 2*math.pi,
         "√π": math.sqrt(math.pi), "√2": math.sqrt(2), "√3": math.sqrt(3), "e": math.e,
         "φ": (1+math.sqrt(5))/2, "1,15": 1.15, "1,30": 1.30, "1,41": 1.41, "16(3/2−√2)": tf,
         "π²": math.pi**2, "π²/45": K_S}
werte = {}
for (a, x), (b, y) in itertools.product(basis.items(), repeat=2):
    for op, v in [("+", x+y), ("−", x-y), ("×", x*y), ("÷", x/y)]:
        if v > 0: werte.setdefault(round(v, 10), f"{a}{op}{b}")
for a, x in basis.items(): werte.setdefault(round(x, 10), a)
vals = sorted(werte)
n_rng = sum(NU/2 <= v <= 2*NU for v in vals)
erw = n_rng * math.log(HI/LO) / math.log(4)
treffer = [(v, werte[v]) for v in vals if LO <= v <= HI]
P(f"  {len(vals)} verschiedene Werte, {n_rng} in [ν/2, 2ν]; Treffer im 1σ-Band {len(treffer)}, erwartet {erw:.1f}")
for v, e_ in treffer: P(f"    {e_:>16s} = {v:.4f} ({zs(v):+.2f}σ)")
n2 = sum(abs(zs(v)) <= abs(zs(2*math.pi)) for v in vals)
P(f"  gleich nah wie 2π oder näher: {n2} Ausdrücke")

# ---------------------------------------------------------------- 6 Beobachterfehler
P("\n6 Blickrichtung je Punkt")
for z in ["1 Nome des flachen Torus: G (Zelle von außen), Form; nur Formverhältnis, kein Längenwert.",
          "2 eingebetteter Torus: G, Form; R/r ist die Einbettung, der flache Torus hat kein R/r (Einbettung = Zusatzsetzung).",
          "3 Schluss r*/Λ: r* ist L-Länge, Λ ist I (innere Richtung); ihr Verhältnis √(6/ν) ist eine Form. Kreispaar-Grenze: "
          "G mit Größe a₀ und d_LSS (Zustand: heute, mitatmend gleich) → Grenze hängt am Wert a₀, nicht rein Form.",
          "4 R_Q/r*, χ: L (am Ort der Galaxie), Form. ħ je Radiant: L, Zählweise.",
          "5 Zufallsprobe: nur Zahlen, ohne Blickrichtung."]:
    P("  " + z)

# ---------------------------------------------------------------- 7 von Hand
P("\n7 Von Hand")
P("  ln 0,84287 = ln 0,84 + ln 1,003417 = −0,174353 + 0,003411 = −0,170942; /2π = 0,027206 (Rechnung "
  f"{im_tau_aus_nome(q_KOSM):.6f}; q hier {q_KOSM:.5f})")
P("  innen: 2π/6,3643 = 0,987255; 1 − 0,987255 = 0,012745; 2 sin ½ = 0,958851; r/R = 0,013292; R/r = 75,23 "
  f"(Rechnung {Rr_innen(NU):.2f})")
P(f"  σ(Q) = 5,3643·√(0,01² + 0,006705²) = 5,3643·0,012040 = 0,06459 (Rechnung {SIG:.5f})")
out.close()
