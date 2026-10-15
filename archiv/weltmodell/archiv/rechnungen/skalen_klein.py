# -*- coding: utf-8 -*-
"""
FTM - Fraktales Torus-Faltungs-Modell
Rechenskript: Skalenebenen von ganz klein bis zu den Sternen.

Formeln des Autors:
  (A1) g(x,y)  = x / y^2 * g_E          x = Masse in Erdmassen, y = Radius in Erdradien
  (A2) Ds(x)   ~ G*M / (3 c^2)  ~ x     globale Faltungstiefe ("Mehr Raum im Raum")
  (A3) dtau/dt = sqrt(1 - 2GM/(r c^2))  Zeit als lokales Faltungs-Differential

Direkt daraus folgende Beziehungen (hier ausgerechnet):
  (F1) G*M_E / R_E^2 = 9.8195 m/s^2 (SI);  Bezugswert g_E = 9.81 m/s^2 in (A1)
  (F2) Ds       = x * Ds_E ,  Ds_E = G*M_E/(3c^2)
  (F3) g * R^2  = 3 c^2 * Ds                      (Kopplung von g und Ds ueber R)
  (F4) Ds/R     = g*R/(3c^2) = (x/y) * (Ds_E/R_E)
  (F5) dtau/dt  = sqrt(1 - 6*Ds/r) = sqrt(1 - 2*G*M/(r*c^2))
  (F6) delta    = 1 - dtau/dt  ~ 3*Ds/r           (Faltungs-Differential)
  (F7) r0       = 6*Ds                            (Nullstelle von dtau/dt)

Aufruf:  python skalen_klein.py   -> schreibt skalen_klein.txt daneben
"""
import math
import os
import sys

# ---------------------------------------------------------------- Konstanten
# Einheitliche Konstanten des Weltmodells (01-08)
c = 299_792_458.0            # m/s
G = 6.674e-11                # m^3 kg^-1 s^-2
M_E = 5.972e24               # kg, Erdmasse
R_E = 6.3710e6               # m, mittlerer Erdradius
g_E = 9.81                   # m/s^2, Bezugswert g_Erde in (A1)
M_S = 1.989e30               # kg, Sonnenmasse
R_S = 6.957e8                # m, Sonnenradius (695 700 km)
GM_E = G * M_E               # m^3/s^2
GM_S = G * M_S               # m^3/s^2
AE = 1.495978707e11          # m, Abstand Erde-Sonne
JAHR = 365.25 * 86400.0      # s
TAG = 86400.0

gE_SI = GM_E / R_E**2        # (F1) G M_E / R_E^2 in SI (zum Vergleich mit g_E)
Ds_E = GM_E / (3 * c**2)     # (F2) Faltungstiefe der Erde
k_E = Ds_E / R_E             # Ds/R der Erde
x_S = M_S / M_E              # Sonnenmasse in Erdmassen
y_S = R_S / R_E              # Sonnenradius in Erdradien

# elektrische Bindung im Wasserstoff (fuer die Atom-Ebene)
ke2 = 2.307077e-28           # J*m  (Elementarladung^2 / (4 pi eps0))
m_e = 9.1093837e-31          # kg
m_p = 1.67262192e-27         # kg
a0 = 5.29177211e-11          # m, Radius des Wasserstoff-Grundzustands


# ---------------------------------------------------------------- Formeln
def g_xy(x, y):
    """(A1) Huell-Beschleunigung in m/s^2."""
    return x / y**2 * g_E


def ds_x(x):
    """(A2)/(F2) Faltungstiefe in m."""
    return x * Ds_E


def dtau_dt(x, r_m):
    """(A3)/(F5) Gangrate im Abstand r (m) von einem Zentrum der Masse x."""
    arg = 1.0 - 6.0 * ds_x(x) / r_m
    return math.sqrt(arg) if arg > 0 else float("nan")


def delta(x, r_m):
    """(F6) Faltungs-Differential 1 - dtau/dt (numerisch stabil)."""
    q = 6.0 * ds_x(x) / r_m
    if q < 1e-6:
        return q / 2 + q * q / 8 + q**3 / 16
    return 1.0 - math.sqrt(1.0 - q)


def fmt(v, n=4):
    if v != v:
        return "n/a"
    if v == 0:
        return "0"
    a = abs(v)
    if 1e-3 <= a < 1e5:
        return f"{v:.{n}g}"
    return f"{v:.{n-1}e}"


def len_fmt(m):
    """Laenge mit passender Einheit."""
    a = abs(m)
    if a >= 1e3:
        return f"{fmt(m/1e3)} km"
    if a >= 1:
        return f"{fmt(m)} m"
    if a >= 1e-3:
        return f"{fmt(m*1e3)} mm"
    return f"{fmt(m)} m"


out = []


def p(s=""):
    out.append(s)


def table(header, rows):
    p("| " + " | ".join(header) + " |")
    p("|" + "|".join("---" for _ in header) + "|")
    for r in rows:
        p("| " + " | ".join(str(v) for v in r) + " |")
    p()


# ---------------------------------------------------------------- Koerper
# (Name, Masse kg, Radius m, Ebene)
koerper = [
    ("Elektron", m_e, 2.8179403e-15, "Elementarteilchen-Schwingung"),
    ("Proton", m_p, 0.8409e-15, "Gebundene Kernmasse"),
    ("He-4-Kern", 6.6446573e-27, 1.6755e-15, "Gebundene Kernmasse"),
    ("Fe-56-Kern", 55.9349375 * 1.66053907e-27 - 26 * m_e, 3.7377e-15, "Gebundene Kernmasse"),
    ("H-Atom", m_p + m_e, a0, "Atom (Torus-Wirbel)"),
    ("Mond", 7.342e22, 1.7374e6, "Planetar"),
    ("Merkur", 3.3011e23, 2.4397e6, "Planetar"),
    ("Mars", 6.4171e23, 3.3895e6, "Planetar"),
    ("Venus", 4.8675e24, 6.0518e6, "Planetar"),
    ("Erde", M_E, R_E, "Planetar"),
    ("Uranus", 8.6810e25, 2.5362e7, "Planetar"),
    ("Neptun", 1.02413e26, 2.4622e7, "Planetar"),
    ("Saturn", 5.6834e26, 5.8232e7, "Planetar"),
    ("Jupiter", 1.89813e27, 6.9911e7, "Planetar"),
    ("Sonne", M_S, R_S, "Stern"),
    ("Sirius B (Weisser Zwerg)", 1.018 * M_S, 0.0084 * R_S, "Stern (kompakt)"),
    ("Neutronenstern", 1.4 * M_S, 1.2e4, "Stern (kompakt)"),
    ("Beteigeuze (Roter Riese)", 18.0 * M_S, 764.0 * R_S, "Stern (Roter Riese)"),
]

# ================================================================ Ausgabe
p("FTM - Rechnungen zu den Skalenebenen (klein bis Sterne)")
p("=" * 60)
p()
p("## 0. Erd-Einheiten (Bezugsgroessen)")
p()
table(["Groesse", "Wert"], [
    ["c", f"{c:.0f} m/s"],
    ["G", f"{G:.5e} m^3/(kg s^2)"],
    ["M_E (Erdmasse)", f"{M_E:.5e} kg"],
    ["R_E (Erdradius)", f"{R_E/1e3:.1f} km"],
    ["g_E (Bezugswert in A1)", f"{g_E:.2f} m/s^2"],
    ["G M_E / R_E^2 (SI)", f"{gE_SI:.4f} m/s^2"],
    ["g_E / (G M_E / R_E^2)", f"{g_E/gE_SI:.5f}"],
    ["Ds_E = G M_E / (3 c^2)", f"{Ds_E*1e3:.4f} mm"],
    ["Ds_E / R_E", f"{k_E:.4e}"],
    ["6 Ds_E / R_E = 2 G M_E / (R_E c^2)", f"{6*k_E:.4e}"],
    ["x_Sonne = M_S / M_E", f"{x_S:.1f}"],
    ["Ds / M = G / (3 c^2)", f"{G/(3*c**2):.4e} m/kg"],
    ["y_Sonne = R_S / R_E", f"{y_S:.3f}"],
    ["Ds_Sonne = x_S * Ds_E", f"{ds_x(x_S):.2f} m"],
])

p("## 1. Pruefung der Kopplung G M = 3 c^2 Ds (F3) an der Erde")
p()
p(f"G M_E             = {GM_E:.6e} m^3/s^2")
p(f"3 c^2 * Ds_E      = {3*c**2*Ds_E:.6e} m^3/s^2")
p(f"g_E * R_E^2       = {g_E*R_E**2:.6e} m^3/s^2   (mit Bezugswert g_E = 9.81; Verhaeltnis {g_E*R_E**2/GM_E:.5f})")
p(f"Gangrate Erde, (A3) direkt : {math.sqrt(1-2*GM_E/(R_E*c**2)):.15f}")
p(f"Gangrate Erde, (F5) ueber Ds: {math.sqrt(1-6*Ds_E/R_E):.15f}")
p()

p("## 2. Grosse Tabelle: alle Koerper (Oberflaeche)")
p()
rows = []
rows2 = []
for name, M, R, ebene in koerper:
    x = M / M_E
    y = R / R_E
    g = g_xy(x, y)
    ds = ds_x(x)
    ratio = ds / R
    d = delta(x, R)
    rate = dtau_dt(x, R)
    r0 = 6 * ds
    rows.append([name, fmt(x), fmt(y), fmt(g), fmt(g / g_E), fmt(ds) + " m", fmt(ratio)])
    rate_s = f"{rate:.12f}" if d > 1e-12 else "1 - " + fmt(d)
    rows2.append([name, rate_s, fmt(d), fmt(3 * ratio), fmt(d * JAHR) + " s", fmt(r0) + " m", fmt(R / r0)])
table(["Koerper", "x [M_E]", "y [R_E]", "g [m/s^2]", "g/g_E", "Ds", "Ds/R"], rows)
p("Gangrate und Faltungs-Differential an der Oberflaeche:")
p()
table(["Koerper", "dtau/dt", "delta = 1-dtau/dt", "3 Ds/R", "delta * 1 Jahr", "r0 = 6 Ds", "R / r0"], rows2)

p("## 3. Subatomare Ebene und Atom")
p()
a_el = ke2 / (m_e * a0**2)
v_el = math.sqrt(a_el * a0)
f_el = v_el / (2 * math.pi * a0)
g_H = g_xy((m_p + m_e) / M_E, a0 / R_E)     # (A1)
g_Hp = g_xy(m_p / M_E, a0 / R_E)            # (A1)
x_eq = a_el * (a0 / R_E) ** 2 / g_E
table(["Groesse", "Wert"], [
    ["Radius a0 des H-Atoms", f"{a0:.5e} m"],
    ["g aus (A1), Masse Proton+Elektron, y = a0/R_E", f"{g_H:.4e} m/s^2"],
    ["Ds des H-Atoms", f"{ds_x((m_p+m_e)/M_E):.4e} m"],
    ["Ds / a0", f"{ds_x((m_p+m_e)/M_E)/a0:.4e}"],
    ["elektro-geometrische Bindungs-Steilheit a_el", f"{a_el:.4e} m/s^2"],
    ["a_el / g_E", f"{a_el/g_E:.4e}"],
    ["a_el / g (Faltungs-Steilheit des Protons bei a0, A1)", f"{a_el/g_Hp:.4e}"],
    ["Umlaufgeschwindigkeit v = sqrt(a_el * a0)", f"{v_el:.5e} m/s"],
    ["v / c", f"{v_el/c:.5e}"],
    ["Umlauffrequenz f = v/(2 pi a0)", f"{f_el:.4e} Hz"],
    ["Umlaufzeit 1/f", f"{1/f_el:.4e} s"],
    ["aequivalentes x aus (A1) umgestellt: x = a_el y^2 / g_E", f"{x_eq:.4e} M_E = {x_eq*M_E:.4e} kg"],
    ["Ds_aequiv = x Ds_E", f"{x_eq*Ds_E:.4e} m"],
    ["Ds_aequiv / a0", f"{x_eq*Ds_E/a0:.4e}"],
    ["3 Ds_aequiv / a0", f"{3*x_eq*Ds_E/a0:.4e}"],
    ["(v/c)^2", f"{(v_el/c)**2:.4e}"],
])
p("Additivitaet von Ds in der Kernmasse (Ds ~ x):")
p()
ds_p = ds_x(m_p / M_E)
rowsA = []
for name, M, R, _ in koerper[1:4]:
    ds = ds_x(M / M_E)
    rowsA.append([name, fmt(M) + " kg", fmt(ds) + " m", fmt(ds / ds_p)])
table(["Kern", "Masse", "Ds", "Ds / Ds_Proton"], rowsA)

p("## 4. Sol-System: Planetenbahnen als dynamischer Wirbel (Neigung der Sonne am Bahnradius)")
p()
p("x = x_Sonne + x_Planet (Ds additiv), g(r) = 3 c^2 Ds(x) / r^2 = G M / r^2 (SI),")
p("v = sqrt(g * r) = c sqrt(3 Ds / r), T = 2 pi r / v")
p()
bahnen = [
    ("Merkur", 0.387098, 3.3011e23), ("Venus", 0.723332, 4.8675e24), ("Erde", 1.0000010, M_E),
    ("Mars", 1.523679, 6.4171e23), ("Jupiter", 5.2044, 1.89813e27), ("Saturn", 9.5826, 5.6834e26),
    ("Uranus", 19.2184, 8.6810e25), ("Neptun", 30.110, 1.02413e26),
]
rowsB = []
for name, a_ae, Mp in bahnen:
    r = a_ae * AE
    x = x_S + Mp / M_E
    y = r / R_E
    g = 3 * c**2 * ds_x(x) / r**2
    v = math.sqrt(g * r)
    T = 2 * math.pi * r / v
    rowsB.append([name, f"{a_ae:.4f}", fmt(y), fmt(g), fmt(v / 1e3) + " km/s",
                  fmt(T / JAHR) + " a", fmt(ds_x(x_S) / r), fmt(delta(x_S, r)), fmt((v / c)**2)])
table(["Planet", "r [AE]", "y [R_E]", "g_Sonne(r) [m/s^2]", "v", "T", "Ds_S / r", "delta_S(r)", "(v/c)^2"], rowsB)

r_m = 3.844e8
x_em = 1 + 7.342e22 / M_E
g_m = 3 * c**2 * ds_x(x_em) / r_m**2
v_m = math.sqrt(g_m * r_m)
p("Mond um die Erde (x = x_Erde + x_Mond, r = 384 400 km):")
p(f"  g = {g_m:.4e} m/s^2, v = {v_m:.1f} m/s, T = {2*math.pi*r_m/v_m/TAG:.3f} d")
p()
p("System-Summe der Faltung (Sonne + Planeten + Mond):")
x_pl = sum(Mp for _, _, Mp in bahnen) / M_E + 7.342e22 / M_E
p(f"  Sum x = {x_S:.1f} + {x_pl:.2f} = {x_S + x_pl:.2f}")
p(f"  Ds_System = {ds_x(x_S + x_pl):.3f} m;  Anteil der Sonne = {100 * x_S / (x_S + x_pl):.3f} %")
p()

p("## 5. Interplanetarer Raum: unbewegter Hintergrund um die Sonne")
p()
rowsC = []
for lab, r in [("Sonnenoberflaeche", R_S), ("10 R_S", 10 * R_S), ("0.1 AE", 0.1 * AE),
               ("1 AE", AE), ("5.2 AE", 5.2044 * AE), ("30 AE", 30.11 * AE),
               ("100 AE", 100 * AE), ("1000 AE", 1000 * AE), ("1 Lichtjahr", 9.4607e15)]:
    g = g_xy(x_S, r / R_E)
    d = delta(x_S, r)
    rowsC.append([lab, fmt(r / R_E), fmt(g), fmt(g / g_E), fmt(ds_x(x_S) / r), fmt(d), fmt(d * JAHR) + " s"])
table(["Ort", "y [R_E]", "g [m/s^2]", "g/g_E", "Ds_S/r", "delta", "delta * 1 Jahr"], rowsC)
y_eq = math.sqrt(x_S)
p(f"Radius, bei dem die Sonnenmasse die Erd-Steilheit g_E haette: y = sqrt(x_S) = {y_eq:.1f} R_E"
  f" = {y_eq*R_E/1e3:.4e} km = {y_eq*R_E/R_S:.3f} R_S")
p()

p("## 6. Planetare Kerne: Faltungstiefe Ds im Zentrum")
p()
rowsD = []
for name, M, R, ebene in koerper:
    if ebene != "Planetar":
        continue
    x = M / M_E
    ds = ds_x(x)
    rowsD.append([name, fmt(x), len_fmt(ds), fmt(ds / R), len_fmt(6 * ds)])
table(["Koerper", "x", "Ds", "Ds/R", "r0 = 6 Ds"], rowsD)
R_kern = 3.4855e6
R_ikern = 1.2215e6
p(f"Erde: Ds / R_Kern (3485.5 km)        = {Ds_E/R_kern:.4e}")
p(f"Erde: Ds / R_innerer Kern (1221.5 km) = {Ds_E/R_ikern:.4e}")
p()

p("## 7. Luftballon-Effekt: gleiche Masse, wachsender Radius")
p()
for lab, xs in [("eine Sonnenmasse", x_S), ("Beteigeuze-Masse (18 Sonnenmassen)", 18 * x_S)]:
    p(f"### {lab}: x = {fmt(xs)} M_E, Ds = {ds_x(xs):.2f} m (konstant)")
    p()
    rowsE = []
    for zlab, R in [("12 km (Neutronenstern)", 1.2e4), ("0.0084 R_S (Weisser Zwerg)", 0.0084 * R_S),
                    ("1 R_S", R_S), ("10 R_S", 10 * R_S), ("100 R_S", 100 * R_S),
                    ("200 R_S", 200 * R_S), ("764 R_S", 764 * R_S), ("1000 R_S", 1000 * R_S)]:
        if R <= 6 * ds_x(xs):
            continue          # Radius innerhalb r0 = 6 Ds: nicht in der Reihe
        y = R / R_E
        g = g_xy(xs, y)
        d = delta(xs, R)
        rate = dtau_dt(xs, R)
        rate_s = f"{rate:.6f}" if d > 1e-6 else "1 - " + fmt(d)
        rowsE.append([zlab, fmt(y), fmt(g), fmt(g / g_E), f"{ds_x(xs):.2f} m", fmt(ds_x(xs) / R), rate_s])
    table(["Radius", "y [R_E]", "g [m/s^2]", "g/g_E", "Ds", "Ds/R", "dtau/dt Oberflaeche"], rowsE)
p("Invariante des Luftballon-Effekts: g * y^2 = x * g_E (A1);  SI-Form G M = 3 c^2 Ds")
xs = x_S
for R in (R_S, 100 * R_S):
    y = R / R_E
    p(f"  R = {R/R_S:>5.0f} R_S: g*y^2 = {g_xy(xs, y)*y**2:.6e} m/s^2,  g*R^2 = {g_xy(xs, y)*R**2:.6e} m^3/s^2,"
      f"  G M = 3c^2 Ds = {3*c**2*ds_x(xs):.6e} m^3/s^2")
for lab, xs in [("1 M_S", x_S), ("18 M_S", 18 * x_S)]:
    p(f"  y fuer g = g_E bei {lab}: y = sqrt(x) = {math.sqrt(xs):.1f} R_E = {math.sqrt(xs)*R_E/R_S:.3f} R_S")
p(f"  g(12 km) / g(1000 R_S) bei 1 M_S = {(1000*R_S/1.2e4)**2:.4e}")
p(f"  Sonne / Beteigeuze (Oberflaechen-Neigung) = {g_xy(x_S, y_S)/g_xy(18*x_S, 764*R_S/R_E):.1f}")
p()

p("## 8. Skalenvergleich Ds/R (Mikro bis Makro)")
p()
rowsF = []
for name, M, R, ebene in sorted(koerper, key=lambda k: ds_x(k[1] / M_E) / k[2]):
    x = M / M_E
    rowsF.append([name, ebene, fmt(ds_x(x) / R), fmt(g_xy(x, R / R_E))])
table(["Koerper", "Ebene", "Ds/R", "g [m/s^2]"], rowsF)

p("## 9. Quasi-periodische Wicklung: Faltungsverhaeltnis rho = pi")
p()
p("Bahn auf dem Torus: theta1 = 2 pi * n * rho, theta2 = 2 pi * n  (n = Umlaeufe um theta2).")
p("Rueckkehrabstand nach n Umlaeufen: d(n) = Abstand von n*rho zur naechsten ganzen Zahl")
p("(Bruchteil einer vollen Windung). Aufgelistet: jede neue beste Annaeherung bis n = 10^6.")
p()
rho = math.pi
best = 1.0
rowsW = []
for n in range(1, 1_000_001):
    v = n * rho
    m = round(v)
    d = abs(v - m)
    if d < best:
        best = d
        rowsW.append([n, m, f"{m}/{n}", f"{d:.3e}", f"{d*360:.3e} Grad"])
table(["n (Umlaeufe)", "m (Windungen)", "Naeherung m/n", "d(n)", "Winkelabstand"], rowsW)
p("d(n) wird nie 0: die Bahn schliesst sich fuer kein n.")
p()

text = "\n".join(out)
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "skalen_klein.txt"), "w", encoding="utf-8") as fh:
    fh.write(text + "\n")
sys.stdout.reconfigure(encoding="utf-8")
print(text)
