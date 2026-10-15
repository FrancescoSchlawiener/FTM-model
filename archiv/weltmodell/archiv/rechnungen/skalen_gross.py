# -*- coding: utf-8 -*-
"""
FTM - Fraktales Torus-Faltungs-Modell
Rechnungen fuer die grossen Skalen (Galaxie bis Kosmos) und die vier Bausteine
(Zeit, Verschraenkung, Syntropie, Gedaechtnis).

Formeln des Autors:
    g(x, y)  = x / y^2 * g_Erde          x = Masse in Erdmassen, y = Radius in Erdradien
    ds(x)    ~ G * M / (3 c^2)  ~ x      Faltungstiefe (Zusatzstrecke)
    dtau/dt  = sqrt(1 - 2 G M / (r c^2)) Zeit als lokales Faltungsdifferential

Direkt daraus folgende Umformungen:
    ds(x)    = x * ds_Erde
    2GM/c^2  = 6 * ds        ->  dtau/dt = sqrt(1 - 6 ds / r)
    Zeitdifferential  d(r) = 1 - dtau/dt  ~ 3 ds / r   (fuer ds << r)
    Grenzradius r0 = 6 ds    (dort dtau/dt = 0)
    Summe:   S = sum_i ds_i = G/(3c^2) * sum_i M_i = ds(sum_i x_i)

Ausgabe: skalen_gross.txt im selben Ordner.
Aufruf:  python skalen_gross.py
"""
import math
import os
import sys

# ---------------------------------------------------------------- Konstanten
G = 6.674e-11            # m^3 kg^-1 s^-2  (einheitliche Konstanten 01-08)
c = 2.99792458e8         # m/s
g_E = 9.81               # m/s^2  (g_Erde)
M_E = 5.972e24           # kg     (1 Erdmasse, x = 1)
R_E = 6.371e6            # m      (1 Erdradius, y = 1)
M_S = 1.989e30           # kg     Sonnenmasse
R_S = 6.957e8            # m      Sonnenradius (695 700 km)
AU = 1.495978707e11      # m
LY = 9.4607e15           # m      Lichtjahr
PC = 3.0856775814913673e16  # m   Parsec
YR = 3.15576e7           # s      Jahr
L_S = 3.828e26           # W      Leuchtkraft der Sonne

out_lines = []


def p(s=""):
    out_lines.append(s)


# ---------------------------------------------------------------- Formeln
def x_of(M):
    """Masse in Erdmassen."""
    return M / M_E


def y_of(R):
    """Radius in Erdradien."""
    return R / R_E


def g_xy(x, y):
    """Huell-Beschleunigung g(x,y) = x/y^2 * g_Erde  [m/s^2]."""
    return x / y**2 * g_E


def ds_M(M):
    """Faltungstiefe ds = G M / (3 c^2)  [m]."""
    return G * M / (3 * c**2)


DS_E = ds_M(M_E)          # Faltungstiefe der Erde (x = 1)
DS_S = ds_M(M_S)          # Faltungstiefe der Sonne


def ds_x(x):
    """Faltungstiefe in Erdmassen: ds(x) = x * ds_Erde."""
    return x * DS_E


def dtau_dt(M, r):
    """dtau/dt = sqrt(1 - 2GM/(r c^2)) ; fuer r <= 2GM/c^2 -> 0."""
    q = 2 * G * M / (r * c**2)
    if q >= 1:
        return 0.0, q
    return math.sqrt(1 - q), q


def zeitdiff(M, r):
    """1 - dtau/dt, numerisch stabil fuer sehr kleine Werte."""
    q = 2 * G * M / (r * c**2)
    if q >= 1:
        return 1.0
    return q / (1 + math.sqrt(1 - q))


def fmt(v, n=3):
    if v == 0:
        return "0"
    a = abs(v)
    if 1e-3 <= a < 1e5:
        return f"{v:.{n+1}g}"
    return f"{v:.{n}e}"


def table(header, rows):
    w = [len(h) for h in header]
    for r in rows:
        for i, cell in enumerate(r):
            w[i] = max(w[i], len(cell))
    line = " | ".join(h.ljust(w[i]) for i, h in enumerate(header))
    p(line)
    p("-+-".join("-" * w[i] for i in range(len(header))))
    for r in rows:
        p(" | ".join(cell.ljust(w[i]) for i, cell in enumerate(r)))
    p()


# ================================================================ 0. Basis
p("=" * 78)
p("FTM - Rechnungen grosse Skalen (Galaxie bis Kosmos) und vier Bausteine")
p("=" * 78)
p()
p("0. Bezugsgroessen (x = 1, y = 1)")
p("-" * 78)
p(f"g_Erde                       = {g_E} m/s^2")
p(f"ds_Erde = G M_E/(3c^2)       = {DS_E:.6e} m  = {DS_E*1e3:.4f} mm")
p(f"ds_Erde / R_Erde             = {DS_E/R_E:.6e}")
p(f"6 * ds_Erde / R_Erde         = {6*DS_E/R_E:.6e}")
p(f"ds_Sonne                     = {DS_S:.4f} m")
p(f"x_Sonne (Erdmassen)          = {x_of(M_S):.1f}")
p(f"y_Sonne (Erdradien)          = {y_of(R_S):.2f}")
p(f"Grenzradius Sonne r0 = 6 ds  = {6*DS_S:.2f} m")
p()
p("Umformung: dtau/dt = sqrt(1 - 6*ds/r) = sqrt(1 - 6*(ds_E/R_E)*x/y)")
p(f"           mit 6*ds_E/R_E = {6*DS_E/R_E:.6e}")
p()

# ================================================================ 1. Systeme
p("1. Grosse Systeme: g(x,y), ds, ds/R, dtau/dt am Rand")
p("-" * 78)
p("Eingabewerte (sichtbare Massen, Aussenradius):")
systems = [
    # name, Masse kg, Radius m, Radius-Text
    ("Erde (Bezug)", M_E, R_E, "6371 km"),
    ("Sonne (Bezug)", M_S, R_S, "695 700 km"),
    ("Milchstrasse", 6.0e10 * M_S, 50000 * LY, "50 000 Lj"),
    ("Galaxienhaufen", 1.7e14 * M_S, 2.0e6 * PC, "2 Mpc"),
    ("beobacht. Universum", 1.5e53, 46.5e9 * LY, "46,5 Mrd Lj"),
]
rows = []
for n, M, R, rt in systems:
    rows.append([n, fmt(M), f"{M/M_S:.3e}", rt, fmt(R)])
table(["System", "M [kg]", "M [M_Sonne]", "R", "R [m]"], rows)

rows = []
for n, M, R, rt in systems:
    x, y = x_of(M), y_of(R)
    g = g_xy(x, y)
    ds = ds_M(M)
    dt, q = dtau_dt(M, R)
    rows.append([n, fmt(x), fmt(y), fmt(g), fmt(g / g_E), fmt(ds), fmt(ds / R),
                 f"{dt:.12f}" if dt > 0.999 else f"{dt:.5f}", fmt(zeitdiff(M, R))])
table(["System", "x", "y", "g [m/s^2]", "g/g_E", "ds [m]", "ds/R", "dtau/dt (Rand)",
       "1-dtau/dt"], rows)

p("Faltungstiefe in anschaulichen Einheiten:")
rows = []
for n, M, R, rt in systems:
    ds = ds_M(M)
    rows.append([n, fmt(ds), fmt(ds / AU), fmt(ds / LY), fmt(ds / DS_E), fmt(6 * ds)])
table(["System", "ds [m]", "ds [AE]", "ds [Lj]", "ds/ds_Erde", "r0 = 6ds [m]"], rows)

# ================================================================ 2. Summe ds
p("2. Dunkle Materie als Summeneffekt: S = sum_i ds_i = G/(3c^2) * sum_i M_i")
p("-" * 78)
p("2a. Milchstrasse, Summe ueber Komponenten")
mw = [
    ("Sternscheibe", 4.0e10),
    ("Zentralverdickung (Bulge)", 1.0e10),
    ("Gas", 1.0e10),
    ("Zentrum (kompakte Masse)", 4.3e6),
]
S = 0.0
rows = []
for n, m in mw:
    ds = ds_M(m * M_S)
    S += ds
    rows.append([n, f"{m:.3e}", fmt(ds), fmt(ds / LY)])
rows.append(["SUMME S_MW", f"{sum(m for _, m in mw):.4e}", fmt(S), fmt(S / LY)])
table(["Komponente", "M [M_Sonne]", "ds_i [m]", "ds_i [Lj]"], rows)
S_MW = S
R_MW = 50000 * LY
p(f"Kontrolle Linearitaet: ds(sum x_i) = {ds_M(sum(m for _, m in mw)*M_S):.6e} m = S")
p(f"S_MW / R_MW (R = 50 000 Lj)          = {S_MW/R_MW:.4e}")
p()
p("Summe ueber Einzelsterne: S = N * m_mittel * ds_Sonne")
rows = []
for N, mm in [(1e11, 0.5), (2e11, 0.3), (4e11, 0.15)]:
    Sx = N * mm * DS_S
    rows.append([f"{N:.0e}", f"{mm}", fmt(Sx), fmt(Sx / LY)])
table(["N Sterne", "m_mittel [M_Sonne]", "S [m]", "S [Lj]"], rows)

p("Summe innerhalb eines Radius r (Straffung S(<r)/r), Eingabewerte eingeschl. Masse:")
prof = [
    ("Zentrum, r = 0,1 Lj", 0.1 * LY, 4.3e6),
    ("r = 1 000 Lj", 1000 * LY, 1.0e10),
    ("Sonnenbahn, r = 26 700 Lj", 26700 * LY, 5.5e10),
    ("Rand, r = 50 000 Lj", 50000 * LY, 6.0e10),
]
rows = []
for n, r, m in prof:
    Sx = ds_M(m * M_S)
    dt, q = dtau_dt(m * M_S, r)
    rows.append([n, f"{m:.2e}", fmt(Sx), fmt(Sx / r), fmt(zeitdiff(m * M_S, r))])
table(["Ort", "M(<r) [M_Sonne]", "S(<r) [m]", "S(<r)/r", "1-dtau/dt"], rows)

p("2b. Galaxienhaufen, Summe ueber Komponenten")
cl = [
    ("Sterne in ~1000 Galaxien", 2.0e13),
    ("Heisses Haufengas", 1.5e14),
]
S = 0.0
rows = []
for n, m in cl:
    ds = ds_M(m * M_S)
    S += ds
    rows.append([n, f"{m:.3e}", fmt(ds), fmt(ds / LY)])
rows.append(["SUMME S_H", f"{sum(m for _, m in cl):.4e}", fmt(S), fmt(S / LY)])
table(["Komponente", "M [M_Sonne]", "ds_i [m]", "ds_i [Lj]"], rows)
S_H = S
R_H = 2.0e6 * PC
p(f"S_H / R_H (R = 2 Mpc)   = {S_H/R_H:.4e}")
p(f"S_H / S_MW              = {S_H/S_MW:.4e}")
p(f"Zentrum (4.3e6 M_S): ds = {ds_M(4.3e6*M_S):.4e} m = {ds_M(4.3e6*M_S)/R_S:.3f} R_S, r0 = {6*ds_M(4.3e6*M_S):.4e} m")
p()
p("Netzwerk-Summe fuer eine Galaxie im Haufen: S_Netz = S_eigen + sum_j ds_j (Umgebung)")
rows = []
for n, Sx in [("nur eigene Galaxie (MW-artig)", S_MW),
              ("+ Haufengas + 999 weitere Galaxien", S_H + S_MW)]:
    rows.append([n, fmt(Sx), fmt(Sx / LY), fmt(Sx / S_MW)])
table(["Pfad", "S [m]", "S [Lj]", "S / S_MW"], rows)

# ================================================================ 3. Filamente
p("3. Filamente und Knoten: S_Knoten = sum_k S_Haufen,k ; Linienfaltung lambda = S/L")
p("-" * 78)
rows = []
for nf, nh, Lmpc in [(1, 5, 30), (3, 5, 30), (6, 5, 30)]:
    Sk = nf * nh * S_H
    L = Lmpc * 1e6 * PC
    rows.append([str(nf), str(nh), f"{Lmpc}", fmt(Sk), fmt(Sk / LY),
                 fmt(nh * S_H / L)])
table(["Filamente", "Haufen/Filament", "L [Mpc]", "S_Knoten [m]", "S_Knoten [Lj]",
       "lambda = S_Fil/L"], rows)

# ================================================================ 4. Torus-Schleife
p("4. Quasi-periodische Wicklung: theta_n = 2 pi n w (mod 2 pi), w = pi")
p("-" * 78)
w = math.pi
best = 1.0
rows = []
for n in range(1, 400000):
    d = abs(n * w - round(n * w))
    if d < best:
        best = d
        rows.append([str(n), str(round(n * w)), f"{d:.3e}",
                     f"{d*360:.3e}"])
table(["n Umlaeufe", "naechste ganze Zahl", "Restabstand [Umlauf]", "Restwinkel [Grad]"],
      rows)
p("Kein n liefert Restabstand 0 (w irrational).")
p()

# ================================================================ 5. Expansion / Rotverschiebung
p("5. Expansion (Atmen) und Rotverschiebung: 1 + z = L(t_ob)/L(t_em)")
p("-" * 78)
rows = []
for st in [0.01, 0.1, 0.5, 1.0, 2.0, 10.0]:
    rows.append([f"{st*100:.0f} %", f"{1+st:.2f}", f"{st:.2f}",
                 f"{1/(1+st):.4f}"])
table(["Dehnung des Pfades", "L_ob/L_em", "z", "E_ob/E_em"], rows)
alpha = 70e3 / (1e6 * PC)  # 1/s, Atemrate aus 70 km/s/Mpc
p(f"Atemrate alpha = (1/L) dL/dt, Eingabewert 70 km/s/Mpc = {alpha:.4e} 1/s")
p(f"   relative Pfadverlaengerung pro Jahr = {alpha*YR:.4e}")
p(f"   Zeitskala 1/alpha                     = {1/alpha/YR:.4e} Jahre")
p()
p("Faltungsanteil am Lichtpfad: L_eff = L_flach + sum ds_i")
rows = []
for n, Sx, L in [("Durchgang Milchstrasse (Durchmesser 100 000 Lj)", S_MW, 2 * R_MW),
                 ("Durchgang Haufen (Durchmesser 4 Mpc)", S_H, 2 * R_H),
                 ("Universum, R = 46,5 Mrd Lj", ds_M(1.5e53), 46.5e9 * LY)]:
    rows.append([n, fmt(Sx), fmt(L), fmt(Sx / L)])
table(["Pfad", "sum ds [m]", "L [m]", "sum ds / L"], rows)

# ================================================================ 6. Kaleidoskop
p("6. Optisches Kaleidoskop: Bilder bei d_n = d_0 + n U, Zeitversatz n U / c")
p("-" * 78)
p("Autorbeispiel: Galaxie in 10 Mrd Lj = eigene Region vor 10 Mrd Jahren -> U = 10 Mrd Lj")
rows = []
for n in range(1, 5):
    rows.append([str(n), f"{n*10:.0f}", f"{n*10:.0f}"])
table(["n", "Abstand d_n [Mrd Lj]", "Zeitversatz [Mrd Jahre]"], rows)
p("Anzahl Kopien einer Grundzelle (Kantenlaenge U) in Kugel vom Radius R_B:")
p("   k = (4/3) pi (R_B/U)^3 ;  N_real = N_beob / k ;  N_beob = 2e12 Galaxien")
rows = []
for RB in [13.8, 46.5]:
    for U in [5.0, 10.0, 20.0]:
        k = 4 / 3 * math.pi * (RB / U) ** 3
        rows.append([f"{RB}", f"{U}", f"{k:.1f}", f"{2e12/k:.3e}"])
table(["R_B [Mrd Lj]", "U [Mrd Lj]", "k", "N_real"], rows)

# ================================================================ 7. Dunkle Energie
p("7. Dunkle Energie: Ausgleich der globalen Spannung  dL_aussen = sum ds_innen")
p("-" * 78)
S_U = ds_M(1.5e53)
R_U = 46.5e9 * LY
p(f"sum ds (Universum, M = 1,5e53 kg) = {S_U:.4e} m = {S_U/LY:.4e} Lj")
p(f"relative Rueckstellung dL/R       = {S_U/R_U:.4e}")
p()

# ================================================================ 8. Zeit
p("8. Zeit: dtau/dt = sqrt(1 - 6 ds / r) fuer verschiedene Faltungstiefen")
p("-" * 78)
MSgA = 4.3e6 * M_S
r0SgA = 6 * ds_M(MSgA)
cases = [
    ("flacher Aussenraum (keine Faltung)", 0.0, 1.0),
    ("Erde, Bahn im Sonnentrichter (1 AE)", M_S, AU),
    ("Erde, Oberflaeche", M_E, R_E),
    ("Milchstrasse, Rand", 6.0e10 * M_S, 50000 * LY),
    ("Galaxienhaufen, Rand", 1.7e14 * M_S, 2.0e6 * PC),
    ("Sonne, Oberflaeche", M_S, R_S),
    ("Weisser Zwerg (1 M_S, 5800 km)", M_S, 5.8e6),
    ("Neutronenstern (1,4 M_S, 12 km)", 1.4 * M_S, 1.2e4),
    ("MW-Zentrum, r = 100 r0", MSgA, 100 * r0SgA),
    ("MW-Zentrum, r = 10 r0", MSgA, 10 * r0SgA),
    ("MW-Zentrum, r = 2 r0", MSgA, 2 * r0SgA),
    ("MW-Zentrum, r = 1,1 r0", MSgA, 1.1 * r0SgA),
    ("Universum, r = 46,5 Mrd Lj", 1.5e53, 46.5e9 * LY),
    ("Grenzradius r = r0 = 6 ds", M_S, 6 * DS_S),
]
rows = []
for n, M, r in cases:
    ds = ds_M(M)
    dt, q = dtau_dt(M, r)
    zd = zeitdiff(M, r)
    dtstr = f"{dt:.12f}" if dt > 0.9999 else f"{dt:.6f}"
    rows.append([n, fmt(ds), fmt(r), fmt(ds / r), dtstr, fmt(zd), fmt(zd * YR),
                 fmt(3 * ds / r)])
table(["Ort", "ds [m]", "r [m]", "ds/r", "dtau/dt", "1-dtau/dt", "s pro Jahr",
       "Naeherung 3 ds/r"], rows)

p("Zeitdifferenz zwischen zwei Faltungstiefen (pro Jahr Aussenzeit):")
pairs = [
    ("Erde Oberflaeche vs. flach", zeitdiff(M_E, R_E)),
    ("Sonne Oberflaeche vs. Erde Oberflaeche", zeitdiff(M_S, R_S) - zeitdiff(M_E, R_E)),
    ("Neutronenstern vs. Sonne", zeitdiff(1.4 * M_S, 1.2e4) - zeitdiff(M_S, R_S)),
]
rows = [[n, fmt(v), fmt(v * YR), fmt(v * YR / 86400)] for n, v in pairs]
table(["Vergleich", "Differential", "s / Jahr", "Tage / Jahr"], rows)

p("dtau/dt ueber x und y (Erdeinheiten): sqrt(1 - 6 (ds_E/R_E) x / y)")
rows = []
for n, M, R, rt in systems:
    x, y = x_of(M), y_of(R)
    val = 6 * DS_E / R_E * x / y
    g = g_xy(x, y)
    val2 = 6 * DS_E / R_E * (g / g_E) * y
    rows.append([n, fmt(x / y), fmt(val), fmt(val2)])
table(["System", "x/y", "6(ds_E/R_E) x/y", "6(ds_E/R_E)(g/g_E) y"], rows)

# ================================================================ 9. Verschraenkung
p("9. Verschraenkung: aeusserer Abstand d_aussen vs. innerer Abstand d_innen = 0")
p("-" * 78)
rows = []
for n, d in [("1 m", 1.0), ("1 AE", AU), ("100 Lj (Autorbeispiel)", 100 * LY),
             ("50 000 Lj", 50000 * LY)]:
    rows.append([n, fmt(d), fmt(d / c), "0", "0"])
table(["d_aussen", "[m]", "Lichtlaufzeit aussen [s]", "d_innen [m]", "Laufzeit innen [s]"],
      rows)
m_p = 1.67262192e-27
m_e = 9.1093837e-31
p("Faltungstiefen im Mikro-Faltungsbereich:")
rows = [["Proton", fmt(m_p), fmt(ds_M(m_p))], ["Elektron", fmt(m_e), fmt(ds_M(m_e))]]
table(["Teilchen", "M [kg]", "ds [m]"], rows)

# ================================================================ 10. Syntropie
p("10. Syntropie / Recycling: Massefluss aus Abstrahlung  dM/dt = L / c^2")
p("-" * 78)
Mdot_S = L_S / c**2
p(f"Sonne: L = {L_S:.3e} W  -> dM/dt = {Mdot_S:.4e} kg/s")
p(f"       d(ds)/dt = G/(3c^2) dM/dt = {ds_M(Mdot_S):.4e} m/s = {ds_M(Mdot_S)*YR:.4e} m/Jahr")
p(f"       in 10 Mrd Jahren: dM = {Mdot_S*1e10*YR:.4e} kg = {Mdot_S*1e10*YR/M_S:.4e} M_Sonne,"
  f" d(ds) = {ds_M(Mdot_S*1e10*YR):.4f} m")
L_MW = 5e36
Mdot_MW = L_MW / c**2
p(f"Milchstrasse: L = {L_MW:.1e} W -> dM/dt = {Mdot_MW:.4e} kg/s = "
  f"{Mdot_MW*YR/M_S:.4e} M_Sonne/Jahr")
p(f"       d(ds)/dt = {ds_M(Mdot_MW)*YR:.4e} m/Jahr")
p(f"       Anteil pro Jahr an S_MW = {ds_M(Mdot_MW)*YR/S_MW:.4e}")
p()

# ================================================================ 11. Gedaechtnis
p("11. Gedaechtnis des Raumes: Rillentiefe P(p) = sum_k ds_k (Knickspuren)")
p("-" * 78)
rows = []
for n, M in [("Erde", M_E), ("Sonne", M_S), ("Summe von 1000 Sonnenmassen", 1000 * M_S),
             ("Milchstrasse (sichtbar)", 6.0e10 * M_S), ("Galaxienhaufen", 1.7e14 * M_S)]:
    rows.append([n, fmt(ds_M(M)), fmt(ds_M(M) / DS_E)])
table(["Spur von", "ds [m]", "ds / ds_Erde"], rows)

# ================================================================ 12. Skalen
p("12. Selbstaehnlichkeit: g(l^2 x, l y) = g(x,y) ;  ds(l x)/(l y) = ds(x)/y")
p("-" * 78)
x0, y0 = 1.0, 1.0
rows = []
for lam in [1, 10, 1e3, 1e6]:
    rows.append([fmt(lam), fmt(g_xy(lam**2 * x0, lam * y0)),
                 fmt(ds_x(lam * x0) / (lam * y0 * R_E))])
table(["lambda", "g(l^2 x, l y) [m/s^2]", "ds(l x)/(l y R_E)"], rows)

# ---------------------------------------------------------------- Ausgabe
text = "\n".join(out_lines)
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "skalen_gross.txt"), "w", encoding="utf-8") as f:
    f.write(text + "\n")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
print(text)
