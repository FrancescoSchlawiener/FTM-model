# -*- coding: utf-8 -*-
"""
FTM - Forschungsnotiz 11: Syntropie und die Engstelle (Punkt II.4 des Autors)
Rechnungen zu forschung/11-syntropie-engstelle.md. Nur Python-Standardbibliothek.
Aufruf:  python syntropie.py      (schreibt zusaetzlich syntropie.txt daneben)

Formeln (FTM-Groessen, GM = 3 c^2 Ds, r0 = 6 Ds = 2GM/c^2, l_P^2 = hbar G / c^3):
  S_BH / k   = A/(4 l_P^2) = 36 pi (Ds/l_P)^2           = 4 pi G M^2/(hbar c)
  T_H        = hbar c /(24 pi k Ds)                      = hbar c^3/(8 pi G M k)
  t_ev       = 138240 pi (Ds/l_P)^2 Ds/c = 3840 (S/k) Ds/c  = 5120 pi G^2 M^3/(hbar c^4)
               (Hawking-Abschaetzung nur fuer masselose Strahlung ohne Grauekoerperfaktoren)
  Smarr      : M c^2 = 2 T_H S_BH
  global     : r_H = c/H = 6 Ds_H,  Ds_H = c/(6H),  S_H/k = pi c^2/(l_P^2 H^2) = 36 pi (Ds_H/l_P)^2
               T_GH = hbar H /(2 pi k) = hbar c/(12 pi k Ds_H)
"""
import math, os, sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

_out = []
def P(s=""):
    print(s); _out.append(s)
def H(t):
    P(""); P("=" * 100); P(t); P("=" * 100)

# ---------------------------------------------------------------- Konstanten
G    = 6.674e-11            # m^3 kg^-1 s^-2   (wie im Weltmodell)
c    = 299792458.0          # m/s
hbar = 1.054571817e-34      # J s
h    = 2 * math.pi * hbar
kB   = 1.380649e-23         # J/K
m_u  = 1.66053907e-27       # kg
m_p  = 1.67262192e-27       # kg
m_e  = 9.1093837e-31        # kg
sigT = 6.6524587e-29        # m^2  Thomson
aRad = 7.5657e-16           # J m^-3 K^-4  Strahlungskonstante
M_S  = 1.989e30             # kg  (wie im Weltmodell)
R_S  = 6.957e8              # m
L_S  = 3.828e26             # W
yr   = 3.15576e7            # s
MPC  = 3.0857e22            # m
lP   = math.sqrt(hbar * G / c**3)
tP   = lP / c

# Kosmologie: Eingabewerte wie im Weltmodell (H0 = 70 km/s/Mpc), flaches LambdaCDM
H0   = 70e3 / MPC
Om, OL = 0.3, 0.7
T_CMB = 2.7255              # K
R_obs = 4.4e26              # m  Radius beobachtbares Universum (wie 03: 46,5 Mrd Lj)
M_U   = 1.5e53              # kg Masse des beobachtbaren Universums (wie 03, 05-folgerungen)
H_L   = H0 * math.sqrt(OL)  # asymptotischer de-Sitter-Wert
t0_yr = 13.8e9

def Ds(M):   return G * M / (3 * c**2)
def S_BH(M): return 36 * math.pi * (Ds(M) / lP)**2                 # in k
def T_H(M):  return hbar * c / (24 * math.pi * kB * Ds(M))
def t_ev(M): return 138240 * math.pi * (Ds(M) / lP)**2 * Ds(M) / c  # s
def S_BH_std(M): return 4 * math.pi * G * M**2 / (hbar * c)
def T_H_std(M):  return hbar * c**3 / (8 * math.pi * G * M * kB)
def t_ev_std(M): return 5120 * math.pi * G**2 * M**3 / (hbar * c**4)

def a_of_t(t):   # flaches Materie+Lambda-Universum, a(t0)=1 per Konstruktion der Formel
    return (Om / OL)**(1/3) * math.sinh(1.5 * math.sqrt(OL) * H0 * t)**(2/3)
def t_of_a(a):
    return 2 / (3 * math.sqrt(OL) * H0) * math.asinh(math.sqrt(OL / Om * a**3))
def Hubble(a): return H0 * math.sqrt(Om * a**-3 + OL)

def L(x):   # Laenge lesbar
    AU, LY = 1.496e11, 9.461e15
    if x < 1e3: return f"{x:.3g} m"
    if x < 1e9: return f"{x/1e3:.3g} km"
    if x < 1e15: return f"{x/AU:.3g} AE"
    if x < 1e19: return f"{x/LY:.3g} Lj"
    return f"{x/MPC:.3g} Mpc"

H("0) Konstanten und Eingaben")
P(f"  l_P = sqrt(hbar G/c^3) = {lP:.4e} m,   t_P = {tP:.4e} s")
P(f"  H0 = 70 km/s/Mpc = {H0:.4e} 1/s,  Omega_m = {Om},  Omega_L = {OL},  H_Lambda = H0 sqrt(Omega_L) = {H_L:.4e} 1/s")
t0 = t_of_a(1.0)
P(f"  Alter aus a(t) = (Om/OL)^(1/3) sinh^(2/3)(1.5 sqrt(OL) H0 t):  t0 = {t0/yr:.4e} a  (Kontrolle)")
P(f"  Fuer Verhaeltnisse und Raten wird das gemessene Weltalter t0 = {t0_yr:.3e} a verwendet.")
P(f"  T_CMB = {T_CMB} K,  R_obs = {R_obs:.3e} m,  M_U = {M_U:.3e} kg (FTM-Eingabe)")

# ====================================================================== 1
H("1) Lokale Engstelle: Bekenstein-Hawking-Entropie, Hawking-Temperatur, Verdampfung in Ds")
BHs = [("stellar 10 M_S", 10 * M_S),
       ("Sonne als Falte 1 M_S", M_S),
       ("Sgr A* 4,15e6 M_S", 4.15e6 * M_S),
       ("M87* 6,5e9 M_S", 6.5e9 * M_S),
       ("Milchstrasse 6e10 M_S", 6e10 * M_S),
       ("Universum 1,5e53 kg", M_U)]
P(f"  {'Falte':24s} {'Ds':>11s} {'r0=6Ds':>11s} {'S/k':>10s} {'T_H [K]':>10s} {'t_ev [a]':>10s} {'t_ev/t0':>10s}")
for n, M in BHs:
    P(f"  {n:24s} {L(Ds(M)):>11s} {L(6*Ds(M)):>11s} {S_BH(M):10.3e} {T_H(M):10.3e} {t_ev(M)/yr:10.3e} {t_ev(M)/yr/t0_yr:10.3e}")
P("")
P("  Kontrolle Ds-Form gegen Standardform (Quotienten muessen 1 sein):")
for n, M in BHs[:1] + BHs[3:4]:
    P(f"    {n:22s} S: {S_BH(M)/S_BH_std(M):.12f}   T: {T_H(M)/T_H_std(M):.12f}   t_ev: {t_ev(M)/t_ev_std(M):.12f}"
      f"   Smarr 2TS/(Mc^2): {2*T_H(M)*S_BH(M)*kB/(M*c**2):.12f}")
P("  Identitaet t_ev = 3840 (S/k) (Ds/c):  " +
  f"{t_ev(M_S)/(3840*S_BH(M_S)*Ds(M_S)/c):.12f}")
P("  Hinweis: t_ev ist Hawkings Abschaetzung (nur Photonen, ohne Grauekoerperfaktoren); Page (1976) mit realen")
P("  Teilchenarten aendert den Vorfaktor um eine Groessenordnung O(1), nicht die M^3-Skalierung.")

# --- Temperatur des Aussenraums: wann beginnt die Rueckgabe?
H("1b) Rueckgabe erst, wenn der Aussenraum kaelter ist als die Falte: T_CMB(a) = T_CMB/a < T_H")
M_eq = hbar * c**3 / (8 * math.pi * G * kB * T_CMB)
P(f"  Heute gilt T_H = T_CMB fuer M = {M_eq:.3e} kg (= {M_eq/7.342e22:.2f} Mondmassen).")
P("  Jede astrophysikalische Falte ist heute kaelter als der CMB und waechst netto (Absorption > Emission).")
T_GH = hbar * H_L / (2 * math.pi * kB)
P(f"  Untergrenze des Aussenraums (de Sitter): T_GH = hbar H_L/(2 pi k) = {T_GH:.3e} K")
P(f"  {'Falte':24s} {'a_start=T_CMB/T_H':>18s} {'t_start [a]':>12s} {'T_H/T_GH':>10s}")
for n, M in BHs:
    a_s = T_CMB / T_H(M)
    ts = t_of_a(a_s) / yr
    P(f"  {n:24s} {a_s:18.3e} {ts:12.3e} {T_H(M)/T_GH:10.3e}")
P("  -> Fuer alle Falten bis 6e10 M_S setzt die Nettoverdampfung nach ~1e11-1e12 Jahren ein (Expansion kuehlt den")
P("     Aussenraum), also verschwindend frueh gegen t_ev. Fuer die Universums-Masse gilt T_H < T_GH: keine Rueckgabe.")

# --- Materie-Entropie eines Sterns (Sackur-Tetrode, Mittelwerte)
H("1c) Vergleich: Entropie der Materie der Sonne (Sackur-Tetrode, gemittelt) gegen dieselbe Masse als Falte")
X, Y = 0.73, 0.27
n_H_per_kg  = X / m_p
n_He_per_kg = Y / (4.0026 * m_u)
n_e_per_kg  = n_H_per_kg + 2 * n_He_per_kg
N_tot = M_S * (n_H_per_kg + n_He_per_kg + n_e_per_kg)
mu = M_S / (N_tot * m_u)
# Virial fuer Polytrop n=3: |E_grav| = 1.5 GM^2/R, (3/2) N k T = |E_grav|/2
T_mean = 1.5 * G * M_S**2 / R_S / (3 * N_tot * kB)
V = 4/3 * math.pi * R_S**3
def s_ST(Npart, m, g):
    nQ = (2 * math.pi * m * kB * T_mean / h**2)**1.5
    return Npart * (math.log(g * nQ * V / Npart) + 2.5)
S_p  = s_ST(M_S * n_H_per_kg, m_p, 2)
S_he = s_ST(M_S * n_He_per_kg, 4.0026 * m_u, 1)
S_el = s_ST(M_S * n_e_per_kg, m_e, 2)
S_rad_sun = 4/3 * aRad * T_mean**3 * V / kB
S_sun = S_p + S_he + S_el + S_rad_sun
N_bar = M_S / m_u
P(f"  mu = {mu:.3f},  mittlere Temperatur (Virial, n=3) T = {T_mean:.3e} K,  mittlere Dichte = {M_S/V:.0f} kg/m^3")
P(f"  S(Protonen) = {S_p:.3e} k,  S(He) = {S_he:.3e} k,  S(Elektronen) = {S_el:.3e} k,  S(Strahlung) = {S_rad_sun:.3e} k")
P(f"  S_Sonne ~ {S_sun:.3e} k   = {S_sun/N_bar:.1f} k pro Baryon")
P(f"  S_BH(1 M_S) = {S_BH(M_S):.3e} k   ->  Faktor S_BH/S_Sonne = {S_BH(M_S)/S_sun:.3e}")
P("  -> Der Kollaps zur Falte erhoeht die Entropie um ~19 Groessenordnungen (Penrose).")

# --- CMB
H("1d) Entropie der kosmischen Hintergrundstrahlung")
s_gamma = (4 * math.pi**2 / 45) * (kB * T_CMB / (hbar * c))**3     # k pro m^3
T_nu = T_CMB * (4/11)**(1/3)
s_nu = (7/8) * 6/2 * (4 * math.pi**2 / 45) * (kB * T_nu / (hbar * c))**3
V_obs = 4/3 * math.pi * R_obs**3
n_gamma = 2 * 1.2020569 / math.pi**2 * (kB * T_CMB / (hbar * c))**3
P(f"  n_gamma = {n_gamma/1e6:.1f} /cm^3,  s_gamma = {s_gamma/1e6:.1f} k/cm^3,  s_nu = {s_nu/1e6:.1f} k/cm^3")
P(f"  V_obs = {V_obs:.3e} m^3")
P(f"  S_CMB = {s_gamma*V_obs:.3e} k,   S_Neutrinos = {s_nu*V_obs:.3e} k")
P(f"  Vergleich: eine einzige Falte von 1e6 M_S hat S = {S_BH(1e6*M_S):.3e} k  (> S_CMB)")
M_cmb_eq = math.sqrt(s_gamma * V_obs * hbar * c / (4 * math.pi * G))
P(f"  Eine Falte mit S = S_CMB haette M = {M_cmb_eq/M_S:.3e} M_S")

# --- Syntropie ausserhalb der Falte: Beispiel Erde
H("1e) Syntropie als Ordnungsbildung mit Entropieexport: Beispiel Erde im Strahlungsstrom der Sonne")
AU = 1.496e11; R_E = 6.371e6; albedo = 0.30
P_abs = L_S / (4 * math.pi * AU**2) * math.pi * R_E**2 * (1 - albedo)
T_sun_eff, T_E = 5772.0, 255.0
J_exp = 4/3 * P_abs * (1 / T_E - 1 / T_sun_eff)     # W/K, Strahlungsentropie (Faktor 4/3)
P(f"  absorbierte Leistung P = {P_abs:.3e} W,  T_ein = {T_sun_eff} K,  T_aus = {T_E} K")
P(f"  Entropieexport J = (4/3) P (1/T_aus - 1/T_ein) = {J_exp:.3e} W/K = {J_exp/kB:.3e} k/s")
P("  -> Ordnungsbildung (Syntropie) ist nur moeglich, weil eine heisse Quelle (die Falte Sonne) in einem kalten")
P("     Aussenraum steht. Die Falte erzeugt den Temperaturgradienten, der dissipative Strukturen treibt.")

# ====================================================================== 2
H("2) Informationserhalt: Page-Kurve (einfaches Modell S_BH(t) = S0 (1 - t/t_ev)^(2/3))")
P("  grob gekoernte Strahlungsentropie S_rad = beta (S0 - S_BH);  feinkoernig (unitaer): S_vN = min(S_rad, S_BH)")
for beta in (1.0, 1.4847):
    f = beta / (1 + beta)
    tpage = 1 - f**1.5
    P(f"  beta = {beta:.4f}:  Page-Zeit t_P/t_ev = 1 - (beta/(1+beta))^(3/2) = {tpage:.4f},  S_BH(t_P)/S0 = {f:.4f},"
      f"  M(t_P)/M0 = {math.sqrt(f):.4f}")
P("  (beta = 1,4847 ist Pages Verhaeltnis fuer Photonen+Gravitonen; Page 2013: t_Page = 0,5381 t_ev, 40,25 % von S0 abgegeben)")
P(f"  {'t/t_ev':>8s} {'S_BH/S0':>9s} {'S_rad/S0':>9s} {'S_vN/S0':>9s}")
beta = 1.4847
for x in (0, 0.2, 0.4, 0.5, 0.5381, 0.6, 0.8, 0.95, 1.0):
    sb = (1 - x)**(2/3); sr = beta * (1 - sb)
    P(f"  {x:8.4f} {sb:9.4f} {sr:9.4f} {min(sb, sr):9.4f}")
P("  Page-Zeiten konkret (beta = 1,4847):")
for n, M in BHs[:5]:
    P(f"    {n:24s} t_Page = {0.5381*t_ev(M)/yr:.3e} a")

# ====================================================================== 3
H("3) Globale Engstelle: Horizont als Schliessungsregion")
def S_hor(R): return math.pi * R**2 / lP**2
R_H0 = c / H0
R_L = c / H_L
# kosmologischer Ereignishorizont heute: R_CEH = c * int_0^1 du / H(1/u)
N = 200000
acc = 0.0
for i in range(N):
    u = (i + 0.5) / N
    acc += 1.0 / (H0 * math.sqrt(Om * u**3 + OL))
R_CEH = c * acc / N
P(f"  Hubble-Radius heute  c/H0      = {R_H0:.4e} m = {R_H0/MPC/1e3:.3f} Gpc  -> S = pi R^2/l_P^2 = {S_hor(R_H0):.3e} k")
P(f"  Ereignishorizont heute R_CEH   = {R_CEH:.4e} m = {R_CEH/9.461e15/1e9:.2f} Mrd Lj  -> S_CEH = {S_hor(R_CEH):.3e} k")
P(f"  de-Sitter-Endzustand c/H_L     = {R_L:.4e} m = {R_L/MPC/1e3:.3f} Gpc  -> S_dS  = {S_hor(R_L):.3e} k   (Maximum)")
P(f"  Gibbons-Hawking-Temperatur T_GH = {T_GH:.3e} K")
P("")
P("  Brueckenidentitaet: Masse in der Hubble-Kugel bei kritischer Dichte")
M_H = c**3 / (2 * G * H0)
M_Hsph = 4/3 * math.pi * (3 * H0**2 / (8 * math.pi * G)) * R_H0**3
P(f"    M_H = c^3/(2 G H0) = {M_H:.4e} kg ;  (4pi/3) rho_c (c/H0)^3 = {M_Hsph:.4e} kg  (identisch)")
P(f"    Ds_H = G M_H/(3c^2) = c/(6 H0) = {Ds(M_H):.4e} m ;  6 Ds_H = {6*Ds(M_H):.4e} m = c/H0 = {R_H0:.4e} m")
P(f"    S_Falte(M_H) = 36 pi (Ds_H/l_P)^2 = {S_BH(M_H):.4e} k  =  S_Hubble = {S_hor(R_H0):.4e} k  (Quotient {S_BH(M_H)/S_hor(R_H0):.12f})")
P(f"    Temperatur: Falte T_H(M_H) = hbar H0/(4 pi k) = {T_H(M_H):.3e} K ;  Horizont hbar H0/(2 pi k) = {hbar*H0/(2*math.pi*kB):.3e} K  (Faktor 2)")
P(f"    FTM-Universum: M_U/M_H = {M_U/M_H:.3f},  r0(M_U) = {6*Ds(M_U):.3e} m = {6*Ds(M_U)/R_H0:.3f} c/H0,  S(M_U)/S_dS = {S_BH(M_U)/S_hor(R_L):.3f}")
P("")
P("  Groesste Falte in de Sitter (Nariai): M_N = c^3/(3 sqrt3 G H_L), r = c/(sqrt3 H_L)")
M_N = c**3 / (3 * math.sqrt(3) * G * H_L)
P(f"    M_N = {M_N:.3e} kg = {M_N/M_S:.3e} M_S ;  S_Nariai (ein Horizont) = {S_hor(R_L/math.sqrt(3)):.3e} k = S_dS/3 ;"
  f"  beide Horizonte 2/3 S_dS")
P("")
P("  Entropiehaushalt (Egan & Lineweaver 2010): S_obs = 3,1e104 k (SMBH dominieren), S_CEH = 2,6e122 k")
S_obs = 3.1e104
P(f"    Fuellgrad heute S_obs/S_dS = {S_obs/S_hor(R_L):.2e}")
S_all_matter = S_BH(Om * M_H)
P(f"    Obergrenze Sum S_Falten <= 4 pi G (Sum M)^2/(hbar c) mit Sum M = Omega_m M_H: {S_all_matter:.3e} k = Om^2 S_H")
P(f"    aequivalent: 1e11 Falten mit S_obs zusammen haetten sqrt(<M^2>) = {math.sqrt(S_obs/1e11/S_BH(M_S)):.2e} M_S")
P(f"    Beitrag einer M87*-Falte: {S_BH(6.5e9*M_S)/S_obs:.2e} von S_obs;  Sgr A*: {S_BH(4.15e6*M_S)/S_obs:.2e}")
P(f"    Hierarchie: S_Sonne(Materie) {S_sun:.1e} << S_CMB {s_gamma*V_obs:.1e} << S_obs {S_obs:.1e} << S_CEH {S_hor(R_CEH):.1e} < S_dS {S_hor(R_L):.1e}")

# ====================================================================== 4
H("4) Kreislauf, Zeitskalen, zweiter Hauptsatz")
# Zeitskalen
n_cloud = 1e4 * 1e6; rho_cloud = 2.33 * m_u * n_cloud
t_ff = math.sqrt(3 * math.pi / (32 * G * rho_cloud))
t_KH = G * M_S**2 / (R_S * L_S)
t_Edd = sigT * c / (4 * math.pi * G * m_p)
eps = 0.1
t_Sal = eps / (1 - eps) * t_Edd
P("  Zeitskalen des Kreislaufs")
P(f"    Verdichtung: Freifallzeit Molekuelwolkenkern (n = 1e4 /cm^3)      t_ff  = {t_ff/yr:.3e} a")
P(f"    Verdichtung: Kelvin-Helmholtz-Zeit Sonne GM^2/(R L)            t_KH  = {t_KH/yr:.3e} a")
P(f"    Akkretion:   Salpeter-Zeit (eps = 0,1)                          t_Sal = {t_Sal/yr:.3e} a")
P(f"    Rueckgabe:   Verdampfung 1 M_S ... 6e10 M_S                     t_ev  = {t_ev(M_S)/yr:.2e} ... {t_ev(6e10*M_S)/yr:.2e} a")
P(f"    Rueckgabe:   Verdampfung 1e11 M_S                                t_ev  = {t_ev(1e11*M_S)/yr:.2e} a")
P("")
# Omega-Rate
A_H = hbar * c**4 / (15360 * math.pi * G**2)
P(f"  Omega physikalisch: dM/dt = -A/M^2, A = hbar c^4/(15360 pi G^2) = {A_H:.3e} kg^3/s")
P("  relative Rueckgaberate |dM/dt|/M = A/M^3 = 1/(3 t_ev):  faellt mit der Zentrierung wie M^-3")
for n, M in BHs[:5]:
    P(f"    {n:24s} |dM/dt|/M = {A_H/M**3*yr:.3e} 1/a   Leistung L_H = {A_H/M**2*c**2:.3e} W")
P("")
# Rueckverdichtung
P("  Rueckverdichtung der Strahlung zu einer neuen Falte (zweiter Hauptsatz):")
P("    N Falten der Masse M strahlen S_rad = N beta S0 ab (Energie N M c^2). Bis zur Rueckverdichtung wird die")
P("    Strahlung um den Faktor z_f = a2/a1 rotverschoben (Energie / z_f, Entropie bleibt). Eine neue Falte aus")
P("    dieser Energie hat S = (N/z_f)^2 S0. Erlaubt ist das nur, wenn (N/z_f)^2 S0 >= N beta S0, also")
P("        N >= beta * z_f^2")
for zf in (1, 10, 1e3):
    P(f"    z_f = {zf:g}:  N >= {1.4847*zf**2:.4g}")
P("    -> Rueckverdichtung ist nur hierarchisch aufwaerts erlaubt (mehr Masse in weniger Falten). In einer")
P("       Kontraktionsphase (z_f < 1, Blauverschiebung) wird sie leichter: das ist das Tolman-Problem zyklischer Modelle.")
P("")
# Kasten-Gleichgewicht
P("  Falte im geschlossenen Kasten (Hawking 1976): stabiles Gleichgewicht mit eigener Strahlung nur fuer E_rad < M c^2/4")
for n, M in BHs[1:4]:
    Vmax = M * c**2 / (4 * aRad * T_H(M)**4)
    P(f"    {n:24s} V_max = {Vmax:.3e} m^3,  R_max = {(3*Vmax/(4*math.pi))**(1/3):.3e} m  (R_obs = {R_obs:.1e} m)")
P("    -> Ohne Expansion (fester Kasten) stellt sich ein stationaeres Gleichgewicht Falte + Strahlung ein: das ist")
P("       maximale Entropie, kein Kreislauf. Erst das Atmen (V waechst, T_aussen faellt) macht das Gleichgewicht")
P("       instabil (negative Waermekapazitaet der Falte) und startet die Rueckgabe.")
C_BH = -8 * math.pi * G * M_S**2 * kB / (hbar * c)
P(f"    Waermekapazitaet der Falte C = dE/dT = -8 pi G M^2 k/(hbar c) = {C_BH:.3e} J/K fuer 1 M_S  (negativ)")
P("")
# Stationaritaetsbedingung
P("  Stationaritaet: dS_sys/dt = sigma - J = 0 und Senke mit Kapazitaet S_max = pi k c^2/(l_P^2 H^2)")
P("    dauerhaft nur, wenn dS_max/dt >= sigma_ges   <=>   -dH/dt >= l_P^2 H^3 sigma_ges /(2 pi c^2 k)")
sigma_now = S_obs / (t0_yr * yr)          # k/s, mittlere Produktionsrate (dominiert durch SMBH-Wachstum)
Hdot_now = -1.5 * H0**2 * Om
need = lP**2 * H0**3 * sigma_now / (2 * math.pi * c**2)
P(f"    sigma_ges heute ~ S_obs/t0 = {sigma_now:.3e} k/s")
P(f"    -dH/dt heute (LambdaCDM) = 1,5 H0^2 Om = {-Hdot_now:.3e} 1/s^2 ;  benoetigt: {need:.3e} 1/s^2 ;  Reserve {(-Hdot_now)/need:.2e}")
# in LambdaCDM: -Hdot = 1.5 H0^2 Om a^-3 ; need ~ lP^2 H_L^3 sigma/(2 pi c^2) (sigma konstant angenommen)
need_late = lP**2 * H_L**3 * sigma_now / (2 * math.pi * c**2)
a_fail = (1.5 * H0**2 * Om / need_late)**(1/3)
P(f"    bei konstantem sigma versagt die Bedingung in LambdaCDM ab a = {a_fail:.3e}, t = {t_of_a(a_fail)/yr:.3e} a")
t_fill = (S_hor(R_L) - S_obs) / sigma_now
P(f"    Fuellzeit der Senke bei heutigem sigma: (S_dS - S_obs)/sigma = {t_fill/yr:.3e} a")
P("    -> Mit Lambda = konst. endet jeder Kreislauf im de-Sitter-Gleichgewicht S -> S_dS. Ein dauerhaft")
P("       stationaerer Kreislauf verlangt, dass die globale Engstelle (Horizontflaeche) mitwaechst: H muss weiter fallen.")

# ====================================================================== 5
H("5) Lokal und global: dieselbe Formel, zwei Grenzfaelle")
P("  lokal : r0 = 6 Ds,        S = 36 pi k (Ds/l_P)^2,  T = hbar c/(24 pi k Ds),  kappa = c^2/(2 r0)")
P("  global: r_H = c/H = 6 Ds_H, S = 36 pi k (Ds_H/l_P)^2, T = hbar c/(12 pi k Ds_H), kappa = c^2/r_H")
for n, M in BHs:
    P(f"    {n:24s} S/S_dS = {S_BH(M)/S_hor(R_L):.3e}   r0/(c/H_L) = {6*Ds(M)/R_L:.3e}")

# ====================================================================== Einheiten
H("Einheiten-Kontrolle")
P("  [Ds/l_P] = 1 ; [S] = k ; [hbar c/(k Ds)] = J m /(J/K m) = K ; [Ds^3/(l_P^2 c)] = m/(m/s) = s")
P("  [A] = J s m^4 s^-4 /(m^6 kg^-2 s^-4) = kg^3/s ; [l_P^2 H^3 sigma/c^2] = m^2 s^-3 s^-1 /(m^2 s^-2) = s^-2")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "syntropie.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(_out) + "\n")
