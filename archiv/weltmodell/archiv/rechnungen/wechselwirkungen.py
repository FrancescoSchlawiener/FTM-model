# -*- coding: utf-8 -*-
"""
Wechselwirkungen im Zustandsraum  (Forschungsnotiz 03)
======================================================

Minimalsystem des Zustandsraums  z = (q, I, x, v, T, D, A, g):

  EM-Oszillator     q, I      (LC-Kreis, aktives Element a I - b I^3 aus DC-Quelle,
                               Joule R(T,A) I^2, Abstrahlung R_rad I^2, Seebeck S_b (T - T_a))
  Mech. Oszillator  x, v      (Feder k(T,D,g), Strömungs-Selbsterregung a_m v - b_m v^3,
                               innere Dämpfung c v^2, Schallabstrahlung c_ac v^2)
  Kopplung EM<->Mech          Gyrator Γ (Lorentzkraft Γ I, Bewegungs-EMK Γ v)
  Temperatur        T         (Wärmekapazität C_th, Kühlung G (T - T_a), Peltier, thermoelast. Term)
  Schädigung        D         (Ermüdung, dD/dt >= 0)
  Alterung          A         (Arrhenius, dA/dt >= 0)
  Gravitation       g         (langsamer Parameter, Vorspannung m g / l)

Freie Energie   Ψ = q²/2C + L I²/2 + m v²/2 + k(T,D,g) x²/2 + ψ_th(T)
                k(T,D,g) = k_el (1-D) - κ_T (T - T0) + m g / l
Innere Energie  E = q²/2C + L I²/2 + m v²/2 + K(D,g) x²/2 + C_th T,   K = k - T ∂k/∂T
Entropie        S = C_th ln T + κ_T x²/2

Zerlegung  ż = F_tan + F_trans + ε F_par,  F_trans = F_irr + F_port
  F_tan  : reversibel     (∇E·F_tan = 0, ∇S·F_tan = 0, div F_tan = 0 in (q,I,x,v,σ=C_th ln T))
  F_irr  : innere Dissipation (∇E·F_irr = 0, ∇S·F_irr >= 0)
  F_port : offene Tore    (Quellen, Abstrahlung, Kühlung, Seebeck/Peltier)
  F_par  : langsame Parameter D, A, g

Teile
  1  Identitäten der Zerlegung an Zufallszuständen
  2  Energiebilanz und Entropieproduktion entlang einer Trajektorie
  3  Jacobi-Struktur (Rotation, Kontraktion, Langsamkeit)
  4  Onsager / Casimir
  5  Schleifenverstärkungen und Σ qⁿ
  6  Torus-Attraktor: Lyapunov-Signatur, Rotationszahl, Poincaré-Schnitt
  7  Kollaps T² -> T¹ (Rast-Scan in g, langsame Rampe auf/ab, Alterung)
  8  Gravitation: FTM-Größen, Zeitdilatation, Vorspannungskanal
  9  Fluktuations-Dissipations-Beziehung
Ausgabe zusätzlich nach wechselwirkungen.txt (gleicher Ordner).
"""
import io
import os
import time

import numpy as np

_BUF = io.StringIO()


def out(*args):
    text = " ".join(str(a) for a in args)
    print(text)
    _BUF.write(text + "\n")


def hr(title):
    out("")
    out("=" * 78)
    out(title)
    out("=" * 78)


# ---------------------------------------------------------------------------
# Parameter (dimensionslos: L = C = m = l = 1, ω_EM ≈ 1)
# ---------------------------------------------------------------------------
P = dict(L=1.0, C=1.0, m=1.0, ell=1.0,
         a=0.35, b=0.2, R0=0.10, alphaR=0.5, Rrad=0.05, betaA=3.0,
         kel=0.2, kappaT=0.05, c=0.02, cac=0.01, am=0.10, bm=0.2,
         Gam=0.03,
         Cth=20.0, Gth=1.0, Ta=1.0, T0=1.0, Sb=0.05,
         epsD=0.0, Yref=0.1, mD=2.0, nu0=0.0, Ea=10.0, epsg=0.0)

IQ, II, IX, IV, IT, ID, IA, IG = range(8)


def R_of(T, A, p):
    return p['R0'] * (1 + p['alphaR'] * (T - p['T0'])) * (1 + p['betaA'] * A)


def k_of(T, D, g, p):
    return p['kel'] * (1 - D) - p['kappaT'] * (T - p['T0']) + p['m'] * g / p['ell']


def K_of(D, g, p):
    # K = k - T dk/dT  (temperaturunabhängig)
    return p['kel'] * (1 - D) + p['kappaT'] * p['T0'] + p['m'] * g / p['ell']


def slow_rates(z, p):
    x, T, A = z[IX], z[IT], z[IA]
    Y = 0.5 * p['kel'] * x * x                       # Energiefreisetzungsrate -∂Ψ/∂D
    Dd = p['epsD'] * (Y / p['Yref']) ** p['mD']      # Ermüdung, >= 0
    Ad = p['nu0'] * np.exp(-p['Ea'] / T) * (1 - A)   # Arrhenius-Alterung, >= 0
    gd = p['epsg'] * np.ones_like(x)                 # langsame Gravitationsrampe
    return Y, Dd, Ad, gd


def F_tan(z, p):
    q, I, x, v, T, D, A, g = z
    k = k_of(T, D, g, p)
    f = np.zeros_like(z)
    f[IQ] = I
    f[II] = (-q / p['C'] - p['Gam'] * v) / p['L']
    f[IX] = v
    f[IV] = (-k * x + p['Gam'] * I) / p['m']
    f[IT] = -p['kappaT'] * T * x * v / p['Cth']      # reversibler thermoelast. (Kelvin-)Term
    return f


def F_irr(z, p):
    q, I, x, v, T, D, A, g = z
    R = R_of(T, A, p)
    f = np.zeros_like(z)
    f[II] = (-R * I - p['b'] * I ** 3) / p['L']
    f[IV] = (-p['c'] * v - p['bm'] * v ** 3) / p['m']
    f[IT] = (R * I * I + p['b'] * I ** 4 + p['c'] * v * v + p['bm'] * v ** 4) / p['Cth']
    return f


def F_port(z, p):
    q, I, x, v, T, D, A, g = z
    f = np.zeros_like(z)
    f[II] = ((p['a'] - p['Rrad']) * I + p['Sb'] * (T - p['Ta'])) / p['L']
    f[IV] = (p['am'] - p['cac']) * v / p['m']
    f[IT] = (-p['Gth'] * (T - p['Ta']) - p['Sb'] * T * I) / p['Cth']
    return f


def F_par(z, p):
    Y, Dd, Ad, gd = slow_rates(z, p)
    f = np.zeros_like(z)
    f[ID] = Dd
    f[IA] = Ad
    f[IG] = gd
    f[IT] = Y * Dd / p['Cth']                        # Schädigungsdissipation Y·Ḋ -> Wärme
    return f


def F(z, p):
    return F_tan(z, p) + F_irr(z, p) + F_port(z, p) + F_par(z, p)


def energy(z, p):
    q, I, x, v, T, D, A, g = z
    return (0.5 * q * q / p['C'] + 0.5 * p['L'] * I * I + 0.5 * p['m'] * v * v
            + 0.5 * K_of(D, g, p) * x * x + p['Cth'] * T)


def entropy(z, p):
    return p['Cth'] * np.log(z[IT]) + 0.5 * p['kappaT'] * z[IX] ** 2


def gradE(z, p):
    q, I, x, v, T, D, A, g = z
    return np.array([q / p['C'], p['L'] * I, K_of(D, g, p) * x, p['m'] * v, p['Cth'] * np.ones_like(T),
                     -0.5 * p['kel'] * x * x, 0 * x, 0.5 * p['m'] * x * x / p['ell']])


def gradS(z, p):
    q, I, x, v, T, D, A, g = z
    zero = 0 * x
    return np.array([zero, zero, p['kappaT'] * x, zero, p['Cth'] / T, zero, zero, zero])


def flows(z, p):
    """Alle Leistungsflüsse (Energie pro Zeit) des Minimalsystems."""
    q, I, x, v, T, D, A, g = z
    R = R_of(T, A, p)
    Y, Dd, Ad, gd = slow_rates(z, p)
    d = {}
    d['src_EM'] = p['a'] * I * I                  # DC-Quelle über aktives Element
    d['src_flow'] = p['am'] * v * v               # Strömung (Selbsterregung)
    d['grav'] = 0.5 * p['m'] * gd * x * x / p['ell']   # parametrische Arbeit der g-Änderung
    d['joule'] = R * I * I
    d['sat_EM'] = p['b'] * I ** 4
    d['fric'] = p['c'] * v * v
    d['sat_m'] = p['bm'] * v ** 4
    d['damage'] = Y * Dd
    d['rad_EM'] = p['Rrad'] * I * I
    d['ac'] = p['cac'] * v * v
    d['cool'] = p['Gth'] * (T - p['Ta'])
    d['pelt_env'] = p['Sb'] * p['Ta'] * I         # Peltier-Wärme an kalter Lötstelle (Umgebung)
    d['seebeck'] = p['Sb'] * (T - p['Ta']) * I    # el. Leistung der Seebeck-EMK in den Kreis
    d['peltier_node'] = p['Sb'] * T * I           # dem Wärmeknoten entzogene Peltier-Wärme
    d['lorentz'] = p['Gam'] * I * v               # EM -> Mech (Lorentzkraft · v)
    d['thermoel'] = p['kappaT'] * T * x * v       # reversibler Austausch Mech <-> Wärme
    d['P_in'] = d['src_EM'] + d['src_flow'] + d['grav']
    d['P_out'] = d['rad_EM'] + d['ac'] + d['cool'] + d['pelt_env']
    d['P_diss'] = d['joule'] + d['sat_EM'] + d['fric'] + d['sat_m'] + d['damage']
    # Entropieproduktion (System + Umgebung bei T_a); Arbeitsquellen tragen keine Entropie
    d['Pi'] = (d['P_diss'] / T + p['Gth'] * (T - p['Ta']) ** 2 / (T * p['Ta'])
               + (d['rad_EM'] + d['ac']) / p['Ta'])
    d['S_flux_out'] = (d['cool'] + d['pelt_env'] + d['rad_EM'] + d['ac']) / p['Ta']
    return d


def rk4(fun, z, dt):
    k1 = fun(z)
    k2 = fun(z + 0.5 * dt * k1)
    k3 = fun(z + 0.5 * dt * k2)
    k4 = fun(z + dt * k3)
    return z + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def phases(z, p):
    q, I, x, v, T, D, A, g = z
    wm = np.sqrt(np.maximum(k_of(T, D, g, p) / p['m'], 1e-9))
    we = 1.0 / np.sqrt(p['L'] * p['C'])
    return np.arctan2(-I / we, q), np.arctan2(-v / wm, x)


t_start = time.time()
out("Wechselwirkungen im Zustandsraum – Rechnung zur Forschungsnotiz 03")
out("Parameter:", ", ".join(f"{k}={v}" for k, v in P.items()))

# ===========================================================================
hr("TEIL 1  Exakte Zerlegung  ż = F_tan + F_irr + F_port + ε F_par  (Identitäten)")
# ===========================================================================
rng = np.random.default_rng(1)
p1 = dict(P, epsD=1e-3, nu0=0.05, epsg=1e-3)
Nz = 2000
Z = np.vstack([rng.normal(0, 1, Nz), rng.normal(0, 1, Nz), rng.normal(0, 1, Nz), rng.normal(0, 1, Nz),
               rng.uniform(0.5, 2.0, Nz), rng.uniform(0, 0.5, Nz), rng.uniform(0, 0.5, Nz),
               rng.uniform(0.1, 1.0, Nz)])
gE, gS = gradE(Z, p1), gradS(Z, p1)
fl = flows(Z, p1)
dE_tan = np.sum(gE * F_tan(Z, p1), axis=0)
dS_tan = np.sum(gS * F_tan(Z, p1), axis=0)
dE_irr = np.sum(gE * F_irr(Z, p1), axis=0)
dS_irr = np.sum(gS * F_irr(Z, p1), axis=0)
dE_port = np.sum(gE * F_port(Z, p1), axis=0)
dE_par = np.sum(gE * F_par(Z, p1), axis=0)
dS_par = np.sum(gS * F_par(Z, p1), axis=0)
out(f"Stichprobe: {Nz} Zufallszustände")
out(f"max |∇E·F_tan|                         = {np.max(np.abs(dE_tan)):.3e}   (tangential: energieneutral)")
out(f"max |∇S·F_tan|                         = {np.max(np.abs(dS_tan)):.3e}   (tangential: entropieneutral)")
out(f"max |∇E·F_irr|                         = {np.max(np.abs(dE_irr)):.3e}   (innere Dissipation verschiebt nur)")
out(f"min  ∇S·F_irr                          = {np.min(dS_irr):.3e}   (>= 0)")
out(f"max |∇E·F_port - (P_in-P_grav-P_out)|  = {np.max(np.abs(dE_port - (fl['P_in'] - fl['grav'] - fl['P_out']))):.3e}")
out(f"max |∇E·F_par - P_grav|                = {np.max(np.abs(dE_par - fl['grav'])):.3e}   (Schädigung: E-neutral)")
out(f"min  ∇S·F_par                          = {np.min(dS_par):.3e}   (>= 0, Y·Ḋ/T)")


# Divergenz von F_tan in Koordinaten (q, I, x, v, σ = C_th ln T): numerisch
def div_tan_sigma(z, p, h=1e-6):
    # Koordinatenwechsel T -> σ ; dσ/dt = C_th Ṫ / T
    def G(y):
        zz = y.copy()
        zz[IT] = np.exp(y[IT] / p['Cth'])
        f = F_tan(zz, p)
        f[IT] = p['Cth'] * f[IT] / zz[IT]
        return f
    y = z.copy()
    y[IT] = p['Cth'] * np.log(z[IT])
    s = 0.0
    for i in range(5):
        e = np.zeros_like(y)
        e[i] = h
        s = s + (G(y + e)[i] - G(y - e)[i]) / (2 * h)
    return s


dv = div_tan_sigma(Z, p1)
out(f"max |div F_tan| in (q,I,x,v,σ=C_th ln T) = {np.max(np.abs(dv)):.3e}   (Liouville: tangential kontrahiert nie)")

# ===========================================================================
hr("TEIL 2  Energiebilanz und Entropieproduktion entlang einer Trajektorie")
# ===========================================================================
p2 = dict(P, epsD=1e-5, nu0=0.05, epsg=1e-4)
keys = ['src_EM', 'src_flow', 'grav', 'joule', 'sat_EM', 'fric', 'sat_m', 'damage', 'rad_EM', 'ac',
        'cool', 'pelt_env', 'seebeck', 'peltier_node', 'lorentz', 'thermoel']
nk = len(keys)


def rhs_aug(y):
    z = y[:8]
    d = flows(z, p2)
    dy = np.empty_like(y)
    dy[:8] = F(z, p2)
    dy[8] = d['P_in'] - d['P_out']         # ∫ Netto-Zufluss
    dy[9] = d['Pi']                        # ∫ Entropieproduktion
    dy[10] = d['S_flux_out']               # ∫ Entropiefluss an die Umgebung
    dy[11] = d['P_in'] + d['P_out']        # ∫ Durchsatz (Normierung)
    for j, kk in enumerate(keys):
        dy[12 + j] = d[kk]
    return dy


y = np.zeros(12 + nk)
y[:8] = [0.5, 0.0, 0.3, 0.0, 1.0, 0.0, 0.0, 0.5]
dt2, T2 = 0.01, 1000.0
n2 = int(T2 / dt2)
E0, S0 = energy(y[:8], p2), entropy(y[:8], p2)
maxerrE, maxerrS, minPi = 0.0, 0.0, np.inf
for n in range(n2):
    y = rk4(rhs_aug, y, dt2)
    if n % 10 == 0:
        z = y[:8]
        eE = abs(energy(z, p2) - E0 - y[8])
        eS = abs(entropy(z, p2) - S0 - (y[9] - y[10]))
        maxerrE, maxerrS = max(maxerrE, eE), max(maxerrS, eS)
        minPi = min(minPi, flows(z, p2)['Pi'])
zf = y[:8]
out(f"Integration: RK4, dt = {dt2}, t = {T2:.0f}, langsame Terme aktiv (ε_D={p2['epsD']}, ν0={p2['nu0']}, ġ={p2['epsg']})")
out(f"Endzustand: T = {zf[IT]:.4f}, D = {zf[ID]:.3e}, A = {zf[IA]:.3e}, g = {zf[IG]:.4f}")
out(f"Energiebilanz  E(t) - E(0) - ∫(P_in - P_out)dt : max. Fehler = {maxerrE:.3e}")
out(f"   Durchsatz ∫(P_in + P_out)dt = {y[11]:.3f}  -> relativer Fehler = {maxerrE / y[11]:.3e}")
out(f"Entropiebilanz S(t) - S(0) - ∫(Π - Φ_S)dt       : max. Fehler = {maxerrS:.3e}")
out(f"Entropieproduktion: min Π(t) = {minPi:.4e}  (>= 0),  ∫Π dt = {y[9]:.4f}")
out("Zeitmittel der Leistungsflüsse über t = 0..1000 (Vorzeichen: >0 in Pfeilrichtung):")
names = {'src_EM': 'DC-Quelle -> EM (aktives Element)', 'src_flow': 'Strömung -> Mech (Selbsterregung)',
         'grav': 'Gravitation -> Mech (parametrisch ½ṁg x²/l)', 'joule': 'EM -> Wärme (Joule R I²)',
         'sat_EM': 'EM -> Wärme (Sättigung b I⁴)', 'fric': 'Mech -> Wärme (innere Dämpfung c v²)',
         'sat_m': 'Mech -> Wärme (Sättigung b_m v⁴)', 'damage': 'Mech -> Wärme über Schädigung (Y Ḋ)',
         'rad_EM': 'EM -> außen (Abstrahlung)', 'ac': 'Mech -> außen (Schall)', 'cool': 'Wärme -> Umgebung',
         'pelt_env': 'Peltier -> Umgebung (S T_a I)', 'seebeck': 'Wärme -> EM (Seebeck-EMK · I)',
         'peltier_node': 'Wärmeknoten -> Peltier (S T I)', 'lorentz': 'EM -> Mech (Lorentz Γ I v)',
         'thermoel': 'Mech -> Wärme (reversibel κ T x v)'}
for j, kk in enumerate(keys):
    out(f"   {names[kk]:<46s} {y[12 + j] / T2: .6e}")

# ===========================================================================
hr("TEIL 3  Jacobi-Struktur")
# ===========================================================================
p3 = dict(P, g0=0.5)


def jac(fun, z, h=1e-7):
    n = len(z)
    J = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n)
        e[i] = h
        J[:, i] = (fun(z + e) - fun(z - e)) / (2 * h)
    return J


# Vorhersage der Arbeitstemperatur (Herleitung in Teil 5): ΔT = ΔT0 / (1 - q_S2)
_mue0 = P['a'] - P['R0'] - P['Rrad']
_mum = P['am'] - P['c'] - P['cac']
_q1 = -(P['a'] - P['Rrad']) * 2 / (3 * P['b']) * P['alphaR'] * P['R0'] / P['Gth']
_dT0 = ((P['a'] - P['Rrad']) * 2 * _mue0 / (3 * P['b']) + (P['am'] - P['cac']) * 2 * _mum / (3 * P['bm'])) / P['Gth']
Tstar = 1 + _dT0 / (1 - _q1)
out(f"Arbeitspunkt T* = {Tstar:.5f} (aus Teil 5)")
zc = np.array([0.0, 0.0, 0.0, 0.0, Tstar, 0.0, 0.0, 0.5])
Jt = jac(lambda z: F_tan(z, P), zc)
kst = k_of(Tstar, 0, 0.5, P)
Sc = np.diag([1 / np.sqrt(P['C']), np.sqrt(P['L']), np.sqrt(kst), np.sqrt(P['m'])])
Jn = Sc @ Jt[:4, :4] @ np.linalg.inv(Sc)
out("J_tan im Ruhepunkt (EM, Mech), energienormiert (q/√C, √L I, √k x, √m v):")
for row in Jn:
    out("   [" + "  ".join(f"{v: .4f}" for v in row) + "]")
out(f"   ||J + Jᵀ|| = {np.linalg.norm(Jn + Jn.T):.2e}   (exakt antisymmetrisch)")
ev = np.linalg.eigvals(Jn)
out("   Eigenwerte:", "  ".join(f"{e.real:+.1e}{e.imag:+.5f}i" for e in sorted(ev, key=lambda c: c.imag)))
Om = np.sort(np.abs(ev.imag))[::2]
out(f"   Normalmoden Ω1 = {Om[0]:.5f}, Ω2 = {Om[1]:.5f}   (ω_e = 1, ω_m = √k = {np.sqrt(kst):.5f})")


# Torus-Stichprobe: Spur von J_tan in (q,I,x,v,σ) und Mittel der Divergenz des Gesamtfelds
def J_sigma(fun, z, h=1e-7):
    def G(yy):
        zz = yy.copy()
        zz[IT] = np.exp(yy[IT] / P['Cth'])
        f = fun(zz)
        f[IT] = P['Cth'] * f[IT] / zz[IT]
        return f
    yy = z.copy()
    yy[IT] = P['Cth'] * np.log(z[IT])
    return jac(G, yy, h)[:5, :5]


# Kurzer Lauf auf den Torus (g = 0.5)
z = np.array([0.5, 0, 0.3, 0, 1.0, 0, 0, 0.5])
for n in range(int(600 / 0.02)):
    z = rk4(lambda s: F(s, P), z, 0.02)
Jtor_t = J_sigma(lambda s: F_tan(s, P), z)
Jtor_tr = J_sigma(lambda s: F_irr(s, P) + F_port(s, P), z)
p_eps = dict(P, epsD=1e-5, nu0=0.05, epsg=1e-4)
Jfull = jac(lambda s: F(s, p_eps), z)
out(f"Punkt auf dem Torus (g=0.5): tr J_tan (in σ-Koord.) = {np.trace(Jtor_t):.2e};  tr J_trans = {np.trace(Jtor_tr):+.4f}")
out("Blocknormen der vollen Jacobimatrix (ε-Terme aktiv: ε_D=1e-5, ν0=0.05, ġ=1e-4):")
blocks = {'EM(q,I)': [0, 1], 'Mech(x,v)': [2, 3], 'T': [4], 'langsam(D,A,g)': [5, 6, 7]}
hdr = "   Zeile\\Spalte      " + "".join(f"{b:>16s}" for b in blocks)
out(hdr)
for bn, bi in blocks.items():
    row = []
    for cn, ci in blocks.items():
        row.append(np.linalg.norm(Jfull[np.ix_(bi, ci)]))
    out(f"   {bn:<18s}" + "".join(f"{v:16.3e}" for v in row))

# ===========================================================================
hr("TEIL 4  Onsager (dissipativ, symmetrisch) und Casimir (reaktiv, antisymmetrisch)")
# ===========================================================================
Tn = Tstar
Rn = R_of(Tn, 0, P)
Kth = P['Gth']
S_b = P['Sb']
Lmat = np.array([[Tn / Rn, S_b * Tn ** 2 / Rn],
                 [S_b * Tn ** 2 / Rn, S_b ** 2 * Tn ** 3 / Rn + Kth * Tn ** 2]])
out("Thermoelektrisches Element (Kräfte X_e = ΔV/T, X_q = ΔT/T²; Flüsse I, J_q):")
out(f"   L = [[{Lmat[0, 0]:.5f}, {Lmat[0, 1]:.5f}], [{Lmat[1, 0]:.5f}, {Lmat[1, 1]:.5f}]]")
out(f"   L12 - L21 = {Lmat[0, 1] - Lmat[1, 0]:.1e}  (Kelvin: Π = S T  <=>  Onsager L12 = L21)")
out(f"   Eigenwerte L: {np.linalg.eigvalsh(Lmat)[0]:.5f}, {np.linalg.eigvalsh(Lmat)[1]:.5f}  (> 0  => Π_s >= 0)")
Xs = rng.normal(0, 1, (2, 5000))
Pis = np.einsum('in,ij,jn->n', Xs, Lmat, Xs)
out(f"   min X·L·X über 5000 Zufallskräfte = {Pis.min():.3e}")
# Beitrag Seebeck/Peltier zur Entropieproduktion im dynamischen Modell
zz = Z.copy()
pp = p1
dS_node = -pp['Sb'] * zz[II]                              # Knoten: -S T I / T
dS_env = pp['Sb'] * pp['Ta'] * zz[II] / pp['Ta']          # Umgebung: +S T_a I / T_a
out(f"   Seebeck/Peltier-Beitrag zur Gesamt-Entropieproduktion: max |ΔΠ| = {np.max(np.abs(dS_node + dS_env)):.1e}")
# Casimir: Gyrator- und thermoelastischer Block
Jg = Jn
out(f"Gyrator (Lorentz/Induktion): J[Ĩ,ṽ] = {Jg[1, 3]:+.5f}, J[ṽ,Ĩ] = {Jg[3, 1]:+.5f}  -> antisymmetrisch (B-ungerade, Casimir)")
xs, Ts_ = 0.7, Tstar
Jte = jac(lambda s: F_tan(s, P), np.array([0, 0, xs, 0.0, Ts_, 0, 0, 0.5]))
a_vT = np.sqrt(P['m']) * Jte[IV, IT] * np.sqrt(Ts_ / P['Cth'])
a_Tv = np.sqrt(P['Cth'] / Ts_) * Jte[IT, IV] / np.sqrt(P['m'])
out(f"Thermoelastik (Metrik C_th/T):   J[ṽ,T̃] = {a_vT:+.5f}, J[T̃,ṽ] = {a_Tv:+.5f}  -> antisymmetrisch (reversibel)")

# ===========================================================================
hr("TEIL 5  Schleifenverstärkungen und die Reihe Σ qⁿ")
# ===========================================================================
# S2: EM -> Wärme -> Material(R) -> EM
out("S2  EM -> Wärme -> Material σ(T) -> EM")
for q_ in [0.5, 0.9]:
    # Stromgetrieben: C Ṫ = R0 (1 + α ΔT) I0² - G ΔT,  q = α R0 I0² / G
    Gth, R0, al = 1.0, 1.0, 1.0
    I0sq = q_ * Gth / (al * R0)
    dT0 = R0 * I0sq / Gth
    dT = 0.0
    for n in range(200000):
        dT += 0.001 * (R0 * (1 + al * dT) * I0sq - Gth * dT)
    partial = [dT0 * sum(q_ ** j for j in range(N + 1)) for N in (0, 1, 2, 5, 20)]
    out(f"   stromgetrieben q = {q_}: ΔT_sim = {dT:.5f},  ΔT0/(1-q) = {dT0 / (1 - q_):.5f},"
        f"  Teilsummen ΔT0·Σ_(n<=N) qⁿ (N=0,1,2,5,20): " + ", ".join(f"{s:.4f}" for s in partial))
q_ = 1.2
I0sq = q_ * 1.0
dT, t_run = 0.0, 0.0
while dT < 50 and t_run < 1000:
    dT += 0.001 * ((1 + dT) * I0sq - dT)
    t_run += 0.001
out(f"   stromgetrieben q = {q_}: ΔT erreicht 50·ΔT0 nach t = {t_run:.2f} C/G  (Reihe divergiert: thermisches Durchgehen)")
# selbsterregter Kreis (van der Pol): P_heat,EM = (a - R_rad)·2μ_e/(3b)
dPdT = -(P['a'] - P['Rrad']) * 2 / (3 * P['b']) * P['alphaR'] * P['R0']
q1 = dPdT / P['Gth']
mu_e0 = P['a'] - P['R0'] - P['Rrad']
mu_m = P['am'] - P['c'] - P['cac']
Pe0 = (P['a'] - P['Rrad']) * 2 * mu_e0 / (3 * P['b'])
Pm0 = (P['am'] - P['cac']) * 2 * mu_m / (3 * P['bm'])
dT0 = (Pe0 + Pm0) / P['Gth']
dTpred = dT0 / (1 - q1)
out(f"   selbsterregt (Minimalsystem): q_S2 = ∂P/∂T / G = {q1:+.4f}  (negativ: stabilisierend)")
out(f"   ΔT0 (ohne Rücklauf) = {dT0:.5f},  ΔT = ΔT0/(1-q) = {dTpred:.5f}  -> T* = {1 + dTpred:.5f}")
# S3: Selbsterregung (Barkhausen), Steigung -> λ⊥
A2e = 4 * (mu_e0 - P['R0'] * P['alphaR'] * dTpred) / (3 * P['b'])
mu_e = mu_e0 - P['R0'] * P['alphaR'] * dTpred
A2m = 4 * mu_m / (3 * P['bm'])
out("S3  Selbsterregung (aktives Element bzw. Strömung):  q_self(A) = a / (R_tot + 3bA²/4)")
out(f"   EM:   A*² = {A2e:.5f},  q_self(A*) = {P['a'] / (R_of(1 + dTpred, 0, P) + P['Rrad'] + 0.75 * P['b'] * A2e):.6f},"
    f"  λ⊥,EM = (aA*/2L)·∂_A ln q = -μ_e = {-mu_e:.5f}")
out(f"   Mech: A*² = {A2m:.5f},  q_self(A*) = {P['am'] / (P['c'] + P['cac'] + 0.75 * P['bm'] * A2m):.6f},"
    f"  λ⊥,Mech = -μ_m = {-mu_m:.5f}")
lamT_pred = -(P['Gth'] - dPdT) / P['Cth']
out(f"   Wärme: λ_T = -(G - ∂P/∂T)/C_th = {lamT_pred:.5f}")
# S1: Lorentz/Induktion-Schleife, q(ω) = Γ²ω² / (L m (ω²-ω_e²)(ω²-ω_m²))


def q_LI(w, we2=1.0, wm2=kst):
    return P['Gam'] ** 2 * w ** 2 / (P['L'] * P['m'] * (w ** 2 - we2) * (w ** 2 - wm2))


out("S1  EM -> Lorentz -> Mech -> Bewegungs-EMK -> EM:  q_LI(ω) = Γ²ω² / (L m (ω²-ω_e²)(ω²-ω_m²))")
out(f"   q_LI(Ω1) = {q_LI(Om[0]):.9f},  q_LI(Ω2) = {q_LI(Om[1]):.9f}   (Normalmoden = Barkhausen q = 1)")
for w in [0.5, 0.95, 1.2]:
    out(f"   ω = {w}: q_LI = {q_LI(w):+.5f},  Σqⁿ = 1/(1-q) = {1 / (1 - q_LI(w)):.5f}")
# Netzwerkform: (I - A)^-1 = Σ Aⁿ
Anet = np.array([[0, 0, 1 / P['Gth']],                  # ΔT <- ΔP
                 [P['R0'] * P['alphaR'], 0, 0],          # ΔR <- ΔT
                 [0, -(P['a'] - P['Rrad']) * 2 / (3 * P['b']), 0]])   # ΔP <- ΔR (über μ_e, A²)
rho_sp = max(abs(np.linalg.eigvals(Anet)))
Inv = np.linalg.inv(np.eye(3) - Anet)
Neu = sum(np.linalg.matrix_power(Anet, n) for n in range(60))
out(f"Netzwerk (ΔT, ΔR, ΔP): Spektralradius ρ(A) = {rho_sp:.5f} = |q_S2|^(1/3) = {abs(q1) ** (1 / 3):.5f}")
out(f"   ||(I-A)^-1 - Σ_(n<60) Aⁿ|| = {np.linalg.norm(Inv - Neu):.2e};  Gesamtverstärkung ΔP0 -> ΔT: {Inv[0, 2]:.5f} = (1/G)/(1-q) = {1 / P['Gth'] / (1 - q1):.5f}")
# S5: elektrostatischer Pull-in als statische Schleife
kk_, d_, eA = 1.0, 1.0, 1.0
Vpi = np.sqrt(8 * kk_ * d_ ** 3 / (27 * eA))
xpi = d_ / 3
q_pi = eA * Vpi ** 2 / (kk_ * (d_ - xpi) ** 3)
out(f"S5  Maxwell-Druck -> Spalt -> Feld (Pull-in): q_ES = ε0A V²/(k (d-x)³);  am Pull-in-Punkt q_ES = {q_pi:.6f}")
# S7: Ermüdung als integrierende positive Schleife (Amplitude der Geschwindigkeit fest)
mDexp = P['mD']
out(f"S7  Ermüdung: Ḋ = ε_D (Y/Y_ref)^m, Y ∝ x² ∝ 1/k ∝ 1/(1-D)  ->  Ḋ ∝ (1-D)^(-m), m = {mDexp:.0f}")
out(f"   Lebensdauer t_f = (1-D0)^(m+1) / ((m+1) Ḋ0) = {1 / ((mDexp + 1)):.4f} / Ḋ0  (endliche Zeit statt Konvergenz)")
# S8: Arrhenius-Empfindlichkeit
out(f"S8  Alterung: ∂ln Ȧ/∂T = E_a/T² = {P['Ea'] / (1 + dTpred) ** 2:.4f} (Modell);"
    f" real E_a = 0.7 eV, 300->310 K: Faktor {np.exp(0.7 * 1.602176634e-19 / 1.380649e-23 * (1 / 300 - 1 / 310)):.3f}")

# ===========================================================================
hr("TEIL 6  Torus-Attraktor: Lyapunov-Signatur, Rotationszahl, Poincaré-Schnitt")
# ===========================================================================


def lyap_run(g, Ttot=3000.0, Ttr=500.0, dt=0.02, delta=1e-7, renorm=5):
    p = dict(P)
    z0 = np.array([0.5, 0, 0.3, 0, 1.0, 0, 0, g])
    Zb = np.tile(z0[:, None], (1, 6))
    Q = np.eye(8)[:, :5]
    for j in range(5):
        Zb[:, j + 1] = z0 + delta * Q[:, j]
    fun = lambda s: F(s, p)
    nst, ntr = int(Ttot / dt), int(Ttr / dt)
    sums = np.zeros(5)
    ae, am_ = phases(Zb[:, 0], p)
    th_e = th_m = 0.0
    divsum, lorsum, cnt = 0.0, 0.0, 0
    cpsi = spsi = 0.0
    kprev = None
    occ = np.zeros(72, bool)
    ncross = 0
    for n in range(nst):
        Zb = rk4(fun, Zb, dt)
        z = Zb[:, 0]
        be, bm = phases(z, p)
        if n >= ntr:
            th_e += wrap(be - ae)
            th_m += wrap(bm - am_)
            psi = bm - be
            cpsi += np.cos(psi)
            spsi += np.sin(psi)
            q, I, x, v, T, D, A, gg = z
            R = R_of(T, A, p)
            divsum += ((p['a'] - p['Rrad'] - R - 3 * p['b'] * I * I) / p['L']
                       + (p['am'] - p['cac'] - p['c'] - 3 * p['bm'] * v * v) / p['m']
                       + (-p['kappaT'] * x * v + p['R0'] * p['alphaR'] * I * I - p['Gth'] - p['Sb'] * I) / p['Cth'])
            lorsum += p['Gam'] * I * v
            cnt += 1
            kc = np.floor(th_e / (2 * np.pi))
            if kprev is not None and kc != kprev:
                occ[int((bm % (2 * np.pi)) / (2 * np.pi) * 72) % 72] = True
                ncross += 1
            kprev = kc
        ae, am_ = be, bm
        if (n + 1) % renorm == 0:
            Dm = (Zb[:5, 1:] - Zb[:5, [0]]) / delta
            Qm, Rm = np.linalg.qr(Dm)
            s = np.sign(np.diag(Rm))
            Qm = Qm * s
            if n >= ntr:
                sums += np.log(np.abs(np.diag(Rm)))
            Zb[:, 1:] = Zb[:, [0]]
            Zb[:5, 1:] += delta * Qm
    lam = sums / (Ttot - Ttr)
    return dict(lam=lam, rho=th_m / th_e, Oe=th_e / (Ttot - Ttr), Om=th_m / (Ttot - Ttr),
                div=divsum / cnt, lor=lorsum / cnt, Rpsi=np.hypot(cpsi, spsi) / cnt,
                fill=occ.mean(), ncross=ncross, zend=Zb[:, 0])


for g_, label in [(0.5, "quasiperiodisch"), (0.80, "eingerastet")]:
    r = lyap_run(g_)
    lam = np.sort(r['lam'])[::-1]
    out(f"g = {g_} ({label}):")
    out("   Lyapunov-Spektrum (q,I,x,v,T):  " + "  ".join(f"{l:+.5f}" for l in lam))
    out(f"   Σλ = {lam.sum():+.5f},  <div F> = {r['div']:+.5f}")
    out(f"   Ω_EM = {r['Oe']:.6f}, Ω_mech = {r['Om']:.6f}, Rotationszahl ρ = Ω_mech/Ω_EM = {r['rho']:.6f}")
    out(f"   Ordnungsparameter der Phasendifferenz |<e^iψ>| = {r['Rpsi']:.4f}")
    out(f"   Poincaré-Schnitt θ_EM = 0 mod 2π: {r['ncross']} Durchgänge, belegte θ_mech-Bins (72): {r['fill'] * 100:.1f} %")
    out(f"   mittlerer Lorentz-Energiefluss EM->Mech  Γ<I v> = {r['lor']:+.5f}")
    out(f"   T* = {r['zend'][IT]:.4f}")
    if g_ == 0.5:
        out(f"   Vorhersage (Teil 5): λ⊥,EM = {-mu_e:.4f}, λ⊥,Mech = {-mu_m:.4f}, λ_T = {lamT_pred:.4f}")

# ===========================================================================
hr("TEIL 7  Kollaps T² -> T¹ : Rast-Scan in g, langsame Rampe, Alterung")
# ===========================================================================
gs = np.array([0.20, 0.30, 0.40, 0.50, 0.60, 0.65, 0.70, 0.72, 0.74, 0.75, 0.76, 0.78, 0.80, 0.82, 0.84,
               0.86, 0.87, 0.88, 0.90, 0.95, 1.00])
As = np.array([0.30, 0.50, 0.60, 0.65, 0.70, 0.80])
Ng, Na = len(gs), len(As)
N = Ng + Na
Zs = np.zeros((8, N))
Zs[IQ], Zs[IX], Zs[IT] = 0.5, 0.3, 1.0
Zs[IG, :Ng] = gs
Zs[IG, Ng:] = 0.5
Zs[IA, Ng:] = As
dt7, T7, Ttr7 = 0.02, 2500.0, 1000.0
the = np.zeros(N)
thm = np.zeros(N)
cps = np.zeros(N)
sps = np.zeros(N)
Imax = np.zeros(N)
vmax = np.zeros(N)
cnt = 0
ae, am_ = phases(Zs, P)
fun = lambda s: F(s, P)
for n in range(int(T7 / dt7)):
    Zs = rk4(fun, Zs, dt7)
    be, bm = phases(Zs, P)
    if n >= int(Ttr7 / dt7):
        the += wrap(be - ae)
        thm += wrap(bm - am_)
        cps += np.cos(bm - be)
        sps += np.sin(bm - be)
        Imax = np.maximum(Imax, np.abs(Zs[II]))
        vmax = np.maximum(vmax, np.abs(Zs[IV]))
        cnt += 1
    ae, am_ = be, bm
rho = thm / the
Rpsi = np.hypot(cps, sps) / cnt
out("Rast-Scan (A = 0):  g,  ω_m,0 = √(k(T*,0,g)),  ρ = Ω_m/Ω_e,  |<e^iψ>|,  max|I|,  max|v|")
for j in range(Ng):
    wm0 = np.sqrt(k_of(Zs[IT, j], 0, gs[j], P))
    out(f"   g = {gs[j]:.2f}   ω_m,0 = {wm0:.4f}   ρ = {rho[j]:.5f}   |<e^iψ>| = {Rpsi[j]:.4f}"
        f"   max|I| = {Imax[j]:.3f}   max|v| = {vmax[j]:.3f}")
locked = np.abs(rho[:Ng] - 1) < 2e-4
if locked.any():
    out(f"   Rastfenster (|ρ-1| < 2e-4): g ∈ [{gs[:Ng][locked].min():.2f}, {gs[:Ng][locked].max():.2f}]")
out("Alterungs-Scan (g = 0.5, A fest): Zyklustod des EM-Umlaufs bei μ_e(A) = 0")
for j in range(Na):
    jj = Ng + j
    Tloc = Zs[IT, jj]
    mu = P['a'] - P['Rrad'] - R_of(Tloc, As[j], P)
    out(f"   A = {As[j]:.2f}:  μ_e = a - R_rad - R(T*,A) = {mu:+.4f}   max|I| = {Imax[jj]:.4f}   max|v| = {vmax[jj]:.3f}"
        f"   T* = {Tloc:.4f}")
# kritisches A (selbstkonsistent, nur mech. Heizung wenn EM tot)
Tm = 1 + Pm0 / P['Gth']
Ac = ((P['a'] - P['Rrad']) / (P['R0'] * (1 + P['alphaR'] * (Tm - 1))) - 1) / P['betaA']
out(f"   kritisches A_c (μ_e = 0 bei T = T_a + P_mech/G = {Tm:.4f}): A_c = {Ac:.4f}")

# Langsame Rampe g auf / ab (Hysterese)
eps_r = 5e-5
gA, gB = 0.60, 0.95
Tr = (gB - gA) / eps_r
Zr = np.zeros((8, 2))
Zr[IQ], Zr[IX], Zr[IT] = 0.5, 0.3, 1.0
Zr[IG] = [gA, gB]
pr = dict(P)
dtr = 0.02
# Einschwingen bei festem g
for n in range(int(500 / dtr)):
    Zr = rk4(lambda s: F(s, pr), Zr, dtr)
pr['epsg'] = np.array([eps_r, -eps_r])
win = 250.0
nwin = int(win / dtr)
the_w = np.zeros(2)
thm_w = np.zeros(2)
ae, am_ = phases(Zr, pr)
rows = []
for n in range(int(Tr / dtr)):
    Zr = rk4(lambda s: F(s, pr), Zr, dtr)
    be, bm = phases(Zr, pr)
    the_w += wrap(be - ae)
    thm_w += wrap(bm - am_)
    ae, am_ = be, bm
    if (n + 1) % nwin == 0:
        rows.append((Zr[IG].copy(), thm_w / the_w))
        the_w[:] = 0
        thm_w[:] = 0
out(f"Langsame Rampe ġ = ±{eps_r} (ε = ġ/ω ≈ {eps_r:.0e}), Fenster Δt = {win:.0f}: Rotationszahl im Fenster")
out("   g(auf)   ρ_auf      |  g(ab)    ρ_ab")
for (gg, rr) in rows:
    out(f"   {gg[0]:.4f}   {rr[0]:.5f}   |  {gg[1]:.4f}   {rr[1]:.5f}")
# Kriterium: im Rastzustand driftet nur die Rastphase ψ*(g) adiabatisch mit
# (|ρ-1| ~ |dψ*/dt|/ω ≲ 3e-3); ein einziger Phasensprung pro Fenster gäbe |ρ-1| >= 1/40.
tol_r = 5e-3
out(f"   Rastkriterium in der Rampe: |ρ-1| < {tol_r} (ein Phasensprung pro Fenster ergäbe >= {2*np.pi/(win):.3f})")
up_lock = [gg[0] for gg, rr in rows if abs(rr[0] - 1) < tol_r]
dn_lock = [gg[1] for gg, rr in rows if abs(rr[1] - 1) < tol_r]
if up_lock:
    out(f"   Rast bei Rampe auf:  g ∈ [{min(up_lock):.4f}, {max(up_lock):.4f}] (Fenstermitten)")
if dn_lock:
    out(f"   Rast bei Rampe ab:   g ∈ [{min(dn_lock):.4f}, {max(dn_lock):.4f}] (Fenstermitten)")

# Bistabilitätstest: statische Läufe, gestartet im eingerasteten Zustand (g = 0.80)
zl = np.array([0.5, 0, 0.3, 0, 1.0, 0, 0, 0.80])
for n in range(int(1500 / 0.02)):
    zl = rk4(lambda s_: F(s_, P), zl, 0.02)
gb = np.array([0.72, 0.74, 0.87, 0.88, 0.89])
Zb2 = np.tile(zl[:, None], (1, len(gb)))
Zb2[IG] = gb
the_b = np.zeros(len(gb))
thm_b = np.zeros(len(gb))
ae, am_ = phases(Zb2, P)
for n in range(int(3000 / 0.02)):
    Zb2 = rk4(lambda s_: F(s_, P), Zb2, 0.02)
    be, bm = phases(Zb2, P)
    if n >= int(1500 / 0.02):
        the_b += wrap(be - ae)
        thm_b += wrap(bm - am_)
    ae, am_ = be, bm
out("Bistabilitätstest (fester Parameter, Start im Rastzustand von g = 0.80, t = 3000):")
for j, gg in enumerate(gb):
    out(f"   g = {gg:.2f}:  ρ = {thm_b[j] / the_b[j]:.5f}")

# ===========================================================================
hr("TEIL 8  Gravitation als langsamer Parameter – FTM-Größen")
# ===========================================================================
c0 = 299792458.0
bodies = [("Erde", 3.986004418e14, 6.371e6), ("Sonne", 1.32712440018e20, 6.957e8),
          ("Sirius B", 1.018 * 1.32712440018e20, 5.8e6), ("Neutronenstern 1.4 M☉", 1.4 * 1.32712440018e20, 1.2e4)]
out("Körper, g = GM/R², Δs = GM/(3c²), dτ/dt = √(1-6Δs/R), 1-dτ/dt")
for nm, GM, Rb in bodies:
    ds = GM / (3 * c0 ** 2)
    gsurf = GM / Rb ** 2
    dtau = np.sqrt(1 - 6 * ds / Rb)
    out(f"   {nm:<22s} g = {gsurf:.4e} m/s²   Δs = {ds:.5e} m   dτ/dt = {dtau:.12f}   1-dτ/dt = {1 - dtau:.4e}")
ds_E = 3.986004418e14 / (3 * c0 ** 2)
out(f"   Kontrolle Erde: Δs⊕ = {ds_E * 1e3:.4f} mm,  Schwarzschild-Radius r₀ = 6Δs⊕ = {6 * ds_E * 1e3:.3f} mm")
# Invarianz der Rotationszahl unter Zeitdilatation
we_loc, wm_loc = 1.0, 0.8307
for nm, GM, Rb in bodies:
    ds = GM / (3 * c0 ** 2)
    f_ = np.sqrt(1 - 6 * ds / Rb)
    out(f"   {nm:<22s} ω_EM,obs = {we_loc * f_:.12f}, ω_m,obs = {wm_loc * f_:.12f},"
        f" ρ_obs - ρ_loc = {(wm_loc * f_) / (we_loc * f_) - wm_loc / we_loc:.1e}")
# Pound-Rebka-Höhe und Vorspannungskanal
h = 22.5
gE = 3.986004418e14 / 6.371e6 ** 2
out(f"Zeitdilatations-Gradient: Δω/ω = g h / c² bei h = {h} m: {gE * h / c0 ** 2:.4e}")
h2 = 1000.0
dg_g = (6.371e6 / (6.371e6 + h2)) ** 2 - 1
out(f"System um h = {h2:.0f} m angehoben:  δg/g = {dg_g:.4e}")
out(f"   Pendel-Zyklus (ω_m = √(g/l)): δω_m/ω_m = {np.sqrt(1 + dg_g) - 1:.4e};  EM-Zyklus: δω_e/ω_e = 0")
out(f"   gemeinsamer Zeitdilatationsfaktor: δ(dτ/dt) = {np.sqrt(1 - 6 * ds_E / (6.371e6 + h2)) - np.sqrt(1 - 6 * ds_E / 6.371e6):.4e}")
out(f"   Verhältnis differentieller Kanal / gemeinsamer Kanal = "
    f"{abs(np.sqrt(1 + dg_g) - 1) / abs(np.sqrt(1 - 6 * ds_E / (6.371e6 + h2)) - np.sqrt(1 - 6 * ds_E / 6.371e6)):.3e}")
# Modell: Empfindlichkeit der Rotationszahl auf g
Tq = 1 + dTpred
wm05 = np.sqrt(k_of(Tq, 0, 0.5, P))
out(f"Modell: ∂ln ω_m/∂ln g = (m g/l)/(2k) bei g = 0.5: {0.5 / (2 * k_of(Tq, 0, 0.5, P)):.4f}")
gc = P['kel'] * P['ell'] / P['m'] - P['kappaT'] * (Tq - 1) * P['ell'] / P['m']
out(f"Druckvorspannung (k = k_el(1-D) - κΔT - m g/l): Knick-/Zyklusverlust bei g_c = {gc:.4f}")
# g-Wert für ω_m = ω_e (Resonanzmitte) im Modell
g_res = (1.0 - (P['kel'] - P['kappaT'] * (Tq - 1))) * P['ell'] / P['m']
out(f"Resonanzmitte ω_m(g) = ω_e bei g = {g_res:.4f}")

# ===========================================================================
hr("TEIL 9  Fluktuations-Dissipations-Beziehung")
# ===========================================================================
rngf = np.random.default_rng(7)
Mreal, kB_T, gam, kf, mf = 4000, 0.01, 0.1, 1.0, 1.0
dtf, Tf = 0.01, 300.0
xf = np.zeros(Mreal)
vf = np.zeros(Mreal)
acc_v, acc_x, nacc = 0.0, 0.0, 0
c1 = np.exp(-gam * dtf / mf)
c2 = np.sqrt(kB_T / mf * (1 - c1 ** 2))
for n in range(int(Tf / dtf)):
    # BAOAB: exakter Ornstein-Uhlenbeck-Schritt mit <ξξ> = 2 γ k_B T δ
    vf += 0.5 * dtf * (-kf * xf) / mf
    xf += 0.5 * dtf * vf
    vf = c1 * vf + c2 * rngf.normal(size=Mreal)
    xf += 0.5 * dtf * vf
    vf += 0.5 * dtf * (-kf * xf) / mf
    if n * dtf > 100:
        acc_v += np.mean(0.5 * mf * vf ** 2)
        acc_x += np.mean(0.5 * kf * xf ** 2)
        nacc += 1
out(f"Langevin-Oszillator, Rauschen <ξξ> = 2 c k_B T δ mit c = {gam}, k_B T = {kB_T}, {Mreal} Realisierungen:")
out(f"   <½ m v²>/(½ k_B T) = {acc_v / nacc / (0.5 * kB_T):.4f},   <½ k x²>/(½ k_B T) = {acc_x / nacc / (0.5 * kB_T):.4f}")
kB = 1.380649e-23
out(f"Johnson-Nyquist S_V = 4 k_B T R bei R = 1 kΩ, T = 300 K: {4 * kB * 300 * 1e3:.4e} V²/Hz"
    f" = {np.sqrt(4 * kB * 300 * 1e3) * 1e9:.3f} nV/√Hz")

out("")
out(f"Rechenzeit: {time.time() - t_start:.1f} s")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "wechselwirkungen.txt"), "w",
          encoding="utf-8") as fh:
    fh.write(_BUF.getvalue())
