# -*- coding: utf-8 -*-
"""
Durchrechnung zur Zustandsraum-Schicht des Weltmodells (05-zustandsraum.md).

Setzt die Formeln des Autors in Zahlen um:
  A  Asymptotisches Profil a_extra ~ gamma*c1/r fuer v0 = 220 km/s
  B  Wachstum M_DM ~ r und effektive Masse des Bindungsoperators
  C  Skalenfaktor f(l/lambda) aus Phi_eff = Phi_bary (1 + f)
  D  Umlaufzeiten T ~ r^(3/2) (Kepler) und T ~ r (flach)
  E  q-Entwicklung: Gewicht der hoeheren Ordnungen (Ramanujan-Spitzenform Delta)
  F  Modulare Werkzeuge: Transformation, E2 = P, Q, R, Ramanujan-DGLn, Operator D_k
  G  1D-Toy-Model mit Memory-Kernel (Schliessungsstaerke gamma, Gedaechtniszeit tau_K)
  H  FTM-Anschluss: Hull-Beschleunigung g und Faltungstiefe Delta_s = GM/(3c^2)

Nur Standardbibliothek. Ausgabe: zustandsraum.txt im selben Ordner.
"""
import cmath
import math
import os

OUT = []


def p(s=""):
    OUT.append(s)


def hr(title):
    p()
    p("=" * 78)
    p(title)
    p("=" * 78)


# ---------------------------------------------------------------------------
# Konstanten (SI)
# ---------------------------------------------------------------------------
G = 6.674e-11            # m^3 kg^-1 s^-2  (einheitliche Konstanten 01-08)
c = 2.99792458e8         # m/s
M_SUN = 1.989e30         # kg
KPC = 3.085677581e19     # m
AU = 1.495978707e11      # m
MYR = 3.15576e13         # s
G_EARTH = 9.81           # m/s^2 (Bezugswert g_Erde des FTM)

# Eingaben
V0 = 220e3               # m/s, Milchstrasse (Wert aus dem Toy Paper)
M_B = 6.0e10 * M_SUN     # Rechenannahme: baryonische Masse als Punktmasse (frei waehlbar)

# ---------------------------------------------------------------------------
hr("A  Asymptotisches Profil  a_extra ~ gamma*c1/r   (v0 = 220 km/s)")
# ---------------------------------------------------------------------------
# a_total = v0^2/r  (flache Rotation);  Fernfeld: a_extra ~ gamma*c1 * q, q = r_q/r
# => gamma*c1*r_q = v0^2.  Mit q = 1/r (r_q = 1 Laengeneinheit): gamma*c1 = v0^2.
gc1 = V0 ** 2
p(f"v0                         = {V0/1e3:.1f} km/s")
p(f"gamma*c1 (q = r_q/r)       = v0^2 = {gc1:.4e} m^2/s^2 = {gc1/1e6:.0f} (km/s)^2")
p(f"  als Beschleunigung*Laenge = {gc1/KPC:.4e} m/s^2 * kpc")
p(f"  gamma*c1 in Einheiten r_q: gamma*c1 = v0^2 / r_q")
for rq_name, rq in [("1 m", 1.0), ("1 AU", AU), ("1 kpc", KPC)]:
    p(f"     r_q = {rq_name:6s} -> gamma*c1 = {gc1/rq:.4e} m/s^2")

r_x = G * M_B / V0 ** 2
p()
p(f"Rechenannahme M_bary       = {M_B/M_SUN:.2e} M_sun (Punktmasse)")
p(f"Uebergangsradius r_x       = GM/v0^2 = {r_x/KPC:.3f} kpc")
p("  (dort gilt GM/r^2 = gamma*c1/r;  a_extra/a_Newton = r/r_x)")
p()
p(f"{'r [kpc]':>8} {'a_Newton [m/s^2]':>17} {'a_extra [m/s^2]':>16} "
  f"{'a_extra/a_N':>12} {'v_N [km/s]':>11} {'v_tot [km/s]':>13}")
for rk in [0.5, 1, 2, 5, r_x / KPC, 8.2, 10, 20, 30, 50, 100, 200]:
    r = rk * KPC
    aN = G * M_B / r ** 2
    ax = gc1 / r
    vN = math.sqrt(G * M_B / r)
    vt = math.sqrt(r * (aN + ax))      # v^2 = GM/r + gamma*c1
    p(f"{rk:8.2f} {aN:17.4e} {ax:16.4e} {ax/aN:12.4f} {vN/1e3:11.1f} {vt/1e3:13.1f}")
p()
p("Lesart: v_tot^2 = GM/r + gamma*c1  ->  v_tot -> sqrt(gamma*c1) = v0 fuer r >> r_x.")

# ---------------------------------------------------------------------------
hr("B  Masse: M_DM ~ r (Quelle)  vs.  M_eff(r) des Bindungsoperators (Struktur)")
# ---------------------------------------------------------------------------
slope = V0 ** 2 / G
p(f"dM_DM/dr (asymptotisch)    = v0^2/G = {slope:.4e} kg/m = {slope*KPC/M_SUN:.4e} M_sun/kpc")
p("Standardansatz exakt:  v0^2/r = GM/r^2 + a_DM  ->  M_DM(r) = v0^2 r/G - M_bary")
p("Operatoransatz, 1. Ordnung:  a_total = GM/r^2 + gamma*c1/r  ->")
p("   M_eff(r) := r^2 a_total / G = M_bary + (gamma*c1/G) r   (keine Quelle, Operator-Term)")
p()
p(f"{'r [kpc]':>8} {'M_DM [M_sun]':>14} {'M_DM/M_bary':>12} {'M_eff [M_sun]':>15} "
  f"{'M_dyn=v0^2r/G':>15}")
for rk in [r_x / KPC, 10, 20, 50, 100, 200]:
    r = rk * KPC
    m_dm = V0 ** 2 * r / G - M_B
    if abs(m_dm) < 1e-6 * M_B:
        m_dm = 0.0
    m_eff = M_B + gc1 * r / G
    m_dyn = V0 ** 2 * r / G
    p(f"{rk:8.2f} {m_dm/M_SUN:14.4e} {m_dm/M_B:12.3f} {m_eff/M_SUN:15.4e} {m_dyn/M_SUN:15.4e}")
p()
p("M_DM(r) waechst linear mit Steigung v0^2/G; M_eff(r) waechst mit derselben Steigung")
p("gamma*c1/G = v0^2/G.  Beide Bilanzen sind im Fernfeld zahlengleich; der Unterschied")
p("liegt allein darin, ob der lineare Anteil als Quelle (rho_DM) oder als Operator-Term")
p("(gamma*M_k(tau)) gebucht wird.")
p(f"Verdopplungsradius (M_DM = M_bary): r = 2 r_x = {2*r_x/KPC:.2f} kpc")
p()
p("Reifikation als Rechenschritt (Poisson-Inversion des Beobachters):")
p("  rho_X := (1/(4 pi G)) div(a_obs - a_bary) ;  a_extra = v0^2/r  ->  rho_X(r) = v0^2/(4 pi G r^2)")
PC = KPC / 1e3
for rk in [5, 8.2, 20, 50, 100]:
    r = rk * KPC
    rho = V0 ** 2 / (4 * math.pi * G * r ** 2)
    p(f"   r = {rk:6.1f} kpc: rho_X = {rho:.4e} kg/m^3 = {rho*PC**3/M_SUN:.4e} M_sun/pc^3")
p("  Die Dichte ist Ergebnis der Inversion (Ausgabe), nicht Eingabe der Dynamik.")

# ---------------------------------------------------------------------------
hr("C  Skalenfaktor f(l/lambda) aus Phi_eff = Phi_bary (1 + f)")
# ---------------------------------------------------------------------------
p("Forderung: -dPhi_eff/dr = -(GM/r^2 + v0^2/r),  Phi_bary = -GM/r")
p("=> Phi_eff = -GM/r + v0^2 ln(r/r_ref)")
p("=> f(u) = u (C - ln u),  u = l/lambda,  lambda = GM/v0^2 = r_x,  C frei (Potentialkonstante)")
p("   (der Anteil C*u erzeugt keine Kraft: Phi_bary*C*u = -GM C/lambda = const)")
p("   Referenzradius im Logarithmus: r_ref (beliebig, nur Potentialkonstante)")
p()
p(f"{'r [kpc]':>8} {'u=r/lambda':>11} {'f(u), C=0':>11} {'|f|':>9} {'Kontrolle a [m/s^2]':>21} {'v0^2/r+GM/r^2':>15}")
for rk in [1, 2, 5, 10, 20, 50, 100]:
    r = rk * KPC
    u = r / r_x
    f = -u * math.log(u)
    # numerische Kontrolle: a = d/dr [Phi_bary (1+f)]
    h = r * 1e-6

    def phi_eff(rr):
        uu = rr / r_x
        return -G * M_B / rr * (1 - uu * math.log(uu))
    a_num = (phi_eff(r + h) - phi_eff(r - h)) / (2 * h)
    a_soll = V0 ** 2 / r + G * M_B / r ** 2
    p(f"{rk:8.1f} {u:11.4f} {f:11.4f} {abs(f):9.4f} {a_num:21.6e} {a_soll:15.6e}")

# ---------------------------------------------------------------------------
hr("D  Umlaufzeit: Kepler T ~ r^(3/2)  vs.  flache Rotation T ~ r")
# ---------------------------------------------------------------------------
p(f"{'r [kpc]':>8} {'T_Kepler [Myr]':>15} {'T_flach [Myr]':>14} {'T_K/T_f':>9}")
for rk in [1, 5, 8.2, 10, 20, 50, 100]:
    r = rk * KPC
    tk = 2 * math.pi * math.sqrt(r ** 3 / (G * M_B))
    tf = 2 * math.pi * r / V0
    p(f"{rk:8.1f} {tk/MYR:15.1f} {tf/MYR:14.1f} {tk/tf:9.3f}")
p("Skalierung: T_K ~ r^1.5 -> Verdopplung von r: Faktor 2^1.5 = %.3f; T_f: Faktor 2." % 2 ** 1.5)
p("Toy-Paper-Kopplung: lambda ~ T/t_relax, q ~ 1/r.")

# ---------------------------------------------------------------------------
hr("E  q-Entwicklung M_k(tau) = sum c_n q^n  mit c_n = tau(n) (Delta, Gewicht 12)")
# ---------------------------------------------------------------------------
NMAX = 40


def delta_coeffs(n):
    """Koeffizienten von q * prod (1-q^m)^24 bis q^n."""
    poly = [0] * (n + 1)
    poly[0] = 1
    for m in range(1, n + 1):
        for _ in range(24):
            for k in range(n, m - 1, -1):
                poly[k] -= poly[k - m]
    return [0] + poly[:n]     # Index = Potenz von q


tau = delta_coeffs(NMAX)
p("tau(n), n=1..12: " + ", ".join(str(tau[n]) for n in range(1, 13)))
p()
p("Gewicht der hoeheren Ordnungen relativ zum ersten Term (c1 q):")
p(f"{'q':>8} {'c2 q^2/(c1 q)':>15} {'c3 q^3/(c1 q)':>15} {'Summe/(c1 q)':>14}")
for q in [0.1, 0.03, 0.01, 1 / 240, 1e-3, 1e-4]:
    s = sum(tau[n] * q ** n for n in range(1, NMAX + 1))
    p(f"{q:8.5f} {tau[2]*q**2/(q):15.5f} {tau[3]*q**3/(q):15.6f} {s/q:14.6f}")
p("Erste Ordnung auf 1 % genau fuer |q| < 1/2400 (da |c2/c1| = 24).")
p("Mit q = r_q/r: Dominanz der ersten Ordnung (1 %) fuer r > 2400 r_q.")

# ---------------------------------------------------------------------------
hr("F  Modulare Werkzeuge: Transformation, P,Q,R, Ramanujan-DGLn, Operator D_k")
# ---------------------------------------------------------------------------


def sigma(k, n):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


NS = 80


def P(q):
    return 1 - 24 * sum(sigma(1, n) * q ** n for n in range(1, NS))


def Q(q):
    return 1 + 240 * sum(sigma(3, n) * q ** n for n in range(1, NS))


def R(q):
    return 1 - 504 * sum(sigma(5, n) * q ** n for n in range(1, NS))


def dq(fun, q, h=1e-7):
    return (fun(q + h) - fun(q - h)) / (2 * h)


q0 = 0.05
thP, thQ, thR = q0 * dq(P, q0), q0 * dq(Q, q0), q0 * dq(R, q0)
Pv, Qv, Rv = P(q0), Q(q0), R(q0)
p(f"Pruefpunkt q = {q0}")
p(f"  q dP/dq = {thP:.10f}   (P^2-Q)/12  = {(Pv**2-Qv)/12:.10f}")
p(f"  q dQ/dq = {thQ:.10f}   (PQ-R)/3    = {(Pv*Qv-Rv)/3:.10f}")
p(f"  q dR/dq = {thR:.10f}   (PR-Q^2)/2  = {(Pv*Rv-Qv**2)/2:.10f}")


def Delta(q):
    return sum(tau[n] * q ** n for n in range(1, NMAX + 1))


Dv = Delta(q0)
p(f"  Delta = (Q^3-R^2)/1728: {Dv:.12e}  vs. {(Qv**3-Rv**2)/1728:.12e}")
thD = q0 * dq(Delta, q0)
p(f"  q dDelta/dq = {thD:.12e}   P*Delta = {Pv*Dv:.12e}")
p(f"  D_12 Delta = q dDelta/dq - P Delta = {thD-Pv*Dv:.3e}   (= 0)")
p(f"  D_4 Q      = q dQ/dq - (1/3) P Q   = {thQ-Pv*Qv/3:.10f}   -R/3   = {-Rv/3:.10f}")
p(f"  D_6 R      = q dR/dq - (1/2) P R   = {thR-Pv*Rv/2:.10f}   -Q^2/2 = {-Qv**2/2:.10f}")

p()
p("Transformation unter z -> -1/z (a,b,c,d) = (0,-1,1,0), Pruefpunkt z = 0.3 + 1.1 i:")
z = complex(0.3, 1.1)


def qz(zz):
    return cmath.exp(2j * math.pi * zz)


zi = -1 / z
p(f"  E4(-1/z)/(z^4 E4(z))       = {Q(qz(zi))/(z**4*Q(qz(z)))}")
p(f"  Delta(-1/z)/(z^12 Delta(z)) = {Delta(qz(zi))/(z**12*Delta(qz(z)))}")
anom = P(qz(zi)) - z ** 2 * P(qz(z))
p(f"  E2(-1/z) - z^2 E2(z)        = {anom}")
p(f"  12 c (cz+d)/(2 pi i), c=1,d=0 = {12*z/(2j*math.pi)}")
p("  -> E2(gz) = (cz+d)^2 E2(z) + 12c(cz+d)/(2 pi i)")
p("     (dieser Term hebt den Ableitungsterm c k (cz+d)^(k+1) f/(2 pi i) in D_k genau auf)")

# ---------------------------------------------------------------------------
hr("G  1D-Toy-Model:  x'' = -GM/x^2 - gamma * int_{-inf}^t K(t-t') x'(t') dt'")
# ---------------------------------------------------------------------------
p("Einheiten: GM = 1, Start x0 = 1 (vorher Ruhe), Startgeschwindigkeit v_s > v_esc = sqrt(2).")
p("Kern K(s) = exp(-s/tau_K);  y(t) := int K x' dt'  =>  y' = x' - y/tau_K,  x'' = -1/x^2 - gamma y")
p("tau_K -> inf: K = 1, y = x - x0  (Gedaechtnis ohne Vergessen)")
p("Erweiterte Energie H = x'^2/2 - 1/x + gamma y^2/2,   dH/dt = -gamma y^2/tau_K")
p()


def simulate(gamma, tauk, vs=1.6, tmax=5000.0, dt0=2e-3):
    """RK4 mit Kepler-skalierter Schrittweite; Abbruch bei Rueckkehr nach x0 = 1
    (nach Umkehr) oder bei t > tmax."""
    x, v, y = 1.0, vs, 0.0
    inv = 0.0 if tauk is None else 1.0 / tauk

    def f(s):
        xx, vv, yy = s
        return (vv, -1.0 / xx ** 2 - gamma * yy, vv - yy * inv)

    s = (x, v, y)
    xmax = x
    turned = False
    H0 = 0.5 * v * v - 1 / x
    t = 0.0
    while t < tmax:
        dt = dt0 * min(max(s[0], 1.0) ** 1.5, 50.0)
        if tauk is not None:
            dt = min(dt, 0.05 * tauk)
        k1 = f(s)
        k2 = f(tuple(s[i] + 0.5 * dt * k1[i] for i in range(3)))
        k3 = f(tuple(s[i] + 0.5 * dt * k2[i] for i in range(3)))
        k4 = f(tuple(s[i] + dt * k3[i] for i in range(3)))
        s = tuple(s[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(3))
        t += dt
        if s[0] > xmax:
            xmax = s[0]
        if s[1] < 0:
            turned = True
        if turned and s[0] <= 1.0:
            break
    xx, vv, yy = s
    Emech = 0.5 * vv * vv - 1 / xx
    H = Emech + 0.5 * gamma * yy * yy
    return xmax, turned, Emech, H, H0, t


p("Abbruch bei Rueckkehr nach x0 = 1 bzw. bei t = 5000 (ohne Rueckkehr).")
p(f"{'gamma':>6} {'tau_K':>7} {'x_max':>10} {'Rueckkehr':>9} {'t_Ende':>9} {'E_mech':>10} {'H':>10} {'H(0)':>8}")
vs = 1.6
for gamma, tauk in [(0.0, None), (0.01, 1.0), (0.01, 10.0), (0.01, 100.0), (0.01, None),
                    (0.05, 1.0), (0.05, 10.0), (0.05, None)]:
    xmax, turned, Em, H, H0, tend = simulate(gamma, tauk, vs=vs)
    tk = "inf" if tauk is None else f"{tauk:g}"
    p(f"{gamma:6.2f} {tk:>7} {xmax:10.3f} {'ja' if turned else 'nein':>9} {tend:9.1f} "
      f"{Em:10.5f} {H:10.5f} {H0:8.4f}")
p()
p(f"Startgeschwindigkeit v_s = {vs} > sqrt(2) = {math.sqrt(2):.4f}: ohne B (gamma = 0) Flucht.")
p("tau_K = inf: H bleibt erhalten (reine Bindung, kein Energieverlust), Bahn kehrt um.")
p("endliches tau_K: H nimmt um int gamma y^2/tau_K dt ab (Reibungsanteil der Memory-Wirkung).")

# Grenzfall kurzer Kern: int exp(-s/tau) x' ~ tau x'  -> Reibung gamma*tau
p()
p("Grenzfaelle des Kerns:")
p("  tau_K -> 0 :  gamma int K x' dt' -> gamma tau_K x'           (lokale Reibung)")
p("  tau_K -> inf: gamma int K x' dt' -> gamma (x - x_(-inf))       (Rueckstellung, konservativ)")
p(f"  Bindungsradius fuer tau_K = inf, gamma = 0.01, v_s = {vs}:")
# Umkehrpunkt aus H-Erhaltung: v_s^2/2 - 1 = -1/x + gamma (x-1)^2/2
g_ = 0.01
lo, hi = 1.0, 1e4
Hs = 0.5 * vs * vs - 1.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if -1 / mid + 0.5 * g_ * (mid - 1) ** 2 < Hs:
        lo = mid
    else:
        hi = mid
p(f"    x_Umkehr (analytisch aus H) = {lo:.4f}")

# ---------------------------------------------------------------------------
hr("H  FTM-Anschluss: g = x/y^2 * g_Erde  und  Delta_s = GM/(3c^2)")
# ---------------------------------------------------------------------------
ds_b = G * M_B / (3 * c ** 2)
p(f"Delta_s(M_bary)            = {ds_b:.4e} m = {ds_b/AU:.2f} AU  (unabhaengig vom Radius)")
p(f"Delta_s(M_sun)             = {G*M_SUN/(3*c**2):.2f} m")
ratio = V0 ** 2 / (3 * c ** 2)
p(f"Linear wachsende Quelle M_DM(r) ~ v0^2 r/G  ->  Delta_s_DM(r)/r = v0^2/(3c^2) = {ratio:.4e}")
for rk in [10, 50, 100, 200]:
    p(f"   r = {rk:4d} kpc: Delta_s_DM ~ {ratio*rk*KPC:.4e} m = {ratio*rk*KPC/AU:.1f} AU")
p()
p("Hull-Beschleunigung g (FTM) entspricht a_Newton = GM/r^2; in g_Erde-Einheiten:")
for rk in [1, 8.2, 50, 100]:
    r = rk * KPC
    aN = G * M_B / r ** 2
    ax = gc1 / r
    p(f"   r = {rk:6.1f} kpc: g = {aN/G_EARTH:.4e} g_E,  a_extra = {ax/G_EARTH:.4e} g_E,"
      f"  a_extra/g = {ax/aN:.3f}")
p("Fernfeld: g ~ 1/r^2 -> 0 (Hull-Steilheit verschwindet), a_extra ~ 1/r bleibt;")
p("Delta_s(M_bary) bleibt konstant, Delta_s_DM(r) waechst ~ r.")

# ---------------------------------------------------------------------------
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "zustandsraum.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n".join(OUT))
