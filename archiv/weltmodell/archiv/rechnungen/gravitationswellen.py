# -*- coding: utf-8 -*-
"""
FTM - Gravitationswellen und Starkfeld-Regime
Rechnungen zu forschung/07-gravitationswellen-starkfeld.md. Nur Python-Standardbibliothek.
Aufruf:  python gravitationswellen.py > gravitationswellen.txt

Kernformeln des Modells:
  Ds        = GM/(3c^2)                  Faltungstiefe
  g         = 3c^2 Ds / r^2
  dtau/dt   = sqrt(1 - 6 Ds/r),  r0 = 6 Ds   (Takt = 0)
  r*        = sqrt(3c^2 Ds / a0)            Schwellenradius
  nu        = 1/(1-q),  q = exp(-sqrt(a_N/a0)) = exp(-r*/r)
  a0        = c H0 / 2pi                    (Atmen, Vorschlag Notiz 01)
"""
import math, cmath, sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ---------------------------------------------------------------- Konstanten
G     = 6.6743e-11             # m^3 kg^-1 s^-2 (CODATA 2018)
c     = 299792458.0            # m/s
GM_S  = 1.32712440018e20       # m^3/s^2 (IAU, heliozentrische Gravitationskonstante)
M_S   = GM_S / G               # kg
pc    = 3.0856775814913673e16  # m
kpc, Mpc, Gpc = 1e3*pc, 1e6*pc, 1e9*pc
AU    = 1.495978707e11
JAHR  = 365.25*86400.0
a0    = float(__import__("os").environ.get("FTM_A0", "1.2e-10"))                # m/s^2, gemessene Schwelle
eV    = 1.602176634e-19
hbar  = 1.054571817e-34
H0_70 = 70e3/Mpc               # s^-1
muas  = math.pi/180/3600*1e-6  # rad

def Ds(M_sun):            return M_sun*GM_S/(3*c*c)          # m
def GMc3(M_sun):          return M_sun*GM_S/c**3              # s  (= 3 Ds / c)
def log10q(aN, a=a0):     return -math.sqrt(aN/a)/math.log(10) # log10 q
def fmt_log(l):
    """q = 10^l fuer sehr kleine q lesbar"""
    if l > -300: return f"{10**l:.3e}"
    return f"10^({l:.4e})"

def T(title):
    print("\n" + "="*100); print(title); print("="*100)

print("FTM - Gravitationswellen und Starkfeld (alle Zahlen dieser Datei stammen aus diesem Skript)")
print(f"Konstanten: G = {G}, c = {c}, GM_sun = {GM_S:.11e} m^3/s^2, a0 = {a0} m/s^2")
print(f"Ds_sun = GM_sun/(3c^2) = {Ds(1):.4f} m,  GM_sun/c^2 = 3 Ds_sun = {3*Ds(1):.4f} m,  GM_sun/c^3 = 3 Ds_sun/c = {GMc3(1):.6e} s")
for H in (67.4, 70.0, 73.0):
    print(f"   a0 aus Atmen: c*H0/2pi (H0 = {H}) = {c*H*1e3/Mpc/(2*math.pi):.4e} m/s^2")

# =====================================================================================
T("0) Wellengleichung des Faltungsfeldes: Pruefungen der Linearisierung")
# h_00 = 6 Ds / r,  h_ij = 6 Ds / r delta_ij,  hbar_00 = 12 Ds / r ;  Box hbar = -48 pi sigma,  sigma = G T/(3 c^4)
print("Statischer Grenzfall: hbar_00 = 12 Ds/r, h_00 = h_ii = 6 Ds/r  ->  dtau/dt = sqrt(1-h_00) ~ 1 - 3Ds/r")
for name, M, r in [("Sonne, 1 AE", 1, AU), ("Sonne, Oberflaeche", 1, 6.957e8), ("Neutronenstern 1.4, 12 km", 1.4, 1.2e4)]:
    d = Ds(M); h00 = 6*d/r
    print(f"   {name:26s} h_00 = {h00:.6e}   1-sqrt(1-h_00) = {1-math.sqrt(1-h00):.6e}   3Ds/r = {3*d/r:.6e}")
# Poisson-Pruefung numerisch: Laplace von 12 Ds/r in Kugelschale = 0 ausserhalb, Fluss = -48 pi Ds
d = Ds(1); R = 1e9
flux = 4*math.pi*R*R*(-12*d/R**2)
print(f"   Gauss-Fluss von grad(12 Ds/r) durch Kugel: {flux:.6e} m ;  -48 pi Ds = {-48*math.pi*d:.6e} m  (Quelle 48 pi sigma, sigma = G rho/(3c^2))")
alpha_sun = 12*Ds(1)/6.957e8
print(f"Lichtablenkung am Sonnenrand aus h_00 + h_ij: alpha = 12 Ds/b = {alpha_sun:.6e} rad = {alpha_sun/muas/1e6:.4f} arcsec")
print(f"   nur Takt (h_00): 6 Ds/b = {alpha_sun/2/muas/1e6:.4f} arcsec  (Hälfte, Wert von 1911)")
print("Freiheitsgrade: 10 (symm. Tensor) - 4 (Lorenz-Eichung) - 4 (Resteichung) = 2  -> h_+, h_x")
print(f"Quadrupolfaktor: P = (9/5) (c^5/G) <(dddQ_F/c^3)^2>,  c^5/G = {c**5/G:.4e} W")

# =====================================================================================
T("1) Operator im Starkfeld: q = exp(-sqrt(a/a0)) = exp(-r*/r)")
def kerr_isco(j):
    Z1 = 1 + (1-j*j)**(1/3)*((1+j)**(1/3)+(1-j)**(1/3))
    Z2 = math.sqrt(3*j*j+Z1*Z1)
    return 3 + Z2 - math.sqrt((3-Z1)*(3+Z1+2*Z2))   # in GM/c^2, prograd
def r_orb_from_fgw(Mtot, fgw):
    w = math.pi*fgw; return (Mtot*GM_S/w**2)**(1/3)

rows = []
M1, M2, Mf, jf = 36.0, 29.0, 62.0, 0.67
Mt = M1+M2
print(f"GW150914: m1 = {M1}, m2 = {M2}, M = {Mt}, M_f = {Mf}, j_f = {jf} (Quellsystem)")
print(f"   Ds(M) = {Ds(Mt):.1f} m,  Ds(M_f) = {Ds(Mf):.1f} m,  r0 = 6 Ds(M_f) = {6*Ds(Mf)/1e3:.2f} km")
rstar = math.sqrt(Mf*GM_S/a0)
print(f"   r*(M_f) = sqrt(GM/a0) = {rstar:.4e} m = {rstar/AU:.1f} AE = {rstar/pc:.4f} pc")
def line(name, GM, r):
    aN = GM/r**2
    l = log10q(aN)
    print(f"   {name:44s} r = {r/1e3:12.2f} km  a_N = {aN:.3e} m/s^2  a_N/a0 = {aN/a0:.3e}  r*/r = {math.sqrt(aN/a0):.4e}  q = nu-1 = {fmt_log(l)}")
GMf, GMt = Mf*GM_S, Mt*GM_S
rp = (1+math.sqrt(1-jf*jf))*GMf/c**2
line("Horizont Schwarzschild r0 = 6 Ds(M_f)", GMf, 2*GMf/c**2)
line(f"Horizont Kerr r+ (j={jf}) = {rp*c*c/GMf:.4f} GM/c^2", GMf, rp)
line("Photonensphaere 9 Ds (Schwarzschild)", GMf, 3*GMf/c**2)
line("ISCO Schwarzschild 18 Ds(M)", GMt, 6*GMt/c**2)
print(f"   ISCO Kerr prograd, j = {jf}: r_ISCO = {kerr_isco(jf):.4f} GM/c^2 = {kerr_isco(jf)*3:.3f} Ds")
line(f"ISCO Kerr (j={jf}, M_f)", GMf, kerr_isco(jf)*GMf/c**2)
for f in (20, 35, 75, 150):
    line(f"Inspiral, f_gw = {f} Hz (Kepler-Abstand)", GMt, r_orb_from_fgw(Mt, f))
fisco = c**3/(6**1.5*math.pi*GMt)
print(f"   f_gw(ISCO, Schwarzschild) = c/(pi (18Ds)^(3/2)/sqrt(3Ds)) = {fisco:.1f} Hz  (Quellsystem)")

print("\nGW170817: m1 = 1.46, m2 = 1.27 M_sun, R_NS = 12 km")
m1, m2 = 1.46, 1.27; Mbns = m1+m2
line("NS-Oberflaeche (1.46 M_sun, 12 km)", m1*GM_S, 1.2e4)
line("Kontakt, Abstand 24 km", Mbns*GM_S, 2.4e4)
line("ISCO 18 Ds(M)", Mbns*GM_S, 6*Mbns*GM_S/c**2)
for f in (24, 100, 400):
    line(f"Inspiral, f_gw = {f} Hz", Mbns*GM_S, r_orb_from_fgw(Mbns, f))
print(f"   f_gw(ISCO) = {c**3/(6**1.5*math.pi*Mbns*GM_S):.0f} Hz")

print("\nPSR B1913+16 (Hulse-Taylor): mittlere grosse Halbachse")
mp, mc_, Pb, e = 1.438, 1.390, 27906.9795, 0.6171340
Mht = mp+mc_; aht = (Mht*GM_S*(Pb/(2*math.pi))**2)**(1/3)
line("grosse Halbachse", Mht*GM_S, aht)
print(f"\nDetektor LIGO im Erdfeld: a_N = 9.81 m/s^2 -> log10 q = {log10q(9.81):.4e}")
print(f"LISA im Sonnenfeld (1 AE): a_N = {GM_S/AU**2:.4e} -> log10 q = {log10q(GM_S/AU**2):.4e}")
aGW = 0.5*(2*math.pi*150)**2*1e-21*4e3
print(f"Gezeitenbeschleunigung der Welle auf LIGO-Spiegel (f=150 Hz, h=1e-21, L=4 km): {aGW:.3e} m/s^2 = {aGW/a0:.3e} a0")
qw = math.exp(-math.sqrt(aGW/a0))
print(f"   haette die Welle allein das Argument von nu: q = {qw:.4f}, nu = {1/(1-qw):.3f}  (tatsaechlich: Argument = Gesamtfeld |grad Phi_N|)")
print("Ordnung: q = exp(-eps^(-1/2)), eps = a0/a;  d^n q/d eps^n -> 0 fuer eps -> 0 fuer alle n  (jenseits aller PN-Ordnungen)")
for eps in (1e-2, 1e-4, 1e-6):
    print(f"   eps = {eps:.0e}:  q = {fmt_log(-eps**-0.5/math.log(10))}   vgl. eps^10 = {eps**10:.1e}")

# =====================================================================================
T("2) Ringdown: Eigenschwingungen des Faltungstrichters")
# (a) Leaver-Kettenbruch Schwarzschild, s=2, l=2, Einheiten 2M = 1
def leaver(w, l=2, s=2, N=400):
    eps = s*s-1
    def al(n): return n*n + (-2j*w+2)*n - 2j*w + 1
    def be(n): return -(2*n*n + (-8j*w+2)*n - 8*w*w - 4j*w + l*(l+1) - eps)
    def ga(n): return n*n - 4j*w*n - 4*w*w - eps - 1
    x = 0j
    for n in range(N, 0, -1):
        x = ga(n)/(be(n) - al(n)*x)
    return be(0) - al(0)*x
def secant(f, x0, x1, it=80):
    f0, f1 = f(x0), f(x1)
    for _ in range(it):
        x2 = x1 - f1*(x1-x0)/(f1-f0)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(x1-x0) < 1e-14: break
    return x1
w2M = secant(leaver, 0.74-0.18j, 0.75-0.17j)
Mw_S = w2M/2
print(f"(a) Leaver-Kettenbruch, Regge-Wheeler-Gleichung mit Takt (dtau/dt)^2 = 1-6Ds/r:")
print(f"    M omega (l=m=2, n=0, j=0) = {Mw_S.real:.6f} {Mw_S.imag:+.6f} i   (Literatur: 0.373672 - 0.088962 i)")
# (b) Potential-Maximum und WKB 1. Ordnung (Schutz & Will 1985), Einheiten GM/c^2 = 1
l = 2
V  = lambda r: (1-2/r)*(l*(l+1)/r**2 - 6/r**3)
# Maximum
lo, hi = 2.5, 4.5
for _ in range(200):
    m1_, m2_ = lo+(hi-lo)/3, hi-(hi-lo)/3
    if V(m1_) < V(m2_): lo = m1_
    else: hi = m2_
rm = (lo+hi)/2; h = 1e-4
Vrr = (V(rm+h)-2*V(rm)+V(rm-h))/h**2
V0 = V(rm); Vss = (1-2/rm)**2*Vrr
w_wkb = cmath.sqrt(V0 - 0.5j*math.sqrt(-2*Vss))
print(f"(b) Potentialmaximum l=2: r_max = {rm:.4f} GM/c^2 = {3*rm:.3f} Ds (Photonensphaere: 9 Ds)")
print(f"    WKB 1. Ordnung (Schutz & Will): M omega = {w_wkb.real:.4f} {-abs(w_wkb.imag):+.4f} i")
# (c) Eikonal: Photonensphaere
Om = 1/math.sqrt(27)
print(f"(c) Eikonal: Omega_ph = c/(3 sqrt27 Ds) -> M Omega = 1/sqrt27 = {Om:.4f}; Lyapunov lambda = Omega_ph")
print(f"    M omega ~ (l+1/2) Omega - i (n+1/2) lambda = {2.5*Om:.4f} - {0.5*Om:.4f} i  (l=2)")
for L in (2, 3, 5, 10):
    print(f"    l = {L:2d}: Eikonal Re = {(L+0.5)*Om:.4f}")
# (d) Berti-Cardoso-Will 2006 Fit (l=m=2, n=0)
f1, f2, f3 = 1.5251, -1.1568, 0.1292
q1, q2, q3 = 0.7000, 1.4187, -0.4990
def F(j): return f1 + f2*(1-j)**f3
def Q(j): return q1 + q2*(1-j)**q3
print(f"(d) BCW-Fit: M omega_R = f1 + f2 (1-j)^f3, Q = q1 + q2 (1-j)^q3 ;  f = (f1,f2,f3) = ({f1},{f2},{f3}), q = ({q1},{q2},{q3})")
print(f"    Probe j=0: F = {F(0):.4f} (Leaver {Mw_S.real:.4f}),  Q = {Q(0):.4f} (Leaver {Mw_S.real/(2*abs(Mw_S.imag)):.4f})")
for j in (0.0, 0.5, 0.67, 0.9, 0.99):
    print(f"    j = {j:4.2f}: F(j) = M omega_R = {F(j):.4f},  Q = {Q(j):.4f},  M omega_I = F/(2Q) = {F(j)/(2*Q(j)):.4f}")
z150914 = 0.09
print(f"\nGW150914: M_f = {Mf} M_sun, j = {jf}, z = {z150914}")
d_f = Ds(Mf)
tM = 3*d_f/c
fq = c/(2*math.pi*3*d_f)*F(jf)
tau = Q(jf)/(math.pi*fq)
print(f"   Ds(M_f) = {d_f:.1f} m,  3 Ds = GM/c^2 = {3*d_f/1e3:.3f} km,  GM/c^3 = 3 Ds/c = {tM*1e3:.5f} ms")
print(f"   f_QNM = c/(2 pi 3Ds) F(j) = {c/(2*math.pi*3*d_f):.2f} Hz x {F(jf):.4f} = {fq:.1f} Hz (Quellsystem)")
print(f"   tau = Q/(pi f) = {tau*1e3:.3f} ms (Quellsystem)")
print(f"   Detektorsystem (x 1/(1+z) bzw. x (1+z)): f = {fq/(1+z150914):.1f} Hz,  tau = {tau*(1+z150914)*1e3:.3f} ms")
print(f"   Wellenlaenge lambda = c/f = {c/fq/1e3:.0f} km = {c/fq/d_f:.2f} Ds = {c/fq/(6*d_f):.2f} r0")
print(f"   Dämpfung in Takt-Einheiten: tau / (3Ds/c) = {tau/tM:.2f}")
# Schwarzschild-Wert zum Vergleich
print(f"   zum Vergleich j = 0: f = {c/(2*math.pi*3*d_f)*Mw_S.real:.1f} Hz, tau = {1/abs(Mw_S.imag)*tM*1e3:.3f} ms")
# Unsicherheit M_f = 62 +- 4
for Mx in (58, 62, 66):
    dx = Ds(Mx); fx = c/(2*math.pi*3*dx)*F(jf)/(1+z150914)
    print(f"   M_f = {Mx}: f_det = {fx:.1f} Hz, tau_det = {Q(jf)/(math.pi*fx)*1e3:.2f} ms")

# =====================================================================================
T("2b) Quadrupolformel und Chirp in Faltungsgroessen")
def eta(m1, m2): return m1*m2/(m1+m2)**2
def Mchirp(m1, m2): return (m1*m2)**0.6/(m1+m2)**0.2
for name, a, b, z, DL in [("GW150914", 36., 29., 0.09, 410*Mpc), ("GW170817", 1.46, 1.27, 0.0099, 40*Mpc)]:
    Mc = Mchirp(a, b); Mcz = Mc*(1+z)
    dsc, dscz = Ds(Mc), Ds(Mcz)
    print(f"{name}: M_c = {Mc:.3f} M_sun (Quelle), (1+z) M_c = {Mcz:.3f} M_sun (Detektor), eta = {eta(a,b):.4f}")
    print(f"   Chirp-Faltung Ds_c = G M_c/(3c^2) = {dsc:.1f} m (Quelle), {dscz:.1f} m (Detektor)")
    tc = 3*dscz/c
    fl = (20., 35., 150.) if name == "GW150914" else (24., 100., 400.)
    for f in fl:
        fdot = 96/5*math.pi**(8/3)*tc**(5/3)*f**(11/3)
        tau_c = 5/256*tc**(-5/3)*(math.pi*f)**(-8/3)
        Ncyc = 1/(32*math.pi**(8/3))*tc**(-5/3)*f**(-5/3)
        h0 = 4/DL*(3*dscz)**(5/3)*(math.pi*f/c)**(2/3)
        # Rueckrechnung Ds_c aus f, fdot
        ds_back = c/3*((5/96)*math.pi**(-8/3)*f**(-11/3)*fdot)**0.6
        print(f"   f = {f:5.0f} Hz: fdot = {fdot:.4e} Hz/s, t bis Verschmelzung = {tau_c:.4f} s, Zyklen = {Ncyc:.0f}, h0(face-on) = {h0:.3e}, Ds_c aus (f,fdot) = {ds_back:.1f} m")
    # Leistung bei r = 18 Ds (ISCO) und 9 Ds
    Mt_ = a+b
    for k, lab in ((18, "ISCO 18 Ds"), (12, "r = 12 Ds"), (9, "r = 9 Ds")):
        delta = 3/k
        P = 32/5*c**5/G*eta(a,b)**2*delta**5
        print(f"   Leistung P = (32/5)(c^5/G) eta^2 (3Ds/r)^5 bei {lab:10s}: (3Ds/r) = {delta:.4f}, P = {P:.3e} W")
print("Vergleich: GW150914 Spitzenleistung (LIGO) 3.6e49 W;  Energie 3.0 M_sun c^2 =", f"{3.0*M_S*c*c:.3e} J")
# Hulse-Taylor
fe = (1+73/24*e**2+37/96*e**4)/(1-e*e)**3.5
Mc_ht = Mchirp(mp, mc_)
Pbdot = -192*math.pi/5*(2*math.pi*3*Ds(Mc_ht)/c/Pb)**(5/3)*fe
print(f"PSR B1913+16: Ds_c = {Ds(Mc_ht):.2f} m, f(e) = {fe:.4f}, Pbdot = -(192pi/5)(2pi 3Ds_c/(c Pb))^(5/3) f(e) = {Pbdot:.4e} (Literatur GR: -2.40263e-12)")
print(f"Faltungsbilanz GW150914: Ds(65) - Ds(62) = {Ds(65)-Ds(62):.1f} m = Ds(3 M_sun) = {Ds(3):.1f} m")

# =====================================================================================
T("3) Ausbreitung: GW170817 / GRB 170817A")
dt, D_low, D_mid = 1.74, 26*Mpc, 40*Mpc
print(f"dt = 1.74 s, D >= 26 Mpc:  c dt / D = {c*dt/D_low:.3e}  (obere Grenze fuer c_gw > c)")
print(f"   mit angenommener Eigenverzögerung bis 10 s: (10-1.74) c/D = {(10-1.74)*c/D_low:.3e}")
print(f"   bei D = 40 Mpc: c dt/D = {c*dt/D_mid:.3e}")
print("   Literaturwert (Abbott et al. 2017, ApJL 848 L13): -3e-15 <= dv/c <= +7e-16")
# Graviton-Masse durch Atmen: m c^2 = hbar H0/(2pi) ?
wa = H0_70/(2*math.pi)     # = a0/c bei a0 = cH0/2pi
print(f"Frequenzskala des Atmens: H0/2pi = a0/c = {wa:.3e} s^-1  (a0/c mit a0 = 1.2e-10: {a0/c:.3e} s^-1)")
print(f"   als Massenskala: hbar H0/2pi = {hbar*wa/eV:.3e} eV   (LIGO-Grenze: m_g < 1.2e-22 eV (2016), 1.76e-23 eV (GWTC-3))")
for f in (100., 1e-3, 1e-8):
    om = 2*math.pi*f
    print(f"   Dispersion bei f = {f:.0e} Hz: |v_g - c|/c = (1/2)(omega_a/omega)^2 = {0.5*(wa/om)**2:.3e}")
# Shapiro-Unterschied, falls nur Licht (nicht GW) den nu-Anteil saehe
for vf_kms in (175.8, 230.0):
    vf = vf_kms*1e3
    R0, re = 8.2*kpc, 200*kpc
    I = re - R0 - R0*math.log(re/R0)
    dtS = 2*(vf/c)**2*I/c
    print(f"Shapiro-Anteil des logarithmischen Trichters (v_flach = {vf_kms} km/s, 8.2 -> 200 kpc): {dtS:.3e} s = {dtS/86400:.0f} Tage")
MW_Ds = Ds(6e10)
print(f"Newton-Anteil Milchstrasse 6 Ds/c = {6*MW_Ds/c:.3e} s = {6*MW_Ds/c/86400:.1f} Tage (mal Log-Faktor)")

# =====================================================================================
T("4) Gedaechtnis des Raumes = GW-Memory (Christodoulou)")
# h_mem = (2 eta M x / R) sin^2 Th (17 + cos^2 Th)/96 ; E_rad = eta M x/2 -> h_mem = (4 E_rad/R) k(Th)
def kap(th): return math.sin(th)**2*(17+math.cos(th)**2)/96
print(f"Winkelfaktor k(Theta) = sin^2(17+cos^2)/96: Kante (90 Grad) = {kap(math.pi/2):.4f}, 150 Grad = {kap(math.radians(150)):.4f}, 0 Grad = {kap(0):.4f}")
print("In Faltungsgroessen: E_rad = 3 c^4 Ds_rad / G  ->  h_mem = 12 k(Theta) Ds_rad / D_L  (Kante: 2.125 Ds_rad/D_L)")
z, DL = 0.09, 410*Mpc
Erad = 3.0
dsr = Ds(Erad*(1+z))
for th in (90, 150):
    hm = 12*kap(math.radians(th))*dsr/DL
    print(f"GW150914 (E_rad = 3.0 M_sun, D_L = 410 Mpc), Theta = {th}: Ds_rad(1+z) = {dsr:.1f} m, h_mem = {hm:.3e}, Anteil an h_peak = 1.0e-21: {hm/1e-21:.3f}")
hm90 = 12*kap(math.pi/2)*dsr/DL
print(f"   bleibende Armverschiebung LIGO (L = 4 km): dL = h/2 L = {hm90/2*4e3:.2e} m (Kante)")
print(f"   Verhaeltnis Memory/oszillierend bei gleichem x (Kante, fuehrende Ordnung): 17/96 = {17/96:.3f}")
xeff = 2*(Erad/65)/eta(36,29)
print(f"   effektives x = 2 E_rad/(eta M) = {xeff:.3f}")
# Stacking (Lasky et al. 2016): SNR 3 mit ~35, SNR 5 mit ~90 Ereignissen
for N, snr in ((35, 3), (90, 5)):
    print(f"   Lasky et al.: SNR {snr} mit ~{N} Ereignissen -> SNR pro Ereignis ~ {snr/math.sqrt(N):.2f}")
# GW170817 memory
x_c = (2.73*GM_S/c**2)/2.4e4
Ebns = eta(1.46,1.27)*2.73*x_c/2
hb = 12*kap(math.pi/2)*Ds(Ebns)/(40*Mpc)
print(f"GW170817: x(Kontakt 24 km) = {x_c:.3f}, E_rad(Inspiral) ~ {Ebns:.3f} M_sun, h_mem(Kante, 40 Mpc) ~ {hb:.2e}")
# LISA, PTA
def DL_flat(z, H0=67.7e3/Mpc, Om=0.31, N=4000):
    s = 0; dz = z/N
    for i in range(N):
        zz = (i+0.5)*dz; s += dz/math.sqrt(Om*(1+zz)**3+1-Om)
    return (1+z)*c/H0*s
for name, M, z in (("LISA: SMBH 1e6 M_sun (gleich), z = 1", 1e6, 1.0), ("LISA: SMBH 1e7 M_sun, z = 3", 1e7, 3.0),
                   ("PTA: SMBH 1e10 M_sun, z = 0.2", 1e10, 0.2), ("PTA: SMBH 1e9 M_sun, z = 0.05", 1e9, 0.05)):
    D = DL_flat(z); E = 0.048*M
    hm = 12*kap(math.pi/2)*Ds(E*(1+z))/D
    trise = 3*Ds(M*(1+z))/c*100
    print(f"{name:40s} D_L = {D/Gpc:.2f} Gpc, E_rad = 4.8% M, h_mem(Kante) = {hm:.2e}, Anstiegszeit ~100 GM_z/c^3 = {trise:.2e} s")
    if name.startswith("PTA"):
        Tobs = 15*JAHR
        print(f"      Timing-Rest ~ h_mem * T_obs /2 (T = 15 a) = {hm*Tobs/2*1e9:.0f} ns")

# =====================================================================================
T("5) Tief-a0-Regime: wo das FTM abweichen koennte")
print("(a) Schwellenfrequenz eines Binaersystems: Bahnbeschleunigung = a0 bei r = r*:")
print("    f*_gw = (1/pi) sqrt(a0/r*) = a0^(3/4) / (pi (3 c^2 Ds)^(1/4))")
for M in (1, 60, 1e6, 1e9, 1e10):
    GM = M*GM_S; fs = a0**0.75/(math.pi*GM**0.25)
    rs = math.sqrt(GM/a0)
    print(f"    M = {M:8.0e} M_sun: r* = {rs/pc:10.4f} pc, f*_gw = {fs:.3e} Hz, Umlaufzeit = {2/fs/JAHR:.3e} a")
print("    Baender: PTA 1e-9..1e-7 Hz, LISA 1e-4..1e-1 Hz, LIGO 10..1e3 Hz, CMB-B-Moden ~1e-18..1e-16 Hz")
for M, f in ((1e9, 1e-8), (1e10, 1e-9), (1e9, 1e-9)):
    GM = M*GM_S; r = r_orb_from_fgw(M, f); aN = GM/r**2
    print(f"    PTA-Quelle M = {M:.0e}, f_gw = {f:.0e} Hz: Abstand = {r/AU:.0f} AE, a = {aN:.3e} m/s^2 = {aN/a0:.2e} a0, log10 q = {log10q(aN):.3e}")
print("(b) Galaktisches Feld am Ort der Sonne (Pulsar-Beschleunigungen):")
aN_sun = 6e10*GM_S/(8.2*kpc)**2
qs = math.exp(-math.sqrt(aN_sun/a0)); nus = 1/(1-qs)
print(f"    a_N(8.2 kpc, 6e10 M_sun) = {aN_sun:.4e} m/s^2 = {aN_sun/a0:.3f} a0, q = {qs:.4f}, nu = {nus:.4f}, a = nu a_N = {nus*aN_sun:.4e} m/s^2")
v233 = 233e3
print(f"    Vergleich v_c = 233 km/s: v^2/R = {v233**2/(8.2*kpc):.4e} m/s^2")
print(f"    Periodenaenderung Pdot/P = a_los/c: Newton {aN_sun/c:.3e} s^-1, FTM {nus*aN_sun/c:.3e} s^-1")
Sb = 50.0  # M_sun/pc^2 angenommene baryonische Scheibe
print(f"    Vertikal (fuehrende Ordnung): Sigma_eff = nu Sigma_b; fuer angenommenes Sigma_b = {Sb} M_sun/pc^2 -> {nus*Sb:.0f} M_sun/pc^2")
print("(c) Wellen durch kollektive Trichter (Linseneffekt, gleich fuer Licht und GW):")
for name, M in (("Milchstrasse 6e10", 6e10), ("Galaxie 2e11", 2e11), ("Haufen 1.7e14", 1.7e14)):
    GM = M*GM_S; vf2 = math.sqrt(GM*a0)
    al = 2*math.pi*vf2/c**2
    print(f"    {name:18s}: v_flach = {math.sqrt(vf2)/1e3:7.1f} km/s, alpha_inf = 2 pi (v/c)^2 = 6 pi Ds/r* = {al:.3e} rad = {al/muas/1e6:.3f} arcsec")
    # Probe 6 pi Ds / r*
    print(f"        Probe 6 pi Ds/r* = {6*math.pi*(GM/(3*c*c))/math.sqrt(GM/a0):.3e} rad")
    for b in (100*kpc,):
        Meff = b*vf2/G
        fwo = c**3/(8*math.pi*G*Meff)
        print(f"        M_eff(<{b/kpc:.0f} kpc) = b v^2/G = {Meff/M_S:.3e} M_sun, Wellenoptik unterhalb f ~ c^3/(8 pi G M_eff) = {fwo:.2e} Hz")
print("(d) Atmen als Hintergrund:")
for H in (67.4, 70.0, 73.0):
    Hs = H*1e3/Mpc
    print(f"    H0 = {H}: H0 = {Hs:.4e} s^-1, H0/2pi = {Hs/(2*math.pi):.4e} Hz, 1/H0 = {1/Hs/JAHR/1e9:.2f} Gyr")
print(f"    Omega_GW(f) = (2 pi^2/3H0^2) f^2 h_c^2 ; bei f = H0/2pi: Omega = h_c^2/6 = {(2*math.pi**2/3)*(1/(2*math.pi))**2:.6f} h_c^2")
print(f"    Verhaeltnis f_LIGO(100 Hz)/(H0/2pi) = {100/(H0_70/(2*math.pi)):.2e},  f_PTA(1e-8)/(H0/2pi) = {1e-8/(H0_70/(2*math.pi)):.2e}")
for sH, lab in ((1.0e-6, "Tensor regulaer"), (4.7e-11, "Vektor")):
    print(f"    Scherung (sigma/H)0 < {sH:.1e} ({lab}): Omega_sigma = (sigma/H)^2/3 < {sH**2/3:.1e}")
print("    Korrelationsmuster PTA:  Hellings-Downs (Tensor) vs. Atmen-Mode (3 + cos z)/8")
for zdeg in (0, 30, 60, 90, 120, 180):
    zr = math.radians(zdeg); x = (1-math.cos(zr))/2
    hd = 0.5 if x == 0 else 0.5 - x/4 + 1.5*x*math.log(x)
    st = (3+math.cos(zr))/8
    print(f"    zeta = {zdeg:3d} Grad: HD = {hd:+.4f}   Atmen (ST) = {st:+.4f}")
print("    Standard-Sirenen: Xi = d_L^GW/d_L^EM = 1 (Operator ohne Reibungsterm)")

# =====================================================================================
T("6) Schatten: Photonensphaere r = 9 Ds, Schattenradius b_c = 3 sqrt27 Ds")
print(f"b(r) = r/(dtau/dt) minimal bei r = 9 Ds: b_c = 9 Ds/sqrt(1/3) = {9/math.sqrt(1/3):.4f} Ds = 3 sqrt27 Ds = {3*math.sqrt(27):.4f} Ds")
print(f"   dtau/dt an der Photonensphaere = sqrt(1/3) = {math.sqrt(1/3):.4f}")
for name, M, D, obs in (("M87*", 6.5e9, 16.8*Mpc, "Ring 42 +- 3 muas (EHT 2019)"),
                        ("Sgr A* (4.0e6, 8.2 kpc)", 4.0e6, 8.2*kpc, "Ring 51.8 +- 2.3 muas (EHT 2022)"),
                        ("Sgr A* (4.297e6, 8.277 kpc, GRAVITY)", 4.297e6, 8.277*kpc, "")):
    d = Ds(M); thg = 3*d/D
    dia = 2*3*math.sqrt(27)*d/D
    print(f"{name}: Ds = {d:.4e} m = {d/AU:.3f} AE, r0 = 6Ds = {6*d/AU:.3f} AE, Photonensphaere 9Ds = {9*d/AU:.3f} AE")
    print(f"   theta_g = 3Ds/D = {thg/muas:.3f} muas, Schattendurchmesser 2 x 3 sqrt27 Ds/D = {dia/muas:.2f} muas, Ring (alpha = 11 theta_g) = {11*thg/muas:.1f} muas   {obs}")
    aN = M*GM_S/(9*d)**2
    print(f"   an 9 Ds: a_N = {aN:.3e} m/s^2, log10 q = {log10q(aN):.4e}")
print("\nEnde.")
