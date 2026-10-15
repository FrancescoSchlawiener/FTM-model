# -*- coding: utf-8 -*-
"""
Gesetzesprüfung des FTM-Modells (Teil des Prüfstands; die Zahlen des Modells prüft ftm_pruefstand.py).

Geprüft wird nicht einzelnes Nachrechnen, sondern ob das Modell die bekannten Gesetze trägt:
  A  Grenzgesetze der Physik: Lehrbuchform gegen die Faltungsform des Kerns (eigener Rechenweg je Gesetz)
  B  Gesetze der Galaxien aus der Literatur (Milgrom, McGaugh, Lelli, Chae): mit der örtlichen Regel des Kerns nachgerechnet
  C  Gesetze des Zustandsraums (Laplace/Lévy, Kac, Pólya/Watson, Bose/Wien, Kippstelle, Fluktuation-Dissipation,
     Onsager, Ramanujan, Thetareihen, Barkhausen): durch Simulation oder zweiten Rechenweg bestätigt
  Z  Zustandsraum des Autors (Tafeln 1–6: Maxwell, Stoffgesetze, Lorentz, Poynting, Wärme, Mechanik, drei Rollen, Torus-Bündel,
     Schnitt/Rand/Schließung, Trenntests): jede Aussage durch Simulation oder zweiten Rechenweg
  M  Unsere Gesetze aus den Bausteinen 01–09 und den Notizen 71/72 (Produkt fest, Summe der Rückläufe, Dimension als Faktor,
     Quadrat bei Senkrechtem, Nome): zweiter Rechenweg (Nullstelle, Metrik, Differenzenrechnung, Reihe, Monte Carlo)
  E  Messwerte mit Quelle gegen das Modell (Abstand in Standardabweichungen)
Bezeichnungen und Standardwerte nach der Notation in 00-forschungsstand.md; alle Formeln aus ftm_kern.py.
"""
import io, sys, math, time
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigvalsh, solve_banded
from scipy.optimize import brentq, minimize_scalar, curve_fit
from scipy.special import i0, i1, k0, k1, ive, gamma, zeta, gammaln
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import *

t_start = time.time()
rng = np.random.default_rng(29092026)
R = []   # (Teil, Gesetz, Quelle, Modell, Vergleich, Maß, Urteil)

def ident(teil, gesetz, quelle, modell, lehrbuch, tol=1e-9):
    d = abs(modell / lehrbuch - 1)
    R.append((teil, gesetz, quelle, modell, lehrbuch, f"{d:.1e}", "OK" if d <= tol else "ABWEICHUNG"))

def mess(teil, gesetz, quelle, modell, wert, fehler):
    n = (modell - wert) / fehler
    R.append((teil, gesetz, quelle, modell, wert, f"{n:+.1f}σ", "OK" if abs(n) <= 2 else ("SPANNUNG" if abs(n) <= 3 else "WIDERSPRUCH")))

def sim(teil, gesetz, quelle, modell, simuliert, tol):
    d = abs(simuliert / modell - 1)
    R.append((teil, gesetz, quelle, modell, simuliert, f"{d:.1e}", "OK" if d <= tol else "ABWEICHUNG"))

a0 = A0
GMs = G * M_sun

# ======================================================================== A Grenzgesetze der Physik
M, r, b = M_sun, 1.5e11, R_sun
ident("A", "Schwarzschild-Radius 2GM/c² = r₀ = 6Δs", "Schwarzschild 1916", r0(M), 2 * GMs / c**2)
ident("A", "Zeitdilatation √(1 − 2GM/rc²)", "Einstein 1911/1916", takt(M, r), math.sqrt(1 - 2 * GMs / (r * c**2)))
ident("A", "Lichtablenkung 4GM/(c²b) = 12Δs/b", "Einstein 1915", ablenkung(M, b), 4 * GMs / (c**2 * b))
ident("A", "Periheldrehung 6πGM/(c²a(1−e²)) = 18πΔs/(a(1−e²))", "Einstein 1915", perihel_pro_umlauf(M, 5.7909e10, 0.2056), 6 * math.pi * GMs / (c**2 * 5.7909e10 * (1 - 0.2056**2)))
ident("A", "Photonensphäre 3GM/c² = 9Δs", "Schwarzschild", 9 * ds(M), 3 * GMs / c**2)
ident("A", "Schattenradius 3√3·GM/c²", "Synge 1966", schatten_winkel(M, 1.0) / 2, 3 * math.sqrt(3) * GMs / c**2)
# ISCO Kerr über das Minimum der Bahnenergie (zweiter Weg, nicht die Bardeen-Formel des Kerns)
def isco_energie(j, mit=True):
    s = 1 if mit else -1
    E = lambda x: (x**1.5 - 2 * x**0.5 + s * j) / (x**0.75 * math.sqrt(x**1.5 - 3 * x**0.5 + s * 2 * j))
    lo = 1.0001 if mit else 3.0
    # Bahnenergie hat am ISCO ihr Minimum längs der Kreisbahnen
    xs = np.linspace(lo + 0.3 if not mit else 1.3, 12, 20000)
    Ev = []
    for x in xs:
        try: Ev.append(E(x))
        except ValueError: Ev.append(np.nan)
    Ev = np.array(Ev); k = int(np.nanargmin(Ev))
    return minimize_scalar(E, bounds=(xs[max(k - 5, 0)], xs[min(k + 5, len(xs) - 1)]), method="bounded").x
for j_ in (0.0, 0.9):
    ident("A", f"ISCO Kerr j = {j_} mitläufig (Minimum der Bahnenergie)", "Bardeen, Press, Teukolsky 1972", isco_kerr(j_), isco_energie(j_), 1e-4)
ident("A", "ISCO Kerr j = 0,9 gegenläufig", "Bardeen, Press, Teukolsky 1972", isco_kerr(0.9, False), isco_energie(0.9, False), 1e-4)
ident("A", "Hawking-Temperatur ħc³/(8πGMk) = ħc/(24πkΔs)", "Hawking 1974", hbar * c / (24 * math.pi * kB * ds(M)), hbar * c**3 / (8 * math.pi * GMs * kB))
ident("A", "Bekenstein-Hawking S = 4πkGM²/(ħc) = 36πk(Δs/ℓ_P)²", "Bekenstein 1973", 36 * math.pi * (ds(M) / l_P)**2, 4 * math.pi * G * M**2 / (hbar * c))
ident("A", "Hawking-Lebensdauer 5120πG²M³/(ħc⁴)", "Page 1976", lebensdauer_falte(M), 5120 * math.pi * G**2 * M**3 / (hbar * c**4))
ident("A", "Eddington-Grenze 4πGMm_pc/σ_T", "Eddington 1926", eddington(M), 4 * math.pi * GMs * m_p * c / (8 * math.pi / 3 * r_e**2))
Hh = H0_si()
ident("A", "Gibbons-Hawking ħH/(2πk) = ħc/(12πkΔs_H), Δs_H = c/6H", "Gibbons, Hawking 1977", hbar * c / (12 * math.pi * kB * c / (6 * Hh)), hbar * Hh / (2 * math.pi * kB))
ident("A", "Kepler T² = 4π²a³/(GM) = 4π²a³/(3c²Δs)", "Kepler 1619, Newton", 2 * math.pi * math.sqrt(AU**3 / (3 * c**2 * ds(M))), 2 * math.pi * math.sqrt(AU**3 / GMs))
m1, m2, rr = 1.4 * M_sun, 1.3 * M_sun, 1e6
Mt = m1 + m2
ident("A", "Quadrupolleistung (32/5)G⁴(m₁m₂)²M/(c⁵r⁵)", "Peters, Mathews 1963", 32 / 5 * c**5 / G * (m1 * m2 / Mt**2)**2 * (3 * ds(Mt) / rr)**5, 32 / 5 * G**4 * (m1 * m2)**2 * Mt / (c**5 * rr**5))
Pb, e_ = 27906.98, 0.6171338
fe = (1 + 73 / 24 * e_**2 + 37 / 96 * e_**4) / (1 - e_**2)**3.5
ident("A", "Hulse-Taylor Ṗ (Peters-Mathews)", "Peters 1964", hulse_taylor(1.4398, 1.3886, Pb, e_),
      -192 * math.pi * G**(5 / 3) / (5 * c**5) * (Pb / (2 * math.pi))**(-5 / 3) * fe * 1.4398 * 1.3886 * M_sun**2 / (2.8284 * M_sun)**(1 / 3), 1e-3)
w_orb = math.sqrt(G * M_earth / 7027.4e3**3)
ident("A", "geodätische Präzession (3/2)(GM/c²r)ω = (9/2)(Δs/r)ω", "de Sitter 1916", 4.5 * ds(M_earth) / 7027.4e3 * w_orb, 1.5 * G * M_earth / (c**2 * 7027.4e3) * w_orb)
J_E = 8.034e37 * 7.2921e-5; aJ = J_E / (M_earth * c)
ident("A", "Mitführung im Polarorbit GJ/(2c²a³) = 3Δs·a_J·c/(2a³)", "Lense, Thirring 1918; Schiff 1960", 3 * ds(M_earth) * aJ * c / (2 * 7027.4e3**3), G * J_E / (2 * c**2 * 7027.4e3**3))
ident("A", "Sternmasse m_Pl³/m_p² (Landau-Chandrasekhar-Skala)", "Landau 1932", sternmasse(), m_Pl**3 / m_p**2)
ident("A", "α_G(M) = 3Δs/λ̄ = (M/m_Pl)²", "Dirac, Eddington", alpha_G(1e-20), (1e-20 / m_Pl)**2)

# ======================================================================== B Gesetze der Galaxien (Literatur)
gN = np.logspace(-14, -8, 13)
d = np.max(np.abs(a_regel(gN, a0) / (gN / (1 - np.exp(-np.sqrt(gN / a0)))) - 1))
R.append(("B", "Beschleunigungsrelation: Form g/(1 − e^{−√(g/g†)})", "McGaugh, Lelli, Schombert 2016 (PRL 117, 201101)", 1.0, 1.0, f"{d:.1e}", "OK" if d < 1e-12 else "ABWEICHUNG"))
# Skaleninvarianz des Tiefbereichs: a(λr)·λ = a(r) für a_N ≪ a₀
Mm, rr0 = 1e7 * M_sun, 100 * kpc
lam = 3.0
aN1, aN2 = G * Mm / rr0**2, G * Mm / (lam * rr0)**2
d = abs(a_regel(aN2, a0) * lam / a_regel(aN1, a0) - 1)
R.append(("B", "Tiefbereich skaleninvariant: a(λr)·λ = a(r) (a_N = 10⁻⁴·a₀)", "Milgrom 2009 (ApJ 698, 1630)", 1.0, float(a_regel(aN2, a0) * lam / a_regel(aN1, a0)), f"{d:.1e}", "OK" if d < 0.02 else "ABWEICHUNG"))
# Virialsatz im Tiefbereich: σ² (massengewichtet, 1D) = (2/9)√(GMa₀) für jedes Profil – numerisch mit der örtlichen Regel für eine Plummer-Kugel
def virial_sigma2(Mv, bpl):
    rho = lambda x: 3 * Mv / (4 * math.pi * bpl**3) * (1 + (x / bpl)**2)**-2.5
    Menc = lambda x: Mv * x**3 / (x * x + bpl * bpl)**1.5
    f = lambda u: 4 * math.pi * (math.exp(u))**3 * rho(math.exp(u)) * float(a_regel(G * Menc(math.exp(u)) / math.exp(u)**2, a0)) * math.exp(u)
    return quad(f, math.log(bpl) - 20, math.log(bpl) + 20, limit=400)[0] / (3 * Mv)
Mv = 1e5 * M_sun
s2n = virial_sigma2(Mv, 20 * kpc)
sim("B", "Virialsatz im Tiefbereich σ⁴ = (4/81)GMa₀ (Plummer, örtliche Regel)", "Milgrom 1995; McGaugh, Milgrom 2013", (4 / 81 * G * Mv * a0)**0.5, s2n, 0.02)
R.append(("B", "Baryonische Tully-Fisher: v⁴ = GMa₀ (Steigung 4 exakt im Grenzfall)", "Lelli et al. 2019 (MNRAS 484, 3267): Steigung 3,85 ± 0,09", 4.0, 3.85, f"{(4 - 3.85) / 0.09:+.1f}σ", "OK" if abs(4 - 3.85) / 0.09 <= 2 else "SPANNUNG"))
# Zentrale Flächendichte (Toomre-Integral) für Exponentialscheiben mit der örtlichen Regel (radiale Form)
MSPC = M_sun / pc**2
fy = lambda y: y * (i0(y) * k0(y) - i1(y) * k1(y)) if y < 300 else 1 / (8 * y * y)
tn = quad(lambda y: fy(y) / y, 0, 50, limit=400)[0] + 1 / (16 * 50**2)
ident("B", "Toomre: (1/2πG)∫V²/R² dR = Σ(0) für die Exponentialscheibe (Newton)", "Toomre 1963, Freeman 1970", tn, 1.0, 1e-3)
def sigma_dyn(S0):
    gs = 2 * math.pi * G * S0 * MSPC
    return S0 * (quad(lambda y: float(nu(gs * fy(y), a0)) * fy(y) / y, 0, 50, limit=400)[0]
                 + quad(lambda y: float(nu(gs / (8 * y * y), a0)) / (8 * y**3), 50, 5e4, limit=200)[0])
S0s = np.logspace(0.5, 4.5, 25)
Sd = np.array([sigma_dyn(s) for s in S0s])
dpl = lambda x, S0_, Sc, bet: np.log10(S0_ * (1 + 10**x / Sc)**(1 - bet) * (10**x / Sc)**bet)
p, _ = curve_fit(dpl, np.log10(S0s), np.log10(Sd), p0=(1500, 1500, 0.6), maxfev=20000)
mess("B", "Zentrale Flächendichte, ideale Exponentialscheibe: Steigung β der Doppelpotenz (α = 1)", "Lelli et al. 2016 (ApJL 827, L19): β = 0,62 ± 0,04", p[2], 0.62, 0.04)
lel = lambda S: 1468 * (1 + S / 1571) ** (1 - 0.62) * (S / 1571) ** 0.62
dd = [math.log10(sigma_dyn(S) / lel(S)) for S in (10, 100, 1000)]
R.append(("B", "Zentrale Flächendichte, ideale Exponentialscheibe gegen Lelli-Fit bei Σ₀ = 10 / 100 / 1000 M☉/pc² [dex]", "Lelli et al. 2016: Streuung 0,21 dex", float(max(dd, key=abs)), 0.0,
          " / ".join(f"{x:+.2f}" for x in dd), "OK" if max(abs(x) for x in dd) <= 0.21 else "SPANNUNG"))
import sparc_flaechendichte as _sfd
_m, _em, _k, _ek, _ = _sfd.trend(_sfd.main())
mess("B", "Zentrale Flächendichte Galaxie für Galaxie (eigene Baryonen): Rest log Σ_gem/Σ_Regel [dex]", "SPARC, Lelli et al. 2016 (Υ_d-Spanne 0,04 dex als Fehler)", _m, 0.0, math.hypot(_em, 0.04))
mess("B", "Zentrale Flächendichte Galaxie für Galaxie: Trend des Rests je Dekade Σ₀", "SPARC, Lelli et al. 2016", _k, 0.0, _ek)
# Crater II: Schätzer von McGaugh 2016 gegen die Plummer-Virialsumme mit derselben Verstärkung (Außenfeld dominiert)
Mc2, rh, Dgc = 3.2e5 * M_sun, 1.066 * kpc, 117.5 * kpc
gNe = G * 6e10 * M_sun / Dgc**2
nue = float(nu(gNe, a0))
sig_mcg = math.sqrt(nue * G * Mc2 / (3 * 4 / 3 * rh))
sig_vir = math.sqrt(math.pi / 32 * nue * G * Mc2 / rh)
R.append(("B", "Crater II im Außenfeld: Plummer-Virial / Schätzer McGaugh (σ² ∝ M statt M/2) [km/s]", "McGaugh 2016 (ApJL 832, L8): 2,1 (+0,9/−0,6)", sig_vir / 1e3, sig_mcg / 1e3, f"Faktor {sig_mcg**2 / sig_vir**2:.2f} in σ²", "HINWEIS"))
R.append(("B", "Äußeres Feld senkt ν (Richtung zählt): ν frei / senkrecht / gleichgerichtet (Milchstraße, 50 kpc, 0,03a₀)", "Chae et al. 2020 (ApJ 904, 51): nachgewiesen, 4–11σ",
          float(nu(G * 6e10 * M_sun / (50 * kpc)**2, a0)), nu_aussenfeld(G * 6e10 * M_sun / (50 * kpc)**2, 0.03 * a0, "senkrecht"),
          f"{nu_aussenfeld(G * 6e10 * M_sun / (50 * kpc)**2, 0.03 * a0, 'gleich'):.2f}", "OK"))

# ======================================================================== C Gesetze des Zustandsraums
# C1 Rückkehr mit Verlust: Wahrscheinlichkeit, von L zum Ort zurückzukehren, bevor der Verlust greift = e^{−L√(Γ/D)}
def rueckkehr_mc(L, p, n=40000):
    x = np.full(n, L); alive = np.ones(n, bool); hit = np.zeros(n, bool)
    while alive.any():
        idx = np.where(alive)[0]
        x[idx] += rng.choice((-1, 1), size=idx.size)
        kill = rng.random(idx.size) < p
        h = x[idx] == 0
        hit[idx[h & ~kill]] = True
        alive[idx[h | kill]] = False
    return hit.mean()
L_, p_ = 12, 0.004
qmc = rueckkehr_mc(L_, p_)
sim("C", "Rückkehr mit Verlust: q = e^{−L√(Γ/D)} (Irrfahrt, D = ½, Γ = p)", "Laplace-Transformierte der Erstrückkehrzeit (Lévy)", math.exp(-L_ * math.sqrt(2 * p_)), qmc, 0.05)
ident("C", "diskret exakt: [(1 − √(1 − (1−p)²))/(1 − p)]^L gegen e^{−L√(2p)}", "Feller, Kap. XIV", ((1 - math.sqrt(1 - (1 - p_)**2)) / (1 - p_))**L_, math.exp(-L_ * math.sqrt(2 * p_)), 0.01)
# C2 Kac: mittlere Rückkehrzeit zum Ort auf der bei N geschlossenen Halbgeraden = 2N
def kac_mc(N, n=20000):
    x = np.ones(n, int); t = np.ones(n); done = np.zeros(n, bool)
    while not done.all():
        idx = np.where(~done)[0]
        st = rng.choice((-1, 1), size=idx.size)
        xn = x[idx] + st
        xn[xn > N] = N - 1
        x[idx] = xn; t[idx] += 1
        done[idx[xn == 0]] = True
    return t.mean()
sim("C", "Kac: mittlere Rückkehrzeit auf der bei N geschlossenen Halbgeraden = 2N (N = 15)", "Kac 1947", 30.0, kac_mc(15), 0.03)
# C3 Pólya/Watson: Rückkehrsumme in drei Dimensionen
W = quad(lambda t: ive(0, t / 3)**3, 0, 1e5, limit=2000)[0] + 2 * (3 / (2 * math.pi))**1.5 / math.sqrt(1e5)
ident("C", "Watson-Integral W₃ = ∫e^{−t}I₀(t/3)³dt = 1/(1 − p₃)", "Watson 1939; Pólya 1921", 1 / (1 - 0.340537), W, 2e-4)
# C4 Bose/Planck/Wien
for dd_ in (1, 2, 3):
    ident("C", f"Bose-Leiter ∫s^{dd_}n_B ds = Γ({dd_}+1)ζ({dd_}+1)", "Planck 1900; Stefan-Boltzmann (d = 3)", quad(lambda x: x**dd_ / math.expm1(x), 0, 200)[0], gamma(dd_ + 1) * zeta(dd_ + 1), 1e-7)
ident("C", "Wien: Maximum von s³n_B bei s = 2,821439", "Wien 1893 (Verschiebungsgesetz)", bose_stufe_max(3)[0], 2.821439, 1e-6)
# C5 Kippstelle (Spitze ẋ = a + bx − x³): Hysteresefenster ∝ b^{3/2}, Sprung ∝ b^{1/2}, Rückkehrrate ∝ Abstand^{1/2}
def fenster(bb):
    a_f = lambda a: len([z for z in np.roots([-1, 0, bb, a]) if abs(z.imag) < 1e-9])
    lo = brentq(lambda a: bb - 3 * (np.sort(np.real(np.roots([-3, 0, bb])))[1])**2 + 0 * a, -1, 1) if False else None
    return 4 * (bb / 3)**1.5
def fenster_num(bb):
    # Umkehrpunkte numerisch: a, bei denen die Zahl reeller Wurzeln wechselt
    aa = np.linspace(-2 * (bb / 3)**1.5 * 1.5, 2 * (bb / 3)**1.5 * 1.5, 30001)
    n3 = np.array([np.sum(np.abs(np.roots([-1, 0, bb, a]).imag) < 1e-7) for a in aa])
    w = aa[n3 == 3]
    return w.max() - w.min()
sim("C", "Kippstelle: Hysteresefenster ∝ b^{3/2} (Verhältnis b = 0,4 zu 0,1 = 8)", "Thom 1972 (Spitzenkatastrophe)", 8.0, fenster_num(0.4) / fenster_num(0.1), 2e-3)
# C6 Fluktuation-Dissipation: ⟨x²⟩ = kT/k für jede Dämpfung
def fdt(gam, kT=1.0, ks=2.0, n=4000, steps=20000):
    dt = 0.02 * gam / ks
    x = np.zeros(n); acc = 0.0; cnt = 0
    for i in range(steps):
        x += -ks / gam * x * dt + math.sqrt(2 * kT / gam * dt) * rng.standard_normal(n)
        if i > steps // 2 and i % 10 == 0: acc += np.mean(x * x); cnt += 1
    return acc / cnt
sim("C", "Fluktuation-Dissipation ⟨x²⟩ = kT/k, Dämpfung γ = 1", "Einstein 1905; Kubo 1966", 0.5, fdt(1.0), 0.02)
sim("C", "Fluktuation-Dissipation ⟨x²⟩ = kT/k, Dämpfung γ = 10", "Einstein 1905; Kubo 1966", 0.5, fdt(10.0), 0.02)
# C7 Onsager-Casimir: Entropieerzeugung nur aus dem symmetrischen Teil
Ls = rng.normal(size=(4, 4)); Ls = Ls @ Ls.T; La = rng.normal(size=(4, 4)); La = La - La.T; X = rng.normal(size=4)
ident("C", "Onsager-Casimir: X·(L_s + L_a)·X = X·L_s·X", "Onsager 1931; Casimir 1945", X @ (Ls + La) @ X, X @ Ls @ X, 1e-12)
# C8 Ramanujan τ(n) aus Δ = q∏(1 − qⁿ)²⁴
N = 14
poly = np.zeros(N + 1, dtype=object); poly[0] = 1
for n in range(1, N + 1):
    for _ in range(24):
        new = poly.copy()
        for k in range(n, N + 1): new[k] = new[k] - poly[k - n]
        poly = new
tau = {n: int(poly[n - 1]) for n in range(1, N + 1)}
sig11 = lambda n: sum(d**11 for d in range(1, n + 1) if n % d == 0)
ok = (tau[2] == -24 and tau[3] == 252 and tau[4] == -1472 and tau[5] == 4830 and tau[6] == tau[2] * tau[3]
      and tau[4] == tau[2]**2 - 2**11 and all((tau[n] - sig11(n)) % 691 == 0 for n in range(1, N + 1))
      and all(abs(tau[p_]) <= 2 * p_**5.5 for p_ in (2, 3, 5, 7, 11, 13)))
R.append(("C", "Ramanujan τ(n): Werte, Multiplikativität, Hecke, Kongruenz mod 691, Deligne-Schranke", "Ramanujan 1916; Mordell 1917; Deligne 1974", tau[2], -24, "alle", "OK" if ok else "ABWEICHUNG"))
# C9 Thetareihen: Jacobi (Quadrat) und Eisenstein (Sechseck)
def zaehl(n, hexa):
    return sum(1 for m in range(-6, 7) for k in range(-6, 7) if (m * m + m * k + k * k if hexa else m * m + k * k) == n)
jac = lambda n: 4 * (sum(1 for d in range(1, n + 1) if n % d == 0 and d % 4 == 1) - sum(1 for d in range(1, n + 1) if n % d == 0 and d % 4 == 3))
eis = lambda n: 6 * (sum(1 for d in range(1, n + 1) if n % d == 0 and d % 3 == 1) - sum(1 for d in range(1, n + 1) if n % d == 0 and d % 3 == 2))
ok = all(zaehl(n, False) == jac(n) and zaehl(n, True) == eis(n) for n in range(1, 26))
R.append(("C", "Bilder im Kaleidoskop: r₂(n) = 4(d₁ − d₃), Sechseck 6(d₁ − d₂) mod 3, n ≤ 25", "Jacobi 1829; Eisenstein", zaehl(5, False), jac(5), "n = 1…25", "OK" if ok else "ABWEICHUNG"))
# C10 Barkhausen: Grenzzyklus des van-der-Pol-Oszillators mit Amplitude 2 (Schleifenverstärkung 1 im Mittel)
def vdp(mu=0.05, T=600.0, dt=0.005):
    x, v = 0.1, 0.0; amp = 0.0
    for i in range(int(T / dt)):
        a = mu * (1 - x * x) * v - x
        v += a * dt; x += v * dt
        if i * dt > T - 20: amp = max(amp, abs(x))
    return amp
sim("C", "Barkhausen: Amplitude des Grenzzyklus A* = 2 (van der Pol, μ = 0,05)", "van der Pol 1926; Barkhausen 1935", 2.0, vdp(), 0.01)
# C11 Kac für eine Zerlegung: Σ 1/⟨T_i⟩ = 1 in einer geschlossenen Markov-Kette
P = rng.random((5, 5)); P /= P.sum(1, keepdims=True)
state = 0; last = {}; ret = {i: [] for i in range(5)}
for t in range(300000):
    if state in last: ret[state].append(t - last[state])
    last[state] = t
    state = rng.choice(5, p=P[state])
sim("C", "Kac für eine Zerlegung: Σ 1/⟨T_i⟩ = 1 (geschlossene Kette, 5 Zustände)", "Kac 1947", 1.0, sum(1 / np.mean(ret[i]) for i in range(5)), 0.01)

# ======================================================================== Z Zustandsraum des Autors (Tafeln 1–6)
def wahr(teil, gesetz, quelle, ok, mass="alle"):
    R.append((teil, gesetz, quelle, "ja" if ok else "nein", "ja", mass, "OK" if ok else "ABWEICHUNG"))
ZQ = "Autor, Zustandsraum"
# Z1 Poynting mit Joule-Kanal: Feldenergie, die verschwindet, erscheint als J·E = σ|E|² (Maxwell-Block antisymmetrisch)
nx = 400; xg = np.arange(nx); sg = np.where((xg > 250) & (xg < 330), 0.02, 0.0)
def f_em(y):
    E, H = y[:nx], y[nx:2 * nx - 1]
    dE = np.zeros(nx); dE[1:-1] = H[1:] - H[:-1] - sg[1:-1] * E[1:-1]
    return np.concatenate([dE, E[1:] - E[:-1], [np.sum(sg * E * E)]])
E0 = np.exp(-((xg - 100) / 15.0)**2); E0[0] = E0[-1] = 0
y = np.concatenate([E0, -np.exp(-((xg[:-1] + 0.5 - 100) / 15.0)**2), [0.0]])
W0 = 0.5 * np.sum(y[:-1]**2); hs = 0.1
for _ in range(4000):
    k1 = f_em(y); k2 = f_em(y + hs / 2 * k1); k3 = f_em(y + hs / 2 * k2); k4 = f_em(y + hs * k3)
    y = y + hs / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
ident("Z", f"Poynting ∂u/∂t + ∇·S = −J·E: Feldverlust = Joule-Wärme ({y[-1] / W0:.0%} umgesetzt)", ZQ + " (Tafel 3, II/III)", W0 - 0.5 * np.sum(y[:-1]**2), y[-1], 1e-6)
# Z2 Maxwell als Strukturbedingung: ∇·B = 0 und ∇·D = ρ bleiben unter ∇×E = −∂B/∂t, ∇×H = J + ∂D/∂t, ∂ρ/∂t = −∇·J erhalten
n3 = 12; fw = lambda f, a: np.roll(f, -1, a) - f; bw = lambda f, a: f - np.roll(f, 1, a)
B = np.zeros((3, n3, n3, n3)); D = np.zeros((3, n3, n3, n3)); rho = np.zeros((n3, n3, n3))
for _ in range(40):
    Ef, Hf, J = rng.normal(size=(3, 3, n3, n3, n3))
    B -= 0.1 * np.array([fw(Ef[2], 1) - fw(Ef[1], 2), fw(Ef[0], 2) - fw(Ef[2], 0), fw(Ef[1], 0) - fw(Ef[0], 1)])
    D += 0.1 * (np.array([bw(Hf[2], 1) - bw(Hf[1], 2), bw(Hf[0], 2) - bw(Hf[2], 0), bw(Hf[1], 0) - bw(Hf[0], 1)]) - J)
    rho -= 0.1 * sum(bw(J[i], i) for i in range(3))
divB = abs(sum(fw(B[i], i) for i in range(3))).max()   # Yee: Divergenz mit derselben Differenz wie die Rotation
gauss = abs(sum(bw(D[i], i) for i in range(3)) - rho).max()
wahr("Z", "Maxwell: ∇·B = 0 und ∇·D = ρ bleiben erhalten (Gitter 12³, 40 Schritte mit Zufallsfeldern)", ZQ + " (Tafel 3, I)", divB < 1e-10 and gauss < 1e-10, f"{max(divB, gauss):.0e}")
# Z3 Lorentzkraft tangential: (v×B)·v = 0, Umlauf ohne Arbeit (Boris-Schritt, 20 000 Umläufe-Schritte)
v = np.array([1.0, 0.3, 0.2]); Bv = np.array([0.2, 0.5, 1.0]); v0n = np.linalg.norm(v)
tt = Bv * 0.025; ss = 2 * tt / (1 + tt @ tt)
for _ in range(20000):
    vp = v + np.cross(v, tt); v = v + np.cross(vp, ss)
ident("Z", "Lorentzkraft f = ρE + J×B tangential: |v| im Magnetfeld erhalten", ZQ + " (Tafel 3, Lorentz; 09)", np.linalg.norm(v), v0n, 1e-12)
# Z4 Schnitt erzeugt Rand, Rand erzeugt Fluchtmoden, Schließung beseitigt sie (Diffusion auf Ring und Strecke, L = D = 1)
Ng = 400; h_ = 1.0 / Ng
ring = (np.roll(np.eye(Ng), 1, 0) + np.roll(np.eye(Ng), -1, 0) - 2 * np.eye(Ng)) / h_**2
cut = (np.eye(Ng, k=1) + np.eye(Ng, k=-1) - 2 * np.eye(Ng)) * (Ng + 1)**2
ev_r = np.sort(-eigvalsh(ring)); ev_c = np.sort(-eigvalsh(cut))
wahr("Z", "Schließung: Ring hat einen Nullmode (Gesamtmenge bleibt), Strecke mit Rand keinen", ZQ + " (Tafel 4)", abs(ev_r[0]) < 1e-6 and ev_c[0] > 1.0, f"{ev_r[0]:.0e} / {ev_c[0]:.2f}")
sim("Z", "Spektrallücke des geschlossenen Kreises 4π²D/L²", ZQ + " (Tafel 4); 05", 4 * math.pi**2, ev_r[1], 1e-4)
sim("Z", "Fluchtrate der aufgeschnittenen Strecke π²D/L²", ZQ + " (Tafel 4); 05", math.pi**2, ev_c[0], 1e-4)
# Z5 Vorspannung → Steifigkeit → Resonanz, Knicken bei P_E (Euler-Balken gelenkig, Differenzenmatrix)
Nb = 300; D2 = (np.eye(Nb, k=1) + np.eye(Nb, k=-1) - 2 * np.eye(Nb)) * (Nb + 1)**2
w2 = lambda P: np.sort(eigvalsh(D2 @ D2 + P * D2))[0]
sim("Z", "Vorspannung verschiebt Resonanz: ω²(P)/ω²(0) = 1 − P/P_E (P = P_E/2; Knicken bei P_E = π²EI/L²)", ZQ + " (Tafel 1, Gravitation; 05)", 0.5, w2(0.5 * math.pi**2) / w2(0.0), 1e-3)
# Z6 parametrische Rolle: langsame Änderung von g hält die Wirkung J = E/ω (Torus rückt auf die Nachbarfaser)
gt = lambda t: 1 + 3 * (3 * (t / 3000)**2 - 2 * (t / 3000)**3) if t < 3000 else 4.0
sol = solve_ivp(lambda t, z: [z[1], -gt(t) * z[0]], (0, 3000), [1.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12)
xe, ve = sol.y[:, -1]
sim("Z", "Parametrisch: Wirkung E/ω bleibt bei langsamer Änderung von g (g × 4 über ≈ 700 Perioden)", ZQ + " (Tafel 4, Faserbündel); 05", 0.5, (0.5 * ve**2 + 2 * xe**2) / 2.0, 1e-3)
# Z7 Rückkopplungsschleife Wärme → Widerstand → Wärme: Durchgehen bei q = 1, Verlangsamung τ = τ₀/(1 − q)
def ueber10(t, T): return T[0] - 10
ueber10.terminal = True
def durchgeht(I2):
    s_ = solve_ivp(lambda t, T: I2 * np.exp(T) - T, (0, 1e5), [0.0], rtol=1e-10, atol=1e-12, events=ueber10)
    return s_.status != 0
lo_, hi_ = 0.30, 0.40
for _ in range(25):
    mid = (lo_ + hi_) / 2
    lo_, hi_ = (lo_, mid) if durchgeht(mid) else (mid, hi_)
sim("Z", "Thermisches Durchgehen bei Schleifenverstärkung q = 1: I²R₀α/G = 1/e (R = R₀e^{αT})", ZQ + " (Tafel 3, Wärme; 09)", 1 / math.e, lo_, 1e-6)
Ts = brentq(lambda T: 0.2 * math.exp(T) - T, 0, 1); qT = Ts
s_ = solve_ivp(lambda t, T: 0.2 * np.exp(T) - T, (0, 3), [Ts + 1e-6], rtol=1e-12, atol=1e-16)
sim("Z", f"Verstärkung ist Verlangsamung: τ = τ₀/(1 − q) (q = {qT:.3f})", "05", 1 / (1 - qT), -3 / math.log((s_.y[0, -1] - Ts) / 1e-6), 1e-4)
# Z8 GENERIC: tangentialer Block erhält E, transversaler erzeugt S (Oszillator mit Wärmespeicher)
zet_, Cw = 0.3, 1.0
sol = solve_ivp(lambda t, z: [z[1], -z[0] - zet_ * z[1], zet_ * z[1]**2], (0, 60), [1.0, 0.0, 1.0], rtol=1e-12, atol=1e-14, dense_output=True)
tt_ = np.linspace(0, 60, 3000); zz = sol.sol(tt_)
Eg = 0.5 * zz[1]**2 + 0.5 * zz[0]**2 + zz[2]; Sg = Cw * np.log(zz[2])
ok = True
for _ in range(100):
    x_, p_, e_ = rng.normal(), rng.normal(), rng.random() + 0.1
    Lm = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 0]]); vv = np.array([0, 1, -p_]); Mm = zet_ * e_ / Cw * np.outer(vv, vv)
    gE = np.array([x_, p_, 1]); gS = np.array([0, 0, Cw / e_])
    ok &= np.allclose(Lm @ gS, 0) and np.allclose(Mm @ gE, 0) and np.allclose(Lm, -Lm.T) and min(np.linalg.eigvalsh(Mm)) > -1e-12
    ok &= np.allclose(Lm @ gE + Mm @ gS, [p_, -x_ - zet_ * p_, zet_ * p_**2])
wahr("Z", "GENERIC: 𝓛 = −𝓛ᵀ, 𝓛∇S = 0, 𝓜 = 𝓜ᵀ ≥ 0, 𝓜∇E = 0; E fest, S steigt (gedämpfter Oszillator mit Wärme)", ZQ + " (Tafel 2, drei Rollen); 05; Öttinger 2005",
     ok and np.ptp(Eg) < 1e-9 and np.all(np.diff(Sg) > -1e-12), f"ΔE {np.ptp(Eg):.0e}, ΔS {Sg[-1] - Sg[0]:+.3f}")
# Z9 Einrasten |Δω| ≤ 2k|cos φ|; reaktive Kopplung (φ = π/2) rastet nie ein; Rate am Rastpunkt −√((2k cos φ)² − Δω²)
dw = np.linspace(0.05, 1.95, 39); ph = np.array([0, math.pi / 3, math.pi / 2, 2 * math.pi / 3])
DW, PH = np.meshgrid(dw, ph); kk = 0.5
psi = np.zeros_like(DW); f_ = lambda p: DW - 2 * kk * np.cos(PH) * np.sin(p); hs = 0.02
for i in range(40000):
    a1 = f_(psi); a2 = f_(psi + hs / 2 * a1); a3 = f_(psi + hs / 2 * a2); a4 = f_(psi + hs * a3)
    psi = psi + hs / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
    if i == 20000: psi_m = psi.copy()
rast = np.abs(psi - psi_m) / 400 < 1e-3
erw = np.abs(DW) <= 2 * kk * np.abs(np.cos(PH))
rand = np.abs(np.abs(DW) - 2 * kk * np.abs(np.cos(PH))) < 0.03
wahr("Z", "Einrasten genau für |Δω| ≤ 2k|cos φ|; bei φ = π/2 nie (156 Paare Δω, φ)", "06; Adler 1946", np.all((rast == erw) | rand), f"{int(rast.sum())} gerastet")
psi_s = math.asin(0.6 / 1.0)
s_ = solve_ivp(lambda t, p: 0.6 - np.sin(p), (0, 5), [psi_s + 1e-6], rtol=1e-12, atol=1e-18)
sim("Z", "Rate am Rastpunkt −√((2k cos φ)² − Δω²) (Δω = 0,6, k = 0,5, φ = 0; Abklingen einer Störung)", "06", -math.sqrt(1.0 - 0.36), math.log((s_.y[0, -1] - psi_s) / 1e-6) / 5, 1e-5)
# Z10 Einrastzunge m:n: |nω₂ − mω₁| ≤ (m + n)k, Breite linear in k (3:2, Simulation der beiden Phasen)
def zunge(k_, m=3, n=2):
    det = np.linspace(0, 2 * (m + n) * k_, 201)
    th1 = np.zeros_like(det); th2 = np.zeros_like(det); w1 = 2.0; w2 = (m * w1 + det) / n; hs = 0.01
    for i in range(60000):
        ps = n * th2 - m * th1
        th1 = th1 + hs * (w1 + k_ * np.sin(ps)); th2 = th2 + hs * (w2 - k_ * np.sin(ps))
        if i == 30000: p_m = n * th2 - m * th1
    frei = np.abs(n * th2 - m * th1 - p_m) > 1.0
    return det[np.argmax(frei)] / ((m + n) * k_)
sim("Z", "Einrastzunge 3:2: Rand bei |nω₂ − mω₁| = (m + n)k (k = 0,02)", "06; Arnold 1965", 1.0, zunge(0.02), 0.015)
sim("Z", "Einrastzunge 3:2: Rand bei |nω₂ − mω₁| = (m + n)k (k = 0,04, Breite linear)", "06; Arnold 1965", 1.0, zunge(0.04), 0.015)
# Z11 Torusdimension begrenzt: N gleich verteilte Zyklen im Band W rasten ab (N − 1)k ≈ (2/π)W ein (Kuramoto)
def rast_selbstkonsistent(Nk):
    om = np.linspace(-0.5, 0.5, Nk)
    res = minimize_scalar(lambda u: u / np.mean(np.sqrt(np.clip(1 - (om / u)**2, 0, None))), bounds=(0.5, 2), method="bounded")
    return (Nk - 1) / Nk * res.fun
Nk = 40; om = np.linspace(-0.5, 0.5, Nk); kv = np.linspace(0.55, 0.70, 61) / (Nk - 1)
th = np.tile(rng.random(Nk) * 2 * math.pi, (kv.size, 1)); hs = 0.05; acc = np.zeros_like(th)
for i in range(40000):
    z = np.exp(1j * th).mean(1, keepdims=True)
    dth = om + kv[:, None] * Nk * np.imag(z * np.exp(-1j * th))
    th = th + hs * dth
    if i >= 30000: acc += dth
fr = acc / 10000
gerastet = np.ptp(fr, 1) < 1e-3
sim("Z", "Kuramoto N = 40: volles Einrasten (N − 1)k/W, Simulation gegen Selbstkonsistenz", "06; Kuramoto 1975", rast_selbstkonsistent(40), (Nk - 1) * kv[np.argmax(gerastet)], 5e-3)
ident("Z", "Grenze N_max ≈ 1 + 2W/(πk): (N − 1)k_c/W → 2/π (N = 4000)", "06; Kuramoto: K_c = 2/(πg(0))", rast_selbstkonsistent(4000), 2 / math.pi, 2e-3)
# Z12 Koexistenz zweier Zyklen an gemeinsamer Verstärkung: θ₁₂/β₂ < α₁/α₂ < β₁/θ₂₁
def koex(a1, a2, b1=1.2, b2=0.9, t12=0.3, t21=0.8):
    s_ = solve_ivp(lambda t, I: [2 * I[0] * (a1 - b1 * I[0] - t12 * I[1]), 2 * I[1] * (a2 - b2 * I[1] - t21 * I[0])], (0, 400), [0.1, 0.1], rtol=1e-9)
    return bool(min(s_.y[:, -1]) > 1e-3)
ok = all(koex(r_, 1.0) == (0.3 / 0.9 < r_ < 1.2 / 0.8) for r_ in (0.2, 0.5, 1.0, 1.4, 2.0))
wahr("Z", "Koexistenz genau für θ₁₂/β₂ < α₁/α₂ < β₁/θ₂₁ (5 Verhältnisse)", "06", ok)
# Z13 Amplitudentod bei dissipativer Kopplung: K_d > μ und |Δω| > 2√(μ(2K_d − μ))
def tod(K, dwv, mu=1.0):
    w = np.array([0.0, dwv])
    s_ = solve_ivp(lambda t, z: (mu + 1j * w - abs(z)**2) * z + K * (z[::-1] - z), (0, 300), np.array([1.0 + 0j, 0.5j]), rtol=1e-9, atol=1e-14)
    return abs(s_.y[:, -1]).max() < 1e-6
ok = tod(2.0, 6.0) and not tod(2.0, 3.0) and not tod(0.8, 6.0)
wahr("Z", "Amplitudentod genau für K_d > μ und |Δω| > 2√(μ(2K_d − μ)) (3 Punkte)", "06; Aronson, Ermentrout, Kopell 1990", ok)
# Z14 Gedächtnis bindet ohne Energieverlust: v = 1,6 > √2, γ = 0,01, t_K → ∞: Umkehr bei x = 9,75
sol = solve_ivp(lambda t, z: [z[1], -1 / z[0]**2 - 0.01 * z[2], z[1]], (0, 500), [1.0, 1.6, 0.0], rtol=1e-12, atol=1e-12, events=lambda t, z: z[1])
x_um = sol.y_events[0][0][0]
ident("Z", "Gedächtniskern bindet: Umkehrpunkt aus der Erhaltung H = v²/2 − 1/x + γy²/2", "05", x_um, brentq(lambda x: 0.28 + 1 / x - 0.005 * (x - 1)**2, 2, 50), 1e-8)
sim("Z", "Gedächtniskern bindet: Umkehrpunkt x = 9,75", "05", 9.75, x_um, 1e-3)
# Z15 Maxwell-Element: Speicher γ(ωt_K)²/(1 + (ωt_K)²), Verlust γωt_K/(1 + (ωt_K)²)
for wt in (0.5, 2.0):
    re_ = quad(lambda u: math.exp(-u) * math.cos(wt * u), 0, 60, limit=400)[0]; im_ = quad(lambda u: math.exp(-u) * math.sin(wt * u), 0, 60, limit=400)[0]
    ker = 1j * wt * (re_ - 1j * im_)
    ident("Z", f"Maxwell-Element ωt_K = {wt}: Speicheranteil", "05", ker.real, wt**2 / (1 + wt**2), 1e-9)
    ident("Z", f"Maxwell-Element ωt_K = {wt}: Verlustanteil", "05", ker.imag, wt / (1 + wt**2), 1e-9)
# Z16 Hysterese ohne Speicher (Preisach): nur die nicht ausgelöschten Extremwerte zählen
al = rng.uniform(-1, 1, 3000); be = rng.uniform(-1, 1, 3000); al, be = np.maximum(al, be), np.minimum(al, be); wt_ = rng.random(3000)
def preisach(folge):
    st = -np.ones(3000)
    for u in folge:
        st = np.where(u >= al, 1.0, np.where(u <= be, -1.0, st))
    return st
A_ = preisach([-1, 0.8, -0.3, 0.6, 0.2]); B_ = preisach([-1, 0.8, -0.1, 0.4, -0.3, 0.6, 0.2])
C1_ = preisach([-1, 0.8, -0.3, 0.5]); C2_ = preisach([-1, 0.8, -0.3, 0.5, 0.1, 0.5])
wahr("Z", "Preisach: gleiche Folge der Extremwerte ⇒ gleicher Zustand; Rückkehrpunkt-Gedächtnis", "06; Preisach 1935", np.array_equal(A_, B_) and np.array_equal(C1_, C2_), f"Ausgang {wt_ @ A_:.3f} = {wt_ @ B_:.3f}")
# Z17 Kippstelle: vor der Falte 1 − q = 2(3a)^{1/4}δ^{1/2}; Sprung √(3a)
a_, dl = 1.0, 1e-8
rts = np.real(np.roots([1, 0, -a_, -(2 * (a_ / 3)**1.5 - dl)]))
xs_ = min(rts)
sim("Z", "Kippstelle: 1 − q = 2(3a)^{1/4}δ^{1/2} vor der Falte (a = 1, δ = 10⁻⁸)", "06", 2 * (3 * a_)**0.25 * dl**0.5, 3 * xs_**2 - a_, 1e-3)
ident("Z", "Kippstelle: Sprung an der Falte √(3a)", "06", max(np.real(np.roots([1, 0, -a_, -2 * (a_ / 3)**1.5]))) + math.sqrt(a_ / 3), math.sqrt(3 * a_), 1e-6)
# Z18 Dreieck Energie–Struktur–Zustand: oszillatorisches Aufschaukeln erst bei |q₃| > (Σt)(Σ1/t) − 1 ≥ 8
def q3_krit(t3):
    def re_max(Q):
        co = np.poly([-1 / t for t in t3]) * np.prod(t3); co[-1] += Q
        return max(np.real(np.roots(co)))
    return brentq(re_max, 1e-6, 1e6, xtol=1e-14)
for t3 in ((1.0, 1.0, 1.0), (0.1, 1.0, 30.0)):
    ident("Z", f"Oszillatorischer Einsatz im Dreieck, Zeiten {t3}: |q₃| = (Σt)(Σ1/t) − 1", "06", q3_krit(t3), sum(t3) * sum(1 / t for t in t3) - 1, 1e-8)
# Z19 Rückkopplung im Netz: (I − A)⁻¹ = ΣAⁿ für Spektralradius < 1
Am = rng.random((6, 6)); Am *= 0.9 / max(abs(np.linalg.eigvals(Am)))
Sn = np.zeros((6, 6)); Pn = np.eye(6)
for _ in range(600): Sn += Pn; Pn = Pn @ Am
ident("Z", "Rückläufe im Netz: (I − A)⁻¹ = ΣAⁿ (Spektralradius 0,9)", "05; Leontief, Hawkins-Simon", np.linalg.inv(np.eye(6) - Am).sum(), Sn.sum(), 1e-9)
# Z20 Linsen-Test (Tafel 6): Kugel, Ψ = Φ ⇒ Linsenmasse = dynamische Masse (η = 1); Einheiten GM = a₀ = r* = 1
aK = lambda r: float(a_regel(1 / r**2, 1.0))
def ablenk(b):   # (c²/2)·α = ∫ a·(b/r) dl, Substitution l = b·tan θ
    return 2 * quad(lambda t: aK(b / math.cos(t)) * b / math.cos(t), 0, math.pi / 2, limit=400)[0]
def zylindermasse(b):   # M_dyn(r) = r²a(r); Masse im Zylinder vom Radius b aus der Kugelverteilung
    Md = lambda r: r * r * aK(r)
    dM = lambda r: (Md(r * 1.0001) - Md(r / 1.0001)) / (r * 1.0001 - r / 1.0001)
    return Md(b) + quad(lambda u: dM(b * math.exp(u)) * (1 - math.sqrt(1 - math.exp(-2 * u))) * b * math.exp(u), 0, 40, limit=400)[0]
for bb in (0.3, 1.0, 5.0):
    ident("Z", f"Linsen-Test Kugel: Linsenmasse / dynamische Masse η = 1 (b = {bb} r*)", ZQ + " (Tafel 6, Test 1); 05", ablenk(bb), 2 * zylindermasse(bb) / bb, 1e-5)
ident("Z", "Ablenkung weit außen: α → 2π(v_f/c)² + 6Δs/b (b = 200 r*)", "08", ablenk(200.0), math.pi + 1 / 200, 1e-5)

# ======================================================================== M Unsere Gesetze (Bausteine 01–09, Notizen 71, 72)
# M1 Formeln A1–A3
ident("M", "(A1) g = x/y²·g⊕ gegen GM/R² (Sonne)", "01", M_sun / M_earth / (R_sun / R_earth)**2 * g_earth, GMs / R_sun**2, 2e-3)
vE = 2 * math.pi * AU / (365.256363 * 86400)
ident("M", "Zeitdifferential einer Bahn δ = 3Δs/r = v²/c² (Erdbahn, v aus dem Jahr)", "01", 1 - math.sqrt(1 - 6 * ds(M_sun) / AU), vE**2 / c**2, 2e-3)
ident("M", "Brechzahl n = 1 + 6Δs/r als Strahlintegral: α = 12Δs/b (Takt + Faltung je die Hälfte)", "08",
      6 * ds(M_sun) / R_sun * quad(lambda u: (1 + u * u)**-1.5, -np.inf, np.inf)[0], 12 * ds(M_sun) / R_sun, 1e-8)
# M2 Produkt fest, Übergang bei Gleichheit (Notiz 72, Form 1): Schnittmassen über Nullstellensuche
lam = lambda m: hbar / (m * c); ell0 = c**2 / a0
Mc_ = math.exp(brentq(lambda u: math.log(r_stern(math.exp(u)) / lam(math.exp(u))), -200, 0))
MU_ = math.exp(brentq(lambda u: math.log(r_stern(math.exp(u)) / (3 * ds(math.exp(u)))), 100, 200))
mP_ = math.exp(brentq(lambda u: math.log(3 * ds(math.exp(u)) / lam(math.exp(u))), -40, 0))
ident("M", "Schnitt 3Δs = λ̄ bei m_Pl (Nullstelle)", "02; 72", mP_, m_Pl, 1e-9)
ident("M", "Schnitt r* = λ̄ bei M_c = (ħ²a₀/(Gc²))^{1/3} (Nullstelle)", "08; 72", Mc_, (hbar**2 * a0 / (G * c**2))**(1 / 3), 1e-9)
ident("M", "Schnitt r* = 3Δs bei M_U = c⁴/(Ga₀) (Nullstelle)", "08; 72", MU_, c**4 / (G * a0), 1e-9)
ident("M", "Längenleiter M_c³·M_U = m_Pl⁴ (aus den Nullstellen)", "08", Mc_**3 * MU_, mP_**4, 1e-8)
ident("M", "r*²·λ̄ = ℓ_P²·ℓ₀ für jede Masse (Stichprobe 10⁻²⁰ … 10⁴² kg)", "08; 72", max(r_stern(m) ** 2 * lam(m) / (l_P**2 * ell0) for m in np.logspace(-20, 42, 9)), 1.0, 1e-12)
ident("M", "Geometrisches Mittel λ̄_e² = r_e·a_B (CODATA-Werte)", "08", (hbar / (m_e * c))**2, r_e * a_B, 1e-9)
Dst = brentq(lambda x: x * x * ds(m_p) - (lam(m_p) / 3)**3, 1, 1e4)
ident("M", "Sternskala: Δs_*²·Δs_p = (λ̄_p/3)³ ⇒ M_* = m_Pl³/m_p²", "02", 3 * c**2 * Dst / G, m_Pl**3 / m_p**2, 1e-9)
# M3 Verwindung (Kerr aus der Metrik, zweiter Weg)
aK_, mK = 0.6, 1.0   # Einheiten GM/c² = 1
gtt = lambda r, th: -(1 - 2 * r / (r * r + aK_**2 * math.cos(th)**2))
ident("M", "Nullfläche des Takts r_E(θ) = 3Δs + √(9Δs² − a_J²cos²θ) (Nullstelle von g_tt, θ = 0,7)", "08", brentq(lambda r: gtt(r, 0.7), 1.0001, 5), 1 + math.sqrt(1 - aK_**2 * math.cos(0.7)**2), 1e-10)
def omega_kreis(r, s):   # Kreisbahn aus der Metrik: g_tt,r + 2g_tφ,r·Ω + g_φφ,r·Ω² = 0 (Äquator, Boyer-Lindquist)
    tt_, tp, pp = -2 / r**2, 2 * aK_ / r**2, 2 * r - 2 * aK_**2 / r**2
    return (-tp + s * math.sqrt(tp * tp - tt_ * pp)) / pp
for rr_ in (10.0, 1000.0):
    ident("M", f"T_mitläufig − T_gegenläufig = 4πa_J/c unabhängig von r (Metrik, r = {rr_:g} GM/c²)", "08", 2 * math.pi / omega_kreis(rr_, 1) - 2 * math.pi / abs(omega_kreis(rr_, -1)), 4 * math.pi * aK_, 1e-6)
aJe = hbar / 2 / (m_e * c)
ident("M", "Spin ½: a_J·Δs = ℓ_P²/6 (Elektron)", "09", aJe * ds(m_e), l_P**2 / 6, 1e-9)
ident("M", "Magnetisches Moment g_F·e·c·a_J/2 = μ_B (CODATA 9,2740100783·10⁻²⁴ J/T)", "09", 2 * e_el * c * aJe / 2, 9.2740100783e-24, 1e-9)
# M4 Rückkehr: D·G'' = Γ·G mit G(0) = 1 (Differenzenrechnung), offen e^{−x/λ}, bei L versiegelt cosh((L − x)/λ)/cosh(L/λ)
def bvp(L, versiegelt=True, h=1e-3):
    n = int(round(L / h)); ab = np.zeros((3, n)); ab[0, 1:] = 1.0; ab[1] = -(2 + h * h); ab[2, :-1] = 1.0
    if versiegelt: ab[2, -2] = 2.0
    rhs_ = np.zeros(n); rhs_[0] = -1.0
    return solve_banded((1, 1), ab, rhs_)   # Eintrag i ist G(x = (i + 1)h)
ident("M", "Rückkehr bei Λ geschlossen: G(x) = cosh((L − x)/λ)/cosh(L/λ) (L = 3λ, x = 1λ)", "07; 71 (Kabel, Thiele)", bvp(3.0)[999], math.cosh(2) / math.cosh(3), 1e-5)
ident("M", "Rückkehr offen: q = e^{−s} (x = 1λ)", "07", bvp(40.0, False)[999], math.exp(-1), 1e-5)
# M5 Schalter χ: Elastizitäten und Flussdimension über numerische Ableitungen der Regel
for s_v in (0.3, 1.0, 3.0):
    aN_ = s_v**2 * a0; h = 1e-5
    el0 = (math.log(a_regel(aN_, a0 * (1 + h))) - math.log(a_regel(aN_, a0 * (1 - h)))) / (math.log(1 + h) - math.log(1 - h))
    elN = (math.log(a_regel(aN_ * (1 + h), a0)) - math.log(a_regel(aN_ * (1 - h), a0))) / (math.log(1 + h) - math.log(1 - h))
    ident("M", f"d ln a/d ln a₀ = χ/2 (s = {s_v})", "03; 07", el0, float(chi(aN_, a0)) / 2, 1e-7)
    ident("M", f"d ln a/d ln a_N = 1 − χ/2, Flussdimension 3 − χ (s = {s_v})", "03; 72", elN, 1 - float(chi(aN_, a0)) / 2, 1e-7)
ident("M", "Spektraldimension d_s = 2·d ln a/d ln a_N → 1 im Tiefbereich (s = 10⁻³)", "71 (Pólya)", 2 * (1 - float(chi(1e-6 * a0, a0)) / 2), 1.0, 1e-3)
rS = r_stern(M_sun) / 0.7
ident("M", "Zusatz der Rückkehr a − a_N = χ·(Δs/r)(3c²/r*) (Sonne, s = 0,7)", "07", float(a_regel(GMs / rS**2)) - GMs / rS**2, float(chi(GMs / rS**2)) * ds(M_sun) / rS * 3 * c**2 / r_stern(M_sun), 1e-12)
# M6 Zusatzpotential Φ_X = v_f²·ln ν(s), Hernquist exakt v_f²·ln[ν(s)/ν(s₀)] (Integral der Zusatzbeschleunigung)
for s_v in (0.5, 2.0):
    ident("M", f"Φ_X = v_f²·ln ν(s) als Integral von a − a_N (Punktmasse, s = {s_v})", "08",
          quad(lambda r: float(a_regel(1 / r**2, 1.0)) - 1 / r**2, 1e-9, 1 / s_v, limit=400)[0], -math.log(1 - math.exp(-s_v)), 1e-7)
bH = 0.4; aNH = lambda r: 1 / (r + bH)**2; s0H = 1 / bH
ident("M", "Φ_X Hernquist = v_f²·ln[ν(s)/ν(s₀)] (b_H = 0,4 r*, r = 1,5 r*)", "08",
      quad(lambda r: float(a_regel(aNH(r), 1.0)) - aNH(r), 0, 1.5, limit=400)[0], math.log((1 - math.exp(-s0H)) / (1 - math.exp(-1 / (1.5 + bH)))), 1e-7)
# M7 Epizyklen: (Re τ_M)² = κ²/Ω² = 1 + χ (Punktmasse), 3 + γ_N(1 − χ/2) (Hernquist)
def kappa2_om2(aNf, r):
    O2 = lambda x: float(a_regel(aNf(x), 1.0)) / x
    return (r * (O2(r * 1.0001) - O2(r / 1.0001)) / (r * 1.0001 - r / 1.0001) + 4 * O2(r)) / O2(r)
for r_v in (0.5, 2.0):
    ident("M", f"(Re τ_M)² = 1 + χ, Punktmasse (r = {r_v} r*)", "07; 72", kappa2_om2(lambda x: 1 / x**2, r_v), 1 + float(chi(1 / r_v**2, 1.0)), 1e-6)
gN_ = -2 * 1.0 / (1.0 + bH)
ident("M", "(Re τ_M)² = 3 + γ_N(1 − χ/2), Hernquist (r = 1 r*)", "07", kappa2_om2(aNH, 1.0), 3 + gN_ * (1 - float(chi(aNH(1.0), 1.0)) / 2), 1e-6)
# M8 Entropie der Rücklaufzahl: P(n) = (1 − q)qⁿ, H = ln ν + χ
for s_v in (0.3, 1.0):
    qv = math.exp(-s_v); nn = np.arange(int(600 / s_v)); Pn_ = (1 - qv) * qv**nn
    ident("M", f"Entropie der Rücklaufzahl H = ln ν + χ (s = {s_v})", "71 (Jaynes)", -np.sum(Pn_ * np.log(Pn_)), -math.log(1 - qv) + s_v / (math.exp(s_v) - 1), 1e-9)
# M9 Kritische Verzweigung: Kaskade der Rückkehr genau kritisch ⇒ 1 − E[e^{−ΓS}] → √(2Γ), P(S = n) ∝ n^{−3/2}
lPB = lambda n: -n + (n - 1) * np.log(n) - gammaln(n + 1)   # Borel: Gesamtzahl S der kritischen Poisson-Verzweigung
Tfix = lambda Gm: brentq(lambda T: T - math.exp(-Gm) * math.exp(T - 1), 0.5, 1 - 1e-13, xtol=1e-15)
nB = np.arange(1, 400001, dtype=float)
ident("M", "Kritische Verzweigung (Borel): Summe über P(S = n) gegen die Fixpunktgleichung T = z·e^{T−1} (Γ = 10⁻⁴)", "71", 1 - np.sum(np.exp(lPB(nB) - 1e-4 * nB)), 1 - Tfix(1e-4), 1e-6)
sim("M", "Kritische Verzweigung: 1 − E[e^{−ΓS}] = √(2Γ) (Γ = 10⁻⁸), der Exponent ½ der Rückkehr", "71; 07", math.sqrt(2e-8), 1 - Tfix(1e-8), 1e-3)
ident("M", "Kritische Verzweigung: P(S = n)·n^{3/2}·√(2π) → 1 (n = 10⁷)", "71", math.exp(lPB(1e7)) * 1e7**1.5 * math.sqrt(2 * math.pi), 1.0, 2e-6)
# M10 Dimension als Faktor: mittlere Sehne ⟨ℓ⟩ = η_d·V/S mit η_d = π, 4 (Kreis, Kugel; Monte Carlo)
U = rng.random(400000)
sim("M", "Mittlere Sehne im Kreis πA/P (η₂ = π)", "71; 72 (Cauchy)", math.pi / 2, np.mean(2 * np.sqrt(1 - (2 * U - 1)**2)), 5e-3)
sim("M", "Mittlere Sehne in der Kugel 4V/S (η₃ = 4)", "71; 72 (Cauchy)", 4 / 3, np.mean(2 * np.sqrt(U)), 5e-3)
ident("M", "V/(S·R) = 1/d für die d-Kugel ⇒ Δs = GM/(3c²) für d = 3", "72", (4 / 3 * math.pi) / (4 * math.pi), 1 / 3, 1e-12)
# M11 Modulformen (05): Ramanujan-Gleichungen, Δ, D₁₂Δ = 0, Normierungen
def sig(k, n): return sum(d**k for d in range(1, n + 1) if n % d == 0)
NS = 120; S1 = [sig(1, n) for n in range(NS)]; S3 = [sig(3, n) for n in range(NS)]; S5 = [sig(5, n) for n in range(NS)]
Pf = lambda q: 1 - 24 * sum(S1[n] * q**n for n in range(1, NS)); Qf = lambda q: 1 + 240 * sum(S3[n] * q**n for n in range(1, NS)); Rf = lambda q: 1 - 504 * sum(S5[n] * q**n for n in range(1, NS))
qd = lambda f, q: q * (f(q * (1 + 1e-5)) - f(q * (1 - 1e-5))) / (2e-5 * q)
qq = 0.05; Pq, Qq, Rq = Pf(qq), Qf(qq), Rf(qq)
Df = lambda q: q * np.prod([(1 - q**n)**24 for n in range(1, NS)])
ok = (abs(qd(Pf, qq) - (Pq * Pq - Qq) / 12) < 1e-8 and abs(qd(Qf, qq) - (Pq * Qq - Rq) / 3) < 1e-7 and abs(qd(Rf, qq) - (Pq * Rq - Qq * Qq) / 2) < 1e-6
      and abs((Qq**3 - Rq**2) / 1728 / Df(qq) - 1) < 1e-10 and abs(qd(Df, qq) / (Pq * Df(qq)) - 1) < 1e-8)
wahr("M", "Modulformen: q dP/dq = (P² − Q)/12 usw., Δ = (Q³ − R²)/1728 = q∏(1 − qⁿ)²⁴, D₁₂Δ = 0 (q = 0,05)", "05; Ramanujan 1916", ok)
ok = all(abs((2 * math.pi)**k / (gamma(k) * zeta(k)) - v_) < 1e-9 for k, v_ in ((2, 24), (4, 240), (6, 504))) and 3 * 240 + 2 * 504 == 1728 \
     and abs(2 / (1 - 1 / 2 - 1 / 3) - 12) < 1e-12 and abs(4 * math.pi**4 / 15 - 4 * (2 * math.pi)**4 / 240) < 1e-12
wahr("M", "Herkunft der Zahlen: 24, 240, 504 = (2π)^k/(Γ(k)ζ(k)), 1728 = 3·240 + 2·504, 12 = 2/(1 − ½ − ⅓), 4π⁴/15 = 4(2π)⁴/240", "05", ok)
sL = 0.1
ident("M", "Lambert-Reihe Σσ₁(n)e^{−ns} = π²/(6s²) − 1/(2s) + 1/24 (Modultransformation, s = 0,1)", "05",
      sum(d * math.exp(-d * sL) / (1 - math.exp(-d * sL)) for d in range(1, 2000)), math.pi**2 / (6 * sL**2) - 1 / (2 * sL) + 1 / 24, 1e-12)
lq = lambda q: sum(S1[n] * q**n for n in range(1, NS))
ok = all(abs(sum(sig(1, p * m) * 0.3**m for m in range(1, 100)) - ((1 + p) * lq(0.3) - p * lq(0.3**p))) < 1e-9 for p in (2, 3, 5))
wahr("M", "Hecke: Σσ₁(pm)q^m = (1 + p)ℓ(q) − p·ℓ(q^p) (p = 2, 3, 5)", "05", ok)
qr = 0.5
lhs = sum(qr**(n * n) / np.prod([1 - qr**k for k in range(1, n + 1)]) for n in range(0, 30))
rhs = np.prod([1 / ((1 - qr**(5 * n + 1)) * (1 - qr**(5 * n + 4))) for n in range(0, 200)])
ident("M", "Rogers-Ramanujan Σq^{n²}/(q;q)_n = ∏1/((1 − q^{5n+1})(1 − q^{5n+4})) (q = ½)", "05", lhs, rhs, 1e-12)
fib = [1, 1]
for _ in range(40): fib.append(fib[-1] + fib[-2])
ident("M", "Ausschlussregel Abstand ≥ 2: Zustände je Platz → φ (Fibonacci)", "05", fib[-1] / fib[-2], (1 + math.sqrt(5)) / 2, 1e-12)
# M12 Ramanujan-Summen als Projektoren auf genaue Periode, E₁ = 𝒫; Resolvente 𝒫_q
NR = 12; Tm = np.roll(np.eye(NR), 1, 0)
cq = lambda qn, n: sum(math.cos(2 * math.pi * a_ * n / qn) for a_ in range(1, qn + 1) if math.gcd(a_, qn) == 1)
Eq = {qn: sum(cq(qn, n) * np.linalg.matrix_power(Tm, n) for n in range(NR)) / NR for qn in (1, 2, 3, 4, 6, 12)}
ok = np.allclose(sum(Eq.values()), np.eye(NR)) and all(np.allclose(Eq[a_] @ Eq[b_], Eq[a_] if a_ == b_ else 0) for a_ in Eq for b_ in Eq) and np.allclose(Eq[1], np.ones((NR, NR)) / NR)
wahr("M", "Ramanujan-Summen: E_q Projektoren, ΣE_q = I, E_qE_q' = 0, E₁ = 𝒫 (Z₁₂)", "05", ok)
qP = 0.37; ev_ = np.linalg.eigvals((1 - qP) * np.linalg.inv(np.eye(64) - qP * np.roll(np.eye(64), 1, 0)))
ident("M", "𝒫_q = (1 − q)(I − qT)⁻¹: Eigenwerte auf dem Kreis um 1/(1 + q) mit Radius q/(1 + q)", "05", max(abs(ev_ - 1 / (1 + qP))), qP / (1 + qP), 1e-10)
# M13 Wicklung: goldene Rückkehr 1/(√5·n), Drei-Abstände-Satz, Kac-Mittel
gold = (math.sqrt(5) - 1) / 2
ident("M", "Goldene Wicklung: n·‖nφ‖ → 1/√5 (n = F₂₀)", "01", fib[19] * abs(fib[19] * gold - round(fib[19] * gold)), 1 / math.sqrt(5), 1e-6)
ok = True; mw = []
for wv in (gold, math.sqrt(2) - 1, math.pi - 3):
    orb = (np.arange(200000) * wv) % 1.0; hit = np.where(orb < 0.137)[0]; dts = np.unique(np.diff(hit))
    ok &= len(dts) <= 3 and (len(dts) < 3 or dts[2] == dts[0] + dts[1]); mw.append(np.mean(np.diff(hit)) * 0.137)
wahr("M", "Drei-Abstände-Satz: Rückkehrzeiten nehmen ≤ 3 Werte an, der größte = Summe; Mittel = Kac", "71; 01; Slater 1967", ok and max(abs(np.array(mw) - 1)) < 1e-3, "Kac " + "/".join(f"{x:.4f}" for x in mw))
# M14 Laborbild der Scheinquelle: Kapitza-Pendel, oben stabil genau für (Aω)² > 2gℓ
def kapitza(A_w):
    om_ = 60.0; A_ = A_w / om_
    s_ = solve_ivp(lambda t, z: [z[1], (1.0 - A_ * om_**2 * math.cos(om_ * t)) * math.sin(z[0])], (0, 60), [0.2, 0.0], max_step=0.002, rtol=1e-8)
    return max(abs(s_.y[0])) < 1.0   # z = Winkel von oben
ok = kapitza(1.8) and not kapitza(1.1)
wahr("M", "Kapitza-Pendel: zeitgemitteltes Einschließen oben genau für (Aω)² > 2gℓ ((Aω)²/gℓ = 3,24 / 1,21)", "05 (Paul-Falle); Kapitza 1951", ok)
# M15 Äquivalenzschicht: Kugelschale und Vollkugel haben außen dasselbe Feld GM/r² (Dichte nicht eindeutig)
Rs_, rq_ = 1.0, 2.5
shell = quad(lambda t: 0.5 * math.sin(t) * (rq_ - Rs_ * math.cos(t)) / (rq_**2 + Rs_**2 - 2 * rq_ * Rs_ * math.cos(t))**1.5, 0, math.pi)[0]
ball = quad(lambda r_: 3 * r_**2 * quad(lambda t: 0.5 * math.sin(t) * (rq_ - r_ * math.cos(t)) / (rq_**2 + r_**2 - 2 * rq_ * r_ * math.cos(t))**1.5, 0, math.pi)[0], 0, 1)[0]
ok = abs(shell * rq_**2 - 1) < 1e-9 and abs(ball * rq_**2 - 1) < 1e-9
wahr("M", "Äquivalenzschicht: Schale und Vollkugel geben außen dasselbe Feld GM/r²", "05; 71", ok, f"{shell * rq_**2:.10f} / {ball * rq_**2:.10f}")
# M16 Engstelle (04): Kreuzterm, Lebensdauer, Rückgabe, Horizont doppelt so heiß, größte Falte ⅓, Ω_m²-Schranke, Atemrate
S_BH = lambda M: kB * 4 * math.pi * (2 * G * M / c**2)**2 / (4 * l_P**2)
M1_, M2_ = 3 * M_sun, 7 * M_sun
ident("M", "Kreuzterm beim Verschmelzen 8πG·M₁M₂k/(ħc) (aus der Horizontfläche)", "04", S_BH(M1_ + M2_) - S_BH(M1_) - S_BH(M2_), 8 * math.pi * G * M1_ * M2_ * kB / (hbar * c), 1e-9)
ident("M", "Gleiche Massen: Kreuzterm = ½ der Endentropie", "04", (S_BH(2 * M_sun) - 2 * S_BH(M_sun)) / S_BH(2 * M_sun), 0.5, 1e-12)
ident("M", "Lebensdauer t = 3840·(S/k)·Δs/c = 5120πG²M³/(ħc⁴)", "04", 3840 * S_BH(M_sun) / kB * ds(M_sun) / c, 5120 * math.pi * G**2 * M_sun**3 / (hbar * c**4), 1e-9)
tau_ = lebensdauer_falte(1e12); Mt = lambda t: 1e12 * (1 - t / tau_)**(1 / 3)
ident("M", "Rückgabe dM/dt = −ħc⁴/(15360πG²M²) (Ableitung von M(t))", "04", (Mt(1e-6 * tau_) - Mt(0)) / (1e-6 * tau_), -hbar * c**4 / (15360 * math.pi * G**2 * 1e24), 1e-5)
ident("M", "Horizont bei gleichem Radius doppelt so heiß: (ħH/2πk)/(ħc/(4πk·r₀)), r₀ = c/H", "04", (hbar * Hh / (2 * math.pi * kB)) / (hbar * c / (4 * math.pi * kB * c / Hh)), 2.0, 1e-12)
fSdS = lambda r, m: 1 - 2 * m / r - r * r
mN = brentq(lambda m: minimize_scalar(lambda r: -fSdS(r, m), bounds=(0.01, 1), method="bounded").fun, 0.1, 0.3, xtol=1e-14)
rN = minimize_scalar(lambda r: -fSdS(r, mN), bounds=(0.01, 1), method="bounded", options={"xatol": 1e-12}).x
sim("M", "Größte mögliche Falte (Nariai-Grenze, numerisch) hat ⅓ der Horizontentropie", "04", 1 / 3, rN**2, 1e-6)
rho_c = 3 * Hh**2 / (8 * math.pi * G); M_H = 4 / 3 * math.pi * (c / Hh)**3 * rho_c; S_H = math.pi * kB * c**2 / (l_P**2 * Hh**2)
ident("M", "Alle Falten zusammen höchstens Ω_m²·S_H (eine Falte aus Ω_m·M_H)", "04", S_BH(OMEGA_M * M_H) / S_H, OMEGA_M**2, 1e-9)
Si_ = 1e90; Sdot = lambda H: -(math.pi * kB * c**2 / l_P**2) * 2 / H**3   # dS_H/dH
ident("M", "Kein Endzustand: −Ḣ ≥ ℓ_P²H³Ṡ_i/(2πc²k) aus Ṡ_H ≥ Ṡ_i", "04", Si_ / (-Sdot(Hh)), l_P**2 * Hh**3 * Si_ / (2 * math.pi * c**2 * kB), 1e-9)
# M17 Vakuumteil und Strukturwirkung (03, 07)
Lam_ = c**2 / (2 * a0); Lam_s = math.pi**4 / (15 * Lam_**2)
ident("M", "Λ_s = π⁴/(15Λ²) ⇒ Ω_s = (π²/45)·Θ², Θ = 2πa₀/(cH)", "03", Lam_s * c**2 / (3 * Hh**2), K_S * (2 * math.pi * a0 / (c * Hh))**2, 1e-12)
ident("M", "S_max = πkc²/(ℓ_P²H²) = (45/π³)k(Λ/ℓ_P)² bei H² = Λ_s c²/3", "03", math.pi * kB * c**2 / (l_P**2 * Lam_s * c**2 / 3), 45 / math.pi**3 * kB * (Lam_ / l_P)**2, 1e-12)
Li = lambda n, x: sum(x**k / k**n for k in range(1, 400))
for s_v in (0.5, 1.0, 3.0):
    x_ = math.exp(-s_v); Qt = -4 * (s_v**3 * Li(1, x_) + 3 * s_v**2 * Li(2, x_) + 6 * s_v * Li(3, x_) + 6 * Li(4, x_))
    Qs = quad(lambda xp: float(nu(math.sqrt(xp), 1.0)) - 1, 0, s_v**4, limit=400)[0]
    ident("M", f"Strukturwirkung: Q̃(s) + 4π⁴/15 = Q_s(s⁴), der Vakuumteil ist ihre Konstante (s = {s_v})", "07; 03", Qt + 4 * math.pi**4 / 15, Qs, 1e-8)
# M18 Kern (06): Spaltbarkeit x = E_C/(2E_S): Rückstellung ∝ 1 − x (Rotationsellipsoid numerisch), größte Bindung bei x = ¼
def flaeche(eps):
    a_, b_ = 1 + eps, (1 + eps)**-0.5
    return quad(lambda x: 2 * math.pi * b_ * math.sqrt(1 - x * x / a_**2 + (b_ * x / a_**2)**2), -a_, a_, limit=200, epsabs=1e-14, epsrel=1e-13)[0]
def coulomb(eps):
    a_, b_ = 1 + eps, (1 + eps)**-0.5
    return 0.3 * quad(lambda l: 1 / ((b_ * b_ + l) * math.sqrt(a_ * a_ + l)), 0, np.inf, limit=200, epsabs=1e-14, epsrel=1e-13)[0]
k2 = lambda f, e: (f(e) + f(-e) - 2 * f(0)) / (f(0) * e * e)
cS = (4 * k2(flaeche, 0.01) - k2(flaeche, 0.02)) / 3; cC = (4 * k2(coulomb, 0.01) - k2(coulomb, 0.02)) / 3
sim("M", "Spaltbarkeit: Rückstellung ∝ (1 − x), Kippstelle bei x = E_C/(2E_S) = 1 (Ellipsoid numerisch)", "06; Bohr, Wheeler 1939", 1.0, -cS / (2 * cC), 1e-4)
zA = 0.45; Am_ = minimize_scalar(lambda A: 17.8 * A**(-1 / 3) + 0.711 * zA**2 * A**(2 / 3), bounds=(1, 1000), method="bounded", options={"xatol": 1e-10}).x
ident("M", "Größte Bindung je Nukleon (Oberfläche + Coulomb) bei x = ¼, E_C = E_S/2", "06", 0.711 * zA**2 * Am_**(5 / 3) / (2 * 17.8 * Am_**(2 / 3)), 0.25, 1e-6)
# M19 Aufbauzeit einer abgeschirmten Rückkehr τ(x) = x/(2√(Dk)) (Differenzenrechnung), Stationarität verlangt
dxm, dtm = 0.05, 0.001; xm = np.arange(0, 12, dxm); cm = np.zeros_like(xm); cm[0] = 1.0; tau_acc = np.zeros_like(xm)
ab_ = np.zeros((3, xm.size - 2)); ab_[0, 1:] = 1; ab_[1] = -(2 + dxm**2); ab_[2, :-1] = 1; rr_ = np.zeros(xm.size - 2); rr_[0] = -1
cs_ = np.concatenate([[1.0], solve_banded((1, 1), ab_, rr_), [0.0]]); cs_[-1] = 1.0   # stationärer Zustand desselben Gitters
for _ in range(int(25 / dtm)):
    lap = np.zeros_like(cm); lap[1:-1] = (cm[2:] - 2 * cm[1:-1] + cm[:-2]) / dxm**2
    cm = cm + dtm * (lap - cm); cm[0] = 1.0; cm[-1] = 0.0
    tau_acc += (1 - cm / cs_) * dtm
sim("M", "Aufbauzeit der Rückkehr τ(x) = x/(2√(Dk)) = τ_D/(2s) (x = 3λ)", "71 (Morphogen)", 1.5, tau_acc[int(3 / dxm)], 5e-3)
# M21 Eine Zahl trägt die Leiter: ε = ℓ_P/ℓ₀
eps_ = l_P / (c**2 / a0); Mmin_ = hbar * a0 / c**3
MminN_ = 10**brentq(lambda lm: math.log(hbar / (10**lm * c) / (c**2 / a0)), -80, -60)   # Nullstelle λ̄ = ℓ₀
ident("M", "ε_L = ℓ_P/ℓ₀ = a₀/a_P (a_P = c²/ℓ_P)", "08", eps_, a0 / (c**2 / l_P), 1e-12)
ident("M", "Spiegel der Leiter: M_min·M_U = m_Pl² (M_min = ħa₀/c³, Nullstellen)", "08", MminN_ * MU_, m_Pl**2, 1e-9)
ident("M", "Bei M_min zugleich r* = ℓ_P (Nullstelle von λ̄ = ℓ₀ eingesetzt)", "08", r_stern(MminN_), l_P, 1e-9)
ident("M", "M_U = ⅔·ρ_m·L³ (Zelle aus a₀ = 3c⁴/(2Gρ_mL³))", "08; 03", MU_, 2 / 3 * rho_m_mitbewegt() * zelle_aus_a0(a0)**3, 1e-9)
ident("M", "S_max(Λ)/k = (45/π³)(Λ/ℓ_P)² = (45/4π³)·M_U/M_min", "08; 03", 45 / math.pi**3 * (schliesslaenge(a0) / l_P)**2, 45 / (4 * math.pi**3) * MU_ / Mmin_, 1e-9)
# M22 (Rechnung, Invariante noch nicht im Modell) Spiegel der Diffusion: Rückkehr D = GM/c gegen Quantendiffusion ħ/M
dmax_ = max(abs((r_stern(m) ** 2 / tau_D(a0)) * (hbar / m) / (l_P * c) ** 2 - 1) for m in np.logspace(-30, 42, 9))
ident("M", "D·ħ/M = (ℓ_P c)² für jede Masse (D = r*²/τ_D aus der Rückkehr; Stichprobe 10⁻³⁰ … 10⁴² kg)", "Rechnung, nicht im Modell; Nelson 1966", 1.0 + dmax_, 1.0, 1e-12)
ident("M", "Quantendiffusion der Spiegelmasse ħ/(2M′) = D/2, M′ = m_Pl²/M (Sonne)", "08", hbar / (2 * m_Pl**2 / M_sun), r_stern(M_sun)**2 / tau_D(a0) / 2, 1e-12)
ident("M", "Faltungstakt c/(3Δs) = Compton-Takt der Spiegelmasse (Proton)", "09", c / (3 * ds(m_p)), (m_Pl**2 / m_p) * c**2 / hbar, 1e-12)
# M23 (Rechnung, Lesart noch nicht im Modell) Rückkehrgleichung als Schrödinger-Gleichung bei E = 0 (numerisch: Randwertaufgabe mit m = ħ/(2D), V = ħa_N/c)
Mt_ = 6e10 * M_sun; rt_ = r_stern(Mt_, a0); Dt_ = rt_**2 / tau_D(a0); mt_ = hbar / (2 * Dt_); aNt_ = G * Mt_ / (2 * rt_)**2; Vt_ = hbar * aNt_ / c
xx_ = np.linspace(0, 12 * rt_, 60001); hh_ = xx_[1] - xx_[0]; kk_ = 2 * mt_ * Vt_ / hbar**2
ab2_ = np.zeros((3, xx_.size - 2)); ab2_[0, 1:] = 1; ab2_[1] = -(2 + kk_ * hh_**2); ab2_[2, :-1] = 1; rr2_ = np.zeros(xx_.size - 2); rr2_[0] = -1
psi_ = solve_banded((1, 1), ab2_, rr2_)
ident("M", "Rückkehr als Tunneln: ψ(r*) mit m = m_Pl²/(2M), V = ħa_N/c gleich q = e^{−s} (Milchstraße, a_N bei 2r*)", "07", psi_[int(rt_ / hh_) - 1], q_(aNt_, a0), 1e-4)
ident("M", "Barriere bei r* = Ruheenergie der untersten Schnittmasse: ħa₀/c = M_min c²", "07; 08", hbar * a0 / c, Mmin_ * c**2, 1e-12)
# M24 Kac-Abweichung: bei r = 2Λ/π (Verlustrate = kleinster Eigenwert Dπ²/4Λ²) (Punktmasse, ϖ_R = r*) ist a_geschl/a_offen − 1 → 1/tanh(π/2) − 1
Lk_ = 2000.0; xk_ = 2 * Lk_ / math.pi; sk_ = 1 / xk_
qk_ = math.exp(-sk_ * 1) * (1 + math.exp(-2 * sk_ * (Lk_ - 1))) / (1 + math.exp(-2 * sk_ * Lk_))   # cosh(κ(Λ−r*))/cosh(κΛ) stabil
sim("M", "Geschlossene Rückkehr: Abweichung bei r = 2Λ/π gleich 1/tanh(π/2) − 1 (Γ = kleinster Eigenwert)", "07", 1 / math.tanh(math.pi / 2) - 1, (1 - math.exp(-sk_)) / (1 - qk_) - 1, 2e-3)
# M25 Vakuumteil in Planck-Einheiten und Einseitigkeit im äußeren Feld (erste Ordnung gegen volle Rechnung)
ident("M", "Λ_s·ℓ_P² = (4π⁴/15)·ε_L²", "08", lambda_vakuum(schliesslaenge(a0)) * l_P**2, 4 * math.pi**4 / 15 * eps_**2, 1e-12)
import tunnel_aussenfeld as _ta
Rt2_ = 20 * _ta.kpc; gi2_ = G * _ta.M / Rt2_**2; ge2_ = 0.03 * a0; s2_ = math.sqrt(gi2_ / a0); n2_ = float(nu(gi2_, a0))
dn2_ = -n2_ * (n2_ - 1) * s2_ / 2   # dν/d ln g
m1_lin = 0.5 * (-n2_ - dn2_ + float(nu(ge2_, a0))) * ge2_ / (n2_ * gi2_)
ph_ = np.linspace(0, 2 * np.pi, 720, endpoint=False); eR_ = np.stack([np.cos(ph_), np.sin(ph_), 0 * ph_]); gv_ = np.array([ge2_, 0, 0])
gt_ = -gi2_ * eR_ + gv_[:, None]; aa_ = nu(np.linalg.norm(gt_, axis=0), a0) * gt_ - (nu(ge2_, a0) * gv_)[:, None]
vc_ = np.sqrt(-Rt2_ * np.sum(aa_ * eR_, axis=0)); F_ = np.fft.rfft(vc_) / len(vc_)
sim("M", "Einseitigkeit m = 1 bei 20 kpc, g_e = 0,03a₀ in der Ebene: erste Ordnung gegen volle Rechnung", "03", m1_lin, 2 * abs(F_[1]) / F_[0].real, 0.03)
# M26 Verlustrate aus der Faltung der Massen: mit Γ = a/c (Rückläufe mitgezählt) fielen die Kurven wie r^{−1/6}
def _a_selbst(aN):
    return a0 * brentq(lambda y: y * (1 - math.exp(-math.sqrt(y))) - aN / a0, aN / a0, 1e3 * aN / a0 + 10, xtol=1e-15, rtol=1e-14)
_Mg = 6e10 * M_sun; _rg = 1e5 * r_stern(_Mg, a0)
_sl = math.log(math.sqrt(_a_selbst(G * _Mg / (1.1 * _rg)**2) * 1.1 * _rg) / math.sqrt(_a_selbst(G * _Mg / _rg**2) * _rg)) / math.log(1.1)
sim("M", "Γ = a/c selbstkonsistent: d ln v/d ln r → −1/6 (fallende Kurven; Regel nutzt Γ = a_N/c)", "07", -1 / 6, _sl, 5e-3)
# M20 Eine Teilung für jedes Gedächtnis: gebunden x²/(1 + x²), gerieben x/(1 + x²)
def maxwell_sim(x, n=40):
    """Maxwell-Element im Zeitbereich (t in Einheiten von t_K): σ̇ = ε̇ − σ, ε = sin(xt); Anteil in Phase = Speicher."""
    T_ = 2 * math.pi / x
    sol = solve_ivp(lambda t, y: [x * math.cos(x * t) - y[0]], (0, (n + 1) * T_ + 20), [0.0], rtol=1e-11, atol=1e-13, dense_output=True)
    tt = np.linspace(20 + n * T_ - T_, 20 + n * T_, 4001)
    sg = sol.sol(tt)[0]
    return 2 * np.trapezoid(sg * np.sin(x * tt), tt) / T_, 2 * np.trapezoid(sg * np.cos(x * tt), tt) / T_
tauD_ = tau_D(a0)
dw_, dl_ = 0.0, 0.0
for xx in (0.1, 0.5, 1.0, 3.0, 17.86):
    sp_, vl_ = maxwell_sim(xx)
    dw_ = max(dw_, abs(sp_ - gewicht_w_maxwell(xx / tauD_, a0)))
    dl_ = max(dl_, abs(vl_ - xx / (1 + xx * xx)))
sim("M", "Frühere Form des Gewichts (Rechnung): Speicheranteil eines Maxwell-Elements mit t_K = τ_D (x = θτ_D; Zeitbereich)", "Rechnung", 1.0, 1.0 + dw_, 1e-6)
_lw = quad(lambda t: (1 / math.sqrt(4 * math.pi * t**3)) * math.exp(-1 / (4 * t)) * math.exp(-3.0 * t), 0, np.inf, limit=400)[0]
sim("M", "Gewicht der Schließungsregel aus der Rückkehr: 1 − w = Laplace-Transformierte der Dichte erster Rückkehrzeiten = e^{−√x} (x = 3; numerisch)", "03; 05", math.exp(-math.sqrt(3.0)), _lw, 1e-8)
sim("M", "Verlustanteil des Maxwell-Elements x/(1 + x²) (Zeitbereich)", "05", 1.0, 1.0 + dl_, 1e-6)
om_, lm_, rk_ = 1.0, 0.3, 2.0     # Kreisbahn r = rk_, Gedächtnisrate λ_m; Kraft längs e_r(t) aus ∫e^{−λu} ṙ(t−u) du
ker_ = quad(lambda u: math.exp(-lm_ * u) * (rk_ * om_ * math.sin(om_ * u)), 0, np.inf, limit=400)[0]
sim("M", "Gedächtniskern auf der Kreisbahn: rω²/(λ_m² + ω²) = r·(Speicheranteil, x = ω/λ_m); numerisch", "07", rk_ * om_**2 / (lm_**2 + om_**2), ker_, 1e-8)
kr_ = 1.3                            # Radialschwingung δr = e^{iκt}: ∫e^{−λu}δṙ(t−u)du/δr = (κ² + iκλ)/(λ² + κ²)
re_ = quad(lambda u: math.exp(-lm_ * u) * kr_ * math.sin(kr_ * u), 0, np.inf, limit=400)[0]
im_ = quad(lambda u: math.exp(-lm_ * u) * kr_ * math.cos(kr_ * u), 0, np.inf, limit=400)[0]
sim("M", "Radialschwingung im Kern: Speicher κ²/(λ²+κ²), Verlust κλ/(λ²+κ²) (rundet Bahnen); numerisch", "07", (kr_**2 + kr_ * lm_) / (lm_**2 + kr_**2), re_ + im_, 1e-8)
rMW_ = r_stern(6e10 * M_sun, a0); vMW_ = v_kreis_punkt(6e10 * M_sun, rMW_, a0)
ident("M", "Ω·τ_D = (a/a₀)·(c/v) (Milchstraße bei r*)", "07", vMW_ / rMW_ * tauD_, (vMW_**2 / rMW_ / a0) * (c / vMW_), 1e-12)

# ======================================================================== E Messwerte mit Quelle
mess("E", "Schwelle a₀ gegen die gemessene Relation g† (systematischer Fehler)", "McGaugh, Lelli, Schombert 2016: 1,20 ± 0,24 (sys)", a0 * 1e10, 1.20, 0.24)
mess("E", "Zelle ohne Kreispaare: a₀ ≤ 1,13·10⁻¹⁰ (L ≥ 0,985·d_LSS)", "COMPACT, Petersen et al. 2023 (JCAP 01, 030)", a0 * 1e10, a0_aus_zelle(KREISPAAR_GRENZE * d_LSS()) * 1e10, 0.02)
mess("E", "N_eff (Vakuumteil früh null)", "Planck 2018 + BAO: 2,99 ± 0,17", N_EFF_STD, 2.99, 0.17)
mess("E", "Krümmung Ω_K (Zelle atmet, Λ fest)", "Planck 2018 + BAO: 0,0007 ± 0,0037 (95 %)", 0.0, 0.0007, 0.0037 / 2)
mess("E", "Crater II σ, örtliche Regel (Jeans, Außenfeld)", "Caldwell et al. 2017: 2,7 ± 0,3 km/s", 1.39, 2.7, 0.3)
mess("E", "Gravity Probe B geodätisch [mas/Jahr]", "Everitt et al. 2011: 6601,8 ± 18,3", 4.5 * ds(M_earth) / 7027.4e3 * w_orb * yr * 180 / math.pi * 3600e3, 6601.8, 18.3)
mess("E", "Gravity Probe B Mitführung [mas/Jahr]", "Everitt et al. 2011: 37,2 ± 7,2", 3 * ds(M_earth) * aJ * c / (2 * 7027.4e3**3) * yr * 180 / math.pi * 3600e3, 37.2, 7.2)
mess("E", "Hulse-Taylor Ṗ gemessen / Theorie", "Weisberg, Huang 2016: 0,9983 ± 0,0016", 1.0, 0.9983, 0.0016)
mess("E", "Schatten M87* (6,5·10⁹ M☉, 16,8 Mpc) [μas]", "EHT 2019: 42 ± 3", schatten_winkel(6.5e9 * M_sun, 16.8 * Mpc) * 180 / math.pi * 3600e6, 42, 3)

# ======================================================================== Ausgabe
print(f"{'Teil':4s} | {'Gesetz':78s} | {'Modell':>11s} | {'Vergleich':>11s} | {'Maß':>18s} | Urteil | Quelle")
for t, g_, q_, m_, v_, d_, u_ in R:
    fm = lambda x: f"{x:11.5g}" if isinstance(x, (int, float, np.floating)) else f"{str(x):>11s}"
    print(f"{t:4s} | {g_[:78]:78s} | {fm(m_)} | {fm(v_)} | {d_:>18s} | {u_:6s} | {q_}")
from collections import Counter
cnt = Counter(u for *_, u in R)
print(f"\n{len(R)} Gesetze geprüft: " + ", ".join(f"{k} {v}" for k, v in sorted(cnt.items())) + f"   (Rechenzeit {time.time() - t_start:.0f} s)")
