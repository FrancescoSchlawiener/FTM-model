# -*- coding: utf-8 -*-
"""
FTM - Forschungsnotiz 26: Freiheitsgrad-Einschraenkung statt Zusatzquelle.

Frage: Kann der Bindungsoperator des Autors (Gedaechtnis bzw. Projektion (I - Pi_perp)),
der KEINE neue Substanz einfuehrt, im CMB das leisten, was kalte Dunkle Materie leistet:
ein Potential, das nicht mit dem Photon-Baryon-Fluid schwingt?

Abschnitte der Ausgabe:
 1   Hintergrund, Schallhorizont, Daempfungsskala
 1a  Analytischer Oszillator (Hu & Sugiyama): Nulllage und Ungerade/Gerade
 1b  Spielzeug-Boltzmann (lineare Stoerungen, Newton-Eichung, enge Kopplung):
     (a) LCDM-artig, (b) nur Baryonen, (c) verstaerkte Baryonenkraft (mitschwingend)
 2a  Gedaechtnis als Tiefpass: Uebertragungsfunktion
 2b  Variante (i) Gedaechtnis (Rueckkopplungsform, Kern exponentiell in e-Folds, Rate beta)
 2c  Variante (ii) Projektion (I - w Pi_perp) auf den Baryonenfluss; (ii') Momentan-Projektion
 2d  Noetige Staerke; q zu 5,36; Vergleich mit y-Werten aus Notiz 08
 3   Galaxien: Operator allein gegen Hybrid (Notiz 17); konstantes q = 0,843 gegen RAR
 4   Wachstum nach der Rekombination (lineare Rechnung bis z = 0)

Nur numpy/math. Aufruf: python freiheitsgrad.py  (schreibt freiheitsgrad.txt daneben)
Konventionen: Ma & Bertschinger (1995), konforme Newton-Eichung, psi = phi (keine
Neutrino-Scherung), Neutrinos als ideales Fluid, enge Photon-Baryon-Kopplung bis z* = 1090.
Das ist ein Spielzeugmodell: absolute Peakhoehen weichen von CLASS/CAMB ab; verglichen
werden die Varianten untereinander mit identischer Methode.
"""
import math
import os
import numpy as np

OUT = []


def P(s=""):
    OUT.append(str(s))


def hdr(s):
    P()
    P("=" * 100)
    P(s)
    P("=" * 100)


# ------------------------------------------------------------------ Eingaben (wie Notiz 08)
h = 0.6736
wb, wc = 0.02237, 0.1200
wg, wr = 2.4729e-5, 4.1835e-5
Ob, Oc, Or = wb / h**2, wc / h**2, wr / h**2
fg = wg / wr                       # Photonenanteil der Strahlung
H0c = h / 2997.92458               # H0/c in 1/Mpc
zs = 1090.0
a_s = 1.0 / (1.0 + zs)
Q536 = wc / wb                     # 5,364
q536 = Q536 / (1.0 + Q536)         # 0,8429
YHe = 0.245
sigT = 6.6524587e-29
mp = 1.67262192e-27
rhoc_h2 = 1.87847e-26              # kg/m^3 pro h^2
MpcM = 3.0856776e22
nH0 = (1 - YHe) * wb * rhoc_h2 / mp            # 1/m^3
tau0 = nH0 * sigT * MpcM                       # n_e sigma_T heute in 1/Mpc (x_e = 1)


class BG:
    """Flacher Hintergrund mit staubartiger Gesamtmaterie Om (Baryonen + ggf. Operatoranteil)."""

    def __init__(self, Om):
        self.Om = Om
        self.OL = 1.0 - Om - Or

    def E(self, a):
        return np.sqrt(Or * a**-4 + self.Om * a**-3 + self.OL)

    def Hc(self, a):          # konformes H in 1/Mpc
        return a * H0c * self.E(a)

    def frac(self, a):        # Anteile an der kritischen Dichte zur Zeit a
        E2 = self.E(a) ** 2
        return Or * a**-4 / E2, Ob * a**-3 / E2, self.Om * a**-3 / E2, self.OL / E2

    def dlnE(self, a):
        E2 = self.E(a) ** 2
        return 0.5 * (-4 * Or * a**-4 - 3 * self.Om * a**-3) / E2

    def eta(self, a, n=4000):
        x = np.linspace(math.log(1e-10), math.log(a), n)
        aa = np.exp(x)
        f = 1.0 / self.Hc(aa)                      # d eta / d ln a
        return float(np.trapezoid(f, x) + 1e-10 / (H0c * math.sqrt(Or)))

    def rs(self, a):
        x = np.linspace(math.log(1e-10), math.log(a), 6000)
        aa = np.exp(x)
        R = 0.75 * wb / wg * aa
        f = 1.0 / np.sqrt(3 * (1 + R)) / self.Hc(aa)
        return float(np.trapezoid(f, x))

    def DM(self, a):
        x = np.linspace(math.log(a), 0.0, 6000)
        aa = np.exp(x)
        return float(np.trapezoid(1.0 / self.Hc(aa), x))

    def kD(self, a):
        x = np.linspace(math.log(1e-8), math.log(a), 6000)
        aa = np.exp(x)
        R = 0.75 * wb / wg * aa
        tdot = tau0 / aa**2                        # n_e sigma_T a, x_e = 1
        f = (R**2 / (1 + R) + 16.0 / 15.0) / (6 * (1 + R) * tdot) / self.Hc(aa)
        return 1.0 / math.sqrt(float(np.trapezoid(f, x)))


# ------------------------------------------------------------------ Spielzeug-Boltzmann
# Zustand: 0 phi, 1 d_g, 2 th (Photon-Baryon), 3 d_nu, 4 th_nu, 5 d_s, 6 th_s, 7 M
NV = 8


def build_A(a, k, bg, cfg):
    """Matrix A mit dy/dlna = A y (fuer alle k)."""
    Hc = bg.Hc(a)
    fr, fb, fm, fL = bg.frac(a)
    fgam, fnu = fr * fg, fr * (1 - fg)
    R = 0.75 * wb / wg * a
    nk = k.size
    A = np.zeros((nk, NV, NV))
    k2 = k * k
    Qenh = cfg.get("enhQ", 0.0)        # (c): Baryonenquelle (1+Q)-fach, schwingt mit
    wS = cfg.get("wS", 0.0)            # Gewicht der Schatten-/CDM-Dichte in Einheiten von fb
    wproj = cfg.get("w", 1.0)          # Projektionsgewicht (1 = Druckrichtung ganz inaktiv)
    mem = cfg.get("mem", False)
    beta, qm = cfg.get("beta", 1.0), cfg.get("q", 0.0)
    Qeq = cfg.get("eqQ", 0.0)          # (ii') Momentan-Projektion auf die Gleichgewichtslage
    # Quelle S = sum Omega_i delta_i
    S = np.zeros((nk, NV))
    S[:, 1] = fgam + 0.75 * fb * (1 + Qenh)
    S[:, 3] = fnu
    S[:, 5] = wS * fb
    if mem:
        S[:, 7] = fb
    S[:, 0] = Qeq * fb * (-3.0 * (1 + R))     # delta_eq = -3(1+R) psi, psi = phi
    # phi-Zeile (00-Gleichung):  dphi/dN = -phi - k^2 phi/(3 Hc^2) - S/2
    rowphi = -0.5 * S
    rowphi[:, 0] += -1.0 - k2 / (3 * Hc**2)
    A[:, 0, :] = rowphi
    # d_g' = -(4/3) th/Hc + 4 phi'
    A[:, 1, :] = 4 * rowphi
    A[:, 1, 2] += -4.0 / 3.0 / Hc
    # th' = -R/(1+R) th + k^2 d_g/(4(1+R)Hc) + k^2 phi/Hc
    A[:, 2, 2] = -R / (1 + R)
    A[:, 2, 1] = k2 / (4 * (1 + R) * Hc)
    A[:, 2, 0] = k2 / Hc
    # Neutrinos (ideales Fluid)
    A[:, 3, :] = 4 * rowphi
    A[:, 3, 4] += -4.0 / 3.0 / Hc
    A[:, 4, 3] = k2 / (4 * Hc)
    A[:, 4, 0] = k2 / Hc
    # Schatten/CDM: d_s' = -th_s/Hc + 3 phi';  th_s' = -th_s + k^2 phi/Hc + (1-w)[th/(1+R) + k^2 d_g/(4(1+R)Hc)]
    A[:, 5, :] = 3 * rowphi
    A[:, 5, 6] += -1.0 / Hc
    A[:, 6, 6] = -1.0
    A[:, 6, 0] = k2 / Hc
    A[:, 6, 2] += (1 - wproj) / (1 + R)
    A[:, 6, 1] += (1 - wproj) * k2 / (4 * (1 + R) * Hc)
    # Gedaechtnis (Rueckkopplungsform): dM/dN = beta [ q (d_b + M) - M ],  d_b = 0.75 d_g
    if mem:
        A[:, 7, 1] = beta * qm * 0.75
        A[:, 7, 7] = beta * (qm - 1.0)
    return A


def run(cfg, k, a_end=a_s, dN=0.0015, ai=1e-7):
    bg = BG(cfg["Om"])
    eta_i = bg.eta(ai)
    y = np.zeros((k.size, NV))
    y[:, 0] = 1.0
    y[:, 1] = -2.0
    y[:, 2] = 0.5 * k**2 * eta_i
    y[:, 3] = -2.0
    y[:, 4] = 0.5 * k**2 * eta_i
    y[:, 5] = -1.5
    y[:, 6] = 0.5 * k**2 * eta_i
    if cfg.get("mem"):
        qm = cfg["q"]
        y[:, 7] = qm / (1 - qm) * (-1.5)          # stationaer: M = Q d_b
    N = np.arange(math.log(ai), math.log(a_end), dN)
    N = np.append(N, math.log(a_end))
    I = np.eye(NV)[None, :, :]
    An = build_A(math.exp(N[0]), k, bg, cfg)
    for j in range(len(N) - 1):
        hN = N[j + 1] - N[j]
        An1 = build_A(math.exp(N[j + 1]), k, bg, cfg)
        rhs = y + 0.5 * hN * np.einsum("kij,kj->ki", An, y)
        y = np.linalg.solve(I - 0.5 * hN * An1, rhs[..., None])[..., 0]
        An = An1
    return y, bg


KD_PLANCK = 0.1409     # 1/Mpc, Planck-2018-Ableitung (aus dem Gedaechtnis) fuer den LCDM-Hintergrund


def spectrum(y, k, bg, kDfac=1.0):
    """Proxy fuer l(l+1)C_l: (Theta0+psi)^2 + Theta1^2 bei z*, gedaempft mit exp(-2 k^2/k_D^2).
    k_D: Planck-Wert fuer LCDM, fuer andere Hintergruende skaliert mit dem Verhaeltnis der x_e=1-Integrale."""
    Th = y[:, 1] / 4 + y[:, 0]
    Th1 = y[:, 2] / (3 * k)
    kD = KD_PLANCK * bg.kD(a_s) / BG(Ob + Oc).kD(a_s) * kDfac
    return (Th**2 + Th1**2) * np.exp(-2 * (k / kD) ** 2), kD


def peaks(Pk, k, n=5, kmin=0.008):
    idx = [i for i in range(1, len(k) - 1) if k[i] > kmin and Pk[i] > Pk[i - 1] and Pk[i] >= Pk[i + 1]]
    idx = idx[:n]
    return [(k[i], Pk[i]) for i in idx]


def ratios(pk):
    H = [p[1] for p in pk]
    while len(H) < 5:
        H.append(float("nan"))
    return H[1] / H[0], H[2] / H[0], H[2] / H[1], H[3] / H[2], H[4] / H[3]


kgrid = np.linspace(0.002, 0.32, 520)
results = {}


def do(name, cfg, show=True):
    y, bg = run(cfg, kgrid)
    Pk, kD = spectrum(y, kgrid, bg)
    pk = peaks(Pk, kgrid)
    results[name] = dict(y=y, bg=bg, Pk=Pk, pk=pk, cfg=cfg, kD=kD)
    return results[name]


def line(name, r, ref=None):
    rr = ratios(r["pk"])
    plateau = r["Pk"][0]
    s = f"  {name:44s} H2/H1 {rr[0]:6.3f}  H3/H1 {rr[1]:6.3f}  H3/H2 {rr[2]:6.3f}  H4/H3 {rr[3]:6.3f}  H5/H4 {rr[4]:6.3f}  H1/Plateau {r['pk'][0][1]/plateau:6.2f}"
    if ref is not None:
        ra = ratios(ref["pk"])
        s += f"  |dH3/H2| {abs(rr[2]-ra[2])/ra[2]*100:5.1f}%"
    P(s)
    return rr


# ==========================================================================================
hdr("1  Hintergrund, Schallhorizont, Daempfung (k_D-Spalte: Integral mit x_e = 1, zu gross; im Spektrum auf Planck k_D = 0,1409 fuer LCDM skaliert)")
for nm, Om in (("LCDM / Operator mit Q = 5,364 (Om = Ob(1+Q))", Ob + Oc), ("nur Baryonen", Ob)):
    bg = BG(Om)
    zeq = Om / Or - 1
    P(f"  {nm:48s} z_eq {zeq:7.0f}  r_s(z*) {bg.rs(a_s):7.2f} Mpc  D_M(z*) {bg.DM(a_s):8.1f} Mpc  "
      f"100 theta* {100*bg.rs(a_s)/bg.DM(a_s):6.4f}  k_D {bg.kD(a_s):6.4f}/Mpc  Hc(a*) {bg.Hc(a_s):.4e}/Mpc")
P(f"  Q = omega_c/omega_b = {Q536:.4f};  q = Q/(1+Q) = {q536:.5f};  1/(1-q) = {1/(1-q536):.4f}")
P("  Befund: Mit Om_eff = Ob(1+Q) ist der Hintergrund identisch mit LCDM. Das gilt fuer JEDE Variante,")
P("  in der der Operator auf die erhaltene mitbewegte Baryonenmasse wirkt (dort ist das Gedaechtnis exakt stationaer).")

# Gedaechtnis auf dem physikalischen Potential statt auf der mitbewegten Masse:
P()
P("  Gegenprobe: Gedaechtnis auf das physikalische Hintergrundpotential Phi_bg ~ 1/a (fester mitbewegter Radius),")
P("  Kern Q*beta*exp(-beta s) in e-Folds: wirksamer Faktor = Q*beta/(beta-1) (divergiert fuer beta <= 1):")
for b in (1.2, 2.0, 5.0, 20.0):
    P(f"    beta = {b:5.1f}:  Q_eff/Q = {b/(b-1):.3f}")
P("  => Ein langes Gedaechtnis muss auf der erhaltenen Groesse (mitbewegte Masse, Kontrast delta) sitzen, nicht auf Phi.")

# ==========================================================================================
hdr("1a  Analytischer Oszillator mit festem Psi (Hu & Sugiyama 1995):  Theta0+Psi = (1/3+R)Psi cos(k r_s) - R Psi")
Rst = 0.75 * wb / wg * a_s
P(f"  R(z*) = {Rst:.4f}")
odd, even = (1/3 + 2*Rst), (1/3)
P(f"  statische Mulde: |ungerade| = 1/3+2R = {odd:.4f}, |gerade| = 1/3 = {even:.4f}, Amplitudenverh. 1+6R = {odd/even:.3f}, Leistung {(odd/even)**2:.2f}")
P("  Mitschwingende Mulde (Psi ~ Theta, zeitlicher Mittelwert 0): Nulllage -R<Psi> = 0 -> Verhaeltnis 1 (keine Ungerade/Gerade-Modulation)")
P("  Das ist die Freiheitsgrad-Aussage: nur ein Potentialanteil, der NICHT Teil der akustischen Mode ist, verschiebt die Nulllage.")

# ==========================================================================================
hdr("1b  Spielzeug-Boltzmann: Peakhoehen im k-Raum (Proxy fuer l(l+1)C_l), Referenz Notiz 08 (Fitformel, Planck-Werte): H2/H1 0,460  H3/H1 0,450  H3/H2 0,978")
ra = do("a", dict(Om=Ob + Oc, wS=Oc / Ob, w=1.0))
rb = do("b", dict(Om=Ob))
rc = do("c", dict(Om=Ob * (1 + Q536), enhQ=Q536))
rrA = line("(a) LCDM-artig (Mulde aus 5,36-fach, statisch)", ra)
line("(b) nur Baryonen", rb, ra)
line("(c) Baryonenkraft x6,36 (Mulde schwingt mit)", rc, ra)
P()
P("  Peaklagen k [1/Mpc] und l ~ k D_M:")
for nm in ("a", "b", "c"):
    r = results[nm]
    DM = r["bg"].DM(a_s)
    P(f"   ({nm}) " + "  ".join(f"k={p[0]:.4f} (l~{p[0]*DM:5.0f})" for p in r["pk"]))
P()
P("  Robustheit (k_D x 0,8, staerkere Daempfung), nur H3/H2:")
for nm in ("a", "b", "c"):
    r = results[nm]
    Pk2, _ = spectrum(r["y"], kgrid, r["bg"], 0.8)
    P(f"   ({nm}) H3/H2 = {ratios(peaks(Pk2, kgrid))[2]:.3f}")

# Potentialanteile bei z*: gebunden (nicht schwingend) gegen frei (Baryonen)
P()
P("  Anteile an der Poisson-Quelle bei z* (Einheiten Ob(a*)*delta; Baryonen: Huellkurve |d_b| ueber +-0,01/Mpc):")
fr, fb, fm, fL = ra["bg"].frac(a_s)
for kk in (0.02, 0.05, 0.1, 0.2):
    i = int(np.argmin(abs(kgrid - kk)))
    win = (abs(kgrid - kk) < 0.01)
    db_env = np.max(np.abs(0.75 * ra["y"][win, 1]))
    P(f"   k = {kk:4.2f}: (a) Q*d_c = {Q536*ra['y'][i,5]:9.2f},  |d_b|max = {db_env:7.2f},  Verh. gebunden/frei = {abs(Q536*ra['y'][i,5])/db_env:6.2f}")

# ==========================================================================================
hdr("2a  Gedaechtnis als Tiefpass: M = int K (x + M), K = q lam e^{-lam s}  =>  M^ = q lam/(lam(1-q) + i w) x^")
P("  stationaer (w = 0): M = q/(1-q) x = Q x  (bei q = 0,8429: Q = 5,364). Schwingung (w >> lam(1-q)): |M^|/x^ -> q lam/w.")
P("  In e-Folds: lam = beta*H, Schwingungsrate w_N = k c_s / Hc(a*):")
bgA = BG(Ob + Oc)
cs = 1 / math.sqrt(3 * (1 + Rst))
P("   k[1/Mpc]  w_N    |H|/Q fuer beta = 0,03   0,1    0,3    1      3      10")
for kk in (0.02, 0.05, 0.1, 0.2):
    wN = kk * cs / bgA.Hc(a_s)
    vals = []
    for b in (0.03, 0.1, 0.3, 1, 3, 10):
        lam1 = b * (1 - q536)
        vals.append(lam1 / math.sqrt(lam1**2 + wN**2))
    P(f"   {kk:6.2f}   {wN:6.2f}   " + "  ".join(f"{v:6.3f}" for v in vals))
P("  => Fuer beta(1-q) << w_N wird die Schwingung herausgefiltert, der statische Anteil bleibt voll (Q).")

# ==========================================================================================
hdr("2b  Variante (i) Gedaechtnis auf den mitbewegten Baryonenkontrast, q = 0,8429, Hintergrund Ob(1+Q)")
for b in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 30.0):
    r = do(f"i{b}", dict(Om=Ob * (1 + Q536), mem=True, beta=b, q=q536))
    line(f"(i) beta = {b:g} (Gedaechtnis {1/(b*(1-q536)):.2f} e-Folds)", r, ra)
P()
P("  Gedaechtnisanteil M bei z* gegen (a) Q*d_c (Einheiten phi_p):")
P("   k      (a) Q d_c   (i) M: beta=0,01   0,1     1      30")
for kk in (0.02, 0.05, 0.1, 0.2):
    i = int(np.argmin(abs(kgrid - kk)))
    P(f"   {kk:4.2f}   {Q536*ra['y'][i,5]:9.2f}      " + "  ".join(f"{results[f'i{b}']['y'][i,7]:7.2f}" for b in (0.01, 0.1, 1.0, 30.0)))

# ==========================================================================================
hdr("2c  Variante (ii) Projektion z' = (I - w Pi_perp) F auf den Baryonenfluss; Pi_perp = Druck-/Streurichtung (Flucht aus der Mulde)")
P("  Der projizierte Fluss hat keine eigenen Anfangsdaten: Start = Baryonenzustand (adiabatisch von selbst). Gewicht Q = 5,364.")
for w in (0.0, 0.25, 0.5, 0.75, 0.9, 0.97, 1.0):
    r = do(f"ii{w}", dict(Om=Ob * (1 + Q536), wS=Q536, w=w))
    line(f"(ii) w = {w:g}", r, ra)
P("  w = 1 ist Zeile fuer Zeile (a): die Gleichungen sind identisch mit CDM der Dichte Q*rho_b.")
P()
r = do("iip", dict(Om=Ob * (1 + Q536), eqQ=Q536))
line("(ii') Momentan-Projektion auf Gleichgewichtslage", r, ra)
P("  (ii') projiziert den momentanen Baryonenzustand auf die nicht schwingende Richtung: delta_eq = -3(1+R) psi.")
i = int(np.argmin(abs(kgrid - 0.05)))
P(f"   Beispiel k = 0,05: |Q delta_eq| = {abs(Q536*3*(1+Rst)*r['y'][i,0]):.3f}  gegen (a) |Q d_c| = {abs(Q536*ra['y'][i,5]):.2f}")
P(f"   Plateau (k = 0,002): (ii') {r['Pk'][0]:.3f} gegen (a) {ra['Pk'][0]:.3f}; Rueckkopplungsfaktor (9/2)(Hc/k)^2 Ob(a*) Q (1+R) bei k = 0,002: "
  f"{4.5*(ra['bg'].Hc(a_s)/0.002)**2*ra['bg'].frac(a_s)[1]*Q536*(1+Rst):.1f}, bei k = 0,05: {4.5*(ra['bg'].Hc(a_s)/0.05)**2*ra['bg'].frac(a_s)[1]*Q536*(1+Rst):.4f}")

# ==========================================================================================
hdr("2d  Noetige Staerke. (1) Gedaechtnis mit Stoerungs-Q frei, Hintergrund fest bei 5,364; (2) q zu 5,36; (3) y-Werte aus Notiz 08")
P("  (1) beta = 0,03 und 0,1; Q_pert gescannt:")
best = {}
for b in (0.03, 0.1):
    for Qp in (5.364, 10, 20, 40, 80):
        qp = Qp / (1 + Qp)
        r = do(f"Q{b}_{Qp}", dict(Om=Ob * (1 + Q536), mem=True, beta=b, q=qp))
        rr = line(f"beta={b:g}, Q_pert={Qp:g}", r, ra)
        best[(b, Qp)] = rr
P()
P("  (2) Rueckkopplungssumme: gebunden/frei = sum_{n>=1} q^n = q/(1-q) = 5,364  =>")
y_need = math.log(1 / q536) ** 2
P(f"      q = {q536:.5f};  mit q = exp(-sqrt(y)):  y = {y_need:.5f},  nu = 1/(1-q) = {1/(1-q536):.4f}")
P("  (3) y-Werte der Stoerungen bei z = 1100 (Formeln aus Notiz 08) und das zugehoerige Q = q/(1-q):")
c_ = 2.99792458e8
H1100 = 100 * h * 1e3 / MpcM * BG(Ob + Oc).E(1 / 1101)
for Rc, Phi in ((147, 1e-5), (50, 1e-5), (10, 1e-5), (1, 1e-4)):
    Rphys = Rc * MpcM / 1101
    yv = 2 * math.pi * Phi * c_ / (Rphys * H1100)
    qv = math.exp(-math.sqrt(yv))
    P(f"      R = {Rc:4d} Mpc, Phi/c^2 = {Phi:.0e}:  y = {yv:.3e},  q = {qv:.5f},  Q = {qv/(1-qv):9.1f}")
fr, fb, fm, fL = BG(Ob + Oc).frac(1 / 1101)
for nm, yv in (("Hubble-Kugel y = pi", math.pi), (f"Hintergrund y = pi*Omega_b(1100) = pi*{fb:.4f}", math.pi * fb)):
    qv = math.exp(-math.sqrt(yv))
    P(f"      {nm:44s}: y = {yv:.4f},  q = {qv:.4f},  Q = {qv/(1-qv):.3f}")
P(f"      Noetige Bezugslaenge fuer y_bg = pi Omega_b(z)/X = {y_need:.4f}:  X = {math.pi*fb/y_need:.2f} (bei z = 1100), X = {math.pi*(1-q536)/y_need:.2f} (Materie-Aera, Omega_b = 1-q)")

# ==========================================================================================
hdr("3  Galaxien: nur der Operator (keine Quelle) gegen Hybrid (Schatten klumpt wie CDM, Notiz 17)")
G = 6.674e-11; Msun = 1.989e30; kpc = 3.0857e19
a0 = 1.2e-10
H0s = 100 * h * 1e3 / MpcM
rho_c0 = 3 * H0s**2 / (8 * math.pi * G)


def nu_ftm(y):
    return 1.0 / (1.0 - math.exp(-math.sqrt(y)))


def nfw_M(r, M200, cc):
    r200 = (3 * M200 / (4 * math.pi * 200 * rho_c0)) ** (1 / 3)
    rs_ = r200 / cc
    m = lambda x: math.log(1 + x) - x / (1 + x)
    return M200 * m(r / rs_) / m(cc)


def hern_M(r, M, a):
    return M * r * r / (r + a) ** 2


P("  Baryonen Hernquist a = 1,5 kpc (M/1e10)^0,3; Schatten = NFW mit kosm. Anteil (wie Notiz 17); a0 = 1,2e-10")
P("  Spalten: g_Op/g_RAR (Gedaechtnis, stationaer, q = exp(-sqrt y) lokal) | g_Hybrid/g_RAR | g(q=0,8429 fest)/g_RAR")
for Mb, M200, cc in ((6e10, 1e12, 10), (1e9, 1e11, 12)):
    ab = 1.5 * kpc * (Mb / 1e10) ** 0.3
    P(f"   M_b = {Mb:.0e}, M200 = {M200:.0e}:")
    for rk in (5, 10, 20, 50, 100):
        r = rk * kpc
        mb = hern_M(r, Mb * Msun, ab)
        mh = nfw_M(r, M200 * Msun * (1 - Ob / (Ob + Oc)), cc)
        gb = G * mb / r**2
        gt = G * (mb + mh) / r**2
        gR = nu_ftm(gb / a0) * gb
        # stationaeres Gedaechtnis: Fixpunkt g = g_N + q g  (Iteration als Kontrolle)
        qloc = math.exp(-math.sqrt(gb / a0))
        g = gb
        for _ in range(2000):
            g = gb + qloc * g
        gH = nu_ftm(gt / a0) * gt
        gfix = gb / (1 - q536)
        P(f"     r = {rk:4d} kpc  y = {gb/a0:8.4f}   {g/gR:7.4f}   {gH/gR:7.3f}   {gfix/gR:7.3f}")

# ==========================================================================================
hdr("4  Wachstum nach z* (Baryonen druckfrei, Strahlung glatt, subhorizontal)  ' = d/dln a")
P("  d_b'' + (2 + dlnH/dlna) d_b' = 3/2 [Ob(a)(d_b + M) + Oc(a) d_c];  M' = beta[q(d_b+M) - M]")
P("  Wachstumsexponent im Materie-Zeitalter, D ~ a^p:  p(p+1/2) = 3/2 (1-q)[1 + q beta/(beta(1-q)+p)]")
for b in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 1e6):
    lo, hi = 0.0, 2.0
    for _ in range(80):
        p = 0.5 * (lo + hi)
        f = p * (p + 0.5) - 1.5 * (1 - q536) * (1 + q536 * b / (b * (1 - q536) + p))
        lo, hi = (p, hi) if f < 0 else (lo, p)
    P(f"   beta = {b:8g}:  p = {p:.4f}")


def grow(r, kval, mode):
    """Lineares Wachstum von a* bis 1 fuer ein k (Index), Anfangswerte aus dem Spielzeug-Boltzmann."""
    y = r["y"]
    bg = r["bg"]
    cfg = r["cfg"]
    i = int(np.argmin(abs(kgrid - kval)))
    Hc = bg.Hc(a_s)
    db = 0.75 * y[i, 1]
    ub = -y[i, 2] / Hc
    dc = y[i, 5]
    uc = -y[i, 6] / Hc
    M = y[i, 7]
    wS = cfg.get("wS", 0.0)
    enh = cfg.get("enhQ", 0.0)
    beta, qm = cfg.get("beta", 0.0), cfg.get("q", 0.0)
    mem = cfg.get("mem", False)
    Ns = np.linspace(math.log(a_s), 0.0, 20001)
    st = np.array([db, ub, dc, uc, M])

    def rhs(N, s):
        a = math.exp(N)
        fr, fb, fm, fL = bg.frac(a)
        src = 1.5 * (fb * ((1 + enh) * s[0] + (s[4] if mem else 0.0)) + wS * fb * s[2])
        # d/dlna von d mit (1/a^2 d/dt a^2 ...) : d'' + (2 + dlnH/dlna) d' = src, dlnH/dlna = dlnE/dlna
        fric = 2 + bg.dlnE(a)
        dM = beta * (qm * (s[0] + s[4]) - s[4]) if mem else 0.0
        return np.array([s[1], -fric * s[1] + src, s[3], -fric * s[3] + src, dM])

    hist = {}
    for j in range(len(Ns) - 1):
        hN = Ns[j + 1] - Ns[j]
        N = Ns[j]
        k1 = rhs(N, st); k2 = rhs(N + hN / 2, st + hN / 2 * k1)
        k3 = rhs(N + hN / 2, st + hN / 2 * k2); k4 = rhs(N + hN, st + hN * k3)
        st = st + hN / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        for zz in (1.0, 0.5, 0.0):
            if zz not in hist and Ns[j + 1] >= math.log(1 / (1 + zz)) - 1e-12:
                hist[zz] = st.copy()
    return (db, dc, M), hist


# Normierung: phi_p = (2/3) zeta, Delta_zeta^2 = 2,1e-9 (Planck, skaleninvariant genaehert)
phip = (2 / 3) * math.sqrt(2.1e-9)
P()
P(f"  Normierung: Delta_phi_p = (2/3) sqrt(2,1e-9) = {phip:.3e} pro ln k (n_s = 1). Tabellen: Delta = |delta|*Delta_phi_p.")
P("  Variante                         k[1/Mpc]  Delta_b(z*)  Delta_well(z*)  Delta_b(z=0)  Delta_eff(z=0)  f_b(z=0,5)  Delta_b(0)/LCDM")
cases = [("(a) LCDM", "a"), ("(b) nur Baryonen", "b"), ("(c) mitschwingend", "c"),
         ("(i) beta=0,01", "i0.01"), ("(i) beta=0,1", "i0.1"), ("(i) beta=1", "i1.0"), ("(i) beta=30", "i30.0"),
         ("(ii) w=1", "ii1.0")]
refD = {}
growD = {}
for kval in (0.05, 0.1, 0.135):
    for nm, key in cases:
        r = results[key]
        (db, dc, M), hist = grow(r, kval, key)
        cfg = r["cfg"]
        s0, s05 = hist[0.0], hist[0.5]
        wS = cfg.get("wS", 0.0)
        mem = cfg.get("mem", False)
        well = wS * dc + (M if mem else 0.0) + cfg.get("enhQ", 0.0) * db
        Om_eff = r["bg"].Om
        eff0 = Ob * ((1 + cfg.get("enhQ", 0.0)) * s0[0] + (s0[4] if mem else 0) + wS * s0[2]) / Om_eff
        fb05 = s05[1] / s05[0]
        if key == "a":
            refD[kval] = s0[0]
        growD[(key, kval)] = s0[0]
        P(f"  {nm:32s} {kval:6.3f}   {abs(db)*phip:10.2e}   {abs(well)*phip:10.2e}    {abs(s0[0])*phip:10.3e}    {abs(eff0)*phip:10.3e}    {fb05:7.3f}    {s0[0]/refD[kval]:8.3f}")
    P()
P("  (k = 0,135/Mpc = 0,2 h/Mpc entspricht grob der sigma_8-Skala R = 8 Mpc/h.)")
P()
P("  Aufholen durch lokale Verstaerkung (Tief-a0-Regime): Operator mit lokalem nu auf die Baryonenquelle,")
P("  p(p+1/2) = 3/2 (1-q) nu_lok; noetige e-Folds n = ln(F)/(p-1) fuer den Rueckstand F von (i) beta = 1 bei k = 0,135:")
Fdef = abs(refD[0.135] / growD[("i1.0", 0.135)])
for nul in (3.3, 28.2, 62.3, 106.5):
    pp = (-0.5 + math.sqrt(0.25 + 4 * 1.5 * (1 - q536) * nul)) / 2
    P(f"   nu_lok = {nul:6.1f}:  p = {pp:.3f},  n = {math.log(Fdef)/(pp-1) if pp > 1 else float('inf'):.2f} e-Folds (F = {Fdef:.1f})")
P("  (nu_lok-Werte = 1+Q aus Abschnitt 2d(3); nu haengt von delta ab, Superposition verletzt, vgl. Notiz 08, Befund 2.)")
P()
H0inv_Gyr = 9.7779 / h
for tmem in (0.170, 0.260):
    P(f"  Bullet-Schranke 1/lam >= {tmem*1000:.0f} Myr heute (Notiz 12), lam = beta H0:  beta <= {H0inv_Gyr/tmem:.0f}  (H0^-1 = {H0inv_Gyr:.2f} Gyr)")

# ------------------------------------------------------------------
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "freiheitsgrad.txt")
with open(path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n".join(OUT))
