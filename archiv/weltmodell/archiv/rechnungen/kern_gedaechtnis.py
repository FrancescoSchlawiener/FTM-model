# -*- coding: utf-8 -*-
"""
Kopfgleichung von 07 auf exzentrischen Bahnen: Ist die örtliche Regel a = ν(|g_N|)g_N der Grenzfall
x = ω/λ → ∞, und wie klein muss die Zerfallsrate λ_K des Kerns sein, damit alte Scheiben ihre Radialstreuung behalten?

Gleichung (07):  r̈ = −∇Φ_bar − χ κ₀ ∫_{−∞}^{t} Δ (e_r⊗e_r) e^{−λ(t−t′)} ṙ(t′) dt′
Lesart: der Projektor e_r⊗e_r steht zur Zeit t (Kraft längs e_r(t)); nur so gibt die Kreisbahn
a_K = κ₀Δ·rω²/(λ²+ω²) nach innen (e_r(t)·e_φ(t′) = sin ω(t−t′)). χ, Δ = q_Δs = Δs/r und κ₀ am jetzigen Ort
(Variante mit χΔκ₀ unter dem Integral am Ort t′ unten, Teil 2c).

Hilfsvariable (Herleitung):
  Y(t) = ∫_{−∞}^{t} e^{−λ(t−t′)} ṙ(t′) dt′   ⇒   Ẏ = ṙ − λY          (Ableitung nach der oberen Grenze + Zerfall)
  partielle Integration:  Y = r(t) − B(t),   B(t) = λ∫_0^∞ e^{−λs} r(t−s) ds   (exponentiell gewichtete Vergangenheitslage)
  ⇒  Ḃ = λ(r − B),   Kernkraft  F_K = −K(r)·(e_r·(r − B))·e_r = −K(r)·(r − e_r·B)·e_r,   K = χ(s)·(Δs/r)·κ₀(r)
Der Kern ist eine Feder mit Stärke K zur gleitend gemittelten eigenen Vergangenheitslage.
Die Kraft ist zentral (längs e_r): der Drehimpuls bleibt exakt erhalten, Energie geht nur über die Radialbewegung.

Normierung (07, Auswertung): κ₀/λ_eff = 3c²/(r*·v),  λ_eff = ω + λ²/ω, auf der Kreisbahn der Regel (v, ω am Ort r):
  K(r) = χ(s)·(Δs/r)·κ₀(r) = χ(s)·(GM/(r*r²))·(1 + λ²/ω_c(r)²)   mit Δs·3c² = GM.
Grenzfall λ → 0: B → 0 (Mittel über viele gedrehte Rosettenblätter), F_K → −χ·GM/(r*r) = −χ√(a_N a₀): genau die örtliche Regel.
Einheiten der Integration: GM = 1, r* = 1 (also a₀ = 1), Zeit in √(r*³/GM).
"""
import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from ftm_kern import *

out = []
def P(s=""):
    out.append(s); print(s)

M = 6e10 * M_sun
RS = r_stern(M, A0)
TS = math.sqrt(RS**3 / (G * M))               # Zeiteinheit
P(f"M = 6·10¹⁰ M☉, a₀ = {A0:.2e} m/s², r* = {RS/kpc:.3f} kpc, Δs = {ds(M):.4e} m, Zeiteinheit √(r*³/GM) = {TS/Gyr*1e3:.2f} Mio. Jahre")

# ---------------------------------------------------------------- dimensionslose Regel (GM = r* = a₀ = 1)
def aN_(r):
    return 1.0 / r**2

def a_rule(r):                                 # Betrag der örtlichen Regel
    return float(a_regel(aN_(r), 1.0))

def om_c(r):                                    # Kreisfrequenz der Regel
    return math.sqrt(a_rule(r) / r)

def K_(r, lam):
    return float(chi(aN_(r), 1.0)) / r**2 * (1 + lam**2 / om_c(r)**2)

# ---------------------------------------------------------------- 1 Normierung und Kreisbahn
P()
P("1) Normierung κ₀ und Kreisbahnprobe (Kernintegral numerisch mit quad, nicht die Formel)")
P("   κ₀(r) = 3c²/(r*·v)·(ω + λ²/ω) mit v, ω der Regel am Ort r, λ = 0,01·ω(r*)")
Om_star = v_kreis_punkt(M, RS) / RS
lamSI = 0.01 * Om_star

def kern_int(kap0, lam, r, w):
    """∫_0^∞ κ₀ e^{−λs} r ω sin(ωs) ds numerisch (Periodenstücke)."""
    T = 2 * math.pi / w
    tot, n = 0.0, 0
    while True:
        piece = quad(lambda s: math.exp(-lam * s) * math.sin(w * s), n * T, (n + 1) * T, limit=200)[0]
        tot += piece; n += 1
        if math.exp(-lam * n * T) < 1e-16:
            break
    return kap0 * r * w * tot

def v_kern(r, kap0, lam):
    aN = G * M / r**2
    D = ds(M) / r
    ch = float(chi(aN, A0))
    f = lambda v: v * v / r - aN - ch * D * kern_int(kap0, lam, r, v / r)
    return brentq(f, math.sqrt(aN * r), 3 * math.sqrt(aN * r) + 3e5)

def kap0_lokal(r, lam):
    v = v_kreis_punkt(M, r); w = v / r
    return 3 * c**2 / (RS * v) * (w + lam**2 / w)

k0_star = kap0_lokal(RS, lamSI)
P("   r/r*   κ₀ lokal [s⁻²]   v Regel [km/s]   v Kern, κ₀ lokal   Abw.        v Kern, κ₀ fest (Wert bei r*)   Abw.")
for x in [0.5, 1, 2, 5]:
    r = x * RS
    vr = v_kreis_punkt(M, r)
    k0 = kap0_lokal(r, lamSI)
    v1 = v_kern(r, k0, lamSI)
    v2 = v_kern(r, k0_star, lamSI)
    P(f"   {x:4.1f}   {k0:.4e}      {vr/1e3:8.3f}        {v1/1e3:8.3f}        {v1/vr-1:+.1e}     {v2/1e3:8.3f}                        {v2/vr-1:+.3f}")
P("   κ₀ ∝ (ω + λ²/ω)/v ≈ 1/r: mit festem κ₀ trifft die Kreisbahn nur am Normierungsradius; κ₀ muss am Ort gelesen werden.")

# ---------------------------------------------------------------- 2 Bahnen
def rhs_regel(t, u):
    x, y, vx, vy = u
    r = math.hypot(x, y)
    a = a_rule(r)
    return [vx, vy, -a * x / r, -a * y / r]

def make_rhs_kern(lam, variante="t"):
    def rhs(t, u):
        x, y, vx, vy, Bx, By = u
        r = math.hypot(x, y); ex, ey = x / r, y / r
        aN = aN_(r)
        if variante == "t":                    # B = gewichtete Vergangenheitslage, Ḃ = λ(r − B)
            proj = r - (ex * Bx + ey * By)
            aK = K_(r, lam) * proj
            dB = [lam * (x - Bx), lam * (y - By)]
        else:                                  # Variante: K am Ort t′ unter dem Integral, Ẏ = K(r)ṙ − λY
            aK = ex * Bx + ey * By
            k = K_(r, lam)
            dB = [k * vx - lam * Bx, k * vy - lam * By]
        a = aN + aK
        return [vx, vy, -a * ex, -a * ey, dB[0], dB[1]]
    return rhs

def Phi_regel(r):
    return quad(a_rule, 1.0, r)[0]

def start(rp, e):
    ra = rp * (1 + e) / (1 - e)
    L2 = 2 * (Phi_regel(ra) - Phi_regel(rp)) / (1 / rp**2 - 1 / ra**2)
    return ra, math.sqrt(L2)

def wendepunkte(sol):
    """Peri/Apo aus Vorzeichenwechseln von v_r (dichte Ausgabe, Brent)."""
    t = sol.t
    x, y, vx, vy = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    vr = (x * vx + y * vy) / np.hypot(x, y)
    ev = []
    for i in range(len(t) - 1):
        if vr[i] == 0 or vr[i] * vr[i + 1] < 0:
            f = lambda tt: (lambda u: (u[0] * u[2] + u[1] * u[3]) / math.hypot(u[0], u[1]))(sol.sol(tt))
            try:
                tt = brentq(f, t[i], t[i + 1])
            except ValueError:
                continue
            u = sol.sol(tt)
            ev.append((tt, math.hypot(u[0], u[1]), math.atan2(u[1], u[0]), "apo" if vr[i] > 0 else "peri"))
    return ev

def auswerten(sol):
    ev = wendepunkte(sol)
    peri = [e for e in ev if e[3] == "peri"]
    apo = [e for e in ev if e[3] == "apo"]
    n = min(len(peri), len(apo))
    hub = np.array([apo[i][1] - peri[i][1] for i in range(n)])
    th = np.array([0.5 * (apo[i][0] + peri[i][0]) for i in range(n)])
    tp = np.array([p[0] for p in peri])
    ph = np.unwrap(np.array([p[2] for p in peri]))
    Tr = np.mean(np.diff(tp[-10:]))                       # Radialperiode aus den letzten 10 Perizentren
    dphi = np.mean(np.diff(ph[-10:])) % (2 * math.pi)
    dphi = dphi if dphi > 0.5 else dphi + 2 * math.pi
    kap = 2 * math.pi / Tr
    Om = dphi / Tr                                         # mittlere Winkelgeschwindigkeit
    return hub, th, kap, Om

def B_vorgeschichte(lam, u0, Tr_guess):
    """B(0) = λ∫_0^∞ e^{−λs} r(−s) ds unter der Regel; Bahn spiegelsymmetrisch zum Perizentrum auf der x-Achse:
       z(−s) = conj z(s); z(s + T_r) = z(s)·e^{iΔφ}  ⇒  ∫_0^∞ = ∫_0^{T_r}(…)/(1 − e^{−λT_r + iΔφ})."""
    sol = solve_ivp(rhs_regel, [0, 3 * Tr_guess], u0, method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
    ev = [e for e in wendepunkte(sol) if e[3] == "peri" and e[0] > 1e-6]
    Tr, dphi = ev[0][0], ev[0][2]
    I = quad(lambda s: math.exp(-lam * s) * sol.sol(s)[0], 0, Tr, limit=500)[0] + \
        1j * quad(lambda s: math.exp(-lam * s) * sol.sol(s)[1], 0, Tr, limit=500)[0]
    Z = lam * I / (1 - cmath_exp(-lam * Tr + 1j * dphi))
    Z = Z.conjugate()
    return Z.real, Z.imag

def cmath_exp(z):
    import cmath
    return cmath.exp(z)

NRAD = 50
LAMS = [0.1, 0.01, 0.001, 1e-4]
ergebnis = {}
P()
P("2) Exzentrische Bahnen, Start im Perizentrum r_p = 1 r*, e = (Apo − Peri)/(Apo + Peri); 50 Radialperioden")
P("   Kern: B(0) aus der Vorgeschichte unter der Regel (exakt über eine Radialperiode aufsummiert)")
for e in [0.1, 0.3]:
    rp = 1.0
    ra, L = start(rp, e)
    rg = brentq(lambda r: r**3 * a_rule(r) - L**2, rp, ra)
    Omg = om_c(rg)
    kap_ep = Omg * math.sqrt(realteil_quadrat_punkt(1 / rg))
    u0 = [rp, 0.0, 0.0, L / rp]
    Tr0 = 2 * math.pi / kap_ep
    T = NRAD * Tr0
    solR = solve_ivp(rhs_regel, [0, T], u0, method="DOP853", rtol=1e-11, atol=1e-12, dense_output=True, max_step=Tr0 / 40)
    hubR, thR, kapR, OmR = auswerten(solR)
    P()
    P(f"   e = {e}: Apo = {ra:.4f} r*, Führungsradius r_g = {rg:.4f} r*, Ω_g = {Omg:.4f}, √(1+χ(r*/r_g)) = {kap_ep/Omg:.5f}")
    P(f"   Regel:  Hub {hubR[0]:.5f} → {hubR[-1]:.5f} r*,  κ/Ω = {kapR/OmR:.5f},  Ω/Ω_g = {OmR/Omg:.5f},  κ/Ω_g = {kapR/Omg:.5f}")
    P("   λ_K/Ω_g   Hub Anfang   Hub nach 50 T_r   Abnahme je T_r   Rate Γ_Hub/λ   κ/Ω       Ω/Ω_g     κ/Ω_g     |L−L₀|/L₀")
    for f in LAMS:
        lam = f * Omg
        Bx, By = B_vorgeschichte(lam, u0, Tr0)
        sol = solve_ivp(make_rhs_kern(lam), [0, T], u0 + [Bx, By], method="DOP853", rtol=1e-11, atol=1e-12,
                        dense_output=True, max_step=Tr0 / 40)
        hub, th, kap, Om = auswerten(sol)
        # Rate aus linearer Anpassung von ln Hub gegen t (ohne die ersten 2 Perioden)
        sel = th > 2 * Tr0
        selR = thR > 2 * Tr0
        rate = -np.polyfit(th[sel], np.log(hub[sel]), 1)[0] + np.polyfit(thR[selR], np.log(hubR[selR]), 1)[0]
        per = 1 - math.exp(-rate * 2 * math.pi / kap)
        Lend = sol.y[0, -1] * sol.y[3, -1] - sol.y[1, -1] * sol.y[2, -1]
        P(f"   {f:7.0e}   {hub[0]:.5f}      {hub[-1]:.5f}          {per:.3e}        {rate/lam:7.4f}     {kap/Om:.5f}   {Om/Omg:.5f}   {kap/Omg:.5f}   {abs(Lend-L)/L:.1e}")
        ergebnis[(e, f)] = dict(rate=rate / lam, kap=kap / Om, rg=rg, Omg=Omg)
        i20 = np.argmin(abs(th - 20 * 2 * math.pi / kapR))
        P(f"             Hub nach 20 T_r / Anfang = {hub[i20]/hub[0]:.3f}")
    if e == 0.3:
        lam = 0.01 * Omg
        sol = solve_ivp(make_rhs_kern(lam, "t'"), [0, T], u0 + [0.0, 0.0],   # Y(0) = 0: Einschwingen, daher Rate erst ab 20 T_r
                        method="DOP853", rtol=1e-11, atol=1e-12, dense_output=True, max_step=Tr0 / 40)
        hub, th, kap, Om = auswerten(sol)
        sel = th > 20 * Tr0
        rate = -np.polyfit(th[sel], np.log(hub[sel]), 1)[0]
        P(f"   2c) Variante χΔκ₀ am Ort t′ (Ẏ = K(r)ṙ − λY), λ = 0,01 Ω_g, Y(0) = 0, Rate aus T_r 20–50: Γ_Hub/λ = {rate/lam:.4f}, κ/Ω = {kap/Om:.5f}")

# ---------------------------------------------------------------- 3 Linear: Eigenwerte und Handformel
P()
P("3) Unabhängige Prüfung ohne Simulation")
P("   a) Handformel (Störungsrechnung um die Kreisbahn r_g, Epizykel x = X cos κt, r_gφ₁ = −(2Ω/κ)X sin κt;")
P("      der Projektor dreht mit: im Inertialsystem schwingt die Bahn mit Ω ± κ, nicht mit κ):")
P("      Γ_X = −⟨P⟩/(κ²X²),  ⟨P⟩/X² = K[−λΩ²/(λ²+Ω²) + ½κλ(−A(Ω+κ)/(λ²+(Ω+κ)²) + B(Ω−κ)/(λ²+(Ω−κ)²))],")
P("      A = ½ − Ω/κ, B = ½ + Ω/κ;  1D-Maxwell (nur Frequenz κ, Verlustanteil κλ/(λ²+κ²)): Γ_X = Kλ/(2(λ²+κ²))")
P("   b) Eigenwerte des linearisierten Systems im mitrotierenden Rahmen (r, v_r, φ₁, B′), B′ = e^{−iΩt}B")

def hand(rg, lam):
    Om = om_c(rg); kap = Om * math.sqrt(realteil_quadrat_punkt(1 / rg)); K = K_(rg, lam)
    A = 0.5 - Om / kap; Bc = 0.5 + Om / kap; w1 = Om + kap; w2 = Om - kap
    Pp = K * (-lam * Om**2 / (lam**2 + Om**2) + 0.5 * kap * lam * (-A * w1 / (lam**2 + w1**2) + Bc * w2 / (lam**2 + w2**2)))
    return -Pp / kap**2, K * lam / (2 * (lam**2 + kap**2)), kap

def eigen(rg, lam):
    Om = om_c(rg); L = rg**2 * Om
    def F(u):
        r, vr, p1, bx, by = u
        Bc = complex(bx, by)
        proj = r - (np.exp(-1j * p1) * Bc).real
        dv = L**2 / r**3 - aN_(r) - K_(r, lam) * proj
        dB = lam * (r * np.exp(1j * p1) - Bc) - 1j * Om * Bc
        return np.array([vr, dv, L / r**2 - Om, dB.real, dB.imag])
    B0 = lam * rg / (lam + 1j * Om)
    u0 = np.array([rg, 0, 0, B0.real, B0.imag])
    assert np.max(np.abs(F(u0))) < 1e-10
    J = np.zeros((5, 5)); h = 1e-7
    for j in range(5):
        du = np.zeros(5); du[j] = h
        J[:, j] = (F(u0 + du) - F(u0 - du)) / (2 * h)
    return np.linalg.eigvals(J) / Om, Om

P("   r_g/r*   λ/Ω     Γ_X/λ Hand   Γ_X/λ Eigenwert   Γ_X/λ Simulation   Γ_X/λ 1D-Maxwell   Eigenfrequenzen/Ω (Im, Re)")
for e in [0.1, 0.3]:
    for f in LAMS:
        rg = ergebnis[(e, f)]["rg"]
        lam = f * om_c(rg)
        gh, g1, kap = hand(rg, lam)
        ev, Om = eigen(rg, lam)
        osc = [z for z in ev if abs(z.imag) > 0.5]
        # Epizykelmode: Imaginärteil nahe κ/Ω
        epi = min(osc, key=lambda z: abs(abs(z.imag) - kap / Om))
        ge = -epi.real * Om / lam
        fr = "  ".join(f"{abs(z.imag):.4f}/{z.real:+.2e}" for z in sorted(ev, key=lambda z: abs(z.imag)) if z.imag >= 0)
        P(f"   {rg:.4f}  {f:6.0e}   {gh/lam:9.4f}     {ge:9.4f}          {ergebnis[(e, f)]['rate']:9.4f}          {g1/lam:9.4f}        {fr}")

for rg in [0.917]:
    for f in [1e-3, 1e-4]:
        lam = f * om_c(rg); gh, g1, kap = hand(rg, lam); ev, Om = eigen(rg, lam)
        P(f"   Sonnenkreis r_g = {rg} r*, λ/Ω = {f:.0e}: Γ_X/λ Hand = {gh/lam:.4f}, κ/Ω = {kap/Om:.4f}")
P("   Eigenfrequenzen: 0 = Drehung (Symmetrie); ≈ 1 Ω im mitrotierenden Rahmen = ruhender Versatz des Bahnmittelpunkts B im Inertialsystem;")
P("   ≈ κ/Ω = Epizykel. Re = Anwachsen (+) oder Abklingen (−) in Einheiten Ω.")

# ---------------------------------------------------------------- 4 Hochrechnung Sonnenkreis
P()
P("4) Hochrechnung auf 10 Mrd. Jahre am Sonnenkreis")
Om_sun = 220e3 / (8 * kpc)
t10 = 10 * Gyr
P(f"   Ω = 220 km/s / 8 kpc = {Om_sun:.4e} s⁻¹,  Ω·10 Mrd. Jahre = {Om_sun*t10:.1f} rad ({Om_sun*t10/(2*math.pi):.1f} Umläufe)")
P(f"   Sonnenkreis im Modell: 8 kpc = {8*kpc/RS:.3f} r* (Punktmasse 6·10¹⁰ M☉), v Regel = {v_kreis_punkt(M, 8*kpc)/1e3:.1f} km/s")
P("   Hub(t) = Hub₀·e^{−Γ t}, Γ = g·λ_K;  weniger als 10 % Abnahme: λ_K/Ω < −ln 0,9/(g·Ω t)")
gmax = 0
for e in [0.1, 0.3]:
    g = ergebnis[(e, 1e-4)]["rate"]
    gmax = max(gmax, g)
    P(f"   e = {e}: g = Γ/λ = {g:.3f} (λ_K = 10⁻⁴ Ω)  →  λ_K/Ω < {-math.log(0.9)/(g*Om_sun*t10):.2e},  1/λ_K > {g*Om_sun*t10/(-math.log(0.9))/Om_sun/Gyr:.0f} Mrd. Jahre")
lam_max = -math.log(0.9) / (gmax * Om_sun * t10)
P("   λ_K/Ω = 0,1: Epizykel nach ~10 T_r abgeklungen; übrig bleibt der Versatz des Bahnmittelpunkts (κ = Ω), Abklingrate 0,0057 Ω (e = 0,1)")
for f in LAMS[1:]:
    P(f"   λ_K/Ω = {f:6.0e}: Hub nach 10 Mrd. Jahren (g = {gmax:.2f}) = {math.exp(-gmax*f*Om_sun*t10):.3e} des Anfangs")

# ---------------------------------------------------------------- 5 Übergaberate der Mode ξ
P()
P("5) Übergaberate λ_m = ω_orb/N*, N* = 0,3–2,4")
for D in [200, 500]:
    w_, lo, hi = uebergabezeit(D)
    P(f"   gebundener Bereich mit mittlerer Dichte {D}·ρ_krit: ω_orb = {1/w_:.3e} s⁻¹, λ_m/Ω_Sonne = {1/(hi*Om_sun):.4f} – {1/(lo*Om_sun):.4f}")
w_loc = math.sqrt(G * M / (8 * kpc)**3)     # (4πG/3)ρ̄ mit ρ̄ = baryonische Masse in 8 kpc (Punktmasse)
P(f"   Bereich innerhalb des Sonnenkreises (Baryonen 6·10¹⁰ M☉ in 8 kpc): ω_orb = {w_loc:.3e} s⁻¹ = {w_loc/Om_sun:.3f} Ω, λ_m/Ω = {w_loc/Om_sun/2.4:.2f} – {w_loc/Om_sun/0.3:.2f}")
P(f"   Grenze aus 4): λ_K/Ω < {lam_max:.1e};  Verhältnis kleinstes λ_m / größtes λ_K = {1/(uebergabezeit(200)[2]*Om_sun)/lam_max:.0f}")

with open("kern_gedaechtnis.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
