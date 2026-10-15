# -*- coding: utf-8 -*-
"""
Torus-Wechselkräfte  (Forschungsnotiz 27)
=========================================

Baut auf 03-wechselwirkungen, 05-zustandsraum, 06-torus-dynamik auf und rechnet
nur Neues:

  1  Schließung eliminiert Fluchtmoden
     1a Spektren: Intervall (Rand, absorbierend) / Kreis (identifiziert) / offene Gerade
     1b Identifikation als Projektor P = Mittel über die Decktranslationen,
        [P, Δ] = 0, Spektrum von PΔ = Kreisspektrum, Π⊥ = I - P
     1c Diffusion: Verbleib im Fenster, Rückkehrdichte, Aufenthaltsanteil
     1d Liouville (ballistisch): Verbleib, Gleichverteilung; Kac-Rückkehrzeiten
     1e Wellenpaket: exakte Wiederkehr auf dem Ring, Zerfließen auf der Geraden
     1f Pólya: Rückkehr Z^3 (transient) vs Torus Z_N^3 (Kac: N^3)
     1g Bindung als Funktion der Schließungslänge L; Toy Model mit K = 1
  2  Die drei Rollen als Kräfte am Torus (zwei Stuart-Landau-Zyklen, Kopplung K e^{iφ})
     Radius-, Dicken- und Phasengleichgewicht; dünn / eingerastet / zerfallen;
     Temperaturpfad; Zahlenbeispiel
  3  Dreieck Energie - Struktur - Zustand als 3x3-Kopplungssystem
     3a Hurwitz über Schleifenverstärkungen (Zufallstest)
     3b Sekantenbedingung für den negativen 3-Zyklus
     3c nichtlineares E-S-Z-Modell: Attraktor, Alterungsdrift, Kippen, kritisches Verlangsamen
     3d Gedächtniskern: Toy-Model-Kern = Maxwell-Element
  4  Wärme als langsames Gedächtnis: eingespannter Kupferdraht (magnetomotorisch),
     Schleifenverstärkung zerlegt in σ(T), E(T), Wärmespannung; Ferrit am Curie-Punkt
  5  Gravitation: Vorspannung -> Steifigkeit -> Resonanz; χ_g, Existenzgrenzen,
     Verhältnis differentieller Kanal / Takt-Kanal in FTM-Größen
  6  Neue Aussagen: Kac-Budget, Koexistenz (Lamb), Höchstzahl freier Zyklen (Kuramoto)

Nur numpy. Ausgabe zusätzlich nach torus_wechselkraefte.txt (gleicher Ordner).
"""
import io
import math
import os
import time

import numpy as np

_BUF = io.StringIO()
t_start = time.time()


def out(*args):
    text = " ".join(str(a) for a in args)
    print(text)
    _BUF.write(text + "\n")


def hr(title):
    out("")
    out("=" * 78)
    out(title)
    out("=" * 78)


def sub(title):
    out("")
    out("--- " + title)


erf = np.vectorize(math.erf)
PI = math.pi

# ===========================================================================
hr("TEIL 1  Schließung eliminiert Fluchtmoden")
# ===========================================================================
D = 1.0
L = 1.0

sub("1a Spektrum von -D∂² (diskret, n = 200 Punkte pro Länge L)")
n = 200
h = L / n
kk = np.arange(1, n + 1)
lam_dir = (4 / h ** 2) * np.sin(kk * PI / (2 * (n + 1))) ** 2 * D   # Intervall [0,L], Rand absorbierend
kp = np.arange(0, n)
lam_per = np.sort((4 / h ** 2) * np.sin(PI * kp / n) ** 2 * D)       # Kreis der Länge L
out(f"Intervall (Dirichlet):  λ1 = {lam_dir[0]:.4f}   (π²D/L² = {PI**2*D/L**2:.4f}),  kein Nullmode")
out(f"Kreis (periodisch):     λ0 = {lam_per[0]:.2e}, λ1 = λ2 = {lam_per[1]:.4f}   (4π²D/L² = {4*PI**2*D/L**2:.4f})")
gap = 4 * PI ** 2 * D / L ** 2
for Lam in [10.0, 100.0, 1000.0]:
    nn = int(n * Lam)
    hh = Lam / nn
    kk2 = np.arange(1, nn + 1)
    lamO = (4 / hh ** 2) * np.sin(kk2 * PI / (2 * (nn + 1))) ** 2 * D
    cnt = int(np.sum(lamO < gap))
    out(f"offener Bereich Λ = {Lam:6.0f} L:  λ1 = {lamO[0]:.3e},  Zahl der Moden unter der Kreislücke 4π²D/L²: {cnt:5d}"
        f"  (Weyl Λ·2/L = {int(2*Lam)})")
out("Kreis: genau 1 Mode (die Konstante) unter der Lücke; offen: Zahl ~ Λ -> Kontinuum ohne Lücke.")

sub("1b Identifikation als Projektor  P = (1/M) Σ_m T_{mL}  auf dem Überlagerungskreis der Länge M·L")
n1, M = 32, 8
N1 = n1 * M
hh = L / n1
Lap = (2 * np.eye(N1) - np.roll(np.eye(N1), 1, axis=0) - np.roll(np.eye(N1), -1, axis=0)) / hh ** 2
S = np.roll(np.eye(N1), n1, axis=0)
Pm = np.zeros((N1, N1))
Sk = np.eye(N1)
for m in range(M):
    Pm += Sk
    Sk = S @ Sk
Pm /= M
comm = np.linalg.norm(Pm @ Lap - Lap @ Pm)
idem = np.linalg.norm(Pm @ Pm - Pm)
w, V = np.linalg.eigh(Pm)
B = V[:, w > 0.5]
Bp = V[:, w < 0.5]
lam_r = np.sort(np.linalg.eigvalsh(B.T @ Lap @ B))
lap_c = (2 * np.eye(n1) - np.roll(np.eye(n1), 1, axis=0) - np.roll(np.eye(n1), -1, axis=0)) / hh ** 2
lam_c = np.sort(np.linalg.eigvalsh(lap_c))
lam_p = np.sort(np.linalg.eigvalsh(Bp.T @ Lap @ Bp))
out(f"‖P² - P‖ = {idem:.1e},  ‖[P, Δ]‖ = {comm:.1e}   (Identifikation vertauscht mit dem Generator)")
out(f"Rang P = {B.shape[1]} (= Punkte auf dem Kreis L),  Rang Π⊥ = I - P = {Bp.shape[1]}")
out(f"max |Spektrum(PΔP auf Bild P) - Spektrum(Δ auf Kreis L)| = {np.max(np.abs(lam_r - lam_c)):.1e}")
out(f"Moden unter der halben Kreislücke 2π²/L² = {gap/2:.3f}:  im Bild von P: {int(np.sum(lam_r < gap/2))} (die Konstante),"
    f"  im Bild von Π⊥: {int(np.sum(lam_p < gap/2))}  (kleinster Π⊥-Eigenwert {lam_p[0]:.4f} ≈ (2π/(M L))² = {(2*PI/(M*L))**2:.4f})")
out(f"   kleinster Π⊥-Eigenwert analytisch (4/h²)sin²(π/N) = {(4/hh**2)*np.sin(PI/N1)**2:.4f}")
out("=> Π⊥ enthält genau die langwelligen Fluchtmoden (k ∉ 2πZ/L); P löscht sie, ohne einen Term zu Δ hinzuzufügen.")

sub("1c Diffusion (D = 1, L = 1, Start in der Fenstermitte)")
out(" t [L²/D] | Gerade: im Fenster L | Intervall absorbierend: Überleben | Kreis: im Raum | Rückkehrdichte Gerade | Kreis")
for t in [0.01, 0.1, 1.0, 10.0, 100.0]:
    p_line = math.erf((L / 2) / math.sqrt(4 * D * t))
    ks = np.arange(1, 400, 2)
    surv = np.sum(4 / (ks * PI) * (-1.0) ** ((ks - 1) // 2) * np.exp(-D * ks ** 2 * PI ** 2 * t / L ** 2))
    rho_line = 1 / math.sqrt(4 * PI * D * t)
    nn_ = np.arange(1, 200)
    rho_circ = (1 / L) * (1 + 2 * np.sum(np.exp(-4 * PI ** 2 * nn_ ** 2 * D * t / L ** 2)))
    out(f" {t:7.2f}  | {p_line:.6f}             | {surv:.6e}                    | 1            | {rho_line:.5f}              | {rho_circ:.5f}")
out(f"Abklingrate des Überlebens im Intervall: π²D/L² = {PI**2:.4f}; Rückkehrdichte Kreis -> 1/L, Gerade ~ t^(-1/2) -> 0")

a_w = 0.1 * L
out(f"Aufenthaltsanteil (Zeitmittel über [0,t]) im Fenster der Breite a = {a_w} um den Start:")
nn_ = np.arange(1, 20001)
for t in [10.0, 1e3, 1e5]:
    kn2 = (2 * PI * nn_ / L) ** 2
    frac_c = a_w / L + np.sum(2 / (PI * nn_) * np.sin(PI * nn_ * a_w / L) * (1 - np.exp(-kn2 * D * t)) / (kn2 * D * t))
    # Gerade: (1/t) ∫0^t erf(a/(4 sqrt(D s))) ds, Substitution s = t u²
    u = np.linspace(1e-6, 1, 200001)
    integrand = erf(a_w / (4 * np.sqrt(D * t) * u)) * 2 * u
    frac_l = np.trapezoid(integrand, u)
    out(f"   t = {t:8.0f}:  Kreis {frac_c:.5f}  (-> a/L = {a_w/L})   Gerade {frac_l:.5f}  (≈ a/√(πDt) = {a_w/math.sqrt(PI*D*t):.5f})")

sub("1d Liouville-Fluss ẋ = v (Gauß-Ensemble σ_v = 1)")
for t in [1.0, 10.0, 100.0]:
    fl = math.erf(L / (2 * math.sqrt(2) * t))
    dev = 2 * sum(math.exp(-2 * PI ** 2 * m_ ** 2 * t ** 2) for m_ in range(1, 50))
    out(f"   t = {t:5.0f}: Gerade im Fenster L: {fl:.5f};  Kreis: 1, max|ρ·L - 1| = {dev:.2e} (Phasenmischung, kein Verlust)")
out("Koopman-Spektrum: Kreis rein punktförmig {i v 2πn/L}; Gerade kontinuierlich {i v k}.")
alpha = (math.sqrt(5) - 1) / 2
Nst = 2_000_000
xs = (np.arange(Nst) * alpha) % 1.0
for a_k in [0.01, 0.001]:
    visits = np.nonzero(xs < a_k)[0]
    mean_ret = np.mean(np.diff(visits))
    out(f"Kac, Kreisdrehung um den goldenen Winkel, Menge A = [0,{a_k}): mittlere Rückkehrzeit {mean_ret:.3f}"
        f"  (Kac 1/|A| = {1/a_k:.0f});  offene Translation x -> x + α auf R: 0 Rückkehrer")

sub("1e Wellenpaket (ħ = m = 1): Ring der Länge L gegen Gerade")
ng = 2048
x = np.arange(ng) * L / ng
sig0 = 0.05
k0 = 2 * PI * 10 / L
psi0 = np.exp(-(x - 0.5) ** 2 / (4 * sig0 ** 2) + 1j * k0 * x)
psi0 /= np.linalg.norm(psi0)
kgrid = 2 * PI * np.fft.fftfreq(ng, d=L / ng)
Trev = L ** 2 / PI
Mbig = 64
ngb = ng * Mbig
xb = np.arange(ngb) * L / ng
psib = np.exp(-(xb - Mbig * L / 2) ** 2 / (4 * sig0 ** 2) + 1j * k0 * xb)
psib /= np.linalg.norm(psib)
kb = 2 * PI * np.fft.fftfreq(ngb, d=L / ng)
out(f"Wiederkehrzeit T_rev = mL²/(πħ) = {Trev:.5f}")
for frac in [0.25, 0.5, 1.0, 2.0]:
    t = frac * Trev
    pr = np.fft.ifft(np.exp(-1j * kgrid ** 2 * t / 2) * np.fft.fft(psi0))
    pb = np.fft.ifft(np.exp(-1j * kb ** 2 * t / 2) * np.fft.fft(psib))
    Fr = abs(np.vdot(psi0, pr)) ** 2
    Fb = abs(np.vdot(psib, pb)) ** 2
    out(f"   t = {frac:4.2f} T_rev:  Ring |<ψ0|ψ(t)>|² = {Fr:.6f}    Gerade (Ring 64L) = {Fb:.3e}")
out(f"   Gerade analytisch (k0 = 0-Anteil): 1/√(1+(t/2σ0²)²) bei T_rev = {1/math.sqrt(1+(Trev/(2*sig0**2))**2):.4e}")

sub("1f Pólya: Irrfahrt auf Z³ gegen Torus Z_N³")
u_w = math.sqrt(6) / (32 * PI ** 3) * math.gamma(1 / 24) * math.gamma(5 / 24) * math.gamma(7 / 24) * math.gamma(11 / 24)
tt = np.concatenate([np.linspace(0, 10, 20001), np.linspace(10, 600, 200001)[1:]])
integ = np.exp(-tt) * np.i0(tt / 3) ** 3
u_num = np.trapezoid(integ, tt) + (3 / (2 * PI)) ** 1.5 * 2 / math.sqrt(600)
out(f"Watson-Integral u = {u_w:.6f} (numerisch {u_num:.6f});  Rückkehrwahrscheinlichkeit Z³: 1 - 1/u = {1-1/u_w:.6f}")
rng = np.random.default_rng(1)
for Nt in [4, 6]:
    W = 20000
    pos = np.zeros((W, 3), dtype=np.int64)
    ret = np.full(W, -1)
    steps_max = 40 * Nt ** 3
    dirs = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    alive = np.ones(W, bool)
    for s_ in range(1, steps_max + 1):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        pos[idx] = (pos[idx] + dirs[rng.integers(0, 6, idx.size)]) % Nt
        back = idx[np.all(pos[idx] == 0, axis=1)]
        ret[back] = s_
        alive[back] = False
    out(f"Torus Z_{Nt}³: zurückgekehrt {np.mean(ret > 0)*100:.2f} % (bis {steps_max} Schritte), "
        f"mittlere Rückkehrzeit {np.mean(ret[ret>0]):.1f}  (Kac N³ = {Nt**3})")

sub("1g Bindung als Funktion der Schließungslänge L (D = 1)")
out("   L    | κ_var = 12D/L² | κ_gap = 4π²D/L² | Fluchtrate Intervall π²D/L² | Kac-Rückkehr Gitter h = 0.01 | Rückkehrdichte t=100")
for Lx in [1.0, 2.0, 5.0, 10.0, 100.0]:
    nn_ = np.arange(1, 400)
    rc = (1 / Lx) * (1 + 2 * np.sum(np.exp(-4 * PI ** 2 * nn_ ** 2 * D * 100 / Lx ** 2)))
    out(f" {Lx:6.1f} | {12*D/Lx**2:13.5f} | {4*PI**2*D/Lx**2:14.5f} | {PI**2*D/Lx**2:26.5f} | {int(round(Lx/0.01)):9d} Schritte             | {rc:.5f}")
out("   (Kac auf dem Gitter Z_{L/h}: mittlere Rückkehrzeit = L/h Schritte; auf Z: unendlich; Rückkehrdichte bei L = 100 noch nicht relaxiert, t < L²/4π²D)")
out("Toy Model mit K ≡ 1 (τ_K = ∞): Umkehrpunkt aus ½γ(x-x0)² - GM/x = H0, GM = 1, x0 = 1, v_s = 1.6, H0 = 0.28")
H0 = 0.5 * 1.6 ** 2 - 1.0
for gam in [0.05, 0.01, 1e-3, 1e-4, 1e-6]:
    lo, hi = 1.0, 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * gam * (mid - 1) ** 2 - 1 / mid - H0 > 0:
            hi = mid
        else:
            lo = mid
    xm = 0.5 * (lo + hi)
    lB = xm - 1
    out(f"   γ = {gam:8.1e}: x_max = {xm:12.4f},  ℓ_B = x_max - x0 = {lB:11.4f},  γ ℓ_B²/(2 H0) = {gam*lB**2/(2*H0):.5f}")
out("=> Schließungsstärke γ ≈ 2 H0/ℓ_B²: dieselbe 1/L²-Form wie κ_eff ≈ c·D/L².")

# ===========================================================================
hr("TEIL 2  Die drei Rollen als Kräfte am Torus (zwei Stuart-Landau-Zyklen)")
# ===========================================================================
out("ż_j = (μ + iω_j) z_j - β|z_j|² z_j + K e^{iφ}(z_k - z_j) [+ √(2D) ξ_j]")
out("Polarform: Ȧ1 = μA1 - βA1³ + K[A2 cos(ψ+φ) - A1 cos φ],  θ̇1 = ω1 + K[(A2/A1) sin(ψ+φ) - sin φ],  ψ = θ2 - θ1")


def sl_run(mu, beta, w1, w2, K, phi, T=400.0, dt=0.01, z0=None, t_rec=200.0):
    mu, beta, w1, w2, K, phi = [np.atleast_1d(np.asarray(v_, float)) for v_ in (mu, beta, w1, w2, K, phi)]
    nc = max(len(mu), len(K), len(w1))
    c = K * np.exp(1j * phi)

    def f(z):
        z1, z2 = z
        d1 = (mu + 1j * w1) * z1 - beta * np.abs(z1) ** 2 * z1 + c * (z2 - z1)
        d2 = (mu + 1j * w2) * z2 - beta * np.abs(z2) ** 2 * z2 + c * (z1 - z2)
        return np.array([d1, d2])
    z = np.array([np.full(nc, 0.9 + 0j), np.full(nc, 0.9 * np.exp(1j * 2.0))]) if z0 is None else z0.copy()
    nsteps = int(T / dt)
    th = np.zeros((2, nc))
    A_acc = np.zeros((2, nc)); A_min = np.full((2, nc), 1e9); A_max = np.zeros((2, nc))
    e_acc = np.zeros(nc, complex)
    cnt = 0
    t_rec_steps = int(t_rec / dt)
    for s_ in range(nsteps):
        k1 = f(z); k2 = f(z + 0.5 * dt * k1); k3 = f(z + 0.5 * dt * k2); k4 = f(z + dt * k3)
        zn = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if s_ >= nsteps - t_rec_steps:
            th += np.angle(zn / z)
            Aa = np.abs(zn)
            A_acc += Aa; A_min = np.minimum(A_min, Aa); A_max = np.maximum(A_max, Aa)
            e_acc += np.exp(1j * (np.angle(zn[1]) - np.angle(zn[0])))
            cnt += 1
        z = zn
    return dict(Omega=th / (cnt * dt), A=A_acc / cnt, Amin=A_min, Amax=A_max, r=np.abs(e_acc / cnt), z=z)


def classify(mu, K, Dw, phi=0.0):
    Kd = K * math.cos(phi)
    if Kd > mu and abs(Dw) > 2 * math.sqrt(max(mu * (2 * Kd - mu), 0.0)):
        return "Zerfall (Amplitudentod, Fixpunkt)"
    if abs(Dw) <= 2 * abs(Kd):
        return "eingerastet T¹"
    return "freier Torus T²"


sub("2a Gleichgewichte (μ = β = 1, ω̄ = 10, φ = 0) – Vorhersage gegen Simulation")
cases = [("A", 0.1, 1.00, 0.0), ("B", 0.1, 0.21, 0.0), ("C", 0.3, 0.50, 0.0),
         ("D", 1.5, 4.00, 0.0), ("E", 1.5, 2.00, 0.0), ("F", 0.3, 0.50, PI / 2)]
mu0, be0, wbar = 1.0, 1.0, 10.0
Ks = np.array([c_[1] for c_ in cases]); Dws = np.array([c_[2] for c_ in cases]); phs = np.array([c_[3] for c_ in cases])
res = sl_run(mu0, be0, wbar - Dws / 2, wbar + Dws / 2, Ks, phs, T=600.0, t_rec=300.0)
out(" Fall |  K   |  Δω  |  φ    | Vorhersage                         | A_pred  | A_sim (min..max)        | Ω2-Ω1 pred | Ω2-Ω1 sim | r pred | r sim")
for i, (nm, K, Dw, ph) in enumerate(cases):
    cl = classify(mu0, K, Dw, ph)
    Kd = K * math.cos(ph)
    if cl.startswith("Zerfall"):
        Ap, beat, rp = 0.0, float('nan'), float('nan')
    elif cl.startswith("einger"):
        s_ = Dw / (2 * Kd); cp = math.sqrt(1 - s_ ** 2)
        Ap = math.sqrt(max(mu0 - Kd * (1 - cp), 0) / be0); beat = 0.0; rp = 1.0
    else:
        Ap = math.sqrt(max(mu0 - Kd, 0) / be0)
        beat = math.sqrt(Dw ** 2 - 4 * Kd ** 2)
        rp = (Dw - beat) / (2 * Kd) if Kd > 0 else 0.0
    beat_s = res['Omega'][1, i] - res['Omega'][0, i]
    out(f"  {nm}   | {K:4.2f} | {Dw:4.2f} | {ph:5.3f} | {cl:34s} | {Ap:7.4f} | {res['A'][0,i]:.4f} ({res['Amin'][0,i]:.3f}..{res['Amax'][0,i]:.3f}) |"
        f" {beat:10.4f} | {beat_s:9.4f} | {rp:6.3f} | {res['r'][i]:.3f}")
out("   A_pred: freier Torus √((μ-K cosφ)/β) (mit <cos ψ> = 0), eingerastet √((μ - K(1-cos ψ*))/β), sin ψ* = Δω/2K")
out("   Schwebung Adler √(Δω² - 4K²cos²φ);  Phasenkonzentration r = |<e^{iψ}>| = (Δω - √(Δω²-4K²))/(2K)")
out(f"   Zerfallsgrenze für K = 1.5: Δω_c = 2√(μ(2K-μ)) = {2*math.sqrt(1*(3-1)):.4f}")

sub("2b Dicke des Torus: Rauschen gegen transversale Rückstellung (stochastisches Heun, 400 Realisierungen,")
out("    mitrotierendes System ω̄ = 0; Rauschen: Re und Im je Varianz D·dt)")
rngn = np.random.default_rng(7)
for nm, K, Dw, Dn in [("A", 0.1, 1.0, 1e-3), ("C", 0.3, 0.5, 1e-3)]:
    Wn = 400
    base = sl_run(mu0, be0, -Dw / 2, Dw / 2, K, 0.0, T=300.0, t_rec=10.0)
    dt = 0.005

    def drift(z):
        z1, z2 = z
        d1 = (mu0 - 1j * Dw / 2) * z1 - be0 * np.abs(z1) ** 2 * z1 + K * (z2 - z1)
        d2 = (mu0 + 1j * Dw / 2) * z2 - be0 * np.abs(z2) ** 2 * z2 + K * (z1 - z2)
        return np.array([d1, d2])
    for Dcur in [0.0, Dn]:
        z = np.repeat(base['z'], Wn if Dcur > 0 else 1, axis=1)
        W_ = z.shape[1]
        Asum = 0.0; A2sum = 0.0; psis = []; nsamp = 0
        nst = int(150 / dt)
        for s_ in range(nst):
            noise = math.sqrt(Dcur * dt) * (rngn.normal(size=(2, W_)) + 1j * rngn.normal(size=(2, W_)))
            f0 = drift(z)
            zp = z + dt * f0 + noise
            z = z + 0.5 * dt * (f0 + drift(zp)) + noise
            if s_ > nst // 3 and s_ % 20 == 0:
                A1 = np.abs(z[0]); Asum += np.sum(A1); A2sum += np.sum(A1 ** 2); nsamp += W_
                psis.append(np.angle(z[1] / z[0]))
        Am = Asum / nsamp
        varA = A2sum / nsamp - Am ** 2
        if Dcur == 0.0:
            varA_det = varA
            continue
        # radiale Diffusion D/2, Rückstellrate |f'(A*)| = 2βA*²  ->  Var = (D/2)/(2βA*²)
        varA_pred = (Dn / 2) / (2 * be0 * Am ** 2)
        out(f"  Fall {nm}: <A> = {Am:.4f},  Var(A) sim = {varA:.3e} (davon deterministische Schwebung {varA_det:.3e}),"
            f"  Rauschanteil {varA-varA_det:.3e},  Vorhersage (D/2)/(2βA*²) = {varA_pred:.3e},"
            f"  relative Dicke = {math.sqrt(max(varA-varA_det,0))/Am:.4f}")
        if nm == "C":
            ps = np.array(psis).ravel()
            sp = Dw / (2 * K); cp = math.sqrt(1 - sp ** 2)
            psis_c = np.angle(np.exp(1j * (ps - math.asin(sp))))
            varpsi = np.var(psis_c)
            # Phasendiffusion je Zyklus D/(2A²), für ψ: D/A²; Rückstellrate 2K cos ψ*
            varpsi_pred = (Dn / Am ** 2) / (2 * K * cp)
            out(f"          Phasendicke der Rastkurve: Var(ψ) sim = {varpsi:.3e}, Vorhersage (D/A²)/(2K cos ψ*) = {varpsi_pred:.3e}")

sub("2c Parametrischer Pfad (Temperatur ϑ): μ(ϑ) = 1 - 0.1ϑ, Δω(ϑ) = 0.6 + 0.05ϑ, K = 0.4, D(ϑ) = 1e-3(1+0.1ϑ)")
Kp = 0.4
thetas = np.array([0.0, 2.0, 3.8, 4.2, 5.0, 5.8, 6.2, 7.0])
mus = 1 - 0.1 * thetas
Dws = 0.6 + 0.05 * thetas
resp = sl_run(mus, be0, wbar - Dws / 2, wbar + Dws / 2, np.full(len(thetas), Kp), np.zeros(len(thetas)), T=800.0, t_rec=300.0)
out("   ϑ   |  μ   |  Δω  | Regime (analytisch)                | A_sim  | Ω2-Ω1 sim | r sim | rel. Dicke √(D/2/(2(μ-K)))/A")
for i, th_ in enumerate(thetas):
    cl = classify(mus[i], Kp, Dws[i])
    Dn = 1e-3 * (1 + 0.1 * th_)
    if cl.startswith("frei") and mus[i] > Kp:
        thick = math.sqrt(Dn / 2 / (2 * (mus[i] - Kp))) / math.sqrt((mus[i] - Kp) / be0)
        ts = f"{thick:.4f}"
    else:
        ts = "  –"
    if cl.startswith("Zerfall"):
        bs, rs_ = "      –  ", "  –  "
    else:
        bs, rs_ = f"{resp['Omega'][1,i]-resp['Omega'][0,i]:9.4f}", f"{resp['r'][i]:.3f}"
    out(f"  {th_:4.1f} | {mus[i]:4.2f} | {Dws[i]:4.2f} | {cl:34s} | {resp['A'][0,i]:.4f} | {bs} | {rs_} | {ts}")
out("   analytische Schwellen: Entriegeln Δω = 2K -> ϑ = 4.0;  Zerfall K > μ -> ϑ = 6.0 (dort Δω = 0.9 > 2√(μ(2K-μ)) = 0.8)")

sub("2d Zahlenbeispiel (Laborskala, Annahmen): zwei selbsterregte MEMS-Oszillatoren bei f0 = 1 MHz")
mu_Hz = 50.0
out(f"   Einheit μ0 = 2π·{mu_Hz:.0f} Hz (Nettoverstärkung), Zeiteinheit 1/μ0 = {1e3/(2*PI*mu_Hz):.3f} ms, 1 ϑ-Einheit = 10 K")
out(f"   K = 0.4 μ0 = 2π·{0.4*mu_Hz:.0f} Hz (dissipative Kopplung über gemeinsame Einspannung)")
out(f"   Rastbereich |Δf| ≤ 2K/2π = {2*0.4*mu_Hz:.0f} Hz;  kalt Δf = {0.6*mu_Hz:.0f} Hz -> eingerastet")
out(f"   differentielle Frequenzdrift {0.05*mu_Hz/10:.2f} Hz/K, Verstärkungsabfall {0.1*mu_Hz/10:.1f} Hz/K (je 2π)")
out(f"   Entriegeln bei ΔT = 40 K (Schwebung setzt mit √-Gesetz ein), Zerfall bei ΔT = 60 K")

# ===========================================================================
hr("TEIL 3  Das Dreieck Energie – Struktur – Zustand als 3x3-Kopplungssystem")
# ===========================================================================
out("J = diag(1/τ_E, 1/τ_S, 1/τ_Z)·(G - I),  G_ij = dimensionslose Verstärkung j -> i (Diagonale 0)")
out("a1 = Σ 1/τ_i,  a2 = Σ_{i<j} (1 - g_ij g_ji)/(τ_i τ_j),  a3 = det(I - G)/(τ_E τ_S τ_Z)")
out("det(I - G) = 1 - q_ES - q_EZ - q_SZ - q_ESZ - q_EZS  (2-Zyklen und beide 3-Zyklen)")


def hurwitz(G, tau):
    Dg = np.diag(1 / tau)
    Jm = Dg @ (G - np.eye(3))
    a1 = np.sum(1 / tau)
    a2 = sum((1 - G[i, j] * G[j, i]) / (tau[i] * tau[j]) for i in range(3) for j in range(i + 1, 3))
    q2 = G[0, 1] * G[1, 0] + G[0, 2] * G[2, 0] + G[1, 2] * G[2, 1]
    q3 = G[1, 0] * G[2, 1] * G[0, 2] + G[2, 0] * G[1, 2] * G[0, 1]
    a3 = (1 - q2 - q3) / np.prod(tau)
    return Jm, a1, a2, a3


sub("3a Zufallstest: Hurwitz aus Schleifenverstärkungen gegen Eigenwerte")
rng3 = np.random.default_rng(3)
Ntest = 20000
agree = 0; n_st = 0; n_real = 0; n_hopf = 0; max_a3err = 0.0
for _ in range(Ntest):
    G = rng3.normal(0, 1.0, (3, 3)); np.fill_diagonal(G, 0)
    tau = 10 ** rng3.uniform(-1, 3, 3)
    Jm, a1, a2, a3 = hurwitz(G, tau)
    ev = np.linalg.eigvals(Jm)
    st_e = np.max(ev.real) < 0
    st_h = (a1 > 0) and (a3 > 0) and (a1 * a2 > a3)
    max_a3err = max(max_a3err, abs(a3 + np.linalg.det(Jm)) / max(abs(a3), 1e-300))
    agree += (st_e == st_h)
    if st_e:
        n_st += 1
    elif a3 <= 0:
        n_real += 1
    else:
        n_hopf += 1
out(f"   {Ntest} Zufallssysteme (g ~ N(0,1), τ log-gleichverteilt 0.1..1000): Übereinstimmung {agree/Ntest*100:.3f} %")
out(f"   max rel. Fehler a3 = -det J: {max_a3err:.1e};  stabil {n_st}, instabil über a3 <= 0 (Durchgehen, Σq >= 1) {n_real},"
    f" instabil oszillatorisch (a1a2 < a3, Aufschaukeln) {n_hopf}")

sub("3b Negativer 3-Zyklus E -> S -> Z -> E (Vorzeichen +,+,-), Schwelle |q3| für Aufschaukeln")


def q3_threshold(tau):
    lo, hi = 0.0, 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        g = mid ** (1 / 3)
        G = np.zeros((3, 3)); G[1, 0] = g; G[2, 1] = g; G[0, 2] = -g
        ev = np.linalg.eigvals(np.diag(1 / tau) @ (G - np.eye(3)))
        if np.max(ev.real) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


for tau in [np.array([1.0, 1.0, 1.0]), np.array([1.0, 3.0, 10.0]), np.array([1.0, 10.0, 100.0]), np.array([1.0, 100.0, 1e4])]:
    thr = q3_threshold(tau)
    a1 = np.sum(1 / tau); a2 = sum(1 / (tau[i] * tau[j]) for i in range(3) for j in range(i + 1, 3))
    thr_an = a1 * a2 * np.prod(tau) - 1
    out(f"   τ = {tau}:  |q3|_krit numerisch = {thr:12.4f},  analytisch a1·a2·Πτ - 1 = {thr_an:12.4f}")
out("   Sekantenbedingung: |q3| < sec³(π/3) = 8 garantiert Stabilität für alle τ; Zeitskalentrennung hebt die Schwelle.")

sub("3c Nichtlineares E-S-Z-Modell")
out("Ė = 2(P·ℒ(s) - γ(Z)E),  ℒ(s) = 1/(1 + (Δ0 - s)²),  γ(Z) = 1 + bZ          (Struktur -> Energie: Resonanzlage)")
out("ṡ = (aE + dZ - s)/τ_S                                                    (Energie -> Struktur; Zustand -> Struktur)")
out("Ż = νE(1 + c_s s²) - hZ                                                  (Energie/Struktur -> Zustand; h = Ausheilen)")
Pp, aa, Delta0, tauS = 1.0, 5.0, 3.0, 10.0


def fixed_points(A):
    ss = np.linspace(-1, 8, 180001)
    fn = A / (1 + (Delta0 - ss) ** 2) - ss
    sg = np.where(fn >= 0, 1, -1)
    idx = np.nonzero(sg[:-1] != sg[1:])[0]
    return ss[idx]


def qloop(s, A):
    return A * 2 * (Delta0 - s) / (1 + (Delta0 - s) ** 2) ** 2


fps = fixed_points(aa * Pp)
out(f"Ohne Zustand (Z = 0): A = aP/γ = {aa*Pp}; Fixpunkte s* = {np.round(fps,4)}, Schleifenverstärkung q = {np.round([qloop(s_, aa*Pp) for s_ in fps],4)}")
# Falten: A(s) = s(1 + (Δ0 - s)²), dA/ds = 0  <=>  u = Δ0 - s = (Δ0 ± √(Δ0² - 3))/3
u_f = np.array([(Delta0 + math.sqrt(Delta0 ** 2 - 3)) / 3, (Delta0 - math.sqrt(Delta0 ** 2 - 3)) / 3])
s_f = Delta0 - u_f
A_f = s_f * (1 + u_f ** 2)
A_fold_lo = float(np.min(A_f))
out(f"Falten (q = 1): s_f = {np.round(s_f,4)}, A_f = {np.round(A_f,4)}; bistabil für {A_f.min():.4f} < A < {A_f.max():.4f} (Bedingung Δ0 > √3)")
out(f"Der heiße Ast (s > s_f) verschwindet, sobald A = aP/γ unter {A_fold_lo:.4f} fällt: Sattel-Knoten, Sprung auf den kalten Ast")


def esz_run(bZ, d, nu, h, c_s=0.0, T=4000.0, dt=0.01, s0=None, rec_every=100):
    E = Pp / 1.0 / (1 + (Delta0 - 3.0) ** 2); s = 3.0 if s0 is None else s0; Z = 0.0
    E = Pp / (1 + (Delta0 - s) ** 2)
    rec = []
    nst = int(T / dt)

    def f(E, s, Z):
        Lr = 1 / (1 + (Delta0 - s) ** 2)
        return (2 * (Pp * Lr - (1 + bZ * Z) * E), (aa * E + d * Z - s) / tauS, nu * E * (1 + c_s * s * s) - h * Z)
    for i in range(nst):
        k1 = f(E, s, Z)
        k2 = f(E + .5 * dt * k1[0], s + .5 * dt * k1[1], Z + .5 * dt * k1[2])
        k3 = f(E + .5 * dt * k2[0], s + .5 * dt * k2[1], Z + .5 * dt * k2[2])
        k4 = f(E + dt * k3[0], s + dt * k3[1], Z + dt * k3[2])
        E += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        s += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Z += dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        if i % rec_every == 0:
            rec.append((i * dt, E, s, Z))
    return np.array(rec)


# Start auf dem heißen Ast
s_hot = fps[-1]
out(f"Start auf dem heißen Ast s = {s_hot:.4f}")
out(" Szenario                                 | b    | d     | ν      | h      | Ende: E      s       Z       | Kippen bei t")
scen = [("(i)  Ausheilen, Attraktor", 0.5, 0.0, 1e-3, 5e-3),
        ("(ii) Alterung ohne Ausheilen", 0.5, 0.0, 1e-3, 0.0),
        ("(iii) Alterung, schneller (ν x 4)", 0.5, 0.0, 4e-3, 0.0),
        ("(iv) Z schiebt Struktur weg (d > 0)", 0.0, 0.5, 1e-3, 0.0),
        ("(v)  Z schiebt Struktur zurück (d < 0)", 0.0, -0.5, 1e-3, 0.0)]
tips = {}
for nm, bZ, d, nu, hh_ in scen:
    rec = esz_run(bZ, d, nu, hh_, s0=s_hot, T=6000.0, dt=0.02, rec_every=50)
    jump = np.nonzero(rec[:, 2] < 2.0)[0]
    ttip = rec[jump[0], 0] if jump.size else float('nan')
    tips[nm] = (ttip, rec)
    out(f" {nm:40s} | {bZ:4.2f} | {d:5.2f} | {nu:6.0e} | {hh_:6.0e} | {rec[-1,1]:.4f}  {rec[-1,2]:7.4f}  {rec[-1,3]:7.4f} | {ttip:8.1f}")
# analytische Kippzeit (ii): γ_c = aP/A_fold, Z_c = (γ_c - 1)/b, t ≈ Z_c/(ν <E>)
gam_c = aa * Pp / A_fold_lo
Zc = (gam_c - 1) / 0.5
rec =tips["(ii) Alterung ohne Ausheilen"][1]
tt_ii = tips["(ii) Alterung ohne Ausheilen"][0]
iz = np.searchsorted(rec[:, 0], tt_ii)
Emean = np.mean(rec[:iz, 1]) if iz > 0 else float('nan')
out(f"Analytisch (ii): γ_c = aP/A_fold = {gam_c:.4f}, Z_c = (γ_c-1)/b = {Zc:.4f}, gemessen Z(t_kipp) = {rec[iz,3]:.4f};"
    f"  t ≈ Z_c/(ν<E>) = {Zc/(1e-3*Emean):.1f} mit <E> = {Emean:.4f}")
t3 = tips["(iii) Alterung, schneller (ν x 4)"][0]
out(f"Kippzeit-Verhältnis (ii)/(iii) = {tt_ii/t3:.3f}  (erwartet ν-Verhältnis 4)")
# Ausheilschwelle (i): Fixpunkt existiert auf dem heißen Ast, wenn aP/(1 + b ν E*/h) > A_fold
Ehot = Pp / (1 + (Delta0 - s_hot) ** 2)
hc = 0.5 * 1e-3 * Ehot / (gam_c - 1)
out(f"Ausheilschwelle (i): h_c ≈ bνE*/(γ_c - 1) = {hc:.2e} (E* ≈ {Ehot:.4f}); gewählt h = 5e-3 > h_c -> Attraktor")
# kritisches Verlangsamen: langsamster Eigenwert des schnellen (E,s)-Teilsystems entlang (ii)
jf = np.argmax(rec[:, 3] >= Zc)
t_fold = rec[jf, 0]
out(f"Faltendurchgang (Z = Z_c) bei t_F = {t_fold:.1f}; Sprung (s < 2) bei t = {tt_ii:.1f}: Geisterverzögerung {tt_ii - t_fold:.1f}"
    f" (vgl. τ_S = {tauS})")
out("Kritisches Verlangsamen vor dem Kippen (ii): langsamster Eigenwert λ_s des (E,s)-Blocks am quasistatischen Punkt")
out("    γ_c - γ  | t_F - t (Lauf) | q = A ℒ'(s*) | λ_s         | λ_s/√(γ_c - γ) | λ_s/√(t_F - t)")
for dg in [0.3, 0.1, 0.03, 0.01, 0.003, 0.001]:
    gq = gam_c - dg
    A = aa * Pp / gq
    fp = fixed_points(A)
    sq = fp[-1]
    Lr = 1 / (1 + (Delta0 - sq) ** 2)
    dL = 2 * (Delta0 - sq) * Lr ** 2
    Jfs = np.array([[-2 * gq, 2 * Pp * dL], [aa / tauS, -1 / tauS]])
    lam_s = np.max(np.linalg.eigvals(Jfs).real)
    Zq = (gq - 1) / 0.5
    jq = np.argmax(rec[:, 3] >= Zq)
    dtq = t_fold - rec[jq, 0]
    out(f"   {dg:9.3f} | {dtq:14.1f} | {qloop(sq, A):11.5f} | {lam_s:.4e} | {lam_s/math.sqrt(dg):14.4e} | {lam_s/math.sqrt(max(dtq,1e-9)):.4e}")
out("   (λ_s/√(γ_c - γ) -> konstant: Wurzelgesetz des Sattel-Knotens; da Z nahezu linear wächst, auch λ_s ∝ √(t_F - t))")

sub("3d Gedächtniskern: Toy Model = Teilchen am Maxwell-Element")
out("Toy Model: ẍ = -ω0² x - γ y,  ẏ = ẋ - y/τ_K   <=>   Kraft = -∫ G(t-t') ẋ(t') dt',  G(t) = γ e^{-t/τ_K}")
out("Harmonische Antwort: Speicheranteil γ(ωτ)²/(1+(ωτ)²) (bindet), Verlustanteil γτ/(1+(ωτ)²) (reibt)")
w0, gm = 1.0, 0.1
out("   ωτ_K   | Δ(ω²)/γ = (ωτ)²/(1+(ωτ)²) | Dämpfung η/(γτ) = 1/(1+(ωτ)²) | sim: ω² - ω0²  | sim: Abklingrate")
for wt in [0.1, 1.0, 10.0, 100.0]:
    tauK = wt / w0
    # Eigenwerte des linearen Systems (x, v, y)
    Jm = np.array([[0, 1, 0], [-w0 ** 2, 0, -gm], [0, 1, -1 / tauK]])
    ev = np.linalg.eigvals(Jm)
    osc = ev[np.argmax(np.abs(ev.imag))]
    w2 = osc.imag ** 2 + osc.real ** 2
    wa = abs(osc.imag)
    out(f"   {wt:6.1f}  | {wt**2/(1+wt**2):25.5f}  | {1/(1+wt**2):28.5f}   | {(w2-w0**2)/gm:12.5f}·γ | {-osc.real:.5f}"
        f" (½γτ/(1+(ωτ)²) mit tatsächlichem ω = {wa:.4f}: {0.5*gm*tauK/(1+(wa*tauK)**2):.5f}; Speicher mit ω: {(wa*tauK)**2/(1+(wa*tauK)**2):.5f})")

# ===========================================================================
hr("TEIL 4  Wärme als langsames Gedächtnis: Rückwirkung auf die EM-Regeln")
# ===========================================================================
sub("4a Eingespannter Kupferdraht im Magnetfeld (magnetomotorisch, Spannungsquelle)")
rho_el, alR = 1.68e-8, 3.93e-3
E0, alE = 128e9, -3.8e-4
alth, rho, cp, kth = 16.5e-6, 8960.0, 385.0, 401.0
a_r, ell, Bf, Qm = 50e-6, 0.02, 0.02, 1000.0
sig0 = 20e6
A_c = PI * a_r ** 2
I_m = PI * a_r ** 4 / 4
R0 = rho_el * ell / A_c
G_th = 12 * kth * A_c / ell
C_th = rho * cp * A_c * ell
m_eff = 1.5 * rho * A_c * ell
Ncr = 4 * PI ** 2 * E0 * I_m / ell ** 2
N0 = sig0 * A_c
Gam = Bf * ell


def mech(dT, use_E=True, use_N=True):
    Et = E0 * (1 + alE * dT) if use_E else E0
    N = N0 - (E0 * A_c * alth * dT if use_N else 0.0)
    k_el = 8 * PI ** 4 * Et * I_m / ell ** 3
    k_geo = 2 * PI ** 2 * N / ell
    return k_el + k_geo


k00 = mech(0.0)
w_m0 = math.sqrt(k00 / m_eff)
c_m = m_eff * w_m0 / Qm
R_m = Gam ** 2 / c_m
out(f"Kupfer: ρ_el = {rho_el:.2e} Ωm, α_R = {alR:.2e}/K, E = {E0/1e9:.0f} GPa, (1/E)dE/dT ≈ {alE:.1e}/K, α = {alth:.2e}/K")
out(f"Draht: Radius {a_r*1e6:.0f} µm, Länge {ell*1e3:.0f} mm, beidseitig eingespannt, Vorspannung {sig0/1e6:.0f} MPa, B = {Bf} T, Q_m = {Qm:.0f}")
out(f"R0 = {R0*1e3:.3f} mΩ,  G_th = 12kA/ℓ = {G_th*1e3:.4f} mW/K,  C_th = {C_th*1e3:.4f} mJ/K,  τ_th = C/G = {C_th/G_th:.4f} s")
out(f"Galerkin (φ = 1 - cos 2πx/ℓ): m_eff = {m_eff*1e6:.4f} mg, N_cr = 4π²EI/ℓ² = {Ncr*1e3:.3f} mN, N0 = {N0*1e3:.2f} mN")
out(f"f_m0 = {w_m0/(2*PI):.2f} Hz,  Γ = Bℓ = {Gam:.4f} T·m,  R_m = Γ²/c = {R_m*1e3:.3f} mΩ,  R_m/R0 = {R_m/R0:.3f}")
dlnw_N = 0.5 * (-(2 * PI ** 2 / ell) * E0 * A_c * alth) / k00
dlnw_E = 0.5 * alE * (8 * PI ** 4 * E0 * I_m / ell ** 3) / k00
out(f"Frequenz-Temperaturkoeffizient: über Wärmespannung {dlnw_N:.4e}/K, über E(T) {dlnw_E:.4e}/K; Leitwert-Koeffizient α_R = {alR:.2e}/K")
Ql = Qm / (1 + R_m / R0)
out(f"belastete Güte Q_l = Q/(1 + R_m/R0) = {Ql:.1f};  Verschiebung pro K in Linienbreiten: {abs(dlnw_N)*Ql:.2f}")
out(f"Knicken durch Wärmespannung bei ΔT_b = (N0 + N_cr)/(EAα) = {(N0+Ncr)/(E0*A_c*alth):.3f} K;  τ_m = 2Q_l/ω = {2*Ql/w_m0:.4f} s")


def heat(dT, wd, V, use_R=True, use_E=True, use_N=True, eta_m=1.0):
    R = R0 * (1 + alR * dT) if use_R else R0
    k = mech(dT, use_E, use_N)
    Zm = Gam ** 2 / (c_m + 1j * (m_eff * wd - k / wd))
    I = V / (R + Zm)
    return 0.5 * abs(I) ** 2 * (R + eta_m * Zm.real)


def loop_parts(dT, wd, V):
    e = 1e-6
    q = (heat(dT + e, wd, V) - heat(dT - e, wd, V)) / (2 * e) / G_th
    parts = {}
    for nm, kw in [("σ(T)", dict(use_E=False, use_N=False)), ("E(T)", dict(use_R=False, use_N=False)),
                   ("Wärmespannung", dict(use_R=False, use_E=False))]:
        # nur dieser Kanal temperaturabhängig, übrige am Arbeitspunkt eingefroren
        def hp(x):
            R = R0 * (1 + alR * (x if kw.get('use_R', True) else dT))
            Et = E0 * (1 + alE * (x if kw.get('use_E', True) else dT))
            N = N0 - E0 * A_c * alth * (x if kw.get('use_N', True) else dT)
            k = 8 * PI ** 4 * Et * I_m / ell ** 3 + 2 * PI ** 2 * N / ell
            Zm = Gam ** 2 / (c_m + 1j * (m_eff * wd - k / wd))
            I = V / (R + Zm)
            return 0.5 * abs(I) ** 2 * (R + Zm.real)
        parts[nm] = (hp(dT + e) - hp(dT - e)) / (2 * e) / G_th
    return q, parts


def steady_states(wd, V, Tmax=3.0):
    xs2 = np.linspace(0.0, Tmax, 100001)
    fnv = heat(xs2, wd, V) - G_th * xs2
    sg = np.where(fnv >= 0, 1, -1)
    idx = np.nonzero(sg[:-1] != sg[1:])[0]
    return [xs2[i] - fnv[i] * (xs2[i + 1] - xs2[i]) / (fnv[i + 1] - fnv[i]) for i in idx]


out("")
out("Stationäre Temperaturen G ΔT = P(ΔT) und Schleifenverstärkung q = P'(ΔT)/G (zerlegt), Anregung f_d relativ zu f_m0:")
out("  V [mV] | f_d - f_m0 [Hz] | Lösungen ΔT [K]              | q je Lösung                 | q_σ(T)      | q_E(T)      | q_Wärmespannung")
for V, dfs in [(2e-3, [-10.0, -2.0, 0.0, 2.0]), (5e-3, [-10.0, -7.0, -4.0, 0.0]), (10e-3, [-40.0, -30.0, -25.0, -20.0, -10.0])]:
    for dfd in dfs:
        wd = w_m0 + 2 * PI * dfd
        roots = steady_states(wd, V)
        qs = [loop_parts(r_, wd, V) for r_ in roots]
        rs = ", ".join(f"{r_:.4f}" for r_ in roots)
        qt = ", ".join(f"{q_[0]:+.3f}" for q_ in qs)
        q0 = qs[0][1]
        out(f"   {V*1e3:4.1f}  | {dfd:+7.1f}         | {rs:28s} | {qt:27s} | {q0['σ(T)']:+.2e} | {q0['E(T)']:+.2e} | {q0['Wärmespannung']:+.3e}")
out("   (zerlegte Anteile an der ersten Lösung; eine Lösung mit q > 1 ist instabil -> Bistabilität, Sprung, Hysterese)")
out("Bistabilitätsfenster (drei stationäre Lösungen), Scan f_d - f_m0 von -80 bis +10 Hz in 0.25 Hz:")
for V in [2e-3, 5e-3, 6e-3, 7e-3, 10e-3]:
    win = []
    qmax = -1e9
    for dfd in np.arange(-80.0, 10.0001, 0.25):
        wd = w_m0 + 2 * PI * dfd
        roots = steady_states(wd, V)
        for r_ in roots:
            e = 1e-6
            qmax = max(qmax, (heat(r_ + e, wd, V) - heat(r_ - e, wd, V)) / (2 * e) / G_th)
        if len(roots) >= 3:
            win.append(dfd)
    ws = f"[{min(win):.2f}, {max(win):.2f}] Hz" if win else "keins"
    out(f"   V = {V*1e3:4.1f} mV: max q über alle Lösungen = {qmax:+.3f};  Fenster: {ws}")
V = 10e-3
out(f"Quasistatischer Sweep bei V = {V*1e3:.0f} mV (C dT/dt = P - G ΔT, je 0.25-Hz-Schritt 10 τ_th):")
dT = steady_states(w_m0 + 2 * PI * (-80.0), V)[0]
freqs = np.arange(-80.0, 10.0001, 0.25)


def relax(dT, wd, V, n=1000):
    h_ = 10 * C_th / G_th / n
    for _ in range(n):
        dT += (heat(dT, wd, V) - G_th * dT) / C_th * h_
    return dT


up = []
for dfd in freqs:
    dT = relax(dT, w_m0 + 2 * PI * dfd, V); up.append(dT)
down = []
for dfd in freqs[::-1]:
    dT = relax(dT, w_m0 + 2 * PI * dfd, V); down.append(dT)
down = np.array(down[::-1]); up = np.array(up)
diffm = np.abs(up - down)
hyst = freqs[diffm > 1e-3]
if hyst.size:
    out(f"   Hysterese (|ΔT_auf - ΔT_ab| > 1 mK): f_d - f_m0 ∈ [{hyst.min():.2f}, {hyst.max():.2f}] Hz, max Sprung {diffm.max()*1e3:.1f} mK")
else:
    out("   keine Hysterese")
jmax = np.argmax(diffm)
out(f"   bei f_d - f_m0 = {freqs[jmax]:.2f} Hz: ΔT_auf = {up[jmax]*1e3:.1f} mK, ΔT_ab = {down[jmax]*1e3:.1f} mK")
# Verhältnis der Kanäle, allgemein
out(f"Kanalverhältnis (Frequenzhub pro K in Linienbreiten) / α_R = {abs(dlnw_N)*Ql/alR:.0f}")

sub("4b Ferrit am Curie-Punkt (Modell): P(T) = P0 (1 - T/T_C)^{2β}, 2β = 0.72, T_C = 493 K (MnZn typisch), T0 = 293 K")
Tc, T0 = 493.0, 293.0
out("  P0/G [K] | T* [K]   | T_C - T* [K] | q = P'(T*)/G | Verstärkung 1/(1-q) | dT*/d ln P0 [K] | Kupfer (P konst.) dT*/d ln P0 [K]")
for P0G in [100.0, 300.0, 1000.0, 3000.0, 10000.0]:
    lo, hi = T0, Tc
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if P0G * (1 - mid / Tc) ** 0.72 - (mid - T0) > 0:
            lo = mid
        else:
            hi = mid
    Ts = 0.5 * (lo + hi)
    q = -P0G * 0.72 * (1 - Ts / Tc) ** (-0.28) / Tc
    dTd = P0G * (1 - Ts / Tc) ** 0.72 / (1 - q)
    out(f"  {P0G:8.0f} | {Ts:8.3f} | {Tc-Ts:12.4f} | {q:12.3f} | {1/(1-q):19.5f} | {dTd:15.4f} | {P0G:10.1f}")

# ===========================================================================
hr("TEIL 5  Gravitation: Vorspannung -> Steifigkeit -> Resonanzverschiebung")
# ===========================================================================
Gn, cl, Msun = 6.674e-11, 2.998e8, 1.989e30
bodies = [("Erde", 5.972e24, 6.371e6), ("Sonne", Msun, 6.957e8), ("Sirius B", 1.018 * Msun, 5.8e6),
          ("Neutronenstern", 1.4 * Msun, 12e3)]
out("Galerkin, Kragbalken (unten frei), φ = 1 - cos(πx/2ℓ): k_el = (π⁴/32) EI/ℓ³, k_geo = (π²/16 - 1/4) ρAg,")
out("   m_eff = (3/2 - 4/π) ρAℓ.  hängend: ω² = (k_el + k_geo)/m, stehend: ω² = (k_el - k_geo)/m")
kel_c = PI ** 4 / 32
kgeo_c = PI ** 2 / 16 - 0.25
out(f"   Koeffizienten: k_el = {kel_c:.4f} EI/ℓ³, k_geo = {kgeo_c:.4f} ρAg (unabhängig von ℓ); Knicklänge stehend ℓ_c³ = {kel_c/kgeo_c:.3f} EI/(ρAg)"
    f"  (exakt Greenhill 7.837, Galerkin +{(kel_c/kgeo_c/7.837-1)*100:.1f} %)")
structs = [("Si-Kragbalken 100 µm x 1 µm", 169e9, 2330.0, 100e-6, (1e-6) ** 2 / 12),
           ("Stahlstab 1 m, d = 10 mm", 200e9, 7850.0, 1.0, (10e-3) ** 2 / 16),
           ("Stahldraht 1 m, d = 0.1 mm", 200e9, 7850.0, 1.0, (0.1e-3) ** 2 / 16)]
out("")
out(" Körper          | g [m/s²]  | Δs [m]     | Takt √(1-6Δs/r) | h_max Stahl σ_y/ρg | ℓ_c Stab d=10mm | χ_g (hängend): Si-Balken | Stahlstab | Stahldraht")
for bn, Mb, rb in bodies:
    g = Gn * Mb / rb ** 2
    ds = Gn * Mb / (3 * cl ** 2)
    takt = math.sqrt(1 - 6 * ds / rb)
    hmax = 250e6 / (7850 * g)
    chis = []
    for sn, Ey, rh, ln, IA in structs:
        kel = kel_c * Ey * IA / ln ** 3
        kgeo = kgeo_c * rh * g
        chis.append(kgeo / (kel + kgeo))
    lc = (7.837 * 200e9 * (10e-3) ** 2 / 16 / (7850 * g)) ** (1 / 3)
    out(f" {bn:15s} | {g:9.3e} | {ds:10.4e} | {takt:.12f} | {hmax:15.4e} m | {lc:12.4e} m | {chis[0]:.3e}               | {chis[1]:.3e} | {chis[2]:.4f}")
out("Hebelgesetz: ∂ln ω/∂ln g = ½ χ_g (hängend), -½ χ_g/(1 - 2χ_g) (stehend);  χ_g = k_geo/(k_el + k_geo)")
out("")
out("Verhältnis differentieller Kanal / Takt-Kanal bei radialer Verschiebung (FTM-Form g = 3c²Δs/r²):")
out("   R = (½ χ_g |d ln g/dr|) / (d ln Takt/dr) = χ_g (r - 6Δs)/(3Δs)")
out("   mit Eigenbeschleunigung a = g/√(1-6Δs/r) (ART, statischer Beobachter): R' = χ_g[(r - 6Δs)/(3Δs) + ½]")
for bn, Mb, rb in bodies:
    ds = Gn * Mb / (3 * cl ** 2)
    Rr = (rb - 6 * ds) / (3 * ds)
    out(f"   {bn:15s}: r/Δs = {rb/ds:12.4e},  R/χ_g = {Rr:12.4e},  R'/χ_g = {Rr+0.5:12.4e}")
out("   R = 1 bei χ_g = 1: FTM-Form r = 9Δs (= 3GM/c², Radius der Photonensphäre), Eigenbeschleunigung r = 7.5Δs")
out("")
out("Beispiel: Pendel/hängender Draht (χ_g ≈ 1) 1 m anheben auf der Erde:")
ds_e = Gn * 5.972e24 / (3 * cl ** 2); r_e = 6.371e6
dlg = -2 * 1.0 / r_e
out(f"   δω/ω differentiell = ½ δln g = {0.5*dlg:.4e};  Takt δln(dτ/dt) = {3*ds_e/r_e**2/(1-6*ds_e/r_e)*1.0:.4e}")
for sn, Ey, rh, ln, IA in structs[:2]:
    g = 9.8203
    kel = kel_c * Ey * IA / ln ** 3; kgeo = kgeo_c * rh * g
    chi = kgeo / (kel + kgeo)
    out(f"   {sn}: χ_g = {chi:.3e}; hängend gegen stehend: ω_h/ω_s - 1 = {math.sqrt((kel+kgeo)/(kel-kgeo))-1:.3e}")

# ===========================================================================
hr("TEIL 6  Neue Aussagen: Auswahlregeln und Grenzen für die Zahl der Zyklen")
# ===========================================================================
sub("6a Koexistenz zweier Zyklen an gemeinsamer Verstärkung (Lamb): İ_i = 2I_i(α_i - β_i I_i - θ_ij I_j)")
for th12, th21 in [(0.5, 0.5), (0.9, 0.9), (1.5, 1.5), (2.0, 0.3)]:
    I = np.array([0.01, 0.02]); al = np.array([1.0, 0.95]); be = np.array([1.0, 1.0])
    dt = 0.01
    for _ in range(40000):
        dI = 2 * I * (al - be * I - np.array([th12 * I[1], th21 * I[0]]))
        I = np.maximum(I + dt * dI, 0)
    Cc = th12 * th21 / (be[0] * be[1])
    win_ok = (th12 / be[1] < al[0] / al[1] < be[0] / th21)
    out(f"   θ12 = {th12}, θ21 = {th21}: C = {Cc:.3f}, Fenster θ12/β2 = {th12/be[1]:.3f} < α1/α2 = {al[0]/al[1]:.4f} < β1/θ21 = {be[0]/th21:.3f}: {win_ok}"
        f" -> I = ({I[0]:.4f}, {I[1]:.4f})  {'Koexistenz (T²)' if min(I) > 1e-3 else 'ein Zyklus verdrängt (T¹)'}")

sub("6c Kac-Budget: Zerlegung des Kreises in 10 gleiche Zellen, Drehung um den goldenen Winkel (2·10⁶ Schritte)")
sums = 0.0
mr = []
for i in range(10):
    vis = np.nonzero((xs >= i / 10) & (xs < (i + 1) / 10))[0]
    mr.append(np.mean(np.diff(vis)))
mr = np.array(mr)
out(f"   mittlere Rückkehrzeiten je Zelle: min {mr.min():.4f}, max {mr.max():.4f};  Σ 1/<τ_i> = {np.sum(1/mr):.6f}  (Kac: = 1)")
out("   offene Translation auf R: jede Zelle wird höchstens einmal besucht, Σ 1/<τ_i> = 0")

sub("6b Höchstzahl freier Zyklen an einem gemeinsamen Medium: θ̇_i = ω_i + k Σ_j sin(θ_j - θ_i), ω_i gleichabständig in Band W = 1")
out("   N  | k1 (erstes Einrasten) | (N-1)·k1 | k_s (vollständig synchron) | N·k_s  | Kuramoto 2W/π = 0.6366")
for Nn in [4, 8, 16, 32]:
    kgrid = np.linspace(0.3, 1.2, 91) / Nn
    om = np.linspace(-0.5, 0.5, Nn)
    th = np.tile(np.linspace(0, 2 * PI, Nn, endpoint=False), (len(kgrid), 1)) * 1.0
    kk_ = kgrid[:, None]

    def fk(th):
        Zm = np.sum(np.exp(1j * th), axis=1, keepdims=True)
        return om[None, :] + kk_ * np.imag(Zm * np.exp(-1j * th))
    dt = 0.05
    Ttot = 3000.0
    nst = int(Ttot / dt)
    th_start = None
    for s_ in range(nst):
        k1 = fk(th); k2 = fk(th + .5 * dt * k1); k3 = fk(th + .5 * dt * k2); k4 = fk(th + dt * k3)
        th = th + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if s_ == nst // 2:
            th_start = th.copy()
    Om = (th - th_start) / (Ttot / 2)
    ncl = []
    for row in Om:
        srt = np.sort(row)
        ncl.append(1 + int(np.sum(np.diff(srt) > 2e-3)))
    ncl = np.array(ncl)
    i1 = np.argmax(ncl < Nn)
    is_ = np.argmax(ncl == 1) if np.any(ncl == 1) else None
    k1v = kgrid[i1]
    ksv = kgrid[is_] if is_ is not None else float('nan')
    out(f"   {Nn:2d} | {k1v:21.4e} | {(Nn-1)*k1v:8.4f} | {ksv:26.4e} | {Nn*ksv:6.4f}")
out("   Gitterauflösung in k: 0.01/N")

out("")
out(f"Rechenzeit: {time.time() - t_start:.1f} s")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "torus_wechselkraefte.txt"), "w",
          encoding="utf-8") as fh:
    fh.write(_BUF.getvalue())
