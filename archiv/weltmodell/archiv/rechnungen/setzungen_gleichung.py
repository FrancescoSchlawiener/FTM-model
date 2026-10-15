"""Setzungen der Rückkehrgleichung auf eine zurückführen: Irrfahrt der Regel mit Schritt r₀ und c.

Frage: Welche der Setzungen (1) parabolisch, (2) ϖ_R² = D·τ_D, (3) D = r₀c/2, (4) eine Uhr für θ und T,
(5) Abschirmung lokal Γ/D, (6) Erneuerung mit einem Sinn folgen aus
"die Rückkehr ist eine Irrfahrt mit Schritt r₀ und Tempo c, Verlust je Schritt = Taktgefälle über einen Schritt"?
Konstanten und Formeln aus ftm_kern. Ausgabe: setzungen_gleichung.txt
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ftm_kern as K  # noqa: E402

rng = np.random.default_rng(20260929)
out = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def kopf(t):
    P("")
    P("=" * 110)
    P(t)
    P("=" * 110)


# ---------------------------------------------------------------- Gitter-Hilfen
def lam_gitter(p):
    """Beschränkte Lösung der Gitterrückkehr u_n = (1−p)(u_{n+1}+u_{n−1})/2: u_n = λⁿ."""
    z = 1 - p
    return (1 - math.sqrt(p * (2 - p))) / z


def q_gitter(m, p):
    """Anteil, der aus Tiefe m Schritten zum Ort zurückkehrt (Verlust p je Schritt, Gitter exakt)."""
    return math.exp(m * math.log(lam_gitter(p)))


def mc_rueckkehr(m, p_takt, n_walk, L=None, theta_p=0.0, renew=False, maxstep=400000):
    """Monte-Carlo auf dem Gitter in ϖ: Start in Tiefe m, ±1 je Takt, Verlust p_takt + theta_p je Takt,
    bei L reflektierend (geschlossen). renew=False: Anteil der Rückläufe; renew=True: Zahl der Rückläufe
    (nach jedem Rücklauf Neustart in Tiefe m) bis zum Verlust."""
    pos = np.full(n_walk, m, dtype=np.int64)
    alive = np.ones(n_walk, bool)
    count = np.zeros(n_walk)
    ploss = p_takt + theta_p
    for _ in range(maxstep):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        dead = rng.random(idx.size) < ploss
        alive[idx[dead]] = False
        idx = idx[~dead]
        step = np.where(rng.random(idx.size) < 0.5, 1, -1)
        pos[idx] += step
        if L is not None:
            over = pos[idx] > L
            pos[idx[over]] = 2 * L - pos[idx[over]]
        hit = idx[pos[idx] <= 0]
        count[hit] += 1
        if renew:
            pos[hit] = m
        else:
            alive[hit] = False
    return count.mean(), count.std() / math.sqrt(n_walk), int(alive.sum())


# ============================================================================================================
kopf("1  Irrfahrt mit Schritt h und Takt h/c: Diffusion, Kac, Wegfall von D")
# (1) parabolisch: Schritt ±h je Takt h/c ⇒ Varianz N·h² ⇒ D = h c/2 je Achse; Welle erst unterhalb h/c
n, N = 20000, 4000
x = (2 * rng.integers(0, 2, size=(n, N), dtype=np.int8) - 1).sum(axis=1, dtype=np.int64)
P(f"   1D, N = {N} Schritte: ⟨x²⟩/(N h²) = {np.mean(x.astype(float)**2)/N:.4f} (Diffusion D = h·c/2: 1)")
P("   Kurtosis der Verteilung nach N Schritten (Gauß: 3):", f"{np.mean(x.astype(float)**4)/np.mean(x.astype(float)**2)**2:.3f}")
P("   ⇒ (1) folgt: eine Irrfahrt ist auf Skalen ≫ h parabolisch (Telegraph-Übergang bei |p| ~ 2c/h).")
P("   Normierung: bewegt sich jede Achse in jedem Takt, ist D = h c/2 je Achse; wählt jeder Takt eine von d Achsen,")
P("   ist D = h c/(2d) (Raum + ϖ, d = 4: h c/8). Die Irrfahrt legt D also nur bis auf die Geometrie fest.")

# Kac auf der bei Λ geschlossenen Halbgeraden: ⟨T⟩ = 2L Takte = 2Λ/c, unabhängig von h
P("")
P("   Kac, Gitter 0…L reflektierend, Rückkehrzeit zum Ort (Monte-Carlo, 20000 Rückläufe):")
for L in (10, 40, 160):
    pos = np.zeros(20000, np.int64)
    T = np.zeros(20000)
    act = np.ones(20000, bool)
    t = 0
    pos[:] = 0
    first = True
    while act.any():
        idx = np.nonzero(act)[0]
        st = np.where(rng.random(idx.size) < 0.5, 1, -1)
        st[pos[idx] == 0] = 1
        st[pos[idx] == L] = -1
        pos[idx] += st
        T[idx] += 1
        back = idx[pos[idx] == 0]
        act[back] = False
    P(f"      L = {L:4d}: ⟨T⟩/(2L) = {T.mean()/(2*L):.4f} ± {T.std()/math.sqrt(T.size)/(2*L):.4f}   (Kac: 1; in Zeit 2L·h/c = 2Λ/c)")

# Wegfall von D: q hängt nur von Γ·⟨T⟩ ab, wenn die Rückkehrtiefe die Diffusionstiefe einer Kac-Zeit ist
P("")
P("   Verlust je Schritt = Taktgefälle über einen Schritt: p = h·|∇ln τ̇| = h·a_N/c²  ⇒  Rate p/(h/c) = a_N/c = Γ, für jedes h.")
P("   Tiefe ϖ_R² = D·⟨T⟩ mit D = hc/2, ⟨T⟩ = 2Λ/c ⇒ m = ϖ_R/h = √(Λ/h). Gitter exakt, Milchstraße 6·10¹⁰ M☉ bei 8 kpc:")
M = 6e10 * K.M_sun
r = 8e3 * K.pc
aN = K.G * M / r**2
Lam = K.schliesslaenge()
for h in (K.r0(M), 10 * K.r0(M), 1e3 * K.r0(M), 1e6 * K.r0(M)):
    p = h * aN / K.c**2
    m = math.sqrt(Lam / h)
    P(f"      h = {h/K.r0(M):8.0e}·r₀:  m = {m:11.4e},  p = {p:.4e},  q = {q_gitter(m, p):.6f}")
P(f"      Regel e^(−s) = {float(K.q_(aN)):.6f}   (s = {float(K.s_(aN)):.6f})")
P("   ⇒ Die Schrittlänge fällt heraus: q, Gewicht und Kern hängen nur an Γ = c|∇ln τ̇| und ⟨T⟩ = 2Λ/c.")
P("     (3) wird nur gebraucht, um die Tiefe als Länge zu lesen (ϖ_R = r* bei h = r₀) und Γ = D k² mit k = 1/r zu schreiben:")
r0 = K.r0(M)
P(f"     Verlust je Schritt r₀ aus dem Takt p = r₀·a_N/c² = {r0*aN/K.c**2:.6e}; Zerfall der Bahnmode je Schritt 1 − cos(r₀/r) = {1-math.cos(r0/r):.6e}")
tk = lambda rr: math.sqrt(1 - r0 / rr)
p_exakt = math.log(tk(r + r0 / 2) / tk(r - r0 / 2))
P(f"     mit dem vollen Takt √(1 − r₀/r) über einen Schritt: p = {p_exakt:.6e}  (gleich: das Taktgefälle über r₀ ist die Glättung von k = 1/r über r₀)")

# (4) eine Uhr: Verluste je Takt addieren sich
kopf("2  Eine Uhr (4) und ein Sinn (6): Monte-Carlo auf dem Gitter")
m = 25
nw = 40000
P(f"   Tiefe m = {m} Schritte, geschlossen bei L = 40·m; Verlust je Takt p_Γ (Takt) und p_θ (Ausdehnung) auf derselben Uhr")
for sG, sT in ((1.0, 0.0), (0.0, 2.0), (1.0, 2.0), (0.5, 1.0)):
    pG = sG**2 / (2 * m * m)
    pT = sT**2 / (2 * m * m)
    qm, qe, _ = mc_rueckkehr(m, pG, nw, L=40 * m, theta_p=pT)
    P(f"      Γτ_D = {sG**2:4.2f}, θτ_D = {sT**2:4.2f}:  MC {qm:.4f} ± {qe:.4f}   exp(−√Σ) {math.exp(-math.sqrt(sG**2+sT**2)):.4f}"
      f"   Gitter {q_gitter(m, pG+pT):.4f}   Produkt {math.exp(-sG-sT):.4f}")
P("   ⇒ (4) folgt: Ein Takt für beide Verluste, sie addieren sich je Schritt und stehen unter einer Wurzel.")
P("")
P("   Zahl der Rückläufe (nach jedem Rücklauf Erneuerung aus der Tiefe) bis zum Verlust, + 1 für den Ursprung:")
for s in (0.5, 1.0, 2.0):
    p = s**2 / (2 * m * m)
    cm, ce, _ = mc_rueckkehr(m, p, nw, L=40 * m, renew=True)
    q = q_gitter(m, p)
    P(f"      s = {s:3.1f}:  MC 1 + ⟨n⟩ = {1+cm:.4f} ± {ce:.4f}   1/(1−q) = {1/(1-q):.4f}   (1+q)/(1−q) = {(1+q)/(1-q):.4f}   Regel ν = {1/(1-math.exp(-s)):.4f}")
P("   ⇒ (6) folgt: Rückläufe sind Ereignisse einer Folge; gezählt wird jedes einmal (Erneuerungszahl 1/(1−q)).")

# (2) bleibt: Ohne Tiefe (Rückkehr des ersten Schritts) keine Schwelle
kopf("3  Warum die Tiefe (2) nicht aus der Irrfahrt folgt")
P("   Rückkehr vom Ort aus (erster Schritt, dann erster Rücklauf), Gitter exakt: Anteil 1 − √(1 − (1−p)²) ⇒ ν = 1/√(1−(1−p)²) ≈ r/r₀")
for nm, MM, rr in (("Milchstraße 8 kpc", M, r), ("Sonne 1000 AE", K.M_sun, 1000 * K.AU)):
    pp = K.r0(MM) * K.G * MM / rr**2 / K.c**2
    nu_ex = 1 / math.sqrt(pp * (2 - pp))
    P(f"      {nm:18s}: ν = {nu_ex:.4e}   r/r₀ = {rr/K.r0(MM):.4e}   Regel ν = {float(K.nu(K.G*MM/rr**2)):.4f}")
P("   ⇒ a = c²/(2r), ohne a₀: Die Irrfahrt allein kennt nur r₀. Die Schwelle kommt erst mit der Tiefe ϖ_R² = D·⟨T⟩,")
P("     d. h. mit der Kac-Zeit der bei Λ geschlossenen Halbgeraden. Kac liefert ⟨T⟩ = 2Λ/c und die Massenfreiheit;")
P("     dass die Erneuerung aus der in einer Kac-Zeit erreichten Tiefe kommt, bleibt gesetzt.")
P("   Andere naheliegende Tiefen aus derselben Irrfahrt:")
P(f"      mittlere Durchgangszeit ϖ→0 auf [0,Λ] gleich ⟨T⟩: ϖ(2Λ − ϖ)/(2D) = 2Λ/c ⇒ ϖ ≈ 2D/c = h (ein Schritt, wie oben)")
P(f"      Anteil der Ausflüge bis ϖ (h/ϖ) gleich Anteil der Linie (ϖ/Λ) ⇒ ϖ = √(hΛ) = r* (bei h = r₀): dieselbe Gleichung, kein neuer Grund")

# (5) lokal: Raum als Verlust gegen Raum als Transport
kopf("4  Warum lokal Γ/D (5): Raumschritt als Taktvergleich gegen Raumschritt als Transport (Punktmasse)")
P("   Einheiten r* = 1, τ_D = 1, D = 1 (ϖ_R = 1). Verlustrate am Ort Γτ_D = a_N/a₀ = 1/ρ², ρ = |x|/r*.")
P("   (R) Rückkehr in ϖ am festen Ort, Verlust = Taktgefälle am Ort;   (T) Rückkehr läuft zugleich im Raum (3D) und verliert")
P("   mit dem Takt am jeweiligen Ort, ∂_t u = D(∇² + ∂_ϖ²)u − Γ(x)u. Monte-Carlo, gleicher Code, (R) = Raumdiffusion aus.")


def mc_4d(rho0, n, raum, Lam=20.0, eps=2e-3):
    X = np.zeros((n, 3))
    X[:, 0] = rho0
    w = np.ones(n)
    z = np.ones(n)
    done = np.zeros(n, bool)
    res = np.zeros(n)
    for _ in range(3_000_000):
        idx = np.nonzero(~done)[0]
        if idx.size == 0:
            break
        rho = np.linalg.norm(X[idx], axis=1)
        dt = np.clip(0.01 * np.minimum(rho, z[idx])**2, 1e-7, 0.02)
        sq = np.sqrt(2 * dt)
        w[idx] *= np.exp(-dt / rho**2)
        if raum:
            X[idx] += sq[:, None] * rng.standard_normal((idx.size, 3))
        z[idx] += sq * rng.standard_normal(idx.size)
        z[idx] = np.where(z[idx] > Lam, 2 * Lam - z[idx], z[idx])
        hit = idx[z[idx] <= eps]
        res[hit] = w[hit]
        done[hit] = True
        low = idx[w[idx] < 1e-6]
        done[low] = True
    return res.mean(), res.std() / math.sqrt(n)


P("      ρ = r/r*   a_N/a₀  |  q (R) MC        e^(−1/ρ)  |  q (T) MC        |  ν (R)    ν (T)")
rows = []
for rho0 in (0.05, 0.1, 0.3, 1.0, 3.0):
    nn = 3000
    qR, eR = mc_4d(rho0, nn, False)
    qT, eT = mc_4d(rho0, nn, True)
    rows.append((rho0, qT))
    P(f"      {rho0:8.2f}  {1/rho0**2:8.2f}  |  {qR:.4f} ± {eR:.4f}   {math.exp(-1/rho0):.4f}  |  {qT:.4f} ± {eT:.4f}  |  {1/(1-qR):.4f}   {1/(1-qT):.4f}")
lr = np.log([x[0] for x in rows[:3]])
lq = np.log([max(x[1], 1e-12) for x in rows[:3]])
P(f"   (T): q ∝ ρ^b bei kleinem ρ, b = {np.polyfit(lr, lq, 1)[0]:.3f};  radial u'' + 2u'/ρ − u/ρ² = 0 gibt b = (√5 − 1)/2 = {(math.sqrt(5)-1)/2:.3f}")
b = (math.sqrt(5) - 1) / 2
rho_E = K.AU / K.r_stern(K.M_sun)
P(f"   Erdbahn: ρ = 1 AE/r*☉ = {rho_E:.3e};  (T) ν − 1 ~ q(ρ=1)·ρ^b = {rows[3][1]*rho_E**b:.1e}  gegen Regel {float(K.nu(K.G*K.M_sun/K.AU**2))-1:.1e}")
P("   ⇒ Läuft die Rückkehr im Raum mit (Transport), wandert sie in einer Kac-Zeit um r* vom Ort weg, liest dort kleineres Gefälle")
P("     und verliert die Newton-Grenze (Potenz statt Exponent). Ohne eigenen Verlust (nur ∇²) ist es die nichtlokale Lesart")
P("     mit v⁴ = (2/π)²GMa₀. Die Regel entsteht nur, wenn der Raumschritt nichts trägt, sondern den Takt über einen Schritt")
P("     vergleicht: Verlust p = r₀|∇ln τ̇| am Ort. Da r₀/r* = 2δ* ≪ 1 (Milchstraße 6,6·10⁻⁷), ist dieser Vergleich lokal.")

# SPARC
kopf("5  Gezählte Rückläufe gegen ν(g_N): Punktmasse und SPARC-Galaxien (Gitter exakt und Monte-Carlo)")
P("   Verlust je Schritt aus dem Taktgefälle des baryonischen Felds, p = h·g_N/c², Tiefe m = √(Λ/h); MC mit h so, dass m = 30.")
P("   Punktmasse Milchstraße 6·10¹⁰ M☉:")
for rk in (2, 8, 20, 50):
    rr = rk * 1e3 * K.pc
    g = K.G * M / rr**2
    hh = K.r0(M)
    q_ex = q_gitter(math.sqrt(Lam / hh), hh * g / K.c**2)
    mm = 30
    p_mc = (mm**2) ** -1 * 0.5 * g / K.A0 * 1.0  # h = Λ/m² ⇒ p = h g/c² = g/(2 a₀ m²)
    qm, qe, _ = mc_rueckkehr(mm, p_mc, 20000, L=40 * mm)
    P(f"      r = {rk:3d} kpc: ν Gitter(r₀) {1/(1-q_ex):.4f}   ν MC {1/(1-qm):.4f} (q {qm:.4f} ± {qe:.4f})   Regel ν {float(K.nu(g)):.4f}")

sp = os.path.join(os.path.dirname(__file__), "..", "daten", "sparc", "Rotmod_LTG")
for nm in ("NGC3198", "DDO154"):
    d = np.loadtxt(os.path.join(sp, nm + "_rotmod.dat"), comments="#")
    R, Vo, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    Vbar2 = np.abs(Vg) * Vg + 0.5 * Vd**2 + 0.7 * Vb**2
    gN = np.maximum(Vbar2, 1e-6) * 1e6 / (R * 1e3 * K.pc)
    Mt = Vbar2[-1] * 1e6 * R[-1] * 1e3 * K.pc / K.G   # Masse für r₀ (fällt heraus)
    hh = K.r0(Mt)
    nu_g = np.array([1 / (1 - q_gitter(math.sqrt(Lam / hh), hh * g / K.c**2)) for g in gN])
    nu_r = K.nu(gN)
    v_g = np.sqrt(nu_g * gN * R * 1e3 * K.pc) / 1e3
    v_r = np.sqrt(nu_r * gN * R * 1e3 * K.pc) / 1e3
    chi_g = np.mean(((v_g - Vo) / eV)**2)
    chi_r = np.mean(((v_r - Vo) / eV)**2)
    P(f"   {nm} (Υ_d = 0,5, Υ_b = 0,7, N = {len(R)}): max|ν_Gitter/ν_Regel − 1| = {np.max(np.abs(nu_g/nu_r-1)):.1e};  χ²/N Gitter {chi_g:.2f}, Regel {chi_r:.2f}")
    sel = np.linspace(0, len(R) - 1, 4).astype(int)
    for i in sel:
        mm = 30
        qm, qe, _ = mc_rueckkehr(mm, 0.5 * gN[i] / K.A0 / mm**2, 80000, L=40 * mm)
        P(f"      R = {R[i]:6.2f} kpc, g_N/a₀ = {gN[i]/K.A0:7.3f}: ν MC {1/(1-qm):.4f} (± {qe/(1-qm)**2:.4f})   Gitter {nu_g[i]:.4f}   Regel {nu_r[i]:.4f}"
              f"   v {v_g[i]:6.1f} (gemessen {Vo[i]:6.1f} ± {eV[i]:4.1f}) km/s")
P("   ⇒ Die gezählten Rückläufe der Irrfahrt mit Taktverlust geben die Regel; MC-Abweichungen nur statistisch und Gitter m = 30.")

kopf("6  Handprüfung (anderer Weg)")
P("   Milchstraße 6·10¹⁰ M☉, r = 8 kpc = 2,4686e20 m: GM = 1,32712e20·6e10 = 7,9627e30 m³/s²")
P("   a_N = 7,9627e30/(2,4686e20)² = 7,9627e30/6,0940e40 = 1,3067e-10 m/s²;  s = √(1,3067/1,1) = √1,18791 = 1,08991")
P("   r₀ = 2GM/c² = 1,5925e31/8,98755e16 = 1,7719e14 m;  r* = √(7,9627e30/1,1e-10) = 2,6905e20 m;  m = r*/r₀ = 1,5184e6")
P("   Verlust je Schritt p = r₀a_N/c² = (r₀/r)²/2 = (7,1778e-7)²/2 = 2,5760e-13;  √(2p)·m = 7,1778e-7·1,5184e6 = 1,0899")
P("   q = e^{−1,0899} = 0,33625;  ν = 1/(1 − 0,33625) = 1,5066;  Γ = p/(r₀/c) = 2,5760e-13·2,99792e8/1,7719e14 = 4,358e-19 s⁻¹")
P(f"   Code: s = {float(K.s_(aN)):.5f}, ν = {float(K.nu(aN)):.4f}, Γ = {K.verlust_gamma(aN):.4e}")

with open(os.path.join(os.path.dirname(__file__), "setzungen_gleichung.txt"), "w") as f:
    f.write("\n".join(out) + "\n")
