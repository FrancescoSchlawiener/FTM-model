"""Kernfrage τ: Welche Rolle spielt Re τ = ω_rad/ω_azim = κ/Ω im Nome q = e^{2πiτ}?
Modul-Symmetrie T: τ → τ+1 lässt q unverändert. Kepler: κ/Ω = 1 ≡ 0 (mod 1)."""
import math, cmath
G = 6.674e-11; MS = 1.989e30; KPC = 3.0857e19; a0 = 1.2e-10
M = 6e10*MS; rs = math.sqrt(G*M/a0)

def nu_rar(y): return 1/(1-math.exp(-math.sqrt(y)))
def v_of(r):
    gN = G*M/r**2; return math.sqrt(nu_rar(gN/a0)*gN*r)
def kappa_over_omega(r, h=1e-4):
    dlnv = (math.log(v_of(r*(1+h))) - math.log(v_of(r*(1-h))))/(2*h)
    return math.sqrt(2*(1+dlnv))

print(f"Milchstraße (Punktmasse 6e10 M☉), r* = {rs/KPC:.2f} kpc\n")
print(f"{'r [kpc]':>8s} {'y=gN/a0':>9s} {'κ/Ω':>7s} {'Re τ mod 1':>10s} {'ν (|q|)':>9s} {'Re ν(q komplex)':>16s} {'Δ [dex]':>8s} {'|ν komplex|':>12s}")
worst = 0
for rk in [0.5, 1, 2, 5, 8.35, 15, 30, 60, 120]:
    r = rk*KPC; gN = G*M/r**2; y = gN/a0
    ko = kappa_over_omega(r); ph = ko % 1.0
    qa = math.exp(-math.sqrt(y))
    qc = qa*cmath.exp(2j*math.pi*ko)
    nu_r = 1/(1-qa); nu_c = (1/(1-qc)).real; nu_abs = abs(1/(1-qc))
    d = math.log10(nu_c/nu_r); worst = max(worst, abs(d))
    print(f"{rk:8.2f} {y:9.3e} {ko:7.4f} {ph:10.4f} {nu_r:9.4f} {nu_c:16.4f} {d:+8.3f} {nu_abs:12.4f}")
print(f"\nGrößte Abweichung mit Phase: {worst:.3f} dex  (Streuung der gemessenen RAR ≈ 0,13 dex, systematisch < 0,05 dex)")
print("Im Kepler-Bereich ist κ/Ω = 1 ≡ 0 (mod 1): die Phase verschwindet exakt (T-Symmetrie).")
print("Im flachen Bereich ist κ/Ω → √2, Re τ mod 1 → 0,414: die Phase dreht q aus der reellen Achse.")
