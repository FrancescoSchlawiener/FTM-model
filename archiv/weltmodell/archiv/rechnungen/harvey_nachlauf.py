"""Harvey et al. 2015 (72 Stoßpaare): Messgröße β = δ_SI/δ_SG gegen den Nachlauf der Mode ξ.

Messgröße (Harvey 2015, aus Suchauszügen, Volltext nicht erreichbar):
  β = Versatz(Linse − Sterne) / Versatz(Gas − Sterne), längs der Stoßachse, ⟨β⟩ = −0,04 ± 0,07 (72 Unterstrukturen, 30 Systeme).
Abbildung im Modell (Schwerpunkt der gelesenen Linse, Sterne bei 0, Gas bei 1):
  gelesene Linse = Baryonen (Newton, am wahren Ort)
                 + Operatorform der Regel (ν − 1) auf die Baryonen, gebildet an der Gedächtnislage:
                   Sterne bei 0, Gas bei φ(x),  φ(x) = 1 − (1 − e^{−x})/x,  x = λ_m·t  (sofortige Abbremsung)
                 + nicht übergegangener Rest der Mode ξ: zieht mit den stoßfreien Bahnen (Sterne, 0).
  β(x) = A + B·φ(x),  A = M_Gas/M_Linse (Newton-Gas; genauso in ΛCDM, weil das Gas dort ebenso linst),
                      B = Anteil der Operatorform am Gas an der Linse = f_Gas·(1 − M_b/M_L)·f_Op.
  Gemessen wird gegen ΛCDM der Überschuss Δβ = B·φ(x).
Rechnet nur mit ftm_kern.
"""
import math
from ftm_kern import uebergabezeit, Gyr, Q_KOSM

def phi_s(x):
    return 1 - (-math.expm1(-x)) / x if x > 1e-9 else x / 2

def phi_inv(p):
    lo, hi = 1e-9, 1e4
    for _ in range(200):
        m = math.sqrt(lo * hi)
        lo, hi = (m, hi) if phi_s(m) < p else (lo, m)
    return math.sqrt(lo * hi)

print("Harvey 2015 gegen den Nachlauf der Mode ξ (Schwerpunkt-Abbildung)\n")

# --- Messung [G: Suchauszüge, Volltext gesperrt] ---
BETA, DBETA = -0.04, 0.07
print(f"Messung: ⟨β⟩ = {BETA} ± {DBETA} (72 Unterstrukturen, 30 Systeme) [G]")
for k in (1, 2):
    print(f"   Obergrenze {k}σ: β ≤ {BETA + k*DBETA:+.2f}")
DB_MAX = BETA + 2 * DBETA
print(f"   verwendet: Überschuss Δβ ≤ {DB_MAX:.2f} (2σ, ΛCDM-Grundwert A abgezogen)\n")

# --- Anteile ---
f_gas = 10.2 / 11.2            # Gas/Sterne 10,2 (Bullet-Rechnung linsen.py); Haufen allgemein 0,85–0,92
ML_Mb = 1 + Q_KOSM            # Linse/Baryonen ≈ 1 + Q (Modell trifft die Linsenmasse auf 1–5 %)
extra = 1 - 1 / ML_Mb
print(f"Anteile: Gas an den Baryonen f_Gas = {f_gas:.3f}, M_L/M_b = 1 + Q = {ML_Mb:.2f}, Zusatz an der Linse = {extra:.3f}")
A = f_gas / ML_Mb
print(f"   A = M_Gas/M_L = {A:.3f} (Newton-Gas; gleich in ΛCDM)")
F_OP = {"100 kpc (Op 20 %)": 0.20, "Bullet (Op 20 % der gelesenen Masse)": None, "1,5 Mpc (Op 68 %)": 0.68}
Bs = {}
for name, fop in F_OP.items():
    if fop is None:
        B = 0.20 * extra      # 03: 20 % der gelesenen Masse folgen dem Gas mit Nachlauf
    else:
        B = f_gas * extra * fop
    Bs[name] = B
    print(f"   B[{name:38s}] = {B:.3f}")
print()

# --- Übergaberate λ_m = ω_orb/N* ---
print("1/ω_orb und 1/λ_m = N*/ω_orb  [Gyr]")
print(f"   {'Δ·ρ_krit':>9s} {'z':>4s} {'1/ω_orb':>8s} {'N*=0,3':>7s} {'N*=2,4':>7s}")
TAB = {}
for D in (200, 500):
    for z in (0.0, 0.3, 0.5):
        w1, lmin, lmax = uebergabezeit(D, z)
        TAB[(D, z)] = w1 / Gyr
        print(f"   {D:9d} {z:4.1f} {w1/Gyr:8.3f} {lmin/Gyr:7.3f} {lmax/Gyr:7.3f}")
print()

# --- β(x) ---
print("Überschuss Δβ(x) = B·φ(x) und gesamtes β = A + B·φ(x)")
print(f"   {'x=λt':>6s} {'φ':>6s} " + " ".join(f"{'Δβ B='+format(b,'.2f'):>11s}" for b in Bs.values()))
for x in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0):
    p = phi_s(x)
    print(f"   {x:6.2f} {p:6.3f} " + " ".join(f"{b*p:11.3f}" for b in Bs.values()))
print()

# --- Schranke ---
print(f"Schranke Δβ ≤ {DB_MAX:.2f}: φ ≤ Δβ/B,  x_max = φ⁻¹, 1/λ_m ≥ t/x_max, N* ≥ ω_orb·t/x_max")
TS = (0.2, 0.5, 1.0)
print(f"   {'B':>6s} {'φ_max':>6s} {'x_max':>7s} " + " ".join(f"{'t='+str(t)+' Gyr: 1/λ≥':>18s}" for t in TS)
      + "   N* nötig (Δ=200…500, z=0,3; t=0,5 Gyr)")
for name, B in Bs.items():
    pm = DB_MAX / B
    if pm >= 1:
        print(f"   {B:6.3f}   ≥ 1     –   jede Nachlaufzeit verträglich")
        continue
    xm = phi_inv(pm)
    lam = [t / xm for t in TS]
    ns = [TAB[(D, 0.3)] ** -1 * 0.5 / xm for D in (200, 500)]
    print(f"   {B:6.3f} {pm:6.3f} {xm:7.3f} " + " ".join(f"{l:15.2f} Gyr" for l in lam)
          + f"   N* ≥ {ns[0]:.2f} … {ns[1]:.2f}")
print()

# --- Modellband N* = 0,3…2,4 gegen die Schranke, je t ---
print("Modellband: Δβ für N* = 0,3 / 1 / 2,4 (Δ = 200 bzw. 500, z = 0,3), B = Bullet-Wert")
Bb = Bs["Bullet (Op 20 % der gelesenen Masse)"]
print(f"   {'t[Gyr]':>6s} " + " ".join(f"{'N*='+str(n)+' Δ='+str(D):>14s}" for n in (0.3, 1.0, 2.4) for D in (200, 500)))
for t in (0.15, 0.3, 0.5, 1.0, 2.0):
    row = []
    for n in (0.3, 1.0, 2.4):
        for D in (200, 500):
            x = t / (n * TAB[(D, 0.3)])
            row.append(Bb * phi_s(x))
    print(f"   {t:6.2f} " + " ".join(f"{v:14.3f}" for v in row))
print()

# --- Handprüfung ---
print("Handprüfung:")
print(f"   φ(1) = 1 − (1 − e⁻¹) = e⁻¹ = {math.exp(-1):.4f}; Rechnung {phi_s(1.0):.4f}")
print(f"   1/ω_orb(Δ=200, z=0) = 1/(10·H₀) = 977,8/67,4/10 Gyr = {977.8/67.4/10:.3f}; Rechnung {TAB[(200,0.0)]:.3f}")
print(f"   B_Bullet = 0,20·(1 − 1/6,36) = {0.20*(1-1/6.364):.3f}; Rechnung {Bb:.3f}")
print("\nFertig.")
