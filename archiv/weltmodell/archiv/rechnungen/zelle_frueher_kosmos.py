# Bruchstelle 1: Welche Zelle setzt a0 = c²/U, geprüft am ganzen Kosmos.
# Vakuumteil der Strukturwirkung: Ω_s = k_s·Θ², Θ = 2πa0/(cH), k_s = π²/45, w_s = −1 − (2/3)·d ln a0/d ln a.
# Lesarten: A Zelle atmet mit dem Horizont (U = 2πc/H, a0 ∝ H)
#           B Zelle mitbewegt wie die Schleifen des Autors (U ∝ a, a0 ∝ 1+z)
#           C Schließlänge fest (a0 konstant); C∞: fest am de-Sitter-Endhorizont
import numpy as np
from scipy.integrate import quad
ks = np.pi**2/45; Om = 0.315; Or = 9.1e-5; H0 = 67.4; c = 2.998e5
a0H = 2.998e8*H0*1e3/3.0857e22/(2*np.pi)          # cH0/2π in m/s²
print(f"cH0/2π = {a0H:.4e} m/s², k_s = {ks:.4f}")

def E2(z, L, Th0=1.0):
    Rm = Om*(1+z)**3 + Or*(1+z)**4
    if L == "A":  return (Rm + (1-Om-Or-ks)*1.0)/(1-ks)          # Θ = 1 immer, R(0) = 1 − k_s
    if L == "B":  return Rm + (1-Om-Or-ks*Th0**2) + ks*Th0**2*(1+z)**2
    if L == "C":  return Rm + (1-Om-Or)                            # Λ_s konstant, Summe = Ω_DE
    if L == "L":  return Rm + (1-Om-Or)                            # ΛCDM
def Oms(z, L, Th0=1.0):
    if L == "A": return ks
    if L == "B": return ks*Th0**2*(1+z)**2/E2(z, L, Th0)
    if L == "C": return ks*Th0**2/E2(z, L, Th0)
def a0rel(z, L, Th0=1.0):
    return {"A": np.sqrt(E2(z, "A")), "B": 1+z, "C": 1.0}[L]
def DM(z, L, Th0=1.0): return c/H0*quad(lambda x: 1/np.sqrt(E2(x, L, Th0)), 0, z, limit=200)[0]
def dNeff(Os):  # konstanter Anteil f am Gesamtbudget im Strahlungszeitalter → ΔN_eff
    return Os/(1-Os)*(1+0.2271*3.046)/0.2271

Einf = np.sqrt(1-Om-Or)                                    # C∞: Θ∞ = 1 ⇒ Θ0 = E∞
print(f"C∞: E∞ = H∞/H0 = {Einf:.4f}, a0 = cH∞/2π = {Einf*a0H:.4e} m/s², U = 2πc/H∞ = {2*np.pi*c/(H0*Einf)/1e3:.1f} Gpc")
DL = DM(1090, "L")
for L, Th0, lab in (("A", 1, "A Horizont"), ("B", 1, "B mitbewegt"), ("C", 1.0, "C fest, Θ0 = 1"), ("C", Einf, "C∞ de Sitter")):
    zs = (1e9, 1090, 2, 1, 0)
    Os = [Oms(z, L, Th0) for z in zs]
    # w_s heute numerisch aus a0(z)
    h = 1e-4; dlna0 = (np.log(a0rel(h, L, Th0)) - np.log(a0rel(0, L, Th0)))/(-h)   # d ln a0/d ln a = −(1+z) d ln a0/dz
    ws = -1 - (2/3)*dlna0
    print(f"{lab:16s} Ω_s(10⁹) = {Os[0]:.2e} (ΔN_eff {dNeff(Os[0]):.2f}), Ω_s(1090) = {Os[1]:.2e}, Ω_s(0) = {Os[4]:.3f}, "
          f"w_s(0) = {ws:+.2f}, a0(1)/a0(0) = {a0rel(1, L, Th0):.2f}, a0(2)/a0(0) = {a0rel(2, L, Th0):.2f}, "
          f"D_M(1090)/ΛCDM = {DM(1090, L, Th0)/DL:.4f}")
