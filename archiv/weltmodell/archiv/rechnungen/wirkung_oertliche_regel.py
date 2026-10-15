# Bruchstelle 2: Wirkung für die örtliche Regel F = ν(|g_N|)·g_N = −ν∇Φ_N.
# F ist flächennormal zu den Newton-Äquipotentialen (F·∇×F = 0), 1/ν integrierender Faktor.
# Ruhendes Gas: ∇P = ρF ⇒ P = P(Φ_N), ρν = −dP/dΦ_N ⇒ kT ∝ P/ρ ∝ ν längs jeder Isobare.
# Abgeplattete Galaxie: Miyamoto-Nagai-Potential, Einheiten G = M = a0 = 1 (r* = 1).
import numpy as np
def nu(g): s = np.sqrt(g); return 1/(-np.expm1(-s))
def gN(R, z, a, b):
    zeta = np.sqrt(z*z + b*b); D = np.sqrt(R*R + (a + zeta)**2)
    return np.hypot(R/D**3, z*(a + zeta)/(zeta*D**3))
# Probe von Hand: Punktmasse (a = b = 0): |g| = 1/r², ν bei r = r* = 1/(1 − e⁻¹) = 1,5820
print("Probe Punktmasse r = 1:", round(nu(gN(1.0, 0.0, 0.0, 1e-9)), 4))
print("MN (a, b) | R_eq | a_N/a0 Ebene / Achse | T_Achse/T_Ebene längs der Isobare | ν(g)·g-Wirbel: ∮F·dl / ∮|F|dl")
for a, b in ((0.25, 0.05), (0.25, 0.15), (0.25, 0.25)):
    for Req in (0.3, 1.0, 3.0):
        D0 = np.sqrt(Req**2 + (a + b)**2)
        zmax = np.sqrt((D0 - a)**2 - b**2)
        zs = np.linspace(0, zmax, 2001)
        Rs = np.sqrt(np.maximum(D0**2 - (a + np.sqrt(zs**2 + b*b))**2, 0))
        g = gN(Rs, zs, a, b); v = nu(g)
        # Arbeit der Regel auf einem geschlossenen Weg: Ebene R_eq → 1,5 R_eq, Äquipotential hinauf zur Achse, Achse hinab
        Req2 = 1.5*Req; D1 = np.sqrt(Req2**2 + (a + b)**2)
        def path_work(D):
            zm = np.sqrt((D - a)**2 - b**2); z = np.linspace(0, zm, 4001); R = np.sqrt(np.maximum(D**2 - (a + np.sqrt(z**2 + b*b))**2, 0))
            return R, z
        # auf Äquipotentialen ist F·dl = 0; Arbeit nur auf den Radialstücken: Ebene (z = 0) und Achse (R = 0)
        Rp = np.linspace(Req, Req2, 4001); Wp = np.trapezoid(nu(gN(Rp, 0*Rp, a, b))*gN(Rp, 0*Rp, a, b), Rp)
        za0, za1 = np.sqrt((D0 - a)**2 - b**2), np.sqrt((D1 - a)**2 - b**2)
        za = np.linspace(za0, za1, 4001); Wa = np.trapezoid(nu(gN(0*za, za, a, b))*gN(0*za, za, a, b), za)
        dPhi = 1/D0 - 1/D1
        print(f"({a}, {b}) | {Req} | {g[0]:.3f} / {g[-1]:.3f} | {v[-1]/v[0]:.3f} (größtes/kleinstes {v.max()/v.min():.3f}) | "
              f"Umlaufarbeit {(Wp - Wa):+.4f} bei Newton-Stufe {dPhi:.4f}: {(Wp - Wa)/(0.5*(Wp + Wa)):+.3f}")


# Dünne Exponentialscheibe (M = 1, G = a0 = 1): Der Wirbel der Regel sitzt in der Schicht.
# In der Ebene |g_N| = g_R (Freeman), knapp darüber √(g_R² + (2πGΣ)²). Umlauf durch die Schicht zwischen R1 und R2:
# ∮F·dl = ∫[ν(g_R) − ν(√(g_R² + (2πGΣ)²))]·g_R dR  (senkrechte Stücke durch die dünne Schicht heben sich auf)
from scipy.special import i0, i1, k0, k1
def gR_freeman(R, Rd):
    y = R/(2*Rd); return (2/Rd)*y*y*(i0(y)*k0(y) - i1(y)*k1(y))/R
print("Probe Freeman: v² bei 2,2 R_d (R_d = 1) =", round(2.2*gR_freeman(2.2, 1.0), 5), "(Handwert 0,3864)")
print("R_d | Abschnitt | 2πGΣ/g_R | Umlaufarbeit / ∫ν(g_R)g_R dR")
for Rd in (0.3, 1.0, 3.0):
    for (r1, r2) in ((0.25, 1.0), (1.0, 3.0), (3.0, 6.0)):
        R = np.linspace(r1*Rd, r2*Rd, 4001); gR = gR_freeman(R, Rd); gz = 2*np.pi*np.exp(-R/Rd)/(2*np.pi*Rd*Rd)
        W = np.trapezoid((nu(gR) - nu(np.hypot(gR, gz)))*gR, R); W0 = np.trapezoid(nu(gR)*gR, R)
        print(f"{Rd} | {r1}–{r2} R_d | {gz[0]/gR[0]:.2f}–{gz[-1]/gR[-1]:.2f} | {W/W0:.3f}")
