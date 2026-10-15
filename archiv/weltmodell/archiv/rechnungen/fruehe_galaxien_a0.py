# -*- coding: utf-8 -*-
"""
Frühe Scheibengalaxien (z ≈ 0,9–2,5) und die Schwelle: festes a₀ (Modell) gegen a₀ ∝ H(z) (verworfene Lesart).
1 Anteil f_DM(R_e) = 1 − a_N/a = q (örtliche Regel) für Exponentialscheiben (+ Hernquist-Bulge),
  bei festem a₀ (Kern a0_frueh_relativ) und bei a₀(z) = A0·E(z).
2 Tully-Fisher-Nullpunkt: tiefer Grenzfall v⁴ = GM_b a₀ und am Ort v(2,2 R_d) derselben Scheibe.
3 Druckstütze (Burkert 2010 [G]): v_c² = v_rot² + 2σ₀²·R/R_d; wie stark f_DM ohne Korrektur verrutscht.
4 Preisach: Abweichung frisch kollabierter Galaxien (Werte aus gewicht_preisach.txt, Abschnitt 3)
  in dex von M_b und v; Stichprobe für 3σ bei 0,1–0,2 dex Fehler je Galaxie.
5 Handprüfung.
Kern: A0, G, M_sun, kpc, nu, chi, q_, E_hintergrund, a0_frueh_relativ, hernquist_masse, v_flach.
Nicht im Kern (Newton, hier gerechnet): Freeman-Scheibe v²_N = 4πGΣ₀R_d y²[I₀K₀ − I₁K₁], y = R/2R_d.
Galaxienwerte (M_b, R_e, σ₀) typisch nach Genzel 2017/2020, Nestor Shachar 2023 [G] (Quellen nicht geholt).
"""
import math, os, sys
import numpy as np
from scipy.special import i0, i1, k0, k1
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ftm_kern import *

OUT = []
def P(s=""):
    OUT.append(s); print(s)


def aN_scheibe(M, Rd, R, fb=0.0, ab=None):
    """Newton-Beschleunigung in der Ebene: Freeman-Scheibe (1−fb)M + Hernquist-Bulge fb·M."""
    Md = (1 - fb) * M
    y = R / (2 * Rd)
    S0 = Md / (2 * math.pi * Rd**2)
    v2 = 4 * math.pi * G * S0 * Rd * y**2 * (i0(y) * k0(y) - i1(y) * k1(y))
    a = v2 / R
    if fb > 0:
        a += G * hernquist_masse(fb * M, ab, R) / R**2
    return a


def a0_fest(z):
    return A0 * a0_frueh_relativ(z)


def a0_H(z):
    return A0 * E_hintergrund(z)


def fDM(aN, a0):
    return 1 - 1 / float(nu(aN, a0))


# ------------------------------------------------------------------ 0
P("0  Eingaben: A0 = %.2e m/s², H(z) aus E_hintergrund (Ω_m = %.4f)" % (A0, OMEGA_M))
for z in (0.9, 1.5, 2.0, 2.3, 2.5):
    P("   z = %.1f: E(z) = %.4f, a0_fest/A0 = %.5f, a0_H/A0 = %.4f" % (z, E_hintergrund(z), a0_frueh_relativ(z), E_hintergrund(z)))
P("   R_e = 1,678 R_d (Exponentialscheibe); Gas in der Scheibe; Bulge Hernquist a_b = R_e,b/1,815 mit R_e,b = 1 kpc [G]")

# ------------------------------------------------------------------ 1
P()
P("1  f_DM(R_e) = 1 − a_N/a = q = e^{−√(a_N/a₀)} (örtliche Regel, Scheibe in der Ebene)")
P("   Typische Scheiben [G]: M_b, R_e, Bulgeanteil B/T")
gal = [("z≈0,9", 0.9, 7e10, 6.0, 0.1), ("z≈1,5", 1.5, 1.0e11, 5.5, 0.2),
       ("z≈2,3", 2.3, 1.2e11, 4.5, 0.3), ("z≈2,3 groß", 2.3, 1.2e11, 7.0, 0.2)]
P("   Galaxie        M_b[M☉]   R_e[kpc]  B/T   a_N(R_e)/A0   f_DM fest   f_DM a0∝H   Σ_b(<R_e)[M☉/pc²]")
for name, z, M, Re, fb in gal:
    Rd = Re / 1.678
    aN = aN_scheibe(M * M_sun, Rd * kpc, Re * kpc, fb, 1.0 / 1.815 * kpc)
    Sig = 0.5 * M / (math.pi * (Re * 1e3)**2)
    P("   %-12s  %.1e   %4.1f     %.1f   %6.3f        %.3f       %.3f       %6.0f"
      % (name, M, Re, fb, aN / A0, fDM(aN, a0_fest(z)), fDM(aN, a0_H(z)), Sig))
P("   f_DM hängt bei festem a₀ nur an a_N(R_e)/a₀, also an der Flächendichte; a₀ ∝ H hebt es bei gleicher Scheibe mit z.")
P()
P("   Raster (reine Scheibe, B/T = 0): f_DM(R_e) fest / a₀∝H bei z = 1 und z = 2")
P("   M_b \\ R_e     3 kpc              5 kpc              7 kpc              9 kpc")
for M in (3e10, 6e10, 1e11, 2e11):
    row = []
    for Re in (3, 5, 7, 9):
        aN = aN_scheibe(M * M_sun, Re / 1.678 * kpc, Re * kpc)
        row.append("%.2f | %.2f/%.2f" % (fDM(aN, a0_fest(1)), fDM(aN, a0_H(1)), fDM(aN, a0_H(2))))
    P("   %.0e   %s" % (M, "   ".join(row)))
P("   (je Feld: fest | a₀∝H bei z=1 / z=2)")
P()
P("   Welche Flächendichte gibt die gemessenen Mediane [S: RC100-Abstract im Suchtreffer] 0,38 (z≈1) und 0,27 (z≈2)?")
for fm, z in ((0.38, 1.0), (0.27, 2.0), (0.28, 2.0)):
    x_f = math.log(1 / fm)**2                     # a_N/a0 aus q = e^{−√x}
    x_h = x_f * E_hintergrund(z)                  # a0 ∝ H: a_N/A0 = x·E
    P("   f_DM = %.2f bei z = %.0f: a_N(R_e)/a₀ = %.3f (fest), bei a₀∝H verlangt a_N(R_e)/A0 = %.3f (%.2f-fach)"
      % (fm, z, x_f, x_h, E_hintergrund(z)))
P("   Steigung an der Stelle: d f_DM/d log a₀ = f·s·ln10/2 (s = √(a_N/a₀)); a₀∝H verschiebt log a₀ um log E(z):")
for fm, z in ((0.38, 1.0), (0.27, 2.0)):
    s = math.log(1 / fm)
    fH = math.exp(-s / math.sqrt(E_hintergrund(z)))
    P("   gleiche Scheibe, f_DM fest = %.2f → a₀∝H bei z = %.0f: %.2f" % (fm, z, fH))

P("   Empfindlichkeit: d f_DM/d log a_N = −f·s·ln10/2; ein Fehler von 0,1 / 0,2 dex in M_b verschiebt f_DM:")
for fm in (0.38, 0.27):
    s_ = math.log(1 / fm); k_ = fm * s_ * math.log(10) / 2
    P("   f = %.2f: %.3f je dex → ±%.3f / ±%.3f; Abstand fest ↔ a₀∝H bei gleicher Scheibe: %.2f"
      % (fm, k_, 0.1 * k_, 0.2 * k_, math.exp(-s_ / math.sqrt(E_hintergrund(1.0 if fm > 0.3 else 2.0))) - fm))
P("   Median aus N Galaxien mit Streuung 0,23 / 0,18 [S]: Fehler 0,23/√50 = %.3f, 0,18/√50 = %.3f" % (0.23 / 50**0.5, 0.18 / 50**0.5))

# ------------------------------------------------------------------ 2
P()
P("2  Tully-Fisher-Nullpunkt")
P("   tiefer Grenzfall v⁴ = G M_b a₀: bei festem v Δlog M_b = −log(a₀(z)/a₀(0))")
for z in (0.9, 1.5, 2.3):
    P("   z = %.1f: fest %+.4f dex, a₀∝H %+.3f dex" % (z, -math.log10(a0_frueh_relativ(z)), -math.log10(E_hintergrund(z))))
P("   am Ort v(2,2 R_d) derselben Scheibe (R_d fest), M_b bei gleichem v gesucht:")


def v_at(M, Rd, a0, fb=0.0):
    R = 2.2 * Rd
    aN = aN_scheibe(M, Rd, R, fb, 1.0 / 1.815 * kpc)
    return math.sqrt(float(nu(aN, a0)) * aN * R)


for name, z, M, Re, fb in gal:
    Rd = Re / 1.678 * kpc
    v0 = v_at(M * M_sun, Rd, a0_fest(z), fb)
    MH = brentq(lambda m: v_at(m * M_sun, Rd, a0_H(z), fb) - v0, 1e8, 1e13)
    vH = v_at(M * M_sun, Rd, a0_H(z), fb)
    # Abstand zur tiefen Relation v⁴ = GMa₀ (fest)
    off = math.log10(M * M_sun / (v0**4 / (G * a0_fest(z))))
    P("   %-12s v(2,2R_d) = %5.1f km/s (fest), %5.1f (a₀∝H);  Δlog M_b bei festem v = %+.3f dex;  log M_b − log(v⁴/Ga₀) fest = %+.3f"
      % (name, v0 / 1e3, vH / 1e3, math.log10(MH / M), off))
P("   Gemessen [S, Übler 2017 Suchtreffer]: bei festem v höhere M_b bei z≈2,3 als bei z≈0,9 (positiver Nullpunktgang);")
P("   a₀∝H verlangt das Gegenteil (niedrigere M_b bei höherem z). Zahlenwerte der Messung nicht geholt.")

# ------------------------------------------------------------------ 3
P()
P("3  Druckstütze (Burkert 2010 [G]): v_c² = v_rot² + 2σ₀² R/R_d, bei R_e: 3,36 σ₀²")
for name, z, M, Re, fb in gal:
    Rd = Re / 1.678
    aN = aN_scheibe(M * M_sun, Rd * kpc, Re * kpc, fb, 1.0 / 1.815 * kpc)
    vc2 = float(nu(aN, a0_fest(z))) * aN * Re * kpc
    vbar2 = aN * Re * kpc
    for sig in (30e3, 50e3):
        vrot2 = vc2 - 3.356 * sig**2
        f_ohne = 1 - vbar2 / vrot2
        P("   %-12s σ₀ = %2.0f km/s: v_c = %5.1f, v_rot = %5.1f km/s; f_DM korrigiert %.3f, ohne Korrektur %+.3f (Δ = %+.3f)"
          % (name, sig / 1e3, math.sqrt(vc2) / 1e3, math.sqrt(vrot2) / 1e3, 1 - vbar2 / vc2, f_ohne, f_ohne - (1 - vbar2 / vc2)))
P("   Die Korrektur gilt für die Bahn wie in Newton (Jeans, radiales Gleichgewicht); die Regel ändert a, nicht den Druckterm.")
P("   Vergleich also mit v_c (korrigiert), wie die Quellen es tun.")

# ------------------------------------------------------------------ 4
P()
P("4  Preisach: frisch kollabierte Galaxien (Hernquist M_b = 6·10¹⁰, a = 2,57 kpc; g/g_RAR − 1 aus gewicht_preisach.txt)")
Mb = 6e10 * M_sun
ab = 1.5 * kpc * 6**0.3
cases = [("z_c=3,4 → z=2", 10, 0.0193), ("z_c=3,4 → z=2", 20, 0.0093),
         ("z_c=1 → z=0,5", 10, -0.0173), ("z_c=1 → z=0,5", 20, -0.0402), ("z_c=1 → z=0,5", 50, -0.0613)]
P("   Fall             R[kpc]  δg/g     a_N/A0   χ      dlng/dlnM   Δlog v [dex]  Δlog M_b [dex]   N(3σ) bei 0,10 / 0,15 / 0,20 dex")
for name, rk, d in cases:
    R = rk * kpc
    aN = G * hernquist_masse(Mb, ab, R) / R**2
    ch = float(chi(aN))
    slope = 1 - ch / 2
    dlv = 0.5 * math.log10(1 + d)
    dlM = math.log10(1 + d) / slope
    Ns = ["%6.0f" % ((3 * s / abs(dlM))**2) for s in (0.10, 0.15, 0.20)]
    P("   %-15s  %3d    %+.4f  %6.3f   %.3f   %.3f       %+.4f       %+.4f          %s" % (name, rk, d, aN / A0, ch, slope, dlv, dlM, " / ".join(Ns)))
P("   Δlog M_b: Masse, die bei gleichem g am selben Ort die Abweichung vortäuscht (M_b bei festem R und g).")
P("   Systematik: ein gemeinsamer Fehler im Masse-zu-Licht (SPARC: Υ_d 0,4–0,6 verschiebt a₀ um 0,9–1,5·10⁻¹⁰) bleibt bei jedem N;")
d10 = math.log10(1.0193) / (1 - float(chi(G * hernquist_masse(Mb, ab, 10 * kpc) / (10 * kpc)**2)) / 2)
P("   für 3σ muss dieser gemeinsame Nullpunkt der Massen auf Δ/3 = %.4f dex stimmen (z=2, 10 kpc)." % (d10 / 3))
P("   in v gemessen (Δlog v = 0,0042 dex): N(3σ) bei 0,03 / 0,05 dex Fehler je Galaxie = %.0f / %.0f"
  % ((3 * 0.03 / 0.0042)**2, (3 * 0.05 / 0.0042)**2))
# Vergleich: a₀∝H bei z=2, gleiche Galaxie 10 kpc
aN10 = G * hernquist_masse(Mb, ab, 10 * kpc) / (10 * kpc)**2
r_H = float(nu(aN10, a0_H(2))) / float(nu(aN10, a0_fest(2)))
P("   Zum Maß: a₀∝H gäbe bei z = 2 an derselben Stelle g/g_RAR − 1 = %+.3f (Faktor %.0f über dem Preisach-Signal)"
  % (r_H - 1, (r_H - 1) / 0.0193))

# ------------------------------------------------------------------ 5
P()
P("5  Handprüfung")
Ez2 = math.sqrt(0.3134 * 27 + 9.14e-5 * 81 + 0.6865)
P("   E(2) = √(0,3134·27 + 9,14e-5·81 + 0,6865) = √(8,4618 + 0,0074 + 0,6865) = √9,1557 = %.4f; Kern %.4f" % (Ez2, E_hintergrund(2)))
P("   → Δlog M_b (tief, a₀∝H, z=2) = −log 3,0258 = −0,4809")
P("   N(3σ), z=2, 10 kpc, σ = 0,15: Δ = log(1,0193)/(1−χ/2); Kern-χ s. o.; (3·0,15/Δ)²")
Rk = 10 * kpc; mh = 6e10 * (10 / (10 + 2.5695))**2
aNh = 6.6743e-11 * mh * 1.98841e30 / (3.0857e20)**2
sh = math.sqrt(aNh / 1.1e-10); chh = sh / (math.exp(sh) - 1)
Dh = math.log10(1.0193) / (1 - chh / 2)
P("   von Hand: M(<10) = 6e10·(10/12,57)² = %.3e M☉; a_N = %.3e m/s²; s = %.4f; χ = %.4f; Δ = %.5f dex; N = %.0f"
  % (mh, aNh, sh, chh, Dh, (0.45 / Dh)**2))
P("   f_DM von Hand: a_N/a₀ = 1 → q = e^{−1} = 0,368; a_N/a₀ = 2 → e^{−1,414} = 0,243; Kern: %.3f, %.3f"
  % (fDM(A0, A0), fDM(2 * A0, A0)))

with open(os.path.splitext(os.path.abspath(__file__))[0] + ".txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
