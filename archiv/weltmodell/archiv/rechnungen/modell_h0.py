# -*- coding: utf-8 -*-
"""
H₀ aus dem eigenen Atmen des Modells und Prüfung am Winkel des Schallhorizonts θ* und an BAO.
Blick: fest sind die physikalischen Dichten (ω_b, ω_c = Q·ω_b, ω_r), a₀, Λ = c²/2a₀, der Vakuumteil π⁴/(15Λ²) (G, Zelle)
und das Fluchtmaß X_an(t) (selbstkonsistent, S2, K1). Anteile Ω = ω/h² hängen an H (Zustand): Folge, keine Eingabe.
1 H₀ aus der Flachheit, zwei Wege: (A) Iteration S2 (E_eff(1) in Einheiten des Bezugs 67,4); (B) absolut in km/s/Mpc:
  H₀ = 100·√(ω_m + ω_r + ω_s) + ⅓R(t₀), ω_s aus Λ_s = π⁴/(15Λ²) direkt, R aus der Faltung.  a₀ = 1,1 und 1,06 (SPARC).
2 θ* = r_s(z*)/D_M(z*) und BAO D_M/r_d, D_H/r_d [G] gegen Standardmodell (Kern-Werte).  3 Was sich ändern müsste.
4 Beobachterfehler.  5 Von Hand.
Aus ftm_kern: Konstanten, OMEGA_*, lambda_vakuum, schliesslaenge, omega_vakuum, E_hintergrund (Kern mit a₀ je Lauf).
Aus de_rest_anstieg (Kopf, Werkzeuge): schritt, iteriere, E_aus (Bildung im eigenen Atmen, S2).
Lokal: r_s-Integral (Schallgeschwindigkeit c/√(3(1+R_b))), z* nach Hu–Sugiyama, Abstände, Skalierung k·R, a₀ für θ*.
"""
import math, os, sys, io, re, types, builtins, contextlib, warnings
import numpy as np
from scipy.optimize import brentq
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

_open = builtins.open
@contextlib.contextmanager
def nur_lesen():
    """Die Hilfsmodule schreiben beim Laden ihre .txt; hier umgeleitet, damit keine Datei im Ordner geändert wird."""
    def op(f, mode="r", *a, **k):
        return _open(os.devnull, mode, *a, **k) if any(x in mode for x in "wax") else _open(f, mode, *a, **k)
    builtins.open = op
    try:
        yield
    finally:
        builtins.open = _open

KERN_SRC = _open(os.path.join(HIER, "ftm_kern.py"), encoding="utf-8").read()
DRA_SRC = _open(os.path.join(HIER, "de_rest_anstieg.py"), encoding="utf-8").read().split("\n")
DRA_STOP = next(i for i, l in enumerate(DRA_SRC) if l.startswith("# ============================================================== 1 Selbstkonsistenz"))

def lauf(a0):
    """Kern mit a₀, dann S2 (Bildung und Umkehr im eigenen Atmen) bis zur Selbstkonsistenz."""
    for m in ("ftm_kern", "de_fluchtmoden"):
        sys.modules.pop(m, None)
    K = types.ModuleType("ftm_kern"); K.__file__ = os.path.join(HIER, "ftm_kern.py")
    exec(compile(re.sub(r"(?m)^A0 = [0-9.e+-]+", "A0 = %r" % a0, KERN_SRC, count=1), K.__file__, "exec"), K.__dict__)
    sys.modules["ftm_kern"] = K
    ns = {"__file__": os.path.join(HIER, "de_rest_anstieg.py"), "__name__": "dra"}
    with nur_lesen(), contextlib.redirect_stdout(io.StringIO()):
        exec("\n".join(DRA_SRC[:DRA_STOP]), ns)
        d0 = ns['schritt'](ns['E_L'], "Λ")
        d1, v1, _ = ns['iteriere']("Λ", d0)
        d2, v2, dm = ns['iteriere']("eigen", d1)
    return dict(K=K, ns=ns, d=d2, E=ns['E_aus'](d2), runden=(len(v1) - 1, len(v2) - 1), dm=dm)

L11 = lauf(1.1e-10)
L106 = lauf(1.06e-10)
K = L11['K']
c, Mpc, Gyr = K.c, K.Mpc, K.Gyr
H_REF = K.H0_KMS                                    # Bezug der Einheiten (67,4), kein Ergebnis
h_ref = H_REF / 100
w_b, w_c = K.OMEGA_B_H2, K.OMEGA_C_H2
w_m = w_b + w_c
w_r = K.OMEGA_R * h_ref**2
w_g = 2.4728e-5 * (K.T_CMB / 2.7255)**4              # Photonen allein (Schallgeschwindigkeit)
KMS_GYR = (1 / Gyr) / (1e3 / Mpc)                    # 1/Gyr in km/s/Mpc
HH = 100e3 / Mpc                                     # 100 km/s/Mpc in 1/s

# ============================================================== 1 H₀
kopf("1  H₀ aus dem eigenen Atmen (Flachheit, kein freier Parameter)")
P("   Fest: ω_b = %.5f, ω_c = Q·ω_b = %.4f (Q = %.3f), ω_r = %.3e, a₀, Λ = c²/2a₀, Λ_s = π⁴/(15Λ²), X_an(t) aus S2."
  % (w_b, w_c, K.Q_KOSM, w_r))
P("   Die Rechnung läuft in absoluten Größen; 67,4 ist nur die Einheit von E (ρ_m, R_ta, σ, t in SI bzw. Gyr).")
ERG = {}
P()
P("   a₀        Runden  E_eff(1)  (A) H₀=67,4·E   ω_s        H_b [km/s/Mpc]  ⅓R(t₀)   (B) H₀     t₀[Gyr]   Ω_m     Ω_s     Ω_X     Σ")
for tag, L in (("1,10", L11), ("1,06", L106)):
    Kl = L['K']; d = L['d']; E1 = L['E'](1.0)
    H_A = H_REF * E1
    Lam = Kl.schliesslaenge(); lam_s = Kl.lambda_vakuum(Lam)             # 1/m²
    w_s = lam_s * c**2 / 3 / HH**2
    H_b = 100 * math.sqrt(w_m + w_r + w_s)
    t0 = d['t0']; R0 = float(np.interp(t0, d['TG'], d['R']))           # 1/Gyr
    H_B = H_b + R0 / 3 * KMS_GYR
    h = H_B / 100
    Om, Os = (w_m + w_r) / h**2, w_s / h**2
    OX = 1 - Om - Os
    ERG[tag] = dict(E1=E1, H_A=H_A, H_B=H_B, w_s=w_s, Om=Om, Os=Os, OX=OX, t0=t0, R0=R0, H_b=H_b, L=L)
    P("   %s·10⁻¹⁰  %d/%d    %.4f    %6.2f          %.4f     %6.2f          %5.2f    %6.2f    %6.2f   %.3f   %.3f   %.3f   %.3f"
      % (tag, *L['runden'], E1, H_A, w_s, H_b, R0 / 3 * KMS_GYR, H_B, t0, Om, Os, OX, Om + Os + OX))
P("   Kontrolle: de_rest_anstieg S2 (a₀ = 1,1) E_eff(1) = 0,9414; Ω_s bei 67,4: omega_vakuum() = %.4f, aus Λ_s: %.4f."
  % (K.omega_vakuum(), ERG["1,10"]['w_s'] / h_ref**2))
P("   Ω_X = 1 − Ω_m − Ω_s = (H₀² − H_b²)/H₀²: Flachheit schließt von selbst, weil H = H_b + ⅓R (Kanal b); kein freier Parameter.")
P("   (A) und (B) verschieden gerechnet (E-Gitter gegen absolute Summe aus Λ_s und R(t₀)); Abweichung %.2f bzw. %.2f km/s/Mpc."
  % (ERG["1,10"]['H_A'] - ERG["1,10"]['H_B'], ERG["1,06"]['H_A'] - ERG["1,06"]['H_B']))

# ============================================================== 2 θ*, BAO
kopf("2  Schallhorizont θ* und BAO [Messwerte G]")
def z_stern(wb, wm):
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763); g2 = 0.560 / (1 + 21.1 * wb**1.81)
    return 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)
def H_frueh(z):                                       # km/s/Mpc; Vakuumteil und Rückkehr vor z ~ 30 vernachlässigbar (gerechnet unten)
    return 100 * math.sqrt(w_m * (1 + z)**3 + w_r * (1 + z)**4)
def r_s(z1):
    lz = np.linspace(math.log(1 + z1), math.log(1 + 1e8), 200001)
    zz = np.exp(lz) - 1
    Rb = 3 * w_b / (4 * w_g) / (1 + zz)
    Hz = 100 * np.sqrt(w_m * (1 + zz)**3 + w_r * (1 + zz)**4)
    return np.trapezoid(c / 1e3 / np.sqrt(3 * (1 + Rb)) / Hz * (1 + zz), lz)   # Mpc
ZS = z_stern(w_b, w_m); ZD = 1059.9                 # z_d [G, Planck 2018]; r_s, r_d für beide Modelle gleich
RS, RD = r_s(ZS), r_s(ZD)
def D_M(Ef, z):                                      # Mpc, Ef(a) = H/H_REF
    zs = np.concatenate([np.linspace(0, min(z, 10), 8001), np.geomspace(10, z, 8001)[1:] if z > 10 else []])
    return c / 1e3 / H_REF * np.trapezoid(1 / np.array([Ef(1 / (1 + x)) for x in zs]), zs)
E_L = lambda a: K.E_hintergrund(1 / a - 1)
dmL = D_M(E_L, ZS)
P("   z* = %.1f (Hu–Sugiyama aus ω_b, ω_m), r_s(z*) = %.2f Mpc, r_d(z_d = %.1f) = %.2f Mpc – gleich im Modell und im" % (ZS, RS, ZD, RD))
Hz_ = H_REF * 1e3 / Mpc * E_L(1 / (1 + ZS))
P("   Standardmodell: vor der Rekombination c·H(z*) = %.1e m/s² ≫ a₀ (s = %.0f, q = e^{−s}), Θ = 2πa₀/(cH) = %.1e,"
  % (c * Hz_, math.sqrt(c * Hz_ / 1.1e-10), 2 * math.pi * 1.1e-10 / (c * Hz_)))
P("   Ω_s/Ω_m(z*) = %.0e, R = 0 vor z ≈ 30. Nur ω_b, ω_c, ω_r setzen r_s." % (K.omega_vakuum() / (K.OMEGA_M * (1 + ZS)**3)))
P()
P("   Fall                  H₀      D_M(z*) [Mpc]  100θ*    gegen 1,0411 [G]   D_M/D_M,Λ")
TH = {}
P("   Standardmodell (Kern) %5.2f   %9.1f      %.4f   %+.2f %%" % (H_REF, dmL, 100 * RS / dmL, 100 * (RS / dmL / 0.010411 - 1)))
for tag in ("1,10", "1,06"):
    e = ERG[tag]; dm = D_M(e['L']['E'], ZS); TH[tag] = dm
    P("   Modell a₀ = %s    %5.2f   %9.1f      %.4f   %+.2f %%            %.4f"
      % (tag, e['H_A'], dm, 100 * RS / dm, 100 * (RS / dm / 0.010411 - 1), dm / dmL))
P("   Planck-Fehler 100θ* ± 0,0003 [G] (0,03 %).")
P()
BAO = [(0.510, 13.588, 0.167, 21.863, 0.425), (0.706, 17.351, 0.177, 19.455, 0.330), (0.934, 21.576, 0.152, 17.641, 0.193),
       (1.321, 27.601, 0.318, 14.176, 0.221), (1.484, 30.512, 0.760, 12.817, 0.516), (2.330, 38.988, 0.531, 8.632, 0.101)]
P("   BAO (DESI DR2 [G]):  z      D_M/r_d: Mess   Λ(Kern)  Modell 1,10  1,06      D_H/r_d: Mess   Λ(Kern)  Modell 1,10  1,06")
chi = {"Λ": 0.0, "1,10": 0.0, "1,06": 0.0}
for z, dmm, sdm, dhm, sdh in BAO:
    row = {}
    for tag, Ef in (("Λ", E_L), ("1,10", ERG["1,10"]['L']['E']), ("1,06", ERG["1,06"]['L']['E'])):
        vm = D_M(Ef, z) / RD; vh = c / 1e3 / (H_REF * Ef(1 / (1 + z))) / RD
        row[tag] = (vm, vh); chi[tag] += ((vm - dmm) / sdm)**2 + ((vh - dhm) / sdh)**2
    P("                      %.3f   %6.2f±%.2f  %6.2f   %6.2f  %6.2f          %6.2f±%.2f  %6.2f   %6.2f  %6.2f"
      % (z, dmm, sdm, row["Λ"][0], row["1,10"][0], row["1,06"][0], dhm, sdh, row["Λ"][1], row["1,10"][1], row["1,06"][1]))
P("   χ² (12 Werte, ohne Korrelation): Standardmodell (Kern) %.1f, Modell 1,10: %.0f, 1,06: %.0f." % (chi["Λ"], chi["1,10"], chi["1,06"]))

# ============================================================== 3 Was sich ändern müsste
kopf("3  Was sich ändern müsste, damit θ* stimmt (nur benennen, keine Anpassung)")
L = ERG["1,10"]['L']; d = L['d']; ns = L['ns']
AGd, Ed = d['AG'], d['E']
EBf = lambda a: math.sqrt(K.OMEGA_M * a**-3 + K.OMEGA_R * a**-4 + K.omega_vakuum())
Rrel = np.array([Ed[i] - EBf(AGd[i]) for i in range(len(AGd))])          # ⅓R/H_REF
def E_k(k, os_fak=1.0):
    def f(a):
        eb = math.sqrt(K.OMEGA_M * a**-3 + K.OMEGA_R * a**-4 + os_fak * K.omega_vakuum())
        return eb + k * (float(np.interp(a, AGd, Rrel)) if a >= AGd[0] else 0.0)
    return f
ziel = RS / 0.010411
k_th = brentq(lambda k: D_M(E_k(k), ZS) - ziel, 0.5, 5.0, xtol=1e-4)
s_th = brentq(lambda s: D_M(E_k(1.0, s), ZS) - ziel, 0.5, 5.0, xtol=1e-4)
P("   Fehlend in D_M(z*): %.2f %% (Modell 1,10). Zwei Größen des Modells tragen die späte Rate:" % (100 * (TH["1,10"] / ziel - 1)))
P("   (i) Fluchtmaß: Rate ⅓R mal k = %.2f ⇒ H₀ = %.2f km/s/Mpc (X_an und R je Baryon größer, Form fest)." % (k_th, H_REF * E_k(k_th)(1.0)))
P("   (ii) Vakuumteil: Ω_s mal %.2f, also a₀ mal %.3f (a₀ = %.2f·10⁻¹⁰, Λ = c²/2a₀) ⇒ H₀ = %.2f (X fest gehalten)."
  % (s_th, math.sqrt(s_th), 1.1 * math.sqrt(s_th), H_REF * E_k(1.0, s_th)(1.0)))
P("   (iii) ω_c = Q·ω_b kleiner: senkt E früh und spät zugleich, ändert auch r_s – nicht gerechnet.")
for nm, Ef in (("(i)", E_k(k_th)), ("(ii)", E_k(1.0, s_th))):
    ch = 0.0
    for z, dmm, sdm, dhm, sdh in BAO:
        ch += ((D_M(Ef, z) / RD - dmm) / sdm)**2 + ((c / 1e3 / (H_REF * Ef(1 / (1 + z))) / RD - dhm) / sdh)**2
    P("   %s mit θ* getroffen: BAO χ² = %.1f (12 Werte); D_H/r_d bei z = 0,51: %.2f, z = 2,33: %.2f"
      % (nm, ch, c / 1e3 / (H_REF * Ef(1 / 1.51)) / RD, c / 1e3 / (H_REF * Ef(1 / 3.33)) / RD))
P("   θ* legt ∫dz/H fest; weil die Rückkehr-Rate spät anhebt, verlangt sie ein höheres H₀ als das Standardmodell.")

# ============================================================== 4 Beobachterfehler
kopf("4  Beobachterfehler")
P("   Fest (G, Zelle): ω_b, ω_c, ω_r, Λ_s = π⁴/(15Λ²), a₀, τ_D = 2Λ/c; X_an(t) je Baryon (Zustandsraum, I → G).")
P("   Zustandsabhängig: H(t) und alle Anteile Ω = ω/h²; die 0,687 ist ein Anteil beim H des Standardmodells.")
P("   Der Test mit θ* und BAO vergleicht Winkel und Verhältnisse (r_s/D_M, D/r_d): unabhängig von einem gesetzten H₀.")
P("   Die Lücke erscheint dort nicht als Anteil, sondern als zu langsame späte Rate ∫dz/H: sie ist kein Bezugsfehler.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand")
e = ERG["1,10"]
P("   H₀ = 67,4·0,9414 = 63,45 km/s/Mpc (Code (A) %.2f).  Ω_s = 0,2440·(67,4/63,45)² = 0,2440·1,1284 = 0,2753 (Code %.4f)."
  % (e['H_A'], e['Os']))
P("   Ω_m = 0,14197/0,6345² = 0,14197/0,40259 = 0,3526 (+ Strahlung ⇒ Code %.4f)." % e['Om'])
P("   θ*: 1,0411 ist gemessen; Modell/Standard = 1/%.4f ⇒ 100θ* ≈ 1,0411/%.4f = %.4f bei gleicher Eichung."
  % (TH["1,10"] / dmL, TH["1,10"] / dmL, 1.0411 / (TH["1,10"] / dmL)))

with _open(os.path.join(HIER, "modell_h0.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
