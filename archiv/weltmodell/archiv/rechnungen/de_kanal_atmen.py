# -*- coding: utf-8 -*-
"""
Kanal des inaktiven Fluchtmaßes: wirkt Maß, das über die Rückkehr ins Ganze kommt, (a) auf die Länge Λ der inneren
Richtung, (b) direkt auf die Schleifenlängen L_i (Atmen), oder (c) auf beides in festem Verhältnis (Schließung Λ ≙ ΔL)?
1 Minimalmodell: Maß je Baryon = Raum (L³) × innere Richtung (Λ) × Impulse; Rückkehr als Irrfahrt auf der bei Λ
  geschlossenen Halbgeraden (Transit), Kac, Λ aus der Faltung (Λ = ν_kosm·ΣΔs_b); Folgen je Kanal für Λ, a₀, Ω_s.
2 Fehlender Faktor 1,33 in X: 3D-Vollständigkeit der Tori (Monte Carlo, anderer Weg), Gas, Mode ξ (1 + Q), Einbettung
  unter 10⁶ M☉ (obere Schranke), Übergabe ab Umkehr gegen I_F.
3 Gesamtbild je Kanal: Ω_DE(z), w₀, w_a, H/H_Λ, D_M(1090), Q-Fenster, fσ₈ [Messwerte G].  4 Beobachterfehler.  5 Von Hand.
Aus ftm_kern: A0, c, Mpc, Gyr, OMEGA_*, omega_vakuum, schliesslaenge, zelle_aus_a0, q_fenster, d_LSS, E_hintergrund.
Aus de_fluchtmoden_mass / de_fluchtmoden_verzoegert (Hilfsfunktionen): mass_bereich, R_ta, Massenfunktion, lauf.
Lokal: Eigenreihe der Irrfahrt auf [0, ℓ], Monte Carlo der Newton-Tori, Kanalintegration (a), (c).
"""
import math, os, sys, io, contextlib, warnings
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp, quad
from scipy.special import erfc
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
from ftm_kern import (A0, c, G, M_sun, Mpc, Gyr, OMEGA_M, OMEGA_R, OMEGA_DE, Q_KOSM, H0_si, E_hintergrund,
                      omega_vakuum, schliesslaenge, zelle_aus_a0, q_fenster, d_LSS, KREISPAAR_GRENZE)

# ---- Hilfsfunktionen aus de_fluchtmoden_mass (Kopf bis vor die Tabellen; enthält den Kopf von de_fluchtmoden_verzoegert)
_src = open(os.path.join(HIER, "de_fluchtmoden_mass.py"), encoding="utf-8").read().split("\n")
_stop = next(i for i, l in enumerate(_src) if l.startswith("# ============================================================== 1 Herleitung"))
M_ = {"__file__": os.path.join(HIER, "de_fluchtmoden_mass.py"), "__name__": "dfmass"}
with contextlib.redirect_stdout(io.StringIO()):
    exec("\n".join(_src[:_stop]), M_)
V = M_['V']; dfm = M_['dfm']; AA = M_['AA']; AG = M_['AG']; i0 = M_['i0']; EBG = M_['EBG']; E_L = M_['E_L']
mass_bereich = M_['mass_bereich']; TG = V['TG']; H0 = V['H0']

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

NK = 1 + Q_KOSM
OM_S = omega_vakuum()
OM_X = OMEGA_DE - OM_S
LAM0 = schliesslaenge()                     # m
LZ = zelle_aus_a0(A0)                       # mitbewegte Zellenlänge, m
LN10 = math.log(10)

# ============================================================== 1 Minimalmodell
kopf("1  Kanal: Minimalmodell (Raumschleife L, innere Richtung Λ, gemeinsames Maß erhalten)")
P("   Maß je Baryon: μ₁ = L³·(2Λ)·(Impulsanteile). Jede Falte hat ihre eigene Halbgerade bis Λ (07), also trägt jedes Baryon")
P("   einen Faktor Λ wie einen Faktor V = L³. Erhaltung: δ ln μ₁ = 3 δ ln L + δ ln Λ = X (mitbewegte Impulse fest).")
P()
P("   (i) Rückkehr als Transit: Irrfahrt in ϖ auf [0, ℓ], Start in der Rückkehrtiefe ϖ_R = 1, Ort (Raum) bei ϖ = 0,")
P("       bei ℓ = Λ/ϖ_R geschlossen (reflektierend); D = 1 ⇒ τ_D = 1. Eigenreihe k_n = (n − ½)π/ℓ.")
def transit(ell, t, nmax=20000):
    k = (np.arange(1, nmax + 1) - 0.5) * math.pi / ell
    b = (2 / ell) * np.sin(k)
    f = np.array([np.sum(b * k * np.exp(-k**2 * x)) for x in np.atleast_1d(t)])          # Ankunft im Raum
    S = np.array([np.sum(b / k * np.exp(-k**2 * x)) for x in np.atleast_1d(t)])          # noch in I
    Tm = float(np.sum(b / k**3))
    return f, S, Tm
levy = lambda t: t**-1.5 * math.exp(-1 / (4 * t)) / math.sqrt(4 * math.pi)
P("       ℓ       ⟨T⟩/τ_D  (ϖ_R(2Λ−ϖ_R)/2D)   in I bei t = 1, 10, 100, 10⁴ τ_D       f(t)/Lévy bei t = 0,2  1  10")
for ell in (10.0, 100.0, 1000.0):
    ts = np.array([1.0, 10.0, 100.0, 1e4])
    f, S, Tm = transit(ell, ts, nmax=int(40 * ell) + 2000)
    fr = [transit(ell, np.array([x]), nmax=int(40 * ell) + 2000)[0][0] / levy(x) for x in (0.2, 1.0, 10.0)]
    P("   %8.0f   %10.1f   (%10.1f)          %s        %s"
      % (ell, Tm, (2 * ell - 1) / 2, "  ".join("%.4f" % s for s in S), "  ".join("%.4f" % r for r in fr)))
P("       Offen (ℓ → ∞): ⟨T⟩ = ∞, in I bei t: erf(1/√(4t)) (t = 10⁴: %.4f); geschlossen: alles kommt an, ⟨T⟩ endlich." % math.erf(1 / math.sqrt(4e4)))
P("       Heute ist Λ/ϖ_R ≈ 10⁶ (Milchstraße: ϖ_R = r* ≈ 8,6 kpc, Λ = 13,2 Gpc): f(T) = Lévy bis T ≈ 10¹¹ τ_D; die Schließung")
P("       ändert nicht die Laufzeit, sondern macht den Transit vollständig: Bei gleichmäßigem Zufluss r liegt in I nur")
P("       r·⟨T⟩ im Durchgang, fest; was zufließt, fließt in den Raum ab. Die Gleichung ist linear im Inhalt u auf festem [0, Λ]:")
P("       Der Inhalt verschiebt den Abschluss nicht (Kac: ⟨T⟩ = 2Λ/c hängt weder an Dynamik noch an Masse noch am Inhalt).")
P()
P("   (ii) Λ aus dem Modell: Λ = ν_kosm·ΣΔs_b (03), Δs = GM/3c² je Masse, unabhängig vom Bahnzustand; die gezählte Faltung")
P("       gebundener Bereiche hebt sich in der Naht auf (03). Bindung ändert ΣΔs_b nicht, ν_kosm = 1 + Q (mitatmend) auch nicht.")
P("       ⇒ δ ln Λ = 0 bei fester Faltung; aus der Erhaltung bleibt δ ln V = 3 δ ln L = X: Kanal (b).")
P()
P("   (iii) Kanäle gegeneinander (arrivierter Anteil X_an(t) aus der erzwungenen Kombination B, 10⁶, oben, θ̃ = 0, I_F):")
P("       (a) δ ln Λ = X_an ⇒ a₀ ∝ e^{−X_an}, Ω_s ∝ Λ⁻² ∝ e^{−2X_an}")
P("       (b) δ ln V = X_an ⇒ H_eff = H_b + ⅓·dX_an/dt, Λ fest")
P("       (c) Schließung Λ ≙ ΔL (03), δΛ = δL_i je Schleife: dΛ·(3/L + 1/Λ) = dX_an, δ ln a = dΛ/L; L = L_Zelle·a.")
P("           Anteil ins Atmen 3Λ/(3Λ + L) = %.3f heute (L_Zelle = %.1f Gpc, Λ = %.2f Gpc)."
  % (3 * LAM0 / (3 * LAM0 + LZ), LZ / (1e3 * Mpc), LAM0 / (1e3 * Mpc)))

# ---- Tabellen der Maße (Definition B) auf dem a-Gitter
TABM = {}
with contextlib.redirect_stdout(io.StringIO()):
    for a in AA:
        TABM[a] = [mass_bereich(10**l * M_sun / dfm.NK, dfm.R_ta(10**l * M_sun / dfm.NK, a)) for l in dfm.LMg]
FAK = {"k": 1.0}
def X_B(a, lm_min, art):
    Fg = erfc(dfm.delta_ta(a) / (math.sqrt(2) * dfm.SIG * dfm.Dg(a)))
    Lv = np.interp(dfm.LM, dfm.LMg, np.array([t["B"] for t in TABM[a]]))
    msk = dfm.LM >= lm_min
    dF = -np.gradient(Fg, dfm.LM * LN10)
    return FAK["k"] * float(np.trapezoid((dF * Lv)[msk], dfm.LM[msk] * LN10))
V['X_zaehlung'] = X_B

def kanal(d, art):
    """E(a) je Kanal aus dem arrivierten Maß X_an(t) = d['Xarr'], Rate d['R'] (1/Gyr)."""
    R = d['R']; Xa = d['Xarr']; Hg = H0 * Gyr
    Em2 = EBG**2 - OM_S                                      # Materie + Strahlung
    if art == "b":
        return EBG + R / (3 * Hg), np.ones_like(AG)
    if art == "a":
        lam = np.exp(Xa)
        return np.sqrt(Em2 + OM_S / lam**2), lam
    # (c): dΛ = dX/(3/L + 1/Λ), zusätzlich δ ln a = dΛ/L
    lam = np.ones_like(AG); Rb = np.zeros_like(AG)
    L0 = LZ / LAM0                                           # in Einheiten Λ₀
    for n in range(len(AG)):
        L = L0 * AG[n]
        g = 1 / (3 / L + 1 / lam[n])
        Rb[n] = 3 * R[n] * g / L                             # dln V/dt aus dem Atmen-Anteil
        if n + 1 < len(AG):
            lam[n + 1] = lam[n] + R[n] * g * (TG[n + 1] - TG[n])
    return np.sqrt(Em2 + OM_S / lam**2) + Rb / (3 * Hg), lam

def ausw(E, lam):
    OX = E**2 - (EBG**2 - OM_S)                              # gesamte Dunkle Energie (Vakuumteil + Übergang)
    a = AG; lna = np.log(a)
    wz = -1 - np.gradient(np.log(np.clip(OX, 1e-300, None)), lna) / 3
    m2 = (a >= 1 / 3) & (a <= 1.0)
    A_ = np.vstack([np.ones(m2.sum()), 1 - a[m2]]).T
    w0, wa = np.linalg.lstsq(A_, wz[m2], rcond=None)[0]
    at = lambda arr, z: float(np.interp(1 / (1 + z), a, arr))
    return dict(OD={z: at(OX, z) / at(E, z)**2 * 1.0 for z in (0, 0.5, 1, 2)}, OD0=float(OX[i0]),
                w=float(wz[i0]), w0=w0, wa=wa, hz={z: at(E, z) / E_L(1 / (1 + z)) for z in (0, 0.5, 1, 2, 3)},
                a0rel={z: 1 / at(lam, z) for z in (0, 1, 3)})

d6 = V['lauf'](6, "oben", "0", "I_F")
P()
P("   Kanal  Ω_DE(0) [0,687 G]  Ω_s(0)   a₀(z=1)/a₀(z=3)  a₀(0)/a₀(z=3)   H/H_Λ z=0    1      2     w(0)    w₀      w_a")
KAN = {}
for art in ("a", "b", "c"):
    E, lam = kanal(d6, art); e = ausw(E, lam); KAN[art] = (E, lam, e)
    P("   (%s)     %7.3f            %6.3f     %7.4f          %7.4f         %6.3f   %6.3f  %6.3f  %6.3f  %6.3f  %7.3f"
      % (art, e['OD0'], OM_S / float(np.interp(1.0, AG, lam))**2, e['a0rel'][1] / e['a0rel'][3], e['a0rel'][0] / e['a0rel'][3],
         e['hz'][0], e['hz'][1], e['hz'][2], e['w'], e['w0'], e['wa']))
P("   X_an(0) = %.3f (arriviert je Baryon), X(0) = %.3f (inaktiv).  Rückwirkung von Λ auf τ_D in (a), (c) nicht mitgeführt"
  % (float(d6['Xarr'][i0]), float(np.interp(1.0, AA, d6['X']))))
P("   (würde den Transit weiter verlangsamen).")
P("   (a) senkt Ω_DE unter den Vakuumteil (Vorzeichen falsch) und ließe a₀ seit z = 1 um %.0f %% fallen, auf Dauer (X → %.1f)"
  % (100 * (1 - KAN['a'][2]['a0rel'][0] / KAN['a'][2]['a0rel'][1]), float(np.interp(1.0, AA, d6['X']))))
P("       um den Faktor e^{−X} ≈ %.3f; das widerspricht Λ = ν_kosm·ΣΔs_b (a₀ seit z ≈ 10 bis auf 0,3 %% fest)." % math.exp(-float(np.interp(1.0, AA, d6['X']))))
P("   (c) teilt: %.0f %% ins Atmen, der Rest verlängert Λ und senkt Ω_s; netto kleiner als (b), a₀ wandert (%.1f %% seit z = 1)."
  % (100 * 3 * LAM0 / (3 * LAM0 + LZ), 100 * (1 - KAN['c'][2]['a0rel'][0] / KAN['c'][2]['a0rel'][1])))
P("       Λ ≙ ΔL (03) gilt für die Faltung (Rückstellung σ = ΣΔs − ΔL = 0), nicht für Maß der Bewegung: Bewegung faltet nicht.")
P("   ⇒ Kanal (b): Weg durch I (Laufzeit f(T), τ_D = 2Λ/c), Ankunft als Dehnung der Raumschleifen (G); Λ bleibt.")

# ============================================================== 2 Faktor
kopf("2  Fehlender Faktor in X (nötig 1,33: X(0) = 6,09 statt 4,56)")
# (i) 3D-Vollständigkeit: Monte Carlo der Newton-Tori innerhalb R (GM = R = 1) gegen Ellipsoid-Formel
rng = np.random.default_rng(7)
N = 4_000_000
r = rng.random(N)**(2 / 3)                                  # Dichte ∝ r²·v_max³ ∝ r^{1/2}
vmax = np.sqrt(2 / r)
u = rng.standard_normal((N, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
v = u * (vmax * rng.random(N)**(1 / 3))[:, None]
xr = np.zeros((N, 3)); xr[:, 0] = r
E = 0.5 * np.sum(v**2, 1) - 1 / r
Lv = np.cross(xr, v); L2 = np.sum(Lv**2, 1)
disc = 1 + 2 * E * L2
rapo = np.where(E < 0, (-1 - np.sqrt(np.clip(disc, 0, None))) / (2 * E), np.inf)
geb = rapo <= 1.0
Vs = (16 * math.pi**2 / 3) * 2**1.5 * (2 / 3)
mu_mc = Vs * geb.mean(); err = Vs * math.sqrt(geb.mean() * (1 - geb.mean()) / N)
mu_an = (16 * math.pi**2 / 3) * 2**1.5 * math.pi * (1.5 - math.sqrt(2))
vesc_R = np.sqrt(np.clip(2 * (1 / r - 1), 0, None))
mu_A = Vs * (np.linalg.norm(v, axis=1) <= vesc_R).mean()
P("   (i) Drehimpuls/Bahnebene: Monte Carlo im vollen 6D-Raum (r, v; alle Bahnebenen, L, L_z), gebunden ⇔ Apozentrum ≤ R:")
P("       μ_b(MC) = %.4f ± %.4f, Ellipsoid analytisch %.4f; μ_a(MC) = %.4f (π²·(16/3)·2^{3/2}·π/16 = %.4f); Verhältnis %.4f (16(3/2 − √2) = %.4f)"
  % (mu_mc, err, mu_an, mu_A, (16 * math.pi**2 / 3) * 2**1.5 * math.pi / 16, mu_mc / mu_A, 16 * (1.5 - math.sqrt(2))))
P("       Das Maß der Tori in 3D ist vollständig (Bahnebene und Drehimpuls stecken in d³v); kein Faktor.")
P("   (ii) Gas: Liouville zählt das Einteilchenmaß des Bereichs; ob ein Baryon auf Sternbahnen oder im Gas (Stöße, Druck)")
P("       läuft, ändert die Menge der Tori in der Naht nicht; F zählt alle Baryonen gebundener Bereiche schon mit. Kein Faktor;")
P("       Druck verhindert Bindung kleiner Bereiche, das senkt X eher (nur Richtung).")
X0 = float(np.interp(1.0, AA, d6['X']))
P("   (iii) Mode ξ: (1 + Q) auf X gäbe den Faktor %.2f (statt 1,33) und läse die Mode innen (w = 0) als bindenden Stoff; als" % NK)
P("       Bezug (D: ln B − 1,5 ln(1 + Q)) senkt sie X. Kein zwingender Anteil.")
# (iv) Einbettung unter 10⁶: obere Schranke (alle übrigen Baryonen heute in Bereichen der Masse M)
Fg0 = erfc(dfm.delta_ta(1.0) / (math.sqrt(2) * dfm.SIG * dfm.Dg(1.0)))
F6 = float(Fg0[0])
P("   (iv) Einbettung: F(>10⁶, z = 0) = %.3f; die übrigen %.3f der Baryonen liegen in keinem Bereich ≥ 10⁶ M☉." % (F6, 1 - F6))
for lm in (5.0, 4.0, 3.0):
    Mb = 10**lm * M_sun / NK
    lb = mass_bereich(Mb, dfm.R_ta(Mb, 1.0))["B"]
    P("       Lägen alle in Bereichen von 10^%.0f M☉ (ln B = %.3f): ΔX ≤ %.3f ⇒ X(0) ≤ %.3f, Faktor ≤ %.3f"
      % (lm, lb, (1 - F6) * lb, X0 + (1 - F6) * lb, (X0 + (1 - F6) * lb) / X0))
P("       Obere Schranke (voller Anteil, jüngste Epoche); die Bindung kleiner Bereiche ist früher und seltener: bleibt unter 1,33.")
# (v) Übergabe ab Umkehr gegen I_F
dU = V['lauf'](6, "oben", "0", "Umkehr"); EU, _ = kanal(dU, "b"); eU = ausw(EU, np.ones_like(AG))
FAK["k"] = 1.0
P("   (v) Übergabe ab Umkehr statt ab I_F: Ω_DE(0) = %.3f (I_F %.3f); X(0) gleich. I_F ist erzwungen (03), Umkehr nur Spannweite."
  % (eU['OD0'], KAN['b'][2]['OD0']))
P("   ⇒ Kein im Modell zwingender Anteil liefert 1,33; die Rechnung bleibt bei X(0) = %.2f." % X0)

# ============================================================== 3 Gesamtbild
kopf("3  Gesamtbild mit dem geklärten Kanal (b), Vergleich (a), (c) [Messwerte nur Test, G]")
def D_M(Efun, z=1090.0):
    zs = np.concatenate([np.linspace(0, 10, 4001), np.geomspace(10, z, 4001)[1:]])
    return np.trapezoid(1 / np.array([Efun(1 / (1 + x)) for x in zs]), zs)
dm_L = D_M(E_L)
def wachs(Efun):
    def f(x, y):
        a = math.exp(x); h = 1e-4
        dl = (math.log(Efun(a * math.exp(h))) - math.log(Efun(a * math.exp(-h)))) / (2 * h)
        return [y[1], -(2 + dl) * y[1] + 1.5 * OMEGA_M * a**-3 / Efun(a)**2 * y[0]]
    s = solve_ivp(f, (math.log(AG[0]), 0.0), [AG[0], AG[0]], dense_output=True, rtol=1e-8, atol=1e-12)
    return lambda a: s.sol(math.log(a))
gL = wachs(E_L)
for art in ("b", "a", "c"):
    E, lam, e = KAN[art]
    Ef = lambda a, E=E: float(np.interp(a, AG, E)) if a >= AG[0] else M_['EB'](a)
    dm = D_M(Ef)
    Lmin = KREISPAAR_GRENZE * d_LSS() * dm / dm_L
    gE = wachs(Ef)
    P("   (%s) Ω_DE(z) = %s; w(0) = %.3f, w₀ = %.3f, w_a = %.3f" % (art, ", ".join("z=%g: %.3f" % (z, e['OD'][z]) for z in (0, 0.5, 1, 2)), e['w'], e['w0'], e['wa']))
    P("       H/H_Λ: %s; D_M(1090)/D_M,Λ = %.3f; Q-Fenster oben %.2f (ΛCDM %.2f), unten 5"
      % (", ".join("z=%g: %.3f" % (z, e['hz'][z]) for z in (0, 0.5, 1, 2, 3)), dm / dm_L, q_fenster(L_min=Lmin)[1], q_fenster()[1]))
    P("       fσ₈/fσ₈,Λ (gleiche frühe Amplitude): " + ", ".join("z=%.1f: %.3f" % (z, gE(1 / (1 + z))[1] / gL(1 / (1 + z))[1]) for z in (0, 0.3, 0.6, 1.0, 1.5)))
P("   Ω_DE(0) gemessen [G] %.3f. (b) trägt 0,56 davon (82 %%); Form: w(0) ≈ −1, w_a ≪ 0 (Anstieg erst spät)." % OMEGA_DE)

# ============================================================== 4 Blick
kopf("4  Beobachterfehler")
P("   X: Maß des Zustandsraums (G des Bereichs, je Baryon über die Zelle G). Λ: Länge der inneren Richtung (I), gesetzt durch")
P("   die Faltung mit ihren Rückläufen, kein Behälter, dessen Größe am Inhalt hängt. Atmen: L_i(t) des Ganzen (G).")
P("   Der Weg läuft in I (Rückkehrzeit T, keine Raumstrecke), die Wirkung kommt in G an (Rate ⅓R, keine Dichte am Ort).")
P("   Kanal (a) läse Maß der Bewegung als Faltung (Stoff in der inneren Richtung): derselbe Fehler wie Regel als Stoff.")
P("   Zustandsabhängig: X, X_an, Ω_DE(z); fest: Λ, τ_D, a₀ (seit z ≈ 10).")

# ============================================================== 5 Von Hand
kopf("5  Von Hand")
P("   Anteil ins Atmen in (c): 3Λ/(3Λ + L) = 3·13,24/(3·13,24 + %.2f) = 39,72/%.2f = %.3f (Code %.3f)."
  % (LZ / (1e3 * Mpc), 39.72 + LZ / (1e3 * Mpc), 39.72 / (39.72 + LZ / (1e3 * Mpc)), 3 * LAM0 / (3 * LAM0 + LZ)))
P("   Transit ℓ = 100: ⟨T⟩ = ϖ_R(2ℓ − ϖ_R)/2D = (200 − 1)/2 = 99,5 τ_D (Eigenreihe %.1f)." % transit(100.0, np.array([1.0]), 6000)[2])
P("   Newton-Tori: μ_b/μ_a = 16(3/2 − √2) = 16·0,085786 = 1,3726 (Monte Carlo %.4f)." % (mu_mc / mu_A))
xa = float(d6['Xarr'][i0])
P("   Kanal (a) heute: Ω_s·e^{−2X_an} = %.3f·e^{−2·%.3f} = %.3f (Code %.3f)." % (OM_S, xa, OM_S * math.exp(-2 * xa), KAN['a'][2]['OD0']))

with open(os.path.join(HIER, "de_kanal_atmen.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
