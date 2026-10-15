# -*- coding: utf-8 -*-
"""
Welche Laufzeit gilt für den Übergang des inaktiven Fluchtmaßes ins Atmen? (Fortsetzung von de_rest_anstieg)
1 Kandidaten aus dem Modell (vor dem Ergebnis): K0 sofort (Ankunft beim Eintritt in I), K1 Erstdurchgang bis ϖ_R
  (= Rückkehr zum Ort über die Erneuerung u(0) = g_N + u(ϖ_R), Lévy τ_D), K2 hin und zurück 0 → ϖ_R → 0 (Lévy 4τ_D),
  K3 Erstdurchgang bis Λ (Lévy (Λ/ϖ_R)²τ_D), K4 Durchquerung mit Reflexion bei Λ zurück zum Ort (Lévy (2Λ/ϖ_R)²τ_D).
2 Selbstkonsistent (wie S2 in de_rest_anstieg) je Kandidat: Ω_X(z), w₀, w_a, H/H_Λ, D_M(1090), Q-Fenster, fσ₈ [Messwerte nur Test, G].
3 Verlust unterwegs: Γ = a_N/c im gebundenen Bereich (je Masse, an der Naht und bei R_ta/2), θ̃ = 0.
4 Beobachterfehler.  5 Von Hand.
Aus ftm_kern: G, c, M_sun, kpc, Gyr, A0, tau_D, schliesslaenge, r_stern, rueckkehr, OMEGA_*.
Aus de_rest_anstieg (Kopf bis Abschnitt 1): schritt, iteriere, gesamt, X_von, Kugel, Tafeln (darin de_fluchtmoden_mass,
de_fluchtmoden_verzoegert, de_fluchtmoden).
"""
import math, os, sys, io, contextlib, warnings
import numpy as np
from scipy.special import erfc
from scipy.integrate import cumulative_trapezoid
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
from ftm_kern import (G, c, M_sun, kpc, Gyr, A0, OMEGA_DE, tau_D, schliesslaenge, r_stern, rueckkehr, omega_vakuum)

_src = open(os.path.join(HIER, "de_rest_anstieg.py"), encoding="utf-8").read().split("\n")
_stop = next(i for i, l in enumerate(_src) if l.startswith("# ============================================================== 1 Selbstkonsistenz"))
NS = {"__file__": os.path.join(HIER, "de_rest_anstieg.py"), "__name__": "dra"}
with contextlib.redirect_stdout(io.StringIO()):
    exec("\n".join(_src[:_stop]), NS)
dfm = NS['dfm']; AA = NS['AA']; fdich = NS['fdich']; T_F = NS['T_F']; NT = NS['NT']; NK = NS['NK']
OM_S = omega_vakuum()

OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)
def kopf(t):
    P(); P("=" * 112); P(t); P("=" * 112)

TD = tau_D() / Gyr
LAM = schliesslaenge()
T0 = 13.8

# ============================================================== 1 Kandidaten
kopf("1  Laufzeit des Übergangs: Kandidaten aus dem Modell")
P("   Rückkehrgleichung (07): ∂_t u + (Γ + θ)u = D∂_ϖ²u auf der Halbgeraden ϖ ≥ 0, Ort = Rand ϖ = 0 (Raum), Erneuerung")
P("   u(0) = g_N + u(ϖ_R), Schluss bei Λ. Erstdurchgang über die Strecke x: Lévy √(x²/4πD)·T^{−3/2}e^{−x²/4DT}, Skala x²/D.")
P("   τ_D = 2Λ/c = %.2f Gyr, Λ = %.3f Gpc." % (TD, LAM / (1e3 * kpc * 1e3)))
P()
P("   Λ/ϖ_R je Masse (ϖ_R = r* = √(GM_b/a₀), M_b = M/(1 + Q)):")
for lm in (8, 10, 12, 14):
    Mb = 10**lm * M_sun / NK
    rs = r_stern(Mb)
    P("     M = 10^%d M☉: ϖ_R = %8.2f kpc, Λ/ϖ_R = %.2e, (Λ/ϖ_R)²τ_D = %.1e Gyr" % (lm, rs / kpc, LAM / rs, (LAM / rs)**2 * TD))
P()
kand = [("K0", "sofort beim Eintritt in I", 0.0),
        ("K1", "Erstdurchgang bis ϖ_R = Rückkehr zum Ort (Erneuerung)", 1.0),
        ("K2", "hin und zurück 0 → ϖ_R → 0", 4.0),
        ("K3", "Erstdurchgang bis Λ (10¹² M☉)", (LAM / r_stern(1e12 * M_sun / NK))**2),
        ("K4", "Durchquerung, Reflexion bei Λ, zurück zum Ort (10¹² M☉)", (2 * LAM / r_stern(1e12 * M_sun / NK))**2)]
P("   Kandidat  Skala/τ_D     Gipfel [Gyr]   angekommen nach 13,8 Gyr   K̂(0; Γ) für die Regel")
for k, nm, sk in kand:
    if sk == 0:
        P("   %-3s  %-58s  0           –              1                          1 (keine Regel, q = 1)" % (k, nm)); continue
    frac = erfc(math.sqrt(sk * TD / (4 * T0)))
    P("   %-3s  %-58s  %-10.3g  %-12.4g   %-12.4g               e^{−%s·s}" % (k, nm, sk, sk * TD / 6, frac,
      "%.0f" % math.sqrt(sk) if sk < 100 else "%.1e" % math.sqrt(sk)))
# Gegenprobe K2 = Faltung zweier K1 (numerisch, anderer Weg als die Formel)
dt = 0.01; tg = np.arange(dt, 60, dt)
f1 = fdich(tg)
conv = np.convolve(f1, f1)[:len(tg)] * dt
lev4 = np.sqrt(4 * TD / (4 * math.pi)) * tg**-1.5 * np.exp(-4 * TD / (4 * tg))
P("   Gegenprobe K2 = K1 * K1 (numerische Faltung der Talbot-Dichte) gegen Lévy 4τ_D bei T = 20, 40, 59 Gyr: %s"
  % ", ".join("%.3e/%.3e" % (float(np.interp(T, tg + dt, conv)), float(np.interp(T, tg, lev4))) for T in (20, 40, 59)))
P()
P("   Zuordnung aus dem Modell:")
P("   – Die L_i sind Schleifen im Raum (01: die Rückführung erscheint als Dehnung der Raumschleifen). Der Raum ist der Rand")
P("     ϖ = 0 (07: Ort als unterste Ebene, G(0) = 1). Solange das Maß bei ϖ > 0 liegt, ist es in I, nicht in G: K0 und K3")
P("     (Ankunft beim Eintritt in I bzw. bei Λ) läsen eine Lage in I als Ankunft im Ganzen – Beobachterfehler.")
P("   – Ein Rücklauf ist ein Durchgang über die Tiefe ϖ_R, danach setzt die Erneuerung u(0) = g_N + u(ϖ_R) wieder am Ort ein;")
P("     'Rückkehr zum Ort' und 'Erstdurchgang bis ϖ_R' sind dasselbe Ereignis. Nur dieser Kern gibt K̂(0; Γ) = e^{−s} = q (07);")
P("     zwei Wege (K2) gäben e^{−2s}, also die Regel mit a₀/4. Das Maß geht in die Rückkehr über (05), also in diesen Kern: K1.")
P("   – Jeder stetige Weg von 0 nach Λ kreuzt ϖ_R vorher; an ϖ_R erneuert er sich am Ort. Λ setzt nicht die Laufzeit, sondern")
P("     die Kac-Zeit ⟨T⟩ = 2Λ/c = τ_D (über ϖ_R² = D·⟨T⟩) und macht den Durchgang vollständig (01: schließt bei Λ;")
P("     07: Halbgerade bei Λ geschlossen). Reflektierte Wege (K4) wiegen erst ab T ~ (Λ/ϖ_R)²τ_D (de_kanal_atmen: f/Lévy = 1).")
P("   – Die Kac-Rückkehr einer Irrfahrt mit Schritt h zum Randpunkt hat viele kurze Rückläufe; die Schrittlänge fällt in der")
P("     Regel heraus (07), gezählt wird der Durchgang bis zur Erneuerungstiefe: wieder K1.")
P("   ⇒ Erzwungen: K1, Lévy mit Skala τ_D (die bisher verwendete Dichte). Die Laufzeit ist kein freier Hebel.")

# ============================================================== 2 Selbstkonsistent je Kandidat
kopf("2  Selbstkonsistent (Bildung, Umkehr und Übergang im eigenen Atmen, wie S2) je Kandidat [Messwerte nur Test, G]")
uebergang_alt = NS['uebergang']
fdich_alt = NS['fdich']

def mache_uebergang(skala):
    def ueb(TA, X):
        TG = np.linspace(TA[0], TA[-1], NT); dtg = TG[1] - TG[0]
        sp = NS['CubicSpline'](TA, X)
        tq = TG / T_F
        S = sp(np.clip(tq, TA[0], None), 1) * (tq >= TA[0]) / T_F
        t_puls = TA[0] * T_F
        if skala == 0:
            R = S.copy()
            Xarr = cumulative_trapezoid(R, TG, initial=0.0) + X[0] * (TG >= t_puls)
            return TG, R, Xarr
        fk = lambda T: fdich_alt(np.asarray(T) / skala) / skala
        Fi = fk(np.arange(NT) * dtg)
        R = np.array([np.trapezoid(Fi[n::-1] * S[:n + 1], dx=dtg) for n in range(NT)])
        R += X[0] * fk(TG - t_puls) * (TG > t_puls)
        return TG, R, cumulative_trapezoid(R, TG, initial=0.0)
    return ueb

ERG = {}; RUN = {}
for k, skala in (("K1", 1.0), ("K2", 4.0), ("K0", 0.0)):
    NS['uebergang'] = mache_uebergang(skala)
    with contextlib.redirect_stdout(io.StringIO()):
        d1, v1, _ = NS['iteriere']("Λ")
        d2, v2, dm_ = NS['iteriere']("eigen", d1)
    RUN[k] = (d2, len(v2) - 1, dm_)
    ERG[k] = NS['gesamt'](d2)
    P("   %s: %d Runden (eigen), letzte Änderung max|ΔE/E| = %.1e" % (k, len(v2) - 1, dm_))
NS['uebergang'] = uebergang_alt

def zeile(k, e):
    return ("   %-8s %6.3f  %6.3f  %6.3f  %6.3f  %6.3f   %6.3f  %6.3f  %6.3f    %6.3f  %6.3f %6.3f %6.3f %6.3f   %6.3f     %5.2f   %6.3f  %6.3f %6.3f %6.3f  %6.3f"
            % (k, e['OX0'], OM_S + e['OX0'], e['OD'][0.5], e['OD'][1], e['OD'][2], e['w'], e['w0'], e['wa'], e['hz'][0],
               e['hz'][0.5], e['hz'][1], e['hz'][2], e['hz'][3], e['DM'], e['Q'], e['fs'][0], e['fs'][0.3], e['fs'][0.6],
               e['fs'][1.0], e['Xan0'] / e['X0']))
P()
P("   Kand Ω_X(0) Ω_DE(0)  z=0,5   z=1     z=2      w(0)    w₀      w_a     H/H_Λ z=0    0,5    1      2      3     D_M(1090)  Q oben  fσ₈/Λ z=0   0,3    0,6    1,0   angek./X")
for k in ("K1", "K2", "K0"):
    P(zeile(k, ERG[k]))
P("   K3, K4: bis heute nichts angekommen (erfc(√(10¹²)) = 0): Ω_DE = Ω_s = %.3f, w = −1, H/H_Λ(0) = %.3f." % (OM_S, math.sqrt(1 - 0.687 + OM_S)))
P("   Gemessen [G]: Ω_DE(0) = 0,687; Standardmodell Ω_DE-Anteil 0,394 / 0,215 / 0,068 bei z = 0,5 / 1 / 2.")
P("   K1 = S2 aus de_rest_anstieg (Ω_X(0) 0,328, z = 1: 0,131) – Gegenprobe der Neuimplementierung.")

# ============================================================== 3 Verlust unterwegs
kopf("3  Verlust unterwegs: Γ = a_N/c im gebundenen Bereich (θ̃ = 0)")
P("   Mit Γ: Dichte f(T)·e^{−ΓT}, gesamt ∫ = K̂(0; Γ) = e^{−s}, s² = Γτ_D = a_N/a₀; der Rest 1 − e^{−s} käme nie an.")
P("   Γ je Masse und Epoche der Bindung aus a_N = GM_b/r² an der Naht r = R_ta bzw. bei r = R_ta/2 (Baryonen, Newton).")
for a in (1 / 4, 1 / 2, 1.0):
    ss = []
    for lm in (8, 10, 12, 14):
        Mb = 10**lm * M_sun / NK; R = dfm.R_ta(Mb, a)
        ss.append((math.sqrt(G * Mb / R**2 / A0), math.sqrt(4 * G * Mb / R**2 / A0)))
    P("   a = %.2f: s (Naht / R_ta/2) für M = 10^8, 10^10, 10^12, 10^14: %s" % (a, ",  ".join("%.3f/%.3f" % x for x in ss)))

XB_STORE = {}
def X_von_bins(Dfun, dfak, rfak):
    X = []; XB = []
    for kk, a in enumerate(AA):
        Fg = erfc(dfm.delta_ta(a) * dfak[kk] / (math.sqrt(2) * dfm.SIG * Dfun(a)))
        Lv = np.interp(dfm.LM, dfm.LMg, NS['TB'][a] + NS['TS'][a] * math.log(rfak[kk]))
        dF = -np.gradient(Fg, dfm.LM * NS['LN10'])
        integ = dF * Lv
        msk = dfm.LM >= 6.0
        X.append(float(np.trapezoid(integ[msk], dfm.LM[msk] * NS['LN10'])))
        row = []
        for lo in np.arange(6.0, 16.5, 0.5):
            m = (dfm.LM >= lo - 1e-9) & (dfm.LM <= lo + 0.5 + 1e-9)
            row.append(float(np.trapezoid(integ[m], dfm.LM[m] * NS['LN10'])))
        XB.append(row)
    XB_STORE['XB'] = np.array(XB); XB_STORE['rfak'] = np.array(rfak)
    return np.array(X)

LMB = np.arange(6.0, 16.5, 0.5) + 0.25
def mache_uebergang_G(faktor_r):
    def ueb(TA, X):
        XB = XB_STORE['XB']; rf = XB_STORE['rfak']
        TG = np.linspace(TA[0], TA[-1], NT); dtg = TG[1] - TG[0]
        tq = TG / T_F
        a_b = np.exp(np.interp(np.clip(tq, TA[0], TA[-1]), TA, np.log(AA)))       # Epoche der Bindung
        lag = TG[:, None] - TG[None, :]
        Fm = np.where(lag > 0, fdich_alt(np.clip(lag, 0, None)), 0.0)
        R = np.zeros(NT); verl = 0.0; ges = 0.0
        for j, lm in enumerate(LMB):
            Mb = 10**lm * M_sun / NK
            sp = NS['CubicSpline'](TA, XB[:, j])
            S = sp(np.clip(tq, TA[0], None), 1) * (tq >= TA[0]) / T_F
            Rt = np.array([dfm.R_ta(Mb, a) for a in AA]) * rf / faktor_r
            aN = G * Mb / np.interp(np.log(a_b), np.log(AA), Rt)**2
            Gam = aN / c * Gyr                                                    # 1/Gyr, je Quellzeit
            Kmat = Fm * np.exp(-Gam[None, :] * np.clip(lag, 0, None))
            R += (Kmat * S[None, :]).sum(axis=1) * dtg
            t_puls = TA[0] * T_F; G0 = Gam[0]
            R += XB[0, j] * fdich_alt(TG - t_puls) * np.exp(-G0 * (TG - t_puls)) * (TG > t_puls)
            s_ = np.sqrt(Gam * TD)
            ges += np.trapezoid(S, TG) + XB[0, j]; verl += np.trapezoid(S * np.exp(-s_), TG) + XB[0, j] * math.exp(-s_[0])
        XB_STORE['anteil'] = verl / ges
        return TG, R, cumulative_trapezoid(R, TG, initial=0.0)
    return ueb

NS['X_von'] = X_von_bins
for k, fr in (("Γ-Naht", 1.0), ("Γ-R_ta/2", 2.0)):
    NS['uebergang'] = mache_uebergang_G(fr)
    with contextlib.redirect_stdout(io.StringIO()):
        d1, v1, _ = NS['iteriere']("Λ")
        d2, v2, dm_ = NS['iteriere']("eigen", d1)
    RUN[k] = (d2, len(v2) - 1, dm_)
    ERG[k] = NS['gesamt'](d2)
    P("   %s: %d Runden, max|ΔE/E| = %.1e; auf Dauer ankommender Anteil ∫S·e^{−s}/∫S = %.3f"
      % (k, len(v2) - 1, dm_, XB_STORE['anteil']))
NS['uebergang'] = uebergang_alt
P()
P("   Kand     Ω_X(0) Ω_DE(0)  z=0,5   z=1     z=2      w(0)    w₀      w_a     H/H_Λ z=0    0,5    1      2      3     D_M(1090)  Q oben  fσ₈/Λ z=0   0,3    0,6    1,0   angek./X")
for k in ("K1", "Γ-Naht", "Γ-R_ta/2"):
    P(zeile(k, ERG[k]))
P("   Zuordnung: 𝒫 erhält das Maß (05, doppelt stochastisch); Γ ist ein Taktverlust der Rückkehr am Ort (Regel: q < 1 heißt,")
P("   ein Anteil kehrt nicht zurück, nicht, dass Maß verschwindet). Ein Γ-Verlust des Fluchtmaßes nähme Maß aus dem geschlossenen")
P("   Ganzen – das verbietet die Erhaltung; der Maß-Übergang läuft ohne Γ (auf Dauer vollständig, 05). Γ dämpft nur die Regel.")
P("   θ: im gebundenen Bereich θ̃ = 0 (07, kleinste Rate der Geschichte); θ = 3H gäbe e^{−√(3H₀τ_D)} = 0,015 (de_fluchtmoden_verzoegert).")

# ============================================================== 4 Beobachterfehler
kopf("4  Beobachterfehler")
P("   Laufzeit: in I (Rückkehrzeit T über die Tiefe ϖ_R, keine Raumstrecke, keine Bahnzeit); τ_D = 2Λ/c fest (Λ Länge in I).")
P("   Ankunft: in G, als Dehnung der Raumschleifen L_i (H_eff = H_b + ⅓·dX_an/dt), erst am Rand ϖ = 0 (Raum).")
P("   Eine Lage in I (Eintritt, Tiefe Λ) als Ankunft in G zu lesen (K0, K3) wäre die Verwechslung der Blickrichtung.")
P("   Γ ist lokal (L, Gefälle am Ort) und dämpft die Rückkehr der Regel; das Maß ist G und erhalten.")
P("   Zustandsabhängig: X, X_an, Ω_X(z), t₀; fest: τ_D, Λ, a₀, Ω_s.")

# ============================================================== 5 Von Hand
kopf("5  Von Hand")
P("   Angekommen nach 13,8 Gyr (Puls bei t = 0): K1 erfc(√(86,36/55,2)) = erfc(1,2508) = %.4f; K2 erfc(2,5016) = %.2e."
  % (erfc(math.sqrt(TD / (4 * T0))), erfc(math.sqrt(4 * TD / (4 * T0)))))
Mb = 1e12 * M_sun / NK; R = dfm.R_ta(Mb, 1.0)
aN = G * Mb / R**2
P("   10¹² M☉ heute: M_b = %.3e kg, R_ta = %.0f kpc, a_N = GM_b/R_ta² = %.3e m/s² ⇒ s = √(a_N/a₀) = %.4f, e^{−s} = %.4f"
  % (Mb, R / kpc, aN, math.sqrt(aN / A0), math.exp(-math.sqrt(aN / A0))))
P("   Kern: Laplace mit Γ: rueckkehr(0, Γ) = %.4f (Kernfunktion, gleicher Wert)." % float(rueckkehr(0.0, aN / c)))
P("   Ω_X(0) K1 = E(1)² − E_b(1)² = %.4f² − 0,7467² = %.4f (Code %.4f)."
  % (ERG['K1']['E'](1.0), ERG['K1']['E'](1.0)**2 - 0.7467**2, ERG['K1']['OX0']))

with open(os.path.join(HIER, "de_laufzeit.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
