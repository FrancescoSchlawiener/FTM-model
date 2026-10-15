# -*- coding: utf-8 -*-
"""
Wogegen verliert ein Bereich in der Rückkehr?  (Rechnung, noch nicht im Modell)
Rückkehr auf der bei Λ geschlossenen Halbgeraden (innere Richtung I, Λ fest), Raumlänge des Bereichs L ∝ a (mitatmend).
Verlust = Anteil der Rückkehr, dessen Ziel im Raum sich während der Rückkehrzeit T gedehnt hat:
  1 − w = E[ D(T) ],  D(T) = exp(−∫_t0^{t0+T} θ_ref dt)
  (a) absolut gegen das feste Λ:   θ_ref = θ̃       (Hintergrund: D = (a(t0)/a(t0+T))³)
  (b) relativ zum Ganzen:          θ_ref = θ̃ − 3H  (Hintergrund: D = 1)
1 Formel und Probe (konstante Rate: e^{−√x})
2 Hintergrund z = 0 … 1090: (a) quasistatisch (Modellformel), (a) mit wirklichem a(t) längs T, (b)
3 Lineare Störungen (θ_lok = 3H(1 − ε), ε = fδ/3): (a) gegen (b)
4 Das Ganze: Kac-Mittel (ungedämpft, I) gegen gedämpften Anteil (Laplace, Rate aus L) auf derselben Dynamik
5 Ergebnis, 6 Handprüfungen
Nur Kern-Formeln (ftm_kern).
"""
import math, os, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftm_kern as K
from ftm_kern import c, Gyr, A0, Q_KOSM, OMEGA_M, OMEGA_R, E_hintergrund, a0_frueh_relativ, gewicht_w

OUT = []
def P(s=""):
    OUT.append(s); print(s)
def hdr(s):
    P(); P("=" * 110); P(s); P("=" * 110)

H0 = K.H0_si()
tauD0 = K.tau_D()
Lam0 = K.schliesslaenge()

def H(a):
    return H0 * E_hintergrund(1 / a - 1)

# Zeit t(a) = ∫ da/(aH); Tabelle in ln a, auch Zukunft (a bis 1e3)
lna = np.linspace(math.log(1e-6), math.log(1e3), 400001)
aa = np.exp(lna)
Ha = H0 * np.sqrt(OMEGA_M * aa**-3 + OMEGA_R * aa**-4 + (1 - OMEGA_M - OMEGA_R))
integrand = 1 / Ha
t_tab = np.concatenate([[aa[0] ** 2 / (2 * H0 * math.sqrt(OMEGA_R))], np.zeros(len(aa) - 1)])
t_tab[1:] = t_tab[0] + np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(lna))

def t_of_a(a):
    return np.interp(math.log(a), lna, t_tab)

def lna_of_t(t):
    return np.interp(t, t_tab, lna)

def f_rueck(T, tau):
    """Dichte erster Rückkehrzeiten aus der Tiefe ϖ (ϖ² = D·τ): √(τ/(4πT³))·e^{−τ/(4T)}."""
    return math.sqrt(tau / (4 * math.pi * T**3)) * math.exp(-tau / (4 * T))

def eins_minus_w_zeitabh(z, tau):
    """(a) mit wirklichem a(t): E[(a(t0)/a(t0+T))³] über die Rückkehrdichte (Integral in ln T)."""
    a0_ = 1 / (1 + z); t0 = t_of_a(a0_); l0 = math.log(a0_)
    def g(u):
        T = math.exp(u)
        l1 = lna_of_t(t0 + T)
        return f_rueck(T, tau) * T * math.exp(-3 * (l1 - l0))
    # log-Integral, Stützstellen dicht um den Sattel
    us = np.linspace(math.log(tau * 1e-6), math.log(tau * 1e4), 20001)
    vals = np.array([g(u) for u in us])
    return np.trapezoid(vals, us), us[np.argmax(vals)]

def eins_minus_w_konst(theta, tau):
    val, _ = quad(lambda u: f_rueck(math.exp(u), tau) * math.exp(u) * math.exp(-theta * math.exp(u)),
                  math.log(tau * 1e-8), math.log(tau * 1e4), limit=400, points=[math.log(tau / (2 * math.sqrt(max(theta * tau, 1e-9))))])
    return val

hdr("0  Grundgrößen")
P(f"   a₀ = {A0:.2e}, Λ = c²/2a₀ = {Lam0/K.Gpc:.4f} Gpc (I, fest seit z ≈ 10), τ_D = 2Λ/c = {tauD0/Gyr:.3f} Gyr, 3H₀τ_D = {3*H0*tauD0:.4f}")
P(f"   Q = {Q_KOSM:.4f}, ν_kosm = 1 + Q = {1+Q_KOSM:.4f}, Kac: ν_kosm = Λ/ΣΔs_b ⇒ ΣΔs_b = {Lam0/(1+Q_KOSM)/K.Gpc:.4f} Gpc (I, fest)")

hdr("1  Formel: Verlust = Dehnung des Ziels im Raum während der Rückkehrzeit in I")
P("   Die Rückkehr läuft in der inneren Richtung (Tiefe ϖ, ϖ² = Dτ_D, geschlossen bei Λ); ihr Ziel ist der Bereich im Raum (L ∝ a).")
P("   Kommt sie nach T an, hat sich das Ziel um exp(∫θ_ref dt) gedehnt; 1 − w = E[exp(−∫θ_ref dt)] ist der Anteil, der das Ziel trifft.")
P("   Konstante Rate θ: 1 − w = e^{−√(θτ_D)} (Laplace der Lévy-Dichte). Probe (Integral gegen Formel, τ = 1):")
for x in (0.01, 1.0, 17.86, 1e3):
    P(f"     x = {x:8.2f}:  Integral {eins_minus_w_konst(x, 1.0):.10e}   e^(−√x) {math.exp(-math.sqrt(x)):.10e}   Kern 1 − gewicht_w {1-gewicht_w(x*A0/c):.10e}")
P("   Abschluss bei Λ: 1 − w = cosh(k(Λ−ϖ))/cosh(kΛ); für ϖ ≪ Λ (jede Falte unter der Zelle) gleich e^{−√x} (Rechnung gewicht_rueckkehr).")

hdr("2  Hintergrund (mitbewegter Bereich, θ̃ = 3H): (a) absolut gegen Λ, (b) relativ zum Ganzen")
P("   (a-q) quasistatisch = Modellformel, θ̃ = 3H(z) konstant während T;  (a-t) wirkliches a(t) längs T (Zukunft aus dem Kern-Hintergrund);")
P("   τ_D mit a₀(z) (Kern a0_frueh_relativ) bzw. mit dem heutigen a₀ (Λ fest); (b) θ̃ − 3H = 0 ⇒ D = 1 ⇒ 1 − w = P(Rückkehr) = 1.")
P("   z        x = 3Hτ_D(z)    1 − w (a-q)          1 − w (a-t), τ(z)     Sattel T [Gyr]   1 − w (a-t), τ heute    w (b)")
res = {}
for z in (0, 0.5, 1, 3, 10, 100, 1090):
    a0z = A0 * a0_frueh_relativ(z)
    tz = c / a0z
    x = 3 * H(1 / (1 + z)) * tz
    wq = math.sqrt(x)
    wt, us = eins_minus_w_zeitabh(z, tz)
    wt0, _ = eins_minus_w_zeitabh(z, tauD0)
    res[z] = (x, wq, wt, wt0)
    P(f"   {z:<7}  {x:12.5e}   e^(−{wq:8.3f}) = 10^{-wq/math.log(10):8.2f}   {wt:10.4e}          {math.exp(us)/Gyr:8.3f}        {wt0:10.4e}            0")
P("   heute: w (a-q) = %.5f, w (a-t) = %.5f" % (1 - math.exp(-res[0][1]), 1 - res[0][2]))

# Hand-Sattel für z = 1090 (Materiezeit, a ∝ t^{2/3}): Exponent 2·ln(T/t0) + τ/(4T), Minimum bei T = τ/8
t1090 = t_of_a(1 / 1091)
P(f"   Hand z = 1090: t₀ = {t1090/Gyr*1e3:.3f} Mio. J.; mit a ∝ t^(2/3) ist D ≈ (t₀/T)²; Exponent 2 ln(T/t₀) + τ/(4T) minimal bei T = τ/8 = {tauD0/8/Gyr:.2f} Gyr")
T8 = tauD0 / 8
expo = 2 * math.log(T8 / t1090) + 2
P(f"     Wert 2·ln({T8/Gyr:.2f} Gyr/{t1090/Gyr*1e3:.3f} Mio. J.) + 2 = {expo:.2f} ⇒ 1 − w ~ e^(−{expo:.1f}) = {math.exp(-expo):.1e} (ohne Vorfaktor; Materiezeit bis T nur grob)")
P("   Bedingung aus der Hintergrundstrahlung: Der nicht mitatmende Anteil 1 − w wirkt mit ν(a_N) (+3–4 % Potential, C_ℓ 6–8 % bei 1 − w = 1),")
P(f"   linear in 1 − w: (a-t) gibt bei z = 1090 {res[1090][2]:.1e}·(6–8 %) = {res[1090][2]*0.06:.0e}–{res[1090][2]*0.08:.0e}: unsichtbar; (b) gibt 1 − w = 1: 6–8 %.")

hdr("3  Lineare Störungen, θ_lok = 3H(1 − ε), ε = f·δ/3 (Werte aus stoerungen_rekombination: θ_lok/3H)")
P("   (a): x = 3H(1 − ε)τ_D;  (b): x_rel = −3Hετ_D, nur für ε < 0 (Unterdichte) positiv, sonst w = 0")
P("   z      λ        ε = 1 − θ_lok/3H    1 − w (a) Über-/Unterdichte         w (b) Überdichte   w (b) Unterdichte")
for z, lam, rat in ((1090, "147 Mpc", 0.999847), (1090, "10 Mpc", 0.999184), (10, "147 Mpc", 0.985928), (0, "147 Mpc", 0.935892)):
    eps = 1 - rat
    a0z = A0 * a0_frueh_relativ(z); tz = c / a0z; x = 3 * H(1 / (1 + z)) * tz
    wo = math.sqrt(x * (1 - eps)); wu = math.sqrt(x * (1 + eps))
    wbu = -math.expm1(-math.sqrt(x * eps))
    P(f"   {z:<5}  {lam:8s} {eps:10.3e}          10^{-wo/math.log(10):8.2f} / 10^{-wu/math.log(10):8.2f}          0                  {wbu:.6f}")
P("   (b) macht das Gewicht vom Vorzeichen von δ abhängig: jede Überdichte (die Hälfte jeder Welle) wäre ab z = 1090 gebunden (w = 0),")
P("   die Regel ν(a_N) wirkte dort (a_N = 2,6–12·a₀(z)), die Mode wäre kein Staub; (a) behandelt beide Hälften gleich (w ≈ 1).")

hdr("4  Das Ganze: Kac-Mittel (I, ungedämpft) und gedämpfter Anteil (Rate aus L) auf derselben Rückkehr")
P("   Kette auf dem Kreis (M Zellen, Maß X = Λ), offener Bereich A (Maß ΣΔs_b = (1 − q)·Λ), Irrfahrt mit Verweilen, exakt (lineare Gleichungen).")
M = 200; nA = 32
A = np.zeros(M, bool); A[:nA] = True
Pm = np.zeros((M, M))
for i in range(M):
    Pm[i, i] = 0.5; Pm[i, (i + 1) % M] += 0.25; Pm[i, (i - 1) % M] += 0.25
def kac_und_laplace(theta_step):
    # h_i = E_i[e^{−θT_A}] für Start außerhalb von A (T_A = erster Eintritt), dann Rückkehr aus A: ein Schritt + h
    B = ~A; idx = np.where(B)[0]
    PB = Pm[np.ix_(idx, idx)]; PBA = Pm[np.ix_(idx, np.where(A)[0])].sum(1)
    d = math.exp(-theta_step)
    h = np.linalg.solve(np.eye(len(idx)) - d * PB, d * PBA)
    m = np.linalg.solve(np.eye(len(idx)) - PB, np.ones(len(idx)))   # mittlere Eintrittszeit
    hfull = np.ones(M); hfull[idx] = h
    mfull = np.zeros(M); mfull[idx] = m
    # Start gleichverteilt in A (stationär auf A)
    sA = np.where(A)[0]
    lap = np.mean([d * (Pm[i] @ hfull) for i in sA])
    mean = np.mean([1 + Pm[i] @ mfull for i in sA])
    return lap, mean
for th in (0.0, 0.001, 0.01, 0.1):
    lap, mean = kac_und_laplace(th)
    P(f"   Rate je Schritt θ = {th:5.3f}:  E[e^(−θT)] = {lap:.5f}   ⟨T⟩ = {mean:.4f}   (Kac μ(X)/μ(A) = {M/nA:.4f})")
P("   θ = 0 (Γ = 0): jede Rückkehr kommt an (q = 1), das Mittel ist das Kac-Maß ν = μ(X)/μ(A); die Dämpfung θ ändert ⟨T⟩ nicht,")
P("   nur den Anteil, der vor der Dehnung ankommt. q = 1 und ν_kosm = 1 + Q sind Aussagen über T in I (Maß, ungedämpft);")
P("   w ist eine Aussage über T gegen eine Rate im Raum. Beides gilt zugleich, für den Bereich wie für den Hintergrund des Ganzen.")
P(f"   Hintergrund des Ganzen, Lesart (a): Γ = 0 ⇒ q_Regel = 1;  θ̃ = 3H ⇒ w = {1-math.exp(-res[0][1]):.4f}: die Mode atmet mit (Staub),")
P("   und genau dieser Träger zählt im Ausdehnungsgesetz Ω_m = (1 + Q)·Ω_b. Lesart (b): w_Ganzes = 0 ⇒ Mode in Operatorform,")
P("   die ohne Gefälle nichts zählt ⇒ ν_kosm hätte im Hintergrund keinen Staubträger.")
P("   Bezugslänge: Λ und ΣΔs_b liegen in I und atmen nicht (seit z ≈ 10 fest), ν_kosm = Λ/ΣΔs_b ist zeitfest.")
P("   Atmete Λ mit (Λ ∝ a), wäre die Rate gegen Λ θ̃ − 3H (Volumen) — das ist Lesart (b); die Faltungsbilanz hält Λ fest ⇒ (a).")

hdr("5  Ergebnis")
P(f"   (a) absolut gegen das feste Λ: w_bg heute {1-math.exp(-res[0][1]):.4f} (quasistatisch) bzw. {1-res[0][2]:.4f} (wirkliches a(t)); z = 1090: 1 − w = e^(−738) bzw. {res[1090][2]:.1e};")
P("       Hintergrund und Störungen beider Vorzeichen atmen mit; ν_kosm = 1 + Q bleibt (Kac in I); Peaks unverändert.")
P("   (b) relativ zum Ganzen: w_bg = 0 in jeder Epoche, Überdichten w = 0 ab z = 1090 ⇒ Regel wirkt ⇒ Peaks 6–8 % daneben; trägt nicht.")
P("   Beobachterfehler: „Das Ganze verliert nichts“ gilt für Γ (Gefälle, Regel): q = 1, ν = Kac-Maß. Gegen Λ atmet auch der Hintergrund")
P("   des Ganzen (w = 0,985). Die Trennung ist nicht G gegen L, sondern Maß der Rückkehr (I, ungedämpft) gegen Dämpfung durch die Rate im Raum.")
P("   θ̃ wird im Raum gemessen (Divergenz der Strömung des Bereichs, L; im Hintergrund 3H für jeden Beobachter), es zählt als Rate")
P("   gegen die Rückkehrzeit in I (τ_D = 2Λ/c).")
P("   Setzung bleibt: ϖ² = Dτ_D (Tiefe = Kac-Diffusionstiefe); dass die Ziel-Dehnung (Volumen, 3H) der Verlust ist; Preisach θ̃ = min θ.")
P("   Die Modellformel (Rate zur Zeit t0 konstant) unterschätzt 1 − w, weil H während T fällt: heute %.4f statt %.4f." % (res[0][2], math.exp(-res[0][1])))

hdr("6  Unabhängige Prüfungen (von Hand)")
x0 = 3 * H0 * tauD0
P(f"   1) 3H₀τ_D: H₀ = 67,4/3,0857e19 = {67.4/3.0857e19:.4e} s⁻¹; τ_D = c/a₀ = 2,9979e8/1,1e-10 = {2.9979e8/1.1e-10:.4e} s; x = {3*67.4/3.0857e19*2.9979e8/1.1e-10:.3f}")
P(f"      √17,86 = 4,2261; e^(−4) = 0,018316, e^(−0,2261) = 0,79765 ⇒ 1 − w = 0,014610 (Rechnung {math.exp(-math.sqrt(x0)):.6f}) ⇒ w = 0,98539")
P(f"   2) z = 1090: x = {res[1090][0]:.4e}, √x = {math.sqrt(res[1090][0]):.2f} (Kontrolle gewicht_rueckkehr: 738,42)")
P("   3) Sattel z = 1090 (Materiezeit, oben): e^(−22…23) ≈ 1e-10; die Zahl der Tabelle muss in dieser Größenordnung liegen (Vorfaktor √(τ/T)… ~ 1).")
P(f"   4) Kette: Kac μ(X)/μ(A) = 200/32 = 6,25 (Tabelle θ = 0: {kac_und_laplace(0.0)[1]:.4f})")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "theta_bezug.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
