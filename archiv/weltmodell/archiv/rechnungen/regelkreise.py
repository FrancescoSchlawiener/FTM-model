"""Regelkreise je Ebene und durch die Epochen (Kybernetik, Systemtheorie, Thermodynamik).
Jede Ebene: Schleifenverstärkung q, Auslass 1 − q, Verstärkung ν = 1/(1 − q), Kac-Budget 1/ν + q = 1.
Alle Formeln und Konstanten aus ftm_kern. Ausgabe regelkreise.txt."""
import math, sys, io
import ftm_kern as K

out = io.StringIO()
def p(*a): print(*a, file=out)
def kopf(t): p("=" * 110); p(t); p("=" * 110)

kpc, Mpc = 1e3 * K.pc, 1e6 * K.pc
Ms = K.M_sun
MW = 6e10 * Ms

# ------------------------------------------------------------------ 1 Tafel je Ebene
kopf("1 Regelkreis je Ebene (L am Ort; Kosmos G/I). q = e^{−√(a_N/a₀)}, |q_Δs| = Δs/r, Takt = Umlauf 2πr/v bzw. τ_D")
def zeile(name, M, r, ein, aus, takt=None):
    aN = K.G * M / r**2
    q = float(K.q_(aN)); nu = float(K.nu(aN))
    v = math.sqrt(K.a_regel(aN) * r)
    T = takt if takt is not None else 2 * math.pi * r / v
    p(f"{name:<22} a_N/a₀={aN/K.A0:9.3e}  q={q:9.3e}  1−q={1-q:9.3e}  ν={nu:9.4g}  "
      f"|q_Δs|={K.ds(M)/r:8.2e}  Takt={T/K.yr:9.3e} a  1/ν+q={1/nu+q:.6f}")
    p(f"{'':<22} ein: {ein};  aus: {aus}")
    return q, nu

zeile("Teilchen (Proton, R_p)", K.m_p, K.R_p, "Faltung Δs (fest), Windung verlustfrei",
      "nichts im Betrag; Kreis offen (Newton), nur Phase läuft um", takt=K.hbar / (K.m_p * K.c**2))
zeile("Stern (Sonne, R☉)", Ms, K.R_sun, "Faltung aus Gas", "Strahlung (Faltung geht rotverschoben weg), Entropie")
zeile("Galaxie innen (5 kpc)", MW, 5 * kpc, "Gas, Faltung", "Sternlicht; Kreis fast offen")
rs = K.r_stern(MW)
zeile("Galaxie außen (r*)", MW, rs, "Rückläufe beginnen", "Auslass e^{−1} offen: Kreis gekoppelt")
RQ = K.r_Q(MW)
zeile("R_Q (q = q_kosm)", MW, RQ, "örtliche = kosmische Schleife", "Übergabe der Mode ξ vollständig")
Mges = (1 + K.Q_KOSM) * MW
Rta = K.umkehrradius_heute(Mges)
zeile("Umkehrradius R_ta", MW, Rta, "Kreis fast geschlossen (nahe Kippstelle q=1)", "Naht: jenseits w→1, Anschluss ans Atmen (1+wQ)")
Mh = 1.7e14 * Ms
rho_c = 3 * K.H0_si()**2 / (8 * math.pi * K.G)
R200 = (3 * (1 + K.Q_KOSM) * Mh / (4 * math.pi * 200 * rho_c))**(1 / 3)
zeile("Haufen (R₂₀₀)", Mh, R200, "Galaxien als Wirbel, Einfall", "Gas heiß (Tiefe als Temperatur), Übergabe nach wenigen Radiant")
p(f"   Milchstraße: r* = {rs/kpc:.2f} kpc, R_Q = {RQ/kpc:.1f} kpc, R_ta = {Rta/Mpc:.3f} Mpc; Haufen R₂₀₀ = {R200/Mpc:.2f} Mpc, M_b,Q = {K.masse_bQ()/Ms:.2e} M☉")
# Kosmos
qk, nuk = K.q_KOSM, 1 + K.Q_KOSM
Lam = K.schliesslaenge(); tD = K.tau_D()
Rh = K.d_LSS(1e8) / 2
rho_b = K.OMEGA_B * rho_c
Mb_h = rho_b * 4 / 3 * math.pi * Rh**3
p(f"{'Kosmos (Zelle)':<22} q={qk:.4f}  1−q={1-qk:.4f} (=Ω_b/Ω_m, offener Teil = Ort)  ν={nuk:.4f}  "
  f"ΣΔs/R(Horizont)={K.ds(Mb_h)/Rh:.3f}  Takt τ_D=2Λ/c={tD/K.Gyr:.1f} Gyr  1/ν+q={1/nuk+qk:.6f}")
p(f"{'':<22} ein: Faltung der Baryonen (Λ = ν·ΣΔs_b = {Lam/(1e9*K.pc):.1f} Gpc, innere Richtung);  aus: Rückstellung (Vakuumteil) + Atmen")
p("   Lesart: Verstärkung wächst von innen nach außen monoton bis zur Naht; dort springt das Gewicht w und der Kreis")
p("   öffnet in den kosmischen Kreis (fester Auslass 1 − q_kosm). Ein Kreis im Kreis mit Begrenzer an der Naht.")

# größte Verstärkung je Masse (am Umkehrradius)
kopf("1b Größte örtliche Verstärkung ν(R_ta) je Masse (Begrenzer Naht): tief ν ≈ 1/s, s ∝ M_b^{1/6}")
for Mb in [1e8, 1e10, 6e10, 1e12, 1e14]:
    R = K.umkehrradius_heute((1 + K.Q_KOSM) * Mb * Ms)
    aN = K.G * Mb * Ms / R**2
    p(f"   M_b = {Mb:8.1e} M☉:  R_ta = {R/Mpc:6.3f} Mpc  a_N/a₀ = {aN/K.A0:8.2e}  ν_max = {float(K.nu(aN)):7.1f}  q = {float(K.q_(aN)):.4f}")

# ------------------------------------------------------------------ 2 Epochen
kopf("2 Epochen: kosmischer Kreis (Form q_kosm fest), Schwelle a₀(z), Θ = 2πa₀/(cH), Hintergrundrest 1 − w (θ̃ = 3H)")
p(f"{'Epoche':<26}{'z':>10}{'a₀(z)/a₀':>10}{'τ_D·H':>12}{'Θ':>10}{'1−w_bg':>11}{'Ω_s':>11}   Kreise")
ep = [("Nukleosynthese", 4e8, "alle örtlichen offen; Schwelle 0 (Λ strahlungsgefaltet); kosmischer Kreis q fest, mitatmend"),
      ("Rekombination", 1090, "offen/mitatmend (w≈1): Mode wie Staub, Auslass → Atmen"),
      ("erste Umkehr", 20, "erste kleine Kreise schließen (10⁸ M☉), Rest mitatmend"),
      ("a₀ fest", 10, "Λ fest: Takt des Ganzen steht; örtliche Kreise koppeln an festen Sollwert"),
      ("Bildung/Übergabe", 2, "Galaxien schließen (w→0), Mode übergeben; Kehrseite → Atmen (verzögert τ_D)"),
      ("DE-Epoche", 0.5, "Haufen schließen; Vakuumteil trägt; Fluchtmaß kommt an"),
      ("heute", 0, "Θ≈1: Horizonttakt = Schwellentakt (Übereinstimmung von heute)")]
for n, z, t in ep:
    H = K.H0_si() * K.E_hintergrund(z)
    a0z = K.A0 * K.a0_frueh_relativ(z)
    Th = K.theta_(a0z, H)
    rest = 1 - K.gewicht_w(3 * H, a0z)
    Om_s = K.lambda_vakuum(K.schliesslaenge(a0z)) * K.c**2 / (3 * H**2)
    p(f"{n:<26}{z:>10.3g}{a0z/K.A0:>10.4f}{K.tau_D(a0z)*H:>12.4g}{Th:>10.3g}{rest:>11.3e}{Om_s:>11.3e}   {t}")
# Zukunft: fester DE-Deckel
p("   Zukunft (H → H_∞, fester Deckel k_DE = n·k_s): Θ_∞ = √(45/n)/π, 1 − w_∞ = e^{−√(6π/Θ_∞)}")
for nn in [1, 3]:
    Thi = math.sqrt(1 / (nn * K.K_S))
    p(f"     n = {nn}: Θ_∞ = {Thi:.4f},  1 − w_∞ = {math.exp(-math.sqrt(6*math.pi/Thi)):.4e}")

# ------------------------------------------------------------------ 3 Thermodynamik
kopf("3 Bilanz von Ordnung und Entropie")
H0 = K.H0_si()
SH = K.horizontentropie(H0)
Hs = K.H_vakuum(Lam)
SL = K.horizontentropie(Hs)
p(f"   Horizontentropie heute S_H = {SH:.3e} k; Deckel Vakuumteil allein S_Λ = {SL:.3e} k (S_H/S_Λ = {SH/SL:.3f})")
p(f"   Deckel mit k_DE = 3k_s: S_∞ = {SL/3:.3e} k (S_H/S_∞ = {SH/(SL/3):.3f}) – S_H läuft auf ihn zu")
p("   Bindung ordnet innen (inaktive Fluchtmoden, 𝒫 erhält das Maß) → Kehrseite muss aus dem gebundenen Bereich:")
p("     Rückstellung: hängt nur an der Faltung ΣΔs, nicht am Maß → trägt die Kehrseite nicht")
p("     Atmen (Fluchtmaß): nimmt das Maß auf, wächst mit der Bindung → nicht zeitfest")
p("   Zeitfest: Λ (seit z≈10), q_kosm (Form), der Deckel k_DE. Wachsend: Entropie, gebundenes Maß, Zahl der Umkehrungen, S_H → Deckel.")
r, S, _ = K.kreislauf_mindestrate(Lam)
p(f"   Ist S_H am Deckel, trägt nur wachsendes Λ die Entropie: Λ̇/Λ ≥ {r:.2e} 1/s = {r*K.yr:.1e} je Jahr")

# ------------------------------------------------------------------ 4 Dunkle Energie: welcher Auslass
kopf("4 Welcher Auslass kann der zeitfeste Rest sein?  nötig n = k_DE/k_s")
k_s = K.K_S * (2 * math.pi)**2 * 3 / (8 * math.pi)
p(f"   k_s = (π²/45)(2π)²·3/(8π) = {k_s:.4f} (π³/30 = {math.pi**3/30:.4f})")
omega_DE = 0.317                       # aus θ* [G] (Rechnung de_invariante)
H100 = 100e3 / (1e6 * K.pc)
def n_noetig(a0): return omega_DE * 3 * H100**2 / (8 * math.pi * K.G) * K.G * K.c**2 / a0**2 / k_s
def a0_fuer(n): return math.sqrt(omega_DE * 3 * H100**2 / (8 * math.pi) * K.c**2 / (n * k_s))
p(f"   nötig: n(1,06) = {n_noetig(1.06e-10):.3f}, n(1,10) = {n_noetig(1.10e-10):.3f}, SPARC 1,04–1,19: {n_noetig(1.19e-10):.2f}–{n_noetig(1.04e-10):.2f}")
p("   (a) Zerlegung derselben Rückkehr in Teilkreise: Σ1/T_i = 1/T fest, Vakuumteil ∝ Σ(1/T_i)² ≤ (Σ1/T_i)² ⇒ n ≤ 1")
for m in [2, 3]:
    p(f"       {m} gleiche Teilkreise: n = {m*(1/m)**2:.3f}")
Q = K.Q_KOSM
p(f"   (b) geschlossener Teil q des kosmischen Kreises mit eigener Rückstellung: ∝ (1/T_q)² → n = 1 + Q² = {1+Q**2:.1f};  linear 1 + Q = {1+Q:.2f}")
p("   (c) Kehrseite örtlicher Bindung: wächst mit der Bildung (nicht zeitfest, BAO-Form verlangt fest)")
p("   (d) Kehrseite des kosmischen Kreises: mitatmend (w ≈ 1) ⇒ keine inaktiven Fluchtmoden ⇒ 0")
p("   (e) Zahlen im Fenster ohne Herleitung (a₀, das sie verlangen):")
for nm, n in [("3 (Raumschleifen)", 3.0), ("π", math.pi), ("ν_kosm/2", nuk / 2), ("1 + ln ν_kosm", 1 + math.log(nuk)), ("√ν_kosm", math.sqrt(nuk))]:
    p(f"       {nm:<20} n = {n:.3f}  →  a₀ = {a0_fuer(n)*1e10:.4f}·10⁻¹⁰")
p("       Fenster n = 2,86–3,08 ist 7 % breit; jede einfache Zahl darin passt. Das Bild erzwingt keine.")

# ------------------------------------------------------------------ 5 Beobachterfehler
kopf("5 Beobachterfehler je Ebene (L am Ort, I innere Richtung, G Ganzes; Zustand)")
for t in ["Teilchen/Stern (L): q ≈ 0 heißt offener Kreis, nicht 'keine Rückkehr'; Strahlung ist Faltungsverlust, keine Wärme des Mediums.",
          "Galaxie (L): ν ist Zählung der Rückläufe am Ort, keine Masse; die gelesene Zusatzmasse ist Ausgabe der Rechnung.",
          "R_Q, R_ta (L/G): 1 + Q ist das Maß des Ganzen, keine Grenze am Ort; ν_max ist ein Zustand (w = 0 ab Umkehr), wächst mit R_ta(t).",
          "Haufen (L, Zustand): Übergabe hängt an der Geschichte (Preisach), nicht am heutigen Feld.",
          "Kosmos (I): Λ und τ_D liegen in der inneren Richtung, keine Raumlänge, keine Laufzeit im Raum; q_kosm ist Form (zeitfest).",
          "Kosmos (G, Epoche): Θ, τ_D·H, S_H/S_Λ sind Größen mit laufender Uhr, gelten je Beobachter und Epoche; zeitfest nur q, Λ, k_DE.",
          "Vakuumteil/Deckel (G): Rückstellung als Energie gelesen, kein Stoff am Ort; der Deckel ist ein Maß des Ganzen."]:
    p("   " + t)

# ------------------------------------------------------------------ 6 von Hand
kopf("6 Von Hand (anderer Weg)")
GM = K.GM_sun * 6e10
rs_h = math.sqrt(GM / 1.1e-10)
lq = math.log(5.364 / 6.364)
p(f"   GM = 1,32712e20·6e10 = {GM:.4e}; r* = √(GM/a₀) = {rs_h:.4e} m = {rs_h/kpc:.3f} kpc")
p(f"   ln(5,364/6,364) = {lq:.5f}; R_Q = r*/0,17094 = {rs_h/(-lq)/kpc:.2f} kpc (Code {RQ/kpc:.2f})")
p(f"   Kac: 1/6,364 + 0,8429 = {1/6.364 + 5.364/6.364:.6f}")
p(f"   Θ_∞(n=1) = √45/π = {math.sqrt(45)/math.pi:.4f} (Code {math.sqrt(1/K.K_S):.4f})")
p(f"   n = 3 ⇒ a₀ = 1,1·√(2,856/3) = {1.1*math.sqrt(2.856/3):.4f} (Code {a0_fuer(3)*1e10:.4f})")

txt = out.getvalue()
open(__file__.replace('.py', '.txt'), 'w').write(txt)
print(txt)
