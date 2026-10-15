# -*- coding: utf-8 -*-
"""
Triade auf allen Ebenen: Hintergrund – Wirbel – gefaltetes Zentrum (01), mit den Größen des Modells.
Nur ftm_kern und weltmodell/daten (SPARC). Annahmen mit [A].
 1 Tabelle je Ebene (Zentrum, Wirbel, Hintergrund, Übergänge, Zustand)
 2 Verhältnisse: Identitäten (auf allen Ebenen gleich), Skalengesetze (Form), Brüche
 3 Unvollständige Triaden (Teilchen unter M_c, Atom, kompakte Ebenen, Haufen, Kosmos)
 4 Vorhersagen: Sternhaufen im Feld der Milchstraße, Planetensysteme, Filamente und Wände, oberste gebundene Ebene
 5 Werte für die Prüfung von Hand
"""
import math, os, sys, io
import numpy as np
from scipy.optimize import brentq
from ftm_kern import (G, c, hbar, kB, m_p, m_e, alpha, r_e, a_B, R_p, M_sun, R_sun, M_earth, R_earth, AU, pc, kpc, Mpc, Gpc,
                      yr, Gyr, l_P, m_Pl, A0, Q_KOSM, q_KOSM, OMEGA_M, H0_si, ds, r0, r_stern, nu, chi, v_flach,
                      tau_D, schliesslaenge, zelle_aus_a0, rho_m_mitbewegt, compton_red, alpha_G, masse_r_stern_gleich_compton,
                      gewicht_w, realteil_quadrat_punkt)

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
out = []
P = lambda s="": out.append(s)
H0 = H0_si()
ell0 = c**2 / A0
LAM = schliesslaenge(A0)
TD = tau_D(A0)
eps_L = l_P / ell0
M_c = masse_r_stern_gleich_compton(A0)
M_U = c**4 / (G * A0)
M_X_STD = None

# ------------------------------------------------------------------ Hilfen
def aN_pm(M, r):
    return G * M / r**2

def a_pm(M, r):
    x = aN_pm(M, r)
    return float(nu(x)) * x

def gezeit_radial(a_fn, R, h=1e-4):
    """T = √(3/2)|E^TF| für ein kugelsymmetrisches Feld a(R) (Betrag, nach innen): Hesse diag(a', a/R, a/R)."""
    ap = (a_fn(R * (1 + h)) - a_fn(R * (1 - h))) / (2 * R * h)
    d = np.array([ap, a_fn(R) / R, a_fn(R) / R])
    E = d - d.sum() / 3
    return math.sqrt(1.5) * math.sqrt((E**2).sum())

M_MW = 6e10 * M_sun                             # Modellwert (03)
aMW = lambda R: a_pm(M_MW, R)                   # Milchstraße als Punktmasse mit der Regel (wie 03)
R_SONNE = 8.2 * kpc                             # [A] Sonnenkreis
T_SONNE_07 = A0 / (642 * pc)                    # T am Sonnenort aus 07 (R× = a₀/T = 642 pc)

def r_umkehr(M_b, eins_plus_delta):
    """Umkehrradius: mittlere Dichte (1+Δ)·ρ̄_m [A: Δ gegen die mittlere Materiedichte] mit der gezählten Masse (1+Q)·M_b."""
    return (3 * (1 + Q_KOSM) * M_b / (4 * math.pi * eins_plus_delta * rho_m_mitbewegt()))**(1 / 3)

def fmt(x, n=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:.{n}e}"

# ------------------------------------------------------------------ SPARC (Daten des Repos)
D_SP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "daten", "sparc", "SPARC_Lelli2016c.mrt")
sp = []
for l in open(D_SP, encoding="latin-1"):
    f = l.split()
    if len(f) < 18:
        continue
    try:
        d = dict(name=f[0], L=float(f[7]), Rd=float(f[11]), SB0=float(f[12]), MHI=float(f[13]), RHI=float(f[14]), Vf=float(f[15]), Q=int(f[17]))
    except ValueError:
        continue
    d['Mb'] = (0.5 * d['L'] + 1.33 * d['MHI']) * 1e9 * M_sun      # Υ_d = 0,5 (Modell), Helium 1,33 [A]
    sp.append(d)
spR = [d for d in sp if d['RHI'] > 0]
Mb_arr = np.array([d['Mb'] for d in spR]); R_arr = np.array([d['RHI'] * kpc for d in spR])
i_med = np.argsort(Mb_arr)[len(Mb_arr) // 2]
M_SP, R_SP = Mb_arr[i_med], R_arr[i_med]

# ------------------------------------------------------------------ 1 Ebenen
M_CELL_L = zelle_aus_a0(A0)
M_CELL = rho_m_mitbewegt() * M_CELL_L**3             # mit den Rückläufen gezählte Masse der Zelle
R_CELL = M_CELL_L * (3 / (4 * math.pi))**(1 / 3)     # Kugel gleichen Volumens

# (Name, M, R, Art, Elterngezeit T oder None, Wirbel: (v, ω, Träger))
def tide_mw(Rgc):
    return gezeit_radial(aMW, Rgc)

EB = []
EB.append(("Proton", m_p, R_p, "kompakt", T_SONNE_07, (c, c / compton_red(m_p), "Compton-Takt (QCD-Bindung)")))
EB.append(("H-Atom", m_p + m_e, a_B, "EM", T_SONNE_07, (alpha * c, alpha * c / a_B, "Hülle EM, v/c = α")))
EB.append(("Erde", M_earth, R_earth, "kompakt", T_SONNE_07, None))
EB.append(("Sonne", M_sun, R_sun, "kompakt", T_SONNE_07, None))
EB.append(("KSH dicht [A]", 1e5 * M_sun, 5 * pc, "kollektiv", tide_mw(10 * kpc), None))          # [A] 10⁵ M☉, 5 pc, R_gc 10 kpc
EB.append(("KSH dünn [A]", 2e4 * M_sun, 25 * pc, "kollektiv", tide_mw(70 * kpc), None))         # [A] 2·10⁴ M☉, 25 pc, R_gc 70 kpc
EB.append(("SPARC-Median", M_SP, R_SP, "kollektiv", None, None))
EB.append(("Milchstraße", M_MW, 15 * kpc, "kollektiv", None, None))                               # R 15 kpc [A]
EB.append(("Haufen", 1.7e14 * M_sun, 1 * Mpc, "kollektiv", None, None))                           # M aus 03, R 1 Mpc [A]
EB.append(("Kosmos-Zelle", M_CELL, R_CELL, "Zelle", None, (H0 * R_CELL, H0, "Atmen H₀")))

P("TRIADE AUF ALLEN EBENEN – Hintergrund / Wirbel / gefaltetes Zentrum")
P(f"Kern: a₀ = {A0:.2e} m/s², ℓ₀ = c²/a₀ = {ell0/Gpc:.2f} Gpc, Λ = {LAM/Gpc:.2f} Gpc, τ_D = c/a₀ = {TD/Gyr:.1f} Mrd. J., "
  f"H₀τ_D = {H0*TD:.3f}, Q = {Q_KOSM:.4f}, q = {q_KOSM:.4f}, Ω_m = {OMEGA_M:.4f}")
P(f"Elterngezeit am Sonnenort T = a₀/(642 pc) = {T_SONNE_07:.3e} s⁻² (aus 07); Milchstraße als Punktmasse 6·10¹⁰ M☉ mit der Regel")
P(f"SPARC ({len(sp)} Galaxien, {len(spR)} mit R_HI; Υ_d = 0,5, Gas 1,33·M_HI [A]): Median M_b = {M_SP/M_sun:.3e} M☉, R_HI = {R_SP/kpc:.2f} kpc")
P("Zelle: L = zelle_aus_a0(a₀); Masse ρ_m·L³ (mit Rückläufen gezählt), R = Radius der volumengleichen Kugel")
P()

rows = []
for (name, M, R, art, T, wir) in EB:
    D = ds(M); R0 = r0(M); rs = r_stern(M); gR = aN_pm(M, R)
    if name == "Kosmos-Zelle":
        gR = 4 * math.pi / 3 * G * rho_m_mitbewegt() * R
    kompakt = art in ("kompakt", "EM")
    nuR = 1.0 if (kompakt or art == "Zelle") else float(nu(gR))   # Zelle: w ≈ 1, im homogenen Hintergrund zählt keine Rückkehr (03)
    aR = nuR * gR
    if wir is None:
        v = math.sqrt(aR * R); om = v / R; trg = "Umlauf an R" + (" (Newton, geglättet)" if kompakt else " (mit Regel)")
    else:
        v, om, trg = wir
    ell = (G * M / T)**(1 / 3) if T else float('nan')
    dR = 3 * D / R
    d_eff = aR * R / c**2
    dstar = (v_flach(M) / c)**2
    RQ = rs / (-math.log(q_KOSM))
    rta = (r_umkehr(M, 10), r_umkehr(M, 12)) if art == "kollektiv" and M > 1e9 * M_sun else (float('nan'), float('nan'))
    w_orb = math.sqrt(G * M / R**3)
    rows.append(dict(name=name, M=M, R=R, art=art, D=D, R0=R0, rs=rs, gR=gR, nuR=nuR, v=v, om=om, trg=trg, ell=ell, dR=dR,
                     d_eff=d_eff, dstar=dstar, RQ=RQ, rta=rta, w_orb=w_orb, T=T))

P("1) DIE DREI GLIEDER JE EBENE (SI; Längen in m)")
P("   Zentrum: Δs = GM/3c², r₀ = 6Δs, Tiefe an der Hülle δ = 3Δs/R")
P(f"   {'Ebene':15s} {'M [kg]':>10s} {'Δs':>10s} {'r₀':>10s} {'δ(R)':>10s}")
for r in rows:
    P(f"   {r['name']:15s} {fmt(r['M'])} {fmt(r['D'])} {fmt(r['R0'])} {fmt(r['dR'])}")
P("   Wirbel: Hülle R, Neigung g(R), Umlauf v, ω an R (kompakt Newton; kollektiv mit ν(g(R)); Atom EM; Proton Compton-Takt; Zelle H₀)")
P(f"   {'Ebene':15s} {'R':>10s} {'g(R)/a₀':>10s} {'ν(R)':>7s} {'v [m/s]':>10s} {'ω [1/s]':>10s}  Träger")
for r in rows:
    P(f"   {r['name']:15s} {fmt(r['R'])} {fmt(r['gR']/A0)} {r['nuR']:7.3f} {fmt(r['v'])} {fmt(r['om'])}  {r['trg']}")
P("   Hintergrund und Übergänge: Schwellenradius r*, Glättungslänge ℓ = (GM/T)^{1/3}, Bilanzradius R_Q = r*/(−ln q), Umkehrradius r_ta (1+Δ = 10–12)")
P(f"   {'Ebene':15s} {'r*':>10s} {'ℓ':>10s} {'R_Q':>10s} {'r_ta(11)':>10s} {'r*/R':>10s} {'ℓ/r*':>10s} {'δ*=(v_f/c)²':>12s}")
for r in rows:
    P(f"   {r['name']:15s} {fmt(r['rs'])} {fmt(r['ell'])} {fmt(r['RQ'])} {fmt(0.5*(r['rta'][0]+r['rta'][1]))} {fmt(r['rs']/r['R'])} "
      f"{fmt(r['ell']/r['rs'])} {fmt(r['dstar'])}")
P("   Zustand: Tiefe δ_eff = a(R)·R/c² (Regel), Temperatur als Tiefe kT = f·δ_eff·μm_pc² (f = 0,5–0,85 aus 02/03, μ = 0,6 [A]);")
P("            Eigenfrequenz ω_orb = √(GM/R³), Übergabe 1/λ_m = N*/ω_orb (N* = 0,3–2,4) nur wo ξ in die Regel übergeht")
P(f"   {'Ebene':15s} {'δ_eff':>10s} {'T_δ(0,5) K':>11s} {'T_δ(0,85) K':>11s} {'ω_orb':>10s} {'1/λ_m [Mrd J]':>16s}  Gedächtnis")
mem = {"Proton": "Windungszahl n (Ladung), Chern-Zahl", "H-Atom": "Schalenschluss (Zählregel)",
       "Erde": "Knickspur (fest: Gitter als Identifikation)", "Sonne": "Zusammensetzung = größte Wärmetiefe δ_W",
       "KSH dicht [A]": "w = 0 (tiefster Kollaps)", "KSH dünn [A]": "w = 0 (tiefster Kollaps)",
       "SPARC-Median": "w = 0, ξ übergeben", "Milchstraße": "w = 0, ξ übergeben bis ~100 kpc",
       "Haufen": "w = 0, ξ nur zu 20–68 % übergeben", "Kosmos-Zelle": "Λ = ν_kosm·ΣΔs_b (gesamte Faltung, fest)"}
for r in rows:
    gas = r['art'] in ("kollektiv",) or r['name'] == "Sonne"
    T5 = 0.5 * r['d_eff'] * 0.6 * m_p * c**2 / kB if gas else float('nan')
    T85 = 0.85 * r['d_eff'] * 0.6 * m_p * c**2 / kB if gas else float('nan')
    lm = f"{0.3/r['w_orb']/Gyr:.3g}–{2.4/r['w_orb']/Gyr:.3g}" if r['art'] == "kollektiv" and r['M'] > 1e9 * M_sun else "–"
    P(f"   {r['name']:15s} {fmt(r['d_eff'])} {fmt(T5,2):>11s} {fmt(T85,2):>11s} {fmt(r['w_orb'])} {lm:>16s}  {mem[r['name']]}")
kTI = hbar * A0 / (c * kB)
P(f"   Zelle: Periode der inneren Richtung 2Λ/c = τ_D, als Temperatur geschrieben ħa₀/(ck) = {kTI:.3e} K (keine Wärme); Sonne: δ(R) ist Oberflächentiefe (Kern nach 02 höher)")
P(f"   Umkehrradius Milchstraße (1+Δ = 10–12): {rows[7]['rta'][0]/Mpc:.3f}–{rows[7]['rta'][1]/Mpc:.3f} Mpc (03: 0,58 Mpc); zwischen R_Q und r_ta ist der Übergang der Zählung offen (03)")
P()

# ------------------------------------------------------------------ 2 Verhältnisse
P("2) VERHÄLTNISSE ÜBER DIE EBENEN")
P("2a) Identitäten (auf jeder Ebene exakt gleich, jede Masse):")
P(f"   {'Ebene':15s} {'3Δs/(r*δ*)':>11s} {'r*/(ℓ₀δ*)':>11s} {'ν(r*)':>8s} {'Reτ(r*)':>8s} {'R_Q/r*':>8s} {'ω_N(r*)τ_D√δ*':>14s} {'r*²λ̄/ℓ_P²ℓ₀':>13s} {'δ_Bahn/(v/c)²':>13s}")
for r in rows:
    M = r['M']; rs = r['rs']; dst = r['dstar']
    omN = math.sqrt(G * M / rs**3)
    lam = compton_red(M)
    vN2 = G * M / r['R']
    P(f"   {r['name']:15s} {3*r['D']/(rs*dst):11.8f} {rs/(ell0*dst):11.8f} {float(nu(A0)):8.5f} {math.sqrt(realteil_quadrat_punkt(1.0)):8.5f} "
      f"{r['RQ']/rs:8.4f} {omN*TD*math.sqrt(dst):14.8f} {rs**2*lam/(l_P**2*ell0):13.8f} {(3*r['D']/r['R'])/(vN2/c**2):13.8f}")
P("   (3Δs : r* : ℓ₀ = δ*² : δ* : 1 – geometrische Leiter mit dem Schritt δ* = (v_flach/c)² auf jeder Ebene;")
P("    Atom in derselben Form: r_e : a_B = α² = (v/c)², Übergangslänge λ̄_e² = r_e·a_B wie r*² = 3Δs·ℓ₀)")
P(f"   Atom: r_e/a_B = {r_e/a_B:.6e}, α² = {alpha**2:.6e};  λ̄_e²/(r_e·a_B) = {compton_red(m_e)**2/(r_e*a_B):.8f}")
P()
P("2b) Skalengesetze (Form gleich, Zahl wandert mit M; Exponenten numerisch über 10⁻²⁷–10⁵³ kg):")
Ms = np.array([r['M'] for r in rows])
def slope(y):
    return np.polyfit(np.log(Ms), np.log(y), 1)[0]
P(f"   d ln r*/d ln M = {slope([r['rs'] for r in rows]):.4f} (½),  d ln v_flach/d ln M = {slope([v_flach(M) for M in Ms]):.4f} (¼),  "
  f"d ln δ*/d ln M = {slope([r['dstar'] for r in rows]):.4f} (½)")
P(f"   d ln ω_N(r*)/d ln M = {slope([math.sqrt(G*M/r_stern(M)**3) for M in Ms]):.4f} (−¼),  d ln ℓ/d ln M (feste Gezeit) = 1/3,  d ln r_ta/d ln M = 1/3,  "
  f"d ln (Δs·λ̄)/d ln M = {slope([ds(M)*compton_red(M) for M in Ms]):.2e} (0)")
P()
P("2c) Brüche: Größen, die zwischen den Ebenen springen")
comp = [r for r in rows if r['art'] in ("kompakt", "EM")]
coll = [r for r in rows if r['art'] == "kollektiv"]
for grp, nm in ((comp, "kompakt/EM"), (coll, "kollektiv")):
    Mg = np.array([r['M'] for r in grp]); Rg = np.array([r['R'] for r in grp])
    k = np.polyfit(np.log(Rg), np.log(Mg), 1)[0]
    P(f"   {nm:11s}: M ∝ R^k, k = {k:.2f};  r*/R = {min(r['rs']/r['R'] for r in grp):.3g} … {max(r['rs']/r['R'] for r in grp):.3g};  "
      f"g(R)/a₀ = {min(r['gR']/A0 for r in grp):.3g} … {max(r['gR']/A0 for r in grp):.3g};  δ(R) = {min(r['dR'] for r in grp):.2e} … {max(r['dR'] for r in grp):.2e}")
kS = np.polyfit(np.log(R_arr), np.log(Mb_arr), 1)[0]
gS = G * Mb_arr / R_arr**2 / A0
rsS = np.sqrt(G * Mb_arr / A0) / R_arr
P(f"   SPARC ({len(spR)} Galaxien, R = R_HI): M_b ∝ R^{kS:.2f};  g(R_HI)/a₀: Median {np.median(gS):.3f}, 16–84 % {np.percentile(gS,16):.3f}–{np.percentile(gS,84):.3f}, "
  f"Streuung {np.std(np.log10(gS)):.2f} dex über {np.log10(Mb_arr.max()/Mb_arr.min()):.1f} Dekaden in M_b")
P(f"                      r*/R_HI: Median {np.median(rsS):.3f}, 16–84 % {np.percentile(rsS,16):.3f}–{np.percentile(rsS,84):.3f}")
vf = np.array([d['Vf'] for d in spR]) * 1e3; okv = vf > 0
ratio = (vf[okv] / c)**2 / np.array([(v_flach(M) / c)**2 for M in Mb_arr[okv]])
P(f"   SPARC V_flat²/√(GM_b a₀): Median {np.median(ratio):.3f}, Streuung {np.std(np.log10(ratio)):.2f} dex ({okv.sum()} mit V_flat; enthält das Randglied ≈ 1 + r*/2r, keine Lesung von a₀)")
S0 = np.array([0.5 * d['SB0'] for d in sp if d['SB0'] > 0])
from ftm_kern import sigma_kollektiv
P(f"   SPARC zentrale Sternflächendichte Σ₀ = 0,5·SB_disk gegen a₀/(πG) = {sigma_kollektiv():.0f} M☉/pc²: Median {np.median(S0)/sigma_kollektiv():.3f}, "
  f"16–84 % {np.percentile(S0,16)/sigma_kollektiv():.3f}–{np.percentile(S0,84)/sigma_kollektiv():.3f}, Streuung {np.std(np.log10(S0)):.2f} dex")
P("   (R_HI ist eine Isophote (1 M☉/pc²): M ∝ R² dort teils durch die Definition; Σ₀ ist davon frei)")
P()
P("   Wirbeltakt gegen Rückkehrzeit: Umläufe je τ_D, ω·τ_D/2π (ω = Wirbelfrequenz aus Tabelle 1)")
for r in rows:
    P(f"   {r['name']:15s} {r['om']*TD/(2*math.pi):.3e}")
P(f"   Zelle: H₀τ_D/2π = {H0*TD/(2*math.pi):.4f} = 1/Θ (Θ = 2πa₀/cH₀ = {2*math.pi*A0/(c*H0):.4f}): ein Atemtakt je Rückkehr, heute; alle kleineren Ebenen laufen vielfach um")
P("   Zwei Ähnlichkeitsgruppen der Formeln (04): M ∝ R² hält g (und g/a₀), M ∝ R hält δ. Beide zugleich nur für λ = 1.")
P("   Kollektive Ebenen (Sternhaufen bis Zelle, 19 Dekaden in M) folgen näherungsweise der ersten: g(R)/a₀ = 0,01–5, r*/R = 0,1–2,3;")
P("   kompakte keiner (g(R)/a₀ = 10⁻⁷–10¹², r*/R = 6·10⁻⁴–1,6·10⁶). Selbstähnlich ist das Regime der kollektiven Ebenen, nicht die Tiefe δ.")
P()

# ------------------------------------------------------------------ 3 Unvollständige Triaden
P("3) UNVOLLSTÄNDIGE ODER ANDERE TRIADEN")
lam_e, lam_p = compton_red(m_e), compton_red(m_p)
P(f"   Schnittmasse M_c (r* = λ̄) = {M_c:.3e} kg;  Proton/M_c = {m_p/M_c:.2f}, Elektron/M_c = {m_e/M_c:.4f}")
P(f"   Elektron: r* = {r_stern(m_e):.3e} m, λ̄_e = {lam_e:.3e} m, r*/λ̄ = {r_stern(m_e)/lam_e:.3e}  → Schwelle innerhalb der eigenen Quantenlänge;")
P(f"             Spinzahl j = (1/2)/α_G = {0.5/alpha_G(m_e):.3e}: Verwindung ≫ Faltung, Triade nur in der inneren Richtung (I)")
P(f"   Proton:   r*/λ̄ = {r_stern(m_p)/lam_p:.1f}, r*/R_p = {r_stern(m_p)/R_p:.1f}, j = {0.5/alpha_G(m_p):.3e}; Wirbel (Quarks) ist QCD, Faltung α_G = {alpha_G(m_p):.3e}")
P(f"   H-Atom:   a_B/r*_p = {a_B/r_stern(m_p):.0f} (Hülle formal im Tiefbereich des Protons), Glättungslänge ℓ = {rows[1]['ell']*1e3:.2f} mm ≫ a_B → Newton;")
P(f"             Bindung EM: α/α_G = {alpha/alpha_G(m_p):.4e}; EM-Tiefe α² = {alpha**2:.3e} gegen Faltungstiefe 3Δs/a_B = {3*ds(m_p)/a_B:.3e}")
P(f"   Kompakte Ebenen (Erde, Sonne): r* weit außen (r*/R = {rows[2]['rs']/rows[2]['R']:.3g}, {rows[3]['rs']/rows[3]['R']:.3g}), ℓ/r* = "
  f"{rows[2]['ell']/rows[2]['rs']:.3g}, {rows[3]['ell']/rows[3]['rs']:.3g}: der Hintergrund ist nicht flach, sondern der Wirbel der Elternebene")
P(f"             Glättungskette in derselben Gezeit: ℓ_Erde = {rows[2]['ell']/AU:.0f} AE < ℓ_Sonne = {rows[3]['ell']/pc:.3f} pc; ℓ ∝ M^(1/3)")
MQ_R200 = (1 + Q_KOSM)**3 * math.log(q_KOSM)**6 * A0**3 / (1e4 * G * H0**4)
P(f"   Haufen:   R_Q/R = {rows[8]['RQ']/rows[8]['R']:.2f}; M_Q (R_Q = R₂₀₀) = {MQ_R200/M_sun:.3e} M☉ → Rückführung (Übergabe ξ) nicht ganz im Bereich")
cz = rows[-1]
P(f"   Kosmos:   M_Zelle/M_U = {cz['M']/M_U:.4f} (M_U = c⁴/Ga₀, oberste Schnittmasse: r* = 3Δs = ℓ₀);  Δs_Zelle/Λ = {cz['D']/LAM:.4f}")
P(f"             in Λ: 3Δs = {3*cz['D']/LAM:.3f}, ℓ₀ = {ell0/LAM:.3f}, r* = {cz['rs']/LAM:.3f} (√6 = {math.sqrt(6):.3f}), L = {M_CELL_L/LAM:.3f}, r₀ = {cz['R0']/LAM:.3f}")
P(f"             Zentrum, Wirbel, Schwelle fallen auf eine Länge (Faktor ≤ {max(3*cz['D'],ell0,cz['rs'],M_CELL_L)/min(3*cz['D'],ell0,cz['rs'],M_CELL_L):.2f}); L < r₀: kein Außenraum")
P(f"             g(R_Zelle)/a₀ = {cz['gR']/A0:.3f}; Randneigung der Hubble-Kugel cH₀/2 = {c*H0/2/A0:.3f}·a₀")
P()

# ------------------------------------------------------------------ 4 Vorhersagen
P("4) VORHERSAGEN AUS DER SELBSTÄHNLICHKEIT")
P("4a) Sternhaufen im Feld der Milchstraße (Punktmasse 6·10¹⁰ M☉ mit Regel). Eigene Triade vollständig (eigener Tiefbereich zwischen r* und")
P("    r_e = √(GM/g_e), innerhalb ℓ = (GM/T)^{1/3}) genau wenn g_e < a₀; Fenster r_e/r* = √(a₀/g_e), unabhängig von M.")
Rk = brentq(lambda R: aMW(R) - A0, 0.5 * kpc, 100 * kpc)
P(f"    g_e = a₀ bei R_gc = {Rk/kpc:.2f} kpc (a_N/a₀ dort {aN_pm(M_MW,Rk)/A0:.3f})")
P(f"    {'R_gc':>6s} {'g_e/a₀':>8s} {'r_e/r*':>7s} {'M×[M☉]':>10s} | M = 10⁴ M☉: {'r*[pc]':>7s} {'r_e[pc]':>8s} {'ℓ[pc]':>7s} | 10⁵: {'r*':>6s} {'r_e':>7s} {'ℓ':>7s} | 10⁶: {'r*':>6s} {'r_e':>7s} {'ℓ':>7s}")
for Rg in (5, 10, 20, 50, 100, 200):
    R = Rg * kpc; ge = aMW(R); T = gezeit_radial(aMW, R); Mx = A0**3 / (G * T**2)
    s = f"    {Rg:6.0f} {ge/A0:8.3f} {math.sqrt(A0/ge):7.3f} {Mx/M_sun:10.3e} |"
    for Mm in (1e4, 1e5, 1e6):
        M = Mm * M_sun; rs = r_stern(M); re = math.sqrt(G * M / ge); el = (G * M / T)**(1 / 3)
        s += f" {rs/pc:6.2f} {re/pc:7.2f} {el/pc:7.2f} |" if Mm > 1e4 else f"          {rs/pc:7.2f} {re/pc:8.2f} {el/pc:7.2f} |"
    P(s)
P("    Folge: Innerhalb R_gc ≈ {:.1f} kpc hat kein Sternhaufen einen eigenen Tiefbereich (Hintergrund = Newton-Feld der Galaxie);".format(Rk / kpc))
P("    außerhalb wächst das Fenster wie √(a₀/g_e). Das ist die Form, die 03 für Pal 14/NGC 2419 rechnet (dort gemessen weniger): keine neue Aussage.")
P()
P("4b) Planetensysteme (kompakte Zentren, Gezeit Sonnenort): Stern-Masse, r*, ℓ, ℓ/r*")
for Mm in (0.1, 1.0, 10.0):
    M = Mm * M_sun; rs = r_stern(M); el = (G * M / T_SONNE_07)**(1 / 3)
    P(f"    {Mm:5.1f} M☉: r* = {rs/AU:8.0f} AE, ℓ = {el/AU:8.0f} AE, ℓ/r* = {el/rs:6.2f}  (M× = {A0**3/(G*T_SONNE_07**2)/M_sun:.3e} M☉)")
P(f"    Mond um Erde: a_Mond/r*_Erde = {3.844e8/r_stern(M_earth):.2e} [A: 3,844·10⁸ m]; jedes Planetensystem liegt tief im Newton-Bereich seines Zentrums")
P("    und bis ℓ geglättet: kein eigener Tiefbereich, Hintergrund = Wirbel der Galaxie. Vorhersage: Newton bis ℓ (wie 03, weite Doppelsterne).")
P()
P("4c) Filamente und Wände (Flussgeometrie d = 3 Kugel, 2 Zylinder, 1 Platte):")
P("    a_N ∝ r^{−(d−1)}; mit der Regel d ln a/d ln r = −(d−1)(1 − χ/2) ⇒ Flussdimension d_eff = 1 + (d−1)(1 − χ/2)")
for d in (3, 2, 1):
    P(f"    d = {d}: Newton d_eff = {d}, tief (χ = 1) d_eff = {1+(d-1)*0.5:.1f}, verlorene Dimension (d−1)χ/2 = {(d-1)/2:.1f}")
P("    Numerisch (Zylinder, a_N = 2Gμ/r): d ln a/d ln r bei s² = a_N/a₀:")
for x in (100, 1, 0.01, 1e-4):
    ch = float(chi(x * A0))
    aa = lambda r_: float(nu(x * A0 / r_)) * x * A0 / r_
    sl = (math.log(aa(1.001)) - math.log(aa(0.999))) / (math.log(1.001) - math.log(0.999))
    P(f"      a_N/a₀ = {x:7.0e}: χ = {ch:.4f}, Steigung {sl:+.4f}, −(1 − χ/2) = {-(1-ch/2):+.4f}")
mu = 1e12 * M_sun / Mpc   # [A] baryonische Liniendichte
rs_f = 2 * G * mu / A0
P(f"    Schwelle des Zylinders r*_fil = 2Gμ/a₀; μ = 10¹² M☉/Mpc [A]: r*_fil = {rs_f/kpc:.2f} kpc; tief v_c² = a·r = √(2Gμa₀r) ∝ r^½:")
for rr in (0.1, 0.5, 1.0, 2.0):
    r_ = rr * Mpc; aN = 2 * G * mu / r_; a = float(nu(aN)) * aN
    P(f"      r = {rr:3.1f} Mpc: a_N/a₀ = {aN/A0:.3e}, ν = {float(nu(aN)):.2f}, v_c = {math.sqrt(a*r_)/1e3:6.1f} km/s (Newton {math.sqrt(aN*r_)/1e3:5.1f})")
thc = A0 / c
P(f"    Gewicht der Schließungsregel: w = ½ bei Ausdehnungsrate θ_c = a₀/c = {thc/H0:.4f}·H₀. Die Regel greift erst, wenn das Volumen langsamer atmet:")
for nm, n_ax in (("Filament (1 Achse mit H₀)", 1), ("Wand (2 Achsen mit H₀)", 2)):
    th_quer = (thc - n_ax * H0) / (3 - n_ax)
    P(f"      {nm}: Querrate je kollabierende Richtung ≤ {th_quer/H0:+.3f}·H₀ für w ≤ ½")
for th in (1.0, 0.5, 0.3, 0.168, 0.1, 0.0):
    P(f"      θ = {th:5.3f}·H₀: w = {gewicht_w(th*H0):.4f}")
P("    Folge: Filamente mit Querschnitt, der langsamer als 0,42·H₀ schrumpft (bei axialer Dehnung H₀), sind mitatmend (w > ½):")
P("    sie tragen die Mode ξ wie Staub mit 1+Q (Standardform); erst danach gilt quer das Zylindergesetz d_eff = 2 − χ/2.")
P()
P("4d) Oberste Ebene mit vollständiger Triade:")
for d1 in (10, 12):
    # R_Q = r_ta
    fn = lambda lM: r_stern(10**lM * M_sun) / (-math.log(q_KOSM)) - r_umkehr(10**lM * M_sun, d1)
    lM1 = brentq(fn, 8, 25)
    fn2 = lambda lM: r_stern(10**lM * M_sun) - r_umkehr(10**lM * M_sun, d1)
    lM2 = brentq(fn2, 8, 30)
    P(f"    1+Δ = {d1}: R_Q = r_ta bei M_b = {10**lM1:.3e} M☉ (M = (1+Q)M_b = {(1+Q_KOSM)*10**lM1:.3e}); r* = r_ta bei M_b = {10**lM2:.3e} M☉")
P(f"    M_Q (R_Q = R₂₀₀, 03) = {MQ_R200/M_sun:.3e} M☉ (gesamt);  Baryonen der Zelle ΣM_b = {M_CELL/(1+Q_KOSM)/M_sun:.3e} M☉")
P("    Folge: Jede gebundene Masse schließt ihre Rückführung (R_Q) innerhalb ihres Umkehrradius; nur über M_Q liegt R_Q außerhalb des")
P(f"    virialisierten Teils (Haufen). r* = r_ta erst bei M_b ≈ 5–7·10²¹ M☉ (r* = {r_stern(6e21*M_sun)/Gpc:.1f} Gpc), weit über jeder gebundenen Struktur;")
P("    nur die Zelle hat ihre Schwelle am eigenen Rand (r*/L = {:.3f}).".format(cz['rs'] / M_CELL_L))
P()

# ------------------------------------------------------------------ 5 Handprüfung
P("5) WERTE FÜR DIE PRÜFUNG VON HAND")
P(f"   r*(Proton) = √(G m_p/a₀) = {r_stern(m_p):.6e} m")
P(f"   ν(r*) = 1/(1 − e⁻¹) = {float(nu(A0)):.6f}")
P(f"   M_Zelle/M_U = {cz['M']/M_U:.6f};  r*_Zelle/Λ = {cz['rs']/LAM:.6f};  H₀τ_D = {H0*TD:.6f};  θ_c/H₀ = {thc/H0:.6f}")
P(f"   R_gc(g_e = a₀) = {Rk/kpc:.4f} kpc;  M_Q = {MQ_R200/M_sun:.4e} M☉;  M_b(R_Q = r_ta, 1+Δ = 10) = {brentq(lambda lM: r_stern(10**lM*M_sun)/(-math.log(q_KOSM)) - r_umkehr(10**lM*M_sun, 10), 8, 25):.4f} (log₁₀ M☉)")
P("   Von Hand (unabhängig, Taschenrechner-Weg): r*_p = √(6,6743e-11·1,67262e-27/1,1e-10) = √(1,01487e-27) = 3,18570e-14 m;")
P("   ν(r*) = 1/(1 − 0,367879) = 1,581977;  M_Zelle = 3c⁴/(2Ga₀) aus Λ = Gρ_mL³/3c² und a₀ = c²/2Λ ⇒ 1,5·M_U, r* = √(3/2)·2Λ = √6·Λ;")
P("   R_gc: x = 0,511, x/(1 − e^{−√x}) = 1,0006, R = r*_MW/√x = 2,6905e20/0,7148 m = 12,20 kpc;  M_Q = 257,8·2,49e-5·1,331e-30/(6,674e-7·2,276e-71) kg = 2,83e14 M☉;")
P("   R_Q = r_ta (1+Δ = 10): GM_b = a₀³ln⁶q·(2(1+Q)/((1+Δ)Ω_mH₀²))² = 1,331e-30·2,4948e-5·16,49/2,2763e-71 ⇒ M_b = 1,813e17 M☉")

txt = "\n".join(out)
print(txt)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "triade_ebenen.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
