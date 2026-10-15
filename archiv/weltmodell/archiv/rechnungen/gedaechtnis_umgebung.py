# -*- coding: utf-8 -*-
"""
Gedächtnis der Rückkehr im zeitlich wechselnden äußeren Feld: Ordnet es die Abweichungen der druckgestützten Systeme?

1  Herleitung (07): Kopfgleichung r̈ = −∇Φ_bar − χκ₀∫Δ(τ_M)(e_r⊗e_r)e^{−λ_K(t−t′)}ṙ(t′)dt′.
   Der Kern trägt Δ ≈ q_Δs, den Nome; mit Im τ_M = s/2π (Punkt 1) ist |q| = e^{−s} der Rückkehranteil q.
   07: „Die Bindung läuft nur über den Betrag des Nomes“, „das Gedächtnis wirkt auf den Zustand der Regel, nicht auf die
   einzelne Bahn“, die Rückkehr ist „eine stationäre Schleife am festen Ort“, die den Betrag des Gefälles dort sieht.
   Der Zustand der Regel am Ort x ist q(x) = e^{−s(x)}, s² = |g_i(x) + g_e(t)|/a₀. Mit Gedächtnis:
       q_K(x,t) = λ_K ∫_{−∞}^{t} e^{−λ_K(t−t′)} q(|g_i(x) + g_e(t′)|) dt′,     ν_eff(x,t) = 1/(1 − q_K(x,t)),
       a(x,t) = ν_eff(x,t)·(g_i(x) + g_e(t))                                     (die Regel wirkt auf das heutige Feld)
   λ_K → ∞: q_K = q(heute), örtliche Regel.  λ_K → 0: q_K = Zeitmittel von q über die Bahn.
   Varianten zum Vergleich: Gedächtnis auf s (ν = 1/(1 − e^{−⟨s⟩})) und auf die Verlustrate Γ ∝ |g| (s² = ⟨|g|⟩/a₀).
   Lesart A (Hauptlesart): Ort im System, Richtung des äußeren Felds fest (heutige), nur sein Betrag g_e(t′) wechselt
            (Bezug mit der Linie zum Wirt mitgeführt). Kreisbahn: g_e konstant → unabhängig von λ_K.
   Lesart B: Ort im nicht drehenden Bezug des Systems; die Richtung zum Wirt dreht mit der Bahn (Rosette, Vorrücken Δφ je T_r).
2  Bahnen im Wirt: Milchstraße 6·10¹⁰ M☉ als Punktmasse mit der Regel bzw. logarithmisches Potential v = 200 km/s [A];
   NGC 1052 10¹¹ M☉ als Punktmasse mit der Regel [A]. Das äußere Feld der Regel ist das Newton-Feld g_e = GM_h/r².
3  σ_los: Plummer (b = R_e), isotrope Jeans-Gleichung, Fenster R < R_e, Monopol des Felds (wie glaettung_stetig.py).
   [A] σ im Gleichgewicht mit der heutigen wirksamen Regel.
4  λ_K in s⁻¹ gegen die Grenze der alten Scheiben λ_K ≲ 3·10⁻⁴ Ω_Sonne (kern_gedaechtnis.py).
5  Prüfungen: λ_K → ∞ gegen die bisherigen Werte; λ_K → 0 von Hand (Bahnzeit mit quad statt Integrator, Quasi-Newton-Formel
   σ² = σ_N²·ν̄ bzw. σ_N²·ν(1 + K/3)); Kreisbahn; Schranken σ_N ≤ σ ≤ max über konstantes g_e.
[A] Annahme, [G] Messwert aus dem Gedächtnis, ungeprüft.
"""
import io, sys, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from ftm_kern import *

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
np.seterr(all="ignore")
out = []
def P(s=""): out.append(s); print(s)
def H(t): P(); P("=" * 110); P(t); P("=" * 110)

# ------------------------------------------------------------------ Einheiten: kpc, km/s; Felder in a₀
A0K = A0 * kpc / 1e6                     # a₀ in (km/s)²/kpc
GK = G * M_sun / kpc / 1e6               # G in kpc (km/s)²/M☉
TU = kpc / 1e3                           # Zeiteinheit kpc/(km/s) in s
GMu = lambda M: GK * M / A0K             # GM in kpc²·a₀
kms = lambda s2: np.sqrt(np.maximum(s2, 0) * A0K)
def qf(x): return np.exp(-np.sqrt(np.maximum(x, 0.0)))          # q(|g|/a₀)
def nuq(q): return 1.0 / (1.0 - q)

# ------------------------------------------------------------------ Systeme (Eingaben wie glaettung_stetig.py)
SYS = [  # Name, M [M☉], R_e [kpc], Wirt, heutiger Abstand [kpc], σ gemessen, Fehler, Richtung verlangt
    ("Pal 14",    2e4,   0.0247, "MW", 71.0,  0.38, 0.12, "kleiner"),
    ("NGC 2419",  9e5,   0.019,  "MW", 90.0,  4.1,  0.5,  "kleiner"),
    ("Crater II", 3.2e5, 1.066,  "MW", 117.5, 2.7,  0.3,  "größer"),
    ("Draco",     5.8e5, 0.221,  "MW", 76.0,  9.1,  1.2,  "größer"),
    ("Fornax",    4.0e7, 0.71,   "MW", 147.0, 11.7, 0.9,  "kleiner"),
    ("DF2",       2.0e8, 2.2,    "N1052", 80.0, 8.5, 2.3, "kleiner"),
    ("DF4",       1.5e8, 1.6,    "N1052", 165.0, 4.2, 2.2, "kleiner")]
REF = {"Pal 14": 1.78, "NGC 2419": 6.49, "Crater II": 1.39, "Draco": 3.33, "Fornax": 13.96, "DF2": 17.70, "DF4": 19.21}
MHOST = {"MW": 6e10, "N1052": 1e11}
# Bahnen [G] (Gaia-Bahnen, aus dem Gedächtnis, grob; Bereich der Literatur als zwei Szenarien)
ORB = {"Pal 14": [(5, 90), (30, 90)], "NGC 2419": [(15, 92), (25, 92)], "Crater II": [(20, 135), (40, 135)],
       "Draco": [(40, 95), (60, 95)], "Fornax": [(60, 155), (100, 155)]}
LAMS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, np.inf]
OM_SUN = 220.0 / 8.0 / TU                # s⁻¹ (wie kern_gedaechtnis.py)
LAM_DISK = 3.07e-4 * OM_SUN              # Grenze der alten Scheiben (e = 0,3)

# ------------------------------------------------------------------ Plummer, Jeans
def plummer(M, Re, n=1500):
    b = Re; GM = GMu(M); r = np.logspace(math.log10(b) - 4, math.log10(b) + 4.5, n)
    GMenc = GM * r**3 / (r**2 + b**2)**1.5; Grho = 3 * GM / (4 * np.pi * b**3) * (1 + (r / b)**2)**-2.5
    return r, GMenc / r**2, Grho
def jeans(r, rho, gr, Re):
    f = rho * gr
    p = np.concatenate([np.cumsum((0.5 * (f[1:] + f[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    w = np.where(r <= Re, 1.0, 1 - np.sqrt(np.clip(1 - (Re / r)**2, 0, 1)))
    return float(kms(np.trapezoid(r**2 * p * w, r) / np.trapezoid(r**2 * rho * w, r)))

MU, WMU = np.polynomial.legendre.leggauss(96)
def geff_A(gi, ge_hist, W, ge0, mode="q"):
    """Lesart A. gi: (nr,), ge_hist: (nt,) Beträge der Geschichte, W: (nl, nt) Gewichte (Summe 1), ge0: heute.
       Rückgabe (nl, nr): Monopol des einwärts gerichteten Felds ⟨ν_eff·(g_i + g_e0 μ)⟩_μ (Richtung fest, μ = cos)."""
    res = np.empty((W.shape[0], gi.size))
    for j0 in range(0, gi.size, 100):
        g = gi[j0:j0 + 100, None, None]
        x = np.sqrt(g**2 + ge_hist[None, None, :]**2 + 2 * g * ge_hist[None, None, :] * MU[None, :, None])  # |g| in a₀
        if mode == "q":   Z = qf(x)
        elif mode == "s": Z = np.sqrt(x)
        else:             Z = x
        Zk = np.einsum("rmt,lt->lrm", Z, W)
        if mode == "q":   nue = nuq(Zk)
        elif mode == "s": nue = 1 / (-np.expm1(-Zk))
        else:             nue = 1 / (-np.expm1(-np.sqrt(Zk)))
        fld = g[:, :, 0] + ge0 * MU[None, :]
        res[:, j0:j0 + 100] = 0.5 * np.sum(WMU[None, None, :] * nue * fld[None], axis=2)
    return res

def geffA_fine(r, gi, geh, W, ge0, mode):
    """geff_A auf 150 Stützradien, Verhältnis g_eff/g_i in ln r auf das feine Gitter interpoliert"""
    ic = np.linspace(0, len(r) - 1, 150).astype(int)
    gc = geff_A(gi[ic], geh, W, ge0, mode)
    return np.array([gi * np.interp(np.log(r), np.log(r[ic]), g_ / gi[ic]) for g_ in gc])

def sigma_const(M, Re, ge):
    r, gi, rho = plummer(M, Re)
    x = np.sqrt(gi[:, None]**2 + ge**2 + 2 * gi[:, None] * ge * MU[None, :])
    g = 0.5 * np.sum(WMU * nuq(qf(x)) * (gi[:, None] + ge * MU[None, :]), axis=1)
    return jeans(r, rho, g, Re)

# ------------------------------------------------------------------ Bahnen
def host_acc(host, pot):
    GMh = GK * MHOST[host]
    if pot == "log": return lambda r: 200.0**2 / r
    return lambda r: float(nu(GMh / r**2 / A0K, 1.0)) * GMh / r**2
def host_phi(acc, r, rref=1.0):
    return quad(acc, rref, r, limit=200)[0]
def orbit(host, pot, rp, ra, n=801):
    """Bahn von Peri nach Peri (eine Radialperiode): t, r, φ. E, L aus rp, ra."""
    acc = host_acc(host, pot)
    if abs(ra - rp) < 1e-9:
        v = math.sqrt(acc(rp) * rp); T = 2 * math.pi * rp / v
        t = np.linspace(0, T, n); return t, np.full(n, rp), v / rp * t, T, 2 * math.pi
    dPhi = host_phi(acc, ra, rp)
    L = math.sqrt(2 * dPhi / (1 / rp**2 - 1 / ra**2))
    def f(t, y):
        r, vr, ph = y; return [vr, L**2 / r**3 - acc(r), L / r**2]
    ev = lambda t, y: y[1]; ev.direction = 1; ev.terminal = False
    Tg = 4 * math.pi * ra / (L / ra + 1e-9)
    sol = solve_ivp(f, (0, 3 * Tg), [rp * (1 + 1e-12), 0.0, 0.0], rtol=1e-10, atol=1e-12, events=ev, dense_output=True, max_step=Tg / 400)
    Tr = [te for te in sol.t_events[0] if te > 1e-3 * Tg][0]
    t = np.linspace(0, Tr, n); y = sol.sol(t)
    return t, y[0], y[2], Tr, y[2][-1]

def history(t, r, ph, T, D, outbound=True):
    """Zeitpunkt heute: r = D (auswärts oder einwärts). Rückgabe τ (Rückblick), r(t0 − τ), φ relativ zu heute."""
    n = len(t)
    if np.ptp(r) < 1e-9 * r[0]: i0 = 0
    else:
        half = np.argmax(r)
        idx = np.arange(0, half + 1) if outbound else np.arange(half, n)
        i0 = idx[np.argmin(np.abs(r[idx] - D))]
    tau = t[:n - 1]
    j = (i0 - np.arange(n - 1)) % (n - 1)
    back = (i0 - np.arange(n - 1)) < 0
    return tau, r[j], ph[j] - ph[i0] - np.where(back, ph[-1], 0.0)

def weights(tau, T, lam_over_Om):
    """Gewichte der Geschichte über eine Radialperiode, periodisch fortgesetzt (exakt): λe^{−λτ}/(1 − e^{−λT})."""
    Om = 2 * math.pi / T; dt = tau[1] - tau[0]
    W = []
    for lo in lam_over_Om:
        if np.isinf(lo):
            w = np.zeros_like(tau); w[0] = 1.0
        else:
            lam = lo * Om; w = np.exp(-lam * tau)
        W.append(w / w.sum())
    return np.array(W)

# ============================================================================ 1
H("1   Wirksame Verstärkung mit Gedächtnis (aus der Kopfgleichung von 07)")
P("   Kern: Δ ≈ q_Δs, Nome; Punkt 1: |q| = e^{−s} = Rückkehranteil. Gedächtnis auf den Zustand der Regel am Ort:")
P("     q_K(x,t) = λ_K∫_{−∞}^t e^{−λ_K(t−t′)} exp(−√(|g_i(x) + g_e(t′)|/a₀)) dt′,   ν_eff = 1/(1 − q_K),   a = ν_eff·(g_i + g_e(t))")
P("   periodische Bahn (Radialperiode T): Gewicht der Geschichte w(τ) = λ_K e^{−λ_Kτ}/(1 − e^{−λ_K T}), 0 ≤ τ < T (exakt).")
P("   Grenzfälle: λ_K → ∞ örtliche Regel (heutiges g_e); λ_K → 0 Zeitmittel ⟨q⟩_T; Kreisbahn (Lesart A) unabhängig von λ_K.")
P("   Schwaches Eigenfeld (g_i ≪ g_e), Lesart A: g_eff/g_i = ν̄ + (ν̄²g_e0/3)⟨q′(g_e(t′))⟩_K, ν̄ = 1/(1 − ⟨q(g_e)⟩_K);")
P("   Lesart B: dasselbe mit ⟨q′·cos(φ(t′) − φ₀)⟩_K; für λ_K → 0 fällt das zweite Glied weg (Richtungen gleichverteilt).")

# ============================================================================ 2 + 3
H("2/3 Bahnen, g_e(t) und σ_los (Lesart A, Gedächtnis auf q); λ_K/Ω_Bahn, Ω_Bahn = 2π/T_r")
P(f"   a₀ = {A0:.2e} m/s² = {A0K:.1f} (km/s)²/kpc;  Grenze der Scheiben λ_K ≲ 3,07·10⁻⁴·Ω_Sonne = {LAM_DISK:.3e} s⁻¹ (1/λ = {1/LAM_DISK/Gyr:.0f} Mrd. Jahre)")
RES = {}          # (name, szenario) -> dict
SCN = {}
for nm, M, Re, host, D, so, er, want in SYS:
    GMh = GMu(MHOST[host])
    r, gi, rho = plummer(M, Re)
    sN = jeans(r, rho, gi, Re)
    ge_now = GMh / D**2
    s_inf = sigma_const(M, Re, ge_now)
    # Schranken über konstantes g_e
    gg = np.logspace(-5, 2, 57); ss = np.array([sigma_const(M, Re, x) for x in gg]); s_iso = sigma_const(M, Re, 0.0)
    smax = max(ss.max(), s_iso)
    P()
    P(f"-- {nm}: M = {M:.2g} M☉, R_e = {Re} kpc, heute D = {D} kpc, g_e heute = {ge_now:.4f} a₀;  σ_N = {sN:.2f}, σ(λ=∞) = {s_inf:.2f} "
      f"(bisher {REF[nm]}), isoliert {s_iso:.2f}, max über konst. g_e {smax:.2f} km/s;  gemessen {so} ± {er} [G] ({want})")
    if so < sN: ger = None; req = "keines (gemessen unter Newton)"
    elif so > smax: ger = None; req = "keines (gemessen über dem größten Wert)"
    else:
        ger = brentq(lambda x: sigma_const(M, Re, x) - so, 1e-6, 1e2, xtol=1e-9)
        req = f"g_e = {ger:.4f} a₀ (Abstand {math.sqrt(GMh / ger):.0f} kpc)"
    P(f"   verlangtes konstantes äußeres Feld: {req}")
    scen = []
    if host == "MW":
        for pot in ("punkt", "log"):
            for rp, ra in ORB[nm]:
                for ob in (True, False):
                    scen.append((pot, rp, ra, ob))
    else:
        for e in (0.0, 0.5, 0.8):
            for pl in (("peri", "apo", "mitte") if e > 0 else ("kreis",)):
                if pl == "kreis": rp = ra = D
                elif pl == "peri": rp, ra = D, D * (1 + e) / (1 - e)
                elif pl == "apo": rp, ra = D * (1 - e) / (1 + e), D
                else: a_ = D; rp, ra = a_ * (1 - e), a_ * (1 + e)
                scen.append(("punkt", rp, ra, True, e, pl))
    P(f"   {'Bahn':32s} {'T_r[Gyr]':>8s} {'g_e peri/heute/apo [a₀]':>24s} {'Ω_Bahn[s⁻¹]':>11s} {'λ_Scheibe/Ω':>11s}  "
      + "  ".join(f"{('λ/Ω=' + ('∞' if np.isinf(l) else f'{l:g}')):>12s}" for l in LAMS))
    for sc in scen:
        pot, rp, ra, ob = sc[:4]
        if rp > D + 1e-9 or ra < D - 1e-9: continue
        t, rr, ph, T, dphi = orbit(host, pot, rp, ra)
        tau, rh, phh = history(t, rr, ph, T, D, ob)
        geh = GMh / rh**2
        W = weights(tau, T, LAMS)
        gq = geffA_fine(r, gi, geh, W, ge_now, "q")
        sig = [jeans(r, rho, gq[k], Re) for k in range(len(LAMS))]
        Om = 2 * math.pi / (T * TU)
        lab = (f"{pot}, {rp:g}–{ra:g} kpc, {'aus' if ob else 'ein'}" if host == "MW"
               else f"e = {sc[4]:.1f}, heute {sc[5]}, {rp:.0f}–{ra:.0f}")
        RES[(nm, lab)] = dict(sig=sig, T=T, Om=Om, geh=geh, W=W, gi=gi, r=r, rho=rho, tau=tau, phh=phh, dphi=dphi, sc=sc)
        P(f"   {lab:32s} {T*TU/Gyr:8.2f} {GMh/rp**2:7.4f}/{ge_now:7.4f}/{GMh/ra**2:7.4f} {Om:11.3e} {LAM_DISK/Om:11.2e}  "
          + "  ".join(f"{s:6.2f} ({(so/s)**2:4.2f})" for s in sig))
    SCN[nm] = dict(sN=sN, s_inf=s_inf, s_iso=s_iso, smax=smax, ger=ger, so=so, er=er, want=want, M=M, Re=Re, ge_now=ge_now)
P()
P("   Einträge: σ_los [km/s] und in Klammern (σ_gem/σ)².")

# ============================================================================ Varianten
H("Varianten des Gedächtnisses (Lesart A): auf q / auf s / auf Γ ∝ |g|; σ [km/s] bei λ/Ω = 10⁻² und 1")
for (nm, lab), d in RES.items():
    if not (lab.startswith("log") or lab.startswith("e = 0.8")): continue
    M, Re = SCN[nm]["M"], SCN[nm]["Re"]
    row = []
    for mode in ("q", "s", "g"):
        gg_ = geffA_fine(d["r"], d["gi"], d["geh"], d["W"][[2, 4]], SCN[nm]["ge_now"], mode)
        row.append([jeans(d["r"], d["rho"], gg_[k], Re) for k in range(2)])
    P(f"   {nm:10s} {lab:32s}  q: {row[0][0]:6.2f} {row[0][1]:6.2f}   s: {row[1][0]:6.2f} {row[1][1]:6.2f}   Γ: {row[2][0]:6.2f} {row[2][1]:6.2f}")

# ============================================================================ Lesart B
H("Lesart B (Richtung zum Wirt dreht mit der Bahn, nicht drehender Bezug des Systems), Gedächtnis auf q")
TH, WTH = np.polynomial.legendre.leggauss(16)          # cos θ gegen die Bahnnormale
PSI = np.linspace(0, 2 * np.pi, 32, endpoint=False)
def geff_B(gi_c, ge_t, phi_t, Wt, ge0):
    """Monopol auf der Kugel: n·u(t′) = sinθ cos(ψ − φ(t′)); heute φ = 0. Wt: (nl, nt)."""
    st = np.sqrt(1 - TH**2)
    cosang = st[:, None, None] * np.cos(PSI[None, :, None] - phi_t[None, None, :])     # (θ, ψ, t)
    cos0 = st[:, None] * np.cos(PSI[None, :])
    res = np.empty((Wt.shape[0], gi_c.size))
    for i, g in enumerate(gi_c):
        x = np.sqrt(g**2 + ge_t[None, None, :]**2 - 2 * g * ge_t[None, None, :] * cosang)
        qk = np.einsum("abt,lt->lab", qf(x), Wt)
        rad = g - ge0 * cos0                                                  # einwärts: −n·(g_i + g_e0)
        res[:, i] = np.einsum("a,lab->l", WTH, nuq(qk) * rad[None]) / (2 * len(PSI))
    return res
def hist_B(d, lam_over_Om, nt=60, nk_max=40000):
    """Geschichte über viele Radialperioden: φ(t0 − τ − kT) = φ(t0 − τ) − kΔφ, Gewicht e^{−λ(τ + kT)}; auf nt Stützstellen je Periode."""
    tau, geh, phh, T, dphi = d["tau"], d["geh"], d["phh"], d["T"], d["dphi"]
    sel = np.linspace(0, len(tau) - 1, nt).astype(int)
    tau_s, ge_s, ph_s = tau[sel], geh[sel], phh[sel]
    Om = 2 * math.pi / T
    ges, phs, ws = [], [], []
    for lo in lam_over_Om:
        if np.isinf(lo):
            ges.append(np.array([ge_s[0]])); phs.append(np.array([0.0])); ws.append(np.array([1.0])); continue
        lam = lo * Om
        K = int(min(nk_max, max(1, math.ceil(12 / (lam * T)))))
        k = np.arange(K)
        w = np.exp(-lam * (tau_s[None, :] + k[:, None] * T))
        if K == nk_max:                                                      # λ → 0: Richtungen gleichverteilt
            pass
        ph = ph_s[None, :] - k[:, None] * dphi
        # auf 36 Winkelklassen zusammenfassen
        nb = 36; b = np.floor((ph % (2 * np.pi)) / (2 * np.pi) * nb).astype(int)
        Wb = np.zeros((nb, nt)); np.add.at(Wb, (b, np.broadcast_to(np.arange(nt), b.shape)), w)
        Wb /= Wb.sum()
        cen = (np.arange(nb) + 0.5) / nb * 2 * np.pi
        ges.append(np.tile(ge_s, nb)); phs.append(np.repeat(cen, nt)); ws.append(Wb.ravel())
    return ges, phs, ws
def sigma_B(d, M, Re, ge0, lams):
    r, rho = d["r"], d["rho"]; gi = d["gi"]
    ic = np.linspace(0, len(r) - 1, 60).astype(int)
    ges, phs, ws = hist_B(d, lams)
    out_ = []
    for ge_t, ph_t, w in zip(ges, phs, ws):
        m = w > 1e-12 * w.max()
        gc = geff_B(gi[ic], ge_t[m], ph_t[m], (w[m] / w[m].sum())[None, :], ge0)[0]
        ratio = np.interp(np.log(r), np.log(r[ic]), gc / gi[ic])
        out_.append(jeans(r, rho, gi * ratio, Re))
    return out_
BSEL = {"Pal 14": "log, 5–90 kpc, aus", "NGC 2419": "log, 15–92 kpc, aus", "Crater II": "log, 20–135 kpc, aus",
        "Draco": "log, 40–95 kpc, aus", "Fornax": "log, 60–155 kpc, aus", "DF2": "e = 0.5, heute apo, 27–80", "DF4": "e = 0.5, heute apo, 55–165"}
P("   " + f"{'System, Bahn':44s}" + "  ".join(f"{('λ/Ω=' + ('∞' if np.isinf(l) else f'{l:g}')):>12s}" for l in LAMS))
for nm, lab in BSEL.items():
    d = RES.get((nm, lab))
    if d is None: continue
    sb = sigma_B(d, SCN[nm]["M"], SCN[nm]["Re"], SCN[nm]["ge_now"], LAMS)
    so = SCN[nm]["so"]
    P(f"   {nm + ', ' + lab:44s}" + "  ".join(f"{s:6.2f} ({(so/s)**2:4.2f})" for s in sb))
P("   Kreisbahn in Lesart B (Crater II, DF4 auf heutigem Abstand):")
for nm in ("Crater II", "DF4"):
    lab = [k[1] for k in RES if k[0] == nm and ("kreis" in k[1])]
    if not lab:
        host = "MW"; D = [s for s in SYS if s[0] == nm][0][4]
        t, rr, ph, T, dphi = orbit(host, "log", D, D)
        tau, rh, phh = history(t, rr, ph, T, D, True)
        r_, gi_, rho_ = plummer(SCN[nm]["M"], SCN[nm]["Re"])
        d = dict(tau=tau, geh=GMu(MHOST[host]) / rh**2, phh=phh, T=T, dphi=dphi, r=r_, gi=gi_, rho=rho_)
    else: d = RES[(nm, lab[0])]
    sb = sigma_B(d, SCN[nm]["M"], SCN[nm]["Re"], SCN[nm]["ge_now"], LAMS)
    P(f"      {nm:10s}: " + "  ".join(f"{s:6.2f}" for s in sb) + "   (Lesart A: für alle λ gleich σ(λ=∞))")

# ============================================================================ 4
H("4   Ordnung der Vorzeichen: Richtung der Änderung gegen λ = ∞ und Erreichbarkeit")
P("   verlangt: 'kleiner' → Gedächtnis muss σ senken (wirksames g_e größer als heute), 'größer' → heben.")
P(f"   {'System':10s} {'verl.':>7s} {'σ_N':>6s} {'σ(∞)':>6s} {'max':>6s} {'gem.':>9s}   Spanne σ über alle Bahnen je λ/Ω (min–max)")
for nm, M, Re, host, D, so, er, want in SYS:
    ks = [k for k in RES if k[0] == nm]
    rng = []
    for il in range(len(LAMS)):
        v = [RES[k]["sig"][il] for k in ks]; rng.append((min(v), max(v)))
    sc = SCN[nm]
    P(f"   {nm:10s} {want:>7s} {sc['sN']:6.2f} {sc['s_inf']:6.2f} {sc['smax']:6.2f} {so:5.2f}±{er:<4.2f} "
      + "  ".join(f"{a:5.2f}–{b:<5.2f}" for a, b in rng))
P()
P("   Je λ/Ω: Zahl der Systeme, bei denen irgendeine Bahn in die verlangte Richtung ändert / den Messwert (±1 Fehler) erreicht:")
for il, l in enumerate(LAMS[:-1]):
    ndir, nhit, names = 0, 0, []
    for nm, M, Re, host, D, so, er, want in SYS:
        ks = [k for k in RES if k[0] == nm]
        v = np.array([RES[k]["sig"][il] for k in ks]); s0 = SCN[nm]["s_inf"]
        okdir = np.any(v < s0 - 1e-3) if want == "kleiner" else np.any(v > s0 + 1e-3)
        hit = np.any(np.abs(v - so) <= er)
        ndir += okdir; nhit += hit; names.append(f"{nm}:{'R' if okdir else '-'}{'T' if hit else '-'}")
    P(f"   λ/Ω = {l:g}: Richtung {ndir}/7, getroffen {nhit}/7   " + " ".join(names))
P("   (R = Richtung stimmt bei mindestens einer Bahn, T = Messwert ±1σ bei mindestens einer Bahn getroffen)")

# ============================================================================ 5
H("5   Unabhängige Prüfung")
P("   a) λ_K → ∞ gegen die bisherigen Werte (glaettung_stetig.py, a₀ = 1,1·10⁻¹⁰): " +
  ", ".join(f"{nm} {SCN[nm]['s_inf']:.2f}/{REF[nm]}" for nm in REF))
# b) Crater II, λ → 0, von Hand: Zeitverteilung über quad (dt = dr/v_r), Quasi-Newton g_i ≪ g_e: σ² = σ_N²·ν̄
nm = "Crater II"; M, Re = SCN[nm]["M"], SCN[nm]["Re"]; GMh = GMu(6e10)
for rp, ra in ORB[nm]:
    v = 200.0
    L = math.sqrt(2 * v**2 * math.log(ra / rp) / (1 / rp**2 - 1 / ra**2)); E = v**2 * math.log(rp) + L**2 / (2 * rp**2)
    vr = lambda x: math.sqrt(max(2 * (E - v**2 * math.log(x)) - L**2 / x**2, 1e-30))
    # Substitution x = rp + (ra − rp)(1 − cos u)/2 hebt die Wurzelsingularitäten
    jac = lambda u: (ra - rp) / 2 * math.sin(u)
    xr = lambda u: rp + (ra - rp) * (1 - math.cos(u)) / 2
    Th = quad(lambda u: jac(u) / vr(xr(u)), 0, math.pi, limit=400)[0]
    qbar = quad(lambda u: jac(u) / vr(xr(u)) * math.exp(-math.sqrt(GMh / xr(u)**2)), 0, math.pi, limit=400)[0] / Th
    dq = lambda x: -math.exp(-math.sqrt(x)) / (2 * math.sqrt(x))
    dqbar = quad(lambda u: jac(u) / vr(xr(u)) * dq(GMh / xr(u)**2), 0, math.pi, limit=400)[0] / Th
    nbar = 1 / (1 - qbar); ge0 = SCN[nm]["ge_now"]
    fA = nbar + nbar**2 * ge0 * dqbar / 3
    sN = SCN[nm]["sN"]
    key = f"log, {rp:g}–{ra:g} kpc, aus"
    s_num = RES[(nm, key)]["sig"][0]
    P(f"   b) Crater II, log, {rp}–{ra} kpc, λ → 0 von Hand (g_i ≪ g_e): T_r = {2*Th*TU/Gyr:.3f} Gyr (Integrator {RES[(nm, key)]['T']*TU/Gyr:.3f}), "
      f"⟨q⟩ = {qbar:.5f}, ν̄ = {nbar:.3f}, ⟨q′⟩ = {dqbar:.3f}")
    P(f"      Lesart A: σ ≈ σ_N√(ν̄ + ν̄²g_e0⟨q′⟩/3) = {sN*math.sqrt(fA):.3f} km/s (Rechnung λ/Ω = 10⁻⁴: {s_num:.3f});"
      f"  Lesart B: σ ≈ σ_N√ν̄ = {sN*math.sqrt(nbar):.3f} km/s (Tabelle Lesart B nur für 20–135 kpc)")
ge0 = SCN[nm]["ge_now"]; s0 = math.sqrt(ge0); nu0 = 1 / (-math.expm1(-s0)); K0 = -(nu0 - 1) * s0 / 2
P(f"   c) Crater II, λ = ∞ von Hand: g_e = {ge0:.5f} a₀, ν = {nu0:.3f}, K = {K0:.4f}, σ ≈ σ_N√(ν(1 + K/3)) = "
  f"{SCN[nm]['sN']*math.sqrt(nu0*(1+K0/3)):.3f} km/s (Rechnung {SCN[nm]['s_inf']:.3f})")
# d) Kreisbahn Lesart A
kc = [k for k in RES if k[0] == "DF2" and "kreis" in k[1]][0]
P(f"   d) Kreisbahn (DF2, e = 0), Lesart A: σ für alle λ = " + ", ".join(f"{s:.4f}" for s in RES[kc]["sig"]))
# e) Gewichte: Summe und Mittelwert τ
d = RES[("Fornax", "log, 60–155 kpc, aus")]
for il, l in enumerate(LAMS[:-1]):
    w = d["W"][il]; lam = l * 2 * math.pi / d["T"]; Tt = d["T"]
    exact = 1 / lam - Tt * math.exp(-lam * Tt) / (1 - math.exp(-lam * Tt))
    if il in (0, 3, 4):
        P(f"   e) Fornax, λ/Ω = {l:g}: ⟨τ⟩ der Gewichte = {np.sum(w*d['tau']):.4f}, exakt 1/λ − T e^{{−λT}}/(1 − e^{{−λT}}) = {exact:.4f} (Einheit kpc/(km/s))")

open(__file__.replace(".py", ".txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
