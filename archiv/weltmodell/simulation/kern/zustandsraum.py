# -*- coding: utf-8 -*-
"""
Zustandsraum (Ebene Z): Form des Ganzen als zeitfeste Verhältnisse. Nome und Im τ des flachen Torus,
Re τ_M der Punktmasse und ihre Resonanzradien, Messfehler σ(Q), Bose-Leiter, Rücklaufzahl (Fisher-Metrik, Entropie),
Rückkopplungssumme und Kippstelle, Einrasten, Kac-Budget, Preisach-Gedächtnis, RC-Kette als Nachbau der Erneuerung,
Eisenstein-Gittersumme; Labor: eingespannter Kupferdraht im Magnetfeld (06).
q ist ein Objekt: der Nome e^{2πiτ}; sein Betrag ist der Rückkehranteil, im Regelkreis die Schleifenverstärkung.
Am Ende die Rollen der Größen als Daten: Triade (01) und Rollen im Zustandsraum (05) je Ebene L/I/Z/G.
"""
from .konstanten import *
from .formeln import formel
from .rueckkehr import chi, n_B


# ------------------------------------------------------------------ Nome und Formverhältnisse
@formel('realteil_quadrat_punkt', '(Re τ_M)² = 1 + χ(s), Punktmasse',
        ('Re τ_M',), ('s',), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'γ_N = d ln a_N/d ln r;     Punktmasse: 1 + χ'))
def realteil_quadrat_punkt(s):
    """(Re τ_M)² = 1 + χ(s) für die Punktmasse, s = r*/r (a_N = s²a₀)."""
    return 1 + chi(s * s * A0, A0)


@formel('radius_resonanz', '1 + χ(s) = (p/k)², r_p/k = r*/s',
        ('r_p/k',), ('p/k', 'r*'), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'örtliche Resonanzen bei rationalem κ/Ω'))
def radius_resonanz(p_zu_q, rstern):
    """Radius, an dem Re τ_M = p/q (Punktmasse): 1 + χ(s) = (p/q)²."""
    from scipy.optimize import brentq
    target = p_zu_q**2 - 1
    s = brentq(lambda x: chi(x * x * A0, A0) - target, 1e-6, 50)
    return rstern / s


@formel('Q_fehler', 'σ_Q = Q·√((σ_c/ω_c)² + (σ_b/ω_b)²)',
        ('σ_Q',), ('ω_b', 'ω_c'), 'G', 'hilfe', 'Baustein 03', None)
def Q_fehler():
    """σ(Q) aus σ(ω_b), σ(ω_c), unkorreliert: Q·√((σ_c/ω_c)² + (σ_b/ω_b)²)."""
    sb, sc = OMEGA_BC_H2_FEHLER
    return Q_KOSM * math.hypot(sc / OMEGA_C_H2, sb / OMEGA_B_H2)


@formel('im_tau_aus_nome', 'Im τ_Z = −ln|q|/2π, q = e^{2πiτ_Z}',
        ('Im τ_Z',), ('q',), 'Z', 'modell', 'Baustein 05', ('05-zustandsraum.md', 'Die Spitze q = e^{2πiτ} → 0'))
def im_tau_aus_nome(q):
    """Flacher Torus ℂ/(ℤ+τℤ), Nome q = e^{2πiτ}: Im τ = −ln|q|/2π."""
    return -math.log(q) / (2 * math.pi)


# ------------------------------------------------------------------ Kac-Zählung (06)
@formel('kac_rueckkehrzeit', '⟨T⟩ = μ(X)/μ(A)',
        ('⟨T⟩',), ('μ(X)', 'μ(A)'), 'Z', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'mittlere Rückkehrzeit ⟨T_A⟩ = μ(X)/μ(A)'))
def kac_rueckkehrzeit(mu_X, mu_A):
    """Kac: mittlere Rückkehrzeit (in Schritten) in einen Bereich A vom Maß μ(A) auf dem geschlossenen Raum X."""
    return mu_X / mu_A


@formel('kac_budget', 'Σ_i 1/⟨T_i⟩ = Σ_i μ(A_i)/μ(X) = 1 für jede Zerlegung',
        ('Σ1/⟨T_i⟩',), ('μ(X)', 'μ(A)'), 'Z', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'Σ_i 1/⟨T_{A_i}⟩ = 1'))
def kac_budget(mu_X, mu_A_liste):
    """Kac-Budget einer Zerlegung von X in Bereiche A_i: Σ 1/⟨T_i⟩; 1 genau dann, wenn die Bereiche X ausschöpfen."""
    return sum(1 / kac_rueckkehrzeit(mu_X, m) for m in mu_A_liste)


# ------------------------------------------------------------------ Labor: eingespannter Kupferdraht im Magnetfeld (06)
# Stoffwerte Kupfer aus dem Messwert-Register, Draht, Vorspannung, B und Q_m aus dem Eingabe-Register
from .messwerte import MESSWERTE as _M
from .eingabe import EINGABEN as _E
CU = dict(rho_el=_M["CU_RHO_EL"].wert, alR=_M["CU_AL_R"].wert, E=_M["CU_E"].wert, alE=_M["CU_AL_E"].wert,
          alth=_M["CU_AL_TH"].wert, rho=_M["CU_RHO"].wert, cp=_M["CU_CP"].wert, kth=_M["CU_K_TH"].wert)
DRAHT = dict(a=_E["DRAHT_A"].wert, ell=_E["DRAHT_ELL"].wert, B=_E["DRAHT_B"].wert, Qm=_E["DRAHT_QM"].wert,
             sig0=_E["DRAHT_SIG0"].wert)


@formel('draht', 'R_el,0, G_th, C_th, m_eff, Γ_B = Bℓ_D, k_st(ΔT), ω₀, γ_m (Galerkin-Mode 1 − cos 2πx/ℓ_D)',
        ('R_el,0', 'G_th', 'C_th', 'm_eff', 'Γ_B', 'k_st', 'ω₀', 'γ_m'), ('ρ_el', 'α_R', 'E_Y', 'α_E', 'α_th', 'ρ_Cu', 'c_p', 'k_th', 'a_D', 'ℓ_D', 'B', 'Q_m', 'σ_0'), 'L', 'hilfe', 'Baustein 06', None)
def draht(cu=CU, d=DRAHT):
    """Größen des Drahts (Galerkin-Mode 1 − cos 2πx/ℓ): R0, G_th, C_th, m_eff, Γ = Bℓ, c_m, k(ΔT)."""
    A = math.pi * d['a']**2
    I = math.pi * d['a']**4 / 4
    ell = d['ell']
    p = dict(R0=cu['rho_el'] * ell / A, G_th=12 * cu['kth'] * A / ell, C_th=cu['rho'] * cu['cp'] * A * ell,
             m=1.5 * cu['rho'] * A * ell, Gam=d['B'] * ell, alR=cu['alR'])

    def k(dT):
        Et = cu['E'] * (1 + cu['alE'] * dT)
        N = d['sig0'] * A - cu['E'] * A * cu['alth'] * dT
        return 8 * math.pi**4 * Et * I / ell**3 + 2 * math.pi**2 * N / ell
    p['k'] = k
    p['w0'] = math.sqrt(k(0.0) / p['m'])
    p['c_m'] = p['m'] * p['w0'] / d['Qm']
    return p


# ------------------------------------------------------------------ Bose-Leiter und Rücklaufzahl (07)
@formel('bose_leiter', 'ℬ_d = s^d·n_B(s)',
        ('ℬ_d',), ('s', 'd'), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'dieselbe Funktion s^d·n_B(s)'))
def bose_leiter(s, d):
    """Stufe d der Leiter: d = 0 ν − 1, d = 1 χ, d = 2 (a − a_N)/a₀, d = 3 ¼·dQ̃/ds; jede Stufe ist die vorige mal s."""
    return s**d * n_B(s)


@formel('bose_leiter_integral', 'ℐ_d = ∫₀^∞ s^d·n_B ds = Γ(d+1)ζ(d+1)',
        ('ℐ_d',), ('d',), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', '∫₀^∞ ds = Γ(d+1)ζ(d+1)'))
def bose_leiter_integral(d):
    """Fläche der Stufe d im Maß der Bose-Variablen: π²/6, 2ζ(3), π⁴/15 für d = 1, 2, 3 (d = 0 divergent)."""
    from scipy.special import gamma, zeta
    if d <= 0:
        return math.inf
    return float(gamma(d + 1) * zeta(d + 1))


@formel('fisher_metrik', 'g_F = −dν/ds = ν(ν − 1)',
        ('g_F',), ('s',), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'Fisher-Metrik g = −dν/ds = ν(ν − 1)'))
def fisher_metrik(s):
    """Fisher-Metrik der geometrischen Verteilung der Rücklaufzahl, ν = 1 + n_B(s)."""
    v = 1 + n_B(s)
    return v * (v - 1)


@formel('entropie_ruecklauf', 'S_n = ln ν + χ',
        ('S_n',), ('s',), 'Z', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'ihre Entropie ist H = ln ν + χ'))
def entropie_ruecklauf(s):
    """Entropie der geometrischen Verteilung der Rücklaufzahl (größte Entropie bei festem Mittel ν − 1)."""
    return math.log(1 + n_B(s)) + chi(s * s * A0, A0)


# ------------------------------------------------------------------ Rückkopplungssumme, Kippstelle, Einrasten (05, 06)
@formel('schleifensumme', 'Σ_{n≥0} qⁿ = 1/(1 − q)',
        ('ν',), ('q',), 'Z', 'modell', 'Baustein 05',
        ('05-zustandsraum.md', 'alle Rückläufe summieren sich zu Σ qⁿ = 1/(1 − q)'), alternative='nu')
def schleifensumme(q):
    """Summe aller Rückläufe einer Schleife mit Schleifenverstärkung q (Betrag des Nomes); divergiert bei q = 1."""
    if q >= 1:
        return math.inf
    return 1 / (1 - q)


@formel('schleife_selbsterregt', 'q = (∂P/∂T)/G_th = −g_E·2/(3b_I)·(dR/dT)/G_th, P = g_E⟨I²⟩, ⟨I²⟩ = 2(g_E − R)/(3b_I)',
        ('q',), ('g_E', 'b_I', 'dR/dT', 'G_th'), 'Z', 'modell', 'Baustein 05',
        ('05-zustandsraum.md', 'stabilisiert selbsterregt (q = −0,05)'), alternative='schleife_thermisch')
def schleife_selbsterregt(g_E, b_I, dR_dT, G_th):
    """Selbsterregter Umlauf L·İ = … + (g_E − R)I − b_I·I³: Die Amplitude stellt sich so ein, dass die Entdämpfung g_E
    die ganze Leistung liefert, P = g_E·⟨I²⟩ mit ⟨I²⟩ = 2(g_E − R)/(3b_I) (Mittelung über den Umlauf). Steigt R mit T,
    sinkt P: Dieselbe Kette, die stromgetrieben durchgeht, stabilisiert selbsterregt."""
    return -g_E * 2 / (3 * b_I) * dR_dT / G_th


@formel('kipp_schleife', 'q = 1 − 2(3a)^{1/4}δ^{1/2} vor der Falte',
        ('q',), ('a_F', 'δ_F'), 'Z', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'vor der Falte gilt 1 − q = 2(3a)^{1/4}δ^{1/2}'), alternative='q_')
def kipp_schleife(a, delta):
    """Schleifenverstärkung im Abstand δ vor der Falte der Normalform; bei δ = 0 ist q = 1 (Kippstelle)."""
    return 1 - 2 * (3 * a)**0.25 * math.sqrt(delta)


@formel('rastbreite', '|Δω| ≤ 2k|cos φ|',
        ('Δω_R',), ('k_K', 'φ_K'), 'Z', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'Die Zyklen rasten genau dann ein, wenn |Δω| ≤ 2k|cos φ|'))
def rastbreite(k, phi):
    """Größter Frequenzabstand zweier Phasen (ψ̇ = Δω − 2k cos φ sin ψ), bei dem sie einrasten; nur der dissipative
    Anteil k cos φ rastet ein, bei φ = π/2 nie."""
    return 2 * k * abs(math.cos(phi))


# ------------------------------------------------------------------ Preisach-Gedächtnis (06)
@formel('preisach_ausgang', 'y_P = Σ μ_P·Zustand(α_P, β_P) nach der Folge der nicht ausgelöschten Extremwerte',
        ('y_P',), ('u_P', 'α_P', 'β_P', 'μ_P'), 'Z', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'Ein Preisach-Zustand ist nur die Folge der nicht ausgelöschten Extremwerte'))
def preisach_ausgang(folge, alpha, beta, mass):
    """Hysteronen (α ≥ β) aus dem Zustand −1 (Geschichte von −∞): auf +1 bei u ≥ α, auf −1 bei u ≤ β.
    Gerechnet wird nur mit der Folge der nicht ausgelöschten Extremwerte (ein Extremum mit größerer Spanne löscht das
    letzte Paar); wenige Zahlen genügen für eine lange Geschichte."""
    st = []
    for u in folge:
        u = float(u)
        if st and u == st[-1]:
            continue
        if len(st) >= 2 and (st[-1] - st[-2]) * (u - st[-1]) > 0:
            st[-1] = u
        else:
            st.append(u)
        while len(st) >= 3 and abs(st[-1] - st[-2]) >= abs(st[-2] - st[-3]):
            del st[-3:-1]
    alpha, beta, mass = np.asarray(alpha, float), np.asarray(beta, float), np.asarray(mass, float)
    z = -np.ones(alpha.shape)
    for u in st:
        z = np.where(u >= alpha, 1.0, np.where(u <= beta, -1.0, z))
    return float(mass @ z)


# ------------------------------------------------------------------ Labor der Rückkehr: RC-Kette (06)
@formel('rc_kette_nu', 'ν = 1/(1 − e^{−s}), s = n_R·arcosh(1 + RG/2)',
        ('ν',), ('n_R', 'RG'), 'Z', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 's = n_R·arcosh(1 + RG/2)'),
        alternative='nu')
def rc_kette_nu(n_R, RG):
    """RC-Kette (Kabelgleichung) mit Rückführung des Abgriffs n_R auf den Eingang, V(0) = V_ein + V(n_R): Nachbau der
    Erneuerung u(0) = g_N + u(ϖ_R). Verstärkung V(0)/V_ein = 1 + n_B(s)."""
    return 1 + float(n_B(n_R * math.acosh(1 + RG / 2)))


# ------------------------------------------------------------------ Modulformen: Gittersumme und markierter Torus (05)
@formel('eisenstein_G', 'G_k(τ) = Σ′(m + nτ)^{−k} = 2ζ(k) + 2(2πi)^k/(k−1)!·Σ σ_{k−1}(n)qⁿ',
        ('G_k',), ('k_w', 'τ'), 'Z', 'modell', 'Baustein 05',
        ('05-zustandsraum.md', 'die Eisenstein-Reihe G_k(τ) = Σ′(m + nτ)^{−k}'))
def eisenstein_G(k, tau):
    """𝒫 auf dem Torus ℂ/(ℤ + τℤ), angewandt auf z^{−k} (k gerade ≥ 4), längs des Umlaufs zuerst summiert: jeder
    Rückkehrschritt n gibt eine Reihe im Nome q = e^{2πiτ}. Reell für τ auf der imaginären Achse."""
    from scipy.special import zeta
    if k != int(k) or k < 4 or int(k) % 2:
        raise ValueError("Gewicht k gerade und ≥ 4")
    k = int(k)
    tau = complex(tau)
    if tau.imag < 0.05:
        raise ValueError("Im τ ≥ 0,05 (sonst Reihe zu lang)")
    q = np.exp(2j * math.pi * tau)
    summe = 0j
    for n in range(1, int(40 / (2 * math.pi * tau.imag)) + 2):
        summe += sum(d**(k - 1) for d in range(1, n + 1) if n % d == 0) * q**n
    w = 2 * float(zeta(k)) + 2 * (-1)**(k // 2) * (2 * math.pi)**k / math.factorial(k - 1) * summe
    return w.real if tau.real == 0 else w


@formel('nu_aus_modul', 'ν = 1/(1 − |e^{2πiτ}|)',
        ('ν',), ('τ',), 'Z', 'modell', 'Baustein 05', ('05-zustandsraum.md', 'τ = i/2 und 2i geben 1,045 und 1,0000035'),
        alternative='nu')
def nu_aus_modul(tau):
    """Verstärkung am markierten Torus: nur der Betrag des Nomes bindet; die Modulgruppe vertauscht Umlauf und
    Rückkehr, ν aber nicht."""
    return 1 / (1 - abs(np.exp(2j * math.pi * complex(tau))))


@formel('delta_reihe', 'Δ = q∏_{n≥1}(1 − qⁿ)²⁴ = Σ c_Δ(n)qⁿ',
        ('c_Δ',), (), 'Z', 'modell', 'Baustein 05', ('05-zustandsraum.md', 'q∏(1 − qⁿ)²⁴ = q − 24q² + 252q³'))
def delta_reihe(n_max):
    """Die ersten n_max Koeffizienten c_Δ(1 … n_max) der Diskriminante als Produkt ausmultipliziert (ganzzahlig)."""
    n_max = int(n_max)
    if not 1 <= n_max <= 200:
        raise ValueError("1 ≤ n_max ≤ 200")
    c = [0] * (n_max + 1)
    c[1] = 1
    for n in range(1, n_max):
        for _ in range(24):
            for j in range(n_max, n - 1, -1):
                c[j] -= c[j - n]
    return np.array(c[1:], dtype=float)


@formel('delta_modul', 'Δ(τ) = q∏(1 − qⁿ)²⁴, q = e^{2πiτ}',
        ('Δ_τ',), ('τ',), 'Z', 'modell', 'Baustein 05', ('07-feldgleichung.md', 'bei r* nur noch 4,7·10⁻⁸ der Reihe'))
def delta_modul(tau):
    """Diskriminante am Punkt τ der oberen Halbebene als Produkt; reell für τ auf der imaginären Achse."""
    tau = complex(tau)
    if tau.imag < 0.02:
        raise ValueError("Im τ ≥ 0,02 (sonst Produkt zu lang)")
    q = np.exp(2j * math.pi * tau)
    n = np.arange(1, int(60 / (2 * math.pi * tau.imag)) + 2)
    w = q * np.exp(24 * np.sum(np.log1p(-q**n)))
    return w.real if tau.real == 0 else w


@formel('resolvente_flucht', 'ν_F = max_{1≤k≤K_M} 1/|1 − q·e^{2πikα_D}|',
        ('ν_F',), ('q', 'K_M', 'α_D'), 'Z', 'modell', 'Baustein 05', ('05-zustandsraum.md', 'Fluchtmoden höchstens 567-fach'))
def resolvente_flucht(q, K_M, alpha=(math.sqrt(5) - 1) / 2):
    """Rückkehr als Translation T (Drehung um α_D): Die Resolvente (I − qT)⁻¹ verstärkt die invariante Mode um
    1/(1 − q), jede Fluchtmode e^{2πikx} nur um 1/|1 − qλ_k|, λ_k = e^{2πikα_D}; gibt das größte davon."""
    if not 0 <= q < 1 or not 1 <= K_M <= 1e6:
        raise ValueError("0 ≤ q < 1, 1 ≤ K_M ≤ 10⁶")
    k = np.arange(1, int(K_M) + 1)
    return float(np.max(1 / np.abs(1 - q * np.exp(2j * math.pi * k * alpha))))


@formel('spiel_umkehrpunkt', 'v_K²/2 − 1 = −1/x_U + γ_K(x_U − 1)²/2 (H = v²/2 − 1/x + γ_K·y²/2, y = x − 1, t_K → ∞)',
        ('x_U',), ('v_K', 'γ_K'), 'Z', 'modell', 'Baustein 05', ('05-zustandsraum.md', 'kehrt bei γ = 0,01 an x = 9,75 um'))
def spiel_umkehrpunkt(v, gamma):
    """Spielbeispiel des Gedächtnisses (Einheiten GM = 1, Start x = 1): Mit unendlichem Gedächtnis (t_K → ∞) wird der
    Maxwell-Spielkern zur Feder γ_K·y um den Startpunkt; die Energie bleibt erhalten, das Teilchen kehrt dort um, wo
    v = 0 wird, auch über der Fluchtgeschwindigkeit √2."""
    from scipy.optimize import brentq
    if not (v > 0 and 0 < gamma < 10):
        raise ValueError("v > 0, 0 < γ < 10")
    E = 0.5 * v * v - 1
    f = lambda x: -1 / x + 0.5 * gamma * (x - 1)**2 - E
    hi = 2.0
    while f(hi) < 0:
        hi *= 2
    return brentq(f, 1.0, hi, xtol=1e-14)


@formel('rate_rastpunkt', 'λ_ψ = −√((2k cos φ)² − Δω²)',
        ('λ_ψ',), ('k_K', 'φ_K', 'Δω'), 'Z', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'ein Exponent wird −0,098'))
def rate_rastpunkt(k, phi, domega):
    """Rate am stabilen Rastpunkt von ψ̇ = Δω − 2k cos φ sin ψ; nur innerhalb der Rastbreite reell."""
    b = 2 * k * math.cos(phi)
    if abs(domega) > abs(b):
        return math.nan
    return -math.sqrt(b * b - domega * domega)


@formel('vo2_hopf', 'Re λ = 0 des elektrothermischen Paars der Jacobi-Matrix am Fixpunkt des VO₂-Elements',
        ('V_H', 'f_H'), ('𝔄_VO₂',), 'L', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'Hopf bei 3,55 V (3,26 MHz)'))
def vo2_hopf(aufbau):
    """VO₂ auf einer Mikrobrücke im Magnetfeld: Quellspannung V_H [V], unterhalb der das elektrothermische Paar
    instabil wird (Geburt des Zyklus), und seine Frequenz f_H [Hz] dort. Fixpunkt über die Kennlinie G(T − T₀) = V²/R
    mit iterierter Auslenkung, Jacobi-Matrix als zentrale Differenzen, Paar mit der größten Frequenz."""
    a = dict(aufbau)

    def _vo2_fluss(Z, Vs):
        """Fluss (V, T, x, v) des VO₂-Elements mit Lastwiderstand, Parallelkapazität und Mikrobrücke im Magnetfeld B
        (Einheiten µs, K, V, kΩ, mA, mW, nF, nJ/K, µm): C·V̇ = (V_s − V)/R_L − I, I = (V − Γ_B·v)/R(T, x),
        C_th·Ṫ = I²R − G_th(T − T₀), ẍ = −(ω/Q)ẋ − ω²y(1 + (y/x_D)²), y = x − D_tr(1 − s(T)) − λ_B·I; die Dehnung
        verschiebt die Übergangstemperatur um κ·x."""
        V, T, x, v = Z
        s = lambda sh: 1 / (1 + np.exp((T - a['T_C'] - sh) / a['breite']))
        R = np.exp(math.log(a['R_M']) + math.log(a['R_I'] / a['R_M']) * s(a['kappa'] * x))
        I = (V - a['Gamma'] * a['B'] * v) / R
        om = 2 * math.pi * a['f_m']
        y = x - (a['D_tr'] * (1 - s(0.0)) + a['lambda_B'] * a['B'] * I)
        return np.array([((Vs - V) / a['R_L'] - I) / a['C'], (I * I * R - a['G_th'] * (T - a['T_0'])) / a['C_th'], v,
                         -om / a['Q'] * v - om * om * y * (1 + (y / a['x_D'])**2)])

    def fixpunkt(Vs):
        Ts = np.linspace(a['T_0'] + 0.01, a['T_0'] + 70.0, 70001)
        x = 0.0
        for _ in range(100):
            R = np.exp(math.log(a['R_M']) + math.log(a['R_I'] / a['R_M'])
                       / (1 + np.exp((Ts - a['T_C'] - a['kappa'] * x) / a['breite'])))
            V = Vs * R / (R + a['R_L'])
            g = V * V / R - a['G_th'] * (Ts - a['T_0'])
            i = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
            if len(i) != 1:
                raise ValueError("kein eindeutiger Fixpunkt")
            k = i[0]
            T = Ts[k] - g[k] * (Ts[k + 1] - Ts[k]) / (g[k + 1] - g[k])
            Rt = math.exp(math.log(a['R_M']) + math.log(a['R_I'] / a['R_M'])
                          / (1 + math.exp((T - a['T_C'] - a['kappa'] * x) / a['breite'])))
            V = Vs * Rt / (Rt + a['R_L'])
            xn = a['D_tr'] / (1 + math.exp(-(T - a['T_C']) / a['breite'])) + a['lambda_B'] * a['B'] * V / Rt
            if abs(xn - x) < 1e-13:
                break
            x = xn
        return np.array([V, T, x, 0.0])

    def paar(Vs):
        z = fixpunkt(Vs)
        J = np.zeros((4, 4))
        for i in range(4):
            dz = np.zeros(4)
            dz[i] = 1e-6 * max(1.0, abs(z[i]))
            J[:, i] = (_vo2_fluss(z + dz, Vs) - _vo2_fluss(z - dz, Vs)) / (2 * dz[i])
        ev = np.linalg.eigvals(J)
        ev = ev[ev.imag > 0]
        return ev[int(np.argmax(ev.imag))]
    lo, hi = a['V_suche']
    for _ in range(40):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if paar(m).real > 0 else (lo, m)
    VH = 0.5 * (lo + hi)
    return VH, paar(VH).imag / (2 * math.pi) * 1e6


_KURAMOTO = {}


@formel('kuramoto_einrasten', 'θ̇_i = ω_i + k_K·Σ_j sin(θ_j − θ_i), ω_i gleichabständig im Band W; '
        '(N − 1)k₁/W beim ersten Einrasten',
        ('(N−1)k₁/W',), ('N_K',), 'Z', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'erstes Einrasten bei (N − 1)k ≈ 0,57–0,63 W'))
def kuramoto_einrasten(N_liste):
    """Je N: kleinste Kopplung k₁ auf dem Gitter k·N = 0,30 … 1,20 (Schritt 0,01, Band W = 1), bei der zwei Zyklen
    dieselbe mittlere Frequenz haben (Abstand ≤ 2·10⁻³ über die zweite Hälfte von t = 3000; DOP853, rtol 10⁻⁹; Start
    gleichverteilt auf dem Kreis); gibt (N − 1)k₁/W je N. Alle N und alle k laufen in einem Zustand."""
    from scipy.integrate import solve_ivp
    roh = np.atleast_1d(np.asarray(N_liste, dtype=float))
    if roh.size == 0 or np.any(roh != np.round(roh)) or np.any((roh < 2) | (roh > 64)):
        raise ValueError("N ganzzahlig, 2 ≤ N ≤ 64")
    N_liste = tuple(int(n) for n in roh)
    if N_liste in _KURAMOTO:
        return _KURAMOTO[N_liste]
    kN = np.round(np.arange(0.30, 1.2001, 0.01), 10)
    teile = [(n, np.linspace(-0.5, 0.5, n), (kN / n)[:, None]) for n in N_liste]
    grenzen = np.concatenate([[0], np.cumsum([kN.size * n for n in N_liste])])

    def f(t, y):
        aus = np.empty_like(y)
        for j, (n, om, kk) in enumerate(teile):
            th = y[grenzen[j]:grenzen[j + 1]].reshape(kN.size, n)
            cs, sn = np.cos(th), np.sin(th)
            aus[grenzen[j]:grenzen[j + 1]] = (om[None, :] + kk * (sn.sum(1, keepdims=True) * cs
                                                                   - cs.sum(1, keepdims=True) * sn)).ravel()
        return aus
    y0 = np.concatenate([np.tile(np.linspace(0, 2 * math.pi, n, endpoint=False), kN.size) for n in N_liste])
    T = 3000.0
    s = solve_ivp(f, (0, T), y0, method='DOP853', rtol=1e-9, atol=1e-9, t_eval=[T / 2, T])
    Om = (s.y[:, 1] - s.y[:, 0]) / (T / 2)
    aus = []
    for j, n in enumerate(N_liste):
        teil = Om[grenzen[j]:grenzen[j + 1]].reshape(kN.size, n)
        gruppen = 1 + np.sum(np.diff(np.sort(teil, axis=1), axis=1) > 2e-3, axis=1)
        rast = np.nonzero(gruppen < n)[0]
        aus.append((n - 1) * kN[rast[0]] / n if rast.size else math.nan)
    _KURAMOTO[N_liste] = np.array(aus)
    return _KURAMOTO[N_liste]


@formel('rate_quer', 'λ_⊥ = −D_Z',
        ('λ_⊥',), ('D_Z',), 'Z', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'λ_⊥ = −0,5'))
def rate_quer(D):
    """Transversale Rate des Torus T² × S: Die Querrichtung S zieht mit der Dissipation D auf den Torus."""
    return -D


# ------------------------------------------------------------------ Labor: Kupferdraht, Hüllkurve und Hopf-Grenze (06)
@formel('kupfer_huellkurve', 'Ẋ = [Γ_B·I − (k_st(ΔT) − mω² + iωγ_m)X]/(2iωm), I = (V − iωΓ_B X)/R(ΔT), '
        'C_thΔṪ = ½R|I|² + ½γ_mω²|X|² − G_thΔT',
        ('ΔT_D',), ('Δf_d', 'V_D', 'Q_m'), 'L', 'hilfe', 'Baustein 06', None)
def kupfer_huellkurve(dfd, V, Qm):
    """Hüllkurve X(t) (x = Re X e^{iωt}) des Drahts bei Anregung V cos ωt, ω = ω₀ + 2πΔf_d: dict mit den stationären
    Punkten ΔT (aufsteigend), dem heißen Ast (größtes ΔT), seinem Zustand (Re X, Im X, ΔT), den Eigenwerten der
    Hüllkurve dort und dem Fluss (Re Ẋ, Im Ẋ, ΔṪ) als Funktion des Zustands."""
    from scipy.optimize import brentq
    d = dict(DRAHT)
    d['Qm'] = Qm
    p = draht(CU, d)
    R0, alR, Gam, m, cm, k, Gth, Cth = (p[x] for x in ('R0', 'alR', 'Gam', 'm', 'c_m', 'k', 'G_th', 'C_th'))
    wd = p['w0'] + 2 * math.pi * dfd

    def rhs(y):
        X, dT = y[0] + 1j * y[1], y[2]
        Rt = R0 * (1 + alR * dT)
        Ic = (V - 1j * wd * Gam * X) / Rt
        dX = (Gam * Ic - (k(dT) - m * wd**2 + 1j * wd * cm) * X) / (2j * wd * m)
        return np.array([dX.real, dX.imag, (0.5 * Rt * abs(Ic)**2 + 0.5 * cm * wd**2 * abs(X)**2 - Gth * dT) / Cth])

    def X_stat(dT):
        Rt = R0 * (1 + alR * dT)
        return Gam * V / Rt / (k(dT) - m * wd**2 + 1j * wd * cm + 1j * wd * Gam**2 / Rt)

    def f_stat(dT):
        Rt = R0 * (1 + alR * dT)
        X = X_stat(dT)
        Ic = (V - 1j * wd * Gam * X) / Rt
        return 0.5 * Rt * np.abs(Ic)**2 + 0.5 * cm * wd**2 * np.abs(X)**2 - Gth * dT

    xs = np.linspace(0.0, 3.0, 30001)
    fv = f_stat(xs)
    idx = np.nonzero(np.sign(fv[:-1]) != np.sign(fv[1:]))[0]
    st = [brentq(f_stat, xs[i], xs[i + 1], xtol=1e-14) for i in idx]
    aus = dict(stationaer=st, heiss=None, zustand=None, eigenwerte=None, fluss=rhs)
    if st:
        X = X_stat(st[-1])
        y = np.array([X.real, X.imag, st[-1]])
        J = np.zeros((3, 3))
        sk = [abs(y[0]) + abs(y[1]) + 1e-12] * 2 + [1e-3]
        for j in range(3):
            e = np.zeros(3)
            e[j] = 1e-6 * sk[j]
            J[:, j] = (rhs(y + e) - rhs(y - e)) / (2 * e[j])
        aus.update(heiss=st[-1], zustand=y, eigenwerte=np.linalg.eigvals(J))
    return aus


_KUPFER = {}


@formel('kupfer_hopf', 'heißer Ast: Falte = Ende der drei stationären Punkte; Hopf-Grenze = Re λ des Paars = 0 davor',
        ('Δf_F',), ('V_D', 'Q_m'), 'L', 'hilfe', 'Baustein 06', None)
def kupfer_hopf(Qm, V):
    """Längs des heißen Asts zur Falte hin (Verstimmung nach unten): gibt (Abstand der Hopf-Grenze vor der Falte [Hz],
    Falte Δf_d [Hz], Ringfrequenz |Im λ|/2π an der Hopf-Grenze [Hz]); ohne Hopf-Grenze (nan, Falte, nan).
    Hopf-Grenze: der am weitesten von der Falte entfernte Punkt, an dem der führende Eigenwert die Achse als komplexes
    Paar kreuzt (davor stabil, dahinter bis zur Falte instabil)."""
    schluessel = (float(Qm), float(V))
    if schluessel in _KUPFER:
        return _KUPFER[schluessel]

    def drei(d):
        return len(kupfer_huellkurve(d, V, Qm)['stationaer']) >= 3

    def fuehrend(d):
        """größter Realteil der Eigenwerte am heißen Ast und der Imaginärteil dieses Eigenwerts"""
        ev = kupfer_huellkurve(d, V, Qm)['eigenwerte']
        e = max(ev, key=lambda x: x.real)
        return e.real, abs(e.imag)
    gitter = np.arange(0.0, -400.0, -1.0)
    ja = [i for i, d in enumerate(gitter) if drei(d)]
    aus = (math.nan, math.nan, math.nan)
    if ja:
        a, b = gitter[ja[-1]], gitter[ja[-1]] - 1.0
        for _ in range(34):
            mitte = 0.5 * (a + b)
            a, b = (mitte, b) if drei(mitte) else (a, mitte)
        falte = a
        abst = np.geomspace(1e-3, 5.0, 40)
        re = [fuehrend(falte + x)[0] for x in abst]
        aus = (math.nan, falte, math.nan)
        instabil = [j for j in range(len(abst)) if re[j] > 0]
        if instabil and instabil[-1] < len(abst) - 1:
            lo, hi = abst[instabil[-1]], abst[instabil[-1] + 1]
            for _ in range(30):
                mitte = 0.5 * (lo + hi)
                lo, hi = (mitte, hi) if fuehrend(falte + mitte)[0] > 0 else (lo, mitte)
            im = fuehrend(falte + lo)[1]
            if im > 0:
                aus = (lo, falte, im / (2 * math.pi))
    _KUPFER[schluessel] = aus
    return aus


@formel('kupfer_hopf_abstand', 'Δf_H = Abstand der Hopf-Grenze vor der Falte des heißen Asts',
        ('Δf_H',), ('V_D', 'Q_m'), 'L', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'endet der heiße Ast 0,0010–0,6 Hz vor seiner Falte in einer unterkritischen Hopf-Verzweigung'))
def kupfer_hopf_abstand(Qm, V):
    return kupfer_hopf(Qm, V)[0]


@formel('kupfer_ringfrequenz', 'f_ring = |Im λ|/2π an der Hopf-Grenze',
        ('f_ring',), ('V_D', 'Q_m'), 'L', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'Ringfrequenz 0,20–1,5 Hz'))
def kupfer_ringfrequenz(Qm, V):
    return kupfer_hopf(Qm, V)[2]


@formel('kupfer_zeitverhaeltnis', 'τ_th/τ_m = (C_th/G_th)/(2Q_l/ω₀), Q_l = Q_m/(1 + Γ_B²/(γ_m R_el,0))',
        ('τ_th/τ_m',), ('Q_m',), 'L', 'modell', 'Baustein 06', ('06-torus-dynamik.md', 'Wärmezeit/Ringzeit ≤ 0,74'))
def kupfer_zeitverhaeltnis(Qm):
    """Wärmezeit gegen Ringzeit der Amplitude (belastete Güte Q_l); die Hopf-Grenze verlangt eine Ringzeit, die länger
    ist als die Wärmezeit."""
    d = dict(DRAHT)
    d['Qm'] = Qm
    p = draht(CU, d)
    Ql = p['m'] * p['w0'] / p['c_m'] / (1 + p['Gam']**2 / p['c_m'] / p['R0'])
    return (p['C_th'] / p['G_th']) / (2 * Ql / p['w0'])


@formel('kupfer_lyapunov', 'ℓ₁ = Re[⟨p, C(q,q,q̄)⟩ − 2⟨p, B(q, J⁻¹B(q,q̄))⟩ + ⟨p, B(q̄, (2iω − J)⁻¹B(q,q))⟩]/2ω '
        'an der Hopf-Grenze der Hüllkurve (Re X, Im X in |X₀|, ΔT in K)',
        ('ℓ₁',), ('V_D', 'Q_m'), 'L', 'modell', 'Baustein 06',
        ('06-torus-dynamik.md', 'in einer unterkritischen Hopf-Verzweigung'))
def kupfer_lyapunov(Qm, V):
    """Erster Lyapunov-Koeffizient der Hüllkurve an der Hopf-Grenze des heißen Asts. ℓ₁ > 0: unterkritisch (ein
    instabiler Zyklus umgibt den stabilen Punkt davor, dahinter springt der Draht vom heißen Ast ab, kein
    dauerhaftes Selbstpulsieren); ℓ₁ < 0: überkritisch (kleiner stabiler Zyklus). Ableitungen als zentrale
    Differenzen in den Koordinaten (Re X/|X₀|, Im X/|X₀|, ΔT/K). Ohne Hopf-Grenze nan."""
    abst, falte, _ = kupfer_hopf(Qm, V)
    if not math.isfinite(abst):
        return math.nan
    a = kupfer_huellkurve(falte + abst, V, Qm)
    y0, fl = a['zustand'], a['fluss']
    sk = np.array([math.hypot(y0[0], y0[1])] * 2 + [1.0])
    F = lambda u: fl(u * sk) / sk
    u0 = y0 / sk
    e = np.eye(3)

    def d1(j, h=1e-7):
        return (F(u0 + h * e[j]) - F(u0 - h * e[j])) / (2 * h)

    def d2(j, k, h=1e-4):
        return (F(u0 + h * (e[j] + e[k])) - F(u0 + h * (e[j] - e[k])) - F(u0 - h * (e[j] - e[k]))
                + F(u0 - h * (e[j] + e[k]))) / (4 * h * h)

    def d3(j, k, l, h=2e-3):
        s = np.zeros(3)
        for sj in (1, -1):
            for sk_ in (1, -1):
                for sl in (1, -1):
                    s = s + sj * sk_ * sl * F(u0 + h * (sj * e[j] + sk_ * e[k] + sl * e[l]))
        return s / (8 * h**3)
    J = np.array([d1(j) for j in range(3)]).T
    B2 = np.zeros((3, 3, 3))
    C3 = np.zeros((3, 3, 3, 3))
    for j in range(3):
        for k in range(j, 3):
            B2[:, j, k] = B2[:, k, j] = d2(j, k)
            for l in range(k, 3):
                w = d3(j, k, l)
                for (x, y_, z) in {(j, k, l), (j, l, k), (k, j, l), (k, l, j), (l, j, k), (l, k, j)}:
                    C3[:, x, y_, z] = w
    ev, vr = np.linalg.eig(J)
    i = int(np.argmax(ev.imag))
    om, q = ev[i].imag, vr[:, i]
    evl, vl = np.linalg.eig(J.T)
    pv = vl[:, int(np.argmin(np.abs(evl - np.conj(ev[i]))))]
    pv = pv / np.conj(np.vdot(pv, q))
    B = lambda x, y_: np.einsum('ijk,j,k->i', B2, x, y_)
    C = lambda x, y_, z: np.einsum('ijkl,j,k,l->i', C3, x, y_, z)
    h11 = np.linalg.solve(J, B(q, np.conj(q)))
    h20 = np.linalg.solve(2j * om * np.eye(3) - J, B(q, q))
    return float((np.vdot(pv, C(q, q, np.conj(q))) - 2 * np.vdot(pv, B(q, h11))
                  + np.vdot(pv, B(np.conj(q), h20))).real / (2 * om))


@formel('draht_schleife', 'q = (∂P/∂T)/G_th am stationären Punkt, zerlegt: q_σ (nur R(T)), q_E (nur E(T)), q_N (nur '
        'Wärmespannung N(T)), die übrigen Glieder am Arbeitspunkt festgehalten',
        ('q', 'q_σ', 'q_E', 'q_N'), ('Δf_d', 'V_D', 'Q_m', 'ΔT_D'), 'L', 'modell', 'Baustein 09',
        ('09-elektromagnetismus-und-frequenzleiter.md', 'über die Wärmespannung in der Struktur 200- bis 1900-mal stärker'),
        alternative='schleife_thermisch')
def draht_schleife(dT, dfd, V, Qm):
    """Schleife Feld → Strom → Wärme → Regel → Feld am eingespannten Draht (quasistatisch, Impedanz der Mode
    Z_m = Γ_B²/(γ_m + i(mω − k_st/ω))): P = ½|V/(R + Z_m)|²·(R + Re Z_m). Gibt (q, q_σ, q_E, q_N); je Kanal ist nur
    diese Regel temperaturabhängig, die übrigen stehen am Arbeitspunkt ΔT."""
    d = dict(DRAHT)
    d['Qm'] = Qm
    p = draht(CU, d)
    A = math.pi * d['a']**2
    I = math.pi * d['a']**4 / 4
    ell = d['ell']
    wd = p['w0'] + 2 * math.pi * dfd

    def P(tR, tE, tN):
        R = p['R0'] * (1 + p['alR'] * tR)
        k = (8 * math.pi**4 * CU['E'] * (1 + CU['alE'] * tE) * I / ell**3
             + 2 * math.pi**2 * (d['sig0'] * A - CU['E'] * A * CU['alth'] * tN) / ell)
        Zm = p['Gam']**2 / (p['c_m'] + 1j * (p['m'] * wd - k / wd))
        return 0.5 * abs(V / (R + Zm))**2 * (R + Zm.real)
    from .teilchen_stern import schleife_thermisch
    h = 1e-6
    ab = lambda f: schleife_thermisch((f(dT + h) - f(dT - h)) / (2 * h), p['G_th'])
    return (ab(lambda x: P(x, x, x)), ab(lambda x: P(x, dT, dT)), ab(lambda x: P(dT, x, dT)),
            ab(lambda x: P(dT, dT, x)))


# ------------------------------------------------------------------ Rollen der Größen als Daten (01, 05, 00)
# Triade (01, Hodge): Hintergrund = harmonischer Teil (auf einem Teilgebiet das Feld der Elternebene), Wirbel = übrige
# Neigung und Rotation (aufgelöste Struktur), Zentrum = Fluss GM des Gradienten (Faltung). Rollen im Zustandsraum (05):
# tangential (verlustfrei, Umläufe und Phasen), transversal (dissipativ, zieht auf den Torus), parametrisch (bewegt den
# Torus). Ebene je Größe steht im Symbol-Register. Jede Zuordnung nennt ihre Stelle (Datei, Stichwort wörtlich).
TRIADE_GLIEDER = {
    "Hintergrund": ("harmonisch: das gleichförmige Feld, auf einem Teilgebiet das Feld der Elternebene",
                    ("01-grundlagen.md", "harmonisch (Hintergrund, das gleichförmige Feld)")),
    "Wirbel": ("übrige Neigung und Rotation, die aufgelöste Struktur um das Zentrum",
               ("01-grundlagen.md", "der Wirbel die aufgelöste Struktur um ihn")),
    "Zentrum": ("Fluss GM des Gradienten, die Faltung",
                ("01-grundlagen.md", "das Zentrum ist der Fluss GM des Gradienten")),
}

ROLLEN_ZUSTANDSRAUM = {
    "tangential": ("θ̇ = ω(p)", "Umläufe und Phasen, verlustfrei",
                   ("05-zustandsraum.md", "θ̇ = ω(p) | Umläufe und Phasen")),
    "transversal": ("ρ̇ = −κ(p)ρ", "zieht auf den Torus, dissipativ",
                    ("05-zustandsraum.md", "ρ̇ = −κ(p)ρ | zieht auf den Torus, stabilisiert")),
    "parametrisch": ("ṗ = εG", "bewegt den Torus selbst",
                     ("05-zustandsraum.md", "ṗ = εG | bewegt den Torus selbst")),
}

# (Zeichen, Glied der Triade oder None, Stelle, Rolle im Zustandsraum oder None, Stelle)
GROESSEN_ROLLEN = (
    ("Δs", "Zentrum", ("01-grundlagen.md", "Faltung (Zusatzstrecke Δs im Zentrum)"),
     "parametrisch", ("06-torus-dynamik.md", "Gravitation: C, langsame Parametrisierung")),
    ("g", "Wirbel", ("01-grundlagen.md", "lokaler Wirbel (Neigung g)"),
     "parametrisch", ("06-torus-dynamik.md", "Gravitation: C, langsame Parametrisierung")),
    ("dτ/dt", None, None,
     "parametrisch", ("05-zustandsraum.md", "sie verschiebt alle Frequenzen gleich und parametrisiert nur die Zeit um")),
    ("𝒯", "Hintergrund", ("01-grundlagen.md", "Das Feld der Elternebene ist dort Hintergrund der Kindebene"),
     None, None),
    ("P_ϖ", "Wirbel", ("01-grundlagen.md", "nur die eigene aufgelöste Struktur klingt ab"), None, None),
    ("M", "Zentrum", ("01-grundlagen.md", "bleiben der Fluss GM (Zentrum)"), None, None),
    ("Λ", "Zentrum", ("03-von-der-galaxie-bis-zum-kosmos.md", "Λ = ν_kosm·ΣΔs_b"), None, None),
    ("L", "Hintergrund", ("01-grundlagen.md", "ihnen steht die gleichförmige Nullmode gegenüber"), None, None),
    # q ist ein Objekt (Nome): der Betrag |q| trägt die transversale Rolle, die Phase arg q die tangentiale (gleicher Stamm q)
    ("|q|", None, None,
     "transversal", ("00-bauplan.md", "Sein Betrag |q| = e^{−s} ist der Verlust")),
    ("arg q", None, None,
     "tangential", ("00-bauplan.md", "Seine Phase ist die Windung")),
    ("Γ", None, None,
     "transversal", ("05-zustandsraum.md", "sind dissipativ und transversal")),
    ("Q_el", None, None,
     "tangential", ("00-bauplan.md", "Seine Phase ist die Windung (verlustfrei, beide Sinne): Ladung")),
    ("Re τ_M", None, None,
     "tangential", ("07-feldgleichung.md", "Realteil = Umlaufphase")),
    ("Δω_R", None, None,
     "transversal", ("05-zustandsraum.md", "Einrasten braucht Dissipation")),
)

# Grenzlängen jeder Ebene (08): geometrische Leiter mit dem Schritt δ*, symmetrisch um r*
GRENZLAENGEN = (("3Δs", "Faltung", "L"), ("r*", "Schwellenradius", "L"), ("ℓ₀", "Wirkungslänge", "I"))
GRENZLAENGEN_ORT = ("08-licht-und-trichter.md", "Die Leiter der Längen λ̄, ℓ_P, 3Δs, r*, ℓ₀")

# Brücken der Ebenen: (von, nach, Gleichung, Formel-Kennung oder None, Stelle)
BRUECKEN = (
    ("L", "I", "r*² = 3Δs·ℓ₀", "r_stern_leiter", ("08-licht-und-trichter.md", "r*² = 3Δs·ℓ₀")),
    ("L", "I", "r*² = r₀·Λ", None, ("07-feldgleichung.md", "r*² = r₀·Λ")),
    ("L", "I", "Γ = a_N/c", "verlust_gamma", ("07-feldgleichung.md", "Γ = c|∇σ| = a_N/c")),
    ("I", "G", "Λ = ν_kosm·ΣΔs_b", None, ("03-von-der-galaxie-bis-zum-kosmos.md", "Λ = ν_kosm·ΣΔs_b")),
    ("L", "G", "R_Q = r*/(−ln q)", "r_Q", ("03-von-der-galaxie-bis-zum-kosmos.md", "R_Q = r*/(−ln q)")),
    # Ebene Z: Zustandsraum ↔ Ort über die Regel (𝒫, L_eff = 𝒫L); Zustandsraum ↔ Ganzes über das Kac-Budget
    ("Z", "L", "L_eff = 𝒫L", None, ("05-zustandsraum.md", "L_eff = 𝒫L")),
    ("Z", "G", "Σ_i 1/⟨T_{A_i}⟩ = 1", "kac_budget", ("06-torus-dynamik.md", "Σ_i 1/⟨T_{A_i}⟩ = 1")),
)
