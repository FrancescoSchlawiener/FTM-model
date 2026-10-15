# -*- coding: utf-8 -*-
"""
Faltung (Ebene L, lokal am Ort): Neigung g, Faltungstiefe Δs, Takt, δ, r0, Schwellradius r*.
Die Leiter r*² = 3Δs·ℓ₀ verbindet Δs (L) mit ℓ₀ (I). Dazu die Hierarchie der Trichter (Glättung in der Gezeit der
Elternebene) und das Flussgesetz mit der Regel (Flussdimension, Epizyklen).
"""
from .konstanten import *
from .formeln import formel
from .rueckkehr import ell0, chi


@formel('neigung_xy', 'g = x/y²·g⊕ (A1)',
        ('g',), ('x', 'y', 'g⊕'), 'L', 'modell', 'Autor (FTM-PDF, A1)', ('01-grundlagen.md', '(A1)  g(x,y) = x/y² · g⊕'), triade='Wirbel')
def neigung_xy(x, y):
    """(A1) g(x, y) = x/y²·g⊕, x = M/M⊕, y = R/R⊕."""
    return x / y**2 * g_earth


@formel('ds', 'Δs = GM/(3c²) (A2)',
        ('Δs',), ('M',), 'L', 'modell', 'Autor (FTM-PDF, A2)', ('01-grundlagen.md', '(A2)  Δs     = GM/(3c²)'), triade='Zentrum')
def ds(M):
    """(A2) Faltungstiefe Δs = GM/(3c²)."""
    return G * M / (3 * c**2)


@formel('neigung', 'g = GM/R² = 3c²Δs/R²',
        ('g',), ('M', 'R'), 'L', 'modell', 'Autor (FTM-PDF)', ('01-grundlagen.md', 'g·R² = GM = 3c²·Δs'), alternative='neigung_xy',
        triade='Wirbel')
def neigung(M, R):
    return G * M / R**2


@formel('takt', 'dτ/dt = √(1 − 6Δs/r) (A3)',
        ('dτ/dt',), ('M', 'r'), 'L', 'modell', 'Autor (FTM-PDF, A3)', ('01-grundlagen.md', '(A3)  dτ/dt  = √(1 − 2GM/(rc²)) = √(1 − 6Δs/r)'))
def takt(M, r):
    """(A3) dτ/dt = √(1 − 6Δs/r)."""
    return math.sqrt(1 - 6 * ds(M) / r)


@formel('delta', 'δ = 1 − dτ/dt',
        ('δ',), ('M', 'r'), 'L', 'modell', 'Autor (FTM-PDF)', ('01-grundlagen.md', 'δ = 1 − dτ/dt ≈ 3Δs/r'))
def delta(M, r):
    """δ = 1 − √(1 − x) mit x = 6Δs/r, ohne Auslöschung als x/(1 + √(1 − x)): genau bis δ ≈ 10⁻³⁹ (1 fm neben dem
    Proton), wo 1 − dτ/dt in Gleitkomma 0 gäbe."""
    x = 6 * ds(M) / r
    return x / (1 + math.sqrt(1 - x))


@formel('r0', 'r₀ = 6Δs',
        ('r₀',), ('M',), 'L', 'modell', 'Autor (FTM-PDF)', ('01-grundlagen.md', 'r₀ = 6Δs:  dτ/dt = 0'), triade='Zentrum')
def r0(M):
    return 6 * ds(M)


@formel('r_stern', 'r* = √(GM/a₀)',
        ('r*',), ('M', 'a₀'), 'L', 'modell', 'Autor (FTM-PDF)', ('08-licht-und-trichter.md', 'r* = √(GM/a₀) = c·√(3Δs/a₀)'))
def r_stern(M, a0=A0):
    return math.sqrt(G * M / a0)


@formel('r_stern_leiter', 'r*² = 3Δs·ℓ₀',
        ('r*',), ('M', 'a₀'), 'LI', 'modell', 'Baustein 08', ('08-licht-und-trichter.md', 'r*² = 3Δs·ℓ₀'), alternative='r_stern')
def r_stern_leiter(M, a0=A0):
    """Leiter: r* als geometrisches Mittel der Faltung 3Δs (L) und der Länge der Wirkung ℓ₀ (I)."""
    return math.sqrt(3 * ds(M) * ell0(a0))


# ------------------------------------------------------------------ Leiter der Ebene: 3Δs : r* : ℓ₀ (08)
@formel('leiter_schritt', 'δ* = 3Δs/r* = r*/ℓ₀ = (v_flach/c)²',
        ('δ*',), ('M', 'a₀'), 'LI', 'modell', 'Baustein 08',
        ('08-licht-und-trichter.md', 'Diese eine Zahl δ* = (v_flach/c)²'))
def leiter_schritt(M, a0=A0):
    """Schritt der geometrischen Leiter Faltung 3Δs (L) : Schwellenradius r* : Wirkungslänge ℓ₀ (I) jeder Ebene."""
    return 3 * ds(M) / r_stern(M, a0)


# ------------------------------------------------------------------ Glättung und Hierarchie der Trichter (07)
@formel('gezeit_betrag', '𝒯 = √(3/2)·|Ē^TF|',
        ('𝒯',), ('Ē^TF',), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'ℓ_i³ = GM_i / (√(3/2)·|⟨Ē^TF⟩_i|)'),
        triade='Hintergrund')
def gezeit_betrag(E):
    """Stärke der Gezeit des geglätteten Elternfelds aus seiner Hesse-Form (3×3), Spurteil entfernt, Frobenius-Norm.
    Das Feld der Elternebene ist der Hintergrund der Kindebene."""
    E = np.asarray(E, float)
    tf = E - np.trace(E) / 3 * np.eye(3)
    return math.sqrt(1.5) * float(np.sqrt(np.sum(tf * tf)))


@formel('glaettungslaenge', 'ℓ³ = GM/𝒯',
        ('ℓ',), ('M', '𝒯'), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'Die Glättungswolke mit ℓ³ = GM/𝒯'))
def glaettungslaenge(M, T):
    """Glättungslänge eines kompakten Körpers in der Gezeit 𝒯 seiner Elternebene: bis ℓ bleibt er Newtonsch."""
    return (G * M / T)**(1 / 3)


@formel('gezeitendichte', 'ρ_J = 3𝒯/(4πG)',
        ('ρ_J',), ('𝒯',), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'ρ_J = 3𝒯/(4πG)'))
def gezeitendichte(T):
    """Schwelle der Glättung als mittlere Dichte der Glättungswolke; geglättet wird nur der Überschuss (ρ − ρ_J)₊.
    Eine Schwelle der Regel, kein Stoff."""
    return 3 * T / (4 * math.pi * G)


@formel('schnittmasse_glaettung', 'M× = a₀³/(G𝒯²)',
        ('M×',), ('𝒯', 'a₀'), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'M× = a₀³/(G𝒯²)'))
def schnittmasse_glaettung(T, a0=A0):
    """Masse, an der die Glättungsgrenze ℓ = R die Schwelle r* = R schneidet."""
    return a0**3 / (G * T**2)


@formel('schnittradius_glaettung', 'R× = a₀/𝒯',
        ('R×',), ('𝒯', 'a₀'), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'R× = ℓ = r* = a₀/𝒯'))
def schnittradius_glaettung(T, a0=A0):
    return a0 / T


@formel('glaettung_zu_schwelle', 'ℓ/r* = (M×/M)^{1/6}',
        ('ℓ/r*',), ('M', '𝒯', 'a₀'), 'L', 'modell', 'Baustein 07', ('07-feldgleichung.md', 'ℓ/r* = (M×/M)^{1/6}'))
def glaettung_zu_schwelle(M, T, a0=A0):
    return (schnittmasse_glaettung(T, a0) / M)**(1 / 6)


# ------------------------------------------------------------------ Flussgesetz mit der Regel (07)
@formel('flussdimension', 'd_eff = 1 + (d − 1)(1 − χ/2); Kugel (d = 3): 1 − d ln a/d ln r = 3 − χ',
        ('d_eff',), ('a_N', 'a₀', 'd'), 'L', 'modell', 'Baustein 03, 07',
        ('03-von-der-galaxie-bis-zum-kosmos.md', 'd_eff = 1 + (d − 1)(1 − χ/2)'))
def flussdimension(aN, a0=A0, d=3):
    """Flussdimension einer Form der Dimension d unter der Regel (Kugel 3 → 2, Filament 2 → 1,5, Wand 1). Newton: Fluss
    durch die Kugelfläche (3); Tiefbereich: durch einen Umfang (2); χ zählt die verlorene Dimension."""
    return 1 + (d - 1) * (1 - chi(aN, a0) / 2)


@formel('epizyklen_summe', '(κ_E² + ω_z²)/Ω² = 2 + χ',
        ('(κ_E²+ω_z²)/Ω²',), ('a_N', 'a₀'), 'L', 'modell', 'Baustein 07',
        ('07-feldgleichung.md', 'mit der Regel wird daraus (2 + χ)Ω²'))
def epizyklen_summe(aN, a0=A0):
    """Kreisbahn um eine Punktmasse mit der örtlichen Regel; leerer Raum: 2. Den Überschuss χΩ² liest ein Beobachter,
    der den leeren Pythagoras erwartet, als Dichte; er ist Regel, kein Stoff."""
    return 2 + chi(aN, a0)


# ------------------------------------------------------------------ Zwei Kanäle der Gravitation (05), Rückgabe der Falte (04)
@formel('kanal_verhaeltnis', 'R_K = f_g·(r − 6Δs)/(3Δs)',
        ('R_K',), ('M', 'r', 'f_g'), 'L', 'modell', 'Baustein 05',
        ('05-zustandsraum.md', 'ist der differentielle Kanal 1,44·10⁹-mal stärker'))
def kanal_verhaeltnis(M, r, f_g):
    """ω_m,beob = √(k_el/m ± 3c²Δs/(ℓr²))·√(1 − 6Δs/r) bei radialer Verschiebung: differentieller Kanal
    ½f_g·|d ln g/dr| = f_g/r gegen Takt-Kanal d ln√(1 − 6Δs/r)/dr; f_g Anteil der Gravitation an der Rückstellung."""
    D = ds(M)
    return f_g * (r - 6 * D) / (3 * D)
