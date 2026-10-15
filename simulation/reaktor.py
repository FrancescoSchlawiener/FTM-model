"""Rechnungen zum Schauer-Faltungs-Fusionsreaktor (SFFR): die Regel im Labor, Energien der Kernreaktionen aus AME2020,
Feld der Torusspule, Schumann-Kopplung. Lauf: python simulation/reaktor.py gibt alle Zahlen des Berichts."""
import math

import kern as K
import daten as D

U_MEV = 931.494_102_42                 # MeV je u (CODATA 2018)


def bindung():
    """Gesamtbindung [MeV] je (Z, A) aus AME2020."""
    return {(Z, A): b * A / 1e3 for Z, N, A, b, _ in D.ame2020()}


def q_wert(vorher, nachher):
    """Q [MeV] einer Kernreaktion aus den Bindungen: Summe nachher minus vorher (das freie Neutron bindet 0)."""
    B = bindung()
    return sum(B.get(k, 0.0) for k in nachher) - sum(B.get(k, 0.0) for k in vorher)


REAKTIONEN = {
    "D(n,γ)T": ([(1, 2)], [(1, 3)]),
    "⁶Li(n,α)T": ([(3, 6)], [(2, 4), (1, 3)]),
    "⁷Li(n,γ)⁸Li": ([(3, 7)], [(3, 8)]),
    "³He(n,p)T": ([(2, 3)], [(1, 1), (1, 3)]),
    "D + T → α + n": ([(1, 2), (1, 3)], [(2, 4)]),
    "D(γ,n)p": ([(1, 2)], [(1, 1)]),
    "⁹Be(γ,n)⁸Be": ([(4, 9)], [(4, 8)]),
}


def regel_im_labor(M=1000.0, r=1.0):
    """Masse M [kg] im Abstand r [m] und das Erdfeld: a_N/a₀, s, ν − 1 (als lg) und r* der Masse [m]."""
    aus = {}
    for name, aN in (("Reaktor", K.G * M / r**2), ("Erde", 9.81)):
        s = math.sqrt(aN / K.a0)
        aus[name] = dict(aN=aN, verh=aN / K.a0, s=s, lg_nu1=-s / math.log(10))
    aus["r_stern"] = math.sqrt(K.G * M / K.a0)
    return aus


def feldenergie(B=10.0, V=1.0):
    """Energie [J] eines Magnetfelds B [T] im Volumen V [m³], ihre Masse E/c² [kg] und deren Neigung in 1 m [m/s²]."""
    E = B**2 / (2 * K.mu0) * V
    m = E / K.c**2
    return E, m, K.G * m


def neutronen_bilanz():
    """Q je Reaktion [MeV]; Kreislauf Photospaltung des Deuterons plus Einfang auf D (= D + D → p + T)."""
    q = {n: q_wert(*r) for n, r in REAKTIONEN.items()}
    q["Kreislauf D(γ,n)p + D(n,γ)T"] = q["D(γ,n)p"] + q["D(n,γ)T"]
    return q


def leistung(phi=1e14, n_li6=4.6e28, sigma_b=940.0, q_mev=None):
    """Leistungsdichte [W/m³] aus Fluss φ [n/(cm²·s)], Kerndichte [1/m³], Wirkungsquerschnitt [b] und Q [MeV]."""
    q_mev = neutronen_bilanz()["⁶Li(n,α)T"] if q_mev is None else q_mev
    return phi * 1e4 * n_li6 * sigma_b * 1e-28 * q_mev * 1e6 * K.e


def schumann(f=7.83, l=1000.0):
    """Wellenlänge [km] bei f [Hz] und der Strahlungswiderstand [Ω] eines kurzen Dipols der Länge l [m], 80π²(l/λ)²."""
    lam = K.c / f
    return lam / 1e3, 80 * math.pi**2 * (l / lam)**2


def neutron_welle(v=2200.0):
    """De-Broglie-Wellenlänge [m] eines Neutrons der Geschwindigkeit v [m/s] (thermisch 2200 m/s) gegen den Kernradius
    λ_H·A^{1/3} des Deuterons [m] (λ_H = 1,25 fm wie in der Vorlage)."""
    lam = 2 * math.pi * K.hbar / (K.m_n * v)
    return lam, 1.25e-15 * 2**(1 / 3)


def compton_takt():
    """Compton-Frequenz des Nukleons m_n·c²/h [Hz] gegen einen Takt der Leistungselektronik von 1 GHz."""
    f = K.m_n * K.c**2 / (2 * math.pi * K.hbar)
    return f, f / 1e9


def neutronen_je_megawatt(q_mev=None):
    """Eingefangene Neutronen je Sekunde für 1 MW bei ⁶Li(n,α)T."""
    q_mev = neutronen_bilanz()["⁶Li(n,α)T"] if q_mev is None else q_mev
    return 1e6 / (q_mev * 1e6 * K.e)


def thomson():
    """Streuquerschnitt eines Photons an einem Elektron, (8π/3)·r_e² [barn]."""
    r_e = K.e**2 / (4 * math.pi * K.eps0 * K.m_e * K.c**2)
    return 8 * math.pi / 3 * r_e**2 / 1e-28


def fermium():
    """Halbwertszeiten [s] und Zerfall der Fermium-Isotope 255 bis 260 (NUBASE2020): wo der Einfang mit β-Zerfällen endet."""
    nb = D.nubase()
    return {A: nb.get((100, A)) for A in range(255, 261)}


def li6_aufteilung():
    """Energie [MeV] von α und Triton aus ⁶Li(n,α)T (thermisch, Impulse gleich): E_α = Q·m_T/(m_α + m_T)."""
    q = neutronen_bilanz()["⁶Li(n,α)T"]
    return q * 3 / 7, q * 4 / 7


def torus_feld(N=1000, I=5000.0, r=1.0):
    """Feld [T] im Innern einer Torusspule mit N Windungen, Strom I [A], Radius r [m]: B = μ₀NI/(2πr), außen null."""
    return K.mu0 * N * I / (2 * math.pi * r)


REAKTIONEN["⁹Be(n,2n)2α"] = ([(4, 9)], [(2, 4), (2, 4)])


def masse(Z, A):
    """Kernmasse [MeV] aus der Bindung (AME2020)."""
    return Z * 938.272_088 + (A - Z) * 939.565_421 - bindung()[(Z, A)]


def huelle(m_teilchen):
    """Größe der Hülle um das Deuteron, ħ/(m_red·α·c) [m], mit der reduzierten Masse des umlaufenden Teilchens."""
    m_d = masse(1, 2) * 1e6 * K.e / K.c**2
    m_red = m_teilchen * m_d / (m_teilchen + m_d)
    return K.hbar / (m_red * K.alpha * K.c)


def tunnel(abstand, E_kev=0.0):
    """Durchgang durch den Coulombwall von D und T aus dem Abstand [m] bis zur Berührung R = λ_H·(2^{1/3} + 3^{1/3}),
    λ_H = 1,25 fm, bei der Energie E [keV] gegen den nackten Wall: lg der Wahrscheinlichkeit nach WKB."""
    from scipy.integrate import quad
    mu = masse(1, 2) * masse(1, 3) / (masse(1, 2) + masse(1, 3)) * 1e6 * K.e / K.c**2
    R = 1.25e-15 * (2**(1 / 3) + 3**(1 / 3))
    E = E_kev * 1e3 * K.e
    k = K.e**2 / (4 * math.pi * K.eps0)
    r_aus = abstand if E <= 0 else min(abstand, k / E)
    S = quad(lambda r: math.sqrt(max(2 * mu * (k / r - E), 0.0)) / K.hbar, R, r_aus, limit=400)[0]
    return -2 * S / math.log(10)


def kreisrate():
    """Aus Y = L/(1 + ω·L), L = λ_c·τ_μ, mit der gemessenen Ausbeute und Haftung: die Kreisrate λ_c [1/s] und die Grenze 1/ω."""
    L = 1 / (1 / K.MCF_FUSIONEN - K.MCF_HAFTEN)
    return L / K.TAU_MYON, 1 / K.MCF_HAFTEN


def energie_je_verschmelzung():
    """Energie [MeV] je Verschmelzung: allein, mit Einfang des Neutrons auf ⁶Li, mit Vermehrung über ⁹Be(n,2n) und zwei
    Einfängen auf ⁶Li; Energie von α und Neutron; geladener Anteil im Mantel aus ⁹Be und ⁶Li."""
    q = neutronen_bilanz()
    allein = q["D + T → α + n"]
    li = allein + q["⁶Li(n,α)T"]
    be = allein + q["⁹Be(n,2n)2α"] + 2 * q["⁶Li(n,α)T"]
    m_n, m_a = 939.565_421, masse(2, 4)
    e_a = allein * m_n / (m_n + m_a)
    return dict(allein=allein, li=li, be=be, e_alpha=e_a, e_neutron=allein - e_a)


def gewinnschwelle(kosten_gev, wirkungsgrad=0.4):
    """Verschmelzungen je Myon, ab denen das Myon der Kosten [GeV, elektrisch] mehr Strom liefert als es kostet (Mantel Be, ⁶Li)."""
    return kosten_gev * 1e3 / (wirkungsgrad * energie_je_verschmelzung()["be"])


def myonen_fuer(leistung_w, fusionen=None):
    """Myonen je Sekunde für eine thermische Leistung [W] (Mantel Be, ⁶Li)."""
    f = K.MCF_FUSIONEN if fusionen is None else fusionen
    return leistung_w / (f * energie_je_verschmelzung()["be"] * 1e6 * K.e)


def abschalten(t_us=10.0):
    """Anteil der Myonen, die t [µs] nach dem Abschalten des Strahls noch leben."""
    return math.exp(-t_us * 1e-6 / K.TAU_MYON)


def larmor(E_mev, B=5.0, Z=2, A=4):
    """Bahnradius [m] eines Kerns (Z, A) der Energie E [MeV] im Feld B [T]."""
    m = masse(Z, A) * 1e6 * K.e / K.c**2
    return math.sqrt(2 * m * E_mev * 1e6 * K.e) / (Z * K.e * B)


def tritium_halbwert():
    """Halbwertszeit des Tritiums [a] (NUBASE2020)."""
    return D.nubase()[(1, 3)][0] / K.yr


def spaltung():
    """Q [MeV] der Spaltung n + ²³³U → ¹⁴⁰Ce + ⁹²Zr + 2n bis zu den stabilen Enden (AME2020), einschließlich der Antineutrinos."""
    return q_wert([(92, 233)], [(58, 140), (40, 92)])


def mantel(k=0.95, e_spalt=180.0):
    """Unterkritischer Thorium-Mantel mit Vermehrung k: Spaltungen je Quellneutron k/(ν(1 − k)), alle Neutronen je Quellneutron
    1/(1 − k), frei für Lithium und Verluste nach Spaltung und Nachbrüten des ²³³U aus ²³²Th, Energie je Verschmelzung [MeV] mit dem Auslegungswert e_spalt [MeV] je Spaltung (ohne Antineutrinos)."""
    F = k / (K.NU_U233 * (1 - k))
    return dict(spaltungen=F, neutronen=1 / (1 - k), frei=1 / (1 - k) - 2 * F, energie=neutronen_bilanz()["D + T → α + n"] + F * e_spalt)


def gewinn(kosten_gev=5.0, wirkungsgrad=0.4, k=0.95, fusionen=None):
    """Strom aus einem Myon geteilt durch den Strom, den es kostet."""
    f = K.MCF_FUSIONEN if fusionen is None else fusionen
    return f * mantel(k)["energie"] * wirkungsgrad / (kosten_gev * 1e3)



def bahn_myon(p_mev=100.0, B=1.0):
    """Bahnradius [m] eines Myons des Impulses p [MeV/c] im Feld B [T]."""
    return p_mev * 1e6 / (K.c * B)


def eta_u233():
    """Neutronen je im ²³³U absorbiertem Neutron, η = ν·σ_f/(σ_f + σ_γ)."""
    return K.NU_U233 * K.U233_SPALT / (K.U233_SPALT + K.U233_FANG)


def kette(k=0.97, e_spalt=180.0, leistung_w=1e8, wirkungsgrad=None, beschleuniger=1 / 3):
    """Kette der Neutronen, angestoßen von einem Proton (Spallation, kern.SPALLATION): Verstärkung ν = 1/(1 − k) = Σkⁿ, Spaltungen
    je Quellneutron und je Proton, Wärme je Proton [GeV], Gewinn Wärme durch Strahlenergie, kleinster Wirkungsgrad des
    Beschleunigers, ab dem der Reaktor seinen Strahl selbst trägt, Strahlstrom [mA] und Strahlleistung [MW] für die Leistung,
    Strom [MW] und Eigenbedarf des Beschleunigers [MW] beim Wirkungsgrad beschleuniger."""
    wirkungsgrad = K.WIRKUNGSGRAD_THERMOAKUSTIK if wirkungsgrad is None else wirkungsgrad
    E_p, Y = K.SPALLATION
    F = k / (K.NU_U233 * (1 - k))
    waerme = Y * F * e_spalt / 1e3
    protonen = leistung_w / (waerme * 1e9 * K.e)
    return dict(nu=1 / (1 - k), spaltungen=F, je_proton=Y * F, waerme=waerme, gewinn=waerme / E_p,
                wirkungsgrad_min=E_p / (wirkungsgrad * waerme), strom_ma=protonen * K.e * 1e3,
                strahl_mw=protonen * E_p * 1e9 * K.e / 1e6, strom_mw=wirkungsgrad * leistung_w / 1e6,
                eigenbedarf_mw=protonen * E_p * 1e9 * K.e / 1e6 / beschleuniger)


def haushalt(k=0.97, anteil_li=0.5):
    """Neutronen je Generation: Anteil, der im ²³³U absorbiert wird (k/η), der ²³²Th zu ²³³U aufbaut (ebenso viel, Brutrate 1),
    der Rest; davon anteil_li im ⁶Li; Einfänge im ⁶Li je Quellneutron."""
    p_f = k / eta_u233()
    rest = 1 - 2 * p_f
    li = anteil_li * rest / (1 - k)
    return dict(uran=p_f, thorium=p_f, rest=rest, li_je_quelle=li,
                li_mev_je_proton=K.SPALLATION[1] * li * neutronen_bilanz()["⁶Li(n,α)T"])


def moderator(A=12, sigma_streu=None, sigma_abs=None, dichte=1.7, E0=2e6, v=2200.0):
    """Bremsen und Halten der Neutronen in einem Moderator der Massenzahl A und Dichte [g/cm³] (Graphit nach Sears): Stöße
    bis thermisch, Bremszeit 2λ_s/(ξ·v) [s], Diffusionslänge L [cm], Albedo der dicken Wand (1 − 2D/L)/(1 + 2D/L), Lebensdauer
    eines thermischen Neutrons im reinen Moderator [s]."""
    sigma_streu = K.QUERSCHNITT["C_streu"] if sigma_streu is None else sigma_streu
    sigma_abs = K.QUERSCHNITT["C_abs"] if sigma_abs is None else sigma_abs
    xi = 1 + (A - 1)**2 / (2 * A) * math.log((A - 1) / (A + 1))
    N = dichte * 6.022_140_76e23 / A
    S_s = N * sigma_streu * (A / (A + 1))**2 * 1e-24
    S_a = N * sigma_abs * 1e-24
    D_ = 1 / (3 * S_s * (1 - 2 / (3 * A)))
    L = math.sqrt(D_ / S_a)
    return dict(xi=xi, stoesse=math.log(E0 / 0.0253) / xi, bremszeit=2 / (S_s * xi * v * 100), L=L,
                albedo=(1 - 2 * D_ / L) / (1 + 2 * D_ / L), leben=1 / (S_a * v * 100))


def einfangflaeche(A=233, sigma=None):
    """Fläche des gehaltenen Trichters π(λ_H·A^{1/3})² [b] gegen den gemessenen Querschnitt [b] und die Wellenlänge ħ/(m_n·v)
    des thermischen Neutrons gegen den Radius."""
    sigma = K.U233_SPALT + K.U233_FANG if sigma is None else sigma
    R = 1.25e-15 * A**(1 / 3)
    return dict(geometrisch=math.pi * R**2 / 1e-28, gemessen=sigma, verhaeltnis=sigma / (math.pi * R**2 / 1e-28),
                welle_durch_radius=K.hbar / (K.m_n * 2200.0) / R)


def paar_einfang():
    """Energie [MeV], die der Einfang eines Neutrons frei macht (AME2020): in ²³³U schließt er ein Paar, in ²³²Th öffnet er eines;
    dazu die Halbwertszeiten von ²³³Th und ²³³Pa [s] (NUBASE2020) auf dem Weg ²³²Th → ²³³U."""
    B = bindung()
    nb = D.nubase()
    return dict(u233=B[(92, 234)] - B[(92, 233)], th232=B[(90, 233)] - B[(90, 232)],
                th233=nb[(90, 233)][0], pa233=nb[(91, 233)][0])


def abklingen(k=0.97, anteil_graphit=0.05):
    """Generationszeit [s] im Kern (Bremszeit plus Leben bis zur Absorption, wenn der Graphit anteil_graphit der Neutronen
    absorbiert) und Abklingzeit der Kette nach dem Abschalten des Strahls, l/(1 − k) [s]; Abklingrate α = (1 − k)/l [1/s]."""
    m = moderator()
    l = m["bremszeit"] + anteil_graphit * m["leben"]
    return dict(generation=l, abklingzeit=l / (1 - k), alpha=(1 - k) / l)


def rueckkopplung(k=0.97, dichte_weniger=0.01):
    """q nach Ausdehnung des flüssigen Brennstoffs um den Anteil dichte_weniger: Uran und Thorium verlieren denselben Anteil,
    der Rest der Absorption bleibt. Gibt (q neu, ν neu)."""
    h = haushalt(k)
    f = 1 - dichte_weniger
    k_neu = k * f / ((h["uran"] + h["thorium"]) * f + h["rest"])
    return k_neu, 1 / (1 - k_neu)


def torus_dipol(I=5000.0, R=1.0, d=10.0):
    """Fernfeld [T] auf der Achse im Abstand d [m] einer einfach gewickelten Torusspule: Ihre Windungen schreiten einmal um
    den Ring fort und bilden eine Schleife mit dem Moment I·πR²; die gegensinnige bifilare Lage hebt es auf."""
    return K.mu0 * I * math.pi * R**2 / (2 * math.pi * d**3)


def kuramoto(sigma_grad=1.0):
    """Ordnung r = |⟨e^{iθ}⟩| vieler Kavitäten mit gaußisch gestreuter Phase σ: r = e^{−σ²/2}, der Anteil des Energiegewinns."""
    return math.exp(-math.radians(sigma_grad)**2 / 2)


def thorium_verbrauch(e_spalt=180.0):
    """Gespaltenes und insgesamt aufgenommenes ²³³U je Gigawattjahr Wärme [kg]; bei Brutrate 1 baut ²³²Th das aufgenommene nach."""
    spaltungen = 1e9 * K.yr / (e_spalt * 1e6 * K.e)
    m = spaltungen * 233 * K.m_u
    return m, m * (K.U233_SPALT + K.U233_FANG) / K.U233_SPALT


if __name__ == "__main__":
    print("Li-6: E_alpha, E_T [MeV]", li6_aufteilung(), "Torus [T]", torus_feld())
    print("Neutronwelle, Kernradius D [m]", neutron_welle())
    print("Compton-Takt [Hz], Verhältnis zu 1 GHz", compton_takt())
    print("Neutronen je MW [1/s]", neutronen_je_megawatt())
    print("Thomson [b]", thomson())
    print("Fermium", fermium())
    for k, v in regel_im_labor().items():
        print(k, v)
    print("Feld 10 T in 1 m³: E, m, g", feldenergie())
    for k, v in neutronen_bilanz().items():
        print(f"{k:32s} {v:8.3f} MeV")
    print("Leistung Li-6 bei 1e14 n/cm²s [W/m³]", leistung())
    print("Schumann", schumann())
    a_e, a_m = huelle(K.m_e), huelle(K.M_MYON)
    print("Hülle Elektron, Myon [m], Verhältnis", a_e, a_m, a_e / a_m)
    print("Tunnel lg P aus 1,4·Hülle: Elektron, Myon", tunnel(1.4 * a_e), tunnel(1.4 * a_m))
    print("Kreisrate [1/s], Grenze 1/ω", kreisrate())
    print("Energie je Verschmelzung", energie_je_verschmelzung())
    for kosten in (5.0, 3.0, 2.0):
        print("Gewinnschwelle", kosten, "GeV, Wirkungsgrad 0,4:", gewinnschwelle(kosten))
    print("Myonen je s für 1 GW", myonen_fuer(1e9), "Strahl [µA]", myonen_fuer(1e9) * K.e * 1e6)
    print("nach 10 µs", abschalten(), "Larmor α 3,5 MeV 5 T [m]", larmor(energie_je_verschmelzung()["e_alpha"]))
    print("Tritium [a]", tritium_halbwert())
    print("Spaltung Q [MeV]", spaltung())
    for k in (0.90, 0.95, 0.97):
        print("Mantel k", k, mantel(k), "Gewinn 5 GeV", gewinn(k=k), "Gewinn 5 GeV, 100 Fusionen", gewinn(k=k, fusionen=100))
    m = mantel()["energie"]
    print("Myonen je s für 1 GW thermisch", 1e9 / (K.MCF_FUSIONEN * m * 1e6 * K.e), "Strahl [mA]", 1e9 / (K.MCF_FUSIONEN * m * 1e6 * K.e) * K.e * 1e3)
    print("Myonbahn 100 MeV/c in 1 T [m]", bahn_myon())
    print("η(²³³U)", eta_u233())
    for k in (0.95, 0.97, 0.98):
        print("Kette", k, kette(k), haushalt(k))
    print("Moderator", moderator(), "Einfangfläche", einfangflaeche(), "⁶Li", einfangflaeche(6, K.QUERSCHNITT["Li6"]))
    print("Paar", paar_einfang(), "Abklingen", abklingen(), "Rückkopplung 1 %, 10 %", rueckkopplung(), rueckkopplung(dichte_weniger=0.1))
    print("Torusdipol [T]", torus_dipol(), "Kuramoto 1°, 10°", kuramoto(), kuramoto(10.0), "Thorium kg/GWa", thorium_verbrauch())
    print("Larmor T 2,73 MeV und α 2,05 MeV bei 1 T [m]", larmor(2.73, 1.0, 1, 3), larmor(2.05, 1.0))
