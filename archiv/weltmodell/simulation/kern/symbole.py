# -*- coding: utf-8 -*-
"""
Symbol-Register: jedes Zeichen des Modells genau einmal, mit Code-Name, Bedeutung, SI-Einheit, Ebene, Zustand
und Herkunft. Kein Zeichen hat zwei Bedeutungen; verschiedene Sachen tragen verschiedene Zeichen (w_Z, Q_el …).
Dieselbe Sache auf verschiedenen Ebenen trägt dasselbe Zeichen mit Ebenen-Index (Selbstähnlichkeit der Triade):
ein Eintrag mit stamm = Grundzeichen, Zeichen = Grundzeichen_Index (q_Δs, q_kosm zu q; τ_M, τ_Z zu τ).

Ebene: L lokal am Ort (Raum), I innere Richtung (Rückkehr, orthogonal zum Raum), Z Zustandsraum (zeitfeste
Formverhältnisse), G das Ganze (Torus von außen, Zelle). Brücken tragen mehrere Buchstaben (z. B. "LI");
Naturkonstanten und Einheiten gelten auf allen Ebenen ("LIZG").
Zustand: fest | seit_z10 (fest seit z ≈ 10, früher anders) | epoche (atmet mit, je Epoche und Beobachter) |
je_system (je Körper, Bereich oder Aufbau).
"""

EBENEN = set("LIZG")
ZUSTAENDE = {"fest", "seit_z10", "epoche", "je_system"}


class Symbol:
    __slots__ = ("zeichen", "code_name", "bedeutung", "einheit_si", "ebene", "zustand", "herkunft", "stamm")

    def __init__(self, zeichen, code_name, bedeutung, einheit_si, ebene, zustand, herkunft, stamm=None):
        if not zeichen or not code_name:
            raise ValueError("Symbol ohne Zeichen oder Code-Name")
        if not ebene or not set(ebene) <= EBENEN:
            raise ValueError(f"{zeichen}: Ebene {ebene!r} nicht aus L/I/Z/G")
        if zustand not in ZUSTAENDE:
            raise ValueError(f"{zeichen}: Zustand {zustand!r} unbekannt")
        if not herkunft:
            raise ValueError(f"{zeichen}: Herkunft fehlt")
        self.zeichen, self.code_name, self.bedeutung = zeichen, code_name, bedeutung
        self.einheit_si, self.ebene, self.zustand, self.herkunft = einheit_si, ebene, zustand, herkunft
        self.stamm = stamm

    def __repr__(self):
        return f"Symbol({self.zeichen!r}, {self.code_name!r})"


class SymbolRegister(dict):
    """Zeichen → Symbol. Doppeltes Zeichen oder doppelter Code-Name wirft einen Fehler; ein Eintrag mit stamm
    verlangt das Grundzeichen im Register und das Zeichen Grundzeichen_Index."""

    def __init__(self):
        super().__init__()
        self.nach_code = {}

    def neu(self, *args, stamm=None):
        s = Symbol(*args, stamm=stamm)
        if stamm is not None and (stamm not in self or not s.zeichen.startswith(stamm + "_")):
            raise ValueError(f"{s.zeichen}: Grundzeichen {stamm!r} fehlt oder Zeichen nicht {stamm}_Index")
        if s.zeichen in self:
            raise ValueError(f"Zeichen doppelt: {s.zeichen}")
        if s.code_name in self.nach_code:
            raise ValueError(f"Code-Name doppelt: {s.code_name}")
        self[s.zeichen] = s
        self.nach_code[s.code_name] = s
        return s


SYMBOLE = SymbolRegister()
_s = SYMBOLE.neu
REGISTERMODUL = True
A, MESS, DEF = "Autor (FTM-PDF)", "Messwert [G]", "Definition SI"

# ------------------------------------------------------------------ Konstanten und Einheiten (alle Ebenen)
_s("G", "G", "Gravitationskonstante", "m³/(kg·s²)", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("c", "c", "Lichtgeschwindigkeit, Skalierungsfaktor der Hülle", "m/s", "LIZG", "fest", DEF)
_s("ħ", "hbar", "reduzierte Planck-Konstante", "J·s", "LIZG", "fest", DEF)
_s("k_B", "kB", "Boltzmann-Konstante", "J/K", "LIZG", "fest", DEF)
_s("e", "e_el", "Elementarladung", "C", "LIZG", "fest", DEF)
_s("N_A", "N_A", "Avogadro-Konstante", "mol⁻¹", "LIZG", "fest", DEF)
_s("ε₀", "eps0", "elektrische Feldkonstante", "F/m", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("m_u", "m_u", "atomare Masseneinheit m(¹²C)/12", "kg", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("m_p", "m_p", "Protonmasse", "kg", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("m_e", "m_e", "Elektronmasse", "kg", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("α", "alpha", "Feinstrukturkonstante r_e/λ̄_e", "1", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("r_e", "r_e", "klassischer Elektronenradius", "m", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("a_B", "a_B", "Bohr-Radius", "m", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("R_p", "R_p", "Ladungsradius des Protons", "m", "LIZG", "fest", "Messwert CODATA 2018 [G]")
_s("M☉", "M_sun", "Sonnenmasse", "kg", "LIZG", "fest", "IAU (GM☉) [G]")
_s("GM☉", "GM_sun", "Gravitationsparameter der Sonne", "m³/s²", "LIZG", "fest", "IAU [G]")
_s("R☉", "R_sun", "Sonnenradius", "m", "LIZG", "fest", "IAU nominal")
_s("P_L☉", "L_sun", "Leuchtkraft der Sonne", "W", "LIZG", "fest", "IAU nominal")
_s("T_eff☉", "TEFF_sun", "Effektivtemperatur der Sonne", "K", "LIZG", "fest", "IAU nominal")
_s("M⊕", "M_earth", "Erdmasse", "kg", "LIZG", "fest", MESS)
_s("GM⊕", "GM_earth", "Gravitationsparameter der Erde", "m³/s²", "LIZG", "fest", MESS)
_s("R⊕", "R_earth", "Erdradius", "m", "LIZG", "fest", MESS)
_s("g⊕", "g_earth", "Neigung an der Erdoberfläche, 9,81 m/s²", "m/s²", "L", "fest", "Autor (FTM-PDF, A1)")
_s("AE", "AU", "Astronomische Einheit", "m", "LIZG", "fest", DEF)
_s("pc", "pc", "Parsec (kpc, Mpc, Gpc Vielfache)", "m", "LIZG", "fest", DEF)
_s("Lj", "ly", "Lichtjahr", "m", "LIZG", "fest", DEF)
_s("Jahr", "yr", "julianisches Jahr (Mrd. Jahre Vielfache)", "s", "LIZG", "fest", DEF)
_s("ℓ_P", "l_P", "Planck-Länge", "m", "LIZG", "fest", "aus ħ, G, c")
_s("m_Pl", "m_Pl", "Planck-Masse", "kg", "LIZG", "fest", "aus ħ, G, c")
_s("T_CMB", "T_CMB", "Temperatur der Hintergrundstrahlung heute", "K", "G", "epoche", MESS)

# ------------------------------------------------------------------ Ebene L: Faltung am Ort (01, 02, 08)
_s("M", "M", "Masse (Faltung eines Körpers oder Systems)", "kg", "L", "je_system", A)
_s("m", "m", "Masse eines Teilchens", "kg", "L", "je_system", "Baustein 02")
_s("R", "R", "Radius der Hülle eines Körpers", "m", "L", "je_system", A)
_s("r", "r", "Abstand vom Zentrum im Raum", "m", "L", "je_system", A)
_s("x", "x_A1", "Masse in Erdmassen, x = M/M⊕ (A1)", "1", "L", "je_system", "Autor (FTM-PDF, A1)")
_s("y", "y_A1", "Radius in Erdradien, y = R/R⊕ (A1)", "1", "L", "je_system", "Autor (FTM-PDF, A1)")
_s("g", "g", "Neigung an der Kante des Trichters, GM/R²", "m/s²", "L", "je_system", "Autor (FTM-PDF, A1)")
_s("Δs", "ds", "Faltungstiefe GM/(3c²)", "m", "L", "je_system", "Autor (FTM-PDF, A2)")
_s("dτ/dt", "takt", "Abwicklungstakt √(1 − 6Δs/r)", "1", "L", "je_system", "Autor (FTM-PDF, A3)")
_s("t_E", "t_eigen", "Eigenzeit am Ort (im Modelltext τ in dτ/dt)", "s", "L", "je_system", "Autor (FTM-PDF, A3)")
_s("δ", "delta", "Faltungsdifferential 1 − dτ/dt", "1", "L", "je_system", A)
_s("r₀", "r0", "Grenzradius 6Δs (Takt steht)", "m", "L", "je_system", A)
_s("r*", "r_stern", "Schwellenradius √(GM/a₀) = √(3Δs·ℓ₀)", "m", "L", "je_system", A)
_s("λ̄", "lambda_C", "reduzierte Compton-Länge ħ/(mc)", "m", "L", "je_system", "Baustein 02")
_s("α_G", "alpha_G", "gravitative Kopplung GM²/(ħc) = 3Δs/λ̄", "1", "L", "je_system", "Baustein 02")
_s("M_*", "M_stern", "Massenskala der Sterne m_Pl³/m_p²", "kg", "L", "fest", "Baustein 02")
_s("P_Edd", "L_Edd", "Eddington-Leuchtkraft (9/2)·m_p·c³·Δs/r_e²", "W", "L", "je_system", "Baustein 02")
_s("T", "T", "Temperatur", "K", "L", "je_system", "Baustein 04")
_s("T_H", "T_H", "Temperatur einer Falte ħc/(24πkΔs)", "K", "L", "je_system", "Baustein 04")
_s("M_T", "M_T", "Masse der Falte mit T_H = T (die Epoche kommt über T, etwa T_CMB)", "kg", "L", "je_system", "Baustein 04")
_s("t_F", "t_falte", "Lebensdauer der Falte 3840·(S/k)·Δs/c", "s", "L", "je_system", "Baustein 04")
_s("ω_C", "omega_C", "Compton-Takt Mc²/ħ", "s⁻¹", "L", "je_system", "Baustein 09")
_s("a_J", "a_J", "Verwindungslänge J/(Mc)", "m", "L", "je_system", "Baustein 08")
_s("j", "j_spin", "Spinzahl a_J/(3Δs)", "1", "L", "je_system", "Baustein 08")
_s("Q_el", "Q_el", "Ladung n·e (Windungszahl)", "C", "L", "fest", "Baustein 02")
_s("α_abl", "alpha_abl", "Ablenkwinkel des Lichts 12Δs/b", "rad", "L", "je_system", "Baustein 08")
_s("α_∞", "alpha_inf", "kollektive Ablenkung 6πΔs/r* = 2π(v_flach/c)²", "rad", "L", "je_system", "Baustein 08")
_s("b", "b_stoss", "Stoßparameter des Lichts", "m", "L", "je_system", "Baustein 08")
_s("Δφ_P", "dphi_perihel", "Periheldrehung je Umlauf", "rad", "L", "je_system", "Baustein 08")
_s("a_b", "a_bahn", "große Halbachse einer Bahn", "m", "L", "je_system", "Baustein 08")
_s("e_B", "e_bahn", "Exzentrizität einer Bahn", "1", "L", "je_system", "Baustein 08")
_s("α_S", "alpha_schatten", "Winkeldurchmesser des Schattens 2b/D_A, b = 3√27·Δs", "rad", "L", "je_system", "Baustein 08")
_s("D_A", "D_abstand", "Abstand des Beobachters", "m", "L", "je_system", "Baustein 08")
_s("r̃_ISCO", "r_isco", "innerste stabile Kreisbahn in Einheiten 3Δs", "1", "L", "je_system", "Baustein 08")
_s("f*_gw", "f_gw_stern", "Schwellenfrequenz der Gravitationswellen Ω(r*)/π", "Hz", "L", "je_system", "Baustein 08")
_s("P_b", "P_b", "Umlaufzeit eines Doppelsterns", "s", "L", "je_system", "Baustein 08")
_s("Ṗ_b", "dPb_dt", "Abnahme der Umlaufzeit", "1", "L", "je_system", "Baustein 08")
_s("m₁", "m1", "Masse der ersten Komponente (im Kern in M☉)", "kg", "L", "je_system", "Baustein 08")
_s("m₂", "m2", "Masse der zweiten Komponente (im Kern in M☉)", "kg", "L", "je_system", "Baustein 08")
_s("η_m", "eta_m", "Massenverhältnis m₁m₂/M²", "1", "L", "je_system", "Baustein 08")
_s("η", "eta", "Massenverhältnis M_dyn/M_lens", "1", "L", "je_system", "Baustein 05")
_s("ρ_J", "rho_J", "Gezeitendichte 3𝒯/(4πG)", "kg/m³", "L", "je_system", "Baustein 07")
_s("ℓ", "ell_glatt", "Glättungslänge, ℓ³ = GM/𝒯", "m", "L", "je_system", "Baustein 07")
_s("𝒯", "T_gezeit", "Gezeit des geglätteten Elternfelds √(3/2)·|Ē^TF|", "s⁻²", "L", "je_system", "Baustein 07")
_s("σ", "sigma_falt", "Faltungspotential (kovariante Fassung)", "1", "L", "je_system", "Baustein 07")
_s("σ_v", "sigma_v", "Geschwindigkeitsstreuung", "m/s", "L", "je_system", "Baustein 03")
_s("σ_el", "sigma_el", "elektrische Leitfähigkeit", "S/m", "L", "je_system", "Baustein 05")
_s("f_g", "f_g", "Anteil der Schwere an der Steifigkeit", "1", "L", "je_system", "Baustein 05")
_s("Λ_QCD", "Lambda_QCD", "Skala der starken Kopplung", "J", "L", "fest", "Baustein 02")

# ------------------------------------------------------------------ Ebene L mit der Regel (03, 07)
_s("a", "a", "Beschleunigung nach der örtlichen Regel ν(a_N)·a_N", "m/s²", "L", "je_system", "Baustein 03")
_s("a_N", "a_N", "Newton-Beschleunigung |g_N| am Ort", "m/s²", "L", "je_system", "Baustein 03")
_s("g_N", "g_N", "Newton-Feld (Vektor) als Funktion des Orts", "m/s²", "L", "je_system", "Baustein 07")
_s("a_i", "a_i", "inneres Feld eines Systems", "m/s²", "L", "je_system", "Baustein 03")
_s("a_e", "a_e", "äußeres Feld am Ort eines Systems", "m/s²", "L", "je_system", "Baustein 03")
_s("v", "v", "Kreisgeschwindigkeit", "m/s", "L", "je_system", "Baustein 03")
_s("v_flach", "v_flach", "flache Rotationsgeschwindigkeit (GMa₀)^{1/4}", "m/s", "L", "je_system", "Baustein 03")
_s("Σ_a₀", "sigma_a0", "Grenzflächendichte a₀/(πG)", "kg/m²", "L", "seit_z10", "Baustein 03")
_s("M(<r)", "M_innen", "eingeschlossene Masse (Hernquist-Profil)", "kg", "L", "je_system", "Baustein 08")
_s("b_H", "b_H", "Skalenlänge des Hernquist-Profils", "m", "L", "je_system", "Baustein 08")
_s("R_zyl", "R_zyl", "Zylinderradius in der Scheibe", "m", "L", "je_system", "Baustein 07")
_s("z_h", "z_hoehe", "Höhe über der Mittelebene", "m", "L", "je_system", "Baustein 07")
_s("v_R", "v_R", "Radialgeschwindigkeit", "m/s", "L", "je_system", "Baustein 07")
_s("v_z", "v_z", "Vertikalgeschwindigkeit", "m/s", "L", "je_system", "Baustein 07")
_s("J_z", "J_z", "Drehimpuls je Masse um die Achse", "m²/s", "L", "je_system", "Baustein 07")
_s("t", "t", "Zeit längs der Bahn", "s", "L", "je_system", "Baustein 07")
_s("ω_orb", "omega_orb", "Umlauffrequenz √((4πG/3)ρ̄)", "s⁻¹", "L", "je_system", "Baustein 03")
_s("λ_m", "lambda_m", "Übergaberate der Mode ξ, ω_orb/N*", "s⁻¹", "L", "je_system", "Baustein 03")
_s("N*", "N_stern", "Radiant der Bahn bis zur Übergabe", "rad", "L", "je_system", "Baustein 03")
_s("λ_K", "lambda_K", "Zerfallsrate des Gedächtniskerns der Kopfgleichung", "s⁻¹", "L", "je_system", "Baustein 05")
_s("Δ_B", "Delta_B", "Bindungsgrad eines Bereichs gegen ρ_krit, mit Rückläufen gezählt", "1", "L", "je_system", "Baustein 03")
_s("θ", "theta", "Ausdehnungsrate eines Bereichs im Raum", "s⁻¹", "L", "je_system", "Baustein 03")
_s("θ̃", "theta_tilde", "kleinste Ausdehnungsrate der Geschichte, max(min θ, 0)", "s⁻¹", "L", "je_system", "Baustein 03")

# ------------------------------------------------------------------ Ebene I: Rückkehr (03, 07)
_s("a₀", "a0", "Schwelle c²/(2Λ), Standardwert 1,1·10⁻¹⁰", "m/s²", "I", "seit_z10", A)
_s("Λ", "Lam", "Schließlänge der Rückkehr (innere Länge, keine Raumlänge)", "m", "I", "seit_z10", A)
_s("ℓ₀", "ell0", "Länge der Wirkung c²/a₀ = 2Λ", "m", "I", "seit_z10", "Baustein 08")
_s("τ_D", "tau_D", "Rückkehrzeit ℓ₀/c = c/a₀", "s", "I", "seit_z10", "Baustein 07")
_s("Γ", "Gamma", "Verlustrate a_N/c (Taktgefälle)", "s⁻¹", "LI", "je_system", "Baustein 07")
_s("θ_c", "theta_c", "kritische Ausdehnungsrate a₀/c = 1/τ_D", "s⁻¹", "I", "seit_z10", "Baustein 03")
_s("ϖ", "varpi", "Glättungstiefe der Rückkehr (Halbgerade)", "m", "I", "je_system", A)
_s("ϖ_R", "varpi_R", "Rückkehrtiefe, Punktmasse ϖ_R = r*", "m", "I", "je_system", "Baustein 07")
_s("D", "D_rueck", "Rückkehrdiffusion ϖ_R²/τ_D (Punktmasse GM/c)", "m²/s", "I", "je_system", "Baustein 07")
_s("u", "u_rueck", "Feld der Rückkehr in (x, ϖ)", "m/s²", "I", "je_system", "Baustein 07")
_s("p", "p_laplace", "Laplace-Variable des Kerns", "s⁻¹", "I", "je_system", "Baustein 07")
_s("s_p", "s_p", "Exponent der einen Rückkehr √((p + Γ + θ)·τ_D)", "1", "I", "je_system", "Baustein 07")
_s("û", "u_dach", "projizierte Rückkehr e^{−s_p}: q, 1 − w, K̂", "1", "I", "je_system", "Baustein 07")
_s("K̂", "K_dach", "Laplace-Bild des Kerns e^{−√((p+Γ)τ_D)}", "1", "I", "je_system", "Baustein 07")
_s("x_θ", "x_theta", "Rate mal Zeitkonstante eines Gedächtnisses, θ₊τ_D", "1", "I", "je_system", "Baustein 03")
_s("s", "s", "Bose-Variable √(a_N/a₀), Punktmasse r*/r", "1", "LI", "je_system", "Baustein 07")
_s("q", "q", "Nome e^{2πiτ}; |q| = e^{−s} ist der Rückkehranteil, als Regelkreis die Schleifenverstärkung", "1", "LIZG",
   "je_system", A)
_s("ν", "nu", "Verstärkung 1/(1 − q) = 1 + n_B(s)", "1", "LI", "je_system", A)
_s("χ", "chi", "Schalter s/(e^s − 1)", "1", "LI", "je_system", A)
_s("n_B", "n_B", "Bose-Besetzung 1/(e^s − 1)", "1", "LI", "je_system", "Baustein 07")
_s("Q̃", "Q_tilde", "Strukturwirkung −4∫_s^∞ x³n_B dx", "1", "LG", "je_system", "Baustein 07")
_s("Q_s", "Q_s", "Strukturteil Q̃ + 4π⁴/15", "1", "L", "je_system", "Baustein 07")
_s("d", "d_dim", "Dimension der Stufe der Bose-Leiter", "1", "Z", "fest", "Baustein 07")
_s("s_d", "s_d", "Maximum von s^d·n_B (Wien-Stelle)", "1", "Z", "fest", "Baustein 07")
_s("w", "w", "Gewicht der Schließungsregel 1 − e^{−√(θ̃/θ_c)}", "1", "G", "je_system", "Baustein 03")
_s("ξ", "xi", "Mode der Schließungsregel, (∂ξ)² = −1", "1", "G", "je_system", A)

# ------------------------------------------------------------------ Ebene Z: Zustandsraum (05, 06, 07)
_s("𝒫", "P_op", "Schließungsoperator (1/N)ΣT_{mL}", "1", "Z", "fest", A)
_s("Π⊥", "Pi_perp", "Projektor auf Flucht- und Instabilitätsrichtungen", "1", "Z", "je_system", A)
_s("L_eff", "L_eff", "Erzeuger der Bewegung mit Schließung, 𝒫L", "s⁻¹", "Z", "je_system", A)
_s("τ", "tau", "Modul des Torus; q = e^{2πiτ}", "1", "Z", "je_system", A)
_s("τ_M", "tau_M", "Modul des Bahntorus, Re = ω_rad/ω_azim, Im = Bindung", "1", "Z", "je_system", A, stamm="τ")
_s("τ_Z", "tau_Z", "Modul des flachen Torus der Zelle ℂ/(ℤ + τ_Zℤ)", "1", "Z", "fest", "Baustein 05", stamm="τ")
_s("Re τ_M", "re_tau_M", "Umlaufphase ω_rad/ω_azim", "1", "Z", "je_system", "Baustein 07")
_s("q_Δs", "q_ds", "Nome am Kern e^{2πiτ_M}, |q_Δs| = Δs/r", "1", "Z", "je_system", A, stamm="q")
_s("p/k", "p_zu_k", "rationales Verhältnis einer Resonanz", "1", "Z", "fest", "Baustein 07")
_s("r_p/k", "r_resonanz", "Radius der Resonanz Re τ_M = p/k (Punktmasse)", "m", "L", "je_system", "Baustein 07")
_s("Im τ_Z", "im_tau_Z", "Imaginärteil des Moduls τ_Z", "1", "Z", "fest", "Baustein 05")
_s("R_T/r_T", "R_zu_r_torus", "Radienverhältnis eines eingebetteten Torus", "1", "Z", "fest", "Baustein 05")
_s("b_S", "b_streifen", "Breite des offenen Streifens im Rohrwinkel", "rad", "Z", "fest", "Baustein 05")
_s("ν_T", "nu_torus", "Kehrwert des offenen Anteils der Torusfläche", "1", "Z", "fest", "Baustein 05")
_s("μ(X)", "mu_X", "Maß des ganzen Zustandsraums", "1", "Z", "fest", "Baustein 06")
_s("μ(A)", "mu_A", "Maß eines Bereichs A", "1", "Z", "je_system", "Baustein 06")
_s("⟨T⟩", "T_kac", "mittlere Rückkehrzeit in A (Kac), in Schritten", "1", "Z", "je_system", "Baustein 06")
_s("Σ1/⟨T_i⟩", "kac_budget", "Kac-Budget einer Zerlegung, 1 bei Schließung", "1", "Z", "fest", "Baustein 06")
_s("ε_L", "eps_L", "Leiterzahl ℓ_P/ℓ₀", "1", "Z", "seit_z10", "Baustein 08")
_s("M_min", "M_min", "untere Schnittmasse ħ/(c²τ_D)", "kg", "Z", "seit_z10", "Baustein 08")
_s("M_c", "M_c", "Schnittmasse r* = λ̄, ε_L^{1/3}·m_Pl", "kg", "Z", "seit_z10", "Baustein 08")
_s("M_U", "M_U", "obere Schnittmasse c³τ_D/G", "kg", "Z", "seit_z10", "Baustein 08")

# ------------------------------------------------------------------ Labor im Zustandsraum: Kupferdraht (06)
_s("ρ_el", "rho_el", "spezifischer Widerstand", "Ω·m", "L", "je_system", "Baustein 06")
_s("α_R", "alpha_R", "Temperaturbeiwert des Widerstands", "K⁻¹", "L", "je_system", "Baustein 06")
_s("E_Y", "E_young", "Elastizitätsmodul", "Pa", "L", "je_system", "Baustein 06")
_s("α_E", "alpha_E", "Temperaturbeiwert des Elastizitätsmoduls", "K⁻¹", "L", "je_system", "Baustein 06")
_s("α_th", "alpha_th", "Wärmeausdehnung", "K⁻¹", "L", "je_system", "Baustein 06")
_s("ρ_Cu", "rho_cu", "Dichte des Drahts", "kg/m³", "L", "je_system", "Baustein 06")
_s("c_p", "c_p", "spezifische Wärme", "J/(kg·K)", "L", "je_system", "Baustein 06")
_s("k_th", "k_th", "Wärmeleitfähigkeit", "W/(m·K)", "L", "je_system", "Baustein 06")
_s("a_D", "a_draht", "Drahtradius", "m", "L", "je_system", "Baustein 06")
_s("ℓ_D", "ell_draht", "Drahtlänge", "m", "L", "je_system", "Baustein 06")
_s("B", "B_feld", "Magnetfeld", "T", "L", "je_system", "Baustein 06")
_s("Q_m", "Q_m", "mechanische Güte", "1", "L", "je_system", "Baustein 06")
_s("σ_0", "sigma_0", "mechanische Vorspannung", "Pa", "L", "je_system", "Baustein 06")
_s("R_el,0", "R0_el", "elektrischer Widerstand bei ΔT = 0", "Ω", "L", "je_system", "Baustein 06")
_s("G_th", "G_th", "Wärmeleitwert", "W/K", "L", "je_system", "Baustein 06")
_s("C_th", "C_th", "Wärmekapazität", "J/K", "L", "je_system", "Baustein 06")
_s("m_eff", "m_eff", "effektive Masse der Mode", "kg", "L", "je_system", "Baustein 06")
_s("Γ_B", "Gamma_B", "Kopplung Bℓ", "T·m", "L", "je_system", "Baustein 06")
_s("k_st", "k_steif", "Steifigkeit der Mode k(ΔT)", "N/m", "L", "je_system", "Baustein 06")
_s("ω₀", "omega0", "Eigenfrequenz der Mode", "s⁻¹", "L", "je_system", "Baustein 06")
_s("γ_m", "c_m", "mechanische Dämpfung m·ω₀/Q_m", "kg/s", "L", "je_system", "Baustein 06")

# ------------------------------------------------------------------ Ebene G: das Ganze (03, 04)
_s("L", "L_zelle", "Zelle, räumliche Schleife des Torus (mitbewegt; fest wie a₀ und ρ_m)", "m", "G", "seit_z10", A)
_s("L_min", "L_min", "kürzeste Schleife ohne Kreispaare, 0,985·d_LSS", "m", "G", "epoche", MESS)
_s("Q", "Q", "kosmisches Gewicht q/(1 − q) = Ω_c/Ω_b", "1", "G", "fest", A)
_s("q_kosm", "q_kosm", "Rückkehranteil des Ganzen Q/(1 + Q)", "1", "G", "fest", "Baustein 03", stamm="q")
_s("ν_kosm", "nu_kosm", "mittlere Zahl der Rückläufe im Kosmos 1 + Q (Kac)", "1", "G", "fest", "Baustein 03")
_s("σ_Q", "sigma_Q", "Messfehler von Q aus σ(ω_b), σ(ω_c)", "1", "G", "fest", MESS)
_s("Q_min", "Q_min", "untere Grenze des Fensters von Q", "1", "G", "fest", "Baustein 03")
_s("Q_max", "Q_max", "obere Grenze des Fensters von Q (hängt an L_min, also am Alter des Beobachters)", "1", "G", "epoche", "Baustein 03")
_s("H", "H", "Atemrate ȧ/a", "s⁻¹", "G", "epoche", "Baustein 03")
_s("H₀", "H0", "Atemrate heute (im Kern in km/s/Mpc)", "s⁻¹", "G", "epoche", "Messwert Planck 2018 [G]")
_s("h_H", "h_", "H₀/(100 km/s/Mpc)", "1", "G", "epoche", MESS)
_s("z", "z", "Rotverschiebung", "1", "G", "epoche", "Baustein 03")
_s("E(z)", "E_z", "H(z)/H₀ im Hintergrund", "1", "G", "epoche", "Baustein 03")
_s("d_C", "d_C", "mitbewegter Abstand", "m", "G", "epoche", "Baustein 03")
_s("d_LSS", "d_LSS", "Durchmesser der letzten Streufläche (mitbewegt)", "m", "G", "epoche", "Baustein 03")
_s("Θ", "Theta", "2πa₀/(cH)", "1", "G", "epoche", "Baustein 03")
_s("Λ_s", "Lambda_s", "Vakuumteil π⁴/(15Λ²)", "m⁻²", "G", "seit_z10", "Baustein 03")
_s("Ω_s", "Omega_s", "Vakuumteil als Anteil (π²/45)Θ²", "1", "G", "epoche", "Baustein 03")
_s("K_s", "K_S", "Vorfaktor π²/45 des Vakuumteils bei Θ = 1", "1", "G", "fest", "Baustein 03")
_s("H_s", "H_s", "Atemrate des Vakuumteils allein, H_s² = Λ_s c²/3", "s⁻¹", "G", "seit_z10", "Baustein 03")
_s("ρ_m", "rho_m", "mittlere Dichte der Materie, mitbewegt, ν_kosm·ρ_b", "kg/m³", "G", "fest", "Baustein 03")
_s("ρ_b", "rho_b", "mittlere Dichte der Baryonen, mitbewegt", "kg/m³", "G", "fest", "Baustein 03")
_s("Ω_b", "Omega_b", "Anteil der Baryonen heute", "1", "G", "epoche", "Messwert Planck 2018 [G]")
_s("Ω_m", "Omega_m", "Anteil der Materie heute, Ω_b(1 + Q)", "1", "G", "epoche", "Baustein 03")
_s("Ω_r", "Omega_r", "Anteil der Strahlung heute", "1", "G", "epoche", MESS)
_s("Ω_K", "Omega_K", "Krümmungsanteil heute", "1", "G", "epoche", MESS)
_s("Ω_DE", "Omega_DE", "Anteil der Dunklen Energie heute", "1", "G", "epoche", MESS)
_s("H_DE", "H_DE", "Atemrate der Dunklen Energie", "s⁻¹", "G", "epoche", "Baustein 03")
_s("ω_b", "omega_b", "Ω_b h²", "1", "G", "fest", "Messwert Planck 2018 [G]")
_s("ω_c", "omega_c", "Ω_c h²", "1", "G", "fest", "Messwert Planck 2018 [G]")
_s("a₀(z)/a₀", "a0_frueh", "Schwelle früh gegen heute", "1", "G", "epoche", "Baustein 03")
_s("f_N", "f_anteil", "fester Anteil am Gesamtbudget im Strahlungszeitalter", "1", "G", "epoche", "Baustein 03")
_s("N_eff", "N_eff", "effektive Zahl der Neutrinoarten", "1", "G", "fest", MESS)
_s("ΔN_eff", "dN_eff", "Zusatz zu N_eff", "1", "G", "epoche", "Baustein 03")
_s("S_H", "S_H", "Horizontentropie πc²/(ℓ_P²H²) in Einheiten von k_B", "1", "G", "epoche", "Baustein 03")
_s("Ṡ_i", "S_i_dot", "Entropieproduktion je Horizont in k_B/s", "s⁻¹", "G", "epoche", "Baustein 04")
_s("Λ̇/Λ", "Lam_rate", "Wachstumsrate der Schließlänge", "s⁻¹", "G", "epoche", "Baustein 03")
_s("w_Z", "w_Z", "Zustandsgleichung p/(ρc²) (w_Z0, w_Za, w_Zs)", "1", "G", "epoche", "Baustein 03")
_s("Δw_Z", "w_Z_abstand", "Abstand der Zustandsgleichung von −1: 1 + w_Z", "1", "G", "epoche", "Baustein 03")
_s("σ_B", "sigma_bilanz", "Faltungsbilanz ΣΔs − ΔL = 0", "m", "G", "fest", A)

# ------------------------------------------------------------------ Brücke L–G: Naht und Haufen (03)
_s("R_Q", "R_Q", "Radius, an dem das örtliche q das kosmische erreicht", "m", "LG", "je_system", "Baustein 03")
_s("M_b,Q", "M_bQ", "Baryonenmasse, bei der R_Q die Grenze R₂₀₀ erreicht", "kg", "LG", "epoche", "Baustein 03")
_s("Δ_ta", "Delta_ta", "Dichte gegen ρ_m einer heute umkehrenden Schale (1 + Δ_ta)", "1", "LG", "epoche", "Baustein 03")
_s("R_ta", "R_ta", "Umkehrradius heute", "m", "LG", "epoche", "Baustein 03")
_s("M_ges", "M_ges", "Maß des Ganzen (1 + Q)·M_b", "kg", "LG", "je_system", "Baustein 03")
from . import symbole_naht_kosmos  # noqa: E402,F401  (Symbole von kosmos und naht_haufen)
from . import symbole_teilchen_licht  # noqa: E402,F401  (Symbole der Gruppe Teilchen, Sterne, Licht)
from . import symbole_rueckkehr_zustand  # noqa: E402,F401  (Symbole der Rückkehr, Faltung und des Zustandsraums)
from . import symbole_galaxie_daten  # noqa: E402,F401  (Symbole der Gruppe Galaxie und Daten)
