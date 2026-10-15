"""Rechnungen zum Wardenclyffe-Turm (Tesla, Shoreham 1901 bis 1917) mit den Formen des Modells: Schließung der Erdschale je
Umlauf und in Kugelfunktionen, Kuppel im Coulombfeld (Phase statt Faltung), Erde als Gegenplatte (ΣQ_el = 0), Turm als
Schwinger und sein Einrasten an der Erdschale, die verstärkende Rückkehr Σqⁿ, was die Erdschale für Energie leisten kann,
die Regel am Turm. Lauf: python simulation/turm.py gibt alle Zahlen des Berichts."""
import functools
import math

import kern as K

# Teslas Baumaße (Eingabe, nach den Bauplänen): Turm 187 ft, Kuppel 68 ft Durchmesser, Schacht 120 ft, Eisenrohre 300 ft
TURM = dict(hoehe=57.0, kuppel=20.7, schacht=36.6, rohre=91.4)
E_LUFT = 3.0e6                          # V/m, Durchschlagsfeld trockener Luft auf Meereshöhe (Eingabe)
P_EIN = 1.0e6                           # W, Einspeisung je Megawatt (Bezugsgröße der Rechnung)
# Auslegung des Aufbaus (Eingaben, gewählt): Einspeisung des Kraftwerks, Wendel aus Kupferrohr im Turm, Erdung, Primärkreis
AUSLEGUNG = dict(P=2.0e5, wendel=6.0, rohr=0.05, gang=2.0, boden=100.0, stab_radius=0.05, staebe=16, guete_primaer=50.0,
                 spannung_primaer=5.0e4)
SIGMA_CU = 5.96e7                       # S/m, Leitfähigkeit von Kupfer bei 20 °C (Eingabe)


def erdschale(n_max=5):
    """Die Erdschale schließt je Umlauf, 2πR = n·λ: f_n = n·c/(2πR); als Fläche einer Kugel schließen ihre stehenden Wellen
    in Kugelfunktionen, f_n = √(n(n + 1))·c/(2πR) mit 2n + 1 Moden je Grad. Gibt (Umlaufzeit [s], [f Umlauf], [f Kugel])."""
    f1 = K.c / (2 * math.pi * K.R_ERDE)
    return 1 / f1, [n * f1 for n in range(1, n_max + 1)], [math.sqrt(n * (n + 1)) * f1 for n in range(1, n_max + 1)]


def kuppel(aus=None):
    """Die Kuppel als Kugel vom Radius r im Coulombfeld (1/r², ohne Schwelle): C = 4πε₀r, größte Spannung bis zum
    Durchschlag V = E_Luft·r, Ladung Q = C·V, Energie W = ½CV²; aus: geänderte Baumaße. Gibt dict."""
    r = {**TURM, **(aus or {})}["kuppel"] / 2
    C = 4 * math.pi * K.eps0 * r
    V = E_LUFT * r
    return dict(r=r, C=C, V=V, Q=C * V, W=C * V * V / 2)


def erde_gegenplatte(Q=None):
    """Die Erde als Gegenplatte: C_E = 4πε₀R; die Ladung Q der Kuppel verschiebt ihr Potential um Q/C_E (ΣQ_el = 0: was
    die Kuppel abgibt, nimmt die Erde auf und gibt es zurück). Gibt (C_E [F], ΔV [V])."""
    Q = kuppel()["Q"] if Q is None else Q
    C_E = 4 * math.pi * K.eps0 * K.R_ERDE
    return C_E, Q / C_E


def nagaoka(r, l):
    """Nagaoka-Zahl einer Stromlage vom Radius r und der Länge l, exakt aus den elliptischen Integralen:
    K_N = 4/(3πk′)·[(k′²/k²)(K − E) + E − k] mit k² = 4r²/(4r² + l²); lange Spule 1, kurze 0."""
    from scipy.special import ellipe, ellipk
    m = 4 * r * r / (4 * r * r + l * l)
    k, kp = math.sqrt(m), math.sqrt(1 - m)
    return 4 / (3 * math.pi * kp) * ((1 - m) / m * (ellipk(m) - ellipe(m)) + ellipe(m) - k)


def wendel(N=None, aus=None):
    """Die Wendel im Turm: N Windungen vom Durchmesser D über die Turmhöhe h (höchstens h/(Gang·Rohr)). Induktivität der
    einlagigen Stromlage L = μ₀πr²N²/h·K_N (Nagaoka, exakt), als zweiter Weg μ₀πr²N²/(h + 0,9r) (Wheeler), Eigenkapazität
    C_W = D·(0,1126·h/D + 0,08 + 0,27/√(h/D)) pF mit D in cm (Medhurst), Drahtlänge NπD. Mit der Kuppel schwingt sie bei f = 1/(2π√(L(C_K + C_W))); als Viertelwelle des Drahts
    schließt die stehende Welle bei c/(4·NπD); aus: geänderte Auslegungswerte. Gibt dict."""
    A = {**TURM, **AUSLEGUNG, **(aus or {})}
    h, D = A["hoehe"], A["wendel"]
    N = int(h / (A["gang"] * A["rohr"])) if N is None else N
    L = K.mu0 * math.pi * (D / 2)**2 * N * N / h * nagaoka(D / 2, h)
    L_wheeler = K.mu0 * math.pi * (D / 2)**2 * N * N / (h + 0.9 * D / 2)
    C_W = D * 100 * (0.1126 * h / D + 0.08 + 0.27 / math.sqrt(h / D)) * 1e-12
    C = kuppel(aus)["C"] + C_W
    lw = N * math.pi * D
    return dict(N=N, L=L, L_wheeler=L_wheeler, C_W=C_W, C=C, draht=lw, f=1 / (2 * math.pi * math.sqrt(L * C)),
                f_viertel=K.c / (4 * lw))


def verluste(N=None, aus=None):
    """Verluste der Wendel bei ihrer Frequenz, alle im transversalen Kanal oder als Abgabe: Joule-Wärme im Kupferrohr bis zur
    Hauttiefe δ = √(2/(ωμ₀σ)), R_D = ℓ/(σ·πd·δ) (ohne Nähewirkung der Windungen), Abstrahlung des oben belasteten Turms
    R_S = 160π²(h/λ)², Erdung aus 16 Rohren R_E = ρ/(2πℓ)·(ln(4ℓ/a) − 1)/16 (ohne gegenseitige Abschirmung, untere Grenze).
    Güte Q_s = ωL/R; aus: geänderte Auslegungswerte, mit supra=True leitet die Wendel ohne Verlust (vollständiges
    Einrasten, der transversale Kanal geschlossen). Gibt dict."""
    A, w_ = {**TURM, **AUSLEGUNG, **(aus or {})}, wendel(N, aus)
    om = 2 * math.pi * w_["f"]
    delta = math.sqrt(2 / (om * K.mu0 * SIGMA_CU))
    R_D = 0.0 if A.get("supra") else w_["draht"] / (SIGMA_CU * math.pi * A["rohr"] * delta)
    R_S = 160 * math.pi**2 * (A["hoehe"] * w_["f"] / K.c)**2
    lr = A["rohre"]
    R_E = A["boden"] / (2 * math.pi * lr) * (math.log(4 * lr / A["stab_radius"]) - 1) / A["staebe"]
    R = R_D + R_S + R_E
    return dict(delta=delta, R_D=R_D, R_S=R_S, R_E=R_E, R=R, Q=om * w_["L"] / R)


def betrieb(N=None, P=None, aus=None):
    """Betriebspunkt bei der Einspeisung P: Strom I = √(P/R) (Effektivwert), Scheitelspannung der Kuppel auf zwei Wegen,
    über den Strom V = √2·I/(ωC) und über die gespeicherte Energie W = Q_s·P/ω, V = √(2W/C); Feld an der Kuppel V/r gegen
    den Durchschlag, Verstärkung ν = Q_s/π, Aufteilung der Leistung; aus: geänderte Auslegungswerte. Gibt dict."""
    P = AUSLEGUNG["P"] if P is None else P
    w_, v = wendel(N, aus), verluste(N, aus)
    om = 2 * math.pi * w_["f"]
    I = math.sqrt(P / v["R"])
    V = math.sqrt(2) * I / (om * w_["C"])
    W = v["Q"] * P / om
    V2 = math.sqrt(2 * W / w_["C"])
    E = V / kuppel(aus)["r"]
    return dict(I=I, V=V, V_energie=V2, W=W, E=E, E_anteil=E / E_LUFT, nu=v["Q"] / math.pi,
                P_D=I * I * v["R_D"], P_S=I * I * v["R_S"], P_E=I * I * v["R_E"])


def primaerkreis(N=None, P=None):
    """Primärkreis, auf die Wendel abgestimmt (Δω = 0, Einrasten): kritische Kopplung k = 1/√(Q_p·Q_s), zwischen zu schwach
    (die Schwingung zerfällt vor der Übergabe) und zu stark (die Moden spalten auf), Übergabezeit 1/(2kf); umlaufende
    Energie W_p = Q_p·P/ω, Kapazität C_p = 2W_p/V_p², Induktivität L_p = 1/(ω²C_p). Gibt dict."""
    A = AUSLEGUNG
    P = A["P"] if P is None else P
    f = wendel(N)["f"]
    om = 2 * math.pi * f
    k = 1 / math.sqrt(A["guete_primaer"] * verluste(N)["Q"])
    W_p = A["guete_primaer"] * P / om
    C_p = 2 * W_p / A["spannung_primaer"]**2
    return dict(k=k, t=1 / (2 * k * f), W_p=W_p, C_p=C_p, L_p=1 / (om * om * C_p))


def signal(d, N=None, aus=None):
    """Feld der Bodenwelle [V/m, effektiv] im Abstand d [m] für den oben belasteten kurzen Turm über leitender Ebene,
    E = √(90·P_S)/d, mit der abgestrahlten Leistung des Betriebspunkts."""
    return math.sqrt(90 * betrieb(N, aus=aus)["P_S"]) / d


def grenzleistung(aus=None):
    """Die Einspeisung, bei der die Kuppel den Durchschlag der Luft erreicht, P = P_A·(V_max/V)² (V ∝ √P), und die dann in
    der Schwingung gespeicherte Energie ½·C·V_max². Gibt (P [W], W [J])."""
    V_max = kuppel(aus)["V"]
    return AUSLEGUNG["P"] * (V_max / betrieb(aus=aus)["V"])**2, wendel(aus=aus)["C"] * V_max**2 / 2


def viertelwelle_strahlung(P=None, d=1e6):
    """Der Turm ohne Wendel als Viertelwelle über der Erdung: Strahlungswiderstand 36,5 Ω (die Hälfte des Halbwellendipols),
    Anteil der abgestrahlten Leistung gegen die Erdung und Feld der Bodenwelle im Abstand d. Gibt (Anteil, E [V/m])."""
    P = AUSLEGUNG["P"] if P is None else P
    eta = 36.5 / (36.5 + verluste()["R_E"])
    return eta, math.sqrt(90 * eta * P) / d


def kopplung(d, aus=None, ringe=24, punkte=96):
    """Kopplung zweier gleicher Wendeln Seite an Seite im Achsabstand d [m], aus dem Neumann-Integral
    M = μ₀/(4π)·∮∮ dl₁·dl₂/|r₁ − r₂| über die Stromlagen (je Wendel ringe Ringe mit N/ringe Windungen, punkte Stützstellen je
    Ring), k = |M|/L. Fern gilt der Dipol Seite an Seite, k → r²(h + 0,9r)/(4d³), als zweiter Weg. Gibt (k, k_Dipol)."""
    import numpy as np
    A = {**TURM, **AUSLEGUNG, **(aus or {})}
    w_ = wendel(aus=aus)
    r, h, N = A["wendel"] / 2, A["hoehe"], w_["N"]
    z = (np.arange(ringe) + 0.5) * h / ringe
    ph = np.arange(punkte) * 2 * math.pi / punkte
    x1, y1 = r * np.cos(ph), r * np.sin(ph)
    tx, ty = -np.sin(ph) * r * 2 * math.pi / punkte, np.cos(ph) * r * 2 * math.pi / punkte
    skalar = tx[:, None] * tx[None, :] + ty[:, None] * ty[None, :]
    abst2 = (x1[:, None] - d - x1[None, :])**2 + (y1[:, None] - y1[None, :])**2
    Mz = [K.mu0 / (4 * math.pi) * np.sum(skalar / np.sqrt(abst2 + (z[j] - z[0])**2)) for j in range(ringe)]
    M = sum(Mz[abs(i - j)] for i in range(ringe) for j in range(ringe)) * (N / ringe)**2
    return abs(M) / w_["L"], r * r * (h + 0.9 * r) / (4 * d**3)


def empfaenger(d, aus=None):
    """Ein gleicher Turm als Empfänger im Abstand d [m], neben dem Sender, beide auf dieselbe Frequenz eingerastet: mit
    U = k·Q_s (k aus kopplung) ist der beste Wirkungsgrad η = U²/(1 + √(1 + U²))², die Grenze des Kopplungsfensters liegt bei
    U = 1. Gibt (k, U, η)."""
    k = kopplung(d, aus)[0]
    U = k * verluste(aus=aus)["Q"]
    return k, U, U * U / (1 + math.sqrt(1 + U * U))**2


def reichweite(eta=0.5, aus=None):
    """Abstand [m], bis zu dem ein gleicher, eingerasteter Turm den Wirkungsgrad η erreicht (η aus empfaenger, Nullstelle
    zwischen dem Durchmesser der Wendel und 100 km)."""
    from scipy.optimize import brentq
    A = {**AUSLEGUNG, **(aus or {})}
    return brentq(lambda d: empfaenger(d, aus)[2] - eta, 1.5 * A["wendel"], 1e5, xtol=0.01)


def anschwingen(aus=None):
    """Zeit, in der die Amplitude der Wendel auf ihren Endwert zuläuft, 2Q_s/ω = Q_s/(πf), und die Breite der Resonanz f/Q_s,
    die eine Phasenregelung (Einrasten des Antriebs an der Wendel) halten muss. Gibt (Zeit [s], Breite [Hz])."""
    f, Q = wendel(aus=aus)["f"], verluste(aus=aus)["Q"]
    return Q / (math.pi * f), f / Q


def kopplung_fluss(d, aus=None, ringe=24, radial=24, winkel=48):
    """Zweiter Weg zur Kopplung: der Fluss des Feldes jedes Rings der einen Wendel durch die Kreisflächen der Ringe der
    anderen, das Feld geschlossen aus den elliptischen Integralen, B_z = μ₀I/(2π)/√((a + ρ)² + z²)·[K(m) + (a² − ρ² − z²)/
    ((a − ρ)² + z²)·E(m)] mit m = 4aρ/((a + ρ)² + z²). Gibt k = |M|/L."""
    import numpy as np
    from scipy.special import ellipe, ellipk
    A = {**TURM, **AUSLEGUNG, **(aus or {})}
    w_ = wendel(aus=aus)
    a, h, N = A["wendel"] / 2, A["hoehe"], w_["N"]
    xs, gs = np.polynomial.legendre.leggauss(radial)
    s_, th = a * (xs + 1) / 2, np.arange(winkel) * 2 * math.pi / winkel
    rho = np.hypot(d + s_[:, None] * np.cos(th)[None, :], s_[:, None] * np.sin(th)[None, :])
    dA = (gs * a / 2 * s_)[:, None] * (2 * math.pi / winkel)
    Mz = []
    for j in range(ringe):
        z = j * h / ringe
        m = 4 * a * rho / ((a + rho)**2 + z * z)
        Bz = K.mu0 / (2 * math.pi) / np.sqrt((a + rho)**2 + z * z) * (
            ellipk(m) + (a * a - rho**2 - z * z) / ((a - rho)**2 + z * z) * ellipe(m))
        Mz.append(np.sum(Bz * dA))
    M = sum(Mz[abs(i - j)] for i in range(ringe) for j in range(ringe)) * (N / ringe)**2
    return abs(M) / w_["L"]


def bauform(x, aus=None, groesse=1.0):
    """Eine Wendel aus derselben Kupfermenge wie in Teslas Maßen (Rohr und Ganghöhe gleich, also D·h fest) mit dem
    Verhältnis x = D/h, die Kuppel so weit wie die Wendel, damit sie deren obere Windungen schirmt; groesse: alle Längen
    außer dem Rohr mal diesem Faktor. Gibt die Auslegung."""
    h = groesse * math.sqrt(TURM["hoehe"] * AUSLEGUNG["wendel"] / x)
    return {**(aus or {}), "wendel": x * h, "hoehe": h, "kuppel": x * h}


def form_guete(aus=None):
    """Das Verhältnis x = D/h, bei dem die Güte der Wendel aus derselben Kupfermenge am größten ist. Gibt (x, Q_s)."""
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda x: -verluste(aus=bauform(x, aus))["Q"], bounds=(0.2, 4.0), method="bounded",
                        options=dict(xatol=1e-3))
    return r.x, -r.fun


SUPRA = dict(supra=True, staebe=64)     # supraleitende Wendel, 64 Erdrohre, Antrieb mit Phasenregelung


@functools.lru_cache(maxsize=None)
def formen():
    """Die beiden Bauformen aus derselben Leitermenge mit der größten Güte: Kupfer (16 Erdrohre) und supraleitend (SUPRA).
    Gibt (Kupfer, supraleitend) als Auslegungen."""
    return bauform(form_guete()[0]), bauform(form_guete(SUPRA)[0], SUPRA)


def guete_geschlossen(aus=None):
    """Zweiter Weg zur Güte des Kupfers: Mit N = h/(g·d) wird ωL/R_D = π·K·D/(2gδ), unabhängig von der Zahl der Windungen;
    K nach Wheeler 1/(1 + 0,9r/h), δ bei der Frequenz der Wendel. Gibt Q_D."""
    A = {**TURM, **AUSLEGUNG, **(aus or {})}
    D, h = A["wendel"], A["hoehe"]
    delta = math.sqrt(2 / (2 * math.pi * wendel(aus=aus)["f"] * K.mu0 * SIGMA_CU))
    return math.pi * D / (1 + 0.9 * D / 2 / h) / (2 * A["gang"] * delta)


def turm_schwinger():
    """Der senkrechte Leiter als Viertelwelle, f = c/(4ℓ): Turm und Schacht, und mit den Eisenrohren. Gibt (f, f mit Rohren)."""
    l1 = TURM["hoehe"] + TURM["schacht"]
    return K.c / (4 * l1), K.c / (4 * (l1 + TURM["rohre"]))


def spule_fuer(f, C=None):
    """Induktivität [H], die mit der Kuppel bei der Frequenz f [Hz] schwingt: L = 1/((2πf)²·C)."""
    C = kuppel()["C"] if C is None else C
    return 1 / ((2 * math.pi * f)**2 * C)


def einrasten(f_turm=None, f_erde=None):
    """Zwei Zyklen rasten ein, wenn |Δω| ≤ 2k|cos φ|: die Kopplung, die der Turm (Wendel mit Kuppel) an der tiefsten
    Kugelmode der Erdschale bräuchte, je Kreisfrequenz des Turms. Gibt (f_Turm/f_Erde, k/ω_Turm mindestens)."""
    f_turm = wendel()["f"] if f_turm is None else f_turm
    f_erde = erdschale()[2][0] if f_erde is None else f_erde
    return f_turm / f_erde, (f_turm - f_erde) / (2 * f_turm)


def rueckkehr(guete):
    """Die verstärkende Rückkehr eines Schwingers der Güte Q_s, angestoßen je Umlauf: Je Umlauf bleibt q = e^{−π/Q_s} der
    Amplitude, die Summe der Rückläufe ist ν = Σqⁿ = 1/(1 − q) ≈ Q_s/π. Gibt (q, ν)."""
    q = math.exp(-math.pi / guete)
    return q, 1 / (1 - q)


def schale_fuer_energie(fluss=1.0, P=P_EIN, guete=None):
    """Was die Erdschale für Energie leisten kann: Die Rückkehr hebt die umlaufende Amplitude um ν = 1/(1 − q), die
    gespeicherte Energie wird W = P·T·ν; abgeben kann die Schale im Gleichgewicht je Umlauf nur, was der Turm einspeist
    (fester Haushalt der Wiederkehr). Gibt (abgebbare Leistung je Fläche der Erde bei P [W/m²], Einspeisung für den Fluss
    überall [W], gespeicherte Energie bei der Güte [J])."""
    T = erdschale()[0]
    A = 4 * math.pi * K.R_ERDE**2
    W = P * T * rueckkehr(guete)[1] if guete else None
    return P / A, fluss * A, W


def regel_am_turm():
    """Die Regel am Ort des Turms: s = √(g/a₀) mit g an der Erdoberfläche, der Rückkehranteil q = e^{−s} (als lg).
    Die Windung trägt keine Regel (Phase statt Faltung): Für das Feld des Turms gibt es keine Schwelle."""
    g = K.G * K.M_earth / K.R_ERDE**2
    s = math.sqrt(g / K.a0)
    return g, s, -s / math.log(10)


if __name__ == "__main__":
    T, fu, fk = erdschale()
    print("Erdschale: Umlaufzeit [s]", T, "f je Umlauf", fu, "f Kugel", fk)
    print("Kuppel", kuppel())
    print("Erde als Gegenplatte C_E, ΔV", erde_gegenplatte())
    print("Wendel", wendel())
    print("Verluste", verluste())
    print("Betrieb", betrieb())
    print("Primärkreis", primaerkreis())
    print("Signal bei 100 und 1000 km [V/m]", signal(1e5), signal(1e6))
    for N in (100, 200, 400):
        print("Vergleich N", N, wendel(N)["f"], verluste(N)["Q"], betrieb(N)["V"], betrieb(N)["P_S"])
    print("Turm als Viertelwelle, mit Rohren [Hz]", turm_schwinger())
    print("Spule für die tiefste Kugelmode [H]", spule_fuer(fk[0]))
    print("Einrasten f_Turm/f_Erde, k/ω", einrasten())
    print("Rückkehr bei Güte 100", rueckkehr(100.0))
    print("Schale: je Fläche bei 200 kW, Einspeisung für 1 W/m²", schale_fuer_energie(P=AUSLEGUNG["P"]))
    print("Regel am Turm g, s, lg q", regel_am_turm())
    F, S = formen()
    for name, a in (("Teslas Maße", None), ("Kupfer D = h", F), ("supraleitend", S)):
        w_, v, b = wendel(aus=a), verluste(aus=a), betrieb(aus=a)
        print(name, a, "f", w_["f"], "N", w_["N"], "Q", v["Q"], "R", v["R_D"], v["R_S"], v["R_E"], "V", b["V"], "W", b["W"],
              "Grenze", grenzleistung(a), "Reichweite 50 %", reichweite(aus=a), "Anschwingen", anschwingen(a))
    print("Viertelwelle: Anteil, Signal 1000 km", viertelwelle_strahlung())
