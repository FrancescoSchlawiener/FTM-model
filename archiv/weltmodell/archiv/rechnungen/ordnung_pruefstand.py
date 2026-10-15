# -*- coding: utf-8 -*-
"""
Ordnung der Prüfzahlen: Wovon hängt jede Zahl des Prüfstands ab?
Jede Grundgröße des Kerns wird um ±1 % verschoben, der ganze Prüfstand neu gerechnet; daraus der Exponent
d ln(Zahl)/d ln(Grundgröße). Ergebnis: Exponentenvektor je Zahl. Gleiche Vektoren = gleiche Struktur,
Vektor null = reine Zahl, lineare Abhängigkeiten = exakte Verknüpfungen zwischen Zahlen.
Massen von Himmelskörpern skalieren gemeinsam (M), Radien einzeln; Zahlen in kpc/Mpc sind feste Daten.
"""
import io, sys, re, math, types, itertools
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = r"C:\DEV\FTM\weltmodell\rechnungen"
KSRC = open(HERE + r"\ftm_kern.py", encoding="utf-8").read()
PSRC = open(HERE + r"\ftm_pruefstand.py", encoding="utf-8").read()
PSRC = PSRC[:PSRC.index("# ---- Ausgabe")].replace("sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding=\"utf-8\")", "")

# Grundgrößen: (Name, Zeile im Kern, Ersatz mit Faktor F['name'])
INPUTS = {
    "G":     ("G = 6.67430e-11", "G = 6.67430e-11 * F['G']"),
    "c":     ("c = 299792458.0", "c = 299792458.0 * F['c']"),
    "ħ":     ("hbar = 6.62607015e-34 / (2 * math.pi)", "hbar = 6.62607015e-34 / (2 * math.pi) * F['ħ']"),
    "k":     ("kB = 1.380649e-23", "kB = 1.380649e-23 * F['k']"),
    "e":     ("e_el = 1.602176634e-19", "e_el = 1.602176634e-19 * F['e']"),
    "m_p":   ("m_p = 1.67262192369e-27", "m_p = 1.67262192369e-27 * F['m_p']"),
    "m_e":   ("m_e = 9.1093837015e-31", "m_e = 9.1093837015e-31 * F['m_e']"),
    "α":     ("alpha = 1 / 137.035999084", "alpha = 1 / 137.035999084 * F['α']"),
    "R_p":   ("R_p = 0.8414e-15", "R_p = 0.8414e-15 * F['R_p']"),
    "M":     ("GM_sun = 1.32712440018e20  # IAU\nM_sun = GM_sun / G", "M_sun = 1.32712440018e20 / 6.67430e-11 * F['M']\nGM_sun = G * M_sun"),
    "R☉":    ("R_sun = 6.957e8", "R_sun = 6.957e8 * F['R☉']"),
    "R⊕":    ("R_earth = 6.371e6", "R_earth = 6.371e6 * F['R⊕']"),
    "AU":    ("AU = 1.495978707e11", "AU = 1.495978707e11 * F['AU']"),
    "a₀":    ("A0_LIT = 1.2e-10", "A0_LIT = 1.2e-10 * F['a₀']"),
    "H₀":    ("H0_KMS = 67.4 ", "H0_KMS = 67.4 * F['H₀'] "),
    "Ω_bh²": ("OMEGA_B_H2 = 0.02237", "OMEGA_B_H2 = 0.02237 * F['Ω_bh²']"),
    "Ω_ch²": ("OMEGA_C_H2 = 0.1200", "OMEGA_C_H2 = 0.1200 * F['Ω_ch²']"),
    "T_CMB": ("T_CMB = 2.7255", "T_CMB = 2.7255 * F['T_CMB']"),
}
EXTRA = [("r_e = 2.8179403262e-15", "r_e = alpha * hbar / (m_e * c)"),
         ("a_B = 5.29177210903e-11", "a_B = hbar / (alpha * m_e * c)"),
         ("GM_earth = 3.986004418e14\nM_earth = GM_earth / G", "M_earth = 3.986004418e14 / 6.67430e-11 * F['M']\nGM_earth = G * M_earth")]

def kern_mit(F):
    src = KSRC
    for k, (old, new) in INPUTS.items():
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    for old, new in EXTRA:
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    mod = types.ModuleType("ftm_kern"); mod.__dict__["F"] = F
    exec(compile(src, "ftm_kern", "exec"), mod.__dict__)
    sys.modules["ftm_kern"] = mod
    ns = {"__name__": "pruef"}
    exec(compile(PSRC, "pruef", "exec"), ns)
    return [(b, t, float(k)) for b, t, m, k, tol in ns["checks"]]

base = kern_mit({k: 1.0 for k in INPUTS})
names = list(INPUTS)
E = np.zeros((len(base), len(names)))
h = 0.01
for j, nm in enumerate(names):
    Fp = {k: 1.0 for k in INPUTS}; Fm = dict(Fp); Fp[nm] = 1 + h; Fm[nm] = 1 / (1 + h)
    vp, vm = kern_mit(Fp), kern_mit(Fm)
    for i in range(len(base)):
        a, b_ = abs(vp[i][2]), abs(vm[i][2])
        E[i, j] = (math.log(a) - math.log(b_)) / (2 * math.log(1 + h)) if a > 1e-95 and b_ > 1e-95 else 0.0

def rund(x):   # auf Zwölftel runden, wenn nah
    r = round(x * 12) / 12
    return r if abs(x - r) < 0.01 else x

def fmt(v):
    parts = []
    for nm, x in zip(names, v):
        if abs(x) > 0.004:
            r = rund(x)
            from fractions import Fraction
            fr = Fraction(r).limit_denominator(12) if abs(r - round(r * 12) / 12) < 1e-9 else None
            parts.append(f"{nm}^{fr}" if fr is not None and fr != 1 else (nm if fr == 1 else f"{nm}^{x:.3f}"))
    return " ".join(parts) if parts else "rein"

print(f"{len(base)} Zahlen, {len(names)} Grundgrößen: " + ", ".join(names))
M = np.array([[rund(x) for x in row] for row in E])
print(f"Rang der Exponentenmatrix: {np.linalg.matrix_rank(M, tol=1e-3)}")

# 1 reine Zahlen
print("\n=== 1 Reine Zahlen (von keiner Grundgröße abhängig)")
rein = [i for i in range(len(base)) if np.all(np.abs(E[i]) < 0.004) and abs(base[i][2]) > 1e-90]
for i in rein:
    print(f"  {base[i][0]} | {base[i][1]:60s} | {base[i][2]:.6g}")

# 2 Klassen gleicher Struktur
print("\n=== 2 Klassen gleicher Abhängigkeit (Exponentenvektor)")
cls = {}
for i in range(len(base)):
    if i in rein or abs(base[i][2]) < 1e-90: continue
    key = tuple(np.round(M[i], 3))
    cls.setdefault(key, []).append(i)
for key, idx in sorted(cls.items(), key=lambda kv: -len(kv[1])):
    print(f"  [{fmt(np.array(key))}]  ({len(idx)})")
    for i in idx:
        print(f"      {base[i][0]} | {base[i][1]:58s} | {base[i][2]:.5g}")

# 3 Verknüpfungen: Paare mit proportionalen Vektoren (Exponent k aus {±1/3, ±1/2, ±2/3, ±1, ±3/2, ±2, ±3, ±4})
print("\n=== 3 Verknüpfungen zwischen Klassen: x_i / x_j^k ist rein")
keys = list(cls)
ks = [1/4, 1/3, 1/2, 2/3, 3/4, 1, 4/3, 3/2, 2, 3, 4]
for a, b in itertools.combinations(range(len(keys)), 2):
    va, vb = np.array(keys[a]), np.array(keys[b])
    for k in ks + [-x for x in ks]:
        if np.allclose(va, k * vb, atol=1e-3) and np.any(np.abs(vb) > 0):
            i, j = cls[keys[a]][0], cls[keys[b]][0]
            val = abs(base[i][2]) / abs(base[j][2])**k
            print(f"  {base[i][1][:40]:40s} = C · ({base[j][1][:40]})^{k:.3g}   [{fmt(va)}]  C = {val:.5g}")
            break

np.save(HERE + r"\ordnung_pruefstand_E.npy", E)
