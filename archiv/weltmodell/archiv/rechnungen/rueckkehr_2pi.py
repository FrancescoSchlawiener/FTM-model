"""Hypothese: kosmische Rückkehrzeit ν = 1/(1−q) = Ω_m/Ω_b = 2π (ein voller Umlauf)."""
import math, itertools
data = [("Planck 2018 TT,TE,EE+lowE+lensing", 0.1200, 0.0012, 0.02237, 0.00015, 0.0),
        ("Planck 2018 + BAO",                 0.11933, 0.00091, 0.02242, 0.00014, 0.0)]
print("1) Messwerte und Abstand zu ν = 2π  (ν = Ω_m/Ω_b = 1 + Ω_c/Ω_b)")
for n, oc, doc, ob, dob, _ in data:
    nu = 1 + oc/ob; dnu = (oc/ob)*math.hypot(doc/oc, dob/ob)
    print(f"   {n:36s} ν = {nu:.4f} ± {dnu:.4f}   2π = {2*math.pi:.4f}   Abstand {(nu-2*math.pi)/dnu:+.2f} σ")
q = 1 - 1/(2*math.pi)
print(f"\n   Folge bei ν = 2π:  q = 1 − 1/2π = {q:.5f},  Ω_c/Ω_b = 2π − 1 = {2*math.pi-1:.4f},  Ω_b/Ω_m = 1/2π = {1/(2*math.pi):.5f}")
print(f"   Vorhersage Ω_c h² = (2π−1)·Ω_b h² = {(2*math.pi-1)*0.02237:.5f} (bei Ω_b h² = 0,02237)")

print("\n2) Zufallsprüfung: 'natürliche' Kandidaten für ν im 2σ-Band des Planck+BAO-Werts [6,216 ; 6,428]")
lo, hi = 6.3224-2*0.0526, 6.3224+2*0.0526
cands = {}
consts = {"π": math.pi, "e": math.e, "φ": (1+5**.5)/2, "√2": 2**.5, "√3": 3**.5, "√5": 5**.5}
for a in range(1, 80):
    for b in range(1, 13):
        if math.gcd(a, b) == 1 and lo <= a/b <= hi: cands[f"{a}/{b}"] = a/b
for n, x in consts.items():
    for k in [1, 2, 3, 4]:
        for m in range(-6, 7):
            for val, s in [(k*x + m, f"{k}{n}{m:+d}"), (x**k + m, f"{n}^{k}{m:+d}")]:
                if lo <= val <= hi: cands[s] = val
    for (n2, x2) in consts.items():
        for val, s in [(x*x2, f"{n}·{n2}"), (x+x2, f"{n}+{n2}")]:
            if lo <= val <= hi: cands[s] = val
print(f"   {len(cands)} Kandidaten, z. B.: {', '.join(sorted(cands, key=lambda s: len(s))[:16])}")
print("   -> 2π ist einer von vielen Zahlentreffern. Es unterscheidet sich nur durch die Deutung (Kac-Rückkehr = ein Umlauf).")

print("\n3) Schärfe des Tests: Wie genau muss Ω_m/Ω_b gemessen werden, um 2π auf 5σ zu prüfen?")
for nu_true in [6.3224, 6.364]:
    d = abs(nu_true - 2*math.pi)
    print(f"   falls ν = {nu_true:.4f}: Abstand zu 2π = {d:.4f} -> nötige Unsicherheit σ ≤ {d/5:.4f} ({d/5/nu_true*100:.2f} %)")
print("   Heute: σ(ν)/ν ≈ 0,8 % (Planck+BAO). CMB-S4/LiteBIRD-Ziel für ω_c, ω_b: je ~0,2–0,3 % (aus dem Gedächtnis).")
