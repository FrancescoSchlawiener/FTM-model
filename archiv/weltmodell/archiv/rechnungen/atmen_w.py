"""Dunkle Energie aus dem Strukturterm: Λ_s = (π²/15)·H²/c² (Notiz 06). Konstant oder mit dem Atmen mitlaufend?"""
import math
Om, Or = 0.3153, 9.14e-5
Os_const = math.pi**2/45            # Ω_s = Λ_s c²/(3H0²) = π²/45
print(f"1) Vakuumwert des Strukturterms: Ω_s = π²/45 = {Os_const:.4f}  (Anteil an Ω_Λ = 0,685: {Os_const/0.685:.2f})")

print("\n2) Variante 'mitlaufend' (a0 ∝ H, also Λ_s ∝ H²): konstanter Anteil Ω_s an jeder Epoche")
print(f"   Wirkung wie eine vergrößerte Gravitationskonstante für die Expansion: G_eff/G = 1/(1−Ω_s) = {1/(1-Os_const):.3f}")
print("   Urknall-Nukleosynthese erlaubt etwa |ΔG/G| ≲ 0,1 (aus dem Gedächtnis) und die CMB einen frühen Anteil ≲ 0,01–0,05")
print("   -> mitlaufende Variante ausgeschlossen.")

print("\n3) Variante 'eingefroren' (Λ_s konstant, heutiger Wert) + Rest Λ_r: w(z) der Summe bleibt −1")
print("   -> reines Λ; der DESI-Hinweis (w₀ > −1, w_a < 0) wird damit nicht erklärt.")

print("\n4) Zwischenform: Λ_s ∝ H(z)^n, 0 < n < 2 (Schwelle folgt dem Atmen nur teilweise). Effektives w(z) des Strukturterms:")
def H2(z, n, OL_rest):
    # iterativ: H² = Om(1+z)³ + Or(1+z)⁴ + OL_rest + Os·(H/H0)^n   (Einheiten H0=1)
    h2 = 1.0
    for _ in range(200):
        h2 = Om*(1+z)**3 + Or*(1+z)**4 + OL_rest + Os_const*h2**(n/2)
    return h2
for n in [0.0, 0.25, 0.5, 1.0]:
    OL_rest = 1 - Om - Or - Os_const
    ws = []
    for z in [0.0, 0.5, 1.0]:
        dz = 1e-4
        rho = lambda zz: Os_const*H2(zz, n, OL_rest)**(n/2)
        dlnrho = (math.log(rho(z+dz)) - math.log(rho(z-dz)))/(2*dz)
        ws.append(-1 + (1+z)*dlnrho/3)
    # Gesamt-w der Dunklen Energie (Λ_rest + Struktur)
    wt = []
    for z in [0.0, 0.5, 1.0]:
        dz = 1e-4
        rde = lambda zz: OL_rest + Os_const*H2(zz, n, OL_rest)**(n/2)
        wt.append(-1 + (1+z)*(math.log(rde(z+dz))-math.log(rde(z-dz)))/(2*dz)/3)
    print(f"   n = {n:4.2f}:  w_s(z=0/0,5/1) = {ws[0]:+.3f}/{ws[1]:+.3f}/{ws[2]:+.3f}   w_DE,gesamt = {wt[0]:+.3f}/{wt[1]:+.3f}/{wt[2]:+.3f}")
print("   Hinweis: n > 0 ergibt w > −1, aber mit der Zeit steigendem |w+1| in die Vergangenheit, also w_a > 0.")
print("   DESI bevorzugt w₀ > −1 UND w_a < 0 (Dunkle Energie war früher stärker 'phantomartig').")
