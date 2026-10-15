"""CAMB: CMB-Spektren für die Schließungsregel (linear = Staub mit Ω_c = Q·Ω_b) mit Q = 5,364 gegen Q = 2π−1."""
import math, numpy as np, camb
ombh2, tau, ns, As, theta = 0.02237, 0.0544, 0.9649, 2.1e-9, 0.0104092
def run(Q):
    p = camb.CAMBparams()
    p.set_cosmology(cosmomc_theta=theta, ombh2=ombh2, omch2=Q*ombh2, tau=tau, mnu=0.06)
    p.InitPower.set_params(As=As, ns=ns); p.set_for_lmax(2600, lens_potential_accuracy=1)
    r = camb.get_results(p)
    cl = r.get_cmb_power_spectra(p, CMB_unit='muK', raw_cl=True)['total']
    d = r.get_derived_params()
    return cl, p.H0, d
cases = {"Q = 5,364 (Messwert)": 5.364, "Q = 2π − 1 = 5,283": 2*math.pi - 1}
res = {k: run(v) for k, v in cases.items()}
(cl0, H0a, d0), (cl1, H0b, d1) = res[list(cases)[0]], res[list(cases)[1]]
print("Hintergrund (θ* fest):")
for k, (cl, H0, d) in res.items():
    print(f"   {k:24s} ω_c = {cases[k]*ombh2:.5f}  H0 = {H0:6.2f}  z* = {d['zstar']:.2f}  r_drag = {d['rdrag']:.3f} Mpc  z_eq = {d['zeq']:.1f}")

l = np.arange(2, 2501)
TT0, TT1 = cl0[2:2501, 0], cl1[2:2501, 0]
EE0, EE1 = cl0[2:2501, 1], cl1[2:2501, 1]
TE0, TE1 = cl0[2:2501, 3], cl1[2:2501, 3]
print("\nRelative Änderung der TT-Leistung (Q = 2π−1 gegen 5,364):")
for L in [2, 10, 30, 220, 540, 810, 1130, 1450, 2000]:
    i = L - 2; print(f"   ℓ = {L:5d}   ΔC_ℓ/C_ℓ = {(TT1[i]-TT0[i])/TT0[i]*100:+.3f} %")

# Δχ² (Fisher, ohne Verschiebung anderer Parameter -> obere Schranke), T+E mit Planck-ähnlichem Rauschen
def noise(sig_uKarcmin, fwhm_arcmin):
    s = sig_uKarcmin*math.pi/10800; b = fwhm_arcmin*math.pi/10800/math.sqrt(8*math.log(2))
    return s**2*np.exp(l*(l+1)*b**2)
def dchi2(fsky, NT, NE):
    tot = 0.0
    for i, L in enumerate(l):
        C = np.array([[TT0[i]+NT[i], TE0[i]], [TE0[i], EE0[i]+NE[i]]])
        dC = np.array([[TT1[i]-TT0[i], TE1[i]-TE0[i]], [TE1[i]-TE0[i], EE1[i]-EE0[i]]])
        Ci = np.linalg.inv(C); tot += fsky*(2*L+1)/2*np.trace(Ci@dC@Ci@dC)
    return tot
planck = dchi2(0.6, noise(33, 7.3), noise(33*math.sqrt(2)*1.0, 7.3))
cv = dchi2(0.6, np.zeros_like(TT0), np.zeros_like(EE0))
s4 = dchi2(0.4, noise(1.0, 1.4), noise(1.4, 1.4))
print(f"\nΔχ² (nur ω_c geändert, übrige Parameter fest, obere Schranke):")
print(f"   Planck-ähnlich (f_sky 0,6, 33 µK', 7,3'):   Δχ² = {planck:8.2f}  ->  {math.sqrt(planck):.2f} σ")
print(f"   ideale Messung (nur kosmische Varianz):     Δχ² = {cv:8.2f}  ->  {math.sqrt(cv):.2f} σ")
print(f"   CMB-S4-ähnlich (f_sky 0,4, 1 µK', 1,4'):    Δχ² = {s4:8.2f}  ->  {math.sqrt(s4):.2f} σ")
print("\nMit Entartungen (A_s, n_s, H0 frei) gilt die marginalisierte Planck-Unsicherheit σ(ω_c) = 0,0012:")
dwc = (5.364 - (2*math.pi-1))*ombh2
print(f"   Δω_c = {dwc:.5f}  ->  {dwc/0.0012:.2f} σ (Planck),  mit BAO σ = 0,00091 -> {dwc/0.00091:.2f} σ")
