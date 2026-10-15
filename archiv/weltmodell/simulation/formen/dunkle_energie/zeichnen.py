"""Formen der Dunklen Energie: Kurven zeichnen, Merkmale rechnen. Formeln nur aus Kern und forschung/kosmos.py.
Setzung: a0 = 1,1 (Standard), n_Z = X0/ln(nu_kosm) (nz/ergebnis), sofortige Verschiebung (M1), Mmin = FLUCHT_M_MIN."""
import sys, math, os, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/user/FTM/weltmodell/simulation'); sys.path.insert(0, '/home/user/FTM/weltmodell/simulation/forschung')
import kern as K, forschung.kosmos as Ko
from kern.eingabe import EINGABEN
from scipy.special import erfc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
OUT = os.path.dirname(os.path.abspath(__file__)); c = K.c
X0 = 4.5757
q = K.q_KOSM; nu = 1/(1-q); nZ = X0/math.log(nu)
Mmin = EINGABEN["FLUCHT_M_MIN"].wert
E_L = lambda a: K.E_hintergrund(1/a-1)                      # LambdaCDM (Kern)
Eb = lambda a: K.hintergrund_modell(1/a-1)/K.H0_si()         # Hintergrund mit Vakuumteil
sM = float(Ko.streuung_masse(np.array([math.log10(Mmin)]))[0]); D = Ko.wachstum_linear(E_L); D1 = D(1.0)[0]
F = lambda a: erfc(1.686/(math.sqrt(2)*sM*D(a)[0]/D1))       # Press-Schechter mit Kernhilfen
Ovac = K.lambda_vakuum(K.schliesslaenge(K.A0))*c**2/3/K.H0_si()**2   # rho_s in H0^2-Einheiten (konstant)
Om = K.OMEGA_M
n_ = lambda a: nZ*F(a)/F(1.0)
rho_DE = lambda a: Ovac*(1+n_(a))
E2 = lambda a: Eb(a)**2 - Ovac + rho_DE(a)                    # Vakuumteil durch rho_DE ersetzt (M1, a0 = 1,1 fest)
lna = np.linspace(-3.0, 1.0, 401); a = np.exp(lna); z = 1/a-1
Fa = np.array([F(x) for x in a]); dF = np.gradient(Fa, lna)
rho = np.array([rho_DE(x) for x in a]); Ef = np.sqrt([E2(x) for x in a]); El = np.array([E_L(x) for x in a])
Th = 2*math.pi*K.A0/(c*Ef*K.H0_si())                         # Theta(z) aus H(z) des Modells
Os = math.pi**2/45*Th**2
w = -1-np.gradient(np.log(rho), lna)/3
ODE = rho/Ef**2
res = {}
# Merkmale
i = np.argmax(dF); res['Wendepunkt lna'] = lna[i]; res['a_w'] = a[i]; res['z_w(lna)'] = z[i]
dFz = np.abs(np.gradient(Fa, z)); mask = (z >= 0) & (z < 6); j = np.argmax(np.where(mask, dFz, 0)); res['z max |dF/dz|'] = z[j]
d2 = np.gradient(dF, lna); ok = np.where(np.diff(np.sign(d2)) != 0)[0]; res['Nulldurchgang d2F'] = [(lna[k], z[k]) for k in ok]
res['a_gl'] = (Ovac/Om)**(-1/3) if False else (Om/Ovac)**(1/3); res['z_gl'] = 1/res['a_gl']-1   # rho_m = rho_s
res['a_gl_DE'] = (Om/rho_DE(1.0))**(1/3)
res['Theta0'] = float(Th[-1] if abs(a[-1]-1) < 1e-9 else np.interp(0, lna, Th)); res['Os0'] = float(np.interp(0, lna, Os)); res['ODE0'] = float(np.interp(0, lna, ODE))
res['w0'] = float(np.interp(0, lna, w)); res['nZ'] = nZ; res['kDE/ks'] = 1+nZ
res['F(0),F(z=1),F(z=3)'] = (F(1.0), F(.5), F(.25))
res['H/H_LCDM (z=0,1,2)'] = [float(np.interp(l, lna, Ef/El)) for l in (0, -math.log(2), -math.log(3))]
res['w(z=0,.5,1,2)'] = [float(np.interp(-math.log(1+zz), lna, w)) for zz in (0, .5, 1, 2)]
zm = np.argmax(np.where(z >= 0, ODE*0+rho*Ef**0, 0)); 
res['z Schnitt F=1/2'] = float(np.interp(0.5, Fa[::-1], z[::-1]))
for k, v in res.items(): print(k, v)
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .25})
col = dict(s='#1f6f8b', m='#b5651d', f='#7a3e9d', de='#c0392b', k='#444444')
zp = (z >= 0) & (z <= 6)
# Bild 1: Formen
fig, ax = plt.subplots(3, 3, figsize=(15, 12)); ax = ax.ravel()
t = np.linspace(0.2, 2.2, 200); ax[0].plot(t, math.pi**2/45*t**2, col['s']); ax[0].axvline(Th[np.argmin(abs(lna))], color='k', ls=':')
ax[0].set(title='Vakuumteil: Parabel  Ω_s = (π²/45)Θ²', xlabel='Θ', ylabel='Ω_s'); ax[0].plot(res['Theta0'], res['Os0'], 'o', color=col['de'])
ax[1].plot(z[zp], Os[zp], col['s'], label='Ω_s(z)'); ax[1].plot(z[zp], Th[zp], col['k'], ls='--', label='Θ(z)')
ax[1].set(title='Ω_s über z: Θ ∝ 1/H fällt, Ω_s ∝ 1/H² (Hyperbel in H)', xlabel='z'); ax[1].legend()
ax[2].loglog(a, rho*0+Ovac, col['s'], label='ρ_s = const (Stufe/Gerade)'); ax[2].loglog(a, Om*a**-3, col['m'], label='ρ_m ∝ a⁻³ (Hyperbel)')
ax[2].loglog(a, rho, col['de'], label='ρ_DE'); ax[2].axvline(res['a_gl'], color='k', ls=':'); ax[2].axvline(res['a_gl_DE'], color='k', ls='--'); ax[2].set(title=f"ρ_m=ρ_s: z={res['z_gl']:.2f}; ρ_m=ρ_DE: z={1/res['a_gl_DE']-1:.2f}", xlabel='a', ylim=(.1, 1e4)); ax[2].legend()
ax[3].plot(lna, Fa, col['f']); ax[3].axvline(res['Wendepunkt lna'], color='k', ls=':'); ax[3].set(title='Bindungsanteil F(a): Sigmoide (Fehlerfunktion)', xlabel='ln a', ylabel='F')
ax[4].plot(lna, dF, col['f']); ax[4].axvline(res['Wendepunkt lna'], color='k', ls=':'); ax[4].set(title=f"dF/dln a: Glocke, Gipfel = Wendepunkt (z={res['z_w(lna)']:.2f})", xlabel='ln a')
ax[5].plot(z[zp], rho[zp]/Ovac, col['de']); ax[5].set(title=f'ρ_DE/ρ_s = 1 + n_Z·F/F(1): gehobene Sigmoide (n_Z={nZ:.2f})', xlabel='z'); ax[5].axhline(1, color='k', ls=':')
ax[6].plot(z[zp], w[zp], col['de']); ax[6].axhline(-1, color='k', ls=':'); ax[6].set(title='w(z) = −1 − dlnρ/(3dln a): flache Stufe unter −1', xlabel='z')
ax[7].plot(z[zp], (Ef/El)[zp], col['k']); ax[7].plot(z[zp], (ODE)[zp], col['de'], label='Ω_DE(z)'); ax[7].legend(['H/H_ΛCDM', 'Ω_DE(z)']); ax[7].set(title='H(z)/H_ΛCDM (Beule) und Ω_DE(z) (Sigmoide, fallend)', xlabel='z')
n = np.linspace(0, 6, 300); ax[8].semilogy(n, (1-q)**n, col['k']); ax[8].axhline(math.exp(-X0), color=col['de'], ls=':'); ax[8].axvline(nZ, color=col['de'], ls=':')
ax[8].set(title=f'Kac-Zählung (1−q)ⁿ = e^(−X): Exponential, n_Z={nZ:.2f}', xlabel='n', ylim=(1e-3, 1))
plt.tight_layout(); plt.savefig(f'{OUT}/formen_kurven.png', dpi=110); plt.close()
# Bild 2: Abhängigkeitsbild
fig, ax = plt.subplots(figsize=(14, 8)); ax.axis('off'); ax.set_xlim(0, 14); ax.set_ylim(0, 8)
def box(x, y, txt, fc, w_=3.2, h=1.0):
    ax.add_patch(FancyBboxPatch((x-w_/2, y-h/2), w_, h, boxstyle='round,pad=0.05', fc=fc, ec='k', alpha=.85)); ax.text(x, y, txt, ha='center', va='center', fontsize=8.5)
def pf(a_, b_, t=''):
    ax.annotate('', xy=b_, xytext=a_, arrowprops=dict(arrowstyle='->', lw=1.3)); 
    if t: ax.text((a_[0]+b_[0])/2, (a_[1]+b_[1])/2+.12, t, fontsize=7.5, ha='center')
P = {'H': (2, 7), 'Th': (2, 5.5), 'Os': (2, 4), 'rs': (2, 2.5), 'Mm': (6.2, 7), 'z': (6.2, 5.5), 'F': (6.2, 3.8), 'kac': (10.5, 7), 'nZ': (10.5, 5.5), 'DE': (8.5, 1.9), 'w': (12.5, 3.7), 'HH': (12.5, 1.9)}
box(*P['H'], 'H(z) (G): Hyperbel 1/H', '#dfe')
box(*P['Th'], 'Θ = 2πa₀/(cH): Hyperbel\n(Λ = c²/2a₀: Höhe, innere Richtung)', '#cde')
box(*P['Os'], 'Ω_s = (π²/45)Θ²: Parabel', '#cde')
box(*P['rs'], 'ρ_s = const: Gerade\n(Ω_s·H² ∝ Θ²H² = const)', '#cde')
box(*P['Mm'], 'Mmin, D(a): Lage in z\n(Wachstum, Streuung σ(M))', '#fed')
box(*P['z'], 'ρ_m ∝ a⁻³: Hyperbel in a\nGleichheit mit ρ_DE: z≈%.2f\n(mit ρ_s allein: z=%.2f)' % (1/res['a_gl_DE']-1, res['z_gl']), '#fed')
box(*P['F'], 'F(a) = erfc(δ_c/√2σD): Sigmoide\nWendepunkt z≈%.2f' % res['z_w(lna)'], '#edc', h=1.1)
box(*P['kac'], 'Kac: (1−q)ⁿ = e^(−X): Exponential', '#eee'); box(*P['nZ'], 'n_Z = X/ln ν = %.2f: Höhe der Stufe' % nZ, '#eee')
box(*P['DE'], 'ρ_DE = ρ_s[1 + n_Z F/F(1)]\ngehobene Sigmoide auf Gerade', '#fcc', w_=3.6, h=1.2)
box(*P['w'], 'w(z): Stufe unter −1', '#fcc', w_=2.8); box(*P['HH'], 'H/H_ΛCDM: Beule, Ω_DE(z)', '#fcc', w_=3.2)
pf((2, 6.5), (2, 6.0)); pf((2, 5.0), (2, 4.5)); pf((2, 3.5), (2, 3.0))
pf((3.6, 2.5), (7, 2.1), 'Höhe (Grundlinie)'); pf((6.2, 3.25), (7.6, 2.5), 'Gestalt'); pf((6.2, 6.5), (6.2, 6.0), 'Lage'); pf((6.2, 5.0), (6.2, 4.4), 'D(a)'); pf((10.5, 5.0), (9.5, 2.5), 'Höhe (Faktor)')
pf((10.5, 6.5), (10.5, 6.0)); pf((10.3, 1.9), (10.9, 1.9)); pf((9.5, 2.4), (11.6, 3.3)); pf((3.6, 7), (4.6, 5.8), 'H, Ω_m')
pf((3.6, 7.0), (3.6, 7.0)); ax.set_title('Abhängigkeit der Formen der Dunklen Energie', fontsize=12); plt.savefig(f'{OUT}/abhaengigkeit.png', dpi=110, bbox_inches='tight'); plt.close()
# Bild 3: Wendepunkt gegen z≈0,4
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(z[zp], dFz[zp]*0+np.abs(np.gradient(Fa, z))[zp], col['f']); ax[0].axvline(0.4, color='k', ls='--', label='gemessen z≈0,4 [G]'); ax[0].axvline(res['z_w(lna)'], color=col['de'], ls=':', label='Wendepunkt in ln a'); ax[0].legend(); ax[0].set(title='|dF/dz|', xlabel='z')
ax[1].plot(z[zp], dF[zp], col['f']); ax[1].axvline(0.4, color='k', ls='--'); ax[1].axvline(res['z_w(lna)'], color=col['de'], ls=':'); ax[1].set(title='dF/dln a (Rate der Bindung)', xlabel='z')
plt.tight_layout(); plt.savefig(f'{OUT}/wendepunkt.png', dpi=110); plt.close()
import json; json.dump({k: (v if not isinstance(v, (np.floating,)) else float(v)) for k, v in res.items()}, open(f'{OUT}/merkmale.json', 'w'), default=float, ensure_ascii=False, indent=1)
