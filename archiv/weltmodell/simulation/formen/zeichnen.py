# -*- coding: utf-8 -*-
"""Zeichnet die tragenden Formeln des Formelblatts als Formen. Formeln und Zahlen nur aus dem Kern."""
import sys, os
sys.path.insert(0, '/home/user/FTM/weltmodell/simulation')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import kern as K

HIER = os.path.dirname(os.path.abspath(__file__))
def speichern(fig, name):
    fig.savefig(os.path.join(HIER, name), dpi=130, bbox_inches='tight'); plt.close(fig)

M_MW = 6e10 * K.M_sun
A0 = K.A0
def T(x):  # Takt, x = r/Δs, nur K.takt-Formel: 1 - 6Δs/r
    return np.sqrt(np.clip(1 - 6.0 / x, 0, None))
# Prüfung gegen Kern: Takt in Kernfunktion bei r = 9Δs
_M = K.M_sun; _r = 9 * K.ds(_M)
assert abs(K.takt(_M, _r) - T(9.0)) < 1e-12

# 1 Trichter -----------------------------------------------------------------
def trichter():
    fig = plt.figure(figsize=(13, 5.5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    r = np.linspace(6, 40, 200); ph = np.linspace(0, 2*np.pi, 90)
    R, P = np.meshgrid(r, ph)
    Z = 2*np.sqrt(6*(R-6))            # Flamm: z' = 1/sqrt(r/6Δs -1)... = ∫√g_rr-Ausdruck, Einheit Δs
    ax.plot_surface(R*np.cos(P), R*np.sin(P), -Z, cmap='viridis', alpha=.6, linewidth=0, rstride=2, cstride=2)
    for rk, f, n in [(6, 'r', '6Δs'), (9, 'orange', '9Δs'), (18, 'k', '18Δs')]:
        z = -2*np.sqrt(6*(rk-6)); t = np.linspace(0, 2*np.pi, 200)
        ax.plot(rk*np.cos(t), rk*np.sin(t), z*np.ones_like(t), color=f, lw=2.5, label=n)
    ax.set_title('Flamm-Paraboloid  z = 2√(6Δs(r−6Δs))'); ax.legend(); ax.set_xlabel('r/Δs'); ax.set_zlabel('z/Δs')
    ax2 = fig.add_subplot(1, 2, 2)
    x = np.linspace(6.0001, 60, 800); ax2.plot(x, T(x), lw=2, label='T = √(1−6Δs/r)')
    ax2.plot(x, T(x)/x*6, '--', label='T/r (skaliert)')
    ax2.plot(x, 1/T(x)/4, ':', label='Ringdichte 1/T (÷4)')
    for rk, c, n in [(6, 'r', '6Δs: T=0'), (9, 'orange', '9Δs: 1/T=√3, T/r max'), (18, 'k', '18Δs: κ=0')]:
        ax2.axvline(rk, color=c, ls=':'); ax2.text(rk+.5, .05, n, rotation=90, color=c, va='bottom')
    ax2.axhline(1, color='gray', lw=.5); ax2.set_xlabel('r/Δs'); ax2.set_ylim(0, 1.2); ax2.legend(loc='center right')
    ax2.set_title('Takt: Wurzel, senkrechte Tangente bei 6Δs, Asymptote 1')
    speichern(fig, 'trichter.png')

# 2 q, nu ---------------------------------------------------------------------
def q_nu():
    s = np.linspace(0.02, 6, 600)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].plot(s, np.exp(-s), label='q = e^{−s}'); ax[0].plot(s, 1-np.exp(-s), '--', label='1−q')
    ax[0].axvline(1, color='gray', ls=':'); ax[0].plot(1, np.exp(-1), 'ro'); ax[0].text(1.1, .4, 's=1 (r*): q=0,368')
    ax[0].set_xlabel('s = √(a_N/a₀) = r*/r'); ax[0].legend(); ax[0].set_title('Exponentialkurve, q(0)=1, Asymptote 0')
    nu = 1/(1-np.exp(-s))
    ax[1].plot(s, nu, label='ν = 1/(1−e^{−s})'); ax[1].plot(s, 1/s, '--', label='1/s (Hyperbel, tief)')
    ax[1].plot(s, 1+1/np.expm1(s), ':', label='1+n_B')
    ax[1].axhline(1, color='gray', lw=.6); ax[1].plot(1, 1/(1-np.exp(-1)), 'ro'); ax[1].text(1.1, 1.9, 's=1: ν=1,582')
    ax[1].set_ylim(0.8, 6); ax[1].set_xlabel('s'); ax[1].legend(); ax[1].set_title('Hyperbel über der Exponentialkurve')
    speichern(fig, 'q_nu.png')

# 3 Regel ---------------------------------------------------------------------
def regel():
    x = np.logspace(-4, 4, 600); a = K.a_regel(x*A0)/A0
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.loglog(x, a, lw=2.5, label='a/a₀ = ν·a_N/a₀'); ax.loglog(x, x, '--', label='Newton: a = a_N (Steigung 1)')
    ax.loglog(x, np.sqrt(x), '--', label='Wurzel: a = √(a_N a₀) (Steigung ½)')
    ax.plot(1, K.a_regel(A0)/A0, 'ro'); ax.text(1.3, 1.9, 'Knick bei a_N=a₀\nSchnitt der Asymptoten; a=1,58 a₀')
    ax.set_xlabel('a_N/a₀'); ax.set_ylabel('a/a₀'); ax.legend(); ax.grid(alpha=.3, which='both'); ax.set_title('Regel: Hyperbel∘Exponential∘Wurzel')
    speichern(fig, 'regel.png')

# 4 Rotationskurve ------------------------------------------------------------
def rotation():
    kpc = K.kpc; r = np.logspace(-1, 2.3, 500)*kpc
    v = np.array([K.v_kreis_punkt(M_MW, x) for x in r])/1e3
    vN = np.sqrt(K.G*M_MW/r)/1e3; nu = K.nu(K.G*M_MW/r**2)
    rs = K.r_stern(M_MW)/kpc; vf = K.v_flach(M_MW)/1e3
    assert np.allclose(v, vN*np.sqrt(nu), rtol=1e-9)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogx(r/kpc, v, lw=2.5, label='v = v_Kepler·√ν'); ax.semilogx(r/kpc, vN, '--', label='Kepler v∝r^{−1/2}')
    ax.axhline(vf, color='g', ls=':', label='flach: v⁴=GMa₀ = %.0f km/s' % vf)
    ax.axvline(rs, color='r', ls=':'); ax.text(rs*1.05, 20, 'r* = %.2f kpc' % rs, color='r')
    i = np.argmin(v); ax.plot(r[i]/kpc, v[i], 'ko'); ax.text(r[i]/kpc*1.1, v[i]+20, 'Minimum %.0f km/s bei %.0f kpc' % (v[i], r[i]/kpc))
    ax.set_ylim(0, 600); ax.set_xlabel('r [kpc] (r/r* = r/%.2f)' % rs); ax.set_ylabel('v [km/s]'); ax.legend(); ax.set_title('Rotationskurve der Punktmasse 6·10¹⁰ M☉')
    speichern(fig, 'rotation.png')
    return rs, vf, v[i], r[i]/kpc

# 5 Rückkehr û über p ---------------------------------------------------------
def rueckkehr():
    tD = K.tau_D(); w = np.logspace(-3, 3, 800)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    for G_, c in [(0, 'C0'), (1, 'C1'), (4, 'C2')]:
        u = K.rueckkehr(1j*w/tD, G_/tD)
        ax[0].loglog(w, np.abs(u), c, label='s²=Γτ_D=%g' % G_); ax[1].semilogx(w, np.unwrap(np.angle(u)), c)
        ax[2].plot(u.real, u.imag, c)
    ax[0].set_xlabel('ωτ_D  (p = iω)'); ax[0].set_ylabel('|û|'); ax[0].legend(); ax[0].set_title('Betrag: e^{−Re√} , Abfall bei ωτ_D ≈ 1')
    ax[1].set_xlabel('ωτ_D'); ax[1].set_ylabel('arg û [rad]'); ax[1].set_title('Phase: −Im√((iω+Γ)τ_D), Windung')
    ax[2].set_aspect('equal'); ax[2].set_title('û in der Ebene: Spirale nach 0'); ax[2].set_xlabel('Re û'); ax[2].set_ylabel('Im û')
    p = np.linspace(0, 10, 200)
    ax[2].plot([], [], ' ', label='p=0,Γ=0: û=1'); ax[2].legend()
    speichern(fig, 'rueckkehr.png')
    # reelles p: Laplace-Kern
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for G_ in (0, 1, 4): ax.semilogy(p, K.rueckkehr(p/tD, G_/tD), label='Γτ_D=%g' % G_)
    ax.set_xlabel('pτ_D'); ax.set_ylabel('û'); ax.legend(); ax.set_title('û über reellem p: e^{−√(pτ_D)}')
    speichern(fig, 'rueckkehr_reell.png')

# 6 chi -----------------------------------------------------------------------
def chi():
    s = np.linspace(1e-4, 10, 1000); c = s/np.expm1(s)
    d2 = np.gradient(np.gradient(c, s), s); iw = s[1:-1][np.where(np.diff(np.sign(d2[1:-1])))[0]]
    fig, ax = plt.subplots(figsize=(6.5, 4.5)); ax.plot(s, c, lw=2.5, label='χ = s/(e^s−1)')
    ax.plot(s, np.exp(-s), '--', label='e^{−s}'); ax.plot(s, 1/(1+s), ':', label='1/(1+s)')
    ax.axhline(1, color='gray', lw=.5); ax.plot(1, 1/(np.e-1), 'ro'); ax.text(1.1, .65, 'χ(1)=0,582')
    ax.set_xlabel('s'); ax.set_ylabel('χ'); ax.legend(); ax.set_title('χ: Schalter 1→0, s·n_B(s), Sigmoide (Bernoulli-Erzeugende)')
    speichern(fig, 'chi.png'); return iw

# 7 Leiter --------------------------------------------------------------------
def leiter():
    M = np.logspace(-30, 45, 400)
    ds = K.G*M/(3*K.c**2); rs = np.sqrt(K.G*M/A0); l0 = K.ell0()
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].loglog(M, rs, lw=2.5, label='r* = √(GM/a₀): Steigung ½'); ax[0].loglog(M, 3*ds, label='3Δs: Steigung 1'); ax[0].loglog(M, l0*np.ones_like(M), 'k--', label='ℓ₀=2Λ=c²/a₀ (waagerecht)')
    for n, m in [('Erde', K.M_earth), ('Sonne', K.M_sun), ('MW', M_MW)]:
        ax[0].plot(m, K.r_stern(m), 'ro'); ax[0].text(m*1.5, K.r_stern(m)*.4, n)
    ax[0].set_xlabel('M [kg]'); ax[0].set_ylabel('Länge [m]'); ax[0].legend(); ax[0].set_title('Leiter: 3Δs : r* : ℓ₀ geometrisch (Verhältnis δ*=3Δs/r*=r*/ℓ₀)')
    x = np.logspace(-8, 12, 10); ax[1].loglog(x, l0*x, lw=2.5, label='r*² = ℓ₀·(3Δs): Gerade, Steigung 1')
    lam = lambda m: K.hbar/(m*K.c)
    ax[1].loglog(x, K.l_P**2/x, '--', label='ℓ_P² = 3Δs·λ̄ als λ̄(3Δs): Steigung −1')
    for n, m in [('Sonne', K.M_sun), ('MW', M_MW)]:
        ax[1].plot(3*K.ds(m), K.r_stern(m)**2, 'ro'); ax[1].text(3*K.ds(m)*1.5, K.r_stern(m)**2*.3, n)
    ax[1].plot(3*K.ds(K.m_p), (3*K.ds(K.m_p))*lam(K.m_p), 'bs')
    ax[1].set_xlabel('3Δs [m]'); ax[1].set_ylabel('Länge² [m²]'); ax[1].legend(); ax[1].set_title('Gerade im log-log; Achsenabschnitt log ℓ₀')
    for m in (K.M_sun, M_MW): assert abs(K.r_stern(m)**2/(3*K.ds(m)*l0) - 1) < 1e-9
    speichern(fig, 'leiter.png')

# 8 Nome ---------------------------------------------------------------------
def nome():
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    th = np.linspace(0, 2*np.pi, 300)
    for y in (0.05, 0.15, 0.3, 0.6): ax[0].plot(np.exp(-2*np.pi*y)*np.cos(th), np.exp(-2*np.pi*y)*np.sin(th), label='Im τ=%g: |q|=%.2f' % (y, np.exp(-2*np.pi*y)))
    for x in np.linspace(0, 1, 9)[:-1]:
        y = np.linspace(0.02, 0.8, 100); qq = np.exp(2j*np.pi*(x+1j*y)); ax[0].plot(qq.real, qq.imag, 'gray', lw=.5)
    ax[0].set_aspect('equal'); ax[0].legend(fontsize=7); ax[0].set_title('q = e^{2πiτ}: Kreisring, |q|=e^{−2π Im τ}, arg q = 2π Re τ')
    t = np.linspace(0, 3, 600); qq = np.exp(2j*np.pi*t*(1+0.25j))
    ax[1].plot(qq.real, qq.imag); ax[1].set_aspect('equal'); ax[1].set_title('logarithmische Spirale: τ = t(1+0,25i)')
    s = np.linspace(0.05, 6, 400); ax[2].plot(s, np.exp(-s)); ax[2].plot(s, np.exp(-s)*np.cos(2*np.pi*0.2*s)); ax[2].set_title('reeller Schnitt (s = 2π Im τ): q(s)=e^{−s}; mit Windung gedämpfter Kosinus')
    ax[2].set_xlabel('s')
    speichern(fig, 'nome.png')

# 9 Gauss-Bonnet --------------------------------------------------------------
def gauss_bonnet():
    x = np.linspace(6.0, 60, 800); gb = 2*np.pi*(1-T(x))
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.plot(x, gb, lw=2.5, label='∫K dA = 2π(1−T)')
    ax.plot(x, 2*np.pi*3/x, '--', label='2π·3Δs/r (schwach)'); ax.axhline(2*np.pi, color='r', ls=':'); ax.text(30, 2*np.pi+.1, '2π: Topologiewechsel (T=0)', color='r')
    for rk in (9, 18): ax.plot(rk, 2*np.pi*(1-T(rk)), 'ko'); ax.text(rk+.5, 2*np.pi*(1-T(rk))+.15, '%gΔs: %.2f' % (rk, 2*np.pi*(1-T(rk))))
    ax.set_xlabel('r/Δs'); ax.set_ylabel('∫K dA'); ax.legend(); ax.set_title('Gauß-Bonnet: 2π minus Wurzel, fällt wie 1/r')
    speichern(fig, 'gauss_bonnet.png')

# 10 Q, nu als Spirale --------------------------------------------------------
def spirale():
    s = np.linspace(0.15, 8, 800); k = 1.6
    q = np.exp(-s*(1 - 1j*k)/1.0); Q = q/(1-q); nu = 1/(1-q)
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    ax[0].plot(q.real, q.imag); ax[0].set_title('q(s) = e^{−s(1−ik)}: logarithmische Spirale'); ax[0].set_aspect('equal')
    ax[1].plot(nu.real, nu.imag); ax[1].plot(1, 0, 'ro'); ax[1].set_title('ν = 1/(1−q): Kreisinversion, Spirale um 1 (Wurzel: ν→1)'); ax[1].set_aspect('equal')
    ax[2].plot(Q.real, Q.imag); ax[2].plot(0, 0, 'ro'); ax[2].set_title('Q = q/(1−q) = ν−1: Spirale um 0'); ax[2].set_aspect('equal')
    Qm = K.OMEGA_C_H2/K.OMEGA_B_H2; qm = Qm/(1+Qm)
    ax[2].plot(Qm, 0, 'gs'); ax[2].text(Qm, .3, 'Q=%.3f\n(s=%.3f)' % (Qm, -np.log(qm)))
    speichern(fig, 'spirale_Q_nu.png')
    return Qm, qm

# 11 Abhängigkeit -------------------------------------------------------------
def abhaengigkeit():
    fig, ax = plt.subplots(figsize=(17, 10)); ax.set_xlim(0, 17); ax.set_ylim(0, 10); ax.axis('off')
    N = {}
    def knoten(k, x, y, t, farbe):
        ax.add_patch(FancyBboxPatch((x-.95, y-.32), 1.9, .64, boxstyle='round,pad=0.03', fc=farbe, ec='k')); ax.text(x, y, t, ha='center', va='center', fontsize=8); N[k] = (x, y)
    def pf(a, b):
        ax.add_patch(FancyArrowPatch((N[a][0], N[a][1]-.32), (N[b][0], N[b][1]+.32), arrowstyle='-|>', mutation_scale=12, lw=1, color='#444', shrinkA=0, shrinkB=0))
    L, I, Z, Gc = '#cfe8ff', '#ffe0b3', '#d9f2d0', '#f2d0ef'
    # Ebene 0 Eingaben
    knoten('M', 2, 9.3, 'M (Masse)\nGerade', L); knoten('r', 5, 9.3, 'r (Ort)', L); knoten('a0', 9, 9.3, 'a₀ (Konstante)', I); knoten('tau', 14.5, 9.3, 'τ = Re+i·Im\nGerade/Halbebene', Z)
    # Ebene 1
    knoten('ds', 2, 8, 'Δs = GM/3c²\nGerade Steigung 1', L); knoten('aN', 5, 8, 'a_N = GM/r²\nHyperbel 1/r²', L)
    knoten('Lam', 9, 8, 'Λ = c²/2a₀\nHyperbel 1/a₀', I); knoten('nome', 14.5, 8, 'q = e^{2πiτ}\nKreisring/Spirale', Z)
    # Ebene 2
    knoten('T', 0.9, 6.6, 'T = √(1−6Δs/r)\nWurzel', L); knoten('s', 5, 6.6, 's = √(a_N/a₀)\nWurzel', L)
    knoten('tD', 9, 6.6, 'τ_D = 2Λ/c\nGerade in Λ', I); knoten('l0', 11.2, 6.6, 'ℓ₀ = 2Λ', I)
    # Ebene 3
    knoten('flamm', 0.9, 5.2, 'Flamm-Trichter\nParaboloid z∝√(r−6Δs)', L); knoten('gb', 3.3, 5.2, 'Gauß-Bonnet\n2π(1−T)', L)
    knoten('q', 5.6, 5.2, 'q = e^{−s}\nExponentialkurve', I); knoten('rs', 8, 5.2, 'r* = √(GM/a₀)\nWurzel, Steigung ½', L)
    knoten('u', 10.8, 5.2, 'û = e^{−√((p+Γ)τ_D)}\nExp∘Wurzel, Spirale', I); knoten('leiter', 13.8, 5.2, 'Leiter r*²=3Δs·ℓ₀\nGerade log-log', L)
    # Ebene 4
    knoten('kreise', 0.9, 3.7, 'Kreise 6/9/18Δs\nStufen', L); knoten('nB', 4.2, 3.7, 'n_B = 1/(e^s−1)\nBose-Hyperbel', I)
    knoten('nu', 7, 3.7, 'ν = 1/(1−q)\nHyperbel∘Exp', I); knoten('Q', 9.5, 3.7, 'Q = q/(1−q) = ν−1\nSpirale', I)
    knoten('vK', 12.2, 3.7, 'v_Kepler = √(GM/r)\nWurzel 1/√r', L)
    # Ebene 5
    knoten('chi', 3.2, 2.2, 'χ = s·n_B\nSigmoide', I); knoten('regel', 6.5, 2.2, 'Regel a = ν·a_N\nKnick, Asymptoten\nNewton & Wurzel', I)
    knoten('Qspir', 9.5, 2.2, 'Q,ν als Spirale\nlog. Spirale', I); knoten('rot', 12.2, 2.2, 'Rotationskurve\nv = v_K·√ν', L)
    kanten = [('M','ds'),('r','aN'),('M','aN'),('a0','Lam'),('tau','nome'),('ds','T'),('r','T'),('aN','s'),('a0','s'),('Lam','tD'),('Lam','l0'),
              ('T','flamm'),('T','gb'),('flamm','kreise'),('s','q'),('ds','rs'),('a0','rs'),('tD','u'),('s','u'),('l0','leiter'),('ds','leiter'),('rs','leiter'),
              ('q','nu'),('q','nB'),('s','nB'),('q','Q'),('M','vK'),('nome','q'),('nB','chi'),('s','chi'),('nu','regel'),('aN','regel'),
              ('nu','rot'),('vK','rot'),('Q','Qspir'),('nu','Qspir'),('u','q')]
    for a, b in kanten: pf(a, b)
    ax.text(.1, .3, 'Farben: blau = Ebene L (lokal, Faltung/Trichter) · orange = I (innen, Rückkehr) · grün = Z (Nome/Windung).  Pfeil = geht ein in.\nDie zwei Rückkehr-Pfeile nome→q, u→q: q ist Sonderfall (p=0,Γ=s²/τ_D) der einen Rückkehrgleichung û und Betrag des Nomes.', fontsize=9)
    speichern(fig, 'abhaengigkeit.png')

if __name__ == '__main__':
    trichter(); q_nu(); regel(); r = rotation(); rueckkehr(); iw = chi(); leiter(); nome(); gauss_bonnet(); Qm = spirale(); abhaengigkeit()
    print('r*=%.3f kpc  v_flach=%.1f  vmax=%.1f bei %.2f kpc' % r); print('chi Wendepunkt s=', iw); print('Q,q=', Qm)
