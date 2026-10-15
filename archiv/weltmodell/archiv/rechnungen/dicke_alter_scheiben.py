# Dicke alter Scheiben unter der örtlichen Regel F = −ν(|g_N|)∇Φ_N (nicht wirbelfrei).
# Vergleich: wirbelfreie Kraft F_c = −ν_m(Φ_N)∇Φ_N, ν_m aus der Mittelebene; in der Ebene gleich, gleiche Rotationskurve.
# Milchstraßen-ähnlich: Exponentialscheibe 5·10¹⁰ M☉ (R_d = 2,6 kpc, exponentiell in z mit h), Hernquist-Bulge 10¹⁰ M☉ (a = 0,6 kpc).
# Einheiten kpc, km/s, M☉; Zeit in kpc/(km/s) = 977,8 Myr. a0 = 1,2·10⁻¹⁰ m/s².
import numpy as np, sys
from scipy.special import j0, j1
from scipy.ndimage import map_coordinates
G = 4.30091e-6; A0 = 1.2e-10*3.0857e19/1e6; TU = 977.8
Md, Rd, Mb, ab = 5e10, 2.6, 1e10, 0.6
def nu(g): s = np.sqrt(np.maximum(g, 1e-30)/A0); return 1/(-np.expm1(-s))

def newton_grid(h, Rg, zg):
    k = np.concatenate([np.linspace(1e-5, 2.0, 6000), np.geomspace(2.0, 150.0, 6000)[1:]])
    w = np.gradient(k)                                    # Quadraturgewichte
    S = Md/(2*np.pi)*(1 + (k*Rd)**2)**-1.5               # Σ̃(k)·(1/2πG)-Normierung: 2πG Σ0 Rd² = G Md
    kh = k*h; den = 1 - kh**2; den = np.where(np.abs(den) < 1e-6, 1e-6, den)
    Zp = (np.exp(-np.outer(zg, k)) - kh*np.exp(-zg/h)[:, None])/den          # Φ-Faktor
    Zz = (-k*np.exp(-np.outer(zg, k)) + k*np.exp(-zg/h)[:, None])/den        # ∂/∂z des Faktors
    gR = np.zeros((len(Rg), len(zg))); gz = np.zeros_like(gR); ph = np.zeros_like(gR)
    for i, R in enumerate(Rg):
        a = 2*np.pi*G*S*w
        ph[i] = -(Zp @ (j0(k*R)*a)); gR[i] = -(Zp @ (k*j1(k*R)*a)); gz[i] = (Zz @ (j0(k*R)*a))
    r = np.sqrt(Rg[:, None]**2 + zg[None, :]**2)
    ph += -G*Mb/(r + ab); gb = -G*Mb/(r + ab)**2
    gR += gb*Rg[:, None]/r; gz += gb*zg[None, :]/r
    return ph, gR, gz

def run(h, sz0, N=400, T=10000.0, dt=0.25, seed=1):
    Rg = np.linspace(0.3, 20.0, 395); zg = np.linspace(0.0, 4.0, 401)
    ph, gR, gz = newton_grid(h, Rg, zg)
    g = np.hypot(gR, gz); nl = nu(g)
    # ν_m(Φ): Mittelebenenwert auf derselben Newton-Äquipotentialfläche
    Pm = ph[:, 0]; nm = nu(np.abs(gR[:, 0]))
    order = np.argsort(Pm); Ps, ns = Pm[order], nm[order]
    nc = np.interp(ph, Ps, ns)
    U1 = np.concatenate([[0], np.cumsum(0.5*(ns[1:] + ns[:-1])*np.diff(Ps))])   # U(Φ) = ∫ν_m dΦ
    U = np.interp(ph, Ps, U1)
    FRl, Fzl = nl*gR, nl*gz; FRc, Fzc = nc*gR, nc*gz
    vc = np.sqrt(np.maximum(-Rg*FRl[:, 0], 0))
    iR = lambda R: (R - Rg[0])/(Rg[1] - Rg[0]); iz = lambda z: np.abs(z)/(zg[1] - zg[0])
    def acc(FR, Fz, R, z, Lz):
        c = np.vstack([iR(R), iz(z)])
        return map_coordinates(FR, c, order=1) + Lz**2/R**3, np.sign(z)*map_coordinates(Fz, c, order=1)
    rng = np.random.default_rng(seed)
    R0 = rng.uniform(6.0, 10.0, N); v0 = np.interp(R0, Rg, vc); Lz = R0*v0
    vR0 = rng.normal(0, 30.0, N); vz0 = rng.normal(0, sz0, N)
    out = {}
    for lab, FR, Fz in (("Regel", FRl, Fzl), ("wirbelfrei", FRc, Fzc)):
        R, z, vR, vz = R0.copy(), np.zeros(N), vR0.copy(), vz0.copy()
        E = lambda R, z, vR, vz: 0.5*(vR**2 + vz**2) + 0.5*Lz**2/R**2 + map_coordinates(U, np.vstack([iR(R), iz(z)]), order=1)
        E0 = E(R, z, vR, vz); steps = int(T/TU/(dt/TU)); dtu = dt/TU
        aR, az = acc(FR, Fz, R, z, Lz); zr = []; tr = []
        for n in range(steps):
            vR += 0.5*dtu*aR; vz += 0.5*dtu*az; R += dtu*vR; z += dtu*vz
            aR, az = acc(FR, Fz, R, z, Lz); vR += 0.5*dtu*aR; vz += 0.5*dtu*az
            if n % 20 == 0: zr.append(np.sqrt(np.mean(z**2))); tr.append(n*dt)
        zr, tr = np.array(zr), np.array(tr)
        w0, w1 = tr < 1000, tr > T - 1000
        dE = E(R, z, vR, vz) - E0
        out[lab] = (np.mean(zr[w0]), np.mean(zr[w1]), np.mean(dE), np.std(dE))
    return vc, Rg, out

if __name__ == "__main__":
    print("Probe Kreisgeschwindigkeit (Regel) und Newton bei 8 kpc:")
    for h, sz in ((0.3, 18.0), (0.3, 30.0), (0.15, 12.0)):
        vc, Rg, out = run(h, sz)
        if h == 0.3 and sz == 18.0:
            print(f"  v_c(8 kpc) = {np.interp(8.0, Rg, vc):.1f} km/s")
        print(f"h = {h} kpc, σ_z0 = {sz} km/s")
        for lab, (z0, z1, mE, sE) in out.items():
            print(f"  {lab:10s} rms z: {z0:.4f} → {z1:.4f} kpc ({(z1/z0 - 1)*100:+.1f} %), ΔE nach 10 Gyr: Mittel {mE:+.1f}, Streuung {sE:.1f} (km/s)², σ_z0² = {sz**2:.0f}")
