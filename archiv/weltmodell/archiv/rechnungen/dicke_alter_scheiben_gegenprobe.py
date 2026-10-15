# Gegenprobe zu dicke_alter_scheiben.py auf anderem Weg: kartesische 3D-Bahnen (x, y, z) statt Meridianebene mit L_z,
# Zeitschritt halbiert, kubische Interpolation auf einem Gitter über beide Seiten der Ebene, Energie E = ½v² + U(Φ_N) wie in der 2D-Rechnung.
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
from dicke_alter_scheiben import newton_grid, nu, TU
def run3d(h, sz0, N=150, T=10000.0, dt=0.125, seed=7):
    Rg = np.linspace(0.3, 20.0, 395); zp = np.linspace(0.0, 4.0, 401)
    ph, gR, gz = newton_grid(h, Rg, zp)
    # volles Gitter z ∈ [−4, 4]: gerade Größen gespiegelt, g_z ungerade
    zg = np.concatenate([-zp[:0:-1], zp]); ev = lambda A: np.concatenate([A[:, :0:-1], A], axis=1)
    phF, gRF, gzF = ev(ph), ev(gR), np.concatenate([-gz[:, :0:-1], gz], axis=1)
    nl = nu(np.hypot(gRF, gzF)); Pm = ph[:, 0]; nm = nu(np.abs(gR[:, 0])); o = np.argsort(Pm)
    nc = np.interp(phF, Pm[o], nm[o])
    U1 = np.concatenate([[0], np.cumsum(0.5*(nm[o][1:] + nm[o][:-1])*np.diff(Pm[o]))]); U = np.interp(phF, Pm[o], U1)
    iR = lambda R: (R - Rg[0])/(Rg[1] - Rg[0]); iz = lambda z: (z - zg[0])/(zg[1] - zg[0])
    vc = np.sqrt(np.maximum(-Rg*nl[:, 400]*gR[:, 0], 0))
    S = {k: spline_filter(v, order=3, mode="nearest") for k, v in (("gR", gRF), ("gz", gzF), ("nl", nl), ("nc", nc), ("U", U))}
    mc = lambda A, c: map_coordinates(A, c, order=3, prefilter=False, mode="nearest")
    def force(key, x, y, z):
        R = np.hypot(x, y); c = np.vstack([iR(R), iz(z)]); n = mc(S[key], c)
        fR = n*mc(S["gR"], c); return fR*x/R, fR*y/R, n*mc(S["gz"], c)
    def E(x, y, z, vx, vy, vz):
        R = np.hypot(x, y); return 0.5*(vx*vx + vy*vy + vz*vz) + mc(S["U"], np.vstack([iR(R), iz(z)]))
    rng = np.random.default_rng(seed); R0 = rng.uniform(6, 10, N); v0 = np.interp(R0, Rg, vc)
    vR0 = rng.normal(0, 30, N); vz0 = rng.normal(0, sz0, N)
    res = {}
    for lab, key in (("Regel", "nl"), ("wirbelfrei", "nc")):
        x, y, z = R0.copy(), np.zeros(N), np.zeros(N); vx, vy, vz = vR0.copy(), v0.copy(), vz0.copy()
        E0 = E(x, y, z, vx, vy, vz); dtu = dt/TU; zr = []; tr = []
        ax, ay, az = force(key, x, y, z)
        for n in range(int(T/dt)):
            vx += 0.5*dtu*ax; vy += 0.5*dtu*ay; vz += 0.5*dtu*az
            x += dtu*vx; y += dtu*vy; z += dtu*vz
            ax, ay, az = force(key, x, y, z)
            vx += 0.5*dtu*ax; vy += 0.5*dtu*ay; vz += 0.5*dtu*az
            if n % 40 == 0: zr.append(np.sqrt(np.mean(z**2))); tr.append(n*dt)
        zr, tr = np.array(zr), np.array(tr); dE = E(x, y, z, vx, vy, vz) - E0
        res[lab] = (zr[tr < 1000].mean(), zr[tr > T - 1000].mean(), dE.mean(), dE.std())
    return res
if __name__ == "__main__":
    for h, sz in ((0.3, 18.0),):
        for lab, (z0, z1, mW, sW) in run3d(h, sz).items():
            print(f"3D h = {h}, σ_z0 = {sz}: {lab:10s} rms z {z0:.4f} → {z1:.4f} kpc ({(z1/z0 - 1)*100:+.1f} %), ΔE nach 10 Gyr: Mittel {mW:+.1f}, Streuung {sW:.1f} (km/s)²")
