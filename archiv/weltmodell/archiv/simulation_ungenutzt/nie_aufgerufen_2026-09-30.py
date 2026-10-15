
# ---- aus kern/herkunft_g.py (nie aufgerufen)
def symbol_g(zeichen):
    return bool(symbol_wurzeln(zeichen))

# ---- aus kern/galaxie.py (nie aufgerufen)
@formel('bahnen_meridian', 'Leapfrog in (R_zyl, z_h) mit J_z fest: a = ν(|g_N|)·g_N + J_z²/R_zyl³ (regel=False: Newton)',
        ('R_zyl', 'z_h', 'v_R', 'v_z'), ('g_N', 'J_z', 't', 'a₀'), 'L', 'hilfe', 'Baustein 07', None)
def bahnen_meridian(gR_fn, gz_fn, R0, z0, vR0, vz0, Lz, T, dt, a0=A0, regel=True):
    """Leapfrog in (R, z) mit L_z fest; gR_fn, gz_fn: Newton-Feld (m/s²) als Funktionen von (R, z) in m.
       regel=True: örtliche Regel a = ν(|g_N|)·g_N; False: Newton."""
    R, z, vR, vz = map(np.array, (R0, z0, vR0, vz0))

    def acc(R, z):
        gR, gz = gR_fn(R, z), gz_fn(R, z)
        f = nu(np.hypot(gR, gz), a0) if regel else 1.0
        return f * gR + Lz**2 / R**3, f * gz
    aR, az = acc(R, z)
    for _ in range(int(T / dt)):
        vR = vR + 0.5 * dt * aR; vz = vz + 0.5 * dt * az
        R = R + dt * vR; z = z + dt * vz
        aR, az = acc(R, z)
        vR = vR + 0.5 * dt * aR; vz = vz + 0.5 * dt * az
    return R, z, vR, vz
@formel('bin_mediane', 'Median je Bin in x = log₁₀(g_N/a₀) mit mindestens 20 Punkten; Fehler 1,2533·σ/√n',
        ('Ȳ_bin',), ('g_N',), 'L', 'hilfe', 'Baustein 07')
def bin_mediane(gN, y, lo, hi, n, a0=A0):
    """Array (Bins × 4): Mitte in log₁₀(g_N/a₀), Median von y, Zahl, Fehler des Medians."""
    e = np.linspace(lo, hi, n + 1)
    x = np.log10(gN / a0)
    aus = []
    for a_, b_ in zip(e[:-1], e[1:]):
        m = (x >= a_) & (x < b_) & np.isfinite(y)
        if m.sum() >= 20:
            aus.append((0.5 * (a_ + b_), np.median(y[m]), m.sum(), 1.2533 * y[m].std() / math.sqrt(m.sum())))
    return np.array(aus)

# ---- aus pruefung/zahlen_galaxie.py (nie aufgerufen)
KMSKPC = 1e3 / K.kpc

# ---- aus pruefung/gesetze_M.py (nie aufgerufen)
LEITER_MASSEN = ("m_p", "M_sun", MW, "M_HAUFEN")
