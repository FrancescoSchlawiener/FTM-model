# -*- coding: utf-8 -*-
"""Die Simulation: Körper bewegen sich in der Zeit.

Ein Szenario ist eine Liste von Körpern (Masse, Ladung, Ort, Geschwindigkeit) und eine Wahl der Wirkungen:
  - Neigung der Faltung (Newton-Form g = GM/r²), immer,
  - Rückkehrregel a = ν(|g_N|)·g_N (kollektive Systeme unter der Schwelle), wahlweise,
  - Takt- und Faltungsglied erster Ordnung (Perihel), wahlweise,
  - Coulomb (Windung), wenn Ladungen da sind.
Ein fester Hintergrund (Punktmasse im Ursprung) kann dazukommen, etwa die Masse einer Galaxie für Testkörper.
Integriert wird mit Runge-Kutta 4. Mitgeführt wird die Eigenzeit τ jedes Körpers im Trichter des Hintergrunds.
"""
import hashlib
import math
import os

import numpy as np

import kern as K


class Koerper:
    def __init__(self, name, masse, ort, geschw, ladung=0.0):
        self.name = name
        self.masse = float(masse)
        self.ladung = float(ladung)
        self.ort = np.array(ort, dtype=float)
        self.geschw = np.array(geschw, dtype=float)


class Welt:
    def __init__(self, koerper, regel=False, takt_glied=False, hintergrund=0.0, a0=K.a0):
        self.koerper = koerper
        self.regel = regel
        self.takt_glied = takt_glied
        self.hintergrund = hintergrund          # Masse im Ursprung, fest
        self.a0 = a0
        self.m = np.array([k.masse for k in koerper])
        self.qe = np.array([k.ladung for k in koerper])

    # ------------------------------------------------------------ Wirkungen
    def beschleunigung(self, x, v):
        n = len(self.koerper)
        gN = np.zeros_like(x)
        coul = np.zeros_like(x)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                d = x[j] - x[i]
                r = math.sqrt(d @ d)
                gN[i] += K.G * self.m[j] * d / r**3
                if self.qe[i] and self.qe[j]:
                    coul[i] -= K.coulomb_konstante() * self.qe[i] * self.qe[j] * d / (r**3 * self.m[i])
            if self.hintergrund:
                r = math.sqrt(x[i] @ x[i])
                gN[i] -= K.G * self.hintergrund * x[i] / r**3
        a = gN.copy()
        if self.regel:
            for i in range(n):
                betrag = math.sqrt(gN[i] @ gN[i])
                a[i] = K.verstaerkung(betrag, self.a0) * gN[i]
        if self.takt_glied and self.hintergrund:
            GM = K.G * self.hintergrund
            for i in range(n):
                r = math.sqrt(x[i] @ x[i])
                a[i] += GM / (K.c**2 * r**3) * ((4 * GM / r - v[i] @ v[i]) * x[i] + 4 * (x[i] @ v[i]) * v[i])
        return a + coul

    # ------------------------------------------------------------ Zeitschritte
    def lauf(self, dauer, schritte, proben=2000):
        """Integriert über die Dauer; gibt Zeiten, Orte, Geschwindigkeiten (proben Zeilen) und Eigenzeiten zurück."""
        x = np.array([k.ort for k in self.koerper])
        v = np.array([k.geschw for k in self.koerper])
        dt = dauer / schritte
        jede = max(1, schritte // proben)
        T, X, V = [0.0], [x.copy()], [v.copy()]
        tau = np.zeros(len(self.koerper))
        f = self.beschleunigung
        for s in range(1, schritte + 1):
            if self.hintergrund:
                r = np.sqrt((x * x).sum(axis=1))
                tau += dt * np.sqrt(1 - 6 * K.faltung(self.hintergrund) / r)
            k1x, k1v = v, f(x, v)
            k2x, k2v = v + 0.5 * dt * k1v, f(x + 0.5 * dt * k1x, v + 0.5 * dt * k1v)
            k3x, k3v = v + 0.5 * dt * k2v, f(x + 0.5 * dt * k2x, v + 0.5 * dt * k2v)
            k4x, k4v = v + dt * k3v, f(x + dt * k3x, v + dt * k3v)
            x = x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
            v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
            if s % jede == 0 or s == schritte:
                T.append(s * dt)
                X.append(x.copy())
                V.append(v.copy())
        return np.array(T), np.array(X), np.array(V), tau


# ---------------------------------------------------------------- viele Körper: Einfall der Mode zum Haufen
import numba as nb


@nb.njit(cache=True)
def _baum_bauen(x, m):
    """Achtelbaum über die Orte x (Barnes-Hut): je Knoten Mitte, halbe Kante, acht Kinder (−1 leer, ≥ 0 Knoten,
    −2 − p Körper p), Masse und Schwerpunkt. Kinder haben stets größere Nummern als ihre Eltern."""
    n = x.shape[0]
    gross = 3 * n + 64
    kind = -np.ones((gross, 8), np.int64)
    mitte = np.zeros((gross, 3))
    halb = np.zeros(gross)
    lo = np.empty(3)
    hi = np.empty(3)
    for d in range(3):
        lo[d] = x[:, d].min()
        hi[d] = x[:, d].max()
    for d in range(3):
        mitte[0, d] = 0.5 * (lo[d] + hi[d])
    halb[0] = 0.5 * max(hi[0] - lo[0], max(hi[1] - lo[1], hi[2] - lo[2])) * 1.0001 + 1e-30
    zahl = 1
    for p in range(n):
        k = 0
        tiefe = 0
        while True:
            o = (x[p, 0] > mitte[k, 0]) * 1 + (x[p, 1] > mitte[k, 1]) * 2 + (x[p, 2] > mitte[k, 2]) * 4
            c = kind[k, o]
            if c == -1:
                kind[k, o] = -2 - p
                break
            if c >= 0:
                k = c
                tiefe += 1
                continue
            if tiefe > 60:                                  # gleiche Orte: freier Platz im selben Knoten
                for o2 in range(8):
                    if kind[k, o2] == -1:
                        kind[k, o2] = -2 - p
                        break
                break
            q = -2 - c
            nk = zahl
            zahl += 1
            halb[nk] = 0.5 * halb[k]
            mitte[nk, 0] = mitte[k, 0] + (1 if o & 1 else -1) * halb[nk]
            mitte[nk, 1] = mitte[k, 1] + (1 if o & 2 else -1) * halb[nk]
            mitte[nk, 2] = mitte[k, 2] + (1 if o & 4 else -1) * halb[nk]
            kind[k, o] = nk
            oq = (x[q, 0] > mitte[nk, 0]) * 1 + (x[q, 1] > mitte[nk, 1]) * 2 + (x[q, 2] > mitte[nk, 2]) * 4
            kind[nk, oq] = -2 - q
            k = nk
            tiefe += 1
    masse = np.zeros(zahl)
    sp = np.zeros((zahl, 3))
    for k in range(zahl - 1, -1, -1):
        for o in range(8):
            c = kind[k, o]
            if c >= 0:
                masse[k] += masse[c]
                for d in range(3):
                    sp[k, d] += sp[c, d]
            elif c <= -2:
                p = -2 - c
                masse[k] += m[p]
                for d in range(3):
                    sp[k, d] += m[p] * x[p, d]
    return kind[:zahl], halb[:zahl], masse, sp


@nb.njit(parallel=True, cache=True)
def _baum_kraft(x, m, kind, halb, masse, sp, eps, theta, g):
    """Beschleunigung je Körper aus dem Baum: ferne Knoten als Punkt im Schwerpunkt (Kante < θ·Abstand), nahe einzeln,
    Plummer-geglättet mit eps."""
    n = x.shape[0]
    a = np.zeros((n, 3))
    e2 = eps * eps
    t2 = theta * theta
    for i in nb.prange(n):
        stapel = np.empty(1024, np.int64)
        stapel[0] = 0
        h = 1
        ax = 0.0
        ay = 0.0
        az = 0.0
        while h > 0:
            h -= 1
            k = stapel[h]
            mk = masse[k]
            dx = sp[k, 0] / mk - x[i, 0]
            dy = sp[k, 1] / mk - x[i, 1]
            dz = sp[k, 2] / mk - x[i, 2]
            r2 = dx * dx + dy * dy + dz * dz
            if 4 * halb[k] * halb[k] < t2 * r2:
                f = g * mk / (r2 + e2)**1.5
                ax += f * dx
                ay += f * dy
                az += f * dz
                continue
            for o in range(8):
                c = kind[k, o]
                if c >= 0:
                    stapel[h] = c
                    h += 1
                elif c <= -2:
                    p = -2 - c
                    if p != i:
                        dx = x[p, 0] - x[i, 0]
                        dy = x[p, 1] - x[i, 1]
                        dz = x[p, 2] - x[i, 2]
                        f = g * m[p] / (dx * dx + dy * dy + dz * dz + e2)**1.5
                        ax += f * dx
                        ay += f * dy
                        az += f * dz
        a[i, 0] = ax
        a[i, 1] = ay
        a[i, 2] = az
    return a


def schwere_viele(x, m, eps, theta=0.5):
    """Newtonsche Neigung der Faltung vieler Körper [m/s²] (Barnes-Hut, Öffnungswinkel theta, Glättung eps [m])."""
    kind, halb, masse, sp = _baum_bauen(x, m)
    return _baum_kraft(x, m, kind, halb, masse, sp, eps, theta, K.G)


def haufen_anfang(punkt, R_L, gitter=128, aussen=1.5, saat=1, a_start=0.02, kasten=4.4):
    """Anfang eines Haufens: Gaußfeld des linearen Spektrums im Punkt (Würfel der Kante kasten·aussen·R_L [Mpc],
    mitbewegt, gitter³ Zellen, Zufall saat), mit der Bedingung, dass die lineare Faltung innerhalb R_L um die Mitte heute
    Δ_L hat (Hoffman-Ribak: δ = δ_Zufall + (Δ_L − δ_Zufall,R_L(0))·ξ_×(r)/σ²(R_L), σ² und ξ_× aus den Moden des Würfels,
    so trifft die Bedingung genau); Δ_L wie bei kern.einfall_haufen, die Schale
    bei R_L erreicht ohne Kreuzung heute die Mitte. Körper auf den Gitterpunkten innerhalb aussen·R_L, verschoben nach
    Zel'dovich, Geschwindigkeit H·x + a·H·f·D·ψ. Gibt (x [m], v [m/s], Masse je Körper [kg], t_start [s], t₀ [s])."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    from scipy.integrate import solve_ivp
    rho_m, rho_r, rho_de, w, H = K.atmen_punkt(punkt)
    D, f = K.wachstum_linear(punkt)
    k_t, P_t = K.spektrum_linear(punkt)
    L = kasten * aussen * R_L
    n = gitter
    zelle = L / n
    kx = 2 * math.pi * np.fft.fftfreq(n, zelle)
    kz = 2 * math.pi * np.fft.rfftfreq(n, zelle)
    kk = np.sqrt(kx[:, None, None]**2 + kx[None, :, None]**2 + kz[None, None, :]**2)
    kk[0, 0, 0] = 1.0
    P = np.exp(np.interp(np.log(kk), np.log(k_t), np.log(P_t)))
    P[0, 0, 0] = 0.0
    zufall = np.random.default_rng(saat).standard_normal((n, n, n))
    dk = np.fft.rfftn(zufall) * np.sqrt(P / zelle**3)
    W = lambda y: 3 * (np.sin(y) - y * np.cos(y)) / y**3
    glatt0 = lambda f: float(np.fft.irfftn(f * W(kk * R_L), s=(n, n, n), axes=(0, 1, 2))[0, 0, 0])
    dR0 = glatt0(dk)
    xi_k = P * W(kk * R_L) / zelle**3                              # ⟨δ(x)·δ_R_L(0)⟩ auf dem Gitter, im Raum der Moden
    sig2 = glatt0(xi_k)
    # Δ_L aus dem Einfall der Schale R_L allein (wie kern.einfall_haufen)
    t_von_a = lambda a: quad(lambda y: 1 / H(math.exp(y)), math.log(1e-8), math.log(a), epsrel=1e-10, limit=200)[0]
    t0, ti = t_von_a(1.0), t_von_a(a_start)
    ah = solve_ivp(lambda t, y: [y[0] * H(y[0])], (ti, 1.05 * t0), [a_start], rtol=1e-11, dense_output=True)
    a_t = lambda t: float(ah.sol(t)[0])
    M_L = 4 * math.pi / 3 * rho_m * (R_L * K.Mpc)**3

    def ankunft(DL):
        d = DL * D(a_start)
        r0 = a_start * R_L * K.Mpc * (1 + d)**(-1 / 3)
        v0 = H(a_start) * r0 * (1 - f(a_start) * d / (3 * (1 + d)))
        ext = lambda t, r: -4 * math.pi * K.G / 3 * (2 * rho_r / a_t(t)**4 + rho_de(a_t(t)) * (1 + 3 * w(a_t(t)))) * r
        ev = lambda t, y: y[0] - 1e-4 * r0 / a_start
        ev.terminal = True
        s = solve_ivp(lambda t, y: [y[1], -K.G * M_L / y[0]**2 + ext(t, y[0])], (ti, 1.04 * t0), [r0, v0],
                      events=ev, rtol=1e-10)
        return s.t_events[0][0] if len(s.t_events[0]) else 2 * t0
    DL = brentq(lambda x_: ankunft(x_) - t0, 1.0, 3.0, xtol=1e-6)
    delta_k = dk + (DL - dR0) * xi_k / sig2
    psi = [np.fft.irfftn(1j * kc * delta_k / kk**2, s=(n, n, n), axes=(0, 1, 2))
           for kc in (kx[:, None, None] + 0 * kk, kx[None, :, None] + 0 * kk, kz[None, None, :] + 0 * kk)]
    q1 = (np.arange(n) * zelle + 0.5 * L) % L - 0.5 * L
    qx, qy, qz = np.meshgrid(q1, q1, q1, indexing="ij")
    innen = qx**2 + qy**2 + qz**2 < (aussen * R_L)**2
    q = np.stack([qx[innen], qy[innen], qz[innen]], axis=1)
    ps = np.stack([p_[innen] for p_ in psi], axis=1)
    Di = D(a_start)
    x = a_start * (q + Di * ps) * K.Mpc
    v = H(a_start) * x + a_start * H(a_start) * f(a_start) * Di * ps * K.Mpc
    return x, v, rho_m * (zelle * K.Mpc)**3, ti, t0


def haufen_lauf(punkt, R_L, gitter=128, aussen=1.5, saat=1, a_start=0.02, weich=10.0, eta=0.2, theta=0.5, kasten=4.4,
                melden=None):
    """Entstehung eines Haufens aus druckfreiem Staub (haufen_anfang), Leapfrog im Atmen des Punkts: Neigung der Faltung
    aller Körper (Baum), dazu glatte Strahlung und Dunkle Energie −(4πG/3)(2ρ_r + ρ_DE(1 + 3w))·x. Glättung weich [kpc],
    mitbewegt bis z = 2, danach fest. Zeitschritt η·min(√(ε/|g|)), höchstens 0,02/H. Gibt die Orte heute [m], die
    Geschwindigkeiten heute [m/s], die Masse je Körper [kg] und die Nummern der Körper, die anfangs innerhalb 0,3·R_L
    um die Mitte lagen."""
    from scipy.integrate import solve_ivp
    rho_m, rho_r, rho_de, w, H = K.atmen_punkt(punkt)
    x, v, m1, t, t0 = haufen_anfang(punkt, R_L, gitter, aussen, saat, a_start, kasten)
    start = np.flatnonzero(np.sqrt((x * x).sum(axis=1)) < 0.3 * a_start * R_L * K.Mpc)
    ah = solve_ivp(lambda t_, y: [y[0] * H(y[0])], (t, 1.05 * t0), [a_start], rtol=1e-11, dense_output=True)
    a_t = lambda t_: float(ah.sol(t_)[0])
    m = np.full(len(x), m1)

    def beschl(x, t_):
        a = a_t(t_)
        eps = weich * K.kpc * min(1.0, 3 * a)
        return schwere_viele(x, m, eps, theta) - 4 * math.pi * K.G / 3 * (
            2 * rho_r / a**4 + rho_de(a) * (1 + 3 * w(a))) * x, eps
    g, eps = beschl(x, t)
    schritt = 0
    while t < t0:
        dt = min(0.02 / H(a_t(t)), eta * float(np.sqrt(eps / np.sqrt((g * g).sum(axis=1)).max())), t0 - t)
        v += 0.5 * dt * g
        x += dt * v
        t += dt
        g, eps = beschl(x, t)
        v += 0.5 * dt * g
        schritt += 1
        if melden and schritt % melden == 0:
            print(schritt, round(t / K.Gyr, 3), round(a_t(t), 4), flush=True)
    return x, v, m, start


def haufen_mitte(x, m, kern_zahl=1000, start=None):
    """Mitte des Haufens: Schwerpunkt der Körper start (die anfangs der Mitte nächsten, der Haufen zieht mit der
    großräumigen Strömung), dann schrumpfende Kugel bis kern_zahl Körper. Gibt den Ort [m]."""
    mitte = np.zeros(3) if start is None else x[start].mean(axis=0)
    R = 3 * K.Mpc
    while True:
        d = np.sqrt(((x - mitte)**2).sum(axis=1))
        w = d < R
        if w.sum() < kern_zahl:
            return mitte
        mitte = (x[w] * m[w, None]).sum(axis=0) / m[w].sum()
        R *= 0.9


def haufen_profil(x, m, kern_zahl=1000, start=None):
    """M(<r) um die Mitte des Haufens (haufen_mitte) auf einem Gitter in r. Gibt (r [m], M [kg])."""
    mitte = haufen_mitte(x, m, kern_zahl, start)
    d = np.sort(np.sqrt(((x - mitte)**2).sum(axis=1)))
    r = np.geomspace(5 * K.kpc, 10 * K.Mpc, 300)
    return r, np.searchsorted(d, r) * m[0]


def haufen_heute(punkt, R_L, gitter=128, saat=1, neu=False, melden=None, **kw):
    """Orte und Geschwindigkeiten heute aus haufen_lauf, abgelegt in zwischen/haufen/, gerechnet nur, wenn sie fehlen
    oder neu. Gibt (Name des Laufs, x [m], v [m/s], Masse je Körper [kg], start)."""
    name = f"punkt={[round(p_, 5) for p_ in punkt]} R_L={R_L} gitter={gitter} saat={saat} {sorted(kw.items())}"
    ort = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zwischen", "haufen")
    datei = os.path.join(ort, hashlib.sha1(name.encode()).hexdigest()[:16] + ".npz")
    if neu or not os.path.exists(datei) or "v" not in np.load(datei).files:
        x, v, m, start = haufen_lauf(punkt, R_L, gitter=gitter, saat=saat, melden=melden, **kw)
        os.makedirs(ort, exist_ok=True)
        np.savez_compressed(datei, x=x.astype(np.float32), v=v.astype(np.float32), m=m[:1], start=start)
    z = np.load(datei)
    return name, z["x"].astype(float), z["v"].astype(float), float(z["m"][0]), z["start"]


def haufen_konzentration(punkt, R_L, gitter=128, saat=1, neu=False, melden=None, **kw):
    """Konzentration des Haufens aus haufen_lauf, gelesen wie die Haufen (kern.konzentration_profil). Abgelegt in
    zwischen/etappe_haufen.json (Orte und Geschwindigkeiten heute in zwischen/haufen/), gerechnet nur, wenn er fehlt
    oder neu. Gibt dict(c, M200 [M☉], R200 [kpc], rest [dex], Körper in R₂₀₀)."""
    import huelle_ablage as HA
    name = f"punkt={[round(p_, 5) for p_ in punkt]} R_L={R_L} gitter={gitter} saat={saat} {sorted(kw.items())}"
    alle = HA.lesen()
    if not neu and name in alle.get("_haufen", {}):
        return alle["_haufen"][name]
    _, x, _, m1, start = haufen_heute(punkt, R_L, gitter, saat, neu, melden, **kw)
    m = np.full(len(x), m1)
    r, M = haufen_profil(x, m, start=start)
    c, R2, M2, rest = K.konzentration_profil(r, M, punkt[2])
    erg = dict(c=c, M200=M2 / K.M_sun, R200=R2 / K.kpc, rest=rest, koerper=int(M2 / m[0]))
    alle = HA.lesen()
    alle.setdefault("_haufen", {})[name] = erg
    HA.schreiben(alle)
    return erg


def haufen_drehung(punkt, R_L, gitter=128, saat=1, linear=False, wirbel=False, **kw):
    """Verwindung des Haufens aus haufen_lauf: Drehimpuls J der Körper innerhalb R₂₀₀ um ihren Schwerpunkt, gelesen als
    Anteil der Drehung λ' = J/(√2·M₂₀₀·V₂₀₀·R₂₀₀), V₂₀₀ = √(GM₂₀₀/R₂₀₀), als Spinzahl j = a_J/(3Δs) = J·c/(G·M₂₀₀²) und
    als Tempo der Drehung gegen V₂₀₀, j·V₂₀₀/c = √2·λ'. Gibt dict(lam, j, M200 [M☉], R200 [kpc], V200 [km/s]).
    linear: dasselbe Gebiet (die Körper, die heute in R₂₀₀ liegen) am Anfang (haufen_anfang, x = a·(q + D·ψ)). Der ganze
    lineare Drehimpuls J = a²Ḋ·Σ m·q × ψ, der Kommutator J_i = a²Ḋ·ε_ijk·(T·I)_jk mit der Form I = Σ m·q·q und der Scherung T
    aus ψ ≈ −T·q (kleinste Quadrate), dazu λ' bei der Umkehr des Gebiets: seine mittlere lineare Faltung δ̄ = −3⟨ψ·q̂⟩/R_M
    auf der Lagrange-Kugel R_M seiner Masse, die Umkehr der Schale mit δ̄ im Atmen des Punkts. Gibt dict(lam, j, cos des
    Kommutators gegen den ganzen, |Kommutator|/|ganz|, a der Umkehr). wirbel: Die Körper in R₂₀₀ als Baryonen mit
    übergebener Mode (Masse m/(1 + Q)), ihr Feld nach Newton g_N (Baum, Glättung weich [kpc]) mit der Regel ν(|g_N|)·g_N,
    das Drehmoment N = Σ m·r × ν·g_N um den Schwerpunkt (nach Newton heben sich die Paare auf, die Regel ist kein
    Gradient, ihr Wirbel ∇×A dreht): Gibt dict(t_N = |J|/|N| [Gyr], cos von N gegen J, |N_Newton|/|N|)."""
    _, x, v, m1, start = haufen_heute(punkt, R_L, gitter, saat, **kw)
    m = np.full(len(x), m1)
    r, M = haufen_profil(x, m, start=start)
    _, R2, M2, _ = K.konzentration_profil(r, M, punkt[2])
    mitte = haufen_mitte(x, m, start=start)
    w = np.sqrt(((x - mitte)**2).sum(axis=1)) < R2
    M2 = m1 * w.sum()
    V2 = math.sqrt(K.G * M2 / R2)
    if wirbel:
        mb = np.full(w.sum(), m1 / (1 + K.q_je_radiant()))
        g = schwere_viele(x[w], mb, kw.get("weich", 10.0) * K.kpc, kw.get("theta", 0.5))
        nu = np.array([K.verstaerkung(b) for b in np.sqrt((g * g).sum(axis=1))])
        dx, dv = x[w] - x[w].mean(axis=0), v[w] - v[w].mean(axis=0)
        J = (mb[:, None] * np.cross(dx, dv)).sum(axis=0)
        N = (mb[:, None] * np.cross(dx, nu[:, None] * g)).sum(axis=0)
        N0 = (mb[:, None] * np.cross(dx, g)).sum(axis=0)
        return dict(t_N=float(np.linalg.norm(J) / np.linalg.norm(N)) / K.Gyr,
                    cos=float(N @ J / np.linalg.norm(N) / np.linalg.norm(J)),
                    newton=float(np.linalg.norm(N0) / np.linalg.norm(N)))
    if not linear:
        dx, dv = x[w] - x[w].mean(axis=0), v[w] - v[w].mean(axis=0)
        J = float(np.linalg.norm(np.cross(dx, dv).sum(axis=0))) * m1
        return dict(lam=J / (math.sqrt(2) * M2 * V2 * R2), j=J * K.c / (K.G * M2**2), M200=M2 / K.M_sun,
                    R200=R2 / K.kpc, V200=V2 / 1e3)
    from scipy.integrate import quad, solve_ivp
    a_s = kw.get("a_start", 0.02)
    rho_m, rho_r, rho_de, w_, H = K.atmen_punkt(punkt)
    D, f = K.wachstum_linear(punkt)
    xa, va, _, ti, t0 = haufen_anfang(punkt, R_L, gitter, kw.get("aussen", 1.5), saat, a_s, kw.get("kasten", 4.4))
    psi = (va - H(a_s) * xa) / (a_s * H(a_s) * f(a_s) * D(a_s))
    q = xa / a_s - D(a_s) * psi
    q0 = q[w].mean(axis=0)
    qi, pi = q[w] - q0, psi[w] - psi[w].mean(axis=0)
    S = np.cross(qi, pi).sum(axis=0)
    A = np.linalg.lstsq(qi, pi, rcond=None)[0]
    T, I = -(A + A.T) / 2, qi.T @ qi
    TI = T @ I
    Sk = np.array([TI[1, 2] - TI[2, 1], TI[2, 0] - TI[0, 2], TI[0, 1] - TI[1, 0]])
    R_M = (3 * M2 / (4 * math.pi * rho_m))**(1 / 3)
    dq = q - q0
    rq = np.sqrt((dq**2).sum(axis=1))
    schale = abs(rq / R_M - 1) < 0.05
    dbar = -3 * float(np.mean((psi[schale] * dq[schale]).sum(axis=1) / rq[schale])) / R_M
    ah = solve_ivp(lambda t_, y: [y[0] * H(y[0])], (ti, 1.05 * t0), [a_s], rtol=1e-11, dense_output=True)
    a_t = lambda t_: float(ah.sol(t_)[0])
    d = dbar * D(a_s)
    r0 = a_s * R_M * (1 + d)**(-1 / 3)
    v0 = H(a_s) * r0 * (1 - f(a_s) * d / (3 * (1 + d)))
    ext = lambda t_, r_: -4 * math.pi * K.G / 3 * (2 * rho_r / a_t(t_)**4 + rho_de(a_t(t_)) * (1 + 3 * w_(a_t(t_)))) * r_
    um = lambda t_, y: y[1]
    um.terminal, um.direction = True, -1
    s = solve_ivp(lambda t_, y: [y[1], -K.G * M2 / y[0]**2 + ext(t_, y[0])], (ti, 1.04 * t0), [r0, v0], events=um,
                  rtol=1e-10)
    a_u = a_t(s.t_events[0][0])
    J = float(np.linalg.norm(S)) * m1 * a_u**2 * H(a_u) * f(a_u) * D(a_u)
    return dict(lam=J / (math.sqrt(2) * M2 * V2 * R2), j=J * K.c / (K.G * M2**2),
                cos=float(S @ Sk / np.linalg.norm(S) / np.linalg.norm(Sk)),
                verhaeltnis=float(np.linalg.norm(Sk) / np.linalg.norm(S)), a_umkehr=a_u)


_GRUEN = {}


def gitter_feld(x, m, N, L, art):
    """Feld vieler Körper [m/s²] auf einem Gitter (Kante L [m], N³ Zellen, Mitte im Ursprung, Masse nach CIC, Newton
    isoliert über FFT auf 2N³), art: "newton", "vektor" (die Regel ν(|g_N|)·g_N am Gitterpunkt) oder "wirbelfrei"
    (−∇Φ mit ∇²Φ = ∇·[ν·g_N], die Quelle −∇·[ν·g_N]/(4πG) wieder isoliert gelöst). Körper außerhalb des Gitters
    lesen die ganze Masse als Punkt im Schwerpunkt, mit der Regel ν(|g_N|)·g_N."""
    h = L / N
    if (N, h) not in _GRUEN:
        i = np.fft.fftfreq(2 * N, 1.0 / (2 * N))
        r = h * np.sqrt(i[:, None, None]**2 + i[None, :, None]**2 + i[None, None, :]**2)
        r[0, 0, 0] = 0.5 * h
        _GRUEN[(N, h)] = np.fft.rfftn(-K.G / r)
    gr = _GRUEN[(N, h)]

    def potential(Mg):
        pad = np.zeros((2 * N,) * 3)
        pad[:N, :N, :N] = Mg
        return np.fft.irfftn(np.fft.rfftn(pad) * gr, s=(2 * N,) * 3)[:N, :N, :N]
    u = (x + L / 2) / h - 0.5
    i0 = np.floor(u).astype(int)
    f = u - i0
    drin = np.all((i0 >= 0) & (i0 < N - 1), axis=1)
    ecken = [(dx, dy, dz) for dx in (0, 1) for dy in (0, 1) for dz in (0, 1)]
    gew = lambda d: (np.where(d[0], f[:, 0], 1 - f[:, 0]) * np.where(d[1], f[:, 1], 1 - f[:, 1])
                     * np.where(d[2], f[:, 2], 1 - f[:, 2]))
    Mg = np.zeros((N, N, N))
    for d in ecken:
        np.add.at(Mg, (i0[drin, 0] + d[0], i0[drin, 1] + d[1], i0[drin, 2] + d[2]), m[drin] * gew(d)[drin])
    g = [-q for q in np.gradient(potential(Mg), h)]
    if art != "newton":
        b = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2) + 1e-30
        nu = 1 / (-np.expm1(-np.sqrt(b / K.a0)))
        g = [nu * q for q in g]
        if art == "wirbelfrei":
            quelle = -sum(np.gradient(g[k], h, axis=k) for k in range(3)) * h**3 / (4 * math.pi * K.G)
            g = [-q for q in np.gradient(potential(quelle), h)]
    aus = np.zeros_like(x)
    for d in ecken:
        w = gew(d)[drin]
        for k in range(3):
            aus[drin, k] += w * g[k][i0[drin, 0] + d[0], i0[drin, 1] + d[1], i0[drin, 2] + d[2]]
    if (~drin).any():
        mitte = (m[:, None] * x).sum(axis=0) / m.sum()
        dd = x[~drin] - mitte
        r = np.sqrt((dd * dd).sum(axis=1))
        gn = K.G * m.sum() / r**2
        fak = gn if art == "newton" else gn / (-np.expm1(-np.sqrt(gn / K.a0)))
        aus[~drin] = -dd / r[:, None] * fak[:, None]
    return aus


def _flaechen(schritt, rand, stufe):
    """Zellgrenzen ab 0: gleich weit im Schritt schritt, bis die Zellen um den Faktor stufe wachsen, dann wachsend bis rand."""
    gleich = np.arange(0, schritt / (stufe - 1), schritt)
    n = int(math.ceil(math.log(rand / gleich[-1]) / math.log(stufe)))
    return np.concatenate([gleich, gleich[-1] * stufe**np.arange(1, n + 1)])


def achsen_feld(scheiben, kugel, R_aus, a0s=(K.a0,), schritt=None, rand=None, stufe=1.03, je_teil=False, hoehen=None):
    """Feld in der Ebene einer achsensymmetrischen, zur Ebene spiegelgleichen Masse [m/s², nach innen], in Finiten Volumen
    auf einem Gitter in R und z ≥ 0 (Zellen im Schritt schritt = (dR, dz) [m], dann um den Faktor stufe wachsend bis rand
    [m]). scheiben: Liste (Σ(R) [kg/m²], h [m]) mit senkrechter Form e^{−|z|/h}/(2h), kugel: eingeschlossene Masse M(<r)
    [kg] oder None. Am Rand das Feld der ganzen Masse im Punkt. Gibt (g_N, [ν·g_N je a₀], [wirbelfrei je a₀]) an den
    Radien R_aus [m]; wirbelfrei ist −∇Φ mit ∇²Φ = ∇·[ν(|g_N|)·g_N], am Rand der Fluss von ν·g_N. je_teil: nur g_N in
    der Ebene, je Scheibe und Kugel einzeln (Liste). hoehen [m]: nur das Newton-Feld (g_R, g_z) in den Höhen über der Ebene,
    je ein Feld (Höhe, Radius), zwischen den Mitten der Zellen linear gelesen."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu
    R_aus = np.asarray(R_aus, dtype=float)
    hs = [h for _, h in scheiben]
    if schritt is None:
        schritt = (R_aus.min() / 2, min(hs) / 8 if hs else R_aus.min() / 2)
    rand = rand or 100 * R_aus.max()
    Rf, zf = _flaechen(schritt[0], rand, stufe), _flaechen(schritt[1], rand, stufe)
    nR, nz = len(Rf) - 1, len(zf) - 1
    Rc, zc, dz = (Rf[1:] + Rf[:-1]) / 2, (zf[1:] + zf[:-1]) / 2, np.diff(zf)
    V = (math.pi * np.diff(Rf**2))[:, None] * dz[None, :]
    teile = [S(Rc)[:, None] * ((np.exp(-zf[:-1] / h) - np.exp(-zf[1:] / h)) / (2 * dz))[None, :] for S, h in scheiben]
    if kugel is not None:
        r = np.hypot(Rc[:, None], zc[None, :])
        teile.append((kugel(r * 1.01) - kugel(r * 0.99)) / (4 * math.pi * r**2 * 0.02 * r))
    AR = 2 * math.pi * Rf[:, None] * dz[None, :]                 # Flächen der Grenzen in R (nR + 1, nz)
    Az = (math.pi * np.diff(Rf**2))[:, None] * np.ones(nz + 1)   # Flächen der Grenzen in z (nR, nz + 1)
    TR = AR[1:-1] / np.diff(Rc)[:, None]
    Tz = Az[:, 1:-1] / np.diff(zc)[None, :]
    idx = np.arange(nR * nz).reshape(nR, nz)
    zeilen, spalten, werte = [], [], []
    for a, b, T in ((idx[:-1], idx[1:], TR), (idx[:, :-1], idx[:, 1:], Tz)):
        a, b, T = a.ravel(), b.ravel(), T.ravel()
        zeilen += [a, b, a, b]; spalten += [a, b, b, a]; werte += [T, T, -T, -T]
    zeilen, spalten, werte = np.concatenate(zeilen), np.concatenate(spalten), np.concatenate(werte)
    fest = idx[-1, -1]
    w = zeilen != fest
    Mx = coo_matrix((np.append(werte[w], 1.0), (np.append(zeilen[w], fest), np.append(spalten[w], fest))),
                    shape=(nR * nz,) * 2).tocsc()
    lu = splu(Mx)

    def loesen(rechts):
        rechts = rechts.ravel().copy()
        rechts[fest] = 0.0
        return lu.solve(rechts).reshape(nR, nz)

    def flaechenfeld(Phi, bR, bz):
        gR = np.zeros((nR + 1, nz)); gz = np.zeros((nR, nz + 1))
        gR[1:-1] = -np.diff(Phi, axis=0) / np.diff(Rc)[:, None]
        gz[:, 1:-1] = -np.diff(Phi, axis=1) / np.diff(zc)[None, :]
        gR[-1], gz[:, -1] = bR[0], bz[1]
        return gR, gz

    def newton(rho):
        # Rand: Feld der Masse im Punkt, auf den Gesamtfluss −4πG·M/2 der Hälfte geeicht
        M = 2 * float((rho * V).sum())
        if M == 0:
            return (np.zeros((nR + 1, nz)), np.zeros((nR, nz + 1))), [np.zeros(nz)] * 2, [np.zeros(nR)] * 2
        rR, rz = np.hypot(Rf[-1], zc), np.hypot(Rc, zf[-1])
        bR = [-K.G * M * Rf[-1] / rR**3, -K.G * M * zc / rR**3]      # (g_R, g_z) an der äußeren Grenze in R
        bz = [-K.G * M * Rc / rz**3, -K.G * M * zf[-1] / rz**3]      # (g_R, g_z) an der oberen Grenze
        s = -2 * math.pi * K.G * M / ((AR[-1] * bR[0]).sum() + (Az[:, -1] * bz[1]).sum())
        bR, bz = [s * q for q in bR], [s * q for q in bz]
        rand_fluss = np.zeros((nR, nz))
        rand_fluss[-1] += AR[-1] * bR[0]
        rand_fluss[:, -1] += Az[:, -1] * bz[1]
        return flaechenfeld(loesen(-4 * math.pi * K.G * rho * V - rand_fluss), bR, bz), bR, bz
    ebene = lambda g: np.interp(R_aus, Rf[1:-1], -g[1:-1, 0])
    if je_teil:
        return [ebene(newton(rho)[0][0]) for rho in teile]
    (gR, gz), bR, bz = newton(sum(teile))
    cR, cz = (gR[:-1] + gR[1:]) / 2, (gz[:, :-1] + gz[:, 1:]) / 2
    if hoehen is not None:
        from scipy.interpolate import RegularGridInterpolator
        p = np.stack(np.broadcast_arrays(np.asarray(R_aus)[None, :], np.asarray(hoehen, dtype=float)[:, None]), -1)
        return tuple(-RegularGridInterpolator((Rc, zc), c, bounds_error=False, fill_value=None)(p) for c in (cR, cz))
    BR = np.hypot(gR, np.vstack([cz[:1], (cz[:-1] + cz[1:]) / 2, bR[1][None, :]]))
    Bz = np.hypot(gz, np.hstack([cR[:, :1], (cR[:, :-1] + cR[:, 1:]) / 2, bz[0][:, None]]))
    gN = ebene(gR)
    vektor, wirbelfrei = [], []
    for a0_ in a0s:
        FR, Fz = gR / -np.expm1(-np.sqrt(BR / a0_)), gz / -np.expm1(-np.sqrt(np.maximum(Bz, 1e-300) / a0_))
        q = np.zeros((nR, nz))
        q[:-1] += AR[1:-1] * FR[1:-1]; q[1:] -= AR[1:-1] * FR[1:-1]
        q[:, :-1] += Az[:, 1:-1] * Fz[:, 1:-1]; q[:, 1:] -= Az[:, 1:-1] * Fz[:, 1:-1]
        vektor.append(gN / -np.expm1(-np.sqrt(gN / a0_)))
        wirbelfrei.append(ebene(flaechenfeld(loesen(q), bR, bz)[0]))
    return gN, vektor, wirbelfrei


def polring(scheiben, kugel, r, n=48, a0_=K.a0, **kw):
    """Polarring vom Radius r [m] um eine achsensymmetrische Masse (achsen_feld): Tempo auf der Polbahn gegen die Kreisbahn
    in der Ebene, mit der Neigung nach innen über den Winkel der Bahn gemittelt, v² = r·⟨a_r⟩. Die Regel liest auf der
    Bahn den Betrag am Ort, a = ν(|g_N|)·g_N. Gibt (v_pol/v_eq nach Newton, mit der Regel)."""
    th = (np.arange(n) + 0.5) / n * math.pi / 2
    R, z = np.append(r * np.sin(th), r), np.append(r * np.cos(th), 0.0)
    gR, gz = achsen_feld(scheiben, kugel, R, hoehen=z, **kw)
    gR, gz = np.diag(gR), np.diag(gz)
    nu = 1 / -np.expm1(-np.sqrt(np.hypot(gR, gz) / a0_))
    ar = gR[:-1] * np.sin(th) + gz[:-1] * np.cos(th)
    return math.sqrt(ar.mean() / gR[-1]), math.sqrt((nu[:-1] * ar).mean() / (nu[-1] * gR[-1]))


def scheibe_lauf(n=20000, regel=True, saat=7, dauer=3.0, schritt=1.0, M_b=6e10, streuung=0.15, dicke=0.1, weich=0.1,
                 theta=0.6, proben=100, neu=False, groesse=1.0, feld="baum", gitter=128, kante=12.0):
    """Drehimpuls einer Scheibe an der Schwelle in der Zeit. Scheibe: Σ₀ = a₀/(πG), R_d = r*/√2 (M_b [M☉]), exponentiell
    bis 6·R_d, Dicke dicke·R_d (Gauß), Umlauf aus dem eigenen Feld, mit regel ν(|g_N|)·g_N (der Regel als Vektor, kein
    Gradient) oder nach Newton, Streuung streuung·v_c in allen Richtungen, Zufall saat. Leapfrog mit festem Schritt
    [Myr], Feld aus dem Baum (Glättung weich [kpc], Öffnungswinkel theta). Abgelegt in zwischen/etappe_scheibe.json
    unter "_scheibe". groesse: Scheibenlänge in Einheiten der Schwelle, R_d = groesse·r*/√2 (Σ₀ = a₀/(πG·groesse²)).
    Gibt dict(t [Gyr], jz = J_z/J_z(0), j0 = J(0)/M [kpc km/s], rh = Radius der halben Masse in der Ebene in Einheiten
    von r*, a2 = Stärke des Balkens |Σe^{2iφ}|/N innerhalb 2·R_d). feld: "baum" (Baum, die Regel als Vektor am Körper)
    oder auf dem Gitter (gitter_feld, gitter³ Zellen, Kante kante·r*) "vektor" oder "wirbelfrei", nach Newton "newton"
    gilt dort mit regel=False."""
    import huelle_ablage as HA
    name = f"n={n} regel={regel} saat={saat} dauer={dauer} schritt={schritt} M_b={M_b} streuung={streuung} dicke={dicke} weich={weich} theta={theta}"
    if groesse != 1.0:
        name += f" groesse={groesse}"
    if feld != "baum":
        name += f" feld={feld} gitter={gitter} kante={kante}"
    alle = HA.lesen()
    if not neu and name in alle.get("_scheibe", {}):
        return alle["_scheibe"][name]
    rng = np.random.default_rng(saat)
    M = M_b * K.M_sun
    rs = math.sqrt(K.G * M / K.a0)
    Rd = groesse * rs / math.sqrt(2)
    m = np.full(n, M / n)
    R = rng.gamma(2.0, Rd, n)
    R = np.where(R > 6 * Rd, rng.uniform(0, 6 * Rd, n), R)
    ph = rng.uniform(0, 2 * math.pi, n)
    x = np.stack([R * np.cos(ph), R * np.sin(ph), rng.normal(0, dicke * Rd, n)], axis=1)

    def beschl(x_):
        if feld != "baum":
            return gitter_feld(x_ - x_.mean(axis=0), m, gitter, kante * rs, feld if regel else "newton")
        g = schwere_viele(x_, m, weich * K.kpc, theta)
        if regel:
            b = np.sqrt((g * g).sum(axis=1))
            return g / (-np.expm1(-np.sqrt(b / K.a0)))[:, None]
        return g
    a = beschl(x)
    eR = np.stack([np.cos(ph), np.sin(ph), np.zeros(n)], axis=1)
    eP = np.stack([-np.sin(ph), np.cos(ph), np.zeros(n)], axis=1)
    vc = np.sqrt(np.maximum(-(a * eR).sum(axis=1) * R, 0))
    s = streuung * vc
    v = eP * np.sqrt(np.maximum(vc**2 - s**2, 0))[:, None] + rng.normal(0, 1, (n, 3)) * s[:, None]
    v -= v.mean(axis=0)
    dt = schritt * 1e-3 * K.Gyr
    J = lambda: (m[:, None] * np.cross(x - x.mean(axis=0), v - v.mean(axis=0))).sum(axis=0)
    J0 = J()

    def form():
        d = x - x.mean(axis=0)
        r = np.hypot(d[:, 0], d[:, 1])
        w = r < 2 * Rd
        return float(np.median(r) / rs), float(abs(np.exp(2j * np.arctan2(d[w, 1], d[w, 0])).mean()))
    f_ = form()
    t, jz, rh, a2 = [0.0], [1.0], [f_[0]], [f_[1]]
    for st in range(1, int(round(dauer * 1e3 / schritt)) + 1):
        v += 0.5 * dt * a
        x += dt * v
        a = beschl(x)
        v += 0.5 * dt * a
        if st % proben == 0:
            t.append(st * dt / K.Gyr)
            jz.append(float(J()[2] / J0[2]))
            f_ = form()
            rh.append(f_[0])
            a2.append(f_[1])
    erg = dict(t=t, jz=jz, j0=float(np.linalg.norm(J0) / M / (K.kpc * 1e3)), rh=rh, a2=a2)
    alle = HA.lesen()
    alle.setdefault("_scheibe", {})[name] = erg
    HA.schreiben(alle)
    return erg
