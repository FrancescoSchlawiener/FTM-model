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
import math
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
