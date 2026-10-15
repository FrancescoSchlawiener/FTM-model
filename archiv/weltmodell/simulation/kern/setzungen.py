# -*- coding: utf-8 -*-
"""Setzungen des Modells (Klasse und Regeln in setzung.py). Neue Setzungen kommen hierher."""
from .setzung import SETZUNGEN

REGISTERMODUL = True
_s = SETZUNGEN.neu
_D03 = "03-von-der-galaxie-bis-zum-kosmos.md"

_s("A0", "a₀ hat den Standardwert 1,1·10⁻¹⁰ m/s² (SPARC mit realistischen Fehlern, Zelle ohne Kreispaare ≤ 1,13)",
   ("00-bauplan.md", "das Modell rechnet mit 1,1·10⁻¹⁰"), 1.1e-10, "m/s²", "a₀",
   weitere_orte=((_D03, "Schwelle a₀ = 1,1·10⁻¹⁰"), (_D03, "mit a₀ = 1,1·10⁻¹⁰"), (_D03, "(1,1·10⁻¹⁰:")))
_s("G_ERDE", "Die Neigung an der Erdoberfläche in der Formel (A1) ist g⊕ = 9,81 m/s²",
   ("01-grundlagen.md", "g⊕ = 9,81 m/s²"), 9.81, "m/s²", "g⊕")

# Übergabe der Mode ξ in gebundenen Bereichen: 1/λ_m = N*/ω_orb, N* in Radiant der Bahn (Spanne des Modells)
_NS = ("07-feldgleichung.md", "N* ≈ 0,3–2,4")
_s("NSTERN_MIN", "Die Übergabe der Mode dauert N* ≈ 0,3–2,4 Radiant der Bahn: untere Grenze", (_D03, "N* ≈ 0,3–2,4"),
   0.3, "rad", "N*", weitere_orte=(_NS,))
_s("NSTERN_MAX", "Die Übergabe der Mode dauert N* ≈ 0,3–2,4 Radiant der Bahn: obere Grenze", (_D03, "N* ≈ 0,3–2,4"),
   2.4, "rad", "N*", weitere_orte=(_NS,))

# Festlegungen des Datensatzes SPARC, wie sie alle Rechnungen an SPARC benutzen (Auswahl, Bulge, Helium, Fehlermodell)
_SP = ("03-von-der-galaxie-bis-zum-kosmos.md", "SPARC")
_s("SPARC_Q_MAX", "Auswahl der SPARC-Kurven mit Güte Q ≤ 2", _SP, 2, "1", "Q_fl")
_s("SPARC_INC_MIN", "Auswahl der SPARC-Galaxien mit Neigung ≥ 30°", _SP, 30.0, "°", "i_inc")
_s("SPARC_UPS_B", "Masse-zu-Licht des Bulges Υ_b = 0,7", _SP, 0.7, "M☉/L☉[3,6]", "Υ_b")
_s("SPARC_HELIUM", "Gasmasse 1,33·M_HI (Helium)", _SP, 1.33, "1", "f_He")
_s("SPARC_UPS_STREUUNG", "Streuung von Υ je Galaxie 0,1 dex im Fehlermodell", _SP, 0.1, "1", "σ_Υ")
