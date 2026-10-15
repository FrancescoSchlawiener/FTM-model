# Formen der Dunklen Energie

Setzung: a₀ = 1,1·10⁻¹⁰ (Standard), n_Z = X/ln ν = 2,47, sofortige Verschiebung (M1), Mmin = 10⁶ M☉ (frei). Rechnung: `zeichnen.py` (Kern + forschung/kosmos.py), Merkmale in `merkmale.json`. Alles G-Größen der Epoche, außer Λ, a₀ (I). Rechnung, noch nicht im Modell.

| Formel | Form (Merkmale) | hängt ab von | Bild |
|---|---|---|---|
| Ω_s = (π²/45)Θ² | Parabel in Θ; Θ₀ = 0,98, Ω_s0 = 0,210 (Modell-H); über z Hyperbel, fällt mit 1/H² | Θ | formen_kurven.png (1, 2) |
| Θ = 2πa₀/(cH) | Hyperbel in H; Höhe durch Λ = c²/2a₀ (I) | H(z), a₀ | formen_kurven.png (2) |
| ρ_s = Ω_s·H² | Gerade (konstant); Θ²H² kürzt H | a₀ | formen_kurven.png (3) |
| ρ_m ∝ a⁻³ | Hyperbel in a; Gleichheit mit ρ_s bei z = −0,08 (Zukunft), mit ρ_DE bei z = 0,39 | Ω_m | formen_kurven.png (3) |
| F(a) = erfc(δ_c/(√2σD)) | Sigmoide in ln a; F(0) = 0,835, F(z=1) = 0,73, F(z=3) = 0,51; Wendepunkt bei ln a = −1,82, z = 5,17 | Mmin (Lage), D(a) aus H (Wachstum) | formen_kurven.png (4) |
| dF/dln a | Glocke, Gipfel 0,48 bei z = 5,17; Gipfel von |dF/dz| bei z = 1,39 | F | formen_kurven.png (5), wendepunkt.png |
| ρ_DE = ρ_s[1 + n_Z F/F(1)] | Sigmoide, auf die Gerade gehoben: Grundlinie ρ_s, Stufe n_Z = 2,47, Endwert 3,47·ρ_s bei z = 0 | ρ_s (Höhe: Λ, a₀), n_Z (Höhe), F (Gestalt, Lage) | formen_kurven.png (6) |
| n_Z = X/ln ν | Kac: (1−q)ⁿ = e^(−X), Exponential in n; X = 4,576, ν = 6,364, n_Z = 2,47 | q_kosm, X | formen_kurven.png (9) |
| w(z) = −1 − dlnρ_DE/(3dln a) | flache Stufe unter −1: −1,024 / −1,047 / −1,070 / −1,117 bei z = 0 / 0,5 / 1 / 2 | F | formen_kurven.png (7) |
| H(z)/H_ΛCDM | flache Beule: 1,078 (z=0), 1,014 (z=1), 1,0003 (z=2) | ρ_DE | formen_kurven.png (8) |
| Ω_DE(z) = ρ_DE/H² | fallende Sigmoide: 0,73 bei z = 0 | ρ_DE, H | formen_kurven.png (8) |
| Abhängigkeit | ρ_DE = gehobene Sigmoide (Gestalt F, Lage Mmin/D, Höhe ρ_s und n_Z) | alle | abhaengigkeit.png |

Wendepunkt: Der Wendepunkt der Sigmoide F liegt bei z = 5,17 (in ln a), das Maximum der Rate in z bei z = 1,39. Beide nicht bei z ≈ 0,4 [G]. Bei z = 0,39 liegt die Gleichheit ρ_m = ρ_DE (Schnitt Hyperbel mit gehobener Sigmoide), nicht der Wendepunkt.
