# Formen der Formeln

Bilder erzeugt `python zeichnen.py` (Zahlen aus dem Kern). Achsen dimensionslos: r/Δs, s, a_N/a₀, r/r*.

| Formel | Form (Merkmale) | hängt ab von | Bild |
|---|---|---|---|
| Δs = GM/(3c²) | Gerade durch 0, Steigung 1 in M (log-log ebenso) | M | leiter.png |
| T = √(1 − 6Δs/r) | Wurzel: senkrechte Tangente bei r = 6Δs (T = 0), Asymptote 1; T/r Maximum bei 9Δs | Gerade Δs, Hyperbel 1/r | trichter.png |
| Trichter z = 2√(6Δs(r − 6Δs)) | Flamm-Paraboloid (Drehfläche der Wurzel); Kreise 6Δs (Hals), 9Δs (Ringdichte 1/T = √3), 18Δs (κ = 0) | Wurzel 1/√(1 − 6Δs/r) | trichter.png |
| 2π(1 − T) | 2π minus Wurzel: 2π bei 6Δs, bei 9Δs 2,66, bei 18Δs 1,15, danach ~3π·Δs/r | T | gauss_bonnet.png |
| q = e^{−s} | Exponentialkurve, q(0) = 1, Asymptote 0, kein Wendepunkt; s = 1: q = 0,368 | Wurzel s(a_N) | q_nu.png |
| ν = 1/(1 − e^{−s}) = 1 + n_B | Hyperbel ~1/s bei s → 0, Asymptote 1; s = 1: 1,582 | Exponential q | q_nu.png |
| a/a₀ über a_N/a₀ | Knick im log-log: Steigung 1 (Newton) → ½ (Wurzel), Asymptoten schneiden sich bei a_N = a₀; dort a = 1,582 a₀ | Hyperbel ν · Gerade a_N | regel.png |
| v(r) Punktmasse | Kepler-Wurzel r^{−1/2} mal √ν; flach → 172 km/s (von oben, Minimum 173,9 km/s bei 200 kpc), r* = 8,72 kpc (s = 1) | Kepler v_K, ν, r* | rotation.png |
| û = e^{−√((p+Γ)τ_D)} | Exp∘Wurzel: Betrag fällt bei ωτ_D ≈ 1 (p = iω), Phase −Im√, in der Ebene Spirale nach 0; p = 0: û = q | s² = Γτ_D, τ_D = 2Λ/c | rueckkehr.png, rueckkehr_reell.png |
| χ = s/(e^s − 1) | Sigmoide 1 → 0 (χ(0) = 1, χ(1) = 0,582), konvex, kein Wendepunkt | s · n_B(s) | chi.png |
| r*² = 3Δs·ℓ₀ | Gerade im log-log: r*(M) Steigung ½; r*² über 3Δs Steigung 1, Achsenabschnitt ℓ₀; ℓ_P² = 3Δs·λ̄: Steigung −1 in λ̄(3Δs) | Gerade Δs, Konstante ℓ₀ = c²/a₀ | leiter.png |
| q = e^{2πiτ} | Kreisring |q| = e^{−2π Im τ}; Spirale bei Im τ ∝ Re τ; Betrag Verlust, Phase Windung | Gerade τ | nome.png |
| Q = q/(1 − q) = ν − 1 | Spirale um 0 in der q-Ebene (Kreisinversion von q); reell Q = 5,364 bei s = 0,171 | q, ν | spirale_Q_nu.png |
| Abhängigkeit aller Formen | Knoten = Form, Pfeil = geht ein in | – | abhaengigkeit.png |

Kette: M → Δs → T → Trichter → Gauß-Bonnet; a_N → s → q → ν → Regel → Rotationskurve (mit v_K); s → χ, n_B; Λ → a₀ → τ_D → û → q; Δs · ℓ₀ → r*; τ → Nome → q → Q.
