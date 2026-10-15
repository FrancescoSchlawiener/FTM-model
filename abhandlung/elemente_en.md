# The Elements: Measured, Computed, Possible

This file places all known elements with their measured data next to the pure calculation of the model and derives from it which elements are possible beyond the known ones. Every number is read from the data or computed with the core of the simulation; the file is written by the simulation and compared with it on every run.

**Data.** PubChem: 118 elements with mass, electron configuration, ionization energy, atomic radius, electronegativity, melting and boiling point, density. AME2020: 3557 nuclei, 2549 of them with measured mass; measured up to Z = 110, above that the masses are estimated.

---

## 1 What Is Computed

The model computes the elements on two levels, shell and nucleus, each with a closing form and a coupling that opens it.

**Shell.** The pure 1/r² slope closes after 2n² per shell. The inner shells open this closure: orbits with small winding dip deep and bind more tightly, the shells fill by n + l, for equal n + l the one with smaller n first. The coupling of the shell is Zα: the speed of the innermost electron relative to c. At Zα = 1 it reaches c, the meeting point of the shell, at Z = 1/α = 137.036, for a point nucleus.

**Nucleus.** The nucleus is a held funnel. Its binding follows the liquid-drop formula, fitted once to all 2509 measured nuclei of AME2020 from A = 12 (estimated masses excluded)

```
B = a_V·A − a_S·A^{2/3} − a_C·Z(Z − 1)/A^{1/3} − a_A·(N − Z)²/A ± a_P/√A
```

with a_V = 15.36, a_S = 16.69, a_C = 0.6916, a_A = 22.36, a_P = 11.67 MeV, scatter 3.39 MeV. The Coulomb term, the winding of the protons, gives the holder length λ_H = (3/5)·αħc/a_C = 1.25 fm. The fissility x = E_C/(2E_S) = a_C·Z²/(2a_S·A) compares winding and holding; at x = 1 the restoring of a deformation vanishes, the tipping point of the nucleus.

**Valley of stability.** For fixed A the nucleus with the smallest fold energy Z·E_H + N·E_n − B is the most stable, with E = 3c⁴·Δs/G for the hydrogen atom and the neutron. Against the lightest measured isobar per A (A = 20 to 238, 219 mass numbers) the valley of the liquid-drop formula, taking for each A the integer Z of smallest fold energy, is accurate to 0.32 protons on average, 76 % exact. In the valley the tipping point x = 1 lies at Z = 133.7, A = 370.3; by a second route, taking for each A the integer Z with the smallest fold energy, at Z = 134, A = 368 to 372.

**Shell term.** The liquid-drop formula knows no closures. With the closure numbers of the model (2, 8, 20, 28, 50, 82, 126, 184, 258) as a shell term, in the form of Myers and Swiatecki, with Wigner and surface-asymmetry terms, likewise fitted once, the formula meets the data better: the binding with 1.72 instead of 3.39 MeV scatter, the lightest measured isobar per A to 0.26 instead of 0.32 protons (79 instead of 76 % exact). Its surface term is larger, a_S = 22.52 MeV, a_C = 0.7379 MeV; the tipping point in the valley thus lies at Z = 181.3, A = 538.3 (integer Z = 180). The tipping point itself cannot be tested on the data: masses are measured up to Z = 110, fissions up to x = 0.942 (x of the liquid-drop formula). The tables therefore give both formulas.

---

## 2 The Ladder of the Periods

The length of period p is the number of states with n + l = p + 1, 2⌈(p + 1)/2⌉². Every step 2k² appears twice, because n + l = 2k − 1 and 2k have the same largest winding l = k − 1.

| Period | measured closure (noble gas) | measured length | 2⌈(p + 1)/2⌉² | shells by n + l |
|---|---|---|---|---|
| 1 | 2 (He) | 2 | 2 | 1s |
| 2 | 10 (Ne) | 8 | 8 | 2s 2p |
| 3 | 18 (Ar) | 8 | 8 | 3s 3p |
| 4 | 36 (Kr) | 18 | 18 | 4s 3d 4p |
| 5 | 54 (Xe) | 18 | 18 | 5s 4d 5p |
| 6 | 86 (Rn) | 32 | 32 | 6s 4f 5d 6p |
| 7 | 118 (Og) | 32 | 32 | 7s 5f 6d 7p |
| 8 | (168, computed) | – | 50 | 8s 5g 6f 7d 8p |

**The ladder as a mirror form.** A field mirrored in log r about a radius R (r ↔ R²/r), in log r v·f(1 − f) with f = 1/(eˣ + 1), binds at the threshold E = 0 exactly at the couplings N(N + 1), N = n + l, each step with 2⌈N/2⌉² states. The states of the mirror field at the threshold are spherical harmonics, u = √r·P_N^{2l+1}((r − R)/(r + R)) (exact; numerically overlap above 0.9999): degree N = n + l, order 2l + 1, poles nucleus and distance, equator R.


---

## 3 The Known Elements

### 3.1 Measured Data: Atom and Substance (PubChem)

Mass in u, ionization energy IE in eV, atomic radius r in pm, electronegativity EN (Pauling), melting and boiling point in K, density ρ in g/cm³. A dash: no measured value in PubChem.

| Z | | Name | Mass | Configuration | IE | r | EN | T_melt | T_boil | ρ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | H | Hydrogen | 1.0080 | 1s1 | 13.598 | 120 | 2.2 | 13.81 | 20.28 | 0.00008988 |
| 2 | He | Helium | 4.00260 | 1s2 | 24.587 | 140 | – | 0.95 | 4.22 | 0.0001785 |
| 3 | Li | Lithium | 7.0 | [He]2s1 | 5.392 | 182 | 0.98 | 453.65 | 1615 | 0.534 |
| 4 | Be | Beryllium | 9.012183 | [He]2s2 | 9.323 | 153 | 1.57 | 1560 | 2744 | 1.85 |
| 5 | B | Boron | 10.81 | [He]2s2 2p1 | 8.298 | 192 | 2.04 | 2348 | 4273 | 2.37 |
| 6 | C | Carbon | 12.011 | [He]2s2 2p2 | 11.260 | 170 | 2.55 | 3823 | 4098 | 2.2670 |
| 7 | N | Nitrogen | 14.007 | [He] 2s2 2p3 | 14.534 | 155 | 3.04 | 63.15 | 77.36 | 0.0012506 |
| 8 | O | Oxygen | 15.999 | [He]2s2 2p4 | 13.618 | 152 | 3.44 | 54.36 | 90.2 | 0.001429 |
| 9 | F | Fluorine | 18.99840316 | [He]2s2 2p5 | 17.423 | 135 | 3.98 | 53.53 | 85.03 | 0.001696 |
| 10 | Ne | Neon | 20.180 | [He]2s2 2p6 | 21.565 | 154 | – | 24.56 | 27.07 | 0.0008999 |
| 11 | Na | Sodium | 22.9897693 | [Ne]3s1 | 5.139 | 227 | 0.93 | 370.95 | 1156 | 0.97 |
| 12 | Mg | Magnesium | 24.305 | [Ne]3s2 | 7.646 | 173 | 1.31 | 923 | 1363 | 1.74 |
| 13 | Al | Aluminum | 26.981538 | [Ne]3s2 3p1 | 5.986 | 184 | 1.61 | 933.437 | 2792 | 2.70 |
| 14 | Si | Silicon | 28.085 | [Ne]3s2 3p2 | 8.152 | 210 | 1.9 | 1687 | 3538 | 2.3296 |
| 15 | P | Phosphorus | 30.97376200 | [Ne]3s2 3p3 | 10.487 | 180 | 2.19 | 317.3 | 553.65 | 1.82 |
| 16 | S | Sulfur | 32.07 | [Ne]3s2 3p4 | 10.360 | 180 | 2.58 | 388.36 | 717.75 | 2.067 |
| 17 | Cl | Chlorine | 35.45 | [Ne]3s2 3p5 | 12.968 | 175 | 3.16 | 171.65 | 239.11 | 0.003214 |
| 18 | Ar | Argon | 39.9 | [Ne]3s2 3p6 | 15.760 | 188 | – | 83.8 | 87.3 | 0.0017837 |
| 19 | K | Potassium | 39.0983 | [Ar]4s1 | 4.341 | 275 | 0.82 | 336.53 | 1032 | 0.89 |
| 20 | Ca | Calcium | 40.08 | [Ar]4s2 | 6.113 | 231 | 1 | 1115 | 1757 | 1.54 |
| 21 | Sc | Scandium | 44.95591 | [Ar]4s2 3d1 | 6.561 | 211 | 1.36 | 1814 | 3109 | 2.99 |
| 22 | Ti | Titanium | 47.867 | [Ar]4s2 3d2 | 6.828 | 187 | 1.54 | 1941 | 3560 | 4.5 |
| 23 | V | Vanadium | 50.9415 | [Ar]4s2 3d3 | 6.746 | 179 | 1.63 | 2183 | 3680 | 6.0 |
| 24 | Cr | Chromium | 51.996 | [Ar]3d5 4s1 | 6.767 | 189 | 1.66 | 2180 | 2944 | 7.15 |
| 25 | Mn | Manganese | 54.93804 | [Ar]4s2 3d5 | 7.434 | 197 | 1.55 | 1519 | 2334 | 7.3 |
| 26 | Fe | Iron | 55.84 | [Ar]4s2 3d6 | 7.902 | 194 | 1.83 | 1811 | 3134 | 7.874 |
| 27 | Co | Cobalt | 58.93319 | [Ar]4s2 3d7 | 7.881 | 192 | 1.88 | 1768 | 3200 | 8.86 |
| 28 | Ni | Nickel | 58.693 | [Ar]4s2 3d8 | 7.640 | 163 | 1.91 | 1728 | 3186 | 8.912 |
| 29 | Cu | Copper | 63.55 | [Ar]4s1 3d10 | 7.726 | 140 | 1.9 | 1357.77 | 2835 | 8.933 |
| 30 | Zn | Zinc | 65.4 | [Ar]4s2 3d10 | 9.394 | 139 | 1.65 | 692.68 | 1180 | 7.134 |
| 31 | Ga | Gallium | 69.723 | [Ar]4s2 3d10 4p1 | 5.999 | 187 | 1.81 | 302.91 | 2477 | 5.91 |
| 32 | Ge | Germanium | 72.63 | [Ar]4s2 3d10 4p2 | 7.900 | 211 | 2.01 | 1211.4 | 3106 | 5.323 |
| 33 | As | Arsenic | 74.92159 | [Ar]4s2 3d10 4p3 | 9.815 | 185 | 2.18 | 1090 | 887 | 5.776 |
| 34 | Se | Selenium | 78.97 | [Ar]4s2 3d10 4p4 | 9.752 | 190 | 2.55 | 493.65 | 958 | 4.809 |
| 35 | Br | Bromine | 79.90 | [Ar]4s2 3d10 4p5 | 11.814 | 183 | 2.96 | 265.95 | 331.95 | 3.11 |
| 36 | Kr | Krypton | 83.80 | [Ar]4s2 3d10 4p6 | 14.000 | 202 | 3 | 115.79 | 119.93 | 0.003733 |
| 37 | Rb | Rubidium | 85.468 | [Kr]5s1 | 4.177 | 303 | 0.82 | 312.46 | 961 | 1.53 |
| 38 | Sr | Strontium | 87.62 | [Kr]5s2 | 5.695 | 249 | 0.95 | 1050 | 1655 | 2.64 |
| 39 | Y | Yttrium | 88.90584 | [Kr]5s2 4d1 | 6.217 | 219 | 1.22 | 1795 | 3618 | 4.47 |
| 40 | Zr | Zirconium | 91.22 | [Kr]5s2 4d2 | 6.634 | 186 | 1.33 | 2128 | 4682 | 6.52 |
| 41 | Nb | Niobium | 92.90637 | [Kr]5s1 4d4 | 6.759 | 207 | 1.6 | 2750 | 5017 | 8.57 |
| 42 | Mo | Molybdenum | 95.95 | [Kr]5s1 4d5 | 7.092 | 209 | 2.16 | 2896 | 4912 | 10.2 |
| 43 | Tc | Technetium | 96.90636 | [Kr]5s2 4d5 | 7.28 | 209 | 1.9 | 2430 | 4538 | 11 |
| 44 | Ru | Ruthenium | 101.1 | [Kr]5s1 4d7 | 7.361 | 207 | 2.2 | 2607 | 4423 | 12.1 |
| 45 | Rh | Rhodium | 102.9055 | [Kr]5s1 4d8 | 7.459 | 195 | 2.28 | 2237 | 3968 | 12.4 |
| 46 | Pd | Palladium | 106.42 | [Kr]4d10 | 8.337 | 202 | 2.2 | 1828.05 | 3236 | 12.0 |
| 47 | Ag | Silver | 107.868 | [Kr]5s1 4d10 | 7.576 | 172 | 1.93 | 1234.93 | 2435 | 10.501 |
| 48 | Cd | Cadmium | 112.41 | [Kr]5s2 4d10 | 8.994 | 158 | 1.69 | 594.22 | 1040 | 8.69 |
| 49 | In | Indium | 114.818 | [Kr]5s2 4d10 5p1 | 5.786 | 193 | 1.78 | 429.75 | 2345 | 7.31 |
| 50 | Sn | Tin | 118.71 | [Kr]5s2 4d10 5p2 | 7.344 | 217 | 1.96 | 505.08 | 2875 | 7.287 |
| 51 | Sb | Antimony | 121.760 | [Kr]5s2 4d10 5p3 | 8.64 | 206 | 2.05 | 903.78 | 1860 | 6.685 |
| 52 | Te | Tellurium | 127.6 | [Kr]5s2 4d10 5p4 | 9.010 | 206 | 2.1 | 722.66 | 1261 | 6.232 |
| 53 | I | Iodine | 126.9045 | [Kr]5s2 4d10 5p5 | 10.451 | 198 | 2.66 | 386.85 | 457.55 | 4.93 |
| 54 | Xe | Xenon | 131.29 | [Kr]5s2 4d10 5p6 | 12.130 | 216 | 2.6 | 161.36 | 165.03 | 0.005887 |
| 55 | Cs | Cesium | 132.9054520 | [Xe]6s1 | 3.894 | 343 | 0.79 | 301.59 | 944 | 1.93 |
| 56 | Ba | Barium | 137.33 | [Xe]6s2 | 5.212 | 268 | 0.89 | 1000 | 2170 | 3.62 |
| 57 | La | Lanthanum | 138.9055 | [Xe]6s2 5d1 | 5.577 | 240 | 1.1 | 1191 | 3737 | 6.15 |
| 58 | Ce | Cerium | 140.116 | [Xe]6s2 4f1 5d1 | 5.539 | 235 | 1.12 | 1071 | 3697 | 6.770 |
| 59 | Pr | Praseodymium | 140.90766 | [Xe]6s2 4f3 | 5.464 | 239 | 1.13 | 1204 | 3793 | 6.77 |
| 60 | Nd | Neodymium | 144.24 | [Xe]6s2 4f4 | 5.525 | 229 | 1.14 | 1294 | 3347 | 7.01 |
| 61 | Pm | Promethium | 144.91276 | [Xe]6s2 4f5 | 5.55 | 236 | – | 1315 | 3273 | 7.26 |
| 62 | Sm | Samarium | 150.4 | [Xe]6s2 4f6 | 5.644 | 229 | 1.17 | 1347 | 2067 | 7.52 |
| 63 | Eu | Europium | 151.964 | [Xe]6s2 4f7 | 5.670 | 233 | – | 1095 | 1802 | 5.24 |
| 64 | Gd | Gadolinium | 157.25 | [Xe]6s2 4f7 5d1 | 6.150 | 237 | 1.2 | 1586 | 3546 | 7.90 |
| 65 | Tb | Terbium | 158.92535 | [Xe]6s2 4f9 | 5.864 | 221 | – | 1629 | 3503 | 8.23 |
| 66 | Dy | Dysprosium | 162.500 | [Xe]6s2 4f10 | 5.939 | 229 | 1.22 | 1685 | 2840 | 8.55 |
| 67 | Ho | Holmium | 164.93033 | [Xe]6s2 4f11 | 6.022 | 216 | 1.23 | 1747 | 2973 | 8.80 |
| 68 | Er | Erbium | 167.26 | [Xe]6s2 4f12 | 6.108 | 235 | 1.24 | 1802 | 3141 | 9.07 |
| 69 | Tm | Thulium | 168.93422 | [Xe]6s2 4f13 | 6.184 | 227 | 1.25 | 1818 | 2223 | 9.32 |
| 70 | Yb | Ytterbium | 173.05 | [Xe]6s2 4f14 | 6.254 | 242 | – | 1092 | 1469 | 6.90 |
| 71 | Lu | Lutetium | 174.9667 | [Xe]6s2 4f14 5d1 | 5.426 | 221 | 1.27 | 1936 | 3675 | 9.84 |
| 72 | Hf | Hafnium | 178.49 | [Xe]6s2 4f14 5d2 | 6.825 | 212 | 1.3 | 2506 | 4876 | 13.3 |
| 73 | Ta | Tantalum | 180.9479 | [Xe]6s2 4f14 5d3 | 7.89 | 217 | 1.5 | 3290 | 5731 | 16.4 |
| 74 | W | Tungsten | 183.84 | [Xe]6s2 4f14 5d4 | 7.98 | 210 | 2.36 | 3695 | 5828 | 19.3 |
| 75 | Re | Rhenium | 186.207 | [Xe]6s2 4f14 5d5 | 7.88 | 217 | 1.9 | 3459 | 5869 | 20.8 |
| 76 | Os | Osmium | 190.2 | [Xe]6s2 4f14 5d6 | 8.7 | 216 | 2.2 | 3306 | 5285 | 22.57 |
| 77 | Ir | Iridium | 192.22 | [Xe]6s2 4f14 5d7 | 9.1 | 202 | 2.2 | 2719 | 4701 | 22.42 |
| 78 | Pt | Platinum | 195.08 | [Xe]6s1 4f14 5d9 | 9 | 209 | 2.28 | 2041.55 | 4098 | 21.46 |
| 79 | Au | Gold | 196.96657 | [Xe]6s1 4f14 5d10 | 9.226 | 166 | 2.54 | 1337.33 | 3129 | 19.282 |
| 80 | Hg | Mercury | 200.59 | [Xe]6s2 4f14 5d10 | 10.438 | 209 | 2 | 234.32 | 629.88 | 13.5336 |
| 81 | Tl | Thallium | 204.383 | [Xe]6s2 4f14 5d10 6p1 | 6.108 | 196 | 1.62 | 577 | 1746 | 11.8 |
| 82 | Pb | Lead | 207 | [Xe]6s2 4f14 5d10 6p2 | 7.417 | 202 | 2.33 | 600.61 | 2022 | 11.342 |
| 83 | Bi | Bismuth | 208.98040 | [Xe]6s2 4f14 5d10 6p3 | 7.289 | 207 | 2.02 | 544.55 | 1837 | 9.807 |
| 84 | Po | Polonium | 208.98243 | [Xe]6s2 4f14 5d10 6p4 | 8.417 | 197 | 2 | 527 | 1235 | 9.32 |
| 85 | At | Astatine | 209.98715 | [Xe]6s2 4f14 5d10 6p5 | 9.5 | 202 | 2.2 | 575 | – | 7 |
| 86 | Rn | Radon | 222.01758 | [Xe]6s2 4f14 5d10 6p6 | 10.745 | 220 | – | 202 | 211.45 | 0.00973 |
| 87 | Fr | Francium | 223.01973 | [Rn]7s1 | 3.9 | 348 | 0.7 | 300 | – | – |
| 88 | Ra | Radium | 226.02541 | [Rn]7s2 | 5.279 | 283 | 0.9 | 973 | 1413 | 5 |
| 89 | Ac | Actinium | 227.02775 | [Rn]7s2 6d1 | 5.17 | 260 | 1.1 | 1324 | 3471 | 10.07 |
| 90 | Th | Thorium | 232.038 | [Rn]7s2 6d2 | 6.08 | 237 | 1.3 | 2023 | 5061 | 11.72 |
| 91 | Pa | Protactinium | 231.03588 | [Rn]7s2 5f2 6d1 | 5.89 | 243 | 1.5 | 1845 | – | 15.37 |
| 92 | U | Uranium | 238.0289 | [Rn]7s2 5f3 6d1 | 6.194 | 240 | 1.38 | 1408 | 4404 | 18.95 |
| 93 | Np | Neptunium | 237.048172 | [Rn]7s2 5f4 6d1 | 6.266 | 221 | 1.36 | 917 | 4175 | 20.25 |
| 94 | Pu | Plutonium | 244.06420 | [Rn]7s2 5f6 | 6.06 | 243 | 1.28 | 913 | 3501 | 19.84 |
| 95 | Am | Americium | 243.061380 | [Rn]7s2 5f7 | 5.993 | 244 | 1.3 | 1449 | 2284 | 13.69 |
| 96 | Cm | Curium | 247.07035 | [Rn]7s2 5f7 6d1 | 6.02 | 245 | 1.3 | 1618 | 3400 | 13.51 |
| 97 | Bk | Berkelium | 247.07031 | [Rn]7s2 5f9 | 6.23 | 244 | 1.3 | 1323 | – | 14 |
| 98 | Cf | Californium | 251.07959 | [Rn]7s2 5f10 | 6.30 | 245 | 1.3 | 1173 | – | – |
| 99 | Es | Einsteinium | 252.0830 | [Rn]7s2 5f11 | 6.42 | 245 | 1.3 | 1133 | – | – |
| 100 | Fm | Fermium | 257.09511 | [Rn] 5f12 7s2 | 6.50 | – | 1.3 | 1800 | – | – |
| 101 | Md | Mendelevium | 258.09843 | [Rn]7s2 5f13 | 6.58 | – | 1.3 | 1100 | – | – |
| 102 | No | Nobelium | 259.10100 | [Rn]7s2 5f14 | 6.65 | – | 1.3 | 1100 | – | – |
| 103 | Lr | Lawrencium | 266.120 | [Rn]7s2 5f14 6d1 | – | – | 1.3 | 1900 | – | – |
| 104 | Rf | Rutherfordium | 267.122 | [Rn]7s2 5f14 6d2 | – | – | – | – | – | – |
| 105 | Db | Dubnium | 268.126 | [Rn]7s2 5f14 6d3 | – | – | – | – | – | – |
| 106 | Sg | Seaborgium | 269.128 | [Rn]7s2 5f14 6d4 | – | – | – | – | – | – |
| 107 | Bh | Bohrium | 270.133 | [Rn]7s2 5f14 6d5 | – | – | – | – | – | – |
| 108 | Hs | Hassium | 269.1336 | [Rn]7s2 5f14 6d6 | – | – | – | – | – | – |
| 109 | Mt | Meitnerium | 277.154 | [Rn]7s2 5f14 6d7 (calculated) | – | – | – | – | – | – |
| 110 | Ds | Darmstadtium | 282.166 | [Rn]7s2 5f14 6d8 (predicted) | – | – | – | – | – | – |
| 111 | Rg | Roentgenium | 282.169 | [Rn]7s2 5f14 6d9 (predicted) | – | – | – | – | – | – |
| 112 | Cn | Copernicium | 286.179 | [Rn]7s2 5f14 6d10 (predicted) | – | – | – | – | – | – |
| 113 | Nh | Nihonium | 286.182 | [Rn]5f14 6d10 7s2 7p1 (predicted) | – | – | – | – | – | – |
| 114 | Fl | Flerovium | 290.192 | [Rn]7s2 7p2 5f14 6d10 (predicted) | – | – | – | – | – | – |
| 115 | Mc | Moscovium | 290.196 | [Rn]7s2 7p3 5f14 6d10 (predicted) | – | – | – | – | – | – |
| 116 | Lv | Livermorium | 293.205 | [Rn]7s2 7p4 5f14 6d10 (predicted) | – | – | – | – | – | – |
| 117 | Ts | Tennessine | 294.211 | [Rn]7s2 7p5 5f14 6d10 (predicted) | – | – | – | – | – | – |
| 118 | Og | Oganesson | 295.216 | [Rn]7s2 7p6 5f14 6d10 (predicted) | – | – | – | – | – | – |

### 3.2 Measured Data: Nucleus (AME2020)

Number of nuclei (measured of total), range of measured mass numbers, most tightly bound measured nucleus with its binding per nucleon B/A in keV and its N/Z.

| Z | | nuclei | measured A | most bound nucleus | B/A | N/Z |
|---|---|---|---|---|---|---|
| 1 | H | 6 of 7 | 1–6 | H-3 | 2827.3 | 2.000 |
| 2 | He | 8 of 8 | 3–10 | He-4 | 7073.9 | 1.000 |
| 3 | Li | 10 of 11 | 4–13 | Li-7 | 5606.4 | 1.333 |
| 4 | Be | 11 of 12 | 6–16 | Be-8 | 7062.4 | 1.000 |
| 5 | B | 15 of 16 | 7–21 | B-11 | 6927.7 | 1.200 |
| 6 | C | 14 of 16 | 8–22 | C-12 | 7680.1 | 1.000 |
| 7 | N | 14 of 16 | 10–23 | N-15 | 7699.5 | 1.143 |
| 8 | O | 16 of 18 | 11–26 | O-16 | 7976.2 | 1.000 |
| 9 | F | 16 of 19 | 14–29 | F-19 | 7779.0 | 1.111 |
| 10 | Ne | 17 of 20 | 15–31 | Ne-22 | 8080.5 | 1.200 |
| 11 | Na | 18 of 23 | 17–34 | Na-23 | 8111.5 | 1.091 |
| 12 | Mg | 19 of 23 | 19–37 | Mg-26 | 8333.9 | 1.167 |
| 13 | Al | 15 of 23 | 23–37 | Al-29 | 8348.5 | 1.231 |
| 14 | Si | 17 of 24 | 24–40 | Si-30 | 8520.7 | 1.143 |
| 15 | P | 16 of 24 | 27–42 | P-33 | 8513.8 | 1.200 |
| 16 | S | 17 of 24 | 28–44 | S-34 | 8583.5 | 1.125 |
| 17 | Cl | 17 of 25 | 30–46 | Cl-37 | 8570.3 | 1.176 |
| 18 | Ar | 17 of 26 | 32–48 | Ar-38 | 8614.3 | 1.111 |
| 19 | K | 19 of 29 | 35–53 | K-43 | 8576.2 | 1.263 |
| 20 | Ca | 21 of 29 | 36–56 | Ca-46 | 8669.0 | 1.300 |
| 21 | Sc | 21 of 29 | 39–59 | Sc-49 | 8686.3 | 1.333 |
| 22 | Ti | 20 of 29 | 40–60 | Ti-50 | 8755.7 | 1.273 |
| 23 | V | 21 of 29 | 43–63 | V-51 | 8742.1 | 1.217 |
| 24 | Cr | 21 of 30 | 44–64 | Cr-54 | 8778.0 | 1.250 |
| 25 | Mn | 21 of 31 | 46–66 | Mn-55 | 8765.0 | 1.200 |
| 26 | Fe | 20 of 32 | 48–67 | Fe-58 | 8792.3 | 1.231 |
| 27 | Co | 22 of 32 | 50–71 | Co-59 | 8768.0 | 1.185 |
| 28 | Ni | 22 of 35 | 52–73 | Ni-62 | 8794.6 | 1.214 |
| 29 | Cu | 25 of 33 | 55–79 | Cu-65 | 8757.1 | 1.241 |
| 30 | Zn | 25 of 33 | 58–82 | Zn-66 | 8759.6 | 1.200 |
| 31 | Ga | 25 of 33 | 61–85 | Ga-69 | 8724.6 | 1.226 |
| 32 | Ge | 24 of 33 | 63–86 | Ge-72 | 8731.7 | 1.250 |
| 33 | As | 23 of 33 | 65–87 | As-75 | 8700.9 | 1.273 |
| 34 | Se | 25 of 33 | 67–91 | Se-78 | 8717.8 | 1.294 |
| 35 | Br | 25 of 34 | 69–93 | Br-81 | 8695.9 | 1.314 |
| 36 | Kr | 27 of 35 | 71–97 | Kr-84 | 8717.4 | 1.333 |
| 37 | Rb | 30 of 34 | 73–102 | Rb-87 | 8711.0 | 1.351 |
| 38 | Sr | 28 of 35 | 75–102 | Sr-88 | 8732.6 | 1.316 |
| 39 | Y | 25 of 35 | 79–103 | Y-89 | 8714.0 | 1.282 |
| 40 | Zr | 25 of 37 | 81–105 | Zr-90 | 8710.0 | 1.250 |
| 41 | Nb | 28 of 38 | 83–110 | Nb-91 | 8670.9 | 1.220 |
| 42 | Mo | 27 of 39 | 85–111 | Mo-94 | 8662.3 | 1.238 |
| 43 | Tc | 28 of 40 | 87–114 | Tc-97 | 8623.7 | 1.256 |
| 44 | Ru | 29 of 41 | 89–117 | Ru-98 | 8620.3 | 1.227 |
| 45 | Rh | 29 of 41 | 92–121 | Rh-101 | 8588.2 | 1.244 |
| 46 | Pd | 32 of 42 | 92–123 | Pd-104 | 8584.8 | 1.261 |
| 47 | Ag | 30 of 42 | 96–125 | Ag-107 | 8553.9 | 1.277 |
| 48 | Cd | 36 of 42 | 97–132 | Cd-110 | 8551.3 | 1.292 |
| 49 | In | 33 of 42 | 100–132 | In-113 | 8522.9 | 1.306 |
| 50 | Sn | 35 of 42 | 100–135 | Sn-116 | 8523.1 | 1.320 |
| 51 | Sb | 33 of 41 | 105–137 | Sb-119 | 8487.9 | 1.333 |
| 52 | Te | 36 of 42 | 104–140 | Te-122 | 8478.1 | 1.346 |
| 53 | I | 34 of 42 | 109–142 | I-125 | 8450.3 | 1.358 |
| 54 | Xe | 38 of 43 | 108–146 | Xe-126 | 8443.5 | 1.333 |
| 55 | Cs | 34 of 42 | 113–148 | Cs-129 | 8416.0 | 1.345 |
| 56 | Ba | 34 of 42 | 114–150 | Ba-132 | 8409.4 | 1.357 |
| 57 | La | 28 of 42 | 124–151 | La-135 | 8382.8 | 1.368 |
| 58 | Ce | 26 of 41 | 126–151 | Ce-138 | 8377.0 | 1.379 |
| 59 | Pr | 30 of 41 | 128–157 | Pr-141 | 8354.0 | 1.390 |
| 60 | Nd | 31 of 40 | 130–160 | Nd-142 | 8346.0 | 1.367 |
| 61 | Pm | 29 of 40 | 133–161 | Pm-143 | 8317.7 | 1.344 |
| 62 | Sm | 30 of 41 | 135–164 | Sm-144 | 8303.7 | 1.323 |
| 63 | Eu | 29 of 41 | 137–165 | Eu-145 | 8269.3 | 1.302 |
| 64 | Gd | 28 of 40 | 140–167 | Gd-146 | 8249.5 | 1.281 |
| 65 | Tb | 29 of 40 | 140–168 | Tb-149 | 8209.8 | 1.292 |
| 66 | Dy | 27 of 39 | 143–169 | Dy-154 | 8193.1 | 1.333 |
| 67 | Ho | 28 of 39 | 144–171 | Ho-159 | 8165.1 | 1.373 |
| 68 | Er | 27 of 39 | 146–172 | Er-162 | 8152.4 | 1.382 |
| 69 | Tm | 28 of 39 | 147–176 | Tm-165 | 8125.5 | 1.391 |
| 70 | Yb | 27 of 38 | 151–178 | Yb-166 | 8112.5 | 1.371 |
| 71 | Lu | 29 of 39 | 153–183 | Lu-169 | 8086.3 | 1.380 |
| 72 | Hf | 30 of 38 | 156–186 | Hf-172 | 8071.7 | 1.389 |
| 73 | Ta | 29 of 40 | 157–187 | Ta-175 | 8044.4 | 1.397 |
| 74 | W | 30 of 41 | 160–191 | W-176 | 8030.1 | 1.378 |
| 75 | Re | 31 of 41 | 161–193 | Re-181 | 8004.1 | 1.413 |
| 76 | Os | 32 of 43 | 164–196 | Os-182 | 7989.7 | 1.395 |
| 77 | Ir | 31 of 43 | 167–199 | Ir-185 | 7963.7 | 1.403 |
| 78 | Pt | 34 of 44 | 168–202 | Pt-188 | 7947.9 | 1.410 |
| 79 | Au | 32 of 43 | 171–203 | Au-191 | 7924.7 | 1.418 |
| 80 | Hg | 36 of 47 | 172–208 | Hg-194 | 7914.6 | 1.425 |
| 81 | Tl | 36 of 43 | 176–213 | Tl-199 | 7893.9 | 1.457 |
| 82 | Pb | 38 of 43 | 178–215 | Pb-202 | 7882.2 | 1.463 |
| 83 | Bi | 33 of 41 | 186–218 | Bi-203 | 7857.5 | 1.446 |
| 84 | Po | 37 of 42 | 186–222 | Po-206 | 7840.6 | 1.452 |
| 85 | At | 34 of 39 | 191–224 | At-209 | 7814.8 | 1.459 |
| 86 | Rn | 37 of 39 | 193–229 | Rn-210 | 7796.7 | 1.442 |
| 87 | Fr | 37 of 37 | 197–233 | Fr-213 | 7768.5 | 1.448 |
| 88 | Ra | 34 of 35 | 201–234 | Ra-214 | 7749.2 | 1.432 |
| 89 | Ac | 32 of 33 | 205–236 | Ac-215 | 7719.4 | 1.416 |
| 90 | Th | 29 of 32 | 208–237 | Th-216 | 7697.7 | 1.400 |
| 91 | Pa | 28 of 31 | 211–238 | Pa-217 | 7664.6 | 1.385 |
| 92 | U | 24 of 29 | 215–240 | U-218 | 7640.7 | 1.370 |
| 93 | Np | 21 of 27 | 219–242 | Np-223 | 7607.6 | 1.398 |
| 94 | Pu | 19 of 27 | 228–246 | Pu-230 | 7591.0 | 1.447 |
| 95 | Am | 10 of 27 | 229–245 | Am-235 | 7565.2 | 1.474 |
| 96 | Cm | 18 of 22 | 233–251 | Cm-234 | 7550.7 | 1.438 |
| 97 | Bk | 9 of 22 | 243–251 | Bk-243 | 7517.5 | 1.505 |
| 98 | Cf | 14 of 20 | 237–254 | Cf-240 | 7510.2 | 1.449 |
| 99 | Es | 7 of 20 | 246–255 | Es-246 | 7480.8 | 1.485 |
| 100 | Fm | 11 of 20 | 246–257 | Fm-246 | 7468.0 | 1.460 |
| 101 | Md | 7 of 19 | 249–258 | Md-249 | 7444.0 | 1.465 |
| 102 | No | 7 of 17 | 252–259 | No-252 | 7425.8 | 1.471 |
| 103 | Lr | 4 of 16 | 253–256 | Lr-253 | 7402.9 | 1.456 |
| 104 | Rf | 4 of 16 | 256–261 | Rf-256 | 7385.4 | 1.462 |
| 105 | Db | 3 of 16 | 257–259 | Db-257 | 7362.0 | 1.448 |
| 106 | Sg | 3 of 16 | 260–262 | Sg-260 | 7342.6 | 1.453 |
| 107 | Bh | 2 of 19 | 261–262 | Bh-261 | 7317.3 | 1.439 |
| 108 | Hs | 3 of 18 | 264–266 | Hs-264 | 7298.4 | 1.444 |
| 109 | Mt | 1 of 18 | 266–266 | Mt-266 | 7270.8 | 1.440 |
| 110 | Ds | 2 of 18 | 269–270 | Ds-270 | 7253.8 | 1.455 |
| 111 | Rg | 0 of 15 | – | – | – | – |
| 112 | Cn | 0 of 13 | – | – | – | – |
| 113 | Nh | 0 of 13 | – | – | – | – |
| 114 | Fl | 0 of 8 | – | – | – | – |
| 115 | Mc | 0 of 6 | – | – | – | – |
| 116 | Lv | 0 of 5 | – | – | – | – |
| 117 | Ts | 0 of 4 | – | – | – | – |
| 118 | Og | 0 of 3 | – | – | – | – |

### 3.3 Calculation of the Model

Period and place, last shell by n + l, comparison with the measured configuration, coupling of the shell Zα, binding per nucleon of the most bound measured nucleus by the liquid-drop formula (measured in 3.2), valley of stability (A, N) and fissility x in the valley, by the liquid-drop formula and with the shell term (with its valley).

Magnetism of the atom: ground term by Hund from the measured configuration (largest S, then largest L per open subshell, summed; J = |L − S| below half filling of the open shell of largest winding, otherwise L + S) and moment g·J in Bohr magnetons (Landé factor, g_s = 2). Our data contain no measured moments.

| Z | | period · place | shell | n + l as measured | formation as measured | two countings as measured | ground term · moment [μ_B] | Zα | B/A computed | A_valley | N_valley | x drop | x with shell |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | H | 1 · 1 of 2 | 1s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.007 | – | 2.0 | 1.0 | 0.010 | 0.008 |
| 2 | He | 1 · 2 of 2 | 1s | yes | yes | yes | ¹S₀ · 0 | 0.015 | – | 4.0 | 2.0 | 0.021 | 0.016 |
| 3 | Li | 2 · 1 of 8 | 2s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.022 | – | 6.1 | 3.1 | 0.031 | 0.025 |
| 4 | Be | 2 · 2 of 8 | 2s | yes | yes | yes | ¹S₀ · 0 | 0.029 | – | 8.1 | 4.1 | 0.041 | 0.033 |
| 5 | B | 2 · 3 of 8 | 2p | yes | yes | yes | ²P₁/₂ · 0.33 | 0.036 | – | 10.2 | 5.2 | 0.051 | 0.041 |
| 6 | C | 2 · 4 of 8 | 2p | yes | yes | yes | ³P₀ · 0.00 | 0.044 | 7596.0 | 12.3 | 6.3 | 0.060 | 0.049 |
| 7 | N | 2 · 5 of 8 | 2p | yes | yes | yes | ⁴S₃/₂ · 3.00 | 0.051 | 7708.3 | 14.5 | 7.5 | 0.070 | 0.056 |
| 8 | O | 2 · 6 of 8 | 2p | yes | yes | yes | ³P₂ · 3.00 | 0.058 | 7958.6 | 16.6 | 8.6 | 0.080 | 0.065 |
| 9 | F | 2 · 7 of 8 | 2p | yes | yes | yes | ²P₃/₂ · 2.00 | 0.066 | 8061.4 | 18.8 | 9.8 | 0.089 | 0.071 |
| 10 | Ne | 2 · 8 of 8 | 2p | yes | yes | yes | ¹S₀ · 0 | 0.073 | 8322.3 | 20.9 | 10.9 | 0.099 | 0.079 |
| 11 | Na | 3 · 1 of 8 | 3s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.080 | 8285.9 | 23.1 | 12.1 | 0.108 | 0.086 |
| 12 | Mg | 3 · 2 of 8 | 3s | yes | yes | yes | ¹S₀ · 0 | 0.088 | 8496.7 | 25.3 | 13.3 | 0.118 | 0.094 |
| 13 | Al | 3 · 3 of 8 | 3p | yes | yes | yes | ²P₁/₂ · 0.33 | 0.095 | 8477.3 | 27.5 | 14.5 | 0.127 | 0.101 |
| 14 | Si | 3 · 4 of 8 | 3p | yes | yes | yes | ³P₀ · 0.00 | 0.102 | 8609.9 | 29.7 | 15.7 | 0.137 | 0.108 |
| 15 | P | 3 · 5 of 8 | 3p | yes | yes | yes | ⁴S₃/₂ · 3.00 | 0.109 | 8599.6 | 32.0 | 17.0 | 0.146 | 0.116 |
| 16 | S | 3 · 6 of 8 | 3p | yes | yes | yes | ³P₂ · 3.00 | 0.117 | 8682.6 | 34.2 | 18.2 | 0.155 | 0.123 |
| 17 | Cl | 3 · 7 of 8 | 3p | yes | yes | yes | ²P₃/₂ · 2.00 | 0.124 | 8678.4 | 36.5 | 19.5 | 0.164 | 0.130 |
| 18 | Ar | 3 · 8 of 8 | 3p | yes | yes | yes | ¹S₀ · 0 | 0.131 | 8726.8 | 38.8 | 20.8 | 0.173 | 0.139 |
| 19 | K | 4 · 1 of 18 | 4s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.139 | 8723.4 | 41.0 | 22.0 | 0.182 | 0.146 |
| 20 | Ca | 4 · 2 of 18 | 4s | yes | yes | yes | ¹S₀ · 0 | 0.146 | 8764.1 | 43.3 | 23.3 | 0.191 | 0.153 |
| 21 | Sc | 4 · 3 of 18 | 3d | yes | yes | yes | ²D₃/₂ · 1.20 | 0.153 | 8722.4 | 45.6 | 24.6 | 0.200 | 0.158 |
| 22 | Ti | 4 · 4 of 18 | 3d | yes | yes | yes | ³F₂ · 1.33 | 0.161 | 8805.7 | 48.0 | 26.0 | 0.209 | 0.165 |
| 23 | V | 4 · 5 of 18 | 3d | yes | yes | yes | ⁴F₃/₂ · 0.60 | 0.168 | 8793.9 | 50.3 | 27.3 | 0.218 | 0.172 |
| 24 | Cr | 4 · 6 of 18 | 3d | closes a half or full shell | yes | yes | ⁷S₃ · 6.00 | 0.175 | 8827.0 | 52.6 | 28.6 | 0.227 | 0.182 |
| 25 | Mn | 4 · 7 of 18 | 3d | yes | yes | yes | ⁶S₅/₂ · 5.00 | 0.182 | 8802.2 | 55.0 | 30.0 | 0.236 | 0.189 |
| 26 | Fe | 4 · 8 of 18 | 3d | yes | yes | yes | ⁵D₄ · 6.00 | 0.190 | 8832.8 | 57.4 | 31.4 | 0.244 | 0.195 |
| 27 | Co | 4 · 9 of 18 | 3d | yes | yes | yes | ⁴F₉/₂ · 6.00 | 0.197 | 8798.1 | 59.7 | 32.7 | 0.253 | 0.202 |
| 28 | Ni | 4 · 10 of 18 | 3d | yes | yes | yes | ³F₄ · 5.00 | 0.204 | 8826.5 | 62.1 | 34.1 | 0.262 | 0.208 |
| 29 | Cu | 4 · 11 of 18 | 3d | closes a half or full shell | yes | yes | ²S₁/₂ · 1.00 | 0.212 | 8800.5 | 64.5 | 35.5 | 0.270 | 0.213 |
| 30 | Zn | 4 · 12 of 18 | 3d | yes | yes | yes | ¹S₀ · 0 | 0.219 | 8810.7 | 66.9 | 36.9 | 0.279 | 0.219 |
| 31 | Ga | 4 · 13 of 18 | 4p | yes | yes | yes | ²P₁/₂ · 0.33 | 0.226 | 8787.6 | 69.4 | 38.4 | 0.287 | 0.226 |
| 32 | Ge | 4 · 14 of 18 | 4p | yes | yes | yes | ³P₀ · 0.00 | 0.234 | 8800.3 | 71.8 | 39.8 | 0.296 | 0.233 |
| 33 | As | 4 · 15 of 18 | 4p | yes | yes | yes | ⁴S₃/₂ · 3.00 | 0.241 | 8770.8 | 74.2 | 41.2 | 0.304 | 0.239 |
| 34 | Se | 4 · 16 of 18 | 4p | yes | yes | yes | ³P₂ · 3.00 | 0.248 | 8774.3 | 76.7 | 42.7 | 0.312 | 0.246 |
| 35 | Br | 4 · 17 of 18 | 4p | yes | yes | yes | ²P₃/₂ · 2.00 | 0.255 | 8741.4 | 79.1 | 44.1 | 0.321 | 0.252 |
| 36 | Kr | 4 · 18 of 18 | 4p | yes | yes | yes | ¹S₀ · 0 | 0.263 | 8738.6 | 81.6 | 45.6 | 0.329 | 0.259 |
| 37 | Rb | 5 · 1 of 18 | 5s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.270 | 8703.9 | 84.1 | 47.1 | 0.337 | 0.265 |
| 38 | Sr | 5 · 2 of 18 | 5s | yes | yes | yes | ¹S₀ · 0 | 0.277 | 8721.3 | 86.6 | 48.6 | 0.346 | 0.272 |
| 39 | Y | 5 · 3 of 18 | 4d | yes | yes | yes | ²D₃/₂ · 1.20 | 0.285 | 8700.3 | 89.1 | 50.1 | 0.354 | 0.280 |
| 40 | Zr | 5 · 4 of 18 | 4d | yes | yes | yes | ³F₂ · 1.33 | 0.292 | 8697.8 | 91.6 | 51.6 | 0.362 | 0.289 |
| 41 | Nb | 5 · 5 of 18 | 4d | s to d | yes | yes | ⁶D₁/₂ · 1.67 | 0.299 | 8659.1 | 94.1 | 53.1 | 0.370 | 0.295 |
| 42 | Mo | 5 · 6 of 18 | 4d | closes a half or full shell | yes | yes | ⁷S₃ · 6.00 | 0.306 | 8662.1 | 96.7 | 54.7 | 0.378 | 0.301 |
| 43 | Tc | 5 · 7 of 18 | 4d | yes | yes | yes | ⁶S₅/₂ · 5.00 | 0.314 | 8636.9 | 99.2 | 56.2 | 0.386 | 0.308 |
| 44 | Ru | 5 · 8 of 18 | 4d | s to d | yes | yes | ⁵F₅ · 7.00 | 0.321 | 8622.5 | 101.8 | 57.8 | 0.394 | 0.314 |
| 45 | Rh | 5 · 9 of 18 | 4d | s to d | yes | yes | ⁴F₉/₂ · 6.00 | 0.328 | 8599.0 | 104.3 | 59.3 | 0.402 | 0.320 |
| 46 | Pd | 5 · 10 of 18 | 4d | closes a half or full shell | yes | yes | ¹S₀ · 0 | 0.336 | 8596.5 | 106.9 | 60.9 | 0.410 | 0.326 |
| 47 | Ag | 5 · 11 of 18 | 4d | closes a half or full shell | yes | yes | ²S₁/₂ · 1.00 | 0.343 | 8570.3 | 109.5 | 62.5 | 0.418 | 0.333 |
| 48 | Cd | 5 · 12 of 18 | 4d | yes | yes | yes | ¹S₀ · 0 | 0.350 | 8563.7 | 112.1 | 64.1 | 0.426 | 0.339 |
| 49 | In | 5 · 13 of 18 | 5p | yes | yes | yes | ²P₁/₂ · 0.33 | 0.358 | 8535.6 | 114.7 | 65.7 | 0.434 | 0.345 |
| 50 | Sn | 5 · 14 of 18 | 5p | yes | yes | yes | ³P₀ · 0.00 | 0.365 | 8526.0 | 117.3 | 67.3 | 0.442 | 0.348 |
| 51 | Sb | 5 · 15 of 18 | 5p | yes | yes | yes | ⁴S₃/₂ · 3.00 | 0.372 | 8496.8 | 119.9 | 68.9 | 0.449 | 0.352 |
| 52 | Te | 5 · 16 of 18 | 5p | yes | yes | yes | ³P₂ · 3.00 | 0.379 | 8484.9 | 122.6 | 70.6 | 0.457 | 0.358 |
| 53 | I | 5 · 17 of 18 | 5p | yes | yes | yes | ²P₃/₂ · 2.00 | 0.387 | 8455.0 | 125.2 | 72.2 | 0.465 | 0.364 |
| 54 | Xe | 5 · 18 of 18 | 5p | yes | yes | yes | ¹S₀ · 0 | 0.394 | 8448.6 | 127.9 | 73.9 | 0.473 | 0.370 |
| 55 | Cs | 6 · 1 of 32 | 6s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.401 | 8419.8 | 130.5 | 75.5 | 0.480 | 0.376 |
| 56 | Ba | 6 · 2 of 32 | 6s | yes | yes | yes | ¹S₀ · 0 | 0.409 | 8406.4 | 133.2 | 77.2 | 0.488 | 0.382 |
| 57 | La | 6 · 3 of 32 | 4f | d before f | yes | yes | ²D₃/₂ · 1.20 | 0.416 | 8377.1 | 135.9 | 78.9 | 0.495 | 0.388 |
| 58 | Ce | 6 · 4 of 32 | 4f | d before f | yes | yes | ³H₄ · 3.20 | 0.423 | 8362.3 | 138.6 | 80.6 | 0.503 | 0.394 |
| 59 | Pr | 6 · 5 of 32 | 4f | yes | [Xe] 6s1 4f2 5d2 | yes | ⁴I₉/₂ · 3.27 | 0.431 | 8332.8 | 141.3 | 82.3 | 0.511 | 0.405 |
| 60 | Nd | 6 · 6 of 32 | 4f | yes | yes | yes | ⁵I₄ · 2.40 | 0.438 | 8325.6 | 144.0 | 84.0 | 0.518 | 0.411 |
| 61 | Pm | 6 · 7 of 32 | 4f | yes | yes | yes | ⁶H₅/₂ · 0.71 | 0.445 | 8300.5 | 146.7 | 85.7 | 0.526 | 0.417 |
| 62 | Sm | 6 · 8 of 32 | 4f | yes | yes | yes | ⁷F₀ · 0.00 | 0.452 | 8285.1 | 149.5 | 87.5 | 0.533 | 0.423 |
| 63 | Eu | 6 · 9 of 32 | 4f | yes | yes | yes | ⁸S₇/₂ · 7.00 | 0.460 | 8252.3 | 152.2 | 89.2 | 0.540 | 0.429 |
| 64 | Gd | 6 · 10 of 32 | 4f | closes a half or full shell | yes | yes | ⁹D₂ · 5.33 | 0.467 | 8229.2 | 155.0 | 91.0 | 0.548 | 0.434 |
| 65 | Tb | 6 · 11 of 32 | 4f | yes | [Xe] 6s2 4f8 5d1 | yes | ⁶H₁₅/₂ · 10.00 | 0.474 | 8205.3 | 157.7 | 92.7 | 0.555 | 0.440 |
| 66 | Dy | 6 · 12 of 32 | 4f | yes | [Xe] 6s2 4f9 5d1 | yes | ⁵I₈ · 10.00 | 0.482 | 8201.0 | 160.5 | 94.5 | 0.562 | 0.446 |
| 67 | Ho | 6 · 13 of 32 | 4f | yes | [Xe] 6s2 4f10 5d1 | yes | ⁴I₁₅/₂ · 9.00 | 0.489 | 8175.3 | 163.3 | 96.3 | 0.570 | 0.452 |
| 68 | Er | 6 · 14 of 32 | 4f | yes | [Xe] 6s2 4f11 5d1 | yes | ³H₆ · 7.00 | 0.496 | 8159.4 | 166.1 | 98.1 | 0.577 | 0.457 |
| 69 | Tm | 6 · 15 of 32 | 4f | yes | [Xe] 6s2 4f12 5d1 | yes | ²F₇/₂ · 4.00 | 0.504 | 8131.9 | 168.9 | 99.9 | 0.584 | 0.463 |
| 70 | Yb | 6 · 16 of 32 | 4f | yes | [Xe] 6s2 4f13 5d1 | yes | ¹S₀ · 0 | 0.511 | 8117.9 | 171.7 | 101.7 | 0.591 | 0.469 |
| 71 | Lu | 6 · 17 of 32 | 5d | yes | [Xe] 6s2 4f14 6p1 | yes | ²D₃/₂ · 1.20 | 0.518 | 8091.1 | 174.5 | 103.5 | 0.599 | 0.474 |
| 72 | Hf | 6 · 18 of 32 | 5d | yes | yes | yes | ³F₂ · 1.33 | 0.525 | 8074.7 | 177.3 | 105.3 | 0.606 | 0.480 |
| 73 | Ta | 6 · 19 of 32 | 5d | yes | yes | yes | ⁴F₃/₂ · 0.60 | 0.533 | 8047.8 | 180.2 | 107.2 | 0.613 | 0.485 |
| 74 | W | 6 · 20 of 32 | 5d | yes | [Xe] 6s1 4f14 5d5 | yes | ⁵D₀ · 0.00 | 0.540 | 8032.3 | 183.0 | 109.0 | 0.620 | 0.491 |
| 75 | Re | 6 · 21 of 32 | 5d | yes | yes | yes | ⁶S₅/₂ · 5.00 | 0.547 | 8003.9 | 185.9 | 110.9 | 0.627 | 0.496 |
| 76 | Os | 6 · 22 of 32 | 5d | yes | yes | yes | ⁵D₄ · 6.00 | 0.555 | 7989.3 | 188.8 | 112.8 | 0.634 | 0.502 |
| 77 | Ir | 6 · 23 of 32 | 5d | yes | yes | yes | ⁴F₉/₂ · 6.00 | 0.562 | 7962.9 | 191.6 | 114.6 | 0.641 | 0.507 |
| 78 | Pt | 6 · 24 of 32 | 5d | s to d | yes | yes | ³D₃ · 4.00 | 0.569 | 7945.7 | 194.5 | 116.5 | 0.648 | 0.513 |
| 79 | Au | 6 · 25 of 32 | 5d | closes a half or full shell | yes | yes | ²S₁/₂ · 1.00 | 0.576 | 7919.3 | 197.4 | 118.4 | 0.655 | 0.518 |
| 80 | Hg | 6 · 26 of 32 | 5d | yes | yes | yes | ¹S₀ · 0 | 0.584 | 7901.7 | 200.3 | 120.3 | 0.662 | 0.523 |
| 81 | Tl | 6 · 27 of 32 | 6p | yes | yes | yes | ²P₁/₂ · 0.33 | 0.591 | 7870.1 | 203.3 | 122.3 | 0.669 | 0.529 |
| 82 | Pb | 6 · 28 of 32 | 6p | yes | yes | yes | ³P₀ · 0.00 | 0.598 | 7851.7 | 206.2 | 124.2 | 0.676 | 0.532 |
| 83 | Bi | 6 · 29 of 32 | 6p | yes | yes | yes | ⁴S₃/₂ · 3.00 | 0.606 | 7831.2 | 209.1 | 126.1 | 0.683 | 0.540 |
| 84 | Po | 6 · 30 of 32 | 6p | yes | yes | yes | ³P₂ · 3.00 | 0.613 | 7813.0 | 212.1 | 128.1 | 0.690 | 0.546 |
| 85 | At | 6 · 31 of 32 | 6p | yes | yes | yes | ²P₃/₂ · 2.00 | 0.620 | 7786.9 | 215.0 | 130.0 | 0.696 | 0.551 |
| 86 | Rn | 6 · 32 of 32 | 6p | yes | yes | yes | ¹S₀ · 0 | 0.628 | 7772.7 | 218.0 | 132.0 | 0.703 | 0.556 |
| 87 | Fr | 7 · 1 of 32 | 7s | yes | yes | yes | ²S₁/₂ · 1.00 | 0.635 | 7747.0 | 221.0 | 134.0 | 0.710 | 0.562 |
| 88 | Ra | 7 · 2 of 32 | 7s | yes | yes | yes | ¹S₀ · 0 | 0.642 | 7731.0 | 223.9 | 135.9 | 0.717 | 0.567 |
| 89 | Ac | 7 · 3 of 32 | 5f | d before f | yes | yes | ²D₃/₂ · 1.20 | 0.649 | 7705.6 | 226.9 | 137.9 | 0.723 | 0.572 |
| 90 | Th | 7 · 4 of 32 | 5f | d before f | yes | yes | ³F₂ · 1.33 | 0.657 | 7685.6 | 229.9 | 139.9 | 0.730 | 0.577 |
| 91 | Pa | 7 · 5 of 32 | 5f | d before f | [Rn] 7s2 5f1 6d2 | yes | ⁴K₁₁/₂ · 4.23 | 0.664 | 7656.4 | 232.9 | 141.9 | 0.737 | 0.583 |
| 92 | U | 7 · 6 of 32 | 5f | d before f | yes | yes | ⁵L₆ · 4.29 | 0.671 | 7632.7 | 236.0 | 144.0 | 0.743 | 0.588 |
| 93 | Np | 7 · 7 of 32 | 5f | d before f | yes | yes | ⁶L₁₁/₂ · 3.38 | 0.679 | 7615.6 | 239.0 | 146.0 | 0.750 | 0.593 |
| 94 | Pu | 7 · 8 of 32 | 5f | yes | yes | yes | ⁷F₀ · 0.00 | 0.686 | 7602.5 | 242.0 | 148.0 | 0.757 | 0.598 |
| 95 | Am | 7 · 9 of 32 | 5f | yes | yes | yes | ⁸S₇/₂ · 7.00 | 0.693 | 7575.9 | 245.1 | 150.1 | 0.763 | 0.603 |
| 96 | Cm | 7 · 10 of 32 | 5f | closes a half or full shell | yes | yes | ⁹D₂ · 5.33 | 0.701 | 7559.5 | 248.1 | 152.1 | 0.770 | 0.608 |
| 97 | Bk | 7 · 11 of 32 | 5f | yes | [Rn] 7s2 5f8 6d1 | yes | ⁶H₁₅/₂ · 10.00 | 0.708 | 7528.3 | 251.2 | 154.2 | 0.776 | 0.613 |
| 98 | Cf | 7 · 12 of 32 | 5f | yes | – | yes | ⁵I₈ · 10.00 | 0.715 | 7517.1 | 254.3 | 156.3 | 0.783 | 0.618 |
| 99 | Es | 7 · 13 of 32 | 5f | yes | – | yes | ⁴I₁₅/₂ · 9.00 | 0.722 | 7488.0 | 257.4 | 158.4 | 0.789 | 0.623 |
| 100 | Fm | 7 · 14 of 32 | 5f | yes | – | yes | ³H₆ · 7.00 | 0.730 | 7474.7 | 260.5 | 160.5 | 0.796 | 0.629 |
| 101 | Md | 7 · 15 of 32 | 5f | yes | – | yes | ²F₇/₂ · 4.00 | 0.737 | 7450.5 | 263.6 | 162.6 | 0.802 | 0.634 |
| 102 | No | 7 · 16 of 32 | 5f | yes | – | yes | ¹S₀ · 0 | 0.744 | 7432.3 | 266.7 | 164.7 | 0.808 | 0.639 |
| 103 | Lr | 7 · 17 of 32 | 6d | yes | – | yes | ²D₃/₂ · 1.20 | 0.752 | 7407.8 | 269.8 | 166.8 | 0.815 | 0.643 |
| 104 | Rf | 7 · 18 of 32 | 6d | yes | – | yes | ³F₂ · 1.33 | 0.759 | 7389.7 | 273.0 | 169.0 | 0.821 | 0.648 |
| 105 | Db | 7 · 19 of 32 | 6d | yes | – | yes | ⁴F₃/₂ · 0.60 | 0.766 | 7364.2 | 276.1 | 171.1 | 0.828 | 0.653 |
| 106 | Sg | 7 · 20 of 32 | 6d | yes | – | yes | ⁵D₀ · 0.00 | 0.774 | 7346.3 | 279.3 | 173.3 | 0.834 | 0.658 |
| 107 | Bh | 7 · 21 of 32 | 6d | yes | – | yes | ⁶S₅/₂ · 5.00 | 0.781 | 7319.6 | 282.4 | 175.4 | 0.840 | 0.663 |
| 108 | Hs | 7 · 22 of 32 | 6d | yes | – | yes | ⁵D₄ · 6.00 | 0.788 | 7301.9 | 285.6 | 177.6 | 0.846 | 0.668 |
| 109 | Mt | 7 · 23 of 32 | 6d | yes | – | yes | ⁴F₉/₂ · 6.00 | 0.795 | 7274.0 | 288.8 | 179.8 | 0.853 | 0.673 |
| 110 | Ds | 7 · 24 of 32 | 6d | yes | – | yes | ³F₄ · 5.00 | 0.803 | 7260.9 | 292.0 | 182.0 | 0.859 | 0.678 |
| 111 | Rg | 7 · 25 of 32 | 6d | yes | – | yes | ²D₅/₂ · 3.00 | 0.810 | – | 295.2 | 184.2 | 0.865 | 0.684 |
| 112 | Cn | 7 · 26 of 32 | 6d | yes | – | yes | ¹S₀ · 0 | 0.817 | – | 298.4 | 186.4 | 0.871 | 0.694 |
| 113 | Nh | 7 · 27 of 32 | 7p | yes | – | yes | ²P₁/₂ · 0.33 | 0.825 | – | 301.6 | 188.6 | 0.877 | 0.700 |
| 114 | Fl | 7 · 28 of 32 | 7p | yes | – | yes | ³P₀ · 0.00 | 0.832 | – | 304.8 | 190.8 | 0.884 | 0.705 |
| 115 | Mc | 7 · 29 of 32 | 7p | yes | – | yes | ⁴S₃/₂ · 3.00 | 0.839 | – | 308.1 | 193.1 | 0.890 | 0.710 |
| 116 | Lv | 7 · 30 of 32 | 7p | yes | – | yes | ³P₂ · 3.00 | 0.846 | – | 311.3 | 195.3 | 0.896 | 0.714 |
| 117 | Ts | 7 · 31 of 32 | 7p | yes | – | yes | ²P₃/₂ · 2.00 | 0.854 | – | 314.6 | 197.6 | 0.902 | 0.719 |
| 118 | Og | 7 · 32 of 32 | 7p | yes | – | yes | ¹S₀ · 0 | 0.861 | – | 317.8 | 199.8 | 0.908 | 0.724 |

### 3.4 The 19 Deviations from n + l

99 of the 118 elements fill exactly by n + l. All deviations lie in d and f shells within the periods; none changes a closure number.

| Z | | measured (PubChem) | by n + l | kind | formation | two countings |
|---|---|---|---|---|---|---|
| 24 | Cr | [Ar]3d5 4s1 | [Ar] 4s2 3d4 | closes a half or full shell | as measured | as measured |
| 29 | Cu | [Ar]4s1 3d10 | [Ar] 4s2 3d9 | closes a half or full shell | as measured | as measured |
| 41 | Nb | [Kr]5s1 4d4 | [Kr] 5s2 4d3 | s to d | as measured | as measured |
| 42 | Mo | [Kr]5s1 4d5 | [Kr] 5s2 4d4 | closes a half or full shell | as measured | as measured |
| 44 | Ru | [Kr]5s1 4d7 | [Kr] 5s2 4d6 | s to d | as measured | as measured |
| 45 | Rh | [Kr]5s1 4d8 | [Kr] 5s2 4d7 | s to d | as measured | as measured |
| 46 | Pd | [Kr]4d10 | [Kr] 5s2 4d8 | closes a half or full shell | as measured | as measured |
| 47 | Ag | [Kr]5s1 4d10 | [Kr] 5s2 4d9 | closes a half or full shell | as measured | as measured |
| 57 | La | [Xe]6s2 5d1 | [Xe] 6s2 4f1 | d before f | as measured | as measured |
| 58 | Ce | [Xe]6s2 4f1 5d1 | [Xe] 6s2 4f2 | d before f | as measured | as measured |
| 64 | Gd | [Xe]6s2 4f7 5d1 | [Xe] 6s2 4f8 | closes a half or full shell | as measured | as measured |
| 78 | Pt | [Xe]6s1 4f14 5d9 | [Xe] 6s2 4f14 5d8 | s to d | as measured | as measured |
| 79 | Au | [Xe]6s1 4f14 5d10 | [Xe] 6s2 4f14 5d9 | closes a half or full shell | as measured | as measured |
| 89 | Ac | [Rn]7s2 6d1 | [Rn] 7s2 5f1 | d before f | as measured | as measured |
| 90 | Th | [Rn]7s2 6d2 | [Rn] 7s2 5f2 | d before f | as measured | as measured |
| 91 | Pa | [Rn]7s2 5f2 6d1 | [Rn] 7s2 5f3 | d before f | [Rn] 7s2 5f1 6d2 | as measured |
| 92 | U | [Rn]7s2 5f3 6d1 | [Rn] 7s2 5f4 | d before f | as measured | as measured |
| 93 | Np | [Rn]7s2 5f4 6d1 | [Rn] 7s2 5f5 | d before f | as measured | as measured |
| 96 | Cm | [Rn]7s2 5f7 6d1 | [Rn] 7s2 5f8 | closes a half or full shell | as measured | as measured |

8 close a half or full shell (d⁵, d¹⁰, f⁷), 7 begin the f shell with a d electron, the smaller winding, 4 pass an s electron to the d shell.

**The shell from its formation.** What is computed is how the shell forms while δ of the whole falls: the nucleus, a held funnel, holding k electrons captures the next one where the folding of the whole shell becomes deepest. Every orbit closes per revolution in the slope of the nucleus and of the other windings, at the clock rate of unwinding set by the centre. Windings of the same sense block each other's state (closure operator 𝒫), counted over the modes of the spherical harmonics. Protected is the count, not the single loop: where at least two windings lock and the counting rule closes (d⁵, d¹⁰, f⁷, full shells), it stays bound; the others move after every capture as long as the folding of the whole falls. No value is fitted; the measurement only checks at the end. Of 97 computed elements with measured configuration (PubChem) the formation meets 86, of the 19 deviations from n + l among them 18. Off are Pr, Tb, Dy, Ho, Er, Tm, Yb, Lu, W, Pa, Bk. Open are the pairing of windings of opposite sense in the same mode and the binding of f against d.

**The shell from two countings.** The coupling of the windings is counted in two ways: locally, from the density at the place, and per sense of revolution, each winding against each of the same sense. Which counting holds is set by the filling of the spin half of the changing shell of largest winding, weighted by the ratio of the tempos of the two orbits between which the electron changes (Fermi form). The distance of two configurations is mirrored about the middle, the mean of both directions of view, at the locked level J of both senses. From n + l one electron after the other moves to the lower neighbour, without knowledge of the measurement. Computed so, 118 of 118 elements with a configuration in PubChem (1–102 measured, 103–118 predicted) meet the configuration, among them 19 of the 19 deviations from n + l. Computed towards zero grid step, only Np lies narrowly, by 0.01 eV, on the side of 5f⁵7s².

### 3.5 What the Calculation Says about the Known Elements

- The noble gases (PubChem) close after 2, 10, 18, 36, 54, 86, 118; computed by n + l, exactly the same numbers.
- The formula meets the binding per nucleon of the most bound measured nuclei (A ≥ 12, 105 elements) on average to 32 keV; the deviations carry the closures, which the formula does not contain. Below A = 12 the formula is not fitted.
- The heaviest known nuclei are neutron-poor: Og-294 (mass estimated) has N = 176, the valley requires N = 200 at Z = 118. This is why Og-294 already reaches the tipping point by the liquid-drop formula, x = 0.98; in the valley, Z = 118 would lie at x = 0.91, with the shell term at 0.72 (Og-294: 0.78).
- Masses are measured up to Z = 110. The most neutron-rich measured nucleus there has N = 160, the valley requires N = 182.

---

## 4 New Elements

Beyond 118 the shell continues the eighth period: first 8s with ununennium (Uue, 119) and unbinilium (Ubn, 120), then by n + l the first g shell (l = 4), a kind of element that does not yet exist. By n + l the period has 50 places and would close at 168. Two limits bound the elements, both computed:

- **Shell:** The binding of the innermost electron reaches twice the fold energy of the electron, 6c⁴·Δs_e/G, at Z = 170 to 171 (4.3). This is robust: it hardly depends on the formula of the nucleus.
- **Nucleus:** The tipping point x = 1 in the valley lies at Z = 134 by the liquid-drop formula, and at Z = 181 to 191 with the shell term of the closure numbers (depending on the next proton closure 126, 120 or 114). The difference comes from the surface term (a_S = 16.7 versus 22.5 MeV). The shell formula meets binding and valley of the measured nuclei better (1), Q_α from Z = 100 worse (4.3); the tipping point itself is not measured (masses up to Z = 110), and both x order the measured fission lifetimes equally well.

Secured, held in the valley by both formulas, are therefore the elements Uue (119) to Utt (133). By the formula that meets the data better, with the shell term, they reach the limit of the shell, Usn (170); then the eighth period would close at 168. The tables therefore run to Z = 170 and give x by both formulas.

The names are the systematic IUPAC names, formed from the digits, until an element is discovered and named.

### 4.1 Shell and Position in the Valley

| Z | Name | Symbol | configuration by n + l | from two countings | ground term · moment [μ_B] | shell by n + l | Zα | A_valley | N_valley | N/Z | x drop | x with shell | neutron gap of fusion |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 119 | ununennium | Uue | [Og] 8s1 | as n + l | ²S₁/₂ · 1.00 | 8s | 0.868 | 321.1 | 202.1 | 1.698 | 0.914 | 0.729 | 20.4 |
| 120 | unbinilium | Ubn | [Og] 8s2 | as n + l | ¹S₀ · 0 | 8s | 0.876 | 324.4 | 204.4 | 1.703 | 0.920 | 0.733 | 20.6 |
| 121 | unbiunium | Ubu | [Og] 8s2 5g1 | [Og] 8s2 7d1 | ²D₃/₂ · 1.20 | 5g | 0.883 | 327.7 | 206.7 | 1.708 | 0.926 | 0.738 | 20.8 |
| 122 | unbibium | Ubb | [Og] 8s2 5g2 | [Og] 8s2 7d2 | ³F₂ · 1.33 | 5g | 0.890 | 331.0 | 209.0 | 1.713 | 0.932 | 0.743 | 21.0 |
| 123 | unbitrium | Ubt | [Og] 8s2 5g3 | [Og] 8s2 6f2 7d1 | ⁴K₁₁/₂ · 4.23 | 5g | 0.898 | 334.3 | 211.3 | 1.718 | 0.938 | 0.747 | 21.2 |
| 124 | unbiquadium | Ubq | [Og] 8s2 5g4 | – | ⁵N₈ · 6.22 | 5g | 0.905 | 337.6 | 213.6 | 1.723 | 0.944 | 0.752 | 21.3 |
| 125 | unbipentium | Ubp | [Og] 8s2 5g5 | – | ⁶N₁₅/₂ · 5.29 | 5g | 0.912 | 341.0 | 216.0 | 1.728 | 0.950 | 0.757 | 21.5 |
| 126 | unbihexium | Ubh | [Og] 8s2 5g6 | – | ⁷M₆ · 3.43 | 5g | 0.919 | 344.3 | 218.3 | 1.733 | 0.956 | 0.760 | 21.7 |
| 127 | unbiseptium | Ubs | [Og] 8s2 5g7 | – | ⁸K₇/₂ · 0.78 | 5g | 0.927 | 347.7 | 220.7 | 1.738 | 0.961 | 0.759 | 21.9 |
| 128 | unbioctium | Ubo | [Og] 8s2 5g8 | – | ⁹G₀ · 0.00 | 5g | 0.934 | 351.0 | 223.0 | 1.743 | 0.967 | 0.763 | 22.1 |
| 129 | unbiennium | Ube | [Og] 8s2 5g9 | – | ¹⁰S₉/₂ · 9.00 | 5g | 0.941 | 354.4 | 225.4 | 1.747 | 0.973 | 0.768 | 22.3 |
| 130 | untrinilium | Utn | [Og] 8s2 5g10 | – | ⁹G₈ · 12.00 | 5g | 0.949 | 357.8 | 227.8 | 1.752 | 0.979 | 0.772 | 22.5 |
| 131 | untriunium | Utu | [Og] 8s2 5g11 | – | ⁸K₂₁/₂ · 14.00 | 5g | 0.956 | 361.2 | 230.2 | 1.757 | 0.985 | 0.777 | 22.7 |
| 132 | untribium | Utb | [Og] 8s2 5g12 | – | ⁷M₁₂ · 15.00 | 5g | 0.963 | 364.6 | 232.6 | 1.762 | 0.990 | 0.781 | 22.9 |
| 133 | untritrium | Utt | [Og] 8s2 5g13 | – | ⁶N₂₅/₂ · 15.00 | 5g | 0.971 | 368.0 | 235.0 | 1.767 | 0.996 | 0.786 | 23.1 |
| 134 | untriquadium | Utq | [Og] 8s2 5g14 | – | ⁵N₁₂ · 14.00 | 5g | 0.978 | 371.4 | 237.4 | 1.772 | 1.002 | 0.790 | 23.2 |
| 135 | untripentium | Utp | [Og] 8s2 5g15 | – | ⁴M₂₁/₂ · 12.00 | 5g | 0.985 | 374.8 | 239.8 | 1.777 | 1.008 | 0.795 | 23.4 |
| 136 | untrihexium | Uth | [Og] 8s2 5g16 | – | ³K₈ · 9.00 | 5g | 0.992 | 378.3 | 242.3 | 1.782 | 1.013 | 0.799 | 23.6 |
| 137 | untriseptium | Uts | [Og] 8s2 5g17 | – | ²G₉/₂ · 5.00 | 5g | 1.000 | 381.7 | 244.7 | 1.786 | 1.019 | 0.804 | 23.8 |
| 138 | untrioctium | Uto | [Og] 8s2 5g18 | – | ¹S₀ · 0 | 5g | 1.007 | 385.2 | 247.2 | 1.791 | 1.025 | 0.808 | 24.0 |
| 139 | untriennium | Ute | [Og] 8s2 5g18 6f1 | – | ²F₅/₂ · 2.14 | 6f | 1.014 | 388.7 | 249.7 | 1.796 | 1.030 | 0.813 | 24.2 |
| 140 | unquadnilium | Uqn | [Og] 8s2 5g18 6f2 | – | ³H₄ · 3.20 | 6f | 1.022 | 392.1 | 252.1 | 1.801 | 1.036 | 0.817 | 24.4 |
| 141 | unquadunium | Uqu | [Og] 8s2 5g18 6f3 | – | ⁴I₉/₂ · 3.27 | 6f | 1.029 | 395.6 | 254.6 | 1.806 | 1.041 | 0.821 | 24.6 |
| 142 | unquadbium | Uqb | [Og] 8s2 5g18 6f4 | – | ⁵I₄ · 2.40 | 6f | 1.036 | 399.1 | 257.1 | 1.811 | 1.047 | 0.826 | 24.8 |
| 143 | unquadtrium | Uqt | [Og] 8s2 5g18 6f5 | – | ⁶H₅/₂ · 0.71 | 6f | 1.044 | 402.6 | 259.6 | 1.816 | 1.053 | 0.836 | 25.0 |
| 144 | unquadquadium | Uqq | [Og] 8s2 5g18 6f6 | – | ⁷F₀ · 0.00 | 6f | 1.051 | 406.1 | 262.1 | 1.820 | 1.058 | 0.843 | 25.2 |
| 145 | unquadpentium | Uqp | [Og] 8s2 5g18 6f7 | – | ⁸S₇/₂ · 7.00 | 6f | 1.058 | 409.7 | 264.7 | 1.825 | 1.064 | 0.848 | 25.3 |
| 146 | unquadhexium | Uqh | [Og] 8s2 5g18 6f8 | – | ⁷F₆ · 9.00 | 6f | 1.065 | 413.2 | 267.2 | 1.830 | 1.069 | 0.852 | 25.5 |
| 147 | unquadseptium | Uqs | [Og] 8s2 5g18 6f9 | – | ⁶H₁₅/₂ · 10.00 | 6f | 1.073 | 416.8 | 269.8 | 1.835 | 1.075 | 0.856 | 25.7 |
| 148 | unquadoctium | Uqo | [Og] 8s2 5g18 6f10 | – | ⁵I₈ · 10.00 | 6f | 1.080 | 420.3 | 272.3 | 1.840 | 1.080 | 0.860 | 25.9 |
| 149 | unquadennium | Uqe | [Og] 8s2 5g18 6f11 | – | ⁴I₁₅/₂ · 9.00 | 6f | 1.087 | 423.9 | 274.9 | 1.845 | 1.085 | 0.865 | 26.1 |
| 150 | unpentnilium | Upn | [Og] 8s2 5g18 6f12 | – | ³H₆ · 7.00 | 6f | 1.095 | 427.4 | 277.4 | 1.850 | 1.091 | 0.869 | 26.3 |
| 151 | unpentunium | Upu | [Og] 8s2 5g18 6f13 | – | ²F₇/₂ · 4.00 | 6f | 1.102 | 431.0 | 280.0 | 1.854 | 1.096 | 0.873 | 26.5 |
| 152 | unpentbium | Upb | [Og] 8s2 5g18 6f14 | – | ¹S₀ · 0 | 6f | 1.109 | 434.6 | 282.6 | 1.859 | 1.102 | 0.877 | 26.7 |
| 153 | unpenttrium | Upt | [Og] 8s2 5g18 6f14 7d1 | – | ²D₃/₂ · 1.20 | 7d | 1.116 | 438.2 | 285.2 | 1.864 | 1.107 | 0.882 | 26.9 |
| 154 | unpentquadium | Upq | [Og] 8s2 5g18 6f14 7d2 | – | ³F₂ · 1.33 | 7d | 1.124 | 441.8 | 287.8 | 1.869 | 1.112 | 0.886 | 27.1 |
| 155 | unpentpentium | Upp | [Og] 8s2 5g18 6f14 7d3 | – | ⁴F₃/₂ · 0.60 | 7d | 1.131 | 445.5 | 290.5 | 1.874 | 1.118 | 0.890 | 27.3 |
| 156 | unpenthexium | Uph | [Og] 8s2 5g18 6f14 7d4 | – | ⁵D₀ · 0.00 | 7d | 1.138 | 449.1 | 293.1 | 1.879 | 1.123 | 0.894 | 27.5 |
| 157 | unpentseptium | Ups | [Og] 8s2 5g18 6f14 7d5 | – | ⁶S₅/₂ · 5.00 | 7d | 1.146 | 452.7 | 295.7 | 1.884 | 1.128 | 0.898 | 27.6 |
| 158 | unpentoctium | Upo | [Og] 8s2 5g18 6f14 7d6 | – | ⁵D₄ · 6.00 | 7d | 1.153 | 456.4 | 298.4 | 1.888 | 1.134 | 0.902 | 27.8 |
| 159 | unpentennium | Upe | [Og] 8s2 5g18 6f14 7d7 | – | ⁴F₉/₂ · 6.00 | 7d | 1.160 | 460.0 | 301.0 | 1.893 | 1.139 | 0.906 | 28.0 |
| 160 | unhexnilium | Uhn | [Og] 8s2 5g18 6f14 7d8 | – | ³F₄ · 5.00 | 7d | 1.168 | 463.7 | 303.7 | 1.898 | 1.144 | 0.911 | 28.2 |
| 161 | unhexunium | Uhu | [Og] 8s2 5g18 6f14 7d9 | – | ²D₅/₂ · 3.00 | 7d | 1.175 | 467.4 | 306.4 | 1.903 | 1.149 | 0.915 | 28.4 |
| 162 | unhexbium | Uhb | [Og] 8s2 5g18 6f14 7d10 | – | ¹S₀ · 0 | 7d | 1.182 | 471.1 | 309.1 | 1.908 | 1.155 | 0.919 | 28.6 |
| 163 | unhextrium | Uht | [Og] 8s2 5g18 6f14 7d10 8p1 | – | ²P₁/₂ · 0.33 | 8p | 1.189 | 474.8 | 311.8 | 1.913 | 1.160 | 0.923 | 28.8 |
| 164 | unhexquadium | Uhq | [Og] 8s2 5g18 6f14 7d10 8p2 | – | ³P₀ · 0.00 | 8p | 1.197 | 478.5 | 314.5 | 1.918 | 1.165 | 0.927 | 29.0 |
| 165 | unhexpentium | Uhp | [Og] 8s2 5g18 6f14 7d10 8p3 | – | ⁴S₃/₂ · 3.00 | 8p | 1.204 | 482.2 | 317.2 | 1.922 | 1.170 | 0.931 | 29.2 |
| 166 | unhexhexium | Uhh | [Og] 8s2 5g18 6f14 7d10 8p4 | – | ³P₂ · 3.00 | 8p | 1.211 | 485.9 | 319.9 | 1.927 | 1.175 | 0.935 | 29.4 |
| 167 | unhexseptium | Uhs | [Og] 8s2 5g18 6f14 7d10 8p5 | – | ²P₃/₂ · 2.00 | 8p | 1.219 | 489.7 | 322.7 | 1.932 | 1.180 | 0.939 | 29.6 |
| 168 | unhexoctium | Uho | [Og] 8s2 5g18 6f14 7d10 8p6 | – | ¹S₀ · 0 | 8p | 1.226 | 493.4 | 325.4 | 1.937 | 1.185 | 0.943 | 29.8 |
| 169 | unhexennium | Uhe | [Og] 8s2 5g18 6f14 7d10 8p6 6g1 | – | ²S₁/₂ · 1.00 | 6g | 1.233 | 497.2 | 328.2 | 1.942 | 1.191 | 0.947 | 30.0 |
| 170 | unseptnilium | Usn | [Og] 8s2 5g18 6f14 7d10 8p6 6g2 | – | ¹S₀ · 0 | 6g | 1.241 | 500.9 | 330.9 | 1.947 | 1.196 | 0.951 | 30.2 |

Configuration by n + l and from the two countings (search without measurement as in 3.3; – not yet computed); ground term from the computed configuration, otherwise by n + l. A_valley and N_valley by the liquid-drop formula; the valley with the shell term lies the same to within a few nucleons. The neutron gap is N_valley(Z) minus the largest neutron number that two nuclei from the valley (each Z ≥ 20) bring along when they fuse.

### 4.2 Properties of the Nucleus in the Valley

For the nucleus in the valley (A rounded to an integer) the formula gives the binding per nucleon B/A, the alpha-decay energy Q_α and the neutron separation energy S_n. Tested on measured nuclei (AME2020): Q_α on 45 decays from Z = 100 off by 0.12 MeV on average, scatter 0.27 MeV; S_n on 187 nuclei from Z = 90 off by 0.31 MeV on average, scatter 0.43 MeV. The lifetime lg T_α (4.3) holds for even-even nuclei; beyond the tipping point of the respective formula the nucleus is not held.

| Z | Symbol | A | B/A [keV] | Q_α [MeV] | Q_α with shell [MeV] | S_n [MeV] | lg T_α [s] |
|---|---|---|---|---|---|---|---|
| 119 | Uue | 321 | 7008 | 9.05 | 9.43 | 5.83 | 3.5 to 4.7 |
| 120 | Ubn | 324 | 6990 | 9.23 | 9.55 | 5.81 | 3.4 to 4.4 |
| 121 | Ubu | 328 | 6960 | 9.26 | 9.49 | 4.43 | 3.9 to 4.6 |
| 122 | Ubb | 331 | 6943 | 9.44 | 9.61 | 4.42 | 3.8 to 4.4 |
| 123 | Ubt | 334 | 6921 | 9.62 | 9.73 | 4.42 | 3.8 to 4.1 |
| 124 | Ubq | 338 | 6899 | 9.67 | 9.70 | 5.60 | 4.2 to 4.3 |
| 125 | Ubp | 341 | 6878 | 9.85 | 9.81 | 5.59 | 4.0 to 4.1 |
| 126 | Ubh | 344 | 6860 | 10.03 | 9.94 | 5.58 | 3.8 to 4.0 |
| 127 | Ubs | 348 | 6832 | 10.07 | 10.55 | 4.25 | 2.5 to 3.9 |
| 128 | Ubo | 351 | 6814 | 10.25 | 11.32 | 4.24 | 0.7 to 3.7 |
| 129 | Ube | 354 | 6793 | 10.43 | 11.43 | 4.24 | 0.7 to 3.4 |
| 130 | Utn | 358 | 6772 | 10.49 | 11.41 | 5.40 | 1.0 to 3.6 |
| 131 | Utu | 361 | 6751 | 10.67 | 11.52 | 5.39 | 1.0 to 3.3 |
| 132 | Utb | 365 | 6727 | 10.72 | 11.50 | 4.09 | 1.3 to 3.5 |
| 133 | Utt | 368 | 6706 | 10.90 | 11.62 | 4.09 | 1.3 to 3.2 |
| 134 | Utq | 371 | 6689 | 11.08 | 11.75 | 4.09 | 1.2 to 3.0 |
| 135 | Utp | 375 | 6665 | 11.14 | 11.73 | 5.23 | 1.5 to 3.1 |
| 136 | Uth | 378 | 6648 | 11.32 | 11.86 | 5.22 | 1.5 to 2.9 |
| 137 | Uts | 382 | 6621 | 11.38 | 11.84 | 3.95 | 1.8 to 3.0 |
| 138 | Uto | 385 | 6604 | 11.56 | 11.98 | 3.95 | 1.7 to 2.8 |
| 139 | Ute | 389 | 6580 | 11.62 | 11.97 | 5.07 | 2.0 to 2.9 |
| 140 | Uqn | 392 | 6563 | 11.81 | 12.11 | 5.07 | 1.9 to 2.6 |
| 141 | Uqu | 396 | 6537 | 11.87 | 12.11 | 3.83 | 2.1 to 2.8 |
| 142 | Uqb | 399 | 6520 | 12.06 | 12.25 | 3.83 | 2.0 to 2.5 |
| 143 | Uqt | 403 | 6497 | 12.12 | 13.75 | 4.94 | −1.2 to 2.6 |
| 144 | Uqq | 406 | 6480 | 12.31 | 13.88 | 4.93 | −1.3 to 2.4 |
| 145 | Uqp | 410 | 6454 | 12.37 | 13.85 | 3.72 | −1.0 to 2.5 |
| 146 | Uqh | 413 | 6438 | 12.56 | 13.98 | 3.72 | −1.0 to 2.3 |
| 147 | Uqs | 417 | 6415 | 12.63 | 13.97 | 4.81 | −0.8 to 2.3 |
| 148 | Uqo | 420 | 6399 | 12.81 | 14.11 | 4.81 | −0.9 to 2.1 |
| 149 | Uqe | 424 | 6373 | 12.88 | 14.10 | 3.62 | −0.6 to 2.2 |
| 150 | Upn | 427 | 6357 | 13.07 | 14.23 | 3.62 | −0.7 to 2.0 |
| 151 | Upu | 431 | 6335 | 13.14 | 14.23 | 4.70 | −0.5 to 2.0 |
| 152 | Upb | 435 | 6312 | 13.21 | 14.23 | 3.52 | −0.2 to 2.1 |
| 153 | Upt | 438 | 6294 | 13.40 | 14.37 | 3.53 | −0.3 to 1.9 |
| 154 | Upq | 442 | 6274 | 13.48 | 14.38 | 4.59 | −0.1 to 1.9 |
| 155 | Upp | 445 | 6256 | 13.66 | 14.52 | 4.59 | −0.2 to 1.7 |
| 156 | Uph | 449 | 6234 | 13.74 | 14.53 | 3.44 | 0.0 to 1.8 |
| 157 | Ups | 453 | 6212 | 13.81 | 14.54 | 4.49 | 0.2 to 1.8 |
| 158 | Upo | 456 | 6196 | 14.00 | 14.69 | 4.49 | 0.1 to 1.6 |
| 159 | Upe | 460 | 6172 | 14.08 | 14.70 | 3.35 | 0.3 to 1.7 |
| 160 | Uhn | 464 | 6153 | 14.16 | 14.73 | 4.39 | 0.4 to 1.7 |
| 161 | Uhu | 467 | 6135 | 14.35 | 14.87 | 4.40 | 0.3 to 1.5 |
| 162 | Uhb | 471 | 6113 | 14.43 | 14.89 | 3.28 | 0.5 to 1.5 |
| 163 | Uht | 475 | 6092 | 14.51 | 14.92 | 4.31 | 0.7 to 1.6 |
| 164 | Uhq | 478 | 6077 | 14.70 | 15.07 | 4.31 | 0.5 to 1.4 |
| 165 | Uhp | 482 | 6053 | 14.78 | 15.10 | 3.20 | 0.7 to 1.4 |
| 166 | Uhh | 486 | 6035 | 14.87 | 15.13 | 4.22 | 0.8 to 1.4 |
| 167 | Uhs | 490 | 6011 | 14.95 | 15.16 | 3.13 | 1.0 to 1.4 |
| 168 | Uho | 493 | 5996 | 15.14 | 15.32 | 3.13 | 0.9 to 1.2 |
| 169 | Uhe | 497 | 5976 | 15.23 | 15.35 | 4.14 | 1.0 to 1.3 |
| 170 | Usn | 501 | 5955 | 15.31 | 15.39 | 3.06 | 1.1 to 1.3 |

In the secured range up to Z = 133, Q_α lies at 9.0 to 10.9 MeV, with the shell term at 9.4 to 11.6 MeV; the heaviest measured nucleus Ds-270 has 11.12 MeV. S_n stays positive up to Z = 170, at 3.1 to 5.8 MeV: the nuclei in the valley hold their neutrons. The alternation of S_n follows the pairing.

### 4.3 Shells, Relativity, Lifetime

**Shell correction in the nucleus.** The liquid-drop formula knows no closures. Inserting the closure numbers of the model (2, 8, 20, 28, 50, 82, 126, 184, 258) as a shell term, in the form of Myers and Swiatecki, with Wigner and surface-asymmetry terms, lowers the scatter of the binding on 2509 measured nuclei from 3.39 to 1.72 MeV. S_n improves (0.43 → 0.31 MeV from Z = 90), Q_α from Z = 84 too (1.43 → 1.07 MeV), but not from Z = 100 (0.27 → 0.47 MeV): the form lacks the deformation of the nuclei between the shells. This is why 4.2 gives both values; their difference is the uncertainty. Whether the next proton closure lies at 114, 120 or 126 is not decided by the measured nuclei: scatter from Z = 90 2.15, 2.22 and 2.28 MeV. The valley nuclei from Uue to Utt, with N ≈ 202–235, lie between the neutron closures 184 and 258.

**Relativity of the shell.** By Dirac, the energy of the innermost electron around a point nucleus is (3c⁴·Δs_e/G)·√(1 − (Zα)²), with the fold Δs_e of the electron; at Zα = 1, Z = 137, it ends. Around the extended nucleus with R = λ_H·A^{1/3} the state stays bound and reaches the binding 6c⁴·Δs_e/G, twice the fold energy of the electron, only at Z = 170.1 to 170.9 (Dirac equation for a uniformly charged sphere, computed by two routes; the range comes from the holder length and the valley without and with the shell term).

**Lifetime.** Alpha decay is a passage through the Coulomb barrier, the tunnelling factor e^{−√(E_G/E)} of the model. From it follows lg T_α = a·Z/√Q_α + b·Z + c. On 82 even-even alpha emitters from Z = 84 (NUBASE2020, AME2020), a = 1.526, b = −0.1831, c = −33.89 (T in s), scatter 0.51 decades; it meets U-238 to 0.06 decades. With the two Q_α of 4.2 it gives the range lg T_α of the table, for even-even nuclei; odd ones live longer. In the secured range, Uue to Utt, the valley nuclei would hold against alpha decay for 10^{0.7} to 10^{4.7} s, 5.2 seconds to 14 hours; the range per element is given in 4.2. Against fission there is no such form: on 43 measured nuclei the fission lifetime does fall with x (correlation −0.89), but scatters by 4.1 decades, and measurements reach only x = 0.942 (liquid drop). How long the new elements hold against fission is predicted neither by the model nor by the data; near x = 1 the restoring vanishes.

### 4.4 How, Where, How Long

**How they can form.** The valley curves: N_valley grows faster than Z. Two nuclei from the valley therefore always bring too few neutrons when they fuse; for Z = 119 to 170 the table shows 20 to 30 missing. Nuclei made by fusion lie neutron-poor and tip earlier, like Og-294. Only the stepwise addition of neutrons to heavy nuclei leads into the valley, until N ≈ 202 is reached for Uue (Z = 119), N ≈ 235 for Utt (Z = 133) and N ≈ 331 for Usn (Z = 170), with conversions n → p that raise Z. Neutron-rich nuclei hold better, x = a_C·Z²/(2a_S·A) falls with A: between the valley and the drip edge (the next neutron pair no longer bound) the nucleus stays held by the liquid-drop formula up to Z = 172 (edge at N = 442, x = 0.999), with the shell term beyond Z = 199; the valley itself only up to 133 and 181. On a neutron-rich path the shell (Z = 170.9) thus sets the end, not the nucleus. The edges lie far beyond the measured masses; they are a continuation of the fitted formulas.

**Where the elements form.** Each holder of the funnel does its part. *Heat:* in the star the nuclei burn as a ladder in the depth as long as the binding per nucleon grows, up to x = ¼; in the valley this is Z = 26.7, A = 58.9 by the liquid-drop formula, between Fe-56 (smallest mass per nucleon) and Ni-62 (largest binding per nucleon, AME2020); with the shell term Z = 34.7. Up to the iron group the elements form by burning. *Counting rule of the electrons:* a funnel whose depth it holds reaches at most kT = 0.0964·(μ/μ_e)·m_ec², for carbon and oxygen 43 keV (5.0·10⁸ K), at the limiting mass 1.456 M☉. *Counting rule of the nucleons:* beyond the iron group only neutron capture leads on (How, Where).

**Macro state: where the shell holds.** The atoms as mass fold the body, δ = GM/(Rc²) = 3Δs/R; the body gives every atom pressure, P·V_A = ½·m_A·c²·δ (centre of the homogeneous sphere), and as a heat-held funnel heat, kT = f·δ·m̄c² (f = 0.8543 from the polytrope n = 3). The shell holds as long as both stay below the binding of the outer electron IE (PubChem): up to δ_P = 2·IE/(m_A·c²) and δ_T = IE/(f·m_A·c²). Against them stand the bodies: Earth δ = 6.96·10^-10, Sun δ = 2.12·10^-6, neutron star 1.4 M☉, 12 km δ = 1.72·10^-1. The thresholds of all elements lie between 3.8·10⁻¹¹ and 2.9·10⁻⁸: where the elements form, in the star and the neutron star, every shell is open; they form as nuclei, and their shell, configuration and chemistry form only in flat regions, in the limit δ → 0 of the gas, where the configurations are measured. Conversely a body of the element with its measured density ρ reaches δ_P in its centre from the mass M_krit = (δ_P·c²/G)^{3/2}·(3/(4πρ))^{1/2} (elements that are solid or liquid). The first column names the range of the holder (burning ladder up to x = ¼, beyond it neutron capture), not the path of each nucleus. For all elements with a measured ionisation energy.

| Z | | range of the holder | IE [eV] | δ_P | δ_T | M_krit [Earth masses] | shell in the Earth (centre) |
|---|---|---|---|---|---|---|---|
| 1 | H | start of the ladder | 13.60 | 2.90·10^-8 | 1.70·10^-8 | – | holds |
| 2 | He | burning | 24.59 | 1.32·10^-8 | 7.72·10^-9 | – | holds |
| 3 | Li | burning | 5.39 | 1.65·10^-9 | 9.68·10^-10 | 11.767 | holds |
| 4 | Be | burning | 9.32 | 2.22·10^-9 | 1.30·10^-9 | 9.839 | holds |
| 5 | B | burning | 8.30 | 1.65·10^-9 | 9.65·10^-10 | 5.556 | holds |
| 6 | C | burning | 11.26 | 2.01·10^-9 | 1.18·10^-9 | 7.668 | holds |
| 7 | N | burning | 14.53 | 2.23·10^-9 | 1.30·10^-9 | – | holds |
| 8 | O | burning | 13.62 | 1.83·10^-9 | 1.07·10^-9 | – | holds |
| 9 | F | burning | 17.42 | 1.97·10^-9 | 1.15·10^-9 | – | holds |
| 10 | Ne | burning | 21.57 | 2.29·10^-9 | 1.34·10^-9 | – | holds |
| 11 | Na | burning | 5.14 | 4.80·10^-10 | 2.81·10^-10 | 1.365 | open |
| 12 | Mg | burning | 7.65 | 6.75·10^-10 | 3.95·10^-10 | 1.701 | open |
| 13 | Al | burning | 5.99 | 4.76·10^-10 | 2.79·10^-10 | 0.809 | open |
| 14 | Si | burning | 8.15 | 6.23·10^-10 | 3.65·10^-10 | 1.303 | open |
| 15 | P | burning | 10.49 | 7.27·10^-10 | 4.25·10^-10 | 1.857 | holds |
| 16 | S | burning | 10.36 | 6.94·10^-10 | 4.06·10^-10 | 1.624 | open |
| 17 | Cl | burning | 12.97 | 7.85·10^-10 | 4.60·10^-10 | – | holds |
| 18 | Ar | burning | 15.76 | 8.48·10^-10 | 4.96·10^-10 | – | holds |
| 19 | K | burning | 4.34 | 2.38·10^-10 | 1.40·10^-10 | 0.499 | open |
| 20 | Ca | burning | 6.11 | 3.27·10^-10 | 1.92·10^-10 | 0.610 | open |
| 21 | Sc | burning | 6.56 | 3.13·10^-10 | 1.83·10^-10 | 0.410 | open |
| 22 | Ti | burning | 6.83 | 3.06·10^-10 | 1.79·10^-10 | 0.323 | open |
| 23 | V | burning | 6.75 | 2.84·10^-10 | 1.66·10^-10 | 0.250 | open |
| 24 | Cr | burning | 6.77 | 2.79·10^-10 | 1.64·10^-10 | 0.223 | open |
| 25 | Mn | burning | 7.43 | 2.91·10^-10 | 1.70·10^-10 | 0.234 | open |
| 26 | Fe | burning | 7.90 | 3.04·10^-10 | 1.78·10^-10 | 0.241 | open |
| 27 | Co | neutron capture | 7.88 | 2.87·10^-10 | 1.68·10^-10 | 0.209 | open |
| 28 | Ni | neutron capture | 7.64 | 2.79·10^-10 | 1.64·10^-10 | 0.200 | open |
| 29 | Cu | neutron capture | 7.73 | 2.61·10^-10 | 1.53·10^-10 | 0.180 | open |
| 30 | Zn | neutron capture | 9.39 | 3.08·10^-10 | 1.80·10^-10 | 0.259 | open |
| 31 | Ga | neutron capture | 6.00 | 1.85·10^-10 | 1.08·10^-10 | 0.132 | open |
| 32 | Ge | neutron capture | 7.90 | 2.34·10^-10 | 1.37·10^-10 | 0.198 | open |
| 33 | As | neutron capture | 9.81 | 2.81·10^-10 | 1.65·10^-10 | 0.251 | open |
| 34 | Se | neutron capture | 9.75 | 2.65·10^-10 | 1.55·10^-10 | 0.252 | open |
| 35 | Br | neutron capture | 11.81 | 3.17·10^-10 | 1.86·10^-10 | 0.410 | open |
| 36 | Kr | neutron capture | 14.00 | 3.59·10^-10 | 2.10·10^-10 | – | open |
| 37 | Rb | neutron capture | 4.18 | 1.05·10^-10 | 6.14·10^-11 | 0.111 | open |
| 38 | Sr | neutron capture | 5.70 | 1.40·10^-10 | 8.17·10^-11 | 0.130 | open |
| 39 | Y | neutron capture | 6.22 | 1.50·10^-10 | 8.79·10^-11 | 0.111 | open |
| 40 | Zr | neutron capture | 6.63 | 1.56·10^-10 | 9.14·10^-11 | 0.098 | open |
| 41 | Nb | neutron capture | 6.76 | 1.56·10^-10 | 9.14·10^-11 | 0.085 | open |
| 42 | Mo | neutron capture | 7.09 | 1.59·10^-10 | 9.29·10^-11 | 0.080 | open |
| 43 | Tc | neutron capture | 7.28 | 1.61·10^-10 | 9.44·10^-11 | 0.079 | open |
| 44 | Ru | neutron capture | 7.36 | 1.56·10^-10 | 9.15·10^-11 | 0.072 | open |
| 45 | Rh | neutron capture | 7.46 | 1.56·10^-10 | 9.11·10^-11 | 0.070 | open |
| 46 | Pd | neutron capture | 8.34 | 1.68·10^-10 | 9.84·10^-11 | 0.081 | open |
| 47 | Ag | neutron capture | 7.58 | 1.51·10^-10 | 8.83·10^-11 | 0.073 | open |
| 48 | Cd | neutron capture | 8.99 | 1.72·10^-10 | 1.01·10^-10 | 0.098 | open |
| 49 | In | neutron capture | 5.79 | 1.08·10^-10 | 6.33·10^-11 | 0.053 | open |
| 50 | Sn | neutron capture | 7.34 | 1.33·10^-10 | 7.77·10^-11 | 0.073 | open |
| 51 | Sb | neutron capture | 8.64 | 1.52·10^-10 | 8.92·10^-11 | 0.093 | open |
| 52 | Te | neutron capture | 9.01 | 1.52·10^-10 | 8.87·10^-11 | 0.096 | open |
| 53 | I | neutron capture | 10.45 | 1.77·10^-10 | 1.03·10^-10 | 0.135 | open |
| 54 | Xe | neutron capture | 12.13 | 1.98·10^-10 | 1.16·10^-10 | – | open |
| 55 | Cs | neutron capture | 3.89 | 6.29·10^-11 | 3.68·10^-11 | 0.046 | open |
| 56 | Ba | neutron capture | 5.21 | 8.15·10^-11 | 4.77·10^-11 | 0.049 | open |
| 57 | La | neutron capture | 5.58 | 8.62·10^-11 | 5.05·10^-11 | 0.041 | open |
| 58 | Ce | neutron capture | 5.54 | 8.49·10^-11 | 4.97·10^-11 | 0.038 | open |
| 59 | Pr | neutron capture | 5.46 | 8.33·10^-11 | 4.87·10^-11 | 0.037 | open |
| 60 | Nd | neutron capture | 5.53 | 8.22·10^-11 | 4.81·10^-11 | 0.036 | open |
| 61 | Pm | neutron capture | 5.55 | 8.22·10^-11 | 4.81·10^-11 | 0.035 | open |
| 62 | Sm | neutron capture | 5.64 | 8.06·10^-11 | 4.72·10^-11 | 0.034 | open |
| 63 | Eu | neutron capture | 5.67 | 8.01·10^-11 | 4.69·10^-11 | 0.040 | open |
| 64 | Gd | neutron capture | 6.15 | 8.40·10^-11 | 4.91·10^-11 | 0.035 | open |
| 65 | Tb | neutron capture | 5.86 | 7.92·10^-11 | 4.64·10^-11 | 0.031 | open |
| 66 | Dy | neutron capture | 5.94 | 7.85·10^-11 | 4.59·10^-11 | 0.030 | open |
| 67 | Ho | neutron capture | 6.02 | 7.84·10^-11 | 4.59·10^-11 | 0.030 | open |
| 68 | Er | neutron capture | 6.11 | 7.84·10^-11 | 4.59·10^-11 | 0.029 | open |
| 69 | Tm | neutron capture | 6.18 | 7.86·10^-11 | 4.60·10^-11 | 0.029 | open |
| 70 | Yb | neutron capture | 6.25 | 7.76·10^-11 | 4.54·10^-11 | 0.033 | open |
| 71 | Lu | neutron capture | 5.43 | 6.66·10^-11 | 3.90·10^-11 | 0.022 | open |
| 72 | Hf | neutron capture | 6.83 | 8.21·10^-11 | 4.80·10^-11 | 0.026 | open |
| 73 | Ta | neutron capture | 7.89 | 9.36·10^-11 | 5.48·10^-11 | 0.029 | open |
| 74 | W | neutron capture | 7.98 | 9.32·10^-11 | 5.45·10^-11 | 0.026 | open |
| 75 | Re | neutron capture | 7.88 | 9.09·10^-11 | 5.32·10^-11 | 0.024 | open |
| 76 | Os | neutron capture | 8.70 | 9.82·10^-11 | 5.75·10^-11 | 0.026 | open |
| 77 | Ir | neutron capture | 9.10 | 1.02·10^-10 | 5.95·10^-11 | 0.028 | open |
| 78 | Pt | neutron capture | 9.00 | 9.91·10^-11 | 5.80·10^-11 | 0.027 | open |
| 79 | Au | neutron capture | 9.23 | 1.01·10^-10 | 5.89·10^-11 | 0.029 | open |
| 80 | Hg | neutron capture | 10.44 | 1.12·10^-10 | 6.54·10^-11 | 0.041 | open |
| 81 | Tl | neutron capture | 6.11 | 6.42·10^-11 | 3.76·10^-11 | 0.019 | open |
| 82 | Pb | neutron capture | 7.42 | 7.69·10^-11 | 4.50·10^-11 | 0.026 | open |
| 83 | Bi | neutron capture | 7.29 | 7.49·10^-11 | 4.38·10^-11 | 0.026 | open |
| 84 | Po | neutron capture | 8.42 | 8.65·10^-11 | 5.06·10^-11 | 0.034 | open |
| 85 | At | neutron capture | 9.50 | 9.71·10^-11 | 5.69·10^-11 | 0.046 | open |
| 86 | Rn | neutron capture | 10.74 | 1.04·10^-10 | 6.08·10^-11 | – | open |
| 87 | Fr | neutron capture | 3.90 | 3.75·10^-11 | 2.20·10^-11 | – | open |
| 88 | Ra | neutron capture | 5.28 | 5.01·10^-11 | 2.93·10^-11 | 0.020 | open |
| 89 | Ac | neutron capture | 5.17 | 4.89·10^-11 | 2.86·10^-11 | 0.014 | open |
| 90 | Th | neutron capture | 6.08 | 5.63·10^-11 | 3.29·10^-11 | 0.016 | open |
| 91 | Pa | neutron capture | 5.89 | 5.47·10^-11 | 3.20·10^-11 | 0.013 | open |
| 92 | U | neutron capture | 6.19 | 5.59·10^-11 | 3.27·10^-11 | 0.012 | open |
| 93 | Np | neutron capture | 6.27 | 5.68·10^-11 | 3.32·10^-11 | 0.012 | open |
| 94 | Pu | neutron capture | 6.06 | 5.33·10^-11 | 3.12·10^-11 | 0.011 | open |
| 95 | Am | neutron capture | 5.99 | 5.29·10^-11 | 3.10·10^-11 | 0.013 | open |
| 96 | Cm | neutron capture | 6.02 | 5.23·10^-11 | 3.06·10^-11 | 0.013 | open |
| 97 | Bk | neutron capture | 6.23 | 5.41·10^-11 | 3.17·10^-11 | 0.014 | open |
| 98 | Cf | neutron capture | 6.30 | 5.39·10^-11 | 3.15·10^-11 | – | open |
| 99 | Es | neutron capture | 6.42 | 5.47·10^-11 | 3.20·10^-11 | – | open |
| 100 | Fm | neutron capture | 6.50 | 5.43·10^-11 | 3.18·10^-11 | – | open |
| 101 | Md | neutron capture | 6.58 | 5.47·10^-11 | 3.20·10^-11 | – | open |
| 102 | No | neutron capture | 6.65 | 5.51·10^-11 | 3.23·10^-11 | – | open |
| 103 | Unt | neutron capture | – | – | – | – | – |
| 104 | Unq | neutron capture | – | – | – | – | – |
| 105 | Unp | neutron capture | – | – | – | – | – |
| 106 | Unh | neutron capture | – | – | – | – | – |
| 107 | Uns | neutron capture | – | – | – | – | – |
| 108 | Uno | neutron capture | – | – | – | – | – |
| 109 | Une | neutron capture | – | – | – | – | – |
| 110 | Uun | neutron capture | – | – | – | – | – |
| 111 | Uuu | neutron capture | – | – | – | – | – |
| 112 | Uub | neutron capture | – | – | – | – | – |
| 113 | Uut | neutron capture | – | – | – | – | – |
| 114 | Uuq | neutron capture | – | – | – | – | – |
| 115 | Uup | neutron capture | – | – | – | – | – |
| 116 | Uuh | neutron capture | – | – | – | – | – |
| 117 | Uus | neutron capture | – | – | – | – | – |
| 118 | Uuo | neutron capture | – | – | – | – | – |
| 119 | Uue | neutron capture | – | – | – | – | – |
| 120 | Ubn | neutron capture | – | – | – | – | – |
| 121 | Ubu | neutron capture | – | – | – | – | – |
| 122 | Ubb | neutron capture | – | – | – | – | – |
| 123 | Ubt | neutron capture | – | – | – | – | – |
| 124 | Ubq | neutron capture | – | – | – | – | – |
| 125 | Ubp | neutron capture | – | – | – | – | – |
| 126 | Ubh | neutron capture | – | – | – | – | – |
| 127 | Ubs | neutron capture | – | – | – | – | – |
| 128 | Ubo | neutron capture | – | – | – | – | – |
| 129 | Ube | neutron capture | – | – | – | – | – |
| 130 | Utn | neutron capture | – | – | – | – | – |
| 131 | Utu | neutron capture | – | – | – | – | – |
| 132 | Utb | neutron capture | – | – | – | – | – |
| 133 | Utt | neutron capture | – | – | – | – | – |
| 134 | Utq | neutron capture | – | – | – | – | – |
| 135 | Utp | neutron capture | – | – | – | – | – |
| 136 | Uth | neutron capture | – | – | – | – | – |
| 137 | Uts | neutron capture | – | – | – | – | – |
| 138 | Uto | neutron capture | – | – | – | – | – |
| 139 | Ute | neutron capture | – | – | – | – | – |
| 140 | Uqn | neutron capture | – | – | – | – | – |
| 141 | Uqu | neutron capture | – | – | – | – | – |
| 142 | Uqb | neutron capture | – | – | – | – | – |
| 143 | Uqt | neutron capture | – | – | – | – | – |
| 144 | Uqq | neutron capture | – | – | – | – | – |
| 145 | Uqp | neutron capture | – | – | – | – | – |
| 146 | Uqh | neutron capture | – | – | – | – | – |
| 147 | Uqs | neutron capture | – | – | – | – | – |
| 148 | Uqo | neutron capture | – | – | – | – | – |
| 149 | Uqe | neutron capture | – | – | – | – | – |
| 150 | Upn | neutron capture | – | – | – | – | – |
| 151 | Upu | neutron capture | – | – | – | – | – |
| 152 | Upb | neutron capture | – | – | – | – | – |
| 153 | Upt | neutron capture | – | – | – | – | – |
| 154 | Upq | neutron capture | – | – | – | – | – |
| 155 | Upp | neutron capture | – | – | – | – | – |
| 156 | Uph | neutron capture | – | – | – | – | – |
| 157 | Ups | neutron capture | – | – | – | – | – |
| 158 | Upo | neutron capture | – | – | – | – | – |
| 159 | Upe | neutron capture | – | – | – | – | – |
| 160 | Uhn | neutron capture | – | – | – | – | – |
| 161 | Uhu | neutron capture | – | – | – | – | – |
| 162 | Uhb | neutron capture | – | – | – | – | – |
| 163 | Uht | neutron capture | – | – | – | – | – |
| 164 | Uhq | neutron capture | – | – | – | – | – |
| 165 | Uhp | neutron capture | – | – | – | – | – |
| 166 | Uhh | neutron capture | – | – | – | – | – |
| 167 | Uhs | neutron capture | – | – | – | – | – |
| 168 | Uho | neutron capture | – | – | – | – | – |
| 169 | Uhe | neutron capture | – | – | – | – | – |
| 170 | Usn | neutron capture | – | – | – | – | – |

**Where.** This requires free neutrons at high density. In the model they exist where the counting rule of the nucleons itself holds the matter, the held funnel with a spread of order one. Whether and how fast nuclei grow there up to A ≈ 370 and beyond the model does not compute.

**How long.** The restoring against deformation is ∝ 1 − x. In the valley, by the liquid-drop formula it is 0.086 for Uue and 0.004 for Utt, zero at the tipping point; with the shell term 0.271 for Uue. For comparison, U-238: 0.263 by the liquid-drop formula, 0.417 with the shell term.

**What remains open and what is not possible.** Above Z = 170 the shell ends. Between 133 and 170 everything depends on the shell term. The double closure unbihexium (Ubh, Z = 126) with N = 184 is neutron-poor: x = 1.06 by the liquid-drop formula, 0.84 with the shell term.

**Limit of the calculation.** The nuclei above A = 295 are a continuation of the fitted formulas, not a measurement.
