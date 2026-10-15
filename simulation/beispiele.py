"""Drei gerechnete Beispiele an unseren Daten: Atom und Kern, eine Galaxie, der Kosmos."""
import math
import numpy as np

import kern as K
import daten as D


def atom_und_kern():
    """Wasserstoff: Umlauf mit v = αc auf a_B, Bindung ½μv² gegen PubChem. Kern: Halter λ_H aus AME2020 (Simulation der Kerne)."""
    import kerne as KE
    mu = K.m_e * K.m_p / (K.m_e + K.m_p)
    E_H = 0.5 * mu * (K.alpha * K.c)**2 / K.e
    aC = KE.a_C(KE.TROPFEN)
    lam_H = KE.halterlaenge(aC)
    lam_p = K.hbar / (K.m_p * K.c) * 1e15
    return dict(E_H=E_H, IE_H=D.pubchem()["H"]["IE"], aC=aC, lam_H=lam_H,
                R_Fe=lam_H * 56**(1 / 3), v_groesse=lam_p / lam_H,
                v_bindung=math.sqrt(2 * KE.bekannte()["b_je_a"] / (K.energie(K.m_u) / K.e / 1e6)),
                lam_spiegel=KE.halterlaenge(D.spiegelkerne()[1]))


def galaxie(name="NGC2403", radien=(1.0, 3.0, 6.0, 10.0, 15.0, 20.0)):
    """Rotationskurve allein aus dem Sichtbaren (Katalogabstand und -neigung, Υ_d aus der Eichung) gegen SPARC."""
    yd = D.sparc_eichung(K.kpc)
    k = D.sparc_kurve(name)
    vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + D.YB * k["Vbul"]**2
    aN = vb2 * 1e6 / (k["R"] * K.kpc)
    v = np.sqrt(np.array([K.verstaerkung(x) for x in aN]) * vb2)
    zeilen = []
    for r in radien:
        i = int(np.argmin(abs(k["R"] - r)))
        zeilen.append((k["R"][i], math.sqrt(vb2[i]), aN[i] / K.a0, v[i], k["V"][i], k["eV"][i]))
    rest = np.log10(v / k["V"])
    return dict(yd=yd, zeilen=zeilen, punkte=len(v), rms=float(np.sqrt(np.mean(rest**2))),
                verhaeltnis=float(np.mean(v / k["V"])))


def kosmos():
    """D_M/D_H aus DESI DR2 gegen das Modell: Ω_m folgt aus a₀ (Atmen), keine freie Zahl; r_d und H₀ kürzen sich."""
    from scipy.integrate import quad
    H0, Ode = K.atmen_heute()
    o_r = K.dichte_heute()[0] / (K.rho_krit100 * (H0 / 100)**2)
    om = 1 - Ode - o_r
    E = lambda z: math.sqrt(o_r * (1 + z)**4 + om * (1 + z)**3 + Ode)
    z, v, g, C = D.desi()
    zeilen, chi = [], 0.0
    for zz in sorted(set(z[np.array(g) == "DM_over_rs"])):
        i = [j for j in range(len(z)) if z[j] == zz and g[j] == "DM_over_rs"][0]
        k = [j for j in range(len(z)) if z[j] == zz and g[j] == "DH_over_rs"][0]
        f = v[i] / v[k]
        df = f * math.sqrt(C[i, i] / v[i]**2 + C[k, k] / v[k]**2 - 2 * C[i, k] / (v[i] * v[k]))
        m = quad(lambda x: 1 / E(x), 0, zz)[0] * E(zz)
        zeilen.append((zz, f, df, m))
        chi += ((f - m) / df)**2
    return dict(H0=H0, Ode=Ode, om=om, zeilen=zeilen, chi2=chi)


def sonnensystem():
    """Erde, Sonne, Erdbahn, Licht am Sonnenrand, Merkur und die Schwelle der Sonne: der Trichter ohne Rückkehr."""
    a_merkur, e_merkur = 0.387098 * K.AU, 0.205630
    umlauf = 2 * math.pi * math.sqrt(a_merkur**3 / K.GM_sun)
    return dict(erde=(K.faltung(K.M_earth) * 1e3, K.grenzradius(K.M_earth) * 1e3, (1 - K.takt(K.M_earth, K.R_ERDE)) * K.yr),
                sonne=(K.faltung(K.M_sun), K.grenzradius(K.M_sun) / 1e3, (1 - K.takt(K.M_sun, K.R_SONNE)) * K.yr),
                v_erde=K.bahntempo(K.M_sun, K.AU) / 1e3, T_erde=2 * math.pi * math.sqrt(K.AU**3 / K.GM_sun) / 86400,
                uhr_erdbahn=3 * K.faltung(K.M_sun) / K.AU * K.yr,
                licht=math.degrees(K.ablenkung(K.M_sun, K.R_SONNE)) * 3600,
                merkur=math.degrees(K.perihel_je_umlauf(K.M_sun, a_merkur, e_merkur)) * 3600 * 100 * K.yr / umlauf,
                r_stern=K.schwellenradius(K.M_sun) / K.AU)


def schliesszahlen():
    """Schließzahlen aus den Zählregeln (Hülle nach n + l, Kern nach O(k) + 2(k + 2)) gegen die stärksten Abfälle der
    Ionisierungsenergie (PubChem) und die stärksten Sprünge der Zwei-Nukleonen-Bindung (AME2020)."""
    import kerne as KE
    huelle = K.schliesszahlen_huelle(7)
    N, Z, atom, he, rn = D.schliesszahlen()
    ie = {v["Z"]: v["IE"] for v in D.pubchem().values() if v["IE"]}
    return dict(huelle=huelle, perioden=tuple(b - a for a, b in zip((0,) + huelle, huelle)),
                rein=tuple(sum(2 * n * n for n in range(1, k + 1)) for k in range(1, 5)),
                gehalten=tuple(sum((n + 1) * (n + 2) for n in range(k)) for k in range(1, 7)), kern=KE.schliesszahlen_kern()[3:8],
                N=N, Z=Z, atom=atom, he=he, rn=rn, abfall={z: ie[z] / ie[z + 1] for z in ie if z + 1 in ie})


def sterne():
    """Massenskala, Grenzmasse aus der Form n = 3, Decke des elektronengehaltenen Trichters, Ende der Brennleiter, Carter."""
    import kerne as KE
    xi1, om3 = K.lane_emden(3.0)
    carter, aG = K.alpha**12 * (K.m_e / K.m_p)**4 * 1e39, K.kopplung_grav(K.m_p) * 1e39
    return dict(M_stern=K.sternmasse() / K.M_sun, xi1=xi1, om3=om3, M_grenz=K.grenzmasse() / K.M_sun,
                decke=K.decke_entartet(1.0)[0] / (K.m_e * K.c**2), T_CO=K.decke_entartet(0.875)[0] / K.k_B / 1e8,
                Z_ende=KE.ende_brennleiter(KE.TROPFEN)[0], carter=carter, aG=aG, faktor=aG / carter)


def falte():
    """Die Falte: 10⁶ M☉ in der Galaxie, Temperatur und Periode der Sonnenfalte, Abgabe beim Fall, Entropie, Treffpunkt M_U/4."""
    M6 = 1e6 * K.M_sun
    fe = 100 * [e for e in D.ame2020() if e[0] == 26 and e[2] == 56][0][3] / (K.m_u * K.c**2 / K.e / 1e3)
    abgabe = 1 - 2 / 3 * math.sqrt(2)
    return dict(s2=K.c**4 / (4 * K.G * M6) / K.a0 / 1e17, r_stern=K.schwellenradius(M6) / K.pc,
                T=K.falte_temperatur(K.M_sun) / 1e-8, periode=K.falte_temperatur(K.M_sun) * K.c * K.k_B / (K.hbar * K.a0) / 1e22,
                abgabe=100 * abgabe, abgabe_j=100 * (1 - 1 / math.sqrt(3)), fe=fe, mal=abgabe / (fe / 100),
                S=K.entropie_falte(K.M_sun) / 1e77, MU4=K.schnittmassen()[3] / 4 / K.M_sun / 1e23)


def milchstrasse(radien=(5.0, None, 15.0, 50.0)):
    """Die Milchstraße als Punktmasse (6·10¹⁰ M☉): Regel gegen Newton, Licht weit außen, Naht."""
    M = 6e10 * K.M_sun
    rs = K.schwellenradius(M) / K.kpc
    r = [rs if x is None else x for x in radien]
    aN50 = K.G * M / (50 * K.kpc)**2
    return dict(ds=K.faltung(M) / K.AU, rs=rs, vflach=K.v_flach(M) / 1e3, nu_rs=K.verstaerkung(K.a0), r=r,
                v_N=[math.sqrt(K.G * M / (x * K.kpc)) / 1e3 for x in r], v=[K.bahntempo_regel(M, x * K.kpc) / 1e3 for x in r],
                nu50=K.verstaerkung(aN50),
                alpha_inf=math.degrees(K.ablenkung_grenz(M)) * 3600, RQ=K.naht_radius(M) / K.kpc,
                aN_Q=math.log(K.Q_kosm / (1 + K.Q_kosm))**2)


def schwelle_sparc():
    """Die Schwelle aus 141 Galaxien (SPARC): Υ_d aus der Eichung, a₀, die drei Formen, Dimension der Rückkehr, Scheiben an der
    Schwelle; Einzeltrichter und kollektive Systeme."""
    yd = D.sparc_eichung(K.kpc)
    chi, n_pkt = D.sparc_formen(K.kpc, yd)
    a0, (lo, hi) = D.sparc_schwelle(K.kpc, yd=yd)
    dim = D.sparc_dimension(K.kpc, yd=yd)
    sch = D.sparc_scheiben(K.kpc, yd, K.M_sun, K.a0, K.G)
    return dict(yd=yd, n_gal=len(D.sparc_auswahl(yd)), n_pkt=n_pkt, a0=a0 / 1e-10, lo=lo / 1e-10, hi=hi / 1e-10,
                d_einfach=chi["einfach"] - chi["Bose"], d_standard=chi["standard"] - chi["Bose"], dim=dim[0], dim_s=dim[1],
                Rd=sch[0], Rd_lo=sch[1][0], Rd_hi=sch[1][1], steig=sch[2], steig_s=sch[3],
                nu_kugel=K.verstaerkung(25 * K.a0), nu_duenn=K.verstaerkung(0.13 * K.a0),
                sig7=K.streuung_kugel(1e7 * K.M_sun) / 1e3, sig8=K.streuung_kugel(2e8 * K.M_sun) / 1e3,
                Sigma=K.flaechendichte_schwelle() * K.pc**2 / K.M_sun)


def leiter():
    """Aus a₀ allein: Schließlänge, Rückkehrzeit, Leiter der Längen, Atemzahl und Vakuumteil, Zelle."""
    Mmin, Mc, mPl, MU = K.schnittmassen()
    H = K.atmen_heute()[0] * 1e3 / K.Mpc
    return dict(Lam=K.schliesslaenge() / K.Gpc, tau=K.rueckkehrzeit() / K.Gyr, eps=K.leiterzahl() / 1e-62,
                Mmin=Mmin / 1e-70, Mc=Mc / 1e-29, mPl=mPl / 1e-8, MU=MU / 1e54, MU_sun=MU / K.M_sun / 1e23,
                theta=K.theta_(H=H), Os=K.vakuumteil(H=H), verh=(1 - K.Omega_b - K.Omega_c) / K.vakuumteil(),
                Ls=4 * math.pi**4 / 15 * K.leiterzahl()**2 / 1e-122, ks=math.pi**3 / 30,
                Hs=K.H0_kms * math.sqrt(K.vakuumteil()), H3=math.sqrt(8 * math.pi * K.G / 3 * K.dichte_dunkle_energie()) * K.Mpc / 1e3,
                form=K.form_zelle(), L=K.zelle_laenge() / K.Gpc)


def zahlen():
    """Jede gerechnete Zahl der Dateien rechnen_de/en.md, geschrieben mit Punkt (z wie im Katalog)."""
    a, g, k = atom_und_kern(), galaxie(), kosmos()
    g2 = galaxie("NGC3198")
    aus = [(a["E_H"], 3), (a["IE_H"], 3), (a["aC"], 2), (a["lam_H"], 2), (a["R_Fe"], 1),
           (a["v_groesse"], 2), (a["v_bindung"], 2), (a["lam_spiegel"], 2),
           (g["yd"], 2), (g["rms"], 2), (g["verhaeltnis"], 2), (g2["verhaeltnis"], 2),
           (k["H0"], 1), (k["Ode"], 3), (k["om"], 3), (k["chi2"], 1)]
    for r, vb, x, v, vo, ev in g["zeilen"]:
        aus += [(r, 2), (vb, 0), (x, 2), (v, 0), (vo, 0), (ev, 1)]
    for z, f, df, m in k["zeilen"]:
        aus += [(f, 3), (df, 3), (m, 3)]
    so = sonnensystem()
    aus += [(so["erde"][0], 3), (so["erde"][1], 2), (so["erde"][2], 3), (so["sonne"][0], 1), (so["sonne"][1], 2), (so["sonne"][2], 0),
            (so["v_erde"], 2), (so["T_erde"], 1), (so["uhr_erdbahn"], 2), (so["licht"], 4), (so["licht"] / 2, 4), (so["merkur"], 2),
            (so["r_stern"], 0)]
    sz = schliesszahlen()
    aus += [(sz["he"], 1), (sz["rn"], 1)]
    st = sterne()
    aus += [(st["M_stern"], 3), (st["xi1"], 3), (st["om3"], 3), (st["M_grenz"], 3), (st["decke"], 3), (st["T_CO"], 1),
            (st["Z_ende"], 1), (st["carter"], 1), (st["aG"], 3), (st["faktor"], 1)]
    f = falte()
    aus += [(f["s2"], 1), (f["r_stern"], 0), (f["T"], 1), (f["periode"], 1), (f["abgabe"], 1), (f["abgabe_j"], 0), (f["fe"], 2),
            (f["mal"], 1), (f["S"], 2), (f["MU4"], 1)]
    mw = milchstrasse()
    aus += [(mw["ds"], 1), (mw["rs"], 2), (mw["vflach"], 1), (mw["nu_rs"], 3), (mw["nu50"], 2),
            (mw["alpha_inf"], 2), (mw["RQ"], 1), (mw["aN_Q"], 3)] + [(v, 0) for v in mw["v_N"] + mw["v"]]
    sp = schwelle_sparc()
    aus += [(sp["yd"], 2), (sp["n_gal"], 0), (sp["n_pkt"], 0), (sp["a0"], 2), (sp["lo"], 2), (sp["hi"], 2), (sp["d_einfach"], 0),
            (sp["d_standard"], 0), (sp["dim"], 2), (sp["dim_s"], 2), (sp["Rd"], 2), (sp["Rd_lo"], 2), (sp["Rd_hi"], 2), (sp["steig"], 2),
            (sp["steig_s"], 2), (sp["nu_kugel"], 3), (sp["nu_duenn"], 1), (sp["sig7"], 1), (sp["sig8"], 1), (sp["Sigma"], 0)]
    le = leiter()
    aus += [(le["Lam"], 1), (le["tau"], 1), (le["eps"], 1), (le["Mmin"], 1), (le["Mc"], 1), (le["mPl"], 1), (le["MU"], 1), (le["MU_sun"], 1),
            (le["theta"], 3), (le["Os"], 3), (le["verh"], 2), (le["Ls"] * 10, 1), (le["ks"], 3), (le["Hs"], 1), (le["H3"], 1),
            (le["form"], 3), (le["L"], 1)]
    return [f"{v:.{n}f}" for v, n in aus] + [f"{z:g}" for z, *_ in k["zeilen"]]


# ---------------------------------------------------------------- Bilder zu den Beispielen
BESCHRIFTUNG = {
    "de": dict(aN="Newton aus dem Sichtbaren  a_N  [m/s²]", a="gemessene Beschleunigung  a  [m/s²]",
               regel="Regel  a = a_N/(1 − q)", newton="Newton  a = a_N", sparc="SPARC, 141 Galaxien",
               A="Massenzahl A", BA="Bindung je Nukleon [MeV]", ame="AME2020, Tal der Stabilität",
               tropfen="gehaltener Trichter (Tröpfchenform)", R="Radius R [kpc]", v="Geschwindigkeit v [km/s]",
               mess="gemessen (SPARC)", sichtbar="Sichtbares allein (Newton)", schwelle="a_N = a₀",
               z="Rotverschiebung z", DMDH="D_M/D_H", desi="DESI DR2", modell="Modell, keine freie Zahl",
               ohne="ohne Rückstellung (Ω_DE = 0)", Z="Kernladung Z", abfall="Abfall der Ionisierungsenergie  IE(Z)/IE(Z + 1)",
               pubchem="PubChem", vflach="v⁴ = GMa₀", masse="Masse M [kg]", laenge="Länge [m]", lam="λ̄ = ħ/(Mc)",
               ds="3Δs = GM/c²", rs="r* = √(GM/a₀)", l0="ℓ₀ = 2Λ", lP="ℓ_P", proton="Proton", sonne="Sonne", mw="Milchstraße"),
    "en": dict(aN="Newton from the visible  a_N  [m/s²]", a="measured acceleration  a  [m/s²]",
               regel="rule  a = a_N/(1 − q)", newton="Newton  a = a_N", sparc="SPARC, 141 galaxies",
               A="mass number A", BA="binding per nucleon [MeV]", ame="AME2020, valley of stability",
               tropfen="held funnel (droplet form)", R="radius R [kpc]", v="velocity v [km/s]",
               mess="measured (SPARC)", sichtbar="visible alone (Newton)", schwelle="a_N = a₀",
               z="redshift z", DMDH="D_M/D_H", desi="DESI DR2", modell="model, no free number",
               ohne="without restoring force (Ω_DE = 0)", Z="nuclear charge Z", abfall="drop of the ionisation energy  IE(Z)/IE(Z + 1)",
               pubchem="PubChem", vflach="v⁴ = GMa₀", masse="mass M [kg]", laenge="length [m]", lam="λ̄ = ħ/(Mc)",
               ds="3Δs = GM/c²", rs="r* = √(GM/a₀)", l0="ℓ₀ = 2Λ", lP="ℓ_P", proton="proton", sonne="Sun", mw="Milky Way")}


def _tex(t):
    """Beschriftung mit Index als Formelsatz, aufrecht wie im Text: a_N → $\\mathrm{a_N}$; m/s² bleibt."""
    import re
    return re.sub(r"\b([a-zA-Z]+|Ω)_([A-Za-z]+)\b", r"$\\mathrm{\1_{\2}}$", t).replace("Ω", "\\Omega")


def bilder(sprache="de", ordner=None):
    """Vier Bilder zu den Beispielen als SVG in abhandlung/bilder/<sprache>/: die Regel an allen SPARC-Punkten, die Bindung
    je Nukleon (AME2020) gegen den gehaltenen Trichter, die Rotationskurve von NGC 2403, D_M/D_H an DESI DR2.
    Nur Katalogwerte, keine Anpassung je Galaxie. Gibt die Pfade."""
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.integrate import quad
    import kerne as KE
    B = BESCHRIFTUNG[sprache]
    ordner = ordner or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "abhandlung", "bilder", sprache)
    os.makedirs(ordner, exist_ok=True)
    plt.rcParams.update({"font.family": "FreeSerif", "font.size": 9.5, "axes.linewidth": 0.6, "svg.fonttype": "path",
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                         "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.fontsize": 9,
                         "mathtext.fontset": "custom", "mathtext.rm": "FreeSerif", "mathtext.it": "FreeSerif:italic",
                         "mathtext.bf": "FreeSerif:bold", "axes.labelsize": 9.5})
    komma = FuncFormatter(lambda v, _: f"{v:g}".replace(".", "," if sprache == "de" else "."))
    wege = []

    def ablage(fig, name):
        weg = os.path.join(ordner, name + ".svg")
        fig.savefig(weg, bbox_inches="tight", pad_inches=0.04)
        plt.close(fig)
        wege.append(weg)

    # Bild 1: die Regel an allen Punkten der 141 Galaxien (Katalogabstand und -neigung, Υ_d aus der Eichung)
    yd = D.sparc_eichung(K.kpc)
    gb, go = [], []
    for n, p, k, eD, eI in D.sparc_auswahl(yd):
        R = k["R"] * K.kpc
        gb += list((k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + D.YB * k["Vbul"]**2) * 1e6 / R)
        go += list(k["V"]**2 * 1e6 / R)
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.loglog(gb, go, ".", color="#9a9a9a", ms=2.2, mew=0, alpha=0.7, label=B["sparc"], rasterized=False)
    x = np.logspace(-12.5, -8.3, 300)
    ax.loglog(x, x, "--", color="#555", lw=0.9, label=_tex(B["newton"]))
    ax.loglog(x, [K.verstaerkung(v) * v for v in x], "-", color="k", lw=1.4, label=_tex(B["regel"]))
    ax.axvline(K.a0, color="#aaa", lw=0.7, ls=":")
    ax.text(K.a0 * 1.15, 10**-12.3, "$\\mathrm{a_0}$", fontsize=10, ha="left", va="bottom")
    ax.set_xlim(10**-12.5, 10**-8.3); ax.set_ylim(10**-12.5, 10**-8.3)
    ax.set_xlabel(_tex(B["aN"])); ax.set_ylabel(_tex(B["a"]))
    ax.legend(loc="upper left")
    ablage(fig, "regel")

    # Bild 2: Bindung je Nukleon im Tal (AME2020) gegen den gehaltenen Trichter
    b = KE.b_gemessen()
    tal = {}
    for (z, a), bb in b.items():
        if bb / a > tal.get(a, (0, 0))[0]:
            tal[a] = (bb / a, z)
    A = np.array(sorted(a for a in tal if a >= 4)); BA = np.array([tal[a][0] for a in A])
    A2 = np.arange(12, 271)
    T = np.array([float(KE.bindung(KE.z_tal(a, KE.TROPFEN), a, KE.TROPFEN)[0]) / a for a in A2])
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.plot(A, BA, ".", color="#9a9a9a", ms=3.2, mew=0, label=B["ame"])
    ax.plot(A2, T, "-", color="k", lw=1.3, label=B["tropfen"])
    for a, name in ((56, "$^{56}$Fe"), (62, "$^{62}$Ni")):
        ax.annotate(name, (a, tal[a][0]), xytext=(a + (14 if a == 62 else -30), tal[a][0] + 0.35),
                    fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.5, color="#555"))
    ax.set_xlim(0, 275); ax.set_ylim(5, 9.3)
    ax.yaxis.set_major_formatter(komma)
    ax.set_xlabel(B["A"]); ax.set_ylabel(B["BA"])
    ax.legend(loc="lower right")
    ablage(fig, "kerne")

    # Bild 3: NGC 2403, nur das Sichtbare mit Katalogwerten
    k = D.sparc_kurve("NGC2403")
    vb2 = k["Vgas"] * np.abs(k["Vgas"]) + yd * k["Vdisk"]**2 + D.YB * k["Vbul"]**2
    aN = vb2 * 1e6 / (k["R"] * K.kpc)
    v = np.sqrt(np.array([K.verstaerkung(x) for x in aN]) * vb2)
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.errorbar(k["R"], k["V"], yerr=k["eV"], fmt="o", color="k", ms=2.8, lw=0.7, capsize=0, label=B["mess"])
    ax.plot(k["R"], np.sqrt(vb2), "--", color="#555", lw=1.0, label=B["sichtbar"])
    ax.plot(k["R"], v, "-", color="k", lw=1.4, label=B["regel"].split("  ")[0])
    i = int(np.max(np.where(aN > K.a0)[0]))                     # äußerer Schnitt mit der Schwelle
    ax.axvline(k["R"][i], color="#aaa", lw=0.7, ls=":")
    ax.text(k["R"][i] + 0.3, 142, _tex(B["schwelle"]), fontsize=9, color="#555", va="top")
    ax.set_xlim(0, 21.5); ax.set_ylim(0, 150)
    ax.xaxis.set_major_formatter(komma)
    ax.set_xlabel(B["R"]); ax.set_ylabel(B["v"])
    ax.legend(loc="lower right")
    ablage(fig, "ngc2403")

    # Bild 4: D_M/D_H an DESI DR2 gegen das Modell, daneben die Ausdehnung ohne Rückstellung
    kos = kosmos()
    o_r = K.dichte_heute()[0] / (K.rho_krit100 * (kos["H0"] / 100)**2)
    def kurve(om, ode):
        E = lambda z: math.sqrt(o_r * (1 + z)**4 + om * (1 + z)**3 + ode)
        return lambda z: quad(lambda x: 1 / E(x), 0, z)[0] * E(z)
    z = np.linspace(0.01, 2.6, 120)
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.plot(z, [kurve(1 - o_r, 0.0)(x) for x in z], "--", color="#555", lw=0.9, label=_tex(B["ohne"]))
    ax.plot(z, [kurve(kos["om"], kos["Ode"])(x) for x in z], "-", color="k", lw=1.4, label=B["modell"])
    zz, f, df, m = zip(*kos["zeilen"])
    ax.errorbar(zz, f, yerr=df, fmt="o", color="k", ms=3.2, lw=0.8, capsize=2, label=B["desi"])
    ax.set_xlim(0, 2.6); ax.set_ylim(0, 6.5)
    ax.xaxis.set_major_formatter(komma); ax.yaxis.set_major_formatter(komma)
    ax.set_xlabel(B["z"]); ax.set_ylabel(_tex(B["DMDH"]))
    ax.legend(loc="upper left")
    ablage(fig, "desi")

    # Bild 5: die Milchstraße als Punktmasse, Newton gegen die Regel
    mw = milchstrasse()
    M = 6e10 * K.M_sun
    r = np.logspace(0, 2, 200)
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.semilogx(r, np.sqrt(K.G * M / (r * K.kpc)) / 1e3, "--", color="#555", lw=0.9, label=B["newton"].split("  ")[0])
    ax.semilogx(r, [K.bahntempo_regel(M, x * K.kpc) / 1e3 for x in r], "-", color="k", lw=1.4, label=B["regel"].split("  ")[0])
    ax.axhline(mw["vflach"], color="#aaa", lw=0.7, ls=":")
    ax.text(1.1, mw["vflach"] + 6, B["vflach"], fontsize=9, color="#555")
    ax.axvline(mw["rs"], color="#aaa", lw=0.7, ls=":")
    ax.text(mw["rs"] * 1.06, 30, "r*", fontsize=10, color="#555")
    ax.plot(mw["r"], mw["v"], "o", color="k", ms=3.5)
    for x, v in zip(mw["r"], mw["v"]):
        ax.annotate(f"{v:.0f}", (x, v), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8.5)
    ax.set_xlim(1, 100); ax.set_ylim(0, 520)
    ax.xaxis.set_major_formatter(komma)
    ax.set_xlabel(B["R"]); ax.set_ylabel(B["v"])
    ax.legend(loc="upper right")
    ablage(fig, "milchstrasse")

    # Bild 6: die Abfälle der Ionisierungsenergie (PubChem) und die Schließzahlen der Hülle
    sz = schliesszahlen()
    sym = {v["Z"]: s for s, v in D.pubchem().items()}
    zs = sorted(sz["abfall"])
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    ax.bar(zs, [sz["abfall"][z] for z in zs], width=0.8, color="#bbb", lw=0)
    ax.bar(list(sz["atom"]), [sz["abfall"][z] for z in sz["atom"]], width=0.8, color="k", lw=0)
    for z in sz["atom"]:
        ax.annotate(f"{sym[z]} {z}", (z, sz["abfall"][z]), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8.5)
    ax.axhline(1, color="#aaa", lw=0.6)
    ax.set_xlim(0, 119); ax.set_ylim(0, 5.6)
    ax.yaxis.set_major_formatter(komma)
    ax.set_xlabel(B["Z"]); ax.set_ylabel(B["abfall"].split("  ")[1])
    ax.text(118, 5.2, B["pubchem"], ha="right", fontsize=9, color="#555")
    ablage(fig, "schliesszahlen")

    # Bild 7: die Leiter der Längen
    from matplotlib.ticker import LogLocator
    Mmin, Mc, mPl, MU = K.schnittmassen()
    M = np.logspace(math.log10(Mmin), math.log10(MU), 400)
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.loglog(M, K.hbar / (M * K.c), "-", color="k", lw=1.0, label=B["lam"])
    ax.loglog(M, K.G * M / K.c**2, "-", color="#555", lw=1.0, label=B["ds"])
    ax.loglog(M, np.sqrt(K.G * M / K.a0), "-", color="#999", lw=1.2, label=B["rs"])
    l0 = 2 * K.schliesslaenge()
    ax.axhline(l0, color="k", lw=0.7, ls="--"); ax.text(mPl, l0 * 1e3, B["l0"], fontsize=9, ha="center")
    ax.axhline(K.l_P, color="k", lw=0.5, ls=":"); ax.text(Mmin * 1e3, K.l_P * 1e3, "$\\ell_{\\mathrm{P}}$", fontsize=10)
    for m, y, name in ((Mmin, l0, "M_min"), (Mc, K.hbar / (Mc * K.c), "M_c"), (mPl, K.l_P, "m_Pl"), (MU, l0, "M_U")):
        ax.plot(m, y, "o", color="k", ms=4)
        ax.annotate(_tex(name), (m, y), xytext=(6, -13) if name != "M_min" else (6, 4), textcoords="offset points", fontsize=9)
    for m, name in ((K.m_p, B["proton"]), (K.M_sun, B["sonne"]), (6e10 * K.M_sun, B["mw"])):
        ax.axvline(m, color="#ccc", lw=0.6)
        ax.text(m * 2, 1e-90, name, fontsize=8, rotation=90, va="bottom", color="#555")
    ax.set_xlim(Mmin / 10, MU * 10); ax.set_ylim(1e-100, 1e30)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=8)); ax.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
    ax.set_xlabel(B["masse"]); ax.set_ylabel(B["laenge"])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8.5)
    ablage(fig, "leiter")
    return wege


if __name__ == "__main__":
    import pprint
    pprint.pprint((atom_und_kern(), galaxie(), kosmos()))
    print(bilder())
