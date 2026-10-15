# B10: Modellwert cH0/2π gegen gemessenes a0 = 1,2e-10 – Rolle des Masse-zu-Licht-Verhältnisses
import numpy as np, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from scipy.optimize import minimize_scalar
c=2.99792458e8; Mpc=3.0857e22
nu=lambda y: 1/(1-np.exp(-np.sqrt(y)))
gobs=lambda gb,a0: gb*nu(gb/a0)
out=[]; P=lambda s="": out.append(s)
a0m=1.2e-10
P("1) Modellwert a0 = cH0/2π und das zum Messwert gehörige H0")
for H in [67.4,70.0,73.0]:
    a=c*H*1e3/Mpc/(2*np.pi); P(f"   H0 = {H:5.1f} km/s/Mpc:  a0 = {a:.4e} m/s²   Abstand zu 1,2e-10 bei σ_sys = 0,24e-10: {(a0m-a)/0.24e-10:.2f} σ")
P(f"   a0 = 1,2e-10 verlangt H0 = {2*np.pi*a0m/c*Mpc/1e3:.1f} km/s/Mpc")
P()
P("2) Stellarer Anteil: Wahres g_bar = f_*·k·g + (1−f_*)·g (Sterne mit Faktor k skaliert, Gas fest).")
P("   Gesucht k, bei dem die Anpassung mit Υ_ann (k = 1) das Modell-a0 als 1,2e-10 erscheinen lässt.")
gb=np.logspace(-12,-9,200)
for H in [67.4,70.0]:
    a0t=c*H*1e3/Mpc/(2*np.pi)
    for fs in [1.0,0.7,0.5]:
        def fit_a0(k):
            gtrue=gb*(fs*k+(1-fs)); go=gobs(gtrue,a0t)
            r=minimize_scalar(lambda la: np.sum((np.log10(gobs(gb,10**la))-np.log10(go))**2),bounds=(-11,-9),method='bounded')
            return 10**r.x
        ks=np.linspace(0.7,1.3,601); af=np.array([fit_a0(k) for k in ks])
        k=ks[np.argmin(abs(af-a0m))]
        P(f"   H0 = {H}: Sternanteil {fs:.1f}: a0_fit = 1,2e-10 bei k = {k:.3f}  → Υ_ann = 0,50 entspricht Υ_wahr = {0.5*k:.3f}")
open("schwelle_modell_mess.txt","w",encoding="utf-8").write("\n".join(out)); print("\n".join(out))

# 3) Äußeres Feld: Daten mit EFE erzeugt, ohne EFE angepasst
out=[]
P("3) Äußeres Feld g_e: g = ν(√(g_N²+g_e²)/a0)·g_N erzeugt, mit ν(g_N/a0)·g_N angepasst (g_N = 1e-12 … 1e-9)")
for H in [67.4,70.0]:
    a0t=c*H*1e3/Mpc/(2*np.pi)
    for ge in [0.0,0.01,0.03,0.1]:
        go=gb*nu(np.hypot(gb,ge*a0t)/a0t)
        r=minimize_scalar(lambda la: np.sum((np.log10(gobs(gb,10**la))-np.log10(go))**2),bounds=(-11,-9),method='bounded')
        P(f"   H0 = {H}: g_e = {ge:4.2f} a0 → angepasstes a0 = {10**r.x:.4e}  ({10**r.x/a0t:.3f} × Modellwert)")
open("schwelle_modell_mess.txt","a",encoding="utf-8").write("\n\n"+"\n".join(out)); print("\n".join(out))
