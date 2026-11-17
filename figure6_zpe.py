"""Figure 6 and Table 1: Arrhenius plot and KIE(T) with ZPE in the barrier."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wkb_core import EV_TO_J, M_PROTON, M_DEUTERON
from rate_sweep import rates, arrhenius_fit, fwhm_angstrom, ZPE, KB, A

T = np.linspace(200, 400, 41)
cases = [(0.9, 0.40), (1.1, 0.30), (1.5, 0.25)]
res = {c: rates(*c, T, zpe=True) for c in cases}
dZPE = ZPE[M_PROTON] - ZPE[M_DEUTERON]

plt.rcParams.update({"font.size": 12})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5), dpi=300)
kH, kD = res[(0.9, 0.40)]
kH_cl = A * np.exp(-(0.9 * EV_TO_J - ZPE[M_PROTON]) / (KB * T))
kD_cl = A * np.exp(-(0.9 * EV_TO_J - ZPE[M_DEUTERON]) / (KB * T))
a1.plot(1000 / T, np.log(kH_cl), "--", color="#1f77b4", lw=1.4, alpha=.7, label="H, no tunneling")
a1.plot(1000 / T, np.log(kD_cl), "--", color="#d62728", lw=1.4, alpha=.7, label="D, no tunneling")
a1.plot(1000 / T, np.log(kH), color="#1f77b4", lw=2.2, label="H, WKB")
a1.plot(1000 / T, np.log(kD), color="#d62728", lw=2.2, label="D, WKB")
a1.set_xlabel(r"$1000/T$ (K$^{-1}$)"); a1.set_ylabel(r"$\ln k$"); a1.invert_xaxis(); a1.grid(alpha=.3); a1.legend()
a1.set_title(r"$V_{max}$ = 0.9 eV, $L$ = 0.40 Å, ZPE in barrier")
for (V, LA), (kH, kD) in res.items():
    a2.plot(T, kH / kD, lw=2.2, label=f"{V} eV, {LA:.2f} Å")
a2.plot(T, np.exp(dZPE / (KB * T)), "--", color="gray", lw=1.4, label="ZPE only, no tunneling")
a2.set_yscale("log"); a2.set_xlabel("Temperature (K)"); a2.set_ylabel(r"KIE ($k_H/k_D$)")
a2.grid(True, which="both", alpha=.3); a2.legend(fontsize=10); a2.set_title("KIE with ZPE in the barrier")
fig.tight_layout(); fig.savefig("figure6_zpe_arrhenius.png", bbox_inches="tight")

i300, i220 = np.argmin(abs(T - 300)), np.argmin(abs(T - 220))
sel = (T >= 280) & (T <= 320)
print("dZPE = %.2f kcal/mol; ZPE-only KIE at 300 K = %.2f, at 220 K = %.2f" % (
    dZPE * 6.02214076e23 / 4184, np.exp(dZPE / (KB * 300)), np.exp(dZPE / (KB * 220))))
print("Table 1: Vmax(eV) L(A) FWHM(A) KIE300 KIE220 dEa(kcal/mol) AH/AD")
for (V, LA), (kH, kD) in res.items():
    EaH, AH = arrhenius_fit(T[sel], kH[sel]); EaD, AD = arrhenius_fit(T[sel], kD[sel])
    print("%.1f %.2f %.2f %.3g %.3g %.2f %.3g" % (V, LA, fwhm_angstrom(LA), kH[i300] / kD[i300], kH[i220] / kD[i220], EaD - EaH, AH / AD))
EaH, AH = arrhenius_fit(T[sel], kH_cl[sel]); EaD, AD = arrhenius_fit(T[sel], kD_cl[sel])
print("no tunneling, 0.9 eV: dEa %.2f AH/AD %.3g" % (EaD - EaH, AH / AD))
