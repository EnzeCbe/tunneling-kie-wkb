import numpy as np
import matplotlib.pyplot as plt
from wkb_core import (
    EV_TO_J, M_PROTON, M_DEUTERON,
    eckart_barrier_height, solve_V2_for_barrier_height, thermal_rate_constant,
)

KB = 1.380649e-23

V1 = 0.0
V_MAX_EV = 1.5
V2 = solve_V2_for_barrier_height(V_MAX_EV * EV_TO_J, V1)
L = 0.25e-10
A = 1e13
Ea = eckart_barrier_height(V1, V2)

T_values = np.linspace(200, 400, 41)   # 5 K steps, includes 220 K and 300 K

k_classical = A * np.exp(-Ea / (KB * T_values))
k_H = np.array([thermal_rate_constant(T, V1, V2, L, M_PROTON, A) for T in T_values])
k_D = np.array([thermal_rate_constant(T, V1, V2, L, M_DEUTERON, A) for T in T_values])
KIE_T = k_H / k_D

inv_T = 1000 / T_values

plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.0})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=300)

ax1.plot(inv_T, np.log(k_classical), color="gray", linewidth=2.0, linestyle="--", label="Classical (no tunneling)")
ax1.plot(inv_T, np.log(k_H), color="#1f77b4", linewidth=2.2, label="WKB\u2013H")
ax1.plot(inv_T, np.log(k_D), color="#d62728", linewidth=2.2, label="WKB\u2013D")
ax1.set_xlabel(r"$1000\,/\,T$  (K$^{-1}$)", fontsize=13)
ax1.set_ylabel(r"$\ln(k)$", fontsize=13)
ax1.set_title("Tunneling-Corrected Arrhenius Plot", fontsize=13)
ax1.legend(fontsize=11, frameon=True)
ax1.grid(True, alpha=0.3)
ax1.invert_xaxis()

ax2.plot(T_values, KIE_T, color="#2ca02c", linewidth=2.2)
ax2.axhline(7, color="gray", linestyle="--", linewidth=1.4, label="Semiclassical limit (~7)")
ax2.set_xlabel("Temperature, $T$ (K)", fontsize=13)
ax2.set_ylabel(r"KIE  ($k_\mathrm{H}\,/\,k_\mathrm{D}$)", fontsize=13)
ax2.set_title("Temperature Dependence of KIE", fontsize=13)
ax2.set_yscale("log")
ax2.legend(fontsize=11, frameon=True)
ax2.grid(True, which="both", alpha=0.3)

fig.suptitle(r"Eckart barrier: $V_\mathrm{max}$ = 1.5 eV, $L$ = 0.25 Å", fontsize=12)
fig.tight_layout()
fig.savefig("figure4_arrhenius_final.png", dpi=300, bbox_inches="tight")
# plt.show()  # uncomment to open the figure window
print(f"KIE at 300K: {np.interp(300, T_values, KIE_T):.1f}")
print(f"KIE at 220K: {np.interp(220, T_values, KIE_T):.1f}")
