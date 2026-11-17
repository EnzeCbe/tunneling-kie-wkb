"""
Figure 3: Hydrogen vs. deuterium tunneling and the resulting kinetic isotope
effect (KIE), computed via WKB transmission through an Eckart barrier.

KIE here is defined at fixed energy E as the ratio of transmission
probabilities, P_H(E) / P_D(E), which is the dominant tunneling contribution
to the overall rate-constant KIE k_H/k_D (the full thermal KIE in Figure 4
additionally folds in the Boltzmann-averaged rate constants).
"""

import numpy as np
import matplotlib.pyplot as plt
from wkb_core import (
    EV_TO_J, M_PROTON, M_DEUTERON,
    solve_V2_for_barrier_height, wkb_transmission,
)

# --- Barrier parameters ------------------------------------------------
V1 = 0.0
V_MAX_EV = 0.8
E_EV = 0.4   # sub-barrier energy, same for H and D
V2 = solve_V2_for_barrier_height(V_MAX_EV * EV_TO_J, V1)
E = E_EV * EV_TO_J

L_values_angstrom = np.linspace(0.1, 1.0, 80)
L_values = L_values_angstrom * 1e-10

T_H = np.array([wkb_transmission(E, V1, V2, L, M_PROTON) for L in L_values])
T_D = np.array([wkb_transmission(E, V1, V2, L, M_DEUTERON) for L in L_values])
KIE = T_H / T_D

# --- Plot ----------------------------------------------------------------
plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.0})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=300)

ax1.plot(L_values_angstrom, T_H, color="#1f77b4", linewidth=2.2, label="Hydrogen (H)")
ax1.plot(L_values_angstrom, T_D, color="#d62728", linewidth=2.2, label="Deuterium (D)")
ax1.set_xlabel("Eckart barrier width parameter, $L$ (\u00c5)", fontsize=13)
ax1.set_ylabel("WKB transmission probability, $P(E)$", fontsize=13)
ax1.set_title("H vs. D Tunneling", fontsize=14)
ax1.set_yscale("log")
ax1.legend(fontsize=12, frameon=True)
ax1.grid(True, which="both", alpha=0.3)

ax2.plot(L_values_angstrom, KIE, color="#2ca02c", linewidth=2.2)
ax2.set_xlabel("Eckart barrier width parameter, $L$ (\u00c5)", fontsize=13)
ax2.set_ylabel(r"Transmission ratio  ($P_H / P_D$)", fontsize=13)
ax2.set_title("H/D Transmission Ratio at One Energy", fontsize=14)
ax2.set_yscale("log")
ax2.grid(True, which="both", alpha=0.3)

fig.suptitle(f"$V_\\mathrm{{max}}$ = {V_MAX_EV} eV, $E$ = {E_EV} eV", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig("figure3_KIE.png", dpi=300, bbox_inches="tight")
# plt.show()  # uncomment to open the figure window

print(f"P_H/P_D at L=0.3A: {np.interp(0.3, L_values_angstrom, KIE):.2f}")
print(f"P_H/P_D at L=0.6A: {np.interp(0.6, L_values_angstrom, KIE):.2f}")
print("(Compare: experimental enzymatic KIEs range from ~7 (semiclassical limit) "
      "up to about 80 for extreme cases like soybean lipoxygenase.)")
