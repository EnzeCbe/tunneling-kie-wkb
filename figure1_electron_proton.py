"""
Figure 1: Electron vs. proton tunneling probability through an Eckart barrier,
as a function of barrier width.

This replaces the square-barrier analytic shortcut (T = exp(-2*kappa*L)) with
a full numerical WKB transmission probability through an asymmetric Eckart
barrier, evaluated at a fixed sub-barrier energy E.
"""

import numpy as np
import matplotlib.pyplot as plt
from wkb_core import (
    EV_TO_J, M_ELECTRON, M_PROTON,
    eckart_barrier_height, solve_V2_for_barrier_height, wkb_transmission,
)

# --- Barrier and energy parameters -----------------------------------------
V1 = 0.0                                  # symmetric barrier
V_MAX_EV = 1.0                            # barrier height (eV)
E_EV = 0.5                                # particle energy (eV), E < V_max
V2 = solve_V2_for_barrier_height(V_MAX_EV * EV_TO_J, V1)
E = E_EV * EV_TO_J

# Barrier width parameter L swept over a physically relevant range (Angstrom)
L_values_angstrom = np.linspace(0.1, 1.5, 80)
L_values = L_values_angstrom * 1e-10

# --- Compute WKB transmission probability at each width ---------------------
T_electron = np.array([wkb_transmission(E, V1, V2, L, M_ELECTRON) for L in L_values])
T_proton = np.array([wkb_transmission(E, V1, V2, L, M_PROTON) for L in L_values])

# --- Plot --------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 12,
    "axes.linewidth": 1.0,
})

fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
ax.plot(L_values_angstrom, T_electron, color="#1f77b4", linewidth=2.2, label="Electron")
ax.plot(L_values_angstrom, T_proton, color="#d62728", linewidth=2.2, label="Proton (C\u2013H coordinate)")

ax.set_xlabel("Eckart barrier width parameter, $L$ (\u00c5)", fontsize=13)
ax.set_ylabel("WKB transmission probability, $P(E)$", fontsize=13)
ax.set_title(f"Electron vs. Proton Tunneling Through an Eckart Barrier\n"
             f"($V_\\mathrm{{max}}$ = {V_MAX_EV} eV, $E$ = {E_EV} eV)", fontsize=13)
ax.set_yscale("log")
ax.legend(fontsize=12, frameon=True)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("figure1_tunneling.png", dpi=300, bbox_inches="tight")
# plt.show()  # uncomment to open the figure window

print(f"Barrier height check: {eckart_barrier_height(V1, V2)/EV_TO_J:.4f} eV (target {V_MAX_EV} eV)")
print(f"P(electron) at L=0.5A: {np.interp(0.5, L_values_angstrom, T_electron):.3e}")
print(f"P(proton)   at L=0.5A: {np.interp(0.5, L_values_angstrom, T_proton):.3e}")
