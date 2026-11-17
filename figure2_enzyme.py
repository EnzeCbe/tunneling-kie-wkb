"""
Figure 2: Effect of the enzyme active site on proton tunneling probability.

Physical picture: the enzyme's role is modeled as modifying the barrier
geometry along the C-H reaction coordinate -- specifically, lowering the
barrier height and/or narrowing its width relative to the uncatalyzed
(solution-phase) reaction. This is a simplified but standard way to
represent pre-organization and electrostatic stabilization of the
transition state by the protein environment (Klinman & Kohen 2013).
"""

import numpy as np
import matplotlib.pyplot as plt
from wkb_core import (
    EV_TO_J, M_PROTON,
    solve_V2_for_barrier_height, wkb_transmission,
)

# --- Barrier parameters: enzyme lowers barrier height AND narrows width ----
V1 = 0.0

V_MAX_NO_ENZYME_EV = 1.2    # higher, uncatalyzed barrier
V_MAX_ENZYME_EV = 0.6       # enzyme-stabilized, lower barrier

E_EV = 0.4                  # particle energy (eV) -- same for both cases

V2_no_enzyme = solve_V2_for_barrier_height(V_MAX_NO_ENZYME_EV * EV_TO_J, V1)
V2_enzyme = solve_V2_for_barrier_height(V_MAX_ENZYME_EV * EV_TO_J, V1)
E = E_EV * EV_TO_J

L_values_angstrom = np.linspace(0.1, 1.5, 80)
L_values = L_values_angstrom * 1e-10

T_no_enzyme = np.array([wkb_transmission(E, V1, V2_no_enzyme, L, M_PROTON) for L in L_values])
T_with_enzyme = np.array([wkb_transmission(E, V1, V2_enzyme, L, M_PROTON) for L in L_values])

# --- Plot ---------------------------------------------------------------
plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.0})

fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
ax.plot(L_values_angstrom, T_no_enzyme, color="#d62728", linewidth=2.2, label=f"Higher barrier ({V_MAX_NO_ENZYME_EV} eV)")
ax.plot(L_values_angstrom, T_with_enzyme, color="#2ca02c", linewidth=2.2, label=f"Lower barrier ({V_MAX_ENZYME_EV} eV)")

ax.set_xlabel("Eckart barrier width parameter, $L$ (\u00c5)", fontsize=13)
ax.set_ylabel("WKB transmission probability, $P(E)$", fontsize=13)
ax.set_title("Effect of Barrier Height on Proton Tunneling\n"
             f"($V_\\mathrm{{max}}$: {V_MAX_NO_ENZYME_EV} eV \u2192 {V_MAX_ENZYME_EV} eV, $E$ = {E_EV} eV)",
             fontsize=13)
ax.set_yscale("log")
ax.legend(fontsize=12, frameon=True)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("figure2_enzyme.png", dpi=300, bbox_inches="tight")
# plt.show()  # uncomment to open the figure window

enhancement = T_with_enzyme / T_no_enzyme
print(f"Tunneling enhancement factor at L=0.5A: {np.interp(0.5, L_values_angstrom, enhancement):.3e}")
