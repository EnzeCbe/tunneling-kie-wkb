"""Rate constants and KIEs for scans over barrier height and width.

Two versions of the model are used:
  * tunneling only: H and D see the same Eckart barrier;
  * ZPE in barrier: each isotope sees an effective barrier lowered by its own
    C-H or C-D zero-point energy, Vmax - (1/2) h c nu, and WKB is applied to it.
"""
import numpy as np
from wkb_core import (EV_TO_J, M_PROTON, M_DEUTERON, solve_V2_for_barrier_height,
                      eckart_barrier_height, wkb_transmission)

KB = 1.380649e-23
H = 6.62607015e-34
C_CM = 2.99792458e10          # speed of light in cm/s
NU_H, NU_D = 3000.0, 2250.0   # cm^-1, typical C-H / C-D stretch (Schaller, Reactivity in Chemistry)
ZPE = {M_PROTON: 0.5 * H * C_CM * NU_H, M_DEUTERON: 0.5 * H * C_CM * NU_D}   # J
A = 1e13                      # s^-1, same for both isotopes
N_E = 800                     # energy points (KIE converged to <1%)


def _transmission_table(Vmax_J, L, mass):
    V1 = 0.0
    V2 = solve_V2_for_barrier_height(Vmax_J, V1)
    Vm = eckart_barrier_height(V1, V2)
    E = np.linspace(1e-6 * Vm, 3 * Vm, N_E)
    P = np.array([wkb_transmission(e, V1, V2, L, mass) for e in E])
    return E, P


def rates(Vmax_eV, L_angstrom, T, zpe=False):
    """Return k_H(T), k_D(T) arrays for a symmetric Eckart barrier."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    L = L_angstrom * 1e-10
    out = []
    for m in (M_PROTON, M_DEUTERON):
        Vm = Vmax_eV * EV_TO_J - (ZPE[m] if zpe else 0.0)
        E, P = _transmission_table(Vm, L, m)
        k = np.array([A * np.trapezoid(P * np.exp(-E / (KB * t)), E) / (KB * t) for t in T])
        out.append(k)
    return out[0], out[1]


def arrhenius_fit(T, k):
    """Fit ln k = ln A - Ea/(kB T); return (Ea in kcal/mol, A)."""
    slope, intercept = np.polyfit(1 / T, np.log(k), 1)
    return -slope * KB * 6.02214076e23 / 4184, np.exp(intercept)


def fwhm_angstrom(L_angstrom):
    """Full width at half maximum of V = Vmax sech^2(x/2L)."""
    return 4 * np.arccosh(np.sqrt(2)) * L_angstrom
