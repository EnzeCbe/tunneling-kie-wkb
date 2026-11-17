"""
Core physics module for WKB tunneling calculations through an Eckart barrier.

This module replaces the square-barrier analytic shortcut (T = exp(-2*kappa*L))
used in earlier drafts with a proper numerical WKB treatment of an asymmetric
Eckart barrier, which is a standard model for reaction-coordinate barriers in
enzymatic C-H activation (see e.g. Klinman & Kohen, Annu. Rev. Biochem. 2013).

All energies are handled internally in Joules; barrier heights/well depths are
supplied in eV for convenience and converted at the function boundary.
"""

import numpy as np
from scipy import integrate, optimize

# ---------------------------------------------------------------------------
# Physical constants (SI units)
# ---------------------------------------------------------------------------
HBAR = 1.0545718e-34      # J*s
EV_TO_J = 1.602176634e-19 # J per eV

# Particle masses (kg)
M_ELECTRON = 9.10938e-31
M_PROTON = 1.67262e-27
M_DEUTERON = 3.34358e-27   # ~2x proton mass (deuteron, for D KIE)


def eckart_potential(x, V1, V2, L):
    """
    Eckart barrier potential.

    V(x) = V1 * y / (1 + y) + V2 * y / (1 + y)^2,   y = exp(x / L)

    Parameters
    ----------
    x : array_like
        Reaction coordinate (meters), defined on (-inf, inf).
    V1 : float
        Asymptotic energy difference V(+inf) - V(-inf) (Joules). Sets the
        reaction exothermicity/endothermicity (0 = symmetric barrier).
    V2 : float
        Controls the barrier height above the higher asymptote (Joules).
        Must satisfy V2 > 0 and, for a physical barrier, V2 > -V1 as well.
    L : float
        Barrier width parameter (meters).

    Returns
    -------
    V : ndarray
        Potential energy at each x (Joules).
    """
    y = np.exp(np.clip(x / L, -700, 700))  # avoid overflow in exp
    return V1 * y / (1 + y) + V2 * y / (1 + y) ** 2


def eckart_barrier_height(V1, V2):
    """
    Maximum height of the Eckart barrier above the reactant (x -> -inf) well.

    Standard result: V_max = (V1 + V2)^2 / (4 * V2)  [for V2 > 0]
    """
    return (V1 + V2) ** 2 / (4 * V2)


def solve_V2_for_barrier_height(V_max, V1):
    """
    Given a desired barrier height V_max and asymmetry V1, solve for the
    shape parameter V2 such that eckart_barrier_height(V1, V2) == V_max.

    From V_max = (V1+V2)^2 / (4*V2), rearranging as a quadratic in V2:
        V2^2 + (2*V1 - 4*V_max)*V2 + V1^2 = 0
    We take the positive root, which is the physically relevant branch.
    """
    a = 1.0
    b = 2 * V1 - 4 * V_max
    c = V1 ** 2
    disc = b ** 2 - 4 * a * c
    if disc < 0:
        raise ValueError("No real solution for V2 with these V_max, V1")
    root1 = (-b + np.sqrt(disc)) / (2 * a)
    root2 = (-b - np.sqrt(disc)) / (2 * a)
    candidates = [r for r in (root1, root2) if r > 0]
    if not candidates:
        raise ValueError("No positive V2 solution found")
    return max(candidates)


def find_turning_points(E, V1, V2, L, x_range=None, n_scan=4000):
    """
    Find the classical turning points x1 < x2 where V(x) = E, i.e. the
    boundaries of the classically forbidden (tunneling) region for energy E.

    Uses a coarse scan to bracket sign changes of (V(x) - E), then refines
    each bracket with Brent's method.
    """
    if x_range is None:
        x_range = (-50 * L, 50 * L)

    x_scan = np.linspace(x_range[0], x_range[1], n_scan)
    f = eckart_potential(x_scan, V1, V2, L) - E

    sign_changes = np.where(np.diff(np.sign(f)) != 0)[0]
    if len(sign_changes) < 2:
        # No classically forbidden region at this energy (E above barrier top,
        # or scan range too narrow). Caller should handle this case.
        return None

    def root_fn(x):
        return eckart_potential(x, V1, V2, L) - E

    x1 = optimize.brentq(root_fn, x_scan[sign_changes[0]], x_scan[sign_changes[0] + 1])
    x2 = optimize.brentq(root_fn, x_scan[sign_changes[-1]], x_scan[sign_changes[-1] + 1])
    return x1, x2


def wkb_transmission(E, V1, V2, L, mass):
    """
    WKB tunneling probability through the Eckart barrier at energy E.

    P(E) = exp( -2 * Integral_{x1}^{x2} sqrt(2*mass*(V(x)-E)) / hbar  dx )

    For E at or above the barrier maximum, returns 1.0 (classical
    over-the-barrier passage; no tunneling integral needed).
    """
    V_max = eckart_barrier_height(V1, V2)
    if E >= V_max:
        return 1.0

    turning_points = find_turning_points(E, V1, V2, L)
    if turning_points is None:
        return 1.0

    x1, x2 = turning_points

    def integrand(x):
        V = eckart_potential(x, V1, V2, L)
        diff = np.maximum(V - E, 0.0)  # guard against tiny negative roundoff
        return np.sqrt(2 * mass * diff) / HBAR

    integral, _ = integrate.quad(integrand, x1, x2, limit=200)
    return np.exp(-2 * integral)


def thermal_rate_constant(T, V1, V2, L, mass, A_prefactor, kB=1.380649e-23,
                           E_max_factor=3.0, n_E=800):
    """
    Tunneling-corrected rate constant at temperature T, obtained by averaging
    the WKB transmission probability over a Boltzmann energy distribution:

        k(T) = A * Integral_0^Emax  P(E) * exp(-E / kB T) dE  /  (kB T)

    The division by kB*T normalizes the Boltzmann weight so that in the
    classical limit (P(E) -> Heaviside(E - V_max)) this reduces to the
    standard Arrhenius form k = A * exp(-V_max / kB T).

    Parameters
    ----------
    A_prefactor : float
        Pre-exponential (attempt frequency) factor, same units as k's output.
    E_max_factor : float
        Integrate energies up to E_max_factor * V_max to capture the full
        thermally-populated range.
    """
    V_max = eckart_barrier_height(V1, V2)
    E_max = E_max_factor * V_max
    E_values = np.linspace(1e-6 * V_max, E_max, n_E)

    P_values = np.array([wkb_transmission(E, V1, V2, L, mass) for E in E_values])
    boltzmann_weight = np.exp(-E_values / (kB * T))

    integral = np.trapezoid(P_values * boltzmann_weight, E_values)
    k = A_prefactor * integral / (kB * T)
    return k
