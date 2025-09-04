"""
Minimal working example for the derivatives of the Epstein zeta function.
If the library is installed, run with `python epstein_zeta_derivatives.py`.
"""

# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
# SPDX-License-Identifier: AGPL-3.0-only


import numpy as np

from epsteinlib import epstein_zeta, set_zeta_der


def check_epstein_zeta_derivative() -> bool:
    """
    Calculate derivative of the epstein zeta function by finite differences
    and compare to the set zeta derivatives as provided by the library.

    Note in particular, that set_zeta_der(nu, a, x, y, alpha) is the partial
    derivative with respect to y and a multi-index alpha of
    exp(2 pi i  x y) epstein_zeta(nu, a, x, y)

    @return true if the relative_error is smaller than 10**(-10)
    """

    nu = 1.0
    dim = 3
    a = np.identity(dim)
    x = np.zeros(dim)
    y = np.full(dim, 0.5)
    alpha = np.array([2, 0, 0], dtype=np.uint32)

    h = 0.00005
    y_minus_h = y + np.array([h, 0, 0])
    y_plus_h = y - np.array([h, 0, 0])

    epstein_zeta_der_fin_diff = np.real(
        epstein_zeta(nu, a, x, y_plus_h)
        - 2 * epstein_zeta(nu, a, x, y)
        + epstein_zeta(nu, a, x, y_minus_h)
    ) / (h * h)

    epstein_zeta_der_ref = np.real(set_zeta_der(nu, a, x, y, alpha))

    relative_error = abs(
        epstein_zeta_der_fin_diff - epstein_zeta_der_ref
    ) / abs(epstein_zeta_der_ref)

    print(
        f"Epstein zeta derivative finite differences:\t {epstein_zeta_der_fin_diff:.16f}"
    )
    print(
        f"Epstein zeta derivative internal routine:  \t {epstein_zeta_der_ref:.16f}"
    )
    print(f"Relative error:                          \t+{relative_error:.2e}")

    return relative_error > 10**-10


if __name__ == "__main__":
    check_epstein_zeta_derivative()
