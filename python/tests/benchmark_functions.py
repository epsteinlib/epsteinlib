# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
#
# SPDX-License-Identifier: AGPL-3.0-only


"""Reference functions for benchmarking purposes."""

import numpy as np
from mpmath import gamma, mp, pi, zeta
from numpy.typing import NDArray

# 50 decimal places precision in mpmath.zeta
# this is needed for precise evaluation of the
# hurwitz zeta function for benchmarking
mp.dps = 50


def dirichlet_eta(s: float) -> float:
    """
    Compute the Dirichlet eta function for a given s.
    Representation in terms of the Riemann Zeta function.
    """
    return (1 - 2 ** (1 - s)) * float(zeta(s))


def dirichlet_beta(s: float) -> float:
    """
    Compute the Dirichlet eta function for a given s.
    Representation in Terms of the Hurwitz zeta function.
    """
    return 4 ** (-s) * (float(zeta(s, 1 / 4)) - float(zeta(s, 3 / 4)))


def dirichlet_lambda(s: float) -> float:
    """
    Compute the Dirichlet lambda function for a given s.
    Representation in Terms of the Riemann zeta function.
    """
    return (1 - 2 ** (-s)) * float(zeta(s))


def singularity_in_id(y: NDArray[np.float64], s: float, dim: int) -> float:
    """
    Compute the Singularity of the Epstein Zeta functions
    when y goes to zero for any vector x, identity matrix a
    and nu ≠ 2 + n for any natural number n including zero.

    epstein_zeta_reg(x, y, id, nu) =
        exp(2 * Pi * i * x . y)
        * epstein_zeta(x, y, id, nu)
        - singularity(y, nu)
    """
    y2 = np.sqrt(np.dot(y, y))
    if y2 == 0:
        return 0
    return float(
        y2 ** (s - dim)
        * pi ** (s - dim / 2)
        * gamma((dim - s) / 2)
        / gamma(s / 2)
    )


def epstein_zeta_00_mhalfmhalf_id(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu ≠ 2 and

    x = [0, 0]
    y = [- 1 / 2, - 1 / 2]
    a = [[1, 0], [0, 1]]

    Representation in terms of the Dirichlet Eta function
    and the Dirichlet Beta function.
    """
    return -4 * dirichlet_eta(nu / 2) * dirichlet_beta(nu / 2)


def epstein_zeta_m1m1_halfhalf_id(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu ≠ 2 and

    x = [-1, -1]
    y = [1 / 2, 1 / 2]
    a = [[1, 0], [0, 1]]

    Representation in terms of the Dirichlet Beta function
    and the Zeta function.
    """
    return (
        -(2 ** (2 - nu / 2))
        * (2 ** (nu / 2) - 2)
        * dirichlet_beta(nu / 2)
        * float(zeta(nu / 2))
    )


def epstein_zeta_m1m1_half0_id(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu ≠ 2 and

    x = [-1, -1]
    y = [1 / 2, 0]
    a = [[1, 0], [0, 1]]

    Representation in terms of the Dirichlet Beta function
    and the Zeta function.
    """
    return (
        2 ** (2 - nu)
        * (2 ** (nu / 2) - 2)
        * dirichlet_beta(nu / 2)
        * float(zeta(nu / 2))
    )


def epstein_zeta_onehalf0sqrt3half_00_00(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu and

    x = [0, 0]
    y = [0, 0]
    a = [[1, 1/2], [0, sqrt(3)/2]]

    Representation in terms of the (hurwitz) Zeta function.
    """
    return (
        3 ** (1 - nu / 2)
        * 2
        * float(zeta(nu / 2) * (zeta(nu / 2, 1 / 3) - zeta(nu / 2, 2 / 3)))
    )


def epstein_zeta_diag2sqrt242_0m1m1_4sqrt2th00(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu and

    x = [0, -1, -1]
    y = [-4*sqrt(2), 0, 0]
    a = [[2sqrt(2), 0, 0], [0, 4, 0], [0, 0, 2]]

    Representation in terms of the Dirichlet Beta function.
    """
    return 2 ** (1 - nu / 2) * dirichlet_beta(nu - 1)


def epstein_zeta_half000_0000_id(nu: float) -> float:
    """
    Compute the Epstein Zeta function for a given nu ≠ 2 and

    x = [1 / 2, 0, 0, 0]
    y = [0, 0, 0, 0]
    a = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

    Representation in terms of the Dirichlet Beta function
    and the Zeta function.
    """
    return 2 ** (nu) * (
        dirichlet_lambda(nu / 2) * dirichlet_lambda(nu / 2 - 1)
        + dirichlet_beta(nu / 2) * dirichlet_beta(nu / 2 - 1)
    )


def min_errors_abs_error_rel(true: complex, approx: complex) -> float:
    """Return the minimum of absolute and relative error."""
    error_abs = (
        (true.real - approx.real) ** 2 + (true.imag - approx.imag) ** 2
    ) ** (1 / 2)
    error_rel = (
        error_abs * ((true.real) ** 2 + (true.imag) ** 2) ** (-1 / 2)
        if true != 0
        else float("inf")
    )
    return float(min(error_abs, error_rel))
