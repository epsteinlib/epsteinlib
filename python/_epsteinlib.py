"""
Python wrapper for the Epstein Zeta function
"""

# SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: AGPL-3.0-only

import cython
import numpy as np
from cython.cimports.epsteinlib import epsteinZeta, epsteinZetaReg
from numpy.typing import NDArray


def validate_inputs(
    nu: float | int,
    A: NDArray[np.float64],  # pylint: disable=invalid-name
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> None:
    """
    Validate the inputs for the Epstein zeta function calculation.

    Raises:
    TypeError: If x or y is not a 1D NumPy array,
    or if nu is not a real number
    ValueError: If x and y have different lengths,
    or if A is not a 2D array with correct dimensions
    """
    if (
        not isinstance(x, np.ndarray)
        or x.ndim != 1
        or not np.issubdtype(x.dtype, np.number)
        or np.issubdtype(x.dtype, complex)
    ):
        raise TypeError("x must be a 1D NumPy array of real numbers")
    if (
        not isinstance(y, np.ndarray)
        or y.ndim != 1
        or not np.issubdtype(y.dtype, np.number)
        or np.issubdtype(y.dtype, complex)
    ):
        raise TypeError("y must be a 1D NumPy array of real numbers")

    if x.shape[0] != y.shape[0]:
        raise ValueError("x and y must have the same length")

    dim = x.shape[0]
    if dim < 1:
        raise ValueError("x and y cannot be empty arrays")

    if (
        not isinstance(A, np.ndarray)
        or A.shape != (dim, dim)
        or not np.issubdtype(A.dtype, np.number)
        or np.issubdtype(A.dtype, complex)
    ):
        raise ValueError(
            f"""A must be a 2D NumPy array with shape"
            ({dim}, {dim}) of real numbers"""
        )

    if not isinstance(nu, (float, int)):
        raise TypeError("nu must be a real number")


def prepare_inputs(
    nu: float | int,
    A: NDArray[np.float64],  # pylint: disable=invalid-name
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> tuple[
    float | int,
    int,
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """
    Prepare the inputs for the Epstein zeta function calculation.
    Makes sure, the input arrays are C_CONTIGUOUS of type float, see
    https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flags.html
    """
    a = A.reshape(-1)
    dim = x.shape[0]

    # Ensure the arrays are of type float64
    A = np.array(A, dtype=np.float64)
    x = np.array(x, dtype=np.float64)
    y = np.array(y, dtype=np.float64)

    # ENSURE arrays are C_CONTIGUOUS
    a = np.ascontiguousarray(a) if not a.flags.c_contiguous else a
    x = np.ascontiguousarray(x) if not x.flags.c_contiguous else x
    y = np.ascontiguousarray(y) if not y.flags.c_contiguous else y
    return nu, dim, a, x, y


def epstein_zeta_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[::1],
    x: cython.double[::1],
    y: cython.double[::1],
) -> complex:
    """
    Call the C function to calculate the Epstein zeta function.
    """
    return epsteinZeta(  # type: ignore [no-any-return]
        nu,
        dim,
        cython.address(a[0]),
        cython.address(x[0]),
        cython.address(y[0]),
    )


def epstein_zeta(
    nu: float | int,
    A: NDArray[np.float64],  # pylint: disable=invalid-name
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> complex:
    """
    Calculate the Epstein zeta function.
    """
    validate_inputs(nu, A, x, y)
    nu, dim, a, x, y = prepare_inputs(nu, A, x, y)
    return epstein_zeta_c_call(nu, dim, a, x, y)


def epstein_zeta_reg_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[::1],
    x: cython.double[::1],
    y: cython.double[::1],
) -> complex:
    """
    Call the C function to calculate the regularized Epstein zeta function.
    """
    return epsteinZetaReg(  # type: ignore [no-any-return]
        nu,
        dim,
        cython.address(a[0]),
        cython.address(x[0]),
        cython.address(y[0]),
    )


def epstein_zeta_reg(
    nu: float | int,
    A: NDArray[np.float64],  # pylint: disable=invalid-name
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> complex:
    """
    Calculate the regularized Epstein zeta function.
    """
    validate_inputs(nu, A, x, y)
    nu, dim, a, x, y = prepare_inputs(nu, A, x, y)
    return epstein_zeta_reg_c_call(nu, dim, a, x, y)
