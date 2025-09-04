"""
Python wrapper for the Epstein Zeta function
"""

# SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: AGPL-3.0-only

from typing import Any, Union

import cython
import numpy as np
from cython.cimports.epsteinlib import epsteinZeta, epsteinZetaReg, setZetaDer
from numpy.typing import NDArray


def validate_inputs(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
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
        or not (
            np.issubdtype(x.dtype, np.integer)
            or np.issubdtype(x.dtype, np.floating)
        )
    ):
        raise TypeError(
            "x must be a 1D NumPy array of real numbers (int or float)"
        )
    if (
        not isinstance(y, np.ndarray)
        or y.ndim != 1
        or not (
            np.issubdtype(y.dtype, np.integer)
            or np.issubdtype(y.dtype, np.floating)
        )
    ):
        raise TypeError(
            "y must be a 1D NumPy array of real numbers (int or float)"
        )

    if x.shape[0] != y.shape[0]:
        raise ValueError("x and y must have the same length")

    dim = x.shape[0]
    if dim < 1:
        raise ValueError("x and y cannot be empty arrays")

    if (
        not isinstance(A, np.ndarray)
        or A.shape != (dim, dim)
        or not (
            np.issubdtype(A.dtype, np.integer)
            or np.issubdtype(A.dtype, np.floating)
        )
    ):
        raise ValueError(
            f"A must be a 2D NumPy array with shape ({dim}, {dim}) of real numbers (int or float)"
        )

    if not isinstance(nu, (float, int)):
        raise TypeError("nu must be a real number (int or float)")


def prepare_inputs(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> tuple[
    np.float64,
    int,
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """
    Prepare the inputs for the Epstein zeta function calculation.
    Makes sure, the input arrays are C_CONTIGUOUS of type float64, see
    https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flags.html
    """

    # Convert (dim, dim) matrix A to 1D array a for cython compatibility
    dim = x.shape[0]
    a = A.reshape(-1)  # pylint: disable=invalid-name

    # Ensure arrays and nu are of type float64
    nu_out: np.float64 = np.float64(nu)
    a_out: NDArray[np.float64] = a.astype(np.float64, copy=False)
    x_out: NDArray[np.float64] = x.astype(np.float64, copy=False)
    y_out: NDArray[np.float64] = y.astype(np.float64, copy=False)

    # Ensure arrays are C_CONTIGUOUS
    a_out = np.ascontiguousarray(a_out)
    x_out = np.ascontiguousarray(x_out)
    y_out = np.ascontiguousarray(y_out)

    return nu_out, dim, a_out, x_out, y_out


def validate_inputs_der(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
    alpha: NDArray[np.integer[Any]],
) -> None:
    """
    Validate the inputs for the Epstein zeta function calculation.

    Raises:
    TypeError: If x, y or alpha is not a 1D NumPy array,
    or if nu is not a real number
    ValueError: If x and y have different lengths,
    if A is not a 2D array with correct dimensions
    or if alpha as negative entries
    """
    if (
        not isinstance(x, np.ndarray)
        or x.ndim != 1
        or not (
            np.issubdtype(x.dtype, np.integer)
            or np.issubdtype(x.dtype, np.floating)
        )
    ):
        raise TypeError(
            "x must be a 1D NumPy array of real numbers (int or float)"
        )
    if (
        not isinstance(y, np.ndarray)
        or y.ndim != 1
        or not (
            np.issubdtype(y.dtype, np.integer)
            or np.issubdtype(y.dtype, np.floating)
        )
    ):
        raise TypeError(
            "y must be a 1D NumPy array of real numbers (int or float)"
        )
    if (
        not isinstance(alpha, np.ndarray)
        or alpha.ndim != 1
        or not (np.issubdtype(alpha.dtype, np.integer))
        or (alpha < 0).any()
    ):
        raise TypeError(
            "alpha must be a 1D NumPy array of non-negative integers"
        )

    if x.shape[0] != y.shape[0] or x.shape[0] != alpha.shape[0]:
        raise ValueError("x, y and alpha must have the same length")

    dim = x.shape[0]
    if dim < 1:
        raise ValueError("x, y and alpha cannot be empty arrays")

    if (
        not isinstance(A, np.ndarray)
        or A.shape != (dim, dim)
        or not (
            np.issubdtype(A.dtype, np.integer)
            or np.issubdtype(A.dtype, np.floating)
        )
    ):
        raise ValueError(
            f"A must be a 2D NumPy array with shape ({dim}, {dim}) of real numbers (int or float)"
        )

    if not isinstance(nu, (float, int)):
        raise TypeError("nu must be a real number (int or float)")


def prepare_inputs_der(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
    alpha: NDArray[np.integer[Any]],
) -> tuple[
    np.float64,
    int,
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.unsignedinteger[Any]],
]:
    """
    Prepare the inputs for the Epstein zeta function calculation.
    Makes sure, the input arrays are C_CONTIGUOUS of type float64, see
    https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flags.html
    """

    # Convert (dim, dim) matrix A to 1D array a for cython compatibility
    dim = x.shape[0]
    a = A.reshape(-1)  # pylint: disable=invalid-name

    # Ensure arrays and nu are of type float64
    nu_out: np.float64 = np.float64(nu)
    a_out: NDArray[np.float64] = a.astype(np.float64, copy=False)
    x_out: NDArray[np.float64] = x.astype(np.float64, copy=False)
    y_out: NDArray[np.float64] = y.astype(np.float64, copy=False)
    alpha_out: NDArray[np.unsignedinteger[Any]] = alpha.astype(
        np.uint32, copy=False
    )

    # Ensure arrays are C_CONTIGUOUS
    a_out = np.ascontiguousarray(a_out)
    x_out = np.ascontiguousarray(x_out)
    y_out = np.ascontiguousarray(y_out)
    alpha_out = np.ascontiguousarray(alpha_out)

    return nu_out, dim, a_out, x_out, y_out, alpha_out


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
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> complex:
    """
    Calculate the Epstein zeta function.
    """
    validate_inputs(nu, A, x, y)
    nu_cython, dim, a_cython, x_cython, y_cython = prepare_inputs(nu, A, x, y)
    return epstein_zeta_c_call(nu_cython, dim, a_cython, x_cython, y_cython)


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
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> complex:
    """
    Calculate the regularized Epstein zeta function.
    """
    validate_inputs(nu, A, x, y)
    nu_cython, dim, a_cython, x_cython, y_cython = prepare_inputs(nu, A, x, y)
    return epstein_zeta_reg_c_call(
        nu_cython, dim, a_cython, x_cython, y_cython
    )


def set_zeta_der_c_call(  # pylint: disable=too-many-arguments, too-many-positional-arguments
    nu: cython.double,
    dim: cython.int,
    a: cython.double[::1],
    x: cython.double[::1],
    y: cython.double[::1],
    alpha: cython.uint[::1],
) -> complex:
    """
    Call the C function to calculate the regularized Epstein zeta function.
    """
    return setZetaDer(  # type: ignore [no-any-return]
        nu,
        dim,
        cython.address(a[0]),
        cython.address(x[0]),
        cython.address(y[0]),
        cython.address(alpha[0]),
    )


def set_zeta_der(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
    alpha: NDArray[np.integer[Any]],
) -> complex:
    """
    Calculate the derivatives of the set zeta function.
    """
    validate_inputs_der(nu, A, x, y, alpha)
    nu_cython, dim, a_cython, x_cython, y_cython, alpha_cython = (
        prepare_inputs_der(nu, A, x, y, alpha)
    )
    return set_zeta_der_c_call(
        nu_cython, dim, a_cython, x_cython, y_cython, alpha_cython
    )
