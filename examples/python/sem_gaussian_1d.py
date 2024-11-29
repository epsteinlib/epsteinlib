"""
Singular Euler-Maclaurin (SEM) expansion for a Gaussian function in 1D.

This script demonstrates the application of the Singular Euler-Maclaurin
expansion to a Gaussian function in one dimension. It calculates and compares
the integral, sum, and SEM approximations for various orders.

Features:
- Calculates and plots integral, sum, and SEM approximations
- Compares SEM approximations of different orders
- Provides error analysis between sum and SEM approximations

Usage:
    python sem_gaussian_1d.py [--nu NU]

Arguments:
    --nu NU    Optional. Set the value of nu for calculations. Default is 1.5.
               Example: python sem_gaussian_1d.py --nu 1

The script generates plots comparing the different approximations and their
errors for the specified nu value.
"""

# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
# SPDX-License-Identifier: AGPL-3.0-only


import argparse
from typing import Callable, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from mpmath import factorial, gamma, hyp1f1
from numpy.typing import NDArray

from epsteinlib import epstein_zeta_reg

EPS = 1e-8  # taylor expansion at nu = 1 - 2 EPS in integral and lattice_contribution


def gaussian(
    y: Union[float, NDArray[np.float64]], sigma: float
) -> Union[float, NDArray[np.float64]]:
    """Gaussian function with input y and standard deviation sigma."""
    return np.exp(-np.pi * (y / sigma) ** 2)


def power_law(
    k: Union[float, NDArray[np.float64]], nu: float
) -> Union[float, NDArray[np.float64]]:
    """Power law function with input k and exponent nu."""
    return 1 / (np.sqrt(k**2) ** nu)


def finite_differences_coefficients(
    num_points: int, derivative_order: int = 1
) -> NDArray[np.float64]:
    """Calculate finite differences coefficients for Taylor methods."""

    def matrix_element(row: int, col: int) -> np.float64:
        if row == 1:
            return np.float64(1)
        center = np.float64(num_points + 1) / 2
        return np.float64(1 / factorial(row - 1)) * (col - center) ** (row - 1)

    # Create the matrix
    matrix = np.array(
        [
            [
                matrix_element(i, j)
                for j in np.arange(1, num_points + 1, dtype=np.float64)
            ]
            for i in np.arange(1, num_points + 1, dtype=np.float64)
        ]
    )

    # Create the right-hand side vector
    rhs = np.zeros(num_points, dtype=np.float64)
    rhs[derivative_order] = np.float64(num_points + 1) * derivative_order

    # Solve the system of equations
    return np.linalg.solve(matrix, rhs)


def finite_differences(
    func: Callable[[float], float],
    x_val: float,
    order: int = 1,
    grid_points: int = 3,
    step_size: float = 1e-5,
) -> float:
    """Calculate finite differences for given function and parameters."""
    if order == 0:
        return func(x_val)

    # Calculate coefficients and scaling factor
    coeffs = finite_differences_coefficients(grid_points, order)
    scaling_factor = (step_size**order) * order * (grid_points + 1)

    # Calculate the center point of the grid
    center = (grid_points + 1) / 2 - 1

    # Sum up the weighted function values
    result = 0.0
    for i in range(grid_points):
        x_i = x_val + step_size * (i - center)
        f_x_i = func(x_i)
        weighted_value = coeffs[i] * f_x_i
        result += weighted_value / scaling_factor

    return result


def binomial(a: int, b: int) -> float:
    """Binomial coefficients"""
    return float(factorial(a) / (factorial(b) * factorial(a - b)))


def hermite_number(n: int) -> int:
    """Value of the n'th Hermite polynomial at 0"""
    if n % 2 == 1:
        return 0

    def dfact(n: int) -> int:
        """Double factorial"""
        result: int = 1 if n <= 0 else n * dfact(n - 2)
        return result

    result: int = (-2) ** (n // 2) * dfact(n - 1)
    return result


def gaussian_derivative(y: float, sigma: float, order: int) -> float:
    """Calculate the Gaussian derivative for given input, sigma, and order."""
    return float(
        (-1) ** order
        * gaussian(y, sigma)
        * sum(
            2**j
            * binomial(order, j)
            * hermite_number(order - j)
            * (np.pi / sigma**2) ** ((order + j) / 2)
            * (1 if (j == 0) else y ** j)
            for j in range(order + 1)
        )
    )


def lattice_contribution(
    x_val: float, nu: float, sigma: float, order: int
) -> float:
    """Calculate the lattice contribution for given parameters."""
    if nu == 1.0:
        nu0 = 1 - 2 * EPS
        return lattice_contribution(
            x_val, nu0, sigma, order
        ) + 2 * EPS * finite_differences(
            lambda nu: lattice_contribution(x_val, nu, sigma, order), nu0
        )

    def epstein_zeta_reg_wrapper(y: float) -> float:
        return float(
            np.real(
                epstein_zeta_reg(
                    nu,
                    np.array([[1.0]]),
                    np.array([0.0]),
                    np.array([np.double(y)]),
                )
            )
        )

    result = float(
        np.real(
            sum(
                (1 / factorial(j))
                * (1j / (2 * np.pi)) ** j
                * finite_differences(
                    epstein_zeta_reg_wrapper,
                    0,
                    j,
                    j + 1,
                )
                * gaussian_derivative(x_val, sigma, j)
                for j in range(0, order + 1, 2)
            )
        )
    )
    return float(result)


def sum_func(x_val: float, nu: float, sigma: float) -> float:
    """Calculate the sum function for given x, nu, and sigma."""
    sum_lim = 10 * sigma
    return float(
        np.sum(
            np.array(
                [
                    (
                        0
                        if np.isclose(y - x_val, 0)
                        else gaussian(y, sigma) * power_law(y - x_val, nu)
                    )
                    for y in np.arange(-sum_lim, sum_lim)
                ]
            )
        )
    )


def integral(x_val: float, nu: float, sigma: float) -> float:
    """Analytic representation of the integral over the summands."""
    if nu == 1.0:
        nu0 = 1.0 - 2 * EPS
        return integral(x_val, nu0, sigma) + 2 * EPS * finite_differences(
            lambda nu: integral(x_val, nu, sigma), nu0
        )
    return float(
        2**nu
        * np.pi ** (-1 + nu / 2)
        * np.sqrt(1 / (sigma**2))
        * (sigma**2) ** (1 - nu / 2)
        * gamma(1 - nu)
        * gamma(nu / 2)
        * hyp1f1(nu / 2, 1 / 2, -np.pi * x_val**2 / (sigma**2))
        * np.sin(np.pi * nu / 2)
    )


def sem(x_val: float, nu: float, sigma: float, order: int) -> float:
    """Calculate the SEM function for given parameters."""
    return integral(x_val, nu, sigma) + lattice_contribution(
        x_val, nu, sigma, order
    )


def min_abs_rel_error(
    ref: Union[float, NDArray[np.float64]],
    approx: Union[float, NDArray[np.float64]],
) -> Union[float, NDArray[np.float64]]:
    """
    Calculate the minimum between the absolute and relative error of two numbers.
    """
    abs_error = np.abs(approx - ref)
    rel_error = abs_error / (1 if np.isclose(ref, 0) else np.abs(ref))
    return np.minimum(abs_error, rel_error)


def plot_left(
    ax: Axes,
    x_values: NDArray[np.float64],
    plot_data_values: dict[str, NDArray[np.float64]],
    nu: float,
) -> None:
    """Plot the left subplot with integral, sum, and SEM (order 0) values."""
    ax.plot(x_values, plot_data_values["integral"], label="Integral")
    ax.plot(x_values, plot_data_values["sum_func"], "ro", label="Sum")
    ax.plot(x_values, plot_data_values["sem_order_0"], label="SEM (order 0)")

    ax.set_xlabel("x")
    ax.set_ylabel("Value")
    if nu == 1.0:
        ax.set_title("Integral, Sum, and SEM Plot in the nu -> 1 limit")
    else:
        ax.set_title(f"Integral, Sum, and SEM Plots in nu = {nu}")
    ax.legend()
    ax.grid()


def plot_right(
    ax: Axes,
    x_values: NDArray[np.float64],
    diff1: NDArray[np.float64],
    diff2: Union[NDArray[np.float64], None],
) -> None:
    """Plot the right subplot with error analysis."""
    ax.semilogy(
        x_values,
        diff1,
        "r-",
        linewidth=2,
        marker="o",
        markersize=4,
        markerfacecolor="none",
        markeredgecolor="r",
        label="|Sum - SEM (order 0)|",
    )
    if diff2 is not None:
        ax.semilogy(
            x_values,
            diff2,
            "b-",
            linewidth=2,
            marker="s",
            markersize=4,
            markerfacecolor="none",
            markeredgecolor="b",
            label="|Sum - SEM (order 2)|",
        )

    ax.set_xlabel("x")
    ax.set_ylabel("Minimum of absolute and relative error")
    ax.set_title("Error Analysis: |Sum - SEM|")
    ax.legend(loc="best", frameon=False)
    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.tick_params(which="minor", length=4, color="k")
    ax.tick_params(which="major", length=10, color="k")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Calculate SEM function with custom nu value."
    )
    parser.add_argument(
        "--nu",
        type=float,
        default=1.5,
        help="Exponent of power law nu (default: 1.5)",
    )
    args = parser.parse_args()

    # Set fixed values
    NU0 = args.nu
    SIGMA0 = 100

    # SEM displacement values
    x = np.linspace(-2 * SIGMA0, 2 * SIGMA0, 51)

    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6))

    # Precompute values
    integral_values = np.array([integral(xx, NU0, SIGMA0) for xx in x])
    sum_func_values = np.array([sum_func(xx, NU0, SIGMA0) for xx in x])
    sem_order_0_values = np.array([sem(xx, NU0, SIGMA0, 0) for xx in x])

    # Calculate differences for order 0
    diff_order_0 = np.array(
        [
            min_abs_rel_error(s, s0)
            for s, s0 in zip(sum_func_values, sem_order_0_values)
        ]
    )

    # Calculate and plot sem_order_2 if NU0 != 1
    if NU0 == 1:
        diff_order_2: Union[NDArray[np.float64], None] = None
    else:
        sem_order_2_values = np.array([sem(xx, NU0, SIGMA0, 2) for xx in x])

        diff_order_2 = np.array(
            [
                min_abs_rel_error(s, s2)
                for s, s2 in zip(sum_func_values, sem_order_2_values)
            ]
        )
    # Plot left and right subplots
    plot_data = {
        "integral": integral_values,
        "sum_func": sum_func_values,
        "sem_order_0": sem_order_0_values,
    }
    plot_left(ax1, x, plot_data, NU0)
    plot_right(ax2, x, diff_order_0, diff_order_2)

    plt.tight_layout()
    plt.show()
