# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
#
# SPDX-License-Identifier: AGPL-3.0-only


"""Tests for python wrapper of the Epstein Zeta function"""

import unittest
from typing import Callable

import benchmark_functions as bf
import numpy as np
from numpy.typing import NDArray

from epsteinlib import (
    epstein_zeta,
    epstein_zeta_reg,
    prepare_inputs,
    validate_inputs,
)


class TestEpsteinZeta(unittest.TestCase):
    """
    Unit test class for testing the Epstein Zeta function
    by comparison to analytic representations in special cases.
    """

    # Class variables
    offset: float = (
        0.01  # problems for python values with 0ffset <= 0.001, no problems with matheamtica values
    )
    stepsize: float = 0.1
    nu_values: NDArray[np.float64] = np.concatenate(
        [
            np.arange(-9 + offset, -8 + offset, stepsize),
            np.arange(-1 + offset, 4 + offset, stepsize),
        ]
    )
    threshold: float = 2 * 10 ** (-13)

    def compare_epstein_zeta_with_ref(
        self,
        a: NDArray[np.float64],
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        epstein_zeta_ref_func: Callable[[float], float],
    ) -> None:
        """
        Test the Epstein Zeta function and the regularized Epstein Zeta
        function for a range of nu values by comparing with a reference
        function.

        Raises:
            AssertionError: If the error between the computed and Reference
            values exceeds the threshold.
        """
        for nu in self.nu_values:
            with self.subTest(nu=nu):
                # test epstein_zeta
                epstein_zeta_num = epstein_zeta(nu, a, x, y)
                epstein_zeta_ref = epstein_zeta_ref_func(nu)
                error = bf.min_errors_abs_error_rel(
                    epstein_zeta_ref, epstein_zeta_num
                )
                self.assertLessEqual(
                    error,
                    self.threshold,
                    f"""
                    epstein_zeta({nu}, {np.ndarray.flatten(a)}, {x}, {y})
                    """,
                )
                # test epstein_zeta_reg
                epstein_zeta_reg_num = epstein_zeta_reg(nu, a, x, y)
                epstein_zeta_reg_ref = np.exp(
                    2 * np.pi * 1j * np.dot(x, y)
                ) * epstein_zeta_ref - bf.singularity_in_id(
                    y, nu, np.size(x)
                ) / np.abs(
                    np.linalg.det(a)
                )
                error = bf.min_errors_abs_error_rel(
                    epstein_zeta_reg_ref, epstein_zeta_reg_num
                )
                self.assertLessEqual(
                    error,
                    self.threshold,
                    f"""
                    epstein_zeta_reg({nu}, {np.ndarray.flatten(a)}, {x}, {y})
                    """,
                )

    def test_epstein_zeta_00_mhalfmhalf_id(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[0,0], y=[-1/2,-1/2], and a is the identity matrix.
        """
        a: NDArray[np.float64] = np.identity(2)
        x: NDArray[np.float64] = np.array([0, 0])
        y: NDArray[np.float64] = np.array([-1 / 2, -1 / 2])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_00_mhalfmhalf_id,
        )

    def test_epstein_zeta_m1m1_halfhalf_id(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[-1,-1], y=[1/2,1/2], and a is the identity matrix.
        """
        a: NDArray[np.float64] = np.identity(2)
        x: NDArray[np.float64] = np.array([-1, -1])
        y: NDArray[np.float64] = np.array([1 / 2, 1 / 2])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_m1m1_halfhalf_id,
        )

    def test_epstein_zeta_m1m1_half0_id(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[-1,-1], y=[1/2,0], and a is the identity matrix.
        """
        a: NDArray[np.float64] = np.identity(2)
        x: NDArray[np.float64] = np.array([-1, -1])
        y: NDArray[np.float64] = np.array([1 / 2, 0])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_m1m1_half0_id,
        )

    def epstein_zeta_onehalf0sqrt3half_00_00(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[1/2,0,0,0], y=[0,0,0,0], and a= [[1, 1/2], [0, sqrt(3)/2]].
        """
        a: NDArray[np.float64] = np.array([[1, 1 / 2], [0, np.sqrt(3) / 2]])
        x: NDArray[np.float64] = np.array([0, 0])
        y: NDArray[np.float64] = np.array([0, 0])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_half000_0000_id,
        )

    def test_epstein_zeta_diag2sqrt242_0m1m1_4sqrt2th00(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[0,-1,-1], y=[1/4*sqrt(2),0,0], and a=diag(2*sqrt(2),4,2).
        """
        a: NDArray[np.float64] = np.diag([2 * np.sqrt(2), 4, 2])
        x: NDArray[np.float64] = np.array([0, -1, -1])
        y: NDArray[np.float64] = np.array([1 / (4 * np.sqrt(2)), 0, 0])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_diag2sqrt242_0m1m1_4sqrt2th00,
        )

    def test_epstein_zeta_half000_0000_id(self) -> None:
        """
        Test the Epstein Zeta function for the case where
        x=[1/2,0,0,0], y=[0,0,0,0], and a is the identity matrix.
        """
        a: NDArray[np.float64] = np.identity(4)
        x: NDArray[np.float64] = np.array([1 / 2, 0, 0, 0])
        y: NDArray[np.float64] = np.array([0, 0, 0, 0])
        self.compare_epstein_zeta_with_ref(
            a,
            x,
            y,
            bf.epstein_zeta_half000_0000_id,
        )

    def test_singularity(self) -> None:
        """
        Test if the Epstein Zeta function returns nan for singularity case.
        """
        nu: float = 1.0
        a: NDArray[np.float64] = np.array([[1.0]])
        x: NDArray[np.float64] = np.array([0.0])
        y: NDArray[np.float64] = np.array([0.0])

        result = epstein_zeta(nu, a, x, y)
        self.assertTrue(np.isnan(result), f"Expected nan, but got {result}")

    def test_no_singularity_in_reg(self) -> None:
        """
        Test if the Epstein Zeta function returns nan for singularity case.
        """
        nu: float = 1.0
        a: NDArray[np.float64] = np.array([[1.0]])
        x: NDArray[np.float64] = np.array([0.0])
        y: NDArray[np.float64] = np.array([0.0])

        result = epstein_zeta_reg(nu, a, x, y)
        self.assertFalse(
            np.isnan(result), f"Expected a number, but got {result}"
        )
        self.assertTrue(
            np.isfinite(result), f"Expected a finite number, but got {result}"
        )


class TestValidateInputs(unittest.TestCase):
    """
    Unit test class for validate_inputs function in epsteinzetalib.
    """

    def test_valid_inputs(self) -> None:
        """
        Ensures that no exceptions are raised for valid inputs.
        """
        np.random.seed(0)  # Set the seed for reproducibility
        for _ in range(10):  # Generate and test 10 random matrices
            length = np.random.randint(1, 10)
            nu = np.random.rand() * 10
            a = np.random.rand(length, length)
            x = np.random.rand(length)
            y = np.random.rand(length)
            try:
                validate_inputs(nu, a, x, y)
            except Exception as e:  # pylint: disable=broad-exception-caught
                self.fail(
                    f"validate_inputs raised {type(e).__name__} unexpectedly!"
                )

    def test_invalid_x_type(self) -> None:
        """
        Ensures that a TypeError is raised for invalid x types.
        """
        np.random.seed(1)
        invalid_x_types = [
            [1, 2, 3],
            np.array([[1], [2], [3]]),
            np.array(["1", 2, 2]),
            "1 2 3",
            {"1": 1, "2": 2, "3": 3},
            np.array([1j, 2, 3]),
        ]  # not non-complex numpy arrays of subtype np.number
        for x in invalid_x_types:
            nu = np.random.rand() * 10
            a = np.random.rand(3, 3)
            y = np.random.rand(3)
            with self.assertRaises(TypeError):
                validate_inputs(nu, a, x, y)  # type: ignore[arg-type]

    def test_invalid_y_type(self) -> None:
        """
        Ensures that a TypeError is raised for invalid y types.
        """
        np.random.seed(2)
        invalid_y_types = [
            [1, 2, 3],
            np.array([[1], [2], [3]]),
            np.array(["1", 2, 2]),
            "1 2 3",
            {"1": 1, "2": 2, "3": 3},
            np.array([1j, 2, 3]),
        ]  # not non-complex numpy arrays of subtype np.number
        for y in invalid_y_types:
            nu = np.random.rand() * 10
            a = np.random.rand(3, 3)
            x = np.random.rand(3)
            with self.assertRaises(TypeError):
                validate_inputs(nu, a, x, y)  # type: ignore[arg-type]

    def test_x_y_different_lengths(self) -> None:
        """
        Ensures that a ValueError is raised when x and y have different lengths.
        """
        np.random.seed(3)
        odd_dimension = 11
        for length_x in range(1, odd_dimension):
            length_y = odd_dimension - length_x  # different length than x
            nu = np.random.rand() * 10
            a = np.random.rand(length_x, length_x)
            x = np.random.rand(length_x)
            y = np.random.rand(length_y)
            with self.assertRaises(ValueError):
                validate_inputs(nu, a, x, y)

    def test_empty_inputs(self) -> None:
        """
        Ensures that a ValueError is raised for empty inputs.
        """
        nu = 0.0
        a = np.array([[]])
        x = np.array([])
        y = np.array([])
        with self.assertRaises(ValueError):
            validate_inputs(nu, a, x, y)

    def test_invalid_a_type(self) -> None:
        """
        Ensures that a ValueError is raised for each invalid matrices a.
        """
        np.random.seed(4)
        invalid_a_types = [
            [[1, 0], [0, 1]],
            1,
            np.array([["1", 0], [0, 1]]),
            "1 2 3",
            {"1": 1, "2": 2, "3": 3},
            np.array([[1, 0], [0, 1j]]),
        ]  # not non-complex 2D numpy arrays of subtype np.number
        for a in invalid_a_types:
            nu = np.random.rand() * 10
            x = np.random.rand(2)
            y = np.random.rand(2)
            with self.assertRaises(ValueError):
                validate_inputs(nu, a, x, y)  # type: ignore[arg-type]

    def test_invalid_a_shape(self) -> None:
        """
        Ensures that a ValueError is raised for invalid shapes of a.
        """
        np.random.seed(5)
        odd_dimension = 11
        for length in range(1, odd_dimension):
            diff_length = odd_dimension - length
            nu = np.random.rand() * 10
            a = np.random.rand(length, diff_length)  # Incorrect shape
            x = np.random.rand(length)
            y = np.random.rand(length)
            with self.assertRaises(ValueError):
                validate_inputs(nu, a, x, y)

    def test_invalid_nu_type(self) -> None:
        """
        Ensures that a TypeError is raised for invalid nu.
        """
        np.random.seed(6)
        invalid_nu_types = [
            "0.5",
            0.5 + 1j,
            [0.5],
            np.array([0.5]),
        ]  # not floats or ints
        for nu in invalid_nu_types:
            length = np.random.randint(1, 10)
            a = np.random.rand(length, length)
            x = np.random.rand(length)
            y = np.random.rand(length)
            with self.assertRaises(TypeError):
                validate_inputs(nu, a, x, y)  # type: ignore[arg-type]


class TestPrepareInputs(unittest.TestCase):
    """
    Unit test class for prepare_inputs function in epsteinzetalib.
    """

    def test_preserve_inputs(self) -> None:
        """
        Ensures that the inputs are preserved and correctly prepared.
        """
        np.random.seed(0)  # Set the seed for reproducibility

        for _ in range(10):  # Generate and test 10 random sets of inputs
            length = np.random.randint(1, 10)
            nu = np.random.rand() * 10
            a = np.random.rand(length, length)
            x = np.random.rand(length)
            y = np.random.rand(length)

            # Prepare inputs and check returned values
            nu_prep, dim_prep, a_prep, x_prep, y_prep = prepare_inputs(
                nu, a, x, y
            )

            # Check if the dimensions are correct
            self.assertEqual(dim_prep, length)

            # Check if the values are correctly assigned
            self.assertEqual(nu_prep, nu)
            np.testing.assert_array_equal(a_prep, a.flatten())
            np.testing.assert_array_equal(x_prep, x)
            np.testing.assert_array_equal(y_prep, y)

    def test_non_contiguous_inputs(self) -> None:
        """
        Ensures that non-contiguous arrays are returned contiguous.
        """
        np.random.seed(1)

        for _ in range(10):
            length = np.random.randint(1, 10)
            nu = np.random.rand() * 10
            a = np.random.rand(length, length)
            x = np.random.rand(length)
            y = np.random.rand(length)

            # Make arrays non-contiguous
            a = np.asfortranarray(a)
            x = np.asfortranarray(x)
            y = np.asfortranarray(y)

            # Prepare inputs and check returned values
            _, _, a_prep, x_prep, y_prep = prepare_inputs(nu, a, x, y)

            # Check if the returned arrays are contiguous
            self.assertTrue(x_prep.flags.c_contiguous)
            self.assertTrue(y_prep.flags.c_contiguous)
            self.assertTrue(a_prep.flags.c_contiguous)


if __name__ == "__main__":
    unittest.main()
