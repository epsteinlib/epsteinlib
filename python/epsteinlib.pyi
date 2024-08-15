import cython
import numpy as np
from numpy.typing import NDArray as NDArray

def validate_inputs(
    nu: float | int,
    A: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> None: ...
def prepare_inputs(
    nu: float | int,
    A: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> tuple[
    float | int,
    int,
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]: ...
def epstein_zeta_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[None],
    x: cython.double[None],
    y: cython.double[None],
) -> complex: ...
def epstein_zeta(
    nu: float | int,
    A: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> complex: ...
def epstein_zeta_reg_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[None],
    x: cython.double[None],
    y: cython.double[None],
) -> complex: ...
def epstein_zeta_reg(
    nu: float | int,
    A: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> complex: ...
