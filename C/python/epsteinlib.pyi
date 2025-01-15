import cython
import numpy as np
from numpy.typing import NDArray as NDArray
from typing import Any, Union

def validate_inputs(
    nu: Union[float, int],
    A: NDArray[Union[np.integer[Any], np.floating[Any]]],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> None: ...
def prepare_inputs(
    nu: Union[float, int],
    A: NDArray[Union[np.integer[Any], np.floating[Any]]],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> tuple[
    np.float64,
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
    nu: Union[float, int],
    A: NDArray[Union[np.integer[Any], np.floating[Any]]],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> complex: ...
def epstein_zeta_reg_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[None],
    x: cython.double[None],
    y: cython.double[None],
) -> complex: ...
def epstein_zeta_reg(
    nu: Union[float, int],
    A: NDArray[Union[np.integer[Any], np.floating[Any]]],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
) -> complex: ...
