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
def validate_inputs_der(
    nu: Union[float, int],
    A: NDArray[  # pylint: disable=invalid-name
        Union[np.integer[Any], np.floating[Any]]
    ],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
    alpha: NDArray[np.integer[Any]],
) -> None: ...
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
    NDArray[np.unsignedinteger],
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
def set_zeta_der_c_call(
    nu: cython.double,
    dim: cython.int,
    a: cython.double[::1],
    x: cython.double[::1],
    y: cython.double[::1],
    alpha: cython.uint[::1],
) -> complex: ...
def set_zeta_der(
    nu: Union[float, int],
    A: NDArray[Union[np.integer[Any], np.floating[Any]]],
    x: NDArray[Union[np.integer[Any], np.floating[Any]]],
    y: NDArray[Union[np.integer[Any], np.floating[Any]]],
    alpha: NDArray[np.integer],
) -> complex: ...
