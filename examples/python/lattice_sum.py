"""
Minimal working example for the Epstein Zeta Library.
If the library is installed, run with `python lattice_sum.py`.
"""

# SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
# SPDX-License-Identifier: AGPL-3.0-only


import numpy as np

from epsteinlib import epstein_zeta


def check_madelung() -> bool:
    """
    Calculate Madelung constant and compare to precomputed value.

    Madelung constant:
    sum_{i, j, k in Z} (-1)**(i + j + k) / sqrt(i**2 + j**2 + k**2)

    @return true if the difference to precomputed value is smaller than 10**(-14)
    """
    # Madelung constant found in literature
    madelung_ref = -1.7475645946331821906362120355443974
    dim = 3
    a = np.identity(dim)  # identity matrix for whole numbers
    x = np.zeros(dim)  # no shift
    y = np.full(dim, 0.5)  # alternating sum
    nu = 1.0
    madelung = np.real(epstein_zeta(nu, a, x, y))
    print(f"Madelung sum in 3 dimensions:\t {madelung:.16f}")
    print(f"Reference value:\t\t {madelung_ref:.16f}")
    print(
        f"Relative error:\t\t\t +{abs(madelung_ref - madelung) / abs(madelung_ref):.2e}"
    )

    return abs(madelung - madelung_ref) > 10**-14


if __name__ == "__main__":
    check_madelung()
