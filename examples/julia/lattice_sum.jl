# SPDX-FileCopyrightText: 2025 Jonathan K. Busse <jonathan@jbusse.de>
# SPDX-FileCopyrightText: 2025 David Gómez-Castro <david.gomezcastro@uam.es>
# SPDX-License-Identifier: AGPL-3.0-only

"""
Minimal working example for the Epstein zeta Library.
If the library is installed by, e.g. `using Pkg; Pkg.add("EpsteinLib")`
run with: `julia lattice_sum.jl`
"""

using EpsteinLib
using Printf

ν = 1.0
A = [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
]
x = [0.0, 0.0, 0.0]
y = [0.5, 0.5, 0.5]

# Calculate Madelung constant
madelung = epsteinzeta(ν, A, x, y)

# Reference value and relative error
madelung_ref = -1.7475645946331821906362120355443974
relerr = abs(madelung_ref - real(madelung)) / abs(madelung_ref)

println("Madelung sum in 3 dimensions:\t", real(madelung))
println("Reference value:\t\t", madelung_ref)
@printf("Relative error:\t\t\t +%.2e\n", relerr)
