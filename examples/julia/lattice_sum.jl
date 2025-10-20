# SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
# SPDX-License-Identifier: AGPL-3.0-only

##
# @file madelung_test.jl
# @brief Calculates the Madelung constant.
#
# Minimal working example for the Epstein Zeta Library (Julia interface).
# Requires `Epsteinlib_jll`.
#
# Usage:
# ```bash
# julia madelung_test.jl
# ```
#
# Reference Madelung constant from literature:
# sum_{i,j,k∈ℤ} (-1)^(i+j+k) / sqrt(i² + j² + k²)
#
# @date 20/10/2025
# @authors
#   Andreas Buchheit
#   Jonathan Busse
#   Ruben Gutendorf
##

using Pkg; Pkg.add("Epsteinlib_jll")
using Epsteinlib_jll
using Printf

##
# @brief Calculate Madelung constant and compare to reference.
#
# Computes:
#   sum_{i,j,k∈ℤ} (-1)^(i+j+k) / sqrt(i² + j² + k²)
# using the Epstein zeta function from the `Epsteinlib_jll` package.
#
# @return 0 if relative error < 1e-14, otherwise 1.
##
function madelung_test()
    madelung_ref = -1.7475645946331821906362120355443974

    nu  = 1.0
    dim = UInt32(3)
    A = [1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0]
    x = [0.0, 0.0, 0.0]
    y = [0.5, 0.5, 0.5]

    madelung = ccall((:epsteinZeta, Epsteinlib_jll.libepstein),
                     ComplexF64,
                     (Cdouble, Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                     nu, dim, A, x, y)

    println("Madelung sum in 3 dimensions:\t", real(madelung))
    println("Reference value:\t\t", madelung_ref)

    relerr = abs(madelung_ref - real(madelung)) / abs(madelung_ref)
    @printf("Relative error:\t\t\t +%.2e\n", relerr)

    return relerr > 1e-14 ? 1 : 0
end

exit(madelung_test())
