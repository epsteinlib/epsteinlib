# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
#
# Original source:
# https://github.com/epsteinlib/epsteinlib/blob/main/examples/c/lattice_sum.c
#
# Original file copyright:
# // SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# // SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
# // SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
# // SPDX-License-Identifier: AGPL-3.0-only
#
# Description:
# Calculates the 3D Madelung constant using the Epstein zeta function


using Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/JuliaBinaryWrappers/Epsteinlib_jll.jl.git")

using Epsteinlib_jll
using Printf

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

exit(relerr > 1e-14 ? 1 : 0)
