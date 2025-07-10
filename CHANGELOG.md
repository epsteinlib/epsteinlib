<!--
SPDX-FileCopyrightText: 2025 Andreas Buchheit <buchheit@num.uni-sb.de>
SPDX-FileCopyrightText: 2025 Jan Schmitz <schmitz@num.uni-sb.de>
SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan.busse@dlr.de>
SPDX-FileCopyrightText: 2025 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>

SPDX-License-Identifier: AGPL-3.0-only
-->

# Changelog

## [0.5.0] - 2025-07-10

### Added
- Tests that ensure that EpsteinZeta can be expressed as EpsteinZetaReg + singularity.
- Standalone example notebook `examples/mathematica/BenchmarkQuick.wls` that compares the (regularized) Epstein zeta function to known formulas in special cases.
- Extensive example notebook `examples/mathematica/BenchmarkAndPaperFigures.wls` that reproduces every figure from "Computation and Properties of the Epstein zeta function".
- Executables `./build/src/tests/epsteinlib_benchmark_gamma` and `./build/src/tests/epsteinlib_benchmark_epstein` that generate benchmark values for the above notebook.
- Improved unit test terminal output formatting and expressive test result reporting.
- Python wrapper now uses `Union` from `typing` to support Python versions $<3.10$
- Python wrapper supports any `np.floating[Any]` or `np.integer[Any]` type arrays for `A`, `x`, and `y`
- Loading message of `<<"EpsteinZeta.wl"` in Mathematica can now be suppressed using `Quiet`, e.g. `Quiet@<<"EpsteinZeta.wl"`

### Fixed
- Increased edge case parameter in `assignzArgBound` in `src/crandall.c` from $10^{16}$ to `DBL_MAX` from `<float.h>` to ensure that the asymptotic method is not used for large $|\nu|$.
- Added correction term for non-diagonal matrixes in $\nu = d+2 k$ for the regularized Epstein zeta function.
- Removed special use of asymptotic expansion around $\nu = 2$ and $\nu=4$, which reduces the error from $<10^{-11}$ to $<10^{-14}$.
- Fixed special cases `NU` element of $\{\ldots,-4,-2,0,2,3,4,\ldots\}$ for evaluation of `python examples/python/sem_gaussian_1d.py --nu [NU]`.
- Reduced cutoff around zero in `crandall_g` from `src/crandall.c` from $2^{-62}$ to $10^{-64}$. This ensures the correct evaluation close to the singularity.
- Reduced cutoff around $\boldsymbol y^2 = 0$ in `zeta.c` from $2^{-32}$ element wise to $10^{-64}$ for the scalar product. This ensures the correct evaluation close to the singularity.

## [0.4.1] - 2024-12-29
### Fixed
- pypi release

## [0.4.0] - 2024-12-29
### Breaking Changes
- Shared library is now called `libepstein.so`
- Moved Python library from `epsteinlib` nix package to `epsteinlib_python`

### Added
- Build support as meson subproject
- Build support on Windows
- Mathematica interface
- Epstein zeta regularized evaluates at $\nu=-2\mathbb N$
- Github CI

### Fixed
- macOS builds with Python 3.8

## [0.3.0] - 2024-08-19

### Breaking Changes
- Changed function signatures: unified type of `dim` to `unsigned int`

### Changed
- Removed lapack dependency
- Added Kahan summation to improve cutoff error

### Added
- Added Python interface

### Fixed
- Increased stability of gamma function for small values of `x` and negative `nu`

## [0.2.0] - 2024-06-26

### Breaking Changes
- Moved `epsteinZeta.h` to `/include`

### Changed
- Renamed `__epsteinZeta` to `epsteinZetaInternal`
- Removed underscores in front of include guards
- Added `[in]`, `[out]` in documentation to indicate whether vectors may be overwritten
- Changed dynamic memory allocation to static memory allocation to improve performance
- Marked constants as `static`
- Marked unchanged vectors as `const` in `tools.h`

### Added
- Added regularized Epstein zeta function `epsteinZetaReg`

### Fixed
- `egf_gammaStar` works for small `x` and negative `a`

## [0.1.0] - 2024-06-14

_First release._
