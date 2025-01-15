<!--
SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>

SPDX-License-Identifier: AGPL-3.0-only
-->

# Changelog

## [0.5.0] - unreleased
### Breaking Changes

### Added

### Fixed

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
