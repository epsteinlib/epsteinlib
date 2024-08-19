# Changelog

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