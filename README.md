<!--
SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>

SPDX-License-Identifier: AGPL-3.0-only
-->

## EpsteinLib

Experimental alpha version 0.1.0

Authors: Andreas A. Buchheit, Jonathan Busse, Ruben Gutendorf, DevOps: Jan Schmitz

Contact: buchheit@num.uni-sb.de

EpsteinLib is a C library for the fast and efficient computation of the Epstein zeta function for arbitary multidimensional lattices. Originally studied by Epstein [1,2], it forms the basis for the computation of general multidimensional lattice sums in classical and quantum physics applications [3]. Together with its regularization, it serves as the central ingredient in the singular Euler-Maclaurin (SEM) expansion, which generalizes the 300-year-old Euler summation formula to lattice sums in higher dimensions with physically relevant power-law interactions [4-5]. An efficiently computable representation of the Epstein zeta function is provided in [6,7].

For a $d$-dimensional lattice $\Lambda=A\mathbb Z^d$, with $A\in \mathbb R^{d\times d}$ regular, $x,y \in \mathbb R^d$, and $\nu \in \mathbb C$, the Epstein zeta function is defined by the Dirichlet series

$$
Z_{\Lambda,\nu}\begin{vmatrix} x \\ y \end{vmatrix}
= \sum_{z \in \Lambda}{}^{'} \frac{e^{-2\pi i  y \cdot  z}}{\left| x-  z\right|^\nu},\quad \mathrm{Re}(\nu)>d,
$$

which can be meromorphically continued to $\nu \in \mathbb C$. Here, the primed sum excludes the case $z = x.$

The Epstein zeta function is implemented in this library as

```c
double complex epsteinZeta(double nu, int dim, double *A, double *x, double *y);
```
and evalutates to full precision over the whole parameter range for $\nu \in (-10,10)$.

In addition, this library includes the regularized Epstein zeta function, which is analytic around $y=0$, and is defined via 

$$
Z_{\Lambda,\nu}^{\mathrm{reg}}\begin{vmatrix} x \\ y \end{vmatrix} =
e^{2\pi i x\cdot y}
Z_{\Lambda,\nu}\left|\begin{aligned} x \\ y \end{aligned}\right| 
-\frac{\hat{s}(y)}{V_{\Lambda}},
$$

where $V_{\Lambda}=|\det A|$ is the volume of the elementary lattice cell, 

$$
\hat{s}(y)=-\pi^{\nu-\frac{d}{2}}
	\frac{\Gamma((d-\nu)/2)}{\Gamma(\nu/2)}|y|^{\nu-d}
$$

is the distributional Fourier transform of $\vert z \vert^{-\nu}$, where $\Gamma$ denotes the gamma function.

In this library, the regularized Epstein zeta function is included as
```c
double complex epsteinZetaReg(double nu, int dim, double *A, double *x, double *y);
```

## Installation with meson


1. Install meson, ninja, pkg-config, lapacke, e.g. with 
```bash
# Archlinux
pacman -S meson ninja pkgconf lapacke

# MacOS
brew install meson ninja pkg-config lapack
```
2. `cd <path/to/repo>`
3. `meson setup build`
4. `meson compile -C build`
5. To test the library, run `meson test -C build`
### For system-wide installation
Meson supports a system-wide installation of the compiled library. After that, you can use `#include <epsteinZeta.h>` and link with `gcc -lepsteinZeta`. This may require superuser rights.

6. To install system-wide: `meson install -C build`.

7. Try to compile the sample program in `test/lattice_sum.c` with the command `gcc -o lattice_sum lattice_sum.c -lm -lepsteinZeta`. You may encounter the problem that the shared library cannot be found. In this case, you need to modify the environment variables. Please continue with the next step. Otherwise, you are done.

8. Update the environment variables to correctly locate the shared library at runtime. You can find `/path/to/library` in the output given by `meson install`.
```bash
# Linux
export $LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/library

# MacOS
export $DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/path/to/library
```

### For local installation
6. Copy the header `include/epsteinZeta.h` and move the compiled library `build/src/libepsteinZeta.so` to places of your choice, e. g.
```bash
cp include/epsteinZeta.h /your/path/to/include
mv build/src/libepsteinZeta.so /your/path/to/library
```
7. To test your library, try to compile `test/lattice_sum.c` with the command `gcc -o lattice_sum lattice_sum.c -lm -L/your/path/to/library -lepsteinZeta -I/your/path/to/include`.

## View api documentation

1. install [doxygen](https://www.doxygen.nl/manual/install.html)
2. `cd <path/to/repo>`
3. `doxygen`
4. open `html/index.html` in browser

## Usage

Minimal working examples for calculating the Madelung constant in $3$ dimensions.
### in c
``` c
// If the library is installed, compile with `gcc -o lattice_sum lattice_sum.c -lm -lepsteinZeta 
// If the library is not installed, compile with `gcc -o lattice_sum lattice_sum.c -lm -L/path/to/library -lepsteinZeta -I/path/to/include

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "epsteinZeta.h"


int main() {

    // Madelung constant found in literature
    double madelungRef = -1.7475645946331821906362120355443974;
    int dim = 3;
    double m[] = {1, 0, 0, 0, 1,
                  0, 0, 0, 1};    // identity matrix for whole numbers
    double x[] = {0, 0, 0};       // no shift
    double y[] = {0.5, 0.5, 0.5}; // alternating sum
    double nu = 1.0;
    double madelung = creal(epsteinZeta(nu, dim, m, x, y));
    printf("Madelung sum in 3 dimensions:\t %.16lf\n", creal(madelung));
    printf("Reference value:\t\t %.16lf\n", madelungRef);
    printf("Relative error:\t\t\t +%.2e\n",
           fabs(madelungRef - madelung) / fabs(madelungRef));

    return fabs(madelung - madelungRef) > pow(10, -14);
}
```


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## References

[1] P. Epstein. “Zur Theorie allgemeiner Zetafunctionen”. Math. Ann. 56 (1903), pp. 615–644.

[2] P. Epstein. “Zur Theorie allgemeiner Zetafunktionen. II”. Math. Ann. 63 (1906), pp. 205–216

[3] Andreas A. Buchheit et al. “Exact Continuum Representation of Long-range Interacting Systems and Emerging Exotic Phases in Unconventional Superconductors”, Phys. Rev. Research 5, 043065 (2023)

[4] Andreas A Buchheit and Torsten Keßler. “On the Efficient Computation of Large Scale Singular Sums with Applications to Long-Range Forces in Crystal Lattices”. J. Sci. Comput. 90.1 (2022), pp. 1–20

[5] Andreas A Buchheit and Torsten Keßler. “Singular Euler–Maclaurin expansion on multidimensional lattices”. Nonlinearity 35.7 (2022), p. 3706

[6] R. Crandall. “Unified algorithms for polylogarithm, L-series, and zeta variants”. Algorithmic Reflections: Selected Works. PSIpress, 2012

[7] Andreas A. Buchheit, Torsten Keßler, and Kirill Serkh. "On the computation of lattice sums without translational invariance". arXiv preprint arXiv:2403.03213.
