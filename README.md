<!--
SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>

SPDX-License-Identifier: AGPL-3.0-only
-->

# EpsteinLib

Authors: Andreas A. Buchheit, Jonathan Busse, Ruben Gutendorf, DevOps: Jan Schmitz

Contact: buchheit@num.uni-sb.de

EpsteinLib is a C library designed for the fast and efficient computation of the Epstein zeta function for arbitrary multidimensional lattices. In addition to the C library, we also offer a Python package, epsteinlib, which can be easily installed via pip.

Originally studied by Epstein [1,2], the Epstein zeta function forms the basis for computing general multidimensional lattice sums in classical and quantum physics applications [3]. Together with its regularization, it serves as the central ingredient in the singular Euler-Maclaurin (SEM) expansion, which generalizes the 300-year-old Euler summation formula to lattice sums in higher dimensions with physically relevant power-law interactions [4-5]. An efficiently computable representation of the Epstein zeta function is provided in [6,7].

For a $d$-dimensional lattice $\Lambda=A\mathbb Z^d$, with $A\in \mathbb R^{d\times d}$ regular, $\boldsymbol x,\boldsymbol y \in \mathbb R^d$, and $\nu \in \mathbb C$, the Epstein zeta function is defined by the Dirichlet series

$$
Z_{\Lambda,\nu}\begin{vmatrix} \boldsymbol x \newline\boldsymbol y \end{vmatrix}
= \sum_{z \in \Lambda}{}^{'} \frac{e^{-2\pi i \boldsymbol y \cdot \boldsymbol z}}{\left| \boldsymbol x- \boldsymbol z\right|^\nu},\quad \mathrm{Re}(\nu)>d,
$$

which can be meromorphically continued to $\nu \in \mathbb C$. Here, the primed sum excludes the case $\boldsymbol z = \boldsymbol x.$

The Epstein zeta function is implemented in this library as

```c
double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
```
In the Python package, it is implemented as
```python
def epstein_zeta(nu: float | int, A: NDArray[np.float64], x: NDArray[np.float64], y: NDArray[np.float64]) -> complex
```
In the Mathematica package, it is implemented as
```mathematica
EpsteinZeta[\[Nu],A,x,y]
```
and evaluates to full precision over the whole parameter range up to ten dimensions.

In addition, this library includes the regularized Epstein zeta function, which is analytic around $\boldsymbol y=0$, and is defined via

$$
Z_{\Lambda,\nu}^{\mathrm{reg}}\begin{vmatrix} \boldsymbol x \newline\boldsymbol y \end{vmatrix} =
e^{2\pi i \boldsymbol x\cdot\boldsymbol y}
Z_{\Lambda,\nu}\left|\begin{aligned} \boldsymbol x \newline\boldsymbol y \end{aligned}\right|
-\frac{\hat{s}(\boldsymbol y)}{V_{\Lambda}},
$$

where $V_{\Lambda}=|\det A|$ is the volume of the elementary lattice cell, and

$$
\hat{s}_\nu(\boldsymbol y) = \frac{\pi^{\nu/2}}{\Gamma(\nu/2)}\Gamma\big((d-\nu)/2\big)  (\pi \boldsymbol y^2)^{(\nu - d)/2},\quad \nu \not\in (d+2\mathbb N_0)
$$

is the distributional Fourier transform of $\vert\boldsymbol z \vert^{-\nu}$, where $\Gamma$ denotes the gamma function and we adopt the choice

$$
\hat s_{d+2k}(\boldsymbol y)= \frac{\pi^{k+d/2}}{\Gamma(k+d/2)}\frac{(-1)^{k+1}}{k!} ( \pi \boldsymbol y^2 )^{k} \log (\pi  \boldsymbol y^{2}),\quad k\in \mathbb N_0.
$$

In the c library, the regularized Epstein zeta function is included as
```c
double complex epsteinZetaReg(double nu, unsigned int dim, const double *A, const double *x, const double *y);
```
in the Python package as
```python
def epstein_zeta_reg(nu: float | int, A: NDArray[np.float64], x: NDArray[np.float64], y: NDArray[np.float64]) -> complex
```
and in the Mathematica package as
```mathematica
EpsteinZetaReg[\[Nu],A,x,y]
```
## Installation
Install our required dependencies: meson, ninja, pkg-config, python3 e.g. with
```bash
# Archlinux
pacman -S meson ninja pkgconf python

# MacOS
brew install meson ninja pkg-config python3
```

Currently, we support native Windows builds only with GCC installed via https://www.msys2.org/. Different environments may or may not work.
However, for the full out of the box development experience
we encourage Windows users to use [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and follow the Linux installation instructions.

### Installing only the Python wrapper with pip
```bash
# Create a virtualenvironment and activate it if you are not already inside one.
python3 -m venv .venv && source .venv/bin/activate
# Install epsteinlib
python -m pip install epsteinlib
```

### Installing the C library and the Python wrapper with meson
1. git clone https://github.com/epsteinlib/epsteinlib.git
2. `cd epsteinlib`
3. `meson setup build`
4. `meson compile -C build`
5. To test the library, run `meson test -C build`

Proceed either with system-wide or local installation

<details><summary>System-wide installation</summary>

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
</details>

<details><summary>Local installation</summary>

6. Copy the header `include/epsteinZeta.h` and move the compiled library `build/src/libepsteinZeta.so` to places of your choice, e. g.
```bash
cp include/epsteinZeta.h /your/path/to/include
mv build/src/libepsteinZeta.so /your/path/to/library
```
7. To test your library, try to compile `test/lattice_sum.c` with the command `gcc -o lattice_sum lattice_sum.c -lm -L/your/path/to/library -lepsteinZeta -I/your/path/to/include`.

</details>

## View api documentation
See https://epsteinlib.github.io/epsteinlib/

## Usage

Minimal working examples for calculating the Madelung constant in $3$ dimensions.
### in C
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
    unsigned int dim = 3;
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
### in Python

```py
import numpy as np
from epsteinlib import epstein_zeta

madelung_ref = -1.7475645946331821906362120355443974
dim = 3
a = np.identity(dim)            # identity matrix for whole numbers
x = np.zeros(dim)               # no shift
y = np.full(dim, 0.5)           # alternating sum
nu = 1.0
madelung = np.real(epstein_zeta(nu, a, x, y))
print(f"Madelung sum in 3 dimensions:\t {madelung:.16f}")
print(f"Reference value:\t\t {madelung_ref:.16f}")
print(f"Relative error:\t\t\t +{abs(madelung_ref - madelung) / abs(madelung_ref):.2e}")
```

In the `examples/python/` folder, you can find two more Python examples:

1. `dispersion_relation_3d.py`: This script demonstrates how to use EpsteinLib to calculate quantum dispersion relations in 3D.
2. `sem_gaussian_1d.py`: This script showcases the Singular Euler-Maclaurin (SEM) expansion for a Gaussian function in 1D. It has an optional argument `--nu` that can be used to set the value of nu. For example, you can run it with `python sem_gaussian_1d.py --nu 1`. If no value is provided, it defaults to nu = 1.5.

These examples, along with the `lattice_sum.py` script, provide a comprehensive overview of how to use EpsteinLib in various scenarios.



### in Mathematica

```mathematica
<<"EpsteinZeta.wl"

madelungRef = -1.7475645946331821906362120355443974;

dim = 3;
A = IdentityMatrix[dim];
x = ConstantArray[0, dim];
y = ConstantArray[0.5, dim];
\[Nu] = 1.0;

madelung = Re[EpsteinZeta[\[Nu], A, x, y]];

Print["Madelung sum in 3 dimensions: ", NumberForm[madelung, 16]];
Print["Reference value:              ", NumberForm[madelungRef, 16]];
Print["Relative error:               +", ScientificForm[Abs[madelungRef - madelung]/Abs[madelungRef], 2]];
```

Executing this code snipped in the same folder as `EpsteinZeta.wl` and setting `SetDirectory[NotebookDirectory[]]` is the easiest way to help mathematica find the package.

## Development environment
We provide a nix devshell to have a reproducible development environment with the same dependencies across different operating systems. Once you have installed and configured nix starting developing is as easy as running `nix develop`.

<details><summary>Nix installation instructions</summary>

### Nix based - recommended
1. Install [nix](https://nixos.org/download/); Follow the [wiki](https://wiki.archlinux.org/title/Nix)
2. Configure nix by executing
```bash
sudo tee -a /etc/nix/nix.conf <<CFG
max-jobs = auto
#max-jobs = 1
experimental-features = nix-command flakes auto-allocate-uids
auto-allocate-uids = true
auto-optimise-store = true
CFG
```
3. `systemctl enable --now nix-daemon.socket`
4. `usermod -a -G nix-users <your username>`
5. Reboot
6. `cd <path/to/repo>`
7. `nix develop` or `nix run -- <your args>`

### Nix-Portable based - if you do not have root rights
1. Install [nix-portable](https://github.com/DavHau/nix-portable):
```bash
mkdir -p ~/.local/bin
cd ~/.local/bin

curl -L https://github.com/DavHau/nix-portable/releases/latest/download/nix-portable-$(uname -m) > ./nix-portable
chmod +x ./nix-portable
cat > ./nix <<NIX
#!/usr/bin/env bash
CURDIR=\$(dirname "\$(readlink -f "\$0")")
NP_RUNTIME=bwrap "\$CURDIR/nix-portable" nix \$@
NIX
chmod +x ./nix

export PATH=~/.local/bin:"$PATH"
cd ~
nix run 'nixpkgs#hello'

```
2. Configure nix.conf by executing
```bash
tee -a ~/.nix-portable/conf/nix.conf <<CFG
max-jobs = auto
#max-jobs = 1
auto-optimise-store = true
CFG
```
3. Add .local/bin permanently to your PATH
```bash
echo 'PATH=$HOME/.local/bin:"$PATH"' >> ~/.env
echo 'export $(envsubst < .env)' | tee -a .bashrc >> .zshrc
```

4. `cd <path/to/repo>`
5. `nix develop` or `nix run -- <your args>`
</details>

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
