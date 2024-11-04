# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

zeta_src += files('zeta.c', 'gamma.c', 'tools.c', 'crandall.c', 'epsteinZeta.c')
epsteinlib = library('epstein', zeta_src, include_directories : incdir, dependencies: deps, install: true, override_options: override_options)

epsteinlib_dep = declare_dependency(include_directories : incdir, link_with : epsteinlib)

# Generate pkgconfig file
pkgconf = import('pkgconfig')
pkgconf.generate(epsteinlib)

# include unittests
subdir('tests')
