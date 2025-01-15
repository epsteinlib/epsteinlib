# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

python_only = not build_C and build_python

zeta_src += files('zeta.c', 'gamma.c', 'tools.c', 'crandall.c', 'epsteinZeta.c')
epsteinlib = both_libraries('epstein', zeta_src, include_directories : incdir, dependencies: deps, install: not python_only, override_options: override_options)

epsteinlib_dep = declare_dependency(include_directories : incdir, link_with : epsteinlib)

if not python_only
  # Generate pkgconfig file
  pkgconf = import('pkgconfig')
  pkgconf.generate(epsteinlib)
endif

# include unittests
subdir('tests')
