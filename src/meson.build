# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

zeta_src = ['zeta.c', 'gamma.c', 'tools.c', 'crandall.c', 'epsteinZeta.c']
epsteinZetalib = library('epsteinZeta', zeta_src, include_directories : incdir, dependencies: deps, install: true, override_options: override_options)
