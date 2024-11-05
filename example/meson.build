# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

testscript = executable('epsteinlib_c-lattice_sum', 'lattice_sum.c', include_directories : incdir, dependencies: deps, install: not python_only, override_options: override_options, link_with : epsteinlib)
#One could include some tests here
test('Minimal working example', testscript)
