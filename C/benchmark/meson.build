# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
#
# SPDX-License-Identifier: CC0-1.0

benchmark_epstein = executable('epsteinlib_benchmark_epstein',
    ['benchmark_epstein.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_gamma = executable('epsteinlib_benchmark_gamma',
    ['benchmark_gamma.c'],
    include_directories : [incdir, '../src'],
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)
