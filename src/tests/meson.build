# SPDX-FileCopyrightText: 2025 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
#
# SPDX-License-Identifier: CC0-1.0

test_epsteinZeta = executable('epsteinlib_test_epsteinZeta',
    ['../tools.c','utils.c','test_epsteinZeta.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_crandall = executable('epsteinlib_test_crandall',
    ['utils.c','test_crandall.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_epstein = executable('epsteinlib_benchmark_epstein',
    ['benchmark_epstein.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_gamma = executable('epsteinlib_benchmark_gamma',
    ['../gamma.c','benchmark_gamma.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test('test_epsteinZeta',
    test_epsteinZeta,
    workdir: meson.current_source_dir()
)

test('test_crandall',
    test_crandall,
    workdir: meson.current_source_dir()
)
