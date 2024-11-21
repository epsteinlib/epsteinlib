# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
#
# SPDX-License-Identifier: CC0-1.0

test_epsteinZeta = executable('epsteinlib_test_epsteinZeta',
    ['utils.c','test_epsteinZeta.c'],
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

test('test_epsteinZeta',
    test_epsteinZeta,
    workdir: meson.current_source_dir()
)

test('test_crandall',
    test_crandall,
    workdir: meson.current_source_dir()
)
