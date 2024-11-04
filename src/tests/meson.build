# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
#
# SPDX-License-Identifier: CC0-1.0

test_epsteinZeta = executable('epsteinlib_test_epsteinZeta',
    'test_epsteinZeta.c',
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test('test_epsteinZeta',
    test_epsteinZeta,
    workdir: meson.current_source_dir()  # This makes the test run from the tests/ directory
)
