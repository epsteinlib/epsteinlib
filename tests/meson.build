# SPDX-FileCopyrightText: 2025 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
# SPDX-FileCopyrightText: 2025-2026 Jonathan Busse <jonathan@jbusse.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-License-Identifier: CC0-1.0

test_epsteinZeta = executable('epsteinlib_test_epsteinZeta',
    ['../src/tools.c','utils.c','test_epsteinZeta.c'],
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

test_setZetaDer = executable('epsteinlib_test_setZetaDer',
    ['wrappers.c','../src/tools.c','utils.c','test_setZetaDer.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_epsteinZetaRegDer = executable('epsteinlib_test_epsteinZetaRegDer',
    ['wrappers.c','../src/tools.c','utils.c','test_epsteinZetaRegDer.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_tools= executable('epsteinlib_test_tools',
    ['../src/tools.c','utils.c','test_tools.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_hpdyad= executable('epsteinlib_test_hpdyad',
    ['../src/hpdyad.c','utils.c','test_hpdyad.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_harmonics= executable('epsteinlib_test_harmonics',
    ['wrappers.c','utils.c','test_harmonics.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

test_crandall_der = executable('epsteinlib_test_crandall_derivatives',
    ['wrappers.c','utils.c','test_crandall_der.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_epstein = executable('epsteinlib_benchmark_epstein',
    ['utils.c', 'benchmark_epstein.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_gamma = executable('epsteinlib_benchmark_gamma',
    ['utils.c','../src/gamma.c','benchmark_gamma.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib
)

benchmark_aniso_timing = executable('epsteinlib_benchmark_aniso_timing',
    ['wrappers.c','utils.c','benchmark_aniso_timing.c'],
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

test('test_setZetaDer',
    test_setZetaDer,
    workdir: meson.current_source_dir()
)

test('test_epsteinZetaRegDer',
    test_epsteinZetaRegDer,
    workdir: meson.current_source_dir()
)

test('test_harmonics',
    test_harmonics,
    workdir: meson.current_source_dir()
)

test('test_crandall_derivatives',
    test_crandall_der,
    workdir: meson.current_source_dir()
)

test('test_tools',
    test_tools,
    workdir: meson.current_source_dir()
)

test('test_hpdyad',
    test_hpdyad,
    workdir: meson.current_source_dir()
)
