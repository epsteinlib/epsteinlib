# SPDX-FileCopyrightText: 2026 Jonathan K. Hartlaub <contact@jkhartlaub.com>
#
# SPDX-License-Identifier: CC0-1.0

# static library, so that gprof/perf see the library functions inside the executable
prof_aniso = executable('epsteinlib_prof_aniso',
    ['prof_aniso.c'],
    include_directories : incdir,
    dependencies: deps,
    install: false,
    override_options: override_options,
    link_with : epsteinlib.get_static_lib()
)
