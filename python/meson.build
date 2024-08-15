# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: CC0-1.0

py_mod = import('python')
py = py_mod.find_installation(pure: false)

py.extension_module(
    'epsteinlib',
    'epsteinlib.pyx',
    zeta_src,
    dependencies: deps,
    include_directories: zeta_include,
    install: true,
)
# Install stub file
py.install_sources('epsteinlib.pyi')

python_module_build_dir = meson.current_build_dir()
subdir('tests')
