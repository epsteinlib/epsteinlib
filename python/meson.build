# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: CC0-1.0

py.extension_module(
    'epsteinlib',
    'epsteinlib.pyx',
    link_whole: epsteinlib.get_static_lib(),
    install: true,
)
# Install stub file
py.install_sources('epsteinlib.pyi')
meson.add_dist_script('remove_symlink.sh')

python_module_build_dir = meson.current_build_dir()
subdir('tests')
