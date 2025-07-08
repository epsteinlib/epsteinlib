# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

project('epsteinlib', 'c', 'cython', version: files('VERSION'), default_options: ['werror=true', 'c_std=gnu11', 'buildtype=release', 'default_library=shared'])

add_project_arguments('-Wno-error=cpp', language: 'c')

cc = meson.get_compiler('c')

deps = []
deps += cc.find_library('m', required : true)

# Initialize source files list
# Populate in subdirectories using zeta_src +=
zeta_src = []

incdir = include_directories('include')

override_options = []
# Enable sanitizer in debug builds
if get_option('debug') == true
    override_options += ['b_sanitize=address']
    add_project_arguments('-DDEBUG', language : 'c')
endif
if get_option('buildtype') == 'release'
    add_project_arguments(['-fno-math-errno'], language: 'c')
endif

build_C = get_option('build_C')
build_python = get_option('build_python')

if not build_C and not build_python
    error('The build_C and build_Python options are disabled. Nothing to build!')
endif

if host_machine.system() == 'windows' and meson.get_compiler('c').get_id() == 'gcc'
    # This avoids having to run 'os.add_dll_directory("/msys64/mingw64/bin")' before importing the python module
    add_project_link_arguments('-static-libgcc', language: 'c')
endif

if build_python
    py_mod = import('python')
    py = py_mod.find_installation(pure: false)

    if host_machine.system() == 'darwin' and py.language_version() == '3.8'
        # This fixes build failure on MacOS with Python 3.8
        add_project_arguments('-Wno-deprecated-declarations', language: 'c')
    endif
endif


subdir('include')
subdir('mathematica')
subdir('src')
subdir('examples/c')
if build_python
    subdir('python')
endif
