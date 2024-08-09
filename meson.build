# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: CC0-1.0

project('epstein', 'c', version:'0.1.0', default_options: ['werror=true', 'c_std=gnu11', 'buildtype=release', 'default_library=shared'])

add_global_arguments('-Wno-error=cpp', language: 'c')

cc = meson.get_compiler('c')

deps = []
deps += cc.find_library('m', required : true)

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

subdir('include')
subdir('src')
subdir('test')

