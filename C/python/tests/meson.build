# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: CC0-1.0

fs = import('fs')

unittests = files(
    'test_epsteinlib.py',
)

test_env = environment()
test_env.prepend('PYTHONPATH', python_module_build_dir)

foreach unittest : unittests
    test(
        fs.stem(unittest),
        py,
        args : [unittest],
        env : test_env,
    )
endforeach
