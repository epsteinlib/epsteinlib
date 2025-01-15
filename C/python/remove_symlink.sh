#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: CC0-1.0

set -eu
cd "$MESON_DIST_ROOT/python"
rm epsteinlib.pyx
cp _epsteinlib.py epsteinlib.pyx
