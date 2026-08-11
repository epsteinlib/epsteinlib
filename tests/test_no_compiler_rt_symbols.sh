#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
#
# SPDX-License-Identifier: AGPL-3.0-only

# Guards against compiler-emitted helper functions unavailable when
# cross-compiling for aarch64-apple-darwin (as Yggdrasil does for the Julia
# bindings). GCC resolves these silently from libgcc, so they otherwise only
# surface as link errors. The cause is complex-by-complex division which
# we avoid.
set -euo pipefail

lib="$1"
forbidden='__divdc3|__divsc3'

found=0
for obj in "$lib".p/*.o; do
    [ -e "$obj" ] || continue
    syms=$(nm "$obj" | grep -Eo "$forbidden" | sort -u || true)
    [ -n "$syms" ] || continue
    found=1
    for sym in $syms; do
        echo "$(basename "$obj"): references $sym"
        locs=$(objdump -dlr "$obj" 2>/dev/null | awk -v s="$sym" '
            /[^ \/]*\.[ch]:[0-9]+/ { loc = $0 }
            index($0, s) && loc != "" { print loc }
        ' | sed 's/ (discriminator [0-9]*)//' | sort -u || true)
        if [ -n "$locs" ]; then
            echo "$locs" | sed 's/^/    /'
        else
            echo "    (no source locations -- rebuild with debug info:"
            echo "     meson setup --reconfigure -Ddebug=true build/)"
        fi
    done
done

if [ "$found" -ne 0 ]; then
    echo "ERROR: complex-by-complex division will not link on darwin."
    exit 1
fi

echo "OK: no compiler-rt helper symbols in $lib"
