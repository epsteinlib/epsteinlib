#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Jonathan K. Hartlaub <contact@jkhartlaub.com>
#
# SPDX-License-Identifier: CC0-1.0

# Profiling epsteinZetaAniso. Run from the repository root.
# case: zero (identity, x = y = 0, alpha = (n,0,..)) | generic (oblique, x, y, alpha spread)
#
#   profiling/prof.sh setup                        once, and after code changes (rebuilds all three)
#   profiling/prof.sh sweep before.txt             timing table 2D/3D x zero/generic + values;
#                                                  run before and after a change, diff the files
#   profiling/prof.sh cg 3 16 generic 3            callgrind: where do the instructions go,
#                                                  deterministic, includes libm -> cg.out
#   profiling/prof.sh line src/harmonics.c         which lines of a file are hot (needs cg.out)
#   callgrind_annotate --tree=both cg.out | less   who calls whom, with costs
#   kcachegrind cg.out                             same as clickable GUI
#   profiling/prof.sh gprof 2 20 generic 1500      gprof flat profile + call graph -> gprof.txt;
#                                                  misses libm time, keep run > ~2 s
#
# Iterations for cg: native run ~0.05 s (valgrind ~50x slower), e.g. 3D n=16 -> 3, 2D n=20 -> 10.

set -euo pipefail
P=profiling/epsteinlib_prof_aniso
COMMON=(-Dprofiling=true -Dbuild_python=false)
NOINL='-fno-inline-functions-called-once -fno-inline-small-functions'

case "${1:-}" in
setup)
    [ -d build-prof ] || meson setup build-prof "${COMMON[@]}" -Dc_args='-g'
    # no inlining, so double and double-double show up as separate functions
    [ -d build-noinl ] || meson setup build-noinl "${COMMON[@]}" -Dc_args="-g $NOINL"
    # gprof needs -pg in both; valgrind cannot run -pg binaries
    [ -d build-gprof ] || meson setup build-gprof "${COMMON[@]}" -Db_lto=false \
        -Dc_args="-pg -g $NOINL" -Dc_link_args='-pg'
    for b in build-prof build-noinl build-gprof; do meson compile -C $b; done
    ;;
sweep)
    ./build-prof/$P sweep 0.3 | tee "${2:-sweep.txt}"
    ;;
cg)
    shift
    valgrind --tool=callgrind --callgrind-out-file=cg.out ./build-noinl/$P "$@" > /dev/null 2>&1
    callgrind_annotate --inclusive=yes cg.out | sed -n '/file:function/,$p' | head -40
    ;;
line)
    callgrind_annotate cg.out "${2:-src/harmonics.c}" | less
    ;;
gprof)
    shift
    ./build-gprof/$P "$@"
    gprof -b ./build-gprof/$P gmon.out > gprof.txt
    head -30 gprof.txt
    echo "full profile in gprof.txt"
    ;;
*)
    sed -n '6,21p' "$0"
    exit 1
    ;;
esac
