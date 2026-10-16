<!--
SPDX-FileCopyrightText: 2026 Jonathan K. Hartlaub <contact@jkhartlaub.com>

SPDX-License-Identifier: AGPL-3.0-only
-->

# Profiling epsteinZetaAniso

`prof_aniso.c` drives `epsteinZetaAniso` in two cases:

- `zero`: identity lattice, x = y = 0, α = (n, 0, …, 0)
  (odd n is a symmetry zero and returns immediately),
- `generic`: oblique lattice (0.3 on the upper off-diagonal), generic x and y,
  α spread over all directions, e.g. n = 11 in 2D gives α = (6, 5).

Modes:

```sh
epsteinlib_prof_aniso <dim> <n> <zero|generic> <iterations> [nu]   # one case, for gprof/callgrind/perf
epsteinlib_prof_aniso sweep [seconds per configuration] [nu]      # timing table d = 2, 3, both cases
```

The library is linked statically, so the profilers see its functions inside the
executable. Run all commands from the repository root; profile output files
land in the current directory.

## Timing table (before/after comparisons)

```sh
meson setup build-prof -Dprofiling=true -Dbuild_python=false -Dc_args='-g -fno-omit-frame-pointer'
meson compile -C build-prof
./build-prof/profiling/epsteinlib_prof_aniso sweep 0.3 > sweep_before.txt
```

The values block at the end of the output catches accuracy changes: diff it
against a reference run.

## callgrind: deterministic instruction counts, line level

Preferred for diagnosis. Counts do not depend on machine load, include libm,
and attribute inlined code to source lines.

```sh
valgrind --tool=callgrind --callgrind-out-file=cg.out \
    ./build-prof/profiling/epsteinlib_prof_aniso 2 20 generic 20
callgrind_annotate --inclusive=yes cg.out | head -40   # functions incl. callees
callgrind_annotate cg.out src/harmonics.c | less       # per source line
```

Cache misses (array access cost), with the same binary:

```sh
valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=cache.out \
    ./build-prof/profiling/epsteinlib_prof_aniso 3 16 generic 5
cg_annotate cache.out src/harmonics.c | less
```

## gprof

Inlining of called-once and small functions is switched off so that
`harmonic_h_outer_double` and `harmonic_h_dd` show up with call counts.
Caveats: libm (`cexp`, `pow`, `fma`, `tgamma`, …) is not instrumented and its
time is missing or attributed to `_init`; the `mcount` overhead inflates small,
frequently called functions such as `hpdyad_*`.

```sh
meson setup build-gprof -Dprofiling=true -Dbuild_python=false -Db_lto=false \
    -Dc_args='-pg -g -fno-inline-functions-called-once -fno-inline-small-functions' \
    -Dc_link_args='-pg'
meson compile -C build-gprof
./build-gprof/profiling/epsteinlib_prof_aniso 2 20 generic 1500   # keep wall time > ~2 s
gprof -b ./build-gprof/profiling/epsteinlib_prof_aniso gmon.out | less
```

## perf (if available)

```sh
perf record -g --call-graph fp -- ./build-prof/profiling/epsteinlib_prof_aniso 3 16 generic 200
perf report --sort=dso,symbol
```
