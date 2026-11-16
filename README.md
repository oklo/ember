# Ember

Ember is a one-dimensional stellar evolution code for very low-mass stars,
written in C++23. The C++ source is written to trace lineage back to the
IBM 7094-compatible [FORTRAN IV reconstruction](https://oklo.github.io/Henyey/)
of the algorithms described in
[Henyey et al. (1959)](https://adsabs.harvard.edu/pdf/1959ApJ...129..628H) and
[Henyey, Forbes & Gould (1964)](https://adsabs.harvard.edu/pdf/1964ApJ...139..306H),
thus providing continuity to the dawn of modern computing.

The solver couples Henyey relaxation to burning, mixing and diffusion, using
composition-dependent EOS, opacity and atmosphere inputs. It includes deuterium,
pp and CN burning, electron conduction, settling and experimental cold-remnant
physics. The [code guide](docs/DEVELOPMENT.md) identifies the main routines.

The code is under development. Low-mass boundary studies use external clear
and cloudy atmosphere tables. Full CNO, element separation during freezing
and extreme-cold evolution are under current development. Some tables
required for the cooling calculations are still local.
See [physics and limitations](docs/COLD_REMNANT.md).

## Build

Requires CMake, Ninja and a C++23 compiler.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build -LE material-data --output-on-failure
```

This runs the tests that do not require external material fixtures. To run
the full suite, install the [test inputs](tests/README.md) and omit
`-LE material-data`. Stellar runs also require the
[physics inputs](docs/DATA_REPRODUCTION.md); a fresh clone does not contain
all inputs used by the development calculations. The solver uses CPU threads;
it does not currently use Metal or the GPU.

[Autonomous control](docs/AUTONOMY.md) ·
[Run configuration](docs/LIFETIME_DRIVER.md) ·
[Physics inputs](docs/DATA_REPRODUCTION.md) ·
[Numerical lineage](docs/LINEAGE.md)

## License

MIT. External physics data and source packages retain their own licenses;
see the READMEs under `data/`.
