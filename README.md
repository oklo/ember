# Ember

Ember is a one-dimensional stellar evolution code for very low-mass stars,
written in C++23. Its numerical lineage is **Henyey FORTRAN IV → Bodenheimer
FORTRAN 77 → Astra/Opus Ember C++**. The [lineage notes](docs/LINEAGE.md)
distinguish the inherited methods from the local historical reconstructions.

The main research products are exploratory tracks across the low-mass spectrum,
especially at 0.1 solar masses, and studies of the hydrogen-burning boundary and
late reheating versus metallicity. Start with the [research record](research/record/README.md)
for retained data, qualifications and reproduction commands. Scientific campaigns
are paused as of November 1, 2026.

The solver couples Henyey relaxation to burning, mixing and diffusion, using
composition-dependent EOS, opacity and atmosphere inputs. It includes deuterium,
pp and CN burning, electron conduction, settling and experimental cold-remnant
physics. The [code guide](docs/DEVELOPMENT.md) identifies the main routines.

The code is under development. Low-mass boundary studies use external clear
and cloudy atmosphere tables. Full CNO, element separation during freezing
and extreme-cold evolution remain unfinished. Some tables
required for the cooling calculations are still local.
See [physics and limitations](docs/COLD_REMNANT.md).

## Build

Requires CMake, Ninja and a C++23 compiler.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Stellar runs and table-dependent tests also require the
[physics inputs](docs/DATA_REPRODUCTION.md); a fresh clone does not contain
all inputs used by the development calculations. The solver uses CPU threads;
it does not currently use Metal or the GPU.

[Run configuration](docs/LIFETIME_DRIVER.md) ·
[Working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) ·
[Paper source](docs/reports/2026-09-27/README.md) ·
[Cluster and field mass functions](research/mass-functions/README.md) ·
[Hydrogen-burning boundary](research/metallicity-boundary/README.md)

## License

MIT. External physics data and source packages retain their own licenses;
see the READMEs under `data/`.
