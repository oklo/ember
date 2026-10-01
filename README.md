# Ember

Ember is a one-dimensional stellar evolution code for very low-mass stars,
written in C++23. One solver follows Hayashi contraction, hydrogen burning
and helium-white-dwarf cooling.

- Henyey relaxation with analytic derivatives and coupled burning and diffusion.
- Composition-dependent EOS and opacity, non-gray atmospheres, electron conduction
  and mixing-length convection with Ledoux stability.
- Deuterium, pp and CN burning, nuclear screening and plasma-neutrino losses.
- H/He/metal settling, composition-dependent energy accounting, adaptive time
  steps and energy and isotope conservation checks.
- Optional quantum-ion thermodynamics, crystallization and latent heat.
- Double precision and CPU parallelism optimized for Apple silicon.

The code is under development. Full CNO, grain-bearing stellar atmospheres,
element separation during freezing and extreme-cold evolution remain unfinished. Some tables
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
[Paper source](docs/reports/2026-09-27/README.md)

## License

MIT. External physics data and source packages retain their own licenses;
see the READMEs under `data/`.
