# Retained Hayashi-start calculation

This local package contains every runtime input for the retained contraction
calculation. `configuration.txt` and the five family manifests use relative
paths. Numerical tables in `tables/` are unchanged copies; their SHA-256 hashes,
sizes, source descriptions, generators and licence notes are recorded in
[`../MANIFEST.json`](../MANIFEST.json). The original inputs remain in place.

Verify a copy before using it:

```sh
python3 scripts/package_runtime_data.py --verify data/MANIFEST.json
```

From the repository root, reproduce the contraction:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 6
build/apps/ember-evolve --lifetime data/production/configuration.txt \
  1000000000 out/contraction --max-steps 2000 --cpu-seconds 1200
```

The output directory must not already exist. A transferred package needs the
whole `production/` directory and `MANIFEST.json`; the manifest alone is not a
download or a complete regeneration recipe. Large numerical payloads remain
local. Their external source terms are not replaced by Ember's code licence.

This package preserves the calculation for reproduction. It does not extend
its physics or supported range. In particular, the retained atmosphere repeats
one composition block and does not support further evolution beyond about
4.34 Gyr. A composition-dependent replacement is being calculated. The current
driver also stops when its whole-star convective approximation ceases to apply.
