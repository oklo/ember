# Testing

Build and run tests that need no external material fixtures:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build -LE material-data --output-on-failure
```

The `material-data` tests additionally use generated EOS and atmosphere tables,
archived source calculations, and the configured production input family.
`data/TEST_DATA_MANIFEST.json` pins the source fixtures by size and SHA-256.
Given a bundle made from those exact inputs:

```sh
python3 scripts/package_test_data.py --install /path/to/test-data.tar.gz
python3 scripts/package_test_data.py --verify
ctest --test-dir build --output-on-failure
```

The lifetime and screened-envelope integration tests also need the inputs
named by `data/production/configuration.txt`. The complete bundle is not
currently distributed; see [physics inputs](../docs/DATA_REPRODUCTION.md).
Missing data is a failed integration test, never a physical pass.

To create the source-fixture bundle from an existing complete checkout:

```sh
python3 scripts/package_test_data.py --pack /path/to/test-data.tar.gz
```

Small numerical fixtures and source-regression records remain in `tests/data/`.
Synthetic fixtures exercise numerical behavior; they are not stellar inputs.
