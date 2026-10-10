# Joint stellar-count and brown-dwarf prediction test

Completed on 2026-10-09. This is a conditional population-model check, not a
recovered birth mass function. No collisional jobs remain running.

Report and figure: `docs/research/independent/CLUSTER_STELLAR_TEST_20261009.*`.
The report defines the populations, limitations, selection assumptions and units.
All associated/tail candidates remain the target. The matched central Pleiades
comparison uses DANCe within 3 degrees; the full candidate union is also saved.

Run from the repository root using the existing science interpreter:

```sh
export PYTHONPATH=out/cluster-stellar-fit-oct9-v1/deps:out/cluster-dissolution-oct9-v1/deps:out/cluster-mass-functions-oct9-v1/deps
export MPLCONFIGDIR=/private/tmp/ember-cluster-mpl
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/selection.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/galactic_diffusion.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/build_diffuse.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/build_templates.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/fit_stars.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/audit_fit.py
out/science-env-sept18-v1/bin/python out/cluster-stellar-fit-oct9-v1/make_results.py
```

Generation scripts skip existing diffuse or count templates. Reuse is appropriate
only if their generating code, inputs and assumptions are unchanged. A changed
physical/selection calculation needs fresh templates, not just a fresh fit.
`build_diffuse.py` records costs in each geometry JSON; on a repeat, its aggregate
execution JSON must continue to retain the completed-generation receipts.

Inputs:

- Existing cluster catalogue SQLite, Ember isochrones and 256 mass draws in
  `out/cluster-mass-functions-oct9-v1/`.
- Observed geometry in `out/cluster-dissolution-oct9-v1/` and completed collisional
  endpoints/calibrations in `out/cluster-dissolution-oct9-v2/`.
- GaiaUnlimited source and BSD-3-Clause license in `raw/`. Map: original Zenodo
  8063930 `allsky_M10_hpx7.hdf5`, CC-BY-4.0, published MD5 verified. This models
  Gaia source detection, not cluster-membership recovery. Citation and SHA256
  are in `selection_receipt.json`.
- Private additional dependencies: h5py and astropy-healpix. Exact versions used
  are listed in `dependency_versions.json`; nothing was installed into shared code.

Outputs:

- `stellar_fits.csv`: all stellar objective evaluations; no fitted bound-cluster
  amplitude, only a nonnegative diffuse population amplitude.
- `frozen_stellar_choices.csv`: choices made without brown-dwarf counts. Its hash
  is recorded before generating `brown_dwarf_predictions.csv`.
- `fitted_profiles.csv`: disjoint angular by mass bins; counts are per bin, not
  densities per square degree or per logarithmic mass interval.
- `diagonal_objective_check.csv`: independent diagonal-variance fit. Its scores
  cannot be compared numerically with covariance-weighted scores.
- `catalogue_attribution.csv`: selected catalogue contributions to the full union.
  Sources overlap: never sum their counts as independent populations.
- `cluster_stellar_test.sqlite`: queryable copies plus metadata and receipts.

No formal likelihood, posterior, contamination correction or complete membership
selection is claimed. The map, orbital convergence, count ranges, fixed birth
ratios, full endpoint ages and neglected-gravity checks are recorded in the
receipts. Best-fit diffuse central masses are negligible for the structure screen.
