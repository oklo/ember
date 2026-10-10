# Cluster dissolution and catalogue selection

Private analysis workspace for the Hyades, Pleiades and Ember field comparison.
Scientific report: `docs/research/independent/CLUSTER_DISSOLUTION_RESULTS_20261009.md`
relative to the Ember root. The PNG/PDF stem is `CLUSTER_DISSOLUTION_20261009`.

Run commands from the Ember root with this environment:

```sh
export PYTHONPATH=out/cluster-dissolution-oct9-v1/deps:out/cluster-mass-functions-oct9-v1/deps
export MPLCONFIGDIR=/private/tmp/ember-cluster-mpl
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v1/analyze_clusters.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v1/plot_results.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v1/archive_results.py
```

Before starting a new run, inspect existing `runs/*/pid`, continuation PIDs and
receipts. The initial-condition driver refuses an existing run directory:

```sh
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v1/run_cluster.py NEW_NAME --cluster hyades --n 3000 --rh 4 --seed 41 --threads 4 --output 25
```

An assortative binary control adds `--binary-fraction .25 --binary-mode
assortative`. That mode actually assigns separate binary fractions of .15 to
brown dwarfs and .30 to stars; the common .25 switch is merely nonzero metadata
for the control. Exact component counts and parameters are in `initial.json`.
Random pairing is a stress test, not the preferred binary model.

The `cluster_dissolution.sqlite` file contains run configurations, every initial
component, surviving final descendants, dynamical histories, mock observations,
the observed aperture counts and execution cost. Initial IDs absent from final
survivors must be checked against the BSE events for mergers/zero-mass removal;
they were not deleted as distant escapers. Full event logs and sparse snapshots
are retained alongside the ledger.

`vendor_receipt.json`, `build_command.json` and `petar_m4_compatibility.patch`
record the local CPU build. PeTar uses scalar ARM gravity with OpenMP, not CUDA.
GSL is local. The GCC 16 OpenMP include order and the Fortran runtime link are
build compatibility changes. Shared Ember source, tables and Opus work are
untouched. The shell initial-data helper is invoked with explicit unit scaling
to avoid its uninitialized AWK variable on macOS.

The compact Hyades model continues from its saved 150-Myr state after a
25-Myr comparison of faster timestep/changeover settings. `continuation.json`
identifies the saved state; its `continuation/` directory stores the later
snapshots and separate cost. The one-time restart script must not be rerun.

The larger-step experiments that set the neighbor radius from the velocity
dispersion of *all* descendants were discarded: fast escapers inflate that
dispersion and cause unnecessarily large encounter groups. The adopted setting
instead keeps the search radius appropriate to the cluster and retains each
particle's own velocity-dependent neighbor search. No force softening or
escaper deletion is introduced.

Photometric thresholds and footprints are approximate sensitivity brackets.
The database must not be presented as a calibrated completeness-corrected
cluster birth mass function. See `survey_selection_scope.json` and the report.
