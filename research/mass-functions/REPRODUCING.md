# Reproducing the mass-function research

The public bundle contains reports, current figures, summary CSVs, a summary
SQLite database, source scripts and provenance. The original raw observations,
model inputs, intermediate numerical draws and complete results remain under
`out/`. They are linked from the ignored `local/` directory.

| Study | Local workspace | Principal local database |
|---|---|---|
| Cluster catalogues and Ember masses | `out/cluster-mass-functions-oct9-v1` | `cluster_catalogues.sqlite` |
| Field population | `out/field-mass-function-oct9-v1` | `field_census.sqlite` |
| Field metallicity tables | `out/field-metallicity-oct9-v1` | — |
| Cluster dynamics and methods | `out/cluster-dissolution-oct9-v1` | `cluster_dissolution.sqlite` |
| Complete 17-model comparison | `out/cluster-dissolution-oct9-v2` | `cluster_followup.sqlite` |
| Stellar fits and held-out BD predictions | `out/cluster-stellar-fit-oct9-v1` | `cluster_stellar_test.sqlite` |
| Membership-selection tests | `out/cluster-membership-oct9-v1` | `cluster_membership_test.sqlite` |

The databases above retain object-level and simulation-level records. Public
`results.sqlite` contains only the compact result tables listed in its `datasets`
table. CSVs are the source for that database. Identifiers remain text when
applicable; raw numerical precision is retained in data files.

## Inspecting or rebuilding the published summary

Only Python's standard library is needed:

```sh
python3 research/mass-functions/verify.py
python3 research/mass-functions/build_database.py
```

The first command checks the file manifest, all public Markdown links, Python
syntax, database integrity, table row counts and equality of database values
with the CSV inputs. It does not rerun the scientific calculations. Rebuilding
the database after changing CSVs requires deliberately refreshing the manifest.

## Running the scientific analyses

The script copies under `analysis/` expect the original sibling workspaces,
which contain their input catalogues and cached tables. They are preserved for
reproduction, rather than exposed as an installed package. Do not launch them
directly from this bundle. Their study READMEs give processing order and the
environment used. Restore missing scripts to their original workspace paths with:

```sh
python3 research/mass-functions/restore_analysis.py --root /path/to/ember
```

This command never overwrites an existing file. `@EMBER_ROOT@` in published
scripts and receipts stands for the chosen checkout path. It is substituted in
restored scripts only. Receipts preserve the original source hashes; packaging
changes have their own checksums in `manifest.json` and `provenance/origins.json`.
Original validation records describe the completed calculations, not a new
validation against the current GitHub solver.

Scientific processing uses NumPy, SciPy, pandas, Matplotlib and Astropy. The
dynamics and selection stages also use galpy, PeTar/FDPS/SDAR/SSE/BSE and
GaiaUnlimited. Version/build records and the Apple-silicon compatibility patch
are in `provenance`. External packages and catalogues retain their own licenses.
No vendor code is redistributed in this folder.

The broad dependency order is: cluster catalogues and tracks → field comparison
and cluster dynamics → ensemble comparison → stellar fit → membership selection.
The field metallicity tables supply the field controls. Several figure scripts
also export reports to the original `docs/research/independent` location; review
their outputs before replacing the curated reports here.

The stellar grid uses the frozen Ember executable based on commit `087351c`
and a private initialization variant documented in the local cluster `physics/`
directory. The selected field histories are listed in
`tables/field/selected_tracks.csv`; the cluster grid is inventoried in
`provenance/clusters/track_inventory.json`. Most tracks begin after deuterium
exhaustion and use an approximate initial-age correction. Atmosphere and
initial-state alternatives are assumptions, not a claim of uniform accuracy
across the whole mass range.

A fresh clone lacks the large source catalogues, atmosphere/EOS inputs, frozen
executables and complete histories needed to rerun everything. Retain or recover
those inputs before using the generation scripts. Existing completed runs should
not be relaunched merely to inspect the figures.

## Retention

The complete catalogue databases, unique simulation states, continuation states,
event logs, full mass-inference draws, numerical checks and executable used for
the dynamics are retained. Duplicate build executables, unused compiler products,
Python caches, stale process markers and superseded presentation copies can be
removed without changing these results. The local cleanup record lives at
`local/cleanup.json`. Old report paths contain pointers to the current reports.
