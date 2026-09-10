# Ember working paper — September 10, 2026

- [Working paper PDF](ember_status_and_future.pdf)
- [LaTeX source](ember_status_and_future.tex)
- [Evolution against age](evolution_history.pdf) and [Ember data](evolution_history.csv)
- [Comparison at equal hydrogen abundance](f77_matched_hydrogen.pdf)
- [F77 figure data](f77_history.csv) and [extraction checks](f77_history_provenance.json)
- [Validation receipt](validation.json)
- [Local recovery manifest](recovery_manifest.json)

The paper describes the completed 512-point, 0.1 Msun calculation through 3.40 trillion
years, refined input families, independent F77 comparisons, performance work,
and the remaining requirements for a helium-remnant cooling track to 100 K.
The star is still fully convective and burning hydrogen.

The trial table of 108 atmospheres still needs composition refinement.
A stellar trial using it began losing full convection near 3.543 trillion
years and stopped near 3.544 trillion years at the program's limit on rejected
time steps. The trial is discussed separately from the plotted 3.40-trillion-year
reference. For later job status, read the newest [handoff entry](../../../HANDOFF.md).

## Rebuild

From this directory:

```sh
python3 build_figures.py
tectonic ember_status_and_future.tex
```

Both figures use the committed CSV files and F77 extraction record, and require
NumPy and Matplotlib. They can be rebuilt without raw stellar archives.
F77 luminosity and radius are converted to the same nominal solar units as
Ember. Thin F77 lines have opacity 0.30; individual segments changing a surface
quantity by more than 5%, or an abundance by more than 0.01, have opacity 0.08.
No data are deleted or smoothed, and no claim is made that every abrupt change
is numerical error. With local recovery files present,
`python3 build_figures.py --from-archives` checks and joins the fresh
0–3.30T history and exact 3.30–3.40T restart, then rebuilds the CSV. The
overlapping physical states must agree exactly; duplicated restart rows are
removed, leaving 1861 states with strictly increasing ages.

With the recorded F77 outputs present, re-extract their accepted models with:

```sh
python3 build_figures.py --f77-work /tmp/ember-f77-core-onset-v1
```

Extraction checks the original file hashes, joins central and surface records
by model number, checks their ages and hydrogen fractions, and verifies row
counts and final ages against the independent F77 summary. It retains the
last accepted Ferguson state before the later convergence failure.

## Writing and daily updates

The paper is updated in place during the day, with one dated paper per day.
Text, tables and figure labels use at most four significant figures. Machine
inputs, calculated data and exact file identifiers retain their required
precision. Write for a reader arriving fresh: explain the calculation and
its meaning, and keep internal workflow terminology out of the discussion.

## Local recovery and publication

The ignored `artifacts/` directory holds deterministic gzip copies of both
completed histories and timing/provenance receipts, the 3.40T native checkpoint,
and the exact macOS arm64 executable. The [manifest](recovery_manifest.json)
records raw and compressed hashes. These files and newly generated bulk physics
tables are not published. A manifest does not provide the numeric inputs it
describes; see [DATA_REPRODUCTION.md](../../DATA_REPRODUCTION.md).

Recover into a new scratch directory and verify all hashes before use.
[Native restart](../../RESTART.md) requires the recorded executable, tables,
mesh, tolerances and physical selections. The selected atmosphere ends at
XH=0.2, close to the current composition; recovering the checkpoint does not
extend its physical support.

September 9 local recovery files remain intact under
`docs/reports/2026-09-09/artifacts/`, including their original manifest. The
earlier paper remains accessible in Git history. Source licenses are unchanged.
