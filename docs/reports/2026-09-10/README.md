# Ember working paper — September 10, 2026

- [Working paper PDF](ember_status_and_future.pdf)
- [LaTeX source](ember_status_and_future.tex)
- [Ember and LBA97 in the HR diagram](lba97_hr.pdf)
- [Hydrogen and helium-3 comparison](lba97_composition.pdf)
- [Ember history](evolution_history.csv) and [later trial state](ember_transition_trial.json)
- [LBA97 figure coordinates and stated values](lba97_figure1_digitization.json)
- [Digitized HR curve](lba97_digitized_hr.csv) and [composition curves](lba97_digitized_composition.csv)
- [Supplementary F77 comparison](f77_matched_hydrogen.pdf), [data](f77_history.csv), and [extraction checks](f77_history_provenance.json)
- [Validation receipt](validation.json)
- [Local recovery manifest](recovery_manifest.json)

The paper compares the ongoing 0.1-solar-mass Ember calculation with the
published model of Laughlin, Bodenheimer and Adams (1997), abbreviated LBA97.
The completed Ember track ends at 3.40 trillion years. Open orange points show
the later saved state of a refined-atmosphere trial at 3.548 trillion years.
That trial develops a stable shell around a convective center and stopped
because of a checkpoint counter limit. Corrected and higher-resolution
calculations are checking the transition. No cooling endpoint has been reached.

LBA97 is represented by 42 selected HR positions, 20 hydrogen positions and
15 helium-3 positions read from Figure 1, plus independently stated values.
The original numerical history has not been recovered. Image coordinates,
axis calibration, original PDF checksum and an illustrative four-pixel reading
scale are recorded explicitly. These approximate curves are not precise
numerical output or statistically independent measurements. The current F77
reconstruction is a separate comparison, not the original LBA97 program.

## Rebuild

From this directory:

```sh
python3 build_lba97_figures.py
python3 build_figures.py
tectonic ember_status_and_future.tex
```

The scripts require NumPy and Matplotlib. The LBA97 figures use the committed
Ember CSV, trial-state JSON and figure-coordinate JSON; they need neither the
original paper image nor raw stellar archives. Diamonds distinguish values
stated by LBA97 from points read from its figure. The stellar ages and physical
quantities have not been adjusted to improve agreement.

The supplementary F77 figures use the committed CSV files and extraction
record. Luminosity and radius use the same nominal solar units as Ember.
Thin F77 lines have opacity 0.30; individual segments changing a surface
quantity by more than 5%, or an abundance by more than 0.01, have opacity 0.08.
No data are deleted or smoothed. No claim is made that every abrupt change is
numerical error.

With local recovery files present, `python3 build_figures.py --from-archives`
checks and joins the fresh 0–3.30T history and exact 3.30–3.40T restart.
Overlapping physical states must agree exactly before the duplicate restart
row is removed, leaving 1861 states with strictly increasing ages.
The later trial point can be re-extracted from its unmodified checkpoint using
`scripts/extract_transition_point.py`; this conversion is only for plotting
and does not create a restart file.

With the recorded F77 outputs present, re-extract their accepted models with:

```sh
python3 build_figures.py --f77-work /tmp/ember-f77-core-onset-v1
```

Extraction checks original file hashes, matches central and surface records
by model number, checks ages and hydrogen fractions, and verifies counts and
final ages against the independent F77 summary. The last accepted Ferguson
state is retained before its later convergence failure.

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
