# Ember working paper — September 10, 2026

- [Working paper PDF](ember_status_and_future.pdf)
- [LaTeX source](ember_status_and_future.tex)
- [Ember and LBA97 in the HR diagram](lba97_hr.pdf)
- [Hydrogen and helium-3 comparison](lba97_composition.pdf)
- [Ember history](evolution_history.csv), [1024-point continuation](evolution_1024_continuation.csv), and [comparison checks](../../results/evolution_transition_3560gyr_v1.json)
- [Convection transition and hydrogen profiles](convection_transition.pdf)
- [LBA97 figure coordinates and stated values](lba97_figure1_digitization.json)
- [Digitized HR curve](lba97_digitized_hr.csv) and [composition curves](lba97_digitized_composition.csv)
- [Supplementary F77 comparison](f77_matched_hydrogen.pdf), [data](f77_history.csv), and [extraction checks](f77_history_provenance.json)
- [Validation receipt](validation.json)
- [Local recovery manifest](recovery_manifest.json)

The paper compares the ongoing 0.1-solar-mass Ember calculation with the
published model of Laughlin, Bodenheimer and Adams (1997), abbreviated LBA97.
Both calculations now reach 3.560 trillion years. The complete 512-point
history contains 2354 states; the 1024-point continuation contains 532 states
from its saved restart near 3.545 trillion years. Both now have stable central
regions comprising about 37–38% of the mass, with convective envelopes outside.
The final luminosities differ by 0.07813%, while central hydrogen differs by
1.752%. Further mesh and time-step checks are needed. Hydrogen burning continues;
no white-dwarf cooling endpoint has been reached.

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
python3 build_transition_figures.py
tectonic ember_status_and_future.tex
```

The scripts require NumPy and Matplotlib. The LBA97 figures use the committed
512-point CSV, 1024-point continuation CSV and figure-coordinate JSON; they need neither the
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
reads the complete new 512-point calculation, rather than joining histories
with different atmosphere inputs. It must contain 2354 states with strictly
increasing ages from zero through 3.560 trillion years. The earlier 3.40-trillion-year
history and the saved 3.548-trillion-year trial remain in local recovery files
and Git history. The old trial-state JSON is retained as evidence for the
previously convective center; it is no longer used as a plotted endpoint.

From the repository root, `scripts/summarize_transition_runs.py --mixing-probe
/tmp/ember-evolution-mixing-probe-v1` rechecks the raw completed calculations,
receipts, profile compositions and actual Ledoux mixing regions, then exports
the published histories, profiles and comparison report. Its report records
the probe and input checksums. The transition figure uses only the published
history and profile CSV files.
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

The ignored `artifacts/` directory holds deterministic gzip copies of the
completed histories and receipts, saved checkpoints and exact macOS arm64
executables. The [manifest](recovery_manifest.json) records raw and compressed
hashes for 14 files, including the older calculation. These files and generated
bulk physics tables remain local; see [DATA_REPRODUCTION.md](../../DATA_REPRODUCTION.md).

Recover into a new scratch directory and verify all hashes before use.
[Native restart](../../RESTART.md) requires identical executable bytes, tables,
mesh, tolerances and physical selections. The new 512-point checkpoint supports
further continuation. The 1024-point continuation used its unchanged earlier
binary without writing another checkpoint; its older saved state is preserved.
Atmosphere composition coverage now reaches hydrogen fraction 0.1 in the
specified warm models, with incomplete coverage of other combinations.
September 9 local recovery files remain intact under
`docs/reports/2026-09-09/artifacts/`, including their original manifest. The
earlier paper remains accessible in Git history. Source licenses are unchanged.
