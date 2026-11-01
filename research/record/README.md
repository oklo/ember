# Ember research record

Scientific campaigns are **paused**. This record was assembled on November 1,
2026 from retained calculations. Replotting it requires no stellar, EOS or
atmosphere calculation.

| Study | Retained result | Qualification |
| --- | --- | --- |
| Exploratory mass spectrum | `exploratory/tracks.json.gz`: 23 tracks, 0.005–0.235 solar masses, including the 0.1-solar-mass evolution into cooling | Cold physics changed along documented continuations. This is not one fixed-physics lifetime calculation. |
| Continuous 0.1-solar-mass study | `exploratory/continuous-0.1Msun.csv.gz` and `continuous-provenance.json`: 16,947 states through 4.005 Tyr | The September 27 snapshot is a distinct, source-limited study. Its recorded cumulative cost is about 31 CPU hours. |
| MS/BD boundary and reheating | `boundary/tracks.json.gz`: the 32 tracks in the October 31 comparison; catalogue, source hashes, ledger and conditional widths | Different mesh sizes, initial states and material approximations must remain distinct. Reported sample brackets are not converged physical widths. |
| Low-metallicity follow-up | `low-metallicity/tracks.json.gz`: five reviewed tracks through the November 1 cutoff | Pop III and Z=0.005 use compact starts. The time coordinate is elapsed integration time, not a calibrated PMS age. |
| Analytic first look | `analytic/map.json.gz`, source manifest and fold review | Extrapolated scenarios and asymptotic theory guide the search; they do not replace the sampled boundary controls. |

The fixed-helium clear/cloudy studies have their own compact data and plots in
[metallicity-boundary](../metallicity-boundary/README.md). The retained
[working paper](../../docs/reports/2026-09-27/README.md) provides the longer
scientific discussion. `paused-state.json` preserves the latest research
qualifications and unfinished controls. Historical status fields inside a
launch record do not indicate a currently running process.

## Verify and redraw

From the repository root:

```sh
python3 research/record/reproduce.py verify
python3 -m venv /tmp/ember-record-env
/tmp/ember-record-env/bin/pip install -r research/record/requirements.txt
/tmp/ember-record-env/bin/python research/record/reproduce.py figures
```

The existing local interpreter at
`out/picaso-pilot-oct11-fable/condaenv/bin/python` also works. The renderer writes
three PNGs under `research/record/figures/`; these are new renderings of the
same saved samples, not pixel-identical recreations of the earlier Desktop
layouts. It does not write to the Desktop. The JSON/CSV data, their original
byte hashes and compressed-file hashes are in `manifest.json`.

`provenance.tar.gz` retains 2,020 source histories, configurations and
review files used by these studies. All 2,020 referenced paths were present
and every supplied historical checksum matched at collection.
`provenance-index.json` lists the archived files. `current-source.tar.gz` preserves the 792-file source/test
snapshot after this maintenance, including uncommitted code. It is separate
from the historical source versions retained under `out/`.

The publication was reconciled with GitHub after that snapshot was made.
The repository source retains the newer upstream physics and tests;
`current-source.tar.gz` records the local cleanup baseline, not a replacement
for the published checkout. See `publication.json` for the publication base.

## Reproduce an integration

Use the executable, configuration, start state and material inputs recorded for
that run. The current checkout is not interchangeable with all historical
executables. For a retained lifetime launch, this command verifies the recorded
executable, configuration and restart hashes and **prints** a command using a
new output directory:

```sh
python3 research/record/reproduce.py command \
  out/adaptive-production-oct18/runs/z2_020_to100000Gyr/launch.json \
  /tmp/ember-replay-z2
```

It does not execute the printed command. Material inputs and their internal
manifests must also be available; see [physics reproduction](../../docs/DATA_REPRODUCTION.md).
Some old workflows require absolute local paths, private input tables or
several source versions. The compact plots are portable; a fresh public clone
is not yet a complete reproduction bundle for all historical integrations.
The obsolete `/tmp/ember-f77-baseline-henyey77.f` reference was already missing
when this cleanup began; the local F77 reconstruction is not evidence that the
missing baseline has been recovered.

## Retention

The user set a recomputation budget of a few hours per study. Unique source
inputs, code, configurations, scientific summaries, complete surface histories,
restart states and expensive or uncertain results remain. The full 0.1-solar-mass
interior archives and the latest paused campaign were excluded from profile
reduction.

Duplicate physics inputs and selected older tabular archives now live once under `data/retained/sha256/`, with
relative symbolic links at their old locations. Keep this directory when
copying the project. These shared inputs are read-only; generate changed inputs
in a new location.

For 615 earlier integrations, full interior dumps were replaced by selected
complete profiles: endpoints, time samples, scalar extrema and radiative-boundary
changes. Each retained launch had matching executable/configuration/restart
hashes and a measured integration cost below three CPU hours. This limit applies
to replay from that restart, not to regenerating an entire campaign or its
material tables. Surface histories and checkpoints remain byte-for-byte intact.
Each affected run has `profile_retention.json`; its old `profile_archive.json`
describes the intentionally removed full dump.

The [maintenance record](maintenance/2026-11-01/README.md) contains the complete
inventory, deletion plans, checksums, receipts and validation results.
