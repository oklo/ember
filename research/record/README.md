# Retained research record

This is the November 1, 2026 snapshot of selected exploratory calculations.
It is not a live progress report or a claim that later physics changes have
validated earlier tracks.

| Study | Retained data | Scope |
| --- | --- | --- |
| Exploratory mass spectrum | `exploratory/tracks.json.gz` | 23 tracks, 0.005–0.235 solar masses; physical choices changed across continuations. |
| Continuous 0.1-solar-mass sequence | `exploratory/continuous-0.1Msun.csv.gz` | 16,947 states through 4.005 Tyr; a distinct, source-limited September calculation. |
| Hydrogen-burning boundary | `boundary/tracks.json.gz`, catalogue and ledger | 32 sampled tracks; mass brackets are not converged physical widths. |
| Low metallicity | `low-metallicity/tracks.json.gz` | Five tracks; compact starts do not supply calibrated PMS ages. |
| Analytic first look | `analytic/map.json.gz`, fold review | Conditional extrapolations used to guide further calculations. |

`manifest.json` records original and compressed hashes. `paused-state.json`
records the qualifications at collection, not the state of current campaigns.
The separate [clear/cloudy comparison](../metallicity-boundary/README.md)
retains its own composition conventions and compact data.

## Verify and redraw

```sh
python3 research/record/reproduce.py verify
python3 -m venv /tmp/ember-record-env
/tmp/ember-record-env/bin/pip install -r research/record/requirements.txt
/tmp/ember-record-env/bin/python research/record/reproduce.py figures
```

The renderer produces three PNGs under `research/record/figures/`, without new
stellar or atmosphere calculations. Generated figures are not versioned.

## Reproduce an integration

Keep each run's executable, configuration, starting checkpoint and material
inputs together. For a local retained lifetime launch, the following verifies
the executable, configuration and restart hashes and prints a replay command:

```sh
python3 research/record/reproduce.py command /path/to/launch.json /tmp/ember-replay
```

The current executable is not interchangeable with every historical one.
Some studies need local tables, absolute paths or several source versions;
a public clone is not a complete replay bundle. See [physics inputs](../../docs/DATA_REPRODUCTION.md).

`provenance-index.json` identifies the historical inputs. Full source/history
archives and maintenance inventories stay with the local research workspace.
The hashes and Git revision of the two previously published archives remain
in `manifest.json` under `historical_archives`; verification checks the compact
files distributed here. Unreviewed papers and repeated output are absent from
the current package and remain recoverable from Git history.
