# Maintenance record — November 1, 2026 UTC

The measured project footprint fell from about **502.5 to 375.4 GiB**, a
reduction of **127.1 GiB**. The measurements count allocated file blocks and
shared hardlinks once. They do not independently measure APFS clone extents
or snapshots, so this is not a claim about an identical increase in Finder's
free-space reading. The scientific campaign remains paused.

## Disposition

| Action | Record |
| --- | --- |
| Consolidated 4,740 identical material paths into 1,257 shared objects; estimated 89.34 GiB reduction | `input-consolidation-plan.jsonl`, matching receipts and verification |
| Consolidated 960 older opacity/tabular/archive paths in 197 groups; estimated 6.42 GiB reduction | `archive-consolidation-plan.jsonl`, matching receipts and verification |
| Replaced 615 full interior archives with selected complete profiles; 30.18 GiB reduction | `profile-retention-plan.json`, receipts and per-run `profile_retention.json` |
| Removed 422 generated-cache directories containing 108,073 inventoried files | `build-cache-plan.json`, `build-cache-receipts.jsonl` |
| Removed the 331 MiB temporary SQL index after exporting its records | `temporary-index-removal.json`, `duplicate-checksums.csv.gz` |

The few-hours limit was applied to measured integrations replayable from retained
restarts. Source-table generation is separate. Every reduced profile archive
had a matching executable, configuration and restart identity, a complete
surface history, and seed/final checkpoints. Every selected model retains all
its cells. The full 0.1-solar-mass interiors and the latest paused campaign were
excluded. Exclusions and their reasons remain in `profile-retention-exclusions.json`.

Unique material inputs, uncertain-cost studies, expensive interiors, historical
source snapshots and root Git history remain. `retained-large-files.csv` lists
435 large objects left after cleanup. The Desktop figures occupy only about
25 MiB; the separate historical `low_mass_stars` project about 28 MiB. Both
were left intact. No scientific results were changed or extended.

## Research record

The [study index](../../README.md) contains the exploratory mass spectrum,
continuous 0.1-solar-mass trajectory, boundary comparison, low-metallicity
follow-up and analytic first look. Twenty-two files and archives are pinned
by hashes. All 2,020 referenced provenance paths were present, and every
supplied historical hash matched; the provenance archive preserves their
contents. A separate 792-file source/test archive includes the current
uncommitted code. Existing material tables and historical executable/source
versions remain local and are still needed for full physical replay.

## Code and documentation

- Restored `opacity_cold_dense.hpp` from two agreeing historical snapshots,
  together with the matching optional abundance-basis argument in the table
  reader. The initial checkout failed to compile without these dependencies.
- Replaced 41 repeated CMake test definitions with one small helper. All 83
  test commands and properties are unchanged.
- Shared the repeated warm-opacity joins and expanded dense Henyey boundary
  and residual-check code into readable statements without changing its
  arithmetic order or tolerances.
- Reworked the README around research products and the Henyey FORTRAN IV →
  Bodenheimer FORTRAN 77 → Astra/Opus Ember C++ lineage. Removed stale solver
  status claims and provided a short code guide and reproduction entry point.

The [review criteria and web sources](../../../../docs/DEVELOPMENT.md) distinguish
concrete maintenance problems from unreliable guesses about authorship.

## Validation and audit files

Both the restored baseline and the cleaned code passed **83/83 tests**.
`tests-before.log` and `tests-after.log` retain the results. The 5,700 shared
paths resolve correctly and all retained object hashes were verified. The
compact data and provenance archive passed byte checks, three PNGs were
redrawn and visually checked, and the replay helper printed a command without
launching it. The paused control and restart-memory checksum are unchanged.

The inventories cover 1,314,746 files before and 1,208,902 at the after-scan
checkpoint, excluding root `.git`. Small maintenance files were added and the
temporary SQL index removed after that scan. `summary.json` accounts for the
index removal. `git-status-before.txt`, `git-index-before.txt` and
`tracked-changes-before.patch` preserve the pre-existing dirty worktree record;
this work did not create a commit or rewrite Git history.

The maintenance scripts and plans describe this completed operation. They are
not supervisors and do not run automatically. To make a new inventory:

```sh
python3 scripts/research_inventory.py --output /tmp/ember-inventory
```

The SQL index used for grouping can be rebuilt from `inventory-before/files.csv.gz`;
it was a disposable working file. Receipt files, rather than the index, are the
authoritative record of changes made here.
