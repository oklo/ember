# Ember scientific and operational checkpoint — 2026-09-09

- [Report PDF](ember_status_and_future.pdf)
- [LaTeX source](ember_status_and_future.tex)
- [Next-session prompt](../../NEXT_SESSION_PROMPT.md)
- [Local recovery artifact manifest](recovery_manifest.json)
- [Validation receipt](validation.json)

This snapshot describes the completed 2.85-trillion-year gas-atmosphere model,
the physics and numerical improvements, the unfinished source grids, and the
requirements for complete lifetimes and helium-remnant cooling. It preserves
the distinction between accepted inputs, failed/diagnostic experiments and
future work. The physics snapshot is 23:09 UTC; validation was performed afterward.

**Publication scope:** code, generation recipes, small metadata and this report
are pushed to GitHub. Bulk tables and everything in `artifacts/` stay local and
ignored. The hash manifest describes those local recovery files; it does not mean
their contents are bundled in a fresh clone. See
[DATA_REPRODUCTION.md](../../DATA_REPRODUCTION.md) for the table-generation pipeline.

## Rebuild the document

From this directory, with Tectonic installed:

```sh
tectonic ember_status_and_future.tex
```

The vector figure and compact history CSV are committed. Regenerate the figure
from the CSV using Python with NumPy and Matplotlib:

```sh
python3 build_figures.py
```

If the local raw histories are available, `python3 build_figures.py --from-archives`
also checks their joins and rebuilds the CSV. Raw histories are not published.

The figure combines the original 0->2.5T checkpoint run, 2.5->2.75T restart,
and 2.75->2.85T restart. It checks matching states at both seams and removes
duplicated initial restart rows. It is one supported-physics trajectory with
stops/restarts, not a claim that a single uninterrupted run has identical steps.

## Recovering the native checkpoint

On the development machine, the ignored `artifacts/` directory contains deterministic gzip archives of all three raw evolution
JSON segments, their SHA-256 provenance receipts, the final native checkpoint,
and the exact original macOS arm64 `ember-evolve` executable. Raw and compressed
hashes appear in `recovery_manifest.json`. The executable is a local, platform-specific
research artifact; it is not published. Rebuild the source for a different platform.

Decompress chosen artifacts into a **new** scratch directory, verify their
uncompressed SHA-256 hashes against the manifest, and give the recovered binary
execute permission there. Preserve the selected data files exactly. Follow
[RESTART.md](../../RESTART.md) for native restart semantics. No recovery helper
silently overwrites an existing model or changes table selections.

The present atmosphere stops at XH=.3, so the 2.85T state has little forward
coverage. Adding different physics requires a new calculation under the current
checkpoint contract; merely recovering this checkpoint does not authorize or
enable a changed-physics continuation. A fresh build may have different bytes
even if the source is unchanged.

## Pending work preserved here

The local `artifacts/pending_tops_refinement/` directory contains four completed raw TOPS source
requests plus the source controller's completion receipt. They have **not** been
imported or independently accepted as a refined runtime family. The filenames
with `x012` and `x018` encode actual X=.125 and .175. Read the requests.

Controller log snapshots establish completion/failure at inspection time, not
the current liveness of any background process. The final X3=0 solar 2800 K,
logg4.9 atmosphere failed the grain-enthalpy audit after converging. No condensate
install receipt existed, and the stellar reference includes no condensates.

The data and source artifacts elsewhere in the repository remain subject to
their documented source licenses. This report does not change those terms.
