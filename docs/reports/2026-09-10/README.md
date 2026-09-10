# Ember working paper — September 10, 2026

- [Working paper PDF](ember_status_and_future.pdf)
- [LaTeX source](ember_status_and_future.tex)
- [Completed evolution figure](evolution_history.pdf) and [history CSV](evolution_history.csv)
- [Validation receipt](validation.json)
- [Local recovery manifest](recovery_manifest.json)

This edition replaces the September 9 paper in the current repository. It
describes the completed 512-point, 0.1 Msun calculation through 3.40 trillion
years, refined input families, independent F77 comparisons, performance work,
and the remaining requirements for a helium-remnant cooling track to 100 K.
The star is still fully convective and burning hydrogen.

The 108-node atmosphere candidate, its composition refinement, and the
separate stellar diagnostic toward 3.60 trillion years are ongoing work.
They are not completed-reference results in this paper. For later job status,
read the newest [handoff entry](../../../HANDOFF.md).

## Rebuild

From this directory:

```sh
python3 build_figures.py
tectonic ember_status_and_future.tex
```

The figure uses the committed CSV and requires NumPy and Matplotlib. It can
be rebuilt without raw stellar archives. With local recovery files present,
`python3 build_figures.py --from-archives` checks and joins the fresh
0–3.30T history and exact 3.30–3.40T restart, then rebuilds the CSV. The
overlapping physical states must agree exactly; duplicated restart rows are
removed, leaving 1861 states with strictly increasing ages.

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
