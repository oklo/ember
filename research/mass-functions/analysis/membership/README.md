# Cluster membership selection test

Completed 2026-10-09. Report and figure:
`docs/research/independent/CLUSTER_MEMBERSHIP_TEST_20261009.*`.
No cluster, injection or download job remains running after completion.

All associated and tidal-tail candidates remain the target. These calculations
diagnose survey selection; they do not infer a corrected birth mass function.
The report distinguishes actual candidate retention from conditional simulation
recovery, whose faint-source Gaia quality losses are still unmodeled.

## Reproduction

From the repository root:

```sh
export PYTHONPATH=out/cluster-stellar-fit-oct9-v1/deps:out/cluster-dissolution-oct9-v1/deps:out/cluster-mass-functions-oct9-v1/deps
export MPLCONFIGDIR=/private/tmp/ember-cluster-mpl
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/source_audit.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/fetch_parent.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/fetch_hyades_extension.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/cp_selection.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/inject_clusters.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/fetch_candidates.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/audit_candidates.py
out/science-env-sept18-v1/bin/python out/cluster-membership-oct9-v1/make_results.py
```

Fetch scripts reuse existing downloads; query text, URLs and SHA256 receipts
are in `raw/`. Broad queries are limited to one million rows and checked against
truncation. An initial 300000-row query was rejected and not used. Final parent
samples contain about 7.455e5 Hyades-field and 3.195e5 Pleiades-field sources. The Hyades
includes both X hemispheres because published members occupy both, despite the
paper's general hemisphere-cut description. The Pleiades uses negative X.

Inputs remain in the three preceding cluster-analysis workspaces. The physical
tracks, mass draws, collisional endpoints and Gaia M10 detection map are reused.
No shared stellar code was edited, and Opus was not redirected from Luhman.

## Outputs and qualifications

- `candidate_selection_audit.csv`: full luminosity-usable union, source overlap,
  temperature/CMD flags and means over mass draws. Never sum overlapping sources.
- `pleiades_group_partition.csv`: nearest published 3D center, not the authors'
  subgroup membership. For comparison with Table 1, add 8122 pc to catalogue xp;
  yp and zp keep the published Galactic frame. The plotted Sun has X=Y=0.
- `*_cp_reproduction.csv`: automatic neighbor/window checks of published members.
  The alternate eroded window is retained separately; it is not the adopted cut.
- `cp_injection_draws.csv`: 288 trials, 2880 count rows, including all escapers.
  Counts obey recovered <= window <= detected <= intrinsic source counts.
- `candidate_recovery_details.csv`: exact Gaia IDs and every cumulative selection
  stage. Corresponding counts are in `candidate_recovery_counts.csv`.
- `conditional_CP_comparison.csv`: compare only with the CP catalogue, not the
  full infrared-plus-Gaia union. Mock quality/manual/SCCP stages remain omitted.
- `temperature_scatter_sensitivity.csv`: assumed 100/200 K Gaussian errors,
  uniform log-mass weighting. This is not a calibrated TIC error model.
- `cluster_membership_test.sqlite`: all compact results and detailed candidate
  and injection tables, with metadata. Source IDs are TEXT, never floating point.
- `validation_receipt.json`: checks, versions and file hashes. The input Gaia
  CSV download hashes are verified independently of derived-table hashes.

The injection field excludes all known candidate Gaia IDs before the simulated
population is inserted. The query data themselves satisfy Gaia quality cuts.
Mock astrometric errors come from that selected field, so they cannot describe
the rejected astrometric-solution population. Mock binary orbital velocities are
replaced by center-of-mass velocities before projected proper motions are formed;
binary photocenter jitter and its effect on RUWE remain unmodeled.

Raw inputs consume roughly 0.3 GB and support exact recovery checks; no new
large dynamical snapshots were generated. Existing source coordinates and model
ages are unchanged. Displayed numbers use at most four significant figures;
machine-readable data retain numerical precision.
