# Ember: the long life of a 0.1-solar-mass star

- [Working paper PDF](ember_status_and_future.pdf) and [LaTeX source](ember_status_and_future.tex)
- [Ember and LBA97 in the HR diagram](lba97_hr.pdf)
- [Hydrogen and helium-3 comparison](lba97_composition.pdf)
- [Convective mass and hydrogen profile](convection_current.pdf)
- [Photospheric evolution over LBA97 and Ember opacity](photospheric_evolution.pdf), [PNG](photospheric_evolution.png), [data](photospheric_evolution.csv), and [source checks](photospheric_evolution_provenance.json)
- [LBA97 Figure 6 geometry and approximate opacity calibration](lba97_figure6_extraction.json)
- [Plotted stellar history](evolution_latest.csv), [final structure](evolution_latest_profile.csv), and [input and validation record](evolution_latest_provenance.json)
- [LBA97 figure coordinates and stated values](lba97_figure1_digitization.json)
- [Digitized HR curve](lba97_digitized_hr.csv) and [composition curves](lba97_digitized_composition.csv)
- [Supplementary F77 comparison](f77_matched_hydrogen.pdf), [data](f77_history.csv), and [extraction checks](f77_history_provenance.json)
- [Validation receipt](validation.json) and [local recovery manifest](recovery_manifest.json)
- [Lifetime timeline source](lifetime_timeline.tex), [derived milestones and decay benchmark](lifetime_benchmarks.json), and [conditional decay figure](lifetime_decay_reference.pdf)

The paper compares the 0.1-solar-mass Ember calculation with Laughlin,
Bodenheimer and Adams (1997), abbreviated LBA97. All Ember curves in the main
paper use one 512-point track through **3.875 trillion years**. The sequence
contains the initial stellar structure and 8455 accepted time steps, each
updating both composition and structure. Central hydrogen reaches
**X = 0.001294** and surface hydrogen **X = 0.1730** at the atmosphere table's
**5000 K** limit. The core carries heat by radiation and conduction; conduction
supplies **84.71% of the local central flux**. The paper also compares the
initial model with observed stars near 0.1 solar masses.

The final four pages give the lifetime timeline: formation, hydrogen burning,
helium-remnant cooling, diffusion and crystallization, environmental and nuclear
heating, and conditional nucleon decay through disappearance. Computed ages,
literature comparisons and future conditions are identified separately. The
decay example assumes independent equal lifetimes for all baryons; its final
particle distribution is not a computed Ember cooling trajectory. The text
describes how physical detail is matched to the uncertainty and the question.

LBA97 is represented by 42 selected HR positions, 20 hydrogen positions and
15 helium-3 positions read from Figure 1, together with stated values.
Image coordinates, axis calibration, original PDF checksum and an illustrative
four-pixel reading scale are recorded. These approximate curves are not precise
numerical output or statistically independent measurements. The F77 code seeks
to reproduce LBA97's assumptions and methods; its results form a separate
comparison because its inputs and quantitative evolution differ.

## Rebuild

From this directory, with NumPy, Matplotlib and Tectonic installed:

```sh
python3 build_lba97_figures.py
python3 build_figures.py
python3 build_current_convection.py
python3 build_photospheric_figure.py
python3 build_lifetime_timeline.py
tectonic ember_status_and_future.tex
```

The figures use committed CSV and JSON data. They need neither the original
paper image nor the raw stellar archives. Diamonds distinguish values stated
by LBA97 from points read from its figure. Ages and physical quantities have
not been adjusted to improve agreement. Figure descriptions are in captions;
the plots have no titles.

`build_lifetime_timeline.py` reads the plotted track, derives the helium-3 peak,
first departure from full convection and endpoint, and checks the analytic
decay identities. It writes the timeline values, benchmark record and figure.
Raw data retain precision; generated manuscript values use at most four
significant figures. No new stellar integration is performed by this script.

Thin F77 lines have opacity 0.30; individual segments changing a surface
quantity by more than 5%, or an abundance by more than 0.01, have opacity 0.08.
No data are deleted or smoothed. Abrupt changes are not assumed to be errors.

The photospheric diagram compares LBA97's grain-inclusive map with Ember's
gas absorption map on identical axes and one opacity color scale. LBA97's
156 shaded rectangles and 897 markers for its 0.1-solar-mass track are recovered
from the original vector figure. Its shading has no colorbar: the numerical
scale is approximate, calibrated against Alexander et al. (1983), Table 2.
Nineteen withheld cells differ by at most 0.12 in log opacity. This is a check
on the diagram calibration, not an error bound for the grain physics.
Both panels show the two tracks, with circle diameters proportional to radius.
LBA97 radii are normalized using its stated main-sequence luminosity and
temperature. **Grains and scattering are absent from the Ember background.**
Hatching marks conditions outside its wavelength table.
The photosphere is extracted at Rosseland optical depth 2/3 from all 196
accepted atmosphere structures, using their molecular equation of state.
This is distinct from the interior's matching boundary at optical depth 100.
The exporter verifies original source hashes, all matching states and 65
comparisons with the runtime atmosphere interpolation. It does not establish
an independent photospheric interpolation error bound.

With local source files restored, regenerate those data from the repository root:

```sh
python3 scripts/export_photospheric_evolution.py docs/reports/2026-09-11 /tmp/ember-nongrey-exhaustion-t5000-v1.manifest.json /tmp/ember-nongrey-extended-v2/plane-004/opacity/fort.63 --probe /tmp/ember-nongrey-grid-domain-probe-v2
```

The independent figure extraction uses NumPy, SciPy and PyMuPDF:

```sh
python3 scripts/extract_lba97_photospheric_figure.py /path/to/lba97.pdf /path/to/ajr83.pdf docs/reports/2026-09-11/lba97_figure6_extraction.json
```

With local recovery files present, `python3 build_figures.py --from-archives`
verifies the plotted CSV against the exact archived raw history. The source
result's `converged: false`, source-domain errors and temperature ceiling are
retained in the validation record.

To archive and check a fresh 512-point run, use
`scripts/archive_paper_evolution.py TRACK PAPER_DIRECTORY --mixing-probe PROBE`
from the repository root. It checks the receipt and executable, finite and
normalized profiles, profile/history equality, actual Ledoux mixing regions,
and exhaustion/cooling milestones. An atmosphere-limited result additionally
requires logged domain errors and an endpoint at the recorded temperature limit.

With the recorded F77 outputs present, re-extract their accepted models with:

```sh
python3 build_figures.py --f77-work /tmp/ember-f77-core-onset-v1
```

Extraction checks original file hashes, matches central and surface records
by model number, checks ages and hydrogen fractions, and verifies counts and
final ages against the independently recorded comparison.

## Writing and daily updates

The September 11 paper is overwritten during the day, with one dated paper
retained per day. Prose, tables and figure labels use at most four significant
figures. Machine inputs, calculated data and file identifiers retain their
required precision. Describe the scientific model and its results for a reader
who has not followed the development process. Keep internal development
chronology out of the paper and use measured computing times where relevant.

## Local recovery and publication

The ignored `artifacts/` directory holds deterministic gzip copies of raw
histories, receipts, checkpoints and exact macOS arm64 executables. The
[manifest](recovery_manifest.json) records raw and compressed hashes. Bulk
physics tables and raw outputs remain local under the
[data reproduction policy](../../DATA_REPRODUCTION.md).

Separate resolution-control data through 3.560 trillion years are retained in
`evolution_history.csv`, `evolution_1024_continuation.csv` and the corresponding
profile CSVs. `build_transition_figures.py` rebuilds their comparison. These
controls are not curves in the current paper and do not test resolution at its
3.875-trillion-year endpoint.

Recover into a fresh scratch directory and verify all hashes.
[Native restart](../../RESTART.md) requires identical executable bytes, tables,
mesh, tolerances and physical selections. Changing physical tables requires an
explicitly checked restart procedure or a consistent fresh calculation.
Source licenses are unchanged. The September 10 paper and its recovery records
remain in their dated directory; September 9 is accessible in Git history.
