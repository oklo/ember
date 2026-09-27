# Ember working paper — September 27, 2026

[Read the PDF](ember_status_and_future.pdf).

The draft is reduced from 63 to 30 pages while retaining the physical
assumptions, comparisons, numerical qualifications and lifetime timeline.
The continuous Hayashi-started calculation reaches **3.213 Tyr** at **3154 K**,
with **X = 0.2385**, still fully convective. It selects composition-dependent
atmospheres at optical depth 100 throughout. Figure 1 now shows its complete
computed HR trajectory and its initial rise to sustained hydrogen burning.
No continuous flash result is established.

The flash figure and table now use the same accepted histories. The largest
sampled nuclear powers on 666, 848 and 1695 mass points are **5.537e36**,
**4.835e36** and **3.281e35 erg/s**. The peak is not converged, and every flash
continuation retains the changed atmosphere normalization.

The observed-star comparison uses the continuous model near **1 Gyr**. Rates,
EOS composition coverage, transport assumptions and grain limitations are
specified separately for each calculation.

Rebuild with `tectonic --keep-logs ember_status_and_future.tex`.
The section files and the 12 figure PDFs actually included are kept here.
The contraction/main-sequence and late-continuation CSVs and plotting scripts
are adjacent. Rebuild Figure 1 with `python plot_pms_main_sequence.py`; it reads
only the preserved CSV. `continuous_values.tex` records the matching text values. Shared
numerical inputs and their checksums for the other figures remain in
[the retained figure-input directory](../2026-09-26/figure_inputs/manifest.json);
they have not been duplicated. Some comparison scripts also need the retained
local pulse query files identified there; the public PDF alone is not a full
runtime dataset. Figure provenance and publication hashes are in `artifacts.json`.
