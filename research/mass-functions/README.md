# Low-mass populations: Hyades, Pleiades and the field

Ember tracks convert observed luminosities into provisional mass distributions.
The cluster sample includes identified associated populations, escaped members
and tidal tails. Collisional simulations and catalogue-selection tests examine
why the clusters appear deficient in brown dwarfs relative to the nearby field.

**A difference in the birth mass functions has not been established.** Faint
surveys cover much less sky than stellar surveys, Gaia quality cuts remove many
brown-dwarf candidates, and the extended Pleiades catalogue selects warmer stars.
The cluster models also fail to reproduce the full stellar spatial distribution.

## Results

| Analysis | Report | Figure | Numerical results |
|---|---|---|---|
| Cluster candidate mass distributions | [Report](reports/cluster-mass-functions.md) | [Normalized comparison](figures/normalized-cluster-mass-functions.png) | [Tables](tables/clusters) |
| Volume-limited field population | [Report](reports/field-mass-function.md) | [Three populations](figures/field-mass-function.png) | [Tables](tables/field) |
| Dissolution and survey coverage | [Report](reports/cluster-dynamics.md) | [17 cluster models](figures/cluster-dynamics.png) | [Tables](tables/ensemble) |
| Fit stars, then predict brown dwarfs | [Report](reports/stellar-fit.md) | [Spatial comparison](figures/stellar-fit.png) | [Tables](tables/stellar-fit) |
| Membership and faint-source selection | [Report](reports/membership-selection.md) | [Selection tests](figures/membership-selection.png) | [Tables](tables/membership) |

Every figure also has a PDF in [figures](figures). The reports give source
citations, model assumptions and limitations. Results are current through
October 9, 2026; the field census is the archived 2024 catalogue, not a live
inventory of discoveries.

For the fixed comparison intervals, define
**R = N(0.03–0.08 M☉) / N(0.10–0.20 M☉)**. These intervals are not an exact
classification of hydrogen burning.

| Sample | Brown-dwarf candidates | Reference-star candidates | R |
|---|---:|---:|---:|
| Hyades and associated populations | 47.02 | 361.4 | 0.1301 |
| Pleiades and associated populations | 201.0 | 594.8 | 0.3379 |

These are means over mass-inference draws, conditional on candidate membership.
They are source counts, with unresolved binaries blended, and have no complete
contamination or completeness correction. Some report tables use medians and
therefore differ slightly. The field **component** ratio has a model-sensitivity
range of about **0.55–0.71**; comparison with source counts requires a binary model.

The 17 cluster simulations reach 650 Myr for the Hyades and 125 Myr for the
Pleiades, retaining escapers. Their birth ratios are **assumed inputs**, not
measurements of the initial population. The simulations used **8.931 CPU-hours**;
this excludes the Ember tracks, installation and catalogue analysis.

The most recent selection test retrieves Gaia attributes for 5415 candidate
counterparts and performs 288 simulated-member injections. Quality cuts retain
26.75% and 47.53% of the Hyades and Pleiades brown-dwarf counterparts, respectively.
These selected-sample retention rates are **not independent completeness estimates**.
The injection predictions still omit those quality losses.

The next useful inference needs independently calibrated faint-source and
infrared selection, plus a description of the associated groups. More cluster
realizations alone would not remove these limitations.

## Files and reuse

- [results.sqlite](results.sqlite): queryable copy of the published summary tables;
  table descriptions and file paths are in its `datasets` table.
- [Reproduction guide](REPRODUCING.md): dependencies, retained inputs and command order.
- [analysis](analysis): analysis and generation scripts, arranged by study.
- [provenance](provenance): input lists, assumptions, checksums and validation records.
- `local/`: ignored links to the seven working datasets on the research machine.

To verify this bundle, run `python3 research/mass-functions/verify.py` from the
repository root. Raw catalogues, full draw tables, stellar histories and cluster
snapshots remain local; this is not a self-contained release of every input.

```sql
SELECT table_name, row_count, description FROM datasets;
SELECT * FROM membership__candidate_recovery_counts;
SELECT * FROM ensemble__execution_cost;
```
