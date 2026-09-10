# AMES-COND non-grey tau=100 boundary

**Storage policy:** new generated tables and raw source archives described below
remain local and are not included in a fresh clone. Generator code, source patches,
small specifications and provenance metadata are versioned. See
[the reproduction guide](../../docs/DATA_REPRODUCTION.md). Existing published
data history is retained. References to installed/archived files describe the
development machine unless explicitly stated otherwise.

The source `tau100_T.data` and `tau100_Pgas.data` files are from
[MESA](https://github.com/MESAHub/mesa/tree/fd396fd73d3f936da8063ffdf9d92361882eb557/atm/atm_data),
commit `fd396fd73d3f936da8063ffdf9d92361882eb557`, retrieved 2026-09-07.
They supply temperatures in K and material pressures in dyn/cm² at
Rosseland optical depth 100. Their header identifies solar GN93, [M/H]=0.
The selected cool-dwarf region comes from AMES-COND-2000 structures
(Allard et al. 2001), as documented by
[MESA's atmosphere documentation](https://docs.mesastar.org/en/latest/atm/table.html).
These scientific data are not relicensed under ember's MIT code license.

The full MESA rectangle includes grey integrations, extrapolation and
transition smoothing. The importer **only selects Teff=1800..3300 K in
100-K steps and log g=3.5..6.0 in .5-dex steps** (96 nodes). Inspection of
[`create_tau100.f90`](https://github.com/MESAHub/mesa/blob/fd396fd73d3f936da8063ffdf9d92361882eb557/atm/preprocessor/src/create_tau100.f90)
establishes that this region retains COND values: grey integration is at
log g≤2, C&K overwrite begins at 3500 K, and transition smoothing alters
3400/3500 K. None of those cells is imported. The converter adds PHOENIX's
electron pressure to its heavy-particle pressure, so the distributed
`Pgas` is material pressure; ember adds radiation once.

The source's detailed mixture and helium abundance are not matched to
ember's X=.7,Z=.02 interior. The output is explicitly a **solar-mixture
proxy**, never labeled an exact-composition atmosphere. Version 2 records
both that approximation and the allowed interior composition. The reader
rejects it unless constructed with `Mixture::allow_documented_proxy`;
the CLI selects this with `--atmosphere cond-solar-proxy`.

```
python3 scripts/import_cond_atmosphere.py data/atmosphere/sources \
  data/atmosphere/cond_gn93_tau100_solar_proxy.dat
```

Pinned source SHA-256:

- Pgas: `77606d68e2f0c02602624e04e31792c5150c7a15940b53cc8a364ba3682f5900`
- T: `e9551fac9e2d8d036ab94c55ed87cb66811d2a0c54e82d27c846c1aabdea4c3a`
- MESA converter inspected: `ab5f7d105bf10983d674b998b641ecaa6319f0a32bc7484a42d5b83e490ae57f`

No spectrum file was used to infer pressure. These are extracted boundary
states from model structures, not a complete atmosphere transfer solver.
The interface retains the thin-atmosphere approximation: the tau=100
radius is identified with the photospheric radius and the omitted mass is
neglected. Resolving that geometric extension and matching modern mixture
and atmosphere physics remain further improvements.

## Historical bounded initial evolution

The underlying table and fixed-composition reader remain strict, including
abundance basis. `ember-evolve` with transport option `early` selects
`FrozenCompositionAtmosphere`: the same GN93 T/P relation with density
recomputed from the actual evolving EOS. It rejects |Xsurface-.7|>.005,
He3>.005, changed metals or unsupported Teff/gravity. These operational
limits do not establish an atmosphere-error bound. A composition-dependent
non-grey source grid is still needed for substantial depletion; see
[EVOLUTION.md](../../docs/EVOLUTION.md).

## Extended evolution

The default evolution boundary now multiplies these original solar COND
T and gas-pressure states by composition ratios from an independently
integrated grey convective column. The column uses the evolving EOS and
elemental opacity, including helium isotope number-density changes. The
reference composition recovers the original COND boundary exactly.

This is a differential approximation, not a helium-enriched non-grey source
grid. It retains the source Teff/gravity bounds and the column's EOS/opacity
bounds. The historical frozen wrapper and its abundance caps remain intact.
Equations, derivative checks and atmosphere sensitivity comparisons are in
[ATMOSPHERE.md](../../docs/ATMOSPHERE.md); the longer track is documented in
[EXTENDED_EVOLUTION.md](../../docs/EXTENDED_EVOLUTION.md).

## Helium-rich non-grey source pipeline

The optional pinned TLUSTY208/SYNSPEC54 pipeline and
`CompositionAtmosphereGrid` backend are described in
[NONGREY.md](../../docs/NONGREY.md). They compute composition-dependent
radiative/convective gas atmospheres from atomic and molecular opacities.
The installed `nongrey_gs98_z020_tau100.dat` contains 48 accepted physical
atmospheres: XH=.45/.7, X3=0/.12, Teff=2600/2800/3000/3200 K and
log g=4.9/5.15/5.4, at fixed Z=.02 and Rosseland matching depth 100.
Every model uses 300 depths and 20,000 transfer frequencies. Runtime
interpolation supplies temperature and gas pressure; Ember computes density
from its own EOS at the evolving composition.

Original atmosphere inputs, outputs, starting structures and source receipts
are archived in `sources/nongrey_gs98_z020/`, with a manifest that reproduces
the installed grid offline. Large opacity binaries are represented by hashes
and their pinned generation recipe. Source patches, line-data metadata and
the specification are in `sources/`; independent numerical controls are in
`sources/nongrey_validation/`. All cells pass the flux, hydrostatic,
chemical-density, convergence and material-domain checks. The family still
omits condensates/rainout and retains the documented metal/isotope/EOS
approximations. It does not cover hydrogen exhaustion.

## Forward extension and condensation experiment

The XH=.3 extension and shared GS98 interior mixture are documented in
[FORWARD_EVOLUTION.md](../../docs/FORWARD_EVOLUTION.md). The 72-cell gas
family is installed separately as `nongrey_gs98_z020_extended_tau100.dat`.
All 72 cells pass original-source checks and offline reimport. Its original
artifacts are in `sources/nongrey_gs98_z020_extended/`; the older 48-cell
table remains available. The extension reaches XH=.3 at fixed Z=.02,
with the same temperature, gravity and helium-3 axes.

Pinned FastChem 4 chemistry is an optional offline dependency. The source
preparation scripts build a separate TLUSTY/SYNSPEC experiment that removes
equilibrium condensates from the gas and sets grain opacity to zero. The
hybrid chemistry retains the original gas partition functions and line data.
It does not include grain enthalpy, grain-size evolution or a cloud model.
Ember's runtime library does not link FastChem. Source licensing and chemistry
receipts are retained in `sources/condensation_gs98/`.

Completed opacity experiments are archived under
`sources/condensation_opacity_x700/` (original 16-temperature solar plane),
`sources/condensation_opacity_x700_extended/` (17 temperatures, 1001–7762 K),
and `sources/condensation_opacity_x300_he3120/` (16-temperature helium-rich
plane). Each archive contains the original frequency tables, input decks,
source logs and checksums. They are opacity planes, not accepted atmosphere
families. Coupled atmosphere validation remains in progress.
The matching 17-temperature gas control is archived in
`sources/condensation_gas_opacity_x700_extended/`.

Three matched-control models now pass source and independent chemistry checks:
`sources/condensation_atmosphere_x700_2600_g515/` and
`sources/condensation_atmosphere_x300_he3120_2800_g515/`, plus
`sources/condensation_atmosphere_x700_2800_g515/`.
Their pressure changes relative to gas controls on the same material grids
are +.8095%, +.4833% and +.12461%, respectively. They place condensation entirely in
layers with zero convective heat flux. These are individual source models;
a complete composition-dependent boundary family is still being computed.
The attempted cold rectangle is not yet supported: XH=.3, Teff=2600 K,
logg=5.15 trials cross the atomic partition-function floor of 1000 K and are
rejected. Warmer continuation nodes are being tested; a production rectangle
must bracket the actual stellar track and pass every original source check.

`prepare_fastchem_sources.py` and `prepare_condensate_sources.py` build the
offline sources; `generate_condensate_opacity.py` computes a plane;
`archive_condensate_opacity.py` verifies and archives its original artifacts.
`run_condensate_atmosphere.py` then solves one atmosphere, requiring the
same strict final diagnostics as the gas family. A completed solution must
also pass `audit_condensate_atmosphere.py`, which independently checks
chemical closure and excludes condensation in heat-carrying convective
layers because grain enthalpy is omitted.
`generate_condensate_grid.py` runs temperature-continuation chains, and
`archive_condensate_grid.py` requires a complete rectangular family before
importing it. The importer rechecks archived inputs, source fingerprints,
elemental numbers, original chemistry outputs and the enthalpy restriction.
The offline checks in `tests/test_condensate_import.py` require numpy.

The optional `prepare_atmosphere_initializer.py` provides Newton-step,
opacity-derivative and native CONREF initialization controls. Numerical
changes affect starting models only. Those outputs are marked as
initializers and require a separate original-source replay. No trial
atmosphere, extrapolated material state or unconverged control supplies an
Ember boundary condition.
