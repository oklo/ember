# Physics inputs

The repository contains the solver, offline generators and importers, source
patches, manifests, compact validation records and test fixtures. Bulk physics
tables, raw source archives, native histories and executable archives remain
local. A fresh clone can build the C++ code but cannot reproduce every stellar
run or table-dependent test without those inputs.

Use the configuration and source manifests for the calculation being reproduced.
They specify compositions, source versions, coverage and hashes. A hash identifies
an input; it cannot reconstruct a missing file. Preserve the original source
output or regenerate and validate it before selecting a new table.

## Source calculations

| Input | Generation, retrieval and import code | Independent checks / detailed instructions |
|---|---|---|
| FreeEOS H/He and GS98 potentials | `build_freeeos_probe.py`, `generate_freeeos_grid.py`, `generate_metal_eos.py`, `import_freeeos_potential.py`, `import_metal_eos.py`, `assemble_metal_eos_family.py` | `audit_freeeos.py`, `audit_freeeos_metals.py`, `audit_metal_eos_family.py`; [FREEEOS.md](FREEEOS.md), [EOS README](../data/eos/README.md) |
| TOPS opacity compositions | `fetch_tops_composition.py`, `import_tops_composition.py`, `import_tops_mixtures.py` | `audit_opacity_extension.py`, `audit_tops_heldout.py`; [opacity README](../data/opacity/README.md) |
| AESOPUS low-temperature opacity | `archive_aesopus_mixtures.py`, `import_aesopus.py`, `import_aesopus_mixtures.py` | Original source hashes and unchanged cells; [opacity README](../data/opacity/README.md) |
| Ioffe conduction | `import_conduction.py`; direct-source reference via `conduction_reference_probe.f90`, `generate_conduction_reference.py` | [conduction README](../data/conduction/README.md), [CONDUCTION.md](CONDUCTION.md) |
| Microscopic collision response | `build_electron_pair_table.py`, `electron_pair_table.py`, `extend_ion_collision_table.py` | `audit_electron_pair_interpolation.py`; [cold coverage and overlap checks](results/electron_pair_eta4096_oct1_v1.json) |
| Non-grey gas atmospheres | `prepare_nongrey_sources.py`, `generate_nongrey_grid.py`, `solve_column.py`, `archive_nongrey_grid.py`, `import_nongrey_grid.py`, `run_nongrey_plan.py`, `assemble_nongrey_grid.py` | `audit_nongrey_family.py`, source/chemistry/flux checks; [NONGREY.md](NONGREY.md) |
| Condensate experiments | `prepare_fastchem_sources.py`, `prepare_condensate_sources.py`, `generate_condensate_opacity.py`, `run_condensate_atmosphere.py`, `generate_condensate_grid.py`, archive/collect scripts | `audit_condensate_atmosphere.py`, `audit_condensate_material.py`, `audit_condensate_interpolation.py`; [FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md) |

All script names in the table are under `scripts/`. The source distributions keep
their own licenses; see the data READMEs. External Fortran/C++ source builds and
line lists are offline-generation dependencies, not links added to Ember's runtime.

## Generate and validate

1. Build the pinned external source, retaining its version, patches and license.
2. Generate the requested compositions and thermal range into a new directory.
3. Import only supported source states; retain failed-state masks and derivative
   stencil exclusions. Check source values, thermodynamic identities and independent
   interpolation points.
4. Compare stellar evolution over a suitable overlap before selecting new inputs.
   Record the executable, configuration, source hashes and numerical differences.

For example, after building the [FreeEOS probe](FREEEOS.md):

```sh
python3 scripts/generate_metal_eos.py /path/to/probe /tmp/ember-eos-source \
  --hydrogen .3 .4 .5 .6 .7 .75 --helium3 0 .12 --jobs 4
```

This is a generation example, not a complete lifetime input family. Follow the
selected manifest for the full set of source compositions and use each script's
`--help` for import and assembly options. Tighter electron integration is available
through `build_freeeos_probe.py --electron-quadrature-error`; source physics and
thermodynamic acceptance checks must remain explicit.

For difficult cold FreeEOS states, `generate_metal_eos.py --start-temperature
200000` first follows temperature at the upper requested density, then samples
each isotherm in decreasing density. This changes the starting guesses, not the
source physics or the import checks; failed target states remain flagged.
`--grid-from` preserves the supplied coordinates, including unequal axis spacing.
Both choices are recorded, and incompatible cached requests are rejected.
[Source and cache checks](results/freeeos_temperature_start_oct1_v1.json).

TOPS retrieval requires the external service or retained original replies.
Non-gray atmosphere generation requires the documented TLUSTY/SYNSPEC sources,
line data and opacity tables. These external source builds are offline dependencies;
Ember does not link to them at runtime.

Electron-pair source integration supports degeneracy parameter eta through
4096. The checked cold table extends its screening coordinate down to 0.005
and retains the classical-ion grid through log strength 5. A matched 100 Myr
stellar comparison changes global quantities by less than 6.893e-13.
This extends coverage of the same screened Born/Pauli model. It does not add
relativistic or strongly correlated mixture physics. Source refinement at the
high-screening corner changes the response by 0.1084%; numerical accuracy is
better in the tested small-screening stellar region. The compact record above
separates source convergence, interpolation, and stellar comparisons.

## Run inputs and archived comparisons

The [lifetime configuration](LIFETIME_DRIVER.md) names all required files.
Keep that configuration with its matching inputs and [checkpoint](RESTART.md).
The integrated cooling envelope and some development inputs are not yet public;
there is no complete downloadable data bundle for the latest calculation.

`docs/results/` retains compact numerical evidence. The working paper includes
its LaTeX and figure PDFs; some underlying figure data remain local. Comparison
CSVs under older report dates are retained because analysis scripts use them.
Superseded papers and figures are available in Git history.
