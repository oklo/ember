# Generating the local physics inputs

GitHub carries the C++ code, offline generators/importers/auditors, source patches,
generation specifications, small provenance manifests, regression fixtures and
reports. **New bulk EOS/opacity/atmosphere/conduction tables, raw source archives,
native checkpoints and executable archives stay local and are ignored by Git.**
Previously published tables remain in the existing Git history; this policy does
not rewrite that history. Do not force-add newly generated tables or archives.

The development machine retains the completed inputs under `data/`, models under
`out/`, and recovery files under `docs/reports/2026-09-10/artifacts/`. Clearing chat
context does not remove those files. A fresh clone does not contain the new
complete GS98/non-grey families. C++ compilation needs no source-service access,
but stellar runs and table-dependent tests require the corresponding local inputs.
The successful full-suite result in the report used those locally installed inputs.

Source manifests record expected hashes and parameter choices. They are metadata,
not a substitute for the missing numeric source files. Exact reconstruction may
depend on the recorded source version, compiler and original service output.
Changing upstream results or numerical source builds requires a new validation
receipt; do not silently replace a pinned hash to make an import pass.

## Pipeline map

| Input | Generation, retrieval and import code | Independent checks / detailed instructions |
|---|---|---|
| FreeEOS H/He and GS98 potentials | `build_freeeos_probe.py`, `generate_freeeos_grid.py`, `generate_metal_eos.py`, `import_freeeos_potential.py`, `import_metal_eos.py`, `assemble_metal_eos_family.py` | `audit_freeeos.py`, `audit_freeeos_metals.py`, `audit_metal_eos_family.py`; [FREEEOS.md](FREEEOS.md), [EOS README](../data/eos/README.md) |
| TOPS opacity compositions | `fetch_tops_composition.py`, `import_tops_composition.py`, `import_tops_mixtures.py` | `audit_opacity_extension.py`, `audit_tops_heldout.py`; [opacity README](../data/opacity/README.md) |
| AESOPUS low-temperature opacity | `archive_aesopus_mixtures.py`, `import_aesopus.py`, `import_aesopus_mixtures.py` | Original source hashes and unchanged cells; [opacity README](../data/opacity/README.md) |
| Ioffe conduction | `import_conduction.py`; direct-source reference via `conduction_reference_probe.f90`, `generate_conduction_reference.py` | [conduction README](../data/conduction/README.md), [CONDUCTION.md](CONDUCTION.md) |
| Non-grey gas atmospheres | `prepare_nongrey_sources.py`, `generate_nongrey_grid.py`, `archive_nongrey_grid.py`, `import_nongrey_grid.py`, `run_nongrey_plan.py`, `assemble_nongrey_grid.py` | `audit_nongrey_family.py`, source/chemistry/flux checks; [NONGREY.md](NONGREY.md) |
| Condensate experiments | `prepare_fastchem_sources.py`, `prepare_condensate_sources.py`, `generate_condensate_opacity.py`, `run_condensate_atmosphere.py`, `generate_condensate_grid.py`, archive/collect scripts | `audit_condensate_atmosphere.py`, `audit_condensate_material.py`, `audit_condensate_interpolation.py`; [FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md) |

All script names in the table are under `scripts/`. The source distributions keep
their own licenses; see the data READMEs. External Fortran/C++ source builds and
line lists are offline-generation dependencies, not links added to Ember's runtime.

## EOS example

Build the pinned FreeEOS probe as described in [FREEEOS.md](FREEEOS.md). Generate
the original GS98 family in a new work directory:

```sh
python3 scripts/generate_metal_eos.py /tmp/ember-freeeos-source-build/probe \
  /tmp/ember-gs98-source --hydrogen .3 .4 .5 .6 .7 .75 --helium3 0 .12 --jobs 4
```

The selected zero-hydrogen composition extension is the 64-plane family
`data/eos/exhaustion_refined_v2/freeeos300_gs98_z020.dat`. Its explicit source
plane manifests and [528-point audit](results/metal_eos_exhaustion_refined_v2_audit.json)
record the refined composition axis and thermal limits. The earlier 24-plane
family uses `.1 .125 .15 .175 .2 .25 .3 .4 .5 .6 .7 .75` at He3=0/.12.
Generate each identified family in a separate directory. Import and assemble the completed
source planes using the commands/formats in the EOS README and each script's
`--help`. `assemble_metal_eos_family.py` consumes source manifests **and the raw
files they reference**; manifests alone cannot reconstruct the EOS. Recheck
potential masks, source responses, heldout compositions and the thermodynamic
identities before selecting a newly built family.

## Opacity example

The TOPS retriever requests one composition at a time and verifies the returned
mixture. For example, a new independent midpoint request is:

```sh
python3 scripts/fetch_tops_composition.py /tmp/ember-tops-midpoint .15 --metallicity .02
```

Use the committed request/manifests for the complete elemental mixture and source
settings. Recreate all requested compositions before importing a family. Read actual
request fractions: rounded filenames `x012` and `x018` denote `.125` and `.175`.
The source service is an external dependency, so retrieval is not equivalent to
a pinned local archive. The selected hydrogen-poor family is `hydrogen_poor_refined_v4`. Its independent
active-domain source comparisons are in
[the v4 audit](results/opacity_hydrogen_poor_refined_v4_audit.json). Preserve
the rejected coarse midpoint audit when reproducing the refinements.

The current hot-core extension adds X=0/.025/.05/.075 at each Z=.01/.02/.03.
Its complete 54-plane source manifest, 12 independent check records and original
requests are under `data/opacity/sources/hydrogen_exhaustion/`. Raw numeric replies
remain local. TOPS omits the hydrogen row at exactly zero H; the importer verifies
that format and every remaining element rather than substituting a small H value.

```sh
python3 scripts/import_tops_mixtures.py \
  data/opacity/sources/hydrogen_exhaustion/family_manifest.json \
  /tmp/ember-tops-source-import
python3 scripts/assemble_hot_opacity_extension.py \
  data/opacity/hydrogen_poor_refined_v4 /tmp/ember-tops-source-import \
  /tmp/ember-hot-opacity-runtime
```

Only the high-temperature tables are selected from the import. The original cool
TOPS and AESOPUS files are copied unchanged. Every old hot entry is checked for
exact agreement before assembly. This reproduces the selected
`data/opacity/hydrogen_exhaustion_hot_v1` numeric files byte for byte. Low-H
opacity remains unsupported where the cool TOPS branch is required. The
[hot-profile checks](results/tops_exhaustion_hot_profile_all_v2.json) and
[derivative checks](results/hot_opacity_runtime_v1.json) record the tested domain.

AESOPUS is retrieved from the authors' distribution, then archived and imported;
these scripts do not implement the authors' opacity engine. The source URL, hashes
and selected composition planes are documented in the opacity README.

## Non-grey atmosphere example

```sh
python3 scripts/prepare_nongrey_sources.py /tmp/ember-atmosphere-source
python3 scripts/generate_nongrey_grid.py \
  /tmp/ember-atmosphere-source/prepared.json \
  data/atmosphere/sources/nongrey_extended_specification.json \
  /tmp/ember-atmosphere-extended --jobs 3
```

The committed small source assets include the patches, atomic-ion input deck,
element metadata and grid specifications consumed by these scripts. Preparation
downloads and checks the external source packages and line lists. Do not use
`--initial-models` with an absent local archive. From-scratch starting structures
can require additional initialization work; the existing continuation/CONREF
tools help generate initial guesses, and canonical replay supplies acceptance.
This is a reproducible pipeline, not a promise that every cold source cell will
converge without further work.

Archive and import only a complete, independently accepted family following
[NONGREY.md](NONGREY.md). The condensate pipeline is an unfinished experiment:
the grain-enthalpy failures remain real, and no command here promotes it to an
accepted runtime atmosphere. Large generated outputs stay ignored/local.

## What a fresh checkout can reproduce immediately

The report PDF can be rebuilt from its LaTeX source and committed vector figure.
The figures can be regenerated from the committed Ember and F77 history CSV
files and F77 extraction record using
`docs/reports/2026-09-10/build_figures.py`. The F77 file contains the quantities
needed for the graphs, extracted from accepted models; the complete printed
outputs remain local. With local raw histories available,
`--from-archives` repeats the history-join checks and rebuilds that CSV.

For native stellar restart, use the local archived binary/checkpoint and the
exact matching local tables; see [RESTART.md](RESTART.md). Rebuilding tables or
changing the executable does not satisfy the existing checkpoint's identity
contract automatically. On a different machine, generate/validate the inputs
and start a consistent new calculation unless an exact compatible local archive
has separately been provided.
