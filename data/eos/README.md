# Equation-of-state data

**Storage policy:** new generated tables and raw source archives described below
remain local and are not included in a fresh clone. Generator code, source patches,
small specifications and provenance metadata are versioned. See
[the reproduction guide](../../docs/DATA_REPRODUCTION.md). Existing published
data history is retained. References to installed/archived files describe the
development machine unless explicitly stated otherwise.

## GS98 baryonic H/He3 family

The checked stellar trajectory through 3.848 trillion years selects the
72-plane family in `low_density_refined_v2/`. It extends the
low-density range and adds four hydrogen fractions at both helium-3 fractions.
All 5.937 million old source states were reused; 2.220 million additional
low-density states and eight new composition planes were calculated.

The first wider family failed a heat-capacity check. After composition
refinement, 1188 independent comparisons in the added range pass, with maximum
heat-capacity difference 0.2952%. Another 792 comparisons in the old range pass.
All 132 existing atmosphere models now have EOS support, including the four
low-density models outside the older family. These are interpolation and
consistency checks, not physical uncertainty bounds. See the
[new-density checks](../../docs/results/metal_eos_low_density_source_v2.json),
[old-density checks](../../docs/results/metal_eos_low_density_retained_v2.json),
and [atmosphere integration](../../docs/results/nongrey_t3600_low_density_eos_v2.json).

With the local numeric files present, reassemble the checked family without
calculating the source again:

```sh
python3 scripts/assemble_metal_eos_family.py \
  /tmp/ember-eos-low-density-reassembled/freeeos300_gs98_z020.dat \
  data/eos/low_density_refined_v2/sources/low_density_parent_00_manifest.json \
  data/eos/low_density_refined_v2/sources/low_density_parent_01_manifest.json
```

The source specifications and installation record are in that same `sources/`
directory. The older EOS inputs and all running calculations remain intact.

### Numerical electron integration

The accepted 72-mixture family in `numerical_electron_base_v1/` uses direct
numerical electron integrals (FreeEOS options `3 223 -2`). The material source
omits radiation; Ember adds it once. Its source temperature range extends
through 22.39 MK before trimming derivative stencils. All 4736 independent
source comparisons, 2240 thermodynamic identities, 512 stellar-profile zones
and 227 atmosphere queries pass. Maximum heat-capacity differences are
0.2952% across the general tests and 0.09778% on the saved stellar profile.
A fresh stellar calculation selects this EOS with the checked 5000 K atmosphere.
[Acceptance](../../docs/results/numerical_electron_base_acceptance_v1.json),
[general comparisons](../../docs/results/numerical_electron_base_general_v4.json),
[profile comparisons](../../docs/results/numerical_electron_base_profile_v4.json).

The nominal source contains 66 states failing thermodynamic identities and one
failing its density-coordinate criterion. Recomputing their complete isotherms
with tighter electron quadrature retains 9.782 million unaffected states
exactly. One state remains inconsistent; its actual source values are retained
and every potential derivative stencil touching it is excluded. The runtime
rejects queries requiring an unsupported stencil. No acceptance criterion is
relaxed. All preserved source and runtime files match the independently checked
family byte for byte.
[Source repair](../../docs/results/numerical_electron_source_repair_v4.json),
[installation](numerical_electron_base_v1/sources/installation_specification.json).

A separate hot density addition reuses the completed source isotherms. One
exchange iteration fails at the nominal electron quadrature accuracy. Tightening
the relative integration target from 1e-9 to 1e-11 resolves that failure with
unchanged physical formulas and Newton tolerance. The 100 converged controls
agree within 6.382e-10 in scaled physical response; a 605-state pilot also
passes. The generator records the actual executable and input hashes for
every isotherm requiring this fallback. No unsuccessful source values are
filled or extrapolated.
[Precision specification](sources/numerical_electron_precision_fallback_v1_specification.json),
[controls](../../docs/results/numerical_electron_precision_controls_v1.json),
[pilot](../../docs/results/numerical_electron_precision_fallback_pilot_v1.json).

At five dense isotherms, independent quadrature targets of 1e-11 and 1e-13
agree within 2.446e-10 across 605 states. A fresh build reproduces all 605
reference responses exactly. The dense family is therefore being recalculated
at 1e-11, retaining any source isotherm already calculated at that accuracy.
Changes larger than 1e-5 require an independent 1e-13 calculation of the entire
isotherm, agreeing within 1e-8; thermodynamic and coordinate criteria remain
unchanged. This source calculation is not yet an accepted dense EOS.
[Convergence comparison](../../docs/results/numerical_electron_dense_quadrature_convergence_v1.json),
[build reproduction](../../docs/results/numerical_electron_quadrature13_builder_v1.json).

The offline builder requires a fresh working directory and accepts
`--electron-quadrature-error 1e-9`, `1e-11` or `1e-13`.

The density merger preserves the original source rows and the numerical
precision record for each added interval. Unrequested cold dense states are
explicitly absent, and every derivative stencil touching them is masked.
The added domain still requires complete source and interpolation checks.

The following paragraphs describe the original smaller families.

`freeeos300_gs98_z020.dat` selects twelve material-potential planes with
XH=.3/.4/.5/.6/.7/.75 and X3=0/.12, fixed Z=.02. Source number densities
match the atmosphere's GS98 representative-isotope inventory except for
FreeEOS's unsupported trace K (4.56e-6 baryonic mass). Original raw sources
and `sources/freeeos300_gs98_manifest.json` record inputs, transformations,
source failures and hashes. Masked derivative stencils are never used.

Use `MetalHelmholtzEos` with explicit documented-approximation selection and
`MetalInventory::gs98`. Radiation and ideal helium isotope entropy are added
once. The source and interpolation checks, remaining isotope approximations,
and evolution selection are in [FORWARD_EVOLUTION.md](../../docs/FORWARD_EVOLUTION.md).
The H/He proxy families below remain available as controls.

### Hydrogen-poor extension

`freeeos300_gs98_hydrogen_poor_z020.dat` is a separate 24-plane family,
extending XH to .1 with new .1/.125/.15/.175/.2/.25 planes at both X3 values.
It retains the original material grid, source options, masks and all original
table bytes. Existing evolution runs still select the original family.

Coarse .1 composition spacing gave up to 1.22% heat-capacity interpolation
error, prompting the additional planes. The installed family passes 336
fresh source comparisons: maximum pressure, energy and heat-capacity
differences are .04985%, .1087% and .1934%. First-law and response checks
are below 7.2e-10; 1536 archived stellar-profile queries reproduce the old
family byte for byte. These are source/interpolation checks, not physical
uncertainty bounds. Reports are in
`docs/results/metal_eos_hydrogen_poor_{audit,coarse_audit,overlap}.json`.

The raw evaluations and three disjoint parent manifests are archived in
`sources/`. Reassemble into a new directory without recomputing the source:

```
python3 scripts/assemble_metal_eos_family.py /tmp/ember-eos-reassembled/freeeos300_gs98_hydrogen_poor_z020.dat \
  data/eos/sources/freeeos300_gs98_hydrogen_poor_parent_00.json \
  data/eos/sources/freeeos300_gs98_hydrogen_poor_parent_01.json \
  data/eos/sources/freeeos300_gs98_hydrogen_poor_parent_02.json
```

The family SHA-256 is
`93cfb5809dd3548b003fe33aa009f1d572b0303eb58e75fd30c2a961d2cd5924`.
`scripts/audit_metal_eos_family.py` compares a runtime probe built from
`scripts/metal_eos_probe.cpp` against the original FreeEOS source probe;
its audit archives the fresh source queries and executable checksums.

## FreeEOS 3.0 material Helmholtz potential

`freeeos300_hhe_x070_potential.dat` is ember's C2 biquintic representation
of a material free-energy potential from FreeEOS 3.0.0 EOS1. It is the
current fixed-composition equilibrium EOS. The numerical source mixture
is H=.7, He=.3 (metals-as-helium approximation); it is not a full-metal,
He3-capable EOS. Original direct numerical evaluations are retained in
`sources/freeeos300_hhe_x070_raw.json.gz`. Normal builds need no Fortran.

The 305×593 source grid has .0125-dex spacing in log T and
log Q=log[rho/(T/1e6)^1.5]. The potential retains 301×589 nodes after
trimming derivative stencils, with 106 masked nodes. Temperature/density
support is strict. Radiation is removed from the source before forming
the material potential and is added analytically at runtime.

See [FREEEOS.md](../../docs/FREEEOS.md) for equations, source provenance,
reproduction from the checksum-pinned external source archive, and the
small source-fit discontinuities regularized by this representation.
Sampled response discrepancies reach 1–2% near those joins; mathematical
consistency does not establish exact source reproduction or physical accuracy.

```
python3 scripts/import_freeeos_potential.py data/eos/sources/freeeos300_hhe_x070_raw.json.gz /tmp/freeeos-potential.dat
```

SHA-256:

- FreeEOS 3.0.0 source archive: `4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09`
- Raw evaluations: `26165748e1493438f4d62d1e04cb2eb67761a51d57c7abbf99cc0f1a01ff1e62`
- Potential: `34f5e83c898ab74273dd3d907b2ebac77c61976b8e1d7c5e1fda8530c69d823a`

All raw evaluations were reproduced byte for byte from a fresh source
build on the development host; rebuilding the potential is byte-identical.
FreeEOS source is GPL-2.0-or-later and remains outside ember's build; its
source is not relicensed or linked into ember. Cite Alan Irwin and the
FreeEOS release for these numerical calculations.

## Composition family

`freeeos300_hhe_composition.dat` selects the four source H mass fractions
.60, .65, .70 and .75. All use the grid and EOS1 options above, with 43, 89,
106 and 107 masked potential nodes respectively. Original raw responses for
every plane are archived; source and runtime SHA-256 hashes are recorded in
`sources/freeeos300_composition_manifest.json`. All four imports reproduce
byte for byte. The new three raw grids were generated with the same probe
and `generate_freeeos_grid.py --hydrogen X`; every source call converged.

```
python3 scripts/verify_composition_data.py
python3 scripts/import_freeeos_potential.py data/eos/sources/freeeos300_hhe_x065_raw.json.gz /tmp/freeeos-x065.dat --hydrogen .65
```

The runtime composition wrapper preserves H and He isotope number densities
under an explicit density/composition transformation and adds ideal isotope
entropy. Metals remain a He4 proxy. Source support and composition bounds
are strict. See [EVOLUTION.md](../../docs/EVOLUTION.md) for the baryonic
abundance convention, independent off-composition tests and limitations.
The fixed .7 table itself has not changed.

## Retained CMS19 H/He tables

`cms19_h_tp.dat` and `cms19_he_tp.dat` contain the original density and
entropy columns from `TABLE_H_TP_v1` and `TABLE_HE_TP_v1` in the authors'
[DirEOS2019.tar.gz distribution](https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz),
retrieved 2026-09-07. Cite Chabrier, Mazevet & Soubiran (2019), ApJ 872, 51,
[A new equation of state for dense hydrogen-helium mixtures](https://arxiv.org/abs/1902.01852).
The archive root is `DirTABLES-EOS2019/`. Its README identifies the pure
component files as February 2019 v1; the June 2021 archive update added a
mixture table. These external scientific data are not relicensed under
ember's MIT code licence.

The authors' June 2021 interacting-mixture README recommends the 2019
additive-volume tables for low-mass stars. We use pure H and He, never the
2021 effective-hydrogen table as pure H. The optional metals-as-helium
approximation must be explicitly selected in `Cms19Eos`.

### Reproduction

```
curl -fLO https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz
python3 scripts/import_cms19.py DirEOS2019.tar.gz data/eos
python3 scripts/audit_cms19_energy.py DirEOS2019.tar.gz
```

The importer and audit verify SHA-256
`736de2a0b26c02b897fbd906504bb6aa144f08a6d0f3cff36da62970054d149b`.
Reimport was checked byte for byte against both versioned files. Python's
standard library is sufficient; tables are read directly by C++ at runtime.

### Format and support

Each file has 121 temperatures by 441 pressures, or 53,361 pairs of original
values. The format is:

```
121 441 CMS19 H: <description>     # HE: for helium
<121 log10 T[K] values>
<441 log10 P[GPa] values>
<log10 rho[g/cm3]> <log10 S[MJ/kg/K]>  # pressure varies fastest
...
```

Axes run from log T=2 to 8 and log P[GPa]=-9 to 13, both in steps of .05.
The importer retains the original printed density/entropy strings, including
unphysical rectangular corners. It neither interpolates nor repairs them.
The runtime supports a much smaller fluid subset: log T=3.2..7.3, with
each entire 4×4 interpolation stencil inside density and phase masks.
Consequently an imported rectangle does **not** imply rectangular physical
coverage. See [CMS19 implementation and audit](../../docs/CMS19.md).

### Static calculations only

Energy columns are deliberately absent. At rho=1 g/cm³, the original He
table's internal energy decreases between log T=5.95 and 6.0 despite a
positive entropy-derived heat capacity. The source H table also violates
the pressure/entropy Maxwell identity by 26.5% at one audited dense
ionization point. These are reproducible diagnostics of this particular
source version, not a claim that every CMS19 state has such errors.

`Cms19Eos::has_internal_energy()` is false. Energy and its density
derivative are NaN; time-dependent structure and central boundary calls
reject this EOS. Density, entropy and their derivatives support an
**experimental static calculation**, with consistency errors exposed in
the equilibrium JSON. Do not manufacture internal energy or enforce one
identity by overwriting an independently derived response.

## Extended evolution family

`freeeos300_hhe_extended.dat` adds source X=.3/.4/.5/.55 to the original .6/.65/.7/.75 family, using the same FreeEOS 3.0 EOS1 options and .0125-dex potential grid. All eight raw source planes are archived; `sources/freeeos300_extended_manifest.json` pins both raw and imported hashes. The old manifest is unchanged for historical reproduction.

```
python3 scripts/verify_composition_data.py --extended
python3 scripts/generate_freeeos_composition_reference.py /tmp/freeeos-probe \
  tests/data/freeeos300_extended_reference.dat --extended
```

The 63 direct off-composition queries, including He3 mass fractions up to .105, differ by at most 0.083% in P/E and 0.421% in thermal responses. These test the stated elemental/isotope approximation; they are not comparisons with an independent physical EOS. The metals-as-He4 approximation remains explicit. A separate [metal sensitivity audit](../../docs/results/extended_eos_metal_sensitivity.json) compares direct FreeEOS evaluations with GS98 metals at 15 sampled states; differences are below .91% in P/E and .37% in cp/adiabatic gradient. That audit omits unsupported potassium and does not provide a consistent new metal mixture grid.

Pressure inversion now honors the composition family's density support even when the initial ideal-gas guess lies outside it. Immutable mask intervals and integer polynomial scale factors are precomputed for speed; the potential, interpolation order, masks and physical derivatives are unchanged.
