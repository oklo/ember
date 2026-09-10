# Grain absorption and scattering

Ember now has an offline calculation of wavelength-dependent grain absorption,
scattering and the mean scattering direction. **Grain opacity is not yet coupled
to the atmosphere or selected by the evolving star.** The current stellar run
still uses gas-only atmospheres. The new calculation is the first part of adding
suspended grains consistently.

The optional source tools are [LX-MIE](https://github.com/NewStrangeWorlds/LX-MIE)
and [Optool](https://github.com/cdominik/optool). Their versions and archive
checksums are fixed by `scripts/prepare_grain_sources.py`; normal Ember builds
neither download nor link them. LX-MIE and its small adapter are GPL licensed;
Optool is MIT licensed. The original optical-data references are retained in
each generated record. The compilation and spherical-particle calculation are
described by [Kitzmann and Heng (2018)](https://doi.org/10.1093/mnras/stx3141).

For each material and explicitly chosen grain size distribution, the calculation
integrates absorption and scattering cross sections over grain number, then
divides by grain mass. Condensate formula-unit densities from FastChem supply
mass fractions using the same elemental mass convention as the atmosphere.
Mixture opacity is per total gas-plus-grain mass. Absorption and scattering
remain separate: scattered light is not counted as absorbed heat.

The calculation also reports scattering multiplied by one minus its mean
direction cosine. This is a transport diagnostic, not a replacement for the
angular scattering calculation in the atmosphere. Material density, grain
sizes and number weights must be explicit. Missing optical materials and
wavelengths outside the supplied data cause an error; there is no automatic
generic silicate substitution or wavelength extrapolation.

## Checks completed

[Eight independent optical controls](results/grain_mie_v1_audit.json) compare
LX-MIE with Optool, covering small and large particles, weak and strong
absorption, and a nonabsorbing material. The small-particle absorption and
scattering approach the analytic Rayleigh limit within a fractional difference
of 5.633e-7. Tests also check extinction as the sum of absorption and scattering,
mass weighting of a size distribution, zero retained grain mass, a grain with
no optical contrast, and rejection of unsupported input.

Optool imposes an absorption floor in its nonabsorbing-material case. That
case therefore compares its scattering cross section and mean direction only;
Ember's zero absorption is checked against the exact lossless identity.
Optool's angular grid was increased from 180 to 1800 points to resolve its
mean scattering direction. The comparison tolerance was retained.

The [first atmosphere-profile control](results/grain_alumina_profile_v1.json)
uses a 2800 K, log g = 5.15, initially solar-hydrogen gas atmosphere. Its
equilibrium chemistry contains only alumina on that fixed structure, with a
maximum condensed mass fraction of 0.0001094. All of those grains are retained
for this control. The optical data assume amorphous alumina and a material
density of 3.2 g/cm3; the assumptions and original references are in the
[material specification](../data/atmosphere/sources/grain_alumina_optics_specification.json).
This is a sensitivity experiment, not a prediction of grain size or phase.

At a wavelength near 1 micron, grain radii of 0.01, 0.1 and 1 micron give
absorption optical depths of 0.002098, 0.002610 and 0.004681. The corresponding
scattering optical depths are 2.743e-6, 0.002829 and 0.007350. These integrate
only the tabulated atmosphere. Their variation shows why a dust mass fraction
alone does not determine its optical effect. No atmospheric temperature or
interior matching condition has yet been recomputed with these grain opacities.

## Work before an atmosphere can use grains

The existing atmosphere frequency range extends to 0.09 microns, whereas this
alumina compilation begins at 0.2 microns. That gap needs physical data or an
explicitly justified treatment with a measured error. Optical constants also
depend on phase, porosity and temperature. Laboratory alumina changes structure
on heating, so an amorphous optical control cannot by itself establish the
appropriate equilibrium grain phase.

Atmospheric feedback changes which materials condense. The existing 2800 K
calculation with gas depletion already forms calcium titanate and gehlenite
in addition to alumina. The first optical compilation includes calcium
titanate but lacks gehlenite. All relevant materials need coverage or a
quantified treatment before that atmosphere is accepted.

Coupling must include both absorption/emission and scattering, their changes
with temperature and pressure, and a check of angular transport. The atmosphere
must converge again with chemistry and grains responding to the new structure.
Condensate enthalpy remains necessary wherever condensing material carries a
convective heat flux. Existing import checks continue to reject unsupported
cases. Grain growth, mixing and settling require a separate physical model;
fully retained grains and completely settled grains are useful explicit limits.

## Reproduce the optical controls

Use a fresh work directory and a Python environment with NumPy:

```sh
python3 scripts/prepare_grain_sources.py /tmp/ember-grain-sources
python3 scripts/audit_grain_mie.py /tmp/ember-grain-sources/prepared.json /tmp/ember-grain-optical-checks
python3 scripts/generate_grain_opacity.py /tmp/ember-grain-sources/prepared.json data/atmosphere/sources/grain_alumina_optics_specification.json /tmp/ember-grain-alumina
python3 scripts/audit_grain_profile.py /tmp/ember-grain-sources/prepared.json data/atmosphere/sources/grain_alumina_optics_specification.json data/atmosphere/sources/nongrey_gs98_z020/manifest.json docs/results/condensation_nongrey_family.json /tmp/ember-grain-profile --coordinates .7 0 2800 5.15 --radii .01 .1 1
```

The last command needs the separately restored local atmosphere and chemistry
archives. Full optical tables and external source packages remain local under
the [data reproduction policy](DATA_REPRODUCTION.md).
