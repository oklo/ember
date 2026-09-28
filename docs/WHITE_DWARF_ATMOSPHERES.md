# White-dwarf atmospheres

The atmosphere supplies pressure and temperature at an optically thick
matching layer. Emergent colors alone cannot provide the boundary needed
for stellar evolution and cooling.

Ember's composition-dependent gas atmospheres match at optical depth 100.
Cooling extensions include mixed hydrogen/helium compositions and molecular
collision-induced absorption. Coverage is defined by each table's axes and
masks. The integrated cooling envelope and some tables remain local research
components; they have not yet been validated in one uninterrupted
Hayashi-to-cooling calculation.
[Configuration](ATMOSPHERE.md) ·
[Boundary checks](results/cooling_boundary_sept28_v1.json).

## Composition and matching

The boundary must follow the calculated surface abundances. Settling can
produce a hydrogen-rich surface; envelope convection can mix hydrogen and
helium again. Metals currently share a common drift velocity, so individual
heavy-element separation is unresolved.

Atmosphere and interior must use compatible thermodynamics at their matching
layer. Comparisons at different depths test the resulting radius and cooling
rate. Cooling also requires the connection between envelope convection and
the conductive interior; this cannot be replaced by an atmosphere color table.

A supplied MESA pure-hydrogen boundary is available for comparisons through
[its importer](../scripts/import_mesa_wd_atmosphere.py). It uses optical depth
25.12, covers 2000–40,000 K and log g = 5.5–9.5, and is restricted in Ember to
helium fraction ≤0.005 and Z ≤1e-20. Ember interpolates it directly without
MESA's transition to a gray boundary below log g = 6. It is not selected for
the continuous Hayashi-origin calculation.

## Remaining physics

Dense, cool atmospheres require supported nonideal ionization and molecular
dissociation, pressure-dependent absorption, the red wing of hydrogen
Lyman-alpha and density-dependent H2 absorption. Refraction and dense-fluid
effects matter where the calculated density requires them. Low gravity makes
uncontrolled extrapolation from ordinary white dwarfs inappropriate.
[Dense-atmosphere calculations](https://arxiv.org/abs/1807.06616) ·
[Collision-induced absorption tests](https://arxiv.org/abs/1809.11122).

Grain chemistry and its coupled radiative effect remain unfinished. Their
importance depends on the metals retained near the surface. The present
wavelength interval ends at 30 microns; a 100 K atmosphere needs far-infrared
coverage and colder chemistry. Molecular data must remain within their
supported temperature and density ranges.
[Grains](GRAINS.md) · [Very cold physics](ULTRACOLD_PHYSICS.md).

## Observational comparisons

Binary evolution produces helium-core white dwarfs within the present age
of the Universe. They provide comparisons for radii, spectra and atmosphere
physics, but their masses and hydrogen layers differ from an isolated
0.1-solar-mass remnant. Their cooling ages are therefore not direct lifetime
tests. [Antoniadis et al. (2013)](https://arxiv.org/abs/1304.6875),
[Calcaferro et al. (2018)](https://arxiv.org/abs/1802.06753).

Published three-dimensional hydrogen atmospheres test convection where
conditions overlap; spectroscopic corrections are not interior boundary
tables. [Tremblay et al. (2015)](https://arxiv.org/abs/1507.01927).
The [Montréal grids](https://www.astro.umontreal.ca/~bergeron/CoolingModels/)
provide further comparisons within their domains, but their carbon/oxygen
cooling sequences do not supply Ember's helium-core evolution.
