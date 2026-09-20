# Atmospheres through contraction and white-dwarf cooling

Ember evolves the structure, composition and thermal energy of an initially
0.1-solar-mass star. Its atmosphere supplies pressure and temperature at a
deep matching layer. Cooling requires this physical boundary as well as an
emergent spectrum; a table of colors alone cannot supply it.

## Conditions to calculate in advance

The checked model reaches 3.881 trillion years, effective temperature 5200 K,
surface hydrogen X = 0.1730 and central hydrogen X = 0.001084. It retains
0.006867 solar masses of hydrogen. These are different measures: central
exhaustion does not determine the composition of the atmosphere or remove
the remaining nuclear fuel.
[Checked endpoint](results/evolution_atmosphere_limit_3881gyr_v1.json).

LBA97 provides a useful estimate of the region that future atmospheres must
cover. The values below are approximate readings from our
[digitized HR curve](reports/2026-09-11/lba97_digitized_hr.csv), not tabulated
stellar structures. We infer radius from luminosity and effective temperature,
then gravity from a mass of 0.1 solar masses. They guide source calculations;
Ember's trajectory remains an independent result.
[LBA97, Figure 1](https://www.astroexplorer.org/details/10_1086_304125_fg1).

| Position on the inferred LBA97 curve | Effective temperature (K) | log g (cgs) |
| --- | ---: | ---: |
| Near the highest effective temperature | 5800 | 5.77 |
| Contracting, cooler surface | 4700 | 6.19 |
| Cooling | 3540 | 6.38 |
| Cooler remnant | 2590 | 6.49 |
| Cool end of the digitized curve | 1760 | 6.54 |

The immediate calculation plan extends the temperature grid to 6400 K and
starts higher-gravity trials before the star requires them. The warm grid
brackets its measured surface composition with X = 0.15 and X = 0.2, and
helium-3 mass fractions 0 and 0.12. Complete interpolation cells and independent
temperature and lower-boundary comparisons are needed before a table is used.
Each calculation reuses compatible wavelength opacities and a converged
starting structure. Independent source models and their controls run in parallel.
The 220-source table through 5600 K has passed its source, interpolation,
lower-boundary, condensation and interior-EOS checks. Its two independent
matching-state comparisons differ by at most 0.4169%. The stellar continuation
using it reaches 3.890 trillion years and 5546 K before an interior opacity
density limit. Its checkpoint and input checks pass; the full endpoint audit
is pending. [Continuation checks](results/evolution_density_limit_3890gyr_initial_v1.json).
[5600 K atmosphere acceptance](results/nongrey_t5600_acceptance_v1.json).

The [warm plan](../data/atmosphere/sources/nongrey_forward_parallel_plan_specification.json)
contains 34 source calculations, including controls and a complete composition
slice at 6000 K and log g = 5.7. Two additional independent models at 5300 K
and 5700 K test the intermediate temperature intervals. Separate
[6000 K, log g = 6 trials](../data/atmosphere/sources/nongrey_high_gravity_column_pilot_plan_specification.json)
and [contraction trials](../data/atmosphere/sources/nongrey_contraction_column_pilots_plan_specification.json)
at 4800 K, log g = 6.3 and 3800 K, log g = 6.6 test numerical and material
coverage ahead of the trajectory. These trials are not complete atmosphere
grids. A longer-term gravity range through about 6.9 would provide margin
around the inferred curve, subject to the actual stellar radius.

## Physics at the turn and during cooling

The HR-diagram turn is not a change of evolution equations. Ember must follow
the balance of nuclear heating, contraction and heat loss while the electron
gas becomes more degenerate. Degenerate electrons and conductive transport
are already present; their density coverage and accuracy need to extend with
the star. The current hydrogen-burning network is incomplete. Carbon and
nitrogen conversion, the other relevant hydrogen reactions, and thermal
neutrino channels need quantitative checks and implementation where required.
The remaining hydrogen distribution must be resolved as burning moves away
from the center. See [nuclear physics](NUCLEAR.md) and
[the remnant work order](COLD_REMNANT.md#work-and-acceptance-order).

Microscopic diffusion and gravitational separation must be coupled to the
existing mixing and burning. Hydrogen can accumulate at the surface of a
helium-core white dwarf; deeper convection can mix hydrogen and helium again.
The atmosphere must use the resulting surface composition, including metal
depletion. The current fixed GS98 metal pattern and X = 0.15–0.2 atmosphere
interval cannot describe all of these possibilities.
[Diffusion in helium white dwarfs](https://academic.oup.com/mnras/article/317/4/952/1039201),
[observed atmospheric composition changes](https://arxiv.org/abs/1905.02174).

The planned cooling boundary uses wavelength-dependent H/He atmospheres in
radiative and convective equilibrium. Their pressure, temperature, density
and thermal responses must agree with the interior equation of state at an
optically thick matching layer. Comparisons at different matching depths
must show that the calculated radius and cooling rate are insensitive to
where atmosphere and interior meet. The treatment must also follow the
connection of envelope convection to the conductive interior during cooling.

The [current gas atmospheres](NONGREY.md#physical-specification) already include
H2 collision-induced absorption. White-dwarf conditions require a consistent
treatment of nonideal ionization and molecular dissociation, pressure-dependent
line and continuum absorption, the red wing of hydrogen Lyman-alpha, and
density-dependent H2 absorption. Refraction and dense-fluid effects must be
included if the actual atmospheric densities require them. The lower gravity
of a 0.1-solar-mass remnant makes extrapolation from ordinary white dwarfs
particularly inappropriate.
[Dense atmosphere physics](https://arxiv.org/abs/1807.06616),
[collision-induced absorption tests](https://arxiv.org/abs/1809.11122),
[published atmosphere prescriptions](https://www.astro.umontreal.ca/~bergeron/CoolingModels/).

Grain formation remains necessary wherever the calculated elemental abundances
and temperature permit condensation. Its importance in the white-dwarf phase
depends on whether metals remain near the surface after gravitational
separation. Coupled grain chemistry and radiation are unfinished. The present
wavelength interval ends at 30 microns and is insufficient for a 100 K
atmosphere; far-infrared coverage and the low-temperature chemistry must also
be extended. Current molecular and collision-induced absorption data have
explicit temperature limits and cannot simply be extrapolated to that endpoint.

Much later cooling additionally needs a consistent helium free energy,
interaction and quantum contributions to heat capacity, and a calculated
phase boundary. Latent heat applies if the material crosses that boundary.
A carbon/oxygen crystallization prescription does not establish helium
crystallization or its cooling delay. These are separate requirements from
the warm atmosphere extension.

## Present-day comparisons

Binary mass loss can expose a helium core before helium ignition, producing
low-mass white dwarfs within the present age of the Universe. The companion
of PSR J0348+0432 has a measured mass of 0.172 ± 0.003 solar masses.
[Antoniadis et al. (2013)](https://arxiv.org/abs/1304.6875).
Such objects provide comparisons for radius, spectra and atmosphere physics.
Their total masses and remaining hydrogen layers differ from Ember's isolated
0.1-solar-mass star, so their ages are not direct lifetime comparisons.
Residual hydrogen burning can substantially delay cooling.
[Calcaferro et al. (2018)](https://arxiv.org/abs/1802.06753).

Published three-dimensional pure-hydrogen models cover effective temperatures
6000–11500 K and log g = 5–6.5. They provide tests of convection and
spectroscopic inferences where conditions overlap. Their spectroscopic
corrections are not themselves an interior boundary table.
[Tremblay et al. (2015)](https://arxiv.org/abs/1507.01927).
The public Montréal atmosphere grids cover log g = 7–9, and their associated
cooling sequences assume carbon/oxygen cores. They are useful comparisons
within their stated domains but do not supply the requested 0.1-solar-mass
helium-core cooling track or the atmosphere down to 100 K.
[Montréal tables](https://www.astro.umontreal.ca/~bergeron/CoolingModels/).
