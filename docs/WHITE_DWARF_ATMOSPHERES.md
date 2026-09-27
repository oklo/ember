# Atmospheres through contraction and white-dwarf cooling

Ember evolves the structure, composition and thermal energy of an initially
0.1-solar-mass star. Its atmosphere supplies pressure and temperature at a
deep matching layer. Cooling requires this physical boundary as well as an
emergent spectrum; a table of colors alone cannot supply it.

The continuous Hayashi-origin calculation uses composition-dependent
atmospheres matched at optical depth **100** throughout. Its nearly pure-H
extension covers **4600–6000 K** and **log g = 5.9–6.2**, with a declared
trace-helium approximation. The extension preserves the retained stellar
trajectory exactly; coverage at still higher gravity is being prepared.
[Current trajectory](reports/2026-09-27/pms_figure_inputs.json),
[boundary checks](results/highg_atmosphere_extension_sept27_v2.json).

## Supplied white-dwarf boundary comparison

The separate flash and cooling comparisons use a supplied MESA hydrogen
white-dwarf table. Its gas pressure and temperature refer to optical depth **25.12**,
cover **2000–40,000 K** and **log g = 5.5–9.5**, and preserve the original
source values. Ember interpolates these values directly; it does not apply
MESA's transition to a gray boundary below **log g = 6**. The interior retains
its actual species abundances, while the atmosphere uses a pure-hydrogen
approximation limited to helium fraction **0.005** and **Z ≤ 1e-20**.
[Importer](../scripts/import_mesa_wd_atmosphere.py),
[source and runtime checks](results/published_wd_boundary_runtime_v1.json).

All source-node comparisons, finite-difference derivatives and sampled
domain-rejection tests pass. The stellar EOS supports **24 of 25** sampled
boundary states; **2500 K** at **log g = 6.75** needs additional EOS coverage.
This is a finite coverage test, not validation of the entire atmosphere grid.
At **4600 K**, **log g = 5.9** and the same optical depth, the supplied
table gives **2.594%** higher matching temperature and **48.69%** lower gas
pressure than the independently converged Ember atmosphere. These are
differences between atmosphere prescriptions; their detailed cause has not
been isolated.

A matched **1 Myr** stellar test changes effective temperature by about
**241 K** when this boundary replaces the selected optical-depth-100
atmosphere. All native stellar solves and conservation checks pass, but the
step exceeds the original time-error allowance. The comparison histories
therefore resolve an initial thermal adjustment with shorter intervals.
Their helium-3 pulse follows this atmosphere change and does not establish
ignition with a consistent boundary history. This table is not selected for
the continuous Hayashi-origin trajectory.
[Finite stellar comparison](results/published_wd_stellar_control_v1.json).

## Conditions to calculate in advance

Source planning follows the calculated surface composition, temperature and
gravity together. The nearly pure-hydrogen surface changes the relevant
opacity and chemistry requirements. The existing composition-dependent
atmosphere families support the hotter and more metal-rich portions of the
track; they are not a substitute for dense, cool hydrogen atmospheres.

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

The checked warm sources include the X = 0.15–0.20 region through 6400 K
and log g = 5.4, and its extension to log g = 5.7/6 at 5400–6400 K.
The first region has 252 source structures and a largest independent
matching-state difference of 0.4336%; the higher-gravity extension adds 48
structures, with six independent comparisons within 0.4860% and 17 passing
depth controls. Those results apply to their stated compositions and domains.
[6400 K checks](results/nongrey_t6400_acceptance_v1.json),
[higher-gravity checks](results/nongrey_warm_gravity_acceptance_v1.json).

Separate trials at 6000 K and log g = 6, 4800 K and log g = 6.3, and 3800 K
and log g = 6.6 pass their local lower-boundary comparisons. They test the
contraction and cooling region inferred from LBA97 but do not constitute a
complete cooling grid or establish the validity of the gas prescription at
every future atmosphere density. A gravity range through about 6.9 would
provide margin around the inferred curve, subject to the calculated radius.
[6000 K trial](results/nongrey_x150_y000_t6000_g600_column_v1.json),
[4800 K trial](results/nongrey_x150_y000_t4800_g630_column_v1.json),
[3800 K trial](results/nongrey_x150_y000_t3800_g660_column_v1.json).

## Physics at the turn and during cooling

The HR-diagram turn is not a change of evolution equations. Ember must follow
the balance of nuclear heating, contraction and heat loss while the electron
gas becomes more degenerate. Degenerate electrons and conductive transport
are already present; their density coverage and accuracy need to extend with
the star. The selected network evolves hydrogen, helium-3, carbon-12, carbon-13 and
nitrogen-14, with helium-4 closing the physical mass fractions. Oxygen
branches, individual catalyst settling and dense-plasma reaction corrections
remain separate questions. Thermal neutrino losses are present; their source
domains still need to be checked as the remnant evolves.
The remaining hydrogen distribution must be resolved as burning moves away
from the center. See [nuclear physics](NUCLEAR.md) and
[the remnant calculation](COLD_REMNANT.md#objective-and-next-milestones).

Conservative microscopic H/He transport and material enthalpy are coupled
to burning and mixing. Carbon, nitrogen and the remaining metals share the metal-group drift velocity.
Hydrogen can accumulate at the surface of a
helium-core white dwarf; deeper convection can mix hydrogen and helium again.
The atmosphere must use the resulting surface composition, including metal
depletion. The common metal-group velocity does not resolve differential separation
of individual heavy elements; atmosphere composition coverage remains finite.
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
