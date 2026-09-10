# Forward evolution with helium enrichment

Current status, September 10: a fresh calculation using the refined 64-plane
EOS, hydrogen-poor TOPS v4 opacity, 96-node gas atmosphere and plasma-neutrino
losses reached 3.30 trillion years, then its exact restart reached **3.40 trillion
years**. It remains fully convective, at XH=.2034 and Teff=3202 K.
See [COLD_REMNANT.md](COLD_REMNANT.md) and the
[September 10 working paper](reports/2026-09-10/ember_status_and_future.pdf)
for current input validation and the path to 100 K remnant cooling.

**The detailed numerical comparisons below describe the preserved earlier
2.85-trillion-year reference and its original input families.** Their
repeatability and two-trillion-year convergence results remain useful controls;
they are not convergence measurements of the later transition. The new
108-node atmosphere candidate is undergoing composition refinement. No
condensate family, hydrogen-exhaustion result or cooling endpoint is accepted.

## Elemental mixture and thermodynamics

`MetalHelmholtzEos` interpolates a material Helmholtz potential in temperature,
density, baryonic hydrogen fraction and helium-3 fraction. The installed
family is `data/eos/freeeos300_gs98_z020.dat`: hydrogen planes
.3/.4/.5/.6/.7/.75 and helium-3 planes 0/.12, with fixed Z=.02. Each source
plane contains the same GS98 representative-element numbers used by the
non-grey atmospheres. FreeEOS 3.0 has no potassium entry; its omitted baryonic
mass fraction, 4.56e-6, is recorded rather than reassigned to another element.

The source computes electronic helium physics with He4. Baryonic mass
conversion preserves the requested H/He/metal number densities; ideal helium
isotope mixing and translational entropy are added analytically. Isotope
energy-level shifts and spin constants are not included. Radiation is added
once, after interpolating material potentials. Pressure, energy, entropy,
heat capacities and thermal derivatives all derive from that one potential.

The source grid is log10(T/K)=3.5..7.1 and
log10[rho/(T/1e6 K)^1.5]=-1.5..2.5, at .0125 spacing. Imported derivative
stencils inset the outer support by .025. Nonconverged source states remain
flagged; every intersecting derivative stencil is masked. Density inversion
uses only the contiguous supported branch attached to low density. There is
no interpolation across an invalid stencil or extrapolation beyond support.

The four composition-weighted node stencils are combined before evaluating
the common Hermite basis. Subtracting a local constant from F/T before
differentiation avoids cancellation of the arbitrary entropy zero; it does
not change the physical potential. In the regression states, independent
first-law and response errors are 3.76e-10 and 4.15e-8. The 144 fresh source
comparisons in `tests/data/freeeos300_gs98_reference.dat` have maximum
differences .06164% in pressure, .1675% in energy and .2595% in heat capacity.
These local interpolation checks are not global physical uncertainty bounds.

A separate hydrogen-poor family is now available as
`data/eos/freeeos300_gs98_hydrogen_poor_z020.dat`. It adds XH=.1/.125/.15/.175/.2/.25,
retains all original XH>=.3 planes byte for byte, and uses the same X3=0/.12
axis and material grid: 24 planes in total. Coarse .1 spacing gave up to
1.22% heat-capacity interpolation error in cool hydrogen-poor gas, so the
additional composition planes were computed before acceptance. In 336 fresh
FreeEOS comparisons, maximum differences are .04985% in pressure, .1087%
in energy and .1934% in heat capacity, including the unchanged XH>.3 interval.
Inside the new XH<.3 region the corresponding maxima are .03491%, .09190%
and .1548%. First-law and response errors are below 7.2e-10. All 1536 states
from the completed 2/2.5/2.75T profiles give byte-identical EOS probe output
under the original and extended families. Original source slices, raw new
planes and assembly manifests are retained. See
`results/metal_eos_hydrogen_poor_audit.json`,
`results/metal_eos_hydrogen_poor_coarse_audit.json` and
`results/metal_eos_hydrogen_poor_overlap.json`. The existing stellar tracks
continue to select the original family explicitly.

At the tau=100 matching states of all 48 existing atmospheres, the new EOS
gives densities .738–1.992% above the atmospheric gas EOS. The independent
comparison is `docs/results/nongrey_gs98_eos_family_audit.json`. Matching
T/Pgas remain the boundary conditions; density is recomputed from the
interior EOS. The gas EOS and partition-function differences remain a
physical limitation despite the improved element matching.

`MetalInventory::gs98` selects the same 19-element inventory in ion/electron
counts, nuclear screening and electron conduction. The existing inert metal
slots carry conserved Z; their names do not become a resolved GS98 isotope
network. The Ioffe conduction imports now retain eleven original charge
planes through Zn, bracketing Ni. Source table values are unchanged. See
[the original conduction tables](https://www.ioffe.ru/astro/conduct/condint.html)
and [CONDUCTION.md](CONDUCTION.md) for their fully ionized approximation and
the neutral-envelope join.

Interior opacity still uses the original AESOPUS/TOPS fixed GS98 pattern.
The wrapper converts the aggregate metal mass and the exact H/He number
counts to source grams. The source's atomic-weight metal pattern differs
slightly from the atmosphere/EOS representative-isotope pattern: individual
metal number ratios differ by up to about 1.13%. This residual is explicit;
the mixture is substantially better matched, but not isotope-resolved.

A separate hydrogen-poor opacity directory now retains the old source
planes and adds TOPS X=.09/.1/.2 at Z=.01/.02/.03. The .09 plane brackets
the atomic-mass mapping of baryonic XH=.1 with He3=.12. Implementation
checks preserve all 1536 archived profile queries byte for byte and verify
new composition responses and source-domain rejection. An independent
X=.15,Z=.02 source calculation exposes up to 11.90% interpolation error in
cool, dense states, motivating further composition refinement before use.
See [the opacity data notes](../data/opacity/README.md) and
`results/opacity_hydrogen_poor_heldout.json`. The current track still uses
the original opacity directory.

The older H/He EOS with metals represented by helium remains an explicit
control. It rejects the GS98 inventory tag. Changing inventory during burning
or mixing is rejected, and atmosphere cache keys include the inventory.

## Convection and composition-gradient transport

Ledoux buoyancy uses two EOS density inversions at the same face temperature
and pressure:

\[
B={\ln\rho(T,P,X_\mathrm{hi})-\ln\rho(T,P,X_\mathrm{lo})
 \over \delta\,[\ln P_\mathrm{hi}-\ln P_\mathrm{lo}]}.
\]

This becomes dln(mu)/dln(P) for an ideal gas. Full-EOS thermal derivatives
enter the analytic structure Jacobian. Both the MLT temperature gradient
and the convective mixing partition use this same buoyancy. A homogeneous
composition gives exactly zero B and recovers the existing Schwarzschild
calculation. Composition discontinuities with unresolved pressure contrast
are rejected rather than hidden by a denominator floor.

Two optional slow mixing prescriptions operate on Ledoux-stable faces. With
thermal diffusivity chi=4acT^3/(3 kappa rho^2 cp), including conduction in the
total diffusive opacity:

- Langer mixing-only semiconvection:
  D=alpha_sc chi/6 (grad_rad-grad_ad)/(grad_ad+B-grad_rad), for
  B>0 and grad_ad<grad_rad<=grad_ad+B.
- Kippenhahn thermohaline mixing:
  D=-1.5 alpha_th chi B/(grad_ad-grad_rad), for an inverted composition
  gradient B<0 that remains Ledoux stable.

These choices follow the corresponding implementations in
[MESA's semiconvection source](https://github.com/MESAHub/mesa/blob/fd396fd73d3f936da8063ffdf9d92361882eb557/turb/private/semiconvection.f90)
and [thermohaline source](https://github.com/MESAHub/mesa/blob/fd396fd73d3f936da8063ffdf9d92361882eb557/turb/private/thermohaline.f90).
The semiconvective option mixes composition without adding a separate heat
flux prescription. Efficiencies are model parameters, not calibrations.

Burning and diffusion are solved together by backward Euler in conserved
baryonic mass coordinates, with zero surface/central diffusion flux. Fully
convective regions collapse to homogeneous abundance unknowns. A block
tridiagonal solve retains the mass term even for extremely strong diffusion;
positivity and integrated reaction balances are checked independently.
Composition and thermal structure are iterated to convergence, recomputing
transport after the structure changes. Zero diffusion recovers the previous
burn/mix solver exactly.

Tests cover composition barriers, conservative fluxes, positivity, the strong
diffusion limit, nonlinear burning, both slow-mixing regimes and the full
analytic structure Jacobian. This does not yet validate a moving radiative
core in an evolved star. Microscopic settling, mesh adaptation and resolved
moving-boundary convergence remain future work.

## Atmospheres and condensates

The separate condensate production rectangle now has 70 of 72 cells accepted.
The two outstanding cells are solar XH=.7, Teff=2800 K, logg=4.9 at X3=0/.12.
One is undergoing numerical recovery; the other converges but fails the
grain-enthalpy guard because condensing layers carry 4.17e-7 of the total
flux by convection, above the unchanged 1e-8 threshold. Two logg=5.0 controls
also fail this physical check. Artificial heat-capacity controls measure
the boundary sensitivity; they do not supply a full grain EOS. No condensate
runtime family has been installed and no condensate stellar feedback has
entered the 2.85-trillion-year track.

The installed atmosphere extension adds XH=.3 to the existing .45/.7 planes,
retaining X3=0/.12, Teff=2600/2800/3000/3200 K and logg=4.9/5.15/5.4.
All six opacity planes and 72 atmosphere cells pass their source checks.
Original inputs and outputs are archived in
`data/atmosphere/sources/nongrey_gs98_z020_extended/`; the separate runtime
file is `data/atmosphere/nongrey_gs98_z020_extended_tau100.dat`, SHA256
`eec2e0553fb17e17fe16b82260a04977de6e0058d03b272992e638a1150d9532`.
The old 48-cell family remains available. The extended family audit is
`docs/results/nongrey_extended_gs98_eos_family_audit.json`; interior EOS
densities at the matching states are .6785–1.991% above the atmosphere EOS.
Final atmosphere cells require the
unchanged correction, flux, hydrostatic, chemical-density and full-profile
source-support checks before installation. Optional native CONREF refinement
and modified trial deep gradients provide starting guesses only; a separate
canonical solve supplies each accepted atmosphere. See [NONGREY.md](NONGREY.md)
for the existing validated 48-cell family.

An independent atmosphere at XH=.375, X3=.06, 3000 K and logg=5.15 now
checks interpolation between the new helium-rich composition corners.
After 14 native-CONREF initialization iterations and two canonical source
iterations, interpolated T differs by +.1958%, Pgas by +.03401%, and source
density by +.1295%. The four corner models and the independent calculation
were revalidated against their original input/output receipts. Original
artifacts and the local comparison are in
`data/atmosphere/sources/nongrey_extended_validation/` and
`docs/results/nongrey_extended_heldout.json`. This is a local comparison,
not a global error bound for the extended family.

The material rectangle for the extension is 1076..7762.30 K. It retains
the original low-temperature isotherms exactly and brackets the converged
cool profiles. A proposed 15000 K rectangle exposed negative native He I
4471 profile extrapolation above its tabulated electron-density range. No
ad hoc broadening correction or negative-opacity clipping was introduced.

Pinned FastChem 4 sources now provide offline equilibrium and bottom-to-top
rainout chemistry. All 48 existing atmospheres were checked in both modes,
including independent element, particle-pressure and nuclei closure. Original
outputs, chemistry files and license are archived under
`data/atmosphere/sources/condensation_gs98/`; the report is
`docs/results/condensation_nongrey_family.json`. Maximum element closure
error is 3.66e-9. Condensation reaches tau=.01356 in 2600 K models and .001784
in 2800 K models; it is absent at 3000/3200 K and at tau>=.1 throughout the
audited family. This fixed-profile result alone does not establish negligible
radiative feedback. The chemistry methods are described in
[FastChem Cond](https://arxiv.org/abs/2309.02337).

A separate source experiment couples local equilibrium gas depletion to
both the atmospheric gas EOS and the opacity calculation. It uses the
efficient-settling endmember: condensed elements are removed from the gas,
and grain opacity is zero, as in the limiting calculation described by
[Allard et al. (2001)](https://arxiv.org/abs/astro-ph/0104256).
The offline GPL adapter is separate from Ember's MIT library. It recomputes
depletion at trial T/P, including derivative calls, and independently checks
FastChem element and pressure closure. It combines FastChem's removal
fractions with TLUSTY/SYNSPEC's original gas partition functions and line
data; this is a hybrid chemistry approximation. Grain enthalpy, grain sizes,
Mie opacity and time-dependent settling are not implemented.

One coupled atmosphere now passes all original-source checks: XH=.7, X3=0,
Teff=2600 K, logg=5.15, with 300 depths and 20000 transfer frequencies.
Its opacity table has 17 temperatures from 1001 to 7762 K, entirely
bracketing the final 1054–4528 K profile. At tau=100 it gives
T=3864 K and Pgas=2.077e+7 dyn/cm2. The maximum flux error is
4.72e-5 and final undamped temperature correction 5.12e-7. Independent
element closure is 2.97e-9; all condensing layers have zero convective heat
flux. Condensation reaches tau=.00631, well above the matching boundary.
The original model and chemistry audit are archived under
`data/atmosphere/sources/condensation_atmosphere_x700_2600_g515/`, with a
compact report in `docs/results/condensation_atmosphere_cold.json`.
A converged gas control on the identical material and frequency grid isolates
the depletion effect: matching pressure increases 0.8095%, temperature falls
0.1349%, and density increases 1.104%. Changing the gas-only material grid
from 16 to 17 temperatures accounts for a separate 0.0740% pressure change.
The control's original artifacts are archived with the condensed atmosphere.
This resolved boundary effect motivates a complete composition-dependent
condensate family; the single atmosphere is not installed as a stellar boundary.

A helium-rich control (XH=.3, X3=.12, 2800 K, logg=5.15) also passes,
using the original 16-temperature material grid. Relative to the gas model
on that grid, matching pressure increases .4833% and temperature falls
.0974%. Its maximum element closure error is 3.56e-9; condensation reaches
tau=.00258 and all condensing layers have zero convective heat flux. Original
artifacts are in `data/atmosphere/sources/condensation_atmosphere_x300_he3120_2800_g515/`;
see `docs/results/condensation_atmosphere_helium.json`.

The solar 2800 K model also passes, after 28 canonical iterations, with
independent chemistry and original artifacts in
`data/atmosphere/sources/condensation_atmosphere_x700_2800_g515/`.
Against the same 17-temperature gas opacity grid, its matching pressure
increases .1246% and temperature falls .02347%. The matched-source audit
and comparison are reproducible with `scripts/compare_condensate_control.py`;
see `docs/results/condensation_atmosphere_warm.json`.

`generate_condensate_grid.py` now runs composition/gravity continuation chains,
and `archive_condensate_grid.py` assembles a complete rectangular family.
The condensate importer independently rechecks original inputs, executable
fingerprints, material support, mixture, chemistry species conservation and
the absence of heat-carrying condensing layers. It rejects initializers and
mislabeled gas controls. Offline archive and altered-input checks are recorded
in `docs/results/condensation_import_validation.json`. All six 17-temperature
opacity planes are computed and archived; atmosphere cells remain in progress.
The attempted 126-cell rectangle includes 2600/2650/2700/2750/2800/3000/3200 K.
Several low-hydrogen, high-gravity 2600 K calculations encounter trial layers
below the source's 1000 K partition-function floor and are rejected. The
explicit forward specification now uses 2750/2800/3000/3200 K, retaining all
composition and gravity axes, for 72 cells. It brackets the star's initial
2768 K. Solved colder nodes remain diagnostic controls; unsupported cold
nodes are never filled or extrapolated. Independent 2775/2900 K models will
check the production grid, with 2725 K as a colder continuation control.
A separate 33-temperature opacity experiment retains the original 17
isotherms exactly and adds their midpoints to check material interpolation.
Only the source array capacity changes for this refinement; two controls on
the original 17-temperature opacity grid agree with the original atmosphere
boundaries to 1.56e-6 relative. The cold 2600 K material refinement passes:
matching T changes +.02652% and pressure +.3605%. The warm 2800 K source
solve converges, but its chemistry/transport audit flags a convective flux of
1.224e-6 of the total in one condensing layer, above the deliberately strict
1e-8 no-grain-enthalpy gate. The fine model remains a diagnostic, explicitly
ineligible for a production grid. Both production import and collection keep
the strict gate. No condensate grid has yet been installed in Ember.

The flagged convective flux does not result from cancellation in the native
mixing-length root: its rationalized form agrees to 4.44e-16 relative. At the
actual native interface state (1769 K, 4062 dyn/cm2), an independent
FastChem formation-enthalpy derivative estimates an additional heat capacity
of 2.89e5 erg/g/K, about .187% of the native gas value. This is an
ideal-neutral-atom reference estimate, not a complete gas/grain EOS or a
pseudoadiabat. Atomic excitation and source thermodynamic differences remain
omitted; see `results/condensation_warm_latent_enthalpy_estimate.json`.

Two separate source experiments add artificial positive heat capacities of
1e6 and 1e7 erg/g/K below 2000 K, tapering smoothly to zero at 2400 K,
with the corresponding adiabatic-gradient change. Matching pressure changes
by 8.39e-7 and 4.39e-6 relative, respectively. They test boundary sensitivity,
not a rigorous bound on missing grain physics. Archived original inputs and
outputs are independently rechecked in `audit_condensate_material.py` before
using the fine model as an interpolation diagnostic: T changes +.03271% and
pressure +.2890%, within the .5% material-resolution criterion. The report
explicitly retains its failed production grain-enthalpy status. Each actual
production cell must separately pass the unmodified strict grain gate. See
`results/condensation_warm_material_refinement.json` and
`results/condensation_warm_heat_capacity_controls.json`.

A sampled CONREF initialization phase spends most of its time in native
molecular-equilibrium element searches. An optional indexed search preserves
the arithmetic expressions and update ordering. In the 2750 K control, its
independent canonical replay agrees with the original final boundary to
1.74e-8 relative in T/P/density. Original equations and final acceptance checks
remain mandatory. This is a numerical initializer, not a different atmosphere
prescription; see `results/condensation_indexed_initializer_validation.json`.

The archived cold, warm and helium-rich controls can be compared in
[the atmosphere profile figure](results/condensation_atmosphere_controls.pdf),
reproduced by `scripts/plot_condensate_controls.py`.

The gas-only control reproduces an original 30000-frequency, 19-density
opacity isotherm byte for byte. Its 300-depth atmosphere converges in one
iteration and changes matching T/Pgas by less than 1.1e-9 relative. The
depletion bridge agrees with an independent 300-layer FastChem calculation
to 2.24e-11 of each element's initial abundance. A 65-pressure scan at each
of 3000/4000/5000/6000 K finds no condensates over 1e-13..1000 bar; the
6000 K boundary supports the explicitly fully vaporized hotter continuation.
For atmospheres containing hotter layers, the independent audit now queries
their actual pressures separately at 6000 K and requires zero appreciable
condensate particles there. Actual-state chemistry is evaluated only inside
the tested domain. The XH=.3, X3=0, Teff=3200 K, logg=4.9 control has
19 layers at 6020–7157 K; its 6000 K join queries contain no condensates and
close the element inventory to 7.56e-13. The hotter continuation is labelled
as the source prescription, never as an extrapolated FastChem result. See
`results/condensation_vaporized_atmosphere_join.json`.
Original controls and the output-filename receipt recovery are archived in
`data/atmosphere/sources/condensation_validation/`; see
`docs/results/condensation_source_validation.json`.

## Current stellar checkpoint and usage

The latest 512-point gas run reaches **2.85 trillion years** with
R=.1488 Rsun, L=.001858 Lsun, Teff=3107 K,
central T=6.879 MK and rho=241.2 g/cm3. XH=.3029 and
X3=.004663; the star remains fully convective and homogeneous.
The native continuation added 37 accepted steps with no rejections, bringing
the total to 1673 accepted steps and one rejection. The
[model report](results/evolution_metal_m010_2850gyr_gas.json) and
[figure](results/evolution_metal_m010_2850gyr_gas.pdf) contain the
2.75-to-2.85T continuation interval. The receipt confirms unchanged data and
restart inputs. This still uses the original gas boundary and EOS family;
XH=.3029 is close to their .3 source boundary. No late-segment refinement or
condensate stellar feedback is claimed.

The preceding 512-point gas run reaches **2.75 trillion years** with
R=.1489 Rsun, L=.001826 Lsun, Teff=3092 K,
central T=6.720 MK and rho=240.4 g/cm3. The homogeneous baryonic
fractions are XH=.3198, X3=.005907, X4=.6543 and Z=.02.
Central pressure is 1.645e+17 dyn/cm2 and surface logg=5.093.
The star remains fully convective; the minimum diffusive/adiabatic gradient
ratio is 3.66 and the maximum estimated local conductive flux fraction is
3.22%. The largest reduced ppII/pp ratio is .002254; plasmon-only losses
are 1.41e-8 of surface luminosity.

The [restart continuation](results/evolution_metal_m010_2750gyr_gas.json)
and its [figure](results/evolution_metal_m010_2750gyr_gas.pdf) contain only
the 2.5-to-2.75-trillion-year interval. This segment took 94 accepted steps
with no rejections, bringing the entire calculation to 1636 accepted steps
and one rejection. Transport and nuclear audits accompany the report;
the receipt confirms unchanged data and restart inputs. Condensate feedback
is still absent. Mesh and timestep refinements have reached two trillion
years; the small late radius turnover is not yet established by refinements
of this later segment. The next composition extension must go below the
current XH=.3 source boundaries to follow substantial further fuel depletion.

The preceding 512-point gas run reaches **2.5 trillion years** with
R=.1487 Rsun, L=.001753 Lsun, Teff=3062 K,
central T=6.375 MK and rho=239.9 g/cm3. XH=.3600 and
X3=.01015; the model remains fully convective. It took 1542 accepted
macrosteps and one rejection. The minimum diffusive/adiabatic gradient ratio
is 4.97 and the largest estimated local conductive flux fraction is 3.09%.
Plasmon-only losses are 1.35e-8 of surface luminosity; the largest reduced
ppII/pp ratio is .001203. Nuclear rest-mass release integrated from the
initial composition is .002392 of baryonic mass. This remains a fixed
baryonic-mass calculation with the nuclear and atmosphere approximations
listed above, without condensate feedback.

The [compact 2.5-trillion-year track](results/evolution_metal_m010_2500gyr_gas.json),
[figure](results/evolution_metal_m010_2500gyr_gas.pdf), transport audit and
nuclear audit record that completed segment. Resolution checks are complete
through the earlier two-trillion-year endpoint.
The new [restart support](RESTART.md) passes exact trajectory and input-hash
checks. A fresh checkpointed 2.5-trillion-year repeat reproduces the original
complete stellar output byte for byte, as recorded in
`results/evolution_metal_2500gyr_checkpoint_repeat.json`. Its saved state
supplies the completed continuation above.

The preceding 512-point run reaches **2 trillion years** with
R=.1478 Rsun, L=.001632 Lsun, Teff=3017 K,
central T=5.834 MK and rho=242.5 g/cm3. XH=.4309 and
X3=.02575; the star remains fully convective. The new XH=.3 atmosphere
plane supports the portion below the previous XH=.45 boundary. There are
1339 accepted macrosteps and one rejection. The minimum diffusive/adiabatic
gradient ratio is 8.24, with maximum estimated local conductive flux fraction
2.92%. Plasmon-only losses are 1.29e-8 of surface luminosity. Condensate
feedback is not yet included. An independent fresh-start repeat gives
byte-identical full history and final profile, as recorded in
`results/evolution_metal_2tyr_gas_repeat.json`. Doubling the mesh from 512 to
1024 points changes final L by +.06163%, R by -.01234%, and XH by -1.214e-4;
see `results/evolution_metal_m010_2tyr_gas_convergence.json`. Fourfold tighter
time tolerances change L by -.003610%, R by -.000594%, and XH by +1.300e-5.
These quantify repeatability and numerical
sensitivity under the recorded inputs and
executable, not physical accuracy or hydrogen-exhaustion convergence.

The [compact track and full provenance](results/evolution_metal_m010_2tyr_gas.json),
[figure](results/evolution_metal_m010_2tyr_gas.pdf), and companion transport and
nuclear audits record this checkpoint. The copied executable used 1632
child CPU seconds and 2365 UTC seconds while source jobs ran concurrently;
these are not isolated-machine benchmarks.
The integrated nuclear rest-mass release is .001810 of the conserved baryonic
mass. Energy bookkeeping includes that release, while the Newtonian gravity
calculation retains the baryonic mass coordinate. A same-mass comparison with
another stellar code must account for this convention and approximation.

The earlier 512-point control reaches **1.3 trillion years** with
R=.1462 Rsun, L=.001465 Lsun, Teff=2953 K,
central T=5.306 MK and rho=248.3 g/cm3. XH=.5057 and
X3=.07072; the star remains fully convective. The minimum
diffusive/adiabatic gradient ratio is12.99 and maximum estimated local
conductive flux fraction2.84%. Maximum pp screening zeta=.07284;
plasmon-only losses are1.30e-8 of surface luminosity. The run uses
`ledoux-diffusive`; no slow transport activates in this homogeneous star.

The complete input/executable hashes, final quantities and conservation
checks are in `docs/results/evolution_metal_m010_1300gyr.json`, with
companion transport/nuclear audits. Python perf_counter measured1418seconds of awake elapsed time while
atmosphere-source jobs ran concurrently. Start and receipt-write timestamps
span about2973seconds; system suspension is excluded by that timer on this
host. Neither interval is a standalone performance benchmark.
The same executable has also completed a 1024-point run and a 512-point run
with fourfold tighter timestep error tolerances. Doubling the mesh changes
final R by -.01022%, L by +.01917% and He3 by -7.51e-5 absolute. Tightening
the timestep tolerances changes R by +.000466%, L by +.000620% and He3 by
+8.59e-6 absolute. All three remain fully convective. These local sensitivity
checks do not establish an asymptotic mesh order or physical uncertainty.
Full hashes and comparisons are in
`docs/results/evolution_metal_1300gyr_convergence.json`. The earlier 1T
checkpoint below retains its original executable/physics provenance.

At one trillion years, the new 512-point calculation has R=.1442 Rsun,
L=.001369 Lsun, Teff=2923 K, central T=5.154 MK and
rho=257.7 g/cm3. Hydrogen and helium-3 fractions are .5331 and
.09106. It remains fully convective. There are 925 accepted macrosteps
and one rejection. This uses the new GS98 interior with the existing gas
atmosphere family. Condensate feedback is not part of this checkpoint.

The compact result, transport audit and nuclear audit are in
`docs/results/{evolution,transport,nuclear}_metal_m010_1tyr.json`; the full
track is `out/evolution-metal-512-1tyr.json`. The executable hash and wall
time were not recorded at that launch and are not retrospectively claimed.

Select the new interior and Ledoux transport explicitly:

```sh
build/apps/ember-evolve 512 1e12 1e8 10 sfii-svh wd \
  nongrey:data/atmosphere/nongrey_gs98_z020_tau100.dat \
  metal:data/eos/freeeos300_gs98_z020.dat ledoux
```

`ledoux-diffusive` adds alpha_sc=.1 and alpha_th=1. `schwarzschild` and
`proxy` select the earlier controls. The default atmosphere/interior choices
are unchanged. All 27 CTest suites pass at this checkpoint. Source-generation
experiments have separate acceptance tests; building the C++ library does
not validate or install an unfinished source grid.
