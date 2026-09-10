# The 0.1-solar-mass calculation through remnant cooling

The user has authorized autonomous implementation and calculation through a
very cold white-dwarf phase. This is the objective, not the current result.
The preserved reference remains the 2.85-trillion-year, fully convective gas
track described in `HANDOFF.md`. Keep its executable, checkpoint and physical
inputs intact; develop in a separate build and use separate output paths.
Generated bulk tables and raw source calculations remain local under the
publication policy in [DATA_REPRODUCTION.md](DATA_REPRODUCTION.md).

## Progress on 2026-09-10

The latest 512-point and 1024-point calculations both completed **3.560 trillion
years**. Both now have central regions stable against convection, containing
37.23% and 37.77% of their mass. Final central hydrogen fractions are 0.1518
and 0.1491; surface temperatures are 3240 K and 3241 K. The luminosity difference
is 0.07813%, but the central hydrogen difference is 1.752%, requiring further
mesh and time-step checks. Hydrogen burning continues. No stellar jobs remained
active at the 16:43 UTC process inspection. The new 512-point checkpoint is
intact; the 1024-point continuation wrote no new checkpoint. See
`results/evolution_transition_3560gyr_v1.json` and the updated working paper.

The paragraphs below preserve the preceding transition diagnostics.

The latest gas-atmosphere table contains 120 accepted source models, including
hydrogen fractions 0.1 and 0.15. Independent temperature checks near fractions
0.125, 0.15 and 0.175 differ by 0.7971%, 0.6468% and 1.089%. A refined
512-point trial has a saved state near 3.548 trillion years with an off-center
stable shell; its center remains convective. Driver and checkpoint limits on
cumulative rejected steps interrupted the trials and have now been corrected.
The expanded restart test and all 32 CTest suites pass. Fresh evolution and an
exact continuation of the earlier 1024-point model are checking the transition.
See `results/evolution_convection_transition_v1.json` and the latest handoff.

[Grain absorption and scattering](GRAINS.md) have been implemented as an offline
optical calculation and checked against independent Mie and small-particle
results. A fixed 2800 K profile control includes its calculated alumina mass.
Atmospheric coupling, complete material coverage, grain phase/size/settling
and condensate thermodynamics remain unfinished. The evolving star still uses
gas-only atmospheres. The [working paper](reports/2026-09-10/ember_status_and_future.pdf)
now compares Ember primarily with values and figure curves from LBA97, while
keeping the reconstructed F77 results as a separate supplementary comparison.

The following paragraphs preserve the earlier completed reference and the
measurements that led to the latest refinement.

The fresh gas track completed3.30 trillion years, followed by an exact restart
to3.40 trillion years. At3.40T the512-zone star is still fully convective:
XH=.2034, X3=.0009703, Teff3202K, Tc8.031MK,
rhoc258.3g/cm3. No exhaustion or cooling milestone has been reached.
The96-node gas atmosphere, refined EOS/TOPSv4 inputs and copied executable
remained fixed. See `results/remnant_endpoints_3400gyr_forward_v1.json`.

The candidate108-node atmosphere now includes all twelve XH=.1 warm nodes.
Reassembly, EOS support and all108 runtime stencils pass;112 queries retain
the prior warm grid's values and derivatives exactly. Two high-gravity3200K
depth comparisons span actual optical depths165/339 and170/350 and pass
at maximum matching-state differences7.13e-6 and1.93e-5. Independent
XH=.15/3300K/g5.1 interpolation differs+2.134% in matching temperature and
-.7868% in gas pressure. Four completed controls separate this into composition-only +1.652% T
and -.7842% Pgas, plus temperature/gravity +.4743% T and -.00254% Pgas.
Those additional source nodes and independent XH=.125/.175 checks are now complete.
The separate stellar trial using this table began losing full convection near
3.543 trillion years. It stopped near 3.544 trillion years after reaching the
program's total rejected-step limit, with convective mass fraction 0.7638.
It did not reach its requested age of 3.60 trillion years. These preliminary
results need numerical checks and are separate from the completed 3.40T reference.
See `results/nongrey_x015_interpolation_separation_v1_audit.json` and the
latest `HANDOFF.md`. This is gas coverage, not a cold-remnant boundary.

The offline current-profile gradient ratio has a minimum1.336 at mass
fraction.211, down from1.588 at3.30T. The central ratio is3.090.
These are diffusive-to-adiabatic nodal gradient ratios, including conduction;
unity is the homogeneous stability threshold. They suggest approaching loss
of full convection without establishing a central radiative-core onset or
an extrapolated transition age. Report:
`results/evolution_convection_margin_3400gyr_v1.json`.

The retained F77 source was reproduced separately with a print-only diagnostic.
Its AJR-opacity option first flags the center nonconvective at3.224T,
XH=.01756; its current Ferguson option does so at2.946T,
XH=1.064e-5 and later loses convergence during cooling. These are explicit
source/prescription comparisons, not the original1997 paper's transition.
The AJR diagnostic and an unmodified-source control have identical original
output bytes. At equal XH~.203, AJR F77 is about75% more luminous than Ember,
with similar radius and central temperature. The model clocks must not be
identified. Details: `results/f77_core_onset_comparison_v1.json`,
`results/f77_ember_matched_hydrogen_v1.json`, and
[COMPUTATIONAL_COST.md](COMPUTATIONAL_COST.md).

## Endpoints and reference assumptions

Continue the initially homogeneous, isolated, fixed-baryonic-mass 0.1 Msun,
XH=.7, X3=0, Z=.02 GS98 reference. The clock starts at the specified static
main-sequence model. Rotation, magnetism, mass loss and external heating are
separate physical sensitivity studies, not silently assigned prescriptions.
Assess the change between baryonic and gravitating mass before quoting a
precision lifetime.

Record central XH crossings at 1e-3, 1e-4 and 1e-5, retaining the bracketing
accepted states. Separately record when deposited hydrogen-burning luminosity
falls below 10%, 1% and 0.1% of the surface photon luminosity on the cooling
branch. A crossing alone does not establish permanent extinction: retain
later recrossings, residual fuel, shell luminosity and integrated energy.
Do not relabel the first core threshold as the end of all hydrogen burning.

Cold-remnant temperature milestones are Teff=1000, 500 and 100 K, on the
descending-temperature cooling branch. The final target is 100 K, conditional
on a supported material EOS and atmospheric boundary. Zero temperature is
not an attainable finite-age endpoint. Neither a Mestel-law continuation nor
an extrapolated atmosphere counts as the requested evolved model.

## Work and acceptance order

1. **Main-sequence input coverage.** Refine and independently check the
   hydrogen-poor TOPS family, including the active AESOPUS/TOPS blend. The
   refined EOS now reaches XH=0 in composition, within its original thermal
   domain. Calculate and check gas atmospheres
   below XH=.3 and at higher Teff, then resolve the condensate enthalpy issue
   with a common thermochemical description or a quantitatively controlled
   alternative. Cold failed corners cannot simply be declared negligible.
   Include explicit input selectors and record exactly which family is used.
2. **Composition through exhaustion.** Extend EOS, opacity and atmosphere
   sources toward XH=0, including the atomic/baryonic isotope mapping. Source
   atmosphere abundances currently normalize to hydrogen; a zero-hydrogen
   endpoint requires a verified alternate normalization, not division by a
   small surrogate abundance. Extend temperature and density domains along
   measured profiles without clipping unsupported source states.
3. **Stratification and shell burning.** Add conservative microscopic
   diffusion/settling with the ambipolar electric field and appropriate
   collision physics, alongside the existing convection and slow mixing.
   Resolve moving boundaries and thin fuel layers with mass/energy-conserving
   remeshing if fixed-mesh convergence fails. Account for diffusion energetics
   and assess semiconvective heat transport. Establish convergence in endpoint
   age, residual hydrogen mass and the cooling delay.
4. **Energy production and losses.** Audit pep, ppIII and CNO contributions
   in newly encountered conditions; an inert GS98 metal inventory is not a
   CNO network. Check screening against the interaction EOS. Add and verify
   the relevant thermal-neutrino channels. The existing plasmon diagnostic
   does not establish a bound on all missing neutrino processes.
5. **Dense remnant thermodynamics.** Use a consistent free-energy treatment
   of degenerate electrons, Coulomb liquid/solid ions, mixing, ion quantum
   effects and latent heat, with validated joins to partially ionized matter.
   Compare independent source implementations. Do not impose a carbon/oxygen
   phase diagram or a classical melting threshold on a helium mixture.
6. **Cooling boundary and transport.** Follow surface H/He stratification
   from the transport solution. Dense cool atmospheres require nonideal
   chemistry, pressure-dependent absorption and collision-induced absorption;
   condensates, refraction and dense-fluid effects need relevance checks.
   Audit conductive transport in the actual quantum/solid regime and model
   the connection of envelope convection to the degenerate interior.
7. **Convergence and uncertainty.** Refine mesh, timesteps, source grids and
   endpoint thresholds; test reproducible fresh starts and exact restarts.
   Separate numerical errors, source interpolation errors and uncertainty in
   physical prescriptions. Archive receipts at every completed milestone.

## Source work under this authorization

- `hydrogen_poor_refinement_specification.json` requests independent TOPS
  eighth points at Z=.01/.02/.03 and the X=.25 midpoint. Previously independent
  quarter points may become candidate nodes only in a separately identified
  family whose acceptance uses fresh source points.
- `nongrey_hydrogen_poor_specification.json` defines a gas pilot at XH=.1/.2,
  X3=0/.12, Teff=3000..3600 K and logg=4.9/5.15/5.4. Its material rectangle
  retains the original 16 isotherms and adds 8855 and 1.010e4 K. These are
  prospective calculations, not installed atmospheres or validated condensate
  physics. The source pipeline checks the hottest rows first.
- These source jobs do not write the old condensate installation receipt or
  trigger its previously queued stellar controller.
- `docs/results/metal_eos_exhaustion_v1_audit.json` records a rejected
  36-plane candidate extending to true XH=0. Its largest cv error is 2.367%
  in cool, trace-H molecular gas. Direct-source curvature checks select 14
  additional hydrogen coordinates (28 planes). The refined64-plane v2
  passed528 fresh heldouts: P/E/cv/cp maxima .05213/.1088/.1841/.1922%.
  All2048 old stellar-profile queries reproduce the old EOS byte for byte.
  This validates composition interpolation, not dense or cold remnant matter.
  The source generator's exact abundance labels distinguish sub-per-mille
  composition planes.
- Both coarse and full-resolution XH=.1, Teff3400, logg5.15 atmosphere
  initializers converge numerically but exceed the 1.010e4 K opacity ceiling
  (deepest temperatures 1.345e4 and 1.342e4 K). Neither is accepted. A 15000 K
  opacity pilot also fails: He I 4471 profile extrapolation yields nonfinite
  log absorption already at rho=4.64e-7, ne=2.38e16. Merely excluding the
  highest-density source row does not solve the physical coverage problem.
- The new `--opacity-directory` driver option hashes the selected manifests
  and all referenced planes for snapshots and exact restarts. History now
  reports deposited nuclear and nuclear-neutrino luminosities, last-halfstep
  gravothermal luminosity and total hydrogen mass. These diagnostics do not
  establish a cooling endpoint. A subsequent explicit `--thermal-neutrinos`
  option now couples the HRW plasma prescription to structure and evolution;
  the original and current2.85T control calculations still select no thermal
  losses. Other channels remain omitted. See [NEUTRINOS.md](NEUTRINOS.md).

## Helium broadening source coverage

Before extending material opacities to accommodate a failed atmosphere,
check its **actual** bottom optical depth. TLUSTY's TAULAS applies to grey
initialization and does not replace the mass grid of a supplied structure.
New continuations reached optical depths 1.108e4–5.790e4 while matching at tau100.
`--initial-bottom-tau` explicitly truncates a verified seed using its measured
depth profile. The final optical depth is recorded separately. Comparisons
at multiple bottom columns must establish matching-state convergence before
adopting this computational-domain change; nominal input depth is insufficient.

The original SYNSPEC He I 4471 profile uses Barnard et al. data terminating
at ne=1e16 cm^-3; its polynomial extrapolation is not a valid extension for
the hotter, dense atmosphere pilot. Preserve the rejected runs and replace
the unsupported physical input before accepting a warmer boundary.

[Tremblay et al. (2026)](https://arxiv.org/abs/2603.04374) and their
[published datasets](https://zenodo.org/records/18722143) provide a useful
modern comparison. The supplied files cover 10000, 20000 and 40000 K through
ne=6e17 cm^-3. The simulation profiles assume He II ion perturbers, omit line
dissolution, and the complete 36-line file includes semi-analytic substitutes
for unsimulated transitions and low densities. A separate semi-analytic file
includes line dissolution. These distinctions matter in a hydrogen-bearing
mixture; this is not a universal plug-in replacement.

[Gigosos & Gonzalez (2009)](https://doi.org/10.1051/0004-6361/200912243)
[CDS tables](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A%2BA/503/293)
offer He I 4471 profiles for multiple ion reduced masses. Coverage is
nonrectangular: at ne >= 10^17.67 cm^-3 only 20000 and 40000 K are provided.
They therefore do not directly cover the pilot's deepest 13400 K,
ne=5.54e17 cm^-3 state. Doppler broadening is excluded from those tables.
Neither dataset has yet been installed in an atmosphere source executable.

## Condensate free-energy correction under investigation

The [FastChem Cond formulation](https://doi.org/10.1093/mnras/stad3515)
permits an equilibrium-condensation correction at fixed bulk composition.
Subtracting the gas-only Gibbs energy from the gas-plus-condensate result
cancels the reference atomic chemical potentials:

```text
Delta phi = Delta(g/T) = kB sum_j (N_j/M) ln(n_j,eq / n_j,gas)
x = ln T; y = ln Pgas
Delta v = (T/Pgas) Delta phi_y
Delta h = -T Delta phi_x
Delta s = -Delta phi - Delta phi_x
Delta cp = -Delta phi_x - Delta phi_xx
```

Here n_j is the neutral atomic number density; N_j/M is the fixed number
of nuclei per gram of bulk baryonic mass. The correction includes changes
to gas molecules and ions as well as condensed material. It assumes ideal
gas chemistry and condensates with negligible volume, retained in local
equilibrium. Rainout is a different thermodynamic/transport problem.

`scripts/audit_condensate_free_energy.py` independently checks the volume
derivative against species densities and the enthalpy/heat-capacity responses
against derivatives of all species' formation constants. Its three-state
pilot is recorded in `results/condensate_free_energy_audit_v1.json`.
At1769K/.004062bar the complete correction to cp is 1.175e5 erg/g/K;
the earlier grain-only estimate was289000. This difference demonstrates why
adding a standalone grain heat capacity would omit gas rearrangement.

No correction is installed yet. Adding it to a native gas EOS requires a
bulk-composition baseline, consistent specific volume and all thermodynamic
responses, and the corresponding opacity mass normalization. Source joins,
phase boundaries and deeper nonideal regimes require separate validation.
The original condensate production guard remains in effect.

The adapter now retains a defensive copy of its fixed bulk reservoir. An
initial synthetic audit fed depleted gas back as the next bulk input and
incorrectly identified this with the native Fortran call path. Inspection of
both prepared sources shows that MOLEQ already keeps `cbase` separate from
mutable `ccomp`/`abndd`; the latter is written back into depth arrays, while
`cbase` is initialized once and retained by `-fno-automatic`.

The corrected `audit_condensate_memory.py` compares both calling conventions
with fresh-process chemistry references. **Both old and new adapters agree
exactly on all twelve fixed-bulk calls.** Only the old adapter fails the
synthetic bulk-feedback experiment. Reports
`results/condensate_memory_{before,after}_v3_audit.json` supersede the earlier
interpretation. No history dependence of actual atmospheres has been
established, and this experiment does not invalidate their prior source
results. It adds no grain enthalpy or opacity. A future source build using the
defensive change still has a new executable identity and requires validation.

## Dense-matter references informing the implementation

Current derivative audits and the cold analytic-electron repair are recorded
in [DENSE_EOS.md](DENSE_EOS.md). They do not install a complete remnant EOS.

[Skye (Jermyn et al. 2021)](https://arxiv.org/abs/2104.00691) provides a
fully ionized free-energy EOS with Coulomb, mixture and quantum contributions
and a self-consistent crystallization calculation. It does not replace a
partially ionized atmospheric EOS.

The [Ioffe electron-ion plasma source](https://www.ioffe.ru/astro/EIP/eipintr.html)
provides another implementation and describes its quantum-liquid and mixture
limitations. Its default classical melting criterion is an explicit source
choice, which requires assessment for a helium-rich remnant.

[Blouin, Dufour & Allard (2018)](https://arxiv.org/abs/1807.06616) demonstrate
why the low-density approximations in ordinary stellar atmospheres are
insufficient for dense, cool helium atmospheres. This motivates a distinct
remnant boundary calculation; it does not establish coverage down to 100 K.
