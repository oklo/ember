# Two-phase thermodynamics

`equilibrate_phases` solves for two densities, both phase compositions and
their mass fractions. It conserves specific volume, H, He3 and metals, with
equal pressure and chemical potentials. He4 is the remaining mass fraction.
Each supplied phase evaluates one free energy and its first two derivatives.

The returned free-energy derivatives include changes in the phase fractions
and compositions. Thus the heat capacity includes phase adjustment; it is not
an average of the two frozen-phase heat capacities. Derivatives follow the
equilibrium equations by implicit differentiation. Unsupported trial states
shorten the Newton step; they do not extend either material model's domain.

This is a library calculation, **not yet selected by stellar evolution**.
The caller must choose the globally stable supported phase pair and provide
a nearby initial split. Phase-diagram search, single-phase boundaries,
zero-abundance limits and separation kinetics are separate responsibilities.
The API uses positive H, He3 and metal fractions, with He4 positive in each
phase. Potential units are F/(R T); coordinates are ln T, ln rho, X_H, X_He3, Z.

The analytic test uses two composition wells with a temperature-dependent
relative energy. Their exact common tangent tests the potential, thermal and
composition derivatives, phase fractions, volume and species conservation.
Nine states pass, including nonzero phase heat. Invalid compositions and
nonfinite material derivatives are rejected.

Independent private checks use four saved cold-core conditions with two
choices of solid mixing entropy. C++ heat capacities agree with the prototype
within 1.011e-11 relative. The largest difference among composition response
entries is 9.196e-6 relative, in a trace-He3 cross derivative; after scaling by
the corresponding diagonal responses the difference is below 1.481e-10.
The four solves take 0.1041 CPU seconds, with seven or eight iterations each.
That is a point-physics measurement, not a stellar runtime estimate.

Those physical probes formally extend the metal-rich phase beyond the current
Z=0.16 material limit. They assess a candidate approximation and do not validate
that extension for evolution. The effective GS98 metal group, quantum mixture
corrections and residual solid entropy still require physical comparisons.

## Material EOS derivatives

`phase_equilibrium_material` returns the material Helmholtz derivatives used by
the EOS: F/T, its thermal and density derivatives through third order, and the
H, He3 and shared-metal composition derivatives. It converts the dimensionless
phase potential to physical units. The supplied phase mixing terms are included;
radiation is not. Do not add mixing or latent heat a second time.

Four nearby coexistence solves differentiate the analytic Hessian in ln T and
ln rho, with a default increment of 1e-4. They start from the central split.
The same two phases must remain supported across that interval. A failed or
unsupported solve is reported, rather than replaced by a single phase. This is
a material evaluation and tabulation interface; phase selection, phase-boundary
handling and efficient reuse are still needed before selecting it in evolution.

Nine analytic controls with moving coexistence endpoints verify all populated
channels and the resulting pressure, energy, heat capacity, adiabatic gradient
and their derivatives. Eight saved-state controls agree with an independent
phase calculation; ultra-trace relative differences are assessed by their
abundance-weighted effect, not mistaken for a uniform relative-accuracy bound.
