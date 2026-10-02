# Phase equilibrium

`equilibrate_phases` conserves volume, H, He3 and metals while finding two or
three phases with equal pressure and chemical potentials. He4 is the remaining
mass fraction. Separate phases can use the same material model, allowing two
solid compositions to coexist with a liquid. Both APIs use the same solver.

Each material supplies F/(R T) and derivatives in ln T, ln rho, X_H, X_He3, Z.
The result differentiates the equilibrium constraints, including changing phase
fractions and compositions. Its heat capacity therefore includes phase heat.
Unsupported trials shorten the Newton step; accepted states retain the material
model's domain restrictions. Nonfinite derivatives are rejected.

This is a library component, **not yet selected by stellar evolution**. The
caller must choose the stable supported phases and provide an initial split.
All species and phase fractions must be positive. Global phase selection,
vanishing species/phases and separation kinetics remain the caller's work.
A converged stationary solution alone does not establish global stability.

## Material derivatives

`phase_equilibrium_material` returns F/T derivatives through third thermal and
density order and second composition order. It includes the supplied mixing
terms and excludes radiation. Do not add mixing or latent heat again.

When every phase supplies third derivatives, the method differentiates the
same equilibrium equations a second time. One matrix solve with 15 right-hand
sides gives the quadratic responses. No displaced coexistence solves or
finite-difference increments are needed. This is useful at narrow transitions,
where differencing large thermal terms can obscure a small transport response.

Otherwise, four nearby coexistence solves differentiate the Hessian. They use
predicted densities, compositions and phase fractions as starting guesses.
The multi-phase interface shortens the requested increment near a phase
boundary; its returned `log_steps` records the increments actually used.
With analytic derivatives both increments are zero and `analytic_third` is true.
Failed or unsupported neighboring solves are reported, not replaced by another
phase set. Finite differences can remain ill-conditioned at narrow transitions.

## Checks and limits

Exact common tangents and a common plane through three moving composition
wells test conservation, phase fractions, potential derivatives and material
responses. Controls include changing phase order, trace species, a stiff thermal
transition and both derivative methods. The stiff third-derivative check uses a
0.1% bound because it subtracts large terms to obtain an order-unity response.

Private cold-interior controls compare C++ with an independent implementation.
The narrow three-phase transitions conserve latent energy when heat capacity
is integrated across them. These controls assess numerical thermodynamics;
they do not validate the formal metal-rich extension of the effective GS98
mixture. Element-specific partitioning and solid entropy remain physical
uncertainties. The numerical summary is in `results/phase_mixture_oct2.json`.
