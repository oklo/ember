# Fingering convection

The library supplies the homogeneous prescription of
[Brown, Garaud & Stellmach (2013)](https://arxiv.org/abs/1212.1688), including
its inward sensible heat flux and extra composition diffusivity. The fastest
linear mode and its parameter derivatives are solved directly. The microscopic
composition flux is retained separately.

`fingering "brown_two_composition"` in the lifetime driver selects an
experimental cold-remnant extension. It requires `transport "screened_core"`;
the default is `"none"`. Hydrogen and helium-3 have separate diffusion
responses. Reducing them to one diffusivity can fail when their opposing
contributions to density nearly cancel.

The extension uses the EOS chemical response at fixed pressure, the native
collision mobility, and [liquid viscosity](VISCOSITY.md). Its two-composition
linear growth calculation recovers the Brown result for equal diffusivities.
Applying Brown's saturation constant to the two-field mode is an uncalibrated
approximation. It does not become a validated transport fit merely because the
linear calculation agrees with an independent eigensolver.

The same mode supplies inward sensible heat and a matrix of composition
transport coefficients. The species Newton solve updates these fluxes as the
composition changes. Each face exchanges equal and opposite amounts between
its cells; each actual CN isotope is mixed. Existing composition-enthalpy
accounting carries that energy once. Ordinary convection remains separate.

The coefficient model covers fully ionized, degenerate liquid H/He with
negligible metal buoyancy. It uses a classical mean-ion viscosity alongside
electron viscosity. It supplies no multicomponent closure outside that material
domain, and refuses a net inverse-composition layer that needs unsupported
coefficients. Opposing buoyancy contributions can support both stationary and
oscillatory modes. The fastest growth determines the classification.

By default an oscillatory fastest mode is refused. The explicit sensitivity
options `fingering_oscillation "growth_squared"` and `"growth_frequency"`
use dimensionless squared velocities `49 alpha²/q` and
`49 alpha sqrt(alpha²+omega²)/q`, respectively. Here `alpha+i omega` is
the fastest growth mode and `q` is its squared wavenumber. Both vanish at
marginal growth and recover the stationary law. Heat and composition use the
real, phase-dependent response to that mode. These are uncalibrated nonlinear
assumptions, not published oscillatory transport fits.

A matched 3 Gyr cooling comparison reaches 853.9 K with surface temperatures
differing by 0.0008866 K and luminosities by 0.0004169%. This establishes a
small sensitivity over that interval, not a general error bound. At 828.2 K,
the weaker choice requires shorter steps through a nonlinear composition
adjustment. Passage through the colder event and spatial convergence remain
unfinished. The model does not treat composition layers, crystal mixing,
rotation or magnetic fields, and is not yet a full-lifetime prescription.
The common executable reproduces the private 6 Myr comparison at the same age:
luminosity differs by 1.020e-9 and abundances by at most 1.121e-9. It took
smaller intervals because the central convective partition changed during
iteration; this comparison does not establish a performance improvement.
[Oscillatory checks](results/fingering_oscillatory_oct1_v1.json) ·
[Stationary and conservation checks](results/fingering_two_composition_oct1_v1.json).
