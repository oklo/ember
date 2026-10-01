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

The dense material range extends to 100,000 K at densities above 300 g/cm³,
subject to the viscosity, EOS and collision limits. Direct source ionization
checks over 100,000–120,000 K find hydrogen ionized and the doubly ionized helium
fraction above 0.9999. In the tested hydrogen-rich mixing layer, omitting or
doubling the classical ion viscosity changes the additional mixing by at most
0.02602%. This is a component sensitivity, not an error bound on quantum
mixture diffusion or nonlinear saturation.

A matched 6 Myr comparison leaves the warmer stellar state unchanged.
Colder projected queries recover the mixing layer without changing previously
supported responses; separate mixed-H/He EOS limits remain. The model does
not treat composition layers, crystal mixing, rotation or magnetic fields.
A final full-lifetime validation remains unfinished.

[Cold material checks](results/fingering_cold100_oct1_v1.json) ·
[Oscillatory checks](results/fingering_oscillatory_oct1_v1.json) ·
[Stationary and conservation checks](results/fingering_two_composition_oct1_v1.json).
