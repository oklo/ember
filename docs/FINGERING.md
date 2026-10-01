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
oscillatory modes. The fastest growth determines the classification; a dominant
oscillatory mode is refused because its transport is not supplied. This is not yet a full-lifetime fingering prescription, a crystal
mixing law, or a treatment of composition layers, rotation or magnetic fields.

A warm-state diagnostic advances 100 Myr from 848.7 K with the original
conservation and time-error checks. It finishes at 847.6 K. The longer warm replay encounters an oscillatory branch in trial states.
Passage through the colder mixing event and spatial convergence remain unproven.
[Checks and limits](results/fingering_two_composition_oct1_v1.json).
