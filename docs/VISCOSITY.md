# Liquid viscosity

`ocp_liquid_viscosity` returns electron and classical-ion shear viscosity and
analytic derivatives for a fully ionized, nonmagnetic one-component liquid.
Kinematic viscosity is their sum divided by density.

The electron-ion term integrates the effective scattering potential of
[Chugunov & Yakovlev (2005)](https://arxiv.org/abs/astro-ph/0511300). The
implementation includes their electron-electron term and restricts use to
nonrelativistic degenerate electrons: momentum below 0.5 electron rest-mass
units, and temperature below 0.05 of the Fermi temperature. The reported fit
accuracy does not bound the underlying plasma-model uncertainty.

The ion term uses the classical liquid fit of
[Daligault, Rasmussen & Baalrud (2014)](https://arxiv.org/abs/1407.3875), through
Coulomb coupling 175. It remains explicitly classical when the returned ion
quantum parameter approaches unity. This interface does not supply an
ionization model, a mixture rule or solid rheology.

At the hydrogen-rich layer tested during cold cooling, electrons contribute
98.11% of the estimated viscosity. Omitting or doubling the classical ion term
changes the total by 1.892%; this is a sensitivity test, not an uncertainty bound.

The EOS also supplies `isobaric_composition_response`: density and chemical
responses at fixed temperature and pressure, including the density adjustment.
Multiplying its chemical matrix by the microscopic mobility and dividing by
density gives the local composition diffusion matrix. A scalar reduction for
fingering convection needs a separate check when composition contributions to
buoyancy nearly cancel.

An optional [cold fingering treatment](FINGERING.md) uses these components;
it remains experimental and is off by default.
