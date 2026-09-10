# Thermal-neutrino cooling

Thermal neutrinos now enter both the central mass element and every shell
through the same discrete first law:

```text
L_photon = L_nuclear_deposited + L_gravothermal - L_thermal_neutrino
L_nuclear_rest_mass_release = L_nuclear_deposited + L_nuclear_neutrino
```

Thermal losses are positive sinks. They do not consume nuclear fuel and do
not enter the nuclear rest-mass audit. Analytic temperature and density
responses enter the Henyey Jacobian. Nodal trapezoidal mass weights agree
between the structure equation and the integrated luminosity diagnostics.

The driver option `--thermal-neutrinos none|plasma-hrw` defaults to `none`,
preserving the explicit historical control. The selection is part of native
restart identity. JSON reports the selected prescription and its limitations,
and appends `thermal_neutrino_Lsun` to the existing history columns. Nuclear
and thermal luminosities are evaluated at the reported state; gravothermal
luminosity still refers to the last accepted implicit half step.

## Plasma prescription and its limits

`plasma-hrw` implements equations23--27 of
[Haft, Raffelt & Weiss (1994)](https://arxiv.org/abs/astro-ph/9309014), using
their numerical coefficients and weak mixing angle, sin²(theta_W)=.23.
The fully ionized electron number follows Ember's selected abundance basis
and metal inventory. Logarithmic evaluation avoids premature underflow of
very small emissivities. The published piecewise correction is retained;
derivatives are the analytic derivatives of the selected branch.

The paper's few-percent fit assessment concerns the **total** neutrino rate
when other processes are held fixed. It does not bound the relative error
of the plasma channel in regions where that channel is negligible. This
implementation is a plasma-only approximation: photo, pair, bremsstrahlung,
recombination and partial-ionization corrections are not supplied. It is
therefore insufficient by itself to certify a complete remnant cooling track.

The preserved2.85-trillion-year gas profile has an offline plasma luminosity
about1.45e-8 of its photon luminosity. That measures this prescription on
that profile; it is neither a total-loss bound nor a cooling-track result.

## Checks

`generate_neutrino_reference.py` uses the earlier independent scalar Python
implementation to produce21 published-fit comparison points. Runtime tests
check those values, temperature/density derivatives, the explicit zero-loss
control, and central/shell energy equations. A separate coupled evolution
test uses an appreciable analytic sink to check the complete energy balance
and the separation from nuclear rest-mass accounting. Driver tests check the
plasma history diagnostics, rejection of a changed loss selection, and exact
restart agreement with the same selection.

These checks validate the port and numerical coupling. Further physical
acceptance requires comparison with independent emissivity calculations and
assessment of all channels along the measured stellar/remnant trajectory.
