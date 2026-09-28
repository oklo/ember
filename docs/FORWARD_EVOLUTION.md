# Composition and transport

The material EOS, atmosphere and opacity use a declared elemental mixture.
The variable-metal EOS interpolates a free-energy potential in temperature,
density, hydrogen, helium-3 and total metal fraction. Thermal responses,
chemical forces and transported enthalpies differentiate that potential.
Helium isotope mixing is included; the electronic source uses a helium-isotope
mapping. Radiation is added once. See [FreeEOS](FREEEOS.md) and the
[dense-EOS extensions](DENSE_EOS.md).

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


The common driver also supports microscopic H/He/metal transport in its assessed
hot ionized domain. Metals share a fixed mixture response. Cool neutral transport
and independent element settling require additional physics. Semiconvection and
thermohaline mixing are optional library prescriptions, not selected in the
current continuous trajectory.

Atmosphere, condensation and grain calculations are described in
[ATMOSPHERE.md](ATMOSPHERE.md), [NONGREY.md](NONGREY.md) and [GRAINS.md](GRAINS.md).
Stellar input choices and numerical controls are in [LIFETIME_DRIVER.md](LIFETIME_DRIVER.md).
