# Electron conduction

The lifetime driver combines radiative opacity with the
[Ioffe conductive tables](../data/conduction/README.md):

```
kappa_cond = 16 sigma_SB T^3 / (3 rho K)
1/kappa_transport = 1/kappa_rad + 1/kappa_cond.
```

Temperature and density interpolation use monotone Hermite polynomials;
composition derivatives differentiate the same mixture evaluation. Source bounds
remain enforced. Optical depth in the atmosphere uses radiative opacity.

Each ion's conductivity is evaluated at the mixture's electron density, then
resistivities are weighted by its share of the electrons. This approximates
additive collision frequencies and counts electron-electron scattering once.
He3 and He4 use the same charge plane, with their different electron counts.
The GS98 option retains the elemental metal pattern. Ion correlations, isotope
transport and partial metal ionization remain approximations. Conductivities
at moderate degeneracy also depend on how the source prescriptions are joined;
see [Cassisi et al. (2021)](https://arxiv.org/abs/2108.11653).

`conduction_envelope "ionized"` uses the source conductivity directly. This
retains conduction in dense, cold matter that remains pressure-ionized. It
assumes ionization is complete where conduction affects the structure; it
supplies no electron-neutral scattering. The default `"temperature_join"`
multiplies conductance by a smooth temperature weight from zero at 300 kK to
one at 1 MK. That approximation can suppress real transport in a cold remnant.

A matched cooling comparison with the direct source gives 1844 K versus
1775 K at the same age. Direct ionization checks support the conducting layers.
Suppressing conduction in the cooler partly ionized layers changes luminosity
by at most 0.02819% over the completed comparison interval; that control then
stops at a convective-boundary iteration failure. Its uncompleted interval is
not counted as evidence. The EOS and cold-atmosphere limitations remain open.
[Calculations and limits](results/cold_conduction_cooling_sept30_v1.json).

The optional radiative-opacity continuation covers conducting H/He layers with
trace metals through `Z = 1e-4`. It uses two supported density anchors and admits
an extension only when its configured opacity uncertainty has a sufficiently
small effect on total heat transport. The selected cooling calculation allows
a factor of 100 in that opacity and at most 0.3% in transport. Direct EOS queries
support full H/He ionization in the tested layers. The change recovers 147
previously unsupported queries; a 100 Myr stellar comparison changes hydrogen
mass fraction by less than 3e-11. These checks do not establish the accuracy of
the conductivity source or the liquid-phase assumption.
[Checks and continued evolution](results/trace_metal_opacity_sept30_v1.json).
