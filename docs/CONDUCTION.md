# Electron conduction in the extended 0.1 Msun model

The evolution driver includes radiative and electron-conductive transport by default. `TabulatedConduction` reads the [Ioffe source tables](../data/conduction/README.md). The implementation uses monotone Hermite interpolation in native log temperature and log density and linear interpolation in log ion charge. Density derivatives propagate through the temperature interpolant's limiter. All thermodynamic and H1/He3 composition derivatives differentiate the actual mixture evaluation.

The conversion from thermal conductivity K to conductive opacity and the combined opacity are

```
kappa_cond = 16 sigma_SB T^3 / (3 rho K)
1/kappa_transport = 1/kappa_rad + 1/kappa_cond.
```

The total enters the interior diffusion gradient and MLT element cooling. The atmosphere uses radiative opacity to define optical depth. This follows the transport framework of [Cassisi et al. (2007)](https://arxiv.org/abs/astro-ph/0703011). Conduction is substantial as a *possible diffusive carrier* in the old 20-Gyr core: kappa_rad≈51 and kappa_cond≈22 cm²/g. A convective star need not change its structure appreciably when conduction is included: efficient convection can supply the remaining flux at practically the same adiabatic gradient.

The conductivity provider exposes the common density support of every mapped
source ion used by its value and composition derivatives. Combined transport
intersects these bounds with the radiative source bounds; it does not imply
conductive support beyond the original temperature/density tables. The cold
zero-conductance branch imposes no additional conductive density restriction.

## Mixture approximation

The source tables describe pure, fully ionized species. For species j define electron fraction `f_j=Z_j n_j/n_e`. Evaluate each pure-ion conductivity at the **same electron density** as the mixture, using its source mass number:

```
rho_source,j = rho * Ye * A_source,j / Z_j
1/K_mix = sum_j f_j / K_j(T, rho_source,j).
```

This corresponds to an additive collision-frequency approximation: electron-ion scattering sums over species, while the electron-electron contribution appears once because the electron fractions sum to one. It preserves both pure-species and fixed-electron-density limits. For charges absent from the source axis, interpolation also holds electron density fixed. He3 and He4 share the helium source conductivity. Their abundance basis and electron counts enter the mapping explicitly.

The approximation transfers pure-ion correlation functions and mass-dependent corrections to a mixture. It is not a calculation of multicomponent structure factors or an isotope-resolved conductivity. Matthiessen addition and the partial-degeneracy continuation are approximate. The current forward calculation uses the shared GS98 elemental inventory in the ion and electron counts, with eleven source charge planes through Zn, including planes that bracket Ni. The retained legacy inventory uses the carried charges, including Ne-like `Zrest`. Neither option resolves partial metal ionization. See [FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md) for the current mixture and its remaining approximations. [Cassisi et al. (2021)](https://arxiv.org/abs/2108.11653) discuss significant uncertainty in the transition between nondegenerate and degenerate H/He plasmas; the weakly damped, uncorrected and undamped variants remain explicit comparison options.

## Neutral-envelope join

Full ionization is inappropriate in the cool envelope. `HotConduction` activates conductance with `w=u²(3−2u)` in log T between 3e5 and 1e6 K, then uses the full table. Below the interval the adopted conductance is zero (`kappa_cond=infinity`). Inside it, `kappa_cond` is divided by w and the derivative of w is included. This is an explicit domain join, not an ionization calculation. The high-temperature endpoint is chosen where H/He are effectively ionized in this star; metal ionization remains approximate.

The static comparisons move that interval to 6e5..2e6 K and also remove conduction or change its source prescription. Their structural effect is below the relaxation accuracy for the sampled fully convective models. This does not justify omitting conduction after a radiative core develops. Full neutrino cooling, crystallization and dense-ion/remnant physics are separate unfinished inputs.
