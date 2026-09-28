# Evolution equations

Evolution uses conserved baryonic abundances `X_i = A_integer,i Y_i`.
Their sum is one. Nuclear energy release uses physical atomic masses:

```
eps_deposited + eps_nuclear_neutrinos
    = -c^2 sum[(A_atomic/A_integer) dX_i/dt]
```

The mesh density and enclosed mass are baryonic. Newtonian gravity uses that
conserved mass; nuclear binding-energy changes do not update gravitating mass.
Composition-dependent thermal, ionization and interaction energies belong to
the material EOS, separately from nuclear rest mass.

The implicit thermal equation uses both old and new compositions. Without
additional losses or transported composition heat, its local form is

```
dL/dm = eps_deposited
    - [u_new-u_old + P_new*(1/rho_new-1/rho_old)]/dt.
```

Selected neutrino losses and composition heat transport add their corresponding
terms. Nuclear neutrino energy is already excluded from deposited heating.
Mass weights are shared by the thermal, nuclear and mixing calculations.

Structure, burning and transport are iterated together. Connected convective
regions mix according to the selected instantaneous or finite-rate prescription;
radiative regions retain local abundances and microscopic transport. Positivity,
isotope inventories and the discrete first law are checked independently.
Conservation accuracy is distinct from spatial and time accuracy.

The controller compares a full step with two half steps. Optional Richardson
extrapolation is accepted only with its physical and accounting checks; otherwise
it retains the two-half-step state. An abundance-change cap supplements the
error estimate. See [configuration and numerical options](LIFETIME_DRIVER.md),
[nuclear physics](NUCLEAR.md), [convection](CONVECTION.md) and
[restarts](RESTART.md).
