# Analytic zone Jacobian

`zone_residual` returns the four structure residuals and their derivatives
with respect to `(ln r, ln rho, ln T, L)` at both zone endpoints. Mass
coordinates, composition, mixing-length parameter, time step, and the
previous model are held fixed. A small forward chain-rule type carries the
eight derivatives through the zone algebra. Physics modules supply analytic
partials; they are not evaluated at perturbed states in this path.

Each endpoint needs one EOS response, one opacity evaluation, and one nuclear
evaluation. `zone_equations` evaluates values only. The independently
differenced `zone_residual_numerical` remains a reference for tests and for
modules that do not yet provide the complete analytic response. It evaluates
the same residual expressions with ordinary EOS calls, without using the
EOS transport derivatives. Its luminosity perturbation uses the zone's flux
and heating scale, with a 1 erg/s floor rather than a fraction of solar
luminosity.

## EOS transport response

`Eos::eval_with_derivatives` returns `EosResponse`: the ordinary state plus
derivatives of `cp`, `delta`, and `grad_ad` with respect to `ln T` and
`ln rho`, and `dE/dln rho`. The energy temperature derivative is already
`cv*T`. These are ordinary derivatives with respect to logarithmic inputs,
not logarithmic derivatives of the outputs.

The default response method throws. A new EOS must explicitly implement it;
missing derivatives must not silently become zeros. Ordinary `eval` remains
available independently. `CompositeEos` adds component Hessians before
forming thermodynamic response functions. Components provide

```
P_TT, P_TR, P_RR, E_TT, E_TR,
```

where `T` and `R` subscripts denote logarithmic temperature and density
partials at fixed composition. From these it differentiates

```
delta   = P_T / P_R
cv      = E_T / T
cp      = cv + P_T*delta/(rho*T)
grad_ad = P*delta/(rho*T*cp).
```

Ions and radiation have explicit power-law Hessians. `IdealEos` shares the
composite transport-response implementation while retaining its independent
ordinary state calculation for value regression tests.

## Degenerate electron derivatives

The electron chemical potential is implicit: electron number density must
match the composition and mass density. Differentiating that constraint
gives the first and second occupation responses. Directly subtracting large
thermal terms under strong degeneracy loses the small response that
convection needs, so the second derivatives use centered number-weighted
moments.

With `g` the kinetic energy in units of electron rest energy,
`beta=kT/(m_e*c^2)`, and occupation `f`, define

```
h = g/beta + d eta/d ln T
s = 1 - 2f
a = d eta/d ln rho.
```

Let angle brackets denote averages over the number kernel times `f(1-f)`,
and put `B=<s*h^2>`, `C=<s*h>`, `D=<s>`. The second occupation responses,
after removing the common factor `f(1-f)`, are

```
TT:  s*h^2 - h - B
TR:  a*(s*h - C)
RR:  a*(1 + a*(s - D)).
```

The pressure and energy kernels are centered at the Fermi surface by
subtracting a constant multiple of the number kernel. The number constraint
makes its TT and TR contributions zero; the RR number contribution is
restored explicitly. Specific energy also requires differentiating its
explicit `1/rho` factor.

The quadrature shell is split on both sides of the Fermi surface. A single
64-node shell gave plausible pressures but caused roughly 0.2% derivative
discrepancies in degenerate tests, which did not shrink with the finite
difference step. Splitting the shell resolved this; increasing the test
tolerance would have hidden the problem. The occupation derivative is
evaluated as `z/(1+z)^2`, with `z=exp(-abs(g/beta-eta))`, to preserve it even
when the occupation rounds to one.

## Heating and transport

The pp heating derivatives use the same nuclide mass defects and escaping
neutrino energies as the heating value. Screening is differentiated too:
below the existing weak-screening cap, its logarithmic temperature and
density derivatives are `-3H/2` and `H/2`; above the cap they vanish.
This corrects the earlier approximate branch weights and omitted screening
response. It does not replace weak screening with a dense-matter model.

The convection chain includes all state dependence of `grad_rad`, `grad_ad`,
and `U`, including pressure scale height, `cp`, `delta`, and opacity. It uses
the partials returned by the bounded MLT cubic. The transport residual stays
in plain gradient form, without multiplication by `grad/grad_rad`.

Positive `dt` requires a previous model on the identical mass mesh. The
existing gravitational heating discretization is retained:

```
eps_grav = -(E - E_old - P/rho^2*(rho-rho_old))/dt.
```

It is evaluated at the left endpoint of the zone. Its analytic derivative
holds the previous state fixed. This is still a backward, left-endpoint
energy difference; Jacobian assembly does not establish its time accuracy
or supply a time integrator. Non-positive `dt` requests a static model.

## Checks and remaining work

The EOS suite compares responses with ordinary EOS evaluations across ideal,
radiation-dominated, low-mass interior, and degenerate regimes, including
`rho=1e8 g/cm^3`. Radiation-limit power laws provide independent checks.
Nuclear tests vary branch mixtures and screening strength. Structure tests
compare the full matrix at both endpoints in burning and contracting
interiors, degenerate helium, finite-efficiency convection, and zero or
inward luminosity. An ideal-gas compression test checks the sign and
magnitude of gravitational heating independently. The efficient-convection
regression still verifies temperature-row conditioning below
`grad/grad_rad=1e-6`.

Central boundary conditions and Henyey block elimination are the next solver
work. Neither an equilibrium stellar solve nor an evolutionary run is
implemented yet. Physical atmosphere grids, additional opacity coverage,
and the other pending physics in `ROADMAP.md` remain necessary for the
0.1 solar-mass end-to-end milestone.
