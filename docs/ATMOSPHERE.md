# Atmosphere boundary conditions

`Atmosphere::eval(Teff, gravity, composition)` supplies temperature, gas and
total pressure, density, matching optical depth, and logarithmic derivatives
of temperature and total pressure with respect to `Teff` and `gravity`.
All dimensional quantities use CGS. Both implementations use an EOS that
includes LTE radiation pressure, `a_rad*T^4/3`.

## Grey fallback

`GreyAtmosphere` implements the plane-parallel radiative Eddington relation
with locally evaluated **radiative** opacity. The defining temperature and
hydrostatic equations are also described in the
[MESA atmosphere documentation](https://docs.mesastar.org/en/latest/atm/t-tau.html):

```
T^4 = (3/4) Teff^4 (tau + 2/3)
dP/dtau = g/kappa.
```

With `Q = sigma_SB*Teff^4/c`, radiation pressure is `Q*(tau+2/3)`.
Consequently the equation for material pressure is

```
dPgas/dtau = g/kappa - Q.
```

The finite radiation pressure at the top must be retained. For constant
opacity the exact result is `P(tau) = (2/3)*Q + g*tau/kappa`; using
`P = g*tau/kappa` as the **total** pressure would be inconsistent with the
EOS. The implementation requires `g/kappa > Q`, excluding layers where
radiation supports or exceeds gravity.

Integration uses `x=ln(tau)`, `u=ln(Pgas)`, and two sensitivity equations for
`du/dln(Teff)` and `du/dln(g)`. The sensitivities include the EOS's density
response, the opacity derivatives, and the starting boundary condition.
Adaptive RK4 with step doubling controls the local error in all three
quantities. No finite differences of whole atmosphere integrations are needed
to obtain these derivatives.

The top condition is the local approximation
`Pgas_top = tau_top*(g/kappa_top - Q)`, solved with a self-consistent opacity.
Where the opacity declares density bounds, they bracket this solve; Newton
must not start outside the table merely because its initial guess was poor.
Ferguson now exposes those bounds. Actual domain violations still throw.

Default options are `tau_match=2/3`, `tau_top=1e-6`, and integration tolerance
`1e-8`. The omitted outer column is an independent approximation: reduce
`tau_top` to check convergence in a new regime. The code does not silently
move the top or extrapolate opacity to make a calculation finish. The
constant-opacity solution and a varying-opacity case with
`kappa proportional to P` test this distinction explicitly.

This fallback is radiative. It has no convective correction to `T(tau)`,
spherical extension, or irradiation, and the present analytic EOS has no
partial ionization or molecular hydrogen. Deep matching is mathematically
available but does not make a grey radiative profile into a convective model
atmosphere.

## Tabulated atmospheres

`TabulatedAtmosphere` loads one rectangular `(log10 Teff, log10 g)` grid at
one composition and matching optical depth. It interpolates `log10 T` and
`log10 Pgas` bilinearly, then adds radiation pressure and inverts the EOS.
Derivatives differentiate that same interpolation and pressure conversion.
They are cellwise, with possible jumps at grid lines. No corner overshoot,
table extrapolation, composition substitution, or implicit grey fallback is
performed.

The whitespace-delimited format is versioned:

```
EMBER_ATMOSPHERE 1
source "Model citation, version, mixture, and extraction provenance"
tau <matching Rosseland optical depth>
composition <H1 He3 He4 C12 C13 N14 O16 Zrest mass fractions>
log_teff <NT> <NT strictly increasing values>
log_g <NG> <NG strictly increasing values>
data
<log10 T> <log10 Pgas>
... NT*NG pairs, gravity varying fastest ...
```

Both axes need at least two values. Composition must sum to one; queries
must match each recorded species to absolute tolerance `1e-10`. Supporting
an evolving surface composition will require appropriately sourced grids and
composition interpolation. Malformed headers, duplicate axis values, missing
nodes, non-finite values, and extra trailing data are errors.

**No production atmosphere grid is bundled yet.** The only table is
`tests/data/synthetic_atmosphere.dat`, a labeled power-law test fixture.
Attempts on 2026-09-07 to retrieve the
[BT-Settl AGSS2009 structures](http://phoenix.ens-lyon.fr/Grids/BT-Settl/AGSS2009/STRUCTURES/)
timed out over HTTP and HTTPS. The public
[PHOENIX spectra description](https://lydu.ens-lyon.fr/phoenix/doc/spectra.html)
describes spectral flux files; those cannot supply an interior pressure
boundary. Importing structures remains work: verify optical-depth definition
(Rosseland versus a reference wavelength), gas versus total pressure,
composition, units, and model provenance before conversion. Version the
resulting physical numbers with the code.

## Surface equations

`surface_residual` evaluates two dimensionless residuals at the outer mesh
point:

```
f_T = ln T_mesh - ln T_atmosphere
f_P = ln P_EOS  - ln P_atmosphere
Teff^4 = L/(4*pi*sigma_SB*r^2)
g = G*M/r^2.
```

The returned `2x4` Jacobian uses the mesh variables `(ln r, ln rho, ln T, L)`
and includes both chains through `Teff` and `g`. Total mass and composition
are held fixed. A positive surface luminosity is required.

The matching point is explicitly the outer boundary, with its radius
identified with the photospheric radius and its enclosed mass with total
mass in the thin-atmosphere approximation. A thick envelope needs a resolved
atmospheric extension; this API does not place a fitting point at a fixed
interior mass fraction. These two equations are ready for assembly, but the
central boundary conditions and Henyey solver remain unimplemented.

## Findings from the integration tests

Connecting the existing modules exposed two errors worth preserving:

- `rho_from_PT` could exhaust its iterations and still return a density. It
  now requires convergence in density as well as pressure, rejects invalid or
  numerically unresolved states, and throws on nonconvergence.
- Ferguson's density derivative was interpolated as a separate field through
  the temperature interpolation. That is not the derivative of a limited
  interpolant: the slope limiter depends on the ordinates. The Hermite routine
  now propagates parameter derivatives through its active limiter branches.
  A direct opacity regression and the complete atmosphere sensitivity check
  both failed before this correction and pass afterward.
