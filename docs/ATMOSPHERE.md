# Atmosphere boundary conditions

`Atmosphere::eval(Teff, gravity, composition)` supplies temperature, gas and
total pressure, density, matching optical depth, and logarithmic derivatives
of temperature and total pressure with respect to `Teff` and `gravity`.
All dimensional quantities use CGS. The implementations use an EOS that
includes LTE radiation pressure, `a_rad*T^4/3`.

The optional `CompositionAtmosphereGrid` now supplies an evolving H/He
boundary from 48 independently calculated non-grey radiative/convective
atmospheres. Its composition mapping, source coverage, interpolation and
physical checks are documented in [NONGREY.md](NONGREY.md). The COND-based
composition correction remains the default comparison backend.

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
table extrapolation, silent composition substitution, or implicit grey fallback is
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
must match each recorded species to absolute tolerance `1e-10`. Evolving
surface composition is handled by the separate composition atmosphere
backends. Malformed headers, duplicate axis values, missing
nodes, non-finite values, and extra trailing data are errors.

Version 2 adds `composition_proxy "description"` between `tau` and
`composition`. In this version `composition` is the permitted interior
query, not a claim about the native atmosphere abundances. It requires
explicit construction with `Mixture::allow_documented_proxy`. Version 1
retains the exact-composition contract and existing synthetic tests.

**A physical AMES-COND-2000 boundary grid is now bundled**, extracted from
checksum-pinned MESA tau=100 material-pressure and temperature data.
The selected Teff=1800..3300 K, log g=3.5..6 rectangle excludes all MESA
grey-fill, extrapolation and transition-smoothing cells. The native GN93
solar atmosphere is an explicitly declared proxy for ember's X=.7,Z=.02
interior, with unmatched detailed mixture and helium abundance. See
[data provenance and reproduction](../data/atmosphere/README.md).

Select it with `--atmosphere cond-solar-proxy`. `--tau-top` applies only
to grey integrations and is rejected with this option. At Teff=2800 K,
log g=5, the original tau=100 state is T=4081.407 K and
Pgas=1.555799e7 dyn/cm². The stellar seed now uses this local matching
temperature; using Teff at a deep boundary was incorrect. Tests check
source values, analytic sensitivities, explicit proxy selection and bounds.
The imported states come from structures, not spectral flux files.

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
interior mass fraction. These equations are now assembled with the central
conditions and interior rows by the Henyey relaxation solver; see `HENYEY.md`.
The first complete continuum benchmark uses an explicitly artificial
polytropic atmosphere, not a physical atmosphere grid.

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

## Composition-dependent convective column

`ConvectiveAtmosphere` extends the grey closure by solving temperature and hydrostatic pressure together, using the evolving EOS and radiative opacity. `CompositionCorrectedAtmosphere` uses its differential composition response to extend the COND boundary. It does not supply a new non-grey composition grid.

The top uses the existing Eddington integration from tau=.001 to .01. Below that point, the independent variable is ln tau and the column solves

```
dlnP/dln(tau) = tau g / (kappa_rad P)
dlnT/dln(tau) = grad * dlnP/dln(tau)
grad_rad = 3 kappa_rad P F / (16 sigma_SB T^4 g)
F = sigma_SB Teff^4.
```

Convection uses the Henyey element-cooling prescription in equations 13–18 of [Gustafsson et al. (2008)](https://arxiv.org/abs/0805.0554). Let `ell=alpha Hp`, `v0=alpha sqrt(g Hp delta/8)` and `tau_ell=rho kappa ell`. With `q²=grad−grad_element`, the positive cubic is

```
grad_rad - grad_ad = q² + B q + C q³
B = 8 sigma_SB T³ tau_ell / [v0 rho cp (1+y tau_ell²)]
C = 3 kappa P rho cp alpha v0 / (32 sigma_SB T³ g).
```

The terms represent the element contrast, element cooling, and convective flux. The adopted alpha=1.9 and y=1/3 recover the interior's Bohm–Vitense coefficients in the optically thick limit. The alternative y=.076 is a sensitivity control. Scaling the cubic bounds its root without subtracting nearly equal gradients. The two Teff/gravity sensitivity equations propagate the EOS, opacity, buoyancy, heat capacity and cubic derivatives analytically.

The default integration tolerance is 2e-8 for both logarithmic states and sensitivities. Sensitivities must also be well resolved: loosening them can allow boundary-value noise that prevents a tight stellar Newton solve from converging. Separate controls permit convergence tests; normal evolution keeps both tight. The EOS inversion now brackets within its declared density interval, so a poor ideal-gas density guess cannot cross the table floor.

For a declared reference mixture c0, with baryonic X=.7, Z=.02 and zero He3, the corrected boundary is

```
T_match(c)    = T_COND * T_column(c) / T_column(c0)
Pgas_match(c) = Pgas_COND * Pgas_column(c) / Pgas_column(c0).
```

Radiation pressure is then added once at the corrected T, and density is inverted using the actual composition. Logarithmic derivatives combine the same three evaluations; the radiation term uses the corrected temperature derivative. Source tau must agree (100 here). At c0, both values and derivatives reproduce COND exactly. Fixed metal abundances, every source composition/density bound, and the original COND Teff/log-g rectangle remain required.

This assumes that the non-grey correction at c0 transfers to other mixtures. Radiative transfer is still represented by a grey closure, and isotope-dependent collision-induced absorption and line broadening are absent. The anchor covers Teff=1800..3300 K and log g=3.5..6; the EOS and column impose additional restrictions, including the temperature of the upper atmosphere. There is no fallback outside their common support.

The independent checks recover the constant-opacity Eddington limit, compare derivatives through convective layers, verify exact reference anchoring and test the top-column truncation. Static comparisons at four homogeneous compositions are in [the boundary sensitivity results](results/extended_boundary_sensitivity.json). Removing the COND anchor changes luminosity by 7–14%; changing y, alpha or the starting optical depth has a smaller effect. These are sensitivity experiments, **not a calibrated physical error bar**. A helium-enriched non-grey atmosphere remains the preferred next replacement.

## Non-grey composition grid

`CompositionAtmosphereGrid` and the pinned source-generation/import pipeline
are documented in [NONGREY.md](NONGREY.md). Runtime values and derivatives
come from the same four-dimensional interpolant in H1, He3, Teff and gravity.
Source support and mixture assumptions are explicit; the reader never
fills holes, extrapolates or substitutes a grey boundary.
